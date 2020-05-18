#pragma once

#include "../datastructure/flow_hypergraph_builder.h"
#include "hypergraph.h"
#include "../datastructure/queue.h"
#include "partition_ca.h"
#include "../datastructure/node_border.h"
#include "config.h"
#include "partition_threadsafe.h"
#include "mock_builder.h"

namespace whfc_rb {
    class FlowHypergraphBuilderExtractor {
    public:
        static constexpr NodeID invalid_node = std::numeric_limits<NodeID>::max();
        whfc::FlowHypergraphBuilder fhgb;

        static constexpr bool doBothBFSInParallel = false;

        FlowHypergraphBuilderExtractor(const size_t maxNumNodes, const size_t maxNumEdges, const size_t maxNumPins, int seed, const PartitionConfig& config) :
                fhgb(2 * config.percentage_bfs_from_cut * maxNumNodes + 2, maxNumEdges),
                visitedNode(maxNumNodes), visitedHyperedge(maxNumEdges),
                globalToLocalID(maxNumNodes),
                localToGlobalID(maxNumNodes + 2),
                mt(seed), config(config),
                mockBuilder_thread_specific(std::ref(fhgb.getNodes()), whfc::invalidNode, whfc::invalidNode),
                queue_thread_specific_1(maxNumNodes), queue_thread_specific_2(maxNumNodes) {}

        struct ExtractorInfo {
            whfc::Node source;
            whfc::Node target;
            whfc::Flow baseCut;
            whfc::Flow cutAtStake;
        };

        template<class PartitionImpl>
        ExtractorInfo
        run(PartitionImpl &partition, const PartitionBase::PartitionID part0, const PartitionBase::PartitionID part1,
            NodeWeight maxW0, NodeWeight maxW1, whfc::DistanceFromCut& distanceFromCut, whfc::TimeReporter& timer) {

            timer.start("Init", "Extraction");
            CSRHypergraph &hg = partition.getGraph();
            initialize(hg.numNodes(), hg.numHyperedges());
            timer.stop("Init");

            whfc::HopDistance delta = config.distancePiercing ? 1 : 0;

            timer.start("Filter", "Extraction");
            auto cut_hes = partition.getCutEdges(part0, part1);
            timer.stop("Filter");

            timer.start("Shuffle", "Extraction");
            cut_hes.shuffle(mt);
            timer.stop("Shuffle");

            timer.start("BFS", "Extraction");

            whfc::NodeWeight w0, w1;

            localToGlobalID[0] = invalid_node;
            fhgb.addNode(whfc::NodeWeight(0));
            result.target = whfc::Node(1);
            fhgb.addNode(whfc::NodeWeight(0));
            localToGlobalID[result.target] = invalid_node;

            timer.start("Scan_hyperedges_and_add_nodes", "BFS");
            tbb::parallel_invoke([&]() {
                w0 = BreadthFirstSearch(hg, cut_hes, partition, part0, part1, maxW0, result.source, -delta, distanceFromCut, timer, queue_thread_specific_1);
            }, [&]() {
                w1 = BreadthFirstSearch(hg, cut_hes, partition, part1, part0, maxW1, result.target, delta, distanceFromCut, timer, queue_thread_specific_2);
            });
            timer.stop("Scan_hyperedges_and_add_nodes");

            timer.start("writeQueuesToBuilder", "BFS");
            writeQueuesToBuilder(fhgb, queue_thread_specific_1, hg, distanceFromCut);
            writeQueuesToBuilder(fhgb, queue_thread_specific_2, hg, distanceFromCut);
            timer.stop("writeQueuesToBuilder");

            visitedHyperedge.reset();

            std::vector<HyperedgeWithSize> hyperedges;
            size_t total_num_hyperedges = 0;
            for (std::vector<HyperedgeWithSize>& hyperedges_local : hyperedges_thread_specific) {
                total_num_hyperedges += hyperedges_local.size();
            }

            hyperedges.reserve(total_num_hyperedges);
            for (std::vector<HyperedgeWithSize>& hyperedges_local : hyperedges_thread_specific) {
                hyperedges.insert(hyperedges.end(), hyperedges_local.begin(), hyperedges_local.end());
                hyperedges_local.clear();
            }

            

            timer.start("Add_hyperedges", "BFS");
            tbb::parallel_invoke([&]() {
                addHyperedges(hg, partition, part0, part1, timer, queue_thread_specific_1);
            }, [&]() {
                addHyperedges(hg, partition, part0, part1, timer, queue_thread_specific_2);
            });
            timer.stop("Add_hyperedges");

            timer.start("addMockBuildersParallel", "BFS");
            fhgb.addMockBuildersParallel(mockBuilder_thread_specific);
            timer.stop("addMockBuildersParallel");

            for (whfc::Flow& baseCut : baseCut_thread_specific) {
                result.baseCut += baseCut;
                baseCut = 0;
            }

            for (whfc::Flow& cutAtStake : cutAtStake_thread_specific) {
                result.cutAtStake += cutAtStake;
                cutAtStake = 0;
            }


            timer.stop("BFS");


            fhgb.nodeWeight(result.source) = partition.partWeight(part0) - w0;
            fhgb.nodeWeight(result.target) = partition.partWeight(part1) - w1;

            timer.start("Finalize", "Extraction");
            fhgb.finalize();
            timer.stop("Finalize");

            //fhgb.printHypergraph(std::cout);

            return result;
        }

        auto localNodeIDs() const {
            return boost::irange<whfc::Node>(whfc::Node(0), whfc::Node(fhgb.numNodes()));
        }

        whfc::Node global2local(const NodeID x) const {
            assert(visitedNode.isSet(x));
            return globalToLocalID[x];
        }

        NodeID local2global(const whfc::Node x) const { return localToGlobalID[x]; }

    private:
        struct NodeWithDistance {
            whfc::Node node;
            whfc::HopDistance d;
        };

        struct HyperedgeWithSize {
            HyperedgeID e;
            size_t pin_count;
        };

        ldc::AtomicTimestampSet<> visitedNode;
        ldc::AtomicTimestampSet<> visitedHyperedge;
        std::vector<whfc::Node> globalToLocalID;
        std::vector<NodeID> localToGlobalID;
        ExtractorInfo result;
        std::mt19937 mt;
        const PartitionConfig& config;

        MockBuilder mock_builder;
        tbb::enumerable_thread_specific<MockBuilder> mockBuilder_thread_specific;
        tbb::enumerable_thread_specific<whfc::Flow> baseCut_thread_specific;
        tbb::enumerable_thread_specific<whfc::Flow> cutAtStake_thread_specific;

        tbb::enumerable_thread_specific<LayeredQueue<NodeWithDistance>> queue_thread_specific_1;
        tbb::enumerable_thread_specific<LayeredQueue<NodeWithDistance>> queue_thread_specific_2;

        tbb::enumerable_thread_specific<std::vector<HyperedgeWithSize>> hyperedges_thread_specific;


        template<class Builder>
        bool writeQueuesToBuilder(Builder& builder, tbb::enumerable_thread_specific<LayeredQueue<NodeWithDistance>>& queues,
                CSRHypergraph& hg, whfc::DistanceFromCut& distanceFromCut) {
            std::atomic<size_t> node_counter = builder.numNodes();
            size_t num_nodes = builder.numNodes();
            std::vector<whfc::FlowHypergraph::NodeData>& nodes = builder.getNodes();

            for (LayeredQueue<NodeWithDistance>& queue : queues) {
                num_nodes += queue.queueEnd();
            }

            bool has_nodes = num_nodes > builder.numNodes();

            // Write node data into builder
            nodes.resize(num_nodes + 1); // maybe use reserve

            tbb::parallel_for_each(queues, [&](LayeredQueue<NodeWithDistance>& queue) {
                size_t begin_index = node_counter.fetch_add(queue.queueEnd());
                tbb::parallel_for(tbb::blocked_range<size_t>(0, queue.queueEnd()), [&](const tbb::blocked_range<size_t>& indices) {
                    for (size_t i = indices.begin(); i < indices.end(); ++i) {
                        const whfc::Node localID(begin_index + i);
                        const NodeWithDistance globalNode = queue.elementAt(i);
                        assert(localID < num_nodes);
                        nodes[localID].weight = whfc::NodeWeight(hg.nodeWeight(globalNode.node));
                        globalToLocalID[globalNode.node] = localID;
                        localToGlobalID[localID] = globalNode.node;
                        distanceFromCut[localID] = globalNode.d;
                    }
                });
            });
            assert(num_nodes == node_counter.load());

            return has_nodes;
        }

        template<typename PartitionImpl>
        inline void processHyperedge(const CSRHypergraph& hg, const whfc_rb::HyperedgeID& e,  PartitionImpl& partition, PartitionBase::PartitionID partID,
                std::atomic<whfc::NodeWeight> &w, NodeWeight& lastSeenValue, NodeWeight maxWeight, whfc::HopDistance d, LayeredQueue<NodeWithDistance>& queue,
                std::vector<HyperedgeWithSize>& hyperedges) {
            size_t num_pins = 0;
            for (NodeID u : hg.pinsOf(e)) {
                if (partition[u] == partID) {
                    lastSeenValue = w.load(std::memory_order_relaxed);
                    if (lastSeenValue + hg.nodeWeight(u) <= maxWeight) {
                        if (visitedNode.set(u)) {
                            w += hg.nodeWeight(u);
                            queue.push({static_cast<whfc::Node>(u), d});
                        }
                        num_pins++;
                    }
                }
            }
            hyperedges.push_back({e, num_pins});
        }

        template<class PartitionImpl, typename CutEdgeRange>
        whfc::NodeWeight BreadthFirstSearch(CSRHypergraph &hg, CutEdgeRange &cut_hes, const PartitionImpl &partition,
                                            PartitionBase::PartitionID partID, PartitionBase::PartitionID otherPartID,
                                            NodeWeight maxWeight, whfc::Node terminal, whfc::HopDistance delta,
                                            whfc::DistanceFromCut& distanceFromCut, whfc::TimeReporter& timer,
                                            tbb::enumerable_thread_specific<LayeredQueue<NodeWithDistance>>& queue_thread_specific) {
            std::atomic<whfc::NodeWeight> w = 0;
            whfc::HopDistance d = delta;


            // Collect boundary vertices
            //timer.start("Collect_boundary_vertices", "BFS");
            tbb::blocked_range<size_t> range(0, cut_hes.size(), 1000);
            tbb::parallel_for(range, [&](const tbb::blocked_range<size_t>& sub_range) {
                NodeWeight lastSeenValue = 0;
                LayeredQueue<NodeWithDistance>& queue = queue_thread_specific.local();
                std::vector<HyperedgeWithSize>& hyperedges = hyperedges_thread_specific.local();
                for (size_t i = sub_range.begin(); i < sub_range.end(); ++i) {
                    const HyperedgeID e = cut_hes[i];
                    for (NodeID v : hg.pinsOf(e)) {
                        if (partition[v] == partID) {
                            lastSeenValue = w.load(std::memory_order_relaxed);
                            if (lastSeenValue + hg.nodeWeight(v) <= maxWeight) {
                                processHyperedge(hg, e, partition, partID, w, lastSeenValue, maxWeight, d, queue, hyperedges);
                            }
                        }
                    }
                    if (lastSeenValue >= maxWeight) break;
                }
            });
            //timer.stop("Collect_boundary_vertices");

            d += delta;

            bool nodes_left = false;

            for (LayeredQueue<NodeWithDistance>& queue : queue_thread_specific) {
                queue.finishNextLayer();
                if (!queue.currentLayerEmpty()) {
                    nodes_left = true;
                }
            }

            while (nodes_left) {
                //timer.start("Scan_hyperedges_and_add_nodes", "BFS");
                // scan hyperedges and add nodes to the next layer
                tbb::parallel_for_each(queue_thread_specific, [&](LayeredQueue<NodeWithDistance>& queue) {
                    tbb::parallel_for(tbb::blocked_range<size_t>(queue.currentLayerStart(), queue.currentLayerEnd()), [&](const tbb::blocked_range<size_t>& indices) {
                        LayeredQueue<NodeWithDistance>& local_queue = queue_thread_specific.local();
                        std::vector<HyperedgeWithSize>& hyperedges = hyperedges_thread_specific.local();
                        NodeWeight lastSeenValue = 0;

                        for (size_t i = indices.begin(); i < indices.end(); ++i) {
                            const NodeID u = queue.elementAt(i).node;

                            for (HyperedgeID e : hg.hyperedgesOf(u)) {
                                if (partition.pinsInPart(otherPartID, e) == 0 && partition.pinsInPart(partID, e) > 1 && visitedHyperedge.set(e)) {
                                    processHyperedge(hg, e, partition, partID, w, lastSeenValue, maxWeight, d, local_queue, hyperedges);
                                }
                                if (lastSeenValue >= maxWeight) break;
                            }
                        }
                    });
                });
                //timer.stop("Scan_hyperedges_and_add_nodes");

                nodes_left = false;

                for (LayeredQueue<NodeWithDistance>& queue : queue_thread_specific) {
                    queue.clearCurrentLayer();
                    queue.finishNextLayer();
                    if (queue.currentLayer().size() > 0) {
                        nodes_left = true;
                    }
                }

                d += delta;

            }

            distanceFromCut[terminal] = d;


            return w.load();
        }

        template<typename PartitionImpl>
        void addHyperedges(CSRHypergraph& hg, const PartitionImpl &partition, PartitionBase::PartitionID part0,
                PartitionBase::PartitionID part1, whfc::TimeReporter& timer, tbb::enumerable_thread_specific<LayeredQueue<NodeWithDistance>>& queue_thread_specific) {

            //timer.start("Add_hyperedges", "BFS");
            tbb::parallel_for_each(queue_thread_specific, [&](LayeredQueue<NodeWithDistance>& queue) {
                tbb::parallel_for(tbb::blocked_range<size_t>(0, queue.queueEnd()), [&](const tbb::blocked_range<size_t>& indices) {
                    auto &local_builder = mockBuilder_thread_specific.local();
                    local_builder.setSource(result.source);
                    local_builder.setTarget(result.target);
                    whfc::Flow &baseCut = baseCut_thread_specific.local();
                    whfc::Flow &cutAtStake = cutAtStake_thread_specific.local();

                    for (size_t i = indices.begin(); i < indices.end(); ++i) {
                        const NodeID u = queue.elementAt(i).node;

                        // skip terminals
                        if (u == invalid_node) {
                            continue;
                        }

                        const PartitionBase::PartitionID partID = partition[u];
                        const PartitionBase::PartitionID otherPartID = partID == part0 ? part1 : part0;
                        const whfc::Node terminal = partID == part0 ? result.source : result.target;

                        for (HyperedgeID e : hg.hyperedgesOf(u)) {
                            if (visitedHyperedge.set(e)) {
                                bool connectToSource = false;
                                bool connectToTarget = false;
                                if (partition.pinsInPart(otherPartID, e) == 0 && partition.pinsInPart(partID, e) > 1) {
                                    local_builder.startHyperedge(hg.hyperedgeWeight(e));
                                    for (NodeID v : hg.pinsOf(e)) {
                                        if (partition[v] == partID) {
                                            if (visitedNode.isSet(v)) {
                                                local_builder.addPin(globalToLocalID[v]);
                                            } else {
                                                if (terminal == result.source) {
                                                    connectToSource = true;
                                                } else {
                                                    connectToTarget = true;
                                                }
                                            }
                                        }
                                    }
                                } else if (partition.pinsInPart(part0, e) > 0 && partition.pinsInPart(part1, e) > 0) {
                                    // This is a cut hyperedge

                                    cutAtStake += hg.hyperedgeWeight(e);
                                    local_builder.startHyperedge(hg.hyperedgeWeight(e));

                                    for (NodeID v : hg.pinsOf(e)) {
                                        if (visitedNode.isSet(v)) {
                                            assert(globalToLocalID[v] < local_builder.numNodes());
                                            local_builder.addPin(globalToLocalID[v]);
                                        } else {
                                            connectToSource |= (partition[v] == part0);
                                            connectToTarget |= (partition[v] == part1);
                                            if (connectToSource && connectToTarget) {
                                                break;
                                            }
                                        }
                                    }
                                }

                                if (connectToSource && connectToTarget) {
                                    local_builder.removeCurrentHyperedge();
                                    baseCut += hg.hyperedgeWeight(e);
                                } else {
                                    if (connectToSource) {
                                        local_builder.addPin(result.source);
                                    }
                                    if (connectToTarget) {
                                        local_builder.addPin(result.target);
                                    }
                                }
                            }

                        }
                    }
                });
            });
            //timer.stop("Add_hyperedges");
        }

        void initialize(uint numNodes, uint numHyperedges) {
            fhgb.clear();
            mock_builder.clear();
            visitedNode.reset();
            visitedHyperedge.reset();
            result = {whfc::Node(0), whfc::Node(0), 0, 0};
            for (MockBuilder& mockBuilder : mockBuilder_thread_specific) {
                mockBuilder.clear();
            }

            for (LayeredQueue<NodeWithDistance>& queue : queue_thread_specific_1) {
                queue.clear();
            }
            for (LayeredQueue<NodeWithDistance>& queue : queue_thread_specific_2) {
                queue.clear();
            }
        }
    };
}