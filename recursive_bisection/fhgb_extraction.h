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

            tryVisitNodeLock = std::vector<std::mutex>(partition.getGraph().numNodes()); //TODO: this has to be replaced

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

            localToGlobalID[0] = invalid_node;
            fhgb.addNode(whfc::NodeWeight(0));
            result.target = whfc::Node(1);
            fhgb.addNode(whfc::NodeWeight(0));
            localToGlobalID[result.target] = invalid_node;

            std::atomic<NodeWeight> w0 = 0;
            std::atomic<NodeWeight> w1 = 0;
            // This fills the first layers for the bfs
            processCutHyperedges(hg, cut_hes, partition, part0, part1, maxW0, maxW1, delta, w0, w1);

            timer.start("Scan_hyperedges_and_add_nodes", "BFS");
            tbb::parallel_invoke([&]() {
                BreadthFirstSearch(hg, cut_hes, partition, part0, part1, maxW0, result.source, -delta, distanceFromCut, timer, queue_thread_specific_1, w0);
            }, [&]() {
                BreadthFirstSearch(hg, cut_hes, partition, part1, part0, maxW1, result.target, delta, distanceFromCut, timer, queue_thread_specific_2, w1);
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

            std::vector<size_t> prefix_sum(hyperedges.size() + 1);

            computePrefixSum(hyperedges, prefix_sum);

            timer.start("Add_hyperedges", "BFS");
            addHyperedges(hg, partition, part0, part1, timer, hyperedges, prefix_sum);
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

        std::vector<std::mutex> tryVisitNodeLock;

        void computePrefixSum(const std::vector<HyperedgeWithSize>& vector, std::vector<size_t>& result) {
            tbb::parallel_scan(tbb::blocked_range<size_t>(0, vector.size()), 0,
                [&](const tbb::blocked_range<size_t>& r, size_t sum, bool is_final_scan) -> size_t {
                size_t temp = sum;
                for (size_t i = r.begin(); i < r.end(); ++i) {
                    temp = temp + vector[i].pin_count;
                    if(is_final_scan)
                        result[i + 1] = temp;
                }
                return temp;
            }, [] (size_t left, size_t right) {
                return left + right;
            });
        }


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

        inline bool visitNode(const CSRHypergraph& hg, const NodeID u,
                std::atomic<whfc::NodeWeight>& w, NodeWeight& lastSeenValue, const NodeWeight maxWeight, const whfc::HopDistance d,
                LayeredQueue<NodeWithDistance>& queue) {
            std::lock_guard<std::mutex> guard(tryVisitNodeLock[u]);
            lastSeenValue = w.load(std::memory_order_relaxed);
            if (lastSeenValue + hg.nodeWeight(u) <= maxWeight) {
                if (visitedNode.set(u)) {
                    w += hg.nodeWeight(u);
                    queue.push({static_cast<whfc::Node>(u), d});
                }
                return true;
            }
            return visitedNode.isSet(u);
        }

        template<class PartitionImpl, typename CutEdgeRange>
        void processCutHyperedges(CSRHypergraph& hg, CutEdgeRange &cut_hes, const PartitionImpl& partition,
                PartitionBase::PartitionID partID, PartitionBase::PartitionID otherPartID,
                NodeWeight maxWeight1, NodeWeight maxWeight2, whfc::HopDistance delta, std::atomic<NodeWeight>& w1,
                std::atomic<NodeWeight>& w2) {
            tbb::blocked_range<size_t> range(0, cut_hes.size());
            tbb::parallel_for(range, [&](const tbb::blocked_range<size_t>& sub_range) {
                NodeWeight lastSeenValue1 = 0;
                NodeWeight lastSeenValue2 = 0;
                LayeredQueue<NodeWithDistance>& queue1 = queue_thread_specific_1.local();
                LayeredQueue<NodeWithDistance>& queue2 = queue_thread_specific_2.local();
                std::vector<HyperedgeWithSize>& hyperedges = hyperedges_thread_specific.local();
                whfc::Flow &baseCut = baseCut_thread_specific.local();
                whfc::Flow &cutAtStake = cutAtStake_thread_specific.local();
                for (size_t i = sub_range.begin(); i < sub_range.end(); ++i) {
                    const HyperedgeID e = cut_hes[i];
                    cutAtStake += hg.hyperedgeWeight(e);
                    size_t num_pins = 0;
                    bool connectToSource = false;
                    bool connectToTarget = false;
                    for (NodeID u : hg.pinsOf(e)) {
                        if (partition[u] == partID) {
                            if (visitNode(hg, u, w1, lastSeenValue1, maxWeight1, delta, queue1)) {
                                num_pins++;
                            } else {
                                connectToSource = true;
                            }
                        } else if (partition[u] == otherPartID) {
                            if (visitNode(hg, u, w2, lastSeenValue2, maxWeight2, -delta, queue2)) {
                                num_pins++;
                            } else {
                                connectToTarget = true;
                            }
                        }
                    }
                    if (!(connectToSource && connectToTarget)) {
                        if (connectToSource || connectToTarget) num_pins++;
                        hyperedges.push_back({e, num_pins});
                    } else {
                        baseCut += hg.hyperedgeWeight(e);
                    }
                    //if (lastSeenValue1 >= maxWeight1 && lastSeenValue2 >= maxWeight2) break;
                }
            });

        }

        template<class PartitionImpl, typename CutEdgeRange>
        void BreadthFirstSearch(CSRHypergraph &hg, CutEdgeRange &cut_hes, const PartitionImpl &partition,
                                            PartitionBase::PartitionID partID, PartitionBase::PartitionID otherPartID,
                                            NodeWeight maxWeight, whfc::Node terminal, whfc::HopDistance delta,
                                            whfc::DistanceFromCut& distanceFromCut, whfc::TimeReporter& timer,
                                            tbb::enumerable_thread_specific<LayeredQueue<NodeWithDistance>>& queue_thread_specific,
                                            std::atomic<NodeWeight>& w) {
            whfc::HopDistance d = delta;

            d += delta;

            bool nodes_left = false;

            for (LayeredQueue<NodeWithDistance>& queue : queue_thread_specific) {
                queue.finishNextLayer();
                if (!queue.currentLayerEmpty()) {
                    nodes_left = true;
                }
            }

            while (nodes_left) {
                // scan hyperedges and add nodes to the next layer
                tbb::parallel_for_each(queue_thread_specific, [&](LayeredQueue<NodeWithDistance>& queue) {
                    tbb::parallel_for(tbb::blocked_range<size_t>(queue.currentLayerStart(), queue.currentLayerEnd()), [&](const tbb::blocked_range<size_t>& indices) {
                        LayeredQueue<NodeWithDistance>& local_queue = queue_thread_specific.local();
                        std::vector<HyperedgeWithSize>& hyperedges = hyperedges_thread_specific.local();
                        NodeWeight lastSeenValue = 0;

                        for (size_t i = indices.begin(); i < indices.end(); ++i) {
                            const NodeID u = queue.elementAt(i).node;

                            //if (lastSeenValue >= maxWeight) break;
                            for (HyperedgeID e : hg.hyperedgesOf(u)) {
                                if (partition.pinsInPart(otherPartID, e) == 0 && partition.pinsInPart(partID, e) > 1 && visitedHyperedge.set(e)) {
                                    size_t num_pins = 0;
                                    bool connectToTerminal = false;
                                    for (NodeID v : hg.pinsOf(e)) {
                                        if (partition[v] == partID) {
                                            if (visitNode(hg, v, w, lastSeenValue, maxWeight, d, local_queue)) {
                                                num_pins++;
                                            } else {
                                                connectToTerminal = true;
                                            }
                                        }
                                    }
                                    if (connectToTerminal) num_pins++;
                                    hyperedges.push_back({e, num_pins});
                                }
                            }
                        }
                    });
                });

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
        }

        using PinIndexRange = mutable_index_range<whfc::PinIndex>;

        template<typename PartitionImpl>
        void addHyperedges(CSRHypergraph& hg, const PartitionImpl &partition, PartitionBase::PartitionID part0,
                PartitionBase::PartitionID part1, whfc::TimeReporter& timer, std::vector<HyperedgeWithSize>& hyperedges,
                const std::vector<size_t>& prefix_sum) {

            fhgb.hyperedges.resize(hyperedges.size() + 1);
            fhgb.pins_sending_flow.resize(hyperedges.size());
            fhgb.pins_receiving_flow.resize(hyperedges.size());
            fhgb.pins.resize(prefix_sum[prefix_sum.size() - 1]);

            fhgb.hyperedges.back().first_out = whfc::PinIndex(prefix_sum[prefix_sum.size() - 1]);
            fhgb.numPinsAtHyperedgeStart = prefix_sum[prefix_sum.size() - 1];

            tbb::parallel_for(tbb::blocked_range<size_t>(0, hyperedges.size()), [&](const tbb::blocked_range<size_t>& indices) {


                for (size_t i = indices.begin(); i < indices.end(); ++i) {
                    const HyperedgeWithSize edge = hyperedges[i];
                    const HyperedgeID e = edge.e;
                    bool connectToSource = false;
                    bool connectToTarget = false;
                    size_t nextPinPosition = prefix_sum[i];

                    fhgb.pins_sending_flow[i] = PinIndexRange(whfc::PinIndex(prefix_sum[i]), whfc::PinIndex(prefix_sum[i]));
                    fhgb.hyperedges[i] = {whfc::PinIndex(prefix_sum[i]), whfc::Flow(0), whfc::Flow(hg.hyperedgeWeight(e))};
                    fhgb.pins_receiving_flow[i] = PinIndexRange(whfc::PinIndex(prefix_sum[i+1]), whfc::PinIndex(prefix_sum[i+1]));

                    if (partition.pinsInPart(part0, e) > 0 && partition.pinsInPart(part1, e) > 0) {
                        // Cut hyperedge
                        for (NodeID v : hg.pinsOf(e)) {
                            if (visitedNode.isSet(v)) {
                                fhgb.pins[nextPinPosition++] = {globalToLocalID[v], whfc::InHeIndex::Invalid()};
                                __sync_fetch_and_add(&fhgb.nodes[globalToLocalID[v]+1].first_out.value(), 1);
                            } else {
                                connectToSource |= (partition[v] == part0);
                                connectToTarget |= (partition[v] == part1);
                                assert(!(connectToSource && connectToTarget));
                            }
                        }
                    } else {
                        const PartitionBase::PartitionID partID = partition.pinsInPart(part0, e) > 0 ? part0 : part1;
                        const whfc::Node terminal = partID == part0 ? result.source : result.target;
                        for (NodeID v : hg.pinsOf(e)) {
                            if (partition[v] == partID) {
                                if (visitedNode.isSet(v)) {
                                    fhgb.pins[nextPinPosition++] = {globalToLocalID[v], whfc::InHeIndex::Invalid()};
                                    __sync_fetch_and_add(&fhgb.nodes[globalToLocalID[v]+1].first_out.value(), 1);
                                } else {
                                    if (terminal == result.source) {
                                        connectToSource = true;
                                    } else {
                                        connectToTarget = true;
                                    }
                                }
                            }
                        }
                    }

                    if (connectToSource) {
                        fhgb.pins[nextPinPosition++] = {result.source, whfc::InHeIndex::Invalid()};
                        __sync_fetch_and_add(&fhgb.nodes[result.source+1].first_out.value(), 1);
                    }
                    if (connectToTarget) {
                        fhgb.pins[nextPinPosition++] = {result.target, whfc::InHeIndex::Invalid()};
                        __sync_fetch_and_add(&fhgb.nodes[result.target+1].first_out.value(), 1);
                    }
                }
            });
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