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
        static constexpr bool processCutHyperedgesInParallel = true;

        FlowHypergraphBuilderExtractor(const size_t maxNumNodes, const size_t maxNumEdges, const size_t maxNumPins, int seed, const PartitionConfig& config) :
                fhgb(2 * config.percentage_bfs_from_cut * maxNumNodes + 2, maxNumEdges),
                visitedNode(maxNumNodes), visitedHyperedge(maxNumEdges),
                globalToLocalID(maxNumNodes),
                localToGlobalID(maxNumNodes + 2),
                mt(seed), config(config),
                thisLayer_thread_specific(new tbb::enumerable_thread_specific<std::vector<whfc::Node>>()),
                nextLayer_thread_specific(new tbb::enumerable_thread_specific<std::vector<whfc::Node>>()) {}

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

            if constexpr (doBothBFSInParallel) {
                /*
                timer.start("Parallel_BFS", "BFS");

                tbb::parallel_invoke([&]() {
                    fhgb.addNode(whfc::NodeWeight(0));
                    layered_queue.push(invalid_node);
                    layered_queue.reinitialize();
                    w0 = BreadthFirstSearch(hg, cut_hes, partition, part0, part1, maxW0, result.source, -delta, distanceFromCut, timer, fhgb, layered_queue);
                }, [&]() {
                    whfc::DistanceFromCut temp_dist(hg.numNodes());
                    mock_builder.addNode(whfc::NodeWeight(0));
                    second_queue.push(invalid_node);
                    second_queue.reinitialize();
                    w1 = BreadthFirstSearch(hg, cut_hes, partition, part1, part0, maxW1, whfc::Node(0), delta, temp_dist, timer, mock_builder, second_queue);
                });

                timer.stop("Parallel_BFS");

                timer.start("Merging", "BFS");
                size_t num_nodes_first_search = fhgb.numNodes();
                result.target = whfc::Node(fhgb.numNodes());

                fhgb.addMockBuilder(mock_builder, true);

                for (NodeID u : second_queue.allElements()) {
                    layered_queue.push(u);
                    if (u != whfc::Node::InvalidValue) {
                        globalToLocalID[u] = globalToLocalID[u] + whfc::Node(num_nodes_first_search);
                    }
                }
                timer.stop("Merging");
                 */
            } else {
                localToGlobalID[0] = invalid_node;
                fhgb.addNode(whfc::NodeWeight(0));
                w0 = BreadthFirstSearch(hg, cut_hes, partition, part0, part1, maxW0, result.source, -delta, distanceFromCut, timer, fhgb);

                result.target = whfc::Node(fhgb.numNodes());
                fhgb.addNode(whfc::NodeWeight(0));
                localToGlobalID[result.target] = invalid_node;
                w1 = BreadthFirstSearch(hg, cut_hes, partition, part1, part0, maxW1, result.target, delta, distanceFromCut, timer, fhgb);
            }


            timer.stop("BFS");

            timer.start("Process Cut Hyperedges", "Extraction");
            tbb::blocked_range<size_t> range(0, cut_hes.size(), 1000);

            if (processCutHyperedgesInParallel && cut_hes.size() > 1000) {
                // TODO: fix mockBuilders to fhgb and don't reinitialize every time
                mockBuilder_thread_specific = tbb::enumerable_thread_specific<MockBuilder>(std::ref(fhgb.getNodes()));

                tbb::parallel_for(range, [=](const tbb::blocked_range<size_t>& sub_range) {
                    MockBuilder& builder = mockBuilder_thread_specific.local();
                    processCutHyperedges(hg, cut_hes, partition, part0, part1, builder, sub_range);
                });

                fhgb.addMockBuildersParallel(mockBuilder_thread_specific);
            } else {
                processCutHyperedges(hg, cut_hes, partition, part0, part1, fhgb, range);
            }

            for (whfc::Flow& baseCut : baseCut_thread_specific) {
                result.baseCut += baseCut;
                baseCut = 0;
            }

            for (whfc::Flow& cutAtStake : cutAtStake_thread_specific) {
                result.cutAtStake += cutAtStake;
                cutAtStake = 0;
            }

            timer.stop("Process Cut Hyperedges");

            std::vector<NodeWeight> totalWeights = partition.partitionWeights();    // REVIEW NOTE where do you need this?

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
        ldc::AtomicTimestampSet<> visitedNode;
        ldc::AtomicTimestampCounter<> visitedHyperedge;
        std::vector<whfc::Node> globalToLocalID;
        std::vector<NodeID> localToGlobalID;
        ExtractorInfo result;
        std::mt19937 mt;
        const PartitionConfig& config;

        MockBuilder mock_builder;
        tbb::enumerable_thread_specific<MockBuilder> mockBuilder_thread_specific;
        tbb::enumerable_thread_specific<whfc::Flow> baseCut_thread_specific;
        tbb::enumerable_thread_specific<whfc::Flow> cutAtStake_thread_specific;
        tbb::enumerable_thread_specific<whfc::NodeWeight> weights_thread_specific;

        std::unique_ptr<tbb::enumerable_thread_specific<std::vector<whfc::Node>>> thisLayer_thread_specific;
        std::unique_ptr<tbb::enumerable_thread_specific<std::vector<whfc::Node>>> nextLayer_thread_specific;

        inline bool tryToVisitNode(const NodeID u, CSRHypergraph &hg, std::atomic<whfc::NodeWeight> &w, std::vector<whfc::Node>& vectorToPush,
                const whfc::NodeWeight& maxWeight, whfc::NodeWeight& localWeight) {
            if (visitedNode.set(u)) {
                if (w.fetch_add(hg.nodeWeight(u)) + hg.nodeWeight(u) <= maxWeight) {
                    vectorToPush.push_back(static_cast<whfc::Node>(u));
                    localWeight += hg.nodeWeight(u);        // REVIEW NOTE why have a thread local weight if you also have an atomic?
                    return true;
                } else {
                    w.fetch_sub(hg.nodeWeight(u));
                    visitedNode.reset(u);   // REVIEW NOTE is this a good idea? this means other threads will try to add u again --> more contention on w
                                            // --> necessary when adding pins
                }
            }
            return false;
        }

        template<class Builder>
        bool writeNodeLayerToBuilder(Builder& builder, tbb::enumerable_thread_specific<std::vector<whfc::Node>>& layer,
                CSRHypergraph& hg, whfc::DistanceFromCut& distanceFromCut, whfc::HopDistance d) {
            std::atomic<size_t> node_counter = builder.numNodes();
            size_t num_nodes = builder.numNodes();
            std::vector<whfc::FlowHypergraph::NodeData>& nodes = builder.getNodes();

            for (std::vector<whfc::Node>& localLayer : layer) {
                num_nodes += localLayer.size();
            }

            bool has_nodes = num_nodes > builder.numNodes();

            // Write node data into builder
            nodes.resize(num_nodes + 1); // maybe use reserve

            tbb::parallel_for_each(layer.range(), [&](const std::vector<whfc::Node>& vector) {
                tbb::parallel_for(tbb::blocked_range<size_t>(0, vector.size()), [&](const tbb::blocked_range<size_t>& indices) {
                    size_t begin_index = node_counter.fetch_add(indices.size());
                    for (size_t i = indices.begin(); i < indices.end(); ++i) {
                        const whfc::Node localID(begin_index);
                        begin_index++;
                        assert(localID < num_nodes);
                        nodes[localID].weight = whfc::NodeWeight(hg.nodeWeight(vector[i]));
                        globalToLocalID[vector[i]] = localID;
                        localToGlobalID[localID] = vector[i];
                        distanceFromCut[localID] = d;
                    }
                });
            });
            assert(num_nodes == node_counter.load());

            return has_nodes;
        }

        template<class PartitionImpl, typename CutEdgeRange, typename Builder>
        whfc::NodeWeight BreadthFirstSearch(CSRHypergraph &hg, CutEdgeRange &cut_hes, const PartitionImpl &partition,
                                            PartitionBase::PartitionID partID, PartitionBase::PartitionID otherPartID,
                                            NodeWeight maxWeight, whfc::Node terminal, whfc::HopDistance delta,
                                            whfc::DistanceFromCut& distanceFromCut, whfc::TimeReporter& timer, Builder& builder) {
            std::atomic<whfc::NodeWeight> w = 0;
            weights_thread_specific.clear();
            whfc::HopDistance d = delta;

            std::vector<whfc::FlowHypergraph::NodeData>& nodes = builder.getNodes();
            mockBuilder_thread_specific = tbb::enumerable_thread_specific<MockBuilder>(std::ref(nodes));

            thisLayer_thread_specific->clear();
            nextLayer_thread_specific->clear();

            visitedHyperedge.nextRound();

            bool nodes_left = false;
            whfc::NodeWeight& localWeight = weights_thread_specific.local();

            // Collect boundary vertices
            timer.start("Collect_boundary_vertices", "BFS");
            tbb::blocked_range<size_t> range(0, cut_hes.size(), 1000);
            tbb::parallel_for(range, [&](const tbb::blocked_range<size_t>& sub_range) {
                auto& thisLayer = thisLayer_thread_specific->local();
                for (size_t i = sub_range.begin(); i < sub_range.end(); ++i) {
                    const HyperedgeID e = cut_hes[i];
                    for (NodeID v : hg.pinsOf(e)) {
                        if (partition[v] == partID) {
                            tryToVisitNode(v, hg, w, thisLayer, maxWeight, localWeight);
                        }
                    }
                }
            });
            timer.stop("Collect_boundary_vertices");

            // Write nodes to the builder, this also ensures globalToLocal-mapping is set correctly
            timer.start("Write_boundary_to_builder", "BFS");
            nodes_left = writeNodeLayerToBuilder(builder, *thisLayer_thread_specific, hg, distanceFromCut, d);
            d += delta;
            timer.stop("Write_boundary_to_builder");

            size_t counterValue = visitedHyperedge.getCounter();

            while (nodes_left) {
                timer.start("Scan_hyperedges_and_add_nodes", "BFS");
                // scan hyperedges and add nodes to the next layer
                tbb::parallel_for(thisLayer_thread_specific->range(), [&](const tbb::blocked_range<tbb::enumerable_thread_specific<std::vector<whfc::Node>>::iterator>& range) {
                    for (auto& vector : range) {
                        tbb::parallel_for(tbb::blocked_range<size_t>(0, vector.size()), [&](const tbb::blocked_range<size_t>& indices) {
                            auto& nextLayer = nextLayer_thread_specific->local();
                            whfc::NodeWeight& localWeight = weights_thread_specific.local();

                            for (size_t i = indices.begin(); i < indices.end(); ++i) {
                                const NodeID u = vector[i];

                                for (HyperedgeID e : hg.hyperedgesOf(u)) {
                                    if (visitedHyperedge.set(e) == 0 && partition.pinsInPart(otherPartID, e) == 0 && partition.pinsInPart(partID, e) > 1) {
                                        for (NodeID v : hg.pinsOf(e)) {
                                            if (partition[v] == partID) {
                                                tryToVisitNode(v, hg, w, nextLayer, maxWeight, localWeight);
                                            }
                                        }
                                    }
                                }
                            }
                        });
                    }
                });
                timer.stop("Scan_hyperedges_and_add_nodes");

                // This ensures a correct globalToLocal-mapping for the next step
                timer.start("Write_layer_to_builder", "BFS");
                nodes_left = writeNodeLayerToBuilder(builder, *nextLayer_thread_specific, hg, distanceFromCut, d);
                timer.stop("Write_layer_to_builder");

                // REVIEW NOTE why do you write the pins to the local builders after every layer
                // I thought we should collect all nodes first --> more work volume in just one call gives better speedups
                // this solves the hyperedge marker issue
                // if we know the number of pins we can also get rid of the mock builders --> sum up pin count in part in the beginning and do resize if necessary?

                visitedHyperedge.nextRound();

                // scan hyperedges again and write them to local builders
                timer.start("Scan_hyperedges_and_add_hyperedges", "BFS");
                tbb::parallel_for(thisLayer_thread_specific->range(), [&](const tbb::blocked_range<tbb::enumerable_thread_specific<std::vector<whfc::Node>>::iterator>& range) {
                    for (auto& vector : range) {
                        tbb::parallel_for(tbb::blocked_range<size_t>(0, vector.size()), [&](const tbb::blocked_range<size_t>& indices) {
                            auto& local_builder = mockBuilder_thread_specific.local();

                            for (size_t i = indices.begin(); i < indices.end(); ++i) {
                                const NodeID u = vector[i];

                                for (HyperedgeID e : hg.hyperedgesOf(u)) {
                                    if (visitedHyperedge.set(e) == counterValue && partition.pinsInPart(otherPartID, e) == 0 &&
                                        partition.pinsInPart(partID, e) > 1) {
                                        local_builder.startHyperedge(hg.hyperedgeWeight(e));
                                        bool connectToTerminal = false;
                                        for (NodeID v : hg.pinsOf(e)) {
                                            if (partition[v] == partID) {
                                                if (visitedNode.isSet(v)) {
                                                    local_builder.addPin(globalToLocalID[v]);
                                                } else {
                                                    connectToTerminal = true;
                                                }
                                            }
                                        }
                                        if (connectToTerminal) {
                                            local_builder.addPin(terminal);
                                        }
                                    }
                                }
                            }
                        });
                    }
                });
                timer.stop("Scan_hyperedges_and_add_hyperedges");

                std::swap(thisLayer_thread_specific, nextLayer_thread_specific);

                for (auto& vector : *nextLayer_thread_specific) {
                    vector.clear();
                }

                d += delta;

            }

            distanceFromCut[terminal] = d;

            timer.start("addMockBuildersParallel", "BFS");
            builder.addMockBuildersParallel(mockBuilder_thread_specific);       // REVIEW NOTE this guy doesn't clear the mock builder --> adds earlier layers multiple times?
            timer.stop("addMockBuildersParallel");

            whfc::NodeWeight totalWeight = 0;
            for (NodeWeight weight : weights_thread_specific) {
                totalWeight += weight;
            }

            return totalWeight;
        }

        template<class PartitionImpl, class CutEdgeRange, class Builder>
        void processCutHyperedges(const CSRHypergraph &hg, CutEdgeRange &cut_hes, const PartitionImpl &partition,
                                  const PartitionBase::PartitionID part0, const PartitionBase::PartitionID part1,
                                  Builder& builder, const tbb::blocked_range<size_t>& range) {
            whfc::Flow& baseCut = baseCut_thread_specific.local();
            whfc::Flow& cutAtStake = cutAtStake_thread_specific.local();

            for (size_t i = range.begin(); i < range.end(); ++i) {
                const HyperedgeID e = cut_hes[i];
                //assert(!visitedHyperedge.isSet(e));
                assert(partition.pinsInPart(part0, e) > 0 && partition.pinsInPart(part1, e) > 0);
                bool connectToSource = false;
                bool connectToTarget = false;
                cutAtStake += hg.hyperedgeWeight(e);
                builder.startHyperedge(hg.hyperedgeWeight(e));

                for (NodeID v : hg.pinsOf(e)) {
                    if (visitedNode.isSet(v)) {
                        assert(globalToLocalID[v] < builder.numNodes());
                        builder.addPin(globalToLocalID[v]);
                    } else {
                        connectToSource |= (partition[v] == part0);
                        connectToTarget |= (partition[v] == part1);
                        if (connectToSource && connectToTarget) {
                            break;
                        }
                    }
                }
                if (connectToSource && connectToTarget) {
                    builder.removeCurrentHyperedge();
                    baseCut += hg.hyperedgeWeight(e);
                } else {
                    if (connectToSource) {
                        builder.addPin(result.source);
                    }
                    if (connectToTarget) {
                        builder.addPin(result.target);
                    }
                }
            }
        }

        void initialize(uint numNodes, uint numHyperedges) {
            fhgb.clear();
            mock_builder.clear();
            visitedNode.reset();
            visitedHyperedge.reset();
            result = {whfc::Node(0), whfc::Node(0), 0, 0};
        }
    };
}