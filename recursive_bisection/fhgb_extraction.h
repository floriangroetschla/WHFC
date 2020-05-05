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

                addHyperedges(fhgb, hg, partition, part0, part1, timer);
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
        tbb::enumerable_thread_specific<whfc::NodeWeight> weights_thread_specific;

        std::unique_ptr<tbb::enumerable_thread_specific<std::vector<whfc::Node>>> thisLayer_thread_specific;
        std::unique_ptr<tbb::enumerable_thread_specific<std::vector<whfc::Node>>> nextLayer_thread_specific;

        // returns last value of w
        inline NodeWeight tryToVisitNode(const NodeID u, CSRHypergraph &hg, std::atomic<whfc::NodeWeight> &w, std::vector<whfc::Node>& vectorToPush,
                const whfc::NodeWeight& maxWeight, whfc::NodeWeight& localWeight) {
            NodeWeight lastSeenValue = 0;
            if (visitedNode.set(u)) {
                lastSeenValue = w.fetch_add(hg.nodeWeight(u));
                if (lastSeenValue + hg.nodeWeight(u) <= maxWeight) {
                    vectorToPush.push_back(static_cast<whfc::Node>(u));
                    localWeight += hg.nodeWeight(u);        // REVIEW NOTE why have a thread local weight if you also have an atomic?
                } else {
                    w.fetch_sub(hg.nodeWeight(u), std::memory_order_relaxed);
                    visitedNode.reset(u);   // REVIEW NOTE is this a good idea? this means other threads will try to add u again --> more contention on w
                                            // --> necessary when adding pins
                }
            }
            return lastSeenValue;
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

            tbb::parallel_for_each(layer, [&](const std::vector<whfc::Node>& vector) {
                size_t begin_index = node_counter.fetch_add(vector.size());
                tbb::parallel_for(tbb::blocked_range<size_t>(0, vector.size()), [&](const tbb::blocked_range<size_t>& indices) {
                    for (size_t i = indices.begin(); i < indices.end(); ++i) {
                        const whfc::Node localID(begin_index + i);
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

            thisLayer_thread_specific->clear();
            nextLayer_thread_specific->clear();

            bool nodes_left = false;
            whfc::NodeWeight& localWeight = weights_thread_specific.local();

            // Collect boundary vertices
            timer.start("Collect_boundary_vertices", "BFS");
            tbb::blocked_range<size_t> range(0, cut_hes.size(), 1000);
            tbb::parallel_for(range, [&](const tbb::blocked_range<size_t>& sub_range) {
                NodeWeight lastSeenValue = 0;
                auto& thisLayer = thisLayer_thread_specific->local();
                for (size_t i = sub_range.begin(); i < sub_range.end(); ++i) {
                    const HyperedgeID e = cut_hes[i];
                    for (NodeID v : hg.pinsOf(e)) {
                        if (partition[v] == partID && lastSeenValue + hg.nodeWeight(v) <= maxWeight) {
                            lastSeenValue = tryToVisitNode(v, hg, w, thisLayer, maxWeight, localWeight);
                        }
                    }
                    if (lastSeenValue == maxWeight) break;
                }
            });
            timer.stop("Collect_boundary_vertices");

            // Write nodes to the builder, this also ensures globalToLocal-mapping is set correctly
            timer.start("Write_boundary_to_builder", "BFS");
            nodes_left = writeNodeLayerToBuilder(builder, *thisLayer_thread_specific, hg, distanceFromCut, d);
            d += delta;
            timer.stop("Write_boundary_to_builder");


            while (nodes_left) {
                timer.start("Scan_hyperedges_and_add_nodes", "BFS");
                // scan hyperedges and add nodes to the next layer
                tbb::parallel_for_each(*thisLayer_thread_specific, [&](const std::vector<whfc::Node>& vector) {
                    tbb::parallel_for(tbb::blocked_range<size_t>(0, vector.size(), 1000), [&](const tbb::blocked_range<size_t>& indices) {
                        auto& nextLayer = nextLayer_thread_specific->local();
                        whfc::NodeWeight& localWeight = weights_thread_specific.local();
                        NodeWeight lastSeenValue = 0;

                        for (size_t i = indices.begin(); i < indices.end(); ++i) {
                            const NodeID u = vector[i];

                            for (HyperedgeID e : hg.hyperedgesOf(u)) {
                                if (partition.pinsInPart(otherPartID, e) == 0 && partition.pinsInPart(partID, e) > 1 && visitedHyperedge.set(e)) {
                                    for (NodeID v : hg.pinsOf(e)) {
                                        if (partition[v] == partID && lastSeenValue + hg.nodeWeight(v) <= maxWeight) {
                                            lastSeenValue = tryToVisitNode(v, hg, w, nextLayer, maxWeight, localWeight);
                                        }
                                    }
                                }
                                if (lastSeenValue == maxWeight) break;
                            }
                        }
                    });
                });
                timer.stop("Scan_hyperedges_and_add_nodes");

                // This ensures a correct globalToLocal-mapping for the next step
                timer.start("Write_layer_to_builder", "BFS");
                nodes_left = writeNodeLayerToBuilder(builder, *nextLayer_thread_specific, hg, distanceFromCut, d);
                timer.stop("Write_layer_to_builder");


                std::swap(thisLayer_thread_specific, nextLayer_thread_specific);

                for (auto& vector : *nextLayer_thread_specific) {
                    vector.clear();
                }

                d += delta;

            }

            distanceFromCut[terminal] = d;

            whfc::NodeWeight totalWeight = 0;
            for (NodeWeight weight : weights_thread_specific) {
                totalWeight += weight;
            }

            return totalWeight;
        }

        template<typename Builder, typename PartitionImpl>
        void addHyperedges(Builder& builder, CSRHypergraph& hg, const PartitionImpl &partition, PartitionBase::PartitionID part0, PartitionBase::PartitionID part1, whfc::TimeReporter& timer) {
            visitedHyperedge.reset();

            timer.start("Add_hyperedges", "BFS");
            tbb::parallel_for(tbb::blocked_range<size_t>(0, builder.numNodes()), [&](const tbb::blocked_range<size_t>& indices) {
                auto& local_builder = mockBuilder_thread_specific.local();
                local_builder.setSource(result.source);
                local_builder.setTarget(result.target);
                whfc::Flow& baseCut = baseCut_thread_specific.local();
                whfc::Flow& cutAtStake = cutAtStake_thread_specific.local();

                for (size_t i = indices.begin(); i < indices.end(); ++i) {
                    const NodeID u = localToGlobalID[i];

                    // skip terminals
                    if (u == invalid_node) {
                        continue;
                    }

                    const PartitionBase::PartitionID partID = i < result.target ? part0 : part1;
                    const PartitionBase::PartitionID otherPartID = i < result.target ? part1 : part0;
                    const whfc::Node terminal = i < result.target ? result.source : result.target;

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
            timer.stop("Add_hyperedges");

            timer.start("addMockBuildersParallel", "BFS");
            builder.addMockBuildersParallel(mockBuilder_thread_specific);
            timer.stop("addMockBuildersParallel");

            for (whfc::Flow& baseCut : baseCut_thread_specific) {
                result.baseCut += baseCut;
                baseCut = 0;
            }

            for (whfc::Flow& cutAtStake : cutAtStake_thread_specific) {
                result.cutAtStake += cutAtStake;
                cutAtStake = 0;
            }
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
        }
    };
}