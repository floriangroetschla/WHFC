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

        FlowHypergraphBuilderExtractor(const size_t maxNumNodes, const size_t maxNumEdges, const size_t maxNumPins, int seed, const PartitionConfig& config) :
                fhgb(2 * config.percentage_bfs_from_cut * maxNumNodes + 2, maxNumEdges),
                layered_queue(maxNumNodes + 2),
                second_queue(maxNumNodes + 2),
                visitedNode(maxNumNodes), visitedHyperedge(maxNumEdges),
                globalToLocalID(maxNumNodes),
                mt(seed), config(config) {}

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

            assert(layered_queue.empty());

            timer.start("BFS", "Extraction");

            timer.start("Parallel_BFS", "BFS");

            whfc::NodeWeight w0, w1;

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

            timer.stop("BFS");

            timer.start("Process Cut Hyperedges", "Extraction");
            tbb::blocked_range<size_t> range(0, cut_hes.size(), 1000);

            if (cut_hes.size() > 1000) {
                mockBuilder_thread_specific = tbb::enumerable_thread_specific<MockBuilder>(std::ref(fhgb.getNodes()));

                tbb::parallel_for(range, [=](const tbb::blocked_range<size_t>& sub_range) {
                    MockBuilder& builder = mockBuilder_thread_specific.local();
                    processCutHyperedges(hg, cut_hes, partition, part0, part1, builder, sub_range);
                });

                for (MockBuilder& builder : mockBuilder_thread_specific) {
                    fhgb.addMockBuilder(builder, false);
                }
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

            std::vector<NodeWeight> totalWeights = partition.partitionWeights();

            fhgb.nodeWeight(result.source) = partition.partWeight(part0) - w0;
            fhgb.nodeWeight(result.target) = partition.partWeight(part1) - w1;

            timer.start("Finalize", "Extraction");
            fhgb.finalize();
            timer.stop("Finalize");

            return result;
        }

        auto localNodeIDs() const {
            return boost::irange<whfc::Node>(whfc::Node(0), whfc::Node(layered_queue.queueEnd()));
        }

        whfc::Node global2local(const NodeID x) const {
            assert(visitedNode.contains(x));
            return globalToLocalID[x];
        }

        NodeID local2global(const whfc::Node x) const { return layered_queue.elementAt(x); }

    private:
        LayeredQueue<NodeID> layered_queue;
        LayeredQueue<NodeID> second_queue;
        ldc::TimestampSet<> visitedNode;
        ldc::TimestampSet<> visitedHyperedge;
        std::vector<whfc::Node> globalToLocalID;
        ExtractorInfo result;
        std::mt19937 mt;
        const PartitionConfig& config;

        MockBuilder mock_builder;
        tbb::enumerable_thread_specific<MockBuilder> mockBuilder_thread_specific;
        tbb::enumerable_thread_specific<whfc::Flow> baseCut_thread_specific;
        tbb::enumerable_thread_specific<whfc::Flow> cutAtStake_thread_specific;

        template<typename Builder>
        void visitNode(const NodeID node, CSRHypergraph &hg, whfc::NodeWeight &w, Builder& builder, LayeredQueue<NodeID>& queue) {
            globalToLocalID[node] = whfc::Node(builder.numNodes());
            queue.push(node);
            visitedNode.add(node);
            builder.addNode(whfc::NodeWeight(hg.nodeWeight(node)));
            w += hg.nodeWeight(node);
        }

        template<class PartitionImpl, typename CutEdgeRange, typename Builder>
        whfc::NodeWeight BreadthFirstSearch(CSRHypergraph &hg, CutEdgeRange &cut_hes, const PartitionImpl &partition,
                                            PartitionBase::PartitionID partID, PartitionBase::PartitionID otherPartID,
                                            NodeWeight maxWeight, whfc::Node terminal, whfc::HopDistance delta,
                                            whfc::DistanceFromCut& distanceFromCut, whfc::TimeReporter& timer, Builder& builder,
                                            LayeredQueue<NodeID>& queue) {
            whfc::NodeWeight w = 0;
            whfc::HopDistance d = delta;

            // Collect boundary vertices
            for (const HyperedgeID e : cut_hes) {
                for (NodeID v : hg.pinsOf(e)) {
                    if (!visitedNode.contains(v) && partition[v] == partID && w + hg.nodeWeight(v) <= maxWeight) {
                        visitNode(v, hg, w, builder, queue);
                        distanceFromCut[globalToLocalID[v]] = d;
                    }
                }
            }
            
            while (!queue.empty()) {
                if (queue.currentLayerEmpty()) {
                    queue.finishNextLayer();
                    d += delta;
                }
                NodeID u = queue.pop();

                for (HyperedgeID e : hg.hyperedgesOf(u)) {
                    if (!visitedHyperedge.contains(e) && partition.pinsInPart(otherPartID, e) == 0 &&
                        partition.pinsInPart(partID, e) > 1) {
                        builder.startHyperedge(hg.hyperedgeWeight(e));
                        bool connectToTerminal = false;
                        for (NodeID v : hg.pinsOf(e)) {
                            if (partition[v] == partID) {
                                if (!visitedNode.contains(v) && w + hg.nodeWeight(v) <= maxWeight) {
                                    visitNode(v, hg, w, builder, queue);
                                    distanceFromCut[globalToLocalID[v]] = d;
                                }

                                if (visitedNode.contains(v)) {
                                    assert(globalToLocalID[v] < builder.numNodes());
                                    builder.addPin(globalToLocalID[v]);
                                } else {
                                    connectToTerminal = true;
                                }
                            }
                        }
                        if (connectToTerminal) {
                            assert(terminal < builder.numNodes());
                            builder.addPin(terminal);
                        }
                        visitedHyperedge.add(e);
                    }
                }
            }

            d += delta;
            distanceFromCut[terminal] = d;

            return w;
        }

        template<class PartitionImpl, class CutEdgeRange, class Builder>
        void processCutHyperedges(const CSRHypergraph &hg, CutEdgeRange &cut_hes, const PartitionImpl &partition,
                                  const PartitionBase::PartitionID part0, const PartitionBase::PartitionID part1,
                                  Builder& builder, const tbb::blocked_range<size_t>& range) {
            whfc::Flow& baseCut = baseCut_thread_specific.local();
            whfc::Flow& cutAtStake = cutAtStake_thread_specific.local();

            for (size_t i = range.begin(); i < range.end(); ++i) {
                const HyperedgeID e = cut_hes[i];
                assert(!visitedHyperedge.contains(e));
                assert(partition.pinsInPart(part0, e) > 0 && partition.pinsInPart(part1, e) > 0);
                bool connectToSource = false;
                bool connectToTarget = false;
                cutAtStake += hg.hyperedgeWeight(e);
                builder.startHyperedge(hg.hyperedgeWeight(e));

                for (NodeID v : hg.pinsOf(e)) {
                    if (visitedNode.contains(v)) {
                        assert(globalToLocalID[v] < fhgb.numNodes());
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
            layered_queue.clear();
            second_queue.clear();
            visitedNode.clear();
            visitedHyperedge.clear();
            result = {whfc::Node(0), whfc::Node(0), 0, 0};
        }
    };
}