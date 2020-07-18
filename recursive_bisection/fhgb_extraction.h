#pragma once

#include "../datastructure/flow_hypergraph_builder.h"
#include "hypergraph.h"
#include "../datastructure/queue.h"
#include "partition_ca.h"
#include "../datastructure/node_border.h"
#include "config.h"
#include "partition_threadsafe.h"

namespace whfc_rb {
    template<class Hypergraph, class PartitionImpl>
    class HypergraphBuilderExtractor {
    public:
        static constexpr NodeID invalid_node = std::numeric_limits<NodeID>::max();
        Hypergraph fhgb;

        HypergraphBuilderExtractor(const size_t maxNumNodes, const size_t maxNumEdges, const size_t maxNumPins, int seed, const PartitionConfig& config) :
                fhgb(2 * config.percentage_bfs_from_cut * maxNumNodes + 2, maxNumEdges),
                queue(maxNumNodes + 2),
                visitedNode(maxNumNodes), visitedHyperedge(maxNumEdges),
                globalToLocalID(maxNumNodes),
                mt(seed), config(config) {}

        struct ExtractorInfo {
            whfc::Node source;
            whfc::Node target;
            whfc::Flow baseCut;
            whfc::Flow cutAtStake;
        };

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

            assert(queue.empty());

            timer.start("BFS", "Extraction");

            // Add source node and run BFS in part0
            fhgb.addNode(whfc::NodeWeight(0));
            queue.push(invalid_node);
            queue.reinitialize();
            whfc::NodeWeight w0 = BreadthFirstSearch(hg, cut_hes, partition, part0, part1, maxW0, result.source, -delta, distanceFromCut, timer);

            // Add target node and run BFS in part1
            result.target = whfc::Node(fhgb.numNodes());
            fhgb.addNode(whfc::NodeWeight(0));
            queue.push(invalid_node);
            queue.reinitialize();
            whfc::NodeWeight w1 = BreadthFirstSearch(hg, cut_hes, partition, part1, part0, maxW1, result.target, delta, distanceFromCut, timer);

            timer.stop("BFS");

            timer.start("Process_Cut_Hyperedges", "Extraction");
            processCutHyperedges(hg, cut_hes, partition, part0, part1);
            timer.stop("Process_Cut_Hyperedges");

            std::vector<NodeWeight> totalWeights = partition.partitionWeights();

            fhgb.nodeWeight(result.source) = partition.partWeight(part0) - w0;
            fhgb.nodeWeight(result.target) = partition.partWeight(part1) - w1;

            timer.start("Finalize", "Extraction");
            fhgb.finalize();
            timer.stop("Finalize");

            return result;
        }

        auto localNodeIDs() const {
            return boost::irange<whfc::Node>(whfc::Node(0), whfc::Node(queue.queueEnd()));
        }

        whfc::Node global2local(const NodeID x) const {
            assert(visitedNode.contains(x));
            return globalToLocalID[x];
        }

        NodeID local2global(const whfc::Node x) const { return queue.elementAt(x); }

    private:
        LayeredQueue<NodeID> queue;
        ldc::TimestampSet<> visitedNode;
        ldc::TimestampSet<> visitedHyperedge;
        std::vector<whfc::Node> globalToLocalID;
        ExtractorInfo result;
        std::mt19937 mt;
        const PartitionConfig& config;

        void visitNode(const NodeID node, CSRHypergraph &hg, whfc::NodeWeight &w) {
            globalToLocalID[node] = whfc::Node(fhgb.numNodes());
            queue.push(node);
            visitedNode.add(node);
            fhgb.addNode(whfc::NodeWeight(hg.nodeWeight(node)));
            w += hg.nodeWeight(node);
        }

        template<typename CutEdgeRange>
        whfc::NodeWeight BreadthFirstSearch(CSRHypergraph &hg, CutEdgeRange &cut_hes, const PartitionImpl &partition,
                                            PartitionBase::PartitionID partID, PartitionBase::PartitionID otherPartID,
                                            NodeWeight maxWeight, whfc::Node terminal, whfc::HopDistance delta,
                                            whfc::DistanceFromCut& distanceFromCut, whfc::TimeReporter& timer) {
            whfc::NodeWeight w = 0;
            whfc::HopDistance d = delta;

            timer.start("Collect_Boundary_Vertices", "BFS");

            // Collect boundary vertices
            for (const HyperedgeID e : cut_hes) {
                for (NodeID v : hg.pinsOf(e)) {
                    if (!visitedNode.contains(v) && partition[v] == partID && w + hg.nodeWeight(v) <= maxWeight) {
                        visitNode(v, hg, w);
                        distanceFromCut[globalToLocalID[v]] = d;
                    }
                }
            }

            timer.stop("Collect_Boundary_Vertices");

            timer.start("Scan_Levels", "BFS");

            // Do the actual breadth first search
            while (!queue.empty()) {
                if (queue.currentLayerEmpty()) {
                    queue.finishNextLayer();
                    d += delta;
                }
                NodeID u = queue.pop();

                for (HyperedgeID e : hg.hyperedgesOf(u)) {
                    if (!visitedHyperedge.contains(e) && partition.pinsInPart(otherPartID, e) == 0 &&
                        partition.pinsInPart(partID, e) > 1) {
                        fhgb.startHyperedge(hg.hyperedgeWeight(e));
                        bool connectToTerminal = false;
                        for (NodeID v : hg.pinsOf(e)) {
                            if (partition[v] == partID) {
                                if (!visitedNode.contains(v) && w + hg.nodeWeight(v) <= maxWeight) {
                                    visitNode(v, hg, w);
                                    distanceFromCut[globalToLocalID[v]] = d;
                                }

                                if (visitedNode.contains(v)) {
                                    assert(globalToLocalID[v] < fhgb.numNodes());
                                    fhgb.addPin(globalToLocalID[v]);
                                } else {
                                    connectToTerminal = true;
                                }
                            }
                        }
                        if (connectToTerminal) {
                            assert(terminal < fhgb.numNodes());
                            fhgb.addPin(terminal);
                        }
                        visitedHyperedge.add(e);
                    }

                }
            }

            timer.stop("Scan_Levels");

            d += delta;
            distanceFromCut[terminal] = d;

            return w;
        }

        template<class CutEdgeRange>
        void processCutHyperedges(CSRHypergraph &hg, CutEdgeRange &cut_hes, const PartitionImpl &partition,
                                  const PartitionBase::PartitionID part0, const PartitionBase::PartitionID part1) {
            for (const HyperedgeID e : cut_hes) {
                assert(!visitedHyperedge.contains(e));
                assert(partition.pinsInPart(part0, e) > 0 && partition.pinsInPart(part1, e) > 0);
                bool connectToSource = false;
                bool connectToTarget = false;
                result.cutAtStake += hg.hyperedgeWeight(e);
                // visitedHyperedge.add(e);
                fhgb.startHyperedge(hg.hyperedgeWeight(e));

                for (NodeID v : hg.pinsOf(e)) {
                    if (visitedNode.contains(v)) {
                        assert(globalToLocalID[v] < fhgb.numNodes());
                        fhgb.addPin(globalToLocalID[v]);
                    } else {
                        connectToSource |= (partition[v] == part0);
                        connectToTarget |= (partition[v] == part1);
                        if (connectToSource && connectToTarget) {
                            break;
                        }
                    }
                }
                if (connectToSource && connectToTarget) {
                    fhgb.removeCurrentHyperedge();
                    result.baseCut += hg.hyperedgeWeight(e);
                } else {
                    if (connectToSource) {
                        fhgb.addPin(result.source);
                    }
                    if (connectToTarget) {
                        fhgb.addPin(result.target);
                    }
                }
            }
        }

        void initialize(uint numNodes, uint numHyperedges) {
            fhgb.clear();
            queue.clear();
            visitedNode.clear();
            visitedHyperedge.clear();
            result = {whfc::Node(0), whfc::Node(0), 0, 0};
        }
    };
}