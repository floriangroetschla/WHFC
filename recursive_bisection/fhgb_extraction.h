#pragma once

#include "../datastructure/flow_hypergraph_builder.h"
#include "hypergraph.h"
#include "../datastructure/queue.h"
#include "partition_ca.h"
#include "../datastructure/node_border.h"
#include "config.h"
#include "partition_threadsafe.h"

namespace whfc_rb {
    class FlowHypergraphBuilderExtractor {
    public:
        static constexpr NodeID invalid_node = std::numeric_limits<NodeID>::max();
        whfc::FlowHypergraphBuilder fhgb;

        FlowHypergraphBuilderExtractor(const size_t maxNumNodes, const size_t maxNumEdges, const size_t maxNumPins, int seed) :
                fhgb(maxNumNodes, maxNumEdges, maxNumPins),
                queue(maxNumNodes + 2),
                visitedNode(maxNumNodes), visitedHyperedge(maxNumEdges),
                globalToLocalID(maxNumNodes),
                mt(seed) {}

        struct ExtractorInfo {
            whfc::Node source;
            whfc::Node target;
            whfc::Flow baseCut;
            whfc::Flow cutAtStake;
        };

        template<class PartitionImpl>
        ExtractorInfo
        run(PartitionImpl &partition, const PartitionBase::PartitionID part0, const PartitionBase::PartitionID part1,
            NodeWeight maxW0, NodeWeight maxW1, const PartitionConfig& config, whfc::DistanceFromCut& distanceFromCut) {
            CSRHypergraph &hg = partition.getGraph();
            initialize(hg.numNodes(), hg.numHyperedges());
            whfc::HopDistance delta = config.distancePiercing ? 1 : 0;

            auto cut_hes = partition.getCutEdges(part0, part1);
            cut_hes.shuffle(mt);

            assert(queue.empty());

            // Add source node and run BFS in part0
            fhgb.addNode(whfc::NodeWeight(0));
            queue.push(invalid_node);
            queue.reinitialize();
            whfc::NodeWeight w0 = BreadthFirstSearch(hg, cut_hes, partition, part0, part1, maxW0, result.source, -delta, distanceFromCut);

            // Add target node and run BFS in part1
            result.target = whfc::Node::fromOtherValueType(fhgb.numNodes());
            fhgb.addNode(whfc::NodeWeight(0));
            queue.push(invalid_node);
            queue.reinitialize();
            whfc::NodeWeight w1 = BreadthFirstSearch(hg, cut_hes, partition, part1, part0, maxW1, result.target, delta, distanceFromCut);

            processCutHyperedges(hg, cut_hes, partition, part0, part1);

            std::vector<NodeWeight> totalWeights = partition.partitionWeights();

            fhgb.nodeWeight(result.source) = totalWeights[part0] - w0;
            fhgb.nodeWeight(result.target) = totalWeights[part1] - w1;

            fhgb.finalize();

            return result;
        }

        auto localNodeIDs() const {
            return boost::irange<whfc::Node>(whfc::Node(0), whfc::Node::fromOtherValueType(queue.queueEnd()));
        }

        whfc::Node global2local(const NodeID x) const {
            assert(visitedNode[x]);
            return globalToLocalID[x];
        }

        NodeID local2global(const whfc::Node x) const { return queue.elementAt(x); }

    private:
        LayeredQueue<NodeID> queue;
        boost::dynamic_bitset<> visitedNode;
        boost::dynamic_bitset<> visitedHyperedge;
        std::vector<whfc::Node> globalToLocalID;
        ExtractorInfo result;
        std::mt19937 mt;

        void visitNode(const NodeID node, CSRHypergraph &hg, whfc::NodeWeight &w) {
            globalToLocalID[node] = whfc::Node::fromOtherValueType(fhgb.numNodes());
            queue.push(node);
            visitedNode[node] = true;
            fhgb.addNode(whfc::NodeWeight(hg.nodeWeight(node)));
            w += hg.nodeWeight(node);
        }

        template<class PartitionImpl, typename CutEdgeRange>
        whfc::NodeWeight BreadthFirstSearch(CSRHypergraph &hg, CutEdgeRange &cut_hes, const PartitionImpl &partition,
                           uint partID, uint otherPartID, NodeWeight maxWeight, whfc::Node terminal, whfc::HopDistance delta, whfc::DistanceFromCut& distanceFromCut) {
            whfc::NodeWeight w = 0;
            whfc::HopDistance d = delta;

            // Collect boundary vertices
            for (const HyperedgeID e : cut_hes) {
                for (NodeID v : hg.pinsOf(e)) {
                    if (!visitedNode[v] && partition[v] == partID && w + hg.nodeWeight(v) <= maxWeight) {
                        visitNode(v, hg, w);
                        distanceFromCut[globalToLocalID[v]] = d;
                    }
                }
            }

            // Do the actual breadth first search
            while (!queue.empty()) {
                if (queue.currentLayerEmpty()) {
                    queue.finishNextLayer();
                    d += delta;
                }
                NodeID u = queue.pop();

                for (HyperedgeID e : hg.hyperedgesOf(u)) {
                    if (!visitedHyperedge[e] && partition.pinsInPart(otherPartID, e) == 0 &&
                        partition.pinsInPart(partID, e) > 1) {
                        fhgb.startHyperedge(hg.hyperedgeWeight(e));
                        bool connectToTerminal = false;
                        for (NodeID v : hg.pinsOf(e)) {
                            if (partition[v] == partID) {
                                if (!visitedNode[v] && w + hg.nodeWeight(v) <= maxWeight) {
                                    visitNode(v, hg, w);
                                    distanceFromCut[globalToLocalID[v]] = d;
                                }

                                if (visitedNode[v]) {
                                    fhgb.addPin(globalToLocalID[v]);
                                } else {
                                    connectToTerminal = true;
                                }
                            }
                        }
                        if (connectToTerminal) {
                            fhgb.addPin(terminal);
                        }
                        visitedHyperedge[e] = true;
                    }

                }
            }

            d += delta;
            distanceFromCut[terminal] = d;

            return w;
        }

        template<class PartitionImpl, class CutEdgeRange>
        void processCutHyperedges(CSRHypergraph &hg, CutEdgeRange &cut_hes, const PartitionImpl &partition,
                             const PartitionBase::PartitionID part0, const PartitionBase::PartitionID part1) {
            for (const HyperedgeID e : cut_hes) {
                assert(!visitedHyperedge[e]);
                assert(partition.pinsInPart(part0, e) > 0 && partition.pinsInPart(part1, e) > 0);
                bool connectToSource = false;
                bool connectToTarget = false;
                result.cutAtStake += hg.hyperedgeWeight(e);
                visitedHyperedge[e] = true;
                fhgb.startHyperedge(hg.hyperedgeWeight(e));

                for (NodeID v : hg.pinsOf(e)) {
                    if (visitedNode[v]) {
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
            visitedNode.reset();
            visitedHyperedge.reset();
            result = {whfc::Node::fromOtherValueType(0), whfc::Node::fromOtherValueType(0), 0, 0};
        }
    };
}