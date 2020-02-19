#pragma once

#include "../datastructure/flow_hypergraph_builder.h"
#include "hypergraph.h"
#include "../datastructure/queue.h"

namespace whfc_rb {
    class FlowHypergraphBuilderExtractor {
    public:
        static constexpr CSRHypergraph::NodeID invalid_node = std::numeric_limits<CSRHypergraph::NodeID>::max();

        whfc::FlowHypergraphBuilder fhgb;

        FlowHypergraphBuilderExtractor(const size_t maxNumNodes, const size_t maxNumEdges, const size_t maxNumPins) :
                fhgb(maxNumNodes, maxNumEdges, maxNumPins), queue(maxNumNodes + 2), globalToLocalID(maxNumNodes) { }

        struct ExtractorInfo {
            whfc::Node source;
            whfc::Node target;
            whfc::Flow baseCut;
            whfc::Flow cutAtStake;
        };

        ExtractorInfo run(CSRHypergraph& hg, const Partition& partition, const Partition::partitionID part0, const Partition::partitionID part1, CSRHypergraph::NodeWeight maxW0, CSRHypergraph::NodeWeight maxW1) {
            initialize(hg.numNodes(), hg.numHyperedges());

            std::vector<CSRHypergraph::HyperedgeID> cut_hes = partition.getCutEdges(hg, part0, part1);

            // shuffle cut edges?

            whfc::NodeWeight w0, w1;

            // Add source node and run BFS in part0
            fhgb.addNode(whfc::NodeWeight(0));
            queue.push(invalid_node);
            queue.reinitialize();
            BreadthFirstSearch(hg, cut_hes, partition, part0, maxW0, result.source, w0);

            // Add target node and run BFS in part1
            result.target = whfc::Node::fromOtherValueType(fhgb.numNodes());
            fhgb.addNode(whfc::NodeWeight(0));
            queue.push(invalid_node);
            queue.reinitialize();
            BreadthFirstSearch(hg, cut_hes, partition, part1, maxW1, result.target, w1);

            processCutHyperedges(hg, cut_hes, partition, part0, part1);

            std::vector<CSRHypergraph::NodeWeight> totalWeights = partition.partitionWeights(hg);

            fhgb.nodeWeight(result.source) = totalWeights[0] - w0;
            fhgb.nodeWeight(result.target) = totalWeights[1] - w1;

            fhgb.finalize();

            return result;
        }

        auto localNodeIDs() const { return boost::irange<whfc::Node>(whfc::Node(0), whfc::Node::fromOtherValueType(queue.queueEnd())); }
        whfc::Node global2local(const CSRHypergraph::NodeID x) const {  assert(visitedNode[x]); return globalToLocalID[x]; }
        CSRHypergraph::NodeID local2global(const whfc::Node x) const { return queue.elementAt(x); }

    private:
        LayeredQueue<CSRHypergraph::NodeID> queue;
        std::vector<bool> visitedNode;
        std::vector<bool> visitedHyperedge;
        std::vector<whfc::Node> globalToLocalID;
        ExtractorInfo result;

        void visitNode(const CSRHypergraph::NodeID node, CSRHypergraph& hg, whfc::NodeWeight& w) {
            globalToLocalID[node] = whfc::Node::fromOtherValueType(fhgb.numNodes());
            queue.push(node);
            visitedNode[node] = true;
            fhgb.addNode(whfc::NodeWeight(hg.nodeWeight(node)));
            w += hg.nodeWeight(node);
        }

        void BreadthFirstSearch(CSRHypergraph& hg, const std::vector<CSRHypergraph::HyperedgeID>& cut_hes, const Partition& partition, uint partID, CSRHypergraph::NodeWeight maxWeight, whfc::Node terminal, whfc::NodeWeight& w) {
            w = 0;

            // Collect boundary vertices
            for (CSRHypergraph::HyperedgeID e : cut_hes) {
                for (CSRHypergraph::NodeID v : hg.pinsOf(e)) {
                    if (!visitedNode[v] && partition[v] == partID && w + hg.nodeWeight(v) <= maxWeight) {
                        visitNode(v, hg, w);
                    }
                }
            }

            // Do the actual breadth first search
            while (!queue.empty()) {
                CSRHypergraph::NodeID u = queue.pop();
                for (CSRHypergraph::HyperedgeID e : hg.hyperedgesOf(u)) {
                    if (hg.pinCount(e) == partition.pinsInPart(hg, partID, e) && hg.pinCount(e) > 1 && !visitedHyperedge[e]) {

                        visitedHyperedge[e] = true;
                        fhgb.startHyperedge(hg.hyperedgeWeight(e));
                        bool connectToTerminal = false;
                        for (CSRHypergraph::NodeID v : hg.pinsOf(e)) {
                            if (partition[v] == partID) {
                                if (!visitedNode[v] && w + hg.nodeWeight(v) <= maxWeight) {
                                    visitNode(v, hg, w);
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
                    }
                }
            }
        }

        void processCutHyperedges(CSRHypergraph& hg, const std::vector<CSRHypergraph::HyperedgeID>& cut_hes, const Partition& partition, const Partition::partitionID part0, const Partition::partitionID part1) {
            for (CSRHypergraph::HyperedgeID e : cut_hes) {
                bool connectToSource = false;
                bool connectToTarget = false;
                result.cutAtStake += hg.hyperedgeWeight(e);
                visitedHyperedge[e] = true;
                fhgb.startHyperedge(hg.hyperedgeWeight(e));

                for (CSRHypergraph::NodeID v : hg.pinsOf(e)) {
                    if (visitedNode[v]) {
                        fhgb.addPin(globalToLocalID[v]);
                    } else {
                        connectToSource |= (partition[v] == 0);
                        connectToTarget |= (partition[v] == 1);
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
            visitedNode.resize(numNodes, false);
            visitedHyperedge.clear();
            visitedHyperedge.resize(numHyperedges, false);
            result = {whfc::Node::fromOtherValueType(0), whfc::Node::fromOtherValueType(0), 0, 0};
        }
    };
}