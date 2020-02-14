#pragma once

#include "../datastructure/flow_hypergraph_builder.h"
#include "hypergraph.h"
#include <queue>


namespace whfc_rb {
    class FlowHypergraphBuilderExtractor {
    public:
        whfc::FlowHypergraphBuilder fhgb;

        FlowHypergraphBuilderExtractor(const size_t maxNumNodes, const size_t maxNumEdges, const size_t maxNumPins) :
            fhgb(maxNumNodes, maxNumEdges, maxNumPins), globalToLocalID(maxNumNodes) { }

        struct ExtractorInfo {
            whfc::Node source;
            whfc::Node target;
        };

        ExtractorInfo run(CSRHypergraph& hg, const std::vector<CSRHypergraph::HyperedgeID>& cut_hes, const std::vector<int>& partition, uint distanceFromCut) {
            initialize(hg.numNodes(), hg.numHyperedges());

            // shuffle cut edges?

            whfc::NodeWeight w0 = whfc::NodeWeight(0), w1 = whfc::NodeWeight(0);

            fhgb.addNode(whfc::NodeWeight(0));
            whfc::Node sourceNode = whfc::Node::fromOtherValueType(0);

            BreadthFirstSearch(hg, cut_hes, partition, 0, distanceFromCut, sourceNode, w0);

            whfc::Node targetNode = whfc::Node::fromOtherValueType(fhgb.numNodes());
            fhgb.addNode(whfc::NodeWeight(0));
            BreadthFirstSearch(hg, cut_hes, partition, 1, distanceFromCut, targetNode, w1);

            // Cut hyperedges
            for (CSRHypergraph::HyperedgeID e : cut_hes) {
                bool connectToSource = false;
                bool connectToTarget = false;
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
                } else {
                    if (connectToSource) {
                        fhgb.addPin(sourceNode);
                    }
                    if (connectToTarget) {
                        fhgb.addPin(targetNode);
                    }
                }
            }

            auto totalWeights = weightsOfPartitions(hg, partition);

            fhgb.nodeWeight(sourceNode) = totalWeights.first - w0;
            fhgb.nodeWeight(targetNode) = totalWeights.second - w1;

            fhgb.finalize();

            return {sourceNode, targetNode};
        }


    private:
        std::queue<CSRHypergraph::NodeID> queue;
        std::vector<bool> visitedNode;
        std::vector<bool> visitedHyperedge;
        uint numVisitedNodes = 0;
        std::vector<whfc::Node> globalToLocalID;

        void visitNode(const CSRHypergraph::NodeID node, CSRHypergraph& hg, whfc::NodeWeight& w) {
            globalToLocalID[node] = whfc::Node::fromOtherValueType(fhgb.numNodes());
            queue.push(node);
            visitedNode[node] = true;
            fhgb.addNode(whfc::NodeWeight(hg.nodeWeight(node)));
            numVisitedNodes++;
            w += hg.nodeWeight(node);
        }

        void BreadthFirstSearch(CSRHypergraph& hg, const std::vector<CSRHypergraph::HyperedgeID>& cut_hes, const std::vector<int>& partition, uint block, uint maxNodeNumber, whfc::Node terminal, whfc::NodeWeight& w) {
            numVisitedNodes = 0;

            // Collect boundary vertices
            for (CSRHypergraph::HyperedgeID e : cut_hes) {
                for (CSRHypergraph::NodeID v : hg.pinsOf(e)) {
                    if (!visitedNode[v] && partition[v] == block && numVisitedNodes < maxNodeNumber) {
                        visitNode(v, hg, w);
                    }
                }
            }

            // Do the actual breadth first search
            while (!queue.empty()) {
                CSRHypergraph::NodeID u = queue.back();
                queue.pop();
                for (CSRHypergraph::HyperedgeID e : hg.hyperedgesOf(u)) {
                    if (hg.pinCount(e) == pinsInPart(block, partition, hg, e) && pinsInPart(block, partition, hg, e) > 1 && !visitedHyperedge[e]) {

                        visitedHyperedge[e] = true;
                        fhgb.startHyperedge(hg.hyperedgeWeight(e));
                        bool connectToTerminal = false;
                        for (CSRHypergraph::NodeID v : hg.pinsOf(e)) {
                            if (partition[v] == block) {
                                if (!visitedNode[v] && numVisitedNodes < maxNodeNumber) {
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

        uint pinsInPart(uint partitionID, const std::vector<int>& partition, const CSRHypergraph& hg, CSRHypergraph::HyperedgeID e) {
            uint count = 0;
            for (CSRHypergraph::NodeID v : hg.pinsOf(e)) {
                if (partition[v] == partitionID) {
                    count++;
                }
            }
            return count;
        }

        void initialize(uint numNodes, uint numHyperedges) {
            fhgb.clear();
            queue.empty();
            visitedNode.clear();
            visitedNode.resize(numNodes, false);
            visitedHyperedge.clear();
            visitedHyperedge.resize(numHyperedges, false);
        }

        std::pair<whfc::NodeWeight, whfc::NodeWeight> weightsOfPartitions(CSRHypergraph& hg, const std::vector<int>& partition) {
            std::array<whfc::NodeWeight,2> weights;
            weights[0] = whfc::NodeWeight(0);
            weights[1] = whfc::NodeWeight(0);
            std::pair<whfc::NodeWeight, whfc::NodeWeight> result = {whfc::NodeWeight(0), whfc::NodeWeight(0)};
            for (CSRHypergraph::NodeID v : hg.nodes()) {
                weights[partition[v]] += hg.nodeWeight(v);
            }
            return {weights[0], weights[1]};
        }
    };
}