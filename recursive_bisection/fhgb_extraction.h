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

        ExtractorInfo run(CSRHypergraph& hg, const Partition& partition, const Partition::partitionID part0, const Partition::partitionID part1, CSRHypergraph::NodeWeight maxW0, CSRHypergraph::NodeWeight maxW1) {
            initialize(hg.numNodes(), hg.numHyperedges());

            std::vector<CSRHypergraph::HyperedgeID> cut_hes = partition.getCutEdges(hg, part0, part1);

            // shuffle cut edges?

            whfc::NodeWeight w0, w1;

            fhgb.addNode(whfc::NodeWeight(0));
            whfc::Node sourceNode = whfc::Node::fromOtherValueType(0);
            BreadthFirstSearch(hg, cut_hes, partition, part0, maxW0, sourceNode, w0);

            whfc::Node targetNode = whfc::Node::fromOtherValueType(fhgb.numNodes());
            fhgb.addNode(whfc::NodeWeight(0));
            BreadthFirstSearch(hg, cut_hes, partition, part1, maxW1, targetNode, w1);

            processCutHyperedges(hg, cut_hes, partition, part0, part1, sourceNode, targetNode);

            std::vector<CSRHypergraph::NodeWeight> totalWeights = partition.partitionWeights(hg);

            fhgb.nodeWeight(sourceNode) = totalWeights[0] - w0;
            fhgb.nodeWeight(targetNode) = totalWeights[1] - w1;

            fhgb.finalize();

            return {sourceNode, targetNode};
        }


    private:
        std::queue<CSRHypergraph::NodeID> queue;
        std::vector<bool> visitedNode;
        std::vector<bool> visitedHyperedge;
        std::vector<whfc::Node> globalToLocalID;

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
                CSRHypergraph::NodeID u = queue.front();
                queue.pop();
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

        void processCutHyperedges(CSRHypergraph& hg, const std::vector<CSRHypergraph::HyperedgeID>& cut_hes, const Partition& partition, const Partition::partitionID part0, const Partition::partitionID part1, whfc::Node sourceNode, whfc::Node targetNode) {
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
        }

        void initialize(uint numNodes, uint numHyperedges) {
            fhgb.clear();
            queue.empty();
            visitedNode.clear();
            visitedNode.resize(numNodes, false);
            visitedHyperedge.clear();
            visitedHyperedge.resize(numHyperedges, false);
        }
    };
}