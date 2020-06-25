
#include "cutter_state.h"
#include "../datastructure/queue.h"
#include "../datastructure/stack.h"
#include "../datastructure/distance_reachable_sets.h"
#include "ford_fulkerson.h"
#include "../recursive_bisection/timestamp_set.hpp"
#include "../datastructure/copyable_atomic.h"
#include "../datastructure/lawler_flow_hypergraph.h"
#include <boost/circular_buffer.hpp>


namespace whfc {
    class PushRelabel {
    public:
        using Type = PushRelabel;
        using ScanList = LayeredQueue<Node>;

        using ReachableNodes = DistanceReachableNodes;
        using ReachableHyperedges = DistanceReachableHyperedges;
        using Hypergraph = LawlerFlowHypergraph;

        using Pin = FlowHypergraph::Pin;
        using InHe = FlowHypergraph::InHe;
        using PinIndexRange = FlowHypergraph::PinIndexRange;
        using PinIterator = FlowHypergraph::PinIterator;
        using DistanceT = DistanceReachableNodes::DistanceT;

        static constexpr size_t ALPHA = 6;
        static constexpr size_t BETA = 12;
        float globalUpdateFreq = 0.5;

        LawlerFlowHypergraph& hg;
        LayeredQueue<Node> queue;
        boost::circular_buffer<Node> nodes;

        int direction = 0;

        Flow upperFlowBound;

        static constexpr bool same_traversal_as_grow_assimilated = false;
        static constexpr bool grow_reachable_marks_flow_sending_pins_when_marking_all_pins = true;
        static constexpr bool log = false;

        TimeReporter& timer;

        Node piercingNode;
        Node target;

        std::vector<bool> inQueue;

        // for nodes
        std::vector<InHeIndex> current_hyperedge;

        // for hyperedges
        std::vector<PinIterator> current_pin_e_in, current_pin_e_out;

        size_t numPushes, numRelabel, workSinceLastRelabel;

        PushRelabel(LawlerFlowHypergraph& hg, TimeReporter& timer, size_t numThreads) : hg(hg), nodes(hg.maxNumLawlerNodes()), timer(timer), inQueue(hg.maxNumLawlerNodes()),
            current_hyperedge(hg.maxNumNodes, InHeIndex::Invalid()), current_pin_e_in(hg.maxNumHyperedges), current_pin_e_out(hg.maxNumHyperedges)
        {

        }

        void reset() {
            queue.clear();
            nodes.clear();
            numPushes = 0;
            numRelabel = 0;
            workSinceLastRelabel = 0;
            globalUpdateFreq = 0.5;
        }

        ScanList& getScanList() {
            return queue;
        }

        Flow recycleDatastructuresFromGrowReachablePhase(CutterState<Type> &cs) {
            throw std::logic_error("Not implemented");
        }

        Flow growFlowOrSourceReachable(CutterState<Type>& cs) {
            throw std::logic_error("Not implemented");
        }

        template<bool pushFlow>
        __attribute__((always_inline)) inline void iterateOverNode(const Node u, CutterState<Type>& cs) {
            if (!pushFlow) workSinceLastRelabel += BETA;
            size_t minLevel = hg.numLawlerNodes() * 2;
            assert(pushFlow || hg.excess(u) > 0);
            InHeIndex he = pushFlow ? current_hyperedge[u] : hg.beginIndexHyperedges(u);
            if (he == InHeIndex::Invalid()) he = hg.beginIndexHyperedges(u);
            for (; he < hg.endIndexHyperedges(u); ++he) {
                InHe& inc_u = hg.getInHe(he);
                const Hyperedge e = inc_u.e;

                if (!cs.h.areAllPinsSourceReachable__unsafe__(e)) {
                    if (!pushFlow) workSinceLastRelabel += 2;
                    const Node e_out = hg.edge_node_out(e);
                    Flow residual = std::min(hg.excess(u), hg.flowReceived(he));
                    if (residual > 0) {
                        if (hg.label(u) == hg.label(e_out) + 1) {
                            assert(pushFlow);
                            numPushes++;
                            hg.push_node_to_edgeOut(u, he, residual);
                            if (!inQueue[e_out]) { assert(hg.excess(e_out) > 0); nodes.push_back(e_out); inQueue[e_out] = true; }
                        } else {
                            assert(hg.label(u) <= hg.label(e_out));
                            if (hg.label(e_out) != 0) minLevel = std::min(minLevel, hg.label(e_out));
                        }
                    }

                    const Node e_in = hg.edge_node_in(e);
                    residual = hg.excess(u);
                    if (residual > 0) {
                        if (hg.label(u) == hg.label(e_in) + 1) {
                            assert(pushFlow);
                            numPushes++;
                            hg.push_node_to_edgeIn(u, he, residual);
                            if (!inQueue[e_in]) { assert(hg.excess(e_in) > 0); nodes.push_back(e_in); inQueue[e_in] = true; }
                        } else {
                            assert(hg.label(u) <= hg.label(e_in));
                            if (hg.label(e_in) != 0) minLevel = std::min(minLevel, hg.label(e_in));
                        }
                    }
                    assert(hg.excess(u) >= 0);
                    if (hg.excess(u) == 0) break;
                }
            }
            current_hyperedge[u] = pushFlow ? he : hg.beginIndexHyperedges(u);
            if (!pushFlow) {
                assert(minLevel + 1 > hg.label(u));
                hg.label(u) = minLevel + 1;
            }
        }

        template<bool pushFlow>
        __attribute__((always_inline)) inline void iterateOverEdgeNodeIn(const Node e_in, CutterState<Type>& cs) {
            assert(hg.is_edge_in(e_in));
            assert(pushFlow || hg.excess(e_in) > 0);
            if (!pushFlow) workSinceLastRelabel += BETA;
            size_t minLevel = hg.numLawlerNodes() * 2;
            const Hyperedge e = hg.edgeFromLawlerNode(e_in);
            PinIterator pin_it = pushFlow ? std::max(current_pin_e_in[e], hg.beginPins(e)) : hg.beginPins(e);

            for (; pin_it < hg.endPins(e); ++pin_it) {
                const Pin pin = *pin_it;
                const Node v = pin.pin;
                if (!cs.n.isSourceReachable(v) || (v == piercingNode) || v == target) {
                    if (!pushFlow) workSinceLastRelabel++;
                    Flow residual = std::min(hg.flowSent(pin.he_inc_iter), hg.excess(e_in));
                    if (residual > 0) {
                        if (hg.label(e_in) == hg.label(v) + 1) {
                            assert(pushFlow);
                            numPushes++;
                            hg.push_edgeIn_to_node(v, pin.he_inc_iter, residual);
                            if (!(v == target) && !(v == piercingNode) && !inQueue[v]) { assert(hg.excess(v) > 0); nodes.push_back(v); inQueue[v] = true; }
                        } else {
                            assert(hg.label(e_in) <= hg.label(v) || v == piercingNode);
                            minLevel = std::min(minLevel, hg.label(v));
                        }
                    }
                }
                if (hg.excess(e_in) == 0) break;
            }

            if (hg.excess(e_in) > 0 && pin_it == hg.endPins(e)) {
                if (!pushFlow) workSinceLastRelabel++;
                const Node e_out = hg.edge_node_out(e);
                Flow residual = std::min(hg.capacity(e) - hg.flow(e), hg.excess(e_in));
                if (residual > 0) {
                    if (hg.label(e_in) == hg.label(e_out) + 1) {
                        assert(pushFlow);
                        numPushes++;
                        hg.push_edgeIn_to_edgeOut(e_in, e_out, residual);
                        if (!inQueue[e_out]) { assert(hg.excess(e_out) > 0); nodes.push_back(e_out); inQueue[e_out] = true; }
                    } else {
                        assert(hg.label(e_in) <= hg.label(e_out));
                        minLevel = std::min(minLevel, hg.label(e_out));
                    }
                }
            }
            current_pin_e_in[e] = pushFlow ? pin_it : hg.beginPins(e);

            if (!pushFlow) {
                assert(minLevel + 1 > hg.label(e_in));
                hg.label(e_in) = minLevel + 1;
            }
        }

        template<bool pushFlow>
        __attribute__((always_inline)) inline void iterateOverEdgeNodeOut(const Node e_out, CutterState<Type>& cs) {
            assert(hg.excess(e_out) > 0);
            if (!pushFlow) workSinceLastRelabel += BETA;
            size_t minLevel = hg.numLawlerNodes() * 2;
            const Hyperedge e = hg.edgeFromLawlerNode(e_out);
            PinIterator pin_it = pushFlow ? std::max(current_pin_e_out[e], hg.beginPins(e)) : hg.beginPins(e);

            for (; pin_it < hg.endPins(e); ++pin_it) {
                const Pin pin = *pin_it;
                const Node v = pin.pin;
                if (!cs.n.isSourceReachable(v) || (v == piercingNode) || (v == target)) {
                    if (!pushFlow) workSinceLastRelabel++;
                    Flow residual = hg.excess(e_out);
                    if (residual > 0) {
                        if (hg.label(e_out) == hg.label(v) + 1) {
                            assert(pushFlow);
                            numPushes++;
                            hg.push_edgeOut_to_node(v, pin.he_inc_iter, residual);
                            if (!(v == target) && !(v == piercingNode) && !inQueue[v]) { assert(hg.excess(v) > 0); nodes.push_back(v); inQueue[v] = true; }
                        } else {
                            assert(hg.label(e_out) <= hg.label(v) || piercingNode == v);
                            minLevel = std::min(minLevel, hg.label(v));
                        }
                    }
                }
                if (hg.excess(e_out) == 0) break;
            }

            if (hg.excess(e_out) > 0 && pin_it == hg.endPins(e)) {
                if (!pushFlow) workSinceLastRelabel++;
                const Node e_in = hg.edge_node_in(e);
                Flow residual = std::min(hg.flow(e), hg.excess(e_out));
                if (residual > 0) {
                    if (hg.label(e_out) == hg.label(e_in) + 1) {
                        assert(pushFlow);
                        numPushes++;
                        hg.push_edgeOut_to_edgeIn(e_out, e_in, residual);
                        if (!inQueue[e_in]) {
                            assert(hg.excess(e_in) > 0);
                            nodes.push_back(e_in);
                            inQueue[e_in] = true;
                        }
                    } else {
                        assert(hg.label(e_out) <= hg.label(e_in));
                        minLevel = std::min(minLevel, hg.label(e_in));
                    }
                }
            }

            current_pin_e_out[e] = pushFlow ? pin_it : hg.beginPins(e);

            if (!pushFlow) {
                assert(minLevel + 1 > hg.label(e_out));
                hg.label(e_out) = minLevel + 1;
            }


        }

        void dischargeNode(const Node u, CutterState<Type>& cs) {
            assert(u < hg.numNodes());
            assert(hg.excess(u) > 0);

            iterateOverNode<true>(u, cs);
            if (hg.excess(u) > 0) {
                // Relabel
                numRelabel++;
                iterateOverNode<false>(u, cs);
            }
        }

        void dischargeEdgeNodeIn(const Node e_in, CutterState<Type>& cs) {
            assert(hg.excess(e_in) > 0);
            assert(hg.is_edge_in(e_in));

            iterateOverEdgeNodeIn<true>(e_in, cs);
            if (hg.excess(e_in) > 0) {
                // Relabel
                numRelabel++;
                iterateOverEdgeNodeIn<false>(e_in, cs);
            }
        }

        void dischargeEdgeNodeOut(const Node e_out, CutterState<Type>& cs) {
            assert(hg.excess(e_out) > 0);

            iterateOverEdgeNodeOut<true>(e_out, cs);
            if (hg.excess(e_out) > 0) {
                // Relabel
                numRelabel++;
                iterateOverEdgeNodeOut<false>(e_out, cs);
            }
        }

        void pushRelabel(CutterState<Type>& cs, const size_t nm, const size_t n, const bool keepNodesWithHighLabel) {
            while (!nodes.empty()) {
                Node u = nodes.front();
                nodes.pop_front();

                if (!keepNodesWithHighLabel && hg.label(u) >= n) {
                    continue;
                }

                if (hg.isNode(u)) {
                    dischargeNode(u, cs);
                } else if (hg.is_edge_in(u)) {
                    dischargeEdgeNodeIn(u, cs);
                } else {
                    assert(hg.is_edge_out(u));
                    dischargeEdgeNodeOut(u, cs);
                }


                if (hg.excess(u) > 0) {
                    nodes.push_back(u);
                } else {
                    inQueue[u] = false;
                }

                if (!keepNodesWithHighLabel && workSinceLastRelabel * globalUpdateFreq > nm) {
                    globalUpdate(cs, n);
                    workSinceLastRelabel = 0;
                }
            }
        }

        Flow exhaustFlow(CutterState<Type>& cs) {
            assert(cs.sourcePiercingNodes.size() == 1);
            timer.start("exhaustFlow");

            hg.alignViewDirection();

            cs.clearForSearch();

            timer.start("initialize", "exhaustFlow");
            hg.initialize_for_push_relabel();
            timer.stop("initialize");

            timer.start("setLabels", "exhaustFlow");
            size_t n = hg.numNodes() - cs.n.numSettledNodes + 2 * hg.numHyperedges();
            std::array<size_t, 2> numScans = globalUpdate(cs, n);
            timer.stop("setLabels");

            size_t nm = ALPHA * numScans[0] + numScans[1];

            target = cs.targetPiercingNodes.begin()->node;

            queue.clear();

            timer.start("mainLoop", "exhaustFlow");
            for (auto& sp : cs.sourcePiercingNodes) {
                for (InHeIndex inc_iter : hg.incidentHyperedgeIndices(sp.node)) {
                    InHe& inc_u = hg.getInHe(inc_iter);
                    const Hyperedge e = inc_u.e;
                    if (!cs.h.areAllPinsSourceReachable__unsafe__(e)) {
                        const Node e_out = hg.edge_node_out(e);
                        if (hg.label(e_out) != n) {
                            Flow residual = hg.flowReceived(inc_iter);
                            if (residual > 0) {
                                hg.push_node_to_edgeOut(sp.node, inc_iter, residual);
                                if (!inQueue[e_out]) { assert(hg.excess(e_out) > 0); nodes.push_back(e_out); inQueue[e_out] = true; }
                            }
                        }

                        const Node e_in = hg.edge_node_in(e);
                        if (hg.label(e_in) != n) {
                            Flow residual = hg.capacity(e);
                            if (residual > 0) {
                                hg.push_node_to_edgeIn(sp.node, inc_iter, residual);
                                if (!inQueue[e_in]) { assert(hg.excess(e_in) > 0); nodes.push_back(e_in); inQueue[e_in] = true; }
                            }
                        }
                    }
                }
                piercingNode = sp.node;
            }

            timer.start("phase1", "mainLoop");
            // Phase 1
            pushRelabel(cs, nm, numScans[0], false);
            timer.stop("phase1");

            assert(nodes.empty());

            timer.start("phase2", "mainLoop");
            // Phase 2
            hg.equalizeLabels(numScans[0]);
            for (Node node_id(0); node_id < hg.numLawlerNodes(); ++node_id) {
                if (hg.excess(node_id) > 0 && node_id != piercingNode && node_id != target) {
                    assert(inQueue[node_id]);
                    nodes.push_back(node_id);
                } else {
                    assert(!inQueue[node_id]);
                }
            }
            pushRelabel(cs, nm, numScans[0], true);
            timer.stop("phase2");
            timer.stop("mainLoop");

            timer.start("sortPins", "exhaustFlow");
            //hg.sortPins();
            hg.writeBackFlow();
            timer.stop("sortPins");

            resetSourcePiercingNodeDistances(cs);

            Flow f = 0;

            for (auto& sp : cs.sourcePiercingNodes) {
                f -= hg.excess(sp.node);
            }
            cs.flowValue += f;

            cs.verifyFlowConstraints();

            timer.start("growReachable", "exhaustFlow");
            bool found_target = growReachable(cs);
            assert(!found_target);
            timer.stop("growReachable");
            timer.stop("exhaustFlow");

            std::cout << "numPushes: " << numPushes << ", numRelabel: " << numRelabel << ", numLawlerNodes: " << hg.numLawlerNodes() << std::endl;

            return f;
        }

        size_t scannedNodesLastIteration = 0;

        std::array<size_t, 2> globalUpdate(CutterState<Type>& cs, size_t n) {
            size_t scannedNodes = 0;
            size_t scannedEdges = 0;
            queue.clear();

            hg.equalizeLabels(n);

            std::fill(current_pin_e_in.begin(), current_pin_e_in.end(), PinIterator(0));
            std::fill(current_pin_e_out.begin(), current_pin_e_out.end(), PinIterator(0));
            std::fill(current_hyperedge.begin(), current_hyperedge.end(), InHeIndex::Invalid());

            // Source and target of the bfs
            const Node source = *cs.targetPiercingNodes.begin()->node;
            const Node target = *cs.sourcePiercingNodes.begin()->node;

            queue.push(source);
            hg.label(source) = 0;

            queue.finishNextLayer();

            size_t currentLabel = 1;

            while (!queue.empty()) {
                while (!queue.currentLayerEmpty()) {
                    scannedNodes++;
                    const Node u = queue.pop();
                    if (hg.isNode(u)) {
                        for (InHe& inc_u : hg.hyperedgesOf(u)) {
                            scannedEdges += 2;
                            const Hyperedge e = inc_u.e;
                            const Node e_in = hg.edge_node_in(e);
                            const Node e_out = hg.edge_node_out(e);
                            if (!cs.h.areAllPinsSourceReachable__unsafe__(e)) {
                                if (hg.flowSent(inc_u) > 0 && hg.label(e_in) == n) {
                                    hg.label(e_in) = currentLabel;
                                    current_pin_e_in[e] = hg.beginPins(e);
                                    queue.push(e_in);
                                }
                                if (hg.label(e_out) == n) {
                                    hg.label(e_out) = currentLabel;
                                    current_pin_e_out[e] = hg.beginPins(e);
                                    queue.push(e_out);
                                }
                            }
                        }
                    } else if (hg.is_edge_out(u)) {
                        const Hyperedge e = hg.edgeFromLawlerNode(u);
                        const Node e_in = hg.edge_node_in(e);

                        scannedEdges++;
                        if (hg.capacity(e) - hg.flow(e) > 0 && hg.label(e_in) == n) {
                            hg.label(e_in) = currentLabel;
                            current_pin_e_in[e] = hg.beginPins(e);
                            queue.push(e_in);
                        }
                        for (Pin& pin : hg.pinsOf(e)) {
                            scannedEdges++;
                            if (hg.label(pin.pin) == n && (hg.flowReceived(pin.he_inc_iter) > 0) && (!cs.n.isSourceReachable__unsafe__(pin.pin))) {
                                hg.label(pin.pin) = currentLabel;
                                current_hyperedge[pin.pin] = hg.beginIndexHyperedges(pin.pin);
                                queue.push(pin.pin);
                            }
                        }
                    } else {
                        const Hyperedge e = hg.edgeFromLawlerNode(u);
                        const Node e_out = hg.edge_node_out(e);

                        scannedEdges++;
                        if (hg.flow(e) > 0 && hg.label(e_out) == n) {
                            hg.label(e_out) = currentLabel;
                            current_pin_e_out[e] = hg.beginPins(e);
                            queue.push(e_out);
                        }
                        for (Pin& pin : hg.pinsOf(e)) {
                            scannedEdges++;
                            if (hg.label(pin.pin) == n && (!cs.n.isSourceReachable__unsafe__(pin.pin))) {
                                hg.label(pin.pin) = currentLabel;
                                current_hyperedge[pin.pin] = hg.beginIndexHyperedges(pin.pin);
                                queue.push(pin.pin);
                            }
                        }
                    }
                }
                currentLabel++;
                queue.finishNextLayer();
            }
            if (scannedNodes == scannedNodesLastIteration) globalUpdateFreq /= 2;
            scannedNodesLastIteration = scannedNodes;
            return {scannedNodes, scannedEdges};
        }

        bool growReachable(CutterState<Type>& cs) {
            hg.alignViewDirection();
            cs.clearForSearch();
            auto& n = cs.n;
            auto& h = cs.h;
            queue.clear();
            bool found_target = false;

            for (auto& sp : cs.sourcePiercingNodes) {
                n.setPiercingNodeDistance(sp.node, false);
                assert(n.isSourceReachable(sp.node));
                queue.push(sp.node);
            }
            n.hop(); h.hop(); queue.finishNextLayer();

            while (!queue.empty()) {
                while (!queue.currentLayerEmpty()) {
                    const Node u = queue.pop();
                    for (InHe& inc_u : hg.hyperedgesOf(u)) {
                        const Hyperedge e = inc_u.e;
                        if (!h.areAllPinsSourceReachable__unsafe__(e)) {
                            const bool scanAllPins = !hg.isSaturated(e) || hg.flowReceived(inc_u) > 0;
                            if (!scanAllPins && h.areFlowSendingPinsSourceReachable__unsafe__(e)) {
                                continue;
                            }

                            if (scanAllPins) {
                                h.reachAllPins(e);
                                assert(n.distance[u] + 1 == h.outDistance[e]);
                            }

                            const bool scanFlowSending = !h.areFlowSendingPinsSourceReachable__unsafe__(e);
                            if (scanFlowSending) {
                                h.reachFlowSendingPins(e);
                                assert(n.distance[u] + 1 == h.inDistance[e]);
                            }

                            auto visit = [&](const Pin& pv, const bool sendsFlow) {
                                const Node v = pv.pin;
                                found_target |= n.isTarget(v);
                                if (!n.isTarget(v) && !n.isSourceReachable__unsafe__(v)) {
                                    n.reach(v);
                                    assert(n.distance[u] + 1 == n.distance[v]);
                                    queue.push(v);
                                }
                            };

                            if (scanFlowSending || scanAllPins) {
                                for (const Pin& pv : hg.pinsOf(e)) {
                                    if (hg.flowSent(pv.he_inc_iter) > 0 && scanFlowSending) {
                                        visit(pv, true);
                                    } else if (hg.flowSent(pv.he_inc_iter) == 0 && scanAllPins) {
                                        visit(pv, false);
                                    }
                                }
                            }
                        }
                    }
                }

                n.hop(); h.hop(); queue.finishNextLayer();
            }

            n.lockInSourceDistance(); h.lockInSourceDistance();
            h.compareDistances(n);

            resetSourcePiercingNodeDistances(cs);
            return found_target;
        }

        void resetSourcePiercingNodeDistances(CutterState<Type>& cs, bool reset = true) {
            for (auto& sp: cs.sourcePiercingNodes)
                cs.n.setPiercingNodeDistance(sp.node, reset);
        }

    };
}