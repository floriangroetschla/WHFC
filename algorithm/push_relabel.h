
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
        using DistanceT = DistanceReachableNodes::DistanceT;


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

        std::vector<std::atomic<bool>> inQueue;

        PushRelabel(LawlerFlowHypergraph& hg, TimeReporter& timer, size_t numThreads) : hg(hg), nodes(hg.maxNumLawlerNodes()), timer(timer), inQueue(hg.maxNumLawlerNodes())
        {

        }

        void reset() {
            queue.clear();
            nodes.clear();
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

        void dischargeNode(const Node u, CutterState<Type>& cs) {
            assert(u < hg.numNodes());
            assert(hg.excess(u) > 0);
            auto& h = cs.h;

            size_t minLevel = hg.numLawlerNodes() * 2;
            for (InHeIndex inc_iter : hg.incidentHyperedgeIndices(u)) {
                InHe& inc_u = hg.getInHe(inc_iter);
                const Hyperedge e = inc_u.e;
                if (!h.areAllPinsTargetReachable__unsafe__(e)) {

                    const Node e_out = hg.edge_node_out(e);
                    Flow residual = std::min(hg.excess(u), hg.flowReceived(inc_iter));
                    if (residual > 0) {
                        if (hg.label(u) == hg.label(e_out) + 1) {
                            hg.push_node_to_edgeOut(u, inc_iter, residual);
                            if (!inQueue[e_out].exchange(true)) { assert(hg.excess(e_out) > 0); nodes.push_back(e_out); }
                        } else {
                            minLevel = std::min(minLevel, hg.label(e_out));
                        }
                    }

                    const Node e_in = hg.edge_node_in(e);
                    residual = hg.excess(u);
                    if (residual > 0) {
                        if (hg.label(u) == hg.label(e_in) + 1) {
                            hg.push_node_to_edgeIn(u, inc_iter, residual);
                            if (!inQueue[e_in].exchange(true)) { assert(hg.excess(e_in) > 0); nodes.push_back(e_in); }
                        } else {
                            minLevel = std::min(minLevel, hg.label(e_in));
                        }
                    }


                }
                assert(hg.excess(u) >= 0);
                if (hg.excess(u) == 0) break;
            }
            if (hg.excess(u) > 0) {
                // Relabel
                assert(minLevel + 1 > hg.label(u));
                hg.label(u) = minLevel + 1;
            }
        }

        void dischargeEdgeNodeIn(const Node e_in, CutterState<Type>& cs) {
            assert(hg.excess(e_in) > 0);
            auto& n = cs.n;

            size_t minLevel = hg.numLawlerNodes() * 2;
            const Hyperedge e = hg.edgeFromLawlerNode(e_in);

            for (Pin& pin : hg.pinsOf(e)) {
                Node v = pin.pin;
                if (hg.label(v) != 0 || (v == piercingNode) || n.isTarget(v)) {
                    Flow residual = std::min(hg.flowSent(pin.he_inc_iter), hg.excess(e_in));
                    if (residual > 0) {
                        if (hg.label(e_in) == hg.label(v) + 1) {
                            hg.push_edgeIn_to_node(v, pin.he_inc_iter, residual);
                            if (!n.isTarget(v) && !n.isSource(v) && !inQueue[v].exchange(true)) { assert(hg.excess(v) > 0); nodes.push_back(v); }
                        } else {
                            minLevel = std::min(minLevel, hg.label(v));
                        }
                    }
                }
                if (hg.excess(e_in) == 0) break;
            }

            Node e_out = hg.edge_node_out(e);
            Flow residual = std::min(hg.capacity(e) - hg.flow(e), hg.excess(e_in));
            if (residual > 0) {
                if (hg.label(e_in) == hg.label(e_out) + 1) {
                    hg.push_edgeIn_to_edgeOut(e_in, e_out, residual);
                    if (!inQueue[e_out].exchange(true)) { assert(hg.excess(e_out) > 0); nodes.push_back(e_out); }
                } else {
                    minLevel = std::min(minLevel, hg.label(e_out));
                }
            }

            if (hg.excess(e_in) > 0) {
                // Relabel
                assert(minLevel + 1 > hg.label(e_in));
                hg.label(e_in) = minLevel + 1;
            }
        }

        void dischargeEdgeNodeOut(const Node e_out, CutterState<Type>& cs) {
            assert(hg.excess(e_out) > 0);
            auto& n = cs.n;

            size_t minLevel = hg.numLawlerNodes() * 2;
            const Hyperedge e = hg.edgeFromLawlerNode(e_out);

            Node e_in = hg.edge_node_in(e);
            Flow residual = std::min(hg.flow(e), hg.excess(e_out));
            if (residual > 0) {
                if (hg.label(e_out) == hg.label(e_in) + 1) {
                    hg.push_edgeOut_to_edgeIn(e_out, e_in, residual);
                    if (!inQueue[e_in].exchange(true)) { assert(hg.excess(e_in) > 0); nodes.push_back(e_in); }
                } else {
                    minLevel = std::min(minLevel, hg.label(e_in));
                }
            }

            for (Pin& pin : hg.pinsOf(e)) {
                Node v = pin.pin;
                if (hg.label(v) != 0 || (v == piercingNode) || n.isTarget(v)) {
                    Flow residual = hg.excess(e_out);
                    if (residual > 0) {
                        if (hg.label(e_out) == hg.label(v) + 1) {
                            hg.push_edgeOut_to_node(v, pin.he_inc_iter, residual);
                            if (!n.isTarget(v) && !n.isSource(v) && !inQueue[v].exchange(true)) { assert(hg.excess(v) > 0); nodes.push_back(v); }
                        } else {
                            minLevel = std::min(minLevel, hg.label(v));
                        }
                    }
                }
                if (hg.excess(e_out) == 0) break;
            }

            if (hg.excess(e_out) > 0) {
                // Relabel
                assert(minLevel + 1 > hg.label(e_out));
                hg.label(e_out) = minLevel + 1;
            }
        }

        Flow exhaustFlow(CutterState<Type>& cs) {
            assert(cs.sourcePiercingNodes.size() == 1);
            auto& h = cs.h;

            hg.initialize_for_push_relabel();

            cs.flipViewDirection();
            hg.alignViewDirection();
            setLabels(cs);
            cs.flipViewDirection();
            hg.alignViewDirection();

            queue.clear();

            hg.printHypergraph(std::cout);
            hg.printExcessAndLabel();

            for (auto& sp : cs.sourcePiercingNodes) {
                size_t maxLevel = 0;
                hg.label(sp.node) = std::min<size_t>(100, hg.numLawlerNodes());  // for debugging purposes
                for (InHeIndex inc_iter : hg.incidentHyperedgeIndices(sp.node)) {
                    InHe& inc_u = hg.getInHe(inc_iter);
                    const Hyperedge e = inc_u.e;
                    if (!h.areAllPinsTargetReachable__unsafe__(e)) {

                        const Node e_out = hg.edge_node_out(e);
                        Flow residual = hg.flowReceived(inc_iter);
                        if (residual > 0) {
                            hg.push_node_to_edgeOut(sp.node, inc_iter, residual);
                            if (!inQueue[e_out].exchange(true)) { assert(hg.excess(e_out) > 0); nodes.push_back(e_out); }
                        }

                        const Node e_in = hg.edge_node_in(e);
                        residual = hg.capacity(e);
                        if (residual > 0) {
                            hg.push_node_to_edgeIn(sp.node, inc_iter, residual);
                            if (!inQueue[e_in].exchange(true)) { assert(hg.excess(e_in) > 0); nodes.push_back(e_in); }
                        }

                        if (hg.label(e_in) > maxLevel) maxLevel = hg.label(e_in);
                        if (hg.label(e_out) > maxLevel) maxLevel = hg.label(e_out);

                    }
                    //hg.label(sp.node) = maxLevel + 1;
                    piercingNode = sp.node;
                }
            }

            //hg.printHypergraph(std::cout);
            //hg.printExcessAndLabel();

            //hg.printHypergraph(std::cout);

            while (!nodes.empty()) {
                Node u = nodes.front();
                nodes.pop_front();

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

                //hg.printHypergraph(std::cout);
                //hg.printExcessAndLabel();
            }

            hg.sortPins();

            resetSourcePiercingNodeDistances(cs);

            Flow f = 0;

            for (auto& sp : cs.sourcePiercingNodes) {
                f -= hg.excess(sp.node);
            }
            cs.flowValue += f;

            //hg.printHypergraph(std::cout);
            //hg.printExcessAndLabel();

            cs.verifyFlowConstraints();

            //growReachable(cs);

            bool found_target = growReachable(cs);
            assert(!found_target);
            return f;
        }

        void setLabels(CutterState<Type>& cs) {
            queue.clear();

            hg.printHypergraph(std::cout);

            for (auto& sp : cs.sourcePiercingNodes) {
                queue.push(sp.node);
            }

            queue.finishNextLayer();

            size_t currentLabel = 1;

            auto edge_in = [&](Hyperedge e) {
                return hg.forwardView() ? hg.edge_node_in(e) : hg.edge_node_out(e);
            };

            auto edge_out = [&](Hyperedge e) {
                return hg.forwardView() ? hg.edge_node_out(e) : hg.edge_node_in(e);
            };

            while (!queue.empty()) {
                while (!queue.currentLayerEmpty()) {
                    const Node u = queue.pop();
                    if (hg.isNode(u)) {
                        for (InHe& inc_u : hg.hyperedgesOf(u)) {
                            const Hyperedge e = inc_u.e;
                            const Node e_in = edge_in(e);
                            const Node e_out = edge_out(e);
                            if (hg.flowReceived(inc_u) > 0 && hg.label(e_out) == 0) {
                                hg.label(e_out) = currentLabel;
                                queue.push(e_out);
                            }
                            if (hg.label(e_in) == 0) {
                                hg.label(e_in) = currentLabel;
                                queue.push(e_in);
                            }
                        }
                    } else if ((hg.is_edge_in(u) && hg.forwardView()) || (hg.is_edge_out(u) && !hg.forwardView())) {
                        const Hyperedge e = hg.edgeFromLawlerNode(u);
                        const Node e_out = edge_out(e);

                        if (hg.capacity(e) - hg.flow(e) > 0 && hg.label(e_out) == 0) {
                            hg.label(e_out) = currentLabel;
                            queue.push(e_out);
                        }
                        for (Pin& pin : hg.pinsSendingFlowInto(e)) {
                            if (hg.label(pin.pin) == 0 && !cs.n.isTargetReachable(pin.pin)) {
                                hg.label(pin.pin) = currentLabel;
                                queue.push(pin.pin);
                            }
                        }
                    } else {
                        const Hyperedge e = hg.edgeFromLawlerNode(u);
                        const Node e_in = edge_in(e);

                        if (hg.flow(e) > 0 && hg.label(e_in) == 0) {
                            hg.label(e_in) = currentLabel;
                            queue.push(e_in);
                        }
                        for (Pin& pin : hg.pinsOf(e)) {
                            if (hg.label(pin.pin) == 0 && !cs.n.isTargetReachable(pin.pin)) {
                                hg.label(pin.pin) = currentLabel;
                                queue.push(pin.pin);
                            }
                        }
                    }
                }
                currentLabel++;
                queue.finishNextLayer();
            }
        }

        bool growReachable(CutterState<Type>& cs, bool setLabel = false) {
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
                if (setLabel) {
                    hg.label(sp.node) = 0;
                }
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
                                if (hg.flowReceived(inc_u) > 0 && hg.label(hg.edge_node_in(e)) > hg.label(u) + 1) {
                                    hg.label(hg.edge_node_in(e)) = hg.label(u) + 1;
                                    for (const Pin& pv : hg.pinsNotSendingFlowInto(e)) {
                                        if (!n.isTarget(pv.pin)) {
                                            hg.label(pv.pin) = hg.label(hg.edge_node_in(e)) + 1;
                                        }
                                    }
                                }
                                continue;
                            }

                            if (scanAllPins) {
                                h.reachAllPins(e);
                                if (setLabel) {
                                    if (hg.flowReceived(inc_u) > 0) {
                                        hg.label(hg.edge_node_in(e)) = hg.label(u) + 1;
                                    } else {
                                        hg.label(hg.edge_node_in(e)) = hg.label(u) + 2;
                                    }
                                }
                                assert(n.distance[u] + 1 == h.outDistance[e]);
                            }

                            const bool scanFlowSending = !h.areFlowSendingPinsSourceReachable__unsafe__(e);
                            if (scanFlowSending) {
                                h.reachFlowSendingPins(e);
                                if (setLabel) hg.label(hg.edge_node_out(e)) = hg.label(u) + 1;
                                assert(n.distance[u] + 1 == h.inDistance[e]);
                            }

                            auto visit = [&](const Pin& pv, const bool sendsFlow) {
                                const Node v = pv.pin;
                                found_target |= n.isTarget(v);
                                if (!n.isTarget(v) && !n.isSourceReachable__unsafe__(v)) {
                                    n.reach(v);
                                    if (setLabel) {
                                        if (sendsFlow) {
                                            hg.label(v) = hg.label(hg.edge_node_out(e)) + 1;
                                        } else {
                                            hg.label(v) = hg.label(hg.edge_node_in(e)) + 1;
                                        }
                                    }
                                    assert(n.distance[u] + 1 == n.distance[v]);
                                    queue.push(v);
                                }
                            };

                            if (scanFlowSending)
                                for (const Pin& pv : hg.pinsSendingFlowInto(e))
                                    visit(pv, true);

                            if (scanAllPins)
                                for (const Pin& pv : hg.pinsNotSendingFlowInto(e))
                                    visit(pv, false);
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