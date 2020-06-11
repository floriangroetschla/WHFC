
#include "cutter_state.h"
#include "../datastructure/queue.h"
#include "../datastructure/stack.h"
#include "../datastructure/distance_reachable_sets.h"
#include "ford_fulkerson.h"
#include "../recursive_bisection/timestamp_set.hpp"
#include "../datastructure/copyable_atomic.h"
#include "../datastructure/LawlerFlowHypergraph.h"


namespace whfc {
    class PushRelabel {
    public:
        using Type = PushRelabel;
        using ScanList = LayeredQueue<Node>;

        using ReachableNodes = DistanceReachableNodes;
        using ReachableHyperedges = DistanceReachableHyperedges;

        using Pin = FlowHypergraph::Pin;
        using InHe = FlowHypergraph::InHe;
        using PinIndexRange = FlowHypergraph::PinIndexRange;
        using DistanceT = DistanceReachableNodes::DistanceT;


        LawlerFlowHypergraph& hg;
        LayeredQueue<Node> queue;

        LayeredQueue<Node> nodes;

        int direction = 0;

        Flow upperFlowBound;

        static constexpr bool same_traversal_as_grow_assimilated = false;
        static constexpr bool grow_reachable_marks_flow_sending_pins_when_marking_all_pins = true;
        static constexpr bool log = false;

        TimeReporter& timer;

        Node piercingNode;

        std::vector<std::atomic<bool>> inQueue;

        PushRelabel(LawlerFlowHypergraph& hg, TimeReporter& timer, size_t numThreads) : timer(timer), hg(hg), inQueue(hg.numLawlerNodes())
        {
            reset();
        }

        void reset() {
            inQueue = std::vector<std::atomic<bool>>(hg.numLawlerNodes());
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
            auto& n = cs.n;
            auto& h = cs.h;

            size_t minLevel = hg.numLawlerNodes() * 2;
            for (InHeIndex inc_iter : hg.incidentHyperedgeIndices(u)) {
                InHe& inc_u = hg.getInHe(inc_iter);
                const Hyperedge e = inc_u.e;
                if (!h.areAllPinsSourceReachable(e)) {

                    const Node e_out = hg.edge_node_out(e);
                    Flow residual = std::min(hg.excess(u), hg.flowReceived(inc_u));
                    if (residual > 0) {
                        if (hg.label(u) == hg.label(e_out) + 1) {
                            hg.pushToEdgeOut(u, inc_u, residual);
                            if (!inQueue[e_out].exchange(true)) { nodes.push(e_out); }
                        } else {
                            minLevel = std::min(minLevel, hg.label(e_out));
                        }
                    }

                    const Node e_in = hg.edge_node_in(e);
                    residual = std::min(hg.excess(u), hg.capacity(e) - hg.flowSent(inc_u));
                    if (residual > 0) {
                        if (hg.label(u) == hg.label(e_in) + 1) {
                            hg.pushToEdgeIn(u, inc_u, residual);
                            if (!inQueue[e_in].exchange(true)) { nodes.push(e_in); }
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
            auto& h = cs.h;

            size_t minLevel = hg.numLawlerNodes() * 2;
            const Hyperedge e = hg.edgeFromLawlerNode(e_in);

            FlowHypergraph::PinRange  pin_range = hg.pinsOf(e);
            std::vector<Pin> pins_to_scan(pin_range.size());
            std::copy(pin_range.begin(), pin_range.end(), pins_to_scan.begin());
            for (Pin& pin : pins_to_scan) {
                InHe& inc = hg.getInHe(pin.he_inc_iter);
                Node v = pin.pin;
                if (!n.isSourceReachable(v) || (v == piercingNode)) {
                    Flow residual = std::min(hg.flowSent(inc), hg.excess(e_in));
                    if (residual > 0) {
                        if (hg.label(e_in) == hg.label(v) + 1) {
                            hg.pushFromEdgeInToNode(v, inc, residual);
                            if (!n.isTarget(v) && !n.isSource(v) && !inQueue[v].exchange(true)) { nodes.push(v); }
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
                    hg.pushFromEdgeInToEdgeOut(e_in, e_out, residual);
                    if (!inQueue[e_out].exchange(true)) { nodes.push(e_out); }
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
            auto& h = cs.h;

            size_t minLevel = hg.numLawlerNodes() * 2;
            const Hyperedge e = hg.edgeFromLawlerNode(e_out);

            Node e_in = hg.edge_node_in(e);
            Flow residual = std::min(hg.flow(e), hg.excess(e_out));
            if (residual > 0) {
                if (hg.label(e_out) == hg.label(e_in) + 1) {
                    hg.pushFromEdgeOutToEdgeIn(e_out, e_in, residual);
                    if (!inQueue[e_in].exchange(true)) { nodes.push(e_in); }
                } else {
                    minLevel = std::min(minLevel, hg.label(e_in));
                }
            }

            FlowHypergraph::PinRange pin_range = hg.pinsOf(e);
            std::vector<Pin> pins_to_scan(pin_range.size());
            std::copy(pin_range.begin(), pin_range.end(), pins_to_scan.begin());
            for (Pin& pin : pins_to_scan) {
                InHe& inc = hg.getInHe(pin.he_inc_iter);
                Node v = pin.pin;
                if (!n.isSourceReachable(v) || (v == piercingNode)) {
                    Flow residual = hg.excess(e_out);
                    if (residual > 0) {
                        if (hg.label(e_out) == hg.label(v) + 1) {
                            hg.pushFromEdgeOutToNode(v, inc, residual);
                            if (!n.isTarget(v) && !n.isSource(v) && !inQueue[v].exchange(true)) { nodes.push(v); }
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
            auto& n = cs.n;
            auto& h = cs.h;

            hg.initialize_for_push_relabel();
            queue.clear();

            for (auto& sp : cs.sourcePiercingNodes) {
                hg.label(sp.node) = hg.numLawlerNodes();
                for (InHeIndex inc_iter : hg.incidentHyperedgeIndices(sp.node)) {
                    InHe& inc_u = hg.getInHe(inc_iter);
                    const Hyperedge e = inc_u.e;
                    if (!h.areAllPinsSourceReachable(e)) {
                        Flow residual = hg.capacity(e) - hg.flowSent(inc_u);
                        if (residual > 0) {
                            hg.pushToEdgeIn(sp.node, inc_u, residual);
                            if (!inQueue[hg.edge_node_in(e)].exchange(true)) { nodes.push(hg.edge_node_in(e)); }
                        }
                    }
                    piercingNode = sp.node;
                }
            }

            while (!nodes.empty()) {
                //hg.printHypergraph(std::cout);
                //hg.printExcessAndLabel();

                Node u = nodes.pop();

                if (hg.isNode(u)) {
                    dischargeNode(u, cs);
                } else if (hg.is_edge_in(u)) {
                    dischargeEdgeNodeIn(u, cs);
                } else {
                    assert(hg.is_edge_out(u));
                    dischargeEdgeNodeOut(u, cs);
                }


                if (hg.excess(u) > 0) {
                    nodes.push(u);
                } else {
                    inQueue[u] = false;
                }
            }

            hg.printHypergraph(std::cout);
            hg.printExcessAndLabel();

            hg.sortPins();

            hg.printHypergraph(std::cout);

            resetSourcePiercingNodeDistances(cs);

            Flow f = 0;

            for (auto& sp : cs.sourcePiercingNodes) {
                f -= hg.excess(sp.node);
            }
            cs.flowValue += f;

            growReachable(cs);
            return f;
        }

        void growReachable(CutterState<Type>& cs) {
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
                            if (!scanAllPins && h.areFlowSendingPinsSourceReachable__unsafe__(e))
                                continue;

                            if (scanAllPins) {
                                h.reachAllPins(e);
                                assert(n.distance[u] + 1 == h.outDistance[e]);
                            }

                            const bool scanFlowSending = !h.areFlowSendingPinsSourceReachable__unsafe__(e);
                            if (scanFlowSending) {
                                h.reachFlowSendingPins(e);
                                assert(n.distance[u] + 1 == h.inDistance[e]);
                            }

                            auto visit = [&](const Pin& pv) {
                                const Node v = pv.pin;
                                found_target |= n.isTarget(v);
                                if (!n.isTarget(v) && !n.isSourceReachable__unsafe__(v)) {
                                    n.reach(v);
                                    assert(n.distance[u] + 1 == n.distance[v]);
                                    queue.push(v);
                                }
                            };

                            if (scanFlowSending)
                                for (const Pin& pv : hg.pinsSendingFlowInto(e))
                                    visit(pv);

                            if (scanAllPins)
                                for (const Pin& pv : hg.pinsNotSendingFlowInto(e))
                                    visit(pv);
                        }
                    }
                }

                n.hop(); h.hop(); queue.finishNextLayer();
            }

            n.lockInSourceDistance(); h.lockInSourceDistance();
            h.compareDistances(n);

            resetSourcePiercingNodeDistances(cs);
        }

        void resetSourcePiercingNodeDistances(CutterState<Type>& cs, bool reset = true) {
            for (auto& sp: cs.sourcePiercingNodes)
                cs.n.setPiercingNodeDistance(sp.node, reset);
        }

    };
}