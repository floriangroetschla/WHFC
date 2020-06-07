
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

            size_t minLevel = hg.numLawlerNodes();
            for (InHeIndex inc_iter : hg.incidentHyperedgeIndices(u)) {
                InHe& inc_u = hg.getInHe(inc_iter);
                const Hyperedge e = inc_u.e;
                if (!h.areAllPinsSourceReachable(e)) {
                    Flow residual = std::min(hg.excess(u), hg.capacity(e) - hg.flowSent(inc_u));
                    if (residual > 0) {
                        if (hg.label(u) == hg.label(hg.edge_node(e)) + 1) {
                            hg.pushToEdge(u, inc_u, residual);
                            if (!inQueue[hg.edge_node(e)].exchange(true)) { queue.push(hg.edge_node(e)); }
                        } else {
                            minLevel = std::min(minLevel, hg.label(hg.edge_node(e)));
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

        void dischargeEdge(const Node u, CutterState<Type>& cs) {
            assert(hg.excess(u) > 0);
            auto& n = cs.n;
            auto& h = cs.h;

            size_t minLevel = hg.numLawlerNodes();
            const Hyperedge e = hg.edgeFromLawlerNode(u);
            for (Pin& pin : hg.pinsOf(e)) {
                InHe& inc = hg.getInHe(pin.he_inc_iter);
                Node v = pin.pin;
                std::cout << "n.isSourceReachable(" << v << ") = " << n.isSourceReachable(v) << ", n.isSource(" << v << ") = " << n.isSource(v) << std::endl;
                if (!n.isSourceReachable(v) || (v == piercingNode)) {
                    Flow residual = std::min({hg.excess(u), hg.capacity(e) - hg.flow(e) + hg.absoluteFlowSent(inc), hg.capacity(e) + hg.flowSent(inc)});
                    if (residual > 0) {
                        if (hg.label(u) == hg.label(v) + 1) {
                            hg.pushToNode(v, inc, residual);
                            if (!n.isTarget(v) && !n.isSource(v) && !inQueue[v].exchange(true)) { queue.push(v); }
                        } else {
                            minLevel = std::min(minLevel, hg.label(v));
                        }
                    }
                }
                if (hg.excess(u) == 0) break;
            }
            if (hg.excess(u) > 0) {
                // Relabel
                assert(minLevel + 1 > hg.label(u));
                hg.label(u) = minLevel + 1;
            }
        }

        Flow exhaustFlow(CutterState<Type>& cs) {
            assert(cs.sourcePiercingNodes.size() == 1);
            std::cout << "Exhaust flow" << std::endl;
            hg.printHypergraph(std::cout);
            auto& n = cs.n;
            auto& h = cs.h;

            hg.initialize_for_push_relabel();
            queue.clear();

            for (auto& sp : cs.sourcePiercingNodes) {
                std::cout << "piercing node: " << sp.node << std::endl;
                std::cout << "n.isSourceReachable(" << sp.node << ") = " << n.isSourceReachable(sp.node) << ", n.isSource(" << sp.node << ") = " << n.isSource(sp.node) << std::endl;
                hg.label(sp.node) = hg.numLawlerNodes();
                for (InHeIndex inc_iter : hg.incidentHyperedgeIndices(sp.node)) {
                    InHe& inc_u = hg.getInHe(inc_iter);
                    const Hyperedge e = inc_u.e;
                    if (!h.areAllPinsSourceReachable(e)) {
                        Flow residual = hg.capacity(e) - hg.flowSent(inc_u);
                        if (residual > 0) {
                            hg.pushToEdge(sp.node, inc_u, residual);
                            if (!inQueue[hg.edge_node(e)].exchange(true)) { queue.push(hg.edge_node(e)); }
                        }
                    }
                    piercingNode = sp.node;
                }
            }

            while (!queue.empty()) {
                hg.printHypergraph(std::cout);
                hg.printExcessAndLabel();
                const Node u = queue.pop();
                if (hg.isNode(u)) {
                    // This is a node in the original graph
                    dischargeNode(u, cs);
                    if (hg.excess(u) > 0) {
                        queue.push(u);
                    } else {
                        inQueue[u] = false;
                    }
                } else {
                    // This node represents an edge
                    dischargeEdge(u, cs);
                    if (hg.excess(u) > 0) {
                        queue.push(u);
                    } else {
                        inQueue[u] = false;
                    }
                }
            }

            hg.printHypergraph(std::cout);
            hg.printExcessAndLabel();
            resetSourcePiercingNodeDistances(cs);

            Flow f = 0;

            for (auto& sp : cs.sourcePiercingNodes) {
                f -= hg.excess(sp.node);
            }
            cs.flowValue += f;
            return f;
        }

        Flow growReachable(CutterState<Type>& cs) {
            cs.clearForSearch();
            ReachableNodes& n = cs.n;
            ReachableHyperedges& h = cs.h;
            queue.clear();
            for (auto& s : cs.sourcePiercingNodes) {
                n.setPiercingNodeDistance(s.node, false);
                queue.push(s.node);
            }

            n.hop(); h.hop();

            while (!queue.empty()) {
                const Node u = queue.pop();
                for (InHeIndex inc_u_iter : hg.incidentHyperedgeIndices(u)) {
                    const InHe& inc_u = hg.getInHe(inc_u_iter);
                    const Hyperedge e = inc_u.e;
                    if (!h.areAllPinsSourceReachable(e)) {
                        const bool scanAllPins = !hg.isSaturated(e) || hg.flowReceived(inc_u) > 0;
                        if (scanAllPins)
                            h.reachAllPins(e);
                        else if (h.areFlowSendingPinsSourceReachable(e))
                            continue;
                        else
                            h.reachFlowSendingPins(e);

                        for (const Pin& pv : scanAllPins ? hg.pinsOf(e) : hg.pinsSendingFlowInto(e)) {
                            const Node v = pv.pin;
                            if (!n.isSourceReachable(v) && !n.isTarget(v)) {
                                n.reach(v);		//this has to be after the return if v is target, to avoid overwriting the targetSettled timestamp with sourceReachable
                                queue.push(v);
                            }
                        }
                    }
                }
            }

            n.hop(); h.hop();
            n.lockInSourceDistance(); h.lockInSourceDistance();
            h.compareDistances(n);
            resetSourcePiercingNodeDistances(cs);
            return 0;
        }

        void resetSourcePiercingNodeDistances(CutterState<Type>& cs, bool reset = true) {
            for (auto& sp: cs.sourcePiercingNodes)
                cs.n.setPiercingNodeDistance(sp.node, reset);
        }

    };
}