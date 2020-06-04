
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

        std::vector<std::atomic<bool>> inQueue;

        PushRelabel(LawlerFlowHypergraph& hg, TimeReporter& timer, size_t numThreads) : timer(timer), hg(hg), inQueue(hg.numLawlerNodes())
        {
            reset();
        }

        void reset() {
            inQueue = std::vector<std::atomic<bool>>(hg.numLawlerNodes(), false);
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
                        if (hg.label(u) == )
                        hg.pushToEdge(u, inc_u, residual);
                        if (inQueue[hg.edge_in_node(e)].exchange(true)) { queue.push(hg.edge_in_node(e)); }
                        if (inQueue[hg.edge_out_node(e)].exchange(true)) { queue.push(hg.edge_out_node(e)); }
                    }
                }
                assert(hg.excess(u) >= 0);
                if (hg.excess(u) == 0) break;
            }

        }

        Flow exhaustFlow(CutterState<Type>& cs) {
            auto& n = cs.n;
            auto& h = cs.h;
            Flow f = 0;

            hg.initialize_for_push_relabel();
            queue.clear();

            for (auto& sp : cs.sourcePiercingNodes) {
                hg.label(sp.node) = hg.numLawlerNodes();
                for (InHeIndex inc_iter : hg.incidentHyperedgeIndices(sp.node)) {
                    InHe& inc_u = hg.getInHe(inc_iter);
                    const Hyperedge e = inc_u.e;
                    if (!h.areAllPinsSourceReachable(e)) {
                        hg.pushToEdge(sp.node, inc_u, hg.capacity(inc_u.e));
                        if (inQueue[hg.edge_in_node(e)].exchange(true)) { queue.push(hg.edge_in_node(e)); }
                    }
                }
            }

            while (!queue.empty()) {
                const Node u = queue.pop();
                dischargeNode(u, cs);
                inQueue[u] = false;
            }

            return f;
        }

        Flow growReachable(CutterState<Type>& cs) {
            cs.clearForSearch();
            ReachableNodes& n = cs.n;
            ReachableHyperedges& h = cs.h;
            queue.clear();
            for (auto& s : cs.sourcePiercingNodes)
                queue.push(s.node);

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
                            if (!n.isSourceReachable(v)) {		//don't do VD label propagation
                                n.reach(v);		//this has to be after the return if v is target, to avoid overwriting the targetSettled timestamp with sourceReachable
                                queue.push(v);
                            }
                        }
                    }
                }
            }
            return 0;
        }

    private:
        void resetSourcePiercingNodeDistances(CutterState<Type>& cs, bool reset = true) {
            for (auto& sp: cs.sourcePiercingNodes)
                cs.n.setPiercingNodeDistance(sp.node, reset);
        }

    };
}