
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

        std::vector<std::atomic<bool>> inQueue;

        // for nodes
        std::vector<InHeIndex> current_hyperedge;

        // for hyperedges
        std::vector<PinIterator> current_pin_e_in, current_pin_e_out;

        size_t numPushes, numRelabel, numEdgeScans;

        std::unique_ptr<tbb::enumerable_thread_specific<std::vector<Node>>> thisLayer_thread_specific;
        std::unique_ptr<tbb::enumerable_thread_specific<std::vector<Node>>> nextLayer_thread_specific;

        PushRelabel(LawlerFlowHypergraph& hg, TimeReporter& timer, size_t numThreads) : hg(hg), nodes(hg.maxNumLawlerNodes()), timer(timer), inQueue(hg.maxNumLawlerNodes()),
            current_hyperedge(hg.maxNumNodes, InHeIndex::Invalid()), current_pin_e_in(hg.maxNumHyperedges), current_pin_e_out(hg.maxNumHyperedges),
            thisLayer_thread_specific(new tbb::enumerable_thread_specific<std::vector<Node>>()),
            nextLayer_thread_specific(new tbb::enumerable_thread_specific<std::vector<Node>>())
        {

        }

        void reset() {
            queue.clear();
            nodes.clear();
            numPushes = 0;
            numRelabel = 0;
            numEdgeScans = 0;
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

        void iterateOverNode(const Node u, CutterState<Type>& cs, const bool pushFlow) {
            size_t minLevel = hg.numLawlerNodes() * 2;
            assert(pushFlow || hg.excess(u) > 0);
            InHeIndex he = pushFlow ? current_hyperedge[u] : hg.beginIndexHyperedges(u);
            if (he == InHeIndex::Invalid()) he = hg.beginIndexHyperedges(u);
            for (; he < hg.endIndexHyperedges(u); ++he) {
                InHe& inc_u = hg.getInHe(he);
                const Hyperedge e = inc_u.e;

                if (!cs.h.areAllPinsSourceReachable__unsafe__(e)) {
                    numEdgeScans += 2;
                    const Node e_out = hg.edge_node_out(e);
                    Flow residual = std::min(hg.excess(u), hg.flowReceived(he));
                    if (residual > 0) {
                        if (hg.label(u) == hg.label(e_out) + 1) {
                            assert(pushFlow);
                            numPushes++;
                            hg.push_node_to_edgeOut(u, he, residual);
                            if (!inQueue[e_out].exchange(true)) { assert(hg.excess(e_out) > 0); nodes.push_back(e_out); }
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
                            if (!inQueue[e_in].exchange(true)) { assert(hg.excess(e_in) > 0); nodes.push_back(e_in); }
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

        void iterateOverEdgeNodeIn(const Node e_in, CutterState<Type>& cs, const bool pushFlow) {
            assert(hg.is_edge_in(e_in));
            assert(pushFlow || hg.excess(e_in) > 0);

            size_t minLevel = hg.numLawlerNodes() * 2;
            const Hyperedge e = hg.edgeFromLawlerNode(e_in);
            PinIterator pin_it = pushFlow ? std::max(current_pin_e_in[e], hg.beginPins(e)) : hg.beginPins(e);

            for (; pin_it < hg.endPins(e); ++pin_it) {
                const Pin pin = *pin_it;
                const Node v = pin.pin;
                if (!cs.n.isSourceReachable(v) || (v == piercingNode) || v == target) {
                    numEdgeScans++;
                    Flow residual = std::min(hg.flowSent(pin.he_inc_iter), hg.excess(e_in));
                    if (residual > 0) {
                        if (hg.label(e_in) == hg.label(v) + 1) {
                            assert(pushFlow);
                            numPushes++;
                            hg.push_edgeIn_to_node(v, pin.he_inc_iter, residual);
                            if (!(v == target) && !(v == piercingNode) && !inQueue[v].exchange(true)) { assert(hg.excess(v) > 0); nodes.push_back(v); }
                        } else {
                            assert(hg.label(e_in) <= hg.label(v) || v == piercingNode);
                            minLevel = std::min(minLevel, hg.label(v));
                        }
                    }
                }
                if (hg.excess(e_in) == 0) break;
            }

            if (hg.excess(e_in) > 0 && pin_it == hg.endPins(e)) {
                numEdgeScans++;
                const Node e_out = hg.edge_node_out(e);
                Flow residual = std::min(hg.capacity(e) - hg.flow(e), hg.excess(e_in));
                if (residual > 0) {
                    if (hg.label(e_in) == hg.label(e_out) + 1) {
                        assert(pushFlow);
                        numPushes++;
                        hg.push_edgeIn_to_edgeOut(e_in, e_out, residual);
                        if (!inQueue[e_out].exchange(true)) { assert(hg.excess(e_out) > 0); nodes.push_back(e_out); }
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

        void iterateOverEdgeNodeOut(const Node e_out, CutterState<Type>& cs, const bool pushFlow) {
            assert(hg.excess(e_out) > 0);

            size_t minLevel = hg.numLawlerNodes() * 2;
            const Hyperedge e = hg.edgeFromLawlerNode(e_out);
            PinIterator pin_it = pushFlow ? std::max(current_pin_e_out[e], hg.beginPins(e)) : hg.beginPins(e);

            for (; pin_it < hg.endPins(e); ++pin_it) {
                const Pin pin = *pin_it;
                const Node v = pin.pin;
                if (!cs.n.isSourceReachable(v) || (v == piercingNode) || (v == target)) {
                    numEdgeScans++;
                    Flow residual = hg.excess(e_out);
                    if (residual > 0) {
                        if (hg.label(e_out) == hg.label(v) + 1) {
                            assert(pushFlow);
                            numPushes++;
                            hg.push_edgeOut_to_node(v, pin.he_inc_iter, residual);
                            if (!(v == target) && !(v == piercingNode) && !inQueue[v].exchange(true)) { assert(hg.excess(v) > 0); nodes.push_back(v); }
                        } else {
                            assert(hg.label(e_out) <= hg.label(v) || piercingNode == v);
                            minLevel = std::min(minLevel, hg.label(v));
                        }
                    }
                }
                if (hg.excess(e_out) == 0) break;
            }

            if (hg.excess(e_out) > 0 && pin_it == hg.endPins(e)) {
                numEdgeScans++;
                const Node e_in = hg.edge_node_in(e);
                Flow residual = std::min(hg.flow(e), hg.excess(e_out));
                if (residual > 0) {
                    if (hg.label(e_out) == hg.label(e_in) + 1) {
                        assert(pushFlow);
                        numPushes++;
                        hg.push_edgeOut_to_edgeIn(e_out, e_in, residual);
                        if (!inQueue[e_in].exchange(true)) {
                            assert(hg.excess(e_in) > 0);
                            nodes.push_back(e_in);
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

            iterateOverNode(u, cs, true);
            if (hg.excess(u) > 0) {
                // Relabel
                numRelabel++;
                iterateOverNode(u, cs, false);
            }
        }

        void dischargeEdgeNodeIn(const Node e_in, CutterState<Type>& cs) {
            assert(hg.excess(e_in) > 0);
            assert(hg.is_edge_in(e_in));

            iterateOverEdgeNodeIn(e_in, cs, true);
            if (hg.excess(e_in) > 0) {
                // Relabel
                numRelabel++;
                iterateOverEdgeNodeIn(e_in, cs, false);
            }
        }

        void dischargeEdgeNodeOut(const Node e_out, CutterState<Type>& cs) {
            assert(hg.excess(e_out) > 0);

            iterateOverEdgeNodeOut(e_out, cs, true);
            if (hg.excess(e_out) > 0) {
                // Relabel
                numRelabel++;
                iterateOverEdgeNodeOut(e_out, cs, false);
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
            setLabels(cs);
            timer.stop("setLabels");

            target = cs.targetPiercingNodes.begin()->node;

            queue.clear();

            timer.start("mainLoop", "exhaustFlow");
            for (auto& sp : cs.sourcePiercingNodes) {
                for (InHeIndex inc_iter : hg.incidentHyperedgeIndices(sp.node)) {
                    InHe& inc_u = hg.getInHe(inc_iter);
                    const Hyperedge e = inc_u.e;
                    if (!cs.h.areAllPinsSourceReachable__unsafe__(e)) {
                        const Node e_out = hg.edge_node_out(e);
                        if (hg.label(e_out) != hg.numLawlerNodes()) {
                            Flow residual = hg.flowReceived(inc_iter);
                            if (residual > 0) {
                                hg.push_node_to_edgeOut(sp.node, inc_iter, residual);
                                if (!inQueue[e_out].exchange(true)) { assert(hg.excess(e_out) > 0); nodes.push_back(e_out); }
                            }
                        }

                        const Node e_in = hg.edge_node_in(e);
                        if (hg.label(e_in) != hg.numLawlerNodes()) {
                            Flow residual = hg.capacity(e);
                            if (residual > 0) {
                                hg.push_node_to_edgeIn(sp.node, inc_iter, residual);
                                if (!inQueue[e_in].exchange(true)) { assert(hg.excess(e_in) > 0); nodes.push_back(e_in); }
                            }
                        }
                    }
                }
                piercingNode = sp.node;
            }

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

                if (numEdgeScans > 12 * hg.numLawlerNodes() + hg.numHyperedges() + 2 * hg.numPins()) {
                    setLabels(cs);
                    numEdgeScans = 0;
                }
            }
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

        inline bool tryToGetNode(const Node u, const size_t currentLabel, const size_t n) {
            return hg.label(u) > currentLabel && __atomic_exchange_n(&hg.label(u), currentLabel, __ATOMIC_ACQ_REL) == n;
        }

        void setLabels(CutterState<Type>& cs) {
            std::vector<Node>& thisLayer = thisLayer_thread_specific->local();
            bool nodes_left = true;

            hg.equalizeLabels();

            std::fill(current_pin_e_in.begin(), current_pin_e_in.end(), PinIterator(0));
            std::fill(current_pin_e_out.begin(), current_pin_e_out.end(), PinIterator(0));
            std::fill(current_hyperedge.begin(), current_hyperedge.end(), InHeIndex::Invalid());

            // Source and target of the bfs, not of the flow instance
            const Node source = *cs.targetPiercingNodes.begin()->node;
            const Node target = *cs.sourcePiercingNodes.begin()->node;

            thisLayer.clear();
            thisLayer.push_back(source);
            hg.label(source) = 0;

            size_t currentLabel = 1;
            const size_t n = hg.numLawlerNodes();

            auto visitNode = [&](const Node u, std::vector<Node>& queue) {
                for (InHe& inc_u : hg.hyperedgesOf(u)) {
                    const Hyperedge e = inc_u.e;
                    const Node e_in = hg.edge_node_in(e);
                    const Node e_out = hg.edge_node_out(e);
                    if (!cs.h.areAllPinsSourceReachable__unsafe__(e)) {
                        if (hg.flowSent(inc_u) > 0 && tryToGetNode(e_in, currentLabel, n)) {
                            current_pin_e_in[e] = hg.beginPins(e);
                            queue.push_back(e_in);
                        }
                        if (tryToGetNode(e_out, currentLabel, n)) {
                            current_pin_e_out[e] = hg.beginPins(e);
                            queue.push_back(e_out);
                        }
                    }
                }
            };

            auto visitEdgeIn = [&](const Node e_in, std::vector<Node>& queue) {
                const Hyperedge e = hg.edgeFromLawlerNode(e_in);
                const Node e_out = hg.edge_node_out(e);

                if (hg.flow(e) > 0 && tryToGetNode(e_out, currentLabel, n)) {
                    current_pin_e_out[e] = hg.beginPins(e);
                    queue.push_back(e_out);
                }
                for (Pin& pin : hg.pinsOf(e)) {
                    if ((!cs.n.isSourceReachable__unsafe__(pin.pin)) && tryToGetNode(pin.pin, currentLabel, n)) {
                        current_hyperedge[pin.pin] = hg.beginIndexHyperedges(pin.pin);
                        queue.push_back(pin.pin);
                    }
                }
            };

            auto visitEdgeOut = [&](const Node e_out, std::vector<Node>& queue) {
                const Hyperedge e = hg.edgeFromLawlerNode(e_out);
                const Node e_in = hg.edge_node_in(e);

                if (hg.capacity(e) - hg.flow(e) > 0 && tryToGetNode(e_in, currentLabel, n)) {
                    current_pin_e_in[e] = hg.beginPins(e);
                    queue.push_back(e_in);
                }
                for (Pin& pin : hg.pinsOf(e)) {
                    if ((!cs.n.isSourceReachable__unsafe__(pin.pin)) && (hg.flowReceived(pin.he_inc_iter) > 0) && tryToGetNode(pin.pin, currentLabel, n)) {
                        current_hyperedge[pin.pin] = hg.beginIndexHyperedges(pin.pin);
                        queue.push_back(pin.pin);
                    }
                }
            };

            while (nodes_left) {
                tbb::parallel_for_each(*thisLayer_thread_specific, [&](const std::vector<Node>& vector) {
                    tbb::parallel_for(tbb::blocked_range<size_t>(0, vector.size()), [&](const tbb::blocked_range<size_t>& nodes) {
                        std::vector<Node>& localQueue = nextLayer_thread_specific->local();
                        for (size_t i = nodes.begin(); i < nodes.end(); ++i) {
                            const Node u = vector[i];
                            if (hg.isNode(u)) {
                                visitNode(u, localQueue);
                            } else if (hg.is_edge_in(u)) {
                                visitEdgeIn(u, localQueue);
                            } else {
                                assert(hg.is_edge_out(u));
                                visitEdgeOut(u, localQueue);
                            }
                        }
                    });
                });

                nodes_left = false;
                for (std::vector<Node>& nextLayer : *nextLayer_thread_specific) {
                    if (nextLayer.size() > 0) {
                        nodes_left = true;
                        break;
                    }
                }

                std::swap(thisLayer_thread_specific, nextLayer_thread_specific);

                for (std::vector<Node>& nextLayer : *nextLayer_thread_specific) {
                    nextLayer.clear();
                }
                currentLabel++;
            }
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