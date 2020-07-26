#pragma once

#include "../hyperflowcutter/cutter_state.h"
#include "../datastructure/queue.h"
#include "../datastructure/stack.h"
#include "../datastructure/distance_reachable_sets.h"
#include "../datastructure/timestamp_set.hpp"
#include "../datastructure/copyable_atomic.h"
#include "../push_relabel/flow_hypergraph.h"
#include "../push_relabel/lawler_flow_hypergraph_parallel.h"
#include <boost/circular_buffer.hpp>


namespace whfc_pr {
    class PushRelabelParallel {
    public:
        using Type = PushRelabelParallel;
        using Node = whfc::Node;
        using InHeIndex = whfc::InHeIndex;
        using Flow = whfc::Flow;
        using ScanList = LayeredQueue<Node>;
        using TimeReporter = whfc::TimeReporter;

        using ReachableNodes = whfc::DistanceReachableNodes<FlowHypergraph>;
        using ReachableHyperedges = whfc::DistanceReachableHyperedges<FlowHypergraph>;
        using Hypergraph = LawlerFlowHypergraphParallel;
        using BaseHypergraph = FlowHypergraph;

        using Pin = FlowHypergraph::Pin;
        using InHe = FlowHypergraph::InHe;
        using PinIndexRange = FlowHypergraph::PinIndexRange;
        using PinIterator = FlowHypergraph::PinIterator;
        using DistanceT = whfc::DistanceReachableNodes<FlowHypergraph>::DistanceT;

        static constexpr size_t ALPHA = 6;
        static constexpr size_t BETA = 12;
        static constexpr float globalUpdateFreq = 10.0;

        Hypergraph& hg;
        LayeredQueue<Node> queue;
        struct StackFrame {
            Node u;
            InHeIndex he_it;
        };
        FixedCapacityStack<StackFrame> stack;

        int direction = 0;

        Flow upperFlowBound;

        static constexpr bool same_traversal_as_grow_assimilated = false;
        static constexpr bool grow_reachable_marks_flow_sending_pins_when_marking_all_pins = true;
        static constexpr bool log = false;

        TimeReporter& timer;

        Node piercingNode;
        Node target;

        size_t n;

        ldc::AtomicTimestampSet<uint16_t> inQueue;

        // for nodes
        std::vector<InHeIndex> current_hyperedge;

        // for hyperedges
        std::vector<PinIterator> current_pin_e_in, current_pin_e_out;

        size_t numPushes, numRelabel, numGlobalUpdate, workSinceLastRelabel;

        std::unique_ptr<tbb::enumerable_thread_specific<std::vector<Node>>> thisLayer_thread_specific;
        std::unique_ptr<tbb::enumerable_thread_specific<std::vector<Node>>> nextLayer_thread_specific;
        std::unique_ptr<tbb::enumerable_thread_specific<std::vector<Node>>> queue_for_global_update;

        tbb::enumerable_thread_specific<size_t> work_added_thread_specific;

        PushRelabelParallel(Hypergraph& hg, TimeReporter& timer, size_t numThreads) : hg(hg), stack(hg.maxNumLawlerNodes()), timer(timer), inQueue(hg.maxNumLawlerNodes()),
            current_hyperedge(hg.maxNumNodes, InHeIndex::Invalid()), current_pin_e_in(hg.maxNumHyperedges), current_pin_e_out(hg.maxNumHyperedges),
            thisLayer_thread_specific(new tbb::enumerable_thread_specific<std::vector<Node>>()), nextLayer_thread_specific(new tbb::enumerable_thread_specific<std::vector<Node>>()),
            queue_for_global_update(new tbb::enumerable_thread_specific<std::vector<Node>>()), work_added_thread_specific(0)
        {
            reset();
        }

        void reset() {
            queue.clear();
            numPushes = 0;
            numRelabel = 0;
            numGlobalUpdate = 0;
            workSinceLastRelabel = 0;
        }

        ScanList& getScanList() {
            return queue;
        }

        Flow recycleDatastructuresFromGrowReachablePhase(whfc::CutterState<Type> &cs) {
            throw std::logic_error("Not implemented");
        }

        Flow growFlowOrSourceReachable(whfc::CutterState<Type>& cs) {
            throw std::logic_error("Not implemented");
        }

        inline bool win_edge(const Node u, const Node v) {
            return (hg.label(u) == hg.label(v) + 1) || (hg.label(u) < hg.label(v) - 1) || ((hg.label(u) == hg.label(v)) && u < v);
        }

        bool iterateOverEdges(const Node u, whfc::CutterState<Type>& cs, std::vector<Node>& queue, size_t& minLevel) {
            bool skipped = false;
            for (InHeIndex he = hg.beginIndexHyperedges(u); he < hg.endIndexHyperedges(u); ++he) {
                InHe& inc_u = hg.getInHe(he);
                const Hyperedge e = inc_u.e;
                if (!cs.h.areAllPinsSourceReachable__unsafe__(e)) {
                    const Node e_in = hg.edge_node_in(e);
                    if (hg.excess(u) == 0) break;
                    bool admissible = hg.label_next_iteration(u) == hg.label(e_in) + 1;
                    if (hg.excess(e_in)) {
                        bool win = win_edge(u, e_in);
                        assert(!win || !win_edge(e_in, u));
                        if (admissible && !win) {
                            skipped = true;
                            continue;
                        }
                    }
                    Flow residual = hg.excess(u);
                    if (admissible && residual > 0) {
                        numPushes++;
                        hg.push_node_to_edgeIn(u, he, residual);
                        if (inQueue.set(e_in)) { queue.push_back(e_in); }
                    }
                    residual = hg.excess(u);
                    if (residual > 0 && hg.label(e_in) >= hg.label_next_iteration(u)) {
                        minLevel = std::min(minLevel, hg.label(e_in));
                    }
                }
            }

            if (hg.excess(u) == 0 || skipped) {
                return false;
            }

            for (InHeIndex he = hg.beginIndexHyperedges(u); he < hg.endIndexHyperedges(u); ++he) {
                InHe& inc_u = hg.getInHe(he);
                const Hyperedge e = inc_u.e;
                if (!cs.h.areAllPinsSourceReachable__unsafe__(e)) {
                    const Node e_out = hg.edge_node_out(e);
                    if (hg.excess(u) == 0) break;
                    bool admissible = hg.label_next_iteration(u) == hg.label(e_out) + 1;
                    if (hg.excess(e_out)) {
                        bool win = win_edge(u, e_out);
                        assert(!win || !win_edge(e_out, u));
                        if (admissible && !win) {
                            skipped = true;
                            continue;
                        }
                    }
                    Flow residual = std::min(hg.excess(u), hg.getPinOut(inc_u).flow);
                    if (admissible && residual > 0) {
                        numPushes++;
                        hg.push_node_to_edgeOut(u, he, residual);
                        if (inQueue.set(e_out)) { queue.push_back(e_out); }
                    }
                    residual = std::min(hg.excess(u), hg.getPinOut(inc_u).flow);
                    if (residual > 0 && hg.label(e_out) >= hg.label_next_iteration(u)) {
                        minLevel = std::min(minLevel, hg.label(e_out));
                    }
                }
            }

            if (hg.excess(u) == 0 || skipped) {
                return false;
            }
            return true;
        }

        bool iterateOverPins(const Node u, whfc::CutterState<Type>& cs, std::vector<Node>& queue, size_t& minLevel) {
            const bool is_edge_in = hg.is_edge_in(u);

            bool skipped = false;

            const Hyperedge e = hg.edgeFromLawlerNode(u);
            for (PinIterator pinIt = hg.beginPinsIn(e); pinIt < hg.endPinsIn(e); ++pinIt) {
                const Pin pin = *pinIt;
                const Node v = pin.pin;
                if (!cs.n.isSourceReachable(v) || v == target) {
                    if (hg.excess(u) == 0) break;
                    bool admissible = hg.label_next_iteration(u) == hg.label(v) + 1;
                    if (hg.excess(v) > 0) {
                        bool win = win_edge(u, v);
                        assert(!win || !win_edge(v, u));
                        if (admissible && !win) {
                            skipped = true;
                            continue;
                        }
                    }
                    Flow residual = is_edge_in ? std::min(pin.flow, hg.excess(u)) : hg.excess(u);
                    if (admissible && residual > 0) {
                        numPushes++;
                        if (is_edge_in) {
                            hg.push_edgeIn_to_node(v, pin.he_inc_iter, residual);
                        } else {
                            hg.push_edgeOut_to_node(v, pin.he_inc_iter, residual);
                        }
                        if (!(v == target) && inQueue.set(v)) { queue.push_back(v); }
                    }
                    residual = is_edge_in ? std::min(pin.flow, hg.excess(u)) : hg.excess(u);
                    if (residual > 0 && hg.label(v) >= hg.label_next_iteration(u)) {
                        minLevel = std::min(minLevel, hg.label(v));
                    }
                }
            }

            if (hg.excess(u) == 0 || skipped) {
                return false;
            }
            return true;
        }

        void dischargeNode(const Node u, whfc::CutterState<Type>& cs, std::vector<Node>& queue, size_t& work_added) {
            assert(u < hg.numNodes());
            assert(hg.excess(u) > 0);

            while (hg.excess(u) > 0) {
                size_t minLevel = n;
                bool should_relabel = iterateOverEdges(u, cs, queue, minLevel);

                if (should_relabel) {
                    assert(minLevel + 1 > hg.label(u));
                    hg.label_next_iteration(u) = minLevel + 1;
                    work_added += (hg.endIndexHyperedges(u) - hg.beginIndexHyperedges(u)) * 2 + BETA;
                } else {
                    break;
                }
                if (hg.label_next_iteration(u) >= n) break;
            }
            if (hg.excess(u) > 0 && hg.label_next_iteration(u) < n && inQueue.set(u)) {
                queue.push_back(u);
            }

        }

        void dischargeEdgeNodeIn(const Node e_in, whfc::CutterState<Type>& cs, std::vector<Node>& queue, size_t& work_added) {
            assert(hg.excess(e_in) > 0);
            assert(hg.is_edge_in(e_in));

            const Hyperedge e = hg.edgeFromLawlerNode(e_in);
            const Node e_out = hg.edge_node_out(e);
            bool should_relabel = true;

            while (hg.excess(e_in) > 0) {
                size_t minLevel = n;

                bool admissible = hg.label_next_iteration(e_in) == hg.label(e_out) + 1;
                if (hg.excess(e_out) > 0) {
                    bool win = win_edge(e_in, e_out);
                    assert(!win || !win_edge(e_out, e_in));
                    if (admissible && !win) {
                        should_relabel = false;
                    }
                }
                Flow residual = std::min(hg.capacity(e) - hg.flow(e), hg.excess(e_in));
                if (should_relabel && admissible && residual > 0) {
                    numPushes++;
                    hg.push_edgeIn_to_edgeOut(e_in, e_out, residual);
                    if (hg.excess(e_in) == 0) should_relabel = false;
                    if (inQueue.set(e_out)) { queue.push_back(e_out); }
                }
                residual = std::min(hg.capacity(e) - hg.flow(e), hg.excess(e_in));
                if (should_relabel && residual > 0 && hg.label(e_out) >= hg.label_next_iteration(e_in)) {
                    minLevel = std::min(minLevel, hg.label(e_out));
                }

                should_relabel = should_relabel && iterateOverPins(e_in, cs, queue, minLevel);

                if (should_relabel) {
                    assert(minLevel + 1 > hg.label_next_iteration(e_in));
                    hg.label_next_iteration(e_in) = minLevel + 1;
                    work_added += (hg.endPinsIn(e) - hg.beginPinsIn(e)) + BETA;
                } else {
                    break;
                }
                if (hg.label_next_iteration(e_in) >= n) break;
            }

            if (hg.excess(e_in) > 0 && hg.label_next_iteration(e_in) < n && inQueue.set(e_in)) {
                queue.push_back(e_in);
            }

        }

        void dischargeEdgeNodeOut(const Node e_out, whfc::CutterState<Type>& cs, std::vector<Node>& queue, size_t& work_added) {
            assert(hg.excess(e_out) > 0);
            assert(hg.is_edge_out(e_out));

            const Hyperedge e = hg.edgeFromLawlerNode(e_out);
            const Node e_in = hg.edge_node_in(e);
            bool should_relabel = true;

            while (hg.excess(e_out) > 0) {
                size_t minLevel = n;

                should_relabel = should_relabel && iterateOverPins(e_out, cs, queue, minLevel);

                bool skipped = false;
                bool admissible = hg.label_next_iteration(e_out) == hg.label(e_in) + 1;
                if (hg.excess(e_in) > 0) {
                    bool win = win_edge(e_out, e_in);
                    assert(!win || !win_edge(e_in, e_out));
                    if (admissible && !win) {
                        skipped = true;
                        should_relabel = false;
                    }
                }
                Flow residual = std::min(hg.flow(e), hg.excess(e_out));
                if (!skipped && admissible && residual > 0) {
                    numPushes++;
                    hg.push_edgeOut_to_edgeIn(e_out, e_in, residual);
                    if (hg.excess(e_out) == 0) should_relabel = false;
                    if (inQueue.set(e_in)) queue.push_back(e_in);
                }
                residual = std::min(hg.flow(e), hg.excess(e_out));
                if (!skipped && residual > 0 && hg.label(e_in) >= hg.label_next_iteration(e_out)) {
                    minLevel = std::min(minLevel, hg.label(e_in));
                }

                if (should_relabel) {
                    assert(minLevel + 1 > hg.label_next_iteration(e_out));
                    hg.label_next_iteration(e_out) = minLevel + 1;
                    work_added += (hg.endPinsIn(e) - hg.beginPinsIn(e)) + BETA;
                } else {
                    break;
                }
                if (hg.label_next_iteration(e_out) >= n) break;
            }

            if (hg.excess(e_out) > 0 && hg.label_next_iteration(e_out) < n && inQueue.set(e_out)) {
                queue.push_back(e_out);
            }

        }

        void pushRelabel(whfc::CutterState<Type>& cs, const size_t nm, const bool keepNodesWithHighLabel) {
            bool nodes_left = true;

            while (nodes_left) {
                inQueue.reset();

                if (workSinceLastRelabel * globalUpdateFreq > nm) {
                    timer.start("globalUpdate", "phase1");
                    globalUpdate(cs, n);
                    numGlobalUpdate++;
                    timer.stop("globalUpdate");
                    workSinceLastRelabel = 0;
                }

                timer.start("discharge", "phase1");
                tbb::parallel_for_each(*thisLayer_thread_specific, [&](const std::vector<Node>& vector) {
                    tbb::parallel_for(tbb::blocked_range<size_t>(0, vector.size()), [&](const tbb::blocked_range<size_t>& nodes) {
                        std::vector<Node>& nextLayer = nextLayer_thread_specific->local();
                        size_t& work_added = work_added_thread_specific.local();
                        for (size_t i = nodes.begin(); i < nodes.end(); ++i) {
                            const Node u = vector[i];
                            assert(hg.excess(u) > 0);
                            hg.label_next_iteration(u) = hg.label(u);
                            if (hg.label(u) >= n) continue;

                            if (hg.isNode(u)) {
                                dischargeNode(u, cs, nextLayer, work_added);
                            } else if (hg.is_edge_in(u)) {
                                dischargeEdgeNodeIn(u, cs, nextLayer, work_added);
                            } else {
                                assert(hg.is_edge_out(u));
                                dischargeEdgeNodeOut(u, cs, nextLayer, work_added);
                            }
                        }
                    });
                });
                timer.stop("discharge");

                timer.start("set_excess_and_label", "phase1");
                tbb::parallel_for_each(*thisLayer_thread_specific, [&](const std::vector<Node>& vector) {
                    tbb::parallel_for(tbb::blocked_range<size_t>(0, vector.size()), [&](const tbb::blocked_range<size_t>& nodes) {
                        for (size_t i = nodes.begin(); i < nodes.end(); ++i) {
                            const Node u = vector[i];
                            hg.excess(u) += hg.excess_change(u);
                            hg.excess_change(u) = 0;
                            hg.label(u) = hg.label_next_iteration(u);
                        }
                    });
                });
                timer.stop("set_excess_and_label");

                std::swap(thisLayer_thread_specific, nextLayer_thread_specific);

                timer.start("set_excess", "phase1");
                tbb::parallel_for_each(*thisLayer_thread_specific, [&](const std::vector<Node>& vector) {
                    tbb::parallel_for(tbb::blocked_range<size_t>(0, vector.size()), [&](const tbb::blocked_range<size_t>& nodes) {
                        for (size_t i = nodes.begin(); i < nodes.end(); ++i) {
                            const Node u = vector[i];
                            hg.excess(u) += hg.excess_change(u);
                            hg.excess_change(u) = 0;
                        }
                    });
                });
                timer.stop("set_excess");

                nodes_left = false;
                for (std::vector<Node>& vector : *thisLayer_thread_specific) {
                    if (vector.size() > 0) nodes_left = true;
                }

                for (std::vector<Node>& vector : *nextLayer_thread_specific) {
                    vector.clear();
                }

                for (size_t& work_added : work_added_thread_specific) {
                    workSinceLastRelabel += work_added;
                    work_added = 0;
                }

            }

            hg.excess(target) += hg.excess_change(target);
            hg.excess_change(target) = 0;
        }

        Flow exhaustFlow(whfc::CutterState<Type>& cs) {
            assert(cs.sourcePiercingNodes.size() == 1);
            timer.start("exhaustFlow");

            hg.alignViewDirection();

            cs.clearForSearch();

            timer.start("initialize", "exhaustFlow");
            hg.initialize_for_push_relabel();
            timer.stop("initialize");

            piercingNode = cs.sourcePiercingNodes.begin()->node;
            target = cs.targetPiercingNodes.begin()->node;

            std::vector<Node>& thisLayer = thisLayer_thread_specific->local();
            thisLayer.clear();

            inQueue.reset();

            bool nodes_left = false;

            for (InHeIndex inc_iter : hg.incidentHyperedgeIndices(piercingNode)) {
                InHe& inc_u = hg.getInHe(inc_iter);
                const Hyperedge e = inc_u.e;
                if (!cs.h.areAllPinsSourceReachable__unsafe__(e)) {
                    const Node e_out = hg.edge_node_out(e);
                    Flow residual = hg.getPinOut(hg.getInHe(inc_iter)).flow;
                    if (residual > 0) {
                        hg.push_node_to_edgeOut(piercingNode, inc_iter, residual);
                        nodes_left = true;
                        if (inQueue.set(e_out)) { thisLayer.push_back(e_out); }
                    }

                    const Node e_in = hg.edge_node_in(e);
                    residual = hg.capacity(e);
                    if (residual > 0) {
                        hg.push_node_to_edgeIn(piercingNode, inc_iter, residual);
                        nodes_left = true;
                        if (inQueue.set(e_in)) { thisLayer.push_back(e_in); }
                    }
                }
            }

            for (Node u : thisLayer) {
                hg.excess(u) += hg.excess_change(u);
                hg.excess_change(u) = 0;
            }

            timer.start("setLabels", "exhaustFlow");
            n = hg.numNodes() + 2 * hg.numHyperedges();
            std::array<size_t, 2> numScans = globalUpdate(cs, n);
            timer.stop("setLabels");

            size_t nm = ALPHA * numScans[0] + numScans[1];
            n = numScans[0];

            timer.start("mainLoop", "exhaustFlow");

            timer.start("phase1", "mainLoop");
            // Phase 1
            if (nodes_left) pushRelabel(cs, nm, false);
            timer.stop("phase1");

            assert(hg.excess_sums_to_zero());

            timer.start("phase2", "mainLoop");
            timer.start("buildResidualNetwork", "phase2");
            hg.buildResidualNetwork();
            timer.stop("buildResidualNetwork");

            Flow augmented_flow = flowDecomposition(cs);
            whfc::unused(augmented_flow);

            timer.stop("phase2");
            timer.stop("mainLoop");

            timer.start("writeBackFlow", "exhaustFlow");
            hg.writeBackFlow();
            timer.stop("writeBackFlow");

            resetSourcePiercingNodeDistances(cs);

            Flow f = 0;

            for (auto& sp : cs.targetPiercingNodes) {
                f += hg.excess(sp.node);
            }
            cs.flowValue += f;
            assert(f == augmented_flow);


            cs.verifyFlowConstraints();

            timer.start("growReachable", "exhaustFlow");
            bool found_target = growReachable(cs);
            whfc::unused(found_target);
            assert(!found_target);
            timer.stop("growReachable");
            timer.stop("exhaustFlow");

            return f;
        }

        Flow flowDecomposition(whfc::CutterState<Type>& cs) {
            Flow total_augmented_flow = 0;

            timer.start("resets", "phase2");
            inQueue.reset();

            std::fill(current_hyperedge.begin(), current_hyperedge.end(), InHeIndex::Invalid());
            std::fill(current_pin_e_in.begin(), current_pin_e_in.end(), PinIterator(0));
            std::fill(current_pin_e_out.begin(), current_pin_e_out.end(), PinIterator(0));
            timer.stop("resets");

            const Node source = *cs.sourcePiercingNodes.begin()->node;
            const Node target = *cs.targetPiercingNodes.begin()->node;

            timer.start("dfs", "phase2");
            assert(stack.empty());
            stack.push({ source, InHeIndex::Invalid() });

            while (!stack.empty()) {
                const Node u = stack.top().u;
                bool found_new_node = false;
                if (hg.isNode(u)) {
                    InHeIndex& he_it = current_hyperedge[u];
                    if (he_it == InHeIndex::Invalid()) he_it = hg.beginIndexHyperedges(u);
                    for (; he_it < hg.endIndexHyperedges(u); ++he_it) {
                        InHe& inc_u = hg.getInHe(he_it);
                        const Hyperedge e = inc_u.e;
                        const Flow residual_to_edge_in = hg.getPinIn(inc_u).flow;
                        const Flow residual_to_edge_out = -hg.getPinOut(inc_u).flow;

                        if (residual_to_edge_in > 0) {
                            stack.push({ hg.edge_node_in(e), he_it });
                            found_new_node = true;
                            break;
                        } else if (residual_to_edge_out > 0) {
                            stack.push({ hg.edge_node_out(e), he_it });
                            found_new_node = true;
                            break;
                        }

                    }
                } else if (hg.is_edge_in(u)) {
                    const Hyperedge e = hg.edgeFromLawlerNode(u);
                    if (hg.flow(e) > 0) {
                        stack.push( { hg.edge_node_out(e), InHeIndex::Invalid() });
                        found_new_node = true;
                    }
                    if (!found_new_node) {
                        PinIterator& pin_it = current_pin_e_in[e];
                        if (pin_it == PinIterator(0)) pin_it = hg.beginPinsIn(e);
                        for ( ; pin_it < hg.endPinsIn(e); ++pin_it) {
                            const Pin pin = *pin_it;
                            const Node v = pin.pin;
                            const Flow residual_to_v = -hg.getPinIn(hg.getInHe(pin.he_inc_iter)).flow;
                            if (residual_to_v > 0) {
                                stack.push( { v, pin.he_inc_iter });
                                found_new_node = true;
                                break;
                            }
                        }
                    }
                } else {
                    assert(hg.is_edge_out(u));
                    const Hyperedge e = hg.edgeFromLawlerNode(u);
                    if (hg.flow(e) < 0) {
                        stack.push( { hg.edge_node_in(e), InHeIndex::Invalid() });
                        found_new_node = true;
                    }
                    if (!found_new_node) {
                        PinIterator& pin_it = current_pin_e_out[e];
                        if (pin_it == PinIterator(0)) pin_it = hg.beginPinsIn(e);
                        for ( ; pin_it < hg.endPinsIn(e); ++pin_it) {
                            const Pin pin = *pin_it;
                            const Node v = pin.pin;
                            const Flow residual_to_v = hg.getPinOut(hg.getInHe(pin.he_inc_iter)).flow;
                            if (residual_to_v > 0) {
                                stack.push( { v, pin.he_inc_iter });
                                found_new_node = true;
                                break;
                            }
                        }
                    }

                }

                if (found_new_node) {
                    if (stack.top().u == target) {
                        //std::cout << "Push back from target" << std::endl;
                        total_augmented_flow += pushBack(maxFlow, source, true);
                    } else if (hg.excess(stack.top().u) > 0) {
                        //std::cout << "Push back excess to source" << std::endl;
                        const Node top = stack.top().u;
                        Flow flow_pushed_back = pushBack(hg.excess(stack.top().u), source, false);
                        hg.excess(top) -= flow_pushed_back;
                    } else if (inQueue.isSet(stack.top().u)) {
                        //std::cout << "Push back looped flow" << std::endl;
                        pushBack(maxFlow, stack.top().u, false);
                    } else {
                        inQueue.set(stack.top().u);
                    }
                } else {
                    inQueue.reset(stack.top().u);
                    stack.pop();
                }

            }
            timer.stop("dfs");
            return total_augmented_flow;
        }

        Flow pushBack(Flow bottleneckCapacity, Node first_node, bool writeFlowToResult) {
            int64_t lowest_bottleneck = stack.size() - 1;
            for (int64_t stack_pointer = stack.size() - 1; stack_pointer > 0; --stack_pointer) {
                const StackFrame& t = stack.at(stack_pointer);
                Flow residual = 0;
                if (hg.isNode(t.u)) {
                    if (hg.is_edge_out(stack.at(stack_pointer-1).u)) {
                        residual = hg.getPinOut(hg.getInHe(t.he_it)).flow;
                        assert(residual > 0);
                    } else {
                        residual = -hg.getPinIn(hg.getInHe(t.he_it)).flow;
                        assert(residual > 0);
                    }
                } else if (hg.is_edge_in(t.u)) {
                    if (t.he_it != InHeIndex::Invalid()) {
                        residual = hg.getPinIn(hg.getInHe(t.he_it)).flow;
                        assert(residual > 0);
                    } else {
                        residual = -hg.flow(hg.edgeFromLawlerNode(t.u));
                        assert(residual > 0);
                    }
                } else {
                    if (t.he_it != InHeIndex::Invalid()) {
                        residual = -hg.getPinOut(hg.getInHe(t.he_it)).flow;
                        assert(residual > 0);
                    } else {
                        residual = hg.flow(hg.edgeFromLawlerNode(t.u));
                        assert(residual > 0);
                    }
                }
                if (residual <= bottleneckCapacity) {
                    bottleneckCapacity = residual;
                    lowest_bottleneck = stack_pointer;
                }
                if (stack.at(stack_pointer - 1).u == first_node) break;
            }
            assert(bottleneckCapacity > 0);
            //std::cout << "Push back " << bottleneckCapacity << std::endl;

            for (int64_t stack_pointer = stack.size() - 1; stack_pointer > 0; --stack_pointer) {
                const StackFrame& t = stack.at(stack_pointer);
                //std::cout << t.u << " ";
                if (hg.isNode(t.u)) {
                    if (hg.is_edge_out(stack.at(stack_pointer-1).u)) {
                        assert(hg.getPinOut(hg.getInHe(t.he_it)).flow >= bottleneckCapacity);
                        hg.getPinOut(hg.getInHe(t.he_it)).flow -= bottleneckCapacity;
                        if (writeFlowToResult) hg.getInHe(t.he_it).flow -= hg.flowSent(bottleneckCapacity);
                    } else {
                        assert(hg.getPinIn(hg.getInHe(t.he_it)).flow <= -bottleneckCapacity);
                        hg.getPinIn(hg.getInHe(t.he_it)).flow += bottleneckCapacity;
                        if (writeFlowToResult) hg.getInHe(t.he_it).flow -= hg.flowSent(bottleneckCapacity);
                    }
                } else if (hg.is_edge_in(t.u)) {
                    if (t.he_it != InHeIndex::Invalid()) {
                        assert(hg.getPinIn(hg.getInHe(t.he_it)).flow >= bottleneckCapacity);
                        hg.getPinIn(hg.getInHe(t.he_it)).flow -= bottleneckCapacity;
                        if (writeFlowToResult) hg.getInHe(t.he_it).flow += hg.flowSent(bottleneckCapacity);
                    } else {
                        assert(hg.flow(hg.edgeFromLawlerNode(t.u)) <= -bottleneckCapacity);
                        hg.flow(hg.edgeFromLawlerNode(t.u)) += bottleneckCapacity;
                    }
                } else {
                    if (t.he_it != InHeIndex::Invalid()) {
                        assert(hg.getPinOut(hg.getInHe(t.he_it)).flow <= -bottleneckCapacity);
                        hg.getPinOut(hg.getInHe(t.he_it)).flow += bottleneckCapacity;
                        if (writeFlowToResult) hg.getInHe(t.he_it).flow += hg.flowSent(bottleneckCapacity);
                    } else {
                        assert(hg.flow(hg.edgeFromLawlerNode(t.u)) >= bottleneckCapacity);
                        hg.flow(hg.edgeFromLawlerNode(t.u)) -= bottleneckCapacity;
                    }
                }
                if (stack_pointer >= lowest_bottleneck) inQueue.reset(t.u);
                if (stack.at(stack_pointer - 1).u == first_node) break;
            }
            //std::cout << std::endl;
            stack.popDownTo(lowest_bottleneck-1);
            return bottleneckCapacity;
        }

        inline bool tryToGetNode(const Node u, const size_t currentLabel, const size_t n) {
            return hg.label(u) > currentLabel && __atomic_exchange_n(&hg.label(u), currentLabel, __ATOMIC_ACQ_REL) == n;
        }

        std::array<size_t, 2> globalUpdate(whfc::CutterState<Type>& cs, size_t n) {
            size_t scannedNodes = 0;
            size_t scannedEdges = 0;
            tbb::enumerable_thread_specific<std::array<size_t, 2>> scanned_nodes_and_edges(std::array<size_t, 2> {0, 0});

            for (std::vector<Node>& vector : *nextLayer_thread_specific) {
                vector.clear();
            }
            for (std::vector<Node>& vector : *queue_for_global_update) {
                vector.clear();
            }
            for (std::vector<Node>& vector : *thisLayer_thread_specific) {
                vector.clear();
            }

            std::vector<Node>& local_queue = queue_for_global_update->local();

            hg.equalizeLabels(n);

            const Node source = *cs.targetPiercingNodes.begin()->node;

            local_queue.push_back(source);
            hg.label(source) = 0;

            size_t currentLabel = 1;

            auto visitNode = [&](const Node u, std::vector<Node>& queue, std::vector<Node>& push_relabel_vector, size_t& scannedEdges) {
                for (InHe& inc_u : hg.hyperedgesOf(u)) {
                    scannedEdges += 2;
                    const Hyperedge e = inc_u.e;
                    const Node e_in = hg.edge_node_in(e);
                    const Node e_out = hg.edge_node_out(e);
                    if (!cs.h.areAllPinsSourceReachable__unsafe__(e)) {
                        if (hg.getPinIn(inc_u).flow > 0 && tryToGetNode(e_in, currentLabel, n)) {
                            hg.label(e_in) = currentLabel;
                            current_pin_e_in[e] = hg.beginPinsIn(e);
                            queue.push_back(e_in);
                            if (hg.excess(e_in) > 0) push_relabel_vector.push_back(e_in);
                        }
                        if (tryToGetNode(e_out, currentLabel, n)) {
                            hg.label(e_out) = currentLabel;
                            current_pin_e_out[e] = hg.beginPinsIn(e);
                            queue.push_back(e_out);
                            if (hg.excess(e_out) > 0) push_relabel_vector.push_back(e_out);
                        }
                    }
                }
            };


            auto visitEdgeIn = [&](const Node e_in, std::vector<Node>& queue, std::vector<Node>& push_relabel_vector, size_t& scannedEdges) {
                const Hyperedge e = hg.edgeFromLawlerNode(e_in);
                const Node e_out = hg.edge_node_out(e);

                scannedEdges++;
                if (hg.flow(e) > 0 && tryToGetNode(e_out, currentLabel, n)) {
                    hg.label(e_out) = currentLabel;
                    current_pin_e_out[e] = hg.beginPinsIn(e);
                    queue.push_back(e_out);
                    if (hg.excess(e_out) > 0) push_relabel_vector.push_back(e_out);
                }
                for (Pin& pin : hg.pinsInOf(e)) {
                    scannedEdges++;
                    if ((!cs.n.isSourceReachable__unsafe__(pin.pin)) && tryToGetNode(pin.pin, currentLabel, n)) {
                        hg.label(pin.pin) = currentLabel;
                        current_hyperedge[pin.pin] = hg.beginIndexHyperedges(pin.pin);
                        queue.push_back(pin.pin);
                        if (hg.excess(pin.pin) > 0) push_relabel_vector.push_back(pin.pin);
                    }
                }
            };

            auto visitEdgeOut = [&](const Node e_out, std::vector<Node>& queue, std::vector<Node>& push_relabel_vector, size_t& scannedEdges) {
                const Hyperedge e = hg.edgeFromLawlerNode(e_out);
                const Node e_in = hg.edge_node_in(e);

                scannedEdges++;
                if (hg.capacity(e) - hg.flow(e) > 0 && tryToGetNode(e_in, currentLabel, n)) {
                    hg.label(e_in) = currentLabel;
                    current_pin_e_in[e] = hg.beginPinsIn(e);
                    queue.push_back(e_in);
                    if (hg.excess(e_in) > 0) push_relabel_vector.push_back(e_in);
                }
                for (Pin& pin : hg.pinsInOf(e)) {
                    scannedEdges++;
                    if ((hg.getPinOut(hg.getInHe(pin.he_inc_iter)).flow > 0) && (!cs.n.isSourceReachable__unsafe__(pin.pin)) && tryToGetNode(pin.pin, currentLabel, n)) {
                        hg.label(pin.pin) = currentLabel;
                        current_hyperedge[pin.pin] = hg.beginIndexHyperedges(pin.pin);
                        queue.push_back(pin.pin);
                        if (hg.excess(pin.pin) > 0) push_relabel_vector.push_back(pin.pin);
                    }
                }
            };

            bool nodes_left = true;

            while (nodes_left) {
                tbb::parallel_for_each(*queue_for_global_update, [&](const std::vector<Node>& vector) {
                    tbb::parallel_for(tbb::blocked_range<size_t>(0, vector.size()), [&](const tbb::blocked_range<size_t>& nodes) {
                        std::vector<Node>& localQueue = nextLayer_thread_specific->local();
                        std::vector<Node>& queue_for_push_relabel = thisLayer_thread_specific->local();
                        std::array<size_t, 2>& scanned_sizes = scanned_nodes_and_edges.local();
                        scanned_sizes[0] += nodes.size();
                        for (size_t i = nodes.begin(); i < nodes.end(); ++i) {
                            const Node u = vector[i];
                            if (hg.isNode(u)) {
                                visitNode(u, localQueue, queue_for_push_relabel, scanned_sizes[1]);
                            } else if (hg.is_edge_in(u)) {
                                visitEdgeIn(u, localQueue, queue_for_push_relabel, scanned_sizes[1]);
                            } else {
                                assert(hg.is_edge_out(u));
                                visitEdgeOut(u, localQueue, queue_for_push_relabel, scanned_sizes[1]);
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

                std::swap(queue_for_global_update, nextLayer_thread_specific);

                for (std::vector<Node>& nextLayer : *nextLayer_thread_specific) {
                    nextLayer.clear();
                }
                currentLabel++;
            }

            for (std::vector<Node>& vector : *nextLayer_thread_specific) {
                vector.clear();
            }

            for (std::array<size_t, 2>& sizes : scanned_nodes_and_edges) {
                scannedNodes += sizes[0];
                scannedEdges += sizes[1];
            }

            return {scannedNodes, scannedEdges};
        }

        bool growReachable(whfc::CutterState<Type>& cs) {
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
                                for (const Pin& pv : hg.pinsInOf(e)) {
                                    if (hg.flowSent(pv) > 0 && scanFlowSending) {
                                        visit(pv, true);
                                    } else if (hg.flowSent(pv) <= 0 && scanAllPins) {
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
            cs.verifySetInvariants();
            return found_target;
        }

        void resetSourcePiercingNodeDistances(whfc::CutterState<Type>& cs, bool reset = true) {
            for (auto& sp: cs.sourcePiercingNodes)
                cs.n.setPiercingNodeDistance(sp.node, reset);
        }

    };
}