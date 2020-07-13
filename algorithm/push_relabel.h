
#include "cutter_state.h"
#include "../datastructure/queue.h"
#include "../datastructure/stack.h"
#include "../datastructure/distance_reachable_sets.h"
//#include "ford_fulkerson.h"
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
        static constexpr float globalUpdateFreq = 1.0;

        LawlerFlowHypergraph& hg;
        LayeredQueue<Node> queue;
        boost::circular_buffer<Node> nodes;
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



        PushRelabel(LawlerFlowHypergraph& hg, TimeReporter& timer, size_t numThreads) : hg(hg), nodes(hg.maxNumLawlerNodes()), stack(hg.maxNumLawlerNodes()), timer(timer), inQueue(hg.maxNumLawlerNodes()),
            current_hyperedge(hg.maxNumNodes, InHeIndex::Invalid()), current_pin_e_in(hg.maxNumHyperedges), current_pin_e_out(hg.maxNumHyperedges)
        {
            reset();
        }

        void reset() {
            queue.clear();
            nodes.clear();
            numPushes = 0;
            numRelabel = 0;
            numGlobalUpdate = 0;
            workSinceLastRelabel = 0;
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

        size_t iterateOverEdges(const Node u, CutterState<Type>& cs, InHeIndex& he, const InHeIndex he_end) {
            size_t minLevel = n;
            for (; he < hg.endIndexHyperedges(u); ++he) {
                InHe& inc_u = hg.getInHe(he);
                const Hyperedge e = inc_u.e;

                if (!cs.h.areAllPinsSourceReachable__unsafe__(e)) {
                    const Node e_out = hg.edge_node_out(e);
                    Flow residual = std::min(hg.excess(u), hg.getPinOut(inc_u).flow);
                    if (residual > 0) {
                        if (hg.label(u) == hg.label(e_out) + 1) {
                            numPushes++;
                            hg.push_node_to_edgeOut(u, he, residual);
                            if (inQueue.set(e_out)) { nodes.push_back(e_out); }
                        } else {
                            assert(hg.label(u) <= hg.label(e_out));
                            if (hg.label(e_out) != 0) minLevel = std::min(minLevel, hg.label(e_out));
                        }
                    }

                    const Node e_in = hg.edge_node_in(e);
                    residual = hg.excess(u);
                    if (residual > 0) {
                        if (hg.label(u) == hg.label(e_in) + 1) {
                            numPushes++;
                            hg.push_node_to_edgeIn(u, he, residual);
                            if (inQueue.set(e_in)) { nodes.push_back(e_in); }
                        } else {
                            assert(hg.label(u) <= hg.label(e_in));
                            if (hg.label(e_in) != 0) minLevel = std::min(minLevel, hg.label(e_in));
                        }
                    }
                    assert(hg.excess(u) >= 0);
                    if (hg.excess(u) == 0) break;
                }
            }
            return minLevel;
        }

        size_t iterateOverPins(const Node u, CutterState<Type>& cs, PinIterator& pinIt, const PinIterator endIndex) {
            const bool is_edge_in = hg.is_edge_in(u);

            size_t minLevel = n;

            for (; pinIt < endIndex; ++pinIt) {
                const Pin pin = *pinIt;
                const Node v = pin.pin;
                if (!cs.n.isSourceReachable(v) || v == target) {
                    Flow residual = is_edge_in ? std::min(pin.flow, hg.excess(u)) : hg.excess(u);
                    if (residual > 0) {
                        if (hg.label(u) == hg.label(v) + 1) {
                            numPushes++;
                            if (is_edge_in) {
                                hg.push_edgeIn_to_node(v, pin.he_inc_iter, residual);
                            } else {
                                hg.push_edgeOut_to_node(v, pin.he_inc_iter, residual);
                            }
                            if (!(v == target) && inQueue.set(v)) { nodes.push_back(v); }
                        } else {
                            assert(hg.label(u) <= hg.label(v));
                            minLevel = std::min(minLevel, hg.label(v));
                        }
                    }
                }
                if (hg.excess(u) == 0) break;
            }
            return minLevel;
        }

        void dischargeNode(const Node u, CutterState<Type>& cs) {
            assert(u < hg.numNodes());
            assert(hg.excess(u) > 0);

            const InHeIndex beginIt = current_hyperedge[u];
            assert(beginIt >= hg.beginIndexHyperedges(u) && beginIt <= hg.endIndexHyperedges(u));
            const InHeIndex endIt = hg.endIndexHyperedges(u);
            InHeIndex edgeIt = beginIt;

            size_t minLevel = iterateOverEdges(u, cs, edgeIt, endIt);

            if (hg.excess(u) > 0) {
                edgeIt = hg.beginIndexHyperedges(u);
                minLevel = std::min(minLevel, iterateOverEdges(u, cs, edgeIt, beginIt));
                assert(minLevel + 1 > hg.label(u));
                hg.label(u) = minLevel + 1;
                current_hyperedge[u] = hg.beginIndexHyperedges(u);
                workSinceLastRelabel += (endIt - hg.beginIndexHyperedges(u)) * 2 + BETA;
            } else {
                current_hyperedge[u] = edgeIt;
            }
        }

        void dischargeEdgeNodeIn(const Node e_in, CutterState<Type>& cs) {
            assert(hg.excess(e_in) > 0);
            assert(hg.is_edge_in(e_in));

            const Hyperedge e = hg.edgeFromLawlerNode(e_in);
            const PinIterator beginIt = current_pin_e_in[e];
            assert(beginIt >= hg.beginPinsIn(e) && beginIt <= hg.endPinsIn(e));
            const PinIterator endIt = hg.endPinsIn(e);
            PinIterator pinIt = beginIt;

            size_t minLevel = n;

            const Node e_out = hg.edge_node_out(e);
            Flow residual = std::min(hg.capacity(e) - hg.flow(e), hg.excess(e_in));
            if (residual > 0) {
                if (hg.label(e_in) == hg.label(e_out) + 1) {
                    numPushes++;
                    hg.push_edgeIn_to_edgeOut(e_in, e_out, residual);
                    if (inQueue.set(e_out)) { nodes.push_back(e_out); }
                } else {
                    assert(hg.label(e_in) <= hg.label(e_out));
                    minLevel = std::min(minLevel, hg.label(e_out));
                }
            }

            // Scan pins
            if (hg.excess(e_in) > 0) {
                minLevel = std::min(iterateOverPins(e_in, cs, pinIt, endIt), minLevel);
            }

            if (hg.excess(e_in) > 0) {
                // relabel
                pinIt = hg.beginPinsIn(e);
                minLevel = std::min(minLevel, iterateOverPins(e_in, cs, pinIt, beginIt));
                numRelabel++;
                assert(minLevel + 1 > hg.label(e_in));
                hg.label(e_in) = minLevel + 1;
                current_pin_e_in[e] = hg.beginPinsIn(e);
                workSinceLastRelabel += endIt - hg.beginPinsIn(e) + BETA;
            } else {
                current_pin_e_in[e] = pinIt;
            }

        }

        void dischargeEdgeNodeOut(const Node e_out, CutterState<Type>& cs) {
            assert(hg.excess(e_out) > 0);
            assert(hg.is_edge_out(e_out));

            const Hyperedge e = hg.edgeFromLawlerNode(e_out);
            const PinIterator beginIt = current_pin_e_out[e];
            assert(beginIt >= hg.beginPinsIn(e) && beginIt <= hg.endPinsIn(e));
            const PinIterator endIt = hg.endPinsIn(e);
            PinIterator pinIt = beginIt;

            // Scan pins
            size_t minLevel = iterateOverPins(e_out, cs, pinIt, endIt);

            if (hg.excess(e_out) > 0) {
                const Node e_in = hg.edge_node_in(e);
                Flow residual = std::min(hg.flow(e), hg.excess(e_out));
                if (residual > 0) {
                    if (hg.label(e_out) == hg.label(e_in) + 1) {
                        numPushes++;
                        hg.push_edgeOut_to_edgeIn(e_out, e_in, residual);
                        if (inQueue.set(e_in)) {
                            nodes.push_back(e_in);
                        }
                    } else {
                        assert(hg.label(e_out) <= hg.label(e_in));
                        minLevel = std::min(minLevel, hg.label(e_in));
                    }
                }
            }

            if (hg.excess(e_out) > 0) {
                // relabel
                pinIt = hg.beginPinsIn(e);
                minLevel = std::min(minLevel, iterateOverPins(e_out, cs, pinIt, beginIt));
                numRelabel++;
                assert(minLevel + 1 > hg.label(e_out));
                hg.label(e_out) = minLevel + 1;
                workSinceLastRelabel += endIt - hg.beginPinsIn(e) + BETA;
                current_pin_e_out[e] = hg.beginPinsIn(e);
            } else {
                current_pin_e_out[e] = pinIt;
            }

        }

        void pushRelabel(CutterState<Type>& cs, const size_t nm) {
            while (!nodes.empty()) {
                Node u = nodes.front();
                nodes.pop_front();

                if (hg.label(u) >= n) {
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


                if (hg.excess(u) > 0 && hg.label(u) < n) {
                    nodes.push_back(u);
                } else {
                    inQueue.reset(u);
                }

                if (workSinceLastRelabel * globalUpdateFreq > nm) {
                    timer.start("globalUpdate", "phase1");
                    globalUpdate(cs, n);
                    numGlobalUpdate++;
                    timer.stop("globalUpdate");
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

            piercingNode = cs.sourcePiercingNodes.begin()->node;
            target = cs.targetPiercingNodes.begin()->node;

            inQueue.reset();


            for (InHeIndex inc_iter : hg.incidentHyperedgeIndices(piercingNode)) {
                InHe& inc_u = hg.getInHe(inc_iter);
                const Hyperedge e = inc_u.e;
                if (!cs.h.areAllPinsSourceReachable__unsafe__(e)) {
                    const Node e_out = hg.edge_node_out(e);
                    Flow residual = hg.getPinOut(hg.getInHe(inc_iter)).flow;
                    if (residual > 0) {
                        hg.push_node_to_edgeOut(piercingNode, inc_iter, residual);
                        if (inQueue.set(e_out)) { nodes.push_back(e_out); }
                    }

                    const Node e_in = hg.edge_node_in(e);
                    residual = hg.capacity(e);
                    if (residual > 0) {
                        hg.push_node_to_edgeIn(piercingNode, inc_iter, residual);
                        if (inQueue.set(e_in)) { nodes.push_back(e_in); }
                    }
                }
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
            pushRelabel(cs, nm);
            timer.stop("phase1");

            timer.start("phase2", "mainLoop");
            timer.start("buildResidualNetwork", "phase2");
            hg.buildResidualNetwork();
            timer.stop("buildResidualNetwork");

            Flow augmented_flow = flowDecomposition(cs);

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
            assert(!found_target);
            timer.stop("growReachable");
            timer.stop("exhaustFlow");

            std::cout << "numPushes: " << numPushes << ", numRelabel: " << numRelabel << ", numGlobalUpdate: " << numGlobalUpdate << ", numLawlerNodes: " << hg.numLawlerNodes() << std::endl;

            return f;
        }

        Flow flowDecomposition(CutterState<Type>& cs) {
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
            int64_t lowest_bottleneck = std::numeric_limits<int64_t>::max();
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
            if (lowest_bottleneck == std::numeric_limits<int64_t>::max()) lowest_bottleneck = 1;

            for (int64_t stack_pointer = stack.size() - 1; stack_pointer > 0; --stack_pointer) {
                const StackFrame& t = stack.at(stack_pointer);
                //std::cout << t.u << " ";
                if (hg.isNode(t.u)) {
                    if (hg.is_edge_out(stack.at(stack_pointer-1).u)) {
                        hg.getPinOut(hg.getInHe(t.he_it)).flow -= bottleneckCapacity;
                        if (writeFlowToResult) hg.getInHe(t.he_it).flow -= hg.flowSent(bottleneckCapacity);
                    } else {
                        hg.getPinIn(hg.getInHe(t.he_it)).flow += bottleneckCapacity;
                        if (writeFlowToResult) hg.getInHe(t.he_it).flow -= hg.flowSent(bottleneckCapacity);
                    }
                } else if (hg.is_edge_in(t.u)) {
                    if (t.he_it != InHeIndex::Invalid()) {
                        hg.getPinIn(hg.getInHe(t.he_it)).flow -= bottleneckCapacity;
                        if (writeFlowToResult) hg.getInHe(t.he_it).flow += hg.flowSent(bottleneckCapacity);
                    } else {
                        hg.flow(hg.edgeFromLawlerNode(t.u)) += bottleneckCapacity;
                    }
                } else {
                    if (t.he_it != InHeIndex::Invalid()) {
                        hg.getPinOut(hg.getInHe(t.he_it)).flow += bottleneckCapacity;
                        if (writeFlowToResult) hg.getInHe(t.he_it).flow += hg.flowSent(bottleneckCapacity);
                    } else {
                        hg.flow(hg.edgeFromLawlerNode(t.u)) -= bottleneckCapacity;
                    }
                }
                if (stack_pointer >= lowest_bottleneck) inQueue.reset(t.u);
                if (stack.at(stack_pointer - 1).u == first_node) break;
            }
            //std::cout << std::endl;
            if (lowest_bottleneck == std::numeric_limits<int64_t>::max()) lowest_bottleneck = 1;
            stack.popDownTo(lowest_bottleneck-1);
            return bottleneckCapacity;
        }

        inline bool tryToGetNode(const Node u, const size_t currentLabel, const size_t n) {
            return hg.label(u) > currentLabel && __atomic_exchange_n(&hg.label(u), currentLabel, __ATOMIC_ACQ_REL) == n;
        }

        std::array<size_t, 2> globalUpdate(CutterState<Type>& cs, size_t n) {
            size_t scannedNodes = 0;
            size_t scannedEdges = 0;
            queue.clear();

            hg.equalizeLabels(n);

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
                                    current_pin_e_in[e] = hg.beginPinsIn(e);
                                    queue.push(e_in);
                                }
                                if (hg.label(e_out) == n) {
                                    hg.label(e_out) = currentLabel;
                                    current_pin_e_out[e] = hg.beginPinsIn(e);
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
                            current_pin_e_in[e] = hg.beginPinsIn(e);
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
                            current_pin_e_out[e] = hg.beginPinsIn(e);
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

        void resetSourcePiercingNodeDistances(CutterState<Type>& cs, bool reset = true) {
            for (auto& sp: cs.sourcePiercingNodes)
                cs.n.setPiercingNodeDistance(sp.node, reset);
        }

    };
}