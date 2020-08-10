#pragma once

#include "../hyperflowcutter/cutter_state.h"
#include "../datastructure/queue.h"
#include "../datastructure/stack.h"
#include "../datastructure/distance_reachable_sets.h"
#include "dinic.h"
#include "ford_fulkerson.h"


namespace whfc {
	class DinicThreadLocalVectors : public DinicBase {
	public:
		using Type = DinicThreadLocalVectors;
		
		static constexpr bool same_traversal_as_grow_assimilated = false;
		static constexpr bool grow_reachable_marks_flow_sending_pins_when_marking_all_pins = true;
		static constexpr bool log = false;


        TimeReporter& timer;
		
		DinicThreadLocalVectors(FlowHypergraph& hg, TimeReporter& timer, size_t numThreads) : DinicBase(hg), timer(timer),
		    thisLayer_thread_specific(new tls_enumerable_thread_specific<std::vector<Node>>()),
		    nextLayer_thread_specific(new tls_enumerable_thread_specific<std::vector<Node>>())
		{
			reset();
		}
		
		void reset() {
		
		}
		
		void alignDirection(CutterState<Type>& cs) {
			if (direction != cs.currentViewDirection()) {
				flipViewDirection();
			}
		}
		
		ScanList& getScanList() {
			return queue;
		}
		
		bool exhaustFlow(CutterState<Type>& cs) {
			cs.flowValue += recycleDatastructuresFromGrowReachablePhase(cs);
			bool hasCut = false;
			while (cs.flowValue <= upperFlowBound) {
				hasCut = !buildLayeredNetwork(cs, true);
				if (hasCut || cs.flowValue >= upperFlowBound) {
					break;
				}
				else {
					cs.flowValue += augmentFlowInLayeredNetwork(cs);
				}
			}
			resetSourcePiercingNodeDistances(cs);
			return hasCut;
		}
		
		Flow growFlowOrSourceReachable(CutterState<Type>& cs) {
			Flow f = 0;
			if (buildLayeredNetwork(cs, true))
				f += augmentFlowInLayeredNetwork(cs);
			resetSourcePiercingNodeDistances(cs);
			return f;
		}
		
		Flow recycleDatastructuresFromGrowReachablePhase(CutterState<Type> &cs) {
			if (!cs.augmentingPathAvailableFromPiercing || std::none_of(cs.sourcePiercingNodes.begin(), cs.sourcePiercingNodes.end(),
																		[](const auto& sp) { return sp.isReachableFromOppositeSide; })) {
				return 0;
			}
			cs.flipViewDirection();
			resetSourcePiercingNodeDistances(cs, false);
			Flow f = augmentFlowInLayeredNetwork(cs);
			resetSourcePiercingNodeDistances(cs);
			cs.flipViewDirection();
			return f;
		}
		
		void growReachable(CutterState<Type>& cs) {
			bool found_target = buildLayeredNetwork(cs, false);
			assert(!found_target); unused(found_target);
			resetSourcePiercingNodeDistances(cs);
		}
		
	private:
        template<typename T>
        using tls_enumerable_thread_specific = tbb::enumerable_thread_specific<T, tbb::cache_aligned_allocator<T>, tbb::ets_key_per_instance>;

        std::unique_ptr<tls_enumerable_thread_specific<std::vector<Node>>> thisLayer_thread_specific;
        std::unique_ptr<tls_enumerable_thread_specific<std::vector<Node>>> nextLayer_thread_specific;
        tls_enumerable_thread_specific<whfc::NodeWeight> sourceReachableWeight_thread_specific;
        bool found_target;

        static constexpr size_t max_write_buffer_size = 128;
        size_t write_buffer_size = max_write_buffer_size;

        static constexpr bool sort_layers = false;

        void resetSourcePiercingNodeDistances(CutterState<Type>& cs, bool reset = true) {
			for (auto& sp: cs.sourcePiercingNodes)
				cs.n.setPiercingNodeDistance(sp.node, reset);
		}

		void searchFromNode(const Node u, CutterState<Type>& cs) {
            auto& n = cs.n;
            auto& h = cs.h;

            auto processIncidences = [&](const Node u, const tbb::blocked_range<std::vector<InHe>::iterator>& hes) {
                std::vector<Node>& nextLayer = nextLayer_thread_specific->local();
                whfc::NodeWeight& weight_of_visited_nodes = sourceReachableWeight_thread_specific.local();

                for (InHe& in_he : hes) {
                    const Hyperedge e = in_he.e;

                    if (!h.areAllPinsSourceReachable__unsafe__(e)) {
                        bool scanAllPins = !hg.isSaturated(e) || hg.flowReceived(in_he) > 0;
                        if (!scanAllPins && h.areFlowSendingPinsSourceReachable__unsafe__(e))
                            continue;

                        scanAllPins = scanAllPins && h.outDistance[e].exchange(h.runningDistance, std::memory_order_relaxed) < h.s.base;
                        if (scanAllPins) {
                            assert(n.distance[u] + 1 == h.outDistance[e]);
                            current_pin[e] = hg.pinsNotSendingFlowIndices(e).begin();
                        }

                        const bool scanFlowSending = !h.areFlowSendingPinsSourceReachable__unsafe__(e)
                                && h.inDistance[e].exchange(h.runningDistance,std::memory_order_relaxed) < h.s.base;
                        if (scanFlowSending) {
                            assert(n.distance[u] + 1 == h.inDistance[e]);
                            current_flow_sending_pin[e] = hg.pinsSendingFlowIndices(e).begin();
                        }

                        auto visit = [&](const Pin &pv) {
                            const Node& v = pv.pin;
                            if (n.isTarget(v)) found_target = true;

                            if (!n.isTarget(v) && !n.isSourceReachable__unsafe__(v) && n.distance[v].exchange(n.runningDistance, std::memory_order_relaxed) < n.s.base) {
                                assert(v < hg.numNodes());
                                weight_of_visited_nodes += hg.nodeWeight(v);
                                assert(n.distance[u] + 1 == n.distance[v]);
                                nextLayer.push_back(v);
                                current_hyperedge[v] = hg.beginIndexHyperedges(v);
                            }
                        };

                        if (scanFlowSending)
                            for (const Pin &pv : hg.pinsSendingFlowInto(e))
                                visit(pv);

                        if (scanAllPins)
                            for (const Pin &pv : hg.pinsNotSendingFlowInto(e))
                                visit(pv);
                    }
                }
            };

            tbb::blocked_range<std::vector<InHe>::iterator> range(hg.beginHyperedges(u), hg.endHyperedges(u));
            if (hg.hyperedgesOf(u).size() > 5000) {
                tbb::parallel_for(range, [=](const tbb::blocked_range<std::vector<InHe>::iterator>& hes) {
                    processIncidences(u, hes);
                });
            } else {
                processIncidences(u, range);
            }
		}

		bool buildLayeredNetwork(CutterState<Type>& cs, const bool augment_flow) {
		    timer.start("buildLayeredNetwork");
		    alignDirection(cs);
		    unused(augment_flow);
		    cs.clearForSearch();
		    auto& n = cs.n;
		    auto& h = cs.h;
		    found_target = false;

            std::vector<Node>& thisLayer = thisLayer_thread_specific->local();
            thisLayer.clear();

            bool nodes_left = false;

            for (auto& sp : cs.sourcePiercingNodes) {
                n.setPiercingNodeDistance(sp.node, false);
                assert(n.isSourceReachable(sp.node));
                thisLayer.push_back(sp.node);
                current_hyperedge[sp.node] = hg.beginIndexHyperedges(sp.node);
                nodes_left = true;
            }
            n.hop(); h.hop();

            while (nodes_left) {
                timer.start("searchFromNodes", "buildLayeredNetwork");
                tbb::parallel_for_each(*thisLayer_thread_specific, [&](const std::vector<Node>& vector) {
                    tbb::parallel_for(tbb::blocked_range<size_t>(0, vector.size()), [&](const tbb::blocked_range<size_t>& nodes) {
                        for (size_t i = nodes.begin(); i < nodes.end(); ++i) {
                            searchFromNode(vector[i], cs);
                        }
                    });
                });
                timer.stop("searchFromNodes");

                n.hop(); h.hop();

                nodes_left = false;
                for (whfc::NodeWeight& sourceReachableWeight : sourceReachableWeight_thread_specific) {
                    if (sourceReachableWeight > 0) nodes_left = true;
                    n.sourceReachableWeight += sourceReachableWeight;
                    sourceReachableWeight = 0;
                }

                std::swap(thisLayer_thread_specific, nextLayer_thread_specific);

                for (std::vector<Node>& nextLayer : *nextLayer_thread_specific) {
                    nextLayer.clear();
                }


            }
            n.lockInSourceDistance(); h.lockInSourceDistance();
            h.compareDistances(n);

            LOGGER_WN << V(found_target) << "#BFS layers =" << (n.s.upper_bound - n.s.base);
            timer.stop("buildLayeredNetwork");
            return found_target;
		}

		
		Flow augmentFlowInLayeredNetwork(CutterState<Type>& cs) {
			alignDirection(cs);
			auto& n = cs.n;
			auto& h = cs.h;
			Flow f = 0;
			
			for (auto& sp : cs.sourcePiercingNodes) {
				assert(stack.empty());
				stack.push({ sp.node, InHeIndex::Invalid() });
				
				while (!stack.empty()) {
					const Node u = stack.top().u;
					Node v = invalidNode;
					InHeIndex inc_v_it = InHeIndex::Invalid();
					DistanceT req_dist = n.distance[u] + 1;
					assert(!n.isDistanceStale(u));
					assert(stack.size() + n.sourceBaseDistance() == req_dist);
					InHeIndex& he_it = current_hyperedge[u];
					for ( ; he_it < hg.endIndexHyperedges(u); he_it++) {
						InHe& inc_u = hg.getInHe(he_it);
						const Hyperedge e = inc_u.e;
						const Flow residual = hg.residualCapacity(e) + hg.absoluteFlowReceived(inc_u);
						assert((residual > 0) == (!hg.isSaturated(e) || hg.absoluteFlowReceived(inc_u) > 0));
						const bool scanAll = req_dist == h.outDistance[e] && residual > 0;
						const bool scanFlowSending = req_dist == h.inDistance[e];
						
						if (scanFlowSending) {
							for (const PinIndex firstInvalid = hg.pinsSendingFlowIndices(e).end(); current_flow_sending_pin[e] < firstInvalid; current_flow_sending_pin[e]++) {
								const Pin& pv = hg.getPin(current_flow_sending_pin[e]);
								if (residual + hg.absoluteFlowSent(pv) > 0 && (n.isTarget(pv.pin) || n.distance[pv.pin] == req_dist)) {
									v = pv.pin;
									inc_v_it = pv.he_inc_iter;
									break;
								}
							}
						}
						
						if (scanAll && v == invalidNode) {
							for (const PinIndex firstInvalid = hg.pinsNotSendingFlowIndices(e).end(); current_pin[e] < firstInvalid; current_pin[e]++) {
								const Pin& pv = hg.getPin(current_pin[e]);
								if (n.isTarget(pv.pin) || n.distance[pv.pin] == req_dist) {
									v = pv.pin;
									inc_v_it = pv.he_inc_iter;
									break;
								}
							}
						}
						
						if (v != invalidNode)
							break;		//don't advance hyperedge iterator
					}
					
					if (v == invalidNode) {
						assert(current_hyperedge[u] == hg.endIndexHyperedges(u));
						stack.pop();
						// Note: the iteration of u's predecessor on the stack still points to u. setting the distance to unreachable prevents the search from pushing u again.
						// It is fine to destroy the reachability datastructures, since we know that this function increases the flow.
						// An alternative method would be to advance the iteration manually, which would be hacky.
						n.distance[u] = ReachableNodes::unreachableDistance;
					}
					else {
						if (n.isTarget(v))
							f += augmentFromTarget(inc_v_it);
						else
							stack.push( { v, inc_v_it } );
					}
					
				}
			}
			assert(f > 0);
			return f;
		}
		
		
		
		Flow augmentFromTarget(InHeIndex inc_target_it) {
			Flow bottleneckCapacity = maxFlow;
			int64_t lowest_bottleneck = std::numeric_limits<int64_t>::max();
			InHeIndex inc_v_it = inc_target_it;
			for (int64_t stack_pointer = stack.size() - 1; stack_pointer >= 0; --stack_pointer) {
				const StackFrame& t = stack.at(stack_pointer);
				const Flow residual = hg.residualCapacity(hg.getInHe(current_hyperedge[t.u]), hg.getInHe(inc_v_it));
				if (residual <= bottleneckCapacity) {
					bottleneckCapacity = residual;
					lowest_bottleneck = stack_pointer;
				}
				inc_v_it = t.parent_he_it;
			}
			assert(bottleneckCapacity > 0);
			inc_v_it = inc_target_it;
			for (int64_t stack_pointer = stack.size() - 1; stack_pointer >= 0; --stack_pointer) {
				const StackFrame& t = stack.at(stack_pointer);
				hg.routeFlow(hg.getInHe(current_hyperedge[t.u]), hg.getInHe(inc_v_it), bottleneckCapacity);
				inc_v_it = t.parent_he_it;
			}
			stack.popDownTo(lowest_bottleneck);
			return bottleneckCapacity;
		}
		
	};
}