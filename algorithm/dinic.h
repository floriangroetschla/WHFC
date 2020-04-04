
#include "cutter_state.h"
#include "../datastructure/queue.h"
#include "../datastructure/stack.h"
#include "../datastructure/distance_reachable_sets.h"
#include "ford_fulkerson.h"
#include "../recursive_bisection/timestamp_set.hpp"
#include "../datastructure/copyable_atomic.h"


namespace whfc {
	class DinicBase {
	public:
		using ScanList = LayeredQueue<Node>;
		
		using ReachableNodes = DistanceReachableNodes;
		using ReachableHyperedges = DistanceReachableHyperedges;
		
		using Pin = FlowHypergraph::Pin;
		using InHe = FlowHypergraph::InHe;
		using PinIndexRange = FlowHypergraph::PinIndexRange;
		using DistanceT = DistanceReachableNodes::DistanceT;
		
		
		FlowHypergraph& hg;
		LayeredQueue<Node> queue;
		struct StackFrame {
			Node u;
			InHeIndex parent_he_it;
		};
		FixedCapacityStack<StackFrame> stack;
		int direction = 0;
		std::vector<PinIndex> current_flow_sending_pin, current_flow_receiving_pin, current_pin;
		std::vector<InHeIndex> current_hyperedge;
		
		Flow upperFlowBound;
		
		DinicBase(FlowHypergraph& hg) : hg(hg), queue(hg.maxNumNodes), stack(hg.maxNumNodes),
										current_flow_sending_pin(hg.maxNumHyperedges, PinIndex::Invalid()),
										current_flow_receiving_pin(hg.maxNumHyperedges, PinIndex::Invalid()),
										current_pin(hg.maxNumHyperedges, PinIndex::Invalid()),
										current_hyperedge(hg.maxNumNodes, InHeIndex::Invalid())
		{
		
		}
		
		void flipViewDirection() {
			std::swap(current_flow_sending_pin, current_flow_receiving_pin);
			direction = 1 - direction;
		}
		
	};
	
	class Dinic : public DinicBase {
	public:
		using Type = Dinic;
		
		static constexpr bool same_traversal_as_grow_assimilated = false;
		static constexpr bool grow_reachable_marks_flow_sending_pins_when_marking_all_pins = true;
		static constexpr bool log = false;
		
        class NodeVectorView {
        public:
            NodeVectorView() = default;

            void addVector(std::vector<Node>* vector) {
                vector_pointers.push_back(vector);
                prefix_sizes.push_back(prefix_sizes.size() == 0 ? vector->size() : prefix_sizes.back() + vector->size());
            }

            size_t size() {
                if (prefix_sizes.size() == 0) return 0;
                return prefix_sizes.back();
            }

            void clear() {
                vector_pointers.clear();
                prefix_sizes.clear();
            }

            Node& get(size_t i) {
                size_t vector_index = 0;
                size_t last_size = 0;
                while (!(i >= last_size && i < prefix_sizes[vector_index])) {
                    last_size = prefix_sizes[vector_index];
                    vector_index++;
                }
                assert(i - last_size < vector_pointers[vector_index]->size());
                return (*vector_pointers[vector_index])[i - last_size];
            }

        private:
            std::vector<std::vector<Node>*> vector_pointers;
            std::vector<size_t> prefix_sizes;
        };

        TimeReporter& timer;
		
		Dinic(FlowHypergraph& hg, TimeReporter& timer) : DinicBase(hg), timer(timer), thisLayer_thread_specific(), nextLayer_thread_specific(), node_visited(hg.maxNumNodes)
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
        tbb::enumerable_thread_specific<std::vector<Node>> thisLayer_thread_specific;
        tbb::enumerable_thread_specific<std::vector<Node>> nextLayer_thread_specific;
		std::vector<Node> currentLayer;
        ldc::AtomicTimestampSet<uint16_t> node_visited;

		
		void resetSourcePiercingNodeDistances(CutterState<Type>& cs, bool reset = true) {
			for (auto& sp: cs.sourcePiercingNodes)
				cs.n.setPiercingNodeDistance(sp.node, reset);
		}


		bool searchFromNode(const Node u, CutterState<Type>& cs, const bool augment_flow, std::vector<Node>& nextLayer) {
            auto& n = cs.n;
            auto& h = cs.h;
            bool found_target;
            for (InHe& inc_u : hg.hyperedgesOf(u)) {
                const Hyperedge e = inc_u.e;
                if (!h.areAllPinsSourceReachable__unsafe__(e)) { // Are there pins that were not already visited
                    const bool scanAllPins = !hg.isSaturated(e) || hg.flowReceived(inc_u) > 0;
                    if (!scanAllPins && h.areFlowSendingPinsSourceReachable__unsafe__(e)) // only sending pins can be pushed back
                        continue;

                    if (scanAllPins) {
                        h.reachAllPins(e);
                        assert(n.distance[u] + 1 == h.outDistance[e]);
                        current_pin[e] = hg.pinsNotSendingFlowIndices(e).begin();
                    }

                    const bool scanFlowSending = !h.areFlowSendingPinsSourceReachable__unsafe__(e);
                    if (scanFlowSending) {
                        h.reachFlowSendingPins(e);
                        assert(n.distance[u] + 1 == h.inDistance[e]);
                        current_flow_sending_pin[e] = hg.pinsSendingFlowIndices(e).begin();
                    }

                    auto visit = [&](const Pin& pv) {
                        const Node v = pv.pin;
                        if (node_visited.set(v)) {
                            assert(augment_flow || !n.isTargetReachable(v));
                            assert(augment_flow || !cs.isIsolated(v) || n.distance[v] == n.s.base);	//checking distance, since the source piercing node is no longer a source at the moment
                            found_target |= n.isTarget(v);
                            if (!n.isTarget(v) && !n.isSourceReachable__unsafe__(v)) {
                                assert(v < hg.numNodes());
                                n.reach(v);
                                assert(n.distance[u] + 1 == n.distance[v]);
                                nextLayer.push_back(v);
                                current_hyperedge[v] = hg.beginIndexHyperedges(v);
                            }
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
            return found_target;
		}

		bool buildLayeredNetwork(CutterState<Type>& cs, const bool augment_flow) {
		    timer.start("buildLayeredNetwork");
		    alignDirection(cs);
		    unused(augment_flow);
		    cs.clearForSearch();
		    auto& n = cs.n;
		    auto& h = cs.h;
		    bool found_target = false;

		    NodeVectorView view;
		    currentLayer.clear();
		    node_visited.reset();

            for (auto& sp : cs.sourcePiercingNodes) {
                n.setPiercingNodeDistance(sp.node, false);
                assert(n.isSourceReachable(sp.node));
                currentLayer.push_back(sp.node);
                current_hyperedge[sp.node] = hg.beginIndexHyperedges(sp.node);
            }
            n.hop(); h.hop();

            view.addVector(&currentLayer);

            while (view.size() > 0) {
                // Only execute in parallel if there are enough nodes left
                if (view.size() > 100) {
                    tbb::this_task_arena::isolate( [&]{
                        tbb::parallel_for(tbb::blocked_range<size_t>(0,view.size()), [&](tbb::blocked_range<size_t> r) {
                            std::vector<Node>& nextLayer = nextLayer_thread_specific.local();
                            for (size_t i = r.begin(); i < r.end(); ++i) {
                                if (searchFromNode(view.get(i), cs, augment_flow, nextLayer)) found_target = true;
                            }
                        });
                    } );
                } else {
                    std::vector<Node>& nextLayer = nextLayer_thread_specific.local();
                    for (size_t i = 0; i < view.size(); ++i) {
                        if (searchFromNode(view.get(i), cs, augment_flow, nextLayer)) found_target = true;
                    }
                }

                n.hop(); h.hop();

                view.clear();
                std::swap(thisLayer_thread_specific, nextLayer_thread_specific);
                for (std::vector<Node>& thisLayer : thisLayer_thread_specific) {
                    view.addVector(&thisLayer);

                }
                for (std::vector<Node>& nextLayer : nextLayer_thread_specific) {
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
	
	class ScalingDinic : public DinicBase {
	public:
		using Type = ScalingDinic;
		
		static constexpr bool same_traversal_as_grow_assimilated = false;
		static constexpr bool grow_reachable_marks_flow_sending_pins_when_marking_all_pins = true;
		static constexpr bool log = false;
		
		FlowCommons::Scaling scaling;
		
		ScalingDinic(FlowHypergraph& hg) : DinicBase(hg)
		{
			reset();
		}
		
		
		void reset() {
			// TODO maybe don't reset to full capacity after piercing. maybe something less
			scaling.initialize(hg.maxHyperedgeCapacity);
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
			
			scaling.reset();
			while (cs.flowValue <= upperFlowBound && scaling.use()) {
				if (buildLayeredNetwork<true>(cs))
					cs.flowValue += augmentFlowInLayeredNetwork(cs);
				else
					scaling.reduceCapacity();
			}
			
			bool hasCut = false;
			while (cs.flowValue <= upperFlowBound) {
				hasCut = !buildLayeredNetwork<true>(cs);
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
			
			while (scaling.use()) {
				if (buildLayeredNetwork<true>(cs)) {
					f += augmentFlowInLayeredNetwork(cs);
					break;
				}
				else {
					scaling.reduceCapacity();
				}
			}
			
			if (f == 0) {
				if (buildLayeredNetwork<true>(cs))
					f += augmentFlowInLayeredNetwork(cs);
				else
					scaling.reset();
			}
			
			resetSourcePiercingNodeDistances(cs);
			return f;
		}
		
		Flow recycleDatastructuresFromGrowReachablePhase(CutterState<Type> &cs) {
			if (!cs.augmentingPathAvailableFromPiercing
				|| std::none_of(cs.sourcePiercingNodes.begin(), cs.sourcePiercingNodes.end(),
								[](const auto& sp) { return sp.isReachableFromOppositeSide; })) {
				return 0;
			}
			cs.flipViewDirection();
			resetSourcePiercingNodeDistances(cs, false);
			scaling.disable();
			Flow f = augmentFlowInLayeredNetwork(cs);
			scaling.enable();
			resetSourcePiercingNodeDistances(cs);
			cs.flipViewDirection();
			return f;
		}
		
		void growReachable(CutterState<Type>& cs) {
			scaling.disable();
			bool found_target = buildLayeredNetwork<false>(cs);
			scaling.enable();
			assert(!found_target); unused(found_target);
			resetSourcePiercingNodeDistances(cs);
		}
	
	private:
		
		void resetSourcePiercingNodeDistances(CutterState<Type>& cs, bool reset = true) {
			for (auto& sp: cs.sourcePiercingNodes)
				cs.n.setPiercingNodeDistance(sp.node, reset);
		}
		
		template<bool augment_flow>
		bool buildLayeredNetwork(CutterState<Type>& cs) {
			alignDirection(cs);
			cs.clearForSearch();
			auto& n = cs.n;
			auto& h = cs.h;
			queue.clear();
			bool found_target = false;
			const Flow scaling_capacity = scaling.getCapacity();
			
			for (auto& sp : cs.sourcePiercingNodes) {
				n.setPiercingNodeDistance(sp.node, false);
				assert(n.isSourceReachable(sp.node));
				queue.push(sp.node);
				current_hyperedge[sp.node] = hg.beginIndexHyperedges(sp.node);
			}
			n.hop(); h.hop(); queue.finishNextLayer();
			
			while (!queue.empty()) {
				while (!queue.currentLayerEmpty()) {
					const Node u = queue.pop();
					for (InHe& inc_u : hg.hyperedgesOf(u)) {
						const Hyperedge e = inc_u.e;
						if (hg.capacity(e) >= scaling_capacity && !h.areAllPinsSourceReachable__unsafe__(e)) {
							
							auto visit = [&](const Pin& pv) {
								const Node v = pv.pin;
								assert(augment_flow || !n.isTargetReachable(v));
								assert(augment_flow || !cs.isIsolated(v) || n.distance[v] == n.s.base);    //checking distance, since the source piercing node is no longer a source at the moment
								found_target |= n.isTarget(v);
								if (!n.isTarget(v) && !n.isSourceReachable__unsafe__(v)) {
									n.reach(v);
									assert(n.distance[u] + 1 == n.distance[v]);
									queue.push(v);
									current_hyperedge[v] = hg.beginIndexHyperedges(v);
								}
							};
							
							Flow residual = hg.residualCapacity(e) + hg.absoluteFlowReceived(inc_u);
							
							if (!h.areFlowSendingPinsSourceReachable__unsafe__(e)) { /* scan flow sending pins */
								h.reachFlowSendingPins(e);
								assert(n.distance[u] + 1 == h.inDistance[e]);
								current_flow_sending_pin[e] = hg.pinsSendingFlowIndices(e).begin();
								for (const Pin& pv : hg.pinsSendingFlowInto(e)) {
									const Flow residual_at_v = residual + hg.absoluteFlowSent(pv);
									if (residual_at_v >= scaling_capacity) {
										visit(pv);
									}
								}
							}
							
							if (residual >= scaling_capacity) { /* scan all pins */
								h.reachAllPins(e);
								assert(n.distance[u] + 1 == h.outDistance[e]);
								current_pin[e] = hg.pinsNotSendingFlowIndices(e).begin();
								for (const Pin& pv : hg.pinsNotSendingFlowInto(e)) {
									visit(pv);
								}
							}
						}
					}
				}
				
				n.hop(); h.hop(); queue.finishNextLayer();
			}
			
			n.lockInSourceDistance(); h.lockInSourceDistance(); h.compareDistances(n);
			return found_target;
		}
		
		Flow augmentFlowInLayeredNetwork(CutterState<Type>& cs) {
			alignDirection(cs);
			auto& n = cs.n;
			auto& h = cs.h;
			Flow f = 0;
			const Flow scaling_capacity = scaling.getCapacity();
			
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
						if (hg.capacity(e) < scaling_capacity)
							continue;
						
						const Flow residual = hg.residualCapacity(e) + hg.absoluteFlowReceived(inc_u);
						
						if (h.inDistance[e] == req_dist) { /* scan flow sending */
							for (const PinIndex firstInvalid = hg.pinsSendingFlowIndices(e).end(); current_flow_sending_pin[e] < firstInvalid; current_flow_sending_pin[e]++) {
								const Pin& pv = hg.getPin(current_flow_sending_pin[e]);
								if (residual + hg.absoluteFlowSent(pv) >= scaling_capacity && (n.isTarget(pv.pin) || n.distance[pv.pin] == req_dist)) {
									v = pv.pin;
									inc_v_it = pv.he_inc_iter;
									break;
								}
							}
						}
						
						if (v == invalidNode && residual >= scaling_capacity && /* scan all pins */ req_dist == h.outDistance[e]) {
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
			assert(f >= scaling_capacity);
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