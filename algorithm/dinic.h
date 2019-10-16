
#include "cutter_state.h"
#include "../datastructure/queue.h"
#include "../datastructure/stack.h"
#include "../datastructure/distance_reachable_sets.h"
#include "ford_fulkerson.h"

namespace whfc {
	
	//TODO scaling
	class Dinic {
	public:
		using Type = Dinic;
		using ScanList = LayeredQueue<Node>;

		using ReachableNodes = DistanceReachableNodes;
		using ReachableHyperedges = DistanceReachableHyperedges;
		
		using Pin = FlowHypergraph::Pin;
		using InHe = FlowHypergraph::InHe;
		using PinIndexRange = FlowHypergraph::PinIndexRange;
		using DistanceT = DistanceReachableNodes::DistanceT;
		
		static constexpr bool same_traversal_as_grow_assimilated = false;
		static constexpr bool grow_reachable_marks_flow_sending_pins_when_marking_all_pins = true;
		static constexpr bool log = false;
		
		FlowHypergraph& hg;
		LayeredQueue<Node> queue;
		struct StackFrame {
			Node u;
			InHeIndex parent_he_it;
		};
		FixedCapacityStack<StackFrame> stack;
		std::vector<PinIndex> current_flow_sending_pin;
		std::vector<PinIndex> current_pin;
		std::vector<InHeIndex> current_hyperedge;
		
		Dinic(FlowHypergraph& hg) : hg(hg), queue(hg.numNodes()), stack(hg.numNodes()),
									current_flow_sending_pin(hg.numHyperedges(), PinIndex::Invalid()),
									current_pin(hg.numHyperedges(), PinIndex::Invalid()),
									current_hyperedge(hg.numNodes(), InHeIndex::Invalid())
		{
			
		}
		
		ScanList& getScanList() {
			return queue;
		}
		
		Flow exhaustFlow(CutterState<Type>& cs) {
			Flow f = 0;
			f += recycleDatastructuresFromGrowReachablePhase(cs);
			while (buildLayeredNetwork<true>(cs)) {
				f += augmentFlowInLayeredNetwork(cs);
			}
			resetSourcePiercingNodeDistances(cs);
			return f;
		}
		
		Flow growFlowOrSourceReachable(CutterState<Type>& cs) {
			Flow f = 0;
			if (buildLayeredNetwork<true>(cs))
				f += augmentFlowInLayeredNetwork(cs);
			resetSourcePiercingNodeDistances(cs);
			return f;
		}
		
		
		Flow recycleDatastructuresFromGrowReachablePhase(CutterState<Type> &cs) {
			if (!cs.augmentingPathAvailableFromPiercing ||
				std::none_of(cs.sourcePiercingNodes.begin(), cs.sourcePiercingNodes.end(), [](const auto& sp) { return sp.isReachableFromOppositeSide; }))
			{
				return 0;
			}
			cs.flipViewDirection();
			resetSourcePiercingNodeDistances(cs, false);
			Flow f = augmentFlowInLayeredNetwork(cs);
			//Flow f = 0;
			resetSourcePiercingNodeDistances(cs);
			cs.flipViewDirection();
			return f;
		}
		
		void growReachable(CutterState<Type>& cs) {
			bool found_target = buildLayeredNetwork<false>(cs);
			Assert(!found_target);
			resetSourcePiercingNodeDistances(cs);
		}
		
	private:
		
		void resetSourcePiercingNodeDistances(CutterState<Type>& cs, bool reset = true) {
			for (auto& sp: cs.sourcePiercingNodes)
				cs.n.setPiercingNodeDistance(sp.node, reset);
		}
		
		template<bool augment_flow>
		bool buildLayeredNetwork(CutterState<Type>& cs) {
			cs.clearForSearch();
			auto& n = cs.n;
			auto& h = cs.h;
			queue.clear();
			bool found_target = false;
			std::set<DistanceT> target_distances;	// for debugging purposes. TODO remove
			
			for (auto& sp : cs.sourcePiercingNodes) {
				n.setPiercingNodeDistance(sp.node, false);
				Assert(n.isSourceReachable(sp.node));
				queue.push(sp.node);
				current_hyperedge[sp.node] = hg.beginIndexHyperedges(sp.node);
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
								Assert(n.distance[u] + 1 == h.outDistance[e]);
								current_pin[e] = hg.pinsNotSendingFlowIndices(e).begin();
							}
							
							const bool scanFlowSending = !h.areFlowSendingPinsSourceReachable__unsafe__(e);
							if (scanFlowSending) {
								h.reachFlowSendingPins(e);
								Assert(n.distance[u] + 1 == h.inDistance[e]);
								current_flow_sending_pin[e] = hg.pinsSendingFlowIndices(e).begin();
							}
							
							auto visit = [&](const Pin& pv) {	// TODO scaling. and don't use hg.absoluteFlowSent(pv) if not sending flow into e. (save the lookup)
								const Node v = pv.pin;
								AssertMsg(augment_flow || !n.isTargetReachable(v), "Not augmenting flow but target side is reachable.");
								Assert(augment_flow || !cs.isIsolated(v) || n.distance[v] == n.s.base);	//checking distance, since the source piercing node is no longer a source at the moment
								found_target |= n.isTarget(v);
								if (n.isTarget(v)) {
									target_distances.insert(n.distance[u] + 1);
								}
								if (!n.isTarget(v) && !n.isSourceReachable__unsafe__(v)) {
									n.reach(v);
									Assert(n.distance[u] + 1 == n.distance[v]);
									queue.push(v);
									current_hyperedge[v] = hg.beginIndexHyperedges(v);
								}
							};
							
							
							//for (const Pin& pv : scanAllPins ? hg.pinsOf(e) : hg.pinsSendingFlowInto(e))	//if you do the variant below: disable same_traversal_as_grow_assimilated !!!
							//	visit(pv);
							
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
			
			LOGGER_WN << V(found_target) << "#BFS layers =" << (n.s.upper_bound - n.s.base) << ". Target distances = ";
			for (const DistanceT& d : target_distances)
				LOGGER_WN << (d - n.s.base);
			LOGGER << " ";
			return found_target;
		}
		
		Flow augmentFlowInLayeredNetwork(CutterState<Type>& cs) {
			auto& n = cs.n;
			auto& h = cs.h;
			Flow f = 0;
			
			for (auto& sp : cs.sourcePiercingNodes) {
				Assert(stack.empty());
				stack.push({ sp.node, InHeIndex::Invalid() });
				
				while (!stack.empty()) {
					const Node u = stack.top().u;
					Node v = invalidNode;
					InHeIndex inc_v_it = InHeIndex::Invalid();
					DistanceT req_dist = n.distance[u] + 1;
					Assert(!n.isDistanceStale(u));
					Assert(stack.size() + n.sourceBaseDistance() == req_dist);
					InHeIndex& he_it = current_hyperedge[u];
					for ( ; he_it < hg.endIndexHyperedges(u); he_it++) {
						InHe& inc_u = hg.getInHe(he_it);
						const Hyperedge e = inc_u.e;
						const Flow residual = hg.residualCapacity(e) + hg.absoluteFlowReceived(inc_u);
						Assert((residual > 0) == (!hg.isSaturated(e) || hg.absoluteFlowReceived(inc_u) > 0));
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
						Assert(current_hyperedge[u] == hg.endIndexHyperedges(u));
						stack.pop();
						// Note: the iteration of u's predecessor on the stack still points to u. setting the distance to unreachable prevents the search from pushing u again.
						// It is fine to destroy the reachability datastructures, since we know that this function increases the flow.
						// An alternative method would be to advance the iteration manually, which would be hacky.
						n.distance[u] = ReachableNodes::unreachableDistance;
					}
					else {
						if (n.isTarget(v))
							f += augmentFromTarget(cs, inc_v_it);
						else
							stack.push( { v, inc_v_it } );
					}
					
				}
			}
			Assert(f > 0);
			return f;
		}
		
		Flow augmentFromTarget(CutterState<Type>& cs, InHeIndex inc_target_it) {
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
			AssertMsg(bottleneckCapacity > 0, "Bottleneck capacity not positive");
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