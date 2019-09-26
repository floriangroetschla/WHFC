#pragma once

#include "cutter_state.h"
#include "../datastructure/stack.h"
#include "../datastructure/queue.h"
#include "../datastructure/bitset_reachable_sets.h"
//#include "../datastructure/timestamp_reachable_sets.h"

namespace whfc {

	//template<class ScanList>
	template<typename ScanList, bool capacityScaling, bool alwaysSetParent = true>
	class FordFulkerson /* : public FlowAlgorithm */ {
	public:
		static constexpr bool debug = false;
		
		using Type = FordFulkerson;
		//using ScanList = FixedCapacityStack<Node>;

		using ReachableNodes = BitsetReachableNodes;
		using ReachableHyperedges = BitsetReachableHyperedges;
		//using ReachableNodes = TimestampReachableNodes;
		//using ReachableHyperedges = TimestampReachableHyperedges;


		using Pin = FlowHypergraph::Pin;
		using InHe = FlowHypergraph::InHe;

		FlowHypergraph& hg;
		ScanList nodes_to_scan;
		
		struct Parent {
			InHeIndex 	parentIncidenceIterator = InHeIndex::Invalid(),
						currentIncidenceIterator = InHeIndex::Invalid();
		};
		std::vector<Parent> parent;

		explicit FordFulkerson(FlowHypergraph& hg) : hg(hg), nodes_to_scan(hg.numNodes()), parent(hg.numNodes()) { }

		static constexpr Flow InitialScalingCapacity = 1 << 24; //NOTE choose sensibly
		static constexpr Flow ScalingCutOff = 4; //NOTE choose sensibly
		Flow scalingCapacity = InitialScalingCapacity;


		Flow recycleDatastructuresFromGrowReachablePhase(CutterState<Type> &cs) {
			Flow flow = 0;
			if constexpr (alwaysSetParent) {
				if (cs.augmentingPathAvailableFromPiercing) {
					for (Node s : cs.sourcePiercingNodes) {
						if (cs.n.isTargetReachable(s)) {
							cs.flipViewDirection();
							flow += augmentFromTarget(cs.n, s);
							cs.flipViewDirection();
							break;	//no VD label propagation --> only one path
						}
					}
				}
			}
			return flow;
		}

		Flow exhaustFlow(CutterState<Type>& cs) {
			Flow flow = 0;
			flow += recycleDatastructuresFromGrowReachablePhase(cs);
			Flow diff = -1;
			if constexpr (capacityScaling) {
				while (scalingCapacity > ScalingCutOff) {
					while (diff != 0) {
						diff = growWithScaling(cs);
						flow += diff;
					}
					scalingCapacity /= 2;
					diff = -1;
				}
			}
			while (diff != 0) {
				diff = growWithoutScaling<true>(cs);
				flow += diff;
			}
			scalingCapacity = InitialScalingCapacity;
			return flow;
		}

		Flow growFlowOrSourceReachable(CutterState<Type>& cs) {
			if constexpr (capacityScaling) {
				while (scalingCapacity > ScalingCutOff) {
					Flow flow = growWithScaling(cs);
					if (flow != 0)
						return flow;
					else
						scalingCapacity /= 2;
				}
			}
			return growWithoutScaling<true>(cs);
		}

		Flow augmentFromTarget(ReachableNodes& n, const Node target) {
			Flow bottleneckCapacity = maxFlow;
			Node v = target;
			while (!n.isSource(v)) {
				Parent p = parent[v];
				const Flow residual = hg.residualCapacity(hg.getInHe(p.parentIncidenceIterator), hg.getInHe(p.currentIncidenceIterator));
				bottleneckCapacity = std::min(bottleneckCapacity, residual);
				v = hg.getPin(hg.getInHe(p.parentIncidenceIterator)).pin;
			}
			AssertMsg(bottleneckCapacity > 0, "Bottleneck capacity not positive");
			v = target;
			LOG << V(v) << V(bottleneckCapacity);
			while (!n.isSource(v)) {
				Parent p = parent[v];
				hg.routeFlow(hg.getInHe(p.parentIncidenceIterator), hg.getInHe(p.currentIncidenceIterator), bottleneckCapacity);
				v = hg.getPin(hg.getInHe(p.parentIncidenceIterator)).pin;
				LOG << V(v);
			}
			return bottleneckCapacity;
		}

		template<bool augment_flow>
		Flow growWithoutScaling(CutterState<Type>& cs) {
			LOG << "grow without scaling";
			cs.clearForSearch();
			ReachableNodes& n = cs.n;
			ReachableHyperedges& h = cs.h;
			nodes_to_scan.clear();
			for (Node s : cs.sourcePiercingNodes)
				nodes_to_scan.push(s);
			
			while (!nodes_to_scan.empty()) {
				const Node u = nodes_to_scan.pop();
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
							AssertMsg(augment_flow || !n.isTargetReachable(v), "Not augmenting flow but target side is reachable.");
							if (!n.isSourceReachable(v)) {		//don't do VD label propagation
								n.reach(v);
								nodes_to_scan.push(v);
								if constexpr (augment_flow || alwaysSetParent)
									parent[v] = { inc_u_iter, pv.he_inc_iter };
								if constexpr (augment_flow)
									if (n.isTarget(v))
										return augmentFromTarget(n, v);
							}
						}
					}
				}
			}
			return 0;
		}


		/*
		 * Note: capacity scaling is implemented separately from search without capacity scaling, as capacity scaling pruning requires more memory accesses than plain search
		 */
		Flow growWithScaling(CutterState<Type>& cs) {
			AssertMsg(scalingCapacity > 1, "Don't call this method with ScalingCapacity <= 1. Use growWithoutScaling instead.");
			LOG << "grow with scaling" << V(scalingCapacity);
			cs.clearForSearch();
			ReachableNodes& n = cs.n;
			ReachableHyperedges& h = cs.h;
			nodes_to_scan.clear();
			for (Node s : cs.sourcePiercingNodes)
				nodes_to_scan.push(s);

			while (!nodes_to_scan.empty()) {
				const Node u = nodes_to_scan.pop();
				for (InHeIndex inc_u_iter : hg.incidentHyperedgeIndices(u)) {
					const InHe& inc_u = hg.getInHe(inc_u_iter);
					const Hyperedge e = inc_u.e;
					//can push at most flow(e) back into flow-sending pin and at most residual(e) = capacity(e) - flow(e) further flow.
					//other pins can receive at most residual(e) <= capacity(e). so checking capacity(e) < scalingCapacity is a good pruning rule
					if (hg.capacity(e) < scalingCapacity)
						continue;

					Flow residualCapacity = hg.absoluteFlowReceived(inc_u) + hg.residualCapacity(e);
					if (!h.areFlowSendingPinsSourceReachable(e)) {
						h.reachFlowSendingPins(e);
						for (const Pin& pv : hg.pinsSendingFlowInto(e)) {
							if (residualCapacity + hg.absoluteFlowSent(pv) >= scalingCapacity) {//residual = flow received by u + residual(e) + flow sent by v
								const Node v = pv.pin;
								if (!n.isSourceReachable(v)) {
									n.reach(v);
									nodes_to_scan.push(v);
									parent[v] = parent[v] = { inc_u_iter, pv.he_inc_iter };
									if (n.isTarget(v))
										return augmentFromTarget(n, v);
								}
							}
						}
					}

					if (residualCapacity >= scalingCapacity && !h.areAllPinsSourceReachable(e) && (!hg.isSaturated(e) || hg.flowReceived(inc_u) > 0)) {
						h.reachAllPins(e);
						for (const Pin& pv : hg.pinsNotSendingFlowInto(e)) {
							const Node v = pv.pin;
							if (!n.isSourceReachable(v)) {
								n.reach(v);
								nodes_to_scan.push(v);
								parent[v] = parent[v] = { inc_u_iter, pv.he_inc_iter };
								if (n.isTarget(v))
									return augmentFromTarget(n, v);
							}
						}
					}
				}
			}
			return 0;
		}

		void growReachable(CutterState<Type>& cs) {
			growWithoutScaling<false>(cs);
		}
	};

	using ScalingFordFulkerson = FordFulkerson<FixedCapacityStack<Node>, true>;
	using BasicFordFulkerson = FordFulkerson<FixedCapacityStack<Node>, false>;

	using ScalingEdmondsKarp = FordFulkerson<LayeredQueue<Node>, true>;
	using BasicEdmondsKarp = FordFulkerson<LayeredQueue<Node>, false>;

	class DepthFirstFordFulkerson {
	public:
		using Type = DepthFirstFordFulkerson;
		using ReachableNodes = BitsetReachableNodes;
		using ReachableHyperedges = BitsetReachableHyperedges;
		//using ReachableNodes = TimestampReachableNodes;
		//using ReachableHyperedges = TimestampReachableHyperedges;

		using Pin = FlowHypergraph::Pin;
		using InHe = FlowHypergraph::InHe;

		FlowHypergraph& hg;

		Flow growWithoutScaling(CutterState<Type>& cs) {
			/*
			for (Node s : cs.sourcePiercingNodes) {
				const InHeIndex e_it = hg.beginIndexHyperedges(s);
				const PinIndex pin_it = hg.getInHe(e_it).pin_iter;
				stack.push( { e_it, pin_it } );
			}
			 */
			return 0;
		}

		struct StackElement {
			InHeIndex e_it;
			PinIndex pin_it;
		};

		FixedCapacityStack<StackElement> stack;


	};
}