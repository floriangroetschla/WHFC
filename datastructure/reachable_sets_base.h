#pragma once

#include "../definitions.h"
#include "../push_relabel/flow_hypergraph.h"
#include "../datastructure/copyable_atomic.h"

namespace whfc {
    template<class Hypergraph>
	class ReachableNodesBase {
	public:
		explicit ReachableNodesBase(const Hypergraph& _hg) : hg(_hg) { }
		
		// boilerplatey constructor and assignment operator stuff because of the const ref flow_hg. which we don't want to ever change.
		
		ReachableNodesBase(const ReachableNodesBase& o) : hg(o.hg) {
			copy_const_ref(o);
		}
		
		ReachableNodesBase& operator=(ReachableNodesBase&& o) {
			copy_const_ref(o);	//this cast works surprisingly.
			return *this;
		}
		
		ReachableNodesBase& operator=(const ReachableNodesBase& o) {
			copy_const_ref(o);
			return *this;
		}

		inline void reach(const Node u) {
			sourceReachableWeight += hg.nodeWeight(u);
		}

		inline void unreachSource(const Node u) {
			sourceReachableWeight -= hg.nodeWeight(u);
		}

		inline void unreachTarget(const Node u) {
			targetReachableWeight -= hg.nodeWeight(u);
		}
		
		inline void unsettleSource(const Node u) {
			sourceWeight -= hg.nodeWeight(u);
			numSettledNodes--;
		}
		
		inline void unsettleTarget(const Node u) {
			targetWeight -= hg.nodeWeight(u);
			numSettledNodes--;
		}

		inline void settle(const Node u) {
			sourceWeight += hg.nodeWeight(u);
			numSettledNodes++;
		}
		
		inline void reachTarget(const Node u) {
			targetReachableWeight += hg.nodeWeight(u);
		}
		
		inline void settleTarget(const Node u) {
			targetWeight += hg.nodeWeight(u);
			numSettledNodes++;
		}

		void resetSourceReachableToSource() {
			sourceReachableWeight = sourceWeight;
		}

		void flipViewDirection() {
			std::swap(sourceReachableWeight, targetReachableWeight);
			std::swap(sourceWeight, targetWeight);
		}
		
		void fullReset() {
			sourceReachableWeight = NodeWeight(0);
			sourceWeight = NodeWeight(0);
			targetReachableWeight = NodeWeight(0);
			targetWeight = NodeWeight(0);
		}

		NodeWeight sourceReachableWeight = NodeWeight(0), sourceWeight = NodeWeight(0), targetReachableWeight = NodeWeight(0), targetWeight = NodeWeight(0);
		size_t numSettledNodes = 0;
		
	protected:
		const Hypergraph& hg;
		
		void copy_const_ref(const ReachableNodesBase& o) {
			sourceReachableWeight = o.sourceReachableWeight;
			targetReachableWeight = o.targetReachableWeight;
			sourceWeight = o.sourceWeight;
			targetWeight = o.targetWeight;
		}

	};
}