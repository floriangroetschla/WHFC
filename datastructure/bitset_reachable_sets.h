#pragma once

#include "reachable_sets_base.h"
#include "bitvector.h"
namespace whfc {



		class BitsetReachableNodes : public ReachableNodesBase {
		public:
			using Base = ReachableNodesBase;
			using Type = BitsetReachableNodes;

			inline size_t capacity() const { return S.size(); }
			inline bool isSource(const Node u) const { return S[u]; }
			inline bool isSourceReachable(const Node u) const { return SR[u]; }
			inline bool isTarget(const Node u) const { return T[u]; }
			inline bool isTargetReachable(const Node u) const { return TR[u]; }
			inline void reach(const Node u, const NodeWeight w) { SR.set(u); Base::reach(u, w); }
			inline void settle(const Node u, const NodeWeight w) { S.set(u); Base::settle(u, w); }


			void flipViewDirection() {
				std::swap(S, T);
				std::swap(SR, TR);
				Base::flipViewDirection();
			}

			void resetSourceReachableToSource() {
				SR = S;
				Base::resetSourceReachableToSource();
			}

			void verifyDisjoint() const {
				assert((SR & TR).none());
				assert((S & T).none());
			}

			void verifySettledIsSubsetOfReachable() const {
				assert(S.is_subset_of(SR));
				assert(T.is_subset_of(TR));
			}

		protected:
			BitVector S, SR, T, TR;
		};

		class BitsetReachableHyperedges {
		public:
			using Type = BitsetReachableHyperedges;

			inline size_t capacity() const { return IN_SETTLED_S.capacity(); }
			inline bool areAllPinsSources(const Hyperedge e) const { return OUT_SETTLED_S[e]; }
			inline bool areAllPinsSourceReachable(const Hyperedge e) const { return OUT_REACHED_S[e]; }
			inline void settleAllPins(const Hyperedge e) { assert(!areAllPinsSources(e)); OUT_SETTLED_S.set(e); IN_SETTLED_S.set(e);}
			inline void reachAllPins(const Hyperedge e) { assert(!areAllPinsSourceReachable(e)); OUT_REACHED_S.set(e); IN_REACHED_S.set(e);}

			inline bool areFlowSendingPinsSources(const Hyperedge e) const { return IN_SETTLED_S[e]; }
			inline bool areFlowSendingPinsSourceReachable(const Hyperedge e) const { return IN_REACHED_S[e]; }
			inline void settleFlowSendingPins(const Hyperedge e) { assert(!areFlowSendingPinsSources(e)); IN_SETTLED_S.set(e); }
			inline void reachFlowSendingPins(const Hyperedge e) { assert(!areFlowSendingPinsSourceReachable(e)); IN_REACHED_S.set(e); }

			void resetSourceReachableToSource() {
				IN_REACHED_S = IN_SETTLED_S;
				OUT_REACHED_S = OUT_SETTLED_S;
			}

			void flipViewDirection() {
				std::swap(IN_SETTLED_S, OUT_SETTLED_T);
				std::swap(OUT_SETTLED_S, IN_SETTLED_T);
				std::swap(IN_REACHED_S, OUT_REACHED_T);
				std::swap(OUT_REACHED_S, IN_REACHED_T);
			}

			void verifyDisjoint() const {
				assert((OUT_REACHED_S & OUT_REACHED_T).none());
				assert((OUT_SETTLED_S & OUT_SETTLED_T).none());
				assert((IN_REACHED_S & IN_REACHED_T).none());
				assert((IN_SETTLED_S & IN_SETTLED_T).none());
			}

			void verifySettledIsSubsetOfReachable() const {
				assert(OUT_SETTLED_S.is_subset_of(OUT_REACHED_S));
				assert(IN_SETTLED_S.is_subset_of(IN_REACHED_S));
				assert(OUT_SETTLED_T.is_subset_of(OUT_REACHED_T));
				assert(IN_SETTLED_T.is_subset_of(IN_REACHED_T));
			}

			protected:
			BitVector IN_SETTLED_S, OUT_SETTLED_S, IN_REACHED_S, OUT_REACHED_S;
			BitVector IN_SETTLED_T, OUT_SETTLED_T, IN_REACHED_T, OUT_REACHED_T;

		};

}