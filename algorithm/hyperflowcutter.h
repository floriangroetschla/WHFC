#pragma once

#include "../datastructure/flow_hypergraph.h"
#include "cutter_state.h"
#include "grow_assimilated.h"
#include "piercing.h"

namespace whfc {

	template<class FlowAlgorithm>
	class HyperFlowCutter {
	public:
		static constexpr const char *algo_name = "HyperFlowCutter";		//this looks weird
		TimeReporter timer;
		FlowHypergraph& hg;
		CutterState<FlowAlgorithm> cs;
		FlowAlgorithm flow_algo;	// = SearchAlgorithm
		Flow upperFlowBound;
		Piercer piercer;
		bool find_most_balanced = true;

		static constexpr bool log = false;
		HyperFlowCutter(FlowHypergraph& hg, NodeWeight maxBlockWeight) : timer(algo_name), hg(hg), cs(hg, maxBlockWeight, timer), flow_algo(hg), upperFlowBound(maxFlow), piercer(hg)
		{
		
		}

		void reset() {
			cs.reset();
			flow_algo.reset();
			upperFlowBound = maxFlow;
			piercer.clear();
			timer.clear();
		}
		
		void flipViewDirection() {
			cs.flipViewDirection();
			piercer.flipViewDirection();
		}

		Node selectPiercingNode() {
			if (cs.notSettledNodeWeight() == 0)
				return invalidNode;
			
			Assert(cs.n.sourceWeight == cs.n.sourceReachableWeight);
			Assert(cs.n.sourceReachableWeight <= cs.n.targetReachableWeight);
			cs.cleanUpCut();
			cs.cleanUpBorder();
			cs.verifyCutPostConditions();
			
			Node piercingNode = piercer.findPiercingNode(cs, cs.borderNodes.sourceSideBorder, cs.maxBlockWeight);
			if (piercingNode == invalidNode) {
				LOGGER << "piercing fallback";
				auto allNodes = hg.nodeIDs();
				piercingNode = piercer.findPiercingNode(cs, allNodes, cs.maxBlockWeight);
			}
			return piercingNode;
		}

		void setPiercingNode(const Node piercingNode) {
			cs.augmentingPathAvailableFromPiercing = cs.n.isTargetReachable(piercingNode);
			cs.sourcePiercingNodes.clear();
			cs.sourcePiercingNodes.emplace_back(piercingNode, cs.n.isTargetReachable(piercingNode));
			cs.settleNode(piercingNode);
			cs.hasCut = false;
		}
		
		bool pierce(bool reject_piercing_if_it_creates_an_augmenting_path = false) {
			timer.start("Piercing");
			Node piercingNode = selectPiercingNode();
			timer.stop("Piercing");
			if (piercingNode == invalidNode)
				return false;
			if (reject_piercing_if_it_creates_an_augmenting_path && cs.n.isTargetReachable(piercingNode))
				return false;
			setPiercingNode(piercingNode);
			LOGGER << "Piercing" << V(piercingNode) << V(cs.augmentingPathAvailableFromPiercing);
			return true;
		}
		
		void exhaustFlowAndGrow() {
			timer.start("Flow");
			if (cs.augmentingPathAvailableFromPiercing) {
				timer.start("Augment", "Flow");
				cs.flowValue += flow_algo.exhaustFlow(cs);
				timer.stop("Augment");
				flipViewDirection();
				timer.start("Grow Backward Reachable", "Flow");
				flow_algo.growReachable(cs);
				timer.stop("Grow Backward Reachable");
			}
			else {
				timer.start("Grow Reachable due to AAP", "Flow");
				flow_algo.growReachable(cs);
				timer.stop("Grow Reachable due to AAP");
			}
			timer.stop("Flow");
			cs.verifyFlowConstraints();
			cs.verifySetInvariants();
			
			cs.hasCut = true;
			if (cs.n.targetReachableWeight <= cs.n.sourceReachableWeight) {
				flipViewDirection();
			}
			timer.start("Grow Assimilated");
			LOGGER << "grow assimilated";
			GrowAssimilated<FlowAlgorithm>::grow(cs, flow_algo.getScanList());
			timer.stop("Grow Assimilated");
			cs.verifyFlowConstraints();
			cs.verifySetInvariants();
			LOGGER << cs.toString();
		}

		//for flow-based interleaving
		bool advanceOneFlowIteration(bool reject_piercing_if_it_creates_an_augmenting_path = false) {
			const bool pierceInThisIteration = cs.hasCut;
			if (pierceInThisIteration)
				if (!pierce(reject_piercing_if_it_creates_an_augmenting_path))
					return false;
			
			timer.start("Flow");
			if (cs.augmentingPathAvailableFromPiercing) {
				timer.start("Augment", "Flow");
				if (pierceInThisIteration)
					cs.flowValue += flow_algo.recycleDatastructuresFromGrowReachablePhase(cs);	//the flow due to recycled datastructures does not matter when deciding whether we have a cut available
				Flow flow_diff = flow_algo.growFlowOrSourceReachable(cs);
				timer.stop("Augment");
				cs.flowValue += flow_diff;
				cs.hasCut = flow_diff == 0;
				if (cs.hasCut) {
					flipViewDirection();
					timer.start("Grow Backward Reachable", "Flow");
					flow_algo.growReachable(cs);
					timer.stop("Grow Backward Reachable");
				}
			}
			else {
				timer.start("Grow Reachable due to AAP", "Flow");
				flow_algo.growReachable(cs);		//don't grow target reachable
				timer.stop("Grow Reachable due to AAP");
				cs.hasCut = true;
			}
			timer.stop("Flow");
			
			cs.verifyFlowConstraints();
			
			if (cs.hasCut) {
				cs.verifySetInvariants();
				if (cs.n.targetReachableWeight <= cs.n.sourceReachableWeight)
					flipViewDirection();
				timer.start("Grow Assimilated");
				GrowAssimilated<FlowAlgorithm>::grow(cs, flow_algo.getScanList());
				timer.stop("Grow Assimilated");
				cs.verifyFlowConstraints();
				cs.verifySetInvariants();
				LOGGER << cs.toString();
			}
			
			return true;
		}

		bool runUntilBalancedOrFlowBoundExceeded(const Node s, const Node t) {
			cs.initialize(s,t);
			bool piercingFailedOrFlowBoundReachedWithNonAAPPiercingNode = false;
			bool has_balanced_cut = false;
			
			while (cs.flowValue <= upperFlowBound && !has_balanced_cut) {
				piercingFailedOrFlowBoundReachedWithNonAAPPiercingNode = !advanceOneFlowIteration();
				if (piercingFailedOrFlowBoundReachedWithNonAAPPiercingNode)
					break;
				has_balanced_cut = cs.hasCut && cs.isBalanced(); //no cut ==> run and don't check for balance.
			}
			
			static constexpr bool log = true;
			
			if (find_most_balanced && has_balanced_cut && cs.flowValue <= upperFlowBound) {
				auto saved_state = cs.saveState();
				NodeWeight bwd = cs.outputMostBalancedPartition();
				LOGGER << "find most balanced and has balanced cut." << V(cs.flowValue) << V(upperFlowBound) << V(bwd);
				if (bwd > 0) {	//TODO replace by a more suitable choice
					auto saved_first_partition = cs.saveState();
					cs.restoreState(saved_state);
					cs.partitionWrittenToNodeSet = false;
					
					bool changed = false;
					while (advanceOneFlowIteration(true /* reject piercing if it creates an augmenting path */))
						changed = true;
					
					LOGGER << V(changed);
					if (!changed || !cs.isBalanced() || cs.outputMostBalancedPartition() > bwd) {
						cs.restoreState(saved_first_partition);
						cs.partitionWrittenToNodeSet = true;
						LOGGER << "restore first partition";
					}
					else {
						LOGGER << "mbmc worked";
					}
				}
			}
			
			Assert(!cs.hasCut || cs.isBalanced() || cs.flowValue > upperFlowBound);
			return !piercingFailedOrFlowBoundReachedWithNonAAPPiercingNode && cs.flowValue <= upperFlowBound && cs.hasCut && has_balanced_cut;
		}
		
		void findCutsUntilBalancedOrFlowBoundExceeded(const Node s, const Node t) {
			cs.initialize(s,t);
			exhaustFlowAndGrow();
			while (cs.flowValue <= upperFlowBound && !cs.isBalanced() && pierce())
				exhaustFlowAndGrow();
		}

	};


	class MultiCutter {
		//implement interleaving and parallelization here, if desired.
	};

}
