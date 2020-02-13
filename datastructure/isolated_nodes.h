#pragma once

#include "../definitions.h"
#include "flow_hypergraph.h"
#include "../datastructure/bitvector.h"

namespace whfc {
	class IsolatedNodes {
	private:
		FlowHypergraph& hg;
		bool useIsolatedNodes = true;
	public:
		NodeWeight weight = NodeWeight(0);
		std::vector<Node> nodes;
		std::vector<InHeIndex> mixedIncidentHyperedges;
		BitVector hasSettledSourcePins, hasSettledTargetPins;

		struct SummableRange {
			NodeWeight from, to;
			SummableRange(NodeWeight _from, NodeWeight _to) : from(_from), to(_to) { AssertMsg(_from != NodeWeight::Invalid() && _to != NodeWeight::Invalid(), "Invalid range"); }
			bool inRange(const NodeWeight w) const { return from <= w && w <= to; }
			bool operator<(const SummableRange& o) const { return std::tie(from, to) < std::tie(o.from, o.to); }
			bool operator==(const SummableRange& o) const { return from == o.from && to == o.to; }
			bool operator!=(const SummableRange& o) const { return !operator==(o); }
		};

	private:
		NodeWeight maxSubsetSumWeight = NodeWeight(0);

		struct TableEntry {
			Node node;
			Index sumsIndex;
			TableEntry() : node(invalidNode), sumsIndex(invalidIndex) { }
			bool summable() const {
						//sumsIndex = 0 is reserved for the first range.
				AssertMsg(sumsIndex == invalidIndex || node != invalidNode || sumsIndex == 0, "Summable but no node");
				return sumsIndex != invalidIndex;
			}
		};

		std::vector<TableEntry> DPTable;
		std::vector<SummableRange> sumRanges;
		std::vector<SummableRange> nextSumRanges;
		bool newSumAvailable = true;
		std::vector<Node> nodesNotInTheDPTable;
		
		// This variant assumes that DPTable[sumRanges[i].from] == DPTable[sumRanges[i].to] == i
		// Pointers in between are considered stale.
		// Updating DP Table entries requires computing all sums available with node u.
		void updateDPTableWithSumRanges() {
			Assert(useIsolatedNodes);
			for (const Node u : nodesNotInTheDPTable) {
				if (newSumAvailable)
					nextSumRanges = sumRanges;
				AssertMsg(nextSumRanges == sumRanges, "nextSumRanges not up to date");
				newSumAvailable = false;

				const NodeWeight wu = hg.nodeWeight(u);
				AssertMsg(wu > 0, "Node has zero weight");
				for (const SummableRange& sr : sumRanges) {
					for (NodeWeight new_sum = sr.from + wu, _end = std::min(sr.to + wu, maxSubsetSumWeight); new_sum <= _end; ++new_sum) {
						if (!isSummable(new_sum)) {
							newSumAvailable = true;
							DPTable[new_sum].node = u;

							NodeWeight left(new_sum - 1), right(new_sum + 1); //new_sum + 1 is valid because of right-ward sentinel.
							Index leftIndex = DPTable[left].sumsIndex, rightIndex = DPTable[right].sumsIndex;
							bool hasLeft = DPTable[left].summable(), hasRight = DPTable[right].summable();

							if (hasLeft && hasRight) {
								SummableRange& leftRange = nextSumRanges[leftIndex], rightRange = nextSumRanges[rightIndex];
								
								//merge ranges. keep left range, and extend it to cover the right range
								AssertMsg(leftRange.to == left, "hasLeft && hasRight: left range does not extend to new_sum-1");
								AssertMsg(rightRange.from == right, "hasLeft && hasRight: right range does not start at new_sum+1");
								DPTable[new_sum].sumsIndex = leftIndex; //bridging cell. the index is stale, but we're setting it anyway to mark it as summable

								//extend leftRange to cover rightRange
								DPTable[rightRange.to].sumsIndex = leftIndex;
								leftRange.to = rightRange.to;

								//delete rightRange
								if (rightIndex != nextSumRanges.size() - 1) {
									//update sumsIndex of the range that was swapped to rightIndex
									nextSumRanges[rightIndex] = nextSumRanges.back();
									DPTable[ nextSumRanges[rightIndex].from ].sumsIndex = rightIndex;
									DPTable[ nextSumRanges[rightIndex].to ].sumsIndex = rightIndex;
								}
								nextSumRanges.pop_back();
							}
							else if (hasLeft) {
								//extend left range's .to by +1
								AssertMsg(nextSumRanges[leftIndex].to == left, "hasLeft: left range does not extend to new_sum-1");
								nextSumRanges[leftIndex].to = new_sum;
								DPTable[new_sum].sumsIndex = leftIndex;
							}
							else if (hasRight) {
								//extend right range's .from by -1
								AssertMsg(nextSumRanges[rightIndex].from == right, "hasRight: right range does not start at new_sum+1");
								nextSumRanges[rightIndex].from = new_sum;
								DPTable[new_sum].sumsIndex = rightIndex;
							}
							else {
								//start new range
								DPTable[new_sum].sumsIndex = static_cast<Index>(nextSumRanges.size());
								nextSumRanges.emplace_back(new_sum, new_sum);
							}
						}
					}
				}
				std::swap(nextSumRanges, sumRanges);
			}
		}

		//This variant assumes that DPTable[sumRanges[i].from] == i and DPTable[sumRanges[i].to] == sumRanges[i].from
		//Pointers in between point to -- potentially former -- left ends (sr.from) of the range sr. Using path compression we can keep "hop distance" to sr.from "short"
		//When a newly computed sum equals an already computed sum, we can jump to the left end of the range, get the corresponding index into sumRanges,
		//thus the entire range, and prune sums that have been previously computed
		void updateDPTableWithSumRangesAndRangePruning() {
			//Implement me if the other approach is too slow.
		}

	public:
		explicit IsolatedNodes(FlowHypergraph& hg, bool useIsolatedNodes, NodeWeight maxBlockWeight = NodeWeight(0)) :
				hg(hg),
				useIsolatedNodes(useIsolatedNodes),
				mixedIncidentHyperedges(useIsolatedNodes ? hg.numNodes() : 0, InHeIndex(0)),
				hasSettledSourcePins(useIsolatedNodes ? hg.numHyperedges() : 0),
				hasSettledTargetPins(useIsolatedNodes ? hg.numHyperedges() : 0),
				maxSubsetSumWeight(maxBlockWeight),
				DPTable(useIsolatedNodes ? maxBlockWeight + 2 : 0U, TableEntry())
		{
			if (useIsolatedNodes) {
				sumRanges.emplace_back(NodeWeight(0), NodeWeight(0));
				DPTable[0].sumsIndex = 0;
			}
		}
		
		void reset() {
			if (useIsolatedNodes) {
				std::fill_n(mixedIncidentHyperedges.begin(), hg.numNodes(), InHeIndex(0));
				std::fill_n(DPTable.begin(), weight + 1, TableEntry());
				sumRanges.clear();
				nextSumRanges.clear();
				newSumAvailable = true;
				nodes.clear();
				nodesNotInTheDPTable.clear();
				weight = NodeWeight(0);
				hasSettledSourcePins.reset(0, hg.numHyperedges());
				hasSettledTargetPins.reset(0, hg.numHyperedges());
				sumRanges.emplace_back(NodeWeight(0), NodeWeight(0));
				DPTable[0].sumsIndex = 0;
			}
		}
		
		void adaptMaxBlockWeight(const NodeWeight mw) {
			if (useIsolatedNodes && mw > maxSubsetSumWeight) {
				maxSubsetSumWeight = mw;
				DPTable.resize(maxSubsetSumWeight + 2, TableEntry());
			}
		}
		
		void flipViewDirection() {
			std::swap(hasSettledSourcePins, hasSettledTargetPins);
		}

		const std::vector<SummableRange>& getSumRanges() const {
			return sumRanges;
		}

		bool isSummable(const NodeWeight w) const {
			Assert(w < DPTable.size());
			return DPTable[w].summable();
		}

		void add(const Node u) {
			// TODO handle nodes with weight 1 explicitly
			
			nodes.push_back(u);
			nodesNotInTheDPTable.push_back(u);
			weight += hg.nodeWeight(u);
		}

		bool isDPTableUpToDate() const {
			return nodesNotInTheDPTable.empty();
		}

		void updateDPTable() {
			updateDPTableWithSumRanges();
			//Internal_UpdateDPTableWithSumRangesAndRangePruning();		NOT IMPLEMENTED YET. Uncertain if faster or even necessary
			nodesNotInTheDPTable.clear();
		}

		std::vector<Node> extractSubset(NodeWeight sum) const {
			AssertMsg(isSummable(sum),
					  std::string("Trying to extract subset for not achieved subset sum. ")
					  + (isDPTableUpToDate() ? "Call updateDPTable() before calling this method. There are nodes that are not included in the table yet." : " ")
			);
			
			std::vector<Node> result;
			while (sum > 0) {
				const Node u = DPTable[sum].node;
				result.push_back(u);
				sum -= hg.nodeWeight(u);
			}
			return result;
		}

		std::pair<std::vector<Node>, std::vector<Node>> extractBipartition(NodeWeight sum) const {
			auto first = extractSubset(sum);
			auto second = setDifference(nodes, first, hg.numNodes());
			return std::make_pair(std::move(first), std::move(second));
		}

		bool isCandidate(const Node u) const {
			Assert(useIsolatedNodes);
			return mixedIncidentHyperedges[u] == hg.degree(u);
		}

	};
}