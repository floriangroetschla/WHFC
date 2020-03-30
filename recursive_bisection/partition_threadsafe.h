#pragma once

#include <tbb/enumerable_thread_specific.h>
#include <boost/dynamic_bitset.hpp>

#include "hypergraph.h"
#include "partition_base.h"

#include "timestamp_set.hpp"

#include "../util/concatenated_range.h"
#include "../util/range.h"


namespace whfc_rb {

    class CutEdgeRange : public concatenated_range<const_range<std::vector<HyperedgeID>>, HyperedgeID> {
    public:
        using Partial = const_range<std::vector<HyperedgeID>>;
        using Base = concatenated_range<Partial, HyperedgeID>;

        CutEdgeRange(std::vector<HyperedgeID>& c0, std::vector<HyperedgeID>& c1) : Base(Partial(c0), Partial(c1)), c0(c0), c1(c1) { }

        template<typename URBG>
        void shuffle(URBG&& urbg) {
            std::shuffle(c0.begin(), c0.end(), urbg);
            std::shuffle(c1.begin(), c1.end(), urbg);
        }

    private:
        std::vector<HyperedgeID>& c0, c1;
    };

    class PartitionThreadsafe : public PartitionBase {
    public:
        static constexpr PartitionID invalidPartition = std::numeric_limits<PartitionID>::max();

        explicit PartitionThreadsafe(PartitionID num_parts, CSRHypergraph &hg) :
                    PartitionBase(num_parts, hg),
                    ets_deduplicator(hg.numHyperedges()), vec_pinsInPart(hg.numHyperedges() * num_parts, 0),
                    vec_partWeights(num_parts, 0), vec_cutEdges(num_parts * num_parts)
        {
            assert(maxID() < num_parts);
            assert(partition.size() == hg.numNodes());
        }

        explicit PartitionThreadsafe(std::vector<PartitionID> vec_partition, PartitionID num_parts, CSRHypergraph &hg) :
                PartitionThreadsafe(num_parts, hg)
        {
            initialize();
        }

        NodeWeight totalWeight() const {
            assert(datastructures_initialized);
            return std::accumulate(vec_partWeights.begin(), vec_partWeights.end(), 0U);
        }

        NodeWeight partWeight(PartitionID id) const {
            assert(datastructures_initialized);
            return vec_partWeights[id];
        }

        std::vector<NodeWeight> partitionWeights() const {
            assert(datastructures_initialized);
            return vec_partWeights;
        }

        void replace(std::vector<PartitionID> vec_part) override {
            assert(vec_part.size() == partition.size());
            assert(*std::max_element(vec_part.begin(), vec_part.end()) < num_parts);

            partition = std::move(vec_part);
            datastructures_initialized = false;
        }

        std::size_t pinsInPart(PartitionID id, HyperedgeID e) const {
            assert(datastructures_initialized);
            return vec_pinsInPart[e * num_parts + id];
        }

        void filter(std::vector<HyperedgeID>& cut_hes, ldc::TimestampSet<>& deduplicator, PartitionID part0, PartitionID part1) {
            for (size_t i = 0; i < cut_hes.size(); ++i) {
                const HyperedgeID e = cut_hes[i];
                if (deduplicator.contains(e) || pinsInPart(part0, e) == 0 || pinsInPart(part1, e) == 0) {
                    cut_hes[i] = cut_hes.back();
                    cut_hes.pop_back();
                    --i;
                }
                deduplicator.add(e);
            }
        }

        CutEdgeRange getCutEdges(PartitionID part0, PartitionID part1) {
            assert(datastructures_initialized);
            ldc::TimestampSet<>& deduplicator = ets_deduplicator.local();
            deduplicator.clear();
            std::vector<HyperedgeID>& cut_hes1 = cutEdges(part0, part1);
            std::vector<HyperedgeID>& cut_hes2 = cutEdges(part1, part0);
            filter(cut_hes1, deduplicator, part0, part1);
            filter(cut_hes2, deduplicator, part0, part1);
            return CutEdgeRange(cut_hes1, cut_hes2);
        }

        void changePart(NodeID u, PartitionID newPart) {
            assert(datastructures_initialized);
            if (partition[u] != newPart) {
                PartitionID oldPart = partition[u];
                partition[u] = newPart;

                for (HyperedgeID e : hg.hyperedgesOf(u)) {
                    // Update cuts
                    if (pinsInPart(newPart, e) == 0) {
                        for (PartitionID otherPart = 0; otherPart < num_parts; ++otherPart) {
                            if (otherPart != newPart && pinsInPart(otherPart, e) > 0) {
                                cutEdges(newPart, otherPart).push_back(e);
                            }
                        }
                    }

                    pinsInPart(oldPart, e)--;
                    pinsInPart(newPart, e)++;
                }
                vec_partWeights[oldPart] -= hg.nodeWeight(u);
                vec_partWeights[newPart] += hg.nodeWeight(u);
            }
        }

        std::vector<HyperedgeID>& cutEdges(PartitionID part0, PartitionID part1) {
            assert(part0 != part1 && part0 < num_parts && part1 < num_parts);
            return (vec_cutEdges[(part0 * num_parts) + part1]);
        }

        double imbalance() override {
            assert(datastructures_initialized);
            auto total_weight = static_cast<double>(totalWeight());
            auto max_part_weight = static_cast<double>(*std::max_element(vec_partWeights.begin(), vec_partWeights.end()));
            return (max_part_weight * num_parts / total_weight) - 1.0;
        }

        void initialize() override {
            // TODO parallelize!

            std::vector<PartitionID> partsOfHyperedge;
            for (HyperedgeID e : hg.hyperedges()) {
                for (NodeID u : hg.pinsOf(e)) {
                    pinsInPart(partition[u], e)++;
                }

                // Note: complexity can be reduced for large k and small hyperedges by using a clearlist in the loop over the pins, if that turns out to be an issue

                for (PartitionID p = 0; p < numParts(); ++p) {
                    if (pinsInPart(p, e) > 0) {
                        partsOfHyperedge.push_back(p);
                    }
                }

                for (size_t i = 0; i < partsOfHyperedge.size() - 1; ++i) {
                    for (size_t j = i + 1; j < partsOfHyperedge.size(); ++j) {
                        cutEdges(partsOfHyperedge[i], partsOfHyperedge[j]).push_back(e);
                    }
                }
                partsOfHyperedge.clear();
            }

            for (NodeID u : hg.nodes()) {
                vec_partWeights[partition[u]] += hg.nodeWeight(u);
            }

            datastructures_initialized = true;
        }

    private:

        tbb::enumerable_thread_specific<ldc::TimestampSet<>> ets_deduplicator;
        std::vector<uint32_t> vec_pinsInPart;
        std::vector<NodeWeight> vec_partWeights;
        std::vector<std::vector<HyperedgeID>> vec_cutEdges;
        bool datastructures_initialized = false;

        uint32_t& pinsInPart(PartitionID id, HyperedgeID e) {
            return vec_pinsInPart[e * num_parts + id];
        }

        PartitionID maxID() {
            return *std::max_element(partition.begin(), partition.end());
        }
    };
}
