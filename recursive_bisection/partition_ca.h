#pragma once

#include <boost/dynamic_bitset.hpp>
#include "hypergraph.h"
#include "partition_base.h"

namespace whfc_rb {

    class PartitionCA : public PartitionBase {
    public:
        static constexpr PartitionID invalidPartition = std::numeric_limits<PartitionID>::max();

        explicit PartitionCA(PartitionID num_parts, CSRHypergraph &hg) : PartitionBase(num_parts, hg),
                                                                         vec_pinsInPart(hg.numHyperedges() * num_parts,
                                                                                        0),
                                                                         vec_partWeights(num_parts, 0) {
            assert(maxID() < num_parts);
            assert(partition.size() == hg.numNodes());
        }

        explicit PartitionCA(std::vector<PartitionID> vec_partition, PartitionID num_parts, CSRHypergraph &hg)
                : PartitionBase(vec_partition, num_parts, hg),
                  vec_pinsInPart(hg.numHyperedges() * num_parts, 0),
                  vec_partWeights(num_parts, 0) {
            assert(maxID() < num_parts);
            assert(partition.size() == hg.numNodes());

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

        const std::vector<NodeWeight> partitionWeights() const {
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
            return vec_pinsInPart[id * hg.numHyperedges() + e];
        }

        std::vector<HyperedgeID> getCutEdges(PartitionID part0, PartitionID part1) const override {
            assert(datastructures_initialized);
            std::vector<HyperedgeID> cut_hes;

            for (HyperedgeID e : hg.hyperedges()) {
                if (pinsInPart(part0, e) > 0 && pinsInPart(part1, e) > 0) {
                    cut_hes.push_back(e);
                }
            }

            return cut_hes;
        }

        void changePart(NodeID u, PartitionID newPart) {
            assert(datastructures_initialized);
            if (partition[u] != newPart) {
                PartitionID oldPart = partition[u];
                partition[u] = newPart;

                for (HyperedgeID e : hg.hyperedgesOf(u)) {
                    if (vec_pinsInPart[oldPart * hg.numHyperedges() + e] == 1 &&
                        vec_pinsInPart[newPart * hg.numHyperedges() + e] > 0) {
                        km1 -= hg.hyperedgeWeight(e);
                    } else if (vec_pinsInPart[oldPart * hg.numHyperedges() + e] > 1 &&
                               vec_pinsInPart[newPart * hg.numHyperedges() + e] == 0) {
                        km1 += hg.hyperedgeWeight(e);
                    }

                    vec_pinsInPart[oldPart * hg.numHyperedges() + e]--;
                    vec_pinsInPart[newPart * hg.numHyperedges() + e]++;
                }
                vec_partWeights[oldPart] -= hg.nodeWeight(u);
                vec_partWeights[newPart] += hg.nodeWeight(u);
            }
        }


        NodeWeight km1AfterChanges(std::vector<PartitionChangeElement> &partChanges) const {
            NodeWeight newkm1 = km1;
            for (PartitionChangeElement partChange : partChanges) {
                km1Update(newkm1, partChange);
            }
            return newkm1;
        }

        size_t km1Objective() const override {
            assert(km1 == PartitionBase::km1Objective());
            return km1;
        }

        double imbalance() override {
            assert(datastructures_initialized);
            NodeWeight totalWeight = std::accumulate(vec_partWeights.begin(), vec_partWeights.end(), 0U);
            NodeWeight maxPartWeight = *std::max_element(vec_partWeights.begin(), vec_partWeights.end());
            return (static_cast<double>(maxPartWeight) * static_cast<double>(num_parts) /
                    static_cast<double>(totalWeight)) - 1.0;
        }

        void initialize() override {
            vec_pinsInPart = std::vector<std::size_t>(hg.numHyperedges() * num_parts, 0);
            vec_partWeights = std::vector<NodeWeight>(num_parts, 0);
            km1 = 0;

            boost::dynamic_bitset<> has_pins_in_part(num_parts);
            for (HyperedgeID e : hg.hyperedges()) {
                for (NodeID u : hg.pinsOf(e)) {
                    has_pins_in_part.set(partition[u]);
                    vec_pinsInPart[partition[u] * hg.numHyperedges() + e]++;
                }
                km1 += (has_pins_in_part.count() - 1) * hg.hyperedgeWeight(e);
                has_pins_in_part.reset();
            }

            for (NodeID u : hg.nodes()) {
                vec_partWeights[partition[u]] += hg.nodeWeight(u);
            }
            datastructures_initialized = true;
        }

    private:
        std::vector<std::size_t> vec_pinsInPart;
        std::vector<NodeWeight> vec_partWeights;
        NodeWeight km1 = 0;
        bool datastructures_initialized = false;


        void km1Update(NodeWeight &value, PartitionChangeElement &partChange) const {
            if (partChange.newPart != partition[partChange.u]) {
                NodeID oldPart = partition[partChange.u];
                NodeID newPart = partChange.newPart;
                for (HyperedgeID e : hg.hyperedgesOf(partChange.u)) {
                    if (vec_pinsInPart[oldPart * hg.numHyperedges() + e] == 1 &&
                        vec_pinsInPart[newPart * hg.numHyperedges() + e] > 0) {
                        value -= hg.hyperedgeWeight(e);
                    } else if (vec_pinsInPart[oldPart * hg.numHyperedges() + e] > 1 &&
                               vec_pinsInPart[newPart * hg.numHyperedges() + e] == 0) {
                        value += hg.hyperedgeWeight(e);
                    }
                }
            }
        }

        PartitionID maxID() {
            return *std::max_element(partition.begin(), partition.end());
        }
    };
}
