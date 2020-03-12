#pragma once

#include <boost/dynamic_bitset.hpp>
#include "hypergraph.h"
#include "partition_base.h"

namespace whfc_rb {

    class PartitionThreadsafe : public PartitionBase {
    public:
        static constexpr PartitionID invalidPartition = std::numeric_limits<PartitionID>::max();
        static constexpr bool precomputeCuts = true;

        explicit PartitionThreadsafe(PartitionID num_parts, CSRHypergraph &hg) : PartitionBase(num_parts, hg),
                                                                         vec_pinsInPart(hg.numHyperedges() * num_parts),
                                                                         vec_partWeights(num_parts), vec_cutEdges(num_parts * (num_parts - 1) / 2),
                                                                         cutEdgeLocks(num_parts * (num_parts - 1) / 2) {
            assert(maxID() < num_parts);
            assert(partition.size() == hg.numNodes());
        }

        explicit PartitionThreadsafe(std::vector<PartitionID> vec_partition, PartitionID num_parts, CSRHypergraph &hg)
                : PartitionBase(vec_partition, num_parts, hg),
                  vec_pinsInPart(hg.numHyperedges() * num_parts),
                  vec_partWeights(num_parts), vec_cutEdges(num_parts * (num_parts - 1) / 2) {
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

            std::vector<NodeWeight> weights(vec_partWeights.size());
            for (uint i = 0; i < vec_partWeights.size(); ++i) {
                weights[i] = vec_partWeights[i].load();
            }

            return weights;
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

        std::vector<HyperedgeID> getCutEdges(PartitionID part0, PartitionID part1) override {
            assert(datastructures_initialized);
            std::vector<HyperedgeID> cut_hes;

            if constexpr (precomputeCuts) {
                for (HyperedgeID e : cutEdges(part0, part1)) {
                    if (pinsInPart(part0, e) > 0 && pinsInPart(part1, e) > 0 && std::find(cut_hes.begin(), cut_hes.end(), e) == cut_hes.end()) {
                        cut_hes.push_back(e);
                    }
                }
                cutEdgeLock(part0, part1).lock();
                cutEdges(part0, part1) = cut_hes;
                cutEdgeLock(part0, part1).unlock();
            } else {
                for (HyperedgeID e : hg.hyperedges()) {
                    if (pinsInPart(part0, e) > 0 && pinsInPart(part1, e) > 0) {
                        cut_hes.push_back(e);
                    }
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
                    // Update cuts
                    if constexpr (precomputeCuts) {
                        if (pinsInPart(newPart, e) == 0) {
                            for (PartitionID pid = 0; pid < num_parts; ++pid) {
                                if (pid != newPart && pinsInPart(pid, e) > 0) {
                                    cutEdgeLock(pid, newPart).lock();
                                    cutEdges(pid, newPart).push_back(e);
                                    cutEdgeLock(pid, newPart).unlock();
                                }
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

            if (part0 < part1) std::swap(part0, part1);
            return (vec_cutEdges[(part0 * (part0 - 1) / 2) + part1]);
        }

        std::mutex& cutEdgeLock(PartitionID part0, PartitionID part1) {
            assert(part0 != part1 && part0 < num_parts && part1 < num_parts);

            if (part0 < part1) std::swap(part0, part1);
            return (cutEdgeLocks[(part0 * (part0 - 1) / 2) + part1]);
        }

        double imbalance() override {
            assert(datastructures_initialized);
            NodeWeight totalWeight = std::accumulate(vec_partWeights.begin(), vec_partWeights.end(), 0U);
            NodeWeight maxPartWeight = *std::max_element(vec_partWeights.begin(), vec_partWeights.end());
            return (static_cast<double>(maxPartWeight) * static_cast<double>(num_parts) /
                    static_cast<double>(totalWeight)) - 1.0;
        }

        void initialize() override {
            for (uint i = 0; i < vec_pinsInPart.size(); ++i) {
                vec_pinsInPart[i] = 0;
            }
            for (uint i = 0; i < vec_partWeights.size(); ++i) {
                vec_partWeights[i] = 0;
            }
            vec_cutEdges = std::vector<std::vector<HyperedgeID>>(num_parts * (num_parts - 1) / 2);

            for (HyperedgeID e : hg.hyperedges()) {
                std::vector<NodeID> partitionsOfEdge;
                for (NodeID u : hg.pinsOf(e)) {
                    pinsInPart(partition[u], e)++;
                    if (std::find(partitionsOfEdge.begin(), partitionsOfEdge.end(), partition[u]) == partitionsOfEdge.end()) {
                        partitionsOfEdge.push_back(partition[u]);
                    }
                }

                if constexpr (precomputeCuts) {
                    for (uint i = 0; i < partitionsOfEdge.size() - 1; ++i) {
                        for (uint j = i + 1; j < partitionsOfEdge.size(); ++j) {
                            cutEdges(partitionsOfEdge[i], partitionsOfEdge[j]).push_back(e);
                        }
                    }
                }
            }

            for (NodeID u : hg.nodes()) {
                vec_partWeights[partition[u]] += hg.nodeWeight(u);
            }
            datastructures_initialized = true;
        }

    private:
        std::vector<std::atomic<std::size_t>> vec_pinsInPart;
        std::vector<std::atomic<NodeWeight>> vec_partWeights;
        std::vector<std::vector<HyperedgeID>> vec_cutEdges;
        std::vector<std::mutex> cutEdgeLocks;
        bool datastructures_initialized = false;

        std::atomic<std::size_t>& pinsInPart(PartitionID id, HyperedgeID e) {
            return vec_pinsInPart[id * hg.numHyperedges() + e];
        }

        PartitionID maxID() {
            return *std::max_element(partition.begin(), partition.end());
        }
    };
}
