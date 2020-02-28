#pragma once

#include <boost/dynamic_bitset.hpp>
#include "hypergraph.h"

namespace whfc_rb {
    class Partition {
    public:
        using PartitionID = uint32_t;
        static constexpr PartitionID invalidPartition = std::numeric_limits<PartitionID>::max();

        explicit Partition(PartitionID num_parts, CSRHypergraph& hg) : partition(hg.numNodes(), 0), hg(hg), num_parts(num_parts), vec_pinsInPart(hg.numHyperedges(), std::vector<std::size_t>(num_parts, 0)), vec_partWeights(num_parts, 0) {

        }
        explicit Partition(std::vector<PartitionID> vec_partition, PartitionID num_parts, CSRHypergraph& hg) : partition(std::move(vec_partition)), hg(hg), num_parts(num_parts), vec_pinsInPart(hg.numHyperedges(), std::vector<std::size_t>(num_parts, 0)), vec_partWeights(num_parts, 0) {
            // insert assertions here
            recompute();
        }

        std::size_t size() const {
            return partition.size();
        }

        NodeWeight totalWeight() const {
            return std::accumulate(vec_partWeights.begin(), vec_partWeights.end(), 0U);
        }

        NodeWeight partWeight(PartitionID id) const {
            return vec_partWeights[id];
        }

        const std::vector<NodeWeight>& partitionWeights() const {
            return vec_partWeights;
        }

        void rebuild(std::vector<PartitionID> vec_part) {
            // TODO: assertions here
            partition = std::move(vec_part);
            recompute();
        }


        const PartitionID& operator[](std::size_t idx) const { return partition[idx]; }

        PartitionID* data() { return partition.data(); }

        /*const_range<PartitionID> entries() {
            return {partition.begin(), partition.end()};
        }*/

        std::size_t pinsInPart(PartitionID id, HyperedgeID e) const {
            return vec_pinsInPart[e][id];
        }

        std::vector<HyperedgeID> getCutEdges(PartitionID part0, PartitionID part1) const {
            std::vector<HyperedgeID> cut_hes;

            for (HyperedgeID e : hg.hyperedges()) {
                if (pinsInPart(part0, e) > 0 && pinsInPart(part1, e) > 0) {
                    cut_hes.push_back(e);
                }
            }

            return cut_hes;
        }

        void changePart(NodeID u, PartitionID newPart) {
            if (partition[u] != newPart) {
                PartitionID oldPart = partition[u];
                partition[u] = newPart;

                for (HyperedgeID e : hg.hyperedgesOf(u)) {
                    vec_pinsInPart[e][oldPart]--;
                    vec_pinsInPart[e][newPart]++;
                }
                vec_partWeights[oldPart] -= hg.nodeWeight(u);
                vec_partWeights[newPart] += hg.nodeWeight(u);
            }
        }

        double imbalance() {
            NodeWeight totalWeight = std::accumulate(vec_partWeights.begin(), vec_partWeights.end(), 0U);
            NodeWeight maxPartWeight = *std::max_element(vec_partWeights.begin(), vec_partWeights.end());
            return (static_cast<double>(maxPartWeight) * static_cast<double>(num_parts) / static_cast<double>(totalWeight)) - 1.0;
        }
        
        size_t km1Objective() const {
            size_t obj = 0;
            boost::dynamic_bitset<> has_pins_in_part(num_parts);
            for (HyperedgeID e : hg.hyperedges()) {
                for (NodeID u : hg.pinsOf(e)) {
                    has_pins_in_part.set( partition[u] );
                }
                obj += (has_pins_in_part.count() - 1) * hg.hyperedgeWeight(e);
                has_pins_in_part.reset();
            }
            return obj;
        }
        
        size_t cutObjective() const {
            size_t obj = 0;
            for (HyperedgeID e : hg.hyperedges()) {
                PartitionID p = invalidPartition;
                for (NodeID u : hg.pinsOf(e)) {
                    if (p == invalidPartition) {
                        p = partition[u];
                    } else if (partition[u] != p) {
                        obj += hg.hyperedgeWeight(e);
                        break;
                    }
                }
            }
            return obj;
        }

        PartitionID numParts() {
            return num_parts;
        }

        CSRHypergraph& getGraph() {
            return hg;
        }

        void print(std::ostream& out) {
            out << "Partition: ";
            for (uint i = 0; i < partition.size(); ++i) {
                out << partition[i] << " ";
            }
            out << std::endl << std::flush;
        }

        void recompute() {
            vec_pinsInPart = std::vector<std::vector<std::size_t>>(hg.numHyperedges(), std::vector<std::size_t>(num_parts, 0));
            vec_partWeights = std::vector<NodeWeight>(num_parts, 0);

            for (HyperedgeID e : hg.hyperedges()) {
                for (NodeID u : hg.pinsOf(e)) {
                    vec_pinsInPart[e][partition[u]]++;
                }
            }

            for (NodeID u : hg.nodes()) {
                vec_partWeights[partition[u]] += hg.nodeWeight(u);
            }
        }

    private:
        std::vector<PartitionID> partition;
        CSRHypergraph& hg;
        PartitionID num_parts;
        std::vector<std::vector<std::size_t>> vec_pinsInPart;
        std::vector<NodeWeight> vec_partWeights;

    };
}
