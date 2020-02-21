#pragma once

#include <boost/dynamic_bitset.hpp>
#include "hypergraph.h"

namespace whfc_rb {
    class Partition {
    public:
        using PartitionID = uint32_t;
        static constexpr PartitionID invalidPartition = std::numeric_limits<PartitionID>::max();

        explicit Partition() : partition(), num_parts(1) {}
        explicit Partition(std::size_t size, PartitionID num_parts) : partition(size, 0), num_parts(num_parts) {}
        explicit Partition(std::vector<PartitionID> vec_partition, PartitionID num_parts) : partition(std::move(vec_partition)), num_parts(num_parts) {}


        std::size_t partSize(CSRHypergraph& hg, PartitionID id) const {
            if (id >= num_parts) throw std::runtime_error("partition id out of range");

            std::size_t size = 0;
            for (CSRHypergraph::NodeID u : hg.nodes()) {
                if (partition[u] == id) { size++; }
            }
            return size;
        }

        std::vector<std::size_t> partSizes(CSRHypergraph& hg) const {
            std::vector<std::size_t> sizes(num_parts, 0);

            for (CSRHypergraph::NodeID u : hg.nodes()) {
                sizes[partition[u]]++;
            }

            return sizes;
        }

        std::size_t size() const {
            return partition.size();
        }

        CSRHypergraph::NodeWeight weight(CSRHypergraph& hg, PartitionID id) const {
            CSRHypergraph::NodeWeight weight(0);

            for (CSRHypergraph::NodeID u : hg.nodes()) {
                if (partition[u] == id) { weight += hg.nodeWeight(u); }
            }

            return weight;
        }

        std::vector<CSRHypergraph::NodeWeight> partitionWeights(CSRHypergraph& hg) const {
            std::vector<CSRHypergraph::NodeWeight> weights(num_parts, 0);

            for (CSRHypergraph::NodeID u : hg.nodes()) {
                weights[partition[u]] += hg.nodeWeight(u);
            }

            return weights;
        }

        std::vector<CSRHypergraph::NodeID> nodesOf(CSRHypergraph& hg, PartitionID id) {
            if (id >= num_parts) throw std::runtime_error("partition id out of range");

            std::vector<CSRHypergraph::NodeID> nodes;
            for (CSRHypergraph::NodeID u : hg.nodes()) {
                if (partition[u] == id) {
                    nodes.push_back(u);
                }
            }

            return nodes;
        }

        PartitionID& operator[](std::size_t idx) { return partition[idx]; }

        const PartitionID& operator[](std::size_t idx) const { return partition[idx]; }

        PartitionID* data() { return partition.data(); }

        /*const_range<PartitionID> entries() {
            return {partition.begin(), partition.end()};
        }*/

        std::size_t pinsInPart(CSRHypergraph& hg, PartitionID id, CSRHypergraph::HyperedgeID e) const {
            std::size_t count = 0;
            for (CSRHypergraph::NodeID u : hg.pinsOf(e)) {
                if (partition[u] == id) {
                    count++;
                }
            }
            return count;
        }

        std::vector<CSRHypergraph::HyperedgeID> getCutEdges(CSRHypergraph& hg, PartitionID part0, PartitionID part1) const {
            std::vector<CSRHypergraph::HyperedgeID> cut_hes;

            for (CSRHypergraph::HyperedgeID e : hg.hyperedges()) {
                std::array<bool, 2> contains_part = {false, false};
                for (CSRHypergraph::NodeID u : hg.pinsOf(e)) {
                    if (partition[u] == part0) { contains_part[0] = true; }
                    if (partition[u] == part1) { contains_part[1] = true; }
                }
                if (contains_part[0] && contains_part[1]) {
                    cut_hes.push_back(e);
                }
            }

            return cut_hes;
        }

        double imbalance(CSRHypergraph& hg) {
            std::vector<CSRHypergraph::NodeWeight> vec_partitionWeights = partitionWeights(hg);
            CSRHypergraph::NodeWeight totalWeight = std::accumulate(vec_partitionWeights.begin(), vec_partitionWeights.end(), 0U);
            CSRHypergraph::NodeWeight maxPartWeight = *std::max_element(vec_partitionWeights.begin(), vec_partitionWeights.end());
            return (static_cast<double>(maxPartWeight) * static_cast<double>(num_parts) / static_cast<double>(totalWeight)) - 1.0;
        }
        
        size_t km1Objective(const CSRHypergraph& hg) const {
            size_t obj = 0;
            boost::dynamic_bitset<> has_pins_in_part(num_parts);
            for (CSRHypergraph::HyperedgeID e : hg.hyperedges()) {
                for (CSRHypergraph::NodeID u : hg.pinsOf(e)) {
                    has_pins_in_part.set( partition[u] );
                }
                obj += (has_pins_in_part.count() - 1) * hg.hyperedgeWeight(e);
                has_pins_in_part.reset();
            }
            return obj;
        }
        
        size_t cutObjective(const CSRHypergraph& hg) const {
            size_t obj = 0;
            for (CSRHypergraph::HyperedgeID e : hg.hyperedges()) {
                PartitionID p = invalidPartition;
                for (CSRHypergraph::NodeID u : hg.pinsOf(e)) {
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

        void setNumParts(PartitionID num) {
            num_parts = num;
        }

        void print(std::ostream& out) {
            out << "Partition: ";
            for (uint i = 0; i < partition.size(); ++i) {
                out << partition[i] << " ";
            }
            out << std::endl << std::flush;
        }

    private:
        std::vector<PartitionID> partition;
        PartitionID num_parts;

    };
}
