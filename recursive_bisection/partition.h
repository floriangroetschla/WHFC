#pragma once

#include "hypergraph.h"

namespace whfc_rb {
    class Partition {
    public:
        using partitionID = uint32_t;

        explicit Partition() : partition(), num_parts(1) {}
        explicit Partition(std::size_t size, partitionID num_parts) : partition(size, 0), num_parts(num_parts) {}
        explicit Partition(std::vector<partitionID> vec_partition, partitionID num_parts) : partition(std::move(vec_partition)), num_parts(num_parts) {}


        std::size_t partSize(CSRHypergraph& hg, partitionID id) const {
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

        CSRHypergraph::NodeWeight weight(CSRHypergraph& hg, partitionID id) const {
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

        std::vector<CSRHypergraph::NodeID> nodesOf(CSRHypergraph& hg, partitionID id) {
            if (id >= num_parts) throw std::runtime_error("partition id out of range");

            std::vector<CSRHypergraph::NodeID> nodes;
            for (CSRHypergraph::NodeID u : hg.nodes()) {
                if (partition[u] == id) {
                    nodes.push_back(u);
                }
            }

            return nodes;
        }

        partitionID& operator[](std::size_t idx) { return partition[idx]; }

        const partitionID& operator[](std::size_t idx) const { return partition[idx]; }

        partitionID* data() { return partition.data(); }

        /*const_range<partitionID> entries() {
            return {partition.begin(), partition.end()};
        }*/

        std::size_t pinsInPart(CSRHypergraph& hg, partitionID id, CSRHypergraph::HyperedgeID e) const {
            std::size_t count = 0;
            for (CSRHypergraph::NodeID u : hg.pinsOf(e)) {
                if (partition[u] == id) {
                    count++;
                }
            }
            return count;
        }

        std::vector<CSRHypergraph::HyperedgeID> getCutEdges(CSRHypergraph& hg, partitionID part0, partitionID part1) const {
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
            CSRHypergraph::NodeWeight totalWeight = 0;
            CSRHypergraph::NodeWeight maxPartWeight = 0;

            for (uint i = 0; i < vec_partitionWeights.size(); ++i) {
                if (vec_partitionWeights[i] > maxPartWeight) {
                    maxPartWeight = vec_partitionWeights[i];
                }
                totalWeight += vec_partitionWeights[i];
            }

            return (static_cast<double>(maxPartWeight) / static_cast<double>(totalWeight)) - (1.0 / static_cast<double>(num_parts));
        }

        void print(std::ostream& out) {
            out << "Partition: ";
            for (uint i = 0; i < partition.size(); ++i) {
                out << partition[i] << " ";
            }
            out << std::endl << std::flush;
        }

    private:
        std::vector<partitionID> partition;
        partitionID num_parts;

    };
}
