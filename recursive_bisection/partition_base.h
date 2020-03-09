#pragma once

namespace whfc_rb {
    class PartitionBase {
    public:
        using PartitionID = uint32_t;

        struct PartitionChangeElement {
            NodeID u;
            PartitionID newPart;
        };

        static constexpr PartitionID invalidPartition = std::numeric_limits<PartitionID>::max();

        PartitionBase(PartitionID num_parts, CSRHypergraph &hg) : partition(hg.numNodes(), 0), hg(hg),
                                                                  num_parts(num_parts) {}

        PartitionBase(std::vector<PartitionID> vec_partition, PartitionID num_parts, CSRHypergraph &hg) : partition(
                std::move(vec_partition)), hg(hg), num_parts(num_parts) {}

        virtual ~PartitionBase() = default;

        std::size_t size() const {
            return partition.size();
        }

        NodeWeight totalWeight() const {    // Note(Lars): that's not really a property of a partition. Also: compute once and store. You can propagate it through the Recursive Bisection
            NodeWeight weight = 0;

            for (NodeID u : hg.nodes()) {
                weight += hg.nodeWeight(u);
            }

            return weight;
        }

        NodeWeight partWeight(PartitionID id) const {
            NodeWeight weight(0);

            for (NodeID u : hg.nodes()) {
                if (partition[u] == id) { weight += hg.nodeWeight(u); }
            }

            return weight;
        }

        const std::vector<NodeWeight> partitionWeights() const {
            std::vector<NodeWeight> weights(num_parts, 0);

            for (NodeID u : hg.nodes()) {
                weights[partition[u]] += hg.nodeWeight(u);
            }

            return weights;
        }

        virtual void replace(std::vector<PartitionID> vec_part) {
            // asserts?
            partition = std::move(vec_part);
        }

        const PartitionID &operator[](std::size_t idx) const { return partition[idx]; }

        PartitionID *data() { return partition.data(); }

        std::size_t pinsInPart(PartitionID id, HyperedgeID e) const {
            std::size_t count = 0;
            for (NodeID u : hg.pinsOf(e)) {
                if (partition[u] == id) {
                    count++;
                }
            }
            return count;
        }

        virtual std::vector<HyperedgeID> getCutEdges(PartitionID part0, PartitionID part1) {
            std::vector<HyperedgeID> cut_hes;

            for (HyperedgeID e : hg.hyperedges()) {
                std::array<bool, 2> contains_part = {false, false};
                for (NodeID u : hg.pinsOf(e)) {
                    if (partition[u] == part0) { contains_part[0] = true; }
                    if (partition[u] == part1) { contains_part[1] = true; }
                }
                if (contains_part[0] && contains_part[1]) {
                    cut_hes.push_back(e);
                }
            }

            return cut_hes;
        }

        void changePart(NodeID u, PartitionID newPart) {
            partition[u] = newPart;
        }

        virtual double imbalance() {
            std::vector<NodeWeight> vec_partitionWeights = partitionWeights();
            NodeWeight totalWeight = std::accumulate(vec_partitionWeights.begin(), vec_partitionWeights.end(), 0U);
            NodeWeight maxPartWeight = *std::max_element(vec_partitionWeights.begin(), vec_partitionWeights.end());
            return (static_cast<double>(maxPartWeight) * static_cast<double>(num_parts) /
                    static_cast<double>(totalWeight)) - 1.0;
        }

        virtual size_t km1Objective() const {
            size_t obj = 0;
            boost::dynamic_bitset<> has_pins_in_part(num_parts);
            for (HyperedgeID e : hg.hyperedges()) {
                for (NodeID u : hg.pinsOf(e)) {
                    has_pins_in_part.set(partition[u]);
                }
                obj += (has_pins_in_part.count() - 1) * hg.hyperedgeWeight(e);
                has_pins_in_part.reset();
            }
            return obj;
        }

        virtual size_t cutObjective() const {
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

        CSRHypergraph &getGraph() {
            return hg;
        }

        // Note (Lars): this should be HyperedgeWeight. occurs in some other places as well.
        NodeWeight km1AfterChanges(std::vector<PartitionChangeElement> &partChanges) const {
            NodeWeight newkm1 = km1Objective();
            for (PartitionChangeElement partChange : partChanges) {
                km1Update(newkm1, partChange);
            }
            return newkm1;
        }

        virtual void print(std::ostream &out) {
            out << "Partition: ";
            for (uint i = 0; i < partition.size(); ++i) {
                out << partition[i] << " ";
            }
            out << std::endl << std::flush;
        }

        virtual void initialize() {}

    protected:
        std::vector<PartitionID> partition;
        CSRHypergraph &hg;
        PartitionID num_parts;

    private:
        void km1Update(NodeWeight &value, PartitionChangeElement &partChange) const {
            if (partChange.newPart != partition[partChange.u]) {
                NodeID oldPart = partition[partChange.u];
                NodeID newPart = partChange.newPart;
                for (HyperedgeID e : hg.hyperedgesOf(partChange.u)) {
                    if (pinsInPart(oldPart, e) == 1 &&
                        pinsInPart(newPart, e) > 0) {
                        value -= hg.hyperedgeWeight(e);
                    } else if (pinsInPart(oldPart, e) > 1 &&
                               pinsInPart(newPart, e) == 0) {
                        value += hg.hyperedgeWeight(e);
                    }
                }
            }
        }

    };
}