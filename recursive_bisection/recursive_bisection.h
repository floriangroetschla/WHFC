#include <stdexcept>
#include <cmath>
#include "../extern/patoh_wrapper.h"
#include "../io/hmetis_io.h"
#include "hypergraph.h"
#include "refinement.h"


namespace whfc_rb {
    class RecursiveBisector {
    public:
        using PartitionID = int;
        static constexpr PartitionID invalidPartition = std::numeric_limits<PartitionID>::max();
        
        using NodeID = CSRHypergraph::NodeID;
        using HyperedgeID = CSRHypergraph::HyperedgeID;

        RecursiveBisector(uint maxNumNodes, uint maxNumEdges, uint maxNumPins, int seed) :
            refiner(maxNumNodes, maxNumEdges, maxNumPins, seed) {

        }

        Partition run(CSRHypergraph& hg, int seed, double epsilon, std::string preset, uint k) {
            epsilon = std::pow(1.0 + epsilon, 1.0 / std::log2(k)) - 1.0;

            return partition_recursively(hg, seed, epsilon, preset, k, true);
        }

    private:
        WHFCRefiner refiner;

        Partition partition_recursively(CSRHypergraph& hg, int seed, double epsilon, std::string preset, uint k, bool alloc) {
            if (k == 1) {
                return Partition(hg.numNodes(), 1);
            }

            std::array<int, 2> numParts;
            Partition partition;

            if (k % 2 == 0) {
                numParts[0] = k / 2;
                numParts[1] = k / 2;
                partition = PaToHInterface::bisectWithPatoh(hg, seed, epsilon, preset, alloc, false);
            } else {
                numParts[0] = k / 2;
                numParts[1] = numParts[0] + 1;
                partition = PaToHInterface::bisectImbalancedWithPatoh(hg, seed, float(numParts[1]) / float(numParts[0]), epsilon, preset, alloc, false);
            }

            partition.print(std::cout);
            bool result = refiner.refine(partition, hg, epsilon);
            std::cout << "Refinement result: " << result << std::endl;
            partition.print(std::cout);

            if (k > 2) {
                std::vector<int> new_ids(partition.size());
                std::vector<int> carries(2, 0);
                for (uint i = 0; i < partition.size(); ++i) {
                    new_ids[i] = carries[partition[i]]++;
                }

                std::array<Partition, 2> sub_partitions;

                for (uint partID = 0; partID < 2; ++partID) {
                    CSRHypergraph partHg;
                    for (HyperedgeID e : hg.hyperedges()) {
                        int hyperedgeSize = 0;
                        for (NodeID pin : hg.pinsOf(e)) {
                            if (partID == partition[pin]) {
                                partHg.pins().push_back(new_ids[pin]);
                                hyperedgeSize++;
                            }
                        }

                        if (hyperedgeSize == 1) {
                            partHg.pins().pop_back();
                        } else if (hyperedgeSize > 1) {
                            partHg.indexPins().push_back(partHg.indexPins().back() + hyperedgeSize);
                            partHg.hyperedgeWeights().push_back(hg.hyperedgeWeight(e));
                        }
                    }

                    for (uint i = 0; i < partition.size(); ++i) {
                        if (partition[i] == partID) {
                            partHg.nodeWeights().push_back(hg.nodeWeight(i));
                        }
                    }

                    partHg.initNodes(carries[partID]);
                    partHg.computeXPins();

                    sub_partitions[partID] = partition_recursively(partHg, seed, epsilon, preset, numParts[partID], false);

                }

                for (uint i = 0; i < partition.size(); ++i) {
                    partition[i] = sub_partitions[partition[i]][new_ids[i]] + partition[i] * numParts[0];
                }
            }
            
			if (alloc) {
				PaToHInterface::freePatoh();
			}
			return partition;
        }
    };
}