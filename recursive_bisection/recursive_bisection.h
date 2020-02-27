#include <stdexcept>
#include <cmath>
#include "../extern/patoh_wrapper.h"
#include "../io/hmetis_io.h"
#include "hypergraph.h"
#include "refinement.h"
#include <random>

namespace whfc_rb {
    class RecursiveBisector {
    public:

        RecursiveBisector(uint maxNumNodes, uint maxNumEdges, uint maxNumPins, std::mt19937& mt, whfc::TimeReporter& timer) :
            refiner(maxNumNodes, maxNumEdges, maxNumPins, mt, timer), mt(mt), timer(timer) {

        }

        Partition run(CSRHypergraph& hg, double epsilon, std::string preset, uint k) {
            epsilon = std::pow(1.0 + epsilon, 1.0 / std::ceil(std::log2(k))) - 1.0;

            Partition partition(k, hg);
            partition_recursively(partition, epsilon, preset, k, true);
            return partition;
        }

    private:
        WHFCRefiner refiner;
        std::mt19937& mt;
        whfc::TimeReporter& timer;

        void partition_recursively(Partition& partition, double epsilon, std::string preset, uint k, bool alloc) {
            // insert assertions here
            if (k == 1) {
                return;
            }

            std::array<int, 2> numParts;
            CSRHypergraph& hg = partition.getGraph();

            timer.start("PaToH", "RecursiveBisector");
            if (k % 2 == 0) {
                numParts[0] = k / 2;
                numParts[1] = k / 2;
                PaToHInterface::bisectWithPatoh(partition, mt(), epsilon, preset, alloc, false);
            } else {
                numParts[0] = k / 2;
                numParts[1] = numParts[0] + 1;
                PaToHInterface::bisectImbalancedWithPatoh(partition, mt(), float(numParts[1]) / float(numParts[0]), epsilon, preset, alloc, false);
            }

            timer.stop("PaToH");

            double maxFractionPart0 = (1.0 + epsilon) * numParts[0] / k;
            double maxFractionPart1 = (1.0 + epsilon) * numParts[1] / k;

            timer.start("Refinement", "RecursiveBisector");
            refiner.refine(partition, maxFractionPart0, maxFractionPart1);
            timer.stop("Refinement");

            if (k > 2) {
                std::vector<int> new_ids(partition.size());
                std::vector<int> carries(2, 0);
                std::vector<Partition::PartitionID> vec_part(hg.numNodes());
                for (uint i = 0; i < partition.size(); ++i) {
                    new_ids[i] = carries[partition[i]]++;
                }

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
                            assert(partHg.indexPins().back() + hyperedgeSize == partHg.pins().size());
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
                    partHg.computeVertexIncidences();

                    Partition subPartition(numParts[partID], partHg);

                    partition_recursively(subPartition, epsilon, preset, numParts[partID], false);

                    for (uint i = 0; i < partition.size(); ++i) {
                        if (partition[i] == partID) {
                            vec_part[i] = subPartition[new_ids[i]] + partID * numParts[0];
                        }
                    }

                }
                partition.rebuild(vec_part);
            }
            
			if (alloc) {
				PaToHInterface::freePatoh();
			}
			return;
        }
    };
}