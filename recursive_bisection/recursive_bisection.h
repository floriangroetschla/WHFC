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
            mt(mt), timer(timer) {

        }

        Partition run(CSRHypergraph& hg, double epsilon, std::string preset, uint k) {
            epsilon = std::pow(1.0 + epsilon, 1.0 / std::ceil(std::log2(k))) - 1.0;

            return partition_recursively(hg, epsilon, preset, k, true);
        }

    private:
        std::mt19937& mt;
        whfc::TimeReporter& timer;

        Partition partition_recursively(CSRHypergraph& hg, double epsilon, std::string preset, uint k, bool alloc) {
            if (k == 1) {
                return Partition(hg.numNodes(), 1);
            }

            std::array<int, 2> numParts;
            Partition partition;

            timer.start("PaToH", "RecursiveBisector");
            if (k % 2 == 0) {
                numParts[0] = k / 2;
                numParts[1] = k / 2;
                partition = PaToHInterface::bisectWithPatoh(hg, mt(), epsilon, preset, alloc, false);
            } else {
                numParts[0] = k / 2;
                numParts[1] = numParts[0] + 1;
                partition = PaToHInterface::bisectImbalancedWithPatoh(hg, mt(), float(numParts[1]) / float(numParts[0]), epsilon, preset, alloc, false);
            }
            timer.stop("PaToH");

            double maxFractionPart0 = (1.0 + epsilon) * numParts[0] / k;
            double maxFractionPart1 = (1.0 + epsilon) * numParts[1] / k;

            if (k > 2) {
                partition.setNumParts(k);
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

                    sub_partitions[partID] = partition_recursively(partHg, epsilon, preset, numParts[partID], false);
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