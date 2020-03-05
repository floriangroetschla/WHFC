#include <stdexcept>
#include <cmath>
#include "../extern/patoh_wrapper.h"
#include "../io/hmetis_io.h"
#include "hypergraph.h"
#include "whfc_refiner_two_way.h"
#include <random>
#include "config.h"

namespace whfc_rb {
    class RecursiveBisector {
    public:

        RecursiveBisector(uint maxNumNodes, uint maxNumEdges, uint maxNumPins, std::mt19937 &mt,
                          whfc::TimeReporter &timer, PartitionConfig &config) :
                refiner(maxNumNodes, maxNumEdges, maxNumPins, mt, timer), mt(mt), timer(timer), config(config) {

        }

        PartitionBase run(CSRHypergraph &hg, double epsilon, uint k) {
            epsilon = std::pow(1.0 + epsilon, 1.0 / std::ceil(std::log2(k))) - 1.0;

            PartitionBase partition(k, hg);
            partition_recursively(partition, epsilon, k, true);
            partition.initialize();
            return partition;
        }

    private:
        WHFCRefinerTwoWay refiner;
        std::mt19937 &mt;
        whfc::TimeReporter &timer;
        PartitionConfig &config;

        void partition_recursively(PartitionBase &partition, double epsilon, uint k, bool alloc) {
            // insert assertions here
            if (k == 1) {
                return;
            }

            std::array<int, 2> numParts;
            CSRHypergraph &hg = partition.getGraph();

            timer.start("PaToH", "Total");
            if (k % 2 == 0) {
                numParts[0] = k / 2;
                numParts[1] = k / 2;
                PaToHInterface::bisectWithPatoh(partition, mt(), epsilon, config.patoh_preset, alloc, false);
            } else {
                numParts[0] = k / 2;
                numParts[1] = numParts[0] + 1;
                PaToHInterface::bisectImbalancedWithPatoh(partition, mt(), float(numParts[1]) / float(numParts[0]),
                                                          epsilon, config.patoh_preset, alloc, false);
            }

            timer.stop("PaToH");

            NodeWeight maxWeight0 = (1.0 + epsilon) * static_cast<double>(numParts[0]) /
                                    static_cast<double>(k) * partition.totalWeight();
            NodeWeight maxWeight1 = (1.0 + epsilon) * static_cast<double>(numParts[1]) /
                                    static_cast<double>(k) * partition.totalWeight();

            if (config.refine) {
                timer.start("Refinement", "Total");
                refiner.refine(partition, 0, 1, maxWeight0, maxWeight1, false);
                timer.stop("Refinement");
            }

            if (k > 2) {
                timer.start("GraphAndPartitionBuilding", "Total");
                std::vector<int> new_ids(partition.size());
                std::vector<int> carries(2, 0);
                std::vector<PartitionBase::PartitionID> vec_part(hg.numNodes());
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

                    if (config.refine) {
                        partHg.computeVertexIncidences();
                    }

                    PartitionBase subPartition(numParts[partID], partHg);
                    timer.stop("GraphAndPartitionBuilding");
                    partition_recursively(subPartition, epsilon, numParts[partID], false);
                    timer.start("GraphAndPartitionBuilding", "Total");

                    for (uint i = 0; i < partition.size(); ++i) {
                        if (partition[i] == partID) {
                            vec_part[i] = subPartition[new_ids[i]] + partID * numParts[0];
                        }
                    }

                }
                partition.replace(vec_part);
                timer.stop("GraphAndPartitionBuilding");
            }

            if (alloc) {
                PaToHInterface::freePatoh();
            }
            return;
        }
    };
}