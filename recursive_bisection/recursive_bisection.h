#include <stdexcept>
#include <cmath>
#include "../extern/patoh_wrapper.h"
#include "../io/hmetis_io.h"
#include "hypergraph.h"
#include "fhgb_extraction.h"
#include "../algorithm/hyperflowcutter.h"
#include "../algorithm/dinic.h"


namespace whfc_rb {
    class RecursiveBisector {
    public:
        using PartitionID = int;
        static constexpr PartitionID invalidPartition = std::numeric_limits<PartitionID>::max();
        
        using NodeID = CSRHypergraph::NodeID;
        using HyperedgeID = CSRHypergraph::HyperedgeID;

        RecursiveBisector(uint maxNumNodes, uint maxNumEdges, uint maxNumPins, int seed) :
            extractor(maxNumNodes, maxNumEdges, maxNumPins),
            hfc(extractor.fhgb, seed) {

        }

        Partition run(CSRHypergraph& hg, int seed, double epsilon, std::string preset, uint k) {
            epsilon = std::pow(1.0 + epsilon, 1.0 / std::log2(k)) - 1.0;

            return partition_recursively(hg, seed, epsilon, preset, k, true);
        }

    private:
        FlowHypergraphBuilderExtractor extractor;
        whfc::HyperFlowCutter<whfc::Dinic> hfc;

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

            std::vector<CSRHypergraph::NodeWeight> partWeights = partition.partitionWeights(hg);
            CSRHypergraph::NodeWeight totalWeight = partWeights[0] + partWeights[1];

            double alpha = 16.0;

            double maxW0 = (1 + alpha * epsilon) * std::ceil(static_cast<double>(totalWeight) / 2.0) - partWeights[0];
            double maxW1 = (1 + alpha * epsilon) * std::ceil(static_cast<double>(totalWeight) / 2.0) - partWeights[1];

            maxW0 = std::min(maxW0, 0.999 * partWeights[0]);
            maxW1 = std::min(maxW1, 0.999 * partWeights[1]);

            whfc::NodeWeight maxBlockWeight = std::ceil((1.0 + epsilon) * totalWeight / 2.0) + 3;

            partition.print(std::cout);
            FlowHypergraphBuilderExtractor::ExtractorInfo extractor_info = extractor.run(hg, partition, 0, 1, maxW0, maxW1);

            extractor.fhgb.printHypergraph(std::cout);

            // call WHFC to improve the bisection
            hfc.reset();
            hfc.cs.setMaxBlockWeight(0, maxBlockWeight);
            hfc.cs.setMaxBlockWeight(1, maxBlockWeight);
            //hfc.upperFlowBound = 10;

            bool result = hfc.runUntilBalancedOrFlowBoundExceeded(extractor_info.source, extractor_info.target);
            std::cout << "WHFC successful: " << result << std::endl;

            if (k > 2) {
                std::vector<int> new_ids(partition.size());
                std::vector<int> carries(2, 0);
                for (uint i = 0; i < partition.size(); ++i) {
                    new_ids[i] = carries[partition[i]]++;
                }

                std::array<Partition, 2> sub_partitions;

                for (uint partID = 0; partID < 2; ++partID) {
                    CSRHypergraph partHg;
                    for (CSRHypergraph::HyperedgeID e : hg.hyperedges()) {
                        int hyperedgeSize = 0;
                        for (CSRHypergraph::NodeID pin : hg.pinsOf(e)) {
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