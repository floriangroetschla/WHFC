#pragma once

#include "fhgb_extraction.h"
#include "../algorithm/hyperflowcutter.h"
#include "../algorithm/dinic.h"
#include "../io/whfc_io.h"
#include <random>
#include "two_way_refiner_interface.h"

namespace whfc_rb {
    class WHFCRefinerTwoWay : public TwoWayRefinerInterface {
    public:
        using PartitionID = PartitionBase::PartitionID;

        WHFCRefinerTwoWay(uint maxNumNodes, uint maxNumEdges, uint maxNumPins, std::mt19937 &mt,
                          whfc::TimeReporter &timer) :
                extractor(maxNumNodes, maxNumEdges, maxNumPins, mt),
                hfc(extractor.fhgb, mt()), mt(mt), timer(timer) {
            hfc.timer.active = false;
        }

        template<class PartitionImpl>
        bool refine(PartitionImpl &partition, PartitionID part0, PartitionID part1, NodeWeight maxBlockWeight0,
                    NodeWeight maxBlockWeight1, bool globalOptimization) {
            std::vector<NodeWeight> partWeights = partition.partitionWeights();

            double maxW0 = 0.2 * partWeights[part0];
            double maxW1 = 0.2 * partWeights[part1];

            double imbalanceBefore = std::max(partWeights[part0] / maxBlockWeight0,
                                              partWeights[part1] / maxBlockWeight1);

            timer.start("Extraction", "Refinement");
            FlowHypergraphBuilderExtractor::ExtractorInfo extractor_info = extractor.run(partition, part0, part1, maxW0,
                                                                                         maxW1);
            timer.stop("Extraction");

            // call WHFC to improve the bisection
            hfc.reset();
            hfc.cs.setMaxBlockWeight(0, maxBlockWeight0);
            hfc.cs.setMaxBlockWeight(1, maxBlockWeight1);
            hfc.upperFlowBound = extractor_info.cutAtStake - extractor_info.baseCut;

            if (extractor_info.cutAtStake - extractor_info.baseCut == 0) return false;

            static constexpr bool write_snapshot = false;
            if constexpr (write_snapshot) {
                writeSnapshot(extractor_info);
            }

            timer.start("WHFC", "Refinement");
            bool hfc_result = hfc.runUntilBalancedOrFlowBoundExceeded(extractor_info.source, extractor_info.target);
            timer.stop("WHFC");

            if (!hfc_result) return false;

            whfc::Flow newCut = extractor_info.baseCut + hfc.cs.flowValue;

            double imbalanceAfter = std::max(hfc.cs.n.sourceReachableWeight / maxBlockWeight0,
                                             hfc.cs.n.targetReachableWeight / maxBlockWeight1);

            assert(hfc.cs.n.sourceReachableWeight <= maxBlockWeight0 &&
                   hfc.cs.n.targetReachableWeight <= maxBlockWeight1);

            if (globalOptimization) {
                NodeWeight km1Before = partition.km1Objective();
                std::vector<PartitionBase::PartitionChangeElement> partChanges;
                for (whfc::Node localID : extractor.localNodeIDs()) {
                    if (localID == extractor_info.source || localID == extractor_info.target) continue;
                    NodeID globalID = extractor.local2global(localID);
                    PartitionID newPart = hfc.cs.n.isSource(localID) ? part0 : part1;
                    if (newPart != partition[globalID]) {
                        partChanges.push_back({globalID, newPart});
                    }
                }
                NodeWeight km1After = partition.km1AfterChanges(partChanges);
                if (km1After < km1Before || (km1After == km1Before && imbalanceAfter < imbalanceBefore)) {
                    // TODO: factor out this function call
                    reassign(partition, extractor_info, part0, part1);
                }
            } else {
                if (newCut < extractor_info.cutAtStake ||
                    (newCut == extractor_info.cutAtStake && imbalanceAfter < imbalanceBefore)) {
                    // TODO: this one too
                    reassign(partition, extractor_info, part0, part1);
                    return true;
                }
            }


            return false;
        }

    private:
        FlowHypergraphBuilderExtractor extractor;
        whfc::HyperFlowCutter<whfc::Dinic> hfc;
        std::mt19937 &mt;
        whfc::TimeReporter &timer;

        size_t instance_counter = 0;

        template<class PartitionImpl>
        void reassign(PartitionImpl &partition, FlowHypergraphBuilderExtractor::ExtractorInfo &info, PartitionID part0,
                      PartitionID part1) {
            for (whfc::Node localID : extractor.localNodeIDs()) {
                if (localID == info.source || localID == info.target) continue;
                NodeID globalID = extractor.local2global(localID);
                PartitionID newPart = hfc.cs.n.isSource(localID) ? part0 : part1;
                partition.changePart(globalID, newPart);
            }
        }

        void writeSnapshot(FlowHypergraphBuilderExtractor::ExtractorInfo &extractor_info) {
            whfc::WHFC_IO::WHFCInformation i = {
                    {hfc.cs.maxBlockWeight(0), hfc.cs.maxBlockWeight(1)},
                    extractor_info.cutAtStake - extractor_info.baseCut,
                    extractor_info.source,
                    extractor_info.target
            };
            std::string hg_filename = "Snapshot" + std::to_string(instance_counter);
            instance_counter++;
            std::cout << "Wrote snapshot: " << hg_filename << std::endl;
            whfc::HMetisIO::writeFlowHypergraph(extractor.fhgb, hg_filename);
            whfc::WHFC_IO::writeAdditionalInformation(hg_filename, i, hfc.cs.rng);
        }

    };
}