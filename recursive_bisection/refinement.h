#pragma once

#include "fhgb_extraction.h"
#include "../algorithm/hyperflowcutter.h"
#include "../algorithm/dinic.h"
#include "../io/whfc_io.h"
#include <random>

namespace whfc_rb {
    class WHFCRefiner {
    public:
        using PartitionID = int;

        WHFCRefiner(uint maxNumNodes, uint maxNumEdges, uint maxNumPins, std::mt19937& mt, whfc::TimeReporter& timer) :
            extractor(maxNumNodes, maxNumEdges, maxNumPins, mt),
            hfc(extractor.fhgb, mt()), mt(mt), timer(timer) {
                hfc.timer.active = false;
        }

        bool refine(Partition& partition, double maxFractionPart0, double maxFractionPart1) {
            std::vector<NodeWeight> partWeights = partition.partitionWeights();
            NodeWeight totalWeight = partWeights[0] + partWeights[1];

            double maxW0 = 0.2 * partWeights[0];
            double maxW1 = 0.2 * partWeights[1];

            whfc::NodeWeight maxBlockWeight0 = std::ceil(maxFractionPart0 * totalWeight);
            whfc::NodeWeight maxBlockWeight1 = std::ceil(maxFractionPart1 * totalWeight);

            double imbalanceBefore = std::max(partWeights[0] / (maxFractionPart0 * totalWeight), partWeights[1] / (maxFractionPart1 + totalWeight));

            timer.start("Extraction", "Refinement");
            FlowHypergraphBuilderExtractor::ExtractorInfo extractor_info = extractor.run(partition, 0, 1, maxW0, maxW1);
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

            double imbalanceAfter = std::max(hfc.cs.n.sourceReachableWeight / (maxFractionPart0 * totalWeight), hfc.cs.n.targetReachableWeight / (maxFractionPart1 + totalWeight));

            if (newCut < extractor_info.cutAtStake || (newCut == extractor_info.cutAtStake && imbalanceAfter < imbalanceBefore)) {
                reassign(partition, extractor_info);
                return true;
            }

            return false;
        }

    private:
        FlowHypergraphBuilderExtractor extractor;
        whfc::HyperFlowCutter<whfc::Dinic> hfc;
        std::mt19937& mt;
        whfc::TimeReporter& timer;

        size_t instance_counter = 0;

        void reassign(Partition& partition, FlowHypergraphBuilderExtractor::ExtractorInfo& info) {
            for (whfc::Node localID : extractor.localNodeIDs()) {
                if (localID == info.source || localID == info.target) continue;
                NodeID globalID = extractor.local2global(localID);
                PartitionID newPart = hfc.cs.n.isSource(localID) ? 0 : 1;
                partition.changePart(globalID, newPart);
            }
        }

        void writeSnapshot(FlowHypergraphBuilderExtractor::ExtractorInfo& extractor_info) {
            whfc::WHFC_IO::WHFCInformation i = {
                    { hfc.cs.maxBlockWeight(0), hfc.cs.maxBlockWeight(1) },
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