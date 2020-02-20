#pragma once

#include "fhgb_extraction.h"
#include "../algorithm/hyperflowcutter.h"
#include "../algorithm/dinic.h"
#include <random>

namespace whfc_rb {
    class WHFCRefiner {
    public:
        using PartitionID = int;
        using NodeID = CSRHypergraph::NodeID;
        using HyperedgeID = CSRHypergraph::HyperedgeID;

        WHFCRefiner(uint maxNumNodes, uint maxNumEdges, uint maxNumPins, std::mt19937& mt) :
            extractor(maxNumNodes, maxNumEdges, maxNumPins, mt),
            hfc(extractor.fhgb, mt()), mt(mt) {

        }

        bool refine(Partition& partition, CSRHypergraph& hg, double maxFractionPart0, double maxFractionPart1) {
            std::vector<CSRHypergraph::NodeWeight> partWeights = partition.partitionWeights(hg);
            CSRHypergraph::NodeWeight totalWeight = partWeights[0] + partWeights[1];

            double maxW0 = 0.2 * partWeights[0];
            double maxW1 = 0.2 * partWeights[1];

            whfc::NodeWeight maxBlockWeight0 = std::ceil(maxFractionPart0 * totalWeight);
            whfc::NodeWeight maxBlockWeight1 = std::ceil(maxFractionPart1 * totalWeight);

            FlowHypergraphBuilderExtractor::ExtractorInfo extractor_info = extractor.run(hg, partition, 0, 1, maxW0, maxW1);


            // call WHFC to improve the bisection
            hfc.reset();
            hfc.cs.setMaxBlockWeight(0, maxBlockWeight0);
            hfc.cs.setMaxBlockWeight(1, maxBlockWeight1);
            hfc.upperFlowBound = extractor_info.cutAtStake - extractor_info.baseCut;

            if (extractor_info.cutAtStake - extractor_info.baseCut == 0) return false;

            bool hfc_result = hfc.runUntilBalancedOrFlowBoundExceeded(extractor_info.source, extractor_info.target);

            if (!hfc_result) return false;

            whfc::Flow newCut = extractor_info.baseCut + hfc.cs.flowValue;

            if (newCut < extractor_info.cutAtStake) {
                reassign(partition, hg, extractor_info);
                return true;
            } else if (newCut == extractor_info.cutAtStake) {

            }

            return false;
        }

    private:
        FlowHypergraphBuilderExtractor extractor;
        whfc::HyperFlowCutter<whfc::Dinic> hfc;
        std::mt19937 mt;

        void reassign(Partition& partition, CSRHypergraph& hg, FlowHypergraphBuilderExtractor::ExtractorInfo& info) {
            for (whfc::Node localID : extractor.localNodeIDs()) {
                if (localID == info.source || localID == info.target) continue;
                NodeID globalID = extractor.local2global(localID);
                PartitionID newPart = hfc.cs.n.isSource(localID) ? 0 : 1;
                partition[globalID] = newPart;
            }
        }
    };
}