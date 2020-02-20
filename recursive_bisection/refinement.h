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

        bool refine(Partition& partition, CSRHypergraph& hg, double epsilon) {
            std::vector<CSRHypergraph::NodeWeight> partWeights = partition.partitionWeights(hg);
            CSRHypergraph::NodeWeight totalWeight = partWeights[0] + partWeights[1];

            /*
            double alpha = 16.0;
            double maxW0 = (1 + alpha * epsilon) * std::ceil(static_cast<double>(totalWeight) / 2.0) - partWeights[0];
            double maxW1 = (1 + alpha * epsilon) * std::ceil(static_cast<double>(totalWeight) / 2.0) - partWeights[1];
            */

            double maxW0 = 0.2 * partWeights[0];
            double maxW1 = 0.2 * partWeights[1];

            maxW0 = std::min(maxW0, 0.999 * partWeights[0]);
            maxW1 = std::min(maxW1, 0.999 * partWeights[1]);



            whfc::NodeWeight maxBlockWeight = std::ceil((1.0 + epsilon) * totalWeight / 2.0);

            FlowHypergraphBuilderExtractor::ExtractorInfo extractor_info = extractor.run(hg, partition, 0, 1, maxW0, maxW1);

            std::cout << "Max block weight: " << maxBlockWeight << std::endl;

            // call WHFC to improve the bisection
            hfc.reset();
            hfc.cs.setMaxBlockWeight(0, maxBlockWeight);
            hfc.cs.setMaxBlockWeight(1, maxBlockWeight);
            hfc.upperFlowBound = extractor_info.cutAtStake - extractor_info.baseCut;

            if (extractor_info.cutAtStake - extractor_info.baseCut == 0) return false;

            std::cout << "WHFC started" << std::endl;
            bool hfc_result = hfc.runUntilBalancedOrFlowBoundExceeded(extractor_info.source, extractor_info.target);
            std::cout << "WHFC finished" << std::endl;

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