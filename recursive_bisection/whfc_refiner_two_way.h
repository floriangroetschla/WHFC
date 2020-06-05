#pragma once

#include "fhgb_extraction.h"
#include "../algorithm/hyperflowcutter.h"
#include "../algorithm/dinic.h"
#include "../io/whfc_io.h"
#include <random>
#include "two_way_refiner_interface.h"
#include "config.h"
#include "../algorithm/push_relabel.h"

namespace whfc_rb {
    class WHFCRefinerTwoWay : public TwoWayRefinerInterface {
    public:
        using PartitionID = PartitionBase::PartitionID;

        WHFCRefinerTwoWay(uint maxNumNodes, uint maxNumEdges, uint maxNumPins, int seed, const PartitionConfig& config) :
                extractor(maxNumNodes, maxNumEdges, maxNumPins, seed, config),
                hfc(extractor.fhgb, seed, config), timer("WHFCRefinerTwoWay"), config(config)
        {
        }

        template<class PartitionImpl>
        bool refine(PartitionImpl &partition, PartitionID part0, PartitionID part1, NodeWeight maxBlockWeight0,
                    NodeWeight maxBlockWeight1, whfc::TimeReporter& timer) {

            double maxW0 = config.percentage_bfs_from_cut * partition.partWeight(part0);
            double maxW1 = config.percentage_bfs_from_cut * partition.partWeight(part1);

            double imbalanceBefore = std::max(partition.partWeight(part0) / maxBlockWeight0,
                                              partition.partWeight(part1) / maxBlockWeight1);

            timer.start("Extraction", "WHFCRefinerTwoWay");
            FlowHypergraphBuilderExtractor<whfc::LawlerFlowHypergraph>::ExtractorInfo extractor_info = extractor.run(partition, part0, part1, maxW0, maxW1, hfc.cs.borderNodes.distance, timer);
            timer.stop("Extraction");

            // call WHFC to improve the bisection
            hfc.reset();
            hfc.cs.setMaxBlockWeight(0, maxBlockWeight0);
            hfc.cs.setMaxBlockWeight(1, maxBlockWeight1);
            hfc.upperFlowBound = extractor_info.cutAtStake - extractor_info.baseCut;

            if (extractor_info.cutAtStake == extractor_info.baseCut) return false;

            static constexpr bool write_snapshot = false;
            if constexpr (write_snapshot) {
                writeSnapshot(extractor_info);
            }

            timer.start("WHFC", "WHFCRefinerTwoWay");
            bool hfc_result = hfc.enumerateCutsUntilBalancedOrFlowBoundExceeded(extractor_info.source, extractor_info.target);
            timer.stop("WHFC");
            timer.merge(hfc.timer, "WHFC", "HyperFlowCutter");

            if (!hfc_result) return false;

            whfc::Flow newCut = extractor_info.baseCut + hfc.cs.flowValue;

            double imbalanceAfter = std::max(hfc.cs.n.sourceReachableWeight / maxBlockWeight0,
                                             hfc.cs.n.targetReachableWeight / maxBlockWeight1);

            assert(hfc.cs.n.sourceReachableWeight <= maxBlockWeight0 &&
                   hfc.cs.n.targetReachableWeight <= maxBlockWeight1);


            if (newCut < extractor_info.cutAtStake ||
                (newCut == extractor_info.cutAtStake && imbalanceAfter < imbalanceBefore)) {
                timer.start("Reassignment", "WHFCRefinerTwoWay");
                reassign(partition, extractor_info, part0, part1);
                timer.stop("Reassignment");

                return true;
            }

            return false;
        }

        whfc::TimeReporter& getTimer() {
            return timer;
        }

    private:
        FlowHypergraphBuilderExtractor<whfc::LawlerFlowHypergraph> extractor;
        whfc::HyperFlowCutter<whfc::PushRelabel, whfc::LawlerFlowHypergraph> hfc;
        whfc::TimeReporter timer;
        const PartitionConfig& config;

        size_t instance_counter = 0;

        template<class PartitionImpl>
        void reassign(PartitionImpl &partition, FlowHypergraphBuilderExtractor<whfc::LawlerFlowHypergraph>::ExtractorInfo &info, PartitionID part0,
                      PartitionID part1) {
            for (whfc::Node localID : extractor.localNodeIDs()) {
                assert(localID < extractor.fhgb.numNodes());
                if (localID == info.source || localID == info.target) continue;
                NodeID globalID = extractor.local2global(localID);
                assert(globalID < partition.getGraph().numNodes());
                PartitionID newPart = hfc.cs.n.isSource(localID) ? part0 : part1;
                partition.changePart(globalID, newPart);
            }
        }

        void writeSnapshot(FlowHypergraphBuilderExtractor<whfc::LawlerFlowHypergraph>::ExtractorInfo &extractor_info) {
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