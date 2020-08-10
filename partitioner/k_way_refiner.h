#pragma once

#include <boost/dynamic_bitset.hpp>


namespace whfc_rb {
    template<class PartitionImpl, class HypergraphImpl, class FlowAlgo, class Extractor>
    class KWayRefiner {
    public:
        explicit KWayRefiner(PartitionImpl &partition, whfc::TimeReporter& timer, std::mt19937 &mt, const PartitionerConfig& config) : partition(
                partition), partActive(partition.numParts()), timer(timer), mt(mt), twoWayRefiner(
                partition.getGraph().numNodes(), partition.getGraph().numHyperedges(), partition.getGraph().numPins(),
                mt(), config), config(config) {}

        uint refine(double epsilon, uint maxIterations) {
            NodeWeight maxWeight = (1.0 + epsilon) * partition.totalWeight() / static_cast<double>(partition.numParts());
            partActive.set();
            uint iterations = 0;
            std::vector<std::pair<PartitionBase::PartitionID, PartitionBase::PartitionID>> partitionPairs;

            while (partActive.count() > 0 && iterations < maxIterations) {
                fillPartitionPairs(partitionPairs);
                partActive.reset();
                std::shuffle(partitionPairs.begin(), partitionPairs.end(), mt);

                for (auto partitionPair : partitionPairs) {
                    if (iterations >= maxIterations) break;
                    PartitionBase::PartitionID part0 = partitionPair.first;
                    PartitionBase::PartitionID part1 = partitionPair.second;
                    timer.start("WHFCRefinerTwoWay");
                    bool refinementResult = twoWayRefiner.refine(partition, part0, part1, maxWeight, maxWeight, timer);
                    timer.stop("WHFCRefinerTwoWay");
                    if (refinementResult) {
                        // Schedule for next round
                        partActive.set(part0);
                        partActive.set(part1);
                    }
                    iterations++;
                }
            }
            return iterations;
        }

    private:
        PartitionImpl &partition;
        boost::dynamic_bitset<> partActive;
        whfc::TimeReporter &timer;
        std::mt19937 &mt;
        WHFCRefinerTwoWay<PartitionImpl, HypergraphImpl, FlowAlgo, Extractor> twoWayRefiner;
        const PartitionerConfig& config;

        void fillPartitionPairs(std::vector<std::pair<PartitionBase::PartitionID, PartitionBase::PartitionID>>& partitionPairs) {
            partitionPairs.clear();
            for (PartitionBase::PartitionID part0 = 0; part0 < partActive.size() - 1; ++part0) {
                if (partActive[part0]) {
                    for (PartitionBase::PartitionID part1 = part0 + 1; part1 < partActive.size(); ++part1) {
                        partitionPairs.push_back(std::pair<PartitionBase::PartitionID, PartitionBase::PartitionID>(part0, part1));
                    }
                }
            }
        }
    };
}