#pragma once

#include <boost/dynamic_bitset.hpp>


namespace whfc_rb {
    class KWayRefiner {
    public:
        explicit KWayRefiner(PartitionCA &partition, whfc::TimeReporter &timer, std::mt19937 &mt) : partition(
                partition), partActive(partition.numParts()), partActiveNextRound(partition.numParts()), timer(timer), mt(mt), twoWayRefiner(
                partition.getGraph().numNodes(), partition.getGraph().numHyperedges(), partition.getGraph().numPins(),
                mt, timer) {}

        void refine(double epsilon, uint maxIterations) {
            NodeWeight maxWeight = (1.0 + epsilon) * partition.totalWeight() / static_cast<double>(partition.numParts());
            partActive.set();
            partActiveNextRound.reset();
            uint iterations = 0;

            // Note(Lars): support random order? the maintenance code is probably not too running time intensive
            while (partActive.count() > 0 && iterations < maxIterations) {
                for (uint part0 = 0; part0 < partActive.size() - 1; ++part0) {
                    if (partActive[part0]) {
                        for (uint part1 = part0 + 1; part1 < partActive.size(); ++part1) {
                            iterations++;
                            bool refinementResult = twoWayRefiner.refine(partition, part0, part1, maxWeight, maxWeight);
                            if (refinementResult) {
                                // Schedule for next round
                                partActiveNextRound.set(part0);
                                partActiveNextRound.set(part1);
                            }
                        }
                    }
                }
                std::swap(partActive, partActiveNextRound);
                partActiveNextRound.reset();
            }
        }

    private:
        PartitionCA &partition;
        boost::dynamic_bitset<> partActive;
        boost::dynamic_bitset<> partActiveNextRound;
        whfc::TimeReporter &timer;
        std::mt19937 &mt;
        WHFCRefinerTwoWay twoWayRefiner;
    };
}