#pragma once

#include <boost/dynamic_bitset.hpp>


namespace whfc_rb {
    class KWayRefiner {
    public:
        explicit KWayRefiner(PartitionCA &partition, whfc::TimeReporter &timer, std::mt19937 &mt) : partition(
                partition), partActive(partition.numParts()), partActiveNextRound(partition.numParts()), timer(timer), mt(mt), twoWayRefiner(
                partition.getGraph().numNodes(), partition.getGraph().numHyperedges(), partition.getGraph().numPins(),
                mt, timer) {}

        void refine(double epsilon) {
            NodeWeight maxWeight = (1.0 + epsilon) * partition.totalWeight() / static_cast<double>(partition.numParts());
            partActive.set();
            partActiveNextRound.reset();
            while (partActive.count() > 0) {
                for (uint part0 = 0; part0 < partActive.size() - 1; ++part0) {
                    if (partActive[part0]) {
                        for (uint part1 = part0 + 1; part1 < partActive.size(); ++part1) {
                            bool refinementResult = twoWayRefiner.refine(partition, part0, part1, maxWeight, maxWeight, true);
                            std::cout << "Refinement of part " << part0 << " and " << part1 << " resulted in "
                                      << refinementResult << std::endl;
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