#pragma once

namespace whfc_rb {
    class KWayRefiner {
    public:
        explicit KWayRefiner(PartitionCA& partition, whfc::TimeReporter& timer, std::mt19937& mt) : partition(partition), partActive(partition.numParts(), true)
        , timer(timer), mt(mt), twoWayRefiner(partition.getGraph().numNodes(), partition.getGraph().numHyperedges(), partition.getGraph().numPins(), mt, timer) {}

        void refine(double epsilon) {
            while (true) {

                uint block0;
                uint block1;
                bool foundPair = false;
                for (block0 = 0; block0 < partActive.size() - 1; ++block0) {
                    if (partActive[block0]) {
                        for (block1 = block0 + 1; block1 < partActive.size(); ++block1) {
                            if (partActive[block1]) {
                                foundPair = true;
                                break;
                            }
                        }
                        break;
                    }
                }

                if (foundPair) {
                    NodeWeight maxWeight = (1.0 + epsilon) * partition.totalWeight() / static_cast<double>(partition.numParts());

                    bool refinementResult = twoWayRefiner.refine(partition, block0, block1, maxWeight, maxWeight);
                    if (!refinementResult) {
                        partActive[block0] = false;
                        partActive[block1] = false;
                    }
                } else {
                    break;
                }
            }
        }

    private:
        PartitionCA& partition;
        std::vector<bool> partActive;
        whfc::TimeReporter& timer;
        std::mt19937 &mt;
        WHFCRefinerTwoWay twoWayRefiner;
    };
}