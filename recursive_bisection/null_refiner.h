#pragma once

namespace whfc_rb {
    class NullRefiner : public TwoWayRefinerInterface {
    public:
        NullRefiner(uint maxNumNodes, uint maxNumEdges, uint maxNumPins, std::mt19937 &mt, whfc::TimeReporter &timer) {}

        bool refine(Partition &partition, Partition::PartitionID part0, Partition::PartitionID part1,
                    NodeWeight maxBlockWeight0, NodeWeight maxBlockWeight1) {
            return false;
        }
    };
}
