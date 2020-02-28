#pragma once

namespace whfc_rb {
    class NullRefiner : public TwoWayRefinerInterface {
    public:
        NullRefiner(uint maxNumNodes, uint maxNumEdges, uint maxNumPins, std::mt19937 &mt, whfc::TimeReporter &timer) {}

        bool refine(PartitionBase &partition, PartitionBase::PartitionID part0, PartitionBase::PartitionID part1,
                    NodeWeight maxBlockWeight0, NodeWeight maxBlockWeight1) {
            return false;
        }
    };
}
