#pragma once

namespace whfc_rb {
    template<class PartitionImpl>
    class NullRefiner : public TwoWayRefinerInterface {
    public:
        NullRefiner(uint maxNumNodes, uint maxNumEdges, uint maxNumPins, int seed, const PartitionerConfig& config) {}

        bool refine(PartitionImpl &partition, PartitionBase::PartitionID part0, PartitionBase::PartitionID part1, NodeWeight maxBlockWeight0,
                    NodeWeight maxBlockWeight1, whfc::TimeReporter& timer) {
            return false;
        }
    };
}
