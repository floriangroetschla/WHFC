#pragma once

namespace whfc_rb {
    class TwoWayRefinerInterface {
    public:
        TwoWayRefinerInterface(const TwoWayRefinerInterface&) = delete;
        TwoWayRefinerInterface(TwoWayRefinerInterface&&) = delete;
        TwoWayRefinerInterface& operator= (const TwoWayRefinerInterface&) = delete;
        TwoWayRefinerInterface& operator= (TwoWayRefinerInterface&&) = delete;

        TwoWayRefinerInterface() = default;
        virtual ~TwoWayRefinerInterface() = default;

        virtual bool refine(Partition& partition, Partition::PartitionID part0, Partition::PartitionID part1, NodeWeight maxBlockWeight0, NodeWeight maxBlockWeight1) = 0;
    };
}
