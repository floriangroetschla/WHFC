#pragma once

namespace whfc_rb {
    class TwoWayRefinerInterface {
    public:
        TwoWayRefinerInterface(const TwoWayRefinerInterface &) = delete;

        TwoWayRefinerInterface(TwoWayRefinerInterface &&) = delete;

        TwoWayRefinerInterface &operator=(const TwoWayRefinerInterface &) = delete;

        TwoWayRefinerInterface &operator=(TwoWayRefinerInterface &&) = delete;

        TwoWayRefinerInterface() = default;

        virtual ~TwoWayRefinerInterface() = default;

        template<class PartitionImpl>
         bool refine(PartitionImpl &partition, PartitionBase::PartitionID part0, PartitionBase::PartitionID part1,
               NodeWeight maxBlockWeight0, NodeWeight maxBlockWeight1);
    };
}
