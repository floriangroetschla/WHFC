#pragma once

namespace whfc_rb {

    struct PartitionConfig {
        bool refine;
        std::string patoh_preset;
        bool precomputeCuts;
        bool distancePiercing;
    };
}
