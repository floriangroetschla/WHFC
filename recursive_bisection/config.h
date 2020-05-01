#pragma once

namespace whfc_rb {

    struct PartitionConfig {
        bool refine;
        std::string patoh_preset;
        bool precomputeCuts;
        bool distancePiercing;
        size_t numThreads;
        double percentage_bfs_from_cut = 0.2;
    };
}
