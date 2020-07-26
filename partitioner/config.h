#pragma once

namespace whfc_rb {

    struct PartitionerConfig {
        bool refine;
        std::string patoh_preset;
        bool precomputeCuts;
        bool distancePiercing;
        size_t numThreads;
        std::string graphName;
        size_t k;
        double percentage_bfs_from_cut = 0.2;
    };
}
