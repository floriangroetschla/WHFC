#pragma once

#include "../recursive_bisection/hypergraph.h"
#include "patoh.h"
#include "../recursive_bisection/partition_ca.h"


class PaToHInterface {
public:

    struct PatohParameters {
        explicit PatohParameters(double epsilon = 0.0) : use_target_weights(false), target_weights(), epsilon(epsilon),
                                                         k(2), alloc(true), free(true) {}

        bool use_target_weights;
        std::vector<float> target_weights;
        double epsilon;
        int k;
        bool alloc;
        bool free;
    };

    template<class PartitionImpl>
    static void bisectImbalancedWithPatoh(PartitionImpl &partition,
                                          int seed,
                                          float imbalanceFactor,
                                          double epsilon = 0.05,
                                          std::string preset = "D", bool alloc = true, bool free = true) {
        PatohParameters p;
        p.use_target_weights = true;
        p.k = 2;
        p.epsilon = epsilon;
        p.target_weights = {float(1),
                            float(imbalanceFactor)};    // Note(Lars): I think these might be just upper bounds on the part weights.
        p.alloc = alloc;
        p.free = free;
        runPatoh(partition, seed, p, preset);
    }

    template<class PartitionImpl>
    static void bisectWithPatoh(PartitionImpl &partition,
                                int seed,
                                double epsilon = 0.0,
                                std::string preset = "D", bool alloc = true, bool free = true) {
        PatohParameters p(epsilon);
        p.alloc = alloc;
        p.free = free;

        runPatoh(partition, seed, p, preset);
    }

    template<class PartitionImpl>
    static void partitionWithPatoh(PartitionImpl &partition, int seed, int numPartitions, double epsilon = 0.05,
                                   std::string preset = "D") {
        PatohParameters p;
        p.use_target_weights = false;
        p.k = numPartitions;
        p.epsilon = epsilon;
        runPatoh(partition, seed, p, preset);
    }

    template<class PartitionImpl>
    static void runPatoh(PartitionImpl &partition, int seed, PatohParameters params, std::string str_preset = "D") {
        whfc_rb::CSRHypergraph &hg = partition.getGraph();
        // For output of PaToH
        //std::vector<whfc_rb::Partition::PartitionID > vec_partition(hg.numNodes());
        std::vector<int> vec_partweights(params.k, 0);

        PaToH_Parameters args;
        int preset = -1;
        if (str_preset == "D") { preset = PATOH_SUGPARAM_DEFAULT; }
        else if (str_preset == "Q") { preset = PATOH_SUGPARAM_QUALITY; }
        else if (str_preset == "S") { preset = PATOH_SUGPARAM_SPEED; }
        else { throw std::runtime_error("Unknown PaToH preset" + str_preset); }

        PaToH_Initialize_Parameters(&args, PATOH_CONPART, preset);

        args._k = params.k;
        args.doinitperm = 0;
        args.final_imbal = params.epsilon;

        if (str_preset == "Q") {
            args.MemMul_Pins = 100;
            args.MemMul_CellNet = 100;
            args.MemMul_General = 100;
        }

        int c, n, nconst, *cwghts, *nwghts, *xpins, *pins, *partvec, cut, *partweights;
        float *targetweights;
        n = hg.numHyperedges();        // Note(Lars): Use static_cast<desired_type>( value )
        c = hg.numNodes();
        nconst = 1;
        partvec = reinterpret_cast<int *>(partition.data());
        partweights = vec_partweights.data();
        cwghts = reinterpret_cast<int *>(hg.nodeWeights().data());
        nwghts = reinterpret_cast<int *>(hg.hyperedgeWeights().data());
        pins = reinterpret_cast<int *>(hg.pins().data());
        xpins = reinterpret_cast<int *>(hg.indexPins().data());
        if (params.use_target_weights) {
            targetweights = params.target_weights.data();
        } else {
            targetweights = nullptr;
        }

        if (params.alloc) {
            PaToH_Alloc(&args, c, n, nconst, cwghts, nwghts, xpins, pins);
        }
        args.seed = seed;
        PaToH_Part(&args, c, n, nconst, 0, cwghts, nwghts, xpins, pins, targetweights, partvec, partweights, &cut);

        if (params.free) {
            PaToH_Free();
        }
        partition.initialize();     // TODO this should be measured separately from PaToH time
    }

    static int freePatoh() {
        return PaToH_Free();
    }
};
