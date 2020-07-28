#include "datastructure/hypergraph.h"
#include "io/hmetis_io.h"
#include <iostream>
#include "extern/patoh_wrapper.h"
#include "partitioner/recursive_bisection.h"
#include <random>
#include "util/timer.h"
#include "datastructure/partition_base.h"
#include "partitioner/fhgb_extraction.h"
#include "partitioner/null_refiner.h"
#include "push_relabel/lawler_fhgb_extraction_parallel.h"
#include "partitioner/tbb_thread_pinning.h"
#include <version.h>
#include <tbb/task_scheduler_init.h>


void printStatistics(whfc_rb::PartitionBase &partition, whfc::TimeReporter &timer) {
    std::cout << "commit: " << GIT_COMMIT_HASH << std::endl;
    timer.report(std::cout);
    std::cout << "Imbalance: " << partition.imbalance() << std::endl;
    std::cout << "Num_parts: " << partition.numParts() << std::endl;
    std::cout << "km1: " << partition.km1Objective() << std::endl;
    std::cout << "cut: " << partition.cutObjective() << std::endl;
}

int main(int argc, const char *argv[]) {

    if (argc != 10) {
        throw std::runtime_error("Usage ./KWayRefinement HypergraphFile epsilon k seed preset numThreads useThreadPinning(0 or 1) distancePiercing(0 or 1) mode");
    }
    whfc_rb::CSRHypergraph hg = whfc::HMetisIO::readCSRHypergraph(argv[1]);
    double epsilon = std::stod(argv[2]);
    uint numParts = std::stoul(argv[3]);
    int seed = std::stoi(argv[4]);
    std::string patoh_preset = argv[5];
    uint numThreads = std::stoi(argv[6]);
    bool useThreadPinning = std::stoi(argv[7]);
    bool distancePiercing = std::stoi(argv[8]);
    std::string mode = argv[9];
    std::mt19937 mt(seed);

    tbb::task_scheduler_init init(numThreads);
    whfc_rb::pinning_observer pinner;
    pinner.observe(useThreadPinning);

    bool precomputeCuts = true;

    whfc_rb::PartitionerConfig config = {true, patoh_preset, precomputeCuts, distancePiercing, numThreads, "unset", numParts};

    whfc::TimeReporter timer("Total");

    if (!mode.compare("RBONLY")) {
        std::cout << "Using mode RBONLY" << std::endl;
        config.refine = false;
        whfc_rb::RecursiveBisector recursive_bisector = whfc_rb::RecursiveBisector<whfc_rb::PartitionBase, whfc_rb::NullRefiner<whfc_rb::PartitionBase>>(hg.numNodes(), hg.numHyperedges(), hg.numPins(), mt, timer, config);
        timer.start("Total");
        whfc_rb::PartitionBase partition = recursive_bisector.run(hg, epsilon, numParts);
        timer.stop("Total");
        printStatistics(partition, timer);
    } else if (!mode.compare("RBWHFC")) {
        std::cout << "Using mode RBWHFC" << std::endl;
        config.refine = true;
        whfc_rb::RecursiveBisector recursive_bisector = whfc_rb::RecursiveBisector<whfc_rb::PartitionThreadsafe, whfc_rb::WHFCRefinerTwoWay<whfc_rb::PartitionThreadsafe, whfc_pr::LawlerFlowHypergraphParallel, whfc_pr::PushRelabelParallel, whfc_pr::LawlerFlowHypergraphBuilderExtractorParallel<whfc_pr::LawlerFlowHypergraphParallel, whfc_rb::PartitionThreadsafe>>>(hg.numNodes(), hg.numHyperedges(),
                                                                                                                                                                                                                                                                        hg.numPins(), mt, timer, config);
        timer.start("Total");
        whfc_rb::PartitionBase partition = recursive_bisector.run(hg, epsilon, numParts);
        timer.stop("Total");
        printStatistics(partition, timer);
    } else {
        throw std::runtime_error("Mode must be one of: RBONlY, RBWHFC");
    }
}

