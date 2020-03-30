#include "recursive_bisection/hypergraph.h"
#include "io/hmetis_io.h"
#include <iostream>
#include "extern/patoh_wrapper.h"
#include <random>
#include "util/timer.h"
#include <tbb/tbb.h>
#include "recursive_bisection/k_way_refiner_parallel.h"
#include "recursive_bisection/partition_threadsafe.h"
#include "recursive_bisection/tbb_thread_pinning.h"
#include "recursive_bisection/config.h"

void printStatistics(whfc_rb::PartitionBase &partition, whfc::TimeReporter &timer) {
    timer.report(std::cout);
    std::cout << "Imbalance: " << partition.imbalance() << std::endl;
    std::cout << "Num_parts: " << partition.numParts() << std::endl;
    std::cout << "km1: " << partition.km1Objective() << std::endl;
    std::cout << "cut: " << partition.cutObjective() << std::endl;
}

int main(int argc, const char *argv[]) {

    if (argc != 9) {
        throw std::runtime_error("Usage ./KWayRefinementParallel HypergraphFile epsilon k seed preset numThreads useThreadPinning(0 or 1) distancePiercing(0 or 1)");
    }
    whfc_rb::CSRHypergraph hg = whfc::HMetisIO::readCSRHypergraph(argv[1]);
    double epsilon = std::stod(argv[2]);
    uint numParts = std::stoul(argv[3]);
    int seed = std::stoi(argv[4]);
    std::string patoh_preset = argv[5];
    uint numThreads = std::stoi(argv[6]);
    bool useThreadPinning = std::stoi(argv[7]);
    bool distancePiercing = std::stoi(argv[8]);
    std::mt19937 mt(seed);

    //uint maxIterations = 20 * numParts * numParts;std::numeric_limits<uint>::max();
    uint maxIterations = std::numeric_limits<uint>::max();

    tbb::task_scheduler_init init(numThreads);
    whfc_rb::pinning_observer pinner;
    pinner.observe(useThreadPinning);

    whfc::TimeReporter timer("Total");

    bool precomputeCuts = true;

    whfc_rb::PartitionConfig config = {true, patoh_preset, precomputeCuts, distancePiercing};
    std::cout << "useThreadPinning: " << useThreadPinning << std::endl;
    std::cout << "precomputeCuts: " << precomputeCuts << std::endl;
    std::cout << "distancePiercing: " << distancePiercing << std::endl;
    std::cout << "numThreads: " << numThreads << std::endl;

    timer.start("Total");
    timer.start("PaToH", "Total");
    whfc_rb::PartitionThreadsafe partition(numParts, hg);
    PaToHInterface::partitionWithPatoh(partition, seed, numParts, epsilon, patoh_preset);
    timer.stop("PaToH");

    printStatistics(partition, timer);

    timer.start("Refinement", "Total");
    whfc_rb::KWayRefinerParallel refiner(partition, timer, mt, config);
    uint iterations = refiner.refine(epsilon, maxIterations, seed);
    timer.stop("Refinement");
    timer.stop("Total");

    printStatistics(partition, timer);
    std::cout << "numThreads: " << numThreads << std::endl;
    std::cout << "maxIterations: " << maxIterations << std::endl;
    std::cout << "iterations_done: " << iterations << std::endl;
}

