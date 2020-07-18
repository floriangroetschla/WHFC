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
#include "algorithm/dinic.h"
#include "algorithm/dinic_parallel.h"
#include "push_relabel/push_relabel.h"
#include "push_relabel/lawler_fhgb_extraction_parallel.h"
#include "recursive_bisection/fhgb_extraction_parallel.h"
#include "recursive_bisection/fhgb_extraction.h"
#include <filesystem>

void printStatistics(whfc_rb::PartitionBase &partition, whfc::TimeReporter &timer) {
    timer.report(std::cout);
    std::cout << "Imbalance: " << partition.imbalance() << std::endl;
    std::cout << "Num_parts: " << partition.numParts() << std::endl;
    std::cout << "km1: " << partition.km1Objective() << std::endl;
    std::cout << "cut: " << partition.cutObjective() << std::endl;
}

int main(int argc, const char *argv[]) {

    if (argc != 10) {
        throw std::runtime_error("Usage ./KWayRefinementParallel HypergraphFile epsilon k seed preset numThreads useThreadPinning(0 or 1) distancePiercing(0 or 1) maxNumIterations");
    }
    whfc_rb::CSRHypergraph hg = whfc::HMetisIO::readCSRHypergraph(argv[1]);
    double epsilon = std::stod(argv[2]);
    uint numParts = std::stoul(argv[3]);
    int seed = std::stoi(argv[4]);
    std::string patoh_preset = argv[5];
    uint numThreads = std::stoi(argv[6]);
    bool useThreadPinning = std::stoi(argv[7]);
    bool distancePiercing = std::stoi(argv[8]);
    int maxNumIterations = std::stoi(argv[9]);
    std::mt19937 mt(seed);

    tbb::task_scheduler_init init(numThreads);
    whfc_rb::pinning_observer pinner;
    pinner.observe(useThreadPinning);

    whfc::TimeReporter timer("Total");

    bool precomputeCuts = true;

    whfc_rb::PartitionerConfig config = {true, patoh_preset, precomputeCuts, distancePiercing, numThreads, std::filesystem::path(argv[1]).filename(), numParts};
    std::cout << "useThreadPinning: " << useThreadPinning << std::endl;
    std::cout << "precomputeCuts: " << precomputeCuts << std::endl;
    std::cout << "distancePiercing: " << distancePiercing << std::endl;
    std::cout << "numThreads: " << numThreads << std::endl;

    timer.start("Total");
    timer.start("PaToH", "Total");
    whfc_rb::PartitionThreadsafe partition(numParts, hg);
    PaToHInterface::partitionWithPatoh(partition, seed, numParts, epsilon, patoh_preset);
    timer.stop("PaToH");

    //printStatistics(partition, timer);

    timer.start("Refinement", "Total");
    whfc_rb::KWayRefinerParallel<whfc_rb::PartitionThreadsafe, whfc::FlowHypergraphBuilder, whfc::Dinic, whfc_rb::HypergraphBuilderExtractor<whfc::FlowHypergraphBuilder, whfc_rb::PartitionThreadsafe>> refiner1(partition, timer, mt, config);
    //whfc_rb::KWayRefinerParallel<whfc_rb::PartitionThreadsafe, whfc_pr::LawlerFlowHypergraphParallel, whfc_pr::PushRelabelParallel, whfc_pr::LawlerFlowHypergraphBuilderExtractorParallel<whfc_pr::LawlerFlowHypergraphParallel, whfc_rb::PartitionThreadsafe>> refiner2(partition, timer, mt, config);
    //whfc_rb::KWayRefinerParallel<whfc_rb::PartitionThreadsafe, whfc_pr::LawlerFlowHypergraph, whfc_pr::PushRelabel, whfc_rb::HypergraphBuilderExtractor<whfc_pr::LawlerFlowHypergraph, whfc_rb::PartitionThreadsafe>> refiner3(partition, timer, mt, config);
    uint iterations = refiner1.refine(epsilon, maxNumIterations);
    timer.stop("Refinement");
    timer.stop("Total");

    printStatistics(partition, timer);
    std::cout << "numThreads: " << numThreads << std::endl;
    std::cout << "maxIterations: " << maxNumIterations << std::endl;
    std::cout << "iterations_done: " << iterations << std::endl;
}

