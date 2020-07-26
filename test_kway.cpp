#include "datastructure/hypergraph.h"
#include "io/hmetis_io.h"
#include <iostream>
#include "extern/patoh_wrapper.h"
#include "partitioner/recursive_bisection.h"
#include <random>
#include "util/timer.h"
#include "datastructure/partition_base.h"
#include "datastructure/partition_ca.h"
#include "partitioner/null_refiner.h"
#include "partitioner/k_way_refiner.h"
#include "push_relabel/lawler_fhgb_extraction_parallel.h"

void printStatistics(whfc_rb::PartitionBase &partition, whfc::TimeReporter &timer) {
    timer.report(std::cout);
    std::cout << "Imbalance: " << partition.imbalance() << std::endl;
    std::cout << "Num_parts: " << partition.numParts() << std::endl;
    std::cout << "km1: " << partition.km1Objective() << std::endl;
    std::cout << "cut: " << partition.cutObjective() << std::endl;
}

int main(int argc, const char *argv[]) {

    if (argc != 7) {
        throw std::runtime_error("Usage ./KWayRefinement HypergraphFile epsilon k seed preset maxIterations");
    }
    whfc_rb::CSRHypergraph hg = whfc::HMetisIO::readCSRHypergraph(argv[1]);
    double epsilon = std::stod(argv[2]);
    uint numParts = std::stoul(argv[3]);
    int seed = std::stoi(argv[4]);
    std::string patoh_preset = argv[5];
    uint maxIterations = std::stoi(argv[6]);
    std::mt19937 mt(seed);

    whfc_rb::PartitionerConfig config = {true, patoh_preset, true, false};
    whfc::TimeReporter timer("Total");

    timer.start("Total");
    timer.start("PaToH", "Total");
    whfc_rb::PartitionThreadsafe partition(numParts, hg);
    PaToHInterface::partitionWithPatoh(partition, seed, numParts, epsilon, patoh_preset);
    timer.stop("PaToH");
    partition.initialize();

    timer.start("Refinement", "Total");
    whfc_rb::KWayRefiner<whfc_rb::PartitionThreadsafe, whfc_pr::LawlerFlowHypergraphParallel, whfc_pr::PushRelabelParallel, whfc_pr::LawlerFlowHypergraphBuilderExtractorParallel<whfc_pr::LawlerFlowHypergraphParallel, whfc_rb::PartitionThreadsafe>> refiner(partition, timer, mt, config);
    refiner.refine(epsilon, maxIterations);
    timer.stop("Refinement");
    timer.stop("Total");

    printStatistics(partition, timer);
}

