#include "recursive_bisection/hypergraph.h"
#include "io/hmetis_io.h"
#include <iostream>
#include "extern/patoh_wrapper.h"
#include "recursive_bisection/recursive_bisection.h"
#include <random>
#include "util/timer.h"
#include "recursive_bisection/partition_base.h"
#include "recursive_bisection/partition_ca.h"
#include "recursive_bisection/null_refiner.h"
#include "recursive_bisection/k_way_refiner.h"

void printStatistics(whfc_rb::PartitionBase &partition, whfc::TimeReporter &timer) {
    timer.report(std::cout);
    std::cout << "Imbalance: " << partition.imbalance() << std::endl;
    std::cout << "Num_parts: " << partition.numParts() << std::endl;
    std::cout << "km1: " << partition.km1Objective() << std::endl;
    std::cout << "cut: " << partition.cutObjective() << std::endl;
}

int main(int argc, const char *argv[]) {

    if (argc != 6) {
        throw std::runtime_error("Usage ./KWayRefinement HypergraphFile epsilon k seed preset");
    }
    whfc_rb::CSRHypergraph hg = whfc::HMetisIO::readCSRHypergraph(argv[1]);
    double epsilon = std::stod(argv[2]);
    uint numParts = std::stoul(argv[3]);
    int seed = std::stoi(argv[4]);
    std::string patoh_preset = argv[5];
    std::mt19937 mt(seed);

    whfc::TimeReporter timer("Total");

    timer.start("Total");
    timer.start("PaToH", "Total");
    whfc_rb::PartitionCA partition(numParts, hg);
    PaToHInterface::partitionWithPatoh(partition, seed, numParts, epsilon, patoh_preset);
    timer.stop("PaToH");

    timer.start("Refinement", "Total");
    whfc_rb::KWayRefiner refiner(partition, timer, mt);
    refiner.refine(epsilon);
    timer.stop("Refinement");
    timer.stop("Total");

    printStatistics(partition, timer);
}

