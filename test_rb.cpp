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

void printStatistics(whfc_rb::PartitionBase& partition, whfc::TimeReporter& timer) {
    timer.report(std::cout);
    std::cout << "Imbalance: " << partition.imbalance() << std::endl;
    std::cout << "Num_parts: " << partition.numParts() << std::endl;
    std::cout << "km1: " << partition.km1Objective() << std::endl;
    std::cout << "cut: " << partition.cutObjective() << std::endl;
}

int main(int argc, const char *argv[]) {

    if (argc != 7) {
        throw std::runtime_error("Usage ./RecursiveBisection HypergraphFile epsilon k seed preset mode");
    }
    whfc_rb::CSRHypergraph hg = whfc::HMetisIO::readCSRHypergraph(argv[1]);
    double epsilon = std::stod(argv[2]);
    uint numParts = std::stoul(argv[3]);
    int seed = std::stoi(argv[4]);
    std::string patoh_preset = argv[5];
    std::string mode = argv[6];
    std::mt19937 mt(seed);

    whfc::TimeReporter timer("RecursiveBisector");

    if (!mode.compare("RBONLY")) {
        std::cout << "Using mode RBONLY" << std::endl;
        whfc_rb::RecursiveBisector recursive_bisector = whfc_rb::RecursiveBisector<whfc_rb::NullRefiner>(
                hg.numNodes(), hg.numHyperedges(), hg.numPins(), mt, timer);
        timer.start("RecursiveBisector");
        whfc_rb::PartitionBase partition = recursive_bisector.run<whfc_rb::PartitionBase>(hg, epsilon, patoh_preset, numParts);
        timer.stop("RecursiveBisector");
        printStatistics(partition, timer);
    } else if (!mode.compare("RBPARTITIONBASE")) {
        std::cout << "Using mode RBPARTITIONBASE" << std::endl;
        whfc_rb::RecursiveBisector recursive_bisector = whfc_rb::RecursiveBisector<whfc_rb::WHFCRefinerTwoWay>(
                hg.numNodes(), hg.numHyperedges(), hg.numPins(), mt, timer);
        timer.start("RecursiveBisector");
        whfc_rb::PartitionBase partition = recursive_bisector.run<whfc_rb::PartitionBase>(hg, epsilon, patoh_preset, numParts);
        timer.stop("RecursiveBisector");
        printStatistics(partition, timer);
    } else if (!mode.compare("RBPARTITIONCA")) {
        std::cout << "Using mode RBPARTITIONCA" << std::endl;
        whfc_rb::RecursiveBisector recursive_bisector = whfc_rb::RecursiveBisector<whfc_rb::WHFCRefinerTwoWay>(
                hg.numNodes(), hg.numHyperedges(), hg.numPins(), mt, timer);
        timer.start("RecursiveBisector");
        whfc_rb::PartitionCA partition = recursive_bisector.run<whfc_rb::PartitionCA>(hg, epsilon, patoh_preset, numParts);
        timer.stop("RecursiveBisector");
        printStatistics(partition, timer);
    } else {
        throw std::runtime_error("Mode must be one of: RBONLY, RBPARTITIONBASE, RBPARTITIONCA");
    }

}

