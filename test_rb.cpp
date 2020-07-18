#include "recursive_bisection/hypergraph.h"
#include "io/hmetis_io.h"
#include <iostream>
#include "extern/patoh_wrapper.h"
#include "recursive_bisection/recursive_bisection.h"
#include <random>
#include "util/timer.h"
#include "recursive_bisection/partition_base.h"
#include "recursive_bisection/fhgb_extraction.h"


void printStatistics(whfc_rb::PartitionBase &partition, whfc::TimeReporter &timer) {
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

    whfc::TimeReporter timer("Total");

    whfc_rb::PartitionerConfig config;
    config.patoh_preset = patoh_preset;

    if (!mode.compare("RBONLY")) {
        std::cout << "Using mode RBONLY" << std::endl;
        config.refine = false;
    } else if (!mode.compare("RBWHFC")) {
        std::cout << "Using mode RBWHFC" << std::endl;
        config.refine = true;
    } else {
        throw std::runtime_error("Mode must be one of: RBONlY, RBWHFC");
    }

    whfc_rb::RecursiveBisector recursive_bisector = whfc_rb::RecursiveBisector<whfc_rb::PartitionThreadsafe, whfc_pr::LawlerFlowHypergraph, whfc_pr::PushRelabel, whfc_rb::HypergraphBuilderExtractor<whfc_pr::LawlerFlowHypergraph, whfc_rb::PartitionThreadsafe>>(hg.numNodes(), hg.numHyperedges(),
                                                                               hg.numPins(), mt, timer, config);
    timer.start("Total");
    whfc_rb::PartitionBase partition = recursive_bisector.run(hg, epsilon, numParts);
    timer.stop("Total");
    printStatistics(partition, timer);
}

