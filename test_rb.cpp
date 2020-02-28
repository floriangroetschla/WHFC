#include "recursive_bisection/hypergraph.h"
#include "io/hmetis_io.h"
#include <iostream>
#include "extern/patoh_wrapper.h"
#include "recursive_bisection/recursive_bisection.h"
#include <random>
#include "util/timer.h"
#include "recursive_bisection/partition_base.h"
#include "recursive_bisection/partition_ca.h"

int main(int argc, const char *argv[]) {

    if (argc != 6) {
        throw std::runtime_error("Usage ./RecursiveBisection HypergraphFile epsilon k seed preset");
    }
    //whfc_rb::CSRHypergraph hg = whfc::HMetisIO::readCSRHypergraph("../test_hypergraphs/testhg.hgr");
    whfc_rb::CSRHypergraph hg = whfc::HMetisIO::readCSRHypergraph(argv[1]);
    double epsilon = std::stod(argv[2]);
    uint numParts = std::stoul(argv[3]);
    int seed = std::stoi(argv[4]);
    std::string patoh_preset = argv[5];
    //whfc_rb::CSRHypergraph hg = whfc::HMetisIO::readCSRHypergraph("../test_hypergraphs/twocenters.hgr");

    std::mt19937 mt(seed);

    whfc::TimeReporter timer("RecursiveBisector");

    //whfc_rb::PartitionBase part(10, hg);

    whfc_rb::RecursiveBisector recursive_bisector = whfc_rb::RecursiveBisector<whfc_rb::WHFCRefinerTwoWay>(
            hg.numNodes(), hg.numHyperedges(), hg.numPins(), mt, timer);
    timer.start("RecursiveBisector");
    whfc_rb::PartitionCA partition = recursive_bisector.run<whfc_rb::PartitionCA>(hg, epsilon, patoh_preset, numParts);
    timer.stop("RecursiveBisector");
    timer.report(std::cout);
    //partition.print(std::cout);
    std::cout << "Imbalance: " << partition.imbalance() << std::endl;
    std::cout << "Num parts: " << partition.numParts() << std::endl;
    std::cout << "km1: " << partition.km1Objective() << std::endl;
    std::cout << "cut: " << partition.cutObjective() << std::endl;
}