#include "recursive_bisection/hypergraph.h"
#include "io/hmetis_io.h"
#include <iostream>
#include "extern/patoh_wrapper.h"
#include "recursive_bisection/recursive_bisection.h"
#include <random>
#include "util/timer.h"

int main(int argc, const char* argv[]) {
    
    if (argc != 4) {
        throw std::runtime_error("Usage ./RecursiveBisection HypergraphFile epsilon k");
    }
    //whfc_rb::CSRHypergraph hg = whfc::HMetisIO::readCSRHypergraph("../test_hypergraphs/testhg.hgr");
    whfc_rb::CSRHypergraph hg = whfc::HMetisIO::readCSRHypergraph(argv[1]);
    double epsilon = std::stod(argv[2]);
    uint numParts = std::stoul(argv[3]);
    std::string patoh_preset = "D";
    //whfc_rb::CSRHypergraph hg = whfc::HMetisIO::readCSRHypergraph("../test_hypergraphs/twocenters.hgr");

    std::mt19937 mt(42);

    whfc::TimeReporter timer("Recursive bisector");

    whfc_rb::RecursiveBisector recursive_bisector = whfc_rb::RecursiveBisector(hg.numNodes(), hg.numHyperedges(), hg.numPins(), mt, timer);
    timer.start("Recursive bisector");
    whfc_rb::Partition partition = recursive_bisector.run(hg, epsilon, patoh_preset, numParts);
    timer.stop("Recursive bisector");
    timer.report(std::cout);
    //partition.print(std::cout);
    std::cout << "Imbalance: " << partition.imbalance(hg) << std::endl;
    std::cout << "Num parts: " << partition.numParts() << std::endl;
    std::cout << "km1: " << partition.km1Objective(hg) << std::endl;
    std::cout << "cut: " << partition.cutObjective(hg) << std::endl;
}