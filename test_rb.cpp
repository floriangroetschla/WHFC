#include "recursive_bisection/hypergraph.h"
#include "io/hmetis_io.h"
#include <iostream>
#include "extern/patoh_wrapper.h"
#include "recursive_bisection/recursive_bisection.h"
#include <random>
#include "util/timer.h"

int main() {
    //whfc_rb::CSRHypergraph hg = whfc::HMetisIO::readCSRHypergraph("../test_hypergraphs/testhg.hgr");
    whfc_rb::CSRHypergraph hg = whfc::HMetisIO::readCSRHypergraph("../../benchmark_set/192bit.mtx.hgr");
    //whfc_rb::CSRHypergraph hg = whfc::HMetisIO::readCSRHypergraph("../test_hypergraphs/twocenters.hgr");

    std::mt19937 mt(42);

    whfc::TimeReporter timer;

    whfc_rb::RecursiveBisector recursive_bisector = whfc_rb::RecursiveBisector(hg.numNodes(), hg.numHyperedges(), hg.numPins(), mt);
    timer.start("Recursive bisector");
    whfc_rb::Partition partition = recursive_bisector.run(hg, 0.1, "D", 5);
    timer.stop("Recursive bisector");
    timer.report(std::cout);
    std::cout << "Imbalance: " << partition.imbalance(hg) << std::endl;

}