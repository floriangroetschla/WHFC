#include "recursive_bisection/hypergraph.h"
#include "io/hmetis_io.h"
#include <iostream>
#include "extern/patoh_wrapper.h"
#include "recursive_bisection/recursive_bisection.h"
#include "recursive_bisection/fhgb_extraction.h"

int main() {
    whfc_rb::CSRHypergraph hg = whfc::HMetisIO::readCSRHypergraph("../test_hypergraphs/testhg.hgr");
    //whfc_rb::CSRHypergraph hg = whfc::HMetisIO::readCSRHypergraph("../../benchmark_set/192bit.mtx.hgr");

    whfc_rb::Partition partition = PaToHInterface::bisectWithPatoh(hg, 42, 0.1, "D", true, true);
    /*
    for (uint i = 0; i < result.size(); ++i) {
        std::cout << result[i] << " ";
    }
    std::cout << std::endl;

    result = PaToHInterface::bisectWithPatoh(hg, 42, 0.1, "D", false, false);
    for (uint i = 0; i < result.size(); ++i) {
        std::cout << result[i] << " ";
    }
    std::cout << std::endl;

    PaToHInterface::freePatoh();*/

    //std::vector<int> partition = whfc_rb::RecursiveBisector::run(hg, 42, 0.1, "D", 5);

    partition.print(std::cout);


    /*whfc_rb::FlowHypergraphBuilderExtractor fhge(hg.numNodes(), hg.numHyperedges(), hg.numPins());
    fhge.run(hg, cut_hes, partition, 100);
    fhge.fhgb.printHypergraph(std::cout);*/

    whfc_rb::RecursiveBisector recursive_bisector = whfc_rb::RecursiveBisector(hg.numNodes(), hg.numHyperedges(), hg.numPins(), 42);
    recursive_bisector.run(hg, 42, 0.1, "D", 4);
}