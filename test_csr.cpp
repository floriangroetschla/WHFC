#include "recursive_bisection/hypergraph.h"
#include "io/hmetis_io.h"
#include <iostream>
#include "extern/patoh_wrapper.h"
#include "recursive_bisection/recursive_bisection.h"

int main() {
	whfc_rb::CSRHypergraph hg = whfc::HMetisIO::readCSRHypergraph("../test_hypergraphs/testhg.hgr");

	/*std::vector<int> result = PaToHInterface::bisectWithPatoh(hg, 42, 0.1, "D", true, false);
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

    std::vector<int> partition = whfc_rb::RecursiveBisector::run(hg, 42, 0.1, "D", 5);


    for (uint i = 0; i < partition.size(); ++i) {
        std::cout << partition[i] << " ";
    }
    std::cout << std::endl;
}