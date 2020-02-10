#include "recursive_bisection/hypergraph.h"
#include "io/hmetis_io.h"
#include <iostream>

int main() {
	whfc_rb::CSRHypergraph hg = whfc::HMetisIO::readCSRHypergraph("../test_hypergraphs/twocenters.hgr");
	std::cout << hg;
}