#include <iostream>
#include "datastructure/hypergraph.h"
#include "util/range.h"
#include "datastructure/flow_hypergraph.h"
#include "algorithm/hyperflowcutter.h"
#include "algorithm/grow_assimilated.h"
#include "algorithm/ford_fulkerson.h"
#include "io/hmetis_io.h"


namespace whfc {
	void run(const std::string& filename, Node s, Node t) {
		FlowHypergraph hg = HMetisIO::readFlowHypergraph(filename);
		NodeWeight mbw(hg.totalNodeWeight() / 2);

		if (s >= hg.numNodes() || t >= hg.numNodes())
			throw std::runtime_error("s or t not within node id range");
		
		HyperFlowCutter<BasicEdmondsKarp> hfc(hg, mbw);
		hfc.initialize(s,t);
		hfc.runUntilBalanced();
	}
}

int main(int argc, const char* argv[]) {
	if (argc != 4)
		throw std::runtime_error("Usage: ./WHFC hypergraphfile s t");
	std::string hgfile = argv[1];
	whfc::Node s(static_cast<unsigned int>(std::stoul(argv[2])));
	whfc::Node t(static_cast<unsigned int>(std::stoul(argv[3])));
	whfc::run(hgfile, s, t);
	return 0;
}