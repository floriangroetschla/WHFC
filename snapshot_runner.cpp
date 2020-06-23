#include "recursive_bisection/hypergraph.h"
#include "io/hmetis_io.h"
#include <iostream>
#include "util/timer.h"
#include "algorithm/hyperflowcutter.h"
#include "algorithm/dinic.h"
#include "io/whfc_io.h"

int main(int argc, const char *argv[]) {

    if (argc != 2) {
        throw std::runtime_error("Usage ./KWayRefinementParallel HypergraphFile");
    }
    whfc::FlowHypergraph hg = whfc::HMetisIO::readFlowHypergraph(argv[1]);


    whfc::HyperFlowCutter<whfc::Dinic> hfc(hg, 0);
    whfc::TimeReporter timer("Total");

    whfc::WHFC_IO::WHFCInformation i = whfc::WHFC_IO::readAdditionalInformation(argv[1]);
    whfc::WHFC_IO::readRandomGeneratorState(argv[1], hfc.cs.rng);
    hfc.cs.setMaxBlockWeight(0, i.maxBlockWeight[0]);
    hfc.cs.setMaxBlockWeight(1, i.maxBlockWeight[1]);
    hfc.upperFlowBound = i.upperFlowBound;

    bool hfc_result = hfc.enumerateCutsUntilBalancedOrFlowBoundExceeded(i.s, i.t);
    std::cout << "num_nodes: " << hg.numNodes() << std::endl;
    std::cout << "num_hyperedges: " << hg.numHyperedges() << std::endl;
    std::cout << "flow_value: " << hfc.cs.flowValue << std::endl;

    timer.merge(hfc.timer, "Total", "HyperFlowCutter");

    timer.report(std::cout);

    return 0;
}

