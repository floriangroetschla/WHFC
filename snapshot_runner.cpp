#include "recursive_bisection/hypergraph.h"
#include "io/hmetis_io.h"
#include <iostream>
#include "extern/patoh_wrapper.h"
#include <random>
#include "util/timer.h"
#include <tbb/tbb.h>
#include "recursive_bisection/k_way_refiner_parallel.h"
#include "recursive_bisection/partition_threadsafe.h"
#include "recursive_bisection/tbb_thread_pinning.h"
#include "recursive_bisection/config.h"
#include <filesystem>

int main(int argc, const char *argv[]) {

    if (argc != 3) {
        throw std::runtime_error("Usage ./KWayRefinementParallel HypergraphFile numThreads");
    }
    uint numThreads = std::stoi(argv[2]);

    tbb::task_scheduler_init init(numThreads);

    whfc_pr::LawlerFlowHypergraph hg = whfc::HMetisIO::readLawlerHypergraph(argv[1]);
    hg.maxNumHyperedges = hg.numHyperedges();
    hg.maxNumNodes = hg.numNodes();

    const whfc_rb::PartitionConfig config = {true, "D", true, true, 1, argv[1], 2};
    whfc::HyperFlowCutter<whfc_pr::PushRelabel, whfc_pr::LawlerFlowHypergraph> hfc(hg, 0, config);
    whfc::TimeReporter timer("Total");

    whfc::WHFC_IO::WHFCInformation i = whfc::WHFC_IO::readAdditionalInformation(argv[1]);
    whfc::WHFC_IO::readRandomGeneratorState(argv[1], hfc.cs.rng);
    hfc.cs.setMaxBlockWeight(0, i.maxBlockWeight[0]);
    hfc.cs.setMaxBlockWeight(1, i.maxBlockWeight[1]);
    hfc.upperFlowBound = i.upperFlowBound;
    hfc.timer.active = false;

    timer.start("Total");
    bool hfc_result = hfc.enumerateCutsUntilBalancedOrFlowBoundExceeded(i.s, i.t);
    timer.stop("Total");
    std::cout << "num_nodes: " << hg.numNodes() << std::endl;
    std::cout << "num_hyperedges: " << hg.numHyperedges() << std::endl;
    if (hfc_result) std::cout << "flow_value: " << hfc.cs.flowValue << std::endl;

    timer.report(std::cout);

    return 0;
}

