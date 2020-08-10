#include "datastructure/hypergraph.h"
#include "io/hmetis_io.h"
#include <iostream>
#include "extern/patoh_wrapper.h"
#include "partitioner/recursive_bisection.h"
#include <random>
#include "util/timer.h"
#include "datastructure/partition_base.h"
#include "datastructure/partition_ca.h"
#include "partitioner/null_refiner.h"
#include "partitioner/k_way_refiner.h"
#include "partitioner/fhgb_extraction_parallel.h"
#include "push_relabel/lawler_fhgb_extraction_parallel.h"
#include "partitioner/tbb_thread_pinning.h"
#include "flow_algorithms/dinic_thread_local_vectors.h"
#include "partitioner/fhgb_extraction.h"
#include "flow_algorithms/dinic_parallel.h"
#include <version.h>
#include <tbb/task_scheduler_init.h>

void printStatistics(whfc_rb::PartitionBase &partition, whfc::TimeReporter &timer) {
    timer.report(std::cout);
    std::cout << "Imbalance: " << partition.imbalance() << std::endl;
    std::cout << "Num_parts: " << partition.numParts() << std::endl;
    std::cout << "km1: " << partition.km1Objective() << std::endl;
    std::cout << "cut: " << partition.cutObjective() << std::endl;
}

int main(int argc, const char *argv[]) {

    if (argc != 11) {
        throw std::runtime_error("Usage ./KWayRefinement HypergraphFile epsilon k seed preset numThreads useThreadPinning(0 or 1) distancePiercing(0 or 1) maxNumIterations mode");
    }
    whfc_rb::CSRHypergraph hg = whfc::HMetisIO::readCSRHypergraph(argv[1]);
    double epsilon = std::stod(argv[2]);
    uint numParts = std::stoul(argv[3]);
    int seed = std::stoi(argv[4]);
    std::string patoh_preset = argv[5];
    uint numThreads = std::stoi(argv[6]);
    bool useThreadPinning = std::stoi(argv[7]);
    bool distancePiercing =     std::stoi(argv[8]);
    int maxNumIterations = std::stoi(argv[9]);
    std::mt19937 mt(seed);
    std::string mode = argv[10];

    tbb::task_scheduler_init init(numThreads);
    whfc_rb::pinning_observer pinner;
    pinner.observe(useThreadPinning);

    bool precomputeCuts = true;

    whfc_rb::PartitionerConfig config = {true, patoh_preset, precomputeCuts, distancePiercing, numThreads, "unset", numParts};
    whfc::TimeReporter timer("Total");

    std::cout << "useThreadPinning: " << useThreadPinning << std::endl;
    std::cout << "precomputeCuts: " << precomputeCuts << std::endl;
    std::cout << "distancePiercing: " << distancePiercing << std::endl;
    std::cout << "numThreads: " << numThreads << std::endl;
    std::cout << "commit: " << GIT_COMMIT_HASH << std::endl;

    timer.start("Total");
    timer.start("PaToH", "Total");
    whfc_rb::PartitionThreadsafe partition(numParts, hg);
    PaToHInterface::partitionWithPatoh(partition, seed, numParts, epsilon, patoh_preset);
    timer.stop("PaToH");
    partition.initialize();

    uint numIterations = 0;

    timer.start("Refinement", "Total");

    if (!mode.compare("seqExtraction_seqDinic")) {
        std::cout << "mode seqExtraction_seqDinic" << std::endl;
        whfc_rb::KWayRefiner<whfc_rb::PartitionThreadsafe, whfc::FlowHypergraphBuilder, whfc::Dinic, whfc_rb::HypergraphBuilderExtractor<whfc::FlowHypergraphBuilder, whfc_rb::PartitionThreadsafe>> refiner(partition, timer, mt, config);
        numIterations = refiner.refine(epsilon, maxNumIterations);
    } else if (!mode.compare("parExtraction_seqDinic")) {
        std::cout << "mode parExtraction_seqDinic" << std::endl;
        whfc_rb::KWayRefiner<whfc_rb::PartitionThreadsafe, whfc::FlowHypergraphBuilder, whfc::Dinic, whfc_rb::FlowHypergraphBuilderExtractorParallel<whfc::FlowHypergraphBuilder, whfc_rb::PartitionThreadsafe>> refiner(partition, timer, mt, config);
        numIterations = refiner.refine(epsilon, maxNumIterations);
    } else if (!mode.compare("seqExtraction_parDinic")) {
        std::cout << "mode seqExtraction_parDinic" << std::endl;
        whfc_rb::KWayRefiner<whfc_rb::PartitionThreadsafe, whfc::FlowHypergraphBuilder, whfc::DinicThreadLocalVectors, whfc_rb::HypergraphBuilderExtractor<whfc::FlowHypergraphBuilder, whfc_rb::PartitionThreadsafe>> refiner(partition, timer, mt, config);
        numIterations = refiner.refine(epsilon, maxNumIterations);
    } else if (!mode.compare("parExtraction_parDinic")) {
        std::cout << "mode parExtraction_parDinic" << std::endl;
        whfc_rb::KWayRefiner<whfc_rb::PartitionThreadsafe, whfc::FlowHypergraphBuilder, whfc::DinicThreadLocalVectors, whfc_rb::FlowHypergraphBuilderExtractorParallel<whfc::FlowHypergraphBuilder, whfc_rb::PartitionThreadsafe>> refiner(partition, timer, mt, config);
        numIterations = refiner.refine(epsilon, maxNumIterations);
    } else if (!mode.compare("seqExtraction_seqPR")) {
        std::cout << "mode seqExtraction_seqPR" << std::endl;
        whfc_rb::KWayRefiner<whfc_rb::PartitionThreadsafe, whfc_pr::LawlerFlowHypergraph, whfc_pr::PushRelabel, whfc_rb::HypergraphBuilderExtractor<whfc_pr::LawlerFlowHypergraph, whfc_rb::PartitionThreadsafe>> refiner(partition, timer, mt, config);
        numIterations = refiner.refine(epsilon, maxNumIterations);
    } else if (!mode.compare("parExtraction_seqPR")) {
        std::cout << "mode parExtraction_seqPR" << std::endl;
        whfc_rb::KWayRefiner<whfc_rb::PartitionThreadsafe, whfc_pr::LawlerFlowHypergraph, whfc_pr::PushRelabel, whfc_pr::LawlerFlowHypergraphBuilderExtractorParallel<whfc_pr::LawlerFlowHypergraph, whfc_rb::PartitionThreadsafe>> refiner(partition, timer, mt, config);
        numIterations = refiner.refine(epsilon, maxNumIterations);
    } else if (!mode.compare("seqExtraction_parPR")) {
        std::cout << "mode seqExtraction_parPR" << std::endl;
        whfc_rb::KWayRefiner<whfc_rb::PartitionThreadsafe, whfc_pr::LawlerFlowHypergraphParallel, whfc_pr::PushRelabelParallel, whfc_rb::HypergraphBuilderExtractor<whfc_pr::LawlerFlowHypergraphParallel, whfc_rb::PartitionThreadsafe>> refiner(partition, timer, mt, config);
        numIterations = refiner.refine(epsilon, maxNumIterations);
    } else if (!mode.compare("parExtraction_parPR")) {
        std::cout << "mode parExtraction_parPR" << std::endl;
        whfc_rb::KWayRefiner<whfc_rb::PartitionThreadsafe, whfc_pr::LawlerFlowHypergraphParallel, whfc_pr::PushRelabelParallel, whfc_pr::LawlerFlowHypergraphBuilderExtractorParallel<whfc_pr::LawlerFlowHypergraphParallel, whfc_rb::PartitionThreadsafe>> refiner(partition, timer, mt, config);
        numIterations = refiner.refine(epsilon, maxNumIterations);
    } else {
        throw std::runtime_error("Mode must be one of: "
                                 "seqExtraction_seqDinic, "
                                 "parExtraction_seqDinic, "
                                 "seqExtraction_parDinic, "
                                 "parExtraction_parDinic, "
                                 "seqExtraction_seqPR, "
                                 "parExtraction_seqPR, "
                                 "seqExtraction_parPR, "
                                 "parExtraction_parPR");
    }

    timer.stop("Refinement");
    timer.stop("Total");

    std::cout << "iterations_done: " << numIterations << std::endl;
    printStatistics(partition, timer);
}

