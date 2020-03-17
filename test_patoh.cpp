#include "recursive_bisection/hypergraph.h"
#include "io/hmetis_io.h"
#include <iostream>
#include "extern/patoh_wrapper.h"
#include "util/timer.h"

int main(int argc, const char *argv[]) {

    if (argc != 7) {
        throw std::runtime_error("Usage ./PaToH HypergraphFile epsilon k seed preset objective");
    }
    std::string graph_path = argv[1];
    double epsilon = std::stod(argv[2]);
    uint numParts = std::stoul(argv[3]);
    int seed = std::stoi(argv[4]);
    std::string patoh_preset = argv[5];
    std::string objective = argv[6];

    whfc_rb::CSRHypergraph hg = whfc::HMetisIO::readCSRHypergraph(graph_path);
    
    whfc::TimeReporter timer("Total");
    whfc_rb::PartitionBase partition(numParts, hg);

    timer.start("Total");
    PaToHInterface::partitionWithPatoh(partition, seed, numParts, epsilon, patoh_preset, objective);
    timer.stop("Total");
//    timer.report(std::cout);
    //partition.print(std::cout);

    std::string graph_name = graph_path.substr(graph_path.find_last_of('/') + 1);
    std::cout   << "PaToH-" << patoh_preset << ","
                << graph_name << ","
                << "no" << ","
                << seed << ","
                << partition.numParts() << ","
                << epsilon << ","
                << partition.imbalance() << ","
                << whfc::inSeconds(timer.get("Total")).count() << ","
                << objective << ","
                << partition.km1Objective() << ","
                << partition.cutObjective() << ","
                << "no" << std::endl;

//    std::cout << "Imbalance: " << partition.imbalance() << std::endl;
//    std::cout << "Num parts: " << partition.numParts() << std::endl;
//    std::cout << "km1: " << partition.km1Objective() << std::endl;
//    std::cout << "cut: " << partition.cutObjective() << std::endl;
}
