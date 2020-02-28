#include "recursive_bisection/hypergraph.h"
#include "io/hmetis_io.h"
#include <iostream>
#include "extern/patoh_wrapper.h"
#include "recursive_bisection/recursive_bisection.h"
#include "util/timer.h"

int main(int argc, const char* argv[]) {

    if (argc != 6) {
        throw std::runtime_error("Usage ./PaToH HypergraphFile epsilon k seed preset");
    }
    whfc_rb::CSRHypergraph hg = whfc::HMetisIO::readCSRHypergraph(argv[1]);
    double epsilon = std::stod(argv[2]);
    uint numParts = std::stoul(argv[3]);
    int seed = std::stoi(argv[4]);
    std::string patoh_preset = argv[5];


    whfc::TimeReporter timer;
    whfc_rb::PartitionBase partition(numParts, hg);

    timer.start("PaToH");
    PaToHInterface::partitionWithPatoh(partition, seed, numParts, epsilon, patoh_preset);
    timer.stop("PaToH");
    timer.report(std::cout);
    //partition.print(std::cout);
    std::cout << "Imbalance: " << partition.imbalance() << std::endl;
    std::cout << "Num parts: " << partition.numParts() << std::endl;
    std::cout << "km1: " << partition.km1Objective() << std::endl;
    std::cout << "cut: " << partition.cutObjective() << std::endl;
}