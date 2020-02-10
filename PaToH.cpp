#include <stdexcept>
#include <iostream>

#include "extern/patoh_interface.h"
#include "io/hmetis_io.h"
#include "util/range.h"
#include "datastructure/flow_hypergraph.h"

int main(int argc, const char* argv[]) {
    if (argc != 6) {
        throw std::runtime_error("Must provide four arguments. Hypergraph, epsilon, seed, PaToH-preset (Q(uality) or D(default)), number of partitions");
    }

    std::string path_hg = argv[1];
    std::string str_eps = argv[2];
    std::string str_seed = argv[3];
    std::string preset = argv[4];
    std::string str_n_partitions = argv[5];
    double epsilon = std::stod(str_eps);
    int seed = std::stoi(str_seed);
    int n_partitions = std::stoi(str_n_partitions);

    whfc::FlowHypergraph hg = whfc::HMetisIO::readFlowHypergraph(path_hg);

    auto t_begin = whfc::time_now();
    std::vector<int> result = PaToHInterface::partitionWithPatoh(hg, seed, n_partitions, epsilon, preset);
    //std::vector<int> result = PaToHInterface::bisectImbalancedWithPatoh(hg, seed, 2.0, epsilon, preset);
    //std::vector<int> result = PaToHInterface::bisectWithPatoh(hg, seed, epsilon, preset);
    auto t_end = whfc::time_now();

    std::cout
        << "Time: "
        << whfc::inSeconds(std::chrono::duration<double, std::micro>(t_end - t_begin)).count()
        << std::endl;

}