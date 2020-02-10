#include <stdexcept>
#include <math.h>       // Note(Lars): I think this guy is outdated
#include "extern/patoh_interface.h"
#include "io/hmetis_io.h"

// Note(Lars): Consider putting this into its own namespace / class so that it's easier to use from different binaries


std::vector<int> partition_recursively(const whfc::FlowHypergraph& hg, int seed, double epsilon, std::string preset, uint n_partitions) {
    struct FlowHypergraphDescription {
        std::vector<whfc::NodeWeight> nodeWeights;
        std::vector<whfc::HyperedgeWeight> hyperedgeWeights;
        std::vector<whfc::Node> pins;
        std::vector<whfc::PinIndex> hyperedgeSizes;
    };

        // rename to n_parts. typically we just say k.
    if (n_partitions == 1) {
        return std::vector<int>(hg.numNodes(), 0);
    }

    std::vector<int> numParts(2);   // Note(Lars): std::array<int, 2> will do the trick
    std::vector<int> partition;

    if (n_partitions % 2 == 0) {
        numParts[0] = n_partitions / 2;
        numParts[1] = n_partitions / 2;
        partition = PaToHInterface::bisectWithPatoh(hg, seed, epsilon, preset);
    } else {
        numParts[0] = n_partitions / 2;
        numParts[1] = numParts[0] + 1;
        partition = PaToHInterface::bisectImbalancedWithPatoh(hg, seed, float(numParts[1]) / float(numParts[0]), epsilon, preset);
    }

    // call WHFC here to improve the bisection
    
    if (n_partitions == 2) {
        return partition;
    } else {
        std::vector<FlowHypergraphDescription> hgds(2);

        std::vector<int> new_ids(partition.size());
        int carry0 = 0, carry1 = 0;
        for (uint i = 0; i < partition.size(); ++i) {
            if (partition[i] == 0) {
                new_ids[i] = carry0;	// Note(Lars): can be one-liner
                carry0++;
            } else {
                new_ids[i] = carry1;
                carry1++;
            }
        }

        for (whfc::Hyperedge e : hg.hyperedgeIDs()) {
            std::vector<whfc::Node> pins;
            int partitionID = partition[hg.getPin(hg.beginIndexPins(e)).pin.value()];	// Note(Lars): dangerous. the hyperedge may have zero pins
            for (auto pin : hg.pinIndices(e)) {
                pins.push_back(hg.getPin(pin).pin);
                if (partitionID != partition[hg.getPin(pin).pin.value()]) {
                    partitionID = -1;
                    // Note(Lars): you could break here
                }
            }
            
            // Note(Lars): We're going for the connectivity objective (not cut) --> throwing away split hyperedges is not okay.
            if (partitionID != -1) {
                hgds[partitionID].hyperedgeSizes.push_back(hg.pinCount(e));
                hgds[partitionID].hyperedgeWeights.push_back(whfc::HyperedgeWeight(hg.capacity(e)));
                for (whfc::Node pin : pins) {
                    hgds[partitionID].pins.push_back(whfc::Node(new_ids[pin.value()]));
                }
            }
        }

        for (whfc::Node node : hg.nodeIDs()) {
            hgds[partition[node.value()]].nodeWeights.push_back(hg.nodeWeight(node));
        }

        std::vector<std::vector<int>> sub_partitions(2);        // std::array< std::vector<int>, 2 > will do the trick as well, but put the member variables of the inner vector on the stack
        for (uint i = 0; i < 2; ++i) {
            whfc::FlowHypergraph flowHG = whfc::FlowHypergraph(hgds[i].nodeWeights, hgds[i].hyperedgeWeights, hgds[i].hyperedgeSizes, hgds[i].pins);
            // Note(Lars): can't we relax the epsilon at each step?
            sub_partitions[i] = partition_recursively(flowHG, seed, epsilon, preset, numParts[i]);
        }

        std::vector<int> result(hg.numNodes()); // Note(Lars): could reuse the partition vector used for PaToH
        for (uint i = 0; i < partition.size(); ++i) {
            result[i] = sub_partitions[partition[i]][new_ids[i]] + partition[i] * numParts[0];
        }

        return result;
    }
}

int main(int argc, const char* argv[]) {
    if (argc != 6) {
        throw std::runtime_error("Must provide four arguments. Hypergraph, epsilon, seed, PaToH-preset (Q(uality) or D(default)), number of partitions.");
    }

    std::string path_hg = argv[1];
    std::string str_eps = argv[2];
    std::string str_seed = argv[3];
    std::string preset = argv[4];
    std::string str_n_partitions = argv[5];
    double epsilon = std::stod(str_eps);
    int seed = std::stoi(str_seed);
    int n_partitions = std::stoi(str_n_partitions);

    epsilon = pow(1.0 + epsilon, 1.0 / log2(n_partitions)) - 1.0;

    whfc::FlowHypergraph hg = whfc::HMetisIO::readFlowHypergraph(path_hg);

    auto t_begin = whfc::time_now();
    std::vector<int> partition = partition_recursively(hg, seed, epsilon, preset, n_partitions);
    auto t_end = whfc::time_now();

    std::cout
            << "Time: "
            << whfc::inSeconds(std::chrono::duration<double, std::micro>(t_end - t_begin)).count()
            << std::endl;
}