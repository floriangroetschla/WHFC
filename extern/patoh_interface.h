#pragma once

#include "../datastructure/flow_hypergraph.h"
#include "patoh.h"

using whfc::FlowHypergraph;

class PaToHInterface {
public:

	struct PatohParameters {
		explicit PatohParameters(double epsilon=0.0) : use_target_weights(false), target_weights(), epsilon(epsilon), k(2) {}
		bool use_target_weights;
		std::vector<float> target_weights;
		double epsilon;
		int k;
	};

    static std::vector<int> bisectImbalancedWithPatoh(const FlowHypergraph& hg,
                                                                              int seed,
                                                                              float imbalanceFactor,
                                                                              double epsilon=0.05,
                                                                              std::string preset = "D") {
        PatohParameters p;
        p.use_target_weights = true;
        p.k = 2;
        p.epsilon = epsilon;
        p.target_weights = { float(1), float(imbalanceFactor) };
        return runPatoh(hg, seed, p, preset);
    }

	static std::vector<int> bisectWithPatoh(const FlowHypergraph &hg,
												  int seed,
												  double epsilon=0.0,
												  std::string preset = "D") {
		return runPatoh(hg, seed, PatohParameters(epsilon), preset);
	}

    static std::vector<int> partitionWithPatoh(const FlowHypergraph& hg, int seed, int numPartitions, double epsilon=0.05,
                                                                                    std::string preset = "D") {
        PatohParameters p;
        p.use_target_weights = false;
        p.k = numPartitions;
        p.epsilon = epsilon;
        return runPatoh(hg, seed, p, preset);
    }


	static std::vector<int> runPatoh(const FlowHypergraph& hg, int seed, PatohParameters params, std::string str_preset = "D") {
		std::vector<int> vec_partition(hg.numNodes());
		std::vector<int> vec_partweights(params.k, 0);

		std::vector<int> vec_cwghts(hg.numNodes(), 1);

		std::vector<int> vec_nwghts(hg.numHyperedges(), 1);

		std::vector<int> vec_pins(hg.numPins());
		std::vector<int> vec_xpins(hg.numHyperedges() + 1);

        {
            for (whfc::Hyperedge e : hg.hyperedgeIDs()) {
                vec_xpins[e.value()] = hg.beginIndexPins(e).value();
                for (auto pin : hg.pinIndices(e)) {
                    vec_pins[pin.value()] = hg.getPin(pin).pin.value();
                }
                vec_nwghts[e.value()] = hg.capacity(e);
            }

            for (whfc::Node node : hg.nodeIDs()) {
                vec_cwghts[node.value()] = hg.nodeWeight(node).value();
            }

            vec_xpins[hg.numHyperedges()] = hg.endIndexPins(whfc::Hyperedge(hg.numHyperedges() - 1)).value();
        }

		PaToH_Parameters args;
		int preset = -1;
		if (str_preset == "D") { preset = PATOH_SUGPARAM_DEFAULT; }
		else if (str_preset == "Q") { preset = PATOH_SUGPARAM_QUALITY; }
		else if (str_preset == "S") { preset = PATOH_SUGPARAM_SPEED; }
		else { throw std::runtime_error("Unknown PaToH preset" + str_preset); }

		PaToH_Initialize_Parameters(&args, PATOH_CUTPART, preset);

		args._k = params.k;
		args.doinitperm = 0;
		args.final_imbal = params.epsilon;

		if (str_preset == "Q") {
			args.MemMul_Pins = 100;
			args.MemMul_CellNet = 100;
			args.MemMul_General = 100;
		}

		int c, n, nconst, *cwghts, *nwghts, *xpins, *pins, *partvec, cut, *partweights;
		float* targetweights;
		n = hg.numHyperedges();
		c = hg.numNodes();
		nconst = 1;
		partvec = vec_partition.data();
		partweights = vec_partweights.data();
		cwghts = vec_cwghts.data();
		nwghts = vec_nwghts.data();
		pins = vec_pins.data();
		xpins = vec_xpins.data();
		if (params.use_target_weights) {
			targetweights = params.target_weights.data();
		}
		else {
			targetweights = NULL;
		}

		PaToH_Alloc(&args, c, n, nconst, cwghts, nwghts, xpins, pins);
		args.seed = seed;
		PaToH_Part(&args, c, n, nconst, 0, cwghts, nwghts, xpins, pins, targetweights, partvec, partweights, &cut);

		PaToH_Free();

		return vec_partition;
	}

};
