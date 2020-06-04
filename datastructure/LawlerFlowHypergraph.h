#pragma once

#include <atomic>
#include <tbb/parallel_for.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for_each.h>
#include "flow_hypergraph_builder.h"
#include "../recursive_bisection/mock_builder.h"

namespace whfc {
    class LawlerFlowHypergraph : public FlowHypergraphBuilder {
    public:
        using Base = FlowHypergraphBuilder;

        LawlerFlowHypergraph() : Base() {
            clear();
        }

        /*
        explicit LawlerFlowHypergraph(size_t nNodes) {
            reinitialize(nNodes);
        }*/

        //use to get rid of any allocations
        LawlerFlowHypergraph(size_t maxNumNodes, size_t maxNumHyperedges) : Base(maxNumNodes, maxNumHyperedges), vec_excess(maxNumNodes + 2 * maxNumHyperedges, 0), vec_label() {
            //don't do clean-up here yet, so that we can use the numbers for allocating the remaining datastructures
        }

        void clear() {
            Base::clear();

            vec_excess.clear();
            vec_label.clear();
        }

        /*
        void reinitialize(size_t numNodes) {
            clear();
            nodes.resize(numNodes + 1);
        }*/


        void shrink_to_fit() {
            Base::shrink_to_fit();

            vec_excess.shrink_to_fit();
            vec_label.shrink_to_fit();
        }

        inline size_t& label(Node u) {
            assert(u < numLawlerNodes());
            return vec_label[u];
        }

        inline Flow& excess(Node u) { assert(u < numLawlerNodes()); return vec_excess[u]; }

        inline void pushToEdge(Node u, InHe& in_he, Flow flow) {
            assert(vec_excess[u] >= flow);
            assert(flow <= capacity(in_he.e) - flowSent(in_he.flow));

            // Push back to e_out
            Flow to_e_out = std::min(absoluteFlowReceived(in_he), flow);
            if (to_e_out > 0) {
                vec_excess[edge_out_node(in_he.e)] += to_e_out;
                in_he.flow += flowSent(to_e_out);
                flow -= to_e_out;
            }

            // Push the rest to e_in
            if (flow > 0) {
                vec_excess[edge_in_node(in_he.e)] += flow;
                in_he.flow += flowSent(flow);
            }

            vec_excess[u] -= flow;
        }

        inline Flow& excess_node(Node u) { assert(u < numNodes()); return vec_excess[u]; }
        inline Flow& excess_edge_in(Hyperedge e) { assert(e < numHyperedges()); return vec_excess[numNodes() + e]; }
        inline Flow& excess_edge_out(Hyperedge e) { assert(e < numHyperedges()); return vec_excess[numNodes() + numHyperedges() + e]; }

        inline Node edge_in_node(Hyperedge e) { assert(e < numHyperedges()); return Node(numNodes() + e); }
        inline Node edge_out_node(Hyperedge e) { assert(e < numHyperedges()); return Node(numNodes() + numHyperedges() + e); }


        inline size_t numLawlerNodes() const {
            return numNodes() + numHyperedges() * 2;
        }

        void initialize_for_push_relabel() {
            vec_excess.clear();
            vec_excess.resize(numLawlerNodes(), 0);

            vec_label.clear();
            vec_label.resize(numLawlerNodes(), 0);
        }

    private:
        std::vector<Flow> vec_excess;
        std::vector<size_t> vec_label;
    };
}