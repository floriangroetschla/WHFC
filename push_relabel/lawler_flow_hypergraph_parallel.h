#pragma once

#include <atomic>
#include <tbb/parallel_for.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for_each.h>
#include "flow_hypergraph_builder.h"
#include "../recursive_bisection/mock_builder.h"
#include "../recursive_bisection/timestamp_set.hpp"
#include "lawler_flow_hypergraph.h"

namespace whfc_pr {
    class LawlerFlowHypergraphParallel : public LawlerFlowHypergraph {
    public:
        using Base = LawlerFlowHypergraph;

        LawlerFlowHypergraphParallel() : Base(), vec_excess_change(), vec_label_next_iteration() {
            clear();
        }

        LawlerFlowHypergraphParallel(size_t maxNumNodes, size_t maxNumHyperedges) : Base(maxNumNodes, maxNumHyperedges),
            vec_excess_change(numLawlerNodes()), vec_label_next_iteration(numLawlerNodes()) {

        }

        void clear() {
            Base::clear();
        }

        void shrink_to_fit() {
            Base::shrink_to_fit();
        }

        inline size_t& label(Node u) { assert(u < numLawlerNodes()); return vec_label[u]; }
        inline size_t& label_next_iteration(Node u) { assert(u < numLawlerNodes()); return vec_label_next_iteration[u]; }

        inline Flow& excess(Node u) { assert(u < numLawlerNodes()); return vec_excess[u]; }
        inline Flow& excess_change(Node u) { assert(u < numLawlerNodes()); return vec_excess_change[u]; }

        inline void push_node_to_edgeIn(const Node u, const InHeIndex inc, const Flow f) {
            assert(f > 0);

            const InHe& in_he = FlowHypergraph::getInHe(inc);

            getPinIn(in_he).flow += f;
            vec_excess[u] -= f;
            __atomic_add_fetch(&vec_excess_change[edge_node_in(in_he.e)], f, __ATOMIC_ACQ_REL);
        }

        inline void push_node_to_edgeOut(const Node u, const InHeIndex inc, const Flow f) {
            assert(f > 0);
            assert(isNode(u));

            const InHe& in_he = FlowHypergraph::getInHe(inc);

            getPinOut(in_he).flow -= f;
            vec_excess[u] -= f;
            __atomic_add_fetch(&vec_excess_change[edge_node_out(in_he.e)], f, __ATOMIC_ACQ_REL);
        }

        inline void push_edgeIn_to_node(const Node u, const InHeIndex inc, const Flow f) {
            assert(f > 0);
            assert(getPinIn(getInHe(inc)).flow >= f);
            assert(isNode(u));

            const InHe& in_he = FlowHypergraph::getInHe(inc);

            getPinIn(in_he).flow -= f;
            assert(vec_excess[edge_node_in(in_he.e)] >= f);
            vec_excess[edge_node_in(in_he.e)] -= f;
            __atomic_add_fetch(&vec_excess_change[u], f, __ATOMIC_ACQ_REL);
        }

        inline void push_edgeOut_to_node(const Node u, const InHeIndex inc, const Flow f) {
            assert(f > 0);
            assert(isNode(u));

            const InHe& in_he = FlowHypergraph::getInHe(inc);

            getPinOut(in_he).flow += f;
            assert(vec_excess[edge_node_out(in_he.e)] >= f);
            vec_excess[edge_node_out(in_he.e)] -= f;
            __atomic_add_fetch(&vec_excess_change[u], f, __ATOMIC_ACQ_REL);
        }

        inline void push_edgeIn_to_edgeOut(const Node e_in, const Node e_out, const Flow f) {
            assert(f > 0);
            assert(f <= excess(e_in));
            assert(is_edge_out(e_out) && is_edge_in(e_in));

            flow(edgeFromLawlerNode(e_in)) += f;
            assert(vec_excess[e_in] >= f);
            vec_excess[e_in] -= f;
            __atomic_add_fetch(&vec_excess_change[e_out], f, __ATOMIC_ACQ_REL);
        }

        inline void push_edgeOut_to_edgeIn(const Node e_out, const Node e_in, const Flow f) {
            assert(f > 0);
            assert(f <= excess(e_out));
            assert(f <= flow(edgeFromLawlerNode(e_in)));
            assert(is_edge_out(e_out) && is_edge_in(e_in));

            flow(edgeFromLawlerNode(e_in)) -= f;
            assert(vec_excess[e_out] >= f);
            vec_excess[e_out] -= f;
            __atomic_add_fetch(&vec_excess_change[e_in], f, __ATOMIC_ACQ_REL);
        }


        void initialize_for_push_relabel() {
            std::fill(vec_excess.begin(), vec_excess.end(), 0);
            std::fill(vec_excess_change.begin(), vec_excess_change.end(), 0);
        }

        void equalizeLabels(size_t n) {
            std::fill(vec_label.begin(), vec_label.end(), n);
            std::fill(vec_label_next_iteration.begin(), vec_label_next_iteration.end(), n);
        }

        void writeBackFlow() {
            for (Hyperedge e : hyperedgeIDs()) {
                Flow flow_on_edge = 0;
                for (Pin& p : pinsInOf(e)) {
                    const InHeIndex inc = p.he_inc_iter;
                    InHe& inc_he = getInHe(p);
                    flow_on_edge += absoluteFlowSent(inc_he);

                    getPinIn(inc_he).flow = absoluteFlowSent(inc_he);
                    getPinOut(inc_he).flow = absoluteFlowReceived(inc_he);
                }
                flow(e) = flow_on_edge;
            }
        }

        void buildResidualNetwork() {
            for (Hyperedge e : hyperedgeIDs()) {
                Flow previous_flow_on_edge = 0;
                for (Pin& p : pinsInOf(e)) {
                    const InHeIndex inc = p.he_inc_iter;
                    InHe& inc_he = getInHe(p);
                    previous_flow_on_edge += absoluteFlowSent(inc_he);

                    getPinIn(inc_he).flow -= absoluteFlowSent(inc_he);
                    getPinOut(inc_he).flow -= absoluteFlowReceived(inc_he);
                }
                flow(e) -= previous_flow_on_edge;
            }
        }

        bool excess_sums_to_zero() {
            Flow totalExcess = 0;
            for (size_t i = 0; i < vec_excess.size(); ++i) {
                totalExcess += vec_excess[i];
            }
            return totalExcess == 0;
        }


        void finalize() {
            Base::finalize();
            vec_excess.resize(numLawlerNodes());
            vec_label.resize(numLawlerNodes());
            vec_excess_change.resize(numLawlerNodes());
            vec_label_next_iteration.resize(numLawlerNodes());
            std::fill(vec_label.begin(), vec_label.end(), 0);
            std::fill(vec_label_next_iteration.begin(), vec_label_next_iteration.end(), 0);
        }

    private:
        std::vector<Flow> vec_excess_change;
        std::vector<size_t> vec_label_next_iteration;
    };
}