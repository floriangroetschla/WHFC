#pragma once

#include <atomic>
#include <tbb/parallel_for.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for_each.h>
#include "flow_hypergraph_builder.h"
#include "../recursive_bisection/mock_builder.h"
#include "../recursive_bisection/timestamp_set.hpp"

namespace whfc {
    class LawlerFlowHypergraph : public FlowHypergraphBuilder {
    public:
        using Base = FlowHypergraphBuilder;

        LawlerFlowHypergraph() : Base(), vec_excess(), vec_label_this_iteration(), vec_excess_change(), vec_label_next_iteration() {
            clear();
        }

        LawlerFlowHypergraph(size_t maxNumNodes, size_t maxNumHyperedges) : Base(maxNumNodes, maxNumHyperedges),
            vec_excess(numLawlerNodes()), vec_label_this_iteration(numLawlerNodes()), vec_excess_change(numLawlerNodes()), vec_label_next_iteration(numLawlerNodes()) {

        }

        void clear() {
            Base::clear();
        }

        void shrink_to_fit() {
            Base::shrink_to_fit();
        }

        inline size_t& label(Node u) { assert(u < numLawlerNodes()); return vec_label_this_iteration[u]; }
        inline size_t& label_next_iteration(Node u) { assert(u < numLawlerNodes()); return vec_label_next_iteration[u]; }

        inline Flow& excess(Node u) { assert(u < numLawlerNodes()); return vec_excess[u]; }
        inline Flow& excess_change(Node u) { assert(u < numLawlerNodes()); return vec_excess_change[u]; }

        inline bool isNode(Node u) { return u < numNodes(); }

        inline void push_node_to_edgeIn(const Node u, const InHeIndex inc, const Flow f) {
            assert(f > 0);

            const InHe& in_he = FlowHypergraph::getInHe(inc);

            //assert(vec_excess[u] >= f);
            vec_excess[u] -= f;
            __atomic_add_fetch(&vec_excess_change[edge_node_in(in_he.e)], f, __ATOMIC_ACQ_REL);
            getPinIn(in_he).flow += f;
        }

        inline void push_node_to_edgeOut(const Node u, const InHeIndex inc, const Flow f) {
            assert(f > 0);
            assert(isNode(u));

            const InHe& in_he = FlowHypergraph::getInHe(inc);

            assert(vec_excess[u] >= f);
            vec_excess[u] -= f;
            __atomic_add_fetch(&vec_excess_change[edge_node_out(in_he.e)], f, __ATOMIC_ACQ_REL);
            getPinOut(in_he).flow -= f;
        }

        inline void push_edgeIn_to_node(const Node u, const InHeIndex inc, const Flow f) {
            assert(f > 0);
            assert(getPinIn(getInHe(inc)).flow >= f);
            assert(isNode(u));

            const InHe& in_he = FlowHypergraph::getInHe(inc);

            assert(vec_excess[edge_node_in(in_he.e)] >= f);
            vec_excess[edge_node_in(in_he.e)] -= f;
            __atomic_add_fetch(&vec_excess_change[u], f, __ATOMIC_ACQ_REL);
            getPinIn(in_he).flow -= f;
        }

        inline void push_edgeOut_to_node(const Node u, const InHeIndex inc, const Flow f) {
            assert(f > 0);
            assert(isNode(u));

            const InHe& in_he = FlowHypergraph::getInHe(inc);

            assert(vec_excess[edge_node_out(in_he.e)] >= f);
            vec_excess[edge_node_out(in_he.e)] -= f;
            __atomic_add_fetch(&vec_excess_change[u], f, __ATOMIC_ACQ_REL);
            getPinOut(in_he).flow += f;
        }

        inline void push_edgeIn_to_edgeOut(const Node e_in, const Node e_out, const Flow f) {
            assert(f > 0);
            assert(f <= excess(e_in));
            assert(is_edge_out(e_out) && is_edge_in(e_in));

            assert(vec_excess[e_in] >= f);
            vec_excess[e_in] -= f;
            __atomic_add_fetch(&vec_excess_change[e_out], f, __ATOMIC_ACQ_REL);

            flow(edgeFromLawlerNode(e_in)) += f;
        }

        inline void push_edgeOut_to_edgeIn(const Node e_out, const Node e_in, const Flow f) {
            assert(f > 0);
            assert(f <= excess(e_out));
            assert(f <= flow(edgeFromLawlerNode(e_in)));
            assert(is_edge_out(e_out) && is_edge_in(e_in));

            assert(vec_excess[e_out] >= f);
            vec_excess[e_out] -= f;
            __atomic_add_fetch(&vec_excess_change[e_in], f, __ATOMIC_ACQ_REL);

            flow(edgeFromLawlerNode(e_in)) -= f;
        }


        inline bool is_node(Node u) { return u < numNodes(); }
        inline bool is_edge_in(Node u) { return u >= numNodes() && u < numNodes() + numHyperedges(); }
        inline bool is_edge_out(Node u) { return u >= numNodes() + numHyperedges() && u < numLawlerNodes(); }

        inline Node edge_node_in(Hyperedge e) { assert(e < numHyperedges()); return Node(numNodes() + e); }
        inline Node edge_node_out(Hyperedge e) { assert(e < numHyperedges()); return Node(numNodes() + numHyperedges() + e); }

        inline Hyperedge edgeFromLawlerNode(Node u) { assert(u >= numNodes()); return Hyperedge(u < numNodes() + numHyperedges() ? u - numNodes() : u - numNodes() - numHyperedges()); }

        inline size_t maxNumLawlerNodes() { return maxNumNodes + 2 * maxNumHyperedges; }

        void alignViewDirection() {
            if (forward != FlowHypergraph::forwardView()) {
                forward = FlowHypergraph::forwardView();
            }
        }


        inline size_t numLawlerNodes() const {
            return numNodes() + 2 * numHyperedges();
        }


        void initialize_for_push_relabel() {
            std::fill(vec_excess.begin(), vec_excess.end(), 0);
            std::fill(vec_excess_change.begin(), vec_excess_change.end(), 0);
        }

        void equalizeLabels(size_t n) {
            std::fill(vec_label_this_iteration.begin(), vec_label_this_iteration.end(), n);
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

        void printExcessAndLabel() {
            std::cout << "Excess: ";
            for (uint i = 0; i < vec_excess.size(); ++i) {
                std::cout << vec_excess[i] << " ";
            }
            std::cout << std::endl;

            std::cout << "Label: ";
            for (uint i = 0; i < vec_label_this_iteration.size(); ++i)  {
                std::cout << vec_label_this_iteration[i] << " ";
            }
            std::cout << std::endl;
        }

        void printHyperedge(Hyperedge e) {
            std::cout << e << " pincount = " << pinCount(e) << " w= " << capacity(e) << " pins (pin,flow,excess,label) = [";
            for (const Pin& u : pinsOf(e)) {
                assert(pin_is_categorized_correctly(getInHe(u)));
                std::cout << "(" << u.pin << "," << getInHe(u).flow << "," << excess(u.pin) << "," << label(u.pin) << ") ";
            }
            std::cout << "]" << "\n";
        }

        void printNodes(std::ostream& out) {
            out << "---Nodes---\n";
            for (const Node u : nodeIDs()) {
                out << u << " deg = " << degree(u) << " w= " << nodeWeight(u) << " inc_hes = [";
                for (const InHe e : hyperedgesOf(u))
                    out << e.e << " ";
                out << "]" << "\n";
            }
            out << std::flush;
        }

        void printHyperedges(std::ostream& out) {
            out << "---Hyperedges---\n";
            for (const Hyperedge e: hyperedgeIDs()) {
                out << e << " pincount = " << pinCount(e) << " w= " << capacity(e) << " pins (pin,flow,flow_in,flow_out) = [";
                for (const Pin& u : pinsOf(e)) {
                    out << "(" << u.pin << "," << getInHe(u).flow << "," << getPinIn(getInHe(u)).flow << "," << getPinOut(getInHe(u)).flow << ") ";
                }
                out << "]" << "\n";
            }
            out << std::flush;
        }

        void printHypergraph(std::ostream& out) {
            printNodes(out);
            printHyperedges(out);
        }

        void finalize() {
            Base::finalize();
            vec_excess.resize(numLawlerNodes());
            vec_label_this_iteration.resize(numLawlerNodes());
            vec_excess_change.resize(numLawlerNodes());
            vec_label_next_iteration.resize(numLawlerNodes());
            //std::fill(vec_excess_this_iteration.begin(), vec_excess_this_iteration.end(), 0);
            std::fill(vec_label_this_iteration.begin(), vec_label_this_iteration.end(), 0);
            //std::fill(vec_excess_next_iteration.begin(), vec_excess_next_iteration.end(), 0);
            std::fill(vec_label_next_iteration.begin(), vec_label_next_iteration.end(), 0);
        }

    private:
        std::vector<Flow> vec_excess;
        std::vector<size_t> vec_label_this_iteration;

        std::vector<Flow> vec_excess_change;
        std::vector<size_t> vec_label_next_iteration;

        bool forward = true;
    };
}