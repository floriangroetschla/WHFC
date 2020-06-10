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
        LawlerFlowHypergraph(size_t maxNumNodes, size_t maxNumHyperedges) : Base(maxNumNodes, maxNumHyperedges), vec_excess(numLawlerNodes(), 0), vec_label() {
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

        inline bool isNode(Node u) { return u < numNodes(); }

        inline void pushToEdgeIn(Node u, InHe& in_he, Flow f) {
            assert(f > 0);
            //assert(vec_excess[u] >= f); not true if u is source
            assert(f <= capacity(in_he.e) - flowSent(in_he.flow));
            //assert(flowReceived(in_he) <= 0);
            assert(isNode(u));

            Flow prevFlow = in_he.flow;

            vec_excess[u] -= f;
            vec_excess[edge_node_in(in_he.e)] += f;
            in_he.flow_in += flowSent(f);
        }

        inline void pushToEdgeOut(Node u, InHe& in_he, Flow f) {
            assert(f > 0);
            //assert(vec_excess[u] >= f); not true if u is source
            assert(f <= capacity(in_he.e) - flowSent(in_he.flow));
            //assert(flowReceived(in_he) > 0);
            assert(isNode(u));

            Flow prevFlow = in_he.flow;

            vec_excess[u] -= f;
            vec_excess[edge_node_out(in_he.e)] += f;
            in_he.flow_out += flowSent(f);
        }

        inline void pushFromEdgeInToNode(Node u, InHe& in_he, Flow f) {
            assert(f > 0);
            assert(flowSent(in_he) >= f);
            assert(isNode(u));

            Flow prevFlow = in_he.flow;

            vec_excess[u] += f;
            vec_excess[edge_node_in(in_he.e)] -= f;
            in_he.flow_in -= flowSent(f);
        }

        inline void pushFromEdgeOutToNode(Node u, InHe& in_he, Flow f) {
            assert(f > 0);
            //assert(flowReceived(in_he) >= 0);
            assert(isNode(u));

            Flow prevFlow = in_he.flow;

            vec_excess[u] += f;
            vec_excess[edge_node_out(in_he.e)] -= f;
            in_he.flow_out -= flowSent(f);
        }

        inline void pushFromEdgeInToEdgeOut(Node e_in, Node e_out, Flow f) {
            assert(f > 0);
            assert(f <= excess(e_in));
            assert(is_edge_out(e_out) && is_edge_in(e_in));

            vec_excess[e_in] -= f;
            vec_excess[e_out] += f;

            flow(edgeFromLawlerNode(e_in)) += f;
        }

        inline void pushFromEdgeOutToEdgeIn(Node e_out, Node e_in, Flow f) {
            assert(f > 0);
            assert(f <= excess(e_out));
            assert(f <= flow(edgeFromLawlerNode(e_in)));
            assert(is_edge_out(e_out) && is_edge_in(e_in));

            vec_excess[e_out] -= f;
            vec_excess[e_in] += f;

            flow(edgeFromLawlerNode(e_in)) -= f;
        }

        inline bool is_node(Node u) { return u < numNodes(); }
        inline bool is_edge_in(Node u) { return u >= numNodes() && u < numNodes() + numHyperedges(); }
        inline bool is_edge_out(Node u) { return u >= numNodes() + numHyperedges() && u < numLawlerNodes(); }

        inline Flow& excess_node(Node u) { assert(u < numNodes()); return vec_excess[u]; }
        inline Flow& excess_edge(Hyperedge e) { assert(e < numHyperedges()); return vec_excess[numNodes() + e]; }

        inline Node edge_node_in(Hyperedge e) { assert(e < numHyperedges()); return Node(numNodes() + e); }
        inline Node edge_node_out(Hyperedge e) { assert(e < numHyperedges()); return Node(numNodes() + numHyperedges() + e); }

        inline Hyperedge edgeFromLawlerNode(Node u) { assert(u >= numNodes()); return Hyperedge(u < numNodes() + numHyperedges() ? u - numNodes() : u - numNodes() - numHyperedges()); }

        inline Flow flowSent(const InHe& inc_u) const { return inc_u.flow_in; }
        inline Flow flowReceived(const InHe& inc_u) const { return inc_u.flow_out; }

        inline Flow flowSent(const Flow f) const { return f * sends_multiplier; }


        inline size_t numLawlerNodes() const {
            return numNodes() + 2 * numHyperedges();
        }

        void initialize_for_push_relabel() {
            vec_excess.clear();
            vec_excess.resize(numLawlerNodes(), 0);

            vec_label.clear();
            vec_label.resize(numLawlerNodes(), 0);
        }

        void printExcessAndLabel() {
            std::cout << "Excess: ";
            for (uint i = 0; i < vec_excess.size(); ++i) {
                std::cout << vec_excess[i] << " ";
            }
            std::cout << std::endl;

            std::cout << "Label: ";
            for (uint i = 0; i < vec_label.size(); ++i)  {
                std::cout << vec_label[i] << " ";
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
                out << e << " pincount = " << pinCount(e) << " w= " << capacity(e) << " pins (pin,flow_in,flow_out) = [";
                for (const Pin& u : pinsOf(e)) {
                    out << "(" << u.pin << "," << getInHe(u).flow_in << "," << getInHe(u).flow_out << ") ";
                }
                out << "]" << "\n";
            }
            out << std::flush;
        }

        void printHypergraph(std::ostream& out) {
            printNodes(out);
            printHyperedges(out);
        }

    private:
        std::vector<Flow> vec_excess;
        std::vector<size_t> vec_label;
    };
}