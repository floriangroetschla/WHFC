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
        LawlerFlowHypergraph(size_t maxNumNodes, size_t maxNumHyperedges) : Base(maxNumNodes, maxNumHyperedges), vec_excess(numLawlerNodes(), 0), vec_label(), in_flow(), out_flow() {
            //don't do clean-up here yet, so that we can use the numbers for allocating the remaining datastructures
        }

        void clear() {
            Base::clear();

            vec_excess.clear();
            vec_label.clear();
            in_flow.clear();
            out_flow.clear();
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
            in_flow.shrink_to_fit();
            out_flow.shrink_to_fit();
        }

        inline size_t& label(Node u) {
            assert(u < numLawlerNodes());
            return vec_label[u];
        }

        inline Flow& excess(Node u) { assert(u < numLawlerNodes()); return vec_excess[u]; }

        inline bool isNode(Node u) { return u < numNodes(); }

        inline Flow& flowIn(InHe in_he) { return in_flow[pins[in_he.pin_iter].he_inc_iter]; }
        inline Flow& flowOut(InHe in_he) { return out_flow[pins[in_he.pin_iter].he_inc_iter]; }

        inline Flow flowIn(InHe in_he) const { return in_flow[pins[in_he.pin_iter].he_inc_iter]; }
        inline Flow flowOut(InHe in_he) const { return out_flow[pins[in_he.pin_iter].he_inc_iter]; }

        inline Flow& flowIn(InHeIndex inc) { return in_flow[inc]; }
        inline Flow& flowOut(InHeIndex inc) { return out_flow[inc]; }
        inline Flow flowIn(InHeIndex inc) const { return in_flow[inc]; }
        inline Flow flowOut(InHeIndex inc) const { return out_flow[inc]; }

        inline void push_node_to_edgeIn(const Node u, const InHeIndex inc, const Flow f) {
            assert(f > 0);

            const InHe& in_he = FlowHypergraph::getInHe(inc);

            vec_excess[u] -= f;
            vec_excess[edge_node_in(in_he.e)] += f;
            flowIn(inc) += f;
        }

        inline void push_node_to_edgeOut(const Node u, const InHeIndex inc, const Flow f) {
            assert(f > 0);
            assert(isNode(u));

            const InHe& in_he = FlowHypergraph::getInHe(inc);

            vec_excess[u] -= f;
            vec_excess[edge_node_out(in_he.e)] += f;
            flowOut(inc) -= f;
        }

        inline void push_edgeIn_to_node(const Node u, const InHeIndex inc, const Flow f) {
            assert(f > 0);
            assert(flowSent(inc) >= f);
            assert(isNode(u));

            const InHe& in_he = FlowHypergraph::getInHe(inc);

            vec_excess[u] += f;
            vec_excess[edge_node_in(in_he.e)] -= f;
            flowIn(inc) -= f;
        }

        inline void push_edgeOut_to_node(const Node u, const InHeIndex inc, const Flow f) {
            assert(f > 0);
            assert(isNode(u));

            const InHe& in_he = FlowHypergraph::getInHe(inc);

            vec_excess[u] += f;
            vec_excess[edge_node_out(in_he.e)] -= f;
            flowOut(inc) += f;
        }

        inline void push_edgeIn_to_edgeOut(const Node e_in, const Node e_out, const Flow f) {
            assert(f > 0);
            assert(f <= excess(e_in));
            assert(is_edge_out(e_out) && is_edge_in(e_in));

            vec_excess[e_in] -= f;
            vec_excess[e_out] += f;

            flow(edgeFromLawlerNode(e_in)) += f;
        }

        inline void push_edgeOut_to_edgeIn(const Node e_out, const Node e_in, const Flow f) {
            assert(f > 0);
            assert(f <= excess(e_out));
            assert(f <= flow(edgeFromLawlerNode(e_in)));
            assert(is_edge_out(e_out) && is_edge_in(e_in));

            vec_excess[e_out] -= f;
            vec_excess[e_in] += f;

            flow(edgeFromLawlerNode(e_in)) -= f;
        }

        void sortPins() {
            for (Hyperedge e : hyperedgeIDs()) {
                PinIndex begin = beginIndexPins(e);
                PinIndex end = PinIndex(endIndexPins(e) - 1);

                for (PinIndex i = begin; i <= end;) {
                    Pin& p = pins[i];
                    const InHeIndex inc = p.he_inc_iter;
                    InHe& inc_he = getInHe(p);

                    if (flowIn(inc) && flowOut(inc)) {
                        Flow loopedFlow = std::min(flowIn(inc), flowOut(inc));
                        flowIn(inc) -= loopedFlow;
                        flowOut(inc) -= loopedFlow;

                        flow(inc_he.e) -= loopedFlow;
                    }

                    assert(!flowIn(inc) || !flowOut(inc));
                    inc_he.flow = flowSent(flowIn(inc) - flowOut(inc));

                    if (inc_he.flow > 0) {
                        InHe& inc_begin = getInHe(pins[begin]);
                        std::swap(inc_he.pin_iter, inc_begin.pin_iter);
                        std::swap(pins[i++], pins[begin++]);
                    } else if (inc_he.flow < 0) {
                        InHe& inc_end = getInHe(pins[end]);
                        std::swap(inc_he.pin_iter, inc_end.pin_iter);
                        std::swap(pins[i], pins[end--]);
                    } else {
                        i++;
                    }
                }

                if (forward) {
                    pins_sending_flow[e] = PinIndexRange(beginIndexPins(e), begin);
                    pins_receiving_flow[e] = PinIndexRange(PinIndex(end + 1), endIndexPins(e));
                } else {
                    pins_receiving_flow[e] = PinIndexRange(beginIndexPins(e), begin);
                    pins_sending_flow[e] = PinIndexRange(PinIndex(end + 1), endIndexPins(e));
                }

            }
        }

        inline bool is_node(Node u) { return u < numNodes(); }
        inline bool is_edge_in(Node u) { return u >= numNodes() && u < numNodes() + numHyperedges(); }
        inline bool is_edge_out(Node u) { return u >= numNodes() + numHyperedges() && u < numLawlerNodes(); }

        inline Flow& excess_node(Node u) { assert(u < numNodes()); return vec_excess[u]; }
        inline Flow& excess_edge(Hyperedge e) { assert(e < numHyperedges()); return vec_excess[numNodes() + e]; }

        inline Node edge_node_in(Hyperedge e) { assert(e < numHyperedges()); return Node(numNodes() + e); }
        inline Node edge_node_out(Hyperedge e) { assert(e < numHyperedges()); return Node(numNodes() + numHyperedges() + e); }

        inline Hyperedge edgeFromLawlerNode(Node u) { assert(u >= numNodes()); return Hyperedge(u < numNodes() + numHyperedges() ? u - numNodes() : u - numNodes() - numHyperedges()); }

        inline Flow flowSent(const InHe& inc_u) const { return flowIn(inc_u); }
        inline Flow flowReceived(const InHe& inc_u) const { return flowOut(inc_u); }

        inline Flow flowSent(const InHeIndex inc) const { return flowIn(inc); }
        inline Flow flowReceived(const InHeIndex inc) const { return flowOut(inc); }

        inline Flow flowSent(const Flow f) const { return f * sends_multiplier; }

        void alignViewDirection() {
            if (forward != FlowHypergraph::forwardView()) {
                std::swap(in_flow, out_flow);
                forward = FlowHypergraph::forwardView();
            }
        }


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
                out << e << " pincount = " << pinCount(e) << " w= " << capacity(e) << " pins (pin,flow,flow_in,flow_out) = [";
                for (const Pin& u : pinsOf(e)) {
                    out << "(" << u.pin << "," << getInHe(u).flow << "," << flowIn(getInHe(u)) << "," << flowOut(getInHe(u)) << ") ";
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
            in_flow = std::vector<Flow>(numPins(), 0);
            out_flow = std::vector<Flow>(numPins(), 0);
        }

    private:
        std::vector<Flow> vec_excess;
        std::vector<size_t> vec_label;

        std::vector<Flow> in_flow;
        std::vector<Flow> out_flow;

        bool forward = true;
    };
}