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

        /*
        explicit LawlerFlowHypergraph(size_t nNodes) {
            reinitialize(nNodes);
        }*/

        LawlerFlowHypergraph(size_t maxNumNodes, size_t maxNumHyperedges) : Base(maxNumNodes, maxNumHyperedges),
            vec_excess(numLawlerNodes()), vec_label_this_iteration(numLawlerNodes()), vec_excess_change(numLawlerNodes()), vec_label_next_iteration(numLawlerNodes()), in_flow(), out_flow() {

        }

        void clear() {
            Base::clear();

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

            in_flow.shrink_to_fit();
            out_flow.shrink_to_fit();
        }

        inline size_t& label(Node u) { assert(u < numLawlerNodes()); return vec_label_this_iteration[u]; }
        inline size_t& label_next_iteration(Node u) { assert(u < numLawlerNodes()); return vec_label_next_iteration[u]; }

        inline Flow& excess(Node u) { assert(u < numLawlerNodes()); return vec_excess[u]; }
        inline Flow& excess_change(Node u) { assert(u < numLawlerNodes()); return vec_excess_change[u]; }

        inline bool isNode(Node u) { return u < numNodes(); }

        inline Flow& flowIn(InHe in_he) { return in_flow[pins[in_he.pin_iter].he_inc_iter]; }
        inline Flow& flowOut(InHe in_he) { return out_flow[pins[in_he.pin_iter].he_inc_iter]; }

        inline Flow flowIn(InHe in_he) const { return in_flow[pins[in_he.pin_iter].he_inc_iter]; }
        inline Flow flowOut(InHe in_he) const { return out_flow[pins[in_he.pin_iter].he_inc_iter]; }

        inline Flow& flowIn(InHeIndex inc) { return in_flow[inc]; }
        inline Flow& flowOut(InHeIndex inc) { return out_flow[inc]; }
        inline Flow flowIn(InHeIndex inc) const { return in_flow[inc]; }
        inline Flow flowOut(InHeIndex inc) const { return out_flow[inc]; }



        void routeFlow(const InHeIndex in_he_u, const InHeIndex in_he_v, Flow flow_delta) {
            InHe& inc_u = getInHe(in_he_u);
            InHe& inc_v = getInHe(in_he_v);

            Flow flow_delta_on_v_eOut_eIn_u = std::min({ absoluteFlowSent(inc_v), absoluteFlowReceived(inc_u), flow_delta });
            inc_u.flow += flowSent(flow_delta_on_v_eOut_eIn_u);
            inc_v.flow -= flowSent(flow_delta_on_v_eOut_eIn_u);
            flow_delta -= flow_delta_on_v_eOut_eIn_u;

            Flow flow_delta_on_v_eOut_u = std::min(flow_delta, absoluteFlowReceived(inc_u));
            inc_u.flow += flowSent(flow_delta_on_v_eOut_u);
            inc_v.flow -= flowSent(flow_delta_on_v_eOut_u);
            //does not influence flow(e)
            flow_delta -= flow_delta_on_v_eOut_u;

            Flow flow_delta_on_v_eIn_u = std::min(flow_delta, absoluteFlowSent(inc_v));
            inc_u.flow += flowSent(flow_delta_on_v_eIn_u);
            inc_v.flow -= flowSent(flow_delta_on_v_eIn_u);
            //does not influence flow(e)
            flow_delta -= flow_delta_on_v_eIn_u;

            Flow flow_delta_on_v_eIn_eOut_u = flow_delta;
            inc_u.flow += flowSent(flow_delta_on_v_eIn_eOut_u);
            inc_v.flow -= flowSent(flow_delta_on_v_eIn_eOut_u);
        }

        inline void push_node_to_edgeIn(const Node u, const InHeIndex inc, const Flow f) {
            assert(f > 0);

            const InHe& in_he = FlowHypergraph::getInHe(inc);

            vec_excess[u] -= f;
            __atomic_add_fetch(&vec_excess_change[edge_node_in(in_he.e)], f, __ATOMIC_ACQ_REL);
            flowIn(inc) += f;
        }

        inline void push_node_to_edgeOut(const Node u, const InHeIndex inc, const Flow f) {
            assert(f > 0);
            assert(isNode(u));

            const InHe& in_he = FlowHypergraph::getInHe(inc);

            vec_excess[u] -= f;
            __atomic_add_fetch(&vec_excess_change[edge_node_out(in_he.e)], f, __ATOMIC_ACQ_REL);
            flowOut(inc) -= f;
        }

        inline void push_edgeIn_to_node(const Node u, const InHeIndex inc, const Flow f) {
            assert(f > 0);
            assert(flowSent(inc) >= f);
            assert(isNode(u));

            const InHe& in_he = FlowHypergraph::getInHe(inc);

            vec_excess[edge_node_in(in_he.e)] -= f;
            __atomic_add_fetch(&vec_excess_change[u], f, __ATOMIC_ACQ_REL);
            flowIn(inc) -= f;
        }

        inline void push_edgeOut_to_node(const Node u, const InHeIndex inc, const Flow f) {
            assert(f > 0);
            assert(isNode(u));

            const InHe& in_he = FlowHypergraph::getInHe(inc);

            vec_excess[edge_node_out(in_he.e)] -= f;
            __atomic_add_fetch(&vec_excess_change[u], f, __ATOMIC_ACQ_REL);
            flowOut(inc) += f;
        }

        inline void push_edgeIn_to_edgeOut(const Node e_in, const Node e_out, const Flow f) {
            assert(f > 0);
            assert(f <= excess(e_in));
            assert(is_edge_out(e_out) && is_edge_in(e_in));

            vec_excess[e_in] -= f;
            __atomic_add_fetch(&vec_excess_change[e_out], f, __ATOMIC_ACQ_REL);

            flow(edgeFromLawlerNode(e_in)) += f;
        }

        inline void push_edgeOut_to_edgeIn(const Node e_out, const Node e_in, const Flow f) {
            assert(f > 0);
            assert(f <= excess(e_out));
            assert(f <= flow(edgeFromLawlerNode(e_in)));
            assert(is_edge_out(e_out) && is_edge_in(e_in));

            vec_excess[e_out] -= f;
            __atomic_add_fetch(&vec_excess_change[e_in], f, __ATOMIC_ACQ_REL);

            flow(edgeFromLawlerNode(e_in)) -= f;
        }

        void writeBackFlow() {
            for (Hyperedge e : hyperedgeIDs()) {
                Flow flow_on_edge = 0;
                for (Pin& p : pinsOf(e)) {
                    const InHeIndex inc = p.he_inc_iter;
                    InHe& inc_he = getInHe(p);
                    flow_on_edge += absoluteFlowSent(inc_he);

                    flowIn(inc) = absoluteFlowSent(inc_he);
                    flowOut(inc) = absoluteFlowReceived(inc_he);
                }
                flow(e) = flow_on_edge;
            }
        }

        void buildResidualNetwork() {
            for (Hyperedge e : hyperedgeIDs()) {
                Flow previous_flow_on_edge = 0;
                for (Pin& p : pinsOf(e)) {
                    const InHeIndex inc = p.he_inc_iter;
                    InHe& inc_he = getInHe(p);
                    previous_flow_on_edge += absoluteFlowSent(inc_he);

                    flowIn(inc) = flowSent(inc) - absoluteFlowSent(inc_he);
                    flowOut(inc) = flowReceived(inc) - absoluteFlowReceived(inc_he);
                }
                flow(e) -= previous_flow_on_edge;
            }
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

        //inline Flow& excess_node(Node u) { assert(u < numNodes()); return vec_excess[u]; }
        //inline Flow& excess_edge(Hyperedge e) { assert(e < numHyperedges()); return vec_excess[numNodes() + e]; }

        inline Node edge_node_in(Hyperedge e) { assert(e < numHyperedges()); return Node(numNodes() + e); }
        inline Node edge_node_out(Hyperedge e) { assert(e < numHyperedges()); return Node(numNodes() + numHyperedges() + e); }

        inline Hyperedge edgeFromLawlerNode(Node u) { assert(u >= numNodes()); return Hyperedge(u < numNodes() + numHyperedges() ? u - numNodes() : u - numNodes() - numHyperedges()); }

        inline Flow flowSent(const InHe& inc_u) const { return flowIn(inc_u); }
        inline Flow flowReceived(const InHe& inc_u) const { return flowOut(inc_u); }

        inline Flow flowSent(const InHeIndex inc) const { return flowIn(inc); }
        inline Flow flowReceived(const InHeIndex inc) const { return flowOut(inc); }

        inline Flow flowSent(const Flow f) const { return f * sends_multiplier; }

        inline size_t maxNumLawlerNodes() { return maxNumNodes + 2 * maxNumHyperedges; }

        void alignViewDirection() {
            if (forward != FlowHypergraph::forwardView()) {
                std::swap(in_flow, out_flow);
                forward = FlowHypergraph::forwardView();
            }
        }


        inline size_t numLawlerNodes() const {
            return numNodes() + 2 * numHyperedges();
        }

        void copy_excess_vectors() {

        }

        void initialize_for_push_relabel() {
            std::fill(vec_excess.begin(), vec_excess.end(), 0);
            std::fill(vec_excess_change.begin(), vec_excess_change.end(), 0);
        }

        void equalizeLabels(size_t n) {
            std::fill(vec_label_this_iteration.begin(), vec_label_this_iteration.end(), n);
            std::fill(vec_label_next_iteration.begin(), vec_label_next_iteration.end(), n);
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
            in_flow.resize(numPins(), 0);
            out_flow.resize(numPins(), 0);
            std::fill(in_flow.begin(), in_flow.end(), 0);
            std::fill(out_flow.begin(), out_flow.end(), 0);
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

        std::vector<Flow> in_flow;
        std::vector<Flow> out_flow;

        bool forward = true;
    };
}