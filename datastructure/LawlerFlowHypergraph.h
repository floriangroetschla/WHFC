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

        inline void pushToEdge(Node u, InHe& in_he, Flow f) {
            assert(f > 0);
            //assert(vec_excess[u] >= f); not true if u is source
            assert(f <= capacity(in_he.e) - flowSent(in_he.flow));

            Flow prevFlow = in_he.flow;

            if (flowSent(in_he) < 0) {
                const Flow flow_pushed_back = std::max(flowSent(in_he), flowSent(-f));
                flow(in_he.e) += flow_pushed_back;
            }

            vec_excess[u] -= f;
            vec_excess[edge_node(in_he.e)] += f;
            in_he.flow += flowSent(f);

            if (flowReceived(prevFlow) > 0 && flowSent(in_he.flow) >= 0)	//u previously received flow and now either has none, or sends flow.
                removePinFromFlowPins(in_he, true);
            if (flowSent(in_he.flow) > 0 && flowSent(prevFlow) <= 0) //u now sends flow and did not previously, thus must be inserted into pins_sending_flow
                insertPinIntoFlowPins(in_he, false);
        }

        inline void pushToNode(Node u, InHe& in_he, Flow f) {
            assert(f > 0);

            Flow prevFlow = in_he.flow;

            if (flowSent(in_he) > 0) {
                const Flow flow_pushed = std::max(0, flowSent(f) - flowSent(in_he));
                flow(in_he.e) += flow_pushed;
            } else {
                const Flow flow_pushed = f;
                flow(in_he.e) += flow_pushed;
            }

            vec_excess[u] += f;
            vec_excess[edge_node(in_he.e)] -= f;
            in_he.flow -= flowSent(f);

            if (flowSent(prevFlow) > 0 && flowReceived(in_he.flow) >= 0) //v previously sent flow and now either has none, or receives flow.
                removePinFromFlowPins(in_he, false);
            if (flowReceived(in_he.flow) > 0 && flowReceived(prevFlow) <= 0)  //v now receives flow and did not previously, thus must be inserted into pins_receiving_flow
                insertPinIntoFlowPins(in_he, true);
        }

        inline Flow& excess_node(Node u) { assert(u < numNodes()); return vec_excess[u]; }
        inline Flow& excess_edge(Hyperedge e) { assert(e < numHyperedges()); return vec_excess[numNodes() + e]; }

        inline Node edge_node(Hyperedge e) { assert(e < numHyperedges()); return Node(numNodes() + e); }

        inline Hyperedge edgeFromLawlerNode(Node u) { assert(u >= numNodes()); return Hyperedge(u - numNodes()); }

        inline size_t numLawlerNodes() const {
            return numNodes() + numHyperedges();
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

    private:
        std::vector<Flow> vec_excess;
        std::vector<size_t> vec_label;
    };
}