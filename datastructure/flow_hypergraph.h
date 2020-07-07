#pragma once

#include "../definitions.h"
#include "../util/unused.h"
#include <boost/dynamic_bitset.hpp>

namespace whfc {

	//TODO. when all flow algorithms are implemented. Reassess whether we want to store the flow at the pins or the incidences. and whether we actually want array of structs
	
	class FlowHypergraph {
	public:
		struct Pin {
			Node pin = invalidNode;
			InHeIndex he_inc_iter;
			Flow flow = Flow(0);
			bool operator==(const Pin& o) const { return o.pin == pin && o.he_inc_iter == he_inc_iter; }
		};

		struct InHe {	//Hyperedge Incidence
			Hyperedge e = invalidHyperedge;
			std::array<PinIndex, 2> pin_iters;
		};

		struct HyperedgeData {
			PinIndex first_out = PinIndex(0);
			Flow flow = Flow(0);
			Flow capacity = Flow(0);
		};

		struct NodeData {
			InHeIndex first_out = InHeIndex(0);
			NodeWeight weight = NodeWeight(0);
		};
		
		using PinRange = mutable_range<std::vector<Pin>>;
		using PinIterator = PinRange::iterator;
		using PinIndexRange = mutable_index_range<PinIndex>;
		using InHeRange = mutable_range<std::vector<InHe>>;
		using InHeIterator = InHeRange::iterator;
		using InHeIndexRange = mutable_index_range<InHeIndex>;

		inline auto nodeIDs() const { return mutable_index_range<Node>(Node(0), Node::fromOtherValueType(numNodes())); }
		inline auto hyperedgeIDs() const { return mutable_index_range<Hyperedge>(Hyperedge(0), Hyperedge::fromOtherValueType(numHyperedges())); }
		inline auto pinIndices() const { return PinIndexRange(PinIndex(0), PinIndex::fromOtherValueType(numPins())); }

		FlowHypergraph() : maxNumNodes(0), maxNumHyperedges(0), nodes(1), hyperedges(1) { }
		
		//use in FlowHypergraphBuilder
		FlowHypergraph(size_t maxNumNodes, size_t maxNumHyperedges) :
                maxNumNodes(maxNumNodes), maxNumHyperedges(maxNumHyperedges), nodes(1), hyperedges(1) { }

		FlowHypergraph(std::vector<NodeWeight>& node_weights, std::vector<HyperedgeWeight>& hyperedge_weights, std::vector<PinIndex>& hyperedge_sizes, std::vector<Node>& _pins) :
                maxHyperedgeCapacity(0),
		        maxNumNodes(node_weights.size()),
                maxNumHyperedges(hyperedge_weights.size()),
				nodes(node_weights.size() + 1),
				hyperedges(hyperedge_weights.size() + 1),
				pins_in(_pins.size()),
				pins_out(_pins.size()),
				incident_hyperedges(_pins.size()),
				pins_in_sending_flow_end(hyperedge_weights.size()),
				pins_out_receiving_flow_end(hyperedge_weights.size()),
				total_node_weight(boost::accumulate(node_weights, NodeWeight(0)))
		{
			size_t i = 0;
			for (const Node p : _pins) {
			    pins_in[i].pin = p;
				pins_out[i++].pin = p;					//copy pins
				nodes[p + 1].first_out++;			//bucket sizes
			}
			
			for (Node u : nodeIDs()) {
				nodes[u + 1].first_out += nodes[u].first_out;			//prefix sum
				nodes[u].weight = node_weights[u];						//copy node weights
			}
			
			for (Hyperedge e : hyperedgeIDs()) {
				hyperedges[e].capacity = hyperedge_weights[e];
				hyperedges[e+1].first_out = hyperedges[e].first_out + hyperedge_sizes[e];		//prefix sum
				for (auto pin_it = beginIndexPins(e); pin_it != endIndexPins(e); pin_it++) {
					Pin& p_in = pins_in[pin_it];
					Pin& p_out = pins_out[pin_it];
					InHeIndex ind_he = nodes[p_in.pin].first_out++;							//destroy first_out temporarily and reset later
					InHe& inc_he = incident_hyperedges[ind_he];
					inc_he.e = e;
					inc_he.pin_iters = { pin_it, pin_it };
					p_in.he_inc_iter = ind_he;					//set iterator for incident hyperedge -> its position in incident_hyperedges of the node
					p_out.he_inc_iter = ind_he;
				}
				maxHyperedgeCapacity = std::max(maxHyperedgeCapacity, hyperedges[e].capacity);
			}
			
			for (Node u(numNodes()-1); u > 0; u--)
				nodes[u].first_out = nodes[u-1].first_out;	//reset temporarily destroyed first_out
			nodes[0].first_out = InHeIndex(0);
			
			PinIndex x = PinIndex(0);
			for (Hyperedge e : hyperedgeIDs()) {
				pins_in_sending_flow_end[e] = x;
				x += pinCount(e);
				pins_out_receiving_flow_end[e] = x;
			}
		}
		

		bool hasNodeWeights() const { return std::any_of(nodes.begin(), nodes.begin() + numNodes(), [](const NodeData& u) { return u.weight > 1; }); }
		bool hasHyperedgeWeights() const { return std::any_of(hyperedges.begin(), hyperedges.begin() + numHyperedges(), [](const HyperedgeData& e) { return e.capacity > 1; }); }
		inline size_t numNodes() const { return nodes.size() - 1 ; }
		inline size_t numHyperedges() const { return hyperedges.size() - 1; }
		inline size_t numPins() const { return pins_in.size(); }
		inline PinIndex pinCount(const Hyperedge e) const { return hyperedges[e+1].first_out - hyperedges[e].first_out; }
		inline InHeIndex degree(const Node u) const { return nodes[u+1].first_out - nodes[u].first_out; }
		inline NodeWeight totalNodeWeight() const { return total_node_weight; }
		inline NodeWeight nodeWeight(const Node u) const { return nodes[u].weight; }
		inline NodeWeight& nodeWeight(const Node u) { return nodes[u].weight; }

		inline InHeIndex beginIndexHyperedges(Node u) const { return nodes[u].first_out; }
		inline InHeIndex endIndexHyperedges(Node u) const { return nodes[u+1].first_out; }
		//interface for irange is front(), drop_front(), empty(), size(), begin(), end()
		inline InHeIndexRange incidentHyperedgeIndices(const Node u) const {
			return InHeIndexRange(beginIndexHyperedges(u), endIndexHyperedges(u));
		}
		inline InHe& getInHe(const InHeIndex ind_e) { return incident_hyperedges[ind_e]; }
		inline InHe& getInHe(const Pin& pin) { return getInHe(pin.he_inc_iter); }
		inline const InHe& getInHe(const InHeIndex ind_e) const { return incident_hyperedges[ind_e]; }
		inline const InHe& getInHe(const Pin& pin) const { return getInHe(pin.he_inc_iter); }

		inline PinIndex beginIndexPins(const Hyperedge e) const { return hyperedges[e].first_out; }
		inline PinIndex endIndexPins(const Hyperedge e) const { return hyperedges[e+1].first_out; }
		inline PinIndexRange pinIndices(const Hyperedge e) const { return PinIndexRange(beginIndexPins(e), endIndexPins(e)); }
		//inline PinIndexRange pinsSendingFlowIndices(const Hyperedge e) const { return pins_sending_flow[e]; }
		//inline Pin& getPinIn(const PinIndex ind_p) { return pins[ind_p]; }
		inline Pin& getPinIn(const InHe& inc_p) { return pins_in[inc_p.pin_iters[sends_index]]; }
		//inline const Pin& getPin(const PinIndex ind_p) const { return pins[ind_p]; }
		inline const Pin& getPinIn(const InHe& inc_p) const { return pins_in[inc_p.pin_iters[sends_index]]; }
		inline Pin& getPinOut(const InHe& inc_p) { return pins_out[inc_p.pin_iters[1-sends_index]]; }
		inline const Pin& getPinOut(const InHe& inc_p) const { return pins_out[inc_p.pin_iters[1-sends_index]]; }


		inline InHeIterator beginHyperedges() { return incident_hyperedges.begin(); }
		inline InHeIterator endHyperedges() { return incident_hyperedges.end(); }
		inline InHeIterator beginHyperedges(const Node u) { return incident_hyperedges.begin() + nodes[u].first_out; }
		inline InHeIterator endHyperedges(const Node u) { return incident_hyperedges.begin() + nodes[u+1].first_out; }
		InHeRange hyperedgesOf(const Node u) { return InHeRange(beginHyperedges(u), endHyperedges(u)); }
		InHeRange hyperedgesInRange(const InHeIndexRange hir) { return InHeRange(beginHyperedges() + hir.begin(), beginHyperedges() + hir.end()); }

		inline PinIterator beginPinsIn() { return pins_in.begin(); }
		inline PinIterator endPinsIn() { return pins_in.end(); }
		inline PinIterator beginPinsOut() { return pins_out.begin(); }
		inline PinIterator endPinsOut() { return pins_out.end(); }
		inline PinIterator beginPinsIn(const Hyperedge e) { return pins_in.begin() + hyperedges[e].first_out; }
		inline PinIterator endPinsIn(const Hyperedge e) { return pins_in.begin() + hyperedges[e+1].first_out; }
        inline PinIterator beginPinsOut(const Hyperedge e) { return pins_out.begin() + hyperedges[e].first_out; }
        inline PinIterator endPinsOut(const Hyperedge e) { return pins_out.begin() + hyperedges[e+1].first_out; }
		PinRange pinsInOf(const Hyperedge e) { return PinRange(beginPinsIn(e), endPinsIn(e)); }
		PinRange pinsOutOf(const Hyperedge e) { return PinRange(beginPinsOut(e), endPinsOut(e)); }
		PinRange pinsInInRange(const PinIndexRange pir) { return PinRange(beginPinsIn() + pir.begin(), beginPinsIn() + pir.end()); }
        PinRange pinsOutInRange(const PinIndexRange pir) { return PinRange(beginPinsOut() + pir.begin(), beginPinsOut() + pir.end()); }

        PinRange pinsSendingFlowInto(const Hyperedge e) { return PinRange(beginPinsIn(e), beginPinsIn(e) + pins_in_sending_flow_end[e]); }
		PinRange pinsReceivingFlowFrom(const Hyperedge e) { return PinRange(beginPinsOut(e), beginPinsOut(e) + pins_out_receiving_flow_end[e]); }

		/*
		PinIndexRange pinsNotSendingFlowIndices(const Hyperedge e) const {
			if (forwardView()) {
				return PinIndexRange(pins_sending_flow[e].end(), endIndexPins(e));
			}
			else {
				return PinIndexRange(beginIndexPins(e), pins_sending_flow[e].begin());
			}
		}*/

		/*
		PinRange pinsNotSendingFlowInto(const Hyperedge e) {
			return pinsInRange(pinsNotSendingFlowIndices(e));
		}*/

		/*
		PinRange pinsInWithoutFlow(const Hyperedge e) {
			return pinsInRange(pins_without_flow(e));
		}*/

		inline bool forwardView() const { return sends_multiplier == 1; }
		void flipViewDirection() {
		    std::swap(pins_in_sending_flow_end, pins_out_receiving_flow_end);
		    std::swap(pins_in, pins_out),
		    std::swap(sends_multiplier, receives_multiplier);
		    sends_index = 1 - sends_index;
		}
		
		
		inline Flow capacity(const Hyperedge e) const { return hyperedges[e].capacity; }
		inline Flow& capacity(const Hyperedge e) { return hyperedges[e].capacity; }
		inline Flow flow(const Hyperedge e) const { return hyperedges[e].flow; }
		inline Flow& flow(const Hyperedge e) { return hyperedges[e].flow; }
		inline Flow residualCapacity(const Hyperedge e) const { return capacity(e) - flow(e); }
		inline bool isSaturated(const Hyperedge e) const { assert(flow(e) <= capacity(e)); return flow(e) == capacity(e); }

		//inline Flow flowSent(const Flow f) const { return f * sends_multiplier; }
		//flow sent from u = getPin(inc_u.pin_iter).pin into e = inc_u.e
		//inline Flow flowSent(const InHe& inc_u) const { return flowSent(inc_u.flow); }
		inline Flow absoluteFlowSent(const InHe& inc_u) const { return getPinIn(inc_u).flow; }
		//inline Flow flowSent(const Pin& pin) const { return flowSent(getInHe(pin)); }
		//inline Flow absoluteFlowSent(const Pin& pin) const { return pin.flow; }

		inline Flow pinFlow(const Pin& pin) const { return pin.flow; }

		//for testing only
		/*
		inline Flow flowSent(const Node u) {
			Flow f = 0;
			for (const InHe& inc_u : hyperedgesOf(u))
				f += flowSent(inc_u.flow);
			return f;
		}*/

		/*
		inline Flow flowReceived(const Node u) {
			return -flowSent(u);
		}*/

		//inline Flow flowReceived(const Flow f) const { return f * receives_multiplier; }
		//flow that u = getPin(inc_u.pin_iter).pin receives from e = inc_u.e
		//inline Flow flowReceived(const InHe& inc_u) const { return flowReceived(inc_u.flow); }
		//inline Flow flowReceived(const Pin& pin) const { return flowReceived(getInHe(pin)); }
		inline Flow absoluteFlowReceived(const InHe& inc_u) const { return getPinOut(inc_u).flow; }
		//inline Flow flowReceived(const Pin& pin) const { return flowReceived(getInHe(pin)); }

		inline Flow residualCapacity(const InHe& inc_u, InHe& inc_v) const {
			return absoluteFlowReceived(inc_u) + absoluteFlowSent(inc_v) + residualCapacity(inc_u.e);
		}
		
		//for testing only
		InHe& findIncidence(const Node u, const Hyperedge e) {
			for (InHe& x : hyperedgesOf(u))
				if (x.e == e)
					return x;
			throw std::out_of_range("e is not in the list of incident hyperedges of u");
		}
		
		//for testing only
		/*
		Pin& findPin(const Hyperedge e, const Node v) {
			for (Pin& x : pinsOf(e))
				if (x.pin == v)
					return x;
			throw std::out_of_range("v is not a pin of e");
		}*/

		/*
		void routeFlow(InHe& inc_u, InHe& inc_v, Flow flow_delta) {
			const Hyperedge e = inc_u.e;
			assert(inc_u.e == inc_v.e && "Routing flow but incident hyperedges are not the same");
			assert(flow_delta > 0 && "Routing <= 0 flow.");
			assert(flow_delta <= residualCapacity(inc_u, inc_v) && "Routing more flow than residual capacity");
			assert(flow(e) <= capacity(e) && "Capacity on e already exceeded");
			assert(std::abs(inc_u.flow) <= capacity(e) && "Pin capacity already violated (u)");
			assert(std::abs(inc_v.flow) <= capacity(e) && "Pin capacity already violated (v)");

			const Flow prevFlowU = inc_u.flow;
			const Flow prevFlowV = inc_v.flow;
			
			Flow flow_delta_on_v_eOut_eIn_u = std::min({ absoluteFlowSent(inc_v), absoluteFlowReceived(inc_u), flow_delta });
			inc_u.flow += flowSent(flow_delta_on_v_eOut_eIn_u);
			inc_v.flow += flowReceived(flow_delta_on_v_eOut_eIn_u);
			flow(e) -= flow_delta_on_v_eOut_eIn_u;
			flow_delta -= flow_delta_on_v_eOut_eIn_u;
			
			Flow flow_delta_on_v_eOut_u = std::min(flow_delta, absoluteFlowReceived(inc_u));
			inc_u.flow += flowSent(flow_delta_on_v_eOut_u);
			inc_v.flow += flowReceived(flow_delta_on_v_eOut_u);
			//does not influence flow(e)
			flow_delta -= flow_delta_on_v_eOut_u;
			
			Flow flow_delta_on_v_eIn_u = std::min(flow_delta, absoluteFlowSent(inc_v));
			inc_u.flow += flowSent(flow_delta_on_v_eIn_u);
			inc_v.flow += flowReceived(flow_delta_on_v_eIn_u);
			//does not influence flow(e)
			flow_delta -= flow_delta_on_v_eIn_u;
			
			Flow flow_delta_on_v_eIn_eOut_u = flow_delta;
			assert(flow_delta_on_v_eIn_eOut_u <= residualCapacity(e));
			inc_u.flow += flowSent(flow_delta_on_v_eIn_eOut_u);
			inc_v.flow += flowReceived(flow_delta_on_v_eIn_eOut_u);
			flow(e) += flow_delta_on_v_eIn_eOut_u;

			
			if (flowReceived(prevFlowU) > 0 && flowSent(inc_u.flow) >= 0)	//u previously received flow and now either has none, or sends flow.
				removePinFromFlowPins(inc_u, true);
			if (flowSent(inc_u.flow) > 0 && flowSent(prevFlowU) <= 0) //u now sends flow and did not previously, thus must be inserted into pins_sending_flow
				insertPinIntoFlowPins(inc_u, false);

			if (flowSent(prevFlowV) > 0 && flowReceived(inc_v.flow) >= 0) //v previously sent flow and now either has none, or receives flow.
				removePinFromFlowPins(inc_v, false);
			if (flowReceived(inc_v.flow) > 0 && flowReceived(prevFlowV) <= 0)  //v now receives flow and did not previously, thus must be inserted into pins_receiving_flow
				insertPinIntoFlowPins(inc_v, true);
			
			assert(check_flow_conservation_locally(e));
			assert(sanity_check_pin_ranges(e));
			assert(assert_backpointers_correct(inc_u));
			assert(assert_backpointers_correct(inc_v));
			assert(pin_is_categorized_correctly(inc_u) && "Pin categorized incorrectly");
			assert(pin_is_categorized_correctly(inc_v) && "Pin categorized incorrectly");
		}
		 */

		std::vector<NodeData>& getNodes() {
		    return nodes;
		}
		
		Flow maxHyperedgeCapacity = maxFlow;

        size_t maxNumNodes;
        size_t maxNumHyperedges;

	public:
		std::vector<NodeData> nodes;
		std::vector<HyperedgeData> hyperedges;
		std::vector<Pin> pins_in;
		std::vector<Pin> pins_out;
		std::vector<InHe> incident_hyperedges;

		//TODO get rid of the range and just store one index, if this turns out to be cache inefficient later on
		std::vector<PinIndex> pins_in_sending_flow_end;	//indexed by hyperedge id. gives range of pin ids/iterators sending flow to that hyperedge. grows right if forwardView = true
		std::vector<PinIndex> pins_out_receiving_flow_end;	//indexed by hyperedge id. gives range of pin ids/iterators receiving flow from that hyperedge. grows left if forwardView = true
		
		NodeWeight total_node_weight = NodeWeight(0);
		int sends_multiplier = 1;						//if forwardView = true, flow entering hyperedge e should be positive and flow exiting e should be negative. reverse, if forwardView = false.
		int receives_multiplier = -1;
		uint sends_index = 0;

		/*PinIndexRange pinsIn_without_flow(const Hyperedge e) const {
			return PinIndexRange();
		}*/


		PinIndex removePinFromFlowPins(InHe& inc_u, bool flow_receiving_pins) {
			const Hyperedge e = inc_u.e;
			PinIndex it_u = inc_u.pin_iter;
			PinIndexRange& flow_pins = flow_receiving_pins ? pins_receiving_flow[e] : pins_sending_flow[e];
			assert(!flow_pins.empty());
			assert(flow_pins.contains(it_u));

			PinIndex it_o = (forwardView() == flow_receiving_pins) ? flow_pins.begin() : PinIndex(flow_pins.end() - 1);
			InHe& inc_o = getInHe(getPin(it_o));
			assert(it_o == it_u || (flow_receiving_pins ? flowReceived(inc_o) > 0 : flowSent(inc_o) > 0));	//ensure it_o, taken from flow_pins, actually receives or sends flow, as appropriate
			if (forwardView() == flow_receiving_pins)
				flow_pins.advance_begin();
			else
				flow_pins.retreat_end();
			std::swap(inc_u.pin_iter, inc_o.pin_iter);
			std::swap(pins[it_u], pins[it_o]);
			assert(pins_without_flow(e).contains(it_o));
			return it_o;
		}

		PinIndex insertPinIntoFlowPins(InHe& inc_u, bool flow_receiving_pins) {
			const Hyperedge e = inc_u.e;
			PinIndex it_u = inc_u.pin_iter;
			PinIndexRange& flow_pins = flow_receiving_pins ? pins_receiving_flow[e] : pins_sending_flow[e];
			assert(pins_without_flow(e).contains(it_u));
			PinIndex it_o = (forwardView() == flow_receiving_pins) ? PinIndex(flow_pins.begin() - 1) : flow_pins.end();
			InHe& inc_o = getInHe(getPin(it_o));
			if (forwardView() == flow_receiving_pins)
				flow_pins.retreat_begin();
			else
				flow_pins.advance_end();
			std::swap(inc_u.pin_iter, inc_o.pin_iter);
			std::swap(pins[it_u], pins[it_o]);
			assert(flow_pins.contains(it_o));
			return it_o;
		}
		
		bool check_flow_conservation_locally(const Hyperedge e) {
			Flow f_in = 0, f_out = 0;
			for (Pin& p : pinsInOf(e)) {
				f_in += pinFlow(p);
				assert(pinFlow(p) <= capacity(e));
			}
			for (Pin& p : pinsOutOf(e)) {
			    f_out += pinFlow(p);
			    assert(pinFlow(p) <= capacity(e));
			}
			Flow f = f_in - f_out;
			assert(f == 0);
			assert(f_in == flow(e));
			return true;
		}

		bool assert_backpointers_correct(const InHe& inhe) const {
			const InHe& doubled_in = getInHe(getPinIn(inhe));
			const InHe& doubled_out = getInHe(getPinOut(inhe));
			unused(inhe, doubled_in, doubled_out);
			assert(doubled_in.pin_iter_in == inhe.pin_iter_in && "Backpointer Pin Iter inconsistent");
			assert(doubled_in.e == inhe.e && "Backpointer hyperedge ID inconsistent");
			return true;
		}

		/*
		bool sanity_check_pin_ranges(const Hyperedge e) const {
			//check left / right end of pin ranges agree with first_out
			const PinIndexRange& s = forwardView() ? pins_sending_flow[e] : pins_receiving_flow[e];
			const PinIndexRange& l = !forwardView() ? pins_sending_flow[e] : pins_receiving_flow[e];
			unused(e, s, l);
			assert(hyperedges[e].first_out == s.begin());
			assert(hyperedges[e+1].first_out == l.end());

			assert(s.begin() <= s.end());
			assert(s.end() <= l.begin());
			assert(l.begin() <= l.end());
			return true;
		}*/

		bool pin_is_categorized_correctly(const InHe& inc_u) {
			const Hyperedge e = inc_u.e;
			const PinIndex it_u = inc_u.pin_iter;
			bool sends = flowSent(inc_u) > 0 && pins_sending_flow[e].contains(it_u);
			bool receives = flowReceived(inc_u.flow) > 0 && pins_receiving_flow[e].contains(it_u);
			bool no_flow = inc_u.flow == 0 && pins_without_flow(e).contains(it_u);
			return 		(sends && !receives && !no_flow)
					|| 	(!sends && receives && !no_flow)
					|| 	(!sends && !receives && no_flow);
		}
		
		static_assert(std::is_trivially_destructible<Pin>::value);
		static_assert(std::is_trivially_destructible<InHe>::value);
		static_assert(std::is_trivially_destructible<HyperedgeData>::value);
		static_assert(std::is_trivially_destructible<NodeData>::value);
		static_assert(std::is_trivially_destructible<PinIndexRange>::value);

	public:

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
				out << e << " pincount = " << pinCount(e) << " w= " << capacity(e) << " pins (pin,flow) = [";
				for (const Pin& u : pinsOf(e)) {
					assert(pin_is_categorized_correctly(getInHe(u)));
					out << "(" << u.pin << "," << getInHe(u).flow << ") ";
				}
				out << "]" << "\n";
			}
			out << std::flush;
		}

		void printHypergraph(std::ostream& out) {
			printNodes(out);
			printHyperedges(out);
		}

		friend std::ostream& operator<<(std::ostream& out, FlowHypergraph& hg) noexcept {
			hg.printHypergraph(out);
			return out;
		}

	};
}