#pragma once

#include "../io/hmetis_io.h"
#include "../logger.h"
#include "../algorithm/ford_fulkerson.h"

namespace whfc {
namespace Test {
	
	class FlowHypergraphTests {
	private:
		using InHe = FlowHypergraph::InHe;
		using Pin = FlowHypergraph::Pin;
		
		template<typename T>
		bool compareRanges(std::vector<T> actual, std::vector<T>& expected) {
			std::sort(actual.begin(), actual.end());
			std::sort(expected.begin(), expected.end());
			return actual == expected;
		}
		
	public:
		static constexpr bool debug = true;
		void routeFlowTests() {
			std::string hgfile = "../test_hypergraphs/twocenters.hgr";
			//FlowHypergraph hg = HMetisIO::readFlowHypergraph(hgfile);
			FlowHypergraphBuilder hg = HMetisIO::readFlowHypergraphWithBuilder(hgfile);
			
			Node u(1), v(11);
			Hyperedge e(3);
			InHe& inc_u = hg.findIncidence(u,e);
			InHe& inc_v = hg.findIncidence(v,e);
			Pin& pu = hg.findPin(e, u);
			Pin& pv = hg.findPin(e, v);
			
			Assert(pu.pin == u);
			Assert(pv.pin == v);
			Assert(hg.getPin(inc_u) == pu);
			Assert(hg.getPin(inc_v) == pv);
			
			hg.routeFlow(inc_u, inc_v, 1);
			
			Assert(hg.flow(e) == 1);
			Assert(hg.flowSent(u) == 1);
			Assert(hg.flowReceived(v) == 1);
			Assert(hg.flowSent(inc_u) == 1);
			Assert(hg.flowReceived(inc_v) == 1);
			Assert(hg.flowSent(hg.findPin(e, u)) == 1);
			Assert(hg.flowReceived(hg.findPin(e, v)) == 1);
			
			Node w(9);
			Hyperedge e2(2);
			InHe& inc_w_e2 = hg.findIncidence(w, e2);
			InHe& inc_u_e2 = hg.findIncidence(u, e2);
			hg.routeFlow(inc_u_e2, inc_w_e2, 1);
			
			Assert(hg.flow(e2) == 1);
			Assert(hg.flowSent(inc_u_e2) == 1);
			Assert(hg.flowReceived(inc_w_e2) == 1);
			Assert(hg.flowSent(u) == 2);
			Assert(hg.flowReceived(w) == 1);
			
			hg.routeFlow(inc_v, inc_u, 1);
			Assert(hg.flow(e) == 0);
			Assert(hg.flowSent(inc_u) == 0);
			Assert(hg.flowReceived(inc_v) == 0);
		}
		
		template<typename FlowAlgo>
		bool tryFlowAlgo(std::string file, Flow expected_flow, Node s, Node t) {
			LOG << V(file);
			FlowHypergraph hg = HMetisIO::readFlowHypergraph(file);
			CutterState<FlowAlgo> cs(hg, NodeWeight(4000));
			cs.initialize(s,t);
			FlowAlgo flow(hg);
			Flow f = flow.exhaustFlow(cs);
			LOG << V(f);
			return f == expected_flow;
		}
		
		void flowAlgoTest(std::string file, Flow expected_flow, Node s, Node t) {
			LOG << "basic ff";
			Assert(tryFlowAlgo<BasicFordFulkerson>(file, expected_flow, s, t));
			LOG << "scaling ff";
			Assert(tryFlowAlgo<ScalingFordFulkerson>(file, expected_flow, s, t));
			LOG << "basic ek";
			Assert(tryFlowAlgo<BasicEdmondsKarp >(file, expected_flow, s, t));
			LOG << "scaling ff";
			Assert(tryFlowAlgo<ScalingEdmondsKarp>(file, expected_flow, s, t));
		}
		
		void run() {
			routeFlowTests();
			flowAlgoTest("../test_hypergraphs/testhg.hgr", Flow(1), Node(14), Node(10));
			flowAlgoTest("../test_hypergraphs/twocenters.hgr", Flow(2), Node(0), Node(2));
			flowAlgoTest("../test_hypergraphs/twocenters.hgr", Flow(2), Node(0), Node(3));
			flowAlgoTest("../test_hypergraphs/push_back.hgr", Flow(6), Node(0), Node(7));
			
		}
	};
}
}