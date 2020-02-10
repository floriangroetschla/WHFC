#pragma once

#include "../util/range.h"

#include <boost/range/irange.hpp>

#include <cstdint>
#include <vector>

namespace whfc_rb {

using Index = uint32_t;

template<typename RowID, typename ColumnID>
class CSR {
public:
	
	explicit CSR(RowID numRows) : first_out(numRows + 1, 0) {
	
	}
	
	std::vector<Index> first_out;
	std::vector<ColumnID> incidence;
	
	Index numRows() const {
		return static_cast<Index>(first_out.size() - 1);
	}
	
	Index degree(const RowID u) const {
		return first_out[u+1] - first_out[u];
	}
	
	const_range< std::vector<ColumnID> > incidentElementsOf(const RowID u) const {
		return const_range(incidence, first_out[u], first_out[u+1]);
	}

	auto rows() const {
		return boost::irange<RowID>(0, numRows());
	}
	
	void finishRow(const RowID u) {
		first_out[u + 1] = incidence.size();
	}
};


class CSRHypergraph {
public:
	using HyperedgeWeight = uint32_t;
	using HyperedgeID = uint32_t;
	using NodeWeight = uint32_t;
	using NodeID = uint32_t;
	
	CSRHypergraph(const NodeID numNodes = 0, const HyperedgeID numHyperedges = 0) :
			Vertices(numNodes + 1), /* +1 is a hack to abuse the first_out array for the prefix sum without reverse sweep reset */
			E(numHyperedges), node_weights(numNodes, 0), hyperedge_weights(numHyperedges, 0)
	{
		// this constructor only prepares the CSRHypergraph object for reading from a hMetis file
	}
	
	const_range< std::vector<HyperedgeID> > hyperedgesOf(const NodeID u) const {
		return Vertices.incidentElementsOf(u);
	}
	
	Index degree(const NodeID u) const {
		return Vertices.degree(u);
	}
	
	const_range< std::vector<NodeID> > pinsOf(const HyperedgeID e) const {
		return E.incidentElementsOf(e);
	}
	
	Index pinCount(const HyperedgeID e) const {
		return E.degree(e);
	}
	
	auto nodes() const {
		return Vertices.rows();
	}
	
	NodeID numNodes() const {
		return Vertices.numRows();
	}
	
	auto hyperedges() const {
		return E.rows();
	}
	
	HyperedgeID numHyperedges() const {
		return E.numRows();
	}
	
	std::vector<NodeID>& pins() {
		return E.incidence;
	}
	
	NodeID numPins() const {
		return static_cast<NodeID>(E.incidence.size());
	}
	
	std::vector<HyperedgeID>& xpins() {
		return Vertices.incidence;
	}
	
	NodeWeight& nodeWeight(const NodeID u) {
		return node_weights[u];
	}
	
	HyperedgeWeight& hyperedgeWeight(const HyperedgeID u) {
		return hyperedge_weights[u];
	}
	
	void finishHyperedge(const HyperedgeID e) {
		E.finishRow(e);
	}
	
	void computeXPins() {
		// still need to compute Vertices.first_out
		// and Vertices.incidence
		
		for (NodeID u : pins()) {
			Vertices.first_out[u + 2]++;
		}
		
		for (NodeID u = 0; u < Vertices.first_out.size() - 2; ++u) {
			Vertices.first_out[u + 2] += Vertices.first_out[u + 1];
		}
		
		Vertices.incidence.resize(numPins());
		for (HyperedgeID e : hyperedges()) {
			for (NodeID pin : pinsOf(e)) {
				Vertices.incidence[ Vertices.first_out[pin + 1]++ ] = e;
			}
		}
		
		Vertices.first_out.pop_back();	// remove the last element
	}
	
	void printNodes(std::ostream& out) {
		out << "---Nodes---\n";
		for (NodeID u : nodes()) {
			out << u + 1 << " deg = " << degree(u) << " w= " << nodeWeight(u) << " inc_hes = [";
			for (const HyperedgeID e : hyperedgesOf(u))
				out << e + 1 << " ";
			out << "]" << "\n";
		}
		out << std::flush;
	}
	
	void printHyperedges(std::ostream& out) {
		out << "---Hyperedges---\n";
		for (HyperedgeID e: hyperedges()) {
			out << e + 1<< " pincount = " << pinCount(e) << " w= " << hyperedgeWeight(e) << " pins (pin,flow) = [";
			for (const NodeID u : pinsOf(e)) {
				out << u + 1 << " ";
			}
			out << "]" << "\n";
		}
		out << std::flush;
	}
	
	void printHypergraph(std::ostream& out) {
		printNodes(out);
		printHyperedges(out);
	}
	
	friend std::ostream& operator<<(std::ostream& out, CSRHypergraph& hg) noexcept {
		hg.printHypergraph(out);
		return out;
	}
	
protected:
	CSR<NodeID, HyperedgeID> Vertices;		// TODO rename. cannot be V
	CSR<HyperedgeID, NodeID> E;
	std::vector<NodeWeight> node_weights;
	std::vector<HyperedgeWeight> hyperedge_weights;
};

}	// namespace whfc_rb