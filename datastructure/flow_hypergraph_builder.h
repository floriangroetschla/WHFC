#pragma once

#include <atomic>
#include <tbb/parallel_for.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for_each.h>
#include "flow_hypergraph.h"
#include "../recursive_bisection/mock_builder.h"

/*
 * hyperedges with zero/one pins are removed automatically during build process
 */

namespace whfc {
	class FlowHypergraphBuilder : public FlowHypergraph {
	public:
		using Base = FlowHypergraph;
		
		FlowHypergraphBuilder() : Base() {
			clear();
		}
		
		explicit FlowHypergraphBuilder(size_t nNodes) {
			reinitialize(nNodes);
		}
		
		//use to get rid of any allocations
		FlowHypergraphBuilder(size_t maxNumNodes, size_t maxNumHyperedges) : Base(maxNumNodes, maxNumHyperedges) {
			//don't do clean-up here yet, so that we can use the numbers for allocating the remaining datastructures
		}
		
		void clear() {
			finalized = false;
			numPinsAtHyperedgeStart = 0;
			maxHyperedgeCapacity = 0;
			
			nodes.clear();
			hyperedges.clear();
			pins.clear();
			incident_hyperedges.clear();
			pins_sending_flow.clear();
			pins_receiving_flow.clear();
			total_node_weight = NodeWeight(0);
			sends_multiplier = 1;
			receives_multiplier = -1;
			
			//sentinels
			nodes.push_back({InHeIndex(0), NodeWeight(0)});
			hyperedges.push_back({PinIndex(0), Flow(0), Flow(0)});
		}
		
		void reinitialize(size_t numNodes) {
			clear();
			nodes.resize(numNodes + 1);
		}
		
		void addNode(const NodeWeight w) {
			nodes.back().weight = w;
			nodes.push_back({InHeIndex(0), NodeWeight(0)});
		}
		
		void startHyperedge(const Flow capacity) {
			finishHyperedge();	//finish last hyperedge
			hyperedges.back().capacity = capacity;	//exploit sentinel
			numPinsAtHyperedgeStart = numPins();
			maxHyperedgeCapacity = std::max(maxHyperedgeCapacity, capacity);
		}
		
		void addPin(const Node u) {
			assert(u < numNodes());
			pins.push_back({u, InHeIndex::Invalid()});
			nodes[u+1].first_out++;
		}
		
		size_t currentHyperedgeSize() const {
			return numPins() - numPinsAtHyperedgeStart;
		}
		
		void removeCurrentHyperedge() {
			while (numPins() > numPinsAtHyperedgeStart)
				removeLastPin();
		}

		void addMockBuildersParallel(tbb::enumerable_thread_specific<MockBuilder>& mockBuilder_thread_specific, std::vector<whfc::Node>& globalToLocalID, whfc::Node source, whfc::Node target) {
		    struct WorkElement {
		        MockBuilder& builder;
		        size_t hyperedgeStartIndex;
		        size_t pinStartIndex;
		    };

		    std::vector<WorkElement> workElements;

            finishHyperedge();

		    size_t numberOfHyperedges = numHyperedges();
		    size_t numberOfPins = numPins();


		    // resize datastructures
		    for (MockBuilder& builder : mockBuilder_thread_specific) {
                builder.finishHyperedge();
                workElements.push_back({builder, numberOfHyperedges, numberOfPins});
		        numberOfHyperedges += builder.numHyperedges();
		        numberOfPins += builder.numPins();
		    }

		    hyperedges.resize(numberOfHyperedges + 1);
		    pins_sending_flow.resize(numberOfHyperedges);
		    pins_receiving_flow.resize(numberOfHyperedges);
		    pins.resize(numberOfPins);

            tbb::enumerable_thread_specific<size_t> sourceOccurrences(0);
            tbb::enumerable_thread_specific<size_t> targetOccurrences(0);

		    tbb::parallel_for_each(workElements, [&](WorkElement& workElement) {
                size_t pinCounter = workElement.pinStartIndex;
                size_t pinOffset = pinCounter;

                MockBuilder& builder = workElement.builder;
                size_t& sourceOcc = sourceOccurrences.local();
                size_t& targetOcc = targetOccurrences.local();

                for (size_t i = 0; i < builder.pins.size(); ++i) {
                    pins[pinCounter++] = {globalToLocalID[builder.pins[i]], InHeIndex::Invalid()};
                    if (builder.pins[i] == source) {
                        sourceOcc++;
                    } else if (builder.pins[i] == target) {
                        targetOcc++;
                    } else {
                        __atomic_fetch_add(&nodes[globalToLocalID[builder.pins[i]]+1].first_out.value(), 1, __ATOMIC_RELAXED);
                    }
                }

                for (size_t i = 0; i < builder.numHyperedges(); ++i) {
                    pins_sending_flow[workElement.hyperedgeStartIndex + i] = PinIndexRange(PinIndex(pinOffset) + builder.hyperedges[i].first_out, PinIndex(pinOffset) + builder.hyperedges[i].first_out);
                    hyperedges[workElement.hyperedgeStartIndex + i] = {PinIndex(pinOffset) + builder.hyperedges[i].first_out, Flow(0), builder.hyperedges[i].capacity};
                    pins_receiving_flow[workElement.hyperedgeStartIndex + i] = PinIndexRange(PinIndex(pinOffset) + builder.hyperedges[i+1].first_out, PinIndex(pinOffset) + builder.hyperedges[i+1].first_out);
                    maxHyperedgeCapacity = std::max(maxHyperedgeCapacity, builder.hyperedges[i].capacity);
                }

		    });
		    hyperedges.back().first_out = PinIndex(numPins());
            numPinsAtHyperedgeStart = numPins();

            for (size_t sourceOcc : sourceOccurrences) {
                nodes[globalToLocalID[source]+1].first_out += sourceOcc;
            }
            for (size_t targetOcc : targetOccurrences) {
                nodes[globalToLocalID[target]+1].first_out += targetOcc;
            }
		}

		void addMockBuilder(MockBuilder& builder, bool add_nodes) {
		    if ( !builder.finishHyperedge() ) {
		        builder.hyperedges.back().capacity = 0;
		    }
		    if ( !finishHyperedge() ) {
		        hyperedges.back().capacity = 0;
		    }

		    size_t numNodesBefore = numNodes();

		    if (add_nodes) {
		        // remove sentinel
		        NodeData node_data = nodes.back();
		        nodes.pop_back();
		        builder.nodes[0].first_out = node_data.first_out;

                for (NodeData u : builder.nodes) {
                    nodes.push_back(u);
                }
		    }

            PinIndex first_out = hyperedges.back().first_out;
            hyperedges.pop_back();

            if (pins_receiving_flow.size() > 0) {
                pins_receiving_flow.back() = PinIndexRange(first_out + builder.hyperedges[0].first_out, first_out + builder.hyperedges[0].first_out);
            }

            hyperedges.push_back({first_out + builder.hyperedges[0].first_out, Flow(0), builder.hyperedges[0].capacity});
            maxHyperedgeCapacity = std::max(maxHyperedgeCapacity, builder.hyperedges[0].capacity);

            for (size_t i = 1; i < builder.hyperedges.size(); ++i) {
                pins_sending_flow.emplace_back(hyperedges.back().first_out, hyperedges.back().first_out);
                hyperedges.push_back({first_out + builder.hyperedges[i].first_out, Flow(0), builder.hyperedges[i].capacity});
                pins_receiving_flow.emplace_back(hyperedges.back().first_out, hyperedges.back().first_out);
                maxHyperedgeCapacity = std::max(maxHyperedgeCapacity, builder.hyperedges[i].capacity);
            }

            for (Node u : builder.pins) {
                if (add_nodes) {
                    pins.push_back({u + Node(numNodesBefore), InHeIndex::Invalid()});
                } else {
                    pins.push_back({u, InHeIndex::Invalid()});
                }
            }

            numPinsAtHyperedgeStart = numPins();
            assert(hyperedges.back().first_out == numPins());
		}
		
		void finalize() {
			if( !finishHyperedge() )	//finish last open hyperedge
				hyperedges.back().capacity = 0;	//maybe the last started hyperedge has zero/one pins and thus we still use the previous sentinel. was never a bug, since that capacity is never read
			
			total_node_weight = NodeWeight(0);
			for (Node u : nodeIDs()) {
				nodes[u+1].first_out += nodes[u].first_out;
				total_node_weight += nodes[u].weight;
			}

			incident_hyperedges.resize(numPins());
			tbb::parallel_for(tbb::blocked_range<size_t>(0, numHyperedges()), [&](const tbb::blocked_range<size_t>& range) {
			    for (size_t i = range.begin(); i < range.end(); ++i) {
			        const Hyperedge e = static_cast<Hyperedge>(i);
                    assert(beginIndexPins(e) < endIndexPins(e));
                    for (auto pin_it = beginIndexPins(e); pin_it != endIndexPins(e); pin_it++) {
                        Pin& p = pins[pin_it];
                        InHeIndex ind_he = static_cast<InHeIndex>(__sync_fetch_and_add(&nodes[p.pin].first_out.value(), 1));	//destroy first_out temporarily and reset later
                        incident_hyperedges[ind_he] = { e, Flow(0), pin_it };
                        p.he_inc_iter = ind_he;					//set iterator for incident hyperedge -> its position in incident_hyperedges of the node
                    }
			    }
			});
			
			for (Node u(numNodes()-1); u > 0; u--)
				nodes[u].first_out = nodes[u-1].first_out;	//reset temporarily destroyed first_out
			nodes[0].first_out = InHeIndex(0);
			
			finalized = true;
		}
		
		void shrink_to_fit() {
			nodes.shrink_to_fit();
			hyperedges.shrink_to_fit();
			pins.shrink_to_fit();
			incident_hyperedges.shrink_to_fit();
			pins_sending_flow.shrink_to_fit();
			pins_receiving_flow.shrink_to_fit();
		}

        size_t numPinsAtHyperedgeStart = 0;

	private:
		
		void removeLastPin() {
			nodes[ pins.back().pin + 1 ].first_out--;
			pins.pop_back();
		}
		
		bool finishHyperedge() {
			if (currentHyperedgeSize() == 1) {
				removeLastPin();
			}
			
			if (currentHyperedgeSize() > 0) {
				pins_sending_flow.emplace_back(hyperedges.back().first_out, hyperedges.back().first_out);
				hyperedges.push_back({PinIndex::fromOtherValueType(numPins()), Flow(0), Flow(0)});//sentinel
				pins_receiving_flow.emplace_back(hyperedges.back().first_out, hyperedges.back().first_out);
				return true;
			}
			return false;
		}
		
		bool finalized = false;
	};
}