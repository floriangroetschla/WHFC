#pragma once

#include "../definitions.h"
#include "../datastructure/flow_hypergraph.h"

class MockBuilder {
public:
    using Flow = whfc::Flow;
    using Node = whfc::Node;
    using NodeData = whfc::FlowHypergraph::NodeData;
    using HyperedgeData = whfc::FlowHypergraph::HyperedgeData;
    using Pin = whfc::FlowHypergraph::Pin;
    using InHeIndex = whfc::InHeIndex;
    using PinIndex = whfc::PinIndex;

    std::vector<NodeData> nodes;
    std::vector<HyperedgeData> hyperedges;
    std::vector<Node> pins;

    size_t numPinsAtHyperedgeStart = 0;
    Flow maxHyperedgeCapacity = 0;

    MockBuilder() {
        clear();
    }

    void clear() {
        nodes.clear();
        hyperedges.clear();
        pins.clear();

        numPinsAtHyperedgeStart = 0;
        maxHyperedgeCapacity = 0;

        nodes.push_back({InHeIndex(0), whfc::NodeWeight(0)});
        hyperedges.push_back({PinIndex(0), Flow(0), Flow(0)});
    }

    void startHyperedge(const Flow capacity) {
        finishHyperedge();	//finish last hyperedge
        hyperedges.back().capacity = capacity;	//exploit sentinel
        numPinsAtHyperedgeStart = numPins();
        maxHyperedgeCapacity = std::max(maxHyperedgeCapacity, capacity);
    }

    void addPin(const Node u) {
        assert(u < numNodes());
        pins.push_back(u);
        nodes[u+1].first_out++;
    }

    void addNode(const whfc::NodeWeight w) {
        nodes.back().weight = w;
        nodes.push_back({InHeIndex(0), whfc::NodeWeight(0)});
    }

    inline size_t numPins() const { return pins.size(); }

    inline size_t numNodes() const { return nodes.size() - 1 ; }

    size_t currentHyperedgeSize() const {
        return numPins() - numPinsAtHyperedgeStart;
    }

    bool finishHyperedge() {
        if (currentHyperedgeSize() == 1) {
            removeLastPin();
        }

        if (currentHyperedgeSize() > 0) {
            //pins_sending_flow.emplace_back(hyperedges.back().first_out, hyperedges.back().first_out);
            hyperedges.push_back({PinIndex::fromOtherValueType(numPins()), Flow(0), Flow(0)});//sentinel
            //pins_receiving_flow.emplace_back(hyperedges.back().first_out, hyperedges.back().first_out);
            return true;
        }
        return false;
    }

private:
    void removeLastPin() {
        nodes[ pins.back() + 1 ].first_out--;
        pins.pop_back();
    }
};
