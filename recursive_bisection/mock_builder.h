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

    std::vector<NodeData>& nodes;
    std::vector<NodeData> nodes_internal;
    std::vector<HyperedgeData> hyperedges;
    std::vector<Node> pins;
    whfc::Node source = whfc::invalidNode;
    whfc::Node target = whfc::invalidNode;
    size_t sourceOccurences = 0;
    size_t targetOccurences = 0;

    size_t numPinsAtHyperedgeStart = 0;
    Flow maxHyperedgeCapacity = 0;

    MockBuilder() : nodes(nodes_internal) {
        clear();
    }

    MockBuilder(std::vector<NodeData>& existing_nodes, Node source, Node target) : nodes(existing_nodes), source(source), target(target) {
        clear(false);
    }

    void addTerminalOccurences() {
        nodes[source+1].first_out += sourceOccurences;
        nodes[target+1].first_out += targetOccurences;
        sourceOccurences = 0;
        targetOccurences = 0;
    }

    void setSource(Node new_source) {
        source = new_source;
    }

    void setTarget(Node new_target) {
        target = new_target;
    }

    void clear(bool clearNodes = true) {
        if (clearNodes) nodes.clear();
        hyperedges.clear();
        pins.clear();

        numPinsAtHyperedgeStart = 0;
        maxHyperedgeCapacity = 0;
        sourceOccurences = 0;
        targetOccurences = 0;

        if (clearNodes) nodes.push_back({InHeIndex(0), whfc::NodeWeight(0)});
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
        if (u != source && u != target) {
            __sync_fetch_and_add(&nodes[u+1].first_out.value(), 1);
        } else if (u == source) {
            sourceOccurences++;
        } else {
            targetOccurences++;
        }
    }

    void addNode(const whfc::NodeWeight w) {
        nodes.back().weight = w;
        nodes.push_back({InHeIndex(0), whfc::NodeWeight(0)});
    }

    inline size_t numPins() const { return pins.size(); }

    inline size_t numNodes() const { return nodes.size() - 1 ; }

    inline size_t numHyperedges() const { return hyperedges.size() - 1 ; }

    size_t currentHyperedgeSize() const {
        return numPins() - numPinsAtHyperedgeStart;
    }

    void removeCurrentHyperedge() {
        while (numPins() > numPinsAtHyperedgeStart)
            removeLastPin();
    }

    bool finishHyperedge() {
        if (currentHyperedgeSize() == 1) {
            removeLastPin();
        }

        if (currentHyperedgeSize() > 0) {
            hyperedges.push_back({PinIndex::fromOtherValueType(numPins()), Flow(0), Flow(0)});//sentinel
            return true;
        }
        return false;
    }

    std::vector<NodeData>& getNodes() {
        return nodes;
    }

private:
    void removeLastPin() {
        if (pins.back() != source && pins.back() != target) {
            __sync_fetch_and_sub(&nodes[pins.back() + 1].first_out.value(), 1);
        } else if (pins.back() == source) {
            sourceOccurences--;
        } else {
            targetOccurences--;
        }

        pins.pop_back();
    }
};
