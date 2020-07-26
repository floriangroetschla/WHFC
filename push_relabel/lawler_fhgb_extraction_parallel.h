#pragma once

#include "../push_relabel/flow_hypergraph_builder.h"
#include "../datastructure/hypergraph.h"
#include "../datastructure/queue.h"
#include "../datastructure/node_border.h"
#include "../partitioner/config.h"
#include <tbb/scalable_allocator.h>
#include "../partitioner/fhgb_extraction_parallel_base.h"
#include <mutex>
#include <tbb/parallel_invoke.h>
#include <tbb/parallel_sort.h>
#include <tbb/parallel_scan.h>

namespace whfc_pr {
    template<class Hypergraph, class PartitionImpl>
    class LawlerFlowHypergraphBuilderExtractorParallel : public whfc_rb::ExtractorParallelBase<Hypergraph, PartitionImpl> {
    public:
        using Base = whfc_rb::ExtractorParallelBase<Hypergraph, PartitionImpl>;

        LawlerFlowHypergraphBuilderExtractorParallel(const size_t maxNumNodes, const size_t maxNumEdges, const size_t maxNumPins, int seed, const whfc_rb::PartitionerConfig& config) :
            Base(maxNumNodes, maxNumEdges, maxNumPins, seed, config)
        {}

    private:
        virtual void addHyperedges(whfc_rb::CSRHypergraph& hg, const PartitionImpl &partition, whfc_rb::PartitionBase::PartitionID part0,
                           whfc_rb::PartitionBase::PartitionID part1, typename Base::template vec<typename Base::HyperedgeWithSize>& hyperedges) {
            tbb::enumerable_thread_specific<size_t> sourceOccurrences(0);
            tbb::enumerable_thread_specific<size_t> targetOccurrences(0);

            if (hyperedges.size() == 0) return;
            Base::fhgb.hyperedges.resize(hyperedges.size() + 1);
            Base::fhgb.pins_in_sending_flow_end.resize(hyperedges.size());
            Base::fhgb.pins_out_receiving_flow_end.resize(hyperedges.size());
            Base::fhgb.pins_in.resize(hyperedges[hyperedges.size() - 1].pin_count);
            Base::fhgb.pins_out.resize(hyperedges[hyperedges.size() - 1].pin_count);

            Base::fhgb.hyperedges.back().first_out = whfc::PinIndex(hyperedges[hyperedges.size() - 1].pin_count);
            Base::fhgb.numPinsAtHyperedgeStart = hyperedges[hyperedges.size() - 1].pin_count;

            tbb::parallel_for(tbb::blocked_range<size_t>(0, hyperedges.size()), [&](const tbb::blocked_range<size_t>& indices) {
                for (size_t i = indices.begin(); i < indices.end(); ++i) {
                    const typename Base::HyperedgeWithSize edge = hyperedges[i];
                    const typename whfc_rb::HyperedgeID e = edge.e;
                    bool connectToSource = false;
                    bool connectToTarget = false;
                    size_t nextPinPosition = i == 0 ? 0 : hyperedges[i-1].pin_count;

                    size_t& sourceOcc = sourceOccurrences.local();
                    size_t& targetOcc = targetOccurrences.local();

                    Base::fhgb.pins_in_sending_flow_end[i] = whfc::PinIndex(nextPinPosition);
                    Base::fhgb.pins_out_receiving_flow_end[i] = whfc::PinIndex(nextPinPosition);
                    Base::fhgb.hyperedges[i] = {whfc::PinIndex(nextPinPosition), whfc::Flow(0), whfc::Flow(hg.hyperedgeWeight(e))};

                    if (partition.pinsInPart(part0, e) > 0 && partition.pinsInPart(part1, e) > 0) {
                        // Cut hyperedge
                        for (whfc_rb::NodeID v : hg.pinsOf(e)) {
                            if (Base::visitedNode.isSet(v)) {
                                Base::fhgb.pins_in[nextPinPosition] = {Base::globalToLocalID[v], whfc::InHeIndex::Invalid()};
                                Base::fhgb.pins_out[nextPinPosition++] = {Base::globalToLocalID[v], whfc::InHeIndex::Invalid()};
                                __atomic_fetch_add(&Base::fhgb.nodes[Base::globalToLocalID[v]+1].first_out.value(), 1, __ATOMIC_RELAXED);
                            } else {
                                connectToSource |= (partition[v] == part0);
                                connectToTarget |= (partition[v] == part1);
                                assert(!(connectToSource && connectToTarget));
                            }
                        }
                    } else {
                        const whfc_rb::PartitionBase::PartitionID partID = partition.pinsInPart(part0, e) > 0 ? part0 : part1;
                        const whfc::Node terminal = partID == part0 ? Base::result.source : Base::result.target;
                        for (whfc_rb::NodeID v : hg.pinsOf(e)) {
                            if (partition[v] == partID) {
                                if (Base::visitedNode.isSet(v)) {
                                    Base::fhgb.pins_in[nextPinPosition] = {Base::globalToLocalID[v], whfc::InHeIndex::Invalid()};
                                    Base::fhgb.pins_out[nextPinPosition++] = {Base::globalToLocalID[v], whfc::InHeIndex::Invalid()};
                                    __atomic_fetch_add(&Base::fhgb.nodes[Base::globalToLocalID[v]+1].first_out.value(), 1, __ATOMIC_RELAXED);
                                } else {
                                    if (terminal == Base::result.source) {
                                        connectToSource = true;
                                    } else {
                                        connectToTarget = true;
                                    }
                                }
                            }
                        }
                    }

                    if (connectToSource) {
                        Base::fhgb.pins_in[nextPinPosition] = {Base::result.source, whfc::InHeIndex::Invalid()};
                        Base::fhgb.pins_out[nextPinPosition++] = {Base::result.source, whfc::InHeIndex::Invalid()};
                        sourceOcc++;
                    }
                    if (connectToTarget) {
                        Base::fhgb.pins_in[nextPinPosition] = {Base::result.target, whfc::InHeIndex::Invalid()};
                        Base::fhgb.pins_out[nextPinPosition++] = {Base::result.target, whfc::InHeIndex::Invalid()};
                        targetOcc++;
                    }
                    assert(nextPinPosition == hyperedges[i].pin_count);
                }
            });
            for (size_t sourceOcc : sourceOccurrences) {
                Base::fhgb.nodes[Base::result.source+1].first_out += sourceOcc;
            }
            for (size_t targetOcc : targetOccurrences) {
                Base::fhgb.nodes[Base::result.target+1].first_out += targetOcc;
            }
        }
    };
}