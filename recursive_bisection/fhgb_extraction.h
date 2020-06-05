#pragma once

#include "../datastructure/flow_hypergraph_builder.h"
#include "hypergraph.h"
#include "../datastructure/queue.h"
#include "partition_ca.h"
#include "../datastructure/node_border.h"
#include "config.h"
#include "partition_threadsafe.h"
#include <tbb/scalable_allocator.h>

namespace whfc_rb {
    template<class Hypergraph>
    class FlowHypergraphBuilderExtractor {
    public:
        template<typename T>
        using vec = std::vector<T, tbb::scalable_allocator<T>>;

        static constexpr NodeID invalid_node = std::numeric_limits<NodeID>::max();
        Hypergraph fhgb;

        static constexpr bool doBothBFSInParallel = false;

        FlowHypergraphBuilderExtractor(const size_t maxNumNodes, const size_t maxNumEdges, const size_t maxNumPins, int seed, const PartitionConfig& config) :
                fhgb(2 * config.percentage_bfs_from_cut * maxNumNodes + 2, maxNumEdges),
                visitedNode(maxNumNodes), visitedHyperedge(maxNumEdges),
                globalToLocalID(maxNumNodes),
                localToGlobalID(maxNumNodes + 2),
                mt(seed), config(config),
                queue_thread_specific_1(maxNumNodes), queue_thread_specific_2(maxNumNodes),
                tryVisitNodeLock(maxNumNodes) //TODO: this has to be replaced
                {}

        struct ExtractorInfo {
            whfc::Node source;
            whfc::Node target;
            whfc::Flow baseCut;
            whfc::Flow cutAtStake;
        };

        template<class PartitionImpl>
        ExtractorInfo
        run(PartitionImpl &partition, const PartitionBase::PartitionID part0, const PartitionBase::PartitionID part1,
            NodeWeight maxW0, NodeWeight maxW1, whfc::DistanceFromCut& distanceFromCut, whfc::TimeReporter& timer) {

            timer.start("Init", "Extraction");
            CSRHypergraph &hg = partition.getGraph();
            initialize(hg.numNodes(), hg.numHyperedges());
            timer.stop("Init");

            whfc::HopDistance delta = config.distancePiercing ? 1 : 0;

            timer.start("Filter", "Extraction");
            auto cut_hes = partition.getCutEdges(part0, part1);
            timer.stop("Filter");

            timer.start("Shuffle", "Extraction");
            cut_hes.shuffle(mt);
            timer.stop("Shuffle");

            timer.start("BFS", "Extraction");

            std::atomic<NodeWeight> w0 = 0;
            std::atomic<NodeWeight> w1 = 0;

            timer.start("Process_cut_hyperedges", "BFS");
            // This fills the first layers for the bfs
            processCutHyperedges(hg, cut_hes, partition, part0, part1, maxW0, maxW1, w0, w1);
            timer.stop("Process_cut_hyperedges");

            timer.start("Scan_hyperedges_and_add_nodes", "BFS");
            tbb::parallel_invoke([&]() {
                BreadthFirstSearch(hg, cut_hes, partition, part0, part1, maxW0, queue_thread_specific_1, w0);
            }, [&]() {
                BreadthFirstSearch(hg, cut_hes, partition, part1, part0, maxW1, queue_thread_specific_2, w1);
            });

            timer.stop("Scan_hyperedges_and_add_nodes");

            timer.start("writeQueuesToBuilder", "BFS");
            writeQueuesToBuilder(fhgb, queue_thread_specific_1, hg, distanceFromCut, -delta);
            result.source = whfc::Node(fhgb.numNodes());
            localToGlobalID[result.source] = invalid_node;
            fhgb.addNode(whfc::NodeWeight(0));
            writeQueuesToBuilder(fhgb, queue_thread_specific_2, hg, distanceFromCut, delta);
            result.target = whfc::Node(fhgb.numNodes());
            localToGlobalID[result.target] = invalid_node;
            fhgb.addNode(whfc::NodeWeight(0));
            timer.stop("writeQueuesToBuilder");

            visitedHyperedge.reset();

            timer.start("Copy_hyperedges", "BFS");
            size_t total_num_hyperedges = 0;
            for (vec<HyperedgeWithSize>& hyperedges_local : hyperedges_thread_specific) {
                total_num_hyperedges += hyperedges_local.size();
            }

            hyperedges.reserve(total_num_hyperedges);
            for (vec<HyperedgeWithSize>& hyperedges_local : hyperedges_thread_specific) {
                hyperedges.insert(hyperedges.end(), hyperedges_local.begin(), hyperedges_local.end());
                hyperedges_local.clear();
            }
            timer.stop("Copy_hyperedges");

            timer.start("Sort_hyperedges", "BFS");
            tbb::parallel_sort(hyperedges.begin(), hyperedges.end(), [](const HyperedgeWithSize& e0, const HyperedgeWithSize& e1) {
                return e0.e < e1.e;
            });
            timer.stop("Sort_hyperedges");

            timer.start("Compute_prefix_sum", "BFS");
            computePrefixSum(hyperedges);
            timer.stop("Compute_prefix_sum");

            timer.start("Add_hyperedges", "BFS");
            addHyperedges(hg, partition, part0, part1, hyperedges);
            timer.stop("Add_hyperedges");

            for (whfc::Flow& baseCut : baseCut_thread_specific) {
                result.baseCut += baseCut;
                baseCut = 0;
            }

            for (whfc::Flow& cutAtStake : cutAtStake_thread_specific) {
                result.cutAtStake += cutAtStake;
                cutAtStake = 0;
            }

            timer.stop("BFS");

            fhgb.nodeWeight(result.source) = partition.partWeight(part0) - w0;
            fhgb.nodeWeight(result.target) = partition.partWeight(part1) - w1;

            timer.start("Finalize", "Extraction");
            fhgb.finalize();
            timer.stop("Finalize");

            //fhgb.printHypergraph(std::cout);

            return result;
        }

        auto localNodeIDs() const {
            return boost::irange<whfc::Node>(whfc::Node(0), whfc::Node(fhgb.numNodes()));
        }

        whfc::Node global2local(const NodeID x) const {
            assert(visitedNode.isSet(x));
            return globalToLocalID[x];
        }

        NodeID local2global(const whfc::Node x) const { return localToGlobalID[x]; }

    private:
        struct HyperedgeWithSize {
            HyperedgeID e;
            size_t pin_count;
        };

        ldc::AtomicTimestampSet<> visitedNode;
        ldc::AtomicTimestampSet<> visitedHyperedge;
        std::vector<whfc::Node> globalToLocalID;
        std::vector<NodeID, tbb::cache_aligned_allocator<NodeID>> localToGlobalID;
        ExtractorInfo result;
        std::mt19937 mt;
        const PartitionConfig& config;

        tbb::enumerable_thread_specific<whfc::Flow> baseCut_thread_specific;
        tbb::enumerable_thread_specific<whfc::Flow> cutAtStake_thread_specific;

        tbb::enumerable_thread_specific<LayeredQueue<whfc::Node>> queue_thread_specific_1;
        tbb::enumerable_thread_specific<LayeredQueue<whfc::Node>> queue_thread_specific_2;

        tbb::enumerable_thread_specific<vec<HyperedgeWithSize>> hyperedges_thread_specific;

        vec<HyperedgeWithSize> hyperedges;

        std::vector<std::mutex> tryVisitNodeLock;

        void computePrefixSum(vec<HyperedgeWithSize>& vector) {
            tbb::parallel_scan(tbb::blocked_range<size_t>(0, vector.size()), 0,
                [&](const tbb::blocked_range<size_t>& r, size_t sum, bool is_final_scan) -> size_t {
                size_t temp = sum;
                for (size_t i = r.begin(); i < r.end(); ++i) {
                    temp = temp + vector[i].pin_count;
                    if(is_final_scan)
                        vector[i].pin_count = temp;
                }
                return temp;
            }, [] (size_t left, size_t right) {
                return left + right;
            });
        }

        template<class Builder>
        void writeQueuesToBuilder(Builder& builder, tbb::enumerable_thread_specific<LayeredQueue<whfc::Node>>& queues,
                CSRHypergraph& hg, whfc::DistanceFromCut& distanceFromCut, whfc::HopDistance delta) {
            struct WorkElement {
                LayeredQueue<whfc::Node>& queue;
                size_t layer;
                size_t indexStart;
            };

            std::vector<whfc::FlowHypergraph::NodeData>& nodes = builder.getNodes();
            if (queues.empty()) return;
            const size_t numLayers = queues.begin()->numLayers();
            std::vector<size_t> layer_start_index(numLayers + 1, 0);
            layer_start_index[0] = builder.numNodes();
            std::vector<WorkElement> work_elements;

            for (LayeredQueue<whfc::Node>& queue : queues) {
                for (size_t i = 0; i < numLayers; ++i) {
                    layer_start_index[i+1] += queue.layerSize(i);
                }
            }

            for (size_t i = 0; i < numLayers; ++i) {
                layer_start_index[i+1] += layer_start_index[i];
            }

            for (LayeredQueue<whfc::Node>& queue : queues) {
                for (size_t i = 0; i < numLayers; ++i) {
                    if (layer_start_index[i+1] - layer_start_index[i] > 0) work_elements.push_back({queue, i, layer_start_index[i]});
                    layer_start_index[i] += queue.layerSize(i);
                }
            }

            // Write node data into builder
            nodes.resize(layer_start_index.back() + 1);

            tbb::parallel_for_each(work_elements, [&](WorkElement& workElement) {
                const size_t layerStart = workElement.queue.layerStart(workElement.layer);
                tbb::parallel_for(tbb::blocked_range<size_t>(0, workElement.queue.layerSize(workElement.layer)), [&](const tbb::blocked_range<size_t>& indices) {
                    for (size_t i = indices.begin(); i < indices.end(); ++i) {
                        const whfc::Node localID(workElement.indexStart + i);
                        const whfc::Node globalNode = workElement.queue.elementAt(layerStart + i);
                        nodes[localID].weight = whfc::NodeWeight(hg.nodeWeight(globalNode));
                        globalToLocalID[globalNode] = localID;
                        localToGlobalID[localID] = globalNode;
                        distanceFromCut[localID] = (workElement.layer + 1) * delta;
                    }
                });
            });
        }

        inline bool visitNode(const CSRHypergraph& hg, const NodeID u,
                std::atomic<whfc::NodeWeight>& w, NodeWeight& lastSeenValue, const NodeWeight maxWeight, LayeredQueue<whfc::Node>& queue) {
            std::lock_guard<std::mutex> guard(tryVisitNodeLock[u]);
            lastSeenValue = w.load(std::memory_order_relaxed);
            if (lastSeenValue + hg.nodeWeight(u) <= maxWeight) {
                if (visitedNode.set(u)) {
                    w += hg.nodeWeight(u);
                    queue.push(static_cast<whfc::Node>(u));
                }
                return true;
            }
            return visitedNode.isSet(u);
        }

        template<class PartitionImpl, typename CutEdgeRange>
        void processCutHyperedges(CSRHypergraph& hg, CutEdgeRange &cut_hes, const PartitionImpl& partition,
                PartitionBase::PartitionID partID, PartitionBase::PartitionID otherPartID,
                NodeWeight maxWeight1, NodeWeight maxWeight2, std::atomic<NodeWeight>& w1,
                std::atomic<NodeWeight>& w2) {
            tbb::blocked_range<size_t> range(0, cut_hes.size());
            tbb::parallel_for(range, [&](const tbb::blocked_range<size_t>& sub_range) {
                NodeWeight lastSeenValue1 = 0;
                NodeWeight lastSeenValue2 = 0;
                LayeredQueue<whfc::Node>& queue1 = queue_thread_specific_1.local();
                LayeredQueue<whfc::Node>& queue2 = queue_thread_specific_2.local();
                vec<HyperedgeWithSize>& hyperedges = hyperedges_thread_specific.local();
                whfc::Flow &baseCut = baseCut_thread_specific.local();
                whfc::Flow &cutAtStake = cutAtStake_thread_specific.local();
                for (size_t i = sub_range.begin(); i < sub_range.end(); ++i) {
                    const HyperedgeID e = cut_hes[i];
                    cutAtStake += hg.hyperedgeWeight(e);
                    size_t num_pins = 0;
                    bool connectToSource = false;
                    bool connectToTarget = false;
                    for (NodeID u : hg.pinsOf(e)) {
                        if (partition[u] == partID) {
                            if (visitNode(hg, u, w1, lastSeenValue1, maxWeight1, queue1)) {
                                num_pins++;
                            } else {
                                connectToSource = true;
                            }
                        } else if (partition[u] == otherPartID) {
                            if (visitNode(hg, u, w2, lastSeenValue2, maxWeight2, queue2)) {
                                num_pins++;
                            } else {
                                connectToTarget = true;
                            }
                        }
                    }
                    if (!(connectToSource && connectToTarget)) {
                        if (connectToSource || connectToTarget) num_pins++;
                        hyperedges.push_back({e, num_pins});
                    } else {
                        baseCut += hg.hyperedgeWeight(e);
                    }
                    //if (lastSeenValue1 >= maxWeight1 && lastSeenValue2 >= maxWeight2) break;
                }
            });

        }

        template<class PartitionImpl, typename CutEdgeRange>
        void BreadthFirstSearch(CSRHypergraph &hg, CutEdgeRange &cut_hes, const PartitionImpl &partition,
                                            PartitionBase::PartitionID partID, PartitionBase::PartitionID otherPartID,
                                            NodeWeight maxWeight,
                                            tbb::enumerable_thread_specific<LayeredQueue<whfc::Node>>& queue_thread_specific,
                                            std::atomic<NodeWeight>& w) {
            bool nodes_left = false;

            for (LayeredQueue<whfc::Node>& queue : queue_thread_specific) {
                queue.finishNextLayer();
                if (!queue.currentLayerEmpty()) {
                    nodes_left = true;
                }
            }

            size_t numLayers = 1;

            while (nodes_left) {
                // scan hyperedges and add nodes to the next layer
                tbb::parallel_for_each(queue_thread_specific, [&](LayeredQueue<whfc::Node>& queue) {
                    tbb::parallel_for(tbb::blocked_range<size_t>(queue.currentLayerStart(), queue.currentLayerEnd()), [&](const tbb::blocked_range<size_t>& indices) {
                        LayeredQueue<whfc::Node>& local_queue = queue_thread_specific.local();
                        vec<HyperedgeWithSize>& hyperedges = hyperedges_thread_specific.local();
                        NodeWeight lastSeenValue = 0;

                        if (local_queue.numLayers() == 0) {
                            // This is a new queue that joined later
                            for (size_t i = 0; i < numLayers; ++i) {
                                local_queue.finishNextLayer();
                            }
                        }

                        for (size_t i = indices.begin(); i < indices.end(); ++i) {
                            const NodeID u = queue.elementAt(i);

                            //if (lastSeenValue >= maxWeight) break;
                            for (HyperedgeID e : hg.hyperedgesOf(u)) {
                                if (partition.pinsInPart(otherPartID, e) == 0 && partition.pinsInPart(partID, e) > 1 && visitedHyperedge.set(e)) {
                                    size_t num_pins = 0;
                                    bool connectToTerminal = false;
                                    for (NodeID v : hg.pinsOf(e)) {
                                        if (partition[v] == partID) {
                                            if (visitNode(hg, v, w, lastSeenValue, maxWeight, local_queue)) {
                                                num_pins++;
                                            } else {
                                                connectToTerminal = true;
                                            }
                                        }
                                    }
                                    if (connectToTerminal) num_pins++;
                                    hyperedges.push_back({e, num_pins});
                                }
                            }
                        }
                    });
                });

                nodes_left = false;

                for (LayeredQueue<whfc::Node>& queue : queue_thread_specific) {
                    queue.clearCurrentLayer();
                    queue.finishNextLayer();
                    if (queue.currentLayer().size() > 0) {
                        nodes_left = true;
                    }
                }
                numLayers++;
            }

            for (LayeredQueue<whfc::Node>& queue : queue_thread_specific) {
                queue.removeLastLayerBound();
            }
        }

        using PinIndexRange = mutable_index_range<whfc::PinIndex>;

        template<typename PartitionImpl>
        void addHyperedges(CSRHypergraph& hg, const PartitionImpl &partition, PartitionBase::PartitionID part0,
                PartitionBase::PartitionID part1, vec<HyperedgeWithSize>& hyperedges) {
            tbb::enumerable_thread_specific<size_t> sourceOccurrences(0);
            tbb::enumerable_thread_specific<size_t> targetOccurrences(0);

            if (hyperedges.size() == 0) return;
            fhgb.hyperedges.resize(hyperedges.size() + 1);
            fhgb.pins_sending_flow.resize(hyperedges.size());
            fhgb.pins_receiving_flow.resize(hyperedges.size());
            fhgb.pins.resize(hyperedges[hyperedges.size() - 1].pin_count);

            fhgb.hyperedges.back().first_out = whfc::PinIndex(hyperedges[hyperedges.size() - 1].pin_count);
            fhgb.numPinsAtHyperedgeStart = hyperedges[hyperedges.size() - 1].pin_count;

            tbb::parallel_for(tbb::blocked_range<size_t>(0, hyperedges.size()), [&](const tbb::blocked_range<size_t>& indices) {
                for (size_t i = indices.begin(); i < indices.end(); ++i) {
                    const HyperedgeWithSize edge = hyperedges[i];
                    const HyperedgeID e = edge.e;
                    bool connectToSource = false;
                    bool connectToTarget = false;
                    size_t nextPinPosition = i == 0 ? 0 : hyperedges[i-1].pin_count;

                    size_t& sourceOcc = sourceOccurrences.local();
                    size_t& targetOcc = targetOccurrences.local();

                    fhgb.pins_sending_flow[i] = PinIndexRange(whfc::PinIndex(nextPinPosition), whfc::PinIndex(nextPinPosition));
                    fhgb.hyperedges[i] = {whfc::PinIndex(nextPinPosition), whfc::Flow(0), whfc::Flow(hg.hyperedgeWeight(e))};
                    fhgb.pins_receiving_flow[i] = PinIndexRange(whfc::PinIndex(hyperedges[i].pin_count), whfc::PinIndex(hyperedges[i].pin_count));

                    if (partition.pinsInPart(part0, e) > 0 && partition.pinsInPart(part1, e) > 0) {
                        // Cut hyperedge
                        for (NodeID v : hg.pinsOf(e)) {
                            if (visitedNode.isSet(v)) {
                                fhgb.pins[nextPinPosition++] = {globalToLocalID[v], whfc::InHeIndex::Invalid()};
                                __atomic_fetch_add(&fhgb.nodes[globalToLocalID[v]+1].first_out.value(), 1, __ATOMIC_RELAXED);
                            } else {
                                connectToSource |= (partition[v] == part0);
                                connectToTarget |= (partition[v] == part1);
                                assert(!(connectToSource && connectToTarget));
                            }
                        }
                    } else {
                        const PartitionBase::PartitionID partID = partition.pinsInPart(part0, e) > 0 ? part0 : part1;
                        const whfc::Node terminal = partID == part0 ? result.source : result.target;
                        for (NodeID v : hg.pinsOf(e)) {
                            if (partition[v] == partID) {
                                if (visitedNode.isSet(v)) {
                                    fhgb.pins[nextPinPosition++] = {globalToLocalID[v], whfc::InHeIndex::Invalid()};
                                    __atomic_fetch_add(&fhgb.nodes[globalToLocalID[v]+1].first_out.value(), 1, __ATOMIC_RELAXED);
                                } else {
                                    if (terminal == result.source) {
                                        connectToSource = true;
                                    } else {
                                        connectToTarget = true;
                                    }
                                }
                            }
                        }
                    }

                    if (connectToSource) {
                        fhgb.pins[nextPinPosition++] = {result.source, whfc::InHeIndex::Invalid()};
                        sourceOcc++;
                    }
                    if (connectToTarget) {
                        fhgb.pins[nextPinPosition++] = {result.target, whfc::InHeIndex::Invalid()};
                        targetOcc++;
                    }
                    assert(nextPinPosition == hyperedges[i].pin_count);
                }
            });
            for (size_t sourceOcc : sourceOccurrences) {
                fhgb.nodes[result.source+1].first_out += sourceOcc;
            }
            for (size_t targetOcc : targetOccurrences) {
                fhgb.nodes[result.target+1].first_out += targetOcc;
            }
        }

        void initialize(uint numNodes, uint numHyperedges) {
            fhgb.clear();
            hyperedges.clear();
            visitedNode.reset();
            visitedHyperedge.reset();
            result = {whfc::Node(0), whfc::Node(0), 0, 0};

            for (LayeredQueue<whfc::Node>& queue : queue_thread_specific_1) {
                queue.clear();
            }
            for (LayeredQueue<whfc::Node>& queue : queue_thread_specific_2) {
                queue.clear();
            }
        }
    };
}