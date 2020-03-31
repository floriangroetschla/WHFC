#pragma once

#include <boost/dynamic_bitset.hpp>
#include <atomic>
#include <tbb/parallel_do.h>
#include <mutex>
#include "whfc_refiner_two_way.h"
#include "partition_threadsafe.h"

namespace whfc_rb {
    class KWayRefinerParallel {
    public:
        using PartitionID = PartitionBase::PartitionID;
        
        explicit KWayRefinerParallel(PartitionThreadsafe &partition, whfc::TimeReporter &timer, std::mt19937 &mt, const PartitionConfig& config) :
                partition(partition), partActive(partition.numParts()), partActiveNextRound(partition.numParts()), block_pair_status(partition.numParts() * (partition.numParts() - 1) / 2), partScheduled(partition.numParts()), timer(timer), mt(mt), config(config) {}

        uint refine(double epsilon, uint maxIterations, int seed) {

            NodeWeight maxWeight = (1.0 + epsilon) * partition.totalWeight() / static_cast<double>(partition.numParts());
            partActive.assign(partActive.size(), 1);
            partActiveNextRound.assign(partActive.size(), 0);
            iterationCounter = 0;

            whfc::TimeReporter timer_dummy; // Timer not useful yet
            timer_dummy.active = true;
            tbb::enumerable_thread_specific<WHFCRefinerTwoWay> localRefiner(partition.getGraph().numNodes(), partition.getGraph().numHyperedges(), partition.getGraph().numPins(), seed, &timer_dummy);
            tbb::enumerable_thread_specific<whfc::TimeReporter> timer_local;

            while (std::any_of(partActive.begin(), partActive.end(), [](auto& x) { return x > 0; }) && iterationCounter < maxIterations) {
                std::vector<WorkElement> tasks = initialBlockPairs();

                tbb::parallel_do(tasks,
                        [&](WorkElement element, tbb::parallel_do_feeder<WorkElement>& feeder) {
                            assert(partScheduled[element.part0]);
                            assert(partScheduled[element.part1]);

                            if (iterationCounter < maxIterations) {
                                WHFCRefinerTwoWay& refiner = localRefiner.local();
                                whfc::TimeReporter& localTimer = timer_local.local();
                                refiner.setTimer(&localTimer);
                                localTimer.start("Refinement");
                                bool refinementResult = refiner.refine(partition, element.part0, element.part1, maxWeight, maxWeight, config);
                                localTimer.stop("Refinement");
                                if (refinementResult) {
                                    // Schedule for next round
                                    partActiveNextRound[element.part0] = 1;
                                    partActiveNextRound[element.part1] = 1;
                                }

                                blockPairStatus(element.part0, element.part1) = TaskStatus::FINISHED;

                                if (!addNewTasks(element.part0, feeder, maxWeight)) {
                                    partScheduled[element.part0] = false;
                                }
                                if (!addNewTasks(element.part1, feeder, maxWeight)) {
                                    partScheduled[element.part1] = false;
                                }
                                iterationCounter++;
                            }
                        }
                );

                //assert(allPairsProcessed()); only if maxIterations allows it
                std::swap(partActive, partActiveNextRound);
                partActiveNextRound.assign(partActive.size(), 0);
            }
            for (auto local_timer = timer_local.begin(); local_timer != timer_local.end(); local_timer++) {
                timer.merge(*local_timer, "Refinement", "Refinement");
            }

            return iterationCounter.load();
        }

    private:
        enum class TaskStatus {UNSCHEDULED, SCHEDULED, FINISHED};

        PartitionThreadsafe &partition;
        std::vector<uint8_t> partActive, partActiveNextRound;
        std::vector<std::atomic<TaskStatus>> block_pair_status;
        std::vector<std::atomic<bool>> partScheduled;
        whfc::TimeReporter &timer;
        std::mt19937 &mt;
        std::atomic<uint> iterationCounter = 0;
        const PartitionConfig& config;

        struct WorkElement {
            PartitionID part0;
            PartitionID part1;
        };

        size_t guessNumCutEdges(PartitionID part0, PartitionID part1) {
            return partition.cutEdges(part0, part1).size() + partition.cutEdges(part1, part0).size();
        }
        
        std::vector<WorkElement> initialBlockPairs() {
            resetBlockPairStatus();

            std::vector<size_t> participations(partition.numParts(), 0);
            std::vector<std::pair<PartitionID, PartitionID>> block_pairs;
            for (PartitionID i = 0; i < partition.numParts() - 1; ++i) {
                for (PartitionID j = i + 1; j < partition.numParts(); ++j) {
                    if (guessNumCutEdges(i, j) > 0) {
                        participations[i]++;
                        participations[j]++;
                        block_pairs.emplace_back(i,j);
                    }
                }
            }

            // schedule block pairs that want to participate in fewer calls first. the intuition is that a later stages
            // we still have the busy ones left
            std::sort(block_pairs.begin(), block_pairs.end(), [&](const auto& lhs, const auto& rhs) {
                return std::max(participations[lhs.first], participations[lhs.second])
                       < std::max(participations[rhs.first], participations[rhs.second]);
            });

            // maybe random shuffle works too
            // std::shuffle(block_pairs.begin(), block_pairs.end(), mt);

            std::vector<WorkElement> tasks;
            for (const auto& bp : block_pairs) {
                const PartitionID p0 = bp.first, p1 = bp.second;
                if (!partScheduled[p0] && !partScheduled[p1] && (partActive[p0] || partActive[p1])) {
                    tasks.push_back({p0, p1});
                    blockPairStatus(p0, p1) = TaskStatus::SCHEDULED;
                    partScheduled[p0] = true;
                    partScheduled[p1] = true;
                }
            }

            return tasks;
        }

        bool addNewTasks(PartitionID part, tbb::parallel_do_feeder<WorkElement>& feeder, NodeWeight maxWeight) {
            bool foundPair = false;
            for (PartitionID pid = 0; pid < partition.numParts(); ++pid) {
                if (pid != part && guessNumCutEdges(part, pid) > 0 && !partScheduled[pid]) {
                    if (!partScheduled[pid].exchange(true)) {
                        if (blockPairStatus(part, pid) == TaskStatus::UNSCHEDULED) {
                            blockPairStatus(part, pid) = TaskStatus::SCHEDULED;
                            feeder.add({part, pid});
                            foundPair = true;
                            break;
                        } else {
                            partScheduled[pid] = false;
                        }
                    }
                }
            }
            return foundPair;
        }

        bool allPairsProcessed() {
            for (uint i = 0; i < block_pair_status.size(); ++i) {
                if (block_pair_status[i] != TaskStatus::FINISHED) return false;
            }
            return true;
        }

        void resetBlockPairStatus() {
            for (uint i = 0; i < block_pair_status.size(); ++i) {
                block_pair_status[i] = TaskStatus::UNSCHEDULED;
            }
        }

        std::atomic<TaskStatus>& blockPairStatus(PartitionID part0, PartitionID part1) {
            if (part0 < part1) std::swap(part0, part1);
            return (block_pair_status[(part0 * (part0 - 1) / 2) + part1]);
        }

    };
}