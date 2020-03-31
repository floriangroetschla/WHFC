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
                partition(partition), partActive(partition.numParts()), partActiveNextRound(partition.numParts()), blockPairStatus(partition.numParts() * (partition.numParts() - 1) / 2), partitionScheduled(partition.numParts()), timer(timer), mt(mt), config(config) {}

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
                            assert(partitionScheduled[element.part0]);
                            assert(partitionScheduled[element.part1]);

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

                                blockPair(element.part0, element.part1) = TaskStatus::FINISHED;

                                if (!addNewTasks(element.part0, feeder, maxWeight)) {
                                    partitionScheduled[element.part0] = false;
                                }
                                if (!this->addNewTasks(element.part1, feeder, maxWeight)) {
                                    partitionScheduled[element.part1] = false;
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
        std::vector<std::atomic<TaskStatus>> blockPairStatus;
        std::vector<std::atomic<bool>> partitionScheduled;
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
            std::vector<WorkElement> tasks;
            PartitionID part0 = 0;
            PartitionID part1 = partition.numParts() - 1;
            resetBlockPairStatus();

            // TODO: count for each block in how many block pairs it wants to be (just take the initial guess)
            // TODO: only schedule block pairs that actually share a cut (partition.cutEdges.size() > 0)
            while (part0 < part1) {
                if (partActive[part0] || partActive[part1]) {
                    WorkElement element = {part0, part1};
                    tasks.push_back(element);
                    blockPair(part0, part1) = TaskStatus::SCHEDULED;
                    partitionScheduled[part0] = true;
                    partitionScheduled[part1] = true;
                }
                part0++;
                part1--;
            }
            return tasks;
        }

        bool addNewTasks(PartitionID part, tbb::parallel_do_feeder<WorkElement>& feeder, NodeWeight maxWeight) {
            bool foundPair = false;
            for (PartitionID pid = 0; pid < partition.numParts(); ++pid) {
                if (pid != part && !partitionScheduled[pid]) {
                    if (!partitionScheduled[pid].exchange(true)) {
                        if (blockPair(part, pid) == TaskStatus::UNSCHEDULED) {
                            blockPair(part, pid) = TaskStatus::SCHEDULED;
                            feeder.add({part, pid});
                            foundPair = true;
                            break;
                        } else {
                            partitionScheduled[pid] = false;
                        }
                    }
                }
            }
            return foundPair;
        }

        bool allPairsProcessed() {
            for (uint i = 0; i < blockPairStatus.size(); ++i) {
                if (blockPairStatus[i] != TaskStatus::FINISHED) return false;
            }
            return true;
        }

        void resetBlockPairStatus() {
            for (uint i = 0; i < blockPairStatus.size(); ++i) {
                blockPairStatus[i] = TaskStatus::UNSCHEDULED;
            }
        }

        std::atomic<TaskStatus>& blockPair(PartitionID part0, PartitionID part1) {
            if (part0 < part1) std::swap(part0, part1);
            return (blockPairStatus[(part0 * (part0 - 1) / 2) + part1]);
        }

    };
}