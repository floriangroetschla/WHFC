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
                partition(partition), partActive(partition.numParts()), partActiveNextRound(partition.numParts()),
                improvement_history(partition.numParts() * (partition.numParts() -1) / 2, 0),
                block_pair_status(partition.numParts() * (partition.numParts() - 1) / 2),
                partScheduled(partition.numParts()),
                timer(timer), mt(mt), config(config) {}

        uint refine(double epsilon, uint maxIterations, int seed) {
            NodeWeight maxWeight = (1.0 + epsilon) * partition.totalWeight() / static_cast<double>(partition.numParts());
            partActive.assign(partActive.size(), 1);
            partActiveNextRound.assign(partActive.size(), 0);
            iterationCounter = 0;

            tbb::enumerable_thread_specific<WHFCRefinerTwoWay> localRefiner(partition.getGraph().numNodes(), partition.getGraph().numHyperedges(), partition.getGraph().numPins(), seed, std::ref(config));

            // TODO maybe also make a row restriction. or terminate if there aren't enough initial tasks. 
            while (std::any_of(partActive.begin(), partActive.end(), [](auto& x) { return x > 0; }) && iterationCounter < maxIterations) {
                std::vector<WorkElement> tasks = initialBlockPairs();
                std::cout << "Round " << round << " initial tasks: " << tasks.size() << std::endl;
                tbb::parallel_do(tasks,
                        [&](WorkElement element, tbb::parallel_do_feeder<WorkElement>& feeder) {
                            assert(partScheduled[element.part0]);
                            assert(partScheduled[element.part1]);

                            if (iterationCounter < maxIterations) {
                                WHFCRefinerTwoWay& refiner = localRefiner.local();
                                bool refinementResult = refiner.refine(partition, element.part0, element.part1, maxWeight, maxWeight);
                                if (refinementResult) {
                                    // Schedule for next round
                                    partActiveNextRound[element.part0] = 1;
                                    partActiveNextRound[element.part1] = 1;
                                    reportImprovement(element.part0, element.part1);
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

                std::cout << "WHFC refiner calls: " << iterationCounter << std::endl;

                //assert(allPairsProcessed()); only if maxIterations allows it
                std::swap(partActive, partActiveNextRound);
                partActiveNextRound.assign(partActive.size(), 0);
                round++;
            }

            for (WHFCRefinerTwoWay& refiner : localRefiner) {
                timer.merge(refiner.getTimer(), "Refinement", "WHFCRefinerTwoWay");
            }

            return iterationCounter.load();
        }

    private:
        enum class TaskStatus {UNSCHEDULED, SCHEDULED, FINISHED};

        size_t round = 0;
        PartitionThreadsafe &partition;
        std::vector<uint8_t> partActive, partActiveNextRound;
        std::vector<uint32_t> improvement_history;
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
                    if (isEligible(i, j)) {
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
            for (PartitionID pid = 0; pid < partition.numParts(); ++pid) {
                if (pid != part && isEligible(part, pid)) {
                    if (!partScheduled[pid].exchange(true)) {
                        if (blockPairStatus(part, pid) == TaskStatus::UNSCHEDULED) {
                            blockPairStatus(part, pid) = TaskStatus::SCHEDULED;
                            feeder.add({part, pid});
                            return true;
                        } else {
                            partScheduled[pid] = false;
                        }
                    }
                }
            }
            return false;
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

        bool isEligible(PartitionID part0, PartitionID part1) {
            if (part0 < part1) std::swap(part0, part1);
            return guessNumCutEdges(part0, part1) > 0
                   && blockPairStatus(part0, part1) == TaskStatus::UNSCHEDULED
                   && (round < 2 || improvement_history[(part0 * (part0 - 1)/2) + part1] > 0);
        }

        void reportImprovement(PartitionID part0, PartitionID part1) {
            if (part0 < part1) std::swap(part0, part1);
            improvement_history[(part0 * (part0 - 1)/2) + part1]++;
        }



    };
}