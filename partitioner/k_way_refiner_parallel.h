#pragma once

#include <boost/dynamic_bitset.hpp>
#include <atomic>
#include <tbb/parallel_do.h>
#include <mutex>
#include "whfc_refiner_two_way.h"
#include "../datastructure/partition_threadsafe.h"

namespace whfc_rb {
    template<class PartitionImpl, class HypergraphImpl, class FlowAlgo, class Extractor>
    class KWayRefinerParallel {
    public:
        using PartitionID = PartitionBase::PartitionID;

        explicit KWayRefinerParallel(PartitionImpl &partition, whfc::TimeReporter &timer, std::mt19937 &mt, const PartitionerConfig& config) :
                partition(partition), partActive(partition.numParts()), partActiveNextRound(partition.numParts()),
                improvement_history(partition.numParts() * (partition.numParts() -1) / 2, 0),
                block_pair_status(partition.numParts() * (partition.numParts() - 1) / 2),
                partScheduled(partition.numParts()),
                timer(timer), mt(mt), config(config),
                refiners_thread_specific(partition.getGraph().numNodes(), partition.getGraph().numHyperedges(), partition.getGraph().numPins(), mt(), std::ref(config)),
                timers_thread_specific(), participations(partition.numParts(), 0), partitionsSortedByParticipations()
        {

        }

        uint refine(double epsilon, uint maxIterations) {
            NodeWeight maxWeight = (1.0 + epsilon) * partition.totalWeight() / static_cast<double>(partition.numParts());
            partActive.assign(partActive.size(), 1);
            partActiveNextRound.assign(partActive.size(), 0);
            iterationCounter = 0;
            std::atomic<size_t> iterations_done = 0;

            while (std::any_of(partActive.begin(), partActive.end(), [&](auto& x) { return x > 0; }) && iterationCounter < maxIterations) {
                std::vector<WorkElement> tasks = initialBlockPairs();
                tbb::parallel_do(tasks,[&](WorkElement element, tbb::parallel_do_feeder<WorkElement>& feeder) {
                    assert(partScheduled[element.part0]);
                    assert(partScheduled[element.part1]);

                    if (iterationCounter.fetch_add(1) < maxIterations) {
                        WHFCRefinerTwoWay<PartitionImpl, HypergraphImpl, FlowAlgo, Extractor>& refiner = refiners_thread_specific.local();
                        whfc::TimeReporter& timer = timers_thread_specific.local();
                        timer.start("WHFCRefinerTwoWay");
                        bool refinementResult;
                        tbb::this_task_arena::isolate( [&] {
                            refinementResult = refiner.refine(partition, element.part0, element.part1,
                                                              maxWeight, maxWeight, timer);
                        });
                        timer.stop("WHFCRefinerTwoWay");
                        if (refinementResult) {
                            // Schedule for next round
                            partActiveNextRound[element.part0] = 1;
                            partActiveNextRound[element.part1] = 1;
                            reportImprovement(element.part0, element.part1);
                        }

                        blockPairStatus(element.part0, element.part1) = TaskStatus::FINISHED;

                        partScheduled[element.part0] = false;
                        partScheduled[element.part1] = false;
                        timer.start("addNewTasks");
                        size_t id_with_more_participations = participations[element.part0] > participations[element.part1] ? element.part0 : element.part1;
                        size_t id_with_less_participations = participations[element.part0] <= participations[element.part1] ? element.part0 : element.part1;
                        WorkElement newElement = addBlockPair(id_with_more_participations);
                        if (newElement.part0 != 0 || newElement.part1 != 0) { feeder.add(newElement); }
                        newElement = addBlockPair(id_with_less_participations);
                        if (newElement.part0 != 0 || newElement.part1 != 0) { feeder.add(newElement); }
                        timer.stop("addNewTasks");
                        iterations_done++;
                    }
                });

                std::swap(partActive, partActiveNextRound);
                partActiveNextRound.assign(partActive.size(), 0);
                round++;
            }

            for (whfc::TimeReporter& local_timer : timers_thread_specific) {
                timer.merge(local_timer, "Refinement", "total");
            }

            return iterations_done.load();
        }

    private:
        enum class TaskStatus {UNSCHEDULED, SCHEDULED, FINISHED};

        struct WorkElement {
            PartitionID part0;
            PartitionID part1;
        };

        size_t round = 0;
        PartitionThreadsafe &partition;
        std::vector<uint8_t> partActive, partActiveNextRound;
        std::vector<uint32_t> improvement_history;
        std::vector<std::atomic<TaskStatus>> block_pair_status;
        std::vector<std::atomic<bool>> partScheduled;
        whfc::TimeReporter &timer;
        std::mt19937 &mt;
        std::atomic<uint> iterationCounter = 0;
        const PartitionerConfig& config;
        tbb::enumerable_thread_specific<WHFCRefinerTwoWay<PartitionImpl, HypergraphImpl, FlowAlgo, Extractor>> refiners_thread_specific;
        tbb::enumerable_thread_specific<whfc::TimeReporter> timers_thread_specific;
        std::vector<size_t> participations;
        std::vector<size_t> partitionsSortedByParticipations;


        size_t guessNumCutEdges(PartitionID part0, PartitionID part1) {
            return partition.cutEdges(part0, part1).size() + partition.cutEdges(part1, part0).size();
        }

        std::vector<WorkElement> initialBlockPairs() {
            timer.start("initialBlockPairs", "Refinement");
            resetBlockPairStatus();
            std::fill(participations.begin(), participations.end(), 0);
            partitionsSortedByParticipations.clear();
            std::vector<WorkElement> tasks;

            for (PartitionID i = 0; i < partition.numParts() - 1; ++i) {
                for (PartitionID j = i + 1; j < partition.numParts(); ++j) {
                    if (isEligible(i, j)) {
                        participations[i]++;
                        participations[j]++;
                    }
                }
            }

            for (PartitionID i = 0; i < partition.numParts(); ++i) {
                if (participations[i] > 0) {
                    partitionsSortedByParticipations.push_back(i);
                }
            }

            auto compareParticipations = [&](size_t i, size_t j) { return participations[i] > participations[j]; };

            std::sort(partitionsSortedByParticipations.begin(), partitionsSortedByParticipations.end(), compareParticipations);

            for (size_t i = 0; i < partitionsSortedByParticipations.size(); ++i) {
                if (partScheduled[partitionsSortedByParticipations[i]]) continue;
                WorkElement element = addBlockPair(partitionsSortedByParticipations[i]);
                if (element.part0 != 0 || element.part1 != 0) {
                    tasks.push_back(element);
                    i--;
                }
                if (tasks.size() == partition.numParts() / 2) { break; }
            }

            timer.stop("initialBlockPairs");
            return tasks;
        }

        std::mutex add_lock;
        WorkElement addBlockPair(PartitionID partID) {
            std::lock_guard<std::mutex> lock_guard(add_lock);
            if (partScheduled[partID]) return {0, 0};
            for (PartitionID i = 0; i < partitionsSortedByParticipations.size(); ++i) {
                if (partID != i && isEligible(partID, i)) {
                    if (!partScheduled[i] && blockPairStatus(partID, i) == TaskStatus::UNSCHEDULED) {
                        partScheduled[i] = true;
                        partScheduled[partID] = true;
                        blockPairStatus(partID, i) = TaskStatus::SCHEDULED;
                        participations[partID]--;
                        fixOrdering(partID);
                        participations[i]--;
                        fixOrdering(i);
                        return {partID, i};
                    }
                }
            }
            return {0, 0};
        }

        void fixOrdering(size_t id) {
            for (size_t i = 0; i < partitionsSortedByParticipations.size(); ++i) {
                if (partitionsSortedByParticipations[i] == id) {
                    size_t swapPosition = i + 1;
                    while (swapPosition < partitionsSortedByParticipations.size() && participations[partitionsSortedByParticipations[swapPosition]] > participations[partitionsSortedByParticipations[i]]) { swapPosition++; }
                    swapPosition--;
                    assert(swapPosition < partitionsSortedByParticipations.size());
                    std::swap(partitionsSortedByParticipations[i], partitionsSortedByParticipations[swapPosition]);
                    break;
                }
            }
            if (participations[partitionsSortedByParticipations[partitionsSortedByParticipations.size() - 1]] == 0) { partitionsSortedByParticipations.pop_back(); }
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
            assert(part0 < partition.numParts() && part1 < partition.numParts());
            if (part0 < part1) std::swap(part0, part1);
            return guessNumCutEdges(part0, part1) > 0
                   && (partActive[part0] || partActive[part1])
                   && blockPairStatus(part0, part1) == TaskStatus::UNSCHEDULED
                   && (round < 2 || improvement_history[(part0 * (part0 - 1)/2) + part1] > 0);
        }

        void reportImprovement(PartitionID part0, PartitionID part1) {
            if (part0 < part1) std::swap(part0, part1);
            improvement_history[(part0 * (part0 - 1)/2) + part1]++;
        }



    };
}