#pragma once

#include <boost/dynamic_bitset.hpp>
#include <atomic>
#include <tbb/parallel_do.h>
#include <mutex>

namespace whfc_rb {
    class KWayRefinerParallel {
    public:
        explicit KWayRefinerParallel(PartitionCA &partition, whfc::TimeReporter &timer, std::mt19937 &mt) : partition(
                partition), partActive(partition.numParts()), partActiveNextRound(partition.numParts()), blockPairStatus(partition.numParts() * (partition.numParts() - 1) / 2), partitionScheduled(partition.numParts()), timer(timer), mt(mt) {}

        void refine(double epsilon, uint maxIterations) {
            NodeWeight maxWeight = (1.0 + epsilon) * partition.totalWeight() / static_cast<double>(partition.numParts());
            partActive.set();
            partActiveNextRound.reset();
            iterationCounter = 0;

            while (partActive.count() > 0 && iterationCounter < maxIterations) {
                std::vector<WorkElement> tasks;
                PartitionBase::PartitionID part0 = 0;
                PartitionBase::PartitionID part1 = partition.numParts() - 1;
                resetBlockPairStatus();
                while (part0 < part1) {
                    if (partActive[part0] || partActive[part1]) {
                        WorkElement element = {*this, part0, part1, maxWeight};
                        tasks.push_back(element);
                        blockPair(part0, part1) = TaskStatus::SCHEDULED;
                        partitionScheduled[part0] = true;
                        partitionScheduled[part1] = true;
                    }
                    part0++;
                    part1--;
                }

                tbb::parallel_do(tasks,
                        [](WorkElement element, tbb::parallel_do_feeder<WorkElement>& feeder)
                        {
                            assert(element.refiner.partitionScheduled[element.part0]);
                            assert(element.refiner.partitionScheduled[element.part1]);

                            std::mt19937 mt(element.refiner.mt());
                            whfc::TimeReporter timer;
                            WHFCRefinerTwoWay refiner(element.refiner.partition.getGraph().numNodes(), element.refiner.partition.getGraph().numHyperedges(), element.refiner.partition.getGraph().numPins(), mt, timer);
                            bool refinementResult = refiner.refine(element.refiner.partition, element.part0, element.part1, element.maxWeight, element.maxWeight);
                            if (refinementResult) {
                                // Schedule for next round
                                element.refiner.partActiveNextRound.set(element.part0);
                                element.refiner.partActiveNextRound.set(element.part1);
                            }

                            element.refiner.blockPair(element.part0, element.part1) = TaskStatus::FINISHED;

                            if (!element.refiner.addNewTasks(element.part0, feeder, element.maxWeight)) {
                                element.refiner.partitionScheduled[element.part0] = false;
                            }
                            if (!element.refiner.addNewTasks(element.part1, feeder, element.maxWeight)) {
                                element.refiner.partitionScheduled[element.part1] = false;
                            }
                            element.refiner.iterationCounter++;
                        });

                std::swap(partActive, partActiveNextRound);
                partActiveNextRound.reset();
            }
        }

    private:
        enum class TaskStatus {UNSCHEDULED, SCHEDULED, FINISHED};

        PartitionCA &partition;
        boost::dynamic_bitset<> partActive;
        boost::dynamic_bitset<> partActiveNextRound;
        std::vector<std::atomic<TaskStatus>> blockPairStatus;
        std::vector<std::atomic<bool>> partitionScheduled;
        whfc::TimeReporter &timer;
        std::mt19937 &mt;
        std::atomic<uint> iterationCounter = 0;

        struct WorkElement {
            KWayRefinerParallel& refiner;
            PartitionBase::PartitionID part0;
            PartitionBase::PartitionID part1;
            NodeWeight maxWeight;
        };

        bool addNewTasks(PartitionBase::PartitionID part, tbb::parallel_do_feeder<WorkElement>& feeder, NodeWeight maxWeight) {
            bool foundPair = false;
            for (PartitionBase::PartitionID pid = 0; pid < partition.numParts(); ++pid) {
                if (pid != part && !partitionScheduled[pid]) {
                    if (!partitionScheduled[pid].exchange(true)) {
                        if (blockPair(part, pid) == TaskStatus::UNSCHEDULED) {
                            blockPair(part, pid) = TaskStatus::SCHEDULED;
                            feeder.add({*this, part, pid, maxWeight});
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

        void resetBlockPairStatus() {
            for (uint i = 0; i < blockPairStatus.size(); ++i) {
                blockPairStatus[i] = TaskStatus::UNSCHEDULED;
            }
        }

        std::atomic<TaskStatus>& blockPair(PartitionBase::PartitionID part0, PartitionBase::PartitionID part1) {
            if (part0 < part1) std::swap(part0, part1);
            return (blockPairStatus[(part0 * (part0 - 1) / 2) + part1]);
        }

    };
}