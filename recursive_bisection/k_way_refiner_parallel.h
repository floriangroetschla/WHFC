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
                timer(timer), mt(mt), config(config),
                refiners_thread_specific(partition.getGraph().numNodes(), partition.getGraph().numHyperedges(), partition.getGraph().numPins(), mt(), std::ref(config)),
                timers_thread_specific(), bucketPQ(partition.numParts())
                {

                }

        uint refine(double epsilon, uint maxIterations) {
            NodeWeight maxWeight = (1.0 + epsilon) * partition.totalWeight() / static_cast<double>(partition.numParts());
            partActive.assign(partActive.size(), 1);
            partActiveNextRound.assign(partActive.size(), 0);
            iterationCounter = 0;

            // TODO maybe also make a row restriction. or terminate if there aren't enough initial tasks. 
            while (std::any_of(partActive.begin(), partActive.end(), [](auto& x) { return x > 0; }) && iterationCounter < maxIterations) {
                std::vector<WorkElement> tasks = initialBlockPairs();
                //std::cout << "Round " << round << " initial tasks: " << tasks.size() << std::endl;
                tbb::parallel_do(tasks,
                        [&](WorkElement element, tbb::parallel_do_feeder<WorkElement>& feeder) {
                            assert(partScheduled[element.part0]);
                            assert(partScheduled[element.part1]);

                            if (iterationCounter < maxIterations) {
                                WHFCRefinerTwoWay& refiner = refiners_thread_specific.local();
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

                                //partScheduled[element.part0] = false;
                                //partScheduled[element.part1] = false;
                                timer.start("addNewTasks");
                                //addNewTasksFromPQ(feeder, maxWeight);

                                if (!addNewTasks(element.part0, feeder, maxWeight)) {
                                    partScheduled[element.part0] = false;
                                }
                                if (!addNewTasks(element.part1, feeder, maxWeight)) {
                                    partScheduled[element.part1] = false;
                                }

                                timer.stop("addNewTasks");
                                iterationCounter++;
                            }
                        }
                );
                //std::cout << "WHFC refiner calls: " << iterationCounter << std::endl;

                //assert(allPairsProcessed()); only if maxIterations allows it
                std::swap(partActive, partActiveNextRound);
                partActiveNextRound.assign(partActive.size(), 0);
                round++;
            }

            for (whfc::TimeReporter& local_timer : timers_thread_specific) {
                timer.merge(local_timer, "Refinement", "total");
            }

            return iterationCounter.load();
        }

        template<typename Element>
        class BucketPriorityQueue {
        public:
            BucketPriorityQueue() = delete;

            BucketPriorityQueue(size_t numBuckets) : buckets(numBuckets),
                first_non_empty_bucket(std::numeric_limits<size_t>::max()), container_size(0) {

            }

            void clear() {
                for (std::vector<Element> bucket : buckets) {
                    bucket.clear();
                }
                first_non_empty_bucket = std::numeric_limits<size_t>::max();
                container_size = 0;
            }

            size_t size() {
                return container_size;
            }

            size_t numBuckets() {
                return buckets.size();
            }

            void addElement(size_t key, Element element) {
                assert(key < buckets.size());
                buckets[key].push_back(element);
                if (key < first_non_empty_bucket) first_non_empty_bucket = key;
                container_size++;
            }

            std::vector<Element> getBucket(size_t i) {
                assert(i < buckets.size());
                return buckets[i];
            }

            size_t firstNonEmptyBucket() {
                return first_non_empty_bucket;
            }

            Element getNextElement() {
                assert(container_size > 0);
                assert(first_non_empty_bucket != std::numeric_limits<size_t>::max());
                container_size--;
                const Element element = buckets[first_non_empty_bucket].pop_back();
                if (container_size != 0) {
                    while (buckets[first_non_empty_bucket].size() == 0) {
                        first_non_empty_bucket++;
                    }
                } else {
                    first_non_empty_bucket = std::numeric_limits<size_t>::max();
                }
                return element;
            }

            void removeElement(size_t bucket, size_t i) {
                assert(bucket < buckets.size() && i < buckets[bucket].size());
                buckets[bucket][i] = buckets[bucket].back();
                buckets[bucket].pop_back();
                container_size--;
            }

        private:
            std::vector<std::vector<Element>> buckets;
            size_t first_non_empty_bucket;
            size_t container_size;
        };

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
        const PartitionConfig& config;
        tbb::enumerable_thread_specific<WHFCRefinerTwoWay> refiners_thread_specific;
        tbb::enumerable_thread_specific<whfc::TimeReporter> timers_thread_specific;
        BucketPriorityQueue<WorkElement> bucketPQ;


        size_t guessNumCutEdges(PartitionID part0, PartitionID part1) {
            return partition.cutEdges(part0, part1).size() + partition.cutEdges(part1, part0).size();
        }
        
        std::vector<WorkElement> initialBlockPairs() {
            timer.start("initialBlockPairs", "Refinement");
            resetBlockPairStatus();
            bucketPQ.clear();

            std::vector<size_t> participations(partition.numParts(), 0);

            for (PartitionID i = 0; i < partition.numParts() - 1; ++i) {
                for (PartitionID j = i + 1; j < partition.numParts(); ++j) {
                    if (isEligible(i, j)) {
                        participations[i]++;
                        participations[j]++;
                    }
                }
            }

            for (PartitionID i = 0; i < partition.numParts() - 1; ++i) {
                for (PartitionID j = i + 1; j < partition.numParts(); ++j) {
                    if (isEligible(i, j)) {
                        bucketPQ.addElement(std::max(participations[i], participations[j]), {i, j});
                    }
                }
            }

            std::vector<WorkElement> tasks;
            size_t bucket = bucketPQ.firstNonEmptyBucket();
            while (bucket < bucketPQ.numBuckets()) {
                std::vector<WorkElement> bucket_elements = bucketPQ.getBucket(bucket);
                size_t current_size = bucket_elements.size();
                for (size_t i = 0; i < current_size; ++i) {
                    const WorkElement element = bucket_elements[i];
                    if (!partScheduled[element.part0] && !partScheduled[element.part1]) {
                        tasks.push_back(element);
                        blockPairStatus(element.part0, element.part1) = TaskStatus::SCHEDULED;
                        partScheduled[element.part0] = true;
                        partScheduled[element.part1] = true;
                        bucketPQ.removeElement(bucket, i);
                        i--;
                        current_size--;
                    }
                }
                bucket++;
            }

            timer.stop("initialBlockPairs");
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

        std::mutex add_lock;
        void addNewTasksFromPQ(tbb::parallel_do_feeder<WorkElement>& feeder, NodeWeight maxWeight) {
            std::lock_guard<std::mutex> lock_guard(add_lock);
            if (bucketPQ.size() > 0) {
                size_t bucket = bucketPQ.firstNonEmptyBucket();
                while (bucket < bucketPQ.numBuckets()) {
                    std::vector<WorkElement> bucket_elements = bucketPQ.getBucket(bucket);
                    size_t current_size = bucket_elements.size();
                    for (size_t i = 0; i < current_size; ++i) {
                        const WorkElement element = bucket_elements[i];
                        if (!partScheduled[element.part0] && !partScheduled[element.part1] && isEligible(element.part0, element.part1)) {
                            feeder.add({element});
                            blockPairStatus(element.part0, element.part1) = TaskStatus::SCHEDULED;
                            partScheduled[element.part0] = true;
                            partScheduled[element.part1] = true;
                            bucketPQ.removeElement(bucket, i);
                            i--;
                            current_size--;
                        }
                    }
                    bucket++;
                }
            }
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