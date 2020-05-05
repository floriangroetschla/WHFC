#pragma once

#include <vector>
#include <cstdint>

namespace ldc {

// TODO: consider making the underlying containers generic to e.g. support different vector allocators?

template<typename TimestampT = uint16_t>
class TimestampSet {
public:
	TimestampSet(size_t n) : timestamps(n, 0), generation(1) {

	}

	bool contains(size_t i) const {
		return timestamps[i] == generation; 
	}

	void add(size_t i) {
		timestamps[i] = generation;
	}

	void remove(size_t i) {
		timestamps[i] = 0;
	}

	void clear() {
		if (generation == std::numeric_limits<TimestampT>::max()) {
			timestamps.assign(timestamps.size(), 0);
			generation = 0;
		}
		generation++;
	}

private:
	std::vector<TimestampT> timestamps;
	TimestampT generation;
};

template<typename T, typename TimestampT = uint16_t>
class TimestampMap {
public:
	TimestampMap(size_t n, T defaultValue = T()) : defaultValue(defaultValue), timestamps(n, 0), generation(0), map(n, defaultValue) {

	}

	T& operator[](size_t i) {
		if (timestamps[i] != generation) {
			timestamps[i] = generation;
			map[i] = defaultValue;
		}
		return map[i];
	}

	void clear() {
	if (generation == std::numeric_limits<TimestampT>::max()) {
			timestamps.assign(timestamps.size(), 0);
			generation = 0;
		}
		generation++;	
	}

	const T defaultValue;
private:
	std::vector<TimestampT> timestamps;
	TimestampT generation;
	std::vector<T> map;
};

template<typename TimestampT = uint16_t>
class AtomicTimestampSet {
public:
    AtomicTimestampSet(size_t n) : timestamps(n), generation(0) {}

    bool set(size_t i) {
        return timestamps[i].exchange(generation, std::memory_order_relaxed) != generation;	// REVIEW NOTE. specify memory order? efficiency of atomic operations on 16 bit ints?
    }

    bool isSet(size_t i) const {
        return timestamps[i].load(std::memory_order_relaxed) == generation;
    }

    void reset(size_t i) {
        timestamps[i] = generation - 1;
    }

    void reset() {
        if (generation == std::numeric_limits<TimestampT>::max()) {
            timestamps = std::vector<std::atomic<TimestampT>>(timestamps.size());
            generation = 0;
        }
        generation++;
    }

private:
    std::vector<std::atomic<TimestampT>> timestamps;
    std::atomic<TimestampT> generation;
};

template<typename TimestampT = uint32_t>
class AtomicTimestampCounter {
public:
    AtomicTimestampCounter(size_t n) : timestamps(n), generation_start(0), generation_current(0) {}

    bool isSet(size_t i) const {
        return timestamps[i] >= generation_start && timestamps[i] <= generation_current;
    }

    TimestampT set(size_t i) {
        return std::max<TimestampT>(timestamps[i].exchange(generation_current) - generation_start, 0);
    }

    TimestampT value(size_t i) {
        assert(isSet(i));
        return timestamps[i] - generation_start;
    }

    TimestampT getCounter() {
        return generation_current - generation_start;
    }

    void nextRound() {
        assert(generation_current != std::numeric_limits<TimestampT>::max());
        generation_current++;
    }

    void reset() {
        generation_start = ++generation_current;
    }

private:
    std::vector<std::atomic<TimestampT>> timestamps;
    std::atomic<TimestampT> generation_start;
    std::atomic<TimestampT> generation_current;
};


}
