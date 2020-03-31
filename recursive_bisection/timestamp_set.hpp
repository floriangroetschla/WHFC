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

}
