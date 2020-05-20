#pragma once

#include <vector>
#include "../util/range.h"
#include <boost/range/irange.hpp>

template<typename T, typename queue_size_type=uint32_t>
class LayeredQueue {
private:
	std::vector<T> queue;
	std::vector<size_t> layer_bounds;
public:
	using size_type = queue_size_type;
	size_type qfront;
	explicit LayeredQueue(const size_type num_elements) : queue(), qfront(0) { layer_bounds.push_back(0); }
	explicit LayeredQueue(const size_t num_elements) : LayeredQueue(static_cast<size_type>(num_elements)) { }
	//Note. Use reinitialize() if you want to keep entries in the underlying vector intact, and ensure these won't be pushed again
	inline void reinitialize(size_type x) { qfront = x; queue.clear(); layer_bounds.clear(); layer_bounds.push_back(x); }
	inline void reinitialize() { reinitialize(queueEnd()); }
	inline void clear() { reinitialize(0); }
	inline bool empty() const { return qfront == queue.size(); }
	inline bool currentLayerEmpty() const { return qfront == layer_bounds.back(); }
	inline T pop() { return queue[qfront++]; }
	//inline T previousLayerPop() { return queue[layerfront++]; }
	inline void finishNextLayer() { layer_bounds.push_back(queue.size()); }
	inline void clearCurrentLayer() { qfront = layer_bounds.back(); }
	inline void push(const T x) { queue.push_back(x); }
	//inline bool previousLayerEmpty() const { return layerfront == layerend; }
	inline T capacity() const { return static_cast<T>(queue.size()); }
	inline std::vector<T>& data() { return queue; }
	template<typename Func> inline void forAllEverContainedElements(Func f) { for (size_type i = 0; i < queue.size(); i++) { f(queue[i]); } }
	inline const_range<std::vector<T>> range(size_type __begin, size_type __end) { return { queue.begin() + __begin, queue.cbegin() + __end }; }
	inline const_range<std::vector<T>> currentLayer() { return range(qfront, layer_bounds.back()); }
	inline const_range<std::vector<T>> allElements() { return range(0, queue.size()); }
	inline size_type queueEnd() const { return queue.size(); }
	inline T popBack() { return queue.pop_back(); }
	inline size_type currentLayerStart() { return qfront; }
	inline size_type currentLayerEnd() { return layer_bounds.back(); }

	inline decltype(auto) currentLayerIndices() { return boost::irange<size_type>(qfront, layer_bounds.back()); }
	/*
	inline void truncateQueue(size_type new_end) {
		qend = std::min(new_end, qend);
		layerend = std::min(new_end, layerend);
		qfront = std::min(new_end, qfront);
	}*/

	inline T elementAt(const size_type pos) const { return queue[pos]; }
	inline void setTo(const size_type pos, T element) { queue[pos] = element; }

	inline T swapFrontToPositionAndPop(size_type pos) {
		std::swap(queue[pos], queue[qfront]);
		return pop();
	}

	std::vector<T> extract() { return std::move(queue); }

	template<typename URBG>
	void shuffleQueue(URBG &&urbg, size_type a, size_type b) {
		std::shuffle(queue.begin() + a, queue.begin() + b, urbg);
	}

	template<typename URBG>
	void shuffleQueue(URBG&& urbg) {
		shuffleQueue(urbg, qfront, queue.size());
	}

	template<typename URBG>
	void shuffleCurrentLayer(URBG &&urbg) {
		shuffleQueue(urbg, qfront, layer_bounds.back());
	}


    /*
	template<bool resize = true>
	void inject(std::vector<T> external_queue, size_type num_elements) {
		queue = std::move(external_queue);
		layerfront = 0;
		layerend = 0;
		qfront = 0;
		qend = queue.size();
		if (resize)
			queue.resize(num_elements);
	}*/
};
