// tensor_shuffle.h — Refined shuffle product, weighted add, and scalar multiply
//
// Adapted from the bootstrap_E6_archive version, with the following changes:
//   - Removed verbose progress monitoring (per-minute stdout, per-thread progress)
//   - Removed tensor_filter_parallel and tensor_remove_parallel (unused, hardcoded threshold)
//   - Simplified the parallel merge (direct concatenation + sort instead of priority queue)
//   - Added tensor_scalar_mul (multiply all entries by a rational scalar)
//   - Cleaned up redundant safety checks
//
// Core functions:
//   tensor_shuffle_product_parallel(A, B, F, pool) — shuffle product of two COO tensors
//   tensor_add_weighted(A, B, w_A, w_B, F)        — A = w_A*A + w_B*B (in-place)
//   tensor_scalar_mul(A, scalar, F)               — returns scalar * A (new tensor)

#ifndef TENSOR_SHUFFLE_H
#define TENSOR_SHUFFLE_H

#include <iostream>
#include <thread>
#include <chrono>
#include <unordered_map>
#include <algorithm>
#include <vector>
#include "SparseRREF/sparse_tensor.h"

namespace SparseRREF {

// Hash function for std::vector<int> (used by unordered_map in shuffle product)
struct VectorHash {
	std::size_t operator()(const std::vector<int>& vec) const {
		std::size_t seed = vec.size();
		for (const auto& i : vec) {
			seed ^= std::hash<int>{}(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		}
		return seed;
	}
};

// Generate all shuffle interleavings of two index sequences using Gosper's hack.
// For each interleaving, accumulate the product value into position_sums.
template <typename index_type, typename T>
void process_shuffles_combinatorial(
	const std::vector<index_type>& seqA,
	const std::vector<index_type>& seqB,
	const T& product_val,
	std::unordered_map<std::vector<index_type>, T, VectorHash>& position_sums) {

	if (seqA.empty() && seqB.empty()) return;
	if (seqA.empty()) { position_sums[seqB] += product_val; return; }
	if (seqB.empty()) { position_sums[seqA] += product_val; return; }

	size_t total_len = seqA.size() + seqB.size();
	size_t a_count = seqA.size();

	// Safety: total_len > 11 would overflow 64-bit in the worst case (2^total_len).
	// In practice, weight <= 11 covers all realistic loop orders.
	if (total_len > 11) {
		throw std::runtime_error("tensor_shuffle_product: total weight " + std::to_string(total_len) + " > 11");
	}

	// Use Gosper's hack to iterate over all combinations of a_count positions
	std::vector<bool> mask(total_len, false);
	std::fill(mask.begin(), mask.begin() + a_count, true);

	do {
		std::vector<index_type> current;
		current.reserve(total_len);

		size_t a_idx = 0, b_idx = 0;
		for (bool use_a : mask) {
			if (use_a) {
				current.push_back(seqA[a_idx++]);
			} else {
				current.push_back(seqB[b_idx++]);
			}
		}
		position_sums[current] += product_val;
	} while (std::prev_permutation(mask.begin(), mask.end()));
}

// Parallel shuffle product: generates all interleavings of positions from two
// sparse tensors. The result has rank = A.rank() + B.rank(), with all dims equal
// to the maximum dimension across both tensors (the "letter" space).
template <typename index_type, typename T>
sparse_tensor<T, index_type, SPARSE_COO> tensor_shuffle_product_parallel(
	const sparse_tensor<T, index_type, SPARSE_COO>& A,
	const sparse_tensor<T, index_type, SPARSE_COO>& B,
	const field_t& F, thread_pool* pool = nullptr) {

	if (A.nnz() == 0 || B.nnz() == 0) {
		throw std::invalid_argument("tensor_shuffle_product_parallel: one or both tensors are empty.");
	}

	// Result dims: rank = A.rank() + B.rank(), all dims = max(dim across A and B).
	// NOTE: B.dims() returns by value (a temporary std::vector<size_t>); calling
	// B.dims().begin() and B.dims().end() separately creates iterators from TWO
	// different temporaries — a dangling-iterator bug that triggers
	// "vector::_M_range_insert" assertions. Capture into a local first.
	std::vector<size_t> dimsCp = A.dims();
	auto bdims = B.dims();
	dimsCp.insert(dimsCp.end(), bdims.begin(), bdims.end());
	size_t maxVal = *max_element(dimsCp.begin(), dimsCp.end());
	std::vector<size_t> dimsC(A.rank() + B.rank(), maxVal);

	sparse_tensor<T, index_type, SPARSE_COO> C(dimsC);

	// --- Sequential version ---
	// Strategy: accumulate ALL shuffle results into a single unordered_map
	// (key=multi-index, value=accumulated sum), then build the result tensor
	// with push_back. This matches the parallel version's approach and avoids
	// the O(n²) ordered insert_add path which was observed to corrupt tensor
	// dims for large results (e.g. E1^3 with 1397 entries).
	if (pool == nullptr) {
		auto permA = A.gen_perm();
		auto permB = B.gen_perm();

		std::unordered_map<std::vector<index_type>, T, VectorHash> accumulated;
		std::unordered_map<std::vector<index_type>, T, VectorHash> position_sums;

		for (auto i : permA) {
			std::vector<index_type> seqA = A.index_vector(i);
			for (auto j : permB) {
				std::vector<index_type> seqB = B.index_vector(j);
				T product_val = scalar_mul(A.val(i), B.val(j), F);

				position_sums.clear();
				process_shuffles_combinatorial(seqA, seqB, product_val, position_sums);

				for (const auto& pair : position_sums) {
					accumulated[pair.first] += pair.second;
				}
			}
		}

		// Build the result tensor from the accumulated map (skip zero entries).
		C.reserve(accumulated.size() + 1);
		for (const auto& pair : accumulated) {
			if (pair.second != (T)0) {
				C.push_back(pair.first, pair.second);
			}
		}
		C.canonicalize();  // remove any zeros created by cancellation
		C.sort_indices();
		return C;
	}

	// --- Parallel version ---
	//
	// Strategy: each thread accumulates results into its own unordered_map
	// (key=multi-index, value=accumulated sum). After all threads finish, the
	// maps are merged into one, then the result tensor is built.
	//
	// This avoids the non-existent sparse_tensor::merge() method and correctly
	// handles duplicate entries (different (i,j) pairs can produce the same
	// shuffle position, and entries from different threads can collide).
	size_t nthread = pool->get_thread_count();
	std::vector<std::unordered_map<std::vector<index_type>, T, VectorHash>> thread_maps(nthread);

	// Permutations are shared across threads (read-only after generation).
	auto permA = A.gen_perm();
	auto permB = B.gen_perm();

	size_t total_work = A.nnz() * B.nnz();

	pool->detach_blocks(0, total_work, [&](size_t start, size_t end) {
		size_t tid = thread_id();
		auto& local_map = thread_maps[tid];

		// Reusable per-(i,j) map (cleared each iteration).
		std::unordered_map<std::vector<index_type>, T, VectorHash> position_sums;

		for (size_t work_idx = start; work_idx < end; work_idx++) {
			size_t i = work_idx / B.nnz();
			size_t j = work_idx % B.nnz();

			auto posA = permA[i];
			auto posB = permB[j];

			std::vector<index_type> seqA = A.index_vector(posA);
			std::vector<index_type> seqB = B.index_vector(posB);

			T product_val = scalar_mul(A.val(posA), B.val(posB), F);

			position_sums.clear();
			process_shuffles_combinatorial(seqA, seqB, product_val, position_sums);

			for (const auto& pair : position_sums) {
				local_map[pair.first] += pair.second;
			}
		}
	}, nthread);

	pool->wait();

	// Merge all thread maps into one, skipping zero entries.
	std::unordered_map<std::vector<index_type>, T, VectorHash> merged;
	for (auto& tm : thread_maps) {
		for (const auto& pair : tm) {
			if (pair.second != (T)0) {
				merged[pair.first] += pair.second;
			}
		}
	}
	thread_maps.clear();

	// Build the result tensor from the merged map.
	C.reserve(merged.size());
	for (const auto& pair : merged) {
		if (pair.second != (T)0) {
			C.push_back(pair.first, pair.second);
		}
	}
	C.canonicalize();  // remove any zeros created by cancellation
	return C;
}

// Weighted tensor addition: A = weight_A * A + weight_B * B (modifies A in-place).
// Handles empty tensors (alloc == 0 or nnz == 0) gracefully.
template <typename index_type, typename T>
void tensor_add_weighted(
	sparse_tensor<T, index_type, SPARSE_COO>& A,
	const sparse_tensor<T, index_type, SPARSE_COO>& B,
	const T& weight_A,
	const T& weight_B,
	const field_t& F) {

	// Both empty: nothing to do
	if (A.alloc() == 0 && B.alloc() == 0) return;

	// A is empty: A = weight_B * B
	if (A.alloc() == 0) {
		A = sparse_tensor<T, index_type, SPARSE_COO>(B.dims());
		A.reserve(B.nnz());
		auto Bperm = B.gen_perm();
		for (auto i : Bperm) {
			auto val = scalar_mul(B.val(i), weight_B, F);
			if (val != 0) A.push_back(B.index_vector(i), val);
		}
		return;
	}

	// B is empty: A = weight_A * A
	if (B.alloc() == 0) {
		sparse_tensor<T, index_type, SPARSE_COO> result(A.dims());
		result.reserve(A.nnz());
		auto Aperm = A.gen_perm();
		for (auto i : Aperm) {
			auto val = scalar_mul(A.val(i), weight_A, F);
			if (val != 0) result.push_back(A.index_vector(i), val);
		}
		A = std::move(result);
		return;
	}

	// Dimension check
	if (A.rank() != B.rank()) {
		throw std::runtime_error("tensor_add_weighted: rank mismatch");
	}
	for (size_t i = 0; i < A.rank(); i++) {
		if (A.dim(i) != B.dim(i)) {
			throw std::runtime_error("tensor_add_weighted: dimension mismatch at axis " + std::to_string(i));
		}
	}

	auto rank = A.rank();

	// Handle zero nnz cases
	if (A.nnz() == 0) {
		sparse_tensor<T, index_type, SPARSE_COO> result(B.dims());
		result.reserve(B.nnz());
		auto Bperm = B.gen_perm();
		for (auto i : Bperm) {
			auto val = scalar_mul(B.val(i), weight_B, F);
			if (val != 0) result.push_back(B.index_vector(i), val);
		}
		A = std::move(result);
		return;
	}
	if (B.nnz() == 0) {
		sparse_tensor<T, index_type, SPARSE_COO> result(A.dims());
		result.reserve(A.nnz());
		auto Aperm = A.gen_perm();
		for (auto i : Aperm) {
			auto val = scalar_mul(A.val(i), weight_A, F);
			if (val != 0) result.push_back(A.index_vector(i), val);
		}
		A = std::move(result);
		return;
	}

	// Merge two sorted tensors with weighted values
	sparse_tensor<T, index_type, SPARSE_COO> C(A.dims());
	C.reserve(A.nnz() + B.nnz());

	auto Aperm = A.gen_perm();
	auto Bperm = B.gen_perm();

	size_t i = 0, j = 0;
	while (i < A.nnz() && j < B.nnz()) {
		auto posA = Aperm[i];
		auto posB = Bperm[j];
		int cmp = lexico_compare(A.index(posA), B.index(posB), rank);

		if (cmp < 0) {
			auto val = scalar_mul(A.val(posA), weight_A, F);
			if (val != 0) C.push_back(A.index_vector(posA), val);
			i++;
		} else if (cmp > 0) {
			auto val = scalar_mul(B.val(posB), weight_B, F);
			if (val != 0) C.push_back(B.index_vector(posB), val);
			j++;
		} else {
			auto valA = scalar_mul(A.val(posA), weight_A, F);
			auto valB = scalar_mul(B.val(posB), weight_B, F);
			auto val = scalar_add(valA, valB, F);
			if (val != 0) C.push_back(A.index_vector(posA), val);
			i++; j++;
		}
	}
	while (i < A.nnz()) {
		auto posA = Aperm[i];
		auto val = scalar_mul(A.val(posA), weight_A, F);
		if (val != 0) C.push_back(A.index_vector(posA), val);
		i++;
	}
	while (j < B.nnz()) {
		auto posB = Bperm[j];
		auto val = scalar_mul(B.val(posB), weight_B, F);
		if (val != 0) C.push_back(B.index_vector(posB), val);
		j++;
	}

	A = std::move(C);
}

// Scalar multiply: returns a new tensor with all entries multiplied by scalar.
// If scalar is 0, returns an empty tensor with the same dims.
template <typename index_type, typename T>
sparse_tensor<T, index_type, SPARSE_COO> tensor_scalar_mul(
	const sparse_tensor<T, index_type, SPARSE_COO>& A,
	const T& scalar,
	const field_t& F) {

	sparse_tensor<T, index_type, SPARSE_COO> C(A.dims());

	if (A.nnz() == 0) return C;

	C.reserve(A.nnz());
	auto perm = A.gen_perm();
	for (auto i : perm) {
		auto val = scalar_mul(A.val(i), scalar, F);
		if (val != 0) {
			C.push_back(A.index_vector(i), val);
		}
	}
	C.canonicalize();
	return C;
}

} // namespace SparseRREF

#endif // TENSOR_SHUFFLE_H
