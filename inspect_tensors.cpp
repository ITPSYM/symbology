// inspect_tensors.cpp — Inspect E1, compute E1^2/2, verify boundary
//
// Usage: ./inspect_tensors [--data-dir <dir>] [--output-dir <dir>]
// Reads <output-dir>/oneloop/E1.wxf and <output-dir>/2loop/boundary_2L.wxf.
// Defaults: --output-dir resolves to <exec_dir>/output (relative paths are
// resolved against the executable directory).
#include "bootstrap.hpp"
#include "projection.hpp"
#include "tensor_shuffle.h"
#include <iostream>
#include <set>
#include <map>

using index_t = int32_t;
using scalar_t = rat_t;

int main(int argc, char* argv[]) {
	std::filesystem::path output_dir;
	std::filesystem::path base = std::filesystem::path(argv[0]).parent_path();
	if (base.empty()) base = std::filesystem::current_path();

	for (int i = 1; i < argc; i++) {
		std::string arg = argv[i];
		if (arg == "--output-dir" && i + 1 < argc) {
			output_dir = argv[++i];
		}
		else if (arg == "--data-dir" && i + 1 < argc) {
			++i;  // accepted for symmetry with the other tools; not used here
		}
		else if (arg == "-h" || arg == "--help") {
			std::cerr << "Usage: " << argv[0] << " [--data-dir <dir>] [--output-dir <dir>]\n";
			return 0;
		}
	}
	if (output_dir.empty()) output_dir = base / "output";
	if (!output_dir.is_absolute()) output_dir = base / output_dir;

	field_t F(FIELD_QQ);
	rref_option_t opt;
	opt->method = 0;
	opt->verbose = false;
	opt->pool.reset();
	thread_pool* pool = &(opt->pool);

	// Read E1
	auto e1_path = output_dir / "oneloop" / "E1.wxf";
	auto E1_csr = projection_read_tensor<scalar_t, index_t>(e1_path, F, pool);
	sparse_tensor<scalar_t, index_t, SPARSE_COO> E1(std::move(E1_csr));
	std::cout << "=== E1 (rank=" << E1.rank() << ", nnz=" << E1.nnz() << ") ===" << std::endl;
	for (size_t i = 0; i < E1.nnz(); i++) {
		auto coords = E1.index(i);
		std::cout << "   E1[" << coords[0] << "," << coords[1] << "] = " << E1.val(i) << std::endl;
	}

	// Compute E1^2 using SEQUENTIAL shuffle product
	constexpr thread_pool* seq_pool = nullptr;
	auto E1_sq = tensor_shuffle_product_parallel(E1, E1, F, seq_pool);
	std::cout << "\n=== E1^2 (shuffle, nnz=" << E1_sq.nnz() << ") ===" << std::endl;

	// Compute boundary = E1^2 / 2
	auto boundary = tensor_scalar_mul(E1_sq, scalar_t(1, 2), F);
	std::cout << "=== boundary = E1^2/2 (nnz=" << boundary.nnz() << ") ===" << std::endl;

	// Dump first 20 boundary entries
	std::cout << "   First 20 boundary entries:" << std::endl;
	size_t count = 0;
	for (auto i : boundary.gen_perm()) {
		auto coords = boundary.index(i);
		std::cout << "   b[";
		for (size_t d = 0; d < boundary.rank(); d++) {
			std::cout << coords[d];
			if (d + 1 < boundary.rank()) std::cout << ",";
		}
		std::cout << "] = " << boundary.val(i) << std::endl;
		count++;
		if (count >= 20) break;
	}

	// Read the boundary file and compare
	auto bfile_path = output_dir / "2loop" / "boundary_2L.wxf";
	auto bfile_csr = projection_read_tensor<scalar_t, index_t>(bfile_path, F, pool);
	sparse_tensor<scalar_t, index_t, SPARSE_COO> bfile(std::move(bfile_csr));
	std::cout << "\n=== boundary file (nnz=" << bfile.nnz() << ") ===" << std::endl;

	// Compare: build map from computed boundary
	std::map<std::vector<index_t>, scalar_t> computed_map;
	for (auto i : boundary.gen_perm()) {
		auto coords = boundary.index(i);
		std::vector<index_t> key(coords, coords + boundary.rank());
		computed_map[key] = boundary.val(i);
	}

	size_t match = 0, mismatch = 0, file_only = 0, computed_only = 0;
	for (auto i : bfile.gen_perm()) {
		auto coords = bfile.index(i);
		std::vector<index_t> key(coords, coords + bfile.rank());
		auto it = computed_map.find(key);
		if (it != computed_map.end()) {
			if (it->second == bfile.val(i)) {
				match++;
			} else {
				mismatch++;
				std::cout << "   MISMATCH at [";
				for (size_t d = 0; d < key.size(); d++) {
					std::cout << key[d];
					if (d + 1 < key.size()) std::cout << ",";
				}
				std::cout << "]: file=" << bfile.val(i) << " computed=" << it->second << std::endl;
			}
		} else {
			file_only++;
		}
	}
	for (const auto& [key, val] : computed_map) {
		// Check if key is in file
		bool found = false;
		for (auto i : bfile.gen_perm()) {
			auto coords = bfile.index(i);
			std::vector<index_t> fkey(coords, coords + bfile.rank());
			if (fkey == key) { found = true; break; }
		}
		if (!found) computed_only++;
	}

	std::cout << "\nComparison: match=" << match << " mismatch=" << mismatch
	          << " file_only=" << file_only << " computed_only=" << computed_only << std::endl;

	// Verify specific values by hand
	std::cout << "\n=== Hand verification ===" << std::endl;
	// b[0,0,0,0] should be 12 (6 shuffles of [0,0]x[0,0], each 4, /2 = 12)
	// b[0,0,1,1] should be 4 (2 shuffles, each 4, /2 = 4)
	// b[0,0,2,6] should be 4 (2 shuffles, each 4, /2 = 4)

	return 0;
}
