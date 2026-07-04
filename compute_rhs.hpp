// compute_rhs.hpp — Compute the RHS (collinear boundary) for collinear constraints
//
// This is a standalone module that calculates the right-hand side of the
// collinear constraint equation c.A = b, where:
//   - A is the expanded SEW tensor after collinear projection
//   - b (the RHS / boundary) is the expected collinear limit
//
// The collinear limit is non-zero because of the normalization difference
// between BDS ansatz and BDS-like ansatz. This difference is an exponential
// that can be expanded loop order by loop order (each loop order L has
// weight 2L, so SEW_5p1 = (5+1)/2 = 3 loop order).
//
// The one-loop seed E1 is provided in data/E1.wxf (rank 2, (11,11)).
// The boundary at L loops is computed recursively:
//   L=2: E1^2 / 2
//   L=3: E1^3 / 6 + E1 * R2
//   L=4: -E1^4 / 12 + E2^2 / 2 + E1 * R3        (not yet implemented)
//   L=5: E1^5 / 20 - E1*E2^2 / 2 + E2*R3 + E1*R4  (not yet implemented)
//
// where E_L is the expanded collinear projection of hepMHV_LL, and
// R_L = E_L - boundary_L (the "remainder" at loop L).
//
// hepMHV_LL = contract(solMHV_LL, projected_SEW_basis, axis 0)
// E_L = expand(hepMHV_LL) to rank 2L, dims (11,...,11)
//
// solMHV_LL is obtained by solving c.A = boundary_L for the SEW_{(2L-1)}p1 target.
// This module invokes the solve-collinear logic for lower loop orders recursively.
//
// Naming convention for outputs:
//   output/oneloop/E1.wxf                           (copy of data/E1.wxf)
//   output/twoloop/{solMHV_2L, hepMHV_2L, E2, R2, boundary_2L}.wxf
//   output/threeloop/{solMHV_3L, hepMHV_3L, E3, R3, boundary_3L}.wxf
//
// SEW-level projections use naming: colprojdiv_{sew_name}.wxf (e.g. colprojdiv_SEW_5p1.wxf)
//
// IMPORTANT (Q3): We project SEW_* to its collinear expression by using the SEW
// collinear-reduced basis (SEW_{name}_basis.wxf, produced by `--project`) directly.
// We do NOT apply colprojdiv to the SEW basis. The colprojdiv operation is only
// conceptual: it ensures all parts contain a divergent component, so the divergent
// part of the collinear limit alone can fix all the parameters. In practice, the
// SEW basis already encodes the collinear expression, and colprojdiv at the SEW
// level is identity, so applying it would be a no-op (and is conceptually wrong
// per Q3: "no need to further project to divergent part").

#ifndef COMPUTE_RHS_HPP
#define COMPUTE_RHS_HPP

#include "bootstrap.hpp"
#include "projection.hpp"
#include "solve_collinear.hpp"
#include "tensor_expand.hpp"
#include "linear_solve.hpp"
#include "tensor_shuffle.h"

#include <set>
#include <map>
#include <utility>

// ========== Utility: write a CSR tensor to a WXF file ==========

template <typename T, typename index_t>
void write_tensor_file(
	const std::filesystem::path& path,
	sparse_tensor<T, index_t, SPARSE_CSR>&& tensor_csr) {

	auto u8arr = sparse_tensor_write_wxf(tensor_csr);
	std::ofstream ofs(path, std::ios::binary);
	ofs.write(reinterpret_cast<const char*>(u8arr.data()), u8arr.size());
	ofs.close();
	std::cout << "   Wrote " << path.filename().string() << std::endl;
}

// ========== Shuffle power: compute E^n via repeated shuffle product ==========
//
// E^n = shuffle(shuffle(...shuffle(E, E)..., E), E)  [n times]
// The result has rank = n * E.rank() and all dims = max dim of E.

template <typename T, typename index_t>
sparse_tensor<T, index_t, SPARSE_COO> shuffle_power(
	const sparse_tensor<T, index_t, SPARSE_COO>& E,
	size_t n,
	const field_t& F, thread_pool* /*pool*/) {

	if (n == 0) {
		throw std::runtime_error("shuffle_power: n=0 is not supported");
	}
	if (n == 1) {
		return E;  // copy
	}

	// Use the SEQUENTIAL version of tensor_shuffle_product (pool=nullptr).
	// The parallel version has been observed to produce incorrect results for
	// the boundary calculation (E1^2/2), so we force sequential here.
	constexpr thread_pool* seq_pool = nullptr;

	// Start with E^1 = E, then repeatedly shuffle with E
	auto result = tensor_shuffle_product_parallel(E, E, F, seq_pool);  // E^2
	std::cout << "      E^2: nnz=" << result.nnz() << std::endl;
	for (size_t i = 3; i <= n; i++) {
		auto next = tensor_shuffle_product_parallel(result, E, F, seq_pool);
		std::cout << "      E^" << i << ": nnz=" << next.nnz() << std::endl;
		result = std::move(next);
	}
	return result;
}

// ========== Expand hepMHV (rank 2) to rank 2L via basis chain ==========
//
// hepMHV has shape (FEC_basis, 11) rank 2. We prepend a dummy dimension
// to make it rank 3: (1, FEC_basis, 11), then use expand_tensor which
// contracts axis 1 with each basis. The result is (1, 11, ..., 11) rank 2L+1.
// Finally, reshape to remove the dummy dimension → (11, ..., 11) rank 2L.
//
// basis_paths: highest weight first, e.g. [first_w3_basis, first_w2_basis] for L=2.

template <typename T, typename index_t>
sparse_tensor<T, index_t, SPARSE_CSR> expand_hepmhv(
	sparse_tensor<T, index_t, SPARSE_CSR>&& hepMHV_csr,
	const std::vector<std::filesystem::path>& basis_paths,
	const field_t& F, thread_pool* pool) {

	// Convert to COO for reshaping (CSR has no reshape method).
	sparse_tensor<T, index_t, SPARSE_COO> hep(std::move(hepMHV_csr));

	// Prepend dummy dimension: (FEC, 11) → (1, FEC, 11)
	std::vector<size_t> new_dims = {1, hep.dim(0), hep.dim(1)};
	hep.reshape(new_dims);

	std::cout << "-- Expand hepMHV --" << std::endl;
	std::cout << "   Starting: rank=" << hep.rank() << " dims=1x" << hep.dim(1) << "x" << hep.dim(2)
	          << " nnz=" << hep.nnz() << std::endl;

	// Convert back to CSR for expand_tensor (which takes CSR input).
	auto hep_csr = sparse_tensor<T, index_t, SPARSE_CSR>(std::move(hep), pool);

	// Use expand_tensor (contracts axis 1)
	auto expanded = expand_tensor<T, index_t>(std::move(hep_csr), basis_paths, F, pool);

	// Remove the dummy dimension: (1, 11, ..., 11) → (11, ..., 11)
	// Convert to COO for reshaping.
	sparse_tensor<T, index_t, SPARSE_COO> expanded_coo(std::move(expanded));
	std::vector<size_t> final_dims;
	for (size_t i = 1; i < expanded_coo.rank(); i++) {
		final_dims.push_back(expanded_coo.dim(i));
	}
	expanded_coo.reshape(final_dims);
	auto result = sparse_tensor<T, index_t, SPARSE_CSR>(std::move(expanded_coo), pool);

	std::cout << "   Final: rank=" << result.rank() << " dims=";
	for (size_t i = 0; i < result.rank(); i++) {
		std::cout << result.dim(i);
		if (i + 1 < result.rank()) std::cout << "x";
	}
	std::cout << " nnz=" << result.nnz() << std::endl;

	return result;
}

// ========== Apply colprojdiv_w1 to each 11-dim letter slot ==========
//
// The collinear constraint requires matching the DIVERGENT part of the boundary.
// Both A (expanded SEW collinear basis) and the boundary live in the full letter
// space (11-dim slots). We project each slot to the 2-dim divergent subspace via
// colprojdiv_w1 (shape (11, 2), (input, output) convention) before solving c.A = b.
//
// After each contraction, the projected axis moves to the end of the tensor, so
// the next slot to project is always at `first_slot_axis`. The final axis order
// is: [non-projected axes..., projected axes...]. For A = (sew_dim, 11^k) this
// gives (sew_dim, 2^k); for boundary = (11^k) this gives (2^k).

template <typename T, typename index_t>
sparse_tensor<T, index_t, SPARSE_COO> apply_colprojdiv_slots(
	sparse_tensor<T, index_t, SPARSE_COO>&& tensor_coo,
	const sparse_tensor<T, index_t, SPARSE_COO>& proj_coo,  // (11, m)
	size_t first_slot_axis,
	size_t n_slots,
	const field_t& F, thread_pool* pool) {

	auto result = std::move(tensor_coo);
	for (size_t s = 0; s < n_slots; s++) {
		result = tensor_contract(result, proj_coo, first_slot_axis, 0, F, pool);
	}
	return result;
}

// ========== Compute the boundary (RHS) for loop order L ==========
//
// Hardcoded formulas:
//   L=2: E1^2 / 2
//   L=3: E1^3 / 6 + E1 * R2
//   L=4: -E1^4 / 12 + E2^2 / 2 + E1 * R3
//   L=5: E1^5 / 20 - E1*E2^2 / 2 + E2*R3 + E1*R4
//
// Inputs: E[1..L-1] and R[2..L-1] as COO tensors.
// Returns the boundary as a COO tensor (rank 2L, dims (11,...,11)).

template <typename T, typename index_t>
sparse_tensor<T, index_t, SPARSE_COO> compute_boundary(
	size_t L,
	const std::map<size_t, sparse_tensor<T, index_t, SPARSE_COO>>& E_list,
	const std::map<size_t, sparse_tensor<T, index_t, SPARSE_COO>>& R_list,
	const field_t& F, thread_pool* pool) {

	std::cout << "== Computing boundary for L=" << L << " ==" << std::endl;

	using coo_t = sparse_tensor<T, index_t, SPARSE_COO>;

	if (L == 2) {
		// boundary = E1^2 / 2
		auto E1_sq = shuffle_power<T, index_t>(E_list.at(1), 2, F, pool);
		std::cout << "   E1^2: nnz=" << E1_sq.nnz() << std::endl;
		auto boundary = tensor_scalar_mul(E1_sq, T(1, 2), F);
		std::cout << "   boundary = E1^2/2: nnz=" << boundary.nnz() << std::endl;
		return boundary;
	}

	if (L == 3) {
		// boundary = E1^3 / 6 + E1 * R2
		auto E1 = E_list.at(1);
		auto R2 = R_list.at(2);

		std::cout << "   Computing E1^3..." << std::endl;
		auto E1_cube = shuffle_power<T, index_t>(E1, 3, F, pool);
		std::cout << "   E1^3: nnz=" << E1_cube.nnz() << std::endl;
		auto term1 = tensor_scalar_mul(E1_cube, T(1, 6), F);
		E1_cube.clear();

		std::cout << "   Computing E1 * R2..." << std::endl;
		auto E1_R2 = tensor_shuffle_product_parallel(E1, R2, F, nullptr);
		std::cout << "   E1*R2: nnz=" << E1_R2.nnz() << std::endl;

		auto boundary = term1;  // start with term1
		tensor_add_weighted(boundary, E1_R2, T(1), T(1), F);  // boundary = term1 + E1*R2
		std::cout << "   boundary = E1^3/6 + E1*R2: nnz=" << boundary.nnz() << std::endl;
		return boundary;
	}

	if (L == 4) {
		// boundary = -E1^4 / 12 + E2^2 / 2 + E1 * R3
		auto E1 = E_list.at(1);
		auto E2 = E_list.at(2);
		auto R3 = R_list.at(3);

		std::cout << "   Computing E1^4..." << std::endl;
		auto E1_4 = shuffle_power<T, index_t>(E1, 4, F, pool);
		auto term1 = tensor_scalar_mul(E1_4, T(-1, 12), F);
		E1_4.clear();

		std::cout << "   Computing E2^2..." << std::endl;
		auto E2_sq = tensor_shuffle_product_parallel(E2, E2, F, nullptr);
		auto term2 = tensor_scalar_mul(E2_sq, T(1, 2), F);
		E2_sq.clear();

		std::cout << "   Computing E1 * R3..." << std::endl;
		auto E1_R3 = tensor_shuffle_product_parallel(E1, R3, F, nullptr);

		// boundary = term1 + term2 + E1*R3
		auto boundary = term1;
		tensor_add_weighted(boundary, term2, T(1), T(1), F);
		tensor_add_weighted(boundary, E1_R3, T(1), T(1), F);
		std::cout << "   boundary = -E1^4/12 + E2^2/2 + E1*R3: nnz=" << boundary.nnz() << std::endl;
		return boundary;
	}

	if (L == 5) {
		// boundary = E1^5 / 20 - E1*E2^2 / 2 + E2*R3 + E1*R4
		auto E1 = E_list.at(1);
		auto E2 = E_list.at(2);
		auto R3 = R_list.at(3);
		auto R4 = R_list.at(4);

		std::cout << "   Computing E1^5..." << std::endl;
		auto E1_5 = shuffle_power<T, index_t>(E1, 5, F, pool);
		auto term1 = tensor_scalar_mul(E1_5, T(1, 20), F);
		E1_5.clear();

		std::cout << "   Computing E1 * E2^2..." << std::endl;
		auto E2_sq = tensor_shuffle_product_parallel(E2, E2, F, nullptr);
		auto E1_E2_sq = tensor_shuffle_product_parallel(E1, E2_sq, F, nullptr);
		E2_sq.clear();
		auto term2 = tensor_scalar_mul(E1_E2_sq, T(-1, 2), F);
		E1_E2_sq.clear();

		std::cout << "   Computing E2 * R3..." << std::endl;
		auto E2_R3 = tensor_shuffle_product_parallel(E2, R3, F, nullptr);

		std::cout << "   Computing E1 * R4..." << std::endl;
		auto E1_R4 = tensor_shuffle_product_parallel(E1, R4, F, nullptr);

		// boundary = term1 + term2 + E2*R3 + E1*R4
		auto boundary = term1;
		tensor_add_weighted(boundary, term2, T(1), T(1), F);
		tensor_add_weighted(boundary, E2_R3, T(1), T(1), F);
		tensor_add_weighted(boundary, E1_R4, T(1), T(1), F);
		std::cout << "   boundary = E1^5/20 - E1*E2^2/2 + E2*R3 + E1*R4: nnz=" << boundary.nnz() << std::endl;
		return boundary;
	}

	throw std::runtime_error("compute_boundary: L=" + std::to_string(L) + " not supported (max 5)");
}

// ========== Helper: locate the dlogmat seed in data_dir ==========
//
// The condition tensor is named dlogmat_<group>.wxf (e.g. dlogmat_E6.wxf).
// For multi-project support we scan data_dir for the first match rather
// than hardcoding the group name.

inline std::filesystem::path find_dlogmat(const std::filesystem::path& data_dir) {
	if (!std::filesystem::exists(data_dir)) {
		throw std::runtime_error("find_dlogmat: data dir not found: " + data_dir.string());
	}
	std::filesystem::path found;
	for (auto& entry : std::filesystem::directory_iterator(data_dir)) {
		auto name = entry.path().filename().string();
		if (name.rfind("dlogmat_", 0) == 0 && entry.path().extension() == ".wxf") {
			found = entry.path();
			break;
		}
	}
	if (found.empty()) {
		throw std::runtime_error("find_dlogmat: no dlogmat_*.wxf found in " + data_dir.string());
	}
	return found;
}

// ========== Helper: run a bootstrap command and log it ==========
//
// All bootstrap subprocess invocations MUST receive --data-dir and --output-dir
// so that --project / --solve-symmetry / --solve-collinear auto-invocations inside
// the subprocess read from / write to the same project directories as compute_rhs.
// The paths must be absolute because the subprocess resolves relative paths
// against its own executable directory (which is the cwd when invoked as
// "./bootstrap").

inline void run_bootstrap_cmd(
	const std::string& cmd,
	const std::filesystem::path& data_dir,
	const std::filesystem::path& output_dir) {

	auto abs_data = std::filesystem::absolute(data_dir);
	auto abs_output = std::filesystem::absolute(output_dir);
	std::string full = cmd
		+ " --data-dir " + abs_data.string()
		+ " --output-dir " + abs_output.string();
	std::cout << "   [bootstrap] " << full << std::endl;
	int ret = std::system(full.c_str());
	if (ret != 0) {
		throw std::runtime_error("bootstrap command failed (code " + std::to_string(ret) + "): " + full);
	}
}

// ========== Ensure FEC tensors exist (bootstrap --extend chain) ==========
//
// FEC_1 is in data/. FEC_{2..F} are generated by bootstrap --extend, stored in output/.
// This generates any missing FEC tensors in the chain.

inline void ensure_fec_tensors(
	size_t F, const std::filesystem::path& data_dir,
	const std::filesystem::path& output_dir) {

	auto dlogmat = find_dlogmat(data_dir);

	auto prev_fec = data_dir / "FEC_1.wxf";
	if (!std::filesystem::exists(prev_fec)) {
		throw std::runtime_error("ensure_fec_tensors: FEC_1 not found: " + prev_fec.string());
	}

	for (size_t w = 2; w <= F; w++) {
		auto curr_fec = output_dir / ("FEC_" + std::to_string(w) + ".wxf");
		if (std::filesystem::exists(curr_fec)) {
			continue;
		}
		std::cout << "   Generating FEC_" << w << " via bootstrap --extend..." << std::endl;
		std::string cmd = "./bootstrap --extend"
			+ std::string(" -c ") + std::filesystem::absolute(dlogmat).string()
			+ std::string(" -f ") + std::filesystem::absolute(prev_fec).string()
			+ std::string(" -o ") + std::filesystem::absolute(curr_fec).string();
		run_bootstrap_cmd(cmd, data_dir, output_dir);
		if (!std::filesystem::exists(curr_fec)) {
			throw std::runtime_error("ensure_fec_tensors: bootstrap did not produce " + curr_fec.string());
		}
		prev_fec = curr_fec;
	}
}

// ========== Ensure SEW tensor and collinear basis exist ==========
//
// The collinear basis (output/collinear/SEW_{name}_basis.wxf) is generated by
// bootstrap --project --symmetry collinear --target SEW_{name}, which also
// generates all FEC collinear bases (first_w{N}_basis.wxf).
// The --project command needs the SEW tensor (output/SEW_{name}.wxf), which
// is generated by bootstrap --sew, which needs FEC_{F} (output/FEC_{F}.wxf).

inline void ensure_sew_basis(
	const std::string& sew_name, size_t F, size_t L,
	const std::filesystem::path& data_dir,
	const std::filesystem::path& output_dir) {

	auto collinear_dir = output_dir / "collinear";
	auto sew_basis_path = collinear_dir / (sew_name + "_basis.wxf");

	if (std::filesystem::exists(sew_basis_path)) {
		std::cout << "   SEW collinear basis already exists: " << sew_basis_path.filename().string() << std::endl;
		return;
	}

	std::cout << "   SEW collinear basis not found, generating pipeline..." << std::endl;

	// Step 1: Ensure FEC_{2..F} exist (via bootstrap --extend)
	ensure_fec_tensors(F, data_dir, output_dir);

	// Step 2: Ensure SEW tensor exists (via bootstrap --sew)
	auto sew_tensor_path = output_dir / (sew_name + ".wxf");
	if (!std::filesystem::exists(sew_tensor_path)) {
		std::cout << "   Generating SEW tensor via bootstrap --sew..." << std::endl;
		auto dlogmat = find_dlogmat(data_dir);
		auto fec_path = output_dir / ("FEC_" + std::to_string(F) + ".wxf");
		auto lec_path = data_dir / ("LEC_" + std::to_string(L) + ".wxf");
		if (!std::filesystem::exists(lec_path)) {
			throw std::runtime_error("ensure_sew_basis: LEC file not found: " + lec_path.string());
		}
		std::string cmd = "./bootstrap --sew"
			+ std::string(" -c ") + std::filesystem::absolute(dlogmat).string()
			+ std::string(" -f ") + std::filesystem::absolute(fec_path).string()
			+ std::string(" -l ") + std::filesystem::absolute(lec_path).string()
			+ std::string(" -o ") + std::filesystem::absolute(sew_tensor_path).string();
		run_bootstrap_cmd(cmd, data_dir, output_dir);
		if (!std::filesystem::exists(sew_tensor_path)) {
			throw std::runtime_error("ensure_sew_basis: bootstrap --sew did not produce " + sew_tensor_path.string());
		}
	}

	// Step 3: Run --project to generate collinear projections and bases
	std::cout << "   Running bootstrap --project --symmetry collinear --target " << sew_name << "..." << std::endl;
	std::string cmd = "./bootstrap --project --symmetry collinear --target " + sew_name;
	run_bootstrap_cmd(cmd, data_dir, output_dir);

	if (!std::filesystem::exists(sew_basis_path)) {
		throw std::runtime_error("ensure_sew_basis: --project did not produce " + sew_basis_path.string());
	}
	std::cout << "   SEW collinear basis generated successfully." << std::endl;
}

// ========== Ensure FEC collinear basis files exist ==========
//
// After ensure_sew_basis runs --project, all FEC collinear bases should exist.
// This is a verification check.

inline void ensure_fec_bases(size_t target_weight, const std::filesystem::path& output_dir) {
	auto collinear_dir = output_dir / "collinear";
	for (size_t w = 2; w < target_weight; w++) {
		auto path = collinear_dir / ("first_w" + std::to_string(w) + "_basis.wxf");
		if (!std::filesystem::exists(path)) {
			throw std::runtime_error("ensure_fec_bases: " + path.string()
				+ " not found (should have been generated by --project).");
		}
	}
}

// ========== Main: compute RHS for a given loop order L ==========
//
// This function recursively computes E[1..L] and R[1..L] (R[1] = E1 - E1 = 0).
// For each loop l from 2 to L:
//   1. Compute boundary_l using E[1..l-1] and R[2..l-1]
//   2. Ensure SEW_{(2l-1)}p1_basis exists (generate if needed)
//   3. Expand SEW_basis → A (rank 2l+1, first dim = sew_dim) [no colprojdiv, per Q3]
//   4. Solve c.A = boundary_l → solMHV_lL (vector of length sew_dim)
//   5. hepMHV_lL = contract(solMHV_lL, SEW_basis, axis 0) [no colprojdiv, per Q3]
//   6. E_l = expand(hepMHV_lL) → rank 2l
//   7. R_l = E_l - boundary_l
//
// Q3: colprojdiv is NOT applied to the SEW basis. The SEW basis IS the collinear
// expression. colprojdiv is only used conceptually (and by the leaf solver
// `bootstrap --solve-collinear` when the user provides a divergent RHS).

template <typename T, typename index_t>
void compute_rhs_for_loop(
	size_t L,
	const std::filesystem::path& data_dir,
	const std::filesystem::path& output_dir,
	const field_t& F, rref_option_t& opt) {

	thread_pool* pool = &(opt->pool);
	auto collinear_dir = output_dir / "collinear";

	std::cout << "========================================" << std::endl;
	std::cout << "Compute RHS for L=" << L << " (SEW_" << (2*L - 1) << "p1)" << std::endl;
	std::cout << "========================================" << std::endl;

	// Ensure E1 exists
	auto oneloop_dir = output_dir / "oneloop";
	std::filesystem::create_directories(oneloop_dir);
	auto E1_path = oneloop_dir / "E1.wxf";
	if (!std::filesystem::exists(E1_path)) {
		auto seed_E1 = data_dir / "E1.wxf";
		if (!std::filesystem::exists(seed_E1)) {
			throw std::runtime_error("compute_rhs: E1 seed not found: " + seed_E1.string());
		}
		std::filesystem::copy_file(seed_E1, E1_path);
		std::cout << "Copied E1 seed from data/E1.wxf" << std::endl;
	}

	// E and R storage (COO tensors, indexed by loop number)
	std::map<size_t, sparse_tensor<T, index_t, SPARSE_COO>> E_list, R_list;

	// Load E1
	{
		auto E1_csr = projection_read_tensor<T, index_t>(E1_path, F, pool);
		sparse_tensor<T, index_t, SPARSE_COO> E1_coo(std::move(E1_csr));
		E_list[1] = std::move(E1_coo);
		std::cout << "E1: rank=" << E_list[1].rank() << " dims=";
		for (size_t i = 0; i < E_list[1].rank(); i++) {
			std::cout << E_list[1].dim(i) << (i + 1 < E_list[1].rank() ? "x" : "");
		}
		std::cout << " nnz=" << E_list[1].nnz() << std::endl;
	}

	// Recursively compute E[2..L-1] and R[2..L-1]
	for (size_t l = 2; l < L; l++) {
		auto l_loop_dir = output_dir / (std::to_string(l) + "loop");
		auto E_l_path = l_loop_dir / ("E" + std::to_string(l) + ".wxf");
		auto R_l_path = l_loop_dir / ("R" + std::to_string(l) + ".wxf");

		if (std::filesystem::exists(E_l_path) && std::filesystem::exists(R_l_path)) {
			std::cout << "L=" << l << ": E" << l << " and R" << l << " already exist, loading." << std::endl;
			auto E_csr = projection_read_tensor<T, index_t>(E_l_path, F, pool);
			auto R_csr = projection_read_tensor<T, index_t>(R_l_path, F, pool);
			E_list[l] = sparse_tensor<T, index_t, SPARSE_COO>(std::move(E_csr));
			R_list[l] = sparse_tensor<T, index_t, SPARSE_COO>(std::move(R_csr));
		} else {
			std::cout << "L=" << l << ": computing recursively..." << std::endl;
			compute_rhs_for_loop<T, index_t>(l, data_dir, output_dir, F, opt);
			// Reload
			auto E_csr = projection_read_tensor<T, index_t>(E_l_path, F, pool);
			auto R_csr = projection_read_tensor<T, index_t>(R_l_path, F, pool);
			E_list[l] = sparse_tensor<T, index_t, SPARSE_COO>(std::move(E_csr));
			R_list[l] = sparse_tensor<T, index_t, SPARSE_COO>(std::move(R_csr));
		}
	}

	// Now compute for loop L
	auto L_loop_dir = output_dir / (std::to_string(L) + "loop");
	std::filesystem::create_directories(L_loop_dir);

	std::string sew_name = "SEW_" + std::to_string(2*L - 1) + "p1";
	size_t F_weight = 2*L - 1;
	size_t L_weight = 1;
	size_t target_weight = F_weight + L_weight;  // = 2L

	// Step 1: Compute boundary_L
	auto boundary = compute_boundary<T, index_t>(L, E_list, R_list, F, pool);
	auto boundary_path = L_loop_dir / ("boundary_" + std::to_string(L) + "L.wxf");
	{
		auto boundary_csr = sparse_tensor<T, index_t, SPARSE_CSR>(std::move(boundary), pool);
		write_tensor_file<T, index_t>(boundary_path, std::move(boundary_csr));
	}

	// Step 2: Ensure SEW basis exists
	std::cout << "== Ensuring SEW basis ==" << std::endl;
	ensure_sew_basis(sew_name, F_weight, L_weight, data_dir, output_dir);

	// Step 3: Ensure FEC bases exist
	ensure_fec_bases(target_weight, output_dir);

	// Step 4: Compute colprojdiv_SEW via collinear projection chain (for the leaf
	// solver `bootstrap --solve-collinear` — NOT applied to A or hepMHV per Q3).
	std::cout << "== Computing collinear projections (for leaf solver) ==" << std::endl;
	auto sew_basis_path = collinear_dir / (sew_name + "_basis.wxf");
	auto colprojdiv_path = collinear_dir / ("colprojdiv_" + sew_name + ".wxf");
	auto colprojfin_path = collinear_dir / ("colprojfin_" + sew_name + ".wxf");

	if (!std::filesystem::exists(colprojdiv_path)) {
		// Build chain base paths (weights 2..target_weight)
		std::vector<std::filesystem::path> chain_bases;
		for (size_t w = 2; w <= target_weight; w++) {
			if (w == target_weight) {
				chain_bases.push_back(sew_basis_path);  // SEW level
			} else {
				chain_bases.push_back(collinear_dir / ("first_w" + std::to_string(w) + "_basis.wxf"));
			}
		}
		run_collinear_proj_chain<T, index_t>(chain_bases, data_dir, output_dir, F, opt, sew_name);
	} else {
		std::cout << "   colprojdiv_" << sew_name << " already exists." << std::endl;
	}

	// Step 5: Load SEW basis (the collinear expression). Per Q3, do NOT apply
	// colprojdiv — the SEW basis IS the collinear expression. colprojdiv at the
	// SEW level is identity anyway, and applying it would be conceptually wrong.
	std::cout << "== Loading SEW basis (collinear expression, no colprojdiv) ==" << std::endl;
	auto sew_basis_csr = projection_read_tensor<T, index_t>(sew_basis_path, F, pool);

	// Step 6: Expand SEW basis → A (rank 2L+1)
	std::cout << "== Expanding SEW tensor ==" << std::endl;
	std::vector<std::filesystem::path> expansion_bases;
	for (long w = static_cast<long>(target_weight) - 1; w >= 2; w--) {
		expansion_bases.push_back(collinear_dir / ("first_w" + std::to_string(w) + "_basis.wxf"));
	}
	auto A_expanded = expand_tensor<T, index_t>(std::move(sew_basis_csr), expansion_bases, F, pool);

	std::cout << "   A: rank=" << A_expanded.rank() << " dims=";
	for (size_t i = 0; i < A_expanded.rank(); i++) {
		std::cout << A_expanded.dim(i) << (i + 1 < A_expanded.rank() ? "x" : "");
	}
	std::cout << std::endl;

	// Step 7: Solve c.A = boundary at matching positions.
	//
	// A (expanded SEW basis) and boundary (E1^L / L!) may have different supports.
	// We solve using only matching positions (both nonzero), where the ratio c =
	// boundary / A should be consistent. Non-matching positions become part of R*.
	// The indicator-vector method (Step 13) verifies R* = c.A - boundary is
	// divergent-free.
	std::cout << "== Solving at matching positions ==" << std::endl;

	// Load colprojdiv_w1 (needed for indicator-vector verification in Step 13)
	auto colprojdiv_w1_path = collinear_dir / "colprojdiv_w1.wxf";
	if (!std::filesystem::exists(colprojdiv_w1_path)) {
		throw std::runtime_error("compute_rhs: colprojdiv_w1 not found: " + colprojdiv_w1_path.string()
			+ " (run the collinear projection chain first)");
	}
	auto colprojdiv_w1_csr = projection_read_tensor<T, index_t>(colprojdiv_w1_path, F, pool);
	sparse_tensor<T, index_t, SPARSE_COO> colprojdiv_w1_coo(std::move(colprojdiv_w1_csr));
	std::cout << "   colprojdiv_w1: rank=" << colprojdiv_w1_coo.rank() << " dims=";
	for (size_t i = 0; i < colprojdiv_w1_coo.rank(); i++) {
		std::cout << colprojdiv_w1_coo.dim(i) << (i + 1 < colprojdiv_w1_coo.rank() ? "x" : "");
	}
	std::cout << " nnz=" << colprojdiv_w1_coo.nnz() << std::endl;

	// Load boundary and A
	auto boundary_csr = projection_read_tensor<T, index_t>(boundary_path, F, pool);
	sparse_tensor<T, index_t, SPARSE_CSR> A_csr(std::move(A_expanded));

	std::cout << "   A: rank=" << A_csr.rank() << " dims=";
	for (size_t i = 0; i < A_csr.rank(); i++) {
		std::cout << A_csr.dim(i) << (i + 1 < A_csr.rank() ? "x" : "");
	}
	std::cout << " nnz=" << A_csr.nnz() << std::endl;
	std::cout << "   boundary: rank=" << boundary_csr.rank() << " dims=";
	for (size_t i = 0; i < boundary_csr.rank(); i++) {
		std::cout << boundary_csr.dim(i) << (i + 1 < boundary_csr.rank() ? "x" : "");
	}
	std::cout << " nnz=" << boundary_csr.nnz() << std::endl;

	// Convert to COO for index inspection
	sparse_tensor<T, index_t, SPARSE_COO> A_coo(std::move(A_csr));
	sparse_tensor<T, index_t, SPARSE_COO> b_coo(std::move(boundary_csr));

	// Project both A and boundary to the DIVERGENT subspace by applying colprojdiv_w1
	// to each letter slot. The collinear constraint only requires c.A = boundary in the
	// divergent subspace; the finite part is R* (the finite remainder). Solving in the
	// full space fails at L>=3 because c.A and boundary legitimately differ at finite
	// positions. After projection, each 11-dim letter slot becomes 2-dim.
	size_t n_letter_slots = b_coo.rank();          // = 2L
	size_t A_first_slot = 1;                       // A has leading sew_dim axis
	std::cout << "== Projecting A and boundary to divergent subspace ==" << std::endl;
	std::cout << "   Projecting A (n_slots=" << n_letter_slots << ", first_slot=" << A_first_slot << ")..." << std::endl;
	A_coo = apply_colprojdiv_slots<T, index_t>(std::move(A_coo), colprojdiv_w1_coo, A_first_slot, n_letter_slots, F, pool);
	std::cout << "   A_proj: rank=" << A_coo.rank() << " dims=";
	for (size_t i = 0; i < A_coo.rank(); i++) std::cout << A_coo.dim(i) << (i + 1 < A_coo.rank() ? "x" : "");
	std::cout << " nnz=" << A_coo.nnz() << std::endl;

	std::cout << "   Projecting boundary (n_slots=" << n_letter_slots << ", first_slot=0)..." << std::endl;
	b_coo = apply_colprojdiv_slots<T, index_t>(std::move(b_coo), colprojdiv_w1_coo, 0, n_letter_slots, F, pool);
	std::cout << "   b_proj: rank=" << b_coo.rank() << " dims=";
	for (size_t i = 0; i < b_coo.rank(); i++) std::cout << b_coo.dim(i) << (i + 1 < b_coo.rank() ? "x" : "");
	std::cout << " nnz=" << b_coo.nnz() << std::endl;

	// Build b_map: letter key -> boundary value (boundary has no sew axis).
	// For A, iterate over all entries and match against b_map. Each A entry
	// keeps its full index (including the sew axis) so that positions with
	// BOTH sew entries contribute two coefficients to the linear system
	// c0*A[0,key] + c1*A[1,key] = b[key]. The old code used std::map<key,T>
	// which overwrote one sew entry — broken for sew_dim > 1.
	std::map<std::vector<index_t>, T> b_map;
	for (auto i : b_coo.gen_perm()) {
		b_map[b_coo.index_vector(i)] = b_coo.val(i);
	}

	// Create A_match and b_match: A_match preserves the sew axis (all matching
	// A entries with their original sew index); b_match has one entry per matching
	// position.
	sparse_tensor<T, index_t, SPARSE_COO> A_match(A_coo.dims());
	sparse_tensor<T, index_t, SPARSE_COO> b_match(b_coo.dims());
	std::set<std::vector<index_t>> matched_keys;
	size_t n_A_entries = 0;
	for (auto i : A_coo.gen_perm()) {
		auto full_idx = A_coo.index_vector(i);              // 7 elements [sew, p1..p6]
		std::vector<index_t> key(full_idx.begin() + 1, full_idx.end());  // [p1..p6]
		auto it = b_map.find(key);
		if (it != b_map.end()) {
			A_match.push_back(full_idx, A_coo.val(i));      // preserve sew index
			if (matched_keys.insert(key).second) {
				b_match.push_back(key, it->second);          // push b once per key
			}
			n_A_entries++;
		}
	}
	A_match.canonicalize();
	b_match.canonicalize();

	size_t n_match = matched_keys.size();
	size_t n_A_only = 0, n_b_only = 0;
	// Count A_only: A entries whose key is not in b_map
	for (auto i : A_coo.gen_perm()) {
		auto full_idx = A_coo.index_vector(i);
		std::vector<index_t> key(full_idx.begin() + 1, full_idx.end());
		if (b_map.find(key) == b_map.end()) n_A_only++;
	}
	n_b_only = b_map.size() - n_match;

	std::cout << "   Matching positions: " << n_match
	          << " (A_entries=" << n_A_entries
	          << ", A_only=" << n_A_only
	          << ", b_only=" << n_b_only << ")" << std::endl;

	// Convert to CSR for solve_linear_system
	auto A_match_csr = sparse_tensor<T, index_t, SPARSE_CSR>(std::move(A_match), pool);
	auto b_match_csr = sparse_tensor<T, index_t, SPARSE_CSR>(std::move(b_match), pool);

	// Step 8: Solve c.A_match = b_match
	auto solve_result = solve_linear_system<T, index_t>(
		std::move(A_match_csr), std::move(b_match_csr), F, opt);

	if (!solve_result.consistent) {
		std::cout << "========================================" << std::endl;
		std::cout << "FAILED: System is inconsistent — no solution exists." << std::endl;
		std::cout << "   Ratios at matching positions are not unique." << std::endl;
		std::cout << "========================================" << std::endl;
		throw std::runtime_error("compute_rhs: collinear constraints inconsistent at L=" + std::to_string(L));
	}

	sparse_vec<T, index_t> solution_vec = std::move(solve_result.solution);
	size_t n_unknowns = solve_result.n_unknowns;

	// Step 9: Build sol_mat (1, n_unknowns) once and reuse for both writing and contracting.
	auto solMHV_path = L_loop_dir / ("solMHV_" + std::to_string(L) + "L.wxf");
	auto hepMHV_path = L_loop_dir / ("hepMHV_" + std::to_string(L) + "L.wxf");
	auto E_L_path = L_loop_dir / ("E" + std::to_string(L) + ".wxf");
	auto R_L_path = L_loop_dir / ("R" + std::to_string(L) + ".wxf");
	{
		sparse_mat<T, index_t> sol_mat(1, n_unknowns);
		for (size_t i = 0; i < solution_vec.nnz(); i++) {
			sol_mat[0].push_back(solution_vec(i), solution_vec[i]);
		}
		sol_mat[0].compress();

		std::cout << "   solMHV_" << L << "L: " << solution_vec.nnz()
		          << " non-zero coefficients (n_unknowns=" << n_unknowns << ")" << std::endl;
		for (size_t i = 0; i < solution_vec.nnz(); i++) {
			std::cout << "      c[" << solution_vec(i) << "] = " << solution_vec[i] << std::endl;
		}

		// Write solMHV_LL as a (1, n_unknowns) rank-2 CSR tensor.
		sparse_tensor<T, index_t, SPARSE_CSR> sol_csr(sol_mat);
		write_tensor_file<T, index_t>(solMHV_path, std::move(sol_csr));

		// Step 10: Compute hepMHV_LL = contract(solMHV_LL, SEW_basis, axis 1, 0)
		// Per Q3: use the SEW basis directly (no colprojdiv). The SEW basis IS the
		// collinear expression; colprojdiv is not applied to hepMHV/E_L.
		// sol_csr has logical dims (1, n_unknowns); we contract axis 1 (the n_unknowns
		// axis, which is the dummy "1" + actual vector) with axis 0 of the SEW basis.
		std::cout << "== Computing hepMHV (no colprojdiv, per Q3) ==" << std::endl;
		// Reload SEW basis (was consumed by expand_tensor in Step 6)
		auto sew_basis_csr2 = projection_read_tensor<T, index_t>(sew_basis_path, F, pool);

		// Rebuild sol_csr for the contraction (write_tensor_file moved it)
		sparse_tensor<T, index_t, SPARSE_CSR> sol_csr2(sol_mat);
		sparse_tensor<T, index_t, SPARSE_COO> sol_coo(std::move(sol_csr2));
		sparse_tensor<T, index_t, SPARSE_COO> sew_coo(std::move(sew_basis_csr2));

		std::cout << "   sol: rank=" << sol_coo.rank() << " dims=";
		for (size_t i = 0; i < sol_coo.rank(); i++) {
			std::cout << sol_coo.dim(i) << (i + 1 < sol_coo.rank() ? "x" : "");
		}
		std::cout << std::endl;
		std::cout << "   SEW basis: rank=" << sew_coo.rank() << " dims=";
		for (size_t i = 0; i < sew_coo.rank(); i++) {
			std::cout << sew_coo.dim(i) << (i + 1 < sew_coo.rank() ? "x" : "");
		}
		std::cout << std::endl;

		// Contract axis 1 of sol (n_unknowns) with axis 0 of SEW basis (sew_dim).
		// Result axes: [sol.0=1, sew.1, sew.2, ...] → rank-3 (1, FEC, LEC).
		// The leading "1" is a dummy dimension from the solution vector's row shape.
		// Reshape to (FEC, LEC) rank-2 so expand_hepmhv can prepend its own dummy "1".
		auto hepMHV = tensor_contract(sol_coo, sew_coo, 1, 0, F, pool);
		std::cout << "   hepMHV (raw): rank=" << hepMHV.rank() << " dims=";
		for (size_t i = 0; i < hepMHV.rank(); i++) {
			std::cout << hepMHV.dim(i) << (i + 1 < hepMHV.rank() ? "x" : "");
		}
		std::cout << " nnz=" << hepMHV.nnz() << std::endl;

		// Remove the leading "1": (1, FEC, LEC) → (FEC, LEC)
		std::vector<size_t> hepMHV_dims;
		for (size_t i = 1; i < hepMHV.rank(); i++) {
			hepMHV_dims.push_back(hepMHV.dim(i));
		}
		hepMHV.reshape(hepMHV_dims);
		auto hepMHV_csr = sparse_tensor<T, index_t, SPARSE_CSR>(std::move(hepMHV), pool);

		std::cout << "   hepMHV: rank=" << hepMHV_csr.rank() << " dims=";
		for (size_t i = 0; i < hepMHV_csr.rank(); i++) {
			std::cout << hepMHV_csr.dim(i) << (i + 1 < hepMHV_csr.rank() ? "x" : "");
		}
		std::cout << " nnz=" << hepMHV_csr.nnz() << std::endl;

		write_tensor_file<T, index_t>(hepMHV_path, std::move(hepMHV_csr));
	}

	// Step 11: Expand hepMHV → E_L
	std::cout << "== Expanding hepMHV to E" << L << " ==" << std::endl;
	auto hepMHV_csr2 = projection_read_tensor<T, index_t>(hepMHV_path, F, pool);
	auto E_L = expand_hepmhv<T, index_t>(std::move(hepMHV_csr2), expansion_bases, F, pool);
	write_tensor_file<T, index_t>(E_L_path, std::move(E_L));

	// Step 12: Compute R_L = E_L - boundary
	std::cout << "== Computing R" << L << " ==" << std::endl;
	auto E_L_csr2 = projection_read_tensor<T, index_t>(E_L_path, F, pool);
	auto boundary_csr2 = projection_read_tensor<T, index_t>(boundary_path, F, pool);
	sparse_tensor<T, index_t, SPARSE_COO> E_L_coo(std::move(E_L_csr2));
	sparse_tensor<T, index_t, SPARSE_COO> boundary_coo(std::move(boundary_csr2));
	tensor_add_weighted(E_L_coo, boundary_coo, T(1), T(-1), F);  // E_L = E_L - boundary
	auto R_L_csr = sparse_tensor<T, index_t, SPARSE_CSR>(std::move(E_L_coo), pool);
	write_tensor_file<T, index_t>(R_L_path, std::move(R_L_csr));

	// Step 13: Verify R* is free of divergent letters (indicator-vector method).
	// R_L = c·A - boundary. The collinear constraint requires the divergent part to
	// match; R_L is the finite remainder R* if and only if no divergent letter appears
	// in its support. Method: collect all distinct letter indices from R_L's nonzero
	// entries, form an 11-dim indicator vector v (v[k]=1 if letter k appears), project
	// v with colprojdiv_w1 (11, m_div). If the projection is zero → R_L = R*.
	std::cout << "== Verifying R* is divergent-free (indicator-vector method) ==" << std::endl;
	auto R_L_verify_csr = projection_read_tensor<T, index_t>(R_L_path, F, pool);
	sparse_tensor<T, index_t, SPARSE_COO> R_L_verify(std::move(R_L_verify_csr));

	// Collect distinct letter indices from R_L's nonzero entries
	std::set<size_t> letter_indices;
	size_t r_rank = R_L_verify.rank();
	for (size_t i = 0; i < R_L_verify.nnz(); i++) {
		auto coords = R_L_verify.index(i);  // pointer to r_rank index_t values
		for (size_t d = 0; d < r_rank; d++) {
			letter_indices.insert(static_cast<size_t>(coords[d]));
		}
	}
	std::cout << "   R_L nnz=" << R_L_verify.nnz()
	          << ", distinct letter indices: " << letter_indices.size() << std::endl;
	std::cout << "   Letters appearing: {";
	{
		bool first = true;
		for (auto idx : letter_indices) {
			if (!first) std::cout << ",";
			std::cout << idx;
			first = false;
		}
	}
	std::cout << "}" << std::endl;

	// Build indicator vector (rank-1 COO tensor, dim 11)
	sparse_tensor<T, index_t, SPARSE_COO> indicator({11});
	indicator.resize(letter_indices.size());
	size_t pos = 0;
	for (auto idx : letter_indices) {
		indicator.index(pos)[0] = static_cast<index_t>(idx);
		indicator.val(pos) = T(1);
		pos++;
	}
	indicator.canonicalize();
	indicator.sort_indices();
	std::cout << "   Indicator vector (11-dim): nnz=" << indicator.nnz() << std::endl;

	// Project: indicator (11) · colprojdiv_w1 (11, m_div) → (m_div,)
	auto div_check = tensor_contract(indicator, colprojdiv_w1_coo, 0, 0, F, pool);
	std::cout << "   Projection to divergent subspace: rank=" << div_check.rank()
	          << " dims=";
	for (size_t i = 0; i < div_check.rank(); i++) {
		std::cout << div_check.dim(i) << (i + 1 < div_check.rank() ? "x" : "");
	}
	std::cout << " nnz=" << div_check.nnz() << std::endl;

	if (div_check.nnz() != 0) {
		std::cout << "========================================" << std::endl;
		std::cout << "FAILED: R* contains divergent letters." << std::endl;
		std::cout << "   The indicator vector projected to a non-zero vector in the" << std::endl;
		std::cout << "   divergent subspace — R_L is NOT free of divergent letters." << std::endl;
		std::cout << "========================================" << std::endl;
		throw std::runtime_error("compute_rhs: R* divergent-letter check failed at L="
			+ std::to_string(L));
	}
	std::cout << "   [OK] R* is free of divergent letters — R_L = R*" << std::endl;

	std::cout << "========================================" << std::endl;
	std::cout << "Done: L=" << L << std::endl;
	std::cout << "   E" << L << ": " << (L_loop_dir / ("E" + std::to_string(L) + ".wxf")).string() << std::endl;
	std::cout << "   R" << L << ": " << (L_loop_dir / ("R" + std::to_string(L) + ".wxf")).string() << std::endl;
	std::cout << "   boundary_" << L << "L: " << boundary_path.string() << std::endl;
	std::cout << "   solMHV_" << L << "L: " << (L_loop_dir / ("solMHV_" + std::to_string(L) + "L.wxf")).string() << std::endl;
	std::cout << "   hepMHV_" << L << "L: " << (L_loop_dir / ("hepMHV_" + std::to_string(L) + "L.wxf")).string() << std::endl;
	std::cout << "========================================" << std::endl;
}

#endif // COMPUTE_RHS_HPP
