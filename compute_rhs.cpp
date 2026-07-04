// compute_rhs.cpp — Standalone executable for computing the RHS (collinear boundary)
//
// This module recursively computes E[L] and R[L] from loop 2 up to the target
// loop order. It invokes the main bootstrap module (--sew, --extend) to generate
// missing SEW basis files, and uses the collinear solver logic internally for
// each loop order.
//
// Usage:
//   ./compute_rhs --target <SEW_FpL> [--data-dir <dir>] [--output-dir <dir>]
//
// Examples:
//   ./compute_rhs --target SEW_3p1    (2-loop: computes E2, R2, boundary_2L)
//   ./compute_rhs --target SEW_5p1    (3-loop: computes E3, R3, boundary_3L, requires E2/R2)
//
// The loop order L is derived from the SEW target: L = (F + L) / 2.
// For SEW_3p1: L = (3+1)/2 = 2. For SEW_5p1: L = (5+1)/2 = 3.

#include "bootstrap.hpp"
#include "projection.hpp"
#include "compute_rhs.hpp"

#include <cstdlib>

using index_t = int32_t;
using scalar_t = rat_t;

struct rhs_args_t {
	std::string target;
	std::filesystem::path data_dir;
	std::filesystem::path output_dir;
	bool help = false;
};

void print_usage(const char* program) {
	std::cerr << "Usage:" << std::endl;
	std::cerr << "  " << program << " --target <SEW_FpL> [--data-dir <dir>] [--output-dir <dir>]" << std::endl;
	std::cerr << std::endl;
	std::cerr << "Options:" << std::endl;
	std::cerr << "  --target <SEW_FpL>  Target SEW name (e.g. SEW_3p1 for 2-loop, SEW_5p1 for 3-loop)" << std::endl;
	std::cerr << "  --data-dir <dir>    Data directory with seed files (default: <exec_dir>/data)" << std::endl;
	std::cerr << "  --output-dir <dir>  Output directory (default: <exec_dir>/output)" << std::endl;
	std::cerr << "  -h, --help          Show this help message" << std::endl;
}

std::string take_value(int& i, int argc, char* argv[], const std::string& flag) {
	if (i + 1 >= argc) {
		throw std::runtime_error("Missing value after " + flag);
	}
	i++;
	return argv[i];
}

rhs_args_t parse_args(int argc, char* argv[]) {
	rhs_args_t args;
	for (int i = 1; i < argc; i++) {
		std::string arg = argv[i];
		if (arg == "--target") {
			args.target = take_value(i, argc, argv, arg);
		}
		else if (arg == "--data-dir") {
			args.data_dir = take_value(i, argc, argv, arg);
		}
		else if (arg == "--output-dir") {
			args.output_dir = take_value(i, argc, argv, arg);
		}
		else if (arg == "-h" || arg == "--help") {
			args.help = true;
		}
		else {
			throw std::runtime_error("Unknown argument: " + arg);
		}
	}
	return args;
}

int main(int argc, char* argv[]) {
	try {
		rhs_args_t args = parse_args(argc, argv);

		if (args.help || args.target.empty()) {
			print_usage(argv[0]);
			return args.help ? 0 : 1;
		}

		// Parse the target to determine the loop order
		auto target = parse_target(args.target);
		if (target.kind != target_kind_t::SEW) {
			throw std::runtime_error("compute_rhs only supports SEW targets, got: " + args.target);
		}

		// L = (F + L) / 2 where F and L are the FEC and LEC weights
		size_t L = (target.fec_weight + target.lec_weight) / 2;
		if ((target.fec_weight + target.lec_weight) % 2 != 0) {
			throw std::runtime_error("compute_rhs: target weight " + args.target
				+ " is odd — expected even (F+L=2*loop_order)");
		}
		if (L < 2) {
			throw std::runtime_error("compute_rhs: loop order L=" + std::to_string(L)
				+ " < 2 — nothing to compute (E1 is the seed)");
		}
		if (L > 5) {
			throw std::runtime_error("compute_rhs: loop order L=" + std::to_string(L)
				+ " > 5 — boundary formulas only implemented up to L=5");
		}

		// Resolve paths relative to the executable directory
		std::filesystem::path base = std::filesystem::path(argv[0]).parent_path();
		if (base.empty()) {
			base = std::filesystem::current_path();
		}

		// Resolve paths relative to the executable directory.
	// Matches bootstrap.cpp/inspect_tensors.cpp: relative --data-dir/--output-dir
	// are resolved against `base` (the executable's parent directory), not cwd.
	// This keeps "./compute_rhs --data-dir data_test" working from any cwd as long
	// as the binary is invoked from the repo root (the common case).
	if (args.data_dir.empty()) {
		args.data_dir = base / "data";
	} else if (!args.data_dir.is_absolute()) {
		args.data_dir = base / args.data_dir;
	}
	if (args.output_dir.empty()) {
		args.output_dir = base / "output";
	} else if (!args.output_dir.is_absolute()) {
		args.output_dir = base / args.output_dir;
	}

		std::cout << "========================================" << std::endl;
		std::cout << "Compute RHS (collinear boundary)" << std::endl;
		std::cout << "   Target: " << args.target << " (loop order L=" << L << ")" << std::endl;
		std::cout << "   Data dir: " << args.data_dir.string() << std::endl;
		std::cout << "   Output dir: " << args.output_dir.string() << std::endl;
		std::cout << "========================================" << std::endl;

		field_t F(FIELD_QQ);

		rref_option_t opt;
		opt->method = 0;
		opt->verbose = true;
		opt->pool.reset();
		thread_pool* pool = &(opt->pool);

		std::cout << "threads: " << pool->get_thread_count() << std::endl;
		auto local_time = std::chrono::zoned_time{
			std::chrono::current_zone(),
			std::chrono::system_clock::now()
		};
		std::cout << "Task begin at: " << std::format("{:%Y-%m-%d %H:%M:%S %z}", local_time) << std::endl;

		Timer total_timer;
		total_timer.start();

		compute_rhs_for_loop<scalar_t, index_t>(
			L, args.data_dir, args.output_dir, F, opt);

		total_timer.stop();
		local_time = std::chrono::zoned_time{
			std::chrono::current_zone(),
			std::chrono::system_clock::now()
		};
		std::cout << "Task end at: " << std::format("{:%Y-%m-%d %H:%M:%S %z}", local_time) << std::endl;
		std::cout << "Done in " << total_timer.milliseconds() << " ms" << std::endl << std::endl;
	}
	catch (const std::exception& e) {
		std::cerr << "Error: " << e.what() << std::endl;
		print_usage(argv[0]);
		return 1;
	}
	return 0;
}
