#include "graphblas/tc.hpp"
#include "graphblas/utils.hpp"
#include "spla/tc.hpp"
#include "spla/utils.hpp"

#include <GraphBLAS.h>
#include <spla.hpp>

#include <chrono>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

enum class Backend { GraphBLAS, SplaCpu, SplaGpu };
enum class Algorithm { Burkhardt, Sandia };

void check_graphblas(GrB_Info status, const char *operation) {
  if (status != GrB_SUCCESS) {
    throw std::runtime_error(std::string(operation) +
                             " failed with GraphBLAS status " +
                             std::to_string(status));
  }
}

int parse_iterations(const char *value) {
  std::size_t parsed = 0;
  const std::string text(value);
  const long count = std::stol(text, &parsed);
  if (parsed != text.size() || count <= 0) {
    throw std::invalid_argument("iterations must be a positive integer");
  }
  return static_cast<int>(count);
}

std::uint64_t run_graphblas(const std::string &dataset, Algorithm algorithm,
                            int iterations) {
  check_graphblas(GrB_init(GrB_NONBLOCKING), "GrB_init");

  GrB_Matrix adjacency = nullptr;
  GrB_Matrix workspace = nullptr;
  try {
    const bool triangular = algorithm == Algorithm::Sandia;
    adjacency = graphblas_utils::load_graph(dataset, triangular);

    GrB_Index n = 0;
    check_graphblas(GrB_Matrix_nrows(&n, adjacency), "GrB_Matrix_nrows");
    check_graphblas(GrB_Matrix_new(&workspace, GrB_UINT64, n, n),
                    "GrB_Matrix_new");

    const auto count_once = [&]() {
      return algorithm == Algorithm::Sandia
                 ? tc_graphblas::sandia(adjacency, workspace)
                 : tc_graphblas::burkhardt(adjacency, workspace);
    };

    // Populate the JIT cache and internal matrix formats before the profiled
    // repetitions. perf still sees this call, but its cost is amortized by the
    // requested repetitions.
    const std::uint64_t expected = count_once();
    check_graphblas(GrB_Matrix_clear(workspace), "GrB_Matrix_clear");

    for (int i = 0; i < iterations; ++i) {
      const std::uint64_t actual = count_once();
      if (actual != expected) {
        throw std::runtime_error("triangle count changed between iterations");
      }
      check_graphblas(GrB_Matrix_clear(workspace), "GrB_Matrix_clear");
    }

    check_graphblas(GrB_Matrix_free(&workspace), "GrB_Matrix_free(workspace)");
    check_graphblas(GrB_Matrix_free(&adjacency), "GrB_Matrix_free(adjacency)");
    check_graphblas(GrB_finalize(), "GrB_finalize");
    return expected;
  } catch (...) {
    if (workspace != nullptr) {
      GrB_Matrix_free(&workspace);
    }
    if (adjacency != nullptr) {
      GrB_Matrix_free(&adjacency);
    }
    GrB_finalize();
    throw;
  }
}

std::uint64_t run_spla(const std::string &dataset, Algorithm algorithm,
                       int iterations, bool accelerated) {
  spla::Library::get()->set_force_no_acceleration(!accelerated);
  const bool triangular = algorithm == Algorithm::Sandia;
  auto adjacency = spla_utils::load_graph(dataset, triangular);
  auto workspace = spla::Matrix::make(adjacency->get_n_rows(),
                                      adjacency->get_n_cols(), spla::INT);

  const auto count_once = [&]() {
    int count = 0;
    if (algorithm == Algorithm::Sandia) {
      tc_spla::sandia(count, adjacency, workspace);
    } else {
      tc_spla::burkhardt(count, adjacency, workspace);
    }
    return static_cast<std::uint64_t>(count);
  };

  const std::uint64_t expected = count_once();
  workspace->clear();

  for (int i = 0; i < iterations; ++i) {
    const std::uint64_t actual = count_once();
    if (actual != expected) {
      throw std::runtime_error("triangle count changed between iterations");
    }
    workspace->clear();
  }

  return expected;
}

} // namespace

int main(int argc, char *argv[]) {
  if (argc != 5) {
    std::cerr << "Usage: " << argv[0]
              << " <gb|spla-cpu|spla-gpu> <burkhardt|sandia> <dataset>"
                 " <iterations>\n";
    return 1;
  }

  try {
    const std::string backend_arg(argv[1]);
    const std::string algorithm_arg(argv[2]);
    const std::string dataset(argv[3]);
    const int iterations = parse_iterations(argv[4]);

    const Backend backend =
        backend_arg == "gb"
            ? Backend::GraphBLAS
            : (backend_arg == "spla" || backend_arg == "spla-cpu")
                  ? Backend::SplaCpu
                  : backend_arg == "spla-gpu"
                        ? Backend::SplaGpu
                        : throw std::invalid_argument(
                              "backend must be 'gb', 'spla-cpu', or "
                              "'spla-gpu'");
    const Algorithm algorithm =
        algorithm_arg == "burkhardt"
            ? Algorithm::Burkhardt
            : algorithm_arg == "sandia"
                  ? Algorithm::Sandia
                  : throw std::invalid_argument(
                        "algorithm must be 'burkhardt' or 'sandia'");

    std::cout << "backend=" << backend_arg << " algorithm=" << algorithm_arg
              << " dataset=" << dataset << " warmup=1"
              << " iterations=" << iterations << '\n';

    const auto start = std::chrono::steady_clock::now();
    const std::uint64_t count =
        backend == Backend::GraphBLAS
            ? run_graphblas(dataset, algorithm, iterations)
            : run_spla(dataset, algorithm, iterations,
                       backend == Backend::SplaGpu);
    const auto end = std::chrono::steady_clock::now();

    std::cout << "triangles=" << count << '\n';
    std::cout << "total_time_s="
              << std::chrono::duration<double>(end - start).count() << '\n';
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "Error: " << error.what() << '\n';
    return 2;
  }
}
