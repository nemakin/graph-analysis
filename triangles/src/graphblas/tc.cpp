#include "tc.hpp"
#include "utils.hpp"
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace tc_graphblas {
uint64_t triangles_counting(GrB_Matrix A, GrB_Matrix workspace,
                            bool triangular) {
  GrB_mxm(workspace, A, nullptr, GxB_PLUS_TIMES_UINT64, A, A, nullptr);

  uint64_t sum = 0;
  GrB_Matrix_reduce_UINT64(&sum, nullptr, GrB_PLUS_MONOID_UINT64, workspace,
                           nullptr);

  return triangular ? sum : sum / 6;
}

uint64_t burkhardt(GrB_Matrix A, GrB_Matrix workspace) {
  return triangles_counting(A, workspace, false);
}

uint64_t sandia(GrB_Matrix A, GrB_Matrix workspace) {
  return triangles_counting(A, workspace, true);
}

std::vector<double> benchmark(const char *filename, bool triangular,
                              const int num_iters) {
  std::vector<double> iteration_times;
  iteration_times.reserve(num_iters);

  GrB_init(GrB_NONBLOCKING);

  GrB_Matrix A;
  A = graphblas_utils::load_graph(filename, triangular);
  GrB_Index n;
  GrB_Matrix_nrows(&n, A);
  GrB_Matrix workspace;
  GrB_Matrix_new(&workspace, GrB_UINT64, n, n);

  for (int i = 0; i < num_iters; ++i) {
    auto start = std::chrono::high_resolution_clock::now();
    uint64_t answer;
    if (triangular) {
      answer = sandia(A, workspace);
    } else {
      answer = burkhardt(A, workspace);
    }
    auto end = std::chrono::high_resolution_clock::now();

    GrB_Matrix_clear(workspace);

    std::chrono::duration<double> elapsed = end - start;

    std::cout << (triangular ? "GB_Sandia" : "GB_Burkhardt") << " Iteration "
              << i + 1 << ": " << elapsed.count() << " s" << std::endl;
    iteration_times.push_back(elapsed.count());
  }

  GrB_Matrix_free(&workspace);
  GrB_Matrix_free(&A);
  GrB_finalize();

  return iteration_times;
}

} // namespace tc_graphblas
