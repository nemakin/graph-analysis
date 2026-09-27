#pragma once

#include <GraphBLAS.h>
#include <string>
#include <vector>

namespace tc_graphblas {
uint64_t triangles_counting(GrB_Matrix A, GrB_Matrix workspace,
                            bool triangular);
uint64_t burkhardt(GrB_Matrix A, GrB_Matrix workspace);
uint64_t sandia(GrB_Matrix A, GrB_Matrix workspace);
std::vector<double> benchmark(const char *filename, bool triangular,
                              const int num_iters);
} // namespace tc_graphblas
