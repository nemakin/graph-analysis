#include "utils.hpp"
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace graphblas_utils {
void print_matrix(GrB_Matrix A) {
  GrB_Index nrows, ncols;
  GrB_Matrix_nrows(&nrows, A);
  GrB_Matrix_ncols(&ncols, A);

  for (GrB_Index i = 0; i < nrows; ++i) {
    for (GrB_Index j = 0; j < ncols; ++j) {
      uint64_t val;
      GrB_Info info = GrB_Matrix_extractElement_UINT64(&val, A, i, j);
      if (info == GrB_SUCCESS) {
        std::cout << val << " ";
      } else if (info == GrB_NO_VALUE) {
        std::cout << "0 ";
      } else {
        throw std::runtime_error("Error extracting element from matrix");
      }
    }
    std::cout << std::endl;
  }
}
// Extract upper triangle of A into *U.
// strict = true → strict upper (j > i)
// strict = false → inclusive upper (j ≥ i)
GrB_Info extract_upper(GrB_Matrix *U, GrB_Matrix A, bool strict) {
  GrB_Index nrows, ncols;

  GrB_Matrix_nrows(&nrows, A);

  GrB_Matrix_ncols(&ncols, A);

  // Allocate U with same type as A
  GrB_Type type;
  GxB_Matrix_type(&type, A);

  GrB_Matrix_new(U, type, nrows, ncols);

  // k = 1 → strict upper; k = 0 → include diagonal
  int64_t k = strict ? 1 : 0;

  return GrB_Matrix_select_INT64(*U,
                                 NULL,     // no mask
                                 NULL,     // no accumulator
                                 GrB_TRIU, // upper-triangular selector
                                 A,
                                 k,   // diagonal offset
                                 NULL // descriptor
  );
}

GrB_Matrix load_graph(const std::string &filepath, bool triangular) {
  std::ifstream infile(filepath);
  if (!infile.is_open())
    throw std::runtime_error("Cannot open file: " + filepath);
  GrB_Index nrows = 0, ncols = 0, n_lines = 0;
  std::string header;
  if (!std::getline(infile, header))
    throw std::runtime_error("File is empty: " + filepath);
  std::istringstream hss(header);
  if (!(hss >> nrows >> ncols >> n_lines))
    throw std::runtime_error("Header format error in file: " + filepath);
  std::vector<std::pair<GrB_Index, GrB_Index>> edges;
  bool zero_based = false;
  std::string line;
  for (GrB_Index i = 0; i < n_lines; ++i) {
    if (!std::getline(infile, line))
      throw std::runtime_error("Not enough edge lines in file: " + filepath);
    std::istringstream iss(line);
    GrB_Index u, v;
    if (!(iss >> u >> v))
      continue;
    zero_based |= (u == 0 || v == 0);
    edges.emplace_back(u, v);
  }
  infile.close();
  GrB_Matrix A;
  GrB_Matrix_new(&A, GrB_UINT64, nrows, ncols);
  for (const auto &[raw_u, raw_v] : edges) {
    GrB_Index u = raw_u - (zero_based ? 0 : 1);
    GrB_Index v = raw_v - (zero_based ? 0 : 1);
    if (u >= nrows || v >= ncols)
      throw std::runtime_error("Edge index outside matrix dimensions: " + filepath);
    GrB_Matrix_setElement_UINT64(A, 1, u, v);
    GrB_Matrix_setElement_UINT64(A, 1, v, u); // undirected
  }
  if (!triangular) {
    return A;
  } else {
    GrB_Matrix U;
    extract_upper(&U, A, true);
    GrB_Matrix_free(&A);
    return U;
  }
}

} // namespace graphblas_utils
