#pragma once

#include <GraphBLAS.h>
#include <vector>
#include <string>

namespace graphblas_utils
{
    GrB_Info extract_upper(GrB_Matrix *U, GrB_Matrix A, bool strict);

    GrB_Matrix load_graph(const std::string &filepath, bool triangular);

    void print_matrix(GrB_Matrix A);
}
