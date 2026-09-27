#pragma once

#include <spla.hpp>
#include <string>
#include <vector>

namespace tc_spla {
using namespace spla;
std::vector<double> benchmark(const char *filename, bool triangular,
                              int num_runs, bool accelerated);

void burkhardt(int &ntrins, const ref_ptr<Matrix> &A, const ref_ptr<Matrix> &B);
void sandia(int &ntrins, const ref_ptr<Matrix> &A, const ref_ptr<Matrix> &B);
} // namespace tc_spla
