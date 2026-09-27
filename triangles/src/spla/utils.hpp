#pragma once
#include <spla.hpp>
#include <string>

namespace spla_utils
{
    spla::ref_ptr<spla::Matrix> load_graph(const std::string &path, bool triangular);

    void print_matrix(const spla::ref_ptr<spla::Matrix> &matrix);
}
