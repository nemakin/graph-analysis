#include "utils.hpp"
#include <algorithm>
#include <fstream>
#include <stdexcept>

namespace spla_utils
{
    spla::ref_ptr<spla::Matrix> load_graph(const std::string &path, bool triangular)
    {
        std::ifstream input(path);
        std::size_t nrows, ncols, edge_count;
        if (!(input >> nrows >> ncols >> edge_count)) {
            throw std::runtime_error("Invalid graph header: " + path);
        }

        bool zero_based = false;
        for (std::size_t k = 0; k < edge_count; ++k) {
            std::size_t u, v;
            if (!(input >> u >> v)) {
                throw std::runtime_error("Invalid graph edge: " + path);
            }
            zero_based |= (u == 0 || v == 0);
        }

        if (zero_based) {
            // SPLA's MtxLoader assumes one-based indices while computing stats.
            // Read zero-based edge lists directly to avoid indexing before row zero.
            input.clear();
            input.seekg(0);
            input >> nrows >> ncols >> edge_count;
            auto A = spla::Matrix::make(nrows, ncols, spla::INT);
            for (std::size_t k = 0; k < edge_count; ++k) {
                std::size_t u, v;
                input >> u >> v;
                if (u >= nrows || v >= ncols) {
                    throw std::runtime_error("Edge index outside matrix dimensions: " + path);
                }
                if (u == v) continue;
                if (triangular) {
                    A->set_int(std::max(u, v), std::min(u, v), 1);
                } else {
                    A->set_int(u, v, 1);
                    A->set_int(v, u, 1);
                }
            }
            return A;
        }

        spla::MtxLoader loader;
        if (!loader.load(path, true, true, true))
        {
            throw std::runtime_error("Failed to load graph: " + path);
        }

        const auto n = loader.get_n_rows();
        auto A = spla::Matrix::make(n, n, spla::INT);
        const auto &Ai = loader.get_Ai();
        const auto &Aj = loader.get_Aj();

        for (std::size_t k = 0; k < loader.get_n_values(); ++k)
        {
            if (!triangular || Ai[k] > Aj[k])
            {
                A->set_int(Ai[k], Aj[k], 1);
            }
        }

        return A;
    }

    void print_matrix(const spla::ref_ptr<spla::Matrix> &matrix)
    {
        std::cout << std::endl;
        for (uint i = 0; i < matrix->get_n_rows(); ++i)
        {
            for (uint j = 0; j < matrix->get_n_cols(); ++j)
            {
                int value;
                matrix->get_int(i, j, value);
                std::cout << value << " ";
            }
            std::cout << std::endl;
        }
    }

}
