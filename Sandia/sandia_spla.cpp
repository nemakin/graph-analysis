#include <spla.hpp>
#include <algorithm.hpp>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <chrono>
using namespace spla;

static spla::ref_ptr<spla::Matrix> a;
static spla::ref_ptr<spla::Matrix> b;

static int el_cnt = 0;
static int edges_count = 0;
// static int n = 0;
int max_node_id = 0;

struct Edge
{
    int u;
    int v;
    int w;
};

void load_graph_mm(const std::string &path)
{
    std::ifstream fin(path);
    if (!fin.is_open())
    {
        throw std::runtime_error("Cannot open file: " + path);
    }

    std::vector<Edge> edges;
    std::string line;
    int rows, cols, nnz;
    int cf = 0;

    while (std::getline(fin, line))
    {
        if (line.empty() || line[0] == '%')
            continue;
        cf += 1;
        std::istringstream iss(line);
        Edge e{};
        Edge e1{};
        if (!(iss >> e.u >> e.v))
            continue;
        if (cf == 1)
        {
            continue;
        };
        max_node_id = std::max(max_node_id, std::max((e.u), (e.v)));
        e.u--, e.v--;
        e.w = 1;
        edges.push_back(e);
    }

    fin.close();
    a = spla::Matrix::make(max_node_id, max_node_id, spla::INT);
    std::cout << "mx_el: " << max_node_id << "\n";
    for (const auto &e : edges)
    {
        a->set_uint(e.v, e.u, e.w); 
        el_cnt += 1;
    }
    std::cout << "loaded elements: " << el_cnt << "\n";
    a->set_format(spla::FormatMatrix::AccCsr);
}
using clock_ = std::chrono::steady_clock;

int main()
{
    try
    {
        std::filesystem::path graph_path = "/Users/nikitalukonenko/Studying/third_course/experiment/sources/gbtl/datasets/graph500-scale18-ef16_adj.mmio";
        load_graph_mm(graph_path);
        b = spla::Matrix::make(max_node_id, max_node_id, spla::INT);
        b->set_format(spla::FormatMatrix::AccCsr);
        int32_t ntrins;
        auto res = tc(ntrins, a, b);
        std::cout << "res: " << static_cast<int>(res) << "\n";
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
