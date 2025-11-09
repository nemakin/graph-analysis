#include <fstream>
#include <set>
#include <spla.hpp>
#include <vector>
#include <chrono>
#include <filesystem>
#include <sstream>
#include <stdexcept>
#include <iostream>

static int n = 0;
static int edges_count = 0;
static unsigned int weight = 0;
static int el_cnt = 0;

static std::vector<unsigned int> buffer1;
static std::vector<unsigned int> buffer2;
static spla::ref_ptr<spla::Matrix> a;
static spla::ref_ptr<spla::Vector> mst;

static const unsigned int INF = 1e9;
static spla::ref_ptr<spla::Scalar> zero_uint = spla::Scalar::make_uint(0);
static spla::ref_ptr<spla::Scalar> inf_uint = spla::Scalar::make_uint(INF);

struct Edge
{
    int u;
    int v;
    int w;
};

void buildMatrixFromDIMACS_lim(const std::string &filename, int max_edges)
{
    std::ifstream fin(filename);
    if (!fin.is_open())
    {
        throw std::runtime_error("Не удалось открыть файл: " + filename);
    }

    std::string line;
    int numNodes = 0, numEdges = 0;

    int max_node_id = 0;

    std::vector<Edge> edges;

    while (std::getline(fin, line) && el_cnt < max_edges)
    {
        if (line.empty() || line[0] == 'c')
            continue;

        std::istringstream iss(line);
        char type;
        iss >> type;

        if (type == 'p')
        {
            std::string sp;
            iss >> sp >> numNodes >> numEdges;
            std::cout << "Original graph: " << numNodes << " nodes, " << numEdges << " edges\n";
        }
        else if (type == 'a')
        {
            int u, v;
            int w;
            iss >> u >> v >> w;

            u -= 1;
            v -= 1;

            max_node_id = std::max(max_node_id, std::max(u, v));
            edges.push_back({u, v, w});
            el_cnt += 1;
        }
    }

    fin.close();
    n = max_edges;
    edges_count = max_edges;

    a = spla::Matrix::make(max_edges, max_edges, spla::UINT);

    for (const auto &e : edges)
    {
        a->set_uint(e.u, e.v, e.w);
    }
    std::cout << "Loaded graph: " << max_node_id << " nodes, " << el_cnt << " edges\n";
}

void load_graph_mm(const std::string &path, int n_loc)
{
    std::ifstream fin(path);
    if (!fin.is_open())
    {
        throw std::runtime_error("Cannot open file: " + path);
    }

    std::vector<Edge> edges;
    std::string line;

    while (std::getline(fin, line))
    {
        if (line.empty() || line[0] == '%')
            continue;

        std::istringstream iss(line);
        Edge e{};
        if (!(iss >> e.u >> e.v >> e.w))
            continue;
        if (e.u > 0 && e.v > 0)
            e.u--, e.v--;

        edges.push_back(e);
    }

    fin.close();
    n = n_loc;
    edges_count = n_loc;
    std::sort(edges.begin(), edges.end(), [](const Edge &a, const Edge &b)
              {
        if (a.u == b.u)
            return a.v < b.v;
        return a.u < b.u; });

    a = spla::Matrix::make(n_loc, n_loc, spla::UINT);
    for (const auto &e : edges)
    {
        if (e.u < n_loc && e.v < n_loc)
        {
            a->set_uint(e.u, e.v, e.w);
            el_cnt += 1;
        }
    }
    std::cout << "loaded elements: " << el_cnt << "\n";
}

void update(std::set<std::pair<unsigned int, unsigned int>> &s, const spla::ref_ptr<spla::Vector> &v)
{
    auto sz = spla::Scalar::make_uint(0);
    spla::exec_v_count_mf(sz, v);
    auto keys_view = spla::MemView::make(buffer1.data(), sz->as_uint(), false);
    auto values_view = spla::MemView::make(buffer2.data(), sz->as_uint(), false);
    v->read(keys_view, values_view);
    auto keys = (unsigned int *)keys_view->get_buffer();
    auto values = (unsigned int *)values_view->get_buffer();
    for (unsigned int i = 0; i < sz->as_uint(); i++)
    {
        s.insert({values[i], keys[i]});
    }
}

using clock_ = std::chrono::steady_clock;

void compute_internal()
{
    mst = spla::Vector::make(n, spla::UINT);

    auto d = spla::Vector::make(n, spla::UINT);
    auto changed = spla::Vector::make(n, spla::UINT);
    auto v_row = spla::Vector::make(n, spla::UINT);
    auto min_v = spla::Scalar::make_uint(INF);

    changed->set_fill_value(zero_uint);
    mst->set_fill_value(inf_uint);
    d->set_fill_value(inf_uint);
    v_row->set_fill_value(inf_uint);

    weight = 0;
    if (n <= 1 || edges_count == 0)
    {
        return;
    }

    std::set<std::pair<unsigned int, unsigned int>> s;
    std::vector<bool> visited(n, false);

    for (int i = 0; i < n; i++)
    {
        if (!visited[i])
        {
            unsigned int v = i;
            d->set_uint(v, 0);
            visited[v] = true;
            spla::exec_m_extract_row(v_row, a, v, spla::IDENTITY_UINT);
            spla::exec_v_eadd_fdb(d, v_row, changed, spla::MIN_UINT);
            spla::exec_v_assign_masked(mst, changed, spla::Scalar::make_uint(v), spla::SECOND_UINT,
                                       spla::NQZERO_UINT);

            update(s, changed);
            while (!s.empty())
            {
                unsigned int w = s.begin()->first;
                v = s.begin()->second;
                s.erase(s.begin());
                if (visited[v])
                    continue;

                weight += w;
                d->set_uint(v, 0);
                visited[v] = true;
                spla::exec_m_extract_row(v_row, a, v, spla::IDENTITY_UINT);
                spla::exec_v_eadd_fdb(d, v_row, changed, spla::MIN_UINT);
                spla::exec_v_assign_masked(mst, changed, spla::Scalar::make_uint(v),
                                           spla::SECOND_UINT,
                                           spla::NQZERO_UINT);

                update(s, changed);
            }
        }
    }
}

std::chrono::seconds compute()
{
    auto start = clock_::now();
    compute_internal();
    auto end = clock_::now();
    return std::chrono::duration_cast<std::chrono::seconds>(end - start);
}

int main()
{
    try
    {
        std::filesystem::path graph_path = "/Users/nikitalukonenko/Studying/third_course/experiment/sources/gbtl/datasets/internet_2/internet.mtx";
        load_graph_mm(graph_path, 50000);
        auto execution_time = compute();
        std::cout << "Algorithm execution time: "
                  << execution_time.count()
                  << " seconds\n";

        std::cout << "MST weight: " << weight << "\n";
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}