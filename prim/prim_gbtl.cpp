
#include <fstream>
#include <sstream>
#include <chrono>

#include <vector>
#include <string>
#include <stdexcept>
#include <graphblas/graphblas.hpp>

#include <iostream>
#include "demo/Timer.hpp"
#include <algorithms/mst.hpp>
#include <graphblas/graphblas.hpp>

using namespace grb;
using namespace algorithms;

#define BOOST_TEST_MAIN
#define BOOST_TEST_MODULE mst_test_suite

//============================================================
// Функция чтения DIMACS и построения матрицы
//============================================================

using namespace grb;

std::tuple<
    Matrix<double>,
    grb::Vector<grb::IndexType>>
load_matrix_mm_lim(const std::string &path, grb::IndexType max_nodes = 0)
{
    std::ifstream file(path);
    if (!file.is_open())
    {
        throw std::runtime_error("Не удалось открыть файл: " + path);
    }

    std::string line;

    // 1️⃣ Пропускаем комментарии и заголовок
    bool header_found = false;
    while (std::getline(file, line))
    {
        if (line.rfind("%%MatrixMarket", 0) == 0)
        {
            header_found = true;
            continue;
        }
        if (line.empty() || line[0] == '%')
        {
            continue;
        }
        else
        {
            break;
        }
    }

    if (!header_found)
    {
        throw std::runtime_error("Неверный формат файла MatrixMarket: отсутствует заголовок %%MatrixMarket");
    }

    // 2️⃣ Читаем размеры и количество ненулевых элементов
    std::istringstream dims(line);
    grb::IndexType nrows, ncols, nnz;
    dims >> nrows >> ncols >> nnz;

    if (nrows == 0 || ncols == 0)
    {
        throw std::runtime_error("Размеры матрицы не могут быть нулевыми");
    }

    // 3️⃣ Читаем все рёбра в вектор
    struct Edge
    {
        grb::IndexType i;
        grb::IndexType j;
        double value;
    };

    std::vector<Edge> edges;
    edges.reserve(nnz);

    grb::IndexType i, j;
    double value;

    while (file >> i >> j >> value)
    {
        // MatrixMarket использует 1-based индексацию
        if (i == 0 || j == 0)
            continue;

        i -= 1;
        j -= 1;

        edges.push_back({i, j, value});
    }

    file.close();

    std::cout << "Прочитано рёбер: " << edges.size() << std::endl;

    // 4️⃣ Сортируем рёбра: сначала по i, потом по j
    std::sort(edges.begin(), edges.end(),
              [](const Edge &a, const Edge &b)
              {
                  if (a.i != b.i)
                      return a.i < b.i;
                  return a.j < b.j;
              });

    std::cout << "Рёбра отсортированы" << std::endl;

    // 5️⃣ Определяем размер матрицы
    grb::IndexType matrix_size = (max_nodes > 0) ? std::min(max_nodes, std::max(nrows, ncols)) : std::max(nrows, ncols);

    grb::Matrix<double> matrix(matrix_size, matrix_size);

    // 6️⃣ Заполняем матрицу отсортированными рёбрами
    grb::IndexType counter = 0;
    for (const auto &edge : edges)
    {
        // Пропускаем рёбра, выходящие за пределы заданного размера
        if (edge.i >= matrix_size || edge.j >= matrix_size)
            continue;

        try
        {
            matrix.setElement(edge.i, edge.j, edge.value);
            ++counter;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Ошибка при добавлении элемента (" << edge.i << "," << edge.j
                      << "): " << e.what() << std::endl;
        }
    }

    std::cout << "Загружено элементов: " << counter
              << " в матрицу размера " << matrix_size << "x" << matrix_size << std::endl;

    grb::Vector<grb::IndexType> parents(matrix_size);

    return std::make_tuple(matrix, parents);
}

std::tuple<
    Matrix<double>,
    grb::Vector<grb::IndexType>>
load_matrix_mm(const std::string &path)
{
    std::ifstream file(path);
    if (!file.is_open())
    {
        throw std::runtime_error("Не удалось открыть файл: " + path);
    }

    std::string line;

    // 1️⃣ Пропускаем комментарии и заголовок
    bool header_found = false;
    while (std::getline(file, line))
    {
        if (line.rfind("%%MatrixMarket", 0) == 0)
        {
            header_found = true;
            continue;
        }
        if (line.empty() || line[0] == '%')
        {
            continue; // комментарии
        }
        else
        {
            // первая строка без % — размеры матрицы
            break;
        }
    }

    if (!header_found)
    {
        throw std::runtime_error("Неверный формат файла MatrixMarket: отсутствует заголовок %%MatrixMarket");
    }

    // 2️⃣ Читаем размеры и количество ненулевых элементов
    std::istringstream dims(line);
    grb::IndexType nrows, ncols, nnz;
    dims >> nrows >> ncols >> nnz;

    if (nrows == 0 || ncols == 0)
    {
        throw std::runtime_error("Размеры матрицы не могут быть нулевыми");
    }

    grb::Matrix<double> matrix(nrows, ncols);

    // 3️⃣ Читаем ненулевые элементы (в формате i j value)
    grb::IndexType i, j;
    double value;

    grb::IndexType counter = 0;
    while (file >> i >> j >> value)
    {
        // MatrixMarket использует 1-based индексацию
        if (i == 0 || j == 0)
            continue;
        i -= 1;
        j -= 1;
        try
        {
            matrix.setElement(i, j, value);
            ++counter;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Ошибка при добавлении элемента (" << i << "," << j
                      << "): " << e.what() << std::endl;
        }
    }

    std::cout << "Загружено элементов: " << counter << " из " << nnz << std::endl;
    grb::Vector<grb::IndexType> parents(ncols);

    return std::make_tuple(matrix, parents);

    // return matrix;
}

std::tuple<
    Matrix<double>,
    grb::Vector<grb::IndexType>>
trim_matrix(const grb::Matrix<double> &matrix, grb::IndexType new_size)
{
    using IndexType = grb::IndexType;

    const IndexType nrows = matrix.nrows();
    const IndexType ncols = matrix.ncols();
    IndexType max_val = 0;

    if (new_size > nrows || new_size > ncols)
    {
        throw std::invalid_argument("Новый размер превышает исходный размер матрицы");
    }

    Matrix<double> trimmed(new_size, new_size);
    grb::Vector<grb::IndexType> parents(new_size);

    // Копируем только элементы в пределах new_size
    for (IndexType i = 0; i < new_size; ++i)
    {
        for (IndexType j = 0; j < new_size; ++j)
        {
            double value;
            try
            {
                value = matrix.extractElement(i, j);
                if (value > max_val)
                {
                    max_val = value;
                }
                trimmed.setElement(i, j, value);
                // }
            }

            catch (const grb::NoValueException &)
            {
                continue;
            }
        }
    }
    std::cout << "Graph: " << max_val << " nodes, " << new_size << " edges\n";
    std::cout << "matrix trimmed" << "\n";
    return std::make_tuple(trimmed, parents);
}

std::tuple<
    Matrix<double>,
    grb::Vector<grb::IndexType>>
buildMatrixFromDIMACS(const std::string &filename)
{
    std::ifstream fin(filename);
    if (!fin.is_open())
    {
        throw std::runtime_error("Не удалось открыть файл: " + filename);
    }

    std::string line;
    size_t numNodes = 0, numEdges = 0;

    std::vector<IndexType> i_indices;
    std::vector<IndexType> j_indices;
    std::vector<double> weights;

    while (std::getline(fin, line))
    {
        if (line.empty() || line[0] == 'c')
            continue; // пропускаем комментарии

        std::istringstream iss(line);
        char type;
        iss >> type;

        if (type == 'p')
        {
            std::string sp;
            iss >> sp >> numNodes >> numEdges;
            std::cout << "Graph: " << numNodes << " nodes, " << numEdges << " edges\n";
        }
        else if (type == 'a')
        {
            size_t u, v;
            double w;
            iss >> u >> v >> w;

            i_indices.push_back(u - 1);
            j_indices.push_back(v - 1);
            weights.push_back(w);
        }
    }

    fin.close();

    Matrix<double> M(numNodes, numNodes);
    M.build(i_indices, j_indices, weights);

    grb::Vector<grb::IndexType> parents(numNodes);

    return std::make_tuple(M, parents);
}

std::tuple<Matrix<double>, grb::Vector<grb::IndexType>>
buildMatrixFromDIMACS_lim(const std::string &filename, size_t max_edges)
{
    std::ifstream fin(filename);
    if (!fin.is_open())
    {
        throw std::runtime_error("Не удалось открыть файл: " + filename);
    }

    std::string line;
    size_t numNodes = 0, numEdges = 0;
    size_t edge_counter = 0;
    size_t max_node_id = 0;

    std::vector<IndexType> i_indices;
    std::vector<IndexType> j_indices;
    std::vector<double> weights;

    while (std::getline(fin, line) && edge_counter < max_edges)
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

            i_indices.reserve(max_edges);
            j_indices.reserve(max_edges);
            weights.reserve(max_edges);
        }
        else if (type == 'a' && edge_counter < max_edges)
        {
            size_t u, v;
            double w;
            iss >> u >> v >> w;

            u -= 1;
            v -= 1;

            max_node_id = std::max(max_node_id, std::max(u, v));

            i_indices.push_back(u);
            j_indices.push_back(v);
            weights.push_back(w);

            edge_counter++;
        }
    }

    fin.close();

    // Используем max_node_id + 1 как размер матрицы
    size_t actual_size = max_node_id + 1;
    std::cout << "Loaded graph: " << actual_size << " nodes, " << edge_counter << " edges\n";

    // Создаём матрицу смежности
    Matrix<double> M(actual_size, actual_size);
    M.build(i_indices, j_indices, weights);

    // Создаём вектор родителей для MST
    grb::Vector<grb::IndexType> parents(actual_size);

    return std::make_tuple(M, parents);
}

using clock_ = std::chrono::steady_clock;

int main()
{
    // auto [graph2, parents2] = load_matrix_mm_lim("../datasets/internet_2/internet.mtx", 50000);
    auto [graph2, parents2] = buildMatrixFromDIMACS_lim("../datasets/USA_road_d_COL.gr", 50000);

    auto start = clock_::now();
    auto result = mst(graph2, parents2);
    auto end = clock_::now();
    std::chrono::seconds execution_time = std::chrono::duration_cast<std::chrono::seconds>(end - start);

    std::cout << "time: " << execution_time.count() << " seconds";

    std::cout << result << "\n";
    return 0;
}
