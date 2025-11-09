#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <GraphBLAS.h>

#define CHECK(x)                                                                                      \
    do                                                                                                \
    {                                                                                                 \
        GrB_Info _info = (x);                                                                         \
        if (_info != GrB_SUCCESS)                                                                     \
        {                                                                                             \
            fprintf(stderr, "GraphBLAS error: %s returned %d (line %d)\n", #x, (int)_info, __LINE__); \
            exit(1);                                                                                  \
        }                                                                                             \
    } while (0)

typedef struct
{
    GrB_Index i;
    GrB_Index j;
    double value;
} MMEdge;

static int mmedge_cmp(const void *a, const void *b)
{
    const MMEdge *A = (const MMEdge *)a;
    const MMEdge *B = (const MMEdge *)b;
    if (A->i < B->i)
        return -1;
    if (A->i > B->i)
        return 1;
    if (A->j < B->j)
        return -1;
    if (A->j > B->j)
        return 1;
    return 0;
}

GrB_Info load_dimacs_lim(const char *path, GrB_Matrix *out_matrix, GrB_Vector *out_parents, GrB_Index max_edges)
{
    FILE *f = fopen(path, "r");
    if (f == NULL)
    {
        fprintf(stderr, "load_dimacs_lim: cannot open file %s\n", path);
        return GrB_INVALID_VALUE;
    }

    char line[1024];
    GrB_Index num_nodes = 0, num_edges = 0;
    GrB_Index edge_counter = 0;
    GrB_Index max_node_id = 0;
    GrB_Index capacity = 16;
    MMEdge *edges = malloc(capacity * sizeof(MMEdge));
    if (edges == NULL)
    {
        fclose(f);
        return GrB_OUT_OF_MEMORY;
    }

    while (fgets(line, sizeof(line), f) != NULL && edge_counter < max_edges)
    {
        if (line[0] == '\n' || line[0] == 'c')
            continue;

        if (line[0] == 'p')
        {
            char sp[32];
            if (sscanf(line, "p %s %llu %llu", sp, &num_nodes, &num_edges) != 3)
            {
                free(edges);
                fclose(f);
                return GrB_INVALID_VALUE;
            }
            fprintf(stderr, "Original graph: %llu nodes, %llu edges\n",
                    (unsigned long long)num_nodes, (unsigned long long)num_edges);
            continue;
        }

        if (line[0] == 'a' && edge_counter < max_edges)
        {
            GrB_Index u, v;
            double w;
            if (sscanf(line, "a %llu %llu %lf", &u, &v, &w) != 3)
            {
                free(edges);
                fclose(f);
                return GrB_INVALID_VALUE;
            }

            u--;
            v--;

            max_node_id = (u > max_node_id) ? u : max_node_id;
            max_node_id = (v > max_node_id) ? v : max_node_id;

            if (edge_counter >= capacity)
            {
                capacity *= 2;
                MMEdge *tmp = realloc(edges, capacity * sizeof(MMEdge));
                if (tmp == NULL)
                {
                    free(edges);
                    fclose(f);
                    return GrB_OUT_OF_MEMORY;
                }
                edges = tmp;
            }

            edges[edge_counter].i = u;
            edges[edge_counter].j = v;
            edges[edge_counter].value = w;
            edge_counter++;
        }
    }

    fclose(f);

    GrB_Index matrix_size = max_node_id + 1;
    fprintf(stderr, "Loaded graph: %llu nodes, %llu edges\n",
            (unsigned long long)matrix_size, (unsigned long long)edge_counter);

    /* Сортируем рёбра */
    if (edge_counter > 1)
    {
        qsort(edges, edge_counter, sizeof(MMEdge), mmedge_cmp);
    }

    GrB_Matrix A = NULL;
    GrB_Info info = GrB_Matrix_new(&A, GrB_FP64, matrix_size, matrix_size);
    if (info != GrB_SUCCESS)
    {
        free(edges);
        return info;
    }

    for (GrB_Index k = 0; k < edge_counter; k++)
    {
        info = GrB_Matrix_setElement_FP64(A, edges[k].value, edges[k].i, edges[k].j);
        if (info != GrB_SUCCESS)
        {
            GrB_Matrix_free(&A);
            free(edges);
            return info;
        }
    }

    GrB_Vector parents = NULL;
    info = GrB_Vector_new(&parents, GrB_UINT64, matrix_size);
    if (info != GrB_SUCCESS)
    {
        GrB_Matrix_free(&A);
        free(edges);
        return info;
    }

    free(edges);

    *out_matrix = A;
    *out_parents = parents;

    return GrB_SUCCESS;
}

GrB_Info load_matrix_mm_lim(const char *path, GrB_Matrix *out_matrix, GrB_Vector *out_parents, GrB_Index max_nodes)
{
    FILE *f = fopen(path, "r");
    if (f == NULL)
    {
        fprintf(stderr, "load_matrix_mm_lim: cannot open file %s\n", path);
        return GrB_INVALID_VALUE;
    }

    char line[1024];
    bool header_found = false;

    while (fgets(line, sizeof(line), f) != NULL)
    {
        if (strncmp(line, "%%MatrixMarket", 14) == 0)
        {
            header_found = true;
            continue;
        }
        if (line[0] == '%' || line[0] == '\n' || line[0] == '\r')
        {
            continue;
        }
        break;
    }

    if (!header_found)
    {
        fprintf(stderr, "load_matrix_mm_lim: not a MatrixMarket file (missing header)\n");
        fclose(f);
        return GrB_INVALID_VALUE;
    }

    GrB_Index nrows = 0, ncols = 0;
    unsigned long long nnz_ull = 0;
    if (sscanf(line, "%llu %llu %llu", &nrows, &ncols, &nnz_ull) != 3)
    {
        fprintf(stderr, "load_matrix_mm_lim: failed to parse dimensions line: %s\n", line);
        fclose(f);
        return GrB_INVALID_VALUE;
    }

    GrB_Index nnz = (GrB_Index)nnz_ull;
    MMEdge *edges = NULL;
    GrB_Index capacity = (nnz > 0) ? nnz : 16;
    edges = malloc(capacity * sizeof(MMEdge));
    if (edges == NULL)
    {
        fprintf(stderr, "load_matrix_mm_lim: malloc failed\n");
        fclose(f);
        return GrB_OUT_OF_MEMORY;
    }

    unsigned long long ii = 0, jj = 0;
    double val = 0.0;
    GrB_Index count = 0;
    while (fscanf(f, "%llu %llu %lf", &ii, &jj, &val) == 3)
    {
        if (ii == 0 || jj == 0)
        {
            continue;
        }

        GrB_Index ri = (GrB_Index)(ii - 1);
        GrB_Index cj = (GrB_Index)(jj - 1);

        if (count >= capacity)
        {
            GrB_Index newcap = capacity * 2;
            MMEdge *tmp = realloc(edges, newcap * sizeof(MMEdge));
            if (tmp == NULL)
            {
                free(edges);
                fclose(f);
                fprintf(stderr, "load_matrix_mm_lim: realloc failed\n");
                return GrB_OUT_OF_MEMORY;
            }
            edges = tmp;
            capacity = newcap;
        }

        edges[count].i = ri;
        edges[count].j = cj;
        edges[count].value = val;
        count++;
    }

    fclose(f);

    if (count > 1)
    {
        qsort(edges, count, sizeof(MMEdge), mmedge_cmp);
    }

    GrB_Index maxdim = (nrows > ncols) ? nrows : ncols;
    GrB_Index matrix_size = (max_nodes > 0) ? ((max_nodes < maxdim) ? max_nodes : maxdim) : maxdim;

    GrB_Matrix A = NULL;
    GrB_Info info = GrB_Matrix_new(&A, GrB_FP64, matrix_size, matrix_size);
    if (info != GrB_SUCCESS)
    {
        free(edges);
        fprintf(stderr, "load_matrix_mm_lim: GrB_Matrix_new failed (%d)\n", (int)info);
        return info;
    }

    GrB_Index loaded = 0;
    for (GrB_Index k = 0; k < count; k++)
    {
        if (edges[k].i >= matrix_size || edges[k].j >= matrix_size)
            continue;
        info = GrB_Matrix_setElement_FP64(A, edges[k].value, edges[k].i, edges[k].j);
        if (info != GrB_SUCCESS)
        {
            GrB_Matrix_free(&A);
            free(edges);
            fprintf(stderr, "load_matrix_mm_lim: GrB_Matrix_setElement_FP64 failed (%d)\n", (int)info);
            return info;
        }
        loaded++;
    }

    GrB_Vector parents = NULL;
    info = GrB_Vector_new(&parents, GrB_UINT64, matrix_size);
    if (info != GrB_SUCCESS)
    {
        GrB_Matrix_free(&A);
        free(edges);
        fprintf(stderr, "load_matrix_mm_lim: GrB_Vector_new failed (%d)\n", (int)info);
        return info;
    }

    free(edges);

    *out_matrix = A;
    *out_parents = parents;

    fprintf(stderr, "load_matrix_mm_lim: read edges=%llu, loaded=%llu, matrix_size=%llu\n",
            (unsigned long long)count, (unsigned long long)loaded, (unsigned long long)matrix_size);

    return GrB_SUCCESS;
}

typedef struct
{
    GrB_Index index;
    double weight;
} MSTType;

void mst_first(void *z, const void *x, const void *y)
{
    MSTType *result = (MSTType *)z;
    const MSTType *first = (const MSTType *)x;
    *result = *first;
}

void mst_weight(void *z, const void *x)
{
    double *result = (double *)z;
    const MSTType *val = (const MSTType *)x;
    *result = val->weight;
}

void mst_min(void *z, const void *x, const void *y)
{
    MSTType *result = (MSTType *)z;
    const MSTType *lhs = (const MSTType *)x;
    const MSTType *rhs = (const MSTType *)y;

    if (rhs->weight > lhs->weight)
    {
        *result = *lhs;
    }
    else
    {
        *result = *rhs;
    }
}
GrB_Index argmin_vector(GrB_Vector v)
{
    GrB_Index nvals;
    CHECK(GrB_Vector_nvals(&nvals, v));

    if (nvals == 0)
    {
        fprintf(stderr, "Empty vector in argmin\n");
        exit(1);
    }

    GrB_Index *indices = malloc(nvals * sizeof(GrB_Index));
    double *values = malloc(nvals * sizeof(double));
    GrB_Index n = nvals;
    CHECK(GrB_Vector_extractTuples_FP64(indices, values, &n, v));

    GrB_Index idx_min = indices[0];
    double val_min = values[0];

    for (GrB_Index i = 1; i < nvals; i++)
    {
        if (values[i] < val_min)
        {
            val_min = values[i];
            idx_min = indices[i];
        }
    }

    free(indices);
    free(values);
    return idx_min;
}

double mst_prim(GrB_Matrix graph, GrB_Vector mst_parents)
{
    GrB_Index rows, cols;
    CHECK(GrB_Matrix_nrows(&rows, graph));
    CHECK(GrB_Matrix_ncols(&cols, graph));

    if (rows != cols)
    {
        fprintf(stderr, "Matrix must be square\n");
        exit(1);
    }

    GrB_Type MSTType_type;
    CHECK(GrB_Type_new(&MSTType_type, sizeof(MSTType)));

    GrB_UnaryOp weight_op;
    CHECK(GrB_UnaryOp_new(&weight_op, mst_weight, GrB_FP64, MSTType_type));

    GrB_BinaryOp MST_min;
    CHECK(GrB_BinaryOp_new(&MST_min, mst_min, MSTType_type, MSTType_type, MSTType_type));

    GrB_Matrix A;
    CHECK(GrB_Matrix_new(&A, MSTType_type, rows, cols));

    GrB_Index nvals;
    CHECK(GrB_Matrix_nvals(&nvals, graph));

    GrB_Index *i = malloc(nvals * sizeof(GrB_Index));
    GrB_Index *j = malloc(nvals * sizeof(GrB_Index));
    double *vals = malloc(nvals * sizeof(double));

    if (!i || !j || !vals)
    {
        fprintf(stderr, "Failed to allocate memory\n");
        exit(1);
    }

    CHECK(GrB_Matrix_extractTuples_FP64(i, j, vals, &nvals, graph));

    MSTType *new_vals = malloc(nvals * sizeof(MSTType));
    if (!new_vals)
    {
        fprintf(stderr, "Failed to allocate memory for new_vals\n");
        exit(1);
    }
    for (GrB_Index ix = 0; ix < nvals; ix++)
    {
        new_vals[ix].index = i[ix];
        new_vals[ix].weight = vals[ix];
    }

    GrB_BinaryOp MST_first;
    CHECK(GrB_BinaryOp_new(&MST_first, mst_first, MSTType_type, MSTType_type, MSTType_type));
    CHECK(GrB_Matrix_build(A, i, j, (void *)new_vals, nvals, MST_first));

    double total_weight = 0.0;

    GrB_Vector mask;
    CHECK(GrB_Vector_new(&mask, GrB_BOOL, rows));

    CHECK(GrB_Vector_clear(mst_parents));

    GrB_Index start = 0;
    CHECK(GrB_Vector_setElement_BOOL(mask, true, start));

    GrB_Vector d;
    CHECK(GrB_Vector_new(&d, MSTType_type, rows));

    GrB_Vector weights;
    CHECK(GrB_Vector_new(&weights, GrB_FP64, rows));

    CHECK(GrB_Col_extract(d, NULL, NULL, A, GrB_ALL, rows, start, NULL));

    GrB_Index visited_count = 1;
    while (visited_count < rows)
    {
        CHECK(GrB_Vector_clear(weights));
        CHECK(GrB_Vector_apply(weights, mask, NULL, weight_op, d, GrB_DESC_RC));

        GrB_Index nvals_weights;
        CHECK(GrB_Vector_nvals(&nvals_weights, weights));
        if (nvals_weights == 0)
            break;

        GrB_Index u = argmin_vector(weights);

        MSTType edge_info;
        CHECK(GrB_Vector_extractElement_UDT(&edge_info, d, u));

        // Добавляем вес в общую сумму
        total_weight += edge_info.weight;

        CHECK(GrB_Vector_setElement_UINT64(mst_parents, edge_info.index, u));

        CHECK(GrB_Vector_setElement_BOOL(mask, true, u));
        visited_count++;

        GrB_Vector new_edges;
        CHECK(GrB_Vector_new(&new_edges, MSTType_type, rows));
        CHECK(GrB_Col_extract(new_edges, NULL, NULL, A, GrB_ALL, rows, u, NULL));

        CHECK(GrB_eWiseAdd(d, mask, NULL, MST_min, d, new_edges, GrB_DESC_R));

        CHECK(GrB_Vector_free(&new_edges));
    }

    free(i);
    free(j);
    free(vals);
    free(new_vals);

    CHECK(GrB_Matrix_free(&A));
    CHECK(GrB_Vector_free(&mask));
    CHECK(GrB_Vector_free(&d));
    CHECK(GrB_Vector_free(&weights));
    CHECK(GrB_Type_free(&MSTType_type));
    CHECK(GrB_BinaryOp_free(&MST_first));
    CHECK(GrB_BinaryOp_free(&MST_min));
    CHECK(GrB_UnaryOp_free(&weight_op));

    return total_weight;
}

int main(int argc, char **argv)
{
    CHECK(GrB_init(GrB_NONBLOCKING));
    printf("=== Пример алгоритма Прима с использованием GraphBLAS ===\n");
    const GrB_Index n = 3;
    GrB_Matrix graph = NULL;
    GrB_Vector mst_parents = NULL;

    if (argc < 2)
    {
        CHECK(GrB_Matrix_new(&graph, GrB_FP64, n, n));
        CHECK(GrB_Matrix_setElement_FP64(graph, 2.0, 0, 1));
        CHECK(GrB_Matrix_setElement_FP64(graph, 2.0, 1, 0));

        CHECK(GrB_Matrix_setElement_FP64(graph, 3.0, 1, 2));
        CHECK(GrB_Matrix_setElement_FP64(graph, 3.0, 2, 1));
        CHECK(GrB_Vector_new(&mst_parents, GrB_UINT64, n));
    }
    else if (argc > 2 && strcmp(argv[1], "--dimacs") == 0)
    {
        CHECK(load_dimacs_lim(argv[2], &graph, &mst_parents, 50000));
    }
    else
    {
        CHECK(load_matrix_mm_lim(argv[1], &graph, &mst_parents, 50000));
    }

    printf("\n--- Запуск алгоритма Прима ---\n");
    double total_weight = mst_prim(graph, mst_parents);

    printf("\n=== Результаты ===\n");
    printf("Общий вес минимального остовного дерева: %.2f\n", total_weight);

    CHECK(GrB_Matrix_free(&graph));
    CHECK(GrB_Vector_free(&mst_parents));
    CHECK(GrB_finalize());

    return 0;
}