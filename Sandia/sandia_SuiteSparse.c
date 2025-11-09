#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <GraphBLAS.h>
#include <LAGraph.h>
#include <sys/resource.h>
#include <stdio.h>
#include <time.h>
//--------------------------------------------------------------------
// Функция чтения MatrixMarket файла и создания GrB_Matrix
//--------------------------------------------------------------------
GrB_Matrix read_matrix_market(const char *filename)
{
    FILE *f = fopen(filename, "r");
    if (!f)
    {
        perror("open file");
        exit(1);
    }

    char line[5096];
    do
    {
        if (!fgets(line, sizeof(line), f))
            exit(1);
    } while (line[0] == '%');

    GrB_Index nrows, ncols, nvals;
    GrB_Index max_el;
    sscanf(line, "%llu %llu %llu", &nrows, &ncols, &nvals);
    max_el = nvals;
    for (GrB_Index k = 0; k < nvals; k++)
    {
        unsigned long long i, j;
        fscanf(f, "%zu %zu", &i, &j);
        if (i - 1 > max_el)
            max_el = i - 1;
        if (j - 1 > max_el)
            max_el = j - 1;
    }

    rewind(f);
    do
    {
        if (!fgets(line, sizeof(line), f))
            exit(1);
    } while (line[0] == '%');
    sscanf(line, "%llu %llu %llu", &nrows, &ncols, &nvals);

    GrB_Matrix A;
    GrB_Matrix_new(&A, GrB_UINT32, max_el, max_el);

    unsigned long long cnt = 0;
    for (GrB_Index k = 0; k < nvals; k++)
    {
        unsigned long long i, j;
        fscanf(f, "%llu %llu", &i, &j);
        i--; 
        j--;

        GrB_Matrix_setElement_BOOL(A, true, i, j);
        GrB_Matrix_setElement_BOOL(A, true, j, i);
        cnt += 2;
    }

    fclose(f);
    printf("  Ребер:  %llu\n", cnt);
    printf("  Ребер:  %llu\n", nvals);
    printf("  Ребер:  %llu\n", max_el);

    GrB_Index nrows1, ncols1, nvals1;

    GrB_Matrix_nrows(&nrows1, A);
    GrB_Matrix_ncols(&ncols1, A);
    GrB_Matrix_nvals(&nvals1, A);

    printf("Размер матрицы смежности:\n");
    printf("  Вершин (строк):  %llu\n", nrows1);
    printf("  Вершин (столбцов): %llu\n", ncols1);
    printf("Количество рёбер (ненулевых элементов): %llu\n", nvals1);
    printf("************************.\n");

    return A;
}

//--------------------------------------------------------------------
// Основная функция
//--------------------------------------------------------------------
int main(int argc, char **argv)
{
    if (argc < 2)
    {
        fprintf(stderr, "Usage: %s input.mtx\n", argv[0]);
        return EXIT_FAILURE;
    }

    const char *filename = argv[1];

    LAGraph_Init(NULL);

    GrB_Matrix A = read_matrix_market(filename);
    if (!A)
    {
        fprintf(stderr, "Error reading MatrixMarket file.\n");
        LAGraph_Finalize(NULL);
        return EXIT_FAILURE;
    }

    GrB_Index nrows, ncols, nvals;

    GrB_Matrix_nrows(&nrows, A);
    GrB_Matrix_ncols(&ncols, A);
    GrB_Matrix_nvals(&nvals, A);

    printf("Размер матрицы смежности:\n");
    printf("  Вершин (строк):  %llu\n", (unsigned long long)nrows);
    printf("  Вершин (столбцов): %llu\n", (unsigned long long)ncols);
    printf("Количество рёбер (ненулевых элементов): %llu\n", (unsigned long long)nvals);
    printf("************************.\n");

    LAGraph_Graph G = NULL;
    char msg[LAGRAPH_MSG_LEN];
    if (LAGraph_New(&G, &A, LAGraph_ADJACENCY_UNDIRECTED, msg) != GrB_SUCCESS)
    {
        fprintf(stderr, "LAGraph_New failed: %s\n", msg);
        LAGraph_Finalize(msg);
        return EXIT_FAILURE;
    }

    if (G->is_symmetric_structure)
        printf("✔ Граф симметричен.\n");
    else
        printf("✖ Граф несимметричен.\n");

    GrB_Matrix matr = G->A; 
    GrB_Matrix_nrows(&nrows, matr);
    GrB_Matrix_ncols(&ncols, matr);
    GrB_Matrix_nvals(&nvals, matr);

    printf("Размер матрицы смежности:\n");
    printf("  Вершин (строк):  %llu\n", (unsigned long long)nrows);
    printf("  Вершин (столбцов): %llu\n", (unsigned long long)ncols);
    printf("Количество рёбер (ненулевых элементов): %llu\n", (unsigned long long)nvals);

    printf("************************.\n");

    LAGraph_DeleteSelfEdges(G, msg);

    LAGraph_Cached_NSelfEdges(G, msg);
    LAGraph_Cached_IsSymmetricStructure(G, msg);
    LAGraph_Cached_OutDegree(G, msg);

    GrB_Matrix_nrows(&nrows, matr);
    GrB_Matrix_ncols(&ncols, matr);
    GrB_Matrix_nvals(&nvals, matr);

    printf("Размер матрицы смежности:\n");
    printf("  Вершин (строк):  %llu\n", (unsigned long long)nrows);
    printf("  Вершин (столбцов): %llu\n", (unsigned long long)ncols);
    printf("Количество рёбер (ненулевых элементов): %llu\n", (unsigned long long)nvals);

    if (G->is_symmetric_structure)
        printf("✔ Граф симметричен.\n");
    else
        printf("✖ Граф несимметричен.\n");
    uint64_t ntri = 0;
    LAGr_TriangleCount_Method method = LAGr_TriangleCount_Sandia_LL;
    LAGr_TriangleCount_Presort presort = LAGr_TriangleCount_AutoSort;
    clock_t c1 = clock();
    double start_time = LAGraph_WallClockTime();

    if (LAGr_TriangleCount(&ntri, G, &method, &presort, msg) != GrB_SUCCESS)
    {
        fprintf(stderr, "TriangleCount failed: %s\n", msg);
        LAGraph_Delete(&G, msg);
        LAGraph_Finalize(msg);
        return EXIT_FAILURE;
    }
    clock_t c2 = clock();
    double end_time = LAGraph_WallClockTime();
    printf("CPU time: %.6f s\n", (double)(c2 - c1) / CLOCKS_PER_SEC);
    printf("Время выполнения: %.6f секунд\n", end_time - start_time);

    printf("Triangle count (Sandia_LL): %llu\n", ntri);
    LAGraph_Delete(&G, msg);
    LAGraph_Finalize(msg);
    return EXIT_SUCCESS;
}
