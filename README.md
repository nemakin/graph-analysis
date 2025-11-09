# Graph Analysis - Algorithm Implementation Comparison

This repository contains implementations and comparisons of graph algorithms (Prim's MST and Sandia's Delta-Stepping SSSP) using different GraphBLAS implementations:

- [GBTL (Graph BLAS Template Library)](https://github.com/cmu-sei/gbtl)
- [SuiteSparse:GraphBLAS](https://github.com/DrTimothyAldenDavis/GraphBLAS)
- [SPLA (Sparse Linear Algebra)](https://github.com/JetBrains-Research/spla)

## Repository Structure

```
.
├── gbtl/               # GBTL library and its implementation
├── prim/              # Prim's algorithm implementations
│   ├── prim_spla.cpp
│   └── prim_SuiteSparse.c
├── Sandia/            # Sandia's Delta-Stepping SSSP implementations
│   ├── sandia_spla.cpp
│   └── sandia_SuiteSparse.c
├── spla/              # SPLA library
└── SuiteSparse/       # SuiteSparse library
```

## Algorithms

### Prim's Minimum Spanning Tree (MST)
- Implementation using GBTL
- Implementation using SPLA
- Implementation using SuiteSparse:GraphBLAS

### Sandia's Delta-Stepping Single-Source Shortest Path (SSSP)
- Implementation using SPLA
- Implementation using SuiteSparse:GraphBLAS

## Libraries Used

### GBTL
The GraphBLAS Template Library (GBTL) is a C++ template library that implements the GraphBLAS standard. It provides a reference implementation for the GraphBLAS specification.

### SuiteSparse:GraphBLAS
SuiteSparse:GraphBLAS is a complete implementation of the GraphBLAS standard, which provides a powerful and expressive framework for creating graph algorithms based on the mathematics of sparse matrix operations.

### SPLA
SPLA (Sparse Linear Algebra) is a high-performance library for sparse linear algebra computations, with support for various backends including CPU and OpenCL.

## Dataset

The repository includes a sample dataset:
- `datasets/USA_road_d_COL.gr`: Road network of Colorado

## Building and Running

Each implementation has its own build requirements. Please refer to the respective library documentation for detailed build instructions:

- [GBTL Build Instructions](https://github.com/cmu-sei/gbtl)
- [SuiteSparse Build Instructions](https://github.com/DrTimothyAldenDavis/GraphBLAS)
- [SPLA Build Instructions](https://github.com/JetBrains-Research/spla)

## License

This project includes components from different libraries, each with its own license:
- GBTL: Apache 2.0
- SuiteSparse: Various (primarily GPL v2+)
- SPLA: MIT License

Please refer to individual library folders for their specific licenses.