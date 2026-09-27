#include "graphblas/tc.hpp"
#include "spla/tc.hpp"
#include <GraphBLAS.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

int main(int argc, char *argv[]) {
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " <datasets_folder> <num_iters>"
              << std::endl;
    return 1;
  }
  std::string folder = argv[1];
  int num_iters = std::stoi(argv[2]);
  std::ofstream csv("bench_tc.csv");
  csv << "algo,dataset,time_of_iter" << std::endl;
  for (const auto &entry : std::filesystem::directory_iterator(folder)) {
    if (entry.is_regular_file() && entry.path().extension() == ".txt") {
      std::string dataset = entry.path().filename().string();
      std::string dataset_path = entry.path().string();
      std::cout << "\nRunning benchmarks for dataset: " << dataset << std::endl;

      // SPLA GPU Burkhardt
      auto times_splagpu_burkhardt =
          tc_spla::benchmark(dataset_path.c_str(), false, num_iters, true);
      for (const auto &t : times_splagpu_burkhardt) {
        csv << "SPLAGPU_Burkhardt," << dataset << "," << t << std::endl;
      }
      // SPLA GPU Sandia
      auto times_splagpu_sandia =
          tc_spla::benchmark(dataset_path.c_str(), true, num_iters, true);
      for (const auto &t : times_splagpu_sandia) {
        csv << "SPLAGPU_Sandia," << dataset << "," << t << std::endl;
      }

      // SPLA Burkhardt
      auto times_spla_burkhardt =
          tc_spla::benchmark(dataset_path.c_str(), false, num_iters, false);
      for (const auto &t : times_spla_burkhardt) {
        csv << "SPLA_Burkhardt," << dataset << "," << t << std::endl;
      }
      // SPLA Sandia
      auto times_spla_sandia =
          tc_spla::benchmark(dataset_path.c_str(), true, num_iters, false);
      for (const auto &t : times_spla_sandia) {
        csv << "SPLA_Sandia," << dataset << "," << t << std::endl;
      }

      // GraphBLAS Burkhardt
      auto times_burkhardt =
          tc_graphblas::benchmark(dataset_path.c_str(), false, num_iters);
      for (const auto &t : times_burkhardt) {
        csv << "GB_Burkhardt," << dataset << "," << t << std::endl;
      }
      // GraphBLAS Sandia
      auto times_sandia =
          tc_graphblas::benchmark(dataset_path.c_str(), true, num_iters);
      for (const auto &t : times_sandia) {
        csv << "GB_Sandia," << dataset << "," << t << std::endl;
      }
    }
  }
  return 0;
}
