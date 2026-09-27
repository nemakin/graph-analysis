#include <gtest/gtest.h>

#include "graphblas/tc.hpp"
#include "graphblas/utils.hpp"
#include "spla/tc.hpp"
#include "spla/utils.hpp"

#include <cstdint>
#include <utility>

namespace {

int count_spla(const spla::ref_ptr<spla::Matrix> &adjacency, bool triangular) {
  auto product = spla::Matrix::make(adjacency->get_n_rows(),
                                    adjacency->get_n_cols(), spla::INT);
  int count = -1;
  if (triangular) {
    tc_spla::sandia(count, adjacency, product);
  } else {
    tc_spla::burkhardt(count, adjacency, product);
  }
  return count;
}

uint64_t count_graphblas(GrB_Matrix adjacency, bool triangular) {
  GrB_Index n;
  EXPECT_EQ(GrB_Matrix_nrows(&n, adjacency), GrB_SUCCESS);
  GrB_Matrix workspace = nullptr;
  EXPECT_EQ(GrB_Matrix_new(&workspace, GrB_UINT64, n, n), GrB_SUCCESS);
  const uint64_t count = triangular
                             ? tc_graphblas::sandia(adjacency, workspace)
                             : tc_graphblas::burkhardt(adjacency, workspace);
  EXPECT_EQ(GrB_Matrix_free(&workspace), GrB_SUCCESS);
  return count;
}

TEST(SplaTriangleCounting, TwoTrianglesFromDataset) {
  spla::Library::get()->set_force_no_acceleration(true);
  EXPECT_EQ(count_spla(spla_utils::load_graph(TC_TEST_DATASET, false), false), 2);
  EXPECT_EQ(count_spla(spla_utils::load_graph(TC_TEST_DATASET, true), true), 2);
}

TEST(SplaTriangleCounting, PathHasNoTriangles) {
  spla::Library::get()->set_force_no_acceleration(true);
  auto full = spla::Matrix::make(4, 4, spla::INT);
  auto lower = spla::Matrix::make(4, 4, spla::INT);
  for (auto [u, v] : {std::pair{0u, 1u}, std::pair{1u, 2u},
                       std::pair{2u, 3u}}) {
    ASSERT_EQ(full->set_int(u, v, 1), spla::Status::Ok);
    ASSERT_EQ(full->set_int(v, u, 1), spla::Status::Ok);
    ASSERT_EQ(lower->set_int(v, u, 1), spla::Status::Ok);
  }
  EXPECT_EQ(count_spla(full, false), 0);
  EXPECT_EQ(count_spla(lower, true), 0);
}

class GraphblasTriangleCounting : public ::testing::Test {
protected:
  static void SetUpTestSuite() { ASSERT_EQ(GrB_init(GrB_BLOCKING), GrB_SUCCESS); }
  static void TearDownTestSuite() { EXPECT_EQ(GrB_finalize(), GrB_SUCCESS); }
};

TEST_F(GraphblasTriangleCounting, TwoTrianglesFromDataset) {
  GrB_Matrix full = graphblas_utils::load_graph(TC_TEST_DATASET, false);
  ASSERT_NE(full, nullptr);
  EXPECT_EQ(count_graphblas(full, false), 2u);
  EXPECT_EQ(GrB_Matrix_free(&full), GrB_SUCCESS);

  GrB_Matrix upper = graphblas_utils::load_graph(TC_TEST_DATASET, true);
  ASSERT_NE(upper, nullptr);
  EXPECT_EQ(count_graphblas(upper, true), 2u);
  EXPECT_EQ(GrB_Matrix_free(&upper), GrB_SUCCESS);
}

TEST_F(GraphblasTriangleCounting, PathHasNoTriangles) {
  GrB_Matrix full = nullptr;
  GrB_Matrix upper = nullptr;
  ASSERT_EQ(GrB_Matrix_new(&full, GrB_UINT64, 4, 4), GrB_SUCCESS);
  ASSERT_EQ(GrB_Matrix_new(&upper, GrB_UINT64, 4, 4), GrB_SUCCESS);
  for (auto [u, v] : {std::pair{0u, 1u}, std::pair{1u, 2u},
                       std::pair{2u, 3u}}) {
    ASSERT_EQ(GrB_Matrix_setElement_UINT64(full, 1, u, v), GrB_SUCCESS);
    ASSERT_EQ(GrB_Matrix_setElement_UINT64(full, 1, v, u), GrB_SUCCESS);
    ASSERT_EQ(GrB_Matrix_setElement_UINT64(upper, 1, u, v), GrB_SUCCESS);
  }
  EXPECT_EQ(count_graphblas(full, false), 0u);
  EXPECT_EQ(count_graphblas(upper, true), 0u);
  EXPECT_EQ(GrB_Matrix_free(&full), GrB_SUCCESS);
  EXPECT_EQ(GrB_Matrix_free(&upper), GrB_SUCCESS);
}

} // namespace
