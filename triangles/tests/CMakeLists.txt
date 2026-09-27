set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../spla/deps/gtest
                 ${CMAKE_CURRENT_BINARY_DIR}/googletest)

add_executable(tc_tests algorithms_test.cpp)
target_link_libraries(tc_tests PRIVATE tc_algorithms gtest_main)
target_compile_definitions(tc_tests PRIVATE
    TC_TEST_DATASET="${CMAKE_CURRENT_SOURCE_DIR}/data/two_triangles.txt")

add_test(NAME tc_algorithms COMMAND tc_tests)
