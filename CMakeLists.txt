cmake_minimum_required(VERSION 3.5)
project(cs133hw56)

set(CS133_HW56_VERSION_MAJOR 1)
set(CS133_HW56_VERSION_MINOR 0)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE Debug)

include(CTest)
include(${PROJECT_SOURCE_DIR}/cmake/misc.cmake)
include_directories(${PROJECT_SOURCE_DIR}/include)

configure_file(${PROJECT_SOURCE_DIR}/cmake/get_eigen.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/eigen-download/CMakeLists.txt)
setup_eigen()

set(CS133_HW56_LIB ${PROJECT_NAME})
set(CS133_HW56_SRC
        ${PROJECT_SOURCE_DIR}/src/AllFunctors.cpp
        ${PROJECT_SOURCE_DIR}/src/functions.cpp)
add_library(${CS133_HW56_LIB} ${CS133_HW56_SRC})
target_link_libraries(${CS133_HW56_LIB} Eigen3::Eigen)

if (NOT DEFINED ENABLE_TESTS)
    set(ENABLE_TESTS "[a-z_]*")
endif ()
file(GLOB CS133_HW56_TEST_CASES ${PROJECT_SOURCE_DIR}/test/*.cpp)
foreach (test_src ${CS133_HW56_TEST_CASES})
    string(REPLACE "${PROJECT_SOURCE_DIR}/test/test_" "" test_src ${test_src})
    string(REPLACE ".cpp" "" test_src ${test_src})
    if (${test_src} MATCHES ${ENABLE_TESTS})
        test_case(${test_src})
    endif ()
endforeach ()
