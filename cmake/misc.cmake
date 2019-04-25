FUNCTION(setup_eigen)
    execute_process(
            COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
            RESULT_VARIABLE result
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/eigen-download
    )
    if (result)
        message(FATAL_ERROR "CMake step for eigen failed: ${result}")
    endif ()

    execute_process(
            COMMAND ${CMAKE_COMMAND} --build .
            RESULT_VARIABLE result
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/eigen-download
    )
    if (result)
        message(FATAL_ERROR "Build step for eigen failed: ${result}")
    endif ()
    list(APPEND ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_BINARY_DIR}/eigen-build)
    find_package(Eigen3 REQUIRED)
ENDFUNCTION()

FUNCTION(test_case name)
    set(exe_name test_${name})
    set(utest_name ${exe_name})
    add_executable(${exe_name} ${PROJECT_SOURCE_DIR}/test/${exe_name}.cpp)
    target_link_libraries(${exe_name} ${CS133_HW56_LIB})
    add_test(${utest_name} ${exe_name})
    set_tests_properties(${utest_name} PROPERTIES PASS_REGULAR_EXPRESSION "converged")
ENDFUNCTION()
