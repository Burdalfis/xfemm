add_executable(cspars_multiply_test cspars_multiply_test.cpp)
target_link_libraries(cspars_multiply_test PRIVATE femm)

add_test(NAME cspars_multiply COMMAND cspars_multiply_test)
set_tests_properties(cspars_multiply PROPERTIES LABELS "libfemm;unit")

add_executable(spars_compact_test spars_compact_test.cpp)
target_link_libraries(spars_compact_test PRIVATE femm)

add_test(NAME spars_compact COMMAND spars_compact_test)
set_tests_properties(spars_compact PROPERTIES LABELS "libfemm;unit")

add_test(NAME spars_compact_mixed16 COMMAND spars_compact_test)
set_tests_properties(spars_compact_mixed16 PROPERTIES
    LABELS "libfemm;unit"
    ENVIRONMENT "XFEMM_PCG_COLUMN_INDEX=mixed16"
    )

add_test(NAME spars_compact_row16 COMMAND spars_compact_test)
set_tests_properties(spars_compact_row16 PROPERTIES
    LABELS "libfemm;unit"
    ENVIRONMENT "XFEMM_PCG_COLUMN_INDEX=row16"
    )

if(XFEMM_ENABLE_OPENMP AND OpenMP_CXX_FOUND)
    add_test(NAME spars_compact_parallel COMMAND spars_compact_test)
    set_tests_properties(spars_compact_parallel PROPERTIES
        LABELS "libfemm;unit"
        ENVIRONMENT "XFEMM_NUM_THREADS=4"
        )
endif()

add_executable(femm_problem_geometry_test femm_problem_geometry_test.cpp)
target_link_libraries(femm_problem_geometry_test PRIVATE femm)

add_test(NAME femm_problem_geometry COMMAND femm_problem_geometry_test)
set_tests_properties(femm_problem_geometry PROPERTIES LABELS "libfemm;unit")

add_executable(material_curve_test material_curve_test.cpp)
target_link_libraries(material_curve_test PRIVATE femm)

add_test(NAME material_curve COMMAND material_curve_test)
set_tests_properties(material_curve PROPERTIES LABELS "libfemm;unit")

if(XFEMM_USE_CUDSS)
    add_executable(cuda_planar_material_test cuda_planar_material_test.cpp)
    target_link_libraries(cuda_planar_material_test PRIVATE femm)
    add_test(NAME cuda_planar_material COMMAND cuda_planar_material_test)
    set_tests_properties(cuda_planar_material PROPERTIES
        LABELS "libfemm;cuda;unit" SKIP_RETURN_CODE 77)
endif()
