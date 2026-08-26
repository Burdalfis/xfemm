# CMake generated Testfile for 
# Source directory: /home/alex/xfemm_cuda/xfemm/cfemm/libfemm/test
# Build directory: /home/alex/xfemm_cuda/xfemm/build-cudss-session/libfemm/test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(cspars_multiply "/home/alex/xfemm_cuda/xfemm/cfemm/bin/cspars_multiply_test")
set_tests_properties(cspars_multiply PROPERTIES  LABELS "libfemm;unit" _BACKTRACE_TRIPLES "/home/alex/xfemm_cuda/xfemm/cfemm/libfemm/test/CMakeLists.txt;4;add_test;/home/alex/xfemm_cuda/xfemm/cfemm/libfemm/test/CMakeLists.txt;0;")
add_test(spars_compact "/home/alex/xfemm_cuda/xfemm/cfemm/bin/spars_compact_test")
set_tests_properties(spars_compact PROPERTIES  LABELS "libfemm;unit" _BACKTRACE_TRIPLES "/home/alex/xfemm_cuda/xfemm/cfemm/libfemm/test/CMakeLists.txt;10;add_test;/home/alex/xfemm_cuda/xfemm/cfemm/libfemm/test/CMakeLists.txt;0;")
add_test(spars_compact_mixed16 "/home/alex/xfemm_cuda/xfemm/cfemm/bin/spars_compact_test")
set_tests_properties(spars_compact_mixed16 PROPERTIES  ENVIRONMENT "XFEMM_PCG_COLUMN_INDEX=mixed16" LABELS "libfemm;unit" _BACKTRACE_TRIPLES "/home/alex/xfemm_cuda/xfemm/cfemm/libfemm/test/CMakeLists.txt;13;add_test;/home/alex/xfemm_cuda/xfemm/cfemm/libfemm/test/CMakeLists.txt;0;")
add_test(spars_compact_row16 "/home/alex/xfemm_cuda/xfemm/cfemm/bin/spars_compact_test")
set_tests_properties(spars_compact_row16 PROPERTIES  ENVIRONMENT "XFEMM_PCG_COLUMN_INDEX=row16" LABELS "libfemm;unit" _BACKTRACE_TRIPLES "/home/alex/xfemm_cuda/xfemm/cfemm/libfemm/test/CMakeLists.txt;19;add_test;/home/alex/xfemm_cuda/xfemm/cfemm/libfemm/test/CMakeLists.txt;0;")
add_test(spars_compact_parallel "/home/alex/xfemm_cuda/xfemm/cfemm/bin/spars_compact_test")
set_tests_properties(spars_compact_parallel PROPERTIES  ENVIRONMENT "XFEMM_NUM_THREADS=4" LABELS "libfemm;unit" _BACKTRACE_TRIPLES "/home/alex/xfemm_cuda/xfemm/cfemm/libfemm/test/CMakeLists.txt;26;add_test;/home/alex/xfemm_cuda/xfemm/cfemm/libfemm/test/CMakeLists.txt;0;")
add_test(femm_problem_geometry "/home/alex/xfemm_cuda/xfemm/cfemm/bin/femm_problem_geometry_test")
set_tests_properties(femm_problem_geometry PROPERTIES  LABELS "libfemm;unit" _BACKTRACE_TRIPLES "/home/alex/xfemm_cuda/xfemm/cfemm/libfemm/test/CMakeLists.txt;36;add_test;/home/alex/xfemm_cuda/xfemm/cfemm/libfemm/test/CMakeLists.txt;0;")
add_test(material_curve "/home/alex/xfemm_cuda/xfemm/cfemm/bin/material_curve_test")
set_tests_properties(material_curve PROPERTIES  LABELS "libfemm;unit" _BACKTRACE_TRIPLES "/home/alex/xfemm_cuda/xfemm/cfemm/libfemm/test/CMakeLists.txt;42;add_test;/home/alex/xfemm_cuda/xfemm/cfemm/libfemm/test/CMakeLists.txt;0;")
