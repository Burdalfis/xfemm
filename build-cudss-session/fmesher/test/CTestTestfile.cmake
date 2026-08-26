# CMake generated Testfile for 
# Source directory: /home/alex/xfemm_cuda/xfemm/cfemm/fmesher/test
# Build directory: /home/alex/xfemm_cuda/xfemm/build-cudss-session/fmesher/test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(fmesher_Temp.fem "/home/alex/xfemm_cuda/xfemm/cfemm/bin/fmesher" "/home/alex/xfemm_cuda/xfemm/build-cudss-session/fmesher/test/Temp.fem.work/Temp.fem")
set_tests_properties(fmesher_Temp.fem PROPERTIES  LABELS "magnetics;mesher" WORKING_DIRECTORY "/home/alex/xfemm_cuda/xfemm/build-cudss-session/fmesher/test/Temp.fem.work" _BACKTRACE_TRIPLES "/home/alex/xfemm_cuda/xfemm/cfemm/fmesher/test/CMakeLists.txt;5;add_test;/home/alex/xfemm_cuda/xfemm/cfemm/fmesher/test/CMakeLists.txt;14;test_fmesher;/home/alex/xfemm_cuda/xfemm/cfemm/fmesher/test/CMakeLists.txt;0;")
add_test(fmesher_split_seg_err_test.fem "/home/alex/xfemm_cuda/xfemm/cfemm/bin/fmesher" "/home/alex/xfemm_cuda/xfemm/build-cudss-session/fmesher/test/split_seg_err_test.fem.work/split_seg_err_test.fem")
set_tests_properties(fmesher_split_seg_err_test.fem PROPERTIES  LABELS "magnetics;mesher" WORKING_DIRECTORY "/home/alex/xfemm_cuda/xfemm/build-cudss-session/fmesher/test/split_seg_err_test.fem.work" _BACKTRACE_TRIPLES "/home/alex/xfemm_cuda/xfemm/cfemm/fmesher/test/CMakeLists.txt;5;add_test;/home/alex/xfemm_cuda/xfemm/cfemm/fmesher/test/CMakeLists.txt;15;test_fmesher;/home/alex/xfemm_cuda/xfemm/cfemm/fmesher/test/CMakeLists.txt;0;")
add_test(fmesher_selected_backend "/home/alex/xfemm_cuda/xfemm/cfemm/bin/fmesher" "--version")
set_tests_properties(fmesher_selected_backend PROPERTIES  LABELS "mesher;Triangle" PASS_REGULAR_EXPRESSION "default mesher backend Triangle" _BACKTRACE_TRIPLES "/home/alex/xfemm_cuda/xfemm/cfemm/fmesher/test/CMakeLists.txt;17;add_test;/home/alex/xfemm_cuda/xfemm/cfemm/fmesher/test/CMakeLists.txt;0;")
add_test(fmesher_tangle_mesh_converter "/home/alex/xfemm_cuda/xfemm/cfemm/bin/tangle_mesh_converter_test")
set_tests_properties(fmesher_tangle_mesh_converter PROPERTIES  LABELS "mesher;Tangle;unit" _BACKTRACE_TRIPLES "/home/alex/xfemm_cuda/xfemm/cfemm/fmesher/test/CMakeLists.txt;29;add_test;/home/alex/xfemm_cuda/xfemm/cfemm/fmesher/test/CMakeLists.txt;0;")
add_test(fmesher_tangle_backend_execution "/home/alex/xfemm_cuda/xfemm/cfemm/bin/tangle_backend_execution_test")
set_tests_properties(fmesher_tangle_backend_execution PROPERTIES  LABELS "mesher;Tangle;unit" _BACKTRACE_TRIPLES "/home/alex/xfemm_cuda/xfemm/cfemm/fmesher/test/CMakeLists.txt;35;add_test;/home/alex/xfemm_cuda/xfemm/cfemm/fmesher/test/CMakeLists.txt;0;")
add_test(fmesher_Triangle_backend_solver_mesh "/home/alex/xfemm_cuda/xfemm/cfemm/bin/fmesher_backend_test" "Triangle" "/home/alex/xfemm_cuda/xfemm/cfemm/fmesher/test/Temp.fem" "/home/alex/xfemm_cuda/xfemm/cfemm/fmesher/test/../../femmcli/test/femmcli_antiperiodicBC_AGE_TorqueBenchmark.fem")
set_tests_properties(fmesher_Triangle_backend_solver_mesh PROPERTIES  LABELS "magnetics;mesher;Triangle" _BACKTRACE_TRIPLES "/home/alex/xfemm_cuda/xfemm/cfemm/fmesher/test/CMakeLists.txt;40;add_test;/home/alex/xfemm_cuda/xfemm/cfemm/fmesher/test/CMakeLists.txt;0;")
add_test(fmesher_Tangle_backend_solver_mesh "/home/alex/xfemm_cuda/xfemm/cfemm/bin/fmesher_backend_test" "Tangle" "/home/alex/xfemm_cuda/xfemm/cfemm/fmesher/test/Temp.fem" "/home/alex/xfemm_cuda/xfemm/cfemm/fmesher/test/../../femmcli/test/femmcli_antiperiodicBC_AGE_TorqueBenchmark.fem")
set_tests_properties(fmesher_Tangle_backend_solver_mesh PROPERTIES  LABELS "magnetics;mesher;Tangle" _BACKTRACE_TRIPLES "/home/alex/xfemm_cuda/xfemm/cfemm/fmesher/test/CMakeLists.txt;40;add_test;/home/alex/xfemm_cuda/xfemm/cfemm/fmesher/test/CMakeLists.txt;0;")
