# CMake generated Testfile for 
# Source directory: /home/alex/xfemm_cuda/xfemm/cfemm/esolver/test
# Build directory: /home/alex/xfemm_cuda/xfemm/build-cudss-session/esolver/test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(esolver_test.mesh "/home/alex/xfemm_cuda/xfemm/cfemm/bin/fmesher" "test.fee")
set_tests_properties(esolver_test.mesh PROPERTIES  LABELS "electrostatics;mesher" _BACKTRACE_TRIPLES "/home/alex/xfemm_cuda/xfemm/cfemm/esolver/test/CMakeLists.txt;15;add_test;/home/alex/xfemm_cuda/xfemm/cfemm/esolver/test/CMakeLists.txt;47;test_esolver;/home/alex/xfemm_cuda/xfemm/cfemm/esolver/test/CMakeLists.txt;0;")
add_test(esolver_test.solve "/home/alex/xfemm_cuda/xfemm/cfemm/bin/esolver" "test")
set_tests_properties(esolver_test.solve PROPERTIES  DEPENDS "esolver_test.mesh" LABELS "electrostatics;solver" _BACKTRACE_TRIPLES "/home/alex/xfemm_cuda/xfemm/cfemm/esolver/test/CMakeLists.txt;21;add_test;/home/alex/xfemm_cuda/xfemm/cfemm/esolver/test/CMakeLists.txt;47;test_esolver;/home/alex/xfemm_cuda/xfemm/cfemm/esolver/test/CMakeLists.txt;0;")
add_test(esolver_test.check "/usr/bin/cmake" "-E" "compare_files" "test.res" "test.res.check")
set_tests_properties(esolver_test.check PROPERTIES  DEPENDS "esolver_test.solve" DISABLED "TRUE" LABELS "electrostatics" _BACKTRACE_TRIPLES "/home/alex/xfemm_cuda/xfemm/cfemm/esolver/test/CMakeLists.txt;32;add_test;/home/alex/xfemm_cuda/xfemm/cfemm/esolver/test/CMakeLists.txt;47;test_esolver;/home/alex/xfemm_cuda/xfemm/cfemm/esolver/test/CMakeLists.txt;0;")
