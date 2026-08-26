# CMake generated Testfile for 
# Source directory: /home/alex/xfemm_cuda/xfemm/cfemm/epproc/test
# Build directory: /home/alex/xfemm_cuda/xfemm/build-cudss-asan/epproc/test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(epproc_test "/home/alex/xfemm_cuda/xfemm/cfemm/bin/epproc-test" "test.res" "test.out")
set_tests_properties(epproc_test PROPERTIES  LABELS "electrostatics;postprocessor" _BACKTRACE_TRIPLES "/home/alex/xfemm_cuda/xfemm/cfemm/epproc/test/CMakeLists.txt;14;add_test;/home/alex/xfemm_cuda/xfemm/cfemm/epproc/test/CMakeLists.txt;31;test_epproc;/home/alex/xfemm_cuda/xfemm/cfemm/epproc/test/CMakeLists.txt;0;")
add_test(epproc_test.out.check "/usr/bin/cmake" "-E" "compare_files" "test.out" "test.out.check")
set_tests_properties(epproc_test.out.check PROPERTIES  DEPENDS "epproc_test" LABELS "electrostatics" _BACKTRACE_TRIPLES "/home/alex/xfemm_cuda/xfemm/cfemm/epproc/test/CMakeLists.txt;21;add_test;/home/alex/xfemm_cuda/xfemm/cfemm/epproc/test/CMakeLists.txt;31;test_epproc;/home/alex/xfemm_cuda/xfemm/cfemm/epproc/test/CMakeLists.txt;0;")
