# CMake generated Testfile for 
# Source directory: /home/alex/xfemm_cuda/xfemm/cfemm/hpproc/test
# Build directory: /home/alex/xfemm_cuda/xfemm/build-cudss-session/hpproc/test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(hpproc_Temp0 "/home/alex/xfemm_cuda/xfemm/cfemm/bin/hpproc-test" "Temp0.anh" "Temp0.out")
set_tests_properties(hpproc_Temp0 PROPERTIES  LABELS "heatflow;postprocessor" _BACKTRACE_TRIPLES "/home/alex/xfemm_cuda/xfemm/cfemm/hpproc/test/CMakeLists.txt;13;add_test;/home/alex/xfemm_cuda/xfemm/cfemm/hpproc/test/CMakeLists.txt;30;test_hpproc;/home/alex/xfemm_cuda/xfemm/cfemm/hpproc/test/CMakeLists.txt;0;")
add_test(hpproc_Temp0.out.check "/usr/bin/cmake" "-E" "compare_files" "Temp0.out" "Temp0.out.check")
set_tests_properties(hpproc_Temp0.out.check PROPERTIES  DEPENDS "hpproc_Temp0" LABELS "heatflow" _BACKTRACE_TRIPLES "/home/alex/xfemm_cuda/xfemm/cfemm/hpproc/test/CMakeLists.txt;20;add_test;/home/alex/xfemm_cuda/xfemm/cfemm/hpproc/test/CMakeLists.txt;30;test_hpproc;/home/alex/xfemm_cuda/xfemm/cfemm/hpproc/test/CMakeLists.txt;0;")
add_test(hpproc_Temp1 "/home/alex/xfemm_cuda/xfemm/cfemm/bin/hpproc-test" "Temp1.anh" "Temp1.out")
set_tests_properties(hpproc_Temp1 PROPERTIES  LABELS "heatflow;postprocessor" _BACKTRACE_TRIPLES "/home/alex/xfemm_cuda/xfemm/cfemm/hpproc/test/CMakeLists.txt;13;add_test;/home/alex/xfemm_cuda/xfemm/cfemm/hpproc/test/CMakeLists.txt;31;test_hpproc;/home/alex/xfemm_cuda/xfemm/cfemm/hpproc/test/CMakeLists.txt;0;")
add_test(hpproc_Temp1.out.check "/usr/bin/cmake" "-E" "compare_files" "Temp1.out" "Temp1.out.check")
set_tests_properties(hpproc_Temp1.out.check PROPERTIES  DEPENDS "hpproc_Temp1" LABELS "heatflow" _BACKTRACE_TRIPLES "/home/alex/xfemm_cuda/xfemm/cfemm/hpproc/test/CMakeLists.txt;20;add_test;/home/alex/xfemm_cuda/xfemm/cfemm/hpproc/test/CMakeLists.txt;31;test_hpproc;/home/alex/xfemm_cuda/xfemm/cfemm/hpproc/test/CMakeLists.txt;0;")
