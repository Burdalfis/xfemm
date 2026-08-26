# CMake generated Testfile for 
# Source directory: /home/alex/xfemm_cuda/xfemm/cfemm/hsolver/test
# Build directory: /home/alex/xfemm_cuda/xfemm/build-cudss-session/hsolver/test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(hsolver_Temp0.mesh "/home/alex/xfemm_cuda/xfemm/cfemm/bin/fmesher" "Temp0.feh")
set_tests_properties(hsolver_Temp0.mesh PROPERTIES  LABELS "heatflow;mesher" _BACKTRACE_TRIPLES "/home/alex/xfemm_cuda/xfemm/cfemm/hsolver/test/CMakeLists.txt;15;add_test;/home/alex/xfemm_cuda/xfemm/cfemm/hsolver/test/CMakeLists.txt;44;test_hsolver;/home/alex/xfemm_cuda/xfemm/cfemm/hsolver/test/CMakeLists.txt;0;")
add_test(hsolver_Temp0.solve "/home/alex/xfemm_cuda/xfemm/cfemm/bin/hsolver" "Temp0")
set_tests_properties(hsolver_Temp0.solve PROPERTIES  DEPENDS "hsolver_Temp0.mesh" LABELS "heatflow;solver" _BACKTRACE_TRIPLES "/home/alex/xfemm_cuda/xfemm/cfemm/hsolver/test/CMakeLists.txt;21;add_test;/home/alex/xfemm_cuda/xfemm/cfemm/hsolver/test/CMakeLists.txt;44;test_hsolver;/home/alex/xfemm_cuda/xfemm/cfemm/hsolver/test/CMakeLists.txt;0;")
add_test(hsolver_Temp0.check "/usr/bin/cmake" "-E" "compare_files" "Temp0.anh" "Temp0.anh.check")
set_tests_properties(hsolver_Temp0.check PROPERTIES  DEPENDS "hsolver_Temp0.solve" DISABLED "TRUE" LABELS "heatflow" _BACKTRACE_TRIPLES "/home/alex/xfemm_cuda/xfemm/cfemm/hsolver/test/CMakeLists.txt;32;add_test;/home/alex/xfemm_cuda/xfemm/cfemm/hsolver/test/CMakeLists.txt;44;test_hsolver;/home/alex/xfemm_cuda/xfemm/cfemm/hsolver/test/CMakeLists.txt;0;")
add_test(hsolver_Temp1.mesh "/home/alex/xfemm_cuda/xfemm/cfemm/bin/fmesher" "Temp1.feh")
set_tests_properties(hsolver_Temp1.mesh PROPERTIES  LABELS "heatflow;mesher" _BACKTRACE_TRIPLES "/home/alex/xfemm_cuda/xfemm/cfemm/hsolver/test/CMakeLists.txt;15;add_test;/home/alex/xfemm_cuda/xfemm/cfemm/hsolver/test/CMakeLists.txt;45;test_hsolver;/home/alex/xfemm_cuda/xfemm/cfemm/hsolver/test/CMakeLists.txt;0;")
add_test(hsolver_Temp1.solve "/home/alex/xfemm_cuda/xfemm/cfemm/bin/hsolver" "Temp1")
set_tests_properties(hsolver_Temp1.solve PROPERTIES  DEPENDS "hsolver_Temp1.mesh" LABELS "heatflow;solver" _BACKTRACE_TRIPLES "/home/alex/xfemm_cuda/xfemm/cfemm/hsolver/test/CMakeLists.txt;21;add_test;/home/alex/xfemm_cuda/xfemm/cfemm/hsolver/test/CMakeLists.txt;45;test_hsolver;/home/alex/xfemm_cuda/xfemm/cfemm/hsolver/test/CMakeLists.txt;0;")
add_test(hsolver_Temp1.check "/usr/bin/cmake" "-E" "compare_files" "Temp1.anh" "Temp1.anh.check")
set_tests_properties(hsolver_Temp1.check PROPERTIES  DEPENDS "hsolver_Temp1.solve" DISABLED "TRUE" LABELS "heatflow" _BACKTRACE_TRIPLES "/home/alex/xfemm_cuda/xfemm/cfemm/hsolver/test/CMakeLists.txt;32;add_test;/home/alex/xfemm_cuda/xfemm/cfemm/hsolver/test/CMakeLists.txt;45;test_hsolver;/home/alex/xfemm_cuda/xfemm/cfemm/hsolver/test/CMakeLists.txt;0;")
