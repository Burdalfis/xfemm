# CMake generated Testfile for 
# Source directory: /home/alex/xfemm_cuda/xfemm/cfemm/fsolver/test
# Build directory: /home/alex/xfemm_cuda/xfemm/build-cudss-session/fsolver/test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(fsolver_Temp.setup "/usr/bin/cmake" "-E" "copy" "/home/alex/xfemm_cuda/xfemm/cfemm/fsolver/test/Temp.edge" "/home/alex/xfemm_cuda/xfemm/cfemm/fsolver/test/Temp.ele" "/home/alex/xfemm_cuda/xfemm/cfemm/fsolver/test/Temp.node" "/home/alex/xfemm_cuda/xfemm/cfemm/fsolver/test/Temp.pbc" ".")
set_tests_properties(fsolver_Temp.setup PROPERTIES  _BACKTRACE_TRIPLES "/home/alex/xfemm_cuda/xfemm/cfemm/fsolver/test/CMakeLists.txt;18;add_test;/home/alex/xfemm_cuda/xfemm/cfemm/fsolver/test/CMakeLists.txt;64;test_fsolver;/home/alex/xfemm_cuda/xfemm/cfemm/fsolver/test/CMakeLists.txt;0;")
add_test(fsolver_Temp.solve "/home/alex/xfemm_cuda/xfemm/cfemm/bin/fsolver" "Temp")
set_tests_properties(fsolver_Temp.solve PROPERTIES  DEPENDS "fsolver_Temp.setup" LABELS "magnetics;solver" _BACKTRACE_TRIPLES "/home/alex/xfemm_cuda/xfemm/cfemm/fsolver/test/CMakeLists.txt;36;add_test;/home/alex/xfemm_cuda/xfemm/cfemm/fsolver/test/CMakeLists.txt;64;test_fsolver;/home/alex/xfemm_cuda/xfemm/cfemm/fsolver/test/CMakeLists.txt;0;")
add_test(fsolver_Temp.check "/usr/bin/cmake" "-E" "compare_files" "Temp.ans" "Temp.ans.check")
set_tests_properties(fsolver_Temp.check PROPERTIES  DEPENDS "fsolver_Temp.solve" DISABLED "TRUE" LABELS "magnetics" _BACKTRACE_TRIPLES "/home/alex/xfemm_cuda/xfemm/cfemm/fsolver/test/CMakeLists.txt;51;add_test;/home/alex/xfemm_cuda/xfemm/cfemm/fsolver/test/CMakeLists.txt;64;test_fsolver;/home/alex/xfemm_cuda/xfemm/cfemm/fsolver/test/CMakeLists.txt;0;")
add_test(fsolver_Temp1.mesh "/home/alex/xfemm_cuda/xfemm/cfemm/bin/fmesher" "Temp1.fem")
set_tests_properties(fsolver_Temp1.mesh PROPERTIES  LABELS "magnetics;mesher" _BACKTRACE_TRIPLES "/home/alex/xfemm_cuda/xfemm/cfemm/fsolver/test/CMakeLists.txt;29;add_test;/home/alex/xfemm_cuda/xfemm/cfemm/fsolver/test/CMakeLists.txt;65;test_fsolver;/home/alex/xfemm_cuda/xfemm/cfemm/fsolver/test/CMakeLists.txt;0;")
add_test(fsolver_Temp1.solve "/home/alex/xfemm_cuda/xfemm/cfemm/bin/fsolver" "Temp1")
set_tests_properties(fsolver_Temp1.solve PROPERTIES  DEPENDS "fsolver_Temp1.mesh" LABELS "magnetics;solver" _BACKTRACE_TRIPLES "/home/alex/xfemm_cuda/xfemm/cfemm/fsolver/test/CMakeLists.txt;36;add_test;/home/alex/xfemm_cuda/xfemm/cfemm/fsolver/test/CMakeLists.txt;65;test_fsolver;/home/alex/xfemm_cuda/xfemm/cfemm/fsolver/test/CMakeLists.txt;0;")
add_test(fsolver_Temp1.check "/usr/bin/cmake" "-E" "compare_files" "Temp1.ans" "Temp1.ans.check")
set_tests_properties(fsolver_Temp1.check PROPERTIES  DEPENDS "fsolver_Temp1.solve" DISABLED "TRUE" LABELS "magnetics" _BACKTRACE_TRIPLES "/home/alex/xfemm_cuda/xfemm/cfemm/fsolver/test/CMakeLists.txt;51;add_test;/home/alex/xfemm_cuda/xfemm/cfemm/fsolver/test/CMakeLists.txt;65;test_fsolver;/home/alex/xfemm_cuda/xfemm/cfemm/fsolver/test/CMakeLists.txt;0;")
add_test(analysis_session "/home/alex/xfemm_cuda/xfemm/cfemm/bin/analysis_session_test")
set_tests_properties(analysis_session PROPERTIES  LABELS "magnetics;solver;unit" _BACKTRACE_TRIPLES "/home/alex/xfemm_cuda/xfemm/cfemm/fsolver/test/CMakeLists.txt;69;add_test;/home/alex/xfemm_cuda/xfemm/cfemm/fsolver/test/CMakeLists.txt;0;")
