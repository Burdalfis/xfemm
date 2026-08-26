# CMake generated Testfile for 
# Source directory: /home/alex/xfemm_cuda/xfemm/cfemm/fpproc
# Build directory: /home/alex/xfemm_cuda/xfemm/build-cudss-session/fpproc
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(fpproc_snapshot_parity "/home/alex/xfemm_cuda/xfemm/cfemm/bin/fpproc-snapshot-test" "/home/alex/xfemm_cuda/xfemm/build-cudss-session/fsolver/test/Temp.ans")
set_tests_properties(fpproc_snapshot_parity PROPERTIES  DEPENDS "fsolver_Temp.solve" LABELS "magnetics;postprocessor;unit" _BACKTRACE_TRIPLES "/home/alex/xfemm_cuda/xfemm/cfemm/fpproc/CMakeLists.txt;21;add_test;/home/alex/xfemm_cuda/xfemm/cfemm/fpproc/CMakeLists.txt;0;")
add_test(fpproc_snapshot_age_parity "/home/alex/xfemm_cuda/xfemm/cfemm/bin/fpproc-snapshot-test" "/home/alex/xfemm_cuda/xfemm/build-cudss-session/femmcli/test/femmcli_antiperiodicBC_AGE_TorqueBenchmark.result.ans")
set_tests_properties(fpproc_snapshot_age_parity PROPERTIES  DEPENDS "femmcli_antiperiodicBC_AGE_TorqueBenchmark.lua" LABELS "magnetics;postprocessor;unit;airgap" _BACKTRACE_TRIPLES "/home/alex/xfemm_cuda/xfemm/cfemm/fpproc/CMakeLists.txt;26;add_test;/home/alex/xfemm_cuda/xfemm/cfemm/fpproc/CMakeLists.txt;0;")
add_test(persistent_motor_session_cudss "/home/alex/xfemm_cuda/xfemm/cfemm/bin/persistent_motor_session_test" "/home/alex/xfemm_cuda/xfemm/cfemm/../mfemm/testing/radial_machine/data/radial_machine_sliding.fem")
set_tests_properties(persistent_motor_session_cudss PROPERTIES  LABELS "magnetics;cuda;cudss;airgap;regression" SKIP_RETURN_CODE "77" _BACKTRACE_TRIPLES "/home/alex/xfemm_cuda/xfemm/cfemm/fpproc/CMakeLists.txt;40;add_test;/home/alex/xfemm_cuda/xfemm/cfemm/fpproc/CMakeLists.txt;0;")
