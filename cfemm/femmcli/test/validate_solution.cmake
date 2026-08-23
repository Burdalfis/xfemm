if(NOT EXISTS "${SOLUTION_FILE}")
    message(FATAL_ERROR "Solution file was not created: ${SOLUTION_FILE}")
endif()
file(SIZE "${SOLUTION_FILE}" _solution_size)
if(_solution_size LESS 100)
    message(FATAL_ERROR "Solution file is unexpectedly small: ${_solution_size} bytes")
endif()
# Search the raw text so FEMM's CRLF output and CMake's line handling do not
# make the section check platform-dependent.
file(READ "${SOLUTION_FILE}" _solution_text)
string(FIND "${_solution_text}" "[Solution]" _solution_marker)
if(_solution_marker EQUAL -1)
    message(FATAL_ERROR "Solution file has no [Solution] section")
endif()
