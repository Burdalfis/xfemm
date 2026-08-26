# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file LICENSE.rst or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION ${CMAKE_VERSION}) # this file comes with cmake

# If CMAKE_DISABLE_SOURCE_CHANGES is set to true and the source directory is an
# existing directory in our source tree, calling file(MAKE_DIRECTORY) on it
# would cause a fatal error, even though it would be a no-op.
if(NOT EXISTS "/home/alex/xfemm_cuda/xfemm/build-cudss-session/_deps/tangle_upstream-src")
  file(MAKE_DIRECTORY "/home/alex/xfemm_cuda/xfemm/build-cudss-session/_deps/tangle_upstream-src")
endif()
file(MAKE_DIRECTORY
  "/home/alex/xfemm_cuda/xfemm/build-cudss-session/_deps/tangle_upstream-build"
  "/home/alex/xfemm_cuda/xfemm/build-cudss-session/_deps/tangle_upstream-subbuild/tangle_upstream-populate-prefix"
  "/home/alex/xfemm_cuda/xfemm/build-cudss-session/_deps/tangle_upstream-subbuild/tangle_upstream-populate-prefix/tmp"
  "/home/alex/xfemm_cuda/xfemm/build-cudss-session/_deps/tangle_upstream-subbuild/tangle_upstream-populate-prefix/src/tangle_upstream-populate-stamp"
  "/home/alex/xfemm_cuda/xfemm/build-cudss-session/_deps/tangle_upstream-subbuild/tangle_upstream-populate-prefix/src"
  "/home/alex/xfemm_cuda/xfemm/build-cudss-session/_deps/tangle_upstream-subbuild/tangle_upstream-populate-prefix/src/tangle_upstream-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/home/alex/xfemm_cuda/xfemm/build-cudss-session/_deps/tangle_upstream-subbuild/tangle_upstream-populate-prefix/src/tangle_upstream-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/home/alex/xfemm_cuda/xfemm/build-cudss-session/_deps/tangle_upstream-subbuild/tangle_upstream-populate-prefix/src/tangle_upstream-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
