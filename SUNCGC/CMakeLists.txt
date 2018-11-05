# Define library. Only source files here!
add_library(SUNCGC
  src/hilbertspace.cpp
  src/symbol3j.cpp
  src/vector.cpp
  )

# Define headers for this library. PUBLIC headers are used for
# compiling the library, and will be added to consumers' build
# paths.
target_include_directories(
  SUNCGC
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>
  PRIVATE src
  )

target_compile_features(SUNCGC PRIVATE cxx_std_17)

target_link_libraries(
  SUNCGC
  PUBLIC TNT
  PRIVATE openblas primme tcl hptt  tbb hdf5
)
#add_definitions(-DTBB_USE_DEBUG)

#tbbmalloc_debug tbbmalloc_proxy_debug
#target_link_libraries(
#  suncgc
#  PUBLIC Boost::system
#  PRIVATE openblas hdf5 Boost::filesystem
#)

# Every library has unit tests, of course
#add_executable(testlib
#    test/main.cpp
#    test/ising/main.cpp
#    )
#target_link_libraries(testlib
#    tns
#    )
#add_test(testlib testlib)
