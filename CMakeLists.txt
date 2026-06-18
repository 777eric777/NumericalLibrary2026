cmake_minimum_required(VERSION 3.20)

# Force 64-bit generator platform when not explicitly provided by the build system
# (prevents mixing x86 system libraries with an x64 target on Windows/CMake+MSVC)
if(WIN32 AND NOT DEFINED CMAKE_GENERATOR_PLATFORM)
  set(CMAKE_GENERATOR_PLATFORM "x64" CACHE STRING "Platform for generator" FORCE)
endif()

project(NumericalLibrary LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(${CMAKE_SOURCE_DIR}/include)

# Library sources
file(GLOB LIB_SOURCES
	src/*.cpp
)

add_library(NumericalLibrary STATIC ${LIB_SOURCES})
target_include_directories(NumericalLibrary PUBLIC ${CMAKE_SOURCE_DIR}/include)

# Example executable
add_executable(usage_example examples/usage_example.cpp)
target_link_libraries(usage_example PRIVATE NumericalLibrary)

# Test runner
add_executable(tests_runner tests/test_runner.cpp)
target_link_libraries(tests_runner PRIVATE NumericalLibrary)

# Examples
add_executable(demo_root_finding examples/demo_root_finding.cpp)
target_link_libraries(demo_root_finding PRIVATE NumericalLibrary)

add_executable(demo_odes examples/demo_odes.cpp)
target_link_libraries(demo_odes PRIVATE NumericalLibrary)

# Install target (optional)
install(TARGETS NumericalLibrary DESTINATION lib)
install(DIRECTORY include/ DESTINATION include)
