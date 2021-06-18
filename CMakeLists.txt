# project
cmake_minimum_required (VERSION 2.6)
project (libelas)

# use sse3 instruction set
SET(CMAKE_CXX_FLAGS "-msse3 -I/home/oalexan1/miniconda3/envs/isis4.4/include")

FIND_PACKAGE(TIFF)
INCLUDE_DIRECTORIES(${TIFF_INCLUDE_DIR})

# directories
set (LIBELAS_SRC_DIR src)

# include directory
include_directories("${LIBELAS_SRC_DIR}")

# sources
FILE(GLOB LIBELAS_SRC_FILES "src/*.cpp" "src/*.c")

# make release version
set(CMAKE_BUILD_TYPE Release)

link_libraries(${TIFF_LIBRARY_RELEASE})

# build demo program
add_executable(elas ${LIBELAS_SRC_FILES})

