cmake_minimum_required(VERSION 3.10)
project(mppnccombine-fast)

# Compiler-specific options
add_compile_options(-Wall -Wextra)
if (${C_COMPILER_ID} MATCHES GNU)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
elseif (${C_COMPILER_ID} MATCHES Clang)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-missing-field-initializers")
elseif (${C_COMPILER_ID} MATCHES Intel)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
endif()

# Locate libraries
if (DEFINED ENV{OPENMPI_ROOT})
        set(MPI_C_ADDITIONAL_INCLUDE_DIRS $ENV{OPENMPI_ROOT}/include)
        set(MPI_CXX_ADDITIONAL_INCLUDE_DIRS $ENV{OPENMPI_ROOT}/include)
endif()
include(FindMPI)

find_path(HDF5_INCLUDE_DIRECTORY
        NAMES hdf5.h
        HINTS $ENV{HDF5_ROOT}/include
        )
find_library(HDF5_LIB
        NAMES hdf5
        HINTS $ENV{HDF5_ROOT}/lib
        )
find_library(HDF5_HL_LIB
        NAMES hdf5_hl
        HINTS $ENV{HDF5_ROOT}/lib
        )
find_path(NETCDF_INCLUDE_DIRECTORY
        NAMES netcdf.h
        HINTS $ENV{NETCDF_ROOT}/include
        )
find_library(NETCDF_LIB
        NAMES netcdf
        HINTS $ENV{NETCDF_ROOT}/lib
        )

add_executable(mppnccombine-fast mppnccombine-fast.c async.c error.c read_chunked.c)
set_target_properties(mppnccombine-fast PROPERTIES
        C_STANDARD 11)
target_link_libraries(mppnccombine-fast PUBLIC
        m
        MPI::MPI_C
        ${HDF5_LIB}
        ${HDF5_HL_LIB}
        ${NETCDF_LIB})
target_include_directories(mppnccombine-fast PRIVATE
        ${HDF5_INCLUDE_DIRECTORY}
        ${NETCDF_INCLUDE_DIRECTORY})

# Enable clang-tidy
find_program(CLANG_TIDY
        NAMES clang-tidy
        )
if (CLANG_TIDY)
        set_target_properties(mppnccombine-fast PROPERTIES
                C_CLANG_TIDY "${CLANG_TIDY};-checks=*,-hicpp-signed-bitwise,-clang-diagnostic-missing-field-initializers")
endif()

# Install
install(TARGETS mppnccombine-fast
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        INCLUDES DESTINATION include
        )
