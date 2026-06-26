cmake_minimum_required(VERSION 3.15)
project(AstroFluid VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(OpenMP REQUIRED)

set(
    SOURCES
    main.cpp
    src/gas_cloud.cpp
    src/gas_cloud.hpp
    src/gravity.cpp
    src/gravity.hpp
    src/spatial_hash.cpp
    src/spatial_hash.hpp
    src/sph.cpp
    src/sph.hpp
    src/vec3.cpp
    src/vec3.hpp
    src/gas_material.hpp
    src/universe.cpp
    src/universe.hpp
)

add_executable(simulation ${SOURCES})

target_include_directories(simulation PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

target_link_libraries(simulation PRIVATE OpenMP::OpenMP_CXX)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|clang")
    target_compile_options(simulation PRIVATE -O3 -march=native)
endif()