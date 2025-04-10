cmake_minimum_required(VERSION 3.28)
project(TFG_PRA C)
find_package(MPI REQUIRED)

set(CMAKE_C_STANDARD 11)

add_executable(TFG_PRA main.c )
target_link_libraries(TFG_PRA PUBLIC MPI::MPI_C m)
