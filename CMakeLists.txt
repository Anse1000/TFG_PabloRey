cmake_minimum_required(VERSION 3.18)
project(TFG_PRA C CXX)
find_package(MPI REQUIRED)
option(ENABLE_AVX512 "Enable AVX-512 optimization" OFF)
option(ENABLE_CUDA "Enable CUDA support" OFF)

# Archivos fuente base
set(SRC
        src/main.c
        src/file_handler.c
        src/simulation.c
)

set(HEADERS
        src/file_handler.h
        src/simulation.h
        src/types.h
)

# Si se activa CUDA, añadimos el .cu a los fuentes antes de crear el ejecutable
if(ENABLE_CUDA)
    enable_language(CUDA)
    list(APPEND SRC src/cuda_functions.cu)
endif()

# Crear el ejecutable
add_executable(TFG_PRA ${SRC} ${HEADERS})

# Opciones comunes
target_compile_options(TFG_PRA PRIVATE -O3 -Wall -Wextra)
target_link_libraries(TFG_PRA PRIVATE PUBLIC MPI::MPI_C m)  # Link a math

# AVX-512
if(ENABLE_AVX512)
    target_compile_options(TFG_PRA PRIVATE -march=cascadelake -mavx512f -mavx512dq)
    target_compile_definitions(TFG_PRA PRIVATE AVX_512)
endif()

# CUDA
if(ENABLE_CUDA)
    target_compile_definitions(TFG_PRA PRIVATE CUDA)
    set_target_properties(TFG_PRA PROPERTIES
            CUDA_SEPARABLE_COMPILATION ON
            CUDA_ARCHITECTURES 75
    )
    target_compile_options(TFG_PRA PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-O3>)
endif()

# (Opcional) Ejecutable en la carpeta build/
set_target_properties(TFG_PRA PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)




