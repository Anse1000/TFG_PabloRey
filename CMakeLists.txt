cmake_minimum_required(VERSION 3.18)

# Configurar compiladores Intel antes de definir el proyecto
find_program(INTEL_COMPILER icx)
if (INTEL_COMPILER)
    set(CMAKE_C_COMPILER "icx")
    set(CMAKE_CXX_COMPILER "icpx")
else ()
    set(CMAKE_C_COMPILER "gcc")
    set(CMAKE_CXX_COMPILER "g++")
endif ()

project(TFG_PRA C CXX)

option(ENABLE_CUDA "Enable CUDA support" OFF)

find_package(OpenMP REQUIRED)

# Archivos fuente base
set(SRC
        src/main.c
        src/file_handler.c
        src/simulation.c
        src/calculations.c
        src/aux_fun.c
        src/octree.c
)

set(HEADERS
        src/file_handler.h
        src/simulation.h
        src/types.h
        src/calculations.h
        src/aux_fun.h
        src/octree.h
)

# CUDA
if(ENABLE_CUDA)
    enable_language(CUDA)
    list(APPEND SRC src/cuda_functions.cu)
    message(STATUS "CUDA enabled, using ${CMAKE_CUDA_COMPILER}")
endif()

add_executable(TFG_PRA ${SRC})

# Opciones de compilación específicas
if(CMAKE_C_COMPILER_ID MATCHES "Intel")
    target_compile_options(TFG_PRA PRIVATE -O3 -Wall -xHost -D_FILE_OFFSET_BITS=64)
else()
    target_compile_options(TFG_PRA PRIVATE -O3 -Wall -Wextra)
endif()

# Enlazar con libm
find_library(MATH_LIBRARY m)
if(MATH_LIBRARY)
    target_link_libraries(TFG_PRA PRIVATE ${MATH_LIBRARY})
endif()

# Añadir flags de OpenMP
target_link_libraries(TFG_PRA PRIVATE OpenMP::OpenMP_C)

# Añadir OpenMP también para archivos CUDA
if(ENABLE_CUDA)
    target_link_libraries(TFG_PRA PRIVATE OpenMP::OpenMP_CXX)
endif()


# Flags de depuración
if(CMAKE_C_COMPILER_ID MATCHES "Intel")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g -O2 -fno-omit-frame-pointer")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O2 -traceback")
else()
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g3 -O1 -fno-omit-frame-pointer")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g3 -O1 -fno-omit-frame-pointer")
endif()

# CUDA opciones y propiedades
if(ENABLE_CUDA)
    target_compile_definitions(TFG_PRA PRIVATE CUDA)

    set_target_properties(TFG_PRA PROPERTIES
            CUDA_RUNTIME_LIBRARY Shared
            CUDA_SEPARABLE_COMPILATION ON
            CUDA_ARCHITECTURES 80
    )
    target_compile_options(TFG_PRA PRIVATE
            $<$<COMPILE_LANGUAGE:CUDA>:-O3>
            $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=-fopenmp>
            # Solo en Debug: info de líneas y verbose de ptxas para profiling
            $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:CUDA>>:-lineinfo -Xptxas=-v>
    )
endif()

# Salida del ejecutable
set_target_properties(TFG_PRA PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)
# === Definir macros según el tipo de build ===
target_compile_definitions(TFG_PRA PRIVATE
        $<$<CONFIG:Debug>:DEBUG_BUILD>
        $<$<CONFIG:Release>:RELEASE_BUILD>
)
