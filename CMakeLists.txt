cmake_minimum_required(VERSION 3.18)

# Configurar compiladores Intel antes de definir el proyecto
set(CMAKE_C_COMPILER "icx")
set(CMAKE_CXX_COMPILER "icpx")

project(TFG_PRA C CXX)

option(ENABLE_AVX512 "Enable AVX-512 optimization" OFF)
option(ENABLE_CUDA "Enable CUDA support" OFF)

find_package(OpenMP REQUIRED)

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

# Opciones de compilación específicas para ICC
if(CMAKE_C_COMPILER_ID MATCHES "Intel")
    # Opciones comunes para ICC
    target_compile_options(TFG_PRA PRIVATE -O3 -Wall -xHost -D_FILE_OFFSET_BITS=64)
    # AVX-512 con sintaxis específica de Intel
    if(ENABLE_AVX512)
        target_compile_options(TFG_PRA PRIVATE -xCASCADELAKE -qopt-zmm-usage=high)
        target_compile_definitions(TFG_PRA PRIVATE AVX_512)
    endif()
else()
    # Opciones para compiladores no Intel (como respaldo)
    target_compile_options(TFG_PRA PRIVATE -O3 -Wall -Wextra)

    # AVX-512 para compiladores GCC/Clang
    if(ENABLE_AVX512)
        target_compile_options(TFG_PRA PRIVATE -march=cascadelake -mavx512f -mavx512dq)
        target_compile_definitions(TFG_PRA PRIVATE AVX_512)
    endif()
endif()

target_link_libraries(TFG_PRA PRIVATE m)  # Link a math

# Añadir flags de OpenMP
target_link_libraries(TFG_PRA PRIVATE OpenMP::OpenMP_C)

# Flags de depuración ajustadas para ICC
if(CMAKE_C_COMPILER_ID MATCHES "Intel")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g -O2 -fno-omit-frame-pointer")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O2 -traceback")
else()
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -fsanitize=thread -g3 -O1 -fno-omit-frame-pointer")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=thread -g3 -O1 -fno-omit-frame-pointer")
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

set_target_properties(TFG_PRA PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)