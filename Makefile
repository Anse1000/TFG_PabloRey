# Nombre del ejecutable
EXEC = TFG_PRA

# Archivos fuente (con las rutas actualizadas)
SRC = src/main.c src/file_handler.c src/simulation.c

# Archivos de encabezado
HEADERS = src/file_handler.h src/simulation.h src/types.h

# Objetos intermedios (con las rutas correctas)
OBJ = main.o file_handler.o simulation.o

# Compiladores
CC = gcc
NVCC = nvcc

# Opciones de compilación generales
CFLAGS = -O3 -lm -Wall -Wextra

# Opciones de AVX-512
ifdef AVX_512
	AVX512_FLAGS = -march=cascadelake -mavx512f -mavx512dq -DAVX_512
endif

# Opciones de CUDA (solo si se define CUDA)
ifdef CUDA
    SRC += src/cuda_functions.cu
    OBJ += cuda_functions.o
    CFLAGS += -DCUDA
    CUDA_ARCH = -gencode arch=compute_75,code=sm_75
    CUDA_FLAGS = -O3 -Xcompiler "-Wall -Wextra" $(CUDA_ARCH)
    CUDA_LIBS =
endif

# Regla por defecto (compilar)
all: $(EXEC)

# Compilar archivos normales (.c en .o)
%.o: src/%.c $(HEADERS)
	$(CC) -c $< -o $@ $(CFLAGS)

# Compilar simulation.c con AVX-512 si se especifica
simulation.o: src/simulation.c src/simulation.h src/types.h
ifdef AVX_512
	$(CC) $(AVX512_FLAGS) -c $< -o $@ $(CFLAGS)
else
	$(CC) -c $< -o $@ $(CFLAGS)
endif

# Compilar CUDA solo si está definido
ifdef CUDA
cuda_functions.o: src/cuda_functions.cu src/cuda_functions.h
	$(NVCC) $(CUDA_FLAGS) -c $< -o $@
endif

# Compilar el ejecutable
$(EXEC): $(OBJ)
ifdef CUDA
	$(NVCC) $^ -o $@ $(CUDA_FLAGS) $(CUDA_LIBS)
else
	$(CC) $^ -o $@ $(CFLAGS)
endif

# Limpiar archivos generados
clean:
	rm -f $(EXEC) *.o




