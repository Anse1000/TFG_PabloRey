# Nombre del ejecutable
EXEC = TFG_PRA

# Archivos fuente
SRC = main.c file_handler.c simulation.c

# Archivos de encabezado
HEADERS = file_handler.h simulation.h types.h

# Objetos intermedios
OBJ = main.o file_handler.o simulation.o

# Compilador
CC = gcc

DIR = reducidos

# Opciones de compilación generales
CFLAGS = -O3 -static -lm -Wall -Wextra

# Opciones de AVX-512 (solo para simulation.c)
AVX512_FLAGS = -march=cascadelake -mavx512f -mavx512dq -DAVX_512

# Regla por defecto (compilar)
all: $(EXEC)

# Compilar el ejecutable
$(EXEC): $(OBJ)
	$(CC) $^ -o $@ $(CFLAGS)

# Compilar archivos normales (.c en .o)
%.o: %.c $(HEADERS)
	$(CC) -c $< -o $@ $(CFLAGS)

# Compilar simulation.c con AVX-512 si se especifica
simulation.o: simulation.c simulation.h types.h
ifdef AVX_512
	$(CC) $(AVX512_FLAGS) -c $< -o $@ $(CFLAGS)
else
	$(CC) -c $< -o $@ $(CFLAGS)
endif

# Ejecutar el programa
run: all
	./$(EXEC) $(DIR)

# Limpiar archivos generados
clean:
	rm -f $(EXEC) $(OBJ)

