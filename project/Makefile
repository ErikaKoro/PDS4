CC = gcc
MPICC = mpicc
CILK = /usr/local/OpenCilk-2.0.0-x86_64-Linux-Ubuntu-22.04/bin/clang

BUILD_DIR := ./make-build-debug-g
MPI_BUILD_DIR := ./make-build-debug-mpi
CILK_BUILD_DIR := ./make-build-debug-cilk-20
PROF_DIR := ./profiling
SRC_DIRS := ./src

# Colors
GREEN = \033[1;32m
RED = \033[1;31m
NC = \033[0m
BOLD = \033[1m

# Directories

DATA_SRC := $(shell find $(SRC_DIRS)/Data -name '*.c')
DATA_SRC := $(DATA_SRC:%=$(BUILD_DIR)/%.o)

SERIAL_SRC := $(shell find $(SRC_DIRS)/sequential -name '*.c')
SERIAL_SRC += $(shell find $(SRC_DIRS) -maxdepth 1 -name '*.c')
SERIAL_SRC := $(SERIAL_SRC:%=$(BUILD_DIR)/%.o)

MPI_SRC := $(shell find $(SRC_DIRS)/mpi -name '*.c')
MPI_SRC += $(shell find $(SRC_DIRS) -maxdepth 1 -name '*.c')
MPI_SRC := $(MPI_SRC:%=$(MPI_BUILD_DIR)/%.o)

CILK_SRC := $(shell find $(SRC_DIRS)/cilk -name '*.c')
CILK_SRC += $(shell find $(SRC_DIRS) -maxdepth 1 -name '*.c')
CILK_SRC := $(CILK_SRC:%=$(CILK_BUILD_DIR)/%.o)

KNN_SRC := $(shell find $(SRC_DIRS)/knn -name '*.c')
KNN_SRC += $(shell find $(SRC_DIRS) -maxdepth 1 -name '*.c')
KNN_SRC := $(KNN_SRC:%=$(BUILD_DIR)/%.o)

# Every folder in ./src will need to be passed to GCC so that it can find header files
INC_DIRS := $(shell find $(SRC_DIRS) -type d)

# Add a prefix to INC_DIRS. So moduleA would become -ImoduleA. GCC understands this -I flag
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

CC_FLAGS := $(INC_FLAGS) -O3 -g
MPI_FLAGS := $(INC_FLAGS) -O3 -g -lm
CILK_FLAGS := $(INC_FLAGS) -O3 -fopencilk # -fsanitize=cilk -Og -g

all: build_sequential build_mpi build_cilk

$(BUILD_DIR)/data.out: $(DATA_SRC)
	@echo
	@echo -e "        $(BOLD)Linking...$(NC)"
	@echo
	@$(CC) $(CC_FLAGS) -o $(BUILD_DIR)/data.out $(DATA_SRC)
	@echo -e "    $(GREEN)Build finished successfully!$(NC)"
	@echo

$(BUILD_DIR)/sequential.out: $(SERIAL_SRC)
	@echo
	@echo -e "        $(BOLD)Linking...$(NC)"
	@echo
	@$(CC) $(CC_FLAGS) -o $(BUILD_DIR)/sequential.out $(SERIAL_SRC)
	@echo -e "    $(GREEN)Build finished successfully!$(NC)"
	@echo

$(BUILD_DIR)/knn.out: $(KNN_SRC)
	@echo
	@echo -e "        $(BOLD)Linking...$(NC)"
	@echo
	@$(CC) $(CC_FLAGS) -o $(BUILD_DIR)/knn.out $(KNN_SRC)
	@echo -e "    $(GREEN)Build finished successfully!$(NC)"
	@echo


$(MPI_BUILD_DIR)/mpi.out: $(MPI_SRC)
	@echo
	@echo -e "        $(BOLD)Linking...$(NC)"
	@echo
	@$(MPICC) $(MPI_FLAGS) -o $(MPI_BUILD_DIR)/mpi.out $(MPI_SRC)
	@echo -e "    $(GREEN)Build finished successfully!$(NC)"
	@echo


$(CILK_BUILD_DIR)/cilk.out: $(CILK_SRC)
	@echo
	@echo -e "        $(BOLD)Linking...$(NC)"
	@echo
	@$(CILK) $(CILK_FLAGS) -o $(CILK_BUILD_DIR)/cilk.out $(CILK_SRC)
	@echo -e "    $(GREEN)Build finished successfully!$(NC)"
	@echo


$(BUILD_DIR)/%.c.o: %.c
	@mkdir -p $(dir $@)
	@echo -e "        $(BOLD)Compiling:$(NC) $(<)..."
	@$(CC) $(CC_FLAGS) -c $< -o $@


$(MPI_BUILD_DIR)/%.c.o: %.c
	@mkdir -p $(dir $@)
	@echo -e "        $(BOLD)Compiling:$(NC) $(<)..."
	@$(MPICC) $(MPI_FLAGS) -c $< -o $@


$(CILK_BUILD_DIR)/%.c.o: %.c
	@mkdir -p $(dir $@)
	@echo -e "        $(BOLD)Compiling:$(NC) $(<)..."
	@$(CILK) $(CILK_FLAGS) -c $< -o $@

%.c:

build_data: $(BUILD_DIR)/data.out
run_data: $(BUILD_DIR)/data.out
	@echo
	@echo
	@$(BUILD_DIR)/data.out 10 4194304 ./src/Data/data
	@echo
	@echo

build_sequential: $(BUILD_DIR)/sequential.out
run_sequential: $(BUILD_DIR)/sequential.out
	@echo
	@echo
	@$(BUILD_DIR)/sequential.out ./src/Data/data
	@echo
	@echo

build_knn: $(BUILD_DIR)/knn.out
run_knn: $(BUILD_DIR)/knn.out
	@echo
	@echo
	@$(BUILD_DIR)/knn.out ./src/Data/data 16
	@echo
	@echo

build_mpi: $(MPI_BUILD_DIR)/mpi.out
run_mpi: $(MPI_BUILD_DIR)/mpi.out
	@echo
	@echo
	@mpirun -hostfile hosts $(MPI_BUILD_DIR)/mpi.out ./src/Data/data 16
	@echo
	@echo


valgrind_mpi: $(MPI_BUILD_DIR)/mpi.out
	@echo
	@echo
	@mpirun -hostfile hosts valgrind --track-origins=yes --leak-check=full --suppressions=openmpi-valgrind.supp $(MPI_BUILD_DIR)/mpi.out ./src/Data/data 5
	@echo
	@echo

build_cilk: $(CILK_BUILD_DIR)/cilk.out
run_cilk: $(CILK_BUILD_DIR)/cilk.out
	@echo
	@echo
	@CILK_NWORKERS=16 $(CILK_BUILD_DIR)/cilk.out ./src/Data/data
	@echo
	@echo


.PHONY: clean
clean:
	@echo -e "$(RED)Clearing build directories...$(NC)"
	@rm -rf $(BUILD_DIR)
	@rm -rf $(CILK_BUILD_DIR)
	@rm -rf $(MPI_BUILD_DIR)
	@echo -e "$(GREEN)Done!$(NC)"