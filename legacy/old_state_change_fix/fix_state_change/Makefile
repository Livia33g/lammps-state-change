# Makefile for LAMMPS fix_state_change

# LAMMPS base directory (adjust this to your LAMMPS installation)
LAMMPS_DIR ?= /work/nvme/bewl/lguttieres/lammps_build/lammps

# Include LAMMPS Makefile
include $(LAMMPS_DIR)/src/MAKE/Makefile.mpi

# Source files
SRC = fix_state_change.cpp
OBJ = $(SRC:.cpp=.o)

# Compiler flags
CXXFLAGS += -O2 -g -fPIC

# Target
all: libfix_state_change.a

libfix_state_change.a: $(OBJ)
	$(AR) $(ARFLAGS) $@ $^
	@echo "Built $@"

clean:
	rm -f *.o libfix_state_change.a

# Installation instructions:
# 1. Copy fix_state_change.h to $(LAMMPS_DIR)/src
# 2. Copy fix_state_change.cpp to $(LAMMPS_DIR)/src
# 3. Add "fix_state_change" to the list in $(LAMMPS_DIR)/src/Makefile
# 4. Rebuild LAMMPS

