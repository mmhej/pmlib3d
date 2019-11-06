#---------------------------------------------------------------------------------#
# pmlib3d makefile
#---------------------------------------------------------------------------------#
TARGET = libpmlib.a

default: $(TARGET)

SRC_DIR := ./source
OBJ_DIR := ./objects
INC_DIR := ./include

#---------------------------------------------------------------------------------#
# Path to libraries
#---------------------------------------------------------------------------------#
MKL_DIR  = /home/mek/mmhej/libs/intel-mkl-2013.1.117
FFTW_DIR = /home/mmh/Programs/fftw-3.3.5_gcc-4.8
#FFTW_DIR = /home/opt/el6/sl230s/fftw-3.3.3-sl230s-tm-intel-2013.1.117-openmpi-1.6.3-1

#---------------------------------------------------------------------------------#
# Flags and libraries for compiling
#---------------------------------------------------------------------------------#
# Fortran compiler
FC      = mpif90

# Fortran compiler flags
#FFLAGS = -O1 -free
FFLAGS  = -O3 -ffree-form -ffree-line-length-none

# Library flags
LDFLAGS = 
#LDFLAGS = -L$(MKL_DIR)/lib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
#LDFLAGS = -L$(FFTW_DIR)/lib -lfftw3

# Include paths
INCL    = -I$(INC_DIR) 
#INCL    = -I$(INC_DIR) -I$(MKL_DIR)/include
#INCL    = -I$(INC_DIR) -I$(FFTW_DIR)/include

# Preprocessor compiler
CPP     = cpp

# Preprocessor compiler flags
CPPFLAGS = -D__double
#CPPFLAGS = -D__double -D__mkl
#CPPFLAGS = -D__double -D__fftw

#---------------------------------------------------------------------------------#
# write scripts
#---------------------------------------------------------------------------------#
SHELL := /bin/sh
LOG   := compile.log
RUN   := $(SHELL) ./make_output.sh

#---------------------------------------------------------------------------------#
# Setup build directories
#---------------------------------------------------------------------------------#
NAMES   := $(notdir $(wildcard $(SRC_DIR)/*.f))
SOURCES := $(NAMES:%.f=$(SRC_DIR)/%.f)
OBJECTS := $(NAMES:%.f=$(OBJ_DIR)/%.o)

$(warning Setup directories...)
$(shell test -d $(OBJ_DIR) || mkdir $(OBJ_DIR))
$(shell test -d $(INC_DIR) || mkdir $(INC_DIR))
$(shell test -e $(LOG) && rm $(LOG))

# Dont delete the given intermediate files
.SECONDARY:

#---------------------------------------------------------------------------------#
# Preproces
#---------------------------------------------------------------------------------#
CPCMD = $(CPP) $(CPPFLAGS) $< $@
$(OBJ_DIR)/%.f : $(SRC_DIR)/%.f
	@printf "  CPP  %-42s" $<; \
	$(RUN) "$(CPCMD)" $(LOG) "Preprocessing Error"

#---------------------------------------------------------------------------------#
# Compile
#---------------------------------------------------------------------------------#
COMPILECMD = $(FC) $(FFLAGS) $(LDFLAGS) $(INCL) -c -o $@ $<
$(OBJ_DIR)/%.o : $(OBJ_DIR)/%.f
	@printf "  FC   %-42s" $<; \
	$(RUN) "$(COMPILECMD)" $(LOG) "Compile Error"

#---------------------------------------------------------------------------------#
# Build library
#---------------------------------------------------------------------------------#
ARCMD = ar crs $@ $(OBJECTS)

$(TARGET): $(OBJECTS)
	@printf "  AR   %-42s" "Creating library archive"; \
	$(RUN) "$(ARCMD)" $(LOG) "Error Creating Archive"; \
	$(shell mv *.mod $(INC_DIR) 2> $(LOG))

#---------------------------------------------------------------------------------#
# Explicit dependencies
#---------------------------------------------------------------------------------#
#A
#B
#C
$(OBJ_DIR)/pmlib_communication.f:\
  $(SRC_DIR)/communication/*.f
$(OBJ_DIR)/pmlib_communication.o:\
  $(OBJ_DIR)/pmlib_topology.o\
  $(OBJ_DIR)/pmlib_parameters.o\
  $(OBJ_DIR)/pmlib_write.o
#D
#E
#F
$(OBJ_DIR)/pmlib_finalise.o:\
  $(OBJ_DIR)/pmlib_parameters.o\
  $(OBJ_DIR)/pmlib_write.o\
  $(OBJ_DIR)/pmlib_particles.o\
  $(OBJ_DIR)/pmlib_mesh.o\
  $(OBJ_DIR)/pmlib_topology.o\
  $(OBJ_DIR)/pmlib_poisson.o
$(OBJ_DIR)/pmlib_fourier.f:\
  $(SRC_DIR)/fourier/*.f
$(OBJ_DIR)/pmlib_fourier.o:\
  $(OBJ_DIR)/pmlib_parameters.o\
  $(OBJ_DIR)/pmlib_write.o\
  $(OBJ_DIR)/pmlib_mesh.o\
  $(OBJ_DIR)/pmlib_topology.o\
  $(OBJ_DIR)/pmlib_communication.o
#G
#H
#I
$(OBJ_DIR)/pmlib_interpolation.f:\
  $(SRC_DIR)/interpolation/*.f
$(OBJ_DIR)/pmlib_interpolation.o:\
  $(OBJ_DIR)/pmlib_particles.o\
  $(OBJ_DIR)/pmlib_mesh.o\
  $(OBJ_DIR)/pmlib_parameters.o\
  $(OBJ_DIR)/pmlib_communication.o\
  $(OBJ_DIR)/pmlib_topology.o
#J
#K
#L
#M
$(OBJ_DIR)/pmlib_mesh.f:\
  $(SRC_DIR)/mesh/*.f
$(OBJ_DIR)/pmlib_mesh.o:\
  $(OBJ_DIR)/pmlib_parameters.o\
  $(OBJ_DIR)/pmlib_communication.o\
  $(OBJ_DIR)/pmlib_topology.o
#N
#O
$(OBJ_DIR)/pmlib_output.f:\
  $(SRC_DIR)/output/*.f
$(OBJ_DIR)/pmlib_output.o:\
  $(OBJ_DIR)/pmlib_parameters.o\
  $(OBJ_DIR)/pmlib_particles.o\
  $(OBJ_DIR)/pmlib_mesh.o\
  $(OBJ_DIR)/pmlib_patch.o\
  $(OBJ_DIR)/pmlib_topology.o\
  $(OBJ_DIR)/pmlib_write.o
#P
$(OBJ_DIR)/pmlib_parameters.f:\
  $(SRC_DIR)/parameters/*.f
$(OBJ_DIR)/pmlib_parameters.o:\
  $(OBJ_DIR)/pmlib_write.o
$(OBJ_DIR)/pmlib_particles.f:\
  $(SRC_DIR)/particles/*.f
$(OBJ_DIR)/pmlib_particles.o:\
  $(OBJ_DIR)/pmlib_topology.o\
  $(OBJ_DIR)/pmlib_patch.o\
  $(OBJ_DIR)/pmlib_communication.o\
  $(OBJ_DIR)/pmlib_parameters.o\
  $(OBJ_DIR)/pmlib_write.o
$(OBJ_DIR)/pmlib_patch.f:\
  $(SRC_DIR)/patch/*.f
$(OBJ_DIR)/pmlib_patch.o:\
  $(OBJ_DIR)/pmlib_parameters.o\
  $(OBJ_DIR)/pmlib_write.o
$(OBJ_DIR)/pmlib_poisson.f:\
  $(SRC_DIR)/poisson/*.f
$(OBJ_DIR)/pmlib_poisson.o:\
  $(OBJ_DIR)/pmlib_parameters.o\
  $(OBJ_DIR)/pmlib_mesh.o\
  $(OBJ_DIR)/pmlib_topology.o\
  $(OBJ_DIR)/pmlib_patch.o\
  $(OBJ_DIR)/pmlib_communication.o\
  $(OBJ_DIR)/pmlib_interpolation.o\
  $(OBJ_DIR)/pmlib_fourier.o\
  $(OBJ_DIR)/pmlib_output.o\
  $(OBJ_DIR)/pmlib_write.o
#Q
#R
$(OBJ_DIR)/pmlib_regularise.f:\
  $(SRC_DIR)/regularise/*.f
$(OBJ_DIR)/pmlib_regularise.o:\
  $(OBJ_DIR)/pmlib_parameters.o\
  $(OBJ_DIR)/pmlib_mesh.o\
  $(OBJ_DIR)/pmlib_topology.o\
  $(OBJ_DIR)/pmlib_patch.o\
  $(OBJ_DIR)/pmlib_communication.o\
  $(OBJ_DIR)/pmlib_fourier.o\
  $(OBJ_DIR)/pmlib_output.o\
  $(OBJ_DIR)/pmlib_write.o
$(OBJ_DIR)/pmlib_remesh.f:\
  $(SRC_DIR)/remesh/*.f
$(OBJ_DIR)/pmlib_remesh.o:\
  $(OBJ_DIR)/pmlib_parameters.o\
  $(OBJ_DIR)/pmlib_mesh.o\
  $(OBJ_DIR)/pmlib_particles.o\
  $(OBJ_DIR)/pmlib_patch.o\
  $(OBJ_DIR)/pmlib_topology.o\
  $(OBJ_DIR)/pmlib_write.o
$(OBJ_DIR)/pmlib_repatch.f:\
  $(SRC_DIR)/repatch/*.f
$(OBJ_DIR)/pmlib_repatch.o:\
  $(OBJ_DIR)/pmlib_parameters.o\
  $(OBJ_DIR)/pmlib_topology.o\
  $(OBJ_DIR)/pmlib_patch.o\
  $(OBJ_DIR)/pmlib_particles.o\
  $(OBJ_DIR)/pmlib_mesh.o\
  $(OBJ_DIR)/pmlib_communication.o\
  $(OBJ_DIR)/pmlib_poisson.o\
  $(OBJ_DIR)/pmlib_write.o
#S
#T
$(OBJ_DIR)/pmlib_topology.f:\
  $(SRC_DIR)/topology/*.f
$(OBJ_DIR)/pmlib_topology.o:\
  $(OBJ_DIR)/pmlib_parameters.o\
  $(OBJ_DIR)/pmlib_patch.o\
  $(OBJ_DIR)/pmlib_write.o
#U
#V
$(OBJ_DIR)/pmlib_visualise.f:\
  $(SRC_DIR)/visualise/*.f
$(OBJ_DIR)/pmlib_visualise.o:\
  $(OBJ_DIR)/pmlib_parameters.o\
  $(OBJ_DIR)/pmlib_write.o\
  $(OBJ_DIR)/pmlib_patch.o\
  $(OBJ_DIR)/pmlib_topology.o
#W
$(OBJ_DIR)/pmlib_write.o:
#X
#Y
#Z


#---------------------------------------------------------------------------------#
# Clean
#---------------------------------------------------------------------------------#
clean:
	rm -fr $(OBJ_DIR)
	rm -fr $(INC_DIR)
	rm -f  *.mod 2> $(LOG)
	rm -f  $(TARGET) 
	rm -f  $(LOG)
