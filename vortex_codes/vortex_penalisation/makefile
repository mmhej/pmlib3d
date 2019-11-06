#---------------------------------------------------------------------------------#
# makefile
#---------------------------------------------------------------------------------#
PROGRAM   = vortex_penalisation

SRC_DIR  := ./source
OBJ_DIR  := ./objects

DEFAULT: $(PROGRAM)

#---------------------------------------------------------------------------------#
# Path to libraries
#---------------------------------------------------------------------------------#
PMLIB_DIR = ./../../
MKL_DIR   = /home/mek/mmhej/libs/intel-mkl-2013.1.117
FFTW_DIR  = /home/mmh/Programs/fftw-3.3.5_gcc-4.8
#FFTW_DIR  = /home/opt/el6/sl230s/fftw-3.3.3-sl230s-tm-intel-2013.1.117-openmpi-1.6.3-1

#---------------------------------------------------------------------------------#
# Flags and libraries for compiling
#---------------------------------------------------------------------------------#
# Fortran compiler
FC      = mpif90

# Fortran compiler flags
#FFLAGS  = -O1 -free
FFLAGS  = -O3 -ffree-form -ffree-line-length-none

# Library flags
LDFLAGS = -L$(PMLIB_DIR) -lpmlib
#LDFLAGS = -L$(PMLIB_DIR) -lpmlib \
#          -L$(MKL_DIR)/lib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
#LDFLAGS = -L$(PMLIB_DIR) -lpmlib -L$(FFTW_DIR)/lib -lfftw3

# Include paths
INCL    = -I./objects -I$(PMLIB_DIR)/include
#INCL    = -I./objects -I$(PMLIB_DIR)/include -I$(MKL_DIR)/include
#INCL    = -I./objects -I$(PMLIB_DIR)/include -I$(FFTW_DIR)/include

# Preprocessor compiler
CPP      = cpp

# Preprocessor compiler flags
CPPFLAGS = 

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
$(shell test -e $(LOG) && rm $(LOG))


# Dont delete the given intermediate files
.SECONDARY:

#---------------------------------------------------------------------------------#
# Preprocess
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
# Link
#---------------------------------------------------------------------------------#
LINKCMD = $(FC) $(OBJECTS) $(LDFLAGS) $(INCL) -o $@
$(PROGRAM): $(OBJECTS)
	@printf " LINK  %-42s" $@; \
	$(RUN) "$(LINKCMD)" $(LOG) "Linking Error"; \
	$(shell mv *.mod $(OBJ_DIR) 2> $(LOG))

#---------------------------------------------------------------------------------#
# Explicit dependencies
#---------------------------------------------------------------------------------#
$(OBJ_DIR)/main.f:\
  $(PMLIB_DIR)/*.a
$(OBJ_DIR)/main.o:\
  $(OBJ_DIR)/mod_parameters.o\
  $(OBJ_DIR)/mod_input.o\
  $(OBJ_DIR)/mod_output.o\
  $(OBJ_DIR)/mod_flowcases.o\
  $(OBJ_DIR)/mod_particle_derivatives.o\
  $(OBJ_DIR)/mod_penalisation.o\
  $(OBJ_DIR)/mod_diagnostics.o
#A
#B
#C
#D
$(OBJ_DIR)/mod_diagnostics.f:\
  $(SRC_DIR)/diagnostics/*.f
$(OBJ_DIR)/mod_diagnostics.o:\
  $(OBJ_DIR)/mod_parameters.o
#E
#F
$(OBJ_DIR)/mod_flowcases.f:\
  $(SRC_DIR)/flowcases/*.f
$(OBJ_DIR)/mod_flowcases.o:\
  $(OBJ_DIR)/mod_parameters.o\
  $(OBJ_DIR)/mod_penalisation.o
#G
#H
#I
$(OBJ_DIR)/mod_input.f:\
  $(SRC_DIR)/input/*.f
$(OBJ_DIR)/mod_input.o:\
  $(OBJ_DIR)/mod_parameters.o
#J
#K
#L
#M
#N
#O
$(OBJ_DIR)/mod_output.f:\
  $(SRC_DIR)/output/*.f
$(OBJ_DIR)/mod_output.o:\
  $(OBJ_DIR)/mod_parameters.o
#P
$(OBJ_DIR)/mod_parameters.o:
$(OBJ_DIR)/mod_penalisation.f:\
  $(SRC_DIR)/penalisation/*.f
$(OBJ_DIR)/mod_penalisation.o:\
  $(OBJ_DIR)/mod_parameters.o
#Q
#R
$(OBJ_DIR)/mod_particle_derivatives.f:\
  $(SRC_DIR)/particle_derivatives/*.f
$(OBJ_DIR)/mod_particle_derivatives.o:\
  $(OBJ_DIR)/mod_parameters.o
#S
#T
#U
#V
#W
#X
#Y
#Z

#---------------------------------------------------------------------------------#
# Clean
#---------------------------------------------------------------------------------#
clean:
	rm -fr $(OBJ_DIR)
	rm -f  *.mod 2> $(LOG)
	rm -f  $(PROGRAM)
	rm -f  $(LOG)

