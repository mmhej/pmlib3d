#---------------------------------------------------------------------------------#
# pmlib_client makefile
#---------------------------------------------------------------------------------#
# <<HEADER>>
#---------------------------------------------------------------------------------#
PROGRAM   = poisson

SRC_DIR  := ./sources
OBJ_DIR  := ./objects

DEFAULT: $(PROGRAM)

#---------------------------------------------------------------------------------#
# Path to libraries
#---------------------------------------------------------------------------------#
PMDIR      = ./../..
POISSONDIR = ./../../../poisson_solver/fortran

#---------------------------------------------------------------------------------#
# Flags and libraries for compiling
#---------------------------------------------------------------------------------#
FC      = mpif90
FFLAGS  = -O3 -ffree-form -ffree-line-length-none -ffree-line-length-none

LDFLAGS = -L$(PMDIR) -lpmlib -L$(POISSONDIR) -lpoissonsolver

INCL    = -I./objects -I$(PMDIR)/include -I$(POISSONDIR)/include

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
CPCMD = cpp $< $@
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
$(OBJ_DIR)/poisson.f:\
  $(PMDIR)/*.a\
  $(POISSONDIR)/*.a
$(OBJ_DIR)/poisson.o:\
  $(OBJ_DIR)/mod_output.o
#A
#B
#C
#D
#E
#F
#G
#H
#I
#J
#K
#L
#M
#N
#O
$(OBJ_DIR)/mod_output.f:\
  ./sources/output/*.f
$(OBJ_DIR)/mod_output.o:
#P
#Q
#R
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
	rm -f  $(LOG)
	rm -fr $(OBJ_DIR)
	rm -f  *.mod 2> /dev/null
	rm -f  $(PROGRAM)

