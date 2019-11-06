#---------------------------------------------------------------------------
# pmlib_client makefile
#---------------------------------------------------------------------------
#                  __    __
#    ___  _ _     / /__ / /
#   | _ \/ ` \ _ / // // _ \
#   | ,_/_,_,_\  \_/\_/\___/
#   |_|                      Mads M. Hejlesen* and Jens H. Walther*^
#                            * Technical University of Denmark
#                            ^ ETH-Zurich
#                            mmhej@mek.dtu.dk
#
#---------------------------------------------------------------------------
PROGRAM   = validate_mesh_mapping

SRC_DIR  := ./sources
OBJ_DIR  := ./objects

DEFAULT: $(PROGRAM)

#---------------------------------------------------------------------------
# Path to libraries
#---------------------------------------------------------------------------
PMDIR    = ./../../lib

#PMDIR    = /home/mmhej/Development/pmlib/lib3D
#PMDIR   = /home/mek/mmhej/pmlib/lib3D

#---------------------------------------------------------------------------
# Flags and libraries for compiling
#---------------------------------------------------------------------------
FC      = mpif90
FFLAGS  = -O3 -ffree-form -ffree-line-length-none -ffree-line-length-none

LDFLAGS =  -llapack -lblas -L$(PMDIR) -lpmlib

INCL    = -I./objects -I$(PMDIR)/include

#---------------------------------------------------------------------------
# write scripts
#---------------------------------------------------------------------------
SHELL := /bin/sh
LOG   := compile.log
RUN   := $(SHELL) ./make_output.sh

#---------------------------------------------------------------------------
# Setup build directories
#---------------------------------------------------------------------------
NAMES   := $(notdir $(wildcard $(SRC_DIR)/*.f))
SOURCES := $(NAMES:%.f=$(SRC_DIR)/%.f)
OBJECTS := $(NAMES:%.f=$(OBJ_DIR)/%.o)

$(warning Setup directories...)
$(shell test -d $(OBJ_DIR) || mkdir $(OBJ_DIR))
$(shell test -e $(LOG) && rm $(LOG))

# Dont delete the given intermediate files
.SECONDARY:

#---------------------------------------------------------------------------
# Preprocess
#---------------------------------------------------------------------------
CPCMD = cpp $< $@
$(OBJ_DIR)/%.f : $(SRC_DIR)/%.f
	@printf "  CPP  %-42s" $<; \
	$(RUN) "$(CPCMD)" $(LOG) "Preprocessing Error"

#---------------------------------------------------------------------------
# Compile
#---------------------------------------------------------------------------
COMPILECMD = $(FC) $(FFLAGS) $(LDFLAGS) $(INCL) -c -o $@ $<
$(OBJ_DIR)/%.o : $(OBJ_DIR)/%.f
	@printf "  FC   %-42s" $<; \
	$(RUN) "$(COMPILECMD)" $(LOG) "Compile Error"

#---------------------------------------------------------------------------
# Link
#---------------------------------------------------------------------------
LINKCMD = $(FC) $(OBJECTS) $(LDFLAGS) $(INCL) -o $@
$(PROGRAM): $(OBJECTS)
	@printf " LINK  %-42s" $@; \
	$(RUN) "$(LINKCMD)" $(LOG) "Linking Error"; \
	mv *.mod $(OBJ_DIR) 2> /dev/null

#---------------------------------------------------------------------------
# Explicit dependencies
#---------------------------------------------------------------------------
$(OBJ_DIR)/validate_mesh_mapping.o:
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

#---------------------------------------------------------------------------
# Clean
#---------------------------------------------------------------------------
clean:
	rm -f  $(LOG)
	rm -fr $(OBJ_DIR)
	rm -f  *.mod 2> /dev/null
	rm -f  $(PROGRAM)

