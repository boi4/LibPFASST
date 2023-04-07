#
# Makefile for advection diffusion example
#

LIBPFASST ?= ../..
BUILDDIR = build

EXE = main.exe
USE_FFT=FALSE
include $(LIBPFASST)/Makefile.defaults

FSRC = src/level.f90 src/sweeper.f90 src/hooks.f90 src/probin.f90 src/main.f90

OBJ  = $(subst src, build,$(FSRC:.f90=.o))

FFLAGS  += -I$(LIBPFASST)/include -I $(MPIDYNRES_PATH)/build/include
LDFLAGS += -L $(MPIDYNRES_PATH)/build/lib -lmpidynres

all: $(EXE)

VPATHS = src 

include $(LIBPFASST)/Makefile.rules

#
# dependencies
#
include .depend
main.exe : $(LIBPFASST)/lib/libpfasst.a
