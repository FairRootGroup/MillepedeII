# #################################################################
# Makefile for MillePede II (Fortran90) with possible input from C
# 
# Works on 64-bit SL5 with gcc version 4.4.4 (or higher).
# See comments about different gcc versions inline to get a
# hint about the necessary adjustments.
# #################################################################
# On 32-bit systems:
#  LARGE_SIZE=4
# #################################################################
#
# All but 'yes' disables support of reading C-binaries:
SUPPORT_READ_C = yes
# If yes (and if SUPPORT_READ_C is yes), uses rfio, i.e. shift library, for IO,
# requires shift library to be installed
SUPPORT_C_RFIO =
# yes
# If yes (and if SUPPORT_READ_C is yes), use zlib to read gzipped binary files
SUPPORT_ZLIB = yes
# Define the size of LARGE integers (4: INTERGER*4, 8: INTEGER*8)
LARGE_SIZE=8
# ompP profiler (http://www.ompp-tool.com)
OMPP = kinst-ompp
#
FCOMP = $(OMPP) gcc44
F_FLAGS = -Wall -fautomatic -fno-backslash -O3 -cpp -DLARGE_SIZE=${LARGE_SIZE}
#
CCOMP = $(OMPP) gcc44 
C_FLAGS = -Wall -O3 -Df2cFortran
C_INCLUDEDIRS =  # e.g. -I .
# gcc4 till gcc44: 
C_LIBS = -lgfortran -lgfortranbegin
# gcc45, gcc46: C_LIBS = -lgfortran -lm
DEBUG =          # e.g. -g
#
# Multithreading with OpenMP (TM)
C_LIBS  += -lgomp
F_FLAGS += -fopenmp
#
LOADER = $(OMPP) gcc44
L_FLAGS = -Wall -O3
#
# objects for this project
#
USER_OBJ_PEDE = mpdef.o mpdalc.o mpmod.o mpbits.o mptest1.o mptest2.o mille.o mpnum.o mptext.o mphistab.o \
	minresblas.o minres.o randoms.o vertpr.o linesrch.o Dbandmatrix.o pede.o
#
# Chose flags/object files for C-binary support:
#
ifeq ($(SUPPORT_READ_C),yes)
  F_FLAGS += -DREAD_C_FILES
  USER_OBJ_PEDE += readc.o
  ifeq ($(SUPPORT_C_RFIO),yes)
    C_FLAGS += -DUSE_SHIFT_RFIO
    C_LIBS += -lshift
  else
    ifeq ($(SUPPORT_ZLIB),yes)
      C_FLAGS += -DUSE_ZLIB
      C_LIBS += -lz
    endif
  endif
endif
#  
#
# Make the executables
EXECUTABLES = pede 
#

all:	$(EXECUTABLES)

pede : 	${USER_OBJ_PEDE} Makefile
	$(LOADER) $(L_FLAGS) \
		-o $@ ${USER_OBJ_PEDE} $(C_LIBS) 
#
clean:
	rm -f *.o *~ */*.o */*~ *.mod */*.mod
#
clobber: clean 
	rm -f $(EXECUTABLES)

install: $(EXECUTABLES) #clean
	mkdir -p bin
	mv $(EXECUTABLES) bin

# Make the object files - depend on source and include file 
#
%.o : %.f90 Makefile
	${FCOMP} ${F_FLAGS} -c $< -o $@ 
%.o: %.c Makefile
	$(CCOMP) -c $(C_FLAGS) $(DEFINES) $(C_INCLUDEDIRS) $(DEBUG) -o $@ $<
#
# ##################################################################
# END
# ##################################################################
