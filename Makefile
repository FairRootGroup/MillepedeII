# #################################################################
# Makefile for MillePede II with possible input from C
# (works on 32-bit SLC3 and and 64-bit SLC4)
# #################################################################
#
# compiler options
#
# To make it link with a BIG static array in dynal.inc,
# but does not work on 32 bit machines:
#LARGE_MEMORY_OPT=-mcmodel=medium
LARGE_MEMORY_OPT=
# All but 'yes' disables support of reading C-binaries:
SUPPORT_READ_C = yes
#
FCOMP = gcc
F_FLAGS = -Wall -fno-automatic -fno-backslash -O3 ${LARGE_MEMORY_OPT}
# Blobel:  -Wunused -fno-automatic -fno-backslash  -O3
#
CCOMP = gcc 
C_FLAGS = -Wall -O3 -Df2cFortran ${LARGE_MEMORY_OPT}
C_INCLUDEDIRS =  # e.g. -I .
DEBUG =          # e.g. -g
#
LOADER = gcc
L_FLAGS = -Wall -O3 ${LARGE_MEMORY_OPT}
#
# objects for this project
#
USER_OBJ_PEDE = pede.o mptest.o mille.o mpnum.o mptext.o mphistab.o \
	dynal.o minresblas.o minres.o vertpr.o linesrch.o
#
# Chose flags/object files for C-binary support:
#
ifeq ($(SUPPORT_READ_C),yes)
  F_FLAGS += -DREAD_C_FILES
  USER_OBJ_PEDE += readc.o
endif
#  
#
# Make the executables
#
pede : 	${USER_OBJ_PEDE} Makefile
	$(LOADER) $(L_FLAGS) -lg2c -lfrtbegin \
		-o $@ ${USER_OBJ_PEDE} 
#
#
#
clean:
	rm -f *.o *~
#
# Make the object files - depend on source and include file 
#
%.o : %.F Makefile
	${FCOMP} ${F_FLAGS} -c $< -o $@ 
#
# these two depend on the included memory:
dynal.o : dynal.F dynal.inc 
pede.o : pede.F dynal.inc 
#
#
%.o: %.c Makefile
	$(CCOMP) -c $(C_FLAGS) $(DEFINES) $(C_INCLUDEDIRS) $(DEBUG) -o $@ $<
#
# ##################################################################
# END
# ##################################################################
