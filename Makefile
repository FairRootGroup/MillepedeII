# #################################################################
# Makefile for MillePede II with possible input from C
# (works on 32-bit SLC4)
# #################################################################
#
# compiler options
#
# To make it link with a BIG (>2 GB) static array in dynal.inc,
# we need '-mcmodel=medium', but does not work on 32 bit machines,
# so here we switch off and restrict to smaller binaries
LARGE_MEMORY_OPT=
# All but 'yes' disables support of reading C-binaries:
SUPPORT_READ_C = yes
# If yes (and if SUPPORT_READ_C is yes), uses rfio, i.e. shift library, for IO,
# requires shift library to be installed
# (valid only for default executables, pede_*rfio*-ones will always get it):
SUPPORT_C_RFIO =
# yes
#
# Define the number of words for matrix/vector storage for the default pede program:
NUMBER_OF_WORDS = 100000000 # =400 MB
#
#
FCOMP = gcc
F_FLAGS = -Wall -fno-automatic -fno-backslash -O3 ${LARGE_MEMORY_OPT}
#
CCOMP = gcc 
C_FLAGS = -Wall -O3 -Df2cFortran ${LARGE_MEMORY_OPT}
C_INCLUDEDIRS =  # e.g. -I .
C_LIBS = -lg2c -lfrtbegin
DEBUG =          # e.g. -g
#
LOADER = gcc
L_FLAGS = -Wall -O3 ${LARGE_MEMORY_OPT}
#
# objects for this project
#
USER_OBJ_PEDE_STATIC = mptest.o mille.o mpnum.o mptext.o mphistab.o \
	minresblas.o minres.o vertpr.o linesrch.o
USER_OBJ_PEDE = ${USER_OBJ_PEDE_STATIC} pede.o dynal.o
#
# Chose flags/object files for C-binary support:
#
ifeq ($(SUPPORT_READ_C),yes)
  F_FLAGS += -DREAD_C_FILES
  USER_OBJ_PEDE += readc.o
  ifeq ($(SUPPORT_C_RFIO),yes)
    C_FLAGS += -DUSE_SHIFT_RFIO
    C_LIBS += -lshift
  endif
endif
#  
#
# Make the executables
# The specific ones (part of hack, see below) + the single:
EXECUTABLES = pede pede_1GB pede_1GB_rfio pede_2GB pede_2GB_rfio

all:	$(EXECUTABLES)

# The single special one:
pede : 	${USER_OBJ_PEDE} Makefile
	$(LOADER) $(L_FLAGS) $(C_LIBS) \
		-o $@ ${USER_OBJ_PEDE} 
#  
#
#
clean:
	rm -f *.o *~
#
clobber: clean 
	rm -f $(EXECUTABLES)

install: $(EXECUTABLES) #clean
	mkdir -p bin
	mv $(EXECUTABLES) bin

# Make the object files - depend on source and include file 
#
%.o : %.F Makefile
	${FCOMP} ${F_FLAGS} -c $< -o $@ 
# These two depend on the memory defined by NUMBER_OF_WORDS, 
# the second also on variable definitions:
dynal.o : dynal.F dynal.inc Makefile
	${FCOMP} ${F_FLAGS} -DNUMBER_OF_WORDS=${NUMBER_OF_WORDS} -c $< -o $@
pede.o : pede.F dynal.inc mpinds.inc Makefile
	${FCOMP} ${F_FLAGS} -DNUMBER_OF_WORDS=${NUMBER_OF_WORDS} -c $< -o $@
#
%.o: %.c Makefile
	$(CCOMP) -c $(C_FLAGS) $(DEFINES) $(C_INCLUDEDIRS) $(DEBUG) -o $@ $<
#
####################################################################
# Here we start the hack for the various executables...

%_1GB.o: %.F dynal.inc Makefile
	${FCOMP} ${F_FLAGS} -DNUMBER_OF_WORDS=250000000 -c $< -o $@
%_2GB.o: %.F dynal.inc Makefile
	${FCOMP} ${F_FLAGS} -DNUMBER_OF_WORDS=500000000 -c $< -o $@
%_rfio.o: %.c Makefile
	$(CCOMP) -c $(C_FLAGS) -DUSE_SHIFT_RFIO $(DEFINES) $(C_INCLUDEDIRS) \
		$(DEBUG) -o $@ $<
#
#
#
pede_1GB: ${USER_OBJ_PEDE_STATIC} pede_1GB.o dynal_1GB.o readc.o Makefile
	$(LOADER) $(L_FLAGS) $(C_LIBS) \
		-o $@ ${USER_OBJ_PEDE_STATIC} pede_1GB.o dynal_1GB.o readc.o  
#
pede_2GB: ${USER_OBJ_PEDE_STATIC} pede_2GB.o dynal_2GB.o readc.o Makefile
	$(LOADER) $(L_FLAGS) $(C_LIBS) \
		-o $@ ${USER_OBJ_PEDE_STATIC} pede_2GB.o dynal_2GB.o readc.o  
#
pede_1GB_rfio: ${USER_OBJ_PEDE_STATIC} pede_1GB.o dynal_1GB.o readc_rfio.o Makefile
	$(LOADER) $(L_FLAGS) $(C_LIBS) -lshift \
		-o $@ ${USER_OBJ_PEDE_STATIC} pede_1GB.o dynal_1GB.o readc_rfio.o  
#
pede_2GB_rfio: ${USER_OBJ_PEDE_STATIC} pede_2GB.o dynal_2GB.o readc_rfio.o Makefile
	$(LOADER) $(L_FLAGS) $(C_LIBS) -lshift \
		-o $@ ${USER_OBJ_PEDE_STATIC} pede_2GB.o dynal_2GB.o readc_rfio.o  
#
# End hack for the various executables...
####################################################################
#
# ##################################################################
# END
# ##################################################################
