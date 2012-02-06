# #################################################################
# Makefile for MillePede II with possible input from C
# (works on 64-bit SL5)
# #################################################################
# On 32-bit systems:
#  LARGE_MEMORY_OPT = 
#  LARGE_SIZE=4
# #################################################################
# If compiler doesn't understands INTEGER(KIND=) (FORTRAN95) use
# LARGE_SIZE=4 and replace INTEGER(KIND=) by INTEGER*4 or
# LARGE_SIZE=8 and replace INTEGER(KIND=) by INTEGER*8
# e.g. change source files with sed:
# for name in pede.F dynal.F dynal.inc mpinds.inc; do sed 's/(KIND=LARGE)/*4/' $name >! tmp; mv tmp $name; done
# #################################################################
#
# compiler options
#
# To make it link with a BIG (>2 GB) static array in dynal.inc,
# we need '-mcmodel=medium', but does not work on 32 bit machines,
LARGE_MEMORY_OPT=-mcmodel=medium
# All but 'yes' disables support of reading C-binaries:
SUPPORT_READ_C = yes
# If yes (and if SUPPORT_READ_C is yes), uses rfio, i.e. shift library, for IO,
# requires shift library to be installed
# (valid only for default executables, pede_*rfio*-ones will always get it):
SUPPORT_C_RFIO =
# yes
# If yes (and if SUPPORT_READ_C is yes), use zlib to read gzipped binary files
SUPPORT_ZLIB = yes
#
# Define the number of words for matrix/vector storage for the default pede program:
NUMBER_OF_WORDS = 100000000 # =400 MB
# Define the size of LARGE integers (4: INTERGER*4, 8: INTEGER*8)
LARGE_SIZE=8
# ompP profiler (http://www.ompp-tool.com)
OMPP = 
#kinst-ompp
#
FCOMP = $(OMPP) gcc
F_FLAGS = -Wall -fautomatic -fno-backslash -O3 ${LARGE_MEMORY_OPT}
#
CCOMP = $(OMPP) gcc 
C_FLAGS = -Wall -O3 -Df2cFortran ${LARGE_MEMORY_OPT}
C_INCLUDEDIRS =  # e.g. -I .
# gcc3: C_LIBS = -lg2c -lfrtbegin
# gcc4: 
C_LIBS = -lgfortran -lgfortranbegin
DEBUG =          # e.g. -g
#
# Multithreading with OpenMP (TM)
C_LIBS  += -lgomp
F_FLAGS += -fopenmp
#
LOADER = $(OMPP) gcc
L_FLAGS = -Wall -O3 ${LARGE_MEMORY_OPT}
#
# objects for this project
#
USER_OBJ_PEDE_STATIC = mptest.o mptst2.o mille.o mpnum.o mptext.o mphistab.o \
	minresblas.o minres.o randoms.o vertpr.o linesrch.o bandmatrix/Dbandmatrix.o
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
  else
    ifeq ($(SUPPORT_ZLIB),yes)
      C_FLAGS += -DUSE_ZLIB
      C_LIBS += -lz
    endif
  endif
endif
#
F_FLAGS += -DLARGE_SIZE=${LARGE_SIZE}
#  
#
# Make the executables
# The specific ones (part of hack, see below) + the single:
EXECUTABLES = pede pede_1GB pede_2GB 
# Or also those supporting rfio-reading:
#EXECUTABLES = pede pede_1GB pede_2GB pede_1GB_rfio pede_2GB_rfio 
# If you need also the executables with more memory, you need a
# 64-bit environment and the '-mcmodel=medium' option (see above).
#	pede_4GB pede_4GB_rfio pede_8GB pede_8GB_rfio
#EXECUTABLES = pede pede_1GB pede_2GB pede_1GB_rfio pede_2GB_rfio \
#	pede_4GB pede_4GB_rfio pede_8GB pede_8GB_rfio

all:	$(EXECUTABLES)

# The single special one:
pede : 	${USER_OBJ_PEDE} Makefile
	$(LOADER) $(L_FLAGS) \
		-o $@ ${USER_OBJ_PEDE} $(C_LIBS)
#  
#
#
clean:
	rm -f *.o *~ */*.o */*~
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
# These two depend on the memory defined by NUMBER_OF_WORDS, LARGE_SIZE
# the second also on variable definitions:
dynal.o : dynal.F dynal.inc largeint.inc Makefile
	${FCOMP} ${F_FLAGS} -DNUMBER_OF_WORDS=${NUMBER_OF_WORDS} -c $< -o $@
pede.o : pede.F dynal.inc largeint.inc mpinds.inc localfit.inc Makefile
	${FCOMP} ${F_FLAGS} -DNUMBER_OF_WORDS=${NUMBER_OF_WORDS} -c $< -o $@
#
%.o: %.c Makefile
	$(CCOMP) -c $(C_FLAGS) $(DEFINES) $(C_INCLUDEDIRS) $(DEBUG) -o $@ $<
#
####################################################################
# Here we start the hack for the various executables...

%_1GB.o: %.F dynal.inc largeint.inc Makefile
	${FCOMP} ${F_FLAGS} -DNUMBER_OF_WORDS=250000000 -c $< -o $@
%_2GB.o: %.F dynal.inc largeint.inc Makefile
	${FCOMP} ${F_FLAGS} -DNUMBER_OF_WORDS=500000000 -c $< -o $@
%_4GB.o: %.F dynal.inc largeint.inc Makefile
	${FCOMP} ${F_FLAGS}  -DNUMBER_OF_WORDS=1000000000 -c $< -o $@
%_8GB.o: %.F dynal.inc largeint.inc Makefile
	${FCOMP} ${F_FLAGS}  -DNUMBER_OF_WORDS=2000000000 -c $< -o $@
%_16GB.o: %.F dynal.inc largeint.inc Makefile
	${FCOMP} ${F_FLAGS}  -DNUMBER_OF_WORDS=4000000000_8 -c $< -o $@
%_24GB.o: %.F dynal.inc largeint.inc Makefile
	${FCOMP} ${F_FLAGS}  -DNUMBER_OF_WORDS=6000000000_8 -c $< -o $@
%_32GB.o: %.F dynal.inc largeint.inc Makefile
	${FCOMP} ${F_FLAGS}  -DNUMBER_OF_WORDS=8000000000_8 -c $< -o $@
%_48GB.o: %.F dynal.inc largeint.inc Makefile
	${FCOMP} ${F_FLAGS}  -DNUMBER_OF_WORDS=12000000000_8 -c $< -o $@
%_64GB.o: %.F dynal.inc largeint.inc Makefile
	${FCOMP} ${F_FLAGS}  -DNUMBER_OF_WORDS=16000000000_8 -c $< -o $@
%_96GB.o: %.F dynal.inc largeint.inc Makefile
	${FCOMP} ${F_FLAGS}  -DNUMBER_OF_WORDS=24000000000_8 -c $< -o $@
%_rfio.o: %.c Makefile
	$(CCOMP) -c $(C_FLAGS) -DUSE_SHIFT_RFIO $(DEFINES) $(C_INCLUDEDIRS) \
		$(DEBUG) -o $@ $<
#
#
#
pede_1GB: ${USER_OBJ_PEDE_STATIC} pede_1GB.o dynal_1GB.o readc.o Makefile
	$(LOADER) $(L_FLAGS) \
		-o $@ ${USER_OBJ_PEDE_STATIC} pede_1GB.o dynal_1GB.o readc.o $(C_LIBS)
#
pede_2GB: ${USER_OBJ_PEDE_STATIC} pede_2GB.o dynal_2GB.o readc.o Makefile
	$(LOADER) $(L_FLAGS) \
		-o $@ ${USER_OBJ_PEDE_STATIC} pede_2GB.o dynal_2GB.o readc.o $(C_LIBS)
#
pede_1GB_rfio: ${USER_OBJ_PEDE_STATIC} pede_1GB.o dynal_1GB.o readc_rfio.o Makefile
	$(LOADER) $(L_FLAGS) \
		-o $@ ${USER_OBJ_PEDE_STATIC} pede_1GB.o dynal_1GB.o readc_rfio.o $(C_LIBS) -lshift
#
pede_2GB_rfio: ${USER_OBJ_PEDE_STATIC} pede_2GB.o dynal_2GB.o readc_rfio.o Makefile
	$(LOADER) $(L_FLAGS) \
		-o $@ ${USER_OBJ_PEDE_STATIC} pede_2GB.o dynal_2GB.o readc_rfio.o $(C_LIBS) -lshift
#
pede_4GB: ${USER_OBJ_PEDE_STATIC} pede_4GB.o dynal_4GB.o readc.o Makefile
	$(LOADER) $(L_FLAGS)  \
		-o $@ ${USER_OBJ_PEDE_STATIC} pede_4GB.o dynal_4GB.o readc.o $(C_LIBS)
#
pede_4GB_rfio: ${USER_OBJ_PEDE_STATIC} pede_4GB.o dynal_4GB.o readc_rfio.o  Makefile
	$(LOADER) $(L_FLAGS) \
		-o $@ ${USER_OBJ_PEDE_STATIC} pede_4GB.o dynal_4GB.o readc_rfio.o $(C_LIBS) -lshift
#
pede_8GB: ${USER_OBJ_PEDE_STATIC} pede_8GB.o dynal_8GB.o readc.o Makefile
	$(LOADER) $(L_FLAGS) \
		-o $@ ${USER_OBJ_PEDE_STATIC} pede_8GB.o dynal_8GB.o readc.o $(C_LIBS)
#
pede_8GB_rfio: ${USER_OBJ_PEDE_STATIC} pede_8GB.o dynal_8GB.o readc_rfio.o  Makefile
	$(LOADER) $(L_FLAGS) \
		-o $@ ${USER_OBJ_PEDE_STATIC} pede_8GB.o dynal_8GB.o readc_rfio.o $(C_LIBS) -lshift
#
pede_16GB: ${USER_OBJ_PEDE_STATIC} pede_16GB.o dynal_16GB.o readc.o Makefile
	$(LOADER) $(L_FLAGS) \
		-o $@ ${USER_OBJ_PEDE_STATIC} pede_16GB.o dynal_16GB.o readc.o $(C_LIBS)
#
pede_24GB: ${USER_OBJ_PEDE_STATIC} pede_24GB.o dynal_24GB.o readc.o Makefile
	$(LOADER) $(L_FLAGS) \
		-o $@ ${USER_OBJ_PEDE_STATIC} pede_24GB.o dynal_24GB.o readc.o $(C_LIBS)
pede_32GB: ${USER_OBJ_PEDE_STATIC} pede_32GB.o dynal_32GB.o readc.o Makefile
	$(LOADER) $(L_FLAGS) \
		-o $@ ${USER_OBJ_PEDE_STATIC} pede_32GB.o dynal_32GB.o readc.o $(C_LIBS)
pede_48GB: ${USER_OBJ_PEDE_STATIC} pede_48GB.o dynal_48GB.o readc.o Makefile
	$(LOADER) $(L_FLAGS) \
		-o $@ ${USER_OBJ_PEDE_STATIC} pede_48GB.o dynal_48GB.o readc.o $(C_LIBS)
pede_64GB: ${USER_OBJ_PEDE_STATIC} pede_64GB.o dynal_64GB.o readc.o Makefile
	$(LOADER) $(L_FLAGS) \
		-o $@ ${USER_OBJ_PEDE_STATIC} pede_64GB.o dynal_64GB.o readc.o $(C_LIBS)
pede_96GB: ${USER_OBJ_PEDE_STATIC} pede_96GB.o dynal_96GB.o readc.o Makefile
	$(LOADER) $(L_FLAGS) \
		-o $@ ${USER_OBJ_PEDE_STATIC} pede_96GB.o dynal_96GB.o readc.o $(C_LIBS)
#
# End hack for the various executables...
####################################################################
#
# ##################################################################
# END
# ##################################################################
