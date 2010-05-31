/* 
   C-methods to handle input of C/C++ binary files as input for
   the fortran pede program.
   This includes macros utilising cfortran.h to allow direct callability
   from fortran.

   initC() has to be called once in the beginning,
   followed by one or several calls to openC(..) to open one or several files.
   readC(..) is then called to read the records sequentially. It internally
   goes through all files as if it were only one, in contrast to the 
   fortran READ used in routine PEREAD of pede.F.

   if compiled with preprocessor macro USE_SHIFT_RFIO, uses libRFIO,
   i.e. includes shift.h instead of stdio.h

   written by Gero Flucke (gero.flucke@cern.ch) in 2006/7,
   update on July 14th, 2008
   last update on October 29th, 2008: return for file number in readC
*/

#ifdef USE_SHIFT_RFIO
#include <shift.h>
// or this??
// // needed?#define _LARGEFILE64_SOURCE
//#include <sys/types.h>
//#include "rfio_api.h"
#else
#include <stdio.h>
#endif
#include "cfortran.h"

/* ________ global variables used for file handling __________ */

#define MAXNUMFILES 490        /* should roughly match MFILES in mpinds.inc */
FILE *files[MAXNUMFILES];      /* pointers to opened binary files */
int fileReadOnce[MAXNUMFILES]; /* flag to printout record number once */
unsigned int numAllFiles;      /* number of opened files */
int fileIndex;                 /* index of current file */
int nRec;                      /* count records per file */


/*______________________________________________________________*/

void initC()
{
  /* initialises the 'global' variables used for file handling */
  { 
    int i = 0;
    for ( ; i < MAXNUMFILES; ++i) {
      files[i] = 0;
      fileReadOnce[i] = 0;
    }
  }
  numAllFiles = 0;
  fileIndex = -1;
  nRec = 0;
}
FCALLSCSUB0(initC,INITC,initc)

/*______________________________________________________________*/

/* void rewinC() */
/* { */
/*   /\* rewind all open files and start again with first file *\/ */

/*   unsigned int i = numAllFiles; */
/*   while (i--) rewind(files[i]); /\* postfix decrement! *\/ */
/*   fileIndex = 0; */
/* } */
/* FCALLSCSUB0(rewinC,REWINC,rewinc) */

/*______________________________________________________________*/

void resetC()
{
  /* start again with first file */
  if (fileIndex < 0) return; /* no file opened at all... */
  /* rewind(files[fileIndex]);  Does not work with rfio, so call: */
  fseek(files[fileIndex], 0L, SEEK_SET);
  clearerr(files[fileIndex]); /* These two should be the same as rewind... */
  if (!fileReadOnce[fileIndex]) {
    printf("readC: %d. file read the first time, read  %d records.\n",
              fileIndex+1, nRec);
    fileReadOnce[fileIndex] = 1;
  } 
  fileIndex = 0;
  nRec = 0; 
}
FCALLSCSUB0(resetC,RESETC,resetc)

/*______________________________________________________________*/

void openC(const char *fileName, int *errorFlag)
{
  /* No return value since to be called as subroutine from fortran, 
     errorFlag:
     * 0: if file opened and OK, 
     * 1: if too many files open,
     * 2: if file could not be opened 
     * 3: if file opened, but with error (can that happen?)
  */

  if (!errorFlag) return; /* 'printout' error? */

  if (numAllFiles >= MAXNUMFILES) {
    *errorFlag = 1;
  } else {
    files[numAllFiles] = fopen(fileName, "rb");
    if (!files[numAllFiles]) {
      *errorFlag = 2;
    } else if (ferror(files[numAllFiles])) {
      fclose(files[numAllFiles]);
      files[numAllFiles] = 0;
      *errorFlag = 3;
    } else {
      if (numAllFiles == 0) fileIndex = 0;
      ++numAllFiles; /* We have one more opened file! */
      *errorFlag = 0;
    }
  }
}
FCALLSCSUB2(openC,OPENC,openc,STRING,PINT)

/*______________________________________________________________*/

 void readC(float *bufferFloat, int *bufferInt, int *lengthBuffers,
	    int *nFileOut, int *errorFlag)
{
   /* No return value since to be called as subroutine from fortran,
      negative *errorFlag are errors, otherwise fine:
      * -1: pointer to a buffer or lengthBuffers/nFileOut are null
      * -2: problem reading record length
      * -4: given buffers too short for record
      * -8: problem with stream or EOF reading floats
      *-16: problem with stream or EOF reading ints
      *  0: reached end of all files (or read empty record?!)
      * >0: number of words (floats + integers) read and stored in buffers
      
      *nFileOut: returns number of the file the record is read from,
                 starting from 1 (not 0), but only if record read succesfully.
   */
   if (!errorFlag) return;
   *errorFlag = 0;
   if (fileIndex < 0) return; /* no file opened at all... */
   if (!bufferFloat || !bufferInt || !lengthBuffers || !nFileOut) {
     *errorFlag = -1;
     return;
   }

   /* read length of 'record' */
   int recordLength = 0; /* becomes number of words following in file */
   size_t nCheckR = fread(&recordLength, sizeof(recordLength), 1,
 			 files[fileIndex]);
   while (feof(files[fileIndex])) {
     /* rewind(files[fileIndex]);  Does not work with rfio, so call: */
     fseek(files[fileIndex], 0L, SEEK_SET);
     clearerr(files[fileIndex]); /* These two should be the same as rewind... */
     if (!fileReadOnce[fileIndex]) {
       printf("readC: %d. file read the first time, found %d records.\n",
	      fileIndex+1, nRec);
       fileReadOnce[fileIndex] = 1;
     }
     nRec = 0;
     if (fileIndex+1 >= numAllFiles) {
       *errorFlag = 0; /* Means EOF of last file. */
       fileIndex = (numAllFiles > 0 ? 0 : -1); /* Start first file, if any. */
       return;
     } else { /* Try next file! */
       ++fileIndex;
       nCheckR = fread(&recordLength, sizeof(recordLength),1,files[fileIndex]);
     }
   }

   if (1 != nCheckR || ferror(files[fileIndex])) {
     printf("readC: problem reading length of record %d, file %d\n",
 	   nRec+1, fileIndex);
     *errorFlag = -2;
     return;
   }

   if (recordLength/2 >= *lengthBuffers) {
     printf("readC: given buffers too short (%d, need > %d)\n", *lengthBuffers,
 	   recordLength/2);
     *errorFlag = -4;
     return;
   } else {
     *lengthBuffers = recordLength/2;
   }

   /* read floats (i.e. derivatives + value + sigma) */
   size_t nCheckF = fread(bufferFloat, sizeof(bufferFloat[0]), *lengthBuffers,
 			 files[fileIndex]);
   if (ferror(files[fileIndex]) || feof(files[fileIndex])
       || nCheckF != *lengthBuffers) {
     printf("readC: problem with stream or EOF reading floats\n");
     *errorFlag = -8;
     return;
   }

   /* read ints (i.e. parameter lables) */
   size_t nCheckI = fread(bufferInt, sizeof(bufferInt[0]), *lengthBuffers,
 			 files[fileIndex]);
   if (ferror(files[fileIndex]) || feof(files[fileIndex])
       || nCheckI != *lengthBuffers) {
     printf("readC: problem with stream or EOF reading ints\n");
     *errorFlag = -16;
     return;
   }

   ++nRec;
   *errorFlag = *lengthBuffers;
   *nFileOut = fileIndex + 1; /* As output, starting from 1. */
 }
FCALLSCSUB5(readC,READC,readc,PFLOAT,PINT,PINT,PINT,PINT)
