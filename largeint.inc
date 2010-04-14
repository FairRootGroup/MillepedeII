* Large Integers are INTEGER(KIND=LARGE)
* if not understood by compiler replace (KIND=LARGE) by *4 or *8 with preprocessor
      INTEGER LARGE,LSCALE,MAXI4
      PARAMETER       (LARGE=LARGE_SIZE) ! size of LARGE integers (*4 or *8) from preprocessor
      PARAMETER       (LSCALE=LARGE/4)   ! number of (32bit) words
      PARAMETER       (MAXI4=2147483647) ! max. INTEGER*4
*
