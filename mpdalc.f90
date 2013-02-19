!> \file
!! Dynamic memory management.
!!
!! \author C. Kleinwort, DESY, 2012 (Claus.Kleinwort@desy.de)

!> (De)Allocate vectors and arrays.
MODULE mpdalc
    USE mpdef
    IMPLICIT NONE
    SAVE
    ! variables
    INTEGER(kind=large) :: numwordsalloc = 0 !< current dynamic memory allocation (words)
    INTEGER(kind=large) :: maxwordsalloc = 0 !< peak dynamic memory allocation (words)
    INTEGER :: nummpalloc = 0     !< number of dynamic allocations
    INTEGER :: nummpdealloc = 0   !< number of dynamic deallocations
    INTEGER :: printflagalloc = 0 !< print flag for dynamic allocations

    !> allocate array
    INTERFACE mpalloc
        MODULE PROCEDURE mpallocdvec, mpallocfvec, mpallocivec, &
            mpallocfarr, mpallociarr, mpalloclarr, mpalloclist, mpalloccvec
    END INTERFACE mpalloc
    !> deallocate array
    INTERFACE mpdealloc
        MODULE PROCEDURE mpdeallocdvec, mpdeallocfvec, mpdeallocivec, &
            mpdeallocfarr, mpdeallociarr, mpdealloclarr, mpdealloclist, mpdealloccvec
    END INTERFACE mpdealloc

CONTAINS
    ! allocate dynamic vector or array
    !> allocate (1D) double precision array
    SUBROUTINE mpallocdvec(array,length,text)
        DOUBLE PRECISION, DIMENSION(:), INTENT(IN OUT), ALLOCATABLE :: array
        INTEGER(kind=large), INTENT(IN) :: length
        CHARACTER (LEN=*), INTENT(IN) :: text

        INTEGER :: ifail
        ALLOCATE (array(length),stat=ifail)
        CALL mpalloccheck(ifail,2*length,text)
    END SUBROUTINE mpallocdvec

    !> allocate (1D) single precision array
    SUBROUTINE mpallocfvec(array,length,text)
        REAL, DIMENSION(:), INTENT(IN OUT), ALLOCATABLE :: array
        INTEGER(kind=large), INTENT(IN) :: length
        CHARACTER (LEN=*), INTENT(IN) :: text

        INTEGER :: ifail
        ALLOCATE (array(length),stat=ifail)
        CALL mpalloccheck(ifail,length,text)
    END SUBROUTINE mpallocfvec

    !> allocate (1D) integer array
    SUBROUTINE mpallocivec(array,length,text)
        INTEGER, DIMENSION(:), INTENT(IN OUT), ALLOCATABLE :: array
        INTEGER(kind=large), INTENT(IN) :: length
        CHARACTER (LEN=*), INTENT(IN) :: text

        INTEGER :: ifail
        ALLOCATE (array(length),stat=ifail)
        CALL mpalloccheck(ifail,length,text)
    END SUBROUTINE mpallocivec

    !> allocate (2D) single precision array
    SUBROUTINE mpallocfarr(array,rows,cols,text)
        REAL, DIMENSION(:,:), INTENT(IN OUT), ALLOCATABLE :: array
        INTEGER(kind=large), INTENT(IN) :: rows
        INTEGER(kind=large), INTENT(IN) :: cols
        CHARACTER (LEN=*), INTENT(IN)  :: text

        INTEGER :: ifail
        ALLOCATE (array(rows,cols),stat=ifail)
        CALL mpalloccheck(ifail,rows*cols,text)
    END SUBROUTINE mpallocfarr

    !> allocate (2D) integer array
    SUBROUTINE mpallociarr(array,rows,cols,text)
        INTEGER, DIMENSION(:,:), INTENT(IN OUT), ALLOCATABLE :: array
        INTEGER(kind=large), INTENT(IN) :: rows
        INTEGER(kind=large), INTENT(IN) :: cols
        CHARACTER (LEN=*), INTENT(IN)  :: text

        INTEGER :: ifail
        ALLOCATE (array(rows,cols),stat=ifail)
        CALL mpalloccheck(ifail,rows*cols,text)
    END SUBROUTINE mpallociarr

    !> allocate (2D) large integer array
    SUBROUTINE mpalloclarr(array,rows,cols,text)
        INTEGER(kind=large), DIMENSION(:,:), INTENT(IN OUT), ALLOCATABLE :: array
        INTEGER(kind=large), INTENT(IN) :: rows
        INTEGER(kind=large), INTENT(IN) :: cols
        CHARACTER (LEN=*), INTENT(IN)  :: text

        INTEGER :: ifail
        ALLOCATE (array(rows,cols),stat=ifail)
        CALL mpalloccheck(ifail,rows*cols*large/4,text)
    END SUBROUTINE mpalloclarr

    !> allocate (1D) list item array
    SUBROUTINE mpalloclist(array,length,text)
        TYPE(listItem), DIMENSION(:), INTENT(IN OUT), ALLOCATABLE :: array
        INTEGER(kind=large), INTENT(IN) :: length
        CHARACTER (LEN=*), INTENT(IN) :: text

        INTEGER :: ifail
        ALLOCATE (array(length),stat=ifail)
        CALL mpalloccheck(ifail,length*2,text)
    END SUBROUTINE mpalloclist

    !> allocate (1D) character array
    SUBROUTINE mpalloccvec(array,length,text)
        CHARACTER, DIMENSION(:), INTENT(IN OUT), ALLOCATABLE :: array
        INTEGER(kind=large), INTENT(IN) :: length
        CHARACTER (LEN=*), INTENT(IN) :: text

        INTEGER :: ifail
        ALLOCATE (array(length),stat=ifail)
        CALL mpalloccheck(ifail,(length+3)/4,text)
    END SUBROUTINE mpalloccvec

    !> check allocation
    SUBROUTINE mpalloccheck(ifail,numwords,text)
        INTEGER, INTENT(IN) :: ifail
        INTEGER(kind=large), INTENT(IN) :: numwords
        CHARACTER (LEN=*), INTENT(IN)  :: text
        IF (ifail == 0) THEN
            nummpalloc=nummpalloc+1
            numwordsalloc = numwordsalloc + numwords
            maxwordsalloc = MAX(maxwordsalloc, numwordsalloc)
            IF (printflagalloc /= 0) THEN
                print *, ' MPALLOC allocated ', numwords, ' words for : ', text
                print *, ' words used ', numwordsalloc, maxwordsalloc
            ENDIF
        ELSE
            print *, ' MPALLOC failed to allocate ', numwords, ' words for : ', text
            print *, ' MPALLOC words used ', numwordsalloc, maxwordsalloc
            print *, ' MPALLOC stat = ', ifail
            STOP
        ENDIF
    END SUBROUTINE mpalloccheck
    ! deallocate dynamic vector or array
    !> deallocate (1D) double precision array
    SUBROUTINE mpdeallocdvec(array)
        DOUBLE PRECISION, DIMENSION(:), INTENT(IN OUT), ALLOCATABLE :: array

        INTEGER :: ifail
        INTEGER(kind=large) :: isize
        isize = 2*size(array,kind=large)
        DEALLOCATE (array,stat=ifail)
        CALL mpdealloccheck(ifail,isize)
    END SUBROUTINE mpdeallocdvec

    !> deallocate (1D) single precision array
    SUBROUTINE mpdeallocfvec(array)
        REAL, DIMENSION(:), INTENT(IN OUT), ALLOCATABLE :: array

        INTEGER :: ifail
        INTEGER(kind=large) :: isize
        isize = size(array,kind=large)
        DEALLOCATE (array,stat=ifail)
        CALL mpdealloccheck(ifail,isize)
    END SUBROUTINE mpdeallocfvec

    !> deallocate (1D) integer array
    SUBROUTINE mpdeallocivec(array)
        INTEGER, DIMENSION(:), INTENT(IN OUT), ALLOCATABLE :: array

        INTEGER :: ifail
        INTEGER(kind=large) :: isize
        isize = size(array,kind=large)
        DEALLOCATE (array,stat=ifail)
        CALL mpdealloccheck(ifail,isize)
    END SUBROUTINE mpdeallocivec

    !> allocate (2D) single precision array
    SUBROUTINE mpdeallocfarr(array)
        REAL, DIMENSION(:,:), INTENT(IN OUT), ALLOCATABLE :: array

        INTEGER :: ifail
        INTEGER(kind=large) :: isize
        isize = size(array,kind=large)
        DEALLOCATE (array,stat=ifail)
        CALL mpdealloccheck(ifail,isize)
    END SUBROUTINE mpdeallocfarr

    !> allocate (2D) integer array
    SUBROUTINE mpdeallociarr(array)
        INTEGER, DIMENSION(:,:), INTENT(IN OUT), ALLOCATABLE :: array

        INTEGER :: ifail
        INTEGER(kind=large) :: isize
        isize = size(array,kind=large)
        DEALLOCATE (array,stat=ifail)
        CALL mpdealloccheck(ifail,isize)
    END SUBROUTINE mpdeallociarr

     !> deallocate (2D) large integer array
    SUBROUTINE mpdealloclarr(array)
        INTEGER(kind=large), DIMENSION(:,:), INTENT(IN OUT), ALLOCATABLE :: array

        INTEGER :: ifail
        INTEGER(kind=large) :: isize
        isize = size(array,kind=large)*large/4
        DEALLOCATE (array,stat=ifail)
        CALL mpdealloccheck(ifail,isize)
    END SUBROUTINE mpdealloclarr

    !> deallocate (1D) list item array
    SUBROUTINE mpdealloclist(array)
        TYPE(listItem), DIMENSION(:), INTENT(IN OUT), ALLOCATABLE :: array

        INTEGER :: ifail
        INTEGER(kind=large) :: isize
        isize = 2*size(array,kind=large)
        DEALLOCATE (array,stat=ifail)
        CALL mpdealloccheck(ifail,isize)
    END SUBROUTINE mpdealloclist

    !> deallocate (1D) character array
    SUBROUTINE mpdealloccvec(array)
        CHARACTER, DIMENSION(:), INTENT(IN OUT), ALLOCATABLE :: array

        INTEGER :: ifail
        INTEGER(kind=large) :: isize
        isize = (size(array,kind=large)+3)/4
        DEALLOCATE (array,stat=ifail)
        CALL mpdealloccheck(ifail,isize)
    END SUBROUTINE mpdealloccvec

    !> check deallocation
    SUBROUTINE mpdealloccheck(ifail,numwords)
        INTEGER, INTENT(IN) :: ifail
        INTEGER(kind=large), INTENT(IN) :: numwords
        IF (ifail == 0) THEN
            numwordsalloc = numwordsalloc - numwords
            nummpdealloc=nummpdealloc+1
            IF (printflagalloc /= 0) THEN
                print *, ' MPDEALLOC deallocated ', numwords, ' words '
                print *, ' words used ', numwordsalloc, maxwordsalloc
            ENDIF
        ELSE
            print *, ' MPDEALLOC failed to deallocate ', numwords, ' words'
            print *, ' MPDEALLOC words used ', numwordsalloc, maxwordsalloc
            print *, ' MPDEALLOC stat = ', ifail
            STOP
        ENDIF
    END SUBROUTINE mpdealloccheck

END MODULE mpdalc
