!> \file
!! Definitions.

!> Definition of constants.
MODULE mpdef
    IMPLICIT NONE
    SAVE
    ! constants
    INTEGER, PARAMETER :: large=LARGE_SIZE !< size of LARGE integers (*4 or *8) from preprocessor
    INTEGER, PARAMETER :: lscale=large/4   !< number of (32bit) words
    INTEGER, PARAMETER :: maxi4=2147483647 !< max. INTEGER*4
    !> list items from steering file
    TYPE listItem
        INTEGER :: label
        REAL :: value
    END TYPE listItem
END MODULE mpdef
