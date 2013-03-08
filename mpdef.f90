!> \file
!! Definitions.

!> Definition of constants.
MODULE mpdef
    IMPLICIT NONE
    SAVE
    ! constants
    INTEGER, PARAMETER :: large=8          !< size of LARGE integers (*8)
    INTEGER, PARAMETER :: lscale=large/4   !< number of (32bit) words
    INTEGER, PARAMETER :: maxi4=2147483647 !< max. INTEGER*4
    !> list items from steering file
    TYPE listItem
        INTEGER :: label
        REAL :: value
    END TYPE listItem
END MODULE mpdef
