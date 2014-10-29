!> \file
!! Definitions.

!> Definition of constants.
MODULE mpdef
    IMPLICIT NONE
    SAVE
    ! precision constants
    INTRINSIC :: selected_real_kind
    INTRINSIC :: selected_int_kind
    INTEGER, PARAMETER :: mpi4  = selected_int_kind(9)         !>  4 byte integer
    INTEGER, PARAMETER :: mpi8  = selected_int_kind(18)        !>  8 byte integer
    INTEGER, PARAMETER :: mpr4  = selected_real_kind(6, 37)    !>  4 byte float
    INTEGER, PARAMETER :: mpr8  = selected_real_kind(15, 307)  !>  8 byte float
    INTEGER, PARAMETER :: mpr16 = selected_real_kind(33, 4931) !> 16 byte float, gcc needs libquadmath    INTEGER, PARAMETER :: mpi = selected_int_kind(9)         !>  4 byte integer
    INTEGER, PARAMETER :: mpi  = mpi4                          !>  integer
    INTEGER, PARAMETER :: mpl  = mpi8                          !>  long integer
    INTEGER, PARAMETER :: mps  = mpr4                          !>  single precision
    INTEGER, PARAMETER :: mpd  = mpr8                          !>  double precision
    !> list items from steering file
    TYPE listItem
        INTEGER(mpi) :: label
        REAL(mpd) :: value
    END TYPE listItem
END MODULE mpdef
