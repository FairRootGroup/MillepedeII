!> \file
!! Bit field counters.
!!
!! Count pairs of global parameters for sparse storage of global matrix,
!! apply pair entries cut and build (compressed) sparsity structure (row offsets, column lists).
!!
!! In sparse storage mode for each row the list of column indices with non zero elements
!! (and those elements) are stored. With compression this list is represented by the
!! first column and their number for continous regions (encoded in single 32bit words).
!! Rare elements may be stored in single precision.
!!

!> Bit field data.
MODULE mpbits
    USE mpdef
    IMPLICIT NONE

    INTEGER(KIND=large) :: ndimb !< dimension for bit (field) array
    INTEGER :: n      !< matrix size
    INTEGER :: ibfw   !< bit field width
    INTEGER :: mxcnt  !< max value for bit field counters
    INTEGER :: nencdm !< max value for column counter
    INTEGER :: nencdb !< number of bits for encoding column counter
    INTEGER :: nthrd  !< number of threads
    INTEGER, DIMENSION(:), ALLOCATABLE :: bitFieldCounters !< fit field counters for global parameters pairs

END MODULE mpbits

!> Fill bit fields (counters).
!!
!! \param [in]    im     first index
!! \param [in]    jm     second index
!! \param [in]    inc    increment (usually 1)
!!
SUBROUTINE inbits(im,jm,inc)        ! include element (I,J)
    USE mpbits

    INTEGER, INTENT(IN) :: im
    INTEGER, INTENT(IN) :: jm
    INTEGER, INTENT(IN) :: inc

    INTEGER(KIND=large) :: l
    INTEGER(KIND=large) :: ll
    INTEGER :: i
    INTEGER :: j
    INTEGER :: noffj
    INTEGER :: m
    INTEGER :: mm
    INTEGER :: icount
    INTEGER :: ib
    INTEGER :: jcount
    INTEGER*8 :: noffi
    LOGICAL :: btest

    IF(im == jm) RETURN  ! diagonal
    j=MIN(im,jm)
    i=MAX(im,jm)
    IF(j <= 0) RETURN    ! out low
    IF(i > n) RETURN    ! out high
    noffi=int8(i-1)*int8(i-2)*int8(ibfw)/2 ! for J=1
    noffj=(j-1)*ibfw
    l=noffi/32+i+noffj/32 ! row offset + column offset
    !     add I instead of 1 to keep bit maps of different rows in different words (openMP !)
    m=MOD(noffj,32)
    IF (ibfw <= 1) THEN
        bitFieldCounters(l)=ibset(bitFieldCounters(l),m)
    ELSE
        !        get counter from bit field
        ll=l
        mm=m
        icount=0
        DO ib=0,ibfw-1
            IF (btest(bitFieldCounters(ll),mm)) icount=ibset(icount,ib)
            mm=mm+1
            IF (mm >= 32) THEN
                ll=ll+1
                mm=mm-32
            END IF
        END DO
        !        increment
        jcount=icount
        icount=MIN(icount+inc,mxcnt)
        !        store counter into bit field
        IF (icount /= jcount) THEN
            ll=l
            mm=m
            DO ib=0,ibfw-1
                IF (btest(icount,ib)) THEN
                    bitFieldCounters(ll)=ibset(bitFieldCounters(ll),mm)
                ELSE
                    bitFieldCounters(ll)=ibclr(bitFieldCounters(ll),mm)
                END IF
                mm=mm+1
                IF (mm >= 32) THEN
                    ll=ll+1
                    mm=mm-32
                END IF
            END DO
        END IF
    END IF
    RETURN

END SUBROUTINE inbits

!> Calculate bit (field) array size, encoding.
!!
!! \param [in]    in        matrix size
!! \param [out]   idimb     dimension for bit (field) array
!! \param [out]   iencdb    number of bits for encoding column counter
!! \param [in]    jbfw      bit field width
!!
SUBROUTINE clbits(in,idimb,iencdb,jbfw)
    USE mpbits
    USE mpdalc

    INTEGER, INTENT(IN) :: in
    INTEGER(KIND=large), INTENT(OUT) :: idimb
    INTEGER, INTENT(OUT) :: iencdb
    INTEGER, INTENT(IN) :: jbfw

    INTEGER*8 :: noffd
    INTEGER :: i
    INTEGER :: mb
    INTEGER :: nbcol
    !$    INTEGER :: OMP_GET_MAX_THREADS

    n=in
    ibfw=jbfw
    mxcnt=2**ibfw-1
    noffd=int8(n)*int8(n-1)*int8(ibfw)/2
    ndimb=noffd/32+n
    idimb=ndimb
    mb=INT(4.0E-6*FLOAT(ndimb))
    WRITE(*,*) ' '
    WRITE(*,*) 'CLBITS: symmetric matrix of dimension',n
    WRITE(*,*) 'CLBITS: off-diagonal elements',noffd
    IF (mb > 0) THEN
        WRITE(*,*) 'CLBITS: dimension of bit-array',ndimb , '(',mb,'MB)'
    ELSE
        WRITE(*,*) 'CLBITS: dimension of bit-array',ndimb , '(< 1 MB)'
    END IF
    CALL mpalloc(bitFieldCounters,ndimb,'INBITS: bit storage')
    bitFieldCounters=0
    !     encoding for compression
    nbcol=16    ! 16 bits for column number, 16 bits for column counter
    DO i=16,30
        IF (btest(n,i)) nbcol=i+1 ! more bits for column number
    END DO
    nencdb=32-nbcol
    iencdb=nencdb
    nencdm=ishft(1,nencdb)-1
    nthrd=1
    !$ NTHRD=OMP_GET_MAX_THREADS()
    RETURN
END SUBROUTINE clbits

!> Analyze bit fields.
!!
!! \param [out]     ndims   (1): (reduced) size of bit array; (2): size of column lists;
!!                          (3/4): number of (double/single precision) off diagonal elements;
!! \param[out]     ncmprs  compression info (per row)
!! \param[out]     nsparr  row offsets
!! \param[in]      mnpair  min. entries for pair
!! \param[in]      ihst    >0: histogram number
!! \param[in]      jcmprs  <>0: compress row information (column indices)
!!
SUBROUTINE ndbits(ndims,ncmprs,nsparr,mnpair,ihst,jcmprs)
    USE mpbits

    INTEGER(KIND=large), DIMENSION(4), INTENT(OUT) :: ndims
    INTEGER, DIMENSION(:), INTENT(OUT) :: ncmprs
    INTEGER(KIND=large), DIMENSION(:,:), INTENT(OUT) :: nsparr
    INTEGER, INTENT(IN) :: mnpair
    INTEGER, INTENT(IN) :: ihst
    INTEGER, INTENT(IN) :: jcmprs

    INTEGER :: nwcp(0:1)
    INTEGER :: irgn(2)
    INTEGER :: inr(2)
    INTEGER :: ichunk
    INTEGER :: i
    INTEGER :: j
    INTEGER :: jb
    INTEGER :: m
    INTEGER :: last
    INTEGER :: lrgn
    INTEGER :: next
    INTEGER :: icp
    INTEGER :: kbfw
    INTEGER :: mm
    INTEGER :: jp
    INTEGER :: icmprs
    INTEGER :: nj
    INTEGER :: ib
    INTEGER :: ir
    INTEGER :: icount
    INTEGER :: iproc
    INTEGER :: jbfw
    INTEGER :: k
    INTEGER :: mb
    INTEGER :: n1
    INTEGER(KIND=large) :: ll
    INTEGER(KIND=large) :: lb
    INTEGER(KIND=large) :: nin
    INTEGER(KIND=large) :: ntot
    INTEGER*8 :: noffi
    REAL :: cpr
    REAL :: fracu
    REAL :: fracz
    LOGICAL :: btest
    !$    INTEGER :: OMP_GET_THREAD_NUM
    ndims(1)=ndimb
    ndims(2)=0
    ndims(3)=0
    ndims(4)=0
    ntot=0
    icmprs=jcmprs
    !     reduce bit field counters to bits
    kbfw=1
    ll=0
    lb=0
    ichunk=MIN((n+nthrd-1)/nthrd/32+1,256)
    IF (ibfw > 1.OR.icmprs /= 0) THEN
        IF (ibfw > 1.AND.icmprs > 0) kbfw=2 ! need to tag single precision entries
        jbfw=kbfw ! temp. bit field width
        IF (nthrd > 1) jbfw=ibfw ! don't mix rows in bitFieldCounters
        ! parallelize row loop
        ! private copy of NDIMS,NTOT for each thread, combined at end, init with 0.
        !$OMP  PARALLEL DO &
        !$OMP  PRIVATE(I,LL,MM,LB,MB,INR,IRGN,LAST,LRGN, &
        !$OMP          J,ICOUNT,NEXT,IB,ICP,NWCP,JB,JP,IR) &
        !$OMP  REDUCTION(+:NDIMS,NTOT) &
        !$OMP  SCHEDULE(DYNAMIC,ICHUNK)
        DO i=1,n
            noffi=int8(i-1)*int8(i-2)*int8(ibfw)/2
            ll=noffi/32+i
            mm=0
            noffi=int8(i-1)*int8(i-2)*int8(jbfw)/2
            lb=noffi/32+i
            mb=0
            inr(1)=0
            inr(2)=0
            irgn(1)=0
            irgn(2)=0
            last=0
            lrgn=0
            iproc=0
            !$ IPROC=OMP_GET_THREAD_NUM()         ! thread number

            DO j=1,i-1

                icount=0
                next=0
                DO ib=0,ibfw-1
                    IF (btest(bitFieldCounters(ll),mm)) icount=ibset(icount,ib)
                    mm=mm+1
                    IF (mm >= 32) THEN
                        ll=ll+1
                        mm=mm-32
                    END IF
                END DO
                DO jb=0,kbfw-1
                    bitFieldCounters(lb)=ibclr(bitFieldCounters(lb),mb+jb)
                END DO

                IF (icount > 0) THEN
                    ntot=ntot+1
                    IF (iproc == 0.AND.ihst > 0) CALL hmpent(ihst,FLOAT(icount))
                END IF
                !              keep pair ?
                IF (icount >= mnpair) THEN
                    next=1 ! double
                    IF (icount <= icmprs.AND.icmprs > 0) next=2 ! single
                    inr(next)=inr(next)+1
                    bitFieldCounters(lb)=ibset(bitFieldCounters(lb),mb+next-1)
                    IF (next /= last.OR.lrgn >= nencdm) THEN
                        irgn(next)=irgn(next)+1
                        lrgn=0
                    END IF
                    lrgn=lrgn+1
                END IF
                mb=mb+kbfw
                IF (mb >= 32) THEN
                    lb=lb+1
                    mb=mb-32
                END IF
                last=next
            END DO

            DO jp=1,kbfw
                icp=0
                nwcp(0)=inr(jp)                       ! list of column indices (default)
                IF (inr(jp) > 0) THEN
                    nwcp(1)=irgn(jp)+(irgn(jp)+7)/8    ! list of regions of consecutive columns
                    !                 compress row ?
                    IF (nwcp(1) < nwcp(0).AND.icmprs /= 0) THEN
                        icp=1
                        ncmprs(i+n*(jp-1))=irgn(jp)
                    END IF
                    ndims(2)   =ndims(2)   +nwcp(icp)
                    ndims(jp+2)=ndims(jp+2)+nwcp(0)
                END IF
                ir=i+(n+1)*(jp-1)
                nsparr(1,ir+1)=nwcp(icp)
                nsparr(2,ir+1)=nwcp(0)
            END DO

        END DO

        !$OMP END PARALLEL DO
        !        sum up, fill row offsets
        lb=1
        n1=0
        ll=n+1
        DO jp=1,kbfw
            DO i=1,n
                n1=n1+1
                nsparr(1,n1)=lb
                nsparr(2,n1)=ll
                lb=lb+nsparr(1,n1+1)
                ll=ll+nsparr(2,n1+1)
            END DO
            n1=n1+1
            nsparr(1,n1)=lb
            nsparr(2,n1)=ll
            ll=1
        END DO

        IF (jbfw /= kbfw) THEN ! move bit fields
            DO i=1,n
                noffi=int8(i-1)*int8(i-2)*int8(jbfw)/2
                ll=noffi/32+i
                noffi=int8(i-1)*int8(i-2)*int8(kbfw)/2
                lb=noffi/32+i
                nj=((i-1)*kbfw)/32
                DO k=0,nj
                    bitFieldCounters(lb+k)=bitFieldCounters(ll+k)
                END DO
            END DO
        END IF

        ibfw=kbfw
        noffi=int8(n)*int8(n-1)*int8(ibfw)/2
        ndimb=noffi/32+n
        ndims(1)=ndimb

    ELSE

        nin=0
        nsparr(1,1)=1
        nsparr(2,1)=n+1
        n1=1
        DO i=1,n
            noffi=int8(i-1)*int8(i-2)/2
            ll=noffi/32+i
            nj=((i-1)*kbfw)/32
            DO k=0,nj
                DO m=0,31
                    IF(btest(bitFieldCounters(ll+k),m)) nin=nin+1
                END DO
            END DO
            n1=n1+1
            nsparr(1,n1)=nsparr(1,1)+nin
            nsparr(2,n1)=nsparr(2,1)+nin
        END DO
        ndims(2)=nin
        ndims(3)=nin
        ntot=nin

    END IF

    nin=ndims(3)+ndims(4)
    fracz=200.0*FLOAT(ntot)/FLOAT(n)/FLOAT(n-1)
    fracu=200.0*FLOAT(nin)/FLOAT(n)/FLOAT(n-1)
    WRITE(*,*) ' '
    WRITE(*,*) 'NDBITS: number of diagonal elements',n
    WRITE(*,*) 'NDBITS: number of used off-diagonal elements',nin
    WRITE(*,1000) 'fraction of non-zero off-diagonal elements', fracz
    WRITE(*,1000) 'fraction of used off-diagonal elements', fracu
    IF (icmprs /= 0) THEN
        cpr=100.0*FLOAT(ndims(2)+2*ndims(3)+ndims(4))/FLOAT(3*nin)
        WRITE(*,1000) 'compression ratio for off-diagonal elements', cpr
    END IF
1000 FORMAT(' NDBITS: ',a,f6.2,' %')
    RETURN
END SUBROUTINE ndbits

!> Check sparsity of matrix.
!!
!! \param [out]    ndims   (1): (reduced) size of bit array; (2): size of column lists;
!!                         (3/4): number of (double/single precision) off diagonal elements;
!! \param[in]      mnpair  min. entries for pair
!! \param[in]      jcmprs  <>0: compress row information (column indices)
!!
SUBROUTINE ckbits(ndims,mnpair,jcmprs)
    USE mpbits

    INTEGER(KIND=large), DIMENSION(4), INTENT(OUT) :: ndims
    INTEGER, INTENT(IN) :: mnpair
    INTEGER, INTENT(IN) :: jcmprs

    INTEGER :: nwcp(0:1)
    INTEGER :: irgn(2)
    INTEGER :: inr(2)
    INTEGER(KIND=large) :: ll
    INTEGER*8 :: noffi
    INTEGER :: i
    INTEGER :: j
    INTEGER :: last
    INTEGER :: lrgn
    INTEGER :: next
    INTEGER :: icp
    INTEGER :: ib
    INTEGER :: icount
    INTEGER :: icmprs
    INTEGER :: kbfw
    INTEGER :: jp
    INTEGER :: mm
    LOGICAL :: btest

    DO i=1,4
        ndims(i)=0
    END DO
    icmprs=jcmprs
    kbfw=1
    IF (ibfw > 1.AND.icmprs > 0) kbfw=2
    ll=0

    DO i=1,n
        noffi=int8(i-1)*int8(i-2)*int8(ibfw)/2
        ll=noffi/32+i
        mm=0
        inr(1)=0
        inr(2)=0
        irgn(1)=0
        irgn(2)=0
        last=0
        lrgn=0
        DO j=1,i-1
            icount=0
            next=0
            DO ib=0,ibfw-1
                IF (btest(bitFieldCounters(ll),mm)) icount=ibset(icount,ib)
                mm=mm+1
                IF (mm >= 32) THEN
                    ll=ll+1
                    mm=mm-32
                END IF
            END DO

            IF (icount > 0) ndims(1)=ndims(1)+1
            !           keep pair ?
            IF (icount >= mnpair) THEN
                next=1 ! double
                IF (icount <= icmprs.AND.icmprs > 0) next=2 ! single
                inr(next)=inr(next)+1
                IF (next /= last.OR.lrgn >= nencdm) THEN
                    irgn(next)=irgn(next)+1
                    lrgn=0
                END IF
                lrgn=lrgn+1
            END IF
            last=next
        END DO

        IF (icmprs /= 0) THEN
            DO jp=1,kbfw
                IF (inr(jp) > 0) THEN
                    icp=0
                    nwcp(0)=inr(jp)                    ! list of column indices (default)
                    nwcp(1)=irgn(jp)+(irgn(jp)+7)/8    ! list of regions of consecutive columns
                    !                 compress row ?
                    IF (nwcp(1) < nwcp(0)) icp=1
                    ndims(2)   =ndims(2)   +nwcp(icp)
                    ndims(jp+2)=ndims(jp+2)+nwcp(0)
                END IF
            END DO
        ELSE
            ndims(2)=ndims(2)+inr(1)
            ndims(3)=ndims(3)+inr(1)
        END IF

    END DO

    RETURN
END SUBROUTINE ckbits

!> Create sparsity information.
!!
!! \param[in ]    nsparr  row offsets
!! \param[out]    nsparc  column indices
!! \param[in]     ncmprs  compression info (per row)
!!
SUBROUTINE spbits(nsparr,nsparc,ncmprs)               ! collect elements
    USE mpbits
    USE mpdalc

    INTEGER(KIND=large), DIMENSION(:,:), INTENT(IN) :: nsparr
    INTEGER, DIMENSION(:), INTENT(OUT) :: nsparc
    INTEGER, DIMENSION(:), INTENT(IN) :: ncmprs

    INTEGER(KIND=large) :: kl
    INTEGER(KIND=large) :: l
    INTEGER(KIND=large) :: ll
    INTEGER(KIND=large) :: l1
    INTEGER(KIND=large) :: k8
    INTEGER(KIND=large) :: n1
    INTEGER*8 :: noffi
    INTEGER :: i
    INTEGER :: j
    INTEGER :: j1
    INTEGER :: jb
    INTEGER :: jn
    INTEGER :: m
    INTEGER :: ichunk
    INTEGER :: next
    INTEGER :: last
    INTEGER :: lrgn
    INTEGER :: nrgn
    INTEGER :: nrgn8
    LOGICAL :: btest

    ichunk=MIN((n+nthrd-1)/nthrd/32+1,256)

    DO jb=0,ibfw-1
        ! parallelize row loop
        !$OMP  PARALLEL DO &
        !$OMP  PRIVATE(I,N1,NOFFI,L,M,KL,L1,NRGN,NRGN8,K8, &
        !$OMP          LAST,LRGN,LL,J1,JN,J,NEXT) &
        !$OMP  SCHEDULE(DYNAMIC,ICHUNK)
        DO i=1,n
            n1=i+jb*(n+1)
            noffi=int8(i-1)*int8(i-2)*int8(ibfw)/2
            l=noffi/32+i
            m=jb
            kl=nsparr(1,n1)-1  ! pointer to row in NSPARC
            l1=nsparr(2,n1)    ! pointer to row in sparse matrix
            nrgn=ncmprs(i+n*jb)! compression  (number of consecutive regions)
            nrgn8=(nrgn+7)/8   ! number of groups (1 offset per group)
            k8=kl
            kl=kl+nrgn8        ! reserve space of offsets
            last=0
            lrgn=0
            ll=l1-1
            j1=0
            jn=0

            DO j=1,i-1         !  loop for off-diagonal elements
                next=0
                IF(bitFieldCounters(l) /= 0) THEN
                    IF(btest(bitFieldCounters(l),m)) THEN
                        ll=ll+1
                        IF (nrgn <= 0) THEN
                            kl=kl+1
                            nsparc(kl)=j ! column index
                        ELSE
                            next=1
                            IF (last == 0.OR.jn >= nencdm) THEN
                                IF (MOD(lrgn,8) == 0) THEN
                                    k8=k8+1
                                    nsparc(k8)=INT(ll-l1)
                                END IF
                                lrgn=lrgn+1
                                kl=kl+1
                                j1=ishft(j,nencdb)
                                jn=0
                            END IF
                            jn=jn+1
                            nsparc(kl)=ior(j1,jn)
                        END IF
                    END IF
                END IF
                last=next
                m=m+ibfw
                IF (m >= 32) THEN
                    m=m-32
                    l=l+1
                END IF

            END DO
        END DO
    !$OMP END PARALLEL DO

    END DO

    n1=(n+1)*ibfw
    WRITE(*,*) ' '
    WRITE(*,*) 'SPBITS: sparse structure constructed ',nsparr(1,n1), ' words'
    WRITE(*,*) 'SPBITS: dimension parameter of matrix',nsparr(2,1)-1
    IF (ibfw <= 1) THEN
        WRITE(*,*) 'SPBITS: index of last used location',nsparr(2,n1)-1
    ELSE
        WRITE(*,*) 'SPBITS: index of last used double',nsparr(2,n1/2)-1
        WRITE(*,*) 'SPBITS: index of last used single',nsparr(2,n1)-1
    END IF
    CALL mpdealloc(bitFieldCounters)
    RETURN
END SUBROUTINE spbits


