!> \file
!! Data structures.

!> Parameters, variables, dynamic arrays.
!!
!! For parameters which can be set from command line or
!! steering files more details are available in: \ref option_page.

MODULE mpmod
    USE mpdef
    IMPLICIT NONE
    SAVE
    ! steering parameters
    INTEGER :: ictest=0  !< test mode '-t'
    INTEGER :: metsol=0  !< solution method (1: inversion, 2: diagonalization, 3: \ref minres "MINRES")
    INTEGER :: matsto=2  !< (global) matrix storage mode (1: full, 2: sparse)
    INTEGER :: mprint=1  !< print flag (0: minimal, 1: normal, >1: more)
    INTEGER :: mdebug=0  !< debug flag (number of records to print)
    INTEGER :: mdebg2=10 !< number of measurements for record debug printout
    INTEGER :: mreqen=10 !< required number of entries (for variable global parameter)
    INTEGER :: mitera=1  !< number of iterations
    INTEGER :: nloopn=0  !< number of data reading, fitting loops
    INTEGER :: mbandw=0  !< band width of preconditioner matrix
    INTEGER :: lunkno=0  !< flag for unkown keywords
    INTEGER :: lhuber=0  !< Huber down-weighting flag
    REAL    :: chicut=0.0  !< cut in terms of 3-sigma cut, first iteration
    REAL    :: chirem=0.0  !< cut in terms of 3-sigma cut, other iterations, approaching 1.
    REAL    :: chhuge=50.0 !< cut in terms of 3-sigma for unreasonable data, all iterations
    INTEGER :: nrecpr=0  !< record number with printout
    INTEGER :: nrecp2=0  !< record number with printout
    INTEGER :: nrec1 =0  !< record number with largest residual
    INTEGER :: nrec2 =0  !< record number with largest chi^2/Ndf
    REAL    :: value1=0.0!< largest residual
    REAL    :: value2=0.0!< largest chi^2/Ndf
    REAL    :: dwcut=0.0 !< down-weight fraction cut
    INTEGER :: isubit=0  !< subito flag '-s'
    REAL    :: wolfc1=0.0!< C_1 of strong Wolfe condition
    REAL    :: wolfc2=0.0!< C_2 of strong Wolfe condition
    DOUBLE PRECISION :: mrestl=1.0D-06 !< tolerance criterion for MINRES
    INTEGER :: nofeas=0  !< flag for skipping making parameters feasible
    INTEGER :: nhistp=0  !< flag for histogram printout
    REAL    :: delfun=0.0!< expected function change
    REAL    :: actfun=0.0!< actual function change
    REAL    :: angras=0.0!< angle between gradient and search direction
    INTEGER :: iterat=0  !< iterations in solution
    INTEGER :: nregul=0  !< regularization flag
    REAL    :: regula=1.0!< regularization parameter, add regula * norm(global par.) to objective function
    REAL    :: regpre=0.0!< default presigma
    INTEGER :: matrit=0  !< matrix calculation up to iteration MATRIT
    INTEGER :: icalcm=0  !< calculation mode (for \ref xloopn "XLOOPN") , >0: calculate matrix
    INTEGER :: numbit=1  !< number of bits for pair counters
    INTEGER :: nbndr =0  !< number of records with bordered band matrix for local fit
    INTEGER :: nbdrx =0  !< max border size for local fit
    INTEGER :: nbndx =0  !< max band width for local fit
    INTEGER :: nrecer=0  !< record with error (rank deficit or Not-a-Number) for printout
    INTEGER :: nrec3 =maxi4 !< (1.) record number with error
    INTEGER :: mreqpe=1  !< min number of pair entries
    INTEGER :: mhispe=0  !< upper bound for pair entry histogrammimg
    INTEGER :: msngpe=0  !< upper bound for pair entry single precision storage
    INTEGER :: mcmprs=0  !< compression flag for sparsity (column indices)
    INTEGER :: mthrd =1  !< number of (OpenMP) threads
    INTEGER :: mxrec =0  !< max number of records
    INTEGER :: matmon=0  !< record interval for monitoring of (sparse) matrix construction
    INTEGER :: lfitnp=maxi4 !< local fit: number of iteration to calculate pulls
    INTEGER :: lfitbb=1  !< local fit: check for bordered band matrix (if >0)
    INTEGER :: mnrsel=0  !< number of MINRES error labels in LBMNRS (calc err, corr with SOLGLO)
    INTEGER :: ncache=-1 !< buffer size for caching (default 100MB per thread)
    REAL, DIMENSION(3) :: fcache = (/ 0.8,  0., 0. /) !< read cache, average fill level; write cache; dynamic size
    INTEGER :: mthrdr=1  !< number of threads for reading binary files
    INTEGER :: mnrsit=0  !< total number of MINRES internal iterations
    INTEGER :: iforce=0  !< switch to SUBITO for (global) rank defects if zero
    INTEGER :: igcorr=0  !< flag for output of global correlations for inversion, =0: none
    INTEGER :: memdbg=0  !< debug flag for memory management
    REAL    :: prange=0.0!< range (-PRANGE..PRANGE) for histograms of pulls, norm. residuals
    INTEGER :: lsearch=2 !< iterations (solutions) with line search:
                         !! >2: all, =2: all with (next) Chi2 cut scaling factor =1., =1: last, <1: none
    ! variables
    INTEGER :: lunlog !< unit for logfile
    INTEGER :: lvllog !< log level
    INTEGER :: ntgb !< total number of global parameters
    INTEGER :: nvgb !< number of variable global parameters
    INTEGER :: nagb !< number of fit parameters (global par. + Lagrange mult.)
    INTEGER :: ncgb !< number of constraints
    INTEGER :: nagbn !< max number of global paramters per record
    INTEGER :: nalcn !< max number of local paramters per record
    INTEGER :: naeqn !< max number of equations (measurements) per record
    INTEGER :: nrec  !< (current) record number
    REAL    :: dflim !< convergence limit
    INTEGER, DIMENSION(0:3) :: nrejec !< rejected events
    REAL, DIMENSION(0:8) :: times !< cpu time counters
    REAL    :: stepl !< step length (line search)
    CHARACTER (LEN=74) :: textl !< name of current MP 'module' (step)
    LOGICAL :: newite !< flag for new iteration
    INTEGER :: ndfsum !< sum(ndf)
    INTEGER :: iitera !< MINRES iterations
    INTEGER :: istopa !< MINRES istop (convergence)
    INTEGER :: lsinfo !< line search: returned information
    REAL    :: rstart !< cpu start time for solution iterations
    REAL    :: deltim !< cpu time difference
    INTEGER :: npresg !< number of pre-sigmas
    INTEGER :: nrecal !< number of records
    INTEGER :: nmiss1 !< rank deficit for constraints
    INTEGER :: lcalcm !< last calclation mode
    INTEGER :: nspc   !< number of precision for sparse global matrix (1=D, 2=D+F)
    INTEGER :: nencdb !< encoding info (number bits for column counter)
    INTEGER, DIMENSION(100) :: lbmnrs !< MINRES error labels
    DOUBLE PRECISION :: fvalue !< function value (chi2 sum) solution
    DOUBLE PRECISION :: flines !< function value line search
    DOUBLE PRECISION :: sumndf !< weighted sum(ndf)
    ! each loop
    INTEGER :: numReadbuffer     !< number of buffers (records) in (read) block
    INTEGER :: numBlocks         !< number of (read) blocks
    INTEGER :: sumRecords        !< sum of records
    INTEGER :: skippedRecords    !< number of skipped records (buffer too small)
    INTEGER :: minRecordsInBlock !< min. records in block
    INTEGER :: maxRecordsInBlock !< max. records in block
    ! accurate sumation
    INTEGER, PARAMETER::nexp20=1048576 ! 2**20
    DOUBLE PRECISION::accurateDsum=0.0D0 !< fractional part of sum
    INTEGER::accurateNsum=0 !< sum mod 2**20
    INTEGER::accurateNexp=0 !< sum  /  2**20
    INTEGER :: lenGlobalVec !< length of global vector 'b' (A*x=b)
    ! dynamic arrays
    !======================================================
    ! global parameters
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: globalParameter !< global parameters (start values + sum(x_i))
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: globalParCopy !< copy of global parameters
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: globalCorrections !< correction x_i (from A*x_i=b_i in iteration i)
    REAL, DIMENSION(:), ALLOCATABLE :: globalParStart     !< start value for global parameters
    REAL, DIMENSION(:), ALLOCATABLE :: globalParPreSigma  !< pre-sigma for global parameters
    REAL, DIMENSION(:), ALLOCATABLE :: globalParPreWeight !< weight from pre-sigma
    ! global matrix, vector
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: globalMatD !< global matrix 'A' (double, full or sparse)
    REAL, DIMENSION(:), ALLOCATABLE :: globalMatF !< global matrix 'A' (float part for compressed sparse)
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: globalVector !< global vector 'x' (in A*x=b)
    ! preconditioning
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: matPreCond !< preconditioner (band) matrix
    INTEGER, DIMENSION(:), ALLOCATABLE :: indPreCond !< preconditioner pointer array
    ! auxiliary vectors
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: workspaceD !< (general) workspace (D)
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: workspaceLinesearch !< workspace line search
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: workspaceDiagonalization !< workspace diag.
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: workspaceEigenValues !< workspace eigen values
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: workspaceEigenVectors !< workspace eigen vectors
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: workspaceMinres !< workspace MINRES
    INTEGER, DIMENSION(:), ALLOCATABLE :: workspaceI !< (general) workspace (I)
    ! constraint matrix, residuals
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: matConsProduct !< product matrix of constraints
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: vecConsResiduals !< residuals of constraints
    ! global parameter mapping
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: globalParLabelIndex !< global parameters label, total -> var. index
    INTEGER, DIMENSION(:), ALLOCATABLE :: globalParHashTable    !< global parameters hash table
    INTEGER, DIMENSION(:), ALLOCATABLE :: globalParVarToTotal   !< global parameters variable -> total index
    INTEGER, DIMENSION(-7:0) :: globalParHeader = 0 !< global parameters (mapping) header
                                                    !!
                                                    !!  0: length of labels/indices; \n
                                                    !! -1: number of stored items; \n
                                                    !! -2: =0 during build-up; \n
                                                    !! -3: next number; \n
                                                    !! -4: (largest) prime number (< length); \n
                                                    !! -5: number of overflows; \n
                                                    !! -6: nr of variable parameters; \n
                                                    !! -7: call counter for build-up;

    ! row information for sparse matrix
    INTEGER, DIMENSION(:), ALLOCATABLE :: sparseMatrixCompression !< compression info (per row)
    INTEGER, DIMENSION(:), ALLOCATABLE :: sparseMatrixColumns     !< (compressed) list of columns for sparse matrix
    INTEGER(kind=large), DIMENSION(:,:), ALLOCATABLE :: sparseMatrixOffsets !< row offsets for column list, sparse matrix elements
    ! read buffer
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: readBufferInfo !< buffer management (per thread)
    INTEGER, DIMENSION(:), ALLOCATABLE :: readBufferPointer !< pointer to used buffers
    INTEGER, DIMENSION(:), ALLOCATABLE :: readBufferDataI !< integer data
    REAL, DIMENSION(:), ALLOCATABLE :: readBufferDataF !< float data
    ! global parameter usage in record
    INTEGER, DIMENSION(:), ALLOCATABLE :: globalIndexUsage !< indices of global par in record
    INTEGER, DIMENSION(:), ALLOCATABLE :: backIndexUsage   !< list of global par in record
    ! local fit
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE::blvec  !< local fit vector 'b' (in A*x=b), replaced by 'x'
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE::clmat  !< local fit matrix 'A' (in A*x=b)
    INTEGER, DIMENSION(:), ALLOCATABLE:: ibandh !< local fit 'band width histogram' (band size autodetection)
    ! scratch arrays for local fit
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE::vbnd !< local fit band part of 'A'
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE::vbdr !< local fit border part of 'A'
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE::aux  !< local fit 'solutions for border rows'
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE::vbk  !< local fit 'matrix for border solution'
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE::vzru !< local fit 'border solution'
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE::scdiag !< local fit workspace (D)
    INTEGER, DIMENSION(:), ALLOCATABLE:: scflag         !< local fit workspace (I)
    REAL, DIMENSION(:), ALLOCATABLE :: localCorrections !< local fit corrections (to residuals)
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: localGlobalMatrix !< matrix correlating local and global par
    ! update of global matrix
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: writeBufferInfo  !< write buffer management (per thread)
    REAL, DIMENSION(:,:), ALLOCATABLE :: writeBufferData     !< write buffer data (largest residual, Chi2/ndf, per thread)
    INTEGER, DIMENSION(:), ALLOCATABLE :: writeBufferIndices !< write buffer for indices
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: writeBufferUpdates !< write buffer for update matrices
    INTEGER, DIMENSION(-6:6) :: writeBufferHeader = 0 !< write buffer header (-6..-1: updates, 1..6: indices)
                                                      !!
                                                      !! +/-1: buffer size (words) per thread; \n
                                                      !! +/-2: min number of free words; \n
                                                      !! +/-3: number of buffer flushes; \n
                                                      !! +/-4: number of buffer overruns; \n
                                                      !! +/-5: average fill level; \n
                                                      !! +/-6: peak fill level;
    !> list items from steering file
    INTEGER :: lenParameters=0   !< length of list of parameters from steering file
    TYPE(listItem), DIMENSION(:), ALLOCATABLE :: listParameters   !< list of parameters from steering file
    INTEGER :: lenPresigmas=0    !< length of list of pre-sigmas from steering file
    TYPE(listItem), DIMENSION(:), ALLOCATABLE :: listPreSigmas    !< list of pre-sgmas from steering file
    INTEGER :: lenConstraints=0  !< length of list of constraints from steering file
    TYPE(listItem), DIMENSION(:), ALLOCATABLE :: listConstraints  !< list of constraints from steering file
    INTEGER :: lenMeasurements=0 !< length of list of measurements from steering file
    TYPE(listItem), DIMENSION(:), ALLOCATABLE :: listMeasurements !< list of measurements from steering file
    !======================================================
    ! file information
    INTEGER, DIMENSION(:), ALLOCATABLE :: mfd   !< file mode: cbinary =1, text =2, fbinary=3
    INTEGER, DIMENSION(:), ALLOCATABLE :: lfd   !< length of file name
    INTEGER, DIMENSION(:), ALLOCATABLE :: nfd   !< index (line) in (steering) file
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: kfd !< (1,.)=  number of records in file, (2,..)= file order
    INTEGER, DIMENSION(:), ALLOCATABLE :: ifd   !< file: integrated record numbers (=offset)
    INTEGER, DIMENSION(:), ALLOCATABLE :: jfd   !< file: number of accepted records
    INTEGER, DIMENSION(:), ALLOCATABLE :: dfd   !< file: ndf sum
    INTEGER, DIMENSION(:), ALLOCATABLE :: xfd   !< file: max. record size
    REAL, DIMENSION(:), ALLOCATABLE :: cfd      !< file: chi2 sum
    REAL, DIMENSION(:), ALLOCATABLE :: ofd      !< file: option
    REAL, DIMENSION(:), ALLOCATABLE :: wfd      !< file: weight
    CHARACTER (LEN=1024) :: filnam !< name of steering file
    INTEGER :: nfnam  !< length of sterring file name
    CHARACTER, DIMENSION(:), ALLOCATABLE :: tfd !< file names (concatenation)
    INTEGER :: ifile  !< current file (index)
    INTEGER :: nfiles !< number of files
    INTEGER :: nfilb  !< number of binary files
    INTEGER :: nfilf  !< number of Fortran binary files
    INTEGER :: nfilc  !< number of C binary files
    INTEGER :: nfilw  !< number of weighted binary files
    INTEGER :: ndimbuf=10000 !< default read buffer size (I/F words, half record length)

END MODULE mpmod
