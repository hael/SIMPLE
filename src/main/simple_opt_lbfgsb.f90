! memory limited bound-constained BFGS optimiser
module simple_opt_lbfgsb
include 'simple_lib.f08'
use simple_optimizer, only: optimizer
use simple_opt_helpers
implicit none

public :: opt_lbfgsb
private

real(dp), parameter :: zero=0.0d0,one=1.0d0,p5=0.5d0,p66=0.66d0
real(dp), parameter :: xtrapl=1.1d0,xtrapu=4.0d0
real(dp), parameter :: big=1.0d+10
real(dp), parameter :: ftol_lbfgsb=1.0d-3,gtol_lbfgsb=0.9d0,xtol_lbfgsb=0.1d0
real(dp), parameter :: two=2.0d0,three=3.0d0
real(dp), parameter :: factr  = 1.0d+7, pgtol  = 1.0d-5

type, extends(optimizer) :: opt_lbfgsb
    private
    character(len=60)      :: task, csave
    logical                :: lsave(4)
    integer                :: isave(44)
    real(dp)               :: f
    real(dp)               :: dsave(24)
    integer,  allocatable  :: nbd(:), iwa(:)
    real(dp), allocatable  ::  l(:), u(:), g(:), wa(:)
    logical                :: exists=.false.
contains
    procedure :: new      => new_opt_lbfgsb
    procedure :: minimize => lbfgsb_minimize
    procedure :: kill     => kill_opt_lbfgsb
end type opt_lbfgsb

contains

    !> \brief  is a constructor
    subroutine new_opt_lbfgsb( self, spec )
        use simple_opt_spec, only: opt_spec
        class(opt_lbfgsb), intent(inout) :: self !< instance
        class(opt_spec),   intent(inout) :: spec !< specification
        integer                          :: n, m
        call self%kill
        n = spec%ndim
        m = spec%lbfgsb_m
        allocate ( self%nbd(n), self%l(n), self%u(n), self%g(n) )
        allocate ( self%iwa(3*n) )
        allocate ( self%wa(2*m*n + 5*n + 11*m*m + 8*m) )
        self%nbd = 2
        self%exists = .true.
    end subroutine new_opt_lbfgsb

    !>  \brief  nonlinear conjugate gradient minimizer
    subroutine lbfgsb_minimize( self, spec, fun_self, lowest_cost )
        use simple_opt_spec, only: opt_spec
        use simple_opt_subs, only: lnsrch
        class(opt_lbfgsb), intent(inout) :: self        !< instance
        class(opt_spec),   intent(inout) :: spec        !< specification
        class(*),          intent(inout) :: fun_self    !< self-pointer for cost function
        real,              intent(out)   :: lowest_cost !< minimum function value
        if ( (.not. associated(spec%costfun) ).and.( .not. associated(spec%costfun_8)) ) then
            stop 'cost function not associated in opt_spec; lbfgsb_minimize; simple_opt_lbfgsb'
        endif
        if ( (.not. associated(spec%gcostfun) ).and.(.not. associated(spec%gcostfun_8)) ) then
            stop 'gradient of cost function not associated in opt_spec; lbfgsb_minimize; simple_opt_lbfgsb'
        endif
        if ( .not. allocated(spec%limits) ) then
            stop 'limits not allocated in opt_spec; lbfgsb_minimize; simple_opt_lbfgsb'
        end if
        spec%x_8 = spec%x
        ! initialise nevals counters
        spec%nevals  = 0
        spec%ngevals = 0

        self%l = spec%limits(:,1)
        self%u = spec%limits(:,2)
        call lbfgsbmin

    contains

        !>  \brief  nonlinear conjugate gradient minimizer
        subroutine lbfgsbmin
            self%task = 'START'
            do while(self%task(1:2) .eq. 'FG' .or. self%task .eq. 'NEW_X' .or. &
                &    self%task .eq. 'START')
                !     This is the call to the L-BFGS-B code.
                call setulb ( spec%ndim, spec%lbfgsb_m, spec%x_8, self%l, self%u, self%nbd, self%f, self%g, &
                    &         factr, pgtol, self%wa, self%iwa, self%task,  &
                    &         self%csave, self%lsave, self%isave, self%dsave      )

                if (self%task(1:2) .eq. 'FG') then
                    call spec%eval_fdf_8( fun_self, spec%x_8, self%f, self%g )
                else if (self%task(1:5) .eq. 'NEW_X') then
                    ! the minimization routine has returned with a new iterate.
                    ! At this point have the opportunity of stopping the iteration
                    ! or observing the values of certain parameters
                    if (associated(spec%opt_callback)) then
                        spec%x = real(spec%x_8)
                        call spec%opt_callback(fun_self)
                    end if
                    if (self%isave(34) .gt. spec%maxits)  then
                        self%task='STOP: TOTAL NO. of f AND g EVALUATIONS EXCEEDS LIMIT'
                    end if
                end if
            end do
            lowest_cost = real(self%f)
            spec%x      = real(spec%x_8)
        end subroutine lbfgsbmin

        !===========   L-BFGS-B (version 3.0.  April 25, 2011  ===================
        !
        !     This is a modified version of L-BFGS-B. Minor changes in the updated
        !     code appear preceded by a line comment as follows
        !
        !     c-jlm-jn
        !
        !     Major changes are described in the accompanying paper:
        !
        !         Jorge Nocedal and Jose Luis Morales, Remark on "Algorithm 778:
        !         L-BFGS-B: Fortran Subroutines for Large-Scale Bound Constrained
        !         Optimization"  (2011). To appear in  ACM Transactions on
        !         Mathematical Software,
        !
        !     The paper describes an improvement and a correction to Algorithm 778.
        !     It is shown that the performance of the algorithm can be improved
        !     significantly by making a relatively simple modication to the subspace
        !     minimization phase. The correction concerns an error caused by the use
        !     of routine dpmeps to estimate machine precision.
        !
        !     The total work space **wa** required by the new version is
        !
        !                  2*m*n + 11m*m + 5*n + 8*m
        !
        !     the old version required
        !
        !                  2*m*n + 12m*m + 4*n + 12*m
        !
        !
        !            J. Nocedal  Department of Electrical Engineering and
        !                        Computer Science.
        !                        Northwestern University. Evanston, IL. USA
        !
        !
        !           J.L Morales  Departamento de Matematicas,
        !                        Instituto Tecnologico Autonomo de Mexico
        !                        Mexico D.F. Mexico.
        !
        !                        March  2011
        !
        !=============================================================================
        subroutine setulb(n, m, x, l, u, nbd, f, g, factr, pgtol, wa, iwa, &
            &             task, csave, lsave, isave, dsave)
            character :: task*60, csave*60
            logical   :: lsave(4)
            integer   :: n, m
            integer   :: nbd(n), iwa(3*n), isave(44)
            real(dp)  :: f, factr, pgtol, x(n), l(n), u(n), g(n)
            real(dp)  :: wa(2*m*n + 5*n + 11*m*m + 8*m), dsave(24)
            !     ************
            !
            !     Subroutine setulb
            !
            !     This subroutine partitions the working arrays wa and iwa, and
            !       then uses the limited memory BFGS method to solve the bound
            !       constrained optimization problem by calling mainlb.
            !       (The direct method will be used in the subspace minimization.)
            !
            !     n is an integer variable.
            !       On entry n is the dimension of the problem.
            !       On exit n is unchanged.
            !
            !     m is an integer variable.
            !       On entry m is the maximum number of variable metric corrections
            !         used to define the limited memory matrix.
            !       On exit m is unchanged.
            !
            !     x is a double precision array of dimension n.
            !       On entry x is an approximation to the solution.
            !       On exit x is the current approximation.
            !
            !     l is a double precision array of dimension n.
            !       On entry l is the lower bound on x.
            !       On exit l is unchanged.
            !
            !     u is a double precision array of dimension n.
            !       On entry u is the upper bound on x.
            !       On exit u is unchanged.
            !
            !     nbd is an integer array of dimension n.
            !       On entry nbd represents the type of bounds imposed on the
            !         variables, and must be specified as follows:
            !         nbd(i)=0 if x(i) is unbounded,
            !                1 if x(i) has only a lower bound,
            !                2 if x(i) has both lower and upper bounds, and
            !                3 if x(i) has only an upper bound.
            !       On exit nbd is unchanged.
            !
            !     f is a double precision variable.
            !       On first entry f is unspecified.
            !       On final exit f is the value of the function at x.
            !
            !     g is a double precision array of dimension n.
            !       On first entry g is unspecified.
            !       On final exit g is the value of the gradient at x.
            !
            !     factr is a double precision variable.
            !       On entry factr >= 0 is specified by the user.  The iteration
            !         will stop when
            !
            !         (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch
            !
            !         where epsmch is the machine precision, which is automatically
            !         generated by the code. Typical values for factr: 1.d+12 for
            !         low accuracy; 1.d+7 for moderate accuracy; 1.d+1 for extremely
            !         high accuracy.
            !       On exit factr is unchanged.
            !
            !     pgtol is a double precision variable.
            !       On entry pgtol >= 0 is specified by the user.  The iteration
            !         will stop when
            !
            !                 max{|proj g_i | i = 1, ..., n} <= pgtol
            !
            !         where pg_i is the ith component of the projected gradient.
            !       On exit pgtol is unchanged.
            !
            !     wa is a double precision working array of length
            !       (2mmax + 5)nmax + 12mmax^2 + 12mmax.
            !
            !     iwa is an integer working array of length 3nmax.
            !
            !     task is a working string of characters of length 60 indicating
            !       the current job when entering and quitting this subroutine.
            !
            !     csave is a working string of characters of length 60.
            !
            !     lsave is a logical working array of dimension 4.
            !       On exit with 'task' = NEW_X, the following information is
            !                                                             available:
            !         If lsave(1) = .true.  then  the initial X has been replaced by
            !                                     its projection in the feasible set;
            !         If lsave(2) = .true.  then  the problem is constrained;
            !         If lsave(3) = .true.  then  each variable has upper and lower
            !                                     bounds;
            !
            !     isave is an integer working array of dimension 44.
            !       On exit with 'task' = NEW_X, the following information is
            !                                                             available:
            !         isave(22) = the total number of intervals explored in the
            !                         search of Cauchy points;
            !         isave(26) = the total number of skipped BFGS updates before
            !                         the current iteration;
            !         isave(30) = the number of current iteration;
            !         isave(31) = the total number of BFGS updates prior the current
            !                         iteration;
            !         isave(33) = the number of intervals explored in the search of
            !                         Cauchy point in the current iteration;
            !         isave(34) = the total number of function and gradient
            !                         evaluations;
            !         isave(36) = the number of function value or gradient
            !                                  evaluations in the current iteration;
            !         if isave(37) = 0  then the subspace argmin is within the box;
            !         if isave(37) = 1  then the subspace argmin is beyond the box;
            !         isave(38) = the number of free variables in the current
            !                         iteration;
            !         isave(39) = the number of active constraints in the current
            !                         iteration;
            !         n + 1 - isave(40) = the number of variables leaving the set of
            !                           active constraints in the current iteration;
            !         isave(41) = the number of variables entering the set of active
            !                         constraints in the current iteration.
            !
            !     dsave is a double precision working array of dimension 29.
            !       On exit with 'task' = NEW_X, the following information is
            !                                                             available:
            !         dsave(1) = current 'theta' in the BFGS matrix;
            !         dsave(2) = f(x) in the previous iteration;
            !         dsave(3) = factr*epsmch;
            !         dsave(4) = 2-norm of the line search direction vector;
            !         dsave(5) = the machine precision epsmch generated by the code;
            !         dsave(6) = the slope of the line search function at
            !                                  the current point of line search;
            !         dsave(7) = the maximum relative step length imposed in
            !                                                           line search;
            !         dsave(8) = the infinity norm of the projected gradient;
            !         dsave(9) = the relative step length in the line search;
            !         dsave(10) = the slope of the line search function at
            !                                 the starting point of the line search;
            !         dsave(11) = the square of the 2-norm of the line search
            !                                                      direction vector.
            !
            !     Subprograms called:
            !
            !       L-BFGS-B Library ... mainlb.
            !
            !
            !     References:
            !
            !       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
            !       memory algorithm for bound constrained optimization'',
            !       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
            !
            !       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: a
            !       limited memory FORTRAN code for solving bound constrained
            !       optimization problems'', Tech. Report, NAM-11, EECS Department,
            !       Northwestern University, 1994.
            !
            !       (Postscript files of these papers are available via anonymous
            !        ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
            !
            !                           *  *  *
            !
            !     NEOS, November 1994. (Latest revision June 1996.)
            !     Optimization Technology Center.
            !     Argonne National Laboratory and Northwestern University.
            !     Written by
            !                        Ciyou Zhu
            !     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
            !
            !
            !     ************
            integer :: lws,lr,lz,lt,ld,lxp,lwa
            integer :: lwy,lsy,lss,lwt,lwn,lsnd
            if (task .eq. 'START') then
                isave(1)  = m*n
                isave(2)  = m**2
                isave(3)  = 4*m**2
                isave(4)  = 1                      ! ws      m*n
                isave(5)  = isave(4)  + isave(1)   ! wy      m*n
                isave(6)  = isave(5)  + isave(1)   ! wsy     m**2
                isave(7)  = isave(6)  + isave(2)   ! wss     m**2
                isave(8)  = isave(7)  + isave(2)   ! wt      m**2
                isave(9)  = isave(8)  + isave(2)   ! wn      4*m**2
                isave(10) = isave(9)  + isave(3)   ! wsnd    4*m**2
                isave(11) = isave(10) + isave(3)   ! wz      n
                isave(12) = isave(11) + n          ! wr      n
                isave(13) = isave(12) + n          ! wd      n
                isave(14) = isave(13) + n          ! wt      n
                isave(15) = isave(14) + n          ! wxp     n
                isave(16) = isave(15) + n          ! wa      8*m
            end if
            lws  = isave(4)
            lwy  = isave(5)
            lsy  = isave(6)
            lss  = isave(7)
            lwt  = isave(8)
            lwn  = isave(9)
            lsnd = isave(10)
            lz   = isave(11)
            lr   = isave(12)
            ld   = isave(13)
            lt   = isave(14)
            lxp  = isave(15)
            lwa  = isave(16)
            call mainlb(n,m,x,l,u,nbd,f,g,factr,pgtol, &
                &  wa(lws),wa(lwy),wa(lsy),wa(lss), wa(lwt), &
                &  wa(lwn),wa(lsnd),wa(lz),wa(lr),wa(ld),wa(lt),wa(lxp), &
                &  wa(lwa), &
                &  iwa(1),iwa(n+1),iwa(2*n+1),task, &
                &  csave,lsave,isave(22),dsave)
        end subroutine setulb

        subroutine mainlb(n, m, x, l, u, nbd, f, g, factr, pgtol, ws, wy, &
            &                  sy, ss, wt, wn, snd, z, r, d, t, xp, wa, &
            &                  index, iwhere, indx2, task, &
            &                  csave, lsave, isave, dsave)
            implicit none
            character(len=60) :: task, csave
            logical   :: lsave(4)
            integer   :: n, m, nbd(n), index(n)
            integer   :: iwhere(n), indx2(n), isave(23)
            real(dp)  :: f, factr, pgtol
            real(dp)  :: x(n), l(n), u(n), g(n), z(n), r(n), d(n), t(n)
            real(dp)  :: xp(n)
            real(dp)  :: wa(8*m)
            real(dp)  :: ws(n, m), wy(n, m), sy(m, m), ss(m, m)
            real(dp)  :: wt(m, m), wn(2*m, 2*m), snd(2*m, 2*m), dsave(24)
            !     ************
            !
            !     Subroutine mainlb
            !
            !     This subroutine solves bound constrained optimization problems by
            !       using the compact formula of the limited memory BFGS updates.
            !
            !     n is an integer variable.
            !       On entry n is the number of variables.
            !       On exit n is unchanged.
            !
            !     m is an integer variable.
            !       On entry m is the maximum number of variable metric
            !          corrections allowed in the limited memory matrix.
            !       On exit m is unchanged.
            !
            !     x is a double precision array of dimension n.
            !       On entry x is an approximation to the solution.
            !       On exit x is the current approximation.
            !
            !     l is a double precision array of dimension n.
            !       On entry l is the lower bound of x.
            !       On exit l is unchanged.
            !
            !     u is a double precision array of dimension n.
            !       On entry u is the upper bound of x.
            !       On exit u is unchanged.
            !
            !     nbd is an integer array of dimension n.
            !       On entry nbd represents the type of bounds imposed on the
            !         variables, and must be specified as follows:
            !         nbd(i)=0 if x(i) is unbounded,
            !                1 if x(i) has only a lower bound,
            !                2 if x(i) has both lower and upper bounds,
            !                3 if x(i) has only an upper bound.
            !       On exit nbd is unchanged.
            !
            !     f is a double precision variable.
            !       On first entry f is unspecified.
            !       On final exit f is the value of the function at x.
            !
            !     g is a double precision array of dimension n.
            !       On first entry g is unspecified.
            !       On final exit g is the value of the gradient at x.
            !
            !     factr is a double precision variable.
            !       On entry factr >= 0 is specified by the user.  The iteration
            !         will stop when
            !
            !         (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch
            !
            !         where epsmch is the machine precision, which is automatically
            !         generated by the code.
            !       On exit factr is unchanged.
            !
            !     pgtol is a double precision variable.
            !       On entry pgtol >= 0 is specified by the user.  The iteration
            !         will stop when
            !
            !                 max{|proj g_i | i = 1, ..., n} <= pgtol
            !
            !         where pg_i is the ith component of the projected gradient.
            !       On exit pgtol is unchanged.
            !
            !     ws, wy, sy, and wt are double precision working arrays used to
            !       store the following information defining the limited memory
            !          BFGS matrix:
            !          ws, of dimension n x m, stores S, the matrix of s-vectors;
            !          wy, of dimension n x m, stores Y, the matrix of y-vectors;
            !          sy, of dimension m x m, stores S'Y;
            !          ss, of dimension m x m, stores S'S;
            !          yy, of dimension m x m, stores Y'Y;
            !          wt, of dimension m x m, stores the Cholesky factorization
            !                                  of (theta*S'S+LD^(-1)L'); see eq.
            !                                  (2.26) in [3].
            !
            !     wn is a double precision working array of dimension 2m x 2m
            !       used to store the LEL^T factorization of the indefinite matrix
            !                 K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
            !                     [L_a -R_z           theta*S'AA'S ]
            !
            !       where     E = [-I  0]
            !                     [ 0  I]
            !
            !     snd is a double precision working array of dimension 2m x 2m
            !       used to store the lower triangular part of
            !                 N = [Y' ZZ'Y   L_a'+R_z']
            !                     [L_a +R_z  S'AA'S   ]
            !
            !     z(n),r(n),d(n),t(n), xp(n),wa(8*m) are double precision working arrays.
            !       z  is used at different times to store the Cauchy point and
            !          the Newton point.
            !       xp is used to safeguard the projected Newton direction
            !
            !     sg(m),sgo(m),yg(m),ygo(m) are double precision working arrays.
            !
            !     index is an integer working array of dimension n.
            !       In subroutine freev, index is used to store the free and fixed
            !          variables at the Generalized Cauchy Point (GCP).
            !
            !     iwhere is an integer working array of dimension n used to record
            !       the status of the vector x for GCP computation.
            !       iwhere(i)=0 or -3 if x(i) is free and has bounds,
            !                 1       if x(i) is fixed at l(i), and l(i) .ne. u(i)
            !                 2       if x(i) is fixed at u(i), and u(i) .ne. l(i)
            !                 3       if x(i) is always fixed, i.e.,  u(i)=x(i)=l(i)
            !                -1       if x(i) is always free, i.e., no bounds on it.
            !
            !     indx2 is an integer working array of dimension n.
            !       Within subroutine cauchy, indx2 corresponds to the array iorder.
            !       In subroutine freev, a list of variables entering and leaving
            !       the free set is stored in indx2, and it is passed on to
            !       subroutine formk with this information.
            !
            !     task is a working string of characters of length 60 indicating
            !       the current job when entering and leaving this subroutine.
            !
            !     csave is a working string of characters of length 60.
            !
            !     lsave is a logical working array of dimension 4.
            !
            !     isave is an integer working array of dimension 23.
            !
            !     dsave is a double precision working array of dimension 29.
            !
            !
            !     Subprograms called
            !
            !       L-BFGS-B Library ... cauchy, subsm, lnsrlb, formk,
            !
            !        errclb, active, projgr,
            !
            !        freev, cmprlb, matupd, formt.
            !
            !       Linpack Library ... dcopy, ddot.
            !
            !
            !     References:
            !
            !       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
            !       memory algorithm for bound constrained optimization'',
            !       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
            !
            !       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: FORTRAN
            !       Subroutines for Large Scale Bound Constrained Optimization''
            !       Tech. Report, NAM-11, EECS Department, Northwestern University,
            !       1994.
            !
            !       [3] R. Byrd, J. Nocedal and R. Schnabel "Representations of
            !       Quasi-Newton Matrices and their use in Limited Memory Methods'',
            !       Mathematical Programming 63 (1994), no. 4, pp. 129-156.
            !
            !       (Postscript files of these papers are available via anonymous
            !        ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
            !
            !                           *  *  *
            !
            !     NEOS, November 1994. (Latest revision June 1996.)
            !     Optimization Technology Center.
            !     Argonne National Laboratory and Northwestern University.
            !     Written by
            !                        Ciyou Zhu
            !     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
            !
            !
            !     ************
            logical   :: prjctd,cnstnd,boxed,updatd,wrk
            character(len=3) :: word
            integer   :: i,k,nintol,itfile,iback,nskip
            integer   :: head,col,iter,itail,iupdat
            integer   :: nseg,nfgv,info,ifun
            integer   :: iword,nfree,nact,ileave,nenter
            real(dp)  :: theta,fold,dr,rr,tol
            real(dp)  :: xstep,sbgnrm,ddum,dnorm,dtd,epsmch
            real(dp)  :: gd,gdold,stp,stpmx
            if (task .eq. 'START') then
                epsmch = epsilon(one)
                ! Initialize counters and scalars when task='START'.
                ! for the limited memory BFGS matrices:
                col    = 0
                head   = 1
                theta  = one
                iupdat = 0
                updatd = .false.
                iback  = 0
                itail  = 0
                iword  = 0
                nact   = 0
                ileave = 0
                nenter = 0
                fold   = zero
                dnorm  = zero
                gd     = zero
                stpmx  = zero
                sbgnrm = zero
                stp    = zero
                gdold  = zero
                dtd    = zero
                ! for operation counts:
                iter   = 0
                nfgv   = 0
                nseg   = 0
                nintol = 0
                nskip  = 0
                nfree  = n
                ifun   = 0
                ! for stopping tolerance:
                tol    = factr*epsmch
                ! 'word' records the status of subspace solutions.
                word   = '---'
                ! 'info' records the termination information.
                info   = 0
                itfile = 8
                !  Check the input arguments for errors.
                call errclb(n,m,factr,l,u,nbd,task,info,k)
                if (task(1:5) .eq. 'ERROR') return
                ! Initialize iwhere & project x onto the feasible set.
                call active(n,l,u,nbd,x,iwhere,prjctd,cnstnd,boxed)
                ! The end of the initialization.
            else
                ! restore local variables.
                prjctd = lsave(1)
                cnstnd = lsave(2)
                boxed  = lsave(3)
                updatd = lsave(4)
                nintol = isave(1)
                itfile = isave(3)
                iback  = isave(4)
                nskip  = isave(5)
                head   = isave(6)
                col    = isave(7)
                itail  = isave(8)
                iter   = isave(9)
                iupdat = isave(10)
                nseg   = isave(12)
                nfgv   = isave(13)
                info   = isave(14)
                ifun   = isave(15)
                iword  = isave(16)
                nfree  = isave(17)
                nact   = isave(18)
                ileave = isave(19)
                nenter = isave(20)
                theta  = dsave(1)
                fold   = dsave(2)
                tol    = dsave(3)
                dnorm  = dsave(4)
                epsmch = dsave(5)
                gd     = dsave(6)
                stpmx  = dsave(7)
                sbgnrm = dsave(8)
                stp    = dsave(9)
                gdold  = dsave(10)
                dtd    = dsave(11)
                ! After returning from the driver go to the point where execution is to resume.
                if (task(1:5) .eq. 'FG_LN') goto 666
                if (task(1:5) .eq. 'NEW_X') goto 777
                if (task(1:5) .eq. 'FG_ST') goto 111
                if (task(1:4) .eq. 'STOP') then
                    if (task(7:9) .eq. 'CPU') then
                        ! restore the previous iterate.
                        call dcopy(n,t,1,x,1)
                        call dcopy(n,r,1,g,1)
                        f = fold
                    end if
                    goto 999
                end if
            end if
            !  Compute f0 and g0.
            task = 'FG_START'
            ! return to the driver to calculate f and g; reenter at 111.
            goto 1000
111         continue
            nfgv = 1
            ! Compute the infinity norm of the (-) projected gradient.
            call projgr(n,l,u,nbd,x,g,sbgnrm)
            if (sbgnrm .le. pgtol) then
                !  terminate the algorithm.
                task = 'CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL'
                goto 999
            end if
            ! ----------------- the beginning of the loop --------------------------
222         continue
            iword = -1
            if (.not. cnstnd .and. col .gt. 0) then
                ! skip the search for GCP.
                call dcopy(n,x,1,z,1)
                wrk = updatd
                nseg = 0
                goto 333
            end if
            !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            !
            !     Compute the Generalized Cauchy Point (GCP).
            !
            !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            call cauchy(n,x,l,u,nbd,g,indx2,iwhere,t,d,z, &
                &            m,wy,ws,sy,wt,theta,col,head, &
                &            wa(1),wa(2*m+1),wa(4*m+1),wa(6*m+1),nseg, &
                &            sbgnrm, info, epsmch)
            if (info .ne. 0) then
                ! singular triangular system detected; refresh the lbfgs memory.
                info   = 0
                col    = 0
                head   = 1
                theta  = one
                iupdat = 0
                updatd = .false.
                goto 222
            end if
            nintol = nintol + nseg
            ! Count the entering and leaving variables for iter > 0;
            ! find the index set of free and active variables at the GCP.
            call freev(n,nfree,index,nenter,ileave,indx2, &
                &           iwhere,wrk,updatd,cnstnd,iter)
            nact = n - nfree
333         continue
            ! If there are no free variables or B=theta*I, then
            ! skip the subspace minimization.
            if (nfree .eq. 0 .or. col .eq. 0) goto 555
            !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            !
            !     Subspace minimization.
            !
            !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            !     Form  the LEL^T factorization of the indefinite
            !       matrix    K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
            !                     [L_a -R_z           theta*S'AA'S ]
            !       where     E = [-I  0]
            !                     [ 0  I]
            if (wrk) call formk(n,nfree,index,nenter,ileave,indx2,iupdat, &
                &                 updatd,wn,snd,m,ws,wy,sy,theta,col,head,info)
            if (info .ne. 0) then
                ! nonpositive definiteness in Cholesky factorization;
                ! refresh the lbfgs memory and restart the iteration.
                info   = 0
                col    = 0
                head   = 1
                theta  = one
                iupdat = 0
                updatd = .false.
                goto 222
            end if
            ! compute r=-Z'B(xcp-xk)-Z'g (using wa(2m+1)=W'(xcp-x) from 'cauchy').
            call cmprlb(n,m,x,g,ws,wy,sy,wt,z,r,wa,index, &
                &           theta,col,head,nfree,cnstnd,info)
            if (info .ne. 0) goto 444
            !-jlm-jn call the direct method.
            call subsm( n, m, nfree, index, l, u, nbd, z, r, xp, ws, wy, &
                &           theta, x, g, col, head, iword, wa, wn, info)
444         continue
            if (info .ne. 0) then
                ! singular triangular system detected;
                ! refresh the lbfgs memory and restart the iteration.
                info   = 0
                col    = 0
                head   = 1
                theta  = one
                iupdat = 0
                updatd = .false.
                goto 222
            end if
555         continue
            !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            !
            !     Line search and optimality tests.
            !
            !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            !     Generate the search direction d:=z-x.
            do i = 1, n
                d(i) = z(i) - x(i)
            end do
666         continue
            call lnsrlb(n,l,u,nbd,x,f,fold,gd,gdold,g,d,r,t,z,stp,dnorm, &
                &            dtd,xstep,stpmx,iter,ifun,iback,nfgv,info,task, &
                &            boxed,cnstnd,csave,isave(22),dsave(12))
            if (info .ne. 0 .or. iback .ge. 20) then
                ! restore the previous iterate.
                call dcopy(n,t,1,x,1)
                call dcopy(n,r,1,g,1)
                f = fold
                if (col .eq. 0) then
                    ! abnormal termination.
                    if (info .eq. 0) then
                        info = -9
                        ! restore the actual number of f and g evaluations etc.
                        nfgv = nfgv - 1
                        ifun = ifun - 1
                        iback = iback - 1
                    end if
                    task = 'ABNORMAL_TERMINATION_IN_LNSRCH'
                    iter = iter + 1
                    goto 999
                else
                    ! refresh the lbfgs memory and restart the iteration.
                    if (info .eq. 0) nfgv = nfgv - 1
                    info   = 0
                    col    = 0
                    head   = 1
                    theta  = one
                    iupdat = 0
                    updatd = .false.
                    task   = 'RESTART_FROM_LNSRCH'
                    goto 222
                end if
            else if (task(1:5) .eq. 'FG_LN') then
                ! return to the driver for calculating f and g; reenter at 666.
                goto 1000
            else
                ! calculate and print out the quantities related to the new X.
                iter = iter + 1
                ! Compute the infinity norm of the projected (-)gradient.
                call projgr(n,l,u,nbd,x,g,sbgnrm)
                ! Print iteration information.
                goto 1000
            end if
777         continue
            ! Test for termination.
            if (sbgnrm .le. pgtol) then
                !  terminate the algorithm.
                task = 'CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL'
                goto 999
            end if
            ddum = max(abs(fold), abs(f), one)
            if ((fold - f) .le. tol*ddum) then
                ! terminate the algorithm.
                task = 'CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH'
                if (iback .ge. 10) info = -5
                ! i.e., to issue a warning if iback>10 in the line search.
                goto 999
            end if
            ! Compute d=newx-oldx, r=newg-oldg, rr=y'y and dr=y's.
            do i = 1, n
                r(i) = g(i) - r(i)
            end do
            rr = ddot_lbfgsb(n,r,1,r,1)
            if (is_equal(stp , one)) then ! stp.eq.one
                dr = gd - gdold
                ddum = -gdold
            else
                dr = (gd - gdold)*stp
                call dscal(n,stp,d,1)
                ddum = -gdold*stp
            end if
            if (dr .le. epsmch*ddum) then
                ! skip the L-BFGS update.
                nskip = nskip + 1
                updatd = .false.
                goto 888
            end if
            !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            !
            !     Update the L-BFGS matrix.
            !
            !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            updatd = .true.
            iupdat = iupdat + 1
            ! Update matrices WS and WY and form the middle matrix in B.
            call matupd(n,m,ws,wy,sy,ss,d,r,itail, &
                &            iupdat,col,head,theta,rr,dr,stp,dtd)
            ! Form the upper half of the pds T = theta*SS + L*D^(-1)*L';
            ! Store T in the upper triangular of the array wt;
            ! Cholesky factorize T to J*J' with
            ! J' stored in the upper triangular of wt.
            call formt(m,wt,sy,ss,col,theta,info)
            if (info .ne. 0) then
                ! nonpositive definiteness in Cholesky factorization;
                ! refresh the lbfgs memory and restart the iteration.
                info = 0
                col = 0
                head = 1
                theta = one
                iupdat = 0
                updatd = .false.
                goto 222
            end if
            !     Now the inverse of the middle matrix in B is
            !       [  D^(1/2)      O ] [ -D^(1/2)  D^(-1/2)*L' ]
            !       [ -L*D^(-1/2)   J ] [  0        J'          ]
888         continue
            ! -------------------- the end of the loop -----------------------------
            goto 222
999         continue
1000        continue
            ! Save local variables.
            lsave(1)  = prjctd
            lsave(2)  = cnstnd
            lsave(3)  = boxed
            lsave(4)  = updatd
            isave(1)  = nintol
            isave(3)  = itfile
            isave(4)  = iback
            isave(5)  = nskip
            isave(6)  = head
            isave(7)  = col
            isave(8)  = itail
            isave(9)  = iter
            isave(10) = iupdat
            isave(12) = nseg
            isave(13) = nfgv
            isave(14) = info
            isave(15) = ifun
            isave(16) = iword
            isave(17) = nfree
            isave(18) = nact
            isave(19) = ileave
            isave(20) = nenter
            dsave(1)  = theta
            dsave(2)  = fold
            dsave(3)  = tol
            dsave(4)  = dnorm
            dsave(5)  = epsmch
            dsave(6)  = gd
            dsave(7)  = stpmx
            dsave(8)  = sbgnrm
            dsave(9)  = stp
            dsave(10) = gdold
            dsave(11) = dtd
        end subroutine mainlb

        subroutine active(n, l, u, nbd, x, iwhere, prjctd, cnstnd, boxed)
            logical  :: prjctd, cnstnd, boxed
            integer  :: n, nbd(n), iwhere(n)
            real(dp) :: x(n), l(n), u(n)
            !     ************
            !
            !     Subroutine active
            !
            !     This subroutine initializes iwhere and projects the initial x to
            !       the feasible set if necessary.
            !
            !     iwhere is an integer array of dimension n.
            !       On entry iwhere is unspecified.
            !       On exit iwhere(i)=-1  if x(i) has no bounds
            !                         3   if l(i)=u(i)
            !                         0   otherwise.
            !       In cauchy, iwhere is given finer gradations.
            !
            !
            !                           *  *  *
            !
            !     NEOS, November 1994. (Latest revision June 1996.)
            !     Optimization Technology Center.
            !     Argonne National Laboratory and Northwestern University.
            !     Written by
            !                        Ciyou Zhu
            !     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
            !
            !
            !     ************
            integer :: nbdd,i
            ! Initialize nbdd, prjctd, cnstnd and boxed.
            nbdd = 0
            prjctd = .false.
            cnstnd = .false.
            boxed = .true.
            ! Project the initial x to the easible set if necessary.
            do i = 1, n
                if (nbd(i) .gt. 0) then
                    if (nbd(i) .le. 2 .and. x(i) .le. l(i)) then
                        if (x(i) .lt. l(i)) then
                            prjctd = .true.
                            x(i) = l(i)
                        end if
                        nbdd = nbdd + 1
                    else if (nbd(i) .ge. 2 .and. x(i) .ge. u(i)) then
                        if (x(i) .gt. u(i)) then
                            prjctd = .true.
                            x(i) = u(i)
                        end if
                        nbdd = nbdd + 1
                    end if
                end if
            end do
            ! Initialize iwhere and assign values to cnstnd and boxed.
            do i = 1, n
                if (nbd(i) .ne. 2) boxed = .false.
                if (nbd(i) .eq. 0) then
                    ! this variable is always free
                    iwhere(i) = -1
                    ! otherwise set x(i)=mid(x(i), u(i), l(i)).
                else
                    cnstnd = .true.
                    if (nbd(i) .eq. 2 .and. u(i) - l(i) .le. zero) then
                        ! this variable is always fixed
                        iwhere(i) = 3
                    else
                        iwhere(i) = 0
                    end if
                end if
            end do
        end subroutine active

        subroutine bmv(m, sy, wt, col, v, p, info)
            integer  :: m, col, info
            real(dp) :: sy(m, m), wt(m, m), v(2*col), p(2*col)
            !     ************
            !
            !     Subroutine bmv
            !
            !     This subroutine computes the product of the 2m x 2m middle matrix
            !       in the compact L-BFGS formula of B and a 2m vector v;
            !       it returns the product in p.
            !
            !     m is an integer variable.
            !       On entry m is the maximum number of variable metric corrections
            !         used to define the limited memory matrix.
            !       On exit m is unchanged.
            !
            !     sy is a double precision array of dimension m x m.
            !       On entry sy specifies the matrix S'Y.
            !       On exit sy is unchanged.
            !
            !     wt is a double precision array of dimension m x m.
            !       On entry wt specifies the upper triangular matrix J' which is
            !         the Cholesky factor of (thetaS'S+LD^(-1)L').
            !       On exit wt is unchanged.
            !
            !     col is an integer variable.
            !       On entry col specifies the number of s-vectors (or y-vectors)
            !         stored in the compact L-BFGS formula.
            !       On exit col is unchanged.
            !
            !     v is a double precision array of dimension 2col.
            !       On entry v specifies vector v.
            !       On exit v is unchanged.
            !
            !     p is a double precision array of dimension 2col.
            !       On entry p is unspecified.
            !       On exit p is the product Mv.
            !
            !     info is an integer variable.
            !       On entry info is unspecified.
            !       On exit info = 0       for normal return,
            !                    = nonzero for abnormal return when the system
            !                                to be solved by dtrsl is singular.
            !
            !     Subprograms called:
            !
            !       Linpack ... dtrsl.
            !
            !
            !                           *  *  *
            !
            !     NEOS, November 1994. (Latest revision June 1996.)
            !     Optimization Technology Center.
            !     Argonne National Laboratory and Northwestern University.
            !     Written by
            !                        Ciyou Zhu
            !     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
            !
            !
            !     ************
            integer  :: i,k,i2
            real(dp) :: sum
            if (col .eq. 0) return
            !     PART I: solve [  D^(1/2)      O ] [ p1 ] = [ v1 ]
            !                   [ -L*D^(-1/2)   J ] [ p2 ]   [ v2 ].
            !       solve Jp2=v2+LD^(-1)v1.
            p(col + 1) = v(col + 1)
            do i = 2, col
                i2 = col + i
                sum = 0.0d0
                do k = 1, i - 1
                    sum = sum + sy(i,k)*v(k)/sy(k,k)
                end do
                p(i2) = v(i2) + sum
            end do
            !     Solve the triangular system
            call dtrsl(wt,m,col,p(col+1),11,info)
            if (info .ne. 0) return
            !       solve D^(1/2)p1=v1.
            do i = 1, col
                p(i) = v(i)/sqrt(sy(i,i))
            end do
            !     PART II: solve [ -D^(1/2)   D^(-1/2)*L'  ] [ p1 ] = [ p1 ]
            !                    [  0         J'           ] [ p2 ]   [ p2 ].
            !       solve J^Tp2=p2.
            call dtrsl(wt,m,col,p(col+1),01,info)
            if (info .ne. 0) return
            !       compute p1=-D^(-1/2)(p1-D^(-1/2)L'p2)
            !                 =-D^(-1/2)p1+D^(-1)L'p2.
            do i = 1, col
                p(i) = -p(i)/sqrt(sy(i,i))
            end do
            do i = 1, col
                sum = 0.d0
                do k = i + 1, col
                    sum = sum + sy(k,i)*p(col+k)/sy(i,i)
                end do
                p(i) = p(i) + sum
            end do
        end subroutine bmv

        subroutine cauchy(n, x, l, u, nbd, g, iorder, iwhere, t, d, xcp, &
            &                  m, wy, ws, sy, wt, theta, col, head, p, c, wbp, &
            &                  v, nseg, sbgnrm, info, epsmch)
            implicit none
            integer  :: n, m, head, col, nseg, info
            integer  :: nbd(n), iorder(n), iwhere(n)
            real(dp) :: theta, epsmch
            real(dp) :: x(n), l(n), u(n), g(n), t(n), d(n), xcp(n)
            real(dp) :: wy(n, col), ws(n, col), sy(m, m)
            real(dp) :: wt(m, m), p(2*m), c(2*m), wbp(2*m), v(2*m)
            !     ************
            !
            !     Subroutine cauchy
            !
            !     For given x, l, u, g (with sbgnrm > 0), and a limited memory
            !       BFGS matrix B defined in terms of matrices WY, WS, WT, and
            !       scalars head, col, and theta, this subroutine computes the
            !       generalized Cauchy point (GCP), defined as the first local
            !       minimizer of the quadratic
            !
            !                  Q(x + s) = g's + 1/2 s'Bs
            !
            !       along the projected gradient direction P(x-tg,l,u).
            !       The routine returns the GCP in xcp.
            !
            !     n is an integer variable.
            !       On entry n is the dimension of the problem.
            !       On exit n is unchanged.
            !
            !     x is a double precision array of dimension n.
            !       On entry x is the starting point for the GCP computation.
            !       On exit x is unchanged.
            !
            !     l is a double precision array of dimension n.
            !       On entry l is the lower bound of x.
            !       On exit l is unchanged.
            !
            !     u is a double precision array of dimension n.
            !       On entry u is the upper bound of x.
            !       On exit u is unchanged.
            !
            !     nbd is an integer array of dimension n.
            !       On entry nbd represents the type of bounds imposed on the
            !         variables, and must be specified as follows:
            !         nbd(i)=0 if x(i) is unbounded,
            !                1 if x(i) has only a lower bound,
            !                2 if x(i) has both lower and upper bounds, and
            !                3 if x(i) has only an upper bound.
            !       On exit nbd is unchanged.
            !
            !     g is a double precision array of dimension n.
            !       On entry g is the gradient of f(x).  g must be a nonzero vector.
            !       On exit g is unchanged.
            !
            !     iorder is an integer working array of dimension n.
            !       iorder will be used to store the breakpoints in the piecewise
            !       linear path and free variables encountered. On exit,
            !         iorder(1),...,iorder(nleft) are indices of breakpoints
            !                                which have not been encountered;
            !         iorder(nleft+1),...,iorder(nbreak) are indices of
            !                                     encountered breakpoints; and
            !         iorder(nfree),...,iorder(n) are indices of variables which
            !                 have no bound constraits along the search direction.
            !
            !     iwhere is an integer array of dimension n.
            !       On entry iwhere indicates only the permanently fixed (iwhere=3)
            !       or free (iwhere= -1) components of x.
            !       On exit iwhere records the status of the current x variables.
            !       iwhere(i)=-3  if x(i) is free and has bounds, but is not moved
            !                 0   if x(i) is free and has bounds, and is moved
            !                 1   if x(i) is fixed at l(i), and l(i) .ne. u(i)
            !                 2   if x(i) is fixed at u(i), and u(i) .ne. l(i)
            !                 3   if x(i) is always fixed, i.e.,  u(i)=x(i)=l(i)
            !                 -1  if x(i) is always free, i.e., it has no bounds.
            !
            !     t is a double precision working array of dimension n.
            !       t will be used to store the break points.
            !
            !     d is a double precision array of dimension n used to store
            !       the Cauchy direction P(x-tg)-x.
            !
            !     xcp is a double precision array of dimension n used to return the
            !       GCP on exit.
            !
            !     m is an integer variable.
            !       On entry m is the maximum number of variable metric corrections
            !         used to define the limited memory matrix.
            !       On exit m is unchanged.
            !
            !     ws, wy, sy, and wt are double precision arrays.
            !       On entry they store information that defines the
            !                             limited memory BFGS matrix:
            !         ws(n,m) stores S, a set of s-vectors;
            !         wy(n,m) stores Y, a set of y-vectors;
            !         sy(m,m) stores S'Y;
            !         wt(m,m) stores the
            !                 Cholesky factorization of (theta*S'S+LD^(-1)L').
            !       On exit these arrays are unchanged.
            !
            !     theta is a double precision variable.
            !       On entry theta is the scaling factor specifying B_0 = theta I.
            !       On exit theta is unchanged.
            !
            !     col is an integer variable.
            !       On entry col is the actual number of variable metric
            !         corrections stored so far.
            !       On exit col is unchanged.
            !
            !     head is an integer variable.
            !       On entry head is the location of the first s-vector (or y-vector)
            !         in S (or Y).
            !       On exit col is unchanged.
            !
            !     p is a double precision working array of dimension 2m.
            !       p will be used to store the vector p = W^(T)d.
            !
            !     c is a double precision working array of dimension 2m.
            !       c will be used to store the vector c = W^(T)(xcp-x).
            !
            !     wbp is a double precision working array of dimension 2m.
            !       wbp will be used to store the row of W corresponding
            !         to a breakpoint.
            !
            !     v is a double precision working array of dimension 2m.
            !
            !     nseg is an integer variable.
            !       On exit nseg records the number of quadratic segments explored
            !         in searching for the GCP.
            !
            !     sg and yg are double precision arrays of dimension m.
            !       On entry sg  and yg store S'g and Y'g correspondingly.
            !       On exit they are unchanged.
            !
            !     sbgnrm is a double precision variable.
            !       On entry sbgnrm is the norm of the projected gradient at x.
            !       On exit sbgnrm is unchanged.
            !
            !     info is an integer variable.
            !       On entry info is 0.
            !       On exit info = 0       for normal return,
            !                    = nonzero for abnormal return when the the system
            !                              used in routine bmv is singular.
            !
            !     Subprograms called:
            !
            !       L-BFGS-B Library ... hpsolb, bmv.
            !
            !       Linpack ... dscal dcopy, daxpy.
            !
            !
            !     References:
            !
            !       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
            !       memory algorithm for bound constrained optimization'',
            !       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
            !
            !       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: FORTRAN
            !       Subroutines for Large Scale Bound Constrained Optimization''
            !       Tech. Report, NAM-11, EECS Department, Northwestern University,
            !       1994.
            !
            !       (Postscript files of these papers are available via anonymous
            !        ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
            !
            !                           *  *  *
            !
            !     NEOS, November 1994. (Latest revision June 1996.)
            !     Optimization Technology Center.
            !     Argonne National Laboratory and Northwestern University.
            !     Written by
            !                        Ciyou Zhu
            !     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
            !
            !
            !     ************
            logical  :: xlower,xupper,bnded
            integer  :: i,j,col2,nfree,nbreak,pointr
            integer  :: ibp,nleft,ibkmin,iter
            real(dp) :: f1,f2,dt,dtm,tsum,dibp,zibp,dibp2,bkmin
            real(dp) :: tu,tl,wmc,wmp,wmw,tj,tj0,neggi,sbgnrm
            real(dp) :: f2_org
            !     Check the status of the variables, reset iwhere(i) if necessary;
            !       compute the Cauchy direction d and the breakpoints t; initialize
            !       the derivative f1 and the vector p = W'd (for theta = 1).
            if (sbgnrm .le. zero) then
                call dcopy(n,x,1,xcp,1)
                return
            end if
            bnded = .true.
            nfree = n + 1
            nbreak = 0
            ibkmin = 0
            bkmin = zero
            col2 = 2*col
            f1 = zero
            !     We set p to zero and build it up as we determine d.
            do i = 1, col2
                p(i) = zero
            end do
            !     In the following loop we determine for each variable its bound
            !        status and its breakpoint, and update p accordingly.
            !        Smallest breakpoint is identified.
            do i = 1, n
                neggi = -g(i)
                if (iwhere(i) .ne. 3 .and. iwhere(i) .ne. -1) then
                    !             if x(i) is not a constant and has bounds,
                    !             compute the difference between x(i) and its bounds.
                    if (nbd(i) .le. 2) tl = x(i) - l(i)
                    if (nbd(i) .ge. 2) tu = u(i) - x(i)
                    !           If a variable is close enough to a bound
                    !             we treat it as at bound.
                    xlower = nbd(i) .le. 2 .and. tl .le. zero
                    xupper = nbd(i) .ge. 2 .and. tu .le. zero
                    !              reset iwhere(i).
                    iwhere(i) = 0
                    if (xlower) then
                        if (neggi .le. zero) iwhere(i) = 1
                    else if (xupper) then
                        if (neggi .ge. zero) iwhere(i) = 2
                    else
                        if (abs(neggi) .le. zero) iwhere(i) = -3
                    end if
                end if
                pointr = head
                if (iwhere(i) .ne. 0 .and. iwhere(i) .ne. -1) then
                    d(i) = zero
                else
                    d(i) = neggi
                    f1 = f1 - neggi*neggi
                    !             calculate p := p - W'e_i* (g_i).
                    do j = 1, col
                        p(j) = p(j) +  wy(i,pointr)* neggi
                        p(col + j) = p(col + j) + ws(i,pointr)*neggi
                        pointr = mod(pointr,m) + 1
                    end do
                    if (nbd(i) .le. 2 .and. nbd(i) .ne. 0 .and. neggi .lt. zero) then
                        ! x(i) + d(i) is bounded; compute t(i).
                        nbreak = nbreak + 1
                        iorder(nbreak) = i
                        t(nbreak) = tl/(-neggi)
                        if (nbreak .eq. 1 .or. t(nbreak) .lt. bkmin) then
                            bkmin = t(nbreak)
                            ibkmin = nbreak
                        end if
                    else if (nbd(i) .ge. 2 .and. neggi .gt. zero) then
                        ! x(i) + d(i) is bounded; compute t(i).
                        nbreak = nbreak + 1
                        iorder(nbreak) = i
                        t(nbreak) = tu/neggi
                        if (nbreak .eq. 1 .or. t(nbreak) .lt. bkmin) then
                            bkmin = t(nbreak)
                            ibkmin = nbreak
                        end if
                    else
                        ! x(i) + d(i) is not bounded.
                        nfree = nfree - 1
                        iorder(nfree) = i
                        if (abs(neggi) .gt. zero) bnded = .false.
                    end if
                end if
            end do
            !     The indices of the nonzero components of d are now stored
            !       in iorder(1),...,iorder(nbreak) and iorder(nfree),...,iorder(n).
            !       The smallest of the nbreak breakpoints is in t(ibkmin)=bkmin.
            if (.not. is_equal(theta , one) ) then  !theta .ne. one
                !                   complete the initialization of p for theta not= one.
                call dscal(col,theta,p(col+1),1)
            end if
            !     Initialize GCP xcp = x.
            call dcopy(n,x,1,xcp,1)
            if (nbreak .eq. 0 .and. nfree .eq. n + 1) then
                ! is a zero vector, return with the initial xcp as GCP.
                return
            end if
            !     Initialize c = W'(xcp - x) = 0.
            do j = 1, col2
                c(j) = zero
            end do
            !     Initialize derivative f2.
            f2 =  -theta*f1
            f2_org  =  f2
            if (col .gt. 0) then
                call bmv(m,sy,wt,col,p,v,info)
                if (info .ne. 0) return
                f2 = f2 - ddot_lbfgsb(col2,v,1,p,1)
            end if
            dtm = -f1/f2
            tsum = zero
            nseg = 1
            !     If there are no breakpoints, locate the GCP and return.
            if (nbreak .eq. 0) goto 888
            nleft = nbreak
            iter = 1
            tj = zero
            !------------------- the beginning of the loop -------------------------
777         continue
            !     Find the next smallest breakpoint;
            !       compute dt = t(nleft) - t(nleft + 1).
            tj0 = tj
            if (iter .eq. 1) then
                !         Since we already have the smallest breakpoint we need not do
                !         heapsort yet. Often only one breakpoint is used and the
                !         cost of heapsort is avoided.
                tj = bkmin
                ibp = iorder(ibkmin)
            else
                if (iter .eq. 2) then
                    !             Replace the already used smallest breakpoint with the
                    !             breakpoint numbered nbreak > nlast, before heapsort call.
                    if (ibkmin .ne. nbreak) then
                        t(ibkmin) = t(nbreak)
                        iorder(ibkmin) = iorder(nbreak)
                    end if
                    !        Update heap structure of breakpoints
                    !           (if iter=2, initialize heap).
                end if
                call hpsolb(nleft,t,iorder,iter-2)
                tj = t(nleft)
                ibp = iorder(nleft)
            end if
            dt = tj - tj0
            !     If a minimizer is within this interval, locate the GCP and return.
            if (dtm .lt. dt) goto 888
            !     Otherwise fix one variable and
            !       reset the corresponding component of d to zero.
            tsum = tsum + dt
            nleft = nleft - 1
            iter = iter + 1
            dibp = d(ibp)
            d(ibp) = zero
            if (dibp .gt. zero) then
                zibp = u(ibp) - x(ibp)
                xcp(ibp) = u(ibp)
                iwhere(ibp) = 2
            else
                zibp = l(ibp) - x(ibp)
                xcp(ibp) = l(ibp)
                iwhere(ibp) = 1
            end if
            if (nleft .eq. 0 .and. nbreak .eq. n) then
                ! all n variables are fixed,
                !  return with xcp as GCP.
                dtm = dt
                goto 999
            end if
            !     Update the derivative information.
            nseg = nseg + 1
            dibp2 = dibp**2
            !     Update f1 and f2.
            !        temporarily set f1 and f2 for col=0.
            f1 = f1 + dt*f2 + dibp2 - theta*dibp*zibp
            f2 = f2 - theta*dibp2
            if (col .gt. 0) then
                !                          update c = c + dt*p.
                call daxpy(col2,dt,p,1,c,1)
                !           choose wbp,
                !           the row of W corresponding to the breakpoint encountered.
                pointr = head
                do j = 1,col
                    wbp(j) = wy(ibp,pointr)
                    wbp(col + j) = theta*ws(ibp,pointr)
                    pointr = mod(pointr,m) + 1
                end do
                ! compute (wbp)Mc, (wbp)Mp, and (wbp)M(wbp)'.
                call bmv(m,sy,wt,col,wbp,v,info)
                if (info .ne. 0) return
                wmc = ddot_lbfgsb(col2,c,1,v,1)
                wmp = ddot_lbfgsb(col2,p,1,v,1)
                wmw = ddot_lbfgsb(col2,wbp,1,v,1)
                ! update p = p - dibp*wbp.
                call daxpy(col2,-dibp,wbp,1,p,1)
                ! complete updating f1 and f2 while col > 0.
                f1 = f1 + dibp*wmc
                f2 = f2 + 2.0d0*dibp*wmp - dibp2*wmw
            end if
            f2 = max(epsmch*f2_org,f2)
            if (nleft .gt. 0) then
                dtm = -f1/f2
                goto 777
                ! to repeat the loop for unsearched intervals.
            else if(bnded) then
                f1 = zero
                f2 = zero
                dtm = zero
            else
                dtm = -f1/f2
            end if
            !------------------- the end of the loop -------------------------------
888         continue
            if (dtm .le. zero) dtm = zero
            tsum = tsum + dtm
            !     Move free variables (i.e., the ones w/o breakpoints) and
            !       the variables whose breakpoints haven't been reached.
            call daxpy(n,tsum,d,1,xcp,1)
999         continue
            !     Update c = c + dtm*p = W'(x^c - x)
            !       which will be used in computing r = Z'(B(x^c - x) + g).
            if (col .gt. 0) call daxpy(col2,dtm,p,1,c,1)
        end subroutine cauchy

        subroutine cmprlb(n, m, x, g, ws, wy, sy, wt, z, r, wa, index, &
            &                 theta, col, head, nfree, cnstnd, info)
            logical  :: cnstnd
            integer  :: n, m, col, head, nfree, info, index(n)
            real(dp) :: theta
            real(dp) :: x(n), g(n), z(n), r(n), wa(4*m)
            real(dp) :: ws(n, m), wy(n, m), sy(m, m), wt(m, m)
            !     ************
            !
            !     Subroutine cmprlb
            !
            !       This subroutine computes r=-Z'B(xcp-xk)-Z'g by using
            !         wa(2m+1)=W'(xcp-x) from subroutine cauchy.
            !
            !     Subprograms called:
            !
            !       L-BFGS-B Library ... bmv.
            !
            !
            !                           *  *  *
            !
            !     NEOS, November 1994. (Latest revision June 1996.)
            !     Optimization Technology Center.
            !     Argonne National Laboratory and Northwestern University.
            !     Written by
            !                        Ciyou Zhu
            !     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
            !
            !
            !     ************
            integer  :: i,j,k,pointr
            real(dp) :: a1,a2
            if (.not. cnstnd .and. col .gt. 0) then
                do i = 1, n
                    r(i) = -g(i)
                end do
            else
                do i = 1, nfree
                    k = index(i)
                    r(i) = -theta*(z(k) - x(k)) - g(k)
                end do
                call bmv(m,sy,wt,col,wa(2*m+1),wa(1),info)
                if (info .ne. 0) then
                    info = -8
                    return
                end if
                pointr = head
                do j = 1, col
                    a1 = wa(j)
                    a2 = theta*wa(col + j)
                    do i = 1, nfree
                        k = index(i)
                        r(i) = r(i) + wy(k,pointr)*a1 + ws(k,pointr)*a2
                    end do
                    pointr = mod(pointr,m) + 1
                end do
            end if
        end subroutine cmprlb

        subroutine errclb(n, m, factr, l, u, nbd, task, info, k)
            character(len=60) :: task
            integer   :: n, m, info, k, nbd(n)
            real(dp)  :: factr, l(n), u(n)
            !     ************
            !
            !     Subroutine errclb
            !
            !     This subroutine checks the validity of the input data.
            !
            !
            !                           *  *  *
            !
            !     NEOS, November 1994. (Latest revision June 1996.)
            !     Optimization Technology Center.
            !     Argonne National Laboratory and Northwestern University.
            !     Written by
            !                        Ciyou Zhu
            !     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
            !
            !
            !     ************
            integer :: i
            !     Check the input arguments for errors.
            if (n .le. 0) task = 'ERROR: N .LE. 0'
            if (m .le. 0) task = 'ERROR: M .LE. 0'
            if (factr .lt. zero) task = 'ERROR: FACTR .LT. 0'
            !     Check the validity of the arrays nbd(i), u(i), and l(i).
            do i = 1, n
                if (nbd(i) .lt. 0 .or. nbd(i) .gt. 3) then
                    ! return
                    task = 'ERROR: INVALID NBD'
                    info = -6
                    k = i
                end if
                if (nbd(i) .eq. 2) then
                    if (l(i) .gt. u(i)) then
                        ! return
                        task = 'ERROR: NO FEASIBLE SOLUTION'
                        info = -7
                        k = i
                    end if
                end if
            end do
        end subroutine errclb

        subroutine formk(n, nsub, ind, nenter, ileave, indx2, iupdat, &
            &            updatd, wn, wn1, m, ws, wy, sy, theta, col, &
            &            head, info)
            integer  :: n, nsub, m, col, head, nenter, ileave, iupdat
            integer  :: info, ind(n), indx2(n)
            real(dp) :: theta, wn(2*m, 2*m), wn1(2*m, 2*m)
            real(dp) :: ws(n, m), wy(n, m), sy(m, m)
            logical  :: updatd
            !     ************
            !
            !     Subroutine formk
            !
            !     This subroutine forms  the LEL^T factorization of the indefinite
            !
            !       matrix    K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
            !                     [L_a -R_z           theta*S'AA'S ]
            !                                                    where E = [-I  0]
            !                                                              [ 0  I]
            !     The matrix K can be shown to be equal to the matrix M^[-1]N
            !       occurring in section 5.1 of [1], as well as to the matrix
            !       Mbar^[-1] Nbar in section 5.3.
            !
            !     n is an integer variable.
            !       On entry n is the dimension of the problem.
            !       On exit n is unchanged.
            !
            !     nsub is an integer variable
            !       On entry nsub is the number of subspace variables in free set.
            !       On exit nsub is not changed.
            !
            !     ind is an integer array of dimension nsub.
            !       On entry ind specifies the indices of subspace variables.
            !       On exit ind is unchanged.
            !
            !     nenter is an integer variable.
            !       On entry nenter is the number of variables entering the
            !         free set.
            !       On exit nenter is unchanged.
            !
            !     ileave is an integer variable.
            !       On entry indx2(ileave),...,indx2(n) are the variables leaving
            !         the free set.
            !       On exit ileave is unchanged.
            !
            !     indx2 is an integer array of dimension n.
            !       On entry indx2(1),...,indx2(nenter) are the variables entering
            !         the free set, while indx2(ileave),...,indx2(n) are the
            !         variables leaving the free set.
            !       On exit indx2 is unchanged.
            !
            !     iupdat is an integer variable.
            !       On entry iupdat is the total number of BFGS updates made so far.
            !       On exit iupdat is unchanged.
            !
            !     updatd is a logical variable.
            !       On entry 'updatd' is true if the L-BFGS matrix is updatd.
            !       On exit 'updatd' is unchanged.
            !
            !     wn is a double precision array of dimension 2m x 2m.
            !       On entry wn is unspecified.
            !       On exit the upper triangle of wn stores the LEL^T factorization
            !         of the 2*col x 2*col indefinite matrix
            !                     [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
            !                     [L_a -R_z           theta*S'AA'S ]
            !
            !     wn1 is a double precision array of dimension 2m x 2m.
            !       On entry wn1 stores the lower triangular part of
            !                     [Y' ZZ'Y   L_a'+R_z']
            !                     [L_a+R_z   S'AA'S   ]
            !         in the previous iteration.
            !       On exit wn1 stores the corresponding updated matrices.
            !       The purpose of wn1 is just to store these inner products
            !       so they can be easily updated and inserted into wn.
            !
            !     m is an integer variable.
            !       On entry m is the maximum number of variable metric corrections
            !         used to define the limited memory matrix.
            !       On exit m is unchanged.
            !
            !     ws, wy, sy, and wtyy are double precision arrays;
            !     theta is a double precision variable;
            !     col is an integer variable;
            !     head is an integer variable.
            !       On entry they store the information defining the
            !                                          limited memory BFGS matrix:
            !         ws(n,m) stores S, a set of s-vectors;
            !         wy(n,m) stores Y, a set of y-vectors;
            !         sy(m,m) stores S'Y;
            !         wtyy(m,m) stores the Cholesky factorization
            !                                   of (theta*S'S+LD^(-1)L')
            !         theta is the scaling factor specifying B_0 = theta I;
            !         col is the number of variable metric corrections stored;
            !         head is the location of the 1st s- (or y-) vector in S (or Y).
            !       On exit they are unchanged.
            !
            !     info is an integer variable.
            !       On entry info is unspecified.
            !       On exit info =  0 for normal return;
            !                    = -1 when the 1st Cholesky factorization failed;
            !                    = -2 when the 2st Cholesky factorization failed.
            !
            !     Subprograms called:
            !
            !       Linpack ... dcopy, dpofa, dtrsl.
            !
            !
            !     References:
            !       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
            !       memory algorithm for bound constrained optimization'',
            !       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
            !
            !       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: a
            !       limited memory FORTRAN code for solving bound constrained
            !       optimization problems'', Tech. Report, NAM-11, EECS Department,
            !       Northwestern University, 1994.
            !
            !       (Postscript files of these papers are available via anonymous
            !        ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
            !
            !                           *  *  *
            !
            !     NEOS, November 1994. (Latest revision June 1996.)
            !     Optimization Technology Center.
            !     Argonne National Laboratory and Northwestern University.
            !     Written by
            !                        Ciyou Zhu
            !     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
            !
            !
            !     ************
            integer  :: m2,ipntr,jpntr,iy,is,jy,js,is1,js1,k1,i,k
            integer  :: col2,pbegin,pend,dbegin,dend,upcl
            real(dp) :: temp1,temp2,temp3,temp4
            !     Form the lower triangular part of
            !               WN1 = [Y' ZZ'Y   L_a'+R_z']
            !                     [L_a+R_z   S'AA'S   ]
            !        where L_a is the strictly lower triangular part of S'AA'Y
            !              R_z is the upper triangular part of S'ZZ'Y.
            if (updatd) then
                if (iupdat .gt. m) then
                    ! shift old part of WN1.
                    do jy = 1, m - 1
                        js = m + jy
                        call dcopy(m-jy,wn1(jy+1,jy+1),1,wn1(jy,jy),1)
                        call dcopy(m-jy,wn1(js+1,js+1),1,wn1(js,js),1)
                        call dcopy(m-1,wn1(m+2,jy+1),1,wn1(m+1,jy),1)
                    end do
                end if
                ! put new rows in blocks (1,1), (2,1) and (2,2).
                pbegin = 1
                pend = nsub
                dbegin = nsub + 1
                dend = n
                iy = col
                is = m + col
                ipntr = head + col - 1
                if (ipntr .gt. m) ipntr = ipntr - m
                jpntr = head
                do jy = 1, col
                    js = m + jy
                    temp1 = zero
                    temp2 = zero
                    temp3 = zero
                    ! compute element jy of row 'col' of Y'ZZ'Y
                    do k = pbegin, pend
                        k1 = ind(k)
                        temp1 = temp1 + wy(k1,ipntr)*wy(k1,jpntr)
                    end do
                    ! compute elements jy of row 'col' of L_a and S'AA'S
                    do k = dbegin, dend
                        k1 = ind(k)
                        temp2 = temp2 + ws(k1,ipntr)*ws(k1,jpntr)
                        temp3 = temp3 + ws(k1,ipntr)*wy(k1,jpntr)
                    end do
                    wn1(iy,jy) = temp1
                    wn1(is,js) = temp2
                    wn1(is,jy) = temp3
                    jpntr = mod(jpntr,m) + 1
                end do
                ! put new column in block (2,1).
                jy = col
                jpntr = head + col - 1
                if (jpntr .gt. m) jpntr = jpntr - m
                ipntr = head
                do i = 1, col
                    is = m + i
                    temp3 = zero
                    ! compute element i of column 'col' of R_z
                    do k = pbegin, pend
                        k1 = ind(k)
                        temp3 = temp3 + ws(k1,ipntr)*wy(k1,jpntr)
                    end do
                    ipntr = mod(ipntr,m) + 1
                    wn1(is,jy) = temp3
                end do
                upcl = col - 1
            else
                upcl = col
            end if
            ! modify the old parts in blocks (1,1) and (2,2) due to changes
            ! in the set of free variables.
            ipntr = head
            do iy = 1, upcl
                is = m + iy
                jpntr = head
                do jy = 1, iy
                    js = m + jy
                    temp1 = zero
                    temp2 = zero
                    temp3 = zero
                    temp4 = zero
                    do k = 1, nenter
                        k1 = indx2(k)
                        temp1 = temp1 + wy(k1,ipntr)*wy(k1,jpntr)
                        temp2 = temp2 + ws(k1,ipntr)*ws(k1,jpntr)
                    end do
                    do k = ileave, n
                        k1 = indx2(k)
                        temp3 = temp3 + wy(k1,ipntr)*wy(k1,jpntr)
                        temp4 = temp4 + ws(k1,ipntr)*ws(k1,jpntr)
                    end do
                    wn1(iy,jy) = wn1(iy,jy) + temp1 - temp3
                    wn1(is,js) = wn1(is,js) - temp2 + temp4
                    jpntr = mod(jpntr,m) + 1
                end do
                ipntr = mod(ipntr,m) + 1
            end do
            ! modify the old parts in block (2,1).
            ipntr = head
            do is = m + 1, m + upcl
                jpntr = head
                do jy = 1, upcl
                    temp1 = zero
                    temp3 = zero
                    do k = 1, nenter
                        k1 = indx2(k)
                        temp1 = temp1 + ws(k1,ipntr)*wy(k1,jpntr)
                    end do
                    do k = ileave, n
                        k1 = indx2(k)
                        temp3 = temp3 + ws(k1,ipntr)*wy(k1,jpntr)
                    end do
                    if (is .le. jy + m) then
                        wn1(is,jy) = wn1(is,jy) + temp1 - temp3
                    else
                        wn1(is,jy) = wn1(is,jy) - temp1 + temp3
                    end if
                    jpntr = mod(jpntr,m) + 1
                end do
                ipntr = mod(ipntr,m) + 1
            end do
            ! Form the upper triangle of WN = [D+Y' ZZ'Y/theta   -L_a'+R_z' ]
            !                                 [-L_a +R_z        S'AA'S*theta]
            m2 = 2*m
            do iy = 1, col
                is = col + iy
                is1 = m + iy
                do jy = 1, iy
                    js = col + jy
                    js1 = m + jy
                    wn(jy,iy) = wn1(iy,jy)/theta
                    wn(js,is) = wn1(is1,js1)*theta
                end do
                do jy = 1, iy - 1
                    wn(jy,is) = -wn1(is1,jy)
                end do
                do jy = iy, col
                    wn(jy,is) = wn1(is1,jy)
                end do
                wn(iy,iy) = wn(iy,iy) + sy(iy,iy)
            end do
            ! Form the upper triangle of WN= [  LL'            L^-1(-L_a'+R_z')]
            !                                [(-L_a +R_z)L'^-1   S'AA'S*theta  ]
            ! first Cholesky factor (1,1) block of wn to get LL'
            !                          with L' stored in the upper triangle of wn.
            call dpofa(wn,m2,col,info)
            if (info .ne. 0) then
                info = -1
                return
            end if
            ! then form L^-1(-L_a'+R_z') in the (1,2) block.
            col2 = 2*col
            do js = col+1 ,col2
                call dtrsl(wn,m2,col,wn(1,js),11,info)
            end do
            ! Form S'AA'S*theta + (L^-1(-L_a'+R_z'))'L^-1(-L_a'+R_z') in the
            ! upper triangle of (2,2) block of wn.
            do is = col+1, col2
                do js = is, col2
                    wn(is,js) = wn(is,js) + ddot_lbfgsb(col,wn(1,is),1,wn(1,js),1)
                end do
            end do
            ! Cholesky factorization of (2,2) block of wn.
            call dpofa(wn(col+1,col+1),m2,col,info)
            if (info .ne. 0) then
                info = -2
                return
            end if
        end subroutine formk

        subroutine formt(m, wt, sy, ss, col, theta, info)
            integer  :: m, col, info
            real(dp) :: theta, wt(m, m), sy(m, m), ss(m, m)
            !     ************
            !
            !     Subroutine formt
            !
            !       This subroutine forms the upper half of the pos. def. and symm.
            !         T = theta*SS + L*D^(-1)*L', stores T in the upper triangle
            !         of the array wt, and performs the Cholesky factorization of T
            !         to produce J*J', with J' stored in the upper triangle of wt.
            !
            !     Subprograms called:
            !
            !       Linpack ... dpofa.
            !
            !
            !                           *  *  *
            !
            !     NEOS, November 1994. (Latest revision June 1996.)
            !     Optimization Technology Center.
            !     Argonne National Laboratory and Northwestern University.
            !     Written by
            !                        Ciyou Zhu
            !     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
            !
            !
            !     ************
            integer  :: i,j,k,k1
            real(dp) :: ddum
            ! Form the upper half of  T = theta*SS + L*D^(-1)*L',
            ! store T in the upper triangle of the array wt.
            do j = 1, col
                wt(1,j) = theta*ss(1,j)
            end do
            do i = 2, col
                do j = i, col
                    k1 = min(i,j) - 1
                    ddum  = zero
                    do k = 1, k1
                        ddum  = ddum + sy(i,k)*sy(j,k)/sy(k,k)
                    end do
                    wt(i,j) = ddum + theta*ss(i,j)
                end do
            end do
            ! Cholesky factorize T to J*J' with
            ! J' stored in the upper triangle of wt.
            call dpofa(wt,m,col,info)
            if (info .ne. 0) then
                info = -3
            end if
        end subroutine formt

        subroutine freev(n, nfree, index, nenter, ileave, indx2, &
            &                 iwhere, wrk, updatd, cnstnd, iter)
            integer :: n, nfree, nenter, ileave, iter
            integer :: index(n), indx2(n), iwhere(n)
            logical :: wrk, updatd, cnstnd
            !     ************
            !
            !     Subroutine freev
            !
            !     This subroutine counts the entering and leaving variables when
            !       iter > 0, and finds the index set of free and active variables
            !       at the GCP.
            !
            !     cnstnd is a logical variable indicating whether bounds are present
            !
            !     index is an integer array of dimension n
            !       for i=1,...,nfree, index(i) are the indices of free variables
            !       for i=nfree+1,...,n, index(i) are the indices of bound variables
            !       On entry after the first iteration, index gives
            !         the free variables at the previous iteration.
            !       On exit it gives the free variables based on the determination
            !         in cauchy using the array iwhere.
            !
            !     indx2 is an integer array of dimension n
            !       On entry indx2 is unspecified.
            !       On exit with iter>0, indx2 indicates which variables
            !          have changed status since the previous iteration.
            !       For i= 1,...,nenter, indx2(i) have changed from bound to free.
            !       For i= ileave+1,...,n, indx2(i) have changed from free to bound.
            !
            !
            !                           *  *  *
            !
            !     NEOS, November 1994. (Latest revision June 1996.)
            !     Optimization Technology Center.
            !     Argonne National Laboratory and Northwestern University.
            !     Written by
            !                        Ciyou Zhu
            !     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
            !
            !
            !     ************
            integer :: iact,i,k
            nenter = 0
            ileave = n + 1
            if (iter .gt. 0 .and. cnstnd) then
                ! count the entering and leaving variables.
                do i = 1, nfree
                    k = index(i)
                    ! write(6,*) ' k  = index(i) ', k
                    ! write(6,*) ' index = ', i
                    if (iwhere(k) .gt. 0) then
                        ileave = ileave - 1
                        indx2(ileave) = k
                    end if
                end do
                do i = 1 + nfree, n
                    k = index(i)
                    if (iwhere(k) .le. 0) then
                        nenter = nenter + 1
                        indx2(nenter) = k
                    end if
                end do
            end if
            wrk = (ileave .lt. n+1) .or. (nenter .gt. 0) .or. updatd
            !     Find the index set of free and active variables at the GCP.
            nfree = 0
            iact = n + 1
            do i = 1, n
                if (iwhere(i) .le. 0) then
                    nfree = nfree + 1
                    index(nfree) = i
                else
                    iact = iact - 1
                    index(iact) = i
                end if
            end do
        end subroutine freev

        subroutine hpsolb(n, t, iorder, iheap)
            integer  :: iheap, n, iorder(n)
            real(dp) :: t(n)
            !     ************
            !
            !     Subroutine hpsolb
            !
            !     This subroutine sorts out the least element of t, and puts the
            !       remaining elements of t in a heap.
            !
            !     n is an integer variable.
            !       On entry n is the dimension of the arrays t and iorder.
            !       On exit n is unchanged.
            !
            !     t is a double precision array of dimension n.
            !       On entry t stores the elements to be sorted,
            !       On exit t(n) stores the least elements of t, and t(1) to t(n-1)
            !         stores the remaining elements in the form of a heap.
            !
            !     iorder is an integer array of dimension n.
            !       On entry iorder(i) is the index of t(i).
            !       On exit iorder(i) is still the index of t(i), but iorder may be
            !         permuted in accordance with t.
            !
            !     iheap is an integer variable specifying the task.
            !       On entry iheap should be set as follows:
            !         iheap .eq. 0 if t(1) to t(n) is not in the form of a heap,
            !         iheap .ne. 0 if otherwise.
            !       On exit iheap is unchanged.
            !
            !
            !     References:
            !       Algorithm 232 of CACM (J. W. J. Williams): HEAPSORT.
            !
            !                           *  *  *
            !
            !     NEOS, November 1994. (Latest revision June 1996.)
            !     Optimization Technology Center.
            !     Argonne National Laboratory and Northwestern University.
            !     Written by
            !                        Ciyou Zhu
            !     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
            !
            !     ************
            integer  :: i,j,k,indxin,indxou
            real(dp) :: ddum,out
            if (iheap .eq. 0) then
                ! Rearrange the elements t(1) to t(n) to form a heap.
                do k = 2, n
                    ddum  = t(k)
                    indxin = iorder(k)
                    ! Add ddum to the heap.
                    i = k
10                  continue
                    if (i.gt.1) then
                        j = i/2
                        if (ddum .lt. t(j)) then
                            t(i) = t(j)
                            iorder(i) = iorder(j)
                            i = j
                            goto 10
                        end if
                    end if
                    t(i) = ddum
                    iorder(i) = indxin
                end do
            end if
            ! Assign to 'out' the value of t(1), the least member of the heap,
            ! and rearrange the remaining members to form a heap as
            ! elements 1 to n-1 of t.
            if (n .gt. 1) then
                i = 1
                out = t(1)
                indxou = iorder(1)
                ddum  = t(n)
                indxin  = iorder(n)
                !        Restore the heap
30              continue
                j = i+i
                if (j .le. n-1) then
                    if (t(j+1) .lt. t(j)) j = j+1
                    if (t(j) .lt. ddum ) then
                        t(i) = t(j)
                        iorder(i) = iorder(j)
                        i = j
                        goto 30
                    end if
                end if
                t(i) = ddum
                iorder(i) = indxin
                ! Put the least member in t(n).
                t(n) = out
                iorder(n) = indxou
            end if
        end subroutine hpsolb

        subroutine lnsrlb(n, l, u, nbd, x, f, fold, gd, gdold, g, d, r, t, &
            &                  z, stp, dnorm, dtd, xstep, stpmx, iter, ifun, &
            &                  iback, nfgv, info, task, boxed, cnstnd, csave, &
            &                  isave, dsave)
            character(len=60) :: task, csave
            logical   :: boxed, cnstnd
            integer   :: n, iter, ifun, iback, nfgv, info
            integer   :: nbd(n), isave(2)
            real(dp)  :: f, fold, gd, gdold, stp, dnorm, dtd, xstep
            real(dp)  :: stpmx, x(n), l(n), u(n), g(n), d(n), r(n), t(n)
            real(dp)  :: z(n), dsave(13)
            !     **********
            !
            !     Subroutine lnsrlb
            !
            !     This subroutine calls subroutine dcsrch from the Minpack2 library
            !       to perform the line search.  Subroutine dscrch is safeguarded so
            !       that all trial points lie within the feasible region.
            !
            !     Subprograms called:
            !
            !       Minpack2 Library ... dcsrch.
            !
            !       Linpack ... dtrsl, ddot.
            !
            !
            !                           *  *  *
            !
            !     NEOS, November 1994. (Latest revision June 1996.)
            !     Optimizatino Technology Center.
            !     Argonne National Laboratory and Northwestern University.
            !     Written by
            !                        Ciyou Zhu
            !     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
            !
            !
            !     **********
            integer  :: i
            real(dp) :: a1,a2
            if (task(1:5) .eq. 'FG_LN') goto 556
            dtd = ddot_lbfgsb(n,d,1,d,1)
            dnorm = sqrt(dtd)
            ! Determine the maximum step length.
            stpmx = big
            if (cnstnd) then
                if (iter .eq. 0) then
                    stpmx = one
                else
                    do i = 1, n
                        a1 = d(i)
                        if (nbd(i) .ne. 0) then
                            if (a1 .lt. zero .and. nbd(i) .le. 2) then
                                a2 = l(i) - x(i)
                                if (a2 .ge. zero) then
                                    stpmx = zero
                                else if (a1*stpmx .lt. a2) then
                                    stpmx = a2/a1
                                end if
                            else if (a1 .gt. zero .and. nbd(i) .ge. 2) then
                                a2 = u(i) - x(i)
                                if (a2 .le. zero) then
                                    stpmx = zero
                                else if (a1*stpmx .gt. a2) then
                                    stpmx = a2/a1
                                end if
                            end if
                        end if
                    end do
                end if
            end if
            if (iter .eq. 0 .and. .not. boxed) then
                stp = min(one/dnorm, stpmx)
            else
                stp = one
            end if
            call dcopy(n,x,1,t,1)
            call dcopy(n,g,1,r,1)
            fold = f
            ifun = 0
            iback = 0
            csave = 'START'
556         continue
            gd = ddot_lbfgsb(n,g,1,d,1)
            if (ifun .eq. 0) then
                gdold=gd
                if (gd .ge. zero) then
                    ! the directional derivative >=0.
                    ! Line search is impossible.
                    info = -4
                    return
                end if
            end if
            call dcsrch(f,gd,stp,real(spec%ftol,dp),real(spec%gtol,dp),xtol_lbfgsb,zero,stpmx,csave,isave,dsave)
            xstep = stp*dnorm
            if (csave(1:4) .ne. 'CONV' .and. csave(1:4) .ne. 'WARN') then
                task = 'FG_LNSRCH'
                ifun = ifun + 1
                nfgv = nfgv + 1
                iback = ifun - 1
                if (is_equal(stp , one)) then   !stp .eq. one
                    call dcopy(n,z,1,x,1)
                else
                    do i = 1, n
                        x(i) = stp*d(i) + t(i)
                    end do
                end if
            else
                task = 'NEW_X'
            end if
        end subroutine lnsrlb

        subroutine matupd(n, m, ws, wy, sy, ss, d, r, itail, &
            &                  iupdat, col, head, theta, rr, dr, stp, dtd)
            integer  :: n, m, itail, iupdat, col, head
            real(dp) :: theta, rr, dr, stp, dtd, d(n), r(n)
            real(dp) :: ws(n, m), wy(n, m), sy(m, m), ss(m, m)
            !     ************
            !
            !     Subroutine matupd
            !
            !       This subroutine updates matrices WS and WY, and forms the
            !         middle matrix in B.
            !
            !     Subprograms called:
            !
            !       Linpack ... dcopy, ddot.
            !
            !
            !                           *  *  *
            !
            !     NEOS, November 1994. (Latest revision June 1996.)
            !     Optimization Technology Center.
            !     Argonne National Laboratory and Northwestern University.
            !     Written by
            !                        Ciyou Zhu
            !     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
            !
            !
            !     ************
            integer  :: j,pointr
            ! Set pointers for matrices WS and WY.
            if (iupdat .le. m) then
                col = iupdat
                itail = mod(head+iupdat-2,m) + 1
            else
                itail = mod(itail,m) + 1
                head = mod(head,m) + 1
            end if
            ! Update matrices WS and WY.
            call dcopy(n,d,1,ws(1,itail),1)
            call dcopy(n,r,1,wy(1,itail),1)
            !     Set theta=yy/ys.
            theta = rr/dr
            ! Form the middle matrix in B.
            ! update the upper triangle of SS,
            ! and the lower triangle of SY:
            if (iupdat .gt. m) then
                ! move old information
                do j = 1, col - 1
                    call dcopy(j,ss(2,j+1),1,ss(1,j),1)
                    call dcopy(col-j,sy(j+1,j+1),1,sy(j,j),1)
                end do
            end if
            ! add new information: the last row of SY
            ! and the last column of SS:
            pointr = head
            do j = 1, col - 1
                sy(col,j) = ddot_lbfgsb(n,d,1,wy(1,pointr),1)
                ss(j,col) = ddot_lbfgsb(n,ws(1,pointr),1,d,1)
                pointr = mod(pointr,m) + 1
            end do
            if  (is_equal(stp , one)) then
                ss(col,col) = dtd
            else
                ss(col,col) = stp*stp*dtd
            end if
            sy(col,col) = dr
        end subroutine matupd

        subroutine projgr(n, l, u, nbd, x, g, sbgnrm)
            integer  :: n, nbd(n)
            real(dp) :: sbgnrm, x(n), l(n), u(n), g(n)
            !     ************
            !
            !     Subroutine projgr
            !
            !     This subroutine computes the infinity norm of the projected
            !       gradient.
            !
            !
            !                           *  *  *
            !
            !     NEOS, November 1994. (Latest revision June 1996.)
            !     Optimization Technology Center.
            !     Argonne National Laboratory and Northwestern University.
            !     Written by
            !                        Ciyou Zhu
            !     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
            !
            !
            !     ************
            integer  :: i
            real(dp) :: gi
            sbgnrm = zero
            do i = 1, n
                gi = g(i)
                if (nbd(i) .ne. 0) then
                    if (gi .lt. zero) then
                        if (nbd(i) .ge. 2) gi = max((x(i)-u(i)),gi)
                    else
                        if (nbd(i) .le. 2) gi = min((x(i)-l(i)),gi)
                    end if
                end if
                sbgnrm = max(sbgnrm,abs(gi))
            end do
        end subroutine projgr

        subroutine subsm ( n, m, nsub, ind, l, u, nbd, x, d, xp, ws, wy, &
            &              theta, xx, gg, &
            &              col, head, iword, wv, wn, info )
            implicit none
            integer  :: n, m, nsub, col, head, iword, info
            integer  :: ind(nsub), nbd(n)
            real(dp) :: theta
            real(dp) :: l(n), u(n), x(n), d(n), xp(n), xx(n), gg(n)
            real(dp) :: ws(n, m), wy(n, m)
            real(dp) :: wv(2*m), wn(2*m, 2*m)
            !     **********************************************************************
            !
            !     This routine contains the major changes in the updated version.
            !     The changes are described in the accompanying paper
            !
            !      Jose Luis Morales, Jorge Nocedal
            !      "Remark On Algorithm 788: L-BFGS-B: Fortran Subroutines for Large-Scale
            !       Bound Constrained Optimization". Decemmber 27, 2010.
            !
            !             J.L. Morales  Departamento de Matematicas,
            !                           Instituto Tecnologico Autonomo de Mexico
            !                           Mexico D.F.
            !
            !             J, Nocedal    Department of Electrical Engineering and
            !                           Computer Science.
            !                           Northwestern University. Evanston, IL. USA
            !
            !                           January 17, 2011
            !
            !      **********************************************************************
            !
            !
            !     Subroutine subsm
            !
            !     Given xcp, l, u, r, an index set that specifies
            !       the active set at xcp, and an l-BFGS matrix B
            !       (in terms of WY, WS, SY, WT, head, col, and theta),
            !       this subroutine computes an approximate solution
            !       of the subspace problem
            !
            !       (P)   min Q(x) = r'(x-xcp) + 1/2 (x-xcp)' B (x-xcp)
            !
            !             subject to l<=x<=u
            !                       x_i=xcp_i for all i in A(xcp)
            !
            !       along the subspace unconstrained Newton direction
            !
            !          d = -(Z'BZ)^(-1) r.
            !
            !       The formula for the Newton direction, given the L-BFGS matrix
            !       and the Sherman-Morrison formula, is
            !
            !          d = (1/theta)r + (1/theta*2) Z'WK^(-1)W'Z r.
            !
            !       where
            !                 K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
            !                     [L_a -R_z           theta*S'AA'S ]
            !
            !     Note that this procedure for computing d differs
            !     from that described in [1]. One can show that the matrix K is
            !     equal to the matrix M^[-1]N in that paper.
            !
            !     n is an integer variable.
            !       On entry n is the dimension of the problem.
            !       On exit n is unchanged.
            !
            !     m is an integer variable.
            !       On entry m is the maximum number of variable metric corrections
            !         used to define the limited memory matrix.
            !       On exit m is unchanged.
            !
            !     nsub is an integer variable.
            !       On entry nsub is the number of free variables.
            !       On exit nsub is unchanged.
            !
            !     ind is an integer array of dimension nsub.
            !       On entry ind specifies the coordinate indices of free variables.
            !       On exit ind is unchanged.
            !
            !     l is a double precision array of dimension n.
            !       On entry l is the lower bound of x.
            !       On exit l is unchanged.
            !
            !     u is a double precision array of dimension n.
            !       On entry u is the upper bound of x.
            !       On exit u is unchanged.
            !
            !     nbd is a integer array of dimension n.
            !       On entry nbd represents the type of bounds imposed on the
            !         variables, and must be specified as follows:
            !         nbd(i)=0 if x(i) is unbounded,
            !                1 if x(i) has only a lower bound,
            !                2 if x(i) has both lower and upper bounds, and
            !                3 if x(i) has only an upper bound.
            !       On exit nbd is unchanged.
            !
            !     x is a double precision array of dimension n.
            !       On entry x specifies the Cauchy point xcp.
            !       On exit x(i) is the minimizer of Q over the subspace of
            !                                                        free variables.
            !
            !     d is a double precision array of dimension n.
            !       On entry d is the reduced gradient of Q at xcp.
            !       On exit d is the Newton direction of Q.
            !
            !    xp is a double precision array of dimension n.
            !       used to safeguard the projected Newton direction
            !
            !    xx is a double precision array of dimension n
            !       On entry it holds the current iterate
            !       On output it is unchanged
            !    gg is a double precision array of dimension n
            !       On entry it holds the gradient at the current iterate
            !       On output it is unchanged
            !
            !     ws and wy are double precision arrays;
            !     theta is a double precision variable;
            !     col is an integer variable;
            !     head is an integer variable.
            !       On entry they store the information defining the
            !                                          limited memory BFGS matrix:
            !         ws(n,m) stores S, a set of s-vectors;
            !         wy(n,m) stores Y, a set of y-vectors;
            !         theta is the scaling factor specifying B_0 = theta I;
            !         col is the number of variable metric corrections stored;
            !         head is the location of the 1st s- (or y-) vector in S (or Y).
            !       On exit they are unchanged.
            !
            !     iword is an integer variable.
            !       On entry iword is unspecified.
            !       On exit iword specifies the status of the subspace solution.
            !         iword = 0 if the solution is in the box,
            !                 1 if some bound is encountered.
            !
            !     wv is a double precision working array of dimension 2m.
            !
            !     wn is a double precision array of dimension 2m x 2m.
            !       On entry the upper triangle of wn stores the LEL^T factorization
            !         of the indefinite matrix
            !
            !              K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
            !                  [L_a -R_z           theta*S'AA'S ]
            !                                                    where E = [-I  0]
            !                                                              [ 0  I]
            !       On exit wn is unchanged.
            !
            !     info is an integer variable.
            !       On entry info is unspecified.
            !       On exit info = 0       for normal return,
            !                    = nonzero for abnormal return
            !                                  when the matrix K is ill-conditioned.
            !
            !     Subprograms called:
            !
            !       Linpack dtrsl.
            !
            !
            !     References:
            !
            !       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
            !       memory algorithm for bound constrained optimization'',
            !       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
            !
            !
            !
            !                           *  *  *
            !
            !     NEOS, November 1994. (Latest revision June 1996.)
            !     Optimization Technology Center.
            !     Argonne National Laboratory and Northwestern University.
            !     Written by
            !                        Ciyou Zhu
            !     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
            !
            !
            !     ************
            integer  :: pointr,m2,col2,ibd,jy,js,i,j,k
            real(dp) :: alpha, xk, dk, temp1, temp2
            real(dp) :: dd_p
            if (nsub .le. 0) return
            ! Compute wv = W'Zd.
            pointr = head
            do i = 1, col
                temp1 = zero
                temp2 = zero
                do j = 1, nsub
                    k = ind(j)
                    temp1 = temp1 + wy(k,pointr)*d(j)
                    temp2 = temp2 + ws(k,pointr)*d(j)
                end do
                wv(i) = temp1
                wv(col + i) = theta*temp2
                pointr = mod(pointr,m) + 1
            end do
            ! Compute wv:=K^(-1)wv.
            m2 = 2*m
            col2 = 2*col
            call dtrsl(wn,m2,col2,wv,11,info)
            if (info .ne. 0) return
            do i = 1, col
                wv(i) = -wv(i)
            end do
            call dtrsl(wn,m2,col2,wv,01,info)
            if (info .ne. 0) return
            ! Compute d = (1/theta)d + (1/theta**2)Z'W wv.
            pointr = head
            do jy = 1, col
                js = col + jy
                do i = 1, nsub
                    k = ind(i)
                    d(i) = d(i) + wy(k,pointr)*wv(jy)/theta + ws(k,pointr)*wv(js)
                end do
                pointr = mod(pointr,m) + 1
            end do
            call dscal( nsub, one/theta, d, 1 )
            !
            !-----------------------------------------------------------------
            !  Let us try the projection, d is the Newton direction
            iword = 0
            call dcopy ( n, x, 1, xp, 1 )
            !
            do i=1, nsub
                k  = ind(i)
                dk = d(i)
                xk = x(k)
                if ( nbd(k) .ne. 0 ) then
                    !
                    if ( nbd(k).eq.1 ) then         ! lower bounds only
                        x(k) = max( l(k), xk + dk )
                        if ( is_equal(x(k),l(k)) ) iword = 1
                    else
                        !
                        if ( nbd(k).eq.2 ) then     ! upper and lower bounds
                            xk   = max( l(k), xk + dk )
                            x(k) = min( u(k), xk )
                            if ( is_equal(x(k),l(k)) .or. is_equal(x(k),u(k)) ) iword = 1
                        else
                            !
                            if ( nbd(k).eq.3 ) then ! upper bounds only
                                x(k) = min( u(k), xk + dk )
                                if ( is_equal(x(k),u(k)) ) iword = 1
                            end if
                        end if
                    end if
                    !
                else                                ! free variables
                    x(k) = xk + dk
                end if
            end do
            !
            if ( iword.eq.0 ) then
                go to 911
            end if
            !
            !     check sign of the directional derivative
            !
            dd_p = zero
            do i=1, n
                dd_p  = dd_p + (x(i) - xx(i))*gg(i)
            end do
            if ( is_gt_zero(dd_p) ) then
                call dcopy( n, xp, 1, x, 1 )
            else
                go to 911
            end if
            !
            !-----------------------------------------------------------------
            !
            alpha = one
            temp1 = alpha
            ibd   = 0
            do i = 1, nsub
                k = ind(i)
                dk = d(i)
                if (nbd(k) .ne. 0) then
                    if (dk .lt. zero .and. nbd(k) .le. 2) then
                        temp2 = l(k) - x(k)
                        if (temp2 .ge. zero) then
                            temp1 = zero
                        else if (dk*alpha .lt. temp2) then
                            temp1 = temp2/dk
                        end if
                    else if (dk .gt. zero .and. nbd(k) .ge. 2) then
                        temp2 = u(k) - x(k)
                        if (temp2 .le. zero) then
                            temp1 = zero
                        else if (dk*alpha .gt. temp2) then
                            temp1 = temp2/dk
                        end if
                    end if
                    if (temp1 .lt. alpha) then
                        alpha = temp1
                        ibd = i
                    end if
                end if
            end do
            if (alpha .lt. one) then
                dk = d(ibd)
                k = ind(ibd)
                if (dk .gt. zero) then
                    x(k) = u(k)
                    d(ibd) = zero
                else if (dk .lt. zero) then
                    x(k) = l(k)
                    d(ibd) = zero
                end if
            end if
            do i = 1, nsub
                k    = ind(i)
                x(k) = x(k) + alpha*d(i)
            end do
911         continue
        end subroutine subsm

        subroutine dcsrch(f,g,stp,ftol,gtol,xtol,stpmin,stpmax, &
            &             task,isave,dsave)
            character(len=*), intent(inout) :: task
            integer       :: isave(2)
            real(dp)      :: f,g,stp,ftol,gtol,xtol,stpmin,stpmax
            real(dp)      :: dsave(13)
            !     **********
            !
            !     Subroutine dcsrch
            !
            !     This subroutine finds a step that satisfies a sufficient
            !     decrease condition and a curvature condition.
            !
            !     Each call of the subroutine updates an interval with
            !     endpoints stx and sty. The interval is initially chosen
            !     so that it contains a minimizer of the modified function
            !
            !           psi(stp) = f(stp) - f(0) - ftol*stp*f'(0).
            !
            !     If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
            !     interval is chosen so that it contains a minimizer of f.
            !
            !     The algorithm is designed to find a step that satisfies
            !     the sufficient decrease condition
            !
            !           f(stp) <= f(0) + ftol*stp*f'(0),
            !
            !     and the curvature condition
            !
            !           abs(f'(stp)) <= gtol*abs(f'(0)).
            !
            !     If ftol is less than gtol and if, for example, the function
            !     is bounded below, then there is always a step which satisfies
            !     both conditions.
            !
            !     If no step can be found that satisfies both conditions, then
            !     the algorithm stops with a warning. In this case stp only
            !     satisfies the sufficient decrease condition.
            !
            !     A typical invocation of dcsrch has the following outline:
            !
            !     task = 'START'
            !  10 continue
            !        call dcsrch( ... )
            !        if (task .eq. 'FG') then
            !           Evaluate the function and the gradient at stp
            !           goto 10
            !           end if
            !
            !     NOTE: The user must no alter work arrays between calls.
            !
            !     The subroutine statement is
            !
            !        subroutine dcsrch(f,g,stp,ftol,gtol,xtol,stpmin,stpmax,
            !                          task,isave,dsave)
            !     where
            !
            !       f is a double precision variable.
            !         On initial entry f is the value of the function at 0.
            !            On subsequent entries f is the value of the
            !            function at stp.
            !         On exit f is the value of the function at stp.
            !
            !       g is a double precision variable.
            !         On initial entry g is the derivative of the function at 0.
            !            On subsequent entries g is the derivative of the
            !            function at stp.
            !         On exit g is the derivative of the function at stp.
            !
            !       stp is a double precision variable.
            !         On entry stp is the current estimate of a satisfactory
            !            step. On initial entry, a positive initial estimate
            !            must be provided.
            !         On exit stp is the current estimate of a satisfactory step
            !            if task = 'FG'. If task = 'CONV' then stp satisfies
            !            the sufficient decrease and curvature condition.
            !
            !       ftol is a double precision variable.
            !         On entry ftol specifies a nonnegative tolerance for the
            !            sufficient decrease condition.
            !         On exit ftol is unchanged.
            !
            !       gtol is a double precision variable.
            !         On entry gtol specifies a nonnegative tolerance for the
            !            curvature condition.
            !         On exit gtol is unchanged.
            !
            !       xtol is a double precision variable.
            !         On entry xtol specifies a nonnegative relative tolerance
            !            for an acceptable step. The subroutine exits with a
            !            warning if the relative difference between sty and stx
            !            is less than xtol.
            !         On exit xtol is unchanged.
            !
            !       stpmin is a double precision variable.
            !         On entry stpmin is a nonnegative lower bound for the step.
            !         On exit stpmin is unchanged.
            !
            !       stpmax is a double precision variable.
            !         On entry stpmax is a nonnegative upper bound for the step.
            !         On exit stpmax is unchanged.
            !
            !       task is a character variable of length at least 60.
            !         On initial entry task must be set to 'START'.
            !         On exit task indicates the required action:
            !
            !            If task(1:2) = 'FG' then evaluate the function and
            !            derivative at stp and call dcsrch again.
            !
            !            If task(1:4) = 'CONV' then the search is successful.
            !
            !            If task(1:4) = 'WARN' then the subroutine is not able
            !            to satisfy the convergence conditions. The exit value of
            !            stp contains the best point found during the search.
            !
            !            If task(1:5) = 'ERROR' then there is an error in the
            !            input arguments.
            !
            !         On exit with convergence, a warning or an error, the
            !            variable task contains additional information.
            !
            !       isave is an integer work array of dimension 2.
            !
            !       dsave is a double precision work array of dimension 13.
            !
            !     Subprograms called
            !
            !       MINPACK-2 ... dcstep
            !
            !     MINPACK-1 Project. June 1983.
            !     Argonne National Laboratory.
            !     Jorge J. More' and David J. Thuente.
            !
            !     MINPACK-2 Project. October 1993.
            !     Argonne National Laboratory and University of Minnesota.
            !     Brett M. Averick, Richard G. Carter, and Jorge J. More'.
            !
            !     **********
            logical  :: brackt
            integer  :: stage
            real(dp) :: finit,ftest,fm,fx,fxm,fy,fym,ginit,gtest
            real(dp) :: gm,gx,gxm,gy,gym,stx,sty,stmin,stmax,width,width1
            !     Initialization block.
            if (task(1:5) .eq. 'START') then
                !        Check the input arguments for errors.
                if (stp .lt. stpmin) task = 'ERROR: STP .LT. STPMIN'
                if (stp .gt. stpmax) task = 'ERROR: STP .GT. STPMAX'
                if (g .ge. zero) task = 'ERROR: INITIAL G .GE. ZERO'
                if (ftol .lt. zero) task = 'ERROR: FTOL .LT. ZERO'
                if (gtol .lt. zero) task = 'ERROR: GTOL .LT. ZERO'
                if (xtol .lt. zero) task = 'ERROR: XTOL .LT. ZERO'
                if (stpmin .lt. zero) task = 'ERROR: STPMIN .LT. ZERO'
                if (stpmax .lt. stpmin) task = 'ERROR: STPMAX .LT. STPMIN'
                !        Exit if there are errors on input.
                if (task(1:5) .eq. 'ERROR') return
                !        Initialize local variables.
                brackt = .false.
                stage = 1
                finit = f
                ginit = g
                gtest = ftol*ginit
                width = stpmax - stpmin
                width1 = width/p5
                ! The variables stx, fx, gx contain the values of the step,
                ! function, and derivative at the best step.
                ! The variables sty, fy, gy contain the value of the step,
                ! function, and derivative at sty.
                ! The variables stp, f, g contain the values of the step,
                ! function, and derivative at stp.
                stx = zero
                fx = finit
                gx = ginit
                sty = zero
                fy = finit
                gy = ginit
                stmin = zero
                stmax = stp + xtrapu*stp
                task = 'FG'
                goto 1000
            else
                !        Restore local variables.
                if (isave(1) .eq. 1) then
                    brackt = .true.
                else
                    brackt = .false.
                end if
                stage = isave(2)
                ginit = dsave(1)
                gtest = dsave(2)
                gx = dsave(3)
                gy = dsave(4)
                finit = dsave(5)
                fx = dsave(6)
                fy = dsave(7)
                stx = dsave(8)
                sty = dsave(9)
                stmin = dsave(10)
                stmax = dsave(11)
                width = dsave(12)
                width1 = dsave(13)
            end if
            ! If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
            ! algorithm enters the second stage.
            ftest = finit + stp*gtest
            if (stage .eq. 1 .and. f .le. ftest .and. g .ge. zero) stage = 2
            ! Test for warnings.
            if (brackt .and. (stp .le. stmin .or. stp .ge. stmax)) task = 'WARNING: ROUNDING ERRORS PREVENT PROGRESS'
            if (brackt .and. stmax - stmin .le. xtol*stmax) task = 'WARNING: XTOL TEST SATISFIED'
            if ((is_equal(stp ,stpmax)) .and. f .le. ftest .and. g .le. gtest) task = 'WARNING: STP = STPMAX'
            if ((is_equal(stp ,stpmin)) .and. (f .gt. ftest .or. g .ge. gtest)) task = 'WARNING: STP = STPMIN'
            ! Test for convergence.
            if (f .le. ftest .and. abs(g) .le. gtol*(-ginit)) task = 'CONVERGENCE'
            ! Test for termination.
            if (task(1:4) .eq. 'WARN' .or. task(1:4) .eq. 'CONV') goto 1000
            ! A modified function is used to predict the step during the
            ! first stage if a lower function value has been obtained but
            ! the decrease is not sufficient.
            if (stage .eq. 1 .and. f .le. fx .and. f .gt. ftest) then
                ! Define the modified function and derivative values.
                fm = f - stp*gtest
                fxm = fx - stx*gtest
                fym = fy - sty*gtest
                gm = g - gtest
                gxm = gx - gtest
                gym = gy - gtest
                ! Call dcstep to update stx, sty, and to compute the new step.
                call dcstep(stx,fxm,gxm,sty,fym,gym,stp,fm,gm,brackt,stmin,stmax)
                ! Reset the function and derivative values for f.
                fx = fxm + stx*gtest
                fy = fym + sty*gtest
                gx = gxm + gtest
                gy = gym + gtest
            else
                ! Call dcstep to update stx, sty, and to compute the new step.
                call dcstep(stx,fx,gx,sty,fy,gy,stp,f,g,brackt,stmin,stmax)
            end if
            ! Decide if a bisection step is needed.
            if (brackt) then
                if (abs(sty-stx) .ge. p66*width1) stp = stx + p5*(sty - stx)
                width1 = width
                width = abs(sty-stx)
            end if
            ! Set the minimum and maximum steps allowed for stp.
            if (brackt) then
                stmin = min(stx,sty)
                stmax = max(stx,sty)
            else
                stmin = stp + xtrapl*(stp - stx)
                stmax = stp + xtrapu*(stp - stx)
            end if
            ! Force the step to be within the bounds stpmax and stpmin.
            stp = max(stp,stpmin)
            stp = min(stp,stpmax)
            ! If further progress is not possible, let stp be the best
            ! point obtained during the search.
            if (brackt .and. (stp .le. stmin .or. stp .ge. stmax) &
                &   .or. (brackt .and. stmax-stmin .le. xtol*stmax)) stp = stx
            ! Obtain another function and derivative.
            task = 'FG'
1000        continue
            ! Save local variables.
            if (brackt) then
                isave(1) = 1
            else
                isave(1) = 0
            end if
            isave(2) = stage
            dsave(1) =  ginit
            dsave(2) =  gtest
            dsave(3) =  gx
            dsave(4) =  gy
            dsave(5) =  finit
            dsave(6) =  fx
            dsave(7) =  fy
            dsave(8) =  stx
            dsave(9) =  sty
            dsave(10) = stmin
            dsave(11) = stmax
            dsave(12) = width
            dsave(13) = width1
        end subroutine dcsrch

        subroutine dcstep(stx,fx,dx,sty,fy,dy,stp,fp,dp_var,brackt,stpmin,stpmax)
            logical  :: brackt
            real(dp) :: stx,fx,dx,sty,fy,dy,stp,fp,dp_var,stpmin,stpmax
            !     **********
            !
            !     Subroutine dcstep
            !
            !     This subroutine computes a safeguarded step for a search
            !     procedure and updates an interval that contains a step that
            !     satisfies a sufficient decrease and a curvature condition.
            !
            !     The parameter stx contains the step with the least function
            !     value. If brackt is set to .true. then a minimizer has
            !     been bracketed in an interval with endpoints stx and sty.
            !     The parameter stp contains the current step.
            !     The subroutine assumes that if brackt is set to .true. then
            !
            !           min(stx,sty) < stp < max(stx,sty),
            !
            !     and that the derivative at stx is negative in the direction
            !     of the step.
            !
            !     The subroutine statement is
            !
            !       subroutine dcstep(stx,fx,dx,sty,fy,dy,stp,fp,dp_var,brackt,
            !                         stpmin,stpmax)
            !
            !     where
            !
            !       stx is a double precision variable.
            !         On entry stx is the best step obtained so far and is an
            !            endpoint of the interval that contains the minimizer.
            !         On exit stx is the updated best step.
            !
            !       fx is a double precision variable.
            !         On entry fx is the function at stx.
            !         On exit fx is the function at stx.
            !
            !       dx is a double precision variable.
            !         On entry dx is the derivative of the function at
            !            stx. The derivative must be negative in the direction of
            !            the step, that is, dx and stp - stx must have opposite
            !            signs.
            !         On exit dx is the derivative of the function at stx.
            !
            !       sty is a double precision variable.
            !         On entry sty is the second endpoint of the interval that
            !            contains the minimizer.
            !         On exit sty is the updated endpoint of the interval that
            !            contains the minimizer.
            !
            !       fy is a double precision variable.
            !         On entry fy is the function at sty.
            !         On exit fy is the function at sty.
            !
            !       dy is a double precision variable.
            !         On entry dy is the derivative of the function at sty.
            !         On exit dy is the derivative of the function at the exit sty.
            !
            !       stp is a double precision variable.
            !         On entry stp is the current step. If brackt is set to .true.
            !            then on input stp must be between stx and sty.
            !         On exit stp is a new trial step.
            !
            !       fp is a double precision variable.
            !         On entry fp is the function at stp
            !         On exit fp is unchanged.
            !
            !       dp_var is a double precision variable.
            !         On entry dp_var is the the derivative of the function at stp.
            !         On exit dp_var is unchanged.
            !
            !       brackt is an logical variable.
            !         On entry brackt specifies if a minimizer has been bracketed.
            !            Initially brackt must be set to .false.
            !         On exit brackt specifies if a minimizer has been bracketed.
            !            When a minimizer is bracketed brackt is set to .true.
            !
            !       stpmin is a double precision variable.
            !         On entry stpmin is a lower bound for the step.
            !         On exit stpmin is unchanged.
            !
            !       stpmax is a double precision variable.
            !         On entry stpmax is an upper bound for the step.
            !         On exit stpmax is unchanged.
            !
            !     MINPACK-1 Project. June 1983
            !     Argonne National Laboratory.
            !     Jorge J. More' and David J. Thuente.
            !
            !     MINPACK-2 Project. October 1993.
            !     Argonne National Laboratory and University of Minnesota.
            !     Brett M. Averick and Jorge J. More'.
            !
            !     **********
            real(dp) :: gamma,p,q,r,s,sgnd,stpc,stpf,stpq,theta
            sgnd = dp_var*(dx/abs(dx))
            ! First case: A higher function value. The minimum is bracketed.
            ! If the cubic step is closer to stx than the quadratic step, the
            ! cubic step is taken, otherwise the average of the cubic and
            ! quadratic steps is taken.
            if (fp .gt. fx) then
                theta = three*(fx - fp)/(stp - stx) + dx + dp_var
                s = max(abs(theta),abs(dx),abs(dp_var))
                gamma = s*sqrt((theta/s)**2 - (dx/s)*(dp_var/s))
                if (stp .lt. stx) gamma = -gamma
                p = (gamma - dx) + theta
                q = ((gamma - dx) + gamma) + dp_var
                r = p/q
                stpc = stx + r*(stp - stx)
                stpq = stx + ((dx/((fx - fp)/(stp - stx) + dx))/two)*(stp - stx)
                if (abs(stpc-stx) .lt. abs(stpq-stx)) then
                    stpf = stpc
                else
                    stpf = stpc + (stpq - stpc)/two
                end if
                brackt = .true.
                ! Second case: A lower function value and derivatives of opposite
                ! sign. The minimum is bracketed. If the cubic step is farther from
                ! stp than the secant step, the cubic step is taken, otherwise the
                ! secant step is taken.
            else if (sgnd .lt. zero) then
                theta = three*(fx - fp)/(stp - stx) + dx + dp_var
                s = max(abs(theta),abs(dx),abs(dp_var))
                gamma = s*sqrt((theta/s)**2 - (dx/s)*(dp_var/s))
                if (stp .gt. stx) gamma = -gamma
                p = (gamma - dp_var) + theta
                q = ((gamma - dp_var) + gamma) + dx
                r = p/q
                stpc = stp + r*(stx - stp)
                stpq = stp + (dp_var/(dp_var - dx))*(stx - stp)
                if (abs(stpc-stp) .gt. abs(stpq-stp)) then
                    stpf = stpc
                else
                    stpf = stpq
                end if
                brackt = .true.
                ! Third case: A lower function value, derivatives of the same sign,
                ! and the magnitude of the derivative decreases.
            else if (abs(dp_var) .lt. abs(dx)) then
                ! The cubic step is computed only if the cubic tends to infinity
                ! in the direction of the step or if the minimum of the cubic
                ! is beyond stp. Otherwise the cubic step is defined to be the
                ! secant step.
                theta = three*(fx - fp)/(stp - stx) + dx + dp_var
                s = max(abs(theta),abs(dx),abs(dp_var))
                ! The case gamma = 0 only arises if the cubic does not tend
                ! to infinity in the direction of the step.
                gamma = s*sqrt(max(zero,(theta/s)**2-(dx/s)*(dp_var/s)))
                if (stp .gt. stx) gamma = -gamma
                p = (gamma - dp_var) + theta
                q = (gamma + (dx - dp_var)) + gamma
                r = p/q
                if (r .lt. zero .and. (.not. is_zero(gamma))) then  !gamma .ne. zero
                    stpc = stp + r*(stx - stp)
                else if (stp .gt. stx) then
                    stpc = stpmax
                else
                    stpc = stpmin
                end if
                stpq = stp + (dp_var/(dp_var - dx))*(stx - stp)
                if (brackt) then
                    ! A minimizer has been bracketed. If the cubic step is
                    ! closer to stp than the secant step, the cubic step is
                    ! taken, otherwise the secant step is taken.
                    if (abs(stpc-stp) .lt. abs(stpq-stp)) then
                        stpf = stpc
                    else
                        stpf = stpq
                    end if
                    if (stp .gt. stx) then
                        stpf = min(stp+p66*(sty-stp),stpf)
                    else
                        stpf = max(stp+p66*(sty-stp),stpf)
                    end if
                else
                    ! A minimizer has not been bracketed. If the cubic step is
                    ! farther from stp than the secant step, the cubic step is
                    ! taken, otherwise the secant step is taken.
                    if (abs(stpc-stp) .gt. abs(stpq-stp)) then
                        stpf = stpc
                    else
                        stpf = stpq
                    end if
                    stpf = min(stpmax,stpf)
                    stpf = max(stpmin,stpf)
                end if
                ! Fourth case: A lower function value, derivatives of the same sign,
                ! and the magnitude of the derivative does not decrease. If the
                ! minimum is not bracketed, the step is either stpmin or stpmax,
                ! otherwise the cubic step is taken.
            else
                if (brackt) then
                    theta = three*(fp - fy)/(sty - stp) + dy + dp_var
                    s = max(abs(theta),abs(dy),abs(dp_var))
                    gamma = s*sqrt((theta/s)**2 - (dy/s)*(dp_var/s))
                    if (stp .gt. sty) gamma = -gamma
                    p = (gamma - dp_var) + theta
                    q = ((gamma - dp_var) + gamma) + dy
                    r = p/q
                    stpc = stp + r*(sty - stp)
                    stpf = stpc
                else if (stp .gt. stx) then
                    stpf = stpmax
                else
                    stpf = stpmin
                end if
            end if
            ! Update the interval which contains a minimizer
            if (fp .gt. fx) then
                sty = stp
                fy = fp
                dy = dp_var
            else
                if (sgnd .lt. zero) then
                    sty = stx
                    fy = fx
                    dy = dx
                end if
                stx = stp
                fx = fp
                dx = dp_var
            end if
            ! Compute the new step
            stp = stpf
        end subroutine dcstep

        function dnrm2(n,x,incx) result(res)
            integer,  intent(in) :: n,incx
            real(dp), intent(in) :: x(n)
            real(dp)             :: res
            !     **********
            !
            !     Function dnrm2
            !
            !     Given a vector x of length n, this function calculates the
            !     Euclidean norm of x with stride incx.
            !
            !     The function statement is
            !
            !       double precision function dnrm2(n,x,incx)
            !
            !     where
            !
            !       n is a positive integer input variable.
            !
            !       x is an input array of length n.
            !
            !       incx is a positive integer variable that specifies the
            !         stride of the vector.
            !
            !     Subprograms called
            !
            !       FORTRAN-supplied ... abs, max, sqrt
            !
            !     MINPACK-2 Project. February 1991.
            !     Argonne National Laboratory.
            !     Brett M. Averick.
            !
            !     **********
            integer  :: i
            real(dp) :: scale
            res = 0.0d0
            scale = 0.0d0
            do i = 1, n, incx
                scale = max(scale, abs(x(i)))
            end do
            if (is_zero(scale)) return
            do i = 1, n, incx
                res = res + (x(i)/scale)**2
            end do
            res = scale*sqrt(res)
        end function dnrm2

        subroutine daxpy(n,da,dx,incx,dy,incy)
            !
            !     constant times a vector plus a vector.
            !     uses unrolled loops for increments equal to one.
            !     jack dongarra, linpack, 3/11/78.
            !
            real(dp), intent(in)  :: dx(*),da
            real(dp), intent(out) :: dy(*)
            integer,  intent(in)  :: incx,incy,n
            integer               :: i,ix,iy,m,mp1
            !
            if (n.le.0) return
            if (is_zero(da) ) return
            if (incx.eq.1.and.incy.eq.1) go to 20
            !
            !        code for unequal increments or equal increments
            !          not equal to 1
            !
            ix = 1
            iy = 1
            if(incx.lt.0)ix = (-n+1)*incx + 1
            if(incy.lt.0)iy = (-n+1)*incy + 1
            do i = 1,n
                dy(iy) = dy(iy) + da*dx(ix)
                ix = ix + incx
                iy = iy + incy
            end do
            return
            !
            !        code for both increments equal to 1
            !        clean-up loop
            !
20          m = mod(n,4)
            if( m .eq. 0 ) go to 40
            do i = 1,m
                dy(i) = dy(i) + da*dx(i)
            end do
            if( n .lt. 4 ) return
40          mp1 = m + 1
            do i = mp1,n,4
                dy(i) = dy(i) + da*dx(i)
                dy(i + 1) = dy(i + 1) + da*dx(i + 1)
                dy(i + 2) = dy(i + 2) + da*dx(i + 2)
                dy(i + 3) = dy(i + 3) + da*dx(i + 3)
            end do
        end subroutine daxpy

        subroutine dcopy(n,dx,incx,dy,incy)
            !
            !     copies a vector, x, to a vector, y.
            !     uses unrolled loops for increments equal to one.
            !     jack dongarra, linpack, 3/11/78.
            !
            real(dp), intent(in)  :: dx(*)
            real(dp), intent(out) :: dy(*)
            integer,  intent(in)  :: incx,incy,n
            integer               :: i,ix,iy,m,mp1
            !
            if (n.le.0) return
            if (incx.eq.1.and.incy.eq.1) go to 20
            !
            !        code for unequal increments or equal increments
            !          not equal to 1
            !
            ix = 1
            iy = 1
            if(incx.lt.0)ix = (-n+1)*incx + 1
            if(incy.lt.0)iy = (-n+1)*incy + 1
            do i = 1,n
                dy(iy) = dx(ix)
                ix = ix + incx
                iy = iy + incy
            end do
            return
            !
            !        code for both increments equal to 1
            !
            !        clean-up loop
            !
20          m = mod(n,7)
            if( m .eq. 0 ) go to 40
            do i = 1,m
                dy(i) = dx(i)
            end do
            if( n .lt. 7 ) return
40          mp1 = m + 1
            do i = mp1,n,7
                dy(i) = dx(i)
                dy(i + 1) = dx(i + 1)
                dy(i + 2) = dx(i + 2)
                dy(i + 3) = dx(i + 3)
                dy(i + 4) = dx(i + 4)
                dy(i + 5) = dx(i + 5)
                dy(i + 6) = dx(i + 6)
            end do
        end subroutine dcopy

        function ddot_lbfgsb(n,dx,incx,dy,incy) result(res)
            !
            !     forms the dot product of two vectors.
            !     uses unrolled loops for increments equal to one.
            !     jack dongarra, linpack, 3/11/78.
            !
            real(dp), intent(in) :: dx(*),dy(*)
            integer,  intent(in) :: n,incx,incy
            real(dp)             :: res
            real(dp)             :: dtemp
            integer              :: i,ix,iy,m,mp1
            !
            res = 0.0d0
            dtemp = 0.0d0
            if (n.le.0) return
            if (incx.eq.1.and.incy.eq.1) go to 20
            !
            !        code for unequal increments or equal increments
            !          not equal to 1
            !
            ix = 1
            iy = 1
            if(incx.lt.0)ix = (-n+1)*incx + 1
            if(incy.lt.0)iy = (-n+1)*incy + 1
            do i = 1,n
                dtemp = dtemp + dx(ix)*dy(iy)
                ix = ix + incx
                iy = iy + incy
            end do
            res = dtemp
            return
            !
            !        code for both increments equal to 1
            !
            !        clean-up loop
            !
20          m = mod(n,5)
            if( m .eq. 0 ) go to 40
            do i = 1,m
                dtemp = dtemp + dx(i)*dy(i)
            end do
            if( n .lt. 5 ) go to 60
40          mp1 = m + 1
            do i = mp1,n,5
                dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) + &
                    &   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
            end do
60          res = dtemp
        end function ddot_lbfgsb

        subroutine  dscal(n,da,dx,incx)
            !
            !     scales a vector by a constant.
            !     uses unrolled loops for increment equal to one.
            !     jack dongarra, linpack, 3/11/78.
            !     modified 3/93 to return if incx .le. 0.
            !
            real(dp), intent(in)    :: da
            real(dp), intent(inout) :: dx(*)
            integer,  intent(in)    :: n,incx
            integer                 :: i,m,mp1,nincx
            !
            if ( n.le.0 .or. incx.le.0 ) return
            if (incx.eq.1) go to 20
            !
            !        code for increment not equal to 1
            !
            nincx = n*incx
            do i = 1,nincx,incx
                dx(i) = da*dx(i)
            end do
            return
            !
            !        code for increment equal to 1
            !
            !        clean-up loop
            !
20          m = mod(n,5)
            if( m .eq. 0 ) go to 40
            do i = 1,m
                dx(i) = da*dx(i)
            end do
            if( n .lt. 5 ) return
40          mp1 = m + 1
            do i = mp1,n,5
                dx(i) = da*dx(i)
                dx(i + 1) = da*dx(i + 1)
                dx(i + 2) = da*dx(i + 2)
                dx(i + 3) = da*dx(i + 3)
                dx(i + 4) = da*dx(i + 4)
            end do
        end subroutine dscal

        subroutine dpofa(a,lda,n,info)
            integer  :: lda,n,info
            real(dp) :: a(lda,*)
            !
            !     dpofa factors a double precision symmetric positive definite
            !     matrix.
            !
            !     dpofa is usually called by dpoco, but it can be called
            !     directly with a saving in time if  rcond  is not needed.
            !     (time for dpoco) = (1 + 18/n)*(time for dpofa) .
            !
            !     on entry
            !
            !        a       double precision(lda, n)
            !                the symmetric matrix to be factored.  only the
            !                diagonal and upper triangle are used.
            !
            !        lda     integer
            !                the leading dimension of the array  a .
            !
            !        n       integer
            !                the order of the matrix  a .
            !
            !     on return
            !
            !        a       an upper triangular matrix  r  so that  a = trans(r)*r
            !                where  trans(r)  is the transpose.
            !                the strict lower triangle is unaltered.
            !                if  info .ne. 0 , the factorization is not complete.
            !
            !        info    integer
            !                = 0  for normal return.
            !                = k  signals an error condition.  the leading minor
            !                     of order  k  is not positive definite.
            !
            !     linpack.  this version dated 08/14/78 .
            !     cleve moler, university of new mexico, argonne national lab.
            !
            !     subroutines and functions
            !
            !     blas ddot
            !     fortran sqrt
            !
            !     internal variables
            !
            real(dp) :: t
            real(dp) :: s
            integer  :: j,jm1,k
            !     begin block with ...exits to 40
            !
            do j = 1, n
                info = j
                s = 0.0d0
                jm1 = j - 1
                if (jm1 .lt. 1) go to 20
                do k = 1, jm1
                    t = a(k,j) - ddot_lbfgsb(k-1,a(1,k),1,a(1,j),1)
                    t = t/a(k,k)
                    a(k,j) = t
                    s = s + t*t
                end do
20              continue
                s = a(j,j) - s
                !     ......exit
                if (s .le. 0.0d0) go to 40
                a(j,j) = sqrt(s)
            end do
            info = 0
40          continue
        end subroutine dpofa

        subroutine dtrsl(t,ldt,n,b,job,info)
            integer  :: ldt,n,job,info
            real(dp) :: t(ldt,*),b(*)
            !
            !     dtrsl solves systems of the form
            !
            !                   t * x = b
            !     or
            !                   trans(t) * x = b
            !
            !     where t is a triangular matrix of order n. here trans(t)
            !     denotes the transpose of the matrix t.
            !
            !     on entry
            !
            !         t         double precision(ldt,n)
            !                   t contains the matrix of the system. the zero
            !                   elements of the matrix are not referenced, and
            !                   the corresponding elements of the array can be
            !                   used to store other information.
            !
            !         ldt       integer
            !                   ldt is the leading dimension of the array t.
            !
            !         n         integer
            !                   n is the order of the system.
            !
            !         b         double precision(n).
            !                   b contains the right hand side of the system.
            !
            !         job       integer
            !                   job specifies what kind of system is to be solved.
            !                   if job is
            !
            !                        00   solve t*x=b, t lower triangular,
            !                        01   solve t*x=b, t upper triangular,
            !                        10   solve trans(t)*x=b, t lower triangular,
            !                        11   solve trans(t)*x=b, t upper triangular.
            !
            !     on return
            !
            !         b         b contains the solution, if info .eq. 0.
            !                   otherwise b is unaltered.
            !
            !         info      integer
            !                   info contains zero if the system is nonsingular.
            !                   otherwise info contains the index of
            !                   the first zero diagonal element of t.
            !
            !     linpack. this version dated 08/14/78 .
            !     g. w. stewart, university of maryland, argonne national lab.
            !
            !     subroutines and functions
            !
            !     blas daxpy,ddot
            !     fortran mod
            !
            !     internal variables
            !
            real(dp) :: temp
            integer  :: cas,j,jj
            !
            !     begin block permitting ...exits to 150
            !
            !        check for zero diagonal elements.
            !
            do info = 1, n
                !     ......exit
                if (is_zero(t(info,info))) return
            end do
            info = 0
            !
            !        determine the task and go to it.
            !
            cas = 1
            if (mod(job,10) .ne. 0) cas = 2
            if (mod(job,100)/10 .ne. 0) cas = cas + 2
            select  case(cas)
            !
            !        solve t*x=b for t lower triangular
            !
            case(1)
                b(1) = b(1)/t(1,1)
                if (n .lt. 2) return
                do j = 2, n
                    temp = -b(j-1)
                    call daxpy(n-j+1,temp,t(j,j-1),1,b(j),1)
                    b(j) = b(j)/t(j,j)
                end do
            !
            !        solve t*x=b for t upper triangular.
            !
            case(2)
                b(n) = b(n)/t(n,n)
                if (n .lt. 2) return
                do jj = 2, n
                    j = n - jj + 1
                    temp = -b(j+1)
                    call daxpy(j,temp,t(1,j+1),1,b(1),1)
                    b(j) = b(j)/t(j,j)
                end do
            !
            !        solve trans(t)*x=b for t lower triangular.
            !
            case(3)
                b(n) = b(n)/t(n,n)
                if (n .lt. 2)return
                do jj = 2, n
                    j = n - jj + 1
                    b(j) = b(j) - ddot_lbfgsb(jj-1,t(j+1,j),1,b(j+1),1)
                    b(j) = b(j)/t(j,j)
                end do
            !
            !        solve trans(t)*x=b for t upper triangular.
            !
            case(4)
                b(1) = b(1)/t(1,1)
                if (n .lt. 2) return
                do j = 2, n
                    b(j) = b(j) - ddot_lbfgsb(j-1,t(1,j),1,b(1),1)
                    b(j) = b(j)/t(j,j)
                end do
            end select

        end subroutine dtrsl

    end subroutine lbfgsb_minimize

    !> \brief  is a destructor
    subroutine kill_opt_lbfgsb( self )
        class(opt_lbfgsb), intent(inout) :: self !< instance
        if( self%exists )then
            deallocate(self%iwa, self%l, self%u, self%g, self%wa )
            self%exists = .false.
        endif
    end subroutine kill_opt_lbfgsb

end module simple_opt_lbfgsb
