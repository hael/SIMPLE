!@descr: atomicity prior of the PCG reconstruction, lambda_A x^T (I - P_G) x with P_G the projector onto atom kernels
! The map is modelled as x = G c + e and only e is penalized. The columns of G are Gaussian atom kernels clipped
! to the solve support; G^T G is sparse (only overlapping footprints couple) and solved by Jacobi-preconditioned CG
! in double precision. Design: doc/implementation_notes/planned/single_pcg_integration.md.
module simple_pcg_atomic_prior
use simple_core_module_api, only: dp, PI, simple_exception
use simple_pcg_solver,      only: pcg_operator
use simple_image,           only: image
implicit none
private
#include "simple_local_flags.inc"

public :: pcg_atomic_prior, pcg_atomic_prior_stats

real,     parameter :: SIGMA_MIN_VOX  = 1.0        !< width of the narrowest kernel (B_min), voxels: the lattice samples it
real,     parameter :: FOOT_NSIGMA    = 3.0        !< footprint half-width in kernel widths
integer,  parameter :: FOOT_RAD_MIN   = 4          !< footprint half-width bounds, voxels
integer,  parameter :: FOOT_RAD_MAX   = 8
real,     parameter :: MIN_SEP_A      = 0.5        !< an atom closer than this (A) to an earlier one is a duplicate
real(dp), parameter :: GRAM_RIDGE_REL = 1.0e-7_dp  !< ridge of the Gram solve relative to each diagonal entry
real(dp), parameter :: GRAM_TOL       = 1.0e-10_dp !< relative residual of the Gram solve
integer,  parameter :: GRAM_MAXITS    = 2000
real,     parameter :: WIDTH_FIT_CORR = 0.5        !< shells where model and map correlate at least this enter the width fit
integer,  parameter :: WIDTH_FIT_MINSHELLS = 3
real,     parameter :: FIT_SIGMA_VOX  = 1.5        !< width of the width fit's model: no aliasing below Nyquist
real,     parameter :: FIT_NSIGMA     = 5.0        !< its footprint half-width in widths: no truncation plateau
real,     parameter :: C2_MAX         = 0.98       !< cap of the explained signal fraction of the strength estimate
real,     parameter :: SIGNAL_FLOOR   = 1.0e-3     !< signal power floor of the strength estimate, in noise units

!> what the prior was built from and the moments of its strength estimate
type :: pcg_atomic_prior_stats
    integer :: natoms_in     = 0   !< atoms given
    integer :: natoms        = 0   !< atoms kept
    integer :: ndup          = 0   !< dropped as duplicates (closer than MIN_SEP_A)
    integer :: noutside      = 0   !< dropped with their centre voxel outside the box or the support
    integer :: nsupport      = 0   !< voxels of the solve support
    integer :: foot_rad      = 0   !< footprint half-width, voxels
    integer :: gram_nnz      = 0   !< non-zeros of the Gram matrix
    integer :: nshells_fit   = 0   !< shells in the width fit
    real    :: bfac_min      = 0.  !< B of the narrowest kernel (A^2)
    real    :: bfac_data     = 0.  !< B added by the width fit or given (A^2)
    real    :: sigma_vox     = 0.  !< kernel width, voxels
    real    :: noise         = 0.  !< N, noise variance per voxel of the merged base map
    real    :: signal        = 0.  !< S, signal power per voxel
    real    :: explained     = 0.  !< E, power of P_G x_b per voxel
    real    :: c2            = 0.  !< C^2, explained fraction of the signal
    real    :: lambda_rel_merged = 0. !< (N/S)/(1-C^2), the ratio for a solve of the merged data
    real    :: lambda_rel_est    = 0. !< the same with the half-map noise 2N, the ratio for a half solve
end type pcg_atomic_prior_stats

type :: pcg_atomic_prior
    private
    integer :: box    = 0
    integer :: natoms = 0
    integer :: rad    = 0      !< footprint half-width
    integer :: nside  = 0      !< footprint cube side, 2*rad+1
    real    :: smpd   = 1.
    real    :: bfac   = 0.     !< B of the current kernels (A^2)
    real    :: sigma  = 0.     !< width of the current kernels (voxels)
    real,     allocatable :: pos(:,:)        !< (3,natoms) atom positions, voxels (1-based)
    real,     allocatable :: amp(:)          !< column scale of each kernel
    integer,  allocatable :: corner(:,:)     !< (3,natoms) lowest voxel of each footprint cube
    real,     allocatable :: fval(:,:)       !< (nside**3,natoms) kernel values, zero off the support and box
    logical,  allocatable :: support(:,:,:)  !< the solve support
    integer,  allocatable :: gptr(:)         !< CSR row pointers of G^T G
    integer,  allocatable :: gcol(:)         !< CSR columns
    real(dp), allocatable :: gval(:)         !< CSR values
    real(dp), allocatable :: gdiag(:)        !< diagonal of G^T G
    type(pcg_atomic_prior_stats) :: stats
    logical :: exists = .false.
  contains
    procedure :: new
    procedure :: build
    procedure :: gather
    procedure :: scatter
    procedure :: solve_gram
    procedure :: project
    procedure :: apply
    procedure :: amplitudes
    procedure :: fit_width
    procedure :: estimate_lambda
    procedure :: get_stats
    procedure :: get_natoms
    procedure :: get_box
    procedure :: get_bfac
    procedure :: get_bfac_min
    procedure :: is_built
    procedure :: kill
    procedure, private :: build_gram
    procedure, private :: gram_matvec
    procedure, private :: overlap
    procedure, private :: scatter_add
    procedure, private :: model_map
end type pcg_atomic_prior

contains

    !> The model: atom positions in voxels (1-based, the solve box), the solve support (> 0) and optional
    !! column scales; duplicates and atoms centred off the support are dropped. Kernels at B_min
    subroutine new( self, box, smpd, pos_vox, support, amps )
        class(pcg_atomic_prior), intent(inout) :: self
        integer,                 intent(in)    :: box
        real,                    intent(in)    :: smpd
        real,                    intent(in)    :: pos_vox(:,:)
        real,                    intent(in)    :: support(:,:,:)
        real,          optional, intent(in)    :: amps(:)
        logical, allocatable :: keep(:)
        integer :: n, a, b, iv(3)
        real    :: minsep2
        call self%kill
        if( size(pos_vox,1) /= 3 ) THROW_HARD('atom positions are (3,natoms); new')
        if( any(shape(support) /= box) ) THROW_HARD('support dimensions differ from the box; new')
        if( smpd <= 0. ) THROW_HARD('sampling distance must be positive; new')
        n = size(pos_vox,2)
        if( present(amps) )then
            if( size(amps) /= n ) THROW_HARD('one column scale per atom; new')
            if( any(amps <= 0.) ) THROW_HARD('column scales must be positive; new')
        endif
        self%box  = box
        self%smpd = smpd
        self%support = support > 0.
        self%stats%natoms_in = n
        self%stats%nsupport  = count(self%support)
        allocate(keep(n), source=.true.)
        minsep2 = (MIN_SEP_A / smpd)**2
        do a = 1, n
            iv = nint(pos_vox(:,a))
            if( any(iv < 1) .or. any(iv > box) )then
                keep(a) = .false.
                self%stats%noutside = self%stats%noutside + 1
                cycle
            endif
            if( .not. self%support(iv(1),iv(2),iv(3)) )then
                keep(a) = .false.
                self%stats%noutside = self%stats%noutside + 1
            endif
        end do
        call drop_duplicates
        self%natoms = count(keep)
        if( self%natoms < 1 ) THROW_HARD('no atom of the model lies on the solve support; new')
        allocate(self%pos(3,self%natoms), self%amp(self%natoms))
        b = 0
        do a = 1, n
            if( .not. keep(a) ) cycle
            b = b + 1
            self%pos(:,b) = pos_vox(:,a)
            self%amp(b)   = 1.
            if( present(amps) ) self%amp(b) = amps(a)
        end do
        self%stats%natoms   = self%natoms
        self%stats%bfac_min = self%get_bfac_min()
        self%exists = .true.
        call self%build(self%stats%bfac_min)

      contains

        !> drops an atom closer than MIN_SEP_A to an earlier kept one, by a cell grid of the separation
        subroutine drop_duplicates
            integer, allocatable :: head(:,:,:), next(:)
            integer :: c(3), ncell, i, j, k, m
            real    :: cellsz
            cellsz = max(1., MIN_SEP_A / smpd)
            ncell  = ceiling(real(box) / cellsz) + 2
            allocate(head(0:ncell,0:ncell,0:ncell), source=0)
            allocate(next(n), source=0)
            do a = 1, n
                if( .not. keep(a) ) cycle
                c = min(ncell-1, max(1, int(pos_vox(:,a) / cellsz) + 1))
                do k = c(3)-1, c(3)+1
                    do j = c(2)-1, c(2)+1
                        do i = c(1)-1, c(1)+1
                            m = head(i,j,k)
                            do while( m > 0 )
                                if( sum((pos_vox(:,m) - pos_vox(:,a))**2) < minsep2 )then
                                    keep(a) = .false.
                                    exit
                                endif
                                m = next(m)
                            end do
                            if( .not. keep(a) ) exit
                        end do
                        if( .not. keep(a) ) exit
                    end do
                    if( .not. keep(a) ) exit
                end do
                if( .not. keep(a) )then
                    self%stats%ndup = self%stats%ndup + 1
                    cycle
                endif
                next(a) = head(c(1),c(2),c(3))
                head(c(1),c(2),c(3)) = a
            end do
            deallocate(head, next)
        end subroutine drop_duplicates

    end subroutine new

    !> Kernels of width B (A^2): exp(-r^2/(2 sigma^2)), B = 8 pi^2 sigma^2, scaled by the column scale and
    !! clipped to the support and the box; then the Gram matrix
    subroutine build( self, bfac )
        class(pcg_atomic_prior), intent(inout) :: self
        real,                    intent(in)    :: bfac
        integer :: a, di, dj, dk, i, j, k, l
        real    :: d2, two_sig2
        if( .not. self%exists ) THROW_HARD('new has not been called; build')
        if( bfac <= 0. ) THROW_HARD('kernel B must be positive; build')
        self%bfac  = bfac
        self%sigma = sqrt(bfac / (8. * PI * PI)) / self%smpd
        self%rad   = min(FOOT_RAD_MAX, max(FOOT_RAD_MIN, ceiling(FOOT_NSIGMA * self%sigma)))
        self%nside = 2 * self%rad + 1
        self%stats%foot_rad  = self%rad
        self%stats%sigma_vox = self%sigma
        self%stats%bfac_data = bfac - self%stats%bfac_min
        if( allocated(self%corner) ) deallocate(self%corner)
        if( allocated(self%fval)   ) deallocate(self%fval)
        allocate(self%corner(3,self%natoms))
        allocate(self%fval(self%nside**3,self%natoms), source=0.)
        two_sig2 = 2. * self%sigma**2
        !$omp parallel do default(shared) private(a,di,dj,dk,i,j,k,l,d2) schedule(static)
        do a = 1, self%natoms
            self%corner(:,a) = nint(self%pos(:,a)) - self%rad
            do dk = 0, self%nside-1
                k = self%corner(3,a) + dk
                if( k < 1 .or. k > self%box ) cycle
                do dj = 0, self%nside-1
                    j = self%corner(2,a) + dj
                    if( j < 1 .or. j > self%box ) cycle
                    do di = 0, self%nside-1
                        i = self%corner(1,a) + di
                        if( i < 1 .or. i > self%box ) cycle
                        if( .not. self%support(i,j,k) ) cycle
                        l  = 1 + di + self%nside * (dj + self%nside * dk)
                        d2 = (real(i) - self%pos(1,a))**2 + (real(j) - self%pos(2,a))**2 + (real(k) - self%pos(3,a))**2
                        self%fval(l,a) = self%amp(a) * exp(-d2 / two_sig2)
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        call self%build_gram
    end subroutine build

    !> G^T G in CSR from a cell grid of the footprint side: footprints of atoms more than one cell apart are disjoint
    subroutine build_gram( self )
        class(pcg_atomic_prior), intent(inout) :: self
        integer, allocatable :: head(:,:,:), next(:), cnt(:)
        integer :: a, b, c(3), ncell, i, j, k, nz
        ncell = self%box / self%nside + 3
        allocate(head(-1:ncell,-1:ncell,-1:ncell), source=0)
        allocate(next(self%natoms), source=0)
        ! reverse insertion keeps every cell list in ascending atom order
        do a = self%natoms, 1, -1
            c = cell_of(a)
            next(a) = head(c(1),c(2),c(3))
            head(c(1),c(2),c(3)) = a
        end do
        allocate(cnt(self%natoms), source=0)
        !$omp parallel do default(shared) private(a,b,c,i,j,k) schedule(static)
        do a = 1, self%natoms
            c = cell_of(a)
            do k = c(3)-1, c(3)+1
                do j = c(2)-1, c(2)+1
                    do i = c(1)-1, c(1)+1
                        b = head(i,j,k)
                        do while( b > 0 )
                            if( all(abs(self%corner(:,b) - self%corner(:,a)) < self%nside) ) cnt(a) = cnt(a) + 1
                            b = next(b)
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        if( allocated(self%gptr)  ) deallocate(self%gptr)
        if( allocated(self%gcol)  ) deallocate(self%gcol)
        if( allocated(self%gval)  ) deallocate(self%gval)
        if( allocated(self%gdiag) ) deallocate(self%gdiag)
        allocate(self%gptr(self%natoms+1))
        self%gptr(1) = 1
        do a = 1, self%natoms
            self%gptr(a+1) = self%gptr(a) + cnt(a)
        end do
        allocate(self%gcol(self%gptr(self%natoms+1)-1), self%gval(self%gptr(self%natoms+1)-1))
        allocate(self%gdiag(self%natoms), source=0._dp)
        !$omp parallel do default(shared) private(a,b,c,i,j,k,nz) schedule(static)
        do a = 1, self%natoms
            c  = cell_of(a)
            nz = self%gptr(a)
            do k = c(3)-1, c(3)+1
                do j = c(2)-1, c(2)+1
                    do i = c(1)-1, c(1)+1
                        b = head(i,j,k)
                        do while( b > 0 )
                            if( all(abs(self%corner(:,b) - self%corner(:,a)) < self%nside) )then
                                self%gcol(nz) = b
                                self%gval(nz) = self%overlap(a, b)
                                if( b == a ) self%gdiag(a) = self%gval(nz)
                                nz = nz + 1
                            endif
                            b = next(b)
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        if( any(self%gdiag <= 0._dp) ) THROW_HARD('an atom kernel is empty on the support; build_gram')
        self%stats%gram_nnz = size(self%gcol)
        deallocate(head, next, cnt)

      contains

        function cell_of( iatom ) result( cc )
            integer, intent(in) :: iatom
            integer :: cc(3)
            cc = min(ncell-1, max(0, (self%corner(:,iatom) + self%rad - 1) / self%nside))
        end function cell_of

    end subroutine build_gram

    !> <g_a, g_b> over the intersection of the two footprint cubes, in a fixed spatial order (exactly symmetric)
    real(dp) function overlap( self, a, b )
        class(pcg_atomic_prior), intent(in) :: self
        integer,                 intent(in) :: a, b
        integer :: off(3), lo(3), hi(3), di, dj, dk, la, lb
        off = self%corner(:,b) - self%corner(:,a)
        lo  = max(0, off)
        hi  = min(self%nside-1, self%nside-1+off)
        overlap = 0._dp
        do dk = lo(3), hi(3)
            do dj = lo(2), hi(2)
                do di = lo(1), hi(1)
                    la = 1 + di + self%nside * (dj + self%nside * dk)
                    lb = 1 + (di-off(1)) + self%nside * ((dj-off(2)) + self%nside * (dk-off(3)))
                    overlap = overlap + real(self%fval(la,a),dp) * real(self%fval(lb,b),dp)
                end do
            end do
        end do
    end function overlap

    !> G^T x
    function gather( self, x ) result( y )
        class(pcg_atomic_prior), intent(in) :: self
        real,                    intent(in) :: x(:,:,:)
        real(dp), allocatable :: y(:)
        real(dp) :: s
        integer  :: a, di, dj, dk, i, j, k
        allocate(y(self%natoms))
        !$omp parallel do default(shared) private(a,di,dj,dk,i,j,k,s) schedule(static)
        do a = 1, self%natoms
            s = 0._dp
            do dk = 0, self%nside-1
                k = self%corner(3,a) + dk
                if( k < 1 .or. k > self%box ) cycle
                do dj = 0, self%nside-1
                    j = self%corner(2,a) + dj
                    if( j < 1 .or. j > self%box ) cycle
                    do di = 0, self%nside-1
                        i = self%corner(1,a) + di
                        if( i < 1 .or. i > self%box ) cycle
                        s = s + real(self%fval(1 + di + self%nside * (dj + self%nside * dk), a),dp) * real(x(i,j,k),dp)
                    end do
                end do
            end do
            y(a) = s
        end do
        !$omp end parallel do
    end function gather

    !> y = y + alpha G c; serial, so concurrent half solves share the prior and the sum order is fixed
    subroutine scatter_add( self, c, alpha, y )
        class(pcg_atomic_prior), intent(in)    :: self
        real(dp),                intent(in)    :: c(:)
        real,                    intent(in)    :: alpha
        real,                    intent(inout) :: y(:,:,:)
        real    :: ca
        integer :: a, di, dj, dk, i, j, k
        do a = 1, self%natoms
            ca = alpha * real(c(a))
            if( ca == 0. ) cycle
            do dk = 0, self%nside-1
                k = self%corner(3,a) + dk
                if( k < 1 .or. k > self%box ) cycle
                do dj = 0, self%nside-1
                    j = self%corner(2,a) + dj
                    if( j < 1 .or. j > self%box ) cycle
                    do di = 0, self%nside-1
                        i = self%corner(1,a) + di
                        if( i < 1 .or. i > self%box ) cycle
                        y(i,j,k) = y(i,j,k) + ca * self%fval(1 + di + self%nside * (dj + self%nside * dk), a)
                    end do
                end do
            end do
        end do
    end subroutine scatter_add

    !> G c on the box
    function scatter( self, c ) result( x )
        class(pcg_atomic_prior), intent(in) :: self
        real(dp),                intent(in) :: c(:)
        real, allocatable :: x(:,:,:)
        allocate(x(self%box,self%box,self%box), source=0.)
        call self%scatter_add(c, 1., x)
    end function scatter

    !> q = (G^T G + eps D) p, D the diagonal of G^T G
    subroutine gram_matvec( self, p, q )
        class(pcg_atomic_prior), intent(in)  :: self
        real(dp),                intent(in)  :: p(:)
        real(dp),                intent(out) :: q(:)
        integer  :: a, nz
        real(dp) :: s
        !$omp parallel do default(shared) private(a,nz,s) schedule(static)
        do a = 1, self%natoms
            s = GRAM_RIDGE_REL * self%gdiag(a) * p(a)
            do nz = self%gptr(a), self%gptr(a+1)-1
                s = s + self%gval(nz) * p(self%gcol(nz))
            end do
            q(a) = s
        end do
        !$omp end parallel do
    end subroutine gram_matvec

    !> c = (G^T G + eps D)^-1 y by Jacobi-preconditioned CG in double precision, serial dots (reproducible)
    subroutine solve_gram( self, y, c, niters, rel_resid )
        class(pcg_atomic_prior), intent(in)  :: self
        real(dp),                intent(in)  :: y(:)
        real(dp), allocatable,   intent(out) :: c(:)
        integer,  optional,      intent(out) :: niters
        real(dp), optional,      intent(out) :: rel_resid
        real(dp), allocatable :: r(:), z(:), p(:), q(:), dg(:)
        real(dp) :: rz, rz_new, alpha, ynorm, rnorm
        integer  :: it
        allocate(c(self%natoms), source=0._dp)
        if( present(niters)    ) niters    = 0
        if( present(rel_resid) ) rel_resid = 0._dp
        ynorm = sqrt(dot_product(y, y))
        if( ynorm <= 0._dp ) return
        allocate(r(self%natoms), source=y)
        allocate(q(self%natoms))
        dg = self%gdiag * (1._dp + GRAM_RIDGE_REL)
        z  = r / dg
        p  = z
        rz = dot_product(r, z)
        rnorm = ynorm
        do it = 1, GRAM_MAXITS
            call self%gram_matvec(p, q)
            alpha = rz / dot_product(p, q)
            c = c + alpha * p
            r = r - alpha * q
            rnorm = sqrt(dot_product(r, r))
            if( rnorm <= GRAM_TOL * ynorm ) exit
            z      = r / dg
            rz_new = dot_product(r, z)
            p      = z + (rz_new / rz) * p
            rz     = rz_new
        end do
        if( present(niters)    ) niters    = min(it, GRAM_MAXITS)
        if( present(rel_resid) ) rel_resid = rnorm / ynorm
    end subroutine solve_gram

    !> P_G x = G (G^T G)^-1 G^T x
    function project( self, x ) result( px )
        class(pcg_atomic_prior), intent(in) :: self
        real,                    intent(in) :: x(:,:,:)
        real, allocatable :: px(:,:,:)
        real(dp), allocatable :: c(:)
        call self%solve_gram(self%gather(x), c)
        px = self%scatter(c)
    end function project

    !> hp = hp + lambda (p - P_G p), the prior's term of the normal operator
    subroutine apply( self, lambda, p, hp )
        class(pcg_atomic_prior), intent(in)    :: self
        real,                    intent(in)    :: lambda
        real,                    intent(in)    :: p(:,:,:)
        real,                    intent(inout) :: hp(:,:,:)
        real(dp), allocatable :: c(:)
        integer :: i, j, k
        call self%solve_gram(self%gather(p), c)
        !$omp parallel do collapse(2) default(shared) private(i,j,k) schedule(static)
        do k = 1, self%box
            do j = 1, self%box
                do i = 1, self%box
                    hp(i,j,k) = hp(i,j,k) + lambda * p(i,j,k)
                end do
            end do
        end do
        !$omp end parallel do
        call self%scatter_add(c, -lambda, hp)
    end subroutine apply

    !> Amplitudes fitted to the data, c = argmin |A G c - b|^2: CG in atom space on G^T H G c = G^T b with the
    !! reconstructor's normal operator H (op) and gathered right-hand side gb = G^T b; Jacobi on diag(G^T G)
    subroutine amplitudes( self, op, gb, c, maxits, rtol, niters, rel_resid )
        class(pcg_atomic_prior), intent(in)    :: self
        class(pcg_operator),     intent(inout) :: op
        real(dp),                intent(in)    :: gb(:)
        real(dp), allocatable,   intent(out)   :: c(:)
        integer,                 intent(in)    :: maxits
        real,                    intent(in)    :: rtol
        integer,                 intent(out)   :: niters
        real,                    intent(out)   :: rel_resid
        real, allocatable, target :: vin(:), vout(:)
        real, pointer         :: v3(:,:,:), w3(:,:,:)
        real(dp), allocatable :: r(:), z(:), p(:), q(:), dg(:)
        real(dp) :: rz, rz_new, alpha, bnorm, rnorm, pq
        integer  :: it, n3
        if( size(gb) /= self%natoms ) THROW_HARD('one right-hand side entry per atom; amplitudes')
        n3 = self%box**3
        if( op%size() /= n3 ) THROW_HARD('operator size differs from the box; amplitudes')
        allocate(c(self%natoms), source=0._dp)
        niters    = 0
        rel_resid = 0.
        bnorm = sqrt(dot_product(gb, gb))
        if( bnorm <= 0._dp ) return
        allocate(vin(n3), vout(n3))
        v3(1:self%box,1:self%box,1:self%box) => vin
        w3(1:self%box,1:self%box,1:self%box) => vout
        dg = self%gdiag
        r  = gb
        z  = r / dg
        p  = z
        rz = dot_product(r, z)
        rnorm = bnorm
        allocate(q(self%natoms))
        do it = 1, maxits
            v3 = 0.
            call self%scatter_add(p, 1., v3)
            call op%apply(vin, vout)
            q  = self%gather(w3)
            pq = dot_product(p, q)
            if( pq <= 0._dp ) exit
            alpha = rz / pq
            c = c + alpha * p
            r = r - alpha * q
            rnorm  = sqrt(dot_product(r, r))
            niters = it
            if( rnorm <= real(rtol,dp) * bnorm ) exit
            z      = r / dg
            rz_new = dot_product(r, z)
            p      = z + (rz_new / rz) * p
            rz     = rz_new
        end do
        rel_resid = real(rnorm / bnorm)
        nullify(v3, w3)
        deallocate(vin, vout)
    end subroutine amplitudes

    !> B_data, the B the map x carries beyond the narrowest kernel: per-shell slope c_s = Re<m_s,x_s>/|m_s|^2 of x on the
    !! model m (unbiased under map noise), log c_s = log K - B s^2/4 fitted over the shells where m and x correlate at
    !! WIDTH_FIT_CORR or more; K is discarded, B_data clamped at zero; 0 with nshells < WIDTH_FIT_MINSHELLS. m is
    !! simulated at FIT_SIGMA_VOX with a FIT_NSIGMA footprint, so neither aliasing nor truncation bends its spectrum
    subroutine fit_width( self, x, bfac_data, nshells )
        class(pcg_atomic_prior), intent(inout) :: self
        real,                    intent(in)    :: x(:,:,:)
        real,                    intent(out)   :: bfac_data
        integer,                 intent(out)   :: nshells
        type(image) :: mimg, ximg
        complex,  allocatable :: mc(:,:,:), xc(:,:,:)
        real,     allocatable :: m(:,:,:)
        real(dp), allocatable :: smx(:), smm(:), sxx(:), s2arr(:), yarr(:)
        real(dp) :: w, s2, sx, sy, sxy, sxx_fit, ys, slope, nfit
        real     :: bfac_fit
        integer  :: cdim(3), logi(3), i, j, k, sh, nyq
        if( any(shape(x) /= self%box) ) THROW_HARD('map dimensions differ from the box; fit_width')
        nyq = self%box / 2
        call self%model_map(FIT_SIGMA_VOX, ceiling(FIT_NSIGMA * FIT_SIGMA_VOX), m)
        bfac_fit = 8. * PI * PI * (FIT_SIGMA_VOX * self%smpd)**2
        call mimg%new([self%box,self%box,self%box], self%smpd)
        call ximg%new([self%box,self%box,self%box], self%smpd)
        call mimg%set_rmat(m, .false.)
        call ximg%set_rmat(x, .false.)
        deallocate(m)
        call mimg%fft()
        call ximg%fft()
        mc   = mimg%get_cmat()
        xc   = ximg%get_cmat()
        cdim = mimg%get_array_shape()
        allocate(smx(0:nyq), smm(0:nyq), sxx(0:nyq), source=0._dp)
        do k = 1, cdim(3)
            do j = 1, cdim(2)
                do i = 1, cdim(1)
                    logi = mimg%comp_addr_logi(i,j,k)
                    sh   = nint(sqrt(real(sum(logi**2))))
                    if( sh > nyq ) cycle
                    ! the packed half holds h >= 0; h = 0 and h = Nyquist have no stored mate
                    w = 2._dp
                    if( logi(1) == 0 .or. logi(1) == nyq ) w = 1._dp
                    smx(sh) = smx(sh) + w * real(real(conjg(mc(i,j,k)) * xc(i,j,k)),dp)
                    smm(sh) = smm(sh) + w * real(abs(mc(i,j,k)),dp)**2
                    sxx(sh) = sxx(sh) + w * real(abs(xc(i,j,k)),dp)**2
                end do
            end do
        end do
        call mimg%kill
        call ximg%kill
        allocate(s2arr(nyq), yarr(nyq))
        nshells = 0
        do sh = 1, nyq - 1
            if( smm(sh) <= 0._dp .or. sxx(sh) <= 0._dp ) cycle
            if( smx(sh) / sqrt(smm(sh) * sxx(sh)) < real(WIDTH_FIT_CORR,dp) ) cycle
            nshells = nshells + 1
            s2 = (real(sh,dp) / (real(self%box,dp) * real(self%smpd,dp)))**2
            s2arr(nshells) = s2
            yarr(nshells)  = log(smx(sh) / smm(sh))
        end do
        bfac_data = 0.
        self%stats%nshells_fit = nshells
        if( nshells < WIDTH_FIT_MINSHELLS ) return
        nfit = real(nshells,dp)
        sx   = sum(s2arr(1:nshells))
        sy   = sum(yarr(1:nshells))
        sxy  = sum(s2arr(1:nshells) * yarr(1:nshells))
        sxx_fit = sum(s2arr(1:nshells)**2)
        ys    = nfit * sxx_fit - sx * sx
        if( ys <= 0._dp ) return
        slope = (nfit * sxy - sx * sy) / ys
        bfac_data = max(0., bfac_fit + real(-4._dp * slope) - self%get_bfac_min())
    end subroutine fit_width

    !> sum of unit Gaussians of width sigma_vox at the atoms, footprint half-width rad, on the support
    subroutine model_map( self, sigma_vox, rad, m )
        class(pcg_atomic_prior), intent(in)  :: self
        real,                    intent(in)  :: sigma_vox
        integer,                 intent(in)  :: rad
        real, allocatable,       intent(out) :: m(:,:,:)
        real    :: two_sig2, d2
        integer :: a, c(3), i, j, k
        allocate(m(self%box,self%box,self%box), source=0.)
        two_sig2 = 2. * sigma_vox**2
        do a = 1, self%natoms
            c = nint(self%pos(:,a))
            do k = max(1, c(3)-rad), min(self%box, c(3)+rad)
                do j = max(1, c(2)-rad), min(self%box, c(2)+rad)
                    do i = max(1, c(1)-rad), min(self%box, c(1)+rad)
                        if( .not. self%support(i,j,k) ) cycle
                        d2 = (real(i) - self%pos(1,a))**2 + (real(j) - self%pos(2,a))**2 + (real(k) - self%pos(3,a))**2
                        m(i,j,k) = m(i,j,k) + exp(-d2 / two_sig2)
                    end do
                end do
            end do
        end do
    end subroutine model_map

    !> Strength from moments of the base pair over the support: N = var((xe-xo)/2), S = max(mean(xb^2) - N, eps N),
    !! E' = max(mean((P_G xb)^2) - N natoms/|V|, 0), C^2 = min(E'/S, C2_MAX); returns (2N/S)/(1-C^2), the ratio for a
    !! half solve (half-map noise 2N); the merged ratio (N/S)/(1-C^2) is kept in the statistics
    real function estimate_lambda( self, xb, xe, xo ) result( lambda_rel )
        class(pcg_atomic_prior), intent(inout) :: self
        real,                    intent(in)    :: xb(:,:,:), xe(:,:,:), xo(:,:,:)
        real, allocatable :: pgx(:,:,:)
        real(dp) :: nv, dmean, d2, xb2, e2, n_noise, s_sig, e_exp, ep, c2
        integer  :: i, j, k
        if( any(shape(xb) /= self%box) .or. any(shape(xe) /= self%box) .or. any(shape(xo) /= self%box) ) &
            &THROW_HARD('map dimensions differ from the box; estimate_lambda')
        pgx   = self%project(xb)
        nv    = real(self%stats%nsupport,dp)
        dmean = 0._dp
        d2    = 0._dp
        xb2   = 0._dp
        e2    = 0._dp
        do k = 1, self%box
            do j = 1, self%box
                do i = 1, self%box
                    if( .not. self%support(i,j,k) ) cycle
                    dmean = dmean + 0.5_dp * (real(xe(i,j,k),dp) - real(xo(i,j,k),dp))
                    d2    = d2    + 0.25_dp * (real(xe(i,j,k),dp) - real(xo(i,j,k),dp))**2
                    xb2   = xb2   + real(xb(i,j,k),dp)**2
                    e2    = e2    + real(pgx(i,j,k),dp)**2
                end do
            end do
        end do
        dmean   = dmean / nv
        n_noise = d2 / nv - dmean**2
        s_sig   = max(xb2 / nv - n_noise, real(SIGNAL_FLOOR,dp) * n_noise)
        e_exp   = e2 / nv
        ep      = max(e_exp - n_noise * real(self%natoms,dp) / nv, 0._dp)
        c2      = min(ep / s_sig, real(C2_MAX,dp))
        self%stats%noise     = real(n_noise)
        self%stats%signal    = real(s_sig)
        self%stats%explained = real(e_exp)
        self%stats%c2        = real(c2)
        self%stats%lambda_rel_merged = real((n_noise / s_sig) / (1._dp - c2))
        self%stats%lambda_rel_est    = real((2._dp * n_noise / s_sig) / (1._dp - c2))
        lambda_rel = self%stats%lambda_rel_est
    end function estimate_lambda

    function get_stats( self ) result( stats )
        class(pcg_atomic_prior), intent(in) :: self
        type(pcg_atomic_prior_stats) :: stats
        stats = self%stats
    end function get_stats

    pure integer function get_natoms( self )
        class(pcg_atomic_prior), intent(in) :: self
        get_natoms = self%natoms
    end function get_natoms

    pure integer function get_box( self )
        class(pcg_atomic_prior), intent(in) :: self
        get_box = self%box
    end function get_box

    !> B of the current kernels (A^2)
    pure real function get_bfac( self )
        class(pcg_atomic_prior), intent(in) :: self
        get_bfac = self%bfac
    end function get_bfac

    !> B of the narrowest kernel, a Gaussian of SIGMA_MIN_VOX voxels
    pure real function get_bfac_min( self )
        class(pcg_atomic_prior), intent(in) :: self
        get_bfac_min = 8. * PI * PI * (SIGMA_MIN_VOX * self%smpd)**2
    end function get_bfac_min

    pure logical function is_built( self )
        class(pcg_atomic_prior), intent(in) :: self
        is_built = self%exists .and. allocated(self%gdiag)
    end function is_built

    subroutine kill( self )
        class(pcg_atomic_prior), intent(inout) :: self
        if( allocated(self%pos)     ) deallocate(self%pos)
        if( allocated(self%amp)     ) deallocate(self%amp)
        if( allocated(self%corner)  ) deallocate(self%corner)
        if( allocated(self%fval)    ) deallocate(self%fval)
        if( allocated(self%support) ) deallocate(self%support)
        if( allocated(self%gptr)    ) deallocate(self%gptr)
        if( allocated(self%gcol)    ) deallocate(self%gcol)
        if( allocated(self%gval)    ) deallocate(self%gval)
        if( allocated(self%gdiag)   ) deallocate(self%gdiag)
        self%box    = 0
        self%natoms = 0
        self%rad    = 0
        self%nside  = 0
        self%smpd   = 1.
        self%bfac   = 0.
        self%sigma  = 0.
        self%stats  = pcg_atomic_prior_stats()
        self%exists = .false.
    end subroutine kill

end module simple_pcg_atomic_prior
