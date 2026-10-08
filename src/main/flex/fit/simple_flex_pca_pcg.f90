!@descr: flex_pca coupled M-step adapter and operator on the shared PCG recurrence
!  Design: doc/implementation_notes/completed/flex_pca_envelope_support.md, 3.4. CG on the native-lattice basis u:
!  b = S^H y and T = S^H S from 2x-lattice KB deposits (scale OSMPL_PAD_FAC**3), Nyquist-ball band limit,
!  hard support P, floored per-voxel coupled divide as preconditioner. maxits<=0 ships solve_coupled_basis_exp;
!  otherwise CG warm-starts from the masked, LS-scaled gridding solution. put_back writes E*u.
module simple_flex_pca_pcg
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
!$ use omp_lib, only: omp_get_thread_num, omp_get_max_threads
use simple_core_module_api, only: dp, dtiny, find_ldim_nptcls, fplane_type, logfhandle, ori, osmpl_pad_fac, pi, &
    &simple_exception, sp, string, sym, tic, timer_int_kind, tiny, toc
use simple_image,                         only: image
use simple_reconstructor,                 only: reconstructor
use simple_pcg_solver,                    only: pcg_operator, pcg_solver_options, pcg_solve, &
    &flex_pcg_outcome_t => pcg_solver_outcome, FLEX_PCG_STOP_INDEFINITE => PCG_STOP_INDEFINITE, &
    &PCG_XTOL, PCG_RESID_REPLACE, PCG_RHO_FLOOR_FRAC
use simple_linalg,                        only: solve_real_spd_complex
use simple_gridding,                      only: kb_stencil_envelope_1d
use simple_math,                          only: ceil_div, floor_div
use simple_flex_reconstructor_latent_ops, only: solve_coupled_basis_exp, pair_index
use simple_pcg_lattice,                   only: pcg_lattice
implicit none

public :: flex_pcg_t, flex_pcg_outcome_t, flex_pcg_environment, test_flex_pcg_operator
private
#include "simple_local_flags.inc"

!> kernel scale relative to the gridding density: the native-lattice normal operator is (1/N^3) T under
!! the forward-normalised FFT, the 2x-lattice kernel path returns (1/(2N)^3) T, so the ratio is padf**3
real,     parameter :: FLEX_PCG_KSCALE             = real(OSMPL_PAD_FAC)**3
!> relative ridge of the per-voxel coupled solve (as COUPLED_MSTEP_RIDGE_REL in the gridding solve)
real(dp), parameter :: FLEX_PCG_RIDGE_REL          = 1.0d-8
!> default Tikhonov term relative to the low-band mean density (production PCG_LAMBDA)
real,     parameter :: FLEX_PCG_LAMBDA_REL_DEFAULT = 1.0e-3
!> resampled support values below this are outside the hard domain (the note's 4.2 item 4; no larger
!! than the solver's PCG_SUPPORT_DIV_MIN = 0.1)
real,     parameter :: FLEX_PCG_SUPPORT_FLOOR      = 0.05


!> A run-owned support/window environment or a fit-owned exact copy. Its optional envelope
!! image follows the caller's lifecycle rather than being cached at module scope.
type :: flex_pcg_environment
    private
    type(image) :: envelope
    logical :: l_active = .false., l_initialized = .false.
  contains
    procedure :: new            => flex_env_new
    procedure :: copy_from      => flex_env_copy_from
    procedure :: kill           => flex_env_kill
    procedure :: active         => flex_env_active
    procedure :: apply          => flex_window_apply
    procedure :: apply_rec      => flex_window_apply_rec
    procedure :: install_window => flex_pcg_install_window
end type flex_pcg_environment

!> the padded-lattice geometry, KB window, wrap table, deposition envelope and solve support come from
!! pcg_lattice
type, extends(pcg_lattice) :: flex_pcg_t
    private
    integer :: ncomp = 0, npairs = 0
    ! persistent apply_operator scratch: this runs once per CG iteration, so allocating it per
    ! call made the solve mmap- and page-fault-bound rather than arithmetic-bound
    complex, allocatable :: op_cq(:,:,:,:)  !< per-component transforms on the padded lattice
    ! ---- band-ball index lists (memory): the accumulators and the packed kernels live only on the
    ! lattice points the band can reach. nband is the plane band in native units (set_band); the
    ! expanded list holds every expanded-lattice point within R = pf*(nband+1/2) + sqrt(3)*iwinsz + 1
    ! of a wrap image of the origin, the physical list every physical point the fold reads into.
    ! Slot maps are 0 outside the lists. eorder/efold replay fold_accum's (m,k,h>=0) loop order, so
    ! the fold is bitwise the fold of the dense arrays.
    integer :: nband = 0, nexp = 0, npk = 0, nfold = 0
    integer, allocatable :: emap(:,:,:)    !< (expanded index triple) -> expanded slot
    integer, allocatable :: pmap(:,:,:)    !< (physical i,j,k) -> physical slot
    integer, allocatable :: eijk(:,:)      !< (3,nexp) expanded slot -> index triple
    integer, allocatable :: pijk(:,:)      !< (3,npk)  physical slot -> (i,j,k)
    integer, allocatable :: eorder(:), efold(:) !< (nfold) fold pairs: kpk(:,efold(t)) += kacc(:,eorder(t))
    logical :: l_lists = .false.
    real,    allocatable :: env(:,:,:)      !< native-lattice KB gather envelope E, unity at the centre
    real,    allocatable :: dep_c(:,:,:)    !< the 2x-lattice deposition envelope cropped to the native box (RHS deapodization)
    real,    allocatable :: khat(:,:,:,:)   !< (npairs, cdim) finalized pair kernels
    real,    allocatable :: ridge(:,:)      !< (ncomp, nshell) optional Fourier shell ridge in the operator
    real,    allocatable :: rhofl(:,:)      !< (ncomp, 0:Rnat) shell-relative preconditioner floors
    real,    allocatable :: lam(:)          !< (ncomp) Tikhonov term of the operator and preconditioner
    real    :: lam_rel = FLEX_PCG_LAMBDA_REL_DEFAULT !< lam = lam_rel x low-band mean of the diagonal density
    integer :: verbose = 0                  !< CG iteration log level (set_verbose; 0 = outcome lines only)
    logical :: l_kernel = .false., l_ridge = .false., l_floor = .false.
    type(image) :: wimg                     !< persistent boxpd^3 work image (keeps its plans)
    type(image) :: nimg                     !< persistent box^3 work image (band limit, ridge, precond)
    ! Per-thread transform pool: independent transforms run one per thread on single-threaded FFTW plans
    ! (the image%construct_thread_safe_tmp_imgs idiom); threaded plans stalled in os_sem_down. Width nthr.
    type(image), allocatable :: wpool(:)    !< per-thread boxpd^3 images
    type(image), allocatable :: npool(:)    !< per-thread box^3 images (band limit)
    real,    allocatable :: tp_work(:,:,:,:)!< (box,box,box,nthr_pool)
    real,    allocatable :: tp_emb(:,:,:,:) !< (boxpd,boxpd,boxpd,nthr_pool)
    complex, allocatable :: tp_acc(:,:,:,:) !< (cdim,nthr_pool)
    integer :: nthr_pool = 0
    complex, allocatable :: pc_cq(:,:,:,:) !< persistent apply_precond accumulator
    complex, allocatable :: op_hq(:,:,:,:) !< coupled-multiply output, one plane set per component
    integer, allocatable :: pidx(:,:)      !< pair_index(min(q,r),max(q,r)) lookup, built once
    logical :: exists = .false.
  contains
    procedure :: new
    procedure :: kill
    procedure :: is_ready
    procedure :: set_window_sphere => flex_set_window_sphere
    procedure :: set_band
    procedure :: get_npk
    procedure :: get_nexp
    procedure :: alloc_accum
    procedure :: alloc_packed
    procedure :: alloc_rhs_accum
    procedure :: alloc_rhs_packed
    procedure :: accumulate
    procedure :: accumulate_rhs
    procedure :: fold_accum
    procedure :: fold_rhs
    procedure :: finalize
    procedure :: finalize_rhs
    procedure :: set_ridge
    procedure :: clear_ridge
    procedure :: set_lambda_relative
    procedure :: set_verbose
    procedure :: solve
    procedure :: cg_core
    procedure :: apply_operator
    procedure :: apply_precond
    procedure :: get_npairs
    procedure :: bytes_accum
    procedure :: bytes_packed
    procedure :: bytes_rhs_accum
    procedure :: bytes_rhs_packed
    procedure, private :: win_wraps
    procedure, private :: build_lists
    procedure, private :: require_lists
    procedure, private :: dot_all
    procedure, private :: bandlimit_img
    procedure, private :: ensure_pool
    procedure, private :: prep_floor
end type flex_pcg_t

!> Rank-1 view of the coupled component operator for the shared PCG engine.
type, extends(pcg_operator) :: flex_pcg_adapter
    class(flex_pcg_t), pointer :: client => null()
    real, pointer, contiguous  :: density(:,:,:,:) => null()
    integer                    :: density_lb(3) = 0
  contains
    procedure :: size    => flex_vector_size
    procedure :: apply   => flex_vector_apply
    procedure :: precond => flex_vector_precond
    procedure :: dot     => flex_vector_dot
end type flex_pcg_adapter

interface
    !> White-box self-test of the operator, checks (A)-(F) (submodule simple_flex_pca_pcg_tester, built
    !! only with BUILD_TESTS=ON)
    module subroutine test_flex_pcg_operator( box, nsamples, l_pass, passes, sweep )
        integer,           intent(in)  :: box, nsamples
        logical,           intent(out) :: l_pass
        logical, optional, intent(out) :: passes(6)
        logical, optional, intent(in)  :: sweep
    end subroutine test_flex_pcg_operator
end interface

contains

    ! ---------------- lifecycle ----------------

    !> geometry of the solve at the covariance box: the native lattice holds the unknown, the padf-times
    !! padded lattice the kernels and the right-hand sides; tables as in reconstructor_pcg%new
    subroutine new( self, box, smpd, ncomp )
        class(flex_pcg_t), intent(inout) :: self
        integer,           intent(in)    :: box, ncomp
        real,              intent(in)    :: smpd
        real, allocatable :: env1d(:)
        integer :: i, j, k
        call self%kill
        call self%new_lattice(box, smpd)
        self%ncomp  = ncomp
        self%npairs = (ncomp*(ncomp+1))/2
        ! the gather envelope of the native lattice (its reciprocal is the tail's gridding correction)
        call kb_stencil_envelope_1d(self%kbwin, box, env1d)
        allocate(self%env(box,box,box))
        do k = 1, box
            do j = 1, box
                do i = 1, box
                    self%env(i,j,k) = env1d(i)*env1d(j)*env1d(k)
                end do
            end do
        end do
        ! the deposition envelope of the doubled-coordinate scatter, cropped to the native box
        self%dep_c = self%native_deposition_envelope()
        deallocate(env1d)
        call self%wimg%new([self%boxpd,self%boxpd,self%boxpd], smpd)
        call self%nimg%new([box,box,box], smpd)
        ! pair_index LUT: the coupled multiply evaluated pair_index(min(q,r),max(q,r)) once per
        ! (component, component, voxel), which at rank 10 and cdim 129x256x256 is 845 million
        ! out-of-line calls per operator application
        allocate(self%pidx(ncomp,ncomp))
        do j = 1, ncomp
            do i = 1, ncomp
                self%pidx(i,j) = pair_index(min(i,j), max(i,j))
            end do
        end do
        self%lam_rel = FLEX_PCG_LAMBDA_REL_DEFAULT
        self%exists = .true.
    end subroutine new


    !> Tikhonov term relative to the low-band mean of each component's diagonal density (set before
    !! the solve; 0 disables it)
    subroutine set_lambda_relative( self, lam_rel )
        class(flex_pcg_t), intent(inout) :: self
        real,              intent(in)    :: lam_rel
        self%lam_rel = max(0.0, lam_rel)
        self%l_floor = .false.   ! lam is derived with the floors
    end subroutine set_lambda_relative

    !> CG iteration log level of every solve on this operator (0: the outcome line only)
    subroutine set_verbose( self, verbose )
        class(flex_pcg_t), intent(inout) :: self
        integer,           intent(in)    :: verbose
        self%verbose = max(0, verbose)
    end subroutine set_verbose

    subroutine kill( self )
        class(flex_pcg_t), intent(inout) :: self
        integer :: i_pool
        call self%kill_lattice
        if( allocated(self%env)    ) deallocate(self%env)
        if( allocated(self%dep_c)  ) deallocate(self%dep_c)
        if( allocated(self%khat)   ) deallocate(self%khat)
        if( allocated(self%ridge)  ) deallocate(self%ridge)
        if( allocated(self%rhofl)  ) deallocate(self%rhofl)
        if( allocated(self%lam)    ) deallocate(self%lam)
        if( allocated(self%op_cq)  ) deallocate(self%op_cq)
        if( allocated(self%tp_work) ) deallocate(self%tp_work)
        if( allocated(self%tp_emb)  ) deallocate(self%tp_emb)
        if( allocated(self%tp_acc)  ) deallocate(self%tp_acc)
        if( allocated(self%pc_cq)   ) deallocate(self%pc_cq)
        if( allocated(self%op_hq)   ) deallocate(self%op_hq)
        if( allocated(self%pidx)    ) deallocate(self%pidx)
        if( allocated(self%emap)    ) deallocate(self%emap)
        if( allocated(self%pmap)    ) deallocate(self%pmap)
        if( allocated(self%eijk)    ) deallocate(self%eijk)
        if( allocated(self%pijk)    ) deallocate(self%pijk)
        if( allocated(self%eorder)  ) deallocate(self%eorder)
        if( allocated(self%efold)   ) deallocate(self%efold)
        self%nband = 0; self%nexp = 0; self%npk = 0; self%nfold = 0; self%l_lists = .false.
        if( allocated(self%wpool) )then
            do i_pool = 1, size(self%wpool)
                call self%wpool(i_pool)%kill
            end do
            deallocate(self%wpool)
        endif
        if( allocated(self%npool) )then
            do i_pool = 1, size(self%npool)
                call self%npool(i_pool)%kill
            end do
            deallocate(self%npool)
        endif
        self%nthr_pool = 0
        if( self%exists )then
            call self%wimg%kill
            call self%nimg%kill
        endif
        self%l_kernel = .false.; self%l_ridge = .false.; self%l_floor = .false.
        self%exists = .false.
    end subroutine kill

    pure logical function is_ready( self )
        class(flex_pcg_t), intent(in) :: self
        is_ready = self%exists .and. self%l_kernel
    end function is_ready

    pure integer function get_npairs( self )
        class(flex_pcg_t), intent(in) :: self
        get_npairs = self%npairs
    end function get_npairs

    !> bytes of one full-range kernel accumulator (all pairs)
    pure real(dp) function bytes_accum( self )
        class(flex_pcg_t), intent(in) :: self
        bytes_accum = 4.d0 * real(self%npairs,dp) * real(self%nexp,dp)
    end function bytes_accum

    !> bytes of one packed (transport) kernel set
    pure real(dp) function bytes_packed( self )
        class(flex_pcg_t), intent(in) :: self
        bytes_packed = 4.d0 * real(self%npairs,dp) * real(self%npk,dp)
    end function bytes_packed

    !> bytes of one full-range right-hand-side accumulator (all components, complex)
    pure real(dp) function bytes_rhs_accum( self )
        class(flex_pcg_t), intent(in) :: self
        bytes_rhs_accum = 8.d0 * real(self%ncomp,dp) * real(self%nexp,dp)
    end function bytes_rhs_accum

    !> bytes of one packed (transport) right-hand-side set
    pure real(dp) function bytes_rhs_packed( self )
        class(flex_pcg_t), intent(in) :: self
        bytes_rhs_packed = 8.d0 * real(self%ncomp,dp) * real(self%npk,dp)
    end function bytes_rhs_packed

    ! ---------------- support ----------------

    !> The lattice's spherical support (pcg_lattice set_window_sphere), logged
    subroutine flex_set_window_sphere( self, mskrad )
        class(flex_pcg_t), intent(inout) :: self
        real,              intent(in)    :: mskrad
        call self%pcg_lattice%set_window_sphere(mskrad)
        if( .not. self%l_mask ) return
        write(logfhandle,'(A,F7.2,A,F6.3,A,F6.3,A,F6.3)') '>>> FLEX_PCA PCG support: sphere radius ', mskrad, &
            &' px, window mean ', sum(self%window)/real(size(self%window)), ', hard support fraction ', &
            &sum(self%mask)/real(size(self%mask)), ', nominal sphere fraction ', &
            &(4.0/3.0)*PI*mskrad**3/real(self%box)**3
    end subroutine flex_set_window_sphere

    ! ---------------- accumulation ----------------

    !> full-range accumulator of the doubled-coordinate kernel scatter, one row per pair;
    !! 1-based on every axis (index = logical - lims3(:,1) + 1): the scatter sees it through an
    !! assumed-shape dummy and fold_accum through an allocatable one, and the two must agree
    subroutine alloc_accum( self, kacc )
        class(flex_pcg_t), intent(in)  :: self
        real, allocatable, intent(out) :: kacc(:,:)
        call self%require_lists('alloc_accum')
        allocate(kacc(self%npairs, self%nexp), source=0.0)
    end subroutine alloc_accum

    !> packed (h >= 0, real) kernel set: the transport and reduction form
    subroutine alloc_packed( self, kpk )
        class(flex_pcg_t), intent(in)  :: self
        real, allocatable, intent(out) :: kpk(:,:)
        call self%require_lists('alloc_packed')
        allocate(kpk(self%npairs, self%npk), source=0.0)
    end subroutine alloc_packed

    !> full-range accumulator of the doubled-coordinate right-hand-side scatter, one row per component
    subroutine alloc_rhs_accum( self, racc )
        class(flex_pcg_t),    intent(in)  :: self
        complex, allocatable, intent(out) :: racc(:,:)
        call self%require_lists('alloc_rhs_accum')
        allocate(racc(self%ncomp, self%nexp), source=cmplx(0.,0.))
    end subroutine alloc_rhs_accum

    !> packed (h >= 0, complex) right-hand-side set: the transport and reduction form
    subroutine alloc_rhs_packed( self, rpk )
        class(flex_pcg_t),    intent(in)  :: self
        complex, allocatable, intent(out) :: rpk(:,:)
        call self%require_lists('alloc_rhs_packed')
        allocate(rpk(self%ncomp, self%npk), source=cmplx(0.,0.))
    end subroutine alloc_rhs_packed

    !> The plane band in native units. Builds the band-ball index lists; every accumulate / packed
    !! allocation needs them. A band of Rnat (the default when nothing is known) covers the whole
    !! lattice reachable by the stencil.
    subroutine set_band( self, nband )
        class(flex_pcg_t), intent(inout) :: self
        integer,           intent(in)    :: nband
        integer :: nb
        nb = max(1, min(self%Rnat, nband))
        if( self%l_lists .and. nb == self%nband ) return
        self%nband = nb
        call self%build_lists
    end subroutine set_band

    subroutine require_lists( self, who )
        class(flex_pcg_t), intent(in) :: self
        character(len=*),  intent(in) :: who
        if( .not. self%l_lists ) THROW_HARD(who//': band lists not built; call set_band first')
    end subroutine require_lists

    pure integer function get_npk( self )
        class(flex_pcg_t), intent(in) :: self
        get_npk = self%npk
    end function get_npk

    pure integer function get_nexp( self )
        class(flex_pcg_t), intent(in) :: self
        get_nexp = self%nexp
    end function get_nexp

    !> Expanded list: every expanded-lattice index within R of a wrap image of the origin, R the
    !! farthest a stencil point can land for a plane point inside the band (|loc| < pf*(nband+1/2),
    !! nint rounding 1/2, stencil half-diagonal sqrt(3)*iwinsz, +1 margin). Physical list and fold
    !! pairs: fold_accum's (m, k, h >= 0) loop, reading kacc at (wrap(h), k, m) into the physical
    !! address of (h,k,m), replayed here once and stored in the same order.
    subroutine build_lists( self )
        class(flex_pcg_t), intent(inout) :: self
        real    :: R, R2
        integer :: n1, n2, n3, ih, ik, im, h, k, m, hh, hm, km, mm, ne, np, nf, phys(3), es, ps, pass
        if( allocated(self%emap)   ) deallocate(self%emap)
        if( allocated(self%pmap)   ) deallocate(self%pmap)
        if( allocated(self%eijk)   ) deallocate(self%eijk)
        if( allocated(self%pijk)   ) deallocate(self%pijk)
        if( allocated(self%eorder) ) deallocate(self%eorder)
        if( allocated(self%efold)  ) deallocate(self%efold)
        R  = real(self%padf) * (real(self%nband) + 0.5) + sqrt(3.0) * real(self%iwinsz) + 1.0
        R2 = R * R
        n1 = self%lims3(1,2) - self%lims3(1,1) + 1
        n2 = self%lims3(2,2) - self%lims3(2,1) + 1
        n3 = self%lims3(3,2) - self%lims3(3,1) + 1
        allocate(self%emap(n1,n2,n3), source=0)
        do pass = 1, 2
            ne = 0
            do im = 1, n3
                m  = im - 1 + self%lims3(3,1)
                mm = min(abs(m), abs(m - self%boxpd), abs(m + self%boxpd))
                do ik = 1, n2
                    k  = ik - 1 + self%lims3(2,1)
                    km = min(abs(k), abs(k - self%boxpd), abs(k + self%boxpd))
                    do ih = 1, n1
                        h  = ih - 1 + self%lims3(1,1)
                        hm = min(abs(h), abs(h - self%boxpd), abs(h + self%boxpd))
                        if( real(hm*hm + km*km + mm*mm) > R2 ) cycle
                        ne = ne + 1
                        if( pass == 2 )then
                            self%emap(ih,ik,im) = ne
                            self%eijk(:,ne)     = [ih,ik,im]
                        endif
                    end do
                end do
            end do
            if( pass == 1 ) allocate(self%eijk(3,max(1,ne)), source=0)
        end do
        self%nexp = ne
        allocate(self%pmap(self%cdim(1),self%cdim(2),self%cdim(3)), source=0)
        do pass = 1, 2
            np = 0; nf = 0
            do m = self%lims3(3,1), self%lims3(3,2)
                do k = self%lims3(2,1), self%lims3(2,2)
                    do h = 0, -self%lims3(1,1)
                        hh = self%wrap(h)
                        ih = hh - self%lims3(1,1) + 1
                        ik = k  - self%lims3(2,1) + 1
                        im = m  - self%lims3(3,1) + 1
                        es = self%emap(ih,ik,im)
                        if( es == 0 ) cycle
                        phys = self%wimg%comp_addr_phys(h,k,m)
                        ps   = self%pmap(phys(1),phys(2),phys(3))
                        if( ps == 0 )then
                            np = np + 1
                            if( pass == 2 )then
                                self%pmap(phys(1),phys(2),phys(3)) = np
                                self%pijk(:,np) = phys
                            endif
                            ps = np
                        endif
                        nf = nf + 1
                        if( pass == 2 )then
                            self%eorder(nf) = es
                            self%efold(nf)  = ps
                        endif
                    end do
                end do
            end do
            if( pass == 1 )then
                allocate(self%pijk(3,max(1,np)), self%eorder(max(1,nf)), self%efold(max(1,nf)), source=0)
            else
                exit
            endif
            ! pass 1 marked nothing in pmap (ps stays 0 for every new point), so the count is exact
            ! only if pmap is reset between passes
            self%pmap = 0
        end do
        self%npk = np; self%nfold = nf
        self%l_lists = .true.
        write(logfhandle,'(A,I0,A,F6.2,A,I0,A,F6.2,A,I0,A)') '>>> FLEX_PCA PCG band lists: nband=', self%nband, &
            &'  expanded ', 100.0*real(self%nexp)/real(n1)/real(n2)/real(n3), ' % (', self%nexp, &
            &' points), physical ', 100.0*real(self%npk)/real(product(self%cdim)), ' % (', self%npk, ' points)'
        call flush(logfhandle)
    end subroutine build_lists

    !> pair-weighted Gram kernels of one batch: |T_i|^2 M_i(q,r) scattered at padf*R_i*[h,k,0] through the
    !! KB window, the same native-subset plane samples and Friedel storage as the coupled insert
    !! (insert_planes_multi); the flex plane lattice is already padf-oversampled
    subroutine accumulate( self, kacc, se, orientations, fpls, density_scales, valid, nrecords )
        class(flex_pcg_t), intent(in)    :: self
        real,              intent(inout) :: kacc(:,:)
        class(sym),        intent(inout) :: se
        type(ori),         intent(inout) :: orientations(:)
        type(fplane_type), intent(in)    :: fpls(:)
        real(dp),          intent(in)    :: density_scales(:,:,:)
        logical,           intent(in)    :: valid(:)
        integer,           intent(in)    :: nrecords
        type(ori) :: o_sym
        real,    allocatable :: rotmats(:,:,:,:), dpack(:,:)
        integer, allocatable :: fpllims(:,:,:), nyq_disks(:)
        real    :: loc(3), w(self%wdim,self%wdim,self%wdim), rot(3,3), ctfsq_raw
        integer :: i, isym, nsym, l, h, k, q, r, pf, i0(3), nyq_eff, fpllims_pd(3,2)
        integer :: h_sq, k_max_h, k_lo, k_hi
        if( nrecords < 1 ) return
        if( size(kacc,1) /= self%npairs ) THROW_HARD('pair kernel accumulator has the wrong leading extent; accumulate')
        call self%require_lists('accumulate')
        if( size(kacc,2) /= self%nexp ) THROW_HARD('pair kernel accumulator is not on the band list; accumulate')
        nsym = se%get_nsym()
        pf   = self%padf
        allocate(rotmats(3,3,nsym,nrecords), source=0.)
        allocate(dpack(self%npairs,nrecords), source=0.)
        allocate(fpllims(3,2,nrecords), nyq_disks(nrecords), source=0)
        do i = 1, nrecords
            if( .not. valid(i) ) cycle
            if( .not. allocated(fpls(i)%ctfsq_plane) ) THROW_HARD('ctfsq plane does not exist; accumulate')
            rotmats(:,:,1,i) = orientations(i)%get_mat()
            do isym = 2, nsym
                call se%apply(orientations(i), isym, o_sym)
                rotmats(:,:,isym,i) = o_sym%get_mat()
            end do
            fpllims_pd     = fpls(i)%frlims
            fpllims(:,:,i) = fpllims_pd
            fpllims(1,1,i) = ceil_div (fpllims_pd(1,1), pf)
            fpllims(1,2,i) = floor_div(fpllims_pd(1,2), pf)
            fpllims(2,1,i) = ceil_div (fpllims_pd(2,1), pf)
            fpllims(2,2,i) = floor_div(fpllims_pd(2,2), pf)
            nyq_eff = self%Rnat
            if( fpls(i)%nyq > 0 ) nyq_eff = min(nyq_eff, max(1, fpls(i)%nyq/pf))
            nyq_disks(i) = nyq_eff * (nyq_eff + 1)
            do r = 1, self%ncomp
                do q = 1, r
                    dpack(pair_index(q,r),i) = real(density_scales(q,r,i))
                end do
            end do
        end do
        call o_sym%kill
        !$omp parallel default(shared) private(i,isym,l,h,k,h_sq,k_max_h,k_lo,k_hi,ctfsq_raw,loc,i0,w,rot) &
        !$omp proc_bind(close)
        do i = 1, nrecords
            if( .not. valid(i) ) cycle
            do isym = 1, nsym
                rot = rotmats(:,:,isym,i)
                do l = 0, self%stride-1
                    !$omp do schedule(static,1)
                    do h = fpllims(1,1,i)+l, fpllims(1,2,i), self%stride
                        h_sq = h*h
                        if( h_sq > nyq_disks(i) ) cycle
                        k_max_h = int(sqrt(real(nyq_disks(i)-h_sq)))
                        k_lo = max(fpllims(2,1,i),-k_max_h)
                        k_hi = min(fpllims(2,2,i), k_max_h)
                        do k = k_lo, k_hi
                            if( k <= 0 )then
                                ctfsq_raw = fpls(i)%ctfsq_plane(h,k)
                            else
                                ctfsq_raw = fpls(i)%ctfsq_plane(-h,-k)
                            endif
                            if( ctfsq_raw <= TINY ) cycle
                            loc = real(pf) * matmul(real([h,k,0]), rot)
                            i0  = nint(loc) - self%iwinsz
                            if( self%win_wraps(i0) ) cycle   ! rim: serial pass below
                            call self%kbwin%apod_mat_3d_fast(loc, self%iwinsz, self%wdim, w)
                            call scatter_pairs_nowrap(self, i0, w, dpack(:,i), ctfsq_raw, kacc)
                        end do
                    end do
                    !$omp end do
                end do
                ! wrapping rim, serialized: the colouring does not survive folding
                !$omp single
                do h = fpllims(1,1,i), fpllims(1,2,i)
                    h_sq = h*h
                    if( h_sq > nyq_disks(i) ) cycle
                    k_max_h = int(sqrt(real(nyq_disks(i)-h_sq)))
                    k_lo = max(fpllims(2,1,i),-k_max_h)
                    k_hi = min(fpllims(2,2,i), k_max_h)
                    do k = k_lo, k_hi
                        if( h_sq + k*k <= self%sq_rim ) cycle   ! provably cannot wrap
                        if( k <= 0 )then
                            ctfsq_raw = fpls(i)%ctfsq_plane(h,k)
                        else
                            ctfsq_raw = fpls(i)%ctfsq_plane(-h,-k)
                        endif
                        if( ctfsq_raw <= TINY ) cycle
                        loc = real(pf) * matmul(real([h,k,0]), rot)
                        i0  = nint(loc) - self%iwinsz
                        if( .not. self%win_wraps(i0) ) cycle
                        call self%kbwin%apod_mat_3d_fast(loc, self%iwinsz, self%wdim, w)
                        call scatter_pairs_wrap(self, i0, w, dpack(:,i), ctfsq_raw, kacc)
                    end do
                end do
                !$omp end single
            end do
        end do
        !$omp end parallel
        deallocate(rotmats, dpack, fpllims, nyq_disks)
    end subroutine accumulate

    !> right-hand sides of one batch on the 2x lattice: the data scales E[z_iq] times the transfer-
    !! conjugated plane samples (with the plane padding factor, as the coupled insert deposits them),
    !! scattered at padf*R_i*[h,k,0] through the KB window; the same samples and Friedel storage as
    !! accumulate, so the kernels and the right-hand sides are one sample set
    subroutine accumulate_rhs( self, racc, se, orientations, fpls, data_scales, valid, nrecords )
        class(flex_pcg_t), intent(in)    :: self
        complex,           intent(inout) :: racc(:,:)
        class(sym),        intent(inout) :: se
        type(ori),         intent(inout) :: orientations(:)
        type(fplane_type), intent(in)    :: fpls(:)
        real(dp),          intent(in)    :: data_scales(:,:)
        logical,           intent(in)    :: valid(:)
        integer,           intent(in)    :: nrecords
        type(ori) :: o_sym
        real,    allocatable :: rotmats(:,:,:,:), dsc(:,:)
        integer, allocatable :: fpllims(:,:,:), nyq_disks(:)
        complex :: cmplx_raw, vals(self%ncomp)
        real    :: loc(3), w(self%wdim,self%wdim,self%wdim), rot(3,3), pf2
        integer :: i, isym, nsym, l, h, k, pf, i0(3), nyq_eff, fpllims_pd(3,2)
        integer :: h_sq, k_max_h, k_lo, k_hi
        if( nrecords < 1 ) return
        if( size(racc,1) /= self%ncomp ) THROW_HARD('rhs accumulator has the wrong leading extent; accumulate_rhs')
        call self%require_lists('accumulate_rhs')
        if( size(racc,2) /= self%nexp ) THROW_HARD('rhs accumulator is not on the band list; accumulate_rhs')
        nsym = se%get_nsym()
        pf   = self%padf
        pf2  = real(pf*pf)
        allocate(rotmats(3,3,nsym,nrecords), source=0.)
        allocate(dsc(self%ncomp,nrecords), source=0.)
        allocate(fpllims(3,2,nrecords), nyq_disks(nrecords), source=0)
        do i = 1, nrecords
            if( .not. valid(i) ) cycle
            if( .not. allocated(fpls(i)%transfer_plane) ) THROW_HARD('transfer plane does not exist; accumulate_rhs')
            rotmats(:,:,1,i) = orientations(i)%get_mat()
            do isym = 2, nsym
                call se%apply(orientations(i), isym, o_sym)
                rotmats(:,:,isym,i) = o_sym%get_mat()
            end do
            fpllims_pd     = fpls(i)%frlims
            fpllims(:,:,i) = fpllims_pd
            fpllims(1,1,i) = ceil_div (fpllims_pd(1,1), pf)
            fpllims(1,2,i) = floor_div(fpllims_pd(1,2), pf)
            fpllims(2,1,i) = ceil_div (fpllims_pd(2,1), pf)
            fpllims(2,2,i) = floor_div(fpllims_pd(2,2), pf)
            nyq_eff = self%Rnat
            if( fpls(i)%nyq > 0 ) nyq_eff = min(nyq_eff, max(1, fpls(i)%nyq/pf))
            nyq_disks(i) = nyq_eff * (nyq_eff + 1)
            dsc(:,i) = real(data_scales(1:self%ncomp,i))
        end do
        call o_sym%kill
        !$omp parallel default(shared) private(i,isym,l,h,k,h_sq,k_max_h,k_lo,k_hi,cmplx_raw,vals,loc,i0,w,rot) &
        !$omp proc_bind(close)
        do i = 1, nrecords
            if( .not. valid(i) ) cycle
            do isym = 1, nsym
                rot = rotmats(:,:,isym,i)
                do l = 0, self%stride-1
                    !$omp do schedule(static,1)
                    do h = fpllims(1,1,i)+l, fpllims(1,2,i), self%stride
                        h_sq = h*h
                        if( h_sq > nyq_disks(i) ) cycle
                        k_max_h = int(sqrt(real(nyq_disks(i)-h_sq)))
                        k_lo = max(fpllims(2,1,i),-k_max_h)
                        k_hi = min(fpllims(2,2,i), k_max_h)
                        do k = k_lo, k_hi
                            if( k <= 0 )then
                                cmplx_raw = conjg(fpls(i)%transfer_plane(h,k))*fpls(i)%cmplx_plane(h,k)
                            else
                                cmplx_raw = conjg(conjg(fpls(i)%transfer_plane(-h,-k))*fpls(i)%cmplx_plane(-h,-k))
                            endif
                            if( abs(real(cmplx_raw))+abs(aimag(cmplx_raw)) <= TINY ) cycle
                            loc = real(pf) * matmul(real([h,k,0]), rot)
                            i0  = nint(loc) - self%iwinsz
                            if( self%win_wraps(i0) ) cycle   ! rim: serial pass below
                            vals = dsc(:,i) * (pf2 * cmplx_raw)
                            call self%kbwin%apod_mat_3d_fast(loc, self%iwinsz, self%wdim, w)
                            call scatter_rhs_nowrap(self, i0, w, vals, racc)
                        end do
                    end do
                    !$omp end do
                end do
                !$omp single
                do h = fpllims(1,1,i), fpllims(1,2,i)
                    h_sq = h*h
                    if( h_sq > nyq_disks(i) ) cycle
                    k_max_h = int(sqrt(real(nyq_disks(i)-h_sq)))
                    k_lo = max(fpllims(2,1,i),-k_max_h)
                    k_hi = min(fpllims(2,2,i), k_max_h)
                    do k = k_lo, k_hi
                        if( h_sq + k*k <= self%sq_rim ) cycle
                        if( k <= 0 )then
                            cmplx_raw = conjg(fpls(i)%transfer_plane(h,k))*fpls(i)%cmplx_plane(h,k)
                        else
                            cmplx_raw = conjg(conjg(fpls(i)%transfer_plane(-h,-k))*fpls(i)%cmplx_plane(-h,-k))
                        endif
                        if( abs(real(cmplx_raw))+abs(aimag(cmplx_raw)) <= TINY ) cycle
                        loc = real(pf) * matmul(real([h,k,0]), rot)
                        i0  = nint(loc) - self%iwinsz
                        if( .not. self%win_wraps(i0) ) cycle
                        vals = dsc(:,i) * (pf2 * cmplx_raw)
                        call self%kbwin%apod_mat_3d_fast(loc, self%iwinsz, self%wdim, w)
                        call scatter_rhs_wrap(self, i0, w, vals, racc)
                    end do
                end do
                !$omp end single
            end do
        end do
        !$omp end parallel
        deallocate(rotmats, dsc, fpllims, nyq_disks)
    end subroutine accumulate_rhs

    pure logical function win_wraps( self, i0 )
        class(flex_pcg_t), intent(in) :: self
        integer,           intent(in) :: i0(3)
        win_wraps = any(i0 < self%wlims(1)) .or. any(i0 + self%wdim - 1 > self%wlims(2))
    end function win_wraps

    !> all pairs of one sample through one non-wrapping KB window; the pair index is the leading
    !! dimension so the innermost update is a contiguous vector operation
    !> stencil deposit into the band-list accumulator (nowrap: the stencil lies inside wlims)
    subroutine scatter_pairs_nowrap( self, i0, w, dpack, val, kacc )
        class(flex_pcg_t), intent(in)    :: self
        integer,           intent(in)    :: i0(3)
        real,              intent(in)    :: w(self%wdim,self%wdim,self%wdim), dpack(:), val
        real,              intent(inout) :: kacc(:,:)
        integer :: di, dj, dk, ih, ik, im, es
        real    :: wv
        do dk = 1, self%wdim
            im = i0(3) + dk - 1 - self%lims3(3,1) + 1
            do dj = 1, self%wdim
                ik = i0(2) + dj - 1 - self%lims3(2,1) + 1
                do di = 1, self%wdim
                    ih = i0(1) + di - 1 - self%lims3(1,1) + 1
                    es = self%emap(ih,ik,im)
                    if( es == 0 ) THROW_HARD('stencil point outside the band list (plane beyond nband); scatter_pairs_nowrap')
                    wv = w(di,dj,dk) * val
                    kacc(:,es) = kacc(:,es) + wv * dpack(:)
                end do
            end do
        end do
    end subroutine scatter_pairs_nowrap

    subroutine scatter_pairs_wrap( self, i0, w, dpack, val, kacc )
        class(flex_pcg_t), intent(in)    :: self
        integer,           intent(in)    :: i0(3)
        real,              intent(in)    :: w(self%wdim,self%wdim,self%wdim), dpack(:), val
        real,              intent(inout) :: kacc(:,:)
        integer :: di, dj, dk, ih, ik, im, es
        real    :: wv
        do dk = 1, self%wdim
            im = self%wrap(i0(3)+dk-1) - self%lims3(3,1) + 1
            do dj = 1, self%wdim
                ik = self%wrap(i0(2)+dj-1) - self%lims3(2,1) + 1
                do di = 1, self%wdim
                    ih = self%wrap(i0(1)+di-1) - self%lims3(1,1) + 1
                    es = self%emap(ih,ik,im)
                    if( es == 0 ) THROW_HARD('stencil point outside the band list (plane beyond nband); scatter_pairs_wrap')
                    wv = w(di,dj,dk) * val
                    kacc(:,es) = kacc(:,es) + wv * dpack(:)
                end do
            end do
        end do
    end subroutine scatter_pairs_wrap

    subroutine scatter_rhs_nowrap( self, i0, w, vals, racc )
        class(flex_pcg_t), intent(in)    :: self
        integer,           intent(in)    :: i0(3)
        real,              intent(in)    :: w(self%wdim,self%wdim,self%wdim)
        complex,           intent(in)    :: vals(:)
        complex,           intent(inout) :: racc(:,:)
        integer :: di, dj, dk, ih, ik, im, es
        do dk = 1, self%wdim
            im = i0(3) + dk - 1 - self%lims3(3,1) + 1
            do dj = 1, self%wdim
                ik = i0(2) + dj - 1 - self%lims3(2,1) + 1
                do di = 1, self%wdim
                    ih = i0(1) + di - 1 - self%lims3(1,1) + 1
                    es = self%emap(ih,ik,im)
                    if( es == 0 ) THROW_HARD('stencil point outside the band list (plane beyond nband); scatter_rhs_nowrap')
                    racc(:,es) = racc(:,es) + w(di,dj,dk) * vals(:)
                end do
            end do
        end do
    end subroutine scatter_rhs_nowrap

    subroutine scatter_rhs_wrap( self, i0, w, vals, racc )
        class(flex_pcg_t), intent(in)    :: self
        integer,           intent(in)    :: i0(3)
        real,              intent(in)    :: w(self%wdim,self%wdim,self%wdim)
        complex,           intent(in)    :: vals(:)
        complex,           intent(inout) :: racc(:,:)
        integer :: di, dj, dk, ih, ik, im, es
        do dk = 1, self%wdim
            im = self%wrap(i0(3)+dk-1) - self%lims3(3,1) + 1
            do dj = 1, self%wdim
                ik = self%wrap(i0(2)+dj-1) - self%lims3(2,1) + 1
                do di = 1, self%wdim
                    ih = self%wrap(i0(1)+di-1) - self%lims3(1,1) + 1
                    es = self%emap(ih,ik,im)
                    if( es == 0 ) THROW_HARD('stencil point outside the band list (plane beyond nband); scatter_rhs_wrap')
                    racc(:,es) = racc(:,es) + w(di,dj,dk) * vals(:)
                end do
            end do
        end do
    end subroutine scatter_rhs_wrap

    !> fold a full-range accumulator into the packed (h >= 0) real set, adding (the accumulator is freed);
    !! h runs to the 2x lattice's Nyquist plane (its accumulator column is the wrapped -Nyquist one)
    !> fold the band-list accumulator into the packed kernels: the same (m, k, h >= 0) visits as the
    !! dense fold, replayed from eorder/efold, so the sums are bitwise those of the dense arrays
    subroutine fold_accum( self, kacc, kpk )
        class(flex_pcg_t), intent(inout) :: self
        real, allocatable, intent(inout) :: kacc(:,:)
        real,              intent(inout) :: kpk(:,:)
        integer :: t
        if( .not. allocated(kacc) ) return
        if( size(kpk,1) /= self%npairs ) THROW_HARD('packed kernel set has the wrong leading extent; fold_accum')
        if( size(kacc,2) /= self%nexp .or. size(kpk,2) /= self%npk ) THROW_HARD('arrays are not on the band lists; fold_accum')
        do t = 1, self%nfold
            kpk(:,self%efold(t)) = kpk(:,self%efold(t)) + kacc(:,self%eorder(t))
        end do
        deallocate(kacc)
    end subroutine fold_accum
    
    !> the same fold for the complex right-hand-side accumulator
    subroutine fold_rhs( self, racc, rpk )
        class(flex_pcg_t),    intent(inout) :: self
        complex, allocatable, intent(inout) :: racc(:,:)
        complex,              intent(inout) :: rpk(:,:)
        integer :: t
        if( .not. allocated(racc) ) return
        if( size(rpk,1) /= self%ncomp ) THROW_HARD('packed rhs set has the wrong leading extent; fold_rhs')
        if( size(racc,2) /= self%nexp .or. size(rpk,2) /= self%npk ) THROW_HARD('arrays are not on the band lists; fold_rhs')
        do t = 1, self%nfold
            rpk(:,self%efold(t)) = rpk(:,self%efold(t)) + racc(:,self%eorder(t))
        end do
        deallocate(racc)
    end subroutine fold_rhs

    !> packed sums -> operator kernels: inverse transform on the 2x lattice, deposition envelope divided
    !! out, forward transform, real part, scaled to the gridding density (FLEX_PCG_KSCALE)
    subroutine finalize( self, kpk )
        class(flex_pcg_t), intent(inout) :: self
        real,              intent(in)    :: kpk(:,:)
        real,    pointer     :: rp(:,:,:)
        integer, parameter :: IBLK = 16
        real,    allocatable :: kt(:,:,:,:)   !< (cdim, npairs): each pair's kernel, written contiguously by its thread
        integer :: ipair, i, j, k, tid, t, ib, i2
        integer(timer_int_kind) :: t_fold
        if( size(kpk,1) /= self%npairs ) THROW_HARD('packed kernel set has the wrong leading extent; finalize')
        if( size(kpk,2) /= self%npk ) THROW_HARD('packed kernel set is not on the band list; finalize')
        t_fold = tic()
        ! One thread per pair writes a contiguous slab; a blocked transpose builds the pair-leading khat
        ! (per-voxel k x k multiply in solve). Not counted in solve()'s seconds.
        if( allocated(self%khat) )then
            if( size(self%khat,1) /= self%npairs .or. size(self%khat,2) /= self%cdim(1) ) deallocate(self%khat)
        endif
        if( .not. allocated(self%khat) ) allocate(self%khat(self%npairs, self%cdim(1), self%cdim(2), self%cdim(3)))
        allocate(kt(self%cdim(1), self%cdim(2), self%cdim(3), self%npairs))
        call self%ensure_pool()
        !$omp parallel do default(shared) private(ipair,i,j,k,t,tid,rp) schedule(static)
        do ipair = 1, self%npairs
            tid = 1
            !$ tid = omp_get_thread_num() + 1
            self%tp_acc(:,:,:,tid) = cmplx(0.,0.)
            do t = 1, self%npk
                self%tp_acc(self%pijk(1,t),self%pijk(2,t),self%pijk(3,t),tid) = cmplx(kpk(ipair,t), 0.)
            end do
            call self%wpool(tid)%set_cmat(self%tp_acc(:,:,:,tid))
            call self%wpool(tid)%ifft()
            ! divide the deposition envelope out IN PLACE on the image's own buffer: ifft leaves
            ! the image in real space, so the get_rmat / set_rmat round trip is pure overhead
            call self%wpool(tid)%get_rmat_ptr(rp)
            do k = 1, self%boxpd
                call self%divide_deposition_plane(rp, k)
            end do
            call self%wpool(tid)%fft()
            call self%wpool(tid)%get_cmat_sub(self%tp_acc(:,:,:,tid))
            do k = 1, self%cdim(3)
                do j = 1, self%cdim(2)
                    do i = 1, self%cdim(1)
                        kt(i,j,k,ipair) = real(self%tp_acc(i,j,k,tid)) * FLEX_PCG_KSCALE
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        ! blocked transpose (cdim, npairs) -> (npairs, cdim): per (j,k) row, IBLK consecutive i share
        ! one cache line of every pair's slab, and the npairs values of one voxel are written contiguously
        !$omp parallel do default(shared) private(k,j,ib,i2,ipair) schedule(static) collapse(2)
        do k = 1, self%cdim(3)
            do j = 1, self%cdim(2)
                do ib = 1, self%cdim(1), IBLK
                    do i2 = ib, min(ib+IBLK-1, self%cdim(1))
                        do ipair = 1, self%npairs
                            self%khat(ipair,i2,j,k) = kt(i2,j,k,ipair)
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        deallocate(kt)
        self%l_kernel = .true.
        write(logfhandle,'(A,I0,A,I0,A,F8.1)') '>>> FLEX_PCA PCG KERNEL FOLD: npairs=', &
            &self%npairs, '  threads=', self%nthr_pool, '  seconds=', real(toc(t_fold))
        call flush(logfhandle)
    end subroutine finalize

    !> packed right-hand sides -> b = P Pi S^H y on the native lattice: inverse transform on the 2x
    !! lattice, central crop, deposition envelope divided out, Nyquist-ball band limit, support
    subroutine finalize_rhs( self, rpk, b )
        class(flex_pcg_t), intent(inout) :: self
        complex,           intent(in)    :: rpk(:,:)
        real,              intent(out)   :: b(:,:,:,:)
        real, pointer :: rp(:,:,:)
        integer :: q, i, j, k, tid, off, t
        if( size(rpk,1) /= self%ncomp ) THROW_HARD('packed rhs set has the wrong leading extent; finalize_rhs')
        if( size(rpk,2) /= self%npk ) THROW_HARD('packed rhs set is not on the band list; finalize_rhs')
        call self%ensure_pool()
        off = (self%boxpd - self%box)/2
        ! parallel over components: independent transforms, per-thread images and scratch. The crop
        ! is inlined against the image buffer (offset (boxpd-box)/2, matching center_crop_real3d)
        ! and the band limit takes the thread's own native image.
        !$omp parallel do default(shared) private(q,i,j,k,t,tid,rp) schedule(static)
        do q = 1, self%ncomp
            tid = 1
            !$ tid = omp_get_thread_num() + 1
            self%tp_acc(:,:,:,tid) = cmplx(0.,0.)
            do t = 1, self%npk
                self%tp_acc(self%pijk(1,t),self%pijk(2,t),self%pijk(3,t),tid) = rpk(q,t)
            end do
            call self%wpool(tid)%set_cmat(self%tp_acc(:,:,:,tid))
            call self%wpool(tid)%ifft()
            call self%wpool(tid)%get_rmat_ptr(rp)
            self%tp_work(:,:,:,tid) = &
                &rp(off+1:off+self%box, off+1:off+self%box, off+1:off+self%box) / self%dep_c
            call self%bandlimit_img(self%tp_work(:,:,:,tid), self%npool(tid))
            if( self%l_mask ) self%tp_work(:,:,:,tid) = self%tp_work(:,:,:,tid) * self%mask
            b(:,:,:,q) = self%tp_work(:,:,:,tid)
        end do
        !$omp end parallel do
    end subroutine finalize_rhs

    subroutine set_ridge( self, ridge )
        class(flex_pcg_t), intent(inout) :: self
        real,              intent(in)    :: ridge(:,:)
        if( size(ridge,1) < self%ncomp ) THROW_HARD('ridge rank mismatch; set_ridge')
        if( allocated(self%ridge) ) deallocate(self%ridge)
        allocate(self%ridge(self%ncomp, size(ridge,2)), source=ridge(1:self%ncomp,:))
        self%l_ridge = .true.
    end subroutine set_ridge

    subroutine clear_ridge( self )
        class(flex_pcg_t), intent(inout) :: self
        if( allocated(self%ridge) ) deallocate(self%ridge)
        self%l_ridge = .false.
    end subroutine clear_ridge

    ! ---------------- operator, preconditioner, solve ----------------

    !> Nyquist-ball band limit on the native lattice (the shell rule of solve_coupled_basis_exp)
    !> allocate the per-thread transform pool once, OUTSIDE any parallel region (FFTW plan
    !! creation is not thread safe; image%new guards it with omp critical, but building the pool
    !! serially is both correct and cheaper). Every pooled image takes wthreads=.false. so that
    !! concurrent transforms never touch FFTW's thread pool.
    subroutine ensure_pool( self )
        class(flex_pcg_t), intent(inout) :: self
        integer :: nthr, t
        if( self%nthr_pool > 0 ) return
        nthr = 1
        !$ nthr = omp_get_max_threads()
        ! width = nthr, not min(nthr,ncomp): finalize folds npairs = ncomp*(ncomp+1)/2 kernels
        ! (55 at rank 10) and can use every thread, while apply_operator uses only the first ncomp
        nthr = max(1, nthr)
        allocate(self%wpool(nthr), self%npool(nthr))
        do t = 1, nthr
            call self%wpool(t)%new([self%boxpd,self%boxpd,self%boxpd], self%smpd, .false.)
            call self%npool(t)%new([self%box,self%box,self%box],       self%smpd, .false.)
        end do
        allocate(self%tp_work(self%box,self%box,self%box,nthr))
        allocate(self%tp_emb(self%boxpd,self%boxpd,self%boxpd,nthr))
        allocate(self%tp_acc(self%cdim(1),self%cdim(2),self%cdim(3),nthr))
        self%nthr_pool = nthr
    end subroutine ensure_pool

    !> band limit against a CALLER-SUPPLIED image, so it can run inside a parallel region over
    !! components. No internal OpenMP: the parallelism is the component loop outside. The shared-image
    !! variant below keeps its own threading for serial callers.
    subroutine bandlimit_img( self, work, img )
        class(flex_pcg_t), intent(inout) :: self
        real,              intent(inout) :: work(:,:,:)
        class(image),      intent(inout) :: img
        real, pointer :: rp(:,:,:)
        integer :: lims(3,2), h, k, m, phys(3)
        call img%set_rmat(work, .false.)
        call img%fft()
        lims = img%loop_lims(2)
        do m = lims(3,1), lims(3,2)
            do k = lims(2,1), lims(2,2)
                do h = lims(1,1), lims(1,2)
                    if( nint(sqrt(real(h*h + k*k + m*m))) > self%Rnat )then
                        phys = img%comp_addr_phys(h,k,m)
                        call img%set_cmat_at(phys(1),phys(2),phys(3), cmplx(0.,0.))
                    endif
                end do
            end do
        end do
        call img%ifft()
        call img%get_rmat_ptr(rp)
        work = rp(1:self%box,1:self%box,1:self%box)
    end subroutine bandlimit_img

    !> B x = P Pi T Pi P x (+ shell ridge): pair kernels on the 2x lattice, band-limited to the Nyquist
    !! ball, bracketed by the hard support
    subroutine apply_operator( self, x, hx )
        class(flex_pcg_t), intent(inout) :: self
        real,              intent(in)    :: x(:,:,:,:)
        real,              intent(out)   :: hx(:,:,:,:)
        real,    pointer :: rp(:,:,:)
        complex :: acc
        integer :: q, r, i, j, k, sh, h, kk, m, phys(3), nsh, off, tid
        if( .not. self%l_kernel ) THROW_HARD('kernels are not finalized; apply_operator')
        ! Persistent scratch: per-call allocation of these component-sized buffers left the solve page-fault bound.
        ! center_embed/center_crop are inlined on it with simple_cartesian_fourier's offset (boxpd-box)/2 and zero fill.
        if( .not. allocated(self%op_cq) ) &
            &allocate(self%op_cq(self%cdim(1),self%cdim(2),self%cdim(3),self%ncomp))
        off = (self%boxpd - self%box)/2
        call self%ensure_pool()
        ! Both component loops are parallel OVER COMPONENTS, each thread owning a padded image, a
        ! native image and its scratch. The transforms dominate: per CG iteration this is 2*ncomp
        ! on the boxpd lattice plus 4*ncomp on the native lattice (bandlimit does two each, twice
        ! per component), all of them independent. The previous code ran them serially through one
        ! shared self%wimg while threading only the k x k voxel loop, which left the dominant term
        ! single-threaded. Threading inside FFTW instead was measured to be a net harm.
        !$omp parallel do default(shared) private(q,tid) schedule(static)
        do q = 1, self%ncomp
            tid = 1
            !$ tid = omp_get_thread_num() + 1
            self%tp_work(:,:,:,tid) = x(:,:,:,q)
            if( self%l_mask ) self%tp_work(:,:,:,tid) = self%tp_work(:,:,:,tid) * self%mask
            call self%bandlimit_img(self%tp_work(:,:,:,tid), self%npool(tid))
            self%tp_emb(:,:,:,tid) = 0.0
            self%tp_emb(off+1:off+self%box, off+1:off+self%box, off+1:off+self%box, tid) = &
                &self%tp_work(:,:,:,tid)
            call self%wpool(tid)%set_rmat(self%tp_emb(:,:,:,tid), .false.)
            call self%wpool(tid)%fft()
            call self%wpool(tid)%get_cmat_sub(self%op_cq(:,:,:,q))
        end do
        !$omp end parallel do
        ! Coupled multiply, VOXEL-OUTER. Component-outer made each of the ncomp threads stream the
        ! whole khat array independently -- 1.86 GB at rank 10, so ~18.6 GB of reads per call -- and
        ! used only ncomp of nthr threads. Voxel-outer reads khat's contiguous npairs triangle once
        ! per voxel, produces every component's output from it, and collapses over the lattice so
        ! all threads work. Bitwise-exact: each output still sums r ascending.
        if( .not. allocated(self%op_hq) ) &
            &allocate(self%op_hq(self%cdim(1),self%cdim(2),self%cdim(3),self%ncomp))
        !$omp parallel do collapse(3) default(shared) private(i,j,k,q,r,acc) schedule(static)
        do k = 1, self%cdim(3)
            do j = 1, self%cdim(2)
                do i = 1, self%cdim(1)
                    do q = 1, self%ncomp
                        acc = cmplx(0.,0.)
                        do r = 1, self%ncomp
                            acc = acc + self%khat(self%pidx(q,r),i,j,k) * self%op_cq(i,j,k,r)
                        end do
                        self%op_hq(i,j,k,q) = acc
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        !$omp parallel do default(shared) private(q,tid,rp) schedule(static)
        do q = 1, self%ncomp
            tid = 1
            !$ tid = omp_get_thread_num() + 1
            call self%wpool(tid)%set_cmat(self%op_hq(:,:,:,q))
            call self%wpool(tid)%ifft()
            call self%wpool(tid)%get_rmat_ptr(rp)
            self%tp_work(:,:,:,tid) = &
                &rp(off+1:off+self%box, off+1:off+self%box, off+1:off+self%box)
            call self%bandlimit_img(self%tp_work(:,:,:,tid), self%npool(tid))
            if( self%l_mask ) self%tp_work(:,:,:,tid) = self%tp_work(:,:,:,tid) * self%mask
            hx(:,:,:,q) = self%tp_work(:,:,:,tid)
            ! Tikhonov term on the support (lam is derived with the floors; zero before prep_floor)
            if( allocated(self%lam) )then
                if( self%lam(q) > 0.0 )then
                    self%tp_work(:,:,:,tid) = x(:,:,:,q)
                    if( self%l_mask ) self%tp_work(:,:,:,tid) = self%tp_work(:,:,:,tid) * self%mask
                    hx(:,:,:,q) = hx(:,:,:,q) + self%lam(q) * self%tp_work(:,:,:,tid)
                endif
            endif
        end do
        !$omp end parallel do
        ! Fourier shell ridge on the native lattice (the cross-FSC precision, optional)
        if( self%l_ridge )then
            nsh = size(self%ridge,2)
            !$omp parallel do default(shared) private(q,tid,h,kk,m,sh,phys,rp) schedule(static)
            do q = 1, self%ncomp
                tid = 1
                !$ tid = omp_get_thread_num() + 1
                self%tp_work(:,:,:,tid) = x(:,:,:,q)
                if( self%l_mask ) self%tp_work(:,:,:,tid) = self%tp_work(:,:,:,tid) * self%mask
                call self%npool(tid)%set_rmat(self%tp_work(:,:,:,tid), .false.)
                call self%npool(tid)%fft()
                block
                    integer :: lims(3,2)
                    complex :: cval
                    lims = self%npool(tid)%loop_lims(2)
                    ! no nested OpenMP: the parallelism is the component loop outside
                    do m = lims(3,1), lims(3,2)
                        do kk = lims(2,1), lims(2,2)
                            do h = lims(1,1), lims(1,2)
                                sh   = nint(sqrt(real(h*h + kk*kk + m*m)))
                                phys = self%npool(tid)%comp_addr_phys(h,kk,m)
                                if( sh < 1 .or. sh > nsh .or. sh > self%Rnat )then
                                    cval = cmplx(0.,0.)
                                else
                                    cval = self%npool(tid)%get_cmat_at(phys(1),phys(2),phys(3)) &
                                        &* self%ridge(q,sh)
                                endif
                                call self%npool(tid)%set_cmat_at(phys(1),phys(2),phys(3), cval)
                            end do
                        end do
                    end do
                end block
                call self%npool(tid)%ifft()
                call self%npool(tid)%get_rmat_ptr(rp)
                self%tp_work(:,:,:,tid) = rp(1:self%box,1:self%box,1:self%box)
                if( self%l_mask ) self%tp_work(:,:,:,tid) = self%tp_work(:,:,:,tid) * self%mask
                hx(:,:,:,q) = hx(:,:,:,q) + self%tp_work(:,:,:,tid)
            end do
            !$omp end parallel do
        endif
    end subroutine apply_operator

    !> shell-relative floors of the preconditioner density: a fraction of the shell mean of each
    !! component's diagonal pair, so the per-voxel divide stays bounded on unsampled voxels
    subroutine prep_floor( self, rho, rho_lb )
        class(flex_pcg_t), intent(inout) :: self
        real,              intent(in)    :: rho(:,:,:,:)
        integer,           intent(in)    :: rho_lb(3)
        real(dp), allocatable :: ssum(:,:)
        integer,  allocatable :: scnt(:)
        integer :: lims(3,2), h, k, m, sh, q, ih, ik, im
        if( size(rho,1) /= self%npairs ) THROW_HARD('preconditioner density must be the full packed triangle; prep_floor')
        if( allocated(self%rhofl) ) deallocate(self%rhofl)
        allocate(self%rhofl(self%ncomp, 0:self%Rnat), source=0.0)
        allocate(ssum(self%ncomp, 0:self%Rnat), source=0.0_dp)
        allocate(scnt(0:self%Rnat), source=0)
        lims = self%nimg%loop_lims(2)
        do m = lims(3,1), lims(3,2)
            do k = lims(2,1), lims(2,2)
                do h = lims(1,1), lims(1,2)
                    sh = nint(sqrt(real(h*h + k*k + m*m)))
                    if( sh > self%Rnat ) cycle
                    ih = h - rho_lb(1) + 1; ik = k - rho_lb(2) + 1; im = m - rho_lb(3) + 1
                    scnt(sh) = scnt(sh) + 1
                    do q = 1, self%ncomp
                        ssum(q,sh) = ssum(q,sh) + real(max(0.0, rho(pair_index(q,q),ih,ik,im)),dp)
                    end do
                end do
            end do
        end do
        do sh = 0, self%Rnat
            if( scnt(sh) < 1 ) cycle
            self%rhofl(:,sh) = PCG_RHO_FLOOR_FRAC * real(ssum(:,sh) / real(scnt(sh),dp))
        end do
        ! Tikhonov term: lam_rel x the mean diagonal density over the low band (shells 1..max(4,Rnat/4)),
        ! the analogue of reconstructor_pcg's data scale
        if( allocated(self%lam) ) deallocate(self%lam)
        allocate(self%lam(self%ncomp), source=0.0)
        if( self%lam_rel > 0.0 )then
            block
                integer  :: shi
                real(dp) :: lsum(self%ncomp)
                integer  :: lcnt
                shi  = max(4, self%Rnat/4)
                lsum = 0.0_dp; lcnt = 0
                do sh = 1, min(shi, self%Rnat)
                    lsum = lsum + ssum(:,sh)
                    lcnt = lcnt + scnt(sh)
                end do
                if( lcnt > 0 ) self%lam = self%lam_rel * real(lsum / real(lcnt,dp))
            end block
        endif
        deallocate(ssum, scnt)
        self%l_floor = .true.
    end subroutine prep_floor

    !> z = P G P r: the per-voxel coupled divide on the gridding density (the Fourier-diagonal Jacobi
    !! inverse of T), diagonal pairs floored per shell, bracketed by the support. rho is the packed
    !! expanded-lattice density of the coupled insert with lower bounds rho_lb on the lattice axes.
    subroutine apply_precond( self, r, rho, rho_lb, z )
        class(flex_pcg_t), intent(inout) :: self
        real,              intent(in)    :: r(:,:,:,:)
        real,              intent(in)    :: rho(:,:,:,:)
        integer,           intent(in)    :: rho_lb(3)
        real,              intent(out)   :: z(:,:,:,:)
        real,    pointer     :: rp(:,:,:)
        complex(dp) :: rhs(self%ncomp), sol(self%ncomp)
        real(dp)    :: amat(self%ncomp,self%ncomp), diag_sum, ridge, denom
        integer :: lims(3,2), h, k, m, sh, q, rr, ih, ik, im, phys(3), flag, tid
        if( .not. self%l_floor ) call self%prep_floor(rho, rho_lb)
        call self%ensure_pool()
        if( .not. allocated(self%pc_cq) )then
            block
                integer :: ndim(3)
                ndim = self%nimg%get_array_shape()
                allocate(self%pc_cq(ndim(1),ndim(2),ndim(3),self%ncomp))
            end block
        endif
        ! forward transforms are independent per component: one per thread on the pool, and the
        ! accumulator is persistent (it was ~42 MB allocated and freed on every CG iteration)
        !$omp parallel do default(shared) private(q,tid) schedule(static)
        do q = 1, self%ncomp
            tid = 1
            !$ tid = omp_get_thread_num() + 1
            self%tp_work(:,:,:,tid) = r(:,:,:,q)
            if( self%l_mask ) self%tp_work(:,:,:,tid) = self%tp_work(:,:,:,tid) * self%mask
            call self%npool(tid)%set_rmat(self%tp_work(:,:,:,tid), .false.)
            call self%npool(tid)%fft()
            call self%npool(tid)%get_cmat_sub(self%pc_cq(:,:,:,q))
        end do
        !$omp end parallel do
        lims = self%nimg%loop_lims(2)
        !$omp parallel do collapse(3) default(shared) schedule(static) &
        !$omp private(h,k,m,sh,q,rr,ih,ik,im,phys,amat,rhs,sol,diag_sum,ridge,denom,flag)
        do m = lims(3,1), lims(3,2)
            do k = lims(2,1), lims(2,2)
                do h = lims(1,1), lims(1,2)
                    phys = self%nimg%comp_addr_phys(h,k,m)
                    sh   = nint(sqrt(real(h*h + k*k + m*m)))
                    if( sh > self%Rnat )then
                        self%pc_cq(phys(1),phys(2),phys(3),:) = cmplx(0.,0.)
                        cycle
                    endif
                    ih = h - rho_lb(1) + 1; ik = k - rho_lb(2) + 1; im = m - rho_lb(3) + 1
                    diag_sum = 0.d0
                    do q = 1, self%ncomp
                        do rr = q, self%ncomp
                            amat(q,rr) = real(rho(pair_index(q,rr),ih,ik,im), dp)
                            amat(rr,q) = amat(q,rr)
                        end do
                        amat(q,q)  = max(amat(q,q), real(self%rhofl(q,sh),dp)) + real(self%lam(q),dp)
                        diag_sum   = diag_sum + amat(q,q)
                        rhs(q)     = cmplx(self%pc_cq(phys(1),phys(2),phys(3),q), kind=dp)
                    end do
                    ridge = FLEX_PCG_RIDGE_REL * diag_sum / real(self%ncomp, dp)
                    do q = 1, self%ncomp
                        amat(q,q) = amat(q,q) + ridge
                    end do
                    call solve_real_spd_complex(amat, rhs, sol, self%ncomp, flag)
                    if( flag /= 0 )then
                        do q = 1, self%ncomp
                            denom = max(abs(amat(q,q)), ridge)
                            if( denom > DTINY )then
                                sol(q) = rhs(q) / denom
                            else
                                sol(q) = (0.d0, 0.d0)
                            endif
                        end do
                    endif
                    do q = 1, self%ncomp
                        self%pc_cq(phys(1),phys(2),phys(3),q) = cmplx(real(sol(q), sp), real(aimag(sol(q)), sp))
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        !$omp parallel do default(shared) private(q,tid,rp) schedule(static)
        do q = 1, self%ncomp
            tid = 1
            !$ tid = omp_get_thread_num() + 1
            call self%npool(tid)%set_cmat(self%pc_cq(:,:,:,q))
            call self%npool(tid)%ifft()
            call self%npool(tid)%get_rmat_ptr(rp)
            self%tp_work(:,:,:,tid) = rp(1:self%box,1:self%box,1:self%box)
            if( self%l_mask ) self%tp_work(:,:,:,tid) = self%tp_work(:,:,:,tid) * self%mask
            z(:,:,:,q) = self%tp_work(:,:,:,tid)
        end do
        !$omp end parallel do
    end subroutine apply_precond

    pure function dot_all( self, a, b ) result( d )
        class(flex_pcg_t), intent(in) :: self
        real,              intent(in) :: a(:,:,:,:), b(:,:,:,:)
        real(dp) :: d
        d = sum(real(a,dp) * real(b,dp))
    end function dot_all

    !> preconditioned CG on B x = b from the initial guess x (warm start), true-residual bookkeeping,
    !! periodic residual replacement, dx/x diminishing-returns stop, one cold restart on indefiniteness
    subroutine cg_core( self, b, x, rho, rho_lb, maxits, rtol, outcome, tag )
        class(flex_pcg_t), target, intent(inout) :: self
        real, contiguous, target,  intent(in)    :: b(:,:,:,:)
        real, contiguous, target,  intent(inout) :: x(:,:,:,:)
        real, contiguous, target,  intent(in)    :: rho(:,:,:,:)
        integer,                  intent(in)    :: rho_lb(3)
        integer,                  intent(in)    :: maxits
        real,                     intent(in)    :: rtol
        type(flex_pcg_outcome_t), intent(inout) :: outcome
        character(len=*),         intent(in)    :: tag
        real, pointer :: b_flat(:), x_flat(:)
        type(flex_pcg_adapter) :: op
        type(pcg_solver_options) :: options
        if( maxits > 0 .and. .not. self%l_kernel ) THROW_HARD('kernels are not finalized; cg_core')
        call self%prep_floor(rho, rho_lb)
        if( self%verbose > 0 ) call audit_terms
        op%client => self
        op%density => rho
        op%density_lb = rho_lb
        b_flat(1:size(b)) => b
        x_flat(1:size(x)) => x
        options%maxits = maxits
        options%rtol = rtol
        options%xtol = PCG_XTOL
        options%residual_replace = PCG_RESID_REPLACE
        options%rescale_start = .true.
        options%cold_restart = .true.
        options%track_preconditioned_residual = .false.
        options%record_history = .false.
        options%verbose = self%verbose
        options%tag = tag
        call pcg_solve(op, b_flat, x_flat, options, outcome)
        ! a vanishing curvature means the search direction has reached the operator's null space (the
        ! unregularised masked problem is only positive semi-definite): the iterate so far is the answer
        if( trim(outcome%stop_reason) == FLEX_PCG_STOP_INDEFINITE ) write(logfhandle,'(A,A,A,I0,A)') '>>> ', tag, &
            &': curvature vanished after ', outcome%iteration_count, &
            &' iterations (null space reached); returning the current iterate'
        if( .not. all(ieee_is_finite(x)) ) THROW_HARD(tag//': non-finite solution')

    contains

        !> How large are the two regularisers inside the operator, relative to the data term, on the
        !! right-hand side as a probe? lam is derived from the GRIDDING density while the data term is
        !! the doubled-lattice Toeplitz Gram, and the cross-FSC ridge arrives on the gridding scale too
        !! (add_invtausq2rho_coupled adds the same numbers to rho). This states the ratios instead of
        !! leaving the relative scaling to be inferred from the source.
        subroutine audit_terms
            real, allocatable :: h0(:,:,:,:), h1(:,:,:,:), lam_save(:)
            logical  :: ridge_save
            real(dp) :: n0, nl, nr
            allocate(h0, mold=b); allocate(h1, mold=b)
            ridge_save = self%l_ridge
            if( allocated(self%lam) ) lam_save = self%lam
            self%l_ridge = .false.
            if( allocated(self%lam) ) self%lam = 0.0
            call self%apply_operator(b, h0)                       ! data term alone
            n0 = sqrt(self%dot_all(h0,h0))
            if( allocated(lam_save) ) self%lam = lam_save
            call self%apply_operator(b, h1)                       ! data + Tikhonov
            h1 = h1 - h0
            nl = sqrt(self%dot_all(h1,h1))
            h1 = h1 + h0                                          ! restore data + Tikhonov
            self%l_ridge = ridge_save
            nr = 0.0_dp
            if( self%l_ridge )then
                call self%apply_operator(b, h0)                   ! data + Tikhonov + FSC prior
                h0 = h0 - h1
                nr = sqrt(self%dot_all(h0,h0))
            endif
            write(logfhandle,'(A,A,A,ES10.3,A,ES10.3,A,L1)') '>>> ', tag, &
                &' operator terms on the rhs probe: Tikhonov/data=', real(nl / max(n0,1.0e-30_dp)), &
                &'  FSC-prior/data=', real(nr / max(n0,1.0e-30_dp)), '  prior_installed=', self%l_ridge
            deallocate(h0, h1)
        end subroutine audit_terms

    end subroutine cg_core

    integer function flex_vector_size( self ) result(n)
        class(flex_pcg_adapter), intent(in) :: self
        n = self%client%box**3 * self%client%ncomp
    end function flex_vector_size

    subroutine flex_vector_apply( self, x, y )
        class(flex_pcg_adapter), intent(inout) :: self
        real, contiguous, target, intent(in)   :: x(:)
        real, contiguous, target, intent(out)  :: y(:)
        real, pointer :: x4(:,:,:,:), y4(:,:,:,:)
        x4(1:self%client%box,1:self%client%box,1:self%client%box,1:self%client%ncomp) => x
        y4(1:self%client%box,1:self%client%box,1:self%client%box,1:self%client%ncomp) => y
        call self%client%apply_operator(x4, y4)
    end subroutine flex_vector_apply

    subroutine flex_vector_precond( self, r, z )
        class(flex_pcg_adapter), intent(inout) :: self
        real, contiguous, target, intent(in)   :: r(:)
        real, contiguous, target, intent(out)  :: z(:)
        real, pointer :: r4(:,:,:,:), z4(:,:,:,:)
        r4(1:self%client%box,1:self%client%box,1:self%client%box,1:self%client%ncomp) => r
        z4(1:self%client%box,1:self%client%box,1:self%client%box,1:self%client%ncomp) => z
        call self%client%apply_precond(r4, self%density, self%density_lb, z4)
    end subroutine flex_vector_precond

    real(dp) function flex_vector_dot( self, a, b ) result(value)
        class(flex_pcg_adapter), intent(in) :: self
        real, contiguous, target, intent(in) :: a(:), b(:)
        real, pointer :: a4(:,:,:,:), b4(:,:,:,:)
        a4(1:self%client%box,1:self%client%box,1:self%client%box,1:self%client%ncomp) => a
        b4(1:self%client%box,1:self%client%box,1:self%client%box,1:self%client%ncomp) => b
        value = self%client%dot_all(a4,b4)
    end function flex_vector_dot

    !> M-step solve on one half. maxits<=0 ships the gridding solution (solve_coupled_basis_exp); otherwise
    !! CG on the 2x right-hand sides starts from it, masked to the support (cg_core rescales it).
    !! Returns E*u on the expanded lattice.
    subroutine solve( self, Y, rho, rpk, maxits, rtol, outcome, tag )
        class(flex_pcg_t),          intent(inout) :: self
        type(reconstructor),        intent(inout) :: Y(:)
        real,                       intent(in)    :: rho(:,:,:,:)
        complex,                    intent(in)    :: rpk(:,:)
        integer,                    intent(in)    :: maxits
        real,                       intent(in)    :: rtol
        type(flex_pcg_outcome_t),   intent(out)   :: outcome
        character(len=*), optional, intent(in)    :: tag
        real, allocatable :: b(:,:,:,:), x(:,:,:,:), xgrid(:,:,:,:)
        real, pointer     :: rp(:,:,:)
        integer  :: q, rho_lb(3), verbose
        real(dp) :: gnorm, dnorm, gdot
        integer(timer_int_kind) :: t0
        character(len=:), allocatable :: ttag
        ttag = 'FLEX_PCA PCG MSTEP'
        if( present(tag) ) ttag = trim(tag)
        if( size(Y) /= self%ncomp ) THROW_HARD('reconstructor count differs from the solver rank; solve')
        t0 = tic()
        outcome%requested_maxits = maxits
        rho_lb = lbound(Y(1)%cmat_exp)
        verbose = self%verbose
        allocate(x(self%box,self%box,self%box,self%ncomp))
        ! Y(q) transforms run one at a time below (ifft, put_back's fft), so use unthreaded FFTW plans:
        ! threaded plans driven serially stalled in os_sem_down. Set here to leave the gridding path alone.
        do q = 1, self%ncomp
            call Y(q)%set_wthreads(.false.)
        end do
        ! the gridding solution, physical units (the tail's inverse envelope undone here)
        call solve_coupled_basis_exp(Y, rho, self%ncomp)
        do q = 1, self%ncomp
            Y(q)%rho_exp = 1.0
            call Y(q)%compress_exp
            call Y(q)%ifft
            call Y(q)%get_rmat_ptr(rp)
            x(:,:,:,q) = rp(1:self%box,1:self%box,1:self%box) / self%env
        end do
        ! budget 0 ships the gridding solution untouched (the tail's soft mask is inside the hard support
        ! anyway, but the half-set FSC that sets the Wiener filter is taken on the unmasked halves and the
        ! support would inflate it); with iterations the solution lives on the support by construction
        if( maxits <= 0 )then
            outcome%stop_reason     = 'block_jacobi'
            outcome%iteration_count = 0
            outcome%start_norm      = real(sqrt(self%dot_all(x,x)))
            call put_back
            outcome%seconds = real(toc(t0))
            deallocate(x)
            return
        endif
        if( .not. self%l_kernel ) THROW_HARD('kernels are not finalized; solve')
        gnorm = 0._dp
        if( verbose > 0 )then
            allocate(xgrid, source=x)
            ! compare on the SUPPORT: the CG solution is zero outside it while the gridding reference is a
            ! full-box object, and with a 2%-occupancy envelope an unmasked comparison reports orthogonality
            ! no matter how close the two solutions are where they are both defined
            if( self%l_mask )then
                do q = 1, self%ncomp
                    xgrid(:,:,:,q) = xgrid(:,:,:,q) * self%mask
                end do
            endif
            gnorm = sqrt(self%dot_all(xgrid,xgrid))
        endif
        ! warm start: the full-box gridding solution masked onto the support
        if( self%l_mask )then
            do q = 1, self%ncomp
                x(:,:,:,q) = x(:,:,:,q) * self%mask
            end do
        endif
        allocate(b(self%box,self%box,self%box,self%ncomp))
        call self%finalize_rhs(rpk, b)
        call self%cg_core(b, x, rho, rho_lb, maxits, rtol, outcome, ttag)
        if( allocated(xgrid) )then
            gdot  = self%dot_all(x,xgrid)
            xgrid = x - xgrid
            dnorm = sqrt(self%dot_all(xgrid,xgrid))
            write(logfhandle,'(A,A,A,I0,A,ES10.3,A,F8.4)') '>>> ', ttag, ' warm-started solve vs the gridding &
                &reference: iters=', outcome%iteration_count, '  ||x-xg||/||xg||=', &
                &real(dnorm / max(gnorm, 1.0e-30_dp)), '  corr(x,xg)=', &
                &real(gdot / max(sqrt(self%dot_all(x,x))*gnorm, 1.0e-30_dp))
            deallocate(xgrid)
        endif
        call put_back
        outcome%seconds = real(toc(t0))
        deallocate(b, x)

    contains

        !> E*u back into the reconstructors' expanded lattice: the tail divides by E and ships u
        subroutine put_back
            integer :: qq
            do qq = 1, self%ncomp
                call Y(qq)%set_rmat(x(:,:,:,qq) * self%env, .false.)
                call Y(qq)%fft
                call Y(qq)%expand_exp
            end do
        end subroutine put_back

    end subroutine solve

    ! ---------------- envelope support (pcg_mskfile) ----------------

    !> pcg_mskfile resampled to a target lattice under the note's 4.2 contract: Fourier crop, ringing
    !! diagnostics (min/max/out-of-range fraction before clamping), clamp to [0,1], floor to exactly 0
    !! below FLEX_PCG_SUPPORT_FLOOR, occupancy of the box and of the mskdiam sphere on the native and
    !! the target lattice
    subroutine flex_pcg_support_volume( pcg_mskfile, box, smpd, mskdiam, box_t, smpd_t, img, tag )
        type(string),     intent(in)    :: pcg_mskfile
        integer,          intent(in)    :: box, box_t
        real,             intent(in)    :: smpd, mskdiam, smpd_t
        type(image),      intent(inout) :: img
        character(len=*), intent(in)    :: tag
        type(image)   :: nat, rc
        real, pointer :: rp(:,:,:)
        real    :: vmin, vmax, frac_out, occ_box_n, occ_sph_n, occ_box_t, occ_sph_t
        integer :: ldim(3), ifoo, n
        if( .not. flex_mskfile_set(pcg_mskfile) ) THROW_HARD('pcg_mskfile is not set; flex_pcg_support_volume')
        call find_ldim_nptcls(pcg_mskfile, ldim, ifoo)
        if( ldim(1) /= ldim(2) .or. ldim(1) /= ldim(3) ) THROW_HARD('pcg_mskfile must be a cubic volume')
        if( ldim(1) /= box ) THROW_HARD('pcg_mskfile box differs from the project box')
        ! native occupancy (the reference for the resampling diagnostic)
        call nat%new(ldim, smpd)
        call nat%read(pcg_mskfile)
        call nat%get_rmat_ptr(rp)
        call occupancy(rp, ldim(1), 0.5*mskdiam/smpd, occ_box_n, occ_sph_n)
        call nat%kill
        ! Fourier crop to the target lattice; the result is copied into a fresh image of the target box so
        ! that its buffers are exactly the target's (the in-place clip keeps the native allocation)
        call rc%read_and_crop(pcg_mskfile, smpd, box_t, smpd_t)
        if( rc%is_ft() ) call rc%ifft
        call rc%get_rmat_ptr(rp)
        n        = box_t**3
        vmin     = minval(rp(1:box_t,1:box_t,1:box_t))
        vmax     = maxval(rp(1:box_t,1:box_t,1:box_t))
        frac_out = real(count(rp(1:box_t,1:box_t,1:box_t) < 0.0 .or. rp(1:box_t,1:box_t,1:box_t) > 1.0)) / real(n)
        rp(1:box_t,1:box_t,1:box_t) = min(1.0, max(0.0, rp(1:box_t,1:box_t,1:box_t)))
        where( rp(1:box_t,1:box_t,1:box_t) < FLEX_PCG_SUPPORT_FLOOR ) rp(1:box_t,1:box_t,1:box_t) = 0.0
        if( .not. any(rp(1:box_t,1:box_t,1:box_t) > 0.0) ) THROW_HARD('pcg_mskfile is empty after resampling')
        call occupancy(rp, box_t, 0.5*mskdiam/smpd_t, occ_box_t, occ_sph_t)
        call img%new([box_t,box_t,box_t], smpd_t)
        call img%set_rmat(rp(1:box_t,1:box_t,1:box_t), .false.)
        call rc%kill
        write(logfhandle,'(A,A,A,I0,A,F7.3,A,F6.3,A,F6.3,A,F7.4)') '>>> FLEX_PCA PCG SUPPORT (', trim(tag), &
            &'): pcg_mskfile at box ', box_t, ' smpd ', smpd_t, ': before clamp min=', vmin, ' max=', vmax, &
            &' out-of-[0,1] fraction=', frac_out
        write(logfhandle,'(A,F6.3,A,F6.3,A,F6.3,A,F6.3)') '>>> FLEX_PCA PCG SUPPORT occupancy (window>0): native box ', &
            &occ_box_n, ' / sphere ', occ_sph_n, '  ->  target box ', occ_box_t, ' / sphere ', occ_sph_t
        call flush(logfhandle)

    contains

        subroutine occupancy( v, nb, rad, f_box, f_sph )
            real,    intent(in)  :: v(:,:,:)
            integer, intent(in)  :: nb
            real,    intent(in)  :: rad
            real,    intent(out) :: f_box, f_sph
            integer :: i, j, k, c, nin, nsph
            real    :: r2, rad2
            c = nb/2 + 1; rad2 = rad*rad
            nin = 0; nsph = 0
            do k = 1, nb
                do j = 1, nb
                    do i = 1, nb
                        r2 = real((i-c)**2 + (j-c)**2 + (k-c)**2)
                        if( r2 <= rad2 )then
                            nsph = nsph + 1
                            if( v(i,j,k) > 0.0 ) nin = nin + 1
                        endif
                    end do
                end do
            end do
            f_box = real(count(v(1:nb,1:nb,1:nb) > 0.0)) / real(nb)**3
            f_sph = 0.0
            if( nsph > 0 ) f_sph = real(nin) / real(nsph)
        end subroutine occupancy

    end subroutine flex_pcg_support_volume

    !> whether pcg_mskfile was given (the string is unallocated when absent)
    logical function flex_mskfile_set( pcg_mskfile )
        type(string), intent(in) :: pcg_mskfile
        flex_mskfile_set = .false.
        if( pcg_mskfile%is_allocated() ) flex_mskfile_set = len_trim(pcg_mskfile%to_char()) > 0
    end function flex_mskfile_set

    !> idempotent: the envelope at the covariance box when pcg_mskfile is set on the PCG backend
    subroutine flex_env_new( self, rec_backend, pcg_mskfile, box, smpd, mskdiam, box_crop, smpd_crop )
        class(flex_pcg_environment), intent(inout) :: self
        character(len=*), intent(in) :: rec_backend
        type(string),     intent(in) :: pcg_mskfile
        integer,          intent(in) :: box, box_crop
        real,             intent(in) :: smpd, mskdiam, smpd_crop
        call self%kill
        self%l_initialized = .true.
        self%l_active = trim(rec_backend) == 'pcg' .and. flex_mskfile_set(pcg_mskfile)
        if( .not. self%l_active ) return
        call flex_pcg_support_volume(pcg_mskfile, box, smpd, mskdiam, box_crop, smpd_crop, self%envelope, &
            &'covariance box')
    end subroutine flex_env_new

    !> Clone an already-resampled support environment without reading the mask again.
    subroutine flex_env_copy_from( self, source )
        class(flex_pcg_environment), intent(inout) :: self
        class(flex_pcg_environment), intent(in)    :: source
        if( .not. source%l_initialized ) THROW_HARD('flex_env_copy_from: source environment is not initialized')
        call self%kill
        self%l_initialized = .true.
        self%l_active      = source%l_active
        if( self%l_active ) call self%envelope%copy(source%envelope)
    end subroutine flex_env_copy_from

    subroutine flex_env_kill( self )
        class(flex_pcg_environment), intent(inout) :: self
        call self%envelope%kill
        self%l_active = .false.
        self%l_initialized = .false.
    end subroutine flex_env_kill

    logical function flex_env_active( self )
        class(flex_pcg_environment), intent(in) :: self
        flex_env_active = self%l_initialized .and. self%l_active
    end function flex_env_active


    !> the window every basis-shaped volume is multiplied by: the envelope when it is set, else the
    !! soft sphere of the gridding path (mask3D_soft at msk_crop). Typed entry points on purpose: a
    !! class(image) dummy segfaulted at -O3 on the first call of the M-step tail (gfortran codegen).
    subroutine flex_window_apply( self, img, box_crop, msk_crop )
        class(flex_pcg_environment), intent(inout) :: self
        type(image), intent(inout) :: img   !< at the covariance box
        integer,     intent(in)    :: box_crop
        real,        intent(in)    :: msk_crop
        real, pointer :: rp(:,:,:)
        if( .not. self%l_initialized ) THROW_HARD('flex_window_apply: flex environment is not initialized')
        if( self%l_active )then
            if( any(img%get_ldim() /= box_crop) ) THROW_HARD('flex_window_apply: volume is not at the covariance box')
            if( img%is_ft() ) call img%ifft
            call img%get_rmat_ptr(rp)
            call window_product(self, rp, box_crop)
        else
            if( msk_crop > TINY ) call img%mask3D_soft(msk_crop, backgr=0.)
        endif
    end subroutine flex_window_apply

    !> the same for a reconstructor (the initial basis is built on one)
    subroutine flex_window_apply_rec( self, rec, box_crop, msk_crop )
        class(flex_pcg_environment), intent(inout) :: self
        type(reconstructor), intent(inout) :: rec
        integer,             intent(in)    :: box_crop
        real,                intent(in)    :: msk_crop
        real, pointer :: rp(:,:,:)
        if( .not. self%l_initialized ) THROW_HARD('flex_window_apply_rec: flex environment is not initialized')
        if( self%l_active )then
            if( any(rec%get_ldim() /= box_crop) ) THROW_HARD('flex_window_apply_rec: volume is not at the covariance box')
            if( rec%is_ft() ) call rec%ifft
            call rec%get_rmat_ptr(rp)
            call window_product(self, rp, box_crop)
        else
            if( msk_crop > TINY ) call rec%mask3D_soft(msk_crop, backgr=0.)
        endif
    end subroutine flex_window_apply_rec

    !> elementwise product with the envelope on the logical box (loops: no array temporary)
    subroutine window_product( self, rp, b )
        class(flex_pcg_environment), intent(inout) :: self
        real, pointer, intent(inout) :: rp(:,:,:)
        integer,       intent(in)    :: b
        real, pointer :: wp(:,:,:)
        integer :: i, j, k
        if( .not. self%envelope%exists() ) THROW_HARD('window_product: the envelope window is not loaded')
        call self%envelope%get_rmat_ptr(wp)
        do k = 1, b
            do j = 1, b
                do i = 1, b
                    rp(i,j,k) = rp(i,j,k) * wp(i,j,k)
                end do
            end do
        end do
    end subroutine window_product

    !> the solver's support at the covariance box: the envelope when set, else the sphere
    subroutine flex_pcg_install_window( self, op, khi, msk_crop )
        class(flex_pcg_environment), intent(inout) :: self
        type(flex_pcg_t), intent(inout) :: op
        integer,          intent(in)    :: khi
        real,             intent(in)    :: msk_crop
        if( .not. self%l_initialized ) THROW_HARD('flex_pcg_install_window: flex environment is not initialized')
        if( self%l_active )then
            call op%set_window_volume(self%envelope)
            write(logfhandle,'(A,F6.3)') '>>> FLEX_PCA PCG support: envelope window (pcg_mskfile), hard support fraction ', &
                &sum(op%mask)/real(size(op%mask))
        else
            if( msk_crop > TINY ) call op%set_window_sphere(msk_crop)
        endif
        ! the band the planes carry (projected_model_kfromto) plus two shells of margin: the
        ! accumulators and packed kernels of every operator live only on the points this band reaches
        call op%set_band(khi + 2)
    end subroutine flex_pcg_install_window

end module simple_flex_pca_pcg
