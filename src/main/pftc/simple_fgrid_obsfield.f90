!@descr: dense Fourier-grid observation field for restoration experiments
module simple_fgrid_obsfield
use simple_core_module_api
use simple_image, only: image
implicit none

public :: fgrid_obs_field, fgrid_obsfield_eo, unsampled_floor
private
#include "simple_local_flags.inc"

! Two h-colors make nearest-cell writes race-free without atomics: all samples
! for one h are handled by one thread, and h values in the same color are at
! least two native grid units apart before and after rotation.
integer,  parameter :: OBSFIELD_NN_OMP_STRIDE     = 2
integer,  parameter :: OBSFIELD_HEADER_SIZE      = 21
integer,  parameter :: OBSFIELD_GATE_FORMAT_SIGN = -2
real(dp), parameter :: OBSFIELD_UNSAMPLED_FLOOR  = 1.d3

! Part-local Fourier-grid observation field. This is intentionally volume-like:
! particle Fourier components are accumulated into dense expanded-grid
! numerator/density arrays, and requested polar central sections are gathered
! directly from those arrays. Insertion uses weighted nearest-cell assignment:
! each sample updates one native grid cell with a center-normalized KB factor
! from its sub-cell offset. This avoids full insertion-side splatting and the
! hash-table overhead that made the first obsfield prototype slow.
type :: fgrid_obs_field
    private
    type(kbinterpol)       :: kb
    integer                :: pf            = OSMPL_PAD_FAC
    integer                :: nyq           = 0
    integer                :: iwinsz        = 0
    integer                :: wdim          = 0
    integer                :: lims(3,2)     = 0
    integer                :: grid_lims(3,2)= 0
    integer                :: grid_shape(3) = 0
    integer                :: nobs          = 0
    integer                :: ncells        = 0
    logical                :: initialized   = .false.
    complex(dp), allocatable :: grid_num(:,:,:)
    real(dp),    allocatable :: grid_den(:,:,:)
    logical,     allocatable :: grid_obs(:,:,:)
  contains
    procedure, public  :: new                    => obsfield_new
    procedure, public  :: reset                  => obsfield_reset
    procedure, public  :: kill                   => obsfield_kill
    procedure, public  :: insert_plane_oversamp
    procedure, public  :: append_field           => obsfield_append_field
    procedure, public  :: blend_field            => obsfield_blend_field
    procedure, public  :: extract_polar          => obsfield_extract_polar
    procedure, public  :: calc_invtau2           => obsfield_calc_invtau2
    procedure, public  :: log_shell_budget_stats => obsfield_log_shell_budget_stats
    procedure, public  :: count_shell_cells      => obsfield_count_shell_cells
    procedure, public  :: count_shell_cell_counts => obsfield_count_shell_cell_counts
    procedure, public  :: get_nobs               => obsfield_get_nobs
    procedure, public  :: get_ncells             => obsfield_get_ncells
    procedure, private :: compatible_with        => obsfield_compatible_with
end type fgrid_obs_field

! Even/odd wrapper matching the reconstruction convention:
! eo=-1 or eo=0 contributes to even, eo=1 contributes to odd.
type :: fgrid_obsfield_eo
    type(fgrid_obs_field) :: even
    type(fgrid_obs_field) :: odd
  contains
    procedure, public :: new           => obsfield_eo_new
    procedure, public :: reset         => obsfield_eo_reset
    procedure, public :: kill          => obsfield_eo_kill
    procedure, public :: insert_plane  => obsfield_eo_insert_plane
    procedure, public :: append_field  => obsfield_eo_append_field
    procedure, public :: blend_field   => obsfield_eo_blend_field
    procedure, public :: extract_polar => obsfield_eo_extract_polar
    procedure, public :: log_shell_budget_stats => obsfield_eo_log_shell_budget_stats
    procedure, public :: write_merged_volume => obsfield_eo_write_merged_volume
    procedure, public :: write         => obsfield_eo_write
    procedure, public :: read          => obsfield_eo_read
end type fgrid_obsfield_eo

contains

    pure real(dp) function unsampled_floor( grid_den )
        real(dp), intent(in) :: grid_den
        ! Match the Cartesian reconstructor's high-invtau2 fallback for
        ! sampled cells in shells without an FSC-derived ML prior.
        unsampled_floor = min(OBSFIELD_UNSAMPLED_FLOOR, OBSFIELD_UNSAMPLED_FLOOR * grid_den)
    end function unsampled_floor

    subroutine obsfield_new( self, lims, nyq, zero_init )
        class(fgrid_obs_field), intent(inout) :: self
        integer,                intent(in)    :: lims(3,2)
        integer,                intent(in)    :: nyq
        logical, optional,      intent(in)    :: zero_init
        real(dp) :: ncells_dp
        integer  :: dim
        logical  :: l_zero_init
        call self%kill
        l_zero_init = .true.
        if( present(zero_init) ) l_zero_init = zero_init
        if( nyq < 1 ) THROW_HARD('invalid nyq; obsfield_new')
        if( any(lims(:,2) < lims(:,1)) ) THROW_HARD('invalid limits; obsfield_new')
        ! KBALPHA selects the standard kernel shape; obsfield coordinates
        ! passed to the gather remain native grid-cell coordinates.
        self%kb        = kbinterpol(KBWINSZ, KBALPHA)
        self%pf        = OSMPL_PAD_FAC
        self%nyq       = nyq
        self%iwinsz    = ceiling(KBWINSZ - 0.5)
        self%wdim      = self%kb%get_wdim()
        self%lims      = lims
        dim = maxval(abs(lims)) + ceiling(KBWINSZ)
        self%grid_lims(1,:) = [lims(1,1)-self%wdim, dim]
        self%grid_lims(2,:) = [-dim, dim]
        self%grid_lims(3,:) = [-dim, dim]
        self%grid_shape = self%grid_lims(:,2) - self%grid_lims(:,1) + 1
        if( any(self%grid_shape < 1) ) THROW_HARD('invalid grid dimensions; obsfield_new')
        ncells_dp = real(self%grid_shape(1),dp) * real(self%grid_shape(2),dp) * real(self%grid_shape(3),dp)
        if( ncells_dp > real(huge(0),dp) ) THROW_HARD('grid exceeds default integer range; obsfield_new')
        self%ncells = int(ncells_dp)
        if( l_zero_init )then
            allocate(self%grid_num(self%grid_lims(1,1):self%grid_lims(1,2), &
                self%grid_lims(2,1):self%grid_lims(2,2), self%grid_lims(3,1):self%grid_lims(3,2)), &
                source=DCMPLX_ZERO)
            allocate(self%grid_den(self%grid_lims(1,1):self%grid_lims(1,2), &
                self%grid_lims(2,1):self%grid_lims(2,2), self%grid_lims(3,1):self%grid_lims(3,2)), &
                source=0.d0)
            allocate(self%grid_obs(self%grid_lims(1,1):self%grid_lims(1,2), &
                self%grid_lims(2,1):self%grid_lims(2,2), self%grid_lims(3,1):self%grid_lims(3,2)), &
                source=.false.)
        else
            allocate(self%grid_num(self%grid_lims(1,1):self%grid_lims(1,2), &
                self%grid_lims(2,1):self%grid_lims(2,2), self%grid_lims(3,1):self%grid_lims(3,2)))
            allocate(self%grid_den(self%grid_lims(1,1):self%grid_lims(1,2), &
                self%grid_lims(2,1):self%grid_lims(2,2), self%grid_lims(3,1):self%grid_lims(3,2)))
            allocate(self%grid_obs(self%grid_lims(1,1):self%grid_lims(1,2), &
                self%grid_lims(2,1):self%grid_lims(2,2), self%grid_lims(3,1):self%grid_lims(3,2)))
        endif
        self%initialized = .true.
    end subroutine obsfield_new

    subroutine obsfield_reset( self )
        class(fgrid_obs_field), intent(inout) :: self
        integer :: h, k, l
        self%nobs      = 0
        if( allocated(self%grid_num) .and. allocated(self%grid_den) .and. allocated(self%grid_obs) )then
            !$omp parallel do collapse(3) default(shared) schedule(static) private(h,k,l) proc_bind(close)
            do l = lbound(self%grid_num,3), ubound(self%grid_num,3)
                do k = lbound(self%grid_num,2), ubound(self%grid_num,2)
                    do h = lbound(self%grid_num,1), ubound(self%grid_num,1)
                        self%grid_num(h,k,l) = DCMPLX_ZERO
                        self%grid_den(h,k,l) = 0.d0
                        self%grid_obs(h,k,l) = .false.
                    enddo
                enddo
            enddo
            !$omp end parallel do
        else
            if( allocated(self%grid_num) ) self%grid_num = DCMPLX_ZERO
            if( allocated(self%grid_den) ) self%grid_den = 0.d0
            if( allocated(self%grid_obs) ) self%grid_obs = .false.
        endif
    end subroutine obsfield_reset

    subroutine obsfield_kill( self )
        class(fgrid_obs_field), intent(inout) :: self
        if( allocated(self%grid_num) ) deallocate(self%grid_num)
        if( allocated(self%grid_den) ) deallocate(self%grid_den)
        if( allocated(self%grid_obs) ) deallocate(self%grid_obs)
        self%pf          = OSMPL_PAD_FAC
        self%nyq         = 0
        self%iwinsz      = 0
        self%wdim        = 0
        self%lims        = 0
        self%grid_lims   = 0
        self%grid_shape  = 0
        self%nobs        = 0
        self%ncells      = 0
        self%initialized = .false.
    end subroutine obsfield_kill

    ! Insert one particle Fourier plane by weighted nearest-cell assignment.
    ! Coordinates are computed exactly as in reconstructor%insert_plane_oversamp,
    ! but each experimental component updates only the closest native 3D grid cell.
    !
    ! Nearest-cell insertion touches one destination cell per Fourier sample.
    ! Use two h-colors to avoid write collisions without OpenMP atomics.
    subroutine insert_plane_oversamp( self, se, o, fpl, pwght, shift_crop )
        use simple_math,    only: ceil_div, floor_div
        use simple_math_ft, only: fplane_get_cmplx, fplane_get_ctfsq
        class(fgrid_obs_field), intent(inout) :: self
        class(sym),             intent(inout) :: se
        class(ori),             intent(inout) :: o
        class(fplane_type),     intent(in)    :: fpl
        real,                   intent(in)    :: pwght
        real, optional,         intent(in)    :: shift_crop(2)
        type(ori)   :: o_sym
        complex(sp) :: cmplx_raw
        complex(sp) :: phase_h, phase_k, phase_shift, w_k
        complex(dp) :: comp
        real(dp)    :: ctfval, pwght_dp, pwght_pf2_dp, cell_w, cell_w_norm, pshift_pf(2)
        real(sp)    :: ctfsq_raw
        real(sp)    :: loc(3), delta(3), R(3,3)
        real        :: rotmats(se%get_nsym(),3,3)
        integer     :: fpllims_pd(3,2), fpllims(3,2), coord(3)
        integer     :: nsym, isym, h, k, hp, kp, pf_local, l
        integer     :: nyq_disk, h_sq, k_max_h, k_lo, k_hi
        integer     :: lim1_lo, gl_lo1, gl_lo2, gl_lo3, gl_hi1, gl_hi2, gl_hi3
        integer     :: nobs_add
        logical     :: l_shift
        if( .not. self%initialized ) THROW_HARD('obsfield not initialized; insert_plane_oversamp')
        if( pwght < TINY ) return
        nobs_add      = 0
        ! rotation matrices (one per sym) in native Fourier-grid units
        nsym = se%get_nsym()
        rotmats(1,:,:) = o%get_mat()
        if( nsym > 1 )then
            do isym = 2, nsym
                call se%apply(o, isym, o_sym)
                rotmats(isym,:,:) = o_sym%get_mat()
            enddo
        endif
        ! native iteration limits so hp=h*pf and kp=k*pf fit the padded plane
        pf_local     = self%pf
        fpllims_pd   = fpl%frlims
        fpllims      = fpllims_pd
        fpllims(1,1) = ceil_div (fpllims_pd(1,1), pf_local)
        fpllims(1,2) = floor_div(fpllims_pd(1,2), pf_local)
        fpllims(2,1) = ceil_div (fpllims_pd(2,1), pf_local)
        fpllims(2,2) = floor_div(fpllims_pd(2,2), pf_local)
        ! hoisted constants
        pwght_dp     = real(pwght, dp)
        pwght_pf2_dp = pwght_dp * real(pf_local*pf_local, dp)
        cell_w_norm  = real(self%kb%apod_fast(0._sp), dp)
        cell_w_norm  = 1.d0 / (cell_w_norm * cell_w_norm * cell_w_norm)
        l_shift      = .false.
        if( present(shift_crop) ) l_shift = any(abs(shift_crop) > 1.e-6)
        if( l_shift )then
            pshift_pf = real(-shift_crop * fpl%shconst(1:2), dp) * real(pf_local, dp)
            w_k       = cmplx(real(cos(pshift_pf(2)),sp), real(sin(pshift_pf(2)),sp), kind=sp)
        endif
        ! integer disk gate: bit-exact equivalent of original nint(sqrt(h^2+k^2)) > nyq
        ! since for non-negative integer n,  n > nyq*(nyq+1)  <=>  nint(sqrt(n)) > nyq
        nyq_disk     = self%nyq * (self%nyq + 1)
        lim1_lo      = self%lims(1,1)
        gl_lo1 = self%grid_lims(1,1); gl_hi1 = self%grid_lims(1,2)
        gl_lo2 = self%grid_lims(2,1); gl_hi2 = self%grid_lims(2,2)
        gl_lo3 = self%grid_lims(3,1); gl_hi3 = self%grid_lims(3,2)
        !$omp parallel default(shared) private(isym,R,l,h,k,h_sq,k_max_h,k_lo,k_hi,&
        !$omp& cmplx_raw,ctfsq_raw,comp,ctfval,loc,delta,coord,hp,kp,cell_w)&
        !$omp& private(phase_h,phase_k,phase_shift)&
        !$omp& reduction(+:nobs_add) proc_bind(close)
        do isym = 1, nsym
            R = rotmats(isym,:,:)
            do l = 0, OBSFIELD_NN_OMP_STRIDE-1
                !$omp do schedule(static)
                do h = fpllims(1,1)+l, fpllims(1,2), OBSFIELD_NN_OMP_STRIDE
                    h_sq = h*h
                    if( h_sq > nyq_disk ) cycle
                    ! tightest k-range satisfying the integer disk gate
                    k_max_h = int(sqrt(real(nyq_disk - h_sq, sp)))
                    k_lo    = max(fpllims(2,1), -k_max_h)
                    k_hi    = min(fpllims(2,2),  k_max_h)
                    hp = h * pf_local
                    if( l_shift )then
                        phase_h = cmplx(real(cos(real(h,dp)*pshift_pf(1)),sp), &
                            real(sin(real(h,dp)*pshift_pf(1)),sp), kind=sp)
                        phase_k = cmplx(real(cos(real(k_lo,dp)*pshift_pf(2)),sp), &
                            real(sin(real(k_lo,dp)*pshift_pf(2)),sp), kind=sp)
                    endif
                    do k = k_lo, k_hi
                        kp = k * pf_local
                        if( l_shift )then
                            phase_shift = phase_h * phase_k
                            phase_k     = phase_k * w_k
                        endif
                        ! raw padded-plane samples (single precision); cheap zero-skip
                        cmplx_raw = fplane_get_cmplx(fpl, hp, kp)
                        ctfsq_raw = fplane_get_ctfsq(fpl, hp, kp)
                        if( abs(real(cmplx_raw)) + abs(aimag(cmplx_raw)) <= TINY .and. &
                            ctfsq_raw <= TINY ) cycle
                        if( l_shift ) cmplx_raw = cmplx_raw * phase_shift
                        ! rotated location on the native lattice; third h-component
                        ! is zero, so 6 muls instead of matmul's 9
                        loc(1) = real(h,sp)*R(1,1) + real(k,sp)*R(2,1)
                        loc(2) = real(h,sp)*R(1,2) + real(k,sp)*R(2,2)
                        loc(3) = real(h,sp)*R(1,3) + real(k,sp)*R(2,3)
                        coord(1) = nint(loc(1))
                        ! Friedel-mate skip: redundant negative first-axis cells are
                        ! recovered inside obsfield_extract_polar (addr -> -addr,
                        ! conjugate grid_num).
                        if( coord(1) < lim1_lo ) cycle
                        if( coord(1) > gl_hi1 ) cycle
                        coord(2) = nint(loc(2))
                        if( coord(2) < gl_lo2 .or. coord(2) > gl_hi2 ) cycle
                        coord(3) = nint(loc(3))
                        if( coord(3) < gl_lo3 .or. coord(3) > gl_hi3 ) cycle
                        ! Center-normalized one-cell KB confidence. This weights the
                        ! nearest-cell deposit without spreading into neighboring cells.
                        delta  = loc - real(coord,sp)
                        cell_w = cell_w_norm * real(self%kb%apod_fast(delta(1)), dp) &
                            * real(self%kb%apod_fast(delta(2)), dp) &
                            * real(self%kb%apod_fast(delta(3)), dp)
                        ! promote to dp once and accumulate
                        comp   = cell_w * pwght_pf2_dp * cmplx(cmplx_raw, kind=dp)
                        ctfval = cell_w * pwght_dp * real(ctfsq_raw, dp)
                        self%grid_num(coord(1),coord(2),coord(3)) = &
                            self%grid_num(coord(1),coord(2),coord(3)) + comp
                        self%grid_den(coord(1),coord(2),coord(3)) = &
                            self%grid_den(coord(1),coord(2),coord(3)) + ctfval
                        self%grid_obs(coord(1),coord(2),coord(3)) = .true.
                        nobs_add = nobs_add + 1
                    enddo
                enddo
                !$omp end do
            enddo
        enddo
        !$omp end parallel
        self%nobs = self%nobs + nobs_add
        if( nsym > 1 ) call o_sym%kill
    end subroutine insert_plane_oversamp

    subroutine obsfield_append_field( self, src )
        class(fgrid_obs_field), intent(inout) :: self
        class(fgrid_obs_field), intent(in)    :: src
        integer :: h, k, l, nobs_new
        if( .not. src%initialized ) return
        if( .not. self%initialized ) THROW_HARD('destination not initialized; obsfield_append_field')
        if( .not. self%compatible_with(src) ) THROW_HARD('incompatible observation fields; obsfield_append_field')
        nobs_new = 0
        !$omp parallel do collapse(3) default(shared) schedule(static) private(h,k,l) reduction(+:nobs_new) proc_bind(close)
        do l = lbound(self%grid_num,3), ubound(self%grid_num,3)
            do k = lbound(self%grid_num,2), ubound(self%grid_num,2)
                do h = lbound(self%grid_num,1), ubound(self%grid_num,1)
                    self%grid_num(h,k,l) = self%grid_num(h,k,l) + src%grid_num(h,k,l)
                    self%grid_den(h,k,l) = self%grid_den(h,k,l) + src%grid_den(h,k,l)
                    self%grid_obs(h,k,l) = self%grid_obs(h,k,l) .or. src%grid_obs(h,k,l)
                    if( self%grid_obs(h,k,l) ) nobs_new = nobs_new + 1
                enddo
            enddo
        enddo
        !$omp end parallel do
        self%nobs = nobs_new
    end subroutine obsfield_append_field

    subroutine obsfield_blend_field( self, prev, update_frac )
        class(fgrid_obs_field), intent(inout) :: self
        class(fgrid_obs_field), intent(in)    :: prev
        real(dp),                intent(in)    :: update_frac
        real(dp) :: weight_prev
        integer  :: lo(3), hi(3)
        integer  :: h, k, l, nobs_new
        logical  :: l_intersection_subset
        if( .not. self%initialized ) THROW_HARD('destination not initialized; obsfield_blend_field')
        if( .not. prev%initialized ) return
        if( self%pf /= prev%pf ) THROW_HARD('incompatible observation field oversampling; obsfield_blend_field')
        if( self%nyq < 1 .or. prev%nyq < 1 ) THROW_HARD('invalid observation field nyquist; obsfield_blend_field')
        weight_prev = 1.d0 - update_frac
        lo = max(self%grid_lims(:,1), prev%grid_lims(:,1))
        hi = min(self%grid_lims(:,2), prev%grid_lims(:,2))
        if( any(hi < lo) )then
            self%nobs = count(self%grid_obs)
            return
        endif
        l_intersection_subset = any(lo /= self%grid_lims(:,1)) .or. any(hi /= self%grid_lims(:,2)) .or. &
            &any(lo /= prev%grid_lims(:,1)) .or. any(hi /= prev%grid_lims(:,2))
        if( l_intersection_subset .or. self%nyq /= prev%nyq .or. any(self%lims /= prev%lims) )then
            write(logfhandle,'(A,2(A,I0),A,6(I0,1X),A,6(I0,1X),A,6(I0,1X))') &
                &'WARNING: obsfield trailing blend uses grid intersection;', &
                &' current_nyq=', self%nyq, ' prev_nyq=', prev%nyq, &
                &' current_grid_lims=', self%grid_lims, ' prev_grid_lims=', prev%grid_lims, ' overlap=', lo, hi
        endif
        nobs_new = 0
        !$omp parallel do collapse(3) default(shared) schedule(static) private(h,k,l) reduction(+:nobs_new) proc_bind(close)
        do l = lo(3), hi(3)
            do k = lo(2), hi(2)
                do h = lo(1), hi(1)
                    self%grid_num(h,k,l) = update_frac * self%grid_num(h,k,l) + weight_prev * prev%grid_num(h,k,l)
                    self%grid_den(h,k,l) = update_frac * self%grid_den(h,k,l) + weight_prev * prev%grid_den(h,k,l)
                    self%grid_obs(h,k,l) = ((update_frac > DTINY) .and. self%grid_obs(h,k,l)) .or. &
                        &((weight_prev > DTINY) .and. prev%grid_obs(h,k,l))
                    if( self%grid_obs(h,k,l) ) nobs_new = nobs_new + 1
                enddo
            enddo
        enddo
        !$omp end parallel do
        if( l_intersection_subset )then
            self%nobs = count(self%grid_obs)
        else
            self%nobs = nobs_new
        endif
    end subroutine obsfield_blend_field

    subroutine obsfield_calc_invtau2( self, kfromto, fsc, fudge, prior_start, invtau2 )
        class(fgrid_obs_field), intent(in)  :: self
        integer,                intent(in)  :: kfromto(2), prior_start
        real(dp),               intent(in)  :: fsc(kfromto(1):kfromto(2))
        real(dp),               intent(in)  :: fudge
        real(dp),               intent(out) :: invtau2(kfromto(1):kfromto(2))
        real(dp), allocatable :: rsum(:), sig2(:), ssnr(:), tau2(:)
        integer,  allocatable :: cnt(:)
        real(dp) :: cc
        integer  :: h, k, l, shell
        invtau2 = 0.d0
        if( .not. self%initialized ) return
        if( kfromto(1) > kfromto(2) ) THROW_HARD('invalid k-range; obsfield_calc_invtau2')
        allocate(rsum(kfromto(1):kfromto(2)), sig2(kfromto(1):kfromto(2)), &
            &ssnr(kfromto(1):kfromto(2)), tau2(kfromto(1):kfromto(2)), source=0.d0)
        allocate(cnt(kfromto(1):kfromto(2)), source=0)
        !$omp parallel do collapse(3) default(shared) schedule(static) proc_bind(close) private(h,k,l,shell) reduction(+:cnt,rsum)
        do h = self%lims(1,1), self%lims(1,2)
            do k = self%lims(2,1), self%lims(2,2)
                do l = self%lims(3,1), self%lims(3,2)
                    shell = nint(sqrt(real(h*h + k*k + l*l, dp)))
                    if( shell < kfromto(1) .or. shell > kfromto(2) ) cycle
                    cnt(shell)  = cnt(shell) + 1
                    rsum(shell) = rsum(shell) + self%grid_den(h,k,l)
                enddo
            enddo
        enddo
        !$omp end parallel do
        do shell = kfromto(1), kfromto(2)
            cc = max(0.001d0, min(0.999d0, fsc(shell)))
            ssnr(shell) = cc / (1.d0 - cc)
        enddo
        where( rsum > DTINY )
            sig2 = real(cnt,dp) / rsum
        elsewhere
            sig2 = 0.d0
        end where
        tau2 = ssnr * sig2
        do shell = max(kfromto(1),prior_start), kfromto(2)
            if( tau2(shell) > DTINY )then
                invtau2(shell) = 1.d0 / (fudge * tau2(shell))
            endif
        enddo
        deallocate(rsum, sig2, ssnr, tau2, cnt)
    end subroutine obsfield_calc_invtau2

    subroutine obsfield_log_shell_budget_stats( self, kfromto, shell_nodes, label )
        class(fgrid_obs_field), intent(in) :: self
        integer,                intent(in) :: kfromto(2), shell_nodes(:)
        character(len=*),       intent(in) :: label
        call log_single_field_shell_budget(self, kfromto, shell_nodes, label, 'single')
    end subroutine obsfield_log_shell_budget_stats

    integer function obsfield_count_shell_cells( self, kfromto )
        class(fgrid_obs_field), intent(in) :: self
        integer,                intent(in) :: kfromto(2)
        integer, allocatable :: shell_cells(:)
        obsfield_count_shell_cells = 0
        call self%count_shell_cell_counts(kfromto, shell_cells)
        if( allocated(shell_cells) )then
            obsfield_count_shell_cells = sum(shell_cells)
            deallocate(shell_cells)
        endif
    end function obsfield_count_shell_cells

    subroutine obsfield_count_shell_cell_counts( self, kfromto, shell_cells )
        class(fgrid_obs_field), intent(in)  :: self
        integer,                intent(in)  :: kfromto(2)
        integer, allocatable,   intent(out) :: shell_cells(:)
        integer :: h, k, l, shell, ik, nk
        nk = kfromto(2) - kfromto(1) + 1
        if( nk < 1 )then
            allocate(shell_cells(0))
            return
        endif
        allocate(shell_cells(nk), source=0)
        if( .not. self%initialized )then
            return
        endif
        !$omp parallel do collapse(3) default(shared) private(h,k,l,shell,ik) reduction(+:shell_cells) schedule(static) proc_bind(close)
        do l = self%lims(3,1), self%lims(3,2)
            do k = self%lims(2,1), self%lims(2,2)
                do h = self%lims(1,1), self%lims(1,2)
                    shell = nint(sqrt(real(h*h + k*k + l*l, dp)))
                    if( shell < kfromto(1) .or. shell > kfromto(2) ) cycle
                    ik = shell - kfromto(1) + 1
                    shell_cells(ik) = shell_cells(ik) + 1
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine obsfield_count_shell_cell_counts

    ! Reconstructor-like obsfield extraction: restore each Cartesian cell as
    ! grid_num / (grid_den + prior(shell)) before KB gathering. This avoids the
    ! non-equivalent polar-space pattern gather(num)/(gather(den)+prior).
    subroutine obsfield_extract_polar( self, eulspace, nrefs, kfromto, polar_x, polar_y, invtau2, prior_start, pfts, denw )
        class(fgrid_obs_field), intent(in)    :: self
        class(oris),            intent(in)    :: eulspace
        integer,                intent(in)    :: nrefs, kfromto(2)
        real(sp),               intent(in)    :: polar_x(:,:), polar_y(:,:)
        real(dp),               intent(in)    :: invtau2(kfromto(1):kfromto(2))
        integer,                intent(in)    :: prior_start
        complex(dp),            intent(inout) :: pfts(:,:,:)
        real(dp), optional,     intent(inout) :: denw(:,:,:)
        complex(dp) :: acc, cell_num
        real(dp)    :: denom, denacc, prior, wd
        real(sp)    :: R(3,3), loc(3), px, py
        real(sp), allocatable :: w(:,:,:)
        integer     :: iproj, irot, k, kloc, pftsz, kfromto1
        integer     :: win1_1, win1_2, win1_3, c1, c2, c3, a1, a2, a3
        integer     :: l, m, n, iwinsz_l, wdim_l, lim1_lo, shell
        integer     :: gl_lo1, gl_lo2, gl_lo3, gl_hi1, gl_hi2, gl_hi3
        logical     :: l_conj
        pftsz = size(polar_x,1)
        pfts(:,:,1:nrefs) = DCMPLX_ZERO
        if( present(denw) ) denw(:,:,1:nrefs) = 0.d0
        if( .not. any(self%grid_obs) ) return
        iwinsz_l = self%iwinsz
        wdim_l   = self%wdim
        lim1_lo  = self%lims(1,1)
        kfromto1 = kfromto(1)
        gl_lo1   = self%grid_lims(1,1); gl_hi1 = self%grid_lims(1,2)
        gl_lo2   = self%grid_lims(2,1); gl_hi2 = self%grid_lims(2,2)
        gl_lo3   = self%grid_lims(3,1); gl_hi3 = self%grid_lims(3,2)
        !$omp parallel default(shared) private(iproj,R,k,kloc,irot,px,py,loc,w,acc,denacc,wd,denom,prior,&
        !$omp& win1_1,win1_2,win1_3,c1,c2,c3,a1,a2,a3,l,m,n,l_conj,cell_num,shell) proc_bind(close)
        allocate(w(wdim_l,wdim_l,wdim_l))
        !$omp do schedule(static)
        do iproj = 1, nrefs
            R = eulspace%get_mat(iproj)
            do k = kfromto1, kfromto(2)
                kloc = k - kfromto1 + 1
                do irot = 1, pftsz
                    px     = polar_x(irot,kloc)
                    py     = polar_y(irot,kloc)
                    loc(1) = px*R(1,1) + py*R(2,1)
                    loc(2) = px*R(1,2) + py*R(2,2)
                    loc(3) = px*R(1,3) + py*R(2,3)
                    win1_1 = nint(loc(1)) - iwinsz_l
                    win1_2 = nint(loc(2)) - iwinsz_l
                    win1_3 = nint(loc(3)) - iwinsz_l
                    call self%kb%apod_mat_3d_fast(loc, iwinsz_l, wdim_l, w)
                    acc = DCMPLX_ZERO
                    denacc = 0.d0
                    do n = 1, wdim_l
                        c3 = win1_3 + n - 1
                        if( c3 < gl_lo3 .or. c3 > gl_hi3 ) cycle
                        do m = 1, wdim_l
                            c2 = win1_2 + m - 1
                            if( c2 < gl_lo2 .or. c2 > gl_hi2 ) cycle
                            do l = 1, wdim_l
                                c1     = win1_1 + l - 1
                                wd     = real(w(l,m,n), dp)
                                l_conj = c1 < lim1_lo
                                if( l_conj )then
                                    a1 = -c1; a2 = -c2; a3 = -c3
                                    if( a1 < gl_lo1 .or. a1 > gl_hi1 ) cycle
                                    if( .not. self%grid_obs(a1,a2,a3) ) cycle
                                    shell = nint(sqrt(real(a1*a1 + a2*a2 + a3*a3, dp)))
                                    prior = 0.d0
                                    if( shell >= kfromto(1) .and. shell <= kfromto(2) ) prior = invtau2(shell)
                                    denom = self%grid_den(a1,a2,a3) + prior
                                    if( shell >= kfromto(1) .and. shell <= kfromto(2) .and. &
                                        &shell >= prior_start .and. prior <= DTINY )then
                                        denom = denom + unsampled_floor(self%grid_den(a1,a2,a3))
                                    endif
                                    if( denom <= DTINY ) cycle
                                    cell_num = conjg(self%grid_num(a1,a2,a3))
                                else
                                    if( c1 > gl_hi1 ) cycle
                                    if( .not. self%grid_obs(c1,c2,c3) ) cycle
                                    shell = nint(sqrt(real(c1*c1 + c2*c2 + c3*c3, dp)))
                                    prior = 0.d0
                                    if( shell >= kfromto(1) .and. shell <= kfromto(2) ) prior = invtau2(shell)
                                    denom = self%grid_den(c1,c2,c3) + prior
                                    if( shell >= kfromto(1) .and. shell <= kfromto(2) .and. &
                                        &shell >= prior_start .and. prior <= DTINY )then
                                        denom = denom + unsampled_floor(self%grid_den(c1,c2,c3))
                                    endif
                                    if( denom <= DTINY ) cycle
                                    cell_num = self%grid_num(c1,c2,c3)
                                endif
                                acc    = acc + wd * cell_num / denom
                                denacc = denacc + wd * denom
                            enddo
                        enddo
                    enddo
                    pfts(irot,kloc,iproj) = acc
                    if( present(denw) ) denw(irot,kloc,iproj) = denacc
                enddo
            enddo
        enddo
        !$omp end do
        deallocate(w)
        !$omp end parallel
    end subroutine obsfield_extract_polar

    ! Gather even/odd half-map references through one shared geometry/KB walk,
    ! then derive the merged reference with the same denominator weighting used
    ! by the separate extraction path.
    subroutine obsfield_eo_extract_polar( self, eulspace, nrefs, kfromto, polar_x, polar_y, &
            &invtau2_even, invtau2_odd, prior_start, pfts_even, pfts_odd, pfts_merg )
        class(fgrid_obsfield_eo), intent(in)    :: self
        class(oris),              intent(in)    :: eulspace
        integer,                  intent(in)    :: nrefs, kfromto(2)
        real(sp),                 intent(in)    :: polar_x(:,:), polar_y(:,:)
        real(dp),                 intent(in)    :: invtau2_even(kfromto(1):kfromto(2))
        real(dp),                 intent(in)    :: invtau2_odd( kfromto(1):kfromto(2))
        integer,                  intent(in)    :: prior_start
        complex(dp),              intent(inout) :: pfts_even(:,:,:), pfts_odd(:,:,:), pfts_merg(:,:,:)
        complex(dp) :: acc_even, acc_odd, cell_num
        real(dp)    :: den_even, den_odd, denom, denom_weight, prior, wd
        real(sp)    :: R(3,3), loc(3), px, py
        real(sp), allocatable :: w(:,:,:)
        integer     :: iproj, irot, k, kloc, pftsz, kfromto1
        integer     :: win1_1, win1_2, win1_3, c1, c2, c3, a1, a2, a3
        integer     :: l, m, n, iwinsz_l, wdim_l, lim1_lo, shell
        integer     :: gl_lo1, gl_lo2, gl_lo3, gl_hi1, gl_hi2, gl_hi3
        logical     :: have_even, have_odd, l_conj
        pftsz = size(polar_x,1)
        pfts_even(:,:,1:nrefs) = DCMPLX_ZERO
        pfts_odd( :,:,1:nrefs) = DCMPLX_ZERO
        pfts_merg(:,:,1:nrefs) = DCMPLX_ZERO
        have_even = any(self%even%grid_obs)
        have_odd  = any(self%odd%grid_obs)
        if( .not. have_even .and. .not. have_odd ) return
        iwinsz_l = self%even%iwinsz
        wdim_l   = self%even%wdim
        lim1_lo  = self%even%lims(1,1)
        kfromto1 = kfromto(1)
        gl_lo1   = self%even%grid_lims(1,1); gl_hi1 = self%even%grid_lims(1,2)
        gl_lo2   = self%even%grid_lims(2,1); gl_hi2 = self%even%grid_lims(2,2)
        gl_lo3   = self%even%grid_lims(3,1); gl_hi3 = self%even%grid_lims(3,2)
        !$omp parallel default(shared) private(iproj,R,k,kloc,irot,px,py,loc,w,acc_even,acc_odd,den_even,den_odd,&
        !$omp& denom,denom_weight,prior,wd,win1_1,win1_2,win1_3,c1,c2,c3,a1,a2,a3,l,m,n,l_conj,cell_num,shell) &
        !$omp& proc_bind(close)
        allocate(w(wdim_l,wdim_l,wdim_l))
        !$omp do schedule(static)
        do iproj = 1, nrefs
            R = eulspace%get_mat(iproj)
            do k = kfromto1, kfromto(2)
                kloc = k - kfromto1 + 1
                do irot = 1, pftsz
                    px     = polar_x(irot,kloc)
                    py     = polar_y(irot,kloc)
                    loc(1) = px*R(1,1) + py*R(2,1)
                    loc(2) = px*R(1,2) + py*R(2,2)
                    loc(3) = px*R(1,3) + py*R(2,3)
                    win1_1 = nint(loc(1)) - iwinsz_l
                    win1_2 = nint(loc(2)) - iwinsz_l
                    win1_3 = nint(loc(3)) - iwinsz_l
                    call self%even%kb%apod_mat_3d_fast(loc, iwinsz_l, wdim_l, w)
                    acc_even = DCMPLX_ZERO
                    acc_odd  = DCMPLX_ZERO
                    den_even = 0.d0
                    den_odd  = 0.d0
                    do n = 1, wdim_l
                        c3 = win1_3 + n - 1
                        if( c3 < gl_lo3 .or. c3 > gl_hi3 ) cycle
                        do m = 1, wdim_l
                            c2 = win1_2 + m - 1
                            if( c2 < gl_lo2 .or. c2 > gl_hi2 ) cycle
                            do l = 1, wdim_l
                                c1 = win1_1 + l - 1
                                wd = real(w(l,m,n), dp)
                                l_conj = c1 < lim1_lo
                                if( l_conj )then
                                    a1 = -c1; a2 = -c2; a3 = -c3
                                    if( a1 < gl_lo1 .or. a1 > gl_hi1 ) cycle
                                else
                                    if( c1 > gl_hi1 ) cycle
                                    a1 = c1; a2 = c2; a3 = c3
                                endif
                                shell = nint(sqrt(real(a1*a1 + a2*a2 + a3*a3, dp)))
                                if( have_even .and. self%even%grid_obs(a1,a2,a3) )then
                                    prior = 0.d0
                                    if( shell >= kfromto(1) .and. shell <= kfromto(2) ) prior = invtau2_even(shell)
                                    denom = self%even%grid_den(a1,a2,a3) + prior
                                    if( shell >= kfromto(1) .and. shell <= kfromto(2) .and. &
                                        &shell >= prior_start .and. prior <= DTINY )then
                                        denom = denom + unsampled_floor(self%even%grid_den(a1,a2,a3))
                                    endif
                                    if( denom > DTINY )then
                                        cell_num = self%even%grid_num(a1,a2,a3)
                                        if( l_conj ) cell_num = conjg(cell_num)
                                        acc_even = acc_even + wd * cell_num / denom
                                        den_even = den_even + wd * denom
                                    endif
                                endif
                                if( have_odd .and. self%odd%grid_obs(a1,a2,a3) )then
                                    prior = 0.d0
                                    if( shell >= kfromto(1) .and. shell <= kfromto(2) ) prior = invtau2_odd(shell)
                                    denom = self%odd%grid_den(a1,a2,a3) + prior
                                    if( shell >= kfromto(1) .and. shell <= kfromto(2) .and. &
                                        &shell >= prior_start .and. prior <= DTINY )then
                                        denom = denom + unsampled_floor(self%odd%grid_den(a1,a2,a3))
                                    endif
                                    if( denom > DTINY )then
                                        cell_num = self%odd%grid_num(a1,a2,a3)
                                        if( l_conj ) cell_num = conjg(cell_num)
                                        acc_odd = acc_odd + wd * cell_num / denom
                                        den_odd = den_odd + wd * denom
                                    endif
                                endif
                            enddo
                        enddo
                    enddo
                    pfts_even(irot,kloc,iproj) = acc_even
                    pfts_odd( irot,kloc,iproj) = acc_odd
                    denom_weight = den_even + den_odd
                    if( denom_weight > DTINY )then
                        pfts_merg(irot,kloc,iproj) = (acc_even * den_even + acc_odd * den_odd) / denom_weight
                    endif
                enddo
            enddo
        enddo
        !$omp end do
        deallocate(w)
        !$omp end parallel
    end subroutine obsfield_eo_extract_polar

    subroutine obsfield_eo_log_shell_budget_stats( self, kfromto, shell_nodes, label )
        class(fgrid_obsfield_eo), intent(in) :: self
        integer,                  intent(in) :: kfromto(2), shell_nodes(:)
        character(len=*),         intent(in) :: label
        call log_single_field_shell_budget(self%even, kfromto, shell_nodes, label, 'even')
        call log_single_field_shell_budget(self%odd,  kfromto, shell_nodes, label, 'odd')
        call log_merged_field_shell_budget(self, kfromto, shell_nodes, label, 'merged')
    end subroutine obsfield_eo_log_shell_budget_stats

    subroutine log_single_field_shell_budget( self, kfromto, shell_nodes, label, half )
        class(fgrid_obs_field), intent(in) :: self
        integer,                intent(in) :: kfromto(2), shell_nodes(:)
        character(len=*),       intent(in) :: label, half
        real(dp), allocatable :: den_sum(:), den_min(:), den_max(:)
        integer,  allocatable :: obs_count(:), cell_count(:)
        real(dp) :: coverage, den_mean
        integer  :: h, k, l, shell, ik, nk, node_count
        if( .not. self%initialized ) return
        nk = kfromto(2) - kfromto(1) + 1
        if( nk < 1 .or. size(shell_nodes) /= nk ) return
        allocate(obs_count(nk), cell_count(nk), source=0)
        allocate(den_sum(nk), source=0.d0)
        allocate(den_min(nk), source=huge(0.d0))
        allocate(den_max(nk), source=0.d0)
        !$omp parallel do collapse(3) default(shared) private(h,k,l,shell,ik) &
        !$omp reduction(+:obs_count,cell_count,den_sum) reduction(min:den_min) reduction(max:den_max) schedule(static) proc_bind(close)
        do l = self%lims(3,1), self%lims(3,2)
            do k = self%lims(2,1), self%lims(2,2)
                do h = self%lims(1,1), self%lims(1,2)
                    shell = nint(sqrt(real(h*h + k*k + l*l, dp)))
                    if( shell < kfromto(1) .or. shell > kfromto(2) ) cycle
                    ik = shell - kfromto(1) + 1
                    cell_count(ik) = cell_count(ik) + 1
                    if( .not. self%grid_obs(h,k,l) ) cycle
                    obs_count(ik) = obs_count(ik) + 1
                    den_sum(ik) = den_sum(ik) + self%grid_den(h,k,l)
                    den_min(ik) = min(den_min(ik), self%grid_den(h,k,l))
                    den_max(ik) = max(den_max(ik), self%grid_den(h,k,l))
                enddo
            enddo
        enddo
        !$omp end parallel do
        call log_shell_budget_summary(label, half, kfromto, shell_nodes, cell_count, obs_count, den_sum)
        do ik = 1,nk
            node_count = max(1, shell_nodes(ik))
            coverage = real(obs_count(ik),dp) / real(node_count,dp)
            if( obs_count(ik) > 0 )then
                den_mean = den_sum(ik) / real(obs_count(ik),dp)
            else
                den_mean = 0.d0
                den_min(ik) = 0.d0
            endif
            write(logfhandle,'(A,1X,A,2X,A,1X,A,2X,A,I4,2X,A,I0,2X,A,I0,2X,A,I0,2X,A,F8.3,2X,A,ES12.4,2X,A,ES12.4,2X,A,ES12.4,2X,A,ES12.4)') &
                &'obsfield shell-accum', trim(label), 'half=', trim(half), 'k=', kfromto(1)+ik-1, &
                &'field_nodes=', shell_nodes(ik), 'cart_cells=', cell_count(ik), 'obs_cells=', obs_count(ik), &
                &'obs_per_node=', coverage, 'den_sum=', den_sum(ik), 'den_mean=', den_mean, &
                &'den_min=', den_min(ik), 'den_max=', den_max(ik)
        enddo
        deallocate(obs_count, cell_count, den_sum, den_min, den_max)
    end subroutine log_single_field_shell_budget

    subroutine log_merged_field_shell_budget( self, kfromto, shell_nodes, label, half )
        class(fgrid_obsfield_eo), intent(in) :: self
        integer,                  intent(in) :: kfromto(2), shell_nodes(:)
        character(len=*),         intent(in) :: label, half
        real(dp), allocatable :: den_sum(:), den_min(:), den_max(:)
        integer,  allocatable :: obs_count(:), cell_count(:)
        real(dp) :: coverage, den_here, den_mean
        integer  :: h, k, l, shell, ik, nk, node_count
        if( .not. self%even%initialized .or. .not. self%odd%initialized ) return
        nk = kfromto(2) - kfromto(1) + 1
        if( nk < 1 .or. size(shell_nodes) /= nk ) return
        allocate(obs_count(nk), cell_count(nk), source=0)
        allocate(den_sum(nk), source=0.d0)
        allocate(den_min(nk), source=huge(0.d0))
        allocate(den_max(nk), source=0.d0)
        !$omp parallel do collapse(3) default(shared) private(h,k,l,shell,ik,den_here) &
        !$omp reduction(+:obs_count,cell_count,den_sum) reduction(min:den_min) reduction(max:den_max) schedule(static) proc_bind(close)
        do l = self%even%lims(3,1), self%even%lims(3,2)
            do k = self%even%lims(2,1), self%even%lims(2,2)
                do h = self%even%lims(1,1), self%even%lims(1,2)
                    shell = nint(sqrt(real(h*h + k*k + l*l, dp)))
                    if( shell < kfromto(1) .or. shell > kfromto(2) ) cycle
                    ik = shell - kfromto(1) + 1
                    cell_count(ik) = cell_count(ik) + 1
                    if( .not. (self%even%grid_obs(h,k,l) .or. self%odd%grid_obs(h,k,l)) ) cycle
                    den_here = self%even%grid_den(h,k,l) + self%odd%grid_den(h,k,l)
                    obs_count(ik) = obs_count(ik) + 1
                    den_sum(ik) = den_sum(ik) + den_here
                    den_min(ik) = min(den_min(ik), den_here)
                    den_max(ik) = max(den_max(ik), den_here)
                enddo
            enddo
        enddo
        !$omp end parallel do
        call log_shell_budget_summary(label, half, kfromto, shell_nodes, cell_count, obs_count, den_sum)
        do ik = 1,nk
            node_count = max(1, shell_nodes(ik))
            coverage = real(obs_count(ik),dp) / real(node_count,dp)
            if( obs_count(ik) > 0 )then
                den_mean = den_sum(ik) / real(obs_count(ik),dp)
            else
                den_mean = 0.d0
                den_min(ik) = 0.d0
            endif
            write(logfhandle,'(A,1X,A,2X,A,1X,A,2X,A,I4,2X,A,I0,2X,A,I0,2X,A,I0,2X,A,F8.3,2X,A,ES12.4,2X,A,ES12.4,2X,A,ES12.4,2X,A,ES12.4)') &
                &'obsfield shell-accum', trim(label), 'half=', trim(half), 'k=', kfromto(1)+ik-1, &
                &'field_nodes=', shell_nodes(ik), 'cart_cells=', cell_count(ik), 'obs_cells=', obs_count(ik), &
                &'obs_per_node=', coverage, 'den_sum=', den_sum(ik), 'den_mean=', den_mean, &
                &'den_min=', den_min(ik), 'den_max=', den_max(ik)
        enddo
        deallocate(obs_count, cell_count, den_sum, den_min, den_max)
    end subroutine log_merged_field_shell_budget

    subroutine log_shell_budget_summary( label, half, kfromto, shell_nodes, cell_count, obs_count, den_sum )
        character(len=*), intent(in) :: label, half
        integer,          intent(in) :: kfromto(2)
        integer,          intent(in) :: shell_nodes(:), cell_count(:), obs_count(:)
        real(dp),         intent(in) :: den_sum(:)
        real(dp) :: obs_per_node, obs_per_cell
        integer  :: total_nodes, total_cells, total_obs
        total_nodes = sum(shell_nodes)
        total_cells = sum(cell_count)
        total_obs   = sum(obs_count)
        obs_per_node = real(total_obs,dp) / real(max(1,total_nodes),dp)
        obs_per_cell = real(total_obs,dp) / real(max(1,total_cells),dp)
        write(logfhandle,'(A,1X,A,2X,A,1X,A,2X,A,I0,2X,A,I0,2X,A,I0,2X,A,F8.3,2X,A,F8.3,2X,A,ES12.4)') &
            &'obsfield shell-accum summary', trim(label), 'half=', trim(half), &
            &'field_nodes=', total_nodes, 'cart_cells=', total_cells, 'obs_cells=', total_obs, &
            &'obs_per_node=', obs_per_node, 'obs_per_cell=', obs_per_cell, 'den_sum=', sum(den_sum)
        call log_shell_budget_mismatch(label, half, kfromto, shell_nodes, cell_count)
    end subroutine log_shell_budget_summary

    subroutine log_shell_budget_mismatch( label, half, kfromto, shell_nodes, cell_count )
        character(len=*), intent(in) :: label, half
        integer,          intent(in) :: kfromto(2)
        integer,          intent(in) :: shell_nodes(:), cell_count(:)
        real(dp) :: rel, rel_min, rel_max, rel_sum, rel_sumsq, rel_absmax, total_rel
        real(dp) :: rel_cell_sumsq, rel_cell_rms, cell_weight
        real(dp) :: band_sumsq(3), band_rms(3)
        integer  :: i, nvalid, iband, band_count(3), band_first(3), band_last(3)
        integer  :: n, total_nodes, total_cells, shell
        n = size(shell_nodes)
        if( n < 1 .or. size(cell_count) /= n ) return
        rel_min    =  huge(0.d0)
        rel_max    = -huge(0.d0)
        rel_sum    = 0.d0
        rel_sumsq  = 0.d0
        rel_cell_sumsq = 0.d0
        rel_cell_rms   = 0.d0
        rel_absmax = 0.d0
        band_sumsq = 0.d0
        band_rms   = 0.d0
        band_count = 0
        band_first = huge(0)
        band_last  = -huge(0)
        nvalid     = 0
        do i = 1,n
            if( cell_count(i) < 1 ) cycle
            shell      = kfromto(1) + i - 1
            rel        = real(shell_nodes(i),dp) / real(cell_count(i),dp) - 1.d0
            cell_weight = real(cell_count(i),dp)
            rel_min    = min(rel_min, rel)
            rel_max    = max(rel_max, rel)
            rel_absmax = max(rel_absmax, abs(rel))
            rel_sum    = rel_sum + rel
            rel_sumsq  = rel_sumsq + rel * rel
            rel_cell_sumsq = rel_cell_sumsq + cell_weight * rel * rel
            iband = min(3, max(1, 1 + (3 * (i - 1)) / n))
            band_count(iband) = band_count(iband) + 1
            band_sumsq(iband) = band_sumsq(iband) + rel * rel
            band_first(iband) = min(band_first(iband), shell)
            band_last(iband)  = max(band_last(iband),  shell)
            nvalid = nvalid + 1
        enddo
        if( nvalid < 1 ) return
        do iband = 1,3
            if( band_count(iband) > 0 )then
                band_rms(iband) = sqrt(band_sumsq(iband) / real(band_count(iband),dp))
            else
                band_first(iband) = 0
                band_last(iband)  = 0
            endif
        enddo
        total_nodes = sum(shell_nodes)
        total_cells = sum(cell_count)
        total_rel = real(total_nodes,dp) / real(max(1,total_cells),dp) - 1.d0
        if( total_cells > 0 ) rel_cell_rms = sqrt(rel_cell_sumsq / real(total_cells,dp))
        write(logfhandle,'(A,1X,A,2X,A,1X,A,2X,A,I0,A,I0,2X,A,I0,A,I0,2X,A,I0,A,I0,2X,*(A,ES12.4,2X))') &
            &'obsfield shell-budget mismatch', trim(label), 'half=', trim(half), &
            &'low_k=', band_first(1), ':', band_last(1), &
            &'mid_k=', band_first(2), ':', band_last(2), &
            &'high_k=', band_first(3), ':', band_last(3), &
            &'total_rel=', total_rel, 'rel_mean=', rel_sum / real(nvalid,dp), &
            &'rel_rms=', sqrt(rel_sumsq / real(nvalid,dp)), 'rel_cell_rms=', rel_cell_rms, 'rel_min=', rel_min, &
            &'rel_max=', rel_max, 'rel_absmax=', rel_absmax, &
            &'low_rms=', band_rms(1), 'mid_rms=', band_rms(2), 'high_rms=', band_rms(3)
    end subroutine log_shell_budget_mismatch

    subroutine obsfield_eo_write_merged_volume( self, fname, box, smpd, kfromto, invtau2_even, invtau2_odd, &
            &prior_start, mag_correction )
        class(fgrid_obsfield_eo), intent(in) :: self
        class(string),             intent(in) :: fname
        integer,                   intent(in) :: box, kfromto(2), prior_start
        real,                      intent(in) :: smpd
        real(dp),                  intent(in) :: invtau2_even(kfromto(1):kfromto(2))
        real(dp),                  intent(in) :: invtau2_odd( kfromto(1):kfromto(2))
        real, optional,            intent(in) :: mag_correction
        type(image) :: vol
        complex(dp) :: val_even, val_odd, merged
        real(dp)    :: den_even, den_odd, denom_weight
        real        :: mag_corr
        integer     :: h, k, l, shell, phys(3), nyq
        logical     :: have_even, have_odd
        if( box < 1 ) THROW_HARD('invalid box; obsfield_eo_write_merged_volume')
        if( kfromto(1) > kfromto(2) ) THROW_HARD('invalid k-range; obsfield_eo_write_merged_volume')
        if( .not. self%even%initialized .or. .not. self%odd%initialized )then
            THROW_HARD('obsfield not initialized; obsfield_eo_write_merged_volume')
        endif
        if( .not. self%even%compatible_with(self%odd) )then
            THROW_HARD('incompatible even/odd obsfields; obsfield_eo_write_merged_volume')
        endif
        mag_corr = 1.0
        if( present(mag_correction) ) mag_corr = mag_correction
        nyq = fdim(box) - 1
        call vol%new([box,box,box], smpd, wthreads=.false.)
        call vol%set_cmat(cmplx(0.,0.))
        !$omp parallel do collapse(3) default(shared) schedule(static) proc_bind(close)&
        !$omp private(h,k,l,shell,phys,val_even,val_odd,merged,den_even,den_odd,denom_weight,have_even,have_odd)
        do h = self%even%lims(1,1), self%even%lims(1,2)
            do k = self%even%lims(2,1), self%even%lims(2,2)
                do l = self%even%lims(3,1), self%even%lims(3,2)
                    shell = nint(sqrt(real(h*h + k*k + l*l, dp)))
                    if( shell < kfromto(1) .or. shell > kfromto(2) ) cycle
                    if( shell > nyq ) cycle
                    call obsfield_restored_cell(self%even, h, k, l, kfromto, invtau2_even, prior_start, &
                        &val_even, den_even, have_even)
                    call obsfield_restored_cell(self%odd, h, k, l, kfromto, invtau2_odd, prior_start, &
                        &val_odd, den_odd, have_odd)
                    if( .not. have_even .and. .not. have_odd ) cycle
                    denom_weight = 0.d0
                    merged       = DCMPLX_ZERO
                    if( have_even )then
                        denom_weight = denom_weight + den_even
                        merged       = merged + val_even * den_even
                    endif
                    if( have_odd )then
                        denom_weight = denom_weight + den_odd
                        merged       = merged + val_odd * den_odd
                    endif
                    if( denom_weight <= DTINY ) cycle
                    phys = vol%comp_addr_phys(h,k,l)
                    call vol%set_cmat_at(phys(1), phys(2), phys(3), cmplx(merged / denom_weight))
                enddo
            enddo
        enddo
        !$omp end parallel do
        call vol%ifft
        if( mag_corr > TINY ) call vol%div(mag_corr)
        call vol%write(fname, del_if_exists=.true.)
        call vol%kill
    end subroutine obsfield_eo_write_merged_volume

    subroutine obsfield_restored_cell( self, h, k, l, kfromto, invtau2, prior_start, value, denom, have_cell )
        class(fgrid_obs_field), intent(in)  :: self
        integer,                intent(in)  :: h, k, l, kfromto(2), prior_start
        real(dp),               intent(in)  :: invtau2(kfromto(1):kfromto(2))
        complex(dp),            intent(out) :: value
        real(dp),               intent(out) :: denom
        logical,                intent(out) :: have_cell
        real(dp) :: prior
        integer  :: shell
        value     = DCMPLX_ZERO
        denom     = 0.d0
        have_cell = .false.
        if( .not. self%initialized ) return
        if( h < lbound(self%grid_num,1) .or. h > ubound(self%grid_num,1) ) return
        if( k < lbound(self%grid_num,2) .or. k > ubound(self%grid_num,2) ) return
        if( l < lbound(self%grid_num,3) .or. l > ubound(self%grid_num,3) ) return
        if( .not. self%grid_obs(h,k,l) ) return
        shell = nint(sqrt(real(h*h + k*k + l*l, dp)))
        prior = 0.d0
        if( shell >= kfromto(1) .and. shell <= kfromto(2) ) prior = invtau2(shell)
        denom = self%grid_den(h,k,l) + prior
        if( shell >= kfromto(1) .and. shell <= kfromto(2) .and. shell >= prior_start .and. prior <= DTINY )then
            denom = denom + unsampled_floor(self%grid_den(h,k,l))
        endif
        if( denom <= DTINY ) return
        value     = self%grid_num(h,k,l) / denom
        have_cell = .true.
    end subroutine obsfield_restored_cell

    integer function obsfield_get_nobs( self )
        class(fgrid_obs_field), intent(in) :: self
        obsfield_get_nobs = self%nobs
    end function obsfield_get_nobs

    integer function obsfield_get_ncells( self )
        class(fgrid_obs_field), intent(in) :: self
        obsfield_get_ncells = self%ncells
    end function obsfield_get_ncells

    logical function obsfield_compatible_with( self, other )
        class(fgrid_obs_field), intent(in) :: self, other
        obsfield_compatible_with = self%initialized .and. other%initialized
        if( .not. obsfield_compatible_with ) return
        obsfield_compatible_with = (self%pf == other%pf) .and. &
            (self%nyq == other%nyq) .and. &
            (self%iwinsz == other%iwinsz) .and. &
            (self%wdim == other%wdim) .and. &
            all(self%lims == other%lims) .and. &
            all(self%grid_lims == other%grid_lims) .and. &
            all(self%grid_shape == other%grid_shape)
    end function obsfield_compatible_with

    subroutine obsfield_write_local( self, funit, label )
        class(fgrid_obs_field), intent(in) :: self
        integer,                intent(in) :: funit
        character(len=*),       intent(in) :: label
        integer :: header(OBSFIELD_HEADER_SIZE), io_stat
        if( .not. self%initialized ) THROW_HARD('obsfield not initialized; '//trim(label))
        header(1:4)   = [self%pf, self%nyq, self%iwinsz, self%wdim]
        header(5:10)  = reshape(self%lims,      [6])
        header(11:16) = reshape(self%grid_lims, [6])
        header(17:19) = self%grid_shape
        header(20:21) = [self%nobs, OBSFIELD_GATE_FORMAT_SIGN * self%ncells]
        write(funit, iostat=io_stat) header
        call fileiochk('obsfield_write_local header; '//trim(label), io_stat)
        write(funit, iostat=io_stat) self%grid_num
        call fileiochk('obsfield_write_local numerator; '//trim(label), io_stat)
        write(funit, iostat=io_stat) self%grid_den
        call fileiochk('obsfield_write_local density; '//trim(label), io_stat)
        write(funit, iostat=io_stat) self%grid_obs
        call fileiochk('obsfield_write_local observation gate; '//trim(label), io_stat)
    end subroutine obsfield_write_local

    subroutine obsfield_read_local( self, funit, label )
        class(fgrid_obs_field), intent(inout) :: self
        integer,                intent(in)    :: funit
        character(len=*),       intent(in)    :: label
        integer :: header(OBSFIELD_HEADER_SIZE), io_stat
        integer :: lims(3,2), grid_lims(3,2), grid_shape(3), ncells
        read(funit, iostat=io_stat) header
        call fileiochk('obsfield_read_local header; '//trim(label), io_stat)
        lims       = reshape(header(5:10),  [3,2])
        grid_lims  = reshape(header(11:16), [3,2])
        grid_shape = header(17:19)
        if( header(21) >= 0 .or. mod(header(21), OBSFIELD_GATE_FORMAT_SIGN) /= 0 )then
            THROW_HARD('obsfield file lacks gate-format marker; remove stale obsfield parts: '//trim(label))
        endif
        ncells     = header(21) / OBSFIELD_GATE_FORMAT_SIGN
        call self%new(lims, header(2), zero_init=.false.)
        if( self%pf /= header(1) .or. self%iwinsz /= header(3) .or. self%wdim /= header(4) ) &
            &THROW_HARD('obsfield scalar header mismatch; '//trim(label))
        if( any(self%grid_lims /= grid_lims) .or. any(self%grid_shape /= grid_shape) .or. &
            &self%ncells /= ncells ) THROW_HARD('obsfield grid header mismatch; '//trim(label))
        read(funit, iostat=io_stat) self%grid_num
        call fileiochk('obsfield_read_local numerator; '//trim(label), io_stat)
        read(funit, iostat=io_stat) self%grid_den
        call fileiochk('obsfield_read_local density; '//trim(label), io_stat)
        read(funit, iostat=io_stat) self%grid_obs
        call fileiochk('obsfield_read_local observation gate; '//trim(label), io_stat)
        self%nobs = count(self%grid_obs)
    end subroutine obsfield_read_local

    subroutine obsfield_eo_new( self, lims, nyq )
        class(fgrid_obsfield_eo), intent(inout) :: self
        integer,                   intent(in)    :: lims(3,2), nyq
        call self%even%new(lims, nyq)
        call self%odd%new( lims, nyq)
    end subroutine obsfield_eo_new

    subroutine obsfield_eo_reset( self )
        class(fgrid_obsfield_eo), intent(inout) :: self
        call self%even%reset
        call self%odd%reset
    end subroutine obsfield_eo_reset

    subroutine obsfield_eo_kill( self )
        class(fgrid_obsfield_eo), intent(inout) :: self
        call self%even%kill
        call self%odd%kill
    end subroutine obsfield_eo_kill

    subroutine obsfield_eo_insert_plane( self, se, o, fpl, eo, pwght, shift_crop )
        class(fgrid_obsfield_eo), intent(inout) :: self
        class(sym),                intent(inout) :: se
        class(ori),                intent(inout) :: o
        class(fplane_type),        intent(in)    :: fpl
        integer,                   intent(in)    :: eo
        real,                      intent(in)    :: pwght
        real, optional,            intent(in)    :: shift_crop(2)
        select case(eo)
            case(-1,0)
                if( present(shift_crop) )then
                    call self%even%insert_plane_oversamp(se, o, fpl, pwght, shift_crop=shift_crop)
                else
                    call self%even%insert_plane_oversamp(se, o, fpl, pwght)
                endif
            case(1)
                if( present(shift_crop) )then
                    call self%odd%insert_plane_oversamp(se, o, fpl, pwght, shift_crop=shift_crop)
                else
                    call self%odd%insert_plane_oversamp(se, o, fpl, pwght)
                endif
            case DEFAULT
                THROW_HARD('unsupported eo flag; obsfield_eo_insert_plane')
        end select
    end subroutine obsfield_eo_insert_plane

    subroutine obsfield_eo_append_field( self, src )
        class(fgrid_obsfield_eo), intent(inout) :: self
        class(fgrid_obsfield_eo), intent(in)    :: src
        call self%even%append_field(src%even)
        call self%odd%append_field(src%odd)
    end subroutine obsfield_eo_append_field

    subroutine obsfield_eo_blend_field( self, prev, update_frac )
        class(fgrid_obsfield_eo), intent(inout) :: self
        class(fgrid_obsfield_eo), intent(in)    :: prev
        real(dp),                  intent(in)    :: update_frac
        call self%even%blend_field(prev%even, update_frac)
        call self%odd%blend_field( prev%odd,  update_frac)
    end subroutine obsfield_eo_blend_field

    subroutine obsfield_eo_write( self, fname )
        class(fgrid_obsfield_eo), intent(in) :: self
        class(string),             intent(in) :: fname
        integer :: funit, io_stat
        call fopen(funit, fname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk('obsfield_eo_write open; '//fname%to_char(), io_stat)
        call obsfield_write_local(self%even, funit, 'even '//fname%to_char())
        call obsfield_write_local(self%odd,  funit, 'odd '//fname%to_char())
        call fclose(funit)
    end subroutine obsfield_eo_write

    subroutine obsfield_eo_read( self, fname )
        class(fgrid_obsfield_eo), intent(inout) :: self
        class(string),             intent(in)    :: fname
        integer :: funit, io_stat
        if( .not. file_exists(fname) ) THROW_HARD('obsfield file does not exist: '//fname%to_char())
        call fopen(funit, fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
        call fileiochk('obsfield_eo_read open; '//fname%to_char(), io_stat)
        call obsfield_read_local(self%even, funit, 'even '//fname%to_char())
        call obsfield_read_local(self%odd,  funit, 'odd '//fname%to_char())
        call fclose(funit)
    end subroutine obsfield_eo_read

end module simple_fgrid_obsfield
