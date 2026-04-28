!@descr: dense Fourier-grid observation field for restoration experiments
module simple_fgrid_obsfield
use simple_core_module_api
implicit none

public :: fgrid_obs_field, fgrid_obs_field_eo
private
#include "simple_local_flags.inc"

! Two h-colors make nearest-cell writes race-free without atomics: all samples
! for one h are handled by one thread, and h values in the same color are at
! least two native grid units apart before and after rotation.
integer, parameter :: OBSFIELD_NN_OMP_STRIDE = 2

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
    logical,     allocatable :: grid_assigned(:,:,:)
  contains
    procedure, public :: new               => obsfield_new
    procedure, public :: reset             => obsfield_reset
    procedure, public :: kill              => obsfield_kill
    procedure, public :: insert_plane_oversamp
    procedure, public :: append_field      => obsfield_append_field
    procedure, public :: extract_polar     => obsfield_extract_polar
    procedure, public :: get_nobs          => obsfield_get_nobs
    procedure, public :: get_ncells        => obsfield_get_ncells
    procedure, private :: compatible_with  => obsfield_compatible_with
end type fgrid_obs_field

! Even/odd wrapper matching the reconstruction convention:
! eo=-1 or eo=0 contributes to even, eo=1 contributes to odd.
type :: fgrid_obs_field_eo
    type(fgrid_obs_field) :: even
    type(fgrid_obs_field) :: odd
  contains
    procedure, public :: new           => obsfield_eo_new
    procedure, public :: reset         => obsfield_eo_reset
    procedure, public :: kill          => obsfield_eo_kill
    procedure, public :: insert_plane  => obsfield_eo_insert_plane
    procedure, public :: append_field  => obsfield_eo_append_field
end type fgrid_obs_field_eo

contains

    subroutine obsfield_new( self, lims, nyq )
        class(fgrid_obs_field), intent(inout) :: self
        integer,                intent(in)    :: lims(3,2)
        integer,                intent(in)    :: nyq
        real(dp) :: ncells_dp
        integer  :: dim
        call self%kill
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
        allocate(self%grid_num(self%grid_lims(1,1):self%grid_lims(1,2), &
            self%grid_lims(2,1):self%grid_lims(2,2), self%grid_lims(3,1):self%grid_lims(3,2)), &
            source=DCMPLX_ZERO)
        allocate(self%grid_den(self%grid_lims(1,1):self%grid_lims(1,2), &
            self%grid_lims(2,1):self%grid_lims(2,2), self%grid_lims(3,1):self%grid_lims(3,2)), &
            source=0.d0)
        allocate(self%grid_assigned(self%grid_lims(1,1):self%grid_lims(1,2), &
            self%grid_lims(2,1):self%grid_lims(2,2), self%grid_lims(3,1):self%grid_lims(3,2)), &
            source=.false.)
        self%initialized = .true.
    end subroutine obsfield_new

    subroutine obsfield_reset( self )
        class(fgrid_obs_field), intent(inout) :: self
        self%nobs      = 0
        if( allocated(self%grid_num) ) self%grid_num = DCMPLX_ZERO
        if( allocated(self%grid_den) ) self%grid_den = 0.d0
        if( allocated(self%grid_assigned) ) self%grid_assigned = .false.
    end subroutine obsfield_reset

    subroutine obsfield_kill( self )
        class(fgrid_obs_field), intent(inout) :: self
        if( allocated(self%grid_num) ) deallocate(self%grid_num)
        if( allocated(self%grid_den) ) deallocate(self%grid_den)
        if( allocated(self%grid_assigned) ) deallocate(self%grid_assigned)
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
    subroutine insert_plane_oversamp( self, se, o, fpl, pwght )
        use simple_math,    only: ceil_div, floor_div
        use simple_math_ft, only: fplane_get_cmplx, fplane_get_ctfsq
        class(fgrid_obs_field), intent(inout) :: self
        class(sym),             intent(inout) :: se
        class(ori),             intent(inout) :: o
        class(fplane_type),     intent(in)    :: fpl
        real,                   intent(in)    :: pwght
        type(ori)   :: o_sym
        complex(sp) :: cmplx_raw
        complex(dp) :: comp
        real(dp)    :: ctfval, pwght_dp, pwght_pf2_dp, cell_w, cell_w_norm
        real(sp)    :: ctfsq_raw
        real(sp)    :: loc(3), delta(3), R(3,3)
        real        :: rotmats(se%get_nsym(),3,3)
        integer     :: fpllims_pd(3,2), fpllims(3,2), coord(3)
        integer     :: nsym, isym, h, k, hp, kp, pf_local, l
        integer     :: nyq_disk, h_sq, k_max_h, k_lo, k_hi
        integer     :: lim1_lo, gl_lo1, gl_lo2, gl_lo3, gl_hi1, gl_hi2, gl_hi3
        integer     :: nobs_add
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
        ! integer disk gate: bit-exact equivalent of original nint(sqrt(h^2+k^2)) > nyq
        ! since for non-negative integer n,  n > nyq*(nyq+1)  <=>  nint(sqrt(n)) > nyq
        nyq_disk     = self%nyq * (self%nyq + 1)
        lim1_lo      = self%lims(1,1)
        gl_lo1 = self%grid_lims(1,1); gl_hi1 = self%grid_lims(1,2)
        gl_lo2 = self%grid_lims(2,1); gl_hi2 = self%grid_lims(2,2)
        gl_lo3 = self%grid_lims(3,1); gl_hi3 = self%grid_lims(3,2)
        !$omp parallel default(shared) private(isym,R,l,h,k,h_sq,k_max_h,k_lo,k_hi,&
        !$omp& cmplx_raw,ctfsq_raw,comp,ctfval,loc,delta,coord,hp,kp,cell_w)&
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
                    do k = k_lo, k_hi
                        kp = k * pf_local
                        ! raw padded-plane samples (single precision); cheap zero-skip
                        cmplx_raw = fplane_get_cmplx(fpl, hp, kp)
                        ctfsq_raw = fplane_get_ctfsq(fpl, hp, kp)
                        if( abs(real(cmplx_raw)) + abs(aimag(cmplx_raw)) <= TINY .and. &
                            ctfsq_raw <= TINY ) cycle
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
                        self%grid_assigned(coord(1),coord(2),coord(3)) = .true.
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
        if( .not. src%initialized ) return
        if( .not. self%initialized ) THROW_HARD('destination not initialized; obsfield_append_field')
        if( .not. self%compatible_with(src) ) THROW_HARD('incompatible observation fields; obsfield_append_field')
        self%grid_num = self%grid_num + src%grid_num
        self%grid_den = self%grid_den + src%grid_den
        self%grid_assigned = self%grid_assigned .or. src%grid_assigned
        self%nobs      = self%nobs + src%nobs
    end subroutine obsfield_append_field

    ! Gather central sections from the dense expanded grid. This keeps the
    ! polar=no volume-sampling shape: direct array addressing in the 3D KB stencil.
    ! Friedel recovery: stencil cells with c1 < lims(1,1) (which were never
    ! deposited on insertion) are mapped to (-c1,-c2,-c3); grid_num is
    ! conjugated, grid_den is read as-is. The decision depends only on c1,
    ! so it is hoisted to the l-loop. c2/c3 grid-bound checks depend only on
    ! m,n, so they are hoisted to those loops as cheap defensive guards.
    subroutine obsfield_extract_polar( self, eulspace, nrefs, kfromto, polar_x, polar_y, pfts, ctf2 )
        class(fgrid_obs_field), intent(in)    :: self
        class(oris),            intent(inout) :: eulspace
        integer,                intent(in)    :: nrefs, kfromto(2)
        real(sp),               intent(in)    :: polar_x(:,:), polar_y(:,:)
        complex(dp),            intent(inout) :: pfts(:,:,:)
        real(dp),               intent(inout) :: ctf2(:,:,:)
        complex(dp) :: acc_num, cell_num
        real(dp)    :: acc_den, wd
        real(sp)    :: R(3,3), loc(3), px, py
        real(sp), allocatable :: w(:,:,:)
        integer     :: iproj, irot, k, kloc, pftsz, kfromto1
        integer     :: win1_1, win1_2, win1_3, c1, c2, c3, a1, a2, a3
        integer     :: l, m, n, iwinsz_l, wdim_l, lim1_lo
        integer     :: gl_lo1, gl_lo2, gl_lo3, gl_hi1, gl_hi2, gl_hi3
        logical     :: l_conj
        ! if( .not. self%initialized ) THROW_HARD('obsfield not initialized; extract_polar')
        ! if( nrefs < 1 ) THROW_HARD('invalid nrefs; extract_polar')
        ! if( kfromto(1) < 1 .or. kfromto(1) > kfromto(2) ) THROW_HARD('invalid kfromto; extract_polar')
        pftsz = size(polar_x,1)
        ! if( size(polar_y,1) /= pftsz .or. size(polar_y,2) /= size(polar_x,2) )&
        !     &THROW_HARD('polar coordinate shape mismatch; extract_polar')
        ! if( size(polar_x,2) /= kfromto(2)-kfromto(1)+1 ) THROW_HARD('polar k-span mismatch; extract_polar')
        ! if( size(pfts,1) /= pftsz .or. size(pfts,2) /= size(polar_x,2) .or. size(pfts,3) < nrefs )&
        !     &THROW_HARD('pfts shape mismatch; extract_polar')
        ! if( size(ctf2,1) /= pftsz .or. size(ctf2,2) /= size(polar_x,2) .or. size(ctf2,3) < nrefs )&
        !     &THROW_HARD('ctf2 shape mismatch; extract_polar')
        pfts(:,:,1:nrefs) = DCMPLX_ZERO
        ctf2(:,:,1:nrefs) = 0.d0
        if( self%nobs < 1 ) return
        ! hoisted scalars
        iwinsz_l = self%iwinsz
        wdim_l   = self%wdim
        lim1_lo  = self%lims(1,1)
        kfromto1 = kfromto(1)
        gl_lo1   = self%grid_lims(1,1); gl_hi1 = self%grid_lims(1,2)
        gl_lo2   = self%grid_lims(2,1); gl_hi2 = self%grid_lims(2,2)
        gl_lo3   = self%grid_lims(3,1); gl_hi3 = self%grid_lims(3,2)
        !$omp parallel default(shared) private(iproj,R,k,kloc,irot,px,py,loc,w,acc_num,acc_den,wd,&
        !$omp& win1_1,win1_2,win1_3,c1,c2,c3,a1,a2,a3,l,m,n,l_conj,cell_num) proc_bind(close)
        allocate(w(wdim_l,wdim_l,wdim_l))
        !$omp do schedule(static)
        do iproj = 1, nrefs
            R = eulspace%get_mat(iproj)
            do k = kfromto1, kfromto(2)
                kloc = k - kfromto1 + 1
                do irot = 1, pftsz
                    px     = polar_x(irot,kloc)
                    py     = polar_y(irot,kloc)
                    ! The observation field is indexed by native grid cells.
                    ! Evaluate the KB gather in the same native coordinate
                    ! system; using padded coordinates here narrows and shifts
                    ! the extraction stencil relative to the deposited cells.
                    loc(1) = px*R(1,1) + py*R(2,1)
                    loc(2) = px*R(1,2) + py*R(2,2)
                    loc(3) = px*R(1,3) + py*R(2,3)
                    win1_1 = nint(loc(1)) - iwinsz_l
                    win1_2 = nint(loc(2)) - iwinsz_l
                    win1_3 = nint(loc(3)) - iwinsz_l
                    call self%kb%apod_mat_3d_fast(loc, iwinsz_l, wdim_l, w)
                    acc_num = DCMPLX_ZERO
                    acc_den = 0.d0
                    do n = 1, wdim_l
                        c3 = win1_3 + n - 1
                        ! c3 bound is m,n-invariant within (l,m,n)-triple; hoisted to n-loop
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
                                    if( .not. self%grid_assigned(a1,a2,a3) ) cycle
                                    cell_num = conjg(self%grid_num(a1,a2,a3))
                                    acc_num  = acc_num + wd * cell_num
                                    acc_den  = acc_den + wd * self%grid_den(a1,a2,a3)
                                else
                                    if( c1 > gl_hi1 ) cycle
                                    if( .not. self%grid_assigned(c1,c2,c3) ) cycle
                                    acc_num = acc_num + wd * self%grid_num(c1,c2,c3)
                                    acc_den = acc_den + wd * self%grid_den(c1,c2,c3)
                                endif
                            enddo
                        enddo
                    enddo
                    pfts(irot,kloc,iproj) = acc_num
                    ctf2(irot,kloc,iproj) = acc_den
                enddo
            enddo
        enddo
        !$omp end do
        deallocate(w)
        !$omp end parallel
    end subroutine obsfield_extract_polar

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

    subroutine obsfield_eo_new( self, lims, nyq )
        class(fgrid_obs_field_eo), intent(inout) :: self
        integer,                   intent(in)    :: lims(3,2), nyq
        call self%even%new(lims, nyq)
        call self%odd%new( lims, nyq)
    end subroutine obsfield_eo_new

    subroutine obsfield_eo_reset( self )
        class(fgrid_obs_field_eo), intent(inout) :: self
        call self%even%reset
        call self%odd%reset
    end subroutine obsfield_eo_reset

    subroutine obsfield_eo_kill( self )
        class(fgrid_obs_field_eo), intent(inout) :: self
        call self%even%kill
        call self%odd%kill
    end subroutine obsfield_eo_kill

    subroutine obsfield_eo_insert_plane( self, se, o, fpl, eo, pwght )
        class(fgrid_obs_field_eo), intent(inout) :: self
        class(sym),                intent(inout) :: se
        class(ori),                intent(inout) :: o
        class(fplane_type),        intent(in)    :: fpl
        integer,                   intent(in)    :: eo
        real,                      intent(in)    :: pwght
        select case(eo)
            case(-1,0)
                call self%even%insert_plane_oversamp(se, o, fpl, pwght)
            case(1)
                call self%odd%insert_plane_oversamp(se, o, fpl, pwght)
            case DEFAULT
                THROW_HARD('unsupported eo flag; obsfield_eo_insert_plane')
        end select
    end subroutine obsfield_eo_insert_plane

    subroutine obsfield_eo_append_field( self, src )
        class(fgrid_obs_field_eo), intent(inout) :: self
        class(fgrid_obs_field_eo), intent(in)    :: src
        call self%even%append_field(src%even)
        call self%odd%append_field(src%odd)
    end subroutine obsfield_eo_append_field

end module simple_fgrid_obsfield
