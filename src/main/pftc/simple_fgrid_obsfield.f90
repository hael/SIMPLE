!@descr: dense Fourier-grid observation field for restoration experiments
module simple_fgrid_obsfield
use simple_core_module_api
implicit none

public :: fgrid_obs_field, fgrid_obs_field_eo
private
#include "simple_local_flags.inc"

! Part-local Fourier-grid observation field. This is intentionally volume-like:
! particle Fourier components are accumulated into dense expanded-grid
! numerator/density arrays, and requested polar central sections are gathered
! directly from those arrays. The experimental insertion step uses pure
! nearest-cell assignment to avoid both insertion-side interpolation and the
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
    integer                :: nrejected     = 0
    logical                :: initialized   = .false.
    complex(dp), allocatable :: grid_num(:,:,:)
    real(dp),    allocatable :: grid_den(:,:,:)
  contains
    procedure, public :: new               => obsfield_new
    procedure, public :: reset             => obsfield_reset
    procedure, public :: kill              => obsfield_kill
    procedure, public :: insert_plane_oversamp
    procedure, public :: append_field      => obsfield_append_field
    procedure, public :: extract_polar     => obsfield_extract_polar
    procedure, public :: get_nobs          => obsfield_get_nobs
    procedure, public :: get_ncells        => obsfield_get_ncells
    procedure, public :: get_nrejected     => obsfield_get_nrejected
    procedure, private :: compatible_with  => obsfield_compatible_with
    procedure, private :: apod_mat_3d_fast => obsfield_apod_mat_3d_fast
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
        self%initialized = .true.
    end subroutine obsfield_new

    subroutine obsfield_reset( self )
        class(fgrid_obs_field), intent(inout) :: self
        self%nobs      = 0
        self%nrejected = 0
        if( allocated(self%grid_num) ) self%grid_num = DCMPLX_ZERO
        if( allocated(self%grid_den) ) self%grid_den = 0.d0
    end subroutine obsfield_reset

    subroutine obsfield_kill( self )
        class(fgrid_obs_field), intent(inout) :: self
        if( allocated(self%grid_num) ) deallocate(self%grid_num)
        if( allocated(self%grid_den) ) deallocate(self%grid_den)
        self%pf          = OSMPL_PAD_FAC
        self%nyq         = 0
        self%iwinsz      = 0
        self%wdim        = 0
        self%lims        = 0
        self%grid_lims   = 0
        self%grid_shape  = 0
        self%nobs        = 0
        self%ncells      = 0
        self%nrejected   = 0
        self%initialized = .false.
    end subroutine obsfield_kill

    ! Insert one particle Fourier plane by nearest-cell assignment. Coordinates
    ! are computed exactly as in reconstructor%insert_plane_oversamp, but each
    ! experimental component updates only the closest native 3D grid cell.
    subroutine insert_plane_oversamp( self, se, o, fpl, pwght )
        use simple_math,    only: ceil_div, floor_div
        use simple_math_ft, only: fplane_get_cmplx, fplane_get_ctfsq
        class(fgrid_obs_field), intent(inout) :: self
        class(sym),             intent(inout) :: se
        class(ori),             intent(inout) :: o
        class(fplane_type),     intent(in)    :: fpl
        real,                   intent(in)    :: pwght
        type(ori) :: o_sym
        complex(dp) :: comp
        real(dp)    :: ctfval
        real(sp)    :: loc(3)
        real        :: rotmats(se%get_nsym(),3,3)
        integer     :: fpllims_pd(3,2), fpllims(3,2), coord(3), addr(3)
        integer     :: nsym, isym, h, k, hp, kp, sh, l, stride
        integer     :: nobs_add, nrejected_add
        if( .not. self%initialized ) THROW_HARD('obsfield not initialized; insert_plane_oversamp')
        if( pwght < TINY ) return
        nobs_add      = 0
        nrejected_add = 0
        stride        = self%wdim
        nsym = se%get_nsym()
        rotmats(1,:,:) = o%get_mat()
        if( nsym > 1 )then
            do isym = 2, nsym
                call se%apply(o, isym, o_sym)
                rotmats(isym,:,:) = o_sym%get_mat()
            enddo
        endif
        rotmats = KBALPHA * rotmats
        fpllims_pd = fpl%frlims
        fpllims    = fpllims_pd
        fpllims(1,1) = ceil_div (fpllims_pd(1,1), self%pf)
        fpllims(1,2) = floor_div(fpllims_pd(1,2), self%pf)
        fpllims(2,1) = ceil_div (fpllims_pd(2,1), self%pf)
        fpllims(2,2) = floor_div(fpllims_pd(2,2), self%pf)
        !$omp parallel default(shared) private(isym,h,k,l,sh,comp,ctfval,loc,coord,addr,hp,kp)&
        !$omp& reduction(+:nobs_add,nrejected_add) proc_bind(close)
        do isym = 1, nsym
            do l = 0, stride-1
                !$omp do schedule(static)
                do h = fpllims(1,1)+l, fpllims(1,2), stride
                    hp = h * self%pf
                    do k = fpllims(2,1), fpllims(2,2)
                        kp = k * self%pf
                        sh = nint(sqrt(real(h*h + k*k)))
                        if( sh > self%nyq ) cycle
                        loc   = matmul(real([h,k,0],sp), real(rotmats(isym,:,:),sp))
                        coord = nint(loc / real(self%pf,sp))
                        comp   = real(pwght,dp) * real(self%pf*self%pf,dp) * &
                            cmplx(fplane_get_cmplx(fpl, hp, kp),kind=dp)
                        ctfval = real(pwght,dp) * real(fplane_get_ctfsq(fpl, hp, kp),dp)
                        if( abs(comp) <= DTINY .and. ctfval <= DTINY ) cycle
                        addr = coord
                        ! Match reconstructor insertion: skip the fully redundant
                        ! negative first-axis Friedel mate and recover it on readout.
                        if( addr(1) < self%lims(1,1) ) cycle
                        if( any(addr < self%grid_lims(:,1)) .or. any(addr > self%grid_lims(:,2)) )then
                            nrejected_add = nrejected_add + 1
                            cycle
                        endif
                        self%grid_num(addr(1),addr(2),addr(3)) = self%grid_num(addr(1),addr(2),addr(3)) + comp
                        self%grid_den(addr(1),addr(2),addr(3)) = self%grid_den(addr(1),addr(2),addr(3)) + ctfval
                        nobs_add = nobs_add + 1
                    enddo
                enddo
                !$omp end do
            enddo
        enddo
        !$omp end parallel
        self%nobs      = self%nobs      + nobs_add
        self%nrejected = self%nrejected + nrejected_add
        call o_sym%kill
    end subroutine insert_plane_oversamp

    subroutine obsfield_append_field( self, src )
        class(fgrid_obs_field), intent(inout) :: self
        class(fgrid_obs_field), intent(in)    :: src
        if( .not. src%initialized ) return
        if( .not. self%initialized ) THROW_HARD('destination not initialized; obsfield_append_field')
        if( .not. self%compatible_with(src) ) THROW_HARD('incompatible observation fields; obsfield_append_field')
        self%grid_num = self%grid_num + src%grid_num
        self%grid_den = self%grid_den + src%grid_den
        self%nobs      = self%nobs + src%nobs
        self%nrejected = self%nrejected + src%nrejected
    end subroutine obsfield_append_field

    ! Gather central sections from the dense expanded grid. This keeps the
    ! polar=no volume-sampling shape: direct array addressing in the 3D KB stencil
    subroutine obsfield_extract_polar( self, eulspace, nrefs, kfromto, polar_x, polar_y, pfts, ctf2 )
        class(fgrid_obs_field), intent(in)    :: self
        class(oris),            intent(inout) :: eulspace
        integer,                intent(in)    :: nrefs, kfromto(2)
        real(sp),               intent(in)    :: polar_x(:,:), polar_y(:,:)
        complex(dp),            intent(inout) :: pfts(:,:,:)
        real(dp),               intent(inout) :: ctf2(:,:,:)
        complex(dp) :: acc_num, cell_num
        real(dp)    :: acc_den
        real(sp)    :: R(3,3), loc(3), px, py, w(self%wdim,self%wdim,self%wdim)
        integer     :: iproj, irot, k, kloc, pftsz, win(2,3), coord(3), addr(3)
        integer     :: l, m, n
        logical     :: l_conj
        if( .not. self%initialized ) THROW_HARD('obsfield not initialized; extract_polar')
        if( nrefs < 1 ) THROW_HARD('invalid nrefs; extract_polar')
        if( kfromto(1) < 1 .or. kfromto(1) > kfromto(2) ) THROW_HARD('invalid kfromto; extract_polar')
        pftsz = size(polar_x,1)
        if( size(polar_y,1) /= pftsz .or. size(polar_y,2) /= size(polar_x,2) )&
            &THROW_HARD('polar coordinate shape mismatch; extract_polar')
        if( size(polar_x,2) /= kfromto(2)-kfromto(1)+1 ) THROW_HARD('polar k-span mismatch; extract_polar')
        if( size(pfts,1) /= pftsz .or. size(pfts,2) /= size(polar_x,2) .or. size(pfts,3) < nrefs )&
            &THROW_HARD('pfts shape mismatch; extract_polar')
        if( size(ctf2,1) /= pftsz .or. size(ctf2,2) /= size(polar_x,2) .or. size(ctf2,3) < nrefs )&
            &THROW_HARD('ctf2 shape mismatch; extract_polar')
        pfts(:,:,1:nrefs) = DCMPLX_ZERO
        ctf2(:,:,1:nrefs) = 0.d0
        if( self%nobs < 1 ) return
        !$omp parallel do default(shared) private(iproj,R,k,kloc,irot,px,py,loc,win,w,acc_num,acc_den,coord)&
        !$omp& private(l,m,n,addr,l_conj,cell_num) schedule(static) proc_bind(close)
        do iproj = 1, nrefs
            R = eulspace%get_mat(iproj)
            do k = kfromto(1), kfromto(2)
                kloc = k - kfromto(1) + 1
                do irot = 1, pftsz
                    px = polar_x(irot,kloc)
                    py = polar_y(irot,kloc)
                    loc(1) = KBALPHA * (px*R(1,1) + py*R(2,1))
                    loc(2) = KBALPHA * (px*R(1,2) + py*R(2,2))
                    loc(3) = KBALPHA * (px*R(1,3) + py*R(2,3))
                    win(1,:) = nint(loc / real(self%pf,sp)) - self%iwinsz
                    win(2,:) = win(1,:) + self%wdim - 1
                    call self%apod_mat_3d_fast(loc, w)
                    acc_num = DCMPLX_ZERO
                    acc_den = 0.d0
                    do n = 1, self%wdim
                        coord(3) = win(1,3) + n - 1
                        do m = 1, self%wdim
                            coord(2) = win(1,2) + m - 1
                            do l = 1, self%wdim
                                coord(1) = win(1,1) + l - 1
                                addr   = coord
                                l_conj = .false.
                                if( addr(1) < self%lims(1,1) )then
                                    addr   = -addr
                                    l_conj = .true.
                                endif
                                if( any(addr < self%grid_lims(:,1)) .or. any(addr > self%grid_lims(:,2)) ) cycle
                                cell_num = self%grid_num(addr(1),addr(2),addr(3))
                                if( l_conj ) cell_num = conjg(cell_num)
                                acc_num = acc_num + real(w(l,m,n),dp) * cell_num
                                acc_den = acc_den + real(w(l,m,n),dp) * self%grid_den(addr(1),addr(2),addr(3))
                            enddo
                        enddo
                    enddo
                    pfts(irot,kloc,iproj) = acc_num
                    ctf2(irot,kloc,iproj) = acc_den
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine obsfield_extract_polar

    integer function obsfield_get_nobs( self )
        class(fgrid_obs_field), intent(in) :: self
        obsfield_get_nobs = self%nobs
    end function obsfield_get_nobs

    integer function obsfield_get_ncells( self )
        class(fgrid_obs_field), intent(in) :: self
        obsfield_get_ncells = self%ncells
    end function obsfield_get_ncells

    integer function obsfield_get_nrejected( self )
        class(fgrid_obs_field), intent(in) :: self
        obsfield_get_nrejected = self%nrejected
    end function obsfield_get_nrejected

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

    ! Same normalized separable 3D KB stencil as kbinterpol%apod_mat_3d, but
    ! using apod_fast to avoid repeated Bessel evaluations during extraction.
    subroutine obsfield_apod_mat_3d_fast( self, loc, kbw )
        class(fgrid_obs_field), intent(in)  :: self
        real(sp),               intent(in)  :: loc(3)
        real(sp),               intent(out) :: kbw(self%wdim,self%wdim,self%wdim)
        integer  :: win_lo(3)
        real(sp) :: base(3), ww(3)
        real(sp) :: wx(self%wdim), wy(self%wdim), wz(self%wdim)
        real(sp) :: recip
        integer  :: i, j
        win_lo = nint(loc) - self%iwinsz
        base   = real(win_lo, sp) - loc
        do i = 1, self%wdim
            ww    = self%kb%apod_fast(base + real(i-1,sp))
            wx(i) = ww(1)
            wy(i) = ww(2)
            wz(i) = ww(3)
        enddo
        wx = wx * (1.0_sp / sum(wx))
        wy = wy * (1.0_sp / sum(wy))
        wz = wz * (1.0_sp / sum(wz))
        do j = 1, self%wdim
            do i = 1, self%wdim
                kbw(:,i,j) = wx(:) * (wy(i) * wz(j))
            enddo
        enddo
        recip = 1.0_sp / sum(kbw)
        kbw   = kbw * recip
    end subroutine obsfield_apod_mat_3d_fast

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
