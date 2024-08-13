! rotational origin shift alignment of band-pass limited polar projections in the Fourier domain, gradient based minimizer
module simple_pftcc_shsrch_fm
include 'simple_lib.f08'
use simple_polarft_corrcalc, only: pftcc_glob
use simple_parameters,       only: params_glob
implicit none

public :: pftcc_shsrch_fm
private
#include "simple_local_flags.inc"

type :: pftcc_shsrch_fm
    private
    real, allocatable :: scores(:)              !< Objective function
    real, allocatable :: coords(:)              !< Offset in pixels
    real, allocatable :: grid1(:,:), grid2(:,:) !< 2D grids
    real              :: trslim  = 0.           !< Search limit haf-width (pixels)
    real              :: trsincr = 0.           !< Grid width (pixels)
    integer           :: ref=0, ptcl=0          !< Reference & particle indices
    integer           :: hn      = 0            !< Search limit halfwidth (in units of trsincr)
    integer           :: nrots   = 0            !< # rotations
    integer           :: pftsz   = 0            !< nrots/2
    logical           :: opt_angle = .true.
contains
    procedure          :: new
    procedure          :: minimize
    procedure          :: exhaustive_search
    procedure, private :: interpolate_peak
    procedure          :: kill

end type pftcc_shsrch_fm

contains

    subroutine new( self, trslim, trsstep, opt_angle )
        class(pftcc_shsrch_fm), intent(inout) :: self     !< instance
        real,                   intent(in)    :: trslim   !< limits for barrier constraint
        real,                   intent(in)    :: trsstep  !<
        logical,      optional, intent(in)    :: opt_angle
        integer :: i
        call self%kill
        self%trslim  = trslim
        self%hn      = ceiling(self%trslim/trsstep)
        self%trsincr = self%trslim / real(self%hn)
        self%nrots = pftcc_glob%get_nrots()
        self%pftsz = pftcc_glob%get_pftsz()
        self%opt_angle = .true.
        if( present(opt_angle) ) self%opt_angle = opt_angle
        allocate(self%grid1(-self%hn:self%hn,-self%hn:self%hn), self%grid2(-self%hn:self%hn,-self%hn:self%hn),&
            &self%coords(-self%hn:self%hn), self%scores(self%nrots),source=0.)
        self%coords = (/(real(i)*self%trsincr,i=-self%hn,self%hn)/)
    end subroutine new

    !> minimisation routine based on identification of in-plane rotation via shift invariant metric
    subroutine minimize( self, iref, iptcl, found, irot, score, offset )
        class(pftcc_shsrch_fm), intent(inout) :: self
        integer,                intent(in)    :: iref, iptcl
        logical,                intent(out)   :: found
        integer,                intent(inout) :: irot
        real,                   intent(out)   :: score
        real,                   intent(out)   :: offset(2)
        real     :: rotmat(2,2), shift1(2), shift2(2), score1, score2
        integer  :: irot1, irot2
        self%ref  = iref
        self%ptcl = iptcl
        if( irot > self%pftsz )then
            irot1 = irot - self%pftsz
            irot2 = irot
        else
            irot1 = irot
            irot2 = irot + self%pftsz
        endif
        ! Lower bound
        call pftcc_glob%gencorrs(self%ref, self%ptcl, self%scores, kweight=.true.)
        irot   = maxloc(self%scores, dim=1)
        score  = self%scores(irot)
        offset = 0.
        found = .false.
        ! Grid search of both orientations
        call pftcc_glob%bidirectional_shift_search(self%ref, self%ptcl, irot1, self%hn, self%coords, self%grid1, self%grid2)
        ! Maxima
        call self%interpolate_peak(self%grid1, irot1, shift1, score1)
        call self%interpolate_peak(self%grid2, irot2, shift2, score2)
        if( score1 > score )then
            found  = .true.
            irot   = irot1
            score  = score1
            offset = shift1
        endif
        if( score2 > score )then
            found  = .true.
            irot   = irot2
            score  = score2
            offset = shift2
        endif
        ! Particle shift
        call rotmat2d(pftcc_glob%get_rot(irot), rotmat)
        offset = matmul(offset, rotmat)
    end subroutine minimize

    ! exhaustive coarse grid search fllowed by peak interpolation
    subroutine exhaustive_search( self, iref, iptcl, irot, score, offset )
        class(pftcc_shsrch_fm), intent(inout) :: self
        integer,                intent(in)    :: iref, iptcl
        integer,                intent(inout) :: irot
        real,                   intent(out)   :: score
        real,                   intent(out)   :: offset(2)
        real     :: rotmat(2,2), shift(2)
        integer  :: i,j,loc,ii,jj
        self%ref  = iref
        self%ptcl = iptcl
        offset    = 0.
        score     = -1.
        ii        = 0
        jj        = 0
        if( self%opt_angle )then
            ! 3D search
            self%grid1 = -1.
            irot       = 0
            do i = -self%hn,self%hn
                do j = -self%hn,self%hn
                    shift = [self%coords(i), self%coords(j)]
                    call pftcc_glob%gencorrs(self%ref, self%ptcl, shift, self%scores, kweight=params_glob%l_kweight_rot)
                    loc = maxloc(self%scores, dim=1)
                    if( self%scores(loc) > score )then
                        score  = self%scores(loc)
                        irot   = loc
                        ii     = i
                        jj     = j
                    endif
                enddo
            enddo
            do i = max(-self%hn,ii-1), min(self%hn,ii+1)
                do j = max(-self%hn,jj-1), min(self%hn,jj+1)
                    shift = [self%coords(i), self%coords(j)]
                    self%grid1(i,j) = real(pftcc_glob%gencorr_for_rot_8(self%ref, self%ptcl, real(shift,dp), irot))
                enddo
            enddo
        else
            ! 2D search, fixed rotation
            if( (irot<1) .or. (irot>self%nrots) )then
                THROW_HARD('Invalid totation index '//int2str(irot))
            endif
            do i = -self%hn,self%hn
                do j = -self%hn,self%hn
                    shift = [self%coords(i), self%coords(j)]
                    self%grid1(i,j) = real(pftcc_glob%gencorr_for_rot_8(self%ref, self%ptcl, real(shift,dp), irot))
                enddo
            enddo
        endif
        ! peak interpolation
        call self%interpolate_peak(self%grid1, irot, offset, score)
        ! Particle shift
        call rotmat2d(pftcc_glob%get_rot(irot), rotmat)
        offset = matmul(offset, rotmat)
    end subroutine exhaustive_search

    subroutine interpolate_peak( self, grid, irot, offset, score )
        class(pftcc_shsrch_fm), intent(in) :: self
        real,    intent(in)  :: grid(-self%hn:self%hn,-self%hn:self%hn)
        integer, intent(in)  :: irot
        real,    intent(out) :: offset(2), score
        real    :: alpha, beta, gamma, denom
        integer :: maxpos(2),i,j
        maxpos = maxloc(grid)-self%hn-1
        i      = maxpos(1)
        j      = maxpos(2)
        beta   = grid(i,j)
        offset = [self%coords(i), self%coords(j)]
        score  = beta
        if( abs(i) /= self%hn )then
            alpha = grid(i-1,j)
            gamma = grid(i+1,j)
            if( alpha<beta .and. gamma<beta )then
                denom = alpha + gamma - 2.*beta
                if( abs(denom) > TINY ) offset(1) = offset(1) + self%trsincr * 0.5 * (alpha-gamma) / denom
            endif
        endif
        if( abs(j) /= self%hn )then
            alpha = grid(i,j-1)
            gamma = grid(i,j+1)
            if( alpha<beta .and. gamma<beta )then
                denom = alpha + gamma - 2.*beta
                if( abs(denom) > TINY ) offset(2) = offset(2) + self%trsincr * 0.5 * (alpha-gamma) / denom
            endif
        endif
        score = pftcc_glob%gencorr_for_rot_8(self%ref, self%ptcl, real(offset,dp), irot)
        if( score < beta )then
            ! fallback
            offset = [self%coords(i), self%coords(j)]
            score  = beta
        endif
    end subroutine interpolate_peak

    subroutine kill( self )
        class(pftcc_shsrch_fm), intent(inout) :: self
        if( allocated(self%scores) )then
            deallocate(self%scores,self%coords,self%grid1,self%grid2)
        endif
        self%trslim  = 0.
        self%trsincr = 0.
        self%ref     = 0
        self%ptcl    = 0
        self%hn      = 0
        self%nrots   = 0
        self%pftsz   = 0
        self%opt_angle = .true.
    end subroutine kill

end module simple_pftcc_shsrch_fm
