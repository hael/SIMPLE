! rotational origin shift alignment of band-pass limited polar projections in the Fourier domain, gradient based minimizer
module simple_pftcc_shsrch_fm
include 'simple_lib.f08'
use simple_polarft_corrcalc, only: pftcc_glob
use simple_parameters,       only: params_glob
use simple_image,            only: image, image_ptr
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
    procedure          :: calc_phasecorr
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
    subroutine minimize( self, iref, iptcl, irot, score, offset )
        class(pftcc_shsrch_fm), intent(inout) :: self
        integer,                intent(in)    :: iref, iptcl
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
        ! Grid search of both orientations
        call pftcc_glob%bidirectional_shift_search(self%ref, self%ptcl, irot1, self%hn, self%coords, self%grid1, self%grid2)
        ! Maxima
        call self%interpolate_peak(self%grid1, irot1, shift1, score1)
        call self%interpolate_peak(self%grid2, irot2, shift2, score2)
        irot   = irot1
        score  = score1
        offset = shift1
        if( score2 > score1 )then
            irot   = irot2
            score  = score2
            offset = shift2
        endif
        ! Particle shift
        call rotmat2d(pftcc_glob%get_rot(irot), rotmat)
        offset = matmul(offset, rotmat)
    end subroutine minimize

    ! Fourier-Mellin transform appproach without scale
    subroutine calc_phasecorr( self, ind1, ind2, img1, img2, imgcc1, imgcc2, cc, mirror )
        class(pftcc_shsrch_fm), intent(inout) :: self
        integer,                intent(in)    :: ind1, ind2     ! iref,iptcl indices for pftcc
        class(image),           intent(in)    :: img1, img2     ! ref,ptcl corresponding images
        class(image),           intent(inout) :: imgcc1, imgcc2 ! temporary correlation images
        real,                   intent(out)   :: cc
        logical,      optional, intent(in)    :: mirror
        type(image_ptr) :: p1, p2, pcc1, pcc2
        complex :: fcomp1, fcomp2, fcompl, fcompll
        real(dp) :: var1, var2
        real :: corrs(self%pftsz), rmat(2,2), dist(2), loc(2)
        real :: kw, norm, pftcc_ang, cc1,cc2, angstep, alpha,beta,gamma,denom
        integer :: logilims(3,2), cyclims(3,2), cyclimsR(2,2), win_corner(2), phys(2)
        integer :: irot,h,k,l,ll,m,mm,sh,box,c,shlim
        logical :: l_mirr
        l_mirr = .false.
        if( present(mirror) ) l_mirr = mirror
        ! evaluate best orientation from rotational correlation of image magnitudes
        call pftcc_glob%gencorrs_mag_cc(ind1, ind2, corrs, kweight=.true.)
        irot      = maxloc(corrs, dim=1) ! one solution in [0;pi[
        pftcc_ang = pftcc_glob%get_rot(irot)
        ! interpolate rotational peak
        angstep = pftcc_glob%get_rot(2)
        alpha   = merge(corrs(self%pftsz), corrs(irot-1), irot==1)
        beta    = corrs(irot)
        gamma   = merge(corrs(1), corrs(irot+1), irot==self%pftsz)
        denom   = alpha + gamma - 2.*beta
        if( abs(denom) > TINY ) pftcc_ang = pftcc_ang + angstep * 0.5 * (alpha-gamma) / denom
        ! rotation matrix
        call rotmat2d(pftcc_ang, rmat)
        ! phase correlations of image rotated by irot & irot+pi
        call img1%get_cmat_ptr(p1%cmat)
        call img2%get_cmat_ptr(p2%cmat)
        logilims = img1%loop_lims(2)
        cyclims  = img2%loop_lims(3)
        cyclimsR(:,1) = cyclims(1,:)
        cyclimsR(:,2) = cyclims(2,:)
        box = img1%get_box()
        imgcc1 = cmplx(0.,0.)
        imgcc2 = cmplx(0.,0.)
        var1 = 0.d0
        var2 = 0.d0
        do h = logilims(1,1),logilims(1,2)
            do k = logilims(2,1),logilims(2,2)
                sh = nint(hyp(h,k))
                if( sh < pftcc_glob%kfromto(1) )cycle
                if( sh > pftcc_glob%kfromto(2) )cycle
                ! Rotation img2 & bilinear interpolation
                if( l_mirr )then
                    loc = matmul(real([-h,k]),rmat)
                else
                    loc = matmul(real([h,k]),rmat)
                endif
                win_corner = floor(loc)
                dist       = loc - real(win_corner)
                l     = cyci_1d(cyclimsR(:,1), win_corner(1))
                ll    = cyci_1d(cyclimsR(:,1), win_corner(1)+1)
                m     = cyci_1d(cyclimsR(:,2), win_corner(2))
                mm    = cyci_1d(cyclimsR(:,2), win_corner(2)+1)
                phys   = img2%comp_addr_phys(l,m)
                kw     = (1.-dist(1))*(1.-dist(2))
                fcompl = kw * p2%cmat(phys(1), phys(2),1)
                phys   = img2%comp_addr_phys(l,mm)
                kw     = (1.-dist(1))*dist(2)
                fcompl = fcompl + kw * p2%cmat(phys(1), phys(2),1)
                if( l < 0 ) fcompl = conjg(fcompl)
                phys    = img2%comp_addr_phys(ll,m)
                kw      = dist(1)*(1.-dist(2))
                fcompll = kw * p2%cmat(phys(1), phys(2),1)
                phys    = img2%comp_addr_phys(ll,mm)
                kw      = dist(1)*dist(2)
                fcompll = fcompll + kw * p2%cmat(phys(1), phys(2),1)
                if( ll < 0 ) fcompll = conjg(fcompll)
                fcomp2 = fcompl + fcompll
                ! img1
                phys   = img1%comp_addr_phys(h,k)
                fcomp1 = p1%cmat(phys(1), phys(2),1)
                ! rotation by ang
                call imgcc1%set_cmat_at([phys(1),phys(2),1], fcomp1*conjg(fcomp2))
                ! rotation by ang+pi
                call imgcc2%set_cmat_at([phys(1),phys(2),1], fcomp1*      fcomp2)
                ! variances
                if( h == 0 .and. k > 0 ) cycle
                var1 = var1 + real(fcomp1*conjg(fcomp1),dp)
                var2 = var2 + real(fcomp2*conjg(fcomp2),dp)
            end do
        end do
        ! inverse FT
        call imgcc1%ifft
        call imgcc2%ifft
        ! normalization
        norm = real(2.d0*sqrt(var1*var2))
        call imgcc1%div(norm)
        call imgcc2%div(norm)
        ! offsets & maxima
        call imgcc1%get_rmat_ptr(pcc1%rmat)
        call imgcc2%get_rmat_ptr(pcc2%rmat)
        c     = box/2+1
        shlim = floor(self%trslim)
        cc1   = maxval(pcc1%rmat(c-shlim:c+shlim,c-shlim:c+shlim,1))
        cc2   = maxval(pcc2%rmat(c-shlim:c+shlim,c-shlim:c+shlim,1))
        cc    = max(cc1,cc2)
        ! For when the location of maxima matter
        ! Conventions for angle & influence of mirror are TBD
        ! pos1  = maxloc(pcc1%rmat(c-shlim:c+shlim,c-shlim:c+shlim,1))
        ! pos2  = maxloc(pcc2%rmat(c-shlim:c+shlim,c-shlim:c+shlim,1))
        ! cc1   = pcc1%rmat(pos1(1)+c-shlim-1,pos1(2)+c-shlim-1,1)
        ! cc2   = pcc2%rmat(pos2(1)+c-shlim-1,pos2(2)+c-shlim-1,1)
        ! if( cc1 > cc2 )then
        !     offset = real(pos1-shlim-1)
        !     cc     = cc1
        ! else
        !     offset = real(pos2-shlim-1)
        !     cc     = cc2
        !     irot   = irot + self%pftsz
        ! endif
        ! cc & pftcc_glob%gencorr_for_rot_8(ind1,ind2,real(offset,dp),irot) should be close
    end subroutine calc_phasecorr

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
