module simple_tseries_graphene_subtr
use simple_parameters,   only: parameters, params_glob
use simple_cmdline,  only: cmdline
use simple_polarizer,only: polarizer
use simple_image,    only: image
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_math
use simple_rnd
use simple_defs ! singleton
use simple_atoms,    only: atoms
implicit none

public :: init_graphene_subtr, calc_peaks, remove_lattices, kill_graphene_subtr
private

real, parameter :: angular_threshold = 3.0  ! Degrees
real, parameter :: GRAPHENE_BAND3    = 1.06 ! Angstroms
real, parameter :: REMOVAL_HWIDTH1 = 3.          ! pixels, obscuring half-width 1
real, parameter :: REMOVAL_HWIDTH2 = sqrt(6.)    ! obscuring half-width 2
real, parameter :: REMOVAL_HWIDTH3 = sqrt(3.)    ! obscuring half-width 3

type(polarft_corrcalc) :: pftcc
type(polarizer)        :: pspec_img
type(image)            :: filter_img
type(parameters)       :: params
real                   :: angstep
integer                :: nrots

contains

    subroutine init_graphene_subtr(box, smpd)
        integer,  intent(in) :: box
        real,     intent(in) :: smpd
        real,      parameter :: gdist = 1.42
        real,      parameter :: gx    = 0.5
        real,      parameter :: gy    = sin(PI/3.)
        type(atoms)          :: sheet
        type(image)          :: vol, img
        type(cmdline)        :: cline
        real,    allocatable :: x(:)
        logical, allocatable :: resmask(:)
        real    :: l, cen, offset,y
        integer :: phys2d(3),phys3d(3),lims(3,2),h,k,i,j,nx,ny,cnt,natoms,imgcen,sh,nyq
        call kill_graphene_subtr
        ! sheet
        imgcen = box/2+1
        cen    = smpd*real(box/2+1-1)
        l      = real(box)*smpd
        nx     = ceiling(2.*real(box)*smpd/(3.*gdist))+1
        ny     = ceiling(real(box)*smpd/(gy*gdist))
        allocate(x(nx),source=0.)
        do i = 2,nx
            if( is_even(i) )then
                x(i) = x(i-1) + gdist
            else
                x(i) = x(i-1) + 2.*gdist
            endif
        enddo
        natoms = nx*ny
        call sheet%new(natoms)
        cnt = 0
        do j = 1,ny
            y      = real(j-1)*gy*gdist
            offset = 0.
            if( is_even(j) ) offset = -1.5*gdist
            do i = 1,nx
                cnt = cnt + 1
                call sheet%set_coord(cnt,[x(i)+offset,y,cen])
                call sheet%set_num(cnt,cnt)
                call sheet%set_name(cnt,'C   ')
                call sheet%set_chain(cnt,'A')
                call sheet%set_element(cnt,'C ')
                call sheet%set_resnum(cnt,1)
                call sheet%set_occupancy(cnt,1.)
            enddo
        enddo
        ! convolution
        call vol%new([box,box,box],smpd)
        call vol%zero
        call sheet%convolve(vol, 3.*gdist)
        call vol%mask(real(box/2-5),'soft')
        ! extraction, copy plane so no interpolation
        call vol%mul(real(box))
        call vol%fft
        call img%new([box,box,1],smpd)
        call img%zero_and_flag_ft
        lims = img%loop_lims(2)
        do h = lims(1,1),lims(1,2)
            do k = lims(2,1),lims(2,2)
                phys2d = img%comp_addr_phys([h,k,0])
                phys3d = vol%comp_addr_phys([h,k,0])
                call img%set_fcomp([h,k,0], phys2d, vol%get_fcomp([h,k,0],phys3d))
            enddo
        enddo
        call pspec_img%new([box,box,1],smpd)
        call img%ft2img('sqrt', pspec_img)
        ! resolutions mask & filter
        call filter_img%new([box,box,1],smpd)
        filter_img = 1.
        nyq        = pspec_img%get_nyq()
        resmask    = calc_3bands_mask(box,smpd)
        do i = 1,box
            h = i-imgcen
            do j = 1,box
                k  = j-imgcen
                sh = nint(sqrt(real(h*h+k*k)))
                if( sh > nyq )then
                    call filter_img%set([i,j,1],0.)
                else
                    if( sh == 0 )then
                        call filter_img%set([i,j,1],0.)
                    else
                        if( resmask(sh) )call filter_img%set([i,j,1],0.)
                    endif
                endif
            enddo
        enddo
        call pspec_img%mul(filter_img)
        ! pftcc
        call cline%set('prg',      'dummy')
        call cline%set('smpd',     smpd)
        call cline%set('box',      real(box))
        call cline%set('nthr',     1.)
        call cline%set('match_filt','no')
        call cline%set('ctf',      'no')
        call params%new(cline, silent=.true.)
        params_glob%ring2      = nint(0.4*real(box))
        params_glob%kfromto(1) = 5
        params_glob%kfromto(2) = nyq-1
        params_glob%kstop      = params_glob%kfromto(2)
        call pftcc%new(1,[1,1])
        angstep = abs(pftcc%get_rot(2)-pftcc%get_rot(1))
        nrots   = pftcc%get_nrots()
        ! reference polar coordinates
        call pspec_img%init_polarizer(pftcc, KBALPHA)
        call pspec_img%fft
        call pspec_img%polarize(pftcc,1, .false., .true.)
        call pftcc%memoize_ffts
        ! cleanup
        call sheet%kill
        call cline%kill
        call img%kill
    end subroutine init_graphene_subtr

    subroutine calc_peaks( img_in, ang1, ang2 )
        class(image), intent(inout) :: img_in
        real,         intent(out)   :: ang1, ang2
        real    :: corrs(nrots), ang, corr, rot, diffang
        logical :: corrs_mask(nrots)
        integer :: i,j,il,ir
        ! rotational average
        call img_in%roavg(60,pspec_img)
        ! bands filtering
        call pspec_img%mul(filter_img)
        ! polar coordinates
        call pspec_img%fft()
        call pspec_img%polarize(pftcc,1,.true.,.true.)
        ! rotational correlations
        call pftcc%memoize_ffts
        call pftcc%gencorrs(1,1,corrs)
        ! mask out non peak-shape values
        do i = 1,nrots
            il = i-1
            ir = i+1
            if( il < 1 ) il = il+nrots
            if( ir > nrots ) ir = ir-nrots
            if( corrs(il) < corrs(i) .and. corrs(ir) < corrs(i) )then
                corrs_mask(i) = .true.
            else
                corrs_mask(i) = .false.
            endif
        enddo
        ! First spectral peak interpolation
        call interp_peak(ang1, corr)
        do while( ang1 > 60. )
            ang1 = ang1-60.
        end do
        ! mask out first peak
        do i = 1,nrots
            rot = 360.-pftcc%get_rot(i)
            do j = 0,5
                ang = ang1 + real(j)*60.
                diffang = abs(rot-ang)
                if( diffang > 180. ) diffang = 360.- diffang
                if( diffang < angular_threshold ) corrs_mask(i) = .false.
            enddo
        enddo
        if( count(corrs_mask) == 0 )then
            ! could not find second peak within angular threshold
            call range_convention(ang1)
            ang2 = ang1
            return
        endif
        ! Second spectral peak interpolation
        call interp_peak(ang2, corr)
        ! angular convention
        call range_convention(ang1)
        call range_convention(ang2)
        ! Sort peaks
        if( ang1 > ang2 )then
            ang  = ang1
            ang1 = ang2
            ang2 = ang
        endif
        contains

            subroutine range_convention( ang )
                real, intent(inout) :: ang
                ang = -(ang + 30.)
                do while( ang < 0. )
                    ang = ang+360.
                end do
                do while( ang > 60. )
                    ang = ang-60.
                end do
            end subroutine range_convention

            subroutine interp_peak(ang, corr)
                real, intent(out) :: ang, corr
                real    :: denom, alpha, beta, gamma, drot
                integer :: maxind, ind
                maxind = maxloc(corrs,dim=1,mask=corrs_mask)
                beta   = corrs(maxind)
                ind    = merge(nrots, maxind-1, maxind==1)
                alpha  = corrs(ind)
                ind    = merge(1, maxind+1, maxind==pftcc%get_nrots())
                gamma  = corrs(ind)
                corr   = beta
                ang    = 360.-pftcc%get_rot(maxind)
                denom  = alpha+gamma-2.*beta
                if( abs(denom) < TINY )return
                drot = 0.5 * (alpha-gamma) / denom
                corr = beta -0.25*(alpha-gamma)*drot
                ang  = ang + drot*angstep
            end subroutine interp_peak

    end subroutine calc_peaks

    subroutine remove_lattices(img, ang1, ang2)
        class(image), intent(inout) :: img
        real,         intent(in)    :: ang1, ang2
        real    :: smpd
        integer :: band_inds(3) ! respectively resolutions 2.14, 1.23, 1.06 angstroms
        integer :: cen(2), ldim(3)
        ldim    = img%get_ldim()
        ldim(3) = 1
        smpd    = img%get_smpd()
        cen     = ldim(1:2)/2 + 1
        band_inds(1) = calc_fourier_index(max(2.*smpd,GRAPHENE_BAND1),ldim(1),smpd)
        band_inds(2) = calc_fourier_index(max(2.*smpd,GRAPHENE_BAND2),ldim(1),smpd)
        band_inds(3) = calc_fourier_index(max(2.*smpd,GRAPHENE_BAND3),ldim(1),smpd)
        call img%fft
        call remove_lattice(deg2rad(ang1))
        call remove_lattice(deg2rad(ang2))
        call img%ifft
        contains

            subroutine remove_lattice(ang)
                real, intent(in) :: ang
                real,  parameter :: pionsix   = PI/6.
                real,  parameter :: pionthree = PI/3.
                real :: b1(2), b1b2(2), rp, r, angp, rot, rxy(2), Rm(2,2)
                integer :: rot_ind
                rp   = real(band_inds(2))
                r    = rp / (2.*cos(pionsix))
                angp = ang+pionsix
                b1b2 = rp*[cos(ang),sin(ang)]
                ! lattice parameters (eg, band 1)
                b1   = r*[cos(angp),sin(angp)]
                ! remove peaks
                !$omp parallel do default(shared) private(rot_ind,rot,Rm,rxy) schedule(static) proc_bind(close)
                do rot_ind = 0,6
                    rot = real(rot_ind)*pionthree
                    Rm(1,:) = [cos(rot), -sin(rot)]
                    Rm(2,:) = [-Rm(1,2), Rm(1,1)]
                    ! band 1 (=b1); resolution ~ 2.14 angs
                    rxy = matmul(Rm,b1) + real(cen)
                    call obscure_peak(rxy, REMOVAL_HWIDTH1)
                    ! 2. * band 1 (=2*b1); resolution ~ 1.06; generally the weaker peak
                    rxy = matmul(Rm,2.*b1) + real(cen)
                    call obscure_peak(rxy, REMOVAL_HWIDTH3)
                    ! band 2 (=b1+b2); resolution ~ 1.23
                    rxy = matmul(Rm,b1b2) + real(cen)
                    call obscure_peak(rxy, REMOVAL_HWIDTH2)
                enddo
                !$omp end parallel do
            end subroutine remove_lattice

            subroutine obscure_peak(vec, hwidth)
                real, intent(in) :: vec(2), hwidth
                real    :: dist, w
                integer :: lims(2,2),phys(3),i,j,h,k
                lims(:,1) = floor(vec-hwidth)
                lims(:,2) = ceiling(vec+hwidth)
                do i=lims(1,1),lims(1,2)
                    if( i <  1 .or. i > ldim(1) )cycle
                    do j=lims(2,1),lims(2,2)
                        if( j < 1 .or. j > ldim(2) )cycle
                        if( i < cen(1) ) cycle
                        dist = sqrt(sum((real([i,j])-vec)**2.)) / hwidth
                        if( dist > 1. ) cycle
                        h = i - cen(1)
                        k = j - cen(2)
                        phys = img%comp_addr_phys(h,k,0)
                        ! quartic kernel (Epanechnikov squared)
                        w = 1.-(1.-dist**2.)**2.
                        call img%mul_cmat_at(phys, w)
                    enddo
                enddo
            end subroutine obscure_peak

    end subroutine remove_lattices

    !> \brief calculate logical mask filtering out the Graphene bands
    function calc_3bands_mask( box, smpd ) result( mask )
        integer, intent(in)  :: box
        real,    intent(in)  :: smpd
        real,    allocatable :: res(:), sqdiff_band(:)
        logical, allocatable :: mask(:)
        integer, parameter   :: NBANDS = 3
        integer              :: loc(NBANDS), n, i
        res = get_resarr( box, smpd )
        n   = size(res)
        allocate(mask(n), source=.true.)
        allocate(sqdiff_band(n), source=(res - GRAPHENE_BAND1)**2.0)
        loc = minnloc(sqdiff_band, NBANDS)
        do i=1,NBANDS
            mask(loc(i)) = .false.
        end do
        deallocate(sqdiff_band)
        allocate(sqdiff_band(n), source=(res - GRAPHENE_BAND2)**2.0)
        loc = minnloc(sqdiff_band, NBANDS)
        do i=1,NBANDS
            mask(loc(i)) = .false.
        end do
        deallocate(sqdiff_band)
        allocate(sqdiff_band(n), source=(res - GRAPHENE_BAND3)**2.0)
        loc = minnloc(sqdiff_band, NBANDS)
        do i=1,NBANDS
            mask(loc(i)) = .false.
        end do
    end function calc_3bands_mask

    ! DESTRUCTOR

    subroutine kill_graphene_subtr
        call pftcc%kill
        call pspec_img%kill
        call filter_img%kill
        angstep = 0.
        nrots   = 0
    end subroutine kill_graphene_subtr

end module simple_tseries_graphene_subtr
