!@descr: extension of the image class to provide 2D/3D envelope and adaptive masking
module simple_image_msk
use simple_core_module_api
use simple_image,      only: image
use simple_image_bin,  only: image_bin
use simple_parameters, only: parameters
use simple_segmentation
implicit none

public :: image_msk, automask2D, automask2D_support_pix, density_inoutside_mask
private
#include "simple_local_flags.inc"

logical, parameter :: DEBUG         = .false.
logical, parameter :: L_WRITE       = .true.

type, extends(image_bin) :: image_msk
    private
    real    :: msk       = 0.   !< maximum circular mask
    real    :: amsklp    = 0.   !< low-pass limit
    real    :: pix_thres = 0.   !< binarisation threshold
    integer :: edge      = 6    !< edge width
    integer :: binwidth  = 1    !< additional layers to grow
    integer :: idim(3)   = 0    !< image dimension
  contains
    procedure, private :: automask3D_1, automask3D_2
    generic            :: automask3D => automask3D_1, automask3D_2
    procedure          :: estimate_spher_mask_diam
    procedure          :: automask3D_filter
    procedure, private :: automask3D_binarize
    procedure, private :: automask3D_keep_largest_cc
    procedure, private :: env_rproject
end type image_msk

contains

    subroutine automask3D_1( self, params, vol_even, vol_odd, vol_masked, l_tight, pix_thres )
        class(image_msk),  intent(inout) :: self
        class(parameters), intent(in)    :: params
        class(image),      intent(in)    :: vol_even, vol_odd
        class(image),      intent(inout) :: vol_masked
        logical,           intent(in)    :: l_tight
        real, optional,    intent(in)    :: pix_thres
        ! prepare volume for masking
        call vol_masked%copy(vol_even)
        call vol_masked%add(vol_odd)
        call vol_masked%mul(0.5)
        ! automasking
        call self%automask3D_2(params, vol_even, vol_odd, l_tight, pix_thres)
        ! apply mask to volume
        call vol_masked%zero_env_background(self)
        call vol_masked%mul(self)
    end subroutine automask3D_1

    subroutine automask3D_2( self, params, vol_even, vol_odd, l_tight, pix_thres )
        class(image_msk),  intent(inout) :: self
        class(parameters), intent(in)    :: params
        class(image),      intent(in)    :: vol_even, vol_odd
        logical,           intent(in)    :: l_tight
        real, optional,    intent(in)    :: pix_thres 
        type(image) :: vol_filt
        if( vol_even%is_2d() )THROW_HARD('automask3D is intended for volumes only; automask3D')
        write(logfhandle,'(A)'  ) '>>> AUTOMASKING'
        self%amsklp   = params%amsklp
        self%binwidth = params%binwidth
        self%edge     = params%edge
        ! filter
        call self%automask3D_filter(params, vol_even, vol_odd, vol_filt)
        ! binarization
        call self%automask3D_binarize(params, l_tight, pix_thres)
        call vol_filt%kill
        ! add layers
        call self%grow_bins(self%binwidth)
        ! add volume soft edge
        call self%cos_edge(self%edge)
    end subroutine automask3D_2

    subroutine estimate_spher_mask_diam( self, params, vol, amsklp, msk_in_pix )
        class(image_msk),  intent(inout) :: self
        class(parameters), intent(in)    :: params
        class(image),      intent(inout) :: vol
        real,              intent(in)    :: amsklp
        real,              intent(out)   :: msk_in_pix
        real        :: diam
        type(image) :: vol_filt
        if( vol%is_2d() )THROW_HARD('estimate_spher_mask_diam is intended for volumes only; automask3D')
        self%amsklp = amsklp
        ! filter 
        call vol_filt%copy(vol)
        call vol_filt%bp(0., self%amsklp)
        ! transfer image to binary image
        call self%transfer2bimg(vol_filt)
        ! binarization
        call self%automask3D_binarize(params, l_tight=.false.)
        ! mask diameter estimation
        call self%diameter_cc(1, diam)
        diam = diam / self%get_smpd() ! in pixels now
        msk_in_pix = diam / 2. + COSMSKHALFWIDTH
        ! destruct
        call vol_filt%kill
    end subroutine estimate_spher_mask_diam

    subroutine automask3D_filter( self, params, vol_even, vol_odd, vol_filt )
        class(image_msk),  intent(inout) :: self
        class(parameters), intent(in)    :: params
        class(image),      intent(in)    :: vol_even, vol_odd
        class(image),      intent(inout) :: vol_filt
        type(image) :: vol_even_icm, vol_odd_icm
        real,    parameter :: LAM = 100.
        if( vol_even%is_2d() )THROW_HARD('automask3D_filter is intended for volumes only')
        self%amsklp = params%amsklp
        ! ensure vol_filt is constructed
        if( all(vol_filt%get_ldim() == 0) ) call vol_filt%new(vol_even%get_ldim(), vol_even%get_smpd())
        ! ICM filter on local copies to preserve input volumes
        call vol_even_icm%copy(vol_even)
        call vol_odd_icm%copy(vol_odd)
        call vol_even_icm%ICM3D_eo(vol_odd_icm, LAM)
        call vol_filt%copy(vol_even_icm)
        call vol_filt%add(vol_odd_icm)
        call vol_filt%mul(0.5)
        if( L_WRITE .and. params%part == 1 ) call vol_filt%write(string('ICM_avg.mrc'))
        ! low-pass filter
        call vol_filt%bp(0., self%amsklp)
        if( L_WRITE .and. params%part == 1 ) call vol_filt%write(string('ICM_avg_lp.mrc'))
        ! transfer image to binary image
        call self%transfer2bimg(vol_filt)
        call vol_even_icm%kill
        call vol_odd_icm%kill
    end subroutine automask3D_filter

    subroutine automask3D_binarize( self, params, l_tight, pix_thres )
        class(image_msk),  intent(inout) :: self
        class(parameters), intent(in)    :: params
        logical,           intent(in)    :: l_tight
        real, optional,    intent(in)    :: pix_thres
        ! binarize volume
        if( present(pix_thres) )then
            call self%binarize(pix_thres)
        else
            call otsu_img(self, tight=l_tight)
        endif
        if( L_WRITE .and. params%part == 1 ) call self%write(string('binarized.mrc'))
        call self%automask3D_keep_largest_cc(params)
    end subroutine automask3D_binarize

    subroutine automask3D_keep_largest_cc( self, params )
        class(image_msk),  intent(inout) :: self
        class(parameters), intent(in)    :: params
        real,    allocatable :: ccsizes(:)
        type(image_bin)       :: ccimage
        integer              :: loc(1), sz
        call self%set_imat
        ! identify connected components
        call self%find_ccs(ccimage, update_imat=.true.)
        ! extract all cc sizes (in # pixels)
        ccsizes = ccimage%size_ccs()
        sz      = size(ccsizes)
        write(logfhandle,'(A,I7,A)'  ) '>>> FOUND: ', sz, ' CONNECTED COMPONENT(S)'
        if( sz > 1 )then
            loc = maxloc(ccsizes)
            call ccimage%cc2bin(loc(1))
        else
            call ccimage%cc2bin(1)
        endif
        call self%copy_bimg(ccimage)
        if( L_WRITE .and. params%part == 1 ) call self%write(string('largest_cc.mrc'))
        ! destruct
        call ccimage%kill_bimg
        if( allocated(ccsizes) ) deallocate(ccsizes)
    end subroutine automask3D_keep_largest_cc

    ! CALCULATORS

    !>  \brief  volume mask projector
    subroutine env_rproject(self, e, img)
        class(image_msk), intent(inout) :: self   !< projector instance
        class(ori),       intent(inout) :: e      !< Euler angle
        type(image),      intent(inout) :: img    !< resulting projection image
        real, allocatable :: rmat(:,:,:)
        real    :: out_coos(3), maxrad, rad(3), thresh
        real    :: incr_k(3), rvec(3), rvec_k(3)
        integer :: orig(3), i, j, k, sqmaxrad, vec(3)
        ! init
        thresh   = 0.5
        img      = 0.
        orig     = self%idim/2+1
        maxrad   = min(self%msk, real(minval(self%idim(1:2)))/2.-1.)
        sqmaxrad = nint(maxrad**2)
        out_coos = 0.
        rad      = 0.
        vec      = 0
        rmat     = self%get_rmat()
        if( DEBUG )then
            write(logfhandle,*) 'maxrad:       ', maxrad
            write(logfhandle,*) 'sqmaxrad:     ', sqmaxrad
            write(logfhandle,*) 'maxval(rmat): ', maxval(rmat)
            write(logfhandle,*) 'minval(rmat): ', minval(rmat)
            write(logfhandle,*) 'self%idim:    ', self%idim
        endif
        incr_k   = matmul([0., 0., 1.], e%get_mat())
        !$omp parallel do default(shared) private(j,out_coos,rad,i,k,vec,rvec,rvec_k)&
        !$omp schedule(static) proc_bind(close)
        do j=1,self%idim(2)-1
            out_coos(2) = real(j-orig(2))
            rad(2)      = out_coos(2)**2.
            do i=1,self%idim(1)-1
                out_coos(1) = real(i-orig(1))
                rad(1)      = rad(2)+out_coos(1)**2.
                ! check whether we are within the radius
                if(rad(1) > sqmaxrad)cycle
                rvec = real(orig) + matmul([out_coos(:2), 0.], e%get_mat())
                do k=1,self%idim(3)-1
                    out_coos(3) = real(k-orig(3))
                    rad(3)      = rad(1)+out_coos(3)**2.
                    ! check that we are within the radius
                    if( rad(3) > sqmaxrad )cycle
                    rvec_k = rvec + out_coos(3)*incr_k
                    vec    = floor(rvec_k)
                    if( any( rmat(vec(1):vec(1)+1, vec(2):vec(2)+1, vec(3):vec(3)+1) > thresh))then
                        call img%set([i,j,1], 1.)
                        exit
                    endif
                enddo
            enddo
        enddo
        !$omp end parallel do
        deallocate(rmat)
    end subroutine env_rproject

    subroutine density_inoutside_mask( img, lp, msk, nin, nout, nmsk, cccen )
        use simple_segmentation
        use simple_image_bin, only: image_bin
        class(image), intent(inout) :: img
        real,         intent(in)    :: lp, msk
        integer,      intent(out)   :: nin, nout, nmsk
        real,         intent(out)   :: cccen(2)
        type(image_bin)    :: img_bin, cc_img
        real, allocatable :: ccsizes(:)
        integer :: loc, ldim(3)
        real    :: smpd
        ldim = img%get_ldim()
        smpd = img%get_smpd()
        call img_bin%new_bimg(ldim, smpd, wthreads=.false.)
        call img_bin%copy(img)
        call cc_img%new_bimg( ldim, smpd, wthreads=.false.)
        call img_bin%zero_edgeavg
        ! low-pass filter
        call img_bin%bp(0., lp)
        ! binarize with Otsu
        call otsu_img(img_bin)
        call img_bin%set_imat
        ! find the largest connected component
        call img_bin%find_ccs(cc_img)
        ccsizes = cc_img%size_ccs()
        loc     = maxloc(ccsizes,dim=1)
        ! turn it into a binary image for mask creation
        call cc_img%cc2bin(loc)     ! the cc's label is now "1"
        call cc_img%masscen_cc(1, cccen)
        call cc_img%density_inoutside(msk, nin, nout, nmsk)
        call img_bin%kill_bimg
        call cc_img%kill_bimg
        if( allocated(ccsizes) ) deallocate(ccsizes)
    end subroutine density_inoutside_mask

    subroutine automask2D( params, imgs, ngrow, winsz, edge, diams, shifts, write2disk )
        class(parameters), intent(in)    :: params
        class(image),      intent(inout) :: imgs(:)
        integer,           intent(in)    :: ngrow, winsz, edge
        real, allocatable, intent(inout) :: diams(:), shifts(:,:)
        logical, optional, intent(in)    :: write2disk
        type(image_bin),   allocatable   :: cc_img(:)
        integer :: i, n
        logical :: l_write
        n = size(imgs)
        l_write = .false.
        if( present(write2disk) ) l_write = write2disk
        l_write = l_write .and. params%part.eq.1
        if( allocated(diams)  ) deallocate(diams)
        if( allocated(shifts) ) deallocate(shifts)
        ! allocate
        allocate(diams(n), shifts(n,2), source=0.)
        allocate(cc_img(n))
        if( trim(params%automsk).eq.'tight' )then
            write(logfhandle,'(A)') '>>> 2D AUTOMASKING, TIGHT'
        else
            write(logfhandle,'(A)') '>>> 2D AUTOMASKING'
        endif
        call imgs(1)%memoize_mask_coords
        !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
        do i = 1,n
            call automask2D_binary_one(params, imgs(i), ngrow, winsz, edge, cc_img(i), diams(i), shifts(i,:))
            ! apply cosine egde to soften mask (to avoid Fourier artefacts)
            call imgs(i)%zero_and_unflag_ft
            call cc_img(i)%cos_edge(edge,imgs(i))
        end do
        !$omp end parallel do
        ! destruct
        do i = 1,n
            call cc_img(i)%write(string('binarized_automask2D.mrc'), i)
            call cc_img(i)%kill_bimg
            call imgs(i)%write(string('masks_automask2D.mrc'), i)
        end do
        deallocate(cc_img)
    end subroutine automask2D

    subroutine automask2D_support_pix( params, img, ngrow, winsz, edge, pix, diam, shift )
        class(parameters), intent(in)    :: params
        class(image),      intent(in)    :: img
        integer,           intent(in)    :: ngrow, winsz, edge
        integer, allocatable, intent(inout) :: pix(:,:)
        real, optional,    intent(out)   :: diam, shift(2)
        type(image_bin)    :: bin_mask
        integer, allocatable :: imat(:,:,:)
        real :: diam_local, shift_local(2)
        call automask2D_binary_one(params, img, ngrow, winsz, edge, bin_mask, diam_local, shift_local, &
            &l_fallback_spherical=.false.)
        call bin_mask%get_imat(imat)
        call binary_imat_to_pix(imat, pix)
        if( present(diam)  ) diam  = diam_local
        if( present(shift) ) shift = shift_local
        call bin_mask%kill_bimg
        if( allocated(imat) ) deallocate(imat)
    end subroutine automask2D_support_pix

    subroutine automask2D_binary_one( params, img, ngrow, winsz, edge, bin_mask, diam, shift, l_fallback_spherical )
        class(parameters), intent(in)    :: params
        class(image),      intent(in)    :: img
        integer,           intent(in)    :: ngrow, winsz, edge
        type(image_bin),   intent(inout) :: bin_mask
        real,              intent(out)   :: diam, shift(2)
        logical, optional, intent(in)    :: l_fallback_spherical
        type(image_bin)    :: img_bin
        real, allocatable  :: ccsizes(:)
        integer :: loc, ldim(3)
        real    :: smpd, xyz(3), mskrad
        logical :: l_spherical_fallback
        l_spherical_fallback = .true.
        if( present(l_fallback_spherical) ) l_spherical_fallback = l_fallback_spherical
        ldim = img%get_ldim()
        smpd = img%get_smpd()
        mskrad = automask2D_mskrad(params, ldim)
        call img_bin%new_bimg(ldim, smpd, wthreads=.false.)
        call img_bin%copy(img)
        call bin_mask%new_bimg(ldim, smpd, wthreads=.false.)
        call img_bin%zero_edgeavg
        ! dampens below zero (object positive in class averages/reprojs)
        call img_bin%div_below(0.,10.)
        ! low-pass filter
        call img_bin%bp(0., params%amsklp)
        ! filter with non-local means
        call img_bin%NLmean2D
        ! binarize with Otsu
        call otsu_img(img_bin, mskrad=mskrad, positive=trim(params%automsk).eq.'tight')
        call img_bin%masscen(xyz)
        shift = xyz(:2)
        call img_bin%set_imat
        ! grow ngrow layers
        if( ngrow > 0 ) call img_bin%grow_bins(ngrow)
        ! find the largest connected component
        call img_bin%find_ccs(bin_mask)
        ccsizes = bin_mask%size_ccs()
        loc     = maxloc(ccsizes,dim=1)
        ! estimate its diameter
        call bin_mask%diameter_cc(loc, diam)
        if( diam <= TINY )then
            shift = 0.
            if( l_spherical_fallback )then
                diam  = 2.*max(1., mskrad-real(edge)-COSMSKHALFWIDTH)*smpd
                call bin_mask%disc(ldim, smpd, diam/(2.*smpd))
                call bin_mask%set_imat
            else
                diam = 0.
                call bin_mask%new_bimg(ldim, smpd, wthreads=.false.)
            endif
        else if( l_spherical_fallback .and. diam > 2.*(mskrad+real(ngrow))*smpd )then
            ! incorrect component was chosen, fall back on circular support
            shift = 0.
            diam  = 2.*max(1., mskrad-real(edge)-COSMSKHALFWIDTH)*smpd
            call bin_mask%disc(ldim, smpd, diam/(2.*smpd))
            call bin_mask%set_imat
        else
            ! turn it into a binary image for mask creation
            call bin_mask%cc2bin(loc)
            call bin_mask%masscen(xyz)
            shift = xyz(:2)
            ! median filter to smoothen
            if( winsz > 0 )then
                call bin_mask%real_space_filter(winsz, 'median')
                call bin_mask%set_imat
            endif
            ! fill-in holes
            call bin_mask%set_edgecc2background
        endif
        call img_bin%kill_bimg
        if( allocated(ccsizes) ) deallocate(ccsizes)
    end subroutine automask2D_binary_one

    real function automask2D_mskrad( params, ldim )
        class(parameters), intent(in) :: params
        integer,           intent(in) :: ldim(3)
        real :: max_mskrad
        max_mskrad = max(1., real(ldim(1)) / 2. - COSMSKHALFWIDTH - 1.)
        automask2D_mskrad = params%msk
        if( params%box_crop > 0 .and. ldim(1) == params%box_crop .and. params%msk_crop > TINY )then
            automask2D_mskrad = params%msk_crop
        endif
        if( automask2D_mskrad <= TINY ) automask2D_mskrad = max_mskrad
        automask2D_mskrad = min(automask2D_mskrad, max_mskrad)
    end function automask2D_mskrad

    subroutine binary_imat_to_pix( imat, pix )
        integer, intent(in) :: imat(:,:,:)
        integer, allocatable, intent(inout) :: pix(:,:)
        integer :: i, j, ipix, npix
        if( size(imat,3) /= 1 ) THROW_HARD('2D binary mask expected; binary_imat_to_pix')
        if( allocated(pix) ) deallocate(pix)
        npix = count(imat(:,:,1) > 0)
        allocate(pix(2,npix), source=0)
        if( npix < 1 ) return
        ipix = 0
        do j = 1, size(imat,2)
            do i = 1, size(imat,1)
                if( imat(i,j,1) <= 0 ) cycle
                ipix = ipix + 1
                pix(:,ipix) = [i,j]
            end do
        end do
        if( ipix /= npix ) THROW_HARD('support pixel count mismatch; binary_imat_to_pix')
    end subroutine binary_imat_to_pix

end module simple_image_msk
