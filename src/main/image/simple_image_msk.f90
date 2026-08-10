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
logical, parameter :: L_WRITE       = .false.

type, extends(image_bin) :: image_msk
    private
    real    :: msk       = 0.   !< maximum circular mask
    real    :: amsklp    = 0.   !< low-pass limit
    real    :: pix_thres = 0.   !< binarisation threshold
    integer :: edge      = 6    !< edge width
    integer :: binwidth  = 1    !< additional layers to grow
    integer :: idim(3)   = 0    !< image dimension
  contains
    procedure          :: automask3D
    procedure          :: apply_3Dmask
    procedure          :: estimate_spher_mask_diam
    procedure          :: envmask3D_from_lmask
    procedure, private :: automask3D_binarize
    procedure, private :: env_rproject
end type image_msk

contains

    subroutine automask3D( self, params, vol, l_tight, pix_thres, vol_masked )
        class(image_msk),       intent(inout) :: self
        class(parameters),      intent(in)    :: params
        class(image),           intent(in)    :: vol
        logical,                intent(in)    :: l_tight
        real,         optional, intent(in)    :: pix_thres
        class(image), optional, intent(inout) :: vol_masked
        if( vol%is_2d() )THROW_HARD('automask3D is intended for volumes only')
        write(logfhandle,'(A)') '>>> AUTOMASKING'
        ! parameters
        self%amsklp   = params%amsklp
        self%binwidth = params%binwidth
        self%edge     = params%edge
        ! low-pass smoothing
        call self%new_bimg(vol%get_ldim(), vol%get_smpd())
        call self%copy(vol)
        call self%bp(0., self%amsklp)
        call self%ifft
        if( L_WRITE .and. params%part == 1 ) call self%write(string('automask_lowpass.mrc'))
        ! segmentation
        call self%automask3D_binarize(params, l_tight, pix_thres)
        ! Morphological growth and soft edge
        call self%grow_bins(self%binwidth)
        call self%cos_edge(self%edge)
        ! optionally mask
        if( present(vol_masked) )then
            call vol_masked%copy( vol )
            call self%apply_3Dmask( vol_masked )
        endif
    end subroutine automask3D

    subroutine apply_3Dmask( self, vol )
        class(image_msk), intent(in)    :: self
        class(image),     intent(inout) :: vol
        call vol%zero_env_background(self)
        call vol%mul(self)
    end subroutine apply_3Dmask

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

    !>  Builds a soft envelope mask from an already-segmented binary field. This is
    !!  the topology and morphology tail shared by any segmentation front end: a
    !!  smoothness prior regularizes boundary area but cannot enforce a connected
    !!  result, so component selection and hole filling have to happen here.
    !!  minvol_frac keeps every component at least that fraction of the largest,
    !!  which preserves a second well-resolved domain joined by a disordered linker
    !!  instead of silently discarding it the way largest-component selection does.
    subroutine envmask3D_from_lmask( self, lmask, smpd, binwidth, edge, minvol_frac, l_fill_holes, n_ccs, n_ccs_kept )
        class(image_msk), intent(inout) :: self
        logical,          intent(in)    :: lmask(:,:,:)
        real,             intent(in)    :: smpd
        integer,          intent(in)    :: binwidth, edge
        real,             intent(in)    :: minvol_frac
        logical,          intent(in)    :: l_fill_holes
        integer,          intent(out)   :: n_ccs, n_ccs_kept
        type(image)          :: vol_bin
        type(image_bin)      :: ccimage
        integer, allocatable :: ccsizes(:), imat(:,:,:), keepmat(:,:,:)
        integer              :: ldim_l(3), minsz, icc, i, j, k
        ldim_l = shape(lmask)
        n_ccs      = 0
        n_ccs_kept = 0
        if( ldim_l(3) < 2 ) THROW_HARD('envmask3D_from_lmask is intended for volumes only')
        if( .not.any(lmask) )then
            THROW_WARN('empty binary field, no envelope produced; envmask3D_from_lmask')
            return
        endif
        call vol_bin%new(ldim_l, smpd)
        call vol_bin%set_rmat(merge(1., 0., lmask), .false.)
        call self%transfer2bimg(vol_bin)
        call vol_bin%kill
        call self%find_ccs(ccimage)
        ccsizes = ccimage%size_ccs()
        n_ccs   = size(ccsizes)
        if( maxval(ccsizes) < 1 )then
            THROW_WARN('no connected components found; envmask3D_from_lmask')
            call ccimage%kill_bimg
            return
        endif
        minsz = 1
        if( minvol_frac > 0. ) minsz = max(1, nint(minvol_frac * real(maxval(ccsizes))))
        n_ccs_kept = count(ccsizes >= minsz)
        write(logfhandle,'(A,I7,A,I7,A,I10,A)') '>>> ENVELOPE COMPONENTS: ', n_ccs, ' FOUND, ', &
            &n_ccs_kept, ' KEPT (min size ', minsz, ' voxels)'
        call ccimage%get_imat(imat)
        allocate(keepmat(ldim_l(1),ldim_l(2),ldim_l(3)), source=0)
        !$omp parallel do collapse(3) schedule(static) default(shared) private(i,j,k,icc) proc_bind(close)
        do k = 1, ldim_l(3)
            do j = 1, ldim_l(2)
                do i = 1, ldim_l(1)
                    icc = imat(i,j,k)
                    if( icc < 1 ) cycle
                    if( ccsizes(icc) >= minsz ) keepmat(i,j,k) = 1
                end do
            end do
        end do
        !$omp end parallel do
        call vol_bin%new(ldim_l, smpd)
        call vol_bin%set_rmat(real(keepmat), .false.)
        call self%transfer2bimg(vol_bin)
        call vol_bin%kill
        ! Hole filling precedes dilation so interior voids are closed exactly rather
        ! than partially bridged by the grown layers.
        if( l_fill_holes ) call self%set_edgecc2background
        if( binwidth > 0 ) call self%grow_bins(binwidth)
        call self%cos_edge(edge)
        call ccimage%kill_bimg
        deallocate(ccsizes, imat, keepmat)
    end subroutine envmask3D_from_lmask

    subroutine automask3D_binarize( self, params, l_tight, pix_thres )
        class(image_msk),  intent(inout) :: self
        class(parameters), intent(in)    :: params
        logical,           intent(in)    :: l_tight
        real, optional,    intent(in)    :: pix_thres
        integer, allocatable :: cc_sz(:)
        type(image_bin)      :: vol_ccs
        integer              :: i
        ! binarization
        if( present(pix_thres) )then
            call self%binarize(pix_thres)
        else
            call otsu_img(self, tight=l_tight)
        endif
        call self%set_imat
        if( L_WRITE .and. params%part == 1 ) call self%write(string('automask_binarized.mrc'))
        ! identify connected components
        call self%find_ccs(vol_ccs, update_imat=.true.)
        cc_sz = vol_ccs%size_ccs()
        write(logfhandle,'(A,I7,A)'  ) '>>> FOUND: ', size(cc_sz), ' CONNECTED COMPONENT(S)'
        ! extract largest CC
        i = maxloc(cc_sz, dim=1)
        call vol_ccs%cc2bin(i)
        call self%copy_bimg(vol_ccs)
        if( L_WRITE .and. params%part == 1 ) call self%write(string('automask_cc.mrc'))
        call vol_ccs%kill_bimg
        if( allocated(cc_sz) ) deallocate(cc_sz)
    end subroutine automask3D_binarize

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

    subroutine automask2D( params, imgs, ngrow, winsz, edge, diams, shifts, write2disk, min_diams, verbose )
        class(parameters),              intent(in)    :: params
        class(image),                   intent(inout) :: imgs(:)
        integer,                        intent(in)    :: ngrow, winsz, edge
        real,              allocatable, intent(inout) :: diams(:), shifts(:,:)
        logical, optional,              intent(in)    :: write2disk
        real,    optional, allocatable, intent(inout) :: min_diams(:)
        logical, optional,              intent(in)    :: verbose
        type(image_bin),   allocatable                :: cc_img(:)
        type(image_bin)                               :: cc_min_dist
        integer :: i, n
        real    :: mm(2)
        logical :: l_write, l_verbose
        n = size(imgs)
        l_write = .false.
        if( present(write2disk) ) l_write = write2disk
        l_write = l_write .and. params%part.eq.1
        l_verbose = .true.
        if( present(verbose) ) l_verbose = verbose
        if( allocated(diams)     ) deallocate(diams)
        if( allocated(shifts)    ) deallocate(shifts)
        ! allocate
        allocate(diams(n), shifts(n,2), source=0.)
        allocate(cc_img(n))
        if( l_verbose )then
            if( trim(params%automsk).eq.'tight' )then
                write(logfhandle,'(A)') '>>> 2D AUTOMASKING, TIGHT'
            else
                write(logfhandle,'(A)') '>>> 2D AUTOMASKING'
            endif
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
        if( present(min_diams) ) then
            if( allocated(min_diams) ) deallocate(min_diams)
            allocate(min_diams(n))
            do i = 1,n
                call cc_min_dist%copy(cc_img(i))
                call cc_min_dist%distance_transform()
                mm = cc_min_dist%minmax()
                min_diams(i) = 2.*mm(2)*imgs(i)%get_smpd()
                write(logfhandle,'(A,I0,A,F7.2,A)') '>>> MIN DIAMETER OF MASK FOR IMAGE ', i, ': ', min_diams(i), ' A'
                call cc_min_dist%kill_bimg()
            end do
        end if
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
