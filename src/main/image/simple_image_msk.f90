! 2D/3D envelope and adaptive masking
module simple_image_msk
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image,      only: image
use simple_image_bin,  only: image_bin
use simple_parameters, only: params_glob
use simple_segmentation
implicit none

public :: image_msk, automask2D, density_inoutside_mask
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
    procedure          :: mask_from_pdb
    procedure, private :: env_rproject
end type image_msk

contains

    subroutine automask3D_1( self, vol_even, vol_odd, vol_masked, l_tight, pix_thres )
        class(image_msk),  intent(inout) :: self
        class(image),      intent(inout) :: vol_even, vol_odd, vol_masked
        logical,           intent(in)    :: l_tight
        real, optional,    intent(in)    :: pix_thres
        ! prepare volume for masking
        call vol_masked%copy(vol_even)
        call vol_masked%add(vol_odd)
        call vol_masked%mul(0.5)
        ! automasking
        call self%automask3D_2(vol_even, vol_odd, l_tight, pix_thres)
        ! apply mask to volume
        call vol_masked%zero_env_background(self)
        call vol_masked%mul(self)
    end subroutine automask3D_1

    subroutine automask3D_2( self, vol_even, vol_odd, l_tight, pix_thres )
        class(image_msk),  intent(inout) :: self
        class(image),      intent(inout) :: vol_even, vol_odd
        logical,           intent(in)    :: l_tight
        real, optional,    intent(in)    :: pix_thres 
        type(image) :: vol_filt
        if( vol_even%is_2d() )THROW_HARD('automask3D is intended for volumes only; automask3D')
        write(logfhandle,'(A)'  ) '>>> AUTOMASKING'
        self%amsklp   = params_glob%amsklp
        self%binwidth = params_glob%binwidth
        self%edge     = params_glob%edge
        ! filter 
        call self%automask3D_filter(vol_even, vol_odd, vol_filt)
        ! binarization
        call self%automask3D_binarize(l_tight, pix_thres)
        ! add layers
        call self%grow_bins(self%binwidth)
        ! add volume soft edge
        call self%cos_edge(self%edge)
        ! destruct
        call vol_filt%kill
    end subroutine automask3D_2

    subroutine estimate_spher_mask_diam( self, vol, amsklp, msk_in_pix )
        class(image_msk),  intent(inout) :: self
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
        call self%automask3D_binarize(l_tight=.false.)
        ! mask diameter estimation
        call self%diameter_cc(1, diam)
        diam = diam / self%get_smpd() ! in pixels now
        msk_in_pix = diam / 2. + COSMSKHALFWIDTH
        ! destruct
        call vol_filt%kill
    end subroutine estimate_spher_mask_diam

    subroutine automask3D_filter( self, vol_even, vol_odd, vol_filt )
        class(image_msk), intent(inout) :: self
        class(image),     intent(inout) :: vol_even, vol_odd, vol_filt
        real,    parameter :: LAM = 100.
        if( vol_even%is_2d() )THROW_HARD('automask3D_filter is intended for volumes only')
        self%amsklp = params_glob%amsklp
        ! ICM filter
        call vol_even%ICM3D_eo(vol_odd, LAM)
        call vol_filt%copy(vol_even)
        call vol_filt%add(vol_odd)
        call vol_filt%mul(0.5)
        if( L_WRITE .and. params_glob%part == 1 ) call vol_filt%write(string('ICM_avg.mrc'))
        ! low-pass filter
        call vol_filt%bp(0., self%amsklp)
        if( L_WRITE .and. params_glob%part == 1 ) call vol_filt%write(string('ICM_avg_lp.mrc'))
        ! transfer image to binary image
        call self%transfer2bimg(vol_filt)
    end subroutine automask3D_filter

    subroutine automask3D_binarize( self, l_tight, pix_thres )
        class(image_msk),  intent(inout) :: self
        logical,           intent(in)    :: l_tight
        real, optional,    intent(in)    :: pix_thres
        real,    allocatable :: ccsizes(:)
        type(image_bin)       :: ccimage
        integer              :: loc(1), sz
        ! binarize volume
        if( present(pix_thres) )then
            call self%binarize(pix_thres)
        else
            call otsu_img(self, tight=l_tight)
        endif
        call self%set_imat
        if( L_WRITE .and. params_glob%part == 1 ) call self%write(string('binarized.mrc'))
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
        if( L_WRITE .and. params_glob%part == 1 ) call self%write(string('largest_cc.mrc'))
        ! destruct
        call ccimage%kill_bimg
        if( allocated(ccsizes) ) deallocate(ccsizes)
    end subroutine automask3D_binarize

    subroutine mask_from_pdb( self,  pdb, vol_inout, os, pdbout )
        use simple_atoms, only: atoms
        class(image_msk),        intent(inout) :: self
        type(atoms),             intent(inout) :: pdb
        class(image),            intent(inout) :: vol_inout
        class(oris),   optional, intent(inout) :: os
        class(string), optional, intent(inout) :: pdbout
        type(image) :: distimg, cos_img
        type(atoms) :: shifted_pdb
        real        :: centre(3), shift(3), pdb_center(3), minmax(2), radius, smpd
        integer     :: i
        logical     :: was_ft
        if( vol_inout%is_2d() ) THROW_HARD('intended for volumes only; mask_from_pdb')
        was_ft = vol_inout%is_ft()
        smpd   = vol_inout%get_smpd()
        call self%new(vol_inout%get_ldim(), smpd)
        call distimg%new(vol_inout%get_ldim(), smpd)
        call distimg%cendist()
        pdb_center = pdb%get_geom_center()
        centre     = real(self%get_ldim()-1)/2. * smpd
        shift      = ( pdb_center - centre ) / smpd
        if( params_glob%binwidth == 0 ) then
            radius = smpd
        else
            radius = real(params_glob%binwidth) * smpd
        endif
        if( params_glob%center .eq. 'yes' )then
            ! shift volume & oris
            call vol_inout%shift(shift)
            if( present(os) ) call os%map3dshift22d(-shift)
            shift       = - shift * smpd
            shifted_pdb = pdb
            call shifted_pdb%translate(shift)
            if( present(pdbout) ) call shifted_pdb%writepdb(pdbout) ! pdb output
            do i = 1, shifted_pdb%get_n()
                call self%set_within( shifted_pdb%get_coord(i), radius, 1.)
            enddo
            call shifted_pdb%kill
        else
            ! build mask
            do i = 1, pdb%get_n()
                call self%set_within( pdb%get_coord(i), radius, 1.)
            enddo
        endif
        call self%grow_bins(1)
        call distimg%mul(self) ! for suggested focusmsk
        call self%cos_edge(params_glob%edge,cos_img)
        ! multiply with mask
        !call vol_inout%ifft()
        call vol_inout%mul(cos_img)
        if( was_ft ) call vol_inout%fft()
        ! focusmsk
        minmax = distimg%minmax()
        write(logfhandle,'(A,I4)') '>>> SUGGESTED FOCUSMSK: ', ceiling(minmax(2))+params_glob%edge
        call distimg%kill
        call cos_img%kill
    end subroutine mask_from_pdb

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

    subroutine automask2D( imgs, ngrow, winsz, edge, diams, shifts, write2disk )
        use simple_segmentation
        use simple_image_bin, only: image_bin
        class(image),      intent(inout) :: imgs(:)
        integer,           intent(in)    :: ngrow, winsz, edge
        real, allocatable, intent(inout) :: diams(:), shifts(:,:)
        logical, optional, intent(in)    :: write2disk
        type(image_bin),    allocatable   :: img_bin(:), cc_img(:)
        real,              allocatable   :: ccsizes(:)
        integer :: i, n, loc, ldim(3)
        real    :: smpd, xyz(3)
        logical :: l_write
        n = size(imgs)
        l_write = .false.
        if( present(write2disk) ) l_write = write2disk
        l_write = l_write .and. params_glob%part.eq.1
        if( allocated(diams)  ) deallocate(diams)
        if( allocated(shifts) ) deallocate(shifts)
        ! allocate
        allocate(diams(n), shifts(n,2), source=0.)
        ldim = imgs(1)%get_ldim()
        smpd = imgs(1)%get_smpd()
        allocate(img_bin(n), cc_img(n))
        do i = 1,n
            call img_bin(i)%new_bimg(ldim, smpd, wthreads=.false.)
            call img_bin(i)%copy(imgs(i))
            call cc_img(i)%new_bimg( ldim, smpd, wthreads=.false.)
        end do
        if( trim(params_glob%automsk).eq.'tight' )then
            write(logfhandle,'(A)') '>>> 2D AUTOMASKING, TIGHT'
        else
            write(logfhandle,'(A)') '>>> 2D AUTOMASKING'
        endif
        !$omp parallel do default(shared) private(i,ccsizes,loc,xyz) schedule(static) proc_bind(close)
        do i = 1,n
            call img_bin(i)%zero_edgeavg
            ! dampens below zero (object positive in class averages/reprojs)
            call img_bin(i)%div_below(0.,10.)
            ! low-pass filter
            call img_bin(i)%bp(0., params_glob%amsklp)
            ! filter with non-local means
            call img_bin(i)%NLmean2D
            ! binarize with Otsu
            call otsu_img(img_bin(i), mskrad=params_glob%msk, positive=trim(params_glob%automsk).eq.'tight')
            call img_bin(i)%masscen(xyz)
            shifts(i,:) = xyz(:2)
            call img_bin(i)%set_imat
            ! grow ngrow layers
            if( ngrow > 0 ) call img_bin(i)%grow_bins(ngrow)
            ! find the largest connected component
            call img_bin(i)%find_ccs(cc_img(i))
            ccsizes = cc_img(i)%size_ccs()
            loc     = maxloc(ccsizes,dim=1)
            ! estimate its diameter
            call cc_img(i)%diameter_cc(loc, diams(i))
            if( diams(i) > 2.*(params_glob%msk+real(ngrow))*smpd )then
                ! incorrect component was chosen, fall back on spherical mask
                diams(i)    = 2.*(params_glob%msk-edge-COSMSKHALFWIDTH)*smpd
                shifts(i,:) = 0.
                call cc_img(i)%disc(ldim, smpd, diams(i)/(2.*smpd))
                call cc_img(i)%set_imat
            else
                ! turn it into a binary image for mask creation
                call cc_img(i)%cc2bin(loc)
                call cc_img(i)%masscen(xyz)
                shifts(i,:) = xyz(:2)
                ! median filter to smoothen
                if( winsz > 0 )then
                    call cc_img(i)%real_space_filter(winsz, 'median')
                    call cc_img(i)%set_imat
                endif
                ! fill-in holes
                call cc_img(i)%set_edgecc2background
            endif
            ! apply cosine egde to soften mask (to avoid Fourier artefacts)
            call imgs(i)%zero_and_unflag_ft
            call cc_img(i)%cos_edge(edge,imgs(i))
        end do
        !$omp end parallel do
        ! destruct
        do i = 1,n
            call img_bin(i)%kill_bimg
            call cc_img(i)%write(string('binarized_automask2D.mrc'), i)
            call cc_img(i)%kill_bimg
            call imgs(i)%write(string('masks_automask2D.mrc'), i)
        end do
        deallocate(img_bin, cc_img)
        if( allocated(ccsizes) ) deallocate(ccsizes)
    end subroutine automask2D

end module simple_image_msk
