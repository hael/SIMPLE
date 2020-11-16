! 2D/3D envelope and adaptive masking
module simple_masker
include 'simple_lib.f08'
use simple_image,      only: image
use simple_binimage,   only: binimage
use simple_parameters, only: params_glob
implicit none

public :: masker, automask2D
private
#include "simple_local_flags.inc"

integer, parameter :: WINSZ = 3 !< real-space filter half-width window
logical, parameter :: DEBUG = .false.

type, extends(binimage) :: masker
    private
    real    :: msk       = 0.   !< maximum circular mask
    real    :: amsklp    = 0.   !< maximum circular mask
    real    :: mw        = 0.   !< moleclular weight (in kDa)
    real    :: pix_thres = 0.   !< binarisation threshold
    integer :: edge      = 3    !< edge width
    integer :: binwidth  = 1    !< additional layers to grow
    integer :: n         = 0    !< number of classes
    integer :: idim(3)   = 0    !< image dimension
  contains
    procedure          :: automask3D
    procedure          :: mask_from_pdb
    procedure, private :: bin_vol_thres
    procedure, private :: env_rproject
end type masker

contains

    !>  \brief  is for 3D automasking. On output the parent volume is the envelope mask
    !!          The returned volume is envelope masked.
    subroutine automask3D( self, vol_inout )
        class(masker), intent(inout) :: self
        class(image),  intent(inout) :: vol_inout
        logical     :: was_ft
        if( vol_inout%is_2d() )THROW_HARD('automask3D is intended for volumes only; automask3D')
        self%msk       = params_glob%msk
        self%amsklp    = params_glob%amsklp
        self%mw        = params_glob%mw
        self%binwidth  = params_glob%binwidth
        self%edge      = params_glob%edge
        self%pix_thres = params_glob%thres
        write(logfhandle,'(A,F7.1,A)') '>>> AUTOMASK LOW-PASS:           ', self%amsklp,  ' ANGSTROMS'
        write(logfhandle,'(A,I7,A)'  ) '>>> AUTOMASK SOFT EDGE WIDTH:    ', self%edge,    ' PIXEL(S)'
        write(logfhandle,'(A,I7,A)'  ) '>>> AUTOMASK BINARY LAYERS WIDTH:', self%binwidth,' PIXEL(S)'
        write(logfhandle,'(A,F7.1,A)') '>>> AUTOMASK MOLECULAR WEIGHT:   ', self%mw,      ' kDa'
        was_ft = vol_inout%is_ft()
        if( was_ft ) call vol_inout%ifft()
        call self%transfer2bimg(vol_inout)
        ! binarize volume
        call self%bin_vol_thres
        ! add volume soft edge
        call self%cos_edge(self%edge)
        ! apply mask to volume
        call vol_inout%zero_background()
        call vol_inout%mul(self)
        ! the end
        if( was_ft )call vol_inout%fft()
    end subroutine automask3D

    subroutine mask_from_pdb( self,  pdb, vol_inout, os, pdbout)
        use simple_oris,  only: oris
        use simple_atoms, only: atoms
        class(masker),              intent(inout) :: self
        type(atoms),                intent(inout) :: pdb
        class(image),               intent(inout) :: vol_inout
        class(oris),      optional, intent(inout) :: os
        character(len=*), optional, intent(inout) :: pdbout
        type(image) :: distimg, cos_img
        type(atoms) :: shifted_pdb
        real        :: centre(3), shift(3), pdb_center(3), minmax(2), radius, smpd
        integer     :: i
        logical     :: was_ft
        if( vol_inout%is_2d() )THROW_HARD('intended for volumes only; mask_from_pdb')
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
        if( params_glob%center.eq.'yes' )then
            ! shift volume & oris
            call vol_inout%shift(shift)
            if(present(os)) call os%map3dshift22d(-shift)
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
        call vol_inout%ifft()
        call vol_inout%mul(cos_img)
        if(was_ft) call vol_inout%fft()
        ! focusmsk
        minmax = distimg%minmax()
        write(logfhandle,'(A,I4)') '>>> SUGGESTED FOCUSMSK: ', ceiling(minmax(2))+params_glob%edge
        call distimg%kill
        call cos_img%kill
    end subroutine mask_from_pdb

    ! BINARISATION ROUTINES

    !>  \brief  is for binarizing the 3D image using thresholding
    subroutine bin_vol_thres( self )
        class(masker), intent(inout) :: self
        integer :: nnvox
        call self%zero_below(self%pix_thres)
        call self%real_space_filter( WINSZ, 'average')
        call self%bp(0., self%amsklp)
        ! find nr of voxels corresponding to mw
        nnvox = nvoxfind(self%get_smpd(), self%mw)
        ! binarize again
        call self%binarize(nnvox)
        ! binary layers
        call self%set_imat
        call self%grow_bins(self%binwidth)
    end subroutine bin_vol_thres

    ! CALCULATORS

    !>  \brief  volume mask projector
    subroutine env_rproject(self, e, img)
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_ori,    only: ori
        class(masker), intent(inout) :: self   !< projector instance
        class(ori),    intent(inout) :: e      !< Euler angle
        type(image),   intent(inout) :: img    !< resulting projection image
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

    subroutine automask2D( imgs, mask, ngrow, winsz, edge, diams, write2disk )
        use simple_segmentation
        use simple_binimage, only: binimage
        class(image),      intent(inout) :: imgs(:)
        logical,           intent(in)    :: mask(:)
        integer,           intent(in)    :: ngrow, winsz, edge
        real, allocatable, intent(inout) :: diams(:)
        logical, optional, intent(in)    :: write2disk
        real, allocatable :: ccsizes(:)
        type(binimage) :: img_bin, cc_img
        type(image)    :: cos_img
        integer        :: i, n, loc(1)
        real           :: thresh(3)
        logical        :: l_write
        n = size(imgs)
        if( size(mask) /= n ) THROW_HARD('mask array size does not conform with image array; automask2D')
        l_write = .true.
        if( present(write2disk) ) l_write = write2disk
        if( allocated(diams) ) deallocate(diams)
        allocate(diams(n), source=0.)
        if( l_write ) write(logfhandle,'(A)') '>>> 2D AUTOMASKING'
        do i=1,n
            call progress(i, n)
            call img_bin%copy(imgs(i))
            if( mask(i) )then
                ! filter with non-local means
                call img_bin%NLmean
                ! binarise with robust Otsu
                call otsu_robust_fast(img_bin, is2D=.true., noneg=.true., thresh=thresh)
                if( l_write ) call img_bin%write(BIN_OTSU, i)
                ! find the largest connected component
                call img_bin%find_ccs(cc_img)
                ccsizes = cc_img%size_ccs()
                loc = maxloc(ccsizes)
                ! estimate its diameter
                call cc_img%diameter_cc(loc(1), diams(i))
                ! turn it into a binary image for mask creation
                call cc_img%cc2bin(loc(1))
                ! grow ngrow layers
                call img_bin%grow_bins(ngrow)
                if( l_write ) call cc_img%write(BIN_OTSU_GROW, i)
                ! median filter to smoothen
                if( mask(i) ) call cc_img%real_space_filter(winsz, 'median')
                if( l_write ) call cc_img%write(BIN_OTSU_GROW_MED, i)
                ! apply cosine egde to soften mask (to avoid Fourier artefacts)
                call cc_img%cos_edge(edge,cos_img)
                if( l_write ) call cos_img%write(MSK_OTSU, i)
                ! zero negative values before applyting the mask
                call imgs(i)%zero_neg
                ! call imgs(i)%remove_neg
                ! apply
                call imgs(i)%mul(cos_img)
                if( l_write ) call imgs(i)%write(AMSK_OTSU, i)
            else
                call img_bin%zero
            endif
        end do
        call img_bin%kill
        call cc_img%kill
        call cos_img%kill
        if( allocated(ccsizes) ) deallocate(ccsizes)
    end subroutine automask2D

end module simple_masker
