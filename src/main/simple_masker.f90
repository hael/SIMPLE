! 2D/3D envelope and adaptive masking
module simple_masker
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image,      only: image
use simple_binimage,   only: binimage
use simple_parameters, only: params_glob
use simple_segmentation
implicit none

public :: masker, automask2D
private
#include "simple_local_flags.inc"

logical, parameter :: DEBUG         = .false.
logical, parameter :: L_WRITE       = .true.

type, extends(binimage) :: masker
    private
    real    :: msk           = 0.   !< maximum circular mask
    real    :: amsklp        = 0.   !< low-pass limit
    real    :: amsklp_prelim = 0.   !< low-pass limit
    real    :: pix_thres     = 0.   !< binarisation threshold
    integer :: edge          = 6    !< edge width
    integer :: binwidth      = 1    !< additional layers to grow
    integer :: idim(3)       = 0    !< image dimension
  contains
    procedure          :: automask3D
    procedure          :: automask3D_icm
    procedure          :: automask3D_otsu
    procedure, private :: automask3D_otsu_priv
    procedure          :: mask_from_pdb
    procedure, private :: env_rproject
end type masker

contains

    !>  \brief  is for 3D automasking. On output the parent volume is the envelope mask
    !!          The returned volume is envelope masked.
    subroutine automask3D( self, vol_inout )
        class(masker), intent(inout) :: self
        class(image),  intent(inout) :: vol_inout
        logical        :: was_ft
        type(binimage) :: msk_prelim
        if( vol_inout%is_2d() )THROW_HARD('automask3D is intended for volumes only; automask3D')
        self%amsklp_prelim = params_glob%amsklp_prelim
        self%amsklp        = params_glob%amsklp
        self%binwidth      = params_glob%binwidth
        self%edge          = params_glob%edge
        self%pix_thres     = params_glob%thres
        write(logfhandle,'(A,F7.1,A)') '>>> AUTOMASK LOW-PASS:           ', self%amsklp,  ' ANGSTROMS'
        write(logfhandle,'(A,I7,A)'  ) '>>> AUTOMASK SOFT EDGE WIDTH:    ', self%edge,    ' PIXEL(S)'
        write(logfhandle,'(A,I7,A)'  ) '>>> AUTOMASK BINARY LAYERS WIDTH:', self%binwidth,' PIXEL(S)'
        was_ft = vol_inout%is_ft()
        if( was_ft ) call vol_inout%ifft()
        call self%transfer2bimg(vol_inout)
        ! preliminary masking
        call msk_prelim%copy(self)
        call msk_prelim%binarize(self%pix_thres)
        call msk_prelim%grow_bins(self%binwidth)
        call msk_prelim%cos_edge(self%edge)
        if( L_WRITE ) call msk_prelim%write('msk_prelim.mrc')
        call self%mul(msk_prelim)
        if( L_WRITE ) call self%write('pre-masked.mrc')
        ! automasking
        call self%automask3D_otsu_priv(l_tight=.false., amsklp=self%amsklp)
        ! apply mask to volume
        call vol_inout%zero_background()
        call vol_inout%mul(self)
        call msk_prelim%kill_bimg
        if( was_ft ) call vol_inout%fft()
    end subroutine automask3D

    subroutine automask3D_icm( self, vol_even, vol_odd, vol_masked )
        class(masker), intent(inout) :: self
        class(image),  intent(inout) :: vol_even, vol_odd, vol_masked
        real, allocatable  :: fsc(:), filt(:)
        type(image)        :: vol_avg
        integer            :: ldim(3), filtsz
        real               :: msk
        real,    parameter :: LAM = 100.
        logical, parameter :: L_NLMEAN = .false. 
        if( vol_even%is_2d() )THROW_HARD('automask3D_icm is intended for volumes only; automask3D')
        self%amsklp   = params_glob%amsklp
        self%binwidth = params_glob%binwidth
        self%edge     = params_glob%edge
        ldim   = vol_even%get_ldim()
        msk    = real(ldim(1) / 2) - COSMSKHALFWIDTH - 1.
        filtsz = fdim(ldim(1)) - 1
        call vol_masked%copy(vol_even)
        call vol_masked%add(vol_odd)
        call vol_masked%mul(0.5)
        ! zero mean of outer pixels
        call vol_even%zero_background
        call vol_odd%zero_background
        ! spherical mask
        call vol_even%mask(msk, 'soft', backgr=0.)
        call vol_odd%mask(msk,  'soft', backgr=0.)
        ! calculate FSC
        allocate(fsc(filtsz), filt(filtsz), source=0.)
        call vol_even%fft()
        call vol_odd%fft()
        call vol_even%fsc(vol_odd, fsc)
        ! calculate a filter to be applied to the individual e/o pairs
        where( fsc > 0.        ) filt = fsc / (fsc + 1.)
        where( filt  > 0.99999 ) filt = 0.99999
        call vol_even%apply_filter(filt)
        call vol_odd%apply_filter(filt)
        ! put back in real-space
        call vol_even%ifft()
        call vol_odd%ifft()
        ! ICM filter
        call vol_even%ICM3D_eo(vol_odd, LAM)
        call vol_avg%copy(vol_even)
        call vol_avg%add(vol_odd)
        call vol_avg%mul(0.5)
        if( L_WRITE ) call vol_avg%write('ICM_avg.mrc')
        ! NLMEAN filter
        if( L_NLMEAN )then
            call vol_even%NLmean3D_eo(vol_odd, vol_avg)
            if( L_WRITE ) call vol_avg%write('NLmean3D_eo.mrc')
        else
            call vol_avg%bp(0., self%amsklp)
        endif
        call self%transfer2bimg(vol_avg)
        ! automasking
        call self%automask3D_otsu_priv(l_tight=.true.)
        ! apply mask to volume
        call vol_masked%zero_background()
        call vol_masked%mul(self)
        ! destruct
        call vol_avg%kill
    end subroutine automask3D_icm

    subroutine automask3D_otsu( self, vol_inout, do_apply )
        class(masker),     intent(inout) :: self
        class(image),      intent(inout) :: vol_inout
        logical, optional, intent(in)    :: do_apply
        type(masker) :: msk_prelim
        logical      :: was_ft, ddo_apply
        ddo_apply = .true.
        if( present(do_apply) ) ddo_apply = do_apply
        if( vol_inout%is_2d() )THROW_HARD('automask3D_otsu is intended for volumes only; automask3D')
        self%amsklp_prelim = params_glob%amsklp_prelim
        self%amsklp        = params_glob%amsklp
        self%binwidth      = params_glob%binwidth
        self%edge          = params_glob%edge
        self%pix_thres     = params_glob%thres
        write(logfhandle,'(A,F7.1,A)') '>>> AUTOMASK FIRST  LOW-PASS:    ', self%amsklp_prelim, ' ANGSTROMS'
        write(logfhandle,'(A,F7.1,A)') '>>> AUTOMASK SECOND LOW-PASS:    ', self%amsklp,        ' ANGSTROMS'
        write(logfhandle,'(A,I7,A)'  ) '>>> AUTOMASK SOFT EDGE WIDTH:    ', self%edge,          ' PIXEL(S)'
        write(logfhandle,'(A,I7,A)'  ) '>>> AUTOMASK BINARY LAYERS WIDTH:', self%binwidth,      ' PIXEL(S)'
        was_ft = vol_inout%is_ft()
        if( was_ft ) call vol_inout%ifft()
        call self%transfer2bimg(vol_inout)
        ! preliminary masking
        call msk_prelim%transfer2bimg(vol_inout)
        call msk_prelim%automask3D_otsu_priv(l_tight=.true., amsklp=self%amsklp_prelim)
        call self%mul(msk_prelim)
        ! automasking
        call self%automask3D_otsu_priv(l_tight=.false., amsklp=self%amsklp)
        if( ddo_apply )then
            ! apply mask to volume
            call vol_inout%zero_background()
            call vol_inout%mul(self)
        endif
        call msk_prelim%kill_bimg
        if( was_ft ) call vol_inout%fft()
    end subroutine automask3D_otsu

    subroutine automask3D_otsu_priv( self, l_tight, amsklp )
        class(masker),  intent(inout) :: self
        logical,        intent(in)    :: l_tight
        real, optional, intent(in)    :: amsklp
        real,    allocatable :: ccsizes(:)
        integer, allocatable :: imat_cc(:,:,:)
        type(binimage)       :: ccimage
        integer              :: loc(1), imax, sz, nccs
        if( present(amsklp) )then
            ! low-pass filter volume
            call self%bp(0., amsklp)
            if( L_WRITE ) call self%write('lped.mrc')
        endif
        ! binarize volume
        call otsu_img(self, tight=l_tight)
        call self%set_imat
        if( L_WRITE ) call self%write('binarized.mrc')
        ! identify connected components
        call self%find_ccs(ccimage, update_imat=.true.)
        ! extract all cc sizes (in # pixels)
        ccsizes = ccimage%size_ccs()
        sz      = size(ccsizes)
        write(logfhandle,'(A,I7,A)'  ) '>>> FOUND:                       ', sz,   ' CONNECTED COMPONENT(S)'
        if( sz > 1 )then
            loc = maxloc(ccsizes)
            call ccimage%cc2bin(loc(1))
        else
            call ccimage%cc2bin(1)
        endif
        call self%copy_bimg(ccimage)
        if( L_WRITE ) call self%write('largest_cc.mrc')
        ! add layers
        call self%grow_bins(self%binwidth)
        ! add volume soft edge
        call self%cos_edge(self%edge)
        ! destruct
        call ccimage%kill_bimg
    end subroutine automask3D_otsu_priv

    subroutine mask_from_pdb( self,  pdb, vol_inout, os, pdbout )
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

    ! CALCULATORS

    !>  \brief  volume mask projector
    subroutine env_rproject(self, e, img)
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

    subroutine automask2D( imgs, ngrow, winsz, edge, diams, shifts, write2disk )
        use simple_segmentation
        use simple_binimage, only: binimage
        class(image),      intent(inout) :: imgs(:)
        integer,           intent(in)    :: ngrow, winsz, edge
        real, allocatable, intent(inout) :: diams(:), shifts(:,:)
        logical, optional, intent(in)    :: write2disk
        type(binimage),    allocatable   :: img_bin(:), cc_img(:)
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
            call img_bin(i)%NLmean
            ! if( l_write ) call img_bin(i)%write('NLmean_filtered.mrc', i)
            ! binarize with Otsu
            call otsu_img(img_bin(i), mskrad=params_glob%msk, positive=trim(params_glob%automsk).eq.'tight')
            call img_bin(i)%masscen(xyz)
            shifts(i,:) = xyz(:2)
            call img_bin(i)%set_imat
            ! if( l_write ) call img_bin(i)%write(BIN_OTSU, i)
            ! grow ngrow layers
            if( ngrow > 0 ) call img_bin(i)%grow_bins(ngrow)
            ! if( l_write ) call img_bin(i)%write(BIN_OTSU_GROWN, i)
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
                    ! if( l_write ) call cc_img(i)%write(BIN_OTSU_MED, i)
                endif
                ! fill-in holes
                call cc_img(i)%fill_holes
                ! if( l_write ) call cc_img(i)%write(BIN_OTSU_HOLES_FILL, i)
            endif
            ! apply cosine egde to soften mask (to avoid Fourier artefacts)
            call imgs(i)%zero_and_unflag_ft
            call cc_img(i)%cos_edge(edge,imgs(i))
        end do
        !$omp end parallel do
        ! destruct
        do i = 1,n
            call img_bin(i)%kill_bimg
            call cc_img(i)%write('binarized_automask2D.mrc', i)
            call cc_img(i)%kill_bimg
            call imgs(i)%write('masks_automask2D.mrc', i)
        end do
        deallocate(img_bin, cc_img)
        if( allocated(ccsizes) ) deallocate(ccsizes)
    end subroutine automask2D

end module simple_masker
