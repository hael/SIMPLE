! 2D/3D envelope and adaptive masking
module simple_masker
include 'simple_lib.f08'

use simple_image,  only: image
use simple_ori,    only: ori
use simple_params, only: params
implicit none

public :: masker
private
#include "simple_local_flags.inc"

integer, parameter :: WINSZ = 3 !< real-space filter half-width window

!> image child type
type, extends(image) :: masker
    private
    real              :: msk           = 0.   !< maximum circular mask
    real              :: amsklp        = 0.   !< maximum circular mask
    real              :: mw            = 0.   !< moleclular weight (in kDa)
    real              :: pix_thres     = 0.   !< binarisation threshold
    integer           :: edge          = 3    !< edge width
    integer           :: binwidth      = 1    !< additional layers to grow
    integer           :: n             = 0    !< number of classes
    integer           :: idim(3)       = 0    !< image dimension
  contains
    procedure          :: automask3D
    procedure          :: resmask
    procedure          :: apply_2Denvmask22Dref
    procedure          :: mask_from_pdb
    procedure, private :: bin_cavg
    procedure, private :: bin_vol_thres
    procedure, private :: env_rproject
end type masker

contains

    !>  \brief  is for 3D automasking. On output the parent volume is the envelope mask
    !!          The returned volume is envelope masked.
    subroutine automask3D( self, p, vol_inout )
        class(masker), intent(inout) :: self
        class(params), intent(in)    :: p
        class(image),  intent(inout) :: vol_inout
        logical :: was_ft
        if( vol_inout%is_2d() )stop 'automask3D is intended for volumes only, simple_masker::automask3D'
        self%msk       = p%msk
        self%amsklp    = p%amsklp
        self%mw        = p%mw
        self%binwidth  = p%binwidth
        self%edge      = p%edge
        self%pix_thres = p%thres
        write(*,'(A,F7.1,A)') '>>> AUTOMASK LOW-PASS:           ', self%amsklp,  ' ANGSTROMS'
        write(*,'(A,I7,A)'  ) '>>> AUTOMASK SOFT EDGE WIDTH:    ', self%edge,    ' PIXEL(S)'
        write(*,'(A,I7,A)'  ) '>>> AUTOMASK BINARY LAYERS WIDTH:', self%binwidth,' PIXEL(S)'
        write(*,'(A,F7.1,A)') '>>> AUTOMASK MOLECULAR WEIGHT:   ', self%mw,      ' kDa'
        was_ft = vol_inout%is_ft()
        if( vol_inout%is_ft() )call vol_inout%bwd_ft
        call self%copy(vol_inout)
        ! binarize volume
        call self%bin_vol_thres
        ! add volume soft edge
        call self%cos_edge(self%edge)
        ! apply mask to volume call vol_inout%norm()
        call vol_inout%zero_background()
        call vol_inout%mul(self)
        ! the end
        if( was_ft )call vol_inout%fwd_ft
    end subroutine automask3D

    !>  \brief  envelope mask for resolution estimation
    !!          it is assumed that the envelope mask (auotmask) is set from the start
    subroutine resmask( self, p )
        class(masker), intent(inout) :: self
        class(params), intent(in)    :: p
        type(image) :: mskimg, distimg
        integer     :: winsz
        real        :: ave, sdev, maxv, minv, med
        ! create edge-less envelope
        call self%remove_edge
        ! calculate distance stats
        call distimg%new(self%get_ldim(), p%smpd)
        call distimg%cendist
        call distimg%stats(self, ave, sdev, maxv, minv, med)
        winsz = max(8, nint((p%msk - maxv) / 2.))
        ! soft edge mask
        call self%grow_bins(2) ! in addition to binwidth (4 safety)
        call self%cos_edge(winsz)
        ! mask with spherical sof mask
        call self%mask(p%msk, 'soft')
    end subroutine resmask

    !>  \brief  is for envelope masking of the reference in prime2D
    subroutine apply_2Denvmask22Dref( self, ref )
        class(masker), intent(inout) :: self
        class(image),  intent(inout) :: ref
        type(image) :: img
        ! binarize image
        img = ref
        call self%bin_cavg(img)
        ! soft edge mask
        call img%cos_edge(self%edge)
        ! apply envelope mask to reference
        call ref%mul(img)
        if( DEBUG )write(*,*)'simple_masker::update_cls done'
    end subroutine apply_2Denvmask22Dref

    !>  \brief  is for
    subroutine mask_from_pdb( self, p, pdb, vol_inout, os, pdbout)
        use simple_oris,  only: oris
        use simple_atoms, only: atoms
        class(masker),              intent(inout) :: self
        class(params),              intent(in)    :: p
        type(atoms),                intent(inout) :: pdb
        class(image),               intent(inout) :: vol_inout
        class(oris),      optional, intent(inout) :: os
        character(len=*), optional, intent(inout) :: pdbout
        type(image) :: distimg
        type(atoms) :: shifted_pdb
        real        :: centre(3), shift(3), pdb_center(3), xyz(3), minmax(2), radius, smpd
        integer     :: i
        logical     :: was_ft
        if( vol_inout%is_2d() )stop 'sphere_mask is intended for volumes only, simple_masker::sphere_mask'
        was_ft = vol_inout%is_ft()
        smpd   = vol_inout%get_smpd()
        call self%new(vol_inout%get_ldim(), smpd)
        call distimg%new(vol_inout%get_ldim(), smpd)
        call distimg%cendist()
        pdb_center = pdb%get_geom_center()
        centre     = real(self%get_ldim()-1)/2. * smpd
        shift      = ( pdb_center - centre ) / smpd
        if( p%binwidth == 0 ) then
            radius = smpd
        else
            radius = real(p%binwidth) * smpd
        endif
        if( p%center.eq.'yes' )then
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
        call self%grow_bin()
        call distimg%mul(self) ! for suggested focusmsk
        call self%cos_edge(p%edge)
        ! multiply with mask
        call vol_inout%bwd_ft()
        call vol_inout%mul(self)
        if(was_ft) call vol_inout%fwd_ft
        ! focusmsk
        minmax = distimg%minmax()
        write(*,'(A,I4)') '>>> SUGGESTED FOCUSMSK: ', ceiling(minmax(2))+p%edge
        call distimg%kill
    end subroutine mask_from_pdb

    ! BINARISATION ROUTINES

    !>  \brief  is for binarizing the 2D image
    subroutine bin_cavg( self, img )
        class(masker), intent(inout) :: self
        class(image),  intent(inout) :: img
        ! normalize
        call img%norm()
        ! soft masking
        call img%mask(self%msk, 'soft')
        ! low-pass
        call img%bp(0., self%amsklp)
        ! binarize within mask
        call img%mask(self%msk, 'hard')
        call img%bin_kmeans
        ! add one layer
        call img%grow_bins(self%binwidth)
        if( DEBUG ) write(*,*)'simple_masker::bin_cavg done'
    end subroutine bin_cavg

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
        call self%bin(nnvox)
        ! binary layers
        call self%grow_bins(self%binwidth)
        if( DEBUG )write(*,*)'simple_masker::bin_vol done'
    end subroutine bin_vol_thres

    ! CALCULATORS

    !>  \brief  volume mask projector
    subroutine env_rproject(self, e, img)
        !$ use omp_lib
        !$ use omp_lib_kinds
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
            print *, 'maxrad:       ', maxrad
            print *, 'sqmaxrad:     ', sqmaxrad
            print *, 'maxval(rmat): ', maxval(rmat)
            print *, 'minval(rmat): ', minval(rmat)
            print *, 'self%idim:    ', self%idim
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
        if( DEBUG )write(*,*)'simple_masker::env_rproject done'
    end subroutine env_rproject

end module simple_masker
