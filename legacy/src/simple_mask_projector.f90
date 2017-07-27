!------------------------------------------------------------------------------!
! SIMPLE v2.5         Elmlund & Elmlund Lab          simplecryoem.com          !
!------------------------------------------------------------------------------!
!>  \brief  SIMPLE class  mask projector
module simple_mask_projector
use simple_defs       ! use all in there
use simple_image,     only: image
use simple_projector, only: projector
use simple_ori,       only: ori
use simple_oris,      only: oris
use simple_params,    only: params
use simple_jiffys,    only: alloc_err
implicit none

public :: mask_projector
private
#include "simple_local_flags.inc"

real,    parameter :: MSKWIDTH = 10.                     !< HALF SOFT MASK WIDTH (+/-)

type, extends(image) :: mask_projector
    private
    type(oris)               :: o_msk                    !< reference orientations
    real                     :: smpd_here = 0.           !< maximum circular mask
    real                     :: msk       = 0.           !< maximum circular mask
    real                     :: amsklp    = 0.           !< maximum circular mask
    real, allocatable        :: adamsks(:)               !< evaluated circular masks
    real                     :: mskwidth  = MSKWIDTH     !< +/- soft masking width
    real                     :: dens      = 0.
    real                     :: mw        = 0.
    integer                  :: edge      = 1
    integer                  :: binwidth  = 1
    integer                  :: n         = 0
    integer                  :: idim(3)   = 0            !< image dimension
    logical                  :: mskproj_exists = .false.
  contains
    ! CONSTRUCTORS
    procedure          :: init2D
    procedure          :: automask3D
    ! CALCULATORS
    procedure, private :: calc_adamsk
    ! 2D CALCULATORS
    procedure          :: update_cls
    procedure, private :: bin_cavg
    ! 3D CALCULATORS
    procedure, private :: bin_vol
    procedure, private :: env_rproject 
    ! MODIFIER
    procedure          :: apply_mask2D
    ! GETTERS
    procedure          :: get_adamsk
    procedure          :: get_msk
    procedure          :: get_imgmsk
    ! DESTRUCTOR
    procedure          :: kill_mskproj
end type mask_projector

contains

    ! CONSTRUCTOR

    !>  \brief  is a 2D constructor
    !>          on exit the parent image is untouched
    subroutine init2D(self, p, ncls)
        class(mask_projector),      intent(inout) :: self
        class(params),              intent(in)    :: p
        integer,                    intent(in)    :: ncls
        integer              :: alloc_stat
        if( .not.self%is_2d() )stop 'this routine is intended for 2D images only, simple_mask_projector::mskref'
        if(.not.self%even_dims())stop 'even dimensions assumed; simple_mask_projector::init2D'
        if(self%is_ft())         stop 'real space only; simple_mask_projector::init2D'
        call self%kill_mskproj
        self%idim      = [p%boxmatch, p%boxmatch, 1]
        self%n         = ncls
        self%amsklp    = p%amsklp
        self%msk       = p%msk
        self%edge      = p%edge
        self%dens      = p%dens
        self%binwidth  = p%binwidth
        self%mw        = p%mw
        self%smpd_here = p%smpd
        self%mskwidth  = min(self%mskwidth, real(minval(self%idim(:2)/2))-self%msk)
        if(self%mskwidth < 1.)stop 'incompatible dimensiosn in simple_mask_projector%init_parms'
        allocate(self%adamsks(self%n), stat=alloc_stat )
        call alloc_err('in simple_mask_projector::init2D 1', alloc_stat)
        self%adamsks = 0.
        self%mskproj_exists = .true.
        DebugPrint'simple_mask_projector::init2D done'
    end subroutine init2D

    !>  \brief  is a 3D constructor and modifier
    !>  On output the parent volume is the enveloppe mask
    !>  The returned volume is enveloppe masked.
    subroutine automask3D( self, vol_inout, amsklp, mw, binwidth, edge, dens )
        class(mask_projector), intent(inout) :: self
        class(image),          intent(inout) :: vol_inout
        real,                  intent(in)    :: amsklp, mw, dens
        integer,               intent(in)    :: binwidth, edge
        logical :: was_ft
        if( vol_inout%is_2d() )stop 'automask3D is intended for volumes only, simple_mask_projector%init_mskproj'
        call self%kill_mskproj
        self%amsklp    = amsklp
        self%edge      = edge
        self%dens      = dens
        self%binwidth  = binwidth
        self%mw        = mw
        was_ft = vol_inout%is_ft()
        if( vol_inout%is_ft() )call vol_inout%bwd_ft
        self = vol_inout
        ! binarize volume
        call self%bin_vol
        ! add volume soft edge
        call self%cos_edge(self%edge)
        ! apply mask to volume
        call vol_inout%mul(self)
        ! the end
        if( was_ft )call vol_inout%fwd_ft
        DebugPrint 'simple_mask_projector::automask2D done'
    end subroutine automask3D

    ! CALCULATORS

    !>  \brief  is for getting the adaptive circular mask
    real function calc_adamsk( self, img_msk )result( new_msk )
        class(mask_projector), intent(inout) :: self
        class(image),          intent(inout) :: img_msk
        type(image) :: img_dist, tmp_img
        real        :: minmax(2)
        tmp_img = img_msk
        call img_dist%new(self%idim, self%get_smpd())
        call img_dist%cendist
        ! multiply enveloppe mask with square distance matrix
        call tmp_img%mul(img_dist)
        ! determine circular mask size
        minmax  = tmp_img%minmax()
        new_msk = real( ceiling(minmax(2) + self%mskwidth) )
        new_msk = min(new_msk, self%msk)
        DebugPrint'simple_mask_projector::calc_adamsk done'
    end function calc_adamsk

    ! 2D CALCULATORS

    !>  \brief  is for enveloppe masking the input image
    subroutine update_cls( self, ref, cls )
        class(mask_projector), intent(inout) :: self
        class(image),          intent(inout) :: ref
        integer,               intent(in)    :: cls
        type(image) :: img
        if( .not.self%mskproj_exists )stop 'mask_projector object has not be instantiated'
        if( cls > self%n )stop 'class index out of range'
        ! binarize image
        img = ref
        call self%bin_cavg(img)
        ! adaptive circular mask
        self%adamsks(cls) = self%calc_adamsk(img)
        ! soft edge mask
        call img%cos_edge(self%edge)
        call img%norm('sigm')
        ! apply enveloppe mask to reference
        call ref%mul(img)
        DebugPrint'simple_mask_projector::update_cls done'
    end subroutine update_cls

    !>  \brief  is for binarizing the 2D image
    subroutine bin_cavg( self, img )
        class(mask_projector), intent(inout) :: self
        class(image),          intent(inout) :: img
        type(image) :: img_pad, img_copy
        integer     :: ldim_pad(3)
        ! init
        ldim_pad(1:2) = self%idim(1:2)*2
        ldim_pad(3)   = 1
        call img_pad%new(ldim_pad,  self%get_smpd())
        call img_copy%copy(img)
        ! normalize
        call img_copy%norm()
        ! soft masking
        call img_copy%mask(self%msk, 'soft')
        ! pad
        call img_copy%pad(img_pad) ! need be padded?
        ! low-pass
        call img_pad%fwd_ft
        call img_pad%bp(0., self%amsklp)
        call img_pad%bwd_ft
        ! binarize within mask
        call img_pad%bin('msk', self%msk)
        ! add one layer
        call img_pad%grow_bins(1)
        ! clip
        call img_pad%clip(img)
        ! clean
        call img_copy%kill
        call img_pad%kill
        DebugPrint'simple_mask_projector::bin_cavg done'
    end subroutine bin_cavg

    ! 3D CALCULATORS

    subroutine bin_vol( self )
        use simple_math, only: nvoxfind
        class(mask_projector), intent(inout) :: self
        integer :: nnvox
        call self%bp(0., self%amsklp)
        ! find nr of voxels corresponding to mw
        if( self%dens > 0. )then
            nnvox = nvoxfind(self%get_smpd(), self%mw, self%dens)
        else
            nnvox = nvoxfind(self%get_smpd(), self%mw)
        endif
        nnvox = nint(1.1*real(nnvox))   ! this is to compensate for the low-pass filter
        ! binarize
        call self%bin(nnvox)
        ! binary layers
        call self%grow_bins(self%binwidth)
        call self%norm_bin
        DebugPrint'simple_mask_projector::bin_vol done'
    end subroutine bin_vol

    !>  \brief  volume mask projector
    subroutine env_rproject(self, e, img)
        !$ use omp_lib
        !$ use omp_lib_kinds
        class(mask_projector), intent(inout) :: self   !< projector instance
        class(ori),            intent(inout) :: e      !< Euler angle
        type(image),           intent(inout) :: img    !< resulting projection image
        real, allocatable :: rmat(:,:,:)
        real              :: out_coos(3), maxrad, rad(3), thresh
        real              :: incr_k(3), rvec(3), rvec_k(3)
        integer           :: orig(3), i, j, k, sqmaxrad, vec(3)
        ! init
        thresh   = 0.9999             ! prior soft masking is discarded
        img      = 0.
        orig     = self%idim/2+1
        maxrad   = min(self%msk, real(minval(self%idim(1:2)))/2.-1.)
        sqmaxrad = nint(maxrad**2)
        out_coos = 0.
        rad      = 0.
        vec      = 0
        rmat     = self%get_rmat()
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
                    if(rad(3) > sqmaxrad)cycle
                    !vec    = orig + floor( matmul(out_coos,e%get_mat()) ) ! the old way
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
        DebugPrint'simple_mask_projector::env_rproject done'
    end subroutine env_rproject

    ! MODIFIER

    !> \brief  applies the circular/envelope mask
    subroutine apply_mask2D( self, img, cls )
        class(mask_projector), intent(inout) :: self
        class(image),          intent(inout) :: img
        integer,               intent(in)    :: cls
        integer :: ldim(3)
        ldim = img%get_ldim()
        if( .not.img%is_2d())stop 'only for images; simple_mask_projector::apply_mask2D'
        if( any(ldim(:2)-self%idim(:2).ne.0) )stop 'Incompatible dimensions; simple_mask_projector::apply_mask2D'
        if( .not.self%is_2d() )stop 'erroneous function call; simple_mask_projector::apply_mask2D'
        if( cls > self%n )stop 'class index out of range; simple_mask_projector::apply_mask2D'
        call img%mask(self%adamsks(cls), 'soft')
        DebugPrint'simple_mask_projector::apply_mask2D done'
    end subroutine apply_mask2D

    ! GETTERS

    real function get_msk( self )
        class(mask_projector), intent(inout) :: self
        get_msk = self%msk
    end function

    real function get_adamsk( self, i )
        class(mask_projector), intent(inout) :: self
        integer, intent(in) :: i
        if( .not.allocated(self%adamsks) )stop 'adamask has not been calculated; simple_mask_projector%get_adamsk'
        if( i > self%n )stop 'index out of range; simple_mask_projector%get_adamsk'
        get_adamsk = self%adamsks(i)
    end function

    function get_imgmsk( self, i )result( img_msk )
        class(mask_projector), intent(inout) :: self
        integer,               intent(in)    :: i
        type(image) :: img_msk
        if( i > self%n )stop 'index out of range; simple_mask_projector%get_imgmsk'
        call img_msk%new([self%idim(1),self%idim(2),1], self%smpd_here)
        img_msk = 1.
        call img_msk%mask(self%adamsks(i), 'soft')
    end function

    ! DESTRUCTORS

    !>  \brief  is the destructor
    subroutine kill_mskproj( self )
        class(mask_projector), intent(inout) :: self
        self%mskwidth = MSKWIDTH
        self%n        = 0
        self%msk      = 0.
        self%mskproj_exists = .false.
    end subroutine kill_mskproj

end module simple_mask_projector
