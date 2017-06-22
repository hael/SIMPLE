!>  \brief  SIMPLE projector class
module simple_masker
use simple_defs       ! use all in there
use simple_image,     only: image
use simple_projector, only: projector
use simple_ori,       only: ori
use simple_oris,      only: oris
use simple_params,    only: params
use simple_jiffys,    only: alloc_err
implicit none

public :: masker
private

real,    parameter :: MSKWIDTH = 10.                     !< HALF SOFT MASK WIDTH (+/-)
logical, parameter :: DEBUG    = .false.

type, extends(image) :: masker
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
    logical                  :: masker_exists = .false.
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
    procedure          :: kill_masker
end type masker

contains

    ! CONSTRUCTOR

    !>  \brief  is a 2D constructor
    !>          on exit the parent image is untouched
    subroutine init2D(self, p, ncls)
        class(masker),      intent(inout) :: self
        class(params),              intent(in)    :: p
        integer,                    intent(in)    :: ncls
        integer              :: alloc_stat
        if( .not.self%is_2d() )stop 'this routine is intended for 2D images only, simple_masker::mskref'
        if(.not.self%even_dims())stop 'even dimensions assumed; simple_masker::init2D'
        if(self%is_ft())         stop 'real space only; simple_masker::init2D'
        call self%kill_masker
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
        if(self%mskwidth < 1.)stop 'incompatible dimensiosn in simple_masker%init_parms'
        allocate(self%adamsks(self%n), stat=alloc_stat )
        call alloc_err('in simple_masker::init2D 1', alloc_stat)
        self%adamsks = 0.
        self%masker_exists = .true.
        if( DEBUG )write(*,*)'simple_masker::init2D done'
    end subroutine init2D

    !>  \brief  is a 3D constructor and modifier
    !>  On output the parent volume is the enveloppe mask
    !>  The returned volume is enveloppe masked.
    subroutine automask3D( self, vol_inout, msk, amsklp, mw, binwidth, edge, dens )
        class(masker), intent(inout) :: self
        class(image),          intent(inout) :: vol_inout
        real,                  intent(in)    :: msk, amsklp, mw, dens
        integer,               intent(in)    :: binwidth, edge
        logical :: was_ft
        if( vol_inout%is_2d() )stop 'automask3D is intended for volumes only, simple_masker%init_mskproj'
        call self%kill_masker
        self%msk       = msk
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
        if( DEBUG )write(*,*)'simple_masker::automask3D done'
    end subroutine automask3D

    ! CALCULATORS

    !>  \brief  is for getting the adaptive circular mask
    real function calc_adamsk( self, img_msk )result( new_msk )
        class(masker), intent(inout) :: self
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
        if( DEBUG )write(*,*)'simple_masker::calc_adamsk done'
    end function calc_adamsk

    ! 2D CALCULATORS

    !>  \brief  is for enveloppe masking the input image
    subroutine update_cls( self, ref, cls )
        class(masker), intent(inout) :: self
        class(image),          intent(inout) :: ref
        integer,               intent(in)    :: cls
        type(image) :: img
        if( .not.self%masker_exists )stop 'masker object has not be instantiated'
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
        if( DEBUG )write(*,*)'simple_masker::update_cls done'
    end subroutine update_cls

    !>  \brief  is for binarizing the 2D image
    subroutine bin_cavg( self, img )
        class(masker), intent(inout) :: self
        class(image),          intent(inout) :: img
        type(image) :: img_pad, img_copy
        integer     :: ldim_pad(3)
        ! init
        ldim_pad(1:2) = self%idim(1:2)*2
        ldim_pad(3)   = 1
        call img_pad%new(ldim_pad,  self%get_smpd())
        img_copy = img
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
        if( DEBUG )write(*,*)'simple_masker::bin_cavg done'
    end subroutine bin_cavg

    ! 3D CALCULATORS

    subroutine bin_vol( self )
        use simple_math, only: nvoxfind
        class(masker), intent(inout) :: self
        integer :: nnvox
        call self%mask(self%msk, 'soft')
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
        if( DEBUG )write(*,*)'simple_masker::bin_vol done'
    end subroutine bin_vol

    !>  \brief  volume mask projector
    subroutine env_rproject(self, e, img)
        !$ use omp_lib
        !$ use omp_lib_kinds
        class(masker), intent(inout) :: self   !< projector instance
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
        if( DEBUG )write(*,*)'simple_masker::env_rproject done'
    end subroutine env_rproject

    ! MODIFIER

    !> \brief  applies the circular/envelope mask
    subroutine apply_mask2D( self, img, cls )
        class(masker), intent(inout) :: self
        class(image),          intent(inout) :: img
        integer,               intent(in)    :: cls
        integer :: ldim(3)
        ldim = img%get_ldim()
        if( .not.img%is_2d())stop 'only for images; simple_masker::apply_mask2D'
        if( any(ldim(:2)-self%idim(:2).ne.0) )stop 'Incompatible dimensions; simple_masker::apply_mask2D'
        if( .not.self%is_2d() )stop 'erroneous function call; simple_masker::apply_mask2D'
        if( cls > self%n )stop 'class index out of range; simple_masker::apply_mask2D'
        call img%mask(self%adamsks(cls), 'soft')
        if( DEBUG )write(*,*)'simple_masker::apply_mask2D done'
    end subroutine apply_mask2D

    ! GETTERS

    real function get_msk( self )
        class(masker), intent(inout) :: self
        get_msk = self%msk
    end function

    real function get_adamsk( self, i )
        class(masker), intent(inout) :: self
        integer, intent(in) :: i
        if( .not.allocated(self%adamsks) )stop 'adamask has not been calculated; simple_masker%get_adamsk'
        if( i > self%n )stop 'index out of range; simple_masker%get_adamsk'
        get_adamsk = self%adamsks(i)
    end function

    function get_imgmsk( self, i )result( img_msk )
        class(masker), intent(inout) :: self
        integer,               intent(in)    :: i
        type(image) :: img_msk
        if( i > self%n )stop 'index out of range; simple_masker%get_imgmsk'
        call img_msk%new([self%idim(1),self%idim(2),1], self%smpd_here)
        img_msk = 1.
        call img_msk%mask(self%adamsks(i), 'soft')
    end function

    ! DESTRUCTORS

    !>  \brief  is the destructor
    subroutine kill_masker( self )
        class(masker), intent(inout) :: self
        self%mskwidth = MSKWIDTH
        self%n        = 0
        self%msk      = 0.
        self%masker_exists = .false.
    end subroutine kill_masker

end module simple_masker


! module simple_masker
! implicit none

! !public  :: automask, automask2D
! private

! interface automask
!     module procedure automask_1
!     module procedure automask_2
!     module procedure automask_3
!     module procedure automask_4
! end interface

! contains
    
!     !>  \brief  is for generating a mask for solvent flattening of an image
!     subroutine automask_1(img, p, img_msk, nvox )
!         use simple_image,  only: image
!         use simple_params, only: params
!         use simple_math    ! use all in there
!         class(image), intent(in)      :: img
!         class(params), intent(in)     :: p
!         class(image), intent(inout)   :: img_msk
!         integer, intent(in), optional :: nvox
!         type(image)                   :: img_tmp
!         integer                       :: i, nnvox
!         call img_msk%copy(img)        ! make a copy of the image
!         if( img_msk%is_ft() )then     ! make sure that msk image is real 
!             call img_msk%bwd_ft
!         endif
!         call img_msk%bp(0., p%amsklp) ! low-pass filter the mask image
!         if( present(nvox) )then
!             nnvox = nvox
!         else
!             ! find nr of voxels corresponding to mw
!             if( p%dens > 0. )then
!                 nnvox = nvoxfind(p%smpd, p%mw, p%dens)
!             else
!                 nnvox = nvoxfind(p%smpd, p%mw)     
!             endif
!             nnvox = nint(1.1*real(nnvox))     ! this is to compensate for the low-pass filtering
!         endif
!         call img_msk%bin(nnvox)               ! binarize
!         if( present(nvox) )then
!             ! todo
!         else
!             call img_msk%grow_bins(p%binwidth)
!         endif
!         call img_msk%cos_edge(p%edge) ! real-space filter based softening of the edge
!     end subroutine automask_1
    
!     !>  \brief  is for generating and applying a mask for 
!     !!          solvent flattening of an image
!     subroutine automask_2(img, p, nvox )
!         use simple_image,  only: image
!         use simple_params, only: params
!         class(image), intent(inout)   :: img
!         class(params), intent(in)     :: p
!         integer, intent(in), optional :: nvox
!         logical                       :: didft
!         type(image)                   :: mask
!         didft = .false.
!         if( img%is_ft() )then
!             call img%bwd_ft
!             didft = .true.
!         endif
!         call automask_1(img, p, mask, nvox)
!         call img%mul(mask)
!         call mask%kill
!         if( didft ) call img%fwd_ft
!     end subroutine automask_2
    
!     !>  \brief  is for generating, applying, and writing to file 
!     !!           a mask for solvent flattening of an image
!     subroutine automask_3( b, p, cline, recvol, maskvol, volnam, masknam )
!         use simple_build,   only: build
!         use simple_params,  only: params
!         use simple_cmdline, only: cmdline
!         use simple_image,   only: image
!         class(build),   intent(inout) :: b
!         class(params),  intent(in)    :: p
!         class(cmdline), intent(inout) :: cline
!         class(image),   intent(inout) :: recvol, maskvol
!         character(len=*), intent(in)  :: volnam, masknam
!         if( .not. recvol%exists() )  stop 'recvol not allocated; automask_3; simple_masker'
!         if( .not. maskvol%exists() ) stop 'maskvol not allocated; automask_3; simple_masker'
!         call automask_4( b, p, cline, recvol, maskvol )
!         if( cline%defined('nvox') )then
!             call maskvol%write(masknam, del_if_exists=.true.)
!         else if( cline%defined('mw') )then
!             call maskvol%write(masknam, del_if_exists=.true.)
!         else if( cline%defined('mskfile') )then
!             call maskvol%write(p%mskfile, del_if_exists=.true.)
!         endif
!         call recvol%write(volnam, del_if_exists=.true.)
!     end subroutine automask_3
    
!     !>  \brief  is for generating & applying a mask for solvent flattening of an image
!     subroutine automask_4( b, p, cline, recvol, maskvol )
!         use simple_build,   only: build
!         use simple_params,  only: params
!         use simple_cmdline, only: cmdline
!         use simple_image,   only: image
!         class(build),   intent(inout) :: b
!         class(params),  intent(in)    :: p
!         class(cmdline), intent(inout) :: cline
!         class(image), intent(inout)   :: recvol, maskvol
!         if( .not. recvol%exists() )  stop 'recvol not allocated; automask_3; simple_masker'
!         if( .not. maskvol%exists() ) stop 'maskvol not allocated; automask_3; simple_masker'
!         if( cline%defined('nvox') )then
!             call automask_1(recvol, p, maskvol, p%nvox)
!             !call maskvol%mask(p%msk,'soft') ! for now
!             call recvol%mul(maskvol)         
!         else if( cline%defined('mw') )then
!             call automask_1(recvol, p, maskvol)
!             !call maskvol%mask(p%msk,'soft') ! for now
!             call recvol%mul(maskvol)
!             !call recvol%mask(p%msk,'soft')
!         else if( cline%defined('mskfile') )then
!             !call maskvol%mask(p%msk,'soft') ! for now
!             call recvol%mul(maskvol) 
!         endif
!     end subroutine automask_4
    
!     !>  \brief  is for automasking in 2D
!     subroutine automask2D( img, p, img_msk_out )
!         use simple_params,  only: params
!         use simple_image,   only: image
!         class(image),           intent(inout) :: img
!         class(params),          intent(in)    :: p
!         class(image), optional, intent(out)   :: img_msk_out
!         type(image) :: img_pad, img_msk, img_copy
!         integer     :: ldim(3), ldim_pad(3)
!         real        :: smpd, rslask
!         ldim = img%get_ldim()
!         smpd = img%get_smpd()
!         if( ldim(3) > 1 ) stop 'not for 3D images; simple_masker :: automask_5'
!         ldim_pad(1:2) = ldim(1:2)*2
!         ldim_pad(3)   = 1
!         call img_pad%new(ldim_pad, smpd)
!         call img_msk%new(ldim,     smpd)
!         img_copy = img
!         call img_copy%norm
!         call img_copy%mask(p%msk, 'soft')
!         call img_copy%pad(img_pad)
!         call img_pad%fwd_ft
!         call img_pad%bp(0., p%amsklp)
!         call img_pad%bwd_ft
!         call img_pad%bin('msk', p%msk)
!         call img_pad%grow_bin
!         call img_pad%cos_edge(p%edge)
!         call img_pad%clip(img_msk)
!         call img_msk%norm('sigm')
!         call img%mul(img_msk)
!         if( present(img_msk_out) ) img_msk_out = img_msk
!         call img_copy%kill
!         call img_pad%kill
!         call img_msk%kill
!     end subroutine automask2D

! end module simple_masker
