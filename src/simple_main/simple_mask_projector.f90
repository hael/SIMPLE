!>  \brief  SIMPLE projector class
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

real,    parameter :: MSKWIDTH = 10.                    !< HALF SOFT MASK WIDTH (+/-)
integer, parameter :: NPROJS   = 200
logical, parameter :: DEBUG    = .false.

type, extends(image) :: mask_projector
    private
    type(oris)               :: o_msk                   !< reference orientations
    type(image), allocatable :: img_msks(:)             !< masks
    real                     :: smpd_here = 0.           !< maximum circular mask
    real                     :: msk      = 0.           !< maximum circular mask
    real                     :: amsklp   = 0.           !< maximum circular mask
    real, allocatable        :: adamsks(:)              !< evaluated circular masks
    real                     :: mskwidth = MSKWIDTH     !< +/- soft masking width
    real                     :: dens     = 0.
    real                     :: mw       = 0.
    integer                  :: edge     = 1
    integer                  :: binwidth = 1
    integer                  :: n        = NPROJS       !< number of masks and/or projection directions
    integer                  :: idim(3)  = 0            !< volume dimension
    character(len=STDLEN)    :: mode     = 'circ'       !< circular or enveloppe masking
    logical                  :: mskproj_exists = .false.
  contains
    ! CONSTRUCTORS
    procedure          :: init3D
    procedure          :: init2D
    procedure, private :: init
    ! CALCULATORS
    procedure, private :: distsq_img 
    procedure, private :: calc_adamsk
    ! 2D CALCULATORS
    procedure          :: update_cls
    procedure, private :: bin_cavg
    ! 3D CALCULATORS
    procedure, private :: bin_vol
    procedure, private :: env_rproject 
    procedure, private :: build_3Dmsks
    ! MODIFIER
    procedure, private :: apply_mask2D
    procedure, private :: apply_mask3D
    generic            :: apply_mask => apply_mask2D, apply_mask3D
    ! GETTERS
    procedure          :: get_adamsk
    procedure          :: get_msk
    procedure          :: get_imgmsk
    ! DESTRUCTOR
    procedure          :: kill_mskproj
end type mask_projector

contains

    ! CONSTRUCTOR

    !>  \brief  is a constructor
    subroutine init3D(self, p, img_inout, mode, nprojs, mskwidth)
        class(mask_projector),      intent(inout) :: self
        class(params),              intent(in)    :: p
        class(projector),           intent(inout) :: img_inout
        character(len=*), optional, intent(in)    :: mode
        integer,          optional, intent(in)    :: nprojs
        real,             optional, intent(in)    :: mskwidth
        integer  :: i, alloc_stat
        if( self%is_2d() )stop 'this routine is intended for 3D images only, simple_mask_projector::init_mskproj'
        call self%kill_mskproj
        self      = img_inout
        self%idim = self%get_ldim()
        call self%init(p, mskwidth)
        ! 3D
        self%mode = 'circ'
        if(present(nprojs))self%n = nprojs
        ! if(present(mode))then
        !     if(trim(mode).ne.'circ' .and. trim(mode).ne.'env')&
        !     &stop 'Unknown 3D masking mode; simple_mask_projector::init'
        !     self%mode = trim(mode)
        ! endif
        select case(trim(self%mode))
            case('env')
                allocate(self%img_msks(self%n), stat=alloc_stat)
                call alloc_err('in simple_mask_projector::init3D 1', alloc_stat)
                do i = 1, self%n
                    call self%img_msks(i)%new([self%idim(1), self%idim(2), 1], self%get_smpd())
                    self%img_msks(i) = 1.
                enddo
            case('circ')
                allocate(self%adamsks(self%n), stat=alloc_stat)
                call alloc_err('in simple_mask_projector::init3D 2', alloc_stat)
                self%adamsks = 0.
        end select
        ! reference mask orientations
        call self%o_msk%new(self%n)
        call self%o_msk%spiral(p%nsym, p%eullims)
        ! binarize volume
        call self%bin_vol
        ! produces and stashes masks
        call self%build_3Dmsks
        ! add volume soft edge
        call self%cos_edge(self%edge)
        ! apply mask to volume
        call img_inout%mul(self)
        ! exists
        self%mskproj_exists = .true.
        if( DEBUG )write(*,*)'simple_mask_projector::init done'
    end subroutine init3D

    !>  \brief  is a constructor
    subroutine init2D(self, p, ncls, mode)
        class(mask_projector),      intent(inout) :: self
        class(params),              intent(in)    :: p
        integer,                    intent(in)    :: ncls
        character(len=*), optional, intent(in)    :: mode
        type(ori)            :: o
        integer              :: i, alloc_stat
        if( .not.self%is_2d() )stop 'this routine is intended for 2D images only, simple_mask_projector::mskref'
        call self%kill_mskproj
        self%idim = [p%boxmatch, p%boxmatch, 1]
        self%n    = ncls
        call self%init(p)
        if(present(mode))then
            if(trim(mode).ne.'circ' .and. trim(mode).ne.'cavg')&
            &stop 'Unknown 2D masking mode; simple_mask_projector::init'
            self%mode = trim(mode)
        endif
        select case(trim(self%mode))
            case('circ')
                allocate(self%adamsks(self%n), stat=alloc_stat )
                call alloc_err('in simple_mask_projector::init2D 1', alloc_stat)
                self%adamsks = 0.
            case('cavg')
                allocate(self%img_msks(self%n), stat=alloc_stat )
                call alloc_err('in simple_mask_projector::init2D 2', alloc_stat)
                do i = 1, self%n
                    call self%img_msks(i)%new([p%boxmatch, p%boxmatch, 1], p%smpd)
                enddo
        end select
        self%mskproj_exists = .true.
        if( DEBUG )write(*,*)'simple_mask_projector::init2D'
    end subroutine init2D

    !>  \brief  is a constructor
    subroutine init(self, p, mskwidth)
        class(mask_projector), intent(inout) :: self
        class(params),         intent(in)    :: p
        real,        optional, intent(in)    :: mskwidth
        if(.not.self%even_dims())stop 'even dimensions assumed; simple_mask_projector::init'
        if(self%is_ft())         stop 'real space only; simple_mask_projector::init'
        self%amsklp   = p%amsklp
        self%msk      = p%msk
        self%edge     = p%edge
        self%dens     = p%dens
        self%binwidth = p%binwidth
        self%mw       = p%mw
        self%smpd_here = p%smpd
        if(present(mskwidth))self%mskwidth = mskwidth
        self%mskwidth = min(self%mskwidth, real(minval(self%idim(:2)/2))-self%msk)
        if( DEBUG )write(*,*)'simple_mask_projector::init done'
    end subroutine init

    ! CALCULATORS

    !>  \brief  is for getting the adaptive circular mask
    real function calc_adamsk( self, img_msk )result( new_msk )
        class(mask_projector), intent(inout) :: self
        class(image),          intent(inout) :: img_msk
        type(image) :: img_dsq, tmp_img
        real        :: minmax(2)
        tmp_img = img_msk
        call self%distsq_img(img_dsq)
        ! multiply enveloppe mask with square distance matrix
        call tmp_img%mul(img_dsq)
        ! determine circular mask size
        minmax  = tmp_img%minmax()
        new_msk = real(ceiling(sqrt(minmax(2))+self%mskwidth))
        new_msk = min(new_msk, self%msk)
    end function calc_adamsk

    ! 2D CALCULATORS

    !>  \brief  is for binarize the 2D image
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
        if(trim(self%mode).eq.'circ')then
            ! adaptive circular mask
            self%adamsks(cls) = self%calc_adamsk(img)
        endif
        ! soft edge mask
        call img%cos_edge(self%edge)
        call img%norm('sigm')
        if(trim(self%mode).eq.'cavg')then
            ! stash for particle masking
            self%img_msks(cls) = img
        endif
        call ref%mul(img)
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
        img_copy = img
        ! normalize
        call img_copy%norm
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
    end subroutine bin_cavg

    ! 3D CALCULATORS

    subroutine bin_vol( self )
        use simple_math
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
        if( DEBUG )write(*,*)'simple_mask_projector::bin_vol done'
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
        !$omp parallel default(shared) private(j,out_coos,rad,i,k,vec,rvec,rvec_k)
        !$omp do schedule(auto)
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
        !$omp end parallel
        deallocate(rmat)
        if( DEBUG )write(*,*)'simple_mask_projector::env_rproject done'
    end subroutine env_rproject

    subroutine build_3Dmsks( self )
        class(mask_projector), intent(inout) :: self
        type(image) :: img_msk
        type(ori)   :: o
        real        :: new_msk
        integer     :: i
        call img_msk%new([self%idim(1), self%idim(2), 1], self%get_smpd())
        ! fills mask images from discrete orientations
        do i = 1,self%n
            o = self%o_msk%get_ori(i)
            call self%env_rproject(o, img_msk)
            select case(trim(self%mode))
                case('circ')
                    ! evaluates new mask
                    self%adamsks(self%n) = self%calc_adamsk(img_msk)
                case('env')
                    ! add soft edge
                    call img_msk%cos_edge(self%edge)
                    ! stash
                    self%img_msks(i) = img_msk
            end select
        enddo
        ! cleanup
        call img_msk%kill
        if( DEBUG )write(*,*)'simple_mask_projector::build_3Dmsks done'
    end subroutine build_3Dmsks

    !>  \brief  produces an image with square distance from the centre of the image
    subroutine distsq_img( self, img )
        class(mask_projector), intent(inout) :: self
        class(image),          intent(inout) :: img
        real, allocatable :: rmat(:,:,:)
        real    :: centre(3)
        integer :: i, alloc_stat
        ! Builds square distance image
        call img%new([self%idim(1), self%idim(2), 1], self%get_smpd())
        allocate(rmat(1:self%idim(1),1:self%idim(2),1), stat=alloc_stat)
        call alloc_err('simple_mask_projector::distsq_img', alloc_stat)
        rmat   = 0.
        centre = real(self%idim-1)/2.
        do i=1,self%idim(1)
            rmat(i,:,1) = rmat(i,:,1) + (real(i)-centre(1))**2.
        enddo
        do i=1,self%idim(2)
            rmat(:,i,1) = rmat(:,i,1) + (real(i)-centre(2))**2.
        enddo
        call img%set_rmat(rmat)
        deallocate(rmat)
        if( DEBUG )write(*,*)'simple_mask_projector::distsq_img done'
    end subroutine distsq_img

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
        select case(trim(self%mode))
            case('circ')
                ! adaptive circular mask
                call img%mask(self%adamsks(cls), 'soft')
            case('cavg')
                ! soft edge mask
                call img%mul(self%img_msks(cls))
        end select
        if( DEBUG )write(*,*)'simple_mask_projector::apply_mask2D done'
    end subroutine apply_mask2D

    subroutine apply_mask3D( self, img, o )
        class(mask_projector), intent(inout) :: self
        class(image),          intent(inout) :: img
        class(ori),            intent(in)    :: o
        integer :: ind, ldim(3)
        ldim = img%get_ldim()
        if( .not.img%is_2d() )stop 'only for images; simple_mask_projector::apply_mask3D'
        if( any(ldim(:2)-self%idim(:2).ne.0) )stop 'Incompatible dimensions; simple_mask_projector::apply_mask3D'
        if( self%is_2d() )stop 'erroneous function call; simple_mask_projector::apply_mask3D'
        ! index to mask with closest projection direction
        ind = self%o_msk%find_closest_proj( o )
        ! apply
        select case(trim(self%mode))
            case('circ')
                ! soft addaptive circular masking
                call img%mask(self%adamsks(ind), 'soft')
            case('env')
                ! soft edge masking
                call img%mul(self%img_msks(ind))
        end select
        if( DEBUG )write(*,*)'simple_mask_projector::apply_mask3D done'
    end subroutine apply_mask3D

    ! GETTERS

    real function get_msk( self )
        class(mask_projector), intent(inout) :: self
        get_msk = self%msk
    end function

    real function get_adamsk( self, i )
        class(mask_projector), intent(inout) :: self
        integer, intent(in) :: i
        if( i > self%n )stop 'index out of range; simple_mask_projector%get_adamsk'
        get_adamsk = self%adamsks(i)
    end function

    function get_imgmsk( self, i )result( img_msk )
        class(mask_projector), intent(inout) :: self
        integer,               intent(in)    :: i
        type(image) :: img_msk
        if( i > self%n )stop 'index out of range; simple_mask_projector%get_imgmsk'
        if(trim(self%mode).eq.'circ')then
            call img_msk%new([self%idim(1),self%idim(2),1],self%smpd_here)
            img_msk = 1.
            call img_msk%mask(self%adamsks(i), 'soft')
        else
            img_msk = self%img_msks(i)
        endif
    end function

    ! DESTRUCTORS

    !>  \brief  is the destructor
    subroutine kill_mskproj( self )
        class(mask_projector), intent(inout) :: self
        integer :: i
        if(self%mskproj_exists)then
            call self%o_msk%kill
            do i = 1,size(self%img_msks)
                call self%img_msks(i)%kill
            enddo
            deallocate(self%img_msks, self%adamsks)
        endif
        self%mode     = 'circ'
        self%mskwidth = MSKWIDTH
        self%n        = NPROJS
        self%msk      = 0.
        self%mskproj_exists = .false.
    end subroutine kill_mskproj

end module simple_mask_projector
