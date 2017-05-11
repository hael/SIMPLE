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
    real                     :: msk      = 0.           !< maximum circular mask
    real                     :: amsklp   = 0.           !< maximum circular mask
    real                     :: adamsk   = 0.           !< newly evaluated circular mask
    real                     :: mskwidth = MSKWIDTH     !< +/- soft masking width
    real                     :: dens     = 0.
    real                     :: mw       = 0.
    integer                  :: edge     = 1
    integer                  :: binwidth = 1
    integer                  :: nprojs   = NPROJS       !< number of masks & projection directions
    integer                  :: idim(3)  = 0            !< volume dimension
    character(len=STDLEN)    :: mode     = 'circ'       !< circular or enveloppe masking
    logical                  :: mskproj_exists   = .false.
  contains
    ! CONSTRUCTORS
    procedure          :: init_mskproj
    ! CALCULATORS
    procedure, private :: distsq_img 
    ! 2D CALCULATORS
    procedure, private :: bin_cavg
    procedure, private :: calc_cavgadamsk
    procedure, private :: build_cavgmsk
    ! 3D CALCULATORS
    procedure, private :: bin_vol
    procedure, private :: env_rproject 
    procedure, private :: adamask 
    procedure, private :: build_3Dmsks
    ! MODIFIER
    procedure, private :: apply_mask2D
    procedure, private :: apply_mask3D
    generic            :: apply_mask => apply_mask2D, apply_mask3D
    ! GETTER
    procedure          :: get_adamsk
    procedure          :: get_msk
    ! DESTRUCTOR
    procedure          :: kill_mskproj
end type mask_projector

contains

    ! CONSTRUCTOR

    !>  \brief  is a constructor
    subroutine init_mskproj(self, p, mode, nprojs, mskwidth)
        class(mask_projector),      intent(inout) :: self
        class(params),              intent(in)    :: p
        character(len=*), optional, intent(in)    :: mode
        integer,          optional, intent(in)    :: nprojs
        real,             optional, intent(in)    :: mskwidth
        type(image)          :: img_msk
        type(ori)            :: o
        real                 :: minmax(2), msk
        integer              :: i, alloc_stat
        call self%kill_mskproj
        if(.not.self%even_dims())stop 'even dimensions assumed; simple_mask_projector::init'
        if(self%is_ft())         stop 'real space only; simple_mask_projector::init'
        self%idim     = self%get_ldim()
        self%amsklp   = p%amsklp
        self%msk      = p%msk
        self%edge     = p%edge
        self%dens     = p%dens
        self%binwidth = p%binwidth
        self%mw       = p%mw
        if(present(mskwidth))self%mskwidth = mskwidth
        self%mskwidth = min(self%mskwidth, real(minval(self%idim(:2)/2))-self%msk)
        if(self%is_2d() )then
            ! 2D case
            self%nprojs = 1
            if(present(mode))then
                if(trim(mode).ne.'circ' .and. trim(mode).ne.'cavg')&
                &stop 'Unknown 2D masking mode; simple_mask_projector::init'
                self%mode = trim(mode)
            endif
            ! binarize image
            call self%bin_cavg(img_msk)
            select case(trim(mode))
                case('circ')
                    ! adaptive circular mask
                    call self%calc_cavgadamsk(img_msk)
                case('cavg')
                    ! soft edge mask
                    allocate(self%img_msks(self%nprojs))
                    call self%build_cavgmsk(img_msk)
                    self%img_msks(1) = img_msk
            end select
        else
            ! 3D
            ! init
            if(present(mode))then
                if(trim(mode).ne.'circ' .and. trim(mode).ne.'env')&
                &stop 'Unknown 3D masking mode; simple_mask_projector::init'
                self%mode = trim(mode)
            endif
            if(present(nprojs))self%nprojs = nprojs
            allocate(self%img_msks(self%nprojs), stat=alloc_stat)
            call alloc_err('in simple_mask_projector::init', alloc_stat)
            do i = 1, self%nprojs
                call self%img_msks(i)%new([self%idim(1), self%idim(2), 1], self%get_smpd())
                self%img_msks(i) = 1.
            enddo
            ! reference mask orientations
            call self%o_msk%new(self%nprojs)
            call self%o_msk%spiral(p%nsym, p%eullims)
            ! binarize volume
            call self%bin_vol(img_msk)
            ! produces and stashes masks
            call self%build_3Dmsks(img_msk)
        endif
        ! cleanup
        call img_msk%kill
        ! exists
        self%mskproj_exists = .true.
        if( DEBUG )write(*,*)'simple_mask_projector::init done'
    end subroutine init_mskproj

    ! 2D CALCULATORS

    !>  \brief  is for binarize the 2D image
    subroutine bin_cavg( self, img_msk )
        class(mask_projector), intent(inout) :: self
        class(image),          intent(out)   :: img_msk
        type(image) :: img_pad, img_copy
        integer     :: ldim_pad(3)
        ! init
        ldim_pad(1:2) = self%idim(1:2)*2
        ldim_pad(3)   = 1
        call img_pad%new(ldim_pad,  self%get_smpd())
        call img_msk%new(self%idim, self%get_smpd())
        img_copy = self
        ! normalize
        call img_copy%norm
        ! soft masking
        call img_copy%mask(self%msk, 'soft')
        ! pad
        call img_copy%pad(img_pad)
        ! low-pass
        call img_pad%fwd_ft
        call img_pad%bp(0., self%amsklp)
        call img_pad%bwd_ft
        ! binarize within mask
        call img_pad%bin('msk', self%msk)
        ! add one layer 
        call img_pad%grow_bin
        ! clip
        call img_pad%clip(img_msk)
        ! clean
        call img_copy%kill
        call img_pad%kill
    end subroutine bin_cavg

    !>  \brief  is for getting the 2D adaptive mask
    subroutine calc_cavgadamsk( self, img_msk )
        class(mask_projector), intent(inout) :: self
        class(image),          intent(inout) :: img_msk
        type(image) :: tmp_img
        real        :: new_msk, minmax(2)
        ! build sauared distance image
        call self%distsq_img(tmp_img)
        ! multiply with binary mask
        call img_msk%mul(tmp_img)
        ! set adamask
        minmax  = img_msk%minmax()
        new_msk = real(ceiling(sqrt(minmax(2))+self%mskwidth))
        self%adamsk = min(new_msk, self%msk)
        call tmp_img%kill
    end subroutine calc_cavgadamsk

    !>  \brief  is for automasking in 2D
    subroutine build_cavgmsk( self, img_msk)
        class(mask_projector), intent(inout) :: self
        class(image),          intent(inout) :: img_msk
        call img_msk%cos_edge(self%edge)
        call img_msk%norm('sigm')           ! necessary?
        self%img_msks(1) = img_msk
    end subroutine build_cavgmsk

    ! 3D CALCULATORS

    subroutine bin_vol( self, vol_msk )
        use simple_math
        class(mask_projector), intent(inout) :: self
        type(image),           intent(out) :: vol_msk
        integer :: nnvox
        vol_msk = self
        call vol_msk%bp(0., self%amsklp)
        ! find nr of voxels corresponding to mw
        if( self%dens > 0. )then
            nnvox = nvoxfind(self%get_smpd(), self%mw, self%dens)
        else
            nnvox = nvoxfind(self%get_smpd(), self%mw)     
        endif
        nnvox = nint(1.1*real(nnvox))   ! this is to compensate for the low-pass filter
        ! binarize
        call vol_msk%bin(nnvox)
        ! binary layers
        call vol_msk%grow_bins(self%binwidth)
        call vol_msk%norm_bin
    end subroutine bin_vol

    !>  \brief  volume mask projector
    subroutine env_rproject(self, e, vol, img, maxrad)
        !$ use omp_lib
        !$ use omp_lib_kinds
        class(mask_projector), intent(inout) :: self   !< projector instance
        class(ori),            intent(inout) :: e      !< Euler angle
        type(image),           intent(inout) :: vol    !< input binary volume
        type(image),           intent(inout) :: img    !< resulting projection image
        real,                  intent(in)    :: maxrad !< project inside this radius
        real, allocatable :: rmat(:,:,:)
        real              :: out_coos(3), mmaxrad, rad(3), thresh
        real              :: incr_k(3), rvec(3), rvec_k(3)
        integer           :: orig(3), i, j, k, sqmaxrad, vec(3)
        ! init
        thresh   = 0.9999             ! prior soft masking is discarded
        img      = 0.
        orig     = self%idim/2+1
        mmaxrad  = min(maxrad, real(minval(self%idim(1:2)))/2.-1.)
        sqmaxrad = nint(mmaxrad**2)
        out_coos = 0.
        rad    = 0.
        vec    = 0
        rmat   = vol%get_rmat()
        incr_k = matmul([0., 0., 1.], e%get_mat())
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

    !>  \brief  is for adaptive spherical masking
    subroutine adamask( self, img, mskrad )
        use simple_math, only: cosedge
        class(mask_projector), intent(inout) :: self
        class(image),          intent(inout) :: img
        real,                  intent(in)    :: mskrad
        real, allocatable :: rmat(:,:,:)
        real              :: e, cis(self%idim(1)), cjs(self%idim(2)), rmskbox
        integer           :: i, j, ir, jr, mskbox
        ! centre as origin
        forall(i=1:self%idim(1)) cis(i) = -real(self%idim(1)-1)/2. + real(i-1)
        forall(i=1:self%idim(2)) cjs(i) = -real(self%idim(2)-1)/2. + real(i-1)
        ! determines box size for cosedge to perform a +/- soft masking
        rmskbox = 2.*(mskrad + self%mskwidth)
        mskbox  = min( nint(rmskbox), minval(self%idim(1:2)) )
        ! MASKING
        rmat = img%get_rmat()
        do i=1,self%idim(1)/2
            ir = self%idim(1)+1-i
            do j=1,self%idim(2)/2
                jr = self%idim(2)+1-j
                e = cosedge(cis(i), cjs(j), mskbox, mskrad)
                if(e > 0.9999 )cycle
                rmat(i, j, 1) = rmat(i, j, 1) * e
                rmat(i, jr,1) = rmat(i, jr,1) * e
                rmat(ir,j, 1) = rmat(ir,j, 1) * e
                rmat(ir,jr,1) = rmat(ir,jr,1) * e
            enddo
        enddo
        call img%set_rmat(rmat)
        deallocate(rmat)      
        if( DEBUG )write(*,*)'simple_mask_projector::adamask done'    
    end subroutine adamask

    subroutine build_3Dmsks( self, vol_msk )
        class(mask_projector), intent(inout) :: self
        class(image),          intent(inout) :: vol_msk
        type(image) :: img_msk, img_dsq
        type(ori)   :: o
        real        :: new_msk, minmax(2)
        integer     :: i
        call img_msk%new([self%idim(1), self%idim(2), 1], self%get_smpd())
        if(trim(self%mode).eq.'circ')call self%distsq_img(img_dsq)
        ! fills mask images from discrete orientations
        do i = 1,self%nprojs
            o = self%o_msk%get_ori(i)
            call self%env_rproject(o, vol_msk, img_msk, self%msk)
            select case(trim(self%mode))
                case('circ')
                    ! multiply enveloppe mask with square distance matrix
                    call img_msk%mul(img_dsq)
                    ! determine circular mask size
                    minmax  = img_msk%minmax()
                    new_msk = real(ceiling(sqrt(minmax(2))+self%mskwidth))
                    new_msk = min(new_msk, self%msk)
                    ! produce & stash corresponding circular image mask
                    call self%adamask(self%img_msks(i), new_msk)
                    !call self%img_msks(i)%mask(msk,'soft')
                case('env')
                    ! add soft edge
                    call vol_msk%cos_edge(self%edge)
                    ! stash
                    self%img_msks(i) = img_msk
            end select
        enddo
        ! cleanup
        call img_dsq%kill
        call img_msk%kill
    end subroutine build_3Dmsks

    !>  \brief  produces an image with square distance from the centre of the image
    subroutine distsq_img( self, img )
        ! TO GO IN IMAGE CLASS
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
    subroutine apply_mask2D( self, img )
        class(mask_projector), intent(inout) :: self
        class(image),          intent(inout) :: img
        integer :: ldim(3)
        ldim = img%get_ldim()
        if( .not.img%is_2d())stop 'only for images; simple_mask_projector::apply_mask2D'
        if( any(ldim(:2)-self%idim(:2).ne.0) )stop 'Incompatible dimensions; simple_mask_projector::apply_mask2D'
        if( .not.self%is_2d() )stop 'erroneous function call; simple_mask_projector::apply_mask2D'
        select case(trim(self%mode))
            case('circ')
                ! adaptive circular mask
                call img%mask(self%adamsk, 'soft')
            case('cavg')
                ! soft edge mask
                call img%mul(self%img_msks(1))
        end select
        if( DEBUG )write(*,*)'simple_mask_projector::apply_mask2D done'
    end subroutine apply_mask2D

    subroutine apply_mask3D( self, img, o )
        class(mask_projector), intent(inout) :: self
        class(image),          intent(inout) :: img
        class(ori),            intent(inout) :: o
        integer :: ind, ldim(3)
        ldim = img%get_ldim()
        if( .not.img%is_2d() )stop 'only for images; simple_mask_projector::apply_mask3D'
        if( any(ldim(:2)-self%idim(:2).ne.0) )stop 'Incompatible dimensions; simple_mask_projector::apply_mask3D'
        if( self%is_2d() )stop 'erroneous function call; simple_mask_projector::apply_mask3D'
        ! index to mask with closest projection direction
        ind = self%o_msk%find_closest_proj( o )
        ! apply
        call img%mul(self%img_msks(ind))
        if( DEBUG )write(*,*)'simple_mask_projector::apply_mask3D done'
    end subroutine apply_mask3D

    ! GETTERS

    real function get_msk( self )
        class(mask_projector), intent(inout) :: self
        get_msk = self%msk
    end function

    real function get_adamsk( self )
        class(mask_projector), intent(inout) :: self
        get_adamsk = self%adamsk
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
            deallocate(self%img_msks)
        endif
        self%mode     = 'circ'
        self%mskwidth = MSKWIDTH
        self%nprojs   = NPROJS
        self%msk      = 0.
        self%adamsk   = 0.
        self%mskproj_exists = .false.
    end subroutine kill_mskproj

end module simple_mask_projector
