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

real,    parameter :: MSKWIDTH = 10.     !< HALF SOFT MASK WIDTH (+/-)
logical, parameter :: DEBUG    = .false. !< debug flag
integer, parameter :: WINSZ    = 3       !< real-space filter half-width window

type, extends(image) :: masker
    private
    type(oris)        :: o_msk                     !< reference orientations
    real              :: smpd_here     = 0.        !< maximum circular mask
    real              :: msk           = 0.        !< maximum circular mask
    real              :: amsklp        = 0.        !< maximum circular mask
    real, allocatable :: adamsks(:)                !< evaluated circular masks
    real              :: mskwidth      = MSKWIDTH  !< +/- soft masking width
    real              :: dens          = 0.
    real              :: mw            = 0.
    real              :: frac_outliers = 0.
    real              :: pix_thres     = 0.
    integer           :: edge          = 3
    integer           :: binwidth      = 1
    integer           :: n             = 0
    integer           :: idim(3)       = 0         !< image dimension
    logical           :: masker_exists = .false.
  contains
    ! CONSTRUCTORS
    procedure          :: init2D
    procedure          :: automask3D
    procedure          :: init_envmask2D
    ! CALCULATORS
    procedure, private :: calc_adamsk
    ! 2D CALCULATORS
    procedure          :: update_cls
    procedure, private :: bin_cavg
    ! 3D CALCULATORS
    procedure, private :: bin_vol_kmeans
    procedure, private :: bin_vol_thres
    procedure, private :: env_rproject 
    ! MODIFIERS 
    procedure          :: apply_mask2D
    procedure          :: apply_envmask2D
    ! GETTERS
    procedure          :: get_adamsk
    procedure          :: get_msk
    procedure          :: get_imgmsk
    ! DESTRUCTOR
    procedure          :: kill_masker
end type masker

contains

    ! CONSTRUCTOR

    

    !>  \brief  is a 3D constructor and modifier
    !>  On output the parent volume is the envelope mask
    !>  The returned volume is envelope masked
    subroutine automask3D( self, vol_inout, msk, amsklp, mw, binwidth, edge, dens, frac_outliers, pix_thres )
        class(masker),  intent(inout) :: self
        class(image),   intent(inout) :: vol_inout
        real,           intent(in)    :: msk, amsklp, mw, dens
        integer,        intent(in)    :: binwidth, edge
        real, optional, intent(in)    :: frac_outliers, pix_thres
        logical :: was_ft
        if( vol_inout%is_2d() )stop 'automask3D is intended for volumes only, simple_masker%init_mskproj'
        call self%kill_masker
        self%msk       = msk
        self%amsklp    = amsklp
        self%mw        = mw
        self%binwidth  = binwidth
        self%edge      = edge
        self%dens      = dens
        self%frac_outliers = 0.
        if( present(frac_outliers) ) self%frac_outliers = frac_outliers
        self%pix_thres = 0.
        if( present(pix_thres)     ) self%pix_thres     = pix_thres 
        write(*,'(A,F7.1,A)') '>>> AUTOMASK LOW-PASS:           ', self%amsklp,  ' ANGSTROMS'
        write(*,'(A,I7,A)'  ) '>>> AUTOMASK SOFT EDGE WIDTH:    ', self%edge,    ' PIXEL(S)'
        write(*,'(A,I7,A)'  ) '>>> AUTOMASK BINARY LAYERS WIDTH:', self%binwidth,' PIXEL(S)'
        write(*,'(A,F7.1,A)') '>>> AUTOMASK MOLECULAR WEIGHT:   ', self%mw,      ' kDa'
        was_ft = vol_inout%is_ft()
        if( vol_inout%is_ft() )call vol_inout%bwd_ft
        self = vol_inout
        ! binarize volume
        if( present(pix_thres) )then
            call self%bin_vol_thres
        else
            call self%bin_vol_kmeans 
        endif
        ! add volume soft edge
        call self%cos_edge(self%edge)
        ! apply mask to volume
        call vol_inout%mul(self)
        ! the end
        if( was_ft )call vol_inout%fwd_ft
        if( DEBUG )write(*,*)'simple_masker::automask3D done'
    end subroutine automask3D

    !>  \brief  is a 2D constructor
    !>          on exit the parent image is untouched
    subroutine init2D(self, p, ncls)
        class(masker),      intent(inout) :: self
        class(params),              intent(in)    :: p
        integer,                    intent(in)    :: ncls
        integer              :: alloc_stat
        if( .not.self%is_2d() )  stop 'this routine is intended for 2D images only, simple_masker::mskref'
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

    !>  \brief  is for initialising the envelope mask used to extract 2D masks
    !!          it is assumed that the mask is already set in the image part of the object
    subroutine init_envmask2D( self, msk )
        class(masker), intent(inout) :: self
        real,          intent(in)    :: msk
        self%idim    = self%get_ldim()
        self%msk     = msk
        call self%remove_edge
    end subroutine init_envmask2D

    ! CALCULATORS

    !>  \brief  is for getting the adaptive circular mask
    function calc_adamsk( self, img_msk )result( new_msk )
        class(masker), intent(inout) :: self
        class(image),  intent(inout) :: img_msk
        type(image) :: img_dist, tmp_img
        real        :: minmax(2), new_msk
        tmp_img = img_msk
        call img_dist%new(self%idim, self%get_smpd())
        call img_dist%cendist
        ! multiply envelope mask with square distance matrix
        call tmp_img%mul(img_dist)
        ! determine circular mask size
        minmax  = tmp_img%minmax()
        new_msk = real( ceiling(minmax(2) + self%mskwidth) )
        new_msk = min(new_msk, self%msk)
        if( DEBUG )write(*,*)'simple_masker::calc_adamsk done'
    end function calc_adamsk

    ! 2D CALCULATORS

    !>  \brief  is for envelope masking the input image in prime2D
    subroutine update_cls( self, ref, cls )
        class(masker), intent(inout) :: self
        class(image),  intent(inout) :: ref
        integer,       intent(in)    :: cls
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
        ! apply envelope mask to reference
        call ref%mul(img)
        if( DEBUG )write(*,*)'simple_masker::update_cls done'
    end subroutine update_cls

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
        call img%grow_bins(1)
        if( DEBUG ) write(*,*)'simple_masker::bin_cavg done'
    end subroutine bin_cavg

    ! 3D CALCULATORS

    subroutine bin_vol_kmeans( self )
        use simple_math, only: nvoxfind
        class(masker), intent(inout) :: self
        integer :: nnvox
        ! normalize
        call self%norm()
        ! spherical mask first
        call self%mask(self%msk, 'soft')
        if( self%frac_outliers > 0. )then
            call self%bin_kmeans(self%frac_outliers)
        else
            call self%bin_kmeans
        endif
        ! apply smoothening real-space filter (to be able to FT)
        call self%real_space_filter( WINSZ, 'average')
        call self%bp(0., self%amsklp)
        ! find nr of voxels corresponding to mw
        if( self%dens > 0. )then
            nnvox = nvoxfind(self%get_smpd(), self%mw, self%dens)
        else
            nnvox = nvoxfind(self%get_smpd(), self%mw)     
        endif
        ! binarize again
        call self%bin(nnvox)
        ! binary layers
        call self%grow_bins(self%binwidth)
        if( DEBUG )write(*,*)'simple_masker::bin_vol done'
    end subroutine bin_vol_kmeans

    subroutine bin_vol_thres( self )
        use simple_math, only: nvoxfind
        class(masker), intent(inout) :: self
        integer :: nnvox
        call self%zero_below(self%pix_thres)
        call self%real_space_filter( WINSZ, 'average')
        call self%bp(0., self%amsklp)
        ! find nr of voxels corresponding to mw
        if( self%dens > 0. )then
            nnvox = nvoxfind(self%get_smpd(), self%mw, self%dens)
        else
            nnvox = nvoxfind(self%get_smpd(), self%mw)     
        endif
        ! binarize again
        call self%bin(nnvox)
        ! binary layers
        call self%grow_bins(self%binwidth)
        if( DEBUG )write(*,*)'simple_masker::bin_vol done'
    end subroutine bin_vol_thres

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

    ! MODIFIERS

    !> \brief  applies the circular/envelope mask
    subroutine apply_mask2D( self, img, cls )
        class(masker), intent(inout) :: self
        class(image),  intent(inout) :: img
        integer,       intent(in)    :: cls
        integer :: ldim(3)
        ldim = img%get_ldim()
        if( .not.img%is_2d())stop 'only for images; simple_masker::apply_mask2D'
        if( any(ldim(:2)-self%idim(:2).ne.0) )stop 'Incompatible dimensions; simple_masker::apply_mask2D'
        if( .not.self%is_2d() )stop 'erroneous function call; simple_masker::apply_mask2D'
        if( cls > self%n )stop 'class index out of range; simple_masker::apply_mask2D'
        call img%mask(self%adamsks(cls), 'soft')
        if( DEBUG )write(*,*)'simple_masker::apply_mask2D done'
    end subroutine apply_mask2D

    !> \brief  applies the envelope mask
    subroutine apply_envmask2D( self, o, img, edge )
        class(masker),     intent(inout) :: self
        class(ori),        intent(inout) :: o
        class(image),      intent(inout) :: img
        integer, optional, intent(in)    :: edge
        type(image) :: mask2D
        integer     :: ldim(3), eedge
        real        :: smpd
        eedge = 10
        if( present(edge) ) eedge = edge
        ldim    = self%get_ldim()
        smpd    = self%get_smpd()
        ldim(3) = 1
        call mask2D%new(ldim, smpd)
        call self%env_rproject(o, mask2D)
        ! add soft edge
        call mask2D%cos_edge(eedge)
        call img%mul(mask2D)
    end subroutine apply_envmask2D

    ! GETTERS

    real function get_msk( self )
        class(masker), intent(inout) :: self
        get_msk = self%msk
    end function get_msk

    real function get_adamsk( self, i )
        class(masker), intent(inout) :: self
        integer, intent(in) :: i
        if( .not.allocated(self%adamsks) )stop 'adamask has not been calculated; simple_masker%get_adamsk'
        if( i > self%n )stop 'index out of range; simple_masker%get_adamsk'
        get_adamsk = self%adamsks(i)
    end function get_adamsk

    function get_imgmsk( self, i )result( img_msk )
        class(masker), intent(inout) :: self
        integer,               intent(in)    :: i
        type(image) :: img_msk
        if( i > self%n )stop 'index out of range; simple_masker%get_imgmsk'
        call img_msk%new([self%idim(1),self%idim(2),1], self%smpd_here)
        img_msk = 1.
        call img_msk%mask(self%adamsks(i), 'soft')
    end function get_imgmsk

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
