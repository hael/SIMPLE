! 2D/3D envelope and adaptive masking
module simple_masker
use simple_defs    ! use all in there
use simple_image,  only: image
use simple_ori,    only: ori
use simple_params, only: params
implicit none

public :: masker
private
#include "simple_local_flags.inc"

integer, parameter :: WINSZ = 3       !< real-space filter half-width window

!> image child type
type, extends(image) :: masker
    private
    type(image)       :: img_dist             !< distances to center of image
    real, allocatable :: adamsks(:)           !< evaluated circular masks   
    real              :: msk           = 0.   !< maximum circular mask
    real              :: amsklp        = 0.   !< maximum circular mask
    real              :: mw            = 0.   !< moleclular weight (in kDa)
    real              :: frac_outliers = 0.   !< fraction of outlying voxels
    real              :: pix_thres     = 0.   !< binarisation threshold
    integer           :: edge          = 3    !< edge width
    integer           :: binwidth      = 1    !< additional layers to grow
    integer           :: n             = 0    !< number of classes
    integer           :: idim(3)       = 0    !< image dimension
  contains
    ! CONSTRUCTORS
    procedure          :: automask3D
    procedure          :: init2D
    procedure          :: init_envmask2D
    ! MASKING ROUTINES
    procedure          :: apply_adamask2ptcl_2D
    procedure          :: apply_adamask2ptcl_3D
    procedure          :: apply_3Denvmask2ptcl
    procedure          :: apply_2Denvmask22Dref
    ! BINARISATION ROUTINES
    procedure, private :: bin_cavg
    procedure, private :: bin_vol_kmeans
    procedure, private :: bin_vol_thres
    ! CALCULATORS
    procedure, private :: calc_adamsk
    procedure, private :: env_rproject 
    ! DESTRUCTOR
    procedure          :: kill_masker
end type masker

contains

    ! CONSTRUCTORS

    !>  \brief  is for 3D automasking. On output the parent volume is the envelope mask
    !!          The returned volume is envelope masked.
    subroutine automask3D( self, p, vol_inout )
        class(masker), intent(inout) :: self
        class(params), intent(in)    :: p
        class(image),  intent(inout) :: vol_inout
        logical :: was_ft
        if( vol_inout%is_2d() )stop 'automask3D is intended for volumes only, simple_masker::automask3D'
        call self%kill_masker
        self%msk      = p%msk
        self%amsklp   = p%amsklp
        self%mw       = p%mw
        self%binwidth = p%binwidth
        self%edge     = p%edge
        self%frac_outliers = 0.
        if( p%frac_outliers > 0. ) self%frac_outliers = p%frac_outliers
        self%pix_thres = 0.
        if( p%thres > 0. ) self%pix_thres = p%thres
        if( p%frac_outliers > 0. .and. p%thres > 0. )then
            write(*,*) 'either tresholding or k-means (frac_outliers or thres), not both!'
            stop 'simple_masker :: automask3D'
        endif
        write(*,'(A,F7.1,A)') '>>> AUTOMASK LOW-PASS:           ', self%amsklp,  ' ANGSTROMS'
        write(*,'(A,I7,A)'  ) '>>> AUTOMASK SOFT EDGE WIDTH:    ', self%edge,    ' PIXEL(S)'
        write(*,'(A,I7,A)'  ) '>>> AUTOMASK BINARY LAYERS WIDTH:', self%binwidth,' PIXEL(S)'
        write(*,'(A,F7.1,A)') '>>> AUTOMASK MOLECULAR WEIGHT:   ', self%mw,      ' kDa'
        was_ft = vol_inout%is_ft()
        if( vol_inout%is_ft() )call vol_inout%bwd_ft
        self = vol_inout
        ! binarize volume
        if( p%thres > 0. )then
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
    end subroutine automask3D

    !>  \brief  is a 2D constructor
    !!          on exit the parent image is untouched
    subroutine init2D(self, p, ncls)
        use simple_jiffys, only: alloc_err
        class(masker), intent(inout) :: self
        class(params), intent(in)    :: p
        integer,       intent(in)    :: ncls
        integer :: alloc_stat
        if( .not. self%is_2d()    ) stop 'this routine is intended for 2D images only, simple_masker::init2D'
        if( .not. self%even_dims()) stop 'even dimensions assumed; simple_masker::init2D'
        if(       self%is_ft()    ) stop 'real space only; simple_masker::init2D'
        call self%kill_masker
        self%idim      = [p%boxmatch, p%boxmatch, 1]
        self%n         = ncls
        self%amsklp    = p%amsklp
        self%msk       = p%msk
        self%edge      = p%edge
        self%binwidth  = p%binwidth
        call self%img_dist%new(self%idim, self%get_smpd())
        call self%img_dist%cendist
        allocate(self%adamsks(self%n), stat=alloc_stat )
        call alloc_err('in simple_masker::init2D 1', alloc_stat)
        self%adamsks = 0.
    end subroutine init2D

    !>  \brief  is for initialising the 3D envelope mask used to extract 2D masks
    !!          it is assumed that the mask is already set in the image part of the object
    subroutine init_envmask2D( self, msk )
        class(masker), intent(inout) :: self
        real,          intent(in)    :: msk
        integer :: ldim(3)
        self%idim = self%get_ldim()
        ldim      = self%idim
        ldim(3)   = 1
        self%msk  = msk
        call self%remove_edge
        call self%img_dist%new(ldim, self%get_smpd())
        call self%img_dist%cendist
    end subroutine init_envmask2D

    ! MASKING ROUTINES

    !> \brief  applies the adaptive circular mask, 2D case
    subroutine apply_adamask2ptcl_2D( self, img, cls )
        class(masker), intent(inout) :: self
        class(image),  intent(inout) :: img
        integer,       intent(in)    :: cls
        integer :: ldim(3)
        ldim = img%get_ldim()
        if( .not.img%is_2d() )                stop 'only for 2D images; simple_masker::apply_adamask2ptcl_2D'
        if( any(ldim(:2)-self%idim(:2).ne.0) )stop 'Incompatible dimensions; simple_masker::apply_adamask2ptcl_2D'
        if( .not.self%is_2d() )               stop 'erroneous function call; simple_masker::apply_adamask2ptcl_2D'
        if( cls > self%n )                    stop 'class index out of range; simple_masker::apply_adamask2ptcl_2D'
        call img%mask(self%adamsks(cls), 'soft')
    end subroutine apply_adamask2ptcl_2D

    !> \brief  applies the adaptive circular mask, 3D case
    subroutine apply_adamask2ptcl_3D( self, o, img )
        class(masker),     intent(inout) :: self
        class(ori),        intent(inout) :: o
        class(image),      intent(inout) :: img
        type(image) :: mask2D
        integer     :: ldim(3)
        real        :: smpd, mskrad
        ldim    = self%get_ldim()
        smpd    = self%get_smpd()
        ldim(3) = 1
        call mask2D%new(ldim, smpd)
        call self%env_rproject(o, mask2D)
        mskrad = self%calc_adamsk(mask2D)
        call img%mask(mskrad, 'soft')
    end subroutine apply_adamask2ptcl_3D

    !> \brief  projects & applies the 3D envelope mask to particle image
    subroutine apply_3Denvmask2ptcl( self, o, img, edge )
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
    end subroutine apply_3Denvmask2ptcl

    !>  \brief  is for envelope masking of the referecne in prime2D
    !!          with memoization of the adaptive circular 2D mask radius
    subroutine apply_2Denvmask22Dref( self, ref, cls )
        class(masker), intent(inout) :: self
        class(image),  intent(inout) :: ref
        integer,       intent(in)    :: cls
        type(image) :: img
        if( cls > self%n ) stop 'class index out of range'
        ! binarize image
        img = ref
        call self%bin_cavg(img)
        ! adaptive circular mask
        self%adamsks(cls) = self%calc_adamsk(img)
        ! soft edge mask
        call img%cos_edge(self%edge)
        ! apply envelope mask to reference
        call ref%mul(img)
        if( DEBUG )write(*,*)'simple_masker::update_cls done'
    end subroutine apply_2Denvmask22Dref

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

    !>  \brief  is for binarizing the 3D image using k-means
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
        nnvox = nvoxfind(self%get_smpd(), self%mw)     
        ! binarize again
        call self%bin(nnvox)
        ! binary layers
        call self%grow_bins(self%binwidth)
        if( DEBUG )write(*,*)'simple_masker::bin_vol done'
    end subroutine bin_vol_kmeans

    !>  \brief  is for binarizing the 3D image using thresholding
    subroutine bin_vol_thres( self )
        use simple_math, only: nvoxfind
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

    !>  \brief  is for getting the adaptive circular mask
    function calc_adamsk( self, img_msk )result( new_msk )
        class(masker), intent(inout) :: self
        class(image),  intent(inout) :: img_msk
        type(image) :: img_dist, tmp_img
        real        :: minmax(2), new_msk
        integer     :: ldim(3)
        tmp_img = img_msk
        ! multiply envelope mask with square distance matrix
        call tmp_img%mul(self%img_dist)
        ! determine circular mask size
        minmax  = tmp_img%minmax()
        new_msk = real(ceiling(minmax(2) + COSMSKHALFWIDTH))
        new_msk = min(new_msk, self%msk)
        if( DEBUG )write(*,*)'simple_masker::calc_adamsk done'
    end function calc_adamsk

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

    ! DESTRUCTOR

    !>  \brief  is the destructor
    subroutine kill_masker( self )
        class(masker), intent(inout) :: self
        call self%img_dist%kill
        if( allocated(self%adamsks) ) deallocate(self%adamsks)
    end subroutine kill_masker

end module simple_masker
