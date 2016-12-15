module simple_volpft_corrcalc
use simple_image,     only: image
use simple_projector, only: projector
use simple_sym,       only: sym
use simple_ori,       only: ori
use simple_defs       ! singleton
implicit none

real, parameter :: ico_alpha = 37.38     !< face-vertex angle, leads to ~4% oversampling, from Baldwin & Penczek JSB2007
real, parameter :: delta1    = ico_alpha ! -"-
real, parameter :: delta2    = ico_alpha ! -"-
real, parameter :: delta3    = 36.       !< 2pi/5/2; icosahedral inplane symmetry

type :: volpft_corrcalc
    private
    class(image), pointer :: vol_ref=>null()    !< pointer to reference volume 
    class(image), pointer :: vol_target=>null() !< pointer to target volume
    type(projector)       :: proj               !< projector object
    type(sym)             :: ico                !< defines the icosahedral group
    integer               :: nspace     = 0     !< number of vec:s in representation
    integer               :: ldim(3)    = 0     !< logical dimensions of original cartesian volume
    integer               :: kfromto(2) = 0     !< Fourier index range
    integer               :: nk         = 0     !< number of resolution elements (Fourier components)
    real                  :: hp                 !< high-pass limit
    real                  :: lp                 !< low-pass limit
    real                  :: sqsum_ref          !< memoized square sum 4 corrcalc (ref)
    real                  :: sqsum_ref_sh       !< memoized square sum 4 corrcalc (shifted ref)
    real                  :: sqsum_target       !< memoized square sum 4 corrcalc (target)
    complex, allocatable  :: vpft_ref(:,:)      !< reference lines 4 matching
    complex, allocatable  :: vpft_ref_sh(:,:)   !< reference lines 4 matching
    complex, allocatable  :: vpft_target(:,:)   !< target lines 4 matching
    real, allocatable     :: locs_ref(:,:,:)    !< nspace x nk x 3 matrix of positions (reference)
    logical               :: existence=.false.  !< to indicate existence
  contains
    ! CONSTRUCTOR
    procedure          :: new
    ! FOURIER INDEX MAPPING
    procedure, private :: k_ind
    ! GETTERS
    procedure          :: get_ldim
    procedure          :: get_kfromto
    procedure          :: get_nspace
    procedure          :: get_ori
    ! CHECKUPS
    procedure          :: exists
    ! INTERPOLATION METHODS
    procedure          :: extract_ref
    procedure          :: shift_orig_ref
    procedure, private :: shift_ref
    procedure, private :: extract_target_1
    procedure, private :: extract_target_2
    ! CORRELATORS
    procedure, private :: corr_1
    procedure, private :: corr_2
    procedure, private :: corr_3
    generic            :: corr => corr_1, corr_2, corr_3
    ! procedure          :: corr_6dim
    ! DESTRUCTOR
    procedure          :: kill
end type

contains

    !>  \brief  is a constructor
    subroutine new( self, vol_ref, vol_target, hp, lp )
        use simple_jiffys, only: alloc_err
        use simple_math,   only: is_even
        class(volpft_corrcalc), intent(inout) :: self
        class(image), intent(in), target      :: vol_ref, vol_target
        real, intent(in)                      :: hp, lp
        integer                               :: alloc_stat, isym, k, i
        real                                  :: vec(3)
        type(ori)                             :: e
        call self%kill
        if( vol_ref.eqdims.vol_target )then
            ! all good
        else
            stop 'The volumes to be matched are not of the same dimension; simple_volpft_corrcalc :: new'
        endif
        ! set pointers
        ! we assume that the volumes have been prepared with prep4cgrid
        self%vol_ref    => vol_ref
        self%vol_target => vol_target
        self%hp         =  hp
        self%lp         =  lp
        ! make projector functionality
        self%proj = projector()
        ! make the icosahedral group
        call self%ico%new('ico')
        self%nspace = self%ico%get_nsym()
        ! set other stuff
        self%ldim       = vol_ref%get_ldim()
        self%kfromto(1) = vol_ref%get_find(hp)
        self%kfromto(2) = vol_ref%get_find(lp)
        self%nk         = self%kfromto(2)-self%kfromto(1)+1
        allocate( self%vpft_ref(self%nspace,self%nk),&
                  self%vpft_ref_sh(self%nspace,self%nk),&
                  self%vpft_target(self%nspace,self%nk),&
                  self%locs_ref(self%nspace,self%nk,3), stat=alloc_stat)
        call alloc_err("In: simple_volpft_corrcalc :: new", alloc_stat)
        ! generate sampling space
        do isym=1,self%nspace
            ! get symmetry rotation matrix
            e = self%ico%get_symori(isym)
            ! loop over resolution shells
            do k=self%kfromto(1),self%kfromto(2)
                ! calculate sampling location
                vec(1) = 0.
                vec(2) = 0.
                vec(3) = real(k)
                self%locs_ref(isym,self%k_ind(k),:) = matmul(vec,e%get_mat())
            end do
        end do
        ! extract the reference lines
        call self%extract_ref
        self%existence = .true.
    end subroutine new
    
    ! FOURIER INDEX MAPPING
    
    !>  \brief  for mapping physical Fourier index to implemented
    function k_ind( self, k_in ) result( k_out )
        class(volpft_corrcalc), intent(in) :: self
        integer, intent(in)                 :: k_in
        integer :: k_out
        k_out = k_in-self%kfromto(1)+1
    end function k_ind
    
    ! GETTERS
    
    !>  \brief  for getting the logical dimension of the original
    !!          Cartesian image
    function get_ldim( self ) result( ldim )
        class(volpft_corrcalc), intent(in) :: self
        integer :: ldim(3)
        ldim = self%ldim
    end function get_ldim
    
    !>  \brief  for getting the Fourier index range (hp/lp)
    function get_kfromto( self ) result( lim )
        class(volpft_corrcalc), intent(inout) :: self
        integer :: lim(2)
        lim = self%kfromto
    end function get_kfromto
    
    !>  \brief  for getting nspace variable
    function get_nspace( self ) result( nspace )
        class(volpft_corrcalc), intent(inout) :: self
        integer :: nspace
        nspace = self%nspace
    end function get_nspace
    
    !>  \brief  for getting one of the icosahedral group orientation
    function get_ori( self, isym ) result( o )
        class(volpft_corrcalc), intent(inout) :: self
        integer, intent(in) :: isym
        type(ori) :: o
        o = self%ico%get_symori(isym)
    end function get_ori
    
    ! CHECKUPS

    !>  \brief  checks for existence
    function exists( self ) result( yes )
        class(volpft_corrcalc), intent(in) :: self
        logical :: yes
        yes = self%existence
    end function exists
    
    ! INTERPOLATION METHODS
    
    !>  \brief  extracts the lines defined by the icosahedral group 
    !!          from the reference
    subroutine extract_ref( self, e )
        use simple_math, only: csq
        class(volpft_corrcalc), intent(inout) :: self
        type(ori), intent(in), optional       :: e
        integer :: ispace, k, kind
        real    :: mat(3,3), loc(3)
        if( present(e) )then
            mat = e%get_mat()
            !$omp parallel do schedule(auto) default(shared) private(ispace,k,kind,loc)
            do ispace=1,self%nspace
                do k=self%kfromto(1),self%kfromto(2)
                    kind = self%k_ind(k)
                    loc  = matmul(self%locs_ref(ispace,kind,:),mat)
                    self%vpft_ref(ispace,kind) =&
                    self%proj%extr_gridfcomp(self%vol_ref, loc)
                end do
            end do
            !$omp end parallel do
        else
            !$omp parallel do schedule(auto) default(shared) private(ispace,k,kind)
            do ispace=1,self%nspace
                do k=self%kfromto(1),self%kfromto(2)
                    kind = self%k_ind(k)
                    self%vpft_ref(ispace,kind) =&
                    self%proj%extr_gridfcomp(self%vol_ref, self%locs_ref(ispace,kind,:))
                end do
            end do
            !$omp end parallel do
        endif
        self%sqsum_ref = sum(csq(self%vpft_ref))
    end subroutine extract_ref
    
    !>  \brief  shifts the reference lines
    subroutine shift_ref( self, shvec )
        use simple_math, only: csq
        class(volpft_corrcalc), intent(inout) :: self
        real, intent(in) :: shvec(3)
        integer :: ispace, k, kind
        !$omp parallel do schedule(auto) default(shared) private(ispace,k,kind)
        do ispace=1,self%nspace
            do k=self%kfromto(1),self%kfromto(2)
                kind = self%k_ind(k)
                self%vpft_ref_sh(ispace,kind) =&
                self%vpft_ref(ispace,kind)*&
                self%vol_ref%oshift(self%locs_ref(ispace,kind,:),shvec)
            end do
        end do
        !$omp end parallel do
        self%sqsum_ref_sh = sum(csq(self%vpft_ref))
    end subroutine shift_ref
    
    !>  \brief  shifts the reference lines
    subroutine shift_orig_ref( self, shvec )
        use simple_math, only: csq
        class(volpft_corrcalc), intent(inout) :: self
        real, intent(in) :: shvec(3)
        integer :: ispace, k, kind
        !$omp parallel do schedule(auto) default(shared) private(ispace,k,kind)
        do ispace=1,self%nspace
            do k=self%kfromto(1),self%kfromto(2)
                kind = self%k_ind(k)
                self%vpft_ref(ispace,kind) =&
                self%vpft_ref(ispace,kind)*&
                self%vol_ref%oshift(self%locs_ref(ispace,kind,:),shvec)
            end do
        end do
        !$omp end parallel do
        self%sqsum_ref_sh = sum(csq(self%vpft_ref))
    end subroutine shift_orig_ref
    
    !>  \brief  extracts the lines defined by the icosahedral group 
    !!          from the reference
    subroutine extract_target_1( self, isym )
        use simple_math, only: csq
        class(volpft_corrcalc), intent(inout) :: self
        integer, intent(in), optional         :: isym
        type(ori) :: e
        real      :: loc(3), mat(3,3)
        integer   :: ispace, k, kind
        if( present(isym) )then
            ! get symmetry op rotation matrix (4 discrete search)
            e   = self%ico%get_symori(isym)
            mat = e%get_mat()
            !$omp parallel do schedule(auto) default(shared) private(ispace,k,kind,loc)
            do ispace=1,self%nspace
                do k=self%kfromto(1),self%kfromto(2)
                    kind = self%k_ind(k)
                    loc  = matmul(self%locs_ref(ispace,kind,:),mat)
                    self%vpft_target(ispace,kind) =&
                    self%proj%extr_gridfcomp(self%vol_target, loc)
                end do
            end do
            !$omp end parallel do
        else
            !$omp parallel do schedule(auto) default(shared) private(ispace,k,kind,loc)
            do ispace=1,self%nspace
                do k=self%kfromto(1),self%kfromto(2)
                    kind = self%k_ind(k)
                    self%vpft_target(ispace,kind) =&
                    self%proj%extr_gridfcomp(self%vol_target, self%locs_ref(ispace,kind,:))
                end do
            end do
            !$omp end parallel do
        endif
        self%sqsum_target = sum(csq(self%vpft_target))
    end subroutine extract_target_1
    
    !>  \brief  extracts the lines defined by the icosahedral group 
    !!          from the reference
    subroutine extract_target_2( self, e )
        use simple_math, only: csq
        class(volpft_corrcalc), intent(inout) :: self
        class(ori), intent(in)                :: e
        real    :: loc(3), mat(3,3)
        integer :: ispace, k, kind
        mat  = e%get_mat()
        !$omp parallel do schedule(auto) default(shared) private(ispace,k,kind,loc)
        do ispace=1,self%nspace
            do k=self%kfromto(1),self%kfromto(2)
                kind = self%k_ind(k)
                loc  = matmul(self%locs_ref(ispace,kind,:),mat)
                self%vpft_target(ispace,kind) =&
                self%proj%extr_gridfcomp(self%vol_target, loc)
            end do
        end do
        !$omp end parallel do
        self%sqsum_target = sum(csq(self%vpft_target))
    end subroutine extract_target_2
    
    !>  \brief  discrete correlator
    function corr_1( self, isym ) result( cc )
        class(volpft_corrcalc), intent(inout) :: self
        integer, intent(in)                   :: isym
        real :: cc
        call self%extract_target_1(isym)
        cc = sum(real(self%vpft_ref*conjg(self%vpft_target)))
        cc = cc/sqrt(self%sqsum_target*self%sqsum_ref)
    end function corr_1
    
    !>  \brief  continous rotational correlator
    function corr_2( self, e ) result( cc )
        class(volpft_corrcalc), intent(inout) :: self
        class(ori), intent(in)                :: e
        real :: cc
        call self%extract_target_2(e)
        cc = sum(real(self%vpft_ref*conjg(self%vpft_target)))
        cc = cc/sqrt(self%sqsum_target*self%sqsum_ref)
    end function corr_2
    
    !>  \brief  continous shift correlator
    function corr_3( self, e, shvec ) result( cc )
        class(volpft_corrcalc), intent(inout) :: self
        class(ori), intent(in)                :: e
        real, intent (in)                     :: shvec(3)
        real :: cc
        call self%extract_target_2(e)
        call self%shift_ref(shvec)
        cc = sum(real(self%vpft_ref_sh*conjg(self%vpft_target)))
        cc = cc/sqrt(self%sqsum_target*self%sqsum_ref_sh)
    end function corr_3
    
    !>  \brief  continous rotational correlator
    ! function corr_6dim( self, o ) result( cc )
    !     class(volpft_corrcalc), intent(inout) :: self
    !     class(ori),             intent(inout) :: o
    !     real :: cc, shvec(3)
    !     shvec(1) = o%get('x')
    !     shvec(2) = o%get('y')
    !     shvec(3) = o%get('z')
    !     call self%extract_target_2(o)
    !     call self%shift_ref(shvec)
    !     cc = sum(real(self%vpft_ref_sh*conjg(self%vpft_target)))
    !     cc = cc/sqrt(self%sqsum_target*self%sqsum_ref_sh)
    ! end function corr_6dim

    !>  \brief  is a destructor
    subroutine kill( self )
        class(volpft_corrcalc), intent(inout) :: self
        if( self%existence )then
            self%vol_ref    => null()
            self%vol_target => null()
            call self%ico%kill
            deallocate(self%vpft_ref,self%vpft_ref_sh,self%vpft_target,self%locs_ref)
            self%existence = .false.
        endif
    end subroutine kill

end module simple_volpft_corrcalc
