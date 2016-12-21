!>  \brief  SIMPLE polarft class
module simple_polarft
use simple_defs ! singleton
use simple_jiffys, only: alloc_err
implicit none

public :: polarft, test_polarft
private

type :: polarft
    private
    integer              :: nrots             !< number of radial vectors (angle index) (equal to number of components in each shell)
    integer              :: khp               !< high-pass frequency limit 
    integer              :: klp               !< low-pass frequency limit
    integer              :: ring2             !< radius of molecule
    integer              :: ldim(3)=0         !< logical dimensions of original cartesian image
    complex, allocatable :: pft(:,:)          !< first coord is angle (1:nrots), second is ring
    real, allocatable    :: sqsums(:)         !< memoized square sums for the corrleation calcualtions
    logical              :: ptcl=.true.       !< to indicate whether the polarft is a particle or a reference polarft
    logical              :: existence=.false. !< to indicate existence
  contains
    ! CONSTRUCTORS
    procedure :: new
    procedure :: copy
    ! GETTERS/SETTERS
    procedure :: convert2img
    procedure :: get_nrots
    procedure :: get_klp
    procedure :: get_khp
    procedure :: get_ring2
    procedure :: get_pft
    procedure :: get_rot
    procedure :: get_roind
    procedure :: get_coord
    procedure :: get_dims
    procedure :: get_ldim
    procedure :: get_fcomp
    procedure :: set_fcomp
    procedure :: set_ldim
    procedure :: vis
    ! CHECKUPS
    procedure :: exists
    procedure :: is_ptcl
    procedure, private :: same_dims
    generic :: operator(.eqdims.) => same_dims
    ! MODIFIERS
    procedure :: noise_fill
    procedure :: shift
    procedure :: cpsh2unsh
    ! CALCULATORS
    procedure :: memoize_sqsums
    procedure :: corr_slow
    procedure :: corr
    procedure :: corr_external
    procedure, private :: corr_shifted_1
    procedure, private :: corr_shifted_2
    generic :: corr_shifted => corr_shifted_1, corr_shifted_2
    procedure :: rotsrch_slow
    procedure :: rotsrch
    procedure :: rotsrch_shifted
    procedure :: gencorrs_slow
    procedure :: gencorrs
    procedure :: simulate_cshifted_gencorrs
    procedure :: gencorrs_shifted
    procedure :: spectrum
    ! ARITHMETICS
    procedure, private :: assign
    procedure, private :: assignr
    generic :: assignment(=) => assign, assignr
    procedure, private :: addition
    procedure, private :: subtraction
    procedure, private :: multiplication
    procedure, private :: division
    generic :: operator(+) => addition
    generic :: operator(-) => subtraction
    generic :: operator(*) => multiplication
    generic :: operator(/) => division
    procedure :: add
    procedure :: subtr
    procedure :: mul_1
    procedure :: mul_2
    generic :: mul => mul_1, mul_2
    procedure :: div_1
    procedure :: div_2
    generic :: div => div_1, div_2
    ! DESTRUCTOR
    procedure :: kill
end type

interface polarft
    module procedure constructor
end interface

! CLASS PARAMETERS/VARIABLES
complex, parameter   :: zero=cmplx(0.,0.)   !< just a complex zero
logical, parameter   :: debug=.false.       !< debug indicator
real,    allocatable :: polar_angtab(:)     !< table of angles (in degrees)
real,    allocatable :: polar_coords(:,:,:) !< table of polar coordinates (in Cartesian coordinates)
integer, allocatable :: rlims(:,:)          !< limits for matrix subsection extraction in correlation calculation
real,    allocatable :: argtx(:,:)          !< argument transfer constants for shifting polarfts
real,    allocatable :: argty(:,:)          !< argument transfer constants for shifting polarfts

contains
    
    ! CONSTRUCTORS
    
    !>  \brief  is a polar Fourier transform constructor
    function constructor( kfromto, ring2, ptcl ) result( self )
        integer, intent(in)           :: kfromto(2), ring2
        logical, intent(in), optional :: ptcl
        type(polarft)                 :: self
        call self%new(kfromto, ring2, ptcl)
    end function
    
    !>  \brief  is a polar Fourier transform constructor
    subroutine new( self, kfromto, ring2, ptcl )
        use simple_math, only: gen_polar_coords
        class(polarft), intent(inout) :: self
        integer, intent(in)           :: kfromto(2), ring2
        logical, intent(in), optional :: ptcl
        integer                       :: alloc_stat,r
        call self%kill
        if( kfromto(2)-kfromto(1) < 1 )then
            write(*,*) 'fromto:', kfromto(1), kfromto(2)
            stop 'resolution range to narrow; new; simple_polarft'
        endif
        if( ring2 > 0 )then
        else
            stop 'ring2 must be > 0; new; simple_polarft'
        endif
        self%ptcl=.true.
        if( present(ptcl) ) self%ptcl = ptcl
        call gen_polar_coords(kfromto, ring2, polar_coords, polar_angtab)
        self%khp     = kfromto(1)
        self%klp     = kfromto(2)
        self%ring2   = ring2
        self%nrots = size(polar_angtab)
        if( debug )then
            print *, 'khp:',     self%khp
            print *, 'klp:',     self%klp
            print *, 'ring2:',   self%ring2
            print *, 'nrots:', self%nrots
        endif
        if( self%ptcl )then
            allocate(self%pft(2*self%nrots,self%khp:self%klp), self%sqsums(self%nrots), stat=alloc_stat)
            ! for particles the first half 1:nrots is identical to the second half nrots+1:2*nrots
            ! the number of sqsums is equal to the number of rotations
        else
            allocate(self%pft(self%nrots,self%khp:self%klp), self%sqsums(2), stat=alloc_stat)
            ! for references the first half 1:nrots/2 represents the unshifted polarft whereas the 
            ! second half nrots/2+1:nrots represents the shifted version.
            ! We only need two square sums since the reference represents our fixed frame.
            ! self%sqsums(1) is the sum for the unshifted part (1:nrots/2) and self%sqsums(2)
            ! is the sum for the shifted part (nrots/2+1:nrots)
        endif
        call alloc_err("When alloctaing instance data: new; polarft", alloc_stat)
        allocate( rlims(self%nrots,2), argtx(self%nrots/2,self%khp:self%klp),&
        argty(self%nrots/2,self%khp:self%klp), stat=alloc_stat)
        call alloc_err("When alloctaing class data: new; polarft", alloc_stat)
        self%pft    = cmplx(0.,0.)
        self%sqsums = 0.
        do r=1,self%nrots
            rlims(r,1) = r
            rlims(r,2) = r+self%nrots/2-1
        end do
        argtx = 0.
        argty = 0.
        self%existence = .true.
    end subroutine
    
    !>  \brief  is a constructor that copies the input object
    subroutine copy( self, self_in ) 
        class(polarft), intent(inout) :: self
        class(polarft), intent(in)    :: self_in
        integer                       :: i, j
        call self%new( [self_in%khp,self_in%klp], self_in%ring2, self_in%ptcl )
        self%pft     = self_in%pft
        self%sqsums  = self_in%sqsums
        self%ldim    = self_in%ldim
    end subroutine
    
    ! GETTERS/SETTERS
    
    !>  \brief makes an image out of a polar FT (for visualization)
    function convert2img( self ) result( img )
        use simple_image, only: image
        use simple_math,  only: csq
        class(polarft), intent(in) :: self
        type(image)                :: img
        integer                    :: i, j, pdim(3)
        pdim = self%get_dims()
        call img%new([pdim(1),pdim(3)-pdim(2)+1,1],1.)
        do i=1,pdim(1)
            do j=pdim(2),pdim(3)
                call img%set([i,j,1],csq(self%pft(i,j)))
            end do
        end do
    end function
    
    !>  \brief is a getter
    pure function get_nrots( self ) result( nrots )
        class(polarft), intent(in) :: self
        integer :: nrots
        nrots = self%nrots
    end function
    
    !> \brief getter for the self object
    pure function get_klp( self ) result( klp_out )
      class(polarft), intent(in) :: self
      integer        :: klp_out
      klp_out = self%klp
    end function
    
    !> \brief getter for the self object
    pure function get_khp( self ) result( khp_out )
      class(polarft), intent(in) :: self
      integer :: khp_out
      khp_out = self%khp
    end function
    
    !> \brief getter for the self object
    pure function get_ring2( self ) result( ring2_out )
      class(polarft), intent(in) :: self
      integer :: ring2_out
      ring2_out = self%ring2
    end function
    
    !> \brief getter for the self object
    pure function get_pft( self ) result( pft_out )
      class(polarft), intent(in) :: self
      complex, allocatable :: pft_out(:,:)
      allocate( pft_out(self%nrots,self%khp:self%klp), source=self%pft(:self%nrots,:) )
    end function

    !>  \brief is a getter
    function get_rot( self, ind ) result( rot )
        class(polarft), intent(in) :: self
        integer, intent(in) :: ind
        real :: rot
        if( ind < 1 .or. ind > self%nrots )then
            stop 'i is out of range; get_rot; simple_polarft'
        endif
        rot = polar_angtab(ind)
    end function
    
    !>  \brief is a getter
    function get_roind( self, rot ) result( ind )
        class(polarft), intent(in) :: self
        real, intent(in)           :: rot
        integer                    :: ind, i, alloc_stat, loc(1)
        real, allocatable          :: dists(:)
        allocate( dists(size(polar_angtab)), stat=alloc_stat)
        call alloc_err('get_roind', alloc_stat)
        do i=1,size(polar_angtab)
            dists(i) = sqrt((polar_angtab(i)-rot)**2.)
        end do
        loc = minloc(dists)
        deallocate(dists)
        ind = loc(1)
    end function
    
    !>  \brief returns polar coordinate
    pure function get_coord( self, i, j ) result( xy )
        class(polarft), intent(in) :: self
        integer, intent(in)        :: i, j
        real                       :: xy(2)
        xy(1) = polar_coords(i,j,1)
        xy(2) = polar_coords(i,j,2)
    end function
    
    !>  \brief is a getter
    pure function get_dims( self ) result( dims )
        class(polarft), intent(in) :: self
        integer :: dims(3)
        if( self%ptcl )then
            dims(1) = self%nrots
        else
            dims(1) = self%nrots/2
        endif 
        dims(2) = self%khp 
        dims(3) = self%klp
    end function
    
    !>  \brief is a getter
    pure function get_ldim( self ) result( ldim )
        class(polarft), intent(in) :: self
        integer :: ldim(3)
        ldim = self%ldim
    end function
    
    !>  \brief is a setter
    pure function get_fcomp( self, rot, k ) result( comp )
        class(polarft), intent(in) :: self
        integer, intent(in)        :: rot, k
        complex                    :: comp
        comp = self%pft(rot,k)
    end function
    
    !>  \brief is a setter
    subroutine set_fcomp( self, rot, k, comp )
        class(polarft), intent(inout) :: self
        integer, intent(in)           :: rot, k
        complex, intent(in)           :: comp
        self%pft(rot,k) = comp
        if( self%ptcl ) self%pft(self%nrots+rot,k) = comp
    end subroutine
    
    !>  \brief is a setter
    subroutine set_ldim( self, ldim )
        use simple_math, only: is_even
        class(polarft), intent(inout) :: self
        integer, intent(in)           :: ldim(3)
        logical :: even_dims, test(3)
        integer :: i, j
        test = .false.
        test(1) = is_even(ldim(1))
        test(2) = is_even(ldim(2))
        test(3) = ldim(3) == 1
        even_dims = all(test)
        if( even_dims )then
            self%ldim = ldim
            ! memoize the argument transfer constants for shifting polarfts
            argtx(:,:) = (TWOPI*polar_coords(:self%nrots/2,:,1))/real(self%ldim(1))
            argty(:,:) = (TWOPI*polar_coords(:self%nrots/2,:,2))/real(self%ldim(2))
        else
            if( ldim(3) > 1 )then
                stop '3D polarfts are not yet supported; set_ldim; simple_polarft'
            else
                stop 'only even logical dimensions supported; set_ldim; simple_polarft'
            endif
        endif
    end subroutine
    
    !>  \brief  is for plotting a polar FT
    subroutine vis( self )
        use gnufor2
        class(polarft), intent(in) :: self
        call gnufor_image(real(self%pft),  palette='gray')
        call gnufor_image(aimag(self%pft), palette='gray')
    end subroutine

    ! CHECKUPS
    
    !>  \brief  checks for existence
    pure function exists( self ) result( yes ) 
        class(polarft), intent(in) :: self
        logical :: yes
        yes = self%existence
    end function
    
    !>  \brief  checks for same dimensions, overloaded as (.eqdims.)
    pure function same_dims( self1, self2 ) result( yep ) 
        class(polarft), intent(in) :: self1
        class(polarft), intent(in) :: self2
        logical :: yep, test(4)
        test = .false.
        test(1) = self1%nrots == self2%nrots
        test(2) = self1%khp == self2%khp
        test(3) = self1%klp == self2%klp
        test(4) = all(self1%ldim == self2%ldim)
        yep = all(test)
    end function
    
    !>  \brief  checks if ptcl
    pure function is_ptcl( self ) result( yes ) 
        class(polarft), intent(in) :: self
        logical :: yes
        yes = self%ptcl
    end function

    ! MODIFIERS
    
    !>  \brief  is for filling a polarft with noise
    subroutine noise_fill( self )
        use simple_rnd, only: ran3
        class(polarft), intent(inout) :: self
        integer :: i, j, pdim(3)
        pdim = self%get_dims()
        do i=1,pdim(1)
            do j=pdim(2),pdim(3)
                self%pft(i,j) = cmplx(ran3(),ran3())
            end do
        end do
        if( self%ptcl ) self%pft(self%nrots+1:,:) = self%pft(:self%nrots,:)
    end subroutine
    
    !>  \brief  is for shifting a polar Fourier transform
    subroutine shift( self, shvec )
        use simple_math, only: csq
        class(polarft), intent(inout) :: self     !< polar FT to shift
        real, intent(in)              :: shvec(2) !< origin shift vector 
        real                          :: argmat(self%nrots/2,self%khp:self%klp)
        complex                       :: shmat(self%nrots/2,self%khp:self%klp)
        integer                       :: i, j
        if( self%ptcl )then
            write(*,*) 'The polarft is not allowed to be a particle polarft,'
            write(*,*) 'because the space complexity is doubled for particles.'
            write(*,*) 'Thus, shifting the references is half of the work'
            stop 'shift; simple_polarft'
        endif
        ! generate the argument matrix from memoized components in argtransf
        argmat                        = argtx*shvec(1)+argty*shvec(2)
        ! generate the complex shift transformation matrix
        shmat                         = cmplx(cos(argmat),sin(argmat))
        ! shift
        self%pft(self%nrots/2+1:,:) = self%pft(:self%nrots/2,:)*shmat
        ! memoize the shifted square sums
        self%sqsums(2) = sum(csq(self%pft(self%nrots/2+1:,:)))
    end subroutine
    
    !>  \brief  is for copying the shifted part to the unshifted par (4 testing)
    subroutine cpsh2unsh( self )
         class(polarft), intent(inout) :: self !< polar FT to shift
         self%pft(:self%nrots/2,:) = self%pft(self%nrots/2+1:,:)
         self%sqsums(1) = self%sqsums(2)  
    end subroutine
    
    ! CALCULATORS
    
    !>  \brief  is for FAST memoization the complex square sums reqired for correlation calcualtion
    subroutine memoize_sqsums( self )
        use simple_math, only: csq
        class(polarft), intent(inout) :: self
        integer :: r
        if( self%ptcl )then
            self%sqsums = 0.
            do r=1,self%nrots
                self%sqsums(r) = sum(csq(self%pft(rlims(r,1):rlims(r,2),:)))
            end do
        else
            self%sqsums(1) = sum(csq(self%pft(:self%nrots/2,:)))
        endif
    end subroutine
    
    !>  \brief  is the original PRIME function for correlating two polar images,
    !!          parameterization is over rotations
    function corr_slow( self1, rot1, self2, rot2 ) result( r )
        use simple_math, only: calc_corr, csq
        class(polarft), intent(in) :: self1
        class(polarft), intent(in) :: self2
        integer, intent(in)        :: rot1, rot2
        real                       :: r, sumasq, sumbsq
        integer                    :: i, j, i1, i2
        if( self1.eqdims.self2 )then
            i1     = rot1
            i2     = rot2
            sumasq = 0.
            sumbsq = 0.
            r      = 0.
            do i=1,self1%nrots/2
                do j=self1%khp,self1%klp
                    r = r+real(self1%pft(i1,j))*real(self2%pft(i2,j))+aimag(self1%pft(i1,j))*aimag(self2%pft(i2,j))
                    sumasq = sumasq+csq(self1%pft(i1,j))
                    sumbsq = sumbsq+csq(self2%pft(i2,j))
                end do
                i1 = i1+1
                if( i1 > self1%nrots ) i1 = 1
                i2 = i2+1
                if( i2 > self2%nrots ) i2 = 1
            end do
            r = calc_corr(r,sumasq*sumbsq)
        else
            write(*,*) self1%ring2, self1%khp, self1%klp
            write(*,*) self2%ring2, self2%khp, self2%klp
            stop 'not equal dims; corr'
        endif
    end function
    
    !>  \brief  is for correlating two polar images of which one (self_ref) 
    !!          is fixed at rotation index 1 and one is parameterized (self)
    !!          over in-plane rotation index rot. This routine is optimized to
    !!          remove any index permutations (do it in one hit!)
    function corr( self_ref, self, r ) result( cc )
        class(polarft), intent(in)    :: self_ref
        class(polarft), intent(inout) :: self
        integer, intent(in)           :: r
        real                          :: cc
        cc = sum(real(self_ref%pft(:self%nrots/2,:)*conjg(self%pft(rlims(r,1):rlims(r,2),:))))
        cc = cc/sqrt(self_ref%sqsums(1)*self%sqsums(r))
    end function

    !>  \brief  is for testing
    function corr_external( self, pft ) result( cc )
        use simple_math, only: csq
        class(polarft), intent(in) :: self
        complex, intent(in)        :: pft(1:self%nrots/2,self%khp:self%klp)
        real                       :: cc, sqsum_ext
        if( self%ptcl ) stop 'this routine is only for correlating references; corr_external; simple_polarft'
        cc        = sum(real(self%pft(:self%nrots/2,:)*conjg(pft(:self%nrots/2,:))))
        sqsum_ext = sum(csq(pft(:self%nrots/2,:)))
        cc        = cc/sqrt(self%sqsums(1)*sqsum_ext)
    end function
    
    !>  \brief  as above, but the shifted reference is correlated
    function corr_shifted_1( self_ref, self, r ) result( cc )
        class(polarft), intent(in)    :: self_ref
        class(polarft), intent(inout) :: self
        integer, intent(in)           :: r
        real                          :: cc
        cc = sum(real(self_ref%pft(self%nrots/2+1:,:)*conjg(self%pft(rlims(r,1):rlims(r,2),:))))
        cc = cc/sqrt(self_ref%sqsums(2)*self%sqsums(r))
    end function

    !>  \brief  as above, but the shift vec is thrown in here
    function corr_shifted_2( self_ref, self, r, shvec ) result( cc )
        class(polarft), intent(inout) :: self_ref, self
        integer, intent(in)           :: r
        real, intent(in)              :: shvec(2)
        real                          :: cc
        call self_ref%shift(shvec)
        cc = self_ref%corr_shifted_1(self, r)
    end function
    
    !>  \brief  is the original suboruitne for in-plane rotational correlation search used in PRIME
    subroutine rotsrch_slow( self_ref, self, ang, corr, iang )
        class(polarft), intent(in)     :: self_ref
        class(polarft), intent(inout)  :: self
        real, intent(out)              :: ang, corr
        integer, intent(out), optional :: iang
        real                           :: corrs(self%nrots)
        integer                        :: r, loc(1)
        if( .not. (self_ref.eqdims.self) )then
            stop 'images need to have same dimensions; rotsrch; polarft'
        endif
        call self_ref%gencorrs_slow(self, corrs)
        loc  = maxloc(corrs)
        corr = corrs(loc(1))
        ang  = polar_angtab(loc(1))
        if( present(iang) ) iang = loc(1)
    end subroutine

    !>  \brief  is for FAST in-plane rotational correlation search
    !!          This rouitne is about 10 times faster than the original one
    subroutine rotsrch( self_ref, self, ang, corr, iang )
        class(polarft), intent(in)     :: self_ref
        class(polarft), intent(inout)  :: self
        real, intent(out)              :: ang, corr
        integer, intent(out), optional :: iang
        real                           :: corrs(self%nrots)   
        integer                        :: r, loc(1)
        if( .not. self%ptcl )then
            stop 'the rotating pft needs to be constructed as a particle pft; rotsrch; simple_polarft'
        endif
        if( .not. (self_ref.eqdims.self) )then
            print *, 'refdim: ', self_ref%get_dims()
            print *, 'ptcldim: ', self%get_dims()
            stop 'images need to have same dimensions; rotsrch; polarft'
        endif
        call self_ref%gencorrs(self, corrs)
        loc  = maxloc(corrs)
        corr = corrs(loc(1))
        ang  = polar_angtab(loc(1))
        if( present(iang) ) iang = loc(1)
    end subroutine
    
    !>  \brief  is for FAST in-plane rotational correlation search
    !!          This rouitne is about 10 times faster than the original one
    subroutine rotsrch_shifted( self_ref, self, ang, corr, iang )
        class(polarft), intent(in)     :: self_ref
        class(polarft), intent(inout)  :: self
        real, intent(out)              :: ang, corr
        integer, intent(out), optional :: iang
        real                           :: corrs(self%nrots)   
        integer                        :: r, loc(1)
        if( .not. self%ptcl )then
            stop 'the rotating pft needs to be constructed as a particle pft; rotsrch; simple_polarft'
        endif
        if( .not. (self_ref.eqdims.self) )then
            print *, 'refdim: ', self_ref%get_dims()
            print *, 'ptcldim: ', self%get_dims()
            stop 'images need to have same dimensions; rotsrch; polarft'
        endif
        call self_ref%gencorrs_shifted(self, corrs)
        loc  = maxloc(corrs)
        corr = corrs(loc(1))
        ang  = polar_angtab(loc(1))
        if( present(iang) ) iang = loc(1)
    end subroutine
    
    !>  \brief  is the original subroutine for generating 
    !!          correlation values between polar FTs in PRIME
    subroutine gencorrs_slow( self_ref, self, corrs )
        class(polarft), intent(in)    :: self_ref
        class(polarft), intent(inout) :: self
        real, intent(out)             :: corrs(self%nrots)
        integer                       :: r
        do r=1,self_ref%nrots
            corrs(r) = self_ref%corr_slow(1,self, r)
        end do
    end subroutine
    
    !>  \brief  is for FAST generation of correlation values between polar FTs
    !!          This routine is about 10 times faster than the original one
    !!          because the computations for the square sums are lifted out 
    !!          (by memoization) and the index permutations are removed by 
    !!          including an identical copy of the particle pft right next to
    !!          the original one. Space complexity for the ptcl pft hence increases
    !!          by a factor 2 but that doesn't matter because it is only ONE
    subroutine gencorrs( self_ref, self, corrs )
        class(polarft), intent(in)    :: self_ref
        class(polarft), intent(inout) :: self
        real, intent(out)             :: corrs(self%nrots)
        integer                       :: r
        if( .not. self%ptcl )then
            stop 'the rotating pft needs to be constructed as a particle pft; gencorrs; simple_polarft'
        endif
        if( .not. (self_ref.eqdims.self) )then
            stop 'images need to have same dimensions; gencorrs; polarft'
        endif
        do r=1,self_ref%nrots
            corrs(r) = self_ref%corr(self, r)
        end do
    end subroutine
    
    !>  \brief  is for simulateing gencorrs based on cshift
    subroutine simulate_cshifted_gencorrs( self )
        class(polarft), intent(inout) :: self
        complex :: foopft(self%nrots,self%khp:self%klp)
        real    :: corrs(self%nrots)
        integer :: r, refsz
        refsz = self%nrots/2
        do r=1,self%nrots
            foopft = cshift(foopft,1,1)
            corrs(r) = sum(real(self%pft(:refsz,:)*conjg(foopft(:refsz,:))))
        end do
    end subroutine
    
    !>  \brief  is for FAST generation of correlation values between polar FTs
    !!          This routine is about 10 times faster than the original one
    !!          because the computations for the square sums are lifted out 
    !!          (by memoization) and the index permutations are removed by 
    !!          including an identical copy of the particle pft right next to
    !!          the original one. Space complexity for the ptcl pft hence increases
    !!          by a factor 2 but that doesn't matter because it is only ONE
    subroutine gencorrs_shifted( self_ref, self, corrs )
        class(polarft), intent(in)    :: self_ref
        class(polarft), intent(inout) :: self
        real, intent(out)             :: corrs(self%nrots)
        integer                       :: r
        if( .not. self%ptcl )then
            stop 'the rotating pft needs to be constructed as a particle pft; gencorrs; simple_polarft'
        endif
        if( .not. (self_ref.eqdims.self) )then
            stop 'images need to have same dimensions; gencorrs; polarft'
        endif
        do r=1,self_ref%nrots
            corrs(r) = self_ref%corr_shifted(self, r)
        end do
    end subroutine
    
    !>  \brief generates the rotationally averaged spectrum
    function spectrum( self, which ) result( spec )
        use simple_math, only: csq, phase_angle
        class(polarft), intent(inout) :: self
        character(len=*), intent(in)  :: which
        real, allocatable             :: spec(:)
        integer                       :: i, j, pdim(3)
        integer                       :: alloc_stat
        allocate( spec(self%khp:self%klp), stat=alloc_stat )
        call alloc_err('spectrum; simple_polarft', alloc_stat)
        spec = 0.
        pdim = self%get_dims()
        do i=1,pdim(1)
            do j=pdim(2),pdim(3)
                if( which .eq. 'power' )then
                    spec(j) = spec(j)+csq(self%pft(i,j))
                else if( which .eq. 'absreal' )then
                    spec(j) = spec(j)+abs(real(self%pft(i,j)))
                else if( which .eq. 'absimag' )then 
                    spec(j) = spec(j)+abs(aimag(self%pft(i,j)))
                else if( which .eq. 'abs' )then 
                    spec(j) = spec(j)+cabs(self%pft(i,j))
                else if( which .eq. 'phase' )then
                    spec(j) = spec(j)+phase_angle(self%pft(i,j))
                else
                    stop 'Unsupported spectrum kind; simple_image; spectrum'
                endif
            end do
        end do
        spec = spec/real(self%nrots)
    end function
    
    ! ARITHMETICS
    
    !>  \brief  polymorphic assignment (=)
    subroutine assign( selfout, selfin )
        class(polarft), intent(inout) :: selfout
        class(polarft), intent(in)    :: selfin
        call selfout%copy(selfin)
    end subroutine
    
    !>  \brief  polymorphic assignment (=)
    subroutine assignr( self, realin )
        class(polarft), intent(inout) :: self
        real, intent(in)              :: realin
        self%pft = cmplx(realin,realin)
    end subroutine
    
    !> \brief  overloaded addition(+)
    function addition( self1, self2 ) result( self )
        class(polarft), intent(in) :: self1
        class(polarft), intent(in) :: self2
        type(polarft)              :: self
        call self%copy(self1)
        call self%add(self2)
    end function
    
    !> \brief  overloaded subtraction(-)
    function subtraction( self1, self2 ) result( self )
        class(polarft), intent(in) :: self1
        class(polarft), intent(in) :: self2
        type(polarft)              :: self
        call self%copy(self1)
        call self%subtr(self2)
    end function
    
    !> \brief  overloaded multiplication(*)
    function multiplication( self1, self2 ) result( self )
        class(polarft), intent(in) :: self1
        class(polarft), intent(in) :: self2
        type(polarft)              :: self
        call self%copy(self1)
        call self%mul(self2)
    end function
    
    !> \brief  overloaded division(/)
    function division( self1, self2 ) result( self )
        class(polarft), intent(in) :: self1
        class(polarft), intent(in) :: self2
        type(polarft)              :: self
        call self%copy(self1)
        call self%div(self2)
    end function
    
    !> \brief  summation, not overloaded
    subroutine add( self, self2add )
        class(polarft), intent(inout) :: self
        class(polarft), intent(in)    :: self2add
        if( self.eqdims.self2add )then
            if( self%ptcl.eqv.self2add%ptcl )then
                self%pft = self%pft+self2add%pft
            else
                stop 'polarfts with different ptcl status cannot be added; add; polarft'
            endif
        else
            stop 'cannot add polarfts of different dims; add; polarft'
        endif
    end subroutine
    
    !> \brief  subtraction, not overloaded
    subroutine subtr( self, self2subtr )
        class(polarft), intent(inout) :: self
        class(polarft), intent(in)    :: self2subtr
        if( self.eqdims.self2subtr )then
            if( self%ptcl.eqv.self2subtr%ptcl )then
                self%pft = self%pft-self2subtr%pft
            else
                stop 'polarfts with different ptcl status cannot be subtracted; add; polarft'
            endif
        else
            stop 'cannot subtract polarfts of different dims; subtr; polarft'
        endif
    end subroutine
    
    !> \brief  multiplication, not overloaded
    subroutine mul_1( self, self2mul )
        class(polarft), intent(inout) :: self
        class(polarft), intent(in)    :: self2mul
        if( self.eqdims.self2mul )then
            if( self%ptcl.eqv.self2mul%ptcl )then
                self%pft = self%pft*self2mul%pft
            else
                 stop 'polarfts with different ptcl status cannot be multiplied; add; polarft'
            endif
        else
            stop 'cannot multiply polarfts of different dims; mul_1; polarft'
        endif
    end subroutine
    
    !> \brief  multiplication with real constant, not overloaded
    subroutine mul_2( self, c )
        class(polarft), intent(inout) :: self
        real, intent(in) :: c
        self%pft = self%pft*c
    end subroutine
    
    !> \brief  division, not overloaded
    subroutine div_1( self, self2div )
        class(polarft), intent(inout) :: self
        class(polarft), intent(in)    :: self2div
        if( self.eqdims.self2div )then
            if( self%ptcl.eqv.self2div%ptcl )then
                self%pft = self%pft/self2div%pft
            else
                stop 'polarfts with different ptcl status cannot be divided; add; polarft'
            endif
        else
            stop 'cannot divide polarfts of different dims; div_1; polarft'
        endif
    end subroutine
    
    !> \brief  division with real constant, not overloaded
    subroutine div_2( self, c )
        class(polarft), intent(inout) :: self
        real, intent(in) :: c
        if( abs(c) < 1e-6 )then
            stop 'division with zero; div; polarft'
        else
            self%pft = self%pft/c
        endif
    end subroutine
    
    ! DESTRUCTORS
    
    !>  \brief  is a polar Fourier transform destructor
    subroutine kill( self )
        class(polarft), intent(inout) :: self
        if( allocated(polar_angtab) ) deallocate(polar_angtab)
        if( allocated(polar_coords) ) deallocate(polar_coords)
        if( allocated(rlims)        ) deallocate(rlims)
        if( allocated(argtx)        ) deallocate(argtx)
        if( allocated(argty)        ) deallocate(argty)
        if( allocated(self%pft)     ) deallocate(self%pft)
        if( allocated(self%sqsums)  ) deallocate(self%sqsums)
        self%existence = .false.
    end subroutine
    
    !>  \brief  is a unit test for polarft
    subroutine test_polarft
        use simple_jiffys, only: progress
        type(polarft) :: pft1, pft2
        integer, parameter :: NITER=3000
        real               :: corr, diff, corr_slow
        integer            :: i, j, nrots, r
        real, allocatable  :: corrs(:), corrs2(:)
        call pft1%new([2,11], 45)
        call pft2%new([2,11], 45)
        nrots = pft1%get_nrots()
        allocate(corrs(nrots), corrs2(nrots))
        call pft1%noise_fill
        call pft2%noise_fill
        call pft1%memoize_sqsums
        call pft2%memoize_sqsums
        ! compare the two correlation functions
        diff = 0.
        do r=1,nrots
            corr_slow = pft1%corr_slow(1, pft2, r)
            corr = pft1%corr(pft2, r)
            diff = diff+abs(corr-corr_slow)
        end do
        diff = diff/real(nrots)
        if( diff > 1e-6 )then
            stop 'the fast correlation does not conform with the old (slow) version; test_polarft'
        endif
        write(*,*) ">>> TESTING THE OLD SLOW CORRELATION GENERATOR"
        do i=1,NITER
            call progress(i,NITER)
            call pft1%gencorrs_slow(pft2, corrs )
        end do
        write(*,*) ">>> TESTING THE NEW FAST CORRELATION GENERATOR"
        do i=1,NITER
            call progress(i,NITER)
            call pft1%gencorrs(pft2, corrs)
        end do
        write(*,*) ">>> TESTING THE CSHIFTED CORRELATION GENERATOR"
        do i=1,NITER
            call progress(i,NITER)
            call pft1%simulate_cshifted_gencorrs
        end do
    end subroutine
    
end module simple_polarft
