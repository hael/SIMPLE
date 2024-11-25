!==Class simple_comlin
!
! simple_comlin is the central class for all common lines based alignment methods in _SIMPLE_.
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution or modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund, 2009-06-11.
! 
!==Changes are documented below
!
!* deugged and incorporated in the _SIMPLE_ library, HE 2009-06-25
!* re-implemented with new comlin data struct and truly OOD, HE June 14 2012
!
module simple_comlin
use simple_oris,   only: oris
use simple_math    ! singleton
use simple_defs    ! singleton
use simple_image
implicit none

public :: comlin
private

type comlin
    private
    integer                 :: nptcls=0          !< nr of ptcls
    integer                 :: xdim=0            !< Fourier dim
    complex, allocatable    :: clines(:,:,:)     !< the interpolated common lines
    real, allocatable       :: lines(:,:,:)      !< the algebraic common lines
    logical, allocatable    :: foundline(:)      !< to indicate found line or not
    class(oris), pointer    :: a=>null()         !< orientations pointer
    class(image), pointer   :: fpls(:)=>null()   !< Fourier planes pointer
    logical                 :: existence=.false. !< to indicate object existence
  contains
    ! CONSTRUCTOR
    procedure :: new
    ! PARTICLE COMMON LINE CORRELATORS
    
    ! PRIVATE STUFF
    procedure, private :: calc_comlin
    ! procedure, private :: extr_comlin
    ! DESTRUCTOR
    procedure :: kill
end type

interface comlin
    module procedure constructor
end interface

contains

    ! CONSTRUCTORS
    
    !>  \brief  is a constructor
    function constructor( a, fpls ) result( self )
        use simple_oris,  only: oris
        class(oris), intent(in), target  :: a       !< orientations
        class(image), intent(in), target :: fpls(:) !< Fourier planes
        type(comlin)                     :: self    !< object
        call self%new( a, fpls ) 
    end function

    !>  \brief  is a constructor
    subroutine new( self, a, fpls )
        use simple_oris, only: oris
        class(comlin), intent(inout)     :: self    !< object
        class(oris), intent(in), target  :: a       !< orientations
        class(image), intent(in), target :: fpls(:) !< Fourier planes
        integer :: alloc_stat, j, ld_here(3), fromk, tok
        do j=1,self%nptcls
            if(.not. fpls(j)%square_dims()) stop 'square dims assumed; new; simple_comlin'
            if(.not. fpls(j)%even_dims())   stop 'even dims assumed; new; simple_comlin'
        end do
        self%nptcls = a%get_noris()
        self%a      => a
        self%fpls   => fpls
        ld_here     = fpls(1)%get_ldim()
        self%xdim   = ld_here(1)/2
        fromk       = fpls(1)%get_lhp(1)
        tok         = fpls(1)%get_lfny(1)
        allocate( self%clines(self%nptcls,fromk:tok,2), self%foundline(self%nptcls),&
        self%lines(self%nptcls,2,2), stat=alloc_stat )
        call alloc_err('new; simple_comlin, 1', alloc_stat )
        self%clines    = cmplx(0.,0.)
        self%foundline = .false.
        self%lines     = 0.
        self%existence = .true.
    end subroutine

    ! PRIVATE STUFF
    
    !>  \brief  calculates the 3D intersection between two planes 
    !!          defined by their normals and maps the 3D intersection 
    !!          to the coordinate systems of the respective planes
    subroutine calc_comlin( self, pind, j )
        class(comlin), intent(inout) :: self
        integer, intent(in)          :: pind, j
        real, dimension(3)           :: comlin, tmp1, tmpb1, norm1, norm2
        real                         :: scalprod, abscom
        real                         :: euls1(3), euls2(3) !!!DEBUG
        if( pind == j )then
            ! no self common lines
            self%foundline(j) = .false.
            return
        endif
        norm1    = self%a%get_normal(pind)
        norm2    = self%a%get_normal(j)
        scalprod = dot_product(norm1, norm2)
        if( scalprod > 0.99 ) then
            ! identical planes have no common line
            self%foundline(j) = .false.
            return
        endif
        ! find intersection in 3D
        comlin(1) = norm1(2)*norm2(3)-norm1(3)*norm2(2)
        comlin(2) = norm1(3)*norm2(1)-norm1(1)*norm2(3)
        comlin(3) = norm1(1)*norm2(2)-norm1(2)*norm2(1)
        abscom = sqrt(comlin(1)**2+comlin(2)**2+comlin(3)**2)
        if( abscom >= 0.0001 ) then
            ! normalize
            comlin(1) = comlin(1)/abscom
            comlin(2) = comlin(2)/abscom
            comlin(3) = comlin(3)/abscom
        else
            ! identical planes have no common line
            ! this should never happen
            self%foundline(j) = .false.
            return
        endif
        ! comlin is the intersection in 3D, map to the
        ! respective coordinate systems
        ! first map onto the target
        tmp1 = matmul( self%a%get_mat(pind), comlin )
        call projz( tmp1, self%lines(j,:,1) )
        ! then map onto the reference:
        tmpb1 = matmul( self%a%get_mat(j), comlin )
        call projz( tmpb1, self%lines(j,:,2) )
        self%foundline(j) = .true.
    end subroutine
    
    ! !>  \brief  calculates common line algebra, interpolates the
    ! !!          complex vectors, and calculates corr precursors
    ! subroutine extr_comlin( self, pind, j, lims, corr, sumasq, sumbsq )
    !     class(comlin), intent(inout) :: self
    !     integer,intent(in)           :: pind,j,lims(2)
    !     real, intent(out)            :: corr,sumasq,sumbsq
    !     integer                      :: k
    !     real                         :: h1,k1,h2,k2,px,py,jx,jy
    !     call self%calc_comlin(pind, j)
    !     corr   = 0.
    !     sumasq = 0.
    !     sumbsq = 0.
    !     if( self%foundline(j) )then
    !         do k=lims(1),lims(2)
    !             h1 = real(k)*self%lines(j,1,1)
    !             k1 = real(k)*self%lines(j,2,1)
    !             h2 = real(k)*self%lines(j,1,2)
    !             k2 = real(k)*self%lines(j,2,2)
    !             call self%a%get(pind, x=px, y=py)
    !             call self%a%get(j,    x=jx, y=jy)
    !             self%clines(j,k,1) = self%fpls(pind)%get_fcomp([nint(h1),nint(k1),1],[nint(px),nint(py),1])
    !             self%clines(j,k,2) = self%fpls(j   )%get_fcomp([nint(h2),nint(k2),1],[nint(jx),nint(jy),1])
    !             corr = corr+real(self%clines(j,k,1))*real(self%clines(j,k,2))+&
    !             aimag(self%clines(j,k,1))*aimag(self%clines(j,k,2))
    !             sumasq = sumasq+csqsum(self%clines(j,k,1)) 
    !             sumbsq = sumbsq+csqsum(self%clines(j,k,2))
    !         end do
    !     endif
    ! end subroutine
    
    ! DESTRUCTOR
    
    !>  \brief  is a destructor
    subroutine kill( self )
        class(comlin), intent(inout) :: self
        integer :: i
        if( self%existence )then
            deallocate(self%clines)     
            self%a    => null()
            self%fpls => null()
        endif
    end subroutine
    
end module simple_comlin