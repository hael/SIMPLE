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
use simple_defs
use simple_image,  only: image 
use simple_oris,   only: oris
use simple_math,   only: csq, calc_corr
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
    procedure, private :: pcorr_1
    procedure, private :: pcorr_2
    procedure, private :: pcorr_3
    procedure, private :: pcorr_4
    procedure, private :: pcorr_5
    procedure          :: extr_lines
    generic :: pcorr => pcorr_1, pcorr_2, pcorr_3, pcorr_4, pcorr_5
    ! PRIVATE STUFF
    procedure, private :: calc_comlin
    procedure, private :: extr_comlin
    ! DESTRUCTOR
    procedure :: kill
end type comlin

interface comlin
    module procedure constructor
end interface comlin

contains

    ! CONSTRUCTORS
    
    !>  \brief  is a constructor
    function constructor( a, fpls ) result( self )
        use simple_oris,  only: oris
        use simple_image, only: image
        class(oris), intent(in), target  :: a       !< orientations
        class(image), intent(in), target :: fpls(:) !< Fourier planes
        type(comlin)                     :: self    !< object
        call self%new( a, fpls ) 
    end function constructor

    !>  \brief  is a constructor
    subroutine new( self, a, fpls )
        use simple_oris,   only: oris
        use simple_image,  only: image
        use simple_jiffys, only: alloc_err
        class(comlin), intent(inout)     :: self    !< object
        class(oris), intent(in), target  :: a       !< orientations
        class(image), intent(in), target :: fpls(:) !< Fourier planes
        integer :: alloc_stat, j, ld_here(3), fromk, tok
        call self%kill
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
    end subroutine new
    
    !>  \brief  is for extracting common lines (for debugging purposes)
    subroutine extr_lines( self, pind, lp_dyn )
        use simple_filehandling, only: get_fileunit
        class(comlin), intent(inout) :: self
        integer, intent(in)          :: pind 
        real, intent(in)             :: lp_dyn
        real                         :: corrs(self%nptcls), sums1(self%nptcls)
        real                         :: sums2(self%nptcls)
        integer                      :: j, lims(2), filnum, recsz
        lims = self%fpls(pind)%get_clin_lims(lp_dyn)
        filnum = get_fileunit()
        inquire( iolength=recsz ) self%clines(1,:,1)
        open(unit=filnum, status='replace', file='comlins_new.bin', access='direct', form='unformatted', recl=recsz )
        do j=1,self%nptcls
            call self%extr_comlin( pind, j, lims, corrs(j), sums1(j), sums2(j) )
            write(unit=filnum,rec=j) self%clines(j,:,2)
        end do
        close(unit=filnum)
    end subroutine extr_lines

    !>  \brief  is for interpolating the common lines associated with one particle image 
    !!          and calculating the per-particle common line correlation 
    function pcorr_1( self, pind, lp_dyn ) result( corr )
        class(comlin), intent(inout) :: self
        integer, intent(in)          :: pind 
        real, intent(in)             :: lp_dyn
        real                         :: corrs(self%nptcls), sums1(self%nptcls)
        real                         :: sums2(self%nptcls), corr
        integer                      :: j, lims(2)
        ! determine the dynamic resolution range
        lims = self%fpls(pind)%get_clin_lims(lp_dyn)
        ! init
        self%foundline = .false.
        corrs = 0.
        sums1 = 0.
        sums2 = 0.       
        !$omp parallel do default(shared) private(j) schedule(auto)
        do j=1,self%nptcls
            call self%extr_comlin( pind, j, lims, corrs(j), sums1(j), sums2(j) )
        end do
        !$omp end parallel do
        if( count(self%foundline) > 0 ) then
            corr = calc_corr(sum(corrs),sum(sums1)*sum(sums2))
        else
            corr = -1.
        endif
    end function pcorr_1
    
    !>  \brief  is for interpolating the common lines associated with one particle image 
    !!          and calculating the per-particle common line correlation, with random 
    !!          sampling
    function pcorr_2( self, pind, lp_dyn, bootarr ) result( corr )
        class(comlin), intent(inout) :: self
        integer, intent(in)          :: pind 
        real, intent(in)             :: lp_dyn
        integer, intent(in)          :: bootarr(:)
        real                         :: corrs(self%nptcls), sums1(self%nptcls)
        real                         :: sums2(self%nptcls), corr
        integer                      :: j, lims(2)
        ! determine the dynamic resolution range
        lims = self%fpls(pind)%get_clin_lims(lp_dyn)
        ! init
        self%foundline = .false.
        corrs = 0.
        sums1 = 0.
        sums2 = 0.     
        !$omp parallel do default(shared) private(j) schedule(auto)
        do j=1,size(bootarr)
            call self%extr_comlin( pind, bootarr(j), lims, corrs(j), sums1(j), sums2(j) )
        end do
        !$omp end parallel do
        if( count(self%foundline) > 0 ) then
            corr = calc_corr(sum(corrs),sum(sums1)*sum(sums2))
        else
            corr = -1.
        endif
    end function pcorr_2
    
    !>  \brief  is for interpolating the common lines associated with one particle image 
    !!          and calculating the per-particle common line correlation, with state
    !!          formalism
    function pcorr_3( self, pind, lp_dyn, statearr, state ) result( corr )
        class(comlin), intent(inout) :: self
        integer, intent(in)          :: pind 
        real, intent(in)             :: lp_dyn
        integer, intent(in)          :: statearr(self%nptcls), state
        real                         :: corrs(self%nptcls), sums1(self%nptcls)
        real                         :: sums2(self%nptcls), corr
        integer                      :: j, lims(2)
        ! determine the dynamic resolution range
        lims = self%fpls(pind)%get_clin_lims(lp_dyn)
        ! init
        self%foundline = .false.
        corrs = 0.
        sums1 = 0.
        sums2 = 0.        
        !$omp parallel do default(shared) private(j) schedule(auto)
        do j=1,self%nptcls
            if( statearr(j) == state ) then
                call self%extr_comlin( pind, j, lims, corrs(j), sums1(j), sums2(j) )
            endif
        end do
        !$omp end parallel do
        if( count(self%foundline) > 0 ) then
            corr = calc_corr(sum(corrs),sum(sums1)*sum(sums2))
        else
            corr = -1.
        endif
    end function pcorr_3
    
    !>  \brief  is for interpolating the common lines associated with one particle image 
    !!          and calculating the per-particle common line correlation, with random
    !!          sampling and state formalism
    function pcorr_4( self, pind, lp_dyn, bootarr, statearr, state ) result( corr )
        class(comlin), intent(inout) :: self
        integer, intent(in)          :: pind 
        real, intent(in)             :: lp_dyn
        integer, intent(in)          :: bootarr(:)
        integer, intent(in)          :: statearr(self%nptcls), state
        real                         :: corrs(self%nptcls), sums1(self%nptcls)
        real                         :: sums2(self%nptcls), corr
        integer                      :: j, lims(2)
        ! determine the dynamic resolution range
        lims = self%fpls(pind)%get_clin_lims(lp_dyn)
        ! init
        self%foundline = .false.
        corrs = 0.
        sums1 = 0.
        sums2 = 0.
        !$omp parallel do default(shared) private(j) schedule(auto)       
        do j=1,size(bootarr)
            if( statearr(bootarr(j)) == state ) then
                call self%extr_comlin( pind, bootarr(j), lims, corrs(j), sums1(j), sums2(j) )
            endif
        end do
        !$omp end parallel do
        if( count(self%foundline) > 0 ) then
            corr = calc_corr(sum(corrs),sum(sums1)*sum(sums2))
        else
            corr = -1.
        endif
    end function pcorr_4
    
    !>  \brief  is for interpolating the common line between a pair of images
    !!          and calculating the common line correlation 
    function pcorr_5( self, iptcl, jptcl, lp_dyn ) result( corr )
        class(comlin), intent(inout) :: self
        integer, intent(in)          :: iptcl, jptcl
        real, intent(in)             :: lp_dyn
        real                         :: corr, sums1, sums2
        integer                      :: lims(2)
        ! determine the dynamic resolution range
        lims = self%fpls(iptcl)%get_clin_lims(lp_dyn)
        ! init
        self%foundline(jptcl) = .false.
        corr  = 0.
        sums1 = 0.
        sums2 = 0.
        call self%extr_comlin( iptcl, jptcl, lims, corr, sums1, sums2 )
        if( self%foundline(jptcl) )then
            corr = calc_corr(corr,sums1*sums2)
        else
            corr = -1.
        endif
    end function pcorr_5
    
    ! PRIVATE STUFF
    
    !>  \brief  calculates the 3D intersection between two planes 
    !!          defined by their normals and maps the 3D intersection 
    !!          to the coordinate systems of the respective planes
    subroutine calc_comlin( self, pind, j )
        use simple_math, only: projz
        class(comlin), intent(inout) :: self
        integer, intent(in)          :: pind, j
        real, dimension(3)           :: comlin, tmp1, tmpb1, norm1, norm2
        real                         :: scalprod, abscom
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
            comlin(:) = comlin(:)/abscom
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
    end subroutine calc_comlin
    
    !>  \brief  calculates common line algebra, interpolates the
    !!          complex vectors, and calculates corr precursors
    subroutine extr_comlin( self, pind, j, lims, corr, sumasq, sumbsq )
        class(comlin), intent(inout) :: self
        integer,intent(in)           :: pind,j,lims(2)
        real, intent(out)            :: corr,sumasq,sumbsq
        integer                      :: k
        real                         :: h1,k1,h2,k2,px,py,jx,jy
        call self%calc_comlin(pind, j)
        corr   = 0.
        sumasq = 0.
        sumbsq = 0.
        if( self%foundline(j) )then
            px = self%a%get(pind, 'x')
            py = self%a%get(pind, 'y')
            jx = self%a%get(j, 'x')
            jy = self%a%get(j, 'y')
            do k=lims(1),lims(2)
                h1 = real(k)*self%lines(j,1,1)
                k1 = real(k)*self%lines(j,2,1)
                h2 = real(k)*self%lines(j,1,2)
                k2 = real(k)*self%lines(j,2,2)
                self%clines(j,k,1) = self%fpls(pind)%extr_fcomp(h1,k1,px,py)
                self%clines(j,k,2) = self%fpls(j   )%extr_fcomp(h2,k2,jx,jy)
                corr = corr+real(self%clines(j,k,1))*real(self%clines(j,k,2))+&
                &aimag(self%clines(j,k,1))*aimag(self%clines(j,k,2))
                sumasq = sumasq+csq(self%clines(j,k,1)) 
                sumbsq = sumbsq+csq(self%clines(j,k,2))
            end do
        endif
    end subroutine extr_comlin
    
    ! DESTRUCTOR
    
    !>  \brief  is a destructor
    subroutine kill( self )
        class(comlin), intent(inout) :: self
        if( self%existence )then
            deallocate(self%clines)     
            self%a    => null()
            self%fpls => null()
        endif
    end subroutine kill
    
end module simple_comlin
