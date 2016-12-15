module simple_ft_expanded
use simple_image, only: image
use simple_defs   ! singleton
use simple_jiffys ! singleton
implicit none

type :: ft_expanded
    private
    integer              :: lims(3,2)         !< physical limits for the Fourier transform
    integer              :: ldim(3)=[1,1,1]   !< logical dimension of originating image
    real                 :: shconst(3)        !< shift constant
    integer              :: hpind=2           !< default high-pass index
    real                 :: smpd=0.           !< sampling distance of originating image
    complex, allocatable :: cmat(:,:,:)       !< Fourier components
    logical              :: existence=.false. !< existence
  contains
    ! constructors
    procedure, private :: new_1
    procedure, private :: new_2
    generic :: new => new_1, new_2
    ! setter
    procedure :: set_fcomp
    ! checkers
    procedure, private :: same_dims
    generic   :: operator(.eqdims.) => same_dims
    ! calculators
    procedure :: corr
    procedure :: corr_shifted
    ! destructor
    procedure :: kill
end type

contains
    
    ! CONSTRUCTORS
    
    !>  \brief  is a constructor
    subroutine new_1( self, img, lp, hp )
        use simple_math, only: is_even
        class(ft_expanded), intent(inout) :: self
        class(image), intent(inout)       :: img
        real, intent(in)                  :: lp
        real, intent(in), optional        :: hp
        integer                           :: alloc_stat,h,k,l,i
        call self%kill
        self%ldim = img%get_ldim()
        if( self%ldim(3) > 1 ) stop 'currently only 4 2D images; simple_ft_expanded::new_1'
        self%smpd = img%get_smpd()
        self%lims = img%loop_lims(1,lp)
        if( present(hp) ) self%hpind = img%get_find(hp)
        ! set shift constant (shconst)
        do i=1,3
            if( is_even(self%ldim(i)) )then
                self%shconst(i) = PI/real(self%ldim(i)/2.)
            else
                self%shconst(i) = PI/real((self%ldim(i)-1)/2.)
            endif
        end do
        ! allocate 
        allocate( self%cmat(self%lims(1,1):self%lims(1,2),&
                            self%lims(2,1):self%lims(2,2),&
                            self%lims(3,1):self%lims(3,2)), stat=alloc_stat)
        call alloc_err("In: new_1; simple_ft_expanded", alloc_stat)
        self%cmat = cmplx(0.,0.)
        ! fill-up
        if( .not. img%is_ft() ) call img%fwd_ft
        do h=self%lims(1,1),self%lims(1,2)
            do k=self%lims(2,1),self%lims(2,2)
                do l=self%lims(3,1),self%lims(3,2)
                    self%cmat(h,k,l) = img%get_fcomp([h,k,l])
                end do
            end do
        end do
        self%existence = .true.
    end subroutine
    
    !>  \brief  is a constructor
    subroutine new_2( self, ldim, smpd, lims, hpind )
        class(ft_expanded), intent(inout) :: self
        integer, intent(in)               :: ldim(3), lims(3,2)
        real, intent(in)                  :: smpd             
        integer, intent(in), optional     :: hpind
        integer                           :: alloc_stat
        call self%kill
        self%ldim = ldim
        self%smpd = smpd
        self%lims = lims
        if( present(hpind) ) self%hpind = hpind
        allocate( self%cmat(self%lims(1,1):self%lims(1,2),&
                            self%lims(2,1):self%lims(2,2),&
                            self%lims(3,1):self%lims(3,2)), stat=alloc_stat)
        call alloc_err("In: new_2; simple_ft_expanded", alloc_stat)
        self%cmat = cmplx(0.,0.)
        self%existence = .true.
    end subroutine
    
    ! SETTERS
    
    !>  \brief  is a setter
    subroutine set_fcomp( self, h, k, l, comp )
        class(ft_expanded), intent(inout) :: self
        integer, intent(in)               :: h,k,l
        complex, intent(in)               :: comp
        self%cmat(h,k,l) = comp
    end subroutine
    
    ! CHECKERS
    
    !>  \brief  checks for same dimensions, overloaded as (.eqdims.)
    pure function same_dims( self1, self2 ) result( yep )
        class(ft_expanded), intent(in) :: self1, self2
        logical :: yep
        yep = all(self1%lims == self2%lims)
    end function
    
    ! CALCULATORS
    
    !>  \brief  is a correlation calculator
    function corr( self1, self2 ) result( r )
        use simple_math, only: csq, calc_corr
        class(ft_expanded), intent(in) :: self1, self2
        real :: r,sumasq,sumbsq
        if( self1.eqdims.self2 )then
            ! corr is real part of the complex mult btw 1 and 2*
            r      = sum(real(self1%cmat*conjg(self2%cmat)))
            ! normalisation terms
            sumasq = sum(csq(self1%cmat))
            sumbsq = sum(csq(self2%cmat))
            ! finalise the correlation coefficient
            r      = calc_corr(r,sumasq*sumbsq)
        else
            stop 'cannot correlate expanded_ft:s with different dims; ft_expanded::corr'
        endif
    end function
    
    !>  \brief  is a correlation calculator with origin shift of self2
    function corr_shifted( self1, self2, shvec ) result( r )
        use simple_math, only: csq, calc_corr
        class(ft_expanded), intent(in) :: self1, self2 !< instances
        real, intent(in)               :: shvec(3)
        complex, allocatable           :: shmat(:,:,:)
        real                           :: r,sumasq,sumbsq,arg
        integer                        :: alloc_stat,h,k,l
        if( self1.eqdims.self2 )then
            ! allocate 
            allocate( shmat(self1%lims(1,1):self1%lims(1,2),&
                            self1%lims(2,1):self1%lims(2,2),&
                            self1%lims(3,1):self1%lims(3,2)), stat=alloc_stat)
            call alloc_err("In: corr_shifted; simple_ft_expanded", alloc_stat)
            !$omp parallel do schedule(auto) default(shared) private(h,k,l)
            do h=self1%lims(1,1),self1%lims(1,2)
                do k=self1%lims(2,1),self1%lims(2,2)
                    do l=self1%lims(3,1),self1%lims(3,2)
                        arg          = real(h)*shvec(1)*self1%shconst(1)+&
                                       real(k)*shvec(2)*self1%shconst(2)+&
                                       real(l)*shvec(3)*self1%shconst(3)
                        shmat(h,k,l) = cmplx(cos(arg),sin(arg))
                    end do
                end do
            end do
            !$omp end parallel do
            
            !>>>>> GPU START
            
            ! corr is real part of the complex mult btw 1 and 2*
            r      = sum(real(self1%cmat**conjg(self2%cmat*shmat)))
            ! normalisation terms
            sumasq = sum(csq(self1%cmat))
            sumbsq = sum(csq(self2%cmat))
            
            !GPU END <<<<<
            
            ! finalise the correlation coefficient
            r      = calc_corr(r,sumasq*sumbsq)
            ! deallocate
            deallocate(shmat)
        else
            stop 'cannot correlate expanded_ft:s with different dims; ft_expanded::corr_shifted'
        endif
    end function
    
    ! DESTRUCTOR
    
    !>  \brief  is a destructor
    subroutine kill( self )
        class(ft_expanded), intent(inout) :: self
        if( self%existence )then
            deallocate(self%cmat)
            self%existence = .false.
        endif
    end subroutine
    
end module simple_ft_expanded