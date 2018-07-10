! Fourier index iterator
module simple_ftiter
    !include 'simple_lib.f08'
    use simple_error, only: simple_stop
    use simple_math, only: is_even, fdim
implicit none

public :: ftiter, test_ftiter
private

type :: ftiter
    private
    integer :: rlogi_lbounds(3)=[0,0,0]     !<  In each dimension, the lower bound of the real image's logical addresses
    integer :: rlogi_ubounds(3)=[0,0,0]     !<  In each dimension, the upper bound of the real image's logical addresses
    integer :: clogi_lbounds(3)=[0,0,0]     !<  In each dimension, the lower bound of the complex image's logical addresses
    integer :: clogi_ubounds(3)=[0,0,0]     !<  In each dimension, the upper bound of the complex image's logical addresses
    integer :: clogi_lbounds_all(3)=[0,0,0] !<  In each dimension, the lower bound of the complex image's logical addresses,
                                            !! including redundant Friedel mates in the negative frequencies of the first dimension
    integer :: clogi_ubounds_all(3)=[0,0,0] !<  In each dimension, the upper bound of the complex image's logical addresses,
                                            !! including redundant Friedel mates in the negative frequencies of the first dimension
    integer :: cphys_ubounds(3)=[0,0,0]     !<  In each dimension, the upper bound of the complex image's physical addresses
    integer :: ldim(3)=[1,1,1]              !< logical image dimensions
    integer :: lfnys(3)=0                   !< Nyqvist indices
    integer :: lhps(3)=0                    !< High-pass indices
    integer :: llps(3)=0                    !< Low-pass indices
    integer :: lims(3,2)=0                  !< Fourier index limits
    real    :: dsteps(3)=0.                 !< wavelengths of first components
    real    :: smpd=0.                      !< sampling distance (Angstroms per pixel)
  contains
    ! CONSTRUCTOR
    procedure :: new
    ! SETTER
    procedure :: set_hp
    procedure :: get_lhp
    procedure :: get_llp
    procedure :: get_find
    procedure :: get_lfny
    procedure :: get_lp
    procedure :: get_spat_freq
    procedure :: get_clin_lims
    ! LOOPING LIMITS
    procedure :: loop_lims
    ! LOGICAL<->PHYSICAL ADDRESS CONVERTERS
    procedure :: comp_addr_phys1
    procedure :: comp_addr_phys2
    generic   :: comp_addr_phys =>  comp_addr_phys1, comp_addr_phys2
    procedure :: comp_addr_logi_1
    procedure :: comp_addr_logi_2
    generic   :: comp_addr_logi => comp_addr_logi_1, comp_addr_logi_2
    procedure :: comp_addr_phys_orig
    ! TESTS
    procedure, private :: test_addr
end type ftiter

interface ftiter
    module procedure constructor
end interface

contains

    !>  \brief is a parameterized constructor
    function constructor( ldim, smpd ) result( self )
        real,    intent(in) :: smpd
        integer, intent(in) :: ldim(:)
        type(ftiter) :: self
        call self%new(ldim, smpd)
    end function constructor

    !>  \brief  is a parameterized constructor
    subroutine new( self, ldim, smpd )
        class(ftiter), intent(inout) :: self
        real,          intent(in)    :: smpd
        integer,       intent(in)    :: ldim(:)
        integer :: d
        if( size(ldim) == 2 )then
            self%ldim(1:2) = ldim
            self%ldim(3) = 1
        else
            self%ldim = ldim
        endif
        self%smpd   = smpd
        if( is_even(self%ldim(1)) )then
            self%lfnys  = fdim(self%ldim(1))-1
            self%dsteps = real(self%lfnys)*2.*self%smpd
        else
            self%lfnys  = fdim(self%ldim(1))-1
            self%dsteps = (real(self%lfnys)*2.+1.)*self%smpd
        endif
        self%llps   = self%lfnys ! default low-pass limits
        self%lhps   = 0          ! default high-pass limits
        do d=1,3
            if (is_even(self%ldim(d))) then
                self%rlogi_lbounds(d)     = -self%ldim(d)/2
                self%rlogi_ubounds(d)     =  self%ldim(d)/2-1
                self%clogi_lbounds(d)     = -self%ldim(d)/2
                self%clogi_ubounds(d)     =  self%ldim(d)/2-1
                self%clogi_lbounds_all(d) =  self%clogi_lbounds(d)+1
                self%clogi_ubounds_all(d) =  self%clogi_ubounds(d)
            else
                self%rlogi_lbounds(d)     = -self%ldim(d)/2
                self%rlogi_ubounds(d)     =  self%ldim(d)/2
                self%clogi_lbounds(d)     = -self%ldim(d)/2
                self%clogi_ubounds(d)     =  self%ldim(d)/2
                self%clogi_lbounds_all(d) =  self%clogi_lbounds(d)
                self%clogi_ubounds_all(d) =  self%clogi_ubounds(d)
            endif
            self%cphys_ubounds(d) = self%ldim(d)
        enddo
        ! The first dimension in the complex case is treated specially
        self%clogi_lbounds(1) = 0
        if(is_even(self%ldim(1)))then
            self%clogi_ubounds(1)     =  self%ldim(1)/2
            self%cphys_ubounds(1)     =  self%ldim(1)/2+1
            self%clogi_lbounds_all(1) = -self%clogi_ubounds(1)
            self%clogi_ubounds_all(1) =  self%clogi_ubounds(1)-1
        else
            self%clogi_ubounds(1)     = (self%ldim(1)-1)/2
            self%cphys_ubounds(1)     = (self%ldim(1)+1)/2
            self%clogi_lbounds_all(1) = -self%clogi_ubounds(1)
        endif
        if(self%ldim(3) == 1)then ! if the image is 2D, the 3rd dimension is special
            self%rlogi_lbounds(3) = 0
            self%rlogi_ubounds(3) = 0
            self%clogi_lbounds(3) = 0
            self%clogi_ubounds(3) = 0
        endif
    end subroutine new

    ! SETTERS/GETTERS

    !>  \brief  is a setter
    subroutine set_hp( self, hp )
        class(ftiter), intent(inout) :: self
        real,          intent(in)    :: hp
        self%lhps = int(self%dsteps/hp)
    end subroutine set_hp

    !>  \brief  is a getter
    pure function get_lhp( self, which ) result( hpl )
        class(ftiter), intent(in) :: self
        integer,       intent(in) :: which
        integer :: hpl
        hpl = self%lhps(which)
    end function get_lhp

    !>  \brief  is a getter
    pure function get_llp( self, which ) result( lpl )
        class(ftiter), intent(in) :: self
        integer,       intent(in) :: which
        integer :: lpl
        lpl = self%llps(which)
    end function get_llp

    !>  \brief  is a getter
    pure function get_find( self, which, res ) result( ind )
        class(ftiter), intent(in) :: self
        integer,       intent(in) :: which
        real,          intent(in) :: res
        integer :: ind
        ind = int(self%dsteps(which)/res)
    end function get_find

    !>  \brief  is a getter
    pure function get_lfny( self, which ) result( fnyl )
        class(ftiter), intent(in) :: self
        integer,       intent(in) :: which
        integer :: fnyl
        fnyl = self%lfnys(which)
    end function get_lfny

    !>  \brief  is a getter
    pure function get_lp( self, which, ind ) result( lp )
        class(ftiter), intent(in) :: self
        integer,       intent(in) :: which, ind
        real :: lp
        lp = self%dsteps(which)/real(ind)
    end function get_lp

    !>  \brief  is a getter
    pure function get_spat_freq( self, which, ind ) result( spat_freq )
        class(ftiter), intent(in) :: self
        integer,       intent(in) :: which, ind
        real :: spat_freq
        spat_freq = real(ind)/self%dsteps(which)
    end function get_spat_freq

    !>  \brief  is a getter
    pure function get_clin_lims( self, lp_dyn ) result( lims )
        class(ftiter), intent(in) :: self
        real,          intent(in) :: lp_dyn
        integer :: lims(2)
        lims(2) = dynfind( self%dsteps(1), lp_dyn, self%lfnys(1) )
        lims(1) = self%lhps(1)
    end function get_clin_lims

    ! LOOPING LIMITS

    !>  \brief is for determining loop limits for transforms
    function loop_lims( self, mode, lp_dyn ) result( lims )
        class(ftiter),  intent(in) :: self
        integer,        intent(in) :: mode
        real, optional, intent(in) :: lp_dyn
        integer                    :: lims(3,2)
        if( present(lp_dyn) )then
            select case(mode)
                case(1) ! limits for correlation calculation etc.
                    lims(1,2) = dynfind( self%dsteps(1), lp_dyn, self%lfnys(1) )
                    lims(2,2) = dynfind( self%dsteps(2), lp_dyn, self%lfnys(2) )
                    lims(1,1) = self%lhps(1)
                    lims(2,1) = -lims(2,2)
                    if( self%ldim(3) == 1 )then
                        lims(3,1) = self%clogi_lbounds(3)
                        lims(3,2) = self%clogi_ubounds(3)
                    else
                        lims(3,1) = self%lhps(3)
                        lims(3,2) = dynfind( self%dsteps(3), lp_dyn, self%lfnys(3) )
                    endif
                case DEFAULT
                     call simple_stop('undefined mode; loop_lims; simple_ftiter')
            end select
        else
            select case(mode)
                case(1) ! loop over physical addresses
                    lims(1,1) = 1
                    lims(1,2) = self%cphys_ubounds(1)
                    lims(2,1) = 1
                    lims(2,2) = self%ldim(2)
                    lims(3,1) = 1
                    lims(3,2) = self%ldim(3)
                case(2) ! loop over logical addresses
                ! (exluding redundant Friedel mates in
                ! the negative frequencies of the 1st dimension)
                    lims(1,1) = self%clogi_lbounds(1)
                    lims(1,2) = self%clogi_ubounds(1)
                    lims(2,1) = self%clogi_lbounds(2)
                    lims(2,2) = self%clogi_ubounds(2)
                    lims(3,1) = self%clogi_lbounds(3)
                    lims(3,2) = self%clogi_ubounds(3)
                case(3) ! loop over logical addresses
                ! (including redundant Friedel mates)
                    lims(1,1) = self%clogi_lbounds_all(1)
                    lims(1,2) = self%clogi_ubounds_all(1)
                    lims(2,1) = self%clogi_lbounds_all(2)
                    lims(2,2) = self%clogi_ubounds_all(2)
                    lims(3,1) = self%clogi_lbounds_all(3)
                    lims(3,2) = self%clogi_ubounds_all(3)
                case DEFAULT
                    call simple_stop('undefined mode; loop_lims; simple_ftiter')
            end select
        end if
    end function loop_lims

    ! LOGICAL<->PHYSICAL ADDRESS CONVERTERS

    !>  \brief  Convert logical address to physical address. Complex image.
    function comp_addr_phys_orig(self,logi) result(phys)
        class(ftiter), intent(in) :: self
        integer,       intent(in) :: logi(3) !<  Logical address
        integer :: phys(3) !<  Physical address
        integer :: i
        if (logi(1) .ge. 0) then
            phys = logi + 1
            ! The above is true except when in negative frequencies of
            ! 2nd or 3rd dimension
            do i=2,3
                if (logi(i) .lt. 0) phys(i) = logi(i) + self%ldim(i) + 1
            enddo
        else
            ! We are in the negative frequencies of the first dimensions,
            ! which are not defined by the output of FFTW's fwd FT,
            ! so we need to look for the Friedel mate in the positive frequencies
            ! of the first dimension
            phys = -logi + 1
            ! The above is true except when in negative frequencies of
            ! 2nd or 3rd dimension
            do i=2,3
                if (-logi(i) .lt. 0) phys(i) = -logi(i) + self%ldim(i) + 1
            enddo
        endif
    end function comp_addr_phys_orig

    pure function comp_addr_phys1(self,logi) result(phys)
        class(ftiter), intent(in) :: self
        integer,       intent(in) :: logi(3) !<  Logical address
        integer :: phys(3)                   !<  Physical address
        if (logi(1) .ge. 0) then
            phys(1) = logi(1) + 1
            phys(2) = logi(2) + 1 + MERGE(self%ldim(2),0, logi(2) < 0)
            phys(3) = logi(3) + 1 + MERGE(self%ldim(3),0, logi(3) < 0)
        else
            phys(1) = -logi(1) + 1
            phys(2) = -logi(2) + 1 + MERGE(self%ldim(2),0, -logi(2) < 0)
            phys(3) = -logi(3) + 1 + MERGE(self%ldim(3),0, -logi(3) < 0)
        endif
    end function comp_addr_phys1
    pure function comp_addr_phys2(self,h,k,m) result(phys)
        class(ftiter), intent(in) :: self
        integer,       intent(in) :: h,k,m !<  Logical address
        integer :: phys(3)                 !<  Physical address
        if (h .ge. 0) then
            phys(1) = h + 1
            phys(2) = k + 1 + MERGE(self%ldim(2),0, k < 0)
            phys(3) = m + 1 + MERGE(self%ldim(3),0, m < 0)
        else
            phys(1) = -h + 1
            phys(2) = -k + 1 + MERGE(self%ldim(2),0, -k < 0)
            phys(3) = -m + 1 + MERGE(self%ldim(3),0, -m < 0)
        endif
    end function comp_addr_phys2

    !> \brief Convert physical address to logical address. Complex image.
    function comp_addr_logi_1(self,phys) result(logi)
        class(ftiter), intent(in) :: self
        integer,       intent(in) :: phys(3) !<  Physical address
        integer :: logi(3)                   !<  Logical address
        integer :: i
        logi = phys - 1
        ! The above is true except when in negative frequencies of
        ! 2nd or 3rd dimension
        do i=2,3
            if (phys(i) .gt. self%clogi_ubounds(i)+1) logi(i) = phys(i) - self%ldim(i) - 1
        enddo
    end function comp_addr_logi_1

    !> \brief Convert physical address to logical address. Complex image.
    pure function comp_addr_logi_2(self,i,j,k) result(logi)
        class(ftiter), intent(in) :: self
        integer,       intent(in) :: i,j,k   !<  Physical address
        integer :: logi(3)                   !<  Logical address
        logi(1) = i - 1
        ! The above is true except when in negative frequencies of
        ! 2nd or 3rd dimension
        logi(2) = merge( j - self%ldim(2) - 1, j - 1 , j .gt. self%clogi_ubounds(2)+1)
        logi(3) = merge( k - self%ldim(3) - 1, k - 1 , k .gt. self%clogi_ubounds(3)+1)
    end function comp_addr_logi_2

    ! PRIVATE STUFF

    !> \brief  for finding the dynamic low-pass limit
    pure function dynfind( dstep, lp_dyn, tofny ) result( target_to )
        real,    intent(in) :: dstep, lp_dyn
        integer, intent(in) :: tofny
        integer :: target_to
        target_to = int(dstep/lp_dyn)
        if( target_to > tofny ) then
            target_to = tofny
        else if( target_to < 3 ) then
            target_to = 3
        endif
    end function dynfind

    ! TESTS

    subroutine test_ftiter
        type(ftiter) :: fit
        write(*,'(a)') '**info(simple_ftiter_unit_test): testing square dimensions'
        call fit%new([100,100,100],2.)
        call fit%test_addr
        call fit%new([100,100,1],2.)
        call fit%test_addr
        write(*,'(a)') '**info(simple_ftiter_unit_test): testing non-square dimensions'
        call fit%new([120,90,80],2.)
        call fit%test_addr
        call fit%new([120,90,1],2.)
        call fit%test_addr
         write(*,'(a)') 'SIMPLE_FTITER_UNIT_TEST COMPLETED SUCCESSFULLY ;-)'
    end subroutine test_ftiter

    !>  \brief  Test the addressing of the FT'ed image for consistency
    subroutine test_addr(self)
        class(ftiter), intent(in) :: self
        integer ::  i, j, k, logi(3), phys(3)
        write(*,'(a)') '**info(test_addr): testing phys->logi->phys address conversion'
        do k=1,self%ldim(3)
            do j=1,self%ldim(2) ! this could be: do j=1,self%cphys_ubounds(2)
                do i=1,self%cphys_ubounds(1)
                    logi = self%comp_addr_logi([i,j,k])
                    phys = self%comp_addr_phys(logi)
                    if (any([i,j,k] .ne. phys)) then
                        call simple_stop('failed complex phys->logi->phys address conversion test')
                    endif
                enddo
            enddo
        enddo
        write(*,'(a)') '**info(test_addr): testing phys->logi->phys address conversion (scalar)'
        do k=1,self%ldim(3)
            do j=1,self%ldim(2) ! this could be: do j=1,self%cphys_ubounds(2)
                do i=1,self%cphys_ubounds(1)
                    logi = self%comp_addr_logi(i,j,k)
                    phys = self%comp_addr_phys(logi(1),logi(2),logi(3))
                    if (any([i,j,k] .ne. phys)) then
                        call simple_stop('failed complex phys->logi->phys address conversion test')
                    endif
                enddo
            enddo
        enddo
        write(*,'(a)') '**info(test_addr): testing logi->phys->logi address conversion (no Friedel redundancy)'
        do k=self%clogi_lbounds(3),self%clogi_ubounds(3)
            do j=self%clogi_lbounds(2),self%clogi_ubounds(2)
                do i=self%clogi_lbounds(1),self%clogi_ubounds(1)
                    phys = self%comp_addr_phys([i,j,k])
                    logi = self%comp_addr_logi(phys)
                    if (any([i,j,k] .ne. logi)) then
                        call simple_stop('failed complex logi->phys->logi address conversion test')
                    endif
                enddo
            enddo
        enddo
        write(*,'(a)') '**info(test_addr): testing logi->phys->logi address conversion (no Friedel redundancy, scalar)'
        do k=self%clogi_lbounds(3),self%clogi_ubounds(3)
            do j=self%clogi_lbounds(2),self%clogi_ubounds(2)
                do i=self%clogi_lbounds(1),self%clogi_ubounds(1)
                    phys = self%comp_addr_phys(i,j,k)
                    logi = self%comp_addr_logi(phys(1),phys(2),phys(3))
                    if (any([i,j,k] .ne. logi)) then
                        call simple_stop('failed complex logi->phys->logi address conversion test')
                    endif
                enddo
            enddo
        enddo
        write(*,'(a)') '**info(test_addr): testing logi->phys->logi address conversion (with Friedel redundancy)'
        do k=self%clogi_lbounds_all(3),self%clogi_ubounds_all(3)
            do j=self%clogi_lbounds_all(2),self%clogi_ubounds_all(2)
                do i=self%clogi_lbounds_all(1),self%clogi_ubounds_all(1)
                    phys = self%comp_addr_phys([i,j,k])
                    logi = self%comp_addr_logi(phys)
                    if (any([i,j,k]    .ne. logi) .and. &
                        any([-i,-j,-k] .ne. logi)) then
                        write(*,'(a,3(i0,1x))') '          i,j,k   = ', i,j,k
                        write(*,'(a,3(i0,1x))') '     phys(i,j,k)  = ', phys
                        write(*,'(a,3(i0,1x))') 'logi(phys(i,j,k)) = ', logi
                        call simple_stop('failed complex logi->phys->logi address conversion test (with redundant voxels)')
                    endif
                enddo
            enddo
        enddo
         write(*,'(a)') '**info(test_addr): testing logi->phys->logi address conversion (with Friedel redundancy, scalar)'
        do k=self%clogi_lbounds_all(3),self%clogi_ubounds_all(3)
            do j=self%clogi_lbounds_all(2),self%clogi_ubounds_all(2)
                do i=self%clogi_lbounds_all(1),self%clogi_ubounds_all(1)
                    phys = self%comp_addr_phys(i,j,k)
                    logi = self%comp_addr_logi(phys(1),phys(2),phys(3))
                    if (any([i,j,k]    .ne. logi) .and. &
                        any([-i,-j,-k] .ne. logi)) then
                        write(*,'(a,3(i0,1x))') '          i,j,k   = ', i,j,k
                        write(*,'(a,3(i0,1x))') '     phys(i,j,k)  = ', phys
                        write(*,'(a,3(i0,1x))') 'logi(phys(i,j,k)) = ', logi
                        call simple_stop('failed complex logi->phys->logi address conversion test (with redundant voxels)')
                    endif
                enddo
            enddo
        enddo
    end subroutine test_addr

end module simple_ftiter
