module simple_pspecs
include 'simple_lib.f08'
!$ use omp_lib
!$ use omp_lib_kinds
use simple_image, only: image
implicit none

public :: pspecs
private
#include "simple_local_flags.inc"

type pspecs
    private
         
    real,    allocatable :: pspecs(:,:)       ! matrix of power spectra
    real,    allocatable :: apspecs(:,:)      ! class average pspecs 
    real,    allocatable :: resarr(:)         ! resolution values in A
    real,    allocatable :: dynranges(:)      ! dynamic spectral ranges
    real,    allocatable :: apows(:)          ! average powers  
    integer, allocatable :: order(:)          ! index ordering

    real                 :: hp                ! high-pass limit
    real                 :: lp                ! low-pass limit
    integer              :: box     = 0       ! box size
    integer              :: kfromto(2)        ! Fourier index range
    integer              :: sz      = 0       ! size of spectrum 
    integer              :: nspecs = 0       ! # of spectra
    integer              :: class   = 0       ! 0:     empty 
                                              ! 1-9:   good
                                              ! 10-19: bad
    logical              :: exists  = .false. ! existence flag
contains
    procedure          :: new
    procedure          :: set_pspec
    procedure          :: master
    procedure          :: calc_dynranges
    procedure          :: calc_apows
    procedure, private :: calc_distance_1, calc_distance_2
    generic            :: calc_distance => calc_distance_1, calc_distance_2
    procedure          :: kill

end type pspecs

integer, parameter :: NCLS_PSPEC = 5

contains

    subroutine new( self, nspecs, img_template, hp, lp )
        class(pspecs), intent(inout) :: self
        integer,       intent(in)    :: nspecs
        class(image),  intent(in)    :: img_template 
        real,          intent(in)    :: hp, lp
        real, allocatable :: resarr(:)
        integer :: ldim(3)
        real    :: smpd
        call self%kill
        self%nspecs    = nspecs
        self%hp         = hp
        self%lp         = lp
        ldim            = img_template%get_ldim()
        self%box        = ldim(1)
        smpd            = img_template%get_smpd()
        resarr          = get_resarr(self%box, smpd)
        self%kfromto(1) = calc_fourier_index(self%hp, self%box, smpd)
        self%kfromto(2) = calc_fourier_index(self%lp, self%box, smpd)
        self%sz         = self%kfromto(2) - self%kfromto(1) + 1
        allocate(self%resarr(self%sz), source=resarr(self%kfromto(1):self%kfromto(2)))
        allocate(self%pspecs(self%nspecs,self%sz), self%apows(self%nspecs), source=0.)
        allocate(self%order(self%sz), source=0)
        deallocate(resarr)
        self%exists     = .true.
    end subroutine new

    ! should be parallelized outside of here
    subroutine set_pspec( self, ispec, img, msk )
        class(pspecs), intent(inout) :: self
        integer,       intent(in)    :: ispec
        class(image),  intent(inout) :: img
        real,          intent(in)    :: msk
        real, allocatable :: spec(:)
        if( ispec < 1 .or. ispec > self%nspecs ) THROW_HARD('ispec index out of range')
        call img%norm
        call img%mask(msk, 'soft')
        call img%spectrum('sqrt', spec)
        self%pspecs(ispec,:) = spec(self%kfromto(1):self%kfromto(2))
        deallocate(spec)
    end subroutine set_pspec

    subroutine master( self )
        class(pspecs), intent(inout) :: self
        real    :: med_apow, mad_apow
        integer :: i
        ! assume object created
        ! assume pspecs set



        call self%calc_apows

        self%order = (/(i,i=1,self%nspecs)/)
        call hpsort(self%apows, self%order)
        do i = 1, self%nspecs
            print *, self%apows(i)
        enddo




        med_apow = median(self%apows)
        mad_apow = mad(self%apows, med_apow)

        print *, 'med_apow ', med_apow
        print *, 'mad_apow ', mad_apow

        call self%calc_dynranges

    end subroutine master

    subroutine calc_apows( self )
        class(pspecs), intent(inout) :: self
        integer :: i
        !$omp parallel do default(shared) private(i) proc_bind(close)
        do i = 1, self%nspecs
            self%apows(i) = sum(self%pspecs(i,:)) / real(self%sz)
        end do
        !$omp end parallel do
    end subroutine calc_apows

    subroutine calc_dynranges( self )
        class(pspecs), intent(inout) :: self
        if( allocated(self%dynranges) ) deallocate(self%dynranges)
        allocate(self%dynranges(self%nspecs), source=self%pspecs(:,1) - self%pspecs(:,self%sz))
    end subroutine calc_dynranges

    function calc_distance_1( self, i, j ) result( dist )
        class(pspecs), intent(in) :: self
        integer,       intent(in) :: i, j
        real :: dist
        dist = euclid(self%pspecs(i,:), self%pspecs(j,:))
    end function calc_distance_1

    function calc_distance_2( self, i, spec ) result( dist )
        class(pspecs), intent(in) :: self
        integer,       intent(in) :: i
        real,          intent(in) :: spec(self%sz)
        real :: dist
        dist = euclid(self%pspecs(i,:), spec)
    end function calc_distance_2

    subroutine kill( self )
        class(pspecs), intent(inout) :: self
        if( self%exists )then
            if( allocated(self%pspecs)    ) deallocate(self%pspecs)
            if( allocated(self%apspecs)   ) deallocate(self%apspecs)
            if( allocated(self%resarr)    ) deallocate(self%resarr)
            if( allocated(self%dynranges) ) deallocate(self%dynranges)
            if( allocated(self%apows)     ) deallocate(self%apows)
            if( allocated(self%order)     ) deallocate(self%order)
            self%exists = .false.
        endif
    end subroutine kill

end module simple_pspecs
