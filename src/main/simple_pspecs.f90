module simple_pspecs
include 'simple_lib.f08'
!$ use omp_lib
!$ use omp_lib_kinds
use simple_image,    only: image
use simple_fsc,      only: plot_fsc
use simple_oris,     only: oris
implicit none

public :: pspecs
private
#include "simple_local_flags.inc"

type pspecs
    private
    real, allocatable :: pspecs(:,:)  ! matrix of power spectra
    real    :: smpd       = 0.        ! sampling distance
    real    :: hp         = 0.        ! high-pass limit
    real    :: lp         = 0.        ! low-pass limit
    integer :: box        = 0         ! box size
    integer :: kfromto(2) = [0,0]     ! Fourier index range
    integer :: sz         = 0         ! size of spectrum 
    integer :: nspecs     = 0         ! # of spectra
    logical :: exists     = .false.   ! existence flag
contains
    ! constructor
    procedure          :: new
    ! getters & setters
    procedure          :: get_nspecs
    ! clustering
    procedure, private :: calc_distmat_1, calc_distmat_2
    generic            :: calc_distmat => calc_distmat_1, calc_distmat_2
    ! destructor
    procedure          :: kill
end type pspecs

contains

    ! constructor

    subroutine new( self, imgs, msk, hp, lp )
        class(pspecs), intent(inout) :: self
        class(image),  intent(inout) :: imgs(:)
        real,          intent(in)    :: msk, hp, lp
        type(image), allocatable :: cavg_threads(:) 
        real,        allocatable :: pspec(:)
        integer :: ldim(3), ispec, ithr
        call self%kill
        self%nspecs     = size(imgs)
        ldim            = imgs(1)%get_ldim()
        self%box        = ldim(1)
        self%smpd       = imgs(1)%get_smpd()
        self%hp         = hp
        self%lp         = lp
        self%kfromto(1) = calc_fourier_index(self%hp, self%box, self%smpd)
        self%kfromto(2) = calc_fourier_index(self%lp, self%box, self%smpd)
        self%sz         = self%kfromto(2) - self%kfromto(1) + 1
        ! allocate arrays
        allocate(self%pspecs(self%nspecs,self%sz), source=0.)
        ! fill-up the instance
        allocate(cavg_threads(nthr_glob))
        do ithr = 1, nthr_glob
            call cavg_threads(ithr)%new(ldim, self%smpd)
        end do
        !$omp parallel do default(shared) private(ispec,ithr,pspec) proc_bind(close) schedule(static)
        do ispec = 1, self%nspecs
            ithr = omp_get_thread_num() + 1
            ! power spectrum
            call cavg_threads(ithr)%copy(imgs(ispec))
            call cavg_threads(ithr)%norm
            call cavg_threads(ithr)%mask(msk, 'soft')
            call cavg_threads(ithr)%spectrum('sqrt', pspec)
            self%pspecs(ispec,:) = pspec(self%kfromto(1):self%kfromto(2))
            deallocate(pspec)
        end do
        !$omp end parallel do
        do ithr = 1, nthr_glob
            call cavg_threads(ithr)%kill
        end do
        deallocate(cavg_threads)
        self%exists = .true.
    end subroutine new

    ! getters & setters

    function get_nspecs( self ) result( nspecs )
        class(pspecs), intent(in) :: self
        integer :: nspecs
        nspecs = self%nspecs
    end function get_nspecs

    ! clustering

    function calc_distmat_1( self ) result( distmat )
        class(pspecs), intent(inout) :: self
        real, allocatable :: distmat(:,:)
        integer :: i, j
        allocate(distmat(self%nspecs,self%nspecs), source=0.)
        !$omp parallel do default(shared) private(i,j) proc_bind(close) schedule(dynamic)
        do i = 1, self%nspecs - 1
            do j = i + 1, self%nspecs
                distmat(i,j) = euclid(self%pspecs(i,:),self%pspecs(j,:))
                distmat(j,i) = distmat(i,j)
            end do
        end do
        !$omp end parallel do
    end function calc_distmat_1

    function calc_distmat_2( self1, self2 ) result( distmat )
        class(pspecs), intent(inout) :: self1, self2
        real, allocatable :: distmat(:,:)
        integer :: i, j
        allocate(distmat(self1%nspecs,self2%nspecs), source=0.)      
        !$omp parallel do default(shared) private(i,j) proc_bind(close) schedule(dynamic)
        do i = 1, self1%nspecs
            do j = 1, self2%nspecs
                distmat(i,j) = euclid(self1%pspecs(i,:),self2%pspecs(j,:))
            end do
        end do
        !$omp end parallel do
    end function calc_distmat_2 

    ! destructor

    subroutine kill( self )
        class(pspecs), intent(inout) :: self
        if( self%exists )then
            if( allocated(self%pspecs) ) deallocate(self%pspecs)
            self%exists = .false.
        endif
    end subroutine kill

end module simple_pspecs
