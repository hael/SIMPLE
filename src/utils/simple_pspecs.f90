module simple_pspecs
include 'simple_lib.f08'
!$ use omp_lib
!$ use omp_lib_kinds
use simple_image, only: image
use simple_fsc,   only: plot_fsc
implicit none

public :: pspecs
private
#include "simple_local_flags.inc"

type pspecs
    private  
    real,    allocatable :: pspecs(:,:)          ! matrix of power spectra
    real,    allocatable :: pspec_good(:)        ! 'good' average power spectrum
    real,    allocatable :: pspec_bad(:)         ! 'bad'  average power spectrum
    real,    allocatable :: resarr(:)            ! resolution values in A
    real,    allocatable :: dynranges(:)         ! dynamic spectral ranges
    real,    allocatable :: dists2good(:)        ! distance to good pspeccs
    real,    allocatable :: dists2bad(:)         ! distance to bad pspecs
    real,    allocatable :: distmat(:,:)         ! Euclidean distance matrix
    integer, allocatable :: ranks(:)             ! quality ranking
    integer, allocatable :: clsinds_spec(:)      ! 0: empty, 1: good, 2: bad
    integer, allocatable :: clsinds(:)           ! 2D class assignment
    integer, allocatable :: clspops(:)           ! class populations 
    real                 :: smpd                 ! sampling distance
    real                 :: hp                   ! high-pass limit
    real                 :: lp                   ! low-pass limit
    real                 :: dynrange_t = 0.      ! dynamic range threshold
    integer              :: box        = 0       ! box size
    integer              :: kfromto(2)           ! Fourier index range
    integer              :: sz         = 0       ! size of spectrum 
    integer              :: nspecs     = 0       ! # of spectra
    logical              :: exists     = .false. ! existence flag
contains
    procedure            :: new
    procedure            :: set_class_pspec
    procedure            :: bincls_pspecs
    procedure, private   :: bincls_dynrange
    procedure, private   :: calc_good_bad_pspecs
    procedure, private   :: kmeans_iter
    procedure            :: plot_good_bad
    procedure            :: rank_pspecs
    procedure            :: find_closest
    procedure            :: lookup_distance
    procedure            :: calc_distmat
    procedure            :: calc_avg_pspec
    procedure            :: kill
end type pspecs

integer, parameter :: CLASS_EMPTY = 0
integer, parameter :: CLASS_GOOD  = 1
integer, parameter :: CLASS_BAD   = 2
real,    parameter :: EMPTY_THRES = 1e-6
integer, parameter :: MAXITS      = 10

contains

    ! should be parallelized outside of here
    subroutine new( self, nspecs, img_template, hp, lp )
        class(pspecs), intent(inout) :: self
        integer,       intent(in)    :: nspecs
        class(image),  intent(in)    :: img_template 
        real,          intent(in)    :: hp, lp
        real, allocatable :: resarr(:)
        integer           :: ldim(3)
        call self%kill
        self%nspecs     = nspecs
        self%hp         = hp
        self%lp         = lp
        ldim            = img_template%get_ldim()
        self%box        = ldim(1)
        self%smpd       = img_template%get_smpd()
        resarr          = get_resarr(self%box, self%smpd)
        self%kfromto(1) = calc_fourier_index(self%hp, self%box, self%smpd)
        self%kfromto(2) = calc_fourier_index(self%lp, self%box, self%smpd)
        self%sz         = self%kfromto(2) - self%kfromto(1) + 1
        allocate(self%resarr(self%sz), source=resarr(self%kfromto(1):self%kfromto(2)))
        allocate(self%pspecs(self%nspecs,self%sz), self%dynranges(self%nspecs),&
        &self%dists2good(self%nspecs), self%dists2bad(self%nspecs), source=0.)
        allocate(self%ranks(self%nspecs), self%clsinds_spec(self%nspecs),&
        self%clsinds(self%nspecs), self%clspops(self%nspecs), source=0)
        deallocate(resarr)
        self%exists     = .true.
    end subroutine new

    ! should be parallelized outside of here
    subroutine set_class_pspec( self, ispec, class, pop, img, msk )
        class(pspecs), intent(inout) :: self
        integer,       intent(in)    :: ispec, class, pop
        class(image),  intent(inout) :: img
        real,          intent(in)    :: msk
        real, allocatable :: spec(:)
        if( ispec < 1 .or. ispec > self%nspecs ) THROW_HARD('ispec index out of range')
        if( pop == 0 )                           THROW_HARD('Empty classes not allowed')
        call img%norm
        call img%mask(msk, 'soft')
        call img%spectrum('sqrt', spec)
        self%pspecs(ispec,:)  = spec(self%kfromto(1):self%kfromto(2))
        self%dynranges(ispec) = self%pspecs(ispec,1) - self%pspecs(ispec,self%sz)
        if( self%dynranges(ispec) <= EMPTY_THRES ) THROW_HARD('Empty spectra not allowed')
        deallocate(spec)
        self%clsinds(ispec) = class
        self%clspops(ispec) = pop
    end subroutine set_class_pspec

    ! binary k-means clustering of power spectra
    subroutine bincls_pspecs( self )
        class(pspecs), intent(inout) :: self
        integer :: iter
        logical :: l_converged
        call self%bincls_dynrange
        call self%calc_good_bad_pspecs
        call self%plot_good_bad('pspec_good_dynrange', 'pspec_bad_dynrange')
        iter = 0
        l_converged = .false.
        do
            call self%kmeans_iter(iter, l_converged)
            if( l_converged ) exit
        end do
        call self%plot_good_bad('pspec_good_kmeans', 'pspec_bad_kmeans')
    end subroutine bincls_pspecs

    ! make initial grouping based on binary dynamic range clustering
    subroutine bincls_dynrange( self )
        class(pspecs), intent(inout) :: self
        integer :: ispec
        call otsu(self%nspecs, self%dynranges, self%dynrange_t)
        do ispec = 1, self%nspecs
            if( self%dynranges(ispec) >= self%dynrange_t )then
                self%clsinds_spec(ispec) = CLASS_GOOD
            else
                self%clsinds_spec(ispec) = CLASS_BAD
            endif
        enddo
    end subroutine bincls_dynrange

    subroutine calc_good_bad_pspecs( self )
        class(pspecs), intent(inout) :: self
        real, allocatable :: pspec_good(:), pspec_bad(:)
        integer :: ispec, ngood, nbad
        allocate(pspec_good(self%sz), pspec_bad(self%sz), source=0.)
        !$omp parallel do default(shared) private(ispec) proc_bind(close) reduction(+:pspec_good,pspec_bad)
        do ispec = 1, self%nspecs
            select case(self%clsinds_spec(ispec))
                case(CLASS_GOOD)
                    pspec_good = pspec_good + self%pspecs(ispec,:)
                case(CLASS_BAD)
                    pspec_bad  = pspec_bad  + self%pspecs(ispec,:)
            end select
        enddo
        !$omp end parallel do
        ngood = count(self%clsinds_spec == CLASS_GOOD)
        nbad  = count(self%clsinds_spec == CLASS_BAD)
        if( allocated(self%pspec_good) )then
            self%pspec_good = pspec_good / real(ngood)
        else
            allocate(self%pspec_good(self%sz), source=pspec_good / real(ngood))
        endif
        if( allocated(self%pspec_bad)  )then
            self%pspec_bad  = pspec_bad  / real(nbad)
        else
            allocate(self%pspec_bad(self%sz),  source=pspec_bad  / real(nbad))
        endif
    end subroutine calc_good_bad_pspecs

    subroutine kmeans_iter( self, iter, l_converged )
        class(pspecs), intent(inout) :: self
        integer,       intent(inout) :: iter
        logical,       intent(inout) :: l_converged
        integer :: nchanges, ispec
        nchanges = 0
        ! assign clusters
        !$omp parallel do default(shared) private(ispec) proc_bind(close)
        do ispec = 1,self%nspecs
            select case(self%clsinds_spec(ispec))
                case(CLASS_GOOD,CLASS_BAD)
                    self%dists2good(ispec) = euclid(self%pspec_good,self%pspecs(ispec,:))
                    self%dists2bad(ispec)  = euclid(self%pspec_bad, self%pspecs(ispec,:))
                    if( self%dists2good(ispec) <= self%dists2bad(ispec) )then ! is good
                        if( self%clsinds_spec(ispec) == CLASS_BAD ) nchanges = nchanges + 1
                        self%clsinds_spec(ispec) = CLASS_GOOD
                    else                                                      ! is bad
                        if( self%clsinds_spec(ispec) == CLASS_GOOD ) nchanges = nchanges + 1
                        self%clsinds_spec(ispec) = CLASS_BAD
                    endif
            end select
        end do
        !$omp end parallel do
        ! update averages
        call self%calc_good_bad_pspecs
        ! update iteration counter
        iter = iter + 1
        ! set l_converged flag
        l_converged = .false.
        if( nchanges == 0 .or. iter == MAXITS) l_converged = .true.
    end subroutine kmeans_iter

    subroutine plot_good_bad( self, fbody_good, fbody_bad )
        class(pspecs),    intent(in) :: self
        character(len=*), intent(in) :: fbody_good, fbody_bad
        call plot_fsc(self%sz, self%pspec_good, self%resarr, self%smpd, fbody_good)
        call plot_fsc(self%sz, self%pspec_bad,  self%resarr, self%smpd, fbody_bad)
    end subroutine plot_good_bad

    subroutine rank_pspecs( self )
        class(pspecs), intent(inout) :: self
        real,    allocatable :: dists_good(:), dists_bad(:)
        integer, allocatable :: inds_good(:), inds_bad(:)
        integer :: ngood, nbad, ispec, rank
        ngood = count(self%clsinds_spec == CLASS_GOOD)
        nbad  = count(self%clsinds_spec == CLASS_BAD)
        allocate(dists_good(ngood), dists_bad(nbad), source=0.)
        allocate(inds_good(ngood),  inds_bad(nbad),  source=0)
        ngood = 0
        nbad  = 0
        do ispec = 1, self%nspecs
            select case(self%clsinds_spec(ispec))
                case(CLASS_GOOD)
                    ngood             = ngood + 1
                    dists_good(ngood) = self%dists2good(ispec)
                    inds_good(ngood)  = ispec
                case(CLASS_BAD)
                    nbad              = nbad + 1
                    dists_bad(nbad)   = self%dists2good(ispec)
                    inds_bad(nbad)    = ispec
            end select
        enddo
        call hpsort(dists_good, inds_good)
        call hpsort(dists_bad,  inds_bad)
        rank = 0
        do ispec = 1, ngood
            rank = rank + 1
            self%ranks(inds_good(ispec)) = rank
        enddo
        do ispec = 1, nbad
            rank = rank + 1
            self%ranks(inds_bad(ispec)) = rank
        enddo
    end subroutine rank_pspecs

    function find_closest( self, ispec ) result( closest )
        class(pspecs), intent(inout) :: self
        integer,       intent(in)    :: ispec
        real    :: x, dists(self%nspecs)
        integer :: loc(1), closest, i
        do i = 1, self%nspecs
            if( i == ispec )then
                dists(i) = huge(x)
            else
                dists(i) = self%lookup_distance(i, ispec)
            endif
        end do
        loc     = minloc(dists)
        closest = loc(1)
    end function find_closest

    function lookup_distance( self, i, j ) result( d )
        class(pspecs), intent(inout) :: self
        integer,       intent(in)    :: i, j
        real :: d
        if( allocated(self%distmat) )then
            d = self%distmat(i,j)
        else
            call self%calc_distmat
        endif
    end function lookup_distance

    subroutine calc_distmat( self )
        class(pspecs), intent(inout) :: self
        integer :: i, j
        if( allocated(self%distmat) )then
            self%distmat = 0.
        else
            allocate(self%distmat(self%nspecs,self%nspecs), source=0.)
        endif
        !$omp parallel do default(shared) private(i,j) proc_bind(close) schedule(dynamic)
        do i = 1, self%nspecs - 1
            do j = i + 1, self%nspecs
                self%distmat(i,j) = euclid(self%pspecs(i,:),self%pspecs(j,:))
                self%distmat(j,i) = self%distmat(i,j)
            end do
        end do
        !$omp end parallel do
    end subroutine calc_distmat

    function calc_avg_pspec( self, mask ) result( pspec_avg )
        class(pspecs), intent(in) :: self
        logical,       intent(in) :: mask(self%nspecs)
        integer           :: ispec, cnt
        real, allocatable :: pspec_avg(:)
        cnt = count(mask)
        if( cnt == 0 ) THROW_HARD('mask empty')
        allocate(pspec_avg(self%sz), source=0.)
        !$omp parallel do default(shared) private(ispec) proc_bind(close) reduction(+:pspec_avg) schedule(static)
        do ispec = 1, self%nspecs
            if( mask(ispec) )then
                pspec_avg = pspec_avg + self%pspecs(ispec,:)
            endif
        enddo
        !$omp end parallel do
        pspec_avg = pspec_avg / real(cnt)
    end function calc_avg_pspec

    subroutine kill( self )
        class(pspecs), intent(inout) :: self
        if( self%exists )then
            if( allocated(self%pspecs)       ) deallocate(self%pspecs)
            if( allocated(self%resarr)       ) deallocate(self%resarr)
            if( allocated(self%dynranges)    ) deallocate(self%dynranges)
            if( allocated(self%dists2good)   ) deallocate(self%dists2good)
            if( allocated(self%dists2bad)    ) deallocate(self%dists2bad)
            if( allocated(self%ranks)        ) deallocate(self%ranks)
            if( allocated(self%clsinds_spec) ) deallocate(self%clsinds_spec)
            if( allocated(self%clsinds)      ) deallocate(self%clsinds)
            if( allocated(self%clspops)      ) deallocate(self%clspops)
            self%exists = .false.
        endif
    end subroutine kill

end module simple_pspecs
