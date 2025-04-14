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
    real,    allocatable :: pspecs(:,:)           ! matrix of power spectra
    real,    allocatable :: pspec_good(:)         ! 'good' average power spectrum
    real,    allocatable :: pspec_bad(:)          ! 'bad'  average power spectrum
    real,    allocatable :: resarr(:)             ! resolution values in A
    real,    allocatable :: dynranges(:)          ! dynamic ranges of power spectra
    real,    allocatable :: dists2good(:)         ! distance to good pspeccs
    real,    allocatable :: dists2bad(:)          ! distance to bad  pspecs
    real,    allocatable :: distmat(:,:)          ! Euclidean distance matrix
    integer, allocatable :: ranks(:)              ! quality ranking
    integer, allocatable :: order(:)              ! quality order
    integer, allocatable :: clsinds_spec(:)       ! 0: empty, 1: good, 2: bad
    integer, allocatable :: clsinds(:)            ! 2D class assignment
    integer, allocatable :: clspops(:)            ! class populations 
    real                 :: smpd                  ! sampling distance
    real                 :: hp                    ! high-pass limit
    real                 :: lp                    ! low-pass limit
    real                 :: dynrange_t  = 0.      ! dynamic range threshold
    integer              :: box         = 0       ! box size
    integer              :: kfromto(2)            ! Fourier index range
    integer              :: sz          = 0       ! size of spectrum 
    integer              :: nspecs      = 0       ! # of spectra
    integer              :: medoid_good = 0       ! index of medoid of good class
    integer              :: medoid_bad  = 0       ! index of medoid of bad  class
    integer              :: ngood       = 0       ! # good power spectra
    integer              :: nbad        = 0       ! # bad power spectra
    logical              :: exists      = .false. ! existence flag
contains
    
    procedure, private   :: new_1, new_2
    generic              :: new => new_1, new_2
    procedure, private   :: set_class_pspec_1, set_class_pspec_2
    generic              :: set_class_pspec => set_class_pspec_1, set_class_pspec_2
    procedure            :: kmeans_bincls_pspecs_and_rank
    procedure            :: kmedoids_bincls_pspecs_and_rank
    procedure            :: plot_good_bad_avg
    procedure            :: plot_good_bad_medoid
    procedure            :: plot_all
    procedure            :: plot_all_ranked
    procedure            :: calc_frac_good
    procedure            :: get_ngood
    procedure            :: set_ngood
    procedure            :: get_good_bad_state_arr
    procedure            :: get_ordered_clsind
    procedure, private   :: bincls_dynrange
    procedure, private   :: calc_good_bad_pspec_avgs
    procedure, private   :: find_good_bad_pspec_medoids
    procedure, private   :: kmeans_iter
    procedure, private   :: kmedoids_iter
    procedure, private   :: rank_pspecs
    procedure, private   :: find_closest
    procedure, private   :: lookup_distance
    procedure, private   :: calc_distmat
    procedure, private   :: calc_avg_pspec
    procedure            :: kill
end type pspecs

integer, parameter :: CLASS_EMPTY = 0
integer, parameter :: CLASS_GOOD  = 1
integer, parameter :: CLASS_BAD   = 2
real,    parameter :: EMPTY_THRES = 1e-6
integer, parameter :: MAXITS      = 10

contains

    ! should be parallelized outside of here
    subroutine new_1( self, nspecs, img_template, hp, lp )
        class(pspecs), intent(inout) :: self
        integer,       intent(in)    :: nspecs
        class(image),  intent(in)    :: img_template 
        real,          intent(in)    :: hp, lp
        real, allocatable :: resarr(:)
        integer           :: ldim(3)
        call self%kill
        ldim = img_template%get_ldim()
        call self%new_2(nspecs, ldim(1), img_template%get_smpd(), hp, lp )
    end subroutine new_1

    subroutine new_2( self, nspecs, box, smpd, hp, lp )
        class(pspecs), intent(inout) :: self
        integer,       intent(in)    :: nspecs
        integer,       intent(in)    :: box
        real,          intent(in)    :: smpd, hp, lp
        real, allocatable :: resarr(:)
        call self%kill
        self%nspecs     = nspecs
        self%hp         = hp
        self%lp         = lp
        self%box        = box
        self%smpd       = smpd
        resarr          = get_resarr(self%box, self%smpd)
        self%kfromto(1) = calc_fourier_index(self%hp, self%box, self%smpd)
        self%kfromto(2) = calc_fourier_index(self%lp, self%box, self%smpd)
        self%sz         = self%kfromto(2) - self%kfromto(1) + 1
        allocate(self%resarr(self%sz), source=resarr(self%kfromto(1):self%kfromto(2)))
        allocate(self%pspecs(self%nspecs,self%sz), self%dynranges(self%nspecs),&
        &self%dists2good(self%nspecs), self%dists2bad(self%nspecs), source=0.)
        allocate(self%ranks(self%nspecs), self%order(self%nspecs), self%clsinds_spec(self%nspecs),&
        self%clsinds(self%nspecs), self%clspops(self%nspecs), source=0)
        deallocate(resarr)
        self%exists     = .true.
    end subroutine new_2

    ! should be parallelized outside of here
    subroutine set_class_pspec_1( self, ispec, class, pop, img, msk )
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
    end subroutine set_class_pspec_1

    subroutine set_class_pspec_2( self, ispec, class, pop, spec )
        class(pspecs), intent(inout) :: self
        integer,       intent(in)    :: ispec, class, pop
        real,          intent(in)    :: spec(:)
        if( ispec < 1 .or. ispec > self%nspecs ) THROW_HARD('ispec index out of range')
        if( pop == 0 )                           THROW_HARD('Empty classes not allowed')
        self%pspecs(ispec,:)  = spec(self%kfromto(1):self%kfromto(2))
        self%dynranges(ispec) = self%pspecs(ispec,1) - self%pspecs(ispec,self%sz)
        if( self%dynranges(ispec) <= EMPTY_THRES ) THROW_HARD('Empty spectra not allowed')
        self%clsinds(ispec) = class
        self%clspops(ispec) = pop
    end subroutine set_class_pspec_2

    ! binary k-means clustering of power spectra
    subroutine kmeans_bincls_pspecs_and_rank( self )
        class(pspecs), intent(inout) :: self
        integer :: iter
        logical :: l_converged
        call self%bincls_dynrange
        call self%calc_good_bad_pspec_avgs
        call self%plot_good_bad_avg('pspec_good_dynrange', 'pspec_bad_dynrange')
        iter = 0
        l_converged = .false.
        do
            call self%kmeans_iter(iter, l_converged)
            if( l_converged ) exit
        end do
        call self%plot_good_bad_avg('pspec_good_kmeans', 'pspec_bad_kmeans')
        call self%rank_pspecs
    end subroutine kmeans_bincls_pspecs_and_rank

    ! binary k-medoids clustering of power spectra
    subroutine kmedoids_bincls_pspecs_and_rank( self )
        class(pspecs), intent(inout) :: self
        integer :: iter
        logical :: l_converged
        call self%bincls_dynrange
        call self%calc_distmat
        call self%find_good_bad_pspec_medoids
        call self%plot_good_bad_medoid('pspec_good_dynrange', 'pspec_bad_dynrange')
        iter = 0
        l_converged = .false.
        do
            call self%kmedoids_iter(iter, l_converged)
            if( l_converged ) exit
        end do
        call self%plot_good_bad_medoid('pspec_good_kmedoids', 'pspec_bad_kmedoids')
        call self%rank_pspecs
    end subroutine kmedoids_bincls_pspecs_and_rank

    subroutine plot_good_bad_avg( self, fbody_good, fbody_bad )
        class(pspecs),    intent(in) :: self
        character(len=*), intent(in) :: fbody_good, fbody_bad
        call plot_fsc(self%sz, self%pspec_good, self%resarr, self%smpd, fbody_good)
        call plot_fsc(self%sz, self%pspec_bad,  self%resarr, self%smpd, fbody_bad)
    end subroutine plot_good_bad_avg

    subroutine plot_good_bad_medoid( self, fbody_good, fbody_bad )
        class(pspecs),    intent(in) :: self
        character(len=*), intent(in) :: fbody_good, fbody_bad
        call plot_fsc(self%sz, self%pspecs(self%medoid_good,:), self%resarr, self%smpd, fbody_good)
        call plot_fsc(self%sz, self%pspecs(self%medoid_bad,:),  self%resarr, self%smpd, fbody_bad)
    end subroutine plot_good_bad_medoid

    subroutine plot_all( self )
        class(pspecs),     intent(in) :: self
        character(len=:), allocatable :: fbody
        integer :: ispec
        do ispec = 1, self%nspecs
            fbody = 'power_spectrum'//int2str_pad(ispec,3)
            call plot_fsc(self%sz, self%pspecs(ispec,:), self%resarr, self%smpd, fbody)
        end do
    end subroutine plot_all

    subroutine plot_all_ranked( self )
        class(pspecs),     intent(in) :: self
        character(len=:), allocatable :: fbody
        integer :: ispec
        do ispec = 1, self%nspecs
            fbody = 'power_spectrum_rank'//int2str_pad(self%ranks(ispec),3)
            call plot_fsc(self%sz, self%pspecs(self%order(ispec),:), self%resarr, self%smpd, fbody)
        end do
    end subroutine plot_all_ranked

    function calc_frac_good( self, ngood, nptcls ) result( frac_good )
        class(pspecs), intent(in) :: self
        integer,       intent(in) :: ngood, nptcls
        real :: frac_good
        frac_good = (real(sum(self%clspops, mask=self%ranks <= ngood))/real(nptcls))
    end function calc_frac_good

    function get_ngood( self ) result( ngood )
        class(pspecs), intent(in) :: self
        integer :: ngood
        ngood = self%ngood
    end function get_ngood

    subroutine set_ngood( self, ngood )
        class(pspecs), intent(inout) :: self
        integer,       intent(in)    :: ngood
        self%ngood = ngood
        self%nbad  = self%nspecs - ngood
    end subroutine set_ngood

    function get_good_bad_state_arr( self, ncls ) result( states )
        class(pspecs), intent(in) :: self
        integer,       intent(in) :: ncls
        integer, allocatable :: states(:)
        integer :: ispec
        allocate(states(ncls), source=0)
        do ispec = 1, self%nspecs
            if( self%ranks(ispec) <= self%ngood ) states(self%clsinds(ispec)) = 1
        end do
    end function get_good_bad_state_arr

    function get_ordered_clsind( self, ispec ) result( clsind )
        class(pspecs), intent(in) :: self
        integer,       intent(in) :: ispec
        integer :: clsind
        if( ispec < 1 .or. ispec > self%nspecs ) THROW_HARD('ispec index out of range')
        clsind = self%clsinds(self%order(ispec))
    end function get_ordered_clsind

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

    subroutine calc_good_bad_pspec_avgs( self )
        class(pspecs), intent(inout) :: self
        real, allocatable :: pspec_good(:), pspec_bad(:)
        integer :: ispec, ngood, nbad
        allocate(pspec_good(self%sz), pspec_bad(self%sz), source=0.)
        !$omp parallel do default(shared) private(ispec) proc_bind(close) reduction(+:pspec_good,pspec_bad) schedule(static)
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
    end subroutine calc_good_bad_pspec_avgs

    subroutine find_good_bad_pspec_medoids( self )
        class(pspecs), intent(inout) :: self
        real    :: dists_good(self%nspecs), dists_bad(self%nspecs)
        integer :: i, j, loc(1)
        dists_good = 0.
        dists_bad  = 0.
        !$omp parallel do default(shared) private(i,j) proc_bind(close) schedule(static)
        do i = 1, self%nspecs
            do j = 1, self%nspecs
                if( i /= j )then
                    if(      self%clsinds_spec(i) == CLASS_GOOD .and. self%clsinds_spec(j) == CLASS_GOOD )then
                        dists_good(i) = dists_good(i) + self%lookup_distance(i, j)
                    else if( self%clsinds_spec(i) == CLASS_BAD  .and. self%clsinds_spec(j) == CLASS_BAD  )then
                        dists_bad(i)  = dists_bad(i)  + self%lookup_distance(i, j)
                    endif
                endif
            end do
        end do
        !$omp end parallel do
        loc = minloc(dists_good)
        self%medoid_good = loc(1)
        loc = minloc(dists_bad)
        self%medoid_bad = loc(1)
    end subroutine find_good_bad_pspec_medoids

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
        call self%calc_good_bad_pspec_avgs
        ! update iteration counter
        iter = iter + 1
        ! set l_converged flag
        l_converged = .false.
        if( nchanges == 0 .or. iter == MAXITS) l_converged = .true.
    end subroutine kmeans_iter

    subroutine kmedoids_iter( self, iter, l_converged )
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
                    self%dists2good(ispec) = euclid(self%pspecs(self%medoid_good,:),self%pspecs(ispec,:))
                    self%dists2bad(ispec)  = euclid(self%pspecs(self%medoid_good,:), self%pspecs(ispec,:))
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
        ! update medoids
        call self%calc_distmat
        call self%find_good_bad_pspec_medoids
        ! update iteration counter
        iter = iter + 1
        ! set l_converged flag
        l_converged = .false.
        if( nchanges == 0 .or. iter == MAXITS) l_converged = .true.
    end subroutine kmedoids_iter

    subroutine rank_pspecs( self )
        class(pspecs), intent(inout) :: self
        real,    allocatable :: dists_good(:), dists_bad(:)
        integer, allocatable :: inds_good(:), inds_bad(:)
        integer :: ispec, rank
        self%ngood = count(self%clsinds_spec == CLASS_GOOD)
        self%nbad  = count(self%clsinds_spec == CLASS_BAD)
        allocate(dists_good(self%ngood), dists_bad(self%nbad), source=0.)
        allocate(inds_good(self%ngood),  inds_bad(self%nbad),  source=0)
        self%ngood = 0
        self%nbad  = 0
        do ispec = 1, self%nspecs
            select case(self%clsinds_spec(ispec))
                case(CLASS_GOOD)
                    self%ngood             = self%ngood + 1
                    dists_good(self%ngood) = self%dists2good(ispec)
                    inds_good(self%ngood)  = ispec
                case(CLASS_BAD)
                    self%nbad              = self%nbad + 1
                    dists_bad(self%nbad)   = self%dists2good(ispec)
                    inds_bad(self%nbad)    = ispec
            end select
        enddo
        call hpsort(dists_good, inds_good)
        call hpsort(dists_bad,  inds_bad)
        rank = 0
        do ispec = 1, self%ngood
            rank = rank + 1
            self%ranks(inds_good(ispec)) = rank
            self%order(rank)             = inds_good(ispec)
        enddo
        do ispec = 1, self%nbad
            rank = rank + 1
            self%ranks(inds_bad(ispec)) = rank
            self%order(rank)            = inds_bad(ispec)
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
            d = self%distmat(i,j)
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
            if( allocated(self%order)        ) deallocate(self%order)
            if( allocated(self%clsinds_spec) ) deallocate(self%clsinds_spec)
            if( allocated(self%clsinds)      ) deallocate(self%clsinds)
            if( allocated(self%clspops)      ) deallocate(self%clspops)
            self%exists = .false.
        endif
    end subroutine kill

end module simple_pspecs
