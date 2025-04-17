module simple_pspecs
include 'simple_lib.f08'
!$ use omp_lib
!$ use omp_lib_kinds
use simple_image, only: image
use simple_fsc,   only: plot_fsc
use simple_oris,  only: oris
implicit none

public :: pspecs
private
#include "simple_local_flags.inc"

type pspecs
    private  
    real,               allocatable :: pspecs(:,:)       ! matrix of power spectra
    real,               allocatable :: pspec_good(:)     ! 'good' average power spectrum
    real,               allocatable :: pspec_bad(:)      ! 'bad'  average power spectrum
    real,               allocatable :: resarr(:)         ! resolution values in A
    real,               allocatable :: dynranges(:)      ! dynamic ranges of power spectra
    real,               allocatable :: dists2good(:)     ! distance to good pspeccs
    real,               allocatable :: dists2bad(:)      ! distance to bad  pspecs
    real,               allocatable :: distmat(:,:)      ! Euclidean distance matrix
    integer,            allocatable :: ranks(:)          ! quality ranking
    integer,            allocatable :: order(:)          ! quality order
    integer,            allocatable :: clsinds_spec(:)   ! 0: empty, 1: good, 2: bad
    integer,            allocatable :: clsinds(:)        ! 2D class assignment
    integer,            allocatable :: clspops(:)        ! class populations
    type(stats_struct), allocatable :: clsscore_stats(:) ! class score statistics
    real                 :: smpd                         ! sampling distance
    real                 :: hp                           ! high-pass limit
    real                 :: lp                           ! low-pass limit
    real                 :: dynrange_t  = 0.             ! dynamic range threshold
    integer              :: box         = 0              ! box size
    integer              :: kfromto(2)                   ! Fourier index range
    integer              :: sz          = 0              ! size of spectrum 
    integer              :: nspecs      = 0              ! # of spectra
    integer              :: ncls        = 0              ! # classes
    integer              :: medoid_good = 0              ! index of medoid of good class
    integer              :: medoid_bad  = 0              ! index of medoid of bad  class
    integer              :: ngood       = 0              ! # good power spectra
    integer              :: nbad        = 0              ! # bad power spectra
    logical              :: exists      = .false.        ! existence flag
contains
    ! constructor
    procedure            :: new
    ! getters & setters
    procedure            :: get_nspecs
    procedure            :: get_frac_good
    procedure            :: get_ngood
    procedure            :: set_ngood
    procedure            :: get_good_bad_state_arr
    procedure            :: get_ordered_clsind
    ! plotting
    procedure            :: plot_all
    procedure            :: plot_good_bad
    procedure            :: plot_all_ranked
    ! clustering 
    procedure            :: otsu_bincls_dynrange
    procedure            :: kmeans_bincls_pspecs_and_rank
    procedure            :: kmedoids_bincls_pspecs_and_rank
    procedure            :: greedy_bincls_pspecs_and_rank
    procedure            :: hybrid_bincls_pspecs_and_rank
    ! calculators
    procedure            :: smoothen_spectra
    procedure            :: median_good_clspop
    procedure, private   :: calc_good_bad_pspec_avgs
    procedure, private   :: find_good_bad_pspec_medoids
    procedure, private   :: kcluster_iter
    procedure, private   :: greedy_min_iter
    procedure, private   :: rank_pspecs
    procedure, private   :: clsind_spec_of_closest
    procedure, private   :: find_closest
    procedure, private   :: lookup_distance
    procedure, private   :: calc_distmat
    ! destructor
    procedure            :: kill
end type pspecs

integer, parameter :: CLASS_GOOD     = 1
integer, parameter :: CLASS_BAD      = 2
real,    parameter :: DYNRANGE_THRES = 1e-6
integer, parameter :: MAXITS         = 10
integer, parameter :: MINPOP         = 20

contains

    ! constructor

    subroutine new( self, ncls, imgs, os_ptcl2D, os_cls2D, msk, hp, lp )
        class(pspecs), intent(inout) :: self
        integer,       intent(in)    :: ncls
        class(image),  intent(inout) :: imgs(ncls)
        class(oris),   intent(inout) :: os_ptcl2D
        class(oris),   intent(in)    :: os_cls2D
        real,          intent(in)    :: msk, hp, lp
        real,    allocatable :: pspec(:), resarr(:)
        logical, allocatable :: mask(:)
        integer :: specinds(ncls)
        logical :: l_valid_spectra(ncls)
        integer :: ldim(3), icls, ispec
        real    :: dynrange
        call self%kill
        self%ncls       = ncls
        ldim            = imgs(1)%get_ldim()
        self%box        = ldim(1)
        self%smpd       = imgs(1)%get_smpd()
        self%hp         = hp
        self%lp         = lp
        resarr          = get_resarr(self%box, self%smpd)
        self%kfromto(1) = calc_fourier_index(self%hp, self%box, self%smpd)
        self%kfromto(2) = calc_fourier_index(self%lp, self%box, self%smpd)
        self%sz         = self%kfromto(2) - self%kfromto(1) + 1
        ! count valid spectra
        l_valid_spectra = .false.
        !$omp parallel do default(shared) private(icls, pspec, dynrange) proc_bind(close) schedule(static)
        do icls = 1, ncls
            call imgs(icls)%norm
            call imgs(icls)%mask(msk, 'soft')
            call imgs(icls)%spectrum('sqrt', pspec)
            dynrange = pspec(self%kfromto(1)) - pspec(self%kfromto(2))
            if( dynrange > DYNRANGE_THRES .and. os_cls2D%get_int(icls, 'pop') >= MINPOP  ) l_valid_spectra(icls) = .true.
        enddo
        !$omp end parallel do
        self%nspecs = count(l_valid_spectra)
        ! allocate arrays
        allocate(self%resarr(self%sz), source=resarr(self%kfromto(1):self%kfromto(2)))
        allocate(self%pspecs(self%nspecs,self%sz), self%dynranges(self%nspecs),&
        &self%dists2good(self%nspecs), self%dists2bad(self%nspecs), source=0.)
        allocate(self%ranks(self%nspecs), self%order(self%nspecs), self%clsinds_spec(self%nspecs),&
        self%clsinds(self%nspecs), self%clspops(self%nspecs), source=0)
        ! set spectrum indices
        ispec = 0
        do icls = 1, ncls
            if( l_valid_spectra(icls) )then
                ispec = ispec + 1
                specinds(icls) = ispec
            endif
        end do
        ! fill-up the instance
        !$omp parallel do default(shared) private(icls,pspec,mask) proc_bind(close) schedule(static)
        do icls = 1, ncls
            if( l_valid_spectra(icls) )then
                self%clsinds(specinds(icls))   = icls
                call imgs(icls)%spectrum('sqrt', pspec)
                self%pspecs(specinds(icls),:)  = pspec(self%kfromto(1):self%kfromto(2))
                self%dynranges(specinds(icls)) = self%pspecs(specinds(icls),1) - self%pspecs(specinds(icls),self%sz)
                self%clspops(specinds(icls))   = os_cls2D%get_int(icls, 'pop')
                mask = os_ptcl2D%gen_ptcl_mask('class', icls)
                call os_ptcl2D%stats('corr', clsstats, mask)
                deallocate(mask, pspec)
            endif
        end do
        !$omp end parallel do
        deallocate(resarr)
        self%exists = .true.
    end subroutine new

    ! getters & setters

    function get_nspecs( self ) result( nspecs )
        class(pspecs), intent(in) :: self
        integer :: nspecs
        nspecs = self%nspecs
    end function get_nspecs

    function get_frac_good( self, ngood, nptcls ) result( frac_good )
        class(pspecs), intent(in) :: self
        integer,       intent(in) :: ngood, nptcls
        real :: frac_good
        frac_good = (real(sum(self%clspops, mask=self%ranks <= ngood))/real(nptcls))
    end function get_frac_good

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

    function get_good_bad_state_arr( self ) result( states )
        class(pspecs), intent(in) :: self
        integer, allocatable :: states(:)
        integer :: ispec
        allocate(states(self%ncls), source=0)
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

    ! plotting

    subroutine plot_good_bad( self, fbody_good, fbody_bad )
        class(pspecs),    intent(in) :: self
        character(len=*), intent(in) :: fbody_good, fbody_bad
        call plot_fsc(self%sz, self%pspec_good, self%resarr, self%smpd, fbody_good)
        call plot_fsc(self%sz, self%pspec_bad,  self%resarr, self%smpd, fbody_bad)
    end subroutine plot_good_bad

    subroutine plot_all( self, fbody )
        class(pspecs),               intent(in) :: self
        character(len=*),  optional, intent(in) :: fbody
        character(len=:), allocatable :: ffbody, fffbody
        integer :: ispec
        if( present(fbody) )then
            allocate(ffbody, source=trim(fbody))
        else
            allocate(ffbody, source='power_spectrum')
        endif
        do ispec = 1, self%nspecs
            fffbody = ffbody//int2str_pad(ispec,3)
            call plot_fsc(self%sz, self%pspecs(ispec,:), self%resarr, self%smpd, fffbody)
        end do
    end subroutine plot_all

    subroutine plot_all_ranked( self, fbody )
        class(pspecs),               intent(in) :: self
        character(len=*),  optional, intent(in) :: fbody
        character(len=:), allocatable :: ffbody, fffbody
        integer :: ispec
        if( present(fbody) )then
            allocate(ffbody, source=trim(fbody))
        else
            allocate(ffbody, source='power_spectrum')
        endif
        do ispec = 1, self%nspecs
            fffbody = ffbody//'_rank'//int2str_pad(self%ranks(ispec),3)
            call plot_fsc(self%sz, self%pspecs(self%order(ispec),:), self%resarr, self%smpd, fffbody)
        end do
    end subroutine plot_all_ranked

    ! clustering

    ! make initial grouping based on binary dynamic range clustering
    subroutine otsu_bincls_dynrange( self )
        class(pspecs), intent(inout) :: self
        real, allocatable :: tmp(:)
        integer :: ispec, rank
        tmp = pack(self%dynranges, self%dynranges > DYNRANGE_THRES)
        call otsu(size(tmp), tmp, self%dynrange_t)
        do ispec = 1, self%nspecs
            if( self%dynranges(ispec) >= self%dynrange_t )then
                self%clsinds_spec(ispec) = CLASS_GOOD
            else
                self%clsinds_spec(ispec) = CLASS_BAD
            endif
        enddo
        deallocate(tmp)
    end subroutine otsu_bincls_dynrange

    subroutine hybrid_bincls_pspecs_and_rank( self )
        class(pspecs), intent(inout) :: self
        integer :: iter
        logical :: l_converged
        write(logfhandle,'(A)') 'HYBRID GREEDY/K-MEANS CLUSTERING OF POWERSPECTRA'
        ! start with k-medoids, as it seems better at finding a solid bad group
        call self%otsu_bincls_dynrange
        call self%calc_distmat
        call self%calc_good_bad_pspec_avgs
        call self%plot_good_bad('pspec_good_dynrange', 'pspec_bad_dynrange')
        iter = 0
        l_converged = .false.
        do
            call self%greedy_min_iter(iter, l_converged)
            if( l_converged ) exit
        end do
        call self%calc_good_bad_pspec_avgs
        call self%plot_good_bad('pspec_good_greedy', 'pspec_bad_greedy')
        ! refine using k-means
        iter = 0
        l_converged = .false.
        do
            call self%kcluster_iter(iter, l_converged, l_medoid=.false.)
            if( l_converged ) exit
        end do
        call self%plot_good_bad('pspec_good_kmeans', 'pspec_bad_kmeans')
        ! rank using means
        call self%rank_pspecs
    end subroutine hybrid_bincls_pspecs_and_rank

    ! binary k-means clustering of power spectra
    subroutine kmeans_bincls_pspecs_and_rank( self )
        class(pspecs), intent(inout) :: self
        integer :: iter
        logical :: l_converged
        write(logfhandle,'(A)') 'K-MEANS CLUSTERING OF POWERSPECTRA'
        call self%otsu_bincls_dynrange
        call self%calc_good_bad_pspec_avgs
        call self%plot_good_bad('pspec_good_dynrange', 'pspec_bad_dynrange')
        iter = 0
        l_converged = .false.
        do
            call self%kcluster_iter(iter, l_converged, l_medoid=.false.)
            if( l_converged ) exit
        end do
        call self%plot_good_bad('pspec_good_kmeans', 'pspec_bad_kmeans')
        call self%rank_pspecs
    end subroutine kmeans_bincls_pspecs_and_rank

    ! binary k-medoids clustering of power spectra
    subroutine kmedoids_bincls_pspecs_and_rank( self )
        class(pspecs), intent(inout) :: self
        integer :: iter
        logical :: l_converged
        call self%otsu_bincls_dynrange
        call self%calc_distmat
        call self%find_good_bad_pspec_medoids
        call self%plot_good_bad('pspec_good_dynrange', 'pspec_bad_dynrange')
        iter = 0
        l_converged = .false.
        do
            call self%kcluster_iter(iter, l_converged, l_medoid=.true.)
            if( l_converged ) exit
        end do
        call self%plot_good_bad('pspec_good_kmedoids', 'pspec_bad_kmedoids')
        call self%rank_pspecs
    end subroutine kmedoids_bincls_pspecs_and_rank

    subroutine greedy_bincls_pspecs_and_rank( self )
        class(pspecs), intent(inout) :: self
        integer :: iter
        logical :: l_converged
        call self%otsu_bincls_dynrange
        call self%calc_distmat
        call self%find_good_bad_pspec_medoids
        call self%plot_good_bad('pspec_good_dynrange', 'pspec_bad_dynrange')
        iter = 0
        l_converged = .false.
        do
            call self%greedy_min_iter(iter, l_converged)
            if( l_converged ) exit
        end do
        call self%find_good_bad_pspec_medoids
        call self%plot_good_bad('pspec_good_kmedoids', 'pspec_bad_kmedoids')
        call self%rank_pspecs
    end subroutine greedy_bincls_pspecs_and_rank

    ! calculators

    subroutine smoothen_spectra( self )
        class(pspecs), intent(inout) :: self
        integer :: ispec
        !$omp parallel do default(shared) private(ispec) proc_bind(close) schedule(static)
        do ispec = 1, self%nspecs
            call SavitzkyGolay_filter(self%sz, self%pspecs(ispec,:))
        enddo
        !$omp end parallel do
    end subroutine smoothen_spectra

    function median_good_clspop( self ) result( med )
        class(pspecs), intent(in) :: self
        real, allocatable :: clspops_good(:)
        real :: med
        clspops_good = real(pack(self%clspops, mask=self%clsinds_spec == CLASS_GOOD))
        med = median(clspops_good)
    end function median_good_clspop

    subroutine calc_good_bad_pspec_avgs( self )
        class(pspecs), intent(inout) :: self
        real, allocatable :: pspec_good(:), pspec_bad(:)
        integer :: ispec
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
        self%ngood = count(self%clsinds_spec == CLASS_GOOD)
        self%nbad  = count(self%clsinds_spec == CLASS_BAD)
        if( allocated(self%pspec_good) )then
            self%pspec_good = pspec_good / real(self%ngood)
        else
            allocate(self%pspec_good(self%sz), source=pspec_good / real(self%ngood))
        endif
        if( allocated(self%pspec_bad)  )then
            self%pspec_bad  = pspec_bad  / real(self%nbad)
        else
            allocate(self%pspec_bad(self%sz),  source=pspec_bad  / real(self%nbad))
        endif
    end subroutine calc_good_bad_pspec_avgs

    subroutine find_good_bad_pspec_medoids( self )
        class(pspecs), intent(inout) :: self
        real    :: dists_good(self%nspecs), dists_bad(self%nspecs)
        integer :: i, j, loc(1)
        !$omp parallel do default(shared) private(i,j) proc_bind(close) schedule(static)
        do i = 1, self%nspecs
            dists_good(i) = 0.
            dists_bad(i)  = 0.
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
        self%pspec_good  = self%pspecs(self%medoid_good,:)
        loc = minloc(dists_bad)
        self%medoid_bad  = loc(1)
        self%pspec_bad   = self%pspecs(self%medoid_bad,:)
    end subroutine find_good_bad_pspec_medoids

    subroutine kcluster_iter( self, iter, l_converged, l_medoid )
        class(pspecs), intent(inout) :: self
        integer,       intent(inout) :: iter
        logical,       intent(inout) :: l_converged
        logical,       intent(in)    :: l_medoid
        integer :: nchanges, ispec
        nchanges = 0
        ! assign clusters
        !$omp parallel do default(shared) private(ispec) proc_bind(close)
        do ispec = 1,self%nspecs
            self%dists2good(ispec) = euclid(self%pspec_good,self%pspecs(ispec,:))
            self%dists2bad(ispec)  = euclid(self%pspec_bad, self%pspecs(ispec,:))
            if( self%dists2good(ispec) <= self%dists2bad(ispec) )then ! is good
                if( self%clsinds_spec(ispec) == CLASS_BAD ) nchanges = nchanges + 1
                self%clsinds_spec(ispec) = CLASS_GOOD
            else                                                      ! is bad
                if( self%clsinds_spec(ispec) == CLASS_GOOD ) nchanges = nchanges + 1
                self%clsinds_spec(ispec) = CLASS_BAD
            endif
        end do
        !$omp end parallel do
        if( l_medoid )then
            call self%find_good_bad_pspec_medoids
        else
            call self%calc_good_bad_pspec_avgs
        endif
        ! update iteration counter
        iter = iter + 1
        ! set l_converged flag
        l_converged = .false.
        if( nchanges == 0 .or. iter == MAXITS) l_converged = .true.
    end subroutine kcluster_iter

    subroutine greedy_min_iter( self, iter, l_converged )
        class(pspecs), intent(inout) :: self
        integer,       intent(inout) :: iter
        logical,       intent(inout) :: l_converged
        integer :: nchanges, ispec
        nchanges = 0
        ! assign clusters
        !$omp parallel do default(shared) private(ispec) proc_bind(close)
        do ispec = 1,self%nspecs
            if( self%clsind_spec_of_closest(ispec) == CLASS_GOOD )then ! is good
                if( self%clsinds_spec(ispec) == CLASS_BAD ) nchanges = nchanges + 1
                self%clsinds_spec(ispec) = CLASS_GOOD
            else                                             ! is bad
                if( self%clsinds_spec(ispec) == CLASS_GOOD ) nchanges = nchanges + 1
                self%clsinds_spec(ispec) = CLASS_BAD
            endif
        end do
        !$omp end parallel do
        ! update iteration counter
        iter = iter + 1
        ! set l_converged flag
        l_converged = .false.
        if( nchanges == 0 .or. iter == MAXITS) l_converged = .true.
    end subroutine greedy_min_iter

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

    function clsind_spec_of_closest( self, ispec ) result( clsind )
        class(pspecs), intent(inout) :: self
        integer,       intent(in)    :: ispec
        integer :: closest, clsind
        closest = self%find_closest(ispec)
        clsind  = self%clsinds_spec(closest)
    end function clsind_spec_of_closest

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

    ! destructor

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
