module simple_pspecs
include 'simple_lib.f08'
!$ use omp_lib
!$ use omp_lib_kinds
use simple_image,    only: image
use simple_fsc,      only: plot_fsc
use simple_oris,     only: oris
use simple_aff_prop, only: aff_prop
use simple_kmedoids, only: kmedoids
implicit none

public :: pspecs
private
#include "simple_local_flags.inc"

type pspecs
    private
    real,               allocatable :: resarr(:)                ! resolution values in A
    real,               allocatable :: pspecs(:,:)              ! matrix of power spectra
    real,               allocatable :: pspecs_cen(:,:)          ! power spectrum centers (averages)
    real,               allocatable :: distmat(:,:)             ! distance matrix
    real,               allocatable :: dists2cens(:,:)          ! distances to center pspecs 
    real,               allocatable :: clsres(:)                ! 2D class resolutions
    integer,            allocatable :: order(:)                 ! quality order
    integer,            allocatable :: clsinds(:)               ! 2D class assignments
    integer,            allocatable :: clsinds_spec(:)          ! spectral class assignments, if binray: 1 is good, 2 is bad
    integer,            allocatable :: clspops(:)               ! class populations
    integer,            allocatable :: clspops_spec(:)          ! spectral class populations
    type(stats_struct), allocatable :: cls_spec_clsres_stats(:) ! spectral class 2D class resolution stats
    real               :: smpd       = 0.                       ! sampling distance
    real               :: hp         = 0.                       ! high-pass limit
    real               :: lp         = 0.                       ! low-pass limit
    integer            :: box        = 0                        ! box size
    integer            :: kfromto(2) = [0,0]                    ! Fourier index range
    integer            :: sz         = 0                        ! size of spectrum 
    integer            :: nspecs     = 0                        ! # of spectra
    integer            :: ncls       = 0                        ! # 2D classes
    integer            :: ncls_spec  = 0                        ! # spectral clusters
    logical            :: exists     = .false.                  ! existence flag
contains
    ! constructor
    procedure          :: new
    procedure          :: update_ncls_spec
    ! getters & setters
    procedure          :: get_nspecs
    procedure          :: get_clsind_spec_state_arr
    procedure          :: get_distmat
    ! plotting
    procedure          :: plot_all
    procedure          :: plot_cens
    ! clustering
    procedure          :: calc_distmat
    procedure, private :: medoid_cls_init
    procedure          :: kmedoid_cls_pspecs
    procedure          :: aprop_cls_pspecs
    procedure          :: kmeans_cls_pspecs
    ! calculators
    procedure, private :: calc_pspec_cls_avgs
    procedure, private :: kcluster_iter
    procedure, private :: calc_dists2cens
    procedure, private :: calc_cls_spec_stats
    procedure          :: rank_cls_spec
    procedure          :: select_cls_spec
    ! destructor
    procedure          :: kill
end type pspecs

integer, parameter :: MAXITS = 10

contains

    ! constructor

    subroutine new( self, imgs, os_cls2D, msk, hp, lp, ncls_spec, l_exclude_junk )
        use simple_strategy2D_utils, only: flag_non_junk_cavgs
        class(pspecs),     intent(inout) :: self
        class(image),      intent(inout) :: imgs(:)
        class(oris),       intent(in)    :: os_cls2D
        real,              intent(in)    :: msk, hp, lp
        integer,           intent(in)    :: ncls_spec
        logical, optional, intent(in)    :: l_exclude_junk
        type(image), allocatable :: cavg_threads(:) 
        real,        allocatable :: pspec(:), resarr(:)
        logical,     allocatable :: l_valid_spectra(:)
        integer :: specinds(size(imgs))
        logical :: l_junk_class, ll_exclude_junk
        integer :: ldim(3), icls, ispec, ithr
        real    :: dynrange
        call self%kill
        ll_exclude_junk = .true.
        if( present(l_exclude_junk) ) ll_exclude_junk = l_exclude_junk
        self%ncls       = size(imgs)
        ldim            = imgs(1)%get_ldim()
        self%box        = ldim(1)
        self%smpd       = imgs(1)%get_smpd()
        self%hp         = hp
        self%lp         = lp
        self%ncls_spec  = ncls_spec
        resarr          = get_resarr(self%box, self%smpd)
        self%kfromto(1) = calc_fourier_index(self%hp, self%box, self%smpd)
        self%kfromto(2) = calc_fourier_index(self%lp, self%box, self%smpd)
        self%sz         = self%kfromto(2) - self%kfromto(1) + 1
        if( ll_exclude_junk )then
            call flag_non_junk_cavgs(imgs, hp, msk, l_valid_spectra, os_cls2D)
        else
            allocate(l_valid_spectra(self%ncls), source=.true.)
        endif
        self%nspecs     = count(l_valid_spectra)
        ! allocate arrays
        allocate(self%resarr(self%sz), source=resarr(self%kfromto(1):self%kfromto(2)))
        allocate(self%pspecs(self%nspecs,self%sz), self%pspecs_cen(self%ncls_spec,self%sz),&
        &self%distmat(self%nspecs,self%nspecs),&
        &self%dists2cens(self%nspecs,self%ncls_spec), self%clsres(self%nspecs), source=0.)
        allocate(self%order(self%nspecs), self%clsinds_spec(self%nspecs), self%clsinds(self%nspecs),&
        self%clspops(self%nspecs), self%clspops_spec(self%ncls_spec), source=0)
        allocate(self%cls_spec_clsres_stats(self%ncls_spec))
        ! set spectrum indices
        ispec = 0
        do icls = 1, self%ncls
            if( l_valid_spectra(icls) )then
                ispec = ispec + 1
                specinds(icls) = ispec
            endif
        end do
        ! fill-up the instance
        allocate(cavg_threads(nthr_glob))
        do ithr = 1, nthr_glob
            call cavg_threads(ithr)%new(ldim, self%smpd)
        end do
        !$omp parallel do default(shared) private(icls,ithr,pspec) proc_bind(close) schedule(static)
        do icls = 1, self%ncls
            if( l_valid_spectra(icls) )then
                ithr = omp_get_thread_num() + 1
                ! 2D class index
                self%clsinds(specinds(icls))   = icls
                ! power spectrum
                call cavg_threads(ithr)%copy(imgs(icls))
                call cavg_threads(ithr)%norm
                call cavg_threads(ithr)%mask(msk, 'soft')
                call cavg_threads(ithr)%spectrum('sqrt', pspec)
                self%pspecs(specinds(icls),:) = pspec(self%kfromto(1):self%kfromto(2))
                ! 2D class population
                self%clspops(specinds(icls))  = os_cls2D%get_int(icls, 'pop')
                ! 2D class resolution
                self%clsres(specinds(icls))   = os_cls2D%get_int(icls, 'res')
                deallocate(pspec)
            endif
        end do
        !$omp end parallel do
        do ithr = 1, nthr_glob
            call cavg_threads(ithr)%kill
        end do
        deallocate(cavg_threads, resarr)
        self%exists = .true.
    end subroutine new

    subroutine update_ncls_spec( self, ncls_spec )
        class(pspecs),     intent(inout) :: self
        integer,           intent(in)    :: ncls_spec
        self%ncls_spec  = ncls_spec
        ! deallocate arrays
        if( allocated(self%pspecs_cen)            ) deallocate(self%pspecs_cen)
        if( allocated(self%dists2cens)            ) deallocate(self%dists2cens)
        if( allocated(self%clspops_spec)          ) deallocate(self%clspops_spec)
        if( allocated(self%cls_spec_clsres_stats) ) deallocate(self%cls_spec_clsres_stats)
        ! allocate arrays
        allocate(self%pspecs_cen(self%ncls_spec,self%sz), self%dists2cens(self%nspecs,self%ncls_spec), source=0.)
        allocate(self%clspops_spec(self%ncls_spec), source=0)
        allocate(self%cls_spec_clsres_stats(self%ncls_spec))
    end subroutine update_ncls_spec

    ! getters & setters

    function get_nspecs( self ) result( nspecs )
        class(pspecs), intent(in) :: self
        integer :: nspecs
        nspecs = self%nspecs
    end function get_nspecs

    function get_clsind_spec_state_arr( self ) result( states )
        class(pspecs), intent(in) :: self
        integer, allocatable :: states(:)
        integer :: icls_spec, ispec
        allocate(states(self%ncls), source=0)
        do icls_spec = 1, self%ncls_spec
            do ispec = 1, self%nspecs
                if( self%clsinds_spec(ispec) == icls_spec ) states(self%clsinds(ispec)) = icls_spec
            end do
        enddo
    end function get_clsind_spec_state_arr

    function get_distmat( self ) result( distmat )
        class(pspecs), intent(in) :: self
        real, allocatable :: distmat(:,:)
        allocate(distmat(self%nspecs,self%nspecs), source=self%distmat)
    end function get_distmat

    ! plotting

    subroutine plot_cens( self, fbody )
        class(pspecs),    intent(in)  :: self
        character(len=*), intent(in)  :: fbody
        character(len=:), allocatable :: ffbody
        integer :: icen
        do icen = 1, self%ncls_spec
            ffbody = trim(fbody)//'_cen'//int2str_pad(icen,2)
            call plot_fsc(self%sz, self%pspecs_cen(icen,:), self%resarr, self%smpd, ffbody)
        enddo
    end subroutine plot_cens

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

    ! clustering

    subroutine calc_distmat( self, is_l1 )
        class(pspecs),     intent(inout) :: self
        logical, optional, intent(in)    :: is_l1
        logical :: ll1
        integer :: i, j
        ll1          = .false.
        self%distmat = 0.
        if( present(is_l1) ) ll1 = is_l1
        if( ll1 )then
            !$omp parallel do default(shared) private(i,j) proc_bind(close) schedule(dynamic)
            do i = 1, self%nspecs - 1
                do j = i + 1, self%nspecs
                    self%distmat(i,j) = l1dist(self%pspecs(i,:),self%pspecs(j,:))
                    self%distmat(j,i) = self%distmat(i,j)
                end do
            end do
            !$omp end parallel do
        else
            !$omp parallel do default(shared) private(i,j) proc_bind(close) schedule(dynamic)
            do i = 1, self%nspecs - 1
                do j = i + 1, self%nspecs
                    self%distmat(i,j) = euclid(self%pspecs(i,:),self%pspecs(j,:))
                    self%distmat(j,i) = self%distmat(i,j)
                end do
            end do
            !$omp end parallel do
        endif
    end subroutine calc_distmat

    ! make initial grouping based on dynamic range ranking
    subroutine medoid_cls_init( self )
        class(pspecs), intent(inout) :: self
        integer, allocatable :: parts(:,:)
        real    :: medoid_dists(self%nspecs)
        integer :: ispec, i_medoid, icls_spec
        ! calculate distance matrix
        call self%calc_distmat
        ! identify medoid
        call medoid_from_dmat(self%distmat, i_medoid)
        ! order spectra according to similarity to medoid
        !$omp parallel do default(shared) private(ispec) proc_bind(close)
        do ispec = 1, self%nspecs
            medoid_dists(ispec) = euclid(self%pspecs(ispec,:), self%pspecs(i_medoid,:))
        end do
        !$omp end parallel do
        self%order = (/(ispec,ispec=1,self%nspecs)/)
        call hpsort(medoid_dists, self%order)
        parts = split_nobjs_even(self%nspecs, self%ncls_spec)
        ! assign classes based on similarity to medoid
        do icls_spec = 1, self%ncls_spec
            do ispec = parts(icls_spec,1), parts(icls_spec,2)
                self%clsinds_spec(self%order(ispec)) = icls_spec
            end do
        end do
        deallocate(parts)
    end subroutine medoid_cls_init

    ! k-medoids clustering of powerspectra
    subroutine kmedoid_cls_pspecs( self, states )
        class(pspecs),        intent(inout) :: self
        integer, allocatable, intent(inout) :: states(:)
        type(kmedoids) :: kmed
        call self%calc_distmat
        call kmed%new(self%nspecs, self%distmat, self%ncls_spec)
        call kmed%init
        call kmed%cluster
        call kmed%get_labels(self%clsinds_spec)
        call kmed%kill
        call self%calc_pspec_cls_avgs
        call self%calc_cls_spec_stats(l_print=.true.)
        if( allocated(states) ) deallocate(states)
        states = self%get_clsind_spec_state_arr()
    end subroutine kmedoid_cls_pspecs

    ! affinity propagation clustering of powerspectra
    subroutine aprop_cls_pspecs( self, states )
        class(pspecs),        intent(inout) :: self
        integer, allocatable, intent(inout) :: states(:)
        real ,   allocatable :: smat(:,:)
        integer, allocatable :: centers(:), labels(:)
        type(aff_prop) :: aprop
        real    :: pref, simsum
        integer :: ncls_aff_prop
        ! calculate distance matrix
        call self%calc_distmat
        ! convert to similarity matrix
        smat = dmat2smat(self%distmat)
        ! calculate a preference that generates a small number of clusters
        pref = calc_ap_pref(smat, 'min_minus_med')
        call aprop%new(self%nspecs, smat, pref=pref)
        call aprop%propagate(centers, labels, simsum)
        call aprop%kill
        ncls_aff_prop = size(centers)
        write(logfhandle,'(A,I3)') '>>> # CLUSTERS FOUND BY AFFINITY PROPAGATION (AP): ', ncls_aff_prop
        call self%update_ncls_spec(ncls_aff_prop)
        self%clsinds_spec = labels
        call self%calc_pspec_cls_avgs
        call self%calc_cls_spec_stats(l_print=.true.)
        if( allocated(states) ) deallocate(states)
        states = self%get_clsind_spec_state_arr()
    end subroutine aprop_cls_pspecs

    ! k-means clustering of power spectra
    subroutine kmeans_cls_pspecs( self, states ) 
        class(pspecs),        intent(inout) :: self
        integer, allocatable, intent(inout) :: states(:)
        integer :: iter
        logical :: l_converged
        write(logfhandle,'(A)') 'K-MEANS CLUSTERING OF POWERSPECTRA'
        call self%medoid_cls_init
        call self%calc_pspec_cls_avgs
        iter = 0
        l_converged = .false.
        do
            call self%kcluster_iter(iter, l_converged)
            if( l_converged ) exit
        end do
        call self%calc_cls_spec_stats(l_print=.true.)
        if( allocated(states) ) deallocate(states)
        states = self%get_clsind_spec_state_arr()
    end subroutine kmeans_cls_pspecs

    ! calculators

    subroutine calc_pspec_cls_avgs( self )
        class(pspecs), intent(inout) :: self
        integer :: icls_spec, ispec
        ! accumulate class sums and populations
        !$omp parallel default(shared) private(icls_spec,ispec) proc_bind(close)
        !$omp do schedule(static)
        do icls_spec = 1, self%ncls_spec
            self%pspecs_cen(icls_spec,:) = 0.
            self%clspops_spec(icls_spec) = 0
            do ispec = 1, self%nspecs
                if( self%clsinds_spec(ispec) == icls_spec )then
                    self%pspecs_cen(icls_spec,:) = self%pspecs_cen(icls_spec,:) + self%pspecs(ispec,:)
                    self%clspops_spec(icls_spec) = self%clspops_spec(icls_spec) + 1
                endif
            end do
        end do
        !$omp end do
        ! calculate class averages
        !$omp do schedule(static)
        do icls_spec = 1, self%ncls_spec
            if( self%clspops_spec(icls_spec) > 0 )then
                self%pspecs_cen(icls_spec,:) = self%pspecs_cen(icls_spec,:) / real(self%clspops_spec(icls_spec))
            else
                self%pspecs_cen(icls_spec,:) = 0.
            endif
        enddo
        !$omp end do
        !$omp end parallel 
    end subroutine calc_pspec_cls_avgs

    subroutine kcluster_iter( self, iter, l_converged )
        class(pspecs), intent(inout) :: self
        integer,       intent(inout) :: iter
        logical,       intent(inout) :: l_converged
        integer :: nchanges, ispec, loc(self%nspecs)
        real    :: dist
        ! calculate distances
        call self%calc_dists2cens
        ! assign clusters
        loc = minloc(self%dists2cens, dim=2)
        nchanges = count(self%clsinds_spec /= loc)
        self%clsinds_spec = loc
        ! find cluster centers
        call self%calc_pspec_cls_avgs
        ! update iteration counter
        iter = iter + 1
        ! report joint distance
        dist = 0.
        do ispec = 1, self%nspecs
            dist = dist + self%dists2cens(ispec,self%clsinds_spec(ispec))
        end do
        write(logfhandle,'(a,1x,f8.2)') 'ITER: '//int2str(iter)//' DIST: ', dist
        ! set l_converged flag
        l_converged = .false.
        if( nchanges == 0 .or. iter == MAXITS) l_converged = .true.
    end subroutine kcluster_iter

    subroutine calc_dists2cens( self )
        class(pspecs), intent(inout) :: self
        integer :: ispec, icls_spec
        !$omp parallel do default(shared) private(ispec,icls_spec) proc_bind(close)
        do ispec = 1, self%nspecs
            do icls_spec = 1, self%ncls_spec
                self%dists2cens(ispec,icls_spec) = euclid(self%pspecs(ispec,:), self%pspecs_cen(icls_spec,:))
            end do
        end do
        !$omp end parallel do 
    end subroutine calc_dists2cens

    subroutine calc_cls_spec_stats( self, l_print )
        class(pspecs), intent(inout) :: self
        logical,       intent(in)    :: l_print
        real    :: res(self%nspecs)
        integer :: icls_spec, cnt, ispec
        do icls_spec = 1, self%ncls_spec
            if( self%clspops_spec(icls_spec) == 0 ) cycle 
            if( l_print) write(logfhandle,'(A)') 'RESOLUTION STATSTICS FOR SPECTRAL CLASS '//int2str(icls_spec)
            cnt = 0
            do ispec = 1,self%nspecs
                if( self%clsinds_spec(ispec) == icls_spec )then
                    cnt            = cnt + 1
                    res(cnt)       = self%clsres(ispec)
                endif
            end do
            call calc_stats(res(:cnt), self%cls_spec_clsres_stats(icls_spec))
            if( l_print) write(logfhandle,'(a,1x,f8.2)') 'MINIMUM RES: ', self%cls_spec_clsres_stats(icls_spec)%minv
            if( l_print) write(logfhandle,'(a,1x,f8.2)') 'MAXIMUM RES: ', self%cls_spec_clsres_stats(icls_spec)%maxv
            if( l_print) write(logfhandle,'(a,1x,f8.2)') 'AVERAGE RES: ', self%cls_spec_clsres_stats(icls_spec)%avg
            if( l_print) write(logfhandle,'(a,1x,f8.2)') 'MEDIAN  RES: ', self%cls_spec_clsres_stats(icls_spec)%med
            if( l_print) write(logfhandle,'(a,1x,f8.2)') 'SDEV    RES: ', self%cls_spec_clsres_stats(icls_spec)%sdev
        end do
    end subroutine calc_cls_spec_stats

    subroutine rank_cls_spec( self, states )
        class(pspecs),        intent(inout) :: self
        integer, allocatable, intent(inout) :: states(:)
        real,    allocatable :: tmp(:)
        integer :: cls_spec_order(self%ncls_spec), icls_spec, rank, rank_assign(self%nspecs), ispec
        tmp = self%cls_spec_clsres_stats(:)%med
        cls_spec_order = (/(icls_spec,icls_spec=1,self%ncls_spec)/)
        call hpsort(tmp, cls_spec_order)
        do rank = 1, self%ncls_spec
            do ispec = 1, self%nspecs
                if( self%clsinds_spec(ispec) == cls_spec_order(rank) ) rank_assign(ispec) = rank
            end do
        end do
        self%clsinds_spec = rank_assign
        call self%calc_cls_spec_stats(l_print=.true.)
        if( allocated(states) ) deallocate(states)
        states = self%get_clsind_spec_state_arr()
        deallocate(tmp)
    end subroutine rank_cls_spec

    subroutine select_cls_spec( self, states )
        class(pspecs),        intent(inout) :: self
        integer, allocatable, intent(inout) :: states(:)
        real,    allocatable :: res_good(:), res_bad(:)
        integer :: good_bad_assign(self%nspecs), icls_spec, ispec, nptcls, nptcls_good
        real    :: best_res, worst_res, dist2best, dist2worst, frac_good
        type(stats_struct) :: res_stats
        ! assign good/bad 2D classes based on closesness to best and worst median resolution of spectral groups
        best_res  = minval(self%cls_spec_clsres_stats(:)%med)
        worst_res = maxval(self%cls_spec_clsres_stats(:)%med)
        good_bad_assign = 0
        do icls_spec = 1, self%ncls_spec
            dist2best  = abs(self%cls_spec_clsres_stats(icls_spec)%med - best_res)
            dist2worst = abs(self%cls_spec_clsres_stats(icls_spec)%med - worst_res)
            if( dist2best < dist2worst )then
                where(self%clsinds_spec == icls_spec) good_bad_assign = 1
            endif
        end do
        nptcls      = sum(self%clspops)
        nptcls_good = sum(self%clspops, mask=good_bad_assign == 1 )
        frac_good   = real(nptcls_good) / real(nptcls)
        write(logfhandle,'(a,1x,f8.2)') '% PARTICLES CLASSIFIED AS GOOD: ', frac_good * 100.
        ! calculate resolution statistics for good/bad classes
        res_good    = pack(self%clsres, mask=good_bad_assign == 1)
        res_bad     = pack(self%clsres, mask=good_bad_assign == 0)
        call calc_stats(res_good, res_stats)
        write(logfhandle,'(A)') 'RESOLUTION STATSTICS FOR GOOD PARTITION'
        write(logfhandle,'(a,1x,f8.2)') 'MINIMUM RES: ', res_stats%minv
        write(logfhandle,'(a,1x,f8.2)') 'MAXIMUM RES: ', res_stats%maxv
        write(logfhandle,'(a,1x,f8.2)') 'AVERAGE RES: ', res_stats%avg
        write(logfhandle,'(a,1x,f8.2)') 'MEDIAN  RES: ', res_stats%med
        write(logfhandle,'(a,1x,f8.2)') 'SDEV    RES: ', res_stats%sdev
        call calc_stats(res_bad, res_stats)
        write(logfhandle,'(A)') 'RESOLUTION STATSTICS FOR BAD  PARTITION'
        write(logfhandle,'(a,1x,f8.2)') 'MINIMUM RES: ', res_stats%minv
        write(logfhandle,'(a,1x,f8.2)') 'MAXIMUM RES: ', res_stats%maxv
        write(logfhandle,'(a,1x,f8.2)') 'AVERAGE RES: ', res_stats%avg
        write(logfhandle,'(a,1x,f8.2)') 'MEDIAN  RES: ', res_stats%med
        write(logfhandle,'(a,1x,f8.2)') 'SDEV    RES: ', res_stats%sdev
        ! translate good/bad assignment to state array output
        if( allocated(states) ) deallocate(states)
        allocate(states(self%ncls), source=0)
        do ispec = 1, self%nspecs
            states(self%clsinds(ispec)) = good_bad_assign(ispec)
        end do
    end subroutine select_cls_spec

    ! destructor

    subroutine kill( self )
        class(pspecs), intent(inout) :: self
        if( self%exists )then
            if( allocated(self%resarr)                ) deallocate(self%resarr)
            if( allocated(self%pspecs)                ) deallocate(self%pspecs)
            if( allocated(self%pspecs_cen)            ) deallocate(self%pspecs_cen)
            if( allocated(self%distmat)               ) deallocate(self%distmat)
            if( allocated(self%dists2cens)            ) deallocate(self%dists2cens)
            if( allocated(self%order)                 ) deallocate(self%order)
            if( allocated(self%clsinds_spec)          ) deallocate(self%clsinds_spec)
            if( allocated(self%clsinds)               ) deallocate(self%clsinds)
            if( allocated(self%clspops)               ) deallocate(self%clspops)
            if( allocated(self%clsres)                ) deallocate(self%clsres)
            if( allocated(self%clspops_spec)          ) deallocate(self%clspops_spec)  
            if( allocated(self%cls_spec_clsres_stats) ) deallocate(self%cls_spec_clsres_stats)
            self%exists = .false.
        endif
    end subroutine kill

end module simple_pspecs
