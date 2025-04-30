module simple_pspecs
include 'simple_lib.f08'
!$ use omp_lib
!$ use omp_lib_kinds
use simple_image,  only: image
use simple_fsc,    only: plot_fsc
use simple_oris,   only: oris
use simple_masker, only: density_outside_mask
implicit none

public :: pspecs
private
#include "simple_local_flags.inc"

type pspecs
    private
    real,               allocatable :: resarr(:)                  ! resolution values in A
    real,               allocatable :: pspecs(:,:)                ! matrix of power spectra
    real,               allocatable :: pspecs_cen(:,:)            ! power spectrum centers (averages)
    real,               allocatable :: distmat(:,:)               ! distance matrix
    real,               allocatable :: dists2cens(:,:)            ! distances to center pspecs 
    real,               allocatable :: clsres(:)                  ! 2D class resolutions
    integer,            allocatable :: order(:)                   ! quality order
    integer,            allocatable :: clsinds(:)                 ! 2D class assignments
    integer,            allocatable :: clsinds_spec(:)            ! spectral class assignments, if binray: 1 is good, 2 is bad
    integer,            allocatable :: clspops(:)                 ! class populations
    integer,            allocatable :: clspops_spec(:)            ! spectral class populations
    integer,            allocatable :: ptclpops_spec(:)           ! spectral particle population
    type(stats_struct), allocatable :: clsscore_stats(:)          ! 2D class score stats
    type(stats_struct), allocatable :: cls_spec_clsscore_stats(:) ! spectral class 2D class score stats
    type(stats_struct), allocatable :: cls_spec_clsres_stats(:)   ! spectral class 2D class resolution stats
    real                 :: smpd             = 0.                 ! sampling distance
    real                 :: hp               = 0.                 ! high-pass limit
    real                 :: lp               = 0.                 ! low-pass limit
    integer              :: box              = 0                  ! box size
    integer              :: kfromto(2)       = [0,0]              ! Fourier index range
    integer              :: sz               = 0                  ! size of spectrum 
    integer              :: nspecs           = 0                  ! # of spectra
    integer              :: ncls             = 0                  ! # 2D classes
    integer              :: ncls_spec        = 0                  ! # spectral clusters
    logical              :: exists           = .false.            ! existence flag
contains
    ! constructor
    procedure            :: new
    ! getters & setters
    procedure            :: get_nspecs
    procedure            :: get_clsind_spec_state_arr
    ! plotting
    procedure            :: plot_all
    procedure            :: plot_cens
    ! clustering
    procedure, private   :: calc_distmat
    procedure, private   :: medoid_cen_init
    procedure            :: kmeans_cls_pspecs
    ! calculators
    procedure, private   :: calc_pspec_cls_avgs
    procedure, private   :: kcluster_iter
    procedure, private   :: calc_dists2cens
    procedure, private   :: cls_spec_avg_dist
    procedure, private   :: calc_cls_spec_stats
    ! destructor
    procedure            :: kill
end type pspecs

integer, parameter       :: CLASS_JUNK       = 0
integer, parameter       :: CLASS_GOOD       = 1
integer, parameter       :: CLASS_BAD        = 2
real,    parameter       :: DYNRANGE_THRES   = 1e-6
real,    parameter       :: FRAC_BEST_PTCLS  = 0.25
integer, parameter       :: MAXITS           = 10
integer, parameter       :: MINPOP           = 20

contains

    ! constructor

    subroutine new( self, ncls, imgs, os_ptcl2D, os_cls2D, msk, hp, lp, ncls_spec )
        class(pspecs),     intent(inout) :: self
        integer,           intent(in)    :: ncls
        class(image),      intent(inout) :: imgs(ncls)
        class(oris),       intent(inout) :: os_ptcl2D
        class(oris),       intent(in)    :: os_cls2D
        real,              intent(in)    :: msk, hp, lp
        integer,           intent(in)    :: ncls_spec
        real,    allocatable :: pspec(:), resarr(:)
        logical, allocatable :: mask(:)
        integer :: specinds(ncls)
        logical :: l_valid_spectra(ncls), l_junk_class
        integer :: ldim(3), icls, ispec
        real    :: dynrange
        call self%kill
        self%ncls       = ncls
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
        ! flag valid spectra
        l_valid_spectra = .false.
        !$omp parallel do default(shared) private(icls,pspec,dynrange,l_junk_class) proc_bind(close) schedule(static)
        do icls = 1, ncls
            call imgs(icls)%norm
            l_junk_class = density_outside_mask(imgs(icls), hp, msk)
            call imgs(icls)%mask(msk, 'soft')
            call imgs(icls)%spectrum('sqrt', pspec)
            dynrange = pspec(self%kfromto(1)) - pspec(self%kfromto(2))
            if( dynrange > DYNRANGE_THRES .and. os_cls2D%get_int(icls, 'pop') >= MINPOP  )then
                if( .not. l_junk_class ) l_valid_spectra(icls) = .true.
            endif
        enddo
        !$omp end parallel do
        self%nspecs = count(l_valid_spectra)
        ! allocate arrays
        allocate(self%resarr(self%sz), source=resarr(self%kfromto(1):self%kfromto(2)))
        allocate(self%pspecs(self%nspecs,self%sz), self%pspecs_cen(self%ncls_spec,self%sz),&
        &self%distmat(self%nspecs,self%nspecs),&
        &self%dists2cens(self%nspecs,self%ncls_spec), self%clsres(self%nspecs), source=0.)
        allocate(self%order(self%nspecs), self%clsinds_spec(self%nspecs), self%clsinds(self%nspecs),&
        self%clspops(self%nspecs), self%clspops_spec(self%ncls_spec), self%ptclpops_spec(self%ncls_spec), source=0)
        allocate(self%clsscore_stats(self%nspecs), self%cls_spec_clsscore_stats(self%ncls_spec),&
        &self%cls_spec_clsres_stats(self%ncls_spec))
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
                ! 2D class index
                self%clsinds(specinds(icls))   = icls
                ! power spectrum
                call imgs(icls)%spectrum('sqrt', pspec)
                self%pspecs(specinds(icls),:)  = pspec(self%kfromto(1):self%kfromto(2))
                ! 2D class population
                self%clspops(specinds(icls))   = os_cls2D%get_int(icls, 'pop')
                ! 2D class resolution
                self%clsres(specinds(icls))    = os_cls2D%get_int(icls, 'res')
                ! mask inlcuding the FRAC_BEST_PTCLS fraction of top ranking particles in 2D class
                mask = os_ptcl2D%gen_ptcl_mask('class', icls, FRAC_BEST_PTCLS)
                ! score stats for the fraction of selected particles within the class
                call os_ptcl2D%stats('corr', self%clsscore_stats(specinds(icls)), mask)
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

    subroutine calc_distmat( self )
        class(pspecs), intent(inout) :: self
        integer :: i, j
        self%distmat = 0.
        !$omp parallel do default(shared) private(i,j) proc_bind(close) schedule(dynamic)
        do i = 1, self%nspecs - 1
            do j = i + 1, self%nspecs
                self%distmat(i,j) = euclid(self%pspecs(i,:),self%pspecs(j,:))
                self%distmat(j,i) = self%distmat(i,j)
            end do
        end do
        !$omp end parallel do
    end subroutine calc_distmat

    ! make initial grouping based on dynamic range ranking
    subroutine medoid_cen_init( self )
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
    end subroutine medoid_cen_init

    ! k-means clustering of power spectra
    subroutine kmeans_cls_pspecs( self, states ) 
        class(pspecs),        intent(inout) :: self
        integer, allocatable, intent(inout) :: states(:)
        integer :: iter
        logical :: l_converged
        write(logfhandle,'(A)') 'K-MEANS CLUSTERING OF POWERSPECTRA'
        call self%medoid_cen_init
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

    function cls_spec_avg_dist( self, icls_spec ) result( dist )
        class(pspecs), intent(inout) :: self
        integer,       intent(in)    :: icls_spec
        real    :: dist
        integer :: ispec, cnt
        dist = 0.
        cnt  = 0
        do ispec = 1, self%nspecs
            if( self%clsinds_spec(ispec) == icls_spec )then
                dist = dist + self%dists2cens(ispec,self%clsinds_spec(ispec))
                cnt  = cnt + 1
            endif
        end do
        dist = dist / real(cnt)
    end function cls_spec_avg_dist

    subroutine calc_cls_spec_stats( self, l_print )
        class(pspecs), intent(inout) :: self
        logical,       intent(in)    :: l_print
        real    :: scores(self%nspecs), res(self%nspecs)
        integer :: icls_spec, cnt, ispec
        do icls_spec = 1, self%ncls_spec
            if( self%clspops_spec(icls_spec) == 0 ) cycle 
            if( l_print) write(logfhandle,'(A)') 'STATSTICS FOR SPECTRAL CLASS '//int2str(icls_spec)
            cnt = 0
            do ispec = 1,self%nspecs
                if( self%clsinds_spec(ispec) == icls_spec )then
                    cnt            = cnt + 1
                    scores(cnt)    = self%clsscore_stats(ispec)%med
                    res(cnt)       = self%clsres(ispec)
                endif
            end do
            call calc_stats(scores(:cnt), self%cls_spec_clsscore_stats(icls_spec))
            ! if( l_print) write(logfhandle,'(a,1x,f8.2)') 'MINIMUM SCORE:    ', self%cls_spec_clsscore_stats(icls_spec)%minv
            ! if( l_print) write(logfhandle,'(a,1x,f8.2)') 'MAXIMUM SCORE:    ', self%cls_spec_clsscore_stats(icls_spec)%maxv
            ! if( l_print) write(logfhandle,'(a,1x,f8.2)') 'AVERAGE SCORE:    ', self%cls_spec_clsscore_stats(icls_spec)%avg
            ! if( l_print) write(logfhandle,'(a,1x,f8.2)') 'SDEV    SCORE:    ', self%cls_spec_clsscore_stats(icls_spec)%sdev
            call calc_stats(res(:cnt), self%cls_spec_clsres_stats(icls_spec))
            ! if( l_print) write(logfhandle,'(a,1x,f8.2)') 'MINIMUM RES:      ', self%cls_spec_clsres_stats(icls_spec)%minv
            ! if( l_print) write(logfhandle,'(a,1x,f8.2)') 'MAXIMUM RES:      ', self%cls_spec_clsres_stats(icls_spec)%maxv
            ! if( l_print) write(logfhandle,'(a,1x,f8.2)') 'AVERAGE RES:      ', self%cls_spec_clsres_stats(icls_spec)%avg
            ! if( l_print) write(logfhandle,'(a,1x,f8.2)') 'SDEV    RES:      ', self%cls_spec_clsres_stats(icls_spec)%sdev
            if( l_print) write(logfhandle,'(a,1x,f10.6)') 'SPEC CLS DIST: ', self%cls_spec_avg_dist(icls_spec)
            if( l_print) write(logfhandle,'(a,1x,f10.2)') 'MEDIAN  SCORE: ', self%cls_spec_clsscore_stats(icls_spec)%med
            if( l_print) write(logfhandle,'(a,1x,f10.2)') 'MEDIAN    RES: ', self%cls_spec_clsres_stats(icls_spec)%med
        end do
    end subroutine calc_cls_spec_stats

    ! destructor

    subroutine kill( self )
        class(pspecs), intent(inout) :: self
        if( self%exists )then
            if( allocated(self%resarr)                  ) deallocate(self%resarr)
            if( allocated(self%pspecs)                  ) deallocate(self%pspecs)
            if( allocated(self%pspecs_cen)              ) deallocate(self%pspecs_cen)
            if( allocated(self%distmat)                 ) deallocate(self%distmat)
            if( allocated(self%dists2cens)              ) deallocate(self%dists2cens)
            if( allocated(self%order)                   ) deallocate(self%order)
            if( allocated(self%clsinds_spec)            ) deallocate(self%clsinds_spec)
            if( allocated(self%clsinds)                 ) deallocate(self%clsinds)
            if( allocated(self%clspops)                 ) deallocate(self%clspops)
            if( allocated(self%clsres)                  ) deallocate(self%clsres)
            if( allocated(self%clspops_spec)            ) deallocate(self%clspops_spec)
            if( allocated(self%ptclpops_spec)           ) deallocate(self%ptclpops_spec)  
            if( allocated(self%clsscore_stats)          ) deallocate(self%clsscore_stats)
            if( allocated(self%cls_spec_clsscore_stats) ) deallocate(self%cls_spec_clsscore_stats)
            if( allocated(self%cls_spec_clsres_stats)   ) deallocate(self%cls_spec_clsres_stats)
            self%exists = .false.
        endif
    end subroutine kill

end module simple_pspecs
