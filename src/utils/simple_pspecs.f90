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
    real,               allocatable :: pspecs(:,:)                ! matrix of power spectra
    real,               allocatable :: pspecs_cen(:,:)            ! power spectrum centers (averages)
    real,               allocatable :: pspec_junk(:)              ! 'junk' average power spectrum
    real,               allocatable :: pspec_good(:)              ! 'good' average power spectrum
    real,               allocatable :: pspec_bad(:)               ! 'bad'  average power spectrum
    real,               allocatable :: resarr(:)                  ! resolution values in A
    real,               allocatable :: dynranges(:)               ! dynamic ranges of power spectra
    real,               allocatable :: dists2cens(:,:)            ! distances to center pspecs 
    real,               allocatable :: dists2good(:)              ! distances to good pspecs
    real,               allocatable :: dists2bad(:)               ! distances to bad  pspecs
    real,               allocatable :: distmat(:,:)               ! Euclidean distance matrix
    real,               allocatable :: clsres(:)                  ! 2D class resolutions
    integer,            allocatable :: ranks(:)                   ! quality ranking
    integer,            allocatable :: order(:)                   ! quality order
    integer,            allocatable :: clsinds_spec(:)            ! spectral class assignments, if binray: 1 is good, 2 is bad
    integer,            allocatable :: clsinds(:)                 ! 2D class assignments
    integer,            allocatable :: clspops(:)                 ! class populations
    integer,            allocatable :: clspops_spec(:)            ! spectral class populations
    integer,            allocatable :: ptclpops_spec(:)           ! spectral particle population
    logical,            allocatable :: l_junk_spec(:)             ! junk spectrum assignment
    type(stats_struct), allocatable :: clsscore_stats(:)          ! 2D class score stats
    type(stats_struct), allocatable :: cls_spec_clsscore_stats(:) ! spectral class 2D class score stats
    type(stats_struct), allocatable :: cls_spec_clsres_stats(:)   ! spectral class 2D class resolution stats
    real                 :: smpd             = 0.                 ! sampling distance
    real                 :: hp               = 0.                 ! high-pass limit
    real                 :: lp               = 0.                 ! low-pass limit
    real                 :: dynrange_t       = 0.                 ! dynamic range threshold
    real                 :: cls_spec_score_t = 0.                 ! spectral class score
    integer              :: box              = 0                  ! box size
    integer              :: kfromto(2)                            ! Fourier index range
    integer              :: sz               = 0                  ! size of spectrum 
    integer              :: njunk            = 0                  ! # of junk spectra
    integer              :: nspecs           = 0                  ! # of spectra
    integer              :: ncls             = 0                  ! # 2D classes
    integer              :: ncls_spec        = 0                  ! # spectral clusters
    integer              :: ngood            = 0                  ! # good power spectra
    integer              :: nbad             = 0                  ! # bad power spectra
    logical              :: exists           = .false.            ! existence flag
contains
    ! constructor
    procedure            :: new
    ! getters & setters
    procedure            :: get_nspecs
    procedure            :: get_frac_good
    procedure            :: get_ngood
    procedure            :: set_ngood
    procedure            :: get_good_bad_state_arr
    procedure            :: get_clsind_spec_state_arr
    procedure            :: get_ordered_clsind
    ! plotting
    procedure            :: plot_all
    procedure            :: plot_good_bad
    procedure            :: plot_cens
    procedure            :: plot_all_ranked
    ! clustering 
    procedure            :: dynrange_cen_init
    procedure            :: kmeans_cls_pspecs_and_rank
    ! calculators
    procedure            :: median_good_clspop
    procedure, private   :: calc_pspec_cls_avgs
    procedure, private   :: kcluster_iter
    procedure, private   :: calc_dists2cens
    procedure, private   :: rank_pspecs
    procedure, private   :: calc_cls_spec_stats
    procedure, private   :: rank_spectral_classes
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
        integer, optional, intent(in)    :: ncls_spec
        real,    allocatable :: pspec(:), pspec_junk(:), resarr(:)
        logical, allocatable :: mask(:)
        integer :: specinds(ncls)
        logical :: l_valid_spectra(ncls), l_junk_class(ncls)
        integer :: ldim(3), icls, ispec
        real    :: dynrange
        call self%kill
        self%ncls_spec  = 2 ! binary clustering by default
        if( present(ncls_spec) ) self%ncls_spec = ncls_spec
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
        ! count junk classes
        l_junk_class    = .false.
        ! create junk power spectrum
        allocate(pspec_junk(self%sz), source=0.)
        !$omp parallel do default(shared) private(icls,pspec,dynrange) reduction(+:pspec_junk) proc_bind(close) schedule(static)
        do icls = 1, ncls
            call imgs(icls)%norm
            if( density_outside_mask(imgs(icls), hp, msk) ) l_junk_class(icls) = .true.
            call imgs(icls)%mask(msk, 'soft')
            call imgs(icls)%spectrum('sqrt', pspec)
            dynrange = pspec(self%kfromto(1)) - pspec(self%kfromto(2))
            if( dynrange > DYNRANGE_THRES .and. l_junk_class(icls) )then
                pspec_junk = pspec_junk + pspec(self%kfromto(1):self%kfromto(2))
            endif
            if( dynrange > DYNRANGE_THRES .and. os_cls2D%get_int(icls, 'pop') >= MINPOP  )then
                if( .not. l_junk_class(icls) ) l_valid_spectra(icls) = .true.
            endif
        enddo
        !$omp end parallel do
        self%njunk = count(l_junk_class)
        if( self%njunk > 0 )then
            allocate(self%pspec_junk(self%sz), source=pspec_junk/real(self%njunk))
            deallocate(pspec_junk)
        endif
        self%nspecs = count(l_valid_spectra)
        ! allocate arrays
        allocate(self%resarr(self%sz), source=resarr(self%kfromto(1):self%kfromto(2)))
        allocate(self%pspecs(self%nspecs,self%sz), self%pspecs_cen(self%ncls_spec,self%sz), self%dynranges(self%nspecs),&
        &self%dists2cens(self%nspecs,self%ncls_spec), self%dists2good(self%nspecs), self%dists2bad(self%nspecs),&
        &self%clsres(self%nspecs), source=0.)
        allocate(self%ranks(self%nspecs), self%order(self%nspecs), self%clsinds_spec(self%nspecs), self%clsinds(self%nspecs),&
        self%clspops(self%nspecs), self%clspops_spec(self%ncls_spec), self%ptclpops_spec(self%ncls_spec), source=0)
        allocate(self%l_junk_spec(self%nspecs), source=.false.)
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
                self%clsinds(specinds(icls))      = icls
                ! power spectrum
                call imgs(icls)%spectrum('sqrt', pspec)
                self%pspecs(specinds(icls),:)     = pspec(self%kfromto(1):self%kfromto(2))
                ! spectral dynamic range
                self%dynranges(specinds(icls))    = self%pspecs(specinds(icls),1) - self%pspecs(specinds(icls),self%sz)
                ! 2D class population
                self%clspops(specinds(icls))      = os_cls2D%get_int(icls, 'pop')
                ! 2D class resolution
                self%clsres(specinds(icls))       = os_cls2D%get_int(icls, 'res')
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

    ! make initial grouping based on dynamic range ranking
    subroutine dynrange_cen_init( self )
        class(pspecs), intent(inout) :: self
        real, allocatable :: tmp(:)
        integer           :: ispec
        real              :: dist
        real              :: transl_tab(self%ncls_spec)
        tmp = pack(self%dynranges, mask=self%dynranges > DYNRANGE_THRES)
        if( self%ncls_spec == 2 )then
            call otsu(size(tmp), tmp, self%dynrange_t)
            do ispec = 1, self%nspecs
                if( self%dynranges(ispec) >= self%dynrange_t )then
                    self%clsinds_spec(ispec) = CLASS_GOOD
                else
                    self%clsinds_spec(ispec) = CLASS_BAD
                endif
            enddo
        else
            call quantize_vec_serial(tmp, self%ncls_spec, transl_tab)
            ! assign clusters
            do ispec = 1, self%nspecs
                call find(transl_tab, self%ncls_spec, self%dynranges(ispec), self%clsinds_spec(ispec), dist)
            end do
        endif
        deallocate(tmp)
    end subroutine dynrange_cen_init 

    ! k-means clustering of power spectra
    subroutine kmeans_cls_pspecs_and_rank( self, states ) 
        class(pspecs),        intent(inout) :: self
        integer, allocatable, intent(inout) :: states(:)
        integer :: iter
        logical :: l_converged
        write(logfhandle,'(A)') 'K-MEANS CLUSTERING OF POWERSPECTRA'
        call self%dynrange_cen_init
        call self%calc_pspec_cls_avgs
        if( self%ncls_spec == 2 ) call self%plot_good_bad('pspec_good_dynrange', 'pspec_bad_dynrange')
        iter = 0
        l_converged = .false.
        do
            print *, 'njunk: ', count(self%l_junk_spec)
            call self%kcluster_iter(iter, l_converged)
            if( l_converged ) exit
        end do
        if( self%ncls_spec == 2 )then
            call self%plot_good_bad('pspec_good_kmeans', 'pspec_bad_kmeans')
            call self%rank_pspecs
        else
            call self%calc_cls_spec_stats(l_print=.false.)
            call self%rank_spectral_classes
            call self%calc_cls_spec_stats(l_print=.true.)
        endif
        if( allocated(states) ) deallocate(states)
        states = self%get_clsind_spec_state_arr()
    end subroutine kmeans_cls_pspecs_and_rank

    ! calculators

    function median_good_clspop( self ) result( med )
        class(pspecs), intent(in) :: self
        real, allocatable :: clspops_good(:)
        real :: med
        clspops_good = real(pack(self%clspops, mask=self%clsinds_spec == CLASS_GOOD))
        med = median(clspops_good)
    end function median_good_clspop

    subroutine calc_pspec_cls_avgs( self )
        class(pspecs), intent(inout) :: self
        real, allocatable :: pspec_good(:), pspec_bad(:)
        integer :: icls_spec, ispec
        if( self%ncls_spec == 2 )then
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
        else
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
        endif
    end subroutine calc_pspec_cls_avgs

    subroutine kcluster_iter( self, iter, l_converged )
        class(pspecs), intent(inout) :: self
        integer,       intent(inout) :: iter
        logical,       intent(inout) :: l_converged
        integer :: nchanges, ispec, loc(self%nspecs)
        if( self%ncls_spec == 2 )then
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
        else
            ! calculate distances
            call self%calc_dists2cens
            ! assign clusters
            loc = minloc(self%dists2cens, dim=2)
            where(self%l_junk_spec) loc = 0
            nchanges = count(self%clsinds_spec /= loc)
            self%clsinds_spec = loc
        endif
        ! find cluster centers
        call self%calc_pspec_cls_avgs
        ! update iteration counter
        iter = iter + 1
        ! set l_converged flag
        l_converged = .false.
        if( nchanges == 0 .or. iter == MAXITS) l_converged = .true.
    end subroutine kcluster_iter

    subroutine calc_dists2cens( self )
        class(pspecs), intent(inout) :: self
        integer :: ispec, icls_spec
        real    :: dist2junk
        !$omp parallel do default(shared) private(ispec,icls_spec,dist2junk) proc_bind(close)
        do ispec = 1, self%nspecs
            do icls_spec = 1, self%ncls_spec
                self%dists2cens(ispec,icls_spec) = euclid(self%pspecs(ispec,:), self%pspecs_cen(icls_spec,:))
            end do
            if( self%njunk > 0 )then
                dist2junk = euclid(self%pspecs(ispec,:), self%pspec_junk)
                if( dist2junk <= minval(self%dists2cens(ispec,:)) )then
                    self%l_junk_spec(ispec) = .true.
                else
                    self%l_junk_spec(ispec) = .false.
                endif
            endif
        end do
        !$omp end parallel do 
    end subroutine calc_dists2cens

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
        self%clsinds_spec = self%ranks
    end subroutine rank_pspecs

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
            if( l_print) write(logfhandle,'(a,1x,f8.2)') 'MINIMUM SCORE:    ', self%cls_spec_clsscore_stats(icls_spec)%minv
            if( l_print) write(logfhandle,'(a,1x,f8.2)') 'MAXIMUM SCORE:    ', self%cls_spec_clsscore_stats(icls_spec)%maxv
            if( l_print) write(logfhandle,'(a,1x,f8.2)') 'MEDIAN  SCORE:    ', self%cls_spec_clsscore_stats(icls_spec)%med
            if( l_print) write(logfhandle,'(a,1x,f8.2)') 'AVERAGE SCORE:    ', self%cls_spec_clsscore_stats(icls_spec)%avg
            if( l_print) write(logfhandle,'(a,1x,f8.2)') 'SDEV    SCORE:    ', self%cls_spec_clsscore_stats(icls_spec)%sdev
            call calc_stats(res(:cnt), self%cls_spec_clsres_stats(icls_spec))
            if( l_print) write(logfhandle,'(a,1x,f8.2)') 'MINIMUM RES:      ', self%cls_spec_clsres_stats(icls_spec)%minv
            if( l_print) write(logfhandle,'(a,1x,f8.2)') 'MAXIMUM RES:      ', self%cls_spec_clsres_stats(icls_spec)%maxv
            if( l_print) write(logfhandle,'(a,1x,f8.2)') 'MEDIAN  RES:      ', self%cls_spec_clsres_stats(icls_spec)%med
            if( l_print) write(logfhandle,'(a,1x,f8.2)') 'AVERAGE RES:      ', self%cls_spec_clsres_stats(icls_spec)%avg
            if( l_print) write(logfhandle,'(a,1x,f8.2)') 'SDEV    RES:      ', self%cls_spec_clsres_stats(icls_spec)%sdev
        end do
    end subroutine calc_cls_spec_stats

    subroutine rank_spectral_classes( self )
        class(pspecs), intent(inout) :: self
        integer :: rank, order(self%ncls_spec), ispec, rank_assign(self%nspecs), icls_spec
        order = (/(ispec,ispec=1,self%ncls_spec)/)
        call hpsort(order, ispec_lt_jspec)
        call reverse(order)
        rank_assign = 0
        do rank = 1, self%ncls_spec
            do ispec = 1, self%nspecs
                if( self%clsinds_spec(ispec) == order(rank) ) rank_assign(ispec) = rank
            end do
        enddo
        self%clsinds_spec = rank_assign
        ! update spectral class populations
        do icls_spec = 1, self%ncls_spec
            self%clspops_spec(icls_spec)  = count(self%clsinds_spec == icls_spec)
            self%ptclpops_spec(icls_spec) = sum(self%clspops, mask=self%clsinds_spec == icls_spec)
        enddo
        
        contains

            function ispec_lt_jspec( ispec, jspec ) result( l_worse )
                integer, intent(in) :: ispec, jspec
                real, parameter :: ABS_SCORE_DIFF_THRES = 0.05
                logical :: l_score_worse, l_res_worse, l_worse
                real    :: abs_score_diff
                l_score_worse = .false.
                if( self%cls_spec_clsscore_stats(ispec)%med < self%cls_spec_clsscore_stats(jspec)%med )then
                    l_score_worse  = .true.
                    abs_score_diff = 0.
                else
                    abs_score_diff = abs(self%cls_spec_clsscore_stats(ispec)%med - self%cls_spec_clsscore_stats(jspec)%med)
                endif
                l_res_worse = .false.
                if( self%cls_spec_clsres_stats(ispec)%med > self%cls_spec_clsres_stats(jspec)%med )then
                    l_res_worse = .true.
                endif
                if( l_score_worse .and. l_res_worse )then
                    l_worse = .true.
                else if( abs_score_diff <= ABS_SCORE_DIFF_THRES .and. l_res_worse )then ! worse resolution trumps similar scores
                    l_worse = .true.
                else
                    l_worse = .false.
                endif
            end function ispec_lt_jspec

    end subroutine rank_spectral_classes

    ! destructor

    subroutine kill( self )
        class(pspecs), intent(inout) :: self
        if( self%exists )then
            if( allocated(self%pspecs)                  ) deallocate(self%pspecs)
            if( allocated(self%pspecs_cen)              ) deallocate(self%pspecs_cen)
            if( allocated(self%pspec_junk)              ) deallocate(self%pspec_junk)
            if( allocated(self%pspec_good)              ) deallocate(self%pspec_good)
            if( allocated(self%pspec_bad)               ) deallocate(self%pspec_bad)
            if( allocated(self%resarr)                  ) deallocate(self%resarr)
            if( allocated(self%dynranges)               ) deallocate(self%dynranges)
            if( allocated(self%dists2cens)              ) deallocate(self%dists2cens)
            if( allocated(self%dists2good)              ) deallocate(self%dists2good)
            if( allocated(self%dists2bad)               ) deallocate(self%dists2bad)
            if( allocated(self%ranks)                   ) deallocate(self%ranks)
            if( allocated(self%order)                   ) deallocate(self%order)
            if( allocated(self%clsinds_spec)            ) deallocate(self%clsinds_spec)
            if( allocated(self%clsinds)                 ) deallocate(self%clsinds)
            if( allocated(self%clspops)                 ) deallocate(self%clspops)
            if( allocated(self%clsres)                  ) deallocate(self%clsres)
            if( allocated(self%clspops_spec)            ) deallocate(self%clspops_spec)
            if( allocated(self%ptclpops_spec)           ) deallocate(self%ptclpops_spec)
            if( allocated(self%l_junk_spec)             ) deallocate(self%l_junk_spec)  
            if( allocated(self%clsscore_stats)          ) deallocate(self%clsscore_stats)
            if( allocated(self%cls_spec_clsscore_stats) ) deallocate(self%cls_spec_clsscore_stats)
            if( allocated(self%cls_spec_clsres_stats)   ) deallocate(self%cls_spec_clsres_stats)
            self%exists = .false.
        endif
    end subroutine kill

end module simple_pspecs
