! regularizer of the cluster2D and refine3D
module simple_regularizer
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_parameters,        only: params_glob
use simple_builder,           only: build_glob
use simple_polarft_corrcalc,  only: polarft_corrcalc
use simple_corr_binfile,      only: corr_binfile
use simple_image
implicit none

public :: regularizer, calc_nrefs2sample
private
#include "simple_local_flags.inc"

type reg_params
    integer :: iptcl        !< iptcl index
    integer :: iref         !< iref index
    integer :: loc          !< inpl index
    real    :: prob         !< probability/corr
    real    :: sh(2)        !< shift
    real    :: w            !< weight
end type reg_params

type :: regularizer
    integer              :: nrots
    integer              :: nrefs
    integer              :: inpl_ns, refs_ns
    real,    allocatable :: corr_loc_tab(:,:,:)         !< 2D corr/loc table
    integer, allocatable :: ptcl_ref_map(:)             !< ptcl -> ref assignment map
    real,    allocatable :: inpl_corr(:,:), refs_corr(:,:)
    integer, allocatable :: inpl_inds(:,:), refs_inds(:,:)
    logical, allocatable :: ptcl_avail(:)
    class(polarft_corrcalc), pointer     :: pftcc => null()
    type(reg_params),        allocatable :: ref_ptcl_tab(:,:)
    contains
    ! CONSTRUCTOR
    procedure          :: new
    ! PROCEDURES
    procedure          :: fill_tab_inpl_smpl
    procedure          :: tab_normalize
    procedure          :: tab_align
    procedure          :: normalize_weight
    procedure          :: shift_search
    procedure          :: write_tab, read_tab
    procedure, private :: ref_multinomal, inpl_multinomal
    ! DESTRUCTOR
    procedure          :: kill
end type regularizer

contains

    ! CONSTRUCTORS

    subroutine new( self, pftcc )
        class(regularizer),      target, intent(inout) :: self
        class(polarft_corrcalc), target, intent(inout) :: pftcc
        real,    allocatable :: dist(:), dist_inpl(:)
        logical, allocatable :: states(:)
        integer :: iptcl, iref
        real    :: dist_thres, athres
        self%nrots = pftcc%nrots
        self%nrefs = pftcc%nrefs
        call calc_nrefs2sample(self%nrefs, self%nrots, params_glob%reg_athres, self%refs_ns, self%inpl_ns)
        self%pftcc => pftcc
        allocate(self%corr_loc_tab(self%nrefs,params_glob%fromp:params_glob%top,2),&
                &self%refs_corr(self%nrefs,params_glob%nthr), self%inpl_corr(self%nrots,params_glob%nthr), source=0.)
        allocate(self%inpl_inds(self%nrots,params_glob%nthr), self%refs_inds(self%nrefs,params_glob%nthr), source=0)
        allocate(self%ref_ptcl_tab(self%nrefs,params_glob%fromp:params_glob%top))
        allocate(self%ptcl_ref_map(params_glob%fromp:params_glob%top))
        allocate(self%ptcl_avail(params_glob%fromp:params_glob%top), source=.true.)
        do iptcl = params_glob%fromp,params_glob%top
            self%ptcl_avail(iptcl) = (build_glob%spproj_field%get_state(iptcl) > 0)
            if( self%ptcl_avail(iptcl) )then
                do iref = 1,self%nrefs
                    self%ref_ptcl_tab(iref,iptcl)%iptcl = iptcl
                    self%ref_ptcl_tab(iref,iptcl)%iref  = iref
                    self%ref_ptcl_tab(iref,iptcl)%loc   = 0
                    self%ref_ptcl_tab(iref,iptcl)%prob  = 0.
                    self%ref_ptcl_tab(iref,iptcl)%sh    = 0.
                    self%ref_ptcl_tab(iref,iptcl)%w     = 1.
                enddo
            endif
        enddo
    end subroutine new

    subroutine fill_tab_inpl_smpl( self, glob_pinds )
        class(regularizer), intent(inout) :: self
        integer,            intent(in)    :: glob_pinds(self%pftcc%nptcls)
        integer :: i, iref, iptcl, irnd
        real    :: inpl_corrs(self%nrots)
        call seed_rnd
        !$omp parallel do collapse(2) default(shared) private(i,iref,iptcl,inpl_corrs,irnd) proc_bind(close) schedule(static)
        do iref = 1, self%nrefs
            do i = 1, self%pftcc%nptcls
                iptcl = glob_pinds(i)
                if( self%ptcl_avail(iptcl) )then
                    ! sampling the inpl rotation
                    call self%pftcc%gencorrs( iref, iptcl, inpl_corrs )
                    irnd = self%inpl_multinomal(inpl_corrs)
                    self%corr_loc_tab(iref,iptcl,1) = inpl_corrs(irnd)
                    self%corr_loc_tab(iref,iptcl,2) = real(irnd)
                endif
            enddo
        enddo
        !$omp end parallel do
    end subroutine fill_tab_inpl_smpl

    subroutine tab_normalize( self )
        class(regularizer), intent(inout) :: self
        integer   :: iref, iptcl
        real      :: sum_corr_all, min_corr, max_corr
        ! normalize so prob of each ptcl is between [0,1] for all refs
        if( params_glob%l_reg_norm )then
            !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iptcl,sum_corr_all)
            do iptcl = params_glob%fromp, params_glob%top
                if( self%ptcl_avail(iptcl) )then
                    sum_corr_all = sum(self%corr_loc_tab(:,iptcl,1))
                    if( sum_corr_all < TINY )then
                        self%corr_loc_tab(:,iptcl,1) = 0.
                    else
                        self%corr_loc_tab(:,iptcl,1) = self%corr_loc_tab(:,iptcl,1) / sum_corr_all
                    endif
                endif
            enddo
            !$omp end parallel do
        endif
        max_corr = 0.
        min_corr = huge(min_corr)
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iptcl)&
        !$omp reduction(min:min_corr) reduction(max:max_corr)
        do iptcl = params_glob%fromp,params_glob%top
            if( self%ptcl_avail(iptcl) )then
                max_corr = max(max_corr, maxval(self%corr_loc_tab(:,iptcl,1), dim=1))
                min_corr = min(min_corr, minval(self%corr_loc_tab(:,iptcl,1), dim=1))
            endif
        enddo
        !$omp end parallel do
        if( (max_corr - min_corr) < TINY )then
            self%corr_loc_tab(:,:,1) = 0.
        else
            self%corr_loc_tab(:,:,1) = (self%corr_loc_tab(:,:,1) - min_corr) / (max_corr - min_corr)
        endif
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iref,iptcl)
        do iptcl = params_glob%fromp,params_glob%top
            if( self%ptcl_avail(iptcl) )then
                do iref = 1, self%nrefs
                    self%ref_ptcl_tab(iref,iptcl)%prob = self%corr_loc_tab(iref,iptcl,1)
                    self%ref_ptcl_tab(iref,iptcl)%loc  = self%corr_loc_tab(iref,iptcl,2)
                enddo
            else
                self%corr_loc_tab(:,iptcl,1) = 1.     ! unselected particles are down on the sorted corrs
            endif
        enddo
        !$omp end parallel do
    end subroutine tab_normalize

    subroutine shift_search( self, glob_pinds )
        use simple_pftcc_shsrch_reg, only: pftcc_shsrch_reg
        class(regularizer), intent(inout) :: self
        integer,            intent(in)    :: glob_pinds(self%pftcc%nptcls)
        type(pftcc_shsrch_reg) :: grad_shsrch_obj(params_glob%nthr)
        integer :: iref, iptcl, ithr, irot, i
        real    :: lims(2,2), cxy(3)
        lims(1,1) = -params_glob%trs
        lims(1,2) =  params_glob%trs
        lims(2,1) = -params_glob%trs
        lims(2,2) =  params_glob%trs
        do ithr = 1, params_glob%nthr
            call grad_shsrch_obj(ithr)%new(lims, opt_angle=.true.)
        enddo
        !$omp parallel do default(shared) private(i,iref,iptcl,irot,ithr,cxy) proc_bind(close) schedule(static)
        do i = 1, self%pftcc%nptcls
            iptcl = glob_pinds(i)
            if( self%ptcl_avail(iptcl) )then
                iref  = self%ptcl_ref_map(iptcl)
                ithr  = omp_get_thread_num() + 1
                call grad_shsrch_obj(ithr)%set_indices(iref, iptcl)
                irot = self%ref_ptcl_tab(iref,iptcl)%loc
                cxy  = grad_shsrch_obj(ithr)%minimize(irot)
                if( irot > 0 )then
                    self%ref_ptcl_tab(iref,iptcl)%sh  = cxy(2:3)
                    self%ref_ptcl_tab(iref,iptcl)%loc = irot
                endif
            endif
        enddo
        !$omp end parallel do
    end subroutine shift_search

    subroutine tab_align( self )
        class(regularizer), intent(inout) :: self
        integer :: iref, iptcl, assigned_iref, assigned_ptcl, &
                  &ref_dist_inds(self%nrefs), stab_inds(params_glob%fromp:params_glob%top, self%nrefs)
        real    :: sorted_tab(params_glob%fromp:params_glob%top, self%nrefs), ref_dist(self%nrefs)
        logical :: ptcl_avail(params_glob%fromp:params_glob%top)
        self%ptcl_ref_map = 1
        ! sorting each columns
        sorted_tab = transpose(self%corr_loc_tab(:,:,1))
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iref,iptcl)
        do iref = 1, self%nrefs
            stab_inds(:,iref) = (/(iptcl,iptcl=params_glob%fromp,params_glob%top)/)
            call hpsort(sorted_tab(:,iref), stab_inds(:,iref))
        enddo
        !$omp end parallel do
        ! first row is the current best ref distribution
        ref_dist_inds = params_glob%fromp
        ref_dist      = sorted_tab(params_glob%fromp,:)
        ptcl_avail    = self%ptcl_avail
        do while( any(ptcl_avail) )
            ! sampling the ref distribution to choose next iref to assign corresponding iptcl to
            assigned_iref = self%ref_multinomal(ref_dist)
            assigned_ptcl = stab_inds(ref_dist_inds(assigned_iref), assigned_iref)
            ptcl_avail(assigned_ptcl)        = .false.
            self%ptcl_ref_map(assigned_ptcl) = assigned_iref
            ! update the ref_dist and ref_dist_inds
            do iref = 1, self%nrefs
                do while( ref_dist_inds(iref) < params_glob%top .and. .not.(ptcl_avail(stab_inds(ref_dist_inds(iref), iref))) )
                    ref_dist_inds(iref) = ref_dist_inds(iref) + 1
                    ref_dist(iref)      = sorted_tab(ref_dist_inds(iref), iref)
                enddo
            enddo
        enddo
    end subroutine tab_align

    subroutine normalize_weight( self )
        class(regularizer), intent(inout) :: self
        integer :: iptcl, iref
        real    :: min_w, max_w
        min_w = huge(min_w)
        max_w = 0.
        do iptcl = params_glob%fromp, params_glob%top
            if( self%ptcl_avail(iptcl) )then
                iref = self%ptcl_ref_map(iptcl)
                self%ref_ptcl_tab(iref,iptcl)%w = 1. - self%ref_ptcl_tab(iref,iptcl)%prob
                if( self%ref_ptcl_tab(iref,iptcl)%w < min_w ) min_w = self%ref_ptcl_tab(iref,iptcl)%w
                if( self%ref_ptcl_tab(iref,iptcl)%w > max_w ) max_w = self%ref_ptcl_tab(iref,iptcl)%w
            endif
        enddo
        do iptcl = params_glob%fromp, params_glob%top
            if( self%ptcl_avail(iptcl) )then
                iref = self%ptcl_ref_map(iptcl)
                self%ref_ptcl_tab(iref,iptcl)%w = (self%ref_ptcl_tab(iref,iptcl)%w - min_w) / (max_w - min_w)
            endif
        enddo
    end subroutine normalize_weight

    !>  \brief  generates a multinomal 1-of-K random number according to the
    !!          distribution in pvec
    function ref_multinomal( self, pvec ) result( which )
        class(regularizer), intent(inout) :: self
        real,               intent(in)    :: pvec(:) !< probabilities
        integer :: i, which, ithr
        real    :: rnd, bound, sum_refs_corr
        ithr = omp_get_thread_num() + 1
        rnd  = ran3()
        self%refs_corr(:,ithr) = pvec
        self%refs_inds(:,ithr) = (/(i,i=1,self%nrefs)/)
        call hpsort(self%refs_corr(:,ithr), self%refs_inds(:,ithr) )
        sum_refs_corr = sum(self%refs_corr(1:self%refs_ns,ithr))
        if( sum_refs_corr < TINY )then
            ! uniform sampling
            which = 1 + floor(real(self%refs_ns) * rnd)
        else
            ! normalizing within the hard-limit
            self%refs_corr(1:self%refs_ns,ithr) = self%refs_corr(1:self%refs_ns,ithr) / sum_refs_corr
            bound = 0.
            do which=1,self%refs_ns
                bound = bound + self%refs_corr(which, ithr)
                if( rnd >= bound )exit
            enddo
            which = min(which, self%refs_ns)
        endif
        which = self%refs_inds(which, ithr)
    end function ref_multinomal

    ! inpl multinomal based on unnormalized pvec
    function inpl_multinomal( self, pvec ) result( which )
        class(regularizer), intent(inout) :: self
        real,               intent(in)    :: pvec(:) !< probabilities
        integer :: i, which, ithr
        real    :: rnd, bound, sum_corr
        ithr = omp_get_thread_num() + 1
        self%inpl_corr(:,ithr) = pvec
        self%inpl_inds(:,ithr) = (/(i,i=1,self%nrots)/)
        call hpsort(self%inpl_corr(:,ithr), self%inpl_inds(:,ithr) )
        rnd      = ran3()
        sum_corr = sum(self%inpl_corr(1:self%inpl_ns,ithr))
        if( sum_corr < TINY )then
            ! uniform sampling
            which = 1 + floor(real(self%inpl_ns) * rnd)
        else
            ! normalizing within the hard-limit
            self%inpl_corr(1:self%inpl_ns,ithr) = self%inpl_corr(1:self%inpl_ns,ithr) / sum_corr
            bound = 0.
            do which=1,self%inpl_ns
                bound = bound + self%inpl_corr(which, ithr)
                if( rnd >= bound )exit
            enddo
            which = min(which,self%inpl_ns)
        endif
        which = self%inpl_inds(which, ithr)
    end function inpl_multinomal

    ! FILE IO
    subroutine write_tab( self, binfname )
        class(regularizer), intent(in) :: self
        character(len=*),   intent(in) :: binfname
        type(corr_binfile) :: binfile
        if( file_exists(binfname) )then
            call binfile%new_from_file(binfname)
        else
            call binfile%new(binfname, params_glob%fromp, params_glob%top, self%nrefs)
        endif
        call binfile%write(self%corr_loc_tab)
        call binfile%kill
    end subroutine write_tab

    subroutine read_tab( self, binfname )
        class(regularizer), intent(inout) :: self
        character(len=*),   intent(in)    :: binfname
        type(corr_binfile) :: binfile
        if( file_exists(binfname) )then
            call binfile%new_from_file(binfname)
        else
            call binfile%new(binfname, params_glob%fromp, params_glob%top, self%nrefs)
        endif
        call binfile%read(self%corr_loc_tab)
        call binfile%kill
    end subroutine read_tab

    ! DESTRUCTOR

    subroutine kill( self )
        class(regularizer), intent(inout) :: self
        deallocate(self%corr_loc_tab,self%ref_ptcl_tab,self%ptcl_ref_map,self%inpl_corr,self%refs_corr,&
                  &self%inpl_inds,self%refs_inds,self%ptcl_avail)
    end subroutine kill

    ! PUBLIC UTILITITES

    subroutine calc_nrefs2sample( nrefs, nrots, reg_athres, refs_ns, inpl_ns)
        integer, intent(in)  :: nrefs, nrots
        real,    intent(in)  :: reg_athres
        integer, intent(out) :: refs_ns, inpl_ns
        real,    allocatable :: vals(:)
        logical, allocatable :: ptcl_mask(:)
        real    :: athres, dist_thres
        integer :: n
        ptcl_mask  = nint(build_glob%spproj_field%get_all('state')) == 1
        n          = count(ptcl_mask)
        ! in-planes
        vals       = build_glob%spproj_field%get_all('dist_inpl')
        dist_thres = sum(vals, mask=ptcl_mask) / real(n)
        athres     = reg_athres
        if( dist_thres > TINY ) athres = min(athres, dist_thres)
        inpl_ns    = min(nrots,max(1,int(athres * real(nrots) / 180.)))
        ! projection directions
        vals       = build_glob%spproj_field%get_all('dist')
        dist_thres = sum(vals,mask=ptcl_mask) / real(n)
        athres     = reg_athres
        if( dist_thres > TINY ) athres = min(athres, dist_thres)
        refs_ns    = min(nrefs,max(1,int(athres * real(nrefs) / 180.)))
    end subroutine

end module simple_regularizer
