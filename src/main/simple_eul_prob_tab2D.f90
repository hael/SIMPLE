!@descr: 2D probability table routines for multi-reference class assignment with probabilistic sampling
module simple_eul_prob_tab2D
use, intrinsic :: iso_fortran_env, only: int64
use simple_pftc_srch_api
use simple_builder,            only: builder
use simple_pftc_shsrch_grad,   only: pftc_shsrch_grad
use simple_decay_funs,         only: extremal_decay2D
use simple_eul_prob_tab_utils, only: build_pind_lookup, eulprob_dist_switch, materialize_seed_shift,&
    &read_seed_shift_table, write_seed_shift_table, sample_likelihood_index
use simple_segmentation,       only: detect_peak_thres_fdr
implicit none

public :: eul_prob_tab2D, PRIOR2D_STAGE5_FNAME
private
#include "simple_local_flags.inc"

real,             parameter :: NHOOD_FRAC           = 0.1
real,             parameter :: PRIOR2D_TOPK_FRAC    = 0.3  !< fraction of nclasses stored per particle as prior top-K
real,             parameter :: PRIOR2D_FDR_Q        = 0.25 !< FDR q for dynamic prob_prior neighborhoods
character(len=*), parameter :: PRIOR2D_STAGE5_FNAME = 'posterior_topk_stage05.dat'

type :: eul_prob_tab2D
    class(builder),    pointer  :: b_ptr => null()
    class(parameters), pointer  :: p_ptr => null()
    type(ptcl_ref), allocatable :: loc_tab(:,:)      !< 2D search table (nclasses, nptcls)
    type(ptcl_ref), allocatable :: assgn_map(:)      !< assignment map (nptcls)
    real,           allocatable :: seed_shifts(:,:)  !< per-particle seeded shift (2,nptcls)
    logical,        allocatable :: seed_has_sh(:)    !< per-particle seeded shift flag
    integer                     :: seed_nrots = 0    !< rotation grid used for deferred seeded shifts
    integer,        allocatable :: pinds(:)          !< particle indices for processing
    logical,        allocatable :: class_exists(:)   !< class population filter
    integer,        allocatable :: eval_touched_refs(:,:)
    integer,        allocatable :: eval_touched_counts(:)
    integer                     :: eval_max_touched = 0
    logical                     :: l_sparse_snhc = .false.
    logical                     :: l_prior_mode  = .false.
    logical,        allocatable :: prior_operated(:)
    integer                     :: prior_kmax = 0
    integer                     :: nptcls            !< size of pinds array
    integer                     :: nclasses          !< number of classes
    integer                     :: nhood_sz = 1      !< probabilistic neighborhood size used in fill_tab/ref_assign
    contains
    ! CONSTRUCTOR
    procedure :: new
    ! MAIN PROCEDURES
    procedure :: fill_tab
    procedure :: fill_tab_range
    procedure :: ref_assign => ref_assign_likelihood
    procedure :: write_tab
    procedure :: read_tab_to_glob
    procedure :: write_assignment
    procedure :: read_assignment
    procedure :: write_prior_topk
    ! DESTRUCTOR
    procedure :: kill
    ! PRIVATE
    procedure, private :: fill_tab_prob_snhc
    procedure, private :: fill_tab_prob_snhc_range
    procedure, private :: ref_assign_sparse_likelihood
    procedure, private :: clear_sparse_eval_table_range
    procedure, private :: record_sparse_eval
    procedure, private :: mark_ref_touched
    procedure, private :: write_sparse_tab
    procedure, private :: read_sparse_tab_to_glob
end type eul_prob_tab2D

type :: eval2D_sparse_ws
    integer, allocatable :: direct_srch_order(:,:) ! [nclasses,nthr]
    integer, allocatable :: vec_nrots(:,:)         ! [nrots,nthr]
    integer, allocatable :: eval_cls(:,:)          ! [nclasses,nthr]
    integer, allocatable :: best_locs(:,:)         ! [smpl_ncls,nthr]
    real,    allocatable :: inpl_corrs(:,:)        ! [nrots,nthr]
    real,    allocatable :: eval_dists(:,:)        ! [nclasses,nthr]
contains
    procedure :: init_eval2D_sparse_ws
    procedure :: dealloc_eval2D_sparse_ws
end type eval2D_sparse_ws

contains

    ! CONSTRUCTORS

    subroutine new( self, params, build, pinds )
        class(eul_prob_tab2D),     intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        class(builder),    target, intent(in)    :: build
        integer,                   intent(in)    :: pinds(:)
        integer :: i, icls, iptcl, nactive
        real    :: x
        call self%kill
        self%p_ptr     => params
        self%b_ptr     => build
        self%nptcls    = size(pinds)
        self%nclasses  = params%ncls
        self%l_sparse_snhc = trim(params%refine) == 'prob_snhc' .or. trim(params%refine) == 'prob_prior'
        self%l_prior_mode  = trim(params%refine) == 'prob_prior'
        allocate(self%class_exists(self%nclasses))
        ! In 2D probabilistic assignment classes must be able to recover, otherwise
        ! low-population classes are permanently excluded and the solution collapses.
        self%class_exists = .true.
        nactive = count(self%class_exists)
        if( nactive == 0 ) nactive = self%nclasses
        self%nhood_sz = max(1, ceiling(NHOOD_FRAC * real(nactive)))
        self%nhood_sz = min(self%nhood_sz, nactive)
        self%nhood_sz = min(self%nhood_sz, params%npeaks_inpl)
        allocate(self%pinds(self%nptcls), source=pinds)
        allocate(self%loc_tab(self%nclasses, self%nptcls), self%assgn_map(self%nptcls))
        allocate(self%seed_shifts(2,self%nptcls), source=0.)
        allocate(self%seed_has_sh(self%nptcls), source=.false.)
        if( self%l_sparse_snhc )then
            self%eval_max_touched = max(1, self%nclasses)
            allocate(self%eval_touched_refs(self%eval_max_touched,self%nptcls), source=0)
            allocate(self%eval_touched_counts(self%nptcls), source=0)
        endif
        !$omp parallel do default(shared) private(i,iptcl,icls) proc_bind(close) schedule(static)
        do i = 1, self%nptcls
            iptcl = self%pinds(i)
            self%assgn_map(i)%pind   = iptcl
            self%assgn_map(i)%icls   = 0
            self%assgn_map(i)%istate = 0
            self%assgn_map(i)%inpl   = 0
            self%assgn_map(i)%dist   = huge(x)
            self%assgn_map(i)%x      = 0.
            self%assgn_map(i)%y      = 0.
            self%assgn_map(i)%has_sh = .false.
            self%assgn_map(i)%frac   = 100.
            self%assgn_map(i)%npeaks = 0
            do icls = 1, self%nclasses
                self%loc_tab(icls,i)%pind   = iptcl
                self%loc_tab(icls,i)%icls   = icls
                self%loc_tab(icls,i)%inpl   = 0
                self%loc_tab(icls,i)%dist   = huge(x)
                self%loc_tab(icls,i)%x      = 0.
                self%loc_tab(icls,i)%y      = 0.
                self%loc_tab(icls,i)%has_sh = .false.
                self%loc_tab(icls,i)%frac   = 100.
                self%loc_tab(icls,i)%npeaks = 0
            end do
        end do
        !$omp end parallel do
    end subroutine new

    ! table filling for 2D multi-class assignment
    subroutine fill_tab( self )
        class(eul_prob_tab2D), intent(inout) :: self
        call self%fill_tab_range(1, self%nptcls)
    end subroutine fill_tab

    subroutine fill_tab_range( self, i_first, i_last )
        class(eul_prob_tab2D), intent(inout) :: self
        integer,               intent(in)    :: i_first, i_last
        type(pftc_shsrch_grad) :: grad_shsrch_obj(nthr_glob)  !< shift search object, L-BFGS with gradient
        type(ori)              :: o_prev
        integer :: i, icls, iptcl, ithr, irot, irot0, icls_prev, iref_n, nactive, nhood_sz_loc, ninpl_smpl
        integer :: i_from, i_to
        integer :: active_cls(self%nclasses)
        integer :: loc(1), vec_nrots(self%b_ptr%pftc%get_nrots())
        real    :: lims(2,2), lims_init(2,2), cxy(3), cxy_prob(3)
        real    :: inpl_corrs(self%b_ptr%pftc%get_nrots()), cls_dists(self%nclasses), cls_dists_work(self%nclasses)
        real    :: inpl_corr, inpl_dist, neigh_frac
        logical :: class_active(self%nclasses)
        if( self%l_sparse_snhc )then
            call self%fill_tab_prob_snhc_range(i_first, i_last)
            return
        endif
        i_from = max(1, i_first)
        i_to   = min(self%nptcls, i_last)
        if( i_to < i_from ) return
        call seed_rnd
        self%seed_nrots = self%b_ptr%pftc%get_nrots()
        class_active = self%class_exists
        nactive      = count(class_active)
        if( nactive == 0 )then
            class_active = .true.
            nactive      = self%nclasses
            THROW_WARN('No active classes after population filtering; falling back to all classes in eul_prob_tab2D')
        endif
        nactive = 0
        do icls = 1, self%nclasses
            if( class_active(icls) )then
                nactive = nactive + 1
                active_cls(nactive) = icls
            endif
        end do
        nhood_sz_loc = min(self%nhood_sz, nactive)
        neigh_frac = extremal_decay2D(self%p_ptr%extr_iter, self%p_ptr%extr_lim)
        ninpl_smpl = neighfrac2nsmpl(neigh_frac, self%b_ptr%pftc%get_nrots())
        ninpl_smpl = max(1, min(ninpl_smpl, self%b_ptr%pftc%get_nrots()))
        if( self%p_ptr%l_doshift )then
            ! make shift search objects
            lims(:,1)      = -self%p_ptr%trs
            lims(:,2)      =  self%p_ptr%trs
            lims_init(:,1) = -SHC_INPL_TRSHWDTH
            lims_init(:,2) =  SHC_INPL_TRSHWDTH
            do ithr = 1, nthr_glob
                call grad_shsrch_obj(ithr)%new(self%b_ptr, lims, lims_init=lims_init, shbarrier=self%p_ptr%shbarrier,&
                    &maxits=self%p_ptr%maxits_sh, opt_angle=.true., coarse_init=.true.)
            end do
            ! fill the table
            !$omp parallel do default(shared) private(i,iptcl,ithr,o_prev,irot,irot0,cxy,icls,icls_prev,inpl_corrs,loc,&
            !$omp& cls_dists,cls_dists_work,iref_n,cxy_prob,vec_nrots,inpl_corr,inpl_dist)&
            !$omp proc_bind(close) schedule(static)
            do i = i_from, i_to
                iptcl = self%pinds(i)
                ithr  = omp_get_thread_num() + 1
                ! (1) identify shifts using the previously assigned best class
                call self%b_ptr%spproj_field%get_ori(iptcl, o_prev)      ! previous ori
                irot0 = self%b_ptr%pftc%get_roind(360. - o_prev%e3get()) ! in-plane angle index seed
                icls_prev = nint(self%b_ptr%spproj_field%get(iptcl, 'class'))
                cxy       = 0.
                if( icls_prev >= 1 .and. icls_prev <= self%nclasses )then
                    if( class_active(icls_prev) )then
                        irot = irot0
                        call grad_shsrch_obj(ithr)%set_indices(icls_prev, iptcl)
                        cxy = grad_shsrch_obj(ithr)%minimize(irot=irot, sh_rot=.false.)
                        if( irot == 0 )then
                            cxy(2:3) = 0.
                            irot = irot0
                        endif
                    endif
                endif
                self%seed_shifts(:,i) = cxy(2:3)
                self%seed_has_sh(i)   = .true.
                ! (2) search class references using that shared shift
                cls_dists = huge(1.0)
                do iref_n = 1, nactive
                    icls = active_cls(iref_n)
                    call self%b_ptr%pftc%gen_prob_likelihood_objfun_val(icls, iptcl, cxy(2:3), ninpl_smpl,&
                        &inpl_dist, inpl_corr, irot, inpl_corrs, vec_nrots)
                    self%loc_tab(icls,i)%dist = inpl_dist
                    self%loc_tab(icls,i)%inpl   = irot
                    cls_dists(icls) = self%loc_tab(icls,i)%dist
                end do
                ! (3) refine shifts for a neighborhood of classes around the best one
                cls_dists_work = cls_dists
                do iref_n = 1, nhood_sz_loc
                    loc = minloc(cls_dists_work)
                    icls = loc(1)
                    if( cls_dists_work(icls) >= huge(1.0)/2.0 ) exit
                    call grad_shsrch_obj(ithr)%set_indices(icls, iptcl)
                    irot     = self%loc_tab(icls,i)%inpl
                    cxy_prob = grad_shsrch_obj(ithr)%minimize(irot=irot, sh_rot=.true., xy_in=cxy(2:3))
                    if( irot > 0 )then
                        self%loc_tab(icls,i)%inpl   = irot
                        self%loc_tab(icls,i)%dist   = eulprob_dist_switch(cxy_prob(1), self%p_ptr%cc_objfun)
                        self%loc_tab(icls,i)%x      = cxy_prob(2)
                        self%loc_tab(icls,i)%y      = cxy_prob(3)
                        self%loc_tab(icls,i)%has_sh = .true.
                    endif
                    cls_dists_work(icls) = huge(1.0)
                end do
            end do
            !$omp end parallel do
        else
            ! shift-free path: evaluate classes/in-planes at zero shift
            !$omp parallel do default(shared) private(i,iptcl,ithr,icls,iref_n,inpl_corrs,vec_nrots,irot,&
            !$omp& inpl_corr,inpl_dist) proc_bind(close) schedule(static)
            do i = i_from, i_to
                iptcl = self%pinds(i)
                ithr  = omp_get_thread_num() + 1
                do iref_n = 1, nactive
                    icls = active_cls(iref_n)
                    call self%b_ptr%pftc%gen_prob_likelihood_objfun_val(icls, iptcl, [0.,0.], ninpl_smpl,&
                        &inpl_dist, inpl_corr, irot, inpl_corrs, vec_nrots)
                    self%loc_tab(icls,i)%dist = inpl_dist
                    self%loc_tab(icls,i)%inpl   = irot
                    self%loc_tab(icls,i)%x      = 0.
                    self%loc_tab(icls,i)%y      = 0.
                    self%loc_tab(icls,i)%has_sh = .false.
                end do
            end do
            !$omp end parallel do
        endif
        do ithr = 1, nthr_glob
            call grad_shsrch_obj(ithr)%kill
        end do
        call o_prev%kill
    end subroutine fill_tab_range

    subroutine fill_tab_prob_snhc( self )
        class(eul_prob_tab2D), intent(inout) :: self
        call self%fill_tab_prob_snhc_range(1, self%nptcls)
    end subroutine fill_tab_prob_snhc

    subroutine fill_tab_prob_snhc_range( self, i_first, i_last )
        class(eul_prob_tab2D), intent(inout) :: self
        integer,               intent(in)    :: i_first, i_last
        type(pftc_shsrch_grad) :: grad_shsrch_obj(nthr_glob)
        type(eval2D_sparse_ws) :: eval_work
        type(ran_tabu), allocatable :: direct_rts(:)
        integer, allocatable :: prior_units(:)
        integer :: ithr, i, iptcl, nrefs_bound, smpl_ncls, ninpl_smpl, nrots
        integer :: i_from, i_to
        real    :: lims(2,2), lims_init(2,2), neigh_frac, cxy(3)
        logical :: prior_file_ok
        i_from = max(1, i_first)
        i_to   = min(self%nptcls, i_last)
        if( i_to < i_from ) return
        call seed_rnd
        nrots = self%b_ptr%pftc%get_nrots()
        self%seed_nrots = nrots
        neigh_frac  = extremal_decay2D(self%p_ptr%extr_iter, self%p_ptr%extr_lim)
        nrefs_bound = max(2, min(self%nclasses, nint(real(self%nclasses) * (1. - neigh_frac))))
        nrefs_bound = max(1, min(nrefs_bound, self%nclasses))
        smpl_ncls   = neighfrac2nsmpl(neigh_frac, self%nclasses)
        smpl_ncls   = max(1, min(smpl_ncls, self%nclasses))
        ninpl_smpl  = neighfrac2nsmpl(neigh_frac, nrots)
        ninpl_smpl  = max(1, min(ninpl_smpl, nrots))
        if( .not. allocated(self%eval_touched_refs) )then
            self%eval_max_touched = max(1, self%nclasses)
            allocate(self%eval_touched_refs(self%eval_max_touched,self%nptcls), source=0)
        endif
        if( .not. allocated(self%eval_touched_counts) ) allocate(self%eval_touched_counts(self%nptcls), source=0)
        call self%clear_sparse_eval_table_range(i_from, i_to)
        call eval_work%init_eval2D_sparse_ws(self%nclasses, nrots, smpl_ncls, nthr_glob)
        allocate(direct_rts(nthr_glob))
        do ithr = 1, nthr_glob
            direct_rts(ithr) = ran_tabu(self%nclasses)
        enddo
        prior_file_ok = .false.
        if( self%l_prior_mode )then
            self%prior_kmax = max(1, nint(PRIOR2D_TOPK_FRAC * real(self%nclasses)))
            if( allocated(self%prior_operated) )then
                self%prior_operated = .false.
            else
                allocate(self%prior_operated(self%nptcls), source=.false.)
            endif
            call open_prior_units(prior_units, prior_file_ok)
        endif
        if( self%p_ptr%l_doshift )then
            lims(:,1)      = -self%p_ptr%trs
            lims(:,2)      =  self%p_ptr%trs
            lims_init(:,1) = -SHC_INPL_TRSHWDTH
            lims_init(:,2) =  SHC_INPL_TRSHWDTH
            do ithr = 1, nthr_glob
                call grad_shsrch_obj(ithr)%new(self%b_ptr, lims, lims_init=lims_init, shbarrier=self%p_ptr%shbarrier,&
                    &maxits=self%p_ptr%maxits_sh, opt_angle=.true., coarse_init=.true.)
            end do
            !$omp parallel do default(shared) private(i,iptcl,ithr,cxy) proc_bind(close) schedule(static)
            do i = i_from, i_to
                iptcl = self%pinds(i)
                ithr  = omp_get_thread_num() + 1
                call process_particle(i, iptcl, ithr, .true., cxy)
                self%seed_shifts(:,i) = cxy(2:3)
                self%seed_has_sh(i)   = .true.
            enddo
            !$omp end parallel do
        else
            !$omp parallel do default(shared) private(i,iptcl,ithr,cxy) proc_bind(close) schedule(static)
            do i = i_from, i_to
                iptcl = self%pinds(i)
                ithr  = omp_get_thread_num() + 1
                cxy   = 0.
                call process_particle(i, iptcl, ithr, .false., cxy)
            enddo
            !$omp end parallel do
        endif
        do ithr = 1, nthr_glob
            call grad_shsrch_obj(ithr)%kill
        end do
        if( allocated(prior_units) ) call close_prior_units(prior_units)
        call eval_work%dealloc_eval2D_sparse_ws
        deallocate(direct_rts)

    contains

        subroutine process_particle(i_loc, iptcl_loc, ithr_loc, l_with_shift, cxy_loc)
            integer, intent(in)    :: i_loc, iptcl_loc, ithr_loc
            logical, intent(in)    :: l_with_shift
            real,    intent(inout) :: cxy_loc(3)
            type(ori) :: o_prev_loc
            integer :: icls_prev_loc, irot0_loc, irot_loc, isample_loc, icls_loc, neval_loc
            integer :: prior_nloc, neval_bound_loc
            real    :: dist_loc, corr_loc, sh_loc(2)
            prior_nloc = 0
            call direct_rts(ithr_loc)%ne_ran_iarr(eval_work%direct_srch_order(:,ithr_loc))
            if( self%l_prior_mode .and. prior_file_ok )then
                call apply_prior_order(ithr_loc, iptcl_loc, eval_work%direct_srch_order(:,ithr_loc), prior_nloc)
            endif
            if( self%l_prior_mode .and. allocated(self%prior_operated) ) self%prior_operated(i_loc) = prior_nloc > 0
            call self%b_ptr%spproj_field%get_ori(iptcl_loc, o_prev_loc)
            icls_prev_loc = nint(self%b_ptr%spproj_field%get(iptcl_loc, 'class'))
            ! restore put_last for unsampled particles (prior_nloc==0): they get prob_snhc behaviour
            if( ((.not. self%l_prior_mode) .or. prior_nloc == 0) .and. &
                &icls_prev_loc >= 1 .and. icls_prev_loc <= self%nclasses )then
                call put_last(icls_prev_loc, eval_work%direct_srch_order(:,ithr_loc))
            endif
            cxy_loc = 0.
            if( l_with_shift .and. icls_prev_loc >= 1 .and. icls_prev_loc <= self%nclasses )then
                irot0_loc = self%b_ptr%pftc%get_roind(360. - o_prev_loc%e3get())
                irot_loc  = irot0_loc
                call grad_shsrch_obj(ithr_loc)%set_indices(icls_prev_loc, iptcl_loc)
                cxy_loc = grad_shsrch_obj(ithr_loc)%minimize(irot=irot_loc, sh_rot=.false.)
                if( irot_loc == 0 ) cxy_loc(2:3) = 0.
            endif
            sh_loc = 0.
            if( l_with_shift ) sh_loc = cxy_loc(2:3)
            ! prior particles evaluate exactly the prior-K classes for acceleration;
            ! unsampled particles evaluate the full stochastic nrefs_bound
            if( prior_nloc > 0 )then
                neval_bound_loc = prior_nloc
            else
                neval_bound_loc = nrefs_bound
            endif
            eval_work%eval_cls(:,ithr_loc)   = 0
            eval_work%eval_dists(:,ithr_loc) = huge(1.0)
            neval_loc = 0
            do isample_loc = 1, self%nclasses
                if( isample_loc > neval_bound_loc ) exit
                icls_loc = eval_work%direct_srch_order(isample_loc,ithr_loc)
                call score_class(icls_loc, ithr_loc, iptcl_loc, sh_loc, dist_loc, corr_loc, irot_loc)
                call self%record_sparse_eval(i_loc, icls_loc, dist_loc, irot_loc, 0., 0., .false.)
                neval_loc = neval_loc + 1
                eval_work%eval_cls(neval_loc,ithr_loc)   = icls_loc
                eval_work%eval_dists(neval_loc,ithr_loc) = dist_loc
            enddo
            call refine_best_classes(i_loc, ithr_loc, iptcl_loc, sh_loc, neval_loc, l_with_shift)
            call o_prev_loc%kill
        end subroutine process_particle

        subroutine score_class(icls_loc, ithr_loc, iptcl_loc, sh_loc, dist_loc, corr_loc, irot_loc)
            integer, intent(in)  :: icls_loc, ithr_loc, iptcl_loc
            real,    intent(in)  :: sh_loc(2)
            real,    intent(out) :: dist_loc, corr_loc
            integer, intent(out) :: irot_loc
            call self%b_ptr%pftc%gen_prob_likelihood_objfun_val(icls_loc, iptcl_loc, sh_loc, ninpl_smpl,&
                &dist_loc, corr_loc, irot_loc, eval_work%inpl_corrs(:,ithr_loc), eval_work%vec_nrots(:,ithr_loc))
            if( irot_loc < 1 .or. irot_loc > nrots ) irot_loc = 1
        end subroutine score_class

        subroutine refine_best_classes(i_loc, ithr_loc, iptcl_loc, sh_loc, neval_loc, l_with_shift)
            integer, intent(in) :: i_loc, ithr_loc, iptcl_loc, neval_loc
            real,    intent(in) :: sh_loc(2)
            logical, intent(in) :: l_with_shift
            integer :: nrefine_loc, j_loc, eval_slot_loc, icls_loc, irot_loc
            real    :: cxy_prob_loc(3)
            if( .not. l_with_shift ) return
            if( neval_loc < 1 ) return
            nrefine_loc = min(smpl_ncls, neval_loc)
            eval_work%best_locs(1:nrefine_loc,ithr_loc) = minnloc(eval_work%eval_dists(1:neval_loc,ithr_loc), nrefine_loc)
            do j_loc = 1, nrefine_loc
                eval_slot_loc = eval_work%best_locs(j_loc,ithr_loc)
                if( eval_slot_loc < 1 ) cycle
                icls_loc = eval_work%eval_cls(eval_slot_loc,ithr_loc)
                if( icls_loc < 1 ) cycle
                call grad_shsrch_obj(ithr_loc)%set_indices(icls_loc, iptcl_loc)
                irot_loc = self%loc_tab(icls_loc,i_loc)%inpl
                cxy_prob_loc = grad_shsrch_obj(ithr_loc)%minimize(irot=irot_loc, sh_rot=.true., xy_in=sh_loc)
                if( irot_loc > 0 )then
                    call self%record_sparse_eval(i_loc, icls_loc, eulprob_dist_switch(cxy_prob_loc(1), self%p_ptr%cc_objfun),&
                        &irot_loc, cxy_prob_loc(2), cxy_prob_loc(3), .true.)
                endif
            enddo
        end subroutine refine_best_classes

        subroutine open_prior_units(units, ok)
            integer, allocatable, intent(out) :: units(:)
            logical,              intent(out) :: ok
            integer, allocatable :: rec_buf(:)
            integer :: reclen, ios, j
            ok = .false.
            if( .not. file_exists(string(PRIOR2D_STAGE5_FNAME)) )then
                THROW_WARN('prob_prior: prior file missing, fallback to stochastic sparse search: '//PRIOR2D_STAGE5_FNAME)
                return
            endif
            allocate(rec_buf(self%prior_kmax), source=0)
            inquire(iolength=reclen) rec_buf
            deallocate(rec_buf)
            allocate(units(nthr_glob), source=-1)
            do j = 1, nthr_glob
                open(newunit=units(j), file=PRIOR2D_STAGE5_FNAME, status='OLD', action='READ', &
                    &access='DIRECT', form='UNFORMATTED', recl=reclen, iostat=ios)
                if( ios /= 0 )then
                    call close_prior_units(units)
                    THROW_WARN('prob_prior: failed to open prior direct-access file, fallback to stochastic sparse search')
                    return
                endif
            enddo
            ok = .true.
        end subroutine open_prior_units

        subroutine close_prior_units(units)
            integer, allocatable, intent(inout) :: units(:)
            integer :: j
            if( .not. allocated(units) ) return
            do j = 1, size(units)
                if( units(j) /= -1 ) close(units(j))
            enddo
            deallocate(units)
        end subroutine close_prior_units

        subroutine apply_prior_order(ithr_loc, iptcl_loc, order, prior_nloc)
            integer, intent(in)    :: ithr_loc, iptcl_loc
            integer, intent(inout) :: order(:)
            integer, intent(out)   :: prior_nloc
            integer :: prior_rec(self%prior_kmax)
            integer :: ios, k_prior, icls_loc, j_order, i_tmp
            prior_nloc = 0
            prior_rec  = 0
            if( iptcl_loc < 1 ) return
            if( .not. allocated(prior_units) ) return
            read(prior_units(ithr_loc), rec=iptcl_loc, iostat=ios) prior_rec
            if( ios /= 0 ) return
            j_order = 1
            do k_prior = 1, self%prior_kmax
                icls_loc = prior_rec(k_prior)
                if( icls_loc < 1 .or. icls_loc > self%nclasses ) cycle
                do i_tmp = j_order, self%nclasses
                    if( order(i_tmp) == icls_loc )then
                        order(i_tmp) = order(j_order)
                        order(j_order) = icls_loc
                        j_order = j_order + 1
                        prior_nloc = prior_nloc + 1
                        exit
                    endif
                enddo
                if( j_order > self%nclasses ) exit
            enddo
        end subroutine apply_prior_order

    end subroutine fill_tab_prob_snhc_range

    subroutine clear_sparse_eval_table_range( self, i_from, i_to )
        class(eul_prob_tab2D), intent(inout) :: self
        integer,               intent(in)    :: i_from, i_to
        integer :: i, k, icls
        if( .not. allocated(self%eval_touched_counts) .or. .not. allocated(self%eval_touched_refs) ) return
        !$omp parallel do default(shared) private(i,k,icls) proc_bind(close) schedule(static)
        do i = i_from, i_to
            do k = 1, self%eval_touched_counts(i)
                icls = self%eval_touched_refs(k,i)
                if( icls < 1 .or. icls > self%nclasses ) cycle
                self%loc_tab(icls,i)%dist   = huge(1.0)
                self%loc_tab(icls,i)%inpl   = 0
                self%loc_tab(icls,i)%x      = 0.
                self%loc_tab(icls,i)%y      = 0.
                self%loc_tab(icls,i)%has_sh = .false.
            enddo
            self%eval_touched_counts(i) = 0
        enddo
        !$omp end parallel do
    end subroutine clear_sparse_eval_table_range

    subroutine record_sparse_eval( self, i, icls, dist, irot, x, y, has_sh )
        class(eul_prob_tab2D), intent(inout) :: self
        integer,               intent(in)    :: i, icls, irot
        real,                  intent(in)    :: dist, x, y
        logical,               intent(in)    :: has_sh
        if( icls < 1 .or. icls > self%nclasses )&
            &THROW_HARD('simple_eul_prob_tab2D::record_sparse_eval; class index out of range')
        self%loc_tab(icls,i)%dist   = dist
        self%loc_tab(icls,i)%inpl   = irot
        self%loc_tab(icls,i)%x      = x
        self%loc_tab(icls,i)%y      = y
        self%loc_tab(icls,i)%has_sh = has_sh
        call self%mark_ref_touched(i, icls)
    end subroutine record_sparse_eval

    subroutine mark_ref_touched( self, i, icls )
        class(eul_prob_tab2D), intent(inout) :: self
        integer,               intent(in)    :: i, icls
        integer :: kt
        if( .not. allocated(self%eval_touched_counts) .or. .not. allocated(self%eval_touched_refs) )&
            &THROW_HARD('simple_eul_prob_tab2D::mark_ref_touched; sparse bookkeeping is not allocated')
        do kt = 1, self%eval_touched_counts(i)
            if( self%eval_touched_refs(kt,i) == icls ) return
        enddo
        if( self%eval_touched_counts(i) >= self%eval_max_touched )then
            THROW_HARD('simple_eul_prob_tab2D::mark_ref_touched; eval_touched overflow')
        endif
        self%eval_touched_counts(i) = self%eval_touched_counts(i) + 1
        self%eval_touched_refs(self%eval_touched_counts(i),i) = icls
    end subroutine mark_ref_touched


    subroutine ref_assign_sparse_likelihood( self )
        class(eul_prob_tab2D), intent(inout) :: self

        type :: assign_graph_ws
            integer, allocatable :: class_counts(:), class_offsets(:), class_fill(:), class_list(:), class_pos(:), active_classes(:)
            integer, allocatable :: ptcl_counts(:), ptcl_offsets(:), ptcl_fill(:), ptcl_classes(:)
            real,    allocatable :: class_dists(:)
        end type assign_graph_ws

        type :: assign_frontier_ws
            integer, allocatable :: order(:), work_ptcl(:), sel_classes(:), sel_pos(:), likelihood_order(:)
            integer, allocatable :: ptcl_classes(:)
            real,    allocatable :: class_dist(:), work_d(:), sel_dists(:), likelihood_dists(:)
            logical, allocatable :: ptcl_avail(:)
        end type assign_frontier_ws

        type(assign_graph_ws)    :: graph
        type(assign_frontier_ws) :: frontier
        integer, allocatable :: raw_offsets(:), raw_classes(:)
        logical, allocatable :: class_active(:)
        integer :: i, k, pos, idx, icls, assigned_icls, assigned_ptcl, assigned_idx, nleft, nsel
        integer :: total_raw, total, nactive, maxclass, neligible, nassigned, last_class
        real    :: best_dist, dist_loc, huge_val
        if( .not. allocated(self%eval_touched_counts) .or. .not. allocated(self%eval_touched_refs) )&
            &THROW_HARD('simple_eul_prob_tab2D::ref_assign_sparse_likelihood; sparse bookkeeping is not allocated')
        write(logfhandle,'(A,I0,A,I0)') '>>> PROB_TAB2D_ASSIGN: preparing sparse likelihood ',&
            &self%nptcls, ' particles x ', self%nclasses
        call flush(logfhandle)
        allocate(class_active(self%nclasses))
        class_active = self%class_exists
        if( count(class_active) == 0 )then
            class_active = .true.
            THROW_WARN('No active classes after population filtering; falling back to all classes in eul_prob_tab2D')
        endif
        huge_val  = huge(1.0) / 2.0
        call build_sparse_raw_table()
        call build_sparse_assignment_graph()
        call assign_particles_from_frontier()
        call assign_remaining_particles()
        write(logfhandle,'(A,I0)') '>>> PROB_TAB2D_ASSIGN: sparse likelihood done; assigned ', nassigned
        call flush(logfhandle)
        deallocate(raw_offsets, raw_classes, class_active)

    contains

        subroutine build_sparse_raw_table()
            integer, allocatable :: raw_counts(:)
            allocate(raw_offsets(self%nptcls+1), raw_counts(self%nptcls), source=0)
            do i = 1,self%nptcls
                do k = 1,self%eval_touched_counts(i)
                    icls = self%eval_touched_refs(k,i)
                    if( .not. valid_sparse_class(i, icls) ) cycle
                    raw_counts(i) = raw_counts(i) + 1
                enddo
            enddo
            raw_offsets(1) = 1
            do i = 1,self%nptcls
                raw_offsets(i+1) = raw_offsets(i) + raw_counts(i)
            enddo
            total_raw = raw_offsets(self%nptcls+1) - 1
            if( total_raw < 1 ) THROW_HARD('empty sparse table in eul_prob_tab2D%ref_assign_sparse_likelihood')
            allocate(raw_classes(total_raw))
            raw_counts = 0
            do i = 1,self%nptcls
                do k = 1,self%eval_touched_counts(i)
                    icls = self%eval_touched_refs(k,i)
                    if( .not. valid_sparse_class(i, icls) ) cycle
                    pos = raw_offsets(i) + raw_counts(i)
                    raw_classes(pos) = icls
                    raw_counts(i)    = raw_counts(i) + 1
                enddo
            enddo
            deallocate(raw_counts)
        end subroutine build_sparse_raw_table

        subroutine build_sparse_assignment_graph()
            allocate(graph%class_counts(self%nclasses), graph%class_offsets(self%nclasses+1),&
                &graph%class_fill(self%nclasses), graph%class_pos(self%nclasses),&
                &graph%ptcl_counts(self%nptcls), graph%ptcl_offsets(self%nptcls+1), graph%ptcl_fill(self%nptcls))
            graph%class_counts = 0
            graph%ptcl_counts  = 0
            do i = 1,self%nptcls
                do pos = raw_offsets(i), raw_offsets(i+1)-1
                    icls = raw_classes(pos)
                    graph%class_counts(icls) = graph%class_counts(icls) + 1
                    graph%ptcl_counts(i)     = graph%ptcl_counts(i) + 1
                enddo
            enddo
            nactive = count(graph%class_counts > 0)
            if( nactive < 1 ) THROW_HARD('no active sparse classes in eul_prob_tab2D%ref_assign_sparse_likelihood')
            allocate(graph%active_classes(nactive))
            graph%active_classes = pack((/(icls, icls=1,self%nclasses)/), graph%class_counts > 0)
            graph%class_offsets(1) = 1
            do icls = 1,self%nclasses
                graph%class_offsets(icls+1) = graph%class_offsets(icls) + graph%class_counts(icls)
            enddo
            total = graph%class_offsets(self%nclasses+1) - 1
            allocate(graph%class_list(total), graph%class_dists(total))
            graph%ptcl_offsets(1) = 1
            do i = 1,self%nptcls
                graph%ptcl_offsets(i+1) = graph%ptcl_offsets(i) + graph%ptcl_counts(i)
            enddo
            allocate(graph%ptcl_classes(max(1,graph%ptcl_offsets(self%nptcls+1)-1)))
            graph%class_fill = graph%class_offsets(1:self%nclasses)
            graph%ptcl_fill  = graph%ptcl_offsets(1:self%nptcls)
            do i = 1,self%nptcls
                do pos = raw_offsets(i), raw_offsets(i+1)-1
                    icls = raw_classes(pos)
                    idx = graph%class_fill(icls)
                    graph%class_list(idx)  = i
                    graph%class_dists(idx) = self%loc_tab(icls,i)%dist
                    graph%class_fill(icls) = idx + 1
                    idx = graph%ptcl_fill(i)
                    graph%ptcl_classes(idx) = icls
                    graph%ptcl_fill(i)      = idx + 1
                enddo
            enddo
            maxclass = max(1, maxval(graph%class_counts))
            allocate(frontier%work_d(maxclass), frontier%order(maxclass), frontier%work_ptcl(maxclass))
            do idx = 1,nactive
                icls = graph%active_classes(idx)
                k    = graph%class_counts(icls)
                if( k <= 1 ) cycle
                pos = graph%class_offsets(icls)
                frontier%work_d(1:k)    = graph%class_dists(pos:pos+k-1)
                frontier%work_ptcl(1:k) = graph%class_list(pos:pos+k-1)
                frontier%order(1:k)     = (/(i,i=1,k)/)
                call hpsort(frontier%work_d(1:k), frontier%order(1:k))
                do i = 1,k
                    graph%class_dists(pos+i-1) = frontier%work_d(i)
                    graph%class_list(pos+i-1)  = frontier%work_ptcl(frontier%order(i))
                enddo
            enddo
            allocate(frontier%class_dist(self%nclasses), frontier%sel_pos(self%nclasses),&
                &frontier%sel_classes(nactive), frontier%sel_dists(nactive), frontier%likelihood_dists(nactive),&
                &frontier%likelihood_order(nactive), frontier%ptcl_avail(self%nptcls))
        end subroutine build_sparse_assignment_graph

        subroutine assign_particles_from_frontier()
            frontier%ptcl_avail = .true.
            nleft     = self%nptcls
            nassigned = 0
            call init_frontier()
            do while( nleft > 0 )
                if( nsel == 0 ) exit
                neligible = min(self%nhood_sz, nsel)
                if( neligible <= 1 )then
                    assigned_idx = 1
                else
                    call sample_frontier_likelihood(neligible, assigned_idx)
                endif
                call commit_selected_assignment()
            enddo
        end subroutine assign_particles_from_frontier

        subroutine init_frontier()
            graph%class_pos      = 1
            frontier%class_dist  = huge(1.0)
            frontier%sel_pos     = 0
            nsel                 = 0
            do idx = 1,nactive
                icls = graph%active_classes(idx)
                call advance_class_head(icls)
                call sync_frontier_class(icls)
            enddo
        end subroutine init_frontier

        subroutine sample_frontier_likelihood( neligible_loc, assigned_idx_loc )
            integer, intent(in)  :: neligible_loc
            integer, intent(out) :: assigned_idx_loc
            ! Top-K softmax over the current class frontier; shared implementation.
            call sample_likelihood_index(nsel, frontier%sel_dists, neligible_loc, huge_val,&
                &frontier%likelihood_dists, frontier%likelihood_order, assigned_idx_loc)
        end subroutine sample_frontier_likelihood

        subroutine commit_selected_assignment()
            assigned_icls = frontier%sel_classes(assigned_idx)
            assigned_ptcl = graph%class_list(graph%class_offsets(assigned_icls) + graph%class_pos(assigned_icls) - 1)
            frontier%ptcl_avail(assigned_ptcl) = .false.
            nleft = nleft - 1
            nassigned = nassigned + 1
            self%assgn_map(assigned_ptcl) = self%loc_tab(assigned_icls,assigned_ptcl)
            call materialize_seed_shift(self%assgn_map(assigned_ptcl), self%seed_shifts(:,assigned_ptcl),&
                &self%seed_has_sh(assigned_ptcl), self%p_ptr%l_doshift, self%seed_nrots)
            self%assgn_map(assigned_ptcl)%frac   = search_frac(assigned_ptcl)
            self%assgn_map(assigned_ptcl)%npeaks = search_npeaks(assigned_ptcl)
            call update_frontier_after_assignment(assigned_ptcl)
        end subroutine commit_selected_assignment

        subroutine update_frontier_after_assignment( iptcl_assigned )
            integer, intent(in) :: iptcl_assigned
            do idx = graph%ptcl_offsets(iptcl_assigned), graph%ptcl_offsets(iptcl_assigned+1)-1
                icls = graph%ptcl_classes(idx)
                if( graph%class_pos(icls) > graph%class_counts(icls) ) cycle
                pos = graph%class_offsets(icls)
                if( graph%class_list(pos + graph%class_pos(icls) - 1) /= iptcl_assigned ) cycle
                call advance_class_head(icls)
                call sync_frontier_class(icls)
            enddo
        end subroutine update_frontier_after_assignment

        subroutine assign_remaining_particles()
            do i = 1,self%nptcls
                if( .not. frontier%ptcl_avail(i) ) cycle
                assigned_icls = pick_best_sparse_class(i)
                if( assigned_icls == 0 )&
                    &THROW_HARD('Failed sparse likelihood particle assignment in eul_prob_tab2D%ref_assign_sparse_likelihood')
                self%assgn_map(i) = self%loc_tab(assigned_icls,i)
                call materialize_seed_shift(self%assgn_map(i), self%seed_shifts(:,i),&
                    &self%seed_has_sh(i), self%p_ptr%l_doshift, self%seed_nrots)
                self%assgn_map(i)%frac   = search_frac(i)
                self%assgn_map(i)%npeaks = search_npeaks(i)
                frontier%ptcl_avail(i) = .false.
                nassigned = nassigned + 1
            enddo
        end subroutine assign_remaining_particles

        integer function pick_best_sparse_class( iptcl_loc ) result(icls_best)
            integer, intent(in) :: iptcl_loc
            integer :: pos_loc, icls_loc
            icls_best = 0
            best_dist = huge(1.0)
            do pos_loc = raw_offsets(iptcl_loc), raw_offsets(iptcl_loc+1)-1
                icls_loc = raw_classes(pos_loc)
                dist_loc = self%loc_tab(icls_loc,iptcl_loc)%dist
                if( dist_loc < best_dist )then
                    best_dist = dist_loc
                    icls_best = icls_loc
                endif
            enddo
        end function pick_best_sparse_class

        real function search_frac( iptcl_loc ) result(frac)
            integer, intent(in) :: iptcl_loc
            frac = 100. * real(self%eval_touched_counts(iptcl_loc)) / real(self%nclasses)
        end function search_frac

        integer function search_npeaks( iptcl_loc ) result(npeaks)
            integer, intent(in) :: iptcl_loc
            npeaks = 0
            if( allocated(self%eval_touched_counts) ) npeaks = self%eval_touched_counts(iptcl_loc)
        end function search_npeaks

        subroutine advance_class_head( icls_loc )
            integer, intent(in) :: icls_loc
            integer :: nloc, start, cand
            nloc = graph%class_counts(icls_loc)
            if( graph%class_pos(icls_loc) > nloc )then
                frontier%class_dist(icls_loc) = huge(1.0)
                return
            endif
            start = graph%class_offsets(icls_loc)
            do while( graph%class_pos(icls_loc) <= nloc )
                cand = graph%class_list(start + graph%class_pos(icls_loc) - 1)
                if( frontier%ptcl_avail(cand) )then
                    frontier%class_dist(icls_loc) = graph%class_dists(start + graph%class_pos(icls_loc) - 1)
                    return
                endif
                graph%class_pos(icls_loc) = graph%class_pos(icls_loc) + 1
            enddo
            frontier%class_dist(icls_loc) = huge(1.0)
        end subroutine advance_class_head

        subroutine sync_frontier_class( icls_loc )
            integer, intent(in) :: icls_loc
            if( frontier%sel_pos(icls_loc) > 0 )then
                if( frontier%class_dist(icls_loc) >= huge(1.0) )then
                    pos = frontier%sel_pos(icls_loc)
                    if( pos < nsel )then
                        last_class = frontier%sel_classes(nsel)
                        frontier%sel_classes(pos)  = last_class
                        frontier%sel_dists(pos)    = frontier%sel_dists(nsel)
                        frontier%sel_pos(last_class)= pos
                    endif
                    frontier%sel_pos(icls_loc) = 0
                    nsel = nsel - 1
                else
                    frontier%sel_dists(frontier%sel_pos(icls_loc)) = frontier%class_dist(icls_loc)
                endif
            else
                if( frontier%class_dist(icls_loc) < huge(1.0) )then
                    nsel = nsel + 1
                    frontier%sel_classes(nsel) = icls_loc
                    frontier%sel_dists(nsel)   = frontier%class_dist(icls_loc)
                    frontier%sel_pos(icls_loc) = nsel
                endif
            endif
        end subroutine sync_frontier_class

        logical function valid_sparse_class( iptcl_loc, icls_loc ) result(l_valid)
            integer, intent(in) :: iptcl_loc, icls_loc
            l_valid = .false.
            if( icls_loc < 1 .or. icls_loc > self%nclasses ) return
            if( .not. class_active(icls_loc) ) return
            if( self%loc_tab(icls_loc,iptcl_loc)%inpl <= 0 ) return
            l_valid = .true.
        end function valid_sparse_class

    end subroutine ref_assign_sparse_likelihood

    ! ptcl -> class assignment preserving class coverage across the particle set
    subroutine ref_assign_likelihood( self )
        class(eul_prob_tab2D), intent(inout) :: self
        integer, allocatable :: stab_inds(:,:), icls_dist_inds(:), active_cls(:), eligible_cls(:), inds_sorted(:)
        real,    allocatable :: sorted_tab(:,:), cls_dists(:), corr_proxy(:), dists_raw(:,:)
        logical, allocatable :: ptcl_avail(:), class_active(:)
        integer :: i, icls, assigned_icls, assigned_ptcl
        integer :: nactive, iact, chosen_active, neligible, nhood_sz_loc, nassigned
        real    :: best_dist
        if( self%l_sparse_snhc )then
            call self%ref_assign_sparse_likelihood
            return
        endif
        write(logfhandle,'(A,I0,A,I0)') '>>> PROB_TAB2D_ASSIGN: preparing ', self%nptcls, ' particles x ', self%nclasses
        call flush(logfhandle)
        allocate(class_active(self%nclasses), active_cls(self%nclasses), eligible_cls(self%nclasses),&
            &inds_sorted(self%nclasses), icls_dist_inds(self%nclasses), cls_dists(self%nclasses), corr_proxy(self%nclasses))
        class_active = self%class_exists
        nactive      = count(class_active)
        if( nactive == 0 )then
            class_active = .true.
            nactive      = self%nclasses
            THROW_WARN('No active classes after population filtering; falling back to all classes in eul_prob_tab2D')
        endif
        nactive = 0
        do icls = 1, self%nclasses
            if( class_active(icls) )then
                nactive = nactive + 1
                active_cls(nactive) = icls
            endif
        end do
        nhood_sz_loc = min(self%nhood_sz, nactive)
        allocate(dists_raw(self%nclasses,self%nptcls), sorted_tab(self%nptcls,self%nclasses),&
            &stab_inds(self%nptcls,self%nclasses), ptcl_avail(self%nptcls))
        write(logfhandle,'(A)') '>>> PROB_TAB2D_ASSIGN: preserving raw distances'
        call flush(logfhandle)
        dists_raw = self%loc_tab(:,:)%dist
        write(logfhandle,'(A)') '>>> PROB_TAB2D_ASSIGN: using likelihood weights exp(-dist)'
        call flush(logfhandle)
        ! sort each column
        write(logfhandle,'(A)') '>>> PROB_TAB2D_ASSIGN: sorting per-class particle distances'
        call flush(logfhandle)
        sorted_tab = transpose(dists_raw)
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(icls,i)
        do icls = 1, self%nclasses
            stab_inds(:,icls) = (/(i,i=1,self%nptcls)/)
            if( class_active(icls) )then
                call hpsort(sorted_tab(:,icls), stab_inds(:,icls))
            else
                sorted_tab(:,icls) = huge(sorted_tab(1,icls))
            endif
        end do
        !$omp end parallel do
        icls_dist_inds = 1
        ptcl_avail     = .true.
        nassigned      = 0
        write(logfhandle,'(A)') '>>> PROB_TAB2D_ASSIGN: assigning particles'
        call flush(logfhandle)
        do while( any(ptcl_avail) )
            ! collect the current best distance for each class, skipping exhausted ones
            neligible = 0
            do iact = 1, nactive
                icls = active_cls(iact)
                do while( icls_dist_inds(icls) <= self%nptcls )
                    assigned_ptcl = stab_inds(icls_dist_inds(icls), icls)
                    if( ptcl_avail(assigned_ptcl) )then
                        if( (.not. self%l_sparse_snhc) .or. self%loc_tab(icls,assigned_ptcl)%inpl > 0 ) exit
                    endif
                    icls_dist_inds(icls) = icls_dist_inds(icls) + 1
                end do
                if( icls_dist_inds(icls) <= self%nptcls )then
                    neligible = neligible + 1
                    eligible_cls(neligible) = icls
                    cls_dists(neligible)    = sorted_tab(icls_dist_inds(icls), icls)
                endif
            end do
            if( neligible == 0 ) exit
            ! Top-K softmax over the eligible class frontier; shared implementation.
            call sample_likelihood_index(neligible, cls_dists, min(nhood_sz_loc, neligible), huge(1.0),&
                &corr_proxy, inds_sorted, chosen_active)
            assigned_icls = eligible_cls(chosen_active)
            assigned_ptcl = stab_inds(icls_dist_inds(assigned_icls), assigned_icls)
            ptcl_avail(assigned_ptcl)     = .false.
            self%assgn_map(assigned_ptcl) = self%loc_tab(assigned_icls, assigned_ptcl)
            self%assgn_map(assigned_ptcl)%dist = dists_raw(assigned_icls, assigned_ptcl)
            call materialize_seed_shift(self%assgn_map(assigned_ptcl), self%seed_shifts(:,assigned_ptcl),&
                &self%seed_has_sh(assigned_ptcl), self%p_ptr%l_doshift, self%seed_nrots)
            self%assgn_map(assigned_ptcl)%frac = 100.
            nassigned = nassigned + 1
        end do
        if( any(ptcl_avail) )then
            do i = 1, self%nptcls
                if( .not. ptcl_avail(i) ) cycle
                assigned_icls = 0
                best_dist     = huge(1.0)
                do iact = 1, nactive
                    icls = active_cls(iact)
                    if( self%loc_tab(icls,i)%inpl <= 0 ) cycle
                    if( self%loc_tab(icls,i)%dist < best_dist )then
                        best_dist     = self%loc_tab(icls,i)%dist
                        assigned_icls = icls
                    endif
                end do
                if( assigned_icls == 0 )then
                    do iact = 1, nactive
                        icls = active_cls(iact)
                        if( self%loc_tab(icls,i)%dist < best_dist )then
                            best_dist     = self%loc_tab(icls,i)%dist
                            assigned_icls = icls
                        endif
                    end do
                endif
                if( assigned_icls == 0 ) THROW_HARD('Failed particle assignment in eul_prob_tab2D%ref_assign_likelihood')
                self%assgn_map(i) = self%loc_tab(assigned_icls, i)
                self%assgn_map(i)%dist = dists_raw(assigned_icls, i)
                call materialize_seed_shift(self%assgn_map(i), self%seed_shifts(:,i),&
                    &self%seed_has_sh(i), self%p_ptr%l_doshift, self%seed_nrots)
                self%assgn_map(i)%frac = 100.
                ptcl_avail(i)     = .false.
                nassigned = nassigned + 1
            end do
        endif
        write(logfhandle,'(A)') '>>> PROB_TAB2D_ASSIGN: done'
        call flush(logfhandle)
        deallocate(stab_inds, icls_dist_inds, active_cls, eligible_cls, inds_sorted, sorted_tab, cls_dists, corr_proxy,&
            &dists_raw, ptcl_avail, class_active)
    end subroutine ref_assign_likelihood

    ! DESTRUCTOR

    subroutine kill( self )
        class(eul_prob_tab2D), intent(inout) :: self
        if( allocated(self%loc_tab)        ) deallocate(self%loc_tab)
        if( allocated(self%assgn_map)      ) deallocate(self%assgn_map)
        if( allocated(self%seed_shifts)    ) deallocate(self%seed_shifts)
        if( allocated(self%seed_has_sh)    ) deallocate(self%seed_has_sh)
        if( allocated(self%pinds)          ) deallocate(self%pinds)
        if( allocated(self%class_exists)   ) deallocate(self%class_exists)
        if( allocated(self%eval_touched_refs)   ) deallocate(self%eval_touched_refs)
        if( allocated(self%eval_touched_counts) ) deallocate(self%eval_touched_counts)
        if( allocated(self%prior_operated)      ) deallocate(self%prior_operated)
        self%seed_nrots = 0
        self%eval_max_touched = 0
        self%l_sparse_snhc = .false.
        self%l_prior_mode = .false.
        self%prior_kmax = 0
        self%b_ptr => null()
        self%p_ptr => null()
    end subroutine kill

    ! FILE IO

    subroutine write_tab( self, binfname )
        class(eul_prob_tab2D), intent(in) :: self
        class(string),         intent(in) :: binfname
        integer :: funit, io_stat, file_header(2)
        integer(int64) :: addr
        if( self%l_sparse_snhc )then
            call self%write_sparse_tab(binfname)
            return
        endif
        file_header(1) = self%nclasses
        file_header(2) = self%nptcls
        call fopen(funit, binfname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        write(unit=funit, pos=1) file_header
        addr = sizeof(file_header) + 1
        call write_seed_shift_table(funit, addr, self%seed_nrots, self%seed_shifts, self%seed_has_sh)
        write(funit, pos=addr) self%loc_tab
        call fclose(funit)
    end subroutine write_tab

    subroutine write_sparse_tab( self, binfname )
        class(eul_prob_tab2D), intent(in) :: self
        class(string),         intent(in) :: binfname
        type(ptcl_ref), allocatable :: sparse_tab(:)
        integer,        allocatable :: sparse_refs(:), sparse_counts(:)
        integer :: funit, io_stat, file_header(3), i, k, icls, nnz, pos
        integer(int64) :: addr
        if( .not. allocated(self%eval_touched_counts) .or. .not. allocated(self%eval_touched_refs) )&
            &THROW_HARD('eul_prob_tab2D%write_sparse_tab; sparse bookkeeping is not allocated')
        nnz = 0
        allocate(sparse_counts(self%nptcls), source=0)
        do i = 1, self%nptcls
            do k = 1, self%eval_touched_counts(i)
                icls = self%eval_touched_refs(k,i)
                if( icls < 1 .or. icls > self%nclasses ) cycle
                if( self%loc_tab(icls,i)%inpl <= 0 ) cycle
                sparse_counts(i) = sparse_counts(i) + 1
                nnz = nnz + 1
            enddo
        enddo
        if( nnz < 1 ) THROW_HARD('eul_prob_tab2D%write_sparse_tab; empty sparse table')
        allocate(sparse_refs(nnz), sparse_tab(nnz))
        pos = 0
        do i = 1, self%nptcls
            do k = 1, self%eval_touched_counts(i)
                icls = self%eval_touched_refs(k,i)
                if( icls < 1 .or. icls > self%nclasses ) cycle
                if( self%loc_tab(icls,i)%inpl <= 0 ) cycle
                pos = pos + 1
                sparse_refs(pos) = icls
                sparse_tab(pos)  = self%loc_tab(icls,i)
            enddo
        enddo
        file_header(1) = self%nclasses
        file_header(2) = self%nptcls
        file_header(3) = nnz
        call fopen(funit, binfname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        write(unit=funit, pos=1) file_header
        addr = sizeof(file_header) + 1
        write(funit, pos=addr) self%pinds
        addr = addr + sizeof(self%pinds)
        call write_seed_shift_table(funit, addr, self%seed_nrots, self%seed_shifts, self%seed_has_sh)
        write(funit, pos=addr) sparse_counts
        addr = addr + sizeof(sparse_counts)
        write(funit, pos=addr) sparse_refs
        addr = addr + sizeof(sparse_refs)
        write(funit, pos=addr) sparse_tab
        call fclose(funit)
        deallocate(sparse_tab, sparse_refs, sparse_counts)
    end subroutine write_sparse_tab

    subroutine read_tab_to_glob( self, binfname )
        class(eul_prob_tab2D), intent(inout) :: self
        class(string),         intent(in)    :: binfname
        type(ptcl_ref), allocatable :: mat_loc(:,:)
        real,           allocatable :: seed_shifts_loc(:,:)
        logical,        allocatable :: seed_has_sh_loc(:)
        integer, allocatable :: pind2glob(:), pinds_loc(:)
        integer :: funit, io_stat, file_header(2), nptcls_loc, nclasses_loc, i_loc, i_glob, pind, max_pind, seed_nrots_loc
        integer(int64) :: addr
        if( self%l_sparse_snhc )then
            call self%read_sparse_tab_to_glob(binfname)
            return
        endif
        if( file_exists(binfname) )then
            call fopen(funit, binfname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
            call fileiochk('simple_eul_prob_tab2D; read_tab_to_glob; file: '//binfname%to_char(), io_stat)
        else
            THROW_HARD('dist files of partitions should be ready!')
        endif
        read(unit=funit, pos=1) file_header
        nclasses_loc = file_header(1)
        nptcls_loc   = file_header(2)
        if( nclasses_loc .ne. self%nclasses ) THROW_HARD('nclasses mismatch in read_tab_to_glob!')
        allocate(mat_loc(nclasses_loc, nptcls_loc))
        allocate(seed_shifts_loc(2,nptcls_loc), seed_has_sh_loc(nptcls_loc))
        addr = sizeof(file_header) + 1
        call read_seed_shift_table(funit, addr, seed_nrots_loc, seed_shifts_loc, seed_has_sh_loc)
        read(unit=funit, pos=addr) mat_loc
        call fclose(funit)
        if( self%seed_nrots == 0 ) self%seed_nrots = seed_nrots_loc
        if( self%seed_nrots /= seed_nrots_loc ) THROW_HARD('seed_nrots mismatch in eul_prob_tab2D%read_tab_to_glob')
        pinds_loc = mat_loc(1,:)%pind
        call build_pind_lookup(self%pinds, pinds_loc, pind2glob, max_pind)
        if( max_pind < 1 )then
            deallocate(mat_loc, seed_shifts_loc, seed_has_sh_loc, pind2glob, pinds_loc)
            return
        endif
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i_loc,i_glob,pind)
        do i_loc = 1, nptcls_loc
            pind = mat_loc(1,i_loc)%pind
            if( pind < 1 .or. pind > max_pind ) cycle
            i_glob = pind2glob(pind)
            if( i_glob > 0 )then
                self%loc_tab(:,i_glob)     = mat_loc(:,i_loc)
                self%seed_shifts(:,i_glob) = seed_shifts_loc(:,i_loc)
                self%seed_has_sh(i_glob)   = seed_has_sh_loc(i_loc)
            endif
        end do
        !$omp end parallel do
        deallocate(mat_loc, seed_shifts_loc, seed_has_sh_loc, pind2glob, pinds_loc)
    end subroutine read_tab_to_glob

    subroutine read_sparse_tab_to_glob( self, binfname )
        class(eul_prob_tab2D), intent(inout) :: self
        class(string),         intent(in)    :: binfname
        type(ptcl_ref), allocatable :: sparse_tab(:)
        real,           allocatable :: seed_shifts_loc(:,:)
        logical,        allocatable :: seed_has_sh_loc(:)
        integer,        allocatable :: pinds_loc(:), sparse_counts(:), sparse_refs(:), pind2glob(:)
        integer :: funit, io_stat, file_header(3), nclasses_loc, nptcls_loc, nnz
        integer(int64) :: addr
        integer :: i_loc, i_glob, k, pos, icls, pind, max_pind, seed_nrots_loc
        if( file_exists(binfname) )then
            call fopen(funit, binfname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
            call fileiochk('simple_eul_prob_tab2D; read_sparse_tab_to_glob; file: '//binfname%to_char(), io_stat)
        else
            THROW_HARD('sparse dist files of partitions should be ready!')
        endif
        read(unit=funit, pos=1) file_header
        nclasses_loc = file_header(1)
        nptcls_loc   = file_header(2)
        nnz          = file_header(3)
        if( nclasses_loc .ne. self%nclasses ) THROW_HARD('nclasses mismatch in eul_prob_tab2D%read_sparse_tab_to_glob')
        if( nnz < 1 ) THROW_HARD('empty sparse table in eul_prob_tab2D%read_sparse_tab_to_glob')
        allocate(pinds_loc(nptcls_loc), seed_shifts_loc(2,nptcls_loc), seed_has_sh_loc(nptcls_loc))
        allocate(sparse_counts(nptcls_loc), sparse_refs(nnz), sparse_tab(nnz))
        addr = sizeof(file_header) + 1
        read(funit, pos=addr) pinds_loc
        addr = addr + sizeof(pinds_loc)
        call read_seed_shift_table(funit, addr, seed_nrots_loc, seed_shifts_loc, seed_has_sh_loc)
        read(funit, pos=addr) sparse_counts
        addr = addr + sizeof(sparse_counts)
        read(funit, pos=addr) sparse_refs
        addr = addr + sizeof(sparse_refs)
        read(funit, pos=addr) sparse_tab
        call fclose(funit)
        if( any(sparse_counts < 0) .or. sum(sparse_counts) /= nnz )&
            &THROW_HARD('sparse table count mismatch in eul_prob_tab2D%read_sparse_tab_to_glob')
        if( self%seed_nrots == 0 ) self%seed_nrots = seed_nrots_loc
        if( self%seed_nrots /= seed_nrots_loc ) THROW_HARD('seed_nrots mismatch in eul_prob_tab2D%read_sparse_tab_to_glob')
        if( .not. allocated(self%eval_touched_refs) )then
            self%eval_max_touched = max(1, self%nclasses)
            allocate(self%eval_touched_refs(self%eval_max_touched,self%nptcls), source=0)
        endif
        if( .not. allocated(self%eval_touched_counts) ) allocate(self%eval_touched_counts(self%nptcls), source=0)
        call build_pind_lookup(self%pinds, pinds_loc, pind2glob, max_pind)
        if( max_pind < 1 )then
            deallocate(pinds_loc, seed_shifts_loc, seed_has_sh_loc, sparse_counts, sparse_refs, sparse_tab, pind2glob)
            return
        endif
        pos = 0
        do i_loc = 1, nptcls_loc
            pind = pinds_loc(i_loc)
            i_glob = 0
            if( pind >= 1 .and. pind <= max_pind ) i_glob = pind2glob(pind)
            if( i_glob > 0 )then
                self%seed_shifts(:,i_glob) = seed_shifts_loc(:,i_loc)
                self%seed_has_sh(i_glob)   = seed_has_sh_loc(i_loc)
            endif
            do k = 1, sparse_counts(i_loc)
                pos = pos + 1
                if( i_glob < 1 ) cycle
                icls = sparse_refs(pos)
                if( icls < 1 .or. icls > self%nclasses ) cycle
                self%loc_tab(icls,i_glob) = sparse_tab(pos)
                call self%mark_ref_touched(i_glob, icls)
            enddo
        enddo
        if( pos /= nnz ) THROW_HARD('sparse table count mismatch in eul_prob_tab2D%read_sparse_tab_to_glob')
        deallocate(pinds_loc, seed_shifts_loc, seed_has_sh_loc, sparse_counts, sparse_refs, sparse_tab, pind2glob)
    end subroutine read_sparse_tab_to_glob

    subroutine write_assignment( self, binfname )
        class(eul_prob_tab2D), intent(in) :: self
        class(string),         intent(in) :: binfname
        integer :: funit, io_stat, headsz
        headsz = sizeof(self%nptcls)
        call fopen(funit, binfname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        write(unit=funit, pos=1)          self%nptcls
        write(unit=funit, pos=headsz + 1) self%assgn_map
        call fclose(funit)
    end subroutine write_assignment

    subroutine write_prior_topk( self, binfname )
        class(eul_prob_tab2D), intent(in) :: self
        class(string),         intent(in) :: binfname
        integer, allocatable :: rec_ints(:)
        integer, allocatable :: cand_cls(:), top_ind(:)
        real,    allocatable :: cand_dist(:)
        real    :: dist_thres
        integer :: funit, io_stat, reclen
        integer :: i, k, icls, ncand, pick, k_use, k_this, npeaks_detected
        if( self%nptcls < 1 .or. self%nclasses < 1 ) return
        k_use = max(1, nint(PRIOR2D_TOPK_FRAC * real(self%nclasses)))
        allocate(rec_ints(k_use), source=0)
        allocate(cand_cls(self%nclasses), cand_dist(self%nclasses))
        allocate(top_ind(k_use))
        inquire(iolength=reclen) rec_ints
        open(newunit=funit, file=binfname%to_char(), status='REPLACE', action='WRITE', &
            &access='DIRECT', form='UNFORMATTED', recl=reclen, iostat=io_stat)
        call fileiochk('simple_eul_prob_tab2D; write_prior_topk; file: '//binfname%to_char(), io_stat)
        do i = 1, self%nptcls
            if( self%pinds(i) < 1 ) cycle
            ncand = 0
            if( self%l_sparse_snhc .and. allocated(self%eval_touched_counts) )then
                do k = 1, self%eval_touched_counts(i)
                    icls = self%eval_touched_refs(k,i)
                    if( icls < 1 .or. icls > self%nclasses ) cycle
                    if( self%loc_tab(icls,i)%inpl <= 0 ) cycle
                    ncand = ncand + 1
                    cand_cls(ncand)  = icls
                    cand_dist(ncand) = self%loc_tab(icls,i)%dist
                enddo
            endif
            if( ncand < 1 )then
                do icls = 1, self%nclasses
                    ncand = ncand + 1
                    cand_cls(ncand)  = icls
                    cand_dist(ncand) = self%loc_tab(icls,i)%dist
                enddo
            endif
            ! Fixed-width direct-access records are retained for simplicity:
            ! k_use is the 0.3*ncls upper bound, while the FDR threshold selects
            ! a dynamic lower-tail objective neighborhood per particle.  Unused
            ! record slots stay zero and are ignored by apply_prior_order.
            call detect_peak_thres_fdr(ncand, cand_dist(1:ncand), PRIOR2D_FDR_Q, 1, min(k_use,ncand),&
                &dist_thres, npeaks_detected, lower_tail=.true.)
            k_this = min(k_use, ncand, max(1, npeaks_detected))
            rec_ints = 0
            do k = 1, k_this
                pick = minloc(cand_dist(1:ncand), dim=1)
                top_ind(k) = pick
                cand_dist(pick) = huge(1.0)
            enddo
            do k = 1, k_this
                rec_ints(k) = cand_cls(top_ind(k))
            enddo
            write(unit=funit, rec=self%pinds(i), iostat=io_stat) rec_ints
            call fileiochk('simple_eul_prob_tab2D; write_prior_topk(rec); file: '//binfname%to_char(), io_stat)
        enddo
        close(funit)
        deallocate(rec_ints, cand_cls, cand_dist, top_ind)
    end subroutine write_prior_topk

    subroutine read_assignment( self, binfname )
        class(eul_prob_tab2D), intent(inout) :: self
        class(string),         intent(in)    :: binfname
        type(ptcl_ref), allocatable :: assgn_glob(:)
        integer, allocatable :: pind2glob(:), pinds_glob(:)
        integer :: funit, io_stat, nptcls_glob, headsz, i_loc, i_glob, pind, max_pind
        headsz = sizeof(nptcls_glob)
        if( .not. file_exists(binfname) )then
            THROW_HARD('file '//binfname%to_char()//' does not exist!')
        else
            call fopen(funit, binfname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
        end if
        call fileiochk('simple_eul_prob_tab2D; read_assignment; file: '//binfname%to_char(), io_stat)
        read(unit=funit, pos=1) nptcls_glob
        allocate(assgn_glob(nptcls_glob))
        read(unit=funit, pos=headsz + 1) assgn_glob
        call fclose(funit)
        pinds_glob = assgn_glob(:)%pind
        call build_pind_lookup(pinds_glob, self%pinds, pind2glob, max_pind)
        if( max_pind < 1 )then
            deallocate(assgn_glob, pind2glob, pinds_glob)
            return
        endif
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i_loc,i_glob,pind)
        do i_loc = 1, self%nptcls
            pind = self%pinds(i_loc)
            if( pind < 1 .or. pind > max_pind ) cycle
            i_glob = pind2glob(pind)
            if( i_glob > 0 ) self%assgn_map(i_loc) = assgn_glob(i_glob)
        end do
        !$omp end parallel do
        deallocate(assgn_glob, pind2glob, pinds_glob)
    end subroutine read_assignment

    subroutine init_eval2D_sparse_ws( self, nclasses, nrots, smpl_ncls, nthr )
        class(eval2D_sparse_ws), intent(inout) :: self
        integer,                 intent(in)    :: nclasses, nrots, smpl_ncls, nthr
        call self%dealloc_eval2D_sparse_ws
        allocate(self%direct_srch_order(nclasses,nthr), self%eval_cls(nclasses,nthr), source=0)
        allocate(self%best_locs(max(1,smpl_ncls),nthr), source=0)
        allocate(self%vec_nrots(nrots,nthr), source=0)
        allocate(self%inpl_corrs(nrots,nthr), self%eval_dists(nclasses,nthr), source=huge(1.0))
    end subroutine init_eval2D_sparse_ws

    subroutine dealloc_eval2D_sparse_ws( self )
        class(eval2D_sparse_ws), intent(inout) :: self
        if( allocated(self%direct_srch_order) ) deallocate(self%direct_srch_order)
        if( allocated(self%vec_nrots)         ) deallocate(self%vec_nrots)
        if( allocated(self%eval_cls)          ) deallocate(self%eval_cls)
        if( allocated(self%best_locs)         ) deallocate(self%best_locs)
        if( allocated(self%inpl_corrs)        ) deallocate(self%inpl_corrs)
        if( allocated(self%eval_dists)        ) deallocate(self%eval_dists)
    end subroutine dealloc_eval2D_sparse_ws

end module simple_eul_prob_tab2D
