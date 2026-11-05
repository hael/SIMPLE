!@descr: 2D probability table routines for multi-reference class assignment by per-particle draws from the tempered posterior
module simple_eul_prob_tab2D
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use, intrinsic :: iso_fortran_env, only: int64
use simple_pftc_srch_api
use simple_builder,            only: builder
use simple_pftc_shsrch_grad,   only: pftc_shsrch_grad, shift_ptcl2ref, shift_ref2ptcl
use simple_polarft_calc, only: SHERRSQ, SHIFTMAP_SCREEN, SHIFTMAP_WINDOW, NSHIFTMAPS
use simple_type_defs,          only: OBJFUN_EUCLID
use simple_srch_sort_loc,      only: minnloc
use simple_decay_funs,         only: extremal_decay2D
use simple_eul_prob_tab_utils, only: update_screen_best, eulprob_dist_switch, record_seed_shift,&
    &prob_candidate,&
    &prob_candidate_buffer, prob_candidate_store, prob_fill_profile, prob_dense_entry, prob_refined_entry,&
    &posterior_temperature, draw_posterior_row, prob_neff_readout, PROB_NEFF_ROWS, PROB_NEFF_NMIN,&
    &screen_pool_size, shift_needs_disk, disk_radius_px, prob_rnd_scope, keyed_candidate_order, shift_prior_sigma2,&
    &prob_readout, prob_readout_fname, write_prob_readout, report_prob_readout, PROB_READOUT_2D_DENSE, PROB_READOUT_2D_SPARSE
use simple_rnd,                only: rnd_key, rnd_keyed, rnd_draw_key, RND_FAM_ASSIGN
implicit none

public :: eul_prob_tab2D
private
#include "simple_local_flags.inc"

real,             parameter :: NHOOD_FRAC           = 0.1
type :: eul_prob_tab2D
    class(builder),    pointer  :: b_ptr => null()
    class(parameters), pointer  :: p_ptr => null()
    type(prob_candidate_store)        :: candidate_store !< sparse: the batch's candidates at the assignment
    type(prob_candidate_buffer), allocatable :: candidate_buffers(:) !< per thread: the fill's candidates of the batch
    type(ptcl_ref), allocatable :: assgn_map(:)      !< assignment map (nptcls)
    real,           allocatable :: seed_shifts(:,:)  !< per-particle seeded shift (2,nptcls)
    logical,        allocatable :: seed_has_sh(:)    !< per-particle seeded shift flag
    real,           allocatable :: temps(:)          !< per-particle posterior temperature pi wsqsum / (2 pftsz) (nptcls; the worker's fill writes it, the master reads it and draws at prob_beta times it)
    integer                     :: seed_nrots = 0    !< rotation grid of the table (checked when the parts are read)
    integer,        allocatable :: pinds(:)          !< particle indices for processing
    logical,        allocatable :: class_exists(:)   !< class population filter
    logical                     :: l_sparse_snhc = .false.
    integer                     :: nptcls            !< size of pinds array
    integer                     :: nclasses          !< number of classes
    integer                     :: nhood_sz = 1      !< classes around the best one whose shifts fill_tab_range refines
    type(prob_fill_profile)     :: fill_profile      !< cost record of the fills since new (bench file)
    type(rnd_key)               :: rnd_scope         !< the keyed draws' scope of the iteration (prob_rnd_scope), taken at every fill and assignment
    real                        :: sigma2_prior_model = 0. !< the iteration's shift prior variance per axis (shift_prior_sigma2), set by the fill's caller
    logical                     :: l_model_set = .false.   !< set_iteration_model called
    ! the batch's rows, kept in the thread buffers by the fill for the batch's own selection
    type(prob_dense_entry),   allocatable :: batch_rows(:,:)    !< dense: (nclasses, i_first:i_last) the batch's rows
    type(prob_refined_entry), allocatable :: batch_refined(:,:) !< dense: (maxref, i_first:i_last) their refined entries
    integer,                  allocatable :: batch_nrefined(:)  !< dense: (i_first:i_last)
    ! the assignment across the batches (open_assignment .. close_assignment): the readout and the counts
    type(prob_neff_readout)     :: neff
    integer,        allocatable :: active_cls(:)    !< the active classes in class order
    integer                     :: nactive = 0
    logical                     :: l_assign_open = .false.
    integer                     :: nassigned = 0, nmode = 0, ntemp_missing = 0
    contains
    ! CONSTRUCTOR
    procedure :: new_pipeline
    ! MAIN PROCEDURES
    procedure :: set_iteration_model
    procedure :: set_model_from_field
    procedure :: fill_tab_range
    procedure :: get_fill_profile
    procedure :: open_assignment
    procedure :: assign_batch
    procedure :: close_assignment
    procedure, private :: assign_rows_dense
    procedure, private :: assign_rows_sparse
    procedure, private :: batch_rows_from_buffers
    procedure, private :: batch_store_from_buffers
    procedure, private :: row_candidate
    ! DESTRUCTOR
    procedure :: kill
    ! PRIVATE
    procedure, private :: fill_tab_prob_neigh_range
    procedure, private :: record_sparse_eval
    procedure, private :: new_common
    procedure, private :: initialize_storage
    procedure, private :: finalize_candidate_assignment
end type eul_prob_tab2D

type :: eval2D_sparse_ws
    integer, allocatable :: direct_srch_order(:,:) ! [nclasses,nthr]
    integer, allocatable :: eval_cls(:,:)          ! [nclasses,nthr]
    integer, allocatable :: best_locs(:,:)         ! [smpl_ncls,nthr]
    real,    allocatable :: eval_dists(:,:)        ! [nclasses,nthr]
contains
    procedure :: init_eval2D_sparse_ws
    procedure :: dealloc_eval2D_sparse_ws
end type eval2D_sparse_ws

contains

    ! CONSTRUCTORS

    subroutine new_common( self, params, build, pinds )
        class(eul_prob_tab2D),     intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        class(builder),    target, intent(in)    :: build
        integer,                   intent(in)    :: pinds(:)
        integer :: nactive
        call self%kill
        self%p_ptr     => params
        self%b_ptr     => build
        self%nptcls    = size(pinds)
        self%nclasses  = params%ncls
        self%l_sparse_snhc = trim(params%refine) == 'prob_neigh'
        ! the posterior is defined by the Euclidean objective (the particle temperature is its noise scale)
        if( params%cc_objfun /= OBJFUN_EUCLID )&
            &THROW_HARD('the 2D probability tables draw from the tempered posterior: refine=prob|prob_neigh needs objfun=euclid')
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
        allocate(self%seed_shifts(2,self%nptcls), source=0.)
        allocate(self%seed_has_sh(self%nptcls), source=.false.)
        allocate(self%temps(self%nptcls), source=0.)
    end subroutine new_common

    subroutine initialize_storage( self )
        class(eul_prob_tab2D), intent(inout) :: self
        integer :: i, icls, iptcl
        real    :: x
        !$omp parallel do default(shared) private(i,iptcl,icls) proc_bind(close) schedule(static)
        do i = 1, self%nptcls
            iptcl = self%pinds(i)
            if( allocated(self%assgn_map) )then
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
            endif
        end do
        !$omp end parallel do
    end subroutine initialize_storage

    !> the pipeline's table of a part (the shared-memory matcher): the assignment map and the thread buffers the fills
    !! keep for the batch's own selection (assign_batch)
    subroutine new_pipeline( self, params, build, pinds )
        class(eul_prob_tab2D),     intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        class(builder),    target, intent(in)    :: build
        integer,                   intent(in)    :: pinds(:)
        call self%new_common(params,build,pinds)
        allocate(self%assgn_map(self%nptcls))
        allocate(self%candidate_buffers(nthr_glob))
        call self%initialize_storage
    end subroutine new_pipeline

    !> the iteration's model from the complete incoming field: the shift prior's width (shift_prior=yes)
    subroutine set_model_from_field( self )
        class(eul_prob_tab2D), intent(inout) :: self
        real :: sigma2_prior
        sigma2_prior = 0.
        if( self%p_ptr%l_shift_prior ) sigma2_prior = shift_prior_sigma2(self%b_ptr%spproj_field, self%p_ptr%lp / self%p_ptr%smpd)
        call self%set_iteration_model(sigma2_prior)
    end subroutine set_model_from_field

    !> The iteration's model the fills read instead of the particle field: the shift prior's variance per
    !! axis (shift_prior_sigma2), computed once per iteration by the fill's caller from the complete incoming
    !! field, so that a batch's rows never depend on what earlier batches committed
    subroutine set_iteration_model( self, sigma2_prior )
        class(eul_prob_tab2D), intent(inout) :: self
        real,                  intent(in)    :: sigma2_prior
        self%sigma2_prior_model = sigma2_prior
        self%l_model_set        = .true.
    end subroutine set_iteration_model

    subroutine fill_tab_range( self, i_first, i_last )
        class(eul_prob_tab2D), intent(inout) :: self
        integer,               intent(in)    :: i_first, i_last
        type(pftc_shsrch_grad), target  :: grad_shsrch_obj(nthr_glob)  !< shift search object, L-BFGS with gradient
        type(pftc_shsrch_grad), target  :: grad_shsrch_disk(nthr_glob) !< the disk's (shift_sigma > 0): its limits bound the polish box
        type(pftc_shsrch_grad), pointer :: gso => null()
        type(ori)              :: o_prev
        integer :: i, icls, iptcl, ithr, irot, irot0, icls_prev, iref_n, nactive, nhood_sz_loc, npool, imap
        integer, allocatable :: pool(:)
        integer :: i_from, i_to
        integer :: active_cls(self%nclasses)
        integer :: loc(1)
        real    :: lims(2,2), lims_init(2,2), cxy(3), cxy_prob(3), xy0(2)
        real    :: cls_dists(self%nclasses), cls_dists_work(self%nclasses)
        real    :: inpl_dist, sigma2_prior, step_map(NSHIFTMAPS), shsig_loc
        logical :: class_active(self%nclasses)
        real(dp) :: ts_seed(nthr_glob), ts_score(nthr_glob), ts_refine(nthr_glob), t_region, t0
        integer(int64) :: nscored(nthr_glob), nbank(nthr_glob), nrefined(nthr_glob), ndisk(nthr_glob), nscreens(nthr_glob)
        complex(sp), allocatable :: seed_bank(:,:,:) !< per thread: the seed-shifted particle memo
        real(dp) :: prior_w, prior_w_map, temp_i, temp_s
        real     :: loss_cell, scr_xy(2,self%nclasses), scr_sdev(self%nclasses), width_best(self%nclasses)
        integer  :: irho, ialpha, kscreen, scr_irot(self%nclasses)
        logical  :: l_bank, l_map, l_disk, l_polished
        if( self%l_sparse_snhc )then
            call self%fill_tab_prob_neigh_range(i_first, i_last)
            return
        endif
        i_from = max(1, i_first)
        i_to   = min(self%nptcls, i_last)
        if( i_to < i_from ) return
        if( .not. self%l_model_set ) THROW_HARD('the iteration model is set by the caller before the fill (set_iteration_model); eul_prob_tab2D%fill_tab_range')
        ! the keyed draws of this fill: the in-plane draw of every scored class on the calculator
        self%rnd_scope = prob_rnd_scope(self%p_ptr%which_iter, .true.)
        call self%b_ptr%pftc%set_rnd_scope(self%rnd_scope)
        self%seed_nrots = self%b_ptr%pftc%get_nrots()
        ts_seed   = 0.d0
        ts_score  = 0.d0
        ts_refine = 0.d0
        nscored   = 0_int64
        nbank     = 0_int64
        nrefined  = 0_int64
        ndisk     = 0_int64
        nscreens  = 0_int64
        sigma2_prior = 0.
        l_map      = .false.
        prior_w    = 0.d0
        kscreen    = 0
        step_map   = 0.
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
        t_region = omp_get_wtime()
        if( self%p_ptr%l_doshift )then
            ! make shift search objects
            lims(:,1)      = -self%p_ptr%trs
            lims(:,2)      =  self%p_ptr%trs
            lims_init(:,1) = -SHC_INPL_TRSHWDTH
            lims_init(:,2) =  SHC_INPL_TRSHWDTH
            do ithr = 1, nthr_glob
                call grad_shsrch_obj(ithr)%new_alternating(self%b_ptr, lims, lims_init=lims_init, &
                    &maxits=self%p_ptr%maxits_sh, coarse_init=.true.)
            end do
            ! the shift prior's width of the iteration (set per particle below)
            if( self%p_ptr%l_shift_prior ) sigma2_prior = self%sigma2_prior_model
            ! the polar shift map (shift_sigma > 0): the coarse map over a particle's box at its likely classes, the
            ! polish at full resolution from its optimum. The box is the disk of three sigma under the prior for a
            ! particle not yet centred in this workflow or whose prior width exceeds a third of the window (the
            ! dynamic centring: shift_needs_disk), the window of trs for the others
            l_map = self%p_ptr%l_centring
            if( l_map )then
                call setup_shiftmap2D(self, prior_w, kscreen, step_map)
                call new_disk_objects2D(self, grad_shsrch_disk, lims_init)
            endif
            call alloc_seed_bank2D(self, seed_bank)
            ! fill the table; dynamic: the particles' shift searches differ in cost, and each particle is processed
            ! whole by one thread with its own search object; its draws are keyed on the particle
            !$omp parallel do default(shared) private(i,iptcl,ithr,irot,irot0,cxy,icls,icls_prev,loc,&
            !$omp& cls_dists,cls_dists_work,iref_n,cxy_prob,inpl_dist,t0,l_bank,xy0,temp_i,irho,ialpha,&
            !$omp& loss_cell,temp_s,npool,pool,gso,l_disk,scr_xy,scr_irot,scr_sdev,width_best,imap,prior_w_map,l_polished,shsig_loc)&
            !$omp firstprivate(o_prev) proc_bind(close) schedule(dynamic)
            do i = i_from, i_to
                iptcl = self%pinds(i)
                ithr  = omp_get_thread_num() + 1
                t0    = omp_get_wtime()
                self%temps(i) = posterior_temperature(self%b_ptr%pftc%wsqsum_ptcl_now(iptcl), self%b_ptr%pftc%get_pftsz())
                call self%b_ptr%spproj_field%get_ori(iptcl, o_prev)      ! previous ori
                icls_prev = nint(self%b_ptr%spproj_field%get(iptcl, 'class'))
                cxy    = 0.
                l_disk = .false.
                gso    => grad_shsrch_obj(ithr)
                if( l_map )then
                    ! (1) the seed is the stored shift the sections already carry; the box is the disk for a particle not
                    ! yet centred in this workflow (under the map's prior) or whose prior width exceeds a third of the
                    ! window (the data's alone), the window for the others
                    shsig_loc   = self%b_ptr%spproj_field%get(iptcl, 'shsig')
                    l_disk      = shift_needs_disk(shsig_loc, self%p_ptr%trs)
                    imap        = SHIFTMAP_WINDOW
                    prior_w_map = 0.d0
                    if( l_disk )then
                        gso         => grad_shsrch_disk(ithr)
                        imap        = SHIFTMAP_SCREEN
                        if( shsig_loc < 0. ) prior_w_map = prior_w
                        ndisk(ithr) = ndisk(ithr) + 1_int64
                    endif
                    temp_i = real(self%p_ptr%prob_beta, dp) * real(self%temps(i), dp)
                    temp_s = real(self%p_ptr%prob_beta, dp) * real(posterior_temperature(&
                        &self%b_ptr%pftc%wsqsum_ptcl_band(iptcl, kscreen), self%b_ptr%pftc%get_pftsz()), dp)
                else
                    ! (1) identify shifts using the previously assigned best class
                    if( self%p_ptr%l_shift_prior ) call grad_shsrch_obj(ithr)%set_shift_prior(sigma2_prior,&
                        &posterior_temperature(self%b_ptr%pftc%wsqsum_ptcl_now(iptcl), self%b_ptr%pftc%get_pftsz()))
                    irot0 = self%b_ptr%pftc%get_roind(360. - o_prev%e3get()) ! in-plane angle index seed
                    if( icls_prev >= 1 .and. icls_prev <= self%nclasses )then
                        if( class_active(icls_prev) )then
                            irot = irot0
                            call grad_shsrch_obj(ithr)%set_indices(icls_prev, iptcl)
                            cxy = grad_shsrch_obj(ithr)%minimize(irot=irot, sh_rot=.false.)
                            if( irot == 0 )then
                                cxy(2:3) = 0.
                            else
                                ! one particle-frame shift from here on: the estimate converted at its rotation
                                cxy(2:3) = shift_ref2ptcl(cxy(2:3), self%b_ptr%pftc%get_rot(irot))
                            endif
                        endif
                    endif
                endif
                self%seed_shifts(:,i) = cxy(2:3)
                self%seed_has_sh(i)   = .true.
                call memoize_seed_bank2D(self, iptcl, cxy(2:3), seed_bank(:,:,ithr), l_bank, nbank(ithr))
                ts_seed(ithr) = ts_seed(ithr) + (omp_get_wtime() - t0)
                ! (2) search class references using that shared shift
                t0 = omp_get_wtime()
                cls_dists = huge(1.0)
                do iref_n = 1, nactive
                    icls = active_cls(iref_n)
                    call score_seeded_class(self, icls, iptcl, cxy(2:3), seed_bank(:,:,ithr), l_bank, inpl_dist, irot)
                    call self%record_sparse_eval(i, ithr, icls, inpl_dist, irot, 0., 0., .false.)
                    cls_dists(icls) = inpl_dist
                end do
                ts_score(ithr) = ts_score(ithr) + (omp_get_wtime() - t0)
                ! (3) the shifts of a neighborhood of classes around the best one. With the map: the screen, the coarse
                ! map at the origin's best classes holding SCREEN_QUANTILE of the particle's posterior there (a pick far
                ! off leaves it flat and keeps every class), at least twice the refined count, each class keeping its
                ! coarse optimum; then the polish at full resolution from that optimum. Else the windowed refinement
                t0 = omp_get_wtime()
                cls_dists_work = cls_dists
                if( l_map )then
                    pool  = minnloc(cls_dists, nactive)
                    npool = min(nactive, max(2 * nhood_sz_loc, screen_pool_size(cls_dists(pool), nactive, temp_i)))
                    cls_dists_work = huge(1.0)
                    width_best = huge(1.)
                    do iref_n = 1, npool
                        icls = pool(iref_n)
                        if( icls < 1 ) cycle
                        call self%b_ptr%pftc%gen_shiftmap(imap, icls, iptcl)
                        call self%b_ptr%pftc%shiftmap_mode(imap, temp_s, prior_w_map, &
                            &irot, irho, ialpha, loss_cell, &
                            &temp_phys=temp_s / real(self%p_ptr%prob_beta, dp), sdev=scr_sdev(icls), &
                            &width_cutoff=maxval(width_best(1:nhood_sz_loc)))
                        call update_screen_best(loss_cell, width_best(1:nhood_sz_loc))
                        scr_irot(icls) = irot
                        scr_xy(:,icls) = real(self%b_ptr%pftc%shiftmap_shift(imap, irho, ialpha))
                        cls_dists_work(icls) = loss_cell
                        nscreens(ithr) = nscreens(ithr) + 1_int64
                    enddo
                    deallocate(pool)
                endif
                do iref_n = 1, nhood_sz_loc
                    loc = minloc(cls_dists_work)
                    icls = loc(1)
                    if( cls_dists_work(icls) >= huge(1.0)/2.0 ) exit
                    nrefined(ithr) = nrefined(ithr) + 1_int64
                    call gso%set_indices(icls, iptcl)
                    if( l_map )then
                        ! the polish at full resolution in a box of one coarse step about the coarse optimum: the box's
                        ! optimum if it beat the origin, minimize's baseline; else the coarse optimum, scored at full
                        ! resolution, if it beats the origin's row; else that row stands
                        irot = scr_irot(icls)
                        call gso%minimize_in_cell(irot, scr_xy(:,icls), step_map(imap), cxy_prob, l_polished)
                        inpl_dist = eulprob_dist_switch(cxy_prob(1), self%p_ptr%cc_objfun)
                        if( l_polished .or. inpl_dist < cls_dists(icls) )then
                            call self%record_sparse_eval(i, ithr, icls, inpl_dist, irot, cxy_prob(2), cxy_prob(3), .true.,&
                                &scr_sdev(icls))
                        endif
                    else
                        irot = self%candidate_buffers(ithr)%get_inpl(i, icls)
                        xy0  = seed_start(self, cxy(2:3), irot)
                        cxy_prob = gso%minimize(irot=irot, sh_rot=.true., xy_in=xy0)
                        if( irot > 0 )then
                            call self%record_sparse_eval(i, ithr, icls,&
                                &eulprob_dist_switch(cxy_prob(1), self%p_ptr%cc_objfun), irot,&
                                &cxy_prob(2), cxy_prob(3), .true.)
                        endif
                    endif
                    cls_dists_work(icls) = huge(1.0)
                end do
                ts_refine(ithr) = ts_refine(ithr) + (omp_get_wtime() - t0)
            end do
            !$omp end parallel do
        else
            ! shift-free path: evaluate classes/in-planes at zero shift
            !$omp parallel do default(shared) private(i,iptcl,ithr,icls,iref_n,irot,inpl_dist,t0)&
            !$omp& proc_bind(close) schedule(dynamic)
            do i = i_from, i_to
                iptcl = self%pinds(i)
                ithr  = omp_get_thread_num() + 1
                t0    = omp_get_wtime()
                self%temps(i) = posterior_temperature(self%b_ptr%pftc%wsqsum_ptcl_now(iptcl), self%b_ptr%pftc%get_pftsz())
                do iref_n = 1, nactive
                    icls = active_cls(iref_n)
                    call self%b_ptr%pftc%gen_posterior_val(icls, iptcl, [0.,0.], self%b_ptr%pftc%ptcl_posterior_temp(iptcl),&
                        &inpl_dist, irot)
                    call self%record_sparse_eval(i, ithr, icls, inpl_dist, irot, 0., 0., .false.)
                end do
                ts_score(ithr) = ts_score(ithr) + (omp_get_wtime() - t0)
            end do
            !$omp end parallel do
        endif
        call record_fill_profile2D(self, grad_shsrch_obj, omp_get_wtime() - t_region, ts_seed, ts_score, ts_refine,&
            &nscored, nbank, nrefined, ndisk, nscreens, i_to - i_from + 1, nactive, nhood_sz_loc, grad_disk=grad_shsrch_disk)
        do ithr = 1, nthr_glob
            call grad_shsrch_obj(ithr)%kill
            if( l_map ) call grad_shsrch_disk(ithr)%kill
        end do
        if( allocated(seed_bank) ) deallocate(seed_bank)
        call o_prev%kill
    end subroutine fill_tab_range

    subroutine fill_tab_prob_neigh_range( self, i_first, i_last )
        class(eul_prob_tab2D), intent(inout) :: self
        integer,               intent(in)    :: i_first, i_last
        type(pftc_shsrch_grad), target :: grad_shsrch_obj(nthr_glob)
        type(pftc_shsrch_grad), target :: grad_shsrch_disk(nthr_glob) !< the disk's (shift_sigma > 0)
        type(eval2D_sparse_ws) :: eval_work
        integer :: ithr, i, iptcl, nrefs_bound, smpl_ncls, nrots
        integer :: i_from, i_to
        real    :: lims(2,2), lims_init(2,2), neigh_frac, cxy(3), sigma2_prior
        real(dp) :: ts_seed(nthr_glob), ts_score(nthr_glob), ts_refine(nthr_glob), t_region
        integer(int64) :: nscored(nthr_glob), nbank(nthr_glob), nrefined(nthr_glob), ndisk(nthr_glob), nscreens(nthr_glob)
        complex(sp), allocatable :: seed_bank(:,:,:) !< per thread: the seed-shifted particle memo
        logical  :: l_bank(nthr_glob), l_map
        real(dp) :: prior_w
        real     :: step_map(NSHIFTMAPS)
        integer  :: kscreen
        i_from = max(1, i_first)
        i_to   = min(self%nptcls, i_last)
        if( i_to < i_from ) return
        l_map     = .false.
        prior_w   = 0.d0
        kscreen   = 0
        step_map  = 0.
        ts_seed   = 0.d0
        ts_score  = 0.d0
        ts_refine = 0.d0
        nscored   = 0_int64
        nbank     = 0_int64
        nrefined  = 0_int64
        ndisk     = 0_int64
        nscreens  = 0_int64
        sigma2_prior = 0.
        if( .not. self%l_model_set ) THROW_HARD('the iteration model is set by the caller before the fill (set_iteration_model); eul_prob_tab2D%fill_tab_prob_neigh_range')
        ! the keyed draws of this fill: the candidate order per particle and the in-plane draw per scored class
        self%rnd_scope = prob_rnd_scope(self%p_ptr%which_iter, .true.)
        call self%b_ptr%pftc%set_rnd_scope(self%rnd_scope)
        nrots = self%b_ptr%pftc%get_nrots()
        self%seed_nrots = nrots
        neigh_frac  = extremal_decay2D(self%p_ptr%extr_iter, self%p_ptr%extr_lim)
        nrefs_bound = max(2, min(self%nclasses, nint(real(self%nclasses) * (1. - neigh_frac))))
        nrefs_bound = max(1, min(nrefs_bound, self%nclasses))
        smpl_ncls   = neighfrac2nsmpl(neigh_frac, self%nclasses)
        smpl_ncls   = max(1, min(smpl_ncls, self%nclasses))
        call eval_work%init_eval2D_sparse_ws(self%nclasses, smpl_ncls, nthr_glob)
        l_bank   = .false.
        call alloc_seed_bank2D(self, seed_bank)
        t_region = omp_get_wtime()
        if( self%p_ptr%l_doshift )then
            lims(:,1)      = -self%p_ptr%trs
            lims(:,2)      =  self%p_ptr%trs
            lims_init(:,1) = -SHC_INPL_TRSHWDTH
            lims_init(:,2) =  SHC_INPL_TRSHWDTH
            do ithr = 1, nthr_glob
                call grad_shsrch_obj(ithr)%new_alternating(self%b_ptr, lims, lims_init=lims_init, &
                    &maxits=self%p_ptr%maxits_sh, coarse_init=.true.)
            end do
            ! the shift prior's width of the iteration (set per particle below)
            if( self%p_ptr%l_shift_prior ) sigma2_prior = self%sigma2_prior_model
            ! the polar shift map (shift_sigma > 0): the classes are scored at the stored shift, without a seed, and
            ! the coarse map initialises the best classes' polish (minimize_in_cell): the window map for a centred
            ! particle, the disk map under the prior over every class for one not yet centred in this workflow or
            ! whose prior width exceeds a third of the window (the dynamic centring: shift_needs_disk)
            l_map = self%p_ptr%l_centring
            if( l_map )then
                call setup_shiftmap2D(self, prior_w, kscreen, step_map)
                call new_disk_objects2D(self, grad_shsrch_disk, lims_init)
            endif
            ! dynamic: the particles' shift searches differ in cost; each is processed whole by one thread with its
            ! own search object; its draws are keyed on the particle
            !$omp parallel do default(shared) private(i,iptcl,ithr,cxy) proc_bind(close) schedule(dynamic)
            do i = i_from, i_to
                iptcl = self%pinds(i)
                ithr  = omp_get_thread_num() + 1
                call process_particle(i, iptcl, ithr, .true., cxy)
                self%seed_shifts(:,i) = cxy(2:3)
                self%seed_has_sh(i)   = .true.
            enddo
            !$omp end parallel do
        else
            !$omp parallel do default(shared) private(i,iptcl,ithr,cxy) proc_bind(close) schedule(dynamic)
            do i = i_from, i_to
                iptcl = self%pinds(i)
                ithr  = omp_get_thread_num() + 1
                cxy   = 0.
                call process_particle(i, iptcl, ithr, .false., cxy)
            enddo
            !$omp end parallel do
        endif
        call record_fill_profile2D(self, grad_shsrch_obj, omp_get_wtime() - t_region, ts_seed, ts_score, ts_refine,&
            &nscored, nbank, nrefined, ndisk, nscreens, i_to - i_from + 1, 0, smpl_ncls, grad_disk=grad_shsrch_disk)
        do ithr = 1, nthr_glob
            call grad_shsrch_obj(ithr)%kill
            if( l_map ) call grad_shsrch_disk(ithr)%kill
        end do
        call eval_work%dealloc_eval2D_sparse_ws
        deallocate(seed_bank)

    contains

        subroutine process_particle(i_loc, iptcl_loc, ithr_loc, l_with_shift, cxy_loc)
            integer, intent(in)    :: i_loc, iptcl_loc, ithr_loc
            logical, intent(in)    :: l_with_shift
            real,    intent(inout) :: cxy_loc(3)
            type(ori) :: o_prev_loc
            integer :: icls_prev_loc, irot0_loc, irot_loc, isample_loc, icls_loc, neval_loc
            integer :: neval_bound_loc
            real    :: dist_loc, sh_loc(2)
            real(dp) :: t0
            logical  :: l_disk_loc, l_prior_loc
            real     :: shsig_loc
            t0 = omp_get_wtime()
            call self%b_ptr%spproj_field%get_ori(iptcl_loc, o_prev_loc)
            icls_prev_loc = nint(self%b_ptr%spproj_field%get(iptcl_loc, 'class'))
            ! the candidate order: the previous class last, the others keyed on the particle
            call keyed_candidate_order(self%rnd_scope, iptcl_loc, self%nclasses, icls_prev_loc, eval_work%direct_srch_order(:,ithr_loc))
            cxy_loc    = 0.
            l_disk_loc  = .false.
            l_prior_loc = .false.
            if( l_map .and. l_with_shift )then
                shsig_loc   = self%b_ptr%spproj_field%get(iptcl_loc, 'shsig')
                l_disk_loc  = shift_needs_disk(shsig_loc, self%p_ptr%trs)
                l_prior_loc = shsig_loc < 0.
            endif
            self%temps(i_loc) = posterior_temperature(self%b_ptr%pftc%wsqsum_ptcl_now(iptcl_loc), self%b_ptr%pftc%get_pftsz())
            if( self%p_ptr%l_shift_prior .and. l_with_shift ) call grad_shsrch_obj(ithr_loc)%set_shift_prior(sigma2_prior,&
                &posterior_temperature(self%b_ptr%pftc%wsqsum_ptcl_now(iptcl_loc), self%b_ptr%pftc%get_pftsz()))
            if( l_with_shift .and. .not. l_map .and. icls_prev_loc >= 1 .and. icls_prev_loc <= self%nclasses )then
                irot0_loc = self%b_ptr%pftc%get_roind(360. - o_prev_loc%e3get())
                irot_loc  = irot0_loc
                call grad_shsrch_obj(ithr_loc)%set_indices(icls_prev_loc, iptcl_loc)
                cxy_loc = grad_shsrch_obj(ithr_loc)%minimize(irot=irot_loc, sh_rot=.false.)
                if( irot_loc == 0 )then
                    cxy_loc(2:3) = 0.
                else
                    ! one particle-frame shift from here on: the estimate converted at its rotation
                    cxy_loc(2:3) = shift_ref2ptcl(cxy_loc(2:3), self%b_ptr%pftc%get_rot(irot_loc))
                endif
            endif
            sh_loc = 0.
            if( l_with_shift ) sh_loc = cxy_loc(2:3)
            call memoize_seed_bank2D(self, iptcl_loc, sh_loc, seed_bank(:,:,ithr_loc), l_bank(ithr_loc),&
                &nbank(ithr_loc))
            ts_seed(ithr_loc) = ts_seed(ithr_loc) + (omp_get_wtime() - t0)
            t0 = omp_get_wtime()
            neval_bound_loc = nrefs_bound
            if( l_disk_loc ) neval_bound_loc = self%nclasses ! a particle on the disk scores every class
            eval_work%eval_cls(:,ithr_loc)   = 0
            eval_work%eval_dists(:,ithr_loc) = huge(1.0)
            neval_loc = 0
            do isample_loc = 1, self%nclasses
                if( isample_loc > neval_bound_loc ) exit
                icls_loc = eval_work%direct_srch_order(isample_loc,ithr_loc)
                call score_class(icls_loc, ithr_loc, iptcl_loc, sh_loc, dist_loc, irot_loc)
                call self%record_sparse_eval(i_loc, ithr_loc, icls_loc, dist_loc, irot_loc, 0., 0., .false.)
                neval_loc = neval_loc + 1
                eval_work%eval_cls(neval_loc,ithr_loc)   = icls_loc
                eval_work%eval_dists(neval_loc,ithr_loc) = dist_loc
            enddo
            ts_score(ithr_loc) = ts_score(ithr_loc) + (omp_get_wtime() - t0)
            nscored(ithr_loc)  = nscored(ithr_loc) + int(neval_loc, int64)
            t0 = omp_get_wtime()
            call refine_best_classes(i_loc, ithr_loc, iptcl_loc, sh_loc, neval_loc, l_with_shift, l_disk_loc, l_prior_loc)
            ts_refine(ithr_loc) = ts_refine(ithr_loc) + (omp_get_wtime() - t0)
            call o_prev_loc%kill
        end subroutine process_particle

        !> one class scored at the seed
        subroutine score_class(icls_loc, ithr_loc, iptcl_loc, sh_loc, dist_loc, irot_loc)
            integer,  intent(in)  :: icls_loc, ithr_loc, iptcl_loc
            real,     intent(in)  :: sh_loc(2)
            real,     intent(out) :: dist_loc
            integer,  intent(out) :: irot_loc
            call score_seeded_class(self, icls_loc, iptcl_loc, sh_loc, seed_bank(:,:,ithr_loc), l_bank(ithr_loc),&
                &dist_loc, irot_loc)
            if( irot_loc < 1 .or. irot_loc > nrots ) irot_loc = 1
        end subroutine score_class

        !> the shifts of the best classes: the windowed refinement from the seed, or under the key the polish in a box
        !! about the coarse map's optimum (kept if it beat the origin, or the cell if it beats the class's row): the
        !! window map at the best classes; for a particle on the disk, the screen of the classes holding
        !! SCREEN_QUANTILE of its posterior at the origin with the disk map under the prior, the best of the screen polished
        subroutine refine_best_classes(i_loc, ithr_loc, iptcl_loc, sh_loc, neval_loc, l_with_shift, l_disk_loc, l_prior_loc)
            integer,  intent(in) :: i_loc, ithr_loc, iptcl_loc, neval_loc
            real,     intent(in) :: sh_loc(2)
            logical,  intent(in) :: l_with_shift, l_disk_loc, l_prior_loc
            type(pftc_shsrch_grad), pointer :: gso
            integer  :: nrefine_loc, j_loc, eval_slot_loc, icls_loc, irot_loc, irho_loc, ialpha_loc, imap_loc, npool_loc
            integer  :: pool_loc(neval_loc), scr_irot_loc(self%nclasses), slot_of_loc(self%nclasses), loc_loc(1)
            real     :: cxy_prob_loc(3), xy0_loc(2), loss_cell_loc, dist_loc, sdev_loc, dists_sorted_loc(neval_loc)
            real     :: scr_xy_loc(2,self%nclasses), scr_loss_loc(self%nclasses), scr_sdev_loc(self%nclasses)
            real :: width_best_loc(max(1,smpl_ncls))
            real(dp) :: temp_s_loc, temp_i_loc, prior_w_loc
            logical  :: l_polished_loc
            if( .not. l_with_shift ) return
            if( neval_loc < 1 ) return
            gso => grad_shsrch_obj(ithr_loc)
            if( l_map )then
                temp_s_loc  = real(self%p_ptr%prob_beta, dp) * real(posterior_temperature(&
                    &self%b_ptr%pftc%wsqsum_ptcl_band(iptcl_loc, kscreen), self%b_ptr%pftc%get_pftsz()), dp)
                imap_loc    = SHIFTMAP_WINDOW
                prior_w_loc = 0.d0
                if( l_disk_loc )then
                    imap_loc    = SHIFTMAP_SCREEN
                    if( l_prior_loc ) prior_w_loc = prior_w
                    gso         => grad_shsrch_disk(ithr_loc)
                    ndisk(ithr_loc) = ndisk(ithr_loc) + 1_int64
                endif
            endif
            nrefine_loc = min(smpl_ncls, neval_loc)
            if( l_map .and. l_disk_loc )then
                ! the screen: the disk map at the classes holding SCREEN_QUANTILE of the posterior at the origin (at
                ! least twice the refined count), then the best of the screen polished in a box about their optimum
                pool_loc         = minnloc(eval_work%eval_dists(1:neval_loc,ithr_loc), neval_loc)
                dists_sorted_loc = eval_work%eval_dists(pool_loc,ithr_loc)
                temp_i_loc       = real(self%p_ptr%prob_beta, dp) * real(self%temps(i_loc), dp)
                npool_loc        = min(neval_loc, max(2 * nrefine_loc, screen_pool_size(dists_sorted_loc, neval_loc, temp_i_loc)))
                scr_loss_loc = huge(1.0)
                width_best_loc = huge(1.)
                slot_of_loc  = 0
                do j_loc = 1, npool_loc
                    eval_slot_loc = pool_loc(j_loc)
                    if( eval_slot_loc < 1 ) cycle
                    icls_loc = eval_work%eval_cls(eval_slot_loc,ithr_loc)
                    if( icls_loc < 1 ) cycle
                    slot_of_loc(icls_loc) = eval_slot_loc
                    call self%b_ptr%pftc%gen_shiftmap(imap_loc, icls_loc, iptcl_loc)
                    call self%b_ptr%pftc%shiftmap_mode(imap_loc, temp_s_loc, prior_w_loc, &
                        &irot_loc, irho_loc, ialpha_loc, loss_cell_loc, &
                        &temp_phys=temp_s_loc / real(self%p_ptr%prob_beta, dp), sdev=scr_sdev_loc(icls_loc), &
                        &width_cutoff=maxval(width_best_loc(1:nrefine_loc)))
                    call update_screen_best(loss_cell_loc, width_best_loc(1:nrefine_loc))
                    scr_irot_loc(icls_loc) = irot_loc
                    scr_xy_loc(:,icls_loc) = real(self%b_ptr%pftc%shiftmap_shift(imap_loc, irho_loc, ialpha_loc))
                    scr_loss_loc(icls_loc) = loss_cell_loc
                    nscreens(ithr_loc) = nscreens(ithr_loc) + 1_int64
                enddo
                do j_loc = 1, nrefine_loc
                    loc_loc  = minloc(scr_loss_loc)
                    icls_loc = loc_loc(1)
                    if( scr_loss_loc(icls_loc) >= huge(1.0)/2. ) exit
                    nrefined(ithr_loc) = nrefined(ithr_loc) + 1_int64
                    call gso%set_indices(icls_loc, iptcl_loc)
                    irot_loc = scr_irot_loc(icls_loc)
                    call gso%minimize_in_cell(irot_loc, scr_xy_loc(:,icls_loc), step_map(imap_loc), cxy_prob_loc, l_polished_loc)
                    dist_loc = eulprob_dist_switch(cxy_prob_loc(1), self%p_ptr%cc_objfun)
                    if( l_polished_loc .or. dist_loc < eval_work%eval_dists(slot_of_loc(icls_loc),ithr_loc) )then
                        call self%record_sparse_eval(i_loc, ithr_loc, icls_loc, dist_loc, irot_loc,&
                            &cxy_prob_loc(2), cxy_prob_loc(3), .true., scr_sdev_loc(icls_loc))
                    endif
                    scr_loss_loc(icls_loc) = huge(1.0)
                enddo
                return
            endif
            eval_work%best_locs(1:nrefine_loc,ithr_loc) = minnloc(eval_work%eval_dists(1:neval_loc,ithr_loc), nrefine_loc)
            do j_loc = 1, nrefine_loc
                eval_slot_loc = eval_work%best_locs(j_loc,ithr_loc)
                if( eval_slot_loc < 1 ) cycle
                icls_loc = eval_work%eval_cls(eval_slot_loc,ithr_loc)
                if( icls_loc < 1 ) cycle
                nrefined(ithr_loc) = nrefined(ithr_loc) + 1_int64
                call gso%set_indices(icls_loc, iptcl_loc)
                if( l_map )then
                    call self%b_ptr%pftc%gen_shiftmap(SHIFTMAP_WINDOW, icls_loc, iptcl_loc)
                    call self%b_ptr%pftc%shiftmap_mode(SHIFTMAP_WINDOW, temp_s_loc, 0.d0, &
                        &irot_loc, irho_loc, ialpha_loc, loss_cell_loc, &
                        &temp_phys=temp_s_loc / real(self%p_ptr%prob_beta, dp), sdev=sdev_loc)
                    xy0_loc = real(self%b_ptr%pftc%shiftmap_shift(SHIFTMAP_WINDOW, irho_loc, ialpha_loc))
                    nscreens(ithr_loc) = nscreens(ithr_loc) + 1_int64
                    call gso%minimize_in_cell(irot_loc, xy0_loc, step_map(SHIFTMAP_WINDOW), cxy_prob_loc, l_polished_loc)
                    dist_loc = eulprob_dist_switch(cxy_prob_loc(1), self%p_ptr%cc_objfun)
                    if( l_polished_loc .or. dist_loc < eval_work%eval_dists(eval_slot_loc,ithr_loc) )then
                        call self%record_sparse_eval(i_loc, ithr_loc, icls_loc, dist_loc, irot_loc,&
                            &cxy_prob_loc(2), cxy_prob_loc(3), .true., sdev_loc)
                    endif
                else
                    irot_loc = self%candidate_buffers(ithr_loc)%get_inpl(i_loc, icls_loc)
                    xy0_loc  = seed_start(self, sh_loc, irot_loc)
                    cxy_prob_loc = gso%minimize(irot=irot_loc, sh_rot=.true., xy_in=xy0_loc)
                    if( irot_loc > 0 )then
                        call self%record_sparse_eval(i_loc, ithr_loc, icls_loc,&
                            &eulprob_dist_switch(cxy_prob_loc(1), self%p_ptr%cc_objfun),&
                            &irot_loc, cxy_prob_loc(2), cxy_prob_loc(3), .true.)
                    endif
                endif
            enddo
        end subroutine refine_best_classes

    end subroutine fill_tab_prob_neigh_range

    !> The polar shift maps of a fill (shift_sigma > 0, the centring stage), built on the calculator before the fill's
    !! parallel region: the disk of three sigma about the stored shift bounded by the search range trs, in the pixels
    !! of the shift matrices (the native box), magnitudes a quarter of the matching wavelength apart; prior_w is the
    !! prior's weight 1 / (2 sigma^2). The fill uses the coarse map (SHIFTMAP_SCREEN) alone
    !> the coarse maps of a fill (shift_sigma > 0): the disk of three sigma (disk_radius_px) and the window of trs,
    !! with each map's radial step, the box of the polish about a cell
    subroutine setup_shiftmap2D( self, prior_w, kscreen, step_map )
        class(eul_prob_tab2D), intent(inout) :: self
        real(dp),              intent(out)   :: prior_w
        integer,               intent(out)   :: kscreen
        real,                  intent(out)   :: step_map(NSHIFTMAPS)
        real(dp) :: sigma_px
        if( .not. (self%p_ptr%prob_beta > 0.) ) THROW_HARD('the shift map draws at prob_beta x T: prob_beta must be positive')
        sigma_px = real(self%p_ptr%shift_sigma, dp) / real(self%p_ptr%smpd, dp)
        call self%b_ptr%pftc%setup_shiftmaps(sigma_px, real(disk_radius_px(self%p_ptr%shift_sigma, self%p_ptr%smpd,&
            &self%p_ptr%box), dp), real(self%p_ptr%trs, dp), prior_w, kscreen)
        step_map = 0.
        step_map(SHIFTMAP_SCREEN) = real(norm2(self%b_ptr%pftc%shiftmap_shift(SHIFTMAP_SCREEN, 2, 1) -&
            &self%b_ptr%pftc%shiftmap_shift(SHIFTMAP_SCREEN, 1, 1)))
        step_map(SHIFTMAP_WINDOW) = real(norm2(self%b_ptr%pftc%shiftmap_shift(SHIFTMAP_WINDOW, 2, 1) -&
            &self%b_ptr%pftc%shiftmap_shift(SHIFTMAP_WINDOW, 1, 1)))
    end subroutine setup_shiftmap2D

    !> the disk's own search objects (shift_sigma > 0): the polish box clamps to the object's limits
    subroutine new_disk_objects2D( self, objs, lims_init )
        class(eul_prob_tab2D),  intent(in)    :: self
        type(pftc_shsrch_grad), intent(inout) :: objs(:)
        real,                   intent(in)    :: lims_init(2,2)
        real    :: lims_disk(2,2)
        integer :: ithr
        lims_disk(:,1) = -disk_radius_px(self%p_ptr%shift_sigma, self%p_ptr%smpd, self%p_ptr%box)
        lims_disk(:,2) = -lims_disk(:,1)
        do ithr = 1, size(objs)
            call objs(ithr)%new_alternating(self%b_ptr, lims_disk, lims_init=lims_init, maxits=self%p_ptr%maxits_sh,&
                &coarse_init=.true.)
        end do
    end subroutine new_disk_objects2D

    !> the per-thread memo of the particle shifted by its particle-frame seed, which every class of the particle is
    !! scored against (a seed below SHERRSQ scores against the plain particle memo)
    subroutine alloc_seed_bank2D( self, seed_bank )
        class(eul_prob_tab2D),    intent(in)    :: self
        complex(sp), allocatable, intent(inout) :: seed_bank(:,:,:)
        integer :: kfromto(2)
        kfromto = self%b_ptr%pftc%get_kfromto()
        allocate(seed_bank(self%b_ptr%pftc%get_pftsz()+1, kfromto(2)-kfromto(1)+1, nthr_glob))
    end subroutine alloc_seed_bank2D

    !> shifts the particle once by its seed (l_bank: the bank holds it) and counts the memos built
    subroutine memoize_seed_bank2D( self, iptcl, shift_xy, bank, l_bank, nbank )
        class(eul_prob_tab2D), intent(in)    :: self
        integer,               intent(in)    :: iptcl
        real,                  intent(in)    :: shift_xy(2)
        complex(sp),           intent(inout) :: bank(:,:)
        logical,               intent(out)   :: l_bank
        integer(int64),        intent(inout) :: nbank
        l_bank = sum(shift_xy**2) > SHERRSQ
        if( .not. l_bank ) return
        call self%b_ptr%pftc%memoize_shifted_ptcl(iptcl, shift_xy, bank)
        nbank = nbank + 1_int64
    end subroutine memoize_seed_bank2D

    !> one class scored against the particle shifted by shift_xy (particle frame): from the seed bank when it holds
    !! the shifted particle, with the angular Nyquist bin, otherwise by the shifted-reference kernel. The in-plane
    !! rotation is drawn from the row's tempered posterior at the particle's temperature over every grid rotation
    subroutine score_seeded_class( self, icls, iptcl, shift_xy, bank, l_bank, dist, irot )
        class(eul_prob_tab2D), intent(in)  :: self
        integer,               intent(in)  :: icls, iptcl
        real,                  intent(in)  :: shift_xy(2)
        complex(sp),           intent(in)  :: bank(:,:)
        logical,               intent(in)  :: l_bank
        real,                  intent(out) :: dist
        integer,               intent(out) :: irot
        if( l_bank )then
            call self%b_ptr%pftc%gen_posterior_val_bank(icls, iptcl, bank, self%b_ptr%pftc%ptcl_posterior_temp(iptcl),&
                &dist, irot)
        else
            call self%b_ptr%pftc%gen_posterior_val(icls, iptcl, shift_xy, self%b_ptr%pftc%ptcl_posterior_temp(iptcl),&
                &dist, irot)
        endif
    end subroutine score_seeded_class

    !> the start of a class's shift refinement: the particle-frame seed as the shift of the class at the rotation it was
    !! scored at, clamped into the search box
    function seed_start( self, shift_xy, irot ) result( xy )
        class(eul_prob_tab2D), intent(in) :: self
        real,                  intent(in) :: shift_xy(2)
        integer,               intent(in) :: irot
        real :: xy(2)
        xy = shift_xy
        if( irot > 0 ) xy = shift_ptcl2ref(shift_xy, self%b_ptr%pftc%get_rot(irot))
        xy = max(-self%p_ptr%trs, min(self%p_ptr%trs, xy))
    end function seed_start

    !> the cost record of the fills since new, read before kill
    function get_fill_profile( self ) result( profile )
        class(eul_prob_tab2D), intent(in) :: self
        type(prob_fill_profile) :: profile
        profile = self%fill_profile
        if( associated(self%b_ptr) )then
            call self%b_ptr%pftc%get_shiftmap_profile(profile%map_times, profile%map_counts, profile%map_geometry)
        endif
    end function get_fill_profile

    !> adds one fill to the cost record (bench file): the wall time of the parallel region, thread-seconds of the seed,
    !! the score loop and the refinement, classes scored (nrefs per particle in the plain fill, nscored in the
    !! prob_neigh fill), the most classes refined per particle and their total, shifted particle memos built and the
    !! search objects' counters
    subroutine record_fill_profile2D( self, grad_obj, wall, ts_seed, ts_score, ts_refine, nscored, nbank, nrefined, &
            &ndisk, nscreens, nfilled, nrefs, nrefine , grad_disk)
        class(eul_prob_tab2D),  intent(inout) :: self
        type(pftc_shsrch_grad), intent(inout) :: grad_obj(:)
        type(pftc_shsrch_grad), optional, intent(inout) :: grad_disk(:)
        real(dp),               intent(in)    :: wall, ts_seed(:), ts_score(:), ts_refine(:)
        integer(int64),         intent(in)    :: nscored(:), nbank(:), nrefined(:), ndisk(:), nscreens(:)
        integer,                intent(in)    :: nfilled, nrefs, nrefine
        integer(int64) :: obj_ithr, grad_ithr, alt_ithr, grid_ithr, min_ithr
        integer        :: ithr
        self%fill_profile%wall_region = self%fill_profile%wall_region + wall
        self%fill_profile%ts_seed     = self%fill_profile%ts_seed     + sum(ts_seed)
        self%fill_profile%ts_score    = self%fill_profile%ts_score    + sum(ts_score)
        self%fill_profile%ts_refine   = self%fill_profile%ts_refine   + sum(ts_refine)
        self%fill_profile%nthr        = nthr_glob
        self%fill_profile%nrefs       = nrefs
        self%fill_profile%nptcls      = self%fill_profile%nptcls      + int(nfilled, int64)
        self%fill_profile%nscored     = self%fill_profile%nscored     + sum(nscored)
        self%fill_profile%bank_builds = self%fill_profile%bank_builds + sum(nbank)
        self%fill_profile%nrefined    = self%fill_profile%nrefined    + sum(nrefined)
        self%fill_profile%disk_ptcls  = self%fill_profile%disk_ptcls  + sum(ndisk)
        self%fill_profile%screens     = self%fill_profile%screens     + sum(nscreens)
        if( .not. self%p_ptr%l_doshift ) return
        self%fill_profile%projs_ns = max(self%fill_profile%projs_ns, nrefine)
        do ithr = 1,size(grad_obj)
            call grad_obj(ithr)%get_profile(obj_ithr, grad_ithr, alt_ithr, grid_ithr, min_ithr)
            call add_counters
        enddo
        if( present(grad_disk) )then
            do ithr = 1,size(grad_disk)
                call grad_disk(ithr)%get_profile(obj_ithr, grad_ithr, alt_ithr, grid_ithr, min_ithr)
                call add_counters
            enddo
        endif

        contains

            subroutine add_counters
                self%fill_profile%objective_evals = self%fill_profile%objective_evals + obj_ithr
                self%fill_profile%gradient_evals  = self%fill_profile%gradient_evals  + grad_ithr
                self%fill_profile%alternations    = self%fill_profile%alternations    + alt_ithr
                self%fill_profile%grid_points     = self%fill_profile%grid_points     + grid_ithr
                self%fill_profile%minimize_calls  = self%fill_profile%minimize_calls  + min_ithr
            end subroutine add_counters

    end subroutine record_fill_profile2D

    subroutine record_sparse_eval( self, i, ithr, icls, dist, irot, x, y, has_sh, sdev )
        class(eul_prob_tab2D), intent(inout) :: self
        integer,               intent(in)    :: i, ithr, icls, irot
        real,                  intent(in)    :: dist, x, y
        logical,               intent(in)    :: has_sh
        real,        optional, intent(in)    :: sdev
        type(prob_candidate) :: candidate
        if( icls < 1 .or. icls > self%nclasses )&
            &THROW_HARD('simple_eul_prob_tab2D::record_sparse_eval; class index out of range')
        candidate%iref   = icls
        candidate%dist   = dist
        candidate%inpl   = irot
        candidate%x      = x
        candidate%y      = y
        candidate%has_sh = has_sh
        if( present(sdev) ) candidate%sdev = sdev
        call self%candidate_buffers(ithr)%append_or_replace(i, candidate)
    end subroutine record_sparse_eval


    !> the entry of class icls of particle i's row, and the shift of its refined entry if it has one; from the batch
    !! the fill kept
    function row_candidate( self, icls, i ) result( candidate )
        class(eul_prob_tab2D), intent(in) :: self
        integer,               intent(in) :: icls, i
        type(prob_candidate) :: candidate
        integer :: k
        candidate%iref = icls
        candidate%inpl = self%batch_rows(icls,i)%inpl
        candidate%dist = self%batch_rows(icls,i)%dist
        do k = 1,self%batch_nrefined(i)
            if( self%batch_refined(k,i)%ri /= icls ) cycle
            candidate%x      = self%batch_refined(k,i)%x
            candidate%y      = self%batch_refined(k,i)%y
            candidate%sdev   = self%batch_refined(k,i)%sdev
            candidate%has_sh = .true.
            exit
        enddo
    end function row_candidate

    subroutine finalize_candidate_assignment( self, iptcl_loc, candidate, frac, npeaks )
        class(eul_prob_tab2D), intent(inout) :: self
        integer,               intent(in)    :: iptcl_loc, npeaks
        type(prob_candidate),  intent(in)    :: candidate
        real,                  intent(in)    :: frac
        self%assgn_map(iptcl_loc)%pind   = self%pinds(iptcl_loc)
        self%assgn_map(iptcl_loc)%iproj  = 0
        self%assgn_map(iptcl_loc)%istate = 0
        self%assgn_map(iptcl_loc)%icls   = candidate%iref
        self%assgn_map(iptcl_loc)%inpl   = candidate%inpl
        self%assgn_map(iptcl_loc)%dist   = candidate%dist
        self%assgn_map(iptcl_loc)%x      = candidate%x
        self%assgn_map(iptcl_loc)%y      = candidate%y
        self%assgn_map(iptcl_loc)%sdev   = candidate%sdev
        self%assgn_map(iptcl_loc)%has_sh = candidate%has_sh
        call record_seed_shift(self%assgn_map(iptcl_loc),self%seed_shifts(:,iptcl_loc),&
            &self%seed_has_sh(iptcl_loc),self%p_ptr%l_doshift)
        self%assgn_map(iptcl_loc)%frac   = frac
        self%assgn_map(iptcl_loc)%npeaks = npeaks
    end subroutine finalize_candidate_assignment

    !> the assignment opened (the batches' assign_batch run between it and
    !! close_assignment): the active classes (every class: a class must be able to recover), the n_eff readout over the rows of the global sample (nglob its size,
    !! rank0 the rank before this table's first row: a part's offset), the draws' scope, the counts; by default the
    !! table's own rows are the sample
    subroutine open_assignment( self, nglob, rank0 )
        class(eul_prob_tab2D), intent(inout) :: self
        integer, optional,     intent(in)    :: nglob, rank0
        logical, allocatable :: class_active(:)
        integer :: icls, nglob_here, rank0_here
        if( self%l_assign_open ) THROW_HARD('the assignment is open; eul_prob_tab2D%open_assignment')
        nglob_here = self%nptcls
        rank0_here = 0
        if( present(nglob) ) nglob_here = nglob
        if( present(rank0) ) rank0_here = rank0
        self%rnd_scope = prob_rnd_scope(self%p_ptr%which_iter, .true.)
        allocate(class_active(self%nclasses))
        class_active = self%class_exists
        if( count(class_active) == 0 )then
            class_active = .true.
            THROW_WARN('No active classes after population filtering; falling back to all classes in eul_prob_tab2D')
        endif
        if( allocated(self%active_cls) ) deallocate(self%active_cls)
        allocate(self%active_cls(self%nclasses), source=0)
        self%nactive = 0
        do icls = 1, self%nclasses
            if( class_active(icls) )then
                self%nactive = self%nactive + 1
                self%active_cls(self%nactive) = icls
            endif
        end do
        deallocate(class_active)
        ! the n_eff readout: the mean effective number of candidates over a subsample of the rows the draws are from
        call self%neff%new(nglob_here, PROB_NEFF_ROWS, PROB_NEFF_NMIN, rank0=rank0_here)
        self%nassigned     = 0
        self%nmode         = 0
        self%ntemp_missing = 0
        self%l_assign_open = .true.
    end subroutine open_assignment

    !> the batch's selection (the pipeline): the candidates the fill kept in its thread buffers become the batch's rows
    !! (dense) or its candidate store (sparse), in their append order; assigned and released
    subroutine assign_batch( self, i_first, i_last )
        class(eul_prob_tab2D), intent(inout) :: self
        integer,               intent(in)    :: i_first, i_last
        if( self%l_sparse_snhc )then
            call self%batch_store_from_buffers(i_first, i_last)
            call self%assign_rows_sparse(i_first, i_last)
            call self%candidate_store%kill
        else
            call self%batch_rows_from_buffers(i_first, i_last)
            call self%assign_rows_dense(i_first, i_last)
            deallocate(self%batch_rows, self%batch_refined, self%batch_nrefined)
        endif
    end subroutine assign_batch

    !> the assignment closed: its readouts as a record, published by a distributed part (the master merges the parts'
    !! records in part order and reports) or reported here by the project owner
    subroutine close_assignment( self )
        class(eul_prob_tab2D), intent(inout) :: self
        type(prob_readout) :: rec
        if( .not. self%l_assign_open ) THROW_HARD('no open assignment; eul_prob_tab2D%close_assignment')
        rec%route         = merge(PROB_READOUT_2D_SPARSE, PROB_READOUT_2D_DENSE, self%l_sparse_snhc)
        rec%nptcls        = self%nassigned
        rec%nrefs         = self%nclasses
        rec%nmode         = self%nmode
        rec%ntemp_missing = self%ntemp_missing
        rec%neff_nrows    = self%neff%nkept()
        rec%neff_sum_draw = self%neff%sum_neff(self%p_ptr%prob_beta)
        rec%neff_sum_phys = self%neff%sum_neff(1.)
        call self%neff%kill
        if( self%p_ptr%l_distr_worker )then
            call write_prob_readout(rec, prob_readout_fname(self%p_ptr%part, self%p_ptr%numlen))
        else
            call report_prob_readout(rec, self%p_ptr%prob_beta)
        endif
        if( allocated(self%active_cls) ) deallocate(self%active_cls)
        self%nactive       = 0
        self%l_assign_open = .false.
    end subroutine close_assignment

    !> the thread buffers of a dense batch as its rows and refined entries: every class of a particle's row once (its
    !! thread's buffer holds the particle whole), the refined entries those with a shift
    subroutine batch_rows_from_buffers( self, i_first, i_last )
        class(eul_prob_tab2D), intent(inout) :: self
        integer,               intent(in)    :: i_first, i_last
        type(prob_candidate) :: candidate
        logical, allocatable :: filled(:,:)
        integer :: ithr, j, i_loc, icls, maxref
        if( allocated(self%batch_rows) ) THROW_HARD('the previous batch is not assigned yet; eul_prob_tab2D%batch_rows_from_buffers')
        allocate(self%batch_nrefined(i_first:i_last), source=0)
        do ithr = 1,size(self%candidate_buffers)
            do j = 1,self%candidate_buffers(ithr)%nused
                i_loc = self%candidate_buffers(ithr)%particle_indices(j)
                if( i_loc < i_first .or. i_loc > i_last ) THROW_HARD('a candidate of another batch in the buffers; eul_prob_tab2D%batch_rows_from_buffers')
                if( self%candidate_buffers(ithr)%candidates(j)%has_sh ) self%batch_nrefined(i_loc) = self%batch_nrefined(i_loc) + 1
            enddo
        enddo
        maxref = maxval(self%batch_nrefined)
        allocate(self%batch_rows(self%nclasses,i_first:i_last), self%batch_refined(max(1,maxref),i_first:i_last))
        allocate(filled(self%nclasses,i_first:i_last), source=.false.)
        self%batch_nrefined = 0
        do ithr = 1,size(self%candidate_buffers)
            do j = 1,self%candidate_buffers(ithr)%nused
                i_loc     = self%candidate_buffers(ithr)%particle_indices(j)
                candidate = self%candidate_buffers(ithr)%candidates(j)
                icls      = candidate%iref
                if( icls < 1 .or. icls > self%nclasses ) THROW_HARD('class out of range; eul_prob_tab2D%batch_rows_from_buffers')
                if( filled(icls,i_loc) ) THROW_HARD('class written twice in a dense row; eul_prob_tab2D%batch_rows_from_buffers')
                self%batch_rows(icls,i_loc)%dist = candidate%dist
                self%batch_rows(icls,i_loc)%inpl = candidate%inpl
                filled(icls,i_loc) = .true.
                if( candidate%has_sh )then
                    self%batch_nrefined(i_loc) = self%batch_nrefined(i_loc) + 1
                    self%batch_refined(self%batch_nrefined(i_loc),i_loc) = &
                        &prob_refined_entry(i_loc, icls, candidate%x, candidate%y, candidate%dist, candidate%sdev)
                endif
            enddo
            call self%candidate_buffers(ithr)%kill
        enddo
        if( .not. all(filled) ) THROW_HARD('incomplete dense row; eul_prob_tab2D%batch_rows_from_buffers')
        deallocate(filled)
    end subroutine batch_rows_from_buffers

    !> the thread buffers of a sparse batch as the candidate store of the part: the particles outside the batch hold
    !! no candidate; a particle's candidates come in its thread's append order
    subroutine batch_store_from_buffers( self, i_first, i_last )
        class(eul_prob_tab2D), intent(inout) :: self
        integer,               intent(in)    :: i_first, i_last
        integer, allocatable :: counts(:), filled(:)
        integer :: ithr, j, i_loc
        integer(int64) :: pos
        allocate(counts(self%nptcls), filled(self%nptcls), source=0)
        do ithr = 1,size(self%candidate_buffers)
            do j = 1,self%candidate_buffers(ithr)%nused
                i_loc = self%candidate_buffers(ithr)%particle_indices(j)
                if( i_loc < i_first .or. i_loc > i_last ) THROW_HARD('a candidate of another batch in the buffers; eul_prob_tab2D%batch_store_from_buffers')
                counts(i_loc) = counts(i_loc) + 1
            enddo
        enddo
        call self%candidate_store%new_ragged(self%pinds, counts)
        do ithr = 1,size(self%candidate_buffers)
            do j = 1,self%candidate_buffers(ithr)%nused
                i_loc = self%candidate_buffers(ithr)%particle_indices(j)
                filled(i_loc) = filled(i_loc) + 1
                pos = self%candidate_store%offsets(i_loc) + int(filled(i_loc) - 1,int64)
                self%candidate_store%candidates(pos) = self%candidate_buffers(ithr)%candidates(j)
            enddo
            call self%candidate_buffers(ithr)%kill
        enddo
        do i_loc = i_first,i_last
            self%candidate_store%seed_shifts(:,i_loc) = self%seed_shifts(:,i_loc)
            self%candidate_store%seed_has_sh(i_loc)   = self%seed_has_sh(i_loc)
        enddo
        deallocate(counts, filled)
    end subroutine batch_store_from_buffers

    !> ptcl -> class over the rows i_first..i_last of the batch's dense rows by one
    !! draw per particle from its tempered posterior over the active classes: p_c proportional to exp(-T_i (L_c -
    !! L_min)) over the valid distances, T_i the particle's temperature from the fill (temps) scaled by prob_beta; no
    !! class forced to be populated, the drawn class's entry committed (its in-plane, its refined shift when refined),
    !! frac 100. A particle without a valid distance takes the first active class; one without a temperature (none
    !! expected) its argmin. The draws are keyed on the particle; the n_eff readout takes the rows it wants in order
    subroutine assign_rows_dense( self, i_first, i_last )
        class(eul_prob_tab2D), intent(inout) :: self
        integer,               intent(in)    :: i_first, i_last
        real,    allocatable :: w(:), dists_act(:)
        real    :: temp, dist, best_loc, invalid_loc
        integer :: i_loc, iact_loc, which, iact_best, icls
        if( .not. self%l_assign_open ) THROW_HARD('no open assignment; eul_prob_tab2D%assign_rows_dense')
        invalid_loc = 0.1 * huge(invalid_loc)
        allocate(w(self%nactive), dists_act(self%nactive))
        do i_loc = i_first,i_last
            iact_best = 1
            best_loc  = invalid_loc
            do iact_loc = 1,self%nactive
                icls = self%active_cls(iact_loc)
                dists_act(iact_loc) = self%batch_rows(icls,i_loc)%dist
                if( ieee_is_finite(dists_act(iact_loc)) .and. dists_act(iact_loc) < best_loc )then
                    best_loc  = dists_act(iact_loc)
                    iact_best = iact_loc
                endif
            enddo
            if( self%neff%wanted(i_loc) ) call self%neff%offer(self%nactive, dists_act, self%temps(i_loc))
            temp = self%p_ptr%prob_beta * self%temps(i_loc)
            if( .not. (temp > 0.) )then
                self%ntemp_missing = self%ntemp_missing + 1
                temp = huge(temp)
            endif
            call draw_posterior_row(self%nactive, dists_act, temp, w, &
                &rnd_keyed(rnd_draw_key(self%rnd_scope, self%pinds(i_loc), RND_FAM_ASSIGN)), which, dist)
            if( which < 1 ) which = iact_best
            if( which == iact_best ) self%nmode = self%nmode + 1
            call self%finalize_candidate_assignment(i_loc, self%row_candidate(self%active_cls(which), i_loc), 100., 0)
            self%nassigned = self%nassigned + 1
        enddo
        deallocate(w, dists_act)
    end subroutine assign_rows_dense

    !> ptcl -> class over the rows i_first..i_last of the batch's sparse candidate store (refine=prob_neigh) by one draw per particle from its tempered posterior over its evaluated candidates (an
    !! active class, an in-plane index), the drawn candidate committed with its shift, the peak count the candidates
    !! with an in-plane index and the fraction their share of the classes; a particle without a temperature (none
    !! expected) takes its best candidate. The draws are keyed on the particle
    subroutine assign_rows_sparse( self, i_first, i_last )
        class(eul_prob_tab2D), intent(inout) :: self
        integer,               intent(in)    :: i_first, i_last
        real,    allocatable :: w(:), dists_c(:)
        integer, allocatable :: pos_c(:)
        integer :: i_loc, n_loc, which, ibest, maxcand, npeaks_loc
        real    :: temp, dist
        if( .not. self%l_assign_open ) THROW_HARD('no open assignment; eul_prob_tab2D%assign_rows_sparse')
        maxcand = 0
        do i_loc = i_first,i_last
            maxcand = max(maxcand, int(self%candidate_store%offsets(i_loc+1) - self%candidate_store%offsets(i_loc)))
        enddo
        allocate(w(max(1,maxcand)), dists_c(max(1,maxcand)), pos_c(max(1,maxcand)))
        do i_loc = i_first,i_last
            call gather_row(i_loc, dists_c, pos_c, n_loc, ibest)
            if( n_loc < 1 )&
                &THROW_HARD('Failed sparse posterior particle assignment in eul_prob_tab2D%assign_rows_sparse')
            if( self%neff%wanted(i_loc) ) call self%neff%offer(n_loc, dists_c(1:n_loc), self%temps(i_loc))
            if( ibest < 1 ) ibest = 1 ! no finite distance: any candidate
            temp = self%p_ptr%prob_beta * self%temps(i_loc)
            if( .not. (temp > 0.) )then
                self%ntemp_missing = self%ntemp_missing + 1
                temp = huge(temp)
            endif
            call draw_posterior_row(n_loc, dists_c(1:n_loc), temp, w, &
                &rnd_keyed(rnd_draw_key(self%rnd_scope, self%pinds(i_loc), RND_FAM_ASSIGN)), which, dist)
            if( which < 1 ) which = ibest ! no valid distance: the best candidate
            if( which == ibest ) self%nmode = self%nmode + 1
            npeaks_loc = count_npeaks(i_loc)
            call self%finalize_candidate_assignment(i_loc, self%candidate_store%candidates(pos_c(which)),&
                &100. * real(npeaks_loc) / real(self%nclasses), npeaks_loc)
            self%nassigned = self%nassigned + 1
        enddo
        deallocate(w, dists_c, pos_c)

    contains

        !> the particle's valid candidates into dists_row/pos_row, with the position of the best (0 without a finite one)
        subroutine gather_row( i_row, dists_row, pos_row, n_row, ibest_row )
            integer, intent(in)    :: i_row
            real,    intent(inout) :: dists_row(:)
            integer, intent(inout) :: pos_row(:)
            integer, intent(out)   :: n_row, ibest_row
            integer :: pos_loc
            real    :: best_loc
            n_row     = 0
            ibest_row = 0
            best_loc  = huge(1.0)
            do pos_loc = int(self%candidate_store%offsets(i_row)), int(self%candidate_store%offsets(i_row+1))-1
                if( .not. valid_sparse_candidate(pos_loc) ) cycle
                n_row            = n_row + 1
                dists_row(n_row) = self%candidate_store%candidates(pos_loc)%dist
                pos_row(n_row)   = pos_loc
                if( dists_row(n_row) < best_loc )then
                    best_loc  = dists_row(n_row)
                    ibest_row = n_row
                endif
            enddo
        end subroutine gather_row

        integer function count_npeaks( iptcl_loc ) result(npeaks)
            integer, intent(in) :: iptcl_loc
            integer :: pos_loc, icls_loc
            npeaks = 0
            do pos_loc = int(self%candidate_store%offsets(iptcl_loc)),&
                &int(self%candidate_store%offsets(iptcl_loc+1))-1
                icls_loc = self%candidate_store%candidates(pos_loc)%iref
                if( icls_loc < 1 .or. icls_loc > self%nclasses ) cycle
                if( self%candidate_store%candidates(pos_loc)%inpl > 0 ) npeaks = npeaks + 1
            enddo
        end function count_npeaks

        !> an active class and an in-plane index
        logical function valid_sparse_candidate( candidate_pos ) result(l_valid)
            integer, intent(in) :: candidate_pos
            integer :: icls_loc, iact_loc
            l_valid = .false.
            if( candidate_pos < 1 .or. candidate_pos > size(self%candidate_store%candidates) ) return
            icls_loc = self%candidate_store%candidates(candidate_pos)%iref
            if( icls_loc < 1 .or. icls_loc > self%nclasses ) return
            l_valid = .false.
            do iact_loc = 1,self%nactive
                if( self%active_cls(iact_loc) == icls_loc )then
                    l_valid = .true.
                    exit
                endif
            enddo
            if( .not. l_valid ) return
            if( self%candidate_store%candidates(candidate_pos)%inpl <= 0 ) l_valid = .false.
        end function valid_sparse_candidate

    end subroutine assign_rows_sparse

    ! DESTRUCTOR

    subroutine kill( self )
        class(eul_prob_tab2D), intent(inout) :: self
        integer :: ithr
        self%l_model_set        = .false.
        self%sigma2_prior_model = 0.
        call self%candidate_store%kill
        if( allocated(self%candidate_buffers) )then
            do ithr = 1,size(self%candidate_buffers)
                call self%candidate_buffers(ithr)%kill
            enddo
            deallocate(self%candidate_buffers)
        endif
        if( allocated(self%assgn_map)      ) deallocate(self%assgn_map)
        if( allocated(self%batch_rows)     ) deallocate(self%batch_rows)
        if( allocated(self%batch_refined)  ) deallocate(self%batch_refined)
        if( allocated(self%batch_nrefined) ) deallocate(self%batch_nrefined)
        if( allocated(self%active_cls)     ) deallocate(self%active_cls)
        call self%neff%kill
        self%l_assign_open = .false.
        self%nactive       = 0
        self%nassigned     = 0
        self%nmode         = 0
        self%ntemp_missing = 0
        if( allocated(self%seed_shifts)    ) deallocate(self%seed_shifts)
        if( allocated(self%seed_has_sh)    ) deallocate(self%seed_has_sh)
        if( allocated(self%temps)          ) deallocate(self%temps)
        if( allocated(self%pinds)          ) deallocate(self%pinds)
        if( allocated(self%class_exists)   ) deallocate(self%class_exists)
        self%seed_nrots = 0
        self%l_sparse_snhc = .false.
        self%fill_profile   = prob_fill_profile()
        self%b_ptr => null()
        self%p_ptr => null()
    end subroutine kill

    ! FILE IO

    subroutine init_eval2D_sparse_ws( self, nclasses, smpl_ncls, nthr )
        class(eval2D_sparse_ws), intent(inout) :: self
        integer,                 intent(in)    :: nclasses, smpl_ncls, nthr
        call self%dealloc_eval2D_sparse_ws
        allocate(self%direct_srch_order(nclasses,nthr), self%eval_cls(nclasses,nthr), source=0)
        allocate(self%best_locs(max(1,smpl_ncls),nthr), source=0)
        allocate(self%eval_dists(nclasses,nthr), source=huge(1.0))
    end subroutine init_eval2D_sparse_ws

    subroutine dealloc_eval2D_sparse_ws( self )
        class(eval2D_sparse_ws), intent(inout) :: self
        if( allocated(self%direct_srch_order) ) deallocate(self%direct_srch_order)
        if( allocated(self%eval_cls)          ) deallocate(self%eval_cls)
        if( allocated(self%best_locs)         ) deallocate(self%best_locs)
        if( allocated(self%eval_dists)        ) deallocate(self%eval_dists)
    end subroutine dealloc_eval2D_sparse_ws

end module simple_eul_prob_tab2D
