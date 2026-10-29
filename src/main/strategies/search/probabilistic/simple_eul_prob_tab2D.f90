!@descr: 2D probability table routines for multi-reference class assignment by per-particle draws from the tempered posterior
module simple_eul_prob_tab2D
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use, intrinsic :: iso_fortran_env, only: int64
use simple_pftc_srch_api
use simple_builder,            only: builder
use simple_pftc_shsrch_grad,   only: pftc_shsrch_grad, shift_ptcl2ref, shift_ref2ptcl
use simple_polarft_calc,       only: SHERRSQ
use simple_type_defs,          only: OBJFUN_EUCLID
use simple_decay_funs,         only: extremal_decay2D
use simple_eul_prob_tab_utils, only: build_pind_lookup, eulprob_dist_switch, record_seed_shift,&
    &read_seed_shift_table, write_seed_shift_table, prob_candidate,&
    &prob_candidate_buffer, prob_candidate_store, prob_fill_profile, prob_dense_entry, prob_refined_entry,&
    &posterior_temperature, draw_posterior_row, shift_prior_sigma2, prob_neff_readout, PROB_NEFF_ROWS, PROB_NEFF_NMIN
implicit none

public :: eul_prob_tab2D, prob_assign_profile2D
private
#include "simple_local_flags.inc"

real,             parameter :: NHOOD_FRAC           = 0.1
integer,          parameter :: PROB_TAB2D_IO_CHUNK  = 65536
! format words of the part files (header word 3): the candidate stream (refine=prob_snhc) and the dense rows
! (refine=prob); -6/-7 since the particle temperatures after the seed table (2026-10-23). The 3D dense part file
! is -4
integer(int64),   parameter :: PROB_TAB2D_STREAM_FORMAT = -6_int64
integer(int64),   parameter :: PROB_TAB2D_DENSE_FORMAT  = -7_int64
! status of a part file read
integer,          parameter :: PART_READ_OK = 0, PART_READ_OPEN = 1, PART_READ_IO = 2, PART_READ_FORMAT = 3
integer,          parameter :: PART_READ_NCLASSES = 4, PART_READ_COUNTS = 5, PART_READ_PTCL = 6, PART_READ_CLASS = 7
integer,          parameter :: PART_READ_REFINED = 8, PART_READ_SIZE = 9
! bytes of a dense row entry (distance, in-plane index) and of a refined entry (particle, class, shift)
integer(int64),   parameter :: DENSE_ENTRY_BYTES = 8_int64, REFINED_ENTRY_BYTES = 16_int64

!> Wall seconds of the master's part reads and assignment, for the PROB_ALIGN2D bench record
type :: prob_assign_profile2D
    real(timer_int_kind) :: read_count = 0. !< read of the parts' heads (candidates per particle, seeds)
    real(timer_int_kind) :: read_fill  = 0. !< read of the parts' candidates
    real(timer_int_kind) :: draws      = 0. !< the posterior draws, one per particle
end type prob_assign_profile2D

!> One part file as the master reads it: its head, and the master index of each of its particles
type :: tab2D_part
    integer              :: nptcls      = 0
    integer              :: nchunks     = 0
    integer              :: seed_nrots  = 0
    integer(int64)       :: nnz         = 0_int64
    integer(int64)       :: chunks_addr = 1_int64  !< the candidate chunks (stream) or the rows (dense)
    integer              :: nrefined    = 0        !< dense: refined entries in the file
    integer, allocatable :: glob_of(:)  !< master index of each particle of the part (0: not the master's)
    integer, allocatable :: counts(:)   !< stream: candidates per particle of the part
    type(prob_refined_entry), allocatable :: refined(:) !< dense: the part's refined entries, master indices
end type tab2D_part

type :: eul_prob_tab2D
    class(builder),    pointer  :: b_ptr => null()
    class(parameters), pointer  :: p_ptr => null()
    type(prob_dense_entry),   allocatable :: rows(:,:)          !< dense master table (nclasses,nptcls)
    type(prob_refined_entry), allocatable :: refined(:)         !< dense master table: the refined entries (shifts)
    integer,                  allocatable :: refined_offsets(:) !< refined entries of particle i: offsets(i):offsets(i+1)-1
    type(prob_candidate_store)        :: candidate_store !< sparse global table
    type(prob_candidate_buffer), allocatable :: candidate_buffers(:) !< worker batch storage
    type(ptcl_ref), allocatable :: assgn_map(:)      !< assignment map (nptcls)
    real,           allocatable :: seed_shifts(:,:)  !< per-particle seeded shift (2,nptcls)
    logical,        allocatable :: seed_has_sh(:)    !< per-particle seeded shift flag
    real,           allocatable :: temps(:)          !< per-particle posterior temperature pi wsqsum / (2 pftsz) (nptcls; the worker's fill writes it, the master reads it and draws at prob_beta times it)
    integer                     :: seed_nrots = 0    !< rotation grid of the table (checked when the parts are read)
    integer,        allocatable :: pinds(:)          !< particle indices for processing
    logical,        allocatable :: class_exists(:)   !< class population filter
    integer,        allocatable :: cand_counts(:)    !< candidates written per particle (part file)
    integer                     :: table_unit    = -1
    logical                     :: table_is_open = .false.
    integer(int64)              :: table_addr    = 1_int64
    integer(int64)              :: table_nnz     = 0_int64
    integer(int64)              :: table_nrefined = 0_int64 !< dense part file: refined entries written
    integer(int64)              :: rows_addr     = 1_int64  !< dense part file: first row
    integer                     :: table_nchunks = 0
    logical                     :: l_sparse_snhc = .false.
    integer                     :: nptcls            !< size of pinds array
    integer                     :: nclasses          !< number of classes
    integer                     :: nhood_sz = 1      !< classes around the best one whose shifts fill_tab refines
    type(prob_fill_profile)     :: fill_profile      !< cost record of the fills since new (bench file)
    type(prob_assign_profile2D) :: assign_profile    !< wall time of the master's reads and assignment (bench file)
    contains
    ! CONSTRUCTOR
    procedure :: new => new_global
    procedure :: new_worker
    procedure :: new_assignment
    ! MAIN PROCEDURES
    procedure :: fill_tab
    procedure :: fill_tab_range
    procedure :: get_fill_profile
    procedure :: get_assign_profile
    procedure :: ref_assign
    procedure :: begin_write
    procedure :: flush_candidate_buffers
    procedure :: write_tab
    procedure :: read_tabs_to_glob
    procedure :: write_assignment
    procedure :: read_assignment
    ! DESTRUCTOR
    procedure :: kill
    ! PRIVATE
    procedure, private :: fill_tab_prob_snhc
    procedure, private :: fill_tab_prob_snhc_range
    procedure, private :: ref_assign_sparse
    procedure, private :: record_sparse_eval
    procedure, private :: write_stream_head
    procedure, private :: write_dense_head
    procedure, private :: flush_dense_rows
    procedure, private :: dense_candidate
    procedure, private :: read_part_head
    procedure, private :: read_part_cands
    procedure, private :: new_common
    procedure, private :: initialize_storage
    procedure, private :: finalize_assignment
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

    subroutine new_global( self, params, build, pinds )
        class(eul_prob_tab2D),     intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        class(builder),    target, intent(in)    :: build
        integer,                   intent(in)    :: pinds(:)
        call self%new_common(params,build,pinds)
        allocate(self%assgn_map(self%nptcls))
        if( .not. self%l_sparse_snhc )then
            allocate(self%rows(self%nclasses,self%nptcls), self%refined(0))
            allocate(self%refined_offsets(self%nptcls+1), source=1)
        endif
        call self%initialize_storage
    end subroutine new_global

    subroutine new_worker( self, params, build, pinds )
        class(eul_prob_tab2D),     intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        class(builder),    target, intent(in)    :: build
        integer,                   intent(in)    :: pinds(:)
        call self%new_common(params,build,pinds)
        allocate(self%candidate_buffers(nthr_glob))
    end subroutine new_worker

    subroutine new_assignment( self, params, build, pinds )
        class(eul_prob_tab2D),     intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        class(builder),    target, intent(in)    :: build
        integer,                   intent(in)    :: pinds(:)
        integer :: i
        real    :: x
        call self%kill
        self%p_ptr  => params
        self%b_ptr  => build
        self%nptcls = size(pinds)
        allocate(self%pinds(self%nptcls),source=pinds)
        allocate(self%assgn_map(self%nptcls))
        do i = 1,self%nptcls
            self%assgn_map(i)%pind   = self%pinds(i)
            self%assgn_map(i)%icls   = 0
            self%assgn_map(i)%istate = 0
            self%assgn_map(i)%inpl   = 0
            self%assgn_map(i)%dist   = huge(x)
            self%assgn_map(i)%x      = 0.
            self%assgn_map(i)%y      = 0.
            self%assgn_map(i)%has_sh = .false.
            self%assgn_map(i)%frac   = 100.
            self%assgn_map(i)%npeaks = 0
        enddo
    end subroutine new_assignment

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
        self%l_sparse_snhc = trim(params%refine) == 'prob_snhc'
        ! the posterior is defined by the Euclidean objective (the particle temperature is its noise scale)
        if( params%cc_objfun /= OBJFUN_EUCLID )&
            &THROW_HARD('the 2D probability tables draw from the tempered posterior: refine=prob|prob_snhc needs objfun=euclid')
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
            if( allocated(self%rows) )then
                do icls = 1, self%nclasses
                    self%rows(icls,i)%dist = huge(x)
                    self%rows(icls,i)%inpl = 0
                end do
            endif
        end do
        !$omp end parallel do
    end subroutine initialize_storage

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
        integer :: i, icls, iptcl, ithr, irot, irot0, icls_prev, iref_n, nactive, nhood_sz_loc
        integer :: i_from, i_to
        integer :: active_cls(self%nclasses)
        integer :: loc(1)
        real    :: lims(2,2), lims_init(2,2), cxy(3), cxy_prob(3), xy0(2)
        real    :: cls_dists(self%nclasses), cls_dists_work(self%nclasses)
        real    :: inpl_dist, sigma2_prior
        logical :: class_active(self%nclasses)
        real(dp) :: ts_seed(nthr_glob), ts_score(nthr_glob), ts_refine(nthr_glob), t_region, t0
        integer(int64) :: nscored(nthr_glob), nbank(nthr_glob), nrefined(nthr_glob)
        complex(sp), allocatable :: seed_bank(:,:,:) !< per thread: the seed-shifted particle memo
        logical :: l_bank
        if( self%l_sparse_snhc )then
            call self%fill_tab_prob_snhc_range(i_first, i_last)
            return
        endif
        i_from = max(1, i_first)
        i_to   = min(self%nptcls, i_last)
        if( i_to < i_from ) return
        call seed_rnd
        self%seed_nrots = self%b_ptr%pftc%get_nrots()
        ts_seed   = 0.d0
        ts_score  = 0.d0
        ts_refine = 0.d0
        nscored   = 0_int64
        nbank     = 0_int64
        nrefined  = 0_int64
        sigma2_prior = 0.
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
            ! the shift prior's width, from the sampled particles' last shift increments (set per particle below)
            if( self%p_ptr%l_shift_prior ) sigma2_prior = shift_prior_sigma2(self%b_ptr%spproj_field, self%p_ptr%lp / self%p_ptr%smpd)
            call alloc_seed_bank2D(self, seed_bank)
            ! fill the table; dynamic: the particles' shift searches differ in cost, and each particle is processed
            ! whole by one thread with its own search object, so only its random draws change
            !$omp parallel do default(shared) private(i,iptcl,ithr,irot,irot0,cxy,icls,icls_prev,loc,&
            !$omp& cls_dists,cls_dists_work,iref_n,cxy_prob,inpl_dist,t0,l_bank,xy0)&
            !$omp firstprivate(o_prev) proc_bind(close) schedule(dynamic)
            do i = i_from, i_to
                iptcl = self%pinds(i)
                ithr  = omp_get_thread_num() + 1
                ! (1) identify shifts using the previously assigned best class
                t0 = omp_get_wtime()
                if( self%p_ptr%l_shift_prior ) call grad_shsrch_obj(ithr)%set_shift_prior(sigma2_prior,&
                    &posterior_temperature(self%b_ptr%pftc%wsqsum_ptcl_now(iptcl), self%b_ptr%pftc%get_pftsz()))
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
                        else
                            ! one particle-frame shift from here on: the estimate converted at its rotation
                            cxy(2:3) = shift_ref2ptcl(cxy(2:3), self%b_ptr%pftc%get_rot(irot))
                        endif
                    endif
                endif
                self%seed_shifts(:,i) = cxy(2:3)
                self%seed_has_sh(i)   = .true.
                self%temps(i)         = posterior_temperature(self%b_ptr%pftc%wsqsum_ptcl_now(iptcl), self%b_ptr%pftc%get_pftsz())
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
                ! (3) refine shifts for a neighborhood of classes around the best one
                t0 = omp_get_wtime()
                cls_dists_work = cls_dists
                do iref_n = 1, nhood_sz_loc
                    loc = minloc(cls_dists_work)
                    icls = loc(1)
                    if( cls_dists_work(icls) >= huge(1.0)/2.0 ) exit
                    call grad_shsrch_obj(ithr)%set_indices(icls, iptcl)
                    nrefined(ithr) = nrefined(ithr) + 1_int64
                    irot     = self%candidate_buffers(ithr)%get_inpl(i, icls)
                    xy0      = seed_start(self, cxy(2:3), irot)
                    cxy_prob = grad_shsrch_obj(ithr)%minimize(irot=irot, sh_rot=.true., xy_in=xy0)
                    if( irot > 0 )then
                        call self%record_sparse_eval(i, ithr, icls,&
                            &eulprob_dist_switch(cxy_prob(1), self%p_ptr%cc_objfun), irot,&
                            &cxy_prob(2), cxy_prob(3), .true.)
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
            &nscored, nbank, nrefined, i_to - i_from + 1, nactive, nhood_sz_loc)
        do ithr = 1, nthr_glob
            call grad_shsrch_obj(ithr)%kill
        end do
        if( allocated(seed_bank) ) deallocate(seed_bank)
        call o_prev%kill
        if( self%table_is_open ) call self%flush_candidate_buffers
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
        integer :: ithr, i, iptcl, nrefs_bound, smpl_ncls, nrots
        integer :: i_from, i_to
        real    :: lims(2,2), lims_init(2,2), neigh_frac, cxy(3), sigma2_prior
        real(dp) :: ts_seed(nthr_glob), ts_score(nthr_glob), ts_refine(nthr_glob), t_region
        integer(int64) :: nscored(nthr_glob), nbank(nthr_glob), nrefined(nthr_glob)
        complex(sp), allocatable :: seed_bank(:,:,:) !< per thread: the seed-shifted particle memo
        logical :: l_bank(nthr_glob)
        i_from = max(1, i_first)
        i_to   = min(self%nptcls, i_last)
        if( i_to < i_from ) return
        ts_seed   = 0.d0
        ts_score  = 0.d0
        ts_refine = 0.d0
        nscored   = 0_int64
        nbank     = 0_int64
        nrefined  = 0_int64
        sigma2_prior = 0.
        call seed_rnd
        nrots = self%b_ptr%pftc%get_nrots()
        self%seed_nrots = nrots
        neigh_frac  = extremal_decay2D(self%p_ptr%extr_iter, self%p_ptr%extr_lim)
        nrefs_bound = max(2, min(self%nclasses, nint(real(self%nclasses) * (1. - neigh_frac))))
        nrefs_bound = max(1, min(nrefs_bound, self%nclasses))
        smpl_ncls   = neighfrac2nsmpl(neigh_frac, self%nclasses)
        smpl_ncls   = max(1, min(smpl_ncls, self%nclasses))
        call eval_work%init_eval2D_sparse_ws(self%nclasses, smpl_ncls, nthr_glob)
        allocate(direct_rts(nthr_glob))
        do ithr = 1, nthr_glob
            direct_rts(ithr) = ran_tabu(self%nclasses)
        enddo
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
            ! the shift prior's width, from the sampled particles' last shift increments (set per particle below)
            if( self%p_ptr%l_shift_prior ) sigma2_prior = shift_prior_sigma2(self%b_ptr%spproj_field, self%p_ptr%lp / self%p_ptr%smpd)
            ! dynamic: the particles' shift searches differ in cost; each is processed whole by one thread with its
            ! own search object, so only its random draws change
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
            &nscored, nbank, nrefined, i_to - i_from + 1, 0, smpl_ncls)
        do ithr = 1, nthr_glob
            call grad_shsrch_obj(ithr)%kill
        end do
        call eval_work%dealloc_eval2D_sparse_ws
        deallocate(direct_rts, seed_bank)
        if( self%table_is_open ) call self%flush_candidate_buffers

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
            t0 = omp_get_wtime()
            call direct_rts(ithr_loc)%ne_ran_iarr(eval_work%direct_srch_order(:,ithr_loc))
            call self%b_ptr%spproj_field%get_ori(iptcl_loc, o_prev_loc)
            icls_prev_loc = nint(self%b_ptr%spproj_field%get(iptcl_loc, 'class'))
            if( icls_prev_loc >= 1 .and. icls_prev_loc <= self%nclasses )then
                call put_last(icls_prev_loc, eval_work%direct_srch_order(:,ithr_loc))
            endif
            cxy_loc = 0.
            if( self%p_ptr%l_shift_prior .and. l_with_shift ) call grad_shsrch_obj(ithr_loc)%set_shift_prior(sigma2_prior,&
                &posterior_temperature(self%b_ptr%pftc%wsqsum_ptcl_now(iptcl_loc), self%b_ptr%pftc%get_pftsz()))
            if( l_with_shift .and. icls_prev_loc >= 1 .and. icls_prev_loc <= self%nclasses )then
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
            self%temps(i_loc) = posterior_temperature(self%b_ptr%pftc%wsqsum_ptcl_now(iptcl_loc), self%b_ptr%pftc%get_pftsz())
            call memoize_seed_bank2D(self, iptcl_loc, sh_loc, seed_bank(:,:,ithr_loc), l_bank(ithr_loc),&
                &nbank(ithr_loc))
            ts_seed(ithr_loc) = ts_seed(ithr_loc) + (omp_get_wtime() - t0)
            t0 = omp_get_wtime()
            neval_bound_loc = nrefs_bound
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
            call refine_best_classes(i_loc, ithr_loc, iptcl_loc, sh_loc, neval_loc, l_with_shift)
            ts_refine(ithr_loc) = ts_refine(ithr_loc) + (omp_get_wtime() - t0)
            call o_prev_loc%kill
        end subroutine process_particle

        subroutine score_class(icls_loc, ithr_loc, iptcl_loc, sh_loc, dist_loc, irot_loc)
            integer, intent(in)  :: icls_loc, ithr_loc, iptcl_loc
            real,    intent(in)  :: sh_loc(2)
            real,    intent(out) :: dist_loc
            integer, intent(out) :: irot_loc
            call score_seeded_class(self, icls_loc, iptcl_loc, sh_loc, seed_bank(:,:,ithr_loc), l_bank(ithr_loc),&
                &dist_loc, irot_loc)
            if( irot_loc < 1 .or. irot_loc > nrots ) irot_loc = 1
        end subroutine score_class

        subroutine refine_best_classes(i_loc, ithr_loc, iptcl_loc, sh_loc, neval_loc, l_with_shift)
            integer, intent(in) :: i_loc, ithr_loc, iptcl_loc, neval_loc
            real,    intent(in) :: sh_loc(2)
            logical, intent(in) :: l_with_shift
            integer :: nrefine_loc, j_loc, eval_slot_loc, icls_loc, irot_loc
            real    :: cxy_prob_loc(3), xy0_loc(2)
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
                nrefined(ithr_loc) = nrefined(ithr_loc) + 1_int64
                irot_loc = self%candidate_buffers(ithr_loc)%get_inpl(i_loc, icls_loc)
                xy0_loc  = seed_start(self, sh_loc, irot_loc)
                cxy_prob_loc = grad_shsrch_obj(ithr_loc)%minimize(irot=irot_loc, sh_rot=.true., xy_in=xy0_loc)
                if( irot_loc > 0 )then
                    call self%record_sparse_eval(i_loc, ithr_loc, icls_loc,&
                        &eulprob_dist_switch(cxy_prob_loc(1), self%p_ptr%cc_objfun),&
                        &irot_loc, cxy_prob_loc(2), cxy_prob_loc(3), .true.)
                endif
            enddo
        end subroutine refine_best_classes

    end subroutine fill_tab_prob_snhc_range

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
    end function get_fill_profile

    !> wall time of the master's part reads and assignment since new, read before kill
    function get_assign_profile( self ) result( profile )
        class(eul_prob_tab2D), intent(in) :: self
        type(prob_assign_profile2D) :: profile
        profile = self%assign_profile
    end function get_assign_profile

    !> adds one fill to the cost record (bench file): the wall time of the parallel region, thread-seconds of the seed,
    !! the score loop and the refinement, classes scored (nrefs per particle in the plain fill, nscored in the
    !! prob_snhc fill), the most classes refined per particle and their total, shifted particle memos built and the
    !! search objects' counters
    subroutine record_fill_profile2D( self, grad_obj, wall, ts_seed, ts_score, ts_refine, nscored, nbank, nrefined, &
            &nfilled, nrefs, nrefine )
        class(eul_prob_tab2D),  intent(inout) :: self
        type(pftc_shsrch_grad), intent(inout) :: grad_obj(:)
        real(dp),               intent(in)    :: wall, ts_seed(:), ts_score(:), ts_refine(:)
        integer(int64),         intent(in)    :: nscored(:), nbank(:), nrefined(:)
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
        if( .not. self%p_ptr%l_doshift ) return
        self%fill_profile%projs_ns = max(self%fill_profile%projs_ns, nrefine)
        do ithr = 1,size(grad_obj)
            call grad_obj(ithr)%get_profile(obj_ithr, grad_ithr, alt_ithr, grid_ithr, min_ithr)
            self%fill_profile%objective_evals = self%fill_profile%objective_evals + obj_ithr
            self%fill_profile%gradient_evals  = self%fill_profile%gradient_evals  + grad_ithr
            self%fill_profile%alternations    = self%fill_profile%alternations    + alt_ithr
            self%fill_profile%grid_points     = self%fill_profile%grid_points     + grid_ithr
            self%fill_profile%minimize_calls  = self%fill_profile%minimize_calls  + min_ithr
        enddo
    end subroutine record_fill_profile2D

    subroutine record_sparse_eval( self, i, ithr, icls, dist, irot, x, y, has_sh )
        class(eul_prob_tab2D), intent(inout) :: self
        integer,               intent(in)    :: i, ithr, icls, irot
        real,                  intent(in)    :: dist, x, y
        logical,               intent(in)    :: has_sh
        type(prob_candidate) :: candidate
        if( icls < 1 .or. icls > self%nclasses )&
            &THROW_HARD('simple_eul_prob_tab2D::record_sparse_eval; class index out of range')
        candidate%iref   = icls
        candidate%dist   = dist
        candidate%inpl   = irot
        candidate%x      = x
        candidate%y      = y
        candidate%has_sh = has_sh
        call self%candidate_buffers(ithr)%append_or_replace(i, candidate)
    end subroutine record_sparse_eval


    subroutine finalize_assignment( self, iptcl_loc, icls_loc, frac, npeaks )
        class(eul_prob_tab2D), intent(inout) :: self
        integer,                intent(in)    :: iptcl_loc, icls_loc, npeaks
        real,                   intent(in)    :: frac
        call self%finalize_candidate_assignment(iptcl_loc,self%dense_candidate(icls_loc,iptcl_loc),frac,npeaks)
    end subroutine finalize_assignment

    !> the full candidate of class icls and particle i of the dense master table: its row entry, and the shift of
    !! its refined entry if it has one
    function dense_candidate( self, icls, i ) result( candidate )
        class(eul_prob_tab2D), intent(in) :: self
        integer,               intent(in) :: icls, i
        type(prob_candidate) :: candidate
        integer :: k
        candidate%iref = icls
        candidate%inpl = self%rows(icls,i)%inpl
        candidate%dist = self%rows(icls,i)%dist
        do k = self%refined_offsets(i), self%refined_offsets(i+1)-1
            if( self%refined(k)%ri /= icls ) cycle
            candidate%x      = self%refined(k)%x
            candidate%y      = self%refined(k)%y
            candidate%has_sh = .true.
            exit
        enddo
    end function dense_candidate

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
        self%assgn_map(iptcl_loc)%has_sh = candidate%has_sh
        call record_seed_shift(self%assgn_map(iptcl_loc),self%seed_shifts(:,iptcl_loc),&
            &self%seed_has_sh(iptcl_loc),self%p_ptr%l_doshift)
        self%assgn_map(iptcl_loc)%frac   = frac
        self%assgn_map(iptcl_loc)%npeaks = npeaks
    end subroutine finalize_candidate_assignment

    !> ptcl -> class assignment over the sparse candidate store (refine=prob_snhc): one draw per particle from its
    !! tempered posterior over the candidates it evaluated. The class-balanced frontier was removed 2026-10-24:
    !! if every image shows one view, every image is assigned to that class
    subroutine ref_assign_sparse( self )
        class(eul_prob_tab2D), intent(inout) :: self
        logical, allocatable :: class_active(:)
        integer :: nassigned
        integer(timer_int_kind) :: t_phase
        if( .not. allocated(self%candidate_store%candidates) )&
            &THROW_HARD('simple_eul_prob_tab2D::ref_assign_sparse; sparse candidate store is not allocated')
        allocate(class_active(self%nclasses))
        class_active = self%class_exists
        if( count(class_active) == 0 )then
            class_active = .true.
            THROW_WARN('No active classes after population filtering; falling back to all classes in eul_prob_tab2D')
        endif
        t_phase = tic()
        call assign_particles_posterior()
        self%assign_profile%draws = self%assign_profile%draws + toc(t_phase)
        deallocate(class_active)

    contains

        !> ptcl -> class by one draw per particle from its tempered posterior over its evaluated candidates:
        !! p_c proportional to exp(-T_i (L_c - L_min)) over the valid candidates (an active class, an in-plane
        !! index), T_i the particle's temperature from the fill (temps) scaled by prob_beta; no class forced to be
        !! populated, the drawn candidate committed with its shift, the peak count the candidates with an in-plane
        !! index and the fraction their share of the classes. A particle without a temperature (none expected)
        !! takes its best candidate with a warning. The draws run in particle order on the master's random stream.
        subroutine assign_particles_posterior()
            type(prob_neff_readout) :: neff
            real,    allocatable :: w(:), dists_c(:)
            integer, allocatable :: pos_c(:)
            integer :: i_loc, n_loc, which, ibest, maxcand, nmode, ntemp_missing, npeaks_loc
            real    :: temp, dist
            maxcand = 0
            do i_loc = 1,self%nptcls
                maxcand = max(maxcand, int(self%candidate_store%offsets(i_loc+1) - self%candidate_store%offsets(i_loc)))
            enddo
            allocate(w(max(1,maxcand)), dists_c(max(1,maxcand)), pos_c(max(1,maxcand)))
            ! the n_eff readout: the mean effective number of candidates over a subsample of the rows the draws are from
            call neff%new(self%nptcls, PROB_NEFF_ROWS, PROB_NEFF_NMIN)
            do i_loc = 1,self%nptcls
                if( .not. neff%wanted(i_loc) ) cycle
                call gather_row(i_loc, dists_c, pos_c, n_loc, ibest)
                if( n_loc < 1 ) cycle
                call neff%offer(n_loc, dists_c(1:n_loc), self%temps(i_loc))
            enddo
            call neff%report(self%p_ptr%prob_beta)
            call neff%kill
            nassigned     = 0
            nmode         = 0
            ntemp_missing = 0
            do i_loc = 1,self%nptcls
                call gather_row(i_loc, dists_c, pos_c, n_loc, ibest)
                if( n_loc < 1 )&
                    &THROW_HARD('Failed sparse posterior particle assignment in eul_prob_tab2D%ref_assign_sparse')
                if( ibest < 1 ) ibest = 1 ! no finite distance: any candidate
                temp = self%p_ptr%prob_beta * self%temps(i_loc)
                if( .not. (temp > 0.) )then
                    ntemp_missing = ntemp_missing + 1
                    temp = huge(temp)
                endif
                call draw_posterior_row(n_loc, dists_c(1:n_loc), temp, w, which, dist)
                if( which < 1 ) which = ibest ! no valid distance: the best candidate
                if( which == ibest ) nmode = nmode + 1
                npeaks_loc = count_npeaks(i_loc)
                call self%finalize_candidate_assignment(i_loc, self%candidate_store%candidates(pos_c(which)),&
                    &100. * real(npeaks_loc) / real(self%nclasses), npeaks_loc)
                nassigned = nassigned + 1
            enddo
            write(logfhandle,'(A,I0,A,I0,A,F6.3)') '>>> PROB_TAB2D_ASSIGN: sparse posterior draw at the particle temperature, ',&
                &nassigned, ' particles x ', self%nclasses, ' classes; fraction drawn at the row mode ',&
                &real(nmode) / real(max(1,nassigned))
            if( ntemp_missing > 0 ) write(logfhandle,'(A,I0)') &
                &'WARNING: particles without a temperature in the 2D part files, assigned to their best candidate: ', ntemp_missing
            call flush(logfhandle)
            deallocate(w, dists_c, pos_c)
        end subroutine assign_particles_posterior

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

        logical function valid_sparse_candidate( candidate_pos ) result(l_valid)
            integer, intent(in) :: candidate_pos
            integer :: icls_loc
            l_valid = .false.
            if( candidate_pos < 1 .or. candidate_pos > size(self%candidate_store%candidates) ) return
            icls_loc = self%candidate_store%candidates(candidate_pos)%iref
            if( icls_loc < 1 .or. icls_loc > self%nclasses ) return
            if( .not. class_active(icls_loc) ) return
            if( self%candidate_store%candidates(candidate_pos)%inpl <= 0 ) return
            l_valid = .true.
        end function valid_sparse_candidate

    end subroutine ref_assign_sparse

    !> ptcl -> class assignment over the dense table (refine=prob): one draw per particle from its tempered posterior
    !! over the active classes. The class-balanced frontier was removed 2026-10-24: if every image shows one
    !! view, every image is assigned to that class
    subroutine ref_assign( self )
        class(eul_prob_tab2D), intent(inout) :: self
        integer, allocatable :: active_cls(:)
        logical, allocatable :: class_active(:)
        integer :: icls, nactive
        integer(timer_int_kind) :: t_phase
        if( self%l_sparse_snhc )then
            call self%ref_assign_sparse
            return
        endif
        allocate(class_active(self%nclasses), active_cls(self%nclasses))
        class_active = self%class_exists
        nactive      = count(class_active)
        if( nactive == 0 )then
            class_active = .true.
            THROW_WARN('No active classes after population filtering; falling back to all classes in eul_prob_tab2D')
        endif
        nactive = 0
        do icls = 1, self%nclasses
            if( class_active(icls) )then
                nactive = nactive + 1
                active_cls(nactive) = icls
            endif
        end do
        t_phase = tic()
        call assign_posterior
        self%assign_profile%draws = self%assign_profile%draws + toc(t_phase)
        deallocate(active_cls, class_active)

    contains

        !> ptcl -> class by one draw per particle from its tempered posterior over the active classes:
        !! p_c proportional to exp(-T_i (L_c - L_min)) over the valid distances, T_i the particle's temperature from
        !! the fill (temps) scaled by prob_beta; no class forced to be populated, the drawn class's dense_candidate
        !! committed (its in-plane, its refined shift when refined), frac 100. A particle without a valid distance
        !! takes the first active class; one without a temperature (none expected) its argmin with a warning. The
        !! draws run in particle order on the master's random stream.
        subroutine assign_posterior
            type(prob_neff_readout) :: neff
            real    :: w(nactive), dists_act(nactive), temp, dist, best_loc, invalid_loc
            integer :: i_loc, iact_loc, which, iact_best, assigned_loc, nmode, ntemp_missing
            invalid_loc   = 0.1 * huge(invalid_loc)
            ! the n_eff readout: the mean effective number of candidates over a subsample of the rows the draws are from
            call neff%new(self%nptcls, PROB_NEFF_ROWS, PROB_NEFF_NMIN)
            do i_loc = 1,self%nptcls
                if( .not. neff%wanted(i_loc) ) cycle
                do iact_loc = 1,nactive
                    dists_act(iact_loc) = self%rows(active_cls(iact_loc),i_loc)%dist
                enddo
                call neff%offer(nactive, dists_act, self%temps(i_loc))
            enddo
            call neff%report(self%p_ptr%prob_beta)
            call neff%kill
            nmode         = 0
            ntemp_missing = 0
            do i_loc = 1,self%nptcls
                iact_best = 1
                best_loc  = invalid_loc
                do iact_loc = 1,nactive
                    dists_act(iact_loc) = self%rows(active_cls(iact_loc),i_loc)%dist
                    if( ieee_is_finite(dists_act(iact_loc)) .and. dists_act(iact_loc) < best_loc )then
                        best_loc  = dists_act(iact_loc)
                        iact_best = iact_loc
                    endif
                enddo
                temp = self%p_ptr%prob_beta * self%temps(i_loc)
                if( .not. (temp > 0.) )then
                    ntemp_missing = ntemp_missing + 1
                    temp = huge(temp)
                endif
                call draw_posterior_row(nactive, dists_act, temp, w, which, dist)
                if( which < 1 ) which = iact_best
                if( which == iact_best ) nmode = nmode + 1
                assigned_loc = active_cls(which)
                call self%finalize_assignment(i_loc, assigned_loc, 100., 0)
            enddo
            write(logfhandle,'(A,I0,A,I0,A,F6.3)') '>>> PROB_TAB2D_ASSIGN: dense posterior draw at the particle temperature, ',&
                &self%nptcls, ' particles x ', self%nclasses, ' classes; fraction drawn at the row mode ',&
                &real(nmode) / real(max(1,self%nptcls))
            if( ntemp_missing > 0 ) write(logfhandle,'(A,I0)') &
                &'WARNING: particles without a temperature in the 2D part files, assigned to their argmin: ', ntemp_missing
            call flush(logfhandle)
        end subroutine assign_posterior

    end subroutine ref_assign

    ! DESTRUCTOR

    subroutine kill( self )
        class(eul_prob_tab2D), intent(inout) :: self
        integer :: ithr
        if( self%table_is_open ) call fclose(self%table_unit)
        if( allocated(self%rows)            ) deallocate(self%rows)
        if( allocated(self%refined)         ) deallocate(self%refined)
        if( allocated(self%refined_offsets) ) deallocate(self%refined_offsets)
        call self%candidate_store%kill
        if( allocated(self%candidate_buffers) )then
            do ithr = 1,size(self%candidate_buffers)
                call self%candidate_buffers(ithr)%kill
            enddo
            deallocate(self%candidate_buffers)
        endif
        if( allocated(self%assgn_map)      ) deallocate(self%assgn_map)
        if( allocated(self%seed_shifts)    ) deallocate(self%seed_shifts)
        if( allocated(self%seed_has_sh)    ) deallocate(self%seed_has_sh)
        if( allocated(self%temps)          ) deallocate(self%temps)
        if( allocated(self%pinds)          ) deallocate(self%pinds)
        if( allocated(self%class_exists)   ) deallocate(self%class_exists)
        if( allocated(self%cand_counts) ) deallocate(self%cand_counts)
        self%seed_nrots = 0
        self%l_sparse_snhc = .false.
        self%table_unit    = -1
        self%table_is_open = .false.
        self%table_addr    = 1_int64
        self%table_nnz      = 0_int64
        self%table_nchunks  = 0
        self%table_nrefined = 0_int64
        self%rows_addr      = 1_int64
        self%fill_profile   = prob_fill_profile()
        self%assign_profile = prob_assign_profile2D()
        self%b_ptr => null()
        self%p_ptr => null()
    end subroutine kill

    ! FILE IO

    !> Opens the part file. The candidate stream (prob_snhc): header (classes, particles, format word, candidates,
    !! chunks), particle indices, candidate counts per particle (they size the master's store), seed table, particle
    !! temperatures, then the candidate chunks. The dense rows (prob): header (classes, particles, format word,
    !! refined entries), particle indices, seed table, particle temperatures, one row of (distance, in-plane index)
    !! per particle and class, then the refined entries.
    subroutine begin_write( self, binfname )
        class(eul_prob_tab2D), intent(inout) :: self
        class(string),         intent(in)    :: binfname
        integer :: io_stat
        if( self%table_is_open ) THROW_HARD('2D probability candidate stream is already open')
        if( .not. allocated(self%candidate_buffers) ) THROW_HARD('2D candidate stream requested without worker buffers')
        call fopen(self%table_unit,binfname,access='STREAM',action='WRITE',status='REPLACE',iostat=io_stat)
        call fileiochk('simple_eul_prob_tab2D; begin_write; file: '//binfname%to_char(),io_stat)
        self%table_is_open  = .true.
        self%table_nnz      = 0_int64
        self%table_nchunks  = 0
        self%table_nrefined = 0_int64
        if( self%l_sparse_snhc )then
            allocate(self%cand_counts(self%nptcls), source=0)
            call self%write_stream_head(self%table_addr)
        else
            call self%write_dense_head
            self%table_addr = self%rows_addr + int(self%nptcls,int64) * int(self%nclasses,int64) * DENSE_ENTRY_BYTES
        endif
    end subroutine begin_write

    !> header, particle indices and seed table of the dense part file; the rows follow at rows_addr
    subroutine write_dense_head( self )
        class(eul_prob_tab2D), intent(inout) :: self
        integer(int64) :: file_header(4), addr
        file_header = [int(self%nclasses,int64), int(self%nptcls,int64), PROB_TAB2D_DENSE_FORMAT, self%table_nrefined]
        write(self%table_unit,pos=1) file_header
        addr = sizeof(file_header) + 1
        write(self%table_unit,pos=addr) self%pinds
        addr = addr + sizeof(self%pinds)
        call write_seed_shift_table(self%table_unit,addr,self%seed_nrots,self%seed_shifts,self%seed_has_sh)
        write(self%table_unit,pos=addr) self%temps
        addr = addr + sizeof(self%temps)
        self%rows_addr = addr
    end subroutine write_dense_head

    !> header, particle indices, candidate counts and seed table of the part file; the chunks follow at chunks_addr
    subroutine write_stream_head( self, chunks_addr )
        class(eul_prob_tab2D), intent(inout) :: self
        integer(int64),        intent(out)   :: chunks_addr
        integer(int64) :: file_header(5)
        file_header = [int(self%nclasses,int64), int(self%nptcls,int64), PROB_TAB2D_STREAM_FORMAT, self%table_nnz,&
            &int(self%table_nchunks,int64)]
        write(self%table_unit,pos=1) file_header
        chunks_addr = sizeof(file_header) + 1
        write(self%table_unit,pos=chunks_addr) self%pinds
        chunks_addr = chunks_addr + sizeof(self%pinds)
        write(self%table_unit,pos=chunks_addr) self%cand_counts
        chunks_addr = chunks_addr + sizeof(self%cand_counts)
        call write_seed_shift_table(self%table_unit,chunks_addr,self%seed_nrots,self%seed_shifts,self%seed_has_sh)
        write(self%table_unit,pos=chunks_addr) self%temps
        chunks_addr = chunks_addr + sizeof(self%temps)
    end subroutine write_stream_head

    subroutine flush_candidate_buffers( self )
        class(eul_prob_tab2D), intent(inout) :: self
        integer :: ithr, nused, j, i_loc
        if( .not. self%table_is_open ) return
        if( .not. self%l_sparse_snhc )then
            call self%flush_dense_rows
            return
        endif
        do ithr = 1,size(self%candidate_buffers)
            nused = self%candidate_buffers(ithr)%nused
            if( nused < 1 ) cycle
            do j = 1,nused
                i_loc = self%candidate_buffers(ithr)%particle_indices(j)
                self%cand_counts(i_loc) = self%cand_counts(i_loc) + 1
            enddo
            write(self%table_unit,pos=self%table_addr) nused
            self%table_addr = self%table_addr + sizeof(nused)
            write(self%table_unit,pos=self%table_addr) self%candidate_buffers(ithr)%particle_indices(1:nused)
            self%table_addr = self%table_addr + sizeof(self%candidate_buffers(ithr)%particle_indices(1:nused))
            write(self%table_unit,pos=self%table_addr) self%candidate_buffers(ithr)%candidates(1:nused)
            self%table_addr = self%table_addr + sizeof(self%candidate_buffers(ithr)%candidates(1:nused))
            self%table_nnz     = self%table_nnz + int(nused,int64)
            self%table_nchunks = self%table_nchunks + 1
            call self%candidate_buffers(ithr)%kill
        enddo
    end subroutine flush_candidate_buffers

    !> the buffered candidates as dense rows, each at its particle's offset, and their refined entries (the ones with
    !! a shift) after the rows; a thread processes a particle whole, so its candidates are consecutive in one buffer
    subroutine flush_dense_rows( self )
        class(eul_prob_tab2D), intent(inout) :: self
        type(prob_refined_entry), allocatable :: refined(:)
        type(prob_dense_entry) :: row(self%nclasses)
        type(prob_candidate)   :: candidate
        logical :: filled(self%nclasses)
        integer :: ithr, nused, j, i_loc, icls, nrefined
        do ithr = 1,size(self%candidate_buffers)
            nused = self%candidate_buffers(ithr)%nused
            if( nused < 1 ) cycle
            allocate(refined(nused))
            j = 1
            do while( j <= nused )
                i_loc    = self%candidate_buffers(ithr)%particle_indices(j)
                filled   = .false.
                nrefined = 0
                do while( j <= nused )
                    if( self%candidate_buffers(ithr)%particle_indices(j) /= i_loc ) exit
                    candidate = self%candidate_buffers(ithr)%candidates(j)
                    icls      = candidate%iref
                    if( icls < 1 .or. icls > self%nclasses ) THROW_HARD('class out of range; eul_prob_tab2D%flush_dense_rows')
                    if( filled(icls) ) THROW_HARD('class written twice in a dense row; eul_prob_tab2D%flush_dense_rows')
                    row(icls)%dist = candidate%dist
                    row(icls)%inpl = candidate%inpl
                    filled(icls)   = .true.
                    if( candidate%has_sh )then
                        nrefined = nrefined + 1
                        refined(nrefined) = prob_refined_entry(i_loc, icls, candidate%x, candidate%y)
                    endif
                    j = j + 1
                enddo
                if( .not. all(filled) ) THROW_HARD('incomplete dense row; eul_prob_tab2D%flush_dense_rows')
                write(self%table_unit,pos=self%rows_addr+int(i_loc-1,int64)*int(self%nclasses,int64)*DENSE_ENTRY_BYTES) row
                if( nrefined > 0 )then
                    write(self%table_unit,pos=self%table_addr) refined(1:nrefined)
                    self%table_addr     = self%table_addr + int(nrefined,int64) * REFINED_ENTRY_BYTES
                    self%table_nrefined = self%table_nrefined + int(nrefined,int64)
                endif
                self%table_nnz = self%table_nnz + int(self%nclasses,int64)
            enddo
            deallocate(refined)
            call self%candidate_buffers(ithr)%kill
        enddo
    end subroutine flush_dense_rows

    !> completes the part file: the header with the counts, and the seed table
    subroutine write_tab( self, binfname )
        class(eul_prob_tab2D), intent(inout) :: self
        class(string),         intent(in)    :: binfname
        integer(int64) :: chunks_addr
        if( .not. self%table_is_open ) call self%begin_write(binfname)
        call self%flush_candidate_buffers
        ! a distributed part without sampled particles writes an empty table
        if( self%table_nnz < 1 .and. self%nptcls > 0 ) THROW_HARD('eul_prob_tab2D%write_tab; empty candidate stream')
        if( .not. self%l_sparse_snhc )then
            if( self%table_nnz /= int(self%nclasses,int64)*int(self%nptcls,int64) )&
                &THROW_HARD('eul_prob_tab2D%write_tab; dense candidate stream is incomplete')
        endif
        if( self%l_sparse_snhc )then
            call self%write_stream_head(chunks_addr)
            deallocate(self%cand_counts)
        else
            call self%write_dense_head
        endif
        call fclose(self%table_unit)
        self%table_unit    = -1
        self%table_is_open = .false.
        self%table_addr    = 1_int64
    end subroutine write_tab

    !> Reads the nparts part files into the master table, in parallel over the parts: every part's head (particle
    !! indices, candidate counts, seeds), which sizes the sparse store, then every part's candidates once. Every
    !! particle of the master must be in exactly one part, so the parts write disjoint parts of the table.
    subroutine read_tabs_to_glob( self, fbody, nparts, numlen )
        class(eul_prob_tab2D), intent(inout) :: self
        class(string),         intent(in)    :: fbody
        integer,               intent(in)    :: nparts, numlen
        type(tab2D_part), allocatable :: parts(:)
        type(string),     allocatable :: part_fnames(:)
        integer,          allocatable :: counts_glob(:), part_status(:), nparts_of(:)
        integer :: ipart, i_loc, i_glob
        integer(timer_int_kind) :: t_phase
        t_phase = tic()
        allocate(parts(nparts), part_fnames(nparts), part_status(nparts))
        allocate(nparts_of(self%nptcls), source=0)
        do ipart = 1,nparts
            part_fnames(ipart) = fbody//int2str_pad(ipart,numlen)//'.dat'
            if( part_fnames(ipart)%strlen() > XLONGSTRLEN ) THROW_HARD('2D probability table part file name too long')
        enddo
        !$omp parallel do default(shared) private(ipart) proc_bind(close) schedule(dynamic)
        do ipart = 1,nparts
            call self%read_part_head(part_fnames(ipart), parts(ipart), nparts_of, part_status(ipart))
        enddo
        !$omp end parallel do
        call raise_part_errors
        if( any(nparts_of /= 1) ) THROW_HARD('a particle of the 2D table is in no part file or in several')
        ! a part without particles never set its rotation grid
        allocate(counts_glob(self%nptcls), source=0)
        do ipart = 1,nparts
            if( parts(ipart)%nptcls > 0 )then
                if( self%seed_nrots == 0 ) self%seed_nrots = parts(ipart)%seed_nrots
                if( self%seed_nrots /= parts(ipart)%seed_nrots ) THROW_HARD('seed_nrots mismatch between 2D part files')
            endif
            if( .not. self%l_sparse_snhc ) cycle
            do i_loc = 1,parts(ipart)%nptcls
                i_glob = parts(ipart)%glob_of(i_loc)
                if( i_glob > 0 ) counts_glob(i_glob) = counts_glob(i_glob) + parts(ipart)%counts(i_loc)
            enddo
        enddo
        if( self%l_sparse_snhc )then
            if( any(counts_glob < 1) ) THROW_HARD('sparse 2D candidate stream has particles without candidates')
            call self%candidate_store%new_ragged(self%pinds,counts_glob)
            self%candidate_store%seed_shifts = self%seed_shifts
            self%candidate_store%seed_has_sh = self%seed_has_sh
        endif
        self%assign_profile%read_count = self%assign_profile%read_count + toc(t_phase)
        t_phase = tic()
        !$omp parallel do default(shared) private(ipart) proc_bind(close) schedule(dynamic)
        do ipart = 1,nparts
            call self%read_part_cands(part_fnames(ipart), parts(ipart), part_status(ipart))
        enddo
        !$omp end parallel do
        call raise_part_errors
        if( .not. self%l_sparse_snhc ) call index_refined
        self%assign_profile%read_fill = self%assign_profile%read_fill + toc(t_phase)
        deallocate(parts, part_fnames, part_status, counts_glob, nparts_of)

    contains

        !> the parts' refined entries indexed by particle, in part order; their offsets are default integers
        subroutine index_refined
            integer, allocatable :: next(:)
            integer(int64) :: nrefined_tot
            integer :: jpart, k, i
            nrefined_tot = 0_int64
            do jpart = 1,nparts
                nrefined_tot = nrefined_tot + int(size(parts(jpart)%refined),int64)
            enddo
            if( nrefined_tot >= int(huge(1),int64) ) THROW_HARD('refined entries exceed default-integer bounds; read_tabs_to_glob')
            allocate(next(self%nptcls), source=0)
            do jpart = 1,nparts
                do k = 1,size(parts(jpart)%refined)
                    i = parts(jpart)%refined(k)%i
                    next(i) = next(i) + 1
                enddo
            enddo
            self%refined_offsets(1) = 1
            do i = 1,self%nptcls
                self%refined_offsets(i+1) = self%refined_offsets(i) + next(i)
            enddo
            next = self%refined_offsets(1:self%nptcls)
            deallocate(self%refined)
            allocate(self%refined(self%refined_offsets(self%nptcls+1)-1))
            do jpart = 1,nparts
                do k = 1,size(parts(jpart)%refined)
                    i = parts(jpart)%refined(k)%i
                    self%refined(next(i)) = parts(jpart)%refined(k)
                    next(i) = next(i) + 1
                enddo
            enddo
            deallocate(next)
        end subroutine index_refined

        subroutine raise_part_errors
            character(len=:), allocatable :: err_msg
            integer :: jpart
            do jpart = 1,nparts
                if( part_status(jpart) == PART_READ_OK ) cycle
                err_msg = '2D probability table part file '//part_fnames(jpart)%to_char()//': '//&
                    &trim(part_read_error(part_status(jpart)))
                THROW_HARD(err_msg)
            enddo
        end subroutine raise_part_errors

    end subroutine read_tabs_to_glob

    !> the head of one part file: its particles' master indices and candidate counts; their seeds go to the master.
    !! nparts_of counts the parts holding each master particle. Runs inside the parallel read: no stop and no string
    !! temporaries, the caller raises the status.
    subroutine read_part_head( self, binfname, part, nparts_of, status )
        class(eul_prob_tab2D), intent(inout) :: self
        class(string),         intent(in)    :: binfname
        type(tab2D_part),      intent(inout) :: part
        integer,               intent(inout) :: nparts_of(:)
        integer,               intent(out)   :: status
        real,    allocatable :: seed_shifts_loc(:,:), temps_loc(:)
        logical, allocatable :: seed_has_sh_loc(:)
        integer, allocatable :: pinds_loc(:), pind2glob(:)
        type(prob_candidate) :: candidate_dummy
        integer(int64) :: file_header(4), nchunks, addr, file_bytes, expected_bytes, ptcl_bytes
        integer :: funit, io_stat, i_loc, i_glob, pind, max_pind, previous_parts
        logical :: l_format, flag
        call fopen(funit,binfname,access='STREAM',form='UNFORMATTED',action='READ',status='OLD',iostat=io_stat)
        if( io_stat /= 0 )then
            status = PART_READ_OPEN
            return
        endif
        status = read_head()
        call fclose(funit)

    contains

        integer function read_head() result( status_loc )
            read(funit,pos=1,iostat=io_stat) file_header
            if( io_stat /= 0 )then
                status_loc = PART_READ_IO
                return
            endif
            if( self%l_sparse_snhc )then
                l_format = file_header(3) == PROB_TAB2D_STREAM_FORMAT
            else
                l_format = file_header(3) == PROB_TAB2D_DENSE_FORMAT
            endif
            if( .not. l_format )then
                status_loc = PART_READ_FORMAT
                return
            endif
            if( file_header(1) /= int(self%nclasses,int64) )then
                status_loc = PART_READ_NCLASSES
                return
            endif
            if( file_header(2) < 0_int64 .or. file_header(2) > int(huge(i_loc),int64) .or. file_header(4) < 0_int64&
                &.or. file_header(4) > int(huge(i_loc),int64) )then
                status_loc = PART_READ_COUNTS
                return
            endif
            part%nptcls = int(file_header(2))
            if( self%l_sparse_snhc )then
                ! the stream's fifth header word
                read(funit,pos=sizeof(file_header)+1,iostat=io_stat) nchunks
                if( io_stat /= 0 )then
                    status_loc = PART_READ_IO
                    return
                endif
                if( nchunks < 0_int64 .or. nchunks > int(huge(i_loc),int64) )then
                    status_loc = PART_READ_COUNTS
                    return
                endif
                part%nnz     = file_header(4)
                part%nchunks = int(nchunks)
                if( part%nnz < 1_int64 .and. part%nptcls > 0 )then
                    status_loc = PART_READ_COUNTS
                    return
                endif
                addr = sizeof(file_header) + sizeof(nchunks) + 1
            else
                part%nrefined = int(file_header(4))
                addr = sizeof(file_header) + 1
            endif
            ! every read below lies inside the size the header implies (index, seed shift, seed flag, temperature)
            ptcl_bytes = int(storage_size(pind)/8 + 3*storage_size(1.)/8 + storage_size(flag)/8,int64)
            if( self%l_sparse_snhc )then
                expected_bytes = addr - 1 + int(part%nptcls,int64) * (ptcl_bytes + int(storage_size(pind)/8,int64)) +&
                    &int(storage_size(pind)/8,int64) * int(1 + part%nchunks,int64) +&
                    &part%nnz * int(storage_size(pind)/8 + storage_size(candidate_dummy)/8,int64)
            else
                expected_bytes = addr - 1 + int(storage_size(pind)/8,int64) +&
                    &int(part%nptcls,int64) * (ptcl_bytes + int(self%nclasses,int64) * DENSE_ENTRY_BYTES) +&
                    &int(part%nrefined,int64) * REFINED_ENTRY_BYTES
            endif
            inquire(unit=funit, size=file_bytes)
            if( file_bytes /= expected_bytes )then
                status_loc = PART_READ_SIZE
                return
            endif
            allocate(pinds_loc(part%nptcls), part%glob_of(part%nptcls),&
                &seed_shifts_loc(2,part%nptcls), seed_has_sh_loc(part%nptcls), temps_loc(part%nptcls))
            read(funit,pos=addr,iostat=io_stat) pinds_loc
            if( io_stat /= 0 )then
                status_loc = PART_READ_IO
                return
            endif
            addr = addr + sizeof(pinds_loc)
            if( self%l_sparse_snhc )then
                allocate(part%counts(part%nptcls))
                read(funit,pos=addr,iostat=io_stat) part%counts
                if( io_stat /= 0 )then
                    status_loc = PART_READ_IO
                    return
                endif
                if( any(part%counts < 0) )then
                    status_loc = PART_READ_COUNTS
                    return
                endif
                addr = addr + sizeof(part%counts)
            endif
            call read_seed_shift_table(funit, addr, part%seed_nrots, seed_shifts_loc, seed_has_sh_loc)
            if( part%nptcls > 0 )then
                read(funit,pos=addr,iostat=io_stat) temps_loc
                if( io_stat /= 0 )then
                    status_loc = PART_READ_IO
                    return
                endif
            endif
            addr = addr + sizeof(temps_loc)
            part%chunks_addr = addr
            call build_pind_lookup(self%pinds, pinds_loc, pind2glob, max_pind)
            part%glob_of = 0
            do i_loc = 1,part%nptcls
                pind = pinds_loc(i_loc)
                if( pind < 1 .or. pind > max_pind ) cycle
                i_glob = pind2glob(pind)
                if( i_glob < 1 ) cycle
                !$omp atomic capture
                previous_parts      = nparts_of(i_glob)
                nparts_of(i_glob)   = nparts_of(i_glob) + 1
                !$omp end atomic
                if( previous_parts > 0 )then
                    status_loc = PART_READ_PTCL
                    return
                endif
                part%glob_of(i_loc)        = i_glob
                self%seed_shifts(:,i_glob) = seed_shifts_loc(:,i_loc)
                self%seed_has_sh(i_glob)   = seed_has_sh_loc(i_loc)
                self%temps(i_glob)         = temps_loc(i_loc)
            enddo
            status_loc = PART_READ_OK
        end function read_head

    end subroutine read_part_head

    !> the candidates of one part file, once, to the master's particles: the dense rows and the part's refined
    !! entries, or the sparse store in stream order. No stop and no string temporaries: the caller raises the status.
    subroutine read_part_cands( self, binfname, part, status )
        class(eul_prob_tab2D), intent(inout) :: self
        class(string),         intent(in)    :: binfname
        type(tab2D_part),      intent(inout) :: part
        integer,               intent(out)   :: status
        type(prob_candidate), allocatable :: candidates_loc(:)
        integer,              allocatable :: particle_indices(:), nread_loc(:)
        integer(int64) :: addr, nnz_read, indices_addr, candidates_addr, pos
        integer :: funit, io_stat
        call fopen(funit,binfname,access='STREAM',form='UNFORMATTED',action='READ',status='OLD',iostat=io_stat)
        if( io_stat /= 0 )then
            status = PART_READ_OPEN
            return
        endif
        if( self%l_sparse_snhc )then
            status = read_cands()
        else
            status = read_rows()
        endif
        call fclose(funit)

    contains

        integer function read_rows() result( status_loc )
            type(prob_dense_entry),   allocatable :: rows_loc(:,:)
            type(prob_refined_entry), allocatable :: refined_loc(:)
            integer :: i_loc, i_glob, k, n
            allocate(rows_loc(self%nclasses,part%nptcls), refined_loc(part%nrefined))
            if( part%nptcls > 0 )then
                read(funit,pos=part%chunks_addr,iostat=io_stat) rows_loc
                if( io_stat /= 0 )then
                    status_loc = PART_READ_IO
                    return
                endif
            endif
            do i_loc = 1,part%nptcls
                i_glob = part%glob_of(i_loc)
                if( i_glob > 0 ) self%rows(:,i_glob) = rows_loc(:,i_loc)
            enddo
            if( part%nrefined > 0 )then
                read(funit,pos=part%chunks_addr+int(part%nptcls,int64)*int(self%nclasses,int64)*DENSE_ENTRY_BYTES,&
                    &iostat=io_stat) refined_loc
                if( io_stat /= 0 )then
                    status_loc = PART_READ_IO
                    return
                endif
            endif
            ! the refined entries, by the master's particle index
            n = 0
            do k = 1,part%nrefined
                if( refined_loc(k)%i < 1 .or. refined_loc(k)%i > part%nptcls .or. refined_loc(k)%ri < 1 .or.&
                    &refined_loc(k)%ri > self%nclasses )then
                    status_loc = PART_READ_REFINED
                    return
                endif
                i_glob = part%glob_of(refined_loc(k)%i)
                if( i_glob < 1 ) cycle
                n = n + 1
                refined_loc(n)   = refined_loc(k)
                refined_loc(n)%i = i_glob
            enddo
            part%refined = refined_loc(1:n)
            status_loc = PART_READ_OK
        end function read_rows

        integer function read_cands() result( status_loc )
            integer :: ichunk, chunk_n, first, last, nread, j, i_loc, i_glob, icls
            allocate(particle_indices(PROB_TAB2D_IO_CHUNK), candidates_loc(PROB_TAB2D_IO_CHUNK))
            allocate(nread_loc(part%nptcls), source=0)
            addr     = part%chunks_addr
            nnz_read = 0_int64
            do ichunk = 1,part%nchunks
                read(funit,pos=addr,iostat=io_stat) chunk_n
                if( io_stat /= 0 )then
                    status_loc = PART_READ_IO
                    return
                endif
                if( chunk_n < 1 .or. nnz_read + int(chunk_n,int64) > part%nnz )then
                    status_loc = PART_READ_COUNTS
                    return
                endif
                indices_addr    = addr + int(sizeof(chunk_n),int64)
                candidates_addr = indices_addr + int(chunk_n,int64)*int(sizeof(chunk_n),int64)
                do first = 1,chunk_n,PROB_TAB2D_IO_CHUNK
                    last  = min(chunk_n,first+PROB_TAB2D_IO_CHUNK-1)
                    nread = last-first+1
                    read(funit,pos=indices_addr+int(first-1,int64)*int(sizeof(chunk_n),int64),iostat=io_stat)&
                        &particle_indices(1:nread)
                    if( io_stat /= 0 )then
                        status_loc = PART_READ_IO
                        return
                    endif
                    read(funit,pos=candidates_addr+int(first-1,int64)*int(sizeof(candidates_loc(1)),int64),&
                        &iostat=io_stat) candidates_loc(1:nread)
                    if( io_stat /= 0 )then
                        status_loc = PART_READ_IO
                        return
                    endif
                    do j = 1,nread
                        i_loc = particle_indices(j)
                        if( i_loc < 1 .or. i_loc > part%nptcls )then
                            status_loc = PART_READ_PTCL
                            return
                        endif
                        icls = candidates_loc(j)%iref
                        if( icls < 1 .or. icls > self%nclasses )then
                            status_loc = PART_READ_CLASS
                            return
                        endif
                        nread_loc(i_loc) = nread_loc(i_loc) + 1
                        if( nread_loc(i_loc) > part%counts(i_loc) )then
                            status_loc = PART_READ_COUNTS
                            return
                        endif
                        i_glob = part%glob_of(i_loc)
                        if( i_glob < 1 ) cycle
                        pos = self%candidate_store%offsets(i_glob) + int(nread_loc(i_loc) - 1,int64)
                        self%candidate_store%candidates(pos) = candidates_loc(j)
                    enddo
                enddo
                addr     = candidates_addr + int(chunk_n,int64)*int(sizeof(candidates_loc(1)),int64)
                nnz_read = nnz_read + int(chunk_n,int64)
            enddo
            if( nnz_read /= part%nnz .or. any(nread_loc /= part%counts) )then
                status_loc = PART_READ_COUNTS
                return
            endif
            status_loc = PART_READ_OK
        end function read_cands

    end subroutine read_part_cands

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

    function part_read_error( status ) result( reason )
        integer, intent(in) :: status
        character(len=64)   :: reason
        select case(status)
            case(PART_READ_OPEN);     reason = 'cannot be opened'
            case(PART_READ_IO);       reason = 'read error'
            case(PART_READ_FORMAT);   reason = 'not a 2D part file of this format'
            case(PART_READ_NCLASSES); reason = 'class count differs from the master'
            case(PART_READ_COUNTS);   reason = 'invalid or inconsistent candidate counts'
            case(PART_READ_PTCL);     reason = 'invalid particle index, or a particle in several parts'
            case(PART_READ_CLASS);    reason = 'invalid class index'
            case(PART_READ_REFINED);  reason = 'invalid refined entry'
            case(PART_READ_SIZE);     reason = 'file size differs from the size its header implies'
            case DEFAULT;             reason = 'unknown error'
        end select
    end function part_read_error

end module simple_eul_prob_tab2D
