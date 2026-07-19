!@descr: the core probability table routines used for probabilistic 3D search
module simple_eul_prob_tab
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use, intrinsic :: iso_fortran_env,  only: int64
use simple_pftc_srch_api
use simple_builder,          only: builder
use simple_eul_prob_tab_utils, only: build_pind_lookup, calc_athres, calc_num2sample,&
    &eulprob_dist_switch, materialize_seed_shift, read_seed_shift_table, sample_likelihood_dist,&
    &write_seed_shift_table, prob_candidate, prob_candidate_buffer
use simple_pftc_shsrch_grad, only: pftc_shsrch_grad
use simple_type_defs,        only: OBJFUN_EUCLID
use simple_prob_prior3D,     only: prior3d_writer, PRIOR3D_NREMAP
use simple_segmentation,     only: detect_peak_thres_fdr
implicit none

public :: eul_prob_tab
private
#include "simple_local_flags.inc"

integer, parameter :: PROB_TAB_IO_CHUNK = 1024
real,    parameter :: PRIOR3D_FDR_Q      = 0.25 !< BH-FDR q: expected false-discovery fraction among retained directions

type :: eul_prob_tab
    class(builder),    pointer  :: b_ptr => null()
    class(parameters), pointer  :: p_ptr => null()
    type(prob_candidate), allocatable :: loc_tab(:,:)   !< evaluated 3D candidates (active refs,nptcls)
    type(prob_candidate), allocatable :: state_tab(:,:) !< state-only candidates (active states,nptcls)
    type(ptcl_ref), allocatable :: assgn_map(:)    !< assignment map                  (nptcls)
    real,           allocatable :: seed_shifts(:,:) !< per-particle seeded shift       (2,nptcls)
    logical,        allocatable :: seed_has_sh(:)   !< per-particle seeded shift flag  (nptcls)
    integer                     :: seed_nrots = 0   !< rotation grid used for deferred seeded shifts
    integer,        allocatable :: pinds(:)        !< particle indices for processing
    integer,        allocatable :: ssinds(:)       !< non-empty state indices
    integer,        allocatable :: state_to_active_rank(:) !< requested state -> compact active-state rank
    logical,        allocatable :: state_exists(:)
    type(prob_candidate_buffer), allocatable :: candidate_buffers(:)
    integer                     :: table_unit    = -1
    logical                     :: table_is_open = .false.
    integer(int64)              :: table_addr    = 1_int64
    integer(int64)              :: table_nnz     = 0_int64
    integer                     :: table_nchunks = 0
    integer                     :: nptcls          !< size of pinds array
    integer                     :: nstates         !< states number
    integer                     :: nrefs           !< reference number
    contains
    ! CONSTRUCTOR
    procedure :: new => new_global
    procedure :: new_state
    procedure :: new_worker
    procedure :: new_compact_global
    procedure :: new_assignment
    procedure, private :: new_common
    procedure, private :: initialize_storage
    ! PARTITION-WISE PROCEDURES (used only by partition-wise eul_prob_tab objects)
    procedure :: fill_tab
    procedure :: fill_tab_range
    procedure :: fill_tab_state_only
    procedure :: fill_tab_state_only_range
    procedure :: write_tab
    procedure :: begin_write
    procedure :: flush_candidate_buffers
    procedure :: write_state_tab
    procedure :: read_assignment
    ! GLOBAL PROCEDURES (used only by the global eul_prob_tab object)
    procedure :: read_state_tab
    procedure :: read_tab_to_glob
    procedure :: ref_assign
    procedure :: write_assignment
    procedure :: write_prior3D
    procedure :: state_assign
    ! REFERENCE INDEX MAPPING
    procedure :: ref_state
    procedure :: ref_proj
    procedure :: ref_full
    procedure :: full_to_compact_ref
    procedure :: assign_candidate
    ! DESTRUCTOR
    procedure :: kill
end type eul_prob_tab

contains

    ! CONSTRUCTORS

    subroutine new_global( self, params, build, pinds )
        class(eul_prob_tab),       intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        class(builder),    target, intent(in)    :: build
        integer,                   intent(in)    :: pinds(:)
        call self%new_common(params,build,pinds)
        allocate(self%assgn_map(self%nptcls),self%loc_tab(self%nrefs,self%nptcls))
        call self%initialize_storage
    end subroutine new_global

    subroutine new_state( self, params, build, pinds )
        class(eul_prob_tab),       intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        class(builder),    target, intent(in)    :: build
        integer,                   intent(in)    :: pinds(:)
        call self%new_common(params,build,pinds)
        allocate(self%assgn_map(self%nptcls),self%state_tab(self%nstates,self%nptcls))
        call self%initialize_storage
    end subroutine new_state

    subroutine new_worker( self, params, build, pinds )
        class(eul_prob_tab),       intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        class(builder),    target, intent(in)    :: build
        integer,                   intent(in)    :: pinds(:)
        call self%new_common(params,build,pinds)
        allocate(self%candidate_buffers(nthr_glob))
    end subroutine new_worker

    subroutine new_compact_global( self, params, build, pinds )
        class(eul_prob_tab),       intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        class(builder),    target, intent(in)    :: build
        integer,                   intent(in)    :: pinds(:)
        call self%new_common(params,build,pinds)
        allocate(self%assgn_map(self%nptcls))
        call self%initialize_storage
    end subroutine new_compact_global

    subroutine new_common( self, params, build, pinds )
        class(eul_prob_tab),       intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        class(builder),    target, intent(in)    :: build
        integer,                   intent(in)    :: pinds(:)
        integer, parameter :: MIN_POP = 5   ! ignoring cavgs with less than 5 particles
        integer :: istate, si
        call self%kill
        self%p_ptr => params
        self%b_ptr  => build
        self%nptcls       = size(pinds)
        self%state_exists = self%b_ptr%spproj_field%states_exist(self%p_ptr%nstates, thres=MIN_POP)
        self%nstates      = count(self%state_exists .eqv. .true.)
        self%nrefs = self%p_ptr%nspace * self%nstates
        allocate(self%ssinds(self%nstates), self%state_to_active_rank(self%p_ptr%nstates), source=0)
        si = 0
        do istate = 1, self%p_ptr%nstates
            if( .not. self%state_exists(istate) )cycle
            si                                = si + 1
            self%ssinds(si)                   = istate
            self%state_to_active_rank(istate) = si
        enddo
        allocate(self%pinds(self%nptcls), source=pinds)
        allocate(self%seed_shifts(2,self%nptcls), source=0.)
        allocate(self%seed_has_sh(self%nptcls), source=.false.)
    end subroutine new_common

    subroutine initialize_storage( self )
        class(eul_prob_tab), intent(inout) :: self
        integer :: i, iptcl, si, ri
        real    :: x
        !$omp parallel do default(shared) private(i,iptcl,si,ri) proc_bind(close) schedule(static)
        do i = 1,self%nptcls
            iptcl = self%pinds(i)
            if( allocated(self%assgn_map) )then
                self%assgn_map(i)%pind   = iptcl
                self%assgn_map(i)%istate = 0
                self%assgn_map(i)%iproj  = 0
                self%assgn_map(i)%inpl   = 0
                self%assgn_map(i)%dist   = huge(x)
                self%assgn_map(i)%x      = 0.
                self%assgn_map(i)%y      = 0.
                self%assgn_map(i)%has_sh = .false.
                self%assgn_map(i)%npeaks = 0
            endif
            if( allocated(self%state_tab) )then
                do si = 1,self%nstates
                    self%state_tab(si,i)%iref   = 0
                    self%state_tab(si,i)%inpl   = 0
                    self%state_tab(si,i)%dist   = huge(x)
                    self%state_tab(si,i)%x      = 0.
                    self%state_tab(si,i)%y      = 0.
                    self%state_tab(si,i)%has_sh = .false.
                enddo
            endif
            if( allocated(self%loc_tab) )then
                do ri = 1, self%nrefs
                    self%loc_tab(ri,i)%iref   = self%ref_full(ri)
                    self%loc_tab(ri,i)%inpl   = 0
                    self%loc_tab(ri,i)%dist   = huge(x)
                    self%loc_tab(ri,i)%x      = 0.
                    self%loc_tab(ri,i)%y      = 0.
                    self%loc_tab(ri,i)%has_sh = .false.
                end do
            endif
        end do
        !$omp end parallel do
    end subroutine initialize_storage

    subroutine new_assignment( self, params, build, pinds )
        class(eul_prob_tab),       intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        class(builder),    target, intent(in)    :: build
        integer,                   intent(in)    :: pinds(:)
        integer :: i, iptcl
        real    :: x
        call self%kill
        self%p_ptr  => params
        self%b_ptr  => build
        self%nptcls = size(pinds)
        allocate(self%pinds(self%nptcls), source=pinds)
        allocate(self%assgn_map(self%nptcls))
        do i = 1,self%nptcls
            iptcl = self%pinds(i)
            self%assgn_map(i)%pind   = iptcl
            self%assgn_map(i)%istate = 0
            self%assgn_map(i)%iproj  = 0
            self%assgn_map(i)%inpl   = 0
            self%assgn_map(i)%dist   = huge(x)
            self%assgn_map(i)%x      = 0.
            self%assgn_map(i)%y      = 0.
            self%assgn_map(i)%has_sh = .false.
            self%assgn_map(i)%npeaks = 0
        enddo
    end subroutine new_assignment

    ! partition-wise table filling, used only in shared-memory commander 'exec_prob_tab'
    subroutine fill_tab( self )
        class(eul_prob_tab), intent(inout) :: self
        call self%fill_tab_range(1, self%nptcls)
    end subroutine fill_tab

    subroutine fill_tab_range( self, i_first, i_last )
        class(eul_prob_tab), intent(inout) :: self
        integer,             intent(in)    :: i_first, i_last
        integer, allocatable   :: locn(:,:)
        type(pftc_shsrch_grad) :: grad_shsrch_obj(nthr_glob) !< origin shift search object, L-BFGS with gradient
        integer :: i, si, iptcl, n, projs_ns, ithr, inds_sorted(self%b_ptr%pftc%get_nrots(),nthr_glob), istate
        integer :: inpl_refs(self%nrefs,nthr_glob)
        real    :: lims(2,2), lims_init(2,2), inpl_athres(self%p_ptr%nstates)
        real    :: dists_inpl_sorted(self%b_ptr%pftc%get_nrots(),nthr_glob), dists_refs(self%nrefs,nthr_glob)
        logical :: l_sh_first
        if( i_first < 1 .or. i_last > self%nptcls .or. i_last < i_first )then
            THROW_HARD('invalid particle range in eul_prob_tab%fill_tab_range')
        endif
        self%seed_nrots = self%b_ptr%pftc%get_nrots()
        l_sh_first      = self%p_ptr%l_doshift .and. self%p_ptr%nstates <= 1
        call seed_rnd
        projs_ns = 0
        do si = 1, self%nstates
            istate = self%ssinds(si)
            call calc_num2sample(self%b_ptr%spproj_field, self%p_ptr%nspace, 'dist', n, self%p_ptr%prob_athres, state=istate)
            projs_ns            = max(projs_ns, n)
            inpl_athres(istate) = calc_athres(self%b_ptr%spproj_field, 'dist_inpl', self%p_ptr%prob_athres, state=istate)
        enddo
        if( allocated(locn) ) deallocate(locn)
        allocate(locn(projs_ns,nthr_glob), source=0)
        if( self%p_ptr%l_doshift )then
            ! make shift search objects
            lims(:,1)      = -self%p_ptr%trs
            lims(:,2)      =  self%p_ptr%trs
            lims_init(:,1) = -SHC_INPL_TRSHWDTH
            lims_init(:,2) =  SHC_INPL_TRSHWDTH
            do ithr = 1,nthr_glob
                call grad_shsrch_obj(ithr)%new(self%b_ptr, lims, lims_init=lims_init, shbarrier=self%p_ptr%shbarrier,&
                    &maxits=self%p_ptr%maxits_sh, opt_angle=.true., coarse_init=.true.)
            end do
            ! fill the table
            !$omp parallel do default(shared) private(i,iptcl,ithr) proc_bind(close) schedule(static)
            do i = i_first, i_last
                iptcl = self%pinds(i)
                ithr  = omp_get_thread_num() + 1
                call process_particle_with_shift(i, iptcl, ithr)
            enddo
            !$omp end parallel do
        else
            ! fill the table
            !$omp parallel do default(shared) private(i,iptcl,ithr) proc_bind(close) schedule(static)
            do i = i_first, i_last
                iptcl = self%pinds(i)
                ithr  = omp_get_thread_num() + 1
                call process_particle_without_shift(i, iptcl, ithr)
            enddo
            !$omp end parallel do
        endif
        do ithr = 1,nthr_glob
            call grad_shsrch_obj(ithr)%kill
        end do
        if( self%table_is_open ) call self%flush_candidate_buffers

    contains

        subroutine process_particle_with_shift( i_loc, iptcl_loc, ithr_loc )
            integer, intent(in) :: i_loc, iptcl_loc, ithr_loc
            type(ori) :: o_prev_loc
            integer   :: prev_state, prev_proj, ri_loc, irot_loc
            real      :: shift_seed(3), dist_loc
            call get_particle_context(iptcl_loc, o_prev_loc, prev_state, prev_proj)
            call estimate_shift_seed(iptcl_loc, ithr_loc, prev_state, prev_proj, o_prev_loc, shift_seed)
            self%seed_shifts(:,i_loc) = shift_seed(2:3)
            self%seed_has_sh(i_loc)   = l_sh_first
            do ri_loc = 1,self%nrefs
                call score_ref(ri_loc, iptcl_loc, shift_seed(2:3), ithr_loc, dist_loc, irot_loc)
                call record_ref_eval(i_loc, ithr_loc, ri_loc, dist_loc, irot_loc, 0., 0., .false.)
                dists_refs(ri_loc,ithr_loc) = dist_loc
                inpl_refs(ri_loc,ithr_loc)  = irot_loc
            enddo
            call refine_best_refs(i_loc, iptcl_loc, ithr_loc, shift_seed)
            call o_prev_loc%kill
        end subroutine process_particle_with_shift

        subroutine process_particle_without_shift( i_loc, iptcl_loc, ithr_loc )
            integer, intent(in) :: i_loc, iptcl_loc, ithr_loc
            integer :: ri_loc, irot_loc
            real    :: dist_loc
            do ri_loc = 1,self%nrefs
                call score_ref(ri_loc, iptcl_loc, [0.,0.], ithr_loc, dist_loc, irot_loc)
                call record_ref_eval(i_loc, ithr_loc, ri_loc, dist_loc, irot_loc, 0., 0., .false.)
            enddo
        end subroutine process_particle_without_shift

        subroutine get_particle_context( iptcl_loc, o_prev_loc, prev_state, prev_proj )
            integer,   intent(in)    :: iptcl_loc
            type(ori), intent(inout) :: o_prev_loc
            integer,   intent(out)   :: prev_state, prev_proj
            call self%b_ptr%spproj_field%get_ori(iptcl_loc, o_prev_loc)
            prev_state = o_prev_loc%get_state()
            prev_proj  = self%b_ptr%eulspace%find_closest_proj(o_prev_loc)
        end subroutine get_particle_context

        subroutine estimate_shift_seed( iptcl_loc, ithr_loc, prev_state, prev_proj, o_prev_loc, shift_seed )
            integer,   intent(in)    :: iptcl_loc, ithr_loc, prev_state, prev_proj
            type(ori), intent(inout) :: o_prev_loc
            real,      intent(out)   :: shift_seed(3)
            integer :: irot, iref
            shift_seed = 0.
            if( .not. l_sh_first ) return
            if( prev_state < 1 .or. prev_state > self%p_ptr%nstates ) return
            if( prev_proj  < 1 .or. prev_proj  > self%p_ptr%nspace  ) return
            if( .not. self%state_exists(prev_state) ) return
            iref = (prev_state-1)*self%p_ptr%nspace + prev_proj
            irot = self%b_ptr%pftc%get_roind(360.-o_prev_loc%e3get())
            call grad_shsrch_obj(ithr_loc)%set_indices(iref, iptcl_loc)
            shift_seed = grad_shsrch_obj(ithr_loc)%minimize(irot=irot, sh_rot=.false.)
            if( irot == 0 ) shift_seed(2:3) = 0.
        end subroutine estimate_shift_seed

        subroutine score_ref( ri_loc, iptcl_loc, shift_xy, ithr_loc, dist_loc, irot_loc )
            integer, intent(in)  :: ri_loc, iptcl_loc, ithr_loc
            real,    intent(in)  :: shift_xy(2)
            real,    intent(out) :: dist_loc
            integer, intent(out) :: irot_loc
            integer :: istate_loc, iproj_loc, full_ref
            real    :: corr_loc
            istate_loc = self%ref_state(ri_loc)
            iproj_loc  = self%ref_proj(ri_loc)
            full_ref   = self%ref_full(ri_loc)
            call self%b_ptr%pftc%gen_likelihood_val(full_ref, iptcl_loc, shift_xy,&
                &inpl_likelihood_nsample(istate_loc), dist_loc, corr_loc, irot_loc,&
                &dists_inpl_sorted(:,ithr_loc), inds_sorted(:,ithr_loc))
        end subroutine score_ref

        integer function inpl_likelihood_nsample( istate_loc ) result(nsample)
            integer, intent(in) :: istate_loc
            real :: athres_ub
            athres_ub = min(self%p_ptr%prob_athres, inpl_athres(istate_loc))
            nsample = min(self%b_ptr%pftc%get_nrots(), max(1, int(athres_ub * real(self%b_ptr%pftc%get_nrots()) / 180.)))
        end function inpl_likelihood_nsample

        subroutine refine_best_refs( i_loc, iptcl_loc, ithr_loc, shift_seed )
            integer, intent(in) :: i_loc, iptcl_loc, ithr_loc
            real,    intent(in) :: shift_seed(3)
            integer :: j_loc, ri_loc, istate_loc, iproj_loc, irot_loc
            real    :: refined_shift(3)
            locn(:,ithr_loc) = minnloc(dists_refs(:,ithr_loc), projs_ns)
            do j_loc = 1,projs_ns
                ri_loc     = locn(j_loc,ithr_loc)
                istate_loc = self%ref_state(ri_loc)
                iproj_loc  = self%ref_proj(ri_loc)
                call grad_shsrch_obj(ithr_loc)%set_indices(self%ref_full(ri_loc), iptcl_loc)
                irot_loc = inpl_refs(ri_loc,ithr_loc)
                if( l_sh_first )then
                    refined_shift = grad_shsrch_obj(ithr_loc)%minimize(irot=irot_loc, sh_rot=.true., xy_in=shift_seed(2:3))
                else
                    refined_shift = grad_shsrch_obj(ithr_loc)%minimize(irot=irot_loc, sh_rot=.true.)
                endif
                if( irot_loc > 0 )then
                    call replace_ref_eval(i_loc,ithr_loc,ri_loc,&
                        &eulprob_dist_switch(refined_shift(1),self%p_ptr%cc_objfun),irot_loc,&
                        &refined_shift(2),refined_shift(3),.true.)
                endif
            enddo
        end subroutine refine_best_refs

        subroutine record_ref_eval( i_loc, ithr_loc, ri_loc, dist_loc, irot_loc, x_loc, y_loc, has_sh_loc )
            integer, intent(in) :: i_loc, ithr_loc, ri_loc, irot_loc
            real,    intent(in) :: dist_loc, x_loc, y_loc
            logical, intent(in) :: has_sh_loc
            type(prob_candidate) :: candidate
            candidate = make_ref_candidate(ri_loc,dist_loc,irot_loc,x_loc,y_loc,has_sh_loc)
            call self%candidate_buffers(ithr_loc)%append(i_loc,candidate)
        end subroutine record_ref_eval

        subroutine replace_ref_eval( i_loc, ithr_loc, ri_loc, dist_loc, irot_loc, x_loc, y_loc, has_sh_loc )
            integer, intent(in) :: i_loc, ithr_loc, ri_loc, irot_loc
            real,    intent(in) :: dist_loc, x_loc, y_loc
            logical, intent(in) :: has_sh_loc
            type(prob_candidate) :: candidate
            candidate = make_ref_candidate(ri_loc,dist_loc,irot_loc,x_loc,y_loc,has_sh_loc)
            call self%candidate_buffers(ithr_loc)%append_or_replace(i_loc,candidate)
        end subroutine replace_ref_eval

        function make_ref_candidate( ri_loc, dist_loc, irot_loc, x_loc, y_loc, has_sh_loc ) result(candidate)
            integer, intent(in) :: ri_loc, irot_loc
            real,    intent(in) :: dist_loc, x_loc, y_loc
            logical, intent(in) :: has_sh_loc
            type(prob_candidate) :: candidate
            candidate%iref   = self%ref_full(ri_loc)
            candidate%dist   = dist_loc
            candidate%inpl   = irot_loc
            candidate%x      = x_loc
            candidate%y      = y_loc
            candidate%has_sh = has_sh_loc
        end function make_ref_candidate

    end subroutine fill_tab_range

    subroutine fill_tab_state_only( self )
        class(eul_prob_tab), intent(inout) :: self
        call self%fill_tab_state_only_range(1, self%nptcls)
    end subroutine fill_tab_state_only

    subroutine fill_tab_state_only_range( self, i_first, i_last )
        class(eul_prob_tab), intent(inout) :: self
        integer,             intent(in)    :: i_first, i_last
        type(pftc_shsrch_grad) :: grad_shsrch_obj(nthr_glob)  !< origin shift search object, L-BFGS with gradient
        type(ori)               :: o_prev
        type(prob_candidate)    :: candidate
        integer :: i, iproj, iptcl, ithr, irot, istate, iref, is
        real    :: lims(2,2), lims_init(2,2), cxy(3)
        if( i_first < 1 .or. i_last > self%nptcls .or. i_last < i_first )then
            THROW_HARD('invalid particle range in eul_prob_tab%fill_tab_state_only_range')
        endif
        call seed_rnd
        if( self%p_ptr%l_doshift )then
            ! make shift search objects
            lims(:,1)      = -self%p_ptr%trs
            lims(:,2)      =  self%p_ptr%trs
            lims_init(:,1) = -SHC_INPL_TRSHWDTH
            lims_init(:,2) =  SHC_INPL_TRSHWDTH
            do ithr = 1,nthr_glob
                call grad_shsrch_obj(ithr)%new(self%b_ptr, lims, lims_init=lims_init, shbarrier=self%p_ptr%shbarrier,&
                    &maxits=self%p_ptr%maxits_sh, opt_angle=.true.)
            end do
            ! fill the table
            !$omp parallel do default(shared) private(i,iptcl,ithr,o_prev,iproj,is,istate,irot,iref,cxy,candidate)&
            !$omp proc_bind(close) schedule(static)
            do i = i_first, i_last
                iptcl = self%pinds(i)
                ithr  = omp_get_thread_num() + 1
                ! identify shifts using the previously assigned best reference
                call self%b_ptr%spproj_field%get_ori(iptcl, o_prev)     ! previous ori
                irot   = self%b_ptr%pftc%get_roind(360.-o_prev%e3get()) ! in-plane angle index
                iproj = self%b_ptr%eulspace%find_closest_proj(o_prev)   ! previous projection direction
                do is = 1, self%nstates
                    istate = self%ssinds(is)
                    iref   = (istate-1)*self%p_ptr%nspace + iproj
                    ! BFGS over shifts
                    call grad_shsrch_obj(ithr)%set_indices(iref, iptcl)
                    cxy = grad_shsrch_obj(ithr)%minimize(irot=irot, sh_rot=.true.)
                    if( irot == 0 )then
                        irot     = self%b_ptr%pftc%get_roind(360.-o_prev%e3get())
                        cxy(1)   = real(self%b_ptr%pftc%gen_corr_for_rot_8(iref, iptcl, irot))
                        cxy(2:3) = 0.
                    endif
                    candidate%dist   = eulprob_dist_switch(cxy(1), self%p_ptr%cc_objfun)
                    candidate%iref   = iref
                    candidate%inpl   = irot
                    candidate%x      = cxy(2)
                    candidate%y      = cxy(3)
                    candidate%has_sh = .true.
                    call self%candidate_buffers(ithr)%append(i,candidate)
                enddo
            enddo
            !$omp end parallel do
        else
            ! fill the table
            !$omp parallel do default(shared) private(i,iptcl,ithr,o_prev,irot,iproj,is,istate,iref,candidate)&
            !$omp proc_bind(close) schedule(static)
            do i = i_first, i_last
                iptcl = self%pinds(i)
                ithr  = omp_get_thread_num() + 1
                ! identify shifts using the previously assigned best reference
                call self%b_ptr%spproj_field%get_ori(iptcl, o_prev)   ! previous ori
                irot  = self%b_ptr%pftc%get_roind(360.-o_prev%e3get())          ! in-plane angle index
                iproj = self%b_ptr%eulspace%find_closest_proj(o_prev) ! previous projection direction
                do is = 1, self%nstates
                    istate = self%ssinds(is)
                    iref   = (istate-1)*self%p_ptr%nspace + iproj
                    candidate%dist   = eulprob_dist_switch(&
                        &real(self%b_ptr%pftc%gen_corr_for_rot_8(iref, iptcl, irot)), self%p_ptr%cc_objfun)
                    candidate%iref   = iref
                    candidate%inpl   = irot
                    candidate%x      = 0.
                    candidate%y      = 0.
                    candidate%has_sh = .true.
                    call self%candidate_buffers(ithr)%append(i,candidate)
                enddo
            enddo
            !$omp end parallel do
        endif
        do ithr = 1,nthr_glob
            call grad_shsrch_obj(ithr)%kill
        end do
        if( self%table_is_open ) call self%flush_candidate_buffers
        call o_prev%kill
    end subroutine fill_tab_state_only_range

    ! ptcl -> (proj, state) assignment using calibrated likelihood weights
    subroutine ref_assign( self )
        class(eul_prob_tab), intent(inout) :: self
        integer, allocatable :: stab_inds(:,:), inds_sorted(:), iref_dist_inds(:)
        integer, allocatable :: greedy_state(:), active_refs(:), active_inds(:), score_mode(:)
        real,    allocatable :: score_work(:,:), score_min(:), score_spread(:)
        real,    allocatable :: iref_dist(:), dists_sorted(:)
        real,    allocatable :: state_projs_athres(:), active_dists(:), active_dists_sorted(:)
        logical, allocatable :: ptcl_avail(:)
        integer :: i, iref, assigned_iref, assigned_ptcl, istate, si, active_idx, nactive, ithr
        integer :: alloc_stat
        real    :: projs_athres, dist_tmp, corr_tmp
        character(len=256) :: alloc_msg
        allocate(stab_inds(self%nptcls, self%nrefs), score_work(self%nptcls,nthr_glob),&
            &score_min(self%nptcls), score_spread(self%nptcls), score_mode(self%nptcls),&
            &inds_sorted(self%nrefs), iref_dist_inds(self%nrefs), greedy_state(self%nptcls),&
            &active_refs(self%nrefs), active_inds(self%nrefs), iref_dist(self%nrefs),&
            &dists_sorted(self%nrefs), state_projs_athres(self%p_ptr%nstates), active_dists(self%nrefs),&
            &active_dists_sorted(self%nrefs), ptcl_avail(self%nptcls), stat=alloc_stat, errmsg=alloc_msg)
        if( alloc_stat /= 0 )then
            write(logfhandle,*) 'eul_prob_tab%ref_assign allocation failed for nptcls/nrefs: ', self%nptcls, self%nrefs
            write(logfhandle,*) trim(alloc_msg)
            THROW_HARD('failed allocating probability assignment work arrays')
        endif
        call prepare_ref_score_vectors
        ! Sort each reference using one scratch score column per worker thread.
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iref,i,ithr)
        do iref = 1, self%nrefs
            ithr = omp_get_thread_num() + 1
            do i = 1,self%nptcls
                stab_inds(i,iref) = i
                score_work(i,ithr) = ref_score(iref, i)
            enddo
            call hpsort(score_work(:,ithr), stab_inds(:,iref))
        enddo
        !$omp end parallel do
        projs_athres = 0.
        state_projs_athres = 0.
        do si = 1, self%nstates
            istate = self%ssinds(si)
            state_projs_athres(istate) = calc_athres(self%b_ptr%spproj_field, 'dist',&
                &self%p_ptr%prob_athres, state=istate)
            projs_athres = max(projs_athres, state_projs_athres(istate))
        enddo
        greedy_state = 0
        if( self%nstates > 1 )then
            call assign_greedy_state_labels()
            do si = 1,self%nstates
                call assign_refs_for_state(self%ssinds(si))
            enddo
        else
            call reset_ref_frontier()
            do while( any(ptcl_avail) )
                call sample_likelihood_dist(self%nrefs, frontier_ref_dist, likelihood_nsample(self%nrefs, projs_athres),&
                    &dist_tmp, corr_tmp, assigned_iref, dists_sorted, inds_sorted)
                call assign_current_ref()
                do iref = 1,self%nrefs
                    call advance_ref_head(iref, 0)
                enddo
            enddo
        endif
        if( allocated(stab_inds) ) deallocate(stab_inds, inds_sorted, iref_dist_inds, greedy_state, active_refs, active_inds,&
            &score_work, score_min, score_spread, score_mode, iref_dist, dists_sorted, state_projs_athres, active_dists,&
            &active_dists_sorted, ptcl_avail)

    contains

        subroutine prepare_ref_score_vectors()
            real    :: min_dist, max_dist, dist_val, spread, invalid_dist
            integer :: i_loc, iref_loc, nvalid, nflat, nbad
            invalid_dist = 0.1 * huge(invalid_dist)
            nflat = 0
            nbad  = 0
            score_min    = 0.
            score_spread = 0.
            score_mode   = 0
            !$omp parallel do default(shared) proc_bind(close) schedule(static)&
            !$omp private(i_loc,iref_loc,min_dist,max_dist,dist_val,spread,nvalid) reduction(+:nflat,nbad)
            do i_loc = 1, self%nptcls
                min_dist = huge(min_dist)
                max_dist = -huge(max_dist)
                nvalid   = 0
                do iref_loc = 1,self%nrefs
                    dist_val = self%loc_tab(iref_loc,i_loc)%dist
                    if( ieee_is_finite(dist_val) .and. dist_val < invalid_dist )then
                        min_dist = min(min_dist, dist_val)
                        max_dist = max(max_dist, dist_val)
                        nvalid   = nvalid + 1
                    endif
                enddo
                if( nvalid == 0 )then
                    nbad = nbad + 1
                    cycle
                endif
                spread = max_dist - min_dist
                score_min(i_loc)    = min_dist
                score_spread(i_loc) = spread
                if( spread <= 0. )then
                    score_mode(i_loc) = 1
                    nflat = nflat + 1
                else
                    score_mode(i_loc) = 2
                endif
            enddo
            !$omp end parallel do
            if( nbad == self%nptcls ) THROW_HARD('all probability-table reference distances are invalid')
            if( nflat == self%nptcls ) THROW_HARD('all probability-table reference distances are flat')
            if( nbad > 0 ) write(logfhandle,*) 'WARNING: particles with invalid probability-table distances: ', nbad
            if( nflat > 0 ) write(logfhandle,*) 'WARNING: particles with flat probability-table reference distances: ', nflat
        end subroutine prepare_ref_score_vectors

        real function ref_score( iref_loc, iptcl_loc ) result(score)
            integer, intent(in) :: iref_loc, iptcl_loc
            real :: dist_val, invalid_dist
            invalid_dist = 0.1 * huge(invalid_dist)
            dist_val = self%loc_tab(iref_loc,iptcl_loc)%dist
            score = 1.
            if( .not.(ieee_is_finite(dist_val) .and. dist_val < invalid_dist) )then
                score = invalid_dist
                return
            endif
            score = dist_val
        end function ref_score

        integer function likelihood_nsample( n, athres_ub_in ) result(nsample)
            integer, intent(in) :: n
            real,    intent(in) :: athres_ub_in
            real :: athres_ub
            athres_ub = min(self%p_ptr%prob_athres, athres_ub_in)
            nsample = min(n, max(1, int(athres_ub * real(n) / 180.)))
        end function likelihood_nsample

        real function frontier_ref_dist( iref_loc ) result(dist)
            integer, intent(in) :: iref_loc
            dist = iref_dist(iref_loc)
        end function frontier_ref_dist

        subroutine reset_ref_frontier()
            do iref = 1,self%nrefs
                iref_dist(iref) = ref_score(iref, stab_inds(1,iref))
            enddo
            iref_dist_inds = 1
            ptcl_avail     = .true.
        end subroutine reset_ref_frontier

        subroutine advance_ref_head( iref_current, state_filter )
            integer, intent(in) :: iref_current, state_filter
            integer :: iptcl_current
            do while( iref_dist_inds(iref_current) <= self%nptcls )
                iptcl_current = stab_inds(iref_dist_inds(iref_current), iref_current)
                if( ptcl_avail(iptcl_current) .and. &
                    &(state_filter == 0 .or. greedy_state(iptcl_current) == state_filter) )then
                    iref_dist(iref_current) = ref_score(iref_current, iptcl_current)
                    return
                endif
                iref_dist_inds(iref_current) = iref_dist_inds(iref_current) + 1
            enddo
            iref_dist(iref_current) = huge(iref_dist(iref_current))
        end subroutine advance_ref_head

        subroutine assign_current_ref()
            assigned_ptcl = stab_inds(iref_dist_inds(assigned_iref), assigned_iref)
            ptcl_avail(assigned_ptcl)     = .false.
            call self%assign_candidate(assigned_ptcl, self%loc_tab(assigned_iref,assigned_ptcl))
            call materialize_seed_shift(self%assgn_map(assigned_ptcl), self%seed_shifts(:,assigned_ptcl),&
                &self%seed_has_sh(assigned_ptcl), self%p_ptr%l_doshift, self%seed_nrots)
            self%assgn_map(assigned_ptcl)%frac = 100.
        end subroutine assign_current_ref

        subroutine assign_greedy_state_labels()
            ! Multi-state state labelling is deterministic (argmin distance) for both weighting
            ! schemes; probabilistic exploration is confined to the within-state projection
            ! assignment (assign_refs_for_state). Only refine=prob_state samples the state label.
            call reset_ref_frontier()
            do while( any(ptcl_avail) )
                assigned_iref = minloc(iref_dist, dim=1)
                assigned_ptcl  = stab_inds(iref_dist_inds(assigned_iref), assigned_iref)
                greedy_state(assigned_ptcl) = self%ref_state(assigned_iref)
                ptcl_avail(assigned_ptcl)  = .false.
                do iref = 1,self%nrefs
                    call advance_ref_head(iref, 0)
                enddo
            enddo
        end subroutine assign_greedy_state_labels

        subroutine gather_active_refs( state_filter )
            integer, intent(in) :: state_filter
            nactive = 0
            do iref = 1,self%nrefs
                if( self%ref_state(iref) /= state_filter ) cycle
                if( iref_dist_inds(iref) > self%nptcls ) cycle
                nactive = nactive + 1
                active_refs(nactive)  = iref
                active_dists(nactive) = iref_dist(iref)
            enddo
        end subroutine gather_active_refs

        subroutine assign_refs_for_state( state_filter )
            integer, intent(in) :: state_filter
            call reset_ref_frontier()
            ptcl_avail = greedy_state == state_filter
            do iref = 1,self%nrefs
                if( self%ref_state(iref) == state_filter ) call advance_ref_head(iref, state_filter)
            enddo
            do while( any(ptcl_avail) )
                call gather_active_refs(state_filter)
                if( nactive == 0 ) THROW_HARD('no active refs left in multi-state probability assignment')
                if( nactive == 1 )then
                    active_idx = 1
                else
                    call sample_likelihood_dist(nactive, active_ref_dist, likelihood_nsample(nactive, state_projs_athres(state_filter)),&
                        &dist_tmp, corr_tmp, active_idx, active_dists_sorted, active_inds)
                endif
                assigned_iref = active_refs(active_idx)
                call assign_current_ref()
                do iref = 1,self%nrefs
                    if( self%ref_state(iref) == state_filter ) call advance_ref_head(iref, state_filter)
                enddo
            enddo
        end subroutine assign_refs_for_state

        real function active_ref_dist( active_loc ) result(dist)
            integer, intent(in) :: active_loc
            dist = active_dists(active_loc)
        end function active_ref_dist

    end subroutine ref_assign

    ! ptcl -> state (using assigned iproj or previous iproj) assignment
    subroutine state_assign( self )
        class(eul_prob_tab), intent(inout) :: self
        integer :: i, istate, assigned_istate, assigned_ptcl, state_dist_inds(self%nstates),&
                    &stab_inds(self%nptcls, self%nstates), inds_sorted(self%nstates)
        real    :: sorted_tab(self%nptcls, self%nstates), state_dist(self%nstates), state_dists_sorted(self%nstates)
        real    :: dist_tmp, corr_tmp
        logical :: ptcl_avail(self%nptcls)
        if( self%nstates == 1 )then
            do i = 1,self%nptcls
                call self%assign_candidate(i, self%state_tab(1,i))
                self%assgn_map(i)%frac = 100.
            enddo
            return
        endif
        ! sorting each columns
        sorted_tab = transpose(self%state_tab%dist)
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(istate,i)
        do istate = 1, self%nstates
            stab_inds(:,istate) = (/(i,i=1,self%nptcls)/)
            call hpsort(sorted_tab(:,istate), stab_inds(:,istate))
        enddo
        !$omp end parallel do
        ! first row is the current best state distribution
        state_dist_inds = 1
        state_dist      = sorted_tab(1,:)
        ptcl_avail      = .true.
        do while( any(ptcl_avail) )
            call sample_likelihood_dist(self%nstates, state_frontier_dist, self%nstates, dist_tmp, corr_tmp,&
                &assigned_istate, state_dists_sorted, inds_sorted)
            assigned_ptcl   = stab_inds(state_dist_inds(assigned_istate), assigned_istate)
            ptcl_avail(assigned_ptcl)     = .false.
            call self%assign_candidate(assigned_ptcl, self%state_tab(assigned_istate,assigned_ptcl))
            self%assgn_map(assigned_ptcl)%frac = 100.
            ! update the state_dist and state_dist_inds
            do istate = 1, self%nstates
                call advance_state_head(istate)
            enddo
        enddo
    contains

        subroutine advance_state_head( state_loc )
            integer, intent(in) :: state_loc
            do while( state_dist_inds(state_loc) <= self%nptcls )
                if( ptcl_avail(stab_inds(state_dist_inds(state_loc), state_loc)) )then
                    state_dist(state_loc) = sorted_tab(state_dist_inds(state_loc), state_loc)
                    return
                endif
                state_dist_inds(state_loc) = state_dist_inds(state_loc) + 1
            enddo
            state_dist(state_loc) = huge(state_dist(state_loc))
        end subroutine advance_state_head

        real function state_frontier_dist( state_loc ) result(dist)
            integer, intent(in) :: state_loc
            dist = state_dist(state_loc)
        end function state_frontier_dist

    end subroutine state_assign

    pure integer function ref_state( self, iref ) result( state )
        class(eul_prob_tab), intent(in) :: self
        integer,             intent(in) :: iref
        integer :: active_rank
        active_rank = (iref - 1) / self%p_ptr%nspace + 1
        state = self%ssinds(active_rank)
    end function ref_state

    pure integer function ref_proj( self, iref ) result( proj )
        class(eul_prob_tab), intent(in) :: self
        integer,             intent(in) :: iref
        proj = modulo(iref - 1, self%p_ptr%nspace) + 1
    end function ref_proj

    pure integer function ref_full( self, iref ) result( full_ref )
        class(eul_prob_tab), intent(in) :: self
        integer,             intent(in) :: iref
        full_ref = (self%ref_state(iref) - 1) * self%p_ptr%nspace + self%ref_proj(iref)
    end function ref_full

    pure integer function full_to_compact_ref( self, full_ref ) result( iref )
        class(eul_prob_tab), intent(in) :: self
        integer,             intent(in) :: full_ref
        integer :: state, proj, active_rank
        iref = 0
        if( full_ref < 1 .or. full_ref > self%p_ptr%nspace * self%p_ptr%nstates ) return
        state = (full_ref - 1) / self%p_ptr%nspace + 1
        proj  = modulo(full_ref - 1, self%p_ptr%nspace) + 1
        active_rank = self%state_to_active_rank(state)
        if( active_rank < 1 ) return
        iref = (active_rank - 1) * self%p_ptr%nspace + proj
    end function full_to_compact_ref

    subroutine assign_candidate( self, iptcl_map, candidate )
        class(eul_prob_tab),   intent(inout) :: self
        integer,               intent(in)    :: iptcl_map
        type(prob_candidate),  intent(in)    :: candidate
        integer :: state, proj
        if( candidate%iref < 1 .or. candidate%iref > self%p_ptr%nspace * self%p_ptr%nstates )then
            THROW_HARD('candidate reference is out of range in assign_candidate')
        endif
        state = (candidate%iref - 1) / self%p_ptr%nspace + 1
        proj  = modulo(candidate%iref - 1, self%p_ptr%nspace) + 1
        self%assgn_map(iptcl_map)%pind   = self%pinds(iptcl_map)
        self%assgn_map(iptcl_map)%istate = state
        self%assgn_map(iptcl_map)%iproj  = proj
        self%assgn_map(iptcl_map)%inpl   = candidate%inpl
        self%assgn_map(iptcl_map)%dist   = candidate%dist
        self%assgn_map(iptcl_map)%x      = candidate%x
        self%assgn_map(iptcl_map)%y      = candidate%y
        self%assgn_map(iptcl_map)%has_sh = candidate%has_sh
    end subroutine assign_candidate

    ! FILE IO

    subroutine begin_write( self, binfname )
        class(eul_prob_tab), intent(inout) :: self
        class(string),       intent(in)    :: binfname
        integer :: io_stat
        integer(int64) :: file_header(4)
        integer(int64) :: addr
        if( self%table_is_open ) THROW_HARD('probability candidate stream is already open')
        if( .not. allocated(self%candidate_buffers) ) THROW_HARD('candidate stream requested without worker buffers')
        file_header = [int(self%nrefs,int64),int(self%nptcls,int64),0_int64,0_int64]
        call fopen(self%table_unit,binfname,access='STREAM',action='WRITE',status='REPLACE',iostat=io_stat)
        call fileiochk('simple_eul_prob_tab; begin_write; file: '//binfname%to_char(),io_stat)
        self%table_is_open = .true.
        write(self%table_unit,pos=1) file_header
        addr = sizeof(file_header) + 1
        write(self%table_unit,pos=addr) self%pinds
        addr = addr + sizeof(self%pinds)
        call write_seed_shift_table(self%table_unit,addr,self%seed_nrots,self%seed_shifts,self%seed_has_sh)
        self%table_addr    = addr
        self%table_nnz     = 0
        self%table_nchunks = 0
    end subroutine begin_write

    subroutine flush_candidate_buffers( self )
        class(eul_prob_tab), intent(inout) :: self
        integer :: ithr, nused
        if( .not. self%table_is_open ) return
        do ithr = 1,size(self%candidate_buffers)
            nused = self%candidate_buffers(ithr)%nused
            if( nused < 1 ) cycle
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

    ! Finalises the one-file-per-part candidate stream.
    subroutine write_tab( self, binfname )
        class(eul_prob_tab), intent(inout) :: self
        class(string),       intent(in)    :: binfname
        integer(int64) :: file_header(4), addr
        if( .not. self%table_is_open ) call self%begin_write(binfname)
        call self%flush_candidate_buffers
        if( self%table_nnz < 1 ) THROW_HARD('eul_prob_tab%write_tab; empty candidate stream')
        file_header = [int(self%nrefs,int64),int(self%nptcls,int64),self%table_nnz,int(self%table_nchunks,int64)]
        write(self%table_unit,pos=1) file_header
        addr = sizeof(file_header) + 1
        write(self%table_unit,pos=addr) self%pinds
        addr = addr + sizeof(self%pinds)
        call write_seed_shift_table(self%table_unit,addr,self%seed_nrots,self%seed_shifts,self%seed_has_sh)
        call fclose(self%table_unit)
        self%table_unit = -1
        self%table_is_open = .false.
        self%table_addr = 1
    end subroutine write_tab

    ! read the partition-wise dist value binary file to global reg object's dist value table
    subroutine read_tab_to_glob( self, binfname )
        class(eul_prob_tab), intent(inout) :: self
        class(string),       intent(in)    :: binfname
        type(prob_candidate), allocatable  :: candidates_loc(:)
        real,                allocatable   :: seed_shifts_loc(:,:)
        logical,             allocatable   :: seed_has_sh_loc(:)
        integer, allocatable :: pind2glob(:), pinds_loc(:), particle_indices(:), candidate_counts(:)
        integer :: funit, io_stat, nptcls_loc, nrefs_loc, nchunks
        integer(int64) :: file_header(4), nnz, nnz_read
        integer(int64) :: addr, chunk_indices_addr, chunk_candidates_addr
        integer :: i_loc, i_glob, pind, max_pind, seed_nrots_loc, ichunk, chunk_n
        integer :: first, last, nread, j, ri
        if( file_exists(binfname) )then
            call fopen(funit,binfname,access='STREAM',action='READ',status='OLD', iostat=io_stat)
            call fileiochk('simple_eul_prob_tab; read_tab_to_glob; file: '//binfname%to_char(), io_stat)
        else
            THROW_HARD( 'corr/rot files of partitions should be ready! ' )
        endif
        read(unit=funit,pos=1) file_header
        nrefs_loc  = int(file_header(1))
        nptcls_loc = int(file_header(2))
        nnz        = file_header(3)
        nchunks    = int(file_header(4))
        if( nrefs_loc .ne. self%nrefs ) THROW_HARD( 'nrefs should be the same as nrefs in this partition file!' )
        allocate(pinds_loc(nptcls_loc), seed_shifts_loc(2,nptcls_loc), seed_has_sh_loc(nptcls_loc))
        addr = sizeof(file_header) + 1
        read(funit,pos=addr,iostat=io_stat) pinds_loc
        call fileiochk('simple_eul_prob_tab; read_tab_to_glob pinds; file: '//binfname%to_char(), io_stat)
        addr = addr + sizeof(pinds_loc)
        call read_seed_shift_table(funit, addr, seed_nrots_loc, seed_shifts_loc, seed_has_sh_loc)
        if( self%seed_nrots == 0 ) self%seed_nrots = seed_nrots_loc
        if( self%seed_nrots /= seed_nrots_loc ) THROW_HARD('seed_nrots mismatch in eul_prob_tab%read_tab_to_glob')
        call build_pind_lookup(self%pinds, pinds_loc, pind2glob, max_pind)
        do i_loc = 1,nptcls_loc
            pind = pinds_loc(i_loc)
            if( pind >= 1 .and. pind <= max_pind )then
                i_glob = pind2glob(pind)
                if( i_glob > 0 )then
                    self%seed_shifts(:,i_glob) = seed_shifts_loc(:,i_loc)
                    self%seed_has_sh(i_glob)   = seed_has_sh_loc(i_loc)
                endif
            endif
        enddo
        allocate(particle_indices(PROB_TAB_IO_CHUNK), candidates_loc(PROB_TAB_IO_CHUNK))
        allocate(candidate_counts(nptcls_loc), source=0)
        nnz_read = 0
        do ichunk = 1,nchunks
            read(funit,pos=addr) chunk_n
            addr = addr + sizeof(chunk_n)
            if( chunk_n < 1 ) THROW_HARD('invalid empty candidate chunk')
            chunk_indices_addr    = addr
            chunk_candidates_addr = chunk_indices_addr + chunk_n * sizeof(chunk_n)
            do first = 1,chunk_n,PROB_TAB_IO_CHUNK
                last  = min(chunk_n,first+PROB_TAB_IO_CHUNK-1)
                nread = last-first+1
                read(funit,pos=chunk_indices_addr+(first-1)*sizeof(chunk_n)) particle_indices(1:nread)
                read(funit,pos=chunk_candidates_addr+(first-1)*sizeof(candidates_loc(1))) candidates_loc(1:nread)
                do j = 1,nread
                    i_loc = particle_indices(j)
                    if( i_loc < 1 .or. i_loc > nptcls_loc ) THROW_HARD('invalid particle index in candidate stream')
                    candidate_counts(i_loc) = candidate_counts(i_loc) + 1
                    pind = pinds_loc(i_loc)
                    if( pind < 1 .or. pind > max_pind ) cycle
                    i_glob = pind2glob(pind)
                    if( i_glob < 1 ) cycle
                    ri = self%full_to_compact_ref(candidates_loc(j)%iref)
                    if( ri < 1 ) THROW_HARD('invalid reference in dense candidate stream')
                    self%loc_tab(ri,i_glob) = candidates_loc(j)
                enddo
            enddo
            addr = chunk_candidates_addr + chunk_n * sizeof(candidates_loc(1))
            nnz_read = nnz_read + int(chunk_n,int64)
        enddo
        if( nnz_read /= nnz ) THROW_HARD('candidate stream count mismatch')
        if( any(candidate_counts /= nrefs_loc) ) THROW_HARD('dense candidate stream is incomplete')
        call fclose(funit)
        deallocate(candidates_loc,particle_indices,candidate_counts,pinds_loc,seed_shifts_loc,seed_has_sh_loc,pind2glob)
    end subroutine read_tab_to_glob

    subroutine write_state_tab( self, binfname )
        class(eul_prob_tab), intent(inout) :: self
        class(string),       intent(in)    :: binfname
        call self%write_tab(binfname)
    end subroutine write_state_tab

    subroutine read_state_tab( self, binfname )
        class(eul_prob_tab), intent(inout) :: self
        class(string),       intent(in)    :: binfname
        type(prob_candidate), allocatable :: candidates_loc(:)
        real, allocatable :: seed_shifts_loc(:,:)
        logical, allocatable :: seed_has_sh_loc(:)
        integer, allocatable :: pind2glob(:), pinds_loc(:), particle_indices(:), candidate_counts(:)
        integer :: funit, io_stat, nptcls_loc, nrefs_loc, nchunks, seed_nrots_loc
        integer :: i_loc, i_glob, pind, max_pind, ichunk, chunk_n, first, last, nread, j, state_rank
        integer(int64) :: file_header(4), nnz, nnz_read, addr, chunk_indices_addr, chunk_candidates_addr
        if( .not. file_exists(binfname) )then
            THROW_HARD('file '//binfname%to_char()//' does not exists!')
        else
            call fopen(funit,binfname,access='STREAM',action='READ',status='OLD', iostat=io_stat)
        end if
        call fileiochk('simple_eul_prob_tab; read_state_tab; file: '//binfname%to_char(), io_stat)
        read(unit=funit,pos=1) file_header
        nrefs_loc  = int(file_header(1))
        nptcls_loc = int(file_header(2))
        nnz        = file_header(3)
        nchunks    = int(file_header(4))
        if( nrefs_loc /= self%nrefs ) THROW_HARD('reference count mismatch in read_state_tab')
        allocate(pinds_loc(nptcls_loc),seed_shifts_loc(2,nptcls_loc),seed_has_sh_loc(nptcls_loc))
        addr = sizeof(file_header) + 1
        read(funit,pos=addr) pinds_loc
        addr = addr + sizeof(pinds_loc)
        call read_seed_shift_table(funit,addr,seed_nrots_loc,seed_shifts_loc,seed_has_sh_loc)
        call build_pind_lookup(self%pinds,pinds_loc,pind2glob,max_pind)
        allocate(particle_indices(PROB_TAB_IO_CHUNK),candidates_loc(PROB_TAB_IO_CHUNK))
        allocate(candidate_counts(nptcls_loc),source=0)
        nnz_read = 0
        do ichunk = 1,nchunks
            read(funit,pos=addr) chunk_n
            addr = addr + sizeof(chunk_n)
            chunk_indices_addr    = addr
            chunk_candidates_addr = chunk_indices_addr + chunk_n*sizeof(chunk_n)
            do first = 1,chunk_n,PROB_TAB_IO_CHUNK
                last  = min(chunk_n,first+PROB_TAB_IO_CHUNK-1)
                nread = last-first+1
                read(funit,pos=chunk_indices_addr+(first-1)*sizeof(chunk_n)) particle_indices(1:nread)
                read(funit,pos=chunk_candidates_addr+(first-1)*sizeof(candidates_loc(1))) candidates_loc(1:nread)
                do j = 1,nread
                    i_loc = particle_indices(j)
                    if( i_loc < 1 .or. i_loc > nptcls_loc ) THROW_HARD('invalid particle index in state stream')
                    candidate_counts(i_loc) = candidate_counts(i_loc) + 1
                    pind = pinds_loc(i_loc)
                    if( pind < 1 .or. pind > max_pind ) cycle
                    i_glob = pind2glob(pind)
                    if( i_glob < 1 ) cycle
                    state_rank = self%state_to_active_rank((candidates_loc(j)%iref-1)/self%p_ptr%nspace+1)
                    if( state_rank < 1 ) THROW_HARD('invalid state in state candidate stream')
                    self%state_tab(state_rank,i_glob) = candidates_loc(j)
                enddo
            enddo
            addr = chunk_candidates_addr + chunk_n*sizeof(candidates_loc(1))
            nnz_read = nnz_read + int(chunk_n,int64)
        enddo
        if( nnz_read /= nnz ) THROW_HARD('state candidate stream count mismatch')
        if( any(candidate_counts /= self%nstates) ) THROW_HARD('state candidate stream is incomplete')
        call fclose(funit)
        deallocate(candidates_loc,particle_indices,candidate_counts,pinds_loc,seed_shifts_loc,seed_has_sh_loc,pind2glob)
    end subroutine read_state_tab

    ! write a global assignment map to binary file
    subroutine write_assignment( self, binfname )
        class(eul_prob_tab), intent(in) :: self
        class(string),       intent(in) :: binfname
        integer :: funit, io_stat, headsz
        headsz = sizeof(self%nptcls)
        call fopen(funit,binfname,access='STREAM',action='WRITE',status='REPLACE', iostat=io_stat)
        write(unit=funit,pos=1)          self%nptcls
        write(unit=funit,pos=headsz + 1) self%assgn_map
        call fclose(funit)
    end subroutine write_assignment

    ! Writes the dense calibrated likelihood support used by a later finer-grid
    ! prob_neigh stage.  The artifact is deliberately a truncated support, not
    ! a claim to represent the full joint posterior over pose.
    subroutine write_prior3D( self, binfname )
        class(eul_prob_tab), intent(in) :: self
        class(string),       intent(in) :: binfname
        type(prior3d_writer) :: writer
        integer, allocatable :: cand_iref(:), cand_state(:), cand_proj(:), cand_inpl(:), cand_has_sh(:)
        integer, allocatable :: out_state(:), out_proj(:), out_inpl(:), out_has_sh(:)
        real,    allocatable :: cand_dist(:), work_dist(:), cand_weight(:), out_dist(:), out_weight(:)
        real,    allocatable :: out_x(:), out_y(:), out_euls(:,:)
        real :: dist_thres, dmin, norm, euls(3)
        integer :: i, k, ncand, kmax, k_this, npeaks_detected, pick, iproj
        if( .not. allocated(self%loc_tab) )then
            THROW_WARN('write_prior3D requires a dense probability table; no artifact written')
            return
        endif
        if( self%nptcls < 1 .or. self%nrefs < 1 ) return
        kmax = min(self%nrefs, max(PRIOR3D_NREMAP, min(128, max(1, 8*self%p_ptr%npeaks_inpl))))
        allocate(cand_iref(self%nrefs), cand_state(self%nrefs), cand_proj(self%nrefs), cand_inpl(self%nrefs),&
            &cand_has_sh(self%nrefs), cand_dist(self%nrefs), work_dist(self%nrefs), cand_weight(self%nrefs))
        allocate(out_state(kmax), out_proj(kmax), out_inpl(kmax), out_has_sh(kmax))
        allocate(out_dist(kmax), out_weight(kmax), out_x(kmax), out_y(kmax), out_euls(3,kmax))
        call writer%open(binfname%to_char(), self%nptcls, self%p_ptr%nstates, self%p_ptr%nspace,&
            &self%p_ptr%nspace_sub, kmax, self%p_ptr%pgrp, self%pinds)
        do i = 1, self%nptcls
            ncand = 0
            do k = 1, self%nrefs
                if( self%loc_tab(k,i)%inpl < 1 ) cycle
                if( .not. ieee_is_finite(self%loc_tab(k,i)%dist) ) cycle
                ncand = ncand + 1
                cand_iref(ncand) = k
                cand_state(ncand) = self%ref_state(k)
                cand_proj(ncand) = self%ref_proj(k)
                cand_inpl(ncand) = self%loc_tab(k,i)%inpl
                cand_has_sh(ncand) = merge(1, 0, self%loc_tab(k,i)%has_sh)
                cand_dist(ncand) = self%loc_tab(k,i)%dist
            enddo
            if( ncand < 1 )then
                call writer%write_row(self%pinds(i), 0, 0., out_state, out_proj, out_inpl, out_has_sh,&
                    &out_dist, out_weight, out_x, out_y, out_euls)
                cycle
            endif
            work_dist(1:ncand) = cand_dist(1:ncand)
            call detect_peak_thres_fdr(ncand, work_dist(1:ncand), PRIOR3D_FDR_Q, 1, min(kmax,ncand),&
                &dist_thres, npeaks_detected, lower_tail=.true.)
            k_this = min(kmax, ncand, max(1, npeaks_detected))
            dmin = minval(cand_dist(1:ncand))
            cand_weight(1:ncand) = exp(-(cand_dist(1:ncand)-dmin))
            norm = sum(cand_weight(1:ncand))
            if( norm > 0. ) cand_weight(1:ncand) = cand_weight(1:ncand) / norm
            out_state = 0
            out_proj = 0
            out_inpl = 0
            out_has_sh = 0
            out_dist = 0.
            out_weight = 0.
            out_x = 0.
            out_y = 0.
            out_euls = 0.
            do k = 1, k_this
                pick = minloc(work_dist(1:ncand), dim=1)
                work_dist(pick) = huge(1.0)
                out_state(k) = cand_state(pick)
                out_proj(k) = cand_proj(pick)
                out_inpl(k) = cand_inpl(pick)
                out_has_sh(k) = cand_has_sh(pick)
                out_dist(k) = cand_dist(pick)
                out_weight(k) = cand_weight(pick)
                out_x(k) = self%loc_tab(cand_iref(pick),i)%x
                out_y(k) = self%loc_tab(cand_iref(pick),i)%y
                iproj = out_proj(k)
                euls = self%b_ptr%eulspace%get_euler(iproj)
                out_euls(:,k) = euls
            enddo
            call writer%write_row(self%pinds(i), k_this, sum(out_weight(1:k_this)), out_state, out_proj, out_inpl, out_has_sh,&
                &out_dist, out_weight, out_x, out_y, out_euls)
        enddo
        call writer%close
        call writer%kill
        deallocate(cand_iref, cand_state, cand_proj, cand_inpl, cand_has_sh, cand_dist, work_dist, cand_weight,&
            &out_state, out_proj, out_inpl, out_has_sh, out_dist, out_weight, out_x, out_y, out_euls)
    end subroutine write_prior3D

    ! read from the global assignment map to local partition for shift search and further refinement
    subroutine read_assignment( self, binfname )
        class(eul_prob_tab), intent(inout) :: self
        class(string),       intent(in)    :: binfname
        type(ptcl_ref),      allocatable   :: assgn_glob(:)
        integer, allocatable :: pind2glob(:), pinds_glob(:)
        integer :: funit, io_stat, nptcls_glob, headsz, i_loc, i_glob, pind, max_pind
        headsz = sizeof(nptcls_glob)
        if( .not. file_exists(binfname) )then
            THROW_HARD('file '//binfname%to_char()//' does not exists!')
        else
            call fopen(funit,binfname,access='STREAM',action='READ',status='OLD', iostat=io_stat)
        end if
        call fileiochk('simple_eul_prob_tab; read_assignment; file: '//binfname%to_char(), io_stat)
        read(unit=funit,pos=1) nptcls_glob
        allocate(assgn_glob(nptcls_glob))
        read(unit=funit,pos=headsz + 1) assgn_glob
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

    ! DESTRUCTOR

    subroutine kill( self )
        class(eul_prob_tab), intent(inout) :: self
        integer :: ithr
        if( self%table_is_open ) call fclose(self%table_unit)
        self%table_unit    = -1
        self%table_is_open = .false.
        self%table_addr    = 1
        self%table_nnz     = 0
        self%table_nchunks = 0
        if( allocated(self%candidate_buffers) )then
            do ithr = 1,size(self%candidate_buffers)
                call self%candidate_buffers(ithr)%kill
            enddo
            deallocate(self%candidate_buffers)
        endif
        if( allocated(self%loc_tab)      ) deallocate(self%loc_tab)
        if( allocated(self%state_tab)    ) deallocate(self%state_tab)
        if( allocated(self%assgn_map)    ) deallocate(self%assgn_map)
        if( allocated(self%seed_shifts)  ) deallocate(self%seed_shifts)
        if( allocated(self%seed_has_sh)  ) deallocate(self%seed_has_sh)
        if( allocated(self%pinds)        ) deallocate(self%pinds)
        if( allocated(self%ssinds)       ) deallocate(self%ssinds)
        if( allocated(self%state_to_active_rank) ) deallocate(self%state_to_active_rank)
        if( allocated(self%state_exists) ) deallocate(self%state_exists)
        self%seed_nrots = 0
        self%b_ptr => null()
        self%p_ptr => null()
    end subroutine kill

end module simple_eul_prob_tab
