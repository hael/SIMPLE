!@descr: the core probability table routines used for probabilistic 3D search
module simple_eul_prob_tab
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_pftc_srch_api
use simple_builder,          only: builder
use simple_eul_prob_tab_utils, only: angle_sampling, build_pind_lookup, calc_athres, calc_num2sample,&
    &eulprob_dist_switch, materialize_seed_shift, read_seed_shift_table, sample_likelihood_dist,&
    &write_seed_shift_table
use simple_pftc_shsrch_grad, only: pftc_shsrch_grad
use simple_type_defs,        only: OBJFUN_EUCLID
implicit none

public :: eul_prob_tab
private
#include "simple_local_flags.inc"

integer, parameter :: PROB_TAB_IO_CHUNK = 1024

type :: eul_prob_tab
    class(builder),    pointer  :: b_ptr => null()
    class(parameters), pointer  :: p_ptr => null()
    type(ptcl_ref), allocatable :: loc_tab(:,:)    !< 2D search table (nspace*nstates, nptcls)
    type(ptcl_ref), allocatable :: state_tab(:,:)  !< 2D search table (nstates,        nptcls)
    type(ptcl_ref), allocatable :: assgn_map(:)    !< assignment map                  (nptcls)
    real,           allocatable :: seed_shifts(:,:) !< per-particle seeded shift       (2,nptcls)
    logical,        allocatable :: seed_has_sh(:)   !< per-particle seeded shift flag  (nptcls)
    integer                     :: seed_nrots = 0   !< rotation grid used for deferred seeded shifts
    integer,        allocatable :: pinds(:)        !< particle indices for processing
    integer,        allocatable :: ssinds(:)       !< non-empty state indices
    integer,        allocatable :: sinds(:)        !< non-empty state indices of each ref
    integer,        allocatable :: jinds(:)        !< non-empty proj  indices of each ref
    logical,        allocatable :: proj_exists(:,:)
    logical,        allocatable :: state_exists(:)
    integer                     :: nptcls          !< size of pinds array
    integer                     :: nstates         !< states number
    integer                     :: nrefs           !< reference number
    contains
    ! CONSTRUCTOR
    procedure :: new
    procedure :: new_assignment
    ! PARTITION-WISE PROCEDURES (used only by partition-wise eul_prob_tab objects)
    procedure :: fill_tab
    procedure :: fill_tab_range
    procedure :: fill_tab_state_only
    procedure :: fill_tab_state_only_range
    procedure :: write_tab
    procedure :: write_state_tab
    procedure :: read_assignment
    ! GLOBAL PROCEDURES (used only by the global eul_prob_tab object)
    procedure :: read_state_tab
    procedure :: read_tab_to_glob
    procedure :: ref_assign
    procedure :: write_assignment
    procedure :: state_assign
    ! DESTRUCTOR
    procedure :: kill
    ! PRIVATE
    procedure, private :: ref_normalize
    procedure, private :: state_normalize
end type eul_prob_tab

contains

    ! CONSTRUCTORS

    subroutine new( self, params, build, pinds, empty_okay, state_only )
        class(eul_prob_tab),       intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        class(builder),    target, intent(in)    :: build
        integer,                   intent(in)    :: pinds(:)
        logical, optional,         intent(in)    :: empty_okay
        logical, optional,         intent(in)    :: state_only
        integer, parameter :: MIN_POP = 5   ! ignoring cavgs with less than 5 particles
        integer :: i, iproj, iptcl, istate, si, ri
        real    :: x
        logical :: l_state_only
        call self%kill
        l_state_only = .false.
        if( present(state_only) ) l_state_only = state_only
        self%p_ptr => params
        self%b_ptr  => build
        self%nptcls       = size(pinds)
        self%state_exists = self%b_ptr%spproj_field%states_exist(self%p_ptr%nstates, thres=MIN_POP)
        self%nstates      = count(self%state_exists .eqv. .true.)
        ! In 3D, projection directions are defined by volume reprojections and
        ! should always be considered available for states that exist.
        allocate(self%proj_exists(self%p_ptr%nspace,self%p_ptr%nstates), source=.false.)
        do istate = 1,self%p_ptr%nstates
            if( self%state_exists(istate) ) self%proj_exists(:,istate) = .true.
        enddo
        self%nrefs = count(self%proj_exists .eqv. .true.)
        allocate(self%ssinds(self%nstates),self%jinds(self%nrefs),self%sinds(self%nrefs))
        si = 0
        ri = 0
        do istate = 1, self%p_ptr%nstates
            if( .not. self%state_exists(istate) )cycle
            si              = si + 1
            self%ssinds(si) = istate
            do iproj = 1,self%p_ptr%nspace
                ri             = ri + 1
                self%jinds(ri) = iproj
                self%sinds(ri) = istate
            enddo
        enddo
        allocate(self%pinds(self%nptcls), source=pinds)
        allocate(self%assgn_map(self%nptcls), self%state_tab(self%nstates,self%nptcls))
        if( .not. l_state_only ) allocate(self%loc_tab(self%nrefs,self%nptcls))
        allocate(self%seed_shifts(2,self%nptcls), source=0.)
        allocate(self%seed_has_sh(self%nptcls), source=.false.)
        !$omp parallel do default(shared) private(i,iptcl,si,ri) proc_bind(close) schedule(static)
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
            do si = 1,self%nstates
                self%state_tab(si,i)%pind   = iptcl
                self%state_tab(si,i)%istate = self%ssinds(si)
                self%state_tab(si,i)%iproj  = 0
                self%state_tab(si,i)%inpl   = 0
                self%state_tab(si,i)%dist   = huge(x)
                self%state_tab(si,i)%x      = 0.
                self%state_tab(si,i)%y      = 0.
                self%state_tab(si,i)%has_sh = .false.
            enddo
            if( allocated(self%loc_tab) )then
                do ri = 1, self%nrefs
                    self%loc_tab(ri,i)%pind   = iptcl
                    self%loc_tab(ri,i)%istate = self%sinds(ri)
                    self%loc_tab(ri,i)%iproj  = self%jinds(ri)
                    self%loc_tab(ri,i)%inpl   = 0
                    self%loc_tab(ri,i)%dist   = huge(x)
                    self%loc_tab(ri,i)%x      = 0.
                    self%loc_tab(ri,i)%y      = 0.
                    self%loc_tab(ri,i)%has_sh = .false.
                end do
            endif
        end do
        !$omp end parallel do
    end subroutine new

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
        real    :: lims(2,2), lims_init(2,2), inpl_athres(self%p_ptr%nstates)
        real    :: dists_inpl(self%b_ptr%pftc%get_nrots(),nthr_glob),&
            &dists_inpl_sorted(self%b_ptr%pftc%get_nrots(),nthr_glob), dists_refs(self%nrefs,nthr_glob)
        logical :: l_prob_objfun, l_sh_first, l_likelihood_inpl
        if( i_first < 1 .or. i_last > self%nptcls .or. i_last < i_first )then
            THROW_HARD('invalid particle range in eul_prob_tab%fill_tab_range')
        endif
        self%seed_nrots = self%b_ptr%pftc%get_nrots()
        l_prob_objfun   = (self%p_ptr%cc_objfun == OBJFUN_EUCLID)
        l_likelihood_inpl = trim(self%p_ptr%prob_assign) == 'likelihood' .and. trim(self%p_ptr%refine) /= 'prob_state'
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
                call record_ref_eval(i_loc, ri_loc, dist_loc, irot_loc, 0., 0., .false.)
                dists_refs(ri_loc,ithr_loc) = dist_loc
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
                call record_ref_eval(i_loc, ri_loc, dist_loc, irot_loc, 0., 0., .false.)
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
            istate_loc = self%sinds(ri_loc)
            iproj_loc  = self%jinds(ri_loc)
            full_ref   = (istate_loc-1)*self%p_ptr%nspace + iproj_loc
            if( l_likelihood_inpl )then
                call self%b_ptr%pftc%gen_prob_likelihood_objfun_val(full_ref, iptcl_loc, shift_xy,&
                    &inpl_likelihood_nsample(istate_loc), dist_loc, corr_loc, irot_loc,&
                    &dists_inpl_sorted(:,ithr_loc), inds_sorted(:,ithr_loc))
            else if( l_prob_objfun )then
                call self%b_ptr%pftc%gen_prob_objfun_val(full_ref, iptcl_loc, shift_xy,&
                    &inpl_athres(istate_loc), self%p_ptr%prob_athres, dist_loc, irot_loc,&
                    &dists_inpl_sorted(:,ithr_loc), inds_sorted(:,ithr_loc))
            else
                call self%b_ptr%pftc%gen_objfun_vals(full_ref, iptcl_loc, shift_xy, dists_inpl(:,ithr_loc))
                dists_inpl(:,ithr_loc) = eulprob_dist_switch(dists_inpl(:,ithr_loc), self%p_ptr%cc_objfun)
                irot_loc = angle_sampling(dists_inpl(:,ithr_loc), dists_inpl_sorted(:,ithr_loc),&
                    &inds_sorted(:,ithr_loc), inpl_athres(istate_loc), self%p_ptr%prob_athres)
                dist_loc = dists_inpl(irot_loc,ithr_loc)
            endif
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
                istate_loc = self%sinds(ri_loc)
                iproj_loc  = self%jinds(ri_loc)
                call grad_shsrch_obj(ithr_loc)%set_indices((istate_loc-1)*self%p_ptr%nspace + iproj_loc, iptcl_loc)
                irot_loc = self%loc_tab(ri_loc,i_loc)%inpl
                if( l_sh_first )then
                    refined_shift = grad_shsrch_obj(ithr_loc)%minimize(irot=irot_loc, sh_rot=.true., xy_in=shift_seed(2:3))
                else
                    refined_shift = grad_shsrch_obj(ithr_loc)%minimize(irot=irot_loc, sh_rot=.true.)
                endif
                if( irot_loc > 0 )then
                    call record_ref_eval(i_loc, ri_loc, eulprob_dist_switch(refined_shift(1), self%p_ptr%cc_objfun),&
                        &irot_loc, refined_shift(2), refined_shift(3), .true.)
                endif
            enddo
        end subroutine refine_best_refs

        subroutine record_ref_eval( i_loc, ri_loc, dist_loc, irot_loc, x_loc, y_loc, has_sh_loc )
            integer, intent(in) :: i_loc, ri_loc, irot_loc
            real,    intent(in) :: dist_loc, x_loc, y_loc
            logical, intent(in) :: has_sh_loc
            self%loc_tab(ri_loc,i_loc)%dist   = dist_loc
            self%loc_tab(ri_loc,i_loc)%inpl   = irot_loc
            self%loc_tab(ri_loc,i_loc)%x      = x_loc
            self%loc_tab(ri_loc,i_loc)%y      = y_loc
            self%loc_tab(ri_loc,i_loc)%has_sh = has_sh_loc
        end subroutine record_ref_eval

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
            !$omp parallel do default(shared) private(i,iptcl,ithr,o_prev,iproj,is,istate,irot,iref,cxy)&
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
                    self%state_tab(is,i)%dist   = eulprob_dist_switch(cxy(1), self%p_ptr%cc_objfun)
                    self%state_tab(is,i)%iproj  = iproj
                    self%state_tab(is,i)%inpl   = irot
                    self%state_tab(is,i)%x      = cxy(2)
                    self%state_tab(is,i)%y      = cxy(3)
                    self%state_tab(is,i)%has_sh = .true.
                enddo
            enddo
            !$omp end parallel do
        else
            ! fill the table
            !$omp parallel do default(shared) private(i,iptcl,o_prev,irot,iproj,is,istate,iref)&
            !$omp proc_bind(close) schedule(static)
            do i = i_first, i_last
                iptcl = self%pinds(i)
                ! identify shifts using the previously assigned best reference
                call self%b_ptr%spproj_field%get_ori(iptcl, o_prev)   ! previous ori
                irot  = self%b_ptr%pftc%get_roind(360.-o_prev%e3get())          ! in-plane angle index
                iproj = self%b_ptr%eulspace%find_closest_proj(o_prev) ! previous projection direction
                do is = 1, self%nstates
                    istate = self%ssinds(is)
                    iref   = (istate-1)*self%p_ptr%nspace + iproj
                    self%state_tab(is,i)%dist   = eulprob_dist_switch(&
                        &real(self%b_ptr%pftc%gen_corr_for_rot_8(iref, iptcl, irot)), self%p_ptr%cc_objfun)
                    self%state_tab(is,i)%iproj  = iproj
                    self%state_tab(is,i)%inpl   = irot
                    self%state_tab(is,i)%x      = 0.
                    self%state_tab(is,i)%y      = 0.
                    self%state_tab(is,i)%has_sh = .true.
                enddo
            enddo
            !$omp end parallel do
        endif
        do ithr = 1,nthr_glob
            call grad_shsrch_obj(ithr)%kill
        end do
        call o_prev%kill
    end subroutine fill_tab_state_only_range

    ! Legacy in-place normalization retained for subclasses that override this
    ! binding. Plain prob assignment uses compact score vectors instead,
    ! preserving raw distances for ASSIGNMENT.dat and downstream score reporting.
    subroutine ref_normalize( self )
        class(eul_prob_tab), intent(inout) :: self
        real    :: sum_dist_all, min_dist, max_dist
        integer :: i, iref
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i,sum_dist_all)
        do i = 1, self%nptcls
            sum_dist_all = sum(self%loc_tab(:,i)%dist)
            if( sum_dist_all < TINY )then
                self%loc_tab(:,i)%dist = 0.
            else
                self%loc_tab(:,i)%dist = self%loc_tab(:,i)%dist / sum_dist_all
            endif
        enddo
        !$omp end parallel do
        !$omp parallel workshare proc_bind(close)
        min_dist = minval(self%loc_tab(:,:)%dist)
        max_dist = maxval(self%loc_tab(:,:)%dist)
        !$omp end parallel workshare
        if( (max_dist - min_dist) < TINY )then
            THROW_WARN('WARNING: numerical unstability in eul_prob_tab')
            !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iref,i)
            do iref = 1, self%nrefs
                do i = 1, self%nptcls
                    self%loc_tab(iref,i)%dist = ran3()
                enddo
            enddo
            !$omp end parallel do
        else
            self%loc_tab(:,:)%dist = (self%loc_tab(:,:)%dist - min_dist) / (max_dist - min_dist)
        endif
    end subroutine ref_normalize

    ! ptcl -> (proj, state) assignment using the global normalized dist value table
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
        logical :: l_likelihood
        character(len=256) :: alloc_msg
        l_likelihood = trim(self%p_ptr%prob_assign) == 'likelihood'
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
        if( l_likelihood )then
            call prepare_ref_score_vectors
        else
            call self%ref_normalize
        endif
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
                if( l_likelihood )then
                    call sample_likelihood_dist(self%nrefs, frontier_ref_dist, likelihood_nsample(self%nrefs, projs_athres),&
                        &dist_tmp, corr_tmp, assigned_iref, dists_sorted, inds_sorted)
                else
                    assigned_iref = angle_sampling(iref_dist, dists_sorted, inds_sorted, projs_athres, self%p_ptr%prob_athres)
                endif
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
                if( l_likelihood ) score = invalid_dist
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
            self%assgn_map(assigned_ptcl) = self%loc_tab(assigned_iref,assigned_ptcl)
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
                greedy_state(assigned_ptcl) = self%sinds(assigned_iref)
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
                if( self%sinds(iref) /= state_filter ) cycle
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
                if( self%sinds(iref) == state_filter ) call advance_ref_head(iref, state_filter)
            enddo
            do while( any(ptcl_avail) )
                call gather_active_refs(state_filter)
                if( nactive == 0 ) THROW_HARD('no active refs left in multi-state probability assignment')
                if( nactive == 1 )then
                    active_idx = 1
                else if( l_likelihood )then
                    call sample_likelihood_dist(nactive, active_ref_dist, likelihood_nsample(nactive, state_projs_athres(state_filter)),&
                        &dist_tmp, corr_tmp, active_idx, active_dists_sorted, active_inds)
                else
                    active_idx = angle_sampling(active_dists(1:nactive), active_dists_sorted(1:nactive),&
                        &active_inds(1:nactive), state_projs_athres(state_filter), self%p_ptr%prob_athres)
                endif
                assigned_iref = active_refs(active_idx)
                call assign_current_ref()
                do iref = 1,self%nrefs
                    if( self%sinds(iref) == state_filter ) call advance_ref_head(iref, state_filter)
                enddo
            enddo
        end subroutine assign_refs_for_state

        real function active_ref_dist( active_loc ) result(dist)
            integer, intent(in) :: active_loc
            dist = active_dists(active_loc)
        end function active_ref_dist

    end subroutine ref_assign

    ! state normalization (same energy) of the state_tab
    ! [0,1] normalization of the whole table
    subroutine state_normalize( self )
        class(eul_prob_tab), intent(inout) :: self
        real    :: sum_dist_all, min_dist, max_dist
        integer :: i, istate
        ! normalize so prob of each ptcl is between [0,1] for all states
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i,sum_dist_all)
        do i = 1, self%nptcls
            sum_dist_all = sum(self%state_tab(:,i)%dist)
            if( sum_dist_all < TINY )then
                self%state_tab(:,i)%dist = 0.
            else
                self%state_tab(:,i)%dist = self%state_tab(:,i)%dist / sum_dist_all
            endif
        enddo
        !$omp end parallel do
        ! min/max normalization to obtain values between 0 and 1
        max_dist = 0.
        min_dist = huge(min_dist)
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i)&
        !$omp reduction(min:min_dist) reduction(max:max_dist)
        do i = 1, self%nptcls
            max_dist = max(max_dist, maxval(self%state_tab(:,i)%dist, dim=1))
            min_dist = min(min_dist, minval(self%state_tab(:,i)%dist, dim=1))
        enddo
        !$omp end parallel do
        ! special case of numerical unstability of dist values
        if( (max_dist - min_dist) < TINY )then
            ! randomize dist so the assignment is stochastic
            !$omp parallel do default(shared) proc_bind(close) schedule(static) private(istate,i)
            do istate = 1, self%nstates
                do i = 1, self%nptcls
                    self%state_tab(istate,i)%dist = ran3()
                enddo
            enddo
            !$omp end parallel do
        else
            self%state_tab%dist = (self%state_tab%dist - min_dist) / (max_dist - min_dist)
        endif
    end subroutine state_normalize

    ! ptcl -> state (using assigned iproj or previous iproj) assignment
    subroutine state_assign( self )
        class(eul_prob_tab), intent(inout) :: self
        integer :: i, istate, assigned_istate, assigned_ptcl, state_dist_inds(self%nstates),&
                    &stab_inds(self%nptcls, self%nstates), inds_sorted(self%nstates)
        real    :: sorted_tab(self%nptcls, self%nstates), state_dist(self%nstates), state_dists_sorted(self%nstates)
        real    :: dist_tmp, corr_tmp
        logical :: ptcl_avail(self%nptcls)
        logical :: l_likelihood
        if( self%nstates == 1 )then
            self%assgn_map = self%state_tab(1,:)
            self%assgn_map%frac = 100.
            return
        endif
        l_likelihood = trim(self%p_ptr%prob_assign) == 'likelihood'
        if( .not. l_likelihood ) call self%state_normalize
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
            if( l_likelihood )then
                call sample_likelihood_dist(self%nstates, state_frontier_dist, self%nstates, dist_tmp, corr_tmp,&
                    &assigned_istate, state_dists_sorted, inds_sorted)
            else
                assigned_istate = minloc(state_dist, dim=1)
            endif
            assigned_ptcl   = stab_inds(state_dist_inds(assigned_istate), assigned_istate)
            ptcl_avail(assigned_ptcl)     = .false.
            self%assgn_map(assigned_ptcl) = self%state_tab(assigned_istate,assigned_ptcl)
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

    ! FILE IO

    ! write the partition-wise (or global) dist value table to a binary file
    subroutine write_tab( self, binfname )
        class(eul_prob_tab), intent(in) :: self
        class(string),       intent(in) :: binfname
        type(ptcl_ref) :: ptcl_ref_sample
        integer         :: funit, addr, io_stat, file_header(2)
        integer         :: first_ptcl, last_ptcl
        integer(kind=8) :: addr8, expected_bytes, file_bytes, elem_bytes
        file_header(1) = self%nrefs
        file_header(2) = self%nptcls
        call fopen(funit,binfname,access='STREAM',action='WRITE',status='REPLACE', iostat=io_stat)
        call fileiochk('simple_eul_prob_tab; write_tab; file: '//binfname%to_char(), io_stat)
        write(unit=funit,pos=1,iostat=io_stat) file_header
        call fileiochk('simple_eul_prob_tab; write_tab header; file: '//binfname%to_char(), io_stat)
        addr = sizeof(file_header) + 1
        call write_seed_shift_table(funit, addr, self%seed_nrots, self%seed_shifts, self%seed_has_sh)
        addr8 = int(addr,kind=8)
        elem_bytes = int(sizeof(ptcl_ref_sample), kind=8)
        do first_ptcl = 1, self%nptcls, PROB_TAB_IO_CHUNK
            last_ptcl = min(first_ptcl + PROB_TAB_IO_CHUNK - 1, self%nptcls)
            write(funit,pos=addr8,iostat=io_stat) self%loc_tab(:,first_ptcl:last_ptcl)
            call fileiochk('simple_eul_prob_tab; write_tab loc_tab; file: '//binfname%to_char(), io_stat)
            addr8 = addr8 + int(self%nrefs,kind=8) * int(last_ptcl-first_ptcl+1,kind=8) * elem_bytes
        enddo
        expected_bytes = addr8 - 1_8
        inquire(unit=funit, size=file_bytes, iostat=io_stat)
        if( io_stat == 0 .and. file_bytes < expected_bytes )then
            write(logfhandle,*) 'prob_tab write size check failed: ', binfname%to_char()
            write(logfhandle,*) 'actual/expected bytes: ', file_bytes, expected_bytes
            THROW_HARD('probability table write truncated: '//binfname%to_char())
        endif
        call fclose(funit)
    end subroutine write_tab

    ! read the partition-wise dist value binary file to global reg object's dist value table
    subroutine read_tab_to_glob( self, binfname )
        class(eul_prob_tab), intent(inout) :: self
        class(string),       intent(in)    :: binfname
        type(ptcl_ref),      allocatable   :: mat_chunk(:,:)
        type(ptcl_ref)                     :: ptcl_ref_sample
        real,                allocatable   :: seed_shifts_loc(:,:)
        logical,             allocatable   :: seed_has_sh_loc(:)
        integer, allocatable :: pind2glob(:)
        integer :: funit, addr, io_stat, file_header(2), nptcls_loc, nrefs_loc, i_loc, i_glob, pind, max_pind, seed_nrots_loc
        integer :: first_loc, last_loc, nchunk, i
        integer(kind=8) :: addr8, file_bytes, expected_bytes, tab_bytes, elem_bytes
        if( file_exists(binfname) )then
            call fopen(funit,binfname,access='STREAM',action='READ',status='OLD', iostat=io_stat)
            call fileiochk('simple_eul_prob_tab; read_tab_to_glob; file: '//binfname%to_char(), io_stat)
        else
            THROW_HARD( 'corr/rot files of partitions should be ready! ' )
        endif
        ! reading header and the nprojs/nptcls in this partition file
        read(unit=funit,pos=1) file_header
        nrefs_loc  = file_header(1)
        nptcls_loc = file_header(2)
        if( nrefs_loc .ne. self%nrefs ) THROW_HARD( 'nrefs should be the same as nrefs in this partition file!' )
        allocate(seed_shifts_loc(2,nptcls_loc), seed_has_sh_loc(nptcls_loc))
        ! read partition information
        addr = sizeof(file_header) + 1
        call read_seed_shift_table(funit, addr, seed_nrots_loc, seed_shifts_loc, seed_has_sh_loc)
        addr8 = int(addr,kind=8)
        elem_bytes     = int(sizeof(ptcl_ref_sample),kind=8)
        tab_bytes      = int(nrefs_loc,kind=8) * int(nptcls_loc,kind=8) * elem_bytes
        expected_bytes = addr8 + tab_bytes - 1_8
        file_bytes = -1_8
        inquire(file=binfname%to_char(), size=file_bytes, iostat=io_stat)
        if( io_stat /= 0 .or. file_bytes < expected_bytes )then
            write(logfhandle,*) 'prob_tab read size check failed: ', binfname%to_char()
            write(logfhandle,*) 'nrefs/nptcls/actual/expected bytes: ', nrefs_loc, nptcls_loc, file_bytes, expected_bytes
            THROW_HARD('probability table is missing or truncated: '//binfname%to_char())
        endif
        if( self%seed_nrots == 0 ) self%seed_nrots = seed_nrots_loc
        if( self%seed_nrots /= seed_nrots_loc ) THROW_HARD('seed_nrots mismatch in eul_prob_tab%read_tab_to_glob')
        max_pind = maxval(self%pinds)
        if( max_pind < 1 )then
            call fclose(funit)
            deallocate(seed_shifts_loc, seed_has_sh_loc)
            return
        endif
        allocate(pind2glob(max_pind), source=0)
        do i = 1, size(self%pinds)
            pind = self%pinds(i)
            if( pind > 0 .and. pind <= max_pind ) pind2glob(pind) = i
        enddo
        allocate(mat_chunk(nrefs_loc, min(PROB_TAB_IO_CHUNK, nptcls_loc)))
        do first_loc = 1, nptcls_loc, PROB_TAB_IO_CHUNK
            last_loc = min(first_loc + PROB_TAB_IO_CHUNK - 1, nptcls_loc)
            nchunk   = last_loc - first_loc + 1
            read(unit=funit,pos=addr8,iostat=io_stat) mat_chunk(:,1:nchunk)
            call fileiochk('simple_eul_prob_tab; read_tab_to_glob loc_tab; file: '//binfname%to_char(), io_stat)
            addr8 = addr8 + int(nrefs_loc,kind=8) * int(nchunk,kind=8) * elem_bytes
            !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i_loc,i_glob,pind)
            do i_loc = 1, nchunk
                pind = mat_chunk(1,i_loc)%pind
                if( pind < 1 .or. pind > max_pind ) cycle
                i_glob = pind2glob(pind)
                if( i_glob > 0 )then
                    self%loc_tab(:,i_glob)    = mat_chunk(:,i_loc)
                    self%seed_shifts(:,i_glob)= seed_shifts_loc(:,first_loc+i_loc-1)
                    self%seed_has_sh(i_glob)  = seed_has_sh_loc(first_loc+i_loc-1)
                endif
            end do
            !$omp end parallel do
        enddo
        call fclose(funit)
        deallocate(mat_chunk, seed_shifts_loc, seed_has_sh_loc, pind2glob)
    end subroutine read_tab_to_glob

    subroutine write_state_tab( self, binfname )
        class(eul_prob_tab), intent(in) :: self
        class(string),       intent(in) :: binfname
        integer :: funit, io_stat, headsz
        headsz = sizeof(self%nptcls)
        call fopen(funit,binfname,access='STREAM',action='WRITE',status='REPLACE', iostat=io_stat)
        write(unit=funit,pos=1)          self%nptcls
        write(unit=funit,pos=headsz + 1) self%state_tab
        call fclose(funit)
    end subroutine write_state_tab

    subroutine read_state_tab( self, binfname )
        class(eul_prob_tab), intent(inout) :: self
        class(string),       intent(in)    :: binfname
        type(ptcl_ref),      allocatable   :: state_tab_glob(:,:)
        integer, allocatable :: pind2glob(:), pinds_glob(:)
        integer :: funit, io_stat, nptcls_glob, headsz, i_loc, i_glob, pind, max_pind
        headsz = sizeof(nptcls_glob)
        if( .not. file_exists(binfname) )then
            THROW_HARD('file '//binfname%to_char()//' does not exists!')
        else
            call fopen(funit,binfname,access='STREAM',action='READ',status='OLD', iostat=io_stat)
        end if
        call fileiochk('simple_eul_prob_tab; read_state_tab; file: '//binfname%to_char(), io_stat)
        read(unit=funit,pos=1) nptcls_glob
        allocate(state_tab_glob(self%nstates,nptcls_glob))
        read(unit=funit,pos=headsz + 1) state_tab_glob
        call fclose(funit)
        pinds_glob = state_tab_glob(1,:)%pind
        call build_pind_lookup(pinds_glob, self%pinds, pind2glob, max_pind)
        if( max_pind < 1 )then
            deallocate(state_tab_glob, pind2glob, pinds_glob)
            return
        endif
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i_loc,i_glob,pind)
        do i_loc = 1, self%nptcls
            pind = self%pinds(i_loc)
            if( pind < 1 .or. pind > max_pind ) cycle
            i_glob = pind2glob(pind)
            if( i_glob > 0 ) self%state_tab(:,i_loc) = state_tab_glob(:,i_glob)
        end do
        !$omp end parallel do
        deallocate(state_tab_glob, pind2glob, pinds_glob)
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
        if( allocated(self%loc_tab)      ) deallocate(self%loc_tab)
        if( allocated(self%state_tab)    ) deallocate(self%state_tab)
        if( allocated(self%assgn_map)    ) deallocate(self%assgn_map)
        if( allocated(self%seed_shifts)  ) deallocate(self%seed_shifts)
        if( allocated(self%seed_has_sh)  ) deallocate(self%seed_has_sh)
        if( allocated(self%pinds)        ) deallocate(self%pinds)
        if( allocated(self%ssinds)       ) deallocate(self%ssinds)
        if( allocated(self%sinds)        ) deallocate(self%sinds)
        if( allocated(self%jinds)        ) deallocate(self%jinds)
        if( allocated(self%state_exists) ) deallocate(self%state_exists)
        if( allocated(self%proj_exists)  ) deallocate(self%proj_exists)
        self%seed_nrots = 0
        self%b_ptr => null()
        self%p_ptr => null()
    end subroutine kill

end module simple_eul_prob_tab
