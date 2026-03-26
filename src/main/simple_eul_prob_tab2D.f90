!@descr: 2D probability table routines for multi-reference class assignment with probabilistic sampling
module simple_eul_prob_tab2D
use simple_pftc_srch_api
use simple_builder,          only: builder
use simple_pftc_shsrch_grad, only: pftc_shsrch_grad
use simple_decay_funs,       only: extremal_decay2D
use simple_eul_prob_tab,     only: eulprob_dist_switch, eulprob_corr_switch
implicit none

public :: eul_prob_tab2D
private
#include "simple_local_flags.inc"

real, parameter :: NHOOD_FRAC = 0.1

type :: eul_prob_tab2D
    class(builder),    pointer  :: b_ptr => null()
    class(parameters), pointer  :: p_ptr => null()
    type(ptcl_ref), allocatable :: loc_tab(:,:)      !< 2D search table (nclasses, nptcls)
    type(ptcl_ref), allocatable :: assgn_map(:)      !< assignment map (nptcls)
    integer,        allocatable :: pinds(:)          !< particle indices for processing
    logical,        allocatable :: class_exists(:)   !< class population filter
    integer                     :: nptcls            !< size of pinds array
    integer                     :: nclasses          !< number of classes
    integer                     :: nhood_sz = 1      !< probabilistic neighborhood size used in fill_tab/ref_assign
    contains
    ! CONSTRUCTOR
    procedure :: new
    ! MAIN PROCEDURES
    procedure :: fill_tab
    procedure :: ref_assign => ref_assign_greedy
    procedure :: write_tab
    procedure :: read_tab_to_glob
    procedure :: write_assignment
    procedure :: read_assignment
    ! DESTRUCTOR
    procedure :: kill
    ! PRIVATE
    procedure, private :: ref_normalize
end type eul_prob_tab2D

contains

    ! CONSTRUCTORS

    subroutine new( self, params, build, pinds, empty_okay )
        class(eul_prob_tab2D),     intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        class(builder),    target, intent(in)    :: build
        integer,                   intent(in)    :: pinds(:)
        logical, optional,         intent(in)    :: empty_okay
        integer :: i, icls, iptcl, nactive
        real    :: x
        logical :: l_empty
        call self%kill
        l_empty = .false.
        if( present(empty_okay) ) l_empty = empty_okay
        self%p_ptr     => params
        self%b_ptr     => build
        self%nptcls    = size(pinds)
        self%nclasses  = params%ncls
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
            do icls = 1, self%nclasses
                self%loc_tab(icls,i)%pind   = iptcl
                self%loc_tab(icls,i)%icls   = icls
                self%loc_tab(icls,i)%inpl   = 0
                self%loc_tab(icls,i)%dist   = huge(x)
                self%loc_tab(icls,i)%x      = 0.
                self%loc_tab(icls,i)%y      = 0.
                self%loc_tab(icls,i)%has_sh = .false.
            end do
        end do
        !$omp end parallel do
    end subroutine new

    ! table filling for 2D multi-class assignment
    subroutine fill_tab( self )
        class(eul_prob_tab2D), intent(inout) :: self
        type(pftc_shsrch_grad) :: grad_shsrch_obj(nthr_glob)  !< shift search object, L-BFGS with gradient
        type(ori)              :: o_prev
        integer :: i, icls, iptcl, ithr, irot, irot0, icls_prev, iref_n, nactive, nhood_sz_loc, ninpl_smpl, order_ind
        integer :: active_cls(self%nclasses)
        integer :: loc(1), vec_nrots(self%b_ptr%pftc%get_nrots())
        real    :: lims(2,2), lims_init(2,2), cxy(3), cxy_prob(3), rotmat(2,2), rot_xy(2)
        real    :: inpl_corrs(self%b_ptr%pftc%get_nrots()), cls_dists(self%nclasses), cls_dists_work(self%nclasses)
        real    :: inpl_corr, power, neigh_frac
        logical :: class_active(self%nclasses)
        call seed_rnd
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
        power      = EXTR_POWER
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
            !$omp parallel do default(shared) private(i,iptcl,ithr,o_prev,irot,irot0,cxy,icls,icls_prev,inpl_corrs,loc,cls_dists,cls_dists_work,iref_n,cxy_prob,vec_nrots,order_ind,inpl_corr,rotmat,rot_xy) proc_bind(close) schedule(static)
            do i = 1, self%nptcls
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
                ! (2) search class references using that shared shift
                cls_dists = huge(1.0)
                do iref_n = 1, nactive
                    icls = active_cls(iref_n)
                    call self%b_ptr%pftc%gen_objfun_vals(icls, iptcl, cxy(2:3), inpl_corrs)
                    call power_sampling(power, self%b_ptr%pftc%get_nrots(), inpl_corrs, vec_nrots, ninpl_smpl, irot, order_ind, inpl_corr)
                    rot_xy = 0.
                    if( irot > 0 )then
                        call rotmat2d(self%b_ptr%pftc%get_rot(irot), rotmat)
                        rot_xy = matmul(cxy(2:3), rotmat)
                    endif
                    self%loc_tab(icls,i)%dist   = eulprob_dist_switch(inpl_corr, self%p_ptr%cc_objfun)
                    self%loc_tab(icls,i)%inpl   = irot
                    self%loc_tab(icls,i)%x      = rot_xy(1)
                    self%loc_tab(icls,i)%y      = rot_xy(2)
                    self%loc_tab(icls,i)%has_sh = irot > 0
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
            !$omp parallel do default(shared) private(i,iptcl,ithr,icls,iref_n,inpl_corrs,vec_nrots,irot,order_ind,inpl_corr) proc_bind(close) schedule(static)
            do i = 1, self%nptcls
                iptcl = self%pinds(i)
                ithr  = omp_get_thread_num() + 1
                do iref_n = 1, nactive
                    icls = active_cls(iref_n)
                    call self%b_ptr%pftc%gen_objfun_vals(icls, iptcl, [0.,0.], inpl_corrs)
                    call power_sampling(power, self%b_ptr%pftc%get_nrots(), inpl_corrs, vec_nrots, ninpl_smpl, irot, order_ind, inpl_corr)
                    self%loc_tab(icls,i)%dist   = eulprob_dist_switch(inpl_corr, self%p_ptr%cc_objfun)
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
    end subroutine fill_tab

    ! class normalization (same energy) of the loc_tab, [0,1] normalization
    subroutine ref_normalize( self )
        class(eul_prob_tab2D), intent(inout) :: self
        real    :: sum_dist_all, min_dist, max_dist
        integer :: i, icls, nactive
        logical :: class_active(self%nclasses)
        class_active = self%class_exists
        nactive      = count(class_active)
        if( nactive == 0 )then
            ! fallback to avoid table degeneration when all classes are temporarily low-population
            class_active = .true.
            nactive      = self%nclasses
            THROW_WARN('No active classes after population filtering; falling back to all classes in eul_prob_tab2D')
        endif
        ! normalize so prob of each ptcl is between [0,1] for all classes
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i,icls,sum_dist_all)
        do i = 1, self%nptcls
            sum_dist_all = 0.
            do icls = 1, self%nclasses
                if( class_active(icls) )then
                    sum_dist_all = sum_dist_all + self%loc_tab(icls,i)%dist
                else
                    self%loc_tab(icls,i)%dist = 0.
                endif
            end do
            if( sum_dist_all < TINY )then
                do icls = 1, self%nclasses
                    if( class_active(icls) ) self%loc_tab(icls,i)%dist = 0.
                end do
            else
                do icls = 1, self%nclasses
                    if( class_active(icls) ) self%loc_tab(icls,i)%dist = self%loc_tab(icls,i)%dist / sum_dist_all
                end do
            endif
        end do
        !$omp end parallel do
        ! min/max normalization to obtain values between 0 and 1
        min_dist = huge(min_dist)
        max_dist = -huge(max_dist)
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i,icls) reduction(min:min_dist) reduction(max:max_dist)
        do i = 1, self%nptcls
            do icls = 1, self%nclasses
                if( .not. class_active(icls) ) cycle
                min_dist = min(min_dist, self%loc_tab(icls,i)%dist)
                max_dist = max(max_dist, self%loc_tab(icls,i)%dist)
            end do
        end do
        !$omp end parallel do
        ! special case of numerical unstability of dist values
        if( (max_dist - min_dist) < TINY )then
            THROW_WARN('WARNING: numerical unstability in eul_prob_tab2D')
            ! randomize dist so the assignment is stochastic
            !$omp parallel do default(shared) proc_bind(close) schedule(static) private(icls,i)
            do icls = 1, self%nclasses
                if( .not. class_active(icls) )then
                    self%loc_tab(icls,:)%dist = 0.
                    cycle
                endif
                do i = 1, self%nptcls
                    self%loc_tab(icls,i)%dist = ran3()
                end do
            end do
            !$omp end parallel do
        else
            !$omp parallel do default(shared) proc_bind(close) schedule(static) private(icls,i)
            do icls = 1, self%nclasses
                if( .not. class_active(icls) )then
                    self%loc_tab(icls,:)%dist = 0.
                    cycle
                endif
                do i = 1, self%nptcls
                    self%loc_tab(icls,i)%dist = (self%loc_tab(icls,i)%dist - min_dist) / (max_dist - min_dist)
                end do
            end do
            !$omp end parallel do
        endif
    end subroutine ref_normalize

    ! deterministic ptcl -> class assignment (debug): no probabilistic sampling
    subroutine ref_assign_greedy( self )
        class(eul_prob_tab2D), intent(inout) :: self
        integer :: i, icls, assigned_icls, assigned_ptcl, stab_inds(self%nptcls, self%nclasses)
        integer :: icls_dist_inds(self%nclasses), active_cls(self%nclasses), eligible_cls(self%nclasses)
        real    :: sorted_tab(self%nptcls, self%nclasses), dists_sorted(self%nclasses)
        real    :: dists_raw(self%nclasses, self%nptcls)
        logical :: ptcl_avail(self%nptcls)
        logical :: class_active(self%nclasses)
        integer :: nactive, iact, chosen_active, neligible
        real    :: best_dist
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
        ! keep raw distances for output (objective values), normalize only for deterministic assignment ordering
        dists_raw = self%loc_tab(:,:)%dist
        ! normalization
        call self%ref_normalize
        ! sorting each columns
        sorted_tab = transpose(self%loc_tab(:,:)%dist)
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
        do while( any(ptcl_avail) )
            ! collect current best distance for each class, skipping exhausted ones
            neligible = 0
            do iact = 1, nactive
                icls = active_cls(iact)
                do while( icls_dist_inds(icls) <= self%nptcls )
                    assigned_ptcl = stab_inds(icls_dist_inds(icls), icls)
                    if( ptcl_avail(assigned_ptcl) ) exit
                    icls_dist_inds(icls) = icls_dist_inds(icls) + 1
                end do
                if( icls_dist_inds(icls) <= self%nptcls )then
                    neligible = neligible + 1
                    eligible_cls(neligible) = icls
                    dists_sorted(neligible) = sorted_tab(icls_dist_inds(icls), icls)
                endif
            end do
            if( neligible == 0 ) exit
            chosen_active = minloc(dists_sorted(1:neligible), dim=1)
            assigned_icls = eligible_cls(chosen_active)
            assigned_ptcl = stab_inds(icls_dist_inds(assigned_icls), assigned_icls)
            ptcl_avail(assigned_ptcl)     = .false.
            self%assgn_map(assigned_ptcl) = self%loc_tab(assigned_icls, assigned_ptcl)
            self%assgn_map(assigned_ptcl)%dist = dists_raw(assigned_icls, assigned_ptcl)
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
                if( assigned_icls == 0 ) THROW_HARD('Failed particle assignment in eul_prob_tab2D%ref_assign_greedy')
                self%assgn_map(i) = self%loc_tab(assigned_icls, i)
                self%assgn_map(i)%dist = dists_raw(assigned_icls, i)
                ptcl_avail(i)     = .false.
            end do
        endif
    end subroutine ref_assign_greedy

    ! DESTRUCTOR

    subroutine kill( self )
        class(eul_prob_tab2D), intent(inout) :: self
        if( allocated(self%loc_tab)        ) deallocate(self%loc_tab)
        if( allocated(self%assgn_map)      ) deallocate(self%assgn_map)
        if( allocated(self%pinds)          ) deallocate(self%pinds)
        if( allocated(self%class_exists)   ) deallocate(self%class_exists)
        self%b_ptr => null()
        self%p_ptr => null()
    end subroutine kill

    ! FILE IO

    subroutine write_tab( self, binfname )
        class(eul_prob_tab2D), intent(in) :: self
        class(string),         intent(in) :: binfname
        integer :: funit, addr, io_stat, file_header(2)
        file_header(1) = self%nclasses
        file_header(2) = self%nptcls
        call fopen(funit, binfname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        write(unit=funit, pos=1) file_header
        addr = sizeof(file_header) + 1
        write(funit, pos=addr) self%loc_tab
        call fclose(funit)
    end subroutine write_tab

    subroutine read_tab_to_glob( self, binfname )
        class(eul_prob_tab2D), intent(inout) :: self
        class(string),         intent(in)    :: binfname
        type(ptcl_ref), allocatable :: mat_loc(:,:)
        integer :: funit, addr, io_stat, file_header(2), nptcls_loc, nclasses_loc, i_loc, i_glob
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
        addr = sizeof(file_header) + 1
        read(unit=funit, pos=addr) mat_loc
        call fclose(funit)
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i_loc,i_glob)
        do i_glob = 1, self%nptcls
            do i_loc = 1, nptcls_loc
                if( mat_loc(1,i_loc)%pind == self%loc_tab(1,i_glob)%pind )then
                    self%loc_tab(:,i_glob) = mat_loc(:,i_loc)
                    exit
                endif
            end do
        end do
        !$omp end parallel do
        deallocate(mat_loc)
    end subroutine read_tab_to_glob

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
        integer :: funit, io_stat, nptcls_glob, headsz, i_loc, i_glob
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
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i_loc,i_glob)
        do i_loc = 1, self%nptcls
            do i_glob = 1, nptcls_glob
                if( self%assgn_map(i_loc)%pind == assgn_glob(i_glob)%pind )then
                    self%assgn_map(i_loc) = assgn_glob(i_glob)
                    exit
                endif
            end do
        end do
        !$omp end parallel do
        deallocate(assgn_glob)
    end subroutine read_assignment

end module simple_eul_prob_tab2D
