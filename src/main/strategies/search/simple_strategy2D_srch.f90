!@descr: common strategy2D methods and type specification for polymorphic strategy2D object creation are delegated to this class
module simple_strategy2D_srch
use simple_core_module_api, only: dp
use simple_pftc_srch_api
use simple_strategy2D_alloc, only: prep_strategy2D_thread, s2D, is_fresh_2D_start
use simple_eul_prob_tab2D,   only: eul_prob_tab2D
use simple_pftc_shsrch_grad, only: pftc_shsrch_grad
use simple_builder,          only: builder
implicit none

public :: strategy2D_srch, strategy2D_spec
private
#include "simple_local_flags.inc"

type strategy2D_spec
    type(eul_prob_tab2D), pointer :: eulprob_obj_part2D => null()
    real    :: stoch_bound = 0.
    integer :: iptcl       = 0  ! global particle index
    integer :: iptcl_batch = 0  ! maps to index in batch
    integer :: iptcl_map   = 0  ! index in all contiguous batches
end type strategy2D_spec

type strategy2D_srch
    class(parameters), pointer :: p_ptr => null()     ! pointer to parameters
    class(builder),    pointer :: b_ptr => null()     !< build handle for access to pftc
    type(pftc_shsrch_grad) :: grad_shsrch_obj         !< origin shift search object, L-BFGS with gradient
    type(pftc_shsrch_grad) :: grad_shsrch_obj2        !< origin shift search object, L-BFGS with gradient, no call back
    type(pftc_shsrch_grad) :: grad_shsrch_first_obj   !< origin shift search object, L-BFGS with gradient, used for initial shift search on previous ref
    type(pftc_shsrch_grad) :: grad_shsrch_joint_obj   !< classical Euclidean joint (sx,sy,theta) L-BFGS-B refinement
    integer                 :: nrefs           =  0   !< number of references
    integer                 :: nrots           =  0   !< number of in-plane rotations in polar representation
    integer                 :: nrefs_eval      =  0   !< nr of references evaluated
    integer                 :: nsolns          =  0   !< # distinct refs with a stored solution this search
    integer                 :: prev_class      =  0   !< previous class index
    integer                 :: prev_class_mi   =  0   !< previous class used for overlap reporting
    integer                 :: best_class      =  0   !< best class index found by search
    integer                 :: best_rot        =  0   !< best in-plane rotation found by search
    integer                 :: prev_rot        =  0   !< previous in-plane rotation found by search
    integer                 :: iptcl           =  0   !< global particle index
    integer                 :: iptcl_batch     =  0   !< index in pre-allocated batch array
    integer                 :: iptcl_map       =  0   !< index in all contiguous batches combined
    integer                 :: ithr            =  0   !< current thread
    real                    :: prev_shvec(2)   =  0.  !< previous origin shift vector
    real                    :: best_shvec(2)   =  0.  !< best shift vector found by search
    real                    :: xy_first(2)     =  0.  !< initial shifts identified by searching the previous best reference
    real                    :: xy_first_rot(2) =  0.  !< initial shifts identified by searching the previous best reference, rotated
    real                    :: prev_corr       = -1.  !< previous best correlation
    real                    :: best_corr       = -1.  !< best corr found by search
    real(dp)                :: sgd_objective_initial = 0._dp  !< raw Euclidean objective at stream entry
    real(dp)                :: sgd_objective_final   = 0._dp  !< raw Euclidean objective after stream steps
    real                    :: trs             =  0.  !< shift boundary
    integer                 :: sgd_accepted_steps = 0  !< accepted bounded gradient steps
    logical                 :: sgd_used          = .false. !< stream shift update was attempted
    logical                 :: l_sh_first      = .false. !< Whether to search the shifts on previous best reference
    logical                 :: l_fresh_start   = .false. !< Whether previous alignment parameters are intentionally ignored
    logical                 :: use_joint_angle = .false. !< selected Euclidean candidates use joint angle refinement
    logical                 :: best_inpl_e3_valid = .false. !< whether best_inpl_e3 is continuous
    real                    :: best_inpl_e3    = 0. !< continuous e3 angle for the selected candidate
  contains
    procedure :: new
    procedure :: prep4srch
    procedure :: inpl_srch_first
    procedure :: inpl_srch
    procedure :: inpl_srch_peaks
    procedure :: store_solution
    procedure :: assign_ori
    procedure :: kill
end type strategy2D_srch

contains

    subroutine new( self, params, spec, build )
        class(strategy2D_srch),    intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        class(strategy2D_spec),    intent(in)    :: spec
        class(builder),    target, intent(in)    :: build
        real :: lims(2,2), lims_init(2,2), joint_lims(3,2)
        call self%kill
        ! set pointer to parameters
        self%p_ptr => params
        ! set pointer to builder
        self%b_ptr => build
        ! set constants
        self%iptcl       = spec%iptcl
        self%iptcl_batch = spec%iptcl_batch
        self%iptcl_map   = spec%iptcl_map
        self%nrefs       = self%p_ptr%ncls
        self%nrots       = self%b_ptr%pftc%get_nrots()
        if( self%nrots < 1 ) THROW_HARD('strategy2D_srch constructed before PFTC rotations were initialized')
        self%nrefs_eval  = 0
        self%use_joint_angle = trim(self%p_ptr%inpl_refine) == 'yes' .and. &
            &self%b_ptr%pftc%is_euclid_objfun() .and. &
            &(.not. self%p_ptr%l_sgd_streaming_active) .and. trim(self%p_ptr%tseries) /= 'yes' .and. &
            &(.not. self%p_ptr%l_prob_align_mode)
        self%best_inpl_e3_valid = .false.
        ! construct composites
        self%trs        =  self%p_ptr%trs
        lims(:,1)       = -self%p_ptr%trs
        lims(:,2)       =  self%p_ptr%trs
        lims_init(:,1)  = -SHC_INPL_TRSHWDTH
        lims_init(:,2)  =  SHC_INPL_TRSHWDTH
        if( self%p_ptr%l_sgd_streaming_active )then
            ! The stream has already selected one discrete in-plane rotation.
            ! Construct shift-only objects so minimize_direct differentiates
            ! L(c,r,s) with respect to s=(sx,sy) at fixed (c,r), and do not
            ! allocate the legacy particle-shift L-BFGS-B optimizer.
            call self%grad_shsrch_obj%new(self%b_ptr, lims, lims_init=lims_init,&
            maxits=self%p_ptr%maxits_sh, opt_angle=.false., direct_only=.true.)
            call self%grad_shsrch_first_obj%new(self%b_ptr, lims, lims_init=lims_init,&
            maxits=self%p_ptr%maxits_sh, opt_angle=.false., coarse_init=.true., direct_only=.true.)
            call self%grad_shsrch_obj%set_diagnostic_mode(self%p_ptr%sgd_diagnostic)
            call self%grad_shsrch_first_obj%set_diagnostic_mode(self%p_ptr%sgd_diagnostic)
            if( .not. self%grad_shsrch_obj%is_direct_shift_only() )then
                THROW_HARD('stream shift optimizer is not configured for fixed-angle direct descent')
            endif
            if( .not. self%grad_shsrch_first_obj%is_direct_shift_only() )then
                THROW_HARD('stream seed-shift optimizer is not configured for fixed-angle direct descent')
            endif
        else if( trim(self%p_ptr%tseries).eq.'yes' )then
            ! shift only search
            call self%grad_shsrch_obj%new(self%b_ptr, lims, lims_init=lims_init,&
            maxits=self%p_ptr%maxits_sh, opt_angle=.false.)
            call self%grad_shsrch_first_obj%new(self%b_ptr, lims, lims_init=lims_init,&
            maxits=self%p_ptr%maxits_sh, opt_angle=.false., coarse_init=.true.)
        else
            call self%grad_shsrch_obj%new(self%b_ptr, lims, lims_init=lims_init,&
            maxits=self%p_ptr%maxits_sh)
            call self%grad_shsrch_first_obj%new(self%b_ptr, lims, lims_init=lims_init,&
            maxits=self%p_ptr%maxits_sh, coarse_init=.true.)
        endif
        call self%grad_shsrch_obj2%new(self%b_ptr, lims, lims_init=lims_init,&
        &maxits=self%p_ptr%maxits_sh, opt_angle=.false.)
        if( self%use_joint_angle )then
            joint_lims(1:2,:) = lims
            joint_lims(3,:) = [1.-2., real(self%nrots)+2.]
            call self%grad_shsrch_joint_obj%new(self%b_ptr, joint_lims, &
                &shbarrier=self%p_ptr%shbarrier, maxits=self%p_ptr%maxits_sh, &
                &opt_angle=.false., joint_angle=.true.)
        endif
    end subroutine new

    subroutine prep4srch( self, os )
        class(strategy2D_srch), intent(inout) :: self
        class(oris),            intent(inout) :: os
        real    :: corrs(self%nrots)
        logical :: has_been_searched
        self%nrefs_eval = 0
        self%nsolns     = 0
        self%sgd_objective_initial = 0.
        self%sgd_objective_final   = 0.
        self%sgd_accepted_steps    = 0
        self%sgd_used              = .false.
        self%best_inpl_e3_valid    = .false.
        self%ithr       = omp_get_thread_num() + 1
        ! find previous discrete alignment parameters
        self%l_fresh_start = is_fresh_2D_start(self%p_ptr, self%p_ptr%which_iter)
        if( self%l_fresh_start )then
            self%prev_class    = 0
            self%prev_class_mi = 0
            self%prev_rot      = 1
            self%prev_shvec    = 0.
        else
            self%prev_class_mi = nint(os%get(self%iptcl,'class'))            ! class index before any fallback
            self%prev_class    = self%prev_class_mi
            self%prev_rot      = self%b_ptr%pftc%get_roind(360.-os%e3get(self%iptcl)) ! in-plane angle index
        if( self%prev_rot < 1 .or. self%prev_rot > self%nrots )then
            if( self%p_ptr%sgd_diagnostic )then
                write(logfhandle,'(A,1X,I0,1X,A,I0)') &
                    '>>> SEARCH SAFETY: invalid previous rotation=', self%prev_rot, 'nrots=', self%nrots
            endif
            THROW_HARD('Invalid previous in-plane rotation index')
        endif
            self%prev_shvec = os%get_2Dshift(self%iptcl)                  ! shift vector
        endif
        self%best_shvec = 0.
        if( self%prev_class > 0 )then
            if( s2D%cls_pops(self%prev_class) > 0 )then
                ! all done
            else
                ! for limiting cases
                self%prev_class = irnd_uni(self%nrefs)
                do while( s2D%cls_pops(self%prev_class) <= 0 )
                    self%prev_class = irnd_uni(self%nrefs)
                enddo
            endif
        else
            ! initialization
            self%prev_class = irnd_uni(self%nrefs)
            do while( s2D%cls_pops(self%prev_class) <= 0 )
                self%prev_class = irnd_uni(self%nrefs)
            enddo
        endif
        ! set best to previous best by default
        self%best_class = self%prev_class
        self%best_rot   = self%prev_rot
        ! calculate previous best corr (treshold for better)
        if( self%p_ptr%l_sgd_streaming_active )then
            call self%b_ptr%pftc%gen_raw_euclid_vals(self%prev_class, self%iptcl, [0.,0.], corrs)
            ! Search bookkeeping is historically score-like (larger is better).
            ! Store -L internally while the stream selects argmin L.
            self%prev_corr = -corrs(self%prev_rot)
        else
            call self%b_ptr%pftc%gen_objfun_vals(self%prev_class, self%iptcl, [0.,0.], corrs)
        endif
        if( self%p_ptr%l_sgd_streaming_active )then
            ! already assigned from the finite raw loss above
        else if( self%p_ptr%cc_objfun == OBJFUN_CC )then
            self%prev_corr  = max(0., corrs(self%prev_rot))
        else
            self%prev_corr  = corrs(self%prev_rot)
        endif
        self%best_corr = self%prev_corr
        ! whether to search shifts first
        has_been_searched = (.not. self%l_fresh_start) .and. os%has_been_searched(self%iptcl)
        self%l_sh_first   = s2D%do_inplsrch(self%iptcl_batch) .and. self%p_ptr%l_doshift &
            &.and. has_been_searched
        self%xy_first     =  0.
        self%xy_first_rot =  0.
        ! init per-thread class-space arrays
        call prep_strategy2D_thread(self%ithr)
    end subroutine prep4srch

    subroutine inpl_srch_first( self )
        class(strategy2D_srch), intent(inout) :: self
        real    :: cxy(3), rotmat(2,2)
        integer :: irot
        self%best_shvec = [0.,0.]
        if( .not. self%l_sh_first ) return
        ! Stream mode uses the previous discrete state only as a shift seed;
        ! class and angle remain a discrete search and the two shifts are
        ! refined by bounded analytical-gradient steps.
        irot = 0
        call self%grad_shsrch_first_obj%set_indices(self%prev_class, self%iptcl)
        if( self%p_ptr%l_sgd_streaming_active )then
            irot = self%prev_rot
            cxy = self%grad_shsrch_first_obj%minimize_direct(irot=irot, xy_in=[0.,0.],&
                &step_size=self%p_ptr%sgd_eta_shift, max_steps=self%p_ptr%sgd_shift_its,&
                &sh_rot=.false., raw_euclid=.true.)
        else
            if( .not.self%grad_shsrch_first_obj%does_opt_angle() )then
                ! shift-only optimization
                irot = self%prev_rot
            endif
            cxy = self%grad_shsrch_first_obj%minimize(irot=irot, sh_rot=.false.)
        endif
        if( irot == 0 ) cxy(2:3) = 0.
        self%xy_first = cxy(2:3)
        self%xy_first_rot = 0.
        if( irot > 0 )then
            ! rotate the shift vector to the frame of reference
            call rotmat2d(self%b_ptr%pftc%get_rot(irot), rotmat)
            self%xy_first_rot = matmul(cxy(2:3), rotmat)
            ! update best
            if( self%p_ptr%l_sgd_streaming_active )then
                ! The direct raw-loss minimizer returns the merit -L, matching
                ! the larger-is-better search bookkeeping used below.
                self%best_corr = real(cxy(1))
            else
                self%best_corr = cxy(1)
            endif
            self%best_rot   = irot
            self%best_shvec = self%xy_first_rot
        endif
    end subroutine inpl_srch_first

    subroutine inpl_srch( self )
        class(strategy2D_srch), intent(inout) :: self
        real    :: cxy(3), cxy_joint(3), theta_stage1, rotmat(2,2)
        real    :: joint_lims(3,2)
        real(dp) :: theta_joint
        integer :: irot, irot_stage1, irot_joint, iang
        irot = 0
        self%best_shvec = [0.,0.]
        if( s2D%do_inplsrch(self%iptcl_batch) )then
            ! Stream mode replaces particle-shift L-BFGS-B after the discrete
            ! class/angle winner.  The direct minimizer retains the input state
            ! when no tested bounded trial improves the loss.
            call self%grad_shsrch_obj%set_indices(self%best_class, self%iptcl)
            if( self%p_ptr%l_sgd_streaming_active )then
                irot = self%best_rot
                self%sgd_used = .true.
                if( self%l_sh_first )then
                    ! Keep the no-improvement state in the particle frame,
                    ! matching the legacy minimizer's handoff convention.
                    self%best_shvec = self%xy_first_rot
                    cxy = self%grad_shsrch_obj%minimize_direct(irot=irot, xy_in=self%xy_first,&
                        &step_size=self%p_ptr%sgd_eta_shift, max_steps=self%p_ptr%sgd_shift_its,&
                        &sh_rot=.true., raw_euclid=.true., &
                        &accepted_steps=self%sgd_accepted_steps, &
                        &objective_initial=self%sgd_objective_initial, &
                        &objective_final=self%sgd_objective_final)
                else
                    cxy = self%grad_shsrch_obj%minimize_direct(irot=irot, xy_in=[0.,0.],&
                        &step_size=self%p_ptr%sgd_eta_shift, max_steps=self%p_ptr%sgd_shift_its,&
                        &sh_rot=.true., raw_euclid=.true., &
                        &accepted_steps=self%sgd_accepted_steps, &
                        &objective_initial=self%sgd_objective_initial, &
                        &objective_final=self%sgd_objective_final)
                endif
            else if( self%use_joint_angle )then
                ! Obtain the validated continuous-angle Stage 1 state in the
                ! native polar shift frame, then refine (sx,sy,theta) jointly.
                irot_stage1 = 0
                if( self%l_sh_first )then
                    cxy = self%grad_shsrch_obj%minimize(irot=irot_stage1, &
                        &xy_in=self%xy_first, sh_rot=.false., theta=theta_stage1)
                else
                    cxy = self%grad_shsrch_obj%minimize(irot=irot_stage1, &
                        &sh_rot=.false., theta=theta_stage1)
                endif
                if( irot_stage1 > 0 )then
                    joint_lims(1,:) = [-self%trs, self%trs]
                    joint_lims(2,:) = [-self%trs, self%trs]
                    joint_lims(3,:) = [theta_stage1-2., theta_stage1+2.]
                    call self%grad_shsrch_joint_obj%set_indices(self%best_class, self%iptcl)
                    call self%grad_shsrch_joint_obj%set_limits(joint_lims)
                    irot_joint = irot_stage1
                    cxy_joint = self%grad_shsrch_joint_obj%minimize_joint(irot_joint, &
                        &theta_stage1, cxy(2:3), sh_rot=.true., theta=theta_joint)
                    if( irot_joint > 0 )then
                        irot = irot_joint
                        cxy = cxy_joint
                        iang = modulo(nint(theta_joint)-1,self%nrots)+1
                        self%best_inpl_e3 = 360. - (self%b_ptr%pftc%get_rot(iang) + &
                            &(real(theta_joint)-real(iang))*self%b_ptr%pftc%get_dang())
                        self%best_inpl_e3_valid = .true.
                    else
                        irot = irot_stage1
                        call rotmat2d(self%b_ptr%pftc%get_rot(irot_stage1), rotmat)
                        cxy(2:) = matmul(cxy(2:), rotmat)
                    endif
                endif
            else
                if( .not.self%grad_shsrch_obj%does_opt_angle() )then
                    ! shift-only optimization
                    irot = self%best_rot
                endif
                if( self%l_sh_first )then
                    cxy = self%grad_shsrch_obj%minimize(irot=irot, xy_in=self%xy_first)
                else
                    cxy = self%grad_shsrch_obj%minimize(irot=irot)
                endif
            endif
            if( irot > 0 )then
                if( self%p_ptr%l_sgd_streaming_active )then
                    ! cxy(1) is the raw-loss merit -L; retain that sign until
                    ! the established upstream corr representation is written.
                    self%best_corr = real(cxy(1))
                else
                    self%best_corr = cxy(1)
                endif
                self%best_rot   = irot
                self%best_shvec = cxy(2:3)
            endif
        endif
    end subroutine inpl_srch

    subroutine inpl_srch_peaks( self, npeaks_inpl )
        class(strategy2D_srch), intent(inout) :: self
        integer,                intent(in)    :: npeaks_inpl
        real    :: sorted_cls_corrs(self%nrefs), cxy(3)
        integer :: sorted_cls_inds(self%nrefs), saved_inpl_inds(self%nrefs), iref, isample, inpl_ind
        if( .not. s2D%do_inplsrch(self%iptcl_batch) ) return
        sorted_cls_corrs         = s2D%class_space_corrs(   :, self%ithr)
        saved_inpl_inds          = s2D%class_space_inplinds(:, self%ithr)
        sorted_cls_inds          = (/(iref,iref=1,self%nrefs)/)
        call hpsort(sorted_cls_corrs, sorted_cls_inds)
        if( self%p_ptr%sgd_diagnostic )then
            write(logfhandle,'(A,1X,I0,1X,A,I0,1X,A,I0)') &
                '>>> SEARCH DIAG: inpl_srch_peaks entry; particle=', self%iptcl, &
                'nsolns=', self%nsolns, 'npeaks=', npeaks_inpl
            write(logfhandle,'(A,1X,I0,1X,A,I0)') &
                '>>> SEARCH DIAG: saved valid rotations=', count(saved_inpl_inds > 0), &
                'of=', self%nrefs
        endif
        ! reset class-space arrays so only shift-refined entries will be valid
        s2D%class_space_corrs(   :,self%ithr) = -1.
        s2D%class_space_e3s(     :,self%ithr) = 0.
        s2D%class_space_inplinds(:,self%ithr) = 0
        do isample = self%nrefs-npeaks_inpl+1,self%nrefs
            iref     = sorted_cls_inds(isample)
            inpl_ind = saved_inpl_inds(iref)
            if( s2D%cls_pops(iref) == 0 ) cycle
            if( inpl_ind == 0 ) cycle
            call self%grad_shsrch_obj2%set_indices(iref, self%iptcl)
            if( self%l_sh_first )then
                cxy = self%grad_shsrch_obj2%minimize(irot=inpl_ind, xy_in=self%xy_first)
                if( inpl_ind == 0 )then
                    inpl_ind = saved_inpl_inds(iref)
                    cxy(1)   = real(self%b_ptr%pftc%gen_corr_for_rot_8(iref, self%iptcl, real(self%xy_first,dp), inpl_ind))
                    cxy(2:3) = self%xy_first_rot
                endif
            else
                cxy = self%grad_shsrch_obj2%minimize(irot=inpl_ind)
                if( inpl_ind == 0 )then
                    inpl_ind = saved_inpl_inds(iref)
                    cxy      = [real(self%b_ptr%pftc%gen_corr_for_rot_8(iref, self%iptcl, inpl_ind)), 0.,0.]
                endif
            endif
            call self%store_solution(iref, inpl_ind, cxy(1))
        enddo
        if( self%p_ptr%sgd_diagnostic )then
            write(logfhandle,'(A,1X,I0,1X,A,I0)') &
                '>>> SEARCH DIAG: inpl_srch_peaks exit valid rotations=', &
                count(s2D%class_space_inplinds(:,self%ithr) > 0), 'of=', self%nrefs
        endif
    end subroutine inpl_srch_peaks

    subroutine store_solution( self, ref, inpl_ind, corr )
        class(strategy2D_srch), intent(inout) :: self
        integer,                intent(in)    :: ref, inpl_ind
        real,                   intent(in)    :: corr
        if( inpl_ind < 1 .or. inpl_ind > self%nrots )then
            if( self%p_ptr%sgd_diagnostic )then
                write(logfhandle,'(A,1X,I0,1X,A,I0,1X,A,I0,1X,A,I0)') &
                    '>>> SEARCH SAFETY: invalid stored rotation; particle=', self%iptcl, &
                    'reference=', ref, 'rotation=', inpl_ind, 'nrots=', self%nrots
            endif
            return
        endif
        if( s2D%class_space_corrs(ref, self%ithr) <= -huge(1.0)/2.0 )then
            self%nsolns = self%nsolns + 1
        elseif( corr <= s2D%class_space_corrs(ref, self%ithr) )then
            return
        endif
        s2D%class_space_inplinds(ref,self%ithr) = inpl_ind
        s2D%class_space_e3s(     ref,self%ithr) = 360. - self%b_ptr%pftc%get_rot(inpl_ind)
        s2D%class_space_corrs(   ref,self%ithr) = corr
    end subroutine store_solution

    subroutine assign_ori( self, os )
        class(strategy2D_srch), intent(in)    :: self
        class(oris),            intent(inout) :: os
        real :: dist, mat(2,2), u(2), x1(2), x2(2)
        real :: e3, mi_class, frac
        real :: best_corr_local
        integer :: iref, best_class_local, best_rot_local
        logical :: found_valid
        integer :: ithr
        ithr = self%ithr
        found_valid = .false.
        best_corr_local = -huge(1.0)
        if( self%prev_rot   <= 0 )then
            if( self%p_ptr%sgd_diagnostic )then
                write(logfhandle,'(A,1X,I0)') '>>> SEARCH SAFETY: assign_ori invalid previous rotation=', self%prev_rot
            endif
            THROW_HARD('Previous in-plane rotation index is invalid, cannot assign orientation.')
        endif
        if( self%prev_class <= 0 )then
            if( self%p_ptr%sgd_diagnostic )then
                write(logfhandle,'(A,1X,I0)') '>>> SEARCH SAFETY: assign_ori invalid previous class=', self%prev_class
            endif
            THROW_HARD('Previous in-plane class index is invalid, cannot assign orientation.')
        endif
        best_class_local = self%prev_class
        best_rot_local   = self%prev_rot
        do iref = 1, self%nrefs
            if( s2D%cls_pops(iref) <= 0 ) cycle
            if( s2D%class_space_inplinds(iref,ithr) <= 0 ) cycle
            if( s2D%class_space_corrs(iref,ithr) > best_corr_local )then
                best_corr_local  = s2D%class_space_corrs(iref,ithr)
                best_class_local = iref
                best_rot_local   = s2D%class_space_inplinds(iref,ithr)
                found_valid      = .true.
            endif
        enddo
        if( .not. found_valid )then
            best_class_local = self%prev_class
            best_rot_local   = self%prev_rot
            best_corr_local  = self%prev_corr
        endif
        ! Keep the legacy integer rotation and continuous metadata coupled.
        ! Only the explicit Euclidean opt-in gateway may replace the grid angle
        ! with a continuous value; all other paths reconstruct e3 from inpl.
        if( self%best_inpl_e3_valid .and. best_class_local == self%best_class .and. &
            &best_rot_local == self%best_rot )then
            e3 = self%best_inpl_e3
        else
            e3 = 360. - self%b_ptr%pftc%get_rot(best_rot_local)
        endif
        ! calculate in-plane rot dist (radians)
        if( self%l_fresh_start )then
            dist = 0.
        else
            u(1) = 0.
            u(2) = 1.
            call rotmat2d(e3, mat)
            x1   = matmul(u,mat)
            call rotmat2d(os%e3get(self%iptcl), mat)
            x2   = matmul(u,mat)
            dist = myacos(dot_product(x1,x2))
        endif
        ! calculate overlap between distributions
        mi_class = 0.
        if( (.not. self%l_fresh_start) .and. self%prev_class_mi == best_class_local ) mi_class = 1.
        ! search space explored
        frac = 100.*(real(self%nrefs_eval)/real(self%nrefs))
        ! update parameters
        call os%e3set(self%iptcl, e3)
        call os%set_shift(self%iptcl, self%prev_shvec + self%best_shvec)
        call os%set(self%iptcl, 'shincarg',   arg(self%best_shvec))
        call os%set(self%iptcl, 'inpl',       real(best_rot_local))
        call os%set(self%iptcl, 'class',      real(best_class_local))
        if( self%p_ptr%l_sgd_streaming_active )then
            ! Selection/refinement stores the merit -L so all comparisons keep
            ! their historical "larger is better" direction.  At the project
            ! boundary restore SIMPLE's established Euclidean score S=exp(-L)
            ! used by the legacy gen_euclids path; do not invent a stream-only
            ! score representation.
            call os%set(self%iptcl, 'corr',       exp(best_corr_local))
        else
            call os%set(self%iptcl, 'corr', best_corr_local)
        endif
        call os%set(self%iptcl, 'dist_inpl',  rad2deg(dist))
        call os%set(self%iptcl, 'mi_class',   mi_class)
        call os%set(self%iptcl, 'frac',       frac)
    end subroutine assign_ori

    subroutine kill( self )
        class(strategy2D_srch), intent(inout) :: self
        call self%grad_shsrch_obj%kill
        call self%grad_shsrch_obj2%kill
        call self%grad_shsrch_first_obj%kill
        call self%grad_shsrch_joint_obj%kill
    end subroutine kill

end module simple_strategy2D_srch
