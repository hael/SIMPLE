! common strategy3D methods and type specification for polymorphic strategy3D object creation are delegated to this class
module simple_strategy3D_srch
include 'simple_lib.f08'
use simple_oris,               only: oris
use simple_ori,                only: ori
use simple_sym,                only: sym
use simple_pftcc_shsrch_grad,  only: pftcc_shsrch_grad  ! gradient-based in-plane angle and shift search
use simple_pftcc_orisrch_grad, only: pftcc_orisrch_grad ! gradient-based search over all df:s
use simple_polarft_corrcalc,   only: pftcc_glob, polarft_corrcalc
use simple_parameters,         only: params_glob
use simple_builder,            only: build_glob
use simple_strategy3D_alloc    ! singleton s3D
implicit none

public :: strategy3D_srch, strategy3D_spec, set_ptcl_stats, eval_ptcl
private
#include "simple_local_flags.inc"

logical, parameter :: DOCONTINUOUS = .false.

type strategy3D_spec
    integer, pointer :: symmat(:,:) => null()
    integer :: iptcl=0, szsn=0
    logical :: do_extr=.false.
    real    :: extr_score_thresh=0.
end type strategy3D_spec

type strategy3D_srch
    type(pftcc_shsrch_grad)  :: grad_shsrch_obj           !< origin shift search object, L-BFGS with gradient
    type(pftcc_orisrch_grad) :: grad_orisrch_obj          !< obj 4 search over all df:s, L-BFGS with gradient
    integer                  :: iptcl         = 0         !< global particle index
    integer                  :: ithr          = 0         !< thread index
    integer                  :: nrefs         = 0         !< total # references (nstates*nprojs)
    integer                  :: nrefsmaxinpl  = 0         !< total # references (nstates*nprojs)
    integer                  :: nnnrefs       = 0         !< total # neighboring references (nstates*nnn)
    integer                  :: nstates       = 0         !< # states
    integer                  :: nprojs        = 0         !< # projections
    integer                  :: nrots         = 0         !< # in-plane rotations in polar representation
    integer                  :: npeaks        = 0         !< # peaks (nonzero orientation weights)
    integer                  :: npeaks_eff    = 0         !< effective # peaks
    integer                  :: nsym          = 0         !< symmetry order
    integer                  :: nbetter       = 0         !< # better orientations identified
    integer                  :: nrefs_eval    = 0         !< # references evaluated
    integer                  :: nnn_static    = 0         !< # nearest neighbors (static)
    integer                  :: nnn           = 0         !< # nearest neighbors (dynamic)
    integer                  :: prev_roind    = 0         !< previous in-plane rotation index
    integer                  :: prev_state    = 0         !< previous state index
    integer                  :: prev_ref      = 0         !< previous reference index
    integer                  :: prev_proj     = 0         !< previous projection direction index
    real                     :: prev_corr     = 1.        !< previous best correlation
    real                     :: specscore     = 0.        !< spectral score
    real                     :: prev_shvec(2) = 0.        !< previous origin shift vector
    logical                  :: neigh         = .false.   !< nearest neighbour refinement flag
    logical                  :: doshift       = .true.    !< 2 indicate whether 2 serch shifts
    logical                  :: dowinpl       = .true.    !< 2 indicate weights over in-planes as well as projection dirs
    logical                  :: exists        = .false.   !< 2 indicate existence
  contains
    procedure :: new
    procedure :: prep4srch
    procedure :: inpl_srch
    procedure :: store_solution
    procedure :: kill
end type strategy3D_srch

contains

    ! class method: set_ptcl_stats for filling in stats
    ! of particles not part of update fraction
    subroutine set_ptcl_stats( pftcc, iptcl, pftcc_pind )
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: iptcl, pftcc_pind
        integer   :: prev_state, prev_roind, prev_proj, prev_ref
        type(ori) :: o_prev
        real      :: specscore
        prev_state = build_glob%spproj_field%get_state(iptcl)      ! state index
        if( prev_state == 0 )return
        ! previous parameters
        call build_glob%spproj_field%get_ori(iptcl, o_prev)
        prev_roind = pftcc%get_roind(360.-o_prev%e3get())          ! in-plane angle index
        prev_proj  = build_glob%eulspace%find_closest_proj(o_prev) ! previous projection direction
        prev_ref   = (prev_state-1)*params_glob%nspace + prev_proj ! previous reference
        ! particle prep
        call pftcc%prep_matchfilt(pftcc_pind,prev_ref,prev_roind)
        ! calc specscore
        specscore = pftcc%specscore(prev_ref, pftcc_pind, prev_roind)
        ! update spproj_field
        call build_glob%spproj_field%set(iptcl, 'specscore',  specscore)
        call o_prev%kill
    end subroutine set_ptcl_stats

    ! class method: evaluation of stats and objective function
    subroutine eval_ptcl( pftcc, iptcl, pftcc_pind )
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: iptcl, pftcc_pind
        integer   :: prev_state, prev_roind, prev_proj, prev_ref
        type(ori) :: o_prev
        real      :: specscore, corr, corrs(pftcc_glob%get_nrots())
        prev_state = build_glob%spproj_field%get_state(iptcl)      ! state index
        if( prev_state == 0 )return
        ! previous parameters
        call build_glob%spproj_field%get_ori(iptcl, o_prev)
        prev_roind = pftcc%get_roind(360.-o_prev%e3get())          ! in-plane angle index
        prev_proj  = build_glob%eulspace%find_closest_proj(o_prev) ! previous projection direction
        prev_ref   = (prev_state-1)*params_glob%nspace + prev_proj ! previous reference
        ! particle prep
        call pftcc%prep_matchfilt(pftcc_pind,prev_ref,prev_roind)
        ! calc specscore
        specscore = pftcc%specscore(prev_ref, pftcc_pind, prev_roind)
        ! prep corr
        call pftcc_glob%gencorrs(prev_ref, pftcc_pind, corrs)
        if( params_glob%cc_objfun == OBJFUN_EUCLID )then
            corr = maxval(corrs)
        else
            corr = max(0.,maxval(corrs))
        endif
        ! update spproj_field
        call build_glob%spproj_field%set(iptcl, 'proj',      real(prev_proj))
        call build_glob%spproj_field%set(iptcl, 'corr',      corr)
        call build_glob%spproj_field%set(iptcl, 'specscore', specscore)
        call build_glob%spproj_field%set(iptcl, 'w',         1.)
        call build_glob%spproj_field%set(iptcl, 'ow',        1.)
        call o_prev%kill
    end subroutine eval_ptcl

    subroutine new( self, spec, npeaks )
        class(strategy3D_srch), intent(inout) :: self
        class(strategy3D_spec), intent(in)    :: spec
        integer,                intent(in)    :: npeaks
        integer, parameter :: MAXITS = 60
        real    :: lims(2,2), lims_init(2,2)
        ! set constants
        self%iptcl        = spec%iptcl
        self%nstates      = params_glob%nstates
        self%nprojs       = params_glob%nspace
        self%nrefs        = self%nprojs*self%nstates
        self%nrefsmaxinpl = self%nrefs*params_glob%ninplpeaks
        self%nrots        = pftcc_glob%get_nrots()
        self%npeaks       = npeaks
        self%nbetter      = 0
        self%nrefs_eval   = 0
        self%nsym         = build_glob%pgrpsyms%get_nsym()
        self%doshift      = params_glob%l_doshift
        self%neigh        = params_glob%neigh == 'yes'
        self%nnn_static   = params_glob%nnn
        self%nnn          = params_glob%nnn
        self%nnnrefs      = self%nnn*self%nstates
        self%dowinpl      = npeaks /= 1
        ! create in-plane search object
        lims(:,1)         = -params_glob%trs
        lims(:,2)         =  params_glob%trs
        lims_init(:,1)    = -SHC_INPL_TRSHWDTH
        lims_init(:,2)    =  SHC_INPL_TRSHWDTH
        call self%grad_shsrch_obj%new(lims, lims_init=lims_init,&
            &shbarrier=params_glob%shbarrier, maxits=MAXITS, opt_angle=.not.self%dowinpl)
        ! create all df:s search object
        call self%grad_orisrch_obj%new
        self%exists = .true.
    end subroutine new

    subroutine prep4srch( self )
        class(strategy3D_srch), intent(inout) :: self
        integer   :: i, istate
        type(ori) :: o_prev
        real      :: corrs(self%nrots), corr
        ! previous parameters
        call build_glob%spproj_field%get_ori(self%iptcl, o_prev)
        self%prev_state = o_prev%get_state()                                ! state index
        self%prev_roind = pftcc_glob%get_roind(360.-o_prev%e3get())         ! in-plane angle index
        self%prev_shvec = o_prev%get_2Dshift()                              ! shift vector
        self%prev_proj  = build_glob%eulspace%find_closest_proj(o_prev)     ! previous projection direction
        self%prev_ref   = (self%prev_state-1)*self%nprojs + self%prev_proj  ! previous reference
        ! init threaded search arrays
        call prep_strategy3D_thread(self%ithr)
        ! search order
        if( self%neigh )then
            if( .not. allocated(build_glob%nnmat) )&
                &THROW_HARD('need optional nnmat to be present for refine=neigh modes; prep4srch')
            do istate = 0, self%nstates - 1
                i = istate * self%nnn + 1
                s3D%srch_order(self%ithr,i:i+self%nnn-1) = build_glob%nnmat(self%prev_proj,:) + istate*self%nprojs
            enddo
            call s3D%rts(self%ithr)%shuffle(s3D%srch_order(self%ithr,:))
        else
            call s3D%rts(self%ithr)%ne_ran_iarr(s3D%srch_order(self%ithr,:))
        endif
        call put_last(self%prev_ref, s3D%srch_order(self%ithr,:))
        ! sanity check
        if( self%prev_state > 0 )then
            if( self%prev_state > self%nstates )          THROW_HARD('previous best state outside boundary; prep4srch')
            if( .not. s3D%state_exists(self%prev_state) ) THROW_HARD('empty previous state; prep4srch')
        endif
        ! particle prep
        call pftcc_glob%prep_matchfilt(self%iptcl,self%prev_ref,self%prev_roind)
        ! calc specscore
        self%specscore = pftcc_glob%specscore(self%prev_ref, self%iptcl, self%prev_roind)
        ! prep corr
        call pftcc_glob%gencorrs(self%prev_ref, self%iptcl, corrs)
        if( params_glob%cc_objfun == OBJFUN_EUCLID )then
            corr = maxval(corrs)
        else
            corr = max(0.,maxval(corrs))
        endif
        self%prev_corr = corr
        call o_prev%kill
    end subroutine prep4srch

    subroutine inpl_srch( self )
        class(strategy3D_srch), intent(inout) :: self
        type(ori) :: o
        real      :: cxy(3)
        integer   :: i, j, ref, irot, cnt
        logical   :: found_better
        if( DOCONTINUOUS )then
            ! BFGS over all df:s
            call o%new
            cnt = 0
            do i=self%nrefs,self%nrefs-self%npeaks+1,-1
                cnt = cnt + 1
                ref = s3D%proj_space_refinds_sorted_highest(self%ithr, i)
                if( cnt <= CONTNPEAKS )then
                    ! continuous refinement over all df:s
                    call o%set_euler(s3D%proj_space_euls(self%ithr,ref,1,:))
                    call o%set_shift([0.,0.])
                    call self%grad_orisrch_obj%set_particle(self%iptcl)
                    cxy = self%grad_orisrch_obj%minimize(o, NPEAKSATHRES/2.0, params_glob%trs, found_better)
                    if( found_better )then
                        s3D%proj_space_euls(self%ithr, ref, 1,:) = o%get_euler()
                        s3D%proj_space_corrs(self%ithr,ref, 1)   = cxy(1)
                        s3D%proj_space_shift(self%ithr,ref, 1,:) = cxy(2:3)
                    endif
                else
                    ! refinement of in-plane rotation (discrete) & shift (continuous)
                    call self%grad_shsrch_obj%set_indices(ref, self%iptcl)
                    cxy = self%grad_shsrch_obj%minimize(irot=irot)
                    if( irot > 0 )then
                        ! irot > 0 guarantees improvement found, update solution
                        s3D%proj_space_euls(self%ithr, ref,1, 3) = 360. - pftcc_glob%get_rot(irot)
                        s3D%proj_space_corrs(self%ithr,ref,1)    = cxy(1)
                        s3D%proj_space_shift(self%ithr,ref,1,:)  = cxy(2:3)
                    endif
                endif
            end do
        else
            if( self%doshift )then
                if( self%dowinpl )then
                    ! BFGS over shifts only
                    do i=self%nrefsmaxinpl,self%nrefsmaxinpl-self%npeaks+1,-1
                        ref  = s3D%proj_space_refinds_sorted(self%ithr, i)
                        if( .not. s3D%proj_space_corrs_srchd(self%ithr,ref) ) cycle ! must have seen the reference before
                        call self%grad_shsrch_obj%set_indices(ref, self%iptcl)
                        j    = s3D%proj_space_inplinds_sorted(self%ithr, i)
                        irot = s3D%proj_space_inplinds(self%ithr, ref, j)
                        cxy  = self%grad_shsrch_obj%minimize(irot=irot)
                        if( irot > 0 )then
                            ! irot > 0 guarantees improvement found, update solution
                            s3D%proj_space_euls( self%ithr,ref,j,3) = 360. - pftcc_glob%get_rot(irot)
                            s3D%proj_space_corrs(self%ithr,ref,j)   = cxy(1)
                            s3D%proj_space_shift(self%ithr,ref,j,:) = cxy(2:3)
                        endif
                    end do
                else
                    ! BFGS over shifts with in-plane rot exhaustive callback
                    do i=self%nrefs,self%nrefs-self%npeaks+1,-1
                        ref = s3D%proj_space_refinds_sorted_highest(self%ithr, i)
                        if( .not. s3D%proj_space_corrs_srchd(self%ithr,ref) ) cycle ! must have seen the reference before
                        call self%grad_shsrch_obj%set_indices(ref, self%iptcl)
                        cxy = self%grad_shsrch_obj%minimize(irot=irot)
                        if( irot > 0 )then
                            ! irot > 0 guarantees improvement found, update solution
                            s3D%proj_space_euls( self%ithr,ref,1,3) = 360. - pftcc_glob%get_rot(irot)
                            s3D%proj_space_corrs(self%ithr,ref,1)   = cxy(1)
                            s3D%proj_space_shift(self%ithr,ref,1,:) = cxy(2:3)
                        endif
                    end do
                endif
            endif
        endif
    end subroutine inpl_srch

    subroutine store_solution( self, ref, inpl_inds, corrs, searched )
        class(strategy3D_srch), intent(inout) :: self
        integer,                intent(in)    :: ref, inpl_inds(params_glob%ninplpeaks)
        real,                   intent(in)    :: corrs(params_glob%ninplpeaks)
        logical,                intent(in)    :: searched
        integer :: inpl
        s3D%proj_space_inplinds(self%ithr,ref,:) = inpl_inds
        do inpl=1,params_glob%ninplpeaks
            s3D%proj_space_euls(self%ithr,ref,inpl,3) = 360. - pftcc_glob%get_rot(inpl_inds(inpl))
        end do
        s3D%proj_space_corrs(self%ithr,ref,:) = corrs
        s3D%proj_space_corrs_calcd(self%ithr,ref) = .true.
        if (searched) s3D%proj_space_corrs_srchd(self%ithr,ref) = .true.
    end subroutine store_solution

    subroutine kill( self )
        class(strategy3D_srch), intent(inout) :: self
        call self%grad_shsrch_obj%kill
    end subroutine kill

end module simple_strategy3D_srch
