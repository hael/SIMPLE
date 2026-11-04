!@descr: the polar shift map of the classic (non-table) searches: a particle's box, the screen of its candidates and the polish from the best cells
! The dense table's map path (process_particle_with_map) for the SHC, greedy and neighbourhood strategies: the
! candidates scored at the stored shift, the pool of the best by that score each given a coarse map over the
! particle's box (the disk or the window, shift_needs_disk), the best cells polished with minimize_in_cell.
module simple_shiftmap_srch
use simple_core_module_api,    only: dp
use simple_error,              only: simple_exception
use simple_defs,               only: SHC_INPL_TRSHWDTH
use simple_srch_sort_loc,      only: hpsort, minnloc
use simple_parameters,         only: parameters
use simple_builder,            only: builder
use simple_polarft_calc,       only: SHIFTMAP_SCREEN, SHIFTMAP_WINDOW, NSHIFTMAPS
use simple_pftc_shsrch_grad,   only: pftc_shsrch_grad
use simple_eul_prob_tab_utils, only: update_screen_best, shift_needs_disk, disk_radius_px, posterior_temperature, screen_pool_size,&
    &eulprob_dist_switch
implicit none

public :: shiftmap_srch, setup_shiftmap_srch
private
#include "simple_local_flags.inc"

! the run's maps, built once before the search's parallel region (setup_shiftmap_srch)
real(dp) :: prior_w_run = 0.d0
real     :: step_run(NSHIFTMAPS) = 0.
integer  :: kscreen_run = 0

type :: shiftmap_srch
    private
    class(parameters), pointer :: p_ptr => null()
    class(builder),    pointer :: b_ptr => null()
    type(pftc_shsrch_grad) :: gso_win              !< the polish within the window
    type(pftc_shsrch_grad) :: gso_disk             !< the polish within the disk
    real(dp) :: temp_i     = 0.d0                  !< the tempered temperature of the full band (the pool's size)
    real(dp) :: temp_s     = 0.d0                  !< the tempered temperature of the coarse band (the maps)
    real(dp) :: prior_w    = 0.d0                  !< the prior's weight: the disk at a particle's first centring
    integer  :: iptcl      = 0
    integer  :: imap       = SHIFTMAP_WINDOW
    logical  :: l_disk     = .false.
    integer, allocatable :: res_ref(:)             !< the polished candidates of the particle
    integer, allocatable :: res_irot(:)
    real,    allocatable :: res_corr(:), res_sh(:,:), res_sdev(:)
    logical, allocatable :: res_valid(:)           !< the polish or its cell beat the candidate's stored-shift score
    integer  :: nres       = 0
    logical  :: exists     = .false.
contains
    procedure :: new
    procedure :: prep
    procedure :: on_disk
    procedure :: srch
    procedure :: get_result
    procedure :: kill
end type shiftmap_srch

contains

    !> The two maps the classic searches use, on the calculator, and the run's constants: the screen over the disk of
    !! three sigma and the window of trs, both on the coarse band (setup_shiftmaps)
    subroutine setup_shiftmap_srch( params, build )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        real(dp) :: sigma_px, radius_disk
        if( .not. (params%prob_beta > 0.) ) THROW_HARD('the shift map screens at prob_beta x T: prob_beta must be positive')
        sigma_px    = real(params%shift_sigma, dp) / real(params%smpd, dp)
        radius_disk = real(disk_radius_px(params%shift_sigma, params%smpd, params%box), dp)
        call build%pftc%setup_shiftmaps(sigma_px, radius_disk, real(params%trs, dp), prior_w_run, kscreen_run)
        step_run(SHIFTMAP_SCREEN) = real(norm2(build%pftc%shiftmap_shift(SHIFTMAP_SCREEN, 2, 1) -&
            &build%pftc%shiftmap_shift(SHIFTMAP_SCREEN, 1, 1)))
        step_run(SHIFTMAP_WINDOW) = real(norm2(build%pftc%shiftmap_shift(SHIFTMAP_WINDOW, 2, 1) -&
            &build%pftc%shiftmap_shift(SHIFTMAP_WINDOW, 1, 1)))
    end subroutine setup_shiftmap_srch

    subroutine new( self, params, build )
        class(shiftmap_srch),      intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        class(builder),    target, intent(in)    :: build
        real :: lims(2,2), lims_init(2,2), radius_disk
        call self%kill
        self%p_ptr => params
        self%b_ptr => build
        lims_init(:,1) = -SHC_INPL_TRSHWDTH
        lims_init(:,2) =  SHC_INPL_TRSHWDTH
        lims(:,1)      = -params%trs
        lims(:,2)      =  params%trs
        call self%gso_win%new_alternating(build, lims, lims_init=lims_init, maxits=params%maxits_sh, coarse_init=.true.)
        radius_disk    = disk_radius_px(params%shift_sigma, params%smpd, params%box)
        lims(:,1)      = -radius_disk
        lims(:,2)      =  radius_disk
        call self%gso_disk%new_alternating(build, lims, lims_init=lims_init, maxits=params%maxits_sh, coarse_init=.true.)
        self%exists = .true.
    end subroutine new

    !> the particle's box (the disk when shift_needs_disk, under the prior at its first centring) and the
    !! temperatures of its screen
    subroutine prep( self, iptcl )
        class(shiftmap_srch), intent(inout) :: self
        integer,              intent(in)    :: iptcl
        real    :: shsig
        integer :: pftsz
        if( kscreen_run < 1 ) THROW_HARD('setup_shiftmap_srch first; shiftmap_srch%prep')
        self%iptcl   = iptcl
        shsig        = self%b_ptr%spproj_field%get(iptcl, 'shsig')
        self%l_disk  = shift_needs_disk(shsig, self%p_ptr%trs)
        self%imap    = SHIFTMAP_WINDOW
        self%prior_w = 0.d0
        if( self%l_disk )then
            self%imap = SHIFTMAP_SCREEN
            if( shsig < 0. ) self%prior_w = prior_w_run
        endif
        pftsz       = self%b_ptr%pftc%get_pftsz()
        self%temp_i = real(self%p_ptr%prob_beta, dp) *&
            &real(posterior_temperature(self%b_ptr%pftc%wsqsum_ptcl_now(iptcl), pftsz), dp)
        self%temp_s = real(self%p_ptr%prob_beta, dp) *&
            &real(posterior_temperature(self%b_ptr%pftc%wsqsum_ptcl_band(iptcl, kscreen_run), pftsz), dp)
        self%nres   = 0
    end subroutine prep

    pure logical function on_disk( self )
        class(shiftmap_srch), intent(in) :: self
        on_disk = self%l_disk
    end function on_disk

    !> The screen and the polish of candidates refs scored at the stored shift (corrs0, larger is better): a map each
    !! for the pool of the best by that score (all of them up to nbest, else at least twice nbest, the origin
    !! posterior's mass quantile at the particle's tempered temperature), the nbest best cells polished in a box of
    !! one map step. The polished candidates are kept for get_result
    subroutine srch( self, refs, corrs0, nbest )
        class(shiftmap_srch), intent(inout) :: self
        integer,              intent(in)    :: refs(:), nbest
        real,                 intent(in)    :: corrs0(:)
        real,    allocatable :: dists(:), scr_loss(:), scr_xy(:,:), scr_sdev(:)
        integer, allocatable :: order(:), scr_irot(:), best(:)
        real     :: cxy(3), loss, width_best(max(1,nbest))
        integer  :: n, npool, nb, k, j, ref, irot, irho, ialpha
        logical  :: l_polished
        self%nres = 0
        n = size(refs)
        if( n < 1 .or. nbest < 1 ) return
        if( size(corrs0) /= n ) THROW_HARD('one stored-shift score per candidate; shiftmap_srch%srch')
        ! the candidates by their stored-shift score, best first, and the pool
        allocate(dists(n), order(n))
        dists = eulprob_dist_switch(corrs0, self%p_ptr%cc_objfun)
        order = [(k, k = 1, n)]
        call hpsort(dists, order)
        npool = n
        if( n > nbest ) npool = min(n, max(2 * nbest, screen_pool_size(dists, n, self%temp_i)))
        ! the screen: each candidate's map and its mode under the box's prior
        width_best = huge(1.)
        allocate(scr_loss(npool), scr_xy(2,npool), scr_sdev(npool), scr_irot(npool))
        do k = 1, npool
            ref = refs(order(k))
            call self%b_ptr%pftc%gen_shiftmap(self%imap, ref, self%iptcl)
            call self%b_ptr%pftc%shiftmap_mode(self%imap, self%temp_s, self%prior_w, &
                &irot, irho, ialpha, loss, &
                &temp_phys=self%temp_s / real(self%p_ptr%prob_beta, dp), sdev=scr_sdev(k), &
                &width_cutoff=maxval(width_best))
            call update_screen_best(loss, width_best)
            scr_loss(k)  = loss
            scr_irot(k)  = irot
            scr_xy(:,k)  = real(self%b_ptr%pftc%shiftmap_shift(self%imap, irho, ialpha))
        enddo
        ! the polish of the screen's best
        nb = min(nbest, npool)
        allocate(best(nb))
        best = minnloc(scr_loss, nb)
        call alloc_results(nb)
        do j = 1, nb
            k    = best(j)
            ref  = refs(order(k))
            irot = scr_irot(k)
            if( self%l_disk )then
                call self%gso_disk%set_indices(ref, self%iptcl)
                call self%gso_disk%minimize_in_cell(irot, scr_xy(:,k), step_run(self%imap), cxy, l_polished)
            else
                call self%gso_win%set_indices(ref, self%iptcl)
                call self%gso_win%minimize_in_cell(irot, scr_xy(:,k), step_run(self%imap), cxy, l_polished)
            endif
            self%nres = self%nres + 1
            self%res_ref(self%nres)   = ref
            self%res_irot(self%nres)  = irot
            self%res_corr(self%nres)  = cxy(1)
            self%res_sh(:,self%nres)  = cxy(2:3)
            self%res_sdev(self%nres)  = scr_sdev(k)
            self%res_valid(self%nres) = l_polished .or. cxy(1) > corrs0(order(k))
        enddo
        deallocate(dists, order, scr_loss, scr_xy, scr_sdev, scr_irot, best)

    contains

        subroutine alloc_results( nmax )
            integer, intent(in) :: nmax
            if( allocated(self%res_ref) )then
                if( size(self%res_ref) >= nmax ) return
                deallocate(self%res_ref, self%res_irot, self%res_corr, self%res_sh, self%res_sdev, self%res_valid)
            endif
            allocate(self%res_ref(nmax), self%res_irot(nmax), self%res_corr(nmax), self%res_sh(2,nmax),&
                &self%res_sdev(nmax), self%res_valid(nmax))
        end subroutine alloc_results

    end subroutine srch

    !> the last srch's result for candidate ref: found says it was polished, valid that its pose (rotation, score,
    !! shift increment in the particle's axes) beat its stored-shift score; sdev the posterior width of its map
    subroutine get_result( self, ref, found, valid, irot, corr, sh, sdev )
        class(shiftmap_srch), intent(in)  :: self
        integer,              intent(in)  :: ref
        logical,              intent(out) :: found, valid
        integer,              intent(out) :: irot
        real,                 intent(out) :: corr, sh(2), sdev
        integer :: i
        found = .false.
        valid = .false.
        irot  = 0
        corr  = 0.
        sh    = 0.
        sdev  = -1.
        do i = 1, self%nres
            if( self%res_ref(i) /= ref ) cycle
            found = .true.
            valid = self%res_valid(i)
            irot  = self%res_irot(i)
            corr  = self%res_corr(i)
            sh    = self%res_sh(:,i)
            sdev  = self%res_sdev(i)
            return
        enddo
    end subroutine get_result

    subroutine kill( self )
        class(shiftmap_srch), intent(inout) :: self
        if( .not. self%exists ) return
        call self%gso_win%kill
        call self%gso_disk%kill
        if( allocated(self%res_ref) ) deallocate(self%res_ref, self%res_irot, self%res_corr, self%res_sh,&
            &self%res_sdev, self%res_valid)
        self%nres = 0
        nullify(self%p_ptr, self%b_ptr)
        self%exists = .false.
    end subroutine kill

end module simple_shiftmap_srch
