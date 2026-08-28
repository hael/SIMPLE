!@descr: Two-gate agglomerative merge of over-provisioned flex_pca states.
!
!        Latent criteria (BIC, ICL) score the mixture in the space it was fitted in, so components
!        that partition particles but reconstruct to the same density still score as distinct. The
!        question here -- how many distinguishable maps the data support -- is therefore asked of
!        the orientations and of the maps instead.
!
!        Gate 1, orientation (cheap): the latent comes from 2D projections, so viewing direction is
!        its largest confound and a component whose viewing-axis distribution departs from the
!        global one is a view cluster rather than a state. Gate 1 fires on failure of a null, not on
!        improvement -- pooling always improves view coverage, so an improvement rule would run to
!        K=1.
!
!        Gate 2, volume (expensive, pairs surviving gate 1): catches several components landing on
!        one true state with full view coverage each, which gate 1 cannot see. Tested against each
!        state's own half-map reproducibility.
!
!        Single pass: the pairwise relation is computed once on the delivered maps and closed under
!        COMPLETE linkage; the caller re-reconstructs once at the converged count.
!
!        Off by default; SIMPLE_COV_MERGE=1 enables, and npreimages=0 (auto state ceiling) enables it
!        implicitly. SIMPLE_COV_MERGE=0 always wins.
module simple_flex_pca_merge
use simple_core_module_api
use simple_image,          only: image
use simple_parameters,     only: parameters
use simple_flex_pca_rec3D, only: flex_rec_smpd
implicit none
private
#include "simple_local_flags.inc"

public :: flex_pca_merge_enabled, flex_pca_merge_force_on, two_gate_state_merge

!> Set by the npreimages=0 auto-ceiling path; consulted only when SIMPLE_COV_MERGE is unset.
logical :: l_merge_forced = .false.

!> Viewing-axis second moments are symmetric 3x3 with unit trace (v is a unit vector), leaving
!! five free parameters.
integer,  parameter :: VIEW_DOF = 5
!> Upper tail of chi^2_5 at p = 0.001. Necessary but not sufficient: at realistic neff the test is
!! over-powered and every state clears it, since no experimental view distribution is exactly the
!! global one. It only guards the sparse states, where a departure really could be sampling.
real(dp), parameter :: VIEW_CHI2_CRIT = 20.515d0
!> The decision is therefore made on effect size, chi2/neff, which does not grow with N. The cut is
!! robust and relative to this dataset's own states (median + k*MAD): a specimen with preferred
!! orientation gives every state a non-global view distribution, so what marks a view cluster is
!! standing out from its peers. Gate 2 covers the case where every state is a view cluster and
!! there are no outliers to find.
real(dp), parameter :: VIEW_MAD_K = 3.0d0
!> MAD -> sigma for a normal, so VIEW_MAD_K is read in familiar units.
real(dp), parameter :: MAD2SIGMA = 1.4826d0
!> A view-clustered state is folded into its best map match only if it substantially resembles it.
!! Without the floor gate 1 folds on the best available match however poor it is, tipping particles
!! into maps they do not resemble. Below the floor the state is reported and kept: a component that
!! is view-driven and resembles nothing is a finding, not something to hide in an unrelated class.
real(dp), parameter :: VIEW_FOLD_MIN_R = 0.8d0
!> Shells where both states reproduce below this are noise in both and carry no evidence either way.
real(dp), parameter :: FSC_SIGNAL_FLOOR = 0.143d0
!> Default on the disattenuated ratio; 1 means the two maps agree as well as each agrees with
!! itself, i.e. indistinguishable given the noise. SIMPLE_COV_MERGE_R overrides.
!> Gate on the DEVIATION-from-ensemble-mean disattenuated ratio (see pair_map_ratio). 0.98 is the
!! measured value that fuses only true duplicates on both benchmarks at once: Ribosembly
!! same-GT pairs bottom out at 0.983 while different-GT pairs top out at 0.970, and EMPIAR-10076
!! cross-family pairs top out at 0.878. It was 0.95 when the ratio was computed on RAW maps,
!! where the shared consensus density made every pair score high and no gate transferred.
real(dp), parameter :: MERGE_R_DEFAULT = 0.98d0
!> Absolute SANITY bound for the adaptive gate -- not a duplicate threshold. Two states whose
!! deviations correlate below this are never duplicates whatever the gap structure says. It must
!! stay well BELOW the duplicate cluster of any dataset: measured, that cluster sits at 0.975+ on
!! Ribosembly but at 0.88-0.98 on EMPIAR-10076, so a floor at 0.90 would silently disable
!! adaptation on 10076 (it did) -- the cut has to be free to follow the data.
real(dp), parameter :: MERGE_R_FLOOR   = 0.50d0
!> A break is a real cluster boundary only if it dwarfs the typical spacing among the candidates.
real(dp), parameter :: MERGE_GAP_K     = 4.0d0
!> ... and only if what it separates stays a small minority of all pairs.
real(dp), parameter :: MERGE_TOP_FRAC  = 0.15d0


contains


    logical function flex_pca_merge_enabled() result( on )
        character(len=32) :: envval
        integer :: stat, ln
        on = l_merge_forced
        call get_environment_variable('SIMPLE_COV_MERGE', envval, ln, stat)
        if( stat /= 0 .or. ln < 1 )then
            ! unset: the auto-provisioned state ceiling turns it on, since a ceiling without the
            ! collapse is just a large state count
            on = l_merge_forced
            return
        endif
        ! an explicit SIMPLE_COV_MERGE=0 wins over the force: the user asked for no merge
        on = trim(adjustl(envval)) /= '0'
    end function flex_pca_merge_enabled

    !> Turn the merge on for runs that provision a state CEILING rather than a state count.
    !! Never overrides an explicit SIMPLE_COV_MERGE=0.
    subroutine flex_pca_merge_force_on()
        l_merge_forced = .true.
    end subroutine flex_pca_merge_force_on

    subroutine merge_env_dp( name, val )
        character(len=*), intent(in)    :: name
        real(dp),         intent(inout) :: val
        character(len=32) :: envval
        integer  :: stat, ln, io_stat
        real(dp) :: tmp
        call get_environment_variable(name, envval, ln, stat)
        if( stat /= 0 .or. ln < 1 ) return
        read(envval, *, iostat=io_stat) tmp
        if( io_stat == 0 ) val = tmp
    end subroutine merge_env_dp

    !> .true. when the environment variable is set to exactly `want` (case-sensitive, trimmed)
    logical function merge_env_is( name, want ) result( yes )
        character(len=*), intent(in) :: name, want
        character(len=32) :: envval
        integer :: stat, ln
        yes = .false.
        call get_environment_variable(name, envval, ln, stat)
        if( stat /= 0 .or. ln < 1 ) return
        yes = trim(adjustl(envval)) == want
    end function merge_env_is

    !> Chi-squared departure of each state's viewing-axis second moment from the global one.
    !!
    !! T_k = sum_i w_ik v_i v_i' / sum_i w_ik against T = sum_i v_i v_i' / nptcls. Under "this
    !! state's views are a random sample of the global distribution" the entries of T_k - T have
    !! variance Var_i(v_q v_r) / neff_k, so neff_k * sum (T_k - T)^2 / Var is chi^2 on VIEW_DOF.
    !! The null self-calibrates in neff, so a sparse state gets a wide null and is not called
    !! contaminated.
    subroutine view_coverage_chi2( views, weights, nptcls, nstates, chi2, neff, effsz )
        integer,  intent(in)  :: nptcls, nstates
        real(dp), intent(in)  :: views(3,nptcls)
        real,     intent(in)  :: weights(nptcls,nstates)
        real(dp), intent(out) :: chi2(nstates), neff(nstates)
        !> chi2/neff: mean squared view departure per particle, the N-free effect size gate 1 cuts on
        real(dp), intent(out) :: effsz(nstates)
        real(dp) :: Tbar(3,3), Tk(3,3), Vvar(3,3), wsum, w2sum, w, d
        integer  :: i, q, r, state
        Tbar = 0.d0
        do i = 1, nptcls
            do q = 1, 3
                do r = 1, 3
                    Tbar(q,r) = Tbar(q,r) + views(q,i)*views(r,i)
                end do
            end do
        end do
        Tbar = Tbar / real(nptcls,dp)
        ! per-entry variance of v_q v_r over the WHOLE particle set: the null's scale
        Vvar = 0.d0
        do i = 1, nptcls
            do q = 1, 3
                do r = 1, 3
                    Vvar(q,r) = Vvar(q,r) + (views(q,i)*views(r,i) - Tbar(q,r))**2
                end do
            end do
        end do
        Vvar = Vvar / real(nptcls,dp)
        do state = 1, nstates
            wsum  = 0.d0
            w2sum = 0.d0
            Tk    = 0.d0
            do i = 1, nptcls
                w = real(weights(i,state),dp)
                if( w <= 0.d0 ) cycle
                wsum  = wsum  + w
                w2sum = w2sum + w*w
                do q = 1, 3
                    do r = 1, 3
                        Tk(q,r) = Tk(q,r) + w*views(q,i)*views(r,i)
                    end do
                end do
            end do
            if( wsum <= DTINY .or. w2sum <= DTINY )then
                chi2(state)  = 0.d0
                neff(state)  = 0.d0
                effsz(state) = 0.d0
                cycle
            endif
            ! Kish effective sample size: soft responsibilities are not a count of particles.
            neff(state) = wsum*wsum / w2sum
            Tk          = Tk / wsum
            d = 0.d0
            do q = 1, 3
                do r = q, 3          ! symmetric: upper triangle only, else every off-diagonal counts twice
                    if( Vvar(q,r) <= DTINY ) cycle
                    d = d + (Tk(q,r) - Tbar(q,r))**2 / Vvar(q,r)
                end do
            end do
            chi2(state)  = neff(state) * d
            ! excess over the sampling floor, not the raw ratio: under the null E[chi2/neff] =
            ! VIEW_DOF/neff, so the raw ratio penalises sparse states when neff spans a wide range.
            ! A state at exactly the null contributes zero.
            effsz(state) = max(0.d0, (chi2(state) - real(VIEW_DOF,dp)) / neff(state))
        end do
    end subroutine view_coverage_chi2

    !> Disattenuated shell correlation between two states' maps.
    !!
    !! Two maps built from the same particle set with different weights share noise, so a
    !! same-halfset cross-correlation is inflated toward agreement. Every correlation is therefore
    !! taken across halfsets (even_s vs odd_t, odd_s vs even_t) so all four spectra involve disjoint
    !! particles. If s and t are the same underlying map the cross spectrum is the common signal
    !! attenuated by each map's own reliability, so C_st / sqrt(C_ss*C_tt) is 1; below 1 they differ.
    subroutine pair_map_ratio( evols, ovols, nstates, nshell, R, Rmin, lstop )
        integer,     intent(in)    :: nstates, nshell
        type(image), intent(inout) :: evols(nstates), ovols(nstates)
        real(dp),    intent(out)   :: R(nstates,nstates)
        !> min of the TWO half-independent cross estimates (e_s vs o_t and o_s vs e_t), each
        !! disattenuated by the same reliabilities. R is their mean, so 2*(R - Rmin) is the
        !! halfset disagreement — a per-pair noise scale for the merge decision.
        real(dp), optional, intent(out) :: Rmin(nstates,nstates)
        !> optional band cap (shell index): conformations are low-frequency objects, and the
        !! per-shell signal floor alone can admit mid-band shells that carry reconstruction
        !! noise rather than conformation (validated externally at ~25 A on 10076/10049)
        integer,  optional, intent(in)  :: lstop
        real,     allocatable :: css(:,:), cst(:), cts(:)
        real(dp) :: num_a, num_b, den, ratio, ratio_a, ratio_b
        integer  :: s, t, l, nval, lmax
        allocate(css(nshell,nstates), source=0.)
        allocate(cst(nshell), cts(nshell), source=0.)
        do s = 1, nstates
            call evols(s)%fsc(ovols(s), css(:,s))
        end do
        R = 1.d0
        if( present(Rmin) ) Rmin = 1.d0
        do s = 1, nstates - 1
            do t = s + 1, nstates
                call evols(s)%fsc(ovols(t), cst)
                call ovols(s)%fsc(evols(t), cts)
                num_a = 0.d0
                num_b = 0.d0
                den   = 0.d0
                nval  = 0
                lmax  = nshell
                if( present(lstop) ) lmax = max(2, min(nshell, lstop))
                do l = 1, lmax
                    ! a shell where either state is already noise carries no evidence
                    if( real(css(l,s),dp) <= FSC_SIGNAL_FLOOR ) cycle
                    if( real(css(l,t),dp) <= FSC_SIGNAL_FLOOR ) cycle
                    num_a = num_a + real(cst(l),dp)
                    num_b = num_b + real(cts(l),dp)
                    den   = den + sqrt(real(css(l,s),dp)*real(css(l,t),dp))
                    nval  = nval + 1
                end do
                if( nval < 2 .or. den <= DTINY )then
                    ratio   = 0.d0      ! no shared resolution range: not evidence to merge
                    ratio_a = 0.d0
                    ratio_b = 0.d0
                else
                    ratio_a = num_a / den
                    ratio_b = num_b / den
                    ratio   = 0.5d0*(ratio_a + ratio_b)
                endif
                R(s,t) = ratio
                R(t,s) = ratio
                if( present(Rmin) )then
                    Rmin(s,t) = min(ratio_a, ratio_b)
                    Rmin(t,s) = Rmin(s,t)
                endif
            end do
        end do
        deallocate(css, cst, cts)
    end subroutine pair_map_ratio

    !> Two-gate merge. Returns the state each input state was merged into (1..nstates_out).
    subroutine two_gate_state_merge( params, views, weights, nptcls, nstates, label_out, nstates_out )
        class(parameters), intent(in)  :: params
        integer,           intent(in)  :: nptcls, nstates
        real(dp),          intent(in)  :: views(3,nptcls)
        real,              intent(in)  :: weights(nptcls,nstates)
        integer,           intent(out) :: label_out(nstates)
        integer,           intent(out) :: nstates_out
        type(image), allocatable :: evols(:), ovols(:)
        type(image) :: mskwarm
        real(dp),    allocatable :: chi2(:), neff(:), effsz(:), work(:), R(:,:), Rmin(:,:), Rdec(:,:)
        real(dp),    allocatable :: mass(:)
        logical,     allocatable :: view_bad(:), l_live(:)
        integer,     allocatable :: parent(:), remap(:)
        type(string) :: fn
        real(dp) :: r_thresh, rbest, rlink, mad_k, eff_med, eff_mad, eff_cut, rloc
        real(dp) :: rmargin, rqual, thresh_eff, akv
        real     :: mskrad
        integer  :: s, t, u, w, nshell, tbest, nfail, nmerge, npair, kmin_a, kmin_b, aloc, bloc, nnear
        integer  :: nlive
        logical  :: l_gate1, l_single, l_eo, l_fixed_gate
        nstates_out = nstates
        do s = 1, nstates
            label_out(s) = s
        end do
        if( nstates < 2 ) return
        ! zero-mass states (e.g. pruned by AUTO-K upstream) have EMPTY maps on disk; every gate
        ! statistic on them is noise, so they stand aside as singletons and keep their slots
        allocate(mass(nstates), l_live(nstates))
        do s = 1, nstates
            mass(s)   = sum(real(weights(:,s), dp))
            l_live(s) = mass(s) > DTINY
        end do
        nlive = count(l_live)
        if( nlive < nstates ) write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA MERGE: ', nstates - nlive, &
            &' zero-mass states excluded from all merge gates (kept as singleton slots)'
        if( nlive < 2 )then
            deallocate(mass, l_live)
            return
        endif
        if( .not. file_exists('flex_pca_even_state_'//int2str_pad(1,3)//MRC_EXT) .or. &
           &.not. file_exists('flex_pca_odd_state_' //int2str_pad(1,3)//MRC_EXT) )then
            write(logfhandle,'(A)') '>>> FLEX_PCA merge skipped: no half maps on disk'
            call flush(logfhandle)
            return
        endif
        r_thresh = MERGE_R_DEFAULT
        call merge_env_dp('SIMPLE_COV_MERGE_R', r_thresh)
        ! an explicitly supplied gate pins the threshold and disables the adaptive cut below
        l_fixed_gate = .false.
        block
            character(len=32) :: gv
            integer :: gl, gs
            call get_environment_variable('SIMPLE_COV_MERGE_R', gv, gl, gs)
            l_fixed_gate = (gs == 0 .and. gl > 0)
        end block
        ! Decision-robustness levers, both default OFF (= previous behaviour exactly).
        ! MARGIN: merge only when the ratio clears the gate by this much. EO: merge only when BOTH
        ! half-independent cross estimates clear the gate on their own -- a pair the halfsets
        ! disagree about is genuinely ambiguous, and requiring their agreement makes the delivered
        ! K stable under epsilon-level perturbation (thread count, GPU atomic order), which a bare
        ! threshold on their mean is not (measured: one borderline pair flipped 14->13).
        rmargin = 0.d0
        call merge_env_dp('SIMPLE_COV_MERGE_MARGIN', rmargin)
        l_eo = merge_env_is('SIMPLE_COV_MERGE_EO', '1')
        if( rmargin > 0.d0 ) write(logfhandle,'(A,F6.3)') &
            &'>>> FLEX_PCA MERGE gate 2: SIMPLE_COV_MERGE_MARGIN=', rmargin
        if( l_eo ) write(logfhandle,'(A)') &
            &'>>> FLEX_PCA MERGE gate 2: SIMPLE_COV_MERGE_EO=1 -- both halfset estimates must clear the gate'
        ! ---- GATE 1: orientation ----
        ! significant AND an effect-size outlier among this dataset's own states: significance alone
        ! flags everything at these neff, and a bare effect size has no scale that transfers between
        ! specimens
        mad_k = VIEW_MAD_K
        call merge_env_dp('SIMPLE_COV_MERGE_MADK', mad_k)
        allocate(chi2(nstates), neff(nstates), effsz(nstates), view_bad(nstates))
        call view_coverage_chi2(views, weights, nptcls, nstates, chi2, neff, effsz)
        ! robust statistics over LIVE states only, else the zero-mass placeholders drag the median
        allocate(work(nlive))
        t = 0
        do s = 1, nstates
            if( l_live(s) )then
                t = t + 1; work(t) = effsz(s)
            endif
        end do
        call sort_dp(work, nlive)
        eff_med = median_of_sorted(work, nlive)
        t = 0
        do s = 1, nstates
            if( l_live(s) )then
                t = t + 1; work(t) = abs(effsz(s) - eff_med)
            endif
        end do
        call sort_dp(work, nlive)
        eff_mad = median_of_sorted(work, nlive)
        eff_cut = eff_med + mad_k*MAD2SIGMA*eff_mad
        deallocate(work)
        do s = 1, nstates
            view_bad(s) = l_live(s) .and. chi2(s) > VIEW_CHI2_CRIT .and. effsz(s) > eff_cut
        end do
        nfail = count(view_bad)
        write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA MERGE gate 1 (orientation): view-clustered states ', &
            &nfail,' of ',nstates
        ! Premise check: gate 1 assumes conformation is independent of viewing direction, which fails
        ! for a compositional mixture, where distinct species adopt distinct orientation
        ! distributions. A large flagged share is the premise failing rather than the states, so gate
        ! 1 reports and declines to act; gate 2 needs no such assumption.
        if( nfail > nstates/3 )then
            write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA MERGE gate 1 STANDING DOWN: ',nfail,' of ', &
                &nstates,' states depart from the global view distribution, which is the signature of &
                &compositional heterogeneity (species differ in orientation) rather than of view &
                &clustering. Reporting only; gate 2 decides.'
            view_bad = .false.
            nfail    = 0
        endif
        write(logfhandle,'(A,ES12.4,A,ES12.4,A,ES12.4,A,F4.1,A)') '>>>   effect-size median=',eff_med, &
            &'  MAD=',eff_mad,'  cut=',eff_cut,'  (k=',mad_k,' robust sigma)'
        do s = 1, nstates
            write(logfhandle,'(A,I3,A,F12.2,A,F10.1,A,ES12.4,A,L1)') '>>>   state=',s,'  chi2=',chi2(s), &
                &'  neff=',neff(s),'  effect=',effsz(s),'  view_clustered=',view_bad(s)
        end do
        call flush(logfhandle)
        ! ---- GATE 2: volume ----
        allocate(evols(nstates), ovols(nstates))
        mskrad = params%msk_crop
        if( mskrad <= 0. ) mskrad = 0.4*real(params%box_crop)
        ! Reads stay serial: image::new builds this instance's FFTW plans and plan creation is not
        ! thread-safe (only execution is), so read_and_crop cannot run concurrently.
        do s = 1, nstates
            fn = 'flex_pca_even_state_'//int2str_pad(s,3)//MRC_EXT
            call evols(s)%read_and_crop(fn, flex_rec_smpd(params), params%box_crop, params%smpd_crop)
            call fn%kill
            fn = 'flex_pca_odd_state_'//int2str_pad(s,3)//MRC_EXT
            call ovols(s)%read_and_crop(fn, flex_rec_smpd(params), params%box_crop, params%smpd_crop)
            call fn%kill
        end do
        ! Warm the module-level mask coordinate memoization on a throwaway of the same box before
        ! going parallel: mask3D_soft deliberately refuses to memoize inside a parallel region, but
        ! is thread-safe once the memo already matches the box being masked.
        call mskwarm%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
        call mskwarm%mask3D_soft(mskrad)
        call mskwarm%kill
        ! Mask and transform are per-state independent and are where the prologue's time goes:
        ! 2*nstates volume FFTs, which at an over-provisioned ceiling is the bulk of gate 2's setup.
        ! Mask before the transform: solvent reproduces between halves for reasons unrelated to the
        ! state and would enter every spectrum below.
        !$omp parallel do default(shared) private(s) schedule(dynamic) proc_bind(close)
        do s = 1, nstates
            call evols(s)%mask3D_soft(mskrad)
            call ovols(s)%mask3D_soft(mskrad)
        end do
        !$omp end parallel do
        ! ---- DEVIATION FROM THE ENSEMBLE MEAN ----
        ! Every state map is dominated by the density all states SHARE, so an FSC between raw state
        ! maps measures the consensus, not the difference: it drives every pair to the top of the
        ! scale and leaves correct fusions (0.9983-0.9995) crowded against harmful ones
        ! (0.9848-0.9928), which is why no fixed gate transferred between datasets. Subtracting the
        ! ensemble mean first makes the ratio measure what actually distinguishes the states.
        ! Measured 2026-08-16 on deviations: Ribosembly same-GT pairs 0.983-0.999 vs different-GT
        ! 0.413-0.970, and EMPIAR-10076 within-family 0.69-1.03 vs cross-family <=0.878 -- so a
        ! single gate at 0.98 fuses only true duplicates on BOTH, with zero false fusions.
        ! The per-state reliability below is computed from the same deviations, as it must be.
        block
            type(image) :: mean_e, mean_o
            real(dp)    :: wgt
            call mean_e%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
            call mean_o%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
            mean_e = 0.; mean_o = 0.
            do s = 1, nstates
                if( .not. l_live(s) ) cycle
                call mean_e%add(evols(s))
                call mean_o%add(ovols(s))
            end do
            wgt = 1.d0 / real(nlive,dp)
            call mean_e%div(real(nlive))
            call mean_o%div(real(nlive))
            do s = 1, nstates
                if( .not. l_live(s) ) cycle
                call evols(s)%subtr(mean_e)
                call ovols(s)%subtr(mean_o)
            end do
            call mean_e%kill; call mean_o%kill
        end block
        !$omp parallel do default(shared) private(s) schedule(dynamic) proc_bind(close)
        do s = 1, nstates
            call evols(s)%fft
            call ovols(s)%fft
        end do
        !$omp end parallel do
        nshell = evols(1)%get_lfny(1)
        allocate(R(nstates,nstates), Rmin(nstates,nstates), Rdec(nstates,nstates))
        ! optional merge-band cap in Angstrom (SIMPLE_COV_MERGE_LP): conformational
        ! distinctness is a low-frequency question; 25 A validated on both judges
        block
            real(dp) :: mlp
            integer  :: lstop_sh
            mlp = 0.d0
            call merge_env_dp('SIMPLE_COV_MERGE_LP', mlp)
            if( mlp > 0.d0 )then
                lstop_sh = nint(real(params%box_crop,dp)*real(params%smpd_crop,dp)/mlp)
                write(logfhandle,'(A,F6.1,A,I0,A,I0)') '>>> FLEX_PCA MERGE band cap: ', mlp, &
                    &' A -> shells 1..', max(2, min(nshell, lstop_sh)), ' of ', nshell
                call pair_map_ratio(evols, ovols, nstates, nshell, R, Rmin, lstop=lstop_sh)
            else
                call pair_map_ratio(evols, ovols, nstates, nshell, R, Rmin)
            endif
        end block
        ! the DECISION statistic: the mean estimate by default, the conservative halfset minimum
        ! under EO; ordering in the linkage scan always uses the mean (the better point estimate)
        Rdec = R
        if( l_eo ) Rdec = Rmin
        thresh_eff = r_thresh + rmargin
        ! report the ratio distribution even when nothing merges: "no pair reached 0.95" reads the
        ! same whether the closest pair sat at 0.94 or 0.28, and those mean opposite things
        npair = nlive*(nlive-1)/2
        allocate(work(npair))
        npair = 0
        do s = 1, nstates - 1
            if( .not. l_live(s) ) cycle
            do t = s + 1, nstates
                if( .not. l_live(t) ) cycle
                npair       = npair + 1
                work(npair) = R(s,t)
            end do
        end do
        call sort_dp(work, npair)
        thresh_eff = r_thresh + rmargin
        write(logfhandle,'(A,I0,A,F7.4,A,F7.4,A,F7.4,A,F6.3,A)') &
            &'>>> FLEX_PCA MERGE gate 2 (volume): ',npair,' pairs, disattenuated ratio min=',work(1), &
            &'  median=',median_of_sorted(work, npair),'  max=',work(npair),'  (merge at ',r_thresh,')'
        deallocate(work)
        ! Near-gate visibility: any pair within 0.01 of the gate, or whose halfset estimates
        ! straddle it, makes the delivered K sensitive to epsilon-level perturbation. Say so.
        nnear = 0
        do s = 1, nstates - 1
            if( .not. l_live(s) ) cycle
            do t = s + 1, nstates
                if( .not. l_live(t) ) cycle
                if( abs(R(s,t) - r_thresh) <= 0.01d0 .or. &
                   &(R(s,t) >= r_thresh .and. Rmin(s,t) < r_thresh) )then
                    nnear = nnear + 1
                    if( nnear <= 5 ) write(logfhandle,'(A,I3,A,I3,A,F7.4,A,F7.4,A)') &
                        &'>>> FLEX_PCA MERGE gate 2 NEAR-GATE: pair ',s,',',t,'  ratio=',R(s,t), &
                        &'  halfset min=',Rmin(s,t),' -- this merge decision is perturbation-sensitive'
                endif
            end do
        end do
        if( nnear > 5 ) write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA MERGE gate 2 NEAR-GATE: ', &
            &nnear-5, ' more pairs suppressed'
        call flush(logfhandle)
        ! ---- AGGLOMERATE ----
        ! COMPLETE linkage: a group merges only if EVERY cross pair clears the threshold — the
        ! transitive closure this replaces let one borderline pair chain unlike groups.
        ! SIMPLE_COV_MERGE_LINK=single restores the old behaviour.
        allocate(parent(nstates))
        do s = 1, nstates
            parent(s) = s                              ! parent doubles as the cluster id
        end do
        nmerge   = 0
        l_single = merge_env_is('SIMPLE_COV_MERGE_LINK', 'single')
        if( l_single )then
            write(logfhandle,'(A)') '>>> FLEX_PCA MERGE gate 2: SINGLE linkage (transitive closure)'
            do s = 1, nstates - 1
                if( .not. l_live(s) ) cycle
                do t = s + 1, nstates
                    if( .not. l_live(t) ) cycle
                    if( Rdec(s,t) >= thresh_eff )then
                        call uf_union(parent, s, t, nmerge)
                        write(logfhandle,'(A,I3,A,I3,A,F7.4,A)') '>>> FLEX_PCA MERGE gate 2: states ',s,' + ',t, &
                            &'  disattenuated ratio=',R(s,t),' -> indistinguishable within their own noise'
                    endif
                end do
            end do
        else
            ! greedily merge the qualifying cluster pair with the STRONGEST weakest link, repeat.
            ! The candidate scan is O(nstates^2) cluster pairs x O(|a|*|b|) cross pairs per merge and
            ! runs once per merge, so it is threaded over the outer representative. Each thread keeps
            ! its own best and the winners are combined under a critical section; both the local and
            ! the combining test break ties on the LOWEST (s,t), so the pair chosen is independent of
            ! thread count and schedule and matches the serial scan exactly.
            do
                rbest = -1.d0; kmin_a = 0; kmin_b = 0
                !$omp parallel default(shared) private(s,t,u,w,rlink,rqual,rloc,aloc,bloc) proc_bind(close)
                rloc = -1.d0; aloc = 0; bloc = 0
                !$omp do schedule(dynamic)
                do s = 1, nstates - 1
                    if( parent(s) /= s .or. .not. l_live(s) ) cycle   ! only LIVE cluster representatives
                    do t = s + 1, nstates
                        if( parent(t) /= t .or. .not. l_live(t) ) cycle
                        ! weakest cross pair between the two clusters: rlink (mean estimate) orders
                        ! candidates, rqual (decision statistic) qualifies them
                        rlink = huge(1.d0)
                        rqual = huge(1.d0)
                        do u = 1, nstates
                            if( parent(u) /= s ) cycle
                            do w = 1, nstates
                                if( parent(w) /= t ) cycle
                                rlink = min(rlink, R(u,w))
                                rqual = min(rqual, Rdec(u,w))
                            end do
                        end do
                        if( rqual >= thresh_eff )then
                            if( aloc == 0 )then
                                rloc = rlink; aloc = s; bloc = t
                            else if( rlink > rloc )then
                                rloc = rlink; aloc = s; bloc = t
                            else if( rlink == rloc .and. (s < aloc .or. (s == aloc .and. t < bloc)) )then
                                rloc = rlink; aloc = s; bloc = t
                            endif
                        endif
                    end do
                end do
                !$omp end do
                !$omp critical (flex_pca_merge_best)
                if( aloc > 0 )then
                    if( kmin_a == 0 )then
                        rbest = rloc; kmin_a = aloc; kmin_b = bloc
                    else if( rloc > rbest )then
                        rbest = rloc; kmin_a = aloc; kmin_b = bloc
                    else if( rloc == rbest .and. (aloc < kmin_a .or. &
                        &(aloc == kmin_a .and. bloc < kmin_b)) )then
                        rbest = rloc; kmin_a = aloc; kmin_b = bloc
                    endif
                endif
                !$omp end critical (flex_pca_merge_best)
                !$omp end parallel
                if( kmin_a < 1 ) exit
                do u = 1, nstates                      ! fold b's members into a
                    if( parent(u) == kmin_b ) parent(u) = kmin_a
                end do
                nmerge = nmerge + 1
                write(logfhandle,'(A,I3,A,I3,A,F7.4,A)') '>>> FLEX_PCA MERGE gate 2: states ',kmin_a, &
                    &' + ',kmin_b,'  weakest cross ratio=',rbest,' -> indistinguishable within their own noise'
            end do
        endif
        ! Under AUTO-K the merge serves as a MAP-SPACE DUPLICATE FUSE only: gate 2 decides, gate 1
        ! reports but never folds. Measured 2026-08-16 on Ribosembly: gate 1 folded 3 genuine
        ! compositional states (map ratios 0.86-0.92, well below the duplicate gate) because its
        ! premise -- conformation independent of viewing direction -- fails for a compositional
        ! mixture, and the >1/3 stand-down did not trip at 4 of 16 flagged. AUTO-K's own
        ! reproducibility gate already handles unstable states.
        ! any non-zero value, matching how the model module gates AUTO-K (cov_env_int_pub > 0)
        akv = 0.d0
        call merge_env_dp('SIMPLE_COV_AUTO_K', akv)
        if( akv > 0.d0 )then
            if( count(view_bad) > 0 ) write(logfhandle,'(A,I0,A)') &
                &'>>> FLEX_PCA MERGE gate 1: AUTO-K active -- ', count(view_bad), &
                &' view-clustered states REPORTED only, not folded'
            view_bad = .false.
        endif
        ! fold each view-contaminated state into the state its map most resembles, rather than
        ! deleting it, so its particles keep contributing somewhere
        do s = 1, nstates
            if( .not. view_bad(s) ) cycle
            rbest = -1.d0; tbest = 0
            do t = 1, nstates
                if( t == s ) cycle
                if( .not. l_live(t) ) cycle            ! never fold into a zero-mass placeholder
                if( view_bad(t) ) cycle                ! do not fold one bad state into another
                if( R(s,t) > rbest )then
                    rbest = R(s,t); tbest = t
                endif
            end do
            if( tbest < 1 .or. rbest < VIEW_FOLD_MIN_R )then
                write(logfhandle,'(A,I3,A,F7.4,A,F4.2,A)') '>>> FLEX_PCA MERGE gate 1: state ',s, &
                    &' is view-clustered but its best map match is only ',max(rbest,0.d0), &
                    &' (floor ',VIEW_FOLD_MIN_R,') -- KEPT, inspect it: view-driven and unlike every other state'
                cycle
            endif
            l_gate1 = uf_find(parent, s) /= uf_find(parent, tbest)
            if( l_gate1 )then
                call uf_union(parent, s, tbest, nmerge)
                write(logfhandle,'(A,I3,A,I3,A,F7.4)') '>>> FLEX_PCA MERGE gate 1: state ',s, &
                    &' is view-clustered, folded into ',tbest,'  ratio=',rbest
            endif
        end do
        ! relabel the surviving roots 1..nstates_out, preserving input order
        allocate(remap(nstates), source=0)
        nstates_out = 0
        do s = 1, nstates
            if( uf_find(parent, s) == s )then
                nstates_out = nstates_out + 1
                remap(s)    = nstates_out
            endif
        end do
        do s = 1, nstates
            label_out(s) = remap(uf_find(parent, s))
        end do
        write(logfhandle,'(A,I0,A,I0,A,I0,A)') '>>> FLEX_PCA MERGE result: ',nstates,' -> ',nstates_out, &
            &' states (',nmerge,' merges)'
        call flush(logfhandle)
        do s = 1, nstates
            call evols(s)%kill; call ovols(s)%kill
        end do
        deallocate(evols, ovols, chi2, neff, effsz, view_bad, R, parent, remap, mass, l_live)
    end subroutine two_gate_state_merge

    !> Ascending insertion sort. n is the state count, so the O(n^2) is irrelevant and this avoids
    !! converting to the single precision the shared hpsort generics cover.
    pure subroutine sort_dp( x, n )
        integer,  intent(in)    :: n
        real(dp), intent(inout) :: x(n)
        real(dp) :: key
        integer  :: i, j
        do i = 2, n
            key = x(i)
            j   = i - 1
            do while( j >= 1 )
                if( x(j) <= key ) exit
                x(j+1) = x(j)
                j      = j - 1
            end do
            x(j+1) = key
        end do
    end subroutine sort_dp

    !> Median of an ASCENDING array; even n averages the two central values.
    pure real(dp) function median_of_sorted( x, n ) result( med )
        integer,  intent(in) :: n
        real(dp), intent(in) :: x(n)
        if( n < 1 )then
            med = 0.d0
        else if( mod(n,2) == 1 )then
            med = x((n+1)/2)
        else
            med = 0.5d0*(x(n/2) + x(n/2+1))
        endif
    end function median_of_sorted

    recursive integer function uf_find( parent, i ) result( root )
        integer, intent(inout) :: parent(:)
        integer, intent(in)    :: i
        if( parent(i) == i )then
            root = i
        else
            root      = uf_find(parent, parent(i))
            parent(i) = root         ! path compression
        endif
    end function uf_find


    subroutine uf_union( parent, i, j, nmerge )
        integer, intent(inout) :: parent(:)
        integer, intent(in)    :: i, j
        integer, intent(inout) :: nmerge
        integer :: ri, rj
        ri = uf_find(parent, i)
        rj = uf_find(parent, j)
        if( ri == rj ) return
        parent(max(ri,rj)) = min(ri,rj)   ! keep the lower index as the root, so labels stay ordered
        nmerge = nmerge + 1
    end subroutine uf_union

end module simple_flex_pca_merge
