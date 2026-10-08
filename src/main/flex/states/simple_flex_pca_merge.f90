!@descr: Two-gate agglomerative merge of over-provisioned flex_pca states.
!        Gate 1 (orientation): a state whose viewing-axis distribution stands out from its peers is a view
!        cluster, folded into a sufficiently similar map. Gate 2 (volume): pairs whose deviation maps agree within
!        their own half-map reproducibility fuse under complete linkage (simple_hac on 1 - ratio). Latent
!        distance never merges.
!        Invoked by the application when preimage_auto=yes.
module simple_flex_pca_merge
use simple_core_module_api, only: dp, dtiny, file_exists, int2str_pad, logfhandle, mrc_ext, simple_exception, string
use simple_defs_flex,        only: FLEX_FSC_SIGNAL_THRESHOLD
use simple_flex_pca_records, only: flex_state_set
use simple_image,            only: image
use simple_parameters,       only: parameters
use simple_flex_pca_pcg,     only: flex_pcg_environment
use simple_stat,             only: median, mad
use simple_hac,              only: hac
implicit none
private
#include "simple_local_flags.inc"

public :: two_gate_state_merge

!> Viewing-axis second moments are symmetric 3x3 with unit trace (v is a unit vector), leaving
!! five free parameters.
integer,  parameter :: VIEW_DOF        = 5
!> Upper tail of chi^2_5 at p = 0.001. Necessary but not sufficient: at realistic neff the test is
!! over-powered and every state clears it, since no experimental view distribution is exactly the
!! global one. It only guards the sparse states, where a departure really could be sampling.
real(dp), parameter :: VIEW_CHI2_CRIT  = 20.515d0
!> The decision is therefore made on effect size, chi2/neff, which does not grow with N. The cut is
!! robust and relative to this dataset's own states (median + k*MAD): a specimen with preferred
!! orientation gives every state a non-global view distribution, so what marks a view cluster is
!! standing out from its peers. Gate 2 covers the case where every state is a view cluster and
!! there are no outliers to find.
real(dp), parameter :: VIEW_MAD_K      = 3.0d0
!> MAD -> sigma for a normal, so VIEW_MAD_K is read in familiar units.
real(dp), parameter :: MAD2SIGMA       = 1.4826d0
!> A view-clustered state is folded into its best map match only if it substantially resembles it.
!! Without the floor gate 1 folds on the best available match however poor it is, tipping particles
!! into maps they do not resemble. Below the floor the state is reported and kept: a component that
!! is view-driven and resembles nothing is a finding, not something to hide in an unrelated class.
real(dp), parameter :: VIEW_FOLD_MIN_R = 0.8d0
!> Default gate on the DEVIATION-from-ensemble-mean disattenuated ratio (see pair_map_ratio); 1 means
!! the two maps agree as well as each agrees with itself.
real(dp), parameter :: MERGE_R_DEFAULT = 0.98d0


contains


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
    subroutine pair_map_ratio( evols, ovols, nstates, nshell, R, Rmin )
        integer,     intent(in)    :: nstates, nshell
        type(image), intent(inout) :: evols(nstates), ovols(nstates)
        real(dp),    intent(out)   :: R(nstates,nstates)
        !> min of the TWO half-independent cross estimates (e_s vs o_t and o_s vs e_t), each
        !! disattenuated by the same reliabilities. R is their mean, so 2*(R - Rmin) is the
        !! halfset disagreement — a per-pair noise scale for the merge decision.
        real(dp),    intent(out)   :: Rmin(nstates,nstates)
        real,     allocatable :: css(:,:), cst(:), cts(:)
        real(dp) :: num_a, num_b, den, ratio, ratio_a, ratio_b
        integer  :: s, t, l, nval
        allocate(css(nshell,nstates), source=0.)
        allocate(cst(nshell), cts(nshell), source=0.)
        do s = 1, nstates
            call evols(s)%fsc(ovols(s), css(:,s))
        end do
        R = 1.d0
        Rmin = 1.d0
        do s = 1, nstates - 1
            do t = s + 1, nstates
                call evols(s)%fsc(ovols(t), cst)
                call ovols(s)%fsc(evols(t), cts)
                num_a = 0.d0
                num_b = 0.d0
                den   = 0.d0
                nval  = 0
                do l = 1, nshell
                    ! a shell where either state is already noise carries no evidence
                    if( real(css(l,s),dp) <= FLEX_FSC_SIGNAL_THRESHOLD ) cycle
                    if( real(css(l,t),dp) <= FLEX_FSC_SIGNAL_THRESHOLD ) cycle
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
                Rmin(s,t) = min(ratio_a, ratio_b)
                Rmin(t,s) = Rmin(s,t)
            end do
        end do
        deallocate(css, cst, cts)
    end subroutine pair_map_ratio

    !> Two-gate merge. Returns the state each input state was merged into (1..nstates_out).
    subroutine two_gate_state_merge( params, env, views, states, half_prefix, label_out, nstates_out )
        class(flex_pcg_environment), intent(inout) :: env
        type(flex_state_set),        intent(in)    :: states
        integer :: nptcls
        class(parameters),    intent(in)  :: params
        real(dp),             intent(in)  :: views(:,:)
        character(len=*),     intent(in)  :: half_prefix  !< raw state half maps <half_prefix>_stateNN_{even,odd}.mrc
        integer,              intent(out) :: label_out(:)
        integer,              intent(out) :: nstates_out
        type(image), allocatable :: evols(:), ovols(:)
        type(image) :: mskwarm
        real(dp),    allocatable :: chi2(:), neff(:), effsz(:), work(:), R(:,:), Rmin(:,:)
        real(dp),    allocatable :: mass(:), dmerge(:,:), heights(:)
        logical,     allocatable :: view_bad(:), l_live(:)
        integer,     allocatable :: lab(:), remap(:), medoids(:), pairs(:,:)
        type(hac)    :: hc
        type(string) :: fn
        real(dp) :: rbest, eff_med, eff_mad, eff_cut
        real     :: mskrad
        integer  :: s, t, m, nshell, tbest, nfail, nmerge, npair, nnear, lold, lnew
        integer  :: nlive
        nptcls = size(states%weights,1)
        nstates_out = states%nstates
        do s = 1, states%nstates
            label_out(s) = s
        end do
        if( states%nstates < 2 ) return
        ! zero-mass states have EMPTY maps on disk; every gate
        ! statistic on them is noise, so they stand aside as singletons and keep their slots
        allocate(mass(states%nstates), l_live(states%nstates))
        do s = 1, states%nstates
            mass(s)   = sum(real(states%weights(:,s), dp))
            l_live(s) = mass(s) > DTINY
        end do
        nlive = count(l_live)
        if( nlive < states%nstates ) write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA MERGE: ', states%nstates - nlive, &
            &' zero-mass states excluded from all merge gates (kept as singleton slots)'
        if( nlive < 2 )then
            deallocate(mass, l_live)
            return
        endif
        if( .not. file_exists(half_prefix//'_state'//int2str_pad(1,2)//'_even'//MRC_EXT) .or. &
           &.not. file_exists(half_prefix//'_state'//int2str_pad(1,2)//'_odd' //MRC_EXT) )then
            write(logfhandle,'(A)') '>>> FLEX_PCA merge skipped: no half maps on disk'
            call flush(logfhandle)
            return
        endif
        ! ---- GATE 1: orientation ----
        ! significant AND an effect-size outlier among this dataset's own states: significance alone
        ! flags everything at these neff, and a bare effect size has no scale that transfers between
        ! specimens
        allocate(chi2(states%nstates), neff(states%nstates), effsz(states%nstates), view_bad(states%nstates))
        call view_coverage_chi2(views, states%weights, nptcls, states%nstates, chi2, neff, effsz)
        ! robust statistics over LIVE states only, else the zero-mass placeholders drag the median
        allocate(work(nlive))
        t = 0
        do s = 1, states%nstates
            if( l_live(s) )then
                t = t + 1; work(t) = effsz(s)
            endif
        end do
        eff_med = median(work)
        eff_mad = mad(work, eff_med)
        eff_cut = eff_med + VIEW_MAD_K*MAD2SIGMA*eff_mad
        deallocate(work)
        do s = 1, states%nstates
            view_bad(s) = l_live(s) .and. chi2(s) > VIEW_CHI2_CRIT .and. effsz(s) > eff_cut
        end do
        nfail = count(view_bad)
        write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA MERGE gate 1 (orientation): view-clustered states ', &
            &nfail,' of ',states%nstates
        ! Premise check: gate 1 assumes conformation is independent of viewing direction, which fails
        ! for a compositional mixture, where distinct species adopt distinct orientation
        ! distributions. A large flagged share is the premise failing rather than the states, so gate
        ! 1 reports and declines to act; gate 2 needs no such assumption.
        if( nfail > states%nstates/3 )then
            write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA MERGE gate 1 STANDING DOWN: ',nfail,' of ', &
                &states%nstates,' states depart from the global view distribution, which is the signature of &
                &compositional heterogeneity (species differ in orientation) rather than of view &
                &clustering. Reporting only; gate 2 decides.'
            view_bad = .false.
            nfail    = 0
        endif
        write(logfhandle,'(A,ES12.4,A,ES12.4,A,ES12.4,A,F4.1,A)') '>>>   effect-size median=',eff_med, &
            &'  MAD=',eff_mad,'  cut=',eff_cut,'  (k=',VIEW_MAD_K,' robust sigma)'
        do s = 1, states%nstates
            write(logfhandle,'(A,I3,A,F12.2,A,F10.1,A,ES12.4,A,L1)') '>>>   state=',s,'  chi2=',chi2(s), &
                &'  neff=',neff(s),'  effect=',effsz(s),'  view_clustered=',view_bad(s)
        end do
        call flush(logfhandle)
        ! ---- GATE 2: volume ----
        allocate(evols(states%nstates), ovols(states%nstates))
        mskrad = params%msk_crop
        if( mskrad <= 0. ) mskrad = 0.4*real(params%box_crop)
        ! Reads stay serial: image::new builds this instance's FFTW plans and plan creation is not
        ! thread-safe (only execution is), so read_and_crop cannot run concurrently.
        do s = 1, states%nstates
            fn = half_prefix//'_state'//int2str_pad(s,2)//'_even'//MRC_EXT
            call evols(s)%read_and_crop(fn, params%smpd_crop, params%box_crop, params%smpd_crop)
            call fn%kill
            fn = half_prefix//'_state'//int2str_pad(s,2)//'_odd'//MRC_EXT
            call ovols(s)%read_and_crop(fn, params%smpd_crop, params%box_crop, params%smpd_crop)
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
        do s = 1, states%nstates
            if( env%active() )then
                call env%apply(evols(s), params%box_crop, params%msk_crop)
                call env%apply(ovols(s), params%box_crop, params%msk_crop)
            else
                call evols(s)%mask3D_soft(mskrad)
                call ovols(s)%mask3D_soft(mskrad)
            endif
        end do
        !$omp end parallel do
        ! ---- DEVIATION FROM THE ENSEMBLE MEAN ----
        ! Raw state maps are dominated by the shared density, so their FSC measures the consensus; the ratio
        ! and the per-state reliability are both computed on deviations from the live-state mean.
        block
            type(image) :: mean_e, mean_o
            real(dp)    :: wgt
            call mean_e%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
            call mean_o%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
            mean_e = 0.; mean_o = 0.
            do s = 1, states%nstates
                if( .not. l_live(s) ) cycle
                call mean_e%add(evols(s))
                call mean_o%add(ovols(s))
            end do
            wgt = 1.d0 / real(nlive,dp)
            call mean_e%div(real(nlive))
            call mean_o%div(real(nlive))
            do s = 1, states%nstates
                if( .not. l_live(s) ) cycle
                call evols(s)%subtr(mean_e)
                call ovols(s)%subtr(mean_o)
            end do
            call mean_e%kill; call mean_o%kill
        end block
        !$omp parallel do default(shared) private(s) schedule(dynamic) proc_bind(close)
        do s = 1, states%nstates
            call evols(s)%fft
            call ovols(s)%fft
        end do
        !$omp end parallel do
        nshell = evols(1)%get_lfny(1)
        allocate(R(states%nstates,states%nstates), Rmin(states%nstates,states%nstates))
        call pair_map_ratio(evols, ovols, states%nstates, nshell, R, Rmin)
        ! The halfset mean decides the merge; the minimum is retained for near-gate diagnostics.
        ! report the ratio distribution even when nothing merges: "no pair reached 0.95" reads the
        ! same whether the closest pair sat at 0.94 or 0.28, and those mean opposite things
        npair = nlive*(nlive-1)/2
        allocate(work(npair))
        npair = 0
        do s = 1, states%nstates - 1
            if( .not. l_live(s) ) cycle
            do t = s + 1, states%nstates
                if( .not. l_live(t) ) cycle
                npair       = npair + 1
                work(npair) = R(s,t)
            end do
        end do
        write(logfhandle,'(A,I0,A,F7.4,A,F7.4,A,F7.4,A,F6.3,A)') &
            &'>>> FLEX_PCA MERGE gate 2 (volume): ',npair,' pairs, disattenuated ratio min=',minval(work), &
            &'  median=',median(work),'  max=',maxval(work),'  (merge at ',MERGE_R_DEFAULT,')'
        deallocate(work)
        ! Near-gate visibility: any pair within 0.01 of the gate, or whose halfset estimates
        ! straddle it, makes the delivered K sensitive to epsilon-level perturbation. Say so.
        nnear = 0
        do s = 1, states%nstates - 1
            if( .not. l_live(s) ) cycle
            do t = s + 1, states%nstates
                if( .not. l_live(t) ) cycle
                if( abs(R(s,t) - MERGE_R_DEFAULT) <= 0.01d0 .or. &
                   &(R(s,t) >= MERGE_R_DEFAULT .and. Rmin(s,t) < MERGE_R_DEFAULT) )then
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
        ! COMPLETE linkage on 1 - R: a group merges only if EVERY cross pair clears the gate (the transitive
        ! closure this replaced let one borderline pair chain unlike groups); the closest qualifying pair,
        ! the one with the strongest weakest link, merges first. Zero-mass states stay singletons.
        allocate(dmerge(states%nstates,states%nstates))
        dmerge = 1.d0 - R
        do s = 1, states%nstates
            dmerge(s,s) = 0.d0
        end do
        call hc%new(states%nstates, dmerge, 'complete', thres=1.d0 - MERGE_R_DEFAULT, mask=l_live)
        call hc%cluster(medoids, lab)
        call hc%get_history(pairs, heights)
        call hc%kill
        nmerge = size(heights)
        do m = 1, nmerge
            write(logfhandle,'(A,I3,A,I3,A,F7.4,A)') '>>> FLEX_PCA MERGE gate 2: states ',pairs(1,m), &
                &' + ',pairs(2,m),'  weakest cross ratio=',1.d0 - heights(m),' -> indistinguishable within their own noise'
        end do
        ! fold each view-contaminated state into the state its map most resembles, rather than
        ! deleting it, so its particles keep contributing somewhere
        do s = 1, states%nstates
            if( .not. view_bad(s) ) cycle
            rbest = -1.d0; tbest = 0
            do t = 1, states%nstates
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
            if( lab(s) /= lab(tbest) )then
                lold = max(lab(s), lab(tbest))
                lnew = min(lab(s), lab(tbest))
                where( lab == lold ) lab = lnew
                nmerge = nmerge + 1
                write(logfhandle,'(A,I3,A,I3,A,F7.4)') '>>> FLEX_PCA MERGE gate 1: state ',s, &
                    &' is view-clustered, folded into ',tbest,'  ratio=',rbest
            endif
        end do
        ! number the surviving groups 1..nstates_out in order of their smallest member, preserving input order
        allocate(remap(states%nstates), source=0)
        nstates_out = 0
        do s = 1, states%nstates
            if( remap(lab(s)) == 0 )then
                nstates_out    = nstates_out + 1
                remap(lab(s))  = nstates_out
            endif
            label_out(s) = remap(lab(s))
        end do
        write(logfhandle,'(A,I0,A,I0,A,I0,A)') '>>> FLEX_PCA MERGE result: ',states%nstates,' -> ',nstates_out, &
            &' states (',nmerge,' merges)'
        call flush(logfhandle)
        do s = 1, states%nstates
            call evols(s)%kill; call ovols(s)%kill
        end do
        deallocate(evols, ovols, chi2, neff, effsz, view_bad, R, Rmin, dmerge, lab, remap, medoids, pairs, heights, &
            &mass, l_live)
    end subroutine two_gate_state_merge

end module simple_flex_pca_merge
