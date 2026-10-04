!@descr: gridding vs PCG reconstruction comparison of simple_commanders_test_highlevel
submodule (simple_commanders_test_highlevel) simple_commanders_test_highlevel_backends
implicit none
#include "simple_local_flags.inc"

!> Collated per-run observables of one rec3D_backends comparison
type backends_run_summary
    real    :: g_base05    =  0.   !< gridding base half-pair FSC=0.5 (A)
    real    :: g_base0143  =  0.   !< gridding base half-pair FSC=0.143 (A)
    real    :: p_base05    =  0.   !< pcg base (_unfil) half-pair FSC=0.5 -- the negative control
    real    :: p_base0143  =  0.   !< pcg base (_unfil) half-pair FSC=0.143
    real    :: p_ship05    =  0.   !< pcg shipped (regularized) half-pair FSC=0.5 -- inflation diagnostic
    real    :: p_ship0143  =  0.   !< pcg shipped (regularized) half-pair FSC=0.143
    real    :: band_ratio  = -1.   !< median gated-band amplitude ratio pcg/gridding
    real    :: band_fsc    = -1.   !< median gated-band FSC(gridding,pcg) -- backend agreement
    real    :: rad_min     = -1.   !< min normalised in-mask radial ratio -- erosion indicator
    real    :: rad_max     = -1.   !< max normalised in-mask radial ratio
    real    :: centre_bin  = -1.   !< pcg centre-bin ratio (diagnostic)
    real    :: truth_fsc_g = -1.   !< median gated-band FSC(truth,gridding), when vol1 given
    real    :: truth_fsc_p = -1.   !< median gated-band FSC(truth,pcg), when vol1 given
    integer :: nfail       =  0    !< gate violations in the run
end type backends_run_summary

contains

!> Gridding vs PCG reconstruct3D on the same project/poses/sigma2 (numbered exec dir unless mkdir=no).
!> Hard gates: agreement band, in-band amplitude ratio and FSC, radial-ratio range; with truth vol1
!> also truth FSC (LS flatness only if ml_reg=no). Thresholds: doc/implementation_notes/completed/drop_legacy_box_division.md
module subroutine exec_test_rec3D_backends( self, cline )
    class(commander_test_rec3D_backends), intent(inout) :: self
    class(cmdline),                       intent(inout) :: cline
    type(backends_run_summary) :: summary
    logical :: l_mlreg
    ! prior-capable defaults: a bare invocation (projfile, pgrp, mskdiam,
    ! nthr) runs a single euclid+ml_reg comparison at the
    ! production-representative budget; every default is overridable
    if( .not. cline%defined('objfun')     ) call cline%set('objfun',     'euclid')
    if( .not. cline%defined('ml_reg')     ) call cline%set('ml_reg',        'yes')
    if( .not. cline%defined('maxits_pcg') ) call cline%set('maxits_pcg',       5.)
    if( .not. cline%defined('rtol')       ) call cline%set('rtol',         1.e-3)
    l_mlreg = cline%get_carg('ml_reg') .eq. 'yes'
    call run_rec3D_backends_single(cline, summary, .true.)
    call simple_end('**** SIMPLE_TEST_REC3D_BACKENDS NORMAL STOP ****', print_simple=.false.)
end subroutine exec_test_rec3D_backends

subroutine run_rec3D_backends_single( cline, summary, l_abort_on_fail )
    use simple_commanders_rec,  only: commander_rec3D
    use simple_image,           only: image
    use simple_sp_project,      only: sp_project
    use simple_refine3D_fnames, only: refine3D_state_vol_fname, refine3D_state_vol_fbody, &
        &refine3D_state_halfvol_fname
    class(cmdline),              intent(inout) :: cline
    type(backends_run_summary),  intent(out)   :: summary
    logical,                     intent(in)    :: l_abort_on_fail
    character(len=8), parameter :: BACKENDS(2) = [character(len=8) :: 'gridding', 'pcg']
    integer,          parameter :: NRBINS = 16
    type(commander_rec3D) :: xrec3D
    type(cmdline)         :: cline_rec
    type(sp_project)      :: spproj
    type(image)           :: vols(2)
    type(string)          :: projfile, vol_fname, out_fnames(2)
    real,    allocatable  :: spec(:,:), spec_tmp(:), corrs(:), radprof(:,:), ratios(:)
    real,    allocatable  :: tprof(:,:), fscs(:)
    integer               :: nfail, nfsc, kgate
    real                  :: med_fsc, tg
    integer, allocatable  :: radcnt(:)
    real,    pointer      :: rmat(:,:,:) => null()
    type(image)           :: truth, tvols(2)
    type(string)          :: truth_fname, exec_dir, dirbody
    real    :: smpd, smpd_out, mskrad, rbin_width, r, l2(2), rr, rnorm, rmin, rmax, med, lp_here, bg, hp_here
    real    :: r05_tmp, r0143_tmp
    real,    allocatable  :: tspec(:), tcorr(:,:), tspec_b(:,:)
    type(string)          :: cwd_orig
    integer :: ldim(3), ib, state, nptcls, k, lfny, irb, i, j, l, c(3), nrb_msk, n, box, kagree, nrb_used, it
    logical :: l_truth, l_gate_ls, l_mkdir
    if( .not. cline%defined('trs')     ) call cline%set('trs', 5.)
    if( .not. cline%defined('mskdiam') ) THROW_HARD('mskdiam is required; exec_test_rec3D_backends')
    l_mkdir = .true.
    if( cline%defined('mkdir') ) l_mkdir = cline%get_carg('mkdir') .ne. 'no'
    call cline%set('oritype',     'ptcl3D')
    call cline%set('mkdir',       'no')   ! the children run in THIS process's cwd
    call cline%set('postprocess', 'no')
    call cline%delete('nparts')   ! shared-memory execution of both backends
    call cline%delete('part')
    call cline%delete('rec_backend')
    if( l_mkdir )then
        ! numbered execution directory like production commanders, with the
        ! settings that define the measurement in the name so a sweep reads as
        ! a directory listing. The project file is addressed absolutely so its
        ! updates land in the caller's project, as with any ../-style execution
        ! directory; the euclid sigmas come from the project's registered
        ! canonical state, so nothing is discovered by name in the cwd.
        dirbody = string('rec3D_backends')
        if( cline%defined('pgrp')       ) dirbody = dirbody//('_'//cline_tok('pgrp'))
        if( cline%defined('objfun')     ) dirbody = dirbody//('_'//cline_tok('objfun'))
        if( cline%defined('ml_reg') )then
            if( cline%get_carg('ml_reg') .eq. 'yes' ) dirbody = dirbody//'_mlreg'
        endif
        if( cline%defined('maxits_pcg') ) dirbody = dirbody//('_its'//int2str(cline%get_iarg('maxits_pcg')))
        if( cline%defined('rtol')       ) dirbody = dirbody//('_rtol'//real_tok(cline%get_rarg('rtol')))
        if( cline%defined('lp')         ) dirbody = dirbody//('_lp'//real_tok(cline%get_rarg('lp')))
        do i = 1, 9999
            exec_dir = string(int2str(i)//'_')//dirbody
            if( .not. dir_exists(exec_dir) ) exit
        end do
        call simple_mkdir(exec_dir)
        call cline%set('projfile', simple_abspath(cline%get_carg('projfile')))
        if( cline%defined('vol1') ) call cline%set('vol1', simple_abspath(cline%get_carg('vol1')))
        call simple_getcwd(cwd_orig)   ! the runner returns here after the comparison
        call simple_chdir(exec_dir)
        write(logfhandle,'(a)') '>>> REC3D BACKENDS: EXECUTION DIRECTORY '//exec_dir%to_char()
    endif
    projfile = cline%get_carg('projfile')
    call spproj%read(projfile)
    smpd = spproj%get_smpd()
    box  = spproj%get_box()
    call spproj%kill
    state = 1
    if( cline%defined('state') ) state = cline%get_iarg('state')
    ! reconstruct with both backends
    do ib = 1, 2
        cline_rec = cline
        call cline_rec%set('prg',         'reconstruct3D')
        call cline_rec%set('rec_backend', trim(BACKENDS(ib)))
        call cline_rec%delete('vol1')   ! ground-truth volume is for the comparison only
        call cline_rec%delete('lp')
        call cline_rec%delete('hp')
        write(logfhandle,'(A)') '>>> REC3D BACKENDS: RECONSTRUCTING WITH '//trim(BACKENDS(ib))
        call xrec3D%execute(cline_rec)
        vol_fname = refine3D_state_vol_fname(state)
        if( .not. file_exists(vol_fname) )then
            THROW_HARD('reconstruct3D ('//trim(BACKENDS(ib))//') did not produce '//vol_fname%to_char())
        endif
        ! Diagnostic half-pair FSCs, per backend. The base (_unfil) pair
        ! carries the workflow's resolution meaning; the shipped pair is
        ! regularized (ML, and any priors on the pcg leg) with the SAME
        ! regularizers on both halves, so its FSC is inflated by shared
        ! regularization -- reported to SEE that effect, never as a
        ! resolution claim.
        if( cline%defined('ml_reg') )then
            if( cline%get_carg('ml_reg') .eq. 'yes' )then
                call report_pair_fsc(BACKENDS(ib), .true.,  'base (_unfil)  ', r05_tmp, r0143_tmp)
                if( ib == 1 )then
                    summary%g_base05 = r05_tmp; summary%g_base0143 = r0143_tmp
                else
                    summary%p_base05 = r05_tmp; summary%p_base0143 = r0143_tmp
                endif
                call report_pair_fsc(BACKENDS(ib), .false., 'shipped (regul)', r05_tmp, r0143_tmp)
                if( ib == 2 )then
                    summary%p_ship05 = r05_tmp; summary%p_ship0143 = r0143_tmp
                    write(logfhandle,'(a)') &
                        &'    (shipped-pair FSC shares regularization between halves; diagnostic only)'
                endif
            endif
        endif
        out_fnames(ib) = refine3D_state_vol_fbody(state)//'_'//trim(BACKENDS(ib))//MRC_EXT
        call simple_rename(vol_fname, out_fnames(ib), overwrite=.true.)
        call cline_rec%kill
    enddo
    ! read the two maps
    call find_ldim_nptcls(out_fnames(1), ldim, nptcls)
    smpd_out = smpd * real(box) / real(ldim(1))
    do ib = 1, 2
        call vols(ib)%new(ldim, smpd_out)
        call vols(ib)%read(out_fnames(ib))
    enddo
    ! real-space radial profiles of |rho| (before any FFT)
    c          = ldim/2 + 1
    rbin_width = real(ldim(1)/2) / real(NRBINS)
    allocate(radprof(NRBINS,2), source=0.)
    allocate(radcnt(NRBINS),    source=0)
    do ib = 1, 2
        call vols(ib)%get_rmat_ptr(rmat)
        l2(ib) = sqrt(sum(rmat(1:ldim(1),1:ldim(2),1:ldim(3))**2))
        do l = 1, ldim(3)
            do j = 1, ldim(2)
                do i = 1, ldim(1)
                    r   = sqrt(real((i-c(1))**2 + (j-c(2))**2 + (l-c(3))**2))
                    irb = int(r / rbin_width) + 1
                    if( irb > NRBINS ) cycle
                    radprof(irb,ib) = radprof(irb,ib) + abs(rmat(i,j,l))
                    if( ib == 1 ) radcnt(irb) = radcnt(irb) + 1
                enddo
            enddo
        enddo
        nullify(rmat)
    enddo
    do irb = 1, NRBINS
        if( radcnt(irb) > 0 ) radprof(irb,:) = radprof(irb,:) / real(radcnt(irb))
    enddo
    mskrad  = 0.5 * cline%get_rarg('mskdiam') / smpd_out
    nrb_msk = max(1, min(NRBINS, int(mskrad / rbin_width) + 1))
    ! ground-truth mode: radial |rho| profiles against a known volume (synthetic data),
    ! all maps low-passed identically
    l_truth = cline%defined('vol1')
    if( l_truth )then
        truth_fname = cline%get_carg('vol1')
        call truth%new(ldim, smpd_out)
        call truth%read(truth_fname)
        lp_here = 0.
        if( cline%defined('lp') ) lp_here = cline%get_rarg('lp')
        hp_here = 0.
        if( cline%defined('hp') ) hp_here = cline%get_rarg('hp')
        call tvols(1)%copy(vols(1))
        call tvols(2)%copy(vols(2))
        ! Fourier-shell comparison against the truth: background (mean outside the mask)
        ! removed, then soft-masked, unfiltered. Without the background removal a map
        ! with a non-zero solvent level (gridding; PCG's support-masked map has none)
        ! acquires a sphere-shaped term from the mask whose spectrum sits in shells 1-3.
        call truth%get_rmat_ptr(rmat)
        bg = background_mean(rmat); rmat(1:ldim(1),1:ldim(2),1:ldim(3)) = rmat(1:ldim(1),1:ldim(2),1:ldim(3)) - bg
        nullify(rmat)
        call truth%mask3D_soft(mskrad)
        call truth%fft
        call truth%spectrum('sqrt', tspec)
        allocate(tcorr(size(tspec),2), tspec_b(size(tspec),2), source=0.)
        do it = 1, 2
            call tvols(it)%get_rmat_ptr(rmat)
            bg = background_mean(rmat); rmat(1:ldim(1),1:ldim(2),1:ldim(3)) = rmat(1:ldim(1),1:ldim(2),1:ldim(3)) - bg
            nullify(rmat)
            write(logfhandle,'(A,A,A,ES11.4)') '>>> REC3D BACKENDS: TRUTH background level ', trim(BACKENDS(it)), ': ', bg
            call tvols(it)%mask3D_soft(mskrad)
            call tvols(it)%fft
            call tvols(it)%spectrum('sqrt', spec_tmp)
            tspec_b(:,it) = spec_tmp
            call truth%fsc(tvols(it), tcorr(:,it))
            call tvols(it)%ifft
        enddo
        call truth%ifft
        write(logfhandle,'(A)') '>>> REC3D BACKENDS: TRUTH SHELL TABLE (soft-masked)  k  res(A)  fsc(truth,gridding)  fsc(truth,pcg)  amp_gridding/truth  amp_pcg/truth'
        do k = 1, size(tspec)
            write(logfhandle,'(A,I4,F9.2,2F12.4,2F12.4)') '>>> REC3D BACKENDS: TRUTH SHELL ', k, truth%get_lp(k), &
                &tcorr(k,1), tcorr(k,2), safe_ratio(tspec_b(k,1), tspec(k)), safe_ratio(tspec_b(k,2), tspec(k))
        enddo
        ! re-read the unmasked maps for the radial comparison
        call truth%read(truth_fname)
        call tvols(1)%copy(vols(1))
        call tvols(2)%copy(vols(2))
        if( lp_here > 0. .or. hp_here > 0. )then
            call truth%bp(hp_here, lp_here)
            do it = 1, 2
                call tvols(it)%bp(hp_here, lp_here)
            enddo
        endif
        ! background: mean outside the mask (+4 px), removed from every map before comparison;
        ! per-shell least-squares scale <recon*truth>/<truth*truth> is immune to additive offsets
        call truth%get_rmat_ptr(rmat)
        bg = background_mean(rmat)
        rmat(1:ldim(1),1:ldim(2),1:ldim(3)) = rmat(1:ldim(1),1:ldim(2),1:ldim(3)) - bg
        allocate(tprof(NRBINS,2), source=0.)
        do it = 1, 2
            call tvols(it)%get_rmat_ptr(rmat)
            bg   = background_mean(rmat)
            rmat(1:ldim(1),1:ldim(2),1:ldim(3)) = rmat(1:ldim(1),1:ldim(2),1:ldim(3)) - bg
            nullify(rmat)
        enddo
        call truth%get_rmat_ptr(rmat)
        do it = 1, 2
            call ls_scale_profile(tvols(it), rmat, tprof(:,it))
        enddo
        nullify(rmat)
        write(logfhandle,'(A)') '>>> REC3D BACKENDS: TRUTH TABLE (per-shell LS scale recon/truth after background removal, normalised to bin 2; hp '//&
            &real2str_trim(hp_here)//' lp '//real2str_trim(lp_here)//' A)  bin  r_lo-r_hi(px)  gridding  pcg'
        do irb = 1, NRBINS
            write(logfhandle,'(A,I4,F7.1,A,F6.1,2F12.4,A)') '>>> REC3D BACKENDS: TRUTH ', irb, &
                &real(irb-1)*rbin_width, ' -', real(irb)*rbin_width, &
                &safe_ratio(tprof(irb,1), tprof(2,1)), &
                &safe_ratio(tprof(irb,2), tprof(2,2)), merge('  (inside mask)', '               ', irb <= nrb_msk)
        enddo
        call truth%kill
        do it = 1, 2
            call tvols(it)%kill
        enddo
    endif
    ! Fourier shell amplitudes and FSC between the backends (both soft-masked: the PCG
    ! solution is support-masked by the solver, the gridding map is not)
    do ib = 1, 2
        call vols(ib)%get_rmat_ptr(rmat)
        bg = background_mean(rmat); rmat(1:ldim(1),1:ldim(2),1:ldim(3)) = rmat(1:ldim(1),1:ldim(2),1:ldim(3)) - bg
        nullify(rmat)
        call vols(ib)%mask3D_soft(mskrad)
        call vols(ib)%fft
        call vols(ib)%spectrum('sqrt', spec_tmp)
        if( ib == 1 ) allocate(spec(size(spec_tmp),2), source=0.)
        spec(:,ib) = spec_tmp
    enddo
    lfny = size(spec, dim=1)
    allocate(corrs(lfny), source=0.)
    call vols(1)%fsc(vols(2), corrs)
    ! report
    write(logfhandle,'(A)') ''
    write(logfhandle,'(A,I0,A,I0,A,F0.4,A)') '>>> REC3D BACKENDS: STATE ', state, ' BOX ', ldim(1), ' SMPD ', smpd_out, &
        &'  MAPS: '//out_fnames(1)%to_char()//' '//out_fnames(2)%to_char()
    write(logfhandle,'(A,ES11.4,A,ES11.4,A,F0.4)') '>>> REC3D BACKENDS: L2 gridding ', l2(1), ' pcg ', l2(2), &
        &' pcg/gridding ', safe_ratio(l2(2), l2(1))
    write(logfhandle,'(A)') '>>> REC3D BACKENDS: SHELL TABLE  k  res(A)  amp_gridding  amp_pcg  pcg/gridding  fsc(gridding,pcg)'
    do k = 1, lfny
        write(logfhandle,'(A,I4,F9.2,2ES14.4,F12.4,F10.4)') '>>> REC3D BACKENDS: SHELL ', k, vols(1)%get_lp(k), &
            &spec(k,1), spec(k,2), safe_ratio(spec(k,2), spec(k,1)), corrs(k)
    enddo
    ! normalise the radial ratio to the MEDIAN over the in-mask bins: the centre
    ! bin has few voxels and carries the known PCG centre deficit (S6), so
    ! normalising to it inflates every other bin and falsifies the flatness gate
    allocate(ratios(max(size(spec,dim=1),NRBINS)), source=0.)
    n = 0
    do irb = 1, nrb_msk
        rr = safe_ratio(radprof(irb,2), radprof(irb,1))
        if( rr > 0. )then
            n = n + 1
            ratios(n) = rr
        endif
    enddo
    rnorm = 1.
    if( n > 0 )then
        call hpsort(ratios(1:n))
        rnorm = ratios(max(1,(n+1)/2))
    endif
    write(logfhandle,'(A)') '>>> REC3D BACKENDS: RADIAL TABLE  bin  r_lo-r_hi(px)  |rho|_gridding  |rho|_pcg  pcg/gridding  norm_to_med'
    do irb = 1, NRBINS
        rr = safe_ratio(radprof(irb,2), radprof(irb,1))
        write(logfhandle,'(A,I4,F7.1,A,F6.1,2ES14.4,2F12.4,A)') '>>> REC3D BACKENDS: RADIAL ', irb, &
            &real(irb-1)*rbin_width, ' -', real(irb)*rbin_width, radprof(irb,1), radprof(irb,2), rr, &
            &safe_ratio(rr, rnorm), merge('  (inside mask)', '               ', irb <= nrb_msk)
    enddo
    ! summary
    ! agreement band: contiguous shells from k=2 with FSC(gridding,pcg) > 0.5
    kagree = 1
    do k = 2, lfny
        if( corrs(k) <= 0.5 ) exit
        kagree = k
    enddo
    ratios = 0.
    n = 0
    do k = 2, kagree
        if( spec(k,1) > 0. .and. spec(k,2) > 0. )then
            n = n + 1
            ratios(n) = spec(k,2) / spec(k,1)
        endif
    enddo
    med = -1.
    if( n > 0 )then
        call hpsort(ratios(1:n))
        med = ratios(max(1, (n+1)/2))
    endif
    ! radial flatness: bins lying fully inside 0.85 x mask radius (clear of the PCG
    ! soft support edge); the centre bin is excluded from the range and reported
    ! separately (known PCG centre deficit, S6)
    summary%centre_bin = safe_ratio(safe_ratio(radprof(1,2), radprof(1,1)), rnorm)
    write(logfhandle,'(A,F0.4)') '>>> REC3D BACKENDS: PCG CENTRE-BIN RATIO (diagnostic, not gated): ', &
        &summary%centre_bin
    rmin = huge(rmin); rmax = -huge(rmax); nrb_used = 0
    do irb = 2, nrb_msk
        if( real(irb)*rbin_width > 0.85*mskrad ) exit
        if( radprof(irb,1) < 1.e-3*radprof(1,1) .or. radprof(irb,2) < 1.e-3*radprof(1,2) ) cycle
        rr = safe_ratio(safe_ratio(radprof(irb,2), radprof(irb,1)), rnorm)
        if( rr <= 0. ) cycle
        nrb_used = nrb_used + 1
        rmin = min(rmin, rr); rmax = max(rmax, rr)
    enddo
    if( nrb_used > 0 )then
        summary%rad_min = rmin
        summary%rad_max = rmax
    endif
    write(logfhandle,'(A,I0,A,F0.2,A,F0.4)') '>>> REC3D BACKENDS: SUMMARY agreement band (FSC gridding/pcg > 0.5) k=2-', kagree, &
        &' (', vols(1)%get_lp(kagree), ' A); median shell amplitude ratio pcg/gridding in band: ', med
    write(logfhandle,'(A,F0.4,A,F0.4,A,I0,A)') '>>> REC3D BACKENDS: SUMMARY radial ratio inside 0.85 x mask radius (centre bin excluded), normalised to the in-mask median: min ', rmin, &
        &' max ', rmax, ' (', nrb_used, ' bins)'
    write(logfhandle,'(A)') '>>> REC3D BACKENDS: EXPECTATION with the box division mirrored by PCG: shell ratio ~1, radial ratio rising toward the edge (gridding under-deapodized, S2.1)'
    write(logfhandle,'(A)') '>>> REC3D BACKENDS: EXPECTATION after dropping the division and fixing gridding deapodization: shell ratio ~1 AND radial ratio flat (~1)'
    ! GATES: violations are hard failures (review S4.4); thresholds from the validated
    ! neutral-phantom fixture and the streptavidin reference runs in the plan document.
    ! The gated band is the agreement band capped at the data band (lp, when given):
    ! beyond the data band the backends correlate with each other on shared noise and
    ! PCG's beyond-band behaviour is a known, separately tracked solver item (S6).
    kgate = kagree
    if( cline%defined('lp') )then
        kgate = min(kagree, calc_fourier_index(cline%get_rarg('lp'), ldim(1), smpd_out))
    else if( summary%p_base0143 > 0. )then
        ! no lp given: cap the gated band at the base-pair FSC=0.143
        ! resolution so the agreement gates are never judged on
        ! beyond-band noise where the backends legitimately decorrelate
        kgate = min(kagree, calc_fourier_index(summary%p_base0143, ldim(1), smpd_out))
    endif
    ! median amplitude ratio and FSC over the gated band
    allocate(fscs(lfny), source=0.)
    n = 0
    do k = 2, kgate
        if( spec(k,1) > 0. .and. spec(k,2) > 0. )then
            n = n + 1
            ratios(n) = spec(k,2) / spec(k,1)
        endif
    enddo
    med = -1.
    if( n > 0 )then
        call hpsort(ratios(1:n))
        med = ratios(max(1, (n+1)/2))
    endif
    summary%band_ratio = med
    write(logfhandle,'(A,I0,A,F0.4)') '>>> REC3D BACKENDS: GATED BAND k=2-', kgate, &
        &'; median shell amplitude ratio pcg/gridding in gated band: ', med
    nfail = 0
    if( kgate < 10 ) call gate_fail('gated band (agreement band capped at lp) ends at k='//int2str(kgate)//' (need >= 10)')
    if( n < 5 )      call gate_fail('only '//int2str(n)//' valid shells in the gated band (need >= 5)')
    if( med < 0.67 .or. med > 1.5 ) &
        &call gate_fail('median gated-band amplitude ratio pcg/gridding '//real2str_trim(med)//' outside [0.67,1.5]')
    nfsc = max(0, kgate - 1)
    if( nfsc > 0 )then
        fscs(1:nfsc) = corrs(2:kgate)
        call hpsort(fscs(1:nfsc))
        med_fsc = fscs(max(1,(nfsc+1)/2))
        summary%band_fsc = med_fsc
        if( med_fsc < 0.9 ) call gate_fail('median gated-band FSC(gridding,pcg) '//real2str_trim(med_fsc)//' below 0.9')
    endif
    if( nrb_used < 3 ) call gate_fail('only '//int2str(nrb_used)//' usable radial bins inside 0.85 x mask radius (need >= 3)')
    if( nrb_used >= 3 .and. (rmin < 0.5 .or. rmax > 2.0) ) &
        &call gate_fail('normalised radial ratio range ['//real2str_trim(rmin)//','//real2str_trim(rmax)//'] outside [0.5,2.0]')
    if( l_truth )then
        ! gridding LS profile vs truth must be flat inside the mask (a fading deapodization fails); ml_reg=yes
        ! maps move it legitimately (the gate is calibrated on unregularized maps), so there it is only reported.
        l_gate_ls = .true.
        if( cline%defined('ml_reg') )then
            if( cline%get_carg('ml_reg') .eq. 'yes' ) l_gate_ls = .false.
        endif
        do irb = 2, nrb_msk
            if( real(irb)*rbin_width > 0.85*mskrad ) exit
            tg = safe_ratio(tprof(irb,1), tprof(2,1))
            if( tg < 0.92 .or. tg > 1.08 )then
                if( l_gate_ls )then
                    call gate_fail('gridding/truth LS profile '//real2str_trim(tg)//' at radial bin '//int2str(irb)//' outside [0.92,1.08]')
                else
                    write(logfhandle,'(a)') '>>> REC3D BACKENDS: GRIDDING/TRUTH LS PROFILE '//real2str_trim(tg)//&
                        &' AT RADIAL BIN '//int2str(irb)//' outside [0.92,1.08] (diagnostic, not gated: ml_reg=yes)'
                endif
            endif
        enddo
        nfsc = max(0, kgate - 1)
        if( nfsc > 0 )then
            fscs(1:nfsc) = tcorr(2:kgate,1)
            call hpsort(fscs(1:nfsc))
            med_fsc = fscs(max(1,(nfsc+1)/2))
            summary%truth_fsc_g = med_fsc
            if( med_fsc < 0.8 ) call gate_fail('median FSC(truth,gridding) over the gated band '//real2str_trim(med_fsc)//' below 0.8')
            fscs(1:nfsc) = tcorr(2:kgate,2)
            call hpsort(fscs(1:nfsc))
            summary%truth_fsc_p = fscs(max(1,(nfsc+1)/2))
        endif
    endif
    if( nfail == 0 ) write(logfhandle,'(A)') '>>> REC3D BACKENDS: PASS (all gates)'
    do ib = 1, 2
        call vols(ib)%kill
    enddo
    if( allocated(spec)     ) deallocate(spec)
    if( allocated(spec_tmp) ) deallocate(spec_tmp)
    if( allocated(corrs)    ) deallocate(corrs)
    if( allocated(radprof)  ) deallocate(radprof)
    if( allocated(radcnt)   ) deallocate(radcnt)
    if( allocated(ratios)   ) deallocate(ratios)
    if( allocated(fscs)     ) deallocate(fscs)
    summary%nfail = nfail
    if( l_mkdir ) call simple_chdir(cwd_orig)
    call cwd_orig%kill
    if( nfail > 0 .and. l_abort_on_fail )then
        THROW_HARD('TEST_REC3D_BACKENDS FAILED: '//int2str(nfail)//' gate(s) violated (see >>> REC3D BACKENDS: FAIL lines)')
    endif

    contains

        subroutine gate_fail( msg )
            character(len=*), intent(in) :: msg
            nfail = nfail + 1
            write(logfhandle,'(A)') '>>> REC3D BACKENDS: FAIL -- '//trim(msg)
        end subroutine gate_fail

        pure real function safe_ratio( a, b )
            real, intent(in) :: a, b
            if( abs(b) > 0. )then
                safe_ratio = a / b
            else
                safe_ratio = -1.
            endif
        end function safe_ratio

        real function background_mean( rm )
            real, intent(in) :: rm(:,:,:)
            integer :: ii, jj, ll, cnt
            real    :: rad, acc
            acc = 0.; cnt = 0
            do ll = 1, ldim(3)
                do jj = 1, ldim(2)
                    do ii = 1, ldim(1)
                        rad = sqrt(real((ii-c(1))**2 + (jj-c(2))**2 + (ll-c(3))**2))
                        if( rad < mskrad + 4. .or. rad > real(ldim(1)/2) ) cycle
                        acc = acc + rm(ii,jj,ll); cnt = cnt + 1
                    enddo
                enddo
            enddo
            background_mean = 0.
            if( cnt > 0 ) background_mean = acc / real(cnt)
        end function background_mean

        subroutine ls_scale_profile( vol_in, tm, prof )
            class(image), intent(inout) :: vol_in
            real,         intent(in)    :: tm(:,:,:)
            real,         intent(out)   :: prof(NRBINS)
            real, pointer :: rm(:,:,:) => null()
            real(dp) :: num(NRBINS), den(NRBINS)
            integer  :: ii, jj, ll, ib_loc
            real     :: rad
            call vol_in%get_rmat_ptr(rm)
            num = 0.d0; den = 0.d0
            do ll = 1, ldim(3)
                do jj = 1, ldim(2)
                    do ii = 1, ldim(1)
                        rad    = sqrt(real((ii-c(1))**2 + (jj-c(2))**2 + (ll-c(3))**2))
                        ib_loc = int(rad / rbin_width) + 1
                        if( ib_loc > NRBINS ) cycle
                        num(ib_loc) = num(ib_loc) + real(rm(ii,jj,ll),dp) * real(tm(ii,jj,ll),dp)
                        den(ib_loc) = den(ib_loc) + real(tm(ii,jj,ll),dp)**2
                    enddo
                enddo
            enddo
            nullify(rm)
            prof = 0.
            where( den > 0.d0 ) prof = real(num / den)
        end subroutine ls_scale_profile

        function real2str_trim( x ) result( str )
            real, intent(in) :: x
            character(len=:), allocatable :: str
            character(len=32) :: buf
            write(buf,'(F0.1)') x
            str = trim(adjustl(buf))
        end function real2str_trim

        !> trimmed character token of a cline string argument, for the
        !! execution-directory name (get_carg's string carries padding)
        function cline_tok( key ) result( tok )
            character(len=*), intent(in)  :: key
            character(len=:), allocatable :: tok
            type(string) :: sval
            sval = cline%get_carg(key)
            tok  = trim(adjustl(sval%to_char()))
            call sval%kill
        end function cline_tok

        !> compact real token for the execution-directory name: plain decimal
        !! with trailing zeros stripped in the human range, scientific outside
        !! it (real2str_trim's F0.1 renders 1e-3 as '.0')
        !> Diagnostic FSC of a written half pair (soft spherical mask, same
        !! resolution readout as the production FSC path)
        subroutine report_pair_fsc( backend, l_unfil, label, r05_out, r0143_out )
            character(len=*), intent(in)  :: backend, label
            logical,          intent(in)  :: l_unfil
            real,             intent(out) :: r05_out, r0143_out
            type(image)  :: he, ho
            type(string) :: fe, fo
            real, allocatable :: corrs(:), res_h(:)
            integer :: ldim_h(3), nvols, nyq
            real    :: smpd_h, r05, r0143, mskrad_h
            r05_out   = 0.
            r0143_out = 0.
            fe = refine3D_state_halfvol_fname(state, 'even', unfil=l_unfil)
            fo = refine3D_state_halfvol_fname(state, 'odd',  unfil=l_unfil)
            if( .not. (file_exists(fe) .and. file_exists(fo)) )then
                call fe%kill
                call fo%kill
                return
            endif
            call find_ldim_nptcls(fe, ldim_h, nvols)
            smpd_h = smpd * real(box) / real(ldim_h(1))
            call he%new(ldim_h, smpd_h)
            call he%read(fe)
            call ho%new(ldim_h, smpd_h)
            call ho%read(fo)
            mskrad_h = 0.5 * cline%get_rarg('mskdiam') / smpd_h
            call he%mask3D_soft(mskrad_h, backgr=0.0)
            call ho%mask3D_soft(mskrad_h, backgr=0.0)
            ! fsc reads cmat: transform explicitly (the production path gets
            ! this implicitly from the preceding cFAR computation)
            call he%fft()
            call ho%fft()
            nyq = he%get_filtsz()
            allocate(corrs(nyq), source=0.0)
            call he%fsc(ho, corrs)
            res_h = he%get_res()
            call get_resolution(corrs, res_h, r05, r0143)
            r05   = max(r05,   2.0*smpd_h)
            r0143 = max(r0143, 2.0*smpd_h)
            r05_out   = r05
            r0143_out = r0143
            write(logfhandle,'(a,F8.3,a,F8.3)') '>>> REC3D BACKENDS: '//trim(backend)//' '//label//&
                &' half-pair FSC=0.500 at ', r05, '  FSC=0.143 at ', r0143
            call he%kill
            call ho%kill
            call fe%kill
            call fo%kill
            deallocate(corrs, res_h)
        end subroutine report_pair_fsc

        function real_tok( x ) result( tok )
            real, intent(in) :: x
            character(len=:), allocatable :: tok
            character(len=32) :: buf
            if( x == 0.0 )then
                tok = '0'
            else if( abs(x) >= 0.01 .and. abs(x) < 1000.0 )then
                write(buf,'(F0.3)') x
                tok = trim(adjustl(buf))
                do while( len(tok) > 1 .and. tok(len(tok):len(tok)) == '0' )
                    tok = tok(1:len(tok)-1)
                end do
                if( tok(len(tok):len(tok)) == '.' ) tok = tok(1:len(tok)-1)
            else
                write(buf,'(ES9.1)') x
                tok = trim(adjustl(buf))
            endif
        end function real_tok

end subroutine run_rec3D_backends_single

end submodule simple_commanders_test_highlevel_backends
