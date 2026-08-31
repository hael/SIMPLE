!@descr: flex_pca EM: polar pose refinement, pose perturbation and banded plane projections
submodule (simple_flex_pca_em) simple_flex_pca_em_pose
use simple_matcher_3Drec,   only: init_rec, cleanup_rec_buffers
use simple_matcher_ptcl_io, only: discrete_read_imgbatch, prepimgbatch
use simple_flex_reconstructor_latent_ops, only: latent_projection_weights, weighted_expanded_cmat,&
    &LATENT_WDIM
use simple_flex_projected_latent_model, only: prep_imgs4projected_model
use simple_flex_gpu,        only: flex_gpu_available, flex_gpu_prep_begin_f, flex_gpu_prep_free_f,&
    &flex_gpu_psample_begin_f, flex_gpu_psample_batch_f, flex_gpu_psample_free_f,&
    &flex_gpu_psample_batch_res_f
use simple_flex_pca_polar,  only: polar_grid_build, polar_grid_kill, polar_project_recs,&
    &polar_sample_particle, polar_relative_inplane, polar_assign_directions, polar_sample_at_pose,&
    &polar_apply_shift, polar_dir_neighbours
implicit none
#include "simple_local_flags.inc"

contains



    !> Is the polar (shared-direction) former requested for the reduced solve?
    logical module function cov_polar_enabled()
        integer :: v
        v = 0
        call cov_env_int('SIMPLE_COV_POLAR', v)
        cov_polar_enabled = v > 0
    end function cov_polar_enabled

    !> ... and for the embedding? Defaults to following SIMPLE_COV_POLAR, but is separable so the two
    !! stages can be A/B'd one at a time -- changing both at once makes a regression unattributable.
    logical module function cov_polar_embed_enabled()
        integer :: v
        v = 0
        call cov_env_int('SIMPLE_COV_POLAR_EMBED', v)
        if( cov_env_is_set('SIMPLE_COV_POLAR_EMBED') )then
            cov_polar_embed_enabled = v > 0
        else
            cov_polar_embed_enabled = cov_polar_enabled()
        endif
    end function cov_polar_embed_enabled

    !> Number of bank directions.
    !!
    !! MEASURED (IgG-RL 10k, box_crop 64, lp 15, d_tilde 64): the largest reduced eigenvalue is
    !! 1194.2 / 1195.5 / 1196.2 / 1198.2 at ndir = 1000 / 2000 / 8000 / 32000 and 1198.3 with the
    !! grid removed entirely (SIMPLE_COV_POLAR_EXACT=1), and ground-truth basis capture is 0.5828
    !! at ndir=2000 against 0.5819 exact and 0.5848 Cartesian. A 6 degree direction grid is
    !! indistinguishable from no discretisation at all here, because this stage never uses a
    !! per-particle b on its own -- it accumulates Sbb and sum_i G_i (x) G_i over 10^4-10^5
    !! particles, and the Gram is additionally a sum over ~10^3 plane samples, so direction error
    !! enters suppressed by 1/sqrt(nsamp) rather than as a per-sample decorrelation.
    !!
    !! So the default targets AMORTISATION (~40 particles per direction) rather than resolution,
    !! with a floor so small datasets still get a reasonable grid. Raise it with
    !! SIMPLE_COV_POLAR_NDIR if a dataset ever shows direction sensitivity -- the bank is streamed
    !! direction by direction, so ndir costs no memory, only bank-build time.
    integer module function cov_polar_ndir( nptcls )
        integer, intent(in) :: nptcls
        integer :: v
        v = min(4000, max(1000, nptcls/40))
        v = 2*((v+1)/2)                                  ! build_refspiral needs an even count
        call cov_env_int('SIMPLE_COV_POLAR_NDIR', v)
        v = 2*((v+1)/2)
        cov_polar_ndir = v
    end function cov_polar_ndir




    !> ---------------------------------------------------------------------------------------------
    !> POLAR E-STEP FOR THE EMBEDDING
    !>
    !> Fills exactly the sufficient statistics the Cartesian batch loop of
    !> embed_latents_with_contrast produces -- Gcache, bcache, ccache, contrast, the two residual
    !> energies and the split-half solves zhalf -- so the reliability prior, the re-solve, the
    !> precision matrices and every downstream consumer are untouched.
    !>
    !> Same two-pass shape as the solve's former. The one thing the solve did not need is the
    !> SPLIT-HALF statistics `rho` is built from. Those come out of the grid's ring half-split
    !> (see polar_grid_t%rmid): each ring stores its even angles then its odd ones, so both halves
    !> are contiguous and each half's Gram and b are still a single BLAS call.
    !> ---------------------------------------------------------------------------------------------
    module subroutine embed_accumulate_polar( params, build, mean_rec, basis_recs, ncomp, eigvals, sig2, &
        &pinds, nptcls, Gcache, bcache, ccache, zhalf, contrast, resid_energy, resid_mean_energy, &
        &prior, nvalid )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: mean_rec
        type(reconstructor), intent(inout) :: basis_recs(ncomp)
        integer,             intent(in)    :: ncomp, pinds(:), nptcls
        real(dp),            intent(in)    :: eigvals(ncomp), sig2, prior(ncomp)
        real(dp),            intent(inout) :: Gcache(ncomp,ncomp,nptcls), bcache(ncomp,nptcls)
        real(dp),            intent(inout) :: ccache(ncomp,nptcls), zhalf(nptcls,ncomp,2)
        real(dp),            intent(inout) :: contrast(nptcls), resid_energy(nptcls)
        real(dp),            intent(inout) :: resid_mean_energy(nptcls)
        integer,             intent(out)   :: nvalid
        type(polar_grid_t)             :: pg
        type(oris)                     :: dirs
        type(ori)                      :: o
        type(fplane_type), allocatable :: fpls(:)
        real,    allocatable :: rmat_p(:,:,:), nrm_p(:,:), nrm_b(:,:), rmat_b(:,:,:), cav(:), sav(:)
        real,    allocatable :: eul_b(:,:)
        logical :: l_exact_embed
        integer :: j_exact
        integer, allocatable :: dir_of(:), dcnt(:), dptr(:), dlist(:)
        logical, allocatable :: l_zero(:)
        real,    allocatable :: xws(:,:), wrs(:,:), xws1(:,:), xws2(:,:)
        real(dp),allocatable :: eyy(:)
        integer, allocatable :: dirs_in_chunk(:), nvalid_thr(:)
        ! per-thread pass-B workspace; index 0 of the half dimension is the FULL ring
        complex, allocatable :: Ubank(:,:,:)
        real,    allocatable :: Us(:,:,:), Xs(:,:,:), Xs1(:,:,:), Xs2(:,:,:), Reb(:,:,:,:), Csp(:,:,:)
        real(dp),allocatable :: Cf(:,:,:,:), Wmat(:,:,:), Gall(:,:,:,:), Corr(:,:,:,:), den_v(:,:)
        real(dp),allocatable :: c00(:,:), Cm0(:,:,:,:)
        real(dp),allocatable :: Ath(:,:,:), zth(:,:)
        ! pose-refinement workspace: the bank for EVERY direction, kept while the images stream past
        logical :: l_pose
        integer :: ipose_test, i_tmp
        real    :: pose_rot_amp, pose_sh_amp, pose_rot_step, pose_sh_step
        real(dp):: bankgb
        real,    allocatable :: Usall(:,:,:), pose_rot(:), pose_shx(:), pose_shy(:)
        integer, allocatable :: dnn(:,:), ndmove_thr(:)
        real(dp),allocatable :: dangle_thr(:), dq0_thr(:), dq1_thr(:)
        logical :: l_p2
        integer :: nnn_dir, jdacc
        real(dp),allocatable :: Cfall(:,:,:), Cm0all(:,:,:), c00all(:,:)
        real(dp),allocatable :: pose_e0(:), pose_e1(:)
        ! pose band guard
        real     :: sh_cap_A, sh_cap_px, rot_cap_d, smpd_c, mskrad_A, shmag
        real(dp) :: acc_sh_A, acc_rot_A
        integer  :: icap, nclamp_sh, nclamp_rot, npose_w
        real(dp), allocatable :: pose_s0(:), pose_s1(:), pose_ea(:), pose_sa(:)
        integer, allocatable :: pose_n(:), pose_t(:), pose_r(:)
        logical :: l_guard
        integer :: nthr, ithr, i, j, q, r, ir, ibatch, batchlims(2), batchsz, row
        integer :: ndir, nk, nsamp, nsamp2, kto, kfrom, knfrom, knto, pf, nc1
        integer :: hlo, hhi, klo, ph0, pk0, nyq_band, nyq_rec, mmax, m, id, ic, ndchunk
        integer :: jc, idir, ip, ii, i0, nsl, ih, nrb, ncc
        real    :: rmb(3,3), ca, sa, tazim, e3p, shp(2)
        type(string) :: fn_pose
        real(dp):: pw1, cnt1, cc, aa, e_yy, e_mm, myv, res
        integer(timer_int_kind) :: t_a, t_b, t_es
        real :: sec_eread, sec_eprep, sec_esamp
        ! device polar sampler (SIMPLE_COV_GPU=1); pose refinement forces the CPU path
        logical  :: l_gpu_eps, l_eps_res, l_devprep
        integer  :: vps_e
        real,    allocatable :: pwv_e(:), tazv_e(:)
        logical, allocatable :: vmask_e(:)
        real,   parameter :: A_LO_C = 0.1, A_HI_C = 5.0
        nthr   = nthr_glob
        nvalid = 0
        pf     = OSMPL_PAD_FAC
        call mean_rec%expand_exp
        do q = 1, ncomp
            call basis_recs(q)%expand_exp
        end do
        call init_rec(params, build, MAXIMGBATCHSZ, fpls)
        call prepimgbatch(params, build, MAXIMGBATCHSZ)
        call cov_dev_prep_start(params, build, l_devprep)
        batchlims = [1, min(nptcls, MAXIMGBATCHSZ)]
        call discrete_read_imgbatch(params, build, nptcls, pinds, batchlims)
        call prep_imgs4projected_model(params, build, batchlims(2), build%imgbatch(:batchlims(2)), &
            &pinds(1:batchlims(2)), fpls(:batchlims(2)), mskrad=cov_image_mask_radius(params))
        ph0 = lbound(fpls(1)%cmplx_plane,1); pk0 = lbound(fpls(1)%cmplx_plane,2)
        hlo = ceil_div (lbound(fpls(1)%cmplx_plane,1), pf)
        hhi = floor_div(ubound(fpls(1)%cmplx_plane,1), pf)
        klo = ceil_div (lbound(fpls(1)%cmplx_plane,2), pf)
        nyq_rec  = mean_rec%get_lfny(1)
        nyq_band = nyq_rec
        if( fpls(1)%nyq > 0 ) nyq_band = min(nyq_band, max(1, fpls(1)%nyq / pf))
        kfrom  = 1
        kto    = nyq_band
        knfrom = max(kto+1, nint(0.7*real(nyq_rec)))
        knto   = max(knfrom, nyq_rec - 2)
        call polar_grid_build(pg, kfrom, kto, knfrom, knto, hlo, hhi, klo, ph0, pk0)
        nsamp = pg%nsamp; nsamp2 = 2*nsamp; nk = pg%nk
        ! EXACT-DIRECTION MODE. SIMPLE_COV_POLAR_EXACT was only ever honoured by
        ! reduced_solve_accumulate_polar; setting it on this former was bit-identical, i.e. a dead
        ! knob. It is the control that separates the two things the banked path confounds -- the
        ! polar QUADRATURE from the direction QUANTISATION -- so it is worth having here. Every
        ! particle becomes its own bank direction, which costs the whole shared-direction
        ! amortisation (C_qr(k) then serves one particle instead of many) and is therefore a
        ! diagnostic, not an operating point.
        j_exact = 0
        call cov_env_int('SIMPLE_COV_POLAR_EXACT', j_exact)
        l_exact_embed = j_exact > 0
        ndir  = cov_polar_ndir(nptcls)
        if( l_exact_embed ) ndir = nptcls
        ! ---------------- orientations and direction assignment
        allocate(rmat_p(3,3,nptcls), nrm_p(3,nptcls), dir_of(nptcls), cav(nptcls), sav(nptcls))
        allocate(l_zero(nptcls), source=.false.)
        do i = 1, nptcls
            call build%spproj_field%get_ori(pinds(i), o)
            rmat_p(:,:,i) = o%get_mat()
            nrm_p(:,i)    = rmat_p(3,:,i)
            l_zero(i)     = o%isstatezero()
        end do
        call o%kill
        allocate(nrm_b(3,ndir), rmat_b(3,3,ndir), eul_b(3,ndir))
        if( l_exact_embed )then
            ! each particle is its own direction: the relative in-plane angle is the identity, so
            ! the polar sampler reads the particle plane at its own orientation exactly
            rmat_b = rmat_p
            do i = 1, nptcls
                nrm_b(:,i) = rmat_b(3,:,i)
                call build%spproj_field%get_ori(pinds(i), o)
                eul_b(:,i) = o%get_euler()
                dir_of(i)  = i
                cav(i)     = 1.0
                sav(i)     = 0.0
                if( l_zero(i) ) dir_of(i) = 0
            end do
            call o%kill
        else
        call dirs%new(ndir, is_ptcl=.false.)
        call build%pgrpsyms%build_refspiral(dirs)
        do j = 1, ndir
            rmat_b(:,:,j) = dirs%get_mat(j)
            nrm_b(:,j)    = rmat_b(3,:,j)
            eul_b(:,j)    = dirs%get_euler(j)
        end do
        call dirs%kill
        call polar_assign_directions(nrm_p, nptcls, nrm_b, ndir, dir_of)
        do i = 1, nptcls
            call polar_relative_inplane(rmat_p(:,:,i), rmat_b(:,:,dir_of(i)), ca, sa)
            cav(i) = ca; sav(i) = sa
            if( l_zero(i) ) dir_of(i) = 0
        end do
        endif
        write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> FLEX_PCA POLAR EMBED: samples=',nsamp, &
            &'  rings=',nk,'  directions=',ndir
        call flush(logfhandle)
        ! bank + ring-Gram scratch: allocated here because the pose block below needs them
        ! before pass A, and pass B reuses the same buffers
        allocate(Ubank(nsamp,0:ncomp,nthr), Us(nsamp2,0:ncomp,nthr))
        allocate(Csp(0:ncomp,0:ncomp,nthr))
        ! ---------------- POSE REFINEMENT SETUP (P1: in-plane angle + shift)
        !
        ! Pose search needs the model in hand while the IMAGE is in hand, so when it is on the bank
        ! is built for every direction up front rather than direction-major in pass B. At embedding
        ! rank (ncomp ~ 16-20, not d_tilde ~ 128-250) that bank is a few hundred MB, which is exactly
        ! why this is affordable here and was not in the solve.
        l_pose = cov_pose_mode() > 0
        ipose_test = 0
        call cov_env_int('SIMPLE_COV_POSE_TEST', ipose_test)
        l_p2 = .false.
        if( l_pose ) l_p2 = cov_pose_mode() >= 2
        ! injected-perturbation amplitudes and search steps, all in 0.1 units so they can be swept
        ! from the environment without a rebuild
        i_tmp = ipose_test; call cov_env_int('SIMPLE_COV_POSE_TEST_ROT', i_tmp)
        pose_rot_amp  = 0.1*real(i_tmp)
        i_tmp = ipose_test; call cov_env_int('SIMPLE_COV_POSE_TEST_SH',  i_tmp)
        pose_sh_amp   = 0.1*real(i_tmp)
        i_tmp = 20; call cov_env_int('SIMPLE_COV_POSE_ROTSTEP', i_tmp)
        pose_rot_step = 0.1*real(i_tmp)
        i_tmp = 10; call cov_env_int('SIMPLE_COV_POSE_SHSTEP',  i_tmp)
        pose_sh_step  = 0.1*real(i_tmp)
        ! The split-half acceptance test is ON for P1 and OFF for P2, because that is what the
        ! measurements say. P1's signal (dominated by shift) is strong enough that the guard rejects
        ! mostly noise: it halves the damage to correct shifts at no cost to recovery. P2's direction
        ! signal is weak, and the same guard throws away half the REAL improvements -- direction
        ! recovery 3.11 -> 2.77 deg with it against 3.11 -> 2.25 without, and the refit ends up worse
        ! than not refining at all (capture 0.546 vs 0.573). SIMPLE_COV_POSE_GUARD overrides.
        l_guard = .not. l_p2
        if( cov_env_is_set('SIMPLE_COV_POSE_GUARD') ) &
            &l_guard = .not. cov_env_int_off('SIMPLE_COV_POSE_GUARD')
        if( l_pose )then
            bankgb = 4.d0*real(nsamp2,dp)*real(ncomp+1,dp)*real(ndir,dp)/1.d9
            write(logfhandle,'(A,F7.2,A)') '>>> FLEX_PCA POLAR POSE: full bank ',bankgb,' GB'
            if( bankgb > 6.d0 )then
                write(logfhandle,'(A)') '>>> FLEX_PCA POLAR POSE DISABLED: bank too large; &
                    &lower SIMPLE_COV_POLAR_NDIR'
                l_pose = .false.
            endif
        endif
        nnn_dir = 1
        if( l_p2 )then
            nnn_dir = 12
            call cov_env_int('SIMPLE_COV_POSE_NNN', nnn_dir)
            nnn_dir = max(2, min(nnn_dir, ndir-1))
        endif
        allocate(dnn(nnn_dir,ndir))
        if( l_pose )then
            allocate(Usall(nsamp2,0:ncomp,ndir), Cfall(ncomp*ncomp,nk,ndir))
            allocate(Cm0all(ncomp,nk,ndir), c00all(nk,ndir))
            if( l_p2 )then
                call polar_dir_neighbours(ndir, nnn_dir, nrm_b, dnn)
                write(logfhandle,'(A,I0,A,F6.2,A)') '>>> FLEX_PCA POLAR POSE P2: direction search &
                    &over ',nnn_dir,' neighbours (grid spacing ~', &
                    &2.0*180.0/(PI*sqrt(real(ndir))),' deg)'
            else
                do id = 1, ndir
                    dnn(:,id) = id
                end do
            endif
            !$omp parallel do default(shared) schedule(dynamic) proc_bind(close) &
            !$omp& private(id,ithr,rmb,q,j,ir)
            do id = 1, ndir
                ithr = omp_get_thread_num() + 1
                rmb  = rmat_b(:,:,id)
                call polar_project_recs(mean_rec, basis_recs, ncomp, rmb, pg, Ubank(:,:,ithr))
                do q = 0, ncomp
                    do j = 1, nsamp
                        Usall(2*j-1,q,id) = pg%sqwq(j)*real(Ubank(j,q,ithr))
                        Usall(2*j,  q,id) = pg%sqwq(j)*aimag(Ubank(j,q,ithr))
                    end do
                end do
                do ir = 1, nk
                    call polar_ring_gram(Usall(1,0,id), nsamp2, ncomp, pg%rbeg(ir), &
                        &pg%rend(ir)-pg%rbeg(ir)+1, Csp(0,0,ithr), Cfall(1,ir,id), Cm0all(1,ir,id))
                    c00all(ir,id) = polar_ring_selfpower(Usall(1,0,id), nsamp2, pg%rbeg(ir), &
                        &pg%rend(ir)-pg%rbeg(ir)+1)
                end do
            end do
            !$omp end parallel do
        endif
        ! ---------------- PASS A
        t_a = tic()
        allocate(xws(nsamp2,nptcls), source=0.0)
        ! the two lattice-parity half-fields. Storing them costs 2 x nsamp2 x nptcls sp, which is the
        ! price of a reliability estimate that is not contaminated by shared interpolation support.
        allocate(xws1(nsamp2,nptcls), xws2(nsamp2,nptcls), source=0.0)
        allocate(wrs(nk,nptcls), source=0.0)
        allocate(eyy(nptcls), source=0.d0)
        allocate(pose_rot(nptcls), pose_shx(nptcls), pose_shy(nptcls), source=0.0)
        allocate(pose_e0(nthr), pose_e1(nthr), pose_s0(nthr), pose_s1(nthr), source=0.d0)
        allocate(pose_ea(nthr), pose_sa(nthr), source=0.d0)
        allocate(pose_n(nthr), pose_t(nthr), pose_r(nthr), ndmove_thr(nthr), source=0)
        allocate(dangle_thr(nthr), dq0_thr(nthr), dq1_thr(nthr), source=0.d0)
        sec_eread = 0.; sec_eprep = 0.; sec_esamp = 0.
        vps_e = 0
        call cov_env_int('SIMPLE_COV_GPU_PSAMPLE', vps_e)
        l_gpu_eps = vps_e > 0 .and. flex_gpu_available() .and. .not. l_pose
        ! resident sampler default (see the solve pass-A note); pose refinement forces CPU
        l_eps_res = (.not. l_gpu_eps) .and. flex_gpu_available() .and. .not. l_pose .and. &
            &.not. cov_env_int_off('SIMPLE_COV_GPU_PSAMPLE') .and. &
            &cov_image_mask_radius(params) <= 0.
        if( l_eps_res .and. params%l_ml_reg ) l_eps_res = allocated(build%esig%sigma2_noise)
        if( l_gpu_eps .or. l_eps_res )then
            call flex_gpu_psample_begin_f(pg%rad, pg%cs, pg%sn, pg%sqwq, pg%rbeg, pg%rend, &
                &pg%nsamp, pg%nk, pg%nrad, pg%ncs, pg%nsn, pg%nwq, pg%nsamp_n)
            allocate(pwv_e(MAXIMGBATCHSZ), tazv_e(MAXIMGBATCHSZ), vmask_e(MAXIMGBATCHSZ))
            if( l_eps_res )then
                call flex_gpu_prep_begin_f(build%lmsk, params%box, params%boxpd, &
                    &MAXIMGBATCHSZ, 0.0, .true.)
                write(logfhandle,'(A)') '>>> FLEX_PCA POLAR EMBED PASS-A RESIDENT device sampler ON'
            else
                write(logfhandle,'(A)') '>>> FLEX_PCA POLAR EMBED PASS-A device sampler ON'
            endif
            call flush(logfhandle)
        endif
        do ibatch = 1, nptcls, MAXIMGBATCHSZ
            batchlims = [ibatch, min(nptcls, ibatch + MAXIMGBATCHSZ - 1)]
            batchsz   = batchlims(2) - batchlims(1) + 1
            t_es = tic()
            call discrete_read_imgbatch(params, build, nptcls, pinds, batchlims)
            sec_eread = sec_eread + toc(t_es)
            t_es = tic()
            call prep_imgs4projected_model(params, build, batchsz, build%imgbatch(:batchsz), &
                &pinds(batchlims(1):batchlims(2)), fpls(:batchsz), &
                &mskrad=cov_image_mask_radius(params), resident=l_eps_res)
            sec_eprep = sec_eprep + toc(t_es)
            t_es = tic()
            if( l_gpu_eps .or. l_eps_res )then
                do i = 1, batchsz
                    vmask_e(i) = dir_of(batchlims(1)+i-1) > 0
                end do
                if( l_eps_res )then
                    call flex_gpu_psample_batch_res_f(cav(batchlims(1):batchlims(2)), &
                        &sav(batchlims(1):batchlims(2)), vmask_e(:batchsz), batchsz, .true., &
                        &xws, xws1, xws2, wrs, batchlims(1), pwv_e(:batchsz), tazv_e(:batchsz))
                else
                    call flex_gpu_psample_batch_f(fpls(:batchsz), cav(batchlims(1):batchlims(2)), &
                        &sav(batchlims(1):batchlims(2)), vmask_e(:batchsz), batchsz, .true., &
                        &xws, xws1, xws2, wrs, batchlims(1), pwv_e(:batchsz), tazv_e(:batchsz))
                endif
                !$omp parallel do default(shared) schedule(static) proc_bind(close) private(i,row)
                do i = 1, batchsz
                    row = batchlims(1) + i - 1
                    if( dir_of(row) <= 0 ) cycle
                    eyy(row) = polar_self_energy(xws(:,row), wrs(:,row), pg)
                end do
                !$omp end parallel do
                sec_esamp = sec_esamp + toc(t_es)
                cycle
            endif
            !$omp parallel do default(shared) schedule(static) proc_bind(close) &
            !$omp& private(i,row,ithr,tazim,pw1,cnt1,jdacc)
            do i = 1, batchsz
                row = batchlims(1) + i - 1
                if( dir_of(row) <= 0 ) cycle
                ithr = omp_get_thread_num() + 1
                call polar_sample_particle_packed(fpls(i), pg, cav(row), sav(row), &
                    &xws(:,row), wrs(:,row), pw1, cnt1, tazim, xws1(:,row), xws2(:,row))
                ! <y,y> in the polar measure; the Cartesian path takes it with cov_herm_inner
                eyy(row) = polar_self_energy(xws(:,row), wrs(:,row), pg)
                if( l_pose )then
                    jdacc = dir_of(row)
                    call polar_pose_refine_one(pg, fpls(i), Usall, Cfall, Cm0all, c00all, &
                        &nsamp2, ncomp, nk, ndir, dir_of(row), row, wrs(:,row), prior, sig2, &
                        &ipose_test, pose_rot_amp, pose_sh_amp, pose_rot_step, pose_sh_step, &
                        &l_guard, l_p2, nnn_dir, dnn(:,dir_of(row)), rmat_p(:,:,row), rmat_b, &
                        &jdacc, cav(row), sav(row), pose_rot(row), pose_shx(row), pose_shy(row), &
                        &xws(:,row), pose_e0(ithr), pose_e1(ithr), pose_s0(ithr), pose_s1(ithr), &
                        &pose_ea(ithr), pose_sa(ithr), pose_n(ithr), pose_t(ithr), pose_r(ithr))
                    ! a direction move changes which bank pass B must use for this particle
                    if( jdacc /= dir_of(row) )then
                        ndmove_thr(ithr) = ndmove_thr(ithr) + 1
                        dangle_thr(ithr) = dangle_thr(ithr) + real(acos(max(-1.,min(1., &
                            &dot_product(rmat_b(3,:,jdacc), rmat_b(3,:,dir_of(row))))))*180./PI, dp)
                    endif
                    ! recovery of a KNOWN direction degradation, if one was injected: the residual is
                    ! measured for BOTH the incoming and the accepted grid direction, so the grid's
                    ! own quantisation floor is present in both and cancels from the comparison
                    if( allocated(COV_TRUTH_NRM) )then
                        dq0_thr(ithr) = dq0_thr(ithr) + real(acos(max(-1.,min(1., &
                            &dot_product(rmat_b(3,:,dir_of(row)), COV_TRUTH_NRM(:,row)))))*180./PI, dp)**2
                        dq1_thr(ithr) = dq1_thr(ithr) + real(acos(max(-1.,min(1., &
                            &dot_product(rmat_b(3,:,jdacc), COV_TRUTH_NRM(:,row)))))*180./PI, dp)**2
                    endif
                    dir_of(row) = jdacc
                endif
            end do
            !$omp end parallel do
            sec_esamp = sec_esamp + toc(t_es)
        end do
        ! ---- POSE BAND GUARD, applied before anything is written back.
        !
        ! The pose block sees ONLY frequencies below `params%lp`, so it cannot justify a correction
        ! larger than that band can resolve. Nothing enforced this, and nothing reported the
        ! correction in physical units, which is how the following got through on EMPIAR-10330
        ! (PfCRT, published at 3.2 A): the block ran at the commander's default `box_crop=64` on a
        ! 300-box particle, i.e. `smpd_crop = 4.85 A/px`, and applied a 1.84 px shift -- which reads
        ! as small and is a **8.9 A** displacement. Consensus resolution went 3.83 A -> 7.06 A.
        !
        ! The same 1.84 px is 2.4 A at `smpd_crop = 1.3`. Pixels hide the failure; Angstroms do not.
        ! Everything below is therefore in Angstroms, and the cap is a fraction of the band: a
        ! displacement of `d` puts a phase error of `2*pi*d/lp` on the band edge, so `lp/4` is a
        ! quarter turn there and is already generous.
        !
        ! Rotation is converted to the displacement it causes at the mask radius, so one budget
        ! covers both degrees of freedom -- which is the right way to compare them
        ! (see the shift-vs-rotation result: they are one identifiability problem in displacement).
        if( l_pose )then
            sh_cap_A = 0.25 * params%lp
            icap = 0
            call cov_env_int('SIMPLE_COV_POSE_MAXSH', icap)      ! tenths of an Angstrom
            if( icap > 0 ) sh_cap_A = 0.1 * real(icap)
            smpd_c   = params%smpd_crop
            mskrad_A = 0.5 * params%msk_crop * smpd_c
            if( mskrad_A <= 0. ) mskrad_A = 0.25 * real(params%box_crop) * smpd_c
            sh_cap_px = sh_cap_A / max(smpd_c, TINY)
            rot_cap_d = sh_cap_A / max(mskrad_A, TINY) * 180. / PI
            nclamp_sh = 0; nclamp_rot = 0
            acc_sh_A  = 0.d0; acc_rot_A = 0.d0; npose_w = 0
            do i = 1, nptcls
                if( dir_of(i) <= 0 ) cycle
                npose_w  = npose_w + 1
                shmag    = sqrt(pose_shx(i)**2 + pose_shy(i)**2)
                acc_sh_A = acc_sh_A + real(shmag*smpd_c, dp)**2
                acc_rot_A= acc_rot_A + real(abs(pose_rot(i))*PI/180.*mskrad_A, dp)**2
                if( shmag > sh_cap_px )then
                    pose_shx(i) = pose_shx(i) * (sh_cap_px/shmag)
                    pose_shy(i) = pose_shy(i) * (sh_cap_px/shmag)
                    nclamp_sh   = nclamp_sh + 1
                endif
                if( abs(pose_rot(i)) > rot_cap_d )then
                    pose_rot(i) = sign(rot_cap_d, pose_rot(i))
                    nclamp_rot  = nclamp_rot + 1
                endif
            end do
            write(logfhandle,'(A,F7.2,A,F6.3,A,F7.2,A,F6.2,A)') &
                &'>>> FLEX_PCA POSE BAND GUARD: lp=',params%lp,' A  smpd_crop=',smpd_c, &
                &' A/px  cap=',sh_cap_A,' A (',sh_cap_px,' px)'
            write(logfhandle,'(A,F7.2,A,F7.2,A)') '>>>   accepted correction rms BEFORE cap: shift=', &
                &sqrt(sum_dp_safe(acc_sh_A, npose_w)),' A   in-plane edge displacement=', &
                &sqrt(sum_dp_safe(acc_rot_A, npose_w)),' A'
            write(logfhandle,'(A,I0,A,I0,A,I0,A,F6.1,A)') '>>>   clamped: shift ',nclamp_sh, &
                &'  rotation ',nclamp_rot,'  of ',npose_w,' particles (', &
                &100.*real(nclamp_sh)/real(max(1,npose_w)),' % of shifts)'
            if( npose_w > 0 .and. real(nclamp_sh)/real(npose_w) > 0.2 )then
                write(logfhandle,'(A,F6.1,A)') '>>> FLEX_PCA POSE WARNING: ', &
                    &100.*real(nclamp_sh)/real(npose_w),' % of shift corrections exceeded the band cap.'
                write(logfhandle,'(A)') '>>>   The pose block is proposing moves its own band cannot &
                    &justify. Either the incoming poses are much better than lp (do NOT refine), or &
                    &box_crop/lp are too coarse for this data. Compare a consensus reconstruction at &
                    &the refined poses against one at the incoming poses BEFORE trusting the result.'
            endif
            call flush(logfhandle)
        endif
        ! ---- WRITE THE REFINED POSES BACK.
        !
        ! The pose block lives in the EMBEDDING, which runs after the mean, the covariance columns
        ! and the basis. So refining here cannot change the delivered eigenvolumes at all -- measured:
        ! GT basis capture was identical to 4 decimal places with refinement on and off. To reach the
        ! basis the refined poses have to go back into the project and the fit has to be repeated;
        ! that outer loop is the joint refinement the design doc's section 4 describes, and this is
        ! the hook for it. Serial on purpose: the pose search runs inside an OpenMP region and the
        ! project field is shared.
        if( l_pose .and. cov_env_is_set('SIMPLE_COV_POSE_WRITE') )then
            do i = 1, nptcls
                if( dir_of(i) <= 0 ) cycle
                if( l_p2 )then
                    ! P2 may have moved the particle to a different bank direction. The accepted
                    ! orientation is that direction's (e1,e2) with the accepted in-plane angle --
                    ! taken from the bank's own Euler triplet rather than by inverting a rotation
                    ! matrix, so no m2euler branch/convention can silently corrupt it.
                    call build%spproj_field%set_euler(pinds(i), &
                        &[eul_b(1,dir_of(i)), eul_b(2,dir_of(i)), &
                        & eul_b(3,dir_of(i)) + polar_inplane_deg(cav(i), sav(i))])
                else
                    e3p = build%spproj_field%e3get(pinds(i))
                    call build%spproj_field%e3set(pinds(i), e3p + pose_rot(i))
                endif
                shp = build%spproj_field%get_2Dshift(pinds(i))
                call build%spproj_field%set_shift(pinds(i), shp + [pose_shx(i), pose_shy(i)])
            end do
            fn_pose = 'flex_pca_refined_poses.simple'
            call build%spproj%write(fn_pose)
            call fn_pose%kill
            write(logfhandle,'(A)') '>>> FLEX_PCA REFINED POSES WRITTEN: flex_pca_refined_poses.simple &
                &-- feed it back as projfile= for the next outer iteration'
            call flush(logfhandle)
        endif
        call cov_dev_prep_stop(l_devprep)
        if( l_eps_res ) call flex_gpu_prep_free_f
        write(logfhandle,'(A,F8.1)') '>>> FLEX_PCA POLAR EMBED PASS-A SECONDS: ',toc(t_a)
        write(logfhandle,'(A,F7.1,A,F7.1,A,F7.1)') '>>> FLEX_PCA POLAR EMBED PASS-A SPLIT (seconds): read=', &
            &sec_eread,'  prep=',sec_eprep,'  sampling=',sec_esamp
        if( l_gpu_eps .or. l_eps_res )then
            call flex_gpu_psample_free_f
            deallocate(pwv_e, tazv_e, vmask_e)
        endif
        if( l_pose .and. sum(pose_n) > 0 )then
            write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA POSE REFINEMENT (P1: in-plane + shift) on ', &
                &sum(pose_n),' particles'
            if( ipose_test > 0 )then
                write(logfhandle,'(A,F8.4,A,F8.4)') '>>>   PERTURB-RECOVER angle (deg): injected rms=', &
                    &sqrt(sum(pose_e0)/real(sum(pose_n),dp)),'  residual rms=', &
                    &sqrt(sum(pose_e1)/real(sum(pose_n),dp))
                write(logfhandle,'(A,F8.4,A,F8.4)') '>>>   PERTURB-RECOVER shift (px):  injected rms=', &
                    &sqrt(sum(pose_s0)/real(sum(pose_n),dp)),'  residual rms=', &
                    &sqrt(sum(pose_s1)/real(sum(pose_n),dp))
                write(logfhandle,'(A,F7.3,A)') '>>>   the TRUE pose outscores the search answer for ', &
                    &100.d0*real(sum(pose_t),dp)/real(sum(pose_n),dp), &
                    &' % of particles (0 % = search converged, ~50 % = objective uninformative)'
            else if( COV_PROJ_PERTURB_ROT > 0. .or. COV_PROJ_PERTURB_SH > 0. )then
                write(logfhandle,'(A,F8.4,A,F8.4,A,F8.4)') &
                    &'>>>   RECOVERY vs the PROJECT degradation, angle (deg): injected rms=', &
                    &sqrt(sum(pose_e0)/real(sum(pose_n),dp)),'  residual(+)=', &
                    &sqrt(sum(pose_e1)/real(sum(pose_n),dp)),'  residual(-)=', &
                    &sqrt(sum(pose_ea)/real(sum(pose_n),dp))
                write(logfhandle,'(A,F8.4,A,F8.4,A,F8.4)') &
                    &'>>>   RECOVERY vs the PROJECT degradation, shift (px):  injected rms=', &
                    &sqrt(sum(pose_s0)/real(sum(pose_n),dp)),'  residual(+)=', &
                    &sqrt(sum(pose_s1)/real(sum(pose_n),dp)),'  residual(-)=', &
                    &sqrt(sum(pose_sa)/real(sum(pose_n),dp))
            else
                write(logfhandle,'(A,F8.4,A,F8.4)') '>>>   accepted move rms: angle (deg)=', &
                    &sqrt(sum(pose_e1)/real(sum(pose_n),dp)),'  shift (px)=', &
                    &sqrt(sum(pose_s1)/real(sum(pose_n),dp))
                write(logfhandle,'(A,F7.3,A)') '>>>   the score improved for ', &
                    &100.d0*real(sum(pose_t),dp)/real(sum(pose_n),dp),' % of particles'
            endif
            write(logfhandle,'(A,F7.3,A,L2)') '>>>   split-half guard rejected ', &
                &100.d0*real(sum(pose_r),dp)/real(sum(pose_n),dp),' % of proposed moves; guard=',l_guard
            if( l_p2 ) write(logfhandle,'(A,F7.3,A,F7.3,A)') '>>>   direction changed for ', &
                &100.d0*real(sum(ndmove_thr),dp)/real(sum(pose_n),dp),' % of particles, mean move ', &
                &sum(dangle_thr)/max(real(sum(ndmove_thr),dp),1.d0),' deg'
            if( allocated(COV_TRUTH_NRM) ) write(logfhandle,'(A,F8.4,A,F8.4)') &
                &'>>>   DIRECTION recovery vs truth (deg): before rms=', &
                &sqrt(sum(dq0_thr)/real(sum(pose_n),dp)),'  after rms=', &
                &sqrt(sum(dq1_thr)/real(sum(pose_n),dp))
        endif
        call flush(logfhandle)
        call cleanup_rec_buffers(build, fpls)
        deallocate(nrm_p, nrm_b)
        ! ---------------- direction -> particle CSR
        allocate(dcnt(ndir), source=0)
        do i = 1, nptcls
            if( dir_of(i) > 0 ) dcnt(dir_of(i)) = dcnt(dir_of(i)) + 1
        end do
        allocate(dptr(ndir+1)); dptr(1) = 1
        do j = 1, ndir
            dptr(j+1) = dptr(j) + dcnt(j)
        end do
        allocate(dlist(max(1,dptr(ndir+1)-1)))
        dcnt = 0
        do i = 1, nptcls
            j = dir_of(i)
            if( j <= 0 ) cycle
            dlist(dptr(j)+dcnt(j)) = i
            dcnt(j) = dcnt(j) + 1
        end do
        mmax = 0
        do j = 1, ndir
            mmax = max(mmax, dcnt(j))
        end do
        mmax = min(max(64, mmax), 256)
        ! ---------------- PASS B
        t_b  = tic()
        ncc  = ncomp + 1                    ! bank slot 0 is the mean
        ! half index 0 = the full ring, 1 and 2 = the two interleaved half-sets
        allocate(Cf(ncomp*ncomp,nk,0:2,nthr), c00(nk,nthr), Cm0(ncomp,nk,0:2,nthr))
        allocate(Xs(nsamp2,mmax,nthr), Reb(0:ncomp,mmax,0:2,nthr), Wmat(nk,mmax,nthr))
        allocate(Xs1(nsamp2,mmax,nthr), Xs2(nsamp2,mmax,nthr))
        allocate(Gall(ncomp*ncomp,mmax,0:2,nthr), Corr(ncomp,mmax,0:2,nthr), den_v(mmax,nthr))
        allocate(Ath(ncomp,ncomp,nthr), zth(ncomp,nthr))
        allocate(nvalid_thr(nthr), source=0)
        ndchunk = max(4*nthr, 64)
        allocate(dirs_in_chunk(ndchunk))
        idir = 1
        do while( idir <= ndir )
            ic = 0
            do while( idir <= ndir .and. ic < ndchunk )
                if( dcnt(idir) > 0 )then
                    ic = ic + 1
                    dirs_in_chunk(ic) = idir
                endif
                idir = idir + 1
            end do
            if( ic == 0 ) cycle
            !$omp parallel do default(shared) schedule(dynamic) proc_bind(close) &
            !$omp& private(jc,id,m,ithr,rmb,j,q,r,ir,ip,ii,i0,nsl,ih,nrb,nc1,row) &
            !$omp& private(cc,aa,e_yy,e_mm,myv,res)
            do jc = 1, ic
                id   = dirs_in_chunk(jc)
                m    = dcnt(id)
                ithr = omp_get_thread_num() + 1
                rmb  = rmat_b(:,:,id)
                call polar_project_recs(mean_rec, basis_recs, ncomp, rmb, pg, Ubank(:,:,ithr))
                do q = 0, ncomp
                    do j = 1, nsamp
                        Us(2*j-1,q,ithr) = pg%sqwq(j)*real(Ubank(j,q,ithr))
                        Us(2*j,  q,ithr) = pg%sqwq(j)*aimag(Ubank(j,q,ithr))
                    end do
                end do
                ! Ring Gram tables. The half-set G is HALF the full one, not a subset sum over
                ! polar samples: the split is now by Cartesian lattice parity (see the b-halves
                ! below), each parity carries half the lattice measure, and the basis is noise-free
                ! and smooth so both parities see the same signal. That keeps A_half consistent with
                ! b_half, which is what the reliability solve needs.
                do ir = 1, nk
                    call polar_ring_gram(Us(1,0,ithr), nsamp2, ncomp, pg%rbeg(ir), &
                        &pg%rend(ir) - pg%rbeg(ir) + 1, &
                        &Csp(0,0,ithr), Cf(1,ir,0,ithr), Cm0(1,ir,0,ithr))
                    Cf(:,ir,1,ithr)  = 0.5d0*Cf(:,ir,0,ithr)
                    Cf(:,ir,2,ithr)  = Cf(:,ir,1,ithr)
                    Cm0(:,ir,1,ithr) = 0.5d0*Cm0(:,ir,0,ithr)
                    Cm0(:,ir,2,ithr) = Cm0(:,ir,1,ithr)
                    ! <T mu, T mu> per ring, from the mean's own diagonal entry
                    c00(ir,ithr) = polar_ring_selfpower(Us(1,0,ithr), nsamp2, pg%rbeg(ir), &
                        &pg%rend(ir) - pg%rbeg(ir) + 1)
                end do
                do i0 = 1, m, mmax
                    nsl = min(mmax, m - i0 + 1)
                    do ii = 1, nsl
                        ip = dlist(dptr(id)+i0+ii-2)
                        Xs(:,ii,ithr)   = xws(:,ip)
                        Xs1(:,ii,ithr)  = xws1(:,ip)
                        Xs2(:,ii,ithr)  = xws2(:,ip)
                        Wmat(:,ii,ithr) = real(wrs(:,ip), dp)
                    end do
                    ! b for the full plane and for each half. The halves are not contiguous ACROSS
                    ! rings, so each is accumulated ring by ring; the full one is one GEMM.
                    call sgemm('T','N', ncomp+1, nsl, nsamp2, 1.0, Us(1,0,ithr), nsamp2, &
                        &Xs(1,1,ithr), nsamp2, 0.0, Reb(0,1,0,ithr), ncomp+1)
                    ! Half-b from the two LATTICE-PARITY fields, over the whole sample list. The old
                    ! ring-interleaved split alternated polar samples one lattice unit apart, which
                    ! the 3-tap kernel makes share most of their support: measured, that pushed every
                    ! rho to 0.95-1.00 on EMPIAR-10076 and flattened the reliability prior.
                    call sgemm('T','N', ncomp+1, nsl, nsamp2, 1.0, Us(1,0,ithr), nsamp2, &
                        &Xs1(1,1,ithr), nsamp2, 0.0, Reb(0,1,1,ithr), ncomp+1)
                    call sgemm('T','N', ncomp+1, nsl, nsamp2, 1.0, Us(1,0,ithr), nsamp2, &
                        &Xs2(1,1,ithr), nsamp2, 0.0, Reb(0,1,2,ithr), ncomp+1)
                    ! G and the mean cross term, full and per half, all particles at once
                    do ih = 0, 2
                        call dgemm('N','N', ncomp*ncomp, nsl, nk, 1.d0, Cf(1,1,ih,ithr), ncomp*ncomp, &
                            &Wmat(1,1,ithr), nk, 0.d0, Gall(1,1,ih,ithr), ncomp*ncomp)
                        call dgemm('N','N', ncomp, nsl, nk, 1.d0, Cm0(1,1,ih,ithr), ncomp, &
                            &Wmat(1,1,ithr), nk, 0.d0, Corr(1,1,ih,ithr), ncomp)
                    end do
                    call dgemv('T', nk, nsl, 1.d0, Wmat(1,1,ithr), nk, c00(1,ithr), 1, 0.d0, &
                        &den_v(1,ithr), 1)
                    do ii = 1, nsl
                        row  = dlist(dptr(id)+i0+ii-2)
                        e_mm = den_v(ii,ithr)
                        myv  = real(Reb(0,ii,0,ithr),dp)
                        e_yy = eyy(row)
                        do q = 1, ncomp
                            bcache(q,row) = real(Reb(q,ii,0,ithr),dp)
                            ccache(q,row) = Corr(q,ii,0,ithr)
                            do r = 1, ncomp
                                Gcache(q,r,row) = Gall((r-1)*ncomp+q,ii,0,ithr)
                            end do
                        end do
                        resid_mean_energy(row) = e_yy - 2.d0*myv + e_mm
                        if( COV_UNIT_CONTRAST .and. .not. cov_fit_contrast_rt )then
                            cc = 1.d0
                        else
                            cc = max(real(A_LO_C,dp), min(real(A_HI_C,dp), myv/max(e_mm,DTINY)))
                        endif
                        contrast(row) = cc
                        aa = cc*cc
                        ! MAP solve at the chosen contrast, identical arithmetic to the Cartesian body
                        Ath(:,:,ithr) = (aa/sig2)*Gcache(:,:,row)
                        do q = 1, ncomp
                            Ath(q,q,ithr) = Ath(q,q,ithr) + prior(q)
                            zth(q,ithr)   = (cc*bcache(q,row) - aa*ccache(q,row))/sig2
                        end do
                        call spd_solve_dp(Ath(:,:,ithr), zth(:,ithr), ncomp)
                        res = e_yy + aa*e_mm - 2.d0*cc*myv + quad_form(Gcache(:,:,row), zth(:,ithr), ncomp)*aa
                        do q = 1, ncomp
                            res = res + 2.d0*aa*zth(q,ithr)*ccache(q,row) - 2.d0*cc*zth(q,ithr)*bcache(q,row)
                        end do
                        res = res/sig2
                        do q = 1, ncomp
                            res = res + prior(q)*zth(q,ithr)*zth(q,ithr)
                        end do
                        resid_energy(row) = res
                        ! the two half-data solves the reliability prior is built from
                        do ih = 1, 2
                            do q = 1, ncomp
                                do r = 1, ncomp
                                    Ath(q,r,ithr) = (aa/sig2)*Gall((r-1)*ncomp+q,ii,ih,ithr)
                                end do
                            end do
                            do q = 1, ncomp
                                Ath(q,q,ithr) = Ath(q,q,ithr) + prior(q)
                                zth(q,ithr)   = (cc*real(Reb(q,ii,ih,ithr),dp) &
                                    &- aa*Corr(q,ii,ih,ithr))/sig2
                            end do
                            call spd_solve_dp(Ath(:,:,ithr), zth(:,ithr), ncomp)
                            zhalf(row,:,ih) = zth(:,ithr)
                        end do
                    end do
                end do
                nvalid_thr(ithr) = nvalid_thr(ithr) + m
            end do
            !$omp end parallel do
        end do
        nvalid = sum(nvalid_thr)
        write(logfhandle,'(A,F8.1,A,I0)') '>>> FLEX_PCA POLAR EMBED PASS-B SECONDS: ',toc(t_b), &
            &'   particles=',nvalid
        call flush(logfhandle)
        call polar_grid_kill(pg)
        deallocate(xws, xws1, xws2, wrs, eyy, dir_of, cav, sav, dcnt, dptr, dlist, l_zero, rmat_b, rmat_p, eul_b)
        deallocate(pose_rot, pose_shx, pose_shy, pose_e0, pose_e1, pose_s0, pose_s1, pose_n, pose_t, pose_r)
        deallocate(pose_ea, pose_sa, ndmove_thr, dangle_thr, dq0_thr, dq1_thr)
        if( allocated(dnn) ) deallocate(dnn)
        if( allocated(Usall) ) deallocate(Usall, Cfall, Cm0all, c00all)
        deallocate(Ubank, Us, Cf, c00, Cm0, Xs, Reb, Wmat, Gall, Corr, den_v, Ath, zth, Csp)
        deallocate(nvalid_thr, dirs_in_chunk)
    end subroutine embed_accumulate_polar

    !> Recover the in-plane angle in degrees from the (cos,sin) pair the pose block carries.
    pure real module function polar_inplane_deg( ca, sa ) result( d )
        real, intent(in) :: ca, sa
        d = atan2(sa, ca) * 180.0 / PI
    end function polar_inplane_deg

    integer module function cov_pose_mode()
        integer :: v
        v = 0
        call cov_env_int('SIMPLE_COV_POSE', v)
        cov_pose_mode = v
    end function cov_pose_mode

    !> ---------------------------------------------------------------------------------------------
    !> P1 POSE REFINEMENT FOR ONE PARTICLE: in-plane angle and shift, scored by the MARGINAL
    !> likelihood with the conformation integrated out.
    !>
    !> Why the marginal and not the fitted residual: the failure mode of joint pose+heterogeneity
    !> refinement is that pose absorbs conformation -- a small rotation mimics a low-order
    !> eigenvolume and an alternating optimiser takes it, because the loss does not care which knob
    !> explains the residual. Integrating z out under its prior removes the knob: a pose is preferred
    !> only if it explains the data better AFTER the whole ensemble has been given its best shot at
    !> every competing pose.
    !>
    !> With  A = (a^2/sig2) G + Lambda^-1  and  u = (a b - a^2 c)/sig2,
    !>     -2 log p(y | pose)  =  [e_yy + a^2 e_mm - 2 a my]/sig2  -  u^T A^-1 u  +  log det A + const
    !> and for IN-PLANE and SHIFT moves G is invariant (the noise weight is radial), so `A` and the
    !> log-det are constant across trials. One Cholesky per particle serves the whole search, and
    !> each trial costs one GEMV plus one triangular solve.
    !>
    !> Two further consequences of working in polar coordinates: the in-plane angle is just a
    !> different set of sampling angles, and the shift is a phase multiply that needs no resampling
    !> at all -- so an entire shift grid is scored from one angular resample.
    !>
    !> `itest > 0` injects a deterministic per-particle perturbation first and reports what fraction
    !> of it the search recovers. On IgG-RL and Ribosembly the project poses ARE the ground truth, so
    !> perturb-and-recover is the only honest way to ask whether the score locates the right pose.
    !> ---------------------------------------------------------------------------------------------
    module subroutine polar_pose_refine_one( pg, fpl, Usall, Cfall, Cm0all, c00all, nsamp2, ncomp, nk, &
        &ndir, idir, iptcl, wr, prior, sig2, itest, rot_amp, sh_amp, rot_step0, sh_step0, l_guard, &
        &l_p2, nnn, dnn, rmat_pi, rmat_b, jdir, &
        &ca, sa, drot, dshx, dshy, xws, e0, e1, s0, s1, ealt, salt, npose, ntrue, nrej )
        type(polar_grid_t), intent(in)    :: pg
        type(fplane_type),  intent(in)    :: fpl
        integer,            intent(in)    :: nsamp2, ncomp, nk, ndir, idir, iptcl, itest
        real,               intent(in)    :: Usall(nsamp2,0:ncomp,ndir)
        real(dp),           intent(in)    :: Cfall(ncomp*ncomp,nk,ndir), Cm0all(ncomp,nk,ndir)
        real(dp),           intent(in)    :: c00all(nk,ndir)
        real,               intent(in)    :: wr(nk)
        real(dp),           intent(in)    :: prior(ncomp), sig2
        real,               intent(in)    :: rot_amp, sh_amp, rot_step0, sh_step0
        logical,            intent(in)    :: l_guard, l_p2
        integer,            intent(in)    :: nnn, dnn(nnn)
        real,               intent(in)    :: rmat_pi(3,3), rmat_b(3,3,ndir)
        integer,            intent(out)   :: jdir          !< the ACCEPTED bank direction
        real,               intent(inout) :: ca, sa
        real,               intent(out)   :: drot, dshx, dshy
        real,               intent(inout) :: xws(nsamp2)
        real(dp),           intent(inout) :: e0, e1, s0, s1, ealt, salt
        integer,            intent(inout) :: npose, ntrue, nrej
        integer, parameter :: NROT = 7, NSH = 7, NROUND = 3
        !> in-plane trials per candidate DIRECTION. Coarse on purpose: the direction stage only has
        !! to rank directions; the winner then gets the full multi-scale refinement.
        integer, parameter :: NROT_D = 3
        ! Coarsest step. Rotation and shift are NOT interchangeable units: a 1 degree in-plane
        ! rotation displaces the particle edge by mskrad*delta pixels, so for a 200 A particle at
        ! 6 A/px it is worth ~0.3 px of shift. The rotation grid therefore starts coarser.
        real,    parameter :: ROT_STEP0_DEF = 2.0, SH_STEP0_DEF = 1.0
        complex, allocatable :: xw0(:), xw1(:)
        real,    allocatable :: xtrial(:), bq(:)
        real(dp),allocatable :: Ach(:,:), uvec(:), cvec(:)
        real(dp) :: e_mm, aa, acon, best, sc, sc_start, sc_true, shc
        real     :: ca0, sa0, cad, sad, cdel, sdel, px, py, brot, bshx, bshy, rstep, sstep
        real     :: rot_in, shx_in, shy_in, cur_rot, cur_shx, cur_shy, brot0, bshx0, bshy0
        real     :: pr, px_in, py_in
        real(dp) :: h1s, h2s, h1f, h2f, hfull, best2, logdet, logdet2, bestd, e_mm2
        real(dp), allocatable :: Ach2(:,:), cvec2(:)
        real     :: drot2, dshx2, dshy2, caj, saj
        integer  :: q, r, ir, info, it, ish, jsh, iround, ic, jd, jdir0
        allocate(xw0(pg%nsamp), xw1(pg%nsamp), xtrial(nsamp2), bq(0:ncomp))
        allocate(Ach(ncomp,ncomp), uvec(ncomp), cvec(ncomp))
        allocate(Ach2(ncomp,ncomp), cvec2(ncomp))
        shc  = real(fpl%shconst(1), dp)
        acon = 1.d0                                        ! contrast; COV_UNIT_CONTRAST is the default
        aa   = acon*acon
        jdir = idir
        call polar_dir_tables(Cfall, Cm0all, c00all, ncomp, nk, ndir, jdir, wr, prior, sig2, &
            &acon, Ach, cvec, e_mm, logdet, info)
        if( info /= 0 )then
            drot = 0.; dshx = 0.; dshy = 0.
            deallocate(xw0, xw1, xtrial, bq, Ach, uvec, cvec, Ach2, cvec2)
            return
        endif
        ca0 = ca; sa0 = sa
        ! optional injected perturbation -- deterministic in the particle index so the test is
        ! reproducible and carries no RNG state
        rot_in = 0.; shx_in = 0.; shy_in = 0.
        if( itest > 0 )then
            rot_in = rot_amp*(2.0*pose_hash(iptcl, 1) - 1.0)
            shx_in = sh_amp *(2.0*pose_hash(iptcl, 2) - 1.0)
            shy_in = sh_amp *(2.0*pose_hash(iptcl, 3) - 1.0)
            cdel = cos(rot_in*PI/180.); sdel = sin(rot_in*PI/180.)
            ca0  = ca*cdel - sa*sdel
            sa0  = sa*cdel + ca*sdel
        endif
        ! The trial shift is expressed in UNPADDED lattice units, but shconst is defined against the
        ! PADDED grid gen_fplane4rec writes on, so the phase carries the pad factor. Missing it makes
        ! every trial shift half of what it claims to be.
        shc  = shc * real(OSMPL_PAD_FAC, dp)
        ! the search starts from the (possibly perturbed) incoming pose
        brot = 0.; bshx = shx_in; bshy = shy_in
        drot = brot; dshx = bshx; dshy = bshy
        brot0 = brot; bshx0 = bshx; bshy0 = bshy
        jdir0 = jdir
        ! ---- P2: OUT-OF-PLANE DIRECTION SEARCH over the neighbours of the current direction.
        !
        ! Each candidate has its own bank, hence its own G, Cholesky and log-determinant. The score
        ! is the full marginal INCLUDING log det A -- see polar_dir_tables for why that term cannot
        ! be dropped here even though P1 was free to ignore it. Candidates are scored with a coarse
        ! in-plane scan at the incoming shift; the winner then gets the full P1 refinement below.
        if( l_p2 )then
            bestd = huge(1.d0)
            do ic = 1, nnn
                jd = dnn(ic)
                call polar_dir_tables(Cfall, Cm0all, c00all, ncomp, nk, ndir, jd, wr, prior, sig2, &
                    &acon, Ach2, cvec2, e_mm2, logdet2, info)
                if( info /= 0 ) cycle
                ! the relative in-plane angle is different at every candidate direction
                call polar_relative_inplane(rmat_pi, rmat_b(:,:,jd), caj, saj)
                if( itest > 0 )then
                    cdel = cos(rot_in*PI/180.); sdel = sin(rot_in*PI/180.)
                    cad  = caj*cdel - saj*sdel
                    sad  = saj*cdel + caj*sdel
                    caj  = cad; saj = sad
                endif
                do it = 1, NROT_D
                    cur_rot = real(it - (NROT_D+1)/2)*rot_step0
                    cdel = cos(cur_rot*PI/180.); sdel = sin(cur_rot*PI/180.)
                    cad  = caj*cdel - saj*sdel
                    sad  = saj*cdel + caj*sdel
                    call polar_sample_at_pose(fpl%cmplx_plane, fpl%transfer_plane, pg, cad, sad, &
                        &-bshx*real(shc), -bshy*real(shc), xw0)
                    call polar_pose_score_halves(pg, xw0, Usall(1,0,jd), nsamp2, ncomp, Ach2, &
                        &cvec2, sig2, acon, e_mm2, xtrial, uvec, h1f, h2f, hfull)
                    sc = hfull + logdet2
                    if( sc < bestd )then
                        bestd = sc; jdir = jd; ca0 = caj; sa0 = saj; brot = cur_rot
                        Ach = Ach2; cvec = cvec2; e_mm = e_mm2; logdet = logdet2
                    endif
                end do
            end do
            drot = brot
        endif
        ! The trial shift is expressed in UNPADDED lattice units, but shconst is defined against the
        ! PADDED grid gen_fplane4rec writes on, so the phase carries the pad factor. Missing it makes
        ! every trial shift half of what it claims to be.
        ! score at the starting pose, for the "did the search find anything" diagnostic
        call polar_sample_at_pose(fpl%cmplx_plane, fpl%transfer_plane, pg, ca0, sa0, &
            &-bshx*real(shc), -bshy*real(shc), xw0)
        call polar_pose_score_halves(pg, xw0, Usall(1,0,jdir), nsamp2, ncomp, Ach, cvec, &
            &sig2, acon, e_mm, xtrial, uvec, h1s, h2s, hfull)
        sc_start = hfull
        best     = sc_start
        best2    = h2s
        drot2    = brot; dshx2 = bshx; dshy2 = bshy
        ! multi-scale coordinate descent: angle grid, then shift grid, halving the step each round
        rstep = rot_step0
        sstep = sh_step0
        do iround = 1, NROUND
            ! ---- angle
            do it = 1, NROT
                cur_rot = brot + real(it - (NROT+1)/2)*rstep
                cdel = cos(cur_rot*PI/180.); sdel = sin(cur_rot*PI/180.)
                cad  = ca0*cdel - sa0*sdel
                sad  = sa0*cdel + ca0*sdel
                call polar_sample_at_pose(fpl%cmplx_plane, fpl%transfer_plane, pg, cad, sad, &
                    &-bshx*real(shc), -bshy*real(shc), xw0)
                call polar_pose_score_halves(pg, xw0, Usall(1,0,jdir), nsamp2, ncomp, Ach, cvec, &
                    &sig2, acon, e_mm, xtrial, uvec, h1f, h2f, hfull)
                sc = hfull                             ! selection always uses the full plane
                if( sc < best )then
                    best = sc; drot = cur_rot; dshx = bshx; dshy = bshy
                endif
                if( h2f < best2 )then
                    best2 = h2f; drot2 = cur_rot; dshx2 = bshx; dshy2 = bshy
                endif
            end do
            brot = drot
            ! ---- shift, at the accepted angle: ONE resample, the whole grid from phase multiplies
            cdel = cos(brot*PI/180.); sdel = sin(brot*PI/180.)
            cad  = ca0*cdel - sa0*sdel
            sad  = sa0*cdel + ca0*sdel
            call polar_sample_at_pose(fpl%cmplx_plane, fpl%transfer_plane, pg, cad, sad, 0.0, 0.0, xw0)
            do ish = 1, NSH
                do jsh = 1, NSH
                    cur_shx = bshx + real(ish - (NSH+1)/2)*sstep
                    cur_shy = bshy + real(jsh - (NSH+1)/2)*sstep
                    px = -cur_shx*real(shc); py = -cur_shy*real(shc)
                    call polar_apply_shift(pg, cad, sad, px, py, xw0, xw1)
                    call polar_pose_score_halves(pg, xw1, Usall(1,0,jdir), nsamp2, ncomp, Ach, &
                        &cvec, sig2, acon, e_mm, xtrial, uvec, h1f, h2f, hfull)
                    sc = hfull
                    if( sc < best )then
                        best = sc; drot = brot; dshx = cur_shx; dshy = cur_shy
                    endif
                    if( h2f < best2 )then
                        best2 = h2f; drot2 = brot; dshx2 = cur_shx; dshy2 = cur_shy
                    endif
                end do
            end do
            bshx = dshx; bshy = dshy
            rstep = 0.5*rstep; sstep = 0.5*sstep
        end do
        ! ---- SPLIT-HALF ACCEPTANCE.
        !
        ! The pose is CHOSEN on the full plane -- restricting selection to one half costs more
        ! recovery than the extra independence buys (measured: 3.45 deg -> 1.69 rather than 1.36).
        ! It is then accepted only if BOTH interleaved halves independently score it better than the
        ! incoming pose. Measured on IgG-10k this halves the damage done to already-correct poses
        ! (0.71 -> 0.35 px) at no cost to recovery from real error (1.6 px -> 0.41 either way).
        !
        ! It does NOT fix the angular case, and cannot: the ~1.0 deg of residual movement is below
        ! the 1.2 deg identifiability floor, so no test built on this data can tell it from signal.
        if( l_guard .and. (abs(drot-brot0) > TINY .or. abs(dshx-bshx0) > TINY .or. &
            &abs(dshy-bshy0) > TINY .or. jdir /= jdir0) )then
            cdel = cos(drot*PI/180.); sdel = sin(drot*PI/180.)
            cad  = ca0*cdel - sa0*sdel
            sad  = sa0*cdel + ca0*sdel
            call polar_sample_at_pose(fpl%cmplx_plane, fpl%transfer_plane, pg, cad, sad, &
                &-dshx*real(shc), -dshy*real(shc), xw0)
            call polar_pose_score_halves(pg, xw0, Usall(1,0,jdir), nsamp2, ncomp, Ach, cvec, &
                &sig2, acon, e_mm, xtrial, uvec, h1f, h2f, hfull)
            if( .not. (h1f < h1s .and. h2f < h2s) )then
                drot = brot0; dshx = bshx0; dshy = bshy0     ! rejected: keep the incoming pose
                jdir = jdir0                                 ! ... including the direction
                nrej = nrej + 1
            endif
        endif
        ! In test mode, score the TRUE pose too. If truth does not beat what the search found, the
        ! score itself does not identify pose and no amount of search will fix that -- this is the
        ! number that separates "search too coarse" from "objective uninformative".
        sc_true = 0.d0
        if( itest > 0 )then
            cdel = cos(-rot_in*PI/180.); sdel = sin(-rot_in*PI/180.)
            cad  = ca0*cdel - sa0*sdel
            sad  = sa0*cdel + ca0*sdel
            call polar_sample_at_pose(fpl%cmplx_plane, fpl%transfer_plane, pg, cad, sad, 0.0, 0.0, xw0)
            sc_true = polar_pose_score(pg, xw0, Usall(1,0,jdir), nsamp2, ncomp, Ach, cvec, prior, &
                &sig2, acon, e_mm, xtrial, bq, uvec)
        endif
        ! adopt the accepted pose: resample there and hand the samples back
        cdel = cos(drot*PI/180.); sdel = sin(drot*PI/180.)
        cad  = ca0*cdel - sa0*sdel
        sad  = sa0*cdel + ca0*sdel
        px   = -dshx*real(shc); py = -dshy*real(shc)
        call polar_sample_at_pose(fpl%cmplx_plane, fpl%transfer_plane, pg, cad, sad, px, py, xw0)
        do q = 1, pg%nsamp
            xws(2*q-1) = pg%sqwq(q)*real(xw0(q))
            xws(2*q)   = pg%sqwq(q)*aimag(xw0(q))
        end do
        ca = cad; sa = sad
        ! statistics
        npose = npose + 1
        if( itest > 0 )then
            e0 = e0 + real(rot_in,dp)**2
            e1 = e1 + real(rot_in + drot, dp)**2
            s0 = s0 + real(shx_in,dp)**2 + real(shy_in,dp)**2
            s1 = s1 + real(dshx,dp)**2 + real(dshy,dp)**2
            ! ntrue counts particles where the TRUE pose scores better than the search's answer
            if( sc_true < best ) ntrue = ntrue + 1
        else if( COV_PROJ_PERTURB_ROT > 0. .or. COV_PROJ_PERTURB_SH > 0. )then
            ! the project's poses were degraded by a KNOWN amount; report what is left of it. The
            ! sign of the in-plane convention is not assumed -- both are accumulated and the log
            ! prints whichever the search actually produced.
            pr = COV_PROJ_PERTURB_ROT*(2.0*pose_hash(iptcl, 1) - 1.0)
            px_in = COV_PROJ_PERTURB_SH*(2.0*pose_hash(iptcl, 2) - 1.0)
            py_in = COV_PROJ_PERTURB_SH*(2.0*pose_hash(iptcl, 3) - 1.0)
            e0 = e0 + real(pr,dp)**2
            e1 = e1 + real(pr + drot, dp)**2
            s0 = s0 + real(px_in,dp)**2 + real(py_in,dp)**2
            s1 = s1 + real(px_in + dshx, dp)**2 + real(py_in + dshy, dp)**2
            ealt = ealt + real(pr - drot, dp)**2
            salt = salt + real(px_in - dshx, dp)**2 + real(py_in - dshy, dp)**2
            ntrue = ntrue + 1
        else
            e1 = e1 + real(drot,dp)**2
            s1 = s1 + real(dshx,dp)**2 + real(dshy,dp)**2
            if( best < sc_start ) ntrue = ntrue + 1
        endif
        deallocate(xw0, xw1, xtrial, bq, Ach, uvec, cvec, Ach2, cvec2)
    end subroutine polar_pose_refine_one

    !> Everything in the pose objective that depends on the candidate DIRECTION: the Gram, its
    !! Cholesky, the log-determinant, the mean cross term and <T mu, T mu>.
    !!
    !! For P1 this was computed once per particle because `G` is invariant under in-plane rotation
    !! and shift. It is NOT invariant across directions, so P2 pays this per candidate -- and, more
    !! importantly, must KEEP the log-det term. It cancels between in-plane trials and does not
    !! cancel between directions; dropping it biases the search toward directions where the basis
    !! happens to be more expressive, which is the term a naive implementation silently loses.
    module subroutine polar_dir_tables( Cfall, Cm0all, c00all, ncomp, nk, ndir, jd, wr, prior, sig2, &
        &acon, Ach, cvec, e_mm, logdet, info )
        integer,  intent(in)  :: ncomp, nk, ndir, jd
        real(dp), intent(in)  :: Cfall(ncomp*ncomp,nk,ndir), Cm0all(ncomp,nk,ndir), c00all(nk,ndir)
        real,     intent(in)  :: wr(nk)
        real(dp), intent(in)  :: prior(ncomp), sig2, acon
        real(dp), intent(out) :: Ach(ncomp,ncomp), cvec(ncomp), e_mm, logdet
        integer,  intent(out) :: info
        real(dp) :: aa, w
        integer  :: q, r, ir
        aa   = acon*acon
        e_mm = 0.d0
        cvec = 0.d0
        Ach  = 0.d0
        do ir = 1, nk
            w    = real(wr(ir),dp)
            e_mm = e_mm + w*c00all(ir,jd)
            do q = 1, ncomp
                cvec(q) = cvec(q) + w*Cm0all(q,ir,jd)
                do r = 1, ncomp
                    Ach(q,r) = Ach(q,r) + w*Cfall((r-1)*ncomp+q,ir,jd)
                end do
            end do
        end do
        Ach = (aa/sig2)*Ach
        do q = 1, ncomp
            Ach(q,q) = Ach(q,q) + prior(q)
        end do
        call dpotrf('U', ncomp, Ach, ncomp, info)
        logdet = 0.d0
        if( info == 0 )then
            do q = 1, ncomp
                logdet = logdet + 2.d0*log(max(Ach(q,q), DTINY))
            end do
        endif
    end subroutine polar_dir_tables

    !> The same objective evaluated on each interleaved half of the sample list. Used ONLY to decide
    !! whether to accept a move, never to choose it: a pose picked as the argmax over ~170 trials is
    !! biased to fit noise, and the measured per-particle floor (~0.4 px of displacement at SNR 0.01)
    !! is large enough that refining an already-correct pose makes it WORSE. Requiring both halves to
    !! independently prefer the move is a threshold-free guard against exactly that.
    module subroutine polar_pose_score_halves( pg, xw, Us, nsamp2, ncomp, Ach, cvec, sig2, acon, e_mm, &
        &xtrial, uvec, sc1, sc2, scfull )
        type(polar_grid_t), intent(in)    :: pg
        complex,            intent(in)    :: xw(:)
        integer,            intent(in)    :: nsamp2, ncomp
        real,               intent(in)    :: Us(nsamp2,0:ncomp)
        real(dp),           intent(in)    :: Ach(ncomp,ncomp), cvec(ncomp), sig2, acon, e_mm
        real,               intent(inout) :: xtrial(nsamp2)
        real(dp),           intent(inout) :: uvec(ncomp)
        real(dp),           intent(out)   :: sc1, sc2, scfull
        real(dp) :: b1(0:ncomp), b2(0:ncomp), aa, quad
        integer  :: j, q
        do j = 1, pg%nsamp
            xtrial(2*j-1) = pg%sqwq(j)*real(xw(j))
            xtrial(2*j)   = pg%sqwq(j)*aimag(xw(j))
        end do
        b1 = 0.d0; b2 = 0.d0
        do q = 0, ncomp
            do j = 1, nsamp2
                if( pg%hrow(j) == 1 )then
                    b1(q) = b1(q) + real(Us(j,q),dp)*real(xtrial(j),dp)
                else
                    b2(q) = b2(q) + real(Us(j,q),dp)*real(xtrial(j),dp)
                endif
            end do
        end do
        aa = acon*acon
        do q = 1, ncomp
            uvec(q) = (acon*b1(q) - 0.5d0*aa*cvec(q))/sig2
        end do
        call dtrsv('U', 'T', 'N', ncomp, Ach, ncomp, uvec, 1)
        quad = 0.d0
        do q = 1, ncomp
            quad = quad + uvec(q)*uvec(q)
        end do
        sc1 = (0.5d0*aa*e_mm - 2.d0*acon*b1(0))/sig2 - quad
        do q = 1, ncomp
            uvec(q) = (acon*b2(q) - 0.5d0*aa*cvec(q))/sig2
        end do
        call dtrsv('U', 'T', 'N', ncomp, Ach, ncomp, uvec, 1)
        quad = 0.d0
        do q = 1, ncomp
            quad = quad + uvec(q)*uvec(q)
        end do
        sc2 = (0.5d0*aa*e_mm - 2.d0*acon*b2(0))/sig2 - quad
        ! The FULL-plane marginal is NOT sc1+sc2: the quadratic form is (u1+u2)^T A^-1 (u1+u2), and
        ! dropping the cross term halves it relative to the linear term. Selection must use this one;
        ! sc1/sc2 exist only for the split-half acceptance test.
        do q = 1, ncomp
            uvec(q) = (acon*(b1(q) + b2(q)) - aa*cvec(q))/sig2
        end do
        call dtrsv('U', 'T', 'N', ncomp, Ach, ncomp, uvec, 1)
        quad = 0.d0
        do q = 1, ncomp
            quad = quad + uvec(q)*uvec(q)
        end do
        scfull = (aa*e_mm - 2.d0*acon*(b1(0) + b2(0)))/sig2 - quad
    end subroutine polar_pose_score_halves

    !> -2 log p(y | pose) up to terms constant across trials of the same particle.
    real(dp) module function polar_pose_score( pg, xw, Us, nsamp2, ncomp, Ach, cvec, prior, sig2, acon, &
        &e_mm, xtrial, bq, uvec ) result( sc )
        type(polar_grid_t), intent(in)    :: pg
        complex,            intent(in)    :: xw(:)
        integer,            intent(in)    :: nsamp2, ncomp
        real,               intent(in)    :: Us(nsamp2,0:ncomp)
        real(dp),           intent(in)    :: Ach(ncomp,ncomp), cvec(ncomp), prior(ncomp), sig2
        real(dp),           intent(in)    :: acon, e_mm
        real,               intent(inout) :: xtrial(nsamp2), bq(0:ncomp)
        real(dp),           intent(inout) :: uvec(ncomp)
        integer  :: j, q
        real(dp) :: aa, myv, quad
        do j = 1, pg%nsamp
            xtrial(2*j-1) = pg%sqwq(j)*real(xw(j))
            xtrial(2*j)   = pg%sqwq(j)*aimag(xw(j))
        end do
        call sgemv('T', nsamp2, ncomp+1, 1.0, Us, nsamp2, xtrial, 1, 0.0, bq, 1)
        aa  = acon*acon
        myv = real(bq(0), dp)
        do q = 1, ncomp
            uvec(q) = (acon*real(bq(q),dp) - aa*cvec(q))/sig2
        end do
        call dtrsv('U', 'T', 'N', ncomp, Ach, ncomp, uvec, 1)
        quad = 0.d0
        do q = 1, ncomp
            quad = quad + uvec(q)*uvec(q)
        end do
        sc = (aa*e_mm - 2.d0*acon*myv)/sig2 - quad
    end function polar_pose_score

    !> ---------------------------------------------------------------------------------------------
    !> Degrade the PROJECT's in-plane angles and shifts before anything reads them.
    !>
    !> SIMPLE_COV_POSE_TEST perturbs only inside the embedding, which measures whether the score can
    !> find a pose but leaves the model itself built from correct poses. That is the easy half of the
    !> question. This perturbs `build%spproj_field` at the very start of the run, so the mean, the
    !> covariance columns, the basis AND the embedding are all built from degraded poses -- the
    !> feedback loop that a real ab-initio run has and the embed-only test does not.
    !>
    !> The perturbation is hashed from the particle index (same hash as the embed-only test), so it
    !> is exactly reproducible and its magnitude is known, which is what makes the delivered basis
    !> comparable against an unperturbed run.
    !> ---------------------------------------------------------------------------------------------
    module subroutine cov_perturb_project_poses( build, pinds, nptcls )
        type(builder), intent(inout) :: build
        integer,       intent(in)    :: pinds(:), nptcls
        type(ori) :: o
        integer :: irot, ish, idir, i, j, iptcl
        real    :: ramp, samp, e3, dr, dx, dy, sh(2), dd, psi, cps, sps, cdd, sdd
        real    :: rm(3,3), rm2(3,3), ax(3), v(3), cr(3)
        real(dp):: acc_r, acc_s, acc_d
        irot = 0; call cov_env_int('SIMPLE_COV_POSE_PERTURB_ROT', irot)
        ish  = 0; call cov_env_int('SIMPLE_COV_POSE_PERTURB_SH',  ish)
        idir = 0; call cov_env_int('SIMPLE_COV_POSE_PERTURB_DIR', idir)
        if( irot <= 0 .and. ish <= 0 .and. idir <= 0 ) return
        COV_PROJ_PERTURB_ROT = 0.1*real(irot)
        COV_PROJ_PERTURB_SH  = 0.1*real(ish)
        COV_PROJ_PERTURB_DIR = 0.1*real(idir)
        ramp = 0.1*real(irot)
        samp = 0.1*real(ish)
        acc_r = 0.d0; acc_s = 0.d0; acc_d = 0.d0
        if( COV_PROJ_PERTURB_DIR > 0. )then
            if( allocated(COV_TRUTH_NRM) ) deallocate(COV_TRUTH_NRM)
            allocate(COV_TRUTH_NRM(3,nptcls))
        endif
        do i = 1, nptcls
            iptcl = pinds(i)
            dr = ramp*(2.0*pose_hash(i, 1) - 1.0)
            dx = samp*(2.0*pose_hash(i, 2) - 1.0)
            dy = samp*(2.0*pose_hash(i, 3) - 1.0)
            if( COV_PROJ_PERTURB_DIR > 0. )then
                ! tilt the VIEWING DIRECTION by dd about an axis lying in the projection plane, so
                ! the normal moves by exactly dd. Rodrigues on each row of the frame.
                call build%spproj_field%get_ori(iptcl, o)
                rm = o%get_mat()
                COV_TRUTH_NRM(:,i) = rm(3,:)
                dd  = COV_PROJ_PERTURB_DIR*(2.0*pose_hash(i, 4) - 1.0)
                psi = 360.0*pose_hash(i, 5)
                cps = cos(psi*PI/180.); sps = sin(psi*PI/180.)
                ax  = cps*rm(1,:) + sps*rm(2,:)                    ! axis in the plane
                cdd = cos(dd*PI/180.); sdd = sin(dd*PI/180.)
                do j = 1, 3
                    v          = rm(j,:)
                    cr(1)      = ax(2)*v(3) - ax(3)*v(2)
                    cr(2)      = ax(3)*v(1) - ax(1)*v(3)
                    cr(3)      = ax(1)*v(2) - ax(2)*v(1)
                    rm2(j,:)   = v*cdd + cr*sdd + ax*dot_product(ax,v)*(1.0 - cdd)
                end do
                call o%ori_from_rotmat(rm2, .true.)
                call build%spproj_field%set_euler(iptcl, o%get_euler())
                acc_d = acc_d + real(dd,dp)**2
                call o%kill
            endif
            e3 = build%spproj_field%e3get(iptcl)
            call build%spproj_field%e3set(iptcl, e3 + dr)
            sh = build%spproj_field%get_2Dshift(iptcl)
            call build%spproj_field%set_shift(iptcl, sh + [dx, dy])
            acc_r = acc_r + real(dr,dp)**2
            acc_s = acc_s + real(dx,dp)**2 + real(dy,dp)**2
        end do
        write(logfhandle,'(A,F8.4,A,F8.4,A,F8.4)') '>>> FLEX_PCA PROJECT POSES DEGRADED: in-plane rms=', &
            &sqrt(acc_r/real(nptcls,dp)),' deg   shift rms=',sqrt(acc_s/real(nptcls,dp)), &
            &' px   direction rms=',sqrt(acc_d/real(nptcls,dp))
        write(logfhandle,'(A)') '>>> FLEX_PCA   (every stage below now sees the degraded poses)'
        call flush(logfhandle)
    end subroutine cov_perturb_project_poses

    !> Deterministic per-particle uniform-ish value in [0,1) for the perturbation test. A hash, not
    !! an RNG, so the injected perturbation is identical across runs and across thread counts.
    pure real module function pose_hash( i, k ) result( h )
        integer, intent(in) :: i, k
        integer(kind=8) :: x
        x = int(i, 8)*2654435761_8 + int(k, 8)*40503_8
        x = iand(ieor(x, ishft(x, -13)), 2147483647_8)
        h = real(modulo(x, 100000_8)) / 100000.0
    end function pose_hash

    !> One ring's contribution to the Gram of the basis (columns 1..ncomp of Us) and to the mean
    !! cross term (column 0). `Cout` is the full ncomp x ncomp block flattened column-major, `Mout`
    !! the ncomp-vector <U_q, T mu>.
    module subroutine polar_ring_gram( Us, ldu, ncomp, row0, nrow, Csp, Cout, Mout )
        integer,  intent(in)    :: ldu, ncomp, row0, nrow
        real,     intent(in)    :: Us(ldu,0:ncomp)
        real,     intent(inout) :: Csp(0:ncomp,0:ncomp)      !< caller-owned scratch
        real(dp), intent(out)   :: Cout(ncomp*ncomp), Mout(ncomp)
        integer :: q, r, i0, n2
        i0 = 2*row0 - 1
        n2 = 2*nrow
        if( n2 <= 0 )then
            Cout = 0.d0; Mout = 0.d0
            return
        endif
        call ssyrk('U','T', ncomp+1, n2, 1.0, Us(i0,0), ldu, 0.0, Csp, ncomp+1)
        do r = 1, ncomp
            do q = 1, r
                Cout((r-1)*ncomp+q) = real(Csp(q,r), dp)
                Cout((q-1)*ncomp+r) = real(Csp(q,r), dp)
            end do
            Mout(r) = real(Csp(0,r), dp)
        end do
    end subroutine polar_ring_gram

    real(dp) module function polar_ring_selfpower( Us, ldu, row0, nrow )
        integer, intent(in) :: ldu, row0, nrow
        real,    intent(in) :: Us(ldu,0:*)
        integer :: j
        polar_ring_selfpower = 0.d0
        do j = 2*row0-1, 2*(row0+nrow-1)
            polar_ring_selfpower = polar_ring_selfpower + real(Us(j,0),dp)*real(Us(j,0),dp)
        end do
    end function polar_ring_selfpower

    !> <y,y> in the polar measure. xws already carries sqrt(wq) and the CTF adjoint, so the CTF has
    !! to be divided back out ring by ring to recover the observation's own energy.
    real(dp) module function polar_self_energy( xws, wr, pg ) result( e )
        real,               intent(in) :: xws(:), wr(:)
        type(polar_grid_t), intent(in) :: pg
        integer  :: ir, j
        real(dp) :: acc
        e = 0.d0
        do ir = 1, pg%nk
            if( real(wr(ir),dp) <= DTINY ) cycle
            acc = 0.d0
            do j = 2*pg%rbeg(ir)-1, 2*pg%rend(ir)
                acc = acc + real(xws(j),dp)*real(xws(j),dp)
            end do
            e = e + acc / real(wr(ir),dp)
        end do
    end function polar_self_energy

    pure real(dp) module function sum_dp_safe( acc, n ) result( v )
        real(dp), intent(in) :: acc
        integer,  intent(in) :: n
        v = acc / real(max(1,n), dp)
    end function sum_dp_safe

    !> thin wrapper so the OpenMP body stays readable; folds sqrt(wq) into the stored samples
    !> resample a stored half-plane by an in-plane rotation into the bank frame: unit-tap 2D KB
    !! at pf-multiples with per-tap Friedel, per-axis normalized weights, OOB taps dropped --
    !! the polar former's interpolation scheme (polar_interp_plane) on the Cartesian lattice.
    !! Positions outside the nyq disk come back ZERO (nothing downstream reads them).
    module subroutine align_halfplane_inplane( frlims, nyq_eff, src, ca, sa, dst )
        integer, intent(in)  :: frlims(3,2), nyq_eff
        complex, intent(in)  :: src(frlims(1,1):frlims(1,2), frlims(2,1):0)
        real,    intent(in)  :: ca, sa
        complex, intent(out) :: dst(frlims(1,1):frlims(1,2), frlims(2,1):0)
        type(kbinterpol) :: kbwin
        real    :: hu, ku, w, wx(LATENT_WDIM), wy(LATENT_WDIM), wz(LATENT_WDIM)
        integer :: win(2,3), h, k, hlo2, hhi2, klo2, hx, ky, ix, iy, nd, pf
        complex :: acc, cv
        kbwin = kbinterpol(KBWINSZ, KBALPHA)
        pf    = OSMPL_PAD_FAC
        hlo2  = ceil_div (frlims(1,1), pf); hhi2 = floor_div(frlims(1,2), pf)
        klo2  = ceil_div (frlims(2,1), pf)
        nd    = nyq_eff*(nyq_eff+1)
        dst   = CMPLX_ZERO
        do k = klo2, 0
            do h = hlo2, hhi2
                if( h*h + k*k > nd ) cycle
                hu =  h*ca + k*sa
                ku = -h*sa + k*ca
                call latent_projection_weights(kbwin, [hu, ku, 0.], win, wx, wy, wz)
                acc = CMPLX_ZERO
                do iy = 1, LATENT_WDIM
                    ky = win(1,2) + iy - 1
                    do ix = 1, LATENT_WDIM
                        hx = win(1,1) + ix - 1
                        w  = wx(ix)*wy(iy)
                        if( pf*ky <= 0 )then
                            if( pf*hx < frlims(1,1) .or. pf*hx > frlims(1,2) .or. pf*ky < frlims(2,1) ) cycle
                            cv = src(pf*hx, pf*ky)
                        else
                            if( -pf*hx < frlims(1,1) .or. -pf*hx > frlims(1,2) .or. -pf*ky < frlims(2,1) ) cycle
                            cv = conjg(src(-pf*hx, -pf*ky))
                        endif
                        acc = acc + w*cv
                    end do
                end do
                dst(pf*h, pf*k) = acc
            end do
        end do
    end subroutine align_halfplane_inplane

    module subroutine polar_sample_particle_packed( fpl, pg, ca, sa, xws, wr, hfpw, hfcnt, tazim, xws1, xws2 )
        type(fplane_type),  intent(in)    :: fpl
        type(polar_grid_t), intent(in)    :: pg
        real,               intent(in)    :: ca, sa
        real,               intent(out)   :: xws(:)
        real,               intent(out)   :: wr(:)
        real(dp),           intent(inout) :: hfpw, hfcnt
        real,               intent(out)   :: tazim
        !> packed forms of the two lattice-parity half-fields, for the reliability prior
        real, optional,     intent(out)   :: xws1(:), xws2(:)
        complex, allocatable :: xw(:), xw1(:), xw2(:)
        real(dp) :: pw, cnt
        integer  :: j
        allocate(xw(pg%nsamp), xw1(pg%nsamp), xw2(pg%nsamp))
        if( present(xws1) )then
            call polar_sample_particle(fpl%cmplx_plane, fpl%transfer_plane, pg, ca, sa, xw, wr, &
                &pw, cnt, tazim, xw1, xw2)
        else
            call polar_sample_particle(fpl%cmplx_plane, fpl%transfer_plane, pg, ca, sa, xw, wr, &
                &pw, cnt, tazim)
        endif
        do j = 1, pg%nsamp
            xws(2*j-1) = pg%sqwq(j)*real(xw(j))
            xws(2*j)   = pg%sqwq(j)*aimag(xw(j))
        end do
        if( present(xws1) )then
            do j = 1, pg%nsamp
                xws1(2*j-1) = pg%sqwq(j)*real(xw1(j))
                xws1(2*j)   = pg%sqwq(j)*aimag(xw1(j))
                xws2(2*j-1) = pg%sqwq(j)*real(xw2(j))
                xws2(2*j)   = pg%sqwq(j)*aimag(xw2(j))
            end do
        endif
        hfpw  = hfpw  + pw
        hfcnt = hfcnt + cnt
        deallocate(xw, xw1, xw2)
    end subroutine polar_sample_particle_packed

    !> Banded mean projection for the polar E-step: reconstructor%project_fplane's numerics --
    !! the SAME banded (h,k) sweep, apod_mat_3d interpolation weights (including their final
    !! global renormalization), per-sample Friedel conjugation and transfer multiply -- with the
    !! per-call full-plane work removed. project_fplane zero-fills the whole PADDED plane and
    !! copies the reference ctfsq and transfer planes into the output on EVERY call; at the
    !! native padded lattice that is several MB of memory traffic per particle, which measured
    !! as ~80% of the polar E-step's project bucket, all spent on values the polar branch never
    !! reads (only mean_fpl%cmplx_plane is consumed, by the residual subtraction). Here the
    !! plane is zero-filled once at (re)allocation; every call rewrites exactly the in-band disc
    !! samples, and out-of-disc positions stay zero -- the same invariant the Cartesian former's
    !! ensure_latent_projection_plane relies on. The interpolated values are bit-identical to
    !! project_fplane's (same expressions, same kbwin), so the residual planes the M-step
    !! consumes are unchanged wherever the mean is nonzero and unchanged-because-zero elsewhere.
    module subroutine project_fplane_mean_banded( rec, o, fpl_ref, fpl_out )
        type(reconstructor), intent(in)    :: rec
        class(ori),          intent(inout) :: o
        type(fplane_type),   intent(in)    :: fpl_ref
        type(fplane_type),   intent(inout) :: fpl_out
        type(kbinterpol) :: kbwin
        real    :: rotmat(3,3), loc(3), loc_friedel(3), hrow(3)
        real    :: w3(LATENT_WDIM,LATENT_WDIM,LATENT_WDIM)
        integer :: fpllims_pd(3,2), fpllims(3,2), h, k, hp, kp, pf, iwinsz, win(2,3)
        integer :: h_sq, k_max_h, k_lo, k_hi, nyq_disk, nyq_eff
        logical :: l_conjg, l_realloc
        complex :: comp
        kbwin  = kbinterpol(KBWINSZ, KBALPHA)
        iwinsz = ceiling(KBWINSZ - 0.5)
        fpl_out%frlims  = fpl_ref%frlims
        fpl_out%shconst = fpl_ref%shconst
        fpl_out%nyq     = fpl_ref%nyq
        l_realloc = .not. allocated(fpl_out%cmplx_plane)
        if( .not. l_realloc )then
            l_realloc = any(lbound(fpl_out%cmplx_plane) /= lbound(fpl_ref%cmplx_plane)) .or. &
                &any(ubound(fpl_out%cmplx_plane) /= ubound(fpl_ref%cmplx_plane))
        endif
        if( l_realloc )then
            if( allocated(fpl_out%cmplx_plane) ) deallocate(fpl_out%cmplx_plane)
            allocate(fpl_out%cmplx_plane(lbound(fpl_ref%cmplx_plane,1):ubound(fpl_ref%cmplx_plane,1), &
                &lbound(fpl_ref%cmplx_plane,2):ubound(fpl_ref%cmplx_plane,2)))
            fpl_out%cmplx_plane = CMPLX_ZERO
        endif
        rotmat      = o%get_mat()
        pf          = OSMPL_PAD_FAC
        fpllims_pd  = fpl_ref%frlims
        fpllims     = fpllims_pd
        fpllims(1,1)= ceil_div (fpllims_pd(1,1), pf)
        fpllims(1,2)= floor_div(fpllims_pd(1,2), pf)
        fpllims(2,1)= ceil_div (fpllims_pd(2,1), pf)
        fpllims(2,2)= floor_div(fpllims_pd(2,2), pf)
        nyq_eff = rec%get_lfny(1)
        if( fpl_ref%nyq > 0 ) nyq_eff = min(nyq_eff, max(1, fpl_ref%nyq / pf))
        nyq_disk = nyq_eff * (nyq_eff + 1)
        do h = fpllims(1,1), fpllims(1,2)
            h_sq = h*h
            if( h_sq > nyq_disk ) cycle
            k_max_h = int(sqrt(real(nyq_disk - h_sq)))
            k_lo    = max(fpllims(2,1), -k_max_h)
            k_hi    = min(0, min(fpllims(2,2), k_max_h))
            hp      = h * pf
            hrow(1) = real(h) * rotmat(1,1)
            hrow(2) = real(h) * rotmat(1,2)
            hrow(3) = real(h) * rotmat(1,3)
            do k = k_lo, k_hi
                kp     = k * pf
                loc(1) = hrow(1) + real(k) * rotmat(2,1)
                loc(2) = hrow(2) + real(k) * rotmat(2,2)
                loc(3) = hrow(3) + real(k) * rotmat(2,3)
                ! interp_cmat_exp, verbatim (it is private to simple_reconstructor)
                l_conjg     = loc(1) < 0.
                loc_friedel = loc
                if( l_conjg ) loc_friedel = -loc_friedel
                win(1,:) = nint(loc_friedel)
                win(2,:) = win(1,:) + iwinsz
                win(1,:) = win(1,:) - iwinsz
                call kbwin%apod_mat_3d(loc_friedel, iwinsz, LATENT_WDIM, w3)
                comp = sum(w3 * rec%cmat_exp(win(1,1):win(2,1), win(1,2):win(2,2), win(1,3):win(2,3)))
                if( l_conjg ) comp = conjg(comp)
                ! apply_ctf_amp=.true. semantics of project_fplane
                if( allocated(fpl_ref%transfer_plane) )then
                    fpl_out%cmplx_plane(hp,kp) = fpl_ref%transfer_plane(hp,kp) * comp
                else
                    fpl_out%cmplx_plane(hp,kp) = sqrt(max(0., fpl_ref%ctfsq_plane(hp,kp))) * comp
                endif
            end do
        end do
    end subroutine project_fplane_mean_banded

    !> Exact Cartesian statistics of the low-k shells for the HYBRID polar E-step, added on top
    !! of the ring statistics. Per lattice position the KB window geometry is computed once and
    !! all ncomp+1 volumes are gathered through it (the Cartesian former's hoist); the data value,
    !! CTF/whitening transfer and quadrature weight (1 per lattice point) are exactly the
    !! Cartesian former's, so the shells this covers contribute to G/b/c/e_mm/myv precisely what
    !! project_fplanes_mean_basis + cov_herm_inner would contribute for them -- including the DC
    !! sample. This is what removes the ring quadrature's multiplicative posterior-variance bias:
    !! after whitening the low-k shells still anchor the latent scale, and rings sample them
    !! worst (few samples, steep integrand).
    module subroutine polar_hybrid_exact_accum( rec0, recs, ncomp, o, fpl, hex, kex, npos, &
            &Gd, bd, cd, e_mm, myv )
        type(reconstructor), intent(in)    :: rec0
        type(reconstructor), intent(in)    :: recs(ncomp)
        integer,             intent(in)    :: ncomp, npos
        class(ori),          intent(inout) :: o
        type(fplane_type),   intent(in)    :: fpl
        integer,             intent(in)    :: hex(npos), kex(npos)
        real(dp),            intent(inout) :: Gd(ncomp,ncomp), bd(ncomp), cd(ncomp)
        real(dp),            intent(inout) :: e_mm, myv
        type(kbinterpol) :: kbwin
        real        :: rotmat(3,3), loc(3), wx(LATENT_WDIM), wy(LATENT_WDIM), wz(LATENT_WDIM)
        integer     :: j, q, r, win(2,3), hp, kp, exp_lb(3), exp_ub(3), pf
        logical     :: l_conjg, l_tf
        complex     :: tf, yv, u0, val
        complex     :: uq(ncomp)
        complex(dp) :: u0d, yd
        kbwin  = kbinterpol(KBWINSZ, KBALPHA)
        rotmat = o%get_mat()
        pf     = OSMPL_PAD_FAC
        exp_lb = lbound(rec0%cmat_exp)
        exp_ub = ubound(rec0%cmat_exp)
        l_tf   = allocated(fpl%transfer_plane)
        do j = 1, npos
            loc(1) = real(hex(j))*rotmat(1,1) + real(kex(j))*rotmat(2,1)
            loc(2) = real(hex(j))*rotmat(1,2) + real(kex(j))*rotmat(2,2)
            loc(3) = real(hex(j))*rotmat(1,3) + real(kex(j))*rotmat(2,3)
            l_conjg = loc(1) < 0.
            if( l_conjg ) loc = -loc
            call latent_projection_weights(kbwin, loc, win, wx, wy, wz)
            if( any(win(1,:) < exp_lb) .or. any(win(2,:) > exp_ub) ) cycle
            hp = pf*hex(j)
            kp = pf*kex(j)
            if( l_tf )then
                tf = fpl%transfer_plane(hp,kp)
            else
                tf = cmplx(sqrt(max(0., fpl%ctfsq_plane(hp,kp))), 0.)
            endif
            yv  = fpl%cmplx_plane(hp,kp)
            val = weighted_expanded_cmat(rec0, win, wx, wy, wz)
            if( l_conjg ) val = conjg(val)
            u0 = tf * val
            do q = 1, ncomp
                val = weighted_expanded_cmat(recs(q), win, wx, wy, wz)
                if( l_conjg ) val = conjg(val)
                uq(q) = tf * val
            end do
            u0d  = cmplx(u0, kind=dp)
            yd   = cmplx(yv, kind=dp)
            e_mm = e_mm + real(conjg(u0d)*u0d, dp)
            myv  = myv  + real(conjg(u0d)*yd,  dp)
            do q = 1, ncomp
                bd(q) = bd(q) + real(conjg(cmplx(uq(q),kind=dp))*yd,  dp)
                cd(q) = cd(q) + real(conjg(cmplx(uq(q),kind=dp))*u0d, dp)
                do r = q, ncomp
                    Gd(q,r) = Gd(q,r) + real(conjg(cmplx(uq(q),kind=dp))*cmplx(uq(r),kind=dp), dp)
                end do
            end do
        end do
        ! mirror the accumulated upper triangle (the ring dgemv filled both triangles already;
        ! the exact increments above touched q<=r only)
        do r = 1, ncomp
            do q = r+1, ncomp
                Gd(q,r) = Gd(r,q)
            end do
        end do
    end subroutine polar_hybrid_exact_accum

    !> Banded residual subtraction, fpl = fpl - a*mean over EXACTLY the disc the banded (or any
    !! full-plane) mean projection wrote. Everywhere outside that disc the mean plane is
    !! identically zero, so the full-array statement this replaces only rewrote unchanged values
    !! there -- another few MB of per-particle traffic at the native padded lattice for no effect.
    !! The loop bounds are the same expressions as project_fplane_mean_banded's, so written and
    !! subtracted sample sets coincide by construction.
    module subroutine subtract_mean_banded( fpl, mean_fpl, a, rec_nyq )
        type(fplane_type), intent(inout) :: fpl
        type(fplane_type), intent(in)    :: mean_fpl
        real,              intent(in)    :: a
        integer,           intent(in)    :: rec_nyq
        integer :: fpllims_pd(3,2), fpllims(3,2), h, k, hp, kp, pf
        integer :: h_sq, k_max_h, k_lo, k_hi, nyq_disk, nyq_eff
        pf          = OSMPL_PAD_FAC
        fpllims_pd  = fpl%frlims
        fpllims     = fpllims_pd
        fpllims(1,1)= ceil_div (fpllims_pd(1,1), pf)
        fpllims(1,2)= floor_div(fpllims_pd(1,2), pf)
        fpllims(2,1)= ceil_div (fpllims_pd(2,1), pf)
        fpllims(2,2)= floor_div(fpllims_pd(2,2), pf)
        nyq_eff = rec_nyq
        if( fpl%nyq > 0 ) nyq_eff = min(nyq_eff, max(1, fpl%nyq / pf))
        nyq_disk = nyq_eff * (nyq_eff + 1)
        do h = fpllims(1,1), fpllims(1,2)
            h_sq = h*h
            if( h_sq > nyq_disk ) cycle
            k_max_h = int(sqrt(real(nyq_disk - h_sq)))
            k_lo    = max(fpllims(2,1), -k_max_h)
            k_hi    = min(0, min(fpllims(2,2), k_max_h))
            hp      = h * pf
            do k = k_lo, k_hi
                kp = k * pf
                fpl%cmplx_plane(hp,kp) = fpl%cmplx_plane(hp,kp) - a*mean_fpl%cmplx_plane(hp,kp)
            end do
        end do
    end subroutine subtract_mean_banded

end submodule simple_flex_pca_em_pose
