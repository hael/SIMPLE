!@descr: kernel-weighted state reconstruction for flex_pca
module simple_flex_pca_rec3D
use simple_core_module_api
use simple_builder,                only: builder
use simple_sp_project,             only: sp_project
use simple_gridding,               only: prep3D_inv_kbenvelope4mul
use simple_image,                  only: image
use simple_matcher_3Drec,          only: init_rec, prep_imgs4rec, cleanup_rec_buffers
use simple_matcher_ptcl_io,        only: discrete_read_imgbatch, discrete_read_imgbatch_source, prepimgbatch
use simple_parameters,             only: parameters
use simple_flex_pca_rounds, only: flex_pca_rounds, PCA_STAGE_STATES, FLEX_PCA_PART_MAGIC, flex_pca_part_path
use simple_flex_reconstructor_latent_ops, only: insert_planes_oversamp_multi_scaled_batch
use simple_flex_gpu, only: flex_gpu_available, flex_gpu_insert_begin_f, flex_gpu_insert_batch_f, &
    &flex_gpu_insert_batch_res_f, flex_gpu_insert_end_f, flex_gpu_prep_begin_f, &
    &flex_gpu_prep_batch_f, flex_gpu_prep_fetch_f, flex_gpu_prep_free_f, flex_gpu_prep_ready
use simple_reconstructor,          only: reconstructor
use simple_flex_pca_util,          only: flex_pca_write_state
use simple_flex_pca_rec3D_pcg,     only: reconstruct_flex_weighted_states_pcg
use simple_estimate_ssnr,          only: fsc2optlp_sub, get_resolution
implicit none
character(len=*), parameter :: WEIGHTS_FNAME    = 'flex_pca_round_weights.bin'
private
#include "simple_local_flags.inc"

public :: reconstruct_flex_weighted_states
public :: read_state_weights_round
public :: flex_rec_box, flex_rec_smpd

contains

    !> Kernel-weighted state reconstruction: each state is a weighted backprojection of all particles.
    !! With outvol_even/outvol_odd present the combined, even and odd maps come from one pass: insertion
    !! is linear in the weights and every nonlinear finalisation runs after compress_exp, so
    !! combined = even + odd. Each particle is inserted once, into its own halfset.
    subroutine reconstruct_flex_weighted_states( params, build, pinds, state_weights, nstates, fsc_projfile, &
        &floor_rho, outvol_even, outvol_odd, split_eo , rounds)
        class(flex_pca_rounds), intent(inout) :: rounds
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: pinds(:), nstates
        real,              intent(in)    :: state_weights(:,:)
        type(string), optional, intent(in) :: fsc_projfile
        ! shellwise rho floor before the divide (default off; flex_pca opts in)
        logical,      optional, intent(in) :: floor_rho
        type(string), optional, intent(in) :: outvol_even, outvol_odd
        !! worker-side halfset split flag, read from the round-weights table (a worker has no output
        !! names); both routes must set l_fuse identically or master and workers disagree on the halves
        logical,      optional, intent(in) :: split_eo
        logical :: l_floor_rho, l_reduced, l_fuse, l_state_eofilt, l_state_filt
        type(reconstructor), allocatable :: recs_o(:), recs_c(:), cur(:)
        type(string) :: outvol_bak, state_vol_fname
        integer :: eo_i, iview, nview
        type(reconstructor), allocatable :: state_recs(:)
        type(fplane_type), allocatable :: fpls(:)
        type(image) :: gridcorr_img, state_img
        type(ori) :: orientation
        ! batch buffers for insert_planes_oversamp_multi_scaled_batch; bvalid_c/bvalid_o are disjoint
        type(ori),            allocatable :: borientations(:)
        real(dp),             allocatable :: bscales(:,:)
        logical,              allocatable :: bvalid_c(:), bvalid_o(:)
        real(dp), allocatable :: scales(:)
        real, allocatable :: lowpass_filters(:,:)
        logical, allocatable :: has_lowpass_filter(:)
        integer, allocatable :: lowpass_source_state(:)
        integer :: batchlims(2), batchsz, ibatch, i, iptcl, state, box_rec
        real    :: smpd_rec, smpd_crop_bak
        ! GPU insertion path: one device accumulation over the combined [state_recs, recs_o]
        ! component layout, with halfset routing through the scale slots
        logical :: l_gpu, l_devprep, l_chk, l_fetch_batch
        integer :: ncomp_gpu, envlen, envstat, nyq_unpd
        character(len=8) :: envval
        real(dp), allocatable :: gscales(:,:)
        ! halfset split: two passes over nstates device slots instead of one over 2*nstates
        logical :: l_split, l_forced
        logical, allocatable :: bvalid_p(:)
        integer :: ipass, npass, ios
        integer(kind=8) :: nvox_acc
        real(dp) :: acc_mb_fused, acc_budget_mb
        integer(timer_int_kind) :: t_ins, t_sec
        real(timer_int_kind) :: sec_read, sec_prep, sec_ins
        l_floor_rho = .false.
        if( present(floor_rho) ) l_floor_rho = floor_rho
        l_fuse = present(outvol_even) .and. present(outvol_odd)
        if( present(split_eo) ) l_fuse = l_fuse .or. split_eo
        if( size(pinds)<1 .or. nstates<1 ) THROW_HARD('invalid flex weighted state reconstruction dimensions')
        if( any(shape(state_weights)/=[size(pinds),nstates]) ) THROW_HARD('flex weighted state table mismatch')
        box_rec  = flex_rec_box(params)
        smpd_rec = flex_rec_smpd(params)
        if( box_rec /= params%box_crop )then
            write(logfhandle,'(A,I0,A,F6.3,A,I0,A,F6.3,A)') '>>> FLEX STATE RECONSTRUCTION decoupled box: rec box=',box_rec, &
                &' smpd=',smpd_rec,' A (covariance box=',params%box_crop,' smpd=',params%smpd_crop,' A)'
        endif
        ! rec_states_backend=pcg: the same weighted least-squares problems on reconstructor_pcg with the
        ! support inside the solve; the weights round and the stage fan-out are shared with the gridding
        ! path. This is deliberately NOT rec_backend: the M-step and the state maps are separate
        ! decisions (doc/refactoring_notes/flex_pca_branch_reconciliation_2026_09_15.md 4.3 -- PCG wins
        ! the basis, gridding wins the state maps), so the default here is gridding even under
        ! rec_backend=pcg.
        if( trim(params%rec_states_backend) == 'pcg' )then
            write(logfhandle,'(A)') '>>> FLEX STATE RECONSTRUCTION: kernel PCG backend (rec_states_backend=pcg)'
            call flush(logfhandle)
            if( rounds%is_master() )then
                call write_state_weights_round(pinds, state_weights, size(pinds), nstates, l_fuse)
                call rounds%run_stage(params, PCA_STAGE_STATES, 'state reconstruction')
            endif
            call reconstruct_flex_weighted_states_pcg(params, build, pinds, state_weights, nstates, l_fuse, &
                &box_rec, smpd_rec, outvol_even, outvol_odd, rounds)
            return
        endif
        ! the delivered state maps are always low-passed at their own eo-FSC(0.143) resolution
        ! (user decision 2026-09-16): a poorly determined state must look poorly determined
        l_state_eofilt = .false.
        l_state_filt   = .true.
        write(logfhandle,'(A)') '>>> FLEX_PCA state maps delivered under a per-state low-pass at each state''s own &
            &eo-FSC(0.143) resolution'
        allocate(state_recs(nstates),scales(nstates))
        call prepare_project_fsc_lowpass_filters(params,build,nstates,lowpass_filters,has_lowpass_filter,lowpass_source_state, &
            &fsc_projfile, state_mass=real(sum(state_weights,dim=1)))
        do state=1,nstates
            call init_state_reconstructor(params,build,state_recs(state))
        end do
        if( l_fuse )then
            allocate(recs_o(nstates))
            do state=1,nstates
                call init_state_reconstructor(params,build,recs_o(state))
            end do
        endif
        call init_rec(params,build,MAXIMGBATCHSZ,fpls)
        call prepimgbatch(params,build,MAXIMGBATCHSZ)
        ! prep_imgs4rec builds planes at params%smpd_crop, which fixes the CTF frequencies; the maps
        ! live on the box_rec/smpd_rec lattice, so smpd_crop is pointed at smpd_rec for the batch loop
        ! and restored afterwards (prep_imgs4rec is the only consumer in between; no-op unless decoupled).
        ! Distributed: the master ships the weight table, fans the particle range out and sums the
        ! compressed partials; every nonlinear finalisation then runs once on the global sums.
        l_reduced = .false.
        if( rounds%is_master() )then
            call write_state_weights_round(pinds, state_weights, size(pinds), nstates, l_fuse)
            call rounds%run_stage(params, PCA_STAGE_STATES, 'state reconstruction')
            block
                type(reconstructor) :: rec_read
                type(string) :: pf
                integer :: ipart
                integer(timer_int_kind) :: t_red
                t_red = tic()
                call init_state_reconstructor(params,build,rec_read)
                do ipart = 1, rounds%nparts()
                    do state = 1, nstates
                        ! on a split round each part carries both halfsets; reduce each into its own
                        ! accumulator so combined = even + odd below sums two populated halves
                        do eo_i = 0, merge(1, 0, l_fuse)
                            pf = flex_state_part_fbody(params, ipart, state, eo_i)
                            if( .not. file_exists(pf//MRC_EXT) ) THROW_HARD('missing states part: '//pf%to_char())
                            call rec_read%read(pf//MRC_EXT)
                            call rec_read%read_rho(flex_pca_rho_part_name(pf))
                            if( eo_i == 1 )then
                                call recs_o(state)%sum_reduce(rec_read)
                            else
                                call state_recs(state)%sum_reduce(rec_read)
                            endif
                            call del_file(pf//MRC_EXT)
                            call del_file(flex_pca_rho_part_name(pf))
                            call pf%kill
                        end do
                    end do
                end do
                call rec_read%dealloc_rho; call rec_read%kill
                write(logfhandle,'(A,I0,A,F8.1)') '>>> FLEX_PCA reduced states parts=', &
                    &rounds%nparts(),' seconds=',toc(t_red)
                call flush(logfhandle)
            end block
            l_reduced = .true.
            goto 300
        endif
        smpd_crop_bak    = params%smpd_crop
        params%smpd_crop = smpd_rec
        allocate(borientations(MAXIMGBATCHSZ), bscales(nstates,MAXIMGBATCHSZ))
        allocate(bvalid_c(MAXIMGBATCHSZ), bvalid_o(MAXIMGBATCHSZ), bvalid_p(MAXIMGBATCHSZ))
        ! GPU insertion is on whenever the CUDA build sees a device; SIMPLE_COV_GPU=0 opts out,
        ! =1 warns when no device is available
        l_gpu = flex_gpu_available()
        call get_environment_variable('SIMPLE_COV_GPU', envval, envlen, envstat)
        if( envstat == 0 .and. envlen > 0 )then
            if( trim(adjustl(envval)) == '0' ) l_gpu = .false.
            if( trim(adjustl(envval)) == '1' .and. .not. l_gpu ) write(logfhandle,'(A)') &
                &'>>> FLEX_PCA WARNING: SIMPLE_COV_GPU=1 but no CUDA build/device; CPU insertion'
        endif
        ! halfset split decision: the fused pass keeps 2*nstates device accumulators resident
        ! (2*nstates*nvox*12 bytes), the peak device memory of the program. Since each particle
        ! belongs to one halfset, the halves can be accumulated in two passes over nstates slots for
        ! half the memory and a second read+prep sweep. Engages once the fused footprint exceeds
        ! SIMPLE_FLEX_GPU_ACC_MB (default 2048 MB); SIMPLE_FLEX_GPU_SPLIT_EO=1/0 forces it on/off.
        l_split = .false.
        if( l_gpu .and. l_fuse )then
            nvox_acc = int(ubound(state_recs(1)%cmat_exp,1)-lbound(state_recs(1)%cmat_exp,1)+1,8) &
                &    * int(ubound(state_recs(1)%cmat_exp,2)-lbound(state_recs(1)%cmat_exp,2)+1,8) &
                &    * int(ubound(state_recs(1)%cmat_exp,3)-lbound(state_recs(1)%cmat_exp,3)+1,8)
            acc_mb_fused  = real(2*nstates,dp) * real(nvox_acc,dp) * 12.d0 / 1048576.d0
            acc_budget_mb = 2048.d0
            call get_environment_variable('SIMPLE_FLEX_GPU_ACC_MB', envval, envlen, envstat)
            if( envstat == 0 .and. envlen > 0 )then
                read(envval, *, iostat=ios) acc_budget_mb
                if( ios /= 0 ) acc_budget_mb = 2048.d0
            endif
            l_split   = acc_mb_fused > acc_budget_mb
            l_forced  = .false.
            call get_environment_variable('SIMPLE_FLEX_GPU_SPLIT_EO', envval, envlen, envstat)
            if( envstat == 0 .and. envlen > 0 )then
                if( trim(adjustl(envval)) == '1' )then
                    l_forced = .not. l_split
                    l_split  = .true.
                endif
                if( trim(adjustl(envval)) == '0' ) l_split = .false.
            endif
            if( l_split .and. l_forced )then
                write(logfhandle,'(A,F8.1,A)') &
                    &'>>> FLEX_PCA STATEREC HALFSET SPLIT ON (forced): fused accumulators would be ', &
                    &acc_mb_fused,' MB; two passes over half the slots'
            else if( l_split )then
                write(logfhandle,'(A,F8.1,A,F8.1,A)') &
                    &'>>> FLEX_PCA STATEREC HALFSET SPLIT ON: fused accumulators would be ', &
                    &acc_mb_fused,' MB > budget ',acc_budget_mb,' MB; two passes over half the slots'
            endif
        endif
        npass = 1
        if( l_split ) npass = 2
        if( l_gpu )then
            ncomp_gpu = nstates
            if( l_fuse .and. .not. l_split ) ncomp_gpu = 2*nstates
            write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA STATEREC GPU insertion ON (', ncomp_gpu, &
                &' device accumulators, ',npass,' pass(es))'
            call flush(logfhandle)
            ! under the split each pass re-inits the accumulators for its own halfset
            if( .not. l_split ) call flex_gpu_insert_begin_f(state_recs, ncomp_gpu)
            allocate(gscales(ncomp_gpu, MAXIMGBATCHSZ), source=0.d0)
        endif
        ! device prep lifecycle: same taper geometry as the covariance stages (non-cached path,
        ! boxpd planes); the premultiplied reconstruction pair stays resident on the device
        l_devprep = .false.
        if( flex_gpu_available() .and. .not. flex_gpu_prep_ready() )then
            call get_environment_variable('SIMPLE_COV_GPU_PREP', envval, envlen, envstat)
            l_devprep = .true.
            if( envstat == 0 .and. envlen > 0 )then
                if( trim(adjustl(envval)) == '0' ) l_devprep = .false.
            endif
            if( l_devprep .and. params%l_ml_reg )then
                l_devprep = allocated(build%esig%sigma2_noise)
            endif
            if( l_devprep )then
                call flex_gpu_prep_begin_f(build%lmsk, params%box, params%boxpd, &
                    &MAXIMGBATCHSZ, 0.0, .true.)
                write(logfhandle,'(A)') '>>> FLEX_PCA STATEREC DEVICE PREP ON (resident hand-off)'
                call flush(logfhandle)
            endif
        endif
        l_chk = .false.
        if( l_devprep )then
            call get_environment_variable('SIMPLE_COV_PREP_CHECK', envval, envlen, envstat)
            if( envstat == 0 .and. envlen > 0 )then
                if( trim(adjustl(envval)) == '1' ) l_chk = .true.
            endif
        endif
        nyq_unpd = 0
        t_ins = tic()
        sec_read = 0.; sec_prep = 0.; sec_ins = 0.
        do ipass = 1, npass
            ! pass 1 accumulates the even halfset into state_recs, pass 2 the odd into recs_o, over
            ! the same nstates device slots; without the split one pass fills the 2*nstates layout
            if( l_split )then
                if( ipass == 1 )then
                    call flex_gpu_insert_begin_f(state_recs, ncomp_gpu)
                else
                    call flex_gpu_insert_begin_f(recs_o, ncomp_gpu)
                endif
            endif
            do ibatch=1,size(pinds),MAXIMGBATCHSZ
                batchlims=[ibatch,min(size(pinds),ibatch+MAXIMGBATCHSZ-1)]
                batchsz=batchlims(2)-batchlims(1)+1
                t_sec = tic()
                if( params%l_ptcl_src_den )then
                    call discrete_read_imgbatch_source(params,build,'den',batchsz, &
                        &pinds(batchlims(1):batchlims(2)),[1,batchsz],build%imgbatch(:batchsz))
                else
                    call discrete_read_imgbatch(params,build,size(pinds),pinds,batchlims)
                endif
                sec_read = sec_read + toc(t_sec)
                t_sec = tic()
                if( l_devprep )then
                    ! fetch only when a host consumer needs the planes: CPU insertion, or the
                    ! one-time cross-check; the GPU insert consumes them resident
                    l_fetch_batch = (.not. l_gpu) .or. (ibatch == 1 .and. ipass == 1 .and. l_chk)
                    call prep_imgs4rec_dev(params,build,batchsz,build%imgbatch(:batchsz), &
                        &pinds(batchlims(1):batchlims(2)),fpls(:batchsz), fetch=l_fetch_batch, &
                        &nyq_unpd_out=nyq_unpd)
                    if( ibatch == 1 .and. ipass == 1 .and. l_chk )then
                        ! one-time cross-check of the delivered reconstruction planes vs the CPU prep
                        block
                            type(fplane_type), allocatable :: fpls_chk(:)
                            real    :: chk_ey, chk_eq, chk_den
                            integer :: ichk
                            allocate(fpls_chk(batchsz))
                            call prep_imgs4rec(params,build,batchsz,build%imgbatch(:batchsz), &
                                &pinds(batchlims(1):batchlims(2)),fpls_chk)
                            chk_ey = 0.; chk_eq = 0.; chk_den = 1.e-12
                            do ichk = 1, batchsz
                                chk_ey  = max(chk_ey, maxval(abs(fpls_chk(ichk)%cmplx_plane - &
                                    &fpls(ichk)%cmplx_plane)))
                                chk_eq  = max(chk_eq, maxval(abs(fpls_chk(ichk)%ctfsq_plane - &
                                    &fpls(ichk)%ctfsq_plane)))
                                chk_den = max(chk_den, maxval(abs(fpls_chk(ichk)%cmplx_plane)))
                                deallocate(fpls_chk(ichk)%cmplx_plane, fpls_chk(ichk)%ctfsq_plane)
                            end do
                            write(logfhandle,'(A,ES10.2,A,ES10.2,A,ES10.2)') &
                                &'>>> FLEX_PCA STATEREC PREP CHECK: max|d y|=', chk_ey, &
                                &'  rel=', chk_ey/chk_den, '  max|d ctfsq|=', chk_eq
                            call flush(logfhandle)
                            deallocate(fpls_chk)
                        end block
                    endif
                else
                    call prep_imgs4rec(params,build,batchsz,build%imgbatch(:batchsz), &
                        &pinds(batchlims(1):batchlims(2)),fpls(:batchsz))
                endif
                sec_prep = sec_prep + toc(t_sec)
                ! gather serially (get_ori/get_eo are not guaranteed thread-safe), then one parallel
                ! region per target; under l_fuse each halfset gets its own mask and call
                do i=1,batchsz
                    iptcl=pinds(batchlims(1)+i-1)
                    call build%spproj_field%get_ori(iptcl,orientation)
                    bscales(:,i) = real(state_weights(batchlims(1)+i-1,:),dp)
                    bvalid_c(i)  = .not. orientation%isstatezero()
                    bvalid_o(i)  = .false.
                    if( bvalid_c(i) .and. l_fuse )then
                        ! one insertion per particle, into its own halfset
                        eo_i = build%spproj_field%get_eo(iptcl)
                        if( eo_i == 1 )then
                            bvalid_o(i) = .true.
                            bvalid_c(i) = .false.
                        endif
                    endif
                    call borientations(i)%copy(orientation)
                end do
                t_sec = tic()
                if( l_gpu )then
                    ! halfset routing through the component slots: valid_c rows live in 1..nstates,
                    ! valid_o rows in nstates+1..2*nstates; the batch packer compacts on nonzero scales
                    gscales(:,:batchsz) = 0.d0
                    if( l_split )then
                        ! this pass owns one halfset; the other half contributes nothing to it
                        do i = 1, batchsz
                            if( ipass == 1 )then
                                bvalid_p(i) = bvalid_c(i)
                            else
                                bvalid_p(i) = bvalid_o(i)
                            endif
                            if( bvalid_p(i) ) gscales(1:nstates,i) = bscales(:,i)
                        end do
                    else
                        do i = 1, batchsz
                            bvalid_p(i) = bvalid_c(i) .or. bvalid_o(i)
                            if( bvalid_c(i) )then
                                gscales(1:nstates,i) = bscales(:,i)
                            else if( bvalid_o(i) )then
                                gscales(nstates+1:2*nstates,i) = bscales(:,i)
                            endif
                        end do
                    endif
                    if( l_devprep )then
                        ! planes are resident from the device prep of this batch
                        call flex_gpu_insert_batch_res_f(build%pgrpsyms, borientations(:batchsz), &
                            &gscales(:,:batchsz), gscales(:,:batchsz), &
                            &bvalid_p(:batchsz), batchsz, nyq_unpd)
                    else
                        call flex_gpu_insert_batch_f(build%pgrpsyms, borientations(:batchsz), &
                            &fpls(:batchsz), gscales(:,:batchsz), gscales(:,:batchsz), &
                            &bvalid_p(:batchsz), batchsz)
                    endif
                else
                    if( l_fuse .and. any(bvalid_o(:batchsz)) )then
                        call insert_planes_oversamp_multi_scaled_batch(recs_o, build%pgrpsyms, &
                            &borientations(:batchsz), fpls(:batchsz), bscales(:,:batchsz), &
                            &bscales(:,:batchsz), bvalid_o(:batchsz), batchsz)
                    endif
                    if( any(bvalid_c(:batchsz)) )then
                        call insert_planes_oversamp_multi_scaled_batch(state_recs, build%pgrpsyms, &
                            &borientations(:batchsz), fpls(:batchsz), bscales(:,:batchsz), &
                            &bscales(:,:batchsz), bvalid_c(:batchsz), batchsz)
                    endif
                endif
                sec_ins = sec_ins + toc(t_sec)
            end do
            ! fetch this halfset before the next pass reuses the slots
            if( l_split )then
                if( ipass == 1 )then
                    call flex_gpu_insert_end_f(state_recs)
                else
                    call flex_gpu_insert_end_f(recs_o)
                endif
            endif
        end do
        if( l_devprep ) call flex_gpu_prep_free_f
        if( l_gpu )then
            if( .not. l_split )then
                if( l_fuse )then
                    call flex_gpu_insert_end_f(state_recs, recs_o)
                else
                    call flex_gpu_insert_end_f(state_recs)
                endif
            endif
            deallocate(gscales)
        endif
        write(logfhandle,'(A,F8.1,A,I0)') '>>> FLEX_PCA STATEREC read+prep+insert seconds=', &
            &toc(t_ins), '  gpu=', merge(1, 0, l_gpu)
        write(logfhandle,'(A,F7.1,A,F7.1,A,F7.1)') '>>> FLEX_PCA STATEREC SPLIT (seconds): read=', &
            &sec_read,'  prep=',sec_prep,'  insert=',sec_ins
        call flush(logfhandle)
        do i = 1, MAXIMGBATCHSZ
            call borientations(i)%kill
        end do
        deallocate(borientations, bscales, bvalid_c, bvalid_o, bvalid_p)
        call orientation%kill
        params%smpd_crop = smpd_crop_bak
        call cleanup_rec_buffers(build,fpls)
        if( rounds%is_worker() )then
            block
                type(string) :: pf
                do state=1,nstates
                    call state_recs(state)%compress_exp
                    pf = flex_state_part_fbody(params, params%part, state, 0)
                    call state_recs(state)%write(pf//MRC_EXT, del_if_exists=.true.)
                    call state_recs(state)%write_rho(flex_pca_rho_part_name(pf))
                    call pf%kill
                    if( l_fuse )then
                        call recs_o(state)%compress_exp
                        pf = flex_state_part_fbody(params, params%part, state, 1)
                        call recs_o(state)%write(pf//MRC_EXT, del_if_exists=.true.)
                        call recs_o(state)%write_rho(flex_pca_rho_part_name(pf))
                        call pf%kill
                    endif
                end do
            end block
            do state=1,nstates
                call state_recs(state)%dealloc_rho; call state_recs(state)%kill
                if( l_fuse )then
                    call recs_o(state)%dealloc_rho; call recs_o(state)%kill
                endif
            end do
            deallocate(state_recs,scales,lowpass_filters,has_lowpass_filter,lowpass_source_state)
            if( l_fuse ) deallocate(recs_o)
            return
        endif
        300     continue
        gridcorr_img=prep3D_inv_kbenvelope4mul([box_rec,box_rec,box_rec], smpd_rec)
        outvol_bak = params%outvol
        nview = 1
        if( l_fuse )then
            ! per-state eo-FSC delivery: the even and odd maps of one state carry identical kernel
            ! weights, so their FSC measures that state's own information content (kernel regression
            ! borrows strength across the dataset, so a small state can resolve beyond its bin count).
            ! The project FSC is only a fallback when the state FSC is unmeasurable.
            block
                type(image) :: img_e, img_o, img_c, msk_e, msk_o
                real, allocatable :: fsc_eo(:), res_arr(:), filt_half(:), filt_merged(:)
                real    :: kc_lp
                integer :: k_lp
                real    :: fsc05, fsc0143, mskrad
                integer :: filtsz_del, iv2
                logical :: l_eo_fsc
                do state=1,nstates
                    if( .not. l_reduced )then
                        call state_recs(state)%compress_exp
                        call recs_o(state)%compress_exp
                    endif
                end do
                allocate(recs_c(nstates))
                do state=1,nstates
                    call init_state_reconstructor(params,build,recs_c(state))
                    call recs_c(state)%sum_reduce(state_recs(state))
                    call recs_c(state)%sum_reduce(recs_o(state))
                end do
                l_reduced  = .true.
                filtsz_del = fdim(box_rec) - 1
                allocate(fsc_eo(filtsz_del), filt_half(filtsz_del), filt_merged(filtsz_del))
                ! delivery mask at mskdiam, capped at the box edge; a box/2 mask lets solvent noise
                ! dominate the state FSC
                mskrad = min(real(box_rec/2) - COSMSKHALFWIDTH - 1., 0.5*params%mskdiam/smpd_rec)
                do state=1,nstates
                    ! finalize the three views of this state (destructive on the reconstructors)
                    call finalize_state_rec(state_recs(state), gridcorr_img, l_floor_rho, img_e)
                    call finalize_state_rec(recs_o(state),     gridcorr_img, l_floor_rho, img_o)
                    call finalize_state_rec(recs_c(state),     gridcorr_img, l_floor_rho, img_c)
                    ! per-state eo FSC on masked copies (delivery mask, background zeroed)
                    call msk_e%copy(img_e)
                    call msk_e%zero_background
                    call msk_e%mask3D_soft(mskrad, backgr=0.)
                    call msk_o%copy(img_o)
                    call msk_o%zero_background
                    call msk_o%mask3D_soft(mskrad, backgr=0.)
                    call msk_e%fft
                    call msk_o%fft
                    call msk_e%fsc(msk_o, fsc_eo)
                    call msk_e%kill
                    call msk_o%kill
                    l_eo_fsc = any(fsc_eo > 0.143)
                    if( l_eo_fsc )then
                        res_arr = img_e%get_res()
                        call get_resolution(fsc_eo, res_arr, fsc05, fsc0143)
                        if( l_state_eofilt )then
                            ! SIMPLE_COV_STATE_EOFILT=1: per-state eo-FSC optimal filter
                            call fsc2optlp_sub(filtsz_del, fsc_eo, filt_half,   merged=.false.)
                            call fsc2optlp_sub(filtsz_del, fsc_eo, filt_merged, merged=.true.)
                            write(logfhandle,'(A,I3,A,F7.2,A,F7.2,A)') '>>> FLEX STATE eo-FSC state=', &
                                &state,'  res(0.143)=',fsc0143,' A  res(0.5)=',fsc05,' A -- per-state optimal filter applied'
                        else
                            ! default: 8th-order Butterworth low-pass at this state's eo-FSC(0.143) resolution,
                            ! the same filter for the combined map and both halves
                            kc_lp = real(box_rec) * smpd_rec / fsc0143
                            do k_lp = 1, filtsz_del
                                filt_merged(k_lp) = 1.0 / (1.0 + (real(k_lp)/max(kc_lp,1.0))**8)
                            end do
                            filt_half = filt_merged
                            write(logfhandle,'(A,I3,A,F7.2,A,F7.2,A)') '>>> FLEX STATE eo-FSC state=', &
                                &state,'  res(0.143)=',fsc0143,' A  res(0.5)=',fsc05,' A -- low-pass at the state eo-FSC(0.143) applied'
                        endif
                        deallocate(res_arr)
                    else
                        write(logfhandle,'(A,I3,A)') '>>> FLEX STATE eo-FSC state=',state, &
                            &'  unmeasurable (no shell above 0.143); project-FSC low-pass if the project has one, else unfiltered'
                    endif
                    call flush(logfhandle)
                    do iv2 = 1, 3
                        select case(iv2)
                        case(1); call state_img%copy(img_c); params%outvol = outvol_bak
                        case(2); call state_img%copy(img_e); params%outvol = outvol_even
                        case(3); call state_img%copy(img_o); params%outvol = outvol_odd
                        end select
                        if( .not. l_state_filt )then
                        ! SIMPLE_COV_STATE_FILT=0: no filter on the delivered maps
                        else if( l_eo_fsc )then
                            if( iv2 == 1 )then
                                call state_img%apply_filter(filt_merged)
                            else
                                call state_img%apply_filter(filt_half)
                            endif
                        else if( has_lowpass_filter(state) )then
                            call state_img%apply_filter(lowpass_filters(:,state))
                        endif
                        call state_img%zero_background
                        call state_img%mask3D_soft(mskrad, backgr=0.)
                        call flex_pca_write_state(params, state_img, state, state_vol_fname)
                        if( iv2 == 1 )then
                            call build%spproj%add_vol2os_out(state_vol_fname, state_img%get_smpd(), state, 'vol_flex',&
                                &box=state_img%get_box())
                        endif
                        call state_img%kill
                    end do
                    call img_e%kill
                    call img_o%kill
                    call img_c%kill
                    call state_recs(state)%dealloc_rho; call state_recs(state)%kill
                    call recs_o(state)%dealloc_rho;     call recs_o(state)%kill
                    call recs_c(state)%dealloc_rho;     call recs_c(state)%kill
                end do
                deallocate(fsc_eo, filt_half, filt_merged)
            end block
        else
            do iview = 1, nview
                call move_alloc(state_recs, cur)
                do state=1,nstates
                    if( .not. l_reduced ) call cur(state)%compress_exp
                    ! kernel weights are mostly near zero, so rho is small and an unfloored divide
                    ! amplifies noise wherever occupancy is low
                    if( l_floor_rho ) call cur(state)%floor_rho_shellwise
                    call cur(state)%sampl_dens_correct
                    call cur(state)%ifft
                    call cur(state)%mul(gridcorr_img)
                    call state_img%copy(cur(state))
                    if( has_lowpass_filter(state) )then
                        call state_img%apply_filter(lowpass_filters(:,state))
                        write(logfhandle,'(A,I0,A,I0)') '>>> FLEX PRE-IMAGE applied project-FSC low-pass filter to state=',state, &
                            &' using_source_state=',lowpass_source_state(state)
                    endif
                    ! background removal + soft spherical mask: each state carries a different total kernel
                    ! weight, so without both the states differ by a baseline offset and solvent noise rather
                    ! than by conformation. Radius is capped at the broadest soft-maskable sphere in the box.
                    call state_img%zero_background
                    call state_img%mask3D_soft(min(real(box_rec/2) - COSMSKHALFWIDTH - 1., 0.5*params%mskdiam/smpd_rec), backgr=0.)
                    call flex_pca_write_state(params, state_img, state, state_vol_fname)
                    call build%spproj%add_vol2os_out(state_vol_fname, state_img%get_smpd(), state, 'vol_flex',&
                        &box=state_img%get_box())
                    ! clean up
                    call state_img%kill
                    call cur(state)%dealloc_rho
                    call cur(state)%kill
                end do
            end do
        endif
        call build%spproj%write_segment_inside('out', params%projfile)
        params%outvol = outvol_bak
        ! clean up
        call gridcorr_img%kill
        call state_vol_fname%kill
        deallocate(scales,lowpass_filters,has_lowpass_filter,lowpass_source_state)
    end subroutine reconstruct_flex_weighted_states

    !> Part-file body for one worker's partial reconstruction of one state; eo=0 is the even/single
    !! accumulator, eo=1 the odd one (a split round writes both per state).
    function flex_state_part_fbody( params, part, state, eo ) result( fbody )
        class(parameters), intent(in) :: params
        integer,           intent(in) :: part, state, eo
        type(string) :: fbody
        fbody = string('flex_pca_statepart')//int2str_pad(part,max(1,params%numlen))// &
            &'_'//int2str_pad(state,2)
        if( eo == 1 ) fbody = fbody//'_o'
        fbody = flex_pca_part_path(fbody%to_char())
    end function flex_state_part_fbody

    !> rho companion of a state part: same directory, 'rho_' on the file name only
    function flex_pca_rho_part_name( pf ) result( fn )
        type(string), intent(in) :: pf
        type(string) :: fn
        character(len=:), allocatable :: c
        integer :: k
        c = pf%to_char()
        k = index(c, '/', back=.true.)
        fn = string(c(1:k))//'rho_'//c(k+1:)//MRC_EXT
    end function flex_pca_rho_part_name

    subroutine prepare_project_fsc_lowpass_filters( params, build, nstates, lowpass_filters, has_filter, source_state, &
        &fsc_projfile , state_mass)
        class(parameters), intent(in) :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in) :: nstates
        real, allocatable, intent(out) :: lowpass_filters(:,:)
        logical, allocatable, intent(out) :: has_filter(:)
        integer, allocatable, intent(out) :: source_state(:)
        type(string), optional, intent(in) :: fsc_projfile
        !> per-state effective weight mass (population-scaled fallback filter)
        real, optional,    intent(in) :: state_mass(:)
        type(sp_project) :: spproj
        type(string) :: fsc_fname, imgkind_here, proj_for_fsc
        real, allocatable :: fsc(:)
        integer :: filtsz, state, fsc_box, i, state1_fsc_count
        logical :: out_loaded
        ! sized to the delivered map (box_rec)
        filtsz=fdim(flex_rec_box(params))-1
        allocate(lowpass_filters(filtsz,nstates),has_filter(nstates),source_state(nstates))
        lowpass_filters=0.
        has_filter=.false.
        source_state=0
        if( filtsz<1 ) return
        proj_for_fsc=params%projfile
        if( present(fsc_projfile) )then
            if( len_trim(fsc_projfile%to_char())>0 ) proj_for_fsc=fsc_projfile
        endif
        if( .not.file_exists(proj_for_fsc) )then
            write(logfhandle,'(A)') '>>> FLEX PRE-IMAGE project-FSC low-pass unavailable (projfile not found); states rely on their own eo-FSC'
            has_filter = .false.   ! the state's own eo-FSC decides
            call proj_for_fsc%kill
            return
        endif
        call spproj%read_segment('out',proj_for_fsc)
        out_loaded=spproj%os_out%get_noris()>0
        if( .not.out_loaded )then
            write(logfhandle,'(A)') '>>> FLEX PRE-IMAGE project-FSC low-pass unavailable (empty out segment); states rely on their own eo-FSC'
            has_filter = .false.   ! the state's own eo-FSC decides
            call proj_for_fsc%kill
            call spproj%kill
            return
        endif
        state1_fsc_count=0
        do i=1,spproj%os_out%get_noris()
            if( .not.spproj%os_out%isthere(i,'imgkind') ) cycle
            call spproj%os_out%getter(i,'imgkind',imgkind_here)
            if( imgkind_here%to_char()/='fsc' ) cycle
            if( spproj%os_out%get_state(i)==1 ) state1_fsc_count=state1_fsc_count+1
        end do
        if( state1_fsc_count/=1 )then
            write(logfhandle,'(A,I0)') '>>> FLEX PRE-IMAGE project-FSC low-pass unavailable (state=1 FSC count=', &
                &state1_fsc_count
            write(logfhandle,'(A)') '>>>   ); states rely on their own eo-FSC'
            has_filter = .false.   ! the state's own eo-FSC decides
            call imgkind_here%kill
            call proj_for_fsc%kill
            call spproj%kill
            return
        endif
        call spproj%get_fsc(1,fsc_fname,fsc_box)
        if( .not.file_exists(fsc_fname) )then
            write(logfhandle,'(A)') '>>> FLEX PRE-IMAGE project-FSC low-pass unavailable (state=1 FSC file missing); states rely on their own eo-FSC'
            has_filter = .false.   ! the state's own eo-FSC decides
            call imgkind_here%kill
            call proj_for_fsc%kill
            call spproj%kill
            return
        endif
        fsc=file2rarr(fsc_fname)
        if( size(fsc)/=filtsz )then
            write(logfhandle,'(A,I0,A,I0)') '>>> FLEX PRE-IMAGE project-FSC low-pass unavailable (state=1 FSC size mismatch; fsc_nyq=', &
                &size(fsc),' model_nyq=',filtsz
            write(logfhandle,'(A)') '>>>   ); states rely on their own eo-FSC'
            has_filter = .false.   ! the state's own eo-FSC decides
            deallocate(fsc)
            call fsc_fname%kill
            call imgkind_here%kill
            call proj_for_fsc%kill
            call spproj%kill
            return
        endif
        do state=1,nstates
            call fsc2optlp_sub(filtsz,fsc,lowpass_filters(:,state),merged=.false.)
            has_filter(state)=any(lowpass_filters(:,state)>0.)
            source_state(state)=1
        end do
        deallocate(fsc)
        call fsc_fname%kill
        call imgkind_here%kill
        call proj_for_fsc%kill
        call spproj%kill
    end subroutine prepare_project_fsc_lowpass_filters


    !> Box/sampling of the delivered state maps. Decoupled from box_crop because the embedding is a
    !! low-frequency object while the state maps are plain backprojections that carry signal beyond
    !! the covariance band. With box_rec==box_crop this is a no-op.
    pure integer function flex_rec_box( params ) result( box_rec )
        class(parameters), intent(in) :: params
        box_rec = params%box_crop
        if( params%box_rec >= 1 ) box_rec = params%box_rec
    end function flex_rec_box

    pure real function flex_rec_smpd( params ) result( smpd_rec )
        class(parameters), intent(in) :: params
        smpd_rec = params%smpd_crop
        if( params%box_rec >= 1 .and. params%smpd_rec > 0. ) smpd_rec = params%smpd_rec
    end function flex_rec_smpd

    !> One state view: rho floor (opt-in), density correction, ifft, gridding correction -> image.
    !! Destructive on the reconstructor's Fourier state.
    subroutine finalize_state_rec( rec, gridcorr_img, l_floor_rho, img )
        type(reconstructor), intent(inout) :: rec
        type(image),         intent(in)    :: gridcorr_img
        logical,             intent(in)    :: l_floor_rho
        type(image),         intent(inout) :: img
        if( l_floor_rho ) call rec%floor_rho_shellwise
        call rec%sampl_dens_correct
        call rec%ifft
        call rec%mul(gridcorr_img)
        call img%copy(rec)
    end subroutine finalize_state_rec

    subroutine init_state_reconstructor( params, build, state_rec )
        class(parameters), intent(inout) :: params
        class(builder), intent(inout) :: build
        type(reconstructor), intent(inout) :: state_rec
        integer :: box_rec
        box_rec = flex_rec_box(params)
        call state_rec%new([box_rec,box_rec,box_rec],flex_rec_smpd(params))
        call state_rec%alloc_rho(params,build%spproj,expand=.true.)
        call state_rec%reset
        call state_rec%reset_exp
    end subroutine init_state_reconstructor

    !> Device variant of prep_imgs4rec (non-cached path): the taper->norm->pad->FFT->plane chain runs
    !! on the GPU and the premultiplied pair (conj(T)y/sigma2, CTF^2/sigma2) is unpacked into
    !! fplane_type. Equivalent to the CPU generator: values on the OSMPL_PAD_FAC lattice, zeros elsewhere.
    subroutine prep_imgs4rec_dev( params, build, nptcls, ptcl_imgs, pinds, fplanes, fetch, &
        &nyq_unpd_out )
        use simple_ftiter,  only: ftiter
        use simple_math,    only: ceil_div, floor_div
        use simple_math_ft, only: resample_sigma2
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: nptcls
        class(image),      intent(inout) :: ptcl_imgs(nptcls)
        integer,           intent(in)    :: pinds(nptcls)
        type(fplane_type), intent(inout) :: fplanes(nptcls)
        logical, optional, intent(in)    :: fetch        !< .false. = leave planes resident only
        integer, optional, intent(out)   :: nyq_unpd_out !< shared plane band / OSMPL_PAD_FAC
        type(ctfparams) :: ctfp_arr(nptcls)
        real            :: shf(2,nptcls)
        logical         :: vld(nptcls)
        complex(sp), allocatable :: plc(:,:,:)
        real(sp),    allocatable :: plct(:,:,:)
        real,        allocatable :: sig2_ups(:,:)
        type(ftiter) :: fit_pd, fit_cr
        real    :: shconst_pd(3)
        integer :: kfromto(2), frlims_pd(3,2), i, h, k, hlo, hhi, klo, nyqpd, signyq
        integer :: hmin, hmax, kmin, iptcl, pf
        logical :: l_fresh, l_fetch
        l_fetch = .true.
        if( present(fetch) ) l_fetch = fetch
        pf      = OSMPL_PAD_FAC
        kfromto = build%esig%get_kfromto()
        call fit_pd%new([params%boxpd, params%boxpd, 1], params%smpd_crop)
        frlims_pd = fit_pd%loop_lims(3)
        nyqpd     = fit_pd%get_lfny(1)
        if( present(nyq_unpd_out) ) nyq_unpd_out = max(1, nyqpd / pf)
        if( params%l_ml_reg )then
            ! the CPU generator's sigma_nyq comes from the padded box / OSMPL_PAD_FAC
            call fit_cr%new([params%boxpd/pf, params%boxpd/pf, 1], params%smpd_crop)
            signyq = fit_cr%get_lfny(1)
            allocate(sig2_ups(0:nyqpd, nptcls), source=1.0)
        endif
        vld = .true.
        !$omp parallel do default(shared) private(i,iptcl) schedule(static) proc_bind(close)
        do i = 1, nptcls
            iptcl       = pinds(i)
            ctfp_arr(i) = build%spproj%get_ctfparams(params%oritype, iptcl)
            shf(:,i)    = build%spproj_field%get_2Dshift(iptcl)
            if( params%l_ml_reg )then
                call resample_sigma2(kfromto(1), signyq, &
                    &build%esig%sigma2_noise(kfromto(1):kfromto(2), iptcl), nyqpd, &
                    &real(signyq)/real(nyqpd), sig2_ups(:,i))
            endif
        end do
        !$omp end parallel do
        if( params%l_ml_reg )then
            call flex_gpu_prep_batch_f(ptcl_imgs, ctfp_arr, shf, vld, nptcls, params%box, &
                &frlims_pd, nyqpd, sig2_ups=sig2_ups)
        else
            call flex_gpu_prep_batch_f(ptcl_imgs, ctfp_arr, shf, vld, nptcls, params%box, &
                &frlims_pd, nyqpd)
        endif
        if( .not. l_fetch )then
            ! consumer is the resident insert; planes stay on device
            if( allocated(sig2_ups) ) deallocate(sig2_ups)
            return
        endif
        hlo = ceil_div (frlims_pd(1,1), pf); hhi = floor_div(frlims_pd(1,2), pf)
        klo = ceil_div (frlims_pd(2,1), pf)
        allocate(plc(hlo:hhi, klo:0, nptcls), plct(hlo:hhi, klo:0, nptcls))
        call flex_gpu_prep_fetch_f(plc, plct, nptcls, hlo, hhi, klo)
        hmin = frlims_pd(1,1); hmax = frlims_pd(1,2); kmin = frlims_pd(2,1)
        shconst_pd      = 0.
        shconst_pd(1:2) = PI/real(params%boxpd/2)
        !$omp parallel do default(shared) private(i,h,k,l_fresh) schedule(static) proc_bind(close)
        do i = 1, nptcls
            l_fresh = .not. allocated(fplanes(i)%cmplx_plane)
            if( .not. l_fresh )then
                if( lbound(fplanes(i)%cmplx_plane,1) /= hmin .or. &
                    &ubound(fplanes(i)%cmplx_plane,1) /= hmax .or. &
                    &lbound(fplanes(i)%cmplx_plane,2) /= kmin .or. &
                    &ubound(fplanes(i)%cmplx_plane,2) /= 0    .or. &
                    &allocated(fplanes(i)%transfer_plane) )then
                    if( allocated(fplanes(i)%cmplx_plane) )    deallocate(fplanes(i)%cmplx_plane)
                    if( allocated(fplanes(i)%ctfsq_plane) )    deallocate(fplanes(i)%ctfsq_plane)
                    if( allocated(fplanes(i)%transfer_plane) ) deallocate(fplanes(i)%transfer_plane)
                    l_fresh = .true.
                endif
            endif
            if( l_fresh )then
                allocate(fplanes(i)%cmplx_plane(hmin:hmax, kmin:0), source=cmplx(0.,0.))
                allocate(fplanes(i)%ctfsq_plane(hmin:hmax, kmin:0), source=0.)
            endif
            do k = klo, 0
                do h = hlo, hhi
                    fplanes(i)%cmplx_plane(pf*h, pf*k) = plc(h,k,i)
                    fplanes(i)%ctfsq_plane(pf*h, pf*k) = plct(h,k,i)
                end do
            end do
            fplanes(i)%frlims  = frlims_pd
            fplanes(i)%nyq     = nyqpd
            fplanes(i)%shconst = shconst_pd
        end do
        !$omp end parallel do
        deallocate(plc, plct)
        if( allocated(sig2_ups) ) deallocate(sig2_ups)
    end subroutine prep_imgs4rec_dev


    !> Per-round state weights, written with the global pinds so a worker matches its rows by index.
    !! Rewritten before every state round since each round uses a different weight table.
    subroutine write_state_weights_round( pinds, weights, nptcls, nstates, split_eo )
        integer, intent(in) :: pinds(:), nptcls, nstates
        real,    intent(in) :: weights(:,:)
        !! .true. when this round accumulates even and odd into separate reconstructors; a worker
        !! cannot infer it from params%stage (every state round arrives as PCA_STAGE_STATES)
        logical, intent(in) :: split_eo
        type(string) :: fname, tmp_fname
        integer :: funit, io_stat, eo_flag
        fname     = string(WEIGHTS_FNAME)
        tmp_fname = fname//'.tmp'
        eo_flag   = merge(1, 0, split_eo)
        call fopen(funit, file=tmp_fname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk('write_state_weights_round; open', io_stat)
        write(funit, iostat=io_stat) FLEX_PCA_PART_MAGIC, nptcls, nstates, eo_flag
        call fileiochk('write_state_weights_round; header', io_stat)
        write(funit, iostat=io_stat) pinds(1:nptcls)
        call fileiochk('write_state_weights_round; pinds', io_stat)
        write(funit, iostat=io_stat) weights(1:nptcls,1:nstates)
        call fileiochk('write_state_weights_round; weights', io_stat)
        call fclose(funit)
        call simple_rename(tmp_fname, fname)
        call fname%kill; call tmp_fname%kill
    end subroutine write_state_weights_round

    !> Worker side: return the weight rows for this part's pinds, in this part's order.
    subroutine read_state_weights_round( my_pinds, my_nptcls, weights_out, nstates, split_eo )
        integer,           intent(in)  :: my_pinds(:), my_nptcls
        real, allocatable, intent(out) :: weights_out(:,:)
        integer,           intent(out) :: nstates
        logical,           intent(out) :: split_eo !< see write_state_weights_round
        type(string) :: fname
        integer, allocatable :: gpinds(:)
        real,    allocatable :: gw(:,:)
        integer :: funit, io_stat, magic, gn, i, j, hit, eo_flag
        logical :: l_sorted
        fname = string(WEIGHTS_FNAME)
        if( .not. file_exists(fname) ) THROW_HARD('flex_pca worker found no '//WEIGHTS_FNAME)
        call fopen(funit, file=fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
        call fileiochk('read_state_weights_round; open', io_stat)
        read(funit, iostat=io_stat) magic, gn, nstates, eo_flag
        call fileiochk('read_state_weights_round; header', io_stat)
        if( magic /= FLEX_PCA_PART_MAGIC ) THROW_HARD('bad round-weights magic')
        split_eo = eo_flag == 1
        allocate(gpinds(gn), gw(gn,nstates))
        read(funit, iostat=io_stat) gpinds
        call fileiochk('read_state_weights_round; pinds', io_stat)
        read(funit, iostat=io_stat) gw
        call fileiochk('read_state_weights_round; weights', io_stat)
        call fclose(funit)
        allocate(weights_out(my_nptcls,nstates), source=0.)
        ! both lists are ascending project rows, so a merge walk replaces the O(N_local x N_global)
        ! scan; falls back to the scan if either list is unordered
        l_sorted = .true.
        do i = 2, my_nptcls
            if( my_pinds(i) <= my_pinds(i-1) ) l_sorted = .false.
        end do
        do j = 2, gn
            if( gpinds(j) <= gpinds(j-1) ) l_sorted = .false.
        end do
        j = 1
        do i = 1, my_nptcls
            hit = 0
            if( l_sorted )then
                do while( j <= gn )
                    if( gpinds(j) >= my_pinds(i) ) exit
                    j = j + 1
                end do
                if( j <= gn )then
                    if( gpinds(j) == my_pinds(i) ) hit = j
                endif
            else
                do j = 1, gn
                    if( gpinds(j) == my_pinds(i) )then
                        hit = j
                        exit
                    endif
                end do
            endif
            if( hit == 0 ) THROW_HARD('flex_pca worker particle absent from the master weight table')
            weights_out(i,:) = gw(hit,:)
        end do
        deallocate(gpinds, gw)
        call fname%kill
    end subroutine read_state_weights_round

end module simple_flex_pca_rec3D
