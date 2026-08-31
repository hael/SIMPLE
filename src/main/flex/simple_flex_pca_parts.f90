!@descr: versioned part-file contracts for distributed flex_pca stage reductions
module simple_flex_pca_parts
use simple_core_module_api
use simple_parameters, only: parameters
implicit none
private
#include "simple_local_flags.inc"

public :: flex_pca_part_fname
public :: write_mean_scale, read_mean_scale
public :: write_sigma_state, check_sigma_state
public :: write_state_weights_round, read_state_weights_round
public :: write_embed_stats_part, reduce_embed_zhalf_parts, read_embed_stats_part
public :: write_probe_part,  reduce_probe_parts

character(len=*), parameter :: MEAN_SCALE_FNAME = 'flex_pca_mean_scale.bin'
character(len=*), parameter :: SIGMA_STATE_FNAME= 'flex_pca_sigma_state.txt'
character(len=*), parameter :: WEIGHTS_FNAME    = 'flex_pca_round_weights.bin'
integer, parameter :: EMBED_STATS_VERSION = 1
integer, parameter :: PROBE_PART_VERSION  = 4   ! v4: rho rows are always the full packed triangle (no diagonal option)

! Every part file is magic + version + shape header + payload, written to a .tmp and renamed, so a
! master that finds the final name is guaranteed a complete file. The master reduces STREAMING --
! one part resident at a time, folded in and deleted -- so its peak memory is one accumulator plus
! one buffer regardless of nparts. Same shape as write_projected_latent_mstep_stats.
integer, parameter :: FLEX_PCA_PART_MAGIC = 1180053590

contains

    function flex_pca_part_fname( prefix, part, numlen ) result( fname )
        character(len=*), intent(in) :: prefix
        integer,          intent(in) :: part, numlen
        type(string) :: fname
        fname = string('flex_pca_')//prefix//'_part'//int2str_pad(part, numlen)//'.bin'
    end function flex_pca_part_fname

    ! ---- master -> worker handoffs ----

    !> The mean's radial amplitude scale. A worker MUST NOT re-fit this: estimate_mean_scale derives
    !! stride = max(1, nptcls/NSAMPLE) from the particle count it is given, so a worker holding a
    !! fraction of the particles would sample a different subset and fit a different scale. Shipping
    !! the nyq-length filter instead of the scaled volume keeps the handoff exact and tiny -- the
    !! worker rebuilds the mean deterministically from vol1 and applies the same array.
    subroutine write_mean_scale( nyq, filt )
        integer, intent(in) :: nyq
        real,    intent(in) :: filt(nyq)
        type(string) :: fname, tmp_fname
        integer :: funit, io_stat
        fname     = string(MEAN_SCALE_FNAME)
        tmp_fname = fname//'.tmp'
        call fopen(funit, file=tmp_fname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk('write_mean_scale; open', io_stat)
        write(funit, iostat=io_stat) FLEX_PCA_PART_MAGIC, nyq
        call fileiochk('write_mean_scale; header', io_stat)
        write(funit, iostat=io_stat) filt
        call fileiochk('write_mean_scale; payload', io_stat)
        call fclose(funit)
        call simple_rename(tmp_fname, fname)
        call fname%kill; call tmp_fname%kill
    end subroutine write_mean_scale

    subroutine read_mean_scale( nyq, filt, ok )
        integer,           intent(in)  :: nyq
        real,              intent(out) :: filt(nyq)
        logical,           intent(out) :: ok
        type(string) :: fname
        integer :: funit, io_stat, magic, nyq_in
        ok    = .false.
        fname = string(MEAN_SCALE_FNAME)
        if( .not. file_exists(fname) )then
            call fname%kill
            return
        endif
        call fopen(funit, file=fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
        call fileiochk('read_mean_scale; open', io_stat)
        read(funit, iostat=io_stat) magic, nyq_in
        call fileiochk('read_mean_scale; header', io_stat)
        if( magic /= FLEX_PCA_PART_MAGIC ) THROW_HARD('bad mean-scale magic')
        if( nyq_in /= nyq ) THROW_HARD('mean-scale band mismatch; master and worker disagree on box_crop')
        read(funit, iostat=io_stat) filt
        call fileiochk('read_mean_scale; payload', io_stat)
        call fclose(funit)
        ok = .true.
        call fname%kill
    end subroutine read_mean_scale

    !> Whether the master resolved real sigma2 spectra or fell back to the unit white-noise spectrum.
    !! Discovery is re-run independently in every worker, and a worker that resolves it differently
    !! whitens against a different noise model -- which changes every inner product without any error
    !! being raised, because both outcomes are individually legitimate. Pin it and make workers assert.
    subroutine write_sigma_state( loaded )
        logical, intent(in) :: loaded
        integer :: funit, io_stat
        call fopen(funit, file=string(SIGMA_STATE_FNAME), action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk('write_sigma_state; open', io_stat)
        write(funit,'(L1)') loaded
        call fclose(funit)
    end subroutine write_sigma_state

    subroutine check_sigma_state( loaded )
        logical, intent(in) :: loaded
        integer :: funit, io_stat
        logical :: loaded_master
        if( .not. file_exists(string(SIGMA_STATE_FNAME)) )then
            THROW_HARD('flex_pca worker found no flex_pca_sigma_state.txt from the master')
        endif
        call fopen(funit, file=string(SIGMA_STATE_FNAME), action='READ', status='OLD', iostat=io_stat)
        call fileiochk('check_sigma_state; open', io_stat)
        read(funit,'(L1)') loaded_master
        call fclose(funit)
        if( loaded_master .neqv. loaded )then
            THROW_HARD('flex_pca worker resolved sigma2 differently from the master; the parts would be &
                &whitened against different noise models')
        endif
    end subroutine check_sigma_state

    !> Per-round state weights. Written with the GLOBAL pinds so a worker can match its own particle
    !! rows by index rather than assuming its part is a contiguous slice of the master's selection.
    !! Rewritten before every state round, because the combined/even/odd and CV passes each use a
    !! different weight table.
    subroutine write_state_weights_round( pinds, weights, nptcls, nstates, split_eo )
        integer, intent(in) :: pinds(:), nptcls, nstates
        real,    intent(in) :: weights(:,:)
        !! .true. when this round accumulates the even and odd halfsets into SEPARATE reconstructors.
        !! A worker cannot infer this from params%stage -- combined/even-odd and CV rounds all arrive
        !! as PCA_STAGE_STATES -- and without it the master's odd halfset reduces to zero.
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

    !> Worker side: return the weight rows for THIS part's pinds, in this part's order.
    subroutine read_state_weights_round( my_pinds, my_nptcls, weights_out, nstates, split_eo )
        integer,           intent(in)  :: my_pinds(:), my_nptcls
        real, allocatable, intent(out) :: weights_out(:,:)
        integer,           intent(out) :: nstates
        logical,           intent(out) :: split_eo !< see write_state_weights_round
        type(string) :: fname
        integer, allocatable :: gpinds(:)
        real,    allocatable :: gw(:,:)
        integer :: funit, io_stat, magic, gn, i, j, hit, eo_flag
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
        do i = 1, my_nptcls
            hit = 0
            do j = 1, gn
                if( gpinds(j) == my_pinds(i) )then
                    hit = j
                    exit
                endif
            end do
            if( hit == 0 ) THROW_HARD('flex_pca worker particle absent from the master weight table')
            weights_out(i,:) = gw(hit,:)
        end do
        deallocate(gpinds, gw)
        call fname%kill
    end subroutine read_state_weights_round

    ! ---- probe stage: a SUM of every accumulator one EM half-pass produces ----

    !> One probe iteration's accumulators from this part's particle range.
    !!
    !! Payload, in order: for each of the ncomp components the even and odd Y_q reconstructor
    !! (cmat then rho), then the two COUPLED normal-matrix arrays rho_e / rho_o -- one entry per
    !! (q,r) pair rather than one shared density, because the M-step solves the components together
    !! at every grid point -- then the EM Gamma numerator and the valid-particle count.
    !! Gamma is shipped as a SUM, never a mean: the master divides by the reduced nval, and dividing
    !! per part would weight small parts equally with large ones.
    !> v3: the four MCFA accumulators are pure ADDITIVE sufficient statistics (already summed
    !! over threads by the caller), so distributing them is the same operation as the existing
    !! thread reduction. header(9) carries kmix, 0 when absent, which keeps v3 readable by the
    !! version check without a second format.
    subroutine write_probe_part( fname, cmat_e, rho_ex, cmat_o, rho_ox, rho_e, rho_o, &
        &gam_sum, nll_sum, nval, ncomp, mix_sr, mix_sm, mix_smm, mix_sainv, z_sub )
        class(string), intent(in) :: fname
        complex,       intent(in) :: cmat_e(:,:,:,:), cmat_o(:,:,:,:)
        real,          intent(in) :: rho_ex(:,:,:,:), rho_ox(:,:,:,:)
        real,          intent(in) :: rho_e(:,:,:,:),  rho_o(:,:,:,:)
        real(dp),      intent(in) :: gam_sum(:)
        !> sum over this part's particles of log det A_i - h_i' A_i^-1 h_i, the particle-dependent
        !! part of -2 log p(y_i). Shipped as a SUM for the same reason gam_sum is.
        real(dp),      intent(in) :: nll_sum
        integer,       intent(in) :: nval, ncomp
        real(dp), optional, intent(in) :: mix_sr(:), mix_sm(:,:), mix_smm(:,:,:), mix_sainv(:,:)
        !> bounded latent subsample: mcfa_init seeds k-means from z, which no single rank holds
        !! under sharding, so each part ships a deterministic slice and the master pools them.
        real(dp), optional, intent(in) :: z_sub(:,:)
        type(string) :: tmp_fname
        integer :: funit, io_stat, header(9), q, kmix_w
        kmix_w = 0
        if( present(mix_sr) ) kmix_w = size(mix_sr)
        header = [FLEX_PCA_PART_MAGIC, PROBE_PART_VERSION, ncomp, &
            &size(cmat_e,1), size(cmat_e,2), size(cmat_e,3), size(rho_e,1), nval, kmix_w]
        tmp_fname = fname//'.tmp'
        call fopen(funit, file=tmp_fname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk('write_probe_part; open', io_stat)
        write(funit, iostat=io_stat) header
        call fileiochk('write_probe_part; header', io_stat)
        do q = 1, ncomp
            write(funit, iostat=io_stat) cmat_e(:,:,:,q), rho_ex(:,:,:,q), cmat_o(:,:,:,q), rho_ox(:,:,:,q)
            call fileiochk('write_probe_part; basis payload', io_stat)
        end do
        write(funit, iostat=io_stat) rho_e, rho_o, gam_sum, nll_sum
        call fileiochk('write_probe_part; coupled payload', io_stat)
        if( kmix_w > 0 )then
            write(funit, iostat=io_stat) mix_sr, mix_sm, mix_smm, mix_sainv
            call fileiochk('write_probe_part; mixture payload', io_stat)
            if( present(z_sub) )then
                write(funit, iostat=io_stat) size(z_sub,1)
                write(funit, iostat=io_stat) z_sub
            else
                write(funit, iostat=io_stat) 0
            endif
            call fileiochk('write_probe_part; z subsample', io_stat)
        endif
        call fclose(funit)
        call simple_rename(tmp_fname, fname)
        call tmp_fname%kill
    end subroutine write_probe_part

    !> Streaming sum of every part into the master's accumulators; one part resident at a time.
    subroutine reduce_probe_parts( params, nparts, cmat_e, rho_ex, cmat_o, rho_ox, rho_e, rho_o, &
        &gam_sum, nll_sum, nval, ncomp, mix_sr, mix_sm, mix_smm, mix_sainv, z_pool, nz_pool )
        class(parameters), intent(in)    :: params
        integer,           intent(in)    :: nparts, ncomp
        complex,           intent(inout) :: cmat_e(:,:,:,:), cmat_o(:,:,:,:)
        real,              intent(inout) :: rho_ex(:,:,:,:), rho_ox(:,:,:,:)
        real,              intent(inout) :: rho_e(:,:,:,:),  rho_o(:,:,:,:)
        real(dp),          intent(inout) :: gam_sum(:)
        real(dp),          intent(inout) :: nll_sum
        integer,           intent(inout) :: nval
        real(dp), optional, intent(inout) :: mix_sr(:), mix_sm(:,:), mix_smm(:,:,:), mix_sainv(:,:)
        real(dp), optional, intent(inout) :: z_pool(:,:)
        integer,  optional, intent(inout) :: nz_pool
        real(dp), allocatable :: zbuf(:,:)
        integer :: nz_part
        real(dp), allocatable :: msr(:), msm(:,:), msmm(:,:,:), msai(:,:)
        complex, allocatable :: cbuf(:,:,:)
        real,    allocatable :: rbuf(:,:,:), rcbuf(:,:,:,:)
        real(dp), allocatable :: gbuf(:)
        real(dp) :: nllbuf
        type(string) :: fname
        integer :: ipart, funit, io_stat, header(9), q
        integer(timer_int_kind) :: t_red
        t_red = tic()
        allocate(cbuf(size(cmat_e,1),size(cmat_e,2),size(cmat_e,3)))
        allocate(rbuf(size(rho_ex,1),size(rho_ex,2),size(rho_ex,3)))
        allocate(rcbuf(size(rho_e,1),size(rho_e,2),size(rho_e,3),size(rho_e,4)))
        allocate(gbuf(size(gam_sum)))
        if( present(mix_sr) ) allocate(msr(size(mix_sr)), msm(size(mix_sm,1),size(mix_sm,2)), &
            &msmm(size(mix_smm,1),size(mix_smm,2),size(mix_smm,3)), &
            &msai(size(mix_sainv,1),size(mix_sainv,2)))
        do ipart = 1, nparts
            fname = flex_pca_part_fname('probe', ipart, params%numlen)
            if( .not. file_exists(fname) ) THROW_HARD('missing probe part: '//fname%to_char())
            call fopen(funit, file=fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
            call fileiochk('reduce_probe_parts; open '//fname%to_char(), io_stat)
            read(funit, iostat=io_stat) header
            call fileiochk('reduce_probe_parts; header', io_stat)
            if( header(1) /= FLEX_PCA_PART_MAGIC ) THROW_HARD('bad probe part magic')
            if( header(2) /= PROBE_PART_VERSION  ) THROW_HARD('bad probe part version')
            if( header(3) /= ncomp               ) THROW_HARD('probe part ncomp mismatch')
            if( header(7) /= size(rho_e,1)       ) THROW_HARD('probe part npairs mismatch')
            do q = 1, ncomp
                read(funit, iostat=io_stat) cbuf
                call fileiochk('reduce_probe_parts; cmat_e', io_stat)
                cmat_e(:,:,:,q) = cmat_e(:,:,:,q) + cbuf
                read(funit, iostat=io_stat) rbuf
                call fileiochk('reduce_probe_parts; rho_ex', io_stat)
                rho_ex(:,:,:,q) = rho_ex(:,:,:,q) + rbuf
                read(funit, iostat=io_stat) cbuf
                call fileiochk('reduce_probe_parts; cmat_o', io_stat)
                cmat_o(:,:,:,q) = cmat_o(:,:,:,q) + cbuf
                read(funit, iostat=io_stat) rbuf
                call fileiochk('reduce_probe_parts; rho_ox', io_stat)
                rho_ox(:,:,:,q) = rho_ox(:,:,:,q) + rbuf
            end do
            read(funit, iostat=io_stat) rcbuf
            call fileiochk('reduce_probe_parts; rho_e', io_stat)
            rho_e = rho_e + rcbuf
            read(funit, iostat=io_stat) rcbuf
            call fileiochk('reduce_probe_parts; rho_o', io_stat)
            rho_o = rho_o + rcbuf
            read(funit, iostat=io_stat) gbuf
            call fileiochk('reduce_probe_parts; gamma', io_stat)
            gam_sum = gam_sum + gbuf
            read(funit, iostat=io_stat) nllbuf
            call fileiochk('reduce_probe_parts; loglik', io_stat)
            nll_sum = nll_sum + nllbuf
            nval    = nval + header(8)
            if( present(mix_sr) )then
                if( header(9) /= size(mix_sr) ) THROW_HARD('probe part kmix mismatch')
                read(funit, iostat=io_stat) msr, msm, msmm, msai
                call fileiochk('reduce_probe_parts; mixture', io_stat)
                mix_sr    = mix_sr    + msr
                mix_sm    = mix_sm    + msm
                mix_smm   = mix_smm   + msmm
                mix_sainv = mix_sainv + msai
                read(funit, iostat=io_stat) nz_part
                call fileiochk('reduce_probe_parts; z subsample count', io_stat)
                if( nz_part > 0 )then
                    allocate(zbuf(nz_part, size(mix_sm,1)))
                    read(funit, iostat=io_stat) zbuf
                    call fileiochk('reduce_probe_parts; z subsample', io_stat)
                    if( present(z_pool) .and. present(nz_pool) )then
                        do q = 1, nz_part
                            if( nz_pool >= size(z_pool,1) ) exit
                            nz_pool = nz_pool + 1
                            z_pool(nz_pool,:) = zbuf(q,:)
                        end do
                    endif
                    deallocate(zbuf)
                endif
            endif
            call fclose(funit)
            call del_file(fname)
            call fname%kill
        end do
        deallocate(cbuf, rbuf, rcbuf, gbuf)
        write(logfhandle,'(A,I0,A,I0,A,F8.1)') '>>> FLEX_PCA reduced probe parts=',nparts, &
            &'  valid particles=',nval,'  seconds=',toc(t_red)
        call flush(logfhandle)
    end subroutine reduce_probe_parts

    ! ---- embedding stage: a PARTITION of the per-particle rows, not a sum ----

    !> Index of key in an ascending array, 0 if absent. pinds arrive in project order, so the scatter
    !! below can binary search rather than scan linearly, which is O(nptcls) per row.
    pure integer function binsrch_int( arr, n, key ) result( pos )
        integer, intent(in) :: n, arr(n), key
        integer :: lo, hi, mid
        pos = 0
        lo  = 1
        hi  = n
        do while( lo <= hi )
            mid = (lo + hi)/2
            if( arr(mid) == key )then
                pos = mid
                return
            else if( arr(mid) < key )then
                lo = mid + 1
            else
                hi = mid - 1
            endif
        end do
    end function binsrch_int

    !> One part's embedding sufficient statistics.
    !!
    !! The embedding is not a clean partition, so a part cannot ship finished latents.
    !! The reliability prior comes from
    !! rho(q) = corr(zhalf(:,q,1), zhalf(:,q,2)) over every particle, and each particle's final z is
    !! re-solved against it; per-part rho would solve the parts against different priors.
    !!
    !! So a part ships what it can compute independently -- the per-particle sufficient statistics
    !! from the image pass plus its own rows of the split-half latents -- and the master does the
    !! coupled arithmetic: reduce zhalf, form rho and the prior once, then re-solve. The re-solve
    !! touches no images, so the stage that actually costs is the part that distributes.
    subroutine write_embed_stats_part( fname, pinds, contrast, resid_energy, resid_mean_energy, &
        &Gcache, bcache, ccache, zhalf, nptcls, ncomp )
        class(string), intent(in) :: fname
        integer,       intent(in) :: pinds(:), nptcls, ncomp
        real(dp),      intent(in) :: contrast(:), resid_energy(:), resid_mean_energy(:)
        real(dp),      intent(in) :: Gcache(:,:,:), bcache(:,:), ccache(:,:), zhalf(:,:,:)
        type(string) :: tmp_fname
        integer :: funit, io_stat, header(4)
        header = [FLEX_PCA_PART_MAGIC, EMBED_STATS_VERSION, nptcls, ncomp]
        tmp_fname = fname//'.tmp'
        call fopen(funit, file=tmp_fname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk('write_embed_stats_part; open', io_stat)
        write(funit, iostat=io_stat) header
        call fileiochk('write_embed_stats_part; header', io_stat)
        ! zhalf before Gcache, deliberately: the master forms rho from the split-half latents before
        ! it can solve anything, and consumes the much larger Gram blocks one part at a time. Small
        ! arrays first lets pass 1 stop reading at zhalf.
        write(funit, iostat=io_stat) pinds(1:nptcls)
        write(funit, iostat=io_stat) contrast(1:nptcls), resid_energy(1:nptcls), resid_mean_energy(1:nptcls)
        write(funit, iostat=io_stat) zhalf(1:nptcls,1:ncomp,1:2)
        write(funit, iostat=io_stat) Gcache(1:ncomp,1:ncomp,1:nptcls)
        write(funit, iostat=io_stat) bcache(1:ncomp,1:nptcls), ccache(1:ncomp,1:nptcls)
        call fileiochk('write_embed_stats_part; payload', io_stat)
        call fclose(funit)
        call simple_rename(tmp_fname, fname)
        call tmp_fname%kill
    end subroutine write_embed_stats_part

    !> Pass 1: gather what the master needs before it can solve anything -- the split-half latents
    !! that form rho, plus the per-particle scalars. Reads no Gram blocks and deletes nothing; the
    !! files are consumed by read_embed_stats_part below.
    !!
    !! Splitting the reduce in two keeps the master's footprint flat in dataset size: holding every
    !! part's Gcache at once costs ncomp^2 doubles per particle, one part at a time costs that
    !! divided by nparts.
    subroutine reduce_embed_zhalf_parts( params, nparts, gpinds, contrast, resid_energy, &
        &resid_mean_energy, zhalf, nptcls, ncomp )
        class(parameters), intent(in)    :: params
        integer,           intent(in)    :: nparts, gpinds(:), nptcls, ncomp
        real(dp),          intent(inout) :: contrast(:), resid_energy(:), resid_mean_energy(:)
        real(dp),          intent(inout) :: zhalf(:,:,:)
        integer,  allocatable :: ppinds(:)
        real(dp), allocatable :: pc(:), pre(:), prme(:), pzh(:,:,:)
        type(string) :: fname
        integer :: ipart, funit, io_stat, header(4), pn, i, hit, nfilled
        integer(timer_int_kind) :: t_red
        t_red   = tic()
        nfilled = 0
        do ipart = 1, nparts
            fname = flex_pca_part_fname('embedstats', ipart, params%numlen)
            if( .not. file_exists(fname) ) THROW_HARD('missing embed-stats part: '//fname%to_char())
            call fopen(funit, file=fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
            call fileiochk('reduce_embed_zhalf_parts; open '//fname%to_char(), io_stat)
            read(funit, iostat=io_stat) header
            call fileiochk('reduce_embed_zhalf_parts; header', io_stat)
            if( header(1) /= FLEX_PCA_PART_MAGIC ) THROW_HARD('bad embed-stats part magic')
            if( header(2) /= EMBED_STATS_VERSION ) THROW_HARD('bad embed-stats part version')
            if( header(4) /= ncomp               ) THROW_HARD('embed-stats part ncomp mismatch')
            pn = header(3)
            allocate(ppinds(pn), pc(pn), pre(pn), prme(pn), pzh(pn,ncomp,2))
            read(funit, iostat=io_stat) ppinds
            read(funit, iostat=io_stat) pc, pre, prme
            read(funit, iostat=io_stat) pzh
            call fileiochk('reduce_embed_zhalf_parts; payload', io_stat)
            call fclose(funit)
            ! match on pinds rather than assuming a contiguous layout, so a part boundary that does
            ! not line up cannot silently misplace rows
            do i = 1, pn
                hit = binsrch_int(gpinds, nptcls, ppinds(i))
                if( hit < 1 ) THROW_HARD('embed-stats part carries a particle not in the global set')
                contrast(hit)          = pc(i)
                resid_energy(hit)      = pre(i)
                resid_mean_energy(hit) = prme(i)
                zhalf(hit,:,:)         = pzh(i,:,:)
                nfilled = nfilled + 1
            end do
            deallocate(ppinds, pc, pre, prme, pzh)
            call fname%kill
        end do
        if( nfilled /= nptcls ) THROW_HARD('embed-stats parts did not cover every particle')
        write(logfhandle,'(A,I0,A,I0,A,F8.1)') '>>> FLEX_PCA reduced embed-stats (zhalf) parts=',nparts, &
            &'  particles=',nfilled,'  seconds=',toc(t_red)
        call flush(logfhandle)
    end subroutine reduce_embed_zhalf_parts

    !> Pass 2: one part's Gram blocks and its global row indices, so the caller can re-solve just
    !! those particles and free the buffer before reading the next. Deletes the part file.
    subroutine read_embed_stats_part( params, ipart, gpinds, rows, Gpart, bpart, cpart, pn, nptcls, ncomp )
        class(parameters),     intent(in)  :: params
        integer,               intent(in)  :: ipart, gpinds(:), nptcls, ncomp
        integer,  allocatable, intent(out) :: rows(:)
        real(dp), allocatable, intent(out) :: Gpart(:,:,:), bpart(:,:), cpart(:,:)
        integer,               intent(out) :: pn
        integer,  allocatable :: ppinds(:)
        real(dp), allocatable :: skip3(:), pzh(:,:,:)
        type(string) :: fname
        integer :: funit, io_stat, header(4), i, hit
        fname = flex_pca_part_fname('embedstats', ipart, params%numlen)
        if( .not. file_exists(fname) ) THROW_HARD('missing embed-stats part: '//fname%to_char())
        call fopen(funit, file=fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
        call fileiochk('read_embed_stats_part; open '//fname%to_char(), io_stat)
        read(funit, iostat=io_stat) header
        if( header(1) /= FLEX_PCA_PART_MAGIC ) THROW_HARD('bad embed-stats part magic')
        if( header(2) /= EMBED_STATS_VERSION ) THROW_HARD('bad embed-stats part version')
        if( header(4) /= ncomp               ) THROW_HARD('embed-stats part ncomp mismatch')
        pn = header(3)
        allocate(ppinds(pn), skip3(3*pn), pzh(pn,ncomp,2))
        allocate(Gpart(ncomp,ncomp,pn), bpart(ncomp,pn), cpart(ncomp,pn), rows(pn))
        read(funit, iostat=io_stat) ppinds
        read(funit, iostat=io_stat) skip3          ! contrast, resid_energy, resid_mean_energy: pass 1
        read(funit, iostat=io_stat) pzh            ! zhalf: pass 1
        read(funit, iostat=io_stat) Gpart
        read(funit, iostat=io_stat) bpart, cpart
        call fileiochk('read_embed_stats_part; payload', io_stat)
        call fclose(funit)
        do i = 1, pn
            hit = binsrch_int(gpinds, nptcls, ppinds(i))
            if( hit < 1 ) THROW_HARD('embed-stats part carries a particle not in the global set')
            rows(i) = hit
        end do
        deallocate(ppinds, skip3, pzh)
        call del_file(fname)
        call fname%kill
    end subroutine read_embed_stats_part

end module simple_flex_pca_parts
