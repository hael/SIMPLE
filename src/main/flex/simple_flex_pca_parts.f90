!@descr: versioned part-file contracts for distributed flex_pca stage reductions
module simple_flex_pca_parts
use simple_core_module_api
use simple_parameters, only: parameters
implicit none
private
#include "simple_local_flags.inc"

public :: flex_pca_part_fname
public :: write_snr_part,   reduce_snr_parts
public :: write_cols_part,  reduce_cols_parts
public :: write_solve_part, reduce_solve_parts
public :: write_mean_scale, read_mean_scale
public :: write_columns_hkl, read_columns_hkl
public :: write_sigma_state, check_sigma_state
public :: write_state_weights_round, read_state_weights_round
public :: write_states_part, reduce_states_parts
public :: write_embed_part,  reduce_embed_parts
public :: write_embed_stats_part, reduce_embed_zhalf_parts, read_embed_stats_part
public :: write_probe_part,  reduce_probe_parts
public :: cleanup_flex_pca_parts

character(len=*), parameter :: MEAN_SCALE_FNAME = 'flex_pca_mean_scale.bin'
character(len=*), parameter :: COLUMNS_FNAME    = 'flex_pca_columns.bin'
character(len=*), parameter :: SIGMA_STATE_FNAME= 'flex_pca_sigma_state.txt'
character(len=*), parameter :: WEIGHTS_FNAME    = 'flex_pca_round_weights.bin'
integer, parameter :: STATES_PART_VERSION = 1
integer, parameter :: EMBED_PART_VERSION  = 1
integer, parameter :: EMBED_STATS_VERSION = 1
integer, parameter :: PROBE_PART_VERSION  = 1

! Every part file is magic + version + shape header + payload, written to a .tmp and renamed, so a
! master that finds the final name is guaranteed a complete file. The master reduces STREAMING --
! one part resident at a time, folded in and deleted -- so its peak memory is one accumulator plus
! one buffer regardless of nparts. Same shape as write_projected_latent_mstep_stats.
integer, parameter :: FLEX_PCA_PART_MAGIC = 1180053590
integer, parameter :: SNR_PART_VERSION    = 1
integer, parameter :: COLS_PART_VERSION   = 1
integer, parameter :: SOLVE_PART_VERSION  = 1

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

    !> Selected covariance columns, so workers do not re-run the SNR-greedy selection (which would
    !! need the reduced SNR volume anyway).
    subroutine write_columns_hkl( col_hkl, ncol )
        integer, intent(in) :: col_hkl(:,:), ncol
        type(string) :: fname, tmp_fname
        integer :: funit, io_stat
        fname     = string(COLUMNS_FNAME)
        tmp_fname = fname//'.tmp'
        call fopen(funit, file=tmp_fname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk('write_columns_hkl; open', io_stat)
        write(funit, iostat=io_stat) FLEX_PCA_PART_MAGIC, ncol
        call fileiochk('write_columns_hkl; header', io_stat)
        write(funit, iostat=io_stat) col_hkl(1:3,1:ncol)
        call fileiochk('write_columns_hkl; payload', io_stat)
        call fclose(funit)
        call simple_rename(tmp_fname, fname)
        call fname%kill; call tmp_fname%kill
    end subroutine write_columns_hkl

    subroutine read_columns_hkl( col_hkl, ncol, ok )
        integer, allocatable, intent(out) :: col_hkl(:,:)
        integer,              intent(out) :: ncol
        logical,              intent(out) :: ok
        type(string) :: fname
        integer :: funit, io_stat, magic
        ok   = .false.
        ncol = 0
        fname = string(COLUMNS_FNAME)
        if( .not. file_exists(fname) )then
            call fname%kill
            return
        endif
        call fopen(funit, file=fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
        call fileiochk('read_columns_hkl; open', io_stat)
        read(funit, iostat=io_stat) magic, ncol
        call fileiochk('read_columns_hkl; header', io_stat)
        if( magic /= FLEX_PCA_PART_MAGIC ) THROW_HARD('bad column-list magic')
        if( ncol < 1 ) THROW_HARD('empty column list')
        allocate(col_hkl(3,ncol))
        read(funit, iostat=io_stat) col_hkl
        call fileiochk('read_columns_hkl; payload', io_stat)
        call fclose(funit)
        ok = .true.
        call fname%kill
    end subroutine read_columns_hkl

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

    ! ---- state reconstruction stage: nstates COMPRESSED reconstructors (cmat + rho) ----

    subroutine write_states_part( fname, cmats, rhos, nstates )
        class(string), intent(in) :: fname
        complex,       intent(in) :: cmats(:,:,:,:)
        real,          intent(in) :: rhos(:,:,:,:)
        integer,       intent(in) :: nstates
        type(string) :: tmp_fname
        integer :: funit, io_stat, header(6), s
        header = [FLEX_PCA_PART_MAGIC, STATES_PART_VERSION, nstates, &
            &size(cmats,1), size(cmats,2), size(cmats,3)]
        tmp_fname = fname//'.tmp'
        call fopen(funit, file=tmp_fname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk('write_states_part; open', io_stat)
        write(funit, iostat=io_stat) header
        call fileiochk('write_states_part; header', io_stat)
        do s = 1, nstates
            write(funit, iostat=io_stat) cmats(:,:,:,s), rhos(:,:,:,s)
            call fileiochk('write_states_part; payload', io_stat)
        end do
        call fclose(funit)
        call simple_rename(tmp_fname, fname)
        call tmp_fname%kill
    end subroutine write_states_part

    subroutine reduce_states_parts( params, nparts, cmats, rhos, nstates )
        class(parameters), intent(in)    :: params
        integer,           intent(in)    :: nparts, nstates
        complex,           intent(inout) :: cmats(:,:,:,:)
        real,              intent(inout) :: rhos(:,:,:,:)
        complex, allocatable :: cbuf(:,:,:)
        real,    allocatable :: rbuf(:,:,:)
        type(string) :: fname
        integer :: ipart, funit, io_stat, header(6), s
        integer(timer_int_kind) :: t_red
        t_red = tic()
        allocate(cbuf(size(cmats,1),size(cmats,2),size(cmats,3)))
        allocate(rbuf(size(rhos,1), size(rhos,2), size(rhos,3)))
        do ipart = 1, nparts
            fname = flex_pca_part_fname('states', ipart, params%numlen)
            if( .not. file_exists(fname) ) THROW_HARD('missing states part: '//fname%to_char())
            call fopen(funit, file=fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
            call fileiochk('reduce_states_parts; open '//fname%to_char(), io_stat)
            read(funit, iostat=io_stat) header
            call fileiochk('reduce_states_parts; header', io_stat)
            if( header(1) /= FLEX_PCA_PART_MAGIC )  THROW_HARD('bad states part magic')
            if( header(2) /= STATES_PART_VERSION )  THROW_HARD('bad states part version')
            if( header(3) /= nstates )              THROW_HARD('states part nstates mismatch')
            do s = 1, nstates
                read(funit, iostat=io_stat) cbuf, rbuf
                call fileiochk('reduce_states_parts; payload', io_stat)
                cmats(:,:,:,s) = cmats(:,:,:,s) + cbuf
                rhos(:,:,:,s)  = rhos(:,:,:,s)  + rbuf
            end do
            call fclose(funit)
            call del_file(fname)
            call fname%kill
        end do
        deallocate(cbuf, rbuf)
        write(logfhandle,'(A,I0,A,F8.1)') '>>> FLEX_PCA reduced states parts=',nparts,' seconds=',toc(t_red)
        call flush(logfhandle)
    end subroutine reduce_states_parts


    ! ---- probe stage: a SUM of every accumulator one EM half-pass produces ----

    !> One probe iteration's accumulators from this part's particle range.
    !!
    !! Payload, in order: for each of the ncomp components the even and odd Y_q reconstructor
    !! (cmat then rho), then the two COUPLED normal-matrix arrays rho_e / rho_o -- one entry per
    !! (q,r) pair rather than one shared density, because the M-step solves the components together
    !! at every grid point -- then the EM Gamma numerator and the valid-particle count.
    !! Gamma is shipped as a SUM, never a mean: the master divides by the reduced nval, and dividing
    !! per part would weight small parts equally with large ones.
    subroutine write_probe_part( fname, cmat_e, rho_ex, cmat_o, rho_ox, rho_e, rho_o, &
        &gam_sum, nval, ncomp )
        class(string), intent(in) :: fname
        complex,       intent(in) :: cmat_e(:,:,:,:), cmat_o(:,:,:,:)
        real,          intent(in) :: rho_ex(:,:,:,:), rho_ox(:,:,:,:)
        real,          intent(in) :: rho_e(:,:,:,:),  rho_o(:,:,:,:)
        real(dp),      intent(in) :: gam_sum(:)
        integer,       intent(in) :: nval, ncomp
        type(string) :: tmp_fname
        integer :: funit, io_stat, header(9), q
        header = [FLEX_PCA_PART_MAGIC, PROBE_PART_VERSION, ncomp, &
            &size(cmat_e,1), size(cmat_e,2), size(cmat_e,3), size(rho_e,1), nval, 0]
        tmp_fname = fname//'.tmp'
        call fopen(funit, file=tmp_fname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk('write_probe_part; open', io_stat)
        write(funit, iostat=io_stat) header
        call fileiochk('write_probe_part; header', io_stat)
        do q = 1, ncomp
            write(funit, iostat=io_stat) cmat_e(:,:,:,q), rho_ex(:,:,:,q), cmat_o(:,:,:,q), rho_ox(:,:,:,q)
            call fileiochk('write_probe_part; basis payload', io_stat)
        end do
        write(funit, iostat=io_stat) rho_e, rho_o, gam_sum
        call fileiochk('write_probe_part; coupled payload', io_stat)
        call fclose(funit)
        call simple_rename(tmp_fname, fname)
        call tmp_fname%kill
    end subroutine write_probe_part

    !> Streaming sum of every part into the master's accumulators; one part resident at a time.
    subroutine reduce_probe_parts( params, nparts, cmat_e, rho_ex, cmat_o, rho_ox, rho_e, rho_o, &
        &gam_sum, nval, ncomp )
        class(parameters), intent(in)    :: params
        integer,           intent(in)    :: nparts, ncomp
        complex,           intent(inout) :: cmat_e(:,:,:,:), cmat_o(:,:,:,:)
        real,              intent(inout) :: rho_ex(:,:,:,:), rho_ox(:,:,:,:)
        real,              intent(inout) :: rho_e(:,:,:,:),  rho_o(:,:,:,:)
        real(dp),          intent(inout) :: gam_sum(:)
        integer,           intent(inout) :: nval
        complex, allocatable :: cbuf(:,:,:)
        real,    allocatable :: rbuf(:,:,:), rcbuf(:,:,:,:)
        real(dp), allocatable :: gbuf(:)
        type(string) :: fname
        integer :: ipart, funit, io_stat, header(9), q
        integer(timer_int_kind) :: t_red
        t_red = tic()
        allocate(cbuf(size(cmat_e,1),size(cmat_e,2),size(cmat_e,3)))
        allocate(rbuf(size(rho_ex,1),size(rho_ex,2),size(rho_ex,3)))
        allocate(rcbuf(size(rho_e,1),size(rho_e,2),size(rho_e,3),size(rho_e,4)))
        allocate(gbuf(size(gam_sum)))
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
            nval    = nval + header(8)
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

    subroutine write_embed_part( fname, pinds, z, contrast, precision, resid_energy, &
        &resid_mean_energy, nptcls, ncomp )
        class(string), intent(in) :: fname
        integer,       intent(in) :: pinds(:), nptcls, ncomp
        real(dp),      intent(in) :: z(:,:), contrast(:), precision(:,:,:)
        real(dp),      intent(in) :: resid_energy(:), resid_mean_energy(:)
        type(string) :: tmp_fname
        integer :: funit, io_stat, header(4)
        header = [FLEX_PCA_PART_MAGIC, EMBED_PART_VERSION, nptcls, ncomp]
        tmp_fname = fname//'.tmp'
        call fopen(funit, file=tmp_fname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk('write_embed_part; open', io_stat)
        write(funit, iostat=io_stat) header
        call fileiochk('write_embed_part; header', io_stat)
        write(funit, iostat=io_stat) pinds(1:nptcls)
        write(funit, iostat=io_stat) z(1:nptcls,1:ncomp)
        write(funit, iostat=io_stat) contrast(1:nptcls)
        write(funit, iostat=io_stat) precision(1:ncomp,1:ncomp,1:nptcls)
        write(funit, iostat=io_stat) resid_energy(1:nptcls), resid_mean_energy(1:nptcls)
        call fileiochk('write_embed_part; payload', io_stat)
        call fclose(funit)
        call simple_rename(tmp_fname, fname)
        call tmp_fname%kill
    end subroutine write_embed_part

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
    !! The embedding is not a clean partition, so a part cannot ship finished latents (which is why
    !! write_embed_part above was never wired up). The reliability prior comes from
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

    !> Scatter each part's rows into the global arrays by matching pinds. The embedding is a partition
    !! of the particles, so this is a gather, not a sum -- and matching on pinds rather than assuming a
    !! contiguous layout means a part boundary that does not line up cannot silently misplace rows.
    subroutine reduce_embed_parts( params, nparts, gpinds, z, contrast, precision, resid_energy, &
        &resid_mean_energy, nptcls, ncomp )
        class(parameters), intent(in)    :: params
        integer,           intent(in)    :: nparts, gpinds(:), nptcls, ncomp
        real(dp),          intent(inout) :: z(:,:), contrast(:), precision(:,:,:)
        real(dp),          intent(inout) :: resid_energy(:), resid_mean_energy(:)
        integer,  allocatable :: ppinds(:)
        real(dp), allocatable :: pz(:,:), pc(:), pp(:,:,:), pre(:), prme(:)
        type(string) :: fname
        integer :: ipart, funit, io_stat, header(4), pn, i, j, hit, nfilled
        integer(timer_int_kind) :: t_red
        t_red   = tic()
        nfilled = 0
        do ipart = 1, nparts
            fname = flex_pca_part_fname('embed', ipart, params%numlen)
            if( .not. file_exists(fname) ) THROW_HARD('missing embed part: '//fname%to_char())
            call fopen(funit, file=fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
            call fileiochk('reduce_embed_parts; open '//fname%to_char(), io_stat)
            read(funit, iostat=io_stat) header
            call fileiochk('reduce_embed_parts; header', io_stat)
            if( header(1) /= FLEX_PCA_PART_MAGIC ) THROW_HARD('bad embed part magic')
            if( header(2) /= EMBED_PART_VERSION )  THROW_HARD('bad embed part version')
            if( header(4) /= ncomp )               THROW_HARD('embed part component count mismatch')
            pn = header(3)
            allocate(ppinds(pn), pz(pn,ncomp), pc(pn), pp(ncomp,ncomp,pn), pre(pn), prme(pn))
            read(funit, iostat=io_stat) ppinds
            read(funit, iostat=io_stat) pz
            read(funit, iostat=io_stat) pc
            read(funit, iostat=io_stat) pp
            read(funit, iostat=io_stat) pre, prme
            call fileiochk('reduce_embed_parts; payload', io_stat)
            call fclose(funit)
            do i = 1, pn
                hit = 0
                do j = 1, nptcls
                    if( gpinds(j) == ppinds(i) )then
                        hit = j
                        exit
                    endif
                end do
                if( hit == 0 ) THROW_HARD('embed part carries a particle absent from the global selection')
                z(hit,:)             = pz(i,:)
                contrast(hit)        = pc(i)
                precision(:,:,hit)   = pp(:,:,i)
                resid_energy(hit)    = pre(i)
                resid_mean_energy(hit) = prme(i)
                nfilled = nfilled + 1
            end do
            deallocate(ppinds, pz, pc, pp, pre, prme)
            call del_file(fname)
            call fname%kill
        end do
        if( nfilled /= nptcls )then
            THROW_HARD('embed reduction filled '//int2str(nfilled)//' of '//int2str(nptcls)//' particles')
        endif
        write(logfhandle,'(A,I0,A,I0,A,F8.1)') '>>> FLEX_PCA reduced embed parts=',nparts, &
            &' particles=',nfilled,' seconds=',toc(t_red)
        call flush(logfhandle)
    end subroutine reduce_embed_parts

    ! ---- SNR variance stage: two real accumulators over the expanded lattice ----

    subroutine write_snr_part( fname, var_acc, dens_acc, lb, ub )
        class(string), intent(in) :: fname
        real,          intent(in) :: var_acc(:,:,:), dens_acc(:,:,:)
        integer,       intent(in) :: lb(3), ub(3)
        type(string) :: tmp_fname
        integer :: funit, io_stat, header(9)
        header = [FLEX_PCA_PART_MAGIC, SNR_PART_VERSION, lb(1),lb(2),lb(3), ub(1),ub(2),ub(3), 0]
        tmp_fname = fname//'.tmp'
        call fopen(funit, file=tmp_fname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk('write_snr_part; open '//tmp_fname%to_char(), io_stat)
        write(funit, iostat=io_stat) header
        call fileiochk('write_snr_part; header', io_stat)
        write(funit, iostat=io_stat) var_acc
        call fileiochk('write_snr_part; var_acc', io_stat)
        write(funit, iostat=io_stat) dens_acc
        call fileiochk('write_snr_part; dens_acc', io_stat)
        call fclose(funit)
        call simple_rename(tmp_fname, fname)
        call tmp_fname%kill
    end subroutine write_snr_part

    !> Streaming sum over parts. var_acc/dens_acc must arrive zeroed; each part is deleted once folded in.
    subroutine reduce_snr_parts( params, nparts, var_acc, dens_acc, lb, ub )
        class(parameters), intent(in)    :: params
        integer,           intent(in)    :: nparts, lb(3), ub(3)
        real,              intent(inout) :: var_acc(:,:,:), dens_acc(:,:,:)
        real, allocatable :: vbuf(:,:,:), dbuf(:,:,:)
        type(string) :: fname
        integer :: ipart, funit, io_stat, header(9)
        integer(timer_int_kind) :: t_red
        t_red = tic()
        allocate(vbuf(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)))
        allocate(dbuf(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)))
        do ipart = 1, nparts
            fname = flex_pca_part_fname('snr', ipart, params%numlen)
            if( .not. file_exists(fname) ) THROW_HARD('missing SNR part: '//fname%to_char())
            call fopen(funit, file=fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
            call fileiochk('reduce_snr_parts; open '//fname%to_char(), io_stat)
            read(funit, iostat=io_stat) header
            call fileiochk('reduce_snr_parts; header', io_stat)
            if( header(1) /= FLEX_PCA_PART_MAGIC )  THROW_HARD('bad SNR part magic: '//fname%to_char())
            if( header(2) /= SNR_PART_VERSION )     THROW_HARD('bad SNR part version: '//fname%to_char())
            if( any(header(3:5) /= lb) .or. any(header(6:8) /= ub) )then
                THROW_HARD('SNR part lattice mismatch: '//fname%to_char())
            endif
            read(funit, iostat=io_stat) vbuf
            call fileiochk('reduce_snr_parts; var_acc', io_stat)
            read(funit, iostat=io_stat) dbuf
            call fileiochk('reduce_snr_parts; dens_acc', io_stat)
            call fclose(funit)
            var_acc  = var_acc  + vbuf
            dens_acc = dens_acc + dbuf
            call del_file(fname)
            call fname%kill
        end do
        deallocate(vbuf, dbuf)
        write(logfhandle,'(A,I0,A,F8.1)') '>>> FLEX_PCA reduced SNR parts=',nparts,' seconds=',toc(t_red)
        call flush(logfhandle)
    end subroutine reduce_snr_parts

    ! ---- Column accumulation stage: two complex and two real accumulators, ncol deep ----

    subroutine write_cols_part( fname, Bcol_e, Hcol_e, Bcol_o, Hcol_o, lb, ub, ncol )
        class(string), intent(in) :: fname
        complex,       intent(in) :: Bcol_e(:,:,:,:), Bcol_o(:,:,:,:)
        real,          intent(in) :: Hcol_e(:,:,:,:), Hcol_o(:,:,:,:)
        integer,       intent(in) :: lb(3), ub(3), ncol
        type(string) :: tmp_fname
        integer :: funit, io_stat, header(9), s
        header = [FLEX_PCA_PART_MAGIC, COLS_PART_VERSION, lb(1),lb(2),lb(3), ub(1),ub(2),ub(3), ncol]
        tmp_fname = fname//'.tmp'
        call fopen(funit, file=tmp_fname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk('write_cols_part; open '//tmp_fname%to_char(), io_stat)
        write(funit, iostat=io_stat) header
        call fileiochk('write_cols_part; header', io_stat)
        ! written column by column so the master can stream a column at a time if it ever needs to
        do s = 1, ncol
            write(funit, iostat=io_stat) Bcol_e(:,:,:,s), Hcol_e(:,:,:,s), Bcol_o(:,:,:,s), Hcol_o(:,:,:,s)
            call fileiochk('write_cols_part; column payload', io_stat)
        end do
        call fclose(funit)
        call simple_rename(tmp_fname, fname)
        call tmp_fname%kill
    end subroutine write_cols_part

    subroutine reduce_cols_parts( params, nparts, Bcol_e, Hcol_e, Bcol_o, Hcol_o, lb, ub, ncol )
        class(parameters), intent(in)    :: params
        integer,           intent(in)    :: nparts, lb(3), ub(3), ncol
        complex,           intent(inout) :: Bcol_e(:,:,:,:), Bcol_o(:,:,:,:)
        real,              intent(inout) :: Hcol_e(:,:,:,:), Hcol_o(:,:,:,:)
        complex, allocatable :: cbuf_e(:,:,:), cbuf_o(:,:,:)
        real,    allocatable :: rbuf_e(:,:,:), rbuf_o(:,:,:)
        type(string) :: fname
        integer :: ipart, funit, io_stat, header(9), s
        integer(timer_int_kind) :: t_red
        t_red = tic()
        allocate(cbuf_e(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)), cbuf_o(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)))
        allocate(rbuf_e(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)), rbuf_o(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)))
        do ipart = 1, nparts
            fname = flex_pca_part_fname('cols', ipart, params%numlen)
            if( .not. file_exists(fname) ) THROW_HARD('missing column part: '//fname%to_char())
            call fopen(funit, file=fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
            call fileiochk('reduce_cols_parts; open '//fname%to_char(), io_stat)
            read(funit, iostat=io_stat) header
            call fileiochk('reduce_cols_parts; header', io_stat)
            if( header(1) /= FLEX_PCA_PART_MAGIC ) THROW_HARD('bad column part magic: '//fname%to_char())
            if( header(2) /= COLS_PART_VERSION )   THROW_HARD('bad column part version: '//fname%to_char())
            if( any(header(3:5) /= lb) .or. any(header(6:8) /= ub) .or. header(9) /= ncol )then
                THROW_HARD('column part shape mismatch: '//fname%to_char())
            endif
            do s = 1, ncol
                read(funit, iostat=io_stat) cbuf_e, rbuf_e, cbuf_o, rbuf_o
                call fileiochk('reduce_cols_parts; column payload', io_stat)
                Bcol_e(:,:,:,s) = Bcol_e(:,:,:,s) + cbuf_e
                Hcol_e(:,:,:,s) = Hcol_e(:,:,:,s) + rbuf_e
                Bcol_o(:,:,:,s) = Bcol_o(:,:,:,s) + cbuf_o
                Hcol_o(:,:,:,s) = Hcol_o(:,:,:,s) + rbuf_o
            end do
            call fclose(funit)
            call del_file(fname)
            call fname%kill
        end do
        deallocate(cbuf_e, cbuf_o, rbuf_e, rbuf_o)
        write(logfhandle,'(A,I0,A,F8.1)') '>>> FLEX_PCA reduced column parts=',nparts,' seconds=',toc(t_red)
        call flush(logfhandle)
    end subroutine reduce_cols_parts

    ! ---- reduced-solve stage ----
    ! Payload is dominated by the rank-k accumulator: Mspk(npk,npk) under the default packed/CG path
    ! (npk = d(d+1)/2 = 8256 at d_tilde=128, 545 MB) or A(dd,dd) under SIMPLE_COV_CGSOLVE=0 (2.1 GB).
    ! Only the UPPER triangle is populated by dsyrk in each part; summing upper triangles and mirroring
    ! once on the master is equivalent to mirroring per part and summing.
    subroutine write_solve_part( fname, acc, nacc, Sbb, SG, d_tilde, pwsh, cntsh, nyq, scal, nvalid )
        class(string), intent(in) :: fname
        integer,       intent(in) :: nacc, d_tilde, nyq, nvalid
        real(dp),      intent(in) :: acc(nacc,nacc), Sbb(d_tilde,d_tilde), SG(d_tilde,d_tilde)
        real(dp),      intent(in) :: pwsh(0:nyq), cntsh(0:nyq), scal(8)
        type(string) :: tmp_fname
        integer :: funit, io_stat, header(6)
        header = [FLEX_PCA_PART_MAGIC, SOLVE_PART_VERSION, nacc, d_tilde, nyq, nvalid]
        tmp_fname = fname//'.tmp'
        call fopen(funit, file=tmp_fname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk('write_solve_part; open '//tmp_fname%to_char(), io_stat)
        write(funit, iostat=io_stat) header
        call fileiochk('write_solve_part; header', io_stat)
        write(funit, iostat=io_stat) acc
        call fileiochk('write_solve_part; accumulator', io_stat)
        write(funit, iostat=io_stat) Sbb, SG, pwsh, cntsh, scal
        call fileiochk('write_solve_part; moments', io_stat)
        call fclose(funit)
        call simple_rename(tmp_fname, fname)
        call tmp_fname%kill
    end subroutine write_solve_part

    subroutine reduce_solve_parts( params, nparts, acc, nacc, Sbb, SG, d_tilde, pwsh, cntsh, nyq, scal, nvalid )
        class(parameters), intent(in)    :: params
        integer,           intent(in)    :: nparts, nacc, d_tilde, nyq
        real(dp),          intent(inout) :: acc(nacc,nacc), Sbb(d_tilde,d_tilde), SG(d_tilde,d_tilde)
        real(dp),          intent(inout) :: pwsh(0:nyq), cntsh(0:nyq), scal(8)
        integer,           intent(out)   :: nvalid
        real(dp), allocatable :: abuf(:,:), sbuf(:,:), gbuf2(:,:), pbuf(:), cbuf(:), scbuf(:)
        type(string) :: fname
        integer :: ipart, funit, io_stat, header(6)
        integer(timer_int_kind) :: t_red
        t_red  = tic()
        nvalid = 0
        allocate(abuf(nacc,nacc), sbuf(d_tilde,d_tilde), gbuf2(d_tilde,d_tilde))
        allocate(pbuf(0:nyq), cbuf(0:nyq), scbuf(8))
        do ipart = 1, nparts
            fname = flex_pca_part_fname('solve', ipart, params%numlen)
            if( .not. file_exists(fname) ) THROW_HARD('missing solve part: '//fname%to_char())
            call fopen(funit, file=fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
            call fileiochk('reduce_solve_parts; open '//fname%to_char(), io_stat)
            read(funit, iostat=io_stat) header
            call fileiochk('reduce_solve_parts; header', io_stat)
            if( header(1) /= FLEX_PCA_PART_MAGIC ) THROW_HARD('bad solve part magic: '//fname%to_char())
            if( header(2) /= SOLVE_PART_VERSION )  THROW_HARD('bad solve part version: '//fname%to_char())
            if( header(3) /= nacc .or. header(4) /= d_tilde .or. header(5) /= nyq )then
                THROW_HARD('solve part shape mismatch: '//fname%to_char())
            endif
            read(funit, iostat=io_stat) abuf
            call fileiochk('reduce_solve_parts; accumulator', io_stat)
            read(funit, iostat=io_stat) sbuf, gbuf2, pbuf, cbuf, scbuf
            call fileiochk('reduce_solve_parts; moments', io_stat)
            call fclose(funit)
            acc    = acc   + abuf
            Sbb    = Sbb   + sbuf
            SG     = SG    + gbuf2
            pwsh   = pwsh  + pbuf
            cntsh  = cntsh + cbuf
            scal   = scal  + scbuf
            nvalid = nvalid + header(6)
            call del_file(fname)
            call fname%kill
        end do
        deallocate(abuf, sbuf, gbuf2, pbuf, cbuf, scbuf)
        write(logfhandle,'(A,I0,A,I0,A,F8.1)') '>>> FLEX_PCA reduced solve parts=',nparts, &
            &' particles=',nvalid,' seconds=',toc(t_red)
        call flush(logfhandle)
    end subroutine reduce_solve_parts

    subroutine cleanup_flex_pca_parts( nparts, numlen )
        integer, intent(in) :: nparts, numlen
        type(string) :: fname
        integer :: ipart
        character(len=5), parameter :: PREFIXES(3) = ['snr  ', 'cols ', 'solve']
        integer :: ip
        do ip = 1, size(PREFIXES)
            do ipart = 1, nparts
                fname = flex_pca_part_fname(trim(PREFIXES(ip)), ipart, numlen)
                call del_file(fname)
                call del_file(fname//'.tmp')
                call fname%kill
            end do
        end do
    end subroutine cleanup_flex_pca_parts

end module simple_flex_pca_parts
