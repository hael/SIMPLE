!@descr: flex_pca: the fit engine: one master loop for single, paired and worker fits, the distributed round, the probe part codec and its reduces
submodule (simple_flex_probe_fit) simple_flex_probe_fit_engine
use simple_core_module_api
use simple_flex_pca_records, only: flex_fit_model, flex_selection
use simple_builder, only: builder
use simple_parameters, only: parameters
use simple_reconstructor, only: reconstructor
use simple_flex_pca_rounds, only: flex_pca_rounds
use simple_flex_pca_stages, only: flex_stage_request, PCA_STAGE_PROBE, PCA_STAGE_POLISH
use simple_flex_pca_artifacts, only: flex_pca_part_fname, FLEX_PCA_PART_MAGIC
use simple_flex_pca_run_types, only: flex_run_settings
use simple_flex_pca_posterior, only: mcfa_init
use simple_flex_pca_basis, only: save_probe_state, COV_PROBE_META
use simple_flex_pca_util, only: cov_stage_subsample
use simple_flex_pca_fit_types, only: flex_probe_part, probe_part_borrow, probe_part_restore, xfsc_ctx_t
use simple_flex_pca_plane_cache, only: plane_cache_in_use
implicit none
#include "simple_local_flags.inc"

contains

    !> Master-side v5 reduce: fold every worker's paired part into BOTH fits' accumulators,
    !! streaming -- one part resident at a time, per-fit blocks in fit order, file deleted after.
    module subroutine paired_reduce_parts_v5( params, fits , rounds)
        class(flex_pca_rounds), intent(inout) :: rounds
        class(parameters), intent(in)    :: params
        type(flex_probe_fit),    intent(inout) :: fits(2)
        !> keep equal to the engine's MIX_ZSUB_MAX (per-part latent subsample, em_iter)
        integer, parameter :: MIX_ZSUB_MAX5 = 2000
        type(flex_probe_part) :: acc(2)
        type(string) :: fname
        integer :: ipart, funit, f, q
        integer(timer_int_kind) :: t_red
        t_red = tic()
        do f = 1, 2
            associate( es => fits(f)%mstep%es, nc => fits(f)%model%ncomp )
            allocate(acc(f)%cmat_e(es(1),es(2),es(3),nc), acc(f)%cmat_o(es(1),es(2),es(3),nc), source=(0.,0.))
            allocate(acc(f)%rho_ex(es(1),es(2),es(3),nc), acc(f)%rho_ox(es(1),es(2),es(3),nc), source=0.)
            end associate
            if( fits(f)%spec%l_mix_req )then
                ! a rank change between iterations leaves these at the old ncomp; the part files carry
                ! the new one, and buffers sized from stale arrays misalign every read that follows
                if( allocated(fits(f)%history%dm_sm) )then
                    if( size(fits(f)%history%dm_sm,1) /= fits(f)%model%ncomp .or. size(fits(f)%history%dm_sr) /= fits(f)%spec%kmix )then
                        deallocate(fits(f)%history%dm_sr, fits(f)%history%dm_sm, fits(f)%history%dm_smm, fits(f)%history%dm_sai, fits(f)%history%dm_z)
                    endif
                endif
                if( .not. allocated(fits(f)%history%dm_sr) )then
                    allocate(fits(f)%history%dm_sr(fits(f)%spec%kmix), fits(f)%history%dm_sm(fits(f)%model%ncomp,fits(f)%spec%kmix), &
                        &fits(f)%history%dm_smm(fits(f)%model%ncomp,fits(f)%model%ncomp,fits(f)%spec%kmix), &
                        &fits(f)%history%dm_sai(fits(f)%model%ncomp,fits(f)%model%ncomp))
                    allocate(fits(f)%history%dm_z(MIX_ZSUB_MAX5*rounds%nparts(), fits(f)%model%ncomp))
                endif
                fits(f)%history%dm_sr = 0.d0; fits(f)%history%dm_sm = 0.d0; fits(f)%history%dm_smm = 0.d0
                fits(f)%history%dm_sai = 0.d0; fits(f)%history%dm_nz = 0
            endif
            ! the fit's own accumulators fold in place: borrowed, never copied
            call probe_part_borrow(acc(f), fits(f), l_mix_buffers=fits(f)%spec%l_mix_req)
        end do
        do ipart = 1, rounds%nparts()
            fname = flex_pca_part_fname('probe', ipart, params%numlen)
            call open_probe_part_v5_read(fname, 2, funit)
            do f = 1, 2
                call fold_probe_part_v5_fit(funit, acc(f))
            end do
            call close_probe_part_v5_read(funit, fname)
            call fname%kill
        end do
        do f = 1, 2
            call probe_part_restore(acc(f), fits(f))
            do q = 1, fits(f)%model%ncomp
                fits(f)%mstep%Yeven(q)%cmat_exp = acc(f)%cmat_e(:,:,:,q); fits(f)%mstep%Yeven(q)%rho_exp = acc(f)%rho_ex(:,:,:,q)
                fits(f)%mstep%Yodd(q)%cmat_exp  = acc(f)%cmat_o(:,:,:,q); fits(f)%mstep%Yodd(q)%rho_exp  = acc(f)%rho_ox(:,:,:,q)
            end do
            call acc(f)%kill
        end do
        write(logfhandle,'(A,I0,A,I0,A,I0,A,F8.1)') '>>> FLEX_PCA reduced v5 paired parts=', &
            &rounds%nparts(), '  valid A=', fits(1)%iter%nval, '  valid B=', fits(2)%iter%nval, &
            &'  seconds=', toc(t_red)
        call flush(logfhandle)
    end subroutine paired_reduce_parts_v5

    !> One probe iteration's accumulators from this part's particle range. Payload, in order: per
    !! component the even and odd Y_q reconstructor (cmat then rho), the coupled normal-matrix arrays
    !! rho_e / rho_o (one entry per (q,r) pair), the EM Gamma numerator and the valid-particle count.
    !! Gamma is shipped as a sum (the master divides by the reduced nval). The MCFA accumulators are
    !! additive sufficient statistics; header(9) carries kmix, 0 when absent.
    !> One probe iteration's accumulators from this part's particle range. Payload, in order: per
    !! component the even and odd Y_q reconstructor (cmat then rho), the coupled normal-matrix arrays
    !! rho_e / rho_o (one entry per (q,r) pair), the EM Gamma numerator and the valid-particle count.
    !! Gamma is shipped as a sum (the master divides by the reduced nval). The MCFA accumulators are
    !! additive sufficient statistics; header(9) carries kmix, 0 when absent.
    module subroutine write_probe_part( fname, part )
        class(string),         intent(in) :: fname
        type(flex_probe_part), intent(in) :: part
        type(string) :: tmp_fname
        integer :: funit, io_stat, header(9), q, kmix_w
        kmix_w = 0
        if( allocated(part%mix_sr) ) kmix_w = size(part%mix_sr)
        header = [FLEX_PCA_PART_MAGIC, PROBE_PART_VERSION, part%ncomp, &
            &size(part%cmat_e,1), size(part%cmat_e,2), size(part%cmat_e,3), size(part%rho_e,1), part%nval, kmix_w]
        tmp_fname = fname//'.tmp'
        call fopen(funit, file=tmp_fname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk('write_probe_part; open', io_stat)
        write(funit, iostat=io_stat) header
        call fileiochk('write_probe_part; header', io_stat)
        do q = 1, part%ncomp
            write(funit, iostat=io_stat) part%cmat_e(:,:,:,q), part%rho_ex(:,:,:,q), part%cmat_o(:,:,:,q), part%rho_ox(:,:,:,q)
            call fileiochk('write_probe_part; basis payload', io_stat)
        end do
        write(funit, iostat=io_stat) part%rho_e, part%rho_o, part%gam_sum, part%nll_sum
        call fileiochk('write_probe_part; coupled payload', io_stat)
        if( kmix_w > 0 )then
            write(funit, iostat=io_stat) part%mix_sr, part%mix_sm, part%mix_smm, part%mix_sainv
            call fileiochk('write_probe_part; mixture payload', io_stat)
            if( allocated(part%z_sub) )then
                write(funit, iostat=io_stat) size(part%z_sub,1)
                write(funit, iostat=io_stat) part%z_sub
            else
                write(funit, iostat=io_stat) 0
            endif
            call fileiochk('write_probe_part; z subsample', io_stat)
        endif
        call write_probe_part_kernels(funit, part%kpk_e, part%kpk_o, part%rpk_e, part%rpk_o, 'write_probe_part')
        call fclose(funit)
        call simple_rename(tmp_fname, fname)
        call tmp_fname%kill
    end subroutine write_probe_part

    !> Streaming sum of every part into the master's accumulators; one part resident at a time.
    !> Streaming sum of every part into the master's accumulators (the borrowed fit arrays inside
    !! `part`); one part resident at a time.
    module subroutine reduce_probe_parts( params, nparts, part )
        class(parameters),     intent(in)    :: params
        integer,               intent(in)    :: nparts
        type(flex_probe_part), intent(inout) :: part
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
        allocate(cbuf(size(part%cmat_e,1),size(part%cmat_e,2),size(part%cmat_e,3)))
        allocate(rbuf(size(part%rho_ex,1),size(part%rho_ex,2),size(part%rho_ex,3)))
        allocate(rcbuf(size(part%rho_e,1),size(part%rho_e,2),size(part%rho_e,3),size(part%rho_e,4)))
        allocate(gbuf(size(part%gam_sum)))
        if( allocated(part%mix_sr) ) allocate(msr(size(part%mix_sr)), msm(size(part%mix_sm,1),size(part%mix_sm,2)), &
            &msmm(size(part%mix_smm,1),size(part%mix_smm,2),size(part%mix_smm,3)), &
            &msai(size(part%mix_sainv,1),size(part%mix_sainv,2)))
        do ipart = 1, nparts
            fname = flex_pca_part_fname('probe', ipart, params%numlen)
            if( .not. file_exists(fname) ) THROW_HARD('missing probe part: '//fname%to_char())
            call fopen(funit, file=fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
            call fileiochk('reduce_probe_parts; open '//fname%to_char(), io_stat)
            read(funit, iostat=io_stat) header
            call fileiochk('reduce_probe_parts; header', io_stat)
            if( header(1) /= FLEX_PCA_PART_MAGIC ) THROW_HARD('bad probe part magic')
            if( header(2) /= PROBE_PART_VERSION  ) THROW_HARD('bad probe part version')
            if( header(3) /= part%ncomp               ) THROW_HARD('probe part ncomp mismatch')
            if( header(7) /= size(part%rho_e,1)       ) THROW_HARD('probe part npairs mismatch')
            do q = 1, part%ncomp
                read(funit, iostat=io_stat) cbuf
                call fileiochk('reduce_probe_parts; cmat_e', io_stat)
                !$omp parallel workshare
                part%cmat_e(:,:,:,q) = part%cmat_e(:,:,:,q) + cbuf
                !$omp end parallel workshare
                read(funit, iostat=io_stat) rbuf
                call fileiochk('reduce_probe_parts; rho_ex', io_stat)
                !$omp parallel workshare
                part%rho_ex(:,:,:,q) = part%rho_ex(:,:,:,q) + rbuf
                !$omp end parallel workshare
                read(funit, iostat=io_stat) cbuf
                call fileiochk('reduce_probe_parts; cmat_o', io_stat)
                !$omp parallel workshare
                part%cmat_o(:,:,:,q) = part%cmat_o(:,:,:,q) + cbuf
                !$omp end parallel workshare
                read(funit, iostat=io_stat) rbuf
                call fileiochk('reduce_probe_parts; rho_ox', io_stat)
                !$omp parallel workshare
                part%rho_ox(:,:,:,q) = part%rho_ox(:,:,:,q) + rbuf
                !$omp end parallel workshare
            end do
            read(funit, iostat=io_stat) rcbuf
            call fileiochk('reduce_probe_parts; rho_e', io_stat)
            !$omp parallel workshare
            part%rho_e = part%rho_e + rcbuf
            !$omp end parallel workshare
            read(funit, iostat=io_stat) rcbuf
            call fileiochk('reduce_probe_parts; rho_o', io_stat)
            !$omp parallel workshare
            part%rho_o = part%rho_o + rcbuf
            !$omp end parallel workshare
            read(funit, iostat=io_stat) gbuf
            call fileiochk('reduce_probe_parts; gamma', io_stat)
            part%gam_sum = part%gam_sum + gbuf
            read(funit, iostat=io_stat) nllbuf
            call fileiochk('reduce_probe_parts; loglik', io_stat)
            part%nll_sum = part%nll_sum + nllbuf
            part%nval    = part%nval + header(8)
            if( allocated(part%mix_sr) )then
                if( header(9) /= size(part%mix_sr) ) THROW_HARD('probe part kmix mismatch')
                read(funit, iostat=io_stat) msr, msm, msmm, msai
                call fileiochk('reduce_probe_parts; mixture', io_stat)
                part%mix_sr    = part%mix_sr    + msr
                part%mix_sm    = part%mix_sm    + msm
                part%mix_smm   = part%mix_smm   + msmm
                part%mix_sainv = part%mix_sainv + msai
                read(funit, iostat=io_stat) nz_part
            call fileiochk('fold_probe_part_v5_fit; latent subsample count', io_stat)
            if( nz_part < 0 .or. nz_part > 1000000 ) THROW_HARD('v5 probe part reader: corrupt latent subsample count (stream misaligned)')
                call fileiochk('reduce_probe_parts; z subsample count', io_stat)
                if( nz_part > 0 )then
                    allocate(zbuf(nz_part, size(part%mix_sm,1)))
                    read(funit, iostat=io_stat) zbuf
                    call fileiochk('reduce_probe_parts; z subsample', io_stat)
                    if( allocated(part%z_sub) )then
                        do q = 1, nz_part
                            if( part%nz >= size(part%z_sub,1) ) exit
                            part%nz = part%nz + 1
                            part%z_sub(part%nz,:) = zbuf(q,:)
                        end do
                    endif
                    deallocate(zbuf)
                endif
            endif
            call fold_probe_part_kernels(funit, part%kpk_e, part%kpk_o, part%rpk_e, part%rpk_o, 'reduce_probe_parts')
            call fclose(funit)
            call del_file(fname)
            call fname%kill
        end do
        deallocate(cbuf, rbuf, rcbuf, gbuf)
        write(logfhandle,'(A,I0,A,I0,A,F8.1)') '>>> FLEX_PCA reduced probe parts=',nparts, &
            &'  valid particles=',part%nval,'  seconds=',toc(t_red)
        call flush(logfhandle)
    end subroutine reduce_probe_parts

    !> Open one v5 part for writing and emit the file header. The caller then writes one
    !! block per fit with write_probe_part_v5_fit, IN FIT ORDER, and closes with
    !! close_probe_part_v5_write (which renames the .tmp so the master only ever sees
    !! complete files -- same contract as write_probe_part).
    module subroutine open_probe_part_v5_write( fname, nfits, funit, tmp_fname )
        class(string), intent(in)  :: fname
        integer,       intent(in)  :: nfits
        integer,       intent(out) :: funit
        type(string),  intent(out) :: tmp_fname
        integer :: io_stat
        tmp_fname = fname//'.tmp'
        call fopen(funit, file=tmp_fname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk('open_probe_part_v5_write; open', io_stat)
        write(funit, iostat=io_stat) FLEX_PCA_PART_MAGIC, PROBE_PART_VERSION5, nfits
        call fileiochk('open_probe_part_v5_write; header', io_stat)
    end subroutine open_probe_part_v5_write

    !> One fit's block: sub-header + the write_probe_part payload.
    !> One fit's block: sub-header + the write_probe_part payload.
    module subroutine write_probe_part_v5_fit( funit, part )
        integer,               intent(in) :: funit
        type(flex_probe_part), intent(in) :: part
        integer :: io_stat, subhdr(7), q, kmix_w, t, pnnz
        logical, allocatable :: pm(:,:,:)
        integer, allocatable :: pii(:), pjj(:), pkk(:)
        complex, allocatable :: gce(:), gco(:)
        real,    allocatable :: gre(:), gro(:), gpe(:,:), gpo(:,:)
        kmix_w = 0
        if( allocated(part%mix_sr) ) kmix_w = size(part%mix_sr)
        subhdr = [part%ncomp, size(part%cmat_e,1), size(part%cmat_e,2), size(part%cmat_e,3), size(part%rho_e,1), part%nval, kmix_w]
        write(funit, iostat=io_stat) subhdr
        call fileiochk('write_probe_part_v5_fit; sub-header', io_stat)
        ! the coupled rho rows must live on the SAME lattice as the basis accumulators, otherwise
        ! one band box cannot describe both and the boxed payload would be mis-scattered
        if( size(part%rho_e,2) /= size(part%cmat_e,1) .or. size(part%rho_e,3) /= size(part%cmat_e,2) .or. &
            &size(part%rho_e,4) /= size(part%cmat_e,3) ) &
            &THROW_HARD('write_probe_part_v5_fit: rho_e lattice differs from the basis lattice')
        ! one nonzero bounding box for the whole crop lattice (see "band boxing" above)
        ! index list of every populated lattice point, unioned over the six shipped arrays
        allocate(pm(size(part%cmat_e,1),size(part%cmat_e,2),size(part%cmat_e,3)), source=.false.)
        call pk_mask_c1(part%cmat_e, pm)
        call pk_mask_r1(part%rho_ex, pm)
        call pk_mask_c1(part%cmat_o, pm)
        call pk_mask_r1(part%rho_ox, pm)
        call pk_mask_r2(part%rho_e,  pm)
        call pk_mask_r2(part%rho_o,  pm)
        call pk_mask_to_idx(pm, pii, pjj, pkk, pnnz)
        deallocate(pm)
        write(funit, iostat=io_stat) pnnz
        call fileiochk('write_probe_part_v5_fit; index count', io_stat)
        write(funit, iostat=io_stat) pii, pjj, pkk
        call fileiochk('write_probe_part_v5_fit; index list', io_stat)
        allocate(gce(pnnz), gco(pnnz), gre(pnnz), gro(pnnz))
        do q = 1, part%ncomp
            do t = 1, pnnz
                gce(t) = part%cmat_e(pii(t),pjj(t),pkk(t),q)
                gre(t) = part%rho_ex(pii(t),pjj(t),pkk(t),q)
                gco(t) = part%cmat_o(pii(t),pjj(t),pkk(t),q)
                gro(t) = part%rho_ox(pii(t),pjj(t),pkk(t),q)
            end do
            write(funit, iostat=io_stat) gce, gre, gco, gro
            call fileiochk('write_probe_part_v5_fit; basis payload', io_stat)
        end do
        deallocate(gce, gco, gre, gro)
        allocate(gpe(size(part%rho_e,1),pnnz), gpo(size(part%rho_o,1),pnnz))
        do t = 1, pnnz
            gpe(:,t) = part%rho_e(:,pii(t),pjj(t),pkk(t))
            gpo(:,t) = part%rho_o(:,pii(t),pjj(t),pkk(t))
        end do
        write(funit, iostat=io_stat) gpe, gpo, part%gam_sum, part%nll_sum
        call fileiochk('write_probe_part_v5_fit; coupled payload', io_stat)
        deallocate(gpe, gpo, pii, pjj, pkk)
        if( kmix_w > 0 )then
            write(funit, iostat=io_stat) part%mix_sr, part%mix_sm, part%mix_smm, part%mix_sainv
            call fileiochk('write_probe_part_v5_fit; mixture payload', io_stat)
            if( allocated(part%z_sub) )then
                write(funit, iostat=io_stat) size(part%z_sub,1)
                write(funit, iostat=io_stat) part%z_sub
            else
                write(funit, iostat=io_stat) 0
            endif
            call fileiochk('write_probe_part_v5_fit; z subsample', io_stat)
        endif
        call write_probe_part_kernels(funit, part%kpk_e, part%kpk_o, part%rpk_e, part%rpk_o, 'write_probe_part_v5_fit')
    end subroutine write_probe_part_v5_fit

    subroutine pk_mask_r2( a, m )
        real,    intent(in)    :: a(:,:,:,:)
        logical, intent(inout) :: m(:,:,:)
        integer :: i, j, k
        do k = 1, size(a,4)
            do j = 1, size(a,3)
                do i = 1, size(a,2)
                    if( any(a(:,i,j,k) /= 0.) ) m(i,j,k) = .true.
                end do
            end do
        end do
    end subroutine pk_mask_r2

    subroutine pk_mask_c1( a, m )
        complex, intent(in)    :: a(:,:,:,:)
        logical, intent(inout) :: m(:,:,:)
        integer :: i, j, k
        !$omp parallel do collapse(3) default(shared) private(i,j,k) schedule(static)
        do k = 1, size(a,3)
            do j = 1, size(a,2)
                do i = 1, size(a,1)
                    if( any(a(i,j,k,:) /= (0.,0.)) ) m(i,j,k) = .true.
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine pk_mask_c1

    subroutine pk_mask_r1( a, m )
        real,    intent(in)    :: a(:,:,:,:)
        logical, intent(inout) :: m(:,:,:)
        integer :: i, j, k
        !$omp parallel do collapse(3) default(shared) private(i,j,k) schedule(static)
        do k = 1, size(a,3)
            do j = 1, size(a,2)
                do i = 1, size(a,1)
                    if( any(a(i,j,k,:) /= 0.) ) m(i,j,k) = .true.
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine pk_mask_r1

    !> flatten the union mask into the three coordinate lists the payload is gathered on.
    !! An empty mask keeps ONE voxel so every extent stays positive; folding one zero is a no-op.
    subroutine pk_mask_to_idx( m, ii, jj, kk, nnz )
        logical,              intent(in)  :: m(:,:,:)
        integer, allocatable, intent(out) :: ii(:), jj(:), kk(:)
        integer,              intent(out) :: nnz
        integer :: i, j, k, t
        nnz = count(m)
        if( nnz == 0 )then
            allocate(ii(1), jj(1), kk(1)); ii = 1; jj = 1; kk = 1; nnz = 1
            return
        endif
        allocate(ii(nnz), jj(nnz), kk(nnz))
        t = 0
        do k = 1, size(m,3)
            do j = 1, size(m,2)
                do i = 1, size(m,1)
                    if( m(i,j,k) )then
                        t = t + 1; ii(t) = i; jj(t) = j; kk(t) = k
                    endif
                end do
            end do
        end do
    end subroutine pk_mask_to_idx

    !> trailing PCG kernel block of one fit: the pair count (0 when the fit does not solve by PCG), then
    !! the packed even and odd kernel sums
    !> PCG kernel and rhs payload of a part: the packed arrays on the band list, slot for slot. Every
    !! process derives the same list from (box, band, rank), so the reader adds by slot; the list
    !! length travels in the shape and is checked. (Format 12: no per-part index list.)
    subroutine write_probe_part_kernels( funit, kpk_e, kpk_o, rpk_e, rpk_o, who )
        integer,           intent(in) :: funit
        real, allocatable, intent(in) :: kpk_e(:,:), kpk_o(:,:)
        complex, allocatable, intent(in) :: rpk_e(:,:), rpk_o(:,:)
        character(len=*),  intent(in) :: who
        integer :: io_stat, nkp, nrp
        nkp = 0
        if( allocated(kpk_e) .and. allocated(kpk_o) ) nkp = size(kpk_e,1)
        write(funit, iostat=io_stat) nkp
        call fileiochk(who//'; PCG kernel count', io_stat)
        if( nkp > 0 )then
            write(funit, iostat=io_stat) shape(kpk_e)
            write(funit, iostat=io_stat) kpk_e, kpk_o
            call fileiochk(who//'; PCG kernel payload', io_stat)
        endif
        nrp = 0
        if( allocated(rpk_e) .and. allocated(rpk_o) ) nrp = size(rpk_e,1)
        write(funit, iostat=io_stat) nrp
        call fileiochk(who//'; PCG rhs count', io_stat)
        if( nrp > 0 )then
            write(funit, iostat=io_stat) shape(rpk_e)
            write(funit, iostat=io_stat) rpk_e, rpk_o
            call fileiochk(who//'; PCG rhs payload', io_stat)
        endif
    end subroutine write_probe_part_kernels

    !> read one fit's trailing PCG kernel block and ADD it into the packed sums (a part written by a
    !! non-PCG fit carries a zero count; a PCG part folded into a non-PCG reduce is a mismatch)
    subroutine fold_probe_part_kernels( funit, kpk_e, kpk_o, rpk_e, rpk_o, who )
        integer,           intent(in)    :: funit
        real, allocatable, intent(inout) :: kpk_e(:,:), kpk_o(:,:)
        complex, allocatable, intent(inout) :: rpk_e(:,:), rpk_o(:,:)
        character(len=*),  intent(in)    :: who
        integer :: io_stat, nkp, nrp, kshape(2)
        real,    allocatable :: gke(:,:), gko(:,:)
        complex, allocatable :: gqe(:,:), gqo(:,:)
        read(funit, iostat=io_stat) nkp
        call fileiochk(who//'; PCG kernel count', io_stat)
        if( nkp == 0 )then
            if( allocated(kpk_e) ) THROW_HARD(who//': part carries no PCG kernels but the reduce expects them')
        else
            if( .not. (allocated(kpk_e) .and. allocated(kpk_o)) ) &
                &THROW_HARD(who//': part carries PCG kernels the reduce did not expect')
            read(funit, iostat=io_stat) kshape
            call fileiochk(who//'; PCG kernel shape', io_stat)
            if( any(kshape /= shape(kpk_e)) ) THROW_HARD(who//': PCG kernel band-list shape mismatch (part from another lattice or band)')
            allocate(gke(kshape(1),kshape(2)), gko(kshape(1),kshape(2)))
            read(funit, iostat=io_stat) gke, gko
            call fileiochk(who//'; PCG kernel payload', io_stat)
            !$omp parallel workshare default(shared)
            kpk_e = kpk_e + gke
            kpk_o = kpk_o + gko
            !$omp end parallel workshare
            deallocate(gke, gko)
        endif
        read(funit, iostat=io_stat) nrp
        call fileiochk(who//'; PCG rhs count', io_stat)
        if( nrp == 0 )then
            if( allocated(rpk_e) ) THROW_HARD(who//': part carries no PCG right-hand sides but the reduce expects them')
            return
        endif
        if( .not. (allocated(rpk_e) .and. allocated(rpk_o)) ) &
            &THROW_HARD(who//': part carries PCG right-hand sides the reduce did not expect')
        read(funit, iostat=io_stat) kshape
        call fileiochk(who//'; PCG rhs shape', io_stat)
        if( any(kshape /= shape(rpk_e)) ) THROW_HARD(who//': PCG rhs band-list shape mismatch (part from another lattice or band)')
        allocate(gqe(kshape(1),kshape(2)), gqo(kshape(1),kshape(2)))
        read(funit, iostat=io_stat) gqe, gqo
        call fileiochk(who//'; PCG rhs payload', io_stat)
        !$omp parallel workshare default(shared)
        rpk_e = rpk_e + gqe
        rpk_o = rpk_o + gqo
        !$omp end parallel workshare
        deallocate(gqe, gqo)
    end subroutine fold_probe_part_kernels

    module subroutine close_probe_part_v5_write( funit, tmp_fname, fname )
        integer,       intent(in)    :: funit
        type(string),  intent(inout) :: tmp_fname
        class(string), intent(in)    :: fname
        call fclose(funit)
        call simple_rename(tmp_fname, fname)
        call tmp_fname%kill
    end subroutine close_probe_part_v5_write

    !> Open one v5 part for reduction and validate the file header. The caller folds one
    !! block per fit with fold_probe_part_v5_fit, IN FIT ORDER (the file is streamed).
    module subroutine open_probe_part_v5_read( fname, nfits, funit )
        class(string), intent(in)  :: fname
        integer,       intent(in)  :: nfits
        integer,       intent(out) :: funit
        integer :: io_stat, magic, ver, nfits_in
        if( .not. file_exists(fname) ) THROW_HARD('missing v5 probe part: '//fname%to_char())
        call fopen(funit, file=fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
        call fileiochk('open_probe_part_v5_read; open '//fname%to_char(), io_stat)
        read(funit, iostat=io_stat) magic, ver, nfits_in
        call fileiochk('open_probe_part_v5_read; header', io_stat)
        if( magic /= FLEX_PCA_PART_MAGIC  ) THROW_HARD('bad v5 probe part magic')
        if( ver   /= PROBE_PART_VERSION5  ) THROW_HARD('bad v5 probe part version')
        if( nfits_in /= nfits             ) THROW_HARD('v5 probe part nfits mismatch')
    end subroutine open_probe_part_v5_read

    !> Fold one fit's block from an open v5 part into that fit's accumulators; ncomp is checked
    !! against that fit's own expectation.
    !> Fold one fit's block from an open v5 part into that fit's accumulators (borrowed inside
    !! `part`); ncomp is checked against that fit's own expectation.
    module subroutine fold_probe_part_v5_fit( funit, part )
        integer,               intent(in)    :: funit
        type(flex_probe_part), intent(inout) :: part
        real(dp), allocatable :: gbuf(:), zbuf(:,:)
        real(dp), allocatable :: msr(:), msm(:,:), msmm(:,:,:), msai(:,:)
        real(dp) :: nllbuf
        integer  :: io_stat, subhdr(7), q, nz_part, t, pnnz
        integer, allocatable :: pii(:), pjj(:), pkk(:)
        complex, allocatable :: gce(:), gco(:)
        real,    allocatable :: gre(:), gro(:), gpe(:,:), gpo(:,:)
        read(funit, iostat=io_stat) subhdr
        call fileiochk('fold_probe_part_v5_fit; sub-header', io_stat)
        if( subhdr(1) /= part%ncomp          ) THROW_HARD('v5 probe part per-fit ncomp mismatch')
        if( subhdr(5) /= size(part%rho_e,1)  ) THROW_HARD('v5 probe part per-fit npairs mismatch')
        ! the lattice dims MUST be validated, not just ncomp/npairs: a part written on a different
        ! expanded lattice yields an in-range band box whose voxels land at the wrong addresses,
        ! i.e. silently misplaced mass instead of a loud failure
        if( subhdr(2) /= size(part%cmat_e,1) .or. subhdr(3) /= size(part%cmat_e,2) .or. &
            &subhdr(4) /= size(part%cmat_e,3) ) &
            &THROW_HARD('v5 probe part per-fit lattice mismatch (part written on another lattice)')
        ! index list of the populated lattice points (see "index-list packing"): the part carries
        ! only these voxels; every omitted one was an exact zero on the write side, so the fold is
        ! bit-identical to adding the full array
        read(funit, iostat=io_stat) pnnz
        call fileiochk('fold_probe_part_v5_fit; index count', io_stat)
        if( pnnz < 1 ) THROW_HARD('v5 probe part index count is not positive')
        allocate(pii(pnnz), pjj(pnnz), pkk(pnnz))
        read(funit, iostat=io_stat) pii, pjj, pkk
        call fileiochk('fold_probe_part_v5_fit; index list', io_stat)
        if( minval(pii) < 1 .or. maxval(pii) > size(part%cmat_e,1) .or. &
            &minval(pjj) < 1 .or. maxval(pjj) > size(part%cmat_e,2) .or. &
            &minval(pkk) < 1 .or. maxval(pkk) > size(part%cmat_e,3) ) &
            &THROW_HARD('v5 probe part index list outside the basis lattice')
        allocate(gce(pnnz), gco(pnnz), gre(pnnz), gro(pnnz))
        allocate(gbuf(size(part%gam_sum)))
        do q = 1, part%ncomp
            read(funit, iostat=io_stat) gce, gre, gco, gro
            call fileiochk('fold_probe_part_v5_fit; basis payload', io_stat)
            !$omp parallel do default(shared) private(t) schedule(static)
            do t = 1, pnnz
                part%cmat_e(pii(t),pjj(t),pkk(t),q) = part%cmat_e(pii(t),pjj(t),pkk(t),q) + gce(t)
                part%rho_ex(pii(t),pjj(t),pkk(t),q) = part%rho_ex(pii(t),pjj(t),pkk(t),q) + gre(t)
                part%cmat_o(pii(t),pjj(t),pkk(t),q) = part%cmat_o(pii(t),pjj(t),pkk(t),q) + gco(t)
                part%rho_ox(pii(t),pjj(t),pkk(t),q) = part%rho_ox(pii(t),pjj(t),pkk(t),q) + gro(t)
            end do
            !$omp end parallel do
        end do
        deallocate(gce, gco, gre, gro)
        allocate(gpe(size(part%rho_e,1),pnnz), gpo(size(part%rho_o,1),pnnz))
        read(funit, iostat=io_stat) gpe, gpo
        call fileiochk('fold_probe_part_v5_fit; coupled payload', io_stat)
        !$omp parallel do default(shared) private(t) schedule(static)
        do t = 1, pnnz
            part%rho_e(:,pii(t),pjj(t),pkk(t)) = part%rho_e(:,pii(t),pjj(t),pkk(t)) + gpe(:,t)
            part%rho_o(:,pii(t),pjj(t),pkk(t)) = part%rho_o(:,pii(t),pjj(t),pkk(t)) + gpo(:,t)
        end do
        !$omp end parallel do
        deallocate(gpe, gpo)
        read(funit, iostat=io_stat) gbuf
        call fileiochk('fold_probe_part_v5_fit; gamma', io_stat)
        part%gam_sum = part%gam_sum + gbuf
        read(funit, iostat=io_stat) nllbuf
        call fileiochk('fold_probe_part_v5_fit; loglik', io_stat)
        part%nll_sum = part%nll_sum + nllbuf
        part%nval    = part%nval + subhdr(6)
        if( allocated(part%mix_sr) )then
            if( subhdr(7) /= size(part%mix_sr) ) THROW_HARD('v5 probe part per-fit kmix mismatch')
            if( size(part%mix_sm,1) /= part%ncomp .or. size(part%mix_smm,1) /= part%ncomp .or. size(part%mix_sainv,1) /= part%ncomp ) &
                &THROW_HARD('v5 probe part reader: mixture accumulators not sized to ncomp (stale after a rank change)')
            allocate(msr(size(part%mix_sr)), msm(size(part%mix_sm,1),size(part%mix_sm,2)), &
                &msmm(size(part%mix_smm,1),size(part%mix_smm,2),size(part%mix_smm,3)), &
                &msai(size(part%mix_sainv,1),size(part%mix_sainv,2)))
            read(funit, iostat=io_stat) msr, msm, msmm, msai
            call fileiochk('fold_probe_part_v5_fit; mixture', io_stat)
            part%mix_sr    = part%mix_sr    + msr
            part%mix_sm    = part%mix_sm    + msm
            part%mix_smm   = part%mix_smm   + msmm
            part%mix_sainv = part%mix_sainv + msai
            deallocate(msr, msm, msmm, msai)
            read(funit, iostat=io_stat) nz_part
            call fileiochk('fold_probe_part_v5_fit; latent subsample count', io_stat)
            if( nz_part < 0 .or. nz_part > 1000000 ) THROW_HARD('v5 probe part reader: corrupt latent subsample count (stream misaligned)')
            call fileiochk('fold_probe_part_v5_fit; z subsample count', io_stat)
            if( nz_part > 0 )then
                allocate(zbuf(nz_part, size(part%mix_sm,1)))
                read(funit, iostat=io_stat) zbuf
                call fileiochk('fold_probe_part_v5_fit; z subsample', io_stat)
                if( allocated(part%z_sub) )then
                    do q = 1, nz_part
                        if( part%nz >= size(part%z_sub,1) ) exit
                        part%nz = part%nz + 1
                        part%z_sub(part%nz,:) = zbuf(q,:)
                    end do
                endif
                deallocate(zbuf)
            endif
        else
            if( subhdr(7) /= 0 ) THROW_HARD('v5 probe part carries a mixture block this fit did not expect')
        endif
        call fold_probe_part_kernels(funit, part%kpk_e, part%kpk_o, part%rpk_e, part%rpk_o, 'fold_probe_part_v5_fit')
        deallocate(gbuf, pii, pjj, pkk)
    end subroutine fold_probe_part_v5_fit

    module subroutine close_probe_part_v5_read( funit, fname )
        integer,       intent(in) :: funit
        class(string), intent(in) :: fname
        call fclose(funit)
        call del_file(fname)
    end subroutine close_probe_part_v5_read

    !> Probe-based subspace iteration: alternate a Wiener E-step (per-particle latents in the current
    !! basis) with a weighted-backprojection M-step (Y_q += sum_i z_iq * backproject(r_i)), then
    !! orthonormalize the refined probe volumes into the next basis. The single-fit entry: the
    !! caller's mean and basis handles are hoisted into one fit, iterated by the shared engine
    !! (fit_engine_iterate), and handed back.
    module subroutine probe_subspace_iteration( params, cfg, build, model, sel, niters, it_glob, niters_glob, fprefix, meta_fname, rounds)
        type(flex_fit_model), intent(inout), target :: model   !< mean in; basis, prior variances, rank refined in place
        type(flex_selection), intent(in)            :: sel
        class(flex_pca_rounds), intent(inout) :: rounds
        class(parameters),   intent(inout) :: params
        type(flex_run_settings), intent(in)    :: cfg
        type(builder),       intent(inout) :: build
        integer,             intent(in) :: niters
        !> master's global-iteration stamp, passed by a distributed probe worker whose own loop
        !! runs once per relaunch; absent (or 0) in shared-memory runs and on the master itself
        integer, optional,   intent(in)    :: it_glob, niters_glob
        !> file namespace override (merge-polish pass: never clobber fit A's delivery)
        character(len=*), optional, intent(in) :: fprefix, meta_fname
        type(flex_probe_fit), target :: fits(1)
        type(flex_mean_ref) :: means(1)
        character(len=:), allocatable :: pfx_l, meta_l
        integer :: nparts_sub, it_glob_l, niters_glob_l
        pfx_l  = 'flex_pca_pc'
        meta_l = COV_PROBE_META
        if( present(fprefix) )    pfx_l  = trim(fprefix)
        if( present(meta_fname) ) meta_l = trim(meta_fname)
        call fits(1)%new(cfg, 0, pfx_l, meta_l, sel)
        call move_alloc(model%basis_recs, fits(1)%model%basis_recs)
        call move_alloc(model%eigvals,    fits(1)%model%eigvals)
        fits(1)%model%ncomp    = model%ncomp
        fits(1)%model%sig2_eff = model%sig2_eff
        fits(1)%model%sig2     = max(model%sig2_eff, DTINY)
        means(1)%p => model%mean_rec
        ! ---- optional STRIDE subsample, for the basis refinement only ----
        ! The probe refines ncomp band-limited, FSC-regularised volumes, and every iteration costs a
        ! full pass over the data -- far more particles than that many parameters need. A stride keeps
        ! both halfsets and every state proportionally represented; fromp/top would NOT, because
        ! particles are commonly ordered by state (on Ribosembly a contiguous window selects whole
        ! states). The embedding stage that follows still uses every particle: only the basis
        ! refinement is subsampled.
        ! The stride MUST be applied within each halfset, not across the particle list. `eo` alternates
        ! strictly by particle index (0,1,0,1,...), so a plain stride of 2 selects one halfset entirely
        ! and leaves the other empty -- and every probe M-step is regularised by an even/odd FSC, which
        ! is then computed against nothing. Measured: the Wiener filter kills the basis and the run dies
        ! at the "embedding collapsed" guard. Striding per halfset keeps both populated at any stride.
        ! ---- absolute cap, not a fixed ratio ----
        ! The probe refines ncomp band-limited, FSC-regularised volumes, and that parameter count
        ! does not grow with the dataset, so the particles needed to determine it do not either. A
        ! constant stride would leave the probe scaling linearly and dominating the run; capping the
        ! count makes it O(1) in dataset size.
        !
        ! COV_PROBE_MAX_PTCLS is the total across all processes, so each takes its share -- a worker
        ! sees only its own partition and would otherwise take the whole budget nparts times over.
        ! Only a WORKER divides the total by nparts: it holds one fromp/top partition. The master
        ! holds every particle, so passing nparts there divides twice and inflates the stride by
        ! exactly nparts. See cov_stage_subsample, which the initialiser shares.
        nparts_sub = 1
        if( rounds%is_worker() ) nparts_sub = params%nparts
        call cov_stage_subsample(build, fits(1)%spec%sel%pinds, fits(1)%spec%sel%nptcls, nparts_sub, cfg%probe_max, 'PROBE', &
            &fits(1)%spec%ppinds, fits(1)%spec%npp)
        allocate(fits(1)%model%z(fits(1)%spec%npp,fits(1)%model%ncomp))
        it_glob_l = 0; niters_glob_l = 0
        if( present(it_glob) )     it_glob_l     = it_glob
        if( present(niters_glob) ) niters_glob_l = niters_glob
        call fit_engine_iterate(params, build, fits, means, 1, niters, it_glob_l, niters_glob_l, .false., rounds)
        ! hoisted model handles back to the caller's arguments; everything else the fit holds is freed
        call move_alloc(fits(1)%model%basis_recs, model%basis_recs)
        call move_alloc(fits(1)%model%eigvals,    model%eigvals)
        model%ncomp = fits(1)%model%ncomp
        call fits(1)%kill
    end subroutine probe_subspace_iteration

    !> THE FIT ENGINE: one master loop advancing nfits resident fits -- 1 for the single-fit probe
    !! and the joint polish pass, 2 for the paired engine over the mod-4 halves (plan par.3.3: one
    !! it_eff advances BOTH fits). Per iteration every fit begins its iteration, the E-step runs
    !! ONCE over the fits' merged read list (fit_estep_pass; distributed: one qsys round per
    !! iteration, because the basis the E-step projects changes every iteration, then the parts'
    !! reduce; a worker writes its part and returns), then each fit's master tail runs: the
    !! cross-fit-FSC prep, the optional merge stash and the coupled M-step solve. Everything below
    !! the reduce is master-only, so the shared-memory result remains the reference. The loop
    !! exits when every fit has converged.
    !!
    !! EFFECTIVE (GLOBAL) ITERATION NUMBERING: a distributed worker is relaunched once per master
    !! EM iteration with niters=1, so its local counter is pinned at 1 in every round. Every
    !! iteration-keyed schedule (mixture warm-up, one-time checks) and every iteration log keys
    !! off it_eff/niters_eff, never it/niters (a schedule keyed off the local counter once left
    !! ~24% of every shard out of the M-step under nparts>1). The master stamps the true iteration
    !! through the stage request; shared memory reduces to it_eff == it.
    module subroutine fit_engine_iterate( params, build, fits, means, nfits, niters, it_glob, niters_glob, l_merge_stash, rounds )
        class(flex_pca_rounds), intent(inout) :: rounds
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        integer,             intent(in)    :: nfits
        type(flex_probe_fit), intent(inout) :: fits(nfits)
        type(flex_mean_ref), intent(in)    :: means(nfits)
        integer,             intent(in)    :: niters
        !> the master's global-iteration stamp on a distributed worker; 0 elsewhere
        integer,             intent(in)    :: it_glob, niters_glob
        !> paired final stage (par.7 "merge, don't refit"): snapshot each fit's raw M-step
        !! sufficient statistics + entry frame every iteration, pre-ridge pre-solve
        logical,             intent(in)    :: l_merge_stash
        !> cross-fit-FSC driver context: the paired master writes honest paired=1 records every
        !! iteration; the single-fit engine writes none; the ridge consumers act per their gates
        type(xfsc_ctx_t) :: xfctx
        character(len=:), allocatable :: xtag
        integer  :: it, it_eff, niters_eff, f, nthr
        logical  :: l_distr
        integer(timer_int_kind) :: t_it
        nthr = omp_get_max_threads()
        ! ---- per-fit stage config: the E-step formulation and every per-fit policy value ----
        do f = 1, nfits
            call fits(f)%estep_begin_stage(params, nthr)
        end do
        ! ---- CROSS-FIT-FSC setup: the ridge defaults OFF; only the paired master is a writer ----
        call xfsc_setup(xfctx, params, fits(1)%spec%cfg, fits(1)%spec%kfr_ann, nfits == 2, .not. rounds%is_worker())
        l_distr = rounds%distributed()
        ! ---- DOWNSCALED-PARTICLE CACHE ----
        ! Every EM iteration re-reads and re-preps the SAME particles, and one qsys round is
        ! launched per iteration, so the workers are fresh processes each time and nothing held in
        ! memory survives them. The on-disk cache does survive: it stores the iteration-independent
        ! prefix (noise normalisation, FFT, crop to box_crop), which at box=360/box_crop=64 is where
        ! the measured 48% of probe time goes. Masking keeps the full box.
        if( plane_cache_in_use(params, build) )then
            write(logfhandle,'(A)') '>>> FLEX_PCA PROBE: particles served from the downscaled cache'
            call flush(logfhandle)
        endif
        if( nfits == 2 )then
            if( l_distr )then
                write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA PAIRED DISTRIBUTED: E-step over ', &
                    &rounds%nparts(), ' parts, one v5 part per worker, per-fit reduce'
                call flush(logfhandle)
            endif
            write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA PAIRED: basis dim A=',fits(1)%model%ncomp, &
                &' B=',fits(2)%model%ncomp,'  CPU+POLAR E-step'
            call flush(logfhandle)
        endif
        niters_eff = niters
        if( niters_glob > 0 ) niters_eff = niters_glob
        do it = 1, niters
            it_eff = it
            if( it_glob > 0 ) it_eff = it_glob
            t_it = tic()
            do f = 1, nfits
                call fits(f)%iter_begin(params, build, means(f)%p, it_eff, niters_eff, nthr)
                fits(f)%model%z = 0.d0
            end do
            ! ---- ACCUMULATE: distributed when the master has parts, in-process otherwise ----
            if( l_distr )then
                call dispatch_round(params, fits, nfits, it_eff, niters_eff, rounds)
            else
                call fit_estep_pass(params, build, fits, means, nfits, it_eff, nthr)
                ! fold the thread accumulators (the distributed reduce above sums every worker's
                ! contribution into the fit accumulators instead)
                do f = 1, nfits
                    call fits(f)%iter_reduce(it_eff, nthr, rounds=rounds)
                end do
                if( rounds%is_worker() )then
                    ! a worker ships its accumulators and stops: no master tail runs here
                    call write_worker_parts(params, fits, nfits, nthr)
                    call xfsc_teardown(xfctx)
                    return
                endif
            endif
            ! ---- per-fit master tails ----
            do f = 1, nfits
                if( nfits == 2 )then
                    write(logfhandle,'(A,A,A,I0)') '>>> FLEX_PCA PAIRED FIT ', merge('A','B',f==1), &
                        &' master tail, it=', it_eff
                    call flush(logfhandle)
                    xtag = merge('  fit=A','  fit=B', f == 1)
                else
                    xtag = ''
                endif
                ! crossfsc per-iteration prep: the harvest flag for the writer payloads (H from the
                ! rho pair diagonals pre-ridge, internal FSC + Gamma) and the per-fit SSNR ridge
                ! from record t-1 (applied by iter_finish just before the coupled solves; fit X's
                ! ridge maps the record through fit X's side of the signed permutation with fit X's
                ! own H -- the spec par.2.3 invariant, keyed off fit%spec%id; on scaffolding
                ! records the arms degrade to arm 0 and log it)
                if( xfctx%l_any ) call xfsc_prep_iter(xfctx, params, fits(f), it_eff, xtag)
                ! par.7 merge stash: raw last-iteration statistics, BEFORE iter_finish ridges rho /
                ! solves the numerators in place / frees the accumulators. Every iteration
                ! overwrites -- any iteration can turn out to be the last.
                if( l_merge_stash ) call fits(f)%merge_stash
                call fits(f)%iter_finish(params, build, it_eff, nthr)
                if( nfits == 2 )then
                    write(logfhandle,'(A,A,A,I0,A,I0,A,ES12.4)') '>>> FLEX_PCA PAIRED FIT ', &
                        &merge('A','B',f==1), ' it=', it_eff, '  refined dim=', fits(f)%model%ncomp, &
                        &'  max var=', maxval(fits(f)%model%eigvals)
                    call flush(logfhandle)
                endif
            end do
            if( nfits == 2 )then
                ! the cross-fit comparison + artifact record + stopping series: one honest paired=1
                ! record per iteration from the two fits' realized bases and stashed payloads
                if( xfctx%l_writer ) call xfsc_paired_record(xfctx, params, fits, it_eff)
                write(logfhandle,'(A,I0,A,I0,A,F8.1)') '>>> FLEX_PCA PAIRED ITER ', it_eff, ' / ', &
                    &niters_eff, '  seconds=', toc(t_it)
            else
                write(logfhandle,'(A,I0,A,ES12.4,A,ES12.4,A,F8.1)') '>>> FLEX_PCA PROBE ITER ',it_eff, &
                    &' refined dim=',real(fits(1)%model%ncomp),' max var=',maxval(fits(1)%model%eigvals),' seconds=',toc(t_it)
            endif
            call flush(logfhandle)
            ! exit when EVERY fit has converged (a converged fit keeps updating until then)
            if( all(fits(1:nfits)%history%l_converged) )then
                if( nfits == 2 )then
                    write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA PAIRED converged after ', it_eff, &
                        &' iterations (both fits)'
                else
                    write(logfhandle,'(A,I0,A,F9.6)') '>>> FLEX_PCA PROBE converged after ',it_eff, &
                        &' iterations: rank-1 basis cosine vs previous >= ',fits(1)%spec%conv_thresh
                endif
                call flush(logfhandle)
                exit
            endif
        end do
        call xfsc_teardown(xfctx)
    end subroutine fit_engine_iterate

    !> One distributed E-step round on the master: stamp the probe-state file(s) a relaunched
    !! worker loads (the basis volumes are already on disk under each fit's prefix from the
    !! previous iter_finish / the data-free init), run ONE qsys round carrying the iteration, the
    !! budget and the fit count in job_descr, then fold every worker's part into the fits'
    !! accumulators. A single fit always stamps the default meta (the worker reads that one; the
    !! POLISH stage selects the polished basis namespace); the paired fits stamp their own.
    subroutine dispatch_round( params, fits, nfits, it_eff, niters_eff, rounds )
        class(flex_pca_rounds), intent(inout) :: rounds
        class(parameters),      intent(inout) :: params
        integer,                intent(in)    :: nfits, it_eff, niters_eff
        type(flex_probe_fit),   intent(inout) :: fits(nfits)
        integer :: f
        if( nfits == 2 )then
            do f = 1, 2
                call save_probe_state(fits(f)%model%ncomp, fits(f)%model%eigvals, fits(f)%model%sig2_eff, &
                    &fname=fits(f)%spec%meta_fname%to_char())
            end do
            call rounds%run_stage(params, flex_stage_request(stage=PCA_STAGE_PROBE, &
                &label='paired probe iteration', which_iter=it_eff, maxits=niters_eff, nfits=2))
            call paired_reduce_parts_v5(params, fits, rounds=rounds)
        else
            call save_probe_state(fits(1)%model%ncomp, fits(1)%model%eigvals, fits(1)%model%sig2_eff)
            ! the joint fit after the merge refines the polished basis, not fit A's
            if( fits(1)%spec%fprefix%to_char() == 'flex_pca_polished_pc' )then
                call rounds%run_stage(params, flex_stage_request(stage=PCA_STAGE_POLISH, &
                    &label='joint-fit iteration', which_iter=it_eff, maxits=niters_eff, nfits=1))
            else
                call rounds%run_stage(params, flex_stage_request(stage=PCA_STAGE_PROBE, &
                    &label='probe iteration', which_iter=it_eff, maxits=niters_eff, nfits=1))
            endif
            call reduce_single_fit_parts(params, fits(1), rounds%nparts())
        endif
    end subroutine dispatch_round

    !> Distributed master, single fit: fold every worker's v12 part into the fit's accumulators
    !! (the MCFA reduce buffers sized to the current rank first; the fit's own accumulators fold
    !! in place -- borrowed into the payload value, never copied).
    subroutine reduce_single_fit_parts( params, fit, nparts )
        class(parameters),    intent(in)    :: params
        type(flex_probe_fit), intent(inout) :: fit
        integer,              intent(in)    :: nparts
        type(flex_probe_part) :: part
        integer :: q
        allocate(part%cmat_e(fit%mstep%es(1),fit%mstep%es(2),fit%mstep%es(3),fit%model%ncomp), &
            &part%cmat_o(fit%mstep%es(1),fit%mstep%es(2),fit%mstep%es(3),fit%model%ncomp), source=(0.,0.))
        allocate(part%rho_ex(fit%mstep%es(1),fit%mstep%es(2),fit%mstep%es(3),fit%model%ncomp), &
            &part%rho_ox(fit%mstep%es(1),fit%mstep%es(2),fit%mstep%es(3),fit%model%ncomp), source=0.)
        if( fit%spec%l_mix_req )then
            if( allocated(fit%history%dm_sm) )then
                if( size(fit%history%dm_sm,1) /= fit%model%ncomp .or. size(fit%history%dm_sr) /= fit%spec%kmix ) &
                    &deallocate(fit%history%dm_sr, fit%history%dm_sm, fit%history%dm_smm, fit%history%dm_sai, fit%history%dm_z)
            endif
            if( .not. allocated(fit%history%dm_sr) )then
                allocate(fit%history%dm_sr(fit%spec%kmix), fit%history%dm_sm(fit%model%ncomp,fit%spec%kmix), &
                    &fit%history%dm_smm(fit%model%ncomp,fit%model%ncomp,fit%spec%kmix), fit%history%dm_sai(fit%model%ncomp,fit%model%ncomp))
                allocate(fit%history%dm_z(MIX_ZSUB_MAX*nparts, fit%model%ncomp))
            endif
            fit%history%dm_sr = 0.d0; fit%history%dm_sm = 0.d0; fit%history%dm_smm = 0.d0; fit%history%dm_sai = 0.d0; fit%history%dm_nz = 0
        endif
        call probe_part_borrow(part, fit, l_mix_buffers=fit%spec%l_mix_req)
        call reduce_probe_parts(params, nparts, part)
        call probe_part_restore(part, fit)
        if( fit%spec%l_mix_req )then
            write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA MIX distributed reduce: kmix=', &
                &fit%spec%kmix,'  pooled latent subsample rows=',fit%history%dm_nz
            call flush(logfhandle)
        endif
        do q = 1, fit%model%ncomp
            fit%mstep%Yeven(q)%cmat_exp = part%cmat_e(:,:,:,q); fit%mstep%Yeven(q)%rho_exp = part%rho_ex(:,:,:,q)
            fit%mstep%Yodd(q)%cmat_exp  = part%cmat_o(:,:,:,q); fit%mstep%Yodd(q)%rho_exp  = part%rho_ox(:,:,:,q)
        end do
        call part%kill
    end subroutine reduce_single_fit_parts

    !> A worker's part(s) for this iteration: one v12 part for a single fit, ONE v5 part carrying
    !! both fits' accumulator blocks (in fit order) for the paired engine.
    subroutine write_worker_parts( params, fits, nfits, nthr )
        class(parameters),    intent(in)    :: params
        integer,              intent(in)    :: nfits, nthr
        type(flex_probe_fit), intent(inout) :: fits(nfits)
        type(flex_probe_part) :: part
        type(string) :: pfname, tmp_fname
        integer :: f, q, funit
        pfname = flex_pca_part_fname('probe', params%part, params%numlen)
        if( nfits == 2 ) call open_probe_part_v5_write(pfname, 2, funit, tmp_fname)
        do f = 1, nfits
            allocate(part%cmat_e(fits(f)%mstep%es(1),fits(f)%mstep%es(2),fits(f)%mstep%es(3),fits(f)%model%ncomp), &
                &part%cmat_o(fits(f)%mstep%es(1),fits(f)%mstep%es(2),fits(f)%mstep%es(3),fits(f)%model%ncomp))
            allocate(part%rho_ex(fits(f)%mstep%es(1),fits(f)%mstep%es(2),fits(f)%mstep%es(3),fits(f)%model%ncomp), &
                &part%rho_ox(fits(f)%mstep%es(1),fits(f)%mstep%es(2),fits(f)%mstep%es(3),fits(f)%model%ncomp))
            do q = 1, fits(f)%model%ncomp
                part%cmat_e(:,:,:,q) = fits(f)%mstep%Yeven(q)%cmat_exp; part%rho_ex(:,:,:,q) = fits(f)%mstep%Yeven(q)%rho_exp
                part%cmat_o(:,:,:,q) = fits(f)%mstep%Yodd(q)%cmat_exp;  part%rho_ox(:,:,:,q) = fits(f)%mstep%Yodd(q)%rho_exp
            end do
            fits(f)%iter%nll_tot = sum(fits(f)%iter%nll_thr)
            if( fits(f)%spec%l_mix_req ) call probe_part_mixture_from_threads(part, fits(f), nthr, MIX_ZSUB_MAX)
            call probe_part_borrow(part, fits(f), l_mix_buffers=.false.)
            if( nfits == 2 )then
                call write_probe_part_v5_fit(funit, part)
            else
                call write_probe_part(pfname, part)
            endif
            call probe_part_restore(part, fits(f))
            call part%kill
        end do
        if( nfits == 2 )then
            call close_probe_part_v5_write(funit, tmp_fname, pfname)
            write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> FLEX_PCA PAIRED WORKER part=', params%part, &
                &'  wrote v5 part: valid A=', fits(1)%iter%nval, '  B=', fits(2)%iter%nval
            call flush(logfhandle)
        endif
        call pfname%kill
    end subroutine write_worker_parts

    !> A worker's MCFA payload: the thread-summed additive statistics and a bounded deterministic
    !! stride subsample of this part's latents (the master seeds mcfa_init from the pool).
    subroutine probe_part_mixture_from_threads( part, fit, nthr, nzsub_max )
        type(flex_probe_part), intent(inout) :: part
        type(flex_probe_fit),        intent(in)    :: fit
        integer,               intent(in)    :: nthr, nzsub_max
        integer :: tt2, kk4, nzs, izs, istep
        allocate(part%mix_sr(fit%spec%kmix), part%mix_sm(fit%model%ncomp,fit%spec%kmix), &
            &part%mix_smm(fit%model%ncomp,fit%model%ncomp,fit%spec%kmix), part%mix_sainv(fit%model%ncomp,fit%model%ncomp))
        part%mix_sr = 0.d0; part%mix_sm = 0.d0; part%mix_smm = 0.d0; part%mix_sainv = 0.d0
        do tt2 = 1, nthr
            part%mix_sr    = part%mix_sr    + fit%history%mxa_sr(:,tt2)
            part%mix_sm    = part%mix_sm    + fit%history%mxa_sm(:,:,tt2)
            part%mix_sainv = part%mix_sainv + fit%history%mxa_sainv(:,:,tt2)
            do kk4 = 1, fit%spec%kmix
                part%mix_smm(:,:,kk4) = part%mix_smm(:,:,kk4) + fit%history%mxa_smm(:,:,kk4,tt2)
            end do
        end do
        ! deterministic stride subsample of this part's latents (bounded)
        nzs   = min(nzsub_max, fit%spec%npp)
        istep = max(1, fit%spec%npp / max(1,nzs))
        nzs   = min(nzs, (fit%spec%npp + istep - 1)/istep)
        allocate(part%z_sub(nzs,fit%model%ncomp))
        do izs = 1, nzs
            part%z_sub(izs,:) = fit%model%z(min(fit%spec%npp, 1 + (izs-1)*istep), :)
        end do
        part%nz = nzs
    end subroutine probe_part_mixture_from_threads

end submodule simple_flex_probe_fit_engine
