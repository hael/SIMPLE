!@descr: flex_pca: the fit engine: one master loop for single, paired and worker fits, the distributed round, the probe part codec and its reduces
submodule (simple_flex_probe_fit) simple_flex_probe_fit_engine
!$ use omp_lib, only: omp_get_max_threads
use simple_core_module_api, only: del_file, dtiny, fclose, file_exists, fileiochk, fopen, logfhandle, &
    &simple_rename, tic, timer_int_kind, toc
use simple_flex_pca_stages,      only: flex_stage_request, PCA_STAGE_PROBE, PCA_STAGE_POLISH
use simple_flex_pca_artifacts,   only: FLEX_PCA_PART_MAGIC
use simple_flex_pca_posterior,   only: mcfa_init
use simple_flex_pca_basis,       only: save_probe_state, COV_PROBE_META
use simple_flex_pca_util,        only: cov_stage_subsample
use simple_flex_pca_fit_types,   only: probe_part_borrow, probe_part_restore
implicit none
#include "simple_local_flags.inc"

contains

    !> Master-side reduce: fold every worker's paired part into both fits' accumulators,
    !! streaming -- one part resident at a time, per-fit blocks in fit order, file deleted after.
    module subroutine paired_reduce_parts( params, fits, rounds )
        class(flex_pca_rounds), intent(inout) :: rounds
        class(parameters),      intent(in)    :: params
        type(flex_probe_fit),   intent(inout) :: fits(2)
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
                    allocate(fits(f)%history%dm_z(MIX_ZSUB_MAX*rounds%nparts(), fits(f)%model%ncomp))
                endif
                fits(f)%history%dm_sr = 0.d0; fits(f)%history%dm_sm = 0.d0; fits(f)%history%dm_smm = 0.d0
                fits(f)%history%dm_sai = 0.d0; fits(f)%history%dm_nz = 0
            endif
            ! the fit's own accumulators fold in place: borrowed, never copied
            call probe_part_borrow(acc(f), fits(f), l_mix_buffers=fits(f)%spec%l_mix_req)
        end do
        do ipart = 1, rounds%nparts()
            fname = rounds%part_fname('probe', ipart, params%numlen)
            call open_probe_part_read(fname, 2, funit)
            do f = 1, 2
                call fold_probe_part_fit(funit, acc(f))
            end do
            call close_probe_part_read(funit, fname)
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
        write(logfhandle,'(A,I0,A,I0,A,I0,A,F8.1)') '>>> FLEX_PCA reduced paired parts=', &
            &rounds%nparts(), '  valid A=', fits(1)%iter%nval, '  valid B=', fits(2)%iter%nval, &
            &'  seconds=', toc(t_red)
        call flush(logfhandle)
    end subroutine paired_reduce_parts

    !> Write one fit through the same container and payload codec used by paired parts.
    module subroutine write_probe_part( fname, part )
        class(string),         intent(in) :: fname
        type(flex_probe_part), intent(in) :: part
        type(string) :: tmp_fname
        integer :: funit
        call open_probe_part_write(fname, 1, funit, tmp_fname)
        call write_probe_part_fit(funit, part)
        call close_probe_part_write(funit, tmp_fname, fname)
    end subroutine write_probe_part

    !> Fold single-fit parts through the shared streaming codec.
    module subroutine reduce_probe_parts( params, rounds, part )
        class(parameters),     intent(in)    :: params
        class(flex_pca_rounds), intent(inout) :: rounds
        type(flex_probe_part), intent(inout) :: part
        type(string) :: fname
        integer :: ipart, funit
        integer(timer_int_kind) :: t_red
        t_red = tic()
        do ipart = 1, rounds%nparts()
            fname = rounds%part_fname('probe', ipart, params%numlen)
            call open_probe_part_read(fname, 1, funit)
            call fold_probe_part_fit(funit, part)
            call close_probe_part_read(funit, fname)
            call fname%kill
        end do
        write(logfhandle,'(A,I0,A,I0,A,F8.1)') '>>> FLEX_PCA reduced probe parts=',rounds%nparts(), &
            &'  valid particles=',part%nval,'  seconds=',toc(t_red)
        call flush(logfhandle)
    end subroutine reduce_probe_parts

    !> Open one part for writing and emit the file header. The caller then writes one
    !! block per fit with write_probe_part_fit, in fit order, and closes with
    !! close_probe_part_write (which renames the .tmp so the master only ever sees
    !! complete files -- same contract as write_probe_part).
    module subroutine open_probe_part_write( fname, nfits, funit, tmp_fname )
        class(string), intent(in)  :: fname
        integer,       intent(in)  :: nfits
        integer,       intent(out) :: funit
        type(string),  intent(out) :: tmp_fname
        integer :: io_stat
        tmp_fname = fname//'.tmp'
        call fopen(funit, file=tmp_fname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk('open_probe_part_write; open', io_stat)
        write(funit, iostat=io_stat) FLEX_PCA_PART_MAGIC, PROBE_PART_VERSION, nfits
        call fileiochk('open_probe_part_write; header', io_stat)
    end subroutine open_probe_part_write

    !> One fit's block: sub-header, band-packed accumulators and packed PCG arrays.
    module subroutine write_probe_part_fit( funit, part )
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
        call fileiochk('write_probe_part_fit; sub-header', io_stat)
        ! the coupled rho rows must live on the SAME lattice as the basis accumulators, otherwise
        ! one index list cannot describe both and the gathered payload would be mis-scattered
        if( size(part%rho_e,2) /= size(part%cmat_e,1) .or. size(part%rho_e,3) /= size(part%cmat_e,2) .or. &
            &size(part%rho_e,4) /= size(part%cmat_e,3) ) &
            &THROW_HARD('probe part coupled and basis lattices differ')
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
        call fileiochk('write_probe_part_fit; index count', io_stat)
        write(funit, iostat=io_stat) pii, pjj, pkk
        call fileiochk('write_probe_part_fit; index list', io_stat)
        allocate(gce(pnnz), gco(pnnz), gre(pnnz), gro(pnnz))
        do q = 1, part%ncomp
            do t = 1, pnnz
                gce(t) = part%cmat_e(pii(t),pjj(t),pkk(t),q)
                gre(t) = part%rho_ex(pii(t),pjj(t),pkk(t),q)
                gco(t) = part%cmat_o(pii(t),pjj(t),pkk(t),q)
                gro(t) = part%rho_ox(pii(t),pjj(t),pkk(t),q)
            end do
            write(funit, iostat=io_stat) gce, gre, gco, gro
            call fileiochk('write_probe_part_fit; basis payload', io_stat)
        end do
        deallocate(gce, gco, gre, gro)
        allocate(gpe(size(part%rho_e,1),pnnz), gpo(size(part%rho_o,1),pnnz))
        do t = 1, pnnz
            gpe(:,t) = part%rho_e(:,pii(t),pjj(t),pkk(t))
            gpo(:,t) = part%rho_o(:,pii(t),pjj(t),pkk(t))
        end do
        write(funit, iostat=io_stat) gpe, gpo, part%gam_sum, part%nll_sum
        call fileiochk('write_probe_part_fit; coupled payload', io_stat)
        deallocate(gpe, gpo, pii, pjj, pkk)
        if( kmix_w > 0 )then
            write(funit, iostat=io_stat) part%mix_sr, part%mix_sm, part%mix_smm, part%mix_sainv
            call fileiochk('write_probe_part_fit; mixture payload', io_stat)
            if( allocated(part%z_sub) )then
                write(funit, iostat=io_stat) size(part%z_sub,1)
                write(funit, iostat=io_stat) part%z_sub
            else
                write(funit, iostat=io_stat) 0
            endif
            call fileiochk('write_probe_part_fit; z subsample', io_stat)
        endif
        call write_probe_part_kernels(funit, part%kpk_e, part%kpk_o, part%rpk_e, part%rpk_o, 'write_probe_part_fit')
    end subroutine write_probe_part_fit

    ! ---- index-list packing --------------------------------------------------------------------
    ! Part payloads are nonzero only inside the covariance band, so only the populated lattice points
    ! ship: one index list per fit from the union of the shipped arrays' nonzeros (omitted voxels are
    ! exact zeros, so the fold is bit-identical). The PCG kernels and rhs ship packed on the band list.

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
    !! length travels in the shape and is checked.
    subroutine write_probe_part_kernels( funit, kpk_e, kpk_o, rpk_e, rpk_o, who )
        integer,              intent(in) :: funit
        real,    allocatable, intent(in) :: kpk_e(:,:), kpk_o(:,:)
        complex, allocatable, intent(in) :: rpk_e(:,:), rpk_o(:,:)
        character(len=*),     intent(in) :: who
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
        integer,              intent(in)    :: funit
        real,    allocatable, intent(inout) :: kpk_e(:,:), kpk_o(:,:)
        complex, allocatable, intent(inout) :: rpk_e(:,:), rpk_o(:,:)
        character(len=*),     intent(in)    :: who
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

    module subroutine close_probe_part_write( funit, tmp_fname, fname )
        integer,       intent(in)    :: funit
        type(string),  intent(inout) :: tmp_fname
        class(string), intent(in)    :: fname
        call fclose(funit)
        call simple_rename(tmp_fname, fname)
        call tmp_fname%kill
    end subroutine close_probe_part_write

    !> Open one part for reduction and validate its fit count before streaming the fit blocks.
    module subroutine open_probe_part_read( fname, nfits, funit )
        class(string), intent(in)  :: fname
        integer,       intent(in)  :: nfits
        integer,       intent(out) :: funit
        integer :: io_stat, magic, ver, nfits_in
        if( .not. file_exists(fname) ) THROW_HARD('missing probe part: '//fname%to_char())
        call fopen(funit, file=fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
        call fileiochk('open_probe_part_read; open '//fname%to_char(), io_stat)
        read(funit, iostat=io_stat) magic, ver, nfits_in
        call fileiochk('open_probe_part_read; header', io_stat)
        if( magic /= FLEX_PCA_PART_MAGIC ) THROW_HARD('bad probe part magic')
        if( ver /= PROBE_PART_VERSION ) THROW_HARD('bad probe part version')
        if( nfits_in /= nfits ) THROW_HARD('probe part nfits mismatch')
    end subroutine open_probe_part_read

    !> Fold one fit's block from an open part into that fit's accumulators (borrowed inside
    !! `part`); ncomp is checked against that fit's own expectation.
    module subroutine fold_probe_part_fit( funit, part )
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
        call fileiochk('fold_probe_part_fit; sub-header', io_stat)
        if( subhdr(1) /= part%ncomp ) THROW_HARD('probe part per-fit ncomp mismatch')
        if( subhdr(5) /= size(part%rho_e,1) ) THROW_HARD('probe part per-fit npairs mismatch')
        ! the lattice dims MUST be validated, not just ncomp/npairs: a part written on a different
        ! expanded lattice yields an in-range index list whose voxels land at the wrong addresses,
        ! i.e. silently misplaced mass instead of a loud failure
        if( subhdr(2) /= size(part%cmat_e,1) .or. subhdr(3) /= size(part%cmat_e,2) .or. &
            &subhdr(4) /= size(part%cmat_e,3) ) &
            &THROW_HARD('probe part per-fit lattice mismatch')
        ! index list of the populated lattice points (see "index-list packing"): the part carries
        ! only these voxels; every omitted one was an exact zero on the write side, so the fold is
        ! bit-identical to adding the full array
        read(funit, iostat=io_stat) pnnz
        call fileiochk('fold_probe_part_fit; index count', io_stat)
        if( pnnz < 1 .or. pnnz > product(subhdr(2:4)) ) &
            &THROW_HARD('probe part index count is outside the lattice extent')
        allocate(pii(pnnz), pjj(pnnz), pkk(pnnz))
        read(funit, iostat=io_stat) pii, pjj, pkk
        call fileiochk('fold_probe_part_fit; index list', io_stat)
        if( minval(pii) < 1 .or. maxval(pii) > size(part%cmat_e,1) .or. &
            &minval(pjj) < 1 .or. maxval(pjj) > size(part%cmat_e,2) .or. &
            &minval(pkk) < 1 .or. maxval(pkk) > size(part%cmat_e,3) ) &
            &THROW_HARD('probe part index list outside the basis lattice')
        allocate(gce(pnnz), gco(pnnz), gre(pnnz), gro(pnnz))
        allocate(gbuf(size(part%gam_sum)))
        do q = 1, part%ncomp
            read(funit, iostat=io_stat) gce, gre, gco, gro
            call fileiochk('fold_probe_part_fit; basis payload', io_stat)
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
        call fileiochk('fold_probe_part_fit; coupled payload', io_stat)
        !$omp parallel do default(shared) private(t) schedule(static)
        do t = 1, pnnz
            part%rho_e(:,pii(t),pjj(t),pkk(t)) = part%rho_e(:,pii(t),pjj(t),pkk(t)) + gpe(:,t)
            part%rho_o(:,pii(t),pjj(t),pkk(t)) = part%rho_o(:,pii(t),pjj(t),pkk(t)) + gpo(:,t)
        end do
        !$omp end parallel do
        deallocate(gpe, gpo)
        read(funit, iostat=io_stat) gbuf
        call fileiochk('fold_probe_part_fit; gamma', io_stat)
        part%gam_sum = part%gam_sum + gbuf
        read(funit, iostat=io_stat) nllbuf
        call fileiochk('fold_probe_part_fit; loglik', io_stat)
        part%nll_sum = part%nll_sum + nllbuf
        part%nval    = part%nval + subhdr(6)
        if( allocated(part%mix_sr) )then
            if( subhdr(7) /= size(part%mix_sr) ) THROW_HARD('probe part per-fit kmix mismatch')
            if( size(part%mix_sm,1) /= part%ncomp .or. size(part%mix_smm,1) /= part%ncomp .or. size(part%mix_sainv,1) /= part%ncomp ) &
                &THROW_HARD('probe part mixture accumulators have stale rank')
            allocate(msr(size(part%mix_sr)), msm(size(part%mix_sm,1),size(part%mix_sm,2)), &
                &msmm(size(part%mix_smm,1),size(part%mix_smm,2),size(part%mix_smm,3)), &
                &msai(size(part%mix_sainv,1),size(part%mix_sainv,2)))
            read(funit, iostat=io_stat) msr, msm, msmm, msai
            call fileiochk('fold_probe_part_fit; mixture', io_stat)
            part%mix_sr    = part%mix_sr    + msr
            part%mix_sm    = part%mix_sm    + msm
            part%mix_smm   = part%mix_smm   + msmm
            part%mix_sainv = part%mix_sainv + msai
            deallocate(msr, msm, msmm, msai)
            read(funit, iostat=io_stat) nz_part
            call fileiochk('fold_probe_part_fit; z subsample count', io_stat)
            if( nz_part < 0 .or. nz_part > 1000000 ) THROW_HARD('corrupt probe part latent count')
            if( nz_part > 0 )then
                allocate(zbuf(nz_part, size(part%mix_sm,1)))
                read(funit, iostat=io_stat) zbuf
                call fileiochk('fold_probe_part_fit; z subsample', io_stat)
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
            if( subhdr(7) /= 0 ) THROW_HARD('probe part carries an unexpected mixture block')
        endif
        call fold_probe_part_kernels(funit, part%kpk_e, part%kpk_o, part%rpk_e, part%rpk_o, 'fold_probe_part_fit')
        deallocate(gbuf, pii, pjj, pkk)
    end subroutine fold_probe_part_fit

    module subroutine close_probe_part_read( funit, fname )
        integer,       intent(in) :: funit
        class(string), intent(in) :: fname
        call fclose(funit)
        call del_file(fname)
    end subroutine close_probe_part_read

    !> Probe-based subspace iteration: alternate a Wiener E-step (per-particle latents in the current
    !! basis) with a weighted-backprojection M-step (Y_q += sum_i z_iq * backproject(r_i)), then
    !! orthonormalize the refined probe volumes into the next basis. The single-fit entry: the
    !! caller's mean and basis handles are hoisted into one fit, iterated by the shared engine
    !! (fit_engine_iterate), and handed back.
    module subroutine probe_subspace_iteration( params, build, plane_store, pcg_env, model, sel, niters, it_glob, &
        &niters_glob, fprefix, meta_fname, rounds )
        type(flex_fit_model),       intent(inout), target :: model   !< mean in; basis, prior variances, rank refined in place
        type(flex_selection),       intent(in)            :: sel
        class(flex_pca_rounds),     intent(inout)         :: rounds
        class(parameters),          intent(inout)         :: params
        type(builder),              intent(inout)         :: build
        class(flex_plane_store),    intent(inout)         :: plane_store
        class(flex_pcg_environment), intent(in)           :: pcg_env
        integer,                    intent(in)            :: niters
        !> master's global-iteration stamp, passed by a distributed probe worker whose own loop
        !! runs once per relaunch; absent (or 0) in shared-memory runs and on the master itself
        integer,          optional, intent(in)            :: it_glob, niters_glob
        !> file namespace override (merge-polish pass: never clobber fit A's delivery)
        character(len=*), optional, intent(in)            :: fprefix, meta_fname
        type(flex_probe_fit), target :: fits(1)
        type(flex_mean_ref) :: means(1)
        character(len=:), allocatable :: pfx_l, meta_l
        integer :: it_glob_l, niters_glob_l
        pfx_l  = 'flex_pca_pc'
        meta_l = COV_PROBE_META
        if( present(fprefix) )    pfx_l  = trim(fprefix)
        if( present(meta_fname) ) meta_l = trim(meta_fname)
        call fits(1)%new(0, pfx_l, meta_l, sel)
        allocate(fits(1)%mstep%env)
        call fits(1)%mstep%env%copy_from(pcg_env)
        call move_alloc(model%basis_recs, fits(1)%model%basis_recs)
        call move_alloc(model%eigvals,    fits(1)%model%eigvals)
        fits(1)%model%ncomp    = model%ncomp
        fits(1)%model%sig2_eff = model%sig2_eff
        fits(1)%model%sig2     = max(model%sig2_eff, DTINY)
        means(1)%p => model%mean_rec
        ! The probe pass uses every selected particle.
        call cov_stage_subsample(build%spproj_field, fits(1)%spec%sel%pinds, fits(1)%spec%sel%nptcls, 1, 0, 'PROBE', &
            &fits(1)%spec%ppinds, fits(1)%spec%npp)
        allocate(fits(1)%model%z(fits(1)%spec%npp,fits(1)%model%ncomp))
        it_glob_l = 0; niters_glob_l = 0
        if( present(it_glob) )     it_glob_l     = it_glob
        if( present(niters_glob) ) niters_glob_l = niters_glob
        call fit_engine_iterate(params, build, plane_store, fits, means, 1, niters, it_glob_l, niters_glob_l, .false., rounds)
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
    module subroutine fit_engine_iterate( params, build, plane_store, fits, means, nfits, niters, &
        &it_glob, niters_glob, l_merge_stash, rounds )
        class(flex_pca_rounds),  intent(inout) :: rounds
        class(parameters),       intent(inout) :: params
        type(builder),           intent(inout) :: build
        class(flex_plane_store), intent(inout) :: plane_store
        integer,                 intent(in)    :: nfits
        type(flex_probe_fit),    intent(inout) :: fits(nfits)
        type(flex_mean_ref),     intent(in)    :: means(nfits)
        integer,                 intent(in)    :: niters
        !> the master's global-iteration stamp on a distributed worker; 0 elsewhere
        integer,                intent(in)    :: it_glob, niters_glob
        !> paired final stage (par.7 "merge, don't refit"): snapshot each fit's raw M-step
        !! sufficient statistics + entry frame every iteration, pre-ridge pre-solve
        logical,                intent(in)    :: l_merge_stash
        !> cross-fit-FSC driver context: the paired master writes paired=1 records every
        !! iteration (the single-fit engine writes none); the ridge is their only consumer
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
        call xfsc_setup(xfctx, params, fits(1)%spec%kfr_ann, nfits == 2, .not. rounds%is_worker())
        l_distr = rounds%distributed()
        ! ---- DOWNSCALED-PARTICLE CACHE ----
        ! Every EM iteration re-reads and re-preps the SAME particles, and one qsys round is
        ! launched per iteration, so the workers are fresh processes each time and nothing held in
        ! memory survives them. The on-disk cache does survive: it stores the iteration-independent
        ! prefix (noise normalisation, FFT, crop to box_crop), which at box=360/box_crop=64 is where
        ! the measured 48% of probe time goes. Masking keeps the full box.
        if( plane_store%cache_in_use() )then
            write(logfhandle,'(A)') '>>> FLEX_PCA PROBE: particles served from the downscaled cache'
            call flush(logfhandle)
        endif
        if( nfits == 2 )then
            if( l_distr )then
                write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA PAIRED DISTRIBUTED: E-step over ', &
                    &rounds%nparts(), ' parts, one shared-codec part per worker, per-fit reduce'
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
                call fit_estep_pass(params, build, plane_store, fits, means, nfits, it_eff, nthr)
                ! fold the thread accumulators (the distributed reduce above sums every worker's
                ! contribution into the fit accumulators instead)
                do f = 1, nfits
                    call fits(f)%iter_reduce(it_eff, nthr, rounds=rounds)
                end do
                if( rounds%is_worker() )then
                    ! a worker ships its accumulators and stops: no master tail runs here
                    call write_worker_parts(params, fits, nfits, nthr, rounds)
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
            call paired_reduce_parts(params, fits, rounds=rounds)
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
            call reduce_single_fit_parts(params, fits(1), rounds)
        endif
    end subroutine dispatch_round

    !> Distributed master, single fit: fold every worker's part into the fit's accumulators
    !! (the MCFA reduce buffers sized to the current rank first; the fit's own accumulators fold
    !! in place -- borrowed into the payload value, never copied).
    subroutine reduce_single_fit_parts( params, fit, rounds )
        class(parameters),    intent(in)    :: params
        type(flex_probe_fit), intent(inout) :: fit
        class(flex_pca_rounds), intent(inout) :: rounds
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
                allocate(fit%history%dm_z(MIX_ZSUB_MAX*rounds%nparts(), fit%model%ncomp))
            endif
            fit%history%dm_sr = 0.d0; fit%history%dm_sm = 0.d0; fit%history%dm_smm = 0.d0; fit%history%dm_sai = 0.d0; fit%history%dm_nz = 0
        endif
        call probe_part_borrow(part, fit, l_mix_buffers=fit%spec%l_mix_req)
        call reduce_probe_parts(params, rounds, part)
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

    !> A worker's part for this iteration: one container holding one fit, or both fits'
    !! accumulator blocks in fit order for the paired engine.
    subroutine write_worker_parts( params, fits, nfits, nthr, rounds )
        class(parameters),    intent(in)    :: params
        integer,              intent(in)    :: nfits, nthr
        type(flex_probe_fit), intent(inout) :: fits(nfits)
        class(flex_pca_rounds), intent(inout) :: rounds
        type(flex_probe_part) :: part
        type(string) :: pfname, tmp_fname
        integer :: f, q, funit
        pfname = rounds%part_fname('probe', params%part, params%numlen)
        call open_probe_part_write(pfname, nfits, funit, tmp_fname)
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
            call write_probe_part_fit(funit, part)
            call probe_part_restore(part, fits(f))
            call part%kill
        end do
        call close_probe_part_write(funit, tmp_fname, pfname)
        if( nfits == 2 )then
            write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> FLEX_PCA PAIRED WORKER part=', params%part, &
                &'  wrote paired part: valid A=', fits(1)%iter%nval, '  B=', fits(2)%iter%nval
            call flush(logfhandle)
        endif
        call pfname%kill
    end subroutine write_worker_parts

    !> A worker's MCFA payload: the thread-summed additive statistics and a bounded deterministic
    !! stride subsample of this part's latents (the master seeds mcfa_init from the pool).
    subroutine probe_part_mixture_from_threads( part, fit, nthr, nzsub_max )
        type(flex_probe_part), intent(inout) :: part
        type(flex_probe_fit),  intent(in)    :: fit
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
