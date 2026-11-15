!@descr: unit test routines for the rec3D backend selector (simple_rec3D_strategy, simple_rec3D_service)
! Pins the parameter defaults the selector relies on (rec_backend=gridding, maxits_pcg=2),
! backend-name resolution and wiring, and the dynamic type the factory returns for every
! branch: shared-memory gridding, distributed gridding (nparts without part), shared-memory
! PCG, distributed PCG (nparts > 1) and the worker part of either backend; the PCG worker's entries of the
! weighted insertion (topk > 1).
submodule (simple_rec3D_strategy) simple_rec3D_strategy_tester
use simple_rec3D_service,  only: rec3D_backend_is_wired, REC3D_BACKEND_INVALID
use simple_reconstructor_pcg, only: reconstructor_pcg, PCG_LAMBDA
use simple_rec3D_pcg_strategy, only: pcg_extra_entries, accumulate_pcg_extras, pcg_inpl_particle_entries, &
    &accumulate_pcg_inpl_batch
use simple_eul_prob_tab_utils, only: prob_insert_poses, inpl_entry_pack, inpl_entry_rot, inpl_entry_count
use simple_matcher_3Drec,      only: insert_inpl_entries, count_inserted_poses
use simple_reconstructor,      only: reconstructor
use simple_image,              only: image
use simple_sp_project,         only: sp_project
use simple_fplane_grid,        only: fplane_grid
use simple_memoize_ft_maps,    only: forget_ft_maps, memoize_ft_maps
use simple_math_ft,            only: fplane_shift
!$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
use simple_type_defs,          only: CTFFLAG_NO
use simple_rnd,                only: seed_rnd_fixed
use simple_test_utils
implicit none

integer, parameter :: KIND_INMEM = 1, KIND_DISTR = 2, KIND_PCG_INMEM = 3, KIND_OTHER = 0

contains

    module subroutine run_all_rec3D_strategy_tests()
        write(*,'(A)') '**** running all rec3D_strategy tests ****'
        call test_defaults()
        call test_backend_resolution()
        call test_factory_branches()
        call test_pcg_extra_entries()
        call test_inpl_gridding_entries()
        call test_pcg_inpl_chunks()
    end subroutine run_all_rec3D_strategy_tests

    !> The PCG worker's weighted insertion: the extra poses as entries after the committed ones, each from its
    !! particle's observation. Extras at the committed pose with weights summing to one reproduce the plain
    !! accumulation, in one batch and in two; an extra moved off the committed pose changes it
    subroutine test_pcg_extra_entries()
        integer, parameter :: NP = 4, BOX = 32
        real,    parameter :: SMPD = 2.
        type(reconstructor_pcg) :: op
        type(oris)              :: sel, entries
        type(ori)               :: o
        type(ctfparams)         :: ctfp
        type(prob_insert_poses) :: extra
        complex, allocatable    :: planes(:,:,:), y_x(:,:,:), b(:,:,:), b0(:,:,:), b2(:,:,:)
        real,    allocatable    :: d(:,:,:), d0(:,:,:), d2(:,:,:), wts(:)
        integer, allocatable    :: xptcl(:), xfirst(:)
        integer :: pinds(NP), lims2(2,2), i, h, k
        write(*,'(A)') 'test_pcg_extra_entries'
        call seed_rnd_fixed(20261103)
        pinds        = [11, 12, 13, 14]
        ctfp%ctfflag = CTFFLAG_NO
        ctfp%smpd    = SMPD
        call sel%new(NP, .true.)
        call o%new(.true.)
        do i = 1, NP
            call o%set_euler([30.*real(i), 40. + 10.*real(i), 15.*real(i)])
            call o%set_ctfvars(ctfp)
            call o%set_shift([0.5*real(i), -0.3*real(i)])
            call sel%set_ori(i, o)
        enddo
        call op%new(BOX, SMPD, PCG_LAMBDA)
        lims2 = op%get_lims2()
        allocate(planes(lims2(1,1):lims2(1,2), lims2(2,1):lims2(2,2), NP))
        do i = 1, NP
            do k = lims2(2,1), lims2(2,2)
                do h = lims2(1,1), lims2(1,2)
                    planes(h,k,i) = cmplx(ran3() - 0.5, ran3() - 0.5)
                enddo
            enddo
        enddo
        call op%kill
        ! the plain accumulation
        call accumulate_plain(b0, d0)
        ! particle 12: one extra at its committed pose (0.6 / 0.4); 14: two (0.5 / 0.3 / 0.2); 11 and 13 none
        extra%nptcls = NP
        extra%npose  = 3
        allocate(extra%loc(maxval(pinds)), extra%nextra(NP), extra%first(NP), source=0)
        allocate(extra%wc(NP), source=1.)
        allocate(extra%euler(3,3), extra%shift(2,3), extra%w(3), source=0.)
        extra%loc(pinds) = [1, 2, 3, 4]
        extra%first      = [1, 1, 2, 2]
        extra%nextra     = [0, 1, 0, 2]
        extra%wc         = [1., 0.6, 1., 0.5]
        extra%w          = [0.4, 0.3, 0.2]
        call sel%get_ori(2, o)
        extra%euler(:,1) = o%get_euler()
        extra%shift(:,1) = o%get_2Dshift()
        call sel%get_ori(4, o)
        extra%euler(:,2) = o%get_euler()
        extra%euler(:,3) = o%get_euler()
        extra%shift(:,2) = o%get_2Dshift()
        extra%shift(:,3) = o%get_2Dshift()
        call pcg_extra_entries(sel, pinds, [(1., i = 1, NP)], extra, 1., entries, wts, xptcl, xfirst)
        call assert_int(NP + 3, entries%get_noris(), 'one entry per committed and per extra pose')
        call assert_true(all(xptcl == [2, 4, 4]) .and. all(xfirst == [1, 1, 2, 2, 4]), 'the extras follow in particle order')
        call assert_true(all(abs(wts - [1., 0.6, 1., 0.5, 0.4, 0.3, 0.2]) < 1.e-6), 'the entries carry the poses'' weights')
        allocate(y_x(lims2(1,1):lims2(1,2), lims2(2,1):lims2(2,2), 3))
        call accumulate_entries([1, NP], b, d)
        call assert_true(maxval(abs(b - b0)) <= 1.e-4 * maxval(abs(b0)), 'extras at the committed pose reproduce the plain right-hand side')
        call assert_true(maxval(abs(d - d0)) <= 1.e-4 * maxval(abs(d0)), 'extras at the committed pose reproduce the plain density')
        call accumulate_entries([2, NP], b2, d2)
        call assert_true(maxval(abs(b2 - b)) <= 1.e-5 * maxval(abs(b)) .and. maxval(abs(d2 - d)) <= 1.e-5 * maxval(abs(d)),&
            &'two batches accumulate as one')
        ! an extra moved off the committed pose enters at its own pose
        extra%euler(:,3) = extra%euler(:,3) + [20., 15., 0.]
        call entries%kill
        call pcg_extra_entries(sel, pinds, [(1., i = 1, NP)], extra, 1., entries, wts, xptcl, xfirst)
        call accumulate_entries([1, NP], b2, d2)
        call assert_true(maxval(abs(b2 - b)) > 1.e-3 * maxval(abs(b)), 'an extra pose of its own changes the right-hand side')
        call o%kill
        call sel%kill
        call entries%kill
        call extra%kill

    contains

        subroutine accumulate_plain( bp, dp_ )
            complex, allocatable, intent(out) :: bp(:,:,:)
            real,    allocatable, intent(out) :: dp_(:,:,:)
            call op%new(BOX, SMPD, PCG_LAMBDA)
            call op%prep_particles(sel, use_ctf=.true.)
            call op%begin_accum
            call op%accumulate_batch(planes, NP, 1)
            call op%get_raw_accum(bp, dp_)
            call op%kill
        end subroutine accumulate_plain

        !> the committed and the extra entries in batches of the given size (bsz(1)) over the NP particles
        subroutine accumulate_entries( bsz, be, de )
            integer,              intent(in)  :: bsz(2)
            complex, allocatable, intent(out) :: be(:,:,:)
            real,    allocatable, intent(out) :: de(:,:,:)
            integer :: ib, ie
            call op%new(BOX, SMPD, PCG_LAMBDA)
            call op%prep_particles(entries, use_ctf=.true., weights=wts)
            call op%begin_accum
            do ib = 1, bsz(2), bsz(1)
                ie = min(bsz(2), ib + bsz(1) - 1)
                call op%accumulate_batch(planes(:,:,ib:ie), ie - ib + 1, ib)
                call accumulate_pcg_extras(op, planes(:,:,ib:ie), [ib, ie], NP, xptcl, xfirst, y_x)
            enddo
            call op%get_raw_accum(be, de)
            call op%kill
        end subroutine accumulate_entries

    end subroutine test_pcg_extra_entries

    !> The key-on gridding loop (insert_inpl_entries) on a committed pose and two extras equals every entry inserted on
    !! its own from the committed plane re-shifted for it, bit for bit; entries at one pose whose weights sum to 1
    !! equal one insertion at weight 1
    subroutine test_inpl_gridding_entries()
        integer, parameter :: BOXG = 16, KQ = 4, NR = 16
        real,    parameter :: SMPDG = 2.0, SH_C(2) = [0.7, -1.3], WBASE = 0.7
        type(parameters), target :: params
        type(sp_project)         :: project
        type(reconstructor)      :: ra, rb
        type(prob_insert_poses)  :: extra
        type(ori)                :: o, oe
        type(sym)                :: c1sym
        type(ctfparams)          :: ctfparms
        type(fplane_type)        :: fpl, fpl_x, fpl_scratch
        type(fplane_grid)        :: grid
        type(image)              :: obs, obs_pad
        complex, allocatable     :: ca(:,:,:), cb(:,:,:)
        real,    allocatable     :: img(:,:,:), da(:,:,:), db(:,:,:)
        real    :: rots(NR), r, euler(3), wpose, e3, cerr, derr
        integer :: i, j, k, iblk, irot, cnt
        write(*,'(A)') 'test_inpl_gridding_entries'
        params%box        = BOXG
        params%box_crop   = BOXG
        params%box_croppd = OSMPL_PAD_FAC * BOXG
        params%smpd_crop  = SMPDG
        params%nstates    = 1
        params%numlen     = 1
        params%oritype    = 'cls3D'
        do i = 1, NR
            rots(i) = real(i - 1) * 360. / real(NR)
        enddo
        call obs%new([BOXG,BOXG,1], SMPDG, wthreads=.false.)
        call obs_pad%new([OSMPL_PAD_FAC*BOXG,OSMPL_PAD_FAC*BOXG,1], SMPDG, wthreads=.false.)
        call grid%new(OSMPL_PAD_FAC*BOXG)
        call c1sym%new('c1')
        call memoize_ft_maps([OSMPL_PAD_FAC*BOXG,OSMPL_PAD_FAC*BOXG,1], SMPDG)
        ctfparms%smpd    = SMPDG
        ctfparms%ctfflag = CTFFLAG_NO
        call seed_rnd_fixed(20261116)
        allocate(img(BOXG,BOXG,1))
        do k = 1, BOXG
            do j = 1, BOXG
                call random_number(r)
                img(j,k,1) = r - 0.5
            enddo
        enddo
        call obs%set_rmat(img, .false.)
        call obs%pad(obs_pad, backgr=0., antialiasing=.false.)
        call obs_pad%fft()
        call obs_pad%gen_fplane4rec(SMPDG, ctfparms, SH_C, fpl, grid)
        call o%new(.false.)
        call o%set_euler([37., 64., 112.5])
        call o%set_shift(SH_C)
        ! one particle (index 7): the committed pose and two extras at other directions and shifts
        extra%nptcls = 1
        extra%npose  = 2
        allocate(extra%loc(7), extra%nextra(1), extra%first(1), source=0)
        allocate(extra%wc(1), extra%euler(3,2), extra%shift(2,2), extra%w(2), source=0.)
        extra%loc(7)     = 1
        extra%nextra(1)  = 2
        extra%first(1)   = 1
        extra%wc(1)      = 0.5
        extra%w          = [0.3, 0.2]
        extra%euler(:,1) = [52., 70., 201.25]
        extra%euler(:,2) = [41., 59., 333.]
        extra%shift(:,1) = [1.6, 0.4]
        extra%shift(:,2) = [-0.9, -2.2]
        call extra%reserve_inpl(KQ, 3, rots)
        call set_block(0, 3, [3, 9, 0], [3, 1, 0])
        call set_block(1, 5, [5, 12, 0], [2, 2, 0])
        call set_block(2, 7, [13, 0, 0], [4, 0, 0])
        ! the loop
        call ra%new_accumulator(params, project, expand=.true., wthreads=.false.)
        call insert_inpl_entries(c1sym, o, fpl, ra, extra, 1, WBASE, 1., fpl_scratch)
        ! every entry on its own, its plane re-shifted from the committed plane for the entry
        call rb%new_accumulator(params, project, expand=.true., wthreads=.false.)
        call oe%new(.false.)
        oe = o
        do k = 0, 2
            if( k == 0 )then
                wpose = extra%wc(1)
                euler = o%get_euler()
            else
                wpose = extra%w(k)
                euler = extra%euler(:,k)
            endif
            iblk = extra%block(1, k)
            do j = 1, extra%blk_nent(iblk)
                irot = inpl_entry_rot(extra%blk_ent(j,iblk))
                cnt  = inpl_entry_count(extra%blk_ent(j,iblk))
                e3   = euler(3)
                if( irot /= extra%blk_bin(iblk) ) e3 = 360. - rots(irot)
                call oe%set_euler([euler(1), euler(2), e3])
                if( k == 0 )then
                    call rb%insert_plane_oversamp(c1sym, oe, fpl, w=WBASE * wpose * real(cnt) / real(KQ))
                else
                    call fplane_shift(fpl, extra%shift(:,k) - SH_C, OSMPL_PAD_FAC, fpl_x)
                    call rb%insert_plane_oversamp(c1sym, oe, fpl_x, w=WBASE * wpose * real(cnt) / real(KQ))
                endif
            enddo
        enddo
        ca = ra%cmat_exp
        cb = rb%cmat_exp
        da = ra%rho_exp
        db = rb%rho_exp
        call assert_true(all(ca == cb) .and. all(da == db), &
            &'the key-on loop equals each entry inserted on its own from its re-shifted plane, bit for bit')
        call ra%kill
        call rb%kill
        ! entries at one orientation and shift: the committed pose collapsed at 0.3 and an extra at the committed pose
        ! collapsed at 0.7 give one insertion at weight 1
        extra%wc(1)      = 0.3
        extra%w          = [0.7, 0.]
        extra%euler(:,1) = o%get_euler()
        extra%shift(:,1) = SH_C
        call set_block(0, 3, [3, 0, 0], [4, 0, 0])
        call set_block(1, 3, [3, 0, 0], [4, 0, 0])
        call ra%new_accumulator(params, project, expand=.true., wthreads=.false.)
        call insert_inpl_entries(c1sym, o, fpl, ra, extra, 1, 1., 1., fpl_scratch)
        call rb%new_accumulator(params, project, expand=.true., wthreads=.false.)
        call rb%insert_plane_oversamp(c1sym, o, fpl, w=1.)
        ca = ra%cmat_exp
        cb = rb%cmat_exp
        da = ra%rho_exp
        db = rb%rho_exp
        cerr = maxval(abs(ca - cb)) / max(maxval(abs(cb)), TINY)
        derr = maxval(abs(da - db)) / max(maxval(abs(db)), TINY)
        call assert_true(cerr < 1.e-6 .and. derr < 1.e-6, 'entries at one pose whose weights sum to 1 equal one insertion at weight 1')
        call ra%kill
        call rb%kill
        call extra%kill
        call forget_ft_maps()
        call c1sym%kill()
        call o%kill()
        call oe%kill()
        call obs%kill()
        call obs_pad%kill()

    contains

        !> pose k's block: the solved bin and up to three entries (rotation, count; count 0 ends the block)
        subroutine set_block( kpose, bin, irots, cnts )
            integer, intent(in) :: kpose, bin, irots(3), cnts(3)
            integer :: jj, ib
            ib = extra%block(1, kpose)
            extra%blk_bin(ib)  = bin
            extra%blk_nent(ib) = count(cnts > 0)
            do jj = 1, extra%blk_nent(ib)
                extra%blk_ent(jj,ib) = inpl_entry_pack(irots(jj), cnts(jj))
            enddo
        end subroutine set_block

    end subroutine test_inpl_gridding_entries

    !> The key-on PCG chunks (accumulate_pcg_inpl_batch, chunks of 3 entries, one and two observation batches) give the
    !! raw B and D of one flat accumulation of the same entries byte for byte, at one thread and at four; the flat list
    !! holds every inserted entry once
    subroutine test_pcg_inpl_chunks()
        integer, parameter :: NP = 4, BOX = 32, KQ = 4, NR = 16, CAP = 3
        real,    parameter :: SMPD = 2.
        type(reconstructor_pcg) :: op
        type(oris)              :: sel, flat
        type(ori)               :: o
        type(ctfparams)         :: ctfp
        type(prob_insert_poses) :: extra
        complex, allocatable    :: planes(:,:,:), yflat(:,:,:), bf(:,:,:), bc(:,:,:)
        real,    allocatable    :: df(:,:,:), dc(:,:,:), inv_sig2(:)
        real    :: rots(NR), eul(3,3*KQ), sh(2,3*KQ), wt(3*KQ), wbase(NP), wflat(64), eflat(3,64), sflat(2,64)
        integer :: pinds(NP), lims2(2,2), i, h, k, np_i, nflat, pflat(64), nthr_saved, iteam
        logical :: l_begun
        write(*,'(A)') 'test_pcg_inpl_chunks'
        call seed_rnd_fixed(20261117)
        pinds        = [11, 12, 13, 14]
        wbase        = [1., 0.8, 1., 0.6]
        ctfp%ctfflag = CTFFLAG_NO
        ctfp%smpd    = SMPD
        do i = 1, NR
            rots(i) = real(i - 1) * 360. / real(NR)
        enddo
        call sel%new(NP, .true.)
        call o%new(.true.)
        do i = 1, NP
            call o%set_euler([30.*real(i), 40. + 10.*real(i), 15.*real(i)])
            call o%set_ctfvars(ctfp)
            call o%set_shift([0.5*real(i), -0.3*real(i)])
            call sel%set_ori(i, o)
        enddo
        call op%new(BOX, SMPD, PCG_LAMBDA)
        lims2 = op%get_lims2()
        call op%kill
        allocate(planes(lims2(1,1):lims2(1,2), lims2(2,1):lims2(2,2), NP))
        do i = 1, NP
            do k = lims2(2,1), lims2(2,2)
                do h = lims2(1,1), lims2(1,2)
                    planes(h,k,i) = cmplx(ran3() - 0.5, ran3() - 0.5)
                enddo
            enddo
        enddo
        ! 11: entries at its committed pose; 12: one extra; 13: no blocks (the weighted insertion's rule); 14: two extras
        extra%nptcls = NP
        extra%npose  = 3
        allocate(extra%loc(maxval(pinds)), extra%nextra(NP), extra%first(NP), source=0)
        allocate(extra%wc(NP), source=1.)
        allocate(extra%euler(3,3), extra%shift(2,3), extra%w(3), source=0.)
        extra%loc(pinds) = [1, 2, 3, 4]
        extra%first      = [1, 1, 2, 2]
        extra%nextra     = [0, 1, 0, 2]
        extra%wc         = [1., 0.6, 1., 0.5]
        extra%w          = [0.4, 0.3, 0.2]
        extra%euler(:,1) = [70., 20., 99.]
        extra%euler(:,2) = [10., 80., 250.]
        extra%euler(:,3) = [120., 33., 17.5]
        extra%shift(:,1) = [1.1, -0.2]
        extra%shift(:,2) = [-0.6, 0.9]
        extra%shift(:,3) = [0.3, 1.7]
        call extra%reserve_inpl(KQ, 3, rots)
        call set_block(1, 0, 3, [3, 9, 14, 0], [2, 1, 1, 0])
        call set_block(2, 0, 2, [2, 7, 0, 0], [3, 1, 0, 0])
        call set_block(2, 1, 5, [13, 0, 0, 0], [4, 0, 0, 0])
        call set_block(4, 0, 1, [1, 4, 8, 12], [1, 1, 1, 1])
        call set_block(4, 1, 6, [6, 10, 0, 0], [2, 2, 0, 0])
        call set_block(4, 2, 2, [2, 3, 11, 0], [1, 2, 1, 0])
        ! the flat list
        nflat = 0
        do i = 1, NP
            call pcg_inpl_particle_entries(sel, pinds, wbase, extra, 1., i, eul, sh, wt, np_i)
            if( i == 2 ) call assert_real(90., eul(3,np_i), 0., &
                &'PCG keeps a collapsed extra at bin 13, away from its solved bin 5')
            eflat(:,nflat+1:nflat+np_i) = eul(:,1:np_i)
            sflat(:,nflat+1:nflat+np_i) = sh(:,1:np_i)
            wflat(nflat+1:nflat+np_i)   = wt(1:np_i)
            pflat(nflat+1:nflat+np_i)   = i
            nflat = nflat + np_i
        enddo
        call assert_true(int(nflat,int64) == count_inserted_poses(pinds, extra), &
            &'the flat list holds one entry per inserted pose entry (count_inserted_poses)')
        call assert_true(nflat > CAP * 2, 'the batch holds more entries than two chunks')
        call flat%new(nflat, .true.)
        allocate(yflat(lims2(1,1):lims2(1,2), lims2(2,1):lims2(2,2), nflat))
        do i = 1, nflat
            call sel%get_ori(pflat(i), o)
            call o%set_euler(eflat(:,i))
            call o%set_shift(sflat(:,i))
            call flat%set_ori(i, o)
            yflat(:,:,i) = planes(:,:,pflat(i))
        enddo
        call op%new(BOX, SMPD, PCG_LAMBDA)
        call op%prep_particles(flat, use_ctf=.true., weights=wflat(1:nflat))
        call op%begin_accum
        call op%accumulate_batch(yflat, nflat, 1)
        call op%get_raw_accum(bf, df)
        call op%kill
        nthr_saved = 1
        !$ nthr_saved = omp_get_max_threads()
        do iteam = 1, 2
            !$ if( iteam == 2 ) call omp_set_num_threads(4)
            ! one observation batch of four particles, chunks of three entries
            call op%new(BOX, SMPD, PCG_LAMBDA)
            l_begun = .false.
            call accumulate_pcg_inpl_batch(op, sel, pinds, wbase, extra, 1., [1, NP], planes, CAP, inv_sig2, l_begun)
            call op%get_raw_accum(bc, dc)
            call op%kill
            call assert_true(all(bc == bf) .and. all(dc == df), &
                &'chunks smaller than the batch accumulate B and D as the flat list, byte for byte')
            ! two observation batches
            call op%new(BOX, SMPD, PCG_LAMBDA)
            l_begun = .false.
            call accumulate_pcg_inpl_batch(op, sel, pinds, wbase, extra, 1., [1, 2], planes(:,:,1:2), CAP, inv_sig2, l_begun)
            call accumulate_pcg_inpl_batch(op, sel, pinds, wbase, extra, 1., [3, NP], planes(:,:,3:NP), CAP, inv_sig2, l_begun)
            call op%get_raw_accum(bc, dc)
            call op%kill
            call assert_true(all(bc == bf) .and. all(dc == df), &
                &'two observation batches in chunks accumulate as the flat list, byte for byte')
        enddo
        !$ call omp_set_num_threads(nthr_saved)
        call o%kill
        call sel%kill
        call flat%kill
        call extra%kill

    contains

        !> particle i's block for pose kpose: the solved bin and up to four entries (count 0 ends the block)
        subroutine set_block( ip, kpose, bin, irots, cnts )
            integer, intent(in) :: ip, kpose, bin, irots(4), cnts(4)
            integer :: jj, ib
            ib = extra%block(ip, kpose)
            extra%blk_bin(ib)  = bin
            extra%blk_nent(ib) = count(cnts > 0)
            do jj = 1, extra%blk_nent(ib)
                extra%blk_ent(jj,ib) = inpl_entry_pack(irots(jj), cnts(jj))
            enddo
        end subroutine set_block

    end subroutine test_pcg_inpl_chunks

    !> the dynamic type of a strategy as a small integer, exact type (pcg_inmem extends inmem)
    integer function strategy_kind( strategy )
        class(rec3D_strategy), intent(in) :: strategy
        select type(strategy)
            type is(rec3D_inmem_strategy)
                strategy_kind = KIND_INMEM
            type is(rec3D_distr_strategy)
                strategy_kind = KIND_DISTR
            type is(rec3D_pcg_inmem_strategy)
                strategy_kind = KIND_PCG_INMEM
            class default
                strategy_kind = KIND_OTHER
        end select
    end function strategy_kind

    subroutine test_defaults()
        type(parameters) :: params
        write(*,'(A)') 'test_defaults'
        call assert_char('gridding', trim(params%rec_backend), 'rec_backend defaults to gridding')
        call assert_int(2, params%maxits_pcg, 'maxits_pcg defaults to 2')
        call assert_true(params%rtol <= 0., 'rtol defaults to <= 0 (exactly maxits_pcg iterations)')
    end subroutine test_defaults

    subroutine test_backend_resolution()
        write(*,'(A)') 'test_backend_resolution'
        call assert_int(REC3D_BACKEND_GRIDDING, rec3D_backend_id('gridding'), 'gridding resolves')
        call assert_int(REC3D_BACKEND_PCG,      rec3D_backend_id('pcg'),      'pcg resolves')
        call assert_int(REC3D_BACKEND_GRIDDING, rec3D_backend_id('gridding   '), 'trailing blanks are ignored')
        call assert_int(REC3D_BACKEND_INVALID,  rec3D_backend_id('invalid'),  'unknown name is invalid')
        call assert_int(REC3D_BACKEND_INVALID,  rec3D_backend_id('PCG'),      'names are case-sensitive')
        call assert_int(REC3D_BACKEND_INVALID,  rec3D_backend_id(''),         'empty name is invalid')
        call assert_true(rec3D_backend_is_wired(REC3D_BACKEND_GRIDDING), 'gridding is wired')
        call assert_true(rec3D_backend_is_wired(REC3D_BACKEND_PCG),      'pcg is wired')
        call assert_false(rec3D_backend_is_wired(REC3D_BACKEND_INVALID), 'invalid is not wired')
    end subroutine test_backend_resolution

    subroutine test_factory_branches()
        type(cmdline) :: cline
        class(rec3D_strategy), allocatable :: strategy
        write(*,'(A)') 'test_factory_branches'
        ! bare command line: shared-memory gridding
        strategy = create_rec3D_strategy(cline)
        call assert_int(KIND_INMEM, strategy_kind(strategy), 'default is the gridding in-memory strategy')
        deallocate(strategy)
        ! nparts without part: the distributed master
        call cline%set('nparts', 4)
        strategy = create_rec3D_strategy(cline)
        call assert_int(KIND_DISTR, strategy_kind(strategy), 'gridding with nparts is the distributed strategy')
        deallocate(strategy)
        ! nparts with part: a worker runs in memory
        call cline%set('part', 1)
        strategy = create_rec3D_strategy(cline)
        call assert_int(KIND_INMEM, strategy_kind(strategy), 'a gridding worker part runs the in-memory strategy')
        deallocate(strategy)
        call cline%delete('part')
        call cline%delete('nparts')
        ! pcg backend, shared memory
        call cline%set('rec_backend', 'pcg')
        strategy = create_rec3D_strategy(cline)
        call assert_int(KIND_PCG_INMEM, strategy_kind(strategy), 'pcg is the PCG in-memory strategy')
        deallocate(strategy)
        ! pcg with nparts > 1: distributed
        call cline%set('nparts', 2)
        strategy = create_rec3D_strategy(cline)
        call assert_int(KIND_DISTR, strategy_kind(strategy), 'pcg with nparts > 1 is the distributed strategy')
        deallocate(strategy)
        ! pcg with nparts = 1: stays in memory
        call cline%set('nparts', 1)
        strategy = create_rec3D_strategy(cline)
        call assert_int(KIND_PCG_INMEM, strategy_kind(strategy), 'pcg with nparts = 1 stays the PCG in-memory strategy')
        deallocate(strategy)
        ! pcg worker part
        call cline%set('nparts', 2)
        call cline%set('part', 2)
        strategy = create_rec3D_strategy(cline)
        call assert_int(KIND_PCG_INMEM, strategy_kind(strategy), 'a pcg worker part runs the PCG in-memory strategy')
        deallocate(strategy)
        call cline%kill
    end subroutine test_factory_branches

end submodule simple_rec3D_strategy_tester
