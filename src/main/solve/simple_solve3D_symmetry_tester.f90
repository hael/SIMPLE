!@descr: library regression of the solve3D symmetry-axis handoff and the next fractional trailing reconstruction
module simple_solve3D_symmetry_tester
use simple_test_utils
!$ use omp_lib,             only: omp_get_max_threads, omp_set_num_threads
use simple_string,           only: string
use simple_cmdline,          only: cmdline
use simple_parameters,       only: parameters
use simple_sp_project,       only: sp_project
use simple_ori,              only: ori
use simple_image,            only: image
use simple_commanders_rec,   only: commander_rec3D
use simple_syslib,           only: simple_chdir, simple_mkdir, file_exists
use simple_fileio,           only: simple_copy_file, files_identical
use simple_type_defs,        only: lp_crop_inf
use simple_refine3D_fnames,  only: refine3D_state_vol_fname, refine3D_trail_manifest_fname, &
    &refine3D_trail_rec_fname, refine3D_trail_rho_fname
use simple_trail_chain_manifest, only: trail_chain_manifest, TRAIL_MANIFEST_OK
use simple_refine3D_provenance,  only: reuse_provenance_matches
use simple_solve3D_utils, only: symmetrize, lpinfo, cline_refine3D, cline_reconstruct3D, cline_symmap, &
    &l_srch4symaxis, l_symran, l_cavgs_mode, nptcls_eff
use simple_test_truth_metrics, only: axis_operator
use simple_test_flex_phantom_fixture, only: create_flex_pca_phantom_fixture, build_flex_pca_phantom_project, &
    &FLEX_PHANTOM_BOX, FLEX_PHANTOM_SMPD, FLEX_PHANTOM_MSKDIAM
implicit none
private
public :: run_all_solve3D_symmetry_tests

integer, parameter :: NPTCLS = 128, BOX = FLEX_PHANTOM_BOX
real, parameter :: SMPD = FLEX_PHANTOM_SMPD, MSKDIAM = FLEX_PHANTOM_MSKDIAM
real, parameter :: MAP_RELTOL = 256. * epsilon(1.) !< single-precision gridding and restoration roundoff

contains

    subroutine run_all_solve3D_symmetry_tests()
        write(*,'(A)') '**** running solve3D symmetry handoff tests ****'
        call test_symmetry_trailing_handoff()
    end subroutine run_all_solve3D_symmetry_tests

    !> A chain made before a nonidentity symmetry alignment must not enter the next fractional blend. Covered:
    !! one state, gridding, equal stage crops; the stale old-frame chain is the negative control, and a run whose
    !! next stage does not trail must leave no chain at all.
    subroutine test_symmetry_trailing_handoff()
        class(parameters), allocatable :: params
        class(sp_project), allocatable :: project
        class(cmdline), allocatable :: saved_refine, saved_rec, saved_symmap
        type(lp_crop_inf), allocatable :: saved_lpinfo(:)
        type(commander_rec3D) :: xrec3D
        type(trail_chain_manifest) :: manifest
        type(ori) :: before, after
        type(string) :: cwd, root, truth_a, truth_b, truth_mean, truth_diff, oritab, stack, projfile, files(5)
        logical :: saved_flags(3), valid
        integer :: nfail0, i, status, nthr_saved, saved_nptcls_eff
        character(len=256) :: message
        write(*,'(A)') 'test_symmetry_trailing_handoff'
        nfail0 = tests_failed
        nthr_saved = 1
        !$ nthr_saved = omp_get_max_threads()
        call enter_fixture('solve3D_symmetry_handoff', cwd, root)
        allocate(params, project, saved_refine, saved_rec, saved_symmap)
        saved_refine = cline_refine3D
        saved_rec = cline_reconstruct3D
        saved_symmap = cline_symmap
        saved_flags = [l_srch4symaxis, l_symran, l_cavgs_mode]
        saved_nptcls_eff = nptcls_eff
        nptcls_eff = 0
        call move_alloc(lpinfo, saved_lpinfo)
        call create_flex_pca_phantom_fixture(root, 1, truth_a, truth_b, truth_mean, truth_diff, oritab, stack, nptcls=NPTCLS)
        projfile = root//'/phantom.simple'
        call build_flex_pca_phantom_project(projfile, stack, oritab, nptcls=NPTCLS)
        call project%read(projfile)
        call project%os_ptcl3D%set_all2single('updatecnt', 1.)
        call project%os_ptcl3D%set_all2single('sampled', 1.)
        call project%write_segment_inside('ptcl3D', projfile)
        call project%os_ptcl3D%get_ori(1, before)
        call reconstruct(projfile, 'c1', .true., .false.)
        call simple_copy_file(refine3D_trail_rec_fname(1, 'even'), string('old_frame_even.mrc'))
        ! the old-frame chain, kept for the negative control
        files = [refine3D_trail_manifest_fname(1), refine3D_trail_rec_fname(1, 'even'), &
            &refine3D_trail_rho_fname(1, 'even'), refine3D_trail_rec_fname(1, 'odd'), refine3D_trail_rho_fname(1, 'odd')]
        call simple_mkdir(string('stale'))
        do i = 1, size(files)
            call simple_copy_file(files(i), string('stale/')//files(i))
        enddo
        call write_turned_c2_map()

        params%box = BOX
        params%box_crop = BOX
        params%smpd = SMPD
        params%smpd_crop = SMPD
        params%mskdiam = MSKDIAM
        params%nstates = 1
        params%pgrp = 'c2'
        params%pgrp_start = 'c1'
        params%projfile = projfile
        allocate(lpinfo(4))
        lpinfo%lp = 8.
        lpinfo%box_crop = BOX
        lpinfo%smpd_crop = SMPD
        l_srch4symaxis = .true.
        l_symran = .false.
        l_cavgs_mode = .false.
        call cline_refine3D%kill
        call cline_refine3D%set('trail_rec', 'yes')
        call cline_refine3D%set('rec_backend', 'gridding')
        call cline_refine3D%set('objfun', 'cc')
        call cline_refine3D%set('ml_reg', 'no')
        call cline_refine3D%set('endit', 10)
        call cline_reconstruct3D%kill
        call reconstruction_cline(cline_reconstruct3D, projfile, 'c1', .false., .false.)
        call cline_symmap%kill
        call cline_symmap%set('prg', 'symmetrize_map')
        call cline_symmap%set('mkdir', 'no')
        call cline_symmap%set('projfile', projfile)
        call cline_symmap%set('oritype', 'ptcl3D')
        call cline_symmap%set('pgrp', 'c2')
        call cline_symmap%set('mskdiam', MSKDIAM)
        call cline_symmap%set('nthr', 1)
        call symmetrize(params, 4, project, projfile, xrec3D)
        call project%read(projfile)
        call project%os_ptcl3D%get_ori(1, after)
        call assert_true(maxval(abs(after%get_mat() - before%get_mat())) > 1.e-3, &
            &'the handoff applies a rotation well above Euler roundoff (precondition: the frame changes)')
        call manifest%read(refine3D_trail_manifest_fname(1), status)
        call assert_int(TRAIL_MANIFEST_OK, status, 'the symmetry reconstruction publishes a trailing seed')
        call assert_int(1, manifest%get_gen(), 'the new frame starts a fresh chain')
        call assert_real(real(NPTCLS), manifest%get_mrep(), 0., 'the new seed represents the full updated population')
        call assert_false(files_identical(string('old_frame_even.mrc'), refine3D_trail_rec_fname(1, 'even')), &
            &'the symmetric reconstruction replaces the old-frame accumulators')
        files = [refine3D_trail_manifest_fname(1), refine3D_trail_rec_fname(1, 'even'), &
            &refine3D_trail_rho_fname(1, 'even'), refine3D_trail_rec_fname(1, 'odd'), refine3D_trail_rho_fname(1, 'odd')]
        valid = reuse_provenance_matches(params, project, project%os_ptcl3D, 'chain', files(1), files, message)
        call assert_true(valid, 'the seed provenance describes the transformed particle poses: '//trim(message))

        ! Select half the poses, with particles from both half sets and conformers.
        do i = 1, NPTCLS
            if( mod((i-1)/4, 2) /= 0 ) cycle
            call project%os_ptcl3D%set(i, 'updatecnt', 2.)
            call project%os_ptcl3D%set(i, 'sampled', 2.)
        enddo
        call project%write_segment_inside('ptcl3D', projfile)
        call simple_mkdir(string('fresh'))
        call simple_copy_file(projfile, string('fresh/phantom.simple'))
        call simple_copy_file(projfile, string('stale/phantom.simple'))
        call reconstruct(projfile, 'c2', .false., .true.)
        call simple_copy_file(refine3D_state_vol_fname(1), string('actual_next.mrc'))
        call manifest%read(refine3D_trail_manifest_fname(1), status)
        call assert_int(TRAIL_MANIFEST_OK, status, 'the next fractional update publishes its chain')
        call assert_int(2, manifest%get_gen(), 'the next fractional update blends with the new-frame seed')

        ! Independent restart: reconstruct all transformed poses without any pre-alignment chain, then update.
        call simple_chdir(string('fresh'))
        call reconstruct(root//'/fresh/phantom.simple', 'c2', .true., .false.)
        call reconstruct(root//'/fresh/phantom.simple', 'c2', .false., .true.)
        call assert_true(map_reldiff(root//'/actual_next.mrc', refine3D_state_vol_fname(1)) <= MAP_RELTOL, &
            &'the next fractional map equals a clean reconstruction and blend in the transformed frame')
        call simple_chdir(root)

        ! Negative control: the old-frame chain, had it survived, is accepted and blended, and gives another map.
        call simple_chdir(string('stale'))
        call reconstruct(root//'/stale/phantom.simple', 'c2', .false., .true.)
        call manifest%read(refine3D_trail_manifest_fname(1), status)
        call assert_int(TRAIL_MANIFEST_OK, status, 'the stale chain is read')
        call assert_int(2, manifest%get_gen(), 'the stale chain is blended, not discarded (the control is live)')
        call assert_true(map_reldiff(refine3D_state_vol_fname(1), root//'/fresh/'//refine3D_state_vol_fname(1)) > &
            &1000. * MAP_RELTOL, 'a blend with the old-frame chain differs from the transformed-frame restart')
        call simple_chdir(root)

        ! A next stage that does not trail (full sampling) leaves no chain: the old one is removed, none is seeded.
        do i = 1, size(files)
            call simple_copy_file(string('stale/')//files(i), files(i))
        enddo
        call write_turned_c2_map()
        nptcls_eff     = NPTCLS
        params%nsample = NPTCLS
        call symmetrize(params, 4, project, projfile, xrec3D)
        do i = 1, size(files)
            call assert_false(file_exists(files(i)), 'a non-trailing next stage leaves no chain: '//files(i)%to_char())
        enddo

        call cline_refine3D%kill
        call cline_reconstruct3D%kill
        call cline_symmap%kill
        cline_refine3D = saved_refine
        cline_reconstruct3D = saved_rec
        cline_symmap = saved_symmap
        deallocate(lpinfo)
        call move_alloc(saved_lpinfo, lpinfo)
        l_srch4symaxis = saved_flags(1)
        l_symran = saved_flags(2)
        l_cavgs_mode = saved_flags(3)
        nptcls_eff   = saved_nptcls_eff
        call project%kill
        call before%kill
        call after%kill
        call manifest%kill
        !$ call omp_set_num_threads(nthr_saved)
        call leave_fixture(cwd, root, nfail0)
    end subroutine test_symmetry_trailing_handoff

    subroutine reconstruction_cline( cline, projfile, pgrp, seed, trail )
        class(cmdline), intent(inout) :: cline
        class(string), intent(in) :: projfile
        character(len=*), intent(in) :: pgrp
        logical, intent(in) :: seed, trail
        call cline%set('prg', 'reconstruct3D')
        call cline%set('mkdir', 'no')
        call cline%set('projfile', projfile)
        call cline%set('nstates', 1)
        call cline%set('pgrp', pgrp)
        call cline%set('mskdiam', MSKDIAM)
        call cline%set('rec_backend', 'gridding')
        call cline%set('objfun', 'cc')
        call cline%set('ml_reg', 'no')
        call cline%set('postprocess', 'no')
        call cline%set('nthr', 1)
        call cline%set('trail_rec', 'no')
        if( seed ) call cline%set('trail_seed', 'yes')
        if( trail )then
            call cline%set('trail_rec', 'yes')
            call cline%set('update_frac', 0.5)
        endif
    end subroutine reconstruction_cline

    subroutine reconstruct( projfile, pgrp, seed, trail )
        class(string), intent(in) :: projfile
        character(len=*), intent(in) :: pgrp
        logical, intent(in) :: seed, trail
        class(cmdline), allocatable :: cline
        type(commander_rec3D) :: xrec3D
        allocate(cline)
        call reconstruction_cline(cline, projfile, pgrp, seed, trail)
        call xrec3D%execute(cline)
        call cline%kill
    end subroutine reconstruct

    !> Analytic C2 blobs about a rotated axis, offset from the box centre.
    subroutine write_turned_c2_map()
        real, parameter :: POS(3,3) = reshape([14.,4.,3., -5.,12.,-10., 9.,-7.,14.], [3,3])
        real, parameter :: AMPS(3) = [1.,0.7,0.9], OFFSET(3) = [4.,-3.,2.], SIGMA = 3.
        type(image) :: vol
        real, allocatable :: density(:,:,:)
        real :: rot(3,3), centers(3,6), point(3)
        integer :: i, j, k, b, pair
        rot = axis_operator([1.,2.,3.] / norm2([1.,2.,3.]), 37.)
        do pair = 0, 1
            do b = 1, 3
                point = POS(:,b)
                if( pair == 1 ) point(1:2) = -point(1:2)
                centers(:,3*pair+b) = matmul(rot, point) + OFFSET
            enddo
        enddo
        allocate(density(BOX,BOX,BOX), source=0.)
        do k = 1, BOX
            do j = 1, BOX
                do i = 1, BOX
                    point = (real([i,j,k]) - real(BOX/2+1)) * SMPD
                    do b = 1, 6
                        density(i,j,k) = density(i,j,k) + AMPS(1+mod(b-1,3)) * &
                            &exp(-sum((point-centers(:,b))**2) / (2.*SIGMA**2))
                    enddo
                enddo
            enddo
        enddo
        call vol%new([BOX,BOX,BOX], SMPD, wthreads=.false.)
        call vol%set_rmat(density, .false.)
        call vol%write(refine3D_state_vol_fname(1), del_if_exists=.true.)
        call vol%kill
    end subroutine write_turned_c2_map

    !> largest voxel difference of two maps relative to the largest magnitude of the second
    real function map_reldiff( actual, expected ) result( reldiff )
        class(string), intent(in) :: actual, expected
        type(image) :: a, b
        real, allocatable :: ra(:,:,:), rb(:,:,:)
        call a%new([BOX,BOX,BOX], SMPD, wthreads=.false.)
        call b%new([BOX,BOX,BOX], SMPD, wthreads=.false.)
        call a%read(actual)
        call b%read(expected)
        ra = a%get_rmat()
        rb = b%get_rmat()
        reldiff = maxval(abs(ra-rb)) / max(maxval(abs(rb)), tiny(1.))
        call a%kill
        call b%kill
    end function map_reldiff

end module simple_solve3D_symmetry_tester
