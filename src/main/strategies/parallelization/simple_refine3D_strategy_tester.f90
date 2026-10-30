!@descr: portable reference, FSC and accumulator staging checks for the refine3D startup strategy
module simple_refine3D_strategy_tester
use simple_string,     only: string
use simple_sp_project, only: sp_project
use simple_image,      only: image
use simple_syslib,     only: file_exists
use simple_test_utils, only: assert_true, assert_false, assert_int, assert_char, enter_fixture, leave_fixture, tests_failed
implicit none
private
public :: run_all_refine3D_staging_tests
integer, parameter :: BOX = 16
real, parameter    :: SMPD = 1.5

contains

    subroutine run_all_refine3D_staging_tests()
        type(string) :: saved, root
        integer      :: nfail0
        write(*,'(A)') '**** running all refine3D staging tests ****'
        nfail0 = tests_failed
        call enter_fixture('local_staging', saved, root)
        call test_startup_staging()
        call test_reference_transition()
        call leave_fixture(saved, root, nfail0)
    end subroutine run_all_refine3D_staging_tests

    ! An impulse has power in every shell: identical/opposite pairs have FSC +1/-1.
    subroutine test_startup_staging()
        use, intrinsic :: iso_fortran_env, only: int64
        use simple_parameters,             only: parameters
        use simple_fileio,                 only: file2rarr
        use simple_halfmap_diagnostics,    only: write_base_halfmap_fsc
        use simple_syslib,                 only: simple_abspath
        use simple_refine3D_strategy,      only: stage_refine3D_local_references
        use simple_refine3D_provenance,    only: refine3D_provenance
        use simple_sigma2_state_file,      only: sigma2_state_digest_file
        type(refine3D_provenance)      :: reference_record
        class(parameters), allocatable :: params
        class(sp_project), allocatable :: project
        class(image), allocatable      :: volume
        type(string)                   :: source(3), target, target_path, fsc_file, missing_fsc, registered
        integer(int64)                 :: digests(3)
        real, allocatable              :: fsc(:)
        real                           :: registered_smpd
        integer                        :: i, registered_box
        character(len=256)             :: message
        logical                        :: reused
        write(*,'(A)') 'test_startup_staging'
        allocate(params, project, volume)
        params%projfile = 'ownership_not_written.simple'
        params%l_nonuniform = .false.
        source = [string('ownership.mrc'), string('ownership_even.mrc'), string('ownership_odd.mrc')]
        fsc_file = 'ownership_fsc.bin'
        missing_fsc = 'ownership_missing_fsc.bin'
        call volume%new([BOX,BOX,BOX], SMPD)
        call volume%set_rmat_at(1, 1, 1, 1.)
        do i = 1,3
            call volume%write(source(i))
        enddo
        call write_base_halfmap_fsc(source(2), source(3), [BOX,BOX,BOX], SMPD, fsc_file)
        fsc = file2rarr(fsc_file)
        call assert_int(BOX/2, size(fsc), 'startup FSC: one value per nonzero shell')
        call assert_true(all(abs(fsc-1.) < 1.e-6), 'startup FSC: identical impulses give one in every shell')
        call volume%set_rmat_at(1, 1, 1, -1.)
        call volume%write(source(3))
        do i = 1,3
            digests(i) = sigma2_state_digest_file(source(i))
        enddo
        ! No project file exists: staging must use and update the supplied project in memory.
        call project%os_ptcl3D%new(1, is_ptcl=.true.)
        call reference_record%new(source(1), 'reference', message)
        call stage_refine3D_local_references(params, project, reference_record, missing_fsc, .false., target, reused)
        call assert_false(reused, 'startup staging: legacy pair takes the base-half fallback')
        call assert_true(project%isthere_in_osout('vol', 1), 'startup staging: caller receives the volume registration')
        call project%get_vol('vol', 1, registered, registered_smpd, registered_box)
        target_path = simple_abspath(target)
        call assert_char(target_path%to_char(), registered%to_char(), 'startup staging: caller volume is the staged volume')
        call assert_int(BOX, registered_box, 'startup staging: registered native box')
        call assert_true(abs(registered_smpd-SMPD) < 1.e-6, 'startup staging: registered native sampling')
        call project%get_fsc(1, fsc_file, registered_box)
        fsc = file2rarr(fsc_file)
        call assert_true(all(abs(fsc+1.) < 1.e-6), 'startup staging: missing FSC is rebuilt from opposite impulses')
        do i = 1,3
            call assert_true(sigma2_state_digest_file(source(i)) == digests(i), 'startup staging: source map bytes unchanged')
        enddo
        call assert_false(file_exists(params%projfile), 'startup staging: project publication remains with the caller')
        call test_chain_staging(params, project)
        params%l_ml_reg = .false.
        call test_reference_derivatives(params, project)
        call reference_record%kill
        call project%kill()
        call volume%kill()
    end subroutine test_startup_staging

    subroutine test_reference_derivatives( params, project )
        use, intrinsic :: iso_fortran_env, only: int64
        use simple_parameters,             only: parameters
        use simple_fileio,                 only: arr2file, file2rarr, add2fbody
        use simple_syslib,                 only: del_file
        use simple_defs_fname,             only: MRC_EXT, NUFILT_SUFFIX, NULOCRES_SUFFIX, AUTOMASK_FBODY, NU_ENVMASK_FBODY
        use simple_refine3D_provenance,    only: begin_reference_set, publish_reference_set, refine3D_provenance
        use simple_refine3D_strategy,      only: stage_reference_derivatives
        use simple_sigma2_state_file,      only: sigma2_state_digest_file
        class(parameters), intent(inout) :: params
        class(sp_project), intent(in)    :: project
        type(refine3D_provenance) :: reference_record
        type(string)              :: merged, halves(2), fsc, source(6), staged(6), target, name
        integer(int64)            :: digests(6)
        real, allocatable         :: values(:)
        real                      :: align_lp
        integer                   :: i
        logical                   :: reused
        character(len=256)        :: message
        write(*,'(A)') 'test_reference_derivatives'
        merged = 'nu_reference.mrc'
        halves = [string('nu_reference_even.mrc'), string('nu_reference_odd.mrc')]
        fsc = 'nu_fsc.bin'
        target = 'staged_reference.mrc'
        call arr2file([1.,0.], merged)
        do i = 1,2
            call arr2file([real(i),0.], halves(i))
        enddo
        call arr2file([1.,0.5], fsc)
        source(1) = add2fbody(merged, MRC_EXT, NUFILT_SUFFIX)
        source(2) = add2fbody(halves(1), MRC_EXT, NUFILT_SUFFIX)
        source(3) = add2fbody(halves(2), MRC_EXT, NUFILT_SUFFIX)
        source(4) = add2fbody(merged, MRC_EXT, NULOCRES_SUFFIX)
        source(5) = AUTOMASK_FBODY//'02'//MRC_EXT
        source(6) = NU_ENVMASK_FBODY//'02'//MRC_EXT
        staged(1) = add2fbody(target, MRC_EXT, NUFILT_SUFFIX)
        name = add2fbody(target, MRC_EXT, '_even')
        staged(2) = add2fbody(name, MRC_EXT, NUFILT_SUFFIX)
        name = add2fbody(target, MRC_EXT, '_odd')
        staged(3) = add2fbody(name, MRC_EXT, NUFILT_SUFFIX)
        staged(4) = add2fbody(target, MRC_EXT, NULOCRES_SUFFIX)
        staged(5) = AUTOMASK_FBODY//'01'//MRC_EXT
        staged(6) = NU_ENVMASK_FBODY//'01'//MRC_EXT
        do i = 1,6
            call arr2file([real(i),0.], source(i))
            digests(i) = sigma2_state_digest_file(source(i))
        enddo
        params%l_nonuniform = .true.
        params%filt_mode = 'nonuniform'
        params%automsk = 'nu'
        params%lp = 9.
        call begin_reference_set(merged)
        call publish_reference_set(params, project, project%os_ptcl3D, 2, merged, halves, fsc, 8.)
        ! The live LP is an output of matching, not a fixed preparation setting.
        params%lp = 20.
        call reference_record%new(merged, 'reference', message)
        reused = stage_reference_derivatives(params, reference_record, target, align_lp, message)
        call assert_true(reused, 'NU bundle: verified derivatives are staged')
        call assert_true(abs(align_lp - 8.) < 1.e-6, 'NU bundle: recorded matching LP is restored')
        do i = 1,6
            call assert_true(file_exists(staged(i)), 'NU bundle: recorded derivative exists at destination')
            if( .not. file_exists(staged(i)) ) cycle
            values = file2rarr(staged(i))
            call assert_true(size(values) == 2, 'NU bundle: destination shape')
            if( size(values) == 2 ) call assert_true(all(values == [real(i),0.]), 'NU bundle: copied independent values')
            call assert_true(sigma2_state_digest_file(source(i)) == digests(i), 'NU bundle: source bytes preserved')
        enddo
        call reference_record%new(merged, 'reference', message)
        reused = stage_reference_derivatives(params, reference_record, merged, align_lp, message)
        call assert_true(reused, 'NU bundle: in-place map copies are accepted')
        do i = 1,6
            call assert_true(sigma2_state_digest_file(source(i)) == digests(i), 'NU bundle: in-place copies preserve bytes')
        enddo
        call assert_true(file_exists(staged(6)), 'NU bundle: old destination envelope exists before staging')
        call del_file(source(6))
        call begin_reference_set(merged)
        call publish_reference_set(params, project, project%os_ptcl3D, 2, merged, halves, fsc, 8.)
        call reference_record%new(merged, 'reference', message)
        reused = stage_reference_derivatives(params, reference_record, target, align_lp, message)
        call assert_true(reused, 'NU bundle: missing recorded envelope still stages verified constituents')
        call assert_false(file_exists(staged(6)), 'NU bundle: unrecorded destination envelope is removed')
        call assert_true(file_exists(staged(5)), 'NU bundle: recorded destination mask is retained')
        call assert_true(sigma2_state_digest_file(source(5)) == digests(5), 'NU bundle: recorded source mask is preserved')
        params%amsklp = params%amsklp + 1.
        call reference_record%new(merged, 'reference', message)
        reused = stage_reference_derivatives(params, reference_record, string('incompatible.mrc'), align_lp, message)
        call assert_false(reused, 'NU bundle: changed preparation settings use the fallback')
        name = add2fbody(string('incompatible.mrc'), MRC_EXT, NUFILT_SUFFIX)
        call assert_false(file_exists(name), 'NU bundle: rejected derivatives are not copied')
        params%amsklp = params%amsklp - 1.
        call begin_reference_set(merged)
        call publish_reference_set(params, project, project%os_ptcl3D, 2, merged, halves, fsc, 0.)
        call reference_record%new(merged, 'reference', message)
        reused = stage_reference_derivatives(params, reference_record, target, align_lp, message)
        call assert_false(reused, 'NU bundle: absent matching LP uses the fallback')
        params%l_nonuniform = .false.
        call begin_reference_set(merged)
        call publish_reference_set(params, project, project%os_ptcl3D, 2, merged, halves, fsc, 8.)
        params%l_nonuniform = .true.
        call reference_record%new(merged, 'reference', message)
        reused = stage_reference_derivatives(params, reference_record, target, align_lp, message)
        call assert_false(reused, 'NU bundle: unrecorded existing derivatives are not trusted')
        call reference_record%new(string('legacy.mrc'), 'reference', message)
        reused = stage_reference_derivatives(params, reference_record, target, align_lp, message)
        call assert_false(reused, 'NU bundle: absent provenance uses the fallback')
        call reference_record%kill
    end subroutine test_reference_derivatives

    subroutine test_chain_staging( params, project )
        use, intrinsic :: iso_fortran_env, only: int64
        use simple_parameters,             only: parameters
        use simple_fileio,                 only: arr2file
        use simple_refine3D_strategy,      only: stage_refine3D_local_chain
        use simple_refine3D_fnames,        only: refine3D_trail_rec_fname, refine3D_trail_rho_fname, &
            &refine3D_trail_manifest_fname, refine3D_pcg_trail_accum_fname, refine3D_pcg_trail_manifest_fname
        use simple_trail_chain_manifest, only: trail_chain_manifest
        use simple_refine3D_provenance,  only: publish_reuse_provenance, reuse_provenance_matches
        use simple_sigma2_state_file,    only: sigma2_state_digest_file
        class(parameters), intent(inout) :: params
        class(sp_project), intent(inout) :: project
        type(trail_chain_manifest) :: manifest
        type(string)               :: files(5), anchor
        integer(int64)             :: sizes(4), digests(5)
        integer                    :: backend, ncomponents, i, status
        logical                    :: valid
        character(len=256)         :: message
        write(*,'(A)') 'test_chain_staging'
        params%l_update_frac = .true.
        params%l_trail_rec = .true.
        params%l_m_estimator_flex = .false.
        params%oritype = 'ptcl3D'
        params%box_crop = BOX
        params%smpd_crop = SMPD
        call project%projinfo%new(1, is_ptcl=.false.)
        call project%projinfo%set(1, 'projname', 'ownership')
        call project%os_stk%new(1, is_ptcl=.false.)
        call project%os_stk%set(1, 'stk', 'ownership_particles.mrc')
        call project%os_ptcl3D%set(1, 'stkind', 1)
        call project%os_ptcl3D%set(1, 'indstk', 1)
        call project%os_ptcl3D%set_state(1, 1)
        call project%os_ptcl3D%set(1, 'eo', 0)
        call project%os_ptcl3D%set(1, 'updatecnt', 1)
        call project%os_ptcl3D%set_euler(1, [0.,0.,0.])
        call project%os_ptcl3D%set_shift(1, [0.,0.])
        do backend = 1,2
            if( backend == 1 )then
                params%rec_backend = 'gridding'
                files(1) = refine3D_trail_manifest_fname(1)
                files(2) = refine3D_trail_rec_fname(1, 'even')
                files(3) = refine3D_trail_rho_fname(1, 'even')
                files(4) = refine3D_trail_rec_fname(1, 'odd')
                files(5) = refine3D_trail_rho_fname(1, 'odd')
                ncomponents = 4
            else
                params%rec_backend = 'pcg'
                files(1) = refine3D_pcg_trail_manifest_fname(1)
                files(2) = refine3D_pcg_trail_accum_fname(1, 'even')
                files(3) = refine3D_pcg_trail_accum_fname(1, 'odd')
                ncomponents = 2
            endif
            sizes = 0_int64
            do i = 1,ncomponents
                call arr2file([real(i),0.], files(i+1))
                inquire(file=files(i+1)%to_char(), size=sizes(i))
            enddo
            call manifest%new(BOX, SMPD, 1, 1, 1, 1, sizes, 1., 1, [0_int64,0_int64])
            call manifest%write(files(1), status)
            call assert_int(0, status, 'startup chain: fixture manifest published')
            call publish_reuse_provenance(params, project, project%os_ptcl3D, 'chain', files(1), files(:ncomponents+1))
            valid = reuse_provenance_matches(params, project, project%os_ptcl3D, 'chain', &
                &files(1), files(:ncomponents+1), message)
            call assert_true(valid, 'startup chain: fixture provenance is consumable: '//trim(message))
            if( .not. valid ) cycle
            do i = 1,ncomponents+1
                digests(i) = sigma2_state_digest_file(files(i))
            enddo
            anchor = 'ownership.mrc'
            valid = stage_refine3D_local_chain(params, project, anchor, BOX, SMPD)
            call assert_true(valid, 'startup chain: compatible in-place chain accepted')
            do i = 1,ncomponents+1
                call assert_true(sigma2_state_digest_file(files(i)) == digests(i), 'startup chain: in-place file bytes unchanged')
            enddo
            call project%os_ptcl3D%set_euler(1, [1.,0.,0.])
            valid = stage_refine3D_local_chain(params, project, anchor, BOX, SMPD)
            call assert_false(valid, 'startup chain: changed producer poses reject reuse')
            call project%os_ptcl3D%set_euler(1, [0.,0.,0.])
            call assert_false(file_exists(params%projfile), 'startup chain: no project file read or write needed')
        enddo
        call manifest%kill()
    end subroutine test_chain_staging

    ! Known impulses distinguish startup maps, assembled raw halves and assembled NU halves.
    subroutine test_reference_transition()
        use simple_parameters,           only: parameters
        use simple_cmdline,              only: cmdline
        use simple_builder,              only: builder
        use simple_fileio,               only: add2fbody
        use simple_defs_fname,           only: MRC_EXT, NUFILT_SUFFIX
        use simple_refine3D_strategy,    only: adopt_refine3D_reference
        use simple_matcher_refvol_utils, only: read_mask_filter_reproject_refvols
        class(parameters), allocatable      :: params
        class(cmdline), allocatable         :: cline
        class(builder), allocatable, target :: build
        class(image), allocatable           :: volume
        type(string) :: merged, even, odd, nu_even, nu_odd, active
        integer      :: generation, center
        real         :: expected(BOX,BOX,BOX), even_value, odd_value
        write(*,'(A)') 'test_reference_transition'
        allocate(params, cline, build, volume)
        params%box               = BOX
        params%box_crop          = BOX
        params%smpd              = SMPD
        params%smpd_crop         = SMPD
        params%nstates           = 1
        params%msk_crop          = real(BOX/2-2)
        params%kfromto           = [1,BOX/2-1]
        params%lp                = 6.
        params%center            = 'no'
        params%automsk           = 'no'
        params%filt_mode         = 'nonuniform'
        params%l_nonuniform      = .true.
        params%l_nonuniform_lpset= .false.
        params%l_cart_refine     = .true.
        params%l_lpset           = .true.
        params%l_noise_reg       = .false.
        params%l_icm             = .false.
        params%l_ml_reg          = .true.
        call build%spproj%os_ptcl3D%new(1, is_ptcl=.true.)
        call build%spproj%os_ptcl3D%set(1, 'state', 1.)
        call build%spproj%os_ptcl3D%set(1, 'res', 6.)
        build%spproj_field => build%spproj%os_ptcl3D
        allocate(build%fsc(1,BOX/2), source=1.)
        call build%new_vols(params)
        call volume%new([BOX,BOX,BOX], SMPD)
        center = BOX/2+1
        merged = 'transition_startvol.mrc'
        call write_triplet(-10., -11.)
        params%vols(1)           = merged
        params%vols_even(1)      = even
        params%vols_odd(1)       = odd
        params%vols(2)           = 'untouched_state.mrc'
        call cline%set('vol_even', even)
        call cline%set('vol_odd', odd)
        do generation = 1,2
            if( generation == 1 )then
                merged = 'transition_recvol.mrc'
            else
                merged = 'transition_next_recvol.mrc'
            endif
            even_value = real(2*generation)
            odd_value  = even_value+1.
            call write_triplet(even_value, odd_value)
            call adopt_refine3D_reference(params, cline, 1, merged)
            call assert_char(merged%to_char(), params%vols(1)%to_char(), 'reference transition: active merged map')
            call assert_char(even%to_char(), params%vols_even(1)%to_char(), 'reference transition: active even map')
            call assert_char(odd%to_char(), params%vols_odd(1)%to_char(), 'reference transition: active odd map')
            active = cline%get_carg('vol1')
            call assert_char(merged%to_char(), active%to_char(), 'reference transition: child volume matches master')
            call assert_false(cline%defined('vol_even'), 'reference transition: obsolete explicit even input removed')
            call assert_false(cline%defined('vol_odd'), 'reference transition: obsolete explicit odd input removed')
            call assert_char('untouched_state.mrc', params%vols(2)%to_char(), 'reference transition: other state unchanged')
            ! Exercise the reader of the next iteration, with no parameter reconstruction.
            call read_mask_filter_reproject_refvols(params, build, cline, map_shift=.false.)
            expected = 0.
            expected(center,center,center) = even_value
            call assert_true(maxval(abs(build%vol%get_rmat()-expected)) < 1.e-5, &
                &'reference transition: next iteration consumes assembled even NU content')
            expected(center,center,center) = odd_value
            call assert_true(maxval(abs(build%vol_odd%get_rmat()-expected)) < 1.e-5, &
                &'reference transition: next iteration consumes assembled odd NU content')
            call build%cftc%kill()
        enddo
        call volume%kill()
        call build%vol%kill()
        call build%vol_odd%kill()
        call build%vol2%kill()
        nullify(build%spproj_field)
        call build%spproj%kill()
        call cline%kill()
    contains
        subroutine write_triplet( even_amplitude, odd_amplitude )
            real, intent(in) :: even_amplitude, odd_amplitude
            even    = add2fbody(merged, MRC_EXT, '_even')
            odd     = add2fbody(merged, MRC_EXT, '_odd')
            nu_even = add2fbody(even, MRC_EXT, NUFILT_SUFFIX)
            nu_odd  = add2fbody(odd, MRC_EXT, NUFILT_SUFFIX)
            call volume%zero()
            call volume%set_rmat_at(center, center, center, 40.)
            call volume%write(merged)
            call volume%write(even)
            call volume%set_rmat_at(center, center, center, 50.)
            call volume%write(odd)
            call volume%set_rmat_at(center, center, center, even_amplitude)
            call volume%write(nu_even)
            call volume%set_rmat_at(center, center, center, odd_amplitude)
            call volume%write(nu_odd)
        end subroutine write_triplet
    end subroutine test_reference_transition

end module simple_refine3D_strategy_tester
