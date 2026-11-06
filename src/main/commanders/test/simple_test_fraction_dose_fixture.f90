!@descr: fraction-dose fixtures shared by testers: a saved motion model with a fraction dose, and the check of a fraction export to its dose-series manifest
module simple_test_fraction_dose_fixture
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_defs, only: dp, CWD_GLOB
use simple_type_defs, only: OBJFUN_CC
use simple_string, only: string
use simple_string_utils, only: int2str, int2str_pad
use simple_syslib, only: simple_chdir, simple_mkdir
use simple_fileio, only: file_exists, simple_copy_file, file2rarr
use simple_parameters, only: parameters
use simple_cmdline, only: cmdline
use simple_sp_project, only: sp_project
use simple_motion_model, only: motion_model
use simple_image, only: image
use simple_dose_series_manifest, only: dose_series_manifest
use simple_refine_motion_model_strategy, only: reconstruct_motion_fractions
use simple_test_utils, only: assert_true, assert_false, assert_int
implicit none
private
public :: check_fraction_dose_export, write_saved_model
public :: MOVIE_NFRAMES, RECON_NPTCLS, SAVED_FRACTION_DOSE

integer, parameter :: MOVIE_NFRAMES = 12, RECON_NPTCLS = 120
real,    parameter :: SAVED_FRACTION_DOSE = 1.25

contains

    !> Persist realistic model metadata with an EER grouping and target distinct from the saved fraction dose.
    subroutine write_saved_model(params, fname, nframes, dose, kv)
        class(parameters), target, intent(inout) :: params
        type(string), intent(in) :: fname
        integer, intent(in) :: nframes
        real, intent(in) :: dose, kv
        class(motion_model), allocatable :: model
        class(image), pointer :: frames(:)
        type(string) :: movie
        integer :: i, funit
        allocate(model, frames(nframes))
        movie = fname//'.eer'
        open(newunit=funit, file=movie%to_char(), status='replace', access='stream', form='unformatted')
        write(funit) 0
        close(funit)
        do i = 1, nframes
            call frames(i)%new([4,4,1], 1.5, wthreads=.false.)
        enddo
        params%fraction_dose_target = 9.
        params%eer_upsampling = 1
        params%total_dose = real(nframes)*merge(dose,0.,ieee_is_finite(dose))
        call model%new(params, movie, [4,4], 1.5, frames, 8*nframes, 1, kv, merge(dose,0.,ieee_is_finite(dose)), 8)
        model%dose_per_frame = dose
        call model%write(fname//'.star', fname, .false.)
        call model%kill
        do i = 1, nframes
            call frames(i)%kill
        enddo
        deallocate(frames)
    end subroutine write_saved_model

    !> Verify success and metadata-only rejection through the production reconstruction loop.
    subroutine check_fraction_dose_export(root, name, stack, oritab, valid_dose, changed_state)
        use simple_test_flex_phantom_fixture, only: build_flex_pca_phantom_project, &
            &FLEX_PHANTOM_BOX, FLEX_PHANTOM_SMPD, FLEX_PHANTOM_MSKDIAM
        type(string), intent(in) :: root, stack, oritab
        character(len=*), intent(in) :: name
        logical, intent(in) :: valid_dose, changed_state
        class(parameters), allocatable :: params
        class(cmdline), allocatable :: cline
        class(sp_project), allocatable :: project
        class(dose_series_manifest), allocatable :: manifest
        type(string), allocatable :: projects(:), stacks(:)
        integer, allocatable :: ranges(:,:)
        type(string) :: arm, directory, model, message, even, odd
        real, allocatable :: fsc(:)
        real(dp) :: dlo, dhi
        integer :: i, j, nwin, nstates, a, b
        logical :: ok
        allocate(params, cline, project, manifest)
        arm = root//'/'//name
        call simple_mkdir(arm)
        call simple_chdir(arm)
        CWD_GLOB = arm%to_char()
        params%box = FLEX_PHANTOM_BOX
        params%smpd = FLEX_PHANTOM_SMPD
        params%mskdiam = FLEX_PHANTOM_MSKDIAM
        params%nparts = 1
        params%nthr = 1
        params%pgrp = 'c1'
        params%cc_objfun = OBJFUN_CC
        params%l_ml_reg = .false.
        params%cc_emit_sigma = 'no'
        model = arm//'/motion.bin'
        call write_saved_model(params, model, MOVIE_NFRAMES, merge(SAVED_FRACTION_DOSE,0.,valid_dose), 300.)
        nwin = 2
        if( valid_dose .and. .not. changed_state ) nwin = 3
        nstates = 1
        if( changed_state ) nstates = 2
        allocate(projects(nwin), stacks(nwin), ranges(2,nwin))
        ranges(:,1) = [1,3]
        ranges(:,2) = [4,6]
        if( nwin == 3 ) ranges(:,3) = [6,8]
        do i = 1, nwin
            directory = arm//'/frames_'//int2str(ranges(1,i))//'-'//int2str(ranges(2,i))
            call simple_mkdir(directory)
            stacks(i) = directory//'/particles.mrcs'
            call simple_copy_file(stack, stacks(i))
            projects(i) = directory//'/fraction.simple'
            call build_flex_pca_phantom_project(projects(i), stacks(i), oritab, nptcls=RECON_NPTCLS)
            call project%read(projects(i))
            call project%os_mic%new(1, is_ptcl=.false.)
            call project%os_mic%set(1, 'mcmodel', model)
            call project%os_mic%set(1, 'nptcls', RECON_NPTCLS)
            call project%os_mic%set_state(1, 1)
            if( changed_state )then
                do j = RECON_NPTCLS-3, RECON_NPTCLS
                    call project%os_ptcl3D%set_state(j, 2)
                enddo
                if( i == 2 )then
                    call assert_int(project%os_ptcl3D%get_eo(1), project%os_ptcl3D%get_eo(RECON_NPTCLS-3), &
                        &'fraction export: state-membership swap keeps the half population')
                    call project%os_ptcl3D%set_state(1, 2)
                    call project%os_ptcl3D%set_state(RECON_NPTCLS-3, 1)
                endif
            endif
            call project%write(projects(i))
            call project%kill
        enddo
        call reconstruct_motion_fractions(params, cline, projects, ranges, nstates)
        call assert_true(file_exists(arm//'/fsc_fractions_state01.tsv'), 'fraction export '//name//': raw FSC matrix retained')
        do i = 1, nwin
            call assert_true(file_exists(stacks(i)), 'fraction export '//name//': extracted stack retained')
            directory = arm//'/rec_frames_'//int2str_pad(ranges(1,i),2)//'-'//int2str_pad(ranges(2,i),2)
            call assert_true(file_exists(directory//'/recvol_state01_even_unfil.mrc'), 'fraction export '//name//': raw even half')
            call assert_true(file_exists(directory//'/recvol_state01_odd_unfil.mrc'), 'fraction export '//name//': raw odd half')
            fsc = file2rarr(directory//'/fsc_state01.bin')
            call assert_int(FLEX_PHANTOM_BOX/2, size(fsc), 'fraction export '//name//': native FSC shell count')
        enddo
        if( .not. valid_dose .or. changed_state )then
            call assert_false(file_exists(arm//'/dose_series.txt'), 'fraction export '//name//': no invalid manifest published')
        else
            ok = manifest%read(arm//'/dose_series.txt', message)
            call assert_true(ok, 'fraction export: generated manifest is readable')
            if( ok )then
                call assert_int(3, manifest%get_nwindows(), 'fraction export: producer retains all actual intervals')
                ok = manifest%select_windows(message)
                call assert_true(ok, 'fraction export: manifest selects disjoint intervals')
                call assert_int(2, manifest%get_nwindows(), 'fraction export: overlapping predecessor removed')
                call manifest%get_window(2, a, b, dlo, dhi, even, odd)
                call assert_int(6, a, 'fraction export: highest-dose interval retained')
                call assert_true(abs(dlo-6.25_dp) < 1.e-12_dp, 'fraction export: final interval physical dose origin')
                call assert_true(abs(dhi-10._dp) < 1.e-12_dp, 'fraction export: final interval physical dose end')

            endif
        endif
        call manifest%kill
        call simple_chdir(root)
        CWD_GLOB = root%to_char()
    end subroutine check_fraction_dose_export

end module simple_test_fraction_dose_fixture
