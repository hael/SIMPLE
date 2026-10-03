!@descr: unit tests for the identity of the 3D trailing chains: gridding manifest and PCG chain header
! Both chains record the population they represent, M(s): the gridding manifest in a versioned
! field, the PCG chain in the particle count of each half's raw header. A chain of any other
! version is refused and re-seeded.
module simple_trail_chain_manifest_tester
use simple_core_module_api
use simple_test_utils
use simple_trail_chain_manifest, only: trail_chain_manifest, TRAIL_MANIFEST_OK, TRAIL_MANIFEST_MISSING, &
                                      &TRAIL_MANIFEST_UNREADABLE
use simple_reconstructor_pcg,    only: reconstructor_pcg, pcg_raw_accum_compatible, read_pcg_raw_accum_header
implicit none
private
public :: run_all_trail_chain_manifest_tests

contains

    subroutine run_all_trail_chain_manifest_tests()
        write(*,'(A)') '**** running all trailing chain identity tests ****'
        call test_gridding_manifest()
        call test_pcg_chain_header()
    end subroutine run_all_trail_chain_manifest_tests

    subroutine test_gridding_manifest()
        type(trail_chain_manifest) :: man, back
        type(string)    :: fname
        integer(kind=8) :: sizes(4)
        integer         :: status, funit
        write(*,'(A)') 'test_gridding_manifest'
        fname = 'trail_chain_manifest_test.txt'
        sizes = [1024_8, 2048_8, 1025_8, 2049_8]
        call man%new(88, 2.9464285, 33016, 3, 2, 7, sizes, 10815.)
        call man%write(fname, status)
        call assert_int(0, status, 'manifest: write status')
        call back%read(fname, status)
        call assert_int(TRAIL_MANIFEST_OK, status, 'manifest: read status')
        call assert_int(88,    back%get_box(),     'manifest round trip: box')
        call assert_real(2.9464285, back%get_smpd(), 1.e-5, 'manifest round trip: sampling')
        call assert_int(33016, back%get_nptcls(),  'manifest round trip: row count')
        call assert_int(3,     back%get_nstates(), 'manifest round trip: nstates')
        call assert_int(2,     back%get_state(),   'manifest round trip: state')
        call assert_int(7,     back%get_gen(),     'manifest round trip: generation')
        call assert_true(back%get_size(3) == 1025_8, 'manifest round trip: component size')
        call assert_real(10815., back%get_mrep(), 1.e-3, 'manifest round trip: represented population M(s)')
        ! a manifest of another version (here the version-1 layout) is unreadable
        open(newunit=funit, file=fname%to_char(), status='replace', action='write')
        write(funit,*) 88, 2.9464285, 33016, 3, 2, 7, sizes
        close(funit)
        call back%read(fname, status)
        call assert_int(TRAIL_MANIFEST_UNREADABLE, status, 'manifest: another version is refused (re-seeded)')
        ! garbage and absence
        open(newunit=funit, file=fname%to_char(), status='replace', action='write')
        write(funit,'(A)') 'not a manifest'
        close(funit)
        call back%read(fname, status)
        call assert_int(TRAIL_MANIFEST_UNREADABLE, status, 'manifest: a corrupt manifest is refused')
        call del_file(fname)
        call back%read(fname, status)
        call assert_int(TRAIL_MANIFEST_MISSING, status, 'manifest: missing')
        call man%kill
        call back%kill
    end subroutine test_gridding_manifest

    ! The PCG chain's represented population is the header particle count of each half; the
    ! chain identity (provenance) carries a version, so a chain of an older version is not
    ! compatible and the master re-seeds it
    subroutine test_pcg_chain_header()
        integer,          parameter :: BOX  = 16
        real,             parameter :: SMPD = 2.0
        integer,          parameter :: MREP_HALF = 37
        character(len=*), parameter :: PROV_OLD = 'pcgtrail-v2|pgrp=c1|objfun=euclid'
        character(len=*), parameter :: PROV_NEW = 'pcgtrail-v3|pgrp=c1|objfun=euclid'
        type(reconstructor_pcg) :: op
        type(oris)              :: os
        type(string)            :: fname
        character(len=256)      :: prov
        real    :: smpd_file
        integer :: state, eo, part, nparts, nptcls, box_file, status
        write(*,'(A)') 'test_pcg_chain_header'
        fname = 'pcg_trail_chain_test.bin'
        call os%new(1, .false.)
        call op%new(BOX, SMPD, 1.e-2)
        call op%prep_particles(os, use_ctf=.false.)
        call op%begin_accum
        call op%write_raw_accum(fname, 1, 0, 1, 1, MREP_HALF, PROV_NEW)
        call read_pcg_raw_accum_header(fname, state, eo, part, nparts, nptcls, box_file, smpd_file, prov, status)
        call assert_int(0, status, 'PCG chain header: read status')
        call assert_int(MREP_HALF, nptcls, 'PCG chain header: represented population round trip')
        call assert_true(pcg_raw_accum_compatible(fname, BOX, SMPD, PROV_NEW), 'PCG chain: same identity is compatible')
        call op%write_raw_accum(fname, 1, 0, 1, 1, MREP_HALF, PROV_OLD)
        call assert_false(pcg_raw_accum_compatible(fname, BOX, SMPD, PROV_NEW), 'PCG chain: an older chain version is re-seeded')
        call del_file(fname)
        call op%kill
        call os%kill
    end subroutine test_pcg_chain_header

end module simple_trail_chain_manifest_tester
