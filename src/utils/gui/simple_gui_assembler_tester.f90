!@descr: Unit tests for gui_assembler — lifecycle, hash suppression, and all assemble_stream_* procedures
!==============================================================================
! MODULE: simple_gui_assembler_tester
!
! PURPOSE:
!   Exercises gui_assembler through a set of unit tests covering object
!   lifecycle (new/kill/reuse/set_stoptime), hash-based change suppression
!   (clear_hashes), and JSON assembly for the stream-preprocess,
!   optics-assignment, initial-picking, reference-picking, and opening-2D stages.
!   Where the assembled JSON is fully deterministic (no live timestamps) the
!   test verifies an FNV-1a hash of the serialised output; otherwise it checks
!   only that the output is non-empty.
!   Note: assemble_stream_heartbeat is not tested here because it requires
!   live forked child processes.
!
! ENTRY POINT:
!   run_all_gui_assembler_tests() — run every test in this module
!
! DEPENDENCIES:
!   simple_gui_metadata_api, simple_gui_assembler, simple_test_utils,
!   simple_string
!==============================================================================
module simple_gui_assembler_tester
  use simple_gui_metadata_api, only: gui_metadata_stream_preprocess,             &
                                     gui_metadata_micrograph,                    &
                                     gui_metadata_histogram,                     &
                                     gui_metadata_timeplot,                      &
                                     gui_metadata_stream_optics_assignment,      &
                                     gui_metadata_optics_group,                  &
                                     gui_metadata_stream_picking,        &
                                     gui_metadata_stream_opening2D,              &
                                     GUI_METADATA_STREAM_PREPROCESS_TYPE,        &
                                     GUI_METADATA_MICROGRAPH_TYPE,               &
                                     GUI_METADATA_HISTOGRAM_TYPE,                &
                                     GUI_METADATA_STREAM_OPTICS_ASSIGNMENT_TYPE, &
                                     GUI_METADATA_OPTICS_GROUP_TYPE,             &
                                     GUI_METADATA_STREAM_INITIAL_PICKING_TYPE,   &
                                     GUI_METADATA_STREAM_REFERENCE_PICKING_TYPE,      &
                                     GUI_METADATA_STREAM_REFERENCE_PICKING_CLS2D_TYPE, &
                                     GUI_METADATA_STREAM_OPENING2D_TYPE,              &
                                     GUI_METADATA_STREAM_OPENING2D_CLS2D_TYPE,        &
                                     gui_metadata_cavg2D,                        &
                                     sprite_sheet_pos
  use simple_gui_assembler,    only: gui_assembler
  use simple_test_utils,       only: assert_true, assert_char
  use simple_string,           only: string
  implicit none

public :: run_all_gui_assembler_tests
private
#include "simple_local_flags.inc"

contains

  ! Run the complete gui_assembler test suite.
  subroutine run_all_gui_assembler_tests()
    write(*,'(A)') '**** running all gui assembler tests ****'
    call test_new_kill()
    call test_new_kill_reuse()
    call test_set_stoptime()
    call test_clear_hashes()
    call test_preprocess()
    call test_optics_assignment()
    call test_initial_picking()
    call test_reference_picking()
    call test_opening2D()
  end subroutine run_all_gui_assembler_tests

  !---------------- lifecycle ----------------

  ! Verify that a single new/kill cycle leaves the assembler in a clean state.
  subroutine test_new_kill()
    type(gui_assembler) :: assembler
    write(*,'(A)') 'test_new_kill'
    call assembler%new(0)
    call assert_true(assembler%is_associated(), 'assembler json associated')
    call assembler%kill()
    call assert_true(.not.assembler%is_associated(), 'assembler json destroyed')
  end subroutine test_new_kill

  ! Verify that the assembler can be re-initialised after being killed.
  subroutine test_new_kill_reuse()
    type(gui_assembler) :: assembler
    write(*,'(A)') 'test_new_kill_reuse'
    call assembler%new(0)
    call assert_true(assembler%is_associated(), 'assembler json associated')
    call assembler%kill()
    call assert_true(.not.assembler%is_associated(), 'assembler json destroyed')
    call assembler%new(1)
    call assert_true(assembler%is_associated(), 'assembler json associated')
    call assembler%kill()
    call assert_true(.not.assembler%is_associated(), 'assembler json destroyed')
  end subroutine test_new_kill_reuse

  ! Verify that set_stoptime can be called without error and the assembler
  ! remains initialised.
  subroutine test_set_stoptime()
    type(gui_assembler) :: assembler
    write(*,'(A)') 'test_set_stoptime'
    call assembler%new(0)
    call assembler%set_stoptime()
    call assert_true(assembler%is_associated(), 'assembler still associated after set_stoptime')
    call assembler%kill()
    call assert_true(.not.assembler%is_associated(), 'assembler json destroyed')
  end subroutine test_set_stoptime

  !---------------- clear_hashes ----------------

  ! Verify that clear_hashes resets content-hash suppression so that an unchanged
  ! section is retransmitted on the next assembly call.
  subroutine test_clear_hashes()
    type(gui_assembler)                              :: assembler
    type(gui_metadata_stream_preprocess)             :: meta_preprocess
    type(gui_metadata_micrograph),       allocatable :: meta_micrographs(:)
    type(gui_metadata_histogram),        allocatable :: meta_histograms(:)
    type(gui_metadata_timeplot),         allocatable :: meta_timeplots(:)
    type(string)                                     :: json_str1, json_str2, json_str3
    write(*,'(A)') 'test_clear_hashes'
    call meta_preprocess%new(GUI_METADATA_STREAM_PREPROCESS_TYPE)
    call meta_preprocess%set(stage=string('test stage'), movies_imported=100, movies_processed=50, &
                             movies_rejected=10, movies_rate=5, average_ctf_res=4.0,               &
                             average_ice_score=0.1, average_astigmatism=0.5,                       &
                             cutoff_ctf_res=8.0, cutoff_ice_score=0.2, cutoff_astigmatism=0.5)
    call assembler%new(0)
    ! first assembly — section included
    call assembler%assemble_stream_preprocess(meta_preprocess, meta_micrographs, meta_histograms, meta_timeplots)
    json_str1 = assembler%to_string()
    ! second identical assembly — section suppressed (hash unchanged)
    call assembler%assemble_stream_preprocess(meta_preprocess, meta_micrographs, meta_histograms, meta_timeplots)
    json_str2 = assembler%to_string()
    call assert_true(json_str2%strlen() < json_str1%strlen(), 'second identical assembly shorter (suppressed)')
    ! after clear_hashes the section is retransmitted
    call assembler%clear_hashes()
    call assembler%assemble_stream_preprocess(meta_preprocess, meta_micrographs, meta_histograms, meta_timeplots)
    json_str3 = assembler%to_string()
    call assert_true(json_str3%strlen() == json_str1%strlen(), 'assembly after clear_hashes same length as first')
    call assembler%kill()
    call assert_true(.not.assembler%is_associated(), 'assembler json destroyed')
  end subroutine test_clear_hashes

  !---------------- preprocess assembly ----------------

  ! Assemble a stream-preprocess JSON payload from synthetic metadata and verify
  ! its structure is stable via an FNV-1a hash of the serialised output.
  subroutine test_preprocess()
    type(gui_assembler)                              :: assembler
    type(gui_metadata_stream_preprocess)             :: meta_preprocess
    type(gui_metadata_micrograph),       allocatable :: meta_micrographs(:)
    type(gui_metadata_histogram),        allocatable :: meta_histograms(:)
    type(gui_metadata_timeplot),         allocatable :: meta_timeplots(:)
    type(string)                                     :: json_str, json_hash
    real,                                allocatable :: labels(:)
    integer,                             allocatable :: data(:)
    integer                                          :: i
    write(*,'(A)') 'test_preprocess'
    call meta_preprocess%new(GUI_METADATA_STREAM_PREPROCESS_TYPE)
    call meta_preprocess%set(stage=string('test stage'), movies_imported=100000, movies_processed=50000, movies_rejected=1000,&
                             movies_rate=400, average_ctf_res=4.5, average_ice_score=0.15, average_astigmatism=0.75, cutoff_ctf_res=8.0,&
                             cutoff_ice_score=0.2, cutoff_astigmatism=0.5)
    ! micrographs
    allocate(meta_micrographs(5))
    do i=1, size(meta_micrographs)
      call meta_micrographs(i)%new(GUI_METADATA_MICROGRAPH_TYPE)
    enddo
    call meta_micrographs(1)%set(path=string('/test/path/to/micrograph1.mrc'), dfx=1.0, dfy=1.1, ctfres=0.1, i_max=5, i=1)
    call meta_micrographs(2)%set(path=string('/test/path/to/micrograph2.mrc'), dfx=2.0, dfy=1.2, ctfres=0.2, i_max=5, i=2)
    call meta_micrographs(3)%set(path=string('/test/path/to/micrograph3.mrc'), dfx=3.0, dfy=1.3, ctfres=0.3, i_max=5, i=3)
    call meta_micrographs(4)%set(path=string('/test/path/to/micrograph4.mrc'), dfx=4.0, dfy=1.4, ctfres=0.4, i_max=5, i=4)
    call meta_micrographs(5)%set(path=string('/test/path/to/micrograph5.mrc'), dfx=5.0, dfy=1.5, ctfres=0.5, i_max=5, i=5)
    ! histograms
    allocate(meta_histograms(5), labels(10), data(10))
    do i=1, 10
      labels(i) = 1.0/i
      data(i)   = i
    enddo
    do i=1, size(meta_histograms)
      call meta_histograms(i)%new(GUI_METADATA_HISTOGRAM_TYPE)
    enddo
    call meta_histograms(1)%set(name=string('histogram1'), labels=labels, data=data)
    call meta_histograms(2)%set(name=string('histogram2'), labels=labels, data=data)
    call meta_histograms(3)%set(name=string('histogram3'), labels=labels, data=data)
    call meta_histograms(4)%set(name=string('histogram4'), labels=labels, data=data)
    call meta_histograms(5)%set(name=string('histogram5'), labels=labels, data=data)
    ! timeplots
    call assembler%new(0)
    call assert_true(assembler%is_associated(), 'assembler json associated')
    call assembler%assemble_stream_preprocess(meta_preprocess, meta_micrographs, meta_histograms, meta_timeplots)
    json_str = assembler%to_string()
    call assert_true(json_str%strlen() > 0, 'json length greater that 0')
    json_hash = json_str%to_fnv1a_hash64()
    call assert_char(json_hash%to_char(), 'DE931F774EBE58C2', 'assembler json hash matches')
    call assembler%kill()
    call assert_true(.not.assembler%is_associated(), 'assembler json destroyed')
    deallocate(meta_micrographs, meta_histograms, labels, data)
  end subroutine test_preprocess

  !---------------- optics assignment assembly ----------------

  ! Assemble a stream optics-assignment JSON payload from synthetic metadata and
  ! verify its structure is stable via an FNV-1a hash of the serialised output.
  ! Assemble a stream optics-assignment JSON payload from synthetic metadata and
  ! verify its structure is stable via an FNV-1a hash of the serialised output.
  ! NOTE: hash recalibration needed — optics groups are now included in the hash
  ! (previously they were added after hash computation and omitted from it).
  subroutine test_optics_assignment()
    type(gui_assembler)                              :: assembler
    type(gui_metadata_stream_optics_assignment)      :: meta_optics_assignment
    type(gui_metadata_optics_group),     allocatable :: meta_optics_groups(:)
    type(string)                                     :: json_str, json_hash
    real,                                allocatable :: xshifts(:), yshifts(:)
    integer                                          :: i
    write(*,'(A)') 'test_optics_assignment'
    allocate(xshifts(10), yshifts(10))
    do i=1, 10
      xshifts(i) = 1.0/i
      yshifts(i) = i
    enddo
    call meta_optics_assignment%new(GUI_METADATA_STREAM_OPTICS_ASSIGNMENT_TYPE)
    call meta_optics_assignment%set(stage=string('test stage'), micrographs_assigned=100, optics_groups_assigned=5, last_micrograph_imported=1234)
    allocate(meta_optics_groups(5))
    do i=1, size(meta_optics_groups)
      call meta_optics_groups(i)%new(GUI_METADATA_OPTICS_GROUP_TYPE)
      call meta_optics_groups(i)%set(i=i, i_max=size(meta_optics_groups), xshifts=xshifts, yshifts=yshifts, n_shifts=size(xshifts))
    enddo
    call assembler%new(0)
    call assert_true(assembler%is_associated(), 'assembler json associated')
    call assembler%assemble_stream_optics_assignment(meta_optics_assignment, meta_optics_groups)
    json_str = assembler%to_string()
    call assert_true(json_str%strlen() > 0, 'json length greater than 0')
    json_hash = json_str%to_fnv1a_hash64()
    call assert_char(json_hash%to_char(), 'BDBE34825830B633', 'assembler json hash matches')
    call assembler%kill()
    call assert_true(.not.assembler%is_associated(), 'assembler json destroyed')
    deallocate(xshifts, yshifts)
  end subroutine test_optics_assignment

  !---------------- initial picking assembly ----------------

  ! Assemble a stream initial-picking JSON payload from synthetic metadata and
  ! verify the section is non-empty.  An exact hash comparison is not possible
  ! because the section embeds a live Unix timestamp (last_micrograph_imported).
  subroutine test_initial_picking()
    type(gui_assembler)                              :: assembler
    type(gui_metadata_stream_picking)        :: meta_initial_picking
    type(gui_metadata_micrograph),       allocatable :: meta_micrographs(:)
    type(string)                                     :: json_str
    integer                                          :: i
    write(*,'(A)') 'test_initial_picking'
    call meta_initial_picking%new(GUI_METADATA_STREAM_INITIAL_PICKING_TYPE)
    call meta_initial_picking%set(stage=string('test stage'), micrographs_imported=200, &
                                  micrographs_accepted=180, particles_extracted=36000, box_size=256)
    allocate(meta_micrographs(3))
    do i=1, size(meta_micrographs)
      call meta_micrographs(i)%new(GUI_METADATA_MICROGRAPH_TYPE)
    enddo
    call meta_micrographs(1)%set(path=string('/test/path/mic1.mrc'), dfx=1.0, dfy=1.1, ctfres=3.5, i_max=3, i=1)
    call meta_micrographs(2)%set(path=string('/test/path/mic2.mrc'), dfx=2.0, dfy=1.2, ctfres=4.0, i_max=3, i=2)
    call meta_micrographs(3)%set(path=string('/test/path/mic3.mrc'), dfx=3.0, dfy=1.3, ctfres=4.5, i_max=3, i=3)
    call assembler%new(0)
    call assert_true(assembler%is_associated(), 'assembler json associated')
    call assembler%assemble_stream_initial_picking(meta_initial_picking, meta_micrographs)
    json_str = assembler%to_string()
    call assert_true(json_str%strlen() > 0, 'json length greater than 0')
    call assembler%kill()
    call assert_true(.not.assembler%is_associated(), 'assembler json destroyed')
    deallocate(meta_micrographs)
  end subroutine test_initial_picking

  !---------------- reference picking assembly ----------------

  ! Assemble a stream reference-picking JSON payload from synthetic metadata and
  ! verify the section is non-empty.  An exact hash comparison is not possible
  ! because the section embeds a live Unix timestamp (last_micrograph_imported).
  subroutine test_reference_picking()
    type(gui_assembler)                              :: assembler
    type(gui_metadata_stream_picking)                :: meta_reference_picking
    type(gui_metadata_micrograph),       allocatable :: meta_micrographs(:)
    type(gui_metadata_cavg2D),           allocatable :: meta_cavgs2D(:)
    type(string)                                     :: json_str
    integer                                          :: i
    write(*,'(A)') 'test_reference_picking'
    call meta_reference_picking%new(GUI_METADATA_STREAM_REFERENCE_PICKING_TYPE)
    call meta_reference_picking%set(stage=string('test stage'), micrographs_imported=200, &
                                    micrographs_accepted=180, particles_extracted=36000, box_size=256)
    allocate(meta_micrographs(3))
    do i=1, size(meta_micrographs)
      call meta_micrographs(i)%new(GUI_METADATA_MICROGRAPH_TYPE)
    enddo
    call meta_micrographs(1)%set(path=string('/test/path/mic1.mrc'), dfx=1.0, dfy=1.1, ctfres=3.5, i_max=3, i=1)
    call meta_micrographs(2)%set(path=string('/test/path/mic2.mrc'), dfx=2.0, dfy=1.2, ctfres=4.0, i_max=3, i=2)
    call meta_micrographs(3)%set(path=string('/test/path/mic3.mrc'), dfx=3.0, dfy=1.3, ctfres=4.5, i_max=3, i=3)
    allocate(meta_cavgs2D(3))
    do i=1, size(meta_cavgs2D)
      call meta_cavgs2D(i)%new(GUI_METADATA_STREAM_REFERENCE_PICKING_CLS2D_TYPE)
    enddo
    call meta_cavgs2D(1)%set(path=string('/test/path/ref.jpeg'), mrcpath=string('/test/path/ref.mrc'), &
                             idx=1, sprite=sprite_sheet_pos(x=0.0,  y=0.0, h=256, w=768), i=1, i_max=3)
    call meta_cavgs2D(2)%set(path=string('/test/path/ref.jpeg'), mrcpath=string('/test/path/ref.mrc'), &
                             idx=2, sprite=sprite_sheet_pos(x=33.3, y=0.0, h=256, w=768), i=2, i_max=3)
    call meta_cavgs2D(3)%set(path=string('/test/path/ref.jpeg'), mrcpath=string('/test/path/ref.mrc'), &
                             idx=3, sprite=sprite_sheet_pos(x=66.6, y=0.0, h=256, w=768), i=3, i_max=3)
    call assembler%new(0)
    call assert_true(assembler%is_associated(), 'assembler json associated')
    call assembler%assemble_stream_reference_picking(meta_reference_picking, meta_micrographs, meta_cavgs2D)
    json_str = assembler%to_string()
    call assert_true(json_str%strlen() > 0, 'json length greater than 0')
    call assembler%kill()
    call assert_true(.not.assembler%is_associated(), 'assembler json destroyed')
    deallocate(meta_micrographs, meta_cavgs2D)
  end subroutine test_reference_picking

  !---------------- opening2D assembly ----------------

  ! Assemble a stream opening-2D JSON payload from synthetic metadata and verify
  ! the section is non-empty.  An exact hash comparison is not possible because
  ! the section embeds a live Unix timestamp (last_particles_imported).
  subroutine test_opening2D()
    type(gui_assembler)                            :: assembler
    type(gui_metadata_stream_opening2D)            :: meta_opening2D
    type(gui_metadata_cavg2D),         allocatable :: meta_cavgs2D(:), meta_final_cavgs2D(:)
    type(string)                                   :: json_str
    integer                                        :: i
    write(*,'(A)') 'test_opening2D'
    call meta_opening2D%new(GUI_METADATA_STREAM_OPENING2D_TYPE)
    call meta_opening2D%set(stage=string('test stage'), particles_imported=50000, &
                            particles_accepted=42000, mask_diam=160, box_size=256, mask_scale=0.75)
    allocate(meta_cavgs2D(3))
    do i=1, size(meta_cavgs2D)
      call meta_cavgs2D(i)%new(GUI_METADATA_STREAM_OPENING2D_CLS2D_TYPE)
    enddo
    call meta_cavgs2D(1)%set(path=string('/test/path/cls.jpeg'), mrcpath=string('/test/path/cls.mrc'), &
                             idx=1, sprite=sprite_sheet_pos(x=0.0,  y=0.0, h=256, w=768), i=1, i_max=3)
    call meta_cavgs2D(2)%set(path=string('/test/path/cls.jpeg'), mrcpath=string('/test/path/cls.mrc'), &
                             idx=2, sprite=sprite_sheet_pos(x=33.3, y=0.0, h=256, w=768), i=2, i_max=3)
    call meta_cavgs2D(3)%set(path=string('/test/path/cls.jpeg'), mrcpath=string('/test/path/cls.mrc'), &
                             idx=3, sprite=sprite_sheet_pos(x=66.6, y=0.0, h=256, w=768), i=3, i_max=3)
    call assembler%new(0)
    call assert_true(assembler%is_associated(), 'assembler json associated')
    call assembler%assemble_stream_opening2D(meta_opening2D, meta_cavgs2D, meta_final_cavgs2D)
    json_str = assembler%to_string()
    call assert_true(json_str%strlen() > 0, 'json length greater than 0')
    call assembler%kill()
    call assert_true(.not.assembler%is_associated(), 'assembler json destroyed')
    deallocate(meta_cavgs2D)
  end subroutine test_opening2D

end module simple_gui_assembler_tester