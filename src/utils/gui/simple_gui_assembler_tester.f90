!@descr: unit tests for gui_assembler
module simple_gui_assembler_tester
use simple_gui_assembler
use simple_test_utils
use simple_string
implicit none

public :: run_all_gui_assembler_tests
private
#include "simple_local_flags.inc"

contains

  subroutine run_all_gui_assembler_tests()
    write(*,'(A)') '**** running all gui assembler tests ****'
    call test_new_kill()
    call test_new_kill_reuse()
    call test_preprocess()
  end subroutine run_all_gui_assembler_tests
    
  !---------------- lifecycle ----------------

  subroutine test_new_kill()
    type(gui_assembler) :: assembler
    write(*,'(A)') 'test_new_kill'
    call assembler%new()
    call assert_true(assembler%is_associated(), 'assembler json associated')
    call assembler%kill()
    call assert_true(.not.assembler%is_associated(), 'assembler json destroyed')
  end subroutine test_new_kill

  subroutine test_new_kill_reuse()
    type(gui_assembler) :: assembler
    write(*,'(A)') 'test_new_kill_reuse'
    call assembler%new()
    call assert_true(assembler%is_associated(), 'assembler json associated')
    call assembler%kill()
    call assert_true(.not.assembler%is_associated(), 'assembler json destroyed')
    call assembler%new()
    call assert_true(assembler%is_associated(), 'assembler json associated')
    call assembler%kill()
    call assert_true(.not.assembler%is_associated(), 'assembler json destroyed')
  end subroutine test_new_kill_reuse

  !---------------- preprocess assembly ----------------

  subroutine test_preprocess()
    type(gui_assembler)                              :: assembler
    type(gui_metadata_stream_preprocess)             :: meta_preprocess
    type(gui_metadata_micrograph),       allocatable :: meta_micrographs(:)
    type(gui_metadata_histogram),        allocatable :: meta_histograms(:)
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
    call meta_micrographs(1)%set(path=string('/test/path/to/micrograph1.mrc'), dfx=1.0, dfy=1.1, ctfres=0.1)
    call meta_micrographs(2)%set(path=string('/test/path/to/micrograph2.mrc'), dfx=2.0, dfy=1.2, ctfres=0.2)
    call meta_micrographs(3)%set(path=string('/test/path/to/micrograph3.mrc'), dfx=3.0, dfy=1.3, ctfres=0.3)
    call meta_micrographs(4)%set(path=string('/test/path/to/micrograph4.mrc'), dfx=4.0, dfy=1.4, ctfres=0.4)
    call meta_micrographs(5)%set(path=string('/test/path/to/micrograph5.mrc'), dfx=5.0, dfy=1.5, ctfres=0.5)
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
    
    call assembler%new()
    call assert_true(assembler%is_associated(), 'assembler json associated')
    call assembler%assemble_stream_preprocess(meta_preprocess, meta_micrographs, meta_histograms)
    json_str = assembler%to_string()
    call assert_true(json_str%strlen() > 0, 'json length greater that 0')
    json_hash = json_str%to_fnv1a_hash64()
    call assert_char(json_hash%to_char(), 'FF7D6F723627167D', 'assembler json hash matches')
    call assembler%kill()
    call assert_true(.not.assembler%is_associated(), 'assembler json destroyed')
    deallocate(meta_micrographs, meta_histograms, labels, data)
  end subroutine test_preprocess

end module simple_gui_assembler_tester