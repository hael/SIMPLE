!@descr: unit tests for gui metadata
module simple_gui_metadata_tester
use simple_test_utils
use simple_gui_metadata_api
implicit none

public :: run_all_gui_metadata_tests
private
#include "simple_local_flags.inc"

contains

  subroutine run_all_gui_metadata_tests()
    write(*,'(A)') '**** running all gui metadata tests ****'
    call test_new_kill()
    call test_serialise()
    call test_set_get_micrograph()
    call test_serialise_micrograph()
    call test_jsonise_micrograph()
    call test_set_get_histogram()
    call test_serialise_histogram()
    call test_jsonise_histogram()
    call test_set_get_stream_preprocess()
    call test_serialise_stream_preprocess()
    call test_jsonise_stream_preprocess()
  end subroutine run_all_gui_metadata_tests
    
  !---------------- base ----------------

  subroutine test_new_kill()
    type(gui_metadata_base) :: meta
    write(*,'(A)') 'test_new_kill'
    call meta%new(GUI_METADATA_MICROGRAPH_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_MICROGRAPH_TYPE, 'type is set correctly')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
  end subroutine test_new_kill

  subroutine test_serialise()
    character(len=:), allocatable :: buffer
    type(gui_metadata_base)       :: meta
    write(*,'(A)') 'test_serialise'
    call meta%new(GUI_METADATA_MICROGRAPH_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_MICROGRAPH_TYPE, 'type is set correctly')
    call meta%serialise(buffer=buffer)
    call assert_true(allocated(buffer), 'buffer allocated')
    call assert_int(len(buffer), int(sizeof(meta), kind=4), 'empty buffer')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
  end subroutine test_serialise

  !---------------- micrograph ----------------

  subroutine test_set_get_micrograph()
    character(len=:),             allocatable :: buffer
    type(gui_metadata_micrograph)             :: meta
    type(string)                              :: path
    real                                      :: dfx, dfy, ctfres
    write(*,'(A)') 'test_set_get_micrograph'
    call meta%new(GUI_METADATA_MICROGRAPH_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_MICROGRAPH_TYPE, 'type is set correctly')
    call meta%set(path=string('/test/path/to/micrograph.mrc'), dfx=0.25, dfy=2.89, ctfres=8.9)
    call assert_true(meta%assigned(), 'metadata object is set')
    call assert_true(meta%get(path=path, dfx=dfx, dfy=dfy, ctfres=ctfres), 'metadata retreived')
    call assert_char(path%to_char(), '/test/path/to/micrograph.mrc', 'path set/get correctly')
    call assert_true(dfx == 0.25,   'dfx set/get correctly' )
    call assert_true(dfy == 2.89,   'dfy set/get correctly')
    call assert_true(ctfres == 8.9, 'ctfres set/get correctly' )
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
  end subroutine test_set_get_micrograph

  subroutine test_serialise_micrograph()
    character(len=:),             allocatable :: buffer
    type(gui_metadata_micrograph)             :: meta
    write(*,'(A)') 'test_serialise_micrograph'
    call meta%new(GUI_METADATA_MICROGRAPH_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_MICROGRAPH_TYPE, 'type is set correctly')
    call meta%set(path=string('/test/path/to/micrograph.mrc'), dfx=0.25, dfy=2.89, ctfres=8.9)
    call assert_true(meta%assigned(), 'metadata object is set')
    call meta%serialise(buffer=buffer)
    call assert_true(allocated(buffer), 'buffer allocated')
    call assert_int(len(buffer), int(sizeof(meta), kind=4), 'buffer correct size')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
  end subroutine test_serialise_micrograph

  subroutine test_jsonise_micrograph()
    character(kind=CK, len=:),    allocatable :: buffer
    type(gui_metadata_micrograph)             :: meta
    type(json_core)                           :: json
    type(json_value),             pointer     :: json_ptr
    type(string)                              :: json_str, json_hash
    write(*,'(A)') 'test_jsonise_micrograph'
    call json%initialize(no_whitespace=.true., compact_reals=.true.)
    call meta%new(GUI_METADATA_MICROGRAPH_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_MICROGRAPH_TYPE, 'type is set correctly')
    call meta%set(path=string('/test/path/to/micrograph.mrc'), dfx=0.25, dfy=2.89, ctfres=8.9)
    call assert_true(meta%assigned(), 'metadata object is set')
    json_ptr => meta%jsonise()
    call json%print_to_string(json_ptr, buffer)
    call assert_true(allocated(buffer), 'buffer allocated')
    call assert_int(len(buffer), 114, 'buffer correct size')
    json_str = buffer
    call assert_int(json_str%strlen(), 114, 'string correct size')
    json_hash = json_str%to_fnv1a_hash64()
    call assert_char(json_hash%to_char(), '1D2D2E9B461DD2EC', 'correct checksum')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
    call json%destroy(json_ptr)
    call assert_true(.not.json%failed(), 'json destroyed')
  end subroutine test_jsonise_micrograph

  !---------------- histogram ----------------

  subroutine test_set_get_histogram()
    type(gui_metadata_histogram)              :: meta
    type(string)                              :: name
    real,                         allocatable :: labels(:)
    integer,                      allocatable :: data(:)
    integer                                   :: i
    write(*,'(A)') 'test_set_get_histogram'
    allocate(labels(10), data(10))
    do i=1, 10
      labels(i) = 1.0/i
      data(i)   = i
    enddo
    call meta%new(GUI_METADATA_HISTOGRAM_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_HISTOGRAM_TYPE, 'type is set correctly')
    call meta%set(name=string('histname'), labels=labels, data=data)
    call assert_true(meta%assigned(), 'metadata object is set')
    deallocate(labels, data)
    call assert_true(meta%get(name=name, labels=labels, data=data), 'metadata retreived')
    call assert_char(name%to_char(), 'histname', 'name set/get correctly')
    call assert_true(size(labels) == 10,   'labels correctly allocated')
    call assert_true(size(data)   == 10,   'data correctly allocated')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
    deallocate(labels, data)
  end subroutine test_set_get_histogram

  subroutine test_serialise_histogram()
    character(len=:),             allocatable :: buffer
    type(gui_metadata_histogram)              :: meta
    type(string)                              :: name
    real,                         allocatable :: labels(:)
    integer,                      allocatable :: data(:)
    integer                                   :: i
    write(*,'(A)') 'test_serialise_histogram'
    allocate(labels(10), data(10))
    do i=1, 10
      labels(i) = 1.0/i
      data(i)   = i
    enddo
    call meta%new(GUI_METADATA_HISTOGRAM_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_HISTOGRAM_TYPE, 'type is set correctly')
    call meta%set(name=string('histname'), labels=labels, data=data)
    call assert_true(meta%assigned(), 'metadata object is set')
    call meta%serialise(buffer=buffer)
    call assert_true(allocated(buffer), 'buffer allocated')
    call assert_int(len(buffer), int(sizeof(meta), kind=4), 'buffer correct size')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
    deallocate(labels, data)
  end subroutine test_serialise_histogram

  subroutine test_jsonise_histogram()
    character(kind=CK, len=:),    allocatable :: buffer
    type(gui_metadata_histogram)              :: meta
    type(json_core)                           :: json
    type(json_value),             pointer     :: json_ptr
    type(string)                              :: json_str, json_hash
    real,                         allocatable :: labels(:)
    integer,                      allocatable :: data(:)
    integer                                   :: i
    write(*,'(A)') 'test_jsonise_histogram'
    allocate(labels(10), data(10))
    do i=1, 10
      labels(i) = 1.0/i
      data(i)   = i
    enddo
    call meta%new(GUI_METADATA_HISTOGRAM_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_HISTOGRAM_TYPE, 'type is set correctly')
    call meta%set(name=string('histname'), labels=labels, data=data)
    call assert_true(meta%assigned(), 'metadata object is set')
    json_ptr => meta%jsonise()
    call json%print_to_string(json_ptr, buffer)
    call assert_true(allocated(buffer), 'buffer allocated')
    call assert_int(len(buffer), 394, 'buffer correct size')
    json_str = buffer
    call assert_int(json_str%strlen(), 394, 'string correct size')
    json_hash = json_str%to_fnv1a_hash64()
    call assert_char(json_hash%to_char(), 'B99ABC60B24FDFBA', 'correct checksum')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
    call json%destroy(json_ptr)
    call assert_true(.not.json%failed(), 'json destroyed')
    deallocate(labels, data)
  end subroutine test_jsonise_histogram

  !---------------- stream preprocess ----------------

  subroutine test_set_get_stream_preprocess()
    type(gui_metadata_stream_preprocess)             :: meta
    type(string)                                     :: stage
    integer                                          :: movies_imported, movies_processed, movies_rejected, movies_rate
    real                                             :: average_ctf_res, average_ice_score, average_astigmatism
    real                                             :: cutoff_ctf_res, cutoff_ice_score, cutoff_astigmatism
    write(*,'(A)') 'test_set_get_stream_preprocess'
    call meta%new(GUI_METADATA_STREAM_PREPROCESS_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_STREAM_PREPROCESS_TYPE, 'type is set correctly')
    call meta%set(stage=string('test stage'), movies_imported=100000, movies_processed=50000, movies_rejected=1000,&
                  movies_rate=400, average_ctf_res=4.5, average_ice_score=0.15, average_astigmatism=0.75, cutoff_ctf_res=8.0,&
                  cutoff_ice_score=0.2, cutoff_astigmatism=0.5)
    call assert_true(meta%assigned(), 'metadata object is set')
    call assert_true(meta%get(stage=stage, movies_imported=movies_imported, movies_processed=movies_processed, movies_rejected=movies_rejected,&
                              movies_rate=movies_rate, average_ctf_res=average_ctf_res, average_ice_score=average_ice_score,&
                              average_astigmatism=average_astigmatism, cutoff_ctf_res=cutoff_ctf_res, cutoff_ice_score=cutoff_ice_score,&
                              cutoff_astigmatism=cutoff_astigmatism), 'metadata retreived')
    call assert_char(stage%to_char(), 'test stage', 'stage set/get correctly')
    call assert_int(movies_imported,  100000, 'movies_imported set/get correctly' )
    call assert_int(movies_processed, 50000,  'movies_processed set/get correctly')
    call assert_int(movies_rejected,  1000,   'movies_rejected set/get correctly' )
    call assert_int(movies_rate,      400,    'movies_rate set/get correctly'     )
    call assert_true(average_ctf_res     == 4.5,  'average_ctf_res set/get correctly'    )
    call assert_true(average_ice_score   == 0.15, 'average_ice_score set/get correctly'  )
    call assert_true(average_astigmatism == 0.75, 'average_astigmatism set/get correctly')
    call assert_true(cutoff_ctf_res      == 8.0,  'cutoff_ctf_res set/get correctly'     )
    call assert_true(cutoff_ice_score    == 0.2,  'cutoff_ice_score set/get correctly'   )
    call assert_true(cutoff_astigmatism  == 0.5,  'cutoff_astigmatism set/get correctly' )
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
  end subroutine test_set_get_stream_preprocess

  subroutine test_serialise_stream_preprocess()
    character(len=:),                    allocatable :: buffer
    type(gui_metadata_stream_preprocess)             :: meta
    write(*,'(A)') 'test_serialise_stream_preprocess'
    call meta%new(GUI_METADATA_STREAM_PREPROCESS_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_STREAM_PREPROCESS_TYPE, 'type is set correctly')
    call meta%set(stage=string('test stage'), movies_imported=100000, movies_processed=50000, movies_rejected=1000,&
                  movies_rate=400, average_ctf_res=4.5, average_ice_score=0.15, average_astigmatism=0.75, cutoff_ctf_res=8.0,&
                  cutoff_ice_score=0.2, cutoff_astigmatism=0.5)
    call assert_true(meta%assigned(), 'metadata object is set')
    call meta%serialise(buffer=buffer)
    call assert_true(allocated(buffer), 'buffer allocated')
    call assert_int(len(buffer), int(sizeof(meta), kind=4), 'buffer correct size')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
    deallocate(buffer)
  end subroutine test_serialise_stream_preprocess

  subroutine test_jsonise_stream_preprocess()
    character(kind=CK, len=:),    allocatable :: buffer
    type(gui_metadata_stream_preprocess)      :: meta
    type(json_core)                           :: json
    type(json_value),             pointer     :: json_ptr
    type(string)                              :: json_str, json_hash
    write(*,'(A)') 'test_jsonise_stream_preprocess'
    call json%initialize(no_whitespace=.true., compact_reals=.true.)
    call meta%new(GUI_METADATA_STREAM_PREPROCESS_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_STREAM_PREPROCESS_TYPE, 'type is set correctly')
    call meta%set(stage=string('test stage'), movies_imported=100000, movies_processed=50000, movies_rejected=1000,&
                  movies_rate=400, average_ctf_res=4.5, average_ice_score=0.15, average_astigmatism=0.75, cutoff_ctf_res=8.0,&
                  cutoff_ice_score=0.2, cutoff_astigmatism=0.5)
    call assert_true(meta%assigned(), 'metadata object is set')
    json_ptr => meta%jsonise()
    call json%print_to_string(json_ptr, buffer)
    call assert_true(allocated(buffer), 'buffer allocated')
    call assert_int(len(buffer), 306, 'buffer correct size')
    json_str = buffer
    call assert_int(json_str%strlen(), 306, 'string correct size')
    json_hash = json_str%to_fnv1a_hash64()
    call assert_char(json_hash%to_char(), 'B57B49274E327BAE', 'correct checksum')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
    call json%destroy(json_ptr)
    call assert_true(.not.json%failed(), 'json destroyed')
    deallocate(buffer)
  end subroutine test_jsonise_stream_preprocess


end module simple_gui_metadata_tester