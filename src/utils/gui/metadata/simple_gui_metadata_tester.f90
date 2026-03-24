!@descr: Unit tests for all gui_metadata types — lifecycle, serialisation, and JSON serialisation
!==============================================================================
! MODULE: simple_gui_metadata_tester
!
! PURPOSE:
!   Exercises every public gui_metadata type through three test tiers:
!     set/get   — verify that all fields round-trip correctly
!     serialise — verify the binary transfer buffer is the expected size
!     jsonise   — verify the JSON output via FNV-1a hash (for deterministic
!                 types) or a non-empty length check (for types that embed
!                 a live Unix timestamp)
!   Types covered:
!     gui_metadata_base, gui_metadata_micrograph, gui_metadata_histogram,
!     gui_metadata_timeplot, gui_metadata_optics_group,
!     gui_metadata_stream_preprocess, gui_metadata_stream_optics_assignment,
!     gui_metadata_stream_update, gui_metadata_stream_initial_picking,
!     gui_metadata_stream_opening2D
!
! ENTRY POINT:
!   run_all_gui_metadata_tests() — run every test in this module
!
! DEPENDENCIES:
!   simple_test_utils, simple_gui_metadata_api
!==============================================================================
module simple_gui_metadata_tester
  use simple_test_utils,    only: assert_true, assert_int, assert_char
  use simple_gui_metadata_api
  implicit none

public :: run_all_gui_metadata_tests
private
#include "simple_local_flags.inc"

contains

  ! Run the complete gui_metadata test suite.
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
    call test_set_get_timeplot()
    call test_serialise_timeplot()
    call test_jsonise_timeplot()
    call test_set_get_optics_group()
    call test_serialise_optics_group()
    call test_jsonise_optics_group()
    call test_set_get_stream_preprocess()
    call test_serialise_stream_preprocess()
    call test_jsonise_stream_preprocess()
    call test_set_get_stream_optics_assignment()
    call test_serialise_stream_optics_assignment()
    call test_jsonise_stream_optics_assignment()
    call test_set_get_stream_update()
    call test_serialise_stream_update()
    call test_jsonise_stream_update()
    call test_set_get_stream_initial_picking()
    call test_serialise_stream_initial_picking()
    call test_jsonise_stream_initial_picking()
    call test_set_get_stream_opening2D()
    call test_serialise_stream_opening2D()
    call test_jsonise_stream_opening2D()
  end subroutine run_all_gui_metadata_tests

  !---------------- base ----------------

  ! Verify that new() initialises the type field and kill() clears it.
  subroutine test_new_kill()
    type(gui_metadata_base) :: meta
    write(*,'(A)') 'test_new_kill'
    call meta%new(GUI_METADATA_MICROGRAPH_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_MICROGRAPH_TYPE, 'type is set correctly')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
  end subroutine test_new_kill

  ! Verify that serialise() produces a buffer of sizeof(meta) bytes.
  subroutine test_serialise()
    character(len=:), allocatable :: buffer
    type(gui_metadata_base)       :: meta
    write(*,'(A)') 'test_serialise'
    call meta%new(GUI_METADATA_MICROGRAPH_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_MICROGRAPH_TYPE, 'type is set correctly')
    call meta%serialise(buffer=buffer)
    call assert_true(allocated(buffer), 'buffer allocated')
    call assert_int(len(buffer), int(sizeof(meta), kind=4), 'buffer correct size')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
  end subroutine test_serialise

  !---------------- micrograph ----------------

  ! Verify that all micrograph fields round-trip through set/get.
  subroutine test_set_get_micrograph()
    type(gui_metadata_micrograph)             :: meta
    type(string)                              :: path
    integer                                   :: i_max, i
    real                                      :: dfx, dfy, ctfres
    write(*,'(A)') 'test_set_get_micrograph'
    call meta%new(GUI_METADATA_MICROGRAPH_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_MICROGRAPH_TYPE, 'type is set correctly')
    call meta%set(path=string('/test/path/to/micrograph.mrc'), dfx=0.25, dfy=2.89, ctfres=8.9, i_max=1, i=1)
    call assert_true(meta%assigned(), 'metadata object is set')
    call assert_true(meta%get(path=path, dfx=dfx, dfy=dfy, ctfres=ctfres, i_max=i_max, i=i), 'metadata retrieved')
    call assert_char(path%to_char(), '/test/path/to/micrograph.mrc', 'path set/get correctly')
    call assert_true(dfx == 0.25,   'dfx set/get correctly' )
    call assert_true(dfy == 2.89,   'dfy set/get correctly')
    call assert_true(ctfres == 8.9, 'ctfres set/get correctly' )
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
  end subroutine test_set_get_micrograph

  ! Verify that the micrograph serialise buffer is the expected size.
  subroutine test_serialise_micrograph()
    character(len=:),             allocatable :: buffer
    type(gui_metadata_micrograph)             :: meta
    write(*,'(A)') 'test_serialise_micrograph'
    call meta%new(GUI_METADATA_MICROGRAPH_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_MICROGRAPH_TYPE, 'type is set correctly')
    call meta%set(path=string('/test/path/to/micrograph.mrc'), dfx=0.25, dfy=2.89, ctfres=8.9, i_max=1, i=1)
    call assert_true(meta%assigned(), 'metadata object is set')
    call meta%serialise(buffer=buffer)
    call assert_true(allocated(buffer), 'buffer allocated')
    call assert_int(len(buffer), int(sizeof(meta), kind=4), 'buffer correct size')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
  end subroutine test_serialise_micrograph

  ! Verify the micrograph JSON output via buffer length and FNV-1a hash.
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
    call meta%set(path=string('/test/path/to/micrograph.mrc'), dfx=0.25, dfy=2.89, ctfres=8.9, i_max=1, i=1)
    call assert_true(meta%assigned(), 'metadata object is set')
    json_ptr => meta%jsonise()
    call json%print_to_string(json_ptr, buffer)
    call assert_true(allocated(buffer), 'buffer allocated')
    call assert_int(len(buffer), 130, 'buffer correct size')
    json_str = buffer
    call assert_int(json_str%strlen(), 130, 'string correct size')
    json_hash = json_str%to_fnv1a_hash64()
    call assert_char(json_hash%to_char(), 'ACEEF0BD08C90E99', 'correct checksum')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
    call json%destroy(json_ptr)
    call assert_true(.not.json%failed(), 'json destroyed')
  end subroutine test_jsonise_micrograph

  !---------------- histogram ----------------

  ! Verify that all histogram fields round-trip through set/get.
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
    call assert_true(meta%get(name=name, labels=labels, data=data), 'metadata retrieved')
    call assert_char(name%to_char(), 'histname', 'name set/get correctly')
    call assert_true(size(labels) == 10,   'labels correctly allocated')
    call assert_true(size(data)   == 10,   'data correctly allocated')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
    deallocate(labels, data)
  end subroutine test_set_get_histogram

  ! Verify that the histogram serialise buffer is the expected size.
  subroutine test_serialise_histogram()
    character(len=:),             allocatable :: buffer
    type(gui_metadata_histogram)              :: meta
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

  ! Verify the histogram JSON output via buffer length and FNV-1a hash.
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
    call json%initialize(no_whitespace=.true., compact_reals=.true.)
    call meta%new(GUI_METADATA_HISTOGRAM_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_HISTOGRAM_TYPE, 'type is set correctly')
    call meta%set(name=string('histname'), labels=labels, data=data)
    call assert_true(meta%assigned(), 'metadata object is set')
    json_ptr => meta%jsonise()
    call json%print_to_string(json_ptr, buffer)
    call assert_true(allocated(buffer), 'buffer allocated')
    call assert_int(len(buffer), 207, 'buffer correct size')
    json_str = buffer
    call assert_int(json_str%strlen(), 207, 'string correct size')
    json_hash = json_str%to_fnv1a_hash64()
    call assert_char(json_hash%to_char(), 'D6E335EB6331B91D', 'correct checksum')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
    call json%destroy(json_ptr)
    call assert_true(.not.json%failed(), 'json destroyed')
    deallocate(labels, data)
  end subroutine test_jsonise_histogram

  !---------------- timeplot ----------------

  ! Verify that all timeplot fields round-trip through set/get.
  subroutine test_set_get_timeplot()
    type(gui_metadata_timeplot)              :: meta
    type(string)                             :: name
    real,                        allocatable :: labels(:), data(:), data2(:)
    integer                                  :: i
    write(*,'(A)') 'test_set_get_timeplot'
    allocate(labels(10), data(10), data2(10))
    do i=1, 10
      labels(i) = 1.0/i
      data(i)   = i
      data2(i)  = i
    enddo
    call meta%new(GUI_METADATA_TIMEPLOT_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_TIMEPLOT_TYPE, 'type is set correctly')
    call meta%set(name=string('timeplotname'), labels=labels, data=data, data2=data2)
    call assert_true(meta%assigned(), 'metadata object is set')
    deallocate(labels, data, data2)
    call assert_true(meta%get(name=name, labels=labels, data=data, data2=data2), 'metadata retrieved')
    call assert_char(name%to_char(), 'timeplotname', 'name set/get correctly')
    call assert_true(size(labels) == 10,   'labels correctly allocated')
    call assert_true(size(data)   == 10,   'data correctly allocated')
    call assert_true(size(data2)  == 10,   'data2 correctly allocated')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
    deallocate(labels, data, data2)
  end subroutine test_set_get_timeplot

  ! Verify that the timeplot serialise buffer is the expected size.
  subroutine test_serialise_timeplot()
    character(len=:),             allocatable :: buffer
    type(gui_metadata_timeplot)               :: meta
    real,                         allocatable :: labels(:), data(:), data2(:)
    integer                                   :: i
    write(*,'(A)') 'test_serialise_timeplot'
    allocate(labels(10), data(10), data2(10))
    do i=1, 10
      labels(i) = 1.0/i
      data(i)   = i
      data2(i)  = i
    enddo
    call meta%new(GUI_METADATA_TIMEPLOT_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_TIMEPLOT_TYPE, 'type is set correctly')
    call meta%set(name=string('timeplotname'), labels=labels, data=data, data2=data2)
    call assert_true(meta%assigned(), 'metadata object is set')
    call meta%serialise(buffer=buffer)
    call assert_true(allocated(buffer), 'buffer allocated')
    call assert_int(len(buffer), int(sizeof(meta), kind=4), 'buffer correct size')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
    deallocate(labels, data, data2)
  end subroutine test_serialise_timeplot

  ! Verify the timeplot JSON output via buffer length and FNV-1a hash.
  subroutine test_jsonise_timeplot()
    character(kind=CK, len=:),    allocatable :: buffer
    type(gui_metadata_timeplot)               :: meta
    type(json_core)                           :: json
    type(json_value),             pointer     :: json_ptr
    type(string)                              :: json_str, json_hash
    real,                         allocatable :: labels(:), data(:), data2(:)
    integer                                   :: i
    write(*,'(A)') 'test_jsonise_timeplot'
    allocate(labels(10), data(10), data2(10))
    do i=1, 10
      labels(i) = 1.0/i
      data(i)   = i
      data2(i)  = i
    enddo
    call json%initialize(no_whitespace=.true., compact_reals=.true.)
    call meta%new(GUI_METADATA_TIMEPLOT_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_TIMEPLOT_TYPE, 'type is set correctly')
    call meta%set(name=string('histname'), labels=labels, data=data, data2=data2)
    call assert_true(meta%assigned(), 'metadata object is set')
    json_ptr => meta%jsonise()
    call json%print_to_string(json_ptr, buffer)
    call assert_true(allocated(buffer), 'buffer allocated')
    call assert_int(len(buffer), 336, 'buffer correct size')
    json_str = buffer
    call assert_int(json_str%strlen(), 336, 'string correct size')
    json_hash = json_str%to_fnv1a_hash64()
    call assert_char(json_hash%to_char(), '153A8CAA0C02A5AD', 'correct checksum')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
    call json%destroy(json_ptr)
    call assert_true(.not.json%failed(), 'json destroyed')
    deallocate(labels, data, data2)
  end subroutine test_jsonise_timeplot

  !---------------- optics group ----------------

  ! Verify that all optics-group fields round-trip through set/get.
  subroutine test_set_get_optics_group()
    type(gui_metadata_optics_group)          :: meta
    real,                        allocatable :: xshifts(:), yshifts(:)
    integer                                  :: i, i_max, n_shifts
    write(*,'(A)') 'test_set_get_optics_group'
    allocate(xshifts(10), yshifts(10))
    do i=1, 10
      xshifts(i) = 1.0/i
      yshifts(i) = i
    enddo
    call meta%new(GUI_METADATA_OPTICS_GROUP_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_OPTICS_GROUP_TYPE, 'type is set correctly')
    call meta%set(i=1, i_max=1, xshifts=xshifts, yshifts=yshifts, n_shifts=size(xshifts))
    call assert_true(meta%assigned(), 'metadata object is set')
    deallocate(xshifts, yshifts)
    call assert_true(meta%get(i=i, i_max=i_max, xshifts=xshifts, yshifts=yshifts, n_shifts=n_shifts), 'metadata retrieved')
    call assert_int(size(xshifts), 10, 'xshifts correctly allocated')
    call assert_int(size(yshifts), 10, 'yshifts correctly allocated')
    call assert_int(n_shifts,      10,  'n_shifts has correct value')
    call assert_int(i,              1,         'i has correct value')
    call assert_int(i_max,          1,     'i_max has correct value')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
    deallocate(xshifts, yshifts)
  end subroutine test_set_get_optics_group

  ! Verify that the optics-group serialise buffer is the expected size.
  subroutine test_serialise_optics_group()
    type(gui_metadata_optics_group)          :: meta
    character(kind=CK, len=:),   allocatable :: buffer
    real,                        allocatable :: xshifts(:), yshifts(:)
    integer                                  :: i
    write(*,'(A)') 'test_serialise_optics_group'
    allocate(xshifts(10), yshifts(10))
    do i=1, 10
      xshifts(i) = 1.0/i
      yshifts(i) = i
    enddo
    call meta%new(GUI_METADATA_OPTICS_GROUP_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_OPTICS_GROUP_TYPE, 'type is set correctly')
    call meta%set(i=1, i_max=1, xshifts=xshifts, yshifts=yshifts, n_shifts=size(xshifts))
    call assert_true(meta%assigned(), 'metadata object is set')
    call meta%serialise(buffer=buffer)
    call assert_true(allocated(buffer), 'buffer allocated')
    call assert_int(len(buffer), int(sizeof(meta), kind=4), 'buffer correct size')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
    deallocate(xshifts, yshifts)
  end subroutine test_serialise_optics_group

  ! Verify the optics-group JSON output via buffer length and FNV-1a hash.
  subroutine test_jsonise_optics_group()
    character(kind=CK, len=:),    allocatable :: buffer
    type(gui_metadata_optics_group)           :: meta
    type(json_core)                           :: json
    type(json_value),             pointer     :: json_ptr
    type(string)                              :: json_str, json_hash
    real,                        allocatable :: xshifts(:), yshifts(:)
    integer                                  :: i
    write(*,'(A)') 'test_jsonise_optics_group'
    allocate(xshifts(10), yshifts(10))
    do i=1, 10
      xshifts(i) = 1.0/i
      yshifts(i) = i
    enddo
    call json%initialize(no_whitespace=.true., compact_reals=.true.)
    call meta%new(GUI_METADATA_OPTICS_GROUP_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_OPTICS_GROUP_TYPE, 'type is set correctly')
    call meta%set(i=1, i_max=1, xshifts=xshifts, yshifts=yshifts, n_shifts=size(xshifts))
    call assert_true(meta%assigned(), 'metadata object is set')
    json_ptr => meta%jsonise()
    call json%print_to_string(json_ptr, buffer)
    call assert_true(allocated(buffer), 'buffer allocated')
    call assert_int(len(buffer), 359, 'buffer correct size')
    json_str = buffer
    call assert_int(json_str%strlen(), 359, 'string correct size')
    json_hash = json_str%to_fnv1a_hash64()
    call assert_char(json_hash%to_char(), '9052A6F0350F6128', 'correct checksum')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
    call json%destroy(json_ptr)
    call assert_true(.not.json%failed(), 'json destroyed')
    deallocate(xshifts, yshifts)
  end subroutine test_jsonise_optics_group

  !---------------- stream preprocess ----------------

  ! Verify that all stream-preprocess fields round-trip through set/get.
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
                              cutoff_astigmatism=cutoff_astigmatism), 'metadata retrieved')
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

  ! Verify that the stream-preprocess serialise buffer is the expected size.
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

  ! Verify the stream-preprocess JSON output via buffer length and FNV-1a hash.
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

  !---------------- stream optics assignment ----------------

  ! Verify that all stream optics-assignment fields round-trip through set/get.
  subroutine test_set_get_stream_optics_assignment()
    type(gui_metadata_stream_optics_assignment)      :: meta
    type(string)                                     :: stage
    integer                                          :: micrographs_assigned, optics_groups_assigned, last_micrograph_imported
    write(*,'(A)') 'test_set_get_stream_optics_assignment'
    call meta%new(GUI_METADATA_STREAM_OPTICS_ASSIGNMENT_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_STREAM_OPTICS_ASSIGNMENT_TYPE, 'type is set correctly')
    call meta%set(stage=string('test stage'), micrographs_assigned=100, optics_groups_assigned=50, last_micrograph_imported=1234)
    call assert_true(meta%assigned(), 'metadata object is set')
    call assert_true(meta%get(stage=stage, micrographs_assigned=micrographs_assigned, &
                     optics_groups_assigned=optics_groups_assigned, last_micrograph_imported=last_micrograph_imported), 'metadata retrieved')
    call assert_char(stage%to_char(), 'test stage', 'stage set/get correctly')
    call assert_int(micrographs_assigned,     100,  'micrographs_assigned set/get correctly')
    call assert_int(optics_groups_assigned,   50,   'optics_groups_assigned set/get correctly')
    call assert_int(last_micrograph_imported, 1234, 'last_micrograph_imported set/get correctly')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
  end subroutine test_set_get_stream_optics_assignment

  ! Verify that the stream optics-assignment serialise buffer is the expected size.
  subroutine test_serialise_stream_optics_assignment()
    type(gui_metadata_stream_optics_assignment) :: meta
    character(len=:),               allocatable :: buffer
    write(*,'(A)') 'test_serialise_stream_optics_assignment'
    call meta%new(GUI_METADATA_STREAM_OPTICS_ASSIGNMENT_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_STREAM_OPTICS_ASSIGNMENT_TYPE, 'type is set correctly')
    call meta%set(stage=string('test stage'), micrographs_assigned=100, optics_groups_assigned=50, last_micrograph_imported=1234)
    call assert_true(meta%assigned(), 'metadata object is set')
    call meta%serialise(buffer=buffer)
    call assert_true(allocated(buffer), 'buffer allocated')
    call assert_int(len(buffer), int(sizeof(meta), kind=4), 'buffer correct size')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
    deallocate(buffer)
  end subroutine test_serialise_stream_optics_assignment

  ! Verify the stream optics-assignment JSON output via buffer length and FNV-1a hash.
  subroutine test_jsonise_stream_optics_assignment()
    character(kind=CK, len=:),    allocatable   :: buffer
    type(gui_metadata_stream_optics_assignment) :: meta
    type(json_core)                             :: json
    type(json_value),             pointer       :: json_ptr
    type(string)                                :: json_str, json_hash
    write(*,'(A)') 'test_jsonise_stream_optics_assignment'
    call json%initialize(no_whitespace=.true., compact_reals=.true.)
    call meta%new(GUI_METADATA_STREAM_OPTICS_ASSIGNMENT_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_STREAM_OPTICS_ASSIGNMENT_TYPE, 'type is set correctly')
    call meta%set(stage=string('test stage'), micrographs_assigned=100, optics_groups_assigned=50, last_micrograph_imported=1234)
    call assert_true(meta%assigned(), 'metadata object is set')
    json_ptr => meta%jsonise()
    call json%print_to_string(json_ptr, buffer)
    call assert_true(allocated(buffer), 'buffer allocated')
    call assert_int(len(buffer), 109, 'buffer correct size')
    json_str = buffer
    call assert_int(json_str%strlen(), 109, 'string correct size')
    json_hash = json_str%to_fnv1a_hash64()
    call assert_char(json_hash%to_char(), '6C8E8AF0F9ACC713', 'correct checksum')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
    call json%destroy(json_ptr)
    call assert_true(.not.json%failed(), 'json destroyed')
    deallocate(buffer)
  end subroutine test_jsonise_stream_optics_assignment

  !---------------- stream update ----------------

  ! Verify that each stream-update threshold field can be set and retrieved independently.
  subroutine test_set_get_stream_update()
    type(gui_metadata_stream_update) :: meta
    write(*,'(A)') 'test_set_get_stream_update'
    call meta%new(GUI_METADATA_STREAM_UPDATE_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_STREAM_UPDATE_TYPE, 'type is set correctly')
    call meta%set_ctfres_update(8.5)
    call assert_true(meta%assigned(), 'metadata object is set after set_ctfres_update')
    call assert_true(meta%get_ctfres_update() == 8.5, 'ctfres_update set/get correctly')
    call meta%set_astigmatism_update(0.4)
    call assert_true(meta%get_astigmatism_update() == 0.4, 'astigmatism_update set/get correctly')
    call meta%set_icescore_update(0.3)
    call assert_true(meta%get_icescore_update() == 0.3, 'icescore_update set/get correctly')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
  end subroutine test_set_get_stream_update

  ! Verify that the stream-update serialise buffer is the expected size.
  subroutine test_serialise_stream_update()
    character(len=:),                allocatable :: buffer
    type(gui_metadata_stream_update)             :: meta
    write(*,'(A)') 'test_serialise_stream_update'
    call meta%new(GUI_METADATA_STREAM_UPDATE_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_STREAM_UPDATE_TYPE, 'type is set correctly')
    call meta%set_ctfres_update(8.5)
    call meta%set_astigmatism_update(0.4)
    call meta%set_icescore_update(0.3)
    call assert_true(meta%assigned(), 'metadata object is set')
    call meta%serialise(buffer=buffer)
    call assert_true(allocated(buffer), 'buffer allocated')
    call assert_int(len(buffer), int(sizeof(meta), kind=4), 'buffer correct size')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
    deallocate(buffer)
  end subroutine test_serialise_stream_update

  ! gui_metadata_stream_update has no jsonise override; the base class returns an
  ! empty JSON object ({}) when assigned.
  subroutine test_jsonise_stream_update()
    character(kind=CK, len=:),       allocatable :: buffer
    type(gui_metadata_stream_update)             :: meta
    type(json_core)                              :: json
    type(json_value),                pointer     :: json_ptr
    write(*,'(A)') 'test_jsonise_stream_update'
    call json%initialize(no_whitespace=.true., compact_reals=.true.)
    call meta%new(GUI_METADATA_STREAM_UPDATE_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call meta%set_ctfres_update(8.5)
    call assert_true(meta%assigned(), 'metadata object is set')
    json_ptr => meta%jsonise()
    call assert_true(associated(json_ptr), 'json pointer is associated')
    call json%print_to_string(json_ptr, buffer)
    call assert_true(allocated(buffer), 'buffer allocated')
    call assert_int(len(buffer), 2, 'base jsonise produces empty object ({})')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
    call json%destroy(json_ptr)
    call assert_true(.not.json%failed(), 'json destroyed')
    deallocate(buffer)
  end subroutine test_jsonise_stream_update

  !---------------- stream initial picking ----------------

  ! Verify that all stream initial-picking fields round-trip, including
  ! derived fields (micrographs_rejected, particles_per_mic) and that
  ! last_micrograph_imported is populated automatically on first set.
  subroutine test_set_get_stream_initial_picking()
    type(gui_metadata_stream_initial_picking) :: meta
    type(string)                              :: stage
    integer :: micrographs_imported, micrographs_accepted, micrographs_rejected
    integer :: particles_extracted, particles_per_mic, last_micrograph_imported
    write(*,'(A)') 'test_set_get_stream_initial_picking'
    call meta%new(GUI_METADATA_STREAM_INITIAL_PICKING_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_STREAM_INITIAL_PICKING_TYPE, 'type is set correctly')
    call meta%set(stage=string('test stage'), micrographs_imported=200, &
                  micrographs_accepted=180, particles_extracted=36000)
    call assert_true(meta%assigned(), 'metadata object is set')
    call assert_true(meta%get(stage=stage, micrographs_imported=micrographs_imported,             &
                              micrographs_accepted=micrographs_accepted,                          &
                              micrographs_rejected=micrographs_rejected,                          &
                              particles_extracted=particles_extracted,                             &
                              particles_per_mic=particles_per_mic,                                &
                              last_micrograph_imported=last_micrograph_imported), 'metadata retrieved')
    call assert_char(stage%to_char(), 'test stage',      'stage set/get correctly')
    call assert_int(micrographs_imported,  200,          'micrographs_imported set/get correctly')
    call assert_int(micrographs_accepted,  180,          'micrographs_accepted set/get correctly')
    call assert_int(micrographs_rejected,  20,           'micrographs_rejected derived correctly')
    call assert_int(particles_extracted,   36000,        'particles_extracted set/get correctly')
    call assert_int(particles_per_mic,     180,          'particles_per_mic derived correctly')
    call assert_true(last_micrograph_imported > 0,       'last_micrograph_imported is set')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
  end subroutine test_set_get_stream_initial_picking

  ! Verify that the stream initial-picking serialise buffer is the expected size.
  subroutine test_serialise_stream_initial_picking()
    character(len=:),                         allocatable :: buffer
    type(gui_metadata_stream_initial_picking)             :: meta
    write(*,'(A)') 'test_serialise_stream_initial_picking'
    call meta%new(GUI_METADATA_STREAM_INITIAL_PICKING_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_STREAM_INITIAL_PICKING_TYPE, 'type is set correctly')
    call meta%set(stage=string('test stage'), micrographs_imported=200, &
                  micrographs_accepted=180, particles_extracted=36000)
    call assert_true(meta%assigned(), 'metadata object is set')
    call meta%serialise(buffer=buffer)
    call assert_true(allocated(buffer), 'buffer allocated')
    call assert_int(len(buffer), int(sizeof(meta), kind=4), 'buffer correct size')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
    deallocate(buffer)
  end subroutine test_serialise_stream_initial_picking

  ! jsonise embeds last_micrograph_imported (live Unix timestamp) so only
  ! non-emptiness is checked — a stable hash comparison is not possible.
  subroutine test_jsonise_stream_initial_picking()
    character(kind=CK, len=:),                allocatable :: buffer
    type(gui_metadata_stream_initial_picking)              :: meta
    type(json_core)                                        :: json
    type(json_value),                          pointer     :: json_ptr
    write(*,'(A)') 'test_jsonise_stream_initial_picking'
    call json%initialize(no_whitespace=.true., compact_reals=.true.)
    call meta%new(GUI_METADATA_STREAM_INITIAL_PICKING_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call meta%set(stage=string('test stage'), micrographs_imported=200, &
                  micrographs_accepted=180, particles_extracted=36000)
    call assert_true(meta%assigned(), 'metadata object is set')
    json_ptr => meta%jsonise()
    call assert_true(associated(json_ptr), 'json pointer is associated')
    call json%print_to_string(json_ptr, buffer)
    call assert_true(allocated(buffer), 'buffer allocated')
    call assert_true(len(buffer) > 0, 'json output is non-empty')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
    call json%destroy(json_ptr)
    call assert_true(.not.json%failed(), 'json destroyed')
    deallocate(buffer)
  end subroutine test_jsonise_stream_initial_picking

  !---------------- stream opening2D ----------------

  ! Verify that all stream opening-2D fields round-trip, including
  ! the separately set user_input flag and the auto-populated
  ! last_particles_imported timestamp.
  subroutine test_set_get_stream_opening2D()
    type(gui_metadata_stream_opening2D) :: meta
    type(string)                        :: stage
    integer :: particles_imported, particles_accepted, last_particles_imported
    write(*,'(A)') 'test_set_get_stream_opening2D'
    call meta%new(GUI_METADATA_STREAM_OPENING2D_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_STREAM_OPENING2D_TYPE, 'type is set correctly')
    call meta%set(stage=string('test stage'), particles_imported=50000, particles_accepted=42000, &
                  mask_diam=160, box_size=256, mask_scale=0.75)
    call meta%set_user_input(.true.)
    call assert_true(meta%assigned(), 'metadata object is set')
    call assert_true(meta%get(stage=stage, particles_imported=particles_imported,           &
                              particles_accepted=particles_accepted,                        &
                              last_particles_imported=last_particles_imported), 'metadata retrieved')
    call assert_char(stage%to_char(), 'test stage', 'stage set/get correctly')
    call assert_int(particles_imported,  50000,     'particles_imported set/get correctly')
    call assert_int(particles_accepted,  42000,     'particles_accepted set/get correctly')
    call assert_true(last_particles_imported > 0,   'last_particles_imported is set')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
  end subroutine test_set_get_stream_opening2D

  ! Verify that the stream opening-2D serialise buffer is the expected size.
  subroutine test_serialise_stream_opening2D()
    character(len=:),                   allocatable :: buffer
    type(gui_metadata_stream_opening2D)             :: meta
    write(*,'(A)') 'test_serialise_stream_opening2D'
    call meta%new(GUI_METADATA_STREAM_OPENING2D_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_STREAM_OPENING2D_TYPE, 'type is set correctly')
    call meta%set(stage=string('test stage'), particles_imported=50000, particles_accepted=42000, &
                  mask_diam=160, box_size=256, mask_scale=0.75)
    call assert_true(meta%assigned(), 'metadata object is set')
    call meta%serialise(buffer=buffer)
    call assert_true(allocated(buffer), 'buffer allocated')
    call assert_int(len(buffer), int(sizeof(meta), kind=4), 'buffer correct size')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
    deallocate(buffer)
  end subroutine test_serialise_stream_opening2D

  ! jsonise embeds last_particles_imported (live Unix timestamp) so only
  ! non-emptiness is checked — a stable hash comparison is not possible.
  subroutine test_jsonise_stream_opening2D()
    character(kind=CK, len=:),          allocatable :: buffer
    type(gui_metadata_stream_opening2D)             :: meta
    type(json_core)                                 :: json
    type(json_value),                   pointer     :: json_ptr
    write(*,'(A)') 'test_jsonise_stream_opening2D'
    call json%initialize(no_whitespace=.true., compact_reals=.true.)
    call meta%new(GUI_METADATA_STREAM_OPENING2D_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call meta%set(stage=string('test stage'), particles_imported=50000, particles_accepted=42000, &
                  mask_diam=160, box_size=256, mask_scale=0.75)
    call assert_true(meta%assigned(), 'metadata object is set')
    json_ptr => meta%jsonise()
    call assert_true(associated(json_ptr), 'json pointer is associated')
    call json%print_to_string(json_ptr, buffer)
    call assert_true(allocated(buffer), 'buffer allocated')
    call assert_true(len(buffer) > 0, 'json output is non-empty')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
    call json%destroy(json_ptr)
    call assert_true(.not.json%failed(), 'json destroyed')
    deallocate(buffer)
  end subroutine test_jsonise_stream_opening2D

end module simple_gui_metadata_tester