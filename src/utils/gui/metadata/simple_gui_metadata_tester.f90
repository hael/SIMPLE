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
!     gui_metadata_timeplot, gui_metadata_optics_group, gui_metadata_cavg2D,
!     gui_metadata_stream_preprocess, gui_metadata_stream_optics_assignment,
!     gui_metadata_stream_update, gui_metadata_stream_picking (initial and
!     reference picking), gui_metadata_stream_opening2D,
!     gui_metadata_stream_particle_sieving, gui_metadata_stream_pool2D,
!     gui_metadata_stream_pool2D_snapshot
!
! ENTRY POINT:
!   run_all_gui_metadata_tests() — run every test in this module
!
! DEPENDENCIES:
!   simple_test_utils, simple_gui_metadata_api
!==============================================================================
module simple_gui_metadata_tester
  use simple_test_utils,       only: assert_true, assert_int, assert_char
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
    call test_set_get_cavg2D()
    call test_serialise_cavg2D()
    call test_jsonise_cavg2D()
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
    call test_set_get_stream_reference_picking()
    call test_serialise_stream_reference_picking()
    call test_jsonise_stream_reference_picking()
    call test_set_get_stream_opening2D()
    call test_serialise_stream_opening2D()
    call test_jsonise_stream_opening2D()
    call test_set_get_stream_particle_sieving()
    call test_serialise_stream_particle_sieving()
    call test_jsonise_stream_particle_sieving()
    call test_set_get_stream_pool2D()
    call test_serialise_stream_pool2D()
    call test_jsonise_stream_pool2D()
    call test_set_get_stream_pool2D_snapshot()
    call test_serialise_stream_pool2D_snapshot()
    call test_jsonise_stream_pool2D_snapshot()
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

  !---------------- cavg2D ----------------

  ! Verify that all cavg2D fields round-trip through set/get, including optional res and pop.
  subroutine test_set_get_cavg2D()
    type(gui_metadata_cavg2D) :: meta
    type(string)              :: path, mrcpath
    integer                   :: idx
    type(sprite_sheet_pos)    :: sprite
    real                      :: res
    integer                   :: pop
    write(*,'(A)') 'test_set_get_cavg2D'
    call meta%new(GUI_METADATA_CAVG2D_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_CAVG2D_TYPE, 'type is set correctly')
    call meta%set(path=string('/test/path/to/cavg.jpeg'), mrcpath=string('/test/path/to/cavg.mrc'), &
                  idx=3, sprite=sprite_sheet_pos(x=25.0, y=50.0, h=800, w=1200), &
                  i=1, i_max=1, res=3.5, pop=1500)
    call assert_true(meta%assigned(), 'metadata object is set')
    call assert_int(meta%get_idx(), 3, 'get_idx returns correct value')
    call assert_true(meta%get(path=path, mrcpath=mrcpath, idx=idx, sprite=sprite, res=res, pop=pop), 'metadata retrieved')
    call assert_char(path%to_char(),    '/test/path/to/cavg.jpeg', 'path set/get correctly')
    call assert_char(mrcpath%to_char(), '/test/path/to/cavg.mrc',  'mrcpath set/get correctly')
    call assert_int(idx,      3,    'idx set/get correctly')
    call assert_true(sprite%x == 25.0, 'spritex set/get correctly')
    call assert_true(sprite%y == 50.0, 'spritey set/get correctly')
    call assert_int(sprite%h, 800,  'spriteh set/get correctly')
    call assert_int(sprite%w, 1200, 'spritew set/get correctly')
    call assert_true(res == 3.5,   'res set/get correctly')
    call assert_int(pop,     1500, 'pop set/get correctly')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
  end subroutine test_set_get_cavg2D

  ! Verify that the cavg2D serialise buffer is the expected size.
  subroutine test_serialise_cavg2D()
    character(len=:),         allocatable :: buffer
    type(gui_metadata_cavg2D)             :: meta
    write(*,'(A)') 'test_serialise_cavg2D'
    call meta%new(GUI_METADATA_CAVG2D_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_CAVG2D_TYPE, 'type is set correctly')
    call meta%set(path=string('/test/path/to/cavg.jpeg'), mrcpath=string('/test/path/to/cavg.mrc'), &
                  idx=3, sprite=sprite_sheet_pos(x=25.0, y=50.0, h=800, w=1200), &
                  i=1, i_max=1, res=3.5, pop=1500)
    call assert_true(meta%assigned(), 'metadata object is set')
    call meta%serialise(buffer=buffer)
    call assert_true(allocated(buffer), 'buffer allocated')
    call assert_int(len(buffer), int(sizeof(meta), kind=4), 'buffer correct size')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
  end subroutine test_serialise_cavg2D

  ! Verify the cavg2D JSON output via buffer length and FNV-1a hash.
  subroutine test_jsonise_cavg2D()
    character(kind=CK, len=:), allocatable :: buffer
    type(gui_metadata_cavg2D)              :: meta
    type(json_core)                        :: json
    type(json_value),          pointer     :: json_ptr
    type(string)                           :: json_str, json_hash
    write(*,'(A)') 'test_jsonise_cavg2D'
    call json%initialize(no_whitespace=.true., compact_reals=.true.)
    call meta%new(GUI_METADATA_CAVG2D_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_CAVG2D_TYPE, 'type is set correctly')
    call meta%set(path=string('/test/path/to/cavg.jpeg'), mrcpath=string('/test/path/to/cavg.mrc'), &
                  idx=3, sprite=sprite_sheet_pos(x=25.0, y=50.0, h=800, w=1200), &
                  i=1, i_max=1, res=3.5, pop=1500)
    call assert_true(meta%assigned(), 'metadata object is set')
    json_ptr => meta%jsonise()
    call json%print_to_string(json_ptr, buffer)
    call assert_true(allocated(buffer), 'buffer allocated')
    call assert_true(len(buffer) > 0, 'json output is non-empty')
    call assert_int(len(buffer), 166, 'buffer correct size')
    json_str = buffer
    call assert_int(json_str%strlen(), 166, 'string correct size')
    json_hash = json_str%to_fnv1a_hash64()
    call assert_char(json_hash%to_char(), 'BDAB38FA961B3B76', 'correct checksum')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
    call json%destroy(json_ptr)
    call assert_true(.not.json%failed(), 'json destroyed')
  end subroutine test_jsonise_cavg2D

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

  ! Verify that each stream-update threshold field can be set and retrieved independently,
  ! including pickrefs_selection, sieverefs_selection, and the snapshot2D compound field.
  subroutine test_set_get_stream_update()
    type(gui_metadata_stream_update) :: meta
    integer,          allocatable    :: sel_out(:)
    type(string)                     :: fname_out
    integer                          :: snap_id_out, iter_out
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
    call meta%set_mskdiam2D_update(180.0)
    call assert_true(meta%get_mskdiam2D_update() == 180.0, 'mskdiam2D_update set/get correctly')
    ! pickrefs_selection
    call assert_int(meta%get_pickrefs_selection_length(), 0, 'pickrefs_selection_length zero before set')
    call meta%set_pickrefs_selection([1, 0, 1, 1, 0, 1])
    call assert_int(meta%get_pickrefs_selection_length(), 6, 'pickrefs_selection_length set/get correctly')
    sel_out = meta%get_pickrefs_selection()
    call assert_int(size(sel_out), 6,   'pickrefs_selection size correct')
    call assert_int(sel_out(1),    1,   'pickrefs_selection(1) correct')
    call assert_int(sel_out(2),    0,   'pickrefs_selection(2) correct')
    call assert_int(sel_out(6),    1,   'pickrefs_selection(6) correct')
    deallocate(sel_out)
    ! sieverefs_selection
    call assert_int(meta%get_sieverefs_selection_length(), 0, 'sieverefs_selection_length zero before set')
    call meta%set_sieverefs_selection([3, 7, 42, 100])
    call assert_int(meta%get_sieverefs_selection_length(), 4, 'sieverefs_selection_length set/get correctly')
    sel_out = meta%get_sieverefs_selection()
    call assert_int(size(sel_out), 4,   'sieverefs_selection size correct')
    call assert_int(sel_out(1),    3,   'sieverefs_selection(1) correct')
    call assert_int(sel_out(2),    7,   'sieverefs_selection(2) correct')
    call assert_int(sel_out(4),    100, 'sieverefs_selection(4) correct')
    deallocate(sel_out)
    call assert_true(.not.meta%has_snapshot2D_update(), 'has_snapshot2D_update false before set')
    call meta%set_snapshot2D_update(snapshot_id=3, iteration=5, &
                                    selection=[23,171,200,46,142], filename=string('snapshot_3.simple'))
    call assert_true(meta%has_snapshot2D_update(), 'has_snapshot2D_update true after set')
    call meta%get_snapshot2D_update(snap_id_out, iter_out, sel_out, fname_out)
    call assert_int(snap_id_out,      3,                   'snapshot2D_id set/get correctly')
    call assert_int(iter_out,         5,                   'snapshot2D_iteration set/get correctly')
    call assert_int(size(sel_out),    5,                   'snapshot2D_selection size correct')
    call assert_int(sel_out(1),       23,                  'snapshot2D_selection(1) correct')
    call assert_int(sel_out(5),       142,                 'snapshot2D_selection(5) correct')
    call assert_char(fname_out%to_char(),       'snapshot_3.simple', 'snapshot2D_filename set/get correctly')
    deallocate(sel_out)
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
  end subroutine test_set_get_stream_update

  ! Verify that the stream-update serialise buffer is the expected size, including
  ! pickrefs_selection, sieverefs_selection, and the snapshot2D fields.
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
    call meta%set_mskdiam2D_update(180.0)
    call meta%set_pickrefs_selection([1, 0, 1, 1, 0])
    call meta%set_sieverefs_selection([3, 7, 42])
    call meta%set_snapshot2D_update(snapshot_id=3, iteration=5, &
                                    selection=[23,171,200,46,142], filename=string('snapshot_3.simple'))
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
    type(gui_metadata_stream_picking) :: meta
    type(string)                              :: stage
    integer :: micrographs_imported, micrographs_accepted, micrographs_rejected
    integer :: particles_extracted, particles_per_mic, last_micrograph_imported, box_size
    write(*,'(A)') 'test_set_get_stream_initial_picking'
    call meta%new(GUI_METADATA_STREAM_INITIAL_PICKING_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_STREAM_INITIAL_PICKING_TYPE, 'type is set correctly')
    call meta%set(stage=string('test stage'), micrographs_imported=200, &
                  micrographs_accepted=180, particles_extracted=36000, box_size=256)
    call assert_true(meta%assigned(), 'metadata object is set')
    call assert_true(meta%get(stage=stage, micrographs_imported=micrographs_imported,             &
                              micrographs_accepted=micrographs_accepted,                          &
                              micrographs_rejected=micrographs_rejected,                          &
                              particles_extracted=particles_extracted,                             &
                              particles_per_mic=particles_per_mic,                                &
                              last_micrograph_imported=last_micrograph_imported,                  &
                              box_size=box_size), 'metadata retrieved')
    call assert_char(stage%to_char(), 'test stage',      'stage set/get correctly')
    call assert_int(micrographs_imported,  200,          'micrographs_imported set/get correctly')
    call assert_int(micrographs_accepted,  180,          'micrographs_accepted set/get correctly')
    call assert_int(micrographs_rejected,  20,           'micrographs_rejected derived correctly')
    call assert_int(particles_extracted,   36000,        'particles_extracted set/get correctly')
    call assert_int(particles_per_mic,     180,          'particles_per_mic derived correctly')
    call assert_true(last_micrograph_imported > 0,       'last_micrograph_imported is set')
    call assert_int(box_size,              256,          'box_size set/get correctly')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
  end subroutine test_set_get_stream_initial_picking

  ! Verify that the stream initial-picking serialise buffer is the expected size.
  subroutine test_serialise_stream_initial_picking()
    character(len=:),                         allocatable :: buffer
    type(gui_metadata_stream_picking)             :: meta
    write(*,'(A)') 'test_serialise_stream_initial_picking'
    call meta%new(GUI_METADATA_STREAM_INITIAL_PICKING_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_STREAM_INITIAL_PICKING_TYPE, 'type is set correctly')
    call meta%set(stage=string('test stage'), micrographs_imported=200, &
                  micrographs_accepted=180, particles_extracted=36000, box_size=256)
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
    type(gui_metadata_stream_picking)             :: meta
    type(json_core)                                       :: json
    type(json_value),                          pointer    :: json_ptr
    type(string)                                          :: json_str, json_hash
    logical                                               :: found
    write(*,'(A)') 'test_jsonise_stream_initial_picking'
    call json%initialize(no_whitespace=.true., compact_reals=.true.)
    call meta%new(GUI_METADATA_STREAM_INITIAL_PICKING_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call meta%set(stage=string('test stage'), micrographs_imported=200, &
                  micrographs_accepted=180, particles_extracted=36000, box_size=256)
    call assert_true(meta%assigned(), 'metadata object is set')
    json_ptr => meta%jsonise()
    call assert_true(associated(json_ptr), 'json pointer is associated')
    call json%update(json_ptr, 'last_micrograph_imported', 0, found)
    call json%print_to_string(json_ptr, buffer)
    call assert_true(allocated(buffer), 'buffer allocated')
    call assert_true(len(buffer) > 0, 'json output is non-empty')
    call assert_int(len(buffer), 198, 'buffer correct size')
    json_str = buffer
    call assert_int(json_str%strlen(), 198, 'string correct size')
    json_hash = json_str%to_fnv1a_hash64()
    call assert_char(json_hash%to_char(), 'E00E9BC45CFB31B4', 'correct checksum')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
    call json%destroy(json_ptr)
    call assert_true(.not.json%failed(), 'json destroyed')
    deallocate(buffer)
  end subroutine test_jsonise_stream_initial_picking

  !---------------- stream reference picking ----------------

  ! Verify that all stream reference-picking fields round-trip, including
  ! derived fields (micrographs_rejected, particles_per_mic) and that
  ! last_micrograph_imported is populated automatically on first set.
  subroutine test_set_get_stream_reference_picking()
    type(gui_metadata_stream_picking) :: meta
    type(string)                              :: stage
    integer :: micrographs_imported, micrographs_accepted, micrographs_rejected
    integer :: particles_extracted, particles_per_mic, last_micrograph_imported, box_size
    write(*,'(A)') 'test_set_get_stream_reference_picking'
    call meta%new(GUI_METADATA_STREAM_REFERENCE_PICKING_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_STREAM_REFERENCE_PICKING_TYPE, 'type is set correctly')
    call meta%set(stage=string('test stage'), micrographs_imported=200, &
                  micrographs_accepted=180, particles_extracted=36000, box_size=256)
    call assert_true(meta%assigned(), 'metadata object is set')
    call assert_true(meta%get(stage=stage, micrographs_imported=micrographs_imported,             &
                              micrographs_accepted=micrographs_accepted,                          &
                              micrographs_rejected=micrographs_rejected,                          &
                              particles_extracted=particles_extracted,                             &
                              particles_per_mic=particles_per_mic,                                &
                              last_micrograph_imported=last_micrograph_imported,                  &
                              box_size=box_size), 'metadata retrieved')
    call assert_char(stage%to_char(), 'test stage',      'stage set/get correctly')
    call assert_int(micrographs_imported,  200,          'micrographs_imported set/get correctly')
    call assert_int(micrographs_accepted,  180,          'micrographs_accepted set/get correctly')
    call assert_int(micrographs_rejected,  20,           'micrographs_rejected derived correctly')
    call assert_int(particles_extracted,   36000,        'particles_extracted set/get correctly')
    call assert_int(particles_per_mic,     180,          'particles_per_mic derived correctly')
    call assert_true(last_micrograph_imported > 0,       'last_micrograph_imported is set')
    call assert_int(box_size,              256,          'box_size set/get correctly')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
  end subroutine test_set_get_stream_reference_picking

  ! Verify that the stream reference-picking serialise buffer is the expected size.
  subroutine test_serialise_stream_reference_picking()
    character(len=:),                         allocatable :: buffer
    type(gui_metadata_stream_picking)                     :: meta
    write(*,'(A)') 'test_serialise_stream_reference_picking'
    call meta%new(GUI_METADATA_STREAM_REFERENCE_PICKING_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_STREAM_REFERENCE_PICKING_TYPE, 'type is set correctly')
    call meta%set(stage=string('test stage'), micrographs_imported=200, &
                  micrographs_accepted=180, particles_extracted=36000, box_size=256)
    call assert_true(meta%assigned(), 'metadata object is set')
    call meta%serialise(buffer=buffer)
    call assert_true(allocated(buffer), 'buffer allocated')
    call assert_int(len(buffer), int(sizeof(meta), kind=4), 'buffer correct size')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
    deallocate(buffer)
  end subroutine test_serialise_stream_reference_picking

  ! jsonise embeds last_micrograph_imported (live Unix timestamp) so only
  ! non-emptiness is checked — a stable hash comparison is not possible.
  subroutine test_jsonise_stream_reference_picking()
    character(kind=CK, len=:),                allocatable :: buffer
    type(gui_metadata_stream_picking)                     :: meta
    type(json_core)                                       :: json
    type(json_value),                          pointer    :: json_ptr
    type(string)                                          :: json_str, json_hash
    logical                                               :: found
    write(*,'(A)') 'test_jsonise_stream_reference_picking'
    call json%initialize(no_whitespace=.true., compact_reals=.true.)
    call meta%new(GUI_METADATA_STREAM_REFERENCE_PICKING_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call meta%set(stage=string('test stage'), micrographs_imported=200, &
                  micrographs_accepted=180, particles_extracted=36000, box_size=256)
    call assert_true(meta%assigned(), 'metadata object is set')
    json_ptr => meta%jsonise()
    call assert_true(associated(json_ptr), 'json pointer is associated')
    call json%update(json_ptr, 'last_micrograph_imported', 0, found)
    call json%print_to_string(json_ptr, buffer)
    call assert_true(allocated(buffer), 'buffer allocated')
    call assert_true(len(buffer) > 0, 'json output is non-empty')
    call assert_int(len(buffer), 198, 'buffer correct size')
    json_str = buffer
    call assert_int(json_str%strlen(), 198, 'string correct size')
    json_hash = json_str%to_fnv1a_hash64()
    call assert_char(json_hash%to_char(), 'E00E9BC45CFB31B4', 'correct checksum')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
    call json%destroy(json_ptr)
    call assert_true(.not.json%failed(), 'json destroyed')
    deallocate(buffer)
  end subroutine test_jsonise_stream_reference_picking

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
    type(string)                                    :: json_str, json_hash
    logical                                         :: found
    write(*,'(A)') 'test_jsonise_stream_opening2D'
    call json%initialize(no_whitespace=.true., compact_reals=.true.)
    call meta%new(GUI_METADATA_STREAM_OPENING2D_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call meta%set(stage=string('test stage'), particles_imported=50000, particles_accepted=42000, &
                  mask_diam=160, box_size=256, mask_scale=0.75)
    call assert_true(meta%assigned(), 'metadata object is set')
    json_ptr => meta%jsonise()
    call assert_true(associated(json_ptr), 'json pointer is associated')
    call json%update(json_ptr, 'last_particles_imported', 0, found)
    call json%print_to_string(json_ptr, buffer)
    call assert_true(allocated(buffer), 'buffer allocated')
    call assert_true(len(buffer) > 0, 'json output is non-empty')
    call assert_int(len(buffer), 199, 'buffer correct size')
    json_str = buffer
    call assert_int(json_str%strlen(), 199, 'string correct size')
    json_hash = json_str%to_fnv1a_hash64()
    call assert_char(json_hash%to_char(), 'D9735541046FBFA7', 'correct checksum')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
    call json%destroy(json_ptr)
    call assert_true(.not.json%failed(), 'json destroyed')
    deallocate(buffer)
  end subroutine test_jsonise_stream_opening2D

  !---------------- stream particle sieving ----------------

  ! Verify that all particle-sieving fields round-trip through set/get,
  ! including the separately set user_input flag and the ref-selection array.
  subroutine test_set_get_stream_particle_sieving()
    type(gui_metadata_stream_particle_sieving) :: meta
    type(string) :: stage
    integer      :: particles_imported, particles_accepted, particles_rejected, last_import_time
    logical      :: user_input
    write(*,'(A)') 'test_set_get_stream_particle_sieving'
    call meta%new(GUI_METADATA_STREAM_PARTICLE_SIEVING_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_STREAM_PARTICLE_SIEVING_TYPE, 'type is set correctly')
    call meta%set(stage=string('sieving'), particles_imported=10000, &
                  particles_accepted=8000, particles_rejected=2000)
    call meta%set_user_input(.true.)
    call meta%set_initial_ref_selection(3)
    call meta%set_initial_ref_selection(7)
    call assert_true(meta%assigned(), 'metadata object is set')
    call assert_true(meta%get(stage=stage, particles_imported=particles_imported,     &
                              particles_accepted=particles_accepted,                  &
                              particles_rejected=particles_rejected,                  &
                              last_import_time=last_import_time,                      &
                              user_input=user_input), 'metadata retrieved')
    call assert_char(stage%to_char(), 'sieving', 'stage set/get correctly')
    call assert_int(particles_imported,  10000,  'particles_imported set/get correctly')
    call assert_int(particles_accepted,   8000,  'particles_accepted set/get correctly')
    call assert_int(particles_rejected,   2000,  'particles_rejected set/get correctly')
    call assert_true(last_import_time > 0,       'last_import_time is set')
    call assert_true(user_input,                 'user_input set/get correctly')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
  end subroutine test_set_get_stream_particle_sieving

  ! Verify that the particle-sieving serialise buffer is the expected size.
  subroutine test_serialise_stream_particle_sieving()
    character(len=:),                          allocatable :: buffer
    type(gui_metadata_stream_particle_sieving)             :: meta
    write(*,'(A)') 'test_serialise_stream_particle_sieving'
    call meta%new(GUI_METADATA_STREAM_PARTICLE_SIEVING_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_STREAM_PARTICLE_SIEVING_TYPE, 'type is set correctly')
    call meta%set(stage=string('sieving'), particles_imported=10000, &
                  particles_accepted=8000, particles_rejected=2000)
    call assert_true(meta%assigned(), 'metadata object is set')
    call meta%serialise(buffer=buffer)
    call assert_true(allocated(buffer), 'buffer allocated')
    call assert_int(len(buffer), int(sizeof(meta), kind=4), 'buffer correct size')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
    deallocate(buffer)
  end subroutine test_serialise_stream_particle_sieving

  ! jsonise embeds last_import_time (live Unix timestamp) so the timestamp
  ! field is zeroed before hashing.
  subroutine test_jsonise_stream_particle_sieving()
    character(kind=CK, len=:),                 allocatable :: buffer
    type(gui_metadata_stream_particle_sieving)             :: meta
    type(json_core)                                        :: json
    type(json_value),                          pointer     :: json_ptr
    type(string)                                           :: json_str, json_hash
    logical                                                :: found
    write(*,'(A)') 'test_jsonise_stream_particle_sieving'
    call json%initialize(no_whitespace=.true., compact_reals=.true.)
    call meta%new(GUI_METADATA_STREAM_PARTICLE_SIEVING_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call meta%set(stage=string('sieving'), particles_imported=10000, &
                  particles_accepted=8000, particles_rejected=2000)
    call assert_true(meta%assigned(), 'metadata object is set')
    json_ptr => meta%jsonise()
    call assert_true(associated(json_ptr), 'json pointer is associated')
    call json%update(json_ptr, 'last_import_time', 0, found)
    call json%print_to_string(json_ptr, buffer)
    call assert_true(allocated(buffer), 'buffer allocated')
    call assert_true(len(buffer) > 0, 'json output is non-empty')
    json_str = buffer
    json_hash = json_str%to_fnv1a_hash64()
    call assert_char(json_hash%to_char(), '743600493199B0C3', 'json is stable')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
    call json%destroy(json_ptr)
    call assert_true(.not.json%failed(), 'json destroyed')
    deallocate(buffer)
  end subroutine test_jsonise_stream_particle_sieving

  !---------------- stream pool2D ----------------

  ! Verify that all pool-2D fields round-trip through set/get,
  ! including iteration, mask geometry, and the user_input flag.
  subroutine test_set_get_stream_pool2D()
    type(gui_metadata_stream_pool2D) :: meta
    type(string) :: stage
    integer      :: iteration, particles_imported, particles_accepted, particles_rejected
    integer      :: last_import_time, mskdiam
    real         :: mskscale
    logical      :: user_input
    write(*,'(A)') 'test_set_get_stream_pool2D'
    call meta%new(GUI_METADATA_STREAM_POOL2D_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_STREAM_POOL2D_TYPE, 'type is set correctly')
    call meta%set(stage=string('pool2D'), iteration=3, particles_imported=20000, &
                  particles_accepted=16000, particles_rejected=4000, mskdiam=180, mskscale=0.5)
    call meta%set_user_input(.true.)
    call meta%set_initial_ref_selection(2)
    call assert_true(meta%assigned(), 'metadata object is set')
    call assert_true(meta%get(stage=stage, iteration=iteration,                      &
                              particles_imported=particles_imported,                  &
                              particles_accepted=particles_accepted,                  &
                              particles_rejected=particles_rejected,                  &
                              last_import_time=last_import_time,                      &
                              user_input=user_input, mskdiam=mskdiam,                &
                              mskscale=mskscale), 'metadata retrieved')
    call assert_char(stage%to_char(), 'pool2D', 'stage set/get correctly')
    call assert_int(iteration,            3,     'iteration set/get correctly')
    call assert_int(particles_imported,  20000,  'particles_imported set/get correctly')
    call assert_int(particles_accepted,  16000,  'particles_accepted set/get correctly')
    call assert_int(particles_rejected,   4000,  'particles_rejected set/get correctly')
    call assert_true(last_import_time > 0,       'last_import_time is set')
    call assert_true(user_input,                 'user_input set/get correctly')
    call assert_int(mskdiam,              180,   'mskdiam set/get correctly')
    call assert_true(mskscale == 0.5,            'mskscale set/get correctly')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
  end subroutine test_set_get_stream_pool2D

  ! Verify that the pool-2D serialise buffer is the expected size.
  subroutine test_serialise_stream_pool2D()
    character(len=:),              allocatable :: buffer
    type(gui_metadata_stream_pool2D)           :: meta
    write(*,'(A)') 'test_serialise_stream_pool2D'
    call meta%new(GUI_METADATA_STREAM_POOL2D_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_STREAM_POOL2D_TYPE, 'type is set correctly')
    call meta%set(stage=string('pool2D'), iteration=3, particles_imported=20000, &
                  particles_accepted=16000, particles_rejected=4000, mskdiam=180, mskscale=0.5)
    call assert_true(meta%assigned(), 'metadata object is set')
    call meta%serialise(buffer=buffer)
    call assert_true(allocated(buffer), 'buffer allocated')
    call assert_int(len(buffer), int(sizeof(meta), kind=4), 'buffer correct size')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
    deallocate(buffer)
  end subroutine test_serialise_stream_pool2D

  ! jsonise embeds last_import_time (live Unix timestamp) so the timestamp
  ! field is zeroed before hashing.
  subroutine test_jsonise_stream_pool2D()
    character(kind=CK, len=:),     allocatable :: buffer
    type(gui_metadata_stream_pool2D)           :: meta
    type(json_core)                            :: json
    type(json_value),              pointer     :: json_ptr
    type(string)                               :: json_str, json_hash
    logical                                    :: found
    write(*,'(A)') 'test_jsonise_stream_pool2D'
    call json%initialize(no_whitespace=.true., compact_reals=.true.)
    call meta%new(GUI_METADATA_STREAM_POOL2D_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call meta%set(stage=string('pool2D'), iteration=3, particles_imported=20000, &
                  particles_accepted=16000, particles_rejected=4000, mskdiam=180, mskscale=0.5)
    call assert_true(meta%assigned(), 'metadata object is set')
    json_ptr => meta%jsonise()
    call assert_true(associated(json_ptr), 'json pointer is associated')
    call json%update(json_ptr, 'last_import_time', 0, found)
    call json%print_to_string(json_ptr, buffer)
    call assert_true(allocated(buffer), 'buffer allocated')
    call assert_true(len(buffer) > 0, 'json output is non-empty')
    json_str = buffer
    json_hash = json_str%to_fnv1a_hash64()
    call assert_char(json_hash%to_char(), '6504A9531CD8C97D', 'json is stable')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
    call json%destroy(json_ptr)
    call assert_true(.not.json%failed(), 'json destroyed')
    deallocate(buffer)
  end subroutine test_jsonise_stream_pool2D

  !---------------- stream pool2D snapshot ----------------

  ! Verify that all pool-2D snapshot fields round-trip through set/get,
  ! and that snapshot_time is populated automatically on set.
  subroutine test_set_get_stream_pool2D_snapshot()
    type(gui_metadata_stream_pool2D_snapshot) :: meta
    type(string)                              :: fname_out
    integer                                   :: id_out, nptcls_out, time_out
    write(*,'(A)') 'test_set_get_stream_pool2D_snapshot'
    call meta%new(GUI_METADATA_STREAM_POOL2D_SNAPSHOT_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_STREAM_POOL2D_SNAPSHOT_TYPE, 'type is set correctly')
    call meta%set(id=3, snapshot_filename=string('snapshot_3.simple'), snapshot_nptcls=12500)
    call assert_true(meta%assigned(), 'metadata object is set')
    call assert_true(meta%get(id_out, fname_out, nptcls_out, time_out), 'metadata retrieved')
    call assert_int(id_out,       3,                  'id set/get correctly')
    call assert_char(fname_out%to_char(),  'snapshot_3.simple', 'snapshot_filename set/get correctly')
    call assert_int(nptcls_out,  12500,               'snapshot_nptcls set/get correctly')
    call assert_true(time_out > 0,                    'snapshot_time is set automatically')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
  end subroutine test_set_get_stream_pool2D_snapshot

  ! Verify that the pool-2D snapshot serialise buffer is the expected size.
  subroutine test_serialise_stream_pool2D_snapshot()
    character(len=:),                         allocatable :: buffer
    type(gui_metadata_stream_pool2D_snapshot)             :: meta
    write(*,'(A)') 'test_serialise_stream_pool2D_snapshot'
    call meta%new(GUI_METADATA_STREAM_POOL2D_SNAPSHOT_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call assert_int(meta%type(), GUI_METADATA_STREAM_POOL2D_SNAPSHOT_TYPE, 'type is set correctly')
    call meta%set(id=3, snapshot_filename=string('snapshot_3.simple'), snapshot_nptcls=12500)
    call assert_true(meta%assigned(), 'metadata object is set')
    call meta%serialise(buffer=buffer)
    call assert_true(allocated(buffer), 'buffer allocated')
    call assert_int(len(buffer), int(sizeof(meta), kind=4), 'buffer correct size')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
    deallocate(buffer)
  end subroutine test_serialise_stream_pool2D_snapshot

  ! jsonise embeds snapshot_time (live Unix timestamp) so the timestamp
  ! field is zeroed before hashing to make the output deterministic.
  subroutine test_jsonise_stream_pool2D_snapshot()
    character(kind=CK, len=:),                allocatable :: buffer
    type(gui_metadata_stream_pool2D_snapshot)             :: meta
    type(json_core)                                       :: json
    type(json_value),                         pointer     :: json_ptr
    type(string)                                          :: json_str, json_hash
    logical                                               :: found
    write(*,'(A)') 'test_jsonise_stream_pool2D_snapshot'
    call json%initialize(no_whitespace=.true., compact_reals=.true.)
    call meta%new(GUI_METADATA_STREAM_POOL2D_SNAPSHOT_TYPE)
    call assert_true(meta%initialized(), 'type is initialised')
    call meta%set(id=3, snapshot_filename=string('snapshot_3.simple'), snapshot_nptcls=12500)
    call assert_true(meta%assigned(), 'metadata object is set')
    json_ptr => meta%jsonise()
    call assert_true(associated(json_ptr), 'json pointer is associated')
    call json%update(json_ptr, 'snapshot_time', 0, found)
    call json%print_to_string(json_ptr, buffer)
    call assert_true(allocated(buffer), 'buffer allocated')
    call assert_true(len(buffer) > 0, 'json output is non-empty')
    json_str = buffer
    json_hash = json_str%to_fnv1a_hash64()
    call assert_char(json_hash%to_char(), '7A42F85D4D807B14', 'json is stable')
    call meta%kill()
    call assert_true(.not.meta%initialized(), 'type is not initialised')
    call json%destroy(json_ptr)
    call assert_true(.not.json%failed(), 'json destroyed')
    deallocate(buffer)
  end subroutine test_jsonise_stream_pool2D_snapshot

end module simple_gui_metadata_tester