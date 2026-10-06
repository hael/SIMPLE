!@descr: the latest GUI metadata the stream stages have sent the master, by message type, for the heartbeat
!==============================================================================
! MODULE: simple_stream_master_meta_store
!
! PURPOSE:
!   The master's copy of what each stage last reported. store() takes one
!   message (a serialised GUI metadata object, its type first) and keeps it:
!   a stage status replaces the previous one; an item of a list (a micrograph,
!   an optics group, a class average, a volume, a reprojection tile) goes to
!   slot i of a list of i_max, which is remade when i_max changes. assemble()
!   hands everything to the GUI assembler for the next heartbeat. clear_stage()
!   drops a stage's lists when the stage is restarted, so the GUI shows only
!   what the new process sends. The master holds its metadata lock around
!   all three: its listener thread stores, its main loop assembles and clears.
!
! TESTS:
!   simple_stream_master_tester
!==============================================================================
module simple_stream_master_meta_store
use simple_defs,             only: logfhandle
use simple_string_utils,     only: int2str
use simple_error,            only: simple_exception
use simple_gui_assembler,    only: gui_assembler
use simple_stream_master_stage_ids, only: STAGE_PREPROCESS, STAGE_ASSIGN_OPTICS, STAGE_INITIAL_ANALYSIS,&
    &STAGE_REFERENCE_PICKING, STAGE_PARTICLE_SIEVING, STAGE_POOL2D, STAGE_SOLVE3D
use simple_gui_metadata_micrograph,                only: gui_metadata_micrograph
use simple_gui_metadata_histogram,                 only: gui_metadata_histogram
use simple_gui_metadata_timeplot,                  only: gui_metadata_timeplot
use simple_gui_metadata_optics_group,              only: gui_metadata_optics_group
use simple_gui_metadata_cavg2D,                    only: gui_metadata_cavg2D
use simple_gui_metadata_vol3D,                     only: gui_metadata_vol3D
use simple_gui_metadata_stream_preprocess,         only: gui_metadata_stream_preprocess
use simple_gui_metadata_stream_optics_assignment,  only: gui_metadata_stream_optics_assignment
use simple_gui_metadata_stream_picking,            only: gui_metadata_stream_picking
use simple_gui_metadata_stream_initial_analysis,   only: gui_metadata_stream_initial_analysis
use simple_gui_metadata_stream_particle_sieving,   only: gui_metadata_stream_particle_sieving
use simple_gui_metadata_stream_pool2D,             only: gui_metadata_stream_pool2D
use simple_gui_metadata_stream_pool2D_snapshot,    only: gui_metadata_stream_pool2D_snapshot
use simple_gui_metadata_stream_solve3D_multistate, only: gui_metadata_stream_solve3D_multistate
use simple_gui_metadata_types, only: &
    &GUI_METADATA_STREAM_PREPROCESS_TYPE, GUI_METADATA_STREAM_PREPROCESS_HISTOGRAM_ASTIG_TYPE,&
    &GUI_METADATA_STREAM_PREPROCESS_HISTOGRAM_CTFRES_TYPE, GUI_METADATA_STREAM_PREPROCESS_HISTOGRAM_ICEFRAC_TYPE,&
    &GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_ASTIG_TYPE, GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_CTFRES_TYPE,&
    &GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_DF_TYPE, GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_RATE_TYPE,&
    &GUI_METADATA_STREAM_PREPROCESS_MICROGRAPH_TYPE, GUI_METADATA_STREAM_OPTICS_ASSIGNMENT_TYPE,&
    &GUI_METADATA_STREAM_OPTICS_ASSIGNMENT_OPTICS_GROUP_TYPE, GUI_METADATA_STREAM_INITIAL_PICKING_TYPE,&
    &GUI_METADATA_STREAM_INITIAL_PICKING_MICROGRAPH_TYPE, GUI_METADATA_STREAM_INITIAL_ANALYSIS_TYPE,&
    &GUI_METADATA_STREAM_INITIAL_ANALYSIS_VOL3D_TYPE, GUI_METADATA_STREAM_INITIAL_ANALYSIS_CLS2D_TYPE,&
    &GUI_METADATA_STREAM_INITIAL_ANALYSIS_CLS2D_FINAL_TYPE, GUI_METADATA_STREAM_REFERENCE_PICKING_TYPE,&
    &GUI_METADATA_STREAM_REFERENCE_PICKING_MICROGRAPH_TYPE, GUI_METADATA_STREAM_REFERENCE_PICKING_CLS2D_TYPE,&
    &GUI_METADATA_STREAM_PARTICLE_SIEVING_TYPE, GUI_METADATA_STREAM_PARTICLE_SIEVING_CLS2D_TYPE,&
    &GUI_METADATA_STREAM_POOL2D_TYPE, GUI_METADATA_STREAM_POOL2D_CLS2D_TYPE, GUI_METADATA_STREAM_POOL2D_SNAPSHOT_TYPE,&
    &GUI_METADATA_STREAM_POOL2D_SNAPSHOT_CLS2D_TYPE, GUI_METADATA_STREAM_SOLVE3D_MULTISTATE_TYPE,&
    &GUI_METADATA_STREAM_SOLVE3D_MULTISTATE_REPROJ_TYPE, GUI_METADATA_VOL3D_TYPE
implicit none

public :: stream_master_meta_store
private
#include "simple_local_flags.inc"

type :: stream_master_meta_store
    ! preprocessing
    type(gui_metadata_stream_preprocess)            :: preprocess
    type(gui_metadata_micrograph),      allocatable :: preprocess_micrographs(:)
    type(gui_metadata_histogram),       allocatable :: preprocess_histograms(:) ! astigmatism, CTF resolution, ice fraction
    type(gui_metadata_timeplot),        allocatable :: preprocess_timeplots(:)  ! astigmatism, CTF resolution, defocus, rate
    ! optics assignment
    type(gui_metadata_stream_optics_assignment)     :: optics_assignment
    type(gui_metadata_optics_group),    allocatable :: optics_groups(:)
    ! initial analysis: picking and initial analysis
    type(gui_metadata_stream_picking)               :: initial_picking
    type(gui_metadata_micrograph),      allocatable :: initial_picking_micrographs(:)
    type(gui_metadata_stream_initial_analysis)             :: initial_analysis
    type(gui_metadata_cavg2D),          allocatable :: initial_analysis_cavgs(:), initial_analysis_final_cavgs(:)
    type(gui_metadata_vol3D)                        :: initial_analysis_vol3D
    ! reference picking
    type(gui_metadata_stream_picking)               :: reference_picking
    type(gui_metadata_micrograph),      allocatable :: reference_picking_micrographs(:)
    type(gui_metadata_cavg2D),          allocatable :: reference_picking_cavgs(:)
    ! particle sieving; no stage sends reference class averages, which the assembler takes
    type(gui_metadata_stream_particle_sieving)      :: particle_sieving
    type(gui_metadata_cavg2D),          allocatable :: particle_sieving_cavgs(:), particle_sieving_ref_cavgs(:)
    ! pool 2D
    type(gui_metadata_stream_pool2D)                :: pool2D
    type(gui_metadata_cavg2D),          allocatable :: pool2D_cavgs(:)
    type(gui_metadata_stream_pool2D_snapshot)       :: pool2D_snapshot
    type(gui_metadata_cavg2D),          allocatable :: pool2D_snapshot_cavgs(:)
    ! multistate 3D
    type(gui_metadata_stream_solve3D_multistate) :: solve3D
    type(gui_metadata_vol3D),           allocatable :: solve3D_vols(:)
    type(gui_metadata_cavg2D),          allocatable :: solve3D_reprojtiles(:)
    logical :: l_exists = .false.
contains
    procedure :: new
    procedure :: store
    procedure :: assemble
    procedure :: clear_stage
    procedure :: log_state
    procedure :: kill
end type stream_master_meta_store

contains

    !> The stage statuses and the fixed lists, initialised; the other lists are made on receipt.
    subroutine new( self )
        class(stream_master_meta_store), intent(inout) :: self
        call self%kill()
        call self%preprocess%new(GUI_METADATA_STREAM_PREPROCESS_TYPE)
        allocate(self%preprocess_histograms(3), self%preprocess_timeplots(4))
        call self%preprocess_histograms(1)%new(GUI_METADATA_STREAM_PREPROCESS_HISTOGRAM_ASTIG_TYPE)
        call self%preprocess_histograms(2)%new(GUI_METADATA_STREAM_PREPROCESS_HISTOGRAM_CTFRES_TYPE)
        call self%preprocess_histograms(3)%new(GUI_METADATA_STREAM_PREPROCESS_HISTOGRAM_ICEFRAC_TYPE)
        call self%preprocess_timeplots(1)%new(GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_ASTIG_TYPE)
        call self%preprocess_timeplots(2)%new(GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_CTFRES_TYPE)
        call self%preprocess_timeplots(3)%new(GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_DF_TYPE)
        call self%preprocess_timeplots(4)%new(GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_RATE_TYPE)
        call self%optics_assignment%new(GUI_METADATA_STREAM_OPTICS_ASSIGNMENT_TYPE)
        call self%initial_picking%new(GUI_METADATA_STREAM_INITIAL_PICKING_TYPE)
        call self%initial_analysis%new(GUI_METADATA_STREAM_INITIAL_ANALYSIS_TYPE)
        call self%initial_analysis_vol3D%new(GUI_METADATA_STREAM_INITIAL_ANALYSIS_VOL3D_TYPE)
        call self%reference_picking%new(GUI_METADATA_STREAM_REFERENCE_PICKING_TYPE)
        call self%particle_sieving%new(GUI_METADATA_STREAM_PARTICLE_SIEVING_TYPE)
        call self%pool2D%new(GUI_METADATA_STREAM_POOL2D_TYPE)
        call self%pool2D_snapshot%new(GUI_METADATA_STREAM_POOL2D_SNAPSHOT_TYPE)
        call self%solve3D%new(GUI_METADATA_STREAM_SOLVE3D_MULTISTATE_TYPE)
        self%l_exists = .true.
    end subroutine new

    !> Keeps one message; one of an unknown type is ignored.
    subroutine store( self, buffer )
        class(stream_master_meta_store), intent(inout) :: self
        character(len=*),         intent(in)    :: buffer
        integer :: meta_type
        meta_type = transfer(buffer, meta_type)
        select case(meta_type)
            ! a stage's status, or a fixed plot
            case(GUI_METADATA_STREAM_PREPROCESS_TYPE)
                if( frame_fits(buffer, sizeof(self%preprocess), meta_type) ) self%preprocess = transfer(buffer, self%preprocess)
            case(GUI_METADATA_STREAM_PREPROCESS_HISTOGRAM_ASTIG_TYPE)
                if( frame_fits(buffer, sizeof(self%preprocess_histograms(1)), meta_type) ) self%preprocess_histograms(1) = transfer(buffer, self%preprocess_histograms(1))
            case(GUI_METADATA_STREAM_PREPROCESS_HISTOGRAM_CTFRES_TYPE)
                if( frame_fits(buffer, sizeof(self%preprocess_histograms(2)), meta_type) ) self%preprocess_histograms(2) = transfer(buffer, self%preprocess_histograms(2))
            case(GUI_METADATA_STREAM_PREPROCESS_HISTOGRAM_ICEFRAC_TYPE)
                if( frame_fits(buffer, sizeof(self%preprocess_histograms(3)), meta_type) ) self%preprocess_histograms(3) = transfer(buffer, self%preprocess_histograms(3))
            case(GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_ASTIG_TYPE)
                if( frame_fits(buffer, sizeof(self%preprocess_timeplots(1)), meta_type) ) self%preprocess_timeplots(1) = transfer(buffer, self%preprocess_timeplots(1))
            case(GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_CTFRES_TYPE)
                if( frame_fits(buffer, sizeof(self%preprocess_timeplots(2)), meta_type) ) self%preprocess_timeplots(2) = transfer(buffer, self%preprocess_timeplots(2))
            case(GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_DF_TYPE)
                if( frame_fits(buffer, sizeof(self%preprocess_timeplots(3)), meta_type) ) self%preprocess_timeplots(3) = transfer(buffer, self%preprocess_timeplots(3))
            case(GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_RATE_TYPE)
                if( frame_fits(buffer, sizeof(self%preprocess_timeplots(4)), meta_type) ) self%preprocess_timeplots(4) = transfer(buffer, self%preprocess_timeplots(4))
            case(GUI_METADATA_STREAM_OPTICS_ASSIGNMENT_TYPE)
                if( frame_fits(buffer, sizeof(self%optics_assignment), meta_type) ) self%optics_assignment = transfer(buffer, self%optics_assignment)
            case(GUI_METADATA_STREAM_INITIAL_PICKING_TYPE)
                if( frame_fits(buffer, sizeof(self%initial_picking), meta_type) ) self%initial_picking = transfer(buffer, self%initial_picking)
            case(GUI_METADATA_STREAM_INITIAL_ANALYSIS_TYPE)
                if( frame_fits(buffer, sizeof(self%initial_analysis), meta_type) ) self%initial_analysis = transfer(buffer, self%initial_analysis)
            case(GUI_METADATA_STREAM_INITIAL_ANALYSIS_VOL3D_TYPE)
                if( frame_fits(buffer, sizeof(self%initial_analysis_vol3D), meta_type) ) self%initial_analysis_vol3D = transfer(buffer, self%initial_analysis_vol3D)
            case(GUI_METADATA_STREAM_REFERENCE_PICKING_TYPE)
                if( frame_fits(buffer, sizeof(self%reference_picking), meta_type) ) self%reference_picking = transfer(buffer, self%reference_picking)
            case(GUI_METADATA_STREAM_PARTICLE_SIEVING_TYPE)
                if( frame_fits(buffer, sizeof(self%particle_sieving), meta_type) ) self%particle_sieving = transfer(buffer, self%particle_sieving)
            case(GUI_METADATA_STREAM_POOL2D_TYPE)
                if( frame_fits(buffer, sizeof(self%pool2D), meta_type) ) self%pool2D = transfer(buffer, self%pool2D)
            case(GUI_METADATA_STREAM_POOL2D_SNAPSHOT_TYPE)
                if( frame_fits(buffer, sizeof(self%pool2D_snapshot), meta_type) ) self%pool2D_snapshot = transfer(buffer, self%pool2D_snapshot)
            case(GUI_METADATA_STREAM_SOLVE3D_MULTISTATE_TYPE)
                if( frame_fits(buffer, sizeof(self%solve3D), meta_type) ) self%solve3D = transfer(buffer, self%solve3D)
            ! an item of a list
            case(GUI_METADATA_STREAM_PREPROCESS_MICROGRAPH_TYPE)
                call place_micrograph(self%preprocess_micrographs, buffer, meta_type)
            case(GUI_METADATA_STREAM_INITIAL_PICKING_MICROGRAPH_TYPE)
                call place_micrograph(self%initial_picking_micrographs, buffer, meta_type)
            case(GUI_METADATA_STREAM_REFERENCE_PICKING_MICROGRAPH_TYPE)
                call place_micrograph(self%reference_picking_micrographs, buffer, meta_type)
            case(GUI_METADATA_STREAM_OPTICS_ASSIGNMENT_OPTICS_GROUP_TYPE)
                call place_optics_group(self%optics_groups, buffer, meta_type)
            case(GUI_METADATA_STREAM_INITIAL_ANALYSIS_CLS2D_TYPE)
                call place_cavg2D(self%initial_analysis_cavgs, buffer, meta_type)
            case(GUI_METADATA_STREAM_INITIAL_ANALYSIS_CLS2D_FINAL_TYPE)
                call place_cavg2D(self%initial_analysis_final_cavgs, buffer, meta_type)
            case(GUI_METADATA_STREAM_REFERENCE_PICKING_CLS2D_TYPE)
                call place_cavg2D(self%reference_picking_cavgs, buffer, meta_type)
            case(GUI_METADATA_STREAM_PARTICLE_SIEVING_CLS2D_TYPE)
                call place_cavg2D(self%particle_sieving_cavgs, buffer, meta_type)
            case(GUI_METADATA_STREAM_POOL2D_CLS2D_TYPE)
                call place_cavg2D(self%pool2D_cavgs, buffer, meta_type)
            case(GUI_METADATA_STREAM_POOL2D_SNAPSHOT_CLS2D_TYPE)
                call place_cavg2D(self%pool2D_snapshot_cavgs, buffer, meta_type)
            case(GUI_METADATA_STREAM_SOLVE3D_MULTISTATE_REPROJ_TYPE)
                call place_cavg2D(self%solve3D_reprojtiles, buffer, meta_type)
            case(GUI_METADATA_VOL3D_TYPE)
                call place_vol3D(self%solve3D_vols, buffer, meta_type)
        end select
    end subroutine store

    !> Everything stored, to the assembler for the next heartbeat.
    subroutine assemble( self, assembler )
        class(stream_master_meta_store), intent(inout) :: self
        type(gui_assembler),      intent(inout) :: assembler
        call assembler%assemble_stream_preprocess(self%preprocess, self%preprocess_micrographs,&
            &self%preprocess_histograms, self%preprocess_timeplots)
        call assembler%assemble_stream_optics_assignment(self%optics_assignment, self%optics_groups)
        call assembler%assemble_stream_initial_picking(self%initial_picking, self%initial_picking_micrographs)
        call assembler%assemble_stream_initial_analysis(self%initial_analysis, self%initial_analysis_cavgs, self%initial_analysis_final_cavgs,&
            &self%initial_analysis_vol3D)
        call assembler%assemble_stream_reference_picking(self%reference_picking, self%reference_picking_micrographs,&
            &self%reference_picking_cavgs)
        call assembler%assemble_stream_particle_sieving(self%particle_sieving, self%particle_sieving_cavgs,&
            &self%particle_sieving_ref_cavgs)
        call assembler%assemble_stream_pool2D(self%pool2D, self%pool2D_cavgs, self%pool2D_snapshot,&
            &self%pool2D_snapshot_cavgs)
        call assembler%assemble_stream_solve3D_multistate(self%solve3D, self%solve3D_vols,&
            &self%solve3D_reprojtiles)
    end subroutine assemble

    !> Before stage @p id is started again: its lists (micrographs, optics groups, class averages,
    !! volumes, reprojection tiles) are dropped, so entries the previous process sent and the new
    !! one does not send again (a run's volumes, an older set's class averages) leave the GUI. Its
    !! status and fixed plots are replaced by the new process's first messages.
    subroutine clear_stage( self, id )
        class(stream_master_meta_store), intent(inout) :: self
        integer,                         intent(in)    :: id
        select case(id)
            case(STAGE_PREPROCESS)
                if( allocated(self%preprocess_micrographs) ) deallocate(self%preprocess_micrographs)
            case(STAGE_ASSIGN_OPTICS)
                if( allocated(self%optics_groups) ) deallocate(self%optics_groups)
            case(STAGE_INITIAL_ANALYSIS)
                if( allocated(self%initial_picking_micrographs) ) deallocate(self%initial_picking_micrographs)
                if( allocated(self%initial_analysis_cavgs)             ) deallocate(self%initial_analysis_cavgs)
                if( allocated(self%initial_analysis_final_cavgs)       ) deallocate(self%initial_analysis_final_cavgs)
            case(STAGE_REFERENCE_PICKING)
                if( allocated(self%reference_picking_micrographs) ) deallocate(self%reference_picking_micrographs)
                if( allocated(self%reference_picking_cavgs)       ) deallocate(self%reference_picking_cavgs)
            case(STAGE_PARTICLE_SIEVING)
                if( allocated(self%particle_sieving_cavgs)     ) deallocate(self%particle_sieving_cavgs)
                if( allocated(self%particle_sieving_ref_cavgs) ) deallocate(self%particle_sieving_ref_cavgs)
            case(STAGE_POOL2D)
                if( allocated(self%pool2D_cavgs)          ) deallocate(self%pool2D_cavgs)
                if( allocated(self%pool2D_snapshot_cavgs) ) deallocate(self%pool2D_snapshot_cavgs)
            case(STAGE_SOLVE3D)
                if( allocated(self%solve3D_vols)        ) deallocate(self%solve3D_vols)
                if( allocated(self%solve3D_reprojtiles) ) deallocate(self%solve3D_reprojtiles)
        end select
    end subroutine clear_stage

    !> The list sizes and @p pending_bytes (unparsed bytes in the pipes), for the memory log.
    subroutine log_state( self, pending_bytes )
        class(stream_master_meta_store), intent(in) :: self
        integer,                  intent(in) :: pending_bytes
        write(logfhandle,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)')&
            &'>>> MASTER META STATE: pending_bytes=', pending_bytes,&
            &' preprocess_mics=',   nmic(self%preprocess_micrographs),&
            &' init_pick_mics=',    nmic(self%initial_picking_micrographs),&
            &' ref_pick_mics=',     nmic(self%reference_picking_micrographs),&
            &' optics_groups=',     nog(self%optics_groups),&
            &' initial_analysis_cls=',   ncavg(self%initial_analysis_cavgs),&
            &' initial_analysis_final=', ncavg(self%initial_analysis_final_cavgs),&
            &' ref_pick_cls=',      ncavg(self%reference_picking_cavgs),&
            &' sieve_cls=',         ncavg(self%particle_sieving_cavgs),&
            &' pool2D_cls=',        ncavg(self%pool2D_cavgs),&
            &' states_vol3D=',      nvol(self%solve3D_vols),&
            &' solve3D_multistate_reprojtiles=', ncavg(self%solve3D_reprojtiles)
        call flush(logfhandle)

    contains

        integer function nmic( arr )
            type(gui_metadata_micrograph), allocatable, intent(in) :: arr(:)
            nmic = 0
            if( allocated(arr) ) nmic = size(arr)
        end function nmic

        integer function nog( arr )
            type(gui_metadata_optics_group), allocatable, intent(in) :: arr(:)
            nog = 0
            if( allocated(arr) ) nog = size(arr)
        end function nog

        integer function ncavg( arr )
            type(gui_metadata_cavg2D), allocatable, intent(in) :: arr(:)
            ncavg = 0
            if( allocated(arr) ) ncavg = size(arr)
        end function ncavg

        integer function nvol( arr )
            type(gui_metadata_vol3D), allocatable, intent(in) :: arr(:)
            nvol = 0
            if( allocated(arr) ) nvol = size(arr)
        end function nvol

    end subroutine log_state

    subroutine kill( self )
        class(stream_master_meta_store), intent(inout) :: self
        type(stream_master_meta_store) :: fresh
        if( .not. self%l_exists ) return
        ! every object back to its defaults, every list deallocated
        select type( self )
            type is( stream_master_meta_store )
                self = fresh
            class default
                THROW_HARD('kill: an extension of stream_master_meta_store needs its own kill')
        end select
    end subroutine kill

    !---------------- the lists ----------------

    ! Each of these puts the item in @p buffer into slot i of @p arr, remaking @p arr with i_max
    ! fresh objects of @p meta_type when its size is not i_max. An item whose i is out of range is
    ! reported and dropped.

    subroutine place_micrograph( arr, buffer, meta_type )
        type(gui_metadata_micrograph), allocatable, intent(inout) :: arr(:)
        character(len=*),                           intent(in)    :: buffer
        integer,                                    intent(in)    :: meta_type
        type(gui_metadata_micrograph) :: item
        integer :: i
        if( .not. frame_fits(buffer, sizeof(item), meta_type) ) return
        item = transfer(buffer, item)
        if( .not. fits(item%get_i(), item%get_i_max()) ) return
        if( allocated(arr) )then
            if( size(arr) /= item%get_i_max() ) deallocate(arr)
        endif
        if( .not. allocated(arr) )then
            allocate(arr(item%get_i_max()))
            do i = 1,size(arr)
                call arr(i)%new(meta_type)
            enddo
        endif
        arr(item%get_i()) = item
    end subroutine place_micrograph

    subroutine place_optics_group( arr, buffer, meta_type )
        type(gui_metadata_optics_group), allocatable, intent(inout) :: arr(:)
        character(len=*),                             intent(in)    :: buffer
        integer,                                      intent(in)    :: meta_type
        type(gui_metadata_optics_group) :: item
        integer :: i
        if( .not. frame_fits(buffer, sizeof(item), meta_type) ) return
        item = transfer(buffer, item)
        if( .not. fits(item%get_i(), item%get_i_max()) ) return
        if( allocated(arr) )then
            if( size(arr) /= item%get_i_max() ) deallocate(arr)
        endif
        if( .not. allocated(arr) )then
            allocate(arr(item%get_i_max()))
            do i = 1,size(arr)
                call arr(i)%new(meta_type)
            enddo
        endif
        arr(item%get_i()) = item
    end subroutine place_optics_group

    subroutine place_cavg2D( arr, buffer, meta_type )
        type(gui_metadata_cavg2D), allocatable, intent(inout) :: arr(:)
        character(len=*),                       intent(in)    :: buffer
        integer,                                intent(in)    :: meta_type
        type(gui_metadata_cavg2D) :: item
        integer :: i
        if( .not. frame_fits(buffer, sizeof(item), meta_type) ) return
        item = transfer(buffer, item)
        if( .not. fits(item%get_i(), item%get_i_max()) ) return
        if( allocated(arr) )then
            if( size(arr) /= item%get_i_max() ) deallocate(arr)
        endif
        if( .not. allocated(arr) )then
            allocate(arr(item%get_i_max()))
            do i = 1,size(arr)
                call arr(i)%new(meta_type)
            enddo
        endif
        arr(item%get_i()) = item
    end subroutine place_cavg2D

    ! The volumes' sender never includes their reprojection tiles (gui_metadata_vol3D%serialise),
    ! so the allocatable component arrives unallocated.
    subroutine place_vol3D( arr, buffer, meta_type )
        type(gui_metadata_vol3D), allocatable, intent(inout) :: arr(:)
        character(len=*),                      intent(in)    :: buffer
        integer,                               intent(in)    :: meta_type
        type(gui_metadata_vol3D) :: item
        integer :: i
        if( .not. frame_fits(buffer, sizeof(item), meta_type) ) return
        item = transfer(buffer, item)
        if( .not. fits(item%get_i(), item%get_i_max()) ) return
        if( allocated(arr) )then
            if( size(arr) /= item%get_i_max() ) deallocate(arr)
        endif
        if( .not. allocated(arr) )then
            allocate(arr(item%get_i_max()))
            do i = 1,size(arr)
                call arr(i)%new(meta_type)
            enddo
        endif
        arr(item%get_i()) = item
    end subroutine place_vol3D

    logical function fits( i, i_max )
        integer, intent(in) :: i, i_max
        fits = i_max >= 1 .and. i >= 1 .and. i <= i_max
        if( .not. fits ) THROW_WARN('a GUI list item with an index out of range is dropped')
    end function fits

    ! .true. when @p buffer is as long as the type it is copied into (@p nbytes): a frame read from
    ! a desynchronised pipe (after a resync) can carry any tag and any length, and a copy of the
    ! wrong length into a type with an allocatable component would leave it undefined.
    logical function frame_fits( buffer, nbytes, meta_type )
        use, intrinsic :: iso_c_binding, only: c_size_t
        character(len=*),  intent(in) :: buffer
        integer(c_size_t), intent(in) :: nbytes
        integer,           intent(in) :: meta_type
        character(len=:), allocatable :: warning
        frame_fits = int(len(buffer), c_size_t) == nbytes
        if( frame_fits ) return
        warning = 'a GUI message of type '//int2str(meta_type)//' with '//int2str(len(buffer))//&
            &' bytes instead of '//int2str(int(nbytes))//' is dropped'
        THROW_WARN(warning)
    end function frame_fits

end module simple_stream_master_meta_store
