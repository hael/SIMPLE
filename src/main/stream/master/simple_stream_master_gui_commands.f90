!@descr: the commands the GUI sends the stream master in its heartbeat response: stops, restarts and the stage updates
!==============================================================================
! MODULE: simple_stream_master_gui_commands
!
! PURPOSE:
!   Each heartbeat the master posts to the GUI is answered with a JSON object
!   that may ask to stop the whole stream ('terminate'), to stop or restart
!   one stage (terminate_<key>, restart_<key>, keys from simple_stream_master_stage_ids),
!   and carry updates for the stages (thresholds, picking references, the 2D
!   mask diameter, a 2D snapshot). parse() reads one response into this
!   record; the master applies it. The update holds this response's fields
!   only, so a field is sent to the stages once, in the update that brought
!   it. A selection larger than the update holds, a 3D snapshot of no state or
!   of a state out of range, and a snapshot whose name is not a bare *.simple
!   file name, are dropped with a warning; nothing in a response stops the
!   master.
!==============================================================================
module simple_stream_master_gui_commands
use json_kinds,                        only: CK
use json_module,                       only: json_core, json_value
use simple_defs,                       only: logfhandle, dp
use simple_defs_fname,                 only: METADATA_EXT
use simple_string,                     only: string
use simple_gui_metadata_types,         only: GUI_METADATA_STREAM_UPDATE_TYPE
use simple_gui_metadata_stream_update, only: gui_metadata_stream_update, MAX_PICKREFS_SELECTION, MAX_SNAPSHOT2D_SELECTION,&
                                            &MAX_SNAPSHOT_FNAME_LEN
use simple_gui_metadata_stream_solve3D_multistate, only: MAX_STATES_SOLVE3D_MULTISTATE
use simple_stream_master_stage_ids,    only: NSTAGES, stage_gui_key
implicit none

public :: stream_master_gui_commands
private

type :: stream_master_gui_commands
    logical                          :: l_terminate_all       = .false. ! stop the stream
    logical                          :: l_terminate(NSTAGES)  = .false. ! stop a stage
    logical                          :: l_restart(NSTAGES)    = .false. ! restart a stopped stage
    type(gui_metadata_stream_update) :: update                         ! this response's updates
contains
    procedure :: parse
    procedure :: kill
end type stream_master_gui_commands

contains

    !> Reads one heartbeat response; .false. (and nothing asked) when it is not valid JSON.
    function parse( self, text ) result( l_ok )
        class(stream_master_gui_commands), intent(inout) :: self
        character(len=*),           intent(in)    :: text
        type(json_core)                        :: json
        type(json_value), pointer              :: root, snapshot
        character(kind=CK, len=:), allocatable :: str_val
        integer,                   allocatable :: i_arr(:)
        real(kind=dp) :: r_val
        integer       :: i_val, id, snapshot_id
        logical       :: l_ok, l_found, l_id, l_iter, l_sel, l_file
        call self%kill()
        call self%update%new(GUI_METADATA_STREAM_UPDATE_TYPE)
        nullify(root, snapshot)
        l_ok = .false.
        call json%parse(root, text)
        if( json%failed() )then
            write(logfhandle,'(A,A)') 'FAILED TO PARSE JSON RESPONSE ', text
            call json%clear_exceptions()
            if( associated(root) ) call json%destroy(root)
            return
        endif
        ! stops and restarts
        self%l_terminate_all = flag('terminate')
        do id = 1,NSTAGES
            self%l_terminate(id) = flag('terminate_'//stage_gui_key(id))
            self%l_restart(id)   = flag('restart_'//stage_gui_key(id))
        enddo
        ! updates for the stages
        call json%get(root, 'ctfresthreshold', r_val, l_found)
        if( l_found ) call self%update%set_ctfres_update(real(r_val))
        call json%get(root, 'astigthreshold', r_val, l_found)
        if( l_found ) call self%update%set_astigmatism_update(real(r_val))
        call json%get(root, 'icefracthreshold', r_val, l_found)
        if( l_found ) call self%update%set_icescore_update(real(r_val))
        ! a selection too large for the update is dropped whole: acting on part of it would act on
        ! classes the user did not choose
        call json%get(root, 'pickrefs_selection', i_arr, l_found)
        if( l_found .and. allocated(i_arr) )then
            if( size(i_arr) <= MAX_PICKREFS_SELECTION )then
                call self%update%set_pickrefs_selection(i_arr)
            else
                call warn_dropped('picking reference selection', size(i_arr), MAX_PICKREFS_SELECTION)
            endif
        endif
        call json%get(root, 'pickrefs_cycle', i_val, l_found)
        if( l_found ) call self%update%set_pickrefs_cycle(i_val)
        call json%get(root, 'mskdiam2D', r_val, l_found)
        if( l_found ) call self%update%set_mskdiam2D_update(real(r_val))
        call json%get(root, 'snapshot2D', snapshot, l_found)
        if( l_found .and. associated(snapshot) )then
            if( allocated(i_arr) ) deallocate(i_arr)
            call json%get(snapshot, 'id',        snapshot_id, l_id)
            call json%get(snapshot, 'iteration', i_val,       l_iter)
            call json%get(snapshot, 'selection', i_arr,       l_sel)
            call json%get(snapshot, 'filename',  str_val,     l_file)
            if( l_id .and. l_iter .and. l_sel .and. l_file .and. allocated(i_arr) .and. allocated(str_val) )then
                if( size(i_arr) > MAX_SNAPSHOT2D_SELECTION )then
                    call warn_dropped('2D snapshot selection', size(i_arr), MAX_SNAPSHOT2D_SELECTION)
                else if( .not. snapshot_fname_ok(str_val) )then
                    ! p06 writes the snapshot into snapshots/<name without .simple>/<name>
                    write(logfhandle,'(A,A,A)') '>>> WARNING: GUI 2D snapshot name ', str_val(1:min(len(str_val),64)),&
                        &' is not a bare *.simple file name; ignored'
                else
                    call self%update%set_snapshot2D_update(snapshot_id, i_val, i_arr, string(str_val))
                endif
            endif
        endif
        ! a 3D snapshot: the particles of the selected states of multistate 3D, merged into one
        call json%get(root, 'snapshot3D', snapshot, l_found)
        if( l_found .and. associated(snapshot) )then
            if( allocated(i_arr)   ) deallocate(i_arr)
            if( allocated(str_val) ) deallocate(str_val)
            call json%get(snapshot, 'id',        snapshot_id, l_id)
            call json%get(snapshot, 'selection', i_arr,       l_sel)
            call json%get(snapshot, 'filename',  str_val,     l_file)
            if( l_id .and. l_sel .and. l_file .and. allocated(i_arr) .and. allocated(str_val) )then
                if( size(i_arr) == 0 .or. size(i_arr) > MAX_STATES_SOLVE3D_MULTISTATE )then
                    write(logfhandle,'(A,I0,A,I0,A)') '>>> WARNING: GUI 3D snapshot selection of ', size(i_arr),&
                        &' states is empty or exceeds ', MAX_STATES_SOLVE3D_MULTISTATE, '; ignored'
                else if( any(i_arr < 1) .or. any(i_arr > MAX_STATES_SOLVE3D_MULTISTATE) )then
                    write(logfhandle,'(A,I0,A)') '>>> WARNING: GUI 3D snapshot selection names a state outside 1..',&
                        &MAX_STATES_SOLVE3D_MULTISTATE, '; ignored'
                else if( .not. snapshot_fname_ok(str_val) )then
                    write(logfhandle,'(A,A,A)') '>>> WARNING: GUI 3D snapshot name ', str_val(1:min(len(str_val),64)),&
                        &' is not a bare *.simple file name; ignored'
                else
                    call self%update%set_snapshot3D_update(snapshot_id, i_arr, string(str_val))
                endif
            endif
        endif
        call json%destroy(root)
        l_ok = .true.

    contains

        logical function flag( key )
            character(len=*), intent(in) :: key
            logical :: l_val, l_key
            call json%get(root, key, l_val, l_key)
            flag = l_key .and. l_val
        end function flag

        subroutine warn_dropped( what, n, nmax )
            character(len=*), intent(in) :: what
            integer,          intent(in) :: n, nmax
            write(logfhandle,'(A,A,A,I0,A,I0,A)') '>>> WARNING: GUI ', what, ' of ', n, ' classes exceeds ', nmax,&
                &'; ignored'
        end subroutine warn_dropped

        ! A bare file name ending in METADATA_EXT, of at most MAX_SNAPSHOT_FNAME_LEN characters
        logical function snapshot_fname_ok( fname )
            character(len=*), intent(in) :: fname
            integer :: n, next
            n    = len_trim(fname)
            next = len(METADATA_EXT)
            snapshot_fname_ok = .false.
            if( n <= next .or. n > MAX_SNAPSHOT_FNAME_LEN ) return
            if( index(fname(1:n), '/') > 0 ) return
            snapshot_fname_ok = fname(n-next+1:n) == METADATA_EXT
        end function snapshot_fname_ok

    end function parse

    !> Nothing asked. The update gets its default fields back (its kill and new reset only the
    !! flags), so no earlier response's field is sent again.
    subroutine kill( self )
        class(stream_master_gui_commands), intent(inout) :: self
        type(gui_metadata_stream_update) :: fresh
        self%l_terminate_all = .false.
        self%l_terminate     = .false.
        self%l_restart       = .false.
        self%update          = fresh
    end subroutine kill

end module simple_stream_master_gui_commands
