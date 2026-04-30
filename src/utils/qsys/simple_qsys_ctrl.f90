!@descr: batch-processing manager - script generation, job scheduling, and persistent worker dispatch
module simple_qsys_ctrl
use simple_core_module_api
use simple_qsys_base,                only: qsys_base
use simple_qsys_slurm,               only: qsys_slurm
use simple_qsys_lsf,                 only: qsys_lsf
use simple_qsys_persistent_worker,   only: qsys_persistent_worker
use simple_cmdline,                  only: cmdline
use simple_parameters,               only: parameters
use simple_syslib,                   only: simple_rmfile
use simple_persistent_worker_message_task, only: qsys_persistent_worker_message_task 
use simple_persistent_worker_server,       only: persistent_worker
use simple_mem_estimator
implicit none

public :: qsys_ctrl
private
#include "simple_local_flags.inc"

!> Poll interval (seconds) used by blocking scheduler loops.
integer, parameter :: SHORTTIME = 1

!> Central controller that owns script generation, job-status tracking, and submission
!! for one range of partitions.  Supports local, SLURM, LSF, PBS, SGE, and the
!! persistent TCP worker backend through a polymorphic qsys_base pointer.
type qsys_ctrl
    private
    ! --- job scripts and sentinel files ---
    type(string)                  :: exec_binary                   !< executable launched by every generated script
    type(string),     allocatable :: script_names(:)               !< per-partition script file paths
    type(string),     allocatable :: jobs_done_fnames(:)           !< per-partition completion sentinel file paths
    type(string),     allocatable :: jobs_exit_code_fnames(:)      !< per-partition exit-code file paths (streaming only)
    ! --- scheduler state ---
    class(qsys_base), pointer     :: myqsys     => null()          !< polymorphic scheduler backend (not owned)
    integer,          pointer     :: parts(:,:) => null()          !< fromp/top particle ranges for all partitions (not owned)
    integer                       :: fromto_part(2)         = 0    !< [first, last] partition index managed by this instance
    integer                       :: nparts_tot             = 0    !< total number of partitions across all instances
    integer                       :: ncomputing_units       = 0    !< total computing slots available to this instance
    integer                       :: ncomputing_units_avail = 0    !< slots currently free for new submissions
    integer                       :: nthr_worker            = 1    !< thread slots required per job (persistent worker backend only)
    integer                       :: numlen                 = 0    !< zero-padding width for partition index in file names
    ! --- job tracking ---
    logical,          allocatable :: jobs_done(:)                  !< .true. once a partition job has finished
    logical,          allocatable :: jobs_submitted(:)             !< .true. once a partition job has been submitted
    ! --- streaming state ---
    class(cmdline),   allocatable :: stream_cline_stack(:)         !< pending command lines waiting for a free slot
    class(cmdline),   allocatable :: stream_cline_submitted(:)     !< command lines currently running in each slot
    class(cmdline),   allocatable :: stream_cline_done_stack(:)    !< command lines that completed successfully
    class(cmdline),   allocatable :: stream_cline_fail_stack(:)    !< command lines that exited with an error
    integer                       :: n_stream_updates       = 0    !< number of times update_queue has been called
    integer                       :: cline_stacksz          = 0    !< current number of entries in stream_cline_stack
    logical                       :: stream    = .false.           !< .true. when running in streaming (continuous) mode
    logical                       :: existence = .false.           !< .true. after new(); .false. after kill()
    contains
    ! CONSTRUCTORS
    procedure          :: new
    ! GETTERS
    procedure          :: get_jobs_status
    procedure          :: print_jobs_status
    procedure          :: exists
    ! SETTERS
    procedure          :: free_all_cunits
    procedure          :: set_jobs_status
    procedure          :: clear_stack
    ! SCRIPT GENERATORS
    procedure          :: generate_scripts
    procedure          :: generate_scripts_subprojects
    procedure          :: generate_array_script
    procedure, private :: generate_script_1, generate_script_2, generate_script_3, generate_script_4
    generic            :: generate_script => generate_script_2, generate_script_3, generate_script_4
    ! SUBMISSION TO QSYS
    procedure          :: submit_scripts
    procedure          :: submit_script
    procedure          :: dispatch_task_to_persistent_worker
    ! QUERIES
    procedure, private :: update_queue
    ! THE MASTER SCHEDULERS
    procedure          :: schedule_jobs
    procedure          :: schedule_subproject_jobs
    procedure          :: schedule_array_jobs
    ! STREAMING
    procedure          :: schedule_streaming
    procedure          :: add_to_streaming
    procedure, private :: add_to_stream_stack
    procedure          :: get_stream_done_stack
    procedure          :: get_stream_fail_stack
    procedure          :: get_stacksz, get_stack_range
    procedure          :: get_done_stacksz
    procedure          :: get_failed_stacksz
    procedure          :: get_ncomputing_units_avail
    ! DESTRUCTOR
    procedure          :: kill
end type qsys_ctrl

interface qsys_ctrl
    module procedure constructor
end interface qsys_ctrl

contains

    ! CONSTRUCTORS

    !> Functional constructor: allocates and initialises in one call.
    function constructor( exec_binary, qsys_obj, parts, fromto_part, ncomputing_units, stream ) result( self )
        class(string),            intent(in) :: exec_binary      !< executable launched by generated scripts
        class(qsys_base), target, intent(in) :: qsys_obj         !< scheduler backend (not owned)
        integer,          target, intent(in) :: parts(:,:)       !< fromp/top particle ranges (not owned)
        integer,                  intent(in) :: fromto_part(2)   !< [first, last] partition index for this instance
        integer,                  intent(in) :: ncomputing_units !< total computing slots (<= partitions controlled)
        logical,                  intent(in) :: stream           !< .true. for streaming / continuous mode
        type(qsys_ctrl) :: self
        call self%new(exec_binary, qsys_obj, parts, fromto_part, ncomputing_units, stream)
    end function constructor

    !> Initialise (or reinitialise) the controller: wire scheduler pointers, allocate
    !! tracking arrays, and pre-compute script and sentinel file names.
    subroutine new( self, exec_binary, qsys_obj, parts, fromto_part, ncomputing_units, stream, numlen, nthr_worker )
        class(qsys_ctrl),         intent(inout) :: self             !< this instance
        class(string),            intent(in)    :: exec_binary      !< executable launched by generated scripts
        class(qsys_base), target, intent(in)    :: qsys_obj         !< scheduler backend (not owned)
        integer,          target, intent(in)    :: parts(:,:)       !< fromp/top particle ranges (not owned)
        integer,                  intent(in)    :: fromto_part(2)   !< [first, last] partition index for this instance
        integer,                  intent(in)    :: ncomputing_units !< total computing slots (<= partitions controlled)
        logical,                  intent(in)    :: stream           !< .true. for streaming / continuous mode
        integer, optional,        intent(in)    :: numlen           !< explicit zero-padding width; derived from nparts_tot when absent
        integer, optional,        intent(in)    :: nthr_worker      !< thread slots per job for persistent worker backend (default 1)
        integer :: ipart
        call self%kill
        ! Wire scheduler state.
        self%stream                 =  stream
        self%exec_binary            =  exec_binary
        self%myqsys                 => qsys_obj
        self%parts                  => parts
        self%fromto_part            =  fromto_part
        self%nparts_tot             =  size(parts, 1)
        self%ncomputing_units       =  ncomputing_units
        self%ncomputing_units_avail =  ncomputing_units
        if( present(nthr_worker) ) self%nthr_worker = nthr_worker
        ! Streaming mode uses fixed 5-digit padding; batch mode derives it from nparts_tot.
        if( self%stream ) then
            self%numlen = 5
        else
            if( present(numlen) ) then
                self%numlen = numlen
            else
                self%numlen = len(int2str(self%nparts_tot))
            end if
        end if
        ! Allocate per-partition tracking arrays.
        allocate( self%jobs_done(fromto_part(1):fromto_part(2)),          &
                  self%jobs_submitted(fromto_part(1):fromto_part(2)),      &
                  self%script_names(fromto_part(1):fromto_part(2)),        &
                  self%jobs_done_fnames(fromto_part(1):fromto_part(2)),    &
                  self%jobs_exit_code_fnames(fromto_part(1):fromto_part(2)) )
        if( self%stream ) then
            ! In streaming mode jobs start as 'done' (free slot) and are filled on demand.
            self%jobs_done = .true.
            allocate(self%stream_cline_submitted(fromto_part(1):fromto_part(2)))
        else
            self%jobs_done = .false.
        end if
        self%jobs_submitted = .false.
        ! Pre-compute script file names.
        do ipart = fromto_part(1), fromto_part(2)
            self%script_names(ipart) = 'distr_simple_script_'//int2str_pad(ipart, self%numlen)
        end do
        ! Pre-compute completion sentinel and exit-code file names.
        do ipart = self%fromto_part(1), self%fromto_part(2)
            self%jobs_done_fnames(ipart)      = JOB_FINISHED_FBODY//int2str_pad(ipart, self%numlen)
            self%jobs_exit_code_fnames(ipart) = 'EXIT_CODE_JOB_'//int2str_pad(ipart, self%numlen)
        end do
        self%existence = .true.
    end subroutine new

    ! GETTERS

    !> Return copies of the jobs_done and jobs_submitted status arrays.
    subroutine get_jobs_status( self, jobs_done, jobs_submitted )
        class(qsys_ctrl),     intent(in)    :: self
        logical, allocatable, intent(inout) :: jobs_done(:)     !< copy of completion flags on return
        logical, allocatable, intent(inout) :: jobs_submitted(:) !< copy of submission flags on return
        if( allocated(jobs_done) )      deallocate(jobs_done)
        if( allocated(jobs_submitted) ) deallocate(jobs_submitted)
        allocate(jobs_done(size(self%jobs_done)),           source=self%jobs_done)
        allocate(jobs_submitted(size(self%jobs_submitted)), source=self%jobs_submitted)
    end subroutine get_jobs_status

    !> Write per-partition submission and completion flags to the log handle.
    subroutine print_jobs_status( self )
        class(qsys_ctrl), intent(in) :: self
        integer :: i
        do i = 1, size(self%jobs_submitted)
            write(logfhandle,*) i, 'submitted: ', self%jobs_submitted(i), 'done: ', self%jobs_done(i)
        end do
    end subroutine print_jobs_status

    !> Return .true. when the controller has been initialized.
    logical function exists( self )
        class(qsys_ctrl), intent(in) :: self
        exists = self%existence
    end function exists

    ! SETTERS

    !> Mark all computing slots as available (resets ncomputing_units_avail).
    subroutine free_all_cunits( self )
        class(qsys_ctrl), intent(inout) :: self
        self%ncomputing_units_avail = self%ncomputing_units
    end subroutine free_all_cunits

    !> Overwrite the first size(jobs_done) entries of the internal tracking arrays.
    !! Aborts if either input is larger than the allocated range.
    subroutine set_jobs_status( self, jobs_done, jobs_submitted )
        class(qsys_ctrl), intent(inout) :: self
        logical,          intent(in)    :: jobs_done(:)      !< new completion flags
        logical,          intent(in)    :: jobs_submitted(:) !< new submission flags
        if( size(jobs_done)      > size(self%jobs_done)      ) THROW_HARD('jobs_done array larger than allocated range; set_jobs_status')
        if( size(jobs_submitted) > size(self%jobs_submitted) ) THROW_HARD('jobs_submitted array larger than allocated range; set_jobs_status')
        self%jobs_done(:size(jobs_done))           = jobs_done
        self%jobs_submitted(:size(jobs_submitted)) = jobs_submitted
    end subroutine set_jobs_status

    !> Deallocate the pending streaming command-line stack and reset its size counter.
    subroutine clear_stack( self )
        class(qsys_ctrl), intent(inout) :: self
        if( allocated(self%stream_cline_stack) ) then
            call self%stream_cline_stack(:)%kill
            deallocate(self%stream_cline_stack)
        end if
        self%cline_stacksz = 0
    end subroutine clear_stack

    ! SCRIPT GENERATORS

    !> Generate one bash script per subproject for qsys-agnostic parallel execution.
    !! Each subproject gets its own script indexed by its position in jobs_descr(:).
    !! After this call, invoke schedule_jobs() to submit them in parallel up to
    !! ncomputing_units concurrency.  Compatible with local, SLURM, LSF, PBS, SGE.
    subroutine generate_scripts_subprojects( self, jobs_descr, q_descr, exec_bin, subproj_dirs )
        class(qsys_ctrl),        intent(inout) :: self
        type(chash),             intent(inout) :: jobs_descr(:)   !< one job description per subproject
        class(chash),            intent(in)    :: q_descr         !< queue-system metadata (scheduler directives)
        class(string), optional, intent(in)    :: exec_bin        !< override default executable
        class(string), optional, intent(in)    :: subproj_dirs(:) !< per-subproject working directories; defaults to CWD
        character(len=512) :: io_msg
        type(string)       :: job_str, execution_binary
        integer            :: isub, ios, funit, nsub
        ! nsub is defined by the input array — each element IS a subproject
        nsub = size(jobs_descr)
        if( nsub < 1 )then
            THROW_HARD('need at least 1 subproject; generate_scripts_subprojects')
        endif
        if( present(subproj_dirs) )then
            if( size(subproj_dirs) /= nsub )then
                THROW_HARD('# subproject directories must match # subprojects; generate_scripts_subprojects')
            endif
        endif
        if( present(exec_bin) )then
            execution_binary = exec_bin
        else
            execution_binary = self%exec_binary
        endif
        do isub = 1, nsub
            if( present(subproj_dirs) )then
                self%jobs_done_fnames(isub)      = trim(subproj_dirs(isub)%to_char())//'/'//SUBPROJECT_JOB_FINISHED_FBODY//int2str_pad(isub,self%numlen)
            else
                self%jobs_done_fnames(isub)      = SUBPROJECT_JOB_FINISHED_FBODY//int2str_pad(isub,self%numlen)
            endif
            call fopen(funit, file=self%script_names(isub), iostat=ios, STATUS='REPLACE', action='WRITE', iomsg=io_msg)
            call fileiochk('simple_qsys_ctrl :: generate_scripts_subprojects; Error when opening file for writing: '&
                &//self%script_names(isub)%to_char()//' ; '//trim(io_msg), ios)
            ! specify shell
            write(funit,'(a)') '#!/bin/bash'
            ! write qsys-specific instructions (run-time polymorphic)
            if( q_descr%get('qsys_name') .ne. 'local' )then
                call self%myqsys%write_instr(q_descr, fhandle=funit)
            else
                call self%myqsys%write_instr(jobs_descr(isub), fhandle=funit)
            endif
            ! change to subproject directory if provided, otherwise use global CWD
            if( present(subproj_dirs) )then
                write(funit,'(a)') 'cd '//subproj_dirs(isub)%to_char()
            else
                write(funit,'(a)') 'cd '//trim(CWD_GLOB)
            endif
            write(funit,'(a)') ''
            write(funit,'(a)') 'set -o pipefail'
            write(funit,'(a)') ''
            ! compose the command line from the subproject's job description
            call jobs_descr(isub)%set('part',   int2str(isub))
            call jobs_descr(isub)%set('numlen', int2str(self%numlen))
            job_str = jobs_descr(isub)%chash2str()
            write(funit,'(a)',advance='no') 'if '//execution_binary%to_char()//' '//job_str%to_char()
            ! direct output
            write(funit,'(a)') ' '//STDERR2STDOUT//' | tee -a '//SIMPLE_SUBPROC_OUT
            write(funit,'(a)') 'then'
            write(funit,'(a)') '  touch '//self%jobs_done_fnames(isub)%to_char()
            write(funit,'(a)') '  exit 0'
            write(funit,'(a)') 'else'
            write(funit,'(a)') '  exit 1'
            write(funit,'(a)') 'fi'
            call fclose(funit)
            call jobs_descr(isub)%delete('part')
            call jobs_descr(isub)%delete('numlen')
            if( q_descr%get('qsys_name') .eq. 'local' )then
                ios = simple_chmod(self%script_names(isub), '+x')
                if( ios /= 0 ) THROW_HARD('simple_qsys_ctrl :: generate_scripts_subprojects; Error chmoding submit script '//self%script_names(isub)%to_char())
            endif
            ! reset job tracking flags for this subproject
            self%jobs_done(isub)      = .false.
            self%jobs_submitted(isub) = .false.
            call del_file(self%jobs_done_fnames(isub))
            call wait_for_closure(self%script_names(isub))
        end do
        ! reset available computing units
        if( .not. self%stream ) self%ncomputing_units_avail = self%ncomputing_units
    end subroutine generate_scripts_subprojects

    !> Generate one script per partition in the fromto_part range via generate_script_1().
    !! Augments job_descr with fromp/top/part/nparts keys for each partition before
    !! generating the script, then deletes those keys afterwards.  Resets
    !! ncomputing_units_avail when not in streaming mode.
    subroutine generate_scripts( self, job_descr, ext, q_descr, outfile_body, part_params, extra_params )
        class(qsys_ctrl),           intent(inout) :: self
        class(chash),               intent(inout) :: job_descr         !< base job description (modified in-place, restored on return)
        class(string),              intent(in)    :: ext               !< output file extension written into outfile keys
        class(chash),               intent(inout) :: q_descr           !< queue-system metadata
        class(string),    optional, intent(in)    :: outfile_body      !< base name for per-partition output files
        class(chash),     optional, intent(in)    :: part_params(:)    !< per-partition extra key/value pairs
        type(parameters), optional, intent(in)    :: extra_params      !< extra parameters passed to memory estimator
        type(string) :: outfile_body_local, key, val
        integer      :: ipart, iadd
        logical      :: part_params_present
        if( present(outfile_body) ) outfile_body_local = outfile_body
        part_params_present = present(part_params)
        do ipart=self%fromto_part(1),self%fromto_part(2)
            call job_descr%set('fromp',  int2str(self%parts(ipart,1)))
            call job_descr%set('top',    int2str(self%parts(ipart,2)))
            call job_descr%set('part',   int2str(ipart))
            call job_descr%set('nparts', int2str(self%nparts_tot))
            if( outfile_body_local%is_allocated() )then
                call job_descr%set('outfile', outfile_body_local%to_char()//int2str_pad(ipart,self%numlen)//METADATA_EXT)
            endif
            if( part_params_present  )then
                do iadd=1,part_params(ipart)%size_of()
                    key = part_params(ipart)%get_key(iadd)
                    val = part_params(ipart)%get(iadd)
                    call job_descr%set(key%to_char(), val)
                end do
            endif
            if(L_USE_AUTO_MEM) call estimate_mem_usage(job_descr, q_descr, extra_params)
            call self%generate_script_1(job_descr, ipart, q_descr)
        end do
        call job_descr%delete('fromp')
        call job_descr%delete('top')
        call job_descr%delete('part')
        call job_descr%delete('nparts')
        if( outfile_body_local%is_allocated() )then
            call job_descr%delete('outfile')
            call outfile_body_local%kill
        endif
        if( part_params_present )then
            do iadd=1,part_params(1)%size_of()
                key = part_params(1)%get_key(iadd)
                call job_descr%delete(key%to_char())
            end do
        endif
        ! when we generate the scripts we also reset the number of available computing units
        if( .not. self%stream ) self%ncomputing_units_avail = self%ncomputing_units
    end subroutine generate_scripts

    !> Generate a single SLURM/LSF array-job script covering all partitions.
    !! Each partition command is stored as one element of a bash array; the scheduler
    !! selects the element indexed by SLURM_ARRAY_TASK_ID at runtime.
    !! Only supported by SLURM and LSF — aborts for other backends.
    subroutine generate_array_script( self, job_descr, ext, q_descr, outfile_body, part_params )
        class(qsys_ctrl),           intent(inout) :: self
        class(chash),               intent(inout) :: job_descr         !< base job description (modified in-place, restored on return)
        class(string),              intent(in)    :: ext               !< output file extension
        class(chash),               intent(in)    :: q_descr           !< queue-system metadata
        class(string),    optional, intent(in)    :: outfile_body      !< base name for per-partition output files
        class(chash),     optional, intent(in)    :: part_params(:)    !< per-partition extra key/value pairs
        type(string) :: outfile_body_local, key, val, job_str
        integer :: ipart, iadd, ios, funit
        logical :: part_params_present
        character(len=512) :: io_msg
        select type( pmyqsys => self%myqsys )
            class is(qsys_slurm)
                ! all good
            class is(qsys_lsf)
                ! all good
            class DEFAULT
                THROW_HARD('array submission only supported by SLURM, and LSF')
        end select
        if( present(outfile_body) ) outfile_body_local = outfile_body
        part_params_present = present(part_params)
        call fopen(funit, file=string(ARRAY_SCRIPT), iostat=ios, STATUS='REPLACE', action='WRITE', iomsg=io_msg)
        call fileiochk('simple_qsys_ctrl :: generate_array_script; Error when opening file for writing: '//ARRAY_SCRIPT//' ; '//trim(io_msg), ios)
        ! need to specify shell
        write(funit,'(a)') '#!/bin/bash'
        ! Write array-job scheduler directives; pass nactive when throttling below nparts_tot.
        if( self%nparts_tot /= self%ncomputing_units ) then
            call self%myqsys%write_array_instr(q_descr, self%fromto_part, fhandle=funit, nactive=self%ncomputing_units)
        else
            call self%myqsys%write_array_instr(q_descr, self%fromto_part, fhandle=funit)
        end if
        write(funit,'(a)') 'cd '//trim(CWD_GLOB)
        write(funit,'(a)') ''
        ! start partsarray definition
        write(funit,'(a)') 'partsarray=('
        do ipart=self%fromto_part(1),self%fromto_part(2)
            call job_descr%set('fromp',  int2str(self%parts(ipart,1)))
            call job_descr%set('top',    int2str(self%parts(ipart,2)))
            call job_descr%set('part',   int2str(ipart))
            call job_descr%set('nparts', int2str(self%nparts_tot))
            if( outfile_body_local%is_allocated() )then
                call job_descr%set('outfile', outfile_body_local%to_char()//int2str_pad(ipart,self%numlen)//METADATA_EXT)
            endif
            if( part_params_present  )then
                do iadd=1,part_params(ipart)%size_of()
                    key = part_params(ipart)%get_key(iadd)
                    val = part_params(ipart)%get(iadd)
                    call job_descr%set(key%to_char(), val)
                end do
            endif
            ! compose the command line as array element inside partsarray. achar(39) is apostrophe
            job_str = job_descr%chash2str()
            write(funit,'(a)',advance='no') achar(39)//self%exec_binary%to_char()//' '//job_str%to_char()
            ! direct output
            write(funit,'(a)') ' '//STDERR2STDOUT//' | tee -a '//SIMPLE_SUBPROC_OUT//achar(39)
            write(funit,'(a)') ''
        end do
        ! close partsarray definition
        write(funit,'(a)') ')'
        write(funit,'(a)') ''
        ! execute command in partsarray with index SLURM_ARRAY_TASK_ID
        write(funit,'(a)') '${partsarray[$SLURM_ARRAY_TASK_ID]}'
        write(funit,'(a)') ''
        ! exit shell when done
        write(funit,'(a)') 'exit'
        call fclose(funit)
        call job_descr%delete('fromp')
        call job_descr%delete('top')
        call job_descr%delete('part')
        call job_descr%delete('nparts')
        if( outfile_body_local%is_allocated() )then
            call job_descr%delete('outfile')
            call outfile_body_local%kill
        endif
        if( part_params_present )then
            do iadd=1,part_params(1)%size_of()
                key = part_params(1)%get_key(iadd)
                call job_descr%delete(key%to_char())
            end do
        endif
        ! when we have generated the script we also unflag jobs_submitted and jobs_done
        self%jobs_done(:)      = .false.
        self%jobs_submitted(:) = .false.
        ! and reset the number of available computing units
        if( .not. self%stream ) self%ncomputing_units_avail = self%ncomputing_units
        call wait_for_closure(string(ARRAY_SCRIPT))
    end subroutine generate_array_script

    !> Internal helper: write one bash submission script for partition ipart.
    !! Called in a loop by generate_scripts().  Resets jobs_done(ipart) and
    !! jobs_submitted(ipart) to .false. so the scheduler starts fresh.
    subroutine generate_script_1( self, job_descr, ipart, q_descr )
        class(qsys_ctrl), intent(inout) :: self
        class(chash),     intent(in)    :: job_descr !< fully populated job description for this partition
        integer,          intent(in)    :: ipart     !< partition index
        class(chash),     intent(in)    :: q_descr   !< queue-system metadata
        character(len=512) :: io_msg
        type(string) :: job_str
        integer      :: ios, funit
        call fopen(funit, file=self%script_names(ipart), iostat=ios, STATUS='REPLACE', action='WRITE', iomsg=io_msg)
        call fileiochk('simple_qsys_ctrl :: generate_script_1; Error when opening file for writing: '&
                //self%script_names(ipart)%to_char()//' ; '//trim(io_msg),ios )
        ! need to specify shell
        write(funit,'(a)') '#!/bin/bash'
        ! write (run-time polymorphic) instructions to the qsys
        if( q_descr%get('qsys_name').ne.'local' )then
            call self%myqsys%write_instr(q_descr, fhandle=funit)
        else
            call self%myqsys%write_instr(job_descr, fhandle=funit)
        endif
        write(funit,'(a)') 'cd '//trim(CWD_GLOB)
        write(funit,'(a)') ''
        ! compose the command line
        job_str = job_descr%chash2str()
        write(funit,'(a)',advance='no') self%exec_binary%to_char()//' '//job_str%to_char()
        ! direct output
        write(funit,'(a)') ' '//STDERR2STDOUT//' | tee -a '//SIMPLE_SUBPROC_OUT
        ! exit shell when done
        write(funit,'(a)') ''
        write(funit,'(a)') 'exit'
        call fclose(funit)
        if( q_descr%get('qsys_name').eq.'local' )then
            ios = simple_chmod(self%script_names(ipart),'+x')
            if( ios .ne. 0 ) THROW_HARD('simple_qsys_ctrl :: generate_script_1; Error chmoding submit script '//self%script_names(ipart)%to_char())
        endif
        ! when we generate the script we also unflag jobs_submitted and jobs_done
        self%jobs_done(ipart)      = .false.
        self%jobs_submitted(ipart) = .false.
        call wait_for_closure(self%script_names(ipart))
    end subroutine generate_script_1

    !> Generic overload generate_script: single job with explicit chash description.
    !! When outfile is present output is redirected to that path; otherwise it is
    !! appended to SIMPLE_SUBPROC_OUT.  When exit_code_fname is present the exit
    !! status is written to that file after the executable finishes.
    subroutine generate_script_2( self, job_descr, q_descr, exec_bin, script_name, outfile, exit_code_fname )
        class(qsys_ctrl),        intent(in) :: self
        class(chash),            intent(in) :: job_descr        !< job key/value parameters
        class(chash),            intent(in) :: q_descr          !< queue-system metadata
        class(string),           intent(in) :: exec_bin         !< executable to invoke
        class(string),           intent(in) :: script_name      !< path of the script to write
        class(string), optional, intent(in) :: outfile          !< redirect stdout/stderr here when present
        class(string), optional, intent(in) :: exit_code_fname  !< write $? to this file when present
        character(len=512) :: io_msg
        type(string) :: job_str
        integer      :: ios, funit
        call fopen(funit, file=script_name, iostat=ios, STATUS='REPLACE', action='WRITE', iomsg=io_msg)
        call fileiochk('simple_qsys_ctrl :: generate_script_2; Error when opening file: '//&
            &script_name%to_char()//' ; '//trim(io_msg),ios )
        ! need to specify shell
        write(funit,'(a)') '#!/bin/bash'
        ! write (run-time polymorphic) instructions to the qsys
        if( q_descr%get('qsys_name').ne.'local' )then
            call self%myqsys%write_instr(q_descr, fhandle=funit)
        else
            call self%myqsys%write_instr(job_descr, fhandle=funit)
        endif
        write(funit,'(a)') 'cd '//trim(CWD_GLOB)
        write(funit,'(a)') ''
        ! compose the command line
        job_str = job_descr%chash2str()
        write(funit,'(a)',advance='no') exec_bin%to_char()//' '//job_str%to_char()
        ! direct output
        if( present(outfile) )then
            ! unique output
            write(funit,'(a)') ' > '//outfile%to_char()//' '//STDERR2STDOUT
        else
            ! subprocess, global output
            write(funit,'(a)') ' '//STDERR2STDOUT//' | tee -a '//SIMPLE_SUBPROC_OUT
        endif
        ! exit code
        if( present(exit_code_fname) )then
            write(funit,'(a)') ''
            write(funit,'(a)') 'echo $? > '//exit_code_fname%to_char()
        endif
        ! exit shell when done
        write(funit,'(a)') ''
        write(funit,'(a)') 'exit'
        call fclose(funit)
        call wait_for_closure(script_name)
        if( q_descr%get('qsys_name').eq.'local' )then
            ios=simple_chmod(script_name,'+x')
            if( ios .ne. 0 ) THROW_HARD('simple_qsys_ctrl :: generate_script_2; Error chmoding submit script '//script_name%to_char())
        endif
    end subroutine generate_script_2

    !> Generic overload generate_script: N sequential jobs packed into one script.
    !! All jobs share one output file.  exec_bins(:), when provided, overrides
    !! exec_bin on a per-job basis and must match size(jobs_descr).
    subroutine generate_script_4( self, jobs_descr, q_descr, exec_bin, script_name, outfile, exec_bins )
        class(qsys_ctrl),          intent(in) :: self
        type(chash),  allocatable, intent(in) :: jobs_descr(:)  !< one job description per sequential task
        class(chash),              intent(in) :: q_descr         !< queue-system metadata
        class(string),             intent(in) :: exec_bin        !< default executable
        class(string),             intent(in) :: script_name     !< path of the script to write
        class(string),             intent(in) :: outfile         !< shared stdout/stderr log path
        class(string), optional,   intent(in) :: exec_bins(:)    !< per-job executable overrides
        type(string) :: execution_binary, job_str
        character(len=512) :: io_msg
        integer :: ios, funit, i, njobs
        call fopen(funit, file=script_name, iostat=ios, STATUS='REPLACE', action='WRITE', iomsg=io_msg)
        call fileiochk('simple_qsys_ctrl :: generate_script_4; Error when opening file: '//&
            &script_name%to_char()//' ; '//trim(io_msg),ios )
        ! need to specify shell
        write(funit,'(a)') '#!/bin/bash'
        ! write (run-time polymorphic) instructions to the qsys
        if( q_descr%get('qsys_name').ne.'local' )then
            call self%myqsys%write_instr(q_descr, fhandle=funit)
        else
            call self%myqsys%write_instr(jobs_descr(1), fhandle=funit)
        endif
        write(funit,'(a)') 'cd '//trim(CWD_GLOB)
        write(funit,'(a)') ''
        ! compose the command line
        njobs = size(jobs_descr)
        if( present(exec_bins) )then
            if( size(exec_bins) /= njobs )then
                THROW_HARD('# of jobs must be the same as number of executables!')
            endif
        endif
        ! When njobs == 1 the loop below is skipped; the single job is written directly after.
        if( njobs > 1 )then
            do i = 1,njobs-1
                if( present(exec_bins) )then
                    execution_binary = exec_bins(i)
                else
                    execution_binary = exec_bin
                endif
                job_str = jobs_descr(i)%chash2str()
                write(funit,'(a)',advance='no') execution_binary%to_char()//' '//job_str%to_char()
                write(funit,'(a)') ' '//STDERR2STDOUT//' | tee -a '//outfile%to_char()
                write(funit,'(a)') ''
            enddo
        endif
        if( present(exec_bins) )then
            execution_binary = exec_bins(njobs)
        else
            execution_binary = exec_bin
        endif
        job_str = jobs_descr(njobs)%chash2str()
        write(funit,'(a)',advance='no') execution_binary%to_char()//' '//job_str%to_char()
        write(funit,'(a)') ' '//STDERR2STDOUT//' | tee -a '//outfile%to_char()
        ! exit shell when done
        write(funit,'(a)') ''
        write(funit,'(a)') 'exit'
        call fclose(funit)
        call wait_for_closure(script_name)
        if( q_descr%get('qsys_name').eq.'local' )then
            ios=simple_chmod(script_name,'+x')
            if( ios .ne. 0 ) THROW_HARD('simple_qsys_ctrl :: generate_script_4; Error chmoding submit script '//script_name%to_char())
        endif
    end subroutine generate_script_4

    !> Generic overload generate_script: single job derived from a cmdline object.
    !! Converts cline to a chash job description internally before writing the script.
    subroutine generate_script_3( self, cline, q_descr, script_name, prgoutput )
        class(qsys_ctrl), intent(inout) :: self
        class(cmdline),   intent(in)    :: cline        !< command-line parameters for the job
        class(chash),     intent(in)    :: q_descr      !< queue-system metadata
        class(string),    intent(in)    :: script_name  !< path of the script to write
        class(string),    intent(in)    :: prgoutput    !< stdout/stderr log path
        type(chash)        :: job_descr
        character(len=512) :: io_msg
        integer            :: ios, funit
        type(string)       :: job_str
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! open for write
        call fopen(funit, file=script_name, iostat=ios, STATUS='REPLACE', action='WRITE', iomsg=io_msg)
        call fileiochk('simple_qsys_ctrl :: generate_script_3; Error when opening file for writing: '&
        &//script_name%to_char()//' ; '//trim(io_msg),ios )
        ! need to specify shell
        write(funit,'(a)') '#!/bin/bash'
        ! write (run-time polymorphic) instructions to the qsys
        if( q_descr%get('qsys_name').ne.'local' )then
            call self%myqsys%write_instr(q_descr, fhandle=funit)
        else
            call self%myqsys%write_instr(job_descr, fhandle=funit)
        endif
        write(funit,'(a)') 'cd '//trim(CWD_GLOB)
        write(funit,'(a)') ''
        ! compose the command line
        job_str = job_descr%chash2str()
        write(funit,'(a)', advance='no') self%exec_binary%to_char()//' '//job_str%to_char()
        ! direct output
        write(funit,'(a)') ' '//STDERR2STDOUT//' | tee -a '//prgoutput%to_char()
        ! exit shell when done
        write(funit,'(a)') ''
        write(funit,'(a)') 'exit'
        call fclose(funit)
        if( q_descr%get('qsys_name').eq.'local' )then
            ios = simple_chmod(script_name,'+x')
            if( ios .ne. 0 ) THROW_HARD('simple_qsys_ctrl :: generate_script_3; Error chmoding submit script '//script_name%to_char())
        endif
        call wait_for_closure(script_name)
        call job_descr%kill
    end subroutine generate_script_3

    ! SUBMISSION TO QSYS

    !> Submit all unsubmitted scripts in the partition range, up to ncomputing_units_avail
    !! concurrency.  For the persistent worker backend each script is enqueued via
    !! dispatch_task_to_persistent_worker; for local/SLURM/LSF it is launched directly
    !! via exec_cmdline with up to QSYS_SUBMISSION_RETRY_LIMIT retries on failure.
    subroutine submit_scripts( self )
        use simple_qsys_local,   only: qsys_local
        class(qsys_ctrl), intent(inout) :: self
        type(string) :: qsys_cmd, script_name
        integer      :: ipart, submission_exitstat, submission_retry
        logical      :: submit_or_not(self%fromto_part(1):self%fromto_part(2))
        ! Build a submission mask: mark each unsubmitted job that fits in available slots.
        submit_or_not = .false.
        do ipart = self%fromto_part(1), self%fromto_part(2)
            if( self%jobs_submitted(ipart) ) cycle   ! already in flight
            if( self%ncomputing_units_avail > 0 ) then
                submit_or_not(ipart)       = .true.
                self%jobs_submitted(ipart) = .true.
                self%ncomputing_units_avail = self%ncomputing_units_avail - 1
            end if
        end do
        if( .not. any(submit_or_not) ) return
        ! Submit each marked partition.
        do ipart = self%fromto_part(1), self%fromto_part(2)
            if( .not. submit_or_not(ipart) ) cycle
            script_name = filepath(string(PATH_HERE), self%script_names(ipart))
            if( .not. file_exists(script_name) ) &
                write(logfhandle,'(A,A)') 'FILE DOES NOT EXIST: ', script_name%to_char()
            select type( pmyqsys => self%myqsys )
                class is(qsys_local)
                    qsys_cmd = self%myqsys%submit_cmd()//' '//script_name%to_char()//' '//SUPPRESS_MSG//'&'
                class is(qsys_persistent_worker)
                    ! Worker backend: enqueue and move to next partition without blocking.
                    call self%dispatch_task_to_persistent_worker(script_name, self%nthr_worker)
                    cycle
                class DEFAULT
                    qsys_cmd = self%myqsys%submit_cmd()//' '//script_name%to_char()
            end select
            ! Attempt submission with exponential-backoff retry.
            submission_exitstat = -1
            do submission_retry = 1, QSYS_SUBMISSION_RETRY_LIMIT
                call exec_cmdline(qsys_cmd, exitstat=submission_exitstat)
                if( submission_exitstat == 0 ) exit
                write(logfhandle,'(A,A,A)') 'qsys submission failed. Retrying in ', &
                    int2str(QSYS_SUBMISSION_RETRY_SLEEP * (QSYS_SUBMISSION_RETRY_MULTI ** submission_retry)), ' seconds'
                call sleep(QSYS_SUBMISSION_RETRY_SLEEP * (QSYS_SUBMISSION_RETRY_MULTI ** submission_retry))
                write(logfhandle,'(A,I2,A,I2)') 'Retrying qsys submission', submission_retry, '/', QSYS_SUBMISSION_RETRY_LIMIT
                if( submission_retry == QSYS_SUBMISSION_RETRY_LIMIT ) THROW_HARD('qsys submission failed after multiple retries!')
            end do
        end do
    end subroutine submit_scripts

    !> Submit a single pre-generated script by name.
    !! For the persistent worker backend the script is enqueued and returns immediately.
    !! For local/SLURM/LSF the submission command is run via exec_cmdline with up to
    !! QSYS_SUBMISSION_RETRY_LIMIT retries on failure.
    subroutine submit_script( self, script_name )
        use simple_qsys_local, only: qsys_local
        class(qsys_ctrl), intent(inout) :: self
        class(string),    intent(in)    :: script_name !< path of the script to submit
        type(string) :: cmd
        integer      :: submission_exitstat, submission_retry
        if( .not. file_exists(filepath(string(PATH_HERE), script_name)) ) &
            write(logfhandle,'(A,A)') 'FILE DOES NOT EXIST: ', script_name%to_char()
        select type( pmyqsys => self%myqsys )
            class is (qsys_local)
                cmd = self%myqsys%submit_cmd()//' '//filepath(string(CWD_GLOB),script_name)//' '//SUPPRESS_MSG//'&'
            class is (qsys_persistent_worker)
                call self%dispatch_task_to_persistent_worker(filepath(string(CWD_GLOB),script_name), self%nthr_worker)
                return
            class DEFAULT
                cmd = self%myqsys%submit_cmd()//' '//filepath(string(CWD_GLOB),script_name)
        end select
        ! Attempt submission with exponential-backoff retry.
        submission_exitstat = -1
        do submission_retry = 1, QSYS_SUBMISSION_RETRY_LIMIT
            call exec_cmdline(cmd, exitstat=submission_exitstat)
            if( submission_exitstat == 0 ) return
            write(logfhandle,'(A,A,A)') 'qsys submission failed. Retrying in ', &
                int2str(QSYS_SUBMISSION_RETRY_SLEEP * (QSYS_SUBMISSION_RETRY_MULTI ** submission_retry)), ' seconds'
            call sleep(QSYS_SUBMISSION_RETRY_SLEEP * (QSYS_SUBMISSION_RETRY_MULTI ** submission_retry))
            write(logfhandle,'(A,I2,A,I2)') 'Retrying qsys submission', submission_retry, '/', QSYS_SUBMISSION_RETRY_LIMIT
        end do
        THROW_HARD('qsys submission failed after multiple retries!')
    end subroutine submit_script

    !> Enqueue a script as a normal-priority task on the shared persistent worker server.
    !! nthr is the number of thread slots required by the task.
    !! If the server is not associated the task is silently dropped with a log warning.
    subroutine dispatch_task_to_persistent_worker( self, script_path, nthr )
        class(qsys_ctrl), intent(inout) :: self
        type(string),     intent(in)    :: script_path !< absolute path of the script to enqueue
        integer,          intent(in)    :: nthr         !< thread slots required by this task
        type(qsys_persistent_worker_message_task)   :: task_msg
        if( .not. associated(persistent_worker%server) ) then
            write(logfhandle,'(A)') '>>> QSYS_CTRL dispatch_task_to_persistent_worker: worker server not initialised; task ignored'
            return
        end if
        call task_msg%new()
        task_msg%nthr        = nthr
        task_msg%script_path = script_path%to_char()
        if( .not. persistent_worker%server%queue_task(task_msg, string('norm')) ) then
            write(logfhandle,'(A,A)') '>>> QSYS_CTRL dispatch_task_to_persistent_worker: normal-priority queue full; dropped job_id ', int2str(task_msg%job_id)
        end if
        call task_msg%kill
    end subroutine dispatch_task_to_persistent_worker

    ! QUERIES

    !> Poll sentinel files to refresh job-done/submitted status and available-unit count.
    !! In streaming mode jobs are moved to the done or fail stack based on their exit
    !! code file; in batch mode the done sentinel file is sufficient.
    subroutine update_queue( self )
        class(qsys_ctrl), intent(inout) :: self
        integer :: ipart, njobs_in_queue, exit_code
        logical :: err
        if( self%stream ) then
            do ipart = self%fromto_part(1), self%fromto_part(2)
                if( self%jobs_done(ipart) ) cycle
                if( file_exists(self%jobs_exit_code_fnames(ipart)) ) then
                    ! Exit-code file present: job has completed (successfully or not).
                    call read_exit_code(self%jobs_exit_code_fnames(ipart), exit_code, err)
                    self%jobs_done(ipart) = .true.
                    if( (exit_code == 0) .and. (.not. err) .and. file_exists(self%jobs_done_fnames(ipart)) ) then
                        ! Clean exit: move cmdline to the done stack.
                        call self%add_to_stream_stack(self%stream_cline_submitted(ipart), self%stream_cline_done_stack)
                        call self%stream_cline_submitted(ipart)%delete('prg')
                    else
                        ! Error exit: move cmdline to the fail stack.
                        call self%add_to_stream_stack(self%stream_cline_submitted(ipart), self%stream_cline_fail_stack)
                    end if
                end if
            end do
            ! Available slots = number of free (done) slots, capped at ncomputing_units.
            self%ncomputing_units_avail = min(count(self%jobs_done), self%ncomputing_units)
            self%n_stream_updates = self%n_stream_updates + 1
        else
            do ipart = self%fromto_part(1), self%fromto_part(2)
                ! Sentinel file present → job finished.
                if( .not. self%jobs_done(ipart) ) self%jobs_done(ipart) = file_exists(self%jobs_done_fnames(ipart))
                ! Completed jobs must also be marked submitted (in case of external injection).
                if( self%jobs_done(ipart) ) self%jobs_submitted(ipart) = .true.
            end do
            ! Available slots = total units minus jobs that are submitted-but-not-yet-done.
            njobs_in_queue             = count(self%jobs_submitted .eqv. (.not. self%jobs_done))
            self%ncomputing_units_avail = self%ncomputing_units - njobs_in_queue
        end if
    end subroutine update_queue

    ! THE MASTER SCHEDULERS

    !> Block until all partition jobs have completed, polling at SHORTTIME intervals.
    subroutine schedule_jobs( self )
        class(qsys_ctrl), intent(inout) :: self
        do
            if( all(self%jobs_done) ) exit
            call self%update_queue
            call self%submit_scripts
            call sleep(SHORTTIME)
        end do
    end subroutine schedule_jobs

    !> Block until all subproject jobs have completed, polling at SHORTTIME intervals.
    subroutine schedule_subproject_jobs( self )
        class(qsys_ctrl), intent(inout) :: self
        do
            if( all(self%jobs_done) ) exit
            call self%update_queue
            call self%submit_scripts
            call sleep(SHORTTIME)
        end do
    end subroutine schedule_subproject_jobs

    !> Submit the single array-job script and block until all array elements complete.
    subroutine schedule_array_jobs( self )
        class(qsys_ctrl), intent(inout) :: self
        call self%submit_script(string(ARRAY_SCRIPT))
        do
            if( all(self%jobs_done) ) exit
            call self%update_queue
            call sleep(SHORTTIME)
        end do
    end subroutine schedule_array_jobs

    ! STREAMING

    !> Dispatch one round of streaming jobs: update completion status, then fill every
    !! free slot from the pending command-line stack.  When path is supplied the working
    !! directory is changed for the duration of the call and restored on return.
    subroutine schedule_streaming( self, q_descr, path )
        class(qsys_ctrl),        intent(inout) :: self
        class(chash),            intent(in)    :: q_descr          !< queue-system metadata
        type(string),  optional, intent(in)    :: path             !< working directory override
        type(chash) :: job_descr
        type(string) :: cwd, cwd_old
        integer      :: ipart
        ! Optionally change to the requested working directory.
        if( present(path) ) then
            cwd_old  = trim(CWD_GLOB)
            call simple_chdir(path)
            call simple_getcwd(cwd)
            CWD_GLOB = cwd%to_char()
        end if
        call self%update_queue
        ! Fill free slots from the front of the pending stack.
        if( self%cline_stacksz /= 0 ) then
            if( self%ncomputing_units_avail > 0 ) then
                do ipart = 1, self%ncomputing_units
                    if( self%cline_stacksz == 0 ) exit
                    if( self%jobs_done(ipart) ) then
                        ! Pop the first pending cmdline into this slot.
                        self%stream_cline_submitted(ipart) = self%stream_cline_stack(1)
                        call updatestack
                        ! Assign the computing-unit index and generate/submit the script.
                        call self%stream_cline_submitted(ipart)%set('part', ipart)
                        call self%stream_cline_submitted(ipart)%gen_job_descr(job_descr)
                        self%jobs_submitted(ipart) = .true.
                        self%jobs_done(ipart)      = .false.
                        call simple_rmfile(self%jobs_done_fnames(ipart))
                        call simple_rmfile(self%jobs_exit_code_fnames(ipart))
                        call self%generate_script_2(job_descr, q_descr, self%exec_binary, self%script_names(ipart), &
                            &exit_code_fname=self%jobs_exit_code_fnames(ipart))
                        call self%submit_script(self%script_names(ipart))
                    end if
                end do
            end if
        end if
        ! Restore working directory if changed.
        if( present(path) ) then
            call simple_chdir(cwd_old)
            CWD_GLOB = cwd_old%to_char()
        end if
        call job_descr%kill
        contains

            !> Remove the first element of stream_cline_stack and decrement cline_stacksz.
            subroutine updatestack
                class(cmdline), allocatable :: tmp_stack(:)
                if( self%cline_stacksz > 1 ) then
                    allocate(tmp_stack(self%cline_stacksz - 1), source=self%stream_cline_stack(2:self%cline_stacksz))
                    call self%stream_cline_stack(:)%kill
                    deallocate(self%stream_cline_stack)
                    self%cline_stacksz = self%cline_stacksz - 1
                    call move_alloc(tmp_stack, self%stream_cline_stack)
                else
                    deallocate(self%stream_cline_stack)
                    self%cline_stacksz = 0
                end if
            end subroutine updatestack

    end subroutine schedule_streaming

    !> Append cline to the pending streaming command-line stack (stream_cline_stack).
    subroutine add_to_streaming( self, cline )
        class(qsys_ctrl), intent(inout) :: self
        class(cmdline),   intent(in)    :: cline !< command line to enqueue
        class(cmdline), allocatable :: tmp_stack(:)
        integer :: i
        if( .not. allocated(self%stream_cline_stack) ) then
            ! First entry: allocate a single-element stack.
            allocate(self%stream_cline_stack(1), source=cline)
            self%cline_stacksz = 1
        else
            ! Grow the stack by one and append.
            call move_alloc(self%stream_cline_stack, tmp_stack)
            self%cline_stacksz = self%cline_stacksz + 1
            allocate(self%stream_cline_stack(self%cline_stacksz))
            do i = 1, self%cline_stacksz - 1
                self%stream_cline_stack(i) = tmp_stack(i)
            end do
            self%stream_cline_stack(self%cline_stacksz) = cline
            call tmp_stack(:)%kill
            deallocate(tmp_stack)
        end if
    end subroutine add_to_streaming

    !> Append cline to an arbitrary allocatable cmdline stack (used for done/fail stacks).
    !! No-op when cline has no 'prg' key set.
    subroutine add_to_stream_stack( self, cline, stack )
        class(qsys_ctrl),            intent(inout) :: self
        class(cmdline),              intent(in)    :: cline   !< command line to append
        class(cmdline), allocatable, intent(inout) :: stack(:) !< target stack
        class(cmdline), allocatable :: tmp_stack(:)
        integer :: stacksz, i
        if( .not. cline%defined('prg') ) return  ! guard: only track jobs that have a program
        if( .not. allocated(stack) ) then
            allocate(stack(1), source=cline)
        else
            stacksz = size(stack)
            call move_alloc(stack, tmp_stack)
            stacksz = stacksz + 1
            allocate(stack(stacksz))
            do i = 1, stacksz - 1
                stack(i) = tmp_stack(i)
            end do
            stack(stacksz) = cline
            call tmp_stack(:)%kill
            deallocate(tmp_stack)
        end if
    end subroutine add_to_stream_stack

    !> Transfer the completed-job stack to clines and clear it; no-op when empty.
    subroutine get_stream_done_stack( self, clines )
        class(qsys_ctrl),            intent(inout) :: self
        class(cmdline), allocatable, intent(out)   :: clines(:) !< receives ownership of the done stack
        if( .not. allocated(self%stream_cline_done_stack) ) return
        call move_alloc(self%stream_cline_done_stack, clines)
    end subroutine get_stream_done_stack

    !> Transfer the failed-job stack to clines, set n = size; no-op when empty.
    subroutine get_stream_fail_stack( self, clines, n )
        class(qsys_ctrl),            intent(inout) :: self
        class(cmdline), allocatable, intent(out)   :: clines(:) !< receives ownership of the fail stack
        integer,                     intent(out)   :: n         !< number of failed jobs transferred
        n = 0
        if( .not. allocated(self%stream_cline_fail_stack) ) return
        call move_alloc(self%stream_cline_fail_stack, clines)
        n = size(clines)
    end subroutine get_stream_fail_stack

    !> Return the number of command lines currently in the pending streaming stack.
    integer function get_stacksz( self )
        class(qsys_ctrl), intent(in) :: self
        get_stacksz = self%cline_stacksz
    end function get_stacksz

    !> Return the total particle range covered by the pending streaming stack.
    !! Computed as sum(top - fromp + 1) across all pending command lines.
    integer function get_stack_range( self )
        class(qsys_ctrl), intent(in) :: self
        if( self%cline_stacksz == 0 ) then
            get_stack_range = 0
        else
            get_stack_range =   sum(self%stream_cline_stack(:)%get_iarg('top'))   &
                              - sum(self%stream_cline_stack(:)%get_iarg('fromp')) &
                              + self%cline_stacksz
        end if
    end function get_stack_range

    !> Return the number of successfully completed streaming jobs accumulated so far.
    integer function get_done_stacksz( self )
        class(qsys_ctrl), intent(in) :: self
        if( .not. allocated(self%stream_cline_done_stack) ) then
            get_done_stacksz = 0
        else
            get_done_stacksz = size(self%stream_cline_done_stack)
        end if
    end function get_done_stacksz

    !> Return the number of failed streaming jobs accumulated so far.
    integer function get_failed_stacksz( self )
        class(qsys_ctrl), intent(in) :: self
        if( .not. allocated(self%stream_cline_fail_stack) ) then
            get_failed_stacksz = 0
        else
            get_failed_stacksz = size(self%stream_cline_fail_stack)
        end if
    end function get_failed_stacksz

    !> Return the number of computing units currently available for new submissions.
    integer function get_ncomputing_units_avail( self )
        class(qsys_ctrl), intent(in) :: self
        get_ncomputing_units_avail = self%ncomputing_units_avail
    end function get_ncomputing_units_avail

    ! DESTRUCTOR

    !> Release all allocated state and null scheduler pointers; resets existence to .false.
    subroutine kill( self )
        class(qsys_ctrl), intent(inout) :: self
        if( .not. self%existence ) return
        ! Null non-owned pointers first.
        self%exec_binary            = ''
        self%myqsys                 => null()
        self%parts                  => null()
        ! Reset scalar state.
        self%fromto_part(2)         = 0
        self%ncomputing_units       = 0
        self%ncomputing_units_avail = 0
        self%nthr_worker            = 1
        self%numlen                 = 0
        self%cline_stacksz          = 0
        ! Deallocate tracking arrays.
        deallocate(self%script_names, self%jobs_done, self%jobs_done_fnames, &
                   self%jobs_exit_code_fnames, self%jobs_submitted)
        ! Deallocate streaming stacks (kill cmdlines before dealloc where required).
        if( allocated(self%stream_cline_stack) )      deallocate(self%stream_cline_stack)
        if( allocated(self%stream_cline_submitted) ) then
            call self%stream_cline_submitted(:)%kill
            deallocate(self%stream_cline_submitted)
        end if
        if( allocated(self%stream_cline_done_stack) ) then
            call self%stream_cline_done_stack(:)%kill
            deallocate(self%stream_cline_done_stack)
        end if
        if( allocated(self%stream_cline_fail_stack) ) then
            call self%stream_cline_fail_stack(:)%kill
            deallocate(self%stream_cline_fail_stack)
        end if
        self%existence = .false.
    end subroutine kill

end module simple_qsys_ctrl
