! batch-processing manager - control module
module simple_qsys_ctrl
include 'simple_lib.f08'
use simple_qsys_base,  only: qsys_base
use simple_qsys_slurm, only: qsys_slurm
use simple_cmdline,    only: cmdline
implicit none

public :: qsys_ctrl
private
#include "simple_local_flags.inc"

integer, parameter :: SHORTTIME = 1

type qsys_ctrl
    private
    character(len=STDLEN)          :: exec_binary = ''              !< binary to execute in parallel
    character(len=32), allocatable :: script_names(:)               !< file names of generated scripts
    character(len=32), allocatable :: jobs_done_fnames(:)           !< touch files indicating completion
    character(len=32), allocatable :: jobs_exit_code_fnames(:)      !< touch files containing exit codes
    class(qsys_base),  pointer     :: myqsys     => null()          !< pointer to polymorphic qsys object
    integer, pointer               :: parts(:,:) => null()          !< defines the fromp/top ranges for all partitions
    class(cmdline),    allocatable :: stream_cline_stack(:)         !< stack of command lines, for streaming only
    class(cmdline),    allocatable :: stream_cline_submitted(:)     !< stack of submitted command lines, for streaming only
    class(cmdline),    allocatable :: stream_cline_done_stack(:)    !< stack of completed command lines, for streaming only
    class(cmdline),    allocatable :: stream_cline_fail_stack(:)    !< stack of failed command lines, for streaming only
    integer(kind=8),   allocatable :: stream_submit_time(:)         !< records submission time
    logical,           allocatable :: jobs_done(:)                  !< to indicate completion of distributed scripts
    logical,           allocatable :: jobs_submitted(:)             !< to indicate which jobs have been submitted
    integer                        :: fromto_part(2)         = 0    !< defines the range of partitions controlled by this instance
    integer                        :: nparts_tot             = 0    !< total number of partitions
    integer                        :: ncomputing_units       = 0    !< number of computing units
    integer                        :: ncomputing_units_avail = 0    !< number of available units
    integer                        :: numlen                 = 0    !< length of padded number string
    integer                        :: n_stream_updates       = 0    !< counter, for streaming only
    integer                        :: cline_stacksz          = 0    !< size of stack of command lines, for streaming only
    logical                        :: stream    = .false.           !< stream flag
    logical                        :: existence = .false.           !< indicates existence
    contains
    ! CONSTRUCTORS
    procedure          :: new
    ! GETTERS
    procedure          :: get_exec_bin
    procedure          :: get_jobs_status
    procedure          :: print_jobs_status
    procedure          :: exists
    ! SETTERS
    procedure          :: free_all_cunits
    procedure          :: set_jobs_status
    ! SCRIPT GENERATORS
    procedure          :: generate_scripts
    procedure          :: generate_array_script
    procedure, private :: generate_script_1
    procedure, private :: generate_script_2
    procedure, private :: generate_script_3
    generic            :: generate_script => generate_script_2, generate_script_3
    ! SUBMISSION TO QSYS
    procedure          :: submit_scripts
    procedure          :: submit_script
    ! QUERIES
    procedure, private :: update_queue
    ! THE MASTER SCHEDULERS
    procedure          :: schedule_jobs
    procedure          :: schedule_array_jobs
    ! STREAMING
    procedure          :: schedule_streaming
    procedure          :: add_to_streaming
    procedure, private :: add_to_stream_stack
    procedure          :: get_stream_done_stack
    procedure          :: get_stream_fail_stack
    procedure          :: get_stacksz
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

    function constructor( exec_binary, qsys_obj, parts, fromto_part, ncomputing_units, stream ) result( self )
        character(len=*),         intent(in) :: exec_binary      !< the binary that we want to execute in parallel
        class(qsys_base), target, intent(in) :: qsys_obj         !< the object that defines the qeueuing system
        integer,          target, intent(in) :: parts(:,:)       !< defines the start_ptcl/stop_ptcl ranges
        integer,                  intent(in) :: fromto_part(2)   !< defines the range of partitions controlled by this object
        integer,                  intent(in) :: ncomputing_units !< number of computing units (<= the number of parts controlled)
        logical,                  intent(in) :: stream           !< stream flag
        type(qsys_ctrl) :: self
        call self%new(exec_binary, qsys_obj, parts, fromto_part, ncomputing_units, stream )
    end function constructor

    subroutine new( self, exec_binary, qsys_obj, parts, fromto_part, ncomputing_units, stream, numlen )
        class(qsys_ctrl),         intent(inout) :: self             !< the instance
        character(len=*),         intent(in)    :: exec_binary      !< the binary that we want to execute in parallel
        class(qsys_base), target, intent(in)    :: qsys_obj         !< the object that defines the qeueuing system
        integer,          target, intent(in)    :: parts(:,:)       !< defines the start_ptcl/stop_ptcl ranges
        integer,                  intent(in)    :: fromto_part(2)   !< defines the range of partitions controlled by this object
        integer,                  intent(in)    :: ncomputing_units !< number of computing units (<= the number of parts controlled)
        logical,                  intent(in)    :: stream           !< stream flag
        integer, optional,        intent(in)    :: numlen           !< length of number string
        integer :: ipart
        call self%kill
        self%stream                 =  stream
        self%exec_binary            =  exec_binary
        self%myqsys                 => qsys_obj
        self%parts                  => parts
        self%fromto_part            =  fromto_part
        self%nparts_tot             =  size(parts,1)
        self%ncomputing_units       =  ncomputing_units
        self%ncomputing_units_avail =  ncomputing_units
        if( self%stream )then
            self%numlen = 5
        else
            if( present(numlen) )then
                self%numlen = numlen
            else
                self%numlen = len(int2str(self%nparts_tot))
            endif
        endif
        ! allocate
        allocate(   self%jobs_done(fromto_part(1):fromto_part(2)),&
                    self%jobs_submitted(fromto_part(1):fromto_part(2)),&
                    self%script_names(fromto_part(1):fromto_part(2)),&
                    self%jobs_done_fnames(fromto_part(1):fromto_part(2)),&
                    self%jobs_exit_code_fnames(fromto_part(1):fromto_part(2)))
        if( self%stream )then
            self%jobs_done  = .true.
            allocate(self%stream_cline_submitted(fromto_part(1):fromto_part(2)))
            allocate(self%stream_submit_time(fromto_part(1):fromto_part(2)),source=int(0,kind=8))
        else
            self%jobs_done  = .false.
        endif
        self%jobs_submitted = .false.
        ! create script names
        do ipart=fromto_part(1),fromto_part(2)
            self%script_names(ipart) = 'distr_simple_script_'//int2str_pad(ipart,self%numlen)
        end do
        ! create jobs done flags
        do ipart=self%fromto_part(1),self%fromto_part(2)
            self%jobs_done_fnames(ipart) = trim(JOB_FINISHED_FBODY)//int2str_pad(ipart,self%numlen)
            self%jobs_exit_code_fnames(ipart) = 'EXIT_CODE_JOB_'//int2str_pad(ipart,self%numlen)
        end do
        self%existence = .true.
    end subroutine new

    ! GETTERS

    function get_exec_bin( self ) result( exec_bin )
        class(qsys_ctrl), intent(in) :: self
        character(len=STDLEN)        :: exec_bin
        exec_bin = self%exec_binary
    end function get_exec_bin

    subroutine get_jobs_status( self, jobs_done, jobs_submitted )
        class(qsys_ctrl), intent(in) :: self
        logical, allocatable :: jobs_done(:), jobs_submitted(:)
        if( allocated(jobs_done) )      deallocate(jobs_done)
        if( allocated(jobs_submitted) ) deallocate(jobs_submitted)
        allocate(jobs_done(size(self%jobs_done)), source=self%jobs_done)
        allocate(jobs_submitted(size(self%jobs_submitted)), source=self%jobs_submitted)
    end subroutine get_jobs_status

    subroutine print_jobs_status( self )
        class(qsys_ctrl), intent(in) :: self
        integer :: i
        do i=1,size(self%jobs_submitted)
            write(logfhandle,*) i, 'submitted: ', self%jobs_submitted(i), 'done: ', self%jobs_done(i)
        end do
    end subroutine print_jobs_status

    logical function exists( self )
        class(qsys_ctrl), intent(in) :: self
        exists = self%existence
    end function exists

    ! SETTERS

    subroutine free_all_cunits( self )
        class(qsys_ctrl), intent(inout) :: self
        self%ncomputing_units_avail = self%ncomputing_units
    end subroutine free_all_cunits

    subroutine set_jobs_status( self, jobs_done, jobs_submitted )
        class(qsys_ctrl), intent(inout) :: self
        logical,          intent(in)    :: jobs_done(:), jobs_submitted(:)
        self%jobs_done(:size(jobs_done)) = jobs_done
        self%jobs_submitted(:size(jobs_submitted)) = jobs_submitted
    end subroutine set_jobs_status

    ! SCRIPT GENERATORS

    subroutine generate_scripts( self, job_descr, ext, q_descr, outfile_body, part_params )
        class(qsys_ctrl),           intent(inout) :: self
        class(chash),               intent(inout) :: job_descr
        character(len=4),           intent(in)    :: ext
        class(chash),               intent(in)    :: q_descr
        character(len=*), optional, intent(in)    :: outfile_body
        class(chash),     optional, intent(in)    :: part_params(:)
        character(len=:), allocatable :: outfile_body_local, key, val
        integer :: ipart, iadd
        logical :: part_params_present
        if( present(outfile_body) ) allocate(outfile_body_local, source=trim(outfile_body))
        part_params_present = present(part_params)
        do ipart=self%fromto_part(1),self%fromto_part(2)
            call job_descr%set('fromp',  int2str(self%parts(ipart,1)))
            call job_descr%set('top',    int2str(self%parts(ipart,2)))
            call job_descr%set('part',   int2str(ipart))
            call job_descr%set('nparts', int2str(self%nparts_tot))
            if( allocated(outfile_body_local) )then
                call job_descr%set('outfile', trim(outfile_body_local)//int2str_pad(ipart,self%numlen)//trim(METADATA_EXT))
            endif
            if( part_params_present  )then
                do iadd=1,part_params(ipart)%size_of()
                    key = part_params(ipart)%get_key(iadd)
                    val = part_params(ipart)%get(iadd)
                    call job_descr%set(key, val)
                end do
            endif
            call self%generate_script_1(job_descr, ipart, q_descr)
        end do
        call job_descr%delete('fromp')
        call job_descr%delete('top')
        call job_descr%delete('part')
        call job_descr%delete('nparts')
        if( allocated(outfile_body_local) )then
            call job_descr%delete('outfile')
            deallocate(outfile_body_local)
        endif
        if( part_params_present )then
            do iadd=1,part_params(1)%size_of()
                key = part_params(1)%get_key(iadd)
                call job_descr%delete(key)
            end do
        endif
        ! when we generate the scripts we also reset the number of available computing units
        if( .not. self%stream ) self%ncomputing_units_avail = self%ncomputing_units
    end subroutine generate_scripts

    subroutine generate_array_script( self, job_descr, ext, q_descr, outfile_body, part_params )
        class(qsys_ctrl),           intent(inout) :: self
        class(chash),               intent(inout) :: job_descr
        character(len=4),           intent(in)    :: ext
        class(chash),               intent(in)    :: q_descr
        character(len=*), optional, intent(in)    :: outfile_body
        class(chash),     optional, intent(in)    :: part_params(:)
        character(len=:), allocatable :: outfile_body_local, key, val
        integer :: ipart, iadd, ios, funit
        logical :: part_params_present
        character(len=512) :: io_msg
        select type( pmyqsys => self%myqsys )
            class is(qsys_slurm)
                ! all good
            class DEFAULT
                THROW_HARD('array submission only supported by SLURM')
        end select
        if( present(outfile_body) ) allocate(outfile_body_local, source=trim(outfile_body))
        part_params_present = present(part_params)
        call fopen(funit, file=ARRAY_SCRIPT, iostat=ios, STATUS='REPLACE', action='WRITE', iomsg=io_msg)
        call fileiochk('simple_qsys_ctrl :: generate_array_script; Error when opening file for writing: '//ARRAY_SCRIPT//' ; '//trim(io_msg), ios)
        ! need to specify shell
        write(funit,'(a)') '#!/bin/bash'
        ! write instructions to SLURM
        if( self%nparts_tot /= self%ncomputing_units )then
            call self%myqsys%write_array_instr(q_descr, self%fromto_part, fhandle=funit, nactive=self%ncomputing_units)
        else
            call self%myqsys%write_array_instr(q_descr, self%fromto_part, fhandle=funit)
        endif
        write(funit,'(a)') 'cd '//trim(cwd_glob)
        write(funit,'(a)') ''
        ! start partsarray definition
        write(funit,'(a)') 'partsarray=('
        do ipart=self%fromto_part(1),self%fromto_part(2)
            call job_descr%set('fromp',  int2str(self%parts(ipart,1)))
            call job_descr%set('top',    int2str(self%parts(ipart,2)))
            call job_descr%set('part',   int2str(ipart))
            call job_descr%set('nparts', int2str(self%nparts_tot))
            if( allocated(outfile_body_local) )then
                call job_descr%set('outfile', trim(outfile_body_local)//int2str_pad(ipart,self%numlen)//trim(METADATA_EXT))
            endif
            if( part_params_present  )then
                do iadd=1,part_params(ipart)%size_of()
                    key = part_params(ipart)%get_key(iadd)
                    val = part_params(ipart)%get(iadd)
                    call job_descr%set(key, val)
                end do
            endif
            ! compose the command line as array element inside partsarray. achar(39) is apostrophe
            write(funit,'(a)',advance='no') achar(39)//trim(self%exec_binary)//' '//trim(job_descr%chash2str())
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
        if( allocated(outfile_body_local) )then
            call job_descr%delete('outfile')
            deallocate(outfile_body_local)
        endif
        if( part_params_present )then
            do iadd=1,part_params(1)%size_of()
                key = part_params(1)%get_key(iadd)
                call job_descr%delete(key)
            end do
        endif
        ! when we have generated the script we also unflag jobs_submitted and jobs_done
        self%jobs_done(:)      = .false.
        self%jobs_submitted(:) = .false.
        ! and reset the number of available computing units
        if( .not. self%stream ) self%ncomputing_units_avail = self%ncomputing_units
        call wait_for_closure(ARRAY_SCRIPT)
    end subroutine generate_array_script

    subroutine generate_script_1( self, job_descr, ipart, q_descr )
        class(qsys_ctrl), intent(inout) :: self
        class(chash),     intent(in)    :: job_descr
        integer,          intent(in)    :: ipart
        class(chash),     intent(in)    :: q_descr
        character(len=512) :: io_msg
        integer :: ios, funit
        call fopen(funit, file=self%script_names(ipart), iostat=ios, STATUS='REPLACE', action='WRITE', iomsg=io_msg)
        call fileiochk('simple_qsys_ctrl :: gen_qsys_script; Error when opening file for writing: '&
                //trim(self%script_names(ipart))//' ; '//trim(io_msg),ios )
        ! need to specify shell
        write(funit,'(a)') '#!/bin/bash'
        ! write (run-time polymorphic) instructions to the qsys
        if( q_descr%get('qsys_name').ne.'local' )then
            call self%myqsys%write_instr(q_descr, fhandle=funit)
        else
            call self%myqsys%write_instr(job_descr, fhandle=funit)
        endif
        write(funit,'(a)') 'cd '//trim(cwd_glob)
        write(funit,'(a)') ''
        ! compose the command line
        write(funit,'(a)',advance='no') trim(self%exec_binary)//' '//trim(job_descr%chash2str())
        ! direct output
        write(funit,'(a)') ' '//STDERR2STDOUT//' | tee -a '//SIMPLE_SUBPROC_OUT
        ! exit shell when done
        write(funit,'(a)') ''
        write(funit,'(a)') 'exit'
        call fclose(funit)
        if( q_descr%get('qsys_name').eq.'local' )then
            ios = simple_chmod(trim(self%script_names(ipart)),'+x')
            if( ios .ne. 0 )then
                write(logfhandle,'(a)',advance='no') 'simple_qsys_scripts :: gen_qsys_script; Error'
                write(logfhandle,'(a)') 'chmoding submit script'//trim(self%script_names(ipart))
                stop
            endif
        endif
        ! when we generate the script we also unflag jobs_submitted and jobs_done
        self%jobs_done(ipart)      = .false.
        self%jobs_submitted(ipart) = .false.
        call wait_for_closure(self%script_names(ipart))
    end subroutine generate_script_1

    !>  \brief  public script generator for single jobs
    subroutine generate_script_2( self, job_descr, q_descr, exec_bin, script_name, outfile, job_descr2, exit_code_fname )
        class(qsys_ctrl),           intent(inout) :: self
        class(chash),               intent(in)    :: job_descr
        class(chash),               intent(in)    :: q_descr
        character(len=*),           intent(in)    :: exec_bin, script_name
        character(len=*), optional, intent(in)    :: outfile
        class(chash),     optional, intent(in)    :: job_descr2
        character(len=*), optional, intent(in)    :: exit_code_fname
        character(len=512) :: io_msg
        integer :: ios, funit
        call fopen(funit, file=script_name, iostat=ios, STATUS='REPLACE', action='WRITE', iomsg=io_msg)
        call fileiochk('simple_qsys_ctrl :: generate_script_2; Error when opening file: '//&
            &trim(script_name)//' ; '//trim(io_msg),ios )
        ! need to specify shell
        write(funit,'(a)') '#!/bin/bash'
        ! write (run-time polymorphic) instructions to the qsys
        if( q_descr%get('qsys_name').ne.'local' )then
            call self%myqsys%write_instr(q_descr, fhandle=funit)
        else
            call self%myqsys%write_instr(job_descr, fhandle=funit)
        endif
        write(funit,'(a)') 'cd '//trim(cwd_glob)
        write(funit,'(a)') ''
        ! compose the command line
        if( present(job_descr2) )then
            write(funit,'(a)',advance='no') trim(exec_bin)//' '//trim(job_descr%chash2str())
            write(funit,'(a)') ' '//STDERR2STDOUT//' | tee -a '//SIMPLE_SUBPROC_OUT
            write(funit,'(a)') ''
            write(funit,'(a)',advance='no') trim(exec_bin)//' '//trim(job_descr2%chash2str())
        else
            write(funit,'(a)',advance='no') trim(exec_bin)//' '//trim(job_descr%chash2str())
        endif
        ! direct output
        if( present(outfile) )then
            ! unique output
            write(funit,'(a)') ' > '//trim(outfile)//' '//STDERR2STDOUT
        else
            ! subprocess, global output
            write(funit,'(a)') ' '//STDERR2STDOUT//' | tee -a '//SIMPLE_SUBPROC_OUT
        endif
        ! exit code
        if( present(exit_code_fname) )then
            write(funit,'(a)') ''
            write(funit,'(a)') 'echo $? > '//trim(exit_code_fname)
        endif
        ! exit shell when done
        write(funit,'(a)') ''
        write(funit,'(a)') 'exit'
        call fclose(funit)
        call wait_for_closure(script_name)
        if( trim(q_descr%get('qsys_name')).eq.'local' )then
            ios=simple_chmod(trim(script_name),'+x')
            if( ios .ne. 0 )then
                write(logfhandle,'(a)',advance='no') 'simple_qsys_ctrl :: generate_script_2; Error'
                write(logfhandle,'(a)') 'chmoding submit script'//trim(script_name)
                stop
            end if
        endif
    end subroutine generate_script_2

    subroutine generate_script_3( self, cline, q_descr, script_name, prgoutput )
        class(qsys_ctrl), intent(inout) :: self
        class(cmdline),   intent(in)    :: cline
        class(chash),     intent(in)    :: q_descr
        character(len=*), intent(in)    :: script_name, prgoutput
        type(chash)        :: job_descr
        character(len=512) :: io_msg
        integer            :: ios, funit
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! open for write
        call fopen(funit, file=trim(script_name), iostat=ios, STATUS='REPLACE', action='WRITE', iomsg=io_msg)
        call fileiochk('simple_qsys_ctrl :: gen_qsys_script; Error when opening file for writing: '&
        &//trim(script_name)//' ; '//trim(io_msg),ios )
        ! need to specify shell
        write(funit,'(a)') '#!/bin/bash'
        ! write (run-time polymorphic) instructions to the qsys
        if( q_descr%get('qsys_name').ne.'local' )then
            call self%myqsys%write_instr(q_descr, fhandle=funit)
        else
            call self%myqsys%write_instr(job_descr, fhandle=funit)
        endif
        write(funit,'(a)') 'cd '//trim(cwd_glob)
        write(funit,'(a)') ''
        ! compose the command line
        write(funit,'(a)', advance='no') trim(self%exec_binary)//' '//trim(job_descr%chash2str())
        ! direct output
        write(funit,'(a)') ' '//STDERR2STDOUT//' | tee -a '//trim(prgoutput)
        ! exit shell when done
        write(funit,'(a)') ''
        write(funit,'(a)') 'exit'
        call fclose(funit)
        if( q_descr%get('qsys_name').eq.'local' )then
            ios = simple_chmod(trim(script_name),'+x')
            if( ios .ne. 0 )then
                write(logfhandle,'(a)',advance='no') 'simple_qsys_scripts :: gen_qsys_script; Error'
                write(logfhandle,'(a)') 'chmoding submit script'//trim(script_name)
                stop
            endif
        endif
        call wait_for_closure(trim(script_name))
        call job_descr%kill
    end subroutine generate_script_3

    ! SUBMISSION TO QSYS

    subroutine submit_scripts( self )
        use simple_qsys_local,   only: qsys_local
        class(qsys_ctrl), intent(inout) :: self
        character(len=LONGSTRLEN) :: qsys_cmd
        character(len=STDLEN)     :: script_name
        integer :: ipart, submission_exitstat, submission_retry
        logical :: submit_or_not(self%fromto_part(1):self%fromto_part(2))
        ! make a submission mask
        submit_or_not = .false.
        do ipart=self%fromto_part(1),self%fromto_part(2)
            if( self%jobs_submitted(ipart) )then
                ! do nothing
                cycle
            endif
            if( self%ncomputing_units_avail > 0 )then
                submit_or_not(ipart) = .true.
                ! flag job submitted
                self%jobs_submitted(ipart) = .true.
                ! decrement ncomputing_units_avail
                self%ncomputing_units_avail = self%ncomputing_units_avail - 1
            endif
        end do
        if( .not. any(submit_or_not) ) return
        ! on the fly submission
        do ipart=self%fromto_part(1),self%fromto_part(2)
            if( submit_or_not(ipart) )then
                script_name = filepath(PATH_HERE, trim(adjustl(self%script_names(ipart))))
                if( .not.file_exists(trim(script_name)))then
                    write(logfhandle,'(A,A)')'FILE DOES NOT EXIST:',trim(script_name)
                endif
                select type( pmyqsys => self%myqsys )
                    class is(qsys_local)
                        qsys_cmd = trim(adjustl(self%myqsys%submit_cmd()))//' '//trim(script_name)//' '//SUPPRESS_MSG//'&'
                    class DEFAULT
                        qsys_cmd = trim(adjustl(self%myqsys%submit_cmd()))//' '//trim(script_name)
                end select
                ! attempt to execute the command up to QSYS_SUBMISSION_RETRY_LIMIT
                submission_exitstat = -1
                do submission_retry = 1, QSYS_SUBMISSION_RETRY_LIMIT
                    call exec_cmdline(trim(adjustl(qsys_cmd)), exitstat=submission_exitstat)
                    if(submission_exitstat == 0) then
                        exit
                    else
                        write(logfhandle,'(A,I2,A,I2)')'qsys submission failed. Retrying ', submission_retry, '/', QSYS_SUBMISSION_RETRY_LIMIT
                        call sleep(QSYS_SUBMISSION_RETRY_SLEEP)
                    end if
                    if(submission_retry == QSYS_SUBMISSION_RETRY_LIMIT) THROW_HARD('qsys submission failed after multiple retries!')
                end do
            endif
        end do
    end subroutine submit_scripts

    subroutine submit_script( self, script_name )
        use simple_qsys_local,   only: qsys_local
        class(qsys_ctrl), intent(inout) :: self
        character(len=*), intent(in)    :: script_name
        character(len=STDLEN) :: cmd
        integer :: submission_exitstat, submission_retry
        if( .not.file_exists(filepath(PATH_HERE,trim(script_name))))then
            write(logfhandle,'(A,A)')'FILE DOES NOT EXIST:',trim(script_name)
        endif
        select type( pmyqsys => self%myqsys )
            type is (qsys_local)
                cmd = trim(adjustl(self%myqsys%submit_cmd()))//' '//&
                    &filepath(trim(cwd_glob),trim(script_name))//' '//&
                    &SUPPRESS_MSG//'&'
            class DEFAULT
                cmd = trim(adjustl(self%myqsys%submit_cmd()))//' '//&
                    &filepath(trim(cwd_glob),trim(script_name))
        end select
        ! attempt to execute the command up to QSYS_SUBMISSION_RETRY_LIMIT
        submission_exitstat = -1
        do submission_retry = 1, QSYS_SUBMISSION_RETRY_LIMIT
            call exec_cmdline(trim(cmd), exitstat=submission_exitstat)
            if(submission_exitstat == 0) return    
            write(logfhandle,'(A,I2,A,I2)')'qsys submission failed. Retrying ', submission_retry, '/', QSYS_SUBMISSION_RETRY_LIMIT
            call sleep(QSYS_SUBMISSION_RETRY_SLEEP)
        end do
        THROW_HARD('qsys submission failed after multiple retries!')
    end subroutine submit_script

    ! QUERIES

    subroutine update_queue( self )
        class(qsys_ctrl),  intent(inout) :: self
        integer :: ipart, njobs_in_queue, exit_code
        logical :: err
        if( self%stream )then
            do ipart=self%fromto_part(1),self%fromto_part(2)
                if( .not.self%jobs_done(ipart) )then
                    if( file_exists(self%jobs_exit_code_fnames(ipart)) )then
                        ! job has finished
                        call read_exit_code(self%jobs_exit_code_fnames(ipart), exit_code, err)
                        if( (exit_code == 0) .and. (.not.err) .and. file_exists(self%jobs_done_fnames(ipart)))then
                            ! gracefully
                            self%jobs_done(ipart) = .true.
                            call self%add_to_stream_stack( self%stream_cline_submitted(ipart), self%stream_cline_done_stack )
                            call self%stream_cline_submitted(ipart)%delete('prg')
                        else
                            ! some error occured, resource is being freed
                            self%jobs_done(ipart) = .true.
                            call self%add_to_stream_stack( self%stream_cline_submitted(ipart), self%stream_cline_fail_stack )
                        endif
                    endif
                endif
            end do
            self%ncomputing_units_avail = min(count(self%jobs_done), self%ncomputing_units)
            self%n_stream_updates = self%n_stream_updates + 1
        else
            do ipart=self%fromto_part(1),self%fromto_part(2)
                if( .not. self%jobs_done(ipart) ) self%jobs_done(ipart) = file_exists(self%jobs_done_fnames(ipart))
                if( self%jobs_done(ipart) ) self%jobs_submitted(ipart) = .true.
            end do
            njobs_in_queue = count(self%jobs_submitted .eqv. (.not.self%jobs_done))
            self%ncomputing_units_avail = self%ncomputing_units - njobs_in_queue
        endif
    end subroutine update_queue

    ! THE MASTER SCHEDULERS

    subroutine schedule_jobs( self )
        class(qsys_ctrl), intent(inout) :: self
        do
            if( all(self%jobs_done) ) exit
            call self%update_queue
            call self%submit_scripts
            call sleep(SHORTTIME)
        end do
    end subroutine schedule_jobs

    subroutine schedule_array_jobs( self )
        class(qsys_ctrl), intent(inout) :: self
        call self%submit_script(ARRAY_SCRIPT)
        do
            if( all(self%jobs_done) ) exit
            call self%update_queue
            call sleep(SHORTTIME)
        end do
    end subroutine schedule_array_jobs

    ! STREAMING

    subroutine schedule_streaming( self, q_descr, path )
        class(qsys_ctrl),            intent(inout) :: self
        class(chash),                intent(in)    :: q_descr
        character(len=*),  optional, intent(in)    :: path
        type(cmdline)              :: cline
        type(chash)                :: job_descr
        character(len=XLONGSTRLEN) :: cwd, cwd_old
        integer                    :: ipart
        if( present(path) )then
            cwd_old = trim(cwd_glob)
            call chdir(path)
            call simple_getcwd(cwd)
            cwd_glob = trim(cwd)
        endif
        call self%update_queue
        if( self%cline_stacksz /= 0 )then
            if( self%ncomputing_units_avail > 0 )then
                do ipart = 1, self%ncomputing_units
                    if( self%cline_stacksz .eq. 0 )exit
                    if( self%jobs_done(ipart) )then
                        cline = self%stream_cline_stack(1)
                        call updatestack
                        call cline%set('part', real(ipart))        ! computing unit allocation
                        self%stream_cline_submitted(ipart) = cline ! stash
                        call cline%gen_job_descr(job_descr)
                        self%jobs_submitted(ipart) = .true.
                        self%jobs_done(ipart)      = .false.
                        call del_file(self%jobs_done_fnames(ipart))
                        call del_file(self%jobs_exit_code_fnames(ipart))
                        call self%generate_script_2(job_descr, q_descr, self%exec_binary, self%script_names(ipart),&
                        &exit_code_fname=self%jobs_exit_code_fnames(ipart) )
                        call self%submit_script( self%script_names(ipart) )
                        self%stream_submit_time(ipart) = time8()
                    endif
                enddo
            endif
        endif
        if( present(path) )then
            call chdir(cwd_old)
            cwd_glob = trim(cwd_old)
        endif
        contains

            ! move everything up so '1' is index for next job to be run
            subroutine updatestack
                class(cmdline), allocatable :: tmp_stack(:)
                integer :: i
                if( self%cline_stacksz > 1 )then
                    allocate(tmp_stack(self%cline_stacksz-1))
                    call self%stream_cline_stack(1)%kill
                    do i = 1, self%cline_stacksz-1
                        tmp_stack(i) = self%stream_cline_stack(i+1)
                        call self%stream_cline_stack(i+1)%kill
                    enddo
                    deallocate(self%stream_cline_stack)
                    self%cline_stacksz = self%cline_stacksz-1
                    call move_alloc(tmp_stack,self%stream_cline_stack)
                else
                    deallocate(self%stream_cline_stack)
                    self%cline_stacksz = 0
                endif
            end subroutine updatestack

    end subroutine schedule_streaming

    !>  \brief  is to append to the streaming command-line stack
    subroutine add_to_streaming( self, cline )
        class(qsys_ctrl), intent(inout) :: self
        class(cmdline),   intent(in)    :: cline
        class(cmdline), allocatable :: tmp_stack(:)
        integer :: i
        if( .not. allocated(self%stream_cline_stack) )then
            ! empty stack
            allocate( self%stream_cline_stack(1))
            self%stream_cline_stack(1) = cline
            self%cline_stacksz = 1
        else
            ! append
            call move_alloc(self%stream_cline_stack, tmp_stack)
            self%cline_stacksz = self%cline_stacksz + 1
            allocate( self%stream_cline_stack(self%cline_stacksz) )
            do i = 1, self%cline_stacksz-1
                self%stream_cline_stack(i) = tmp_stack(i)
                call tmp_stack(i)%kill
            enddo
            self%stream_cline_stack(self%cline_stacksz) = cline
            deallocate(tmp_stack)
        endif
    end subroutine add_to_streaming

    !>  \brief  is to append to the streaming command-line stack
    subroutine add_to_stream_stack( self, cline, stack )
        class(qsys_ctrl),            intent(inout) :: self
        class(cmdline),              intent(in)    :: cline
        class(cmdline), allocatable, intent(inout) :: stack(:)
        class(cmdline), allocatable :: tmp_stack(:)
        integer :: i, stacksz
        if( .not.cline%defined('prg') )return
        if( .not. allocated(stack) )then
            ! empty stack
            allocate( stack(1), source=cline)
        else
            ! append
            stacksz = size(stack)
            call move_alloc(stack, tmp_stack)
            stacksz = stacksz + 1
            allocate( stack(stacksz) )
            do i = 1, stacksz-1
                stack(i) = tmp_stack(i)
                call tmp_stack(i)%kill
            enddo
            stack(stacksz) = cline
            deallocate(tmp_stack)
        endif
    end subroutine add_to_stream_stack


    !>  \brief  is to get the streaming command-line stack of DONE jobs, which also deallocates it
    subroutine get_stream_done_stack( self, clines )
        class(qsys_ctrl),            intent(inout) :: self
        class(cmdline), allocatable, intent(out)   :: clines(:)
        if( .not. allocated(self%stream_cline_done_stack) )return
        call move_alloc(self%stream_cline_done_stack, clines)
    end subroutine get_stream_done_stack

    !>  \brief  is to get the streaming command-line stack of DONE jobs, which also deallocates it
    subroutine get_stream_fail_stack( self, clines, n )
        class(qsys_ctrl),            intent(inout) :: self
        class(cmdline), allocatable, intent(out)   :: clines(:)
        integer,                     intent(out)   :: n
        n = 0
        if( .not. allocated(self%stream_cline_fail_stack) )return
        call move_alloc(self%stream_cline_fail_stack, clines)
        n = size(clines)
    end subroutine get_stream_fail_stack

    !>  \brief  returns streaming jobs stack size
    integer function get_stacksz( self )
        class(qsys_ctrl), intent(in) :: self
        get_stacksz = self%cline_stacksz
    end function get_stacksz

    !>  \brief  returns streaming jobs done stack size
    integer function get_done_stacksz( self )
        class(qsys_ctrl), intent(in) :: self
        if( .not. allocated(self%stream_cline_done_stack) )then
            get_done_stacksz = 0
        else
            get_done_stacksz = size(self%stream_cline_done_stack)
        endif
    end function get_done_stacksz

    !>  \brief  returns streaming jobs done stack size
    integer function get_failed_stacksz( self )
        class(qsys_ctrl), intent(in) :: self
        if( .not. allocated(self%stream_cline_fail_stack) )then
            get_failed_stacksz = 0
        else
            get_failed_stacksz = size(self%stream_cline_fail_stack)
        endif
    end function get_failed_stacksz

    integer function get_ncomputing_units_avail( self )
        class(qsys_ctrl), intent(in) :: self
        get_ncomputing_units_avail = count(self%jobs_submitted)
    end function get_ncomputing_units_avail

    ! DESTRUCTOR

    !>  \brief  is a destructor
    subroutine kill( self )
        class(qsys_ctrl), intent(inout) :: self
        if( self%existence )then
            self%exec_binary            =  ''
            self%myqsys                 => null()
            self%parts                  => null()
            self%fromto_part(2)         =  0
            self%ncomputing_units       =  0
            self%ncomputing_units_avail =  0
            self%numlen                 =  0
            self%cline_stacksz          =  0
            deallocate(self%script_names, self%jobs_done, self%jobs_done_fnames,&
            &self%jobs_exit_code_fnames, self%jobs_submitted)
            if(allocated(self%stream_cline_stack))     deallocate(self%stream_cline_stack)
            if(allocated(self%stream_cline_submitted)) deallocate(self%stream_cline_submitted)
            if(allocated(self%stream_cline_done_stack))deallocate(self%stream_cline_done_stack)
            if(allocated(self%stream_cline_fail_stack))deallocate(self%stream_cline_fail_stack)
            if(allocated(self%stream_submit_time))     deallocate(self%stream_submit_time)
            self%existence = .false.
        endif
    end subroutine kill

end module simple_qsys_ctrl
