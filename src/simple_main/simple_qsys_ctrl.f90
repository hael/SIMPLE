!------------------------------------------------------------------------------!
! SIMPLE v2.5         Elmlund & Elmlund Lab          simplecryoem.com          !
!------------------------------------------------------------------------------!
!> Simple batch-processing manager - control module
module simple_qsys_ctrl
use simple_defs
use simple_qsys_base, only: qsys_base
use simple_chash,     only: chash
use simple_strings,   only: int2str, int2str_pad
implicit none

public :: qsys_ctrl
private
#include "simple_local_flags.inc"
integer, parameter :: SHORTTIME = 3

type qsys_ctrl
    private
    character(len=STDLEN)          :: exec_binary      = ''      !< binary to execute in parallel
                                                                 !< trim(simplepath)//'/bin/simple_exec'
    character(len=32), allocatable :: script_names(:)            !< file names of generated scripts
    character(len=32), allocatable :: jobs_done_fnames(:)        !< touch files indicating completion
    character(len=STDLEN)          :: pwd              = ''      !< working directory
    class(qsys_base), pointer      :: myqsys           => null() !< pointer to polymorphic qsys object
    integer, pointer               :: parts(:,:)       => null() !< defines the fromp/top ranges for all partitions
    logical, allocatable           :: jobs_done(:)               !< to indicate completion of distributed scripts
    logical, allocatable           :: jobs_submitted(:)          !< to indicate which jobs have been submitted
    integer                        :: fromto_part(2)         = 0 !< defines the range of partitions controlled by this instance
    integer                        :: nparts_tot             = 0 !< total number of partitions
    integer                        :: ncomputing_units       = 0 !< number of computing units
    integer                        :: ncomputing_units_avail = 0 !< number of available units
    integer                        :: numlen = 0                 !< length of padded number string
    logical                        :: stream = .false.           !< stream flag
    logical                        :: existence = .false.        !< indicates existence
  contains
    ! CONSTRUCTOR
    procedure          :: new
    ! GETTER
    procedure          :: get_exec_bin
    procedure          :: get_jobs_status
    procedure          :: print_jobs_status
    procedure          :: exists
    ! SETTER
    procedure          :: free_all_cunits
    procedure          :: set_jobs_status
    ! SCRIPT GENERATORS
    procedure          :: generate_scripts
    procedure, private :: generate_script_1
    procedure, private :: generate_script_2
    generic            :: generate_script => generate_script_2
    ! SUBMISSION TO QSYS
    procedure          :: submit_scripts
    procedure          :: submit_script
    ! QUERIES
    procedure          :: update_queue
    ! THE MASTER SCHEDULER
    procedure          :: schedule_jobs
    ! DESTRUCTOR
    procedure          :: kill
end type qsys_ctrl

interface qsys_ctrl
    module procedure constructor
end interface qsys_ctrl

contains

    ! CONSTRUCTORS

    !>  \brief  is a constructor
    function constructor( exec_binary, qsys_obj, parts, fromto_part, ncomputing_units, stream ) result( self )
        character(len=*),          intent(in)    :: exec_binary      !< the binary that we want to execute in parallel
        class(qsys_base),  target, intent(in)    :: qsys_obj         !< the object that defines the qeueuing system
        integer,           target, intent(in)    :: parts(:,:)       !< defines the start_ptcl/stop_ptcl ranges
        integer,                   intent(in)    :: fromto_part(2)   !< defines the range of partitions controlled by this object
        integer,                   intent(in)    :: ncomputing_units !< number of computing units (<= the number of parts controlled)
        logical,                   intent(in)    :: stream           !< stream flag
        type(qsys_ctrl) :: self
        call self%new(exec_binary, qsys_obj, parts, fromto_part, ncomputing_units, stream )
    end function constructor

    !>  \brief  is a constructor
    subroutine new( self, exec_binary, qsys_obj, parts, fromto_part, ncomputing_units, stream )
        use simple_jiffys, only: alloc_err
        class(qsys_ctrl),          intent(inout) :: self             !< the instance
        character(len=*),          intent(in)    :: exec_binary      !< the binary that we want to execute in parallel
        class(qsys_base),  target, intent(in)    :: qsys_obj         !< the object that defines the qeueuing system
        integer,           target, intent(in)    :: parts(:,:)       !< defines the start_ptcl/stop_ptcl ranges
        integer,                   intent(in)    :: fromto_part(2)   !< defines the range of partitions controlled by this object
        integer,                   intent(in)    :: ncomputing_units !< number of computing units (<= the number of parts controlled)
        logical,                   intent(in)    :: stream           !< stream flag
        integer :: alloc_stat, ipart
        call self%kill
        self%stream                 = stream
        self%exec_binary            =  exec_binary
        self%myqsys                 => qsys_obj
        self%parts                  => parts
        self%fromto_part            =  fromto_part
        self%nparts_tot             =  size(parts,1)
        self%ncomputing_units       =  ncomputing_units
        self%ncomputing_units_avail =  ncomputing_units
        if( stream )then
            self%numlen = 5
        else
            self%numlen =  len(int2str(self%nparts_tot))
        endif
        ! allocate
        allocate(   self%jobs_done(fromto_part(1):fromto_part(2)),&
                    self%jobs_submitted(fromto_part(1):fromto_part(2)),&
                    self%script_names(fromto_part(1):fromto_part(2)),&
                    self%jobs_done_fnames(fromto_part(1):fromto_part(2)), stat=alloc_stat)
        call alloc_err("In: simple_qsys_ctrl :: new", alloc_stat)
        self%jobs_done      = .false.
        self%jobs_submitted = .false.
        ! create script names
        do ipart=fromto_part(1),fromto_part(2)
            self%script_names(ipart) = 'distr_simple_script_'//int2str_pad(ipart,self%numlen)
        end do
        ! create jobs done flags
        do ipart=self%fromto_part(1),self%fromto_part(2)
            self%jobs_done_fnames(ipart) = 'JOB_FINISHED_'//int2str_pad(ipart,self%numlen)
        end do
        ! get pwd
        call get_environment_variable('PWD', self%pwd)
        self%existence = .true.
    end subroutine new

    ! GETTER

    !>  \brief  for getting the execute simple binary
    function get_exec_bin( self ) result( exec_bin )
        class(qsys_ctrl), intent(in) :: self
        character(len=STDLEN)        :: exec_bin
        exec_bin = self%exec_binary
    end function get_exec_bin

    !>  \brief  for getting jobs status logical flags
    subroutine get_jobs_status( self, jobs_done, jobs_submitted )
        class(qsys_ctrl), intent(in) :: self
        logical, allocatable :: jobs_done(:), jobs_submitted(:)
        if( allocated(jobs_done) )      deallocate(jobs_done)
        if( allocated(jobs_submitted) ) deallocate(jobs_submitted)
        allocate(jobs_done(size(self%jobs_done)), source=self%jobs_done)
        allocate(jobs_submitted(size(self%jobs_submitted)), source=self%jobs_submitted)
    end subroutine get_jobs_status

    !>  \brief  for printing jobs status logical flags
    subroutine print_jobs_status( self )
        class(qsys_ctrl), intent(in) :: self
        integer :: i
        do i=1,size(self%jobs_submitted)
            print *, i, 'submitted: ', self%jobs_submitted(i), 'done: ', self%jobs_done(i)
        end do
    end subroutine print_jobs_status

    !>  \brief  for checking existence
    logical function exists( self )
        class(qsys_ctrl), intent(in) :: self
        exists = self%existence
    end function exists

    ! SETTER

    !> \brief for freeing all available computing units
    subroutine free_all_cunits( self )
        class(qsys_ctrl), intent(inout) :: self
        self%ncomputing_units_avail = self%ncomputing_units
    end subroutine free_all_cunits

    !>  \brief  for setting jobs status logical flags
    subroutine set_jobs_status( self, jobs_done, jobs_submitted )
        class(qsys_ctrl), intent(inout) :: self
        logical,          intent(in)    :: jobs_done(:), jobs_submitted(:)
        self%jobs_done(:size(jobs_done)) = jobs_done
        self%jobs_submitted(:size(jobs_submitted)) = jobs_submitted
    end subroutine set_jobs_status

    ! SCRIPT GENERATORS

    !>  \brief  public script generator
    subroutine generate_scripts( self, job_descr, ext, q_descr, outfile_body, part_params, chunkdistr )
        class(qsys_ctrl),           intent(inout) :: self
        class(chash),               intent(inout) :: job_descr
        character(len=4),           intent(in)    :: ext
        class(chash),               intent(in)    :: q_descr
        character(len=*), optional, intent(in)    :: outfile_body
        class(chash),     optional, intent(in)    :: part_params(:)
        logical,          optional, intent(in)    :: chunkdistr
        character(len=:), allocatable :: outfile_body_local, key, val
        integer :: ipart, iadd
        logical :: part_params_present, cchunkdistr
        cchunkdistr = .false.
        if( present(chunkdistr) ) cchunkdistr = .true.
        if( present(outfile_body) )then
            allocate(outfile_body_local, source=trim(outfile_body))
        endif
        part_params_present = present(part_params)
        do ipart=self%fromto_part(1),self%fromto_part(2)
            if( cchunkdistr )then
                ! don't input part-dependent parameters
            else
                call job_descr%set('fromp',   int2str(self%parts(ipart,1)))
                call job_descr%set('top',     int2str(self%parts(ipart,2)))
                call job_descr%set('part',    int2str(ipart))
                call job_descr%set('nparts',  int2str(self%nparts_tot))
                if( allocated(outfile_body_local) )then
                    call job_descr%set('outfile', trim(trim(outfile_body_local)//int2str_pad(ipart,self%numlen)//'.txt'))
                endif
            endif
            if( part_params_present  )then
                do iadd=1,part_params(ipart)%size_of_chash()
                    key = part_params(ipart)%get_key(iadd)
                    val = part_params(ipart)%get(iadd)
                    call job_descr%set(key, val)
                end do
            endif
            call self%generate_script_1(job_descr, ipart, q_descr)
        end do
        if( .not. cchunkdistr )then
            call job_descr%delete('fromp')
            call job_descr%delete('top')
            call job_descr%delete('part')
            call job_descr%delete('nparts')
            if( allocated(outfile_body_local) )then
                call job_descr%delete('outfile')
                deallocate(outfile_body_local)
            endif
        endif
        if( part_params_present  )then
            do iadd=1,part_params(1)%size_of_chash()
                key = part_params(1)%get_key(iadd)
                call job_descr%delete(key)
            end do
        endif
        ! when we generate the scripts we also reset the number of available computing units
        if( .not. self%stream ) self%ncomputing_units_avail = self%ncomputing_units
    end subroutine generate_scripts

    !>  \brief  private part script generator
    subroutine generate_script_1( self, job_descr, ipart, q_descr )
        use simple_filehandling, only: get_fileunit, file_exists
        class(qsys_ctrl), intent(inout) :: self
        class(chash),     intent(in)    :: job_descr
        integer,          intent(in)    :: ipart
        class(chash),     intent(in)    :: q_descr
        character(len=512) :: io_msg
        integer :: ios, funit, val
        funit = get_fileunit()
        if( self%stream )then
            if( file_exists(self%script_names(ipart)) ) return
        endif
        open(unit=funit, file=self%script_names(ipart), iostat=ios, STATUS='REPLACE', action='WRITE', iomsg=io_msg)
        if( ios .ne. 0 )then
            close(funit)
            write(*,'(a)') 'simple_qsys_ctrl :: gen_qsys_script; Error when opening file for writing: '&
            //trim(self%script_names(ipart))//' ; '//trim(io_msg)
            stop
        endif
        ! need to specify shell
        write(funit,'(a)') '#!/bin/bash'
        ! write (run-time polymorphic) instructions to the qsys
        if( q_descr%get('qsys_name').ne.'local' )then
            call self%myqsys%write_instr(q_descr, fhandle=funit)
        else
            call self%myqsys%write_instr(job_descr, fhandle=funit)
        endif
        write(funit,'(a)',advance='yes') 'cd '//trim(self%pwd)
        write(funit,'(a)',advance='yes') ''
        ! compose the command line
        write(funit,'(a)',advance='no') trim(self%exec_binary)//' '//job_descr%chash2str()
        ! direct output
        write(funit,'(a)',advance='yes') ' > OUT'//int2str_pad(ipart,self%numlen)
        ! exit shell when done
        write(funit,'(a)',advance='yes') ''
        write(funit,'(a)',advance='yes') 'exit'
        flush(funit)
        close(funit)
        if( q_descr%get('qsys_name').eq.'local' )call chmod(trim(self%script_names(ipart)),'+x')
        if( ios .ne. 0 )then
            write(*,'(a)',advance='no') 'simple_qsys_scripts :: gen_qsys_script; Error'
            write(*,'(a)') 'chmoding submit script'//trim(self%script_names(ipart))
            stop
        endif
        ! when we generate the script we also unflag jobs_submitted and jobs_done
        self%jobs_done(ipart)      = .false.
        self%jobs_submitted(ipart) = .false.
    end subroutine generate_script_1

    !>  \brief  public script generator for single jobs
    subroutine generate_script_2( self, job_descr, q_descr, exec_bin, script_name, outfile )
        use simple_filehandling, only: get_fileunit, file_exists
        class(qsys_ctrl), intent(inout) :: self
        class(chash),     intent(in)    :: job_descr
        class(chash),     intent(in)    :: q_descr
        character(len=*), intent(in)    :: exec_bin, script_name, outfile
        character(len=512) :: io_msg
        integer :: ios, funit, val
        funit = get_fileunit()
        open(unit=funit, file=script_name, iostat=ios, STATUS='REPLACE', action='WRITE', iomsg=io_msg)
        if( ios .ne. 0 )then
            close(funit)
            write(*,'(a)') 'simple_qsys_ctrl :: gen_qsys_script; Error when opening file for writing: '&
            //trim(script_name)//' ; '//trim(io_msg)
            stop
        endif
        ! need to specify shell
        write(funit,'(a)') '#!/bin/bash'
        ! write (run-time polymorphic) instructions to the qsys
        if( q_descr%get('qsys_name').ne.'local' )then
            call self%myqsys%write_instr(q_descr, fhandle=funit)
        else
            call self%myqsys%write_instr(job_descr, fhandle=funit)
        endif
        write(funit,'(a)',advance='yes') 'cd '//trim(self%pwd)
        write(funit,'(a)',advance='yes') ''
        ! compose the command line
        write(funit,'(a)',advance='no') trim(exec_bin)//' '//job_descr%chash2str()
        ! direct output
        write(funit,'(a)',advance='yes') ' > '//outfile
        ! exit shell when done
        write(funit,'(a)',advance='yes') ''
        write(funit,'(a)',advance='yes') 'exit'
        flush(funit)
        close(funit)
        !flush(funit)
        if( q_descr%get('qsys_name').eq.'local' )call chmod(trim(script_name),'+x')
        if( ios .ne. 0 )then
            write(*,'(a)',advance='no') 'simple_qsys_ctrl :: generate_script_2; Error'
            write(*,'(a)') 'chmoding submit script'//trim(script_name)
            stop
        endif
    end subroutine generate_script_2

    ! SUBMISSION TO QSYS

    subroutine submit_scripts( self )
        use simple_filehandling, only: get_fileunit, fopen_err, del_file
        use simple_syscalls,     only: exec_cmdline
        use simple_qsys_local,   only: qsys_local
        class(qsys_ctrl),  intent(inout) :: self
        class(qsys_base),      pointer   :: pmyqsys
        character(len=STDLEN)            :: master_submit_script
        integer :: ipart, fnr, file_stat, chmod_stat, ios, val
        logical :: submit_or_not(self%fromto_part(1):self%fromto_part(2)), err
        ! master command line
        master_submit_script = trim(adjustl(self%pwd))//'/qsys_submit_jobs'
        call del_file(trim(master_submit_script))
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
        ! make the master submission script
        fnr = get_fileunit()
        open(unit=fnr, FILE=master_submit_script, STATUS='NEW', action='WRITE', iostat=file_stat)
        call fopen_err('simple_qsys_ctrl :: submit_scripts', file_stat )
        write(fnr,'(a)') '#!/bin/bash'
        do ipart=self%fromto_part(1),self%fromto_part(2)
            if( submit_or_not(ipart) )then
                select type( pmyqsys => self%myqsys )
                    class is(qsys_local)
                        write(fnr,'(a)') self%myqsys%submit_cmd()//' ./'//trim(adjustl(self%script_names(ipart)))//' &'
                    class DEFAULT
                        write(fnr,'(a)') self%myqsys%submit_cmd()//' ./'//trim(adjustl(self%script_names(ipart)))
                end select
            endif
        end do
        write(fnr,'(a)') 'exit'
        close( unit=fnr )
        flush(fnr)
        call chmod(master_submit_script,'+x')
        if( debug.or.global_debug )then
            call exec_cmdline('echo DISTRIBUTED MODE :: submitting scripts:')
            call exec_cmdline('ls -1 distr_simple_script_*')
        endif
        ! execute the master submission script
        call exec_cmdline(master_submit_script)
    end subroutine submit_scripts

    subroutine submit_script( self, script_name )
        use simple_syscalls,     only: exec_cmdline
        use simple_qsys_local,   only: qsys_local
        class(qsys_ctrl), intent(inout) :: self
        character(len=*), intent(in)    :: script_name
        class(qsys_base),      pointer  :: pmyqsys
        integer :: ios
        character(len=STDLEN) :: cmd
        select type( pmyqsys => self%myqsys )
            class is(qsys_local)
                cmd = self%myqsys%submit_cmd()//' '//trim(adjustl(self%pwd))//'/'//trim(adjustl(script_name))//' &'
            class DEFAULT
                cmd = self%myqsys%submit_cmd()//' '//trim(adjustl(self%pwd))//'/'//trim(adjustl(script_name))
        end select
        ! execute the command
        call exec_cmdline(cmd)
    end subroutine submit_script

    ! QUERIES

    subroutine update_queue( self )
        use simple_filehandling, only: file_exists
        class(qsys_ctrl),  intent(inout) :: self
        integer :: ipart, njobs_in_queue
        do ipart=self%fromto_part(1),self%fromto_part(2)
            if( .not. self%jobs_done(ipart) )then
                self%jobs_done(ipart) = file_exists(self%jobs_done_fnames(ipart))
            endif
            ! this one is for streaming
            if( self%jobs_done(ipart) ) self%jobs_submitted(ipart) = .true.
        end do
        njobs_in_queue = count(self%jobs_submitted .eqv. (.not. self%jobs_done))
        self%ncomputing_units_avail = self%ncomputing_units-njobs_in_queue
    end subroutine update_queue

    ! THE MASTER SCHEDULER

    subroutine schedule_jobs( self )
        use simple_syscalls, only: simple_sleep
        class(qsys_ctrl),  intent(inout) :: self
        do
            if( all(self%jobs_done) ) exit
            call self%update_queue
            call self%submit_scripts
            call simple_sleep(SHORTTIME)
        end do
    end subroutine schedule_jobs

    ! DESTRUCTOR

    !>  \brief  is a destructor
    subroutine kill( self )
        class(qsys_ctrl), intent(inout) :: self
        if( self%existence )then
            self%exec_binary            =  ''
            self%pwd                    =  ''
            self%myqsys                 => null()
            self%parts                  => null()
            self%fromto_part(2)         =  0
            self%ncomputing_units       =  0
            self%ncomputing_units_avail =  0
            self%numlen                 =  0
            deallocate(self%script_names, self%jobs_done, self%jobs_done_fnames, self%jobs_submitted)
            self%existence = .false.
        endif
    end subroutine kill

end module simple_qsys_ctrl
