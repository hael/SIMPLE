module simple_qsys_ctrl
use simple_defs
use simple_qsys_base, only: qsys_base
use simple_chash,     only: chash
use simple_strings,   only: int2str, int2str_pad
implicit none

public :: qsys_ctrl
private

integer, parameter :: SHORTTIME = 5

type qsys_ctrl
    private
    character(len=STDLEN)          :: exec_binary      = ''      !< binary to execute in parallel
                                                                 !< trim(simplepath)//'/bin/simple_exec'
    character(len=32), allocatable :: script_names(:)            !< file names of generated scripts
    character(len=STDLEN)          :: pwd              = ''      !< working directory
    class(qsys_base), pointer      :: myqsys           => null() !< pointer to polymorphic qsys object
    integer, pointer               :: parts(:,:)       => null() !< defines the fromp/top ranges for all partitions
    logical, pointer               :: lmask_stream(:)  => null() !< pointer to the logical stream mask
    logical, allocatable           :: jobs_done(:)               !< to indicate completion of distributed scripts
    logical, allocatable           :: jobs_submitted(:)          !< to indicate which jobs have been submitted
    integer                        :: fromto_part(2)         = 0 !< defines the range of partitions controlled by this instance
    integer                        :: nparts_tot             = 0 !< total number of partitions
    integer                        :: ncomputing_units       = 0 !< number of computing units
    integer                        :: ncomputing_units_avail = 0 !< number of available units
    integer                        :: numlen = 0                 !< length of padded number string
    logical                        :: exists = .false.           !< indicates existence
  contains
    ! CONSTRUCTOR
    procedure          :: new
    ! GETTER
    procedure          :: get_exec_bin
    ! SCRIPT GENERATORS
    procedure          :: generate_scripts
    procedure, private :: generate_script
    ! SUBMISSION TO QSYS
    procedure, private :: submit_scripts
    ! QUERIES
    procedure, private :: update_queue
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
    function constructor( exec_binary, qsys_obj, parts, fromto_part, ncomputing_units, lmask_stream ) result( self )
        character(len=*),          intent(in)    :: exec_binary      !< the binary that we want to execute in parallel
        class(qsys_base),  target, intent(in)    :: qsys_obj         !< the object that defines the qeueuing system
        integer,           target, intent(in)    :: parts(:,:)       !< defines the start_ptcl/stop_ptcl ranges  
        integer,                   intent(in)    :: fromto_part(2)   !< defines the range of partitions controlled by this object
        integer,                   intent(in)    :: ncomputing_units !< number of computing units (<= the number of parts controlled)
        logical, optional, target, intent(inout) :: lmask_stream(:)  !< logical mask for stream processing
        type(qsys_ctrl) :: self
        call self%new(exec_binary, qsys_obj, parts, fromto_part, ncomputing_units)
    end function constructor
    
    !>  \brief  is a constructor
    subroutine new( self, exec_binary, qsys_obj, parts, fromto_part, ncomputing_units, lmask_stream )
        use simple_jiffys, only: alloc_err
        class(qsys_ctrl),          intent(inout) :: self             !< the instance
        character(len=*),          intent(in)    :: exec_binary      !< the binary that we want to execute in parallel
        class(qsys_base),  target, intent(in)    :: qsys_obj         !< the object that defines the qeueuing system
        integer,           target, intent(in)    :: parts(:,:)       !< defines the start_ptcl/stop_ptcl ranges  
        integer,                   intent(in)    :: fromto_part(2)   !< defines the range of partitions controlled by this object
        integer,                   intent(in)    :: ncomputing_units !< number of computing units (<= the number of parts controlled)
        logical, optional, target, intent(inout) :: lmask_stream(:)  !< logical mask for stream processing
        integer :: alloc_stat, ipart
        call self%kill
        self%exec_binary            =  exec_binary
        self%myqsys                 => qsys_obj
        self%parts                  => parts
        self%fromto_part            =  fromto_part
        self%nparts_tot             =  size(parts,1)
        self%numlen                 =  len(int2str(self%nparts_tot))
        self%ncomputing_units       =  ncomputing_units
        if( present(lmask_stream) )then
            self%ncomputing_units_avail =  ncomputing_units - count(lmask_stream)
            self%lmask_stream           => lmask_stream
        else
            self%ncomputing_units_avail =  ncomputing_units
            self%lmask_stream           => null()
        endif        
        ! allocate
        allocate(   self%jobs_done(fromto_part(1):fromto_part(2)),&
                    self%jobs_submitted(fromto_part(1):fromto_part(2)),&
                    self%script_names(fromto_part(1):fromto_part(2)), stat=alloc_stat)
        call alloc_err("In: simple_qsys_ctrl :: new", alloc_stat)
        self%jobs_done      = .false.
        self%jobs_submitted = .false.
        self%script_names   = ''
        ! get pwd
        call get_environment_variable('PWD', self%pwd)
        self%exists = .true.
    end subroutine new

    ! GETTER

    !>  \brief  for getting the execute simple binary
    function get_exec_bin( self ) result( exec_bin )
        class(qsys_ctrl), intent(in) :: self
        character(len=STDLEN)        :: exec_bin
        exec_bin = self%exec_binary
    end function get_exec_bin
    
    ! SCRIPT GENERATORS
    
    !>  \brief  public script generator
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
        if( present(outfile_body) )then
            allocate(outfile_body_local, source=trim(outfile_body))
        endif
        part_params_present  = present(part_params)
        do ipart=self%fromto_part(1),self%fromto_part(2)
            if( associated(self%lmask_stream) )then
                if( .not. self%lmask_stream(ipart) ) cycle
            endif
            call job_descr%set('fromp',   int2str(self%parts(ipart,1)))
            call job_descr%set('top',     int2str(self%parts(ipart,2)))
            call job_descr%set('part',    int2str(ipart))
            call job_descr%set('nparts',  int2str(self%nparts_tot))
            if( allocated(outfile_body_local) )then
                call job_descr%set('outfile', trim(trim(outfile_body_local)//int2str_pad(ipart,self%numlen)//'.txt'))
            endif
            if( part_params_present  )then
                do iadd=1,part_params(ipart)%size_of_chash()
                    key = part_params(ipart)%get_key(iadd)
                    val = part_params(ipart)%get(iadd)
                    call job_descr%set(key, val)
                end do
            endif
            call self%generate_script(job_descr, ipart, q_descr)
        end do
        call job_descr%delete('fromp')
        call job_descr%delete('top')
        call job_descr%delete('part')
        call job_descr%delete('nparts')
        if( allocated(outfile_body_local) )then
            call job_descr%delete('outfile')
            deallocate(outfile_body_local)
        endif
        if( part_params_present  )then
            do iadd=1,part_params(1)%size_of_chash()
                key = part_params(1)%get_key(iadd)
                call job_descr%delete(key)
            end do
        endif
        if( .not. associated(self%lmask_stream) )then
            ! when we generate the scripts we also reset the number of available computing units
            self%ncomputing_units_avail = self%ncomputing_units
        endif
    end subroutine generate_scripts

    !>  \brief  private part script generator
    subroutine generate_script( self, job_descr, ipart, q_descr )
        use simple_filehandling, only: get_fileunit
        class(qsys_ctrl), intent(inout) :: self
        class(chash),     intent(in)    :: job_descr
        integer,          intent(in)    :: ipart
        class(chash),     intent(in)    :: q_descr
        character(len=512) :: io_msg
        integer :: ios, funit
        funit = get_fileunit()
        self%script_names(ipart) = 'distr_simple_script_'//int2str_pad(ipart,self%numlen)
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
        close(funit)
        call flush(funit)
        call chmod(self%script_names(ipart),'+x',ios)
        if( ios .ne. 0 )then
            write(*,'(a)',advance='no') 'simple_qsys_scripts :: gen_qsys_script; Error'
            write(*,'(a)') 'chmoding submit script'//trim(self%script_names(ipart))
            stop
        endif
        ! when we generate the script we also unflag jobs_submitted and jobs_done
        self%jobs_done(ipart)      = .false.
        self%jobs_submitted(ipart) = .false.
    end subroutine generate_script

    ! SUBMISSION TO QSYS

    subroutine submit_scripts( self )
        use simple_filehandling, only: get_fileunit, fopen_err
        use simple_syscalls,     only: exec_cmdline
        use simple_qsys_local,   only: qsys_local
        class(qsys_ctrl),  intent(inout) :: self
        class(qsys_base),      pointer   :: pmyqsys
        character(len=STDLEN), parameter :: master_submit_script = './qsys_submit_jobs'
        integer :: ipart, fnr, file_stat, chmod_stat, ios
        logical :: submit_or_not(self%fromto_part(1):self%fromto_part(2)), err
        ! make a submission mask
        submit_or_not = .false.
        do ipart=self%fromto_part(1),self%fromto_part(2)
            if( associated(self%lmask_stream) )then
                if( .not. self%lmask_stream(ipart) ) cycle
            endif
            if( .not. self%jobs_submitted(ipart) .and. self%ncomputing_units_avail > 0 )then
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
        open(unit=fnr, FILE=master_submit_script, STATUS='REPLACE', action='WRITE', iostat=file_stat)
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
        call chmod(master_submit_script,'+x',ios)
        if( ios .ne. 0 )then
            write(*,'(a)',advance='no') 'ERROR, simple_qsys_ctrl :: submit_scripts :: make_master_submit_script, '
            write(*,'(a)') 'chmoding master submit script '//trim(master_submit_script)
            stop
        endif
        call exec_cmdline('echo DISTRIBUTED MODE :: submitting scripts:')
        call exec_cmdline('ls -1 distr_simple_script_*')
        ! execute the master submission script
        call exec_cmdline(master_submit_script)
    end subroutine submit_scripts
    
    ! QUERIES

    subroutine update_queue( self )
        use simple_filehandling, only: file_exists
        class(qsys_ctrl),  intent(inout) :: self
        character(len=:), allocatable    :: job_done_fname
        logical :: ltmp
        integer :: ipart, njobs_in_queue
        do ipart=self%fromto_part(1),self%fromto_part(2)
            allocate(job_done_fname, source='JOB_FINISHED_'//int2str_pad(ipart,self%numlen))
            self%jobs_done(ipart) = file_exists(job_done_fname)
            ltmp = self%jobs_done(ipart)
            if( associated(self%lmask_stream) )then
                ! any .false. element in the streaming mask means that the job is done 
                ! even though the file does not exist
                if( .not. self%lmask_stream(ipart) ) self%jobs_done(ipart) = .true.
                ! if the file REALLY exists, the job is done and the streaming mask
                ! needs to be updated
                if( ltmp ) self%lmask_stream(ipart) = .false.
            endif
            deallocate(job_done_fname)
        end do
        njobs_in_queue = count(self%jobs_submitted .eqv. (.not. self%jobs_done))
        self%ncomputing_units_avail = self%ncomputing_units-njobs_in_queue
    end subroutine update_queue

    ! THE MASTER SCHEDULER

    subroutine schedule_jobs( self )
        class(qsys_ctrl),  intent(inout) :: self
        if( associated(self%lmask_stream) )then
            call self%update_queue
            call self%submit_scripts
        else
            do
                if( all(self%jobs_done) ) exit
                call self%update_queue
                call self%submit_scripts
                call sleep(SHORTTIME)
            end do
        endif
    end subroutine schedule_jobs
    
    ! DESTRUCTOR
    
    !>  \brief  is a destructor
    subroutine kill( self )
        class(qsys_ctrl), intent(inout) :: self
        if( self%exists )then
            self%exec_binary            =  ''
            self%pwd                    =  ''
            self%myqsys                 => null()
            self%parts                  => null()
            self%lmask_stream           => null()
            self%fromto_part(2)         =  0
            self%ncomputing_units       =  0
            self%ncomputing_units_avail =  0
            self%numlen                 =  0
            deallocate(self%script_names, self%jobs_done, self%jobs_submitted)
            self%exists = .false.
        endif
    end subroutine kill

end module simple_qsys_ctrl
