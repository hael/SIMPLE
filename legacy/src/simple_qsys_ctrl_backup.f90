module simple_qsys_ctrl
use simple_jiffys     ! singleton
use simple_defs       ! singleton
use simple_qsys_base, only: qsys_base
use simple_chash,     only: chash
implicit none

public :: qsys_ctrl
private

type qsys_ctrl
    private
    character(len=STDLEN)          :: exec_binary = ''     !< binary to execute in parallel
                                                           !< for example trim(simplepath)//'/bin/simple_exec prg=prime2D'
    character(len=32), allocatable :: script_names(:)      !< file names of generated scripts
    character(len=STDLEN)          :: pwd = ''             !< working directory
    class(qsys_base), pointer      :: myqsys     => null() !< pointer to polymorphic qsys object
    integer, pointer               :: parts(:,:) => null() !< defines the fromp/top ranges for all partitions
    integer, allocatable           :: partszs(:)           !< size of partitions (number of particles in each partition)
    integer                        :: fromto_part(2) = 0   !< defines the range of partitions controlled by this instance
    integer                        :: first_ptcl     = 0   !< first particle index in first partition
    integer                        :: last_ptcl      = 0   !< last particle index in last partition
    integer                        :: npart          = 0   !< number of partitions controlled by this instance
    integer                        :: npart_tot      = 0   !< total number of partitions
    integer                        :: numlen         = 0   !< length of padded number string
    logical                        :: exists = .false.     !< indicates existence
  contains
    ! CONSTRUCTOR
    procedure          :: new
    ! SCRIPT GENERATORS
    procedure          :: generate_scripts
    procedure, private :: generate_script
    ! SUBMISSION TO QSYS
    procedure          :: submit_scripts
    ! QUERIES
    procedure          :: scripts_are_finished
    ! DESTRUCTOR
    procedure          :: kill
end type qsys_ctrl

interface qsys_ctrl
    module procedure constructor
end interface 

contains
    
    ! CONSTRUCTORS
    
    !>  \brief  is a constructor
    function constructor( exec_binary, qsys_obj, parts, fromto_part ) result( self )
        character(len=*), intent(in)         :: exec_binary
        class(qsys_base), intent(in), target :: qsys_obj
        integer,          intent(in), target :: parts(:,:)
        integer,          intent(in)         :: fromto_part(2)
        type(qsys_ctrl)     :: self
        call self%new(exec_binary, qsys_obj, parts, fromto_part)
    end function
    
    !>  \brief  is a constructor
    subroutine new( self, exec_binary, qsys_obj, parts, fromto_part )
        class(qsys_ctrl), intent(inout)      :: self
        character(len=*), intent(in)         :: exec_binary
        class(qsys_base), intent(in), target :: qsys_obj
        integer,          intent(in), target :: parts(:,:)
        integer,          intent(in)         :: fromto_part(2)
        integer :: alloc_stat, ipart
        call self%kill
        self%exec_binary =  exec_binary
        self%myqsys      => qsys_obj
        self%parts       => parts
        self%fromto_part =  fromto_part
        self%npart       =  fromto_part(2)-fromto_part(1)+1
        self%npart_tot   =  size(parts,1)
        self%numlen      =  len(int2str(self%npart_tot))
        self%first_ptcl  =  parts(fromto_part(1),1)
        self%last_ptcl   =  parts(fromto_part(2),2)
        ! allocate
        allocate( self%partszs(fromto_part(1):fromto_part(2)),&
                  self%script_names(fromto_part(1):fromto_part(2)), stat=alloc_stat)
        call alloc_err("In: simple_qsys_ctrl :: new", alloc_stat)
        ! get pwd
        call get_environment_variable('PWD', self%pwd)
        ! work out the sizes of the partitions
        do ipart=fromto_part(1),fromto_part(2)
            self%partszs = parts(ipart,2)-parts(ipart,1)+1
        end do
        self%exists = .true.
    end subroutine new
    
    ! SCRIPT GENERATION
    
    !>  \brief  public script generator
    subroutine generate_scripts( self, job_descr, outfile_body )
        class(qsys_ctrl),           intent(inout) :: self
        class(chash),               intent(inout) :: job_descr
        character(len=*), optional, intent(in)    :: outfile_body
        character(len=:), allocatable :: outfile_body_local
        integer :: ipart
        if( present(outfile_body) )then
            allocate(outfile_body_local, source=trim(outfile_body))
        else
            allocate(outfile_body_local, source='algndoc_')
        endif
        do ipart=self%fromto_part(1),self%fromto_part(2)
            call job_descr%set('fromp',   int2str(self%parts(ipart,1)))
            call job_descr%set('top',     int2str(self%parts(ipart,2)))
            call job_descr%set('part',    int2str(ipart))
            call job_descr%set('outfile', trim(outfile_body_local//int2str_pad(ipart,self%numlen)//'.txt'))
            call self%generate_script(job_descr, ipart)
        end do
        call job_descr%delete('fromp')
        call job_descr%delete('top')
        call job_descr%delete('part')
        call job_descr%delete('outfile')
        deallocate(outfile_body_local)
    end subroutine generate_scripts

    !>  \brief  private part script generator
    subroutine generate_script( self, job_descr, ipart )
        class(qsys_ctrl), intent(inout) :: self
        class(chash),     intent(in)    :: job_descr
        integer,          intent(in)    :: ipart
        character(len=32)  :: keys2print(4)
        character(len=512) :: io_msg
        integer :: ios, funit
        funit = get_fileunit()
        self%script_names(ipart) = 'distr_simple_script_'//int2str_pad(ipart,self%numlen)
        open(unit=funit, file=self%script_names(ipart), iostat=ios, status='replace', iomsg=io_msg)
        if( ios .ne. 0 )then
            close(funit)
            write(*,'(a)') 'simple_qsys_ctrl :: gen_qsys_script; Error when opening file for writing: '&
            //trim(self%script_names(ipart))//' ; '//trim(io_msg)
            stop
        endif
        ! need to specify shell
        write(funit,'(a)') '#!/bin/bash'
        ! write (run-time polymorphic) instructions to the qsys
        call self%myqsys%write_instr(job_descr, fhandle=funit)
        write(funit,'(a)') 'cd '//trim(self%pwd)
        keys2print(1) = 'fromp'
        keys2print(2) = 'top'
        keys2print(3) = 'part'
        keys2print(4) = 'outfile'
        ! compose the command line
        write(funit,'(a)',advance='no') trim(self%exec_binary)//' '//job_descr%chash2str() 
        ! call job_descr%print_key_val_pairs(keys2print, fhandle=funit, oneline=.true.)
        ! direct output
        write(funit,'(a)') '> OUT'//int2str_pad(ipart,self%numlen)
        ! exit shell when done
        write(funit,'(a)') 'exit'
        close(funit)
        call chmod(self%script_names(ipart),'+x',ios)
        if( ios .ne. 0 )then
            write(*,'(a)',advance='no') 'simple_qsys_scripts :: gen_qsys_script; Error'
            write(*,'(a)') 'chmoding submit script'//trim(self%script_names(ipart))
            stop
        endif
    end subroutine generate_script
    
    ! SUBMISSION TO QSYS
    
    subroutine submit_scripts( self )
        use simple_qsys_local, only: qsys_local
        class(qsys_ctrl), intent(in)  :: self
        integer                       :: cstat, estat, ipart
        character(len=100)            :: cmsg
        character(len=:), allocatable :: script_exec_binary
        class(qsys_base), pointer     :: pmyqsys
        do ipart=self%fromto_part(1),self%fromto_part(2)
            select type(pmyqsys => self%myqsys)
                class is(qsys_local)
                    allocate(script_exec_binary, source=self%myqsys%submit_cmd()//' ./'//trim(self%script_names(ipart))//' &')
                class DEFAULT
                    allocate(script_exec_binary, source=self%myqsys%submit_cmd()//' ./'//trim(self%script_names(ipart)))
            end select
            call execute_command_line(script_exec_binary, exitstat=estat, cmdstat=cstat, cmdmsg=cmsg)
            if( cstat > 0 )then
                write(*,'(a,a)') 'simple_qsys_ctrl :: submit_scripts; command execution failed with error ', trim(cmsg)
            elseif( cstat < 0 )then
                write(*,'(a)') 'simple_qsys_ctrl :: submit_scripts; command execution not supported'
            endif
            deallocate(script_exec_binary)
        end do
    end subroutine submit_scripts
    
    ! QUERIES
    
    function scripts_are_finished( self ) result( are_finished )
        class(qsys_ctrl), intent(in)  :: self
        character(len=:), allocatable :: job_done_fname
        integer :: ipart
        logical :: jobs_done(self%fromto_part(1):self%fromto_part(2)), are_finished
        do ipart=self%fromto_part(1),self%fromto_part(2)
            allocate(job_done_fname, source='JOB_FINISHED_'//int2str_pad(ipart,self%numlen))
            jobs_done(ipart) = file_exists(job_done_fname)
            deallocate(job_done_fname)
        end do
        are_finished = all(jobs_done)
    end function scripts_are_finished

    ! DESTRUCTOR
    
    !>  \brief  is a destructor
    subroutine kill( self )
        class(qsys_ctrl), intent(inout) :: self
        if( self%exists )then
            self%exec_binary       =  ''
            self%pwd            =  ''
            self%fromto_part(2) =  0
            self%first_ptcl     =  0
            self%last_ptcl      =  0 
            self%npart          =  0
            self%numlen         =  0
            self%myqsys         => null()
            self%parts          => null()
            deallocate(self%script_names, self%partszs)
            self%exists = .false.
        endif
    end subroutine kill

end module simple_qsys_ctrl
