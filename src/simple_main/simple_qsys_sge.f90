! batch-processing manager - SunGrid Engine
module simple_qsys_sge
use simple_qsys_base, only: qsys_base
use simple_chash,     only: chash
implicit none

public :: qsys_sge
private

integer, parameter :: MAXENVITEMS=100

type, extends(qsys_base) :: qsys_sge
    private
    type(chash) :: env !< defines the SGE environment
  contains
    procedure :: new         => new_sge_env
    procedure :: submit_cmd  => get_sge_submit_cmd
    procedure :: write_instr => write_sge_header
    procedure :: kill        => kill_sge_env
end type qsys_sge

contains
    
    !> \brief  is a constructor
    subroutine new_sge_env( self )
        class(qsys_sge), intent(inout) :: self
        ! make the container
        call self%env%new(MAXENVITEMS)
        ! define the environment:
        ! ### USER PARAMETERS
        call self%env%push('user_account',           '#$ -A')
        call self%env%push('user_email',             '#$ -M')
        call self%env%push('user_project',           '#$ -P')
        ! ### QSYS PARAMETERS
        call self%env%push('qsys_queue',             '#$ -q')
        call self%env%push('qsys_submit_cmd',         'qsub')
        ! ### JOB PARAMETERS
        call self%env%push('job_name',               '#$ -N')
        call self%env%push('job_addon_line',         '')
    end subroutine new_sge_env

    !> \brief  is a getter
    function get_sge_submit_cmd( self ) result( cmd )
        class(qsys_sge), intent(in) :: self
        character(len=:), allocatable :: cmd
        cmd = self%env%get('qsys_submit_cmd')
    end function get_sge_submit_cmd

    !> \brief  writes the header instructions
    subroutine write_sge_header( self, job_descr, fhandle )
        class(qsys_sge), intent(in)   :: self
        class(chash),      intent(in) :: job_descr
        integer, optional, intent(in) :: fhandle
        character(len=:), allocatable :: addon_line
        character(len=:), allocatable :: key, qsub_cmd, qsub_val, bind2socket
        integer :: i, which
        logical :: write2file
        write2file = .false.
        if( present(fhandle) ) write2file = .true.
        do i=1,job_descr%size_of_chash()
            key   = job_descr%get_key(i)

            print *, 'key from sge: ', key

            which = self%env%lookup(key)
            if( which > 0 )then
                qsub_cmd = self%env%get(which)
                qsub_val = job_descr%get(i)
                if( key .eq. 'job_addon_line' )then
                    if( allocated(addon_line) ) deallocate(addon_line)

                    print *, 'qsub val for key(job_addon_line): ', qsub_val
 
                    allocate( addon_line, source=qsub_val )
                    cycle
                else
                    if( write2file )then
                        write(fhandle,'(a)') qsub_cmd//' '//qsub_val
                    else
                        write(*,'(a)') qsub_cmd//' '//qsub_val
                    endif
                endif
                if(key .eq. 'job_cpus_per_task')then
                    if( allocated(bind2socket) ) deallocate(bind2socket)
                    allocate(bind2socket, source='#$ -binding linear:'//trim(qsub_val))
                endif
                deallocate(qsub_cmd,qsub_val)
            endif
            deallocate(key)
        end do
        ! write default instructions
        if( write2file )then
            write(fhandle,'(a)') '#$ -V'
            write(fhandle,'(a)') '#$ -S /bin/bash'
            write(fhandle,'(a)') '#$ -cwd'
            write(fhandle,'(a)') '#$ -o outfile -j y'
!            write(fhandle,'(a)') '#$ -pe mpi 1'
            if( allocated(bind2socket) )then
                 write(fhandle,'(a)') bind2socket
            endif
            if( allocated(addon_line) )then
                 write(fhandle, '(a)') addon_line
            endif
        else
            write(*,'(a)') '#$ -V'
            write(*,'(a)') '#$ -S /bin/bash'
            write(*,'(a)') '#$ -cwd'
            write(*,'(a)') '#$ -o outfile -j y'
 !           write(*,'(a)') '#$ -pe mpi 1'
            if( allocated(bind2socket) )then
                 write(*,'(a)') bind2socket
            endif
            if( allocated(addon_line) )then
                 write(*,'(a)') addon_line
            endif
        endif
    end subroutine write_sge_header
    
    !> \brief  is a destructor
    subroutine kill_sge_env( self )
        class(qsys_sge), intent(inout) :: self
        call self%env%kill
    end subroutine kill_sge_env

end module simple_qsys_sge
