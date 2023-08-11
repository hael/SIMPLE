! batch-processing manager - environment module
module simple_qsys_env
include 'simple_lib.f08'
use simple_qsys_funs,    only: qsys_watcher,qsys_cleanup
use simple_qsys_factory, only: qsys_factory
use simple_qsys_base,    only: qsys_base
use simple_qsys_local,   only: qsys_local
use simple_qsys_slurm,   only: qsys_slurm
use simple_qsys_pbs,     only: qsys_pbs
use simple_qsys_sge,     only: qsys_sge
use simple_qsys_ctrl,    only: qsys_ctrl
use simple_parameters,   only: parameters, params_glob
implicit none

public :: qsys_env
private
#include "simple_local_flags.inc"

type :: qsys_env
    integer, allocatable,      public  :: parts(:,:)
    type(qsys_ctrl),           public  :: qscripts
    type(chash),               public  :: qdescr
    character(len=STDLEN),     private :: simple_exec_bin
    type(qsys_factory),        private :: qsys_fac
    class(qsys_base), pointer, private :: myqsys=>null()
    integer,                   private :: nparts
    logical,                   private :: existence = .false.
    contains
    procedure :: new
    procedure :: exists
    procedure :: gen_script
    procedure :: gen_scripts_and_schedule_jobs
    procedure :: exec_simple_prg_in_queue
    procedure :: exec_simple_prg_in_queue_async
    procedure :: exec_simple_prgs_in_queue_async
    procedure :: get_qsys
    procedure :: get_navail_computing_units
    procedure :: kill
end type qsys_env

contains

    subroutine new( self, nparts, stream, numlen, nptcls, exec_bin, qsys_name, qsys_nthr )
        use simple_sp_project, only: sp_project
        class(qsys_env),             intent(inout) :: self
        integer,                     intent(in)    :: nparts
        logical,           optional, intent(in)    :: stream
        integer,           optional, intent(in)    :: numlen, nptcls, qsys_nthr
        character(len=*),  optional, intent(in)    :: exec_bin
        character(len=*),  optional, intent(in)    :: qsys_name ! to override the qsys read in
        type(ori)                     :: compenv_o
        type(sp_project)              :: spproj
        character(len=:), allocatable :: qsnam, tpi, hrs_str, mins_str, secs_str
        integer                       :: partsz, hrs, mins, secs, nptcls_here
        real                          :: rtpi, tot_time_sec
        logical                       :: sstream
        integer, parameter            :: MAXENVKEYS = 30
        call self%kill
        sstream = .false.
        if( present(stream) ) sstream = stream
        self%nparts = nparts
        nptcls_here = params_glob%nptcls
        if( present(nptcls) ) nptcls_here = nptcls
        select case(trim(params_glob%split_mode))
            case('even')
                self%parts = split_nobjs_even(nptcls_here, self%nparts)
                partsz     = self%parts(1,2) - self%parts(1,1) + 1
            case('singles')
                allocate(self%parts(nptcls_here,2))
                self%parts(:,:) = 1
                partsz          = 1
            case('stream')
                self%nparts = params_glob%ncunits
                allocate(self%parts(nptcls_here,2)) ! unused
                self%parts(:,:) = 1                 ! unused
                partsz          = 1                 ! unused
            case DEFAULT
                THROW_HARD('split_mode: '//trim(params_glob%split_mode)//' is unsupported; new')
        end select
        ! retrieve environment variables from file
        call self%qdescr%new(MAXENVKEYS)
        call spproj%read_segment('compenv', params_glob%projfile)
        call spproj%compenv%get_ori(1, compenv_o)
        self%qdescr = compenv_o%ori2chash()
        ! deal with time
        if( self%qdescr%isthere('time_per_image') )then
            tpi          = self%qdescr%get('time_per_image')
            rtpi         = str2real(tpi)
            tot_time_sec = rtpi*real(partsz)
            if( self%qdescr%isthere('walltime') )then
                tot_time_sec = min(tot_time_sec, str2real(self%qdescr%get('walltime')))
            else
                tot_time_sec = min(tot_time_sec, real(WALLTIME_DEFAULT))
            endif
            tot_time_sec = min(tot_time_sec, real(params_glob%walltime)) ! command line override
            hrs          = int(tot_time_sec/3600.)
            hrs_str      = int2str(hrs)
            mins         = int((tot_time_sec - 3600.*real(hrs))/60.)
            mins_str     = int2str(mins)
            secs         = int(tot_time_sec - 3600.*real(hrs) - 60.*real(mins))
            secs_str     = int2str(secs)
            call self%qdescr%set('job_time','0-'//hrs_str//':'//mins_str//':'//secs_str)
        endif
        if( present(qsys_name) ) call self%qdescr%set('qsys_name', qsys_name)
        qsnam = self%qdescr%get('qsys_name')
        call self%qsys_fac%new(qsnam, self%myqsys)
        ! create the user specific qsys and qsys controller (script generator)
        if(present(exec_bin))then
                self%simple_exec_bin = filepath(trim(self%qdescr%get('simple_path')),'bin',trim(exec_bin), nonalloc=.true.)
        else
                self%simple_exec_bin = filepath(trim(self%qdescr%get('simple_path')),'bin','simple_private_exec', nonalloc=.true.)
        endif
        if( present(numlen) )then
            call self%qscripts%new(self%simple_exec_bin, self%myqsys, self%parts,&
            &[1, self%nparts], params_glob%ncunits, sstream, numlen)
        else
            call self%qscripts%new(self%simple_exec_bin, self%myqsys, self%parts,&
            &[1, self%nparts], params_glob%ncunits, sstream)
        endif
        if(present(qsys_nthr)) then
            call self%qdescr%set('job_cpus_per_task', int2str(qsys_nthr))         ! overrides env file and params_glob
        else
            call self%qdescr%set('job_cpus_per_task', int2str(params_glob%nthr))  ! overrides env file
        endif
        call self%qdescr%set('job_nparts', int2str(params_glob%nparts)) ! overrides env file
        deallocate(qsnam)
        call compenv_o%kill
        call spproj%kill
        self%existence = .true.
    end subroutine new

    function exists( self ) result( is )
        class(qsys_env) :: self
        logical         :: is
        is = self%existence
    end function exists

    subroutine gen_script( self, cline, script_name, prg_output )
        use simple_cmdline, only: cmdline
        class(qsys_env)              :: self
        class(cmdline)               :: cline
        character(len=*), intent(in) :: script_name, prg_output
        call self%qscripts%generate_script(cline, self%qdescr, script_name, prg_output)
    end subroutine gen_script

    subroutine gen_scripts_and_schedule_jobs( self,  job_descr, part_params, algnfbody, array, extra_params)
        class(qsys_env)                        :: self
        class(chash)                           :: job_descr
        class(chash),     optional             :: part_params(self%nparts)
        type(parameters), optional, intent(in) :: extra_params
        character(len=*), optional             :: algnfbody
        logical,          optional             :: array
        logical :: aarray
        aarray = .false.
        if( present(array) ) aarray = array
        ! we only support array execution by SLURM
        select type(pmyqsys => self%myqsys)
            class is(qsys_local)
                aarray = .false.
            class is(qsys_slurm)
                ! keep aarray value
            class is(qsys_sge)
                aarray = .false.
            class is(qsys_pbs)
                aarray = .false.
        end select
        call qsys_cleanup
        if( aarray )then
            call self%qscripts%generate_array_script(job_descr, trim(params_glob%ext), self%qdescr,&
            &outfile_body=algnfbody, part_params=part_params)
            call self%qscripts%schedule_array_jobs
        else
            call self%qscripts%generate_scripts(job_descr, trim(params_glob%ext), self%qdescr,&
            &outfile_body=algnfbody, part_params=part_params, extra_params=extra_params)
            call self%qscripts%schedule_jobs
        endif
    end subroutine gen_scripts_and_schedule_jobs

    subroutine exec_simple_prg_in_queue( self, cline, finish_indicator )
        use simple_cmdline, only: cmdline
        class(qsys_env),  intent(inout) :: self
        class(cmdline),   intent(inout) :: cline
        character(len=*), intent(in)    :: finish_indicator
        character(len=*), parameter     :: SCRIPT_NAME = 'simple_script_single'
        type(chash)                     :: job_descr
        call del_file(finish_indicator)
        call cline%gen_job_descr(job_descr)
        call self%qscripts%generate_script(job_descr, self%qdescr, self%simple_exec_bin, SCRIPT_NAME)
        call wait_for_closure(SCRIPT_NAME)
        call self%qscripts%submit_script(SCRIPT_NAME)
        call qsys_watcher(finish_indicator)
        call del_file(finish_indicator)
        call job_descr%kill
    end subroutine exec_simple_prg_in_queue

    subroutine exec_simple_prg_in_queue_async( self, cline, script_name, outfile )
        use simple_cmdline, only: cmdline
        class(qsys_env),            intent(inout) :: self
        class(cmdline),             intent(in)    :: cline
        character(len=*),           intent(in)    :: script_name, outfile
        type(chash) :: job_descr
        call cline%gen_job_descr(job_descr)
        call self%qscripts%generate_script(job_descr, self%qdescr, self%simple_exec_bin, script_name,&
            & outfile=outfile)
        call wait_for_closure(script_name)
        call self%qscripts%submit_script(script_name)
        call job_descr%kill
    end subroutine exec_simple_prg_in_queue_async

    !>  To submit a list of jobs asynchronously
    subroutine exec_simple_prgs_in_queue_async( self, clines, script_name, outfile )
        use simple_cmdline, only: cmdline
        class(qsys_env),            intent(inout) :: self
        type(cmdline), allocatable, intent(in)    :: clines(:)
        character(len=*),           intent(in)    :: script_name, outfile
        type(chash), allocatable :: jobs_descr(:)
        integer :: i, njobs
        njobs = size(clines)
        allocate(jobs_descr(njobs))
        do i = 1,njobs
            call clines(i)%gen_job_descr(jobs_descr(i))
        enddo
        call self%qscripts%generate_script(jobs_descr, self%qdescr, self%simple_exec_bin, script_name, outfile)
        call wait_for_closure(script_name)
        call self%qscripts%submit_script(script_name)
        do i = 1,njobs
            call jobs_descr(i)%kill
        enddo
        deallocate(jobs_descr)
    end subroutine exec_simple_prgs_in_queue_async

    function get_qsys( self )result( qsys )
        class(qsys_env), intent(in)   :: self
        character(len=:), allocatable :: qsys
        select type(pmyqsys => self%myqsys)
            class is(qsys_local)
                qsys = 'local'
            class is(qsys_slurm)
                qsys = 'slurm'
            class is(qsys_sge)
                qsys = 'sge'
            class is(qsys_pbs)
                qsys = 'pbs'
        end select
    end function get_qsys

    integer function get_navail_computing_units( self )
        class(qsys_env), intent(in) :: self
        get_navail_computing_units = self%qscripts%get_ncomputing_units_avail()
    end function get_navail_computing_units

    subroutine kill( self )
        class(qsys_env) :: self
        if( self%existence )then
            deallocate(self%parts)
            call self%qscripts%kill
            call self%qdescr%kill
            call self%qsys_fac%kill
            self%myqsys => null()
            self%existence = .false.
        endif
    end subroutine kill

end module simple_qsys_env
