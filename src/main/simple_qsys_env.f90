! batch-processing manager - environment module
module simple_qsys_env
include 'simple_lib.f08'
use simple_qsys_funs,    only: qsys_watcher,qsys_cleanup
use simple_qsys_factory, only: qsys_factory
use simple_qsys_base,    only: qsys_base
use simple_qsys_ctrl,    only: qsys_ctrl
use simple_parameters,   only: params_glob
implicit none

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
    procedure :: gen_scripts_and_schedule_jobs
    ! procedure :: gen_shm_scripts_and_schedule_jobs
    procedure :: exec_simple_prg_in_queue
    procedure :: get_qsys
    procedure :: kill
end type qsys_env
#include "simple_local_flags.inc"
contains

    subroutine new( self, nparts, stream, numlen, exec_bin )
        use simple_ori,        only: ori
        use simple_sp_project, only: sp_project
        class(qsys_env),             intent(inout) :: self
        integer,                     intent(in)    :: nparts
        logical,           optional, intent(in)    :: stream
        integer,           optional, intent(in)    :: numlen
        character(len=*),  optional, intent(in)    :: exec_bin
        type(ori)                     :: compenv_o
        type(sp_project)              :: spproj
        character(len=:), allocatable :: qsnam, tpi, hrs_str, mins_str, secs_str
        integer                       :: partsz, hrs, mins, secs
        real                          :: rtpi, tot_time_sec
        logical                       :: sstream
        integer, parameter            :: MAXENVKEYS = 30
        call self%kill
        sstream = .false.
        if( present(stream) ) sstream = stream
        self%nparts = nparts
        select case(trim(params_glob%split_mode))
            case('even')
                self%parts = split_nobjs_even(params_glob%nptcls, self%nparts)
                partsz     = self%parts(1,2) - self%parts(1,1) + 1
            case('singles')
                allocate(self%parts(params_glob%nptcls,2))
                self%parts(:,:) = 1
                partsz          = 1
            case('stream')
                self%nparts = params_glob%ncunits
                allocate(self%parts(params_glob%nptcls,2)) ! unused
                self%parts(:,:) = 1              ! unused
                partsz          = 1              ! unused
            case DEFAULT
                write(*,*) 'split_mode: ', trim(params_glob%split_mode)
                stop 'Unsupported split_mode'
        end select
        ! retrieve environment variables from file
        call self%qdescr%new(MAXENVKEYS)
        call spproj%read_segment('compenv', params_glob%projfile)
        compenv_o   = spproj%compenv%get_ori(1)
        self%qdescr = compenv_o%ori2chash()
        ! deal with time
        if( self%qdescr%isthere('time_per_image') )then
            tpi          = self%qdescr%get('time_per_image')
            rtpi         = str2real(tpi)
            tot_time_sec = rtpi*real(partsz)
            hrs          = int(tot_time_sec/3600.)
            hrs_str      = int2str(hrs)
            mins         = int((tot_time_sec - 3600.*real(hrs))/60.)
            mins_str     = int2str(mins)
            secs         = int(tot_time_sec - 3600.*real(hrs) - 60.*real(mins))
            secs_str     = int2str(secs)
            if( hrs > 23 )then
                call self%qdescr%set('job_time', '1-23:59:0')
            else
                call self%qdescr%set('job_time','0-'//hrs_str//':'//mins_str//':'//secs_str)
            endif
        endif
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
        call self%qdescr%set('job_cpus_per_task', int2str(params_glob%nthr))   ! overrides env file
        call self%qdescr%set('job_nparts',        int2str(params_glob%nparts)) ! overrides env file
        deallocate(qsnam)
        self%existence = .true.
    end subroutine new

    function exists( self ) result( is )
        class(qsys_env) :: self
        logical         :: is
        is = self%existence
    end function exists

    subroutine gen_scripts_and_schedule_jobs( self,  job_descr, part_params, algnfbody )
        class(qsys_env)            :: self
        class(chash)               :: job_descr
        class(chash),     optional :: part_params(self%nparts)
        character(len=*), optional :: algnfbody
        call qsys_cleanup
        call self%qscripts%generate_scripts(job_descr, trim(params_glob%ext), self%qdescr,&
        outfile_body=algnfbody, part_params=part_params)
        call self%qscripts%schedule_jobs
    end subroutine gen_scripts_and_schedule_jobs

    subroutine exec_simple_prg_in_queue( self, cline, outfile, finish_indicator, script_name )
        use simple_cmdline,   only: cmdline
        class(qsys_env)            , intent(inout)   :: self
        class(cmdline)             , intent(inout)   :: cline
        character(len=*)                             :: outfile
        character(len=*), optional , intent(in)   :: finish_indicator, script_name
        character(len=*), parameter   :: script_name_default = 'simple_script_single'
        type(chash)                   :: job_descr
        character(len=:), allocatable :: halt_ind, script_name_here
        if( present(finish_indicator) )then
            allocate(halt_ind, source=trim(finish_indicator))
        endif
        if( present(script_name) )then
            allocate(script_name_here, source=trim(script_name))
        else
            allocate(script_name_here, source=script_name_default)
        endif
        if( allocated(halt_ind) ) call del_file(halt_ind)
        call cline%gen_job_descr(job_descr)
        call self%qscripts%generate_script(job_descr, self%qdescr, self%simple_exec_bin, script_name_here, outfile)
        call wait_for_closure(script_name_here)
        call self%qscripts%submit_script(script_name_here)
        if( allocated(halt_ind) )then
            call qsys_watcher(halt_ind)
            call del_file(halt_ind)
        endif
    end subroutine exec_simple_prg_in_queue

    function get_qsys( self )result( qsys )
        use simple_qsys_local, only: qsys_local
        use simple_qsys_slurm, only: qsys_slurm
        use simple_qsys_pbs,   only: qsys_pbs
        use simple_qsys_sge,   only: qsys_sge
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
