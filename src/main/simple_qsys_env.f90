! batch-processing manager - environment module

module simple_qsys_env
#include "simple_lib.f08"
use simple_qsys_funs,    only: qsys_watcher,qsys_cleanup,parse_env_file
use simple_map_reduce   ! use all in there
use simple_qsys_factory, only: qsys_factory
use simple_qsys_base,    only: qsys_base
use simple_qsys_ctrl,    only: qsys_ctrl
use simple_chash,        only: chash
use simple_params,       only: params

implicit none
#include "simple_local_flags.inc"

type :: qsys_env
    integer, allocatable,      public  :: parts(:,:)
    type(qsys_ctrl),           public  :: qscripts
    type(chash),               public  :: qdescr
    character(len=STDLEN),     private :: simple_exec_bin
    type(qsys_factory),        private :: qsys_fac
    class(qsys_base), pointer, private :: myqsys=>null()
    logical,                   private :: existence = .false.
  contains
    procedure :: new
    procedure :: exists
    procedure :: gen_scripts_and_schedule_jobs
    procedure :: exec_simple_prg_in_queue
    procedure :: kill
end type qsys_env

contains

    subroutine new( self, p_master, stream, numlen )
        class(qsys_env)               :: self
        class(params)                 :: p_master
        logical, optional             :: stream
        integer, optional             :: numlen
        character(len=:), allocatable :: qsnam, tpi, hrs_str, mins_str, secs_str
        integer                       :: partsz, hrs, mins, secs, nparts
        real                          :: rtpi, tot_time_sec
        logical                       :: sstream
        integer, parameter            :: MAXENVKEYS = 30
        call self%kill
        sstream = .false.
        if( present(stream) ) sstream = stream
        nparts = p_master%nparts
        select case(p_master%split_mode)
            case('even')
                self%parts = split_nobjs_even(p_master%nptcls, nparts)
                partsz     = self%parts(1,2) - self%parts(1,1) + 1
            case('chunk')
                self%parts = split_nobjs_in_chunks(p_master%nptcls, p_master%chunksz)
                partsz     = p_master%chunksz
            case('singles')
                allocate(self%parts(p_master%nptcls,2))
                self%parts(:,:) = 1
                partsz          = 1
            case('stream')
                nparts = p_master%ncunits
                allocate(self%parts(p_master%nptcls,2)) ! unused
                self%parts(:,:) = 1                     ! unused
                partsz          = 1                     ! unused
            case DEFAULT
                write(*,*) 'split_mode: ', trim(p_master%split_mode)
                stop 'Unsupported split_mode'
        end select
        ! retrieve environment variables from file
        call self%qdescr%new(MAXENVKEYS)
        call parse_env_file(self%qdescr) ! parse .env file
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
        self%simple_exec_bin = trim(self%qdescr%get('simple_path'))//'/bin/simple_exec'
        if( present(numlen) )then
            call self%qscripts%new(self%simple_exec_bin, self%myqsys, self%parts,&
            &[1, nparts], p_master%ncunits, sstream, numlen)
        else
            call self%qscripts%new(self%simple_exec_bin, self%myqsys, self%parts,&
            &[1, nparts], p_master%ncunits, sstream)
        endif
        call self%qdescr%set('job_cpus_per_task', int2str(p_master%nthr))   ! overrides env file
        call self%qdescr%set('job_nparts',        int2str(p_master%nparts)) ! overrides env file
        deallocate(qsnam)
        self%existence = .true.
    end subroutine new

    function exists( self ) result( is )
        class(qsys_env) :: self
        logical         :: is
        is = self%existence
    end function exists

    subroutine gen_scripts_and_schedule_jobs( self, p_master, job_descr, part_params, algnfbody, ext_meta, chunkdistr )
        class(qsys_env)            :: self
        class(params)              :: p_master
        class(chash)               :: job_descr
        class(chash),     optional :: part_params(p_master%nparts)
        character(len=*), optional :: algnfbody
        character(len=4), optional :: ext_meta
        logical,          optional :: chunkdistr
        call qsys_cleanup(p_master)
        call self%qscripts%generate_scripts(job_descr, p_master%ext, self%qdescr,&
        outfile_body=algnfbody, outfile_ext=ext_meta, part_params=part_params, chunkdistr=chunkdistr)
        call self%qscripts%schedule_jobs
    end subroutine gen_scripts_and_schedule_jobs

    subroutine exec_simple_prg_in_queue( self, cline, outfile, finish_indicator, script_name )
        use simple_cmdline,   only: cmdline
        class(qsys_env)            :: self
        class(cmdline)             :: cline
        character(len=*)           :: outfile
        character(len=*), optional :: finish_indicator, script_name
        character(len=STDLEN), parameter   :: script_name_default = 'simple_script_single'
        type(chash)                        :: job_descr
        character(len=:),      allocatable :: halt_ind, script_name_here
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
        call wait_for_closure(script_name_here) !!!!!
        call self%qscripts%submit_script(script_name_here)
        if( allocated(halt_ind) )then
            call qsys_watcher(halt_ind)
            call del_file(halt_ind)
        endif
    end subroutine exec_simple_prg_in_queue

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
