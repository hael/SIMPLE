! batch-processing manager - functions
module simple_qsys_funs
include 'simple_lib.f08'
implicit none

interface qsys_watcher
    module procedure qsys_watcher_1
    module procedure qsys_watcher_2
end interface

integer, parameter :: SHORTTIME = 3

contains

    subroutine qsys_cleanup( flag  )
        use simple_parameters, only: params_glob
        logical, optional, intent(in) :: flag
        integer, parameter :: NUMLEN_STATE = 2, NUMLEN_ITER = 3
        logical :: nokeep
        nokeep=.true.
        if(present(flag))nokeep=flag
        ! single files
        call del_file('FOO')
        call del_file('fort.0')
        call del_file('ftab_from_sys_find_last_fname.txt')
        call del_file('VOLASSEMBLE')
        call del_file('CAVGASSEMBLE')
        call del_file('SYMSRCH')
        call del_file(SIMPLE_SUBPROC_OUT)
        ! state numbered files
        call del_files('VOLASSEMBLE_FINISHED_STATE',     params_glob%nstates, numlen=NUMLEN_STATE)
        if(nokeep)call del_files('simple_script_state',  params_glob%nstates, numlen=NUMLEN_STATE)
        ! part numbered files
        call del_files('OUT',                            params_glob%nparts)
        call del_files('algndoc_',                       params_glob%nparts, ext=trim(METADATA_EXT))
        call del_files('JOB_FINISHED_',                  params_glob%nparts)
        if(nokeep)call del_files('distr_simple_script_', params_glob%nparts)
    end subroutine qsys_cleanup

    !>  Writes the JOB_FINISHED_* file to mark end of computing unit job
    subroutine qsys_job_finished(  source )
        use simple_parameters, only: params_glob
        ! generation of this file marks completion of the partition
        ! this file is empty 4 now but may contain run stats etc.
        character(len=*), intent(in) :: source
        if( params_glob%l_distr_exec .or. params_glob%stream.eq.'yes' )then
            call simple_touch('JOB_FINISHED_'//int2str_pad(params_glob%part,params_glob%numlen), errmsg="qsys_job_finished")
        endif
    end subroutine qsys_job_finished

    !>  returns when the inputted file exists in cwd
    subroutine qsys_watcher_1( fname, wtime )
        character(len=*),  intent(in) :: fname
        integer, optional, intent(in) :: wtime
        integer :: wwtime
        logical :: there
        wwtime = SHORTTIME
        if( present(wtime) ) wwtime = wtime
        do
            there = file_exists(trim(fname))
            if( there ) exit
            call simple_sleep(wwtime)
        end do
    end subroutine qsys_watcher_1

    !>  returns when the inputted files exist in cwd
    subroutine qsys_watcher_2( fnames, wtime )
        character(len=STDLEN), intent(in)    :: fnames(:)
        integer, optional,     intent(in)    :: wtime
        integer, parameter   :: MAXITS=20000
        integer              :: wwtime, nfiles, ifile, i
        logical              :: doreturn, fexists
        wwtime = SHORTTIME
        if( present(wtime) ) wwtime = wtime
        nfiles = size(fnames)
        do i=1,MAXITS
            doreturn = .true.
            do ifile=1,nfiles
                fexists = file_exists(trim(fnames(ifile)))
                if( .not. fexists ) doreturn = .false.
            end do
            if( doreturn )then
                return
            else
                call simple_sleep(wwtime)
            endif
        end do
    end subroutine qsys_watcher_2

end module simple_qsys_funs
