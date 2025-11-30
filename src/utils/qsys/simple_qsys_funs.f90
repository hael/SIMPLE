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

    subroutine qsys_cleanup( keep2D, nparts  )
        use simple_parameters, only: params_glob
        logical, optional, intent(in) :: keep2D
        integer, optional, intent(in) :: nparts
        integer, parameter :: NUMLEN_STATE = 2, NUMLEN_ITER = 3
        logical :: l_keep2D
        integer :: np
        l_keep2D = .true.
        if(present(keep2D)) l_keep2D = keep2D
        np = params_glob%nparts
        if( present(nparts) ) np = nparts
        ! single files
        call del_file('FOO')
        call del_file('fort.0')
        call del_file('VOLASSEMBLE_FINISHED')
        call del_file('CAVGASSEMBLE_FINISHED')
        call del_file('SYMSRCH')
        call del_file(SIMPLE_SUBPROC_OUT)
        ! state numbered files
        call del_files('VOLASSEMBLE_FINISHED_STATE', params_glob%nstates, numlen=NUMLEN_STATE)
        call del_files('simple_script_state',        params_glob%nstates, numlen=NUMLEN_STATE)
        ! part numbered files
        call del_files('OUT',                  np)
        call del_files('algndoc_',             np, ext=trim(METADATA_EXT))
        call del_files('algndoc_cavgs_',       np, ext=trim(METADATA_EXT))
        call del_files(JOB_FINISHED_FBODY,     np)
        call del_files('distr_simple_script_', np)
        ! optionally deletes 2D analysis temporary files
        if( .not.l_keep2D )then
            call del_files('cavgs_even_part',     np, ext=params_glob%ext%to_char())
            call del_files('cavgs_odd_part',      np, ext=params_glob%ext%to_char())
            call del_files('ctfsqsums_even_part', np, ext=params_glob%ext%to_char())
            call del_files('ctfsqsums_odd_part',  np, ext=params_glob%ext%to_char())
        endif
        ! flush the log filehandle to avoid delayed printing
        call flush(logfhandle)
    end subroutine qsys_cleanup

    !>  Writes the JOB_FINISHED_* file to mark end of computing unit job
    subroutine qsys_job_finished(  source )
        use simple_parameters, only: params_glob
        ! generation of this file marks completion of the partition
        ! this file is empty 4 now but may contain run stats etc.
        class(string), intent(in) :: source
        if( params_glob%l_distr_exec .or. params_glob%stream.eq.'yes' )then
            call simple_touch(JOB_FINISHED_FBODY//int2str_pad(params_glob%part,params_glob%numlen))
        endif
    end subroutine qsys_job_finished

    !>  returns when the inputted file exists in cwd
    subroutine qsys_watcher_1( fname, wtime )
        class(string),     intent(in) :: fname
        integer, optional, intent(in) :: wtime
        integer :: wwtime
        logical :: there
        wwtime = SHORTTIME
        if( present(wtime) ) wwtime = wtime
        do
            there = file_exists(fname)
            if( there ) exit
            call sleep(wwtime)
        end do
    end subroutine qsys_watcher_1

    !>  returns when the inputted files exist in cwd
    subroutine qsys_watcher_2( fnames, wtime )
        class(string),     intent(in) :: fnames(:)
        integer, optional, intent(in) :: wtime
        integer, parameter   :: MAXITS=20000
        integer              :: wwtime, nfiles, ifile, i
        logical              :: doreturn, fexists
        wwtime = SHORTTIME
        if( present(wtime) ) wwtime = wtime
        nfiles = size(fnames)
        do i=1,MAXITS
            doreturn = .true.
            do ifile=1,nfiles
                fexists = file_exists(fnames(ifile))
                if( .not. fexists ) doreturn = .false.
            end do
            if( doreturn )then
                return
            else
                call sleep(wwtime)
            endif
        end do
    end subroutine qsys_watcher_2

end module simple_qsys_funs
