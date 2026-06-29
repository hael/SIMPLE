!@descr: batch-processing manager - functions
module simple_qsys_funs
use simple_core_module_api
use simple_parameters, only: parameters
implicit none

interface qsys_watcher
    module procedure qsys_watcher_1
    module procedure qsys_watcher_2
end interface

integer, parameter :: SHORTTIME = 3

contains

    subroutine qsys_cleanup( params, keep2D  )
        class(parameters), intent(in) :: params
        logical, optional, intent(in) :: keep2D
        integer, parameter :: NUMLEN_STATE = 2, NUMLEN_ITER = 3
        logical :: l_keep2D
        integer :: nparts, nstates
        l_keep2D = .true.
        if(present(keep2D)) l_keep2D = keep2D
        nparts  = params%nparts
        nstates = params%nstates
        ! single files
        call del_file('FOO')
        call del_file('fort.0')
        call del_file('VOLASSEMBLE_FINISHED')
        call del_file('CAVGASSEMBLE_FINISHED')
        call del_file('SYMSRCH')
        call del_file(SIMPLE_SUBPROC_OUT)
        ! state numbered files
        call del_files('VOLASSEMBLE_FINISHED_STATE', nstates, numlen=NUMLEN_STATE)
        call del_files('simple_script_state',        nstates, numlen=NUMLEN_STATE)
        ! part numbered files
        call del_files('OUT',                  nparts)
        call del_files('algndoc_',             nparts, ext=trim(METADATA_EXT))
        call del_files('algndoc_cavgs_',       nparts, ext=trim(METADATA_EXT))
        call del_files(JOB_FINISHED_FBODY,     nparts)
        call del_files(SUBPROJECT_JOB_FINISHED_FBODY, nparts)
        call del_files('distr_simple_script_', nparts)
        ! optionally deletes 2D analysis temporary files
        if( .not.l_keep2D )then
            call del_files('cavgs_even_part',     nparts, ext=MRC_EXT)
            call del_files('cavgs_odd_part',      nparts, ext=MRC_EXT)
            call del_files('ctfsqsums_even_part', nparts, ext=MRC_EXT)
            call del_files('ctfsqsums_odd_part',  nparts, ext=MRC_EXT)
        endif
        ! flush the log filehandle to avoid delayed printing
        call flush(logfhandle)
    end subroutine qsys_cleanup

    !> Declare completion of one partition job.
    !! File-backed qsys backends retain the historical JOB_FINISHED_* sentinel
    !! contract.  Coarray completion is synchronized by the coarray dispatch
    !! driver, so this worker-side declaration deliberately avoids filesystem
    !! completion state for qsys=coarray.
    subroutine qsys_declare_part_finished( params, source )
        class(parameters), intent(in) :: params
        class(string),     intent(in) :: source
        logical :: l_stream
        l_stream = trim(params%stream) .eq. 'yes'
        if( trim(params%qsys_name) == 'coarray' ) return
        if( params%l_distr_worker .or. l_stream )then
            call simple_touch(JOB_FINISHED_FBODY//int2str_pad(params%part,params%numlen))
        endif
    end subroutine qsys_declare_part_finished

    !> Backward-compatible spelling for partition completion declaration.
    subroutine qsys_job_finished( params, source )
        class(parameters), intent(in) :: params
        class(string),     intent(in) :: source
        call qsys_declare_part_finished(params, source)
    end subroutine qsys_job_finished

    subroutine qsys_subproject_job_finished( params, source )
        class(parameters), intent(in) :: params
        class(string),     intent(in) :: source
        call simple_touch(SUBPROJECT_JOB_FINISHED_FBODY//int2str_pad(params%part,params%numlen))
    end subroutine qsys_subproject_job_finished

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

    !> Diagnostic watcher with periodic missing-file reporting and timeout hard-stop.
    !! Orthogonal to qsys_watcher so existing behavior elsewhere is untouched.
    subroutine qsys_watcher_diag( fnames, wtime, maxits, log_every )
        class(string),     intent(in) :: fnames(:)
        integer, optional, intent(in) :: wtime, maxits, log_every
        integer :: wwtime, mmaxits, llog_every
        integer :: nfiles, ifile, i, nmissing
        logical :: doreturn
        wwtime     = SHORTTIME
        mmaxits    = 20000
        llog_every = 20
        if( present(wtime) )     wwtime     = max(1, wtime)
        if( present(maxits) )    mmaxits    = max(1, maxits)
        if( present(log_every) ) llog_every = max(1, log_every)
        nfiles = size(fnames)
        do i=1,mmaxits
            doreturn = .true.
            nmissing = 0
            do ifile=1,nfiles
                if( .not. file_exists(fnames(ifile)) )then
                    doreturn = .false.
                    nmissing = nmissing + 1
                endif
            enddo
            if( doreturn ) return
            if( mod(i, llog_every) == 0 )then
                write(logfhandle,'(A,I0,A,I0,A,I0,A)') '>>> qsys_watcher_diag waiting: ', nmissing, &
                    ' of ', nfiles, ' files missing (elapsed ', i * wwtime, ' s)'
                do ifile=1,nfiles
                    if( .not. file_exists(fnames(ifile)) ) write(logfhandle,'(A,A)') &
                        '>>>   missing: ', fnames(ifile)%to_char()
                enddo
                call flush(logfhandle)
            endif
            call sleep(wwtime)
        enddo
        write(logfhandle,'(A)') '>>> qsys_watcher_diag timeout; remaining missing files:'
        do ifile=1,nfiles
            if( .not. file_exists(fnames(ifile)) ) write(logfhandle,'(A,A)') '>>>   missing: ', fnames(ifile)%to_char()
        enddo
        call flush(logfhandle)
        call simple_exception('qsys_watcher_diag timeout waiting for completion markers', 'simple_qsys_funs.f90', 0)
    end subroutine qsys_watcher_diag

end module simple_qsys_funs
