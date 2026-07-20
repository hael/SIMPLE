!@descr: analysis of class averages
module simple_commanders_sieve
use simple_commanders_api
use simple_rec_list,         only: rec_list
use simple_ptcl_sieve,       only: ptcl_sieve, DEFAULT_COARSE_POP_THRESHOLD, DEFAULT_FINE_POP_THRESHOLD, &
                                      DEFAULT_COARSE_BOX, DEFAULT_FINE_BOX, DEFAULT_COARSE_NSAMPLE, DEFAULT_FINE_NSAMPLE, &
                                      DEFAULT_LPSTART, DEFAULT_COARSE_LP, DEFAULT_FINE_LP, DEFAULT_NCLS
use simple_ptcl_sieve_utils, only: generate_sieve_projects
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_sieve_ptcls
  contains
    procedure :: execute      => exec_sieve_ptcls
end type commander_sieve_ptcls

contains

    subroutine exec_sieve_ptcls( self, cline )
        class(commander_sieve_ptcls), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(ptcl_sieve)                            :: sieve
        type(parameters)                            :: params
        type(sp_project)                            :: spproj
        type(rec_list)                              :: project_list
        type(string)                                :: tmp_projfile
        logical                                     :: l_once
        integer                                     :: nchunks, nmics, nstks, nptcls
        ! set defaults
        if( .not. cline%defined('mkdir')          ) call cline%set('mkdir',                                 'yes')
        if( .not. cline%defined('nmics')          ) call cline%set('nmics',                                    50)
        if( .not. cline%defined('nptcls_coarse')  ) call cline%set('nptcls_coarse',  DEFAULT_COARSE_POP_THRESHOLD)
        if( .not. cline%defined('nptcls_fine')    ) call cline%set('nptcls_fine',      DEFAULT_FINE_POP_THRESHOLD)
        if( .not. cline%defined('maxnchunks')     ) call cline%set('maxnchunks',                                0)
        if( .not. cline%defined('lpstart')        ) call cline%set('lpstart',                     DEFAULT_LPSTART)
        if( .not. cline%defined('lpstop_coarse')  ) call cline%set('lpstop_coarse',             DEFAULT_COARSE_LP)
        if( .not. cline%defined('lpstop_fine')    ) call cline%set('lpstop_fine',                 DEFAULT_FINE_LP)
        if( .not. cline%defined('box_coarse')     ) call cline%set('box_coarse',               DEFAULT_COARSE_BOX)
        if( .not. cline%defined('box_fine')       ) call cline%set('box_fine',                   DEFAULT_FINE_BOX)
        if( .not. cline%defined('nsample_coarse') ) call cline%set('nsample_coarse',       DEFAULT_COARSE_NSAMPLE)
        if( .not. cline%defined('nsample_fine')   ) call cline%set('nsample_fine',           DEFAULT_FINE_NSAMPLE)
        if( .not. cline%defined('ncls_coarse')    ) call cline%set('ncls_coarse',                    DEFAULT_NCLS)
        if( .not. cline%defined('ncls_fine')      ) call cline%set('ncls_fine',                      DEFAULT_NCLS)
        if( .not. cline%defined('use_model')      ) call cline%set('use_model',                             'yes')
        if( .not. cline%defined('single_pass')    ) call cline%set('single_pass',                            'no')
        ! set up output directories
        ! parse parameters
        call params%new(cline)
        ! read project
        call spproj%read(params%projfile)
        ! --- sanity checks ---
        nmics  = spproj%os_mic%get_noris()
        nstks  = spproj%os_stk%get_noris()
        nptcls = spproj%os_ptcl2D%get_noris()
        if( nmics == 0     ) THROW_HARD('No micrographs found in project file: '//params%projfile%to_char())
        if( nstks == 0     ) THROW_HARD('No stacks found in project file: '//params%projfile%to_char())
        if( nptcls == 0    ) THROW_HARD('No particles found in project file: '//params%projfile%to_char())
        if( nstks /= nmics ) THROW_HARD('Inconsistent # of stacks and micrographs, use prune_project.')
        ! --- partition dataset into per-chunk sub-projects ---
        call generate_sieve_projects(project_list, params, spproj, nchunks)
        call spproj%kill()
        if( nchunks == 0 ) call simple_end('**** NO CHUNKS GENERATED — INSUFFICIENT PARTICLES OR MICROGRAPHS ****')
        ! --- drive multi-tier classification loop until all tiers are complete ---
        call simple_mkdir(PATH_HERE // DIR_STREAM_COMPLETED)
        call sieve%new(params, string(PATH_HERE // DIR_STREAM_COMPLETED), pre_chunked=.true.)
        l_once = .true.
        do
            call sieve%cycle(project_list)
            if( l_once ) then
                call sieve%set_final_ingestion()
                l_once = .false.
            endif
            if( sieve%get_finished() ) exit
            call sleep(WAITTIME)
        end do
        ! --- combine results and finalise ---
        write(logfhandle,'(A)') ''
        write(logfhandle,'(A,I8,A,I8)') '>>> COARSE TIER PARTICLES ACCEPTED  : ', sieve%get_n_coarse_accepted_ptcls(), ' / ', &
            sieve%get_n_coarse_accepted_ptcls() + sieve%get_n_coarse_rejected_ptcls()
        write(logfhandle,'(A,I8,A,I8)') '>>> COARSE TIER PARTICLES REJECTED  : ', sieve%get_n_coarse_rejected_ptcls(), ' / ', &
            sieve%get_n_coarse_accepted_ptcls() + sieve%get_n_coarse_rejected_ptcls()
        write(logfhandle,'(A)') '' 
        write(logfhandle,'(A,I8,A,I8)') '>>> FINE TIER PARTICLES ACCEPTED    : ', sieve%get_n_fine_accepted_ptcls(), ' / ', &
            sieve%get_n_coarse_accepted_ptcls()
        write(logfhandle,'(A,I8,A,I8)') '>>> FINE TIER PARTICLES REJECTED    : ', sieve%get_n_coarse_accepted_ptcls() - sieve%get_n_fine_accepted_ptcls(), ' / ', &
            sieve%get_n_coarse_accepted_ptcls()
        write(logfhandle,'(A)') ''
        write(logfhandle,'(A,I8,A,I8)') '>>> TOTAL NUMBER PARTICLES ACCEPTED : ', sieve%get_n_accepted_ptcls(), ' / ', sieve%get_n_total_particles()
        write(logfhandle,'(A,I8,A,I8)') '>>> TOTAL NUMBER PARTICLES REJECTED : ', sieve%get_n_rejected_ptcls(), ' / ', sieve%get_n_total_particles()
        write(logfhandle,'(A)') ''
        write(logfhandle,'(A)') '>>> ALL CHUNKS PROCESSED, COMBINING RESULTS...'
        tmp_projfile = add2fbody(params%projfile, METADATA_EXT, '_tmp_sieve')
        if( file_exists(tmp_projfile) ) call del_file(tmp_projfile)
        call sieve%combine_completed_chunks(tmp_projfile)
        if( .not. file_exists(tmp_projfile) ) THROW_HARD('Failed to create combined sieve project: '//tmp_projfile%to_char())
        call simple_rename(tmp_projfile, params%projfile)
        call sieve%kill()
        call qsys_cleanup(params)
        call simple_end('**** SIMPLE_SIEVE_PTCLS NORMAL STOP ****')
    end subroutine exec_sieve_ptcls

end module simple_commanders_sieve
