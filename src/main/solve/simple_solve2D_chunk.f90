!@descr: one subset of solve2D_chunks: its project, queue environment, command line and solve2D job
module simple_solve2D_chunk
use simple_core_module_api
use simple_defs_environment, only: SIMPLE_STREAM_CHUNK_PARTITION
use simple_cmdline,      only: cmdline
use simple_parameters,   only: parameters
use simple_qsys_env,     only: qsys_env
use simple_sp_project,   only: sp_project
use simple_rec_list,     only: rec_list
implicit none

public :: solve2D_chunk
private
#include "simple_local_flags.inc"

! Type to handle a single chunk
type solve2D_chunk
    private
    class(parameters), pointer :: p_ptr => null()       ! parameters pointer
    ! allocatable (compile-time policy), allocated in init_chunk and released in kill
    type(sp_project), allocatable :: spproj             ! master project
    type(qsys_env),   allocatable :: qenv               ! submission handler
    type(cmdline),    allocatable :: cline              ! command line
    type(string), allocatable  :: orig_stks(:)          ! list of stacks
    type(string)               :: path, projfile_out    ! physical location
    integer                    :: id                    ! unique id
    integer                    :: it                    ! # of iterations performed
    integer                    :: nmics                 ! # of micrographs
    integer                    :: nptcls                ! # number of particles
    logical                    :: toanalyze2D = .true.  ! whether to perform 2D analysis or only calculate sigmas2
    logical                    :: converged   = .false. ! whether 2D analysis is over
    logical                    :: available   = .true.  ! has been initialized but no 2D analysis peformed
  contains
    procedure          :: init_chunk
    procedure, private :: copy
    procedure, private :: assign
    generic            :: assignment(=) => assign
    procedure          :: generate
    procedure          :: is_available
    procedure          :: get_projfile_fname
    procedure          :: calc_sigma2
    procedure          :: analyze2D
    procedure, private :: gen_final_cavgs
    procedure          :: display_iter
    procedure          :: has_converged
    procedure          :: is_finished
    procedure          :: print_info
    procedure          :: kill
end type solve2D_chunk

logical, parameter :: DEBUG_HERE = .false.

contains

    ! Utility
    subroutine debug_print( string )
        character(len=*), intent(in) :: string
        if( DEBUG_HERE )then
            write(logfhandle,*) trim(string)
            call flush(logfhandle)
        endif
    end subroutine debug_print

    !>  Instantiator
    subroutine init_chunk( self, params, cline, id, master_spproj )
        class(solve2D_chunk),       intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        class(cmdline),            intent(in)    :: cline
        integer,                   intent(in)    :: id
        class(sp_project),         intent(in)    :: master_spproj
        type(string)          :: exec
        character(len=STDLEN) :: chunk_part_env
        integer               :: envlen
        call debug_print('in chunk%init '//int2str(id))
        if( .not. allocated(self%spproj) ) allocate(self%spproj)
        if( .not. allocated(self%qenv)   ) allocate(self%qenv)
        if( .not. allocated(self%cline)  ) allocate(self%cline)
        call self%spproj%kill
        self%p_ptr  => params
        self%id     = id
        self%it     = 0
        self%nmics  = 0 ! # of micrographs & stacks in chunk
        self%nptcls = 0
        self%path   = './'//DIR_CHUNK//int2str(id)//'/'
        self%cline  = cline
        if( .not.cline%defined('box') )      THROW_HARD('Missing BOX on chunk command-line')
        if( .not.cline%defined('box_crop') ) THROW_HARD('Missing BOX_CROP on chunk command-line')
        if( .not.cline%defined('smpd') )     THROW_HARD('Missing SMPD on chunk command-line')
        if( .not.cline%defined('smpd_crop') )THROW_HARD('Missing SMPD_CROP on chunk command-line')
        self%projfile_out    = ''
        self%spproj%projinfo = master_spproj%projinfo
        self%spproj%compenv  = master_spproj%compenv
        exec = 'simple_exec'
        if( self%p_ptr%nparts_chunk == 1 )then
            ! shared memory
            call self%qenv%new(self%p_ptr, self%p_ptr%nparts_chunk, exec_bin=exec, qsys_nthr=self%p_ptr%nthr2D)
            call self%spproj%compenv%set(1,'qsys_name','local')
        else
            ! we need to override the qsys_name for non local distributed execution
            call self%qenv%new(self%p_ptr, self%p_ptr%nparts_chunk, exec_bin=exec, qsys_name=string('local'))
            call get_environment_variable(SIMPLE_STREAM_CHUNK_PARTITION, chunk_part_env, envlen)
            if(envlen > 0) call self%spproj%compenv%set(1,'qsys_partition', trim(chunk_part_env))
        endif
        call self%spproj%projinfo%delete_entry('projname')
        call self%spproj%projinfo%delete_entry('projfile')
        ! A chunk is an independent particle layout and must never inherit the
        ! master's canonical-state identity.
        call self%spproj%projinfo%delete_entry('sigma2_state')
        if( allocated(self%orig_stks) )then
            call self%orig_stks%kill
            deallocate(self%orig_stks)
        endif
        self%toanalyze2D = .true.
        self%converged   = .false.
        self%available   = .true.
        call debug_print('end chunk%init '//int2str(id))
    end subroutine init_chunk

    subroutine copy( dest, src )
        class(solve2D_chunk), intent(inout) :: dest
        class(solve2D_chunk), intent(in)    :: src
        call dest%kill
        ! sp_project and cmdline assign through a defined assignment: allocated first
        if( allocated(src%spproj) )then
            allocate(dest%spproj)
            dest%spproj = src%spproj
        endif
        if( allocated(src%qenv) ) dest%qenv = src%qenv
        if( allocated(src%cline) )then
            allocate(dest%cline)
            dest%cline = src%cline
        endif
        dest%orig_stks    = src%orig_stks(:)
        dest%path         = src%path
        dest%projfile_out = src%projfile_out
        dest%id           = src%id
        dest%it           = src%it
        dest%nmics        = src%nmics
        dest%nptcls       = src%nptcls
        dest%toanalyze2D  = src%toanalyze2D
        dest%converged    = src%converged
        dest%available    = src%available
        dest%p_ptr        => src%p_ptr
    end subroutine copy

    !>  \brief  assign, polymorphic assignment (=)
    subroutine assign( selfout, selfin )
        class(solve2D_chunk), intent(inout) :: selfout
        class(solve2D_chunk), intent(in)    :: selfin
        call selfout%copy(selfin)
    end subroutine assign

    subroutine generate( self, project_list )
        class(solve2D_chunk), intent(inout) :: self
        class(rec_list),     intent(inout) :: project_list
        integer :: istk, sz
        if( .not.self%available ) THROW_HARD('chunk unavailable; chunk%generate')
        sz = project_list%size()
        if( sz == 0 ) THROW_HARD('# ptcls == 0; chunk%generate')
        call debug_print('in chunk%generate '//int2str(self%id)//' '//int2str(sz))
        call self%spproj%projrecords2proj(project_list)
        self%nmics  = self%spproj%os_mic%get_noris()
        self%nptcls = self%spproj%os_ptcl2D%get_noris()
        allocate(self%orig_stks(self%nmics))
        do istk = 1,self%nmics
            self%orig_stks(istk) = self%spproj%get_stkname(istk)
        enddo
        call debug_print('end chunk%generate_2 '//int2str(self%id))
    end subroutine generate

    elemental function is_available( self ) result( avail )
        class(solve2D_chunk), intent(in) :: self
        logical :: avail
        avail = self%available
    end function is_available

    function get_projfile_fname( self )result( fname )
        class(solve2D_chunk), intent(in) :: self
        type(string) :: fname
        fname = self%path//self%projfile_out
    end function get_projfile_fname

    ! Initiates 2D analysis
    subroutine analyze2D( self, makecavgs )
        class(solve2D_chunk), intent(inout) :: self
        logical,   optional, intent(in)    :: makecavgs
        type(string),  allocatable :: bins(:)
        type(cmdline), allocatable :: clines(:)
        type(cmdline) :: cline_pspec
        type(string)  :: cwd
        integer       :: nptcls_sel, nclines
        logical       :: l_makecavgs
        call debug_print('in chunk%analyze2D '//int2str(self%id))
        if( .not.self%available ) return
        if( self%nptcls == 0 ) return
        l_makecavgs = .false.
        if( present(makecavgs) )l_makecavgs = makecavgs
        call simple_mkdir(self%path)
        call simple_chdir(self%path)
        call simple_getcwd(cwd)
        CWD_GLOB    = cwd%to_char()
        self%projfile_out = CHUNK_PROJNAME//METADATA_EXT
        call simple_mkdir(STDERROUT_DIR)
        nptcls_sel = self%spproj%os_ptcl2D%get_noris(consider_state=.true.)
        nclines = 1
        allocate(clines(nclines))
        call self%cline%set('projfile', self%projfile_out)
        call self%cline%set('projname', CHUNK_PROJNAME)
        call self%spproj%update_projinfo(self%cline)
        call self%spproj%write()
        ! 2D analysis
        clines(nclines) = self%cline
        if( l_makecavgs )then
            ! optional postprocessing
            call self%gen_final_cavgs( clines )
            ! submission
            nclines = nclines+1
            allocate(bins(nclines))
            bins(1:nclines-1) = self%qenv%get_exec_bin()
            bins(nclines)     = 'simple_exec'
            call self%qenv%exec_simple_prgs_in_queue_async(clines, string('./distr_chunk2D'), string('simple_log_chunk2d'), exec_bins=bins(:))
        else
            ! submission
            call self%qenv%exec_simple_prgs_in_queue_async(clines, string('./distr_chunk2D'), string('simple_log_chunk2d'))
        endif
        call simple_chdir('..')
        call simple_getcwd(cwd)
        CWD_GLOB = cwd%to_char()
        ! cleanup
        call self%spproj%kill
        call cline_pspec%kill
        call clines(:)%kill
        deallocate(clines)
        call self%qenv%kill
        ! chunk is now busy
        self%available = .false.
        self%converged = .false.
        write(logfhandle,'(A,I6,A,I6,A)')'>>> CHUNK ',self%id,' INITIATED 2D ANALYSIS WITH ',nptcls_sel,' PARTICLES'
        call debug_print('end chunk%analyze2D')
    end subroutine analyze2D

    ! To calculate noise power estimates only
    subroutine calc_sigma2( self, cline_analyze2D, need_sigma )
        class(solve2D_chunk), intent(inout) :: self
        class(cmdline),      intent(in)    :: cline_analyze2D
        logical,             intent(in)    :: need_sigma
        type(cmdline) :: cline_pspec
        type(string)  :: cwd
        type(string)  :: exec
        integer :: nptcls_sel
        call debug_print('in calc_sigma2 '//int2str(self%id))
        if( .not.self%available ) return
        if( self%nptcls == 0 ) return
        call simple_mkdir(self%path)
        call simple_chdir(self%path)
        call simple_getcwd(cwd)
        CWD_GLOB    = cwd%to_char()
        self%projfile_out = CHUNK_PROJNAME//METADATA_EXT
        call simple_mkdir(STDERROUT_DIR)
        nptcls_sel = self%spproj%os_ptcl2D%get_noris(consider_state=.true.)
        call cline_pspec%set('prg',      'calc_pspec')
        call cline_pspec%set('oritype',  'ptcl2D')
        call cline_pspec%set('nthr',     cline_analyze2D%get_iarg('nthr'))
        call cline_pspec%set('mkdir',    'yes')
        call cline_pspec%set('nparts',   1)
        if( self%p_ptr%nparts_chunk > 1 ) call cline_pspec%set('nparts',self%p_ptr%nparts_chunk)
        call cline_pspec%set('projfile', self%projfile_out)
        call cline_pspec%set('projname', CHUNK_PROJNAME)
        call self%spproj%update_projinfo(cline_pspec)
        call self%spproj%write()
        self%available  = .false.
        self%toanalyze2D = .false. ! not to be classified
        self%it         = 1
        if( need_sigma )then
            ! making sure the executable is *always* simple_private_exec
            exec = self%qenv%get_exec_bin()
            call exec%substr_replace('/simple_exec','/simple_private_exec',one=.true.,back=.true.)
            ! submission
            self%converged  = .false.
            call self%qenv%exec_simple_prg_in_queue_async(cline_pspec, string('./distr_chunk2D'), string('simple_log_chunk2d'), exec_bin=exec)
            write(logfhandle,'(A,I6,A,I6,A)')'>>> CHUNK ',self%id,' INITIATED SIGMA2 CALCULATION WITH ',nptcls_sel,' PARTICLES'
        else
            self%converged  = .true. ! nothing to calculate
        endif
        call simple_chdir('..')
        call simple_getcwd(cwd)
        CWD_GLOB = cwd%to_char()
        ! cleanup
        call self%spproj%kill
        call cline_pspec%kill
        call self%qenv%kill
        call debug_print('end calc_sigma2')
    end subroutine calc_sigma2

    ! classes generation at original sampling
    subroutine gen_final_cavgs( self, clines )
        class(solve2D_chunk),         intent(in)    :: self
        type(cmdline),  allocatable, intent(inout) :: clines(:)
        type(cmdline),    allocatable :: tmp(:)
        type(cmdline)                 :: cline_make_cavgs
        type(string) :: finalcavgs
        integer :: n, i
        if( .not.allocated(clines) ) THROW_HARD('Fatal error gen_final_cavgs')
        finalcavgs = 'final_cavgs.mrc'
        call cline_make_cavgs%set('prg',        'make_cavgs')
        call cline_make_cavgs%set('mkdir',      'no')
        call cline_make_cavgs%set('refs',       finalcavgs)
        call cline_make_cavgs%set('ncls',       self%cline%get_iarg('ncls'))
        call cline_make_cavgs%set('mskdiam',    self%p_ptr%mskdiam)
        call cline_make_cavgs%set('async',      'yes')
        call cline_make_cavgs%set('nthr',       self%p_ptr%nthr2D)
        if( self%cline%defined('nparts') )then
            call cline_make_cavgs%set('nparts', self%cline%get_iarg('nparts'))
        else
            call cline_make_cavgs%set('nparts', 1)
        endif
        n = size(clines)
        allocate(tmp(1:n),source=clines(1:n))
        call clines(:)%kill; deallocate(clines); allocate(clines(n+1))
        do i = 1, n
            clines(i) = tmp(i)
        end do
        clines(n+1) = cline_make_cavgs
        call tmp(:)%kill; call cline_make_cavgs%kill; deallocate(tmp)
    end subroutine gen_final_cavgs

    ! get & display convergence stats
    subroutine display_iter( self )
        class(solve2D_chunk), intent(inout) :: self
        type(oris)   :: os
        type(string) :: fname
        real         :: mi_class,frac,corr
        integer      :: it
        fname = self%path//STATS_FILE
        call debug_print('in chunk%display_iter '//int2str(self%id)//' '//fname%to_char())
        if( file_exists(fname) )then
            call os%new(1,is_ptcl=.false.)
            call os%read(fname)
            it = nint(os%get(1,'ITERATION'))
            if( it /= self%it )then
                self%it  = it
                mi_class = os%get(1,'CLASS_OVERLAP')
                frac     = os%get(1,'SEARCH_SPACE_SCANNED')
                corr     = os%get(1,'SCORE')
                write(logfhandle,'(A,I6,A,I6,A,F7.3,A,F7.3,A,F7.3)')'>>> CHUNK ',self%id,' ITERATION ',self%it,&
                    &'; CLASS OVERLAP: ',mi_class,'; SEARCH SPACE SCANNED: ',frac,'; SCORE: ',corr
            endif
            call os%kill
        endif
        call debug_print('end chunk%display_iter '//int2str(self%id))
    end subroutine display_iter

    ! Check for convergence of 2D analysis is complete
    logical function has_converged( self )
        class(solve2D_chunk), intent(inout) :: self
        type(string) :: str_prg
        if( .not.self%converged )then
            if( self%toanalyze2D )then
                str_prg = self%cline%get_carg('prg')
                self%converged = file_exists(self%path//SOLVE2D_FINISHED)
                call str_prg%kill
            else
                self%converged = file_exists(self%path//CALCPSPEC_FINISHED)
            endif
        endif
        has_converged  = self%converged
    end function has_converged

    ! Check for convergence of 2D analysis is complete
    elemental logical function is_finished( self )
        class(solve2D_chunk), intent(in) :: self
        is_finished  = self%converged
    end function is_finished

    ! For debugging
    subroutine print_info( self )
        class(solve2D_chunk), intent(inout) :: self
        print *,'self%id           : ',self%id
        print *,'self%path         : ',self%path%to_char()
        print *,'self%it           : ',self%it
        print *,'self%nmics        : ',self%nmics
        print *,'self%nptcls       : ',self%nptcls
        print *,'self%converged    : ',self%converged
        print *,'self%available    : ',self%available
    end subroutine print_info

    subroutine kill( self )
        class(solve2D_chunk), intent(inout) :: self
        self%id        = 0
        self%it        = 0
        if( allocated(self%spproj) )then
            call self%spproj%kill
            deallocate(self%spproj)
        endif
        if( allocated(self%qenv) )then
            call self%qenv%kill
            deallocate(self%qenv)
        endif
        if( allocated(self%cline) )then
            call self%cline%kill
            deallocate(self%cline)
        endif
        self%nmics     = 0
        self%nptcls    = 0
        self%path      = ''
        self%projfile_out = ''
        if( allocated(self%orig_stks) )then
            call self%orig_stks%kill
            deallocate(self%orig_stks)
        endif
        self%toanalyze2D = .true.
        self%converged   = .false.
        self%available  = .false.
        self%p_ptr      => null()
    end subroutine kill

end module simple_solve2D_chunk
