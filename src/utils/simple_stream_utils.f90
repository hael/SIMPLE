module simple_stream_utils
include 'simple_lib.f08'
use simple_builder,        only: builder
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_parameters,     only: params_glob
use simple_sp_project,     only: sp_project
use simple_qsys_env,       only: qsys_env
use simple_image,          only: image
use simple_qsys_funs
use simple_commander_cluster2D
use simple_timer
use simple_imgproc
implicit none

public :: stream_chunk, DIR_CHUNK
public :: projrecord, projrecords2proj, class_rejection
public :: procrecord, append_procrecord, kill_procrecords
private
#include "simple_local_flags.inc"

character(len=STDLEN), parameter   :: DIR_CHUNK          = 'chunk_'
character(len=STDLEN), parameter   :: PROJNAME_CHUNK     = 'chunk'
character(len=STDLEN), parameter   :: CLS_CHUNK_REJECTED = 'cls_rejected_chunks.mrc'
logical,               parameter   :: DEBUG_HERE         = .false.
integer                            :: ncls_rejected_glob = 0 ! counter of rejected classes

! Convenience type to hold information about individual project files
type projrecord
    character(len=:), allocatable :: projname               ! project file name
    integer                       :: micind     = 0         ! index of micrograph in project
    integer                       :: nptcls     = 0         ! # of particles
    integer                       :: nptcls_sel = 0         ! # of particles (state=1)
    logical                       :: included   = .false.   ! whether record has been imported
end type projrecord

! Convenience type to hold information about processes
type procrecord
    character(len=:), allocatable :: id                     ! unique ID
    character(len=:), allocatable :: folder                 ! location
    character(len=:), allocatable :: projfile               ! filename
    character(len=:), allocatable :: volume                 ! volume filename
    character(len=:), allocatable :: alnvolume              ! aligned volume filename
    logical                       :: submitted = .false.    ! process has been submitted (running)
    logical                       :: complete  = .false.    ! is complete
    logical                       :: included  = .false.    ! has been post-processed/analyzed
end type procrecord

! Type to handle a single chunk
type stream_chunk
    type(sp_project)                       :: spproj                ! master project
    type(qsys_env)                         :: qenv                  ! submission handler
    type(cmdline)                          :: cline
    character(len=LONGSTRLEN), allocatable :: orig_stks(:)          ! list of stacks
    character(len=LONGSTRLEN)              :: path, projfile_out    ! physical location
    integer                                :: id                    ! unique id
    integer                                :: it                    ! # of iterations performed
    integer                                :: nmics                 ! # of micrographs
    integer                                :: nptcls                ! # number of particles
    logical                                :: toanalyze2D = .true.  ! whether to perform 2D analysis or only calculate sigmas2
    logical                                :: converged  = .false.  ! whether 2D analysis is over
    logical                                :: available  = .true.   ! has been initialized but no 2D analysis peformed
  contains
    procedure          :: init
    procedure          :: copy
    procedure, private :: assign
    generic            :: assignment(=) => assign
    procedure          :: generate
    procedure          :: calc_sigma2
    procedure          :: analyze2D
    procedure          :: read
    procedure          :: split_sigmas_into
    procedure, private :: gen_final_cavgs
    procedure          :: remove_folder
    procedure          :: display_iter
    procedure          :: reject
    procedure          :: has_converged
    procedure          :: print_info
    procedure          :: terminate
    procedure          :: kill
end type stream_chunk

contains

    ! Utilities

    subroutine debug_print( string )
        character(len=*), intent(in) :: string
        if( DEBUG_HERE )then
            write(logfhandle,*) trim(string)
            call flush(logfhandle)
        endif
    end subroutine debug_print

    ! Chunk type related routines

    !>  Instantiator
    subroutine init( self, id, cline, master_spproj )
        class(stream_chunk), intent(inout) :: self
        class(cmdline),      intent(in)    :: cline
        integer,             intent(in)    :: id
        class(sp_project),   intent(in)    :: master_spproj
        character(len=STDLEN)              :: chunk_part_env, exec
        integer                            :: envlen
        call debug_print('in chunk%init '//int2str(id))
        call self%spproj%kill
        self%id     = id
        self%it     = 0
        self%nmics  = 0 ! # of micrographs & stacks in chunk
        self%nptcls = 0
        self%path   = './'//trim(DIR_CHUNK)//int2str(id)//'/'
        self%cline  = cline
        if( .not.cline%defined('box') )      THROW_HARD('Missing BOX on chunk command-line')
        if( .not.cline%defined('box_crop') ) THROW_HARD('Missing BOX_CROP on chunk command-line')
        if( .not.cline%defined('smpd') )     THROW_HARD('Missing SMPD on chunk command-line')
        if( .not.cline%defined('smpd_crop') )THROW_HARD('Missing SMPD_CROP on chunk command-line')
        self%projfile_out    = ''
        self%spproj%projinfo = master_spproj%projinfo
        self%spproj%compenv  = master_spproj%compenv
        if( str_has_substr(self%cline%get_carg('prg'), 'abinitio') )then
            exec = 'simple_exec'
        else
            exec = 'simple_private_exec'
        endif
        if( params_glob%nparts_chunk == 1 )then
            ! shared memory
            call self%qenv%new(params_glob%nparts_chunk, exec_bin=trim(exec), qsys_nthr=params_glob%nthr2D)
            call self%spproj%compenv%set(1,'qsys_name','local')
        else
            ! we need to override the qsys_name for non local distributed execution
            call self%qenv%new(params_glob%nparts_chunk, exec_bin=trim(exec), qsys_name='local')
            call get_environment_variable(SIMPLE_STREAM_CHUNK_PARTITION, chunk_part_env, envlen)
            if(envlen > 0) call self%spproj%compenv%set(1,'qsys_partition',trim(chunk_part_env))
        endif
        call self%spproj%projinfo%delete_entry('projname')
        call self%spproj%projinfo%delete_entry('projfile')
        if( allocated(self%orig_stks) ) deallocate(self%orig_stks)
        self%toanalyze2D = .true.
        self%converged  = .false.
        self%available  = .true.
        call debug_print('end chunk%init '//int2str(id))
    end subroutine init

    subroutine copy( dest, src )
        class(stream_chunk), intent(inout) :: dest
        class(stream_chunk), intent(in)    :: src
        call dest%kill
        dest%spproj       = src%spproj
        dest%qenv         = src%qenv
        dest%cline        = src%cline
        dest%orig_stks    = src%orig_stks(:)
        dest%path         = trim(src%path)
        dest%projfile_out = trim(src%projfile_out)
        dest%id           = src%id
        dest%it           = src%it
        dest%nmics        = src%nmics
        dest%nptcls       = src%nptcls
        dest%toanalyze2D   = src%toanalyze2D
        dest%converged    = src%converged
        dest%available    = src%available
    end subroutine copy

    !>  \brief  assign, polymorphic assignment (=)
    subroutine assign( selfout, selfin )
        class(stream_chunk), intent(inout) :: selfout
        class(stream_chunk), intent(in)    :: selfin
        call selfout%copy(selfin)
    end subroutine assign

    subroutine generate( self, micproj_records )
        class(stream_chunk), intent(inout) :: self
        type(projrecord),    intent(in)    :: micproj_records(:)
        integer :: istk
        if( .not.self%available ) THROW_HARD('chunk unavailable; chunk%generate')
        if( size(micproj_records(:)) == 0 ) THROW_HARD('# ptcls == 0; chunk%generate')
        call debug_print('in chunk%generate '//int2str(self%id)//' '//int2str(size(micproj_records(:))))
        call projrecords2proj(micproj_records(:), self%spproj)
        self%nmics  = self%spproj%os_mic%get_noris()
        self%nptcls = self%spproj%os_ptcl2D%get_noris()
        allocate(self%orig_stks(self%nmics))
        do istk = 1,self%nmics
            self%orig_stks(istk) = self%spproj%get_stkname(istk)
        enddo
        call debug_print('end chunk%generate_2 '//int2str(self%id))
    end subroutine generate

    ! Initiates 2D analysis
    subroutine analyze2D( self, calc_pspec, makecavgs )
        class(stream_chunk), intent(inout) :: self
        logical,             intent(in)    :: calc_pspec
        logical,   optional, intent(in)    :: makecavgs
        character(len=STDLEN), allocatable :: bins(:)
        type(cmdline),         allocatable :: clines(:)
        type(cmdline)              :: cline_pspec
        character(len=XLONGSTRLEN) :: cwd
        integer                    :: nptcls_sel, nclines
        logical                    :: l_makecavgs
        call debug_print('in chunk%analyze2D '//int2str(self%id))
        if( .not.self%available ) return
        if( self%nptcls == 0 ) return
        l_makecavgs = .false.
        if( present(makecavgs) )l_makecavgs = makecavgs
        call simple_mkdir(self%path)
        call chdir(self%path)
        call simple_getcwd(cwd)
        cwd_glob    = trim(cwd)
        self%projfile_out = trim(PROJNAME_CHUNK)//trim(METADATA_EXT)
        call simple_mkdir(STDERROUT_DIR)
        nptcls_sel = self%spproj%os_ptcl2D%get_noris(consider_state=.true.)
        nclines = 1
        if( str_has_substr(self%cline%get_carg('prg'), 'abinitio') )then
            allocate(clines(nclines))
        else
            if( calc_pspec ) nclines = nclines + 1
            allocate(clines(nclines))
            ! noise estimates
            if( calc_pspec )then
                call cline_pspec%set('prg',      'calc_pspec_distr')
                call cline_pspec%set('oritype',  'ptcl2D')
                call cline_pspec%set('projfile', self%projfile_out)
                call cline_pspec%set('nthr',     self%cline%get_iarg('nthr'))
                call cline_pspec%set('mkdir',    'yes')
                call cline_pspec%set('nparts',   1)
                if( params_glob%nparts_chunk > 1 ) call cline_pspec%set('nparts',params_glob%nparts_chunk)
                if( self%cline%defined('sigma_est') ) call cline_pspec%set('sigma_est', self%cline%get_carg('sigma_est'))
                clines(1) = cline_pspec
            endif
        endif
        call self%cline%set('projfile', self%projfile_out)
        call self%cline%set('projname', trim(PROJNAME_CHUNK))
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
            bins(1:nclines-1) = trim(self%qenv%get_exec_bin())
            bins(nclines)     = 'simple_exec'
            call self%qenv%exec_simple_prgs_in_queue_async(clines, './distr_chunk2D', 'simple_log_chunk2d', exec_bins=bins(:))
        else
            ! submission
            call self%qenv%exec_simple_prgs_in_queue_async(clines, './distr_chunk2D', 'simple_log_chunk2d')
        endif
        call chdir('..')
        call simple_getcwd(cwd_glob)
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
        class(stream_chunk), intent(inout) :: self
        class(cmdline),      intent(in)    :: cline_analyze2D
        logical,             intent(in)    :: need_sigma
        type(cmdline)                 :: cline_pspec
        character(len=XLONGSTRLEN)    :: cwd
        character(len=:), allocatable :: exec
        integer :: nptcls_sel
        call debug_print('in calc_sigma2 '//int2str(self%id))
        if( .not.self%available ) return
        if( self%nptcls == 0 ) return
        call simple_mkdir(self%path)
        call chdir(self%path)
        call simple_getcwd(cwd)
        cwd_glob    = trim(cwd)
        self%projfile_out = trim(PROJNAME_CHUNK)//trim(METADATA_EXT)
        call simple_mkdir(STDERROUT_DIR)
        nptcls_sel = self%spproj%os_ptcl2D%get_noris(consider_state=.true.)
        call cline_pspec%set('prg',      'calc_pspec_distr')
        call cline_pspec%set('oritype',  'ptcl2D')
        call cline_pspec%set('nthr',     cline_analyze2D%get_iarg('nthr'))
        call cline_pspec%set('mkdir',    'yes')
        call cline_pspec%set('nparts',   1)
        if( params_glob%nparts_chunk > 1 ) call cline_pspec%set('nparts',params_glob%nparts_chunk)
        call cline_pspec%set('projfile', self%projfile_out)
        call cline_pspec%set('projname', trim(PROJNAME_CHUNK))
        call self%spproj%update_projinfo(cline_pspec)
        call self%spproj%write()
        self%available  = .false.
        self%toanalyze2D = .false. ! not to be classified
        self%it         = 1
        if( need_sigma )then
            ! making sure the executable is *always* simple_private_exec
            exec = trim(self%qenv%get_exec_bin())
            call replace_substring(exec,'/simple_exec','/simple_private_exec',one=.true.,back=.true.)
            ! submission
            self%converged  = .false.
            call self%qenv%exec_simple_prg_in_queue_async(cline_pspec, './distr_chunk2D', 'simple_log_chunk2d', exec_bin=exec)
            write(logfhandle,'(A,I6,A,I6,A)')'>>> CHUNK ',self%id,' INITIATED SIGMA2 CALCULATION WITH ',nptcls_sel,' PARTICLES'
        else
            self%converged  = .true. ! nothing to calculate
        endif
        call chdir('..')
        call simple_getcwd(cwd_glob)
        ! cleanup
        call self%spproj%kill
        call cline_pspec%kill
        call self%qenv%kill
        call debug_print('end calc_sigma2')
    end subroutine calc_sigma2

    ! read project and deals with arrays used for fractional update
    subroutine read( self, box )
        class(stream_chunk), intent(inout) :: self
        integer,             intent(in)    :: box
        type(image)                :: img, avg
        character(len=XLONGSTRLEN) :: projfile
        call debug_print('in chunk%read '//int2str(self%id))
        if( .not.self%converged )THROW_HARD('cannot read chunk prior to convergence')
        ! doc & parameters
        projfile = trim(self%path)//trim(self%projfile_out)
        call self%spproj%read(projfile)
        ! classes, to account for nparts /= nparts_chunk
        if( self%toanalyze2D )then
            call img%new([box,box,1],1.0)
            call avg%new([box,box,1],1.0)
            call average_into(trim(self%path)//'cavgs_even_part')
            call average_into(trim(self%path)//'cavgs_odd_part')
            call average_into(trim(self%path)//'ctfsqsums_even_part')
            call average_into(trim(self%path)//'ctfsqsums_odd_part')
            call img%kill
            call avg%kill
        endif
        call debug_print('end chunk%read '//int2str(self%id))
        contains

            subroutine average_into(tmpl)
                character(len=*), intent(in) :: tmpl
                character(len=XLONGSTRLEN)   :: fname
                integer                      :: icls, ipart, numlen_chunk, iostat
                if( params_glob%nparts_chunk > 1  )then
                    numlen_chunk = len(int2str(params_glob%nparts_chunk)) ! as per parameters
                    call img%zero_and_flag_ft
                    do icls = 1,params_glob%ncls_start
                        call avg%zero_and_flag_ft
                        do ipart = 1,params_glob%nparts_chunk
                            fname = trim(tmpl)//int2str_pad(ipart,numlen_chunk)//trim(params_glob%ext)
                            call img%read(fname,icls)
                            call avg%add(img)
                        enddo
                        call avg%div(real(params_glob%nparts_chunk))
                        call avg%write(trim(tmpl)//trim(params_glob%ext),icls)
                    enddo
                else
                    fname = trim(tmpl)//'1'//trim(params_glob%ext)
                    iostat = rename(fname,trim(tmpl)//trim(params_glob%ext))
                endif
            end subroutine average_into

    end subroutine read

    ! split sigmas into individually named per stack documents
    subroutine split_sigmas_into( self, folder )
        use simple_euclid_sigma2, only: split_sigma2_into_groups, sigma2_star_from_iter
        class(stream_chunk),        intent(in) :: self
        character(len=*),           intent(in) :: folder
        character(len=LONGSTRLEN), allocatable :: stks(:)
        character(len=:),          allocatable :: ext, fbody, fname, dest
        integer :: i
        if( trim(params_glob%sigma_est).eq.'group' )then
            ! one star file with # of micrograph/stack groups -> # of micrograph groups files
            allocate(stks(self%nmics))
            do i = 1, self%nmics
                fname   = basename(self%orig_stks(i))
                ext     = fname2ext(fname)
                fbody   = get_fbody(fname, ext)
                stks(i) = trim(folder)//'/'//trim(fbody)//'.star'
            enddo
            fname = trim(self%path)//trim(sigma2_star_from_iter(self%it))
            call split_sigma2_into_groups(fname, stks)
            deallocate(stks)
        else
            ! one star file
            fname = trim(self%path)//trim(sigma2_star_from_iter(self%it))
            dest  = trim(folder)//'/chunk_'//int2str(self%id)//'.star'
            call simple_copy_file(fname,dest)
        endif
    end subroutine split_sigmas_into

    ! classes generation at original sampling
    subroutine gen_final_cavgs( self, clines )
        class(stream_chunk),         intent(in)    :: self
        type(cmdline),  allocatable, intent(inout) :: clines(:)
        type(cmdline),    allocatable :: tmp(:)
        type(cmdline)                 :: cline_make_cavgs
        character(len=:), allocatable :: finalcavgs
        integer :: n
        if( .not.allocated(clines) ) THROW_HARD('Fatal error gen_final_cavgs')
        finalcavgs = 'final_cavgs.mrc'
        call cline_make_cavgs%set('prg',        'make_cavgs')
        call cline_make_cavgs%set('mkdir',      'no')
        call cline_make_cavgs%set('refs',       finalcavgs)
        call cline_make_cavgs%set('ncls',       self%cline%get_iarg('ncls'))
        call cline_make_cavgs%set('mskdiam',    params_glob%mskdiam)
        call cline_make_cavgs%set('async',      'yes')
        call cline_make_cavgs%set('nthr',       params_glob%nthr2D)
        if( self%cline%defined('nparts') )then
            call cline_make_cavgs%set('nparts', self%cline%get_iarg('nparts'))
        else
            call cline_make_cavgs%set('nparts', 1)
        endif
        n = size(clines)
        allocate(tmp(1:n),source=clines(1:n))
        call clines(:)%kill; deallocate(clines); allocate(clines(n+1))
        clines(1:n) = tmp(1:n)
        clines(n+1) = cline_make_cavgs
        call tmp(:)%kill; call cline_make_cavgs%kill; deallocate(tmp)
    end subroutine gen_final_cavgs

    ! removes processing folder
    subroutine remove_folder( self )
        class(stream_chunk), intent(inout) :: self
        call debug_print('in chunk%remove_folder '//int2str(self%id))
        if( .not.self%converged )THROW_HARD('cannot remove chunk prior to convergence; remove_folder')
        call simple_rmdir(self%path)
        call debug_print('end chunk%remove_folder '//int2str(self%id))
    end subroutine remove_folder

    ! to interrupt processing
    subroutine terminate( self )
        class(stream_chunk), intent(inout) :: self
        character(len=XLONGSTRLEN)  :: cwd
        integer                     :: ipart, numlen
        if( self%id == 0 )   return
        if( file_exists(self%path) )then
            numlen = len(int2str(params_glob%nparts_chunk))
            call chdir(self%path)
            call simple_getcwd(cwd)
            cwd_glob = trim(cwd)
            call qsys_cleanup(keep2D=.false.)
            do ipart = 1,params_glob%nparts_chunk
                call simple_touch(trim(JOB_FINISHED_FBODY)//int2str_pad(ipart,numlen),errmsg="chunk%terminate")
            enddo
            call simple_touch('CAVGASSEMBLE_FINISHED',errmsg="chunk%terminate")
            call chdir('..')
            call simple_getcwd(cwd_glob)
        endif
    end subroutine terminate

    ! get & display convergence stats
    subroutine display_iter( self )
        class(stream_chunk), intent(inout) :: self
        type(oris)                 :: os
        character(len=XLONGSTRLEN) :: fname
        real                       :: mi_class,frac,corr
        integer                    :: it
        fname = trim(self%path)//trim(STATS_FILE)
        call debug_print('in chunk%display_iter '//int2str(self%id)//' '//trim(fname))
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

    ! Whether 2D analysis is complete
    logical function has_converged( self )
        class(stream_chunk), intent(inout) :: self
        if( .not.self%converged )then
            if( self%toanalyze2D )then
                if( str_has_substr(self%cline%get_carg('prg'), 'abinitio') )then
                    self%converged = file_exists(trim(self%path)//trim(ABINITIO2D_FINISHED))
                else
                    self%converged = file_exists(trim(self%path)//trim(CLUSTER2D_FINISHED))
                endif
            else
                self%converged = file_exists(trim(self%path)//trim(CALCPSPEC_FINISHED))
            endif
        endif
        has_converged  = self%converged
    end function has_converged

    ! Handles automated 2D analysis
    subroutine reject( self, res_thresh, ndev)
        class(stream_chunk), intent(inout) :: self
        real,                intent(in)    :: res_thresh, ndev
        type(image)                   :: img
        logical,          allocatable :: cls_mask(:), moments_mask(:), corres_mask(:)
        integer,          allocatable :: pops(:)
        character(len=:), allocatable :: cavgs
        character(len=XLONGSTRLEN)    :: projfile
        real                  :: smpd_here
        integer               :: nptcls_rejected, ncls_rejected, iptcl, box
        integer               :: icls, ncls, ncls_rejected_populated, ncls_populated
        call debug_print('in chunk%reject '//int2str(self%id))
        projfile        = trim(self%path)//self%projfile_out
        ncls_rejected   = 0
        nptcls_rejected = 0
        call self%spproj%read_segment('cls2D',projfile)
        call self%spproj%read_segment('out',  projfile)
        call self%spproj%get_cavgs_stk(cavgs, ncls, smpd_here)
        cavgs = trim(self%path)//basename(cavgs)
        box   = self%cline%get_iarg('box_crop')
        allocate(cls_mask(ncls),moments_mask(ncls),corres_mask(ncls),source=.true.)
        allocate(pops(ncls),source=nint(self%spproj%os_cls2D%get_all('pop')))
        ! moments & total variation distance
        if( trim(params_glob%reject_cls).ne.'old' )then
            call class_rejection(self%spproj%os_cls2D, moments_mask)
        endif
        ! correlation and resolution
        if( params_glob%lpthres < LOWRES_REJECT_THRESHOLD )then
            call self%spproj%os_cls2D%find_best_classes(box, smpd_here, res_thresh, corres_mask, ndev)
        endif
        ! overall class rejection
        cls_mask      = moments_mask .and. corres_mask
        ncls_rejected = count(.not.cls_mask)
        ncls_rejected_populated = count((.not.cls_mask).and.(pops>0))
        ncls_populated          = count(pops>0)
        if( ncls_rejected == 0 .or.&
            &ncls_rejected_populated >= min(ncls_populated,nint(real(ncls_populated)*FRAC_SKIP_REJECTION)) )then
            ! no or too many classes to reject
        else
            call self%spproj%read_segment('ptcl2D',projfile)
            ! rejects particles 2D
            !$omp parallel do private(iptcl,icls) reduction(+:nptcls_rejected) proc_bind(close)
            do iptcl=1,self%nptcls
                if( self%spproj%os_ptcl2D%get_state(iptcl) == 0 )cycle
                icls = self%spproj%os_ptcl2D%get_class(iptcl)
                if( cls_mask(icls) ) cycle
                nptcls_rejected = nptcls_rejected+1
                call self%spproj%os_ptcl2D%delete_2Dclustering(iptcl)
                call self%spproj%os_ptcl2D%set_state(iptcl,0)
            enddo
            !$omp end parallel do
            call self%spproj%write_segment_inside('ptcl2D',projfile)
            ! updates class averages
            call img%new([box,box,1],smpd_here)
            do icls = 1,ncls
                if( cls_mask(icls) ) cycle
                if( pops(icls) > 0 )then
                    ! update to global counter, does not include empty classes
                    ncls_rejected_glob = ncls_rejected_glob + 1
                    call img%read(cavgs,icls)
                    call img%write(CLS_CHUNK_REJECTED,ncls_rejected_glob)
                endif
            enddo
            call img%kill
            ! updates cls2D field
            do icls=1,ncls
                if( .not.cls_mask(icls) )then
                    call self%spproj%os_cls2D%set(icls,'pop',    0)
                    call self%spproj%os_cls2D%set_state(icls,    0)
                    call self%spproj%os_cls2D%set(icls,'corr', -1.)
                endif
            enddo
            call self%spproj%write_segment_inside('cls2D',projfile)
            write(logfhandle,'(A,I6,A,I6,A,I6,A,I6,A)')'>>> REJECTED FROM CHUNK ',self%id,': ',&
                &nptcls_rejected,' / ',self%nptcls,' PARTICLES IN ',ncls_rejected_populated,' CLUSTERS'
        endif
        call self%spproj%kill
        call debug_print('end chunk%reject '//int2str(self%id))
    end subroutine reject

    ! For debugging
    subroutine print_info( self )
        class(stream_chunk), intent(inout) :: self
        print *,'self%id           : ',self%id
        print *,'self%path         : ',trim(self%path)
        print *,'self%it           : ',self%it
        print *,'self%nmics        : ',self%nmics
        print *,'self%nptcls       : ',self%nptcls
        print *,'self%converged    : ',self%converged
        print *,'self%available    : ',self%available
    end subroutine print_info

    subroutine kill( self )
        class(stream_chunk), intent(inout) :: self
        self%id        = 0
        self%it        = 0
        call self%spproj%kill
        call self%qenv%kill
        call self%cline%kill
        self%nmics     = 0
        self%nptcls    = 0
        self%path      = ''
        self%projfile_out = ''
        if( allocated(self%orig_stks) ) deallocate(self%orig_stks)
        self%toanalyze2D = .true.
        self%converged  = .false.
        self%available  = .false.
    end subroutine kill

    ! Public Utility to handle the type projrecord

    !> convert a list of projects into one project
    !  previous mic/stk/ptcl2D,ptcl3D are wipped, other fields untouched
    subroutine projrecords2proj( records, spproj )
        class(projrecord), intent(in)    :: records(:)
        class(sp_project), intent(inout) :: spproj
        type(sp_project)              :: tmpproj
        character(len=:), allocatable :: stack_name, projname, prev_projname
        integer :: iptcl, fromp, ifromp, itop, jptcl, nptcls_tot
        integer :: nrecs, nmics, nptcls, imic, micind
        logical :: has_ptcl
        call spproj%os_mic%kill
        call spproj%os_stk%kill
        call spproj%os_ptcl2D%kill
        call spproj%os_ptcl3D%kill
        nrecs      = size(records)
        if( nrecs == 0 ) return 
        nmics      = nrecs
        nptcls_tot = sum(records(:)%nptcls)
        has_ptcl   = nptcls_tot > 0
        call spproj%os_mic%new(nmics,is_ptcl=.false.)
        call spproj%os_stk%new(nmics,is_ptcl=.false.)
        if( has_ptcl ) call spproj%os_ptcl2D%new(nptcls_tot,is_ptcl=.true.)
        prev_projname = ''
        jptcl = 0
        fromp = 1
        do imic = 1,nmics
            ! read individual project (up to NMOVS_SET entries)
            projname = trim(records(imic)%projname)
            if( trim(projname) /= trim(prev_projname) )then
                call tmpproj%kill
                call tmpproj%read_mic_stk_ptcl2D_segments(projname)
                prev_projname = trim(projname)
            endif
            ! mic
            micind = records(imic)%micind
            call spproj%os_mic%transfer_ori(imic, tmpproj%os_mic, micind)
            ! stack
            nptcls = records(imic)%nptcls
            if( nptcls == 0 )cycle
            call spproj%os_stk%transfer_ori(imic, tmpproj%os_stk, micind)
            ! update stack path to absolute
            stack_name = trim(spproj%get_stkname(imic))
            if( stack_name(1:1) == '/' )then
                ! already absolute path, should always be the case
            else if( stack_name(1:3) == '../' )then
                stack_name = simple_abspath(trim(stack_name))
                call spproj%os_stk%set(imic, 'stk', stack_name)
            else
                THROW_HARD('Unexpected file path format for: '//trim(stack_name))
            endif
            ! particles
            ifromp = spproj%os_stk%get_fromp(imic)
            itop   = spproj%os_stk%get_top(imic)
            do iptcl = ifromp,itop
                jptcl = jptcl+1 ! global index
                call spproj%os_ptcl2D%transfer_ori(jptcl, tmpproj%os_ptcl2D, iptcl)
                call spproj%os_ptcl2D%set_stkind(jptcl, imic)
            enddo
            call spproj%os_stk%set(imic, 'fromp', fromp)
            call spproj%os_stk%set(imic, 'top',   fromp+nptcls-1)
            fromp = fromp + nptcls
        enddo
        call tmpproj%kill
        if( has_ptcl ) spproj%os_ptcl3D = spproj%os_ptcl2D
    end subroutine projrecords2proj

    ! Public utility to handle the type procrecord

    subroutine append_procrecord( records, id, folder, projfile )
        type(procrecord), allocatable, intent(inout) :: records(:)
        character(len=*),              intent(in)    :: id, folder, projfile
        type(procrecord) :: record
        record%id        = trim(id)
        record%folder    = trim(folder)
        record%projfile  = trim(projfile)
        record%volume    = ''
        record%alnvolume = ''
        record%submitted = .false.
        record%complete  = .false.
        record%included  = .false.
        if( allocated(records) )then
            records = [records(:), record]
        else
            allocate(records(1), source=[record])
        endif
        call kill_procrecord(record)
    end subroutine append_procrecord

    elemental subroutine kill_procrecord( rec )
        type(procrecord), intent(inout) :: rec
        if( allocated(rec%id) ) deallocate(rec%id,rec%folder,rec%projfile,rec%volume,rec%alnvolume)
    end subroutine kill_procrecord

    subroutine kill_procrecords( recs )
        type(procrecord), allocatable, intent(inout) :: recs(:)
        if( allocated(recs) )then
            call kill_procrecord(recs(:))
            deallocate(recs)
        endif
    end subroutine kill_procrecords

    ! Public utility

    ! Class rejection routine based on image moments & Total Variation Distance
    subroutine class_rejection( os, mask, adjust )
        class(oris),    intent(in)    :: os
        logical,        intent(inout) :: mask(:)
        real, optional, intent(in)    :: adjust
        real,    allocatable :: vals(:), x(:)
        logical, allocatable :: msk(:)
        real    :: eff_mean_thresh, eff_rel_var_thresh, eff_abs_var_thresh
        real    :: eff_tvd_thresh, eff_min_thresh, eff_max_thresh
        integer :: icls, i, n
        logical :: has_mean, has_var, has_tvd, has_minmax
        n   = os%get_noris()
        if( size(mask) /= n )THROW_HARD('Incompatible sizes! class rejection')
        msk = mask
        if( os%isthere('pop') )then
            do icls=1,n
                msk(icls) = os%get(icls,'pop') > 0.5
            enddo
        endif
        mask = msk
        if( count(msk) <= 5 )then
            deallocate(msk)
            return
        endif
        ! Effective threshold
        eff_mean_thresh    = params_glob%stream_mean_threshold
        eff_rel_var_thresh = params_glob%stream_rel_var_threshold
        eff_abs_var_thresh = params_glob%stream_abs_var_threshold
        eff_tvd_thresh     = max(0.001,min(0.999,params_glob%stream_tvd_theshold))
        eff_min_thresh     = -params_glob%stream_minmax_threshold
        eff_max_thresh     =  params_glob%stream_minmax_threshold
        if( present(adjust) )then
            eff_mean_thresh    = adjust * eff_mean_thresh
            eff_rel_var_thresh = adjust * eff_rel_var_thresh
            eff_abs_var_thresh = adjust * eff_abs_var_thresh
            eff_tvd_thresh     = min(0.999, adjust * eff_tvd_thresh)
            eff_min_thresh     = adjust * eff_min_thresh
            eff_max_thresh     = adjust * eff_max_thresh
        endif
        ! selection
        has_mean   = os%isthere('mean')
        has_var    = os%isthere('var')
        has_tvd    = os%isthere('tvd')
        has_minmax = os%isthere('min') .and. os%isthere('max')
        ! Mean
        if( has_mean )then
            vals = os%get_all('mean')
            x    = pack(vals, mask=msk)
            call robust_scaling(x)
            i = 0
            do icls = 1,n
                if( msk(icls) )then
                    i = i+1
                    if( mask(icls) ) mask(icls) = x(i) > eff_mean_thresh
                endif
            enddo
        endif
        ! Variance
        if( has_var )then
            vals = os%get_all('var')
            x    = pack(vals, mask=msk)
            call robust_scaling(x)
            i = 0
            do icls = 1,n
                if( msk(icls) )then
                    i = i+1
                    if( mask(icls) ) mask(icls) = x(i)       < eff_rel_var_thresh
                    if( mask(icls) ) mask(icls) = vals(icls) < eff_abs_var_thresh
                endif
            enddo
        endif
        ! Total Variation Distance
        if( has_tvd )then
            vals = os%get_all('tvd')
            do icls = 1,n
                if( mask(icls) ) mask(icls) = vals(icls) < eff_tvd_thresh
            enddo
        endif
        ! Min/max
        if( has_minmax )then
            do icls = 1,n
                if( mask(icls) )then
                    if(  (os%get(icls,'min') < eff_min_thresh).and.&
                        &(os%get(icls,'max') > eff_max_thresh) )then
                        mask(icls) = .false.
                    endif
                endif
            enddo
        endif
        deallocate(msk)
        if(allocated(vals) ) deallocate(vals, x)
    end subroutine class_rejection

end module simple_stream_utils
