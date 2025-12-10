module simple_stream_chunk
include 'simple_lib.f08'
use simple_cmdline,    only: cmdline
use simple_image,      only: image
use simple_parameters, only: params_glob
use simple_qsys_env,   only: qsys_env
use simple_sp_project, only: sp_project
use simple_stream_utils
use simple_rec_list
use simple_qsys_funs
implicit none

public :: stream_chunk
private
#include "simple_local_flags.inc"

! Type to handle a single chunk
type stream_chunk
    private
    type(sp_project)          :: spproj                ! master project
    type(qsys_env)            :: qenv                  ! submission handler
    type(cmdline)             :: cline                 ! command line
    type(string), allocatable :: orig_stks(:)          ! list of stacks
    type(string)              :: path, projfile_out    ! physical location
    integer                   :: id                    ! unique id
    integer                   :: it                    ! # of iterations performed
    integer                   :: nmics                 ! # of micrographs
    integer                   :: nptcls                ! # number of particles
    logical                   :: toanalyze2D = .true.  ! whether to perform 2D analysis or only calculate sigmas2
    logical                   :: converged   = .false. ! whether 2D analysis is over
    logical                   :: available   = .true.  ! has been initialized but no 2D analysis peformed
  contains
    procedure          :: init_chunk
    procedure, private :: copy
    procedure, private :: assign
    generic            :: assignment(=) => assign
    procedure          :: generate
    procedure          :: get_id
    procedure          :: is_available
    procedure          :: to_analyze2D
    procedure          :: get_projfile_fname
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
    procedure          :: terminate_chunk
    procedure          :: kill
end type stream_chunk

integer            :: ncls_rejected_glob = 0 ! counter of rejected classes
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
    subroutine init_chunk( self, id, cline, master_spproj )
        class(stream_chunk), intent(inout) :: self
        class(cmdline),      intent(in)    :: cline
        integer,             intent(in)    :: id
        class(sp_project),   intent(in)    :: master_spproj
        type(string)          :: exec, str_prg
        character(len=STDLEN) :: chunk_part_env
        integer               :: envlen
        call debug_print('in chunk%init '//int2str(id))
        call self%spproj%kill
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
        str_prg              = self%cline%get_carg('prg')
        if( str_prg%has_substr('abinitio') )then
            exec = 'simple_exec'
        else
            exec = 'simple_private_exec'
        endif
        if( params_glob%nparts_chunk == 1 )then
            ! shared memory
            call self%qenv%new(params_glob%nparts_chunk, exec_bin=exec, qsys_nthr=params_glob%nthr2D)
            call self%spproj%compenv%set(1,'qsys_name','local')
        else
            ! we need to override the qsys_name for non local distributed execution
            call self%qenv%new(params_glob%nparts_chunk, exec_bin=exec, qsys_name=string('local'))
            call get_environment_variable(SIMPLE_STREAM_CHUNK_PARTITION, chunk_part_env, envlen)
            if(envlen > 0) call self%spproj%compenv%set(1,'qsys_partition', trim(chunk_part_env))
        endif
        call self%spproj%projinfo%delete_entry('projname')
        call self%spproj%projinfo%delete_entry('projfile')
        if( allocated(self%orig_stks) )then
            call self%orig_stks%kill
            deallocate(self%orig_stks)
        endif
        self%toanalyze2D = .true.
        self%converged  = .false.
        self%available  = .true.
        call str_prg%kill
        call debug_print('end chunk%init '//int2str(id))
    end subroutine init_chunk

    subroutine copy( dest, src )
        class(stream_chunk), intent(inout) :: dest
        class(stream_chunk), intent(in)    :: src
        call dest%kill
        dest%spproj       = src%spproj
        dest%qenv         = src%qenv
        dest%cline        = src%cline
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
    end subroutine copy

    !>  \brief  assign, polymorphic assignment (=)
    subroutine assign( selfout, selfin )
        class(stream_chunk), intent(inout) :: selfout
        class(stream_chunk), intent(in)    :: selfin
        call selfout%copy(selfin)
    end subroutine assign

    subroutine generate( self, project_list )
        class(stream_chunk), intent(inout) :: self
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

    integer function get_id( self )
        class(stream_chunk), intent(in) :: self
        get_id = self%id
    end function get_id

    elemental function is_available( self ) result( avail )
        class(stream_chunk), intent(in) :: self
        logical :: avail
        avail = self%available
    end function is_available

    elemental function to_analyze2D( self ) result( yes )
        class(stream_chunk), intent(in) :: self
        logical :: yes
        yes = self%toanalyze2D
    end function to_analyze2D

    function get_projfile_fname( self )result( fname )
        class(stream_chunk), intent(in) :: self
        type(string) :: fname
        fname = self%path//self%projfile_out
    end function get_projfile_fname

    ! Initiates 2D analysis
    subroutine analyze2D( self, calc_pspec, makecavgs )
        class(stream_chunk), intent(inout) :: self
        logical,             intent(in)    :: calc_pspec
        logical,   optional, intent(in)    :: makecavgs
        type(string),  allocatable :: bins(:)
        type(cmdline), allocatable :: clines(:)
        type(cmdline) :: cline_pspec
        type(string)  :: cwd, str_prg
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
        str_prg = self%cline%get_carg('prg')
        if( str_prg%has_substr('abinitio') )then
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
        call str_prg%kill
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
        call cline_pspec%set('prg',      'calc_pspec_distr')
        call cline_pspec%set('oritype',  'ptcl2D')
        call cline_pspec%set('nthr',     cline_analyze2D%get_iarg('nthr'))
        call cline_pspec%set('mkdir',    'yes')
        call cline_pspec%set('nparts',   1)
        if( params_glob%nparts_chunk > 1 ) call cline_pspec%set('nparts',params_glob%nparts_chunk)
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

    ! read project and deals with arrays used for fractional update
    subroutine read( self, box )
        class(stream_chunk), intent(inout) :: self
        integer,             intent(in)    :: box
        type(image)  :: img, avg
        type(string) :: projfile
        call debug_print('in chunk%read '//int2str(self%id))
        if( .not.self%converged )THROW_HARD('cannot read chunk prior to convergence')
        ! doc & parameters
        projfile = self%path//self%projfile_out
        call self%spproj%read(projfile)
        ! classes, to account for nparts /= nparts_chunk
        if( self%toanalyze2D )then
            call img%new([box,box,1],1.0)
            call avg%new([box,box,1],1.0)
            call average_into(self%path//'cavgs_even_part')
            call average_into(self%path//'cavgs_odd_part')
            call average_into(self%path//'ctfsqsums_even_part')
            call average_into(self%path//'ctfsqsums_odd_part')
            call img%kill
            call avg%kill
        endif
        call debug_print('end chunk%read '//int2str(self%id))
        contains

            subroutine average_into(tmpl)
                class(string), intent(in) :: tmpl
                type(string) :: fname
                integer      :: icls, ipart, numlen_chunk
                if( params_glob%nparts_chunk > 1  )then
                    numlen_chunk = len(int2str(params_glob%nparts_chunk)) ! as per parameters
                    call img%zero_and_flag_ft
                    do icls = 1,params_glob%ncls_start
                        call avg%zero_and_flag_ft
                        do ipart = 1,params_glob%nparts_chunk
                            fname = tmpl//int2str_pad(ipart,numlen_chunk)//params_glob%ext%to_char()
                            call img%read(fname,icls)
                            call avg%add(img)
                        enddo
                        call avg%div(real(params_glob%nparts_chunk))
                        call avg%write(tmpl//params_glob%ext%to_char(),icls)
                    enddo
                else
                    fname = tmpl//'1'//params_glob%ext%to_char()
                    call simple_rename(fname,tmpl//params_glob%ext)
                endif
            end subroutine average_into

    end subroutine read

    ! split sigmas into individually named per stack documents
    subroutine split_sigmas_into( self, folder )
        use simple_euclid_sigma2, only: split_sigma2_into_groups, sigma2_star_from_iter
        class(stream_chunk), intent(in) :: self
        class(string),       intent(in) :: folder
        type(string), allocatable :: stks(:)
        type(string) :: ext, fbody, fname, dest
        integer      :: i
        if( trim(params_glob%sigma_est).eq.'group' )then
            ! one star file with # of micrograph/stack groups -> # of micrograph groups files
            allocate(stks(self%nmics))
            do i = 1, self%nmics
                fname   = basename(self%orig_stks(i))
                ext     = fname2ext(fname)
                fbody   = get_fbody(fname, ext)
                stks(i) = folder//'/'//fbody//STAR_EXT
            enddo
            fname = self%path//sigma2_star_from_iter(self%it)
            call split_sigma2_into_groups(fname, stks)
            deallocate(stks)
        else
            ! one star file
            fname = self%path//sigma2_star_from_iter(self%it)
            dest  = folder//'/chunk_'//int2str(self%id)//STAR_EXT
            call simple_copy_file(fname,dest)
        endif
    end subroutine split_sigmas_into

    ! classes generation at original sampling
    subroutine gen_final_cavgs( self, clines )
        class(stream_chunk),         intent(in)    :: self
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
        do i = 1, n
            clines(i) = tmp(i)
        end do
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
    subroutine terminate_chunk( self )
        class(stream_chunk), intent(inout) :: self
        type(string) :: cwd
        integer      :: ipart, numlen
        if( self%id == 0 )   return
        if( file_exists(self%path) )then
            numlen = len(int2str(params_glob%nparts_chunk))
            call simple_chdir(self%path)
            call simple_getcwd(cwd)
            CWD_GLOB = cwd%to_char()
            call qsys_cleanup(keep2D=.false.)
            do ipart = 1,params_glob%nparts_chunk
                call simple_touch(JOB_FINISHED_FBODY//int2str_pad(ipart,numlen))
            enddo
            call simple_touch('CAVGASSEMBLE_FINISHED')
            call simple_chdir('..')
            call simple_getcwd(cwd)
            CWD_GLOB = cwd%to_char()
        endif
    end subroutine terminate_chunk

    ! get & display convergence stats
    subroutine display_iter( self )
        class(stream_chunk), intent(inout) :: self
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

    ! Whether 2D analysis is complete
    logical function has_converged( self )
        class(stream_chunk), intent(inout) :: self
        type(string) :: str_prg
        if( .not.self%converged )then
            if( self%toanalyze2D )then
                str_prg = self%cline%get_carg('prg')
                if( str_prg%has_substr('abinitio') )then
                    self%converged = file_exists(self%path//ABINITIO2D_FINISHED)
                else
                    self%converged = file_exists(self%path//CLUSTER2D_FINISHED)
                endif
                call str_prg%kill
            else
                self%converged = file_exists(self%path//CALCPSPEC_FINISHED)
            endif
        endif
        has_converged  = self%converged
    end function has_converged

    ! Handles automated 2D analysis
    subroutine reject( self, res_thresh, ndev)
        class(stream_chunk), intent(inout) :: self
        real,                intent(in)    :: res_thresh, ndev
        logical, allocatable :: cls_mask(:), moments_mask(:), corres_mask(:)
        integer, allocatable :: pops(:)
        type(image)  :: img
        type(string) :: cavgs, projfile
        real         :: smpd_here
        integer      :: nptcls_rejected, ncls_rejected, iptcl, box
        integer      :: icls, ncls, ncls_rejected_populated, ncls_populated
        call debug_print('in chunk%reject '//int2str(self%id))
        projfile        = self%path//self%projfile_out
        ncls_rejected   = 0
        nptcls_rejected = 0
        call self%spproj%read_segment('cls2D',projfile)
        call self%spproj%read_segment('out',  projfile)
        call self%spproj%get_cavgs_stk(cavgs, ncls, smpd_here)
        cavgs = self%path//basename(cavgs)
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
                    call img%write(string(CHUNK_CLS_REJECTED),ncls_rejected_glob)
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
        print *,'self%path         : ',self%path%to_char()
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
        if( allocated(self%orig_stks) )then
            call self%orig_stks%kill
            deallocate(self%orig_stks)
        endif
        self%toanalyze2D = .true.
        self%converged  = .false.
        self%available  = .false.
    end subroutine kill

end module simple_stream_chunk