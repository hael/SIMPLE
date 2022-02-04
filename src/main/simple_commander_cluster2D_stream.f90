! concrete commander: cluster2D_stream for streaming 2D alignment and clustering of single-particle images
module simple_commander_cluster2D_stream
include 'simple_lib.f08'
use simple_builder,        only: builder
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_parameters,     only: parameters, params_glob
use simple_sp_project,     only: sp_project
use simple_qsys_env,       only: qsys_env
use simple_image,          only: image
use simple_oris,           only: oris
use simple_qsys_funs
use simple_commander_cluster2D
use simple_timer
implicit none

public :: cluster2D_commander_stream

private
#include "simple_local_flags.inc"

real,                  parameter   :: GREEDY_TARGET_LP    = 15.
integer,               parameter   :: MINBOXSZ            = 72    ! minimum boxsize for scaling
integer,               parameter   :: WAIT_WATCHER        = 5     ! seconds prior to new stack detection
! integer,               parameter   :: ORIGPROJ_WRITEFREQ  = 1800  ! dev settings
integer,               parameter   :: ORIGPROJ_WRITEFREQ  = 7200  ! Frequency at which the original project file should be updated
integer,               parameter   :: FREQ_POOL_REJECTION = 5     !
character(len=STDLEN), parameter   :: USER_PARAMS         = 'stream2D_user_params.txt'
character(len=STDLEN), parameter   :: LAST_SNAPSHOT       = 'last_snapshot.txt'
character(len=STDLEN), parameter   :: SPPROJ_SNAPSHOT     = 'SIMPLE_PROJECT_SNAPSHOT'
character(len=STDLEN), parameter   :: PROJFILE_POOL       = 'cluster2D.simple'
character(len=STDLEN), parameter   :: PROJNAME_CHUNK      = 'chunk'
character(len=STDLEN), parameter   :: SCALE_DIR           = './scaled_stks/'
character(len=STDLEN), parameter   :: DISTR_EXEC_FNAME    = './distr_cluster2D_pool'
logical,               parameter   :: DEBUG_HERE          = .false.

type, extends(commander_base) :: cluster2D_commander_stream
  contains
    procedure :: execute      => exec_cluster2D_stream
end type cluster2D_commander_stream

type chunk
    type(sp_project)                       :: spproj
    type(qsys_env)                         :: qenv
    character(len=LONGSTRLEN), allocatable :: orig_stks(:)
    character(len=LONGSTRLEN)              :: path, projfile_out
    integer                                :: id
    integer                                :: it
    integer                                :: nmics
    integer                                :: nptcls
    logical                                :: converged = .false.
    logical                                :: available = .true.
contains
    procedure :: init
    procedure :: generate
    procedure :: exec_classify
    procedure :: read
    procedure :: remove_folder
    procedure :: display_iter
    procedure :: reject
    procedure :: has_converged
    procedure :: print_info
    procedure :: terminate
    procedure :: kill
end type chunk

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

    subroutine init( self, id, master_spproj )
        class(chunk),      intent(inout) :: self
        integer,           intent(in)    :: id
        class(sp_project), intent(in)    :: master_spproj
        call debug_print('in chunk%init '//int2str(id))
        call self%spproj%kill
        self%id        = id
        self%it        = 0
        self%nmics     = 0
        self%nptcls    = 0
        self%path      = './chunk_'//int2str(id)//'/'
        self%projfile_out = ''
        ! we need to override the qsys_name for non local distributed execution
        call self%qenv%new(params_glob%nparts_chunk, exec_bin='simple_private_exec', qsys_name='local')
        self%spproj%projinfo = master_spproj%projinfo
        self%spproj%compenv  = master_spproj%compenv
        call self%spproj%projinfo%delete_entry('projname')
        call self%spproj%projinfo%delete_entry('projfile')
        if( allocated(self%orig_stks) ) deallocate(self%orig_stks)
        self%converged = .false.
        self%available = .true.
        call debug_print('end chunk%init '//int2str(id))
    end subroutine init

    subroutine generate( self, fnames, nptcls, ind_glob )
        class(chunk),              intent(inout) :: self
        character(len=LONGSTRLEN), intent(in)    :: fnames(:)
        integer,                   intent(in)    :: nptcls, ind_glob
        type(sp_project)     :: spproj
        logical, allocatable :: spproj_mask(:)
        integer :: iproj, iptcl, cnt, inptcls, nptcls_tot, fromp, iiproj, n_in
        if( .not.self%available ) THROW_HARD('chunk unavailable; chunk%generate')
        n_in = size(fnames)
        if( n_in == 0 ) THROW_HARD('# ptcls == 0; chunk%generate')
        call debug_print('in chunk%generate '//int2str(self%id)//' '//int2str(ind_glob)//' '//int2str(nptcls)//' '//int2str(n_in))
        allocate(spproj_mask(n_in),source=.false.)
        nptcls_tot = 0
        do iproj = 1,n_in
            call spproj%read_segment('mic', fnames(iproj))
            inptcls = nint(spproj%os_mic%get(1,'nptcls'))
            spproj_mask(iproj) = inptcls > 0
            if( spproj_mask(iproj) ) nptcls_tot = nptcls_tot + inptcls
        enddo
        call spproj%kill
        if( nptcls /= nptcls_tot ) THROW_HARD('Inconsistent # of particles; chunk%generate')
        self%nmics  = count(spproj_mask)
        self%nptcls = nptcls
        call self%spproj%os_mic%new(self%nmics, is_ptcl=.false.)
        call self%spproj%os_stk%new(self%nmics, is_ptcl=.false.)
        call self%spproj%os_ptcl2D%new(self%nptcls, is_ptcl=.true.)
        allocate(self%orig_stks(self%nmics))
        call debug_print('in chunk%generate '//int2str(self%nmics)//' '//int2str(self%nptcls))
        cnt    = 0
        fromp  = 1
        iiproj = 0
        do iproj = 1,n_in
            if( .not.spproj_mask(iproj) ) cycle
            iiproj = iiproj + 1
            call spproj%read_segment('mic', fnames(iproj))
            inptcls = nint(spproj%os_mic%get(1,'nptcls'))
            call self%spproj%os_mic%transfer_ori(iiproj, spproj%os_mic, 1)
            call spproj%read_segment('stk', fnames(iproj))
            call self%spproj%os_stk%transfer_ori(iiproj, spproj%os_stk, 1)
            self%orig_stks(iiproj) = simple_abspath(trim(spproj%get_stkname(1)))
            call self%spproj%os_stk%set(iiproj, 'stk', self%orig_stks(iiproj))
            call spproj%read_segment('ptcl2D', fnames(iproj))
            do iptcl = 1,inptcls
                cnt = cnt + 1
                call self%spproj%os_ptcl2D%transfer_ori(cnt, spproj%os_ptcl2D, iptcl)
                call self%spproj%os_ptcl2D%set(cnt, 'stkind', real(iiproj))
            enddo
            call self%spproj%os_stk%set(iiproj, 'fromp', real(fromp))
            call self%spproj%os_stk%set(iiproj, 'top',   real(fromp+inptcls-1))
            fromp = fromp + inptcls
            call debug_print('in chunk%generate '//int2str(iproj)//' '//int2str(iiproj)//' '//int2str(inptcls)//' '//int2str(fromp))
        enddo
        call spproj%kill
        self%spproj%os_ptcl3D = self%spproj%os_ptcl2D
        call debug_print('end chunk%generate '//int2str(self%id))
    end subroutine generate

    subroutine exec_classify( self, cline_classify, orig_smpd, orig_box, box )
        class(chunk),   intent(inout) :: self
        class(cmdline), intent(inout) :: cline_classify
        real,           intent(in)    :: orig_smpd
        integer,        intent(in)    :: orig_box, box
        type(cmdline)              :: cline_scale
        character(len=XLONGSTRLEN) :: cwd
        character(len=LONGSTRLEN)  :: projfile_in
        logical :: err
        call debug_print('in chunk%exec_classify '//int2str(self%id))
        if( .not.self%available ) return
        if( self%nptcls == 0 ) return
        err = .false.
        call simple_mkdir(self%path)
        call chdir(self%path)
        call simple_getcwd(cwd)
        cwd_glob    = trim(cwd)
        projfile_in = trim(PROJNAME_CHUNK)//trim(METADATA_EXT)
        ! scaling
        if( box /= orig_box )then
            self%projfile_out = trim(PROJNAME_CHUNK)//SCALE_SUFFIX//trim(METADATA_EXT) ! as per scale_project convention
            call cline_scale%set('prg',        'scale_project_distr')
            call cline_scale%set('projname',   trim(PROJNAME_CHUNK))
            call cline_scale%set('projfile',   projfile_in)
            call cline_scale%set('smpd',       orig_smpd)
            call cline_scale%set('box',        real(orig_box))
            call cline_scale%set('newbox',     real(box))
            call cline_scale%set('nparts',     real(params_glob%nparts_chunk))
            call cline_scale%set('nthr',       real(params_glob%nthr))
            call cline_scale%set('mkdir',      'yes') ! required but not done as it is a simple_private_exec
            call cline_scale%set('dir_target', '../'//trim(SCALE_DIR))
            call self%spproj%update_projinfo(cline_scale)
            call self%spproj%write(projfile_in)
            call cline_classify%delete('projname')
            call cline_classify%set('projfile',self%projfile_out)
            call self%qenv%exec_simple_prg_in_queue_async(cline_scale, './distr_chunk2D', 'simple_log_chunk2d', cline2=cline_classify )
        else
            self%projfile_out = trim(projfile_in)
            call cline_classify%set('projfile', projfile_in)
            call cline_classify%set('projname', trim(PROJNAME_CHUNK))
            call self%spproj%update_projinfo(cline_classify)
            call self%spproj%write()
            call self%qenv%exec_simple_prg_in_queue_async(cline_classify, './distr_chunk2D', 'simple_log_chunk2d')
        endif
        call chdir('..')
        call simple_getcwd(cwd_glob)
        ! cleanup
        call self%spproj%kill
        call cline_scale%kill
        call self%qenv%kill
        ! chunk is now busy
        self%available = .false.
        self%converged = .false.
        write(logfhandle,'(A,I6,A,I6,A)')'>>> CHUNK ',self%id,&
            &' INITIATED CLASSIFICATION WITH ',self%nptcls,' PARTICLES'
        call debug_print('end chunk%exec_classify')
    end subroutine exec_classify

    subroutine read( self, box )
        class(chunk), intent(inout) :: self
        integer,      intent(in)    :: box
        type(image)                :: img, avg
        character(len=XLONGSTRLEN) :: projfile
        call debug_print('in chunk%read '//int2str(self%id))
        if( .not.self%converged )THROW_HARD('cannot read chunk prior to convergence')
        ! doc & parameters
        projfile = trim(self%path)//trim(self%projfile_out)
        call self%spproj%read_segment('mic',   projfile)
        call self%spproj%read_segment('stk',   projfile)
        call self%spproj%read_segment('ptcl2D',projfile)
        call self%spproj%read_segment('cls2D', projfile)
        call self%spproj%read_segment('out',   projfile)
        ! classes, to account for nparts /= nparts_chunk
        call img%new([box,box,1],1.0)
        call avg%new([box,box,1],1.0)
        call average_into(trim(self%path)//'cavgs_even_part')
        call average_into(trim(self%path)//'cavgs_odd_part')
        call average_into(trim(self%path)//'ctfsqsums_even_part')
        call average_into(trim(self%path)//'ctfsqsums_odd_part')
        call img%kill
        call avg%kill
        call debug_print('end chunk%read '//int2str(self%id))

        contains

            subroutine average_into(tmpl)
                character(len=*), intent(in) :: tmpl
                character(len=XLONGSTRLEN) :: fname
                integer                    :: icls, ipart, numlen_chunk
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
            end subroutine average_into

    end subroutine read

    subroutine remove_folder( self )
        class(chunk), intent(inout) :: self
        call debug_print('in chunk%remove_folder '//int2str(self%id))
        if( .not.self%converged )THROW_HARD('cannot remove chunk prior to convergence; remove_folder')
        call simple_rmdir(self%path)
        call debug_print('end chunk%remove_folder '//int2str(self%id))
    end subroutine remove_folder

    subroutine terminate( self )
        class(chunk), intent(inout) :: self
        character(len=XLONGSTRLEN)  :: cwd
        integer                     :: ipart, numlen
        if( self%id == 0 )   return
        if( file_exists(self%path) )then
            numlen = len(int2str(params_glob%nparts_chunk))
            do ipart = 1,params_glob%nparts_chunk
                call chdir(self%path)
                call simple_getcwd(cwd)
                cwd_glob = trim(cwd)
                call qsys_cleanup(keep2D=.false.)
                call simple_touch('JOB_FINISHED_'//int2str_pad(ipart,numlen),errmsg="chunk%terminate")
                call simple_touch('CAVGASSEMBLE_FINISHED',errmsg="chunk%terminate")
                call chdir('..')
                call simple_getcwd(cwd_glob)
            enddo
        endif
    end subroutine terminate

    subroutine display_iter( self )
        class(chunk), intent(inout) :: self
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
                corr     = os%get(1,'CORRELATION')
                write(logfhandle,'(A,I6,A,I6,A,F7.3,A,F7.3,A,F7.3)')'>>> CHUNK ',self%id,' ITERATION ',self%it,&
                    &'; CLASS OVERLAP: ',mi_class,'; SEARCH SPACE SCANNED: ',frac,'; CORRELATION: ',corr
            endif
            call os%kill
        endif
        call debug_print('end chunk%display_iter '//int2str(self%id))
    end subroutine display_iter

    logical function has_converged( self )
        class(chunk), intent(inout) :: self
        call debug_print('in chunk%has_converged '//int2str(self%id))
        self%converged = file_exists(trim(self%path)//trim(CLUSTER2D_FINISHED))
        has_converged  = self%converged
        if( has_converged )then
            call debug_print('end chunk%has_converged T '//int2str(self%id))
        else
            call debug_print('end chunk%has_converged F '//int2str(self%id))
        endif
    end function has_converged

    subroutine reject( self, res_thresh, ndev, box )
        class(chunk), intent(inout) :: self
        real,          intent(in)    :: res_thresh, ndev
        integer,       intent(in)    :: box
        type(image) :: img
        logical,          allocatable :: cls_mask(:)
        character(len=:), allocatable :: cavgs
        character(len=XLONGSTRLEN) :: projfile
        real                  :: smpd_here
        integer               :: nptcls_rejected, ncls_rejected, iptcl
        integer               :: icls, ncls_here
        call debug_print('in chunk%reject '//int2str(self%id))
        projfile = trim(self%path)//self%projfile_out
        ncls_rejected   = 0
        nptcls_rejected = 0
        call self%spproj%read_segment('cls2D',projfile)
        call self%spproj%read_segment('out',  projfile)
        call self%spproj%get_cavgs_stk(cavgs, ncls_here, smpd_here)
        cavgs = trim(self%path)//basename(cavgs)
        allocate(cls_mask(ncls_here), source=.true.)
        call self%spproj%os_cls2D%find_best_classes(box, smpd_here, res_thresh, cls_mask, ndev)
        ncls_rejected = count(.not.cls_mask)
        if( ncls_rejected == 0 .or. ncls_rejected == ncls_here )then
            ! nothing to do
        else
            call self%spproj%read_segment('ptcl2D',projfile)
            ! rejects particles 2D
            do iptcl=1,self%nptcls
                if( self%spproj%os_ptcl2D%get_state(iptcl) == 0 )cycle
                icls = nint(self%spproj%os_ptcl2D%get(iptcl,'class'))
                if( cls_mask(icls) ) cycle
                nptcls_rejected = nptcls_rejected+1
                call self%spproj%os_ptcl2D%set(iptcl,'state',0.)
            enddo
            call debug_print('in chunk%reject '//int2str(self%id)//' '//int2str(nptcls_rejected))
            call self%spproj%write_segment_inside('ptcl2D',projfile)
            ! updates cls2D field
            do icls=1,ncls_here
                if( .not.cls_mask(icls) )then
                    call self%spproj%os_cls2D%set(icls,'pop',   0.)
                    call self%spproj%os_cls2D%set(icls,'state', 0.)
                    call self%spproj%os_cls2D%set(icls,'corr', -1.)
                endif
            enddo
            call self%spproj%write_segment_inside('cls2D',projfile)
            ! updates class averages
            call img%new([box,box,1],smpd_here)
            do icls=1,ncls_here
                if( cls_mask(icls) ) cycle
                img = 0.
                call img%write(cavgs,icls)
            enddo
            call img%read(cavgs, ncls_here)
            call img%write(cavgs,ncls_here)
            write(logfhandle,'(A,I6,A,I6,A,I6,A,I6,A)')'>>> REJECTED FROM CHUNK ',self%id,': ',&
                &nptcls_rejected,' / ',self%nptcls,' PARTICLES IN ',ncls_rejected,' CLUSTERS'
        endif
        call img%kill
        call self%spproj%kill
        call debug_print('end chunk%reject '//int2str(self%id))
    end subroutine reject

    subroutine print_info( self )
        class(chunk), intent(inout) :: self
        print *,'self%id           : ',self%id
        print *,'self%path         : ',trim(self%path)
        print *,'self%it           : ',self%it
        print *,'self%nmics        : ',self%nmics
        print *,'self%nptcls       : ',self%nptcls
        print *,'self%converged    : ',self%converged
        print *,'self%available    : ',self%available
    end subroutine print_info

    subroutine kill( self )
        class(chunk), intent(inout) :: self
        self%id        = 0
        self%it        = 0
        call self%spproj%kill
        call self%qenv%kill
        self%nmics     = 0
        self%nptcls    = 0
        self%path      = ''
        self%projfile_out = ''
        if( allocated(self%orig_stks) ) deallocate(self%orig_stks)
        self%converged = .false.
        self%available = .false.
    end subroutine kill

    ! Main routine

    subroutine exec_cluster2D_stream( self, cline )
        use simple_class_frcs, only: class_frcs
        class(cluster2D_commander_stream), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        type(parameters)                   :: params
        type(rank_cavgs_commander)         :: xrank_cavgs
        type(cmdline)                      :: cline_cluster2D, cline_cluster2D_chunk, cline_rank_cavgs
        ! global pool variables
        type(qsys_env)                     :: qenv_pool
        type(sp_project)                   :: pool_proj
        character(LONGSTRLEN), allocatable :: imported_spprojs(:), imported_stks(:)
        character(len=STDLEN)              :: refs_glob, refs_glob_ranked
        logical,               allocatable :: pool_stack_mask(:)
        logical                            :: pool_converged, pool_available
        ! chunk-related variables
        type(chunk),           allocatable :: chunks(:), converged_chunks(:)
        integer                            :: n_new_chunks, glob_chunk_id, n_converged_chunks
        ! others
        type(sp_project)                   :: orig_proj, transfer_spproj
        character(len=:),      allocatable :: spproj_list_fname, orig_projfile
        character(len=LONGSTRLEN)          :: prev_snapshot_frcs, prev_snapshot_cavgs
        real    :: SMPD_TARGET = 4.       ! target sampling distance
        real    :: orig_smpd, scale_factor, smpd, lp_greedy, lpstart_stoch
        integer :: nptcls_per_chunk, ncls_glob, last_injection, max_ncls,ichunk, time_iter, ipart, iotest
        integer :: iter, orig_box, box, boxpd, n_spprojs, pool_iter, origproj_time, time_start_iter
        logical :: do_autoscale, l_greedy, l_once, l_restart
        if( cline%defined('refine') )then
            if( trim(cline%get_carg('refine')).ne.'greedy' )then
                if( .not.cline%defined('mskdiam') ) THROW_HARD('MSKDIAM must be defined!')
            endif
        else
            if( .not.cline%defined('mskdiam') ) THROW_HARD('MSKDIAM must be defined!')
        endif
        if( .not. cline%defined('mkdir')     ) call cline%set('mkdir',     'yes')
        if( .not. cline%defined('lp')        ) call cline%set('lp',          15.)
        if( .not. cline%defined('cenlp')     ) call cline%set('cenlp',       20.)
        if( .not. cline%defined('center')    ) call cline%set('center',     'no')
        if( .not. cline%defined('autoscale') ) call cline%set('autoscale', 'yes')
        if( .not. cline%defined('lpthresh')  ) call cline%set('lpthresh',    30.)
        if( .not. cline%defined('ndev')      ) call cline%set('ndev',        1.5)
        if( .not. cline%defined('oritype')   ) call cline%set('oritype','ptcl2D')
        call cline%set('ptclw',      'no')
        call cline%set('stream','yes') ! only for parameters determination
        call seed_rnd
        call params%new(cline)
        ! sanity
        if( .not.file_exists(params%projfile) )then
            THROW_HARD('project file: '//trim(params%projfile)//' does not exist!')
        endif
        orig_projfile = trim(params%projfile)
        if( .not.file_exists(params%dir_target) )then
            THROW_HARD('folder: '//trim(params%dir_target)//' does not exist!')
        endif
        call cline%set('stream','no') ! was only for parameters determination
        call cline%set('mkdir','no')
        ! permissions check
        call simple_touch(filepath(trim(params%dir_target),'iotest'),status=iotest)
        if( iotest /= 0 )then
            THROW_WARN(trim(params%dir_target)//' is not reachable!')
            THROW_HARD('Check permissions...')
        else
            call del_file(filepath(trim(params%dir_target),'iotest'))
        endif
        ! init
        do_autoscale      = params%autoscale .eq. 'yes'
        max_ncls          = floor(real(params%ncls)/real(params%ncls_start))*params%ncls_start ! effective maximum # of classes
        nptcls_per_chunk  = params%nptcls_per_cls*params%ncls_start         ! # of particles in each chunk
        spproj_list_fname = filepath(trim(params%dir_target), trim(STREAM_SPPROJFILES))
        ncls_glob         = 0
        l_greedy          = trim(params%refine).eq.'greedy'
        lp_greedy         = GREEDY_TARGET_LP
        if( cline%defined('lp') ) lp_greedy = params%lp
        SMPD_TARGET       = lp_greedy/2.
        lpstart_stoch     = (GREEDY_TARGET_LP+lp_greedy)/2.
        prev_snapshot_cavgs = ''
        call write_user_params
        call simple_mkdir(SCALE_DIR)
        ! call find_previous_snapshot
        ! init command-lines
        call cline%delete('lp')
        call cline%delete('refine')
        call cline%delete('nptcls_per_cls')
        call cline%delete('ncls_start')
        cline_cluster2D        = cline
        cline_cluster2D_chunk  = cline
        ! chunk classification
        ! Greedy optimisation, no match filter, no centering
        ! no incremental learning, objective function default is standard cross-correlation (cc)
        call cline_cluster2D_chunk%set('prg',       'cluster2D_distr')
        call cline_cluster2D_chunk%delete('projfile')
        call cline_cluster2D_chunk%delete('projname')
        call cline_cluster2D_chunk%set('objfun',    'cc')
        call cline_cluster2D_chunk%set('center',    'no')
        call cline_cluster2D_chunk%set('match_filt','no')
        call cline_cluster2D_chunk%set('autoscale', 'no')
        call cline_cluster2D_chunk%set('ptclw',     'no')
        call cline_cluster2D_chunk%set('mkdir',     'no')
        call cline_cluster2D_chunk%set('stream',    'no')
        call cline_cluster2D_chunk%delete('update_frac')
        call cline_cluster2D_chunk%delete('dir_target')
        if( l_greedy )then
            call cline_cluster2D_chunk%set('maxits', 10.)
            call cline_cluster2D_chunk%set('refine', 'greedy')
            call cline_cluster2D_chunk%set('lp',     GREEDY_TARGET_LP)
        else
            call cline_cluster2D_chunk%set('refine', 'snhc')
            call cline_cluster2D_chunk%set('extr_iter', real(MAX_EXTRLIM2D-2))
            call cline_cluster2D_chunk%set('maxits',  12.)
            call cline_cluster2D_chunk%set('lpstart', GREEDY_TARGET_LP)
            call cline_cluster2D_chunk%set('lpstop',  lp_greedy)
        endif
        call cline_cluster2D_chunk%set('startit',  1.)
        call cline_cluster2D_chunk%set('ncls',     real(params%ncls_start))
        call cline_cluster2D_chunk%set('nparts',   real(params%nparts_chunk))
        ! pool classification: optional stochastic optimisation, optional match filter
        ! automated incremental learning, objective function is standard cross-correlation (cc)
        call cline_cluster2D%set('prg',       'cluster2D_distr')
        call cline_cluster2D%set('autoscale', 'no')
        call cline_cluster2D%set('trs',       MINSHIFT)
        call cline_cluster2D%set('projfile',  trim(PROJFILE_POOL))
        call cline_cluster2D%set('projname',  trim(get_fbody(trim(PROJFILE_POOL),trim('simple'))))
        call cline_cluster2D%set('objfun',    'cc')
        call cline_cluster2D%set('ptclw',     'no')
        if( .not.cline%defined('match_filt') ) call cline_cluster2D%set('match_filt','no')
        call cline_cluster2D%set('extr_iter', 100.)
        call cline_cluster2D%set('mkdir',     'no')
        if( l_greedy )then
            call cline_cluster2D%set('center', 'no')
            call cline_cluster2D%set('refine', 'greedy')
            call cline_cluster2D%set('lp',     lp_greedy)
        else
            call cline_cluster2D%set('lpstart', lpstart_stoch)
        endif
        ! transfer project info
        call orig_proj%read(params%projfile)
        pool_proj%projinfo = orig_proj%projinfo
        pool_proj%compenv  = orig_proj%compenv
        if( orig_proj%jobproc%get_noris()>0 ) pool_proj%jobproc = orig_proj%jobproc
        call pool_proj%projinfo%delete_entry('projname')
        call pool_proj%projinfo%delete_entry('projfile')
        call pool_proj%update_projinfo(cline_cluster2D)
        transfer_spproj = pool_proj
        call transfer_spproj%write(PROJFILE_POOL)
        call orig_proj%kill
        ! initalize pool & chunks objects
        allocate(chunks(params_glob%nchunks))
        do ichunk = 1,params%nchunks
            call chunks(ichunk)%init(ichunk, pool_proj)
        enddo
        glob_chunk_id = params%nchunks
        ! we need to override the qsys_name for non local distributed execution
        call qenv_pool%new(params%nparts,exec_bin='simple_private_exec',qsys_name='local')
        ! wait for first stacks
        do
            if( file_exists(spproj_list_fname) )then
                if( .not.is_file_open(spproj_list_fname) )then
                    call generate_new_chunks(n_new_chunks)
                    if( n_new_chunks > 0 )then
                        exit ! Enough particles to initiate classification
                    endif
                endif
            endif
            call sleep(WAIT_WATCHER)
        enddo
        ! getting general parameters from the first sp_project
        orig_box     = chunks(1)%spproj%get_box()
        orig_smpd    = chunks(1)%spproj%get_smpd()
        if( orig_box == 0 ) THROW_HARD('FATAL ERROR')
        if( do_autoscale )then
            if( orig_box < MINBOXSZ )then
                do_autoscale = .false.
            else
                call autoscale(orig_box, orig_smpd, SMPD_TARGET, box, smpd, scale_factor)
                if( box < MINBOXSZ )then
                    SMPD_TARGET = orig_smpd * real(orig_box) / real(MINBOXSZ)
                    call autoscale(orig_box, orig_smpd, SMPD_TARGET, box, smpd, scale_factor)
                endif
                if( box == orig_box ) do_autoscale = .false.
            endif
        endif
        if( .not. do_autoscale )then
            smpd         = orig_smpd
            box          = orig_box
            scale_factor = 1.
        endif
        boxpd = 2 * round2even(params%alpha * real(box/2)) ! logics from parameters
        if( DEBUG_HERE )then
            print *,'box         ' ,box
            print *,'smpd        ' ,smpd
            print *,'scale_factor' ,scale_factor
        endif
        call cline_cluster2D_chunk%set('box',  real(box))
        call cline_cluster2D%set('box',  real(box)) ! strictly required
        call cline_cluster2D%set('smpd', smpd)      ! strictly required
        ! MAIN LOOP
        last_injection = simple_gettime()
        origproj_time  = last_injection
        pool_iter      = 0
        pool_converged = .false.
        pool_available = .true.
        call simple_touch(CLUSTER2D_FINISHED)
        iter = 0
        do
            iter = iter + 1
            time_start_iter = simple_gettime()
            call debug_print('global iter '//int2str(iter)//' '//trim(cast_time_char(time_start_iter)))
            ! update rejection parameters
            call update_user_params
            if( file_exists(TERM_STREAM) )then
                write(logfhandle,'(A,A)')'>>> TERMINATING CLUSTER2D STREAM ',trim(cast_time_char(simple_gettime()))
                exit
            endif
            ! classify chunks
            call debug_print('chunk section global iter '//int2str(iter))
            do ichunk = 1,params%nchunks
                call chunks(ichunk)%exec_classify(cline_cluster2D_chunk, orig_smpd, orig_box, box)
            enddo
            ! deal with chunk completion, rejection, reset
            do ichunk = 1,params%nchunks
                if( chunks(ichunk)%has_converged() )then
                    call chunks(ichunk)%display_iter
                    ! rejection
                    call chunks(ichunk)%reject(params%lpthresh, params%ndev, box)
                    ! updates list of chunks to import
                    if( allocated(converged_chunks) )then
                        converged_chunks = [converged_chunks(:), chunks(ichunk)]
                    else
                        allocate(converged_chunks(1),source=[chunks(ichunk)])
                    endif
                    ! free chunk
                    glob_chunk_id = glob_chunk_id + 1
                    call chunks(ichunk)%init(glob_chunk_id, pool_proj)
                endif
            enddo
            n_converged_chunks = 0
            if( allocated(converged_chunks) ) n_converged_chunks = size(converged_chunks)
            call debug_print('end chunk section global iter '//int2str(iter))
            ! deal with pool completion, rejection, execution
            call update_user_params
            if( .not.pool_available )then
                call debug_print('pool unavailable '//int2str(iter))
                pool_converged = file_exists(CLUSTER2D_FINISHED)
                if( pool_iter > 1 ) refs_glob = trim(CAVGS_ITER_FBODY)//trim(int2str_pad(pool_iter,3))//trim(params%ext)
                pool_available = pool_converged
            endif
            if( pool_available )then
                call debug_print('pool available '//int2str(iter))
                call del_file(CLUSTER2D_FINISHED)
                ! read in previous iteration
                call update_pool
                ! reject
                if( pool_iter > 2*FREQ_POOL_REJECTION .and. mod(pool_iter,FREQ_POOL_REJECTION)==0 )then
                    call reject_from_pool
                endif
                ! writes spnapshot every ORIGPROJ_WRITEFREQ seconds
                if( (simple_gettime()-origproj_time) > ORIGPROJ_WRITEFREQ .and. pool_iter>1)then
                    call write_snapshot(.true.)
                    origproj_time = simple_gettime()
                endif
                if( file_exists(SPPROJ_SNAPSHOT) )then
                    call write_snapshot(.true.)
                    call del_file(SPPROJ_SNAPSHOT)
                endif
                ! append new chunks
                if( n_converged_chunks > 0) call import_chunks_into_pool
                ! execute
                call exec_classify_pool
                ! tidy
                if( pool_iter > 3 )then
                    call del_file(trim(CAVGS_ITER_FBODY)//trim(int2str_pad(pool_iter-3,3))//'_even'//trim(params%ext))
                    call del_file(trim(CAVGS_ITER_FBODY)//trim(int2str_pad(pool_iter-3,3))//'_odd'//trim(params%ext))
                endif
                call debug_print('end pool available '//int2str(iter))
            endif
            call debug_print('new chunk section global iter '//int2str(iter))
            ! generate new chunk projects
            call generate_new_chunks(n_new_chunks)
            if( n_new_chunks > 0 ) last_injection = simple_gettime()
            ! wait a bit if necessary
            time_iter = simple_gettime() - time_start_iter
            if( time_iter < WAIT_WATCHER ) call sleep(WAIT_WATCHER-time_iter)
            ! optionally pause
            l_once = .false.
            do while( file_exists(trim(PAUSE_STREAM)) )
                if( file_exists(trim(TERM_STREAM)) ) exit
                if( .not.l_once )then
                    l_once = .true.
                    call write_singlelineoftext(PAUSE_STREAM, 'PAUSED')
                    write(logfhandle,'(A,A)')'>>> CLUSTER2D STREAM PAUSED ',cast_time_char(simple_gettime())
                endif
                call sleep(WAIT_WATCHER)
            enddo
            call debug_print('end new chunk section global iter '//int2str(iter))
            ! endif
        enddo
        call debug_print('exited global iter '//int2str(iter))
        ! call qsys_cleanup
        do ichunk = 1,params%nchunks
            call chunks(ichunk)%terminate
        enddo
        if( .not.pool_converged )then
            pool_iter = pool_iter-1
            refs_glob = trim(CAVGS_ITER_FBODY)//trim(int2str_pad(pool_iter,3))//trim(params%ext)
            ! tricking the asynchronous master process to come to a hard stop
            do ipart = 1,params%nparts
                call simple_touch('JOB_FINISHED_'//int2str_pad(ipart,params%numlen))
            enddo
        endif
        if( pool_iter >= 1 )then
            ! updates project
            call write_snapshot(.false.)
            ! ranking
            refs_glob_ranked = add2fbody(refs_glob,params%ext,'_ranked')
            call cline_rank_cavgs%set('projfile', orig_projfile)
            call cline_rank_cavgs%set('stk',      refs_glob)
            call cline_rank_cavgs%set('outstk',   trim(refs_glob_ranked))
            call xrank_cavgs%execute(cline_rank_cavgs)
        endif
        ! cleanup
        if( .not.debug_here )then
            ! call qsys_cleanup
            call simple_rmdir(SCALE_DIR)
            call del_file(PROJFILE_POOL)
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_CLUSTER2D_STREAM NORMAL STOP ****')

        contains

            subroutine update_pool
                type(sp_project)      :: spproj
                type(oris)            :: os
                integer, allocatable  :: pops(:)
                character(len=STDLEN) :: fname
                real    :: corr, frac, mi_class
                integer :: i, it, jptcl, iptcl, istk
                call debug_print('in update_pool pool_iter '//int2str(pool_iter))
                ! iteration info
                fname = trim(STATS_FILE)
                if( file_exists(fname) )then
                    call os%new(1,is_ptcl=.false.)
                    call os%read(fname)
                    it = nint(os%get(1,'ITERATION'))
                    if( it == pool_iter )then
                        mi_class = os%get(1,'CLASS_OVERLAP')
                        frac     = os%get(1,'SEARCH_SPACE_SCANNED')
                        corr     = os%get(1,'CORRELATION')
                        write(logfhandle,'(A,I6,A,F7.3,A,F7.3,A,F7.3)')'>>> POOL         ITERATION ',it,&
                            &'; CLASS OVERLAP: ',mi_class,'; SEARCH SPACE SCANNED: ',frac,'; CORRELATION: ',corr
                    endif
                    call os%kill
                endif
                call debug_print('in update_pool 1')
                ! transfer to pool
                call spproj%read_segment('stk',PROJFILE_POOL)
                call spproj%read_segment('ptcl2D',PROJFILE_POOL)
                call spproj%read_segment('cls2D', PROJFILE_POOL)
                if( spproj%os_cls2D%get_noris() == 0 )then
                    call debug_print('in update_pool 2')
                    ! not executed yet, do nothing
                else
                    call debug_print('in update_pool 3')
                    if( .not.allocated(pool_stack_mask) )then
                        ! first time
                        THROW_HARD('Critical ERROR 0')
                    endif
                    i = 0
                    do istk = 1,size(pool_stack_mask)
                        if( pool_stack_mask(istk) )then
                            i = i+1
                            iptcl = nint(pool_proj%os_stk%get(istk,'fromp'))
                            do jptcl = nint(spproj%os_stk%get(i,'fromp')),nint(spproj%os_stk%get(i,'top'))
                                if( spproj%os_ptcl2D%get_state(jptcl) > 0 )then
                                    call pool_proj%os_ptcl2D%transfer_2Dparams(iptcl, spproj%os_ptcl2D, jptcl)
                                endif
                                iptcl = iptcl+1
                            enddo
                        endif
                    enddo
                    call pool_proj%os_ptcl2D%get_pops(pops, 'class', consider_w=.false., maxn=ncls_glob)
                    pool_proj%os_cls2D = spproj%os_cls2D
                    call pool_proj%os_cls2D%set_all('pop', real(pops))
                    call debug_print('in update_pool 4')
                    ! for gui
                    os = pool_proj%os_cls2D
                    call pool_proj%add_cavgs2os_out(refs_glob, smpd, 'cavg')
                    pool_proj%os_cls2D = os
                    call pool_proj%write_segment_inside('out',   orig_projfile)
                    call pool_proj%write_segment_inside('cls2D', orig_projfile)
                    call os%kill
                endif
                call spproj%kill
                call debug_print('end update_pool')
            end subroutine update_pool

            !> returns the list of projects necessary for creating a new chunk
            subroutine generate_new_chunks( nnewchunks )
                integer,                 intent(inout) :: nnewchunks
                type(sp_project)                       :: stream_proj
                character(len=LONGSTRLEN), allocatable :: tmp(:), spprojs_for_chunk(:), spproj_list(:)
                integer,                   allocatable :: spproj_nptcls(:)
                logical,                   allocatable :: spproj_mask(:)
                integer :: nptcls, iproj, jproj, cnt, nmics_imported, iichunk, first_new, n_new, n_avail, n2fill, iostat
                logical :: isnew, enough
                call debug_print('in generate_new_chunks')
                nmics_imported = 0
                if( allocated(imported_spprojs) ) nmics_imported = size(imported_spprojs)
                nnewchunks = 0
                ! whether any chunk is available
                n_avail = 0
                do iichunk =  1,params%nchunks
                    if (chunks(iichunk)%available)then
                        n_avail = n_avail+1
                    endif
                enddo
                call debug_print('in generate_new_chunks 1 '//int2str(n_avail)//' '//int2str(nmics_imported))
                if(n_avail == 0) return
                call debug_print('in generate_new_chunks reading '//trim(spproj_list_fname))
                call simple_copy_file(spproj_list_fname, trim(STREAM_SPPROJFILES), iostat)
                call fileiochk('failure to copy '//spproj_list_fname, iostat)
                call read_filetable(trim(STREAM_SPPROJFILES), spproj_list)
                if( .not.allocated(spproj_list) )return
                n_spprojs = size(spproj_list)
                call debug_print('in generate_new_chunks 1b  '//int2str(n_spprojs))
                ! whether the projects have been processed
                allocate(spproj_mask(n_spprojs),source=.false.)
                do iproj = 1,n_spprojs
                    isnew = .true.
                    do jproj = 1,nmics_imported
                        if( trim(spproj_list(iproj)) .eq. trim(imported_spprojs(jproj)) )then
                            isnew = .false.
                            exit
                        endif
                    enddo
                    spproj_mask(iproj) = isnew
                enddo
                call debug_print('in generate_new_chunks 2 '//int2str(count(spproj_mask)))
                ! first new project index
                first_new = 0
                do iproj = 1,n_spprojs
                    if( spproj_mask(iproj) )then
                        first_new = iproj
                        exit
                    endif
                enddo
                call debug_print('in generate_new_chunks 3 '//int2str(first_new))
                ! no new data to process
                if( first_new == 0 ) return
                ! gather number of particles
                allocate(spproj_nptcls(n_spprojs),source=0)
                nptcls = 0
                n2fill = 0
                do iproj = first_new,n_spprojs
                    if( spproj_mask(iproj) )then
                        call stream_proj%read_segment('mic',spproj_list(iproj))
                        spproj_nptcls(iproj) = nint(stream_proj%os_mic%get(1,'nptcls'))
                        if( spproj_nptcls(iproj)==0 ) THROW_HARD('zero particles project detected!')
                        nptcls = nptcls + spproj_nptcls(iproj)
                        if( nptcls > nptcls_per_chunk )then
                            n2fill = n2fill + 1
                            if( n2fill >= n_avail )exit
                            nptcls = 0
                        endif
                    endif
                enddo
                call stream_proj%kill
                call debug_print('in generate_new_chunks 4 '//int2str(n2fill))
                ! enough ?
                if( n2fill == 0 )then
                    return
                endif
                do iichunk =  1,params%nchunks
                    if( .not.chunks(iichunk)%available ) cycle
                    do iproj = first_new,n_spprojs
                        nptcls = sum(spproj_nptcls(first_new:iproj))
                        enough = nptcls >= nptcls_per_chunk
                        if( enough )exit
                    enddo
                    if( .not.enough ) exit
                    ! generate new chunk
                    n_new = iproj - first_new + 1
                    allocate(spprojs_for_chunk(n_new))
                    cnt = 0
                    do jproj = first_new,iproj
                        cnt = cnt + 1
                        spprojs_for_chunk(cnt) = trim(spproj_list(jproj))
                    enddo
                    call chunks(iichunk)%generate(spprojs_for_chunk, nptcls, first_new)
                    nnewchunks = nnewchunks + 1
                    ! update list of imported projects
                    if( nmics_imported == 0 )then
                        allocate(imported_spprojs(n_new))
                    else
                        call move_alloc(imported_spprojs, tmp)
                        allocate(imported_spprojs(nmics_imported+n_new))
                        imported_spprojs(1:nmics_imported) = tmp(:)
                        deallocate(tmp)
                    endif
                    imported_spprojs(nmics_imported+1:nmics_imported+n_new) = spprojs_for_chunk(:)
                    nmics_imported = nmics_imported + n_new
                    deallocate(spprojs_for_chunk)
                    ! update for next chunk
                    spproj_nptcls(first_new:iproj) = 0
                    spproj_mask(first_new:iproj)   = .false.
                    first_new = min(iproj+1,n_spprojs)
                    call debug_print('in generate_new_chunks 5 '//int2str(iichunk))
                enddo
                call debug_print('end generate_new_chunks '//int2str(nmics_imported))
            end subroutine generate_new_chunks

            subroutine import_chunks_into_pool
                type(class_frcs)              :: frcs_glob, frcs_chunk, frcs_prev
                character(LONGSTRLEN), allocatable :: tmp(:)
                character(len=:),      allocatable :: cavgs_chunk, dir_chunk
                integer,               allocatable :: cls_pop(:), cls_chunk_pop(:), pinds(:), states(:)
                type(oris) :: os_backup
                character(len=LONGSTRLEN) :: stk_relpath
                real    :: smpd_here
                integer :: nptcls_imported, nptcls2import
                integer :: iptcl, ind, ncls_here, icls, i, poolind, n_remap, pop, iichunk, nmics2import, nmics_imported
                integer :: nptcls, imic, iproj, ncls_tmp, fromp_prev, nchunks2import, fromp, ii, jptcl
                logical :: l_maxed
                call debug_print('in import_chunks_into_pool')
                nptcls2import   = 0
                nmics2import    = 0
                nchunks2import = size(converged_chunks)
                do iichunk = 1,nchunks2import
                    nptcls2import = nptcls2import + converged_chunks(iichunk)%nptcls
                    nmics2import  = nmics2import  + converged_chunks(iichunk)%nmics
                enddo
                if( nptcls2import == 0 ) return
                call debug_print('in import_chunks_into_pool 1 '//int2str(nptcls2import))
                ! reallocations
                nmics_imported  = pool_proj%os_mic%get_noris()
                nptcls_imported = pool_proj%os_ptcl2D%get_noris()
                if( nmics_imported == 0 )then
                    call pool_proj%os_mic%new(nmics2import, is_ptcl=.false.)
                    call pool_proj%os_stk%new(nmics2import, is_ptcl=.false.)
                    call pool_proj%os_ptcl2D%new(nptcls2import, is_ptcl=.true.)
                    allocate(imported_stks(nmics2import))
                    fromp = 1
                else
                    call pool_proj%os_mic%reallocate(nmics_imported+nmics2import)
                    call pool_proj%os_stk%reallocate(nmics_imported+nmics2import)
                    call pool_proj%os_ptcl2D%reallocate(nptcls_imported+nptcls2import)
                    fromp = nint(pool_proj%os_stk%get(nmics_imported,'top'))+1
                    call move_alloc(imported_stks, tmp)
                    allocate(imported_stks(nmics_imported+nmics2import))
                    imported_stks(1:nmics_imported) = tmp(:)
                    deallocate(tmp)
                endif
                imic  = nmics_imported
                iptcl = nptcls_imported
                call debug_print('in import_chunks_into_pool 2 '//int2str(imic)//' '//int2str(iptcl))
                do iichunk = 1,nchunks2import
                    fromp_prev = fromp
                    call converged_chunks(iichunk)%read(boxpd)
                    ! transfer micrographs, stacks & particles parameters
                    jptcl = 0
                    do iproj=1,converged_chunks(iichunk)%nmics
                        imic  = imic+1
                        ! micrographs
                        call pool_proj%os_mic%transfer_ori(imic, converged_chunks(iichunk)%spproj%os_mic, iproj)
                        ! stacks
                        call pool_proj%os_stk%transfer_ori(imic, converged_chunks(iichunk)%spproj%os_stk, iproj)
                        nptcls = nint(converged_chunks(iichunk)%spproj%os_stk%get(iproj,'nptcls'))
                        call pool_proj%os_stk%set(imic, 'fromp', real(fromp))
                        call pool_proj%os_stk%set(imic, 'top',   real(fromp+nptcls-1))
                        call make_relativepath(cwd_glob, converged_chunks(iichunk)%orig_stks(iproj), stk_relpath)
                        imported_stks(imic) = trim(stk_relpath)
                        ! particles
                        do i = 1,nptcls
                            iptcl = iptcl + 1
                            jptcl = jptcl + 1
                            call pool_proj%os_ptcl2D%transfer_ori(iptcl, converged_chunks(iichunk)%spproj%os_ptcl2D, jptcl)
                            call pool_proj%os_ptcl2D%set(iptcl, 'stkind',    real(imic))
                            call pool_proj%os_ptcl2D%set(iptcl, 'updatecnt', 0.) ! new particle
                            call pool_proj%os_ptcl2D%set(iptcl, 'frac',      0.) ! new particle
                        enddo
                        fromp = fromp + nptcls
                    enddo
                    write(logfhandle,'(A,I6,A,I6,A)')'>>> TRANSFERRED ',fromp-fromp_prev,' PARTICLES FROM CHUNK ',&
                    &converged_chunks(iichunk)%id,' TO POOL'
                    call flush(logfhandle)
                    ! transfer classes
                    l_maxed   = ncls_glob >= max_ncls ! max # of classes reached ?
                    ncls_here = ncls_glob
                    if( .not.l_maxed ) ncls_here = ncls_glob+params%ncls_start
                    states = nint(converged_chunks(iichunk)%spproj%os_ptcl2D%get_all('state'))
                    call converged_chunks(iichunk)%spproj%get_cavgs_stk(cavgs_chunk, ncls_tmp, smpd_here)
                    dir_chunk   = trim(converged_chunks(iichunk)%path)
                    cavgs_chunk = trim(dir_chunk)//basename(cavgs_chunk)
                    call frcs_chunk%new(params%ncls_start, box, smpd, nstates=1)
                    call frcs_chunk%read(trim(dir_chunk)//trim(FRCS_FILE))
                    if( l_maxed )then
                        call debug_print('in import_chunks_into_pool 3 '//int2str(iichunk))
                        ! transfer all others
                        cls_pop = nint(pool_proj%os_cls2D%get_all('pop'))
                        n_remap = 0
                        if( any(cls_pop==0) )then
                            if( all(cls_pop==0) ) THROW_HARD('Empty os_cls2D!')
                            ! remapping
                            cls_chunk_pop = nint(converged_chunks(iichunk)%spproj%os_cls2D%get_all('pop'))
                            call frcs_glob%new(ncls_glob, box, smpd, nstates=1)
                            call frcs_glob%read(FRCS_FILE)
                            do icls=1,ncls_glob
                                if( cls_pop(icls)>0 ) cycle          ! class already filled
                                if( all(cls_chunk_pop == 0 ) ) exit ! no more chunk class available
                                ind = irnd_uni(params%ncls_start)    ! selects chunk class stochastically
                                do while( cls_chunk_pop(ind) == 0 )
                                    ind = irnd_uni(params%ncls_start)
                                enddo
                                cls_chunk_pop(ind) = 0             ! excludes from being picked again
                                n_remap = n_remap+1
                                ! class average
                                call transfer_cavg(cavgs_chunk, dir_chunk, ind, refs_glob, icls)
                                ! frcs
                                call frcs_glob%set_frc(icls,frcs_chunk%get_frc(ind, box, 1), 1)
                                ! class parameters transfer
                                call pool_proj%os_cls2D%transfer_ori(icls, converged_chunks(iichunk)%spproj%os_cls2D, ind)
                                call pool_proj%os_cls2D%set(icls, 'class', real(icls))
                                ! particles
                                call converged_chunks(iichunk)%spproj%os_ptcl2D%get_pinds(ind,'class',pinds,consider_w=.false.)
                                pop = size(pinds)
                                do i=1,pop
                                    ii      = pinds(i)               ! in chunk
                                    poolind = fromp_prev + ii - 1 ! in pool
                                    call pool_proj%os_ptcl2D%set(poolind,'class',real(icls))
                                enddo
                                cls_pop(icls) = cls_pop(icls) + pop ! updates class populations
                            enddo
                            ! now transfer particles that were not remapped
                            do icls = 1,params%ncls_start
                                if( cls_chunk_pop(icls) == 0 ) cycle
                                call converged_chunks(iichunk)%spproj%os_ptcl2D%get_pinds(icls,'class',pinds,consider_w=.false.)
                                pop = size(pinds)
                                do i = 1,pop
                                    ii      = pinds(i)            ! in chunk
                                    poolind = fromp_prev + ii - 1 ! in pool
                                    ind     = irnd_uni(ncls_glob) ! stochastic labelling followed by greedy search
                                    do while( cls_pop(ind) == 0 )
                                        ind = irnd_uni(ncls_glob)
                                    enddo
                                    call pool_proj%os_ptcl2D%set(poolind,'class',real(ind))
                                    cls_pop(ind) = cls_pop(ind) + 1 ! updates class populations
                                enddo
                            enddo
                            if( n_remap > 0 )then
                                call transfer_cavg(refs_glob,dir_chunk,ncls_glob,refs_glob,ncls_glob, self_transfer=.true.)
                                call frcs_glob%write(FRCS_FILE)
                                write(logfhandle,'(A,I4)')'>>> # OF RE-MAPPED CLASS AVERAGES: ',n_remap
                            endif
                        else
                            call debug_print('in import_chunks_into_pool 4 '//int2str(iichunk))
                            ! no remapping, just transfer particles & updates 2D population
                            do ii = 1,converged_chunks(iichunk)%nptcls
                                if( states(ii) /= 0 )then
                                    icls          = irnd_uni(ncls_glob) ! stochastic labelling followed by greedy search
                                    cls_pop(icls) = cls_pop(icls) + 1   ! updates populations
                                    poolind       = fromp_prev + ii - 1
                                    call pool_proj%os_ptcl2D%set(poolind, 'class', real(icls))
                                endif
                            enddo
                        endif
                        ! updates class populations
                        call pool_proj%os_cls2D%set_all('pop',real(cls_pop))
                        call debug_print('in import_chunks_into_pool 5 '//int2str(iichunk))
                    else
                        call debug_print('in import_chunks_into_pool 6 '//' '//int2str(iichunk)//' '//int2str(ncls_glob))
                        ! all new classes can be imported, no remapping
                        if( ncls_glob == 0 )then
                            ! first transfer : copy classes, frcs & class parameters
                            refs_glob = 'start_cavgs'//params%ext
                            do icls= 1,params%ncls_start
                                call transfer_cavg(cavgs_chunk,dir_chunk,icls,refs_glob,icls)
                            enddo
                            call simple_copy_file(trim(dir_chunk)//trim(FRCS_FILE),FRCS_FILE)
                            pool_proj%os_cls2D = converged_chunks(iichunk)%spproj%os_cls2D
                        else
                            ! append new classes
                            do icls=1,params%ncls_start
                                call transfer_cavg(cavgs_chunk, dir_chunk, icls, refs_glob, ncls_glob+icls)
                            enddo
                            ! FRCs
                            call frcs_glob%new(ncls_here, box, smpd, nstates=1)
                            call frcs_prev%new(ncls_glob, box, smpd, nstates=1)
                            call frcs_prev%read(FRCS_FILE)
                            do icls=1,ncls_glob
                                call frcs_glob%set_frc(icls,frcs_prev%get_frc(icls, box, 1), 1)
                            enddo
                            do icls=1,params%ncls_start
                                call frcs_glob%set_frc(ncls_glob+icls,frcs_chunk%get_frc(icls, box, 1), 1)
                            enddo
                            call frcs_glob%write(FRCS_FILE)
                            ! class parameters
                            call pool_proj%os_cls2D%reallocate(ncls_here)
                            do icls = 1,params%ncls_start
                                ind = ncls_glob+icls
                                call pool_proj%os_cls2D%transfer_ori(ind, converged_chunks(iichunk)%spproj%os_cls2D, icls)
                                call pool_proj%os_cls2D%set(ind, 'class', real(ind))
                            enddo
                        endif
                        call debug_print('in import_chunks_into_pool 7'//' '//int2str(iichunk))
                        ! particles 2D
                        do ii = 1,converged_chunks(iichunk)%nptcls
                            if( states(ii) /= 0 )then
                                poolind = fromp_prev+ii-1
                                icls    = ncls_glob+nint(converged_chunks(iichunk)%spproj%os_ptcl2D%get(ii,'class'))
                                call pool_proj%os_ptcl2D%set(poolind, 'class', real(icls))
                            endif
                        enddo
                        ! global # of classes
                        ncls_glob = ncls_here
                        call debug_print('in import_chunks_into_pool 8'//' '//int2str(iichunk))
                    endif
                    ! remove chunk
                    call converged_chunks(iichunk)%remove_folder
                    call converged_chunks(iichunk)%kill
                    call debug_print('in import_chunks_into_pool 9'//' '//int2str(iichunk))
                enddo
                ! for gui
                os_backup = pool_proj%os_cls2D
                call pool_proj%add_cavgs2os_out(refs_glob, smpd, 'cavg')
                pool_proj%os_cls2D = os_backup
                call pool_proj%write_segment_inside('out',   orig_projfile)
                call pool_proj%write_segment_inside('cls2D', orig_projfile)
                call os_backup%kill
                call debug_print('in import_chunks_into_pool 10'//' '//int2str(iichunk))
                ! cleanup
                deallocate(converged_chunks)
                call frcs_prev%kill
                call frcs_glob%kill
                call frcs_chunk%kill
                call debug_print('end import_chunks_into_pool')
            end subroutine import_chunks_into_pool

            subroutine exec_classify_pool
                use simple_ran_tabu
                type(ran_tabu)          :: random_generator
                logical, parameter      :: L_BENCH = .false.
                integer, allocatable    :: prev_eo_pops(:,:), min_update_cnts_per_stk(:), nptcls_per_stk(:), stk_order(:)
                real                    :: frac_update
                integer                 :: iptcl,i, nptcls_tot, nptcls_old, fromp, top, nstks_tot, jptcl
                integer                 :: eo, icls, nptcls_sel, istk, nptcls2update, nstks2update
                integer(timer_int_kind) :: t_tot
                if( L_BENCH )then
                    t_tot  = tic()
                endif
                nptcls_tot = pool_proj%os_ptcl2D%get_noris()
                if( nptcls_tot == 0 ) return
                call debug_print('in exec_classify_pool '//int2str(nptcls_tot))
                pool_iter = pool_iter + 1
                call cline_cluster2D%set('refs',    refs_glob)
                call cline_cluster2D%set('ncls',    real(ncls_glob))
                call cline_cluster2D%set('startit', real(pool_iter))
                call cline_cluster2D%set('maxits',  real(pool_iter))
                call cline_cluster2D%set('frcs',    trim(FRCS_FILE))
                transfer_spproj%projinfo = pool_proj%projinfo
                transfer_spproj%compenv  = pool_proj%compenv
                call transfer_spproj%projinfo%delete_entry('projname')
                call transfer_spproj%projinfo%delete_entry('projfile')
                call transfer_spproj%update_projinfo( cline_cluster2D )
                ! counting number of stacks & particles (old and newer)
                nstks_tot  = pool_proj%os_stk%get_noris()
                allocate(nptcls_per_stk(nstks_tot), min_update_cnts_per_stk(nstks_tot), source=0)
                nptcls_old = 0
                do istk = 1,nstks_tot
                    fromp = nint(pool_proj%os_stk%get(istk,'fromp'))
                    top   = nint(pool_proj%os_stk%get(istk,'top'))
                    min_update_cnts_per_stk(istk) = huge(istk)
                    do iptcl = fromp,top
                        if( pool_proj%os_ptcl2D%get(iptcl,'state') > 0.5 )then
                            nptcls_per_stk(istk)          = nptcls_per_stk(istk) + 1 ! # ptcls with state=1
                            min_update_cnts_per_stk(istk) = min(min_update_cnts_per_stk(istk), nint(pool_proj%os_ptcl2D%get(iptcl,'updatecnt')))
                        endif
                    enddo
                    if( min_update_cnts_per_stk(istk) >= STREAM_SRCHLIM )then
                        nptcls_old = nptcls_old + nptcls_per_stk(istk)
                    endif
                enddo
                ! flagging stacks to be skipped
                if( allocated(pool_stack_mask) ) deallocate(pool_stack_mask)
                allocate(pool_stack_mask(nstks_tot), source=.false.)
                allocate(prev_eo_pops(ncls_glob,2),source=0)
                if( nptcls_old > params%ncls_start*params%nptcls_per_cls )then
                    allocate(stk_order(nstks_tot))
                    do istk = 1,nstks_tot
                        stk_order(istk) = istk
                    enddo
                    random_generator = ran_tabu(nstks_tot)
                    call random_generator%shuffle(stk_order)
                    nptcls2update = 0 ! # of ptcls including state=0 within selected stacks
                    nptcls_sel    = 0 ! # of ptcls excluding state=0 within selected stacks
                    do i = 1,nstks_tot
                        istk = stk_order(i)
                        if( (min_update_cnts_per_stk(istk) > STREAM_SRCHLIM) .and. (nptcls_sel > MAX_STREAM_NPTCLS) ) cycle
                        nptcls_sel    = nptcls_sel + nptcls_per_stk(istk)
                        nptcls2update = nptcls2update + nint(pool_proj%os_stk%get(istk,'nptcls'))
                        pool_stack_mask(istk) = .true.
                    enddo
                    call random_generator%kill
                else
                    nptcls2update   = nptcls_tot
                    nptcls_sel      = sum(nptcls_per_stk)
                    do istk = 1,nstks_tot
                        pool_stack_mask(istk) = nptcls_per_stk(istk) > 0
                    enddo
                endif
                nstks2update = count(pool_stack_mask)
                ! transfer stacks and particles
                call transfer_spproj%os_stk%new(nstks2update, is_ptcl=.false.)
                call transfer_spproj%os_ptcl2D%new(nptcls2update, is_ptcl=.true.)
                i     = 0
                jptcl = 0
                do istk = 1,nstks_tot
                    fromp = nint(pool_proj%os_stk%get(istk,'fromp'))
                    top   = nint(pool_proj%os_stk%get(istk,'top'))
                    if( pool_stack_mask(istk) )then
                        ! transfer alignement parameters for selected particles
                        i = i + 1 ! stack index in transfer_spproj
                        call transfer_spproj%os_stk%transfer_ori(i, pool_proj%os_stk, istk)
                        call transfer_spproj%os_stk%set(i, 'fromp', real(jptcl+1))
                        do iptcl = fromp,top
                            jptcl = jptcl+1
                            call transfer_spproj%os_ptcl2D%transfer_ori(jptcl, pool_proj%os_ptcl2D, iptcl)
                            call transfer_spproj%os_ptcl2D%set(jptcl, 'stkind', real(i))
                        enddo
                        call transfer_spproj%os_stk%set(i, 'top', real(jptcl))
                    else
                        ! keeps track of skipped particles
                        do iptcl = fromp,top
                            icls = nint(pool_proj%os_ptcl2D%get(iptcl,'class'))
                            eo   = nint(pool_proj%os_ptcl2D%get(iptcl,'eo')) + 1
                            prev_eo_pops(icls,eo) = prev_eo_pops(icls,eo) + 1
                        enddo
                    endif
                enddo
                call transfer_spproj%os_ptcl3D%new(nptcls2update, is_ptcl=.true.)
                transfer_spproj%os_cls2D = pool_proj%os_cls2D
                ! execution
                if( sum(prev_eo_pops) == 0 )then
                    call cline_cluster2D%set('stream','no')
                    call cline_cluster2D%delete('update_frac')
                else
                    frac_update = real(nptcls_old-sum(prev_eo_pops)) / real(nptcls_old)
                    call cline_cluster2D%set('stream',      'yes')
                    call cline_cluster2D%set('update_frac', frac_update)
                    call cline_cluster2D%set('center',      'no')
                    do icls = 1,ncls_glob
                        call transfer_spproj%os_cls2D%set(icls,'prev_pop_even',real(prev_eo_pops(icls,1)))
                        call transfer_spproj%os_cls2D%set(icls,'prev_pop_odd', real(prev_eo_pops(icls,2)))
                    enddo
                endif
                call transfer_spproj%write(PROJFILE_POOL)
                call transfer_spproj%kill
                call debug_print('in exec_classify_pool 4')
                call qenv_pool%exec_simple_prg_in_queue_async(cline_cluster2D, DISTR_EXEC_FNAME, 'simple_log_cluster2D_pool')
                pool_available = .false.
                pool_converged = .false.
                write(logfhandle,'(A,I6,A,I8,A3,I8,A)')'>>> POOL         INITIATED ITERATION ',pool_iter,' WITH ',nptcls_sel,&
                &' / ', sum(nptcls_per_stk),' PARTICLES'
                if( L_BENCH ) print *,'timer exec_classify_pool tot : ',toc(t_tot)
                call debug_print('end exec_classify_pool')
            end subroutine exec_classify_pool

            subroutine rescale_cavgs( src, dest )
                character(len=*), intent(in) :: src, dest
                type(image)          :: img, img_pad
                integer, allocatable :: cls_pop(:)
                integer              :: ldim(3),icls, ncls_here
                call debug_print('in rescale_cavgs '//trim(src))
                call debug_print('in rescale_cavgs '//trim(dest))
                call img%new([box,box,1],smpd)
                call img_pad%new([orig_box,orig_box,1],params%smpd)
                cls_pop = nint(pool_proj%os_cls2D%get_all('pop'))
                call find_ldim_nptcls(src,ldim,ncls_here)
                do icls = 1,ncls_here
                    if( cls_pop(icls) > 0 )then
                        call img%zero_and_unflag_ft
                        call img%read(src,icls)
                        call img%fft
                        call img%pad(img_pad, backgr=0.)
                        call img_pad%ifft
                    else
                        img_pad = 0.
                    endif
                    call img_pad%write('tmp_cavgs.mrc',icls)
                enddo
                call simple_rename('tmp_cavgs.mrc',dest)
                call img%kill
                call img_pad%kill
                call debug_print('end rescale_cavgs')
            end subroutine rescale_cavgs

            subroutine reject_from_pool
                type(image)          :: img
                logical, allocatable :: cls_mask(:)
                real                 :: ndev_here
                integer              :: nptcls_rejected, ncls_rejected, iptcl
                integer              :: icls, cnt
                call debug_print('in reject from_pool')
                if( pool_proj%os_cls2D%get_noris() == 0 ) return
                ncls_rejected   = 0
                nptcls_rejected = 0
                allocate(cls_mask(ncls_glob), source=.true.)
                ! correlation & resolution
                ndev_here = 1.5*params%ndev ! less stringent rejection
                call pool_proj%os_cls2D%find_best_classes(box,smpd,params%lpthresh,cls_mask,ndev_here)
                if( count(cls_mask) > 1 .and. count(cls_mask) < ncls_glob )then
                    ncls_rejected = 0
                    do iptcl=1,pool_proj%os_ptcl2D%get_noris()
                        if( pool_proj%os_ptcl2D%get_state(iptcl) == 0 )cycle
                        icls = nint(pool_proj%os_ptcl2D%get(iptcl,'class'))
                        if( cls_mask(icls) ) cycle
                        nptcls_rejected = nptcls_rejected+1
                        call pool_proj%os_ptcl2D%set(iptcl,'state',0.)
                    enddo
                    call debug_print('in reject from_pool 1')
                    if( nptcls_rejected > 0 )then
                        do icls=1,ncls_glob
                            if( .not.cls_mask(icls) )then
                                if( pool_proj%os_cls2D%get(icls,'pop') > 0.5 ) ncls_rejected = ncls_rejected+1
                                call pool_proj%os_cls2D%set(icls,'pop',0.)
                                call pool_proj%os_cls2D%set(icls,'corr',-1.)
                            endif
                        enddo
                        call debug_print('in reject from_pool 2')
                        cnt = 0
                        call img%new([box,box,1],smpd)
                        do icls=1,ncls_glob
                            if( cls_mask(icls) ) cycle
                            cnt = cnt+1
                            if( debug_here )then
                                call img%read(refs_glob,icls)
                                call img%write('rejected_pool_'//int2str(pool_iter)//'.mrc',cnt)
                            endif
                            img = 0.
                            call img%write(refs_glob,icls)
                        enddo
                        call debug_print('in reject from_pool 3')
                        call img%read(refs_glob, ncls_glob)
                        call img%write(refs_glob, ncls_glob)
                        deallocate(cls_mask)
                        write(logfhandle,'(A,I4,A,I6,A)')'>>> REJECTED FROM POOL: ',nptcls_rejected,' PARTICLES IN ',ncls_rejected,' CLUSTER(S)'
                    endif
                else
                    write(logfhandle,'(A,I4,A,I6,A)')'>>> NO PARTICLES FLAGGED FOR REJECTION FROM POOL'
                endif
                call img%kill
                call debug_print('end reject from_pool')
            end subroutine reject_from_pool

            !>  Convenience function for transfer_chunk_to_pool
            subroutine transfer_cavg( refs_in, dir, indin, refs_out, indout, self_transfer )
                character(len=*),  intent(in) :: refs_in, dir, refs_out
                integer,           intent(in) :: indin, indout
                logical, optional, intent(in) :: self_transfer
                type(image)                   :: img
                character(len=:), allocatable :: stkout, stkin
                integer :: iipart
                logical :: l_self
                call debug_print('in transfer_cavg '//int2str(indin)//' '//int2str(indout))
                l_self = .false.
                if( present(self_transfer) ) l_self = self_transfer
                if( l_self ) call debug_print('in transfer_cavg self_transfer')
                call img%new([box,box,1],smpd)
                call img%read( refs_in, indin)
                call img%write(refs_out,indout)
                stkin  = add2fbody(refs_in, params%ext,'_even')
                stkout = add2fbody(refs_out,params%ext,'_even')
                call img%read( stkin, indin)
                call img%write(stkout,indout)
                stkin  = add2fbody(refs_in,params%ext,'_odd')
                stkout = add2fbody(refs_out,params%ext,'_odd')
                call img%read( stkin, indin)
                call img%write(stkout,indout)
                ! temporary matrices, logics from chunk%read
                call img%new([boxpd,boxpd,1],smpd)
                call img%zero_and_flag_ft
                if( l_self )then
                    do iipart = 1,params%nparts
                        stkin = 'cavgs_even_part'//int2str_pad(iipart,params%numlen)//trim(params%ext)
                        call img%read(stkin, indin)
                        call img%write(stkin,indout)
                    enddo
                    do iipart = 1,params%nparts
                        stkin = 'cavgs_odd_part'//int2str_pad(iipart,params%numlen)//trim(params%ext)
                        call img%read(stkin, indin)
                        call img%write(stkin,indout)
                    enddo
                    do iipart = 1,params%nparts
                        stkin = 'ctfsqsums_even_part'//int2str_pad(iipart,params%numlen)//trim(params%ext)
                        call img%read(stkin, indin)
                        call img%write(stkin,indout)
                    enddo
                    do iipart = 1,params%nparts
                        stkin = 'ctfsqsums_odd_part'//int2str_pad(iipart,params%numlen)//trim(params%ext)
                        call img%read(stkin, indin)
                        call img%write(stkin,indout)
                    enddo
                else
                    stkin  = trim(dir)//'/cavgs_even_part'//trim(params%ext)
                    call img%read(stkin, indin)
                    do iipart = 1,params%nparts
                        stkout = 'cavgs_even_part'//int2str_pad(iipart,params%numlen)//trim(params%ext)
                        call img%write(stkout,indout)
                    enddo
                    stkin  = trim(dir)//'/cavgs_odd_part'//trim(params%ext)
                    call img%read(stkin, indin)
                    do iipart = 1,params%nparts
                        stkout = 'cavgs_odd_part'//int2str_pad(iipart,params%numlen)//trim(params%ext)
                        call img%write(stkout,indout)
                    enddo
                    stkin  = trim(dir)//'/ctfsqsums_even_part'//trim(params%ext)
                    call img%read(stkin, indin)
                    do iipart = 1,params%nparts
                        stkout = 'ctfsqsums_even_part'//int2str_pad(iipart,params%numlen)//trim(params%ext)
                        call img%write(stkout,indout)
                    enddo
                    stkin  = trim(dir)//'/ctfsqsums_odd_part'//trim(params%ext)
                    call img%read(stkin, indin)
                    do iipart = 1,params%nparts
                        stkout = 'ctfsqsums_odd_part'//int2str_pad(iipart,params%numlen)//trim(params%ext)
                        call img%write(stkout,indout)
                    enddo
                endif
                ! cleanup
                call img%kill
                call debug_print('end transfer_cavg')
            end subroutine transfer_cavg

            !> produces consolidated project at original scale
            subroutine write_snapshot( add_suffix )
                logical,           intent(in) :: add_suffix
                type(class_frcs)              :: frcs, frcs_sc
                type(oris)                    :: os_backup2, os_backup3
                character(len=:), allocatable :: projfile,projfname, cavgsfname, suffix, frcsfname, src, dest
                integer :: istk
                if( add_suffix )then
                    call debug_print('in write_snapshot suffix T')
                else
                    call debug_print('in write_snapshot suffix F')
                endif
                projfname  = get_fbody(orig_projfile, METADATA_EXT, separator=.false.)
                cavgsfname = get_fbody(refs_glob, params%ext, separator=.false.)
                frcsfname  = get_fbody(FRCS_FILE, BIN_EXT, separator=.false.)
                if( add_suffix )then
                    suffix     = '_'//trim(int2str(pool_proj%os_ptcl2D%get_noris()))//'nptcls'
                    cavgsfname = trim(cavgsfname)//trim(suffix)
                    frcsfname  = trim(frcsfname)//trim(suffix)
                endif
                call pool_proj%projinfo%set(1,'projname', projfname)
                projfile   = trim(projfname)//trim(METADATA_EXT)
                cavgsfname = trim(cavgsfname)//trim(params%ext)
                frcsfname  = trim(frcsfname)//trim(BIN_EXT)
                call pool_proj%projinfo%set(1,'projfile', projfile)
                call debug_print('in write_snapshot 1')
                if( add_suffix )then
                    if( trim(prev_snapshot_cavgs) /= '' )then
                        call del_file(prev_snapshot_frcs)
                        call del_file(prev_snapshot_cavgs)
                        src = add2fbody(prev_snapshot_cavgs, params%ext,'_even')
                        call del_file(src)
                        src = add2fbody(prev_snapshot_cavgs, params%ext,'_odd')
                        call del_file(src)
                    endif
                endif
                write(logfhandle,'(A,A,A,A)')'>>> GENERATING PROJECT SNAPSHOT ',trim(projfile), ' AT: ',cast_time_char(simple_gettime())
                if( do_autoscale )then
                    os_backup3 = pool_proj%os_cls2D
                    os_backup2 = pool_proj%os_stk
                    ! rescale classes
                    call rescale_cavgs(refs_glob, cavgsfname)
                    src  = add2fbody(refs_glob, params%ext,'_even')
                    dest = add2fbody(cavgsfname,params%ext,'_even')
                    call rescale_cavgs(src, dest)
                    src  = add2fbody(refs_glob, params%ext,'_odd')
                    dest = add2fbody(cavgsfname,params%ext,'_odd')
                    call rescale_cavgs(src, dest)
                    call pool_proj%os_out%kill
                    call pool_proj%add_cavgs2os_out(cavgsfname, orig_smpd, 'cavg')
                    pool_proj%os_cls2D = os_backup3
                    call os_backup3%kill
                    call debug_print('in write_snapshot 2')
                    ! rescale frcs
                    call frcs_sc%read(FRCS_FILE)
                    call frcs_sc%upsample(orig_smpd, orig_box, frcs)
                    call frcs%write(frcsfname)
                    call frcs%kill
                    call frcs_sc%kill
                    call pool_proj%add_frcs2os_out(frcsfname, 'frc2D')
                    call debug_print('in write_snapshot 3')
                    ! project updates to original scale
                    call pool_proj%os_stk%set_all2single('box', real(orig_box))
                    call pool_proj%os_stk%set_all2single('smpd',orig_smpd)
                    call pool_proj%os_ptcl2D%mul_shifts( 1./scale_factor )
                    do istk = 1,pool_proj%os_stk%get_noris()
                        call pool_proj%os_stk%set(istk,'stk',imported_stks(istk))
                    enddo
                    call debug_print('in write_snapshot 4')
                    ! write
                    pool_proj%os_ptcl3D = pool_proj%os_ptcl2D
                    call pool_proj%os_ptcl3D%delete_2Dclustering
                    call pool_proj%write(projfile)
                    call pool_proj%os_ptcl3D%kill
                    call debug_print('in write_snapshot 5')
                    ! preserve down-scaling
                    call pool_proj%os_ptcl2D%mul_shifts( scale_factor )
                    pool_proj%os_stk   = os_backup2
                    call os_backup2%kill
                    call debug_print('in write_snapshot 6')
                else
                    if( add_suffix )then
                        call simple_copy_file(FRCS_FILE, frcsfname)
                        call simple_copy_file(refs_glob, cavgsfname)
                        src  = add2fbody(refs_glob, params%ext,'_even')
                        dest = add2fbody(cavgsfname,params%ext,'_even')
                        call simple_copy_file(src, dest)
                        src  = add2fbody(refs_glob, params%ext,'_odd')
                        dest = add2fbody(cavgsfname,params%ext,'_odd')
                        call simple_copy_file(src, dest)
                    endif
                    call pool_proj%os_out%kill
                    call pool_proj%add_cavgs2os_out(cavgsfname, orig_smpd, 'cavg')
                    call pool_proj%add_frcs2os_out(frcsfname, 'frc2D')
                    ! write
                    pool_proj%os_ptcl3D = pool_proj%os_ptcl2D
                    call pool_proj%os_ptcl3D%delete_2Dclustering
                    call pool_proj%write(projfile)
                    call pool_proj%os_ptcl3D%kill
                endif
                ! we save the last snapshot name to simplify restarting
                call write_singlelineoftext(LAST_SNAPSHOT, projfile)
                ! cleanup previous snapshot
                prev_snapshot_frcs  = trim(frcsfname)
                prev_snapshot_cavgs = trim(cavgsfname)
                call debug_print('end write_snapshot')
            end subroutine write_snapshot

            !> for initial write of set of user adjustable parameters
            subroutine write_user_params
                type(oris) :: os
                call os%new(1, is_ptcl=.false.)
                call os%set(1,'lpthresh',params%lpthresh)
                call os%set(1,'ndev',    params%ndev)
                call os%write(USER_PARAMS)
                call os%kill
            end subroutine write_user_params

            !> updates current parameters with user input
            subroutine update_user_params
                type(oris) :: os
                real       :: lpthresh, ndev
                if( .not.file_exists(USER_PARAMS) ) return ! use of default/last update
                call debug_print('in update_user_params')
                lpthresh = params%lpthresh
                ndev     = params%ndev
                call os%new(1, is_ptcl=.false.)
                call os%read(USER_PARAMS)
                if( os%isthere(1,'lpthresh') )then
                    lpthresh = os%get(1,'lpthresh')
                    if( abs(lpthresh-params%lpthresh) > 0.001 )then
                        params%lpthresh = lpthresh
                        write(logfhandle,'(A,F8.2)')'>>> REJECTION LPTHRESH UPDATED TO: ',params%lpthresh
                    endif
                endif
                if( os%isthere(1,'ndev') )then
                    ndev = os%get(1,'ndev')
                    if( abs(ndev-params%ndev) > 0.001 )then
                        params%ndev = ndev
                        write(logfhandle,'(A,F8.2)')'>>> REJECTION NDEV    UPDATED TO: ',params%ndev
                    endif
                endif
                call os%kill
                call debug_print('end update_user_params')
            end subroutine update_user_params

            !> Sniffs for most recent snpashot in previous folder
            subroutine find_previous_snapshot()
                character(len=LONGSTRLEN), allocatable :: restart_snapshot(:)
                l_restart = .false.
                if( .not. cline%defined('dir_prev') ) return
                if( .not.file_exists(trim(params%dir_prev)//trim(LAST_SNAPSHOT)) ) then
                    THROW_HARD('Could not find previous project file for restart!')
                endif
                call read_filetable( trim(params%dir_prev)//trim(LAST_SNAPSHOT), restart_snapshot)
                if( size(restart_snapshot) /= 1 ) THROW_HARD('Invalid format for: '//trim(LAST_SNAPSHOT))
            end subroutine find_previous_snapshot

    end subroutine exec_cluster2D_stream

end module simple_commander_cluster2D_stream
