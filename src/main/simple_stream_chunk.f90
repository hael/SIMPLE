module simple_stream_chunk
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
implicit none

public :: stream_chunk

private
#include "simple_local_flags.inc"

character(len=STDLEN), parameter   :: PROJNAME_CHUNK      = 'chunk'
logical,               parameter   :: DEBUG_HERE          = .false.

type stream_chunk
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

    subroutine init( self, id, master_spproj )
        class(stream_chunk), intent(inout) :: self
        integer,             intent(in)    :: id
        class(sp_project),   intent(in)    :: master_spproj
        call debug_print('in chunk%init '//int2str(id))
        call self%spproj%kill
        self%id        = id
        self%it        = 0
        self%nmics     = 0
        self%nptcls    = 0
        self%path      = './chunk_'//int2str(id)//'/'
        self%projfile_out = ''
        self%spproj%projinfo = master_spproj%projinfo
        self%spproj%compenv  = master_spproj%compenv
        if( params_glob%nparts_chunk == 1 )then
            ! shared memory
            call self%qenv%new(params_glob%nparts_chunk, exec_bin='simple_private_exec')
            call self%spproj%compenv%set(1,'qsys_name','local')
        else
            ! we need to override the qsys_name for non local distributed execution
            call self%qenv%new(params_glob%nparts_chunk, exec_bin='simple_private_exec', qsys_name='local')
        endif
        call self%spproj%projinfo%delete_entry('projname')
        call self%spproj%projinfo%delete_entry('projfile')
        if( allocated(self%orig_stks) ) deallocate(self%orig_stks)
        self%converged = .false.
        self%available = .true.
        call debug_print('end chunk%init '//int2str(id))
    end subroutine init

    subroutine generate( self, fnames, nptcls, ind_glob )
        class(stream_chunk),       intent(inout) :: self
        character(len=LONGSTRLEN), intent(in)    :: fnames(:)
        integer,                   intent(in)    :: nptcls, ind_glob
        type(sp_project)              :: spproj
        character(len=:), allocatable :: stack_name
        logical, allocatable          :: spproj_mask(:)
        integer :: iproj, iptcl, cnt, inptcls, nptcls_tot, fromp, iiproj, n_in, inmics, instks
        if( .not.self%available ) THROW_HARD('chunk unavailable; chunk%generate')
        n_in = size(fnames)
        if( n_in == 0 ) THROW_HARD('# ptcls == 0; chunk%generate')
        call debug_print('in chunk%generate '//int2str(self%id)//' '//int2str(ind_glob)//' '//int2str(nptcls)//' '//int2str(n_in))
        allocate(spproj_mask(n_in),source=.false.)
        nptcls_tot = 0
        do iproj = 1,n_in
            call spproj%read_data_info(fnames(iproj), inmics, instks, inptcls)
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
        cnt    = 0
        fromp  = 1
        iiproj = 0
        do iproj = 1,n_in
            if( .not.spproj_mask(iproj) ) cycle
            iiproj = iiproj + 1
            call spproj%read_mic_stk_ptcl2D_segments(fnames(iproj))
            inptcls = nint(spproj%os_mic%get(1,'nptcls'))
            call self%spproj%os_mic%transfer_ori(iiproj, spproj%os_mic, 1)
            call self%spproj%os_stk%transfer_ori(iiproj, spproj%os_stk, 1)
            ! update to stored stack file name because we are in the root folder
            stack_name = trim(spproj%get_stkname(1))
            if( .not.file_exists(stack_name) )then
                ! for cluster2D_stream, 4 is for '../'
                self%orig_stks(iiproj) = simple_abspath(trim(stack_name(4:)))
            else
                ! for cluster2d_subsets
                self%orig_stks(iiproj) = simple_abspath(trim(stack_name))
            endif
            call self%spproj%os_stk%set(iiproj, 'stk', self%orig_stks(iiproj))
            do iptcl = 1,inptcls
                cnt = cnt + 1
                call self%spproj%os_ptcl2D%transfer_ori(cnt, spproj%os_ptcl2D, iptcl)
                call self%spproj%os_ptcl2D%set(cnt, 'stkind', real(iiproj))
            enddo
            call self%spproj%os_stk%set(iiproj, 'fromp', real(fromp))
            call self%spproj%os_stk%set(iiproj, 'top',   real(fromp+inptcls-1))
            fromp = fromp + inptcls
        enddo
        call spproj%kill
        self%spproj%os_ptcl3D = self%spproj%os_ptcl2D
        call debug_print('end chunk%generate '//int2str(self%id))
    end subroutine generate

    subroutine exec_classify( self, cline_classify, orig_smpd, orig_box, box )
        class(stream_chunk), intent(inout) :: self
        class(cmdline),      intent(inout) :: cline_classify
        real,                intent(in)    :: orig_smpd
        integer,             intent(in)    :: orig_box, box
        character(len=XLONGSTRLEN) :: cwd
        real    :: scale
        integer :: nptcls_sel
        call debug_print('in chunk%exec_classify '//int2str(self%id))
        if( .not.self%available ) return
        if( self%nptcls == 0 ) return
        call simple_mkdir(self%path)
        call chdir(self%path)
        call simple_getcwd(cwd)
        cwd_glob    = trim(cwd)
        self%projfile_out = trim(PROJNAME_CHUNK)//trim(METADATA_EXT)
        call simple_mkdir(STDERROUT_DIR)
        nptcls_sel = self%spproj%os_ptcl2D%get_noris(consider_state=.true.)
        if( box < orig_box )then
            scale = real(box) / real(orig_box)
            call cline_classify%set('smpd',      orig_smpd)
            call cline_classify%set('box',       real(orig_box))
            call cline_classify%set('smpd_crop', orig_smpd / scale)
            call cline_classify%set('box_crop',  real(box))
        endif
        call cline_classify%set('projfile', self%projfile_out)
        call cline_classify%set('projname', trim(PROJNAME_CHUNK))
        call self%spproj%update_projinfo(cline_classify)
        call self%spproj%write()
        call self%qenv%exec_simple_prg_in_queue_async(cline_classify, './distr_chunk2D', 'simple_log_chunk2d')
        call chdir('..')
        call simple_getcwd(cwd_glob)
        ! cleanup
        call self%spproj%kill
        call self%qenv%kill
        ! chunk is now busy
        self%available = .false.
        self%converged = .false.
        write(logfhandle,'(A,I6,A,I6,A)')'>>> CHUNK ',self%id,' INITIATED CLASSIFICATION WITH ',nptcls_sel,' PARTICLES'
        call debug_print('end chunk%exec_classify')
    end subroutine exec_classify

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
                integer                    :: icls, ipart, numlen_chunk, iostat
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

    subroutine remove_folder( self )
        class(stream_chunk), intent(inout) :: self
        call debug_print('in chunk%remove_folder '//int2str(self%id))
        if( .not.self%converged )THROW_HARD('cannot remove chunk prior to convergence; remove_folder')
        call simple_rmdir(self%path)
        call debug_print('end chunk%remove_folder '//int2str(self%id))
    end subroutine remove_folder

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

    logical function has_converged( self )
        class(stream_chunk), intent(inout) :: self
        self%converged = file_exists(trim(self%path)//trim(CLUSTER2D_FINISHED))
        has_converged  = self%converged
    end function has_converged

    subroutine reject( self, res_thresh, ndev, box )
        class(stream_chunk), intent(inout) :: self
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
        self%nmics     = 0
        self%nptcls    = 0
        self%path      = ''
        self%projfile_out = ''
        if( allocated(self%orig_stks) ) deallocate(self%orig_stks)
        self%converged = .false.
        self%available = .false.
    end subroutine kill

end module simple_stream_chunk
