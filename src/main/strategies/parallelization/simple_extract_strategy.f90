! Strategy pattern for extract
!
! Hook-less, strategy-heavy, no separate common module.
!
module simple_extract_strategy
use simple_commanders_api
use simple_parameters, only: parameters
use simple_cmdline,    only: cmdline
use simple_qsys_env,   only: qsys_env
use simple_sp_project, only: sp_project
implicit none

public :: extract_strategy
public :: extract_inmem_strategy
public :: extract_distr_strategy
public :: create_extract_strategy
private

#include "simple_local_flags.inc"

! --------------------------------------------------------------------
! Strategy interface
! --------------------------------------------------------------------
type, abstract :: extract_strategy
contains
    procedure(apply_defaults_interface), deferred :: apply_defaults
    procedure(init_interface),           deferred :: initialize
    procedure(exec_interface),           deferred :: execute
    procedure(finalize_interface),       deferred :: finalize_run
    procedure(cleanup_interface),        deferred :: cleanup
    procedure(endmsg_interface),         deferred :: end_message
end type extract_strategy

! Shared-memory / distributed-worker strategy
type, extends(extract_strategy) :: extract_inmem_strategy
contains
    procedure :: apply_defaults => inmem_apply_defaults
    procedure :: initialize     => inmem_initialize
    procedure :: execute        => inmem_execute
    procedure :: finalize_run   => inmem_finalize_run
    procedure :: cleanup        => inmem_cleanup
    procedure :: end_message    => inmem_end_message
end type extract_inmem_strategy

! Distributed-master strategy
type, extends(extract_strategy) :: extract_distr_strategy
    type(sp_project)            :: spproj
    type(qsys_env)              :: qenv
    type(chash)                 :: job_descr
    type(string), allocatable   :: boxfiles(:)
    integer                     :: nmics_tot = 0
contains
    procedure :: apply_defaults => distr_apply_defaults
    procedure :: initialize     => distr_initialize
    procedure :: execute        => distr_execute
    procedure :: finalize_run   => distr_finalize_run
    procedure :: cleanup        => distr_cleanup
    procedure :: end_message    => distr_end_message
end type extract_distr_strategy

abstract interface
    subroutine apply_defaults_interface(self, cline)
        import :: extract_strategy, cmdline
        class(extract_strategy), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
    end subroutine apply_defaults_interface

    subroutine init_interface(self, params, cline)
        import :: extract_strategy, parameters, cmdline
        class(extract_strategy), intent(inout) :: self
        type(parameters),        intent(inout) :: params
        class(cmdline),          intent(inout) :: cline
    end subroutine init_interface

    subroutine exec_interface(self, params, cline)
        import :: extract_strategy, parameters, cmdline
        class(extract_strategy), intent(inout) :: self
        type(parameters),        intent(inout) :: params
        class(cmdline),          intent(inout) :: cline
    end subroutine exec_interface

    subroutine finalize_interface(self, params, cline)
        import :: extract_strategy, parameters, cmdline
        class(extract_strategy), intent(inout) :: self
        type(parameters),        intent(in)    :: params
        class(cmdline),          intent(inout) :: cline
    end subroutine finalize_interface

    subroutine cleanup_interface(self, params, cline)
        import :: extract_strategy, parameters, cmdline
        class(extract_strategy), intent(inout) :: self
        type(parameters),        intent(in)    :: params
        class(cmdline),          intent(inout) :: cline
    end subroutine cleanup_interface

    function endmsg_interface(self) result(msg)
        import :: extract_strategy
        class(extract_strategy), intent(in) :: self
        character(len=:), allocatable :: msg
    end function endmsg_interface
end interface

contains

    ! --------------------------------------------------------------------
    ! Factory
    ! --------------------------------------------------------------------
    function create_extract_strategy(cline) result(strategy)
        class(cmdline), intent(in) :: cline
        class(extract_strategy), allocatable :: strategy
        logical :: is_master
        ! Same heuristic as other refactors:
        ! distributed master when nparts defined but no explicit worker range/part.
        is_master = cline%defined('nparts') .and. (.not.cline%defined('part')) &
                   .and. (.not.cline%defined('fromp')) .and. (.not.cline%defined('top'))
        if( is_master )then
            allocate(extract_distr_strategy :: strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> DISTRIBUTED EXTRACT (MASTER)'
        else
            allocate(extract_inmem_strategy :: strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> EXTRACT (SHARED-MEMORY / WORKER)'
        endif
    end function create_extract_strategy

    ! --------------------------------------------------------------------
    ! Shared defaults (private; no separate common module)
    ! --------------------------------------------------------------------
    subroutine set_extract_defaults(cline)
        class(cmdline), intent(inout) :: cline
        if( .not. cline%defined('mkdir')          ) call cline%set('mkdir',          'yes')
        if( .not. cline%defined('outside')        ) call cline%set('outside',         'no')
        if( .not. cline%defined('pcontrast')      ) call cline%set('pcontrast',    'black')
        if( .not. cline%defined('stream')         ) call cline%set('stream',          'no')
        if( .not. cline%defined('extractfrommov') ) call cline%set('extractfrommov',  'no')
        if( .not. cline%defined('backgr_subtr')   ) call cline%set('backgr_subtr',    'no')
        if( cline%defined('ctf') )then
            if( cline%get_carg('ctf').ne.'flip' .and. cline%get_carg('ctf').ne.'no' )then
                THROW_HARD('Only CTF=NO/FLIP are allowed')
            endif
        endif
        call cline%set('oritype', 'mic')
    end subroutine set_extract_defaults

    ! Helper used by the distributed master for mapping mic->box from dir_box listing
    type(string) function boxfile_from_mic_impl(mic, boxfiles) result(boxfile)
        class(string), intent(in) :: mic
        type(string),  intent(in) :: boxfiles(:)
        type(string) :: box_from_mic
        integer      :: ibox
        boxfile = NIL
        if( size(boxfiles) == 0 ) return
        box_from_mic = fname_new_ext(basename(mic), 'box')
        do ibox = 1, size(boxfiles)
            if( basename(boxfiles(ibox)).eq.box_from_mic )then
                boxfile = boxfiles(ibox)
                return
            endif
        end do
    end function boxfile_from_mic_impl

    ! ====================================================================
    ! EXTRACT (SHARED-MEMORY / WORKER)
    ! ====================================================================
    subroutine inmem_apply_defaults(self, cline)
        class(extract_inmem_strategy), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        call set_extract_defaults(cline)
        call cline%set('oritype', 'mic')
        call cline%set('mkdir',   'no')  ! keep worker behavior (avoid nested dirs)
    end subroutine inmem_apply_defaults

    subroutine inmem_initialize(self, params, cline)
        class(extract_inmem_strategy), intent(inout) :: self
        type(parameters),              intent(inout) :: params
        class(cmdline),                intent(inout) :: cline
        call params%new(cline)
    end subroutine inmem_initialize

    subroutine inmem_execute(self, params, cline)
        use simple_ctf_estimate_fit,    only: ctf_estimate_fit
        use simple_particle_extractor,  only: ptcl_extractor
        class(extract_inmem_strategy), intent(inout) :: self
        type(parameters),              intent(inout) :: params
        class(cmdline),                intent(inout) :: cline
        type(image),                allocatable :: imgs(:)
        type(ptcl_extractor)                    :: extractor
        type(sp_project)                        :: spproj_in, spproj
        type(nrtxtfile)                         :: boxfile
        type(image)                             :: micrograph, micrograph_pad
        type(oris)                              :: os_mic
        type(ori)                               :: o_mic, o_tmp
        type(ctfparams)                         :: ctfparms
        type(ctf_estimate_fit)                  :: ctffit
        type(stack_io)                          :: stkio_w
        type(string)                            :: output_dir, mic_name, imgkind
        real,                       allocatable :: boxdata(:,:)
        integer,                    allocatable :: ptcl_inds(:)
        logical,                    allocatable :: oris_mask(:), mics_mask(:)
        type(string)                            :: stack, boxfile_name, box_fname, ctfdoc, ext
        real                                   :: ptcl_pos(2), stk_mean,stk_sdev,stk_max,stk_min,dfx,dfy,prog
        integer                                :: ldim(3), lfoo(3), fromto(2)
        integer                                :: nframes, imic, iptcl, nptcls,nmics,nmics_here,box, i, iptcl_g
        integer                                :: cnt, nmics_tot, ifoo, state, iptcl_glob, nptcls2extract
        logical                                :: l_ctfpatch, l_gid_present, l_ogid_present, prog_write, prog_part
        ! init
        output_dir = PATH_HERE
        fromto(:)  = [params%fromp, params%top]
        nmics_here = fromto(2)-fromto(1)+1
        prog_write = .false.
        prog_part  = .false.
        if( params%stream.eq.'yes' )then
            output_dir = DIR_EXTRACT
            if( cline%defined('dir') ) output_dir = params%dir//'/'
            call spproj%read(params%projfile)
            nmics_tot = spproj%os_mic%get_noris()
            fromto(:)  = [1,1]
            nmics_here = 1
            if( cline%defined('fromp') .and. cline%defined('top') )then
                fromto(:)  = [params%fromp, params%top]
                nmics_here = params%top-params%fromp+1
            endif
            if( nmics_tot /= nmics_here ) THROW_HARD('Incompatible # of integrated micrograph to process!')
        else
            call spproj_in%read_segment(params%oritype, params%projfile)
            nmics_tot = spproj_in%os_mic%get_noris()
            if( spproj_in%get_nintgs() == 0 ) THROW_HARD('No integrated micrograph to process!')
            call spproj%read_non_data_segments(params%projfile)
            call spproj%projinfo%set(1,'projname', get_fbody(params%outfile,METADATA_EXT,separator=.false.))
            call spproj%projinfo%set(1,'projfile', params%outfile)
            params%projfile = params%outfile
            call spproj%os_mic%new(nmics_here, is_ptcl=.false.)
            cnt = 0
            do imic = fromto(1),fromto(2)
                cnt = cnt + 1
                call spproj_in%os_mic%get_ori(imic, o_tmp)
                call spproj%os_mic%set_ori(cnt, o_tmp)
            enddo
            prog_write = .true.
            if( cline%defined('part') ) then
                prog_part = .true.
                ! call progressfile_init_part(cline%get_iarg('part'))
            else
                ! call progressfile_init()
            endif
            call spproj_in%kill
        endif
        ! input boxes
        if( cline%defined('dir_box') )then
            if( .not.file_exists(params%dir_box) )then
                write(logfhandle,*)'Directory does not exist: ', params%dir_box%to_char(), 'simple_commanders_pick::exec_extract'
                THROW_HARD('box directory does not exist; exec_extract')
            endif
        endif
        ! sanity checks
        allocate(mics_mask(1:nmics_here), source=.false.)
        nmics  = 0
        nptcls = 0
        do imic = 1,nmics_here
            call spproj%os_mic%get_ori(imic, o_mic)
            state = 1
            if( o_mic%isthere('state') ) state = o_mic%get_state()
            if( state == 0 ) cycle
            if( .not. o_mic%isthere('imgkind') )cycle
            if( .not. o_mic%isthere('intg')    )cycle
            call o_mic%getter('imgkind', imgkind)
            if( imgkind.ne.'mic') cycle
            call o_mic%getter('intg', mic_name)
            if( .not.file_exists(mic_name) )cycle
            ! box input
            if( cline%defined('dir_box') )then
                box_fname = params%dir_box//'/'//fname_new_ext(basename(mic_name),'box')
                if( .not.file_exists(box_fname) )cycle
                boxfile_name = simple_abspath(box_fname, check_exists=.false.)
                call spproj%set_boxfile(imic, boxfile_name)
            else
                boxfile_name = o_mic%get_str('boxfile')
                if( .not.file_exists(boxfile_name) )cycle
            endif
            call find_ldim_nptcls(mic_name, lfoo, nframes )
            if( nframes > 1 ) THROW_HARD('multi-frame extraction not supported; exec_extract')
            mics_mask(imic) = .true.
            nmics = nmics + 1
            if( nmics == 1 )call find_ldim_nptcls(mic_name, ldim, ifoo)
            if( nptcls == 0 .and. .not.cline%defined('box') )then
                if( nlines(boxfile_name) > 0 )then
                    call boxfile%new(boxfile_name, 1)
                    nptcls = boxfile%get_ndatalines()
                endif
                if( nptcls == 0 )then
                    call spproj%os_mic%set(imic, 'nptcls', 0)
                    cycle
                endif
                allocate( boxdata(boxfile%get_nrecs_per_line(),nptcls) )
                call boxfile%readNextDataLine(boxdata(:,1))
                call boxfile%kill
                params%box = nint(boxdata(3,1))
            endif
        enddo
        ! actual extraction
        if( nmics == 0 )then
            ! done
        else
            if( params%box == 0 )THROW_HARD('box cannot be zero; exec_extract')
            call micrograph%new([ldim(1),ldim(2),1], params%smpd)
            if( trim(params%extractfrommov).ne.'yes' ) call extractor%init_mic(params%box, (params%pcontrast .eq. 'black'))
            iptcl_glob = 0
            prog = 0.0
            do imic = 1,nmics_here
                if( .not.mics_mask(imic) )then
                    call spproj%os_mic%set(imic, 'nptcls', 0)
                    call spproj%os_mic%set_state(imic, 0)
                    cycle
                endif
                call spproj%os_mic%get_ori(imic, o_mic)
                boxfile_name = o_mic%get_str('boxfile')
                nptcls = 0
                if( nlines(boxfile_name) > 0 )then
                    call boxfile%new(boxfile_name, 1)
                    nptcls = boxfile%get_ndatalines()
                endif
                if( nptcls == 0 ) cycle
                call progress(imic,nmics_tot)
                if(allocated(oris_mask))deallocate(oris_mask)
                allocate(oris_mask(nptcls), source=.false.)
                if(allocated(boxdata))deallocate(boxdata)
                allocate( boxdata(boxfile%get_nrecs_per_line(),nptcls))
                do iptcl=1,nptcls
                    call boxfile%readNextDataLine(boxdata(:,iptcl))
                    box = nint(boxdata(3,iptcl))
                    if( nint(boxdata(3,iptcl)) /= nint(boxdata(4,iptcl)) )then
                        THROW_HARD('only square windows allowed; exec_extract')
                    endif
                    if( box /= params%box ) boxdata(1:2,iptcl) = boxdata(1:2,iptcl) - real(params%box-box)/2.
                    if( .not.cline%defined('box') .and. nint(boxdata(3,iptcl)) /= params%box )then
                        write(logfhandle,*) 'box_current: ', nint(boxdata(3,iptcl)), 'box in params: ', params%box
                        THROW_HARD('inconsistent box sizes in box files; exec_extract')
                    endif
                    oris_mask(iptcl) = (trim(params%outside).eq.'yes') .or. box_inside(ldim, nint(boxdata(1:2,iptcl)), params%box)
                end do
                nptcls2extract = count(oris_mask)
                call spproj%os_mic%set(imic, 'nptcls', nptcls2extract)

                if( nptcls2extract == 0 )then
                    mics_mask(imic) = .false.
                    cycle
                endif
                if(allocated(ptcl_inds))deallocate(ptcl_inds)
                allocate(ptcl_inds(nptcls2extract), source=0)
                cnt = 0
                do iptcl=1,nptcls
                    if( oris_mask(iptcl) )then
                        cnt = cnt + 1
                        ptcl_inds(cnt) = iptcl
                    endif
                enddo
                ctfparms      = o_mic%get_ctfvars()
                ctfparms%smpd = params%smpd
                if( o_mic%isthere('dfx') )then
                    if( .not.o_mic%isthere('cs') .or. .not.o_mic%isthere('kv') .or. .not.o_mic%isthere('fraca') )then
                        THROW_HARD('input lacks at least cs, kv or fraca; exec_extract')
                    endif
                endif
                call o_mic%getter('intg', mic_name)
                ext   = fname2ext(basename(mic_name))
                stack = output_dir//EXTRACT_STK_FBODY//get_fbody(basename(mic_name), ext)//STK_EXT
                call prepimgbatch(nptcls2extract)
                if( trim(params%extractfrommov).eq.'yes' )then
                    if( trim(params%ctf).eq.'flip' .and. o_mic%isthere('dfx') )then
                        THROW_HARD('extractfrommov=yes does not support ctf=flip yet')
                    endif
                    call extractor%init_mov(o_mic, params%box, (params%pcontrast .eq. 'black'))
                    call extractor%extract_particles(ptcl_inds, nint(boxdata), imgs, stk_min,stk_max,stk_mean,stk_sdev)
                else
                    call micrograph%read(mic_name, 1)
                    if( trim(params%backgr_subtr).eq.'yes') call micrograph%subtract_background(HP_BACKGR_SUBTR)

                    call extractor%extract_particles_from_mic(micrograph, ptcl_inds, nint(boxdata), imgs,&
                        &stk_min,stk_max,stk_mean,stk_sdev)
                endif
                call stkio_w%open(stack, params%smpd, 'write', box=params%box)
                do i = 1,nptcls2extract
                    call stkio_w%write(i, imgs(i))
                enddo
                call stkio_w%close
                call imgs(1)%update_header_stats(stack, [stk_min, stk_max, stk_mean, stk_sdev])
                call spproj%add_stk(stack, ctfparms)
                l_ctfpatch = .false.
                if( o_mic%isthere('ctfdoc') )then
                    ctfdoc = o_mic%get_str('ctfdoc')
                    if( file_exists(ctfdoc) )then
                        call ctffit%read_doc(ctfdoc)
                        l_ctfpatch = .true.
                    endif
                endif
                l_ogid_present = o_mic%isthere('ogid')
                l_gid_present  = o_mic%isthere('gid')
                !$omp parallel do schedule(static) default(shared) proc_bind(close)&
                !$omp private(i,iptcl,iptcl_g,ptcl_pos,dfx,dfy)
                do i = 1,nptcls2extract
                    iptcl    = ptcl_inds(i)
                    iptcl_g  = iptcl_glob + i
                    ptcl_pos = boxdata(1:2,iptcl)
                    call spproj%set_boxcoords(iptcl_g, nint(ptcl_pos))
                    if( l_ctfpatch )then
                        ptcl_pos = ptcl_pos+1.+real(params%box/2)
                        call ctffit%pix2polyvals(ptcl_pos(1),ptcl_pos(2), dfx,dfy)
                        call spproj%os_ptcl2D%set_dfx(iptcl_g,dfx)
                        call spproj%os_ptcl2D%set_dfy(iptcl_g,dfy)
                        call spproj%os_ptcl3D%set_dfx(iptcl_g,dfx)
                        call spproj%os_ptcl3D%set_dfy(iptcl_g,dfy)
                    endif

                    if( l_ogid_present )then
                        call spproj%os_ptcl2D%set(iptcl_g,'ogid',o_mic%get('ogid'))
                        call spproj%os_ptcl3D%set(iptcl_g,'ogid',o_mic%get('ogid'))
                    endif

                    if( l_gid_present )then
                        call spproj%os_ptcl2D%set(iptcl_g,'gid',o_mic%get('gid'))
                        call spproj%os_ptcl3D%set(iptcl_g,'gid',o_mic%get('gid'))
                    endif
                end do
                !$omp end parallel do
                iptcl_glob = iptcl_glob + nptcls2extract
                call boxfile%kill
                call ctffit%kill
                if(prog_write) then
                    if( (real(imic) / real(nmics_here)) > prog + 0.05 ) then
                        prog = real(imic) / real(nmics_here)
                        if(prog_part) then
                            ! call progressfile_update_part(cline%get_iarg('part'), prog)
                        else
                            ! call progressfile_update(prog)
                        endif
                        write(logfhandle,'(f4.0,1x,a)') 100.*prog, 'percent of the micrographs processed'
                    endif
                endif
            enddo
            call killimgbatch
        endif
        if( trim(params%stream).eq.'yes' )then
            nmics = count(mics_mask)
            if( nmics == 0 )then
                call spproj%os_mic%kill
                call spproj%os_stk%kill
                call spproj%os_ptcl2D%kill
                call spproj%os_ptcl3D%kill
            else
                if( nmics < nmics_here )then
                    call os_mic%new(nmics, is_ptcl=.false.)
                    cnt = 0
                    do imic = 1, nmics_here
                        if( mics_mask(imic) )then
                            cnt = cnt+1
                            call os_mic%transfer_ori(cnt, spproj%os_mic, imic)
                        endif
                    enddo
                    spproj%os_mic = os_mic
                endif
            endif
        endif
        call spproj%write(params%projfile)
        ! cleanup resources (kept here; finalize_run only reports job finished)
        call extractor%kill
        call micrograph%kill
        call micrograph_pad%kill
        call o_mic%kill
        call o_tmp%kill
        call os_mic%kill
        call spproj%kill

        contains

            subroutine prepimgbatch( batchsz )
                integer,           intent(in) :: batchsz
                integer :: i
                logical :: doprep
                doprep = .false.
                if( .not. allocated(imgs) )then
                    doprep = .true.
                else
                    if( batchsz > size(imgs) ) doprep = .true.
                    if( doprep ) call killimgbatch
                endif
                if( doprep )then
                    allocate(imgs(batchsz))
                    !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
                    do i = 1,batchsz
                        call imgs(i)%new([params%box, params%box, 1], params%smpd, wthreads=.false.)
                    end do
                    !$omp end parallel do
                endif
            end subroutine prepimgbatch

            subroutine killimgbatch
                integer :: i
                if( allocated(imgs) )then
                    do i = 1,size(imgs)
                        call imgs(i)%kill
                    end do
                    deallocate(imgs)
                endif
            end subroutine killimgbatch

            logical function box_inside(ildim, coord, box)
                integer, intent(in) :: ildim(3), coord(2), box
                integer             :: fromc(2), toc(2)
                fromc = coord + 1
                toc   = fromc + (box - 1)
                box_inside = .true.
                if( any(fromc < 1) .or. toc(1) > ildim(1) .or. toc(2) > ildim(2) ) box_inside = .false.
            end function box_inside

    end subroutine inmem_execute

    subroutine inmem_finalize_run(self, params, cline)
        class(extract_inmem_strategy), intent(inout) :: self
        type(parameters),              intent(in)    :: params
        class(cmdline),                intent(inout) :: cline
        call qsys_job_finished(params, string('simple_commanders_pick :: exec_extract'))
    end subroutine inmem_finalize_run

    subroutine inmem_cleanup(self, params, cline)
        class(extract_inmem_strategy), intent(inout) :: self
        type(parameters),              intent(in)    :: params
        class(cmdline),                intent(inout) :: cline
        ! No-op
    end subroutine inmem_cleanup

    function inmem_end_message(self) result(msg)
        class(extract_inmem_strategy), intent(in) :: self
        character(len=:), allocatable :: msg
        msg = '**** SIMPLE_EXTRACT NORMAL STOP ****'
    end function inmem_end_message

    ! ====================================================================
    ! DISTRIBUTED EXTRACT (MASTER)
    ! ====================================================================
    subroutine distr_apply_defaults(self, cline)
        class(extract_distr_strategy), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        call set_extract_defaults(cline)
        call cline%set('oritype', 'mic')
    end subroutine distr_apply_defaults

    subroutine distr_initialize(self, params, cline)
        class(extract_distr_strategy), intent(inout) :: self
        type(parameters),              intent(inout) :: params
        class(cmdline),                intent(inout) :: cline
        type(ori)      :: o_mic
        type(string)   :: mic_name, imgkind, boxfile_name
        integer        :: lfoo(3)
        integer        :: nframes, imic, nmics_tot, nmics, state
        type(string)   :: boxfile_from_dir
        call params%new(cline)
        ! mirror master behavior: avoid nested directory structure in workers
        call cline%set('mkdir', 'no')
        call self%spproj%read(params%projfile)
        if( self%spproj%get_nintgs() == 0 ) THROW_HARD('No integrated micrograph to process!')
        nmics_tot = self%spproj%os_mic%get_noris()
        self%nmics_tot = nmics_tot
        if( nmics_tot < params%nparts ) params%nparts = nmics_tot
        ! wipe previous stacks & particles
        call self%spproj%os_stk%kill
        call self%spproj%os_ptcl2D%kill
        call self%spproj%os_ptcl3D%kill
        call self%spproj%os_cls2D%kill
        call self%spproj%os_cls3D%kill
        call self%spproj%os_out%kill
        ! input directory for boxfiles (optional)
        if( cline%defined('dir_box') )then
            if( trim(params%mkdir).eq.'yes' .and. params%dir_box%to_char([1,1]).ne.'/')then
                params%dir_box = filepath(PATH_PARENT,params%dir_box)
            endif
            params%dir_box = simple_abspath(params%dir_box)
            if( file_exists(params%dir_box) )then
                call simple_list_files_regexp(params%dir_box,'\.box$', self%boxfiles)
                if(.not.allocated(self%boxfiles))then
                    write(logfhandle,*)'No box file found in ', params%dir_box%to_char(), '; exec_extract_distr'
                    THROW_HARD('No box file found; exec_extract_distr, 1')
                endif
                if(size(self%boxfiles)==0)then
                    write(logfhandle,*)'No box file found in ', params%dir_box%to_char(), '; exec_extract_distr'
                    THROW_HARD('No box file found; exec_extract_distr 2')
                endif
            else
                write(logfhandle,*)'Directory does not exist: ', params%dir_box%to_char(), 'exec_extract_distr'
                THROW_HARD('box directory does not exist; exec_extract_distr')
            endif
            call cline%set('dir_box', params%dir_box)
        endif
        call self%spproj%write(params%projfile)
        ! sanity checks (count eligible micrographs)
        nmics  = 0
        do imic = 1, nmics_tot
            call self%spproj%os_mic%get_ori(imic, o_mic)
            state = 1
            if( o_mic%isthere('state') ) state = o_mic%get_state()
            if( state == 0 ) cycle
            if( .not. o_mic%isthere('imgkind') )cycle
            if( .not. o_mic%isthere('intg')    )cycle
            call o_mic%getter('imgkind', imgkind)
            if( imgkind.ne.'mic') cycle
            call o_mic%getter('intg', mic_name)
            if( .not.file_exists(mic_name) )cycle
            ! box input
            if( cline%defined('dir_box') )then
                boxfile_from_dir = boxfile_from_mic_impl(mic_name, self%boxfiles)
                if( boxfile_from_dir.eq.NIL )cycle
            else
                call o_mic%getter('boxfile', boxfile_name)
                if( .not.file_exists(boxfile_name) )cycle
            endif
            call find_ldim_nptcls(mic_name, lfoo, nframes )
            if( nframes > 1 ) THROW_HARD('multi-frame extraction not supported; exec_extract_distr')
            nmics = nmics + 1
        enddo
        if( nmics == 0 ) THROW_HARD('No particles to extract! exec_extract_distr')
        call o_mic%kill
        ! setup distributed environment + job description
        call self%qenv%new(params, params%nparts)
        call cline%gen_job_descr(self%job_descr)
    end subroutine distr_initialize

    subroutine distr_execute(self, params, cline)
        class(extract_distr_strategy), intent(inout) :: self
        type(parameters),              intent(inout) :: params
        class(cmdline),                intent(inout) :: cline
        type(sp_project)           :: spproj_part
        type(ori)                  :: o_mic, o_tmp
        type(oris)                 :: os_stk
        type(string), allocatable  :: stktab(:), parts_fname(:)
        integer                    :: boxcoords(2)
        type(string)               :: partsfile
        real                       :: dfx, dfy, ogid, gid
        integer                    :: imic, i, nmics_tot, cnt, istk, nstks, ipart, numlen
        ! schedule jobs
        call self%qenv%gen_scripts_and_schedule_jobs( self%job_descr, algnfbody=string(ALGN_FBODY), &
            &array=L_USE_SLURM_ARR, extra_params=params)
        ! ASSEMBLY
        allocate(parts_fname(params%nparts))
        numlen = len(int2str(params%nparts))
        do ipart = 1, params%nparts
            parts_fname(ipart) = ALGN_FBODY//int2str_pad(ipart,numlen)//METADATA_EXT
        enddo
        nmics_tot = self%nmics_tot
        ! copy updated micrographs
        cnt   = 0
        nstks = 0
        do ipart = 1, params%nparts
            partsfile = parts_fname(ipart)
            call spproj_part%read_segment('mic', partsfile)
            do imic = 1, spproj_part%os_mic%get_noris()
                cnt = cnt + 1
                call spproj_part%os_mic%get_ori(imic, o_mic)
                call self%spproj%os_mic%set_ori(cnt, o_mic)
                if( o_mic%isthere('nptcls') )then
                    if( o_mic%get_int('nptcls') > 0 ) nstks = nstks + 1
                endif
            enddo
            call spproj_part%kill
        enddo
        if( cnt /= nmics_tot ) THROW_HARD('Inconstistent number of micrographs in individual projects')
        ! fetch stacks table
        if( nstks > 0 )then
            call os_stk%new(nstks, is_ptcl=.false.)
            allocate(stktab(nstks))
            cnt = 0
            do ipart = 1, params%nparts
                call spproj_part%read_segment('stk', parts_fname(ipart))
                do istk = 1, spproj_part%os_stk%get_noris()
                    cnt = cnt + 1
                    call spproj_part%os_stk%get_ori(istk, o_tmp)
                    call os_stk%set_ori(cnt, o_tmp)
                    stktab(cnt) = os_stk%get_str(cnt,'stk')
                enddo
                call spproj_part%kill
            enddo
            call self%spproj%add_stktab(stktab, os_stk)
            ! transfer particle locations + defocus + group IDs
            cnt = 0
            do ipart = 1, params%nparts
                call spproj_part%read_segment('ptcl2D', parts_fname(ipart))
                do i = 1, spproj_part%os_ptcl2D%get_noris()
                    cnt = cnt + 1
                    call spproj_part%get_boxcoords(i, boxcoords)
                    call self%spproj%set_boxcoords(cnt, boxcoords)
                    if( spproj_part%os_ptcl2D%isthere(i,'dfx') )then
                        dfx = spproj_part%os_ptcl2D%get_dfx(i)
                        dfy = spproj_part%os_ptcl2D%get_dfy(i)
                        call self%spproj%os_ptcl2D%set_dfx(cnt, dfx)
                        call self%spproj%os_ptcl2D%set_dfy(cnt, dfy)
                        call self%spproj%os_ptcl3D%set_dfx(cnt, dfx)
                        call self%spproj%os_ptcl3D%set_dfy(cnt, dfy)
                    endif
                    if( spproj_part%os_ptcl2D%isthere(i,'ogid') )then
                        ogid = spproj_part%os_ptcl2D%get(i, 'ogid')
                        call self%spproj%os_ptcl2D%set(cnt,'ogid',ogid)
                        call self%spproj%os_ptcl3D%set(cnt,'ogid',ogid)
                    endif
                    if( spproj_part%os_ptcl2D%isthere(i,'gid') )then
                        gid = spproj_part%os_ptcl2D%get(i, 'gid')
                        call self%spproj%os_ptcl2D%set(cnt,'gid',gid)
                        call self%spproj%os_ptcl3D%set(cnt,'gid',gid)
                    endif
                enddo
                call spproj_part%kill
            enddo
            call os_stk%kill
            if( allocated(stktab) ) deallocate(stktab)
        endif
        ! final write
        call self%spproj%write(params%projfile)
        ! cleanup locals
        if( allocated(parts_fname) ) deallocate(parts_fname)
        call o_mic%kill
        call o_tmp%kill
    end subroutine distr_execute

    subroutine distr_finalize_run(self, params, cline)
        class(extract_distr_strategy), intent(inout) :: self
        type(parameters),              intent(in)    :: params
        class(cmdline),                intent(inout) :: cline
        ! No-op (master ends via workflow)
    end subroutine distr_finalize_run

    subroutine distr_cleanup(self, params, cline)
        class(extract_distr_strategy), intent(inout) :: self
        type(parameters),              intent(in)    :: params
        class(cmdline),                intent(inout) :: cline
        call self%spproj%kill
        if( allocated(self%boxfiles) ) deallocate(self%boxfiles)
        call qsys_cleanup(params)
        call self%qenv%kill
        call self%job_descr%kill
    end subroutine distr_cleanup

    function distr_end_message(self) result(msg)
        class(extract_distr_strategy), intent(in) :: self
        character(len=:), allocatable :: msg
        msg = '**** SIMPLE_EXTRACT_DISTR NORMAL STOP ****'
    end function distr_end_message

end module simple_extract_strategy