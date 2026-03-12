!@descr: for picking, extraction, and making picking references
module simple_commanders_pick
use simple_commanders_api
! use simple_progress
implicit none

public :: commander_pick_distr
public :: commander_pick
public :: commander_extract_distr
public :: commander_extract
public :: commander_reextract_distr
public :: commander_reextract
public :: commander_pick_extract
public :: commander_make_pickrefs
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_pick_distr
  contains
    procedure :: execute      => exec_pick_distr
end type commander_pick_distr

type, extends(commander_base) :: commander_pick
  contains
    procedure :: execute      => exec_pick
end type commander_pick

type, extends(commander_base) :: commander_extract_distr
  contains
    procedure :: execute      => exec_extract_distr
end type commander_extract_distr

type, extends(commander_base) :: commander_extract
  contains
    procedure :: execute      => exec_extract
end type commander_extract

type, extends(commander_base) :: commander_reextract_distr
  contains
    procedure :: execute      => exec_reextract_distr
end type commander_reextract_distr

type, extends(commander_base) :: commander_reextract
  contains
    procedure :: execute      => exec_reextract
end type commander_reextract

type, extends(commander_base) :: commander_pick_extract
  contains
    procedure :: execute      => exec_pick_extract
end type commander_pick_extract

type, extends(commander_base) :: commander_make_pickrefs
  contains
    procedure :: execute      => exec_make_pickrefs
end type commander_make_pickrefs

contains

    subroutine exec_pick_distr( self, cline )
        class(commander_pick_distr), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        call run_pick_workflow(cline)
    end subroutine exec_pick_distr

    subroutine exec_pick( self, cline )
        class(commander_pick), intent(inout) :: self
        class(cmdline),        intent(inout) :: cline
        call run_pick_workflow(cline)
    end subroutine exec_pick

    subroutine run_pick_workflow( cline )
        use simple_core_module_api, only: simple_end
        use simple_pick_strategy,   only: pick_strategy, create_pick_strategy
        use simple_cmdline,         only: cmdline
        use simple_parameters,      only: parameters
        class(cmdline), intent(inout) :: cline
        class(pick_strategy), allocatable :: strategy
        type(parameters) :: params
        call cline%set('prg', 'pick')
        strategy = create_pick_strategy(cline)
        call strategy%initialize(params, cline)
        call strategy%execute(params, cline)
        call strategy%finalize_run(params, cline)
        call strategy%cleanup(params, cline)
        call simple_end(strategy%end_message())
        if( allocated(strategy) ) deallocate(strategy)
    end subroutine run_pick_workflow

    ! Unified extract workflow (worker/shared-memory + distributed master)
    ! driven by runtime polymorphism (strategy pattern).

    subroutine exec_extract_distr( self, cline )
        class(commander_extract_distr), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        call run_extract_workflow(cline)
    end subroutine exec_extract_distr

    subroutine exec_extract( self, cline )
        class(commander_extract), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        call run_extract_workflow(cline)
    end subroutine exec_extract

    subroutine run_extract_workflow( cline )
        use simple_core_module_api,   only: simple_end
        use simple_extract_strategy,  only: extract_strategy, create_extract_strategy
        use simple_cmdline,           only: cmdline
        use simple_parameters,        only: parameters
        class(cmdline), intent(inout) :: cline
        class(extract_strategy), allocatable :: strategy
        type(parameters) :: params
        ! Helps distributed job script generation if it relies on 'prg'
        call cline%set('prg', 'extract')
        strategy = create_extract_strategy(cline)
        call strategy%apply_defaults(cline)
        call strategy%initialize(params, cline)
        call strategy%execute(params, cline)
        call strategy%finalize_run(params, cline)
        call strategy%cleanup(params, cline)
        call simple_end(strategy%end_message())
        if( allocated(strategy) ) deallocate(strategy)
    end subroutine run_extract_workflow

    subroutine exec_reextract_distr( self, cline )
        class(commander_reextract_distr), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline !< command line input
        type(parameters)                        :: params
        type(sp_project)                        :: spproj
        type(sp_project),           allocatable :: spproj_parts(:)
        type(qsys_env)                          :: qenv
        type(chash)                             :: job_descr
        type(ori)                               :: o_mic, o
        type(oris)                              :: os_stk
        type(chash),                allocatable :: part_params(:)
        type(string),               allocatable :: boxfiles(:), stktab(:), parts_fname(:)
        type(string)                            :: mic_name, imgkind
        integer,                    allocatable :: parts(:,:)
        integer :: imic,i,nmics_tot,numlen,nmics,cnt,state,istk,nstks,ipart,stkind,nptcls
        if( cline%defined('ctf') )then
            if( cline%get_carg('ctf').ne.'flip' .and. cline%get_carg('ctf').ne.'no' )then
                THROW_HARD('Only CTF=NO/FLIP are allowed')
            endif
        endif
        if( cline%defined('osmpd') )then
            if( .not.cline%defined('box') ) THROW_HARD('BOX must be defined with OSMPD!')
        endif
        if( .not. cline%defined('mkdir')     )     call cline%set('mkdir',          'yes')
        if( .not. cline%defined('pcontrast') )     call cline%set('pcontrast',    'black')
        if( .not. cline%defined('oritype')   )     call cline%set('oritype',     'ptcl3D')
        if( .not. cline%defined('extractfrommov')) call cline%set('extractfrommov',  'no')
        if( .not. cline%defined('backgr_subtr'))   call cline%set('backgr_subtr',    'no')
        call params%new(cline)
        call cline%set('mkdir', 'no')
        ! read in integrated movies
        call spproj%read( params%projfile )
        if( spproj%get_nintgs() == 0 ) THROW_HARD('No integrated micrograph to process!')
        if( spproj%get_nstks() == 0 ) THROW_HARD('This project file does not contain stacks!')
        nmics_tot = spproj%os_mic%get_noris()
        if( nmics_tot < params%nparts )then
            params%nparts = nmics_tot
        endif
        ! sanity checks
        nmics  = 0
        do imic = 1, nmics_tot
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
            ! update counter
            nmics = nmics + 1
            ! removes boxfile from micrographs
            call spproj%os_mic%delete_entry(imic,'boxfile')
        enddo
        if( nmics == 0 )then
            THROW_WARN('No particles to re-extract! exec_reextract')
            return
        endif
        call spproj%os_mic%kill
        call spproj%os_stk%kill
        call spproj%os_ptcl2D%kill
        call spproj%os_ptcl3D%kill
        ! DISTRIBUTED EXTRACTION
        ! setup the environment for distributed execution
        parts = split_nobjs_even(nmics_tot, params%nparts)
        allocate(part_params(params%nparts))
        do ipart=1,params%nparts
            call part_params(ipart)%new(2)
            call part_params(ipart)%set('fromp',int2str(parts(ipart,1)))
            call part_params(ipart)%set('top',  int2str(parts(ipart,2)))
        end do
        call qenv%new(params, params%nparts)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! schedule & clean
        call qenv%gen_scripts_and_schedule_jobs( job_descr, algnfbody=string(ALGN_FBODY), part_params=part_params, array=L_USE_SLURM_ARR, extra_params=params)
        ! ASSEMBLY
        allocate(spproj_parts(params%nparts),parts_fname(params%nparts))
        numlen = len(int2str(params%nparts))
        do ipart = 1,params%nparts
            parts_fname(ipart) = ALGN_FBODY//int2str_pad(ipart,numlen)//METADATA_EXT
        enddo
        ! copy updated micrographs
        cnt   = 0
        nmics = 0
        do ipart = 1,params%nparts
            call spproj_parts(ipart)%read_segment('mic',parts_fname(ipart))
            nmics = nmics + spproj_parts(ipart)%os_mic%get_noris()
        enddo
        if( nmics > 0 )then
            call spproj%os_mic%new(nmics, is_ptcl=.false.)
            ! transfer stacks
            cnt   = 0
            nstks = 0
            do ipart = 1,params%nparts
                do imic = 1,spproj_parts(ipart)%os_mic%get_noris()
                    cnt = cnt + 1
                    call spproj%os_mic%transfer_ori(cnt, spproj_parts(ipart)%os_mic, imic)
                enddo
                call spproj_parts(ipart)%kill
                call spproj_parts(ipart)%read_segment('stk',parts_fname(ipart))
                nstks = nstks + spproj_parts(ipart)%os_stk%get_noris()
            enddo
            if( nstks /= nmics ) THROW_HARD('Inconstistent number of stacks in individual projects')
            ! generates stacks table
            call os_stk%new(nstks, is_ptcl=.false.)
            allocate(stktab(nstks))
            cnt = 0
            do ipart = 1,params%nparts
                do istk = 1,spproj_parts(ipart)%os_stk%get_noris()
                    cnt = cnt + 1
                    call os_stk%transfer_ori(cnt, spproj_parts(ipart)%os_stk, istk)
                    stktab(cnt) = os_stk%get_str(cnt,'stk')
                enddo
                call spproj_parts(ipart)%kill
            enddo
            ! import stacks into project
            call spproj%add_stktab(stktab,os_stk)
            call os_stk%kill
            ! 2D/3D parameters, transfer everything but stack index
            cnt = 0
            do ipart = 1,params%nparts
                call spproj_parts(ipart)%read_segment('ptcl2D',parts_fname(ipart))
                call spproj_parts(ipart)%read_segment('ptcl3D',parts_fname(ipart))
                nptcls = spproj_parts(ipart)%os_ptcl2D%get_noris()
                if( nptcls /= spproj_parts(ipart)%os_ptcl3D%get_noris())then
                    THROW_HARD('Inconsistent number of particles')
                endif
                do i = 1,nptcls
                    cnt    = cnt + 1
                    stkind = spproj%os_ptcl2D%get_int(cnt,'stkind')
                    call spproj%os_ptcl2D%transfer_ori(cnt, spproj_parts(ipart)%os_ptcl2D, i)
                    call spproj%os_ptcl3D%transfer_ori(cnt, spproj_parts(ipart)%os_ptcl3D, i)
                    call spproj%os_ptcl2D%set(cnt,'stkind',stkind)
                    call spproj%os_ptcl3D%set(cnt,'stkind',stkind)
                enddo
                call spproj_parts(ipart)%kill
            enddo
        endif
        ! final write
        call spproj%write( params%projfile )
        ! clean-up
        call qsys_cleanup(params)
        call spproj%kill
        deallocate(spproj_parts,part_params)
        call o_mic%kill
        call o%kill
        ! end gracefully
        call simple_end('**** SIMPLE_REEXTRACT_DISTR NORMAL STOP ****')
        contains

            type(string) function boxfile_from_mic(mic)
                class(string), intent(in) :: mic
                type(string) :: box_from_mic
                integer      :: ibox
                box_from_mic     = fname_new_ext(basename(mic),'box')
                boxfile_from_mic = NIL
                do ibox=1,size(boxfiles)
                    if(basename(boxfiles(ibox)).eq.box_from_mic)then
                        boxfile_from_mic = boxfiles(ibox)
                        return
                    endif
                enddo
            end function boxfile_from_mic

    end subroutine exec_reextract_distr

    subroutine exec_reextract( self, cline )
        use simple_strategy2D3D_common, only: prepimgbatch, killimgbatch
        use simple_particle_extractor,  only: ptcl_extractor
        class(commander_reextract), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline !< command line input
        type(parameters)              :: params
        type(sp_project)              :: spproj, spproj_in
        type(builder)                 :: build
        type(image)                   :: micrograph, micrograph_sc
        type(ori)                     :: o_mic, o_stk
        type(ctfparams)               :: ctfparms
        type(stack_io)                :: stkio_w
        type(ptcl_extractor)          :: extractor
        type(string)                  :: mic_name, imgkind, ext
        logical,          allocatable :: mic_mask(:), ptcl_mask(:)
        integer,          allocatable :: mic2stk_inds(:), boxcoords(:,:), ptcl_inds(:), mic_dims(:,:)
        type(string)                  :: stack
        real    :: prev_shift(2),shift2d(2),shift3d(2),prev_shift_sc(2), translation(2), prev_center_sc(2)
        real    :: stk_min, stk_max, stk_mean, stk_sdev, scale
        integer :: prev_pos(2), new_pos(2), ishift(2), ldim(3), prev_center(2), new_center(2), ldim_sc(3)
        integer :: i, nframes, imic, iptcl, nmics, prev_box, box_foo, cnt, nmics_tot, stk_ind
        integer :: fromp, top, istk, nptcls2extract, nptcls
        logical :: l_3d, l_scale_particles, l_movie_frames
        call cline%set('mkdir','no')
        call params%new(cline)
        l_movie_frames    = trim(params%extractfrommov).eq.'yes'
        l_scale_particles = cline%defined('osmpd')
        if( l_scale_particles )then
            if( .not.cline%defined('box') ) THROW_HARD('BOX must be defined with OSMPD!')
            if( l_movie_frames ) THROW_HARD('Particle scaling and extraction of movie frames is not supported!')
        endif
        if( l_movie_frames .and. (trim(params%ctf).eq.'flip') )then
            THROW_HARD('extractfrommov=yes does not support ctf=flip!')
        endif
        ! set normalization radius
        params%msk = RADFRAC_NORM_EXTRACT * real(params%box/2)
        ! whether to use shifts from 2D or 3D
        l_3d = .true.
        if(cline%defined('oritype')) l_3d = trim(params%oritype)=='ptcl3D'
        ! read in integrated movies
        call spproj_in%read_segment('mic', params%projfile)
        nmics_tot = spproj_in%os_mic%get_noris()
        if( spproj_in%get_nintgs() == 0 ) THROW_HARD('No integrated micrograph to process!')
        call spproj_in%read_segment('stk', params%projfile)
        ! sanity checks, dimensions & indexing
        box_foo  = 0
        prev_box = 0
        ldim     = 0
        allocate(mic2stk_inds(nmics_tot), source=0)
        allocate(mic_mask(nmics_tot),     source=.false.)
        stk_ind = 0
        do imic = 1,nmics_tot
            if( imic > params%top ) exit
            call spproj_in%os_mic%get_ori(imic, o_mic)
            if( o_mic%isthere('state') )then
                if( o_mic%get_state() == 0 )cycle
            endif
            if( .not. o_mic%isthere('imgkind') )cycle
            if( .not. o_mic%isthere('intg')    )cycle
            call o_mic%getter('imgkind', imgkind)
            if( imgkind.ne.'mic') cycle
            ! find next selected stack
            do istk=stk_ind,spproj_in%os_stk%get_noris()
                stk_ind = stk_ind+1
                if( spproj_in%os_stk%isthere(stk_ind,'state') )then
                    if( spproj_in%os_stk%get_state(stk_ind) == 1 ) exit
                else
                    exit
                endif
            enddo
            ! update index & mask
            if( imic>=params%fromp .and. imic<=params%top )then
                mic_mask(imic) = .true.
                mic2stk_inds(imic) = stk_ind ! index to os_stk
            endif
        enddo
        nmics = count(mic_mask)
        if( nmics > 0 )then
            call build%build_general_tbox(params, cline, do3d=.false.)
            allocate(mic_dims(3,nmics_tot),source=0)
            ! sanity checks
            do imic = 1,nmics_tot
                if( .not.mic_mask(imic) )cycle
                ! sanity checks
                call spproj_in%os_mic%get_ori(imic, o_mic)
                call o_mic%getter('intg', mic_name)
                if( .not.file_exists(mic_name) )cycle
                 ! micrograph dimensions
                call find_ldim_nptcls(mic_name, ldim, nframes )
                if( nframes > 1 ) THROW_HARD('multi-frame extraction not supported; exec_reextract')
                mic_dims(:,imic) = [ldim(1),ldim(2),1]
                if( l_scale_particles )then
                    ! the following checks are not performed
                else
                    stk_ind = mic2stk_inds(imic)
                    call spproj_in%os_stk%get_ori(stk_ind, o_stk)
                    box_foo = o_stk%get_int('box')
                    if( prev_box == 0 ) prev_box = box_foo
                    if( prev_box /= box_foo ) THROW_HARD('Inconsistent box size; exec_reextract')
                endif
            enddo
            if( .not.cline%defined('box') ) params%box = prev_box
            if( is_odd(params%box) ) THROW_HARD('Box size must be of even dimension! exec_extract')
            ! extraction
            write(logfhandle,'(A)')'>>> EXTRACTING... '
            call spproj_in%read_segment('ptcl2D', params%projfile)
            call spproj_in%read_segment('ptcl3D', params%projfile)
            allocate(ptcl_mask(spproj_in%os_ptcl2D%get_noris()),source=.false.)
            ldim = 0
            do imic = params%fromp,params%top
                if( .not.mic_mask(imic) ) cycle
                ! init micrograph object and extractor
                call spproj_in%os_mic%get_ori(imic, o_mic)
                call o_mic%getter('intg', mic_name)
                ctfparms = o_mic%get_ctfvars()
                if( any(ldim /= mic_dims(:,imic)) )then
                    ! first iteration or different micrograph size
                    ldim = mic_dims(:,imic)
                    call micrograph%new(ldim, ctfparms%smpd)
                    if( .not.l_movie_frames ) call extractor%init_mic(params%box, (params%pcontrast .eq. 'black'))
                    if( l_scale_particles )then
                        scale        = ctfparms%smpd / params%osmpd
                        ldim_sc(1:2) = round2even(scale*real(ldim(1:2)))
                        ldim_sc(3)   = 1
                        call micrograph_sc%new(ldim_sc, params%osmpd)
                    endif
                endif
                ! stack
                stk_ind = mic2stk_inds(imic)
                call spproj_in%os_stk%get_ori(stk_ind, o_stk)
                prev_box = o_stk%get_int('box')
                fromp    = o_stk%get_fromp()
                top      = o_stk%get_top()
                ext      = fname2ext(basename(mic_name))
                stack    = string(EXTRACT_STK_FBODY)//get_fbody(basename(mic_name), ext)//STK_EXT
                ! updating shifts, positions, states and doc
                if( allocated(boxcoords) ) deallocate(boxcoords)
                allocate(boxcoords(2,fromp:top),source=0)
                !$omp parallel do default(shared) proc_bind(close) schedule(static)&
                !$omp private(iptcl,prev_pos,prev_shift,prev_center,prev_center_sc,prev_shift_sc,new_center)&
                !$omp private(new_pos,translation,shift2d,shift3d,ishift)
                do iptcl = fromp,top
                    if( spproj_in%os_ptcl2D%get_state(iptcl) == 0 ) cycle
                    if( spproj_in%os_ptcl3D%get_state(iptcl) == 0 ) cycle
                    ! previous position & shift
                    call spproj_in%get_boxcoords(iptcl, prev_pos)
                    if( l_3d )then
                        prev_shift = spproj_in%os_ptcl3D%get_2Dshift(iptcl)
                    else
                        prev_shift = spproj_in%os_ptcl2D%get_2Dshift(iptcl)
                    endif
                    if( l_scale_particles )then
                        ! scale center, shift & positions
                        prev_center      = prev_pos + prev_box/2
                        prev_center_sc   = scale * real(prev_center)
                        prev_shift_sc    = scale * real(prev_shift)
                        new_center       = nint(prev_center_sc - prev_shift_sc)
                        new_pos          = new_center - params%box/2
                        translation      = -(prev_center_sc - real(new_center))
                        ptcl_mask(iptcl) = box_inside(ldim, new_pos, params%box)
                        if( ptcl_mask(iptcl) )then
                            ! updates shifts
                            if( l_3d )then
                                shift2d = scale * spproj_in%os_ptcl2D%get_2Dshift(iptcl) + translation
                                shift3d = scale * prev_shift                             + translation
                            else
                                shift2d = scale * prev_shift                             + translation
                                shift3d = scale * spproj_in%os_ptcl3D%get_2Dshift(iptcl) + translation
                            endif
                        endif
                    else
                        ! calc new position & shift
                        ishift      = nint(prev_shift)
                        new_pos     = prev_pos - ishift
                        translation = -real(ishift)
                        if( prev_box /= params%box ) new_pos = new_pos + (prev_box-params%box)/2
                        ptcl_mask(iptcl) = box_inside(ldim, new_pos, params%box)
                        if( ptcl_mask(iptcl) )then
                            ! updates shifts
                            if( l_3d )then
                                shift2d = spproj_in%os_ptcl2D%get_2Dshift(iptcl) + translation
                                shift3d = prev_shift                             + translation
                            else
                                shift2d = prev_shift                             + translation
                                shift3d = spproj_in%os_ptcl3D%get_2Dshift(iptcl) + translation
                            endif
                        endif
                    endif
                    ! updates document
                    if( ptcl_mask(iptcl) )then
                        ! updates picking position
                        call spproj_in%set_boxcoords(iptcl, new_pos)
                        ! updates shifts
                        call spproj_in%os_ptcl2D%set_shift(iptcl, shift2d)
                        call spproj_in%os_ptcl3D%set_shift(iptcl, shift3d)
                    else
                        ! excluded
                        call spproj_in%os_ptcl2D%set_state(iptcl, 0)
                        call spproj_in%os_ptcl3D%set_state(iptcl, 0)
                    endif
                    ! for actual extraction
                    boxcoords(:,iptcl) = new_pos
                enddo
                !$omp end parallel do
                nptcls2extract = count(ptcl_mask(fromp:top))
                if( nptcls2extract > 0 )then
                    if( allocated(ptcl_inds) ) deallocate(ptcl_inds)
                    allocate(ptcl_inds(nptcls2extract),source=0)
                    cnt = 0
                    do iptcl = fromp,top
                        if( .not.ptcl_mask(iptcl) ) cycle
                        cnt = cnt + 1
                        ptcl_inds(cnt) = iptcl
                        ! updating index of particle in stack
                        call spproj_in%os_ptcl2D%set(iptcl, 'indstk', cnt)
                        call spproj_in%os_ptcl3D%set(iptcl, 'indstk', cnt)
                    enddo
                    ptcl_inds = ptcl_inds -fromp+1 ! because indexing range lost when passed to extractor
                    call prepimgbatch(params, build, nptcls2extract)
                    if( l_movie_frames )then
                        ! extraction from movie
                        call extractor%init_mov(o_mic, params%box, (params%pcontrast .eq. 'black'))
                        call extractor%extract_particles(ptcl_inds, boxcoords, build%imgbatch, stk_min,stk_max,stk_mean,stk_sdev)
                    else
                        ! read micrograph
                        call micrograph%read(mic_name)
                        ! preprocess micrograph
                        if( trim(params%backgr_subtr).eq.'yes') call micrograph%subtract_background(HP_BACKGR_SUBTR)
                        ! Actual extraction
                        if( l_scale_particles )then
                            ! scale
                            if( all(ldim_sc == ldim) )then
                                call micrograph_sc%copy_fast(micrograph)
                            else
                                call micrograph%fft
                                if( any(ldim_sc < ldim) )then
                                    call micrograph%clip(micrograph_sc)
                                else
                                    call micrograph%pad(micrograph_sc, antialiasing=.false.)
                                endif
                                call micrograph_sc%ifft
                                call micrograph%set_ft(.false.)
                            endif
                            call micrograph_sc%set_smpd(params%osmpd) ! safety
                            ! extract
                            call extractor%extract_particles_from_mic(micrograph_sc, ptcl_inds, boxcoords, build%imgbatch,&
                            &stk_min,stk_max,stk_mean,stk_sdev)
                        else
                            ! extract
                            call extractor%extract_particles_from_mic(micrograph, ptcl_inds, boxcoords, build%imgbatch,&
                                &stk_min,stk_max,stk_mean,stk_sdev)
                        endif
                    endif
                    ! write stack
                    if( l_scale_particles )then
                        call stkio_w%open(stack, params%osmpd, 'write', box=params%box)
                    else
                        call stkio_w%open(stack, params%smpd, 'write', box=params%box)
                    endif
                    do i = 1,nptcls2extract
                        call stkio_w%write(i, build%imgbatch(i))
                    enddo
                    call stkio_w%close
                    if( l_scale_particles )then
                        call spproj_in%os_stk%set(stk_ind,'smpd',params%osmpd) !!
                        call micrograph_sc%update_header_stats(stack, [stk_min, stk_max, stk_mean, stk_sdev])
                    else
                        call micrograph%update_header_stats(stack, [stk_min, stk_max, stk_mean, stk_sdev])
                    endif
                    call spproj_in%os_stk%set(stk_ind,'stk',   simple_abspath(stack,check_exists=.false.))
                    call spproj_in%os_stk%set(stk_ind,'box',   params%box)
                    call spproj_in%os_stk%set(stk_ind,'nptcls',nptcls2extract)
                    call spproj_in%os_mic%set(imic,   'nptcls',nptcls2extract)
                    call spproj_in%os_mic%delete_entry(imic,'boxfile')
                else
                    ! all particles in this micrograph excluded
                    call spproj_in%os_stk%set(stk_ind,'state',0)
                    call spproj_in%os_mic%set(imic,'state',0)
                    mic_mask(imic) = .false.
                    mic2stk_inds(imic) = 0
                endif
            enddo
        endif
        call micrograph%kill
        call micrograph_sc%kill
        call extractor%kill
        call killimgbatch(build)
        ! OUTPUT
        call spproj%read_non_data_segments(params%projfile)
        call spproj%projinfo%set(1,'projname', get_fbody(params%outfile,METADATA_EXT,separator=.false.))
        call spproj%projinfo%set(1,'projfile', params%outfile)
        nmics = count(mic_mask)
        ! transfer mics & stk
        call spproj%os_mic%new(nmics, is_ptcl=.false.)
        call spproj%os_stk%new(nmics, is_ptcl=.false.)
        nptcls = count(ptcl_mask)
        cnt = 0
        do imic = params%fromp,params%top
            if( .not.mic_mask(imic) )cycle
            cnt = cnt+1
            call spproj%os_mic%transfer_ori(cnt, spproj_in%os_mic, imic)
            stk_ind = mic2stk_inds(imic)
            call spproj%os_stk%transfer_ori(cnt, spproj_in%os_stk, stk_ind)
        enddo
        ! transfer particles
        nptcls = count(ptcl_mask)
        call spproj%os_ptcl2D%new(nptcls, is_ptcl=.true.)
        call spproj%os_ptcl3D%new(nptcls, is_ptcl=.true.)
        cnt = 0
        do iptcl = 1,size(ptcl_mask)
            if( .not.ptcl_mask(iptcl) )cycle
            cnt = cnt+1
            call spproj%os_ptcl2D%transfer_ori(cnt, spproj_in%os_ptcl2D, iptcl)
            call spproj%os_ptcl3D%transfer_ori(cnt, spproj_in%os_ptcl3D, iptcl)
        enddo
        call spproj_in%kill
        ! final write
        call spproj%write(params%outfile)
        write(logfhandle,'(A,I8)')'>>> RE-EXTRACTED  PARTICLES: ', nptcls
        ! end gracefully
        call qsys_job_finished(params, string('simple_commanders_pick :: exec_reextract'))
        call build%kill_general_tbox
        call o_mic%kill
        call o_stk%kill
        call simple_end('**** SIMPLE_REEXTRACT NORMAL STOP ****')
    end subroutine exec_reextract

    ! Stream only application
    subroutine exec_pick_extract( self, cline )
        use simple_picker_iter, only: picker_iter
        class(commander_pick_extract), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(parameters)              :: params
        type(oris)                    :: os_mic
        type(ori)                     :: o_mic
        type(picker_iter)             :: piter
        type(commander_extract)       :: xextract
        type(cmdline)                 :: cline_extract
        type(sp_project)              :: spproj
        type(string) :: micname, output_dir_picker, fbody, output_dir_extract
        type(string) :: boxfile, thumb_den
        integer :: fromto(2), imic, ntot, state, nvalid, i, nptcls
        logical :: l_extract
        ! set oritype
        call cline%set('oritype', 'mic')
        ! parse parameters
        call params%new(cline)
        ! if( params%stream.ne.'yes' ) THROW_HARD('new streaming only application')
        l_extract   = trim(params%extract).eq.'yes'
        ! read in movies
        call spproj%read( params%projfile )
        if( spproj%get_nintgs() == 0 ) THROW_HARD('no micrograph to process!')
        params%smpd = spproj%os_mic%get(1,'smpd')
        call cline%set('smpd',params%smpd)
        ! output directories
        output_dir_picker  = DIR_PICKER
        if( l_extract ) output_dir_extract = DIR_EXTRACT
        if( cline%defined('dir') )then
            output_dir_picker  = filepath(params%dir,output_dir_picker)//'/'
            if( l_extract) output_dir_extract = filepath(params%dir,output_dir_extract)//'/'
        endif
        call simple_mkdir(output_dir_picker)
        if( l_extract ) call simple_mkdir(output_dir_extract)
        ! picker specs
        select case(trim(params%picker))
            case('new')
                if(cline%defined('pickrefs'))then
                else
                    if( .not.cline%defined('moldiam') )then
                        THROW_HARD('MOLDIAM required for picker=new')
                    endif
                endif
            case('segdiam')
                if( .not.cline%defined('moldiam_max') )then
                    THROW_WARN('MOLDIAM_MAX not set on command line, falling back on default value: '//int2str(int(params%moldiam_max))//' A')
                endif
            case DEFAULT
                THROW_HARD('Unsupported PICKER: '//trim(params%picker))
        end select
        ! command lines
        if( l_extract )then
            cline_extract = cline
            call cline_extract%set('dir', output_dir_extract)
            call cline_extract%set('pcontrast', params%pcontrast)
            if( cline%defined('box_extract') ) call cline_extract%set('box', params%box_extract)
            call cline%delete('box')
            call cline_extract%delete('box_extract')
        endif
        ! file name
        if( cline%defined('fbody') )then
            fbody = params%fbody
        else
            fbody = ''
        endif
        ! range
        fromto(:) = 1
        if( cline%defined('fromp') .and. cline%defined('top') )then
            fromto(:) = [params%fromp, params%top]
        endif
        ntot   = fromto(2) - fromto(1) + 1
        nvalid = 0
        ! main loop
        do imic = fromto(1),fromto(2)
            ! fetch movie orientation
            call spproj%os_mic%get_ori(imic, o_mic)
            ! sanity check
            state = 1
            if( o_mic%isthere('state') ) state = o_mic%get_state()
            if( state == 0 ) cycle
            if( .not.o_mic%isthere('intg')   )cycle
            call o_mic%getter('intg', micname)
            if( .not.file_exists(micname)) cycle
            ! picker
            params%lp = max(params%fny, params%lp_pick)
            call piter%iterate(params, cline, params%smpd, micname, output_dir_picker, boxfile, thumb_den, nptcls)
            call o_mic%set('nptcls', nptcls)
            if( nptcls > 0 )then
                call o_mic%set('boxfile', boxfile)
                call o_mic%set('thumb_den', thumb_den)
            else
                call o_mic%set_state(0)
            endif
            ! update project
            call spproj%os_mic%set_ori(imic, o_mic)
            nvalid = nvalid+1
        enddo
        ! extract particles
        if( l_extract )then
            call spproj%write_segment_inside(params%oritype, params%projfile)
            call xextract%execute(cline_extract)
            ! nothing to write, done by extract
        else
            if( ntot > 1 )then
                ! purging state=0 and nptcls=0 mics such that all mics (nmics>1)
                ! can be assumed valid
                call os_mic%new(nvalid, is_ptcl=.false.)
                i = 0
                do imic = fromto(1),fromto(2)
                    state  = spproj%os_mic%get_state(imic)
                    nptcls = spproj%os_mic%get_int(imic,'nptcls')
                    if( (state == 1) .and. (nptcls > 0) )then
                        i = i+1
                        call os_mic%transfer_ori(i, spproj%os_mic, imic)
                    endif
                enddo
                spproj%os_mic = os_mic
                call os_mic%kill
            endif
            call spproj%write_segment_inside(params%oritype, params%projfile)
        endif
        ! end gracefully
        call qsys_job_finished(params, string('simple_commanders_pick :: exec_pick_extract'))
        call o_mic%kill
        call piter%kill
        call simple_end('**** SIMPLE_PICK_EXTRACT NORMAL STOP ****')
    end subroutine exec_pick_extract

    subroutine exec_make_pickrefs( self, cline )
        class(commander_make_pickrefs), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(parameters)         :: params
        type(stack_io)           :: stkio_r
        type(oris)               :: moldiamori
        type(image)              :: ref2D, ref2D_clip
        type(image), allocatable :: projs(:), masks(:)
        real,        allocatable :: diams(:), shifts(:,:)
        real,    parameter :: MSKDIAM2LP = 0.15, lP_LB = 30., LP_UB = 15.
        integer, parameter :: NREFS=100
        real    :: ang, rot, lp, diam_max, maxdiam, moldiam, mskdiam
        integer :: nrots, iref, irot, ldim_clip(3), ldim(3), ncavgs, icavg
        integer :: cnt, norefs, box_for_pick, box_for_extract
        ! error check
        if( cline%defined('vol1')          ) THROW_HARD('vol1 input no longer supported, use prg=reproject to generate 20 2D references')
        if( .not.cline%defined('pickrefs') ) THROW_HARD('PICKREFS must be informed!')
        ! set defaults
        call set_automask2D_defaults(cline)
        call cline%set('oritype', 'mic')
        if( .not.cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        ! parse parameters
        call params%new(cline)
        if( params%stream.eq.'yes' ) THROW_HARD('not a streaming application')
        ! read selected cavgs
        call find_ldim_nptcls(params%pickrefs, ldim, ncavgs)
        ldim(3) = 1
        params%msk = real(ldim(1)/2) - COSMSKHALFWIDTH ! for automasking
        ! read
        allocate( projs(ncavgs), masks(ncavgs) )
        call stkio_r%open(params%pickrefs, params%smpd, 'read', bufsz=ncavgs)
        do icavg=1,ncavgs
            call projs(icavg)%new(ldim, params%smpd)
            call stkio_r%read(icavg, projs(icavg))
            call masks(icavg)%copy(projs(icavg))
        end do
        call stkio_r%close
        ! Automasking
        call automask2D(params, masks, params%ngrow, nint(params%winsz), params%edge, diams, shifts)
        do icavg=1,ncavgs
            call projs(icavg)%div_below(0.,10.)
            call projs(icavg)%mul(masks(icavg))
            call projs(icavg)%shift([shifts(icavg,1),shifts(icavg,2),0.])
        end do
        ! estimate new box size and clip
        diam_max        = maxval(diams)
        lp              = min(max(LP_LB,MSKDIAM2LP * diam_max),LP_UB)
        box_for_pick    = min(round2even(diam_max / params%smpd + 2. * COSMSKHALFWIDTH), ldim(1))
        moldiam         = params%smpd * box_for_pick
        mskdiam         = moldiam * MSK_EXP_FAC
        maxdiam         = moldiam + moldiam * BOX_EXP_FAC
        box_for_extract = find_larger_magic_box(round2even(maxdiam / params%smpd))
        write(logfhandle,'(A,1X,I4)') 'ESTIMATED BOX SIZE: ', box_for_pick
        ! set mskdiam in cline
        call cline%set('mskdiam', mskdiam)
        ! write info to file
        call moldiamori%new(1,                    .false.)
        call moldiamori%set(1, 'diam_max',        diam_max)
        call moldiamori%set(1, 'lp',              lp)
        call moldiamori%set(1, 'box_for_pick',    box_for_pick)
        call moldiamori%set(1, 'moldiam',         moldiam)
        call moldiamori%set(1, 'mskdiam',         mskdiam)
        call moldiamori%set(1, 'box_for_extract', box_for_extract)
        if(file_exists(STREAM_MOLDIAM)) call del_file(STREAM_MOLDIAM)
        call moldiamori%write(1, string(STREAM_MOLDIAM))
        call moldiamori%kill
        ! expand in in-plane rotation, clip and write to file
        ldim_clip = [box_for_pick, box_for_pick, 1]
        do icavg=1,ncavgs
            call projs(icavg)%bp(0.,lp)
        end do
        nrots  = nint( real(NREFS)/real(ncavgs) )
        norefs = ncavgs
        call ref2D_clip%new([ldim_clip(1),ldim_clip(2),1], params%smpd)
        if( nrots > 1 )then
            call ref2D%new([ldim(1),ldim(2),1], params%smpd)
            ang = 360./real(nrots)
            cnt = 0
            do iref=1,norefs
                rot = 0.
                do irot=1,nrots
                    cnt = cnt + 1
                    call projs(iref)%rtsq(rot, 0., 0., ref2D)
                    call ref2D%clip(ref2D_clip)
                    call ref2D_clip%write(string(PICKREFS_FBODY)//params%ext, cnt)
                    rot = rot + ang
                end do
            end do
        else
            ! should never happen
            do iref=1,norefs
                call projs(iref)%clip(ref2D_clip)
                call ref2D_clip%write(string(PICKREFS_FBODY)//params%ext, iref)
            end do
        endif
        ! cleanup
        do icavg = 1,ncavgs
            call masks(icavg)%kill
            call projs(icavg)%kill
        enddo
        deallocate(masks,projs)
        call ref2D%kill
        call ref2D_clip%kill
        ! end gracefully
        call simple_touch('MAKE_PICKREFS_FINISHED')
        call simple_end('**** SIMPLE_MAKE_PICKREFS NORMAL STOP ****')
    end subroutine exec_make_pickrefs

    logical function box_inside( ildim, coord, box )
        integer, intent(in) :: ildim(3), coord(2), box
        integer             :: fromc(2), toc(2)
        fromc  = coord+1       ! compensate for the c-range that starts at 0
        toc    = fromc+(box-1) ! the lower left corner is 1,1
        box_inside = .true.    ! box is inside
        if( any(fromc < 1) .or. toc(1) > ildim(1) .or. toc(2) > ildim(2) ) box_inside = .false.
    end function box_inside

end module simple_commanders_pick
