!@descr: project file utilities
module simple_projfile_utils
use simple_core_module_api
use simple_image,         only: image
use simple_sp_project,    only: sp_project
use simple_oris,          only: oris
use simple_euclid_sigma2, only: average_sigma2_groups
use simple_class_frcs
implicit none
#include "simple_local_flags.inc"

contains

    subroutine merge_chunk_projfiles( chunk_fnames, folder, merged_proj, projname_out, write_proj, cavgs_out, cavgs_replace, sigma2_out, update_classno )
        class(string),           intent(in)    :: chunk_fnames(:) ! List of project files
        class(string),           intent(in)    :: folder          ! output folder
        class(sp_project),       intent(inout) :: merged_proj     ! output project, assumed to have compuational env info
        class(string), optional, intent(in)    :: projname_out    ! name for output project file
        logical,       optional, intent(in)    :: write_proj      ! write project file
        logical,       optional, intent(in)    :: cavgs_replace   ! replace cavgs
        logical,       optional, intent(in)    :: update_classno  ! update the ptcl class numbers
        class(string), optional, intent(in)    :: cavgs_out       ! name for output cls2D stack
        class(string), optional, intent(in)    :: sigma2_out      ! name for combined sigma2 file
        type(sp_project), allocatable :: chunks(:)
        type(string),     allocatable :: chunks_sigma2(:)
        real,             allocatable :: states(:)
        integer,          allocatable :: clsmap(:)
        type(class_frcs) :: frcs, frcs_chunk
        type(image)      :: img
        type(string)     :: projname, stkname, evenname, oddname, frc_fname, projfile_out, dir, cavgs
        character(len=STDLEN) :: imgkind_here
        type(string)     :: cavgs_tmp, evenname_tmp, oddname_tmp, sigma2_fname
        real             :: smpd
        integer          :: ldim(3), i, ic, icls, ncls, nchunks, nallmics, nallstks, nallptcls, ncls_tot, box4frc
        integer          :: fromp, fromp_glob, top, top_glob, j, iptcl_glob, nstks, nmics, nptcls, istk
        logical          :: l_write_proj, l_cavgs_replace, l_update_classno, l_merge_evenodd, l_merge_frcs, l_merge_sigma2
        logical          :: frcs_initialised
        l_write_proj     = .true.
        l_cavgs_replace  = .false.
        l_update_classno = .true.
        frcs_initialised = .false.
        if( present( write_proj )    ) l_write_proj     = write_proj
        if( present( cavgs_replace  )) l_cavgs_replace  = cavgs_replace
        if( present( update_classno )) l_update_classno = update_classno
        nchunks = size(chunk_fnames)
        allocate(chunks(nchunks))
        allocate(chunks_sigma2(nchunks))
        dir = folder%to_char()//'/'
        if( present(projname_out) )then
            projfile_out = dir%to_char()//projname_out%to_char()//trim(METADATA_EXT)
        else
            projfile_out = dir%to_char()//'set'//METADATA_EXT
        endif
        if( present(sigma2_out) )then
            sigma2_fname   = dir%to_char()//sigma2_out%to_char()//trim(STAR_EXT)
        else
            sigma2_fname   = dir%to_char()//'sigma2_combined'//trim(STAR_EXT)
        endif
        call merged_proj%os_mic%kill
        call merged_proj%os_stk%kill
        call merged_proj%os_ptcl2D%kill
        call merged_proj%os_ptcl3D%kill
        call merged_proj%os_cls2D%kill
        call merged_proj%os_cls3D%kill
        call merged_proj%os_out%kill
        if( present(cavgs_out) )then
            cavgs = dir//cavgs_out
        else
            cavgs = dir%to_char()//'cavgs'//MRC_EXT
        endif
        if( l_cavgs_replace ) then
            cavgs_tmp    = cavgs
            evenname_tmp = dir//get_fbody(basename(cavgs), fname2ext(cavgs))//'_even'//MRC_EXT
            oddname_tmp  = dir//get_fbody(basename(cavgs), fname2ext(cavgs))//'_odd'//MRC_EXT
            cavgs        = dir%to_char()//'cavgs_tmp'//MRC_EXT
        endif
        if( file_exists(cavgs) )THROW_WARN('ouput stack already exists: '//cavgs%to_char())
        nallptcls = 0
        nallstks  = 0
        nallmics  = 0
        icls      = 0
        l_merge_evenodd = .true.
        l_merge_frcs    = .true.
        l_merge_sigma2  = .true.
        do ic = 1,nchunks
            projname = chunk_fnames(ic)
            call chunks(ic)%read_data_info(projname, nmics, nstks, nptcls)
            nallmics  = nallmics  + nmics
            nallstks  = nallstks  + nstks
            nallptcls = nallptcls + nptcls
            call chunks(ic)%read_segment('out',  projname)
            call chunks(ic)%read_segment('cls2D',  projname)
            if( chunks(ic)%os_out%get_noris() == 0 .or. chunks(ic)%os_cls2D%get_noris() == 0 ) then
                l_merge_evenodd = .false.
                l_merge_frcs    = .false.
                l_merge_sigma2  = .false.
                cycle
            end if
            call chunks(ic)%get_cavgs_stk(stkname, ncls, smpd, imgkind='cavg')
            call find_ldim_nptcls(stkname, ldim, ncls)
            ldim(3) = 1
            ncls    = chunks(ic)%os_cls2D%get_noris()
            call img%new(ldim, smpd)
            evenname = add2fbody(stkname, MRC_EXT, '_even')
            oddname  = add2fbody(stkname, MRC_EXT, '_odd')
            l_merge_evenodd = l_merge_evenodd .and. file_exists(evenname) .and. file_exists(oddname)
            call chunks(ic)%get_frcs(frc_fname, 'frc2D', fail=.false.)
            l_merge_frcs = l_merge_frcs .and. frc_fname /= NIL .and. file_exists(frc_fname)
            do i = 1,chunks(ic)%os_out%get_noris()
                if( chunks(ic)%os_out%isthere(i,'imgkind') )then
                    call chunks(ic)%os_out%get_static(i, 'imgkind', imgkind_here)
                    if( trim(imgkind_here) == 'sigma2' ) exit
                endif
            end do
            if( i > chunks(ic)%os_out%get_noris() ) l_merge_sigma2 = .false.
            do i = 1,ncls
                icls = icls+1
                call img%read(stkname,i)
                call img%write(cavgs,icls)
                if( l_merge_evenodd )then
                    call img%read(evenname,i)
                    call img%write(dir//get_fbody(basename(cavgs), fname2ext(cavgs))//'_even'//MRC_EXT,icls)
                    call img%read(oddname,i)
                    call img%write(dir//get_fbody(basename(cavgs), fname2ext(cavgs))//'_odd'//MRC_EXT,icls)
                endif
            enddo
        enddo
        call img%kill
        if( .not. l_merge_evenodd ) THROW_WARN('merge_chunk_projfiles: missing even/odd class-average stacks; skipping even/odd merge')
        if( .not. l_merge_frcs )    THROW_WARN('merge_chunk_projfiles: missing frc2D data; skipping FRC merge')
        if( .not. l_merge_sigma2 )  THROW_WARN('merge_chunk_projfiles: missing sigma2 data; skipping sigma2 merge')
        ncls_tot = icls
        ! micrographs
        if( nallmics > 0 )then
            call merged_proj%os_mic%new(nallmics,.false.)
            j = 0
            do ic = 1,nchunks
                projname = chunk_fnames(ic)
                call chunks(ic)%read_segment('mic', projname)
                do i = 1,chunks(ic)%os_mic%get_noris()
                    j = j+1
                    call merged_proj%os_mic%transfer_ori(j, chunks(ic)%os_mic, i)
                enddo
                call chunks(ic)%os_mic%kill
            enddo
        endif
        ! particles, stacks and classes frcs & metadata
        call merged_proj%os_cls2D%new(ncls_tot,  .false.)
        call merged_proj%os_ptcl2D%new(nallptcls,.true.)
        call merged_proj%os_stk%new(nallstks,    .false.)
        icls       = 0
        istk       = 0
        iptcl_glob = 0
        fromp_glob = 1
        do ic = 1,nchunks
            projname = chunk_fnames(ic)
            call chunks(ic)%read_segment('stk', projname)
            call absolutize_project_stack_paths(chunks(ic), simple_abspath(projname, check_exists=.false.))
            call chunks(ic)%read_segment('ptcl2D',projname)
            ! classes frcs & info
            ncls = chunks(ic)%os_cls2D%get_noris()
            if( l_merge_frcs )then
                call chunks(ic)%get_frcs(frc_fname, 'frc2D')
                call frcs_chunk%read(frc_fname)
                if( .not. frcs_initialised )then
                    box4frc = frcs_chunk%get_box()
                    call frcs%new(ncls_tot, box4frc, smpd)
                    frcs_initialised = .true.
                endif
            endif
            allocate(clsmap(ncls),source=0)
            do i = 1,ncls
                icls      = icls+1
                clsmap(i) = icls
                call merged_proj%os_cls2D%transfer_ori(icls,   chunks(ic)%os_cls2D, i)
                call merged_proj%os_cls2D%set_class(icls, icls)
                call merged_proj%os_cls2D%set(icls,'origclass',i)
                call merged_proj%os_cls2D%set(icls,'chunk',    projname)
                if( l_merge_frcs .and. chunks(ic)%os_cls2D%get_state(i) > 0 )then
                    call frcs%set_frc(icls, frcs_chunk%get_frc(i,  box4frc))
                endif
            enddo
            ! particles and stacks
            nstks  = chunks(ic)%os_stk%get_noris()
            do i = 1,nstks
                istk  = istk + 1
                fromp = chunks(ic)%os_stk%get_fromp(i)
                top   = chunks(ic)%os_stk%get_top(i)
                do j = fromp,top
                    iptcl_glob = iptcl_glob + 1
                    if(l_update_classno .and. chunks(ic)%os_ptcl2D%get_class(j) > 0) call chunks(ic)%os_ptcl2D%set_class(j, clsmap(chunks(ic)%os_ptcl2D%get_class(j)))
                    call chunks(ic)%os_ptcl2D%set_stkind(j, istk)
                    call merged_proj%os_ptcl2D%transfer_ori(iptcl_glob, chunks(ic)%os_ptcl2D, j)
                    if( chunks(ic)%os_ptcl2D%get_state(j) == 0 ) call merged_proj%os_ptcl2D%reject(iptcl_glob)
                enddo
                top_glob = fromp_glob + top - fromp
                call chunks(ic)%os_stk%set(i, 'fromp', fromp_glob)
                call chunks(ic)%os_stk%set(i, 'top',   top_glob)
                fromp_glob = top_glob+1
                call merged_proj%os_stk%transfer_ori(istk, chunks(ic)%os_stk, i)
            enddo
            deallocate(clsmap)
            ! sigma2
            if( l_merge_sigma2 )then
                call chunks(ic)%get_sigma2(chunks_sigma2(ic))
            else
                chunks_sigma2(ic) = NIL
            endif
            ! making sure the compenv is informed
            if( (ic == 1) .and. (merged_proj%compenv%get_noris() == 0) )then
                call chunks(1)%read_non_data_segments(chunk_fnames(1))
                merged_proj%compenv = chunks(1)%compenv
                call merged_proj%update_projinfo(projfile_out)
            endif
            ! cleanup
            call chunks(ic)%kill
        enddo
        deallocate(chunks)
        ! add classes, frcs
        if( l_merge_frcs )then
            call frcs%write(dir//trim(FRCS_FILE))
            call merged_proj%add_frcs2os_out(dir//trim(FRCS_FILE), 'frc2D')
        endif
        if( l_cavgs_replace ) then
            call simple_rename( cavgs,    cavgs_tmp,    overwrite=.true. )
            if( l_merge_evenodd )then
                evenname = dir//get_fbody(basename(cavgs), fname2ext(cavgs))//'_even'//MRC_EXT
                oddname  = dir//get_fbody(basename(cavgs), fname2ext(cavgs))//'_odd'//MRC_EXT
                call simple_rename( evenname, evenname_tmp, overwrite=.true. )
                call simple_rename( oddname,  oddname_tmp,  overwrite=.true. )
            endif
            cavgs = cavgs_tmp
        endif
        call merged_proj%add_cavgs2os_out(cavgs, smpd, imgkind='cavg')
        ! merge and add sigmas
        if( l_merge_sigma2 )then
            call average_sigma2_groups(sigma2_fname, chunks_sigma2)
            if( file_exists(sigma2_fname) ) call merged_proj%add_sigma22os_out(sigma2_fname)
        endif
        deallocate(chunks_sigma2)
        ! propagate 2D states to 3D
        states = merged_proj%os_cls2D%get_all('state')
        call merged_proj%os_cls3D%new(ncls_tot, .false.)
        call merged_proj%os_cls3D%set_all('state', states)
        merged_proj%os_ptcl3D = merged_proj%os_ptcl2D
        call merged_proj%os_ptcl3D%delete_2Dclustering
        ! write
        if(l_write_proj) call merged_proj%write(projfile_out)
        ! cleanup
        if( l_merge_frcs )then
            call frcs%kill
            call frcs_chunk%kill
        endif
    end subroutine merge_chunk_projfiles

    subroutine merge_selected_project_files( project_fnames, projfile_out, merged_proj, write_proj )
        class(string),     intent(in)    :: project_fnames(:) ! SIMPLE project files to merge
        class(string),     intent(in)    :: projfile_out      ! output project file
        class(sp_project), intent(inout) :: merged_proj       ! output project
        logical, optional, intent(in)    :: write_proj        ! write project file
        real, parameter :: SMPD_TOL = 0.001
        type(sp_project), allocatable :: projects(:)
        type(string) :: projfile_abs, projdir, stage_dir, projfile_stage
        type(binoris_seginfo), allocatable :: seginfos(:), hint_infos(:)
        integer, allocatable :: nmics(:), nstks(:), nptcl2Ds(:), nptcl3Ds(:), noptics(:)
        integer, allocatable :: mic_offsets(:), stk_offsets(:), ptcl2D_offsets(:), ptcl3D_offsets(:)
        integer, allocatable :: opt_offsets(:)
        integer, allocatable :: src_ogid_maxs(:), ogid_offsets(:)
        integer, allocatable :: seginds(:), hint_inds(:)
        integer :: nprojs, iproj, imic, istk, iptcl, iopt, iseg, idir_dummy
        integer :: imic_glob, istk_glob, iptcl2D_glob, iptcl3D_glob
        integer :: iopt_glob
        integer :: stk_offset, ptcl2D_offset, ptcl3D_offset
        integer :: fromp, top, range_offset
        integer :: ogid_offset, ogid_glob_max, ogid
        logical :: l_write_proj, l_has_mics, l_has_stks, l_has_ptcl2D, l_has_ptcl3D
        logical :: l_has_optics, l_has_any_data
        nprojs = size(project_fnames)
        if( nprojs < 2 ) THROW_HARD('merge_selected_project_files requires at least two input projects')
        if( fname2format(projfile_out) /= 'O' )then
            THROW_HARD('output file must be a SIMPLE project (*.simple): '//trim(projfile_out%to_char()))
        endif
        l_write_proj = .false.
        if( present(write_proj) ) l_write_proj = write_proj
        allocate(projects(nprojs), nmics(nprojs), nstks(nprojs), nptcl2Ds(nprojs), nptcl3Ds(nprojs), noptics(nprojs))
        do iproj = 1,nprojs
            projfile_abs = simple_abspath(project_fnames(iproj))
            call projects(iproj)%read_segments_info(projfile_abs, seginds, seginfos)
            if( .not.data_segments_present(seginds) )then
                write(logfhandle,*) 'projtab entry: ', projfile_abs%to_char()
                write(logfhandle,*) 'This SIMPLE project file has no mergeable mic/stk/ptcl/optics data segments.'
                write(logfhandle,*) 'Analysis products cls2D, cls3D, and out are ignored by merge_projects.'
                projdir = get_fpath(projfile_abs)
                stage_dir  = ''
                idir_dummy = find_next_int_dir_prefix(projdir, stage_dir)
                if( idir_dummy > 1 .and. stage_dir /= '' )then
                    projfile_stage = projdir//stage_dir//'/'//basename(projfile_abs)
                    if( file_exists(projfile_stage) )then
                        call projects(iproj)%read_segments_info(projfile_stage, hint_inds, hint_infos)
                        if( data_segments_present(hint_inds) )then
                            write(logfhandle,*) 'Did you mean to list this data-bearing stage project?'
                            write(logfhandle,*) projfile_stage%to_char()
                        endif
                        if( allocated(hint_inds)  ) deallocate(hint_inds)
                        if( allocated(hint_infos) ) deallocate(hint_infos)
                    endif
                endif
                THROW_HARD('merge_projects projtab entry has no project data segments')
            endif
            nmics(iproj)    = 0
            nstks(iproj)    = 0
            nptcl2Ds(iproj) = 0
            nptcl3Ds(iproj) = 0
            noptics(iproj)  = 0
            if( allocated(seginds) )then
                do iseg = 1,size(seginds)
                    select case(seginds(iseg))
                        case(MIC_SEG)
                            call projects(iproj)%read_segment('mic', projfile_abs)
                            nmics(iproj) = projects(iproj)%os_mic%get_noris()
                        case(STK_SEG)
                            call projects(iproj)%read_segment('stk', projfile_abs)
                            call absolutize_project_stack_paths(projects(iproj), projfile_abs)
                            nstks(iproj) = projects(iproj)%os_stk%get_noris()
                        case(PTCL2D_SEG)
                            call projects(iproj)%read_segment('ptcl2D', projfile_abs)
                            nptcl2Ds(iproj) = projects(iproj)%os_ptcl2D%get_noris()
                        case(PTCL3D_SEG)
                            call projects(iproj)%read_segment('ptcl3D', projfile_abs)
                            nptcl3Ds(iproj) = projects(iproj)%os_ptcl3D%get_noris()
                        case(CLS2D_SEG, CLS3D_SEG, OUT_SEG)
                            ! Analysis products are intentionally dropped by merge_projects.
                        case(OPTICS_SEG)
                            call projects(iproj)%read_segment('optics', projfile_abs)
                            noptics(iproj) = projects(iproj)%os_optics%get_noris()
                        case(PROJINFO_SEG)
                            call projects(iproj)%read_segment('projinfo', projfile_abs)
                        case(JOBPROC_SEG)
                            call projects(iproj)%read_segment('jobproc', projfile_abs)
                        case(COMPENV_SEG)
                            call projects(iproj)%read_segment('compenv', projfile_abs)
                    end select
                enddo
                deallocate(seginds, seginfos)
            endif
        enddo
        call validate_field_presence(nmics,   'mic')
        call validate_field_presence(nstks,   'stk')
        call validate_field_presence(nptcl2Ds,'ptcl2D')
        call validate_field_presence(nptcl3Ds,'ptcl3D')
        call validate_field_presence(noptics, 'optics')
        l_has_mics    = nmics(1)    > 0
        l_has_stks    = nstks(1)    > 0
        l_has_ptcl2D  = nptcl2Ds(1) > 0
        l_has_ptcl3D  = nptcl3Ds(1) > 0
        l_has_optics  = noptics(1)  > 0
        l_has_any_data = l_has_mics .or. l_has_stks .or. l_has_ptcl2D .or. l_has_ptcl3D .or. l_has_optics
        if( .not.l_has_any_data )then
            do iproj = 1,nprojs
                write(logfhandle,*) 'input project ', iproj, ': ', project_fnames(iproj)%to_char()
                write(logfhandle,*) 'mic/stk/ptcl2D/ptcl3D/optics counts: ', &
                    nmics(iproj), nstks(iproj), nptcl2Ds(iproj), nptcl3Ds(iproj), noptics(iproj)
            enddo
            THROW_HARD('input projects contain no mergeable data fields beyond project metadata')
        endif
        do iproj = 1,nprojs
            call validate_source_project(projects(iproj), iproj)
        enddo
        if( l_has_stks ) call validate_stack_dimensions
        if( l_has_mics .and. .not.l_has_stks ) call validate_mic_sampling
        call merged_proj%kill
        call merged_proj%projinfo%copy(projects(1)%projinfo)
        call merged_proj%jobproc%copy(projects(1)%jobproc)
        call merged_proj%compenv%copy(projects(1)%compenv)
        call merged_proj%update_projinfo(projfile_out)
        if( l_has_mics   ) call merged_proj%os_mic%new(   sum(nmics),    is_ptcl=.false.)
        if( l_has_stks   ) call merged_proj%os_stk%new(   sum(nstks),    is_ptcl=.false.)
        if( l_has_ptcl2D ) call merged_proj%os_ptcl2D%new(sum(nptcl2Ds), is_ptcl=.true.)
        if( l_has_ptcl3D ) call merged_proj%os_ptcl3D%new(sum(nptcl3Ds), is_ptcl=.true.)
        if( l_has_optics ) call merged_proj%os_optics%new(sum(noptics),  is_ptcl=.false.)
        allocate(mic_offsets(nprojs), stk_offsets(nprojs), ptcl2D_offsets(nprojs), ptcl3D_offsets(nprojs))
        allocate(opt_offsets(nprojs))
        allocate(src_ogid_maxs(nprojs), ogid_offsets(nprojs))
        call make_prefix_offsets(nmics,    mic_offsets)
        call make_prefix_offsets(nstks,    stk_offsets)
        call make_prefix_offsets(nptcl2Ds, ptcl2D_offsets)
        call make_prefix_offsets(nptcl3Ds, ptcl3D_offsets)
        call make_prefix_offsets(noptics,  opt_offsets)
        do iproj = 1,nprojs
            src_ogid_maxs(iproj) = source_max_ogid(projects(iproj))
        enddo
        ogid_glob_max = 0
        do iproj = 1,nprojs
            if( src_ogid_maxs(iproj) > 0 )then
                ogid_offsets(iproj) = ogid_glob_max
                ogid_glob_max       = ogid_glob_max + src_ogid_maxs(iproj)
            else
                ogid_offsets(iproj) = 0
            endif
        enddo

        do iproj = 1,nprojs
            stk_offset    = stk_offsets(iproj)
            ptcl2D_offset = ptcl2D_offsets(iproj)
            ptcl3D_offset = ptcl3D_offsets(iproj)
            ogid_offset   = ogid_offsets(iproj)
            if( l_has_ptcl2D )then
                range_offset = ptcl2D_offset
            else if( l_has_ptcl3D )then
                range_offset = ptcl3D_offset
            else
                range_offset = 0
            endif

            !$omp parallel do default(shared) private(imic, imic_glob) schedule(static) if(nmics(iproj) > 1000)
            do imic = 1,nmics(iproj)
                imic_glob = mic_offsets(iproj) + imic
                call merged_proj%os_mic%transfer_ori(imic_glob, projects(iproj)%os_mic, imic)
                call remap_row_ogid(merged_proj%os_mic, imic_glob, ogid_offset)
            enddo
            !$omp end parallel do

            !$omp parallel do default(shared) private(istk, istk_glob, fromp, top) schedule(static) if(nstks(iproj) > 1000)
            do istk = 1,nstks(iproj)
                istk_glob = stk_offset + istk
                call merged_proj%os_stk%transfer_ori(istk_glob, projects(iproj)%os_stk, istk)
                if( (l_has_ptcl2D .or. l_has_ptcl3D) .and. &
                    merged_proj%os_stk%isthere(istk_glob, 'fromp') .and. &
                    merged_proj%os_stk%isthere(istk_glob, 'top') )then
                    fromp = projects(iproj)%os_stk%get_fromp(istk)
                    top   = projects(iproj)%os_stk%get_top(istk)
                    call merged_proj%os_stk%set(istk_glob, 'fromp', fromp + range_offset)
                    call merged_proj%os_stk%set(istk_glob, 'top',   top   + range_offset)
                endif
                call remap_row_ogid(merged_proj%os_stk, istk_glob, ogid_offset)
            enddo
            !$omp end parallel do

            !$omp parallel do default(shared) private(iptcl, iptcl2D_glob) schedule(static) if(nptcl2Ds(iproj) > 10000)
            do iptcl = 1,nptcl2Ds(iproj)
                iptcl2D_glob = ptcl2D_offset + iptcl
                call copy_particle_row(projects(iproj)%os_stk, projects(iproj)%os_ptcl2D, merged_proj%os_ptcl2D, iptcl, &
                    iptcl2D_glob, stk_offset, 0, ogid_offset, l_has_stks, .false.)
                call merged_proj%os_ptcl2D%delete_2Dclustering(iptcl2D_glob)
            enddo
            !$omp end parallel do
            !$omp parallel do default(shared) private(iptcl, iptcl3D_glob) schedule(static) if(nptcl3Ds(iproj) > 10000)
            do iptcl = 1,nptcl3Ds(iproj)
                iptcl3D_glob = ptcl3D_offset + iptcl
                call copy_particle_row(projects(iproj)%os_stk, projects(iproj)%os_ptcl3D, merged_proj%os_ptcl3D, iptcl, &
                    iptcl3D_glob, stk_offset, 0, ogid_offset, l_has_stks, .false.)
            enddo
            !$omp end parallel do

            !$omp parallel do default(shared) private(iopt, iopt_glob, ogid) schedule(static) if(noptics(iproj) > 1000)
            do iopt = 1,noptics(iproj)
                iopt_glob = opt_offsets(iproj) + iopt
                call merged_proj%os_optics%transfer_ori(iopt_glob, projects(iproj)%os_optics, iopt)
                call remap_row_ogid(merged_proj%os_optics, iopt_glob, ogid_offset)
                if( merged_proj%os_optics%isthere(iopt_glob, 'ogid') )then
                    ogid = merged_proj%os_optics%get_int(iopt_glob, 'ogid')
                    if( ogid > 0 ) call merged_proj%os_optics%set(iopt_glob, 'ogname', 'opticsgroup'//int2str(ogid))
                endif
            enddo
            !$omp end parallel do

            call projects(iproj)%kill
        enddo

        if( l_has_ptcl2D .and. (.not.l_has_ptcl3D) )then
            call merged_proj%os_ptcl3D%new(sum(nptcl2Ds), is_ptcl=.true.)
            !$omp parallel do default(shared) private(iptcl) schedule(static) if(sum(nptcl2Ds) > 10000)
            do iptcl = 1,sum(nptcl2Ds)
                call merged_proj%os_ptcl3D%transfer_ori(iptcl, merged_proj%os_ptcl2D, iptcl)
                call merged_proj%os_ptcl3D%delete_2Dclustering(iptcl)
            enddo
            !$omp end parallel do
        endif
        call validate_source_project(merged_proj, 0)
        if( l_write_proj ) call merged_proj%write(projfile_out)

        contains

            logical function data_segments_present( segments )
                integer, allocatable, intent(in) :: segments(:)
                if( .not.allocated(segments) )then
                    data_segments_present = .false.
                else
                    data_segments_present = any(segments == MIC_SEG)    .or. any(segments == STK_SEG)   .or. &
                                            any(segments == PTCL2D_SEG) .or. any(segments == PTCL3D_SEG).or. &
                                            any(segments == OPTICS_SEG)
                endif
            end function data_segments_present

            subroutine make_prefix_offsets( counts, offsets )
                integer, intent(in)  :: counts(:)
                integer, intent(out) :: offsets(:)
                integer :: ip
                offsets(1) = 0
                do ip = 2,size(counts)
                    offsets(ip) = offsets(ip - 1) + counts(ip - 1)
                enddo
            end subroutine make_prefix_offsets

            subroutine validate_field_presence( counts, segment )
                integer,          intent(in) :: counts(:)
                character(len=*), intent(in) :: segment
                integer :: ip
                if( any(counts > 0) .and. any(counts == 0) )then
                    write(logfhandle,*) 'project field mismatch for segment: ', trim(segment)
                    do ip = 1,size(counts)
                        write(logfhandle,*) 'input project, count: ', ip, counts(ip)
                        write(logfhandle,*) project_fnames(ip)%to_char()
                    enddo
                    THROW_HARD('project field mismatch: '//trim(segment)//' is populated in only some input projects')
                endif
            end subroutine validate_field_presence

            subroutine validate_source_project( proj, iproj )
                class(sp_project), intent(inout) :: proj
                integer,           intent(in)    :: iproj
                integer :: imic, istk, iopt, nrange
                logical :: lctf
                if( l_has_mics )then
                    do imic = 1,proj%os_mic%get_noris()
                        if( proj%os_mic%isthere(imic, 'ctf') )then
                            lctf = ctf_enabled_for_row(proj%os_mic, imic, iproj, 'os_mic')
                            if( lctf )then
                                call require_row_field(proj%os_mic, imic, 'smpd',       'os_mic', iproj)
                                call require_row_field(proj%os_mic, imic, 'kv',         'os_mic', iproj)
                                call require_row_field(proj%os_mic, imic, 'cs',         'os_mic', iproj)
                                call require_row_field(proj%os_mic, imic, 'fraca',      'os_mic', iproj)
                                call require_row_field(proj%os_mic, imic, 'phaseplate', 'os_mic', iproj)
                            endif
                        endif
                    enddo
                endif
                if( l_has_stks )then
                    nrange = count_stack_particles(proj, iproj, l_has_ptcl2D .or. l_has_ptcl3D)
                    if( l_has_ptcl2D .and. nrange > 0 .and. nrange /= proj%os_ptcl2D%get_noris() )then
                        THROW_HARD('stack particle ranges do not match os_ptcl2D size in input project '//int2str(iproj))
                    endif
                    if( l_has_ptcl3D .and. nrange > 0 .and. nrange /= proj%os_ptcl3D%get_noris() )then
                        THROW_HARD('stack particle ranges do not match os_ptcl3D size in input project '//int2str(iproj))
                    endif
                    do istk = 1,proj%os_stk%get_noris()
                        if( proj%os_stk%isthere(istk, 'ctf') )then
                            lctf = ctf_enabled_for_row(proj%os_stk, istk, iproj, 'os_stk')
                            if( lctf )then
                                call require_row_field(proj%os_stk, istk, 'smpd',       'os_stk', iproj)
                                call require_row_field(proj%os_stk, istk, 'kv',         'os_stk', iproj)
                                call require_row_field(proj%os_stk, istk, 'cs',         'os_stk', iproj)
                                call require_row_field(proj%os_stk, istk, 'fraca',      'os_stk', iproj)
                                call require_row_field(proj%os_stk, istk, 'phaseplate', 'os_stk', iproj)
                            endif
                        endif
                    enddo
                endif
                if( l_has_ptcl2D )then
                    call validate_particle_field(proj%os_ptcl2D, proj, iproj, 'os_ptcl2D')
                    if( l_has_stks ) call validate_stack_particle_ranges(proj, proj%os_ptcl2D, iproj, 'os_ptcl2D')
                endif
                if( l_has_ptcl3D )then
                    call validate_particle_field(proj%os_ptcl3D, proj, iproj, 'os_ptcl3D')
                    if( l_has_stks ) call validate_stack_particle_ranges(proj, proj%os_ptcl3D, iproj, 'os_ptcl3D')
                endif
                if( l_has_optics )then
                    do iopt = 1,proj%os_optics%get_noris()
                        call require_row_field(proj%os_optics, iopt, 'ogid', 'os_optics', iproj)
                    enddo
                endif
            end subroutine validate_source_project

            subroutine validate_particle_field( os, proj, iproj, segment )
                class(oris),       intent(in)    :: os
                class(sp_project), intent(inout) :: proj
                integer,           intent(in)    :: iproj
                character(len=*),  intent(in)    :: segment
                integer, parameter :: NO_BAD = huge(1)
                integer :: iptcl, stkind, nptcls, nstks_here
                integer :: bad_stkind_missing, bad_stkind_range
                integer :: bad_dfx, bad_angast, bad_phshift
                logical, allocatable :: stack_ctf(:), stack_phaseplate(:)
                if( .not.l_has_stks ) return
                nptcls = os%get_noris()
                nstks_here = proj%os_stk%get_noris()
                allocate(stack_ctf(nstks_here), stack_phaseplate(nstks_here))
                stack_ctf        = .false.
                stack_phaseplate = .false.
                do stkind = 1,nstks_here
                    if( proj%os_stk%isthere(stkind, 'ctf') )then
                        stack_ctf(stkind) = ctf_enabled_for_row(proj%os_stk, stkind, iproj, 'os_stk')
                        if( stack_ctf(stkind) )then
                            stack_phaseplate(stkind) = phaseplate_enabled_for_row(proj%os_stk, stkind, iproj, 'os_stk')
                        endif
                    endif
                enddo
                bad_stkind_missing = NO_BAD
                bad_stkind_range   = NO_BAD
                bad_dfx            = NO_BAD
                bad_angast         = NO_BAD
                bad_phshift        = NO_BAD
                !$omp parallel do default(shared) private(iptcl, stkind) schedule(static) &
                !$omp& reduction(min:bad_stkind_missing,bad_stkind_range,bad_dfx,bad_angast,bad_phshift) &
                !$omp& if(nptcls > 10000)
                do iptcl = 1,nptcls
                    if( .not.os%isthere(iptcl, 'stkind') )then
                        bad_stkind_missing = min(bad_stkind_missing, iptcl)
                    else
                        stkind = os%get_int(iptcl, 'stkind')
                        if( stkind < 1 .or. stkind > nstks_here )then
                            bad_stkind_range = min(bad_stkind_range, iptcl)
                        else
                            if( stack_ctf(stkind) )then
                                if( .not.os%isthere(iptcl, 'dfx') ) bad_dfx = min(bad_dfx, iptcl)
                                if( os%isthere(iptcl, 'dfy') .and. (.not.os%isthere(iptcl, 'angast')) )then
                                    bad_angast = min(bad_angast, iptcl)
                                endif
                                if( stack_phaseplate(stkind) .and. (.not.os%isthere(iptcl, 'phshift')) )then
                                    bad_phshift = min(bad_phshift, iptcl)
                                endif
                            endif
                        endif
                    endif
                enddo
                !$omp end parallel do
                if( bad_stkind_missing /= NO_BAD )then
                    call require_row_field(os, bad_stkind_missing, 'stkind', trim(segment), iproj)
                endif
                if( bad_stkind_range /= NO_BAD )then
                    write(logfhandle,*) 'segment, row, input project: ', trim(segment), bad_stkind_range, iproj
                    THROW_HARD('particle stkind out of range in input project '//int2str(iproj))
                endif
                if( bad_dfx /= NO_BAD ) call require_row_field(os, bad_dfx, 'dfx', trim(segment), iproj)
                if( bad_angast /= NO_BAD ) call require_row_field(os, bad_angast, 'angast', trim(segment), iproj)
                if( bad_phshift /= NO_BAD ) call require_row_field(os, bad_phshift, 'phshift', trim(segment), iproj)
            end subroutine validate_particle_field

            subroutine validate_stack_particle_ranges( proj, os, iproj, segment )
                class(sp_project), intent(inout) :: proj
                class(oris),       intent(in)    :: os
                integer,           intent(in)    :: iproj
                character(len=*),  intent(in)    :: segment
                integer :: istk, iptcl, stkind, nptcls, nstks_here, fromp, top, nptcls_project
                integer :: expected_fromp, bad_iptcl
                nptcls       = os%get_noris()
                nstks_here   = proj%os_stk%get_noris()
                expected_fromp = 1
                do istk = 1,nstks_here
                    call require_row_field(proj%os_stk, istk, 'fromp', 'os_stk', iproj)
                    call require_row_field(proj%os_stk, istk, 'top',   'os_stk', iproj)
                    fromp = proj%os_stk%get_fromp(istk)
                    top   = proj%os_stk%get_top(istk)
                    if( proj%os_stk%isthere(istk, 'nptcls') )then
                        nptcls_project = proj%os_stk%get_int(istk, 'nptcls')
                    else
                        nptcls_project = top - fromp + 1
                    endif
                    if( nptcls_project < 0 )then
                        THROW_HARD('negative stack nptcls in input project '//int2str(iproj))
                    endif
                    if( fromp /= expected_fromp )then
                        write(logfhandle,*) 'segment, stack, input project: ', trim(segment), istk, iproj
                        write(logfhandle,*) 'expected/fromp: ', expected_fromp, fromp
                        THROW_HARD('non-contiguous stack particle ranges in input project '//int2str(iproj))
                    endif
                    if( nptcls_project == 0 )then
                        if( top >= fromp )then
                            THROW_HARD('zero-particle stack has non-empty range in input project '//int2str(iproj))
                        endif
                    else
                        if( top < fromp .or. top > nptcls )then
                            write(logfhandle,*) 'segment, stack, input project: ', trim(segment), istk, iproj
                            write(logfhandle,*) 'fromp/top/nptcls: ', fromp, top, nptcls
                            THROW_HARD('stack particle range out of bounds in input project '//int2str(iproj))
                        endif
                        if( top - fromp + 1 /= nptcls_project )then
                            write(logfhandle,*) 'segment, stack, input project: ', trim(segment), istk, iproj
                            write(logfhandle,*) 'fromp/top/nptcls: ', fromp, top, nptcls_project
                            THROW_HARD('stack nptcls inconsistent with fromp/top in input project '//int2str(iproj))
                        endif
                        bad_iptcl = 0
                        do iptcl = fromp,top
                            call require_row_field(os, iptcl, 'stkind', trim(segment), iproj)
                            stkind = os%get_int(iptcl, 'stkind')
                            if( stkind /= istk )then
                                bad_iptcl = iptcl
                                exit
                            endif
                        enddo
                        if( bad_iptcl > 0 )then
                            write(logfhandle,*) 'segment, particle, input project: ', trim(segment), bad_iptcl, iproj
                            write(logfhandle,*) 'expected/stkind: ', istk, os%get_int(bad_iptcl, 'stkind')
                            THROW_HARD('particle stkind does not match stack range in input project '//int2str(iproj))
                        endif
                    endif
                    expected_fromp = expected_fromp + nptcls_project
                enddo
                if( expected_fromp /= nptcls + 1 )then
                    write(logfhandle,*) 'segment, input project: ', trim(segment), iproj
                    write(logfhandle,*) 'stack ranges cover/project particles: ', expected_fromp - 1, nptcls
                    THROW_HARD('stack ranges do not cover all particles in input project '//int2str(iproj))
                endif
            end subroutine validate_stack_particle_ranges

            ! count project particle rows regardless of selection
            integer function count_stack_particles( proj, iproj, require_ranges )
                class(sp_project), intent(in) :: proj
                integer,           intent(in) :: iproj
                logical,           intent(in) :: require_ranges
                integer :: istk
                count_stack_particles = 0
                do istk = 1,proj%os_stk%get_noris()
                    if( proj%os_stk%isthere(istk, 'nptcls') )then
                        count_stack_particles = count_stack_particles + proj%os_stk%get_int(istk, 'nptcls')
                    else if( proj%os_stk%isthere(istk, 'fromp') .and. proj%os_stk%isthere(istk, 'top') )then
                        count_stack_particles = count_stack_particles + &
                            proj%os_stk%get_top(istk) - proj%os_stk%get_fromp(istk) + 1
                    else if( require_ranges )then
                        THROW_HARD('missing stack nptcls/fromp/top for input project '//int2str(iproj))
                    endif
                enddo
            end function count_stack_particles

            subroutine validate_stack_dimensions
                real    :: smpd_ref, smpd_here
                integer :: box_ref, box_here, ip, is
                logical :: check_box, check_smpd
                check_box  = projects(1)%os_stk%isthere(1, 'box')
                check_smpd = projects(1)%os_stk%isthere(1, 'smpd')
                if( check_box  ) box_ref  = projects(1)%os_stk%get_int(1, 'box')
                if( check_smpd ) smpd_ref = projects(1)%os_stk%get(1, 'smpd')
                do ip = 1,nprojs
                    do is = 1,projects(ip)%os_stk%get_noris()
                        if( check_box )then
                            call require_row_field(projects(ip)%os_stk, is, 'box', 'os_stk', ip)
                            box_here = projects(ip)%os_stk%get_int(is, 'box')
                            if( box_here /= box_ref ) THROW_HARD('project merge requires identical stack boxes')
                        endif
                        if( check_smpd )then
                            call require_row_field(projects(ip)%os_stk, is, 'smpd', 'os_stk', ip)
                            smpd_here = projects(ip)%os_stk%get(is, 'smpd')
                            if( abs(smpd_here - smpd_ref) > SMPD_TOL )then
                                THROW_HARD('project merge requires identical stack sampling distance')
                            endif
                        endif
                    enddo
                enddo
            end subroutine validate_stack_dimensions

            subroutine validate_mic_sampling
                real :: smpd_ref, smpd_here
                integer :: ip, im
                if( .not.projects(1)%os_mic%isthere(1, 'smpd') ) return
                smpd_ref = projects(1)%os_mic%get(1, 'smpd')
                do ip = 1,nprojs
                    do im = 1,projects(ip)%os_mic%get_noris()
                        call require_row_field(projects(ip)%os_mic, im, 'smpd', 'os_mic', ip)
                        smpd_here = projects(ip)%os_mic%get(im, 'smpd')
                        if( abs(smpd_here - smpd_ref) > SMPD_TOL )then
                            THROW_HARD('project merge requires identical micrograph sampling distance')
                        endif
                    enddo
                enddo
            end subroutine validate_mic_sampling

            subroutine copy_particle_row( os_stk_src, os_src, os_dst, i_src, i_dst, stk_off, cls_off, ogid_off, &
                remap_stk, remap_cls )
                class(oris), intent(in)    :: os_stk_src, os_src
                class(oris), intent(inout) :: os_dst
                integer,     intent(in)    :: i_src, i_dst, stk_off, cls_off, ogid_off
                logical,     intent(in)    :: remap_stk, remap_cls
                integer :: val
                call os_dst%transfer_ori(i_dst, os_src, i_src)
                if( remap_stk .and. os_dst%isthere(i_dst, 'stkind') )then
                    call os_dst%set(i_dst, 'indstk', resolved_particle_indstk(os_stk_src, os_src, i_src))
                    val = os_dst%get_int(i_dst, 'stkind')
                    if( val > 0 ) call os_dst%set_stkind(i_dst, val + stk_off)
                endif
                if( remap_cls .and. os_dst%isthere(i_dst, 'class') )then
                    val = os_dst%get_int(i_dst, 'class')
                    if( val > 0 ) call os_dst%set_class(i_dst, val + cls_off)
                endif
                call remap_row_ogid(os_dst, i_dst, ogid_off)
            end subroutine copy_particle_row

            integer function resolved_particle_indstk( os_stk_src, os_src, i_src )
                class(oris), intent(in) :: os_stk_src, os_src
                integer,     intent(in) :: i_src
                integer :: stkind, fromp, top, indstk, nptcls_stk
                if( .not.os_src%isthere(i_src, 'stkind') )then
                    write(logfhandle,*) 'particle row: ', i_src
                    THROW_HARD('missing particle stkind while resolving indstk during project merge')
                endif
                stkind = os_src%get_int(i_src, 'stkind')
                if( stkind < 1 .or. stkind > os_stk_src%get_noris() )then
                    write(logfhandle,*) 'particle row/stkind/nstks: ', i_src, stkind, os_stk_src%get_noris()
                    THROW_HARD('particle stkind out of range while resolving indstk during project merge')
                endif
                if( .not.(os_stk_src%isthere(stkind, 'fromp') .and. os_stk_src%isthere(stkind, 'top')) )then
                    write(logfhandle,*) 'particle row/stkind: ', i_src, stkind
                    THROW_HARD('missing stack fromp/top while resolving indstk during project merge')
                endif
                fromp = os_stk_src%get_fromp(stkind)
                top   = os_stk_src%get_top(stkind)
                if( i_src < fromp .or. i_src > top )then
                    write(logfhandle,*) 'particle row/stkind/fromp/top: ', i_src, stkind, fromp, top
                    THROW_HARD('particle outside stack range while resolving indstk during project merge')
                endif
                resolved_particle_indstk = i_src - fromp + 1
                if( os_stk_src%isthere(stkind, 'nptcls_stk') )then
                    nptcls_stk = os_stk_src%get_int(stkind, 'nptcls_stk')
                    if( os_src%isthere(i_src, 'indstk') )then
                        indstk = os_src%get_int(i_src, 'indstk')
                        if( indstk > 0 .and. indstk <= nptcls_stk ) resolved_particle_indstk = indstk
                    endif
                endif
            end function resolved_particle_indstk

            subroutine require_row_field( os, irow, key, segment, iproj )
                class(oris),      intent(in) :: os
                integer,          intent(in) :: irow, iproj
                character(len=*), intent(in) :: key, segment
                if( .not. os%isthere(irow, key) )then
                    write(logfhandle,*) 'missing field: ', trim(key)
                    write(logfhandle,*) 'segment, row, input project: ', trim(segment), irow, iproj
                    THROW_HARD('missing required row field during project merge')
                endif
            end subroutine require_row_field

            logical function ctf_enabled_for_row( os, irow, iproj, segment )
                class(oris),       intent(in) :: os
                integer,           intent(in) :: irow, iproj
                character(len=*),  intent(in) :: segment
                character(len=STDLEN) :: ctf_flag
                call os%get_static(irow, 'ctf', ctf_flag)
                select case(trim(ctf_flag))
                    case('no')
                        ctf_enabled_for_row = .false.
                    case('yes','flip')
                        ctf_enabled_for_row = .true.
                    case DEFAULT
                        THROW_HARD('unsupported ctf flag in '//trim(segment)//' for input project '//int2str(iproj))
                end select
            end function ctf_enabled_for_row

            logical function phaseplate_enabled_for_row( os, irow, iproj, segment )
                class(oris),       intent(in) :: os
                integer,           intent(in) :: irow, iproj
                character(len=*),  intent(in) :: segment
                character(len=STDLEN) :: phaseplate
                call os%get_static(irow, 'phaseplate', phaseplate)
                select case(trim(phaseplate))
                    case('yes')
                        phaseplate_enabled_for_row = .true.
                    case('no')
                        phaseplate_enabled_for_row = .false.
                    case DEFAULT
                        THROW_HARD('unsupported phaseplate flag in '//trim(segment)//' for input project '//int2str(iproj))
                end select
            end function phaseplate_enabled_for_row

            integer function source_max_ogid( proj )
                class(sp_project), intent(inout) :: proj
                source_max_ogid = max_ogid_in_oris(proj%os_mic)
                source_max_ogid = max(source_max_ogid, max_ogid_in_oris(proj%os_stk))
                source_max_ogid = max(source_max_ogid, max_ogid_in_oris(proj%os_ptcl2D))
                source_max_ogid = max(source_max_ogid, max_ogid_in_oris(proj%os_ptcl3D))
                source_max_ogid = max(source_max_ogid, max_ogid_in_oris(proj%os_optics))
            end function source_max_ogid

            integer function max_ogid_in_oris( os ) result(max_ogid)
                class(oris), intent(in) :: os
                integer :: i, ogid, noris
                max_ogid = 0
                noris = os%get_noris()
                !$omp parallel do default(shared) private(i, ogid) reduction(max:max_ogid) schedule(static) if(noris > 10000)
                do i = 1,noris
                    if( os%isthere(i, 'ogid') )then
                        ogid = os%get_int(i, 'ogid')
                        if( ogid > max_ogid ) max_ogid = ogid
                    endif
                enddo
                !$omp end parallel do
            end function max_ogid_in_oris

            subroutine remap_row_ogid( os, irow, offset )
                class(oris), intent(inout) :: os
                integer,     intent(in)    :: irow, offset
                integer :: ogid
                if( os%isthere(irow, 'ogid') )then
                    ogid = os%get_int(irow, 'ogid')
                    if( ogid > 0 ) call os%set_ogid(irow, ogid + offset)
                endif
            end subroutine remap_row_ogid

    end subroutine merge_selected_project_files

    subroutine absolutize_project_stack_paths( proj, projfile )
        class(sp_project), intent(inout) :: proj
        class(string),     intent(in)    :: projfile
        type(string) :: projdir, fname, absfname
        integer :: istk, nstks
        projdir = get_fpath(projfile)
        nstks = proj%os_stk%get_noris()
        do istk = 1,nstks
            if( proj%os_stk%isthere(istk, 'stk') )then
                call proj%os_stk%getter(istk, 'stk', fname)
                if( fname%to_char([1,1]) /= '/' ) fname = filepath(projdir, fname)
                absfname = simple_abspath(fname, check_exists=.false.)
                call proj%os_stk%set(istk, 'stk', absfname)
            endif
            if( proj%os_stk%isthere(istk, 'boxfile') )then
                call proj%os_stk%getter(istk, 'boxfile', fname)
                if( fname%to_char([1,1]) /= '/' ) fname = filepath(projdir, fname)
                absfname = simple_abspath(fname, check_exists=.false.)
                call proj%os_stk%set(istk, 'boxfile', absfname)
            endif
        enddo
        call fname%kill
        call absfname%kill
        call projdir%kill
    end subroutine absolutize_project_stack_paths

    subroutine validate_and_repair_project_file( projfile_in, projfile_out )
        class(string), intent(in)  :: projfile_in
        type(string),  intent(out) :: projfile_out
        type(sp_project) :: proj
        integer, allocatable :: stack_counts(:)
        logical, allocatable :: trusted_nptcls_stk(:)
        integer :: nstks, nptcls_ref, nptcls2D, nptcls3D
        integer :: nwarns, nerrors, nrepairs
        if( fname2format(projfile_in) /= 'O' )then
            THROW_HARD('validate_projfile requires a SIMPLE project file (*.simple)')
        endif
        projfile_out = swap_suffix(projfile_in, string(METADATA_EXT), string('_validated')//METADATA_EXT)
        nwarns   = 0
        nerrors  = 0
        nrepairs = 0
        write(logfhandle,'(A)') '>>> VALIDATING PROJECT FILE: '//projfile_in%to_char()
        write(logfhandle,'(A)') '>>> VALIDATED OUTPUT FILE:  '//projfile_out%to_char()
        call proj%read(projfile_in)
        nstks    = proj%os_stk%get_noris()
        nptcls2D = proj%os_ptcl2D%get_noris()
        nptcls3D = proj%os_ptcl3D%get_noris()
        if( nstks == 0 )then
            call warn('project has no stack rows; stack-index policy checks skipped')
        endif
        if( nptcls2D == 0 .and. nptcls3D == 0 )then
            call warn('project has no particle rows; particle-index policy checks skipped')
        endif
        if( nstks > 0 )then
            allocate(stack_counts(nstks), trusted_nptcls_stk(nstks))
            stack_counts = 0
            trusted_nptcls_stk = .false.
            if( nptcls2D > 0 )then
                nptcls_ref = nptcls2D
                call repair_stack_ranges(proj%os_ptcl2D, 'ptcl2D', nptcls_ref)
            else if( nptcls3D > 0 )then
                nptcls_ref = nptcls3D
                call repair_stack_ranges(proj%os_ptcl3D, 'ptcl3D', nptcls_ref)
            else
                nptcls_ref = 0
                call repair_stack_counts_without_particles
            endif
            call repair_stack_nptcls_stk
            if( nptcls2D > 0 ) call repair_particle_segment(proj%os_ptcl2D, 'ptcl2D')
            if( nptcls3D > 0 ) call repair_particle_segment(proj%os_ptcl3D, 'ptcl3D')
            if( nptcls2D > 0 .and. nptcls3D > 0 .and. nptcls2D /= nptcls3D )then
                call err('ptcl2D/ptcl3D row-count mismatch: '//int2str(nptcls2D)//' / '//int2str(nptcls3D))
            endif
        endif
        call proj%update_projinfo(projfile_out)
        call proj%write(projfile_out)
        write(logfhandle,'(A)') '>>> VALIDATE_PROJFILE SUMMARY'
        write(logfhandle,'(A,I0)') '    repairs : ', nrepairs
        write(logfhandle,'(A,I0)') '    warnings: ', nwarns
        write(logfhandle,'(A,I0)') '    errors  : ', nerrors
        write(logfhandle,'(A)') '>>> WROTE VALIDATED PROJECT FILE: '//projfile_out%to_char()
        call proj%kill
        if( allocated(stack_counts) ) deallocate(stack_counts)
        if( allocated(trusted_nptcls_stk) ) deallocate(trusted_nptcls_stk)

        contains

            subroutine warn( msg )
                character(len=*), intent(in) :: msg
                nwarns = nwarns + 1
                write(logfhandle,'(A)') 'WARNING validate_projfile: '//trim(msg)
            end subroutine warn

            subroutine err( msg )
                character(len=*), intent(in) :: msg
                nerrors = nerrors + 1
                write(logfhandle,'(A)') 'ERROR validate_projfile: '//trim(msg)
            end subroutine err

            subroutine repair( msg )
                character(len=*), intent(in) :: msg
                nrepairs = nrepairs + 1
                write(logfhandle,'(A)') 'REPAIR validate_projfile: '//trim(msg)
            end subroutine repair

            subroutine repair_stack_counts_without_particles
                integer :: istk, count_here
                do istk = 1,nstks
                    count_here = 0
                    if( proj%os_stk%isthere(istk, 'nptcls') )then
                        count_here = max(0, proj%os_stk%get_int(istk, 'nptcls'))
                    else if( proj%os_stk%isthere(istk, 'fromp') .and. proj%os_stk%isthere(istk, 'top') )then
                        count_here = max(0, proj%os_stk%get_top(istk) - proj%os_stk%get_fromp(istk) + 1)
                    endif
                    stack_counts(istk) = count_here
                enddo
            end subroutine repair_stack_counts_without_particles

            subroutine repair_stack_ranges( os, segment, nptcls )
                class(oris),      intent(inout) :: os
                character(len=*), intent(in)    :: segment
                integer,          intent(in)    :: nptcls
                integer, allocatable :: stkind_counts(:)
                integer :: istk, iptcl, stkind, fromp, top, nptcls_stk, total_counts, expected_fromp
                integer :: count_from_range, count_from_nptcls, excess, trim_now
                allocate(stkind_counts(nstks))
                stkind_counts = 0
                do iptcl = 1,nptcls
                    if( os%isthere(iptcl, 'stkind') )then
                        stkind = os%get_int(iptcl, 'stkind')
                        if( stkind >= 1 .and. stkind <= nstks )then
                            stkind_counts(stkind) = stkind_counts(stkind) + 1
                        else
                            call err(trim(segment)//' row '//int2str(iptcl)//' has out-of-range stkind '//int2str(stkind))
                        endif
                    else
                        call err(trim(segment)//' row '//int2str(iptcl)//' is missing stkind')
                    endif
                enddo
                stack_counts = 0
                do istk = 1,nstks
                    count_from_range = -1
                    count_from_nptcls = -1
                    if( proj%os_stk%isthere(istk, 'fromp') .and. proj%os_stk%isthere(istk, 'top') )then
                        fromp = proj%os_stk%get_fromp(istk)
                        top   = proj%os_stk%get_top(istk)
                        if( fromp >= 1 .and. top >= fromp .and. top <= nptcls )then
                            count_from_range = top - fromp + 1
                        else if( top < fromp )then
                            count_from_range = 0
                            call warn('stack '//int2str(istk)//' has empty or reversed fromp/top range')
                        else
                            call err('stack '//int2str(istk)//' has out-of-bounds fromp/top range')
                        endif
                    endif
                    if( proj%os_stk%isthere(istk, 'nptcls') )then
                        nptcls_stk = proj%os_stk%get_int(istk, 'nptcls')
                        if( nptcls_stk >= 0 )then
                            count_from_nptcls = nptcls_stk
                        else
                            call err('stack '//int2str(istk)//' has negative project nptcls')
                        endif
                    endif
                    if( count_from_range >= 0 )then
                        stack_counts(istk) = count_from_range
                    else if( count_from_nptcls >= 0 )then
                        stack_counts(istk) = count_from_nptcls
                    else
                        stack_counts(istk) = stkind_counts(istk)
                        call warn('stack '//int2str(istk)//' missing fromp/top/nptcls; using particle stkind count')
                    endif
                    if( count_from_range >= 0 .and. count_from_nptcls >= 0 )then
                        if( count_from_range /= count_from_nptcls )then
                            call warn('stack '//int2str(istk)//' nptcls inconsistent with fromp/top; using fromp/top')
                        endif
                    endif
                enddo
                total_counts = sum(stack_counts)
                if( total_counts /= nptcls )then
                    if( sum(stkind_counts) == nptcls )then
                        stack_counts = stkind_counts
                        call warn('stack ranges do not cover particles; using stkind counts')
                    else
                        call err('stack ranges/counts cover '//int2str(total_counts)//' project rows, expected '//int2str(nptcls))
                        if( nstks > 0 )then
                            if( total_counts < nptcls )then
                                stack_counts(nstks) = stack_counts(nstks) + (nptcls - total_counts)
                                call repair('expanded final stack project range to cover all particles')
                            else
                                excess = total_counts - nptcls
                                do istk = nstks,1,-1
                                    if( excess <= 0 ) exit
                                    trim_now = min(stack_counts(istk), excess)
                                    stack_counts(istk) = stack_counts(istk) - trim_now
                                    excess = excess - trim_now
                                enddo
                                call repair('trimmed stack project ranges to fit particle count')
                            endif
                        endif
                    endif
                endif
                expected_fromp = 1
                do istk = 1,nstks
                    fromp = expected_fromp
                    top   = expected_fromp + stack_counts(istk) - 1
                    call set_stack_int_if_changed(istk, 'fromp', fromp)
                    call set_stack_int_if_changed(istk, 'top',   top)
                    call set_stack_int_if_changed(istk, 'nptcls', stack_counts(istk))
                    if( stack_counts(istk) > 0 )then
                        do iptcl = fromp,top
                            if( iptcl < 1 .or. iptcl > nptcls ) cycle
                            if( .not.os%isthere(iptcl, 'stkind') )then
                                call os%set_stkind(iptcl, istk)
                                call repair(trim(segment)//' row '//int2str(iptcl)//' stkind set from stack range')
                            else if( os%get_int(iptcl, 'stkind') /= istk )then
                                call os%set_stkind(iptcl, istk)
                                call repair(trim(segment)//' row '//int2str(iptcl)//' stkind set from stack range')
                            endif
                        enddo
                    endif
                    expected_fromp = top + 1
                enddo
                deallocate(stkind_counts)
            end subroutine repair_stack_ranges

            subroutine repair_stack_nptcls_stk
                integer :: istk, nptcls_stk
                do istk = 1,nstks
                    trusted_nptcls_stk(istk) = .false.
                    if( proj%os_stk%isthere(istk, 'nptcls_stk') )then
                        nptcls_stk = proj%os_stk%get_int(istk, 'nptcls_stk')
                        if( nptcls_stk >= stack_counts(istk) .and. nptcls_stk >= 0 )then
                            trusted_nptcls_stk(istk) = .true.
                        else
                            call warn('stack '//int2str(istk)//' has invalid nptcls_stk; using project range count')
                            call proj%os_stk%set(istk, 'nptcls_stk', stack_counts(istk))
                            call repair('stack '//int2str(istk)//' nptcls_stk set to '//int2str(stack_counts(istk)))
                        endif
                    else
                        call proj%os_stk%set(istk, 'nptcls_stk', stack_counts(istk))
                        call repair('stack '//int2str(istk)//' missing nptcls_stk; set to project range count')
                    endif
                enddo
            end subroutine repair_stack_nptcls_stk

            subroutine repair_particle_segment( os, segment )
                class(oris),      intent(inout) :: os
                character(len=*), intent(in)    :: segment
                integer :: iptcl, nptcls, stkind, fallback_indstk, indstk, nptcls_stk
                logical :: use_existing
                nptcls = os%get_noris()
                do iptcl = 1,nptcls
                    stkind = stkind_for_project_row(os, segment, iptcl)
                    if( stkind < 1 .or. stkind > nstks ) cycle
                    fallback_indstk = iptcl - proj%os_stk%get_fromp(stkind) + 1
                    if( fallback_indstk < 1 )then
                        call err(trim(segment)//' row '//int2str(iptcl)//' cannot be mapped into stack range')
                        cycle
                    endif
                    nptcls_stk = proj%os_stk%get_int(stkind, 'nptcls_stk')
                    if( fallback_indstk > nptcls_stk )then
                        call proj%os_stk%set(stkind, 'nptcls_stk', fallback_indstk)
                        nptcls_stk = fallback_indstk
                        trusted_nptcls_stk(stkind) = .false.
                        call repair('stack '//int2str(stkind)//' nptcls_stk expanded to cover fallback indstk')
                    endif
                    use_existing = .false.
                    if( trusted_nptcls_stk(stkind) .and. os%isthere(iptcl, 'indstk') )then
                        indstk = os%get_int(iptcl, 'indstk')
                        use_existing = indstk > 0 .and. indstk <= nptcls_stk
                    endif
                    if( .not.use_existing )then
                        if( os%isthere(iptcl, 'indstk') )then
                            indstk = os%get_int(iptcl, 'indstk')
                        else
                            indstk = 0
                        endif
                        if( indstk /= fallback_indstk )then
                            call os%set(iptcl, 'indstk', fallback_indstk)
                            call repair(trim(segment)//' row '//int2str(iptcl)//' indstk set from fromp/top fallback')
                        else if( .not.os%isthere(iptcl, 'indstk') )then
                            call os%set(iptcl, 'indstk', fallback_indstk)
                            call repair(trim(segment)//' row '//int2str(iptcl)//' missing indstk set from fromp/top fallback')
                        endif
                    endif
                enddo
            end subroutine repair_particle_segment

            integer function stkind_for_project_row( os, segment, iptcl ) result(stkind)
                class(oris),      intent(inout) :: os
                character(len=*), intent(in)    :: segment
                integer,          intent(in)    :: iptcl
                integer :: istk
                stkind = 0
                if( os%isthere(iptcl, 'stkind') ) stkind = os%get_int(iptcl, 'stkind')
                if( stkind >= 1 .and. stkind <= nstks ) return
                do istk = 1,nstks
                    if( iptcl >= proj%os_stk%get_fromp(istk) .and. iptcl <= proj%os_stk%get_top(istk) )then
                        stkind = istk
                        call os%set_stkind(iptcl, stkind)
                        call repair(trim(segment)//' row '//int2str(iptcl)//' stkind repaired from stack range')
                        return
                    endif
                enddo
                call err(trim(segment)//' row '//int2str(iptcl)//' has no valid stack range')
            end function stkind_for_project_row

            subroutine set_stack_int_if_changed( istk, key, val )
                integer,          intent(in) :: istk, val
                character(len=*), intent(in) :: key
                if( .not.proj%os_stk%isthere(istk, key) )then
                    call proj%os_stk%set(istk, key, val)
                    call repair('stack '//int2str(istk)//' '//trim(key)//' set to '//int2str(val))
                else if( proj%os_stk%get_int(istk, key) /= val )then
                    call proj%os_stk%set(istk, key, val)
                    call repair('stack '//int2str(istk)//' '//trim(key)//' set to '//int2str(val))
                endif
            end subroutine set_stack_int_if_changed

    end subroutine validate_and_repair_project_file

end module simple_projfile_utils
