!@descr: fills the GUI's project metadata from an in-memory project, writing the previews it shows
! build_project_metadata does for gui_communicator%add_metadata what the record cannot: it reads
! the sp_project segments chosen by oritype (mov, mic, ptcl, cls2D, cls3D; all of them for 'all'),
! writes the preview files (a movie thumbnail, a particle montage, per-stage copies of the 3D
! reprojection and heat-map JPEGs) and reads the box files, FSC curves and MRC headers.
! gui_metadata_project only records the result, so the metadata types need nothing from src/main.
! stage selects the 2D/3D stage slot (stage 1 starts a new run), 0 the final one. A reporting gap
! never fails the job: a section that cannot be filled is left out with a warning.
module simple_gui_project_builder
use simple_defs,                    only: logfhandle, GUI_PSPECSZ
use simple_defs_fname,              only: MRC_EXT, JPG_EXT, MOVTHUMB_FBODY, PPROC_SUFFIX, LP_SUFFIX, MIRR_SUFFIX
use simple_error,                   only: simple_exception
use simple_string,                  only: string
use simple_string_utils,            only: int2str, int2str_pad
use simple_fileio,                  only: swap_suffix, file_exists, file2rarr, add2fbody, get_fpath, simple_copy_file
use simple_syslib,                  only: del_file, simple_abspath, get_process_id
use simple_sp_project,              only: sp_project
use simple_image,                   only: image
use simple_imghead,                 only: get_mrc_minmax
use simple_nrtxtfile,               only: nrtxtfile
use simple_math,                    only: round2even
use simple_math_ft,                 only: get_resarr
use simple_estimate_ssnr,           only: get_resolution
use simple_motion_gain_helpers,     only: read_movies_and_sum_frames
use simple_procimgstk,              only: random_selection_from_imgfile, bp_imgfile
use simple_refine3D_fnames,         only: refine3D_oris_heatmap_fname, refine3D_reprojs_fname
use simple_oris_utils,              only: oridist_from_oris
use simple_gui_utils,               only: mrc2jpeg_tiled
use simple_gui_metadata_types,      only: GUI_METADATA_MICROGRAPH_TYPE, GUI_METADATA_CAVG2D_TYPE,&
                                         &GUI_METADATA_PTCL_TYPE, GUI_METADATA_VOL3D_TYPE
use simple_gui_metadata_micrograph, only: gui_metadata_micrograph, MAX_MIC_COORDINATES
use simple_gui_metadata_ptcl,       only: gui_metadata_ptcl
use simple_gui_metadata_cavg2D,     only: gui_metadata_cavg2D, sprite_sheet_pos
use simple_gui_metadata_vol3D,      only: gui_metadata_vol3D, MAX_FSC_VOL3D, ORIDIST_NBINS_X, ORIDIST_NBINS_Y
use simple_gui_metadata_project,    only: gui_metadata_project
implicit none

public :: build_project_metadata
private
#include "simple_local_flags.inc"

contains

    !> Updates @p meta with the sections of @p spproj that @p oritype selects, for @p stage (0 the
    !! final one); @p selection adds the selected counts.
    subroutine build_project_metadata( meta, spproj, oritype, stage, selection )
        type(gui_metadata_project), intent(inout) :: meta
        type(sp_project),           intent(inout) :: spproj
        character(len=*),           intent(in)    :: oritype
        integer,                    intent(in)    :: stage
        logical,                    intent(in)    :: selection
        type(string) :: projname, projfile
        logical      :: l_all
        if( stage < 0 ) THROW_HARD('stage must be >= 1, or 0 for the final output')
        call spproj%projinfo%getter(1, 'projname', projname)
        call spproj%projinfo%getter(1, 'projfile', projfile)
        call meta%set_project(projname, projfile)
        call meta%start_update(stage)
        l_all = oritype == 'all'
        if( l_all .or. oritype == 'mov' )                        call add_movies(meta, spproj)
        if( l_all .or. oritype == 'mic' .or. oritype == 'ptcl' ) call add_micrographs(meta, spproj, oritype == 'ptcl', selection)
        if( l_all .or. oritype == 'ptcl' )                       call add_particles(meta, spproj)
        if( l_all .or. oritype == 'cls2D' )                      call add_cavgs2D(meta, spproj, stage, selection)
        if( l_all .or. oritype == 'cls3D' )                      call add_vols3D(meta, spproj, stage)
    end subroutine build_project_metadata

    ! A thumbnail of the first movie: its frames summed and Fourier-cropped to GUI_PSPECSZ pixels
    ! across, written as movthumb<i>.jpg in the working directory
    subroutine add_movies( meta, spproj )
        type(gui_metadata_project), intent(inout) :: meta
        type(sp_project),           intent(inout) :: spproj
        integer, parameter :: N_MOV_THUMBS = 1
        type(gui_metadata_micrograph) :: movies(N_MOV_THUMBS)
        type(image)  :: movsum, movthumb
        type(string) :: movfname, thumbfname
        integer      :: i, nmics, nthumbs, n_movies_sum, n_frames_sum, ldim_mov(3), ldim_thumb(3)
        real         :: smpd, scale_thumb
        nmics = spproj%os_mic%get_noris()
        call meta%set_summary(nmics=nmics)
        if( nmics == 0 ) return
        nthumbs = 0
        do i = 1,nmics
            if( nthumbs >= N_MOV_THUMBS ) exit
            if( .not. spproj%os_mic%isthere(i, 'imgkind') ) cycle
            if( spproj%os_mic%get_str(i, 'imgkind') /= 'movie' ) cycle
            nthumbs  = nthumbs + 1
            smpd     = spproj%os_mic%get(i, 'smpd')
            movfname = spproj%os_mic%get_str(i, 'movie')
            call read_movies_and_sum_frames([movfname], smpd, movsum, n_movies_sum, n_frames_sum)
            ldim_mov      = movsum%get_ldim()
            scale_thumb   = real(GUI_PSPECSZ) / real(ldim_mov(1))
            ldim_thumb(1) = round2even(real(ldim_mov(1)) * scale_thumb)
            ldim_thumb(2) = round2even(real(ldim_mov(2)) * scale_thumb)
            ldim_thumb(3) = 1
            call movthumb%new(ldim_thumb, smpd)
            call movsum%fft()
            call movsum%clip(movthumb)
            call movthumb%ifft()
            thumbfname = MOVTHUMB_FBODY // int2str(i) // JPG_EXT
            call movthumb%write_jpg(thumbfname, norm=.true., quality=90)
            thumbfname = simple_abspath(thumbfname)
            call movies(nthumbs)%new(GUI_METADATA_MICROGRAPH_TYPE)
            call movies(nthumbs)%set(path=thumbfname, i_max=N_MOV_THUMBS, i=i)
            call movsum%kill()
            call movthumb%kill()
        end do
        call meta%set_movies(movies(1:nthumbs))
    end subroutine add_movies

    ! The first micrographs with their CTF fit and their picks (from their box files), when the
    ! micrographs have thumbnails; @p l_ptcl adds the stack and particle counts
    subroutine add_micrographs( meta, spproj, l_ptcl, l_selection )
        type(gui_metadata_project), intent(inout) :: meta
        type(sp_project),           intent(inout) :: spproj
        logical,                    intent(in)    :: l_ptcl, l_selection
        integer, parameter :: N_MICS_SHOWN = 50
        type(gui_metadata_micrograph), allocatable :: mics(:)
        real,                          allocatable :: boxdata(:)
        type(nrtxtfile) :: boxfile
        type(string)    :: boxpath
        integer :: i, j, x, y, nmics, nshown, nvalid, nrecs, nlines, xdim, ydim
        logical :: l_ctf
        nmics = spproj%os_mic%get_noris()
        call meta%set_summary(nmics=nmics, pspec_size=GUI_PSPECSZ)
        if( l_ptcl      ) call meta%set_summary(nstks=spproj%os_stk%get_noris(), nptcls=spproj%os_ptcl2D%get_noris())
        if( l_selection ) call meta%set_summary(nmics_selected=spproj%os_mic%count_state_gt_zero())
        if( .not. spproj%os_mic%isthere('thumb') ) return
        nshown = min(N_MICS_SHOWN, nmics)
        allocate(mics(nshown))
        xdim   = nint(spproj%os_mic%get(1, 'xdim'))
        ydim   = nint(spproj%os_mic%get(1, 'ydim'))
        l_ctf  = spproj%os_mic%isthere('ctfjpg')
        nvalid = 0
        do i = 1,nshown
            if( spproj%os_mic%get_state(i) == 0 ) cycle ! needs improvement to work with pagination
            nvalid = nvalid + 1
            call mics(nvalid)%new(GUI_METADATA_MICROGRAPH_TYPE)
            if( l_ctf )then
                call mics(nvalid)%set(path=spproj%os_mic%get_str(i, 'thumb'), dfx=spproj%os_mic%get(i, 'dfx'),&
                    &dfy=spproj%os_mic%get(i, 'dfy'), ctfres=spproj%os_mic%get(i, 'ctfres'),&
                    &ctfimg=spproj%os_mic%get_str(i, 'ctfjpg'), i_max=nshown, i=i)
            else
                call mics(nvalid)%set(path=spproj%os_mic%get_str(i, 'thumb'), dfx=spproj%os_mic%get(i, 'dfx'),&
                    &dfy=spproj%os_mic%get(i, 'dfy'), ctfres=spproj%os_mic%get(i, 'ctfres'), i_max=nshown, i=i)
            endif
            boxpath = spproj%os_mic%get_str(i, 'boxfile')
            if( boxpath%strlen() > 0 .and. file_exists(boxpath) )then
                call boxfile%new(boxpath, 1)
                nrecs  = boxfile%get_nrecs_per_line()
                nlines = boxfile%get_ndatalines()
                if( nrecs >= 4 )then
                    allocate(boxdata(nrecs))
                    ! at most the picks the micrograph's metadata holds
                    do j = 1,min(nlines, MAX_MIC_COORDINATES)
                        call boxfile%readNextDataLine(boxdata)
                        x = nint(boxdata(1) + boxdata(3)/2)
                        y = nint(boxdata(2) + boxdata(4)/2)
                        call mics(nvalid)%set_coordinate(j, x, y, xdim, ydim)
                    enddo
                    deallocate(boxdata)
                endif
                call boxfile%kill()
            endif
        enddo
        call meta%set_micrographs(mics(1:nvalid), xdim, ydim, spproj%os_mic%get(1, 'smpd'))
    end subroutine add_micrographs

    ! A JPEG montage of a random sample of the selected particles, with its low-pass twin, named
    ! after this process: concurrent batch jobs share the execution directory
    subroutine add_particles( meta, spproj )
        type(gui_metadata_project), intent(inout) :: meta
        type(sp_project),           intent(inout) :: spproj
        integer, parameter :: N_PTCLS_SAMPLE = 100
        type(gui_metadata_ptcl), allocatable :: ptcls(:)
        integer,                 allocatable :: ptcl_inds(:)
        type(string) :: stem, stk, jpg, lpstk, lpjpg
        integer :: i, nptcls_all, nptcls_valid, nsample, nshown, box, xtiles, ytiles, xtile, ytile
        real    :: smpd, df
        nptcls_all = spproj%os_ptcl2D%get_noris()
        if( nptcls_all == 0 ) return
        if( spproj%os_ptcl2D%isthere('state') )then
            nptcls_valid = spproj%os_ptcl2D%count_state_gt_zero()
        else
            nptcls_valid = nptcls_all
        endif
        nsample = min(N_PTCLS_SAMPLE, nptcls_valid)
        if( nsample == 0 ) return
        box   = nint(spproj%os_stk%get(1, 'box'))
        smpd  = spproj%os_stk%get(1, 'smpd')
        stem  = 'ptcls_sample_' // int2str(get_process_id())
        stk   = stem // MRC_EXT
        jpg   = stem // JPG_EXT
        lpstk = stem // '_lp' // MRC_EXT
        lpjpg = stem // '_lp' // JPG_EXT
        call random_selection_from_imgfile(spproj, stk, box, nsample, pinds=ptcl_inds)
        call mrc2jpeg_tiled(stk, jpg, ntiles=nshown)
        call bp_imgfile(stk, lpstk, smpd, 0., 10.)
        call mrc2jpeg_tiled(lpstk, lpjpg, ntiles=nshown)
        call del_file(lpstk)
        call del_file(stk)
        jpg    = simple_abspath(jpg)
        lpjpg  = simple_abspath(lpjpg)
        xtiles = floor(sqrt(real(nshown)))
        ytiles = ceiling(real(nshown) / real(xtiles))
        allocate(ptcls(nshown))
        do i = 1,nshown
            xtile = mod(i-1, xtiles)
            ytile = (i-1) / xtiles
            df    = (spproj%os_ptcl2D%get(ptcl_inds(i), 'dfx') + spproj%os_ptcl2D%get(ptcl_inds(i), 'dfy')) / 2.0
            call ptcls(i)%new(GUI_METADATA_PTCL_TYPE)
            call ptcls(i)%set(path=jpg, pathlp=lpjpg, i=i, i_max=nshown, df=df, box=box, idx=ptcl_inds(i),&
                &sprite=sprite_sheet_pos(x=xtile * (100.0 / max(1, xtiles - 1)), y=ytile * (100.0 / max(1, ytiles - 1)),&
                &h=100 * ytiles, w=100 * xtiles))
        enddo
        call meta%set_particles(jpg, ptcls)
    end subroutine add_particles

    ! The selected class averages of the project, as tiles of the stack's JPEG sprite sheet
    subroutine add_cavgs2D( meta, spproj, stage, l_selection )
        type(gui_metadata_project), intent(inout) :: meta
        type(sp_project),           intent(inout) :: spproj
        integer,                    intent(in)    :: stage
        logical,                    intent(in)    :: l_selection
        type(gui_metadata_cavg2D), allocatable :: cavgs(:)
        type(string) :: cavgsstk, cavgsjpg
        integer :: i, ncls, ncls_stk, nvalid, out_ind, xtiles, ytiles, xtile, ytile
        real    :: smpd, box, mskdiam
        ncls = spproj%os_cls2D%get_noris()
        call meta%set_summary(nstks=spproj%os_stk%get_noris(), nptcls=spproj%os_ptcl2D%get_noris(), ncls2D=ncls)
        if( l_selection ) call meta%set_summary(nptcls_selected=nint(spproj%os_cls2D%get_sum('pop')),&
            &ncls2D_selected=spproj%os_cls2D%count_state_gt_zero())
        if( ncls == 0 ) return
        box     = 0.
        out_ind = 0
        call spproj%get_cavgs_stk(cavgsstk, ncls_stk, smpd, fail=.false., out_ind=out_ind, box=box)
        if( ncls_stk /= ncls )then
            write(logfhandle,'(A,I0,A,I0,A)') '>>> WARNING: GUI project metadata: ', ncls_stk,&
                &' class averages for ', ncls, ' classes; the 2D classes are left out'
            return
        endif
        mskdiam = 0.
        if( out_ind > 0 )then
            if( spproj%os_out%isthere(out_ind, 'mskdiam') ) mskdiam = spproj%os_out%get(out_ind, 'mskdiam')
        endif
        cavgsjpg = swap_suffix(cavgsstk, JPG_EXT, MRC_EXT)
        xtiles   = floor(sqrt(real(ncls)))
        ytiles   = ceiling(real(ncls) / real(xtiles))
        allocate(cavgs(ncls))
        nvalid = 0
        do i = 1,ncls
            if( spproj%os_cls2D%get_state(i) == 0 ) cycle
            nvalid = nvalid + 1
            xtile  = mod(i-1, xtiles)
            ytile  = (i-1) / xtiles
            call cavgs(nvalid)%new(GUI_METADATA_CAVG2D_TYPE)
            call cavgs(nvalid)%set(path=cavgsjpg, mrcpath=cavgsstk, i=i, i_max=ncls, res=spproj%os_cls2D%get(i, 'res'),&
                &pop=spproj%os_cls2D%get_int(i, 'pop'), idx=i,&
                &sprite=sprite_sheet_pos(x=xtile * (100.0 / max(1, xtiles - 1)), y=ytile * (100.0 / max(1, ytiles - 1)),&
                &h=100 * ytiles, w=100 * xtiles))
        enddo
        call meta%set_cavgs2D(stage, cavgs(1:nvalid), dim_cavgs=nint(box), mskdiam=mskdiam, mskscale=box * smpd)
    end subroutine add_cavgs2D

    ! The state volumes of the project: their products, FSC, MRC header minimum and maximum,
    ! orthogonal reprojections and orientation histogram
    subroutine add_vols3D( meta, spproj, stage )
        type(gui_metadata_project), intent(inout) :: meta
        type(sp_project),           intent(inout) :: spproj
        integer,                    intent(in)    :: stage
        integer, parameter :: NTILES3D = 3
        type(gui_metadata_vol3D), allocatable :: vols(:)
        real,                     allocatable :: fsc_arr(:), res_arr(:)
        type(gui_metadata_cavg2D) :: tiles(NTILES3D)
        type(string) :: volpath, volpath_out, fsc_fname, pprocpath, lppath, pprocmirrpath, reprojpath, oridistpath
        real         :: smpd, res0143, res05, cfar, vmin, vmax
        integer      :: hist(ORIDIST_NBINS_X, ORIDIST_NBINS_Y)
        integer      :: nstates, istate, nvalid, box, pop, fsc_box, nfsc, itile, iptcl
        logical      :: l_final, l_have_fsc
        l_final = stage == 0
        nstates = spproj%os_ptcl3D%get_n('state')
        call meta%set_summary(nstates3D=nstates)
        if( nstates == 0 ) return
        allocate(vols(nstates))
        nvalid = 0
        do istate = 1,nstates
            if( .not. spproj%isthere_in_osout('vol', istate) ) cycle
            call spproj%get_vol('vol', istate, volpath, smpd, box)
            if( volpath%strlen() == 0 ) cycle
            nvalid = nvalid + 1
            pop    = spproj%os_ptcl3D%get_pop(istate, 'state')
            if( l_final )then
                volpath_out   = volpath
                lppath        = add2fbody(volpath, MRC_EXT, LP_SUFFIX)
                pprocpath     = add2fbody(volpath, MRC_EXT, PPROC_SUFFIX)
                if( .not. file_exists(pprocpath) ) pprocpath = string('')
                pprocmirrpath = string('')
                if( pprocpath%strlen() > 0 )then
                    pprocmirrpath = add2fbody(pprocpath, MRC_EXT, MIRR_SUFFIX)
                    if( .not. file_exists(pprocmirrpath) ) pprocmirrpath = string('')
                endif
            else
                ! non-final stages only ever get a per-stage lowpass snapshot,
                ! e.g. recvol_state01_stage03_lp.mrc (simple_solve3D_utils exec_refine3D);
                ! volpath/pprocpath/pprocmirrpath are final-only products, withheld here
                volpath_out   = string('')
                lppath        = add2fbody(volpath, MRC_EXT, '_stage'//int2str_pad(stage,2)//LP_SUFFIX)
                pprocpath     = string('')
                pprocmirrpath = string('')
            endif
            if( .not. file_exists(lppath) ) lppath = string('')
            ! the JPEGs the 3D writes beside the volume
            reprojpath  = stage_jpeg(get_fpath(volpath)//refine3D_reprojs_fname(istate), stage)
            oridistpath = stage_jpeg(get_fpath(volpath)//refine3D_oris_heatmap_fname(istate), stage)
            call vols(nvalid)%new(GUI_METADATA_VOL3D_TYPE)
            ! the FSC curve and the resolutions; cfar from the first particle of the state holding it
            l_have_fsc = .false.
            res0143    = 0.
            res05      = 0.
            cfar       = 0.
            if( spproj%isthere_in_osout('fsc', istate) )then
                call spproj%get_fsc(istate, fsc_fname, fsc_box)
                if( fsc_fname%strlen() > 0 )then
                    if( file_exists(fsc_fname) )then
                        fsc_arr = file2rarr(fsc_fname)
                        res_arr = get_resarr(fsc_box, smpd)
                        call get_resolution(fsc_arr, res_arr, res05, res0143)
                        do iptcl = 1,spproj%os_ptcl3D%get_noris()
                            if( spproj%os_ptcl3D%get_state(iptcl) /= istate ) cycle
                            if( .not. spproj%os_ptcl3D%isthere(iptcl, 'cfar') ) cycle
                            cfar = spproj%os_ptcl3D%get(iptcl, 'cfar')
                            exit
                        enddo
                        l_have_fsc = .true.
                    endif
                endif
            endif
            if( l_have_fsc )then
                call vols(nvalid)%set(reprojpath, volpath_out, pprocpath, lppath, pprocmirrpath, istate, box, smpd,&
                    &nvalid, nstates, res0143=res0143, res05=res05, cfar=cfar, pop=pop, oridistpath=oridistpath)
                nfsc = min(size(fsc_arr), size(res_arr), MAX_FSC_VOL3D)
                call vols(nvalid)%set_fsc(1.0 / res_arr(1:nfsc), fsc_arr(1:nfsc))
            else
                call vols(nvalid)%set(reprojpath, volpath_out, pprocpath, lppath, pprocmirrpath, istate, box, smpd,&
                    &nvalid, nstates, pop=pop, oridistpath=oridistpath)
            endif
            ! header minimum and maximum, read once here so the GUI need not open the volumes
            if( volpath_out%strlen() > 0 )then
                call get_mrc_minmax(volpath_out, vmin, vmax)
                call vols(nvalid)%set_minmax('volpath', vmin, vmax)
            endif
            if( pprocpath%strlen() > 0 )then
                call get_mrc_minmax(pprocpath, vmin, vmax)
                call vols(nvalid)%set_minmax('pprocpath', vmin, vmax)
            endif
            if( lppath%strlen() > 0 )then
                call get_mrc_minmax(lppath, vmin, vmax)
                call vols(nvalid)%set_minmax('lppath', vmin, vmax)
            endif
            if( pprocmirrpath%strlen() > 0 )then
                call get_mrc_minmax(pprocmirrpath, vmin, vmax)
                call vols(nvalid)%set_minmax('pprocmirrpath', vmin, vmax)
            endif
            ! the three orthogonal reprojections as tiles of their JPEG
            if( reprojpath%strlen() > 0 )then
                do itile = 1,NTILES3D
                    call tiles(itile)%new(GUI_METADATA_CAVG2D_TYPE)
                    call tiles(itile)%set(path=reprojpath, mrcpath=volpath_out, idx=istate,&
                        &sprite=sprite_sheet_pos(x=real(itile-1)*(100.0/real(NTILES3D-1)), y=0.0, h=100, w=100*NTILES3D),&
                        &i=itile, i_max=NTILES3D, pop=pop)
                enddo
                call vols(nvalid)%set_reprojtiles(tiles)
            endif
            ! the state's particle orientations as the azimuth/elevation histogram
            call oridist_from_oris(spproj%os_ptcl3D, istate, hist)
            call vols(nvalid)%set_oridist(hist)
        enddo
        call meta%set_vols3D(stage, vols(1:nvalid))
    end subroutine add_vols3D

    ! @p fname when it exists, '' otherwise. For a non-final @p stage, a copy named after the
    ! stage (e.g. ..._state01_stage03.jpg), since the next stage overwrites the file
    function stage_jpeg( fname, stage ) result( path )
        type(string), intent(in) :: fname
        integer,      intent(in) :: stage
        type(string) :: path
        path = string('')
        if( .not. file_exists(fname) ) return
        if( stage == 0 )then
            path = fname
        else
            path = add2fbody(fname, JPG_EXT, '_stage'//int2str_pad(stage,2))
            call simple_copy_file(fname, path)
        endif
    end function stage_jpeg

end module simple_gui_project_builder
