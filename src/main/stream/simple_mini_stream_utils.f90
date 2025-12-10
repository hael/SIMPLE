module simple_mini_stream_utils
include 'simple_lib.f08'
use simple_sp_project,  only: sp_project
use simple_picksegdiam, only: picksegdiam
use simple_parameters,  only: parameters
use simple_image_bin,   only: image_bin
use simple_image,       only: image
use simple_cmdline,     only: cmdline
use simple_gui_utils
use simple_micproc
use simple_nrtxtfile
implicit none
#include "simple_local_flags.inc"

contains

    subroutine segdiampick_preprocess( spproj, pcontrast, moldiam_max, outdir )
        class(sp_project), intent(inout) :: spproj
        character(len=*),  intent(in)    :: pcontrast
        class(string),     intent(in)    :: outdir
        real,              intent(in)    :: moldiam_max
        real,              parameter     :: SMPD_SHRINK1  = 4.0,  SIGMA_CRIT = 2., SIGMA_CRIT_MSK = 2.5
        integer,           parameter     :: BOXFAC = 3, NQ_DIAMS = 10
        type(picksegdiam)                :: picker
        type(image)                      :: mic_raw, mic_shrink
        type(string)  :: mic_name, mic_den_name, mic_topo_name, mic_bin_name, mic_diam_name
        integer :: ldim_raw(3), ldim(3), imic
        real    :: scale, smpd
        logical :: l_empty
        ! parse project
        smpd = spproj%get_smpd()
        scale = smpd / SMPD_SHRINK1
        ! read the first micrograph
        if( .not. spproj%os_mic%isthere(1, 'intg')) return
        call read_mic_subtr_backgr_shrink(spproj%os_mic%get_str(1,'intg'), smpd, scale, pcontrast, mic_raw, mic_shrink, l_empty)
        ! set logical dimensions
        ldim_raw = mic_raw%get_ldim()
        ldim     = mic_shrink%get_ldim()
        ! diameter estimation
        do imic = 1, spproj%os_mic%get_noris()
            if( .not. spproj%os_mic%isthere(imic, 'intg')) cycle
            if( spproj%os_mic%get(imic, 'state') < 1.0 )   cycle
            ! Segmentation & picking
            call spproj%os_mic%getter(imic, 'intg', mic_name)
            mic_den_name  = filepath(outdir, append2basename(mic_name, DEN_SUFFIX))
            mic_topo_name = filepath(outdir, append2basename(mic_name, TOPO_SUFFIX))
            mic_bin_name  = filepath(outdir, append2basename(mic_name, BIN_SUFFIX))
            mic_diam_name = swap_suffix(mic_name, TXT_EXT, '.mrc')
            mic_diam_name = filepath(outdir, append2basename(mic_diam_name, DIAMS_SUFFIX))
            call picker%pick(mic_name, smpd, moldiam_max, pcontrast, denfname=mic_den_name,&
                topofname=mic_topo_name, binfname=mic_bin_name, empty=l_empty )
            if(l_empty) then
                ! empty micrograph -> state=0
                call spproj%os_mic%set(imic, 'state', 0)
                cycle
            endif
            call spproj%os_mic%set(imic, 'mic_den',  mic_den_name)
            call spproj%os_mic%set(imic, 'mic_topo', mic_topo_name)
            call spproj%os_mic%set(imic, 'mic_bin',  mic_bin_name)   
            if( picker%get_nboxes() == 0 ) cycle
            call picker%write_diameters(mic_diam_name)
            call spproj%os_mic%set(imic, 'mic_diam', mic_diam_name)
        end do
        call picker%kill
        call mic_name%kill
    end subroutine segdiampick_preprocess

    subroutine segdiampick_mics( spproj, pcontrast, mic_to, moldiam_max, box_in_pix, mskdiam )
        class(sp_project), intent(inout) :: spproj
        character(len=*),  intent(in)    :: pcontrast
        integer,           intent(inout) :: mic_to     ! last micrograph to process
        real,              intent(in)    :: moldiam_max
        integer,           intent(out)   :: box_in_pix
        real,              intent(out)   :: mskdiam    ! estimated mask diameter
        logical,           parameter     :: DEBUG = .true.
        real,              parameter     :: SMPD_SHRINK1  = 4.0,  SIGMA_CRIT = 2., SIGMA_CRIT_MSK = 2.5
        integer,           parameter     :: BOXFAC = 3, NQ_DIAMS = 10
        type(string),      allocatable   :: micnames(:), mic_den_names(:), mic_topo_names(:), mic_bin_names(:)
        integer,           allocatable   :: diam_labels(:), orimap(:)
        real,              allocatable   :: diams_arr(:), diams_arr_ts(:), tmp(:)
        real,              allocatable   :: diam_means(:), abs_z_scores(:)
        type(string)       :: boxfile, fbody_here, ext, fname_thumb_den, str_intg
        type(picksegdiam)  :: picker
        type(image)        :: mic_raw, mic_shrink, mic_den
        type(stats_struct) :: diam_stats
        type(nrtxtfile)    :: diams_file
        integer :: nmics, ldim_raw(3), ldim(3), imic, pop, nptcls, i
        real    :: scale, mad, smpd
        logical :: l_empty
        ! parse project
        smpd = spproj%get_smpd()
        call spproj%get_mics_table(micnames, orimap)
        nmics = size(micnames)
        if( mic_to > nmics )then
            THROW_WARN('mic_to out of range, setting current mic_to='//int2str(mic_to)//', to nmics='//int2str(nmics))
            mic_to = nmics
        endif
        ! read the first micrograph
        scale = smpd / SMPD_SHRINK1
        call read_mic_subtr_backgr_shrink(micnames(1), smpd, scale, pcontrast, mic_raw, mic_shrink, l_empty)
        ! set logical dimensions
        ldim_raw = mic_raw%get_ldim()
        ldim     = mic_shrink%get_ldim()
        allocate(mic_den_names(mic_to), mic_topo_names(mic_to), mic_bin_names(mic_to))
        ! diameter estimation
        do imic = 1, mic_to
            ! skip picking if pre-picked in preproc
            if(spproj%os_mic%isthere(orimap(imic), 'mic_den') &
                &.and. spproj%os_mic%isthere(orimap(imic), 'mic_topo') &
                &.and. spproj%os_mic%isthere(orimap(imic), 'mic_bin') ) then
                    mic_den_names(imic)  = spproj%os_mic%get_str(orimap(imic), 'mic_den')
                    mic_topo_names(imic) = spproj%os_mic%get_str(orimap(imic), 'mic_topo')
                    mic_bin_names(imic)  = spproj%os_mic%get_str(orimap(imic), 'mic_bin')
                    if( spproj%os_mic%isthere(orimap(imic), 'mic_diam') ) then
                        if (file_exists(spproj%os_mic%get_str(orimap(imic), 'mic_diam'))) then
                            call diams_file%new(spproj%os_mic%get_str(orimap(imic), 'mic_diam'), 1)
                            allocate(tmp(diams_file%get_nrecs_per_line()))
                            call diams_file%readNextDataLine(tmp)
                            str_intg = spproj%os_mic%get_str(orimap(imic), 'intg')
                            write(logfhandle, *) ">>> FOUND PICK-PREPROCESSING FOR MICROGRAPH "//str_intg%to_char()&
                            &//". IMPORTED "//int2str(size(tmp))// " DIAMETERS"
                            call spproj%os_mic%delete_entry(orimap(imic), 'mic_diam')
                            call str_intg%kill
                        endif
                    else
                        str_intg = spproj%os_mic%get_str(orimap(imic), 'intg')
                        write(logfhandle, *) ">>> FOUND PICK-PREPROCESSING FOR MICROGRAPH "//str_intg%to_char()&
                            &//". IMPORTED 0 DIAMETERS"
                        call str_intg%kill
                    endif
                    call spproj%os_mic%delete_entry(orimap(imic), 'mic_den')
                    call spproj%os_mic%delete_entry(orimap(imic), 'mic_topo')
                    call spproj%os_mic%delete_entry(orimap(imic), 'mic_bin')
            else
                ! Segmentation & picking
                mic_den_names(imic)  = append2basename(micnames(imic), DEN_SUFFIX)
                mic_topo_names(imic) = append2basename(micnames(imic), TOPO_SUFFIX)
                mic_bin_names(imic)  = append2basename(micnames(imic), BIN_SUFFIX)
                call picker%pick(micnames(imic), smpd, moldiam_max, pcontrast, denfname=mic_den_names(imic),&
                    topofname=mic_topo_names(imic), binfname=mic_bin_names(imic) )
                if( picker%get_nboxes() == 0 ) cycle
                call picker%get_diameters(tmp)
            endif
            if(allocated(tmp)) then
                if( allocated(diams_arr) )then
                    diams_arr = [diams_arr(:), tmp(:)]
                else
                    allocate(diams_arr(size(tmp)), source=tmp)
                endif
                deallocate(tmp)
            endif
        end do
        call picker%kill
        call diams_file%kill
        ! diameter stats & box size estimation
        diams_arr = diams_arr + 2. * SMPD_SHRINK1 ! because of the 2X erosion in binarization
        call calc_stats(diams_arr, diam_stats)
        print *, 'CC diameter (in Angs) statistics'
        print *, 'avg diam: ', diam_stats%avg
        print *, 'med diam: ', diam_stats%med
        print *, 'sde diam: ', diam_stats%sdev
        print *, 'min diam: ', diam_stats%minv
        print *, 'max diam: ', diam_stats%maxv
        allocate( diam_means(NQ_DIAMS), diam_labels(NQ_DIAMS) )
        ! quantization of diameters
        call sortmeans(diams_arr, NQ_DIAMS, diam_means, diam_labels)
        ! Z-score of quantas
        mad = mad_gau(diams_arr, diam_stats%med)
        allocate(abs_z_scores(NQ_DIAMS), source=abs((diam_means - diam_stats%med) / mad))
        do i = 1, NQ_DIAMS
            print *, 'diam quanta '//int2str_pad(i,2)//', avg diam: ', diam_means(i),&
            &', % pop: ', 100 * real(count(diam_labels == i)) / real(size(diams_arr)),&
            &', abs(zscore): ', abs_z_scores(i)
        end do
        ! tresholding
        pop = 0
        do i = 2, NQ_DIAMS - 1
            if( abs_z_scores(i) < SIGMA_CRIT )then
                pop = pop + count(diam_labels == i)
                tmp = pack(diams_arr, mask=diam_labels == i)
                if( allocated(diams_arr_ts) )then
                    diams_arr_ts = [diams_arr_ts(:), tmp(:)]
                    deallocate(tmp)
                else
                    allocate(diams_arr_ts(size(tmp)), source=tmp)
                    deallocate(tmp)
                endif
            endif
        end do
        print *, 'thresholding, % of boxes: ', 100. * real(pop)/real(size(diams_arr))
        ! calculate diams stats after hybrid thresholding
        call calc_stats(diams_arr_ts, diam_stats)
        print *, 'CC diameter (in Angs) statistics after thresholding'
        print *, 'avg diam: ', diam_stats%avg
        print *, 'med diam: ', diam_stats%med
        print *, 'sde diam: ', diam_stats%sdev
        print *, 'min diam: ', diam_stats%minv
        print *, 'max diam: ', diam_stats%maxv
        ! box size estimate in pixels
        box_in_pix = find_magic_box(BOXFAC * nint(diam_stats%med/smpd))
        print *, 'box diam: ', box_in_pix * smpd
        ! mskdiam estimate in A
        mskdiam = min((real(box_in_pix) - COSMSKHALFWIDTH) * smpd, diam_stats%maxv)
        mskdiam = min(diam_stats%avg + SIGMA_CRIT_MSK * diam_stats%sdev, mskdiam)
        print *, 'msk diam: ', mskdiam
        ! re-pick with diameter constraints applied
        do imic = 1, mic_to
            boxfile = basename(fname_new_ext(micnames(imic),'box'))
            call picker%pick(ldim_raw, smpd, mic_bin_names(imic), [diam_stats%minv,diam_stats%maxv])
            nptcls = picker%get_nboxes()
            if( nptcls > 0 )then
                call picker%write_pos_and_diams(boxfile, nptcls, box_in_pix)
                call spproj%set_boxfile(orimap(imic), simple_abspath(boxfile), nptcls=nptcls)
            endif
            if( file_exists(mic_den_names(imic)) )then
                fbody_here      = basename(micnames(imic))
                ext             = fname2ext(fbody_here)
                fbody_here      = get_fbody(fbody_here, ext)
                fname_thumb_den = fbody_here%to_char()//DEN_SUFFIX//JPG_EXT
                call read_mic(mic_den_names(imic), mic_den)
                call mic2thumb(mic_den, fname_thumb_den, l_neg=.true.) ! particles black
                call spproj%os_mic%set(orimap(imic), 'thumb_den', simple_abspath(fname_thumb_den))
            endif
        end do
        ! write output to disk
        call spproj%write_segment_inside('mic')
        if( allocated(orimap) ) deallocate(orimap)
    end subroutine segdiampick_mics

end module simple_mini_stream_utils
