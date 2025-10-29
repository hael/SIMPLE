module simple_mini_stream_utils
include 'simple_lib.f08'
use simple_cmdline,     only: cmdline
use simple_parameters,  only: parameters
use simple_sp_project,  only: sp_project
use simple_binimage,    only: binimage
use simple_image,       only: image
use simple_picksegdiam, only: picksegdiam
use simple_micproc
use simple_gui_utils
implicit none
#include "simple_local_flags.inc"

contains

    subroutine segdiampick_mics( spproj, pcontrast, mic_to, moldiam_max, box_in_pix, mskdiam )
        class(sp_project), intent(inout) :: spproj
        character(len=*),  intent(in)    :: pcontrast
        integer,           intent(in)    :: mic_to     ! last micrograph to process
        real,              intent(in)    :: moldiam_max
        integer,           intent(out)   :: box_in_pix
        real,              intent(out)   :: mskdiam    ! estimated mask diameter
        logical, parameter :: DEBUG = .true.
        real,    parameter :: SMPD_SHRINK1  = 4.0,  SIGMA_CRIT = 2., SIGMA_CRIT_MSK = 2.5
        integer, parameter :: BOXFAC = 3, NQ_DIAMS = 10
        character(len=LONGSTRLEN)              :: boxfile, boxfile_out
        character(len=LONGSTRLEN), allocatable :: micnames(:), mic_den_names(:), mic_topo_names(:), mic_bin_names(:)
        character(len=:),          allocatable :: fbody_here, ext, fname_thumb_den
        integer,                   allocatable :: labels(:), diam_labels(:)
        real,                      allocatable :: diams_arr(:), diams_arr_ts(:), tmp(:)
        real,                      allocatable :: diam_means(:), abs_z_scores(:)
        type(picksegdiam)  :: picker
        type(image)        :: mic_raw, mic_shrink, mic_den
        type(binimage)     :: mic_bin
        type(stats_struct) :: diam_stats
        type(stats_struct) :: stats_nboxes
        integer :: nmics, ldim_raw(3), ldim(3), imic, loc(1), pop, nptcls, i, nboxes
        real    :: scale, mad, smpd
        logical :: l_empty
        ! parse project
        smpd = spproj%get_smpd()
        call spproj%get_mics_table(micnames)
        nmics = size(micnames)
        if( mic_to > nmics ) THROW_HARD('mic_to out of range')
        ! read the first micrograph
        scale = smpd / SMPD_SHRINK1
        call read_mic_subtr_backgr_shrink(micnames(1), smpd, scale, pcontrast, mic_raw, mic_shrink, l_empty)
        ! set logical dimensions
        ldim_raw = mic_raw%get_ldim()
        ldim     = mic_shrink%get_ldim()
        allocate(mic_den_names(nmics), mic_topo_names(nmics), mic_bin_names(nmics))
        ! diameter estimation
        do imic = 1, nmics
            ! Segmentation & picking
            mic_den_names(imic)  = append2basename(micnames(imic), DEN_SUFFIX)
            mic_topo_names(imic) = append2basename(micnames(imic), TOPO_SUFFIX)
            mic_bin_names(imic)  = append2basename(micnames(imic), BIN_SUFFIX)
            call picker%pick(micnames(imic), smpd, moldiam_max, pcontrast, denfname=mic_den_names(imic),&
                topofname=mic_topo_names(imic), binfname=mic_bin_names(imic) )
            call picker%get_diameters(tmp)
            if( allocated(diams_arr) )then
                diams_arr = [diams_arr(:), tmp(:)]
                deallocate(tmp)
            else
                allocate(diams_arr(size(tmp)), source=tmp)
                deallocate(tmp)
            endif
        end do
        call picker%kill
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
        do imic = 1, nmics
            boxfile = basename(fname_new_ext(trim(micnames(imic)),'box'))
            call picker%pick(ldim_raw, smpd, mic_bin_names(imic), [diam_stats%minv,diam_stats%maxv])
            call picker%write_pos_and_diams(boxfile, nptcls, box_in_pix)
            if( nptcls == 0 )then
                boxfile_out = ''
            else
                boxfile_out = simple_abspath(boxfile)
            endif
            call spproj%set_boxfile(imic, boxfile_out, nptcls=nptcls)
            if( file_exists(mic_den_names(imic)) )then
                fbody_here      = basename(trim(micnames(imic)))
                ext             = fname2ext(trim(fbody_here))
                fbody_here      = get_fbody(trim(fbody_here), trim(ext))
                fname_thumb_den = trim(adjustl(fbody_here))//DEN_SUFFIX//trim(JPG_EXT)
                call read_mic(trim(mic_den_names(imic)), mic_den)
                call mic2thumb(mic_den, fname_thumb_den)
                call spproj%os_mic%set(imic, 'thumb_den', simple_abspath(fname_thumb_den))
            else
                call spproj%os_mic%set(imic, 'thumb_den', '')
            endif
        end do
        ! write output to disk
        call spproj%write_segment_inside('mic')
    end subroutine segdiampick_mics

end module simple_mini_stream_utils
