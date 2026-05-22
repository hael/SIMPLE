!@descr: utilities for running the mini batch version of the stream
module simple_mini_stream_utils
use simple_core_module_api
use simple_micproc
use simple_sp_project,  only: sp_project
use simple_picksegdiam, only: picksegdiam
use simple_parameters,  only: parameters
use simple_image_bin,   only: image_bin
use simple_image,       only: image
use simple_cmdline,     only: cmdline
use simple_gui_utils,   only: mic2thumb
use simple_nrtxtfile,   only: nrtxtfile
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
        call mic_raw%kill
        call mic_shrink%kill
        call mic_den_name%kill
        call mic_topo_name%kill
        call mic_bin_name%kill
        call mic_diam_name%kill
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
        allocate(diam_means(NQ_DIAMS))
        allocate(diam_labels(size(diams_arr)))
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
        call picker%kill
        ! write output to disk
        call spproj%write_segment_inside('mic')
        if( allocated(orimap) ) deallocate(orimap)
        ! cleanup
        call mic_raw%kill
        call mic_shrink%kill
        call mic_den%kill
    end subroutine segdiampick_mics

    subroutine segdiampick_mics_multi( spproj, projs_clusters, pcontrast, mic_to, moldiam_max, boxes_in_pix, mskdiams, maxmins )
        class(sp_project),         intent(inout) :: spproj
        type(string), allocatable, intent(inout) :: projs_clusters(:)
        integer,      allocatable, intent(inout) :: boxes_in_pix(:)
        real,         allocatable, intent(inout) :: mskdiams(:)
        integer,      allocatable, intent(inout) :: maxmins(:,:)
        integer,                   intent(inout) :: mic_to     ! last micrograph to process
        character(len=*),          intent(in)    :: pcontrast
        real,                      intent(in)    :: moldiam_max
        real,                      parameter     :: SMPD_SHRINK1  = 4.0,  SIGMA_CRIT = 3., SIGMA_CRIT_MSK = 2.5
        real,                      parameter     :: BOXFAC = 1.5, MSKFAC=1.1
        type(string),              allocatable   :: micnames(:), mic_den_names(:), mic_topo_names(:), mic_bin_names(:)
        integer,                   allocatable   :: orimap(:)
        integer,                   allocatable   :: hac_labels(:), hac_pops(:)
        real,                      allocatable   :: diams_arr(:), diams_arr_ts(:), tmp(:)
        real,                      allocatable   :: abs_z_scores(:)
        real,                      allocatable   :: hac_centroids(:), hac_area_pops(:), hac_dmin(:), hac_dmax(:), pop_arr(:), pop_z(:)
        type(string)       :: boxfile, fbody_here, ext, fname_thumb_den, str_intg
        type(string),      allocatable :: boxfiles(:)
        type(picksegdiam)  :: picker
        type(image)        :: mic_raw, mic_shrink, mic_den
        type(stats_struct) :: diam_stats
        type(nrtxtfile)    :: diams_file
        type(sp_project)   :: spproj_tmp
        integer            :: nmics, ldim_raw(3), ldim(3), imic, nptcls, i, nclust_hac, n_accepted, i_acc, box_loc
        integer,           allocatable :: nptcls_by_mic(:)
        real               :: scale, mad, smpd, hac_thresh, n_diams_r, area_sum, pop_scale, area_scale
        real               :: pop_med, pop_mad, msk_loc
        logical            :: l_empty
        ! Parse project and basic dimensions.
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
        ! Diameter estimation pass over micrographs.
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
        if( .not. allocated(diams_arr) ) THROW_HARD('No particle diameters collected across all micrographs')
        n_diams_r = real(size(diams_arr))
        ! diameter stats & box size estimation
        diams_arr = diams_arr + 2. * SMPD_SHRINK1 ! because of the 2X erosion in binarization
        call calc_stats(diams_arr, diam_stats)
        print *, 'CC diameter (in Angs) statistics'
        print *, 'avg diam: ', diam_stats%avg
        print *, 'med diam: ', diam_stats%med
        print *, 'sde diam: ', diam_stats%sdev
        print *, 'min diam: ', diam_stats%minv
        print *, 'max diam: ', diam_stats%maxv
        ! cluster & threshold diameters
        mad = mad_gau(diams_arr, diam_stats%med)
        ! --- hierarchical agglomerative clustering path ---
        hac_thresh = diam_stats%sdev
        allocate(hac_labels(size(diams_arr)))
        call hac_1d_fast(diams_arr, hac_thresh, hac_labels, hac_centroids, hac_pops)
        nclust_hac = size(hac_centroids)
        ! area-weighted cluster populations: sum of d² per cluster; per-cluster min/max
        allocate(hac_area_pops(nclust_hac), source=0.)
        allocate(hac_dmin(nclust_hac), source=huge(1.))
        allocate(hac_dmax(nclust_hac), source=0.)
        do i = 1, size(diams_arr)
            hac_area_pops(hac_labels(i)) = hac_area_pops(hac_labels(i)) + diams_arr(i)**2
            if( diams_arr(i) < hac_dmin(hac_labels(i)) ) hac_dmin(hac_labels(i)) = diams_arr(i)
            if( diams_arr(i) > hac_dmax(hac_labels(i)) ) hac_dmax(hac_labels(i)) = diams_arr(i)
        enddo
        pop_scale = 100. / n_diams_r
        area_sum  = sum(hac_area_pops)
        if( area_sum > 0. )then
            area_scale = 100. / area_sum
        else
            area_scale = 0.
        endif
        ! Z-score on cluster means — guard against zero MAD (all diameters identical)
        if( mad > 0. )then
            allocate(abs_z_scores(nclust_hac), source=abs((hac_centroids - diam_stats%med) / mad))
        else
            allocate(abs_z_scores(nclust_hac), source=0.)
        endif
        ! MAD-based Z-score on area-weighted populations (positive = below-median area = sparse)
        allocate(pop_arr(nclust_hac), source=hac_area_pops)
        pop_med = median_nocopy(pop_arr)
        pop_mad = mad_gau(pop_arr, pop_med)
        if( pop_mad > 0.0 )then
            allocate(pop_z(nclust_hac), source=(pop_med - pop_arr) / pop_mad)
        else
            allocate(pop_z(nclust_hac), source=0.)
        endif
        ! print cluster stats
        do i = 1, nclust_hac
            print *, 'hac cluster '//int2str_pad(i,2)//', avg diam: ', hac_centroids(i),&
            &', min diam: ', hac_dmin(i), ', max diam: ', hac_dmax(i),&
            &', % pop: ', real(hac_pops(i)) * pop_scale,&
            &', % area: ', hac_area_pops(i) * area_scale,&
            &', abs(mean_z): ', abs_z_scores(i), ', pop_z: ', pop_z(i)
        end do
        n_accepted = count( (abs_z_scores < SIGMA_CRIT) .and. (pop_z < SIGMA_CRIT) )
        print *, 'n clusters accepted by Z-score thresholding: ', n_accepted, 'out of ', nclust_hac
        if( n_accepted == 0 ) THROW_HARD('No cluster of diameters passed the Z-score thresholding. Consider increasing SIGMA_CRIT or checking the diameter distribution.')
        if( allocated(projs_clusters) ) deallocate(projs_clusters)
        if( allocated(boxes_in_pix)   ) deallocate(boxes_in_pix)
        if( allocated(mskdiams)       ) deallocate(mskdiams)
        if( allocated(maxmins)        ) deallocate(maxmins)
        allocate(projs_clusters(n_accepted))
        allocate(boxes_in_pix(n_accepted))
        allocate(mskdiams(n_accepted))
        allocate(maxmins(n_accepted, 2))
        ! Add denoised thumbs to project metadata.
        do imic = 1, mic_to
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
        ! Thresholding: reject clusters by outlier mean or anomalously low population.
        i_acc = 0
        do i = 1, nclust_hac
            if( abs_z_scores(i) < SIGMA_CRIT .and. pop_z(i) < SIGMA_CRIT ) then
                i_acc        = i_acc + 1
                spproj_tmp   = spproj
                diams_arr_ts = pack(diams_arr, mask=hac_labels == i)
                call calc_stats(diams_arr_ts, diam_stats)
                print *, 'CC diameter (in Angs) statistics for accepted cluster with centroid ', hac_centroids(i), ':'
                print *, 'avg diam: ', diam_stats%avg
                print *, 'med diam: ', diam_stats%med
                print *, 'sde diam: ', diam_stats%sdev
                print *, 'min diam: ', diam_stats%minv
                print *, 'max diam: ', diam_stats%maxv
                ! box/msk estimate for this accepted cluster; values are stored per cluster in output arrays
                box_loc = find_magic_box(nint(BOXFAC * diam_stats%maxv / smpd))
                msk_loc = (real(box_loc) - COSMSKHALFWIDTH) * smpd
            ! msk_loc = min((real(box_loc) - COSMSKHALFWIDTH) * smpd, diam_stats%maxv * MSKFAC)
            !    msk_loc = min(diam_stats%avg + SIGMA_CRIT_MSK * diam_stats%sdev, msk_loc)
                print *, 'box diam: ', box_loc * smpd
                print *, 'msk diam: ', msk_loc
                ! re-pick with diameter constraints applied
                allocate(boxfiles(mic_to))
                allocate(nptcls_by_mic(mic_to), source=0)
                do imic = 1, mic_to
                    boxfiles(imic) = basename(fname_new_ext(micnames(imic),'box'//int2str(i)))
                enddo
                do imic = 1, mic_to
                    call picker%pick(ldim_raw, smpd, mic_bin_names(imic), [diam_stats%minv,diam_stats%maxv])
                    nptcls = picker%get_nboxes()
                    if( nptcls > 0 )then
                        call picker%write_pos_and_diams(boxfiles(imic), nptcls, box_loc)
                        nptcls_by_mic(imic) = nptcls
                    endif
                end do
                call picker%kill
                do imic = 1, mic_to
                    if( nptcls_by_mic(imic) > 0 )then
                        call spproj_tmp%set_boxfile(orimap(imic), simple_abspath(boxfiles(imic)), nptcls=nptcls_by_mic(imic))
                    endif
                    call boxfiles(imic)%kill
                enddo
                deallocate(boxfiles, nptcls_by_mic)
                ! write output to disk
                projs_clusters(i_acc) = 'cluster_project_'//int2str(i)//METADATA_EXT
                boxes_in_pix(i_acc)   = box_loc
                mskdiams(i_acc)       = msk_loc
                maxmins(i_acc, 1)     = nint(diam_stats%maxv)
                maxmins(i_acc, 2)     = nint(diam_stats%minv)
                call spproj_tmp%write(projs_clusters(i_acc))
                ! tidy
                call spproj_tmp%kill()
                if( allocated(diams_arr_ts) ) deallocate(diams_arr_ts)
            endif
        end do
        deallocate(hac_labels, hac_centroids, hac_pops, hac_area_pops, hac_dmin, hac_dmax, abs_z_scores, pop_arr, pop_z)
        ! cleanup
        call mic_raw%kill
        call mic_shrink%kill
        call mic_den%kill
        if( allocated(orimap) ) deallocate(orimap)
    end subroutine segdiampick_mics_multi

end module simple_mini_stream_utils
