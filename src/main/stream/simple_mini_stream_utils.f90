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

real, parameter :: MOLDIAMS_PICK(6) = [20., 100., 200., 300., 400., 500.] ! in Angstroms, for picksegdiam

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
        integer :: ldim_raw(3), ldim(3), imic, nptcls
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
            call picker%pick(mic_name, smpd, MOLDIAMS_PICK, pcontrast, denfname=mic_den_name,&
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
           ! call picker%write_diameters(mic_diam_name)
            call picker%write_pos_and_diams(mic_diam_name, nptcls)
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
                call picker%pick(micnames(imic), smpd, MOLDIAMS_PICK, pcontrast, denfname=mic_den_names(imic),&
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

    subroutine segdiampick_mics_multi( spproj, pcontrast, mic_to, moldiam_max, box_in_pix, mskdiam )
        ! Build one consolidated picking output across micrographs by:
        ! 1) collecting all candidate diameters,
        ! 2) assigning diameters to fixed bins and rejecting outlier bins via Z-score,
        ! 3) writing one box file per micrograph using a shared final box size.
        class(sp_project),         intent(inout) :: spproj
        integer,                   intent(inout) :: box_in_pix
        real,                      intent(inout) :: mskdiam
        integer,                   intent(inout) :: mic_to     ! last micrograph to process
        character(len=*),          intent(in)    :: pcontrast
        real,                      intent(in)    :: moldiam_max ! kept for API compatibility
        real,                      parameter     :: SMPD_SHRINK1 = 4.0   ! shrink factor used in diameter correction path
        real,                      parameter     :: SIGMA_CRIT   = 3.    ! accepted-bin threshold on absolute Z-score
        real,                      parameter     :: BOXFAC_MAX    = 1.5   ! max expansion factor at very small diameters
        real,                      parameter     :: BOXFAC_MIN    = 1.0   ! floor expansion factor
        real,                      parameter     :: BOXFAC_DECAY_PX = 400. ! BOXFAC reaches BOXFAC_MIN at this diameter in pixels
        logical,                   parameter     :: DEBUG        = .false. ! enables verbose diagnostic printing
        type(string),              allocatable   :: micnames(:), mic_den_names(:), mic_topo_names(:), mic_bin_names(:), mic_diam_names(:)
        integer,                   allocatable   :: orimap(:)
        integer,                   allocatable   :: hac_pops(:)
        real,                      allocatable   :: diams_arr(:), tmp(:), line_data(:)
        real,                      allocatable   :: abs_z_scores(:)
        real,                      allocatable   :: hac_centroids(:), hac_dsum(:), hac_dmin(:), hac_dmax(:)
        type(string)       :: fbody_here, ext, fname_thumb_den, str_intg, mic_diam_name_here, boxfile
        type(picksegdiam)  :: picker
        type(image)        :: mic_den
        type(stats_struct) :: diam_stats
        type(nrtxtfile)    :: diams_file
        integer            :: imic, nptcls, i, nclust_hac, n_accepted, i_acc, box_loc, ibin
        integer            :: nrecs_line, ndatalines, iline
        integer            :: funit, iostat, x_old, y_old, old_box, x_new, y_new
        integer            :: line_data_cap
        integer            :: box_single
        integer            :: diams_len, diams_cap, ncopy
        integer,           allocatable :: cluster_for_bin(:)
        real               :: mad, smpd, pop_scale, diam_here
        real               :: diam_lo, diam_hi, msk_single
        real               :: diam_pix, boxfac_eff
        ! Parse project metadata and clamp the requested micrograph range.
        smpd = spproj%get_smpd()
        call spproj%get_mics_table(micnames, orimap)
        if( mic_to > size(micnames) )then
            THROW_WARN('mic_to out of range, setting current mic_to='//int2str(mic_to)//', to nmics='//int2str(size(micnames)))
            mic_to = size(micnames)
        endif
        ! Initialize dynamic buffers and strict bin boundaries.
        line_data_cap = 0
        diams_len = 0
        diams_cap = 0
        diam_lo = MOLDIAMS_PICK(1)
        diam_hi = MOLDIAMS_PICK(size(MOLDIAMS_PICK))
        allocate(mic_den_names(mic_to), mic_topo_names(mic_to), mic_bin_names(mic_to), mic_diam_names(mic_to))
        ! Pass 1: gather diameters from either preprocessed outputs or fresh picking.
        do imic = 1, mic_to
            ! Reuse preprocessed artifacts when available.
            if(spproj%os_mic%isthere(orimap(imic), 'mic_den') &
                &.and. spproj%os_mic%isthere(orimap(imic), 'mic_topo') &
                &.and. spproj%os_mic%isthere(orimap(imic), 'mic_bin') ) then
                    mic_den_names(imic)  = spproj%os_mic%get_str(orimap(imic), 'mic_den')
                    mic_topo_names(imic) = spproj%os_mic%get_str(orimap(imic), 'mic_topo')
                    mic_bin_names(imic)  = spproj%os_mic%get_str(orimap(imic), 'mic_bin')
                    if( spproj%os_mic%isthere(orimap(imic), 'mic_diam') ) then
                        mic_diam_name_here = spproj%os_mic%get_str(orimap(imic), 'mic_diam')
                        if( file_exists(mic_diam_name_here) ) then
                            mic_diam_names(imic) = mic_diam_name_here
                            call diams_file%new(mic_diam_name_here, 1)
                            nrecs_line = diams_file%get_nrecs_per_line()
                            ndatalines = diams_file%get_ndatalines()
                            if( nrecs_line >= 5 .and. ndatalines > 0 )then
                                allocate(tmp(ndatalines))
                                if( allocated(line_data) .and. line_data_cap < nrecs_line ) deallocate(line_data)
                                if( .not. allocated(line_data) ) then
                                    allocate(line_data(nrecs_line))
                                    line_data_cap = nrecs_line
                                endif
                                do iline = 1, ndatalines
                                    call diams_file%readNextDataLine(line_data)
                                    tmp(iline) = line_data(5)
                                enddo
                            endif
                            call diams_file%kill
                            str_intg = spproj%os_mic%get_str(orimap(imic), 'intg')
                            write(logfhandle, *) ">>> FOUND PICK-PREPROCESSING FOR MICROGRAPH "//str_intg%to_char()&
                            &//". IMPORTED "//int2str(merge(ndatalines, 0, nrecs_line >= 5))// " DIAMETERS"
                            call spproj%os_mic%delete_entry(orimap(imic), 'mic_diam')
                            call str_intg%kill
                        endif
                        call mic_diam_name_here%kill
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
                call picker%pick(micnames(imic), smpd, MOLDIAMS_PICK, pcontrast, denfname=mic_den_names(imic),&
                    topofname=mic_topo_names(imic), binfname=mic_bin_names(imic) )
                if( picker%get_nboxes() == 0 ) cycle
                mic_diam_names(imic) = swap_suffix(micnames(imic), TXT_EXT, '.mrc')
                mic_diam_names(imic) = append2basename(mic_diam_names(imic), DIAMS_SUFFIX)
                call picker%write_pos_and_diams(mic_diam_names(imic), nptcls)
                call picker%get_diameters(tmp)
            endif
            if( allocated(tmp) ) then
                ncopy = size(tmp)
                if( ncopy > 0 )then
                    call ensure_real_capacity(diams_arr, diams_cap, diams_len, diams_len + ncopy)
                    diams_arr(diams_len + 1:diams_len + ncopy) = tmp
                    diams_len = diams_len + ncopy
                endif
                deallocate(tmp)
            endif
        end do
        call picker%kill
        call diams_file%kill
        if( diams_len == 0 ) THROW_HARD('No particle diameters collected across all micrographs')
        if( diams_len < diams_cap )then
            allocate(tmp(diams_len), source=diams_arr(1:diams_len))
            call move_alloc(tmp, diams_arr)
            diams_cap = diams_len
        endif
        ! Build robust global diameter statistics for bin-level filtering.
        diams_arr = diams_arr + 2. * SMPD_SHRINK1 ! because of the 2X erosion in binarization
        ! strict boundary bins: keep only diameters within [MOLDIAMS_PICK(1), MOLDIAMS_PICK(end)]
        tmp = pack(diams_arr, mask=diams_arr >= diam_lo .and. diams_arr <= diam_hi)
        call move_alloc(tmp, diams_arr)
        if( .not. allocated(diams_arr) .or. size(diams_arr) == 0 ) THROW_HARD('No particle diameters within strict boundary bins')
        call calc_stats(diams_arr, diam_stats)
        call print_diam_stats('CC diameter (in Angs) statistics', diam_stats%avg, diam_stats%minv, diam_stats%maxv, med=diam_stats%med, sde=diam_stats%sdev)
        ! Fixed-bin clustering and outlier rejection on cluster means.
        mad = mad_gau(diams_arr, diam_stats%med)
        ! --- fixed-bin path using MOLDIAMS_PICK as bin boundaries ---
        nclust_hac = size(MOLDIAMS_PICK) - 1
        allocate(hac_centroids(nclust_hac), source=0.)
        allocate(hac_dsum(nclust_hac),      source=0.)
        allocate(hac_pops(nclust_hac),      source=0)
        ! cluster populations (unweighted counts); per-cluster min/max
        allocate(hac_dmin(nclust_hac), source=huge(1.))
        allocate(hac_dmax(nclust_hac), source=0.)
        do i = 1, size(diams_arr)
            ibin = bin_index_from_bounds(diams_arr(i), MOLDIAMS_PICK, nclust_hac)
            hac_pops(ibin)   = hac_pops(ibin) + 1
            hac_dsum(ibin)   = hac_dsum(ibin) + diams_arr(i)
            if( diams_arr(i) < hac_dmin(ibin) ) hac_dmin(ibin) = diams_arr(i)
            if( diams_arr(i) > hac_dmax(ibin) ) hac_dmax(ibin) = diams_arr(i)
        enddo
        do i = 1, nclust_hac
            if( hac_pops(i) > 0 )then
                hac_centroids(i) = hac_dsum(i) / real(hac_pops(i))
            else
                ! Keep midpoint for empty bins so reporting remains defined.
                hac_centroids(i) = 0.5 * (MOLDIAMS_PICK(i) + MOLDIAMS_PICK(i + 1))
            endif
        end do
        ! Z-score on cluster means - guard against zero MAD (all diameters identical)
        if( mad > 0. )then
            allocate(abs_z_scores(nclust_hac), source=abs((hac_centroids - diam_stats%med) / mad))
        else
            allocate(abs_z_scores(nclust_hac), source=0.)
        endif
        ! Optional cluster diagnostics.
        if( DEBUG )then
            pop_scale = 100. / real(size(diams_arr))
            do i = 1, nclust_hac
                print *, 'hac cluster '//int2str_pad(i,2)//', avg diam: ', hac_centroids(i),&
                &', min diam: ', hac_dmin(i), ', max diam: ', hac_dmax(i),&
                &', % pop: ', real(hac_pops(i)) * pop_scale,&
                &', abs(mean_z): ', abs_z_scores(i)
            end do
        endif
        n_accepted = count( (abs_z_scores < SIGMA_CRIT) .and. (hac_pops > 0) )
        write(logfhandle, *) ">>> N CLUSTERS ACCEPTED BY Z-SCORE THRESHOLDING: "//int2str(n_accepted)//" OUT OF "//int2str(nclust_hac)
        if( n_accepted == 0 )then
            THROW_WARN('No cluster of diameters passed the Z-score thresholding; returning empty cluster set')
            call spproj%write_segment_inside('mic')
            if( allocated(orimap) ) deallocate(orimap)
            call mic_den%kill
            call mic_diam_name_here%kill
            return
        endif
        box_single = 0
        allocate(cluster_for_bin(nclust_hac), source=0)
        ! Collapse accepted bins into one output box/mask envelope.
        i_acc = 0
        do i = 1, nclust_hac
            if( hac_pops(i) > 0 .and. abs_z_scores(i) < SIGMA_CRIT ) then
                i_acc              = i_acc + 1
                cluster_for_bin(i) = i_acc
                if( DEBUG ) call print_diam_stats('CC diameter (in Angs) statistics for accepted cluster with centroid '//real2str(hac_centroids(i))//':', &
                    hac_centroids(i), hac_dmin(i), hac_dmax(i))
                diam_pix   = hac_dmax(i) / smpd
                boxfac_eff = BOXFAC_MIN + (BOXFAC_MAX - BOXFAC_MIN) * max(0., min(1., (BOXFAC_DECAY_PX - diam_pix) / BOXFAC_DECAY_PX))
                box_loc    = find_magic_box(max(2, nint(boxfac_eff * diam_pix)))
                if( DEBUG )then
                    print *, 'boxfac eff: ', boxfac_eff
                    print *, 'box diam: ', box_loc * smpd
                    print *, 'msk diam: ', (real(box_loc) - COSMSKHALFWIDTH) * smpd
                endif
                box_single = max(box_single, box_loc)
            endif
        end do
        msk_single = (real(box_single) - COSMSKHALFWIDTH) * smpd
        ! Pass 2: write per-micrograph box outputs using accepted bins only.
        do imic = 1, mic_to
            boxfile = basename(fname_new_ext(micnames(imic),'box'))
            nptcls = 0
            ! Add denoised thumb to project metadata.
            if( file_exists(mic_den_names(imic)) )then
                fbody_here      = basename(micnames(imic))
                ext             = fname2ext(fbody_here)
                fbody_here      = get_fbody(fbody_here, ext)
                fname_thumb_den = fbody_here%to_char()//DEN_SUFFIX//JPG_EXT
                call read_mic(mic_den_names(imic), mic_den)
                call mic2thumb(mic_den, fname_thumb_den, l_neg=.true.) ! particles black
                call spproj%os_mic%set(orimap(imic), 'thumb_den', simple_abspath(fname_thumb_den))
            endif
            ! Translate accepted diameters to consolidated box coordinates.
            if( file_exists(mic_diam_names(imic)) )then
                write(logfhandle, *) ">>> MIC_DIAM TO CLUSTER BOXES: "//mic_diam_names(imic)%to_char()
                call fopen(funit, status='REPLACE', action='WRITE', file=boxfile, iostat=iostat)
                if( iostat /= 0 ) THROW_HARD('Failed opening box output file: '//boxfile%to_char())
                call diams_file%new(mic_diam_names(imic), 1)
                nrecs_line = diams_file%get_nrecs_per_line()
                ndatalines = diams_file%get_ndatalines()
                if( nrecs_line >= 5 .and. ndatalines > 0 )then
                    if( allocated(line_data) .and. line_data_cap < nrecs_line ) deallocate(line_data)
                    if( .not. allocated(line_data) )then
                        allocate(line_data(nrecs_line))
                        line_data_cap = nrecs_line
                    endif
                    do iline = 1, ndatalines
                        call diams_file%readNextDataLine(line_data)
                        diam_here = line_data(5)
                        if( diam_here < diam_lo .or. diam_here > diam_hi ) cycle
                        ibin  = bin_index_from_bounds(diam_here, MOLDIAMS_PICK, nclust_hac)
                        i_acc = cluster_for_bin(ibin)
                        if( i_acc <= 0 ) cycle
                        x_old   = nint(line_data(1))
                        y_old   = nint(line_data(2))
                        old_box = nint(line_data(3))
                        x_new   = nint(real(x_old) + 0.5 * real(old_box) - 0.5 * real(box_single))
                        y_new   = nint(real(y_old) + 0.5 * real(old_box) - 0.5 * real(box_single))
                        write(funit,'(4I7,2F8.1)') x_new, y_new, box_single, box_single, diam_here
                        nptcls = nptcls + 1
                    enddo
                endif
                call diams_file%kill
                call fclose(funit)
            endif
            if( nptcls > 0 )then
                call spproj%set_boxfile(orimap(imic), simple_abspath(boxfile), nptcls=nptcls)
            else
                call spproj%set_boxfile(orimap(imic), boxfile, nptcls=0)
            endif
        end do
        ! Remove rejected/invalid micrographs from os_mic:
        ! 1) state=0 (or generally state<1)
        ! 2) missing/invalid boxfile
        do i = spproj%os_mic%get_noris(), 1, -1
            if( spproj%os_mic%get(i, 'state') < 1.0 )then
                call spproj%os_mic%delete(i)
                cycle
            endif
            if( .not. spproj%os_mic%isthere(i, 'boxfile') )then
                call spproj%os_mic%delete(i)
                cycle
            endif
            boxfile = spproj%os_mic%get_str(i, 'boxfile')
            if( boxfile%strlen() == 0 .or. .not. file_exists(boxfile) )then
                call boxfile%kill
                call spproj%os_mic%delete(i)
                cycle
            endif
            call boxfile%kill
        enddo
        ! Persist metadata after box generation for on-disk consistency.
        call spproj%write_segment_inside('mic')
        call spproj%write()
        ! Return a single-output view for callers.
        box_in_pix   = box_single
        mskdiam      = msk_single
        deallocate(hac_centroids, hac_dsum, hac_pops, hac_dmin, hac_dmax, abs_z_scores)
        ! cleanup
        if( allocated(cluster_for_bin) ) deallocate(cluster_for_bin)
        if( allocated(line_data) ) deallocate(line_data)
        if( allocated(tmp) ) deallocate(tmp)
        call mic_den%kill
        call mic_diam_name_here%kill
        if( allocated(orimap) ) deallocate(orimap)
   
    contains

        subroutine ensure_real_capacity(arr, cap, used, needed)
            ! Geometric growth helper for the temporary diameter buffer.
            real,    allocatable, intent(inout) :: arr(:)
            integer,              intent(inout) :: cap
            integer,              intent(in)    :: used, needed
            real, allocatable :: grown(:)
            integer :: new_cap
            if( cap >= needed ) return
            new_cap = max(needed, max(1024, 2 * max(1, cap)))
            if( allocated(arr) )then
                allocate(grown(new_cap))
                if( used > 0 ) grown(1:used) = arr(1:used)
                call move_alloc(grown, arr)
            else
                allocate(arr(new_cap))
            endif
            cap = new_cap
        end subroutine ensure_real_capacity

        subroutine print_diam_stats(label, avgv, minv, maxv, med, sde)
            character(len=*), intent(in) :: label
            real,             intent(in) :: avgv, minv, maxv
            real, optional,   intent(in) :: med, sde
            print *, label
            print *, 'avg diam: ', avgv
            if( present(med) ) print *, 'med diam: ', med
            if( present(sde) ) print *, 'sde diam: ', sde
            print *, 'min diam: ', minv
            print *, 'max diam: ', maxv
        end subroutine print_diam_stats

        integer function bin_index_from_bounds( val, bounds, nbin ) result( ib )
            ! Binary search bin lookup for interval convention: [b1,b2], (b2,b3], ...
            real,    intent(in) :: val
            real,    intent(in) :: bounds(:)
            integer, intent(in) :: nbin
            integer :: l, r, m
            ! Interval convention: [b1,b2], (b2,b3], ..., (b_{N-1},bN].
            l = 1
            r = nbin
            do while( l < r )
                m = (l + r) / 2
                if( val > bounds(m + 1) )then
                    l = m + 1
                else
                    r = m
                endif
            end do
            ib = l
        end function bin_index_from_bounds

    end subroutine segdiampick_mics_multi

end module simple_mini_stream_utils
