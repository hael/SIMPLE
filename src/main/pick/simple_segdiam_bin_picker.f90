!@descr: segmentation-based picking of micrograph batches with one diameter-bin selection and box size for all batches
!==============================================================================
! MODULE: simple_segdiam_bin_picker
!
! PURPOSE:
!   Picks a batch of micrographs with picksegdiam (or reuses the output the
!   preprocessing picked), sorts the particle diameters into the bins bounded
!   by MOLDIAMS_PICK, and writes one box file per micrograph with a single box
!   size. The first batch decides which bins are accepted (bins whose mean
!   diameter lies within SIGMA_CRIT robust z-scores of the batch median) and
!   the box size; later batches reuse both, so every batch of a stream is
!   picked alike. The decision can also be supplied to new().
!
!   segdiampick_mics_multi and segdiampick_mics_multi_fixed_bins
!   (simple_mini_stream_utils) were two 300-line copies of this procedure
!   that differed only in that decision; they are now thin wrappers.
!
! LIFECYCLE:
!   new([accepted_bins, box]) -> { pick(spproj, ...) } -> kill()
!==============================================================================
module simple_segdiam_bin_picker
use simple_core_module_api
use simple_image,       only: image
use simple_micproc,     only: read_mic
use simple_sp_project,  only: sp_project
use simple_picksegdiam, only: picksegdiam
use simple_gui_utils,   only: mic2thumb
use simple_nrtxtfile,   only: nrtxtfile
implicit none

public :: segdiam_bin_picker, MOLDIAMS_PICK
private
#include "simple_local_flags.inc"

real,    parameter :: MOLDIAMS_PICK(6) = [20., 100., 200., 300., 400., 500.] ! diameter bin edges (A), also picksegdiam's diameters
integer, parameter :: NBINS           = size(MOLDIAMS_PICK) - 1
real,    parameter :: SMPD_SHRINK1    = 4.0   ! sampling of the segmentation; its 2x erosion is added back to the diameters
real,    parameter :: SIGMA_CRIT      = 3.    ! accepted-bin threshold on the absolute robust z-score of the bin mean
real,    parameter :: BOXFAC_MAX      = 1.5   ! box expansion factor at very small diameters
real,    parameter :: BOXFAC_MIN      = 1.0   ! floor expansion factor
real,    parameter :: BOXFAC_DECAY_PX = 400.  ! the factor reaches BOXFAC_MIN at this diameter in pixels

type :: segdiam_bin_picker
    private
    logical :: accepted_bins(NBINS) = .false.
    integer :: box                  = 0       ! box size (px) written to the box files
    real    :: mskdiam              = 0.      ! mask diameter (A) for that box; 0 when the box was supplied
    logical :: l_bins_set           = .false.
contains
    procedure :: new
    procedure :: pick
    procedure :: bins_set
    procedure :: get_accepted_bins
    procedure :: get_box
    procedure :: get_mskdiam
    procedure :: kill
end type segdiam_bin_picker

contains

    !> Without arguments the first pick() decides the accepted bins and the box size; with
    !! @p accepted_bins (one per bin of MOLDIAMS_PICK) and @p box every pick() uses those.
    subroutine new( self, accepted_bins, box )
        class(segdiam_bin_picker), intent(inout) :: self
        logical, optional,         intent(in)    :: accepted_bins(:)
        integer, optional,         intent(in)    :: box
        call self%kill()
        if( present(accepted_bins) .neqv. present(box) ) THROW_HARD('accepted_bins and box are supplied together; new')
        if( .not. present(accepted_bins) ) return
        if( size(accepted_bins) /= NBINS ) THROW_HARD('size(accepted_bins) must equal size(MOLDIAMS_PICK)-1; new')
        self%accepted_bins = accepted_bins
        self%box           = box
        self%l_bins_set    = .true.
    end subroutine new

    !> Picks the first @p mic_to micrographs of @p spproj and writes their box files; micrographs
    !! left without particles or box file are removed from the project, which is written.
    subroutine pick( self, spproj, pcontrast, mic_to )
        class(segdiam_bin_picker), intent(inout) :: self
        class(sp_project),         intent(inout) :: spproj
        character(len=*),          intent(in)    :: pcontrast
        integer,                   intent(in)    :: mic_to
        type(string),      allocatable :: micnames(:), mic_den_names(:), mic_topo_names(:), mic_bin_names(:), mic_diam_names(:)
        integer,           allocatable :: orimap(:), bin_pops(:), cluster_for_bin(:)
        real,              allocatable :: diams_arr(:), tmp(:), line_data(:), abs_z_scores(:)
        real,              allocatable :: bin_means(:), bin_dsum(:), bin_dmin(:), bin_dmax(:)
        type(string)       :: fbody_here, ext, fname_thumb_den, str_intg, mic_diam_name_here, boxfile
        type(picksegdiam)  :: picker
        type(image)        :: mic_den
        type(stats_struct) :: diam_stats
        type(nrtxtfile)    :: diams_file
        integer :: nmics, imic, nptcls, i, n_accepted, i_acc, box_loc, ibin, nrecs_line, ndatalines, iline
        integer :: funit, iostat, x_old, y_old, old_box, x_new, y_new, line_data_cap, box_single
        integer :: diams_len, diams_cap, ncopy
        real    :: mad, smpd, diam_here, diam_lo, diam_hi, msk_single, diam_pix, boxfac_eff
        ! project metadata and the micrograph range
        smpd = spproj%get_smpd()
        call spproj%get_mics_table(micnames, orimap)
        nmics = mic_to
        if( nmics > size(micnames) )then
            THROW_WARN('mic_to out of range, setting current mic_to='//int2str(nmics)//', to nmics='//int2str(size(micnames)))
            nmics = size(micnames)
        endif
        line_data_cap = 0
        diams_len     = 0
        diams_cap     = 0
        diam_lo       = MOLDIAMS_PICK(1)
        diam_hi       = MOLDIAMS_PICK(size(MOLDIAMS_PICK))
        allocate(mic_den_names(nmics), mic_topo_names(nmics), mic_bin_names(nmics), mic_diam_names(nmics))
        ! pass 1: diameters, from the preprocessing's picking output when present, else picked here
        do imic = 1, nmics
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
        ! global diameter statistics, within the outer bin edges
        diams_arr = diams_arr + 2. * SMPD_SHRINK1 ! because of the 2X erosion in binarization
        tmp = pack(diams_arr, mask=diams_arr >= diam_lo .and. diams_arr <= diam_hi)
        call move_alloc(tmp, diams_arr)
        if( .not. allocated(diams_arr) .or. size(diams_arr) == 0 ) THROW_HARD('No particle diameters within strict boundary bins')
        call calc_stats(diams_arr, diam_stats)
        call print_diam_stats('CC diameter (in Angs) statistics', diam_stats%avg, diam_stats%minv, diam_stats%maxv,&
            &med=diam_stats%med, sde=diam_stats%sdev)
        ! bin populations, means and ranges
        allocate(bin_means(NBINS), bin_dsum(NBINS), source=0.)
        allocate(bin_pops(NBINS), source=0)
        allocate(bin_dmin(NBINS), source=huge(1.))
        allocate(bin_dmax(NBINS), source=0.)
        do i = 1, size(diams_arr)
            ibin = bin_index_from_bounds(diams_arr(i), MOLDIAMS_PICK, NBINS)
            bin_pops(ibin) = bin_pops(ibin) + 1
            bin_dsum(ibin) = bin_dsum(ibin) + diams_arr(i)
            if( diams_arr(i) < bin_dmin(ibin) ) bin_dmin(ibin) = diams_arr(i)
            if( diams_arr(i) > bin_dmax(ibin) ) bin_dmax(ibin) = diams_arr(i)
        enddo
        do i = 1, NBINS
            if( bin_pops(i) > 0 )then
                bin_means(i) = bin_dsum(i) / real(bin_pops(i))
            else
                bin_means(i) = 0.5 * (MOLDIAMS_PICK(i) + MOLDIAMS_PICK(i + 1)) ! midpoint keeps reporting defined
            endif
        end do
        ! accepted bins: decided by this batch the first time, then reused
        if( self%l_bins_set )then
            n_accepted = count(self%accepted_bins .and. (bin_pops > 0))
            write(logfhandle, *) ">>> N CLUSTERS ACCEPTED (CALLER-SUPPLIED): "//int2str(n_accepted)//" OUT OF "//int2str(NBINS)
        else
            mad = mad_gau(diams_arr, diam_stats%med)
            if( mad > 0. )then
                allocate(abs_z_scores(NBINS), source=abs((bin_means - diam_stats%med) / mad))
            else
                allocate(abs_z_scores(NBINS), source=0.) ! all diameters identical
            endif
            self%accepted_bins = (abs_z_scores < SIGMA_CRIT) .and. (bin_pops > 0)
            self%l_bins_set    = .true.
            n_accepted         = count(self%accepted_bins)
            write(logfhandle, *) ">>> N CLUSTERS ACCEPTED BY Z-SCORE THRESHOLDING: "//int2str(n_accepted)//" OUT OF "//int2str(NBINS)
        endif
        if( n_accepted == 0 )then
            THROW_WARN('No accepted cluster of diameters has any population; returning empty cluster set')
            call spproj%write_segment_inside('mic')
            return
        endif
        ! one box for every accepted bin: decided with the bins, then reused
        allocate(cluster_for_bin(NBINS), source=0)
        box_single = self%box
        i_acc      = 0
        do i = 1, NBINS
            if( bin_pops(i) > 0 .and. self%accepted_bins(i) ) then
                i_acc              = i_acc + 1
                cluster_for_bin(i) = i_acc
                if( self%box > 0 ) cycle
                diam_pix   = bin_dmax(i) / smpd
                boxfac_eff = BOXFAC_MIN + (BOXFAC_MAX - BOXFAC_MIN) * max(0., min(1., (BOXFAC_DECAY_PX - diam_pix) / BOXFAC_DECAY_PX))
                box_loc    = find_magic_box(max(2, nint(boxfac_eff * diam_pix)))
                box_single = max(box_single, box_loc)
            endif
        end do
        msk_single = (real(box_single) - COSMSKHALFWIDTH) * smpd
        ! pass 2: one box file per micrograph with the accepted diameters, re-centred on the common box
        do imic = 1, nmics
            boxfile = basename(fname_new_ext(micnames(imic),'box'))
            nptcls = 0
            ! denoised thumbnail
            if( file_exists(mic_den_names(imic)) )then
                fbody_here      = basename(micnames(imic))
                ext             = fname2ext(fbody_here)
                fbody_here      = get_fbody(fbody_here, ext)
                fname_thumb_den = fbody_here%to_char()//DEN_SUFFIX//JPG_EXT
                call read_mic(mic_den_names(imic), mic_den)
                call mic2thumb(mic_den, fname_thumb_den, l_neg=.true.) ! particles black
                call spproj%os_mic%set(orimap(imic), 'thumb_den', simple_abspath(fname_thumb_den))
            endif
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
                        ibin  = bin_index_from_bounds(diam_here, MOLDIAMS_PICK, NBINS)
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
        ! remove micrographs that are rejected, have no particles, or no readable box file
        do i = spproj%os_mic%get_noris(), 1, -1
            if( spproj%os_mic%get(i, 'state') < 1.0 )then
                call spproj%os_mic%delete(i)
                cycle
            endif
            if( .not. spproj%os_mic%isthere(i, 'nptcls') )then
                call spproj%os_mic%delete(i)
                cycle
            endif
            if( spproj%os_mic%get(i, 'nptcls') <= 0.0 )then
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
        call spproj%write_segment_inside('mic')
        call spproj%write()
        if( self%box == 0 )then
            self%box     = box_single
            self%mskdiam = msk_single
        endif
        call mic_den%kill
    end subroutine pick

    !> .true. once the accepted bins are known, decided by a pick() or supplied to new().
    logical function bins_set( self )
        class(segdiam_bin_picker), intent(in) :: self
        bins_set = self%l_bins_set
    end function bins_set

    function get_accepted_bins( self ) result( accepted_bins )
        class(segdiam_bin_picker), intent(in) :: self
        logical :: accepted_bins(NBINS)
        accepted_bins = self%accepted_bins
    end function get_accepted_bins

    !> The box size (px), 0 until a pick() has decided it.
    integer function get_box( self )
        class(segdiam_bin_picker), intent(in) :: self
        get_box = self%box
    end function get_box

    !> The mask diameter (A) for the decided box, 0 when the box was supplied.
    real function get_mskdiam( self )
        class(segdiam_bin_picker), intent(in) :: self
        get_mskdiam = self%mskdiam
    end function get_mskdiam

    subroutine kill( self )
        class(segdiam_bin_picker), intent(inout) :: self
        self%accepted_bins = .false.
        self%box           = 0
        self%mskdiam       = 0.
        self%l_bins_set    = .false.
    end subroutine kill

    ! geometric growth of the diameter buffer
    subroutine ensure_real_capacity( arr, cap, used, needed )
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

    subroutine print_diam_stats( label, avgv, minv, maxv, med, sde )
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

    ! binary search for the bin of val; intervals [b1,b2], (b2,b3], ..., (b_{N-1},bN]
    integer function bin_index_from_bounds( val, bounds, nbin ) result( ib )
        real,    intent(in) :: val
        real,    intent(in) :: bounds(:)
        integer, intent(in) :: nbin
        integer :: l, r, m
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

end module simple_segdiam_bin_picker
