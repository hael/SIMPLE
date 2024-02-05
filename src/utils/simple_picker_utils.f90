module simple_picker_utils 
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_parameters, only: params_glob
use simple_image,      only: image
use simple_pickgau
implicit none

public :: exec_gaupick, calc_multipick_avgs
private
#include "simple_local_flags.inc"

real,    parameter :: SMPD_SHRINK1 = 4.0, SMPD_SHRINK2 = 2.0, BOX_EXP_FAC = 0.5
integer, parameter :: OFFSET       = 3
logical, parameter :: L_WRITE      = .false.
logical, parameter :: L_DEBUG      = .false.

contains

    subroutine exec_gaupick( micname, boxfile_out, smpd, nptcls, pickrefs, mic_stats, dir_out )
        use simple_strings, only: str2real, parsestr
        character(len=*),           intent(in)    :: micname
        character(len=LONGSTRLEN),  intent(out)   :: boxfile_out
        real,                       intent(in)    :: smpd    !< sampling distance in A
        integer,                    intent(out)   :: nptcls
        character(len=*), optional, intent(in)    :: dir_out
        class(image),     optional, intent(inout) :: pickrefs(:)
        real,             optional, intent(out)   :: mic_stats(params_glob%nmoldiams,5)
        type(pickgau)             :: gaup, gaup_refine
        real,         allocatable :: moldiams(:)
        integer,      allocatable :: pos(:,:)
        character(len=LONGSTRLEN) :: boxfile
        character(len=STDLEN)   :: multi_moldiams
        character(len=4), allocatable :: moldiams_str(:)
        real    :: maxdiam, stepsz, moldiam_cur, pick_stats(5), moldiam_entry
        integer :: box, idiam, istr, num_entries, num_moldiams
        logical :: l_roi, l_backgr_subtr
        boxfile = basename(fname_new_ext(trim(micname),'box'))
        if( present(dir_out) ) boxfile = trim(dir_out)//'/'//trim(boxfile)
        l_roi          = trim(params_glob%pick_roi).eq.'yes'
        l_backgr_subtr = l_roi .or. (trim(params_glob%backgr_subtr).eq.'yes')
        call read_mic_raw(micname, smpd, subtr_backgr=l_backgr_subtr)
        if (params_glob%nmoldiams > 1) then
            ! multiple moldiam pick to assess best molecular diameters, does not generate .box files
            allocate(moldiams(params_glob%nmoldiams))
            ! create moldiams array
            stepsz = (params_glob%moldiam_max - params_glob%moldiam) / (params_glob%nmoldiams - 1)
            moldiam_cur = params_glob%moldiam
            do idiam = 1, params_glob%nmoldiams
                moldiams(idiam) = moldiam_cur
                moldiam_cur = moldiam_cur + stepsz
            enddo
            call gaupick_multi(params_glob%pcontrast, SMPD_SHRINK1, moldiams, offset=OFFSET, mic_stats=mic_stats)
            deallocate(moldiams)
        else if(.not. (params_glob%multi_moldiams  .eq. '')) then
            ! multiple moldiam pick that uses multiple gaussians, generates .box file outputs
            istr=1
            num_entries=0
            ! find number of molecular diameters
            do 
                if(params_glob%multi_moldiams(istr:istr) .eq. ' ') then
                    num_entries = num_entries + 1
                    exit
                end if
                if(istr > len(params_glob%multi_moldiams)) exit
                if(params_glob%multi_moldiams(istr:istr) .eq. ',') then
                    num_entries = num_entries + 1
                end if
                istr = istr + 1
            end do
            ! parse moldiams from comma-separated string of numbers into real-valued array
            allocate(moldiams_str(num_entries))
            call parsestr(params_glob%multi_moldiams,',',moldiams_str,num_moldiams)
            allocate(moldiams(num_moldiams))
            do idiam = 1, num_moldiams
                moldiams(idiam) = str2real(moldiams_str(idiam))
            end do
            deallocate(moldiams_str)
            ! execute multiple gaussian pick
            call gaup%new_gaupicker_multi(       params_glob%pcontrast, SMPD_SHRINK1, moldiams, offset=OFFSET, roi=l_roi)
            call gaup_refine%new_gaupicker_multi(params_glob%pcontrast, SMPD_SHRINK2, moldiams, offset=1)
            call gaup%gaupick(gaup_refine)
            ! gather summary statistics so averages can be calculated
            call gaup%get_stats(pick_stats)
            mic_stats(1,1) = moldiams(1)
            mic_stats(1,2) = pick_stats(1)
            mic_stats(1,3) = pick_stats(2)
            mic_stats(1,4) = pick_stats(3)
            mic_stats(1,5) = pick_stats(4)
            maxdiam = maxval(moldiams) + maxval(moldiams) * BOX_EXP_FAC
            box     = find_larger_magic_box(round2even(maxdiam / smpd))
            call gaup_refine%report_boxfile(box, smpd, boxfile, nptcls)
            if( nptcls == 0 )then
                boxfile_out = ''
            else
                call make_relativepath(CWD_GLOB, boxfile, boxfile_out)
            endif
            call gaup%kill
            call gaup_refine%kill
            deallocate(moldiams)
        else
            ! single moldiam pick
            if( present(pickrefs) )then
                call gaup%new_refpicker(       params_glob%pcontrast, SMPD_SHRINK1, params_glob%mskdiam, pickrefs, offset=OFFSET, roi=l_roi)
                call gaup_refine%new_refpicker(params_glob%pcontrast, SMPD_SHRINK2, params_glob%mskdiam, pickrefs, offset=1)
            else
                call gaup%new_gaupicker(       params_glob%pcontrast, SMPD_SHRINK1, params_glob%moldiam, params_glob%moldiam, offset=OFFSET, roi=l_roi)
                call gaup_refine%new_gaupicker(params_glob%pcontrast, SMPD_SHRINK2, params_glob%moldiam, params_glob%moldiam, offset=1)
            endif
            call gaup%gaupick(gaup_refine)
            ! write
            maxdiam = params_glob%moldiam + params_glob%moldiam * BOX_EXP_FAC
            box     = find_larger_magic_box(round2even(maxdiam / smpd))
            call gaup_refine%report_boxfile(box, smpd, boxfile, nptcls)
            if( nptcls == 0 )then
                boxfile_out = ''
            else
                call make_relativepath(CWD_GLOB, boxfile, boxfile_out)
            endif
            call gaup%kill
            call gaup_refine%kill
        endif
    end subroutine exec_gaupick

    ! calc_multipick_avgs finds the average smd, ksstat, and a_peak for each moldiam specified in multi_pick, across all micrographs
    ! it creates a csv file called 'stats.csv' which contains this information
    ! and reports the local maxima of each statistic and at which moldiam they occur (prints and writes to binary file)
    subroutine calc_multipick_avgs(spproj, nmoldiams)
        use simple_sp_project
        use simple_math
        use simple_strings, only: int2str
        type(sp_project),            intent(inout) :: spproj
        integer,                     intent(in)    :: nmoldiams
        real, allocatable :: mic_stats(:,:) 
        real, allocatable :: avg_stats(:,:)
        real, allocatable :: all_stats(:,:,:)
        real, allocatable :: max_smds(:,:)
        real, allocatable :: max_ksstats(:,:)
        real, allocatable :: max_a_peaks(:,:)
        integer           :: nmicrographs, imic, imic_3, idiam, idiam_2, idiam_3, iapeak, ismd, iksstat, loc_max(1)
   
        nmicrographs = spproj%os_mic%get_noris()
        allocate(mic_stats(nmoldiams,5))
        allocate(all_stats(nmicrographs,nmoldiams,5))
        allocate(avg_stats(nmoldiams,5))

        do imic = 1, nmicrographs 
            do idiam = 1, nmoldiams
                mic_stats(idiam,1) = spproj%os_mic%get(imic,'m_'//int2str(idiam))
                mic_stats(idiam,2) = spproj%os_mic%get(imic,'s_'//int2str(idiam))
                mic_stats(idiam,3) = spproj%os_mic%get(imic,'k_'//int2str(idiam))
                mic_stats(idiam,4) = spproj%os_mic%get(imic,'a_'//int2str(idiam))
                mic_stats(idiam,5) = spproj%os_mic%get(imic,'sp_'//int2str(idiam))
            end do 
            all_stats(imic,:,:) = mic_stats
        end do
        
        ! csv to save picker data
        open(unit=1, file='stats.csv', status='unknown')
        write(1, '(9(a,1x))') 'moldiam', ',', 'smd', ',', 'ksstat', ',', 'a_peak', ',', 's_peak'
        do idiam_2 = 1, nmoldiams
            avg_stats(idiam_2,1) = all_stats(1,idiam_2,1)
            avg_stats(idiam_2,2) = sum(all_stats(:,idiam_2,2))/size(all_stats(:,idiam_2,2))
            avg_stats(idiam_2,3) = sum(all_stats(:,idiam_2,3))/size(all_stats(:,idiam_2,3))
            avg_stats(idiam_2,4) = sum(all_stats(:,idiam_2,4))/size(all_stats(:,idiam_2,4))
            avg_stats(idiam_2,5) = sum(all_stats(:,idiam_2,5))/size(all_stats(:,idiam_2,5))
            write(1, '(5(f7.3,a))') avg_stats(idiam_2,1), ',', avg_stats(idiam_2,2), ',', avg_stats(idiam_2,3), ',', avg_stats(idiam_2,4), ',', avg_stats(idiam_2,5)
        end do
        close(1)

        ! find peak values of statistics, if executing multipick
        if (nmoldiams > 1) then
            max_smds = find_local_maxima(avg_stats(:,1),avg_stats(:,2),nmoldiams)
            max_ksstats = find_local_maxima(avg_stats(:,1),avg_stats(:,3),nmoldiams)
            max_a_peaks = find_local_maxima(avg_stats(:,1),avg_stats(:,4),nmoldiams)

            loc_max = maxloc(avg_stats(:,2))
            do ismd=1, size(max_smds(:,1))
                write(logfhandle,'(1x,2(a,f7.3))') 'Local max smd of ', max_smds(ismd,2), ' occurs at moldiam ', max_smds(ismd,1)
            end do
            write(logfhandle,'(1x,2(a,f7.3))') 'Absolute max smd of ', avg_stats(loc_max(1),2), ' occurs at moldiam ', avg_stats(loc_max(1),1)
            
            do iksstat=1, size(max_ksstats(:,1))
                write(logfhandle,'(1x,2(a,f7.3))') 'Local max ksstat of ', max_ksstats(iksstat,2), ' occurs at moldiam ', max_ksstats(iksstat,1)
            end do
            loc_max = maxloc(avg_stats(:,3))
            write(logfhandle,'(1x,2(a,f7.3))') 'Absolute max ksstat of ', avg_stats(loc_max(1),3), ' occurs at moldiam ', avg_stats(loc_max(1),1)

            do iapeak=1, size(max_a_peaks(:,1))
                write(logfhandle,'(1x,2(a,f7.3))') 'Local max a_peak of ', max_a_peaks(iapeak,2), ' occurs at moldiam ', max_a_peaks(iapeak,1)
            end do 
            loc_max = maxloc(avg_stats(:,4))
            write(logfhandle,'(1x,2(a,f7.3))') 'Absolute max a_peak of ', avg_stats(loc_max(1),4), ' occurs at moldiam ', avg_stats(loc_max(1),1)


            deallocate(max_smds)
            deallocate(max_ksstats)
            deallocate(max_a_peaks)
            
        end if
    
        ! deallocate arrays created to process picker statistics
        deallocate(mic_stats)
        deallocate(all_stats)
        deallocate(avg_stats)

        ! cleaning project file
        do imic_3 = 1, nmicrographs
            do idiam_3 = 1, nmoldiams
                call spproj%os_mic%delete_entry(imic_3,'m_'//int2str(idiam_3))
                call spproj%os_mic%delete_entry(imic_3,'s_'//int2str(idiam_3))
                call spproj%os_mic%delete_entry(imic_3,'k_'//int2str(idiam_3))
                call spproj%os_mic%delete_entry(imic_3,'a_'//int2str(idiam_3))
                call spproj%os_mic%delete_entry(imic_3,'p_'//int2str(idiam_3))
            end do
        end do

    end subroutine calc_multipick_avgs

end module simple_picker_utils
