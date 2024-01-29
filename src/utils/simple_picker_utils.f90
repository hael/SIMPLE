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

    subroutine exec_gaupick( micname, boxfile_out, smpd, nptcls, pickrefs, mic_stats )
        character(len=*),          intent(in)    :: micname
        character(len=LONGSTRLEN), intent(out)   :: boxfile_out
        real,                      intent(in)    :: smpd    !< sampling distance in A
        integer,                   intent(out)   :: nptcls
        class(image), optional,    intent(inout) :: pickrefs(:)
        real,         optional,    intent(out)   :: mic_stats(params_glob%nmoldiams,5)
        type(pickgau)             :: gaup, gaup_refine
        real,         allocatable :: moldiams(:)
        integer,      allocatable :: pos(:,:)
        character(len=LONGSTRLEN) :: boxfile
        real    :: maxdiam, stepsz, moldiam_cur
        integer :: box, idiam
        logical :: l_roi
        boxfile = basename(fname_new_ext(trim(micname),'box'))
        l_roi   = trim(params_glob%pick_roi).eq.'yes'
        call read_mic_raw(micname, smpd, subtr_backgr=l_roi)
        if (params_glob%nmoldiams > 1) then
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
        else
            if( present(pickrefs) )then
                call gaup%new_refpicker(       params_glob%pcontrast, SMPD_SHRINK1, params_glob%mskdiam, pickrefs, offset=OFFSET, roi=l_roi)
                call gaup_refine%new_refpicker(params_glob%pcontrast, SMPD_SHRINK2, params_glob%mskdiam, pickrefs, offset=1)
            else
                call gaup%new_gaupicker(       params_glob%pcontrast, SMPD_SHRINK1, params_glob%moldiam, params_glob%moldiam, offset=OFFSET, roi=l_roi)
                call gaup_refine%new_gaupicker(params_glob%pcontrast, SMPD_SHRINK2, params_glob%moldiam, params_glob%moldiam, offset=1)
            endif
            call gaup%gaupick(gaup_refine)
            nptcls = gaup_refine%get_nboxes()
            if( nptcls > 0 )then
                ! write coordinates
                call gaup_refine%get_positions(pos, smpd_new=smpd)
                maxdiam = params_glob%moldiam + params_glob%moldiam * BOX_EXP_FAC
                box     = find_larger_magic_box(round2even(maxdiam / smpd))
                call write_boxfile(nptcls, pos, box, boxfile)
                call make_relativepath(CWD_GLOB, boxfile, boxfile_out) ! returns absolute path
            else
                ! no particles found
                boxfile_out = ''
            endif
            call gaup%kill
            call gaup_refine%kill
        endif
        
    end subroutine exec_gaupick

    subroutine write_boxfile( n, coordinates, box, fname )
        integer,          intent(in) :: n, coordinates(n,2), box
        character(len=*), intent(in) :: fname
        integer :: funit, ibox, iostat
        call fopen(funit, status='REPLACE', action='WRITE', file=trim(adjustl(fname)), iostat=iostat)
        call fileiochk('simple_picker_utils; write_boxfile ', iostat)
        do ibox = 1,n
            write(funit,'(I7,I7,I7,I7,I7)') coordinates(ibox,1), coordinates(ibox,2), box, box, -3
        end do
        call fclose(funit)
    end subroutine write_boxfile

    ! calc_multipick_avgs finds the average smd, ksstat, and a_peak for each moldiam specified in multi_pick, across all micrographs
    ! it creates a csv file called 'stats.csv' which contains this information
    ! and reports the local maxima of each statistic and at which moldiam they occur (prints and writes to binary file)
    subroutine calc_multipick_avgs(spproj, nmoldiams)
        use simple_sp_project
        use simple_math
        use simple_strings, only: int2str
        type(sp_project),  intent(inout) :: spproj
        integer,           intent(in)    :: nmoldiams
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

            do ismd=1, size(max_smds(:,1))
                write(logfhandle,'(1x,2(a,f7.3))') 'Local max smd of ', max_smds(ismd,2), ' occurs at moldiam ', max_smds(ismd,1)
            end do
            loc_max = maxloc(avg_stats(:,2))
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
