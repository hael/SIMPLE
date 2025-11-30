module simple_volanalyzer
include 'simple_lib.f08'
use simple_dock_vols,        only: dock_vols
use simple_parameters,       only: params_glob
use simple_image,            only: image
use simple_simple_volinterp, only: rotvol
use simple_ori,              only: ori
implicit none

public :: init_volanalyzer, dock_compare_volumes
private
#include "simple_local_flags.inc"

type(string), allocatable :: volnames(:)      !< names of volumes to analyze
type(string), allocatable :: volnames_mirr(:) !< names of mirrored volumes to analyze
real,         allocatable :: corrmat(:,:)     !< similarity matrix
type(dock_vols) :: dvols
integer         :: nvols = 0, ldim(3) = 0

contains

    subroutine init_volanalyzer( voltab )
        class(string), intent(in) :: voltab
        type(image) :: vol_mirr
        integer     :: i, j, ifoo, ldim_read(3), ipair, npairs
        real        :: eul(3), eul_mirr(3), shift(3), shift_mirr(3), cc, cc_mirr
        real        :: smpd_here
        call read_filetable(voltab, volnames)
        nvols = size(volnames)
        if( nvols < 3 ) THROW_HARD('Need at least 3 volumes for analysis')
        ! check that dimensions are consistent across volumes
        call find_ldim_nptcls(volnames(1), ldim, ifoo, smpd=smpd_here)
        do i = 2, nvols
            call find_ldim_nptcls(volnames(i), ldim_read, ifoo, smpd=smpd_here)
            if( any(ldim /= ldim_read) )then
                print *, 'ldim      ', ldim
                print *, 'ldim_read ', ldim_read
                THROW_HARD('Inconsistent volume dimensions!')
            endif
        end do
        ! create mirrored volumes
        allocate(volnames_mirr(nvols))
        do i = 1, nvols
            volnames_mirr(i) = 'vol2analyze'//int2str_pad(i,2)//'_mirr.mrc'
            call vol_mirr%new(ldim, params_glob%smpd)
            call vol_mirr%read(volnames(i))
            call vol_mirr%mirror('x')
            call vol_mirr%write(volnames_mirr(i))
        end do
        allocate(corrmat(nvols,nvols), source=1.)
        ! loop over volume pairs
        write(logfhandle,'(A)') '>>> DOCKING ALL PAIRS OF VOLUMES'
        npairs = (nvols * (nvols - 1)) / 2
        ipair  = 0
        do i = 1, nvols - 1
            do j = i + 1, nvols
                ipair = ipair + 1
                call progress_gfortran(ipair, npairs)
                ! evaluate un-mirrored version
                call dvols%new(volnames(i), volnames(j), params_glob%smpd,&
                &params_glob%hp, params_glob%lp, params_glob%mskdiam)
                call dvols%srch()
                call dvols%get_dock_info(eul, shift, cc)
                call dvols%kill
                ! evaluate mirrored version
                call dvols%new(volnames(i), volnames_mirr(j), params_glob%smpd,&
                &params_glob%hp, params_glob%lp, params_glob%mskdiam)
                call dvols%srch()
                call dvols%get_dock_info(eul_mirr, shift_mirr, cc_mirr)
                call dvols%kill
                ! update corrmat
                corrmat(i,j) = max(cc,cc_mirr)
                corrmat(j,i) = corrmat(i,j)
            end do
        end do
        ! destruct
        call vol_mirr%kill
    end subroutine init_volanalyzer

    subroutine dock_compare_volumes
        type(string) :: volfname
        integer      :: i_medoid, i, funit
        real         :: eul(3), eul_mirr(3), shift(3), shift_mirr(3), cc, cc_mirr
        type(image)  :: vol_medoid, vol_docked
        call medoid_from_smat(corrmat, i_medoid)
        ! write ranked volumes
        write(logfhandle,'(A)') '>>> DOCKING VOLUMES WITH RESPECT TO THE MEDOID'
        do i = 1, nvols
            call progress_gfortran(i, nvols)
            ! evaluate un-mirrored version
            call dvols%new(volnames(i_medoid), volnames(i), params_glob%smpd,&
            &params_glob%hp, params_glob%lp, params_glob%mskdiam)
            call dvols%srch()
            call dvols%get_dock_info(eul, shift, cc)
            call dvols%rotate_target(volnames(i), string('vol_tmp1.mrc'))
            call dvols%kill
            ! evaluate mirrored version
            call dvols%new(volnames(i_medoid), volnames_mirr(i), params_glob%smpd,&
            &params_glob%hp, params_glob%lp, params_glob%mskdiam)
            call dvols%srch()
            call dvols%get_dock_info(eul_mirr, shift_mirr, cc_mirr)
            call dvols%rotate_target(volnames_mirr(i), string('vol_tmp2.mrc'))
            call dvols%kill
            ! file management
            volfname = 'vol_docked'//int2str_pad(i,2)//'.mrc'
            if( cc > cc_mirr )then
                call simple_rename('vol_tmp1.mrc', volfname)
                call del_file('vol_tmp2.mrc')
            else
                call simple_rename('vol_tmp2.mrc', volfname)
                call del_file('vol_tmp1.mrc')
            endif
        end do
        ! calculate volume correlations
        call vol_medoid%new(ldim, params_glob%smpd)
        call vol_medoid%read(volnames(i_medoid))
        call vol_medoid%mask(params_glob%msk, 'soft')
        call vol_docked%new(ldim, params_glob%smpd)
        call fopen(funit, string('volanayze_stats.txt'), 'replace', 'unknown')
        do i = 1, nvols
            volfname = 'vol_docked'//int2str_pad(i,2)//'.mrc'
            call vol_docked%read(volfname)
            call vol_docked%mask(params_glob%msk, 'soft')
            cc = vol_medoid%corr(vol_docked, params_glob%lp, params_glob%hp)
            if( i == i_medoid )then
                write(funit,'(A,1X,I2,1X,A,1X,F7.3,1X,A)') '>>> VOLUME:', i, '>>> CORRELATION:', cc, '***MEDOID***'
            else
                write(funit,'(A,1X,I2,1X,A,1X,F7.3)')      '>>> VOLUME:', i, '>>> CORRELATION:', cc
            endif
        end do
        call fclose(funit)
        ! delete mirrored volumes
        call del_files(volnames_mirr)
        ! destruct class vars
        deallocate(corrmat)
        call volnames%kill
        call volnames_mirr%kill
        call dvols%kill
        ! destruct local vars
        call vol_medoid%kill
        call vol_docked%kill
    end subroutine dock_compare_volumes

end module simple_volanalyzer
