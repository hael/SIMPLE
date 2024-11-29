module simple_volanalyzer
include 'simple_lib.f08'
use simple_dock_vols,       only: dock_vols
use simple_parameters,      only: params_glob
use simple_image,           only: image
use simple_projector_hlev,  only: rotvol
use simple_ori,             only: ori
implicit none

public :: init_volanalyzer, rank_dock_compare_volumes
private
#include "simple_local_flags.inc"

character(len=LONGSTRLEN), allocatable :: volnames(:)        !< names of volumes to analyze
character(len=LONGSTRLEN), allocatable :: volnames_mirr(:)   !< names of mirrored volumes to analyze
character(len=LONGSTRLEN), allocatable :: volnames_ranked(:) !< names in ranked order 
character(len=LONGSTRLEN), allocatable :: volnames_docked(:) !< names of docked volumes, indexed according to rank

type dock_inf
    real    :: eul(3)=0., shift(3)=0., cc=0.
    logical :: l_mirr = .false.
end type dock_inf

type(dock_inf), allocatable :: dock_mat(:,:)
real,           allocatable :: corrmat(:,:)
type(dock_vols)             :: dvols
integer                     :: nvols = 0, ldim(3) = 0

contains

    subroutine init_volanalyzer( voltab )
        character(len=*), intent(in)  :: voltab
        type(image) :: vol_mirr
        integer     :: i, j, ifoo, ldim_read(3)
        real        :: eul(3), eul_mirr(3), shift(3), shift_mirr(3), cc, cc_mirr, smpd_here
        call read_filetable(voltab, volnames)
        nvols = size(volnames)
        if( nvols < 3 ) THROW_HARD('Need at least 3 volumes for analysis')
        ! check that dimensions are consistent across volumes
        call find_ldim_nptcls(trim(volnames(1)), ldim, ifoo, smpd=smpd_here)
        do i = 2, nvols
            call find_ldim_nptcls(trim(volnames(i)), ldim_read, ifoo, smpd=smpd_here)
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
            call find_ldim_nptcls(trim(volnames(i)), ldim, ifoo, smpd=smpd_here)
            call vol_mirr%new(ldim, params_glob%smpd)
            call vol_mirr%read(trim(volnames(i)))
            call vol_mirr%mirror('x')
            call vol_mirr%write(trim(volnames_mirr(i)))
        end do
        allocate(dock_mat(nvols,nvols))
        ! loop over volume pairs
        do i = 1, nvols - 1
            do j = i + 1, nvols
                ! evaluate un-mirrored version
                call dvols%new(trim(volnames(i)), trim(volnames(j)), params_glob%smpd,&
                &params_glob%hp, params_glob%lp, params_glob%mskdiam)
                call dvols%srch()
                call dvols%get_dock_info(eul, shift, cc)
                call dvols%kill
                ! evaluate mirrored version
                call dvols%new(trim(volnames(i)), trim(volnames_mirr(j)), params_glob%smpd,&
                &params_glob%hp, params_glob%lp, params_glob%mskdiam)
                call dvols%srch()
                call dvols%get_dock_info(eul_mirr, shift_mirr, cc_mirr)
                call dvols%kill
                ! update matrix
                if( cc_mirr > cc )then
                    dock_mat(i,j)%eul    = eul_mirr
                    dock_mat(i,j)%shift  = shift_mirr
                    dock_mat(i,j)%cc     = cc_mirr
                    dock_mat(i,j)%l_mirr = .true.
                else
                    dock_mat(i,j)%eul    = eul
                    dock_mat(i,j)%shift  = shift
                    dock_mat(i,j)%cc     = cc
                    dock_mat(i,j)%l_mirr = .false.
                endif
                dock_mat(j,i) = dock_mat(i,j)
            end do
        end do
        ! set diagonal elements
        do i = 1, nvols
            dock_mat(i,i)%eul    = 0.
            dock_mat(i,i)%shift  = 0.
            dock_mat(i,i)%cc     = 1.
            dock_mat(i,i)%l_mirr = .false.
        end do
        ! construct the correlation matrix
        allocate(corrmat(nvols,nvols), source=dock_mat(:,:)%cc)
        ! destruct
        call vol_mirr%kill
    end subroutine init_volanalyzer

    subroutine rank_dock_compare_volumes
        use simple_opt_filter, only: estimate_lplim
        integer,          allocatable :: rank(:)
        character(len=:), allocatable :: volfname
        type(image) :: vol, vol_docked, mskvol
        integer     :: i_medoid, i, j, npix, fhandle
        real        :: lpopt
        ! rank according to similarity to medoid
        call medoid_ranking_from_smat(corrmat, i_medoid, rank)
        ! write ranked filetable
        allocate(volnames_ranked(nvols), volnames_docked(nvols))
        do i = 1, nvols
            volnames_ranked(i) = volnames(rank(i))
        end do
        call write_filetable('vols_ranked.txt', volnames_ranked)
        ! write ranked volumes
        do i = 1, nvols
            if( dock_mat(i_medoid,rank(i))%l_mirr )then
                ! the mirrored version is the one
                volfname = trim(volnames_mirr(rank(i)))
                call dvols%new(trim(volnames(i_medoid)), volfname, params_glob%smpd,&
                &params_glob%hp, params_glob%lp, params_glob%mskdiam)
            else
                ! the un-mirrored version is the one
                volfname = trim(volnames(rank(i)))
                call dvols%new(trim(volnames(i_medoid)), volfname, params_glob%smpd,&
                &params_glob%hp, params_glob%lp, params_glob%mskdiam)
            endif
            call dvols%set_dock_info(dock_mat(i_medoid,rank(i))%eul, dock_mat(i_medoid,rank(i))%shift)
            volnames_docked(i) = 'vol_docked'//int2str_pad(i,2)//'.mrc'
            call dvols%rotate_target(volfname, trim(volnames_docked(i)))
            
            call dvols%kill
        end do
        call vol%new(       ldim, params_glob%smpd)
        call vol_docked%new(ldim, params_glob%smpd)
        call mskvol%disc(   ldim, params_glob%smpd, real(min(params_glob%box/2, int(params_glob%msk + COSMSKHALFWIDTH))))
        call vol%read(trim(volnames_docked(1)))
        call fopen(fhandle, 'volanalyze_stats.txt', 'replace', 'unknown') 
        write(fhandle,'(A,1X,I3,1X,A,1X,F6.2,1X,A,1X,F7.3)') '>>> RANK:', 1, '>>> RESOLUTION:', params_glob%lpstop, '>>> CORRELATION:', dock_mat(i_medoid,rank(1))%cc
        do i = 2, nvols
            call vol_docked%read(trim(volnames_docked(i)))
            call estimate_lplim(vol, vol_docked, mskvol, [params_glob%lpstart,params_glob%lpstop], lpopt)
            write(fhandle,'(A,1X,I3,1X,A,1X,F6.2,1X,A,1X,F7.3)') '>>> RANK:', i, '>>> RESOLUTION:', lpopt, '>>> CORRELATION:', dock_mat(i_medoid,rank(i))%cc
        end do
        call fclose(fhandle)
        ! destruct
        call vol%kill
        call vol_docked%kill
        call mskvol%kill
    end subroutine rank_dock_compare_volumes

end module simple_volanalyzer
