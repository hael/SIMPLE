module simple_volanalyzer
include 'simple_lib.f08'
use simple_dock_vols,       only: dock_vols
use simple_parameters,      only: params_glob
use simple_image,           only: image
use simple_projector_hlev,  only: rotvol
use simple_ori,             only: ori
implicit none

public :: init_volanalyzer
private
#include "simple_local_flags.inc"

character(len=LONGSTRLEN), allocatable :: volnames(:)      !< names of volumes to analyze
character(len=LONGSTRLEN), allocatable :: volnames_mirr(:) !< names of mirrored volumes to analyze

type dock_inf
    real    :: eul(3)=0., shift(3)=0., cc=0.
    logical :: l_mirr = .false.
end type dock_inf

type(dock_inf), allocatable :: dock_mat(:,:)
type(dock_vols)             :: dvols
integer                     :: nvols = 0

contains

    subroutine init_volanalyzer( voltab )
        character(len=*), intent(in)  :: voltab
        type(image) :: vol_mirr
        integer     :: i, j, ldim(3), ifoo
        real        :: eul(3), eul_mirr(3), shift(3), shift_mirr(3), cc, cc_mirr, smpd_here
        call read_filetable(voltab, volnames)
        nvols = size(volnames)
        if( nvols < 3 ) THROW_HARD('Need at least 3 volumes for analysis')
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
            end do
        end do
        ! destruct
        call vol_mirr%kill
    end subroutine init_volanalyzer

end module simple_volanalyzer
