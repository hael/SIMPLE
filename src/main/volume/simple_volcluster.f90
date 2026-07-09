!@descr: clustering of pre-docked volumes from Fourier-shell correlations
module simple_volcluster
use simple_core_module_api
use simple_image, only: image
use simple_parameters, only: parameters
implicit none

public :: read_volcluster_volumes, calc_volcluster_dmat, write_volcluster_report
private
#include "simple_local_flags.inc"

contains

    subroutine read_volcluster_volumes( params, volnames, vols )
        class(parameters),        intent(in)    :: params
        type(string),             intent(in)    :: volnames(:)
        type(image), allocatable, intent(inout) :: vols(:)
        integer :: i, ifoo, nvols, ldim(3), ldim_read(3)
        nvols = size(volnames)
        if( nvols < 2 ) THROW_HARD('Need at least 2 volumes for clustering')
        if( params%msk < 0.1 ) THROW_HARD('Invalid mask radius derived from mskdiam')
        do i = 1, nvols
            if( .not. file_exists(volnames(i)) ) THROW_HARD('Missing volume file: '//volnames(i)%to_char())
        end do
        call find_ldim_nptcls(volnames(1), ldim, ifoo)
        if( ldim(3) <= 1 ) THROW_HARD('volcluster expects 3D volumes')
        do i = 2, nvols
            call find_ldim_nptcls(volnames(i), ldim_read, ifoo)
            if( any(ldim /= ldim_read) )then
                write(logfhandle,*) 'ldim first/read: ', ldim, ldim_read
                THROW_HARD('Inconsistent volume dimensions; volcluster')
            endif
        end do
        if( allocated(vols) ) call kill_volumes(vols)
        allocate(vols(nvols))
        write(logfhandle,'(A)') '>>> READING AND MASKING DOCKED VOLUMES'
        do i = 1, nvols
            call vols(i)%new(ldim, params%smpd, wthreads=.false.)
            call vols(i)%read(volnames(i))
            call vols(i)%mask3D_soft(params%msk, backgr=0.)
            call vols(i)%fft
        end do
    end subroutine read_volcluster_volumes

    function calc_volcluster_dmat( params, vols ) result( dmat )
        class(parameters), intent(in)    :: params
        type(image),       intent(inout) :: vols(:)
        real, allocatable :: dmat(:,:), ccmat(:,:)
        integer :: i, j, nvols, npairs, ipair
        nvols = size(vols)
        allocate(ccmat(nvols,nvols), source=1.)
        write(logfhandle,'(A)') '>>> CALCULATING VOLUME CORRELATION DISTANCE MATRIX'
        npairs = (nvols * (nvols - 1)) / 2
        ipair  = 0
        do i = 1, nvols - 1
            do j = i + 1, nvols
                ipair = ipair + 1
                call progress_gfortran(ipair, npairs)
                ccmat(i,j) = vols(i)%corr(vols(j), params%lp, params%hp)
                ccmat(j,i) = ccmat(i,j)
            end do
        end do
        dmat = smat2dmat(ccmat)
        call normalize_minmax(dmat)
        call rmat2file(dmat, string('volcluster_dmat.bin'))
        deallocate(ccmat)
    end function calc_volcluster_dmat

    subroutine write_volcluster_report( volnames, dmat, i_medoids, labels, outfile )
        type(string), intent(in) :: volnames(:)
        real,         intent(in) :: dmat(:,:)
        integer,      intent(in) :: i_medoids(:), labels(:)
        class(string), intent(in) :: outfile
        integer :: i, imed, funit
        logical :: l_medoid
        call fopen(funit, outfile, status='replace', action='write')
        write(funit,'(A)') '# volcluster report'
        write(funit,'(A,I0)') '# nvols=', size(volnames)
        write(funit,'(A,I0)') '# nclusters=', size(i_medoids)
        write(funit,'(A)') '# index cluster medoid distance_to_medoid filename'
        do i = 1, size(volnames)
            imed     = i_medoids(labels(i))
            l_medoid = any(i_medoids == i)
            write(funit,'(I6,1X,I6,1X,L1,1X,ES14.6,1X,A)') &
                i, labels(i), l_medoid, dmat(i,imed), volnames(i)%to_char()
        end do
        call fclose(funit)
    end subroutine write_volcluster_report

    subroutine kill_volumes( vols )
        type(image), allocatable, intent(inout) :: vols(:)
        integer :: i
        if( allocated(vols) )then
            do i = 1, size(vols)
                call vols(i)%kill
            end do
            deallocate(vols)
        endif
    end subroutine kill_volumes

end module simple_volcluster
