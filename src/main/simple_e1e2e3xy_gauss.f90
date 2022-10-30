module simple_e1e2e3xy_gauss
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_online_var
implicit none

public :: e1e2e3xy_gauss
private
#include "simple_local_flags.inc"

type e1e2e3xy_gauss
    private
    type(online_var), allocatable :: gaumod(:,:)
    integer :: nptcls = 0, recsz = 0
contains
    procedure :: new
    procedure :: update
    procedure :: get_gauss_params
    ! procedure :: get_eul_sdev
    ! procedure :: get_xy_sdev
    procedure :: write
    procedure :: read
    procedure :: kill
end type e1e2e3xy_gauss

contains

    subroutine new( self, nptcls )
        class(e1e2e3xy_gauss), intent(inout) :: self
        integer,               intent(in)    :: nptcls
        real(dp) :: arr(5,4) ! five parameters: e1, e2, e3, x, y & four entries: sumw, mean, var, cnt
        call self%kill
        inquire(iolength=self%recsz) arr(:,:)
        self%nptcls = nptcls
        allocate(self%gaumod(5,self%nptcls))
    end subroutine new

    subroutine update( self, which, vals )
        class(e1e2e3xy_gauss), target, intent(inout) :: self
        character(len=*),              intent(in)    :: which
        real(sp),                      intent(in)    :: vals(:)
        type(online_var), pointer :: param_ptr(:) => null()
        integer :: iptcl
        if( size(vals) .ne. self%nptcls ) THROW_HARD('nonconforming input array of values')
        select case(trim(which))
            case('e1')
                param_ptr => self%gaumod(1,:)
            case('e2')
                param_ptr => self%gaumod(2,:)
            case('e3')
                param_ptr => self%gaumod(3,:)
            case('x')
                param_ptr => self%gaumod(4,:)
            case('y')
                param_ptr => self%gaumod(5,:)
        case DEFAULT
            THROW_HARD('unsupported parameter: '//trim(which))
        end select
        !$omp parallel do default(shared) private(iptcl) schedule(static) proc_bind(close)
        do iptcl = 1,self%nptcls
            call param_ptr(iptcl)%add(vals(iptcl))
        end do
        !$omp end parallel do
    end subroutine update

    subroutine get_gauss_params( self, iptcl, means, sdevs )
        class(e1e2e3xy_gauss), intent(inout) :: self
        integer,               intent(in)    :: iptcl
        real(sp),              intent(inout) :: means(5), sdevs(5)
        means(1) = self%gaumod(1,iptcl)%get_mean()
        means(2) = self%gaumod(2,iptcl)%get_mean()
        means(3) = self%gaumod(3,iptcl)%get_mean()
        means(4) = self%gaumod(4,iptcl)%get_mean()
        means(5) = self%gaumod(5,iptcl)%get_mean()
        sdevs(1) = self%gaumod(1,iptcl)%get_var()
        sdevs(2) = self%gaumod(2,iptcl)%get_var()
        sdevs(3) = self%gaumod(3,iptcl)%get_var()
        sdevs(4) = self%gaumod(4,iptcl)%get_var()
        sdevs(5) = self%gaumod(5,iptcl)%get_var()
        where(sdevs > 1.e-12 )
            sdevs(:) = sqrt(sdevs(:))
        elsewhere
            sdevs(:) = 0.
        end where
    end subroutine get_gauss_params

    subroutine write( self, fname )
        class(e1e2e3xy_gauss), intent(inout) :: self
        character(len=*),      intent(in)    :: fname
        integer  :: funit, io_stat, iptcl
        real(dp) :: arr(5,4)
        if( file_exists(trim(fname)) )then
            call fopen(funit, trim(fname), status='replace', action='write', iostat=io_stat, access='direct', form='unformatted', recl=self%recsz)
            call fileiochk("simple_e1e2e3xy_gauss :: write, fopen failed "//trim(fname), io_stat)
            do iptcl = 1,self%nptcls
                arr(1,:) = self%gaumod(1,iptcl)%serialize()
                arr(2,:) = self%gaumod(2,iptcl)%serialize()
                arr(3,:) = self%gaumod(3,iptcl)%serialize()
                arr(4,:) = self%gaumod(4,iptcl)%serialize()
                arr(5,:) = self%gaumod(5,iptcl)%serialize()
                write(funit, rec=iptcl) arr(:,:)
            end do
            call fclose(funit)
        else
            THROW_HARD('file '//trim(fname)//' does not exist')
        endif
    end subroutine write

    subroutine read( self, fname )
        class(e1e2e3xy_gauss), intent(inout) :: self
        character(len=*),      intent(in)    :: fname
        integer  :: funit, io_stat, iptcl
        real(dp) :: arr(5,4)
        if( file_exists(trim(fname)) )then
            call fopen(funit, fname, status='old', action='read', iostat=io_stat, access='direct', form='unformatted', recl=self%recsz)
            call fileiochk("simple_e1e2e3xy_gauss :: read, fopen failed "//trim(fname), io_stat)
            do iptcl = 1,self%nptcls
                read(funit, rec=iptcl) arr(:,:)
                call self%gaumod(1,iptcl)%unserialize(arr(1,:))
                call self%gaumod(2,iptcl)%unserialize(arr(2,:))
                call self%gaumod(3,iptcl)%unserialize(arr(3,:))
                call self%gaumod(4,iptcl)%unserialize(arr(4,:))
                call self%gaumod(5,iptcl)%unserialize(arr(5,:))
            end do
            call fclose(funit)
        else
            THROW_HARD('file '//trim(fname)//' does not exist')
        endif
    end subroutine read

    subroutine kill( self )
        class(e1e2e3xy_gauss), intent(inout) :: self
        if( allocated(self%gaumod) ) deallocate(self%gaumod)
        self%nptcls = 0
        self%recsz  = 0
    end subroutine kill

end module simple_e1e2e3xy_gauss
