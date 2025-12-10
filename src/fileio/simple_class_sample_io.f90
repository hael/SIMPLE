module simple_class_sample_io
use simple_defs
use simple_string_utils
use simple_error
use simple_syslib
use simple_fileio
implicit none

public :: print_class_sample, class_samples_same, write_class_samples, read_class_samples, deallocate_class_samples
private
#include "simple_local_flags.inc"

contains
    
    subroutine print_class_sample( cs_entry )
        type(class_sample), intent(in)  :: cs_entry
        print *, 'clsind             ', cs_entry%clsind
        print *, 'pop                ', cs_entry%pop
        print *, 'nsample            ', cs_entry%nsample
        print *, 'size(pinds)        ', size(cs_entry%pinds)
        print *, 'size(cs_entry%ccs) ', size(cs_entry%ccs)
    end subroutine print_class_sample

    function class_samples_same( cs1, cs2 ) result( l_same )
        type(class_sample), intent(inout) :: cs1, cs2
        real, allocatable :: rarr1(:), rarr2(:)
        integer :: sz1, sz2, sz_pinds, sz_ints
        logical :: l_same
        rarr1  = serialize_class_sample(cs1)
        rarr2  = serialize_class_sample(cs2)
        sz1    = size(rarr1)
        sz2    = size(rarr2)
        l_same = .false.
        if( sz1 == sz2 )then
            sz_pinds = (sz1 - 3) / 2
            if( sz_pinds > 0 )then
                sz_ints = 3 + sz_pinds
                if( all(nint(rarr1(:sz_ints)) == nint(rarr2(:sz_ints))) ) l_same = .true.
            endif
        endif
    end function class_samples_same

    function serialize_class_sample( cs_entry ) result( rarr )
        type(class_sample), intent(in)  :: cs_entry
        real, allocatable :: rarr(:)
        integer :: sz_pinds, sz_rarr, cnt, i
        sz_pinds = 0
        if( allocated(cs_entry%pinds) ) sz_pinds = size(cs_entry%pinds)
        sz_rarr  = 3 + 2 * sz_pinds
        allocate(rarr(sz_rarr), source=0.)
        rarr(1)  = real(cs_entry%clsind)
        rarr(2)  = real(cs_entry%pop)
        rarr(3)  = real(cs_entry%nsample)
        if( sz_rarr > 3 )then
            cnt = 3
            do i = 1,sz_pinds
                cnt = cnt + 1
                rarr(cnt) = real(cs_entry%pinds(i))
            end do
            do i = 1,sz_pinds
                cnt = cnt + 1
                rarr(cnt) = cs_entry%ccs(i)
            end do
        endif
    end function serialize_class_sample

    function unserialize_class_sample( rarr ) result( cs_entry )
        real, allocatable, intent(in) :: rarr(:)
        type(class_sample) :: cs_entry
        integer :: sz_pinds, sz_rarr, cnt, i
        if( .not. allocated(rarr) ) THROW_HARD('Input array not allocated')
        sz_rarr          = size(rarr)
        sz_pinds         = (sz_rarr - 3) / 2
        cs_entry%clsind  = nint(rarr(1))
        cs_entry%pop     = nint(rarr(2))
        cs_entry%nsample = nint(rarr(3))
        if( sz_pinds > 0 )then
            allocate(cs_entry%pinds(sz_pinds), source=0 )
            allocate(cs_entry%ccs(sz_pinds),   source=0.)
            cnt = 3
            do i = 1,sz_pinds
                cnt = cnt + 1
                cs_entry%pinds(i) = nint(rarr(cnt)) 
            end do
            do i = 1,sz_pinds
                cnt = cnt + 1
                cs_entry%ccs(i) = rarr(cnt)
            end do
        endif
    end function unserialize_class_sample

    subroutine write_class_samples( csarr, fname )
        type(class_sample), intent(in) :: csarr(:)
        class(string),      intent(in) :: fname
        real, allocatable :: rarr(:), rmat(:,:)
        integer :: i, nx, ny, sz_rarr
        ! turn data into 2D matrix
        ny = size(csarr)
        nx = 0
        do i = 1, ny 
            rarr    = serialize_class_sample(csarr(i))
            sz_rarr = size(rarr)
            if( sz_rarr > nx ) nx = sz_rarr
        end do
        allocate(rmat(ny,nx), source=0.)
        do i = 1, ny 
            rarr    = serialize_class_sample(csarr(i))
            sz_rarr = size(rarr)
            rmat(i,:sz_rarr) = rarr
        end do
        ! write matrix
        call rmat2file(rmat, fname)
    end subroutine write_class_samples

    subroutine read_class_samples( csarr, fname ) 
        type(class_sample), allocatable, intent(inout) :: csarr(:)
        class(string),                   intent(in)    :: fname
        real, allocatable :: rmat(:,:)
        integer :: i, j, nx, ny, cnt
        ! read matrix
        call file2rmat(fname, rmat)
        ny = size(rmat, dim=1)
        nx = size(rmat, dim=2)
        ! deallocate class sample array
        call deallocate_class_samples(csarr)
        ! fill up class sample array
        allocate(csarr(ny))
        do i = 1, ny 
            csarr(i)%clsind  = nint(rmat(i,1))
            csarr(i)%pop     = nint(rmat(i,2))
            csarr(i)%nsample = nint(rmat(i,3))
            allocate(csarr(i)%pinds(csarr(i)%pop), source=0)
            allocate(csarr(i)%ccs(csarr(i)%pop),   source=0.)
            cnt = 3
            do j = 1,csarr(i)%pop
                cnt = cnt + 1
                csarr(i)%pinds(j) = nint(rmat(i,cnt)) 
            end do
            do j = 1,csarr(i)%pop
                cnt = cnt + 1
                csarr(i)%ccs(j) = rmat(i,cnt) 
            end do
        end do
    end subroutine read_class_samples

    subroutine deallocate_class_samples( csarr )
        type(class_sample), allocatable, intent(inout) :: csarr(:)
        integer :: sz, i
        if( allocated(csarr) )then
            sz = size(csarr)
            do i = 1, sz
                if( allocated(csarr(i)%pinds) ) deallocate(csarr(i)%pinds)
                if( allocated(csarr(i)%ccs)   ) deallocate(csarr(i)%ccs)
            end do
            deallocate(csarr)
        endif
    end subroutine deallocate_class_samples

end module simple_class_sample_io
