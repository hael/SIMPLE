module simple_stksplit
use simple_defs     ! singleton
use simple_image,   only: image
use simple_imgfile, only: imgfile
implicit none

public :: splitstk_in_parts, splitstk_in_chunks
private

contains

    !> \brief  is splitting a stack in balanced partitions (substacks that are processed one per node)
    !!         this is required to reduce the I/O load in parallel PRIME execution
    subroutine splitstk_in_parts( stkname, ldim, nparts )
        use simple_jiffys, only: get_fileunit, fname2ext, int2str, progress
        character(len=*), intent(in) :: stkname !< name of stack to split
        integer, intent(in)          :: ldim(3) !< logical dimension of stack to split
        integer, intent(in)          :: nparts  !< number of partitions (nodes) for parallel execution
        type(imgfile)                :: imghandle
        type(image)                  :: img
        integer                      :: istop, istart, ptcls_per_part
        integer                      :: leftover, iptcl, ipart, funit
        integer                      :: cnt, nimgs 
        character(len=4)             :: ext
        ext = '.'//fname2ext(stkname)
        call imghandle%open(stkname)
        nimgs = imghandle%getStackSz()
        call imghandle%close
        call img%new(ldim,1.)
        ptcls_per_part = nimgs/nparts
        leftover = nimgs-ptcls_per_part*nparts
        istop    = 0
        istart   = 0
        funit    = get_fileunit()
        open(unit=funit, status='replace', file='stack_partitions.txt')
        do ipart=1,nparts
            call progress(ipart,nparts)
            if( ipart == nparts )then
                istart = istop+1
                istop  = nimgs
            else
                if( leftover == 0 )then
                    istart = istop+1
                    istop  = istart+ptcls_per_part-1
                else
                    istop    = ipart*(ptcls_per_part+1)
                    istart   = istop-(ptcls_per_part+1)+1
                    leftover = leftover-1
                endif
            endif
            write(funit,*) [istart,istop]
            cnt = 0
            do iptcl=istart,istop
                cnt = cnt+1
                call img%read(stkname, iptcl)
                call img%write('stack_part'//int2str(ipart)//ext, cnt)
            end do
        end do
        close(funit)
        call img%kill
    end subroutine
    
    subroutine splitstk_in_chunks( stkname, chunksz )
        character(len=*), intent(in) :: stkname
        integer, intent(in)          :: chunksz

    end subroutine
    
end module simple_stksplit
