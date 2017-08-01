module simple_mapanalyzer
use simple_defs   ! singleton
use simple_image, only: image
implicit none

contains

    subroutine peakfind( density, binmask, excl_radii, winsz, peaklist)
        use simple_sll,    only: sll
        use simple_math,   only: hpsort, euclid
        use simple_jiffys, only: alloc_err
        class(image), intent(inout) :: density, binmask
        real, intent(in)            :: excl_radii
        integer, intent(in)         :: winsz
        class(sll), intent(inout)   :: peaklist
        integer, allocatable        :: inds(:), coords(:,:)
        real, allocatable           :: voxels(:), winvoxels(:) 
        real    :: excl_dist
        integer :: nvox, cnt, ldim(3), alloc_stat, i, j, k
        if( density%is_3d() .and. binmask%is_3d() )then
            ! alles ok
        else
            stop 'only 4 3D vols; peakfind; simple_mapanalyzer'
        endif
        if( .not. (density.eqdims.binmask) ) stop 'ne dims; peakfind; simple_mapanalyzer'
        excl_dist = 2.*excl_radii
        nvox = binmask%nforeground()
        print *, 'nvox:', nvox
        allocate( voxels(nvox), inds(nvox), coords(nvox,3), stat=alloc_stat)
        call alloc_err("In: peakfind; simple_mapanalyzer", alloc_stat)
        ldim = density%get_ldim()
        ! extract voxel & coordinate info
        cnt = 0
        do i=1,ldim(1)
            do j=1,ldim(2)
                do k=1,ldim(3)
                    if( binmask%get([i,j,k]) > 0.5 )then
                        cnt = cnt+1
                        winvoxels     = density%win2arr(i,j,k,winsz)
                        voxels(cnt)   = sum(winvoxels)/real(size(winvoxels))
                        coords(cnt,:) = [i,j,k]
                        inds(cnt)     = cnt
                    endif
                end do
            end do
        end do
        ! sort to find the peaks
        call hpsort(nvox, voxels, inds)
        ! exclude all peaks to are too close to each other
        peaklist = sll()
        binmask = 0.
        do i=nvox,1,-1
            if( inds(i) == 0 )then
                ! this is not a peak
                cycle
            else
                ! store the peak
                call peaklist%add(iarr=coords(inds(i),:), rarr=[voxels(inds(i))])
                call binmask%set(coords(inds(i),:),1.)
                if( i > 1 )then
                    do j=i-1,1,-1
                        if( inds(j) == 0 )then
                            ! this is already excluded
                            cycle 
                        else
                            if( euclid(real(coords(inds(i),:)),real(coords(inds(j),:))) <= excl_dist ) inds(j) = 0
                        endif
                    end do
                end if
            endif
        end do
        deallocate(voxels, inds, coords)
    end subroutine


end module simple_mapanalyzer


