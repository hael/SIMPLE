module simple_imgcorrections
implicit none

contains

    subroutine replace_deadhot_pixels( movieframes, nsig, deadhot ) 
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_image,   only: image
        use simple_stat,    only: moment
        use simple_math,    only: median
        class(image),   intent(inout) :: movieframes(:)
        integer,        intent(out)   :: deadhot(2)
        real, optional, intent(in)    :: nsig
        integer                       :: nmovies,nframes,frame,i,j,ncured,ldim(3),hwinsz,winsz
        real,    allocatable          :: rmat(:,:,:), rmat_pad(:,:),win (:,:)
        logical, allocatable          :: outliers(:,:)
        real                          :: smpd, nsig_here
        type(image)                   :: img_sum
        logical                       :: isft
        logical, parameter            :: DEBUG = .false.
        if( movieframes(1)%is_ft() )then
            stop 'cannot replace dead/hot pixels in FT; simple_jiffys :: replace_deadhot_pixels'
        endif
        ! init
        nframes   = size(movieframes)
        smpd      = movieframes(1)%get_smpd()
        ldim      = movieframes(1)%get_ldim()
        ldim(3)   = 1
        nsig_here = 6.
        if( present(nsig) ) nsig_here = nsig
        ! create sum
        call img_sum%new([ldim(1),ldim(2),1], smpd)
        img_sum = 0.
        ! allocate padded matrix
        hwinsz  = 6
        winsz   = 2*hwinsz+1
        allocate(rmat_pad(1-hwinsz:ldim(1)+hwinsz, 1-hwinsz:ldim(2)+hwinsz), win(winsz,winsz))
        do frame=1,nframes
            call img_sum%add(movieframes(frame))
        end do
        call img_sum%cure_outliers(ncured, nsig_here, deadhot, outliers)  
        if( any(outliers) )then
            !$omp parallel do schedule(static,1) default(shared) private(frame,rmat,rmat_pad,i,j,win)
            do frame=1,nframes
                rmat = movieframes(frame)%get_rmat()
                rmat_pad = median( reshape(rmat(:,:,1),(/(ldim(1)*ldim(2))/)) )
                rmat_pad(1:ldim(1),1:ldim(2)) = rmat(1:ldim(1),1:ldim(2),1)
                do i=1,ldim(1)
                    do j=1,ldim(2)
                        if( outliers(i,j))then
                            win = rmat_pad( i-hwinsz:i+hwinsz, j-hwinsz:j+hwinsz )
                            rmat(i,j,1) = median( reshape(win,(/winsz**2/)) )
                        endif
                    enddo
                enddo
                call movieframes(frame)%set_rmat(rmat)
            enddo
            !$omp end parallel do
        endif
    end subroutine replace_deadhot_pixels

end module simple_imgcorrections
