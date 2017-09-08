! spectral signal-to-noise ratio estimation routines
module simple_estimate_ssnr
use simple_defs   ! use all in there
use simple_image, only: image
implicit none

contains

    !>  \brief generate projection FRCs from even/odd pairs
    subroutine gen_projection_frcs( p, ename, oname, state, projfrcs )
        use simple_params,          only: params
        use simple_oris,            only: oris
        use simple_projector_hlev,  only: projvol
        use simple_projection_frcs, only: projection_frcs
        class(params),          intent(inout) :: p
        character(len=*),       intent(in)    :: ename, oname
        integer,                intent(in)    :: state
        class(projection_frcs), intent(inout) :: projfrcs
        type(oris)               :: e_space
        type(image)              :: even, odd
        type(image), allocatable :: even_imgs(:), odd_imgs(:)
        real,        allocatable :: frc(:), res(:)
        integer :: iproj
        ! read even/odd pair
        call even%new([p%box,p%box,p%box], p%smpd)
        call odd%new([p%box,p%box,p%box], p%smpd)
        call even%read(ename)
        call odd%read(oname)
        ! create e_space
        call e_space%new(p%nspace)
        call e_space%spiral(p%nsym, p%eullims)
        ! generate even/odd projections
        even_imgs = projvol(even, e_space, p)
        odd_imgs  = projvol(even, e_space, p)
        ! calculate FRCs and fill-in projfrcs object
        do iproj=1,p%nspace
            call even_imgs(iproj)%fwd_ft
            call odd_imgs(iproj)%fwd_ft
            call even_imgs(iproj)%fsc(odd_imgs(iproj), res, frc)
            call projfrcs%set_frc(iproj, frc, state)
            call even_imgs(iproj)%kill
            call odd_imgs(iproj)%kill
        end do
        deallocate(even_imgs, odd_imgs)
        call even%kill
        call odd%kill
    end subroutine gen_projection_frcs

    !> \brief  converts the FSC to SSNR (the 2.* is because of the division of the data)
    function fsc2ssnr( corrs ) result( ssnr )
        real, intent(in)  :: corrs(:) !<  instrument FSC 
        real, allocatable :: ssnr(:) !<  instrument SSNR
        integer :: nyq, k
        real    :: fsc
        nyq = size(corrs)
        allocate( ssnr(nyq) )
        do k=1,nyq
            fsc = min(abs(corrs(k)),0.999)
            ssnr(k) = (2.*fsc)/(1.-fsc)
        end do
    end function fsc2ssnr

    !> \brief  converts the FSC to the optimal low-pass filter
    function fsc2optlp( corrs ) result( filt )
        real, intent(in)  :: corrs(:) !< fsc plot (correlations)
        real, allocatable :: filt(:)  !< output filter coefficients
        integer :: nyq, k
        nyq = size(corrs)
        allocate( filt(nyq) )
        filt = 0.
        where( corrs > 0. )     filt = sqrt( 2. * corrs / (corrs + 1.) )
        where( filt  > 0.9999 ) filt = 0.99999
    end function fsc2optlp

    !> \brief  converts the SSNR to FSC
    function ssnr2fsc( ssnr ) result( corrs )
        real, intent(in)  :: ssnr(:)  !< input SSNR array
        real, allocatable :: corrs(:) !< output FSC result
        integer :: nyq, k
        nyq = size(ssnr)
        allocate( corrs(nyq) )
        do k=1,nyq
            corrs(k) = ssnr(k)/(ssnr(k)+1.)
        end do
    end function ssnr2fsc

    !> \brief  converts the SSNR 2 the optimal low-pass filter
    function ssnr2optlp( ssnr ) result( w )
        real, intent(in)  :: ssnr(:) !<  instrument SSNR
        real, allocatable :: w(:) !<  FIR low-pass filter
        integer :: nyq, k
        nyq = size(ssnr)
        allocate( w(nyq) )
        do k=1,nyq
            w(k) = ssnr(k)/(ssnr(k)+1.)
        end do
    end function ssnr2optlp

end module simple_estimate_ssnr
