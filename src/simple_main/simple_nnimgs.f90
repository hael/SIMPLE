module simple_nnimgs
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_jiffys            ! singleton
implicit none

public :: img_nearest_neighbors, init_nn_srch, conduct_nn_srch, kill_nn_srch
private

interface conduct_nn_srch
    module procedure :: conduct_nn_srch_1
    module procedure :: conduct_nn_srch_2
    module procedure :: conduct_nn_srch_3
end interface conduct_nn_srch

! module global variables
type(polarft_corrcalc) :: pftcc     ! polarFT corr calculator
integer                :: nnn   = 0 ! number of nearest neighbors 
integer                :: nimgs = 0 ! number of images

contains

    function img_nearest_neighbors( p, imgs ) result( nnmat )
        use simple_image,  only: image
        use simple_params, only: params
        class(params), intent(in)    :: p
        class(image),  intent(inout) :: imgs(:)
        integer, allocatable :: nnmat(:,:)
        integer              :: alloc_stat
        real                 :: nncorr
        call init_nn_srch( p, imgs )
        allocate( nnmat(nimgs,nnn), stat=alloc_stat )
        call alloc_err("In: simple_nnimgs :: img_nearest_neighbors", alloc_stat)
        call conduct_nn_srch_2( nnmat, nncorr )
        write(*,'(a,7x,f7.4)') '>>> NEAREST NEIGHBOUR CORRELATION:', nncorr
        call kill_nn_srch
    end function img_nearest_neighbors
    
    subroutine init_nn_srch( p, imgs )
        use simple_image,     only: image
        use simple_params,    only: params
        use simple_projector, only: projector
        class(params), intent(in)    :: p
        class(image),  intent(inout) :: imgs(:)
        type(projector) :: proj
        integer         :: iimg
        nimgs = size(imgs)
        nnn   = p%nnn
        call pftcc%new(nimgs, [1,nimgs], [p%box,p%box,1], p%kfromto, p%ring2, 'no')
        ! build projector 
        proj = projector(p%wfun,imgkind=p%imgkind)
        ! set the "particles" first
        do iimg=1,nimgs
            ! move to Fourier space
            call imgs(iimg)%fwd_ft
            ! transfer to polar coordinates
            call proj%img2polarft(iimg, imgs(iimg), pftcc, isptcl=.true.)
        end do
        ! copy in the references
        call pftcc%cp_ptcls2refs
    end subroutine init_nn_srch

    subroutine conduct_nn_srch_3( fromto, nnmat_part )
        integer, intent(in)  :: fromto(2)
        integer, intent(out) :: nnmat_part(fromto(1):fromto(2),nnn)
        integer :: iref, cnt, ntot
        real    :: nncorr
        write(*,'(A)') '>>> CONDUCTING NEAREST NEIGHBOUR SEARCH'
        ntot = fromto(2)-fromto(1)+1
        cnt  = 0
        do iref=fromto(1),fromto(2)
            cnt = cnt+1
            call progress(cnt, ntot)
            call conduct_nn_srch_1(iref, nnmat_part(iref,:), nncorr)
        end do
    end subroutine conduct_nn_srch_3

    subroutine conduct_nn_srch_2( nnmat, nncorr )
        integer, intent(out) :: nnmat(nimgs,nnn)
        real,    intent(out) :: nncorr
        integer :: iref
        write(*,'(A)') '>>> CONDUCTING NEAREST NEIGHBOUR SEARCH'
        do iref=1,nimgs
            call progress(iref,nimgs)
            call conduct_nn_srch_1(iref,nnmat(iref,:),nncorr)
        end do
    end subroutine conduct_nn_srch_2

    subroutine conduct_nn_srch_1( iref, nnarr, nncorr )
        use simple_math, only: corr2dist, hpsort
        integer, intent(in)  :: iref
        integer, intent(out) :: nnarr(nnn)
        real,    intent(out) :: nncorr
        integer :: iimg, inds(nimgs)
        real    :: imgcorrs(nimgs), inplcorrs(pftcc%get_nrots())
        do iimg=1,nimgs
            inds(iimg) = iimg
            if( iref == iimg )then
                imgcorrs(iimg) = 1.
            else
                inplcorrs      = pftcc%gencorrs(iref,iimg)
                imgcorrs(iimg) = maxval(inplcorrs)
            endif
        end do
        call hpsort(nimgs, imgcorrs, inds)
        nnarr(:) = inds(nimgs-nnn+1:)
        nncorr   = sum(imgcorrs(nimgs-nnn+1:))/real(nnn)
    end subroutine conduct_nn_srch_1

    subroutine kill_nn_srch
        call pftcc%kill
    end subroutine kill_nn_srch

end module simple_nnimgs