module simple_masker
implicit none

public  :: automask, automask2D
private

interface automask
    module procedure automask_1
    module procedure automask_2
    module procedure automask_3
    module procedure automask_4
end interface

contains
    
    !>  \brief  is for generating a mask for solvent flattening of an image
    subroutine automask_1(img, p, img_msk, nvox )
        use simple_image,  only: image
        use simple_params, only: params
        use simple_math    ! use all in there
        class(image), intent(in)      :: img
        class(params), intent(in)     :: p
        class(image), intent(inout)   :: img_msk
        integer, intent(in), optional :: nvox
        type(image)                   :: img_tmp
        integer                       :: i, nnvox
        call img_msk%copy(img)        ! make a copy of the image
        if( img_msk%is_ft() )then     ! make sure that msk image is real 
            call img_msk%bwd_ft
        endif
        call img_msk%bp(0., p%amsklp) ! low-pass filter the mask image
        if( present(nvox) )then
            nnvox = nvox
        else
            ! find nr of voxels corresponding to mw
            if( p%dens > 0. )then
                nnvox = nvoxfind(p%smpd, p%mw, p%dens)
            else
                nnvox = nvoxfind(p%smpd, p%mw)     
            endif
            nnvox = nint(1.1*real(nnvox))     ! this is to compensate for the low-pass filtering
        endif
        call img_msk%bin(nnvox)               ! binarize
        if( present(nvox) )then
            ! todo
        else
            ! binary layer growth
            do i=1,p%binwidth
                call img_msk%grow_bin
            end do
        endif
        call img_msk%cos_edge(p%edge) ! real-space filter based softening of the edge
    end subroutine automask_1
    
    !>  \brief  is for generating and applying a mask for 
    !!          solvent flattening of an image
    subroutine automask_2(img, p, nvox )
        use simple_image,  only: image
        use simple_params, only: params
        class(image), intent(inout)   :: img
        class(params), intent(in)     :: p
        integer, intent(in), optional :: nvox
        logical                       :: didft
        type(image)                   :: mask
        didft = .false.
        if( img%is_ft() )then
            call img%bwd_ft
            didft = .true.
        endif
        call automask_1(img, p, mask, nvox)
        call img%mul(mask)
        call mask%kill
        if( didft ) call img%fwd_ft
    end subroutine automask_2
    
    !>  \brief  is for generating, applying, and writing to file 
    !!           a mask for solvent flattening of an image
    subroutine automask_3( b, p, cline, recvol, maskvol, volnam, masknam )
        use simple_build,   only: build
        use simple_params,  only: params
        use simple_cmdline, only: cmdline
        use simple_image,   only: image
        class(build),   intent(inout) :: b
        class(params),  intent(in)    :: p
        class(cmdline), intent(inout) :: cline
        class(image),   intent(inout) :: recvol, maskvol
        character(len=*), intent(in)  :: volnam, masknam
        if( .not. recvol%exists() )  stop 'recvol not allocated; automask_3; simple_masker'
        if( .not. maskvol%exists() ) stop 'maskvol not allocated; automask_3; simple_masker'
        call automask_4( b, p, cline, recvol, maskvol )
        if( cline%defined('nvox') )then
            call maskvol%write(masknam, del_if_exists=.true.)
        else if( cline%defined('mw') )then
            call maskvol%write(masknam, del_if_exists=.true.)
        else if( cline%defined('mskfile') )then
            call maskvol%write(p%mskfile, del_if_exists=.true.)
        endif
        call recvol%write(volnam, del_if_exists=.true.)
    end subroutine automask_3
    
    !>  \brief  is for generating & applying a mask for solvent flattening of an image
    subroutine automask_4( b, p, cline, recvol, maskvol )
        use simple_build,   only: build
        use simple_params,  only: params
        use simple_cmdline, only: cmdline
        use simple_image,   only: image
        class(build),   intent(inout) :: b
        class(params),  intent(in)    :: p
        class(cmdline), intent(inout) :: cline
        class(image), intent(inout)   :: recvol, maskvol
        if( .not. recvol%exists() )  stop 'recvol not allocated; automask_3; simple_masker'
        if( .not. maskvol%exists() ) stop 'maskvol not allocated; automask_3; simple_masker'
        if( cline%defined('nvox') )then
            call automask_1(recvol, p, maskvol, p%nvox)
            call maskvol%mask(p%msk,'soft')
            call recvol%mul(maskvol)         
        else if( cline%defined('mw') )then
            call automask_1(recvol, p, maskvol)
            call maskvol%mask(p%msk,'soft')
            call recvol%mul(maskvol)
        else if( cline%defined('mskfile') )then
            call maskvol%mask(p%msk,'soft')
            call recvol%mul(maskvol) 
        endif
    end subroutine automask_4
    
    !>  \brief  is for automasking in 2D
    subroutine automask2D( img, p, img_msk_out )
        use simple_params,  only: params
        use simple_image,   only: image
        class(image),           intent(inout) :: img
        class(params),          intent(in)    :: p
        class(image), optional, intent(out)   :: img_msk_out
        type(image) :: img_pad, img_msk, img_copy
        integer     :: ldim(3), ldim_pad(3)
        real        :: smpd, rslask
        ldim = img%get_ldim()
        smpd = img%get_smpd()
        if( ldim(3) > 1 ) stop 'not for 3D images; simple_masker :: automask_5'
        ldim_pad(1:2) = ldim(1:2)*2
        ldim_pad(3)   = 1
        call img_pad%new(ldim_pad, smpd)
        call img_msk%new(ldim,     smpd)
        img_copy = img
        call img_copy%norm
        call img_copy%mask(p%msk, 'soft')
        call img_copy%pad(img_pad)
        call img_pad%fwd_ft
        call img_pad%bp(0., p%amsklp)
        call img_pad%bwd_ft
        call img_pad%bin
        call img_pad%grow_bin
        call img_pad%cos_edge(p%edge)
        call img_pad%clip(img_msk)
        call img_msk%norm('sigm')
        call img%mul(img_msk)
        if( present(img_msk_out) ) img_msk_out = img_msk
        call img_copy%kill
        call img_pad%kill
        call img_msk%kill
    end subroutine automask2D

end module simple_masker
