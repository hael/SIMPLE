submodule (simple_image) simple_image_vis
!$ use omp_lib
!$ use omp_lib_kinds
include  'simple_lib.f08'
#include "simple_local_flags.inc"
implicit none
contains

    module subroutine print_cmat( self )
        class(image), intent(in) :: self
        write(logfhandle,*) self%cmat
    end subroutine print_cmat

    module subroutine print_rmat( self )
        class(image), intent(in) :: self
        write(logfhandle,*) self%rmat
    end subroutine print_rmat

    module subroutine vis( self, sect, geomorsphr )
        class(image),      intent(in) :: self
        integer, optional, intent(in) :: sect
        logical, optional, intent(in) :: geomorsphr !< geometrical or spherical complex format
        complex, allocatable :: fplane(:,:)
        integer              :: sect_here
        logical              :: geomorsphr_here
        sect_here = 1
        if( present(sect) ) sect_here = sect
        geomorsphr_here=.true.
        if (present(geomorsphr))geomorsphr_here=geomorsphr
        if( self%ft )then
            if( self%ldim(3) == 1 ) sect_here = 0
            fplane = self%expand_ft()
            if(geomorsphr_here)then
                call gnufor_image(real(fplane), palette='gray')
                call gnufor_image(aimag(fplane), palette='gray')
            else
                call gnufor_image(cabs(fplane), palette='gray')
                call gnufor_image(atan2(real(fplane),aimag(fplane)), palette='gray')
            endif
            deallocate(fplane)
        else
            if( self%ldim(3) == 1 ) sect_here = 1
            call gnufor_image(self%rmat(:self%ldim(1),:self%ldim(2),sect_here), palette='gray')
        endif
    end subroutine vis

    module subroutine before_after( left, right, ba, mask )
        class(image),      intent(in)    :: left, right
        type(image),       intent(inout) :: ba
        logical, optional, intent(in)    :: mask(left%ldim(1),left%ldim(2),left%ldim(3))
        integer     :: ldim(3), i, j
        if( left.eqdims.right )then
            if( left.eqsmpd.right )then
                if( left%ft .or. right%ft ) THROW_HARD('not for FTs; before_after')
                if( left%is_3d() .or. right%is_3d() ) THROW_HARD('not for 3D imgs; before_after')
                ldim = left%ldim
                ba = left
                ba%rmat(:ldim(1)/2,:ldim(2),1)   = left%rmat(:ldim(1)/2,:ldim(2),1)
                ba%rmat(ldim(1)/2+1:,:ldim(2),1) = 0.
                if( present(mask) )then
                    do i=ldim(1)/2+1,ldim(1)
                        do j=1,ldim(2)
                            if( mask(i,j,1) )then
                                ba%rmat(i,j,1) = right%rmat(i,j,1)
                            endif
                        end do
                    end do
                else
                    ba%rmat(ldim(1)/2+1:,:ldim(2),1) = right%rmat(ldim(1)/2+1:,:ldim(2),1)
                endif
                ba%rmat(ldim(1)/2:ldim(1)/2+1,:,1) = 0.
            else
                THROW_HARD('before (left) and after (right) not of same smpd; before_after')
            endif
        else
            THROW_HARD('before (left) and after (right) not of same dim; before_after')
        endif
    end subroutine before_after

    module subroutine scale_pspec4viz( self, rsmpd4viz )
        class(image),   intent(inout) :: self
        real, optional, intent(in)    :: rsmpd4viz
        type(image) :: tmp
        real        :: scale4viz, smpd4viz_here
        integer     :: box4clip
        if( self%ft )          THROW_HARD('pspec input assumed to be in real-space; scale_pspec4viz')
        if( self%ldim(3) > 1 ) THROW_HARD('pspec input assumed to be 2D; scale_pspec4viz')
        smpd4viz_here = SMPD4VIZ
        if( present(rsmpd4viz) ) smpd4viz_here = rsmpd4viz
        scale4viz = min(self%smpd / smpd4viz_here, 1.)
        if( scale4viz < 1. )then
            box4clip = round2even(scale4viz * real(self%ldim(1)))
        else
            return
        endif
        call tmp%new([box4clip,box4clip,1], smpd4viz_here)
        call self%clip(tmp)
        call tmp%zero_edgeavg
        call tmp%fft
        call tmp%pad(self)
        call self%ifft
        call tmp%kill
    end subroutine scale_pspec4viz

    module subroutine generate_orthogonal_reprojs( self, reprojs )
        class(image), intent(in)    :: self
        class(image), intent(inout) :: reprojs
        integer, parameter :: b=3
        type(image) :: reproj
        integer     :: ldim_reproj(3),ldim_reprojs(3)
        if( self%is_ft() )      THROW_HARD('Real space only; generate_orthogonal_reprojs')
        if( .not.self%is_3d() ) THROW_HARD('Volumes only; generate_orthogonal_reprojs')
        ldim_reproj(1:2) = self%ldim(1:2)
        ldim_reproj(3)   = 1
        ldim_reprojs(1)  = 3*ldim_reproj(1) + 4*b
        ldim_reprojs(2)  = ldim_reproj(2) + 2*b
        ldim_reprojs(3)  = 1
        call reproj%new(ldim_reproj,   self%smpd)
        call reprojs%new(ldim_reprojs, self%smpd)
        !$omp parallel workshare
        reproj%rmat(1:ldim_reproj(1),1:ldim_reproj(2),1) = sum(self%rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)),dim=3)
        !$omp end parallel workshare
        call reproj%norm
        reprojs%rmat(b+1:b+ldim_reproj(1),b+1:b+ldim_reproj(2),1) = reproj%rmat(1:ldim_reproj(1),1:ldim_reproj(2),1)
        !$omp parallel workshare
        reproj%rmat(1:ldim_reproj(1),1:ldim_reproj(2),1) = sum(self%rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)),dim=2)
        !$omp end parallel workshare
        call reproj%norm
        reprojs%rmat(ldim_reproj(1)+2*b+1:2*(b+ldim_reproj(1)),b+1:b+ldim_reproj(2),1) = reproj%rmat(1:ldim_reproj(1),1:ldim_reproj(2),1)
        !$omp parallel workshare
        reproj%rmat(1:ldim_reproj(1),1:ldim_reproj(2),1) = sum(self%rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)),dim=1)
        !$omp end parallel workshare
        call reproj%norm
        reprojs%rmat(2*ldim_reproj(1)+3*b+1:3*(b+ldim_reproj(1)),b+1:b+ldim_reproj(2),1) = reproj%rmat(1:ldim_reproj(1),1:ldim_reproj(2),1)
        call reproj%kill
    end subroutine generate_orthogonal_reprojs

    module subroutine collage( self1, self2, img_out )
        class(image), intent(inout) :: self1, self2, img_out
        real, parameter :: background = 128. ! taken as centre of [0.255] for jpegs
        type(image)     :: img_pad
        integer         :: ldim(3), ldim_col(3), border
        if( .not.self1%is_2d() ) THROW_HARD('2D only; collage')
        if( self1%is_ft() )      THROW_HARD('Real space only; collage')
        if( .not.self2%is_2d() ) THROW_HARD('2D only; collage')
        if( self2%is_ft() )      THROW_HARD('Real space only; collage')
        border   = 1
        ldim(1)  = max(self1%ldim(1),self2%ldim(1))
        ldim(2)  = max(self1%ldim(2),self2%ldim(2))
        ldim(1)  = max(ldim(1), ldim(2))
        ldim(2)  = ldim(1)
        ldim(3)  = 1
        ldim_col = [2*ldim(1)+border, ldim(2), 1]
        call img_out%new(ldim_col,1.)
        img_out%rmat = background
        ! pad & copy left image
        call img_pad%new(ldim,self1%get_smpd())
        img_pad%rmat = background
        call self1%norm4viz
        call self1%pad(img_pad, backgr=background)
        img_out%rmat(:ldim(1),:ldim(2),1) = img_pad%rmat(:ldim(1),:ldim(2),1)
        ! pad & copy right image
        img_pad%rmat = background
        call self2%norm4viz
        call img_pad%set_smpd(self2%get_smpd())
        call self2%pad(img_pad, backgr=background)
        img_out%rmat(ldim(1)+border+1:ldim_col(1),:ldim_col(2),1) = img_pad%rmat(:ldim(1),:ldim(2),1)
        call img_pad%kill()
    end subroutine collage

    module subroutine tile( self, stkimg, x, y)
        class(image), intent(inout) :: self
        class(image), intent(inout) :: stkimg
        integer,      intent(in)    :: x, y
        integer                     :: stkimg_ldim(3), x_start, y_start, x_end, y_end
        stkimg_ldim = stkimg%get_ldim()
        x_start = (x - 1) * stkimg_ldim(1) + 1
        y_start = (y - 1) * stkimg_ldim(2) + 1
        x_end = x * stkimg_ldim(1)
        y_end = y * stkimg_ldim(2)
        if(x_start .lt. 1 .or. y_start .lt. 1) THROW_HARD('tile: out of bounds')
        if(x_end .gt. self%ldim(1) .or. y_end .gt. self%ldim(2)) THROW_HARD('tile: out of bounds')
        call stkimg%norm4viz(brightness=80.0, maxmin=.true.)
        self%rmat(x_start:x_end, y_start:y_end, 1) = stkimg%rmat(:stkimg_ldim(1), :stkimg_ldim(2), 1)
    end subroutine tile

end submodule simple_image_vis
