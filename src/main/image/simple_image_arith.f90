submodule (simple_image) simple_image_arith
!$ use omp_lib
!$ use omp_lib_kinds
include  'simple_lib.f08'
#include "simple_local_flags.inc"
implicit none
contains

    !===============================
    ! assignment(=)
    !===============================
    module subroutine assign( selfout, selfin )
        class(image), intent(inout) :: selfout
        class(image), intent(in)    :: selfin
        call selfout%copy(selfin)
    end subroutine assign

    module subroutine assign_r2img( self, realin )
        class(image), intent(inout) :: self
        real,         intent(in)    :: realin
        self%rmat = realin
        self%ft = .false.
    end subroutine assign_r2img

    module subroutine assign_c2img( self, compin )
        class(image), intent(inout) :: self
        complex,      intent(in)    :: compin
        self%cmat = compin
        self%ft = .true.
    end subroutine assign_c2img

    !===============================
    ! operator(-)
    !===============================
    module function subtraction( self_from, self_to ) result( self )
        class(image), intent(in) :: self_from, self_to
        type(image) :: self
        if( self_from.eqdims.self_to )then
            call self%new(self_from%ldim, self_from%smpd)
            if( self_from%ft .neqv. self_to%ft )then
                THROW_HARD('cannot subtract images of different FT state; subtraction(+)')
            endif
            self%rmat = self_from%rmat-self_to%rmat
        else
            THROW_HARD('cannot subtract images of different dims; subtraction(-)')
        endif
    end function subtraction

    !===============================
    ! operator(+)
    !===============================
    module function addition( self1, self2 ) result( self )
        class(image), intent(in) :: self1, self2
        type(image) :: self
        if( self1.eqdims.self2 )then
            call self%new(self1%ldim, self1%smpd)
            if( self1%ft .neqv. self2%ft )then
                THROW_HARD('cannot add images of different FT state; addition(+)')
            endif
            self%rmat = self1%rmat+self2%rmat
        else
            THROW_HARD('cannot add images of different dims; addition(+)')
        endif
        self%ft = self1%ft
    end function addition

    module function addition_const_real( self1, rconst ) result( self )
        class(image), intent(in) :: self1
        real,         intent(in) :: rconst
        type(image) :: self
        call self%new(self1%ldim, self1%smpd)
        self%rmat = self1%rmat+rconst
        self%ft = self1%ft
    end function addition_const_real

    !===============================
    ! operator(*)
    !===============================
    module function multiplication( self1, self2 ) result( self )
        class(image), intent(in) :: self1, self2
        type(image) :: self
        if( self1.eqdims.self2 )then
            call self%new(self1%ldim, self1%smpd)
            if( self1%ft .and. self2%ft )then
                self%cmat = self1%cmat*self2%cmat
                self%ft = .true.
            else if( self1%ft .eqv. self2%ft )then
                self%rmat = self1%rmat*self2%rmat
                self%ft = .false.
            else if(self1%ft)then
                self%cmat = self1%cmat*self2%rmat
                self%ft = .true.
            else
                self%cmat = self1%rmat*self2%cmat
                self%ft = .true.
            endif
        else
            THROW_HARD('cannot multiply images of different dims; multiplication(*)')
        endif
    end function multiplication

    module function multiplication_const_real( self1, rconst ) result( self )
        class(image), intent(in) :: self1
        real,         intent(in) :: rconst
        type(image) :: self
        call self%new(self1%ldim, self1%smpd)
        self%ft = self1%ft
        if(self1%ft)then
            self%cmat = self1%cmat*rconst
        else
            self%rmat = self1%rmat*rconst
        endif
    end function multiplication_const_real

    module function multiplication_const_int( self1, iconst ) result( self )
        class(image), intent(in) :: self1
        integer,      intent(in) :: iconst
        type(image) :: self
        call self%new(self1%ldim, self1%smpd)
        self%ft = self1%ft
        if(self1%ft)then
            self%cmat = self1%cmat*iconst
        else
            self%rmat = self1%rmat*iconst
        endif
    end function multiplication_const_int

    !===============================
    ! operator(/)
    !===============================
    module function division( self1, self2 ) result( self )
        class(image), intent(in) :: self1, self2
        type(image) :: self
        if( self1.eqdims.self2 )then
            call self%new(self1%ldim, self1%smpd)
            if( self1%ft .and. self2%ft )then
                self%cmat = self1%cmat/self2%cmat
                self%ft = .true.
            else if( self1%ft .eqv. self2%ft )then
                self%rmat = self1%rmat/self2%rmat
                self%ft = .false.
            else if(self1%ft)then
                self%cmat = self1%cmat/self2%rmat
                self%ft = .true.
            else
                self%cmat = self1%rmat/self2%cmat
                self%ft = .true.
            endif
        else
            THROW_HARD('cannot divide images of different dims; division(/)')
        endif
    end function division

    !===============================
    ! add_* family
    !===============================
    module subroutine add_1( self, self_to_add, w )
        class(image),   intent(inout) :: self
        class(image),   intent(in)    :: self_to_add
        real, optional, intent(in)    :: w
        real :: ww
        ww = 1.
        if( present(w) ) ww = w
        if( self%ft )then
            self%cmat = self%cmat+ww*self_to_add%cmat
        else
            self%rmat = self%rmat+ww*self_to_add%rmat
        endif
    end subroutine add_1

    module subroutine add_2( self, logi, comp, phys_in, phys_out )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: logi(3)
        complex,           intent(in)    :: comp
        integer, optional, intent(in)   :: phys_in(3)
        integer, optional, intent(out)   :: phys_out(3)
        integer :: phys(3)
        complex :: comp_here
        if( .not. self%ft ) THROW_HARD('cannot add complex number to real image; add_2')
        if( present(phys_in) )then
            phys = phys_in
        else
            phys = self%fit%comp_addr_phys(logi)
        endif
        if( logi(1) < 0 )then
            comp_here = conjg(comp)
        else
            comp_here = comp
        endif
        self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3))+comp_here
        if( present(phys_out) ) phys_out = phys
    end subroutine add_2

    module subroutine add_3( self, rcomp, i, j, k )
        class(image), intent(inout) :: self
        real,         intent(in)    :: rcomp
        integer,      intent(in)    :: i, j, k
        if(  self%ft ) THROW_HARD('cannot add real number to transform; add_3')
        self%rmat(i,j,k) = self%rmat(i,j,k)+rcomp
    end subroutine add_3

    module subroutine add_4( self, logi, comp, w, k )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: logi(3)
        complex,      intent(in)    :: comp
        real,         intent(in)    :: w, k(:,:,:)
        integer :: phys(3)
        complex :: comp_here
        if( .not. self%ft ) THROW_HARD('cannot add complex number to real image; add_2')
        phys = self%fit%comp_addr_phys(logi)
        if( logi(1) < 0 )then
            comp_here = conjg(comp)
        else
            comp_here = comp
        endif
        if( abs(k(phys(1),phys(2),phys(3))) > 1e-6 )then
            self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3))+(comp_here/k(phys(1),phys(2),phys(3)))*w
        endif
    end subroutine add_4

    module subroutine add_5( self, c )
        class(image), intent(inout) :: self
        real,         intent(in)    :: c
        if( self%ft )then
            self%cmat = self%cmat+cmplx(c,0.)
        else
            self%rmat = self%rmat+c
        endif
    end subroutine add_5

    module subroutine add_workshare( self, self_to_add )
        class(image),   intent(inout) :: self
        class(image),   intent(in)    :: self_to_add
        if( self%ft )then
            !$omp parallel workshare proc_bind(close)
            self%cmat = self%cmat + self_to_add%cmat
            !$omp end parallel workshare
        else
            !$omp parallel workshare proc_bind(close)
            self%rmat = self%rmat + self_to_add%rmat
            !$omp end parallel workshare
        endif
    end subroutine add_workshare

    !===============================
    ! subtr_* family
    !===============================
    module subroutine subtr_1( self, self_to_subtr, w )
        class(image),   intent(inout) :: self
        class(image),   intent(in)    :: self_to_subtr
        real, optional, intent(in)    :: w
        real :: ww
        ww = 1.0
        if( present(w) ) ww = w
        if( self%ft )then
            self%cmat = self%cmat-ww*self_to_subtr%cmat
        else
            self%rmat = self%rmat-ww*self_to_subtr%rmat
        endif
    end subroutine subtr_1

    module subroutine subtr_2( self, logi, comp, phys_in, phys_out )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: logi(3)
        complex,           intent(in)    :: comp
        integer, optional, intent(out)   :: phys_in(3)
        integer, optional, intent(out)   :: phys_out(3)
        integer :: phys(3)
        complex :: comp_here
        if( .not. self%ft ) THROW_HARD('cannot subtract complex number from real image; subtr_2')
        if( present(phys_in) )then
            phys = phys_in
        else
            phys = self%fit%comp_addr_phys(logi)
        endif
        if( logi(1) < 0 )then
            comp_here = conjg(comp)
        else
            comp_here = comp
        endif
        self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3))-comp_here
        if( present(phys_out) ) phys_out = phys
    end subroutine subtr_2

    module subroutine subtr_3( self, logi, comp, w, k )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: logi(3)
        complex,      intent(in)    :: comp
        real,         intent(in)    :: w, k(:,:,:)
        integer :: phys(3)
        complex :: comp_here
        if( .not. self%ft ) THROW_HARD('cannot subtract complex number from real image; subtr_3')
        phys = self%fit%comp_addr_phys(logi)
        if( logi(1) < 0 )then
            comp_here = conjg(comp)
        else
            comp_here = comp
        endif
        if( .not. is_zero(k(phys(1),phys(2),phys(3))) )then
            self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3)) -&
                &(comp_here/k(phys(1),phys(2),phys(3)))*w
        endif
    end subroutine subtr_3

    module subroutine subtr_4( self, c )
        class(image), intent(inout) :: self
        real,         intent(in)    :: c
        self%rmat = self%rmat-c
    end subroutine subtr_4

    !===============================
    ! div_* family
    !===============================
    module subroutine div_1( self, c )
        class(image), intent(inout) :: self
        real,         intent(in)    :: c
        if( abs(c) < 1e-6 )then
            THROW_HARD('division with zero; div_1')
        else
            if( self%ft )then
                self%cmat = self%cmat/c
            else
                self%rmat = self%rmat/c
            endif
        endif
    end subroutine div_1

    module subroutine div_2( self, logi, k, square )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: logi(3)
        real,         intent(in)    :: k(:,:,:)
        logical,      intent(in)    :: square
        integer :: phys(3)
        if( self%ft )then
            phys = self%fit%comp_addr_phys(logi)
            if( abs(k(phys(1),phys(2),phys(3))) > 1e-6 )then
                self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3))/k(phys(1),phys(2),phys(3))
            else
                self%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.) ! this is desirable for kernel division
            endif
        else
            THROW_HARD('Image need to be Fourier transformed; div_2')
        endif
    end subroutine div_2

    module subroutine div_3( self, logi, k, phys_in )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: logi(3)
        real,              intent(in)    :: k
        integer, optional, intent(in)    :: phys_in(3)
        integer :: phys(3)
        if( self%ft )then
            if( present(phys_in) )then
                phys = phys_in
            else
                phys = self%fit%comp_addr_phys(logi)
            endif
            if( abs(k) > 1e-6 )then
                self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3))/k
            else
                self%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.) ! this is desirable for kernel division
            endif
        else
            if( abs(k) > 1e-6 )then
                self%rmat(logi(1),logi(2),logi(3)) = self%rmat(logi(1),logi(2),logi(3))/k
            endif
        endif
    end subroutine div_3

    module subroutine div_4( self, self2div )
        class(image), intent(inout) :: self
        class(image), intent(in)    :: self2div
        if( self.eqdims.self2div )then
            if( self%ft .and. self2div%ft )then
                self%cmat = self%cmat/self2div%cmat
            else if( self%ft .eqv. self2div%ft )then
                where(abs(self2div%rmat) > 1.e-6) self%rmat = self%rmat/self2div%rmat
                self%ft = .false.
            else if(self%ft)then
              where(abs(self2div%rmat) > 1.e-6) self%cmat = self%cmat/self2div%rmat
            else
                self%cmat = self%rmat/self2div%cmat
                self%ft = .true.
            endif
        else
            THROW_HARD('cannot divide images of different dims; div_4')
        endif
    end subroutine div_4

    !===============================
    ! mul_* family
    !===============================
    module subroutine mul_1( self, logi, rc, phys_in, phys_out )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: logi(3)
        real,              intent(in)    :: rc
        integer, optional, intent(in)    :: phys_in(3)
        integer, optional, intent(out)   :: phys_out(3)
        integer :: phys(3)
        if( self%is_ft() )then
            if( present(phys_in) )then
                phys = phys_in
            else
                phys = self%fit%comp_addr_phys(logi)
            endif
            if( present(phys_out) ) phys_out = phys
            self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3))*rc
        else
            self%rmat(logi(1),logi(2),logi(3)) = self%rmat(logi(1),logi(2),logi(3))*rc
        endif
    end subroutine mul_1

    module subroutine mul_2( self, rc )
        class(image), intent(inout) :: self
        real,         intent(in)    :: rc
        if( self%ft )then
            self%cmat = self%cmat*rc
        else
            self%rmat = self%rmat*rc
        endif
    end subroutine mul_2

    module subroutine mul_3( self, self2mul )
        class(image), intent(inout) :: self
        class(image), intent(in)    :: self2mul
        if( self.eqdims.self2mul )then
            if( self%ft .and. self2mul%ft )then
                self%cmat = self%cmat*self2mul%cmat
            else if( self%ft .eqv. self2mul%ft )then
                self%rmat = self%rmat*self2mul%rmat
                self%ft = .false.
            else if(self%ft)then
                self%cmat = self%cmat*self2mul%rmat
            else
                self%cmat = self%rmat*self2mul%cmat
                self%ft = .true.
            endif
        else
            THROW_HARD('cannot multiply images of different dims; mul_3')
        endif
    end subroutine mul_3

    module subroutine mul_4( self, self2mul, lp )
        class(image), intent(inout) :: self
        class(image), intent(in)    :: self2mul
        real,         intent(in)    :: lp
        integer                     :: lims(3,2),sqlim,h,k,l,phys(3)
        if( .not. self%is_ft() )     THROW_HARD('low-pass limited multiplication requires self to be FT')
        if( .not. self2mul%is_ft() ) THROW_HARD('low-pass limited multiplication requires self2mul to be FT')
        if( self.eqdims.self2mul )then
            lims = self%fit%loop_lims(1,lp)
            sqlim = (maxval(lims(:,2)))**2
            do h=lims(1,1),lims(1,2)
                do k=lims(2,1),lims(2,2)
                    do l=lims(3,1),lims(3,2)
                        if( h * h + k * k + l * l <= sqlim )then
                            phys = self%fit%comp_addr_phys([h,k,l])
                            self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3))*&
                                self2mul%cmat(phys(1),phys(2),phys(3))
                        endif
                    end do
                end do
            end do
        else
            THROW_HARD('cannot multiply images of different dims; mul_3')
        endif
    end subroutine mul_4

    module subroutine mul_5( self, logi, c)
        class(image), intent(inout) :: self
        integer,      intent(in)    :: logi(3)
        complex,      intent(in)    :: c
        integer :: phys(3)
        phys = self%fit%comp_addr_phys(logi)
        self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3)) * c
    end subroutine mul_5

    module function conjugate( self ) result ( self_out )
        class(image), intent(in) :: self
        type(image) :: self_out
        if( self%ft )then
            call self_out%copy(self)
            self_out%cmat = conjg(self%cmat)
        else
            THROW_WARN('cannot conjugate real image')
        endif
    end function conjugate

    !===============================
    ! index-based per-element ops
    !===============================
    module subroutine mul_rmat_at_1( self, logi, rval )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: logi(3)
        real,         intent(in)    :: rval
        self%rmat(logi(1),logi(2),logi(3)) = self%rmat(logi(1),logi(2),logi(3))*rval
    end subroutine mul_rmat_at_1

    module elemental pure subroutine mul_rmat_at_2( self,i, j, k, rval )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: i,j,k
        real,         intent(in)    :: rval
        self%rmat(i,j,k) = self%rmat(i,j,k)*rval
    end subroutine mul_rmat_at_2

    module subroutine div_rmat_at_1( self, logi, rval )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: logi(3)
        real,         intent(in)    :: rval
        if( abs(rval) > 1e-6 )then
            self%rmat(logi(1),logi(2),logi(3)) = self%rmat(logi(1),logi(2),logi(3))/rval
        endif
    end subroutine div_rmat_at_1

    module subroutine div_rmat_at_2( self, i, j, k, rval )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: i,j,k
        real,         intent(in)    :: rval
        if( abs(rval) > 1e-6 )then
            self%rmat(i,j,k) = self%rmat(i,j,k)/rval
        endif
    end subroutine div_rmat_at_2

    module pure subroutine add_cmat_at_1( self , phys , comp)
        class(image), intent(inout) :: self
        integer,      intent(in)    :: phys(3)
        complex,      intent(in)    :: comp
        self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3)) + comp
    end subroutine add_cmat_at_1

    module pure subroutine add_cmat_at_2( self, h, k, l, comp)
        class(image), intent(inout) :: self
        integer,      intent(in) :: h,k,l
        complex,      intent(in) :: comp
        self%cmat(h,k,l) = self%cmat(h,k,l) + comp
    end subroutine add_cmat_at_2

    module pure subroutine mul_cmat_at_1( self, phys, rval)
        class(image), intent(inout) :: self
        integer,      intent(in)    :: phys(3)
        real,         intent(in)    :: rval
        self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3)) * rval
    end subroutine mul_cmat_at_1

    module pure subroutine mul_cmat_at_2( self, phys, cval )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: phys(3)
        complex,      intent(in)    :: cval
        self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3)) * cval
    end subroutine mul_cmat_at_2

    module pure subroutine mul_cmat_at_3( self, h,k,l,rval )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: h,k,l
        real,         intent(in)    :: rval
        self%cmat(h,k,l) = self%cmat(h,k,l) * rval
    end subroutine mul_cmat_at_3

    module pure subroutine mul_cmat_at_4( self, h,k,l, cval )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: h,k,l
        complex,      intent(in)    :: cval
        self%cmat(h,k,l) = self%cmat(h,k,l) * cval
    end subroutine mul_cmat_at_4

    module subroutine mul_cmat_1( self, rmat )
        class(image), intent(inout) :: self
        real,         intent(in)    :: rmat(self%array_shape(1),self%array_shape(2),self%array_shape(3))
        self%cmat = self%cmat * rmat
    end subroutine mul_cmat_1

    module subroutine mul_cmat_2( self, rmat, resmsk )
        class(image), intent(inout) :: self
        real,         intent(in)    :: rmat(self%array_shape(1),self%array_shape(2),self%array_shape(3))
        logical,      intent(in)    :: resmsk(self%array_shape(1),self%array_shape(2),self%array_shape(3))
        where( resmsk )
            self%cmat = self%cmat * rmat
        end where
    end subroutine mul_cmat_2

    module pure subroutine div_cmat_at_1( self, phys, rval )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: phys(3)
        real,              intent(in)    :: rval
        if( abs(rval) > 1.e-6 )then
            self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3)) / rval
        else
            self%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.) ! this is desirable for kernel division
        endif
    end subroutine div_cmat_at_1

    module pure subroutine div_cmat_at_2( self,h,k,l, rval)
        class(image), intent(inout) :: self
        integer,      intent(in) :: h,k,l
        real,         intent(in) :: rval
        if( abs(rval) > 1.e-6 )then
            self%cmat(h,k,l) = self%cmat(h,k,l) / rval
        else
            self%cmat(h,k,l) =cmplx(0.,0.)
        end if
    end subroutine div_cmat_at_2

    module pure subroutine div_cmat_at_3( self, phys, cval )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: phys(3)
        complex,           intent(in)    :: cval
        if( abs(cval) > 1.e-6 )then
            self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3)) / cval
        else
            self%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.) ! this is desirable for kernel division
        endif
    end subroutine div_cmat_at_3

    module pure subroutine div_cmat_at_4( self,h,k,l, cval)
        class(image), intent(inout) :: self
        integer,      intent(in)    :: h,k,l
        complex,      intent(in)    :: cval
        if( abs(cval) > 1.e-6 )then
            self%cmat(h,k,l) = self%cmat(h,k,l) / cval
        else
            self%cmat(h,k,l) = cmplx(0.,0.)
        end if
    end subroutine div_cmat_at_4

    !===============================
    ! sq_rt
    !===============================
    module subroutine sq_rt( self )
        class(image), intent(inout) :: self
        if( self%ft ) THROW_HARD('Real space only; sq_rt')
        where( self%rmat > 0. )
            self%rmat = sqrt(self%rmat)
        else where
            self%rmat = 0.
        end where
    end subroutine sq_rt

    !>  adds complex matrices from images & arrays. Specialized routine for simple_classaverager
    module subroutine add_cmats_to_cmats( self1 , self2 , self3, self4, self2set1, self2set2, lims, expcmat3, expcmat4)
        class(image), intent(in)    :: self1, self2,self3,self4
        class(image), intent(inout) :: self2set1, self2set2
        integer,      intent(in)    :: lims(3,2)
        real,         intent(inout) :: expcmat3(lims(1,1):lims(1,2),lims(2,1):lims(2,2))
        real,         intent(inout) :: expcmat4(lims(1,1):lims(1,2),lims(2,1):lims(2,2))
        integer :: h, k, logi(3), phys(3)
        !$omp parallel default(shared) private(h,k,logi,phys) proc_bind(close)
        !$omp workshare
        self1%cmat = self1%cmat + self2set1%cmat
        self2%cmat = self2%cmat + self2set2%cmat
        !$omp end workshare
        !$omp do collapse(2) schedule(static)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                logi = [h,k,0]
                phys = self1%comp_addr_phys(logi)
                self3%cmat(phys(1),phys(2),phys(3)) = self3%cmat(phys(1),phys(2),phys(3)) + cmplx(expcmat3(h,k),0.)
                self4%cmat(phys(1),phys(2),phys(3)) = self4%cmat(phys(1),phys(2),phys(3)) + cmplx(expcmat4(h,k),0.)
            enddo
        enddo
        !$omp end do nowait
        !$omp end parallel
    end subroutine add_cmats_to_cmats

end submodule simple_image_arith
