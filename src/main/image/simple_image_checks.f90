submodule (simple_image) simple_image_checks
  implicit none
contains

    module pure function exists( self ) result( is )
        class(image), intent(in) :: self
        logical :: is
        is = self%existence
    end function exists

    module pure logical function is_2d(self)
        class(image), intent(in)  ::  self
        is_2d = count(self%ldim .eq. 1) .eq. 1
    end function is_2d

    module pure logical function is_3d(self)
        class(image), intent(in)  ::  self
        is_3d = .not. any(self%ldim .eq. 1)
    end function is_3d

    module pure function even_dims( self ) result( yep )
        class(image), intent(in) :: self
        logical :: yep, test(2)
        test = .false.
        test(1) = is_even(self%ldim(1))
        test(2) = is_even(self%ldim(2))
        yep = all(test)
    end function even_dims

    module pure function square_dims( self ) result( yep )
        class(image), intent(in) :: self
        logical :: yep
        yep = self%ldim(1) == self%ldim(2)
        if( self%ldim(3) == 1 .and. yep )then
        else
            yep = self%ldim(3) == self%ldim(1)
        endif
    end function square_dims

    module pure logical function same_dims_1( self1, self2 )
        class(image), intent(in) :: self1, self2
        same_dims_1 = all(self1%ldim == self2%ldim)
    end function same_dims_1

    module pure logical function same_dims( self, ldim )
        class(image), intent(in) :: self
        integer,      intent(in) :: ldim(3) !< dimensions
        same_dims = all(self%ldim == ldim)
    end function same_dims

    module logical pure function same_smpd( self1, self2 )
        class(image), intent(in) :: self1, self2
        same_smpd = abs(self1%smpd-self2%smpd) < 0.0001
    end function same_smpd

    module pure function is_ft( self ) result( is )
        class(image), intent(in) :: self
        logical :: is
        is = self%ft
    end function is_ft

    module function is_empty( self ) result( is )
        class(image), intent(in) :: self
        real    :: minmax(2)
        logical :: is
        minmax = self%minmax()
        is     = is_equal(minmax(2)-minmax(1),0.) ! empty image
    end function is_empty

end submodule simple_image_checks
