!------------------------------------------------------------------------------!
! SIMPLE v2.5         Elmlund & Elmlund Lab          simplecryoem.com          !
!------------------------------------------------------------------------------!
!> Simple module for converting atomic and volume modes
module simple_AVratios
use simple_defs
implicit none

type AVratios 
    private
    real :: Vbox=0., Vmsk=0., Vptcl=0.
    real :: Abox=0., Amsk=0., Aptcl=0.
  contains
    procedure :: Vmsk_o_Vptcl
    procedure :: Vptcl_o_Vbox
    procedure :: Vbox_o_Vptcl
    procedure :: Vmsk_o_Vbox
    procedure :: Abox_o_Amsk
    procedure :: Amsk_o_Abox
    procedure :: Abox_o_Aptcl
    procedure :: Aptcl_o_Abox
    procedure :: Amsk_o_Aptcl
end type AVratios

interface AVratios
    module procedure constructor
end interface

contains

    !> \brief  is a constructor
    function constructor( vol_ldim, img_ldim, msk, smpd, mw ) result( self )
        integer, intent(in)        :: vol_ldim(3)              !< volume logical dimensions
        integer, intent(in)        :: img_ldim(3)              !< image logical dimensions
        real, intent(in)           :: msk                      !< mask radius (pixels)
        real, intent(in)           :: smpd                     !< sampling distance (in Angstrom)
        real, intent(in), optional :: mw                       !< Molecular weight (in kDa)
        type(AVratios) :: self
        self%Vbox = real(product(vol_ldim))
        self%Abox = real(product(img_ldim))
        self%Vmsk = (4.*pi*msk*msk*msk)/3.
        self%Amsk = pi*msk*msk
        if( present(mw) )then
            self%Vptcl = (1000.*mw)/(0.81*smpd*smpd*smpd)
            self%Aptcl = pi*((3.*self%Vptcl)/(4.*pi)**(2./3.))
        else
            self%Vptcl = self%Vmsk
            self%Aptcl = self%Amsk
        endif
    end function

    !> \brief  is a getter
    real function Vmsk_o_Vptcl( self )
        class(AVratios), intent(in) :: self
        Vmsk_o_Vptcl = self%Vmsk/self%Vptcl
    end function

    !> \brief  is a getter
    real function Vptcl_o_Vbox( self )
        class(AVratios), intent(in) :: self
        Vptcl_o_Vbox = self%Vptcl/self%Vbox
    end function

    !> \brief  is a getter
    real function Vbox_o_Vptcl( self )
        class(AVratios), intent(in) :: self
        Vbox_o_Vptcl = self%Vbox/self%Vptcl
    end function

    !> \brief  is a getter
    real function Vmsk_o_Vbox( self )
        class(AVratios), intent(in) :: self
        Vmsk_o_Vbox = self%Vmsk/self%Vbox
    end function

    !> \brief  is a getter
    real function Abox_o_Amsk( self )
        class(AVratios), intent(in) :: self
        Abox_o_Amsk = self%Abox/self%Amsk
    end function

    !> \brief  is a getter
    real function Amsk_o_Abox( self )
        class(AVratios), intent(in) :: self
        Amsk_o_Abox = self%Amsk/self%Abox
    end function

    !> \brief  is a getter
    real function Abox_o_Aptcl( self )
        class(AVratios), intent(in) :: self
        Abox_o_Aptcl = self%Abox/self%Aptcl
    end function

    !> \brief  is a getter
    real function Aptcl_o_Abox( self )
        class(AVratios), intent(in) :: self
        Aptcl_o_Abox = self%Aptcl/self%Abox
    end function

    !> \brief  is a getter
    real function Amsk_o_Aptcl( self )
        class(AVratios), intent(in) :: self
        Amsk_o_Aptcl = self%Amsk/self%Aptcl
    end function

end module simple_AVratios
