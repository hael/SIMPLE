!>  \brief Units class
!!
!! This class is based on a class used in CTFFIND4, developed by Alexis Rohou
!! and Nikolaus Grigorieff at Janelia Farm. The below copyright statement therefore 
!! needs to be included here:
!! Copyright 2014 Howard Hughes Medical Institute
!! All rights reserved
!! Use is subject to Janelia Farm Research Campus Software Copyright 1.1
!! license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
!!
module simple_units
	use simple_defs
	
	public :: unit2str, unit_conversion, convert_unit
    private
   
    integer,          parameter, public  ::  millimeters              = 0
    integer,          parameter, public  ::  microns                  = 1
    integer,          parameter, public  ::  angstroms                = 2
    integer,          parameter, public  ::  pixels                   = 3
    integer,          parameter, public  ::  degrees                  = 4
    integer,          parameter, public  ::  radians                  = 5
    integer,          parameter, public  ::  reciprocal_angstroms     = 6
    integer,          parameter, public  ::  reciprocal_pixels        = 7
    character(len=2), parameter          ::  millimeters_str          = 'mm'
    character(len=2), parameter          ::  microns_str              = 'um'
    character(len=1), parameter          ::  angstroms_str            = 'A'
    character(len=6), parameter          ::  pixels_str               = 'pixels'
    character(len=7), parameter          ::  degrees_str              = 'degrees'
    character(len=7), parameter          ::  radians_str              = 'radians'
    character(len=3), parameter          ::  reciprocal_angstroms_str = '1/A'
    character(len=8), parameter          ::  reciprocal_pixels_str    = '1/pixels'

    contains

	    !>  \brief  Return converted value
	    real elemental function convert_unit(current_value,current_unit,new_unit,smpd) result(output_value)
	        real,           intent(in) ::  current_value
	        integer,        intent(in) ::  current_unit
	        integer,        intent(in) ::  new_unit
	        real, optional, intent(in) ::  smpd
	        integer ::  temp_unit
	        ! start work
	        output_value = current_value
	        temp_unit = current_unit
	        call unit_conversion(output_value,temp_unit,new_unit,smpd)
	    end function

	    !>  \brief  Generate a converted value
	    elemental subroutine unit_conversion(current_value,current_unit,new_unit,smpd)
	        real,           intent(inout) :: current_value
	        integer,        intent(inout) :: current_unit   !< Will be changed by the routine to the desired new unit
	        integer,        intent(in)    :: new_unit
	        real, optional, intent(in)    :: smpd           !< Size of a pixel in Angstroms. May be required for the conversion.
	        real    :: ssmpd
	        logical :: conversion_success
	        logical :: smpd_missing
	        ssmpd = 1.0
	        if( present(smpd) ) ssmpd = smpd

	        conversion_success = .true.
	        smpd_missing = .false.

	        ! If we need the pixel size, check so that we have it
	        if( current_unit .ne. new_unit )then
	            if (     current_unit   .eq. pixels             &
	                .or. new_unit       .eq. pixels             &
	                .or. current_unit   .eq. reciprocal_pixels  &
	                .or. new_unit       .eq. reciprocal_pixels) then
	                ! we need to know the pixel size
	                smpd_missing = .not. present(smpd)
	                if( smpd_missing ) conversion_success = .false.
	            endif
	        endif

	        if( current_unit .ne. new_unit .and. .not. smpd_missing )then
	            select case (current_unit)
		            case( microns )
		                select case( new_unit )
			                case( pixels )
			                    ! microns to pixels
			                    current_value = current_value*1.0e4
			                    if( abs(current_value) <= 1e-6 )then
			                    	current_value = 0.0
			                    else
			                    	current_value = current_value/ssmpd
			                    endif
			                case default
			                    conversion_success = .false.
			                end select
		            case( millimeters )
		                select case( new_unit )
			                case( pixels )
			                    ! mm to pixels
			                    current_value = current_value*1.0e7/ssmpd
			                case default
			                    conversion_success = .false.
		                end select
		            case( angstroms )
		                select case( new_unit )
			                case( pixels )
			                    ! A to pixels
			                    current_value = current_value/ssmpd
			                case default
			                    conversion_success = .false.
		                end select
		            case( reciprocal_angstroms )
		                select case( new_unit )
			                case( reciprocal_pixels )
			                    ! 1/A to 1/pixels
			                    current_value = current_value*ssmpd
			                case default
			                    conversion_success = .false.
		                end select
		            case( pixels )
		                select case( new_unit )
			                case( angstroms )
			                    ! pix to A
			                    current_value = current_value*ssmpd
			                case default
			                    conversion_success = .false.
		                end select
		            case( degrees )
		                select case (new_unit)
			                case( radians )
			                    ! degrees to radians
			                    current_value = current_value/180.0e0 * pi
			                case default
			                    conversion_success = .false.
		                end select
		            case( radians )
		                select case( new_unit ) 
			                case( degrees )
			                    ! rad to deg
			                    current_value = current_value/pi*180.0e0
			                case default
			                    conversion_success = .false.
		                end select
		            case default
		                conversion_success = .false.
	            end select
	        endif
	        if (conversion_success) current_unit = new_unit
	    end subroutine

	    !>  \brief  Return a short string description for the units
	    function unit2str(unit_identifier) result(string)
	        integer,          intent(in)  ::  unit_identifier
	        character(len=:), allocatable ::  string
	        select case (unit_identifier)
	            case (millimeters)
	                string = millimeters_str
	            case (microns)
	                string = microns_str
	            case (angstroms)
	                string = angstroms_str
	            case (pixels)
	                string = pixels_str
	            case (degrees)
	                string = degrees_str
	            case (radians)
	                string = radians_str
	            case (reciprocal_angstroms)
	                string = reciprocal_angstroms_str
	            case (reciprocal_pixels)
	                string = reciprocal_pixels_str
	            case default
					stop 'unknown identifier; unit2str; simple_units'
	        end select
	    end function

end module simple_units
