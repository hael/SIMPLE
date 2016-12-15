!>  \brief  Class that defines the Contrast Transfer Function (CTF) of the electron microscope
!!
!! This class is based on a class used in CTFFIND4, developed by Alexis Rohou
!! and Nikolaus Grigorieff at Janelia Farm. The below copyright statement therefore 
!! needs to be included here:
!! Copyright 2014 Howard Hughes Medical Institute
!! All rights reserved
!! Use is subject to Janelia Farm Research Campus Software Copyright 1.1
!! license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
!!
module simple_ctf
use simple_defs  ! singleton
use simple_units ! singleton
implicit none

public :: ctf, test_ctf
private

type ctf
    private
    real    :: smpd        = 0.        
    real    :: kV          = 0.         !< acceleration voltage of the electron microscope (kilo Violts)  
    real    :: Cs          = 0.         !< spherical aberration
    integer :: Cs_unit     = millimeters
    real    :: wl          = 0.         !< wavelength
    integer :: wl_unit     = angstroms
    real    :: amp_contr   = 0.07       !< fraction of amplitude contrast (from 0. to 1., usually btw .07 and .14, see Mindell 03)
    real    :: dfx         = 0.         !< underfocus along first axis  (positive values for underfocus; larger value = more underfocus)
    real    :: dfy         = 0.         !< underfocus along second axis
    integer :: dfx_unit    = microns
    integer :: dfy_unit    = microns
    real    :: angast      = 0.         !< azimuth of first axis. 0.0 means axis is at 3 o'clock
    integer :: angast_unit = degrees
    logical :: unity       =.false.     !< a unity CTF model (4 pre-corrected particles or absence of CTF info)
    real    :: phaseq      = 0.         !< phase constrast weight
    real    :: ampliq      = 0.         !< amplitude contrast weight
  contains
    procedure, private :: init
    procedure          :: reinit
    procedure, private :: convert2PixelsAndRadians
    ! SETTERS
    procedure          :: setDefScalars
    procedure          :: setDefArray
    ! GETTERS
    procedure          :: getDefParams
    procedure          :: getDefParamsInAngstromsAndDegs
    procedure          :: getDfxInAngstroms
    procedure          :: getDfyInAngstroms
    procedure          :: getAstigAngInDegrees
    procedure          :: getAstigAngInRadians
    procedure          :: print
    procedure          :: hasPixUnits
    procedure          :: getAstig
    ! CTF EVALUATORS
    procedure          :: eval
    procedure          :: evalPhSh
    procedure          :: eval_df
    ! CTF APPLICATION ROUTINES
    procedure          :: ctf2img
    procedure          :: ctf2spec
    procedure          :: apply
    ! CALCULATORS
    procedure          :: freqOfAZero
    procedure          :: cntNrExtremaBeforeSqSpatFreq
    procedure, private :: solve4PhSh
    procedure, private :: sqFreq4PhSh
    procedure, private :: sqSf4PhSh
end type

interface ctf
    module procedure constructor
end interface

contains
    
    !>  \brief  is a constructor
    function constructor( smpd, kV, Cs, amp_contr, unity ) result( self )
        real,              intent(in) :: smpd, kV, Cs, amp_contr
        logical, optional, intent(in) :: unity
        type(ctf) :: self
        call self%init(kV, Cs, amp_contr, 0., 0., 0., smpd)
        if( present(unity) ) self%unity = unity ! unity CTF model (4 pre-corrected particles or absence of CTF info)
    end function constructor
    
    !>  \brief  initialise a CTF object (from CTFFIND4)
    subroutine init(self, kV, Cs, amp_contr, dfx, dfy, angast, smpd )
        use simple_math, only: fdim, kV2wl
        class(ctf), intent(inout) :: self      !< instance
        real,       intent(in)    :: kV        !< Kilovolts
        real,       intent(in)    :: Cs        !< mm
        real,       intent(in)    :: amp_contr !< fraction
        real,       intent(in)    :: dfx       !< um
        real,       intent(in)    :: dfy       !< um
        real,       intent(in)    :: angast    !< degrees
        real,       intent(in)    :: smpd      !< A
        self%wl          = kV2wl(kV)
        self%wl_unit     = angstroms
        self%Cs          = Cs
        self%Cs_unit     = millimeters
        self%amp_contr   = amp_contr
        self%dfx         = dfx
        self%dfx_unit    = microns
        self%dfy         = dfy
        self%dfy_unit    = microns
        self%angast      = angast
        self%angast_unit = degrees
        self%smpd        = smpd
        call self%convert2PixelsAndRadians(smpd)
        ! Compute phase and amplitude contrast weights
        self%phaseq = sqrt(1.-self%amp_contr**2.)
        self%ampliq = self%amp_contr
    end subroutine init
    
    !>  \brief  reinitialise a CTF object with new defocus/astigmatism params
    subroutine reinit( self, dfx, dfy, angast )
        class(ctf), intent(inout) :: self   !< instance
        real,       intent(in)    :: dfx    !< um
        real,       intent(in)    :: dfy    !< um
        real,       intent(in)    :: angast !< degrees
        self%dfx         = dfx
        self%dfy         = dfy
        self%dfx_unit    = microns
        self%dfy_unit    = microns
        self%angast      = angast
        self%angast_unit = degrees
        call unit_conversion(self%dfx,self%dfx_unit,pixels,self%smpd)
        call unit_conversion(self%dfy,self%dfy_unit,pixels,self%smpd)
        call unit_conversion(self%angast,self%angast_unit,radians)
    end subroutine reinit
    
    !>  \brief  Convert all values to pixel/radians unit
    subroutine convert2PixelsAndRadians( self, smpd )
        class(ctf), intent(inout) :: self
        real,       intent(in)    :: smpd
        call unit_conversion(self%Cs,self%Cs_unit,pixels,smpd)
        call unit_conversion(self%wl,self%wl_unit,pixels,smpd)
        call unit_conversion(self%dfx,self%dfx_unit,pixels,smpd)
        call unit_conversion(self%dfy,self%dfy_unit,pixels,smpd)
        call unit_conversion(self%angast,self%angast_unit,radians)
    end subroutine convert2PixelsAndRadians
    
    ! SETTERS
    
    !>  \brief  Set defocus and astigmatism in pixels/radians
    subroutine setDefScalars( self, dfx, dfy, angast )
        class(ctf), intent(inout) :: self   !<  instance
        real,       intent(in)    :: dfx    !<  pixels
        real,       intent(in)    :: dfy    !<  pixels
        real,       intent(in)    :: angast !<  radians
        self%dfx         = dfx
        self%dfx_unit    = pixels
        self%dfy         = dfy
        self%dfy_unit    = pixels
        self%angast      = angast
        self%angast_unit = radians
    end subroutine setDefScalars
    
    subroutine setDefArray( self, my_array )
        class(ctf), intent(inout) :: self
        real,       intent(in)    :: my_array(:)
        call self%setDefScalars(my_array(1),my_array(2),my_array(3))
    end subroutine setDefArray
    
    ! GETTERS
    
    function getDefParams( self ) result( defocus )
        class(ctf), intent(in) :: self
        real ::  defocus(3)
        defocus = [self%dfx,self%dfy,self%angast]
    end function getDefParams
    
    function getDefParamsInAngstromsAndDegs( self ) result( defocus )
        class(ctf), intent(in) :: self
        real :: defocus(3)
        defocus(1) = convert_unit(self%dfx,self%dfx_unit,angstroms,self%smpd)
        defocus(2) = convert_unit(self%dfy,self%dfy_unit,angstroms,self%smpd)
        defocus(3) = convert_unit(self%angast,self%angast_unit,degrees)
    end function getDefParamsInAngstromsAndDegs
    
    real function getDfxInAngstroms( self )
        class(ctf), intent(in) :: self
        getDfxInAngstroms = convert_unit(self%dfx, self%dfx_unit, angstroms, self%smpd)
    end function getDfxInAngstroms

    real function getDfyInAngstroms( self )
        class(ctf), intent(in) :: self
        getDfyInAngstroms = convert_unit(self%dfy,self%dfy_unit,angstroms,self%smpd)
    end function getDfyInAngstroms
    
    real function getAstigAngInDegrees( self ) result( azimuth )
        class(ctf), intent(in) :: self
        azimuth = convert_unit(self%angast,self%angast_unit,degrees)
        azimuth = azimuth-180.*nint(azimuth/180.)
    end function getAstigAngInDegrees
    
    real function getAstigAngInRadians( self ) result( azimuth )
        class(ctf), intent(in) :: self
        azimuth = convert_unit(self%angast,self%angast_unit,radians)
    end function getAstigAngInRadians

    subroutine print( self )
        class(ctf), intent(in) :: self
        write(*,'(a)')           '** ctf **'
        write(*,'(a,f0.3,1x,a)') 'Wavelength = ',             self%wl, unit2str(self%wl_unit)
        write(*,'(a,f0.3,1x,a)') 'Spherical aberration = ',   self%Cs, unit2str(self%Cs_unit)
        write(*,'(a,f0.3)')      'Amplitude contrast = ',     self%amp_contr
        write(*,'(a,f0.3,1x,a)') 'Defocus x = ',              self%dfx, unit2str(self%dfx_unit)
        write(*,'(a,f0.3,1x,a)') 'Defocus y = ',              self%dfy, unit2str(self%dfy_unit)
        write(*,'(a,f0.3,1x,a)') 'Azimuth of astigmatism = ', self%angast, unit2str(self%angast_unit)
    end subroutine print

    !>  \brief  to check that all units are in pixels & radians
    logical function hasPixUnits(self)
        class(ctf), intent(in) :: self
        hasPixUnits =  self%Cs_unit == pixels &
               .and.   self%wl_unit == pixels &
               .and.   self%dfx_unit == pixels &
               .and.   self%dfy_unit == pixels &
               .and.   self%angast_unit == radians
    end function hasPixUnits

    real function getAstig( self ) result( astigmatism )
        class(ctf), intent(in) :: self
        astigmatism = self%dfx-self%dfy
    end function getAstig
    
    ! CTF EVALUATORS
    
    !>  \brief Returns the CTF, based on CTFFIND3 subroutine (see Mindell 2003)
    !
    ! How to use eval:
    !
    !            lFreqSq = tfun%getLowFreq4Fit**2
    !            hFreqSq = tfun%getHighFreq4Fit**2
    !            do k=-ydim,ydim
    !                kinv = inv_logical_kdim*real(k)
    !                kinvsq = kinv*kinv
    !                do h=-xdim,xdim
    !                   hinv = inv_logical_hdim*real(h)
    !                    hinvsq = hinv*hinv
    !                    spaFreqSq = kinvsq+hinvsq
    !                    if( spaFreqSq .gt. lFreqSq .and. spaFreqSq .le. hFreqSq )then
    !                        if( spaFreqSq .gt. 0. )then
    !                            ang = atan2(k,h)
    !                        else
    !                            ang = 0.
    !                        endif
    function eval( self, spaFreqSq, dfx, dfy, angast, ang ) result( val )
        class(ctf), intent(inout) :: self      !< instance
        real,       intent(in)    :: spaFreqSq !< squared reciprocal pixels
        real,       intent(in)    :: dfx       !< Defocus along first axis (micrometers)
        real,       intent(in)    :: dfy       !< Defocus along second axis (for astigmatic CTF, dfx .ne. dfy) (micrometers)
        real,       intent(in)    :: angast    !< Azimuth of first axis. 0.0 means axis is at 3 o'clock. (radians)
        real,       intent(in)    :: ang       !< Angle at which to compute the CTF (radians)
        real :: val, phase_shift
        if( self%unity )then
            val = 1.
            return
        endif
        ! Reinitialize the CTF object, using the input parameters
        call self%reinit(dfx, dfy, angast)
        ! Compute phase shift
        phase_shift = self%evalPhSh(spaFreqSq,ang)
        ! Compute value of CTF, assuming white particles
        val = (self%phaseq*sin(phase_shift)+self%ampliq*cos(phase_shift))
    end function eval
    
    !>  \brief returns the argument (radians) to the sine and cosine terms of the ctf
    !!
    !!  We follow the convention, like the rest of the cryo-EM/3DEM field, that underfocusing the objective lens
    !!  gives rise to a positive phase shift of scattered electrons, whereas the spherical aberration gives a
    !!  negative phase shift of scattered electrons
    pure elemental real function evalPhSh( self, spaFreqSq, ang ) result( phase_shift )
        class(ctf), intent(in) :: self      !< instance
        real,       intent(in) :: spaFreqSq !< square of spatial frequency at which to compute the ctf (1/pixels^2)
        real,       intent(in) :: ang       !< angle at which to compute the ctf (radians)
        real :: df ! defocus at point at which we're evaluating the ctf
        ! compute the defocus
        df = self%eval_df(ang)
        ! compute the ctf argument
        phase_shift = pi*self%wl*spaFreqSq*(df-0.5*self%wl**2*spaFreqSq*self%Cs)
    end function evalPhSh
    
    !>  \brief  Return the effective defocus given the pre-set CTF parameters (from CTFFIND4)
    pure elemental real function eval_df( self, ang ) result (df)
        class(ctf), intent(in) :: self !< instance
        real,       intent(in) :: ang  !< angle at which to compute the defocus (radians)
        real :: defdiff, angdiff, cosarg, defsum
        defdiff = self%dfx-self%dfy
        defsum  = self%dfx+self%dfy
        angdiff = ang-self%angast
        cosarg  = 2.0*angdiff*defdiff
        if( abs(cosarg) <= TINY )then
            df = 0.5*(defsum + 1.0)         ! effective defocus
        else
            df = 0.5*(defsum + cos(cosarg)) ! effective defocus
        endif        
    end function eval_df
    
    ! CTF APPLICATION ROUTINES
    
    !>  \brief  is for making a CTF image
    !!          modes: abs, ctf, flip, flipneg, neg, square 
    subroutine ctf2img( self, img, dfx, mode, dfy, angast, bfac )
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_image, only: image
        class(ctf),       intent(inout) :: self              !< instance
        class(image),     intent(inout) :: img               !< image (output)
        real,             intent(in)    :: dfx               !< defocus x-axis
        character(len=*), intent(in)    :: mode              !< abs, ctf, flip, flipneg, neg, square 
        real, optional,   intent(in)    :: dfy, angast, bfac !< defocus y-axis, angle of astigmatism, bfactor
        integer :: lims(3,2),h,k,phys(3),hh,kk,ldim(3)
        real    :: ang,tval,ddfy,aangast,spaFreqSq,hinv,kinv,hinvsq,kinvsq,inv_ldim(3)
        if( img%is_3d() )then
            print *, 'ldim: ', img%get_ldim()
            stop 'Only 4 2D images; ctf2img; simple_ctf'
        endif
        ! set CTF params
        if( present(dfy) )then
            ddfy = dfy
        else
            ddfy = dfx
        endif
        aangast = 0.
        if( present(angast) ) aangast = angast
        ! reinit object
        call self%reinit(dfx, ddfy, aangast)
        ! initialize
        img      = cmplx(0.,0.)
        lims     = img%loop_lims(2)
        ldim     = img%get_ldim()
        inv_ldim = 1./real(ldim)
        !$omp parallel do default(shared) private(h,hinv,hinvsq,k,kinv,kinvsq,spaFreqSq,ang,tval,phys) schedule(auto)
        do h=lims(1,1),lims(1,2)
            hinv = real(h)*inv_ldim(1)
            hinvsq = hinv*hinv
            do k=lims(2,1),lims(2,2)
                kinv = real(k)*inv_ldim(2)
                kinvsq = kinv*kinv
                spaFreqSq = hinvsq+kinvsq
                ang = atan2(real(k),real(h))
                tval = self%eval(spaFreqSq, dfx, ddfy, aangast, ang)
                select case(mode)
                    case('abs')
                        tval = abs(tval)
                    case('ctf')
                        ! tval = tval
                    case('flip')
                        tval = sign(1.,tval)
                    case('flipneg')
                        tval = -sign(1.,tval)
                    case('neg')
                        tval = -tval
                    case('square')
                        tval = min(1.,max(tval**2.,0.001))
                    case DEFAULT
                        write(*,*) 'mode:', mode
                        stop 'unsupported in ctf2img; simple_ctf'
                end select
                phys = img%comp_addr_phys([h,k,0])
                call img%set_fcomp([h,k,0], cmplx(tval,0.), phys_in=phys)
            end do
        end do
        !$omp end parallel do 
        if( present(bfac) ) call img%apply_bfac(bfac)       
    end subroutine ctf2img
    
    !>  \brief  is for making a CTF spectrum (1d)
    !!          modes: abs, ctf, flip, flipneg, neg, square 
    function ctf2spec( self, kfromto, box, dfx, mode, bfac ) result(spec)
        use simple_image,  only: image
        use simple_jiffys, only: alloc_err
        class(ctf),       intent(inout) :: self       !< instance
        integer,          intent(in)    :: kfromto(2) !< Fourier index range
        integer,          intent(in)    :: box        !< box size
        real,             intent(in)    :: dfx        !< defocus
        character(len=*), intent(in)    :: mode       !< abs, ctf, flip, flipneg, neg, square 
        real, optional,   intent(in)    :: bfac       !< bfactor
        real, allocatable :: spec(:)
        integer           :: k,alloc_stat
        real              :: tval,res,wght,spaFreqSq,kinv
        ! reinit object
        call self%reinit(dfx, dfx, 0.)
        ! allocate spectrum
        allocate( spec(kfromto(1):kfromto(2)), stat=alloc_stat )
        call alloc_err('In: ctf2spec; simple_ctf', alloc_stat)
        ! do the work
        do k=kfromto(1),kfromto(2) ! loop over resolution range
            kinv = real(k)/real(box)
            spaFreqSq = kinv*kinv
            tval = self%eval(spaFreqSq, dfx, dfx, 0., 0.) ! CFT-value at k
            select case(mode)
                case('abs')
                    tval = abs(tval)
                case('ctf')
                    ! tval = tval
                case('flip')
                    tval = sign(1.,tval)
                case('flipneg')
                    tval = -sign(1.,tval)
                case('neg')
                    tval = -tval
                case('square')
                    tval = min(1.,max(tval**2.,0.001))
                case DEFAULT
                    write(*,*) 'mode:', mode
                    stop 'unsupported in ctf2img; simple_ctf'
            end select
            spec(k) = tval
        end do
        if( present(bfac) )then
            do k=kfromto(1),kfromto(2)                 ! loop over resolution range
                res = real(k)/(real(box)*self%smpd)    ! assuming square dimensions
                wght = max(0.,exp(-(bfac/4.)*res*res)) ! b-factor weight
                spec(k) = spec(k)*wght                 ! apply it
            end do
        endif      
    end function ctf2spec

    !>  \brief  is for applying CTF to an image
    subroutine apply( self, img, dfx, mode, dfy, angast, bfac )
        use simple_image, only: image
        class(ctf),       intent(inout) :: self              !< instance
        class(image),     intent(inout) :: img               !< image (output)
        real,             intent(in)    :: dfx               !< defocus x-axis
        character(len=*), intent(in)    :: mode              !< abs, ctf, flip, flipneg, neg, square
        real, optional,   intent(in)    :: dfy, angast, bfac !< defocus y-axis, angle of astigmatism, bfactor
        integer     :: ldim(3), ldim_pd(3)
        logical     :: didbwdft, ppad, didfwdft
        type(image) :: ctfimg, img_pd
        ldim = img%get_ldim()
        if( img%is_3d() )then
            print *, 'ldim: ', ldim
            stop 'Only 4 2D images; apply; simple_ctf'
        endif
        call ctfimg%new(ldim, self%smpd)
        call self%ctf2img(ctfimg, dfx, mode, dfy, angast, bfac)
        if( img%is_ft() )then
            call img%mul(ctfimg)
        else
            call img%fwd_ft
            call img%mul(ctfimg)
            call img%bwd_ft
        endif
        call ctfimg%kill
    end subroutine apply
    
    ! CALCULATORS        
    
    !>  \brief  Return the spatial frequency (reciprocal pixels) of the n-th zero, where n is given as which_zero
    real function freqOfAZero( self, dfx, dfy, angast, ang, which_zero )
        class(ctf), intent(inout) :: self       !< instance
        real,       intent(in)    :: dfx        !< Defocus along first axis (micrometers)
        real,       intent(in)    :: dfy        !< Defocus along second axis (for astigmatic CTF, dfx .ne. dfy) (micrometers)
        real,       intent(in)    :: angast     !< Azimuth of first axis. 0.0 means axis is at 3 o'clock. (radians)
        real,       intent(in)    :: ang        !< Radians. Angle at which to evaluate the CTF
        integer,    intent(in)    :: which_zero !< Number (starting at 1) of the zero you are interested in
        real, allocatable         :: phase_shifts_of_zeroes(:)
        real, allocatable         :: sq_sfs_of_zeroes(:)
        integer                   :: num_freqs
        ! Reinitialize the CTF object, using the input parameters
        call self%reinit(dfx, dfy, angast)
        ! allocate
        allocate(phase_shifts_of_zeroes(which_zero*8))
        allocate(sq_sfs_of_zeroes(2*size(phase_shifts_of_zeroes)))
        ! Solve the CTF for the phase shift - we will get lots of phase shifts, some of which will not be possible with current CTF
        call self%solve4PhSh(phase_shifts_of_zeroes,0.)
        call self%sqFreq4PhSh(phase_shifts_of_zeroes,ang,sq_sfs_of_zeroes,num_freqs)
        freqOfAZero = sqrt(sq_sfs_of_zeroes(which_zero))
    end function freqOfAZero
    
    !>  \brief  Count how many extrema the CTF goes through before reaching the given spatial frequency at the given ang
    function cntNrExtremaBeforeSqSpatFreq( self, sqSpatFreq, dfx, dfy, angast, ang ) result( number_of_extrema )
        class(ctf),intent(inout) :: self                      !< instance
        real,      intent(in)    :: sqSpatFreq                !< squared reciprocal pixels
        real,      intent(in)    :: dfx                       !< defocus along first axis (micrometers)
        real,      intent(in)    :: dfy                       !< defocus along second axis (for astigmatic CTF, dfx .ne. dfy) (micrometers)
        real,      intent(in)    :: angast                    !< azimuth of first axis. 0.0 means axis is at 3 o'clock. (radians)
        real,      intent(in)    :: ang                       !< radians. ang at which to evalulate the ctf
        integer, parameter       :: maximum_number_of_rings = 128
        real, allocatable        :: phase_shifts_of_minima(:) ! phase shifts at which minima (-1.0) occur
        real, allocatable        :: phase_shifts_of_maxima(:) ! phase shifts at which maxima (+1.0) occur
        real, allocatable        :: sq_sfs_of_minima(:)       ! squared spa freqs at which minima (-1.0) occur
        real, allocatable        :: sq_sfs_of_maxima(:)       ! squared spa freqs at which maxima (+1.0) occur
        integer                  :: num_minima, num_maxima, number_of_extrema, i
        ! Reinitialize the CTF object, using the input parameters
        call self%reinit(dfx, dfy, angast)
        ! Allocate memory for temporary results
        allocate(phase_shifts_of_minima(maximum_number_of_rings))
        allocate(phase_shifts_of_maxima(maximum_number_of_rings))
        allocate(sq_sfs_of_maxima(maximum_number_of_rings*2))
        allocate(sq_sfs_of_minima(maximum_number_of_rings*2))
        ! Solve the CTF equation for phase shifts
        call self%solve4PhSh(phase_shifts_of_minima,-1.)
        call self%solve4PhSh(phase_shifts_of_maxima,+1.)
        ! Convert phase shifts to squared spatial frequencies
        call self%sqFreq4PhSh(phase_shifts_of_minima,ang,sq_sfs_of_minima,num_minima)
        call self%sqFreq4PhSh(phase_shifts_of_maxima,ang,sq_sfs_of_maxima,num_maxima)
        ! Let's count (note that now, phase_shifts_of_minima actually store squared spatial frequencies)
        number_of_extrema = count(sq_sfs_of_minima(1:num_minima) .le. sqSpatFreq &
                            .and. sq_sfs_of_minima(1:num_minima) .gt. 0.)
        number_of_extrema = count(sq_sfs_of_maxima(1:num_maxima) .le. sqSpatFreq &
                            .and. sq_sfs_of_maxima(1:num_maxima) .gt. 0.) &
                            + number_of_extrema
        ! We get pairs of solutions for sq sp freq, so let's divide by 2
        number_of_extrema = number_of_extrema/2
    end function cntNrExtremaBeforeSqSpatFreq
    
    !>  \brief  Find the ctf argument (phase shift) values at which CTF=ctf_value
    !!
    !!  According to Wolfram Alpha, the solutions to a*cos(t) + b*sin(t) = c are:
    !!  1)  t = 2 (atan((b-sqrt(a^2+b^2-c^2))/(a+c))+pi*n),   with n an integer
    !!  2)  t = 2 (atan((b+sqrt(a^2+b^2-c^2))/(a+c))+pi*n),   with n an integer
    !!
    !!  The above two solutions only hold if:
    !!                  a+c != 0
    !!                 -b*sqrt(a**2+b**2-c**2)+a**2+ac+b**2 != 0   [for solution 1]
    !!                  b*sqrt(a**2+b**2-c**2)+a**2+ac+b**2 != 0   [for solution 2]
    !!
    !!  In our case, a = - ampl_const
    !!           and b = - sqrt(1-ampl_const**2)
    !!               c = ctf_value
    !!  and t is the "argument" to the ctf (i.e. the phase shift), which can be "converted" to a spatial frequency
    subroutine solve4PhSh( self, sols, ctf_value )
        class(ctf),        intent(in)    :: self      !< instance
        real, allocatable, intent(inout) :: sols(:)   !< solutions
        real, optional,    intent(in)    :: ctf_value !< value of the ctf at the solutions
        integer :: i,j
        real    :: a,b,c,sqrt_arg,sqrt_val,sol_first_term_1,sol_first_term_2
        if( .not. allocated(sols) ) stop 'array of solutions is not allocated; solve4PhSh; simple_ctf'
        c = 0.
        if (present(ctf_value)) c = ctf_value
        a = -self%amp_contr 
        b = -sqrt(1-self%amp_contr**2)
        ! Note that since 0.0 <= ampl_cont <= 1.0:
        ! a**2+b**2 = 1.0
        ! sqrt_arg = a**2+b**2-c**2
        sqrt_arg = max(1.0-c**2,0.0)
        sqrt_val = sqrt(sqrt_arg)
        sol_first_term_1 = atan((b-sqrt_val)/(a+c))
        sol_first_term_2 = atan((b+sqrt_val)/(a+c))
        ! Loop over the zeroes
        do i=1,size(sols)/2
            j = (i-1)*2+1
            sols(j) = 2.*(sol_first_term_1+pi*real(i-1-size(sols)/4))
            j = (i-1)*2+2
            sols(j) = 2.*(sol_first_term_2+pi*real(i-size(sols)/4))
        enddo
    end subroutine solve4PhSh
    
    !>  \brief  Compute set of squared spatial frequencies at which the given phase shift is obtained by the CTF
    subroutine sqFreq4PhSh( self, phase_shifts, ang, spaFreqSq, nsols )
        use simple_math, only: hpsort
        class(ctf), intent(in)    :: self
        real,       intent(in)    :: phase_shifts(:)
        real,       intent(in)    :: ang
        real,       intent(inout) :: spaFreqSq(:) ! assumed to be 2x size of phase_shift array
        integer,    intent(out)   :: nsols
        integer :: i
        real    :: curr_sols(2)
        integer :: curr_nsols
        nsols = 0
        do i=1,size(phase_shifts)
            call self%sqSf4PhSh(phase_shifts(i), ang, curr_sols, curr_nsols)
            if( curr_nsols .gt. 0 )then
                spaFreqSq(nsols+1:nsols+curr_nsols) = curr_sols(1:curr_nsols)
                nsols = nsols+curr_nsols
            endif
        enddo
        call hpsort(nsols,spaFreqSq(1:nsols))
    end subroutine sqFreq4PhSh
    
    !>  \brief  Given the argument to the ctf (phase shift, in radians), return the squared spatial frequency
    subroutine sqSf4PhSh( self, phase_shift, ang, sq_sf, nsols )
        class(ctf), intent(in)  :: self        !< instance
        real,       intent(in)  :: phase_shift !< phase shift (radians)
        real,       intent(in)  :: ang         !< angle at which to compute the ctf (radians)
        real,       intent(out) :: sq_sf(2)    !< squared spatial frequency (pixels^-2)
        integer,    intent(out) :: nsols       !< there may be 0, 1 or 2 solutions
        real :: df
        real :: a,b,c
        real :: det
        ! compute the defocus
        df = self%eval_df(ang)
        ! solve
        a = pi*self%wl*df
        b = 0.5*pi*self%wl**3.*self%Cs
        c = phase_shift
        det = a**2.-4.*b*c
        if( det .ge. 0.0 )then
            sq_sf(1) = (a+sqrt(det))/(-2.*b)
            sq_sf(2) = (a-sqrt(det))/(-2.*b)
            if( det .eq. 0.0 )then
                nsols = 1
            else
                nsols = 2
            endif
        else
            ! no solutions
            sq_sf = 0.
            nsols = 0
        endif
        ! squared spatial frequencies must be >=0.; let's remove bad solutions
        select case( nsols )
            case( 1 )
                if( sq_sf(1) .lt. 0. )then
                    nsols = 0
                endif
            case (2)
                if( sq_sf(2) .lt. 0. .and. sq_sf(1) .ge. 0. )then
                    nsols = 1
                else if( sq_sf(1) .lt. 0. .and. sq_sf(2) .ge. 0. )then
                    nsols = 1
                    sq_sf(1) = sq_sf(2)
                else if( sq_sf(1) .lt. 0. .and. sq_sf(2) .lt. 0. )then
                    nsols = 0
                endif
        end select
    end subroutine sqSf4PhSh

    !>  \brief  ctf class unit test
    subroutine test_ctf( nthr )
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_image,  only: image
        use simple_rnd,    only: ran3
        use simple_jiffys, only: progress
        integer, intent(in) :: nthr
        integer, parameter :: BOX_SIZES(4) = [64,128,256,512]
        real,    parameter :: SMPDS(4)     = [3.7,2.2,1.07,0.6]
        real,    parameter :: MAX_DF       = 10.
        real,    parameter :: CS           = 2.7
        real,    parameter :: FRACA        = 0.07
        real,    parameter :: KVS(3)       = [120.,200.,300.] 
        integer, parameter :: NTESTS       = 10000
        real               :: dfx, dfy, angast
        integer            :: ibox, itst, ldim(3), imode, ikv, ntot, cnt
        type(image)        :: img
        type(ctf)          :: tfun
        write(*,'(a)') '**info(simple_ctf_unit_test): testing the ctf2img routine'
!$      call omp_set_num_threads(nthr)
        ntot = size(KVS)*size(BOX_SIZES)*NTESTS
        cnt  = 0
        do ikv=1,size(KVS)
            do ibox=1,size(BOX_SIZES)
                ldim = [BOX_SIZES(ibox),BOX_SIZES(ibox),1]
                call img%new(ldim,SMPDS(ibox))
                tfun = ctf(SMPDS(ibox), KVS(ikv), CS, FRACA)
                do itst=1,NTESTS
                    cnt    = cnt + 1
                    call progress(cnt,ntot)
                    dfx    = ran3()*MAX_DF+0.1
                    dfy    = ran3()*MAX_DF+0.1
                    angast = ran3()*360.
                    call tfun%ctf2img(img, dfx, 'ctf', dfy, angast)
                end do
            end do
        end do
        write(*,'(a)') 'SIMPLE_CTF_UNIT_TEST COMPLETED SUCCESSFULLY ;-)'
    end subroutine test_ctf

end module simple_ctf
