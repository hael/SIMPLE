!! Dictionary of STAR format
module simple_star_dict
include 'simple_lib.f08'
implicit none

type star_dict
    type(chash),              public :: keywords_filename
    type(chash),              public :: keywords_general
    type(chash),              public :: keywords_class3D

    type(hash),               public :: dict
contains
    procedure :: new
    procedure :: exists
    procedure :: kill
end type star_dict

interface star_dict
    module procedure constructor
end interface

contains


    !>  \brief  is a constructor
    function constructor( ) result( self )
        type(ori) :: self
        call self%new_star_dict
    end function constructor

    !>  \brief  is a constructor
    subroutine new_star_dict( self )
        class(star_dict), intent(inout) :: self
        call self%kill
        self%htab = hash()

        !! File name keywords and their formats
        self%keywords_filename = chash()
        call self%keywords_filename%new(5)
        call self%keywords_filename%push('ImageName',          '([0-9]{05})@(Extract|Polish)/job[0-9]{03}/filename.mrc[s]?')
        call self%keywords_filename%push('MicrographName',     'MotionCorr/job[0-9]{03}/filename.mrc[s]?')
        call self%keywords_filename%push('MicrographNameNoDW', 'MotionCorr/job[0-9]{03}/filename_noDW.mrc[s]?')
        call self%keywords_filename%push('OriginalParticleName','Extract/job[0-9]{03}/filename.mrc[s]?')
        call self%keywords_filename%push('CtfImage',            'Ctffind/job[0-9]{03}/filename.ctf:mrc')



        !! General keywords
        self%keywords_general = chash()
        call self%keywords_general%new(33)
        call self%keywords_general%push('AmplitudeContrast',                 'fraca')   !!cmd_dict/params::fraca
        call self%keywords_general%push('AnglePsi',                          'e3')   !! imghead: gamma or psi2
        call self%keywords_general%push('AngleRot',                          'e1')   !! imghead: theta
        call self%keywords_general%push('AngleTilt',                         'e2')     !! imghead phi
        call self%keywords_general%push('AutopickFigureOfMerit',             '?')   !!
        call self%keywords_general%push('AverageNrOfFrames',                 '?')   !!
        call self%keywords_general%push('ClassNumber',                       '?')   !!
        call self%keywords_general%push('CoordinateX',                       'x')   !!  imghead
        call self%keywords_general%push('CoordinateY',                       'y')   !!  imghead
        call self%keywords_general%push('CtfBfactor',                        '?')     !! cmd_dict
        call self%keywords_general%push('CtfFigureOfMerit',                  '?')     !!
        call self%keywords_general%push('CtfScalefactor',                    'scale') !!cmd_dict/
        call self%keywords_general%push('DefocusAngle',                      'angast')
        call self%keywords_general%push('DefocusU',                          'dfx')  !! dfmax
        call self%keywords_general%push('DefocusV',                          'dfy') !! dfmin
        call self%keywords_general%push('DetectorPixelSize',                 'smpd') !! ??
        call self%keywords_general%push('GroupName',                         '')
        call self%keywords_general%push('GroupNumber',                       '')
        call self%keywords_general%push('LogLikeliContribution',             '')
        call self%keywords_general%push('Magnification',                     '')
        call self%keywords_general%push('MaxValueProbDistribution',          '')
        call self%keywords_general%push('MovieFramesRunningAverage',         '')
        call self%keywords_general%push('NormCorrection',                    'norm')  !! cmd_dict::norm
        call self%keywords_general%push('NrOfFrames',                        '')
        call self%keywords_general%push('NrOfSignificantSamples',            '?')
        call self%keywords_general%push('OriginX',                           'x')
        call self%keywords_general%push('OriginY',                           'y')
        call self%keywords_general%push('PhaseShift',                        '?')
        call self%keywords_general%push('RandomSubset',                      '?')
        call self%keywords_general%push('SphericalAberration',               'cs') !! cmd_dict
        call self%keywords_general%push('Voltage',                           'kv') !! cmd_dict


!! Class3D data_model general
        call self%keywords_class3D%push('ReferenceDimensionality ',          '')
        call self%keywords_class3D%push('DataDimensionality ',               '')
        call self%keywords_class3D%push('OriginalImageSize',                 '')
        call self%keywords_class3D%push('CurrentResolution',                 '')
        call self%keywords_class3D%push('CurrentImageSize',                  '')
        call self%keywords_class3D%push('PaddingFactor',                     '')
        call self%keywords_class3D%push('IsHelix',                           '')
        call self%keywords_class3D%push('FourierSpaceInterpolator',          '')
        call self%keywords_class3D%push('MinRadiusNnInterpolation',          '')
        call self%keywords_class3D%push('PixelSize',                         '')
        call self%keywords_class3D%push('NrClasses',                         '')
        call self%keywords_class3D%push('NrBodies',                          '')
        call self%keywords_class3D%push('NrGroups',                          '')
        call self%keywords_class3D%push('Tau2FudgeFactor',                   '')
        call self%keywords_class3D%push('NormCorrectionAverage',             '')
        call self%keywords_class3D%push('SigmaOffsets',                      '')
        call self%keywords_class3D%push('OrientationalPriorMode',            '')
        call self%keywords_class3D%push('SigmaPriorRotAngle',                '')
        call self%keywords_class3D%push('SigmaPriorTiltAngle',               '')
        call self%keywords_class3D%push('SigmaPriorPsiAngle',                '')
        call self%keywords_class3D%push('LogLikelihood',                     '')
        call self%keywords_class3D%push('AveragePmax',                       '')
        !! data_model_classes head
        call self%keywords_class3D%push('ReferenceImage',                    '')
        call self%keywords_class3D%push('ClassDistribution',                 '')
        call self%keywords_class3D%push('AccuracyRotations',                 '')
        call self%keywords_class3D%push('AccuracyTranslations',              '')
        call self%keywords_class3D%push('EstimatedResolution',               '')
        call self%keywords_class3D%push('OverallFourierCompleteness',        '')
        !! data_model_class
        call self%keywords_class3D%push('SpectralIndex',                     '')
        call self%keywords_class3D%push('Resolution',                        '')
        call self%keywords_class3D%push('AngstromResolution',                'smpd')
        call self%keywords_class3D%push('SsnrMap',                           '')
        call self%keywords_class3D%push('GoldStandardFsc',                   '')
        call self%keywords_class3D%push('FourierCompleteness',               '')
        call self%keywords_class3D%push('ReferenceSigma2',                   '')
        call self%keywords_class3D%push('ReferenceTau2',                     '')
        call self%keywords_class3D%push('SpectralOrientabilityContribution', '')
        !! data_model_groups
        call self%keywords_class3D%push('GroupNumber',                       '')
        call self%keywords_class3D%push('GroupName',                         '')
        call self%keywords_class3D%push('GroupNrParticles',                  '')
        call self%keywords_class3D%push('GroupScaleCorrection',              '')
        !! data_model_pdf_orient_class_1
        call self%keywords_class3D%push('OrientationDistribution',           '')


   self%existence = .true.
    end subroutine new_star_dict

     !>  \brief  is a getter
    pure function exists( self ) result( t )
        class(star_dict), intent(in) :: self
        logical :: t
        t = self%existence
    end function exists
        !>  \brief  is a getter
    function get( self, key ) result( val )
        class(star_dict), intent(inout)    :: self
        character(len=*), intent(in) :: key
        real :: val
        val = self%htab%get(key)
    end function get
    !>  \brief  is a getter
    subroutine getter_1( self, key, val )
        class(star_dict),                    intent(inout) :: self
        character(len=*),              intent(in)    :: key
        character(len=:), allocatable, intent(inout) :: val
        if( allocated(val) ) deallocate(val)
        val = self%chtab%get(key)
    end subroutine getter_1

    !>  \brief  is a getter
    subroutine getter_2( self, key, val )
        class(star_dict),       intent(inout) :: self
        character(len=*), intent(in)    :: key
        real,             intent(inout) :: val
        val = self%htab%get(key)
    end subroutine getter_2
   !>  \brief  returns size of hash
    function hash_size( self ) result( sz )
        class(star_dict), intent(in) :: self
        integer :: sz
        sz = self%htab%size_of()
    end function hash_size

    !>  \brief  returns the keys of the hash
    function hash_vals( self ) result( vals )
        class(star_dict), intent(inout) :: self
        real(kind=4), allocatable :: vals(:)
        vals = self%htab%get_values()
    end function hash_vals
   !>  \brief  check for presence of key in the star_dict hash
    function isthere( self, key ) result( found )
        class(star_dict),       intent(inout) :: self
        character(len=*), intent(in)    :: key
        logical :: hash_found, chash_found, found
        hash_found  = self%htab%isthere(key)
        chash_found = self%chtab%isthere(key)
        found = .false.
        if( hash_found .or. chash_found ) found = .true.
    end function isthere

    !>  \brief  writes orientation info
    subroutine write( self, fhandle )
        class(star_dict), intent(inout) :: self
        integer,    intent(in)    :: fhandle
        character(len=:), allocatable :: str
        str = self%star_dict2str()
        if( allocated(str) ) write(fhandle,'(a)') str
    end subroutine write

   !>  \brief  reads all orientation info (by line) into the hash-tables
    subroutine read( self, fhandle )
        use simple_sauron, only: sauron_line_parser
        class(star_dict), intent(inout) :: self
        integer,    intent(in)    :: fhandle
        character(len=2048) :: line
        logical :: isthere(3)
        read(fhandle,fmt='(A)') line
        call sauron_line_parser( line, self%htab, self%chtab )
        isthere(1) = self%htab%isthere('e1')
        isthere(2) = self%htab%isthere('e2')
        isthere(3) = self%htab%isthere('e3')
        if( any(isthere) )&
          &call self%set_euler([self%htab%get('e1'),self%htab%get('e2'),self%htab%get('e3')])
    end subroutine read


    subroutine print_star_dict( self )
        class(star_dict), intent(inout) :: self
        call self%htab%print()
        call self%chtab%print_key_val_pairs
    end subroutine print_star_dict
        !>  \brief  is a destructor
    subroutine kill( self )
        class(star_dict), intent(inout) :: self
        if( self%existence )then
            call self%htab%kill
            call self%keyword%kill
            self%existence = .false.
        endif
    end subroutine kill


end module simple_star_map
