! downscaling of image stacks
module simple_scaler
#include "simple_lib.f08"
use simple_cmdline, only: cmdline
implicit none

public :: scaler
private

type :: scaler
    private
    type(cmdline)         :: cline_scale
    character(len=STDLEN) :: original_stk, stk_sc
    real                  :: original_smpd, original_msk, smpd_sc, scale, msk_sc
    logical               :: stkname_changed
    integer               :: original_box, box_sc, nptcls
  contains
    ! init/uninit
    procedure :: init
    procedure :: uninit
    ! exec
    procedure :: scale_exec
    procedure :: scale_distr_exec
    ! getters
    procedure :: update_smpd_msk
    procedure :: update_stk_smpd_msk
    procedure :: get_scaled_var
    procedure :: get_original_var
end type scaler

contains

    subroutine init( self, p_master, cline, smpd_target, stkscaledbody )
        use simple_params,      only: params
        class(scaler)              :: self
        class(params)              :: p_master
        class(cmdline)             :: cline
        real                       :: smpd_target
        character(len=*), optional :: stkscaledbody
        self%stkname_changed = present(stkscaledbody)
        self%original_stk    = p_master%stk
        self%original_smpd   = p_master%smpd
        self%original_msk    = p_master%msk
        self%original_box    = p_master%box_original
        self%nptcls          = p_master%nptcls
        self%cline_scale     = cline
        call autoscale(self%original_box, p_master%smpd,&
        &smpd_target, self%box_sc, self%smpd_sc, self%scale)
        self%msk_sc = self%scale * p_master%msk
        if( self%stkname_changed )then
            self%stk_sc = trim(stkscaledbody)//p_master%ext
            call self%cline_scale%set('outstk', trim(self%stk_sc))
            call cline%set('stk',  trim(self%stk_sc))
        endif
        call self%cline_scale%set('newbox', real(self%box_sc))
        call cline%set('smpd', self%smpd_sc)
        call cline%set('msk',  self%msk_sc)
    end subroutine init

    subroutine uninit( self, cline )
        class(scaler)  :: self
        class(cmdline) :: cline
        if( self%stkname_changed )  call cline%set('stk',  trim(self%original_stk))
        call cline%set('smpd', self%original_smpd)
        call cline%set('msk',  self%original_msk)
    end subroutine uninit

    subroutine scale_exec( self )
        use simple_commander_imgproc, only: scale_commander
        class(scaler)         :: self
        type(scale_commander) :: xscale
        write(*,'(A)') '>>>'
        write(*,'(A)') '>>> AUTO-SCALING IMAGES'
        write(*,'(A)') '>>>'
        call xscale%execute(self%cline_scale)
    end subroutine scale_exec

    subroutine scale_distr_exec( self )
        use simple_commander_distr_wflows, only: scale_stk_parts_commander
        class(scaler)                   :: self
        type(scale_stk_parts_commander) :: xscale_distr
        write(*,'(A)') '>>>'
        write(*,'(A)') '>>> AUTO-SCALING IMAGES'
        write(*,'(A)') '>>>'
        call xscale_distr%execute(self%cline_scale)
    end subroutine scale_distr_exec

    subroutine update_smpd_msk( self, cline, which )
        class(scaler)    :: self
        class(cmdline)   :: cline
        character(len=*) :: which
        select case(which)
            case('scaled')
                call cline%set('smpd', self%smpd_sc)
                call cline%set('msk',  self%msk_sc)
            case('original')
                call cline%set('smpd', self%original_smpd)
                call cline%set('msk',  self%original_msk)
            case DEFAULT
                 write(*,*) 'flag ', trim(which), ' is unsupported'
                stop 'simple_scaler :: update_smpd_msk'
        end select
    end subroutine update_smpd_msk

    subroutine update_stk_smpd_msk( self, cline, which )
        class(scaler)    :: self
        class(cmdline)   :: cline
        character(len=*) :: which
        select case(which)
            case('scaled')
                if( self%stkname_changed ) call cline%set('stk',  self%stk_sc)
                call cline%set('smpd', self%smpd_sc)
                call cline%set('msk',  self%msk_sc)
            case('original')
                if( self%stkname_changed ) call cline%set('stk',  self%original_stk)
                call cline%set('smpd', self%original_smpd)
                call cline%set('msk',  self%original_msk)
            case DEFAULT
                 write(*,*) 'flag ', trim(which), ' is unsupported'
                stop 'simple_scaler :: update_stk_smpd_msk'
        end select
    end subroutine update_stk_smpd_msk

    real function get_scaled_var( self, which )
        class(scaler)    :: self
        character(len=*) :: which
        select case(which)
            case('smpd')
                get_scaled_var = self%smpd_sc
            case('scale')
                get_scaled_var = self%scale
            case('msk')
                get_scaled_var = self%msk_sc
            case('box')
                get_scaled_var = real(self%box_sc)
            case DEFAULT
                write(*,*) 'flag ', trim(which), ' is unsupported'
                stop 'simple_scaler :: get_scaled_var'
        end select
    end function get_scaled_var

    real function get_original_var( self, which )
        class(scaler)    :: self
        character(len=*) :: which
        select case(which)
            case('smpd')
                get_original_var = self%original_smpd
            case('msk')
                get_original_var = self%original_msk
            case('box')
                get_original_var = real(self%original_box)
            case DEFAULT
                write(*,*) 'flag ', trim(which), ' is unsupported'
                stop 'simple_scaler :: get_original_var'
        end select
    end function get_original_var

end module simple_scaler
