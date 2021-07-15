module simple_defs_ori
implicit none

enum, bind(c)
    enumerator :: I_ANGAST    = 1
    enumerator :: I_CLASS     = 2
    enumerator :: I_CORR      = 3
    enumerator :: I_DFX       = 4
    enumerator :: I_DFY       = 5
    enumerator :: I_DIST      = 6
    enumerator :: I_DIST_INPL = 7
    enumerator :: I_E1        = 8
    enumerator :: I_E2        = 9
    enumerator :: I_E3        = 10
    enumerator :: I_EO        = 11
    enumerator :: I_FRAC      = 12
    enumerator :: I_INDSTK    = 13
    enumerator :: I_INPL      = 14
    enumerator :: I_LP        = 15
    enumerator :: I_MI_CLASS  = 16
    enumerator :: I_MI_PROJ   = 17
    enumerator :: I_MI_STATE  = 18
    enumerator :: I_NPEAKS    = 19
    enumerator :: I_OW        = 20
    enumerator :: I_PHSHIFT   = 21
    enumerator :: I_PROJ      = 22
    enumerator :: I_SHWMEAN   = 23
    enumerator :: I_SHWSTDEV  = 24
    enumerator :: I_SPECSCORE = 25
    enumerator :: I_SPREAD    = 26
    enumerator :: I_STATE     = 27
    enumerator :: I_STKIND    = 28
    enumerator :: I_UPDATECNT = 29
    enumerator :: I_W         = 30
    enumerator :: I_X         = 31
    enumerator :: I_XINCR     = 32
    enumerator :: I_XPOS      = 33
    enumerator :: I_Y         = 34
    enumerator :: I_YINCR     = 35
    enumerator :: I_YPOS      = 36
end enum

integer, parameter :: N_PTCL_ORIPARAMS = 36

contains

    pure integer function get_oriparam_ind( flag )
        character(len=*), intent(in) :: flag
        get_oriparam_ind = 0
        select case(trim(adjustl(flag)))
            case('angast')
                get_oriparam_ind = I_ANGAST
            case('class')
                get_oriparam_ind = I_CLASS
            case('corr')
                get_oriparam_ind = I_CORR
            case('dfx')
                get_oriparam_ind = I_DFX
            case('dfy')
                get_oriparam_ind = I_DFY
            case('dist')
                get_oriparam_ind = I_DIST
            case('dist_inpl')
                get_oriparam_ind = I_DIST_INPL
            case('e1')
                get_oriparam_ind = I_E1
            case('e2')
                get_oriparam_ind = I_E2
            case('e3')
                get_oriparam_ind = I_E3
            case('eo')
                get_oriparam_ind = I_EO
            case('frac')
                get_oriparam_ind = I_FRAC
            case('indstk')
                get_oriparam_ind = I_INDSTK
            case('inpl')
                get_oriparam_ind = I_INPL
            case('lp')
                get_oriparam_ind = I_LP
            case('mi_class')
                get_oriparam_ind = I_MI_CLASS
            case('mi_proj')
                get_oriparam_ind = I_MI_PROJ
            case('mi_state')
                get_oriparam_ind = I_MI_STATE
            case('npeaks')
                get_oriparam_ind = I_NPEAKS
            case('ow')
                get_oriparam_ind = I_OW
            case('phshift')
                get_oriparam_ind = I_PHSHIFT
            case('proj')
                get_oriparam_ind = I_PROJ
            case('shwmean')
                get_oriparam_ind = I_SHWMEAN
            case('shwstdev')
                get_oriparam_ind = I_SHWSTDEV
            case('specscore')
                get_oriparam_ind = I_SPECSCORE
            case('spread')
                get_oriparam_ind = I_SPREAD
            case('state')
                get_oriparam_ind = I_STATE
            case('stkind')
                get_oriparam_ind = I_STKIND
            case('updatecnt')
                get_oriparam_ind = I_UPDATECNT
            case('w')
                get_oriparam_ind = I_W
            case('x')
                get_oriparam_ind = I_X
            case('xincr')
                get_oriparam_ind = I_XINCR
            case('xpos')
                get_oriparam_ind = I_XPOS
            case('y')
                get_oriparam_ind = I_Y
            case('yincr')
                get_oriparam_ind = I_YINCR
            case('ypos')
                get_oriparam_ind = I_YPOS
        end select
    end function get_oriparam_ind

    pure function get_oriparam_flag( ind ) result( flag )
        integer,   intent(in) :: ind
        character(len=32) :: flag
        select case(ind)
            case(I_ANGAST)
                flag ='angast'
            case(I_CLASS)
                flag ='class'
            case(I_CORR)
                flag ='corr'
            case(I_DFX)
                flag ='dfx'
            case(I_DFY)
                flag ='dfy'
            case(I_DIST)
                flag ='dist'
            case(I_DIST_INPL)
                flag ='dist_inpl'
            case(I_E1)
                flag ='e1'
            case(I_E2)
                flag ='e2'
            case(I_E3)
                flag ='e3'
            case(I_EO)
                flag ='eo'
            case(I_FRAC)
                flag ='frac'
            case(I_INDSTK)
                flag ='indstk'
            case(I_INPL)
                flag ='inpl'
            case(I_LP)
                flag ='lp'
            case(I_MI_CLASS)
                flag ='mi_class'
            case(I_MI_PROJ)
                flag ='mi_proj'
            case(I_MI_STATE)
                flag ='mi_state'
            case(I_NPEAKS)
                flag ='npeaks'
            case(I_OW)
                flag ='ow'
            case(I_PHSHIFT)
                flag ='phshift'
            case(I_PROJ)
                flag ='proj'
            case(I_SHWMEAN)
                flag ='shwmean'
            case(I_SHWSTDEV)
                flag ='shwstdev'
            case(I_SPECSCORE)
                flag ='specscore'
            case(I_SPREAD)
                flag ='spread'
            case(I_STATE)
                flag ='state'
            case(I_STKIND)
                flag ='stkind'
            case(I_UPDATECNT)
                flag ='updatecnt'
            case(I_W)
                flag ='w'
            case(I_X)
                flag ='x'
            case(I_XINCR)
                flag ='xincr'
            case(I_XPOS)
                flag ='xpos'
            case(I_Y)
                flag ='y'
            case(I_YINCR)
                flag ='yincr'
            case(I_YPOS)
                flag ='ypos'
        end select
    end function get_oriparam_flag

    pure logical function oriparam_isthere( ind, val )
        integer, intent(in) :: ind
        real,    intent(in) :: val
        real, parameter :: TINY = 1e-10
        logical :: is_zero
        is_zero = .false.
        select case(ind)
            case(I_CLASS)
                is_zero = abs(val) < TINY
            case(I_CORR)
                is_zero = abs(val) < TINY
            case(I_DFX)
                is_zero = abs(val) < TINY
            case(I_DFY)
                is_zero = abs(val) < TINY
            case(I_EO)
                is_zero = abs(val) < TINY
            case(I_FRAC)
                is_zero = abs(val) < TINY
            case(I_INDSTK)
                is_zero = abs(val) < TINY
            case(I_INPL)
                is_zero = abs(val) < TINY
            case(I_LP)
                is_zero = abs(val) < TINY
            case(I_NPEAKS)
                is_zero = abs(val) < TINY
            case(I_OW)
                is_zero = abs(val) < TINY
            case(I_PHSHIFT)
                is_zero = abs(val) < TINY
            case(I_PROJ)
                is_zero = abs(val) < TINY
            case(I_SHWMEAN)
                is_zero = abs(val) < TINY
            case(I_SHWSTDEV)
                is_zero = abs(val) < TINY
            case(I_SPECSCORE)
                is_zero = abs(val) < TINY
            case(I_SPREAD)
                is_zero = abs(val) < TINY
            case(I_STKIND)
                is_zero = abs(val) < TINY
            case(I_UPDATECNT)
                is_zero = abs(val) < TINY
            case(I_W)
                is_zero = abs(val) < TINY
            case(I_X)
                is_zero = abs(val) < TINY
            case(I_XINCR)
                is_zero = abs(val) < TINY
            case(I_XPOS)
                is_zero = abs(val) < TINY
            case(I_Y)
                is_zero = abs(val) < TINY
            case(I_YINCR)
                is_zero = abs(val) < TINY
            case(I_YPOS)
                is_zero = abs(val) < TINY
        end select
        oriparam_isthere = .not. is_zero
    end function oriparam_isthere

end module simple_defs_ori
