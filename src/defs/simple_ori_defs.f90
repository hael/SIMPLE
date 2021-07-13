module simple_ori_defs
implicit none

enum, bind(c)
    enumerator :: I_ANGAST    = 1
    enumerator :: I_CLASS     = 2
    enumerator :: I_CORR      = 3
    enumerator :: I_CS        = 4
    enumerator :: I_DFX       = 5
    enumerator :: I_DFY       = 6
    enumerator :: I_DIST      = 7
    enumerator :: I_DIST_INPL = 8
    enumerator :: I_E1        = 9
    enumerator :: I_E2        = 10
    enumerator :: I_E3        = 11
    enumerator :: I_EO        = 12
    enumerator :: I_FRAC      = 13
    enumerator :: I_FRACA     = 14
    enumerator :: I_INDINSTK  = 15
    enumerator :: I_INPL      = 16
    enumerator :: I_KV        = 17
    enumerator :: I_LP        = 18
    enumerator :: I_MI_CLASS  = 19
    enumerator :: I_MI_PROJ   = 20
    enumerator :: I_MI_STATE  = 21
    enumerator :: I_NPEAKS    = 22
    enumerator :: I_OW        = 23
    enumerator :: I_PHSHIFT   = 24
    enumerator :: I_PROJ      = 25
    enumerator :: I_SHWMEAN   = 26
    enumerator :: I_SHWSTDEV  = 27
    enumerator :: I_SPECSCORE = 28
    enumerator :: I_SPREAD    = 29
    enumerator :: I_STATE     = 30
    enumerator :: I_STKIND    = 31
    enumerator :: I_UPDATECNT = 32
    enumerator :: I_W         = 33
    enumerator :: I_X         = 34
    enumerator :: I_XINCR     = 35
    enumerator :: I_XPOS      = 36
    enumerator :: I_Y         = 37
    enumerator :: I_YINCR     = 38
    enumerator :: I_YPOS      = 39
end enum

integer, parameter :: N_PTCL_ORIPARAMS = 39
integer, parameter :: N_MIC_STK_PARAMS = 8
integer, parameter :: MIC_STK_PARAMS2INCLUDE(8) = [I_ANGAST,I_CS,I_DFX,I_DFY,I_FRACA,I_KV,I_PHSHIFT,I_STATE]
integer, parameter :: CLSDOC_PARAMS2INCLUDE(4)  = [I_CLASS,I_CORR,I_W,I_SPECSCORE]

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
            case('cs')
                get_oriparam_ind = I_CS
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
            case('fraca')
                get_oriparam_ind = I_FRACA
            case('indinstk')
                get_oriparam_ind = I_INDINSTK
            case('inpl')
                get_oriparam_ind = I_INPL
            case('kv')
                get_oriparam_ind = I_KV
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
            case(I_CS)
                flag ='cs'
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
            case(I_FRACA)
                flag ='fraca'
            case(I_INDINSTK)
                flag ='indinstk'
            case(I_INPL)
                flag ='inpl'
            case(I_KV)
                flag ='kv'
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

end module simple_ori_defs
