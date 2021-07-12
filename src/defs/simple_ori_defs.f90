module simple_ori_defs
implicit none

integer, parameter :: N_PTCL_ORIPARAMS = 33

contains

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
        enumerator :: I_INPL      = 15
        enumerator :: I_KV        = 16
        enumerator :: I_LP        = 17
        enumerator :: I_MI_CLASS  = 18
        enumerator :: I_MI_PROJ   = 19
        enumerator :: I_NPEAKS    = 20
        enumerator :: I_OW        = 21
        enumerator :: I_PHSHIFT   = 22
        enumerator :: I_PROJ      = 23
        enumerator :: I_SHWMEAN   = 24
        enumerator :: I_SHWSTDEV  = 25
        enumerator :: I_SPECSCORE = 26
        enumerator :: I_SPREAD    = 27
        enumerator :: I_STATE     = 28
        enumerator :: I_STKIND    = 29
        enumerator :: I_UPDATECNT = 30
        enumerator :: I_W         = 31
        enumerator :: I_X         = 32
        enumerator :: I_Y         = 33
    end enum

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
            case('y')
                get_oriparam_ind = I_Y
        end select
    end function get_oriparam_ind

end module simple_ori_defs
