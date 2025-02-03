module simple_defs_ori
implicit none

enum, bind(c)
    enumerator :: I_ANGAST      = 1
    enumerator :: I_CLASS       = 2
    enumerator :: I_CORR        = 3
    enumerator :: I_DFX         = 4
    enumerator :: I_DFY         = 5
    enumerator :: I_DIST        = 6
    enumerator :: I_DIST_INPL   = 7
    enumerator :: I_E1          = 8
    enumerator :: I_E2          = 9
    enumerator :: I_E3          = 10
    enumerator :: I_EO          = 11
    enumerator :: I_FRAC        = 12
    enumerator :: I_INDSTK      = 13
    enumerator :: I_INPL        = 14
    enumerator :: I_LP          = 15
    enumerator :: I_MI_CLASS    = 16
    enumerator :: I_MI_PROJ     = 17
    enumerator :: I_MI_STATE    = 18
    enumerator :: I_PHSHIFT     = 19
    enumerator :: I_PROJ        = 20
    enumerator :: I_SHINCARG    = 21
    enumerator :: I_RES         = 22
    enumerator :: I_STATE       = 23
    enumerator :: I_STKIND      = 24
    enumerator :: I_UPDATECNT   = 25
    enumerator :: I_W           = 26
    enumerator :: I_X           = 27
    enumerator :: I_XINCR       = 28
    enumerator :: I_XPOS        = 29
    enumerator :: I_Y           = 30
    enumerator :: I_YINCR       = 31
    enumerator :: I_YPOS        = 32
    enumerator :: I_GID         = 33
    enumerator :: I_OGID        = 34
    enumerator :: I_PIND        = 35
    enumerator :: I_NEVALS      = 36
    enumerator :: I_NGEVALS     = 37
    enumerator :: I_BETTER      = 38
    enumerator :: I_NPEAKS      = 39
    enumerator :: I_LP_EST      = 40
    enumerator :: I_PIND_PREV   = 41
    enumerator :: I_CC_NONPEAK  = 42 ! unused
    enumerator :: I_FRAC_GREEDY = 43
    enumerator :: I_BETTER_L    = 44
    enumerator :: I_SAMPLED     = 45
    enumerator :: I_CLUSTER     = 46
    ! empties
    enumerator :: I_EMPTY7      = 47
    enumerator :: I_EMPTY8      = 48
    enumerator :: I_EMPTY9      = 49
    enumerator :: I_EMPTY10     = 50
end enum

integer, parameter :: N_PTCL_ORIPARAMS = 50
integer, parameter :: N_NON_EMPTY      = 46

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
            case('phshift')
                get_oriparam_ind = I_PHSHIFT
            case('proj')
                get_oriparam_ind = I_PROJ
            case('shincarg')
                get_oriparam_ind = I_SHINCARG
            case('res')
                get_oriparam_ind = I_RES
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
            case('gid')
                get_oriparam_ind = I_GID
            case('ogid')
                get_oriparam_ind = I_OGID
            case('pind')
                get_oriparam_ind = I_PIND
            case('nevals')
                get_oriparam_ind = I_NEVALS
            case('ngevals')
                get_oriparam_ind = I_NGEVALS
            case('better')
                get_oriparam_ind = I_BETTER
            case('npeaks')
                get_oriparam_ind = I_NPEAKS
            case('lp_est')
                get_oriparam_ind = I_LP_EST
            case('pind_prev')
                get_oriparam_ind = I_PIND_PREV
            case('cc_nonpeak')
                get_oriparam_ind = I_CC_NONPEAK ! unused
            case('frac_greedy')
                get_oriparam_ind = I_FRAC_GREEDY
            case('better_l')
                get_oriparam_ind = I_BETTER_L
            case('sampled')
                get_oriparam_ind = I_SAMPLED
            case('cluster')
                get_oriparam_ind = I_CLUSTER
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
            case(I_PHSHIFT)
                flag ='phshift'
            case(I_PROJ)
                flag ='proj'
            case(I_SHINCARG)
                flag ='shincarg'
            case(I_RES)
                flag ='res'
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
            case(I_GID)
                flag ='gid'
            case(I_OGID)
                flag ='ogid'
            case(I_PIND)
                flag ='pind'
            case(I_NEVALS)
                flag ='nevals'
            case(I_NGEVALS)
                flag ='ngevals'
            case(I_BETTER)
                flag ='better'
            case(I_NPEAKS)
                flag = 'npeaks'
            case(I_LP_EST)
                flag = 'lp_est'
            case(I_PIND_PREV)
                flag = 'pind_prev'
            case(I_CC_NONPEAK)
                flag = 'cc_nonpeak' ! unused
            case(I_FRAC_GREEDY)
                flag = 'frac_greedy'
            case(I_BETTER_L)
                flag ='better_l'
            case(I_SAMPLED)
                flag ='sampled'
            case(I_CLUSTER)
                flag ='cluster'
            case DEFAULT
                flag = 'empty'
        end select
    end function get_oriparam_flag

    pure logical function oriparam_isthere( ind, val )
        integer, intent(in) :: ind
        real,    intent(in) :: val
        real, parameter :: TINY = 1e-10
        oriparam_isthere = .false.
        if( ind < 1 .or. ind > N_NON_EMPTY ) return
        select case(ind)
            ! these variables cannot be zero if defined
            case(I_CLASS)
                oriparam_isthere = abs(val) > TINY
            case(I_CORR)
                oriparam_isthere = abs(val) > TINY
            case(I_DFX)
                oriparam_isthere = abs(val) > TINY
            case(I_DFY)
                oriparam_isthere = abs(val) > TINY
            case(I_EO)
                oriparam_isthere = abs(val) > TINY
            case(I_FRAC)
                oriparam_isthere = abs(val) > TINY
            case(I_INDSTK)
                oriparam_isthere = abs(val) > TINY
            case(I_INPL)
                oriparam_isthere = abs(val) > TINY
            case(I_LP)
                oriparam_isthere = abs(val) > TINY
            case(I_PROJ)
                oriparam_isthere = abs(val) > TINY
            case(I_RES)
                oriparam_isthere = abs(val) > TINY
            case(I_STKIND)
                oriparam_isthere = abs(val) > TINY
            case(I_XPOS)
                oriparam_isthere = abs(val) > TINY
            case(I_YPOS)
                oriparam_isthere = abs(val) > TINY
            case(I_GID)
                oriparam_isthere = abs(val) > TINY
            case(I_OGID)
                oriparam_isthere = abs(val) > TINY
            case(I_PIND)
                oriparam_isthere = abs(val) > TINY
            case(I_LP_EST)
                oriparam_isthere = abs(val) > TINY
            case(I_PIND_PREV)
                oriparam_isthere = abs(val) > TINY
            case(I_CC_NONPEAK)
                oriparam_isthere = abs(val) > TINY
            case(I_CLUSTER)
                oriparam_isthere = abs(val) > TINY
            case DEFAULT
                ! default case is defined
                oriparam_isthere = .true.
        end select
    end function oriparam_isthere

end module simple_defs_ori
