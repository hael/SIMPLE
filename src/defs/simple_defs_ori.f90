!@descr: static orientation record value enumerators and flag conversions
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
    enumerator :: I_RETIRED_W   = 26
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
    ! slot 42 is spare (see N_PTCL_ORIPARAMS)
    enumerator :: I_FRAC_GREEDY = 43
    enumerator :: I_BETTER_L    = 44
    enumerator :: I_SAMPLED     = 45
    enumerator :: I_CLUSTER     = 46
    enumerator :: I_CLASS_MATCH = 47
    enumerator :: I_CONT_INPL_ATTEMPTED = 48
    enumerator :: I_CONT_INPL_IMPROVED  = 49
    enumerator :: I_RES05               = 50 ! resolution @ FSC=0.5 (formerly the spare slot I_EMPTY10)
    enumerator :: I_CORR_CART           = 51 ! score of the last Cartesian pass (cc or exp(-L); C14, O8)
    enumerator :: I_POSE_CONT_IMPROVED  = 52 ! 1 when the last Cartesian pass improved the pose (O8)
    enumerator :: I_CFAR                = 53 ! conical FSC area ratio of the state's latest half-map pair
end enum

! A particle carries the named slots 1-53 in memory (ori%pparms; slot 42 is spare).
! On disk its record is N_PTCL_RECORD_REALS wide: slots 54-64 are zero padding reserved
! for later fields, so adding a field does not change the project file's record width,
! while resident particles pay only for the named slots. Projects written by earlier
! releases (narrower records) still read, with zeros in the missing slots
! (simple_binoris). A spare slot has no flag and is never "there".
integer, parameter :: N_PTCL_ORIPARAMS      = 53
integer, parameter :: N_PTCL_RECORD_REALS   = 64
integer, parameter :: I_LAST_NAMED_ORIPARAM = I_CFAR
integer, parameter :: I_SPARE_ORIPARAM42    = 42

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
            case('res05')
                get_oriparam_ind = I_RES05
            case('cfar')
                get_oriparam_ind = I_CFAR
            case('state')
                get_oriparam_ind = I_STATE
            case('stkind')
                get_oriparam_ind = I_STKIND
            case('updatecnt')
                get_oriparam_ind = I_UPDATECNT
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
            case('frac_greedy')
                get_oriparam_ind = I_FRAC_GREEDY
            case('better_l')
                get_oriparam_ind = I_BETTER_L
            case('sampled')
                get_oriparam_ind = I_SAMPLED
            case('cluster')
                get_oriparam_ind = I_CLUSTER
            case('class_match')
                get_oriparam_ind = I_CLASS_MATCH
            case('cont_inpl_attempted')
                get_oriparam_ind = I_CONT_INPL_ATTEMPTED
            case('cont_inpl_improved')
                get_oriparam_ind = I_CONT_INPL_IMPROVED
            case('corr_cart')
                get_oriparam_ind = I_CORR_CART
            case('pose_cont_improved')
                get_oriparam_ind = I_POSE_CONT_IMPROVED
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
            case(I_FRAC_GREEDY)
                flag = 'frac_greedy'
            case(I_BETTER_L)
                flag ='better_l'
            case(I_SAMPLED)
                flag ='sampled'
            case(I_CLUSTER)
                flag ='cluster'
            case(I_CLASS_MATCH)
                flag = 'class_match'
            case(I_CONT_INPL_ATTEMPTED)
                flag = 'cont_inpl_attempted'
            case(I_CONT_INPL_IMPROVED)
                flag = 'cont_inpl_improved'
            case(I_CORR_CART)
                flag = 'corr_cart'
            case(I_POSE_CONT_IMPROVED)
                flag = 'pose_cont_improved'
            case(I_RES05)
                flag = 'res05'
            case(I_CFAR)
                flag = 'cfar'
            case DEFAULT
                flag = 'unknown'
        end select
    end function get_oriparam_flag

    pure logical function oriparam_isthere( ind, val )
        integer, intent(in) :: ind
        real,    intent(in) :: val
        real, parameter :: TINY = 1e-10
        oriparam_isthere = .false.
        if( ind < 1 .or. ind > N_PTCL_ORIPARAMS ) return
        if( oriparam_is_spare(ind) ) return
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
            case(I_CLUSTER)
                oriparam_isthere = abs(val) > TINY
            case(I_CLASS_MATCH)
                oriparam_isthere = abs(val) > TINY    
            case(I_RETIRED_W)
                oriparam_isthere = .false.
            case(I_RES05)
                oriparam_isthere = abs(val) > TINY
            case(I_CFAR)
                oriparam_isthere = abs(val) > TINY
            case(I_CORR_CART)
                oriparam_isthere = abs(val) > TINY
            case(I_POSE_CONT_IMPROVED)
                oriparam_isthere = abs(val) > TINY
            case DEFAULT
                ! default case is defined
                oriparam_isthere = .true.
        end select
    end function oriparam_isthere

    !> a spare slot of the particle record (no field yet)
    pure logical function oriparam_is_spare( ind )
        integer, intent(in) :: ind
        oriparam_is_spare = ind == I_SPARE_ORIPARAM42 .or. ind > I_LAST_NAMED_ORIPARAM
    end function oriparam_is_spare

end module simple_defs_ori
