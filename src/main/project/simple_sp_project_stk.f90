submodule(simple_sp_project) simple_sp_project_stk
implicit none
#include "simple_local_flags.inc"
contains

    module subroutine add_stk( self, stk, ctfvars )
        class(sp_project), intent(inout) :: self
        class(string),     intent(in)    :: stk
        type(ctfparams),   intent(in)    :: ctfvars ! All CTF parameters associated with stk
        type(ori)                     :: o
        type(string)    :: stk_relpath, cwd
        integer :: ldim(3), nptcls, n_os_stk, n_os_ptcl2D, n_os_ptcl3D
        integer :: i, fromp, top, pind
        ! full path and existence check
        call simple_getcwd(cwd)
        stk_relpath = simple_abspath(stk)
        ! find dimension of inputted stack
        call find_ldim_nptcls(stk_relpath, ldim, nptcls)
        if( ldim(1) /= ldim(2) )then
            write(logfhandle,*) 'xdim: ', ldim(1)
            write(logfhandle,*) 'ydim: ', ldim(2)
            THROW_HARD('nonsquare particle images not supported; add_stk')
        endif
        ! updates_fields
        n_os_stk    = self%os_stk%get_noris() + 1
        n_os_ptcl2D = self%os_ptcl2D%get_noris()
        n_os_ptcl3D = self%os_ptcl3D%get_noris()
        if( n_os_stk == 1 )then
            call self%os_stk%new(1,         is_ptcl=.false.)
            call self%os_ptcl2D%new(nptcls, is_ptcl=.true.)
            call self%os_ptcl3D%new(nptcls, is_ptcl=.true.)
            fromp = 1
            top   = nptcls
        else
            ! stk
            if( .not.self%os_stk%isthere(n_os_stk-1,'top') )then
                THROW_HARD('FROMP/TOP keys should always be informed; add_stk')
            endif
            call self%os_stk%reallocate(n_os_stk)
            ! 2d
            call self%os_ptcl2D%reallocate(n_os_ptcl2D + nptcls)
            ! 3d
            call self%os_ptcl3D%reallocate(n_os_ptcl3D + nptcls)
            fromp = self%os_stk%get_top(n_os_stk-1) + 1
            top   = fromp + nptcls - 1
        endif
        ! updates oris_objects
        call self%os_stk%set(n_os_stk, 'stk',     stk_relpath)
        call self%os_stk%set(n_os_stk, 'box',     ldim(1))
        call self%os_stk%set(n_os_stk, 'nptcls',  nptcls)
        call self%os_stk%set(n_os_stk, 'fromp',   fromp)
        call self%os_stk%set(n_os_stk, 'top',     top)
        call self%os_stk%set(n_os_stk, 'stkkind', 'split')
        call self%os_stk%set(n_os_stk, 'imgkind', 'ptcl')
        call self%os_stk%set(n_os_stk, 'state',   1) ! default on import
        select case(ctfvars%ctfflag)
            case(CTFFLAG_NO,CTFFLAG_YES,CTFFLAG_FLIP)
                call self%os_stk%set_ctfvars(n_os_stk, ctfvars)
            case DEFAULT
                THROW_HARD('unsupported ctfflag: '//int2str(ctfvars%ctfflag)//'; add_stk')
        end select
        call self%os_stk%set_ctfvars(n_os_stk, ctfvars)
        ! update particle oris objects
        pind = fromp
        do i = 1, nptcls
            call o%new(is_ptcl=.true.)
            call o%set_dfx(      ctfvars%dfx)
            call o%set_dfy(      ctfvars%dfy)
            call o%set('angast', ctfvars%angast)
            if( ctfvars%l_phaseplate ) call o%set('phshift', ctfvars%phshift)
            call o%set('stkind', n_os_stk)
            call o%set('state',  1)         ! default on import
            call o%set('pind',   pind)      ! to keep track of particle indices
            call self%os_ptcl2D%set_ori(n_os_ptcl2D+i, o)
            call self%os_ptcl3D%set_ori(n_os_ptcl3D+i, o)
            pind = pind + 1
        enddo
        call o%kill
    end subroutine add_stk

    module subroutine add_single_stk( self, stk, ctfvars, os )
        class(sp_project), intent(inout) :: self
        class(string),     intent(in)    :: stk
        type(ctfparams),   intent(in)    :: ctfvars ! CTF parameters associated with stk (smpd,kv,cs,fraca,phaseplate)
        class(oris),       intent(inout) :: os      ! parameters associated with stk (dfx,dfy,angast,phshift)
        type(string) :: stk_abspath, projname, fbody
        integer      :: n_os_stk, n_os_ptcl2D, n_os_ptcl3D, ldim(3), nptcls, pind
        call self%projinfo%getter(1, 'projname', projname)
        if( stk%has_substr('mrc') )then
            fbody = get_fbody(basename(stk), string('mrc'))
        else if( stk%has_substr('mrcs') )then
            fbody = get_fbody(basename(stk), string('mrcs'))
        else
            THROW_HARD('Unsupported stack format; use *.mrc or *.mrcs for import')
        endif
        if( projname%has_substr(fbody) ) THROW_HARD('stack for import('//stk%to_char()//') not allowed to have same name as project')
        ! check that stk field is empty
        n_os_stk = self%os_stk%get_noris()
        if( n_os_stk > 0 )then
            write(logfhandle,*) 'stack field (self%os_stk) already populated with # entries: ', n_os_stk
            THROW_HARD('add_single_stk')
        endif
        ! check that particle fields are empty
        n_os_ptcl2D = self%os_ptcl2D%get_noris()
        if( n_os_ptcl2D > 0 )then
            write(logfhandle,*) 'ptcl2D field (self%os_ptcl2D) already populated with # entries: ', n_os_ptcl2D
            THROW_HARD('empty particle fields in project file assumed; add_single_stk')
        endif
        n_os_ptcl3D = self%os_ptcl3D%get_noris()
        if( n_os_ptcl3D > 0 )then
            write(logfhandle,*) 'ptcl3D field (self%os_ptcl3D) already populated with # entries: ', n_os_ptcl3D
            THROW_HARD('empty particle fields in project file assumed; add_single_stk')
        endif
        ! set particle indices
        do pind = 1,os%get_noris()
            call os%set(pind, 'pind', pind)
        end do
        ! copy os
        call self%os_ptcl2D%copy(os, is_ptcl=.true.)
        call self%os_ptcl3D%copy(os, is_ptcl=.true.)
        call self%os_ptcl2D%set_all2single('stkind', 1)
        call self%os_ptcl3D%set_all2single('stkind', 1)
        if( .not. self%os_ptcl2D%isthere('state') ) call self%os_ptcl2D%set_all2single('state',  1) ! default on import
        if( .not. self%os_ptcl3D%isthere('state') ) call self%os_ptcl3D%set_all2single('state',  1) ! default on import
        ! full path and existence check
        stk_abspath = simple_abspath(stk)
        ! find dimension of inputted stack
        call find_ldim_nptcls(stk_abspath, ldim, nptcls)
        if( ldim(1) /= ldim(2) )then
            write(logfhandle,*) 'xdim: ', ldim(1)
            write(logfhandle,*) 'ydim: ', ldim(2)
            THROW_HARD('nonsquare particle images not supported; add_single_stk')
        endif
        ! records
        call self%os_stk%new(1, is_ptcl=.false.)
        call self%os_stk%set(1, 'stk',     stk_abspath)
        call self%os_stk%set(1, 'box',     ldim(1))
        call self%os_stk%set(1, 'nptcls',  nptcls)
        call self%os_stk%set(1, 'fromp',   1)
        call self%os_stk%set(1, 'top',     nptcls)
        call self%os_stk%set(1, 'stkkind', 'single')
        call self%os_stk%set(1, 'imgkind', 'ptcl')
        call self%os_stk%set(1, 'smpd',    ctfvars%smpd)
        call self%os_stk%set(1, 'kv',      ctfvars%kv)
        call self%os_stk%set(1, 'cs',      ctfvars%cs)
        call self%os_stk%set(1, 'fraca',   ctfvars%fraca)
        call self%os_stk%set(1, 'state',   1) ! default on import
        if( ctfvars%l_phaseplate )then
            if( .not. os%isthere('phshift') ) THROW_HARD('phaseplate=yes & input oris lack phshift; add_single_stk')
            call self%os_stk%set(1, 'phaseplate', 'yes')
        else
            call self%os_stk%set(1, 'phaseplate', 'no')
        endif
        select case(ctfvars%ctfflag)
            case(CTFFLAG_NO)
                call self%os_stk%set(1, 'ctf', 'no')
            case(CTFFLAG_YES)
                call self%os_stk%set(1, 'ctf', 'yes')
            case(CTFFLAG_FLIP)
                call self%os_stk%set(1, 'ctf', 'flip')
            case DEFAULT
                THROW_HARD('unsupported ctfflag: '//int2str(ctfvars%ctfflag)//'; add_single_stk')
        end select
    end subroutine add_single_stk

    ! adds stktab given per-stk parameters
    module subroutine add_stktab_1( self, stkfnames, os )
        class(sp_project), intent(inout) :: self
        class(string),     intent(inout) :: stkfnames(:)
        class(oris),       intent(inout) :: os ! parameters associated with stktab
        integer,  allocatable :: nptcls_arr(:)
        type(ori)             :: o_ptcl, o_stk
        integer   :: ldim(3), ldim_here(3), n_os_ptcl2D, n_os_ptcl3D, n_os_stk, istate
        integer   :: i, istk, fromp, top, nptcls, n_os, nstks, nptcls_tot, stk_ind, pind
        nstks = size(stkfnames)
        ! check that inputs are of conforming sizes
        n_os = os%get_noris()
        if( n_os /= nstks )then
            write(logfhandle,*) '# input oris      : ', n_os
            write(logfhandle,*) '# stacks in stktab: ', nstks
            THROW_HARD('nonconforming sizes of inputs; add_stktab')
        endif
        ! first pass for sanity check and determining dimensions
        allocate(nptcls_arr(nstks),source=0)
        do istk=1,nstks
            ! full path and existence check
            stkfnames(istk) = simple_abspath(stkfnames(istk))
            call os%get_ori(istk, o_stk)
            ! logical dimension management
            call find_ldim_nptcls(stkfnames(istk), ldim, nptcls)
            ldim(3) = 1
            if( istk == 1 )then
                ldim_here = ldim
            else
                if( .not. all(ldim_here == ldim) )then
                    write(logfhandle,*) 'micrograph stack #  : ', istk
                    write(logfhandle,*) 'stk name            : ', stkfnames(istk)%to_char()
                    write(logfhandle,*) 'ldim in object      : ', ldim_here
                    write(logfhandle,*) 'ldim read from stack: ', ldim
                    THROW_HARD('inconsistent logical dimensions; add_stktab')
                endif
            endif
            if( ldim(1) /= ldim(2) )then
                write(logfhandle,*) 'stk name: ', stkfnames(istk)%to_char()
                write(logfhandle,*) 'xdim:     ', ldim(1)
                write(logfhandle,*) 'ydim:     ', ldim(2)
                THROW_HARD('nonsquare particle images not supported; add_stktab')
            endif
            ! check variable presence
            if( .not. o_stk%isthere('ctf') )    THROW_HARD('ERROR! ctf flag missing in os input; add_stktab')
            if( .not. o_stk%isthere('smpd') )   THROW_HARD('ERROR! smpd missing in os input; add_stktab')
            if( .not. o_stk%isthere('kv') )     THROW_HARD('ERROR! kv missing in os input; add_stktab')
            if( .not. o_stk%isthere('cs') )     THROW_HARD('ERROR! cs missing in os input; add_stktab')
            if( .not. o_stk%isthere('fraca') )  THROW_HARD('ERROR! fraca missing in os input; add_stktab')
            if( .not. o_stk%isthere('dfx') )    THROW_HARD('ERROR! dfx missing in os input; add_stktab')
            if( .not. o_stk%isthere('dfy') )    THROW_HARD('ERROR! dfy missing in os input; add_stktab')
            if( .not. o_stk%isthere('angast') ) THROW_HARD('ERROR! angast missing in os input; add_stktab')
            ! stash number of images
            nptcls_arr(istk) = nptcls
        enddo
        ! oris allocation
        nptcls_tot  = sum(nptcls_arr)
        n_os_stk    = self%os_stk%get_noris()
        n_os_ptcl2D = self%os_ptcl2D%get_noris()
        n_os_ptcl3D = self%os_ptcl3D%get_noris()
        if( n_os_stk == 0 )then
            call self%os_stk%new(nstks,         is_ptcl=.false.)
            call self%os_ptcl2D%new(nptcls_tot, is_ptcl=.true.)
            call self%os_ptcl3D%new(nptcls_tot, is_ptcl=.true.)
            fromp = 1
        else
            if( .not.self%os_stk%isthere(n_os_stk,'top') )then
                THROW_HARD('FROMP/TOP keys should always be informed; add_stk')
            endif
            call self%os_stk%reallocate(n_os_stk + nstks)
            call self%os_ptcl2D%reallocate(n_os_ptcl2D + nptcls_tot)
            call self%os_ptcl3D%reallocate(n_os_ptcl3D + nptcls_tot)
            fromp = self%os_stk%get_top(n_os_stk) + 1
        endif
        ! parameters transfer
        do istk=1,nstks
            call os%get_ori(istk, o_stk)
            top     = fromp + nptcls_arr(istk) - 1 ! global index
            stk_ind = n_os_stk + istk
            ! updates stk segment
            call self%os_stk%set_ori(stk_ind, o_stk)
            call self%os_stk%set(stk_ind, 'stk',     stkfnames(istk))
            call self%os_stk%set(stk_ind, 'box',     ldim(1))
            call self%os_stk%set(stk_ind, 'nptcls',  nptcls_arr(istk))
            call self%os_stk%set(stk_ind, 'fromp',   fromp)
            call self%os_stk%set(stk_ind, 'top',     top)
            call self%os_stk%set(stk_ind, 'stkkind', 'split')
            call self%os_stk%set(stk_ind, 'imgkind', 'ptcl')
            istate = 1 ! default on import
            if( o_stk%isthere('state') ) istate = o_stk%get_state()
            call self%os_stk%set(stk_ind, 'state', istate)
            ! updates particles segment
            call o_ptcl%new(is_ptcl=.true.)
            call o_ptcl%set_dfx(      o_stk%get_dfx())
            call o_ptcl%set_dfy(      o_stk%get_dfy())
            call o_ptcl%set('angast', o_stk%get('angast'))
            if( o_stk%isthere('phshift') ) call o_ptcl%set('phshift', o_stk%get('phshift'))
            call o_ptcl%set_stkind(stk_ind)
            call o_ptcl%set_state(istate)
            pind = fromp
            do i=1,nptcls_arr(istk)
                call o_ptcl%set('pind', pind) ! to keep track of particle indices
                call self%os_ptcl2D%set_ori(fromp+i-1, o_ptcl)
                call self%os_ptcl3D%set_ori(fromp+i-1, o_ptcl)
                pind = pind + 1
            enddo
            ! update
            fromp = top + 1 ! global index
        enddo
        call o_ptcl%kill
        call o_stk%kill
    end subroutine add_stktab_1

    ! adds stktab given per-particle parameters
    module subroutine add_stktab_2( self, stkfnames, ctfvars, os )
        class(sp_project), intent(inout) :: self
        class(string),     intent(inout) :: stkfnames(:)
        type(ctfparams),   intent(in)    :: ctfvars
        class(oris),       intent(inout) :: os ! parameters associated with stktab
        integer,  allocatable :: nptcls_arr(:)
        type(ori)             :: o_ptcl, o_stk
        integer   :: ldim(3), ldim_here(3), n_os_ptcl2D, n_os_ptcl3D, n_os_stk, istate
        integer   :: i, istk, fromp, top, nptcls, n_os, nstks, nptcls_tot, stk_ind, pind
        nstks = size(stkfnames)
        n_os  = os%get_noris()
        ! first pass for sanity check and determining dimensions
        allocate(nptcls_arr(nstks),source=0)
        do istk=1,nstks
            ! full path and existence check
            stkfnames(istk) = simple_abspath(stkfnames(istk))
            ! logical dimension management
            call find_ldim_nptcls(stkfnames(istk), ldim, nptcls)
            ldim(3) = 1
            if( istk == 1 )then
                ldim_here = ldim
            else
                if( .not. all(ldim_here == ldim) )then
                    write(logfhandle,*) 'micrograph stack #  : ', istk
                    write(logfhandle,*) 'stk name            : ', stkfnames(istk)%to_char()
                    write(logfhandle,*) 'ldim in object      : ', ldim_here
                    write(logfhandle,*) 'ldim read from stack: ', ldim
                    THROW_HARD('inconsistent logical dimensions; add_stktab')
                endif
            endif
            if( ldim(1) /= ldim(2) )then
                write(logfhandle,*) 'stk name: ', stkfnames(istk)%to_char()
                write(logfhandle,*) 'xdim:     ', ldim(1)
                write(logfhandle,*) 'ydim:     ', ldim(2)
                THROW_HARD('nonsquare particle images not supported; add_stktab')
            endif
            ! stash number of images
            nptcls_arr(istk) = nptcls
        enddo
        nptcls_tot  = sum(nptcls_arr)
        if( n_os /= nptcls_tot )then
            write(logfhandle,*) '# input oris               : ', n_os
            write(logfhandle,*) '# ptcls in stacks in stktab: ', nptcls_tot
            THROW_HARD('nonconforming sizes of inputs; add_stktab_2')
        endif
        ! oris allocation
        n_os_stk    = self%os_stk%get_noris()
        n_os_ptcl2D = self%os_ptcl2D%get_noris()
        n_os_ptcl3D = self%os_ptcl3D%get_noris()
        if( n_os_stk == 0 )then
            call self%os_stk%new(nstks,         is_ptcl=.false.)
            call self%os_ptcl2D%new(nptcls_tot, is_ptcl=.true.)
            call self%os_ptcl3D%new(nptcls_tot, is_ptcl=.true.)
            fromp = 1
        else
            if( .not.self%os_stk%isthere(n_os_stk,'top') )then
                THROW_HARD('FROMP/TOP keys should always be informed; add_stktab_2')
            endif
            call self%os_stk%reallocate(n_os_stk + nstks)
            call self%os_ptcl2D%reallocate(n_os_ptcl2D + nptcls_tot)
            call self%os_ptcl3D%reallocate(n_os_ptcl3D + nptcls_tot)
            fromp = self%os_stk%get_top(n_os_stk) + 1
        endif
        ! parameters transfer
        do istk=1,nstks
            top     = fromp + nptcls_arr(istk) - 1 ! global index
            stk_ind = n_os_stk + istk
            ! updates stk segment
            call o_stk%new(is_ptcl=.false.)
            call o_stk%set('stk',     stkfnames(istk))
            call o_stk%set('box',     ldim(1))
            call o_stk%set('nptcls',  nptcls_arr(istk))
            call o_stk%set('fromp',   fromp)
            call o_stk%set('top',     top)
            call o_stk%set('stkkind', 'split')
            call o_stk%set('imgkind', 'ptcl')
            call o_stk%set('smpd',    ctfvars%smpd)
            call o_stk%set('kv',      ctfvars%kv)
            call o_stk%set('cs',      ctfvars%cs)
            call o_stk%set('fraca',   ctfvars%fraca)
            call o_stk%set('state',   1.0) ! default on import
            if( ctfvars%l_phaseplate )then
                call o_stk%set('phaseplate', 'yes')
                call o_stk%set('phshift',    ctfvars%phshift)
            else
                call o_stk%set('phaseplate', 'no')
            endif
            select case(ctfvars%ctfflag)
                case(CTFFLAG_NO)
                    call o_stk%set('ctf', 'no')
                case(CTFFLAG_YES)
                    call o_stk%set('ctf', 'yes')
                case(CTFFLAG_FLIP)
                    call o_stk%set('ctf', 'flip')
                case DEFAULT
                    THROW_HARD('unsupported ctfflag: '//int2str(ctfvars%ctfflag)//'; add_stktab_2')
            end select
            call self%os_stk%set_ori(stk_ind, o_stk)
            ! updates particles segment
            do i=1,nptcls_arr(istk)
                pind = fromp+i-1
                call os%get_ori(pind, o_ptcl)
                call o_ptcl%set_stkind(stk_ind)
                call o_ptcl%set('pind',   pind) ! to keep track of particle indices
                istate = 1
                if( o_ptcl%isthere('state') ) istate = o_ptcl%get_state()
                call o_ptcl%set_state(istate)
                select case(ctfvars%ctfflag)
                    case(CTFFLAG_YES,CTFFLAG_FLIP)
                        if( .not.o_ptcl%isthere('dfx') )then
                            call o_ptcl%print_ori
                            THROW_HARD('Missing defocus parameter(s) for particle: '//int2str(pind))
                        endif
                    case DEFAULT
                        ! all good
                end select
                if( ctfvars%l_phaseplate )then
                    if( .not.o_ptcl%isthere('phshift') )then
                        call o_ptcl%print_ori
                        THROW_HARD('Missing phase-shift parameter for particle: '//int2str(pind))
                    endif
                endif
                call self%os_ptcl2D%set_ori(pind, o_ptcl)
                call self%os_ptcl3D%set_ori(pind, o_ptcl)
            enddo
            ! update
            fromp = top + 1 ! global index
        enddo
        call o_ptcl%kill
        call o_stk%kill
    end subroutine add_stktab_2
    
    !>  Only commits to disk when a change to the project is made
    module subroutine split_stk( self, nparts, dir )
        class(sp_project),       intent(inout) :: self
        integer,                 intent(in)    :: nparts
        class(string), optional, intent(in)    :: dir
        character(len=*), parameter   :: EXT = '.mrc'
        type(image)                   :: img
        type(ori)                     :: orig_stk
        type(stack_io)                :: stkio_w
        type(dstack_io)               :: dstkio_r
        type(string) :: tmp_dir, stkpart, stkkind, dest_stkpart, stk_relpath, cwd, stk
        real         :: smpd
        integer      :: parts(nparts,2), ind_in_stk, iptcl, cnt, istk, box, n_os_stk
        integer      :: nptcls, nptcls_part, numlen
        logical      :: l_set_ind_in_stk
        if( nparts < 2 )return
        ! check that stk field is not empty
        n_os_stk = self%os_stk%get_noris()
        if( n_os_stk == 0 )then
            THROW_HARD('No stack to split! split_stk')
        else if( n_os_stk > 1 )then ! re-splitting not supported
            return
        endif
        call self%os_stk%getter(1, 'stkkind', stkkind)
        if( stkkind%to_char() == 'split' ) return
        ! get original simple_parameters
        call self%os_stk%get_ori(1, orig_stk)
        ! copy prep
        nptcls  = self%get_nptcls()
        parts   = split_nobjs_even( nptcls, nparts )
        numlen  = len_trim(int2str(nparts))
        ! images copy
        smpd = orig_stk%get('smpd')
        box  = orig_stk%get_int('box')
        call img%new([box,box,1], smpd)
        call simple_getcwd(cwd)
        if( present(dir) )then
            tmp_dir = filepath(dir,string('tmp_stacks'))
        else
            tmp_dir = filepath(cwd,string('tmp_stacks'))
        endif
        call simple_mkdir(tmp_dir)
        write(logfhandle,'(a)') '>>> SPLITTING STACK INTO PARTS'
        ! just to get the name of the stack to read from
        call self%get_stkname_and_ind('ptcl2D', 1, stk, ind_in_stk)
        call dstkio_r%new(smpd, box)
        do istk = 1,nparts
            call progress(istk,nparts)
            stkpart = filepath(tmp_dir,string('stack_part'//int2str_pad(istk,numlen)//EXT))
            call stkio_w%open(stkpart, smpd, 'write', box=box, is_ft=.false.)
            cnt = 0
            do iptcl = parts(istk,1), parts(istk,2)
                cnt = cnt + 1
                call self%get_stkname_and_ind('ptcl2D', iptcl, stk, ind_in_stk)
                call dstkio_r%read(stk, ind_in_stk, img)
                call stkio_w%write(cnt, img)
            enddo
            call stkpart%kill
            call stkio_w%close
        enddo
        call dstkio_r%kill
        call img%kill
        call self%os_stk%new(nparts, is_ptcl=.false.)
        if( present(dir) )then
           call simple_mkdir(filepath(dir,STKPARTSDIR))
        else
           call simple_mkdir(STKPARTSDIR)
        endif
        l_set_ind_in_stk = self%os_ptcl2D%isthere('indstk')
        do istk = 1,nparts
            ! file stuff
            stkpart = filepath(tmp_dir,string('stack_part'//int2str_pad(istk,numlen)//EXT))
            if( present(dir) )then
                dest_stkpart = filepath(dir,string(STKPARTFBODY//int2str_pad(istk,numlen)//EXT))
            else
                dest_stkpart = STKPARTFBODY//int2str_pad(istk,numlen)//EXT
            endif
            call simple_rename(stkpart, dest_stkpart)
            stk_relpath = simple_abspath(dest_stkpart)
            nptcls_part = parts(istk,2) - parts(istk,1) + 1
            ! set original before overriding
            call self%os_stk%set_ori(istk, orig_stk)
            ! override
            call self%os_stk%set(istk, 'stk',     stk_relpath)
            call self%os_stk%set(istk, 'nptcls',  nptcls_part)
            call self%os_stk%set(istk, 'fromp',   parts(istk,1))
            call self%os_stk%set(istk, 'top',     parts(istk,2))
            call self%os_stk%set(istk, 'stkkind', 'split')
            ind_in_stk = 0
            do iptcl=parts(istk,1),parts(istk,2)
                call self%os_ptcl2D%set(iptcl,'stkind', istk)
                call self%os_ptcl3D%set(iptcl,'stkind', istk)
                if( l_set_ind_in_stk )then
                    ind_in_stk = ind_in_stk + 1
                    call self%os_ptcl2D%set(iptcl,'indstk', ind_in_stk)
                    call self%os_ptcl3D%set(iptcl,'indstk', ind_in_stk)
                endif
            enddo
            call stkpart%kill
            call dest_stkpart%kill
        enddo
        call self%write
        call simple_rmdir(tmp_dir)
        call orig_stk%kill
    end subroutine split_stk

    module subroutine write_substk( self, fromto, stkout )
        class(sp_project), intent(inout) :: self
        integer,           intent(in)    :: fromto(2)
        class(string),     intent(in)    :: stkout
        integer         :: ind_in_stk, iptcl, n_os_stk, nptcls, box, cnt, ffromto(2)
        real            :: smpd
        type(image)     :: img
        type(string)    :: stk
        type(ori)       :: orig_stk
        type(stack_io)  :: stkio_w
        type(dstack_io) :: dstkio_r
        ! check that stk field is not empty
        n_os_stk = self%os_stk%get_noris()
        if( n_os_stk == 0 )then
            THROW_HARD('No stack(s) to extract from! write_substk')
        endif
        ! get original simple_parameters
        call self%os_stk%get_ori(1, orig_stk)
        ! copy prep
        nptcls  = self%get_nptcls()
        ffromto = fromto
        if( ffromto(1) < 1 ) ffromto(1) = 1
        if( ffromto(2) < 1 ) ffromto(2) = nptcls
        ! images copy
        smpd = orig_stk%get('smpd')
        box  = orig_stk%get_int('box')
        call img%new([box,box,1], smpd)
        call dstkio_r%new(smpd, box)
        call stkio_w%open(stkout, smpd, 'write', box=box, is_ft=.false.)
        cnt = 0
        do iptcl = ffromto(1),ffromto(2)
            cnt = cnt + 1
            if( iptcl < 1 .or. iptcl > nptcls ) THROW_HARD('index '//int2str(iptcl)//' out of range')
            call self%get_stkname_and_ind('ptcl2D', iptcl, stk, ind_in_stk)
            call dstkio_r%read(stk, ind_in_stk, img)
            call stkio_w%write(cnt, img)
        end do
        call stkio_w%close
        call dstkio_r%kill
        call img%kill
        call orig_stk%kill
    end subroutine write_substk

    module subroutine add_scale_tag( self, dir )
        class(sp_project),       intent(inout) :: self
        class(string), optional, intent(in)    :: dir
        type(string) :: newname, stkname, abs_dir, nametmp, ext
        integer      :: imic, nmics
        nmics = self%os_stk%get_noris()
        if( present(dir) )then
            call simple_mkdir(dir)
            abs_dir = simple_abspath(dir)
        endif
        do imic=1,nmics
            call self%os_stk%getter(imic, 'stk', stkname)
            ext = fname2ext(stkname)
            if(present(dir))then
                nametmp = basename(add2fbody(stkname, '.'//ext%to_char(), SCALE_SUFFIX))
                newname = filepath(abs_dir, nametmp)
            else
                newname = add2fbody(stkname, '.'//ext%to_char(), SCALE_SUFFIX)
            endif
            newname = fname_new_ext(newname, ext)
            call self%os_stk%set(imic, 'stk', newname)
        end do
    end subroutine add_scale_tag

    ! report state selection to os_stk & os_ptcl2D/3D
    module subroutine report_state2stk( self, states )
        class(sp_project), intent(inout) :: self
        integer,           intent(in)    :: states(:)
        integer :: iptcl, noris_ptcl3D, noris_ptcl2D, istk, fromp, top, nstks, nptcls
        nstks = self%get_nstks()
        if( nstks == 0 )then
            THROW_WARN('empty STK field. Nothing to do; report_state2stk')
            return
        endif
        if(size(states) /= nstks)then
            THROW_WARN('Inconsistent # number of states & stacks; report_state2stk')
            return
        endif
        ! update stacks
        do istk=1,nstks
            call self%os_stk%set(istk, 'state', states(istk))
        enddo
        ! ensure ptcl fields congruent
        noris_ptcl2D = self%os_ptcl2D%get_noris()
        noris_ptcl3D = self%os_ptcl3D%get_noris()
        if( noris_ptcl3D /= noris_ptcl2D )then
            THROW_HARD('Inconsistent # number of 2D/3D particles; report_state2stk')
        else
            do istk=1,nstks
                fromp  = self%os_stk%get_fromp(istk)
                top    = self%os_stk%get_top(istk)
                nptcls = self%os_stk%get_int(istk,'nptcls')
                if(top-fromp+1 /= nptcls)then
                    call self%os_stk%print(istk)
                    THROW_HARD('Incorrect # number of particles in stack '//int2str(istk)//'; report_state2stk')
                endif
                if( states(istk) > 0 )then
                    ! preserve existing states
                    cycle
                else
                    ! de-select
                    do iptcl=fromp,top
                        call self%os_ptcl2D%set(iptcl, 'state', 0.)
                        call self%os_ptcl3D%set(iptcl, 'state', 0.)
                    enddo
                endif
            enddo
        endif
    end subroutine report_state2stk

    ! Getters

    module logical function has_phaseplate( self, oritype )
        class(sp_project), target, intent(inout) :: self
        character(len=*),          intent(in)    :: oritype
        type(string) :: phaseplate
        has_phaseplate = .false.
        if( trim(oritype) .eq. 'cls3D' ) return
        ! get info
        if( self%os_stk%isthere(1, 'phaseplate') )then
            phaseplate = self%os_stk%get_str(1, 'phaseplate')
        else
            phaseplate = 'no'
        endif
        has_phaseplate = phaseplate%to_char().eq.'yes'
    end function has_phaseplate

    module integer function get_box( self )
        class(sp_project), target, intent(in) :: self
        integer :: n_os_stk
        get_box  = 0
        n_os_stk = self%os_stk%get_noris()
        if( n_os_stk == 0 )then
            THROW_HARD('empty os_stk field! get_box')
        endif
        get_box = self%os_stk%get_int(1,'box')
    end function get_box

    module function get_stkname( self, imic ) result( stkname )
        class(sp_project), intent(inout) :: self
        integer,           intent(in)    :: imic
        type(string) :: stkname
        integer :: nmics
        nmics = self%os_stk%get_noris()
        if( imic < 1 .or. imic > nmics )then
            write(logfhandle,*) 'imic : ', imic
            write(logfhandle,*) 'nmics: ', nmics
            THROW_HARD('imic index out of range; get_stkname')
        endif
        stkname = self%os_stk%get_str(imic, 'stk')
    end function get_stkname

    module subroutine get_stkname_and_ind( self, oritype, iptcl, stkname, ind_in_stk )
        class(sp_project), target, intent(inout) :: self
        character(len=*),          intent(in)    :: oritype
        integer,                   intent(in)    :: iptcl
        class(string),             intent(out)   :: stkname
        integer,                   intent(out)   :: ind_in_stk
        real    :: smpd
        integer :: stkind, ncls
        ! do the index mapping
        call self%map_ptcl_ind2stk_ind(oritype, iptcl, stkind, ind_in_stk)
        ! output name
        if( trim(oritype) .eq. 'cls3D' ) then
            call self%get_cavgs_stk(stkname, ncls, smpd)
        else
            stkname = self%os_stk%get_str(stkind, 'stk')
        endif
    end subroutine get_stkname_and_ind

    module integer function get_nstks( self )
        class(sp_project), target, intent(in) :: self
        get_nstks = self%os_stk%get_noris()
    end function get_nstks

end submodule simple_sp_project_stk
