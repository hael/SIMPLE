submodule(simple_sp_project) simple_sp_project_out
implicit none
#include "simple_local_flags.inc"
contains

    module subroutine add_cavgs2os_out( self, stk, smpd, imgkind, clspath )
        class(sp_project),          intent(inout) :: self
        class(string),              intent(in)    :: stk
        real,                       intent(in)    :: smpd ! sampling distance of images in stk
        character(len=*), optional, intent(in)    :: imgkind
        logical,          optional, intent(in)    :: clspath
        type(string) :: iimgkind, abspath
        integer      :: ldim(3), nptcls, ind
        if( present(imgkind) )then
            iimgkind = imgkind
        else
            iimgkind = 'cavg'
        endif
        ! path and existence check
        abspath = simple_abspath(stk)
        ! find dimension of inputted stack
        call find_ldim_nptcls(abspath, ldim, nptcls)
        ! add os_out entry
        call self%add_entry2os_out(iimgkind%to_char(), ind)
        ! fill-in field
        call self%os_out%set(ind, 'stk',     abspath%to_char())
        call self%os_out%set(ind, 'box',     ldim(1))
        call self%os_out%set(ind, 'nptcls',  nptcls)
        call self%os_out%set(ind, 'fromp',   1)
        call self%os_out%set(ind, 'top',     nptcls)
        call self%os_out%set(ind, 'smpd',    smpd)
        call self%os_out%set(ind, 'stkkind', 'single')
        call self%os_out%set(ind, 'imgkind', iimgkind%to_char())
        call self%os_out%set(ind, 'ctf',     'no')
        if( present(clspath) )then
            if( clspath ) call self%os_out%set(ind, 'stkpath', CWD_GLOB)
        else
            if(self%os_out%isthere(ind, 'stkpath')) call self%os_out%delete_entry(ind, 'stkpath')
        endif
        ! add congruent os_cls2D & os_cls3D
        if( self%os_cls2D%get_noris() /= nptcls )then
            call self%os_cls2D%new(nptcls, is_ptcl=.false.)
            call self%os_cls2D%set_all2single('state',1.)
        endif
        if( self%os_cls3D%get_noris() /= nptcls ) self%os_cls3D = self%os_cls2D
    end subroutine add_cavgs2os_out

    module subroutine add_frcs2os_out( self, frc, which_imgkind )
        class(sp_project), intent(inout) :: self
        class(string),     intent(in)    :: frc
        character(len=*),  intent(in)    :: which_imgkind
        type(string) :: path
        integer      :: ind
        select case(trim(which_imgkind))
            case('frc2D','frc3D')
                ! all good
            case DEFAULT
                THROW_HARD('invalid FRC kind: '//trim(which_imgkind)//'; add_frcs2os_out')
        end select
        ! full path and existence check
        path = simple_abspath(frc)
        ! add os_out entry
        call self%add_entry2os_out(which_imgkind, ind)
        ! fill-in field
        call self%os_out%set(ind, 'frcs',    path)
        call self%os_out%set(ind, 'imgkind', which_imgkind)
    end subroutine add_frcs2os_out

    module subroutine add_fsc2os_out( self, fsc, state, box )
        class(sp_project), intent(inout) :: self
        class(string),     intent(in)    :: fsc
        integer,           intent(in)    :: state, box
        type(string) :: imgkind, abspath
        integer      :: i, ind, n_os_out
        ! full path and existence check
        abspath = simple_abspath(fsc)
        ! add os_out entry
        ! check if field is empty
        n_os_out = self%os_out%get_noris()
        if( n_os_out == 0 )then
            n_os_out = 1
            ind      = 1
            call self%os_out%new(n_os_out, is_ptcl=.false.)
        else
            ind = 0
            do i=1,n_os_out
                if( self%os_out%isthere(i,'imgkind') )then
                    call self%os_out%getter(i,'imgkind',imgkind)
                    if(imgkind.eq.'fsc')then
                        if( self%os_out%get_state(i) == state )then
                            ind = i
                            exit
                        endif
                    endif
                endif
            end do
            if( ind == 0 )then
                n_os_out = n_os_out + 1
                call self%os_out%reallocate(n_os_out)
                ind = n_os_out
            endif
        endif
        ! fill-in field
        call self%os_out%set(ind, 'fsc',     abspath)
        call self%os_out%set(ind, 'imgkind', 'fsc')
        call self%os_out%set(ind, 'state',   state)
        call self%os_out%set(ind, 'box',     box)
    end subroutine add_fsc2os_out

    module subroutine add_vol2os_out( self, vol, smpd, state, which_imgkind, box, pop )
        class(sp_project), intent(inout) :: self
        class(string),     intent(in)    :: vol
        character(len=*),  intent(in)    :: which_imgkind
        real,              intent(in)    :: smpd
        integer,           intent(in)    :: state
        integer, optional, intent(in)    :: box, pop
        type(string) :: imgkind, abspath
        integer                       :: n_os_out, ind, i, ldim(3), ifoo
        select case(trim(which_imgkind))
            case('vol_cavg','vol','vol_msk')
                ! find_dimension of inputted volume
                call find_ldim_nptcls(vol, ldim, ifoo)
                if(present(box))then
                    if( ldim(1) /= box )then
                        THROW_HARD('invalid dimensions for volume: '//vol%to_char()//'; add_vol2os_out 1')
                    endif
                endif
            case DEFAULT
                THROW_HARD('invalid VOL kind: '//trim(which_imgkind)//'; add_vol2os_out 3')
        end select
        ! path and existence check
        abspath = simple_abspath(vol)
        ! check if field is empty
        n_os_out = self%os_out%get_noris()
        if( n_os_out == 0 )then
            n_os_out = 1
            ind      = 1
            call self%os_out%new(n_os_out, is_ptcl=.false.)
        else
            select case(trim(which_imgkind))
                case('vol_msk')
                    ! one volume type for all states
                    ind = 0
                    do i=1,n_os_out
                        if( self%os_out%isthere(i,'imgkind') )then
                            call self%os_out%getter(i,'imgkind',imgkind)
                            if(imgkind.eq.which_imgkind)then
                                ind = i
                                exit
                            endif
                        endif
                    end do
                    if( ind == 0 )then
                        n_os_out = n_os_out + 1
                        call self%os_out%reallocate(n_os_out)
                        ind = n_os_out
                    endif
                case DEFAULT
                    ! one volume per state
                    ind = 0
                    do i=1,n_os_out
                        if( self%os_out%isthere(i,'imgkind') )then
                            call self%os_out%getter(i,'imgkind',imgkind)
                            if(imgkind.eq.which_imgkind)then
                                if( self%os_out%get_state(i) == state )then
                                    ind = i
                                    exit
                                endif
                            endif
                        endif
                    end do
                    if( ind == 0 )then
                        n_os_out = n_os_out + 1
                        call self%os_out%reallocate(n_os_out)
                        ind = n_os_out
                    endif
            end select
        endif
        ! fill-in field
        call self%os_out%set(ind, 'vol',     abspath)
        call self%os_out%set(ind, 'box',     ldim(1))
        call self%os_out%set(ind, 'smpd',    smpd)
        call self%os_out%set(ind, 'imgkind', which_imgkind)
        call self%os_out%set(ind, 'state',   state)
        if(present(pop)) call self%os_out%set(ind, 'pop', pop)
    end subroutine add_vol2os_out

    module subroutine add_entry2os_out( self, which_imgkind, ind )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: which_imgkind
        integer,           intent(out)   :: ind
        type(string) :: imgkind
        integer :: n_os_out, i
        ! check if field is empty
        n_os_out = self%os_out%get_noris()
        if( n_os_out == 0 )then
            n_os_out = 1
            ind      = 1
            call self%os_out%new(n_os_out, is_ptcl=.false.)
        else
            ind = 0
            do i=1,n_os_out
                if( self%os_out%isthere(i,'imgkind') )then
                    imgkind = self%os_out%get_str(i,'imgkind')
                    if(imgkind.eq.trim(which_imgkind))then
                        ind = i
                        exit
                    endif
                endif
            end do
            if( ind == 0 )then
                n_os_out = n_os_out + 1
                call self%os_out%reallocate(n_os_out)
                ind = n_os_out
            endif
        endif
    end subroutine add_entry2os_out

    module subroutine remove_entry_from_osout( self, which_imgkind, state )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: which_imgkind
        integer,           intent(in)    :: state
        type(string) :: imgkind
        integer :: n_os_out, i, ind
        n_os_out = self%os_out%get_noris()
        if( n_os_out == 0 ) return
        ind = 0
        do i = 1,n_os_out
            if( self%os_out%isthere(i,'imgkind') )then
                imgkind = self%os_out%get_str(i,'imgkind')
                if(imgkind.eq.trim(which_imgkind))then
                    if( self%os_out%isthere(i, 'state') )then
                        if( self%os_out%get_state(i) == state )then
                            ind = i
                            exit
                        endif
                    endif
                endif
            endif
        enddo
        if( ind == 0 ) return ! entry not found
        call self%os_out%delete(ind)
    end subroutine remove_entry_from_osout

    ! Getters

    module logical function isthere_in_osout( self, which_imgkind, state )
        class(sp_project), intent(in) :: self
        character(len=*),  intent(in) :: which_imgkind
        integer,           intent(in) :: state
        type(string) :: imgkind
        integer :: n_os_out, i, ind
        isthere_in_osout = .false.
        n_os_out = self%os_out%get_noris()
        if( n_os_out == 0 ) return
        ind = 0
        do i = 1,n_os_out
            if( self%os_out%isthere(i,'imgkind') )then
                imgkind = self%os_out%get_str(i,'imgkind')
                if(imgkind.eq.trim(which_imgkind))then
                    if( self%os_out%isthere(i, 'state') )then
                        if( self%os_out%get_state(i) == state )then
                            ind = i
                            exit
                        endif
                    endif
                endif
            endif
        enddo
        if( ind > 0 ) isthere_in_osout = .true.
    end function isthere_in_osout

    module subroutine get_vol( self, imgkind, state, vol_fname, smpd, box )
        class(sp_project), intent(in)    :: self
        character(len=*),  intent(in)    :: imgkind
        integer,           intent(in)    :: state
        type(string),      intent(inout) :: vol_fname
        real,              intent(out)   :: smpd
        integer,           intent(out)   :: box
        type(string) :: imgkind_here
        integer      :: i, ind, cnt
        select case(trim(imgkind))
            case('vol_cavg','vol','vol_msk')
                ! all good
            case DEFAULT
                THROW_HARD('invalid VOL kind: '//trim(imgkind)//'; get_vol')
        end select
        ! defaults
        call vol_fname%kill
        vol_fname = ''
        box  = 0
        smpd = 0.
        ! fetch index
        ind = 0
        cnt = 0
        select case(trim(imgkind))
            case('vol_msk')
                do i=1,self%os_out%get_noris()
                    if( self%os_out%isthere(i,'imgkind') )then
                        call self%os_out%getter(i,'imgkind',imgkind_here)
                        if(imgkind_here.eq.imgkind)then
                            ind = i
                            cnt = cnt + 1
                        endif
                    endif
                enddo
            case DEFAULT
                do i=1,self%os_out%get_noris()
                    if( self%os_out%isthere(i,'imgkind') )then
                        call self%os_out%getter(i,'imgkind',imgkind_here)
                        if(imgkind_here.eq.imgkind)then
                            if( self%os_out%get_state(i).eq.state )then
                                ind = i
                                cnt = cnt + 1
                            endif
                        endif
                    endif
                enddo
        end select
        if( cnt == 0 )then
            if( trim(imgkind).eq.'vol_msk')then
                ! we do not fall over if the volume mask is absent
                return
            else
                THROW_HARD('no os_out entry with imgkind=volXXX identified, aborting...; get_vol')
            endif
        endif
        if( cnt > 1 )  THROW_HARD('multiple os_out entries with imgkind=volXXX, aborting...; get_vol')
        ! set output
        call vol_fname%kill
        call self%os_out%getter(ind, 'vol', vol_fname)
        smpd = self%os_out%get(ind, 'smpd')
        box  = self%os_out%get_int(ind, 'box')
    end subroutine get_vol

    module subroutine get_all_vols( self, orisout )
        class(sp_project), intent(in)    :: self
        type(oris),        intent(inout) :: orisout
        type(string) :: imgkind_here
        integer      :: i, nvols
        nvols = 0
        call orisout%new(0, .false.)
        do i=1, self%os_out%get_noris()
            if( self%os_out%isthere(i, 'imgkind') ) then
                call self%os_out%getter(i, 'imgkind', imgkind_here)
                if(imgkind_here%to_char() .eq. 'vol') then
                    nvols = nvols + 1
                    if(nvols .eq. 1) then
                        call orisout%new(1, .false.)
                    else
                        call orisout%reallocate(nvols)
                    end if
                    call orisout%transfer_ori(nvols, self%os_out, i)
                endif
            endif
        enddo
        call imgkind_here%kill
    end subroutine get_all_vols

    module subroutine get_fsc( self, state, fsc_fname, box )
        class(sp_project), intent(in)    :: self
        integer,           intent(in)    :: state
        class(string),     intent(inout) :: fsc_fname
        integer,           intent(out)   :: box
        type(string) :: imgkind_here
        integer :: i, ind, cnt
        ! defaults
        call fsc_fname%kill
        fsc_fname =NIL
        ! fetch index
        ind   = 0
        cnt   = 0
        do i=1,self%os_out%get_noris()
            if( self%os_out%isthere(i,'imgkind') )then
                call self%os_out%getter(i,'imgkind',imgkind_here)
                if(imgkind_here%to_char().eq.'fsc')then
                    if( self%os_out%get_state(i).eq.state )then
                        ind = i
                        cnt = cnt + 1
                    endif
                endif
            endif
        enddo
        if( cnt == 0 )THROW_HARD('no os_out entry with imgkind=fsc identified, aborting...; get_fsc')
        if( cnt > 1 ) THROW_HARD('multiple os_out entries with imgkind=fsc, aborting...; get_fsc')
        ! set output
        call fsc_fname%kill
        call self%os_out%getter(ind, 'fsc', fsc_fname)
        box = self%os_out%get_int(ind, 'box')
    end subroutine get_fsc

    module subroutine get_all_fscs( self, orisout )
        class(sp_project), intent(in)    :: self
        type(oris),        intent(inout) :: orisout
        type(string) :: imgkind_here
        integer      :: i, nvols
        nvols = 0
        call orisout%new(0, .false.)
        do i=1, self%os_out%get_noris()
            if( self%os_out%isthere(i, 'imgkind') ) then
                call self%os_out%getter(i, 'imgkind', imgkind_here)
                if(imgkind_here%to_char().eq.'fsc') then
                    nvols = nvols + 1
                    if(nvols .eq. 1) then
                        call orisout%new(1, .false.)
                    else
                        call orisout%reallocate(nvols)
                    end if
                    call orisout%transfer_ori(nvols, self%os_out, i)
                endif
            endif
        enddo
        call imgkind_here%kill
    end subroutine get_all_fscs

    module subroutine get_frcs( self, frcs, which_imgkind, fail )
        class(sp_project), intent(in)    :: self
        class(string),     intent(inout) :: frcs
        character(len=*),  intent(in)    :: which_imgkind
        logical, optional, intent(in)    :: fail
        type(string) :: imgkind
        integer      :: n_os_out, ind, i, cnt
        logical      :: fail_here, found
        select case(trim(which_imgkind))
            case('frc2D','frc3D')
                ! all good
            case DEFAULT
                THROW_HARD('invalid FRC kind: '//trim(which_imgkind)//'; get_frcs')
        end select
        fail_here = .true.
        if( present(fail) )fail_here = fail
        ! check if field is empty
        n_os_out = self%os_out%get_noris()
        if( n_os_out == 0 ) THROW_HARD('trying to fetch from empty os_out field; get_frcs')
        ! look for cavgs
        ind   = 0
        cnt   = 0
        found = .false.
        do i=1,n_os_out
            if( self%os_out%isthere(i,'imgkind') )then
                call self%os_out%getter(i,'imgkind',imgkind)
                if(imgkind.eq.which_imgkind)then
                    ind   = i
                    cnt   = cnt + 1
                    found = .true.
                endif
            endif
        end do
        call frcs%kill
        if( fail_here )then
            if( cnt > 1 )  THROW_HARD('multiple os_out entries with imgkind=frcXD, aborting...; get_frcs')
            if( cnt == 0 ) THROW_HARD('no os_out entry with imgkind=frcsXD identified, aborting...; get_frcs')
        endif
        ! set return values
        if( found )then
            call self%os_out%getter(ind,'frcs',frcs)
        else
            frcs = NIL
        endif
    end subroutine get_frcs

    ! static for OpenMP safety
    module subroutine get_imginfo_from_osout( self, smpd, box, nptcls )
        class(sp_project), intent(inout) :: self
        real,              intent(out)   :: smpd
        integer,           intent(out)   :: box, nptcls
        character(len=STDLEN) :: imgkind
        integer :: i, nos
        nptcls = 0
        smpd   = 0.
        box    = 0
        nos    = self%os_out%get_noris()
        if( nos == 0 )return
        ! nptcls: look for cavgs, defaults to zero for other entries (frcs, volumes)
        do i=1,nos
            if( self%os_out%isthere(i,'imgkind') )then
                call self%os_out%get_static(i,'imgkind',imgkind)
                if(trim(imgkind).eq.'cavg')then
                    if( self%os_out%isthere(i,'fromp').and.self%os_out%isthere(i,'top') )then
                        nptcls = nptcls + self%os_out%get_top(i) - self%os_out%get_fromp(i) + 1
                    else
                        THROW_HARD('Missing fromp and top entries in cavg; get_imginfo_from_osout')
                    endif
                endif
            endif
        end do
        ! box/smpd: first in
        do i=1,nos
            if( self%os_out%isthere(i,'smpd').and.self%os_out%isthere(i,'box') )then
                smpd = self%os_out%get(i,'smpd')
                box  = self%os_out%get_int(i,'box')
                return
            endif
        end do
    end subroutine get_imginfo_from_osout

    ! static for OpenMP safety
    module subroutine get_imgdims_from_osout( self, iseg, smpd, box )
        class(sp_project), intent(inout) :: self
        integer,           intent(in)    :: iseg
        real,              intent(out)   :: smpd
        integer,           intent(out)   :: box
        character(len=STDLEN) :: imgkind_target, imgkind
        integer :: i, nos
        smpd   = 0.
        box    = 0
        nos    = self%os_out%get_noris()
        if( nos == 0 )return
        select case(iseg)
            case(PTCL2D_SEG)
                imgkind_target = 'cavg'
            case(CLS3D_SEG)
                imgkind_target = 'vol_cavg'
            case(PTCL3D_SEG)
                imgkind_target = 'vol'
            case DEFAULT
                return
        end select
        ! last record first
        do i = nos,1,-1
            if( self%os_out%isthere(i,'imgkind') )then
                call self%os_out%get_static(i,'imgkind',imgkind)
                if( trim(imgkind) .eq. trim(imgkind_target) )then
                    if( self%os_out%isthere(i,'smpd').and.self%os_out%isthere(i,'box') )then
                        smpd = self%os_out%get(i,'smpd')
                        box  = self%os_out%get_int(i,'box')
                        return
                    endif
                endif
            endif
        end do
    end subroutine get_imgdims_from_osout

end submodule simple_sp_project_out
