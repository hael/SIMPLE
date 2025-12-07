submodule(simple_sp_project) simple_sp_project_mic
implicit none
#include "simple_local_flags.inc"
contains

    module subroutine add_single_movie( self, moviename, ctfvars )
        class(sp_project), target, intent(inout) :: self
        class(string),             intent(in)    :: moviename
        type(ctfparams),           intent(in)    :: ctfvars
        class(oris), pointer :: os_ptr
        type(string)         :: abs_fname
        integer :: n_os_mic, ldim(3), nframes
        ! oris object pointer
        os_ptr => self%os_mic
        ! check that os_mic field is empty
        n_os_mic = os_ptr%get_noris()
        if( n_os_mic > 0 )then
            write(logfhandle,*) 'mic field (self%os_mic) already populated with # entries: ', n_os_mic
            THROW_HARD('add_single_movie')
        endif
        ! update ori
        call os_ptr%new(1, is_ptcl=.false.)
        abs_fname = simple_abspath(moviename)
        call find_ldim_nptcls(abs_fname, ldim, nframes)
        if( nframes <= 0 )then
            THROW_WARN('# frames in movie: '//abs_fname%to_char()//' <= zero, omitting')
        else if( nframes > 1 )then
            call os_ptr%set(1, 'movie',   abs_fname)
            call os_ptr%set(1, 'imgkind', 'movie')
            call os_ptr%set(1, 'nframes', nframes)
        else
            call os_ptr%set(1, 'intg',    abs_fname)
            call os_ptr%set(1, 'imgkind', 'mic')
        endif
        ! updates segment
        call os_ptr%set(1, 'xdim',       ldim(1))
        call os_ptr%set(1, 'ydim',       ldim(2))
        call os_ptr%set(1, 'smpd',       ctfvars%smpd)
        call os_ptr%set(1, 'kv',         ctfvars%kv)
        call os_ptr%set(1, 'cs',         ctfvars%cs)
        call os_ptr%set(1, 'fraca',      ctfvars%fraca)
        if( ctfvars%l_phaseplate )then
            call os_ptr%set(1, 'phaseplate', 'yes')
        else
            call os_ptr%set(1, 'phaseplate', 'no')
        endif
        select case(ctfvars%ctfflag)
            case(CTFFLAG_NO)
                call os_ptr%set(1, 'ctf', 'no')
            case(CTFFLAG_YES)
                call os_ptr%set(1, 'ctf', 'yes')
            case(CTFFLAG_FLIP)
                call os_ptr%set(1, 'ctf', 'flip')
            case DEFAULT
                THROW_HARD('ctfflag: '//int2str(ctfvars%ctfflag)//' unsupported; add_single_movie')
        end select
        call os_ptr%set(1,'state',1.) ! default on import
    end subroutine add_single_movie

    !> Add/append movies or micrographs without ctf parameters
    module subroutine add_movies( self, movies_array, ctfvars, singleframe, verbose )
        class(sp_project), target, intent(inout) :: self
        type(string),              intent(in)    :: movies_array(:)
        type(ctfparams),           intent(in)    :: ctfvars
        logical,         optional, intent(in)    :: singleframe
        logical,         optional, intent(in)    :: verbose
        class(oris), pointer :: os_ptr
        type(ctfparams)      :: prev_ctfvars
        type(string)         :: name, abs_moviename
        integer :: ldim_orig(3), imic, ldim(3), nframes, nmics, nprev_mics, cnt, ntot, nframes_first
        logical :: is_movie, l_singleframe, l_verbose
        l_verbose     = .true.
        if( present(verbose) ) l_verbose = verbose
        is_movie      = .true.
        l_singleframe = .false.
        if( present(singleframe) ) l_singleframe = singleframe
        ! oris object pointer
        os_ptr => self%os_mic
        ! read movie names
        nmics = size(movies_array)
        ! update oris
        nprev_mics = os_ptr%get_noris()
        ntot       = nmics + nprev_mics
        if( nprev_mics == 0 )then
            call os_ptr%new(ntot, is_ptcl=.false.)
        else
            prev_ctfvars = self%get_micparams(1)
            if( ctfvars%ctfflag /= prev_ctfvars%ctfflag ) THROW_HARD('CTF infos do not match! add_movies')
            if( .not.is_equal(ctfvars%smpd, prev_ctfvars%smpd )) THROW_HARD('The sampling distances do not match! add_movies')
            if( .not.is_equal(ctfvars%cs,   prev_ctfvars%cs   )) THROW_HARD('The spherical aberrations do not match! add_movies')
            if( .not.is_equal(ctfvars%kv,   prev_ctfvars%kv   )) THROW_HARD('The voltages do not match! add_movies')
            if( .not.is_equal(ctfvars%fraca,prev_ctfvars%fraca)) THROW_HARD('The amplitude contrasts do not match! add_movies')
            if( ctfvars%l_phaseplate.neqv.prev_ctfvars%l_phaseplate ) THROW_HARD('Phaseplate infos do not match! add_movies')
            call os_ptr%reallocate(ntot)
        endif
        cnt = 0
        nframes_first = 0
        do imic=nprev_mics + 1,ntot
            cnt = cnt + 1
            abs_moviename = simple_abspath(movies_array(cnt))
            call find_ldim_nptcls(abs_moviename, ldim, nframes)
            if( cnt == 1 )then
                ldim_orig = ldim
            else
                if( ldim(1) /= ldim_orig(1) .or. ldim(2) /= ldim_orig(2) )then
                    write(logfhandle,*)'Inconsistent size for file: ',movies_array(cnt)%to_char()
                    write(logfhandle,*)'Dimensions: ', ldim(1),'x ',ldim(2), ' vs. previous dimensions: ', ldim_orig(1),'x ',ldim_orig(2)
                    THROW_HARD('All files imported must have identical dimensions!')
                endif
            endif
            if( nframes <= 0 )then
                THROW_WARN('# frames in movie: '//movies_array(imic)%to_char()//' <= zero, omitting')
                cycle
            else
                if( nframes > 1 )then
                    call os_ptr%set(imic, 'movie',   abs_moviename)
                    call os_ptr%set(imic, 'imgkind', 'movie')
                    is_movie = .true.
                else
                    if( l_singleframe )then
                        call os_ptr%set(imic, 'frame',   abs_moviename)
                        call os_ptr%set(imic, 'imgkind', 'frame')
                    else
                        call os_ptr%set(imic, 'intg',    abs_moviename)
                        call os_ptr%set(imic, 'imgkind', 'mic')
                    endif
                    is_movie = .false.
                endif
                if( nframes_first == 0 )then
                    nframes_first = nframes
                else
                    if( nframes /= nframes_first )then
                        write(logfhandle,*) abs_moviename%to_char(), ' has ', nframes, ' frame(s)'
                        write(logfhandle,*) 'Previous import have ', nframes_first, ' frame(s)'
                        THROW_HARD('You cannot import both micrographs and movies at the same time! add_movies')
                    endif
                endif
            endif
            ! updates segment
            call os_ptr%set(imic, 'xdim',    ldim(1))
            call os_ptr%set(imic, 'ydim',    ldim(2))
            call os_ptr%set(imic, 'nframes', nframes)
            call os_ptr%set(imic, 'smpd',    ctfvars%smpd)
            call os_ptr%set(imic, 'kv',      ctfvars%kv)
            call os_ptr%set(imic, 'cs',      ctfvars%cs)
            call os_ptr%set(imic, 'fraca',   ctfvars%fraca)
            call os_ptr%set(imic, 'state',   1.0) ! default on import
            if( ctfvars%l_phaseplate )then
                call os_ptr%set(imic, 'phaseplate', 'yes')
            else
                call os_ptr%set(imic, 'phaseplate', 'no')
            endif
            select case(ctfvars%ctfflag)
                case(CTFFLAG_NO)
                    call os_ptr%set(imic, 'ctf', 'no')
                case(CTFFLAG_YES)
                    call os_ptr%set(imic, 'ctf', 'yes')
                case(CTFFLAG_FLIP)
                    call os_ptr%set(imic, 'ctf', 'flip')
                case DEFAULT
                    THROW_HARD('unsupported ctfflag: '//int2str(ctfvars%ctfflag)//'; add_movies')
            end select
        enddo
        if( is_movie )then
            name = 'MOVIE(S)'
        else
            if( l_singleframe )then
                name = 'FRAME(S)'
            else
                name = 'MICROGRAPH(S)'
            endif
        endif
        if( l_verbose )then
            write(logfhandle,'(A13,I6,A1,A)')'>>> IMPORTED ', nmics,' ', name%to_char()
            write(logfhandle,'(A20,A,A1,I6)')'>>> TOTAL NUMBER OF ', name%to_char(),':',ntot
        endif
    end subroutine add_movies

    !> Add/append micrographs with ctf parameters
    module subroutine add_intgs( self, intgs_array, os, ctfvars )
        class(sp_project), target, intent(inout) :: self
        class(string),             intent(in)    :: intgs_array(:)
        class(oris),               intent(in)    :: os
        type(ctfparams),           intent(in)    :: ctfvars
        type(ctfparams)           :: prev_ctfvars, ctfparms
        type(string):: rel_micname
        real                      :: intg_smpd
        integer                   :: imic,ldim(3),nframes,nintgs,nprev_intgs,nprev_mics,cnt,ntot
        nprev_mics  = self%os_mic%get_noris()
        nprev_intgs = self%get_nintgs()
        if( nprev_mics > 0 )then
            if( nprev_mics /= nprev_intgs )then
                THROW_HARD('Cannot add lone micrographs to a project with movies; add_intgs')
            endif
            if( nprev_intgs == 0 )then
                THROW_HARD('Cannot add micrographs to a project with movies only; add_intgs')
            endif
            ! previous micrographs parameters
            prev_ctfvars = self%os_mic%get_ctfvars(1)
            if(.not.is_equal(ctfvars%smpd, prev_ctfvars%smpd )) THROW_HARD('Inconsistent sampling distance; add_intgs')
            if(.not.is_equal(ctfvars%cs,   prev_ctfvars%cs   )) THROW_HARD('Inconsistent spherical aberration; add_intgs')
            if(.not.is_equal(ctfvars%kv,   prev_ctfvars%kv   )) THROW_HARD('Inconsistent voltage; add_intgs')
            if(.not.is_equal(ctfvars%fraca,prev_ctfvars%fraca)) THROW_HARD('Inconsistent amplituce contrast; add_intgs')
            if(ctfvars%ctfflag /= prev_ctfvars%ctfflag) THROW_HARD('Incompatible CTF flag; add_intgs')
            if(ctfvars%l_phaseplate .neqv. prev_ctfvars%l_phaseplate ) THROW_HARD('Incompatible phaseplate info; add_intgs')
        endif
        ! read movie names
        nintgs = size(intgs_array)
        if( nintgs /= os%get_noris() )then
            THROW_HARD('Inconsistent # of mics & ctf parameters; add_intgs')
        endif
        ! update oris
        if( nprev_intgs == 0 )then
            ! first import
            call self%os_mic%new(nintgs, is_ptcl=.false.)
            ntot = nintgs
        else
            ! append
            ntot = nintgs+nprev_intgs
            call self%os_mic%reallocate(ntot)
        endif
        cnt = 0
        do imic=nprev_intgs+1,ntot
            cnt = cnt + 1
            rel_micname = simple_abspath(intgs_array(cnt))
            call find_ldim_nptcls(rel_micname, ldim, nframes, smpd=intg_smpd)
            if( nframes <= 0 )then
                THROW_HARD('# frames in movie: '//intgs_array(cnt)%to_char()//' <= zero; add_intgs')
            else if( nframes > 1 )then
                THROW_HARD('Not the interface for adding movies; add_intgs')
            endif
            if( nprev_intgs > 0 )then
                if( .not.is_equal(intg_smpd,prev_ctfvars%smpd) )then
                    THROW_HARD('Incompatible sampling distance: '//intgs_array(cnt)%to_char()//'; add_intgs')
                endif
            endif
            ! updates segment
            ctfparms = os%get_ctfvars(cnt)
            call self%os_mic%set(imic, 'intg',    rel_micname)
            call self%os_mic%set(imic, 'imgkind', 'mic')
            call self%os_mic%set(imic, 'xdim',    ldim(1))
            call self%os_mic%set(imic, 'ydim',    ldim(2))
            call self%os_mic%set(imic, 'smpd',    ctfvars%smpd)
            call self%os_mic%set(imic, 'kv',      ctfvars%kv)
            call self%os_mic%set(imic, 'cs',      ctfvars%cs)
            call self%os_mic%set(imic, 'fraca',   ctfvars%fraca)
            if( os%isthere(cnt,'state') )then
                call self%os_mic%set(imic, 'state', os%get(cnt,'state'))
            else
                call self%os_mic%set(imic, 'state', 1)
            endif
            if( ctfvars%l_phaseplate )then
                call self%os_mic%set(imic, 'phaseplate', 'yes')
            else
                call self%os_mic%set(imic, 'phaseplate', 'no')
            endif
            select case(ctfvars%ctfflag)
                case(CTFFLAG_NO)
                    call self%os_mic%set(imic, 'ctf', 'no')
                case(CTFFLAG_YES)
                    call self%os_mic%set(imic, 'ctf',    'yes')
                    call self%os_mic%set_dfx(imic,       ctfparms%dfx)
                    call self%os_mic%set_dfy(imic,       ctfparms%dfy)
                    call self%os_mic%set(imic, 'angast', ctfparms%angast)
                    call self%os_mic%set(imic, 'phshift',ctfparms%phshift)
                case(CTFFLAG_FLIP)
                    call self%os_mic%set(imic, 'ctf', 'flip')
            end select
        enddo
        write(logfhandle,'(A,I6,A)')'>>> IMPORTED ', nintgs,' INTEGRATED MOVIES'
        write(logfhandle,'(A,I6)')'>>> TOTAL NUMBER OF MICROGRAPHS:',ntot
    end subroutine add_intgs

    ! report state selection to os_stk & os_ptcl2D/3D
    module subroutine report_state2mic( self, states )
        class(sp_project), intent(inout) :: self
        integer,           intent(in)    :: states(:)
        integer    :: imic, nmics, cnt, nsel
        type(oris) :: tmp
        type(ori)  :: o
        nmics = self%os_mic%get_noris()
        if( nmics == 0 )then
            THROW_WARN('empty MIC field. Nothing to do; report_state2mic')
            return
        endif
        if(size(states) /= nmics)then
            THROW_WARN('Inconsistent # number of states & mics; report_state2mic')
            return
        endif
        nsel = count(states == 1)
        call tmp%new(nsel, is_ptcl=.false.)
        cnt = 0
        do imic=1,nmics
            if( states(imic) == 1 )then
                cnt = cnt + 1
                call self%os_mic%get_ori(imic, o)
                call tmp%set_ori(cnt, o)
            endif
        enddo
        call self%os_mic%copy(tmp, is_ptcl=.false.)
        call tmp%kill
    end subroutine report_state2mic

    ! Getters

    module logical function has_boxfile( self )
        class(sp_project), target, intent(in) :: self
        has_boxfile = self%os_mic%isthere('boxfile')
    end function has_boxfile

    ! returns list of movies regardless of 'imgkind' key as it overrides movies
    module subroutine get_movies_table( self, moviestab )
        class(sp_project),         intent(inout) :: self
        type(string), allocatable, intent(out)   :: moviestab(:)
        integer :: i,n,cnt
        if(allocated(moviestab))deallocate(moviestab)
        n = 0
        do i=1,self%os_mic%get_noris()
            if(self%os_mic%isthere(i,'movie'))then
                if( self%os_mic%get_state(i) > 0 ) n = n+1
            endif
        enddo
        if( n==0 )return
        allocate(moviestab(n))
        cnt = 0
        do i=1,self%os_mic%get_noris()
            if(self%os_mic%isthere(i,'movie'))then
                if( self%os_mic%get_state(i) > 0 )then
                    cnt = cnt + 1
                    call self%os_mic%getter(i,'movie',moviestab(cnt))
                endif
            endif
        enddo
    end subroutine get_movies_table

    module subroutine get_mics_table( self, micstab, orimap)
        class(sp_project),         intent(inout) :: self
        type(string), allocatable, intent(out)   :: micstab(:)
        integer,      allocatable, intent(out)   :: orimap(:)
        integer :: i,n,cnt
        if(allocated(micstab)) deallocate(micstab)
        if(allocated(orimap))  deallocate(orimap)
        n = 0
        do i=1,self%os_mic%get_noris()
            if(self%os_mic%isthere(i,'intg'))then
                if( self%os_mic%get_state(i) > 0 ) n = n + 1
            endif
        enddo
        if( n==0 )return
        allocate(micstab(n), orimap(n))
        cnt = 0
        do i=1,self%os_mic%get_noris()
            if(self%os_mic%isthere(i,'intg'))then
                if( self%os_mic%get_state(i) > 0 )then
                    cnt = cnt + 1
                    call self%os_mic%getter(i,'intg',micstab(cnt))
                    orimap(cnt)  = i
                endif
            endif
        enddo
    end subroutine get_mics_table

    module function get_micparams( self, imic ) result( ctfvars )
        class(sp_project), intent(inout) :: self
        integer,           intent(in)    :: imic
        type(ctfparams) :: ctfvars
        ctfvars = self%os_mic%get_ctfvars(imic)
    end function get_micparams

    ! static for OpenMP safety
    module function get_micname( self, iptcl ) result( micname )
        class(sp_project), intent(inout) :: self
        integer,           intent(in)    :: iptcl
        character(len=XLONGSTRLEN) :: micname
        integer :: imic
        if(iptcl < 1 .or. iptcl > self%os_ptcl2D%get_noris()) then
            write(logfhandle,*) 'iptcl : ',  iptcl
            write(logfhandle,*) 'nptcl2D: ', self%os_ptcl2D%get_noris()
            THROW_HARD('iptcl index out of range; get_micname')
        end if
        imic = self%os_ptcl2D%get_int(iptcl, 'stkind')
        if(imic < 1 .or. imic > self%os_mic%get_noris()) then
            write(logfhandle,*) 'imic : ', imic
            write(logfhandle,*) 'nmics: ', self%os_mic%get_noris()
            THROW_HARD('imic index out of range; get_micname')
        end if
        call self%os_mic%get_static(imic, 'intg', micname)
    end function get_micname

    module function get_mic_kind( self, imic ) result( mickind )
        class(sp_project), intent(inout) :: self
        integer,           intent(in)    :: imic
        type(string) :: mickind
        mickind = NIL
        if( self%os_mic%isthere(imic, 'imgkind') )then
            mickind = self%os_mic%get_str(imic, 'imgkind')
            return
        else
            if( self%os_mic%isthere(imic, 'movie') )then
                mickind = 'movie'
                return
            endif
        endif
    end function get_mic_kind

    module subroutine set_boxfile( self, i, boxfname, nptcls )
        class(sp_project), intent(inout) :: self
        integer,           intent(in)    :: i
        class(string),     intent(in)    :: boxfname
        integer, optional, intent(in)    :: nptcls
        type(nrtxtfile) :: boxfile
        integer         :: nptcls_here
        if( present(nptcls) )then
            nptcls_here = nptcls
            if( nptcls_here == 0 )then
                call self%os_mic%set(i, 'nptcls', 0)
                return
            endif
        else
            call boxfile%new(boxfname, 1)
            nptcls_here = boxfile%get_ndatalines()
            call boxfile%kill
        endif
        call self%os_mic%set(i, 'boxfile', boxfname)
        call self%os_mic%set(i, 'nptcls',  nptcls_here)
    end subroutine set_boxfile

    module integer function get_nmovies( self )
        class(sp_project), target, intent(inout) :: self
        type(string) :: imgkind
        integer :: i
        get_nmovies = 0
        do i=1,self%os_mic%get_noris()
            call self%os_mic%getter(i,'imgkind',imgkind)
            if( imgkind%to_char().eq.'movie' .or. imgkind%to_char().eq.'mic' ) get_nmovies = get_nmovies + 1
        enddo
    end function get_nmovies

    module integer function get_nintgs( self )
        class(sp_project), target, intent(inout) :: self
        integer :: i
        get_nintgs = 0
        do i=1,self%os_mic%get_noris()
            if( self%os_mic%isthere(i,'intg') )get_nintgs = get_nintgs + 1
        enddo
    end function get_nintgs

    module integer function get_nframes( self )
        class(sp_project), target, intent(inout) :: self
        integer :: i
        get_nframes = 0
        do i=1,self%os_mic%get_noris()
            if( self%os_mic%isthere(i,'frame') )get_nframes = get_nframes + 1
        enddo
    end function get_nframes

    module subroutine get_mic2stk_inds( self, mic2stk_inds, stk2mic_inds )
        class(sp_project),    intent(inout) :: self
        integer, allocatable, intent(inout) :: mic2stk_inds(:), stk2mic_inds(:)
        integer :: imic,istk,nmics,nstks,nptcls_mic,nptcls_stk,state_mic,state_stk
        if(allocated(mic2stk_inds))deallocate(mic2stk_inds)
        if(allocated(stk2mic_inds))deallocate(stk2mic_inds)
        nmics = self%os_mic%get_noris()
        nstks = self%os_stk%get_noris()
        if( nmics==0 .or. nstks==0 )then
            THROW_WARN('Empty fields! Fields need be populated; get_mic2stk_inds')
        endif
        if( nmics < nstks )then
            THROW_HARD('MIC & STK fileds indexing error 1! get_mic2stk_inds')
        endif
        allocate(mic2stk_inds(nmics),stk2mic_inds(nstks),source=0)
        if( nmics == nstks )then
            do imic = 1,nmics
                mic2stk_inds(imic) = imic
                stk2mic_inds(imic) = imic
            enddo
        else
            istk = 0
            do imic = 1,nmics
                nptcls_mic = self%os_mic%get_int(imic, 'nptcls')
                if( nptcls_mic > 0 )then
                    istk = istk+1
                    if( istk > nstks ) THROW_HARD('Too many stacks!  get_mic2stk_inds')
                else
                    ! micrographs without particles have no stack
                    cycle
                endif
                mic2stk_inds(imic) = istk
                stk2mic_inds(istk) = imic
            enddo
        endif
        ! state consistency
        do imic = 1,nmics
            istk = mic2stk_inds(imic)
            if( istk == 0 ) cycle
            nptcls_mic = self%os_mic%get_int(imic,'nptcls')
            nptcls_stk = self%os_stk%get_int(istk,'nptcls')
            if( nptcls_mic /= nptcls_stk )then
                print *, 'nptcls_mic ', nptcls_mic
                print *, 'nptcls_stk ', nptcls_stk
                THROW_WARN('Inconsistent number of particles!  get_mic2stk_inds')
            endif
            state_mic = self%os_mic%get_state(imic)
            state_stk = self%os_stk%get_state(istk)
            if( state_mic /= state_stk )then
                THROW_WARN('Inconsistent state!  get_mic2stk_inds')
            endif
        enddo
    end subroutine get_mic2stk_inds

end submodule simple_sp_project_mic
