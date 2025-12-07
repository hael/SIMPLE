submodule(simple_sp_project) simple_sp_project_ptcl
implicit none
#include "simple_local_flags.inc"
contains

    ! static for OpenMP safety
    module subroutine map_ptcl_ind2stk_ind( self, oritype, iptcl, stkind, ind_in_stk )
        class(sp_project), target, intent(inout) :: self
        character(len=*),          intent(in)    :: oritype
        integer,                   intent(in)    :: iptcl
        integer,                   intent(out)   :: stkind
        integer,                   intent(out)   :: ind_in_stk
        class(oris), pointer                     :: ptcl_field
        real    :: smpd
        integer :: nptcls, fromp, top, box
        nullify(ptcl_field)
        ! set field pointer
        select case(trim(oritype))
            case('cls3D')
                call self%get_imginfo_from_osout(smpd, box, nptcls)
                if( iptcl < 1 .or. iptcl > nptcls )then
                    write(logfhandle,*) 'iptcl : ', iptcl
                    write(logfhandle,*) 'ncls: ',   nptcls
                    THROW_HARD('iptcl index out of range 1; map_ptcl_ind2stk_ind')
                endif
                stkind     = 1
                ind_in_stk = iptcl
                return
            case('ptcl2D')
                ptcl_field => self%os_ptcl2D
            case('ptcl3D')
                ptcl_field => self%os_ptcl3D
            case DEFAULT
                THROW_HARD('oritype: '//trim(oritype)//' not supported by map_ptcl_ind2stk_ind')
        end select
        nptcls = ptcl_field%get_noris()
        ! first sanity check, range
        if( iptcl < 1 .or. iptcl > nptcls )then
            write(logfhandle,*) 'iptcl : ', iptcl
            write(logfhandle,*) 'nptcls: ', nptcls
            THROW_HARD('iptcl index out of range 2; map_ptcl_ind2stk_ind')
        endif
        ! second sanity check, stack index present in ptcl_field
        if( .not. ptcl_field%isthere(iptcl, 'stkind') )then
            write(logfhandle,*) 'iptcl: ', iptcl
            THROW_HARD('stkind not present in field: '//trim(oritype)//'; map_ptcl_ind2stk_ind')
        endif
        stkind = ptcl_field%get_int(iptcl, 'stkind')
        if( ptcl_field%isthere(iptcl, 'indstk') )then
            ind_in_stk = ptcl_field%get_int(iptcl, 'indstk')
        else
            ! third sanity check, particle index in range
            fromp = self%os_stk%get_fromp(stkind)
            top   = self%os_stk%get_top(stkind)
            if( iptcl < fromp .or. iptcl > top )then
                write(logfhandle,*) 'iptcl            : ', iptcl
                write(logfhandle,*) 'stkind           : ', stkind
                write(logfhandle,*) 'prange for micstk: ', fromp, top
                THROW_HARD('iptcl index out of micstk range; map_ptcl_ind2stk_ind')
            endif
            ! output index in stack
            ind_in_stk = iptcl - fromp + 1
        endif
        ! cleanup
        nullify(ptcl_field)
    end subroutine map_ptcl_ind2stk_ind

    module subroutine map_cavgs_selection( self, states )
        class(sp_project), intent(inout) :: self
        integer,           intent(in)    :: states(:)
        integer, allocatable :: pinds(:)
        integer :: icls, sz_cls2D, sz_cls3D, sz_states, ncls, i, s
        sz_states = size(states)       
        sz_cls2D  = self%os_cls2D%get_noris()
        if( sz_cls2D /= sz_states )then
            write(logfhandle,*) 'size(cls2D): ', sz_cls2D
            write(logfhandle,*) 'sz_states  : ', sz_states
            THROW_HARD('size(cls2D) not consistent with size(states) in map_cavgs_selection, aborting...')
        endif
        ! map selection to self%os_cls2D
        do icls=1,sz_cls2D
            call self%os_cls2D%set(icls, 'state', states(icls))
        end do
        ! map selection to self%os_cls3D
        sz_cls3D = self%os_cls3D%get_noris()
        if( sz_cls3D /= sz_cls2D ) call self%os_cls3D%new(sz_cls2D, is_ptcl=.false.)
        sz_cls3D = sz_cls2D
        do icls=1,sz_cls3D
            call self%os_cls3D%set(icls, 'state', states(icls))
        end do
        ! map selection to self%os_ptcl2D & os_ptcl3D
        ncls = sz_states
        if( self%os_ptcl2D%get_noris() > 0 .and. self%os_ptcl3D%get_noris() > 0)then
            do icls=1,ncls
                call self%os_ptcl2D%get_pinds(icls, 'class', pinds)
                if( allocated(pinds) )then
                    s = states(icls)
                    do i=1,size(pinds)
                        call self%os_ptcl2D%set(pinds(i), 'state', s)
                        call self%os_ptcl3D%set(pinds(i), 'state', s)
                    end do
                    deallocate(pinds)
                endif
            end do
        endif
    end subroutine map_cavgs_selection

    ! map shifts obtained by cluster2D to the 3D field for cases where an external
    ! starting model is used for initializing 3D refinement (nanoparticles)
    module subroutine map2Dshifts23D( self )
        class(sp_project), intent(inout) :: self
        integer :: noris_ptcl3D, noris_ptcl2D
        real, allocatable :: shifts(:)
        if( self%is_virgin_field('ptcl2D') )then
            return
        endif
        ! ensure ptcl3D field congruent with ptcl2D field
        noris_ptcl3D = self%os_ptcl3D%get_noris()
        noris_ptcl2D = self%os_ptcl2D%get_noris()
        if( noris_ptcl3D /= noris_ptcl2D )then
            ! preserve defocus parameters, stack indices
            self%os_ptcl3D = self%os_ptcl2D
            call self%os_ptcl3D%delete_2Dclustering(keepshifts=.true., keepcls=.true.)
        else
            ! transfer shifts
            shifts = self%os_ptcl2D%get_all('x')
            call self%os_ptcl3D%set_all('x', shifts)
            deallocate(shifts)
            shifts = self%os_ptcl2D%get_all('y')
            call self%os_ptcl3D%set_all('y', shifts)
            deallocate(shifts)
        endif
    end subroutine map2Dshifts23D

    ! this map2ptcls routine assumes that any selection of class averages is done
    ! exclusively by state=0 flagging without any physical deletion
    module subroutine map2ptcls( self )
        class(sp_project), intent(inout) :: self
        integer, allocatable :: particles(:)
        type(ori) :: ori2d, ori_comp, o
        integer   :: ncls, icls, iptcl, pind, noris_ptcl3D, noris_ptcl2D
        real      :: corr, rproj, rstate
        if( self%is_virgin_field('cls3D') )then
            THROW_HARD('os_cls3D is virgin field; nothing to map back; map2ptcls')
        endif
        ! in case 2D was not done with SIMPLE but class averages imported from elsewhere
        if( self%os_ptcl2D%get_noris() == 0 ) return
        if( self%is_virgin_field('ptcl2D')  ) return
        ! ensure ptcl3D field congruent with ptcl2D field
        noris_ptcl3D = self%os_ptcl3D%get_noris()
        noris_ptcl2D = self%os_ptcl2D%get_noris()
        if( noris_ptcl3D /= noris_ptcl2D )then
            ! preserve defocus parameters, stack indices
            self%os_ptcl3D = self%os_ptcl2D
            call self%os_ptcl3D%delete_2Dclustering(keepcls=.true.)
        endif
        ! do the mapping
        ncls = self%os_cls3D%get_noris()
        do icls=1,ncls
            ! get particle indices
            call self%os_ptcl2D%get_pinds(icls, 'class', particles)
            if( allocated(particles) )then
                ! get 3d ori info
                call self%os_cls3D%get_ori(icls, o)
                rproj  = o%get('proj')
                rstate = o%get('state')
                corr   = o%get('corr')
                do iptcl=1,size(particles)
                    ! get particle index
                    pind  = particles(iptcl)
                    ! get 2d ori
                    call self%os_ptcl2D%get_ori(pind, ori2d)
                    ! transfer original parameters in self%os_ptcl2D
                    call self%os_ptcl2D%get_ori(pind, ori_comp)
                    ! compose ori3d and ori2d
                    call o%compose3d2d(ori2d, ori_comp)
                    ! update state in self%os_ptcl2D
                    call self%os_ptcl2D%set(pind, 'state', rstate)
                    ! set 3D orientation in self%os_ptcl3D
                    call self%os_ptcl3D%set_ori(pind, ori_comp)
                    ! set proj/state/corr
                    call self%os_ptcl3D%set(pind, 'proj',  rproj)
                    call self%os_ptcl3D%set(pind, 'state', rstate)
                    call self%os_ptcl3D%set(pind, 'corr',  corr)
                end do
                deallocate(particles)
            endif
        end do
        ! state = 0 all entries that don't have a state label
        do iptcl=1,noris_ptcl2D
            if( .not. self%os_ptcl2D%isthere(iptcl, 'state') ) call self%os_ptcl2D%set(iptcl, 'state', 0.)
            if( .not. self%os_ptcl3D%isthere(iptcl, 'state') ) call self%os_ptcl3D%set(iptcl, 'state', 0.)
        end do
        call ori2d%kill
        call ori_comp%kill
        call o%kill
    end subroutine map2ptcls

    ! this map2ptcls routine assumes that any selection of class averages is done
    ! exclusively by state=0 flagging without any physical deletion
    module subroutine map2ptcls_state( self, append, maxpop )
        class(sp_project), intent(inout) :: self
        logical, optional, intent(in)    :: append
        integer, optional, intent(in)    :: maxpop
        integer, allocatable             :: particles(:)
        integer :: ncls, icls, iptcl, i, nptcls, noris_ptcl3D, noris_ptcl2D
        integer :: maxnptcls, istate
        logical :: l_append
        noris_ptcl2D = self%os_ptcl2D%get_noris()
        if( noris_ptcl2D == 0 )then
            THROW_WARN('empty PTCL2D field. Nothing to do; map2ptcls_state')
            return
        endif
        l_append = .false.
        if(present(append)) l_append = append
        maxnptcls = huge(maxnptcls)
        if(present(maxpop)) maxnptcls = maxpop
        ! ensure ptcl3D field congruent with ptcl2D field
        noris_ptcl3D = self%os_ptcl3D%get_noris()
        if( noris_ptcl3D /= noris_ptcl2D )then
            ! preserve defocus parameters, stack indices
            self%os_ptcl3D = self%os_ptcl2D
            call self%os_ptcl3D%delete_2Dclustering(keepcls=.true.)
        endif
        ! undo previous selection if append is false & excludes non classified particles
        !$omp parallel do proc_bind(close) default(shared) private(iptcl)
        do iptcl=1,noris_ptcl2D
            if( .not.self%os_ptcl2D%isthere(iptcl, 'class') )then
                call self%os_ptcl2D%set_state(iptcl, 0)
                call self%os_ptcl3D%set_state(iptcl, 0)
            else if(.not. l_append) then
                call self%os_ptcl2D%set_state(iptcl, 1)
                call self%os_ptcl3D%set_state(iptcl, 1)
            endif
        end do
        !$omp end parallel do
        ! do the mapping
        ncls = self%os_cls2D%get_noris()
        do icls=1,ncls
            ! get particle indices
            call self%os_ptcl2D%get_pinds(icls, 'class', particles, l_shuffle=.true.)
            if( allocated(particles) )then
                istate = self%os_cls2D%get_state(icls)
                nptcls = size(particles)
                !$omp parallel do proc_bind(close) default(shared) private(i,iptcl)
                do i = 1,min(nptcls,maxnptcls)
                    iptcl = particles(i)
                    call self%os_ptcl2D%set_state(iptcl, istate)
                    call self%os_ptcl3D%set_state(iptcl, istate)
                enddo
                !$omp end parallel do
                if( maxnptcls < nptcls)then
                    !$omp parallel do proc_bind(close) default(shared) private(i,iptcl)
                    do i = maxnptcls+1,nptcls
                        iptcl = particles(i)
                        call self%os_ptcl2D%set_state(iptcl, 0)
                        call self%os_ptcl3D%set_state(iptcl, 0)
                    enddo
                    !$omp end parallel do
                endif
                deallocate(particles)
                if(present(maxpop)) call self%os_cls2D%set(icls, 'pop', nptcls)
            endif
        end do
        ! cls3D mirrors cls2D
        if( self%os_cls3D%get_noris() == ncls)then
            do icls=1,ncls
                call self%os_cls3D%set_state(icls, self%os_cls2D%get_state(icls))
            enddo
        else if( self%os_cls3D%get_noris() > 0 )then
            THROW_WARN('Inconsistent number of classes in cls2D & cls3D segments')
        endif
        ! state = 0 all entries that don't have a state/class label
        !$omp parallel do proc_bind(close) default(shared) private(iptcl)
        do iptcl = 1,noris_ptcl2D
            if( .not.self%os_ptcl2D%isthere(iptcl, 'state') .or.&
                &.not.self%os_ptcl3D%isthere(iptcl, 'state') )then
                 call self%os_ptcl2D%set_state(iptcl, 0)
                 call self%os_ptcl3D%set_state(iptcl, 0)
            endif
        end do
        !$omp end parallel do
    end subroutine map2ptcls_state

    ! report a real flag from the cls2D field back to the corresponding particles
    module subroutine map_cls2D_flag_to_ptcls( self, flag )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: flag
        integer, allocatable             :: particles(:)
        real    :: val
        integer :: icls, iptcl, nptcl3D, nptcl2D, ncls2D, pind
        nptcl2D = self%os_ptcl2D%get_noris()
        nptcl3D = self%os_ptcl3D%get_noris()
        ncls2D  = self%os_cls2D%get_noris()
        ! sanity
        if( nptcl2D == 0 )then
            THROW_WARN('empty PTCL2D field. Nothing to do; map_cls_flag_to_ptcls')
            return
        endif
        if( nptcl3D == 0 )then
            THROW_WARN('empty PTCL3D field. Nothing to do; map_cls_flag_to_ptcls')
            return
        endif
        if( ncls2D == 0 )then
            THROW_WARN('empty CLS2D field. Nothing to do; map_cls_flag_to_ptcls')
            return
        endif
        if( .not.self%os_cls2D%isthere(flag) )then
            THROW_WARN('flag is missing from the CLS2D field. Nothing to do; map_cls_flag_to_ptcls')
            return
        endif
        if( get_oriparam_ind(flag) == 0 )then
            THROW_WARN('flag is missing from the PTCL field. Nothing to do; map_cls_flag_to_ptcls')
            return
        endif
        ! do the class to particles mapping
        do icls = 1,ncls2D
            if( self%os_cls2D%get_state(icls)==0 )then
                ! making sure particles are turned off
                call self%os_ptcl2D%get_pinds(icls, 'class', particles)
                if( allocated(particles) )then
                    !$omp parallel do proc_bind(close) default(shared) private(iptcl,pind)
                    do iptcl=1,size(particles)
                        pind = particles(iptcl)
                        call self%os_ptcl2D%set_state(pind, 0)
                        call self%os_ptcl3D%set_state(pind, 0)
                    end do
                    !$omp end parallel do
                    deallocate(particles)
                endif
            else
                val = self%os_cls2D%get(icls, flag)
                ! get particle indices
                call self%os_ptcl2D%get_pinds(icls, 'class', particles)
                if( allocated(particles) )then
                    !$omp parallel do proc_bind(close) default(shared) private(iptcl,pind)
                    do iptcl=1,size(particles)
                        pind = particles(iptcl)
                        call self%os_ptcl2D%set(pind, flag, val)
                        call self%os_ptcl3D%set(pind, flag, val)
                    end do
                    !$omp end parallel do
                    deallocate(particles)
                endif
            endif
        end do
    end subroutine map_cls2D_flag_to_ptcls

    ! this updates cls fields with respect to ptcl2D/3D states
    module subroutine map_ptcls_state_to_cls( self )
        class(sp_project), intent(inout) :: self
        integer, allocatable :: cls_states(:), cls_pops(:)
        integer :: icls, iptcl, noris_ptcl2D, ncls2D, ncls3D
        noris_ptcl2D = self%os_ptcl2D%get_noris()
        if( noris_ptcl2D == 0 ) return
        ncls2D = self%os_cls2D%get_noris()
        if( ncls2D == 0 ) return
        ncls3D = self%os_cls2D%get_noris()
        ! cls2D
        cls_states = nint(self%os_cls2D%get_all('state'))
        cls_pops   = nint(self%os_cls2D%get_all('pop'))
        do iptcl=1,noris_ptcl2D
            if( self%os_ptcl2D%get_state(iptcl) == 1) cycle
            icls = self%os_ptcl2D%get_class(iptcl)
            if( (icls == 0) .or. (icls > ncls2D) ) cycle
            cls_pops(icls) = cls_pops(icls) - 1
        end do
        where( cls_pops < 1 )
            cls_pops   = 0
            cls_states = 0
        end where
        call self%os_cls2D%set_all('pop', real(cls_pops))
        call self%os_cls2D%set_all('state',real(cls_states))
        ! cls3D should be congruent
        if( ncls3D /= ncls2D ) call self%os_cls3D%new(ncls2D, is_ptcl=.false.)
        call self%os_cls3D%set_all('state',real(cls_states))
    end subroutine map_ptcls_state_to_cls

        ! Removes in place mics, stacks and particles with state=0
    ! new images are not generated & the indstk field is updated
    module subroutine prune_particles( self )
        class(sp_project), target, intent(inout) :: self
        type(oris)                :: os_ptcl2d, os_ptcl3d, os_stk, os_mic
        logical,      allocatable :: stks_mask(:), ptcls_mask(:)
        integer,      allocatable :: stkinds(:), stk2mic_inds(:), mic2stk_inds(:)
        integer                   :: iptcl, istk, stk_cnt, nptcls_tot, ptcl_cnt
        integer                   :: nstks, nstks_tot, fromp, top, fromp_glob, top_glob, nmics_tot
        integer                   :: stkind, ptcl_glob, nptcls_eff, indstk
        nstks_tot  = self%get_nstks()
        if( nstks_tot == 0 ) THROW_HARD('No particles to operate on!')
        ! particles reverse indexing
        nptcls_tot = self%os_ptcl2D%get_noris()
        allocate(ptcls_mask(nptcls_tot), stkinds(nptcls_tot))
        nptcls_eff = 0
        stkinds    = 0
        !$omp parallel do proc_bind(close) default(shared) private(iptcl) reduction(+:nptcls_eff)
        do iptcl = 1,nptcls_tot
            ptcls_mask(iptcl) = self%os_ptcl2D%get_state(iptcl) > 0
            if( ptcls_mask(iptcl) )then
                stkinds(iptcl) = self%os_ptcl2D%get_int(iptcl,'stkind')
                nptcls_eff     = nptcls_eff+1
            endif
        enddo
        !$omp end parallel do
        ! stacks
        allocate(stks_mask(nstks_tot))
        do istk = 1,nstks_tot
            stks_mask(istk) = self%os_stk%get_state(istk) > 0
            if( count(stkinds==istk) == 0 ) stks_mask(istk) = .false.
        enddo
        nstks = count(stks_mask)
        call os_stk%new(nstks, is_ptcl=.false.)
        ! mics
        nmics_tot = self%os_mic%get_noris()
        if( nmics_tot > 0 )then
            call self%get_mic2stk_inds(mic2stk_inds, stk2mic_inds)
            call os_mic%new(nstks, is_ptcl=.false.)
        endif
        ! removing deselected particles
        call os_ptcl2d%new(nptcls_eff, is_ptcl=.true.)
        call os_ptcl3d%new(nptcls_eff, is_ptcl=.true.)
        stkind     = 0
        stk_cnt    = 0
        top_glob   = 0
        ptcl_glob  = 0
        do istk = 1,nstks_tot
            if( .not.stks_mask(istk) ) cycle
            stk_cnt    = stk_cnt + 1
            stkind     = stkind  + 1
            fromp      = self%os_stk%get_fromp(istk)
            top        = self%os_stk%get_top(istk)
            fromp_glob = top_glob+1
            ptcl_cnt   = 0
            do iptcl = fromp,top
                if( .not.ptcls_mask(iptcl) )cycle
                ptcl_glob = ptcl_glob + 1
                top_glob  = top_glob+1
                ptcl_cnt  = ptcl_cnt+1
                indstk    = iptcl-fromp+1
                ! update orientations
                call os_ptcl2D%transfer_ori(ptcl_glob, self%os_ptcl2D, iptcl)
                call os_ptcl3D%transfer_ori(ptcl_glob, self%os_ptcl3D, iptcl)
                call os_ptcl2D%set(ptcl_glob,'stkind',stkind)
                call os_ptcl3D%set(ptcl_glob,'stkind',stkind)
                call os_ptcl2D%set(ptcl_glob,'indstk',indstk)
                call os_ptcl3D%set(ptcl_glob,'indstk',indstk)
            enddo
            ! update stack
            call os_stk%transfer_ori(stk_cnt, self%os_stk, istk)
            call os_stk%set(stk_cnt, 'fromp',  fromp_glob)
            call os_stk%set(stk_cnt, 'top',    top_glob)
            call os_stk%set(stk_cnt, 'nptcls', ptcl_cnt)
            ! update micrograph
            if( nmics_tot > 0 ) then
                call os_mic%transfer_ori(stk_cnt, self%os_mic, stk2mic_inds(istk))
                call os_mic%set(stk_cnt,'nptcls', ptcl_cnt)
            endif
        enddo
        self%os_stk    = os_stk
        self%os_mic    = os_mic
        self%os_ptcl2d = os_ptcl2D
        self%os_ptcl3d = os_ptcl3D
        call os_stk%kill
        call os_mic%kill
        call os_ptcl2d%kill
        call os_ptcl3d%kill
    end subroutine prune_particles

    module subroutine scale_projfile( self, smpd_target, new_projfile, cline, cline_scale, dir )
        ! this probably needs an oritype input for dealing with scale class averages
        class(sp_project),       intent(inout) :: self
        real,                    intent(inout) :: smpd_target
        class(string),           intent(out)   :: new_projfile
        class(cmdline),          intent(inout) :: cline
        class(cmdline),          intent(out)   :: cline_scale
        class(string), optional, intent(in)    :: dir
        type(string) :: projfile, projname, new_projname, str_tmp
        real         :: scale_factor, smpd_sc, smpd
        integer      :: box, box_sc, istk, n_os_stk
        n_os_stk = self%os_stk%get_noris()
        if( n_os_stk == 0 ) THROW_HARD('Empty stack object! scale_projfile')
        call self%projinfo%getter(1, 'projfile', projfile)
        call self%projinfo%getter(1, 'projname', projname)
        if( projname%is_blank() )then
            projname = get_fbody(projfile, string(METADATA_EXT), separator=.false.)
        endif
        ! dimensions
        smpd = self%get_smpd()
        box  = self%get_box()
        call autoscale(box, smpd, smpd_target, box_sc, smpd_sc, scale_factor)
        call cline_scale%set('prg',      'scale_project')
        call cline_scale%set('scale',    scale_factor)
        call cline_scale%set('projfile', projfile%to_char())
        call cline_scale%set('smpd',     smpd_sc)
        str_tmp = cline%get_carg('mkdir')
        if( cline%defined('mkdir') )     call cline_scale%set('mkdir',str_tmp%to_char())
        if( present(dir) )               call cline_scale%set('dir_target',dir%to_char()//path_separator)
        if( box == box_sc )then
            ! no scaling
            new_projfile = projfile
            return
        endif
        ! parameter updates
        do istk = 1,n_os_stk
            call self%os_stk%set(istk, 'smpd', smpd_sc)
            call self%os_stk%set(istk, 'box',  box_sc)
        enddo
        call self%os_ptcl2D%mul_shifts(scale_factor)
        call self%os_ptcl3D%mul_shifts(scale_factor)
        ! name changes and list for scaling job
        new_projname = projname%to_char()//SCALE_SUFFIX
        new_projfile = new_projname%to_char()//METADATA_EXT
        call cline%set('projname', new_projname)
        call cline%delete('projfile')
        call self%update_projinfo(cline)
        if(present(dir))then
            call self%add_scale_tag(dir=dir//path_separator)
        else
            call self%add_scale_tag
        endif
        ! save
        call self%write
        ! command line for scaling
        call cline_scale%set('newbox', box_sc)
        if( cline%defined('state') )  call cline_scale%set('state',  cline%get_iarg('state'))
        if( cline%defined('nthr') )   call cline_scale%set('nthr',   cline%get_iarg('nthr'))
        if( cline%defined('nparts') ) call cline_scale%set('nparts', cline%get_iarg('nparts'))
    end subroutine scale_projfile

    ! Getters/Setters

    module subroutine set_ptcl2D_thumb( self, projfile, indices, boxsize )
        class(sp_project), intent(inout) :: self
        class(string),     intent(in)    :: projfile
        integer,           intent(in)    :: indices(:)
        integer,           intent(out)   :: boxsize
        type(string)          :: thumbfile, tmpfile
        type(string)          :: stkname
        integer, allocatable  :: arr(:)
        real,    allocatable  :: sort_vals(:)
        type(image)           :: stkimg
        type(stack_io)        :: stkio_w
        integer               :: iori, iidx, idx, ldim_stk(3), nptcls
        real                  :: smpd
        logical               :: first = .true.
        thumbfile = stemname(projfile) // "/thumbptcl2D.jpeg"
        ! always recreate
        if(file_exists(thumbfile)) call del_file(thumbfile)
        tmpfile = stemname(projfile) // "/thumbptcl2D.mrcs"
        if(file_exists(tmpfile)) call del_file(tmpfile)
        do iidx=1, size(indices)
            iori = indices(iidx)
            call self%get_stkname_and_ind('ptcl2D', iori, stkname, idx)
            call find_ldim_nptcls(stkname, ldim_stk, nptcls, smpd)
            call stkimg%new([ldim_stk(1), ldim_stk(1), 1], smpd)
            if(first) then
                call stkio_w%open(tmpfile, smpd, 'write', box=ldim_stk(1))
                boxsize = ldim_stk(1)
                first = .false.
            end if
            call stkimg%read(stkname, idx)
            call stkimg%fft()
            call stkimg%bpgau2D(0.0, 4.0)
            call stkimg%ifft()
            call stkio_w%write(iidx, stkimg)
            call stkimg%kill()
            call self%os_ptcl2D%set(iori, "thumb",    thumbfile)
        end do
        call stkio_w%close()
        call mrc2jpeg_tiled(tmpfile, thumbfile)
        call del_file(tmpfile)
        if(allocated(sort_vals)) deallocate(sort_vals)
        if(allocated(arr))       deallocate(arr)
        call thumbfile%kill
        call tmpfile%kill
    end subroutine set_ptcl2D_thumb

    module logical function is_virgin_field( self, oritype )
        class(sp_project), target, intent(inout) :: self
        character(len=*),          intent(in)    :: oritype
        class(oris), pointer :: os
        integer :: i, n
        nullify(os)
        is_virgin_field = .false.
        ! set field pointer
        select case(trim(oritype))
            case('ptcl2D')
                os => self%os_ptcl2D
            case('ptcl3D')
                os => self%os_ptcl3D
            case('cls3D')
                os => self%os_cls3D
            case DEFAULT
                THROW_HARD('oritype: '//trim(oritype)//' not supported by is_virgin_field')
        end select
        n = os%get_noris()
        if( n == 0 )then
            THROW_WARN('cannot check virginity of non-existent field (touched for the very first time???)')
            return
        endif
        do i=1,n
            if( os%has_been_searched(i) )return
        enddo
        is_virgin_field = .true.
    end function is_virgin_field

    module integer function count_state_gt_zero( self )
        class(sp_project), target, intent(inout) :: self
        integer :: iori, cnt_s_gt_zero_ptcl2D, cnt_s_gt_zero_ptcl3D
        ! # ptcl2D = # ptcl3D
        if( self%os_ptcl2D%get_noris() /= self%os_ptcl3D%get_noris() )then
            THROW_HARD('Inconsistent number of particles in PTCL2D/PTCL3D segments; count_state_gt_zero')
        endif
        ! check ptcl2D/ptcl3D fields
        cnt_s_gt_zero_ptcl2D = 0
        cnt_s_gt_zero_ptcl3D = 0
        do iori = 1,self%os_ptcl2D%get_noris()
            if( .not. self%os_ptcl2D%isthere(iori,'state') )then
                THROW_HARD('state flag missing from self%os_ptcl2D; count_state_gt_zero')
            endif
             if( .not. self%os_ptcl3D%isthere(iori,'state') )then
                THROW_HARD('state flag missing from self%os_ptcl3D; count_state_gt_zero')
            endif
            if( self%os_ptcl2D%get_state(iori) > 0 )then
                cnt_s_gt_zero_ptcl2D = cnt_s_gt_zero_ptcl2D + 1
            endif
            if( self%os_ptcl3D%get_state(iori) > 0 )then
                cnt_s_gt_zero_ptcl3D = cnt_s_gt_zero_ptcl3D + 1
            endif
        enddo
        if( cnt_s_gt_zero_ptcl2D == cnt_s_gt_zero_ptcl3D )then
            count_state_gt_zero = cnt_s_gt_zero_ptcl2D
        else
            THROW_HARD('state labelling incosistent between PTCL2D/PTCL3D segments')
        endif 
    end function count_state_gt_zero

    module integer function get_nptcls( self )
        class(sp_project), target, intent(in) :: self
        integer :: i, nos
        get_nptcls = 0
        nos        = self%os_stk%get_noris()
        if( nos == 0 )return
        do i=1,nos
            get_nptcls = get_nptcls + self%os_stk%get_int(i,'nptcls')
        enddo
        ! sanity check
        if( self%os_stk%isthere(nos,'top') )then
            if( self%os_stk%get_top(nos) /=  get_nptcls )then
                write(logfhandle,*) 'nptcls from ptcls', get_nptcls
                write(logfhandle,*) 'nptcls from top  ', self%os_stk%get_top(nos)
                THROW_HARD('total # particles .ne. last top index; get_nptcls')
            endif
        endif
    end function get_nptcls

    module subroutine get_boxcoords( self, iptcl, coords )
        class(sp_project), target, intent(in)  :: self
        integer,                   intent(in)  :: iptcl
        integer,                   intent(out) :: coords(2)
        coords(1) = self%os_ptcl2D%get_int(iptcl,'xpos')
        coords(2) = self%os_ptcl2D%get_int(iptcl,'ypos')
    end subroutine get_boxcoords

    subroutine set_boxcoords( self, iptcl, coords )
        class(sp_project), intent(inout) :: self
        integer,           intent(in)    :: iptcl, coords(2)
        call self%os_ptcl2D%set(iptcl,'xpos',coords(1))
        call self%os_ptcl2D%set(iptcl,'ypos',coords(2))
    end subroutine set_boxcoords

end submodule simple_sp_project_ptcl
