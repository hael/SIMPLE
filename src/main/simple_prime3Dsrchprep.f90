! projection-matching based on Hadamard products, high-level search routines for PRIME3D
module simple_prime3Dsrchprep
!$ use omp_lib
!$ use omp_lib_kinds
use simple_math
use simple_oris,              only: oris
use simple_build,            only: build
use simple_params,           only: params
implicit none

type(oris), allocatable, target :: o_refs(:)
type(oris), allocatable, target :: o_peaks(:)
integer,    allocatable, target :: srch_order(:,:)
integer,    allocatable, target :: proj_space_inds(:,:)
logical,    allocatable, target :: state_exists(:)

contains

    subroutine prime3D_prepsrch( b, p )
        use simple_ran_tabu,  only: ran_tabu
        class(build),   intent(inout) :: b
        class(params),  intent(inout) :: p
        type(ran_tabu) :: rt
        real    :: euldists(p%nspace)
        integer :: i, istate, iproj, iptcl, prev_state, nnnrefs, cnt, prev_proj, prev_ref
        call prime3D_prepsrch_clean
        ! reference projections indices, here to avoid online allocation in prime3d_srch
        allocate(o_refs(p%fromp:p%top))
        do iptcl = p%fromp, p%top
            call o_refs(iptcl)%new_clean(p%nstates*p%nspace)
            cnt = 0
            do istate=1,p%nstates
                do iproj=1,p%nspace
                    cnt = cnt + 1
                    call o_refs(iptcl)%set( cnt, 'state', real(istate) ) ! Updates state
                    call o_refs(iptcl)%set( cnt, 'proj', real(iproj) )   ! Updates proj
                    call o_refs(iptcl)%set_euler( cnt, b%e%get_euler(iproj) )
                enddo
            enddo
        enddo
        ! projection direction peaks
        allocate(o_peaks(p%fromp:p%top))
        do iptcl = p%fromp, p%top
            call o_peaks(iptcl)%new_clean(p%npeaks)
            if( p%ctf.ne.'no' )then
                call o_peaks(iptcl)%set_all2single('kv', b%a%get(iptcl,'kv'))
                call o_peaks(iptcl)%set_all2single('cs', b%a%get(iptcl,'cs'))
                call o_peaks(iptcl)%set_all2single('fraca', b%a%get(iptcl,'fraca'))
                call o_peaks(iptcl)%set_all2single('dfx', b%a%get(iptcl,'dfx'))
                if( p%tfplan%mode .eq. 'astig' )then
                    call o_peaks(iptcl)%set_all2single('dfy', b%a%get(iptcl,'dfy'))
                    call o_peaks(iptcl)%set_all2single('angast', b%a%get(iptcl,'angast'))
                else
                    call o_peaks(iptcl)%set_all2single('dfy', b%a%get(iptcl,'dfx'))
                    call o_peaks(iptcl)%set_all2single('angast', 0.)
                endif
                if( p%tfplan%l_phaseplate )then
                    call o_peaks(iptcl)%set_all2single('phshift', b%a%get(iptcl,'phshift'))
                endif
            endif
        enddo
        ! reference projection directions indices
        allocate(proj_space_inds(p%fromp:p%top, p%nspace*p%nstates), source=0)
        ! nearest neighbours
        select case( trim(p%refine) )
            case( 'states' )
                allocate(srch_order(p%fromp:p%top, p%nnn), source=0)
                rt = ran_tabu(p%nnn)
                do iptcl = p%fromp, p%top
                    prev_state = nint(b%a%get(iptcl, 'state'))
                    srch_order(iptcl,:) = b%nnmat(iptcl,:) + (prev_state-1)*p%nspace
                    call rt%shuffle( srch_order(iptcl,:) )
                    prev_ref = (prev_state-1)*p%nspace + b%e%find_closest_proj(b%a%get_ori(iptcl))
                    call put_last(prev_ref, srch_order(iptcl,:))
                enddo
            case( 'neigh','shcneigh', 'greedyneigh' )
                nnnrefs =  p%nnn * p%nstates
                allocate(srch_order(p%fromp:p%top, nnnrefs), source=0)
                rt = ran_tabu(nnnrefs)
                do iptcl = p%fromp, p%top
                    prev_state = nint(b%a%get(iptcl, 'state'))
                    prev_proj  = b%e%find_closest_proj(b%a%get_ori(iptcl))
                    prev_ref   = (prev_state-1)*p%nspace + prev_proj
                    do istate = 0, p%nstates-1
                        i = istate*p%nnn+1
                        srch_order(iptcl,i:i+p%nnn-1) = b%nnmat(prev_proj,:) + istate*p%nspace
                    enddo
                    call rt%shuffle( srch_order(iptcl,:) )
                    call put_last(prev_ref, srch_order(iptcl,:))
                enddo                
            case( 'yes' )
                allocate(srch_order(p%fromp:p%top,p%nspace*p%nstates), source=0)
                do iptcl = p%fromp, p%top
                    prev_state = nint(b%a%get(iptcl, 'state'))
                    call b%e%calc_euldists(b%a%get_ori(iptcl), euldists)
                    call hpsort( p%nspace, euldists, srch_order(iptcl,:p%nspace) )
                    prev_ref = (prev_state-1)*p%nspace + srch_order(iptcl,1)
                    do istate = 0, p%nstates-1
                        i = istate*p%nspace+1
                        srch_order(iptcl,i:i+p%nspace-1) = srch_order(iptcl,1:p%nspace) + istate*p%nspace
                    enddo
                    call put_last(prev_ref, srch_order(iptcl,:))
                enddo
            case('no','shc','snhc','greedy')
                allocate(srch_order(p%fromp:p%top,p%nspace*p%nstates), source=0)
                rt = ran_tabu(p%nspace*p%nstates)
                do iptcl = p%fromp, p%top
                    prev_state = nint(b%a%get(iptcl, 'state'))
                    call rt%ne_ran_iarr( srch_order(iptcl,:) )
                    prev_ref = (prev_state-1)*p%nspace + b%e%find_closest_proj(b%a%get_ori(iptcl))
                    call put_last(prev_ref, srch_order(iptcl,:))
                enddo
            case DEFAULT
                stop 'Unknown refinement mode; simple_prime3D_srch; prep4srch'
        end select
        if( any(srch_order == 0) ) stop 'Invalid index in srch_order; simple_prime3d_srch::prep_ref_oris'
        call rt%kill
        ! states existence
        if( p%oritab.ne.'' )then
            state_exists = b%a%states_exist(p%nstates)
        else
            allocate(state_exists(p%nstates), source=.true.)
        endif
    end subroutine prime3D_prepsrch

    subroutine prime3D_prepsrch_clean
        integer :: i
        if(allocated(o_refs))then
            do i = 1, size(o_refs)
                call o_refs(i)%kill
            enddo
        endif
        if(allocated(o_peaks))then
            do i = 1, size(o_peaks)
                call o_peaks(i)%kill
            enddo
        endif
        if( allocated(srch_order) )deallocate(srch_order)
        if( allocated(proj_space_inds) )deallocate(proj_space_inds)
        if( allocated(state_exists) )deallocate(state_exists)
    end subroutine prime3D_prepsrch_clean

end module simple_prime3Dsrchprep
