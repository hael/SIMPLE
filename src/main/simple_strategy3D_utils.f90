module simple_strategy3D_utils
include 'simple_lib.f08'
use simple_strategy3D_alloc  ! singleton class s3D
use simple_strategy3D_srch,  only: strategy3D_srch
use simple_builder,          only: build_glob
use simple_parameters,       only: params_glob
use simple_polarft_corrcalc, only: pftcc_glob
implicit none

public :: extract_peak_ori, sort_corrs
private
#include "simple_local_flags.inc"

contains

    subroutine extract_peak_ori( s )
        use simple_ori, only: ori
        class(strategy3D_srch), intent(inout) :: s
        type(ori) :: osym, o_prev, o_new
        integer   :: ref, inpl, state, neff_states, proj
        real      :: shvec(2), shvec_incr(2), mi_state, euldist, dist_inpl, corr, mi_proj, frac
        logical   :: l_multistates
        ! stash previous ori
        call build_glob%spproj_field%get_ori(s%iptcl, o_prev)
        ! reference (proj)
        ref = s3D%proj_space_refinds_sorted(s%ithr, s%nrefsmaxinpl)
        if( ref < 1 .or. ref > s%nrefs ) THROW_HARD('ref index: '//int2str(ref)//' out of bound; extract_peak_ori')
        call build_glob%spproj_field%set(s%iptcl, 'proj', real(s3D%proj_space_proj(ref)))
        ! in-plane (inpl)
        inpl = s3D%proj_space_inplinds_sorted(s%ithr, s%nrefsmaxinpl)
        call build_glob%spproj_field%set(s%iptcl, 'inpl', real(inpl))
        ! Euler angle
        call build_glob%spproj_field%set_euler(s%iptcl, s3D%proj_space_euls(s%ithr,ref,inpl,1:3))
        ! shift
        shvec      = s%prev_shvec
        shvec_incr = 0.
        if( s%doshift ) then
            shvec_incr = s3D%proj_space_shift(s%ithr,ref,inpl,1:2)
            shvec      = shvec + shvec_incr
        end if
        where( abs(shvec) < 1e-6 ) shvec = 0.
        call build_glob%spproj_field%set_shift(s%iptcl, shvec)
        call build_glob%spproj_field%set_shift_incr(s%iptcl, shvec_incr)
        call build_glob%spproj_field%set(s%iptcl, 'shincarg', arg(shvec_incr))
        ! state
        state = 1
        l_multistates = s%nstates > 1
        if( l_multistates )then
            state = s3D%proj_space_state(ref)
            if( .not. s3D%state_exists(state) ) THROW_HARD('empty state: '//int2str(state)//'; extract_peak_ori')
        endif
        mi_state = 0.
        if( s%prev_state == state ) mi_state = 1.
        if( l_multistates )then
            call build_glob%spproj_field%set(s%iptcl, 'state',  real(state))
            call build_glob%spproj_field%set(s%iptcl, 'mi_state', mi_state)
        else
            call build_glob%spproj_field%set(s%iptcl, 'state',    1.)
            call build_glob%spproj_field%set(s%iptcl, 'mi_state', 1.)
        endif
        ! correlation
        corr = s3D%proj_space_corrs(s%ithr,ref,inpl)
        if( params_glob%cc_objfun /= OBJFUN_EUCLID )then
            if( corr < 0. ) corr = 0.
        end if
        call build_glob%spproj_field%set(s%iptcl, 'corr', corr)
        ! specscore
        call build_glob%spproj_field%set(s%iptcl, 'specscore', s%specscore)
        ! angular distances
        call build_glob%spproj_field%get_ori(s%iptcl, o_new)
        call build_glob%pgrpsyms%sym_dists(o_prev, o_new, osym, euldist, dist_inpl)
        if( build_glob%spproj_field%isthere(s%iptcl,'dist') )then
            call build_glob%spproj_field%set(s%iptcl, 'dist', 0.5*euldist + 0.5*build_glob%spproj_field%get(s%iptcl,'dist'))
        else
            call build_glob%spproj_field%set(s%iptcl, 'dist', euldist)
        endif
        call build_glob%spproj_field%set(s%iptcl, 'dist_inpl', dist_inpl)
        ! CONVERGENCE STATS
        ! projection direction overlap
        mi_proj  = 0.
        if( euldist < 0.5 ) mi_proj  = 1.
        call build_glob%spproj_field%set(s%iptcl, 'mi_proj', mi_proj)
        ! fraction of search space scanned
        neff_states = 1
        if( l_multistates ) neff_states = count(s3D%state_exists)
        if( s%neigh )then
            frac = 100.*real(s%nrefs_eval) / real(s%nnn * neff_states)
        else
            frac = 100.*real(s%nrefs_eval) / real(s%nprojs * neff_states)
        endif
        call build_glob%spproj_field%set(s%iptcl, 'frac',      frac)
        ! destruct
        call osym%kill
        call o_prev%kill
        call o_new%kill
    end subroutine extract_peak_ori

    subroutine sort_corrs( s )
        class(strategy3D_srch), intent(inout) :: s
        real    :: corrs(s%nrefs*NINPLPEAKS), corrs_highest(s%nrefs)
        integer :: proj_space_tmp(s%nrefs*NINPLPEAKS)
        integer :: i, j, arange(2), idx_array(s%nrefs*NINPLPEAKS)
        integer :: asequence(NINPLPEAKS)
        real    :: areal
        asequence = (/(j, j=1,NINPLPEAKS)/)
        do i = 1,s%nrefs
            arange(1) = (i-1)*NINPLPEAKS+1
            arange(2) = i*NINPLPEAKS
            s3D%proj_space_refinds_sorted        (s%ithr,arange(1):arange(2)) = i
            s3D%proj_space_inplinds_sorted       (s%ithr,arange(1):arange(2)) = asequence(:)
            s3D%proj_space_refinds_sorted_highest(s%ithr,                  i) = i
            if (.not. s3D%proj_space_corrs_srchd(s%ithr,i)) then
                corrs(arange(1):arange(2)) = -HUGE(areal)
                corrs_highest(i)           = -HUGE(areal)
            else
                corrs(arange(1):arange(2)) = s3D%proj_space_corrs(s%ithr,i,:)
                corrs_highest(i)           = s3D%proj_space_corrs(s%ithr,i,1)
            end if
        end do
        do j = 1,s%nrefs*NINPLPEAKS
            idx_array(j) = j
        end do
        call hpsort(corrs, idx_array)
        proj_space_tmp(:) = s3D%proj_space_refinds_sorted(s%ithr,:)
        do j = 1,s%nrefs*NINPLPEAKS
            s3D%proj_space_refinds_sorted(s%ithr,j) = proj_space_tmp(idx_array(j))
        end do
        proj_space_tmp(:) = s3D%proj_space_inplinds_sorted(s%ithr,:)
        do j = 1,s%nrefs*NINPLPEAKS
            s3D%proj_space_inplinds_sorted(s%ithr,j) = proj_space_tmp(idx_array(j))
        end do
        call hpsort(corrs_highest, s3D%proj_space_refinds_sorted_highest(s%ithr, :))
    end subroutine sort_corrs

end module simple_strategy3D_utils
