!@descr: mirrored polar references, state/bank isolation and FFT scores against analytic and brute-force oracles
submodule (simple_polarft_calc) simple_polarft_mirr_tester
use simple_sym,          only: sym
use simple_test_utils,   only: assert_real, assert_int, assert_true, set_fixed_seed
implicit none
character(len=2), parameter :: GROUPS(*) = [character(len=2) :: 'c1','c3','d2','d7','t','o','i']

contains

    module subroutine run_all_polarft_mirr_tests()
        integer :: saved_team, saved_global, nteam, actual_team, i
        saved_team = 1
        nteam = 1
        !$ saved_team = omp_get_max_threads()
        !$ nteam = min(omp_get_num_procs(),max(2,saved_team))
        saved_global = nthr_glob
        !$ call omp_set_num_threads(nteam)
        actual_team = 1
        !$omp parallel default(shared)
        !$omp single
        !$ actual_team = omp_get_num_threads()
        !$omp end single
        !$omp end parallel
        call assert_int(nteam,actual_team,'mirroring regression uses the requested team')
        call check_frame_angles(12,nteam)
        call check_frame_angles(13,nteam)
        do i = 1,size(GROUPS)
            call check_memo_banks(GROUPS(i),2,1,nteam)
            call check_memo_banks(GROUPS(i),3,1,nteam)
            call check_memo_banks(GROUPS(i),3,2,nteam)
            call check_projections(GROUPS(i),nteam)
        enddo
        nthr_glob = saved_global
        !$ call omp_set_num_threads(saved_team)
    end subroutine run_all_polarft_mirr_tests

    subroutine init_params(p,cline,nspace,nstates,nteam)
        class(parameters), intent(inout) :: p
        class(cmdline),    intent(inout) :: cline
        integer, intent(in) :: nspace,nstates,nteam
        call cline%set('box',32)
        call cline%set('smpd',1.)
        call cline%set('mskdiam',24.)
        call cline%set('nspace',nspace)
        call cline%set('nstates',nstates)
        call cline%set('nptcls',2)
        call cline%set('pftsz',12)
        call cline%set('nthr',nteam)
        call cline%set('ctf','no')
        call cline%set('objfun','euclid')
        call p%new(cline,silent=.true.)
        call set_fixed_seed(20261108)
    end subroutine init_params

    !> Angular functions have the real-image Hermitian half-turn relation.
    complex(sp) function angular_field(theta,k,bank,idir,half) result(value)
        real(dp), intent(in) :: theta
        integer, intent(in) :: k,bank,idir,half
        real(dp) :: amplitude,phase
        amplitude = (0.8_dp+0.11_dp*bank+0.03_dp*idir+0.07_dp*half)*exp(-0.08_dp*k) &
            *(1._dp+0.15_dp*cos(2*theta)+0.09_dp*sin(2*theta))
        phase = (0.13_dp*bank+0.04_dp*half*k)*sin(theta)+0.07_dp*idir*cos(theta)
        value = cmplx(amplitude*cmplx(cos(phase),sin(phase),dp),kind=sp)
    end function angular_field

    subroutine check_memo_banks(pgrp,nstates,nmodels,nteam)
        character(len=*), intent(in) :: pgrp
        integer, intent(in) :: nstates,nmodels,nteam
        integer, parameter :: NSPACE=12, KS(2)=[2,6], NP=12, NR=2*NP
        class(parameters), allocatable, target :: p
        class(cmdline), allocatable :: cline
        class(polarft_calc), allocatable :: calc
        class(oris), allocatable :: dirs
        class(sym), allocatable :: group
        complex(sp), allocatable :: refs(:,:,:,:),particles(:,:,:)
        real(sp), allocatable :: ctfs(:,:,:), expected(:,:,:,:,:), actual(:,:,:,:,:), serial(:,:,:,:,:)
        real(sp), allocatable, target :: sigma2(:,:)
        real(sp) :: shifts(2,2),psi,factor
        real(dp) :: theta,expected_euclid(NR),cc(NR)
        integer :: nrefs,bank,local,source,half,ir,k,ref,pdim(3),pass,s,obj
        logical :: mate
        allocate(p,cline,calc,dirs,group)
        call init_params(p,cline,NSPACE,nstates,nteam)
        call group%new(pgrp)
        call dirs%new(NSPACE,is_ptcl=.false.)
        call group%build_refspiral(dirs)
        nrefs = NSPACE*nstates*nmodels
        call calc%new(p,nrefs,[1,2],KS)
        pdim = calc%get_pdim_interp()
        allocate(refs(NP,KS(1):pdim(3),nrefs,2),particles(NP,KS(1):pdim(3),2),ctfs(NP,KS(1):pdim(3),2))
        allocate(sigma2(KS(1):KS(2),2),expected(NR,nrefs,2,2,2),actual(NR,nrefs,2,2,2),serial(NR,nrefs,2,2,2))
        do half = 1,2
            do k = KS(1),pdim(3)
                do ir = 1,NP
                    theta = DPI*real(ir-1,dp)/NP
                    particles(ir,k,half) = angular_field(theta,k,0,2,half)
                    ctfs(ir,k,half) = real(1._dp+0.1_dp*cos(2*theta)+0.03_dp*half*sin(2*theta),sp)
                enddo
            enddo
            sigma2(:,half) = [(1.+0.12*half*k,k=KS(1),KS(2))]
            call calc%set_ptcl_pft(half,particles(:,:,half))
            call calc%set_ptcl_ctfmat(half,ctfs(:,:,half))
            call calc%set_eo(half,half==1)
            do ref = 1,nrefs
                bank = (ref-1)/NSPACE+1
                local = mod(ref-1,NSPACE)+1
                mate = local > NSPACE/2
                source = local
                psi = 0.
                if( mate )then
                    source = dirs%get_int(local,'mirr')
                    psi = dirs%get(local,'psi')
                endif
                do k = KS(1),pdim(3)
                    do ir = 1,NP
                        theta = DPI*real(ir-1,dp)/NP
                        if( mate ) theta = DPI-theta+real(psi,dp)*DPI/180._dp
                        refs(ir,k,ref,half) = angular_field(theta,k,bank,source,half)
                    enddo
                enddo
                call calc%set_ref_pft(ref,refs(:,:,ref,half),half==1)
            enddo
        enddo
        call calc%assign_sigma2_noise(sigma2)
        call calc%memoize_ptcls()
        shifts(:,1) = 0.
        shifts(:,2) = [0.31,-0.27]
        ! Ordinary FFTs, generic no-geometry fallback, serial mirroring, and repeated forced teams.
        do pass = 0,4
            !$ call omp_set_num_threads(1)
            if( pass >= 3 )then
                !$ call omp_set_num_threads(nteam)
            endif
            p%mirr_proj = 'yes'
            if( pass == 0 ) p%mirr_proj = 'no'
            ! Change the bank between routes: an untouched later-state memo must be stale.
            factor = 1.+0.19*min(pass,2)
            do half = 1,2
                do ref = 1,nrefs
                    call calc%set_ref_pft(ref,factor*refs(:,:,ref,half),half==1)
                    do s = 1,2
                        call brute_scores(factor*refs(:,:,ref,half),particles(:,:,half),ctfs(:,:,half), &
                            sigma2(:,half),shifts(:,s),expected_euclid,cc)
                        expected(:,ref,half,s,1) = real(expected_euclid,sp)
                        expected(:,ref,half,s,2) = real(cc,sp)
                    enddo
                enddo
            enddo
            if( pass <= 1 )then
                call calc%memoize_refs()
            else
                call calc%memoize_refs(dirs)
            endif
            do s = 1,2
                p%cc_objfun = OBJFUN_EUCLID
                !$omp parallel do collapse(2) default(shared) private(half,ref) schedule(static)
                do half = 1,2
                    do ref = 1,nrefs
                        call calc%gen_euclid_dists(ref,half,shifts(:,s),actual(:,ref,half,s,1))
                    enddo
                enddo
                !$omp end parallel do
                p%cc_objfun = OBJFUN_CC
                !$omp parallel do collapse(2) default(shared) private(half,ref) schedule(static)
                do half = 1,2
                    do ref = 1,nrefs
                        call calc%gen_objfun_vals(ref,half,shifts(:,s),actual(:,ref,half,s,2))
                    enddo
                enddo
                !$omp end parallel do
            enddo
            do obj = 1,2
                call assert_real(0.,maxval(abs(actual(:,:,:,:,obj)-expected(:,:,:,:,obj))),2.e-5, &
                    'mirror banks match brute angular scores '//pgrp)
            enddo
            if( pass == 2 ) serial = actual
            if( pass >= 3 ) call assert_true(all(actual==serial),'mirror banks preserve serial/team score identity '//pgrp)
        enddo
        if( pgrp == 'c1' .and. nstates == 2 .and. nmodels == 1 ) call check_bank_particle_frame(calc,refs,particles,ctfs,sigma2,1.+0.19*2)
        call calc%kill()
        call dirs%kill()
        call group%kill()
        call cline%kill()
    end subroutine check_memo_banks

    !> The rows of the shifted-particle bank (the 2D in-plane insertion's row: the particle shifted once by its shift
    !! in its own frame, then every rotation of the reference) against two oracles, with a CTF that depends on the
    !! angle (the astigmatic form 1 + a cos 2 theta + b sin 2 theta) and per-ring noise weights: the direct angular
    !! sums of a particle multiplied by exp(-i phase) at its own samples' angles and the reference unshifted; and the
    !! shifted-reference kernel at the shift R(theta_r) u of each rotation r (the frame rule of the polish), whose
    !! row at r is the bank's at r. The tolerance is that of the existing kernel against its oracle (single-precision
    !! FFT sums). Run on the references as the last scoring pass left them (scaled by factor)
    subroutine check_bank_particle_frame(calc,refs,particles,ctfs,sigma2,factor)
        class(polarft_calc), intent(inout) :: calc
        complex(sp), intent(in) :: refs(:,:,:,:),particles(:,:,:)
        real(sp), intent(in) :: ctfs(:,:,:),sigma2(:,:),factor
        integer, parameter :: NP=12, NR=2*NP, NREFS_CHECK=3, IREFS(NREFS_CHECK)=[1,7,13]
        real(sp), parameter :: TOL=5.e-5
        complex(sp), allocatable :: bank(:,:)
        real(sp) :: us(2,2),u(2),row(NR),rowref(NR),sref(2),mat(2,2),worst_oracle,worst_frame,dang
        real(dp) :: oracle(NR)
        integer :: half,iref,ius,r,kfromto(2)
        calc%p_ptr%cc_objfun = OBJFUN_EUCLID
        kfromto = calc%get_kfromto()
        allocate(bank(NP+1,kfromto(2)-kfromto(1)+1))
        us(:,1) = [0.31,-0.27]
        us(:,2) = [1.3,-0.9]
        dang = 360./real(NR)
        worst_oracle = 0.
        worst_frame  = 0.
        do half = 1,2
            do ius = 1,2
                u = us(:,ius)
                call calc%memoize_shifted_ptcl(half,u,bank)
                do iref = 1,NREFS_CHECK
                    call calc%gen_euclid_dists_bank(IREFS(iref),half,bank,row)
                    call brute_particle_frame(factor*refs(:,:,IREFS(iref),half),particles(:,:,half),ctfs(:,:,half), &
                        sigma2(:,half),u,oracle)
                    worst_oracle = max(worst_oracle,maxval(abs(row-real(oracle,sp))))
                    do r = 1,NR
                        call rotmat2d(real(r-1)*dang,mat)
                        sref(1) = mat(1,1)*u(1)+mat(1,2)*u(2)
                        sref(2) = mat(2,1)*u(1)+mat(2,2)*u(2)
                        call calc%gen_euclid_dists(IREFS(iref),half,sref,rowref)
                        worst_frame = max(worst_frame,abs(row(r)-rowref(r)))
                    enddo
                enddo
            enddo
        enddo
        call assert_real(0.,worst_oracle,TOL,'the shifted-particle bank equals the direct sums of the particle shifted in its own frame')
        call assert_real(0.,worst_frame,TOL,'the shifted-particle bank at u equals the shifted-reference kernel at R(theta_r) u at every rotation')
        deallocate(bank)
    end subroutine check_bank_particle_frame

    !> The loss of the particle shifted by shift in its own frame (its samples multiplied by exp(-i phase) at their own
    !! angles) against the unshifted reference at every rotation, as brute_scores sums it
    subroutine brute_particle_frame(ref,particle,ctf,sigma2,shift,expected_euclid)
        complex(sp), intent(in) :: ref(:,2:),particle(:,2:)
        real(sp), intent(in) :: ctf(:,2:),sigma2(2:),shift(2)
        real(dp), intent(out) :: expected_euclid(:)
        complex(dp) :: r,p
        real(dp) :: theta,phase,weight,loss,pnorm
        integer :: np,nr,rot,ir,j,k
        np = size(ref,1)
        nr = 2*np
        do rot = 0,nr-1
            loss=0; pnorm=0
            do k = 2,6
                weight = real(k,dp)/real(sigma2(k),dp)
                do ir = 1,np
                    j = modulo(ir-1-rot,nr)+1
                    if( j <= np )then
                        r = cmplx(ref(j,k),kind=dp)
                    else
                        r = conjg(cmplx(ref(j-np,k),kind=dp))
                    endif
                    r = r*real(ctf(ir,k),dp)
                    theta = DPI*real(ir-1,dp)/np
                    phase = 2*DPI*k/32*(sin(theta)*shift(1)-cos(theta)*shift(2))
                    p = cmplx(particle(ir,k),kind=dp)*cmplx(cos(phase),-sin(phase),dp)
                    loss = loss+weight*abs(p-r)**2
                    pnorm = pnorm+weight*abs(p)**2
                enddo
            enddo
            expected_euclid(rot+1) = loss/pnorm
        enddo
    end subroutine brute_particle_frame

    !> Direct angular sums, without FFTs or memoised reference terms.
    subroutine brute_scores(ref,particle,ctf,sigma2,shift,expected_euclid,cc)
        complex(sp), intent(in) :: ref(:,2:),particle(:,2:)
        real(sp), intent(in) :: ctf(:,2:),sigma2(2:),shift(2)
        real(dp), intent(out) :: expected_euclid(:),cc(:)
        complex(dp) :: r,p
        real(dp) :: theta,phase,weight,loss,pnorm,cross,rpow,ppow
        integer :: np,nr,rot,ir,j,k
        np = size(ref,1)
        nr = 2*np
        do rot = 0,nr-1
            loss=0; pnorm=0; cross=0; rpow=0; ppow=0
            do k = 2,6
                weight = real(k,dp)/real(sigma2(k),dp)
                do ir = 1,np
                    j = modulo(ir-1-rot,nr)+1
                    if( j <= np )then
                        r = cmplx(ref(j,k),kind=dp)
                    else
                        r = conjg(cmplx(ref(j-np,k),kind=dp))
                    endif
                    theta = DPI*real(j-1,dp)/np
                    phase = 2*DPI*k/32*(sin(theta)*shift(1)-cos(theta)*shift(2))
                    r = r*cmplx(cos(phase),sin(phase),dp)*real(ctf(ir,k),dp)
                    p = cmplx(particle(ir,k),kind=dp)
                    loss = loss+weight*abs(p-r)**2
                    pnorm = pnorm+weight*abs(p)**2
                    cross = cross+real(p*conjg(r),dp)
                    rpow = rpow+abs(r)**2
                    ppow = ppow+abs(p)**2
                enddo
            enddo
            expected_euclid(rot+1) = loss/pnorm
            cc(rot+1) = cross/sqrt(rpow*ppow)
        enddo
    end subroutine brute_scores

    !> Finite angular harmonics, including the even/odd half-ring Nyquist mode.
    complex(dp) function finite_field(theta,np) result(value)
        real(dp), intent(in) :: theta
        integer, intent(in) :: np
        value = cmplx(1._dp+0.15_dp*cos(2*theta)+0.09_dp*sin(2*theta), &
            0.22_dp*sin(theta)+0.13_dp*cos(3*theta),dp)
        if( mod(np,2) == 0 )then
            value = value+0.02_dp*cos(np*theta)
        else
            value = value+cmplx(0._dp,0.02_dp*cos(np*theta),dp)
        endif
    end function finite_field

    !> General-angle resampling agrees with closed forms without changing the source directions.
    subroutine check_frame_angles(np,nteam)
        integer, intent(in) :: np,nteam
        integer, parameter :: NSPACE=12, KS(2)=[2,6]
        integer, parameter :: PARTNERS(6)=[12,10,8,11,7,9]
        real, parameter :: ANGLES(6)=[0.,180.,17.3,73.7429,273.5132,359.95]
        class(parameters), allocatable, target :: p
        class(cmdline), allocatable :: cline
        class(polarft_calc), allocatable :: calc
        complex(sp), allocatable :: refs(:,:,:),expected(:,:,:),serial(:,:,:)
        real(dp) :: theta,scale
        integer :: pass,i,k,ir
        allocate(p,cline,calc)
        call init_params(p,cline,NSPACE,1,nteam)
        p%pftsz=np
        call calc%new(p,NSPACE,[1,2],KS)
        allocate(refs(np,5,NSPACE),expected(np,5,NSPACE),serial(np,5,NSPACE))
        do i = 1,6
            do k = 1,5
                scale = 1._dp+0.07_dp*i+0.03_dp*k
                do ir = 1,np
                    theta = DPI*real(ir-1,dp)/np
                    expected(ir,k,i)=cmplx(scale*finite_field(theta,np),kind=sp)
                    theta = DPI-theta+real(ANGLES(i),dp)*DPI/180._dp
                    expected(ir,k,PARTNERS(i))=cmplx(scale*finite_field(theta,np),kind=sp)
                enddo
            enddo
        enddo
        do pass = 0,2
            !$ call omp_set_num_threads(1)
            if( pass > 0 )then
                !$ call omp_set_num_threads(nteam)
            endif
            refs=cmplx(99.,-71.,sp)
            refs(:,:,:6)=expected(:,:,:6)
            call calc%mirror_ref_bank(refs,PARTNERS,ANGLES)
            call assert_real(0.,maxval(abs(refs-expected)),2.e-6,'general mirror frames match angular harmonics and Nyquist mode')
            if( pass == 0 ) serial=refs
            if( pass > 0 ) call assert_true(all(refs==serial),'general mirror frames preserve serial/team identity')
        enddo
        call calc%kill()
        call cline%kill()
    end subroutine check_frame_angles

    subroutine check_projections(pgrp,nteam)
        use simple_projector, only: projector
        character(len=*), intent(in) :: pgrp
        integer, intent(in) :: nteam
        integer, parameter :: BOX=32, NSPACE=12, NSTATES=3, KS(2)=[2,7]
        class(parameters), allocatable, target :: p
        class(cmdline), allocatable :: cline
        class(polarft_calc), allocatable :: calc
        class(projector), allocatable :: padded
        class(image), allocatable :: volume
        class(oris), allocatable :: dirs
        class(sym), allocatable :: group
        complex(sp), allocatable :: direct(:,:,:,:),serial(:,:,:,:),section(:,:),poison(:,:),base_exp(:,:,:)
        real(dp), allocatable :: centres(:,:),widths(:),amplitudes(:)
        real, pointer :: vox(:,:,:)
        real :: base(3,3),mat(3,3),scale
        real(dp) :: xyz(3),q(3),theta,phase,a,gap,interp,refnorm
        complex(dp) :: exact
        integer :: nsym,ig,ib,n,i,j,k,ir,ref,state,half,pass,pdim(3)
        allocate(p,cline,calc,padded,volume,dirs,group)
        call init_params(p,cline,NSPACE,NSTATES,nteam)
        call group%new(pgrp)
        call dirs%new(NSPACE,is_ptcl=.false.)
        call group%build_refspiral(dirs)
        call calc%new(p,NSPACE*NSTATES,[1,2],KS)
        pdim = calc%get_pdim_interp()
        allocate(section(12,KS(1):KS(2)),poison(12,KS(1):pdim(3)), &
            direct(12,KS(1):KS(2),NSPACE*NSTATES,2),serial(12,KS(1):KS(2),NSPACE*NSTATES,2))
        nsym = group%get_nsym()
        allocate(centres(3,3*nsym),widths(3*nsym),amplitudes(3*nsym))
        base(:,1) = BOX*[0.11,0.04,-0.08]
        base(:,2) = BOX*[-0.10,0.07,0.03]
        base(:,3) = BOX*[0.02,-0.13,0.09]
        do ig = 1,nsym
            call group%get_sym_rmat(ig,mat)
            do ib = 1,3
                n = 3*(ig-1)+ib
                centres(:,n) = real(matmul(mat,base(:,ib)),dp)
                widths(n) = 1.3_dp+0.2_dp*ib
                amplitudes(n) = (1.2_dp-0.2_dp*ib)/nsym
            enddo
        enddo
        call volume%new([BOX,BOX,BOX],1.,wthreads=.false.)
        call volume%get_rmat_ptr(vox)
        do k = 1,BOX
            do j = 1,BOX
                do i = 1,BOX
                    xyz = real([i,j,k]-BOX/2-1,dp)
                    a=0
                    do n = 1,size(widths)
                        a = a+amplitudes(n)*exp(-sum((xyz-centres(:,n))**2)/(2*widths(n)**2))
                    enddo
                    vox(i,j,k)=real(a,sp)
                enddo
            enddo
        enddo
        call padded%new([BOX,BOX,BOX]*OSMPL_PAD_FAC,1.,wthreads=.false.)
        call volume%pad_fft(padded)
        call padded%expand_cmat()
        base_exp = padded%cmat_exp
        poison = cmplx(99.,-71.,sp)
        do pass = 0,3
            !$ call omp_set_num_threads(1)
            if( pass >= 2 )then
                !$ call omp_set_num_threads(nteam)
            endif
            p%mirr_proj='yes'
            if( pass == 0 ) p%mirr_proj='no'
            do half = 1,2
                do ref = 1,NSPACE*NSTATES
                    call calc%set_ref_pft(ref,poison,half==1)
                enddo
            enddo
            do state = 1,NSTATES
                do half = 1,2
                    scale = 1.+0.2*(state-1)+0.07*(half-1)
                    padded%cmat_exp = base_exp*scale
                    call vol_pad2ref_pfts_opt(calc,padded,dirs,state,half==1)
                    do i = 1,NSPACE
                        ref = (state-1)*NSPACE+i
                        section = calc%pfts_refs(:,KS(1):KS(2),ref,half)
                        if( pass == 0 )then
                            direct(:,:,ref,half)=section
                            cycle
                        endif
                        if( pass == 1 ) serial(:,:,ref,half)=section
                        if( pass >= 2 ) call assert_true(all(section==serial(:,:,ref,half)), 'mirrored projection serial/team identity '//pgrp)
                        mat = dirs%get_mat(i)
                        interp=0; gap=0; refnorm=0
                        do k = KS(1),KS(2)
                            do ir = 1,12
                                theta = DPI*real(ir-1,dp)/12
                                q = k*(sin(theta)*real(mat(1,:),dp)-cos(theta)*real(mat(2,:),dp))
                                exact=cmplx(0._dp,0._dp,dp)
                                do n = 1,size(widths)
                                    a = amplitudes(n)*(sqrt(2*DPI)*widths(n))**3/BOX**3 &
                                        *exp(-2*DPI**2*widths(n)**2*sum(q*q)/BOX**2)
                                    phase = -2*DPI*dot_product(q,centres(:,n))/BOX
                                    exact = exact+a*cmplx(cos(phase),sin(phase),dp)
                                enddo
                                exact=exact*real(scale,dp)
                                interp = interp+abs(cmplx(direct(ir,k,ref,half),kind=dp)-exact)**2
                                gap = gap+abs(cmplx(direct(ir,k,ref,half),kind=dp)-cmplx(section(ir,k),kind=dp))**2
                                refnorm = refnorm+abs(exact)**2
                            enddo
                        enddo
                        call assert_true(interp/refnorm < 0.01_dp,'direct projection agrees with analytic Gaussian transform '//pgrp)
                        call assert_true(gap <= interp+1.e-10_dp*refnorm,'mirror difference is within interpolation error '//pgrp)
                    enddo
                enddo
            enddo
        enddo
        call calc%kill()
        call padded%kill_expanded()
        call padded%kill()
        call volume%kill()
        call dirs%kill()
        call group%kill()
        call cline%kill()
    end subroutine check_projections
end submodule simple_polarft_mirr_tester
