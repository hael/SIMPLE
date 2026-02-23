!@descr: for all fft tests
module simple_commanders_test_fft
use simple_commanders_api
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_test_corrs2weights_test
  contains
    procedure :: execute      => exec_test_corrs2weights_test
end type commander_test_corrs2weights_test

type, extends(commander_base) :: commander_test_eval_polarftcc
  contains
    procedure :: execute      => exec_test_eval_polarftcc
end type commander_test_eval_polarftcc

type, extends(commander_base) :: commander_test_ft_expanded
  contains
    procedure :: execute      => exec_test_ft_expanded
end type commander_test_ft_expanded

type, extends(commander_base) :: commander_test_gencorrs_fft
  contains
    procedure :: execute      => exec_test_gencorrs_fft
end type commander_test_gencorrs_fft

type, extends(commander_base) :: commander_test_order_corr
  contains
    procedure :: execute      => exec_test_order_corr
end type commander_test_order_corr

type, extends(commander_base) :: commander_test_phasecorr
  contains
    procedure :: execute      => exec_test_phasecorr
end type commander_test_phasecorr

type, extends(commander_base) :: commander_test_polarops
  contains
    procedure :: execute      => exec_test_polarops
end type commander_test_polarops

type, extends(commander_base) :: commander_test_rank_weights
  contains
    procedure :: execute      => exec_test_rank_weights
end type commander_test_rank_weights

type, extends(commander_base) :: commander_test_rotate_ref
  contains
    procedure :: execute      => exec_test_rotate_ref
end type commander_test_rotate_ref

contains

    subroutine exec_test_corrs2weights_test( self, cline )
        use simple_core_module_api
        class(commander_test_corrs2weights_test), intent(inout) :: self
        class(cmdline),                           intent(inout) :: cline
        real    :: corrs(12), weights(12)
        integer :: i
        corrs(1)  = -1.
        corrs(2)  = 0.0
        corrs(3)  = 0.005
        corrs(4)  = 0.1
        corrs(5)  = 0.2
        corrs(6)  = 0.3
        corrs(7)  = 0.4
        corrs(8)  = 0.5
        corrs(9)  = 0.51
        corrs(10) = 0.52
        corrs(11) = 0.53
        corrs(12) = 0.6
        weights = corrs2weights(corrs, CORRW_CRIT)
        do i=1,size(corrs)
            print *, 'corr/weight: ', corrs(i), weights(i)
        end do
        call simple_end('**** SIMPLE_TEST_CORRS2WEIGHTS_TEST_WORKFLOW NORMAL STOP ****')
    end subroutine exec_test_corrs2weights_test

subroutine exec_test_eval_polarftcc( self, cline )
    use simple_pftc_srch_api
    use simple_strategy2D3D_common, only: set_bp_range
    use simple_builder,             only: builder
    use simple_pftc_shsrch_grad,    only: pftc_shsrch_grad
    class(commander_test_eval_polarftcc), intent(inout) :: self
    class(cmdline),                       intent(inout) :: cline
    type(parameters)         :: p
    type(polarft_calc)       :: pftc
    type(builder)            :: b
    type(ori)                :: o
    real                     :: shvec(2), shift_err, ang_err, lims(2,2), cxy(3)
    real, allocatable        :: cc_fft(:)
    integer(timer_int_kind)  :: tfft
    integer                  :: loc
    type(pftc_shsrch_grad)  :: grad_shsrch_obj
    if( command_argument_count() < 4 )then
        write(logfhandle,'(a)',advance='no') 'simple_test_eval_polarftcc vol1=xx mskdiam=xx lp=xx'
        write(logfhandle,'(a)') ' smpd=xx>'
        stop
    endif
    call cline%parse_oldschool
    call cline%checkvar('vol1', 1)
    call cline%checkvar('mskdiam',  2)
    call cline%checkvar('smpd', 3)
    call cline%checkvar('lp',   4)
    call cline%set('nptcls',1.0)
    call cline%set('ctf','no')
    call cline%check
    call b%init_params_and_build_strategy3D_tbox(cline,p)
    call set_bp_range(p, b, cline)
    ang_err   = 16.
    shift_err = 8.
    call b%eulspace%get_ori(irnd_uni(p%nspace), o)
    print *,'Particle orientation:'
    call o%print_ori
    print *,'Shift= 0.0 0.0'
    print *,'---------------------'
    call pftc%new(p, p%nptcls, [1, p%nptcls], p%kfromto)
    call b%vol%read(p%vols(1))
    call b%vol%mask3D_soft(p%msk)
    call b%vol%fft()
    call b%vol%expand_cmat(p%box)
    call b%vol%fproject_polar(1, o, pftc,       iseven=.true., mask=b%l_resmsk)
    call pftc%cp_even_ref2ptcl(1,1)
    call pftc%set_eo(1, .true. )
    if( o%e3get() < 0.)then
        call o%e3set(o%e3get() - 29.5)
    else
        call o%e3set(o%e3get() + 29.5)
    endif
    call b%vol%fproject_polar(1, o, pftc,       iseven=.true., mask=b%l_resmsk)
    shvec(1) = -2.
    shvec(2) =  2.
    print *,'Ref orientation:'
    call o%print_ori
    print *,'Shift= ',shvec
    print *,'---------------------'
    call pftc%shift_ptcl(1,shvec)
    call pftc%memoize_ptcls
    !### TIMING
    allocate(cc_fft(pftc%get_nrots()))
    tfft = tic()
    call pftc%gen_objfun_vals(1, 1, [0.,0.], cc_fft)
    print *, 'time of gen_corrs (no cache): ', toc(tfft)
    loc = maxloc(cc_fft, dim=1)
    print *, pftc%get_rot(loc)
    ! searching
    lims(:,1) = -5.
    lims(:,2) =  5.
    call grad_shsrch_obj%new(lims)
    call grad_shsrch_obj%set_indices(1, 1)
    loc = 1
    tfft = tic()
    cxy = grad_shsrch_obj%minimize(irot=loc)
    print *, 'time of shift_search: ', toc(tfft)
    print *, cxy, pftc%get_rot(loc)
    call simple_end('**** SIMPLE_TEST_EVAL_POLARFTCC_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_eval_polarftcc

subroutine exec_test_ft_expanded( self, cline )
    use simple_ftexp_shsrch
    class(commander_test_ft_expanded), intent(inout) :: self
    class(cmdline),                    intent(inout) :: cline
    call seed_rnd
    call test_ftexp_shsrch
    call test_ftexp_shsrch2
    call simple_end('**** SIMPLE_TEST_FT_EXPANDED_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_ft_expanded

subroutine exec_test_gencorrs_fft( self, cline )
    use simple_pftc_srch_api
    use simple_timer
    use simple_builder, only: builder
    class(commander_test_gencorrs_fft), intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    type(parameters)        :: p
    type(polarft_calc)      :: pftc
    type(builder)           :: b
    real,    allocatable    :: cc(:), cc_fft(:)
    complex, allocatable    :: pft(:,:)
    integer                 :: iptcl, jptcl
    integer(timer_int_kind) :: tfft
    if( command_argument_count() < 3 )then
        write(logfhandle,'(a)',advance='no') 'simple_test_srch stk=<particles.mrc> msk=<mask radius(in pixels)>'
        write(logfhandle,'(a)') ' smpd=<sampling distance(in A)>'
        stop
    endif
    call cline%parse_oldschool
    call cline%checkvar('stk',  1)
    call cline%checkvar('msk',  2)
    call cline%checkvar('smpd', 3)
    call cline%check
    call p%new(cline)
    p%kfromto(1) = 2
    p%kfromto(2) = 100
    call b%build_general_tbox(p, cline)
    call pftc%new(p, p%nptcls, [1, p%nptcls], p%kfromto)
    call b%img_crop%memoize4polarize(pftc%get_pdim())
    pft = pftc%allocate_pft()
    do iptcl=1,p%nptcls
        call b%img_crop%read(p%stk, iptcl)
        call b%img_crop%fft()
        ! transfer to polar coordinates
        call b%img_crop%polarize(pft, mask=b%l_resmsk)
        call pftc%set_ref_pft(iptcl, pft, iseven=.true.)
        call b%img_crop%polarize(pft, mask=b%l_resmsk)
        call pftc%set_ptcl_pft(iptcl, pft)
    end do
    allocate(cc(pftc%get_nrots()), cc_fft(pftc%get_nrots()))

    !### TIMING

    tfft = tic()
    do iptcl=1,p%nptcls - 1
        do jptcl=iptcl + 1, p%nptcls
            call pftc%gen_objfun_vals(iptcl, jptcl, [0.,0.], cc_fft)
        end do
    end do
    print *, 'time of fft_mod: ', toc(tfft)
    call simple_end('**** SIMPLE_TEST_GENCORRS_FFT_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_gencorrs_fft

subroutine exec_test_order_corr( self, cline )
    class(commander_test_order_corr), intent(inout) :: self
    class(cmdline),                   intent(inout) :: cline
    type(oris)           :: os
    type(ori)            :: o
    integer, allocatable :: order(:)
    integer              :: i
    call seed_rnd
    call os%new(11, is_ptcl=.false.)
    call os%set_all2single('state',1.)
    do i=1,11
        call os%set(i, 'corr', ran3())
    end do
    call os%set(7, 'state', 0.)
    order = os%order()
    call os%calc_hard_weights(0.80)
    do i=1,11
        call os%get_ori(order(i), o)
        call o%print_ori()
    end do
    call o%kill
    if( count(os%get_all('w')>0.5) == 8 )then
        print *,'PASSED'
    else
        print *,'FAILED'
    endif
    call simple_end('**** SIMPLE_TEST_ORDER_CORR_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_order_corr

subroutine exec_test_phasecorr( self, cline )
    class(commander_test_phasecorr), intent(inout) :: self
    class(cmdline),                  intent(inout) :: cline
!    use mod_phasecorr
!    class(commander_test_phasecorr), intent(inout) :: self
!    class(cmdline),                  intent(inout) :: cline
!    type(t_phasecorr) :: tphasecorr
!    call tphasecorr%new
!    call tphasecorr%run
!    call tphasecorr%kill
!    call simple_end('**** SIMPLE_TEST_PHASECORR_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_phasecorr

subroutine exec_test_polarops( self, cline )
    use simple_pftc_srch_api
    use simple_builder, only: builder
    class(commander_test_polarops), intent(inout) :: self
    class(cmdline),                 intent(inout) :: cline
    complex,   allocatable :: pft(:,:)
    integer,     parameter :: N=128
    integer,     parameter :: NIMGS=200
    integer,     parameter :: NCLS=5
    type(image)            :: tmpl_img, img, cavgs(NCLS)
    type(polarft_calc)     :: pftc
    type(parameters)       :: p
    type(builder)          :: b
    real    :: ang, shift(2), shifts(2,NIMGS)
    integer :: pinds(NIMGS), i, eo, icls
    ! dummy structure
    call tmpl_img%soft_ring([N,N,1], 1., 8.)
    call tmpl_img%fft
    call tmpl_img%shift2Dserial([ 8.,-16.])
    call img%soft_ring([N,N,1], 1., 12.)
    call img%fft
    call img%shift2Dserial([ 32., 0.])
    call tmpl_img%add(img)
    call img%soft_ring([N,N,1], 1., 16.)
    call img%fft
    call img%shift2Dserial([ -16., 8.])
    call tmpl_img%add(img)
    call img%soft_ring([N,N,1], 1., 32.)
    call img%fft
    call tmpl_img%add(img)
    call tmpl_img%ifft
    call tmpl_img%write(string('template.mrc'))
    ! init of options & parameters
    call cline%set('prg',    'xxx')
    call cline%set('objfun', 'cc')
    call cline%set('smpd',   1.0)
    call cline%set('box',    N)
    call cline%set('ctf',    'no')
    call cline%set('oritype','ptcl2D')
    call cline%set('ncls',    NCLS)
    call cline%set('nptcls',  NIMGs)
    call cline%set('lp',      3.)
    call cline%set('nthr',    8)
    call cline%set('mskdiam', real(N)/2-10.)
    call cline%set('ref_type', 'polar_cavg')
    ! Calculators
    call b%init_params_and_build_strategy2D_tbox(cline, p)
    call pftc%new(p, NCLS, [1,NIMGS], p%kfromto)
    pinds = (/(i,i=1,NIMGS)/)
    call b%img_crop%memoize4polarize(pftc%get_pdim())
    pft = pftc%allocate_pft()
    do i = 1,NIMGS
        shift = 10.*[ran3(), ran3()] - 5.
        ! ang   = 360. * ran3()
        ang   = 0.
        eo    = 0
        if( .not.is_even(i) ) eo = 1
        icls  = ceiling(ran3()*4.)
        call img%copy_fast(tmpl_img)
        call img%fft
        call img%shift2Dserial(-shift)
        call img%ifft
        call img%rtsq(ang, 0.,0.)
        call img%add_gauran(2.)
        call img%write(string('rotimgs.mrc'), i)
        call img%fft
        call b%spproj_field%set_euler(i, [0.,0.,ang])
        call b%spproj_field%set_shift(i, shift)
        call b%spproj_field%set(i,'w',1.0)
        call b%spproj_field%set(i,'state',1)
        call b%spproj_field%set(i,'class', icls)
        call b%spproj_field%set(i,'eo',eo)
        shifts(:,i) = -shift
        call img%polarize(pft, mask=b%l_resmsk)
        call pftc%set_ptcl_pft(i, pft)
    enddo
    call pftc%polar_cavger_new(.false.)
    call pftc%polar_cavger_update_sums(NIMGS, pinds, b%spproj, b%esig%sigma2_noise, shifts)
    call pftc%polar_cavger_merge_eos_and_norm2D
    call pftc%polar_cavger_calc_and_write_frcs_and_eoavg(b%clsfrcs, b%spproj_field%get_update_frac(), string(FRCS_FILE), cline)
    ! write
    call pftc%polar_cavger_write(string('cavgs_even.bin'), 'even')
    call pftc%polar_cavger_write(string('cavgs_odd.bin'),  'odd')
    call pftc%polar_cavger_write(string('cavgs.bin'),      'merged')
    call pftc%polar_cavger_refs2cartesian(cavgs, 'even')
    call write_imgarr(cavgs, string('cavgs_even.mrc'))
    call pftc%polar_cavger_refs2cartesian(cavgs, 'odd')
    call write_imgarr(cavgs, string('cavgs_odd.mrc'))
    call pftc%polar_cavger_refs2cartesian(cavgs, 'merged')
    call write_imgarr(cavgs, string('cavgs_merged.mrc'))
    call pftc%polar_cavger_kill
    ! read & write again
    call pftc%polar_cavger_new(.false.)
    call pftc%polar_cavger_read(string('cavgs_even.bin'), 'even')
    call pftc%polar_cavger_read(string('cavgs_odd.bin'),  'odd')
    call pftc%polar_cavger_read(string('cavgs.bin'),      'merged')
    call pftc%polar_cavger_refs2cartesian(cavgs, 'even')
    call write_imgarr(cavgs, string('cavgs2_even.mrc'))
    call pftc%polar_cavger_refs2cartesian(cavgs, 'odd')
    call write_imgarr(cavgs, string('cavgs2_odd.mrc'))
    call pftc%polar_cavger_refs2cartesian(cavgs, 'merged')
    call write_imgarr(cavgs, string('cavgs2_merged.mrc'))
    call pftc%polar_cavger_kill
    if( allocated(pft) ) deallocate(pft)
    call simple_end('**** SIMPLE_TEST_POLAROPS_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_polarops

subroutine exec_test_rank_weights( self, cline )
    use simple_core_module_api
    use gnufor2, only: plot
    class(commander_test_rank_weights), intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    real    :: ranks(200), weights(200)
    integer :: i
    do i=1,200
        ranks(i) = real(i)
    end do
    call rank_sum_weights(200, weights)
    call plot(ranks, weights)
    call rank_inverse_weights(200, weights)
    call plot(ranks, weights)
    call rank_centroid_weights(200, weights)
    call plot(ranks, weights)
    call rank_exponent_weights(200, 10.0, weights)
    call plot(ranks, weights)
    ! do i=1,100
    !     weights(i) = 101.0 - real(i)
    ! end do
    ! do i=101,200
    !     weights(i) = real(i) - 99.5
    ! end do
    ! call plot(ranks, weights)
    ! call conv2rank_weights(200, weights, RANK_SUM_CRIT)
    ! call plot(ranks, weights)
    ! call conv2rank_weights(200, weights, RANK_CEN_CRIT)
    ! ! call plot(ranks, weights)
    ! call conv2rank_weights(200, weights, RANK_EXP_CRIT, p=2.0)
    ! ! call plot(ranks, weights)
    ! call conv2rank_weights(200, weights, RANK_INV_CRIT)
    ! ! call plot(ranks, weights)
    call simple_end('**** SIMPLE_TEST_RANK_WEIGHTS_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_rank_weights

subroutine exec_test_rotate_ref( self, cline )
    class(commander_test_rotate_ref), intent(inout) :: self
    class(cmdline),                   intent(inout) :: cline
    integer, parameter :: NP = 100, NK = 300, NR = 200
    complex :: ref(NP, NK), ref_rot(NP, NK), fast_ref_rot(NP, NK)
    real    :: a(NP, NK), b(NP, NK), start, finish, norm_all, fast_all
    integer :: i
    call random_number(a)
    call random_number(b)
    ref      = cmplx(a,b)
    norm_all = 0.
    fast_all = 0.
    do i = 1,NR
        call cpu_time(start)
        call rotate_ref(ref, i, ref_rot)
        call cpu_time(finish)
        norm_all = norm_all + (finish - start)
        call cpu_time(start)
        call fast_rotate_ref(ref, i, fast_ref_rot)
        call cpu_time(finish)
        fast_all = fast_all + (finish - start)
        if( .not.(all(ref_rot .eq. fast_ref_rot)) )then
            print *, 'FAILED'
            stop
        endif
    enddo
    print *, 'current timing = ', norm_all
    print *, '-----------'
    print *, 'improved timing = ', fast_all
    print *, '-----------'
    print *, 'PASSED'
    call simple_end('**** SIMPLE_TEST_ROTATE_REF_WORKFLOW NORMAL STOP ****')

    contains

        subroutine rotate_ref( ref_in, irot, ref_rot_out )
            complex, intent(in)  :: ref_in(NP, NK)
            integer, intent(in)  :: irot
            complex, intent(out) :: ref_rot_out(NP, NK)
            integer :: rot, jrot
            do jrot = 1,NP
                rot = jrot - (irot - 1) ! reverse rotation
                if( rot < 1 ) rot = rot + NR
                if( rot > NP )then
                    ref_rot_out(jrot,:) = conjg(ref_in(rot-NP,:))
                else
                    ref_rot_out(jrot,:) = ref_in(rot,:)
                endif
            enddo
        end subroutine rotate_ref

        subroutine rotate_ref_2( ref_in, irot, ref_rot_out )
            complex, intent(in)  :: ref_in(NP, NK)
            integer, intent(in)  :: irot
            complex, intent(out) :: ref_rot_out(NP, NK)
            integer :: rot, jrot
            do jrot = 1,irot-1
                rot = jrot - (irot - 1) + NR ! reverse rotation
                if( rot > NP )then
                    ref_rot_out(jrot,:) = conjg(ref_in(rot-NP,:))
                else
                    ref_rot_out(jrot,:) = ref_in(rot,:)
                endif
            enddo
            do jrot = irot,NP
                rot = jrot - (irot - 1) ! reverse rotation
                if( rot > NP )then
                    ref_rot_out(jrot,:) = conjg(ref_in(rot-NP,:))
                else
                    ref_rot_out(jrot,:) = ref_in(rot,:)
                endif
            enddo
        end subroutine rotate_ref_2

        subroutine rotate_ref_t( ref_in, irot, ref_rot_out )
            complex, intent(in)  :: ref_in(NK, NP)
            integer, intent(in)  :: irot
            complex, intent(out) :: ref_rot_out(NK, NP)
            integer :: rot, jrot
            do jrot = 1,NP
                rot = jrot - (irot - 1) ! reverse rotation
                if( rot < 1 ) rot = rot + NR
                if( rot > NP )then
                    ref_rot_out(:,jrot) = conjg(ref_in(:,rot-NP))
                else
                    ref_rot_out(:,jrot) = ref_in(:,rot)
                endif
            enddo
        end subroutine rotate_ref_t

        subroutine fast_rotate_ref( ref_in, irot, ref_rot_out )
            complex, intent(in)  :: ref_in(NP, NK)
            integer, intent(in)  :: irot
            complex, intent(out) :: ref_rot_out(NP, NK)
            integer :: mid
            if( irot == 1 )then
                ref_rot_out = ref_in
            elseif( irot >= 2 .and. irot <= NP )then
                mid = NP - irot + 1
                ref_rot_out(   1:irot-1,:) = conjg(ref_in(mid+1:NP, :))
                ref_rot_out(irot:NP,    :) =       ref_in(    1:mid,:)
            elseif( irot == NP + 1 )then
                ref_rot_out = conjg(ref_in)
            else
                mid = NR - irot + 1
                ref_rot_out(irot-NP:NP,       :)  = conjg(ref_in(    1:mid,:))
                ref_rot_out(      1:irot-NP-1,:) =        ref_in(mid+1:NP, :)
            endif
        end subroutine fast_rotate_ref

end subroutine exec_test_rotate_ref

end module simple_commanders_test_fft
