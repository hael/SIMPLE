program simple_test_eval_polarftcc
    include 'simple_lib.f08'
    use simple_polarft_corrcalc,    only: polarft_corrcalc
    use simple_cmdline,             only: cmdline
    use simple_builder,             only: builder, build_glob
    use simple_parameters,          only: parameters
    use simple_pftcc_shsrch_grad,   only: pftcc_shsrch_grad
    use simple_strategy2D3D_common, only: set_bp_range
    use simple_ori,                 only: ori
    implicit none
    type(parameters)         :: p
    type(polarft_corrcalc)   :: pftcc
    type(cmdline)            :: cline
    type(builder)            :: b
    type(ori)                :: o
    real                     :: shvec(2), shift_err, ang_err, lims(2,2), cxy(3)
    real, allocatable        :: cc_fft(:)
    integer(timer_int_kind)  :: tfft
    integer                  :: loc
    type(pftcc_shsrch_grad)  :: grad_shsrch_obj
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
    call set_bp_range(cline)
    
    ang_err   = 16.
    shift_err = 8.
    call b%eulspace%get_ori(irnd_uni(p%nspace), o)
    print *,'Particle orientation:'
    call o%print_ori
    print *,'Shift= 0.0 0.0'
    print *,'---------------------'
    
    call pftcc%new(p%nptcls, [1, p%nptcls], .false.)
    call b%vol%read(p%vols(1))
    call b%vol%mask(p%msk,'soft')
    if( p%gridding.eq.'yes' ) call b%vol%div_w_instrfun(p%interpfun, alpha=p%alpha)
    call b%vol%fft()
    call b%vol%expand_cmat(p%alpha,norm4proj=.true.)
    call b%vol%fproject_polar(1, o, pftcc, iseven=.true., mask=b%l_resmsk)
    call pftcc%cp_even_ref2ptcl(1,1)
    call pftcc%set_eo(1, .true. )
        
    if( o%e3get() < 0.)then
        call o%e3set(o%e3get() - 29.5)
    else
        call o%e3set(o%e3get() + 29.5)
    endif
    call b%vol%fproject_polar(1, o, pftcc, iseven=.true., mask=b%l_resmsk)
    shvec(1) = -2.
    shvec(2) =  2.
    print *,'Ref orientation:'
    call o%print_ori
    print *,'Shift= ',shvec
    print *,'---------------------'
    
    call pftcc%shift_ptcl(1,shvec)
    call pftcc%memoize_ffts
    
    lims(:,1) = -5.
    lims(:,2) =  5.
    call grad_shsrch_obj%new(lims)
    call grad_shsrch_obj%set_indices(1, 1)
    loc = 1
    tfft = tic()
    cxy = grad_shsrch_obj%minimize(irot=loc)
    print *, 'time of shift_search: ', toc(tfft)
    print *, cxy, pftcc%get_rot(loc)

    !### TIMING
    allocate(cc_fft(pftcc%get_nrots()))
    tfft = tic()
    call pftcc%gencorrs(1, 1, cc_fft)
    print *, 'time of gencorrs: ', toc(tfft)
    loc = maxloc(cc_fft, dim=1)
    print *, pftcc%get_rot(loc)
end program simple_test_eval_polarftcc
