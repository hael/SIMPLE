program simple_test_eval_polarftcc
    include 'simple_lib.f08'
    use simple_polarft_corrcalc, only: polarft_corrcalc
    use simple_cmdline,           only: cmdline
    use simple_builder,           only: builder
    use simple_parameters,        only: parameters
    use simple_timer
    implicit none
    type(parameters)        :: p
    type(polarft_corrcalc)  :: pftcc
    type(cmdline)           :: cline
    type(builder)           :: b
    real, allocatable       :: cc(:), cc_fft(:)
    integer                 :: iptcl, jptcl, irot, loc_cc(1), loc_cc_fft(1), nerrors, cnt
    integer(timer_int_kind) :: torig, tfft
    real                    :: err, erravg, errmax
    if( command_argument_count() < 3 )then
        write(logfhandle,'(a)',advance='no') 'simple_test_eval_polarftcc stk=<particles.mrc> mskdiam=<mask radius(in pixels)>'
        write(logfhandle,'(a)') ' smpd=<sampling distance(in A)>'
        stop
    endif
    call cline%parse_oldschool
    call cline%checkvar('stk',  1)
    call cline%checkvar('mskdiam',  2)
    call cline%checkvar('smpd', 3)
    call cline%check
    call p%new(cline)
    call b%build_general_tbox(p, cline)
    call pftcc%new(p%nptcls, [1, p%nptcls], .false.)
    call b%img_match%init_polarizer(pftcc, p%alpha)
    do iptcl=1,p%nptcls
        call b%img_match%read(p%stk, iptcl)
        call b%img_match%fft()
        ! transfer to polar coordinates
        call b%img_match%polarize(pftcc, iptcl, isptcl=.false., iseven=.true.)
        call b%img_match%polarize(pftcc, iptcl, isptcl=.true.,  iseven=.true.)
    end do
    call pftcc%memoize_ffts()
    allocate(cc(pftcc%get_nrots()), cc_fft(pftcc%get_nrots()))

    !### TIMING

    p%kfromto(1) = 2
    p%kfromto(2) = 100
    tfft = tic()
    do iptcl=1,p%nptcls
        call pftcc%gencorrs(iptcl, iptcl, cc_fft)
    end do
    print *, 'time of fft_mod: ', toc(tfft)
end program simple_test_eval_polarftcc
