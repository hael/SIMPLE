program simple_test_uni_filt2D
    include 'simple_lib.f08'
    use simple_cmdline,            only: cmdline
    use simple_builder,            only: builder
    use simple_parameters,         only: parameters
    use simple_commander_resolest, only: uniform_filter2D_commander
    use simple_image,              only: image
    use simple_opt_filter,         only: uni_filt2D
    implicit none
    type(parameters)                  :: p, p_opt
    type(cmdline)                     :: cline, cline_opt
    type(image)                       :: img_clean, img_noisy, noise_img
    type(uniform_filter2D_commander) :: opt_commander
    integer                           :: nptcls, iptcl, noise_n, noise_i
    real                              :: ave, sdev, maxv, minv, noise_lvl
    real,    parameter                :: LP_LOWRES_PHASE = 7., NOISE_MIN = .3, NOISE_MAX = 1., NOISE_DEL = 0.1
    integer, parameter                :: NSEARCH = 100, SMOOTH_EXT = 8
    character(len=LONGSTRLEN)         :: cwd=''
    if( command_argument_count() < 3 )then
        write(logfhandle,'(a)') 'Usage: simple_test_uni_filt2D smpd=xx nthr=yy stk=stk.mrc'
        write(logfhandle,'(a)') 'Example: projections of https://www.rcsb.org/structure/1jyx with smpd=1. mskdiam=180'
        stop
    else
        call cline%parse_oldschool
    endif
    call cline%checkvar('smpd',    1)
    call cline%checkvar('nthr',    2)
    call cline%checkvar('stk' ,    3)
    call cline%check
    call p%new(cline)
    call find_ldim_nptcls(p%stk, p%ldim, nptcls)
    ! setting the opt_2D filter cline
    cline_opt = cline
    call cline_opt%delete('stk')
    call cline_opt%set('frcs',       'temp.bin')
    call cline_opt%set('smooth_ext', 8.)
    call cline_opt%set('prg',        'uniform_filter2D')
    call cline_opt%set('mkdir',      'yes')
    p%ldim(3) = 1 ! because we operate on stacks
    call img_noisy%new(p%ldim, p%smpd)
    call img_clean%new(p%ldim, p%smpd)
    call noise_img%new(p%ldim, p%smpd)
    noise_n  = int((NOISE_MAX - NOISE_MIN)/NOISE_DEL + 1.)
    do iptcl = 1, min(1, p%nptcls)
        write(*, *) 'Particle # ', iptcl
        call img_clean%read(p%stk, iptcl)
        call img_clean%write('stk_img_clean.mrc', iptcl)
        call img_clean%mask(p%msk, 'soft')
        call img_clean%stats('foreground', ave, sdev, maxv, minv)
        ! comparing the nonuniform result with the original data
        do noise_i = 1, noise_n
            call cline_opt%set('mkdir', 'yes')
            call p_opt%new(cline_opt)
            noise_lvl = NOISE_MIN + (noise_i - 1)*NOISE_DEL
            ! adding noise
            call noise_img%gauran(0., noise_lvl * sdev)
            call noise_img%mask(1.5 * p%msk, 'soft')
            call img_noisy%copy(img_clean)
            call img_noisy%add(noise_img)
            call img_noisy%write("stk_img_noisy_"//int2str(noise_i)//".mrc", iptcl)
            call cline_opt%delete('mkdir')
            call cline_opt%set('mkdir', 'no')
            call cline_opt%delete('stk')
            call cline_opt%set('stk',  "stk_img_noisy_"//int2str(noise_i)//".mrc")
            call cline_opt%set('stk2', '../stk_img_clean.mrc')
            ! filtering by calling the commander
            call opt_commander%execute(cline_opt)
            ! spherical masking
            call img_noisy%zero_and_unflag_ft
            call noise_img%zero_and_unflag_ft
            call simple_getcwd(cwd)
            call simple_chdir( trim(cwd)//"/../", errmsg="")
        enddo
    enddo
    call img_noisy%kill()
    call img_clean%kill()
    call noise_img%kill()
    end program simple_test_uni_filt2D