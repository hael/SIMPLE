program simple_test_uni_filt3D
    include 'simple_lib.f08'
    use simple_cmdline,            only: cmdline
    use simple_builder,            only: builder
    use simple_parameters,         only: parameters
    use simple_commander_resolest, only: uniform_filter3D_commander
    use simple_image,              only: image
    use simple_opt_filter,         only: uni_filt3D
    implicit none
    type(parameters)                  :: p, p_opt
    type(cmdline)                     :: cline, cline_opt
    type(image)                       :: vol_clean, vol_noisy, noise_vol
    type(uniform_filter3D_commander)  :: opt_commander
    integer                           :: noise_n, noise_i, ifoo
    real                              :: ave, sdev, maxv, minv, noise_lvl
    real,    parameter                :: NOISE_MIN = 2., NOISE_MAX = 6., NOISE_DEL = 1.
    character(len=LONGSTRLEN)         :: cwd=''
    if( command_argument_count() < 3 )then
        write(logfhandle,'(a)') 'Usage: simple_test_uni_filt3D smpd=xx nthr=yy vol1=vol1.mrc'
        write(logfhandle,'(a)') 'Example: projections of https://www.rcsb.org/structure/1jyx with smpd=1. mskdiam=180'
        stop
    else
        call cline%parse_oldschool
    endif
    call cline%checkvar('smpd',    1)
    call cline%checkvar('nthr',    2)
    call cline%checkvar('vol1',    3)
    call cline%check
    call p%new(cline)
    call find_ldim_nptcls(p%vols(1), p%ldim, ifoo)
    ! setting the opt_2D filter cline
    cline_opt = cline
    call cline_opt%delete('vol')
    call cline_opt%set('frc',       'temp.bin')
    call cline_opt%set('smooth_ext', 8.)
    call cline_opt%set('prg',       'uniform_filter3D')
    call cline_opt%set('mkdir',     'yes')
    call vol_noisy%new(p%ldim, p%smpd)
    call vol_clean%new(p%ldim, p%smpd)
    call noise_vol%new(p%ldim, p%smpd)
    noise_n  = int((NOISE_MAX - NOISE_MIN)/NOISE_DEL + 1.)
    call vol_clean%read(p%vols(1))
    call vol_clean%write('vol_clean.mrc')
    call vol_clean%mask(p%msk, 'soft')
    call vol_clean%stats('foreground', ave, sdev, maxv, minv)
    ! comparing the nonuniform result with the original data
    do noise_i = 1, noise_n
        call cline_opt%set('mkdir', 'yes')
        call p_opt%new(cline_opt)
        noise_lvl = NOISE_MIN + (noise_i - 1)*NOISE_DEL
        ! adding noise
        call noise_vol%gauran(0., noise_lvl * sdev)
        call noise_vol%mask(1.5 * p%msk, 'soft')
        call vol_noisy%copy(vol_clean)
        call vol_noisy%add(noise_vol)
        call vol_noisy%write("vol_noisy_"//int2str(noise_i)//".mrc")
        call cline_opt%delete('mkdir')
        call cline_opt%set('mkdir', 'no')
        call cline_opt%delete('stk')
        call cline_opt%set('vol1', "vol_noisy_"//int2str(noise_i)//".mrc")
        call cline_opt%set('vol2', '../vol_clean.mrc')
        ! filtering by calling the commander
        call opt_commander%execute(cline_opt)
        ! spherical masking
        call vol_noisy%zero_and_unflag_ft
        call noise_vol%zero_and_unflag_ft
        call simple_getcwd(cwd)
        call simple_chdir( trim(cwd)//"/../", errmsg="")
    enddo
    call vol_noisy%kill()
    call vol_clean%kill()
    call noise_vol%kill()
    end program simple_test_uni_filt3D