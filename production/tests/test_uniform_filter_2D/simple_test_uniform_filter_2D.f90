program simple_test_uni_filt2D
    include 'simple_lib.f08'
    use simple_cmdline,            only: cmdline
    use simple_builder,            only: builder
    use simple_parameters,         only: parameters
    use simple_commander_resolest, only: uniform_filter2D_commander
    use simple_commander_volops,   only: reproject_commander
    use simple_image,              only: image
    use simple_opt_filter,         only: uni_filt2D
    implicit none
    type(parameters)                 :: p
    type(cmdline)                    :: cline, cline_opt, cline_projection
    type(image)                      :: img_clean, img_noisy, noise_img
    type(uniform_filter2D_commander) :: opt_commander
    type(reproject_commander)        :: xreproject
    character(len=:), allocatable    :: cmd
    integer                          :: nptcls, iptcl, noise_n, noise_i, rc
    real                             :: ave, sdev, maxv, minv, noise_lvl
    real,    parameter               :: LP_LOWRES_PHASE = 7., NOISE_MIN = .3, NOISE_MAX = 1., NOISE_DEL = 0.1
    integer, parameter               :: NSEARCH = 100, SMOOTH_EXT = 8
    character(len=LONGSTRLEN)        :: cwd=''
    logical                          :: mrc_exists
    if( command_argument_count() < 3 )then
        write(logfhandle,'(a)') 'ERROR! Usage: simple_test_uniform_filt2D smpd=xx nthr=yy stk=stk.mrc'
        write(logfhandle,'(a)') 'Example: projections of https://www.rcsb.org/structure/1jyx with smpd=1. mskdiam=180'
        write(logfhandle,'(a)') 'DEFAULT TEST (example above) is running now...'
        inquire(file="1JYX.mrc", exist=mrc_exists)
        if( .not. mrc_exists )then
            write(*, *) 'Downloading the example dataset...'
            cmd = 'curl -s -o 1JYX.pdb https://files.rcsb.org/download/1JYX.pdb'
            call execute_command_line(cmd, exitstat=rc)
            write(*, *) 'Converting .pdb to .mrc...'
            cmd = 'e2pdb2mrc.py 1JYX.pdb 1JYX.mrc'
            call execute_command_line(cmd, exitstat=rc)
            cmd = 'rm 1JYX.pdb'
            call execute_command_line(cmd, exitstat=rc)
            write(*, *) 'Projecting 1JYX.mrc...'
            call cline_projection%set('vol1'      , '1JYX.mrc')
            call cline_projection%set('smpd'      , 1.)
            call cline_projection%set('pgrp'      , 'c1')
            call cline_projection%set('mskdiam'   , 180.)
            call cline_projection%set('nspace'    , 6.)
            call cline_projection%set('nthr'      , 16.)
            call xreproject%execute(cline_projection)
        endif
        call cline%set('smpd'   , 1.)
        call cline%set('nthr'   , 16.)
        call cline%set('stk'    , 'reprojs.mrcs')
        call cline%set('mskdiam', 180.)
    else
        call cline%parse_oldschool
    endif
    call cline%checkvar('smpd',    1)
    call cline%checkvar('nthr',    2)
    call cline%checkvar('stk' ,    3)
    call cline%check
    call p%new(cline)
    call find_ldim_nptcls(p%stk, p%ldim, nptcls)
    p%ldim(3) = 1 ! because we operate on stacks
    call img_noisy%new(p%ldim, p%smpd)
    call img_clean%new(p%ldim, p%smpd)
    call noise_img%new(p%ldim, p%smpd)
    noise_n  = int((NOISE_MAX - NOISE_MIN)/NOISE_DEL + 1.)
    do iptcl = 1, p%nptcls
        write(*, *) 'Particle # ', iptcl
        call img_clean%read(p%stk, iptcl)
        call img_clean%write('stk_img_clean.mrc', iptcl)
        call img_clean%mask(p%msk, 'soft')
        call img_clean%stats('foreground', ave, sdev, maxv, minv)
        ! comparing the nonuniform result with the original data
        do noise_i = 1, noise_n
            noise_lvl = NOISE_MIN + (noise_i - 1)*NOISE_DEL
            ! adding noise
            call noise_img%gauran(0., noise_lvl * sdev)
            call noise_img%mask(1.5 * p%msk, 'soft')
            call img_noisy%copy(img_clean)
            call img_noisy%add(noise_img)
            call img_noisy%write("stk_img_noisy_"//int2str(noise_i)//".mrc", iptcl)
            ! spherical masking
            call img_noisy%zero_and_unflag_ft
            call noise_img%zero_and_unflag_ft
        enddo
    enddo
    call img_noisy%kill()
    call img_clean%kill()
    call noise_img%kill()
    ! setting the opt_2D filter cline
    cline_opt = cline
    call cline_opt%set('frcs',       'temp.bin')
    call cline_opt%set('smooth_ext', 8.)
    call cline_opt%set('mkdir',      'yes')
    do noise_i = 1, noise_n
        call cline_opt%delete('stk')
        call cline_opt%set('stk',  'stk_img_noisy_'//int2str(noise_i)//'.mrc')
        call cline_opt%set('stk2', 'stk_img_clean.mrc')
        ! filtering by calling the commander
        call opt_commander%execute(cline_opt)
    enddo
end program simple_test_uni_filt2D
