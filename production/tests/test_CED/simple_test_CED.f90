program simple_test_CED
    include 'simple_lib.f08'
    use simple_cmdline,            only: cmdline
    use simple_parameters,         only: parameters
    use simple_commander_volops,   only: reproject_commander
    use simple_image,              only: image
    implicit none
    type(parameters)              :: p
    type(cmdline)                 :: cline, cline_projection
    type(image)                   :: img, ker
    type(reproject_commander)     :: xreproject
    integer                       :: k, l, nptcls, iptcl, rc, gaussian_ext, sh
    character(len=:), allocatable :: cmd
    logical                       :: mrc_exists
    real,    parameter            :: SIGMA = 0.7
    real,    allocatable          :: grad(:,:,:), Dc(:,:,:), Dr(:,:,:)
    complex, pointer              :: img_cmat(:,:,:)=>null(), ker_cmat(:,:,:)=>null()

    character(len=15), parameter, dimension(3) :: FIL_ARR = [character(len=15) :: "tv", "butterworth", "lp"]
    if( command_argument_count() < 4 )then
        write(logfhandle,'(a)') 'Usage: simple_test_CED smpd=xx nthr=yy stk=stk.mrc, mskdiam=zz'
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
    call cline%checkvar('mskdiam', 4)
    call cline%check
    call p%new(cline)
    call find_ldim_nptcls(p%stk, p%ldim, nptcls)
    p%ldim(3) = 1 ! because we operate on stacks
    call img%new(p%ldim, p%smpd)
    call ker%new(p%ldim, p%smpd)
    allocate(grad(p%ldim(1),p%ldim(2),1), Dc(p%ldim(1),p%ldim(2),1),&
            &Dr(p%ldim(1),p%ldim(2),1), source=0.)
    do iptcl = 1, p%nptcls
        write(*, *) 'Particle # ', iptcl
        call img%read(p%stk, iptcl)
        ! build the Gaussian kernel
        gaussian_ext = ceiling(2*SIGMA)
        do k = -gaussian_ext, gaussian_ext
        do l = -gaussian_ext, gaussian_ext
            sh = nint(hyp(real(k),real(l)))
            call ker%set_rmat_at(p%ldim(1)/2 + k, p%ldim(2)/2 + l, 1, exp(-(sh**2/(2*SIGMA**2)))) 
        enddo
        enddo
        call img%fft()
        call ker%fft()
        call img%get_cmat_ptr(img_cmat)
        call ker%get_cmat_ptr(ker_cmat)
        img_cmat = img_cmat*ker_cmat
        call img%ifft()
        call img%write('CED_temp_out.mrc', iptcl)
        call img%calc_gradient(grad, Dc, Dr)
    enddo
end program simple_test_CED