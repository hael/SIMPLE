program simple_test_cont_opt_filt
    include 'simple_lib.f08'
    use simple_cmdline,            only: cmdline
    use simple_parameters,         only: parameters
    use simple_commander_volops,   only: reproject_commander
    use simple_image,              only: image
    use simple_ced_filter,         only: ced_filter_2D
    implicit none
    type(parameters)              :: p
    type(cmdline)                 :: cline, cline_projection
    type(reproject_commander)     :: xreproject
    type(image)                   :: img, noise, ker, img_ker, img_pad
    integer                       :: nptcls, iptcl, rc, gaussian_ext, k, l
    character(len=:), allocatable :: cmd
    real,             allocatable :: pat_mat(:, :)
    logical                       :: mrc_exists
    real                          :: ave, sdev, maxv, minv, med, sh
    real, pointer                 :: img_rmat(:, :, :)
    if( command_argument_count() < 5 )then
        write(logfhandle,'(a)') 'Usage: simple_test_CED smpd=xx nthr=yy stk=stk.mrc mskdiam=zz sigma=tt'
        write(logfhandle,'(a)') 'Example: projections of https://www.rcsb.org/structure/1jyx with smpd=1. mskdiam=180 sigma=0.7'
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
        call cline%set('sigma'  , 0.7)
    else
        call cline%parse_oldschool
    endif
    call cline%checkvar('smpd',    1)
    call cline%checkvar('nthr',    2)
    call cline%checkvar('stk' ,    3)
    call cline%checkvar('mskdiam', 4)
    call cline%checkvar('sigma',   5)
    call cline%check
    call p%new(cline)
    call find_ldim_nptcls(p%stk, p%ldim, nptcls)
    p%ldim(3) = 1 ! because we operate on stacks
    call     img%new(p%ldim, p%smpd)
    call   noise%new(p%ldim, p%smpd)
    call     ker%new(p%ldim, p%smpd)
    call img_ker%new(p%ldim, p%smpd)
    ! build the Gaussian kernel with sigma
    gaussian_ext = ceiling(2*p%sigma)
    allocate(pat_mat(2*gaussian_ext, 2*gaussian_ext), source = 0.)
    do k = -gaussian_ext, gaussian_ext-1
    do l = -gaussian_ext, gaussian_ext-1
        sh = hyp(real(k),real(l))
        call ker%set_rmat_at(p%ldim(1)/2 + k, p%ldim(2)/2 + l, 1, exp(-(sh**2/(2*p%sigma**2))))
        pat_mat(gaussian_ext + k + 1, gaussian_ext + l + 1) = exp(-(sh**2/(2*p%sigma**2)))
    enddo
    enddo
    call ker%fft()
    do iptcl = 1, p%nptcls
        write(*, *) 'Particle # ', iptcl
        call img%read(p%stk, iptcl)
        call img_ker%copy_fast(img)
        ! convolving img with the kernel
        call img_ker%fft()
        img_ker = img_ker*ker
        call img_ker%ifft()
        call img_ker%write('test_img_ker.mrc', iptcl)
        ! voxel-wise convolution
        call img_ker%zero_and_unflag_ft()
        call img_pad%new([p%ldim(1)+2*gaussian_ext,p%ldim(2)+2*gaussian_ext,1],1.)
        call img%pad(img_pad)
        call img_pad%get_rmat_ptr(img_rmat)
        do l = 1, p%ldim(2)
        do k = 1, p%ldim(1)
            call img_ker%set_rmat_at(k, l, 1, sum(pat_mat*img_rmat(k:k+2*gaussian_ext,&
                                                                  &l:l+2*gaussian_ext,1)))
        enddo
        enddo
        call img_ker%write('test_img_ker_patching.mrc', iptcl)
    enddo
end program simple_test_cont_opt_filt