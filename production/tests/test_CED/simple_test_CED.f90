program simple_test_CED
    include 'simple_lib.f08'
    use simple_cmdline,            only: cmdline
    use simple_parameters,         only: parameters
    use simple_commander_volops,   only: reproject_commander
    use simple_image,              only: image
    implicit none
    type(parameters)              :: p
    type(cmdline)                 :: cline, cline_projection
    type(image)                   :: img, ker, J11, J12, J22
    type(reproject_commander)     :: xreproject
    integer                       :: k, l, nptcls, iptcl, rc, gaussian_ext, sh
    character(len=:), allocatable :: cmd
    logical                       :: mrc_exists
    real,    parameter            :: SIGMA = 0.7, RHO = 4., C1 = 0.001, C2 = 1.
    real,    allocatable          :: grad(:,:,:), D1(:,:,:), D2(:,:,:), eig_val(:,:,:,:), eig_vec_1(:,:,:,:),&
                                    &lambda(:,:,:,:), a(:,:,:), b(:,:,:), c(:,:,:)
    complex, pointer              :: img_cmat(:,:,:)=>null(), ker_cmat(:,:,:)=>null()
    real,    pointer              :: J11_rmat(:,:,:)=>null(), J12_rmat(:,:,:)=>null(), J22_rmat(:,:,:)=>null()

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
    allocate(grad(p%ldim(1),p%ldim(2),1), D1(p%ldim(1),p%ldim(2),1),&
              &D2(p%ldim(1),p%ldim(2),1), eig_val(p%ldim(1),p%ldim(2),1,2),&
              &eig_vec_1(p%ldim(1),p%ldim(2),1,2),lambda(p%ldim(1),p%ldim(2),1,2),&
              &a(p%ldim(1),p%ldim(2),1), b(p%ldim(1),p%ldim(2),1), c(p%ldim(1),p%ldim(2),1), source=0.)
    do iptcl = 1, p%nptcls
        write(*, *) 'Particle # ', iptcl
        call img%read(p%stk, iptcl)
        ! build the Gaussian kernel with sigma
        gaussian_ext = ceiling(2*SIGMA)
        do k = -gaussian_ext, gaussian_ext
        do l = -gaussian_ext, gaussian_ext
            sh = nint(hyp(real(k),real(l)))
            call ker%set_rmat_at(p%ldim(1)/2 + k, p%ldim(2)/2 + l, 1, exp(-(sh**2/(2*SIGMA**2)))) 
        enddo
        enddo
        ! convolving img with the kernel
        call img%fft()
        call ker%fft()
        call img%get_cmat_ptr(img_cmat)
        call ker%get_cmat_ptr(ker_cmat)
        img_cmat = img_cmat*ker_cmat
        call img%ifft()
        call img%write('CED_temp_out.mrc', iptcl)
        call img%calc_gradient(grad, D1, D2)
        ! build the Gaussian kernel with rho
        gaussian_ext = ceiling(3*RHO)
        call ker%zero_and_unflag_ft()
        do k = -gaussian_ext, gaussian_ext
        do l = -gaussian_ext, gaussian_ext
            sh = nint(hyp(real(k),real(l)))
            call ker%set_rmat_at(p%ldim(1)/2 + k, p%ldim(2)/2 + l, 1, exp(-(sh**2/(2*RHO**2)))) 
        enddo
        enddo
        ! construct the structure tensor
        call J11%new(p%ldim, p%smpd)
        call J12%new(p%ldim, p%smpd)
        call J22%new(p%ldim, p%smpd)
        call J11%set_rmat(D1**2, .false.)
        call J12%set_rmat(D1*D2, .false.)
        call J22%set_rmat(D2**2, .false.)
        call J11%fft()
        call J12%fft()
        call J22%fft()
        call ker%fft()
        call ker%get_cmat_ptr(ker_cmat)
        call J11%get_cmat_ptr(img_cmat)
        img_cmat = img_cmat*ker_cmat
        call J11%ifft()
        call J11%get_rmat_ptr(J11_rmat)
        call J12%get_cmat_ptr(img_cmat)
        img_cmat = img_cmat*ker_cmat
        call J12%ifft()
        call J12%get_rmat_ptr(J12_rmat)
        call J22%get_cmat_ptr(img_cmat)
        img_cmat = img_cmat*ker_cmat
        call J22%ifft()
        call J22%get_rmat_ptr(J22_rmat)
        ! computing eigenvalues/eigenvectors of the structure tensor
        eig_val(  :,:,:,1) = (J11_rmat + J22_rmat + sqrt((J11_rmat - J22_rmat)**2 + 4*J12_rmat**2))/2
        eig_val(  :,:,:,2) = (J11_rmat + J22_rmat - sqrt((J11_rmat - J22_rmat)**2 + 4*J12_rmat**2))/2
        eig_vec_1(:,:,:,1) = 2*J12_rmat
        eig_vec_1(:,:,:,2) = J22_rmat - J11_rmat - sqrt((J11_rmat - J22_rmat)**2 + 4*J12_rmat**2)
        lambda(   :,:,:,1) = C1
        lambda(   :,:,:,2) = C1 + (1-C1)*exp(-C2/(eig_val(:,:,:,1) - eig_val(:,:,:,2))**2)
        a = lambda(:,:,:,1)*eig_vec_1(:,:,:,1)**2 + lambda(:,:,:,2)*eig_vec_1(:,:,:,2)**2
        b = (lambda(:,:,:,1) - lambda(:,:,:,2))*eig_vec_1(:,:,:,1)*eig_vec_1(:,:,:,2)
        c = lambda(:,:,:,1)*eig_vec_1(:,:,:,2)**2 + lambda(:,:,:,2)*eig_vec_1(:,:,:,1)**2
        ! solving the diffusion equations
    enddo
end program simple_test_CED