program simple_test_real_butterworth
    include 'simple_lib.f08'
    use simple_cmdline,            only: cmdline
    use simple_parameters,         only: parameters
    use simple_commander_volops,   only: reproject_commander
    use simple_image,              only: image
    use simple_opt_filter,         only: butterworth_filter
    implicit none
    type(parameters)              :: p
    type(cmdline)                 :: cline, cline_projection
    type(reproject_commander)     :: xreproject
    type(image)                   :: img, filt_img, img_real, img_pad
    integer                       :: nptcls, iptcl, rc, k, l, sh
    character(len=:), allocatable :: cmd
    logical                       :: mrc_exists
    real,             allocatable :: cur_filt(:), pat_mat(:, :)
    integer,          parameter   :: BW_ORDER = 8
    real,             pointer     :: img_rmat(:,:,:)
    if( command_argument_count() < 5 )then
        write(logfhandle,'(a)') 'Usage: simple_test_real_butterworth smpd=xx nthr=yy stk=stk.mrc mskdiam=zz sigma=tt'
        write(logfhandle,'(a)') 'Example: projections of https://www.rcsb.org/structure/1jyx with smpd=1. mskdiam=180 sigma=50.'
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
        call cline%set('sigma'  , 50.)
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
    call img%new(p%ldim, p%smpd)
    allocate(pat_mat(p%ldim(1), p%ldim(2)), cur_filt(p%ldim(1)), source=0.)
    call butterworth_filter(cur_filt, BW_ORDER, p%sigma)
    call filt_img%new([p%ldim(1), p%ldim(2), 1], p%smpd)
    call filt_img%set_ft(.true.)
    call filt_img%set_cmat((1., 0.))
    do l = -p%ldim(2)/2, p%ldim(2)/2-1
        do k = 0, p%ldim(1)/2
            sh = nint(hyp(real(k),real(l),0.))
            if( sh == 0 )then 
                call filt_img%mul([k,l,0], maxval(cur_filt))
            elseif( sh <= p%ldim(1)/2 )then
                call filt_img%mul([k,l,0], cur_filt(sh))
            else
                call filt_img%mul([k,l,0], 0.)
            endif
        enddo
    enddo
    call filt_img%ifft()
    call filt_img%write('filt_2D.mrc')
    call filt_img%get_rmat_sub(pat_mat)
    do iptcl = 1, p%nptcls
        write(*, *) 'Particle # ', iptcl
        call img%read(p%stk, iptcl)
        call img_real%copy(img)
        call img%fft()
        call img%apply_filter(cur_filt)
        call img%ifft()
        call img%write('img_butterworth.mrc', iptcl)
        call img_pad%new([p%ldim(1)*2,p%ldim(2)*2,1], p%smpd)
        call img%copy(img_real)
        call img%pad(img_pad)
        call img_pad%get_rmat_ptr(img_rmat)
        !$omp parallel do collapse(2) default(shared) private(k,l) schedule(static) proc_bind(close)
        do l = 1, p%ldim(2)
        do k = 1, p%ldim(1)
            call img_real%set_rmat_at(k, l, 1, sum(pat_mat*img_rmat(k:k+p%ldim(1)-1,&
                                                                   &l:l+p%ldim(2)-1,1)))
        enddo
        enddo
        !$omp end parallel do
        call img_real%write('img_real_butterworth.mrc', iptcl)
    enddo
end program