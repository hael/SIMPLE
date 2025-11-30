program simple_test_cavg_kpca
! include 'simple_lib.f08'
! use simple_cmdline,            only: cmdline
! use simple_builder,            only: builder
! use simple_parameters,         only: parameters
! use simple_image,              only: image
! use simple_strategy2D3D_common, only: discrete_read_imgbatch, prepimgbatch
! use simple_ctf,                 only: ctf
! use simple_fsc,                 only: plot_fsc
! implicit none
! #include "simple_local_flags.inc"
! type(builder)        :: build
! type(parameters)     :: p
! type(cmdline)        :: cline
! type(ctf)            :: tfun
! type(ctfparams)      :: ctfvars
! type(image)          :: ctfimg, rotimg, rotctfimg, cavg, rho, img
! type(image)          :: even_cavg, even_rho, odd_cavg, odd_rho
! real,    allocatable :: res(:), frc(:)
! integer, allocatable :: pinds(:)
! integer, parameter   :: MAX_CLS = 200
! character(len=:), allocatable :: last_prev_dir
! complex :: fcompl, fcompll
! real    :: shift(2), mat(2,2), dist(2), loc(2), e3, kw, tval
! integer :: logi_lims(3,2), cyc_lims(3,2), cyc_limsR(2,2), win_corner(2), phys(2), classes(MAX_CLS), cur_class, n_cls
! integer :: i,pop,h,k,l,ll,m,mm, iptcl, eo, filtsz, neven, nodd, idir
! logical :: l_ctf, l_phaseplate
! if( command_argument_count() < 3 )then
!     write(logfhandle,'(a)') 'Usage: simple_test_cavg_kpca mskdiam=xxx nthr=yy projfile=zzz.simple (stk=stk.mrc)'
!     stop
! else
!     call cline%parse_oldschool
! endif
! call cline%checkvar('projfile',1)
! call cline%checkvar('nthr',    2)
! call cline%checkvar('mskdiam', 3)
! call cline%set('mkdir',  'yes')
! call cline%set('oritype','ptcl2D')
! call cline%check
! !!!!!!! Hacky, only for test program convenience
! idir = find_next_int_dir_prefix(PATH_HERE, last_prev_dir)
! call cline%set('exec_dir', int2str(idir)//'_test_cavg_kpca')
! call cline%set('projfile',trim(PATH_PARENT)//trim(cline%get_carg('projfile')))
! if(cline%defined('stk')) call cline%set('stk',trim(PATH_PARENT)//trim(cline%get_carg('stk')))
! call simple_mkdir( filepath(string(PATH_HERE), cline%get_carg('exec_dir')))
! call simple_chdir( filepath(string(PATH_HERE), cline%get_carg('exec_dir')))
! !!!!!!!
! call build%init_params_and_build_general_tbox(cline, p, do3d=.false.)
! classes = int(build%spproj_field%get_all('class'))
! n_cls   = maxval(classes)
! do cur_class = 1, n_cls
!     call build%spproj_field%get_pinds(cur_class, 'class', pinds)
!     if( .not.(allocated(pinds)) ) cycle
!     pop = size(pinds)
!     if( pop == 0 ) cycle
!     write(logfhandle,'(A,I3)')'>>> PROCESSING CLASS',cur_class
!     l_ctf        = build%spproj%get_ctfflag(p%oritype,iptcl=pinds(1)).ne.'no'
!     l_phaseplate = .false.
!     if( l_ctf ) l_phaseplate = build%spproj%has_phaseplate(p%oritype)
!     call prepimgbatch(pop)
!     call img%new([p%boxpd,p%boxpd,1],p%smpd)
!     call rotimg%new([p%boxpd,p%boxpd,1],p%smpd)
!     call ctfimg%new([p%boxpd,p%boxpd,1],p%smpd)
!     call rotctfimg%new([p%boxpd,p%boxpd,1],p%smpd)
!     call even_cavg%new([p%boxpd,p%boxpd,1],p%smpd)
!     call odd_cavg%new([p%boxpd,p%boxpd,1],p%smpd)
!     call even_rho%new([p%boxpd,p%boxpd,1],p%smpd)
!     call odd_rho%new([p%boxpd,p%boxpd,1],p%smpd)
!     logi_lims      = img%loop_lims(2)
!     cyc_lims       = img%loop_lims(3)
!     cyc_limsR(:,1) = cyc_lims(1,:)
!     cyc_limsR(:,2) = cyc_lims(2,:)
!     call even_cavg%zero_and_flag_ft
!     call odd_cavg%zero_and_flag_ft
!     call even_rho%zero_and_flag_ft
!     call odd_rho%zero_and_flag_ft
!     neven = 0
!     nodd  = 0
!     if( cline%defined('stk') )then
!         ! reading particles from inputted stk
!         do i = 1,pop
!             call progress_gfortran(i, pop)
!             iptcl = pinds(i)
!             e3    = build%spproj_field%e3get(iptcl)
!             eo    = build%spproj_field%get_eo(iptcl)
!             if( eo ==0 )then
!                 neven = neven + 1
!             else
!                 nodd = nodd + 1
!             endif
!             call rotimg%zero_and_unflag_ft
!             call rotctfimg%zero_and_flag_ft
!             call build%imgbatch(i)%read(int2str(cur_class)//'_'//p%stk,i)
!             call build%imgbatch(i)%pad(rotimg)
!             call rotimg%fft
!             ! CTF rotation
!             if( l_ctf )then
!                 ctfvars     = build%spproj%get_ctfparams(p%oritype,iptcl)
!                 tfun        = ctf(p%smpd, ctfvars%kv, ctfvars%cs, ctfvars%fraca)
!                 if( l_phaseplate )then
!                     call tfun%ctf2img(ctfimg, ctfvars%dfx, ctfvars%dfy, ctfvars%angast, ctfvars%phshift )
!                 else
!                     call tfun%ctf2img(ctfimg, ctfvars%dfx, ctfvars%dfy, ctfvars%angast)
!                 endif
!                 call rotmat2d(-e3, mat)
!                 do h = logi_lims(1,1),logi_lims(1,2)
!                     do k = logi_lims(2,1),logi_lims(2,2)
!                         ! Rotation
!                         loc        = matmul(real([h,k]),mat)
!                         win_corner = floor(loc) ! bottom left corner
!                         dist       = loc - real(win_corner)
!                         ! Bi-linear interpolation
!                         l     = cyci_1d(cyc_limsR(:,1), win_corner(1))
!                         ll    = cyci_1d(cyc_limsR(:,1), win_corner(1)+1)
!                         m     = cyci_1d(cyc_limsR(:,2), win_corner(2))
!                         mm    = cyci_1d(cyc_limsR(:,2), win_corner(2)+1)
!                         ! l, bottom left corner
!                         phys   = ctfimg%comp_addr_phys(l,m)
!                         kw     = (1.-dist(1))*(1.-dist(2))   ! interpolation kernel weight
!                         tval   = kw * real(ctfimg%get_cmat_at(phys(1), phys(2),1))
!                         ! l, bottom right corner
!                         phys   = ctfimg%comp_addr_phys(l,mm)
!                         kw     = (1.-dist(1))*dist(2)
!                         tval   = tval   + kw * real(ctfimg%get_cmat_at(phys(1), phys(2),1))
!                         ! ll, upper left corner
!                         phys    = ctfimg%comp_addr_phys(ll,m)
!                         kw      = dist(1)*(1.-dist(2))
!                         tval    = tval  + kw * real(ctfimg%get_cmat_at(phys(1), phys(2),1))
!                         ! ll, upper right corner
!                         phys    = ctfimg%comp_addr_phys(ll,mm)
!                         kw      = dist(1)*dist(2)
!                         tval    = tval    + kw * real(ctfimg%get_cmat_at(phys(1), phys(2),1))
!                         ! update with interpolated values
!                         phys = ctfimg%comp_addr_phys(h,k)
!                         call rotctfimg%set_cmat_at(phys(1),phys(2),1, cmplx(tval*tval,0.))
!                     end do
!                 end do
!             else
!                 rotctfimg = cmplx(1.,0.)
!             endif
!             if( eo == 0 )then
!                 call even_cavg%add(rotimg)
!                 call even_rho%add(rotctfimg)
!             else
!                 call odd_cavg%add(rotimg)
!                 call odd_rho%add(rotctfimg)
!             endif
!         enddo
!     else
!         call discrete_read_imgbatch(pop, pinds(:), [1,pop] )
!         do i = 1,pop
!             call progress_gfortran(i,pop)
!             iptcl = pinds(i)
!             shift = build%spproj_field%get_2Dshift(iptcl)
!             e3    = build%spproj_field%e3get(iptcl)
!             eo    = build%spproj_field%get_eo(iptcl)
!             call img%zero_and_flag_ft
!             call rotimg%zero_and_flag_ft
!             call rotctfimg%zero_and_flag_ft
!             if( eo ==0 )then
!                 neven = neven + 1
!             else
!                 nodd = nodd + 1
!             endif
!             ! normalisation
!             call build%imgbatch(i)%norm_noise_pad_fft(build%lmsk,img)
!             ! call build%imgbatch(i)%norm_noise(build%lmsk, sdev)
!             ! shift
!             ! call build%imgbatch(i)%fft
!             call img%shift2Dserial(-shift)
!             ! ctf
!             if( l_ctf )then
!                 ctfvars     = build%spproj%get_ctfparams(p%oritype,iptcl)
!                 tfun        = ctf(p%smpd, ctfvars%kv, ctfvars%cs, ctfvars%fraca)
!                 if( l_phaseplate )then
!                     call tfun%ctf2img(ctfimg, ctfvars%dfx, ctfvars%dfy, ctfvars%angast, ctfvars%phshift )
!                 else
!                     call tfun%ctf2img(ctfimg, ctfvars%dfx, ctfvars%dfy, ctfvars%angast)
!                 endif
!                 call img%mul(ctfimg)
!             else
!                 rotctfimg = cmplx(1.,0.)
!             endif
!             ! particle & ctf rotations
!             call rotmat2d(-e3, mat)
!             do h = logi_lims(1,1),logi_lims(1,2)
!                 do k = logi_lims(2,1),logi_lims(2,2)
!                     ! Rotation
!                     loc        = matmul(real([h,k]),mat)
!                     win_corner = floor(loc) ! bottom left corner
!                     dist       = loc - real(win_corner)
!                     ! Bi-linear interpolation
!                     l     = cyci_1d(cyc_limsR(:,1), win_corner(1))
!                     ll    = cyci_1d(cyc_limsR(:,1), win_corner(1)+1)
!                     m     = cyci_1d(cyc_limsR(:,2), win_corner(2))
!                     mm    = cyci_1d(cyc_limsR(:,2), win_corner(2)+1)
!                     ! l, bottom left corner
!                     phys   = img%comp_addr_phys(l,m)
!                     kw     = (1.-dist(1))*(1.-dist(2))   ! interpolation kernel weight
!                     fcompl = kw * img%get_cmat_at(phys(1), phys(2),1)
!                     if( l_ctf ) tval   = kw * real(ctfimg%get_cmat_at(phys(1), phys(2),1))
!                     ! l, bottom right corner
!                     phys   = img%comp_addr_phys(l,mm)
!                     kw     = (1.-dist(1))*dist(2)
!                     fcompl = fcompl + kw * img%get_cmat_at(phys(1), phys(2),1)
!                     if( l_ctf ) tval   = tval   + kw * real(ctfimg%get_cmat_at(phys(1), phys(2),1))
!                     if( l < 0 ) fcompl = conjg(fcompl) ! conjugation when required!
!                     ! ll, upper left corner
!                     phys    = img%comp_addr_phys(ll,m)
!                     kw      = dist(1)*(1.-dist(2))
!                     fcompll = kw * img%get_cmat_at(phys(1), phys(2),1)
!                     if( l_ctf ) tval = tval  + kw * real(ctfimg%get_cmat_at(phys(1), phys(2),1))
!                     ! ll, upper right corner
!                     phys    = img%comp_addr_phys(ll,mm)
!                     kw      = dist(1)*dist(2)
!                     fcompll = fcompll + kw * img%get_cmat_at(phys(1), phys(2),1)
!                     if( l_ctf ) tval = tval    + kw * real(ctfimg%get_cmat_at(phys(1), phys(2),1))
!                     if( ll < 0 ) fcompll = conjg(fcompll) ! conjugation when required!
!                     ! update with interpolated values
!                     phys = img%comp_addr_phys(h,k)
!                     call rotimg%set_cmat_at(phys(1),phys(2),1, fcompl + fcompll)
!                     if( l_ctf ) call rotctfimg%set_cmat_at(phys(1),phys(2),1, cmplx(tval*tval,0.))
!                 end do
!             end do
!             if( eo == 0 )then
!                 call even_cavg%add(rotimg)
!                 call even_rho%add(rotctfimg)
!             else
!                 call odd_cavg%add(rotimg)
!                 call odd_rho%add(rotctfimg)
!             endif
!             ! write
!             call rotimg%ifft
!             call build%imgbatch(i)%set_ft(.false.)
!             call rotimg%clip(build%imgbatch(i))
!             call build%imgbatch(i)%write(int2str(cur_class) // '_aligned_ptcls_stk.mrc',i)
!             ! write ctf
!             if( l_ctf )then
!                 call rotctfimg%ifft
!                 call rotctfimg%clip(build%imgbatch(i))
!                 call build%imgbatch(i)%write(int2str(cur_class) // '_aligned_ctfs_stk.mrc',i)
!             endif
!         enddo
!     endif
!     ! deconvolutions
!     call cavg%copy(even_cavg)
!     call cavg%add(odd_cavg)
!     call rho%copy(even_rho)
!     call rho%add(odd_rho)
!     if( neven > 0 ) call even_cavg%ctf_dens_correct(even_rho)
!     if( nodd > 0 )  call odd_cavg%ctf_dens_correct(odd_rho)
!     call cavg%ctf_dens_correct(rho)
!     ! write w/o drift correction (e/o merging)
!     call even_cavg%ifft
!     call odd_cavg%ifft
!     call cavg%ifft
!     call even_cavg%clip_inplace([p%box,p%box,1])
!     call odd_cavg%clip_inplace([p%box,p%box,1])
!     call cavg%clip_inplace([p%box,p%box,1])
!     call even_cavg%write(int2str(cur_class) // '_cavg_even.mrc')
!     call odd_cavg%write(int2str(cur_class) // '_cavg_odd.mrc')
!     call cavg%write(int2str(cur_class) // '_cavg.mrc')
!     ! frc
!     if( neven>0 .and. nodd>0 )then
!         filtsz = even_cavg%get_filtsz()
!         if( allocated(frc) ) deallocate(frc)
!         allocate(frc(filtsz))
!         res = even_cavg%get_res()
!         call even_cavg%mask(p%msk,'soft',backgr=0.)
!         call odd_cavg%mask(p%msk,'soft',backgr=0.)
!         call even_cavg%fft()
!         call odd_cavg%fft()
!         call even_cavg%fsc(odd_cavg, frc)
!         call plot_fsc(filtsz, frc, res, p%smpd, int2str(cur_class)//'_frc')
!     endif
! enddo
end program simple_test_cavg_kpca
