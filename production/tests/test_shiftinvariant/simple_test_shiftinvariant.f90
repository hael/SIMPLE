program simple_test_shiftinvariant
! !$ use omp_lib
! !$ use omp_lib_kinds
! include 'simple_lib.f08'
! use simple_polarft_calc,  only: polarft_calc
! use simple_cmdline,           only: cmdline
! use simple_builder,           only: builder
! use simple_image,             only: image
! use simple_parameters,        only: parameters
! use simple_polarizer,         only: polarizer
! use simple_pftc_shsrch_grad, only: pftc_shsrch_grad  ! gradient-based in-plane angle and shift search
! use simple_pftc_shsrch_fm
! use simple_strategy2D3D_common
! use simple_simulator
! use simple_ctf
! use simple_ori
! use simple_classaverager
! use simple_euclid_sigma2
! implicit none
! type(cmdline)            :: cline
! type(builder)            :: b
! type(parameters)         :: p
! type(polarft_calc)   :: pftc
! type(polarizer)          :: img_copy
! type(pftc_shsrch_grad)  :: grad_shsrch_obj
! type(pftc_shsrch_fm)    :: grad_shsrch_fm_obj
! type(ori)                :: o
! type(oris)               :: os
! type(ctfparams)          :: ctfparms
! type(euclid_sigma2)      :: eucl
! type(ctf)                :: tfun
! type(image), allocatable :: match_imgs(:), ptcl_match_imgs(:)
! type(image)              :: img, img_rot_pad
! integer(timer_int_kind)  :: t
! real(timer_int_kind)     :: rt
! logical                :: be_verbose=.false.
! logical                :: found
! real,    parameter     :: SHMAG=8.
! real,    parameter     :: SNR  =0.001
! real,    parameter     :: BFAC =20.
! integer, parameter     :: N_PTCLS = 50
! integer, allocatable   :: pinds(:)
! real, allocatable      :: sigma2_group(:,:,:), orishifts(:,:), scores(:), scores2(:), scores3(:), vals(:)
! real                   :: cxy(3), lims(2,2), lims_init(2,2), sh(2), rsh(2), e3, c, s, angerr, cenerr, aerr
! integer                :: i, irot, irot0, ne, no, iptcl, nptcls2update, ithr, nfound, iref
! if( command_argument_count() < 4 )then
!     write(logfhandle,'(a)',advance='no') 'ERROR! required arguments: '
! endif
! call cline%parse_oldschool
! call cline%checkvar('stk',      1)
! call cline%checkvar('mskdiam',  2)
! call cline%checkvar('smpd',     3)
! call cline%checkvar('lp',       4)
! call cline%check
! call cline%set('oritype','ptcl2D')
! if( .not.cline%defined('objfun') ) call cline%set('objfun', 'euclid')
! call cline%set('ml_reg', 'no')
! call cline%set('ncls',   1.)
! call cline%set('ctf',    'yes')
! call cline%set('kv',     300)
! call cline%set('cs',     2.7)
! call cline%set('fraca',  0.1)
! call cline%set('nptcls',  N_PTCLS)
! be_verbose = .false.

! ! general input
! call b%init_params_and_build_strategy2D_tbox(cline, p)
! call b%spproj%projinfo%new(1,is_ptcl=.false.)
! call b%spproj%projinfo%set(1,'projname', 'test')
! call b%spproj%os_ptcl2D%kill
! call b%spproj%os_ptcl3D%kill
! p%fromp  = 1
! p%top    = p%nptcls
! p%frcs   = trim(FRCS_FILE)
! ctfparms%smpd  = p%smpd
! ctfparms%kv    = p%kv
! ctfparms%cs    = p%cs
! ctfparms%fraca = p%fraca
! tfun = ctf(p%smpd, p%kv, p%cs, p%fraca)

! ! generate particles
! call b%img%read(p%stk, p%iptcl)
! call img%copy(b%img)
! call prepimgbatch(p%nptcls)
! call os%new(p%nptcls,is_ptcl=.true.)
! call os%rnd_ctf(p%kv, 2.7, 0.1, 2.5, 1.5, 0.001)
! call img_rot_pad%copy(b%img_pad)
! allocate(orishifts(p%nptcls,2))
! do iptcl = p%fromp,p%top
!     call os%set(iptcl,'state',1.)
!     call os%set(iptcl,'w',    1.)
!     call os%set(iptcl,'class',1.)
!     call os%set(iptcl,'proj', real(p%iptcl))
!     sh  = [gasdev( 0., SHMAG), gasdev( 0., SHMAG)]
!     call os%set(iptcl,'x',sh(1))
!     call os%set(iptcl,'y',sh(2))
!     orishifts(iptcl,:) = sh
!     e3 = ran3()*359.999
!     call os%set(iptcl,'e3',e3)
!     call os%get_ori(iptcl, o)
!     call b%img%pad(b%img_pad)
!     call simimg(b%img_pad, o, tfun, p%ctf, SNR, bfac=BFAC)
!     call b%img_pad%clip(img)
!     call img%fft
!     call img%pad(b%img_pad)
!     call b%img_pad%ifft
!     call b%img_pad%rtsq(e3,0.,0., img_rot_pad)
!     call img_rot_pad%fft
!     call img_rot_pad%shift2Dserial(p%alpha*sh)
!     call img_rot_pad%clip(b%imgbatch(iptcl))
!     call b%imgbatch(iptcl)%ifft
!     call b%imgbatch(iptcl)%write('particles.mrc',iptcl)
!     ! also works
!     ! call b%img%pad(b%img_pad)
!     ! call simimg(b%img_pad, o, tfun, p%ctf, SNR, bfac=BFAC)
!     ! call b%img_pad%rtsq(e3,0.,0., img_rot_pad)
!     ! call img_rot_pad%fft
!     ! call img_rot_pad%shift2Dserial(sh)
!     ! call img_rot_pad%ifft
!     ! call img_rot_pad%clip(b%imgbatch(iptcl))
!     ! call b%imgbatch(iptcl)%write('particles.mrc',iptcl)
! enddo
! call b%spproj%add_single_stk('particles.mrc', ctfparms, os)
! call b%spproj_field%partition_eo
! call b%spproj_field%sample4update_all([p%fromp,p%top],nptcls2update, pinds, .true.)

! ! pftc
! call pftc%new(p%nptcls, [1,p%nptcls], p%kfromto)
! call eucl%new('dummy.dat', p%box)
! call eucl%allocate_ptcls
! allocate(match_imgs(p%ncls),ptcl_match_imgs(nthr_glob))
! call pftc%reallocate_ptcls(p%nptcls, pinds)
! do ithr = 1,nthr_glob
!     call ptcl_match_imgs(ithr)%new([p%box_crop, p%box_crop, 1], p%smpd_crop, wthreads=.false.)
! enddo
! allocate(scores3(pftc%get_nrots()),scores(pftc%get_nrots()),scores2(pftc%get_pftsz()),source=0.)

! ! references
! ! call restore_read_polarize_cavgs(0)
! call read_polarize_refs(0)

! ! particles
! !$omp parallel do default(shared) private(iptcl,ithr)&
! !$omp schedule(static) proc_bind(close)
! do iptcl = 1,p%nptcls
!     ithr  = omp_get_thread_num() + 1
!     call prepimg4align(iptcl, b%imgbatch(iptcl), ptcl_match_imgs(ithr))
!     call b%img_crop_polarizer%polarize(pftc, ptcl_match_imgs(ithr), iptcl, .true., .true.)
!     call pftc%set_eo(iptcl, nint(b%spproj_field%get(iptcl,'eo'))<=0 )
! end do
! !$omp end parallel do
! call pftc%create_polar_absctfmats(b%spproj, 'ptcl2D')
! call pftc%memoize_ptcls

! ! initial sigma2
! allocate( sigma2_group(2,1,1:fdim(p%box)-1), source=0. )
! ne = 0
! no = 0
! do iptcl = p%fromp,p%top
!     call b%spproj_field%get_ori(iptcl, o)
!     ! call eucl%calc_sigma2(pftc, iptcl, o, 'class')
!     call eucl%calc_sigma2(pftc, iptcl, o, 'proj')
!     if( o%get_eo() == 0 )then
!         ne = ne+1
!         sigma2_group(1,1,:) = sigma2_group(1,1,:) + eucl%sigma2_part(:,iptcl)
!     else
!         no = no+1
!         sigma2_group(2,1,:) = sigma2_group(2,1,:) + eucl%sigma2_part(:,iptcl)
!     endif
! enddo
! sigma2_group(1,:,:) = sigma2_group(1,:,:) / real(ne)
! sigma2_group(2,:,:) = sigma2_group(2,:,:) / real(no)
! call write_groups_starfile(sigma2_star_from_iter(0), sigma2_group, 1)
! call eucl%read_groups(b%spproj_field)

! ! perturb images
! do iptcl = p%fromp,p%top
!     ithr  = omp_get_thread_num() + 1
!     call b%spproj_field%set_shift(iptcl, [0.,0.]) !!!
!     call b%imgbatch(iptcl)%read('particles.mrc',iptcl)
!     call prepimg4align(iptcl, b%imgbatch(iptcl), ptcl_match_imgs(ithr))
!     call ptcl_match_imgs(ithr)%ifft
!     call ptcl_match_imgs(ithr)%write('prepped_particles.mrc', iptcl)
!     call ptcl_match_imgs(ithr)%fft
!     call b%img_crop_polarizer%polarize(pftc, ptcl_match_imgs(ithr), iptcl, .true., .true.)
! enddo
! ! memoize again
! call pftc%create_polar_absctfmats(b%spproj, 'ptcl2D')
! call pftc%memoize_ptcls


! ! nfound = 0
! ! allocate(vals(p%fromp:p%top),source=0.)
! ! t = tic()
! ! do iptcl = p%fromp,p%top
! !     do iref = p%fromp,p%top
! !         if( trim(p%sh_inv).eq.'yes')then
! !             scores = -1.0
! !             call pftc%gen_corrs_mag_cc(iref,iptcl,scores3(1:pftc%get_pftsz()),.true.)
! !         else
! !             call pftc%gen_corrs(iref,iptcl,scores3)
! !         endif
! !         ! call pftc%gen_corrs_mag(iptcl,p%iptcl,scores(1:pftc%pftsz),kweight=.true.)
! !         vals(iref) = maxval(scores3)
! !         ! print *,iptcl,maxval(scores3), maxval(scores), orishifts(p%iptcl,:)
! !     enddo
! !     print *,iptcl,maxloc(vals,dim=1),maxval(vals)
! !     if( maxloc(vals,dim=1) == p%iptcl ) nfound = nfound+1
! ! enddo
! ! print *,nfound, toc(t)
! ! stop

! ! iptcl = 1
! ! call pftc%gen_corrs(p%iptcl,iptcl,scores3)
! ! call pftc%gen_corrs_abs(p%iptcl,iptcl,scores(1:pftc%pftsz),kweight=.false.)
! ! call pftc%gen_corrs_mag(p%iptcl,iptcl,scores2,kweight=.false.)
! ! do irot =1,pftc%get_pftsz()
! !     print *,irot,scores3(irot), pftc%calc_corr_rot_shift(p%iptcl,iptcl,[0.0,0.0],irot,.false.),scores(irot), pftc%calc_abscorr_rot(p%iptcl,iptcl,irot,.false.), scores2(irot), pftc%calc_magcorr_rot(p%iptcl,iptcl,irot,.false.)
! ! enddo
! ! ! print *, arg(orishifts(iptcl,:)),pftc%get_roind(360.-b%spproj_field%e3get(iptcl)), maxloc(scores,dim=1), maxloc(scores2,dim=1),maxloc(scores2,dim=1)+pftc%get_pftsz()
! ! stop

! ! shift search
! lims(:,1)       = -p%trs
! lims(:,2)       =  p%trs
! lims_init(:,1)  = -SHC_INPL_TRSHWDTH
! lims_init(:,2)  =  SHC_INPL_TRSHWDTH
! ! lims_init(:,1)       = -p%trs
! ! lims_init(:,2)       =  p%trs
! angerr  = 0.
! cenerr  = 0.
! do iptcl = p%fromp,p%top
!     if( trim(p%sh_inv).eq.'yes')then
!             call grad_shsrch_fm_obj%new(p%trs, 1.)
!             call pftc%gen_corrs_mag_cc(p%iptcl,iptcl,scores2,kweight=.false.)
!             irot = maxloc(scores2,dim=1)
!             call grad_shsrch_fm_obj%minimize(p%iptcl, iptcl, irot, cxy(1), cxy(2:3))
!             e3 = 360. - pftc%get_rot(irot)
!             aerr = abs(b%spproj_field%e3get(iptcl)-e3)
!             if( aerr > 180. ) aerr = 360.-aerr
!             angerr = angerr + aerr
!             cenerr = cenerr + arg(cxy(2:3)-orishifts(iptcl,:))
!             call b%spproj_field%set(iptcl, 'e3', e3)
!             call b%spproj_field%set_shift(iptcl, b%spproj_field%get_2Dshift(iptcl)+cxy(2:3))
!             print *,iptcl,'found', aerr,arg(cxy(2:3)-orishifts(iptcl,:)),cxy(2:3),orishifts(iptcl,:)
!     else
!         call grad_shsrch_obj%new(lims, lims_init=lims_init, maxits=p%maxits_sh, opt_angle=.true., coarse_init=.false.)
!         call grad_shsrch_obj%set_indices(p%iptcl, iptcl)
!         call pftc%gen_corrs(p%iptcl,iptcl,scores)
!         irot0 = maxloc(scores,dim=1)
!         cxy = grad_shsrch_obj%minimize(irot)
!         if( irot > 0 )then
!             e3 = 360. - pftc%get_rot(irot)
!             aerr = abs(b%spproj_field%e3get(iptcl)-e3)
!             if( aerr > 180. ) aerr = 360.-aerr
!             angerr = angerr + aerr
!             cenerr = cenerr + arg(cxy(2:3)-orishifts(iptcl,:))
!             call b%spproj_field%set(iptcl, 'e3', e3)
!             call b%spproj_field%set_shift(iptcl, b%spproj_field%get_2Dshift(iptcl)+cxy(2:3))
!             print *,iptcl,'found', aerr,arg(cxy(2:3)-orishifts(iptcl,:)),cxy(2:3),orishifts(iptcl,:)
!         else
!             e3 = 360. - pftc%get_rot(irot0)
!             aerr = abs(b%spproj_field%e3get(iptcl)-e3)
!             if( aerr > 180. ) aerr = 360.-aerr
!             angerr = angerr + aerr
!             print *,iptcl,aerr,orishifts(iptcl,:),cxy
!         endif
!     endif
! enddo
! angerr = angerr / real(p%nptcls)
! cenerr = cenerr / real(p%nptcls)
! print *,'angularr error: ',angerr
! print *,'shift    error: ',cenerr
! call restore_read_polarize_cavgs(1)


! contains

!     subroutine read_polarize_refs( iter )
!         integer, intent(in) :: iter
!         p%which_iter = iter
!         ! call match_imgs(1)%new([p%box_crop, p%box_crop, 1], p%smpd_crop, wthreads=.false.)
!         call b%img_crop_polarizer%init_polarizer(pftc, p%alpha)
!         do iptcl = p%fromp,p%top
!             call b%img%read(p%stk, iptcl)
!             call b%img%div(2.)
!             call b%img%mask(p%msk, 'soft')
!             call b%img%write('prepped_refs.mrc',iptcl)
!             call b%img%fft
!             call b%img_crop_polarizer%polarize(pftc, b%img, iptcl, isptcl=.false., iseven=.true.)
!             call b%img_crop_polarizer%polarize(pftc, b%img, iptcl, isptcl=.false., iseven=.false.)
!         enddo
!         call pftc%memoize_refs
!     end subroutine read_polarize_refs

!     subroutine restore_read_polarize_cavgs( iter )
!         integer, intent(in) :: iter
!         p%which_iter = iter
!         call cavger_kill()
!         call cavger_new
!         call cavger_transf_oridat( b%spproj )
!         call cavger_assemble_sums( .false. )
!         call cavger_merge_eos_and_norm
!         call cavger_calc_and_write_frcs_and_eoavg(p%frcs, p%which_iter)
!         p%refs      = trim(CAVGS_ITER_FBODY)//int2str_pad(p%which_iter,3)//p%ext
!         p%refs_even = trim(CAVGS_ITER_FBODY)//int2str_pad(p%which_iter,3)//'_even'//p%ext
!         p%refs_odd  = trim(CAVGS_ITER_FBODY)//int2str_pad(p%which_iter,3)//'_odd'//p%ext
!         call cavger_write(trim(p%refs),      'merged')
!         call cavger_write(trim(p%refs_even), 'even'  )
!         call cavger_write(trim(p%refs_odd),  'odd'   )
!         call b%clsfrcs%read(FRCS_FILE)
!         call cavger_read(trim(p%refs_even), 'even' )
!         call cavger_read(trim(p%refs_even), 'odd' )
!         call b%img_crop_polarizer%init_polarizer(pftc, p%alpha)
!         call match_imgs(1)%new([p%box_crop, p%box_crop, 1], p%smpd_crop, wthreads=.false.)
!         call prep2Dref(cavgs_even(1), match_imgs(1), 1, center=.false.)
!         call b%img_crop_polarizer%polarize(pftc, match_imgs(1), 1, isptcl=.false., iseven=.true.)
!         call prep2Dref(cavgs_odd(1), match_imgs(1), 1, center=.false.)
!         call b%img_crop_polarizer%polarize(pftc, match_imgs(1), 1, isptcl=.false., iseven=.false.)
!         call pftc%memoize_refs
!     end subroutine restore_read_polarize_cavgs

end program simple_test_shiftinvariant