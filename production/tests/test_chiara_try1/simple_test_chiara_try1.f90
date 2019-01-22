program simple_test_chiara_try1
  include 'simple_lib.f08'
  use simple_picker_chiara, only : extract_particles, discard_borders
  use simple_micops
  use simple_image,         only : image
  use simple_stackops
  use simple_math
  use simple_segmentation,  only : sobel, otsu_img
  use simple_parameters,    only : parameters
  use simple_cmdline,       only : cmdline
  use simple_sp_project,    only : sp_project
  use simple_ctf,           only : ctf
  implicit none
  type(cmdline)      :: cline
  type(parameters)   :: params
  type(image)        :: mic_shrunken, mic_draw
  type(image)        :: imgcc, mic_reduced, mic_small
  type(ctfparams)    :: ctfparms
  type(ctf)          :: tfun
  type(sp_project)   :: spproj
  character(len=STDLEN) ::  micname, tmpl
  integer            :: ldim_shrunken(3), box_shrunken, min_sz, max_sz, nmics, imic
  real               :: smpd_shrunken, part_radius
  real,    parameter :: SHRINK = 4.
  real               :: lp, thresh(1), ave, sdev, maxv, minv
  real, allocatable  :: grad(:,:,:)
  integer            :: reduced_ldim(3), box, b_fac
  logical            :: outside

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>PARTICLE PICKING WORKFLOW>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! Image processing steps: 1) Shrink and high pass filtering
  !                         2) WienerLike restoration (--->Total Variation??)
  !                         3) LowPass Filtering
  !                         4) Border discarding
  !                         5) Standardization by histogram stretching
  !                         6) Binarization (sobel|otsu|bin)
  !                         7) Connected components (cc) identification
  !                         8) cc size filtering
  !                         9) Particles extraction

  if( command_argument_count() < 3 )then
      write(logfhandle,'(a)',advance='no') 'simple_test_chiara_try projfile=<project file> [part_radius = <radius of the particle (# pixels)] [detector= <binarisation method> (sobel|otsu|bin)]'
      stop
  endif
  call cline%parse_oldschool
  call cline%checkvar('projfile',    1)
  call cline%checkvar('part_radius', 2)
  call cline%checkvar('detector',    3)
  !Set defaults
  call params%new(cline)  !<read cline parameters
  call spproj%read(params%projfile)
  nmics = spproj%get_nintgs()
  do imic = 1,min(params%top,nmics)
      micname  = spproj%os_mic%get_static(imic,'intg')          !all the path to the name of the file
      tmpl     = get_fbody(basename(micname),trim(params%ext))  !just the name of the file with no .mrc
      write(logfhandle,*)'micname = ', trim(micname)
      write(logfhandle,*)'tmpl = ', trim(tmpl)   !TO FIX
      ctfparms = spproj%os_mic%get_ctfvars(imic)
      call progress(imic,nmics)

      ! 0) Reading and saving original micrograph
      call read_micrograph(micfname = micname, smpd = ctfparms%smpd)

      !>>>>>>>>>>>>>>>>>>PRE PROCESSING ROUTINES>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! 1) Shrink and high pass filtering
      part_radius = params%part_radius/SHRINK
      call shrink_micrograph(SHRINK, ldim_shrunken, smpd_shrunken)
      call set_box(int(SHRINK*( 4*part_radius+part_radius )), box_shrunken) !here HP filt
      call mic_shrunken%new(ldim_shrunken, smpd_shrunken)
      call mic_shrunken%read('shrunken_hpassfiltered.mrc')

      ! 2) Wiener-Like restoration
      b_fac = -30
      call mic_shrunken%fft
      tfun = ctf(smpd_shrunken,ctfparms%kv,ctfparms%cs,ctfparms%fraca)
      call tfun%wienerlike_restoration(mic_shrunken, ctfparms, b_fac)
      call mic_shrunken%ifft
      call mic_shrunken%write('wienerlike_restored.mrc')

      ! 3) Non local mean filtering
      ! call mic_shrunken%NLmean()
      ! call mic_shrunken%write(tmpl//'NLmean_filtered.mrc')

      ! 4) Low pass filtering
      lp = (params%part_radius)/6 !lp ~ 1/6 rad of the particle
      write(logfhandle,*)'lp = ', lp
      call mic_shrunken%bp(0.,lp)
      call mic_shrunken%ifft()
      call mic_shrunken%write(tmpl//'bp_filtered.mrc')

      ! 5) Border discarding
      ! In order to have consistent statistics I have to get rid of the borders.
      ! I esimate they are about 4% of the dimension of the image
      call discard_borders(mic_shrunken, mic_reduced, reduced_ldim)

      ! 6) Standardization by Histogram stretching
      call mic_reduced%hist_stretching(mic_small)
      call mic_small%write('HistogramStretched.mrc')
      !call mic_small%copy(mic_reduced)  !to get rid of the histogram stretching
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!just to draw on that version!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      box = minval(ldim_shrunken(:2))- 2*reduced_ldim(1)
      call mic_draw%new([box,box,1], smpd_shrunken) !work on a copy not to modify the original mic
      call mic_shrunken%window_slim(reduced_ldim(:2),box, mic_draw, outside)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>BINARIZATION>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      allocate(grad(reduced_ldim(1),reduced_ldim(2),reduced_ldim(3)), source = 0.)
      ! 7) Binarization
      if(params%detector .eq. 'sobel') then
        call mic_small%calc_gradient(grad)
        call mic_reduced%set_rmat(grad) !now its is the gradient
        call mic_reduced%stats( ave, sdev, maxv, minv )
        thresh(1) = ave+.7*sdev
        write(logfhandle,*)'gradient threshold = ', ave+.7*sdev
        call sobel(mic_small,thresh)
      else if (params%detector .eq. 'bin') then
        call mic_small%stats( ave, sdev, maxv, minv )
        call mic_small%bin(ave+.7*sdev)
        write(logfhandle,*)'mic threshold = ', ave+.7*sdev
      else if (params%detector .eq. 'otsu') then
        call otsu_img(mic_small)
      else
        write(logfhandle,*) 'Invalid detector parameter! simple_test_chiara_try1'
        stop
      endif
      call mic_small%write(tmpl//'Bin.mrc')
      if(params%detector .eq. 'sobel') then
         call mic_small%real_space_filter(1,'median') !median filtering allows easy calculation of cc. winsz=5??
      else
          call mic_small%real_space_filter(3,'median')
      endif
      call mic_small%write('BinMedian.mrc')

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>CONNECTED COMPONENTS METHODS>>>>>>>>>>>>>>>>>>>>>>
      ! 8) Connected components (cc) identification
      call imgcc%new(mic_small%get_ldim(), smpd_shrunken)
      call mic_small%find_connected_comps(imgcc)

      ! 9) cc filtering (to fix according to Adiga s paper?)
      min_sz = 3*int(part_radius)
      max_sz = 6*int(part_radius*part_radius)   !~ 2times area of a circle
      call imgcc%elim_cc([min_sz,max_sz])
      call imgcc%write(tmpl//'ConnectedComponentsElimin.mrc')

      ! 10) Particles extraction
      call extract_particles(mic_draw, imgcc, int(part_radius))
  enddo
  open(unit = 17, file = "PickerInfo.txt")
  write(unit = 17, fmt = '(a)') '>>>>>>>>>>>>>>>>>>>>PARTICLE PICKING>>>>>>>>>>>>>>>>>>'
  write(unit = 17, fmt = '(a)') ''
  write(unit = 17, fmt = "(a,f0.0)")  'Mic Shrunken by factor ', SHRINK
  write(unit = 17, fmt = "(a,i0,tr1,i0,tr1,i0)") 'Dim after shrink ', ldim_shrunken
  write(unit = 17, fmt = "(a,f0.0)")  'Smpd after shrink ', smpd_shrunken
  write(unit = 17, fmt = "(a,i0)")  'Hp box ', int(SHRINK*( 4*part_radius+part_radius ))
  write(unit = 17, fmt = "(a,f0.0)")  'Lp  ', lp
  write(unit = 17, fmt = "(a,i0,tr1, i0)")  'Connected Components size filtering  ', min_sz, max_sz
  write(unit = 17, fmt = '(a)') ''
  write(unit = 17, fmt = "(a)")  'SELECTED PARAMETERS '
  write(unit = 17, fmt = '(a)') ''
  write(unit = 17, fmt = "(a,tr1,f0.0)")  'smpd ', params%smpd
  write(unit = 17, fmt = "(a,tr1,f0.0)")  'part_radius  ', params%part_radius
  write(unit = 17, fmt = "(a,a)")  'detector  ', params%detector
  close(17, status = "keep")
end program simple_test_chiara_try1
