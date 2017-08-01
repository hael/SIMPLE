!==Program simple_stackops
!
! <stackops/begin> is a program that provides standard single-particle image processing routines that are applied to MRC or SPIDER stacks.
! <comment/begin> You can do many things with \prgname{simple\_stackops}. Inputting two stacks of the same size results in the calculation of
! the joint Fourier Ring  Correlation (FRC) between the images. Inputting no stacks, but setting \texttt{nptcls}, results in production of
! \texttt{nptcls} pure noise images, unless \texttt{ctf=yes}, then CTF images are produced. Filtering is controlled by the \texttt{hp}
! and \texttt{lp} arguments.  If you want to extract a particular state, give an alignment document (\texttt{oritab}) and set \texttt{state}
! to the state that you want to extract. If you want to select the fraction of best particles (according to the goal function), input an
! alignment doc (\texttt{oritab}) and set \texttt{frac}. You can combine the \texttt{state} and \texttt{frac} options. If you
! want to apply noise to images, give the desired signal-to-noise ratio via \texttt{snr}. If you want to mask your images with a spherical
! mask with a soft falloff, set \texttt{msk} to the radius in pixels. If you want to calculate the autocorrelation function of your images set \texttt{acf=yes}. If you
! want to randomise the phases of the Fourier transforms of your images, set \texttt{phrand=yes} and \texttt{lp} to the desired low-pass
! limit. If you want to extract a contiguous subset of particle images from the stack, set \texttt{fromp} and \texttt{top}. If you want
! to fish out a number of particle images from your stack at random, set \texttt{nran} to some nonzero integer number less than \texttt{nptcls}.
! If you want to resize you images, set the desired box to \texttt{newbox} or use the \texttt{scale} option. It is often
! convenient to use \texttt{scale} in combination with \texttt{clip} to resize images.  If you want to normalise your images, set
! \texttt{norm=yes}. \texttt{hfun} controls the normalisation function. With \texttt{avg=yes} the global average of the inputted stack
! is calculated. With \texttt{ctf=flip} the contrast inversions due to the CTF are corrected by the
! infamous (but effective) phase-flipping heuristic. This requires additional input of CTF-related parameters (\texttt{kv}, \texttt{fraca}
! and \texttt{cs}) in addition to the defocus and astigmatism angle values, communicated either via \texttt{oritab} or via \texttt{deftab}. 
! Even if you do initially
! phase-flip the images, which you should do for initial model production with PRIME, you can turn on the Wiener restoration later anyway,
! to accomplish correct weighting of information around the CTF zeroes and maximal noise reduction. \texttt{ft2img=yes} produces images of
! the square power spectrum of the images in \texttt{stk}. If you define \texttt{frameavg} to some integer number larger than one averages 
! with chunk sizes of \texttt{frameavg} are produced, which may be useful for analysis of dose-fractionated image series. \texttt{clip} can 
! be used to re-window or pad the images to a different box size. When \texttt{compare=yes}, the two inputted stacks are Fourier ring correlated.
! \texttt{neg} inverts the contrast of the images. \texttt{ctfsq} applies the squared CTF to
! the inputted images. \texttt{inner} is for applying an inner mask with fall-off width \texttt{width}. Finally, \texttt{append} is for
! appending stack \texttt{stk2} with stack \texttt{stk}, so that the \texttt{stk2} images occur last in the series and the stk name is 
! preserved.<comment/end> <stackops/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2015
!
program simple_stackops
use simple_jiffys,     ! singleton
use simple_cmdline,    only: cmdline
use simple_build,      only: build
use simple_params,     only: params
use simple_image,      only: image
use simple_oris,       only: oris
use simple_ran_tabu,   only: ran_tabu
use simple_math,       only: get_resolution
use simple_imgfile,    only: imgfile
use simple_ori,        only: ori
use simple_stat,       only: moment
use simple_map_reduce, only: splitstk_in_parts
use simple_procimgfile
use simple_timing
implicit none
type(params)                  :: p
type(build)                   :: b
type(cmdline)                 :: cline
type(image)                   :: img, img2, mskimg
type(oris)                    :: o_here
type(ori)                     :: o
type(ran_tabu)                :: rt
integer, allocatable          :: pinds(:)
real, allocatable             :: corrs(:), res(:), corrs_sum(:)
integer                       :: i, j, nincl, s, cnt, alloc_stat, ldim_scaled(3), iptcl
integer                       :: np1, np2, npix, funit_stacks, ldim(3), nstacks, npatterns
integer                       :: funit_patterns, nimgs, ndefs, ntot, lfoo(3), ifoo
real                          :: fsc0143, fsc05, ave, sdev, var, corr
character(len=STDLEN)         :: stackname, patternname
character(len=:), allocatable :: fname
character(len=1)              :: imgformat
logical                       :: err, rawproc
logical, parameter            :: debug=.false.
!Start of the execution commands
call timestamp()
call start_Alltimers_cpu()
!parameter
if( command_argument_count() < 1 )then
    write(*,'(a)', advance='no') 'SIMPLE_STACKOPS [stk=<stack.ext>] [stk2=<stack2.ext>] [nptcls=<nr of imgs>]'
    write(*,'(a)', advance='no') ' [smpd=<sampling distance(in A)>] [outstk=<outstk.ext>] [split=<nr of partitions to'
    write(*,'(a)', advance='no') ' split the stack into>] [oritab=<SIMPLE alignment doc>] [hp=<high-pass limit(in A)>]'
    write(*,'(a)', advance='no') '(in A)>]  [mul=<shift multiplication factor{1}>] [trs=<origin'
    write(*,'(a)', advance='no') ' shift halfwidth(in pixels){0}] [lp=<low-pass limit(in A){20}>]'
    write(*,'(a)', advance='no') ' [state=<state to extract>] [frac=<fraction of ptcls to extract{1}>]'
    write(*,'(a)', advance='no') ' [class=<cluster2extract>]'
    write(*,'(a)', advance='no') ' [snr=<signal2noise ratio>] [msk=<mask radius(in pixels){box/2}>] [vis=<yes|no>]'
    write(*,'(a)', advance='no') ' [acf=<yes|no{no}>] [phrand=<yes|no{no}>] [fromp=<start'
    write(*,'(a)', advance='no') ' ptcl>] [top=<stop ptcl>] [nran=<number of random images to select>] [newbox=<scaled box>]'
    write(*,'(a)', advance='no') ' [scale=<scale factor{1}>] [hfun=<sigm|tanh|lin{sigm}>] [norm=<yes|no{no}>]'
    write(*,'(a)', advance='no') ' [nthr=<nr of openMP threads{1}>] [avg=<yes|no>]'
    write(*,'(a)', advance='no') ' [filetab=<filenames.txt>] [stats=<yes|no{yes}>]'
    write(*,'(a)', advance='no') ' [ctf=<yes|no|flip|mul|abs{no}>] [kv=<acceleration voltage(in kV){300.}>]'
    write(*,'(a)', advance='no') ' [fraca=<frac amp contrast{0.07}>] [cs=<spherical aberration constant(in mm){2.7}>]'
    write(*,'(a)', advance='no') ' [deftab=<text file with defocus values>] [ft2img=<yes|no{no}>] [frameavg=<nr of frames'
    write(*,'(a)', advance='no') ' to average{0}>] [clip=<clipped box size{box}>] [compare=<yes|no{no}>] [mirr=<no|x|y{no}>]'
    write(*,'(a)', advance='no') ' [neg=<yes|no{no}>] [box=<image size(in pixels)>] [outfile=<output_params.txt>]'
    write(*,'(a)', advance='no') ' [ctfsq=<yes|no{no}>] '
    write(*,'(a)') ' [inner=<inner mask radius(in pixels)>] [width=<pixels falloff inner mask{10}>] [append=<yes|no{no}>]'
    stop
endif
call cline%parse
if( .not. cline%defined('thres') )then
    call cline%set('thres', 0.6)
endif
call cline%set('prg', 'stackops')
p = params(cline,checkdistr=.false.)            ! constants & derived constants produced
call b%build_general_tbox(p, cline)             ! general objects built
call img%new([p%box,p%box,1],p%smpd,p%imgkind)  ! image created
call img2%new([p%box,p%box,1],p%smpd,p%imgkind) ! image created

! STACK THE INDIVDUAL PATTERNS IN FILETAB
if( cline%defined('filetab') .and. p%xfel .eq. 'yes')then
    if( .not. cline%defined('box') )then
        stop 'need box to be defined to know the size of the patterns'
    endif
    npatterns = nlines(p%filetab)
    funit_patterns = get_fileunit()
    open(unit=funit_patterns, status='old', action='read', file=p%filetab)
    do i=1,npatterns
        read(funit_patterns,'(a256)') patternname
        if( fname2format(patternname) .eq. 'D' )then
            call b%img%read(patternname, isxfel=.true.)
            call b%img%write(p%outstk, i)
        else
            stop 'cannot recognize the inputted format of the xfel patterns; simple_stackops'
        endif
    end do
    goto 999
endif

! STACK THE SUBSTACKS IN FILETAB
if( cline%defined('filetab') )then
    nstacks = nlines(p%filetab)
    funit_stacks = get_fileunit()
    open(unit=funit_stacks, status='old',  action='read', file=p%filetab)
    if( cline%defined('deftab') )then
        ndefs = nlines(p%deftab)
        if( nstacks /= ndefs ) stop 'number of entries in inputted files (filetab/deftab) do not match!'
        call b%a%new(nstacks)
        call b%a%read(p%deftab)
        ! count the number of images
        cnt = 0
        do i=1,nstacks
            read(funit_stacks,'(a256)') stackname
            call find_ldim_nptcls(stackname,lfoo,nimgs)
            cnt = cnt+nimgs
        end do
        rewind(funit_stacks)
        o_here = oris(cnt)
    endif
    read(funit_stacks,'(a256)') stackname
    rewind(funit_stacks)
    imgformat = fname2format(stackname)
    rawproc = .false.
    if( imgformat .eq. 'B' .or. imgformat .eq. 'D' )then
        if( .not. cline%defined('box') )then
            stop 'need inputted box size for conversion of raws'
        endif
        rawproc = .true.
    else
        call find_ldim_nptcls(stackname,ldim,ifoo)
        ! rebuild b%img
        call b%img%new([ldim(1),ldim(2),1],1.)
    endif
    cnt = 0
    do i=1,nstacks
        read(funit_stacks,'(a256)') stackname
        if( .not. rawproc )then
            call find_ldim_nptcls(stackname,lfoo,nimgs)
            if( cline%defined('deftab') ) o = b%a%get_ori(i)
            do j=1,nimgs
                cnt = cnt+1
                if( cline%defined('deftab') ) call o_here%set_ori(cnt, o)
                call b%img%read(stackname, j)
                call b%img%write(p%outstk, cnt)
            end do
        else
            print *, 'trying to read: ', stackname
            call b%img%new([p%box,p%box,1],1.)
            call b%img%read(stackname,1)
            print *, 'trying to write image:', i, 'to: ', p%outstk
            call b%img%write(p%outstk, i)
            print *, 'did write'
        endif
        call progress(i,nstacks)
    end do
    if( cline%defined('deftab') ) call o_here%write(p%outfile)
    goto 999
endif

! MAKE NOISE IMAGES
if( .not. cline%defined('stk') .and. p%ctf .eq. 'no' )then
    do i=1,p%nptcls
        call img%gauran(0., 1.)
        call img%write(p%outstk, i)
    end do
    goto 999
endif

! APPEND STK2 TO STK WHILE PRESERVING THE NAME OF STK
if( p%append .eq. 'yes' )then
    if( cline%defined('stk') .and. cline%defined('stk2') )then
        ! find out image dimension and number of particles
        call find_ldim_nptcls(p%stk,lfoo,np1)
        call find_ldim_nptcls(p%stk2,lfoo,np2)
        ntot = np1+np2
        cnt = 0
        do i=np1+1,ntot
            cnt = cnt+1
            call progress(cnt,np2)
            call img%read(p%stk2,cnt)
            call img%write(p%stk,i)
        end do
    else
        stop 'need two stacks (stk & stk2) to append; simple_stackops'
    endif
    goto 999
endif

! SPLIT THE STACK IN TO p%split PARTITIONS FOR PARALLEL EXECUTION
if( cline%defined('split') )then
    call splitstk_in_parts(p%stk,p%split)
endif

! CREATE FRAME AVERAGES
if( p%frameavg > 0 )then
    call frameavg_imgfile(p%stk,p%outstk,p%frameavg)
    goto 999
endif

! VISUALIZE
if( p%vis .eq. 'yes' )then
    do i=1,p%nptcls
        call img%read(p%stk, i)
        call img%vis
    end do
    goto 999
endif

! AVERAGE
if( p%avg .eq. 'yes' )then
    call make_avg_imgfile(p%stk, p%outstk)
    goto 999
endif

! RANDOM SELECTION
if( cline%defined('nran') )then
    write(*,'(a)') '>>> RANDOMLY SELECTING IMAGES'
    allocate( pinds(p%nran), stat=alloc_stat )
    call alloc_err('In: simple_stackops', alloc_stat)
    rt = ran_tabu(p%nptcls)
    call rt%ne_ran_iarr(pinds)
    do i=1,p%nran
        call progress(i, p%nran)
        call img%read(p%stk, pinds(i))
        call img%write(p%outstk, i)
    end do
    goto 999
endif

! COPYING
if( cline%defined('top') .and. .not. cline%defined('part') )then
!    call copy_imgfile(p%stk, p%outstk, fromto=[p%fromp,p%top])
    call copy_imgfile(p%stk, p%outstk, fromto=[p%fromp,p%top], smpd_in=p%smpd)
    goto 999
endif

! FISHING EXPEDITIONS
! frac only
if( cline%defined('frac') )then
    if( p%oritab == '' ) stop 'need input orientation doc for fishing expedition; simple_stackops'
    ! determine how many particles to include
    if( p%frac < 0.99 )then
        nincl = nint(real(p%nptcls)*p%frac)
    else
        nincl = p%nptcls
    endif
    ! order the particles
    pinds = b%a%order()
    ! fish the best ones out
    if( cline%defined('state') )then
        cnt = 0
        do i=1,nincl
            call progress(i, nincl)
            s = nint(b%a%get(pinds(i), 'state'))
            if( s == p%state )then
                cnt = cnt+1
                call img%read(p%stk, pinds(i))
                call img%write(p%outstk, cnt)
            endif
        end do
        ! make orientation structure for the best ones
        o_here = oris(cnt)
        cnt = 0
        do i=1,nincl
            call progress(i, nincl)
            s = nint(b%a%get(pinds(i), 'state'))
            if( s == p%state )then
                cnt = cnt+1
                call o_here%set_ori(cnt, b%a%get_ori(pinds(i)))
            endif
        end do
        allocate(fname, source='extracted_oris_state'//int2str_pad(p%state,2)//'.txt')
    else if( cline%defined('class') )then
        cnt = 0
        do i=1,nincl
            call progress(i, nincl)
            s = nint(b%a%get(pinds(i), 'class'))
            if( s == p%class )then
                cnt = cnt+1
                call img%read(p%stk, pinds(i))
                call img%write(p%outstk, cnt)
            endif
        end do
        ! make orientation structure for the best ones
        o_here = oris(cnt)
        cnt = 0
        do i=1,nincl
            call progress(i, nincl)
            s = nint(b%a%get(pinds(i), 'class'))
            if( s == p%class )then
                cnt = cnt+1
                call o_here%set_ori(cnt, b%a%get_ori(pinds(i)))
            endif
        end do
        allocate(fname, source='extracted_oris_class'//int2str_pad(p%state,2)//'.txt')
    else
        o_here = oris(nincl)
        do i=1,nincl
            call progress(i, nincl)
            call img%read(p%stk, pinds(i))
            call img%write(p%outstk, i)
            call o_here%set_ori(i, b%a%get_ori(pinds(i)))
        end do
        allocate(fname, source='extracted_oris.txt')
    endif
    call o_here%write(fname)
    deallocate(fname)
    goto 999
endif
! state/class + frac
if( (cline%defined('state') .or. cline%defined('class')) .and. .not.  cline%defined('frac') )then
    if( p%oritab == '' ) stop 'need input orientation doc for fishing expedition; simple_stackops'
    if( cline%defined('state') )then
        cnt = 0
        do i=1,p%nptcls
            call progress(i, p%nptcls)
            s = nint(b%a%get(i, 'state'))
            if( s == p%state )then
                cnt = cnt+1
                call img%read(p%stk, i)
                call img%write(p%outstk, cnt)
            endif
        end do
        ! make orientation structure for the extracted ones
        o_here = oris(cnt)
        cnt = 0
        do i=1,p%nptcls
            call progress(i, p%nptcls)
            s = nint(b%a%get(i, 'state'))
            if( s == p%state )then
                cnt = cnt+1
                call o_here%set_ori(cnt, b%a%get_ori(i))
            endif
        end do
        allocate(fname, source='extracted_oris_state'//int2str_pad(p%state,2)//'.txt')
    else if( cline%defined('class') )then
        cnt = 0
        do i=1,p%nptcls
            call progress(i, p%nptcls)
            s = nint(b%a%get(i, 'class'))
            if( s == p%class )then
                cnt = cnt+1
                call img%read(p%stk, i)
                call img%write(p%outstk, cnt)
            endif
        end do
        ! make orientation structure for the best ones
        o_here = oris(cnt)
        cnt = 0
        do i=1,p%nptcls
            call progress(i, p%nptcls)
            s = nint(b%a%get(i, 'class'))
            if( s == p%class )then
                cnt = cnt+1
                call o_here%set_ori(cnt, b%a%get_ori(i))
            endif
        end do
        allocate(fname, source='extracted_oris_class'//int2str_pad(p%state,2)//'.txt')
    endif
    call o_here%write(fname)
    deallocate(fname)
    goto 999
endif

! ADD NOISE
if( cline%defined('snr') )then
    call add_noise_imgfile(p%stk, p%outstk, p%snr)
    goto 999
endif

! FILTERING
if( (cline%defined('lp') .and. p%phrand .eq. 'no') .or. cline%defined('hp') )then
    if( .not. cline%defined('smpd') ) stop 'need smpd 4  filtering!'
    if( cline%defined('lp') .and. cline%defined('hp') )then
        call bp_imgfile(p%stk, p%outstk, p%smpd, p%hp, p%lp)
    else if( cline%defined('lp') )then
        call bp_imgfile(p%stk, p%outstk, p%smpd, 0., p%lp)
    endif
    goto 999
endif

! VISUALIZE FT
if( p%ft2img .eq. 'yes' )then
    call ft2img_imgfile(p%stk, p%outstk, 'power')
    goto 999
endif

! PHASE RANDOMIZATION
if( p%phrand .eq. 'yes' )then
    if( .not. cline%defined('smpd') ) stop 'need smpd 4 phase randomization!'
    if( .not. cline%defined('lp') )   stop 'low-pass limit needed 4 phase randomization'
    call phase_rand_imgfile(p%stk, p%outstk, p%smpd, p%lp)
    goto 999
endif

! PHASE FLIPPING
if( p%ctf .ne. 'no' )then
    if( debug )then
        write(*,*) 'CTF parameters used in simple_stackops'
        write(*,*) 'kv = ', p%kv
        write(*,*) 'cs = ', p%cs
        write(*,*) 'fraca = ', p%fraca
    endif
    if( .not. cline%defined('smpd') ) stop 'need smpd 2 apply CTF!'
    if( cline%defined('oritab') .or. cline%defined('deftab') )then
        ! all good
    else
        stop 'oritab/deftab with CTF info needed for phase flipping/multiplication'
    endif
    if( .not. b%a%isthere('dfx') ) stop 'need at least dfx to be set if ctf is going to be used!'
endif
if( p%ctf .eq. 'flip' )then
    if( p%neg .eq. 'yes' )then
        call apply_ctf_imgfile(p%stk, p%outstk, b%a, p%smpd, b%tfun, 'flipneg')
    else
        call apply_ctf_imgfile(p%stk, p%outstk, b%a, p%smpd, b%tfun, 'flip')
    endif
    goto 999
else if( p%ctf .eq. 'mul' )then
    if( p%neg .eq. 'yes' )then
        call apply_ctf_imgfile(p%stk, p%outstk, b%a, p%smpd, b%tfun, 'neg')
    else
        call apply_ctf_imgfile(p%stk, p%outstk, b%a, p%smpd, b%tfun, 'ctf')
    endif
    goto 999
else if( p%ctf .eq. 'abs' )then
    call apply_ctf_imgfile(p%stk, p%outstk, b%a, p%smpd, b%tfun, 'abs')
    goto 999  
else if( p%ctf .eq. 'wiener' )then
    
    print *, 'DOING THE WIENER THING'
    
    call apply_wiener_imgfile(p%stk, p%outstk, b%a, p%smpd, b%tfun)    
    goto 999
endif

! APPLY CFTSQ
if( p%ctfsq .eq. 'yes' )then
    if( .not. cline%defined('smpd') ) stop 'need smpd 2 apply ctf**2!'
    if( .not. cline%defined('oritab') ) stop 'oritab with CTF info needed'
    if( .not. b%a%isthere('dfx') ) stop 'need at least dfx to be set if ctf is going to be used!'
    if( cline%defined('bfac') )then
        call apply_ctf_imgfile(p%stk, p%outstk, b%a, p%smpd, b%tfun, 'square', bfac=p%bfac)
    else
        call apply_ctf_imgfile(p%stk, p%outstk, b%a, p%smpd, b%tfun, 'square')
    endif
    goto 999
endif

! INVERT CONTRAST
if( p%neg .eq. 'yes' )then
    call neg_imgfile(p%stk, p%outstk)
    goto 999
endif

! NOISE NORMALIZATION
if( p%noise_norm .eq. 'yes' )then
    if( cline%defined('msk') )then
        call noise_norm_imgfile(p%stk, p%msk, p%outstk)
    else
        stop 'need msk parameter for noise normalization'
    endif
    goto 999
endif

! MASKING
if( cline%defined('msk') .or. cline%defined('inner') )then
    if( cline%defined('inner') )then
        if( cline%defined('width') )then
            call mask_imgfile(p%stk, p%outstk, p%msk, inner=p%inner, width=p%width)
        else
            call mask_imgfile(p%stk, p%outstk, p%msk, inner=p%inner)
        endif
    else
        call mask_imgfile(p%stk, p%outstk, p%msk)
    endif
    goto 999
endif

! NORMALIZATION
if( cline%defined('norm') )then
    if( cline%defined('hfun') )then
        call norm_imgfile(p%stk, p%outstk, hfun=p%hfun)
    else
        call norm_imgfile(p%stk, p%outstk)
    endif
    goto 999
endif

! CURING
if( cline%defined('cure') )then
    call cure_imgfile(p%stk, p%outstk)
    goto 999
endif

! AUTO CORRELATION FUNCTION
if( p%acf .eq. 'yes' )then
    call acf_imgfile(p%stk, p%outstk)
    goto 999
endif

! SCALING
if( cline%defined('newbox') .or. cline%defined('scale') )then
    ldim_scaled = [p%newbox,p%newbox,1] ! dimension of scaled
    if( cline%defined('clip') )then
        if( cline%defined('part') )then
            p%outstk = 'outstk_part'//int2str_pad(p%part, p%numlen)//p%ext
            call resize_and_clip_imgfile(p%stk,p%outstk,ldim_scaled,[p%clip,p%clip,1],[p%fromp,p%top])
        else
            call resize_and_clip_imgfile(p%stk,p%outstk,ldim_scaled,[p%clip,p%clip,1])
        endif
        goto 999
    else
        if( cline%defined('part') )then
            p%outstk = 'outstk_part'//int2str_pad(p%part, p%numlen)//p%ext
            call resize_imgfile(p%stk,p%outstk,ldim_scaled,[p%fromp,p%top])
        else
            call resize_imgfile(p%stk,p%outstk,ldim_scaled)
        endif
        goto 999
    endif
endif

! MIRRORING
if( cline%defined('mirr') )then
    call mirror_imgfile(p%stk,p%outstk,p%mirr)
    goto 999
endif

! CLIPPING
if( cline%defined('clip') )then
    call clip_imgfile(p%stk,p%outstk,[p%clip,p%clip,1])
    goto 999
endif

999 call simple_end('**** SIMPLE_STACKOPS NORMAL STOP ****')
call del_binfile('tmp'//p%ext)

!*******************************************************************************
!    Environment shutdown
!
!*******************************************************************************
!shutting down the timers
call stop_Alltimers_cpu()

end program simple_stackops
