!==Class simple_params
!
! simple_params provides global distribution of constants and derived constants used in the SIMPLE library.
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_. Redistribution
! or modification is regulated by the GNU General Public License. *Author:* Hans Elmlund, 2011-08-18.
!
!==Changes are documented below
!
!* rewrote as a proper class, HE June 7 2012
!
module simple_params
use simple_defs
use simple_ori,         only: ori
use simple_AVratios,    only: AVratios
use simple_cmdline,     only: cmdline
use simple_strings      ! use all in there
use simple_filehandling ! use all in there
use simple_jiffys,      only: find_ldim_nptcls
use simple_magic_boxes, only: find_magic_box
implicit none

public :: params
private

logical, parameter :: debug=.false.

type :: params
    ! global objects
    type(ori)        :: ori_glob
    type(AVratios)   :: avr
    type(ctfplan)    :: tfplan
    ! yes/no decision variables in ascending alphabetical order
    character(len=3) :: acf='no'
    character(len=3) :: append='no'
    character(len=3) :: async='no'
    character(len=3) :: automsk='no'   
    character(len=3) :: avg='no'
    character(len=3) :: bench_gpu='no'
    character(len=3) :: bin='no'
    character(len=3) :: center='no'
    character(len=3) :: clustvalid='no'
    character(len=3) :: compare='no'
    character(len=3) :: countvox='no'
    character(len=3) :: ctfstats='no'
    character(len=3) :: cure='no'
    character(len=3) :: debug='no'
    character(len=3) :: discrete='no'
    character(len=3) :: diverse='no'
    character(len=3) :: doalign='yes'
    character(len=3) :: dopca='yes'
    character(len=3) :: dopick='yes'
    character(len=3) :: doprint='no'
    character(len=3) :: dynlp='yes'
    character(len=3) :: eo='yes'
    character(len=3) :: errify='no'
    character(len=3) :: even='no'
    character(len=3) :: fix_gpu='no'
    character(len=3) :: ft2img='no'
    character(len=3) :: guinier='no'
    character(len=3) :: kmeans='yes'
    character(len=3) :: local='no'
    character(len=3) :: masscen='no'
    character(len=3) :: merge='no'
    character(len=3) :: mirr='no'
    character(len=3) :: neg='no'
    character(len=3) :: noise_norm ='no'
    character(len=3) :: noise='no'
    character(len=3) :: norec='no'
    character(len=3) :: norm='no'
    character(len=3) :: odd='no'
    character(len=3) :: order='no'
    character(len=3) :: outside='no'
    character(len=3) :: pad='no'
    character(len=3) :: phaseplate='no'
    character(len=3) :: phrand='no'
    character(len=3) :: plot='no'
    character(len=3) :: readwrite='no'
    character(len=3) :: remap='no'
    character(len=3) :: restart='no'
    character(len=3) :: rnd='no'
    character(len=3) :: rm_outliers='yes'
    character(len=3) :: roalgn='no'
    character(len=3) :: round='no'
    character(len=3) :: shalgn='no'
    character(len=3) :: shellnorm='no'
    character(len=3) :: shellw='no'
    character(len=3) :: shbarrier='yes'
    character(len=3) :: single='no'
    character(len=3) :: soften='no'
    character(len=3) :: srch_inpl='yes'
    character(len=3) :: stats='no'
    character(len=3) :: stream='no'
    character(len=3) :: swap='no'
    character(len=3) :: test='no'
    character(len=3) :: tomo='no'
    character(len=3) :: time='no'
    character(len=3) :: trsstats='no'
    character(len=3) :: tseries='no'
    character(len=3) :: use_gpu='no'
    character(len=3) :: verbose='no'
    character(len=3) :: vis='no'
    character(len=3) :: xfel='no'
    character(len=3) :: zero='no'
    ! other fixed length character variables in ascending alphabetical order
    character(len=STDLEN) :: angastunit='degrees'
    character(len=STDLEN) :: boxfile=''
    character(len=STDLEN) :: boxtab=''
    character(len=STDLEN) :: boxtype='eman'
    character(len=STDLEN) :: clsdoc=''
    character(len=STDLEN) :: comlindoc=''
    character(len=STDLEN) :: ctf='no'
    character(len=STDLEN) :: cwd=''
    character(len=STDLEN) :: deftab=''
    character(len=STDLEN) :: dfunit='microns'
    character(len=STDLEN) :: dir=''
    character(len=STDLEN) :: dir_movies=''
    character(len=STDLEN) :: dir_reject='rejected'
    character(len=STDLEN) :: dir_select='selected'
    character(len=STDLEN) :: dir_target=''
    character(len=STDLEN) :: dir_ptcls=''
    character(len=STDLEN) :: doclist=''
    character(len=STDLEN) :: endian='native'
    character(len=STDLEN) :: exp_doc=''
    character(len=4)      :: ext='.mrc'
    character(len=STDLEN) :: extrmode='all'
    character(len=STDLEN) :: fbody=''
    character(len=STDLEN) :: featstk='expecstk.bin'
    character(len=STDLEN) :: filetab=''
    character(len=STDLEN) :: fname=''
    character(len=STDLEN) :: fsc='fsc_state01.bin'
    character(len=STDLEN) :: hfun='sigm'
    character(len=STDLEN) :: hist='corr'
    character(len=STDLEN) :: imgkind='em'
    character(len=STDLEN) :: infile='infile.txt'
    character(len=STDLEN) :: label='class'
    character(len=STDLEN) :: masks(MAXS)=''
    character(len=STDLEN) :: mskfile=''
    character(len=STDLEN) :: msktype='soft'
    character(len=STDLEN) :: opt='simplex'
    character(len=STDLEN) :: oritab=''
    character(len=STDLEN) :: oritab2=''
    character(len=STDLEN) :: outfile='outfile.txt'
    character(len=STDLEN) :: outstk=''
    character(len=STDLEN) :: outstk2=''
    character(len=STDLEN) :: outvol=''
    character(len=STDLEN) :: ctffind_doc=''
    character(len=STDLEN) :: pcastk='pcavecinstk.bin'
    character(len=STDLEN) :: pdfile='pdfile.bin'
    character(len=STDLEN) :: pgrp='c1'
    character(len=STDLEN) :: plaintexttab=''
    character(len=STDLEN) :: prg=''
    character(len=STDLEN) :: refine='no'
    character(len=STDLEN) :: refs_msk=''
    character(len=STDLEN) :: refs=''
    character(len=STDLEN) :: speckind='sqrt'
    character(len=STDLEN) :: split_mode='even'
    character(len=STDLEN) :: stk_part=''
    character(len=STDLEN) :: stk=''
    character(len=STDLEN) :: stk2=''
    character(len=STDLEN) :: stk3=''
    character(len=STDLEN) :: tomoseries=''
    character(len=STDLEN) :: vol=''
    character(len=STDLEN) :: vollist=''
    character(len=STDLEN) :: vols_msk(MAXS)=''
    character(len=STDLEN) :: vols(MAXS)=''
    character(len=STDLEN) :: voltab=''
    character(len=STDLEN) :: voltab2=''
    character(len=STDLEN) :: wfun='kb'
    ! integer variables in ascending alphabetical order
    integer :: astep=1
    integer :: avgsz=0
    integer :: binwidth=1
    integer :: box=0
    integer :: boxconvsz=256
    integer :: boxmatch=0
    integer :: boxpd=0
    integer :: chunksz=0
    integer :: class=1
    integer :: clip=0
    integer :: corner=0
    integer :: cube=0
    integer :: edge=14
    integer :: filtsz_pad=0
    integer :: find=1
    integer :: frameavg=0
    integer :: fromf=1
    integer :: fromp=1
    integer :: froms=1
    integer :: fstep=1
    integer :: grow=0
    integer :: iares=10
    integer :: ind=0
    integer :: iptcl=1
    integer :: jptcl=1
    integer :: jumpsz=0
    integer :: kfromto(2)
    integer :: ldim(3)=0
    integer :: maxits=500
    integer :: maxp=0
    integer :: minp=10
    integer :: mrcmode=2
    integer :: navgs=1
    integer :: ncunits=0
    integer :: nbest=100
    integer :: nboot=0
    integer :: ncls=500
    integer :: ncomps=0
    integer :: ndiscrete=0
    integer :: ndocs=0
    integer :: newbox=0
    integer :: newbox2=0
    integer :: nframes=0
    integer :: nmembers=0
    integer :: nnn=500
    integer :: noris=0
    integer :: nparts=1
    integer :: npeaks=1
    integer :: npix=0
    integer :: nptcls=1
    integer :: nran=0
    integer :: nrefs=100
    integer :: nrestarts=1
    integer :: nrots=0
    integer :: nspace=1000
    integer :: nstates=1
    integer :: nsym=1
    integer :: nthr=1
    integer :: nthr_master=1
    integer :: numlen=0
    integer :: nvalid=0
    integer :: nvars=30
    integer :: nvox=0
    integer :: offset=7
    integer :: part=1
    integer :: pcasz=0
    integer :: ppca=0
    integer :: pspecsz=512
    integer :: pspecsz_unblur=512
    integer :: pspecsz_ctffind=1024
    integer :: ptcl=1
    integer :: ring1=2
    integer :: ring2=0
    integer :: set_gpu=0
    integer :: spec=0
    integer :: startit=1
    integer :: state=1
    integer :: state2split=0
    integer :: stepsz=1
    integer :: tofny=0
    integer :: tof=1
    integer :: top=1
    integer :: tos=1
    integer :: trsstep=1
    integer :: update=1000
    integer :: which_iter=0
    integer :: xcoord=0
    integer :: ycoord=0
    integer :: xdim=0
    integer :: xdimpd=0
    integer :: ydim=0
    ! real variables in ascending alphabetical order  
    real    :: alpha=2.
    real    :: amsklp=20.
    real    :: angerr=0.
    real    :: ares=7.
    real    :: astigerr=0.
    real    :: astigstep=0.05
    real    :: athres=0.
    real    :: bfac=200
    real    :: bfacerr=50.
    real    :: cenlp=50.
    real    :: cs=2.7
    real    :: ctfreslim=8.
    real    :: dcrit_rel=0.5
    real    :: deflim=4.
    real    :: defocus=3.
    real    :: dens=0.
    real    :: dferr=1.
    real    :: dfmax=7.0
    real    :: dfmin=0.5
    real    :: dfsdev=0.1
    real    :: dose_rate=30.0
    real    :: dstep=0.
    real    :: dsteppd=0.
    real    :: e1=0.
    real    :: e2=0.
    real    :: e3=0.
    real    :: eps=0.003
    real    :: extr_thresh=EXTRINITHRESH
    real    :: eullims(3,2)=0.
    real    :: expastig=0.1
    real    :: exp_time=2.0
    real    :: filwidth=0.
    real    :: fny=0.
    real    :: frac=1.
    real    :: fraca=0.07
    real    :: fracdeadhot=0.05
    real    :: fraczero=0.
    real    :: ftol=1e-6
    real    :: gw=0.5
    real    :: het_thresh=HETINITTHRESH
    real    :: hp=100.
    real    :: hp_ctffind=30.
    real    :: inner=0.
    real    :: kv=300.
    real    :: lam=0.5
    real    :: lp_dyn=20.
    real    :: lp=20.
    real    :: lp_ctffind=5.0
    real    :: lp_pick=20.
    real    :: lpmed=20.
    real    :: lpstart=0.
    real    :: lpstop=7.0
    real    :: lpvalid=20.
    real    :: moldiam=140.
    real    :: moment=0.
    real    :: msk=0.
    real    :: mul=1.
    real    :: mw=0.
    real    :: neigh=0.2
    real    :: optlims(7,2)=0.
    real    :: outer=0.
    real    :: phranlp=35.
    real    :: power=2.
    real    :: rrate=0.8
    real    :: scale=1.
    real    :: scale2=1.
    real    :: sherr=0.
    real    :: smpd=2.
    real    :: snr
    real    :: thres=0.
    real    :: time_per_image=200.
    real    :: time_per_frame=0.
    real    :: trs=0.
    real    :: var=1.
    real    :: width=10.
    real    :: winsz=1.
    real    :: xsh=0.
    real    :: ysh=0.
    real    :: zsh=0.
    ! logical variables in ascending alphabetical order
    logical :: cyclic(7)     = .false.
    logical :: l_distr_exec  = .false.
    logical :: doautomsk     = .false.
    logical :: doshift       = .false.
    logical :: l_automsk     = .false.
    logical :: l_dose_weight = .false. 
    logical :: l_innermsk    = .false. 
    logical :: l_pick        = .false. 
    logical :: l_shellw      = .false.
    logical :: l_xfel        = .false.
  contains
    procedure :: new
end type params

interface params
    module procedure constructor
end interface

contains
    
    !> \brief  is a constructor
    function constructor( cline, checkdistr, allow_mix ) result( self )
        class(cmdline), intent(inout) :: cline
        logical, intent(in), optional :: checkdistr, allow_mix
        type(params) :: self
        call self%new( cline, checkdistr, allow_mix )    
    end function constructor

    !> \brief  is a constructor
    subroutine new( self, cline, checkdistr, allow_mix )
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_math, only: round2even
        use simple_map_reduce  ! use all in there
        class(params),     intent(inout) :: self
        class(cmdline),    intent(inout) :: cline
        logical, optional, intent(in)    :: checkdistr, allow_mix
        integer                          :: i, s, ncls, ifoo, lfoo(3), cntfile=0, nthr_local
        logical                          :: here, ccheckdistr, aamix
        character(len=STDLEN)            :: cwd_local, debug_local
        character(len=1)                 :: checkupfile(50)
        character(len=:), allocatable    :: conv
        integer, allocatable             :: parts(:,:)
        logical                          :: nparts_set=.false.
        logical                          :: vol_defined(MAXS)=.false.
        ! take care of optionals
        ccheckdistr = .true.
        if( present(checkdistr) ) ccheckdistr = checkdistr
        aamix = .false.
        if( present(allow_mix) ) aamix = allow_mix
        ! file counter
        cntfile = 0
        ! make global ori
        call self%ori_glob%new
        ! get cwd
        call getcwd(self%cwd)
        cwd_local = self%cwd
        ! checkers in ascending alphabetical order
        call check_carg('acf',            self%acf)
        call check_carg('angastunit',     self%angastunit)
        call check_carg('append',         self%append)
        call check_carg('async',          self%async)
        call check_carg('automsk',        self%automsk)
        call check_carg('avg',            self%avg)
        call check_carg('bench_gpu',      self%bench_gpu)
        call check_carg('bin',            self%bin)
        call check_carg('boxtype',        self%boxtype)
        call check_carg('center',         self%center)
        call check_carg('clustvalid',     self%clustvalid)
        call check_carg('compare',        self%compare)
        call check_carg('countvox',       self%countvox)
        call check_carg('ctf',            self%ctf)
        call check_carg('ctfstats',       self%ctfstats)
        call check_carg('cure',           self%cure)
        call check_carg('debug',          self%debug)
        debug_local = self%debug
        call check_carg('deftab',         self%deftab)
        call check_carg('dfunit',         self%dfunit)
        call check_carg('dir',            self%dir)
        call check_carg('dir_movies',     self%dir_movies)
        call check_carg('dir_ptcls',      self%dir_ptcls)
        call check_carg('dir_reject',     self%dir_reject)
        call check_carg('dir_select',     self%dir_select)
        call check_carg('dir_target',     self%dir_target)
        call check_carg('discrete',       self%discrete)
        call check_carg('diverse',        self%diverse)
        call check_carg('doalign',        self%doalign)
        call check_carg('dopca',          self%dopca)
        call check_carg('dopick',         self%dopick)
        call check_carg('doprint',        self%doprint)
        call check_carg('dynlp',          self%dynlp)
        call check_carg('endian',         self%endian)
        call check_carg('eo',             self%eo)
        call check_carg('errify',         self%errify)
        call check_carg('even',           self%even)
        call check_carg('exp_doc',        self%exp_doc)
        call check_carg('extrmode',       self%extrmode)
        call check_carg('fbody',          self%fbody)
        call check_carg('fix_gpu',        self%fix_gpu)
        call check_carg('ft2img',         self%ft2img)
        call check_carg('guinier',        self%guinier)
        call check_carg('hfun',           self%hfun)
        call check_carg('hist',           self%hist)
        call check_carg('kmeans',         self%kmeans)
        call check_carg('label',          self%label)
        call check_carg('local',          self%local)
        call check_carg('masscen',        self%masscen)
        call check_carg('merge',          self%merge)
        call check_carg('mirr',           self%mirr)
        call check_carg('msktype',        self%msktype)
        call check_carg('neg',            self%neg)
        call check_carg('noise_norm',     self%noise_norm)
        call check_carg('noise',          self%noise)
        call check_carg('norec',          self%norec)
        call check_carg('norm',           self%norm)
        call check_carg('odd',            self%odd)
        call check_carg('opt',            self%opt)
        call check_carg('order',          self%order)
        call check_carg('outfile',        self%outfile)
        call check_carg('outside',        self%outside)
        call check_carg('pad',            self%pad)
        call check_carg('ctffind_doc',    self%ctffind_doc)
        call check_carg('pgrp',           self%pgrp)
        call check_carg('phaseplate',     self%phaseplate)
        call check_carg('phrand',         self%phrand)
        call check_carg('prg',            self%prg)
        call check_carg('readwrite',      self%readwrite)
        call check_carg('refine',         self%refine)
        call check_carg('refs',           self%refs)
        call check_carg('remap',          self%remap)
        call check_carg('restart',        self%restart)
        call check_carg('rm_outliers',    self%rm_outliers)
        call check_carg('rnd',            self%rnd)
        call check_carg('roalgn',         self%roalgn)
        call check_carg('round',          self%round)
        call check_carg('shalgn',         self%shalgn)
        call check_carg('shbarrier',      self%shbarrier)
        call check_carg('shellnorm',      self%shellnorm)
        call check_carg('shellw',         self%shellw)
        call check_carg('single',         self%single)
        call check_carg('soften',         self%soften)
        call check_carg('speckind',       self%speckind)
        call check_carg('split_mode',     self%split_mode)
        call check_carg('srch_inpl',      self%srch_inpl)
        call check_carg('stats',          self%stats)
        call check_carg('stream',         self%stream)
        call check_carg('swap',           self%swap)
        call check_carg('test',           self%test)
        call check_carg('time',           self%time)
        call check_carg('tomo',           self%tomo)
        call check_carg('tomoseries',     self%tomoseries)
        call check_carg('trsstats',       self%trsstats)
        call check_carg('tseries',        self%tseries)
        call check_carg('use_gpu',        self%use_gpu)
        call check_carg('verbose',        self%verbose)
        call check_carg('vis',            self%vis)
        call check_carg('vol',            self%vol)
        call check_carg('wfun',           self%wfun)
        call check_carg('xfel',           self%xfel)
        call check_carg('zero',           self%zero)
        call check_file('boxfile',        self%boxfile,'T')
        call check_file('boxtab',         self%boxtab,'T')
        call check_file('clsdoc',         self%clsdoc,'S','T')
        call check_file('comlindoc',      self%comlindoc,'T')
        call check_file('doclist',        self%doclist,'T')
        call check_file('ext',            self%ext,  notAllowed='T')
        call check_file('filetab',        self%filetab,'T')
        call check_file('fname',          self%fname)
        call check_file('fsc',            self%fsc,'B')
        call check_file('infile',         self%infile)
        call check_file('mskfile',        self%mskfile,  notAllowed='T')
        call check_file('oritab',         self%oritab, 'T')
        call check_file('oritab2',        self%oritab2,'T')
        call check_file('outstk',         self%outstk,   notAllowed='T')
        call check_file('outstk2',        self%outstk2,  notAllowed='T')
        call check_file('outvol',         self%outvol,   notAllowed='T')
        call check_file('plaintexttab',   self%plaintexttab,'T')
        call check_file('stk',            self%stk,  notAllowed='T')
        call check_file('stk2',           self%stk2, notAllowed='T')
        call check_file('stk3',           self%stk3, notAllowed='T')
        call check_file('vollist',        self%vollist, 'T')
        call check_file('voltab',         self%voltab, 'T')
        call check_file('voltab2',        self%voltab2,'T')
        call check_iarg('astep',          self%astep)
        call check_iarg('avgsz',          self%avgsz)
        call check_iarg('binwidth',       self%binwidth)
        call check_iarg('box',            self%box)
        call check_iarg('boxconvsz',      self%boxconvsz)
        call check_iarg('chunksz',        self%chunksz)
        call check_iarg('clip',           self%clip)
        call check_iarg('corner',         self%corner)
        call check_iarg('cube',           self%cube)
        call check_iarg('edge',           self%edge)
        call check_iarg('find',           self%find)
        call check_iarg('frameavg',       self%frameavg)
        call check_iarg('fromf',          self%fromf)
        call check_iarg('fromp',          self%fromp)
        call check_iarg('froms',          self%froms)
        call check_iarg('fstep',          self%fstep)
        call check_iarg('grow',           self%grow)
        call check_iarg('iares',          self%iares)
        call check_iarg('ind',            self%ind)
        call check_iarg('jumpsz',         self%jumpsz)
        call check_iarg('maxits',         self%maxits)
        call check_iarg('maxp',           self%maxp)
        call check_iarg('minp',           self%minp)
        call check_iarg('mrcmode',        self%mrcmode)
        call check_iarg('navgs',          self%navgs)
        call check_iarg('nbest',          self%nbest)
        call check_iarg('nboot',          self%nboot)
        call check_iarg('ncls',           self%ncls)
        call check_iarg('ncunits',        self%ncunits)
        call check_iarg('ndiscrete',      self%ndiscrete)
        call check_iarg('ndocs',          self%ndocs)
        call check_iarg('nframes',        self%nframes)
        call check_iarg('nmembers',       self%nmembers)
        call check_iarg('nnn',            self%nnn)
        call check_iarg('noris',          self%noris)
        call check_iarg('nran',           self%nran)
        call check_iarg('nrefs',          self%nrefs)
        call check_iarg('nrestarts',      self%nrestarts)
        call check_iarg('nspace',         self%nspace)
        call check_iarg('nstates',        self%nstates)
        call check_iarg('set_gpu',        self%set_gpu)
        call check_iarg('class',          self%class)
        call check_iarg('nparts',         self%nparts)
        call check_iarg('npeaks',         self%npeaks)
        call check_iarg('npix',           self%npix)
        call check_iarg('nptcls',         self%nptcls)
        call check_iarg('nthr',           self%nthr)
        call check_iarg('nthr_master',    self%nthr_master)
        call check_iarg('numlen',         self%numlen)
        call check_iarg('nvars',          self%nvars)
        call check_iarg('nvox',           self%nvox)
        call check_iarg('offset',         self%offset)
        call check_iarg('part',           self%part)
        call check_iarg('ppca',           self%ppca)
        call check_iarg('pspecsz',        self%pspecsz)
        call check_iarg('pspecsz_unblur', self%pspecsz_unblur)
        call check_iarg('pspecsz_ctffind', self%pspecsz_ctffind)
        call check_iarg('ring1',          self%ring1)
        call check_iarg('ring2',          self%ring2)
        call check_iarg('startit',        self%startit)
        call check_iarg('state',          self%state)
        call check_iarg('state2split',    self%state2split)
        call check_iarg('stepsz',         self%stepsz)
        call check_iarg('tof',            self%tof)
        call check_iarg('top',            self%top)
        call check_iarg('tos',            self%tos)
        call check_iarg('trsstep',        self%trsstep)
        call check_iarg('update',         self%update)
        call check_iarg('which_iter',     self%which_iter)
        call check_iarg('xdim',           self%xdim)
        call check_iarg('xcoord',         self%xcoord)
        call check_iarg('ycoord',         self%ycoord)
        call check_iarg('ydim',           self%ydim)
        call check_rarg('alpha',          self%alpha)
        call check_rarg('amsklp',         self%amsklp)
        call check_rarg('angerr',         self%angerr)
        call check_rarg('ares',           self%ares)
        call check_rarg('astigerr',       self%astigerr)
        call check_rarg('astigstep',      self%astigstep)
        call check_rarg('athres',         self%athres)
        call check_rarg('bfac',           self%bfac)
        call check_rarg('bfacerr',        self%bfacerr)
        call check_rarg('cenlp',          self%cenlp)
        call check_rarg('cs',             self%cs)
        call check_rarg('ctfreslim',      self%ctfreslim)
        call check_rarg('dcrit_rel',      self%dcrit_rel)
        call check_rarg('deflim',         self%deflim)
        call check_rarg('defocus',        self%defocus)
        call check_rarg('dens',           self%dens)
        call check_rarg('dferr',          self%dferr)
        call check_rarg('dfmax',          self%dfmax)
        call check_rarg('dfmin',          self%dfmin)
        call check_rarg('dfsdev',         self%dfsdev)
        call check_rarg('dose_rate',      self%dose_rate)
        call check_rarg('e1',             self%e1)
        call check_rarg('e2',             self%e2)
        call check_rarg('e3',             self%e3)
        call check_rarg('eps',            self%eps)
        call check_rarg('expastig',       self%expastig)
        call check_rarg('exp_time',       self%exp_time)
        call check_rarg('extr_thresh',    self%extr_thresh)
        call check_rarg('filwidth',       self%filwidth)
        call check_rarg('frac',           self%frac)
        call check_rarg('fraca',          self%fraca)
        call check_rarg('fracdeadhot',    self%fracdeadhot)
        call check_rarg('fraczero',       self%fraczero)
        call check_rarg('ftol',           self%ftol)
        call check_rarg('gw',             self%gw)
        call check_rarg('het_thresh',     self%het_thresh)
        call check_rarg('hp',             self%hp)
        call check_rarg('hp_ctffind',     self%hp_ctffind)
        call check_rarg('inner',          self%inner)
        call check_rarg('kv',             self%kv)
        call check_rarg('lam',            self%lam)
        call check_rarg('lp',             self%lp)
        call check_rarg('lp_ctffind',     self%lp_ctffind)
        call check_rarg('lp_pick',        self%lp_pick)
        call check_rarg('lpstart',        self%lpstart)
        call check_rarg('lpstop',         self%lpstop)
        call check_rarg('lpvalid',        self%lpvalid)
        call check_rarg('moldiam',        self%moldiam)
        call check_rarg('msk',            self%msk)
        call check_rarg('mul',            self%mul)
        call check_rarg('mw',             self%mw)
        call check_rarg('neigh',          self%neigh)
        call check_rarg('outer',          self%outer)
        call check_rarg('phranlp',        self%phranlp)
        call check_rarg('power',          self%power)
        call check_rarg('rrate',          self%rrate)
        call check_rarg('scale',          self%scale)
        call check_rarg('scale2',         self%scale2)
        call check_rarg('sherr',          self%sherr)
        call check_rarg('smpd',           self%smpd)
        call check_rarg('snr',            self%snr)
        call check_rarg('thres',          self%thres)
        call check_rarg('time_per_image', self%time_per_image)
        call check_rarg('trs',            self%trs)
        call check_rarg('var',            self%var)
        call check_rarg('width',          self%width)
        call check_rarg('winsz',          self%winsz)
        call check_rarg('xsh',            self%xsh)
        call check_rarg('ysh',            self%ysh)
        call check_rarg('zsh',            self%zsh)
        ! put ctffind_doc (if defined) as oritab
        if( cline%defined('ctffind_doc') )then
            if( .not. cline%defined('oritab') )then
                call cline%set('oritab', self%ctffind_doc)
                self%oritab = self%ctffind_doc
            endif
        endif
        ! make all programs have the simple_ prefix
        if( cline%defined('prg') )then
            if( .not. str_has_substr(self%prg, 'simple_') ) self%prg = 'simple_'//trim(self%prg)
        endif
        ! check nr of states
        if( cline%defined('nstates') )then
            self%nstates = nint(cline%get_rarg('nstates'))
            if( self%nstates > MAXS )then
                write(*,'(a)') 'ERROR, the number of states is limited to 20!'
                write(*,'(a)') 'In: constructor, module: simple_params.f90'
                stop
            endif
        endif
        ! determines whether at least one volume is on the cmdline
        do i=1,self%nstates
            if( cline%defined( trim('vol')//int2str(i) ))vol_defined(i) = .true.
        enddo
        ! check inputted vols
        if( cline%defined('vollist') )then
            if( nlines(self%vollist)< MAXS )then
                call read_vols
            endif
        else
            if( cline%defined('vol') )then
                self%vols(1) = self%vol
            endif
            if( cline%defined('vol') .or. any( vol_defined) )then
                do i=1,MAXS
                    call check_vol( i )
                end do
            endif
        endif
        if( self%stk .eq. '' .and. self%vols(1) .ne. '' )then
            call find_ldim_nptcls(self%vols(1), self%ldim, ifoo, endconv=conv)
            self%box  = self%ldim(1)
            if( debug ) print *, 'found logical dimension of volume: ', self%ldim
        endif
        if( self%stk .ne. '' )then
            inquire(FILE=self%stk, EXIST=here)
            if( here )then
                if( cline%defined('box') )then
                else
                    call find_ldim_nptcls(self%stk, self%ldim, ifoo, endconv=conv)
                    self%ldim(3) = 1
                    if( debug ) print *, 'found logical dimension of stack: ', self%ldim
                    self%box     = self%ldim(1)
                endif
                if( .not. cline%defined('nptcls') )then
                    ! get number of particles from stack
                     call find_ldim_nptcls(self%stk, self%ldim, self%nptcls, endconv=conv)
                     if( debug ) print *, 'found logical dimension of stack: ', self%ldim
                     if( debug ) print *, 'found nr of ptcls from stack: ', self%nptcls
                     self%ldim(3) = 1
                endif
            else
                write(*,'(a,1x,a)') 'Inputted stack file does not exist!', trim(self%stk)
                stop
            endif
        else if( self%oritab .ne. '' )then
            ! get nr of particles from oritab
            if( .not. cline%defined('nptcls') ) self%nptcls = nlines(self%oritab)
        else if( self%refs .ne. '' )then
            inquire(FILE=self%refs, EXIST=here)
            if( here )then
                if( cline%defined('box') )then
                else
                    call find_ldim_nptcls(self%refs, self%ldim, ifoo, endconv=conv)
                    self%ldim(3) = 1
                    if( debug ) print *, 'found logical dimension of refs: ', self%ldim
                    self%box = self%ldim(1)
                endif
            else
                write(*,'(a,1x,a)') 'Inputted stack file does not exist!', trim(self%refs)
                stop
            endif
        endif
        ! check parallelisation mode (even/chunk)
        nparts_set = .false.
        if( self%split_mode .eq. 'chunk' )then
            if( cline%defined('ncls') .and. cline%defined('nspace') )then
                write(*,*) 'ncls (nr of classes) and nspace (nr of references) cannot'&
                &' simlutaneously be part of the command line!'
                stop 'simple_params :: new'
            else if( cline%defined('ncls') )then
                self%chunksz = self%ncls
            else if( cline%defined('nspace') )then
                self%chunksz = self%nspace
            else 
                write(*,*) 'either ncls (nr of classes) and nspace (nr of references) need'&
                &' to be part of the command line for this mode of job distribution'
                stop 'simple_params :: new'
            endif
            parts       = split_nobjs_in_chunks(self%nptcls, self%chunksz)
            self%nparts = size(parts,1)
            nparts_set  = .true.
            deallocate(parts)
        else
            self%split_mode = 'even'
            nparts_set      = .false.
            if( cline%defined('nparts') ) nparts_set  = .true.
        endif
        if( .not. cline%defined('ncunits') )then
            ! we assume that the number of computing units is equal to the number of partitions
            self%ncunits = self%nparts
        endif
        ! check file formats
        call check_file_formats(aamix)
        call double_check_file_formats
        ! make file names
        call mkfnames
        ! check box
        if( self%box > 0 .and. self%box < 26 ) stop 'box size need to be larger than 26; simple_params'
        if( cline%defined('nthr') )then
!$          call omp_set_num_threads(self%nthr)       
        else
!$          self%nthr = omp_get_max_threads()
!$          call omp_set_num_threads(self%nthr)   
        endif
        nthr_glob = self%nthr
        if( .not. cline%defined('xdim') ) self%xdim = self%box/2
        self%xdimpd = round2even(self%alpha*real(self%box/2))
        self%boxpd  = 2*self%xdimpd
        ! set derived Fourier related variables
        self%dstep = real(self%box-1)*self%smpd                   ! first wavelength of FT
        self%dsteppd = real(self%boxpd-1)*self%smpd               ! first wavelength of padded FT
        if( .not. cline%defined('hp') ) self%hp = 0.7*self%dstep  ! high-pass limit
        self%fny = 2.*self%smpd                                   ! Nyqvist limit
        if( .not. cline%defined('lpstop') )then                   ! default highest resolution lp
            self%lpstop = self%fny                                ! deafult lpstop
        endif
        if( self%fny > 0. ) self%tofny = int(self%dstep/self%fny) ! Nyqvist Fourier index
        if( cline%defined('lp') )then                             ! override dynlp=yes and lpstop
            self%dynlp = 'no'
            ! TAKEN OUT BECAUSE THE PREPROC IS UNHAPPY
            ! if( self%lp < self%lpstop )  self%lpstop  = self%lp
            ! if( self%lpstart > self%lp ) self%lpstart = self%lp
        endif
        ! set default ring2 value
        if( .not. cline%defined('ring2') )then
            if( cline%defined('msk') )then
                self%ring2 = round2even(self%msk)
            else
                self%ring2 = round2even(real(self%box/2-2))
            endif
        else
            self%ring2 = round2even(real(min(self%box/2-2,self%ring2)))
        endif
        ! set default msk value
        if( .not. cline%defined('msk') )then
            if( cline%defined('ring2') )then
                self%msk = self%ring2
            else
                self%msk = self%box/2
            endif
        endif
        ! set mode of masking
        if( cline%defined('inner') )then
            self%l_innermsk = .true.
        else
            self%l_innermsk = .false.
        endif
        ! set nr of rotations
        self%nrots = round2even(twopi*real(self%ring2))
        ! set boxmatch (for clipped imags) and xfel logical
        self%l_xfel = .false.
        if( self%xfel .eq. 'yes')then
            self%boxmatch = self%box
            self%l_xfel    = .true.
        else if( self%box > 2*(self%msk+10) )then
            self%boxmatch = find_magic_box(2*(nint(self%msk)+10))
            if( self%boxmatch > self%box ) self%boxmatch = self%box
        else
            self%boxmatch = self%box
        endif
        ! set default outer mask value
        if( .not. cline%defined('outer') ) self%outer = self%msk
        ! checks automask related values
        self%l_automsk = .false.
        self%doautomsk = .false.
        if( self%automsk.eq.'yes' )                          self%l_automsk = .true.
        if( cline%defined('mw') .and. self%automsk.ne.'no' ) self%doautomsk = .true.
        if( .not.cline%defined('mw') .and. self%automsk.eq.'yes') &
            write(*,*) 'WARNING! MW argument not provided in conjunction with AUTOMSK'
        if( self%doautomsk )then
            if( self%edge <= 0    ) stop 'Invalid value for edge' 
            if( self%binwidth < 0 ) stop 'Invalid value for binwidth'
        endif
        ! set newbox if scale is defined
        if( cline%defined('scale') )then
            self%newbox = find_magic_box(nint(self%scale*real(self%box)))
        endif
         ! set newbox if scale is defined
        if( cline%defined('scale2') )then
            self%newbox2 = find_magic_box(nint(self%scale2*real(self%box)))
        endif
        self%kfromto(1) = max(2,int(self%dstep/self%hp)) ! high-pass Fourier index set according to hp
        self%kfromto(2) = int(self%dstep/self%lp)        ! low-pass Fourier index set according to lp
        self%lp         = max(self%fny,self%lp)          ! lowpass limit
        self%lp_dyn     = self%lp                        ! dynamic lowpass limit
        self%lpmed      = self%lp                        ! median lp
        if( .not. cline%defined('ydim') ) self%ydim = self%xdim
        ! set ldim
        if( cline%defined('xdim') ) self%ldim = [self%xdim,self%ydim,1]
        ! if CTF refinement
        if( self%ctf .eq. 'refine' ) stop 'CTF refinement not yet implemented'
        ! set eullims
        self%eullims(:,1) = 0.
        self%eullims(:,2) = 359.99
        self%eullims(2,2) = 180.
        ! set default size of random sample
        if( .not. cline%defined('nran') )then
            self%nran = self%nptcls
        endif
        ! define clipped box if not given
        if( .not. cline%defined('clip') )then
            self%clip = self%box
        endif
        ! error check ncls
        if( file_exists(self%refs) )then
            ! get number of particles from stack
            call find_ldim_nptcls(self%refs, lfoo, ncls, endconv=conv)
            if( debug ) print *, 'found ncls from refs: ', ncls
            if( cline%defined('ncls') )then
                if( ncls /= self%ncls ) stop 'inputtend number of clusters (ncls) not&
                &consistent with the number of references in stack (p%refs)'
            else
                self%ncls = ncls
            endif
        endif
        ! set endconv it simple_defs
        if( allocated(endconv) ) deallocate(endconv)
        if( cline%defined('endian') )then
            select case(self%endian)
                case('big')
                    allocate(endconv, source='BIG_ENDIAN')
                case('little')
                    allocate(endconv, source='LITTLE_ENDIAN')
                case('native')
                    allocate(endconv, source='NATIVE')
                case DEFAULT
                    stop 'unsupported endianness flag; simple_params :: constructor'
            end select
        else if( allocated(conv) )then
            allocate(endconv, source=conv)
        else
            allocate(endconv, source='NATIVE')
        endif
        if( allocated(conv) ) deallocate(conv)        
        ! set to particle index if not defined in cmdlin
        if( .not. cline%defined('top') ) self%top = self%nptcls
        ! set the number of input orientations
        if( .not. cline%defined('noris') ) self%noris = self%nptcls
        ! fix translation param
        self%trs = abs(self%trs)
        if( .not. cline%defined('trs') )then
            if( self%refine.eq.'no' .or. self%refine.eq.'shc' .or. self%refine.eq.'qcont' )then
                self%trs = 0.
            else
                self%trs = 1.
            endif
        endif
        self%doshift = .true.
        if( self%trs < 0.5 )then
            self%trs = 0.00001
            self%doshift = .false.
        endif
        ! set optlims
        self%optlims(:3,:)  = self%eullims
        self%optlims(4:5,1) = -self%trs
        self%optlims(4:5,2) = self%trs
        self%optlims(6:7,1) = -self%dfsdev
        self%optlims(6:7,2) = self%dfsdev
        ! Set first three in cyclic true
        self%cyclic(1) = .true.
        self%cyclic(2) = .true.
        self%cyclic(3) = .true.
        ! Set molecular diameter
        if( .not. cline%defined('moldiam') )then
            self%moldiam = 0.7*real(self%box)*self%smpd
        endif
        ! set execution mode (shmem or distr)
        if( cline%defined('part') )then
            self%l_distr_exec = .true.
        else
            self%l_distr_exec = .false.
        endif
        ! set global distributed execution flag
        l_distr_exec_glob = self%l_distr_exec
        ! set lenght of number string for zero padding
        if( .not. cline%defined('numlen') )then
            if( nparts_set ) self%numlen = len(int2str(self%nparts))
        endif
        ! set name of partial stack in parallel execution
        if( self%numlen > 0 )then
            self%stk_part = 'stack_part'//int2str_pad(self%part,self%numlen)//self%ext
        else
            self%stk_part = 'stack_part'//int2str(self%part)//self%ext
        endif
        ! Check for the existance of this file if part is defined on the command line
        if( cline%defined('part') )then
            if( ccheckdistr )then
                if( .not. file_exists(self%stk_part) )then
                    write(*,*) 'Need partial stacks to be generated for parallel execution'
                    write(*,*) 'Use simple_exec prg=split'
                    stop
                endif
            endif
        endif
        ! set imgkind and check so that ctf and xfel parameters are congruent
        self%imgkind = 'em'
        if( self%xfel .eq. 'yes' ) then
            self%imgkind = 'xfel'
            if( self%ctf .eq. 'no') then
                ! all good
            else
                stop 'when xfel .eq. yes, ctf .ne. no is not allowed; simple_params'
            endif
        endif
        ! create are/volume ratios object for SSNR re-weighting
        if( cline%defined('mw') )then
            self%avr = AVratios([self%box,self%box,self%box], [self%box,self%box,1], self%msk, self%smpd, self%mw)
        else
            self%avr = AVratios([self%box,self%box,self%box], [self%box,self%box,1], self%msk, self%smpd)
        endif
        ! take care of nspace value for refine .eq. 'qcontneigh' mode
        if( self%refine .eq. 'qcontneigh' ) self%nspace = self%nnn
        ! check if we are dose-weighting or not
        self%l_dose_weight = .false.
        if( cline%defined('exp_time') .or. cline%defined('dose_rate') )then
            if( cline%defined('exp_time') .and. .not. cline%defined('dose_rate') )then
                stop 'need dose_rate to be part of the command line! simple_params :: new'
            endif
            if( .not. cline%defined('exp_time') .and. cline%defined('dose_rate') )then
                stop 'need exp_time to be part of the command line! simple_params :: new'
            endif
            self%l_dose_weight = .true.
        endif
        ! prepare CTF plan
        self%tfplan%flag = self%ctf
        ! set logical shellw flag
        self%l_shellw = .true.
        if( cline%defined('shellw') )then
            if( self%shellw .eq. 'no' ) self%l_shellw = .false.
        endif
        ! set logical pick flag
        self%l_pick = .false.
        if( self%dopick .eq. 'yes' ) self%l_pick = .true.

        write(*,'(A)') '>>> DONE PROCESSING PARAMETERS'

      contains
          
          subroutine check_vol( i )
              integer, intent(in) :: i
              character(len=STDLEN) :: nam
              logical :: here
              nam = 'vol'//int2str(i)
              if( cline%defined(nam) )then
                  call check_file(nam, self%vols(i), notAllowed='T')
                  inquire(FILE=self%vols(i), EXIST=here)
                  if( .not. here )then
                      write(*,*) 'Input volume:', self%vols(i), 'does not exist!'
                      stop
                  endif
                  !self%nstates = i
                  if( self%debug == 'yes' ) write(*,*) nam, '=', self%vols(i)
              endif
          end subroutine check_vol

          subroutine read_vols
              character(len=STDLEN) :: filenam, nam
              integer :: nl, fnr, i
              filenam = cline%get_carg('vollist')
              nl = nlines(filenam)
              fnr = get_fileunit( )
              open(fnr, file=filenam)
              do i=1,nl
                  read(fnr,*) nam
                  if( nam .ne. '' )then
                      self%vols(i) = nam
                      !self%nstates = i
                  endif
              end do
              close(fnr)
          end subroutine read_vols

          subroutine check_file( file, var, allowed1, allowed2, notAllowed )
              character(len=*),           intent(in)  :: file
              character(len=*),           intent(out) :: var
              character(len=1), optional, intent(in)  :: allowed1, allowed2, notAllowed
              character(len=1) :: file_descr
              logical          :: raise_exception
              if( cline%defined(file) )then
                  var = cline%get_carg(file)
                  if( debug ) write(*,*) 'var = ', var
                  file_descr = fname2format(var)
                  if( debug ) write(*,*) 'file_descr = ', file_descr
                  raise_exception = .false.
                  if( present(allowed1) )then
                      if( allowed1 == file_descr ) then
                          if( present(allowed2) )then
                              if( allowed2 == file_descr ) then
                                  ! all good
                              else
                                  raise_exception = .true.
                              endif
                          endif
                      else
                          raise_exception = .true.
                      endif
                  endif
                  if( present(notAllowed) )then
                      if( notAllowed == file_descr ) raise_exception = .true.
                  endif
                  if( raise_exception )then
                      write(*,*) 'This format: ', file_descr, ' is not allowed for this file: ', var
                      stop
                  endif
                  select case(file_descr)
                      case ('I')
                          stop 'Support for IMAGIC files is not yet implemented!'
                      case ('M')
                          ! MRC files are supported
                          cntfile = cntfile+1
                          checkupfile(cntfile) = 'M'
                      case ('S')
                          ! SPIDER files are supported
                          cntfile = cntfile+1
                          checkupfile(cntfile) = 'S'
                      case ('N')
                          stop 'This file format is not supported by SIMPLE; simple_params::check_file'
                      case ('T')
                          ! text files are supported
                      case ('B')
                          ! binary files are supported
                      case DEFAULT
                          stop 'This file format is not supported by SIMPLE; simple_params::check_file'
                  end select
                  if( debug_local == 'yes' ) write(*,*) file, '=', var
              endif
          end subroutine check_file

          subroutine check_file_formats( allow_mix )
              logical, intent(in) :: allow_mix
              integer :: i, j
              if( cntfile > 0 )then
                  do i=1,cntfile
                      do j=1,cntfile
                          if( i == j ) cycle
                          if( checkupfile(i) == checkupfile(j) )then
                              ! all ok
                          else
                              if( .not. allow_mix )then
                                  stop 'The inputted file names have nonconforming format (mixed formats not yet allowed)!'
                              endif
                          endif
                      end do
                  end do
                  select case(checkupfile(1))
                      case('M')
                          self%ext = '.mrc'
                      case('S')
                          self%ext = '.spi'
                      case('D')
                          self%ext = '.mrc'
                      case('B')
                          self%ext = '.mrc'
                      case DEFAULT
                          stop 'This file format is not supported by SIMPLE; check_file_formats; simple_params'
                  end select
              endif
          end subroutine check_file_formats

          subroutine double_check_file_formats
              character(len=STDLEN) :: fname
              character(len=1)      :: form
              integer :: funit
              if( cntfile == 0 )then
                  if( cline%defined('filetab') )then
                      funit = get_fileunit()
                      open(unit=funit, status='old', file=self%filetab)
                      read(funit,'(a256)') fname
                      form = fname2format(fname)
                      close(funit)
                      select case(form)
                          case('M')
                              self%ext = '.mrc'
                          case('S')
                              self%ext = '.spi'
                          case('D')
                              self%ext = '.mrc'
                          case('B')
                              self%ext = '.mrc'
                          case DEFAULT
                              write(*,*) 'format string is ', form
                              stop 'This file format is not supported by SIMPLE; double_check_file_formats; simple_params'
                      end select
                  endif
              endif
          end subroutine double_check_file_formats

          subroutine mkfnames
              integer :: i
              do i=1,self%nstates
                  self%vols_msk(i) = add2fbody(self%vols(i), self%ext, 'msk')
                  self%masks(i)    = 'automask_state'//int2str_pad(i,2)//self%ext
              end do
              if( .not. cline%defined('outstk')  ) self%outstk  = 'outstk'//self%ext
              if( .not. cline%defined('outstk2') ) self%outstk2 = 'outstk2'//self%ext
              if( .not. cline%defined('outvol')  ) self%outvol  = 'outvol'//self%ext
          end subroutine mkfnames

          subroutine check_carg( carg, var )
              character(len=*), intent(in)  :: carg
              character(len=*), intent(out) :: var
              if( cline%defined(carg) )then
                  var = cline%get_carg(carg)
                  if( debug_local == 'yes' ) write(*,*) carg, '=', var
              endif
          end subroutine check_carg

          subroutine check_iarg( iarg, var )
              character(len=*), intent(in)  :: iarg
              integer, intent(out) :: var
              if( cline%defined(iarg) )then
                  var = nint(cline%get_rarg(iarg))
                  if( debug_local == 'yes' ) write(*,*) iarg, '=', var
              endif
          end subroutine check_iarg

          subroutine check_rarg( rarg, var )
              character(len=*), intent(in)  :: rarg
              real, intent(out) :: var
              if( cline%defined(rarg) )then
                  var = cline%get_rarg(rarg)
                  if( debug_local == 'yes' ) write(*,*) rarg, '=', var
              endif
          end subroutine check_rarg
        
    end subroutine new

end module simple_params
