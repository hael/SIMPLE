!@descr: constructor and semantic execution, source, derivation, and validation phases for SIMPLE parameters
submodule(simple_parameters) simple_parameters_phases
use simple_sp_project, only: sp_project
implicit none
#include "simple_local_flags.inc"

contains

    module subroutine new(self, cline, silent)
        class(parameters), intent(inout) :: self
        class(cmdline),    intent(inout) :: cline
        logical, optional, intent(in)    :: silent
        character(len=1) :: checkupfile(50)
        integer          :: cntfile
        logical          :: ssilent
        select type(self)
            type is(parameters)
                self = parameters(box_extract=0, hpind_fsc=0, kfromto=[0,0], lplims2D=[0.,0.,0.], &
                    &smpd_pickrefs=0., smpd_targets2D=[0.,0.], total_dose=0.)
            class default
                THROW_HARD('Unsupported parameters dynamic type; simple_parameters_ctor')
        end select
        call self%init_dynamic_defaults()
        ssilent = .false.
        if( present(silent) ) ssilent = silent
        call seed_rnd
        cntfile = 0
        call self%parse_inputs(cline, cntfile, checkupfile)
        call self%setup_execution_context(cline, ssilent)
        call self%resolve_parameter_sources(cline, cntfile, checkupfile)
        call self%derive_sampling_settings(cline)
        call self%derive_parallel_settings(cline)
        call self%derive_image_settings(cline)
        call self%validate_parameter_consistency(cline)
        if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> DONE PROCESSING PARAMETERS'
    end subroutine new

    module subroutine setup_execution_context(self, cline, silent)
        class(parameters), intent(inout) :: self
        class(cmdline),    intent(inout) :: cline
        logical,           intent(in)    :: silent
        type(string), allocatable :: sp_files(:)
        character(len=STDLEN)     :: str_static
        type(string)              :: str_tmp, absname
        integer :: i, idir, nsp_files
        call simple_getcwd(self%cwd)
        if( .not. allocated(CWD_GLOB_ORIG) ) allocate(CWD_GLOB_ORIG, source=self%cwd%to_char())
        CWD_GLOB = self%cwd%to_char()
        self%pid = get_process_id()
        call get_command_argument(0, str_static)
        self%executable = trim(adjustl(str_static))
        if( self%executable%strlen_trim() == 0 )then
            THROW_HARD('get_command_argument failed; setup_execution_context')
        endif
        if( cline%defined('part') )then
            self%l_distr_worker = .true.
            part_glob = self%part
        else
            self%l_distr_worker = .false.
            part_glob = 1
        endif
        l_distr_worker_glob = self%l_distr_worker
        call get_prg_ptr(self%prg, self%ptr2prg)
        call self%last_prev_dir%kill
        idir = find_next_int_dir_prefix(self%cwd, self%last_prev_dir)
        if( .not. cline%defined('projfile') )then
            if( self%last_prev_dir%is_allocated() )then
                call simple_list_files(self%last_prev_dir%to_char()//'/*.simple', sp_files)
            else
                call simple_list_files('*.simple', sp_files)
            endif
            if( allocated(sp_files) )then
                nsp_files = size(sp_files)
            else
                nsp_files = 0
            endif
        else
            nsp_files = 0
        endif
        if( associated(self%ptr2prg) .and. .not. self%executable%has_substr('private') )then
            self%sp_required = self%ptr2prg%requires_sp_project()
            if( trim(self%prg%to_char()) == 'model_cavgs_rejection' .and. &
                (trim(self%quality_mode) == 'learn' .or. trim(self%quality_mode) == 'promote' .or. &
                 (trim(self%quality_mode) == 'evaluate' .and. cline%defined('filetab'))) ) &
                self%sp_required = .false.
            if( .not. cline%defined('projfile') .and. self%sp_required )then
                if( nsp_files > 1 )then
                    write(logfhandle,*) 'Multiple *simple project files detected'
                    do i=1,nsp_files
                        write(logfhandle,*) sp_files(i)%to_char()
                    end do
                    THROW_HARD('a unique *.simple project could NOT be identified; setup_execution_context')
                endif
            endif
        else
            self%sp_required = .false.
        endif
        if( nsp_files == 0 .and. self%sp_required .and. .not. cline%defined('projfile') )then
            write(logfhandle,*) 'program: ', self%prg%to_char(), ' requires a project file!'
            write(logfhandle,*) 'cwd:     ', self%cwd%to_char()
            THROW_HARD('no *.simple project file identified; setup_execution_context')
        endif
        if( nsp_files == 1 .and. self%sp_required )then
            if( .not. cline%defined('projfile') ) self%projfile = sp_files(1)%to_char()
            str_tmp       = get_fbody(basename(self%projfile), string('simple'))
            self%projname = str_tmp%to_char()
            if( .not. cline%defined('projfile') ) call cline%set('projfile', self%projname//'.simple')
            if( .not. cline%defined('projname') ) call cline%set('projname', self%projname)
        endif
        if( self%projfile .ne. '' )then
            if( file_exists(self%projfile) )then
                absname       = simple_abspath(self%projfile)
                self%projfile = absname%to_char()
                str_tmp       = get_fbody(basename(absname), string('simple'))
                self%projname = str_tmp%to_char()
            endif
        endif
        if( self%mkdir .eq. 'yes' )then
            if( associated(self%ptr2prg) .and. .not. self%executable%has_substr('simple_private') )then
                if( self%prg .eq. 'mkdir' )then
                    self%exec_dir = int2str(idir)//'_'//self%dir%to_char()
                else if( cline%defined('dir_exec') )then
                    if( self%executable%has_substr('simple_stream') )then
                        self%exec_dir = self%dir_exec
                    else
                        self%exec_dir = int2str(idir)//'_'//self%dir_exec%to_char()
                    endif
                else
                    if( self%outdir .eq. '' )then
                        self%exec_dir = int2str(idir)//'_'//self%prg%to_char()
                    else
                        self%exec_dir = self%outdir
                    endif
                endif
                call simple_mkdir(filepath(string(PATH_HERE), self%exec_dir))
                if( self%prg .eq. 'mkdir' ) return
                if( .not. silent ) write(logfhandle,'(a)') '>>> EXECUTION DIRECTORY: '//self%exec_dir%to_char()
                call simple_chdir(filepath(string(PATH_HERE), self%exec_dir))
                if( self%sp_required )then
                    call simple_copy_file(self%projfile, filepath(string(PATH_HERE), basename(self%projfile)))
                    self%projfile = filepath(string(PATH_HERE), basename(self%projfile))
                    absname       = simple_abspath(self%projfile)
                    self%projfile = absname
                    str_tmp       = get_fbody(basename(self%projfile), string('simple'))
                    self%projname = str_tmp%to_char()
                    call cline%set('projfile', self%projfile)
                endif
                call simple_getcwd(self%cwd)
                CWD_GLOB = self%cwd%to_char()
            endif
        endif
        if( .not. self%executable%has_substr('private') .and. STDOUT2LOG )then
            call fopen(logfhandle, status='replace', file=string(LOGFNAME), action='write')
        endif
    end subroutine setup_execution_context

    module subroutine resolve_parameter_sources(self, cline, cntfile, checkupfile)
        class(parameters), intent(inout) :: self
        class(cmdline),    intent(inout) :: cline
        integer,           intent(in)    :: cntfile
        character(len=1),  intent(in)    :: checkupfile(:)
        logical                   :: vol_defined(MAXS)
        type(binoris)             :: bos
        type(sp_project)          :: spproj
        type(ori)                 :: o
        type(string)              :: phaseplate, ctfflag
        real                      :: smpd
        integer                   :: i, ifoo, lfoo(3), box, nptcls, istate, cntfile_loc
        logical                   :: def_vol1, def_even, def_odd
        character(len=1)          :: checkupfile_loc(size(checkupfile))
        cntfile_loc = cntfile
        checkupfile_loc = checkupfile
        vol_defined = .false.
        do i=1,size(vol_defined)
            if( cline%defined(trim('vol')//int2str(i)) ) vol_defined(i) = .true.
        end do
        if( any(vol_defined) )then
            do i=1,size(vol_defined)
                if( vol_defined(i) ) call check_vol(self%vols(i))
            end do
        endif
        if( cline%defined('msklist') )then
            if( nlines(self%msklist) < MAXS ) call read_masks
        endif
        if( self%stk .eq. '' .and. cline%defined('vol_even') )then
            call find_ldim_nptcls(self%vol_even, self%ldim, ifoo)
            if( .not. cline%defined('box')      ) self%box      = self%ldim(1)
            if( .not. cline%defined('box_crop') ) self%box_crop = self%box
        else if( self%stk .eq. '' .and. vol_defined(1) )then
            call find_ldim_nptcls(self%vols(1), self%ldim, ifoo)
            if( .not. cline%defined('box')      ) self%box      = self%ldim(1)
            if( .not. cline%defined('box_crop') ) self%box_crop = self%box
        endif
        if( self%mkdir.eq.'yes' )then
            if( self%dir_meta%to_char([1,1])   .ne. PATH_SEPARATOR ) self%dir_meta   = PATH_PARENT//self%dir_meta%to_char()
            if( self%dir_movies%to_char([1,1]) .ne. PATH_SEPARATOR ) self%dir_movies = PATH_PARENT//self%dir_movies%to_char()
            if( self%dir_target%to_char([1,1]) .ne. PATH_SEPARATOR ) self%dir_target = PATH_PARENT//self%dir_target%to_char()
        endif
        if( cline%defined('oritype') )then
            select case(trim(self%oritype))
                case('mic')
                    self%spproj_iseg = MIC_SEG
                case('stk')
                    self%spproj_iseg = STK_SEG
                case('ptcl2D')
                    self%spproj_iseg = PTCL2D_SEG
                case('cls2D')
                    self%spproj_iseg = CLS2D_SEG
                case('cls3D')
                    self%spproj_iseg = CLS3D_SEG
                case('ptcl3D')
                    self%spproj_iseg = PTCL3D_SEG
                case('out')
                    self%spproj_iseg = OUT_SEG
                case('optics')
                    self%spproj_iseg = OPTICS_SEG
                case('projinfo')
                    self%spproj_iseg = PROJINFO_SEG
                case('jobproc')
                    self%spproj_iseg = JOBPROC_SEG
                case('compenv')
                    self%spproj_iseg = COMPENV_SEG
                case DEFAULT
                    write(logfhandle,*) 'oritype: ', trim(self%oritype)
                    THROW_HARD('unsupported oritype; resolve_parameter_sources')
            end select
        else
            self%oritype     = 'ptcl3D'
            self%spproj_iseg = PTCL3D_SEG
        endif
        if( file_exists(self%projfile) )then
            if( self%stream.eq.'no' )then
                if( self%spproj_iseg==OUT_SEG .or. self%spproj_iseg==CLS3D_SEG )then
                    call spproj%read_segment('out', self%projfile)
                    call spproj%get_imginfo_from_osout(smpd, box, nptcls)
                    call spproj%kill
                    if( .not. cline%defined('smpd') )   self%smpd   = smpd
                    if( .not. cline%defined('box') )    self%box    = box
                    if( .not. cline%defined('nptcls') ) self%nptcls = nptcls
                else
                    call bos%open(self%projfile)
                    if( .not. cline%defined('nptcls') ) self%nptcls = bos%get_n_records(self%spproj_iseg)
                    call o%new(is_ptcl=.false.)
                    select case(self%spproj_iseg)
                        case(MIC_SEG)
                            call bos%read_first_segment_record(MIC_SEG, o)
                        case DEFAULT
                            call bos%read_first_segment_record(STK_SEG, o)
                    end select
                    if( o%isthere('smpd') .and. .not. cline%defined('smpd') ) self%smpd = o%get('smpd')
                    if( o%isthere('box')  .and. .not. cline%defined('box')  ) self%box  = nint(o%get('box'))
                    call o%kill
                endif
            endif
            if( .not. bos%is_opened() ) call bos%open(self%projfile)
            select case(trim(self%oritype))
                case('ptcl2D', 'ptcl3D')
                    call bos%read_first_segment_record(STK_SEG, o)
            end select
            if( .not. cline%defined('ctf') )then
                if( o%exists() )then
                    if( o%isthere('ctf') )then
                        call o%getter('ctf', ctfflag)
                        self%ctf = ctfflag%to_char()
                    endif
                endif
            endif
            if( .not. cline%defined('phaseplate') )then
                if( o%isthere('phaseplate') )then
                    call o%getter('phaseplate', phaseplate)
                    self%phaseplate   = phaseplate%to_char()
                    self%l_phaseplate = self%phaseplate.eq.'yes'
                else
                    self%phaseplate   = 'no'
                    self%l_phaseplate = .false.
                endif
            else
                self%l_phaseplate = self%phaseplate.eq.'yes'
            endif
            call bos%close
        else if( self%stk .ne. '' )then
            if( file_exists(self%stk) )then
                if( .not. cline%defined('nptcls') ) call find_ldim_nptcls(self%stk, lfoo, self%nptcls)
            else
                write(logfhandle,'(a,1x,a)') 'Inputted stack (stk) file does not exist!', self%stk%to_char()
                stop
            endif
        else if( self%oritab .ne. '' )then
            if( .not. cline%defined('nptcls') .and. self%oritab%has_substr('.txt') ) self%nptcls = nlines(self%oritab)
        else if( self%refs .ne. '' )then
            if( file_exists(self%refs) )then
                if( .not. cline%defined('box') )then
                    call find_ldim_nptcls(self%refs, self%ldim, ifoo)
                    self%ldim(3) = 1
                    self%box = self%ldim(1)
                endif
                if( .not. cline%defined('box_crop') )then
                    call find_ldim_nptcls(self%refs, self%ldim, ifoo)
                    self%ldim(3) = 1
                    self%box_crop = self%ldim(1)
                endif
            endif
        endif
        call set_ldim_box_from_stk
        self%scale_movies = 1.0
        if( trim(self%downscale).eq.'yes' )then
            self%scale_movies = self%smpd / self%smpd_downscale
            if( self%scale_movies > 1.0 )then
                self%smpd_downscale = self%smpd
                self%scale_movies   = 1.0
            endif
        endif
        if( .not. cline%defined('box_crop') ) self%box_crop = self%box
        if( .not. cline%defined('smpd_crop') )then
            if( cline%defined('box_crop') )then
                self%smpd_crop = real(self%box)/real(self%box_crop) * self%smpd
            else
                self%smpd_crop = self%smpd
            endif
        endif
        call check_file_formats
        call double_check_file_formats
        call mkfnames
        if( self%box > 0 .and. self%box < 26 )then
            THROW_HARD('box needs to be larger than 26')
        endif
        if( self%box_crop > 0 .and. self%box_crop < 26 )then
            THROW_HARD('box_crop needs to be larger than 26')
        endif
        if( cline%defined('refs') )then
            self%refs_even = add2fbody(self%refs, self%ext, '_even')
            self%refs_odd  = add2fbody(self%refs, self%ext, '_odd')
        endif
        def_vol1 = cline%defined('vol1')
        def_even = cline%defined('vol_even')
        def_odd  = cline%defined('vol_odd')
        if( def_vol1 )then
            if( def_even )then
                THROW_HARD('vol1 and vol_even cannot be simultaneously part of the command line')
            endif
            if( def_odd )then
                THROW_HARD('vol1 and vol_odd cannot be simultaneously part of the command line')
            endif
        else
            if( def_even .and. def_odd )then
                if( def_vol1 )then
                    THROW_HARD('vol1 cannot be part of command line when vol_even and vol_odd are')
                endif
            else
                if( (def_even .eqv. .false.) .and. (def_odd .eqv. .false.) )then
                else
                    THROW_HARD('vol_even and vol_odd must be simultaneously part of the command line')
                endif
            endif
        endif
        if( def_even .or. def_odd .or. def_vol1 )then
            if( def_even .and. def_odd )then
                if( def_vol1 )then
                    THROW_HARD('vol1 cannot be part of command line when vol_even and vol_odd are')
                endif
                self%vols_even(1) = self%vol_even
                self%vols_odd(1)  = self%vol_odd
                self%vols(1)      = self%vol_even
            else
                if( def_vol1 )then
                    do istate=1,self%nstates
                        self%vols_even(istate) = add2fbody(self%vols(istate), self%ext, '_even')
                        self%vols_odd(istate)  = add2fbody(self%vols(istate), self%ext, '_odd')
                        if( self%prg%has_substr('refine') )then
                            if( .not. file_exists(self%vols_even(istate)) ) call simple_copy_file(self%vols(istate), self%vols_even(istate))
                            if( .not. file_exists(self%vols_odd(istate))  ) call simple_copy_file(self%vols(istate), self%vols_odd(istate))
                        endif
                    end do
                else
                    THROW_HARD('both vol_even and vol_odd need to be part of the command line')
                endif
            endif
        endif

    contains

        subroutine check_vol(volname, is_even)
            class(string),     intent(inout) :: volname
            logical, optional, intent(in)    :: is_even
            type(string) :: vol, key
            if( present(is_even) )then
                if( is_even )then
                    key = 'vol_even'
                else
                    key = 'vol_odd'
                endif
            else
                key = 'vol'//int2str(i)
            endif
            if( cline%defined(key%to_char()) )then
                vol = cline%get_carg(key%to_char())
                if( vol%to_char([1,1]).eq.PATH_SEPARATOR )then
                    volname = vol
                    if( .not. file_exists(volname) )then
                        write(logfhandle,*) 'Input volume:', volname%to_char(), ' does not exist! 1'
                        stop
                    endif
                else
                    if( self%mkdir .eq. 'yes' ) call cline%set(key%to_char(), PATH_PARENT//vol%to_char())
                    volname = cline%get_carg(key%to_char())
                    if( .not. file_exists(volname) )then
                        write(logfhandle,*) 'Input volume:', volname%to_char(), ' does not exist! 2'
                        stop
                    else
                        call cline%set(key%to_char(), simple_abspath(volname, check_exists=.false.))
                    endif
                endif
            endif
        end subroutine check_vol

        subroutine read_masks
            type(string) :: filename, name
            integer      :: nl, fnr, i, io_stat, ios
            filename = cline%get_carg('msklist')
            if( filename%to_char([1,1]).ne.PATH_SEPARATOR )then
                if( self%mkdir.eq.'yes' ) filename = PATH_PARENT//filename%to_char()
            endif
            nl = nlines(filename)
            call fopen(fnr, file=filename, iostat=io_stat)
            if(io_stat /= 0) call fileiochk("parameters ; read_masks error opening "//filename%to_char(), io_stat)
            do i=1,nl
                call name%readline(fnr, ios)
                if( ios /= 0 ) exit
                self%mskvols(i) = trim(name%to_char())
            end do
            call fclose(fnr)
        end subroutine read_masks

        subroutine check_file_formats
            integer :: ii, jj
            if( cntfile_loc > 0 )then
                do ii=1,cntfile_loc
                    do jj=1,cntfile_loc
                        if( ii == jj ) cycle
                        if( checkupfile_loc(ii) == checkupfile_loc(jj) ) cycle
                    end do
                end do
                call self%set_img_format(checkupfile_loc(1))
            endif
        end subroutine check_file_formats

        subroutine double_check_file_formats
            type(string) :: fname
            character(len=1) :: form
            integer :: funit, io_stat, ios
            if( cntfile_loc == 0 )then
                if( cline%defined('filetab') )then
                    if( trim(self%prg%to_char()) == 'model_cavgs_rejection' .and. &
                        (trim(self%quality_mode) == 'learn' .or. trim(self%quality_mode) == 'evaluate') ) return
                    call fopen(funit, status='old', file=self%filetab, iostat=io_stat)
                    call fileiochk("In parameters::double_check_file_formats fopen failed "//self%filetab%to_char(), io_stat)
                    call fname%readline(funit, ios)
                    form = fname2format(fname)
                    call fclose(funit)
                    call self%set_img_format(form)
                endif
            endif
        end subroutine double_check_file_formats

        subroutine mkfnames
            if( .not. cline%defined('outstk') ) self%outstk = 'outstk'//self%ext%to_char()
            if( .not. cline%defined('outvol') ) self%outvol = 'outvol'//self%ext%to_char()
        end subroutine mkfnames

        subroutine set_ldim_box_from_stk
            integer :: ifoo
            if( cline%defined('stk') )then
                if( file_exists(self%stk) )then
                    if( .not. cline%defined('box') )then
                        call find_ldim_nptcls(self%stk, self%ldim, ifoo)
                        self%ldim(3) = 1
                        self%box     = self%ldim(1)
                        if( .not. cline%defined('box_crop') ) self%box_crop = self%box
                    endif
                else
                    write(logfhandle,'(a)')      'simple_parameters :: set_ldim_box_from_stk'
                    write(logfhandle,'(a,1x,a)') 'Stack file does not exist!', self%stk%to_char()
                    THROW_HARD("set_ldim_box_from_stk")
                endif
            endif
        end subroutine set_ldim_box_from_stk

    end subroutine resolve_parameter_sources

    module subroutine derive_sampling_settings(self, cline)
        class(parameters), intent(inout) :: self
        class(cmdline),    intent(inout) :: cline
        if( self%update_frac <= .99 )then
            self%l_update_frac = .true.
            self%l_trail_rec   = trim(self%trail_rec).eq.'yes'
            self%l_fillin      = trim(self%fillin).eq.'yes'
        else
            self%update_frac   = 1.0
            self%l_update_frac = .false.
            self%l_trail_rec   = .false.
            self%l_fillin      = .false.
        endif
        self%l_update_missing = trim(self%update_missing).eq.'yes'
        self%l_frac_best   = self%frac_best  <= 0.99
        self%l_frac_worst  = self%frac_worst <= 0.99
        self%l_greedy_smpl = trim(self%greedy_sampling).eq.'yes'
    end subroutine derive_sampling_settings

    module subroutine derive_parallel_settings(self, cline)
        use simple_gpu_utils, only: set_offload_device
        class(parameters), intent(inout) :: self
        class(cmdline),    intent(inout) :: cline
        logical :: nparts_set
        integer :: nthr
        self%split_mode = 'even'
        nparts_set      = .false.
        if( cline%defined('nparts') ) nparts_set = .true.
        if( .not. cline%defined('numlen') )then
            if( nparts_set ) self%numlen = len(int2str(self%nparts))
        endif
        if( .not. cline%defined('ncunits') ) self%ncunits = self%nparts
        if( cline%defined('nthr') )then
            nthr = cline%get_iarg('nthr')
            if( nthr == 0 )then
                !$ self%nthr = omp_get_max_threads()
                !$ call omp_set_num_threads(self%nthr)
            else
                !$ call omp_set_num_threads(self%nthr)
            endif
        else
            !$ self%nthr = omp_get_max_threads()
            !$ call omp_set_num_threads(self%nthr)
        endif
        nthr_glob = self%nthr
        call set_offload_device(cline, self%device)
    end subroutine derive_parallel_settings

    module subroutine derive_image_settings(self, cline)
        class(parameters), intent(inout) :: self
        class(cmdline),    intent(inout) :: cline
        real    :: mskdiam_default, msk_default
        integer :: lfoo(3), ncls
        if( .not. cline%defined('xdim') ) self%xdim = self%box / 2
        self%xdimpd     = round2even(KBALPHA*real(self%box/2))
        self%boxpd      = 2*self%xdimpd
        self%box_croppd = 2*round2even(KBALPHA*real(self%box_crop/2))
        self%dstep   = real(self%box-1)   * self%smpd
        self%dsteppd = real(self%boxpd-1) * self%smpd
        if( .not. cline%defined('hp') ) self%hp = 0.5*self%dstep
        self%fny = 2.*self%smpd
        if( .not. cline%defined('lpstop') ) self%lpstop = self%fny
        if( self%fny > 0. ) self%tofny = nint(self%dstep/self%fny)
        self%lpstart = max(self%lpstart, self%fny)
        self%lpstop  = max(self%lpstop,  self%fny)
        self%lplims2D(1)       = max(self%fny, self%lpstart)
        self%lplims2D(2)       = max(self%fny, self%lplims2D(1) - (self%lpstart - self%lpstop)/2.)
        self%lplims2D(3)       = max(self%fny, self%lpstop)
        self%smpd_targets2D(1) = min(MAX_SMPD, self%lplims2D(2)*LP2SMPDFAC2D)
        self%smpd_targets2D(2) = min(MAX_SMPD, self%lplims2D(3)*LP2SMPDFAC2D)
        mskdiam_default = (real(self%box) - COSMSKHALFWIDTH) * self%smpd
        msk_default     = round2even((real(self%box) - COSMSKHALFWIDTH) / 2.)
        if( cline%defined('mskdiam') )then
            self%msk = round2even((self%mskdiam / self%smpd) / 2.)
            if( self%msk > msk_default )then
                if( msk_default > 0. )then
                    THROW_WARN('Mask diameter too large, falling back on default value')
                    self%mskdiam = mskdiam_default
                    self%msk     = msk_default
                endif
            endif
            if( self%msk < 0.1 )then
                if( msk_default > 0. )then
                    THROW_WARN('Mask diameter zero, falling back on default value')
                    self%mskdiam = mskdiam_default
                    self%msk     = msk_default
                endif
            endif
        else
            self%mskdiam = mskdiam_default
            self%msk     = msk_default
        endif
        if( self%box > 0 )then
            if( .not. cline%defined('msk_crop') )then
                self%msk_crop = min(round2even((real(self%box_crop)-COSMSKHALFWIDTH)/2.), &
                    &round2even(self%msk * real(self%box_crop) / real(self%box)))
            endif
        endif
        if( self%box > 0 )then
            if( .not. cline%defined('pftsz') ) self%pftsz = magic_pftsz(self%msk_crop)
        endif
        self%l_lpset  = cline%defined('lp')
        self%l_envfsc = self%envfsc .ne. 'no'
        if( self%nstates > 1 .and. self%l_envfsc )then
            THROW_WARN('envfsc disabled for nstates > 1')
            self%envfsc  = 'no'
            self%l_envfsc = .false.
        endif
        if( cline%defined('icm') )    self%l_icm    = (trim(self%icm).eq.'yes')
        if( cline%defined('gauref') ) self%l_gauref = (trim(self%gauref).eq.'yes')
        self%l_corrw = self%wcrit .ne. 'no'
        if( self%l_corrw )then
            select case(trim(self%wcrit))
                case('softmax')
                    self%wcrit_enum = CORRW_CRIT
                case('zscore')
                    self%wcrit_enum = CORRW_ZSCORE_CRIT
                case('sum')
                    self%wcrit_enum = RANK_SUM_CRIT
                case('cen')
                    self%wcrit_enum = RANK_CEN_CRIT
                case('exp')
                    self%wcrit_enum = RANK_EXP_CRIT
                case('inv')
                    self%wcrit_enum = RANK_INV_CRIT
                case('uniform')
                    self%wcrit_enum = UNIFORM_CRIT
                case DEFAULT
                    THROW_HARD('unsupported correlation weighting method')
            end select
        endif
        self%l_graphene       = self%graphene_filt .ne. 'no'
        self%l_nu_refine      = trim(self%nu_refine).eq.'yes'
        self%l_autoscale      = self%autoscale .eq. 'yes'
        if( .not. cline%defined('newbox') )then
            if( cline%defined('scale') ) self%newbox = find_magic_box(nint(self%scale*real(self%box)))
        endif
        self%kfromto = 2
        if( cline%defined('hp') ) self%kfromto(1) = max(2, int(self%dstep/self%hp))
        self%lp         = max(self%fny, self%lp)
        self%kfromto(2) = int(self%dstep/self%lp)
        if( .not. cline%defined('ydim') ) self%ydim = self%xdim
        if( cline%defined('xdim') ) self%ldim = [self%xdim, self%ydim, 1]
        if( cline%defined('box') )  self%ldim = [self%box, self%box, 1]
        self%eullims(:,1) = 0.
        self%eullims(:,2) = 359.99
        self%eullims(2,2) = 180.
        if( .not. cline%defined('nran') ) self%nran = self%nptcls
        if( .not. cline%defined('clip') ) self%clip = self%box
        if( file_exists(self%refs) )then
            call find_ldim_nptcls(self%refs, lfoo, ncls)
            if( cline%defined('ncls') )then
                if( ncls /= self%ncls )then
                    write(logfhandle,*)'ncls in ', self%refs%to_char(), ' : ', ncls
                    write(logfhandle,*)'self%ncls : ', self%ncls
                    THROW_HARD('input number of clusters (ncls) not consistent with the number of references in stack (p%refs)')
                endif
            else
                self%ncls = ncls
            endif
        endif
        self%l_remap_cls = self%remap_cls .eq. 'yes'
        if( .not. cline%defined('top') )   self%top   = self%nptcls
        if( .not. cline%defined('noris') ) self%noris = self%nptcls
        if( .not. cline%defined('moldiam') ) self%moldiam = 2. * self%msk * self%smpd
    end subroutine derive_image_settings

    module subroutine validate_parameter_consistency(self, cline)
        class(parameters), intent(inout) :: self
        class(cmdline),    intent(inout) :: cline
        type(atoms) :: atoms_obj
        if( self%walltime <= 0 )then
            THROW_HARD('Walltime cannot be negative!')
        else if( self%walltime > WALLTIME_DEFAULT )then
            THROW_WARN('Walltime is superior to maximum allowed of 23h59mins, resetting it to this value')
            self%walltime = WALLTIME_DEFAULT
        endif
        if( self%scale < 0.00001 )then
            THROW_HARD('scale out if range, should be > 0')
        endif
        if( cline%defined('msk') )then
            THROW_HARD('msk (mask radius in pixels) is deprecated! Use mskdiam (mask diameter in A)')
        endif
        if( cline%defined('mskfile') )then
            THROW_HARD('mskfile is no longer supported on command line; masks are internal and per-state')
        endif
        select case(trim(self%automsk))
            case('yes','tight','no')
            case DEFAULT
                THROW_HARD('Unsupported automsk mode: '//trim(self%automsk))
        end select
        select case(trim(self%objfun))
            case('cc')
                self%cc_objfun = OBJFUN_CC
            case('euclid')
                self%cc_objfun = OBJFUN_EUCLID
            case DEFAULT
                write(logfhandle,*) 'objfun flag: ', trim(self%objfun)
                THROW_HARD('unsupported objective function')
        end select
        self%l_dose_weight = cline%defined('total_dose')
        if( self%fraction_dose_target < 0.01 )then
            THROW_HARD('Invalid : fraction_dose_target'//real2str(self%fraction_dose_target))
        endif
        self%l_eer_fraction = cline%defined('eer_fraction')
        select case(self%eer_upsampling)
            case(1,2)
            case DEFAULT
                THROW_HARD('eer_upsampling not supported: '//int2str(self%eer_upsampling))
        end select
        select case(trim(self%center_type))
            case('mass','seg','params')
            case DEFAULT
                THROW_HARD('Unsupported centering scheme: '//trim(self%center_type))
        end select
        select case(trim(self%sigma_est))
            case('group')
                self%l_sigma_glob = .false.
            case('global')
                self%l_sigma_glob = .true.
            case DEFAULT
                THROW_HARD(trim(self%sigma_est)//' is not a supported sigma estimation approach')
        end select
        self%l_lpauto           = .false.
        self%l_nonuniform       = .false.
        self%l_nonuniform_lpset = .false.
        select case(trim(self%filt_mode))
            case('none')
            case('uniform')
                self%l_lpauto = .true.
                self%l_lpset = .true.
            case('fsc')
                self%l_lpauto = .true.
                self%l_lpset = .true.
            case('nonuniform')
                self%l_nonuniform = .true.
            case('nonuniform_lpset')
                self%l_nonuniform       = .true.
                self%l_nonuniform_lpset = .true.
            case DEFAULT
                THROW_HARD('unsupported filt_mode flag')
        end select
        self%l_noise_reg = cline%defined('snr_noise_reg')
        if( self%l_noise_reg )then
            self%eps_bounds(2) = self%snr_noise_reg
            self%eps_bounds(1) = self%eps_bounds(2) / 2.
            self%eps           = self%eps_bounds(1)
        endif
        self%l_lam_anneal = trim(self%lam_anneal).eq.'yes'
        self%l_ml_reg     = trim(self%ml_reg).eq.'yes'
        if( self%l_ml_reg ) self%l_ml_reg = self%cc_objfun == OBJFUN_EUCLID
        self%l_incrreslim = trim(self%incrreslim) == 'yes' .and. .not. self%l_lpset
        self%l_bfac       = cline%defined('bfac')
        if( cline%defined('element') )then
            if( .not. atoms_obj%element_exists(self%element) )then
                THROW_HARD('Element: '//trim(self%element)//' unsupported for now')
            endif
        endif
        select case(trim(self%imgkind))
            case('movie','mic','ptcl','cavg','cavg3D','vol','vol_cavg')
            case DEFAULT
                write(logfhandle,*) 'imgkind: ', trim(self%imgkind)
                THROW_HARD('unsupported imgkind')
        end select
        self%l_prob_align_mode = .false.
        select case(trim(self%refine))
            case('prob','prob_state','prob_neigh','prob_snhc')
                self%l_prob_align_mode = .true.
            case DEFAULT
        end select
        select case(trim(self%prob_neigh_mode))
            case('state','geom','sum','shc','snhc')
            case DEFAULT
                THROW_HARD('unsupported prob_neigh_mode; expected state|geom|sum|shc|snhc')
        end select
        self%l_neigh = .false.
        if( str_has_substr(self%refine, 'neigh') )then
            if( .not. cline%defined('nspace_sub') ) self%nspace_sub = 500
            if( .not. cline%defined('nspace') )     self%nspace     = 20000
            if( .not. cline%defined('athres') )     self%athres     = 10.
            self%nspace_sub = round2even(real(self%nspace_sub))
            self%nspace     = round2even(real(self%nspace))
            self%l_neigh    = .true.
        endif
        if( .not. cline%defined('trs') )then
            select case(trim(self%refine))
                case('snhc','snhc_smpl','snhc_smpl_many','prob_snhc')
                    self%trs = 0.
                case DEFAULT
                    if( trim(self%refine) == 'prob_neigh' .and. trim(self%prob_neigh_mode) == 'snhc' )then
                        self%trs = 0.
                    else
                        self%trs = MINSHIFT
                    endif
            end select
        endif
        self%trs = abs(self%trs)
        self%l_doshift = .true.
        if( self%trs < 0.1 )then
            self%trs = 0.
            self%l_doshift = .false.
        endif
        self%l_prob_inpl = trim(self%prob_inpl).eq.'yes'
        if( self%mcpatch.eq.'yes' .and. self%nxpatch*self%nypatch <= 1 ) self%mcpatch = 'no'
        if( self%mcpatch.eq.'no' )then
            self%nxpatch = 0
            self%nypatch = 0
        endif
        select case(trim(self%ptcl_src))
            case('raw','den')
            case DEFAULT
                THROW_HARD('Unsupported ptcl_src='//trim(self%ptcl_src)//'; expected raw|den')
        end select
        if( trim(self%ptcl_src) == 'den' .and. trim(self%oritype) /= 'ptcl3D' )then
            THROW_HARD('Denoised particle sources are supported only for oritype=ptcl3D')
        endif
        select case(trim(self%mcconvention))
            case('simple','unblur','motioncorr','relion','first','central','cryosparc','cs')
            case DEFAULT
                THROW_HARD('Invalid entry for MCCONVENTION='//trim(self%mcconvention))
        end select
        if( cline%defined('cls_init') )then
            select case(trim(self%cls_init))
                case('ptcl','rand','randcls')
                case DEFAULT
                    THROW_HARD('Unsupported mode of initial class generation CLS_INIT='//trim(self%cls_init))
            end select
        endif
    end subroutine validate_parameter_consistency

end submodule simple_parameters_phases
