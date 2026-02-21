!@descr: 3D refinement and reconstruction commanders used in SINGLE for nanoparticle processing
module single_commanders_nano3D
use simple_commanders_api
use simple_commanders_ori,   only: commander_vizoris
use simple_commanders_rec,    only: commander_reconstruct3D
use simple_commanders_volops, only: commander_reproject
use simple_commanders_cluster2D
use simple_nanoparticle
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_autorefine3D_nano
  contains
    procedure :: execute      => exec_autorefine3D_nano
end type commander_autorefine3D_nano

type, extends(commander_base) :: commander_refine3D_nano
  contains
    procedure :: execute      => exec_refine3D_nano
end type commander_refine3D_nano

type, extends(commander_base) :: commander_trajectory_reconstruct3D_distr
  contains
    procedure :: execute      => exec_commander_trajectory_reconstruct3D_distr
end type commander_trajectory_reconstruct3D_distr

contains

    subroutine exec_autorefine3D_nano( self, cline )
        use simple_commanders_atoms, only: commander_detect_atoms
        class(commander_autorefine3D_nano), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        type(parameters)              :: params
        type(commander_refine3D_nano) :: xrefine3D_nano
        type(commander_detect_atoms)  :: xdetect_atms
        type(commander_reproject)     :: xreproject
        type(commander_vizoris)       :: xvizoris
        type(commander_make_cavgs)    :: xmake_cavgs
        type(cmdline)                 :: cline_refine3D_nano, cline_detect_atms, cline_reproject
        type(cmdline)                 :: cline_vizoris, cline_make_cavgs
        type(image), allocatable      :: imgs(:)
        type(sp_project)              :: spproj
        character(len=*), parameter   :: RECVOL     = 'recvol_state01.mrc'
        character(len=*), parameter   :: EVEN       = 'recvol_state01_even.mrc'
        character(len=*), parameter   :: ODD        = 'recvol_state01_odd.mrc'
        character(len=*), parameter   :: SIMVOL     = 'recvol_state01_SIM.mrc'
        character(len=*), parameter   :: ATOMS      = 'recvol_state01_ATMS.pdb'
        character(len=*), parameter   :: BINARY     = 'recvol_state01_BIN.mrc'
        character(len=*), parameter   :: CCS        = 'recvol_state01_CC.mrc'
        character(len=*), parameter   :: SPLITTED   = 'split_ccs.mrc'
        character(len=*), parameter   :: FINAL_MAPS = './final_results/'
        character(len=*), parameter   :: TAG        = 'xxx' ! for checking command lines
        integer,          parameter   :: NSPACE_CLS3D = 500
        type(string)                  :: iter_dir, cavgs_stk, fname
        real,             allocatable :: rstates(:), corrs(:)
        logical,          allocatable :: state_mask(:)
        type(string) :: fbody, fbody_split, fname_reprojs, fname_reprojs_sim, fname_cvags_vs_reprojs
        integer      :: i, iter, cnt, cnt2, funit, io_stat, endit
        real         :: smpd
        logical      :: fall_over
        fbody       = get_fbody(RECVOL,   'mrc')
        fbody_split = get_fbody(SPLITTED, 'mrc')
        if(       cline%defined('nparts')         ) call cline%delete('nparts') ! shared-memory workflow
        if( .not. cline%defined('maxits')         ) call cline%set('maxits',          5)
        if( .not. cline%defined('maxits_between') ) call cline%set('maxits_between', 10)
        if( .not. cline%defined('overlap')        ) call cline%set('overlap',      0.98)
        if( .not. cline%defined('fracsrch')       ) call cline%set('fracsrch',      0.9)
        if( .not. cline%defined('objfun')         ) call cline%set('objfun',       'cc') ! needs to be here to avoid ERROR! file sigma2_it_10.star does not exist; simple_fileio.f90; line:   932
        if( .not. cline%defined('trail_rec')      ) call cline%set('trail_rec',   'yes') 
        if( .not. cline%defined('ufrac_trec')     ) call cline%set('ufrac_trec',    0.5)
        call cline%set('mkdir', 'yes') ! because we want to create the directory X_autorefine3D_nano & copy the project file
        call params%new(cline)         ! because the parameters class manages directory creation and project file copying, mkdir = yes
        params%mkdir = 'no'            ! to prevent the input vol to be appended with ../
        call cline%set('mkdir', 'no')  ! because we do not want a nested directory structure in the execution directory
        ! read the project file
        call spproj%read(params%projfile)
        call spproj%write_segment_inside('projinfo')
        ! sanity checks
        rstates = spproj%os_ptcl2D%get_all('state')
        fall_over = .false.
        if( any(rstates < 0.5 ) ) fall_over = .true.
        deallocate(rstates)
        rstates = spproj%os_ptcl3D%get_all('state')
        if( any(rstates < 0.5 ) ) fall_over = .true.
        if( fall_over ) THROW_HARD('There are state=0s in the ptcl2D/3D fields of the project, which is not allowed. Use simple_exec prg=prune_project before executing autorefine3D_nano')
        ! copy the input command line as templates for the refine3D_nano/detect_atoms command line
        cline_refine3D_nano = cline
        cline_detect_atms   = cline
        ! then update cline_refine3D_nano accordingly
        call cline_refine3D_nano%set('prg',     'refine3D_nano')
        call cline_refine3D_nano%set('projfile', params%projfile) ! since we are not making directories (non-standard execution) we need to keep track of project file
        call cline_refine3D_nano%set('keepvol',  'yes')
        call cline_refine3D_nano%set('maxits',   params%maxits_between) ! turn maxits_between into maxits (max # iterations between model building)
        call cline_refine3D_nano%delete('maxits_between')
        ! then update cline_detect_atoms accordingly
        call cline_detect_atms%set('prg', 'detect_atoms')
        call cline_detect_atms%set('vol1', RECVOL)               ! this is ALWYAS going to be the input volume to detect_atoms
        iter = 0
        do i = 1, params%maxits
            ! first refinement pass on the initial volume uses the low-pass limit defined by the user
            call xrefine3D_nano%execute_safe(cline_refine3D_nano)
            call cline_refine3D_nano%set('vol1', SIMVOL)         ! the reference volume is ALWAYS SIMVOL
            call cline_refine3D_nano%delete('lp')                ! uses the default 1.0 A low-pass limit
            endit = cline_refine3D_nano%get_iarg('endit')        ! last iteration executed by refine3D_nano
            call cline_refine3D_nano%delete('endit')             ! used internally but not technically allowed
            call cline_refine3D_nano%set('prg', 'refine3D_nano') ! because the command line is modified refine3D_nano -> refine3D internally
            ! model building
            call xdetect_atms%execute_safe(cline_detect_atms)
            ! copy critical output
            iter_dir = 'iteration_'//int2str_pad(i,2)//'/'
            call simple_mkdir(iter_dir)
            call simple_copy_file(string(RECVOL),   iter_dir//fbody//'_iter'//int2str_pad(i,3)//'.mrc')
            call simple_copy_file(string(EVEN),     iter_dir//fbody//'_iter'//int2str_pad(i,3)//'_even.mrc')
            call simple_copy_file(string(ODD),      iter_dir//fbody//'_iter'//int2str_pad(i,3)//'_odd.mrc')
            call simple_copy_file(string(SIMVOL),   iter_dir//fbody//'_iter'//int2str_pad(i,3)//'_SIM.mrc')
            call simple_copy_file(string(ATOMS),    iter_dir//fbody//'_iter'//int2str_pad(i,3)//'_ATMS.pdb')
            call simple_copy_file(string(BINARY),   iter_dir//fbody//'_iter'//int2str_pad(i,3)//'_BIN.mrc')
            call simple_copy_file(string(CCS),      iter_dir//fbody//'_iter'//int2str_pad(i,3)//'_CC.mrc')
            call simple_copy_file(string(SPLITTED), iter_dir//fbody_split//'_iter'//int2str_pad(i,3)//'.mrc')
            if( params%l_needs_sigma )then
                call simple_copy_file(string(SIGMA2_GROUP_FBODY)//int2str(endit)//STAR_EXT,&
                    &iter_dir//SIGMA2_GROUP_FBODY//int2str_pad(i,3)//STAR_EXT)
            endif
            ! clean
            call exec_cmdline('rm -f recvol_state01_iter*')
            if( params%l_needs_sigma ) call exec_cmdline('rm -f '//SIGMA2_GROUP_FBODY//'*'//STAR_EXT)
            call del_file(ATOMS)
            call del_file(BINARY)
            call del_file(CCS)
            call del_file(SPLITTED)
            iter = iter + 1
        end do
        call xdetect_atms%execute_safe(cline_detect_atms)
        call simple_mkdir(FINAL_MAPS)
        call simple_copy_file(string(RECVOL),   string(FINAL_MAPS)//fbody//'_iter'//int2str_pad(iter,3)//'.mrc')
        call simple_copy_file(string(EVEN),     string(FINAL_MAPS)//fbody//'_iter'//int2str_pad(iter,3)//'_even.mrc')
        call simple_copy_file(string(ODD),      string(FINAL_MAPS)//fbody//'_iter'//int2str_pad(iter,3)//'_odd.mrc')
        call simple_copy_file(string(SIMVOL),   string(FINAL_MAPS)//fbody//'_iter'//int2str_pad(iter,3)//'_SIM.mrc')
        call simple_copy_file(string(ATOMS),    string(FINAL_MAPS)//fbody//'_iter'//int2str_pad(iter,3)//'_ATMS.pdb')
        call simple_copy_file(string(BINARY),   string(FINAL_MAPS)//fbody//'_iter'//int2str_pad(iter,3)//'_BIN.mrc')
        call simple_copy_file(string(CCS),      string(FINAL_MAPS)//fbody//'_iter'//int2str_pad(iter,3)//'_CC.mrc')
        call simple_copy_file(string(SPLITTED), string(FINAL_MAPS)//fbody_split//'_iter'//int2str_pad(iter,3)//'.mrc')
        ! clean
        call del_file(SIMVOL)
        call del_file(ATOMS)
        call del_file(BINARY)
        call del_file(CCS)
        call del_file(SPLITTED)
        ! generate 3d classes
        cavgs_stk = 'cavgs3D.mrc'
        call cline_make_cavgs%set('prg',      'make_cavgs')
        call cline_make_cavgs%set('nspace',   NSPACE_CLS3D)
        call cline_make_cavgs%set('pgrp',     params%pgrp)
        call cline_make_cavgs%set('projfile', params%projfile)
        call cline_make_cavgs%set('nthr',     params%nthr)
        call cline_make_cavgs%set('mkdir',    'no')
        call cline_make_cavgs%set('refs',     cavgs_stk)
        call cline_make_cavgs%set('outfile',  'cavgs_oris.txt')
        call cline_make_cavgs%set('ml_reg',   'no')
        call xmake_cavgs%execute_safe(cline_make_cavgs)
        call spproj%os_cls3D%new(NSPACE_CLS3D, is_ptcl=.false.)
        call spproj%os_cls3D%read(string('cavgs_oris.txt')) ! will not be written as part of document
        if( allocated(rstates) ) deallocate(rstates)
        rstates = spproj%os_cls3D%get_all('state')
        ! prepare for re-projection
        call cline_reproject%set('vol1',   string(FINAL_MAPS)//fbody//'_iter'//int2str_pad(iter,3)//'.mrc')
        call cline_reproject%set('outstk', 'reprojs_recvol.mrc')
        call cline_reproject%set('smpd',   params%smpd)
        call cline_reproject%set('oritab', 'cavgs_oris.txt')
        call cline_reproject%set('pgrp',   params%pgrp)
        call cline_reproject%set('nthr',   params%nthr)
        call xreproject%execute_safe(cline_reproject)
        call cline_reproject%set('vol1',   string(FINAL_MAPS)//fbody//'_iter'//int2str_pad(iter,3)//'_SIM.mrc')
        call cline_reproject%set('outstk', 'reprojs_SIM.mrc')
        ! re-project
        call xreproject%execute_safe(cline_reproject)
        ! write cavgs & reprojections in triplets
        allocate(imgs(3), state_mask(NSPACE_CLS3D))
        call imgs(1)%new([params%box,params%box,1], smpd)
        call imgs(2)%new([params%box,params%box,1], smpd)
        call imgs(3)%new([params%box,params%box,1], smpd)
        cnt  = 0
        cnt2 = 1
        fname_reprojs          = 'reprojs_recvol.mrc'
        fname_reprojs_sim      = 'reprojs_SIM.mrc'
        fname_cvags_vs_reprojs = 'cavgs_vs_reprojections_rec_and_sim.mrc'
        do i = 1,3*NSPACE_CLS3D,3
            cnt = cnt + 1
            if( rstates(cnt) > 0.5 )then
                state_mask(cnt) = .true.
                call imgs(1)%read(cavgs_stk,         cnt)
                call imgs(2)%read(fname_reprojs,     cnt)
                call imgs(3)%read(fname_reprojs_sim, cnt)
                call imgs(1)%norm
                call imgs(2)%norm
                call imgs(3)%norm
                call imgs(1)%write(fname_cvags_vs_reprojs, cnt2    )
                call imgs(2)%write(fname_cvags_vs_reprojs, cnt2 + 1)
                call imgs(3)%write(fname_cvags_vs_reprojs, cnt2 + 2)
                cnt2 = cnt2 + 3
            else
                state_mask(cnt) = .false.
            endif
        end do
        call imgs(1)%kill
        call imgs(2)%kill
        call imgs(3)%kill
        deallocate(imgs)
        call exec_cmdline('rm -rf fsc* fft* recvol* RES* reprojs_recvol* reprojs* cavgs3D*mrc reproject_oris.txt cavgs_oris.txt stderrout')
        ! visualization of particle orientations
        ! read the ptcl3D segment first to make sure that we are using the latest information
        call spproj%read_segment('ptcl3D', params%projfile)
        ! extract ptcls oritab
        call spproj%os_ptcl3D%write(string('ptcls_oris.txt'))
        call cline_vizoris%set('oritab', 'ptcls_oris.txt')
        call cline_vizoris%set('pgrp',        params%pgrp)
        call cline_vizoris%set('nspace',     NSPACE_CLS3D)
        call cline_vizoris%set('tseries',           'yes')
        call xvizoris%execute_safe(cline_vizoris)
        ! print CSV file of correlation vs particle number
        corrs = spproj%os_ptcl3D%get_all('corr')
        fname = 'ptcls_vs_reprojs_corrs.csv'
        call fopen(funit, fname, 'replace', 'unknown', iostat=io_stat, form='formatted')
        call fileiochk('autorefine3D_nano fopen failed'//fname%to_char(), io_stat)
        write(funit,*) 'PTCL_INDEX'//CSV_DELIM//'CORR'
        do i = 1,size(corrs)
            write(funit,*) int2str(i)//CSV_DELIM//real2str(corrs(i))
        end do
        call fclose(funit)
        ! deallocate
        call iter_dir%kill
        call cavgs_stk%kill
        ! end gracefully
        call simple_end('**** AUTOREFINE3D_NANO NORMAL STOP ****')
    end subroutine exec_autorefine3D_nano

    subroutine exec_refine3D_nano( self, cline )
        use simple_commanders_refine3D, only: commander_refine3D_distr
        class(commander_refine3D_nano), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        ! commander
        type(commander_refine3D_distr) :: xrefine3D_distr
        ! static parameters
        call cline%set('prg',           'refine3D')
        call cline%set('dir_exec', 'refine3D_nano')
        ! dynamic parameters
        if( .not. cline%defined('cenlp')          ) call cline%set('cenlp',            5.)
        if( .not. cline%defined('graphene_filt')  ) call cline%set('graphene_filt',  'no') ! since Graphene subtraction is part of the workflow
        if( .not. cline%defined('keepvol')        ) call cline%set('keepvol',       'yes')
        if( .not. cline%defined('hp')             ) call cline%set('hp',              3.0)
        if( .not. cline%defined('lp')             ) call cline%set('lp',              1.0)
        if( .not. cline%defined('lpstart_nonuni') ) call cline%set('lpstart_nonuni',  2.5)
        if( .not. cline%defined('lpstop')         ) call cline%set('lpstop',          0.5)
        if( .not. cline%defined('maxits')         ) call cline%set('maxits',           30)
        if( .not. cline%defined('refine')         ) call cline%set('refine',      'neigh')
        if( .not. cline%defined('oritype')        ) call cline%set('oritype',    'ptcl3D')
        if( .not. cline%defined('trs')            ) call cline%set('trs',             5.0)
        if( .not. cline%defined('objfun')         ) call cline%set('objfun',         'cc') ! best objfun for this kind of data
        if( .not. cline%defined('ml_reg')         ) call cline%set('ml_reg',         'no') ! ml_reg=yes -> too few atoms 
        if( .not. cline%defined('sigma_est')      ) call cline%set('sigma_est',  'global') ! only sensible option for this kind of data
        if( .not. cline%defined('icm')            ) call cline%set('icm',           'yes') ! ICM regualrization works 
        if( .not. cline%defined('lambda')         ) call cline%set('lambda',          0.1) ! this is an empirically determined regularization parameter
        call xrefine3D_distr%execute_safe(cline)
    end subroutine exec_refine3D_nano

    subroutine exec_commander_trajectory_reconstruct3D_distr( self, cline )
        use simple_commanders_rec, only: commander_volassemble
        real, parameter :: LP_LIST(4) = [1.5,2.0,2.5,3.0]
        real, parameter :: HP_LIM = 5.0 ! no information at lower res for these kind of data
        class(commander_trajectory_reconstruct3D_distr), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        type(string),          allocatable :: vol_fnames(:)
        real,                  allocatable :: ccs(:,:,:), fsc(:), rstates(:), rad_cc(:), rad_dists(:)
        integer,               allocatable :: parts(:,:)
        type(string)                  :: recname, fname, str_state, recname_even, res_fname
        type(string)                  :: recname_odd, vol_fname_even, vol_fname_odd
        type(commander_reconstruct3D) :: xrec3D_shmem
        type(parameters)              :: params
        type(sp_project)              :: spproj
        type(cmdline)                 :: cline_rec
        type(image)                   :: vol1, vol2
        integer :: state, ipart, istate, nptcls, frame_start, frame_end
        integer :: funit, nparts, i, ind, nlps, ilp, iostat, hp_ind, lifetime
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',      'yes')
        if( .not. cline%defined('trs')     ) call cline%set('trs',           5.) ! to assure that shifts are being used
        if( .not. cline%defined('stepsz')  ) call cline%set('stepsz',      500.)
        if( .not. cline%defined('objfun')  ) call cline%set('objfun',      'cc') ! best objfun
        if( .not. cline%defined('ml_reg')  ) call cline%set('ml_reg',      'no') ! ml_reg=yes -> too few atoms 
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        call cline%delete('refine')
        call params%new(cline)
        call spproj%read(params%projfile)
        if( cline%defined('fromp') .and. cline%defined('top') )then
            call cline%delete('nparts')   ! shared-memory implementation
            call cline%set('mkdir', 'no') ! to avoid nested directory structure
            call xrec3D_shmem%execute(cline)
            return
        endif
        ! state exception
        rstates = spproj%os_ptcl3D%get_all('state')
        if( any(rstates < 0.5) ) THROW_HARD('state=0 entries not allowed, prune project beforehand')
        ! states/stepz
        nptcls = size(rstates)
        if( cline%defined('nparts') )then
            nparts = params%nparts
        else
            nparts = nint(real(nptcls)/real(params%stepsz))
        endif
        parts = split_nobjs_even(nptcls, nparts)
        allocate(vol_fnames(nparts), rad_cc(params%box/2), rad_dists(params%box/2))
        recname      = VOL_FBODY//int2str_pad(1,2)//params%ext%to_char()
        recname_even = VOL_FBODY//int2str_pad(1,2)//'_even'//params%ext%to_char()
        recname_odd  = VOL_FBODY//int2str_pad(1,2)//'_odd'//params%ext%to_char()
        fname        = 'lifetimes.csv'
        call fopen(funit, status='REPLACE', action='WRITE', file=fname, iostat=iostat)
        write(funit,*) 'PARTITION, ', 'FRAME_START, ', 'FRAME_END, ', 'LIFETIME'
        do ipart = 1,nparts
            str_state         = int2str_pad(ipart,2)
            vol_fnames(ipart) = 'partvol'//str_state%to_char()//params%ext%to_char()
            vol_fname_even    = 'partvol'//str_state%to_char()//'_even'//params%ext%to_char()
            vol_fname_odd     = 'partvol'//str_state%to_char()//'_odd'//params%ext%to_char()
            res_fname         = 'RESOLUTION_STATE'//str_state%to_char()
            ! prep 3D rec command line
            cline_rec = cline
            call cline_rec%delete('nparts') ! shared-memory implementation
            call cline_rec%set('fromp', parts(ipart,1))
            call cline_rec%set('top',   parts(ipart,2))
            frame_start = spproj%os_ptcl3D%get(parts(ipart,1), 'pind')
            frame_end   = spproj%os_ptcl3D%get(parts(ipart,2), 'pind')
            lifetime    = frame_end - frame_start + 1
            write(funit,'(I6,I6,I6,I6)') ipart, frame_start, frame_end, lifetime
            call cline_rec%set('mkdir', 'no')
            ! rec
            call xrec3D_shmem%execute(cline_rec)
            ! rename volumes and resolution files
            call simple_rename(recname,      vol_fnames(ipart))
            call simple_rename(recname_even, vol_fname_even)
            call simple_rename(recname_odd,  vol_fname_odd)
            if( ipart == 1 )then
                call simple_rename('RESOLUTION_STATE01', 'RESOLUTION_FILE_FIRST')
            else
                call simple_rename('RESOLUTION_STATE01',  res_fname)
            endif
        end do
        call fclose(funit)
        call simple_rename('RESOLUTION_FILE_FIRST', 'RESOLUTION_STATE01')
        ! Calculate correlation matrices
        nlps   = size(LP_LIST)
        hp_ind = calc_fourier_index(HP_LIM, params%box, params%smpd)
        allocate(fsc(fdim(params%box)-1),ccs(nlps,nparts,nparts))
        call vol1%new([params%box,params%box,params%box],params%smpd)
        call vol2%new([params%box,params%box,params%box],params%smpd)
        ccs = 1.
        do state = 1, nparts - 1
            call vol1%zero_and_unflag_ft
            call vol1%read(vol_fnames(state))
            call vol1%fft
            do istate = state + 1, nparts
                call vol2%zero_and_unflag_ft
                call vol2%read(vol_fnames(istate))
                call vol2%fft
                call vol1%fsc(vol2,fsc)
                do ilp = 1, nlps
                    ind = calc_fourier_index(LP_LIST(ilp), params%box, params%smpd)
                    ccs(ilp,state,istate) = sum(fsc(hp_ind:ind)) / real(ind - hp_ind + 1)
                    ccs(ilp,istate,state) = ccs(ilp,state,istate)
                enddo
            enddo
        enddo
        ! replace the diagonal elements with the maximum corr value in the columns
        ! for improved plotting / graphing
        do ilp = 1, nlps
            do istate = 1, nparts
                ccs(ilp,istate,istate) = maxval(ccs(ilp,istate,:), mask=ccs(ilp,istate,:) < 0.99)
            end do
        end do
        do ilp = 1,nlps
            fname = 'ccmat_lp'//trim(real2str(LP_LIST(ilp)))//'.csv'
            call fopen(funit, status='REPLACE', action='WRITE', file=fname, iostat=iostat)
            do istate = 1, nparts
                do i = 1, nparts - 1
                    write(funit,'(F8.3,A2)', advance='no') ccs(ilp,istate,i), ', '
                enddo
                write(funit,'(F8.3)', advance='yes') ccs(ilp,istate,nparts)
            enddo
            call fclose(funit)
            fname = 'lp'//trim(real2str(LP_LIST(ilp)))//'.txt'
            fname = 'ccneigh_lp'//trim(real2str(LP_LIST(ilp)))//'.csv'
            call fopen(funit, status='REPLACE', action='WRITE', file=fname, iostat=iostat)
            do istate = 1, nparts - 1
                write(funit,'(I3,A3,F8.3)',advance='yes') istate, ', ', ccs(ilp,istate,istate + 1)
            end do
            call fclose(funit)
        enddo
        call simple_end('**** SIMPLE_trajectory_reconstruct3D NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_commander_trajectory_reconstruct3D_distr

  end module single_commanders_nano3D
