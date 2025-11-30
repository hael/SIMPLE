! concrete commander: operations on orientations
module simple_commanders_oris
include 'simple_lib.f08'
use simple_binoris_io
use simple_cmdline,        only: cmdline
use simple_sp_project,     only: sp_project
use simple_commander_base, only: commander_base
use simple_parameters,     only: parameters
use simple_builder,        only: builder
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_make_oris
  contains
    procedure :: execute      => exec_make_oris
end type commander_make_oris

type, extends(commander_base) :: commander_orisops
  contains
    procedure :: execute      => exec_orisops
end type commander_orisops

type, extends(commander_base) :: commander_oristats
  contains
    procedure :: execute      => exec_oristats
end type commander_oristats

type, extends(commander_base) :: commander_oriconsensus
  contains
    procedure :: execute      => exec_oriconsensus
end type commander_oriconsensus

type, extends(commander_base) :: commander_rotmats2oris
  contains
    procedure :: execute      => exec_rotmats2oris
end type commander_rotmats2oris

type, extends(commander_base) :: commander_vizoris
  contains
    procedure :: execute      => exec_vizoris
end type commander_vizoris

type, extends(commander_base) :: commander_check_states
  contains
    procedure :: execute      => check_states
end type commander_check_states

contains

    !> for making SIMPLE orientation/parameter files
    subroutine exec_make_oris( self, cline )
        class(commander_make_oris), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        type(ori)        :: orientation
        type(oris)       :: os_even, spiral
        real             :: e3, x, y
        integer          :: i, class
        call build%init_params_and_build_general_tbox(cline,params,do3d=.false.)
        if( cline%defined('ncls') )then
            os_even = oris(params%ncls, is_ptcl=.false.)
            call build%pgrpsyms%build_refspiral(os_even)
            call build%spproj_field%new(params%nptcls, is_ptcl=.false.)
            do i=1,params%nptcls
                class = irnd_uni(params%ncls)
                call os_even%get_ori(class, orientation)
                call build%spproj_field%set_ori(i, orientation)
                e3 = ran3()*2.*params%angerr-params%angerr
                x  = ran3()*2.0*params%sherr-params%sherr
                y  = ran3()*2.0*params%sherr-params%sherr
                call build%spproj_field%set(i, 'x', x)
                call build%spproj_field%set(i, 'y', y)
                call build%spproj_field%e3set(i, e3)
                call build%spproj_field%set(i, 'class', class)
            end do
        else if( cline%defined('ndiscrete') )then
            if( params%ndiscrete > 0 )then
                call spiral%new(params%ndiscrete, is_ptcl=.false.)
                call build%pgrpsyms%build_refspiral(spiral)
            endif
            call build%spproj_field%rnd_inpls(params%sherr)
        else if( params%even .eq. 'yes' )then
            call build%pgrpsyms%build_refspiral(build%spproj_field)
            call build%spproj_field%rnd_inpls(params%sherr)
        else
            call build%spproj_field%rnd_oris(params%sherr)
            call build%pgrpsyms%rotall_to_asym(build%spproj_field)
            if( params%doprint .eq. 'yes' )then
                call build%spproj_field%print_matrices
            endif
        endif
        if( cline%defined('defocus') )then
            call build%spproj_field%rnd_ctf(0., 0., 0., params%defocus, params%dferr, params%astigerr)
            call build%spproj_field%delete_entry('cs')
            call build%spproj_field%delete_entry('kv')
            call build%spproj_field%delete_entry('fraca')
        endif
        call build%spproj_field%set_all2single('state',1.)
        if( params%nstates > 1 ) call build%spproj_field%rnd_states(params%nstates)
        call build%spproj_field%set_all2single('w',1.)
        call binwrite_oritab(params%outfile, build%spproj, build%spproj_field, [1,build%spproj_field%get_noris()])
        call orientation%kill
        call os_even%kill
        call spiral%kill
        ! end gracefully
        call simple_end('**** SIMPLE_MAKE_ORIS NORMAL STOP ****')
    end subroutine exec_make_oris

    subroutine exec_orisops( self, cline )
        class(commander_orisops), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        type(ori)        :: orientation
        integer          :: s, i
        call build%init_params_and_build_spproj(cline,params)
        if( cline%defined('angerr') )then
            ! introduce error in input orientations
            call build%spproj_field%introd_alig_err(params%angerr, params%sherr)
        endif
        if( cline%defined('dferr') )then
            ! introduce error in defocus parameters
            call build%spproj_field%introd_ctf_err(params%dferr)
        endif
        if( cline%defined('sherr') )then
            ! introduce gaussian error to 2D shifts
            call build%spproj_field%gau_rnd_shifts(params%sherr)
        endif
        if( cline%defined('e1') .or.&
            cline%defined('e2') .or.&
            cline%defined('e3') )then
            ! rotate input Eulers
            call orientation%new(is_ptcl=.false.)
            call orientation%set_euler([params%e1,params%e2,params%e3])
            if( cline%defined('state') )then
                do i=1,build%spproj_field%get_noris()
                    s = build%spproj_field%get_state(i)
                    if( s == params%state )then
                        call build%spproj_field%rot_transp(i,orientation)
                    endif
                end do
            else
                call build%spproj_field%rot_transp(orientation)
            endif
        endif
        if( cline%defined('mul') )       call build%spproj_field%mul_shifts(params%mul)
        if( params%zero .eq. 'yes' )     call build%spproj_field%zero_shifts
        if( cline%defined('ndiscrete') ) call build%spproj_field%discretize(params%ndiscrete)
        if( params%symrnd .eq. 'yes' )   call build%pgrpsyms%symrandomize(build%spproj_field)
        if( cline%defined('nstates') )   call build%spproj_field%rnd_states(params%nstates)
        if( cline%defined('mirr') )then
            select case(trim(params%mirr))
                case('2d')
                    call build%spproj_field%mirror2d()
                case('3d')
                    call build%spproj_field%mirror3d()
                case('no')
                    ! nothing to do
                case DEFAULT
                    THROW_HARD('mirr flag: '//trim(params%mirr)//' is unsupported; exec_orisops')
            end select
        endif
        call binwrite_oritab(params%outfile, build%spproj, build%spproj_field, [1,build%spproj_field%get_noris()])
        call simple_end('**** SIMPLE_ORISOPS NORMAL STOP ****')
    end subroutine exec_orisops

    !> for analyzing SIMPLE orientation/parameter files
    subroutine exec_oristats( self, cline )
        class(commander_oristats), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        type(parameters)     :: params
        type(builder)        :: build
        type(sp_project)     :: spproj
        class(oris), pointer :: o => null()
        type(oris)           :: osubspace
        type(ori)            :: o_single
        real                 :: mind, maxd, avgd, sdevd, sumd, vard, scale
        real                 :: mind2, maxd2, avgd2, sdevd2, vard2
        real                 :: popmin, popmax, popmed, popave, popsdev, popvar, frac_populated
        real                 :: resmin, resmax, resmed, resave, ressdev, resvar, nr, nr_tot
        integer              :: nprojs, iptcl, icls, j, noris, ncls, pop
        real,    allocatable :: clustszs(:), pops(:), res(:), df(:), vals(:)
        integer, allocatable :: clustering(:), tmp(:), states(:), cls(:)
        logical, allocatable :: ptcl_mask(:)
        integer, parameter   :: hlen=50, NSPACE_REDUCED = 600
        logical              :: err
        call build%init_params_and_build_general_tbox(cline,params,do3d=.false.)
        if( cline%defined('oritab2') )then
            ! Comparison
            if( .not. cline%defined('oritab') ) THROW_HARD('need oritab for comparison')
            if( binread_nlines( params%oritab) .ne. binread_nlines( params%oritab2) )then
                THROW_HARD('inconsistent number of lines in the two oritabs!')
            endif
            call spproj%new_seg_with_ptr(params%nptcls, params%oritype, o)
            call binread_oritab(params%oritab2, spproj, o, [1,params%nptcls])
            call build%spproj_field%diststat(o, sumd, avgd, sdevd, mind, maxd)
            write(logfhandle,'(a,1x,f15.6)') 'SUM OF ANGULAR DISTANCE BETWEEN ORIENTATIONS  :', sumd
            write(logfhandle,'(a,1x,f15.6)') 'AVERAGE ANGULAR DISTANCE BETWEEN ORIENTATIONS :', avgd
            write(logfhandle,'(a,1x,f15.6)') 'STANDARD DEVIATION OF ANGULAR DISTANCES       :', sdevd
            write(logfhandle,'(a,1x,f15.6)') 'MINIMUM ANGULAR DISTANCE                      :', mind
            write(logfhandle,'(a,1x,f15.6)') 'MAXIMUM ANGULAR DISTANCE                      :', maxd
        else if( cline%defined('oritab') )then
            if( params%ctfstats .eq. 'yes' )then
                call build%spproj_field%stats('ctfres', avgd, sdevd, vard, err )
                call build%spproj_field%minmax('ctfres', mind, maxd)
                write(logfhandle,'(a,1x,f8.2)') 'AVERAGE CTF RESOLUTION               :', avgd
                write(logfhandle,'(a,1x,f8.2)') 'STANDARD DEVIATION OF CTF RESOLUTION :', sdevd
                write(logfhandle,'(a,1x,f8.2)') 'MINIMUM CTF RESOLUTION (BEST)        :', mind
                write(logfhandle,'(a,1x,f8.2)') 'MAXIMUM CTF RESOLUTION (WORST)       :', maxd
                call build%spproj_field%stats('dfx', avgd, sdevd, vard, err )
                call build%spproj_field%minmax('dfx', mind, maxd)
                call build%spproj_field%stats('dfy', avgd2, sdevd2, vard2, err )
                call build%spproj_field%minmax('dfy', mind2, maxd2)
                write(logfhandle,'(a,1x,f8.2)') 'AVERAGE DF                           :', (avgd+avgd2)/2.
                write(logfhandle,'(a,1x,f8.2)') 'STANDARD DEVIATION OF DF             :', (sdevd+sdevd2)/2.
                write(logfhandle,'(a,1x,f8.2)') 'MINIMUM DF                           :', (mind+mind2)/2.
                write(logfhandle,'(a,1x,f8.2)') 'MAXIMUM DF                           :', (maxd+maxd2)/2.
                if( trim(params%oritype) .eq. 'ptcl2D' )then
                    states = build%spproj_field%get_all_asint('state')
                    cls    = build%spproj_field%get_all_asint('class')
                    ncls   = maxval(cls)
                    df     = (build%spproj_field%get_all('dfx') + build%spproj_field%get_all('dfy')) / 2.0
                    do icls = 1,ncls
                        vals = pack(df, mask=(cls==icls).and.(states>0))
                        if( .not.allocated(vals) ) cycle
                        pop = size(vals)
                        if( pop == 0 ) cycle
                        avgd  = sum(vals) / real(pop)
                        sdevd = sqrt(sum(vals) / real(pop))
                        write(logfhandle,'(a,1x,I4,2f8.2)')'CLASS, AVERAGE DEFOCUS & STANDARD DEV:', icls, avgd, sdevd
                    enddo
                endif
            endif
            if( params%classtats .eq. 'yes' )then
                noris  = build%spproj_field%get_noris()
                ncls   = build%spproj_field%get_n('class')
                ! generate pop stats
                pops   = build%spproj_field%get_all('pop')
                nr_tot = sum(pops)
                popmin = minval(pops)
                popmax = maxval(pops)
                popmed = median(real(pops))
                call moment(pops, popave, popsdev, popvar, err)
                write(logfhandle,'(a)') 'CLASS POPULATION STATS'
                write(logfhandle,'(a,1x,f8.2)') 'MINIMUM POPULATION :', popmin
                write(logfhandle,'(a,1x,f8.2)') 'MAXIMUM POPULATION :', popmax
                write(logfhandle,'(a,1x,f8.2)') 'MEDIAN  POPULATION :', popmed
                write(logfhandle,'(a,1x,f8.2)') 'AVERAGE POPULATION :', popave
                write(logfhandle,'(a,1x,f8.2)') 'SDEV OF POPULATION :', popsdev
                ! generate res stats
                res    = build%spproj_field%get_all('res')
                resmin = minval(res, mask=pops > 0.1)
                resmax = maxval(res, mask=pops > 0.1)
                resmed = median(res)
                call moment(res, resave, ressdev, resvar, err)
                write(logfhandle,'(a)') ''
                write(logfhandle,'(a)') 'CLASS RESOLUTION STATS'
                write(logfhandle,'(a,1x,f8.2)') 'MINIMUM RESOLUTION :', resmin
                write(logfhandle,'(a,1x,f8.2)') 'MAXIMUM RESOLUTION :', resmax
                write(logfhandle,'(a,1x,f8.2)') 'MEDIAN  RESOLUTION :', resmed
                write(logfhandle,'(a,1x,f8.2)') 'AVERAGE RESOLUTION :', resave
                write(logfhandle,'(a,1x,f8.2)') 'SDEV OF RESOLUTION :', ressdev
                nr  = sum(pack(pops, mask=(res < 10. .and. pops > 0.1)))
                res = pack(res, mask=(res < 10. .and. pops > 0.1))
                write(logfhandle,'(a)') ''
                write(logfhandle,'(a)') 'CLASS RESOLUTION STATS FOR CLASSES BETTER THAN 10 A'
                resmin = minval(res)
                resmax = maxval(res)
                resmed = median(res)
                call moment(res, resave, ressdev, resvar, err)
                write(logfhandle,'(a,1x,f8.2)') '% PARTICLES        :', (nr / nr_tot) * 100. 
                write(logfhandle,'(a,1x,f8.2)') 'MEDIAN  RESOLUTION :', resmed
                write(logfhandle,'(a,1x,f8.2)') 'AVERAGE RESOLUTION :', resave
                write(logfhandle,'(a,1x,f8.2)') 'SDEV OF RESOLUTION :', ressdev
                ! produce a histogram of class populations
                ! szmax = maxval(pops)
                ! ! scale to max 50 *:s
                ! scale = 1.0
                ! do while( nint(scale*szmax) > hlen )
                !     scale = scale - 0.001
                ! end do
                ! write(logfhandle,'(a)') '>>> HISTOGRAM OF CLASS POPULATIONS'
                ! do icls=1,ncls
                !     write(logfhandle,*) pops(icls),"|",('*', j=1,nint(real(pops(icls)*scale)))
                ! end do
            endif
            if( params%projstats .eq. 'yes' )then
                if( .not. cline%defined('nspace') ) THROW_HARD('need nspace command line arg to provide projstats')
                noris = build%spproj_field%get_noris()
                ! setup weights
                call build%spproj_field%calc_hard_weights(params%frac)
                ! generate population stats
                call build%spproj_field%get_pops(tmp, 'proj')
                nprojs         = size(tmp)
                pops           = pack(tmp, tmp > 0.5)                   !! realloc warning
                frac_populated = real(size(pops))/real(params%nspace)
                popmin         = minval(pops)
                popmax         = maxval(pops)
                popmed         = median(real(pops))
                call moment(real(pops), popave, popsdev, popvar, err)
                write(logfhandle,'(a)') '>>> STATISTICS BEFORE CLUSTERING'
                write(logfhandle,'(a,1x,f8.2)') 'FRAC POPULATED DIRECTIONS :', frac_populated
                write(logfhandle,'(a,1x,f8.2)') 'MINIMUM POPULATION        :', popmin
                write(logfhandle,'(a,1x,f8.2)') 'MAXIMUM POPULATION        :', popmax
                write(logfhandle,'(a,1x,f8.2)') 'MEDIAN  POPULATION        :', popmed
                write(logfhandle,'(a,1x,f8.2)') 'AVERAGE POPULATION        :', popave
                write(logfhandle,'(a,1x,f8.2)') 'SDEV OF POPULATION        :', popsdev
                ! produce a histogram based on clustering into NSPACE_REDUCED even directions
                ! first, generate a mask based on state flag and w
                ptcl_mask = build%spproj_field%included()
                allocate(clustering(noris), clustszs(NSPACE_REDUCED))
                call osubspace%new(NSPACE_REDUCED, is_ptcl=.false.)
                call build%pgrpsyms%build_refspiral(osubspace)
                call osubspace%write(string('even_pdirs'//TXT_EXT), [1,NSPACE_REDUCED])
                do iptcl=1,build%spproj_field%get_noris()
                    if( ptcl_mask(iptcl) )then
                        call build%spproj_field%get_ori(iptcl, o_single)
                        clustering(iptcl) = osubspace%find_closest_proj(o_single)
                    else
                        clustering(iptcl) = 0
                    endif
                end do
                ! determine cluster sizes
                do icls=1,NSPACE_REDUCED
                    clustszs(icls) = real(count(clustering == icls))
                end do
                frac_populated = real(count(clustszs > 0.5))/real(NSPACE_REDUCED)
                popmin         = minval(clustszs)
                popmax         = maxval(clustszs)
                popmed         = median_nocopy(clustszs)
                call moment(clustszs, popave, popsdev, popvar, err)
                write(logfhandle,'(a)') '>>> STATISTICS AFTER CLUSTERING'
                write(logfhandle,'(a,1x,f8.2)') 'FRAC POPULATED DIRECTIONS :', frac_populated
                write(logfhandle,'(a,1x,f8.2)') 'MINIMUM POPULATION        :', popmin
                write(logfhandle,'(a,1x,f8.2)') 'MAXIMUM POPULATION        :', popmax
                write(logfhandle,'(a,1x,f8.2)') 'MEDIAN  POPULATION        :', popmed
                write(logfhandle,'(a,1x,f8.2)') 'AVERAGE POPULATION        :', popave
                write(logfhandle,'(a,1x,f8.2)') 'SDEV OF POPULATION        :', popsdev
                ! scale to max 50 *:s
                scale = 1.0
                do while( nint(scale * popmax) > hlen )
                    scale = scale - 0.001
                end do
                write(logfhandle,'(a)') '>>> HISTOGRAM OF SUBSPACE POPULATIONS (FROM NORTH TO SOUTH)'
                do icls=1,NSPACE_REDUCED
                    write(logfhandle,*) nint(clustszs(icls)),"|",('*', j=1,nint(clustszs(icls)*scale))
                end do
            endif
            if( params%trsstats .eq. 'yes' )then
                call build%spproj_field%stats('x', avgd, sdevd, vard, err )
                call build%spproj_field%minmax('x', mind, maxd)
                call build%spproj_field%stats('y', avgd2, sdevd2, vard2, err )
                call build%spproj_field%minmax('y', mind2, maxd2)
                write(logfhandle,'(a,1x,f8.2)') 'AVERAGE TRS               :', (avgd+avgd2)/2.
                write(logfhandle,'(a,1x,f8.2)') 'STANDARD DEVIATION OF TRS :', (sdevd+sdevd2)/2.
                write(logfhandle,'(a,1x,f8.2)') 'MINIMUM TRS               :', (mind+mind2)/2.
                write(logfhandle,'(a,1x,f8.2)') 'MAXIMUM TRS               :', (maxd+maxd2)/2.
            endif
        endif
        call osubspace%kill
        call o_single%kill
        call simple_end('**** SIMPLE_ORISTATS NORMAL STOP ****')
    end subroutine exec_oristats

    subroutine exec_oriconsensus( self, cline )
        class(commander_oriconsensus), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(sp_project),          allocatable :: spprojs(:)
        type(string), allocatable :: projfnames(:)
        integer,      allocatable :: states(:,:), states_consensus(:)
        type(parameters) :: params
        integer          :: nprojs, iproj, noris1, noris2, ngood, nbad, iori
        real             :: overlap
        if( .not. cline%defined('oritab') ) THROW_HARD('Need list of project files (oritab) to analyze!')
        call params%new(cline)
        call read_filetable(params%oritab, projfnames)
        nprojs = size(projfnames)
        allocate(spprojs(nprojs))
        do iproj = 1, nprojs
            call spprojs(iproj)%read(projfnames(iproj))
            if( iproj == 1 )then
                noris1 = spprojs(iproj)%os_ptcl2D%get_noris()
            else
                noris2 = spprojs(iproj)%os_ptcl2D%get_noris()
                if( noris1 /= noris2 ) THROW_HARD('Nonconforming # entries in ptcl2D field!')
            endif
        end do

        print *, 'read project files'

        allocate(states(nprojs,noris1), states_consensus(noris1), source=0)

        !$omp parallel do default(shared) collapse(2) private(iproj,iori) proc_bind(close)
        do iproj = 1, nprojs
            do iori = 1, noris1
                states(iproj,iori) = spprojs(iproj)%os_ptcl2D%get_state(iori)
            end do
        end do
        !$omp end parallel do
        
        print *, 'extracted states'

        !$omp parallel do default(shared) private(iori,ngood,nbad) proc_bind(close)
        do iori = 1, noris1
            ngood = count(states(:,iori) >  0)
            nbad  = count(states(:,iori) == 0)
            if( ngood > nbad )then
                states_consensus(iori) = 1
            else
                states_consensus(iori) = 0
            endif
        end do
        !$omp end parallel do

        print *, 'calculated consensus'
        
        do iproj = 1, nprojs
            ngood = count(states(iproj,:) >  0)
            print *, iproj, ' ngood ', ngood
            overlap = real(count(states(iproj,:) == states_consensus)) / real(noris1)
            print *, iproj, ' consensus overlap ', overlap

        end do
        print *, ' consensus ngood ',  count(states_consensus > 0)

        call simple_end('**** SIMPLE_ORICONSENSUS NORMAL STOP ****')
    end subroutine exec_oriconsensus

    !> convert rotation matrix to orientation oris class
    subroutine exec_rotmats2oris( self, cline )
        class(commander_rotmats2oris),  intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(parameters) :: params
        type(nrtxtfile) :: rotmats
        type(ori)       :: o
        type(oris)      :: os_out
        integer         :: nrecs_per_line, iline, ndatlines
        real            :: rline(9), rmat(3,3)
        if( .not. cline%defined('outfile') ) call cline%set('outfile', 'outfile.txt')
        call params%new(cline)
        if( cline%defined('infile') )then
            call rotmats%new(params%infile, 1)
            ndatlines = rotmats%get_ndatalines()
        else
            THROW_HARD('Need infile defined on command line: text file with 9 records per line defining a rotation matrix (11) (12) (13) (21) etc.')
        endif
        if( fname2format(params%outfile) .eq. '.simple' )then
            THROW_HARD('*.simple outfile not supported; rotmats2oris')
        endif
        nrecs_per_line = rotmats%get_nrecs_per_line()
        if( nrecs_per_line /= 9 ) THROW_HARD('need 9 records (real nrs) per line of file (infile) describing rotation matrices')
        call os_out%new(ndatlines, is_ptcl=.false.)
        do iline=1,ndatlines
            call rotmats%readNextDataLine(rline)
            rmat(1,1) = rline(1)
            rmat(1,2) = rline(2)
            rmat(1,3) = rline(3)
            rmat(2,1) = rline(4)
            rmat(2,2) = rline(5)
            rmat(2,3) = rline(6)
            rmat(3,1) = rline(7)
            rmat(3,2) = rline(8)
            rmat(3,3) = rline(9)
            rmat = transpose(rmat)
            call o%ori_from_rotmat(rmat, is_ptcl=.false.)
            call os_out%set_ori(iline,o)
        end do
        call os_out%swape1e3
        call os_out%write(params%outfile, [1,ndatlines])
        call rotmats%kill
        call simple_end('**** ROTMATS2ORIS NORMAL STOP ****')
    end subroutine exec_rotmats2oris

    subroutine exec_vizoris( self, cline )
        class(commander_vizoris), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(parameters)      :: params
        type(builder)         :: build
        type(ori)             :: o, o_prev
        real,    allocatable  :: euldists(:)
        integer, allocatable  :: pops(:)
        type(string)          :: fname, ext
        integer               :: i, n, maxpop, funit, closest,io_stat
        real                  :: radius, maxradius, ang, scale, col, avg_geodist,avg_euldist,geodist
        real                  :: xyz(3), xyz_end(3), xyz_start(3), vec(3)
        call build%init_params_and_build_general_tbox(cline,params,do3d=.true.)
        n = build%spproj_field%get_noris()
        if( .not.cline%defined('fbody') )then
            fname = basename(params%oritab)
            ext   = fname2ext(fname)
            params%fbody = get_fbody(fname, ext)
        endif
        if( params%tseries.eq.'no' )then
            ! Discretization of the projection directions
            ! init
            allocate(pops(params%nspace), source=0)
            ang = 3.6 / sqrt(real(params%nsym*params%nspace))
            maxradius = 0.75 * sqrt( (1.-cos(ang))**2. + sin(ang)**2. )
            ! projection direction attribution
            n = build%spproj_field%get_noris()
            do i = 1, n
                call build%spproj_field%get_ori(i, o)
                if( o%isstatezero() )cycle
                call progress(i, n)
                closest = build%eulspace%find_closest_proj(o)
                pops(closest) = pops(closest) + 1
            enddo
            maxpop = maxval(pops)
            write(logfhandle,'(A,I6)')'>>> NUMBER OF POPULATED PROJECTION DIRECTIONS:', count(pops>0)
            write(logfhandle,'(A,I6)')'>>> NUMBER OF EMPTY     PROJECTION DIRECTIONS:', count(pops==0)
            ! output
            fname = params%fbody//'.bild'
            call fopen(funit, status='REPLACE', action='WRITE', file=fname,iostat=io_stat)
             if(io_stat/=0)call fileiochk("simple_commanders_oris::exec_vizoris fopen failed "//fname%to_char(), io_stat)
            ! header
            write(funit,'(A)')".translate 0.0 0.0 0.0"
            write(funit,'(A)')".scale 10"
            write(funit,'(A)')".comment -- unit sphere --"
            write(funit,'(A)')".color 0.8 0.8 0.8"
            write(funit,'(A)')".sphere 0 0 0 1.0"
            write(funit,'(A)')".comment -- planes --"
            write(funit,'(A)')".color 0.3 0.3 0.3"
            write(funit,'(A)')".cylinder -0.02 0 0 0.02 0 0 1.02"
            write(funit,'(A)')".cylinder 0 -0.02 0 0 0.02 0 1.02"
            write(funit,'(A)')".cylinder 0 0 -0.02 0 0 0.02 1.02"
            write(funit,'(A)')".comment -- x-axis --"
            write(funit,'(A)')".color 1 0 0"
            write(funit,'(A)')".cylinder -1.5 0 0 1.5 0 0 0.02"
            write(funit,'(A)')".comment -- y-axis --"
            write(funit,'(A)')".color 0 1 0"
            write(funit,'(A)')".cylinder 0 -1.5 0 0 1.5 0 0.02"
            write(funit,'(A)')".comment -- z-axis --"
            write(funit,'(A)')".color 0 0 1"
            write(funit,'(A)')".cylinder 0 0 -1.5 0 0 1.5 0.02"
            ! body
            write(funit,'(A)')".comment -- projection firections --"
            write(funit,'(A)')".color 0.4 0.4 0.4"
            do i = 1, params%nspace
                if( pops(i) == 0 )cycle
                scale     = real(pops(i)) / real(maxpop)
                xyz_start = build%eulspace%get_normal(i)
                xyz_end   = (1.05 + scale/4.) * xyz_start
                radius    = max(maxradius * scale, 0.002)
                write(funit,'(A,F7.3,F7.3,F7.3,F7.3,F7.3,F7.3,F6.3)')&
                &'.cylinder ', xyz_start, xyz_end, radius
            enddo
            call fclose(funit)
        else
            ! time series
            ! unit sphere tracking
            fname  = params%fbody//'_motion.bild'
            radius = 0.02
            call fopen(funit, status='REPLACE', action='WRITE', file=fname, iostat=io_stat)
             if(io_stat/=0)call fileiochk("simple_commanders_oris::exec_vizoris fopen failed ", io_stat)
            write(funit,'(A)')".translate 0.0 0.0 0.0"
            write(funit,'(A)')".scale 1"
            write(funit,'(A)')".sphere 0 0 0 0.1"
            write(funit,'(A)')".color 1 0 0"
            write(funit,'(A)')".arrow 1.2 0 0 1.3 0 0 0.05 0.1 0.3"
            write(funit,'(A)')".color 1 1 0"
            write(funit,'(A)')".arrow 0 1.2 0 0 1.3 0 0.05 0.1 0.3"
            write(funit,'(A)')".color 0 0 1"
            write(funit,'(A)')".arrow 0 0 1.2 0 0 1.3 0.05 0.1 0.3" ! North pole in blue
            do i = 1, n
                call build%spproj_field%get_ori(i, o)
                if( o%e2get() > 90. ) call o%mirror2d()
                xyz = o%get_normal()
                col = real(i-1)/real(n-1)
                write(funit,'(A,F6.2,F6.2)')".color 1.0 ", col, col
                if( i==1 )then
                    write(funit,'(A,F7.3,F7.3,F7.3,A)')".sphere ", xyz, " 0.08"
                else
                    vec = xyz - xyz_start
                    if( sqrt(dot_product(vec, vec)) > 0.01 )then
                        write(funit,'(A,F7.3,F7.3,F7.3,A)')".sphere ", xyz, " 0.02"
                    endif
                endif
                xyz_start = xyz
            enddo
            write(funit,'(A,F7.3,F7.3,F7.3,A)')".sphere ", xyz, " 0.08"
            call fclose(funit)
            ! distance output
            avg_geodist = 0.
            avg_euldist = 0.
            allocate(euldists(n))
            fname  = params%fbody//'_motion.csv'
            call fopen(funit, status='REPLACE', action='WRITE', file=fname, iostat=io_stat)
            if(io_stat/=0)call fileiochk("simple_commanders_oris::exec_vizoris fopen failed "//fname%to_char(), io_stat)
            do i = 1, n
                call build%spproj_field%get_ori(i, o)
                if( i==1 )then
                    ang     = 0.
                    geodist = 0.
                else
                    ang     = rad2deg(o_prev.euldist.o)
                    geodist = o_prev.geod.o
                    call o_prev%mirror2d
                    ang     = min(ang, rad2deg(o_prev.euldist.o))
                    geodist = min(geodist, o_prev.geod.o)
                    avg_euldist = avg_euldist + ang
                    avg_geodist = avg_geodist + geodist
                endif
                euldists(i) = ang
                write(funit,'(I7,A1,F8.3,A1,F8.3)')i, ',', ang, ',', geodist
                o_prev = o
            enddo
            call fclose(funit)
            avg_geodist = avg_geodist / real(n-1)
            avg_euldist = avg_euldist / real(n-1)
            write(logfhandle,'(A,F8.3)')'>>> AVERAGE EULER    DISTANCE: ',avg_euldist
            write(logfhandle,'(A,F8.3)')'>>> AVERAGE GEODESIC DISTANCE: ',avg_geodist
        endif
        call o%kill
        call o_prev%kill
        call simple_end('**** VIZORIS NORMAL STOP ****')
    end subroutine exec_vizoris

    subroutine check_states( self, cline )
        class(commander_check_states), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(parameters)     :: params
        type(builder)        :: build
        type(oris)           :: o_truth
        real,    allocatable :: states(:), truth_states(:), diluted(:,:)
        integer, allocatable :: truth_nptcls(:), correct_states(:), state_order(:), cur_perm(:), all_perms(:,:), est_states(:)
        integer              :: nptcls, nstates, iptcl, istate, perm_cnt, max_perm, iperm, truth_state, est_state
        real                 :: cur_sum, max_sum
        call build%init_params_and_build_general_tbox(cline, params, do3d=(trim(params%oritype) .eq. 'ptcl3D'))
        nptcls  = build%spproj_field%get_noris()
        nstates = build%spproj_field%get_n('state')
        allocate(states(nptcls), truth_states(nptcls), diluted(nstates, nstates), source=0.)
        allocate(truth_nptcls(nstates), correct_states(nstates), state_order(nstates), cur_perm(nstates),&
                &all_perms(nstates, 2**nstates), est_states(nstates), source=0)
        states = build%spproj_field%get_all('state')
        call o_truth%new(nptcls, .true.)
        call o_truth%read(params%oritab, [1, nptcls])
        truth_states = o_truth%get_all('state')
        do istate = 1, nstates
            do iptcl = 1, nptcls
                if( int(truth_states(iptcl)) == istate ) truth_nptcls(istate) = truth_nptcls(istate) + 1
            enddo
        enddo
        if( cline%defined('infile') )then
            call read_state_order(params%infile, state_order)
        else
            perm_cnt = 0
            call generate_perm(1)
            max_sum  = 0.
            max_perm = 1
            do iperm = 1, perm_cnt
                state_order    = all_perms(:, iperm)
                correct_states = 0
                do iptcl = 1, nptcls
                    istate = int(truth_states(iptcl))
                    if( int(states(iptcl)) == state_order(istate) ) correct_states(istate) = correct_states(istate) + 1
                enddo
                cur_sum = sum(correct_states * 100. / real(truth_nptcls))
                if( cur_sum > max_sum )then
                    max_perm = iperm
                    max_sum  = cur_sum
                endif
            enddo
            state_order = all_perms(:, max_perm)
            print *, 'BEST STATE PERMUTATION:'
            do istate = 1, nstates
                print *, 'Reconstructed vol ', int2str(state_order(istate)), ' looks mostly like Truth volume ', int2str(istate)
            enddo
        endif
        ! counting est_states
        est_states = 0
        diluted    = 0.
        do istate = 1, nstates
            est_state = state_order(istate)
            do iptcl = 1, nptcls
                if( int(states(iptcl)) == est_state )then
                    est_states(istate) = est_states(istate) + 1
                    truth_state = int(truth_states(iptcl))
                    diluted(truth_state, est_state) = diluted(truth_state, est_state) + 1
                endif
            enddo
        enddo
        print *, 'TRUTH COMPOSITION TABLE: (from Truth vol 1 to Truth vol ', int2str(nstates), ') '
        do istate = 1, nstates
            print *, 'Reconstructed vol ', int2str(state_order(istate)), ': ', (diluted(:, state_order(istate)) * 100. / real(truth_nptcls(istate)))
        enddo
        print *, 'TRUTH COMPOSITION TABLE (ref paper metrics): (from Truth vol 1 to Truth vol ', int2str(nstates), ') '
        do istate = 1, nstates
            print *, 'Reconstructed vol ', int2str(state_order(istate)), ': ', (diluted(:, state_order(istate)) * 100. / real(est_states(istate)))
        enddo
        ! cleanup
        call build%kill_general_tbox
        ! end gracefully
        call simple_end('**** SIMPLE_CHECK_STATES NORMAL STOP ****')

      contains

        subroutine read_state_order( fname, order )
            class(string), intent(in)  :: fname
            integer,       intent(out) :: order(nstates)
            character(len=100)  :: io_message
            type(string) :: line
            integer :: file_stat, fnr, i, io_stat
            if( .not. file_exists(fname) )then
                THROW_HARD("the file you are trying to read: "//fname%to_char()//' does not exist in cwd' )
            endif
            if( fname2ext(fname) == 'bin' )then
                THROW_HARD('this method does not support binary files; read')
            endif
            io_message='No error'
            call fopen(fnr, FILE=fname, STATUS='OLD', action='READ', iostat=file_stat, iomsg=io_message)
            call fileiochk("read_state_order ; read ,Error when opening file for reading: "//fname%to_char()//':'//trim(io_message), file_stat)
            do i = 1, nstates
                call line%readline(fnr, io_stat)
                order(i) = line%to_int()
            enddo
            call fclose(fnr)
        end subroutine read_state_order

        ! generate all permutations of states
        recursive subroutine generate_perm( position )
            integer, intent(in) :: position
            integer :: value
            if( position > nstates )then
                perm_cnt              = perm_cnt + 1
                all_perms(:,perm_cnt) = cur_perm
            else
                do value = 1, nstates
                    if( .not. any(cur_perm(:position-1) == value) )then
                        cur_perm(position) = value
                        call generate_perm(position + 1)
                    endif
                enddo
            endif
        end subroutine generate_perm

    end subroutine check_states

end module simple_commanders_oris
