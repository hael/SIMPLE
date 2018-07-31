program simple_test_star_export
include 'simple_lib.f08'
use, intrinsic :: iso_fortran_env
use, intrinsic :: iso_c_binding
use simple_star
use simple_cmdline,        only: cmdline
use simple_parameters
!use simple_builder,          only: builder
use simple_sp_project,     only: sp_project
use simple_oris,       only: oris
use simple_sp_project, only: sp_project
use simple_binoris,    only: binoris
use simple_stardoc
use simple_star_dict,  only:  star_dict
implicit none
#include "simple_local_flags.inc"
character(len=STDLEN) :: oldCWDfolder,curDir,datestr
character(len=STDLEN) :: simple_testbench,timestr,folder, fileref, filecompare
character(len=:),allocatable:: tmpfile,line
integer(8) :: count1
integer :: io_stat,funit,ier
integer :: num_data_elements, num_data_lines
type(sp_project)     :: myproject
type(binoris)        :: bos
integer, allocatable :: strlens(:)
type(star_project) :: s
type(parameters) :: p
type(stardoc) :: sdoc
type(star_dict) :: sdict
logical :: isopened
CHARACTER(len=STDLEN) :: argtmp

if( command_argument_count() == 1 )then
    call get_command_argument(1,argtmp)
    print *, trim(argtmp)
    allocate(tmpfile,source=trim(argtmp))
else
    allocate(tmpfile,source="/scratch/el85/stars_from_matt/micrographs_all_gctf.star")
!tmpfile="/scratch/el85/stars_from_matt/micrographs_all_gctf.star"
!io_stat = simple_getenv('SIMPLE_TESTBENCH_DATA', simple_testbench)
!tmpfile=filepath(trim(simple_testbench),"star_test/stars_from_matt/micrographs_all_gctf.star")
endif

print *, trim(tmpfile)
global_debug=.true.
debug=.true.
isopened=.false.
call date_and_time(date=datestr)
folder = 'SIMPLE_TEST_STAR_'//trim(datestr)
call simple_mkdir( trim(folder) , status=io_stat)
if(io_stat/=0) call simple_stop("simple_mkdir failed")
print *," Changing directory to ", folder
call simple_chdir( trim(folder),  oldCWDfolder , status=io_stat)
if(io_stat/=0) call simple_stop("simple_chdir failed")
call simple_getcwd(curDir)
print *," Current working directory ", curDir
count1=tic()

 print *,' Testing star_dict'
 if( .not. sdict%exists()) then
     print *,"Star dictionary constructor did not call new() "
     call sdict%new()
 endif

 call sdict%print_star_dict()



print *,"(info)**simple_test_star_export:  Testing star formatted file ", trim(tmpfile)
!tmpfile= trim(adjustl(tmpfile))
!! Testing starformat
call openstar(trim(tmpfile),funit)
print *,"(info)**simple_test_star_export:  Testing star formatted file is_open ", is_open(funit)
if(.not. is_open(funit)) then
    HALT_NOW("readline isopened failed")
endif
 ier=0
 print *,"(info)**simple_test_star_export: Readline testing"
 do while (ier==0)
     call readline(funit,line,ier)
     if(allocated(line))print *,line
 enddo
 print *,"(info)**simple_test_star_export: Reading header"
 call read_header(funit)
! write(*,*)
! write(*,*)
! write(*,*)
! write(*,*)
! print *,"(info)**simple_test_star_export: Testing star formatted file, ", funit, is_open(funit)
! print *,"(info)**simple_test_star_export: Reading body"
! call read_data_lines(funit)
! print *,"(info)**simple_test_star_export: Testing star formatted file, ", funit, is_open(funit)
 call fclose(funit,io_stat,errmsg='stardoc ; close ')
 if( is_open(funit)) then
     HALT_NOW(" star file still open after closing")
 endif

  print *, 'Testing star module'
  call s%prepareimport(myproject,p,tmpfile)
 ! ! call s%import_ctf_estimation(myproject,tmpfile)
  call s%kill(keeptabs=.true.)
 ! print *, '     star module imported successfully'

!  call exec_cmdline("simple_exec prg=new_project projname=test && cd test && "//&
!      &"simple_exec prg=import_movies filetab=../filetab-stardoc.txt  cs=2.7 ctf=yes "//&
!      &"fraca=0.1 kv=300 smpd=14 deftab=../oritab-stardoc.txt && "//&
!      &"simple_exec prg=print_project_info")




! call test_stardoc

print *,' Testing directory star_test'
!call system('ls '//trim(adjustl(simple_testbench))//PATH_SEPARATOR//'star_test')

!! Motion Correction
! call s%export_motion_corrected_micrographs (trim('tmp_mc.star'))

! call create_relion_starfile
! call sp_project_setup

contains

    subroutine create_relion_starfile
        call exec_cmdline( 'relion_star_loopheader rlnMicrographNameNoDW rlnMicrographName > tmp_mc.star')

        ! Generate STAR files from separate stacks for each micrograph
        ! If the input images are in a separate stack for each micrograph, then one could use the following commands to generate the input STAR file:
        call exec_cmdline( ' relion_star_loopheader rlnImageName rlnMicrographName rlnDefocusU rlnDefocusV rlnDefocusAngle rlnVoltage rlnSphericalAberration rlnAmplitudeContrast > my_images.star')
        call exec_cmdline( ' relion_star_datablock_stack 4 mic1.mrcs mic1.mrcs 10000 10500 30 200 2 0.1  >> my_images.star')
        call exec_cmdline( 'relion_star_datablock_stack 3 mic2.mrcs mic2.mrcs 21000 20500 25 200 2 0.1  >> my_images.star')
        call exec_cmdline( ' relion_star_datablock_stack 2 mic3.mrcs mic3.mrcs 16000 15000 35 200 2 0.1  >> my_images.star')



        ! Generate STAR files from particles in single-file format
        !call exec_cmdline( ' relion_star_loopheader rlnImageName rlnMicrographName rlnDefocusU rlnDefocusV rlnDefocusAngle rlnVoltage rlnSphericalAberration rlnAmplitudeContrast > my_images.star
        ! call exec_cmdline( 'relion_star_datablock_singlefiles "mic1/*.spi" mic1 16000 15000 35 200 2 0.1  >> my_images.star
        !call exec_cmdline( ' relion_star_datablock_singlefiles "mic2/*.spi" mic2 16000 15000 35 200 2 0.1  >> my_images.star
        !call exec_cmdline( ' relion_star_datablock_singlefiles "mic3/*.spi" mic3 16000 15000 35 200 2 0.1  >> my_images.star

        ! Generate STAR files from XMIPP-style CTFDAT files
        ! To generate a STAR file from an XMIPP-style ctfdat file, one could use:

        !call exec_cmdline( ' relion_star_loopheader rlnImageName rlnMicrographName rlnDefocusU rlnDefocusV rlnDefocusAngle rlnVoltage rlnSphericalAberration rlnAmplitudeContrast > all_images.star
        !call exec_cmdline( ' relion_star_datablock_ctfdat all_images.ctfdat>>  all_images.star

        ! Generate STAR files from FREALIGN-style .par files
        ! To generate a STAR file from a FREALIGN-style .par file, one could use:

        ! call exec_cmdline( 'relion_star_loopheader rlnImageName rlnMicrographName rlnDefocusU rlnDefocusV rlnDefocusAngle rlnVoltage rlnSphericalAberration rlnAmplitudeContrast > all_images.star
        ! call exec_cmdline( 'awk '{if ($1!="C") {print $1"@./my/abs/path/bigstack.mrcs", $8, $9, $10, $11, " 80 2.0 0.1"}  }' < frealign.par >> all_images.star ')
        ! Assuming the voltage is 80kV, the spherical aberration is 2.0 and the amplitude contrast is 0.1. Also, a single stack is assumed called: /my/abs/path/bigstack.mrcs.
    end subroutine create_relion_starfile



    subroutine sp_project_setup

        ! prepare stack oris in project
        call myproject%os_stk%new(2)
        ! motion_correct
        call myproject%os_stk%set(1, 'movie',       'movie1.mrc')
        call myproject%os_stk%set(1, 'intg',        'movie1_intg.mrc')
        call myproject%os_stk%set(1, 'forctf',      'movie1_forctf.mrc')
        call myproject%os_stk%set(1, 'pspec',       'movie1_pspec.mrc')
        call myproject%os_stk%set(1, 'thumb',       'movie1_thumb.mrc')
        call myproject%os_stk%set(1, 'intg_frames', 'movie1_intg_frames1-14.mrc')
        call myproject%os_stk%set(1, 'smpd',        1.3)
        call myproject%os_stk%set(2, 'movie',       'movie2.mrc')
        call myproject%os_stk%set(2, 'intg',        'movie2_intg.mrc')
        call myproject%os_stk%set(2, 'forctf',      'movie2_forctf.mrc')
        call myproject%os_stk%set(2, 'pspec',       'movie2_pspec.mrc')
        call myproject%os_stk%set(2, 'thumb',       'movie2_thumb.mrc')
        call myproject%os_stk%set(2, 'intg_frames', 'movie2_intg_frames1-14.mrc')
        call myproject%os_stk%set(2, 'smpd',        1.3)
        ! write/read
        call myproject%write('myproject.simple')
        call myproject%read('myproject.simple')
        ! call myproject%print_header
        call myproject%write_segment2txt('stk', 'myproject_os_stk_1.txt')
        ! ctf_estimate
        call myproject%os_stk%set(1, 'kv',       300.)
        call myproject%os_stk%set(1, 'cs',       2.7 )
        call myproject%os_stk%set(1, 'fraca',    0.1 )
        call myproject%os_stk%set(1, 'dfx',      1.2 )
        call myproject%os_stk%set(1, 'dfy',      1.3 )
        call myproject%os_stk%set(1, 'angast',   30. )
        call myproject%os_stk%set(1, 'phshift',  0.  )
        call myproject%os_stk%set(1, 'ctf_estimatecc', 0.8 )
        call myproject%os_stk%set(1, 'ctfres',   5.2 )
        call myproject%os_stk%set(2, 'kv',       300.)
        call myproject%os_stk%set(2, 'cs',       2.7 )
        call myproject%os_stk%set(2, 'fraca',    0.1 )
        call myproject%os_stk%set(2, 'dfx',      1.5 )
        call myproject%os_stk%set(2, 'dfy',      1.7 )
        call myproject%os_stk%set(2, 'angast',   60. )
        call myproject%os_stk%set(2, 'phshift',  0.  )
        call myproject%os_stk%set(2, 'ctf_estimatecc', 0.75)
        call myproject%os_stk%set(2, 'ctfres',   3.8 )
        ! write/read
        call myproject%write('myproject.simple')
        call myproject%read('myproject.simple')
        call myproject%print_info
        call myproject%write_segment2txt('stk', 'myproject_os_stk_2.txt')

        ! generate 3 algndocs for testing merging
        call myproject%os_ptcl3D%new(9)
        call myproject%os_ptcl3D%set_euler(1, [1.,1.,1.])
        call myproject%os_ptcl3D%set_euler(2, [1.,1.,1.])
        call myproject%os_ptcl3D%set_euler(3, [1.,1.,1.])
        print *, 'writing doc1'
        call myproject%write('doc1.simple', [1,3])
        call myproject%os_ptcl3D%new(9)
        call myproject%os_ptcl3D%set_euler(4, [2.,2.,2.])
        call myproject%os_ptcl3D%set_euler(5, [2.,2.,2.])
        call myproject%os_ptcl3D%set_euler(6, [2.,2.,2.])
        print *, 'writing doc2'
        call myproject%write('doc2.simple', [4,6])
        call myproject%os_ptcl3D%new(9)
        call myproject%os_ptcl3D%set_euler(7, [3.,3.,3.])
        call myproject%os_ptcl3D%set_euler(8, [3.,3.,3.])
        call myproject%os_ptcl3D%set_euler(9, [3.,3.,3.])
        print *, 'writing doc3'
        call myproject%write('doc3.simple', [7,9])


    end subroutine sp_project_setup

    subroutine openstar( filename, funit)
        character(len=*),intent(in) :: filename
        integer,intent(out)::funit
        integer :: io_stat, tmpunit
        integer(8):: filesz
        if(.not. file_exists(trim(filename) ))&
            call simple_stop("simple_stardoc::open ERROR file does not exist "//trim(filename) )
        call fopen(tmpunit, file=trim(filename), action='READ',&
            & iostat=io_stat)
        if(io_stat/=0)call fileiochk('star_doc ; open '//trim(filename), io_stat)
        funit  = tmpunit

        ! check size
        filesz = funit_size(funit)
        if( filesz == -1 )then
            stop 'file_size cannot be inquired; stardoc :: open'
        else if (filesz < 10) then
            write(*,*) 'file: ', trim(filename)
            stop 'file size too small to contain a header; stardoc :: open'
        endif
    end subroutine openstar




    subroutine readline(funit, line, ier)
        use iso_fortran_env
        implicit none
        integer, intent(in)                      :: funit
        character(len=:),allocatable,intent(out) :: line
        integer,intent(out)                      :: ier

        integer,parameter                     :: buflen=1024
        character(len=buflen)                 :: buffer
        integer                               :: last
        integer                               :: isize
        logical :: isopened
        line=''
        ier=0
        inquire(unit=funit,opened=isopened,iostat=ier)
        if(ier/=0) call fileiochk("readline isopened failed", ier)
        if(.not. isopened )then
            HALT_NOW("readline isopened failed")
        endif
        ! read characters from line and append to result
        do
            ! read next buffer (an improvement might be to use stream I/O
            ! for files other than stdin so system line limit is not
            ! limiting)
            read(funit,fmt='(a)',advance='no',size=isize,iostat=ier) buffer
            ! append what was read to result
            !isize=len_trim(buffer)
            if(isize.gt.0)line=line//" "//buffer(:isize)
            ! if hit EOR reading is complete unless backslash ends the line
            !print *, ier == iostat_eor, ier == iostat_end, buffer
            !call fileiochk("readline isopened failed", ier)
           if(is_iostat_eor(ier))then
                last=len(line)
                ! if(last.ne.0)then
                !     ! if line ends in backslash it is assumed a continued line
                !     if(line(last:last).eq.'\')then
                !         ! remove backslash
                !         line=line(:last-1)
                !         ! continue on and read next line and append to result
                !         cycle
                !     endif
                ! endif
                ! hitting end of record is not an error for this routine
                ier=0
                ! end of reading line
                exit
                ! end of file or error
            elseif(ier.ne.0)then
                exit
            endif
        enddo
        call compact(line)  ! remove extra spaces and preceeding spaces

    end subroutine readline

    !> First pass of reading in STAR file
    subroutine read_header(funit)
        integer,intent(inout)  :: funit
        integer          :: cnt,n,ios,lenstr,isize,pos1,pos2, nargsline
        character(len=:), allocatable :: line, starlabel, eqvsimplelabel
        !character(len=:),allocatable :: argline(:)
        character(len=STDLEN) :: argline(128)
        logical :: inData, inHeader
        num_data_elements=0
        num_data_lines=0
        inData=.false.;inHeader=.false.
        ios=0;cnt=1
        !! Make sure we are at the start of the file
        rewind( funit, IOSTAT=ios)
        if(ios/=0)call fileiochk('star_doc ; read_header - rewind failed ', ios)
        ios=0
        do
            call readline(funit, line, ios)
            if(ios /= 0) exit
            line = trim(adjustl(line))
            print*, "STAR>> line #",cnt,":",line
            cnt=cnt+1
            !! Parse the start of the STAR file
            lenstr=len_trim(line)
            if (lenstr == 0 )cycle ! empty line

            !! Count number of fields in header
            if(inHeader)then
                if ( .not. (index(trim(line), "_rln") == 0) )then
                    num_data_elements=num_data_elements+1
                    pos1 = firstNonBlank(trim(line))
                    pos2 = firstBlank(trim(line))
                    DebugPrint " Found STAR field line ", trim(line), " VAR= ",pos2, line(pos1+4:pos2)
                    starlabel = line(pos1+4:pos2)
                    if( sdict%isthere(starlabel) )then
                        print *," STAR field found in dictionary ", num_data_elements
                        eqvsimplelabel = sdict%star2simple(starlabel)
                        print *," STAR field: ", trim(starlabel), " equivalent to ", trim(eqvsimplelabel)
                    else
                        print *," STAR field: ", trim(starlabel), " star field not in dictionary"
                    endif
                    cycle
                else
                    inHeader=.false.
                    DebugPrint " End of STAR data field lines "
                    inData = .true.
                endif
            endif

            !! Count number of data lines
            if(inData)then !! continue through to end of file
                num_data_lines=num_data_lines+1
                DebugPrint " Found STAR data line ", line
                !call parsestr(line,' ',argline,nargsline)
                nargsline = cntRecsPerLine(trim(line))
                if(nargsline /=  num_data_elements)then
                    print *, " Records on line mismatch ", nargsline,  num_data_elements
                    print *, "Line number ",num_data_lines, ":: ", line
                    stop " line has insufficient elements "
                endif
                cycle
            end if

            if (  .not. (index(trim(line), "data_") == 0))then !! does line contain ''data_''
                DebugPrint " Found STAR 'data_*' in header ", line
                !! Quick string length comparison
                if(lenstr == len_trim("data_pipeline_general"))then
                    print *," Found STAR 'data_pipeline_general' header -- Not supported "
                    exit
                else if(lenstr == len_trim("data_pipeline_processes"))then
                    print *," Found STAR 'data_pipeline_processes' header -- Not supported "
                    exit
                else if(lenstr == len_trim("data_sampling_general"))then
                    print *," Found STAR 'data_sampling_general' header -- Not supported "
                    exit
                else if(lenstr == len_trim("data_optimiser_general"))then
                    print *," Found STAR 'data_optimiser_general' header -- Not supported "
                    exit
                else if(lenstr == len_trim("data_model_general"))then
                    print *," Found STAR 'data_model_general' header -- Not supported "
                    exit
                end if
                !!otherwise
                cycle
            endif
            if ( .not. (index(trim(line), "loop_") == 0))then
                inHeader=.true.
                DebugPrint "Begin STAR field lines ", line
            endif

        end do
        rewind( funit,IOSTAT=ios)
        if(ios/=0)call fileiochk('star_doc ; read_header - rewind failed ', ios)
        DebugPrint " STAR field lines: ", num_data_elements
        DebugPrint " STAR data lines : ", num_data_lines

    end subroutine read_header

    !> Read STAR file and store information
    subroutine read_data_lines(funit)
        integer,intent(inout) :: funit
        integer          :: n,cnt,ios,lenstr, pos,i, nargsOnDataline, nDatalines
        character(len=:),allocatable :: line ! 8192
        logical :: inData, inHeader
        inHeader=.false.;inData=.false.
        n=1;cnt=1
        nDatalines=0
        write (*,*) " in Read data lines "
        do
            call readline (funit,line, ios)
            if(ios /= 0) exit
            line=trim(adjustl(line))
            print *,cnt,":",trim(line)
            cnt=cnt+1
            lenstr=len_trim(line)
            if (lenstr == 0 )cycle ! empty line

            !! Count number of fields in header
            if(inHeader)then
                if ( .not. (index(trim(line),"_rln") ==0) )then
                    ! ignore header
                    n = n+1
                else
                    inHeader=.false.
                    print *, " Beginning data lines "
                    inData = .true.
                endif
            endif
            if(inData)then
                nDatalines=nDatalines+1
                print *, " Found STAR data line ", line
                !call parsestr(line,' ',argline,nargsline)
                nargsOnDataline = cntRecsPerLine(trim(line))
                if(nargsOnDataline /=  num_data_elements) then
                    print *, " Records on line mismatch ", nargsOnDataline,  num_data_elements
                    print *, line
                    stop " line has insufficient elements "
                endif
                cycle
            endif

            !! Parse the start of the STAR file
            if ( .not. (index(trim(line), "data_") == 0)) cycle
            if ( .not. (index(trim(line) , "loop_") == 0))then
                inHeader=.true.
                print*, "Begin STAR field lines "
            endif
        end do
        if(nDatalines /= num_data_lines)then
            print *," Num data lines mismatch in read_data_lines and read_header"
            stop
        endif

        rewind( funit,IOSTAT=ios)
        if(ios/=0)call fileiochk('star_doc ; read_header - rewind failed ', ios)

    end subroutine read_data_lines



end program simple_test_star_export
