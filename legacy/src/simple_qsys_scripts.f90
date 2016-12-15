module simple_qsys_scripts
use simple_jiffys     ! singleton
use simple_defs       ! singleton
use simple_qsys_base, only: qsys_base
implicit none

public :: init_qsys_scripts, gen_qsys_scripts, submit_qsys_scripts, qsys_scripts_are_finished
public :: qsys_cleanup_iter, qsys_cleanup_all, qsys_stack_is_split
private

character(len=:), allocatable  :: exec_cmd         !< command to execute in parallel
class(qsys_base), pointer      :: myqsys=>null()   !< pointer to polymorphic qsys object
integer                        :: nptcls = 0       !< number of particles
integer                        :: npart  = 0       !< number of partitions
integer                        :: numlen = 0       !< length of padded number string
integer, allocatable           :: parts(:,:)       !< defines the partitions (fromp:top)
character(len=STDLEN)          :: pwd              !< working directory
character(len=32), allocatable :: script_names(:)  !< file-names of generated scripts

contains
    
    ! INITIALISATION

    subroutine init_qsys_scripts( exec_cmd_in, qsys_obj, nptcls_in, npart_or_chunksz, split_mode )
        use simple_map_reduce ! singleton
        character(len=*), intent(in)         :: exec_cmd_in
        class(qsys_base), intent(in), target :: qsys_obj
        integer, intent(in)                  :: nptcls_in, npart_or_chunksz
        character(len=*), intent(in)         :: split_mode
        call kill_qsys_script_gen
        allocate(exec_cmd, source=exec_cmd_in)
        allocate(script_names(npart_or_chunksz))
        myqsys => qsys_obj
        nptcls = nptcls_in
        numlen = len(int2str(npart))
        select case(split_mode)
            case('chunk')
                parts = split_nobjs_in_chunks(nptcls, npart_or_chunksz)
            case('even')
                parts = split_nobjs_even(nptcls, npart_or_chunksz)
            case DEFAULT
                write(*,*) 'split mode: ', trim(split_mode)
                stop 'unsupported split_mode; simple_qsys_scripts :: init_qsys_scripts'
        end select
        npart = size(parts,1)
        call get_environment_variable('PWD', pwd)
    end subroutine init_qsys_scripts
    
    ! SCRIPT GENERATION
    
    subroutine gen_qsys_scripts( job_descr, outfile_body )
        use simple_chash, only: chash
        class(chash),               intent(inout) :: job_descr
        character(len=*), optional, intent(in)    :: outfile_body
        character(len=:), allocatable :: outfile_body_local
        integer :: ipart
        if( present(outfile_body) )then
            allocate(outfile_body_local, source=trim(outfile_body))
        else
            allocate(outfile_body_local, source='algndoc_')
        endif
        do ipart=1,npart
            call job_descr%set('fromp',   int2str(parts(ipart,1)))
            call job_descr%set('top',     int2str(parts(ipart,2)))
            call job_descr%set('part',    int2str(ipart))
            call job_descr%set('outfile', trim(outfile_body_local)//int2str_pad(ipart,numlen)//'.txt')
            script_names(ipart) = gen_qsys_script(job_descr, ipart)
        end do
        call job_descr%delete('fromp')
        call job_descr%delete('top')
        call job_descr%delete('part')
        call job_descr%delete('outfile')
        deallocate(outfile_body_local)
    end subroutine gen_qsys_scripts
    
    function gen_qsys_script( job_descr, ipart ) result( fname )
        use simple_chash, only: chash
        class(chash),     intent(in)  :: job_descr
        integer,          intent(in)  :: ipart
        character(len=32)  :: keys2print(4), fname
        character(len=512) :: io_msg
        integer :: ios, funit
        funit = get_fileunit()
        fname = 'distr_simple_script_'//int2str_pad(ipart,numlen)
        open(unit=funit, file=fname, iostat=ios, status='replace', iomsg=io_msg)
        if( ios .ne. 0 )then
            close(funit)
            write(*,'(a)') 'simple_qsys_scripts :: gen_qsys_script; Error when opening file for writing: '&
            //trim(fname)//' ; '//trim(io_msg)
            stop
        endif
        ! need to specify shell
        write(funit,'(a)') '#!/bin/bash'
        ! write (run-time polymorphic) instructions to the qsys
        call myqsys%write_instr(job_descr, fhandle=funit)
        write(funit,'(a)') 'cd '//trim(pwd)
        keys2print(1) = 'fromp'
        keys2print(2) = 'top'
        keys2print(3) = 'part'
        keys2print(4) = 'outfile'
        ! compose the command line
        write(funit,'(a)',advance='no') trim(exec_cmd)//' '
        call job_descr%print_key_val_pairs(keys2print, fhandle=funit, oneline=.true.)
        ! direct output
        write(funit,'(a)') '> OUT'//int2str_pad(ipart,numlen)
        ! exit shell when done
        write(funit,'(a)') 'exit'
        close(funit)
        call chmod(fname,'+x',ios)
        if( ios .ne. 0 )then
            write(*,'(a)') 'simple_qsys_scripts :: gen_qsys_script; Error chmoding submit script'//trim(fname)
            stop
        endif
    end function gen_qsys_script
    
    ! SUBMISSION TO QSYS
    
    subroutine submit_qsys_scripts
        use simple_qsys_local, only: qsys_local
        integer                       :: cstat, estat, ipart
        character(len=100)            :: cmsg
        character(len=:), allocatable :: script_exec_cmd
        do ipart=1,npart
            select type(myqsys)
                class is(qsys_local)
                    allocate(script_exec_cmd, source=myqsys%submit_cmd()//' ./'//trim(script_names(ipart))//' &')
                class DEFAULT
                    allocate(script_exec_cmd, source=myqsys%submit_cmd()//' ./'//trim(script_names(ipart)))
            end select
            call execute_command_line(script_exec_cmd, exitstat=estat, cmdstat=cstat, cmdmsg=cmsg)
            if( cstat > 0 )then
                print *, 'submit_qsys_scripts; command execution failed with error ', trim(cmsg)
            elseif( cstat < 0 )then
                print *, 'submit_qsys_scripts; command execution not supported'
            endif
            deallocate(script_exec_cmd)
        end do
    end subroutine submit_qsys_scripts
    
    ! QUERIES
    
    function qsys_scripts_are_finished() result( are_finished )
        character(len=:), allocatable :: job_done_fname
        integer :: ipart
        logical :: jobs_done(npart), are_finished
        do ipart=1,npart
            allocate(job_done_fname, source='JOB_FINISHED_'//int2str_pad(ipart,numlen))
            jobs_done(ipart) = file_exists(job_done_fname)
            deallocate(job_done_fname)
        end do
        are_finished = all(jobs_done)
    end function qsys_scripts_are_finished
    
    function qsys_stack_is_split( stkext, num_parts ) result( is_split )
        character(len=4),  intent(in) :: stkext
        integer, optional, intent(in) :: num_parts
        character(len=:), allocatable :: stack_part_fname
        logical, allocatable          :: stack_parts_exist(:) 
        integer :: inpart, inumlen, ipart
        logical :: is_split
        if( present(num_parts) )then
            inpart  = num_parts
            inumlen = len(int2str(inpart))
        else
            inpart  = npart
            inumlen = numlen
        endif
        allocate( stack_parts_exist(inpart) )
        do ipart=1,inpart
            allocate(stack_part_fname, source='stack_part'//int2str_pad(ipart,inumlen)//stkext)
            stack_parts_exist(ipart) = file_exists(stack_part_fname)
            deallocate(stack_part_fname)
        end do
        is_split = all(stack_parts_exist)
    end function qsys_stack_is_split
    
    ! CLEANUP
    
    !> \brief  is for deleting EVERYTHING generated. Use with EXTREME caution
    subroutine qsys_cleanup_all
        integer                       :: cstat, estat
        character(len=100)            :: cmsg
        character(len=:), allocatable :: cleanup_exec_cmd, tmp1, tmp2, tmp3
        allocate(tmp1,source='rm -rf FOO OUT* algndoc_* distr_script_* rho* fort.0 primedoc_* recvol*' )
        allocate(tmp2,source=' JOB_FINISHED_* errfile.* outfile.* shdistr_script SHMEMJOBOUT shmemerrfile.* shmemoutfile.*')
        allocate(tmp3,source=' ctfsqsums_part* noisespecs_part*')
        allocate(cleanup_exec_cmd, source=tmp1//tmp2//tmp3)
        call execute_command_line(cleanup_exec_cmd, exitstat=estat, cmdstat=cstat, cmdmsg=cmsg)
        if( cstat > 0 )then
            print *, 'qsys_cleanup_all; command execution failed with error ', trim(cmsg)
        elseif( cstat < 0 )then
            print *, 'qsys_cleanup_all; command execution not supported'
        endif
        deallocate(cleanup_exec_cmd, tmp1, tmp2, tmp3)
    end subroutine qsys_cleanup_all

    !> \brief  is for cleaning up after an iteration
    subroutine qsys_cleanup_iter
        integer                       :: cstat, estat
        character(len=100)            :: cmsg
        character(len=:), allocatable :: cleanup_exec_cmd, tmp1, tmp2, tmp3
        allocate(tmp1,source='rm -rf FOO OUT* algndoc_* distr_*script_* recvol_state*_part* rho* fort.0' )
        allocate(tmp2,source=' JOB_FINISHED_* errfile.* outfile.* shdistr_script SHMEMJOBOUT shmemerrfile.* shmemoutfile.*')
        allocate(tmp3,source=' ctfsqsums_part* noisespecs_part* cavgs_part*')
        allocate(cleanup_exec_cmd, source=tmp1//tmp2//tmp3)
        call execute_command_line(cleanup_exec_cmd, exitstat=estat, cmdstat=cstat, cmdmsg=cmsg)
        if( cstat > 0 )then
            print *, 'qsys_cleanup_iter; command execution failed with error ', trim(cmsg)
        elseif( cstat < 0 )then
            print *, 'qsys_cleanup_iter; command execution not supported'
        endif
        deallocate(cleanup_exec_cmd, tmp1, tmp2, tmp3)
    end subroutine qsys_cleanup_iter
        
    ! DESTRUCTOR
    
    subroutine kill_qsys_script_gen
        if( allocated(exec_cmd))     deallocate(exec_cmd)
        if( allocated(parts))        deallocate(parts)
        if( allocated(script_names)) deallocate(script_names)
        myqsys=>null()
        nptcls=0 
        numlen=0
    end subroutine kill_qsys_script_gen

end module simple_qsys_scripts
