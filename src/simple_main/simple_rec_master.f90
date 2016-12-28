module simple_rec_master
use simple_build,   only: build
use simple_params,  only: params
use simple_cmdline, only: cmdline
use simple_jiffys   ! singleton
use simple_defs     ! singleton
use simple_math     ! singleton
use simple_filterer ! singleton
implicit none

public :: exec_rec_master
private

logical :: debug=.false.

contains

    subroutine exec_rec_master( b, p, cline, fbody_in, wmat )
        class(build),               intent(inout) :: b
        class(params),              intent(inout) :: p
        class(cmdline),             intent(inout) :: cline
        character(len=*), optional, intent(in)    :: fbody_in
        real,             optional, intent(in)    :: wmat(:,:)
        select case(p%eo)
            case( 'yes' )
                call exec_eorec( b, p, cline, fbody_in, wmat )
            case( 'no' )
                call exec_rec( b, p, cline, fbody_in, wmat )
            case DEFAULT
                stop 'unknonw eo flag; simple_rec_master :: exec_rec_master'
        end select
    end subroutine exec_rec_master
    
    subroutine exec_rec( b, p, cline, fbody_in, wmat )
        class(build),               intent(inout) :: b
        class(params),              intent(inout) :: p
        class(cmdline),             intent(inout) :: cline
        character(len=*), optional, intent(in)    :: fbody_in
        real,             optional, intent(in)    :: wmat(:,:)
        character(len=:), allocatable :: fbody
        integer :: s, fri, toi, file_stat, fnr, nstates_oritab
        ! rebuild b%vol according to box size (beacuse it is otherwise boxmatch)
        call b%vol%new([p%box,p%box,p%box], p%smpd, p%imgkind)
        if( p%mirr .eq. 'yes' ) call b%a%mirror3d
        if( cline%defined('state') )then ! setting iteration from/to state
            fri = p%state
            toi = p%state
        else
            fri = 1
            toi = p%nstates
        endif
        do s=fri,toi ! state loop
            if( debug ) write(*,*) 'processing state: ', s
            if( p%l_distr_exec )then ! embarrasingly parallel rec
                if( present(fbody_in) )then
                    allocate(fbody, source=trim(adjustl(fbody_in))//&
                    &'_state'//int2str_pad(s,2)//'_part'//int2str_pad(p%part,p%numlen))
                else
                    allocate(fbody, source='recvol_state'//int2str_pad(s,2)//&
                    &'_part'//int2str_pad(p%part,p%numlen))
                endif
                if( p%even .eq. 'yes' .and. p%odd .eq. 'no' )then
                    p%vols(s) = fbody//'_even'//p%ext
                    p%masks(s) = 'rho_'//fbody//'_even'//p%ext
                    call b%recvol%rec(p%stk, p, b%a, b%se, s, mul=p%mul, eo=2, part=p%part, wmat=wmat)
                else if( p%odd .eq. 'yes' .and. p%even .eq. 'no' )then
                    p%vols(s) = fbody//'_odd'//p%ext
                    p%masks(s) = 'rho_'//fbody//'_odd'//p%ext
                    call b%recvol%rec(p%stk, p, b%a, b%se, s, mul=p%mul, eo=1, part=p%part, wmat=wmat)
                else if( p%even .eq. 'yes' .and. p%odd .eq. 'yes' )then
                    stop 'ERROR! even and odd cannot both be yes!'
                else
                    p%vols(s)  = fbody//p%ext
                    p%masks(s) = 'rho_'//fbody//p%ext
                    call b%recvol%rec(p%stk, p, b%a, b%se, s, mul=p%mul, part=p%part, wmat=wmat)
                endif
                call b%recvol%write(p%vols(s), del_if_exists=.true.)
                call b%recvol%write_rho(p%masks(s))
            else ! shared-mem parallel rec
                if( present(fbody_in) )then
                    allocate(fbody, source=trim(adjustl(fbody_in))//'_state')
                else
                    allocate(fbody, source='recvol_state')
                endif
                p%vols(s) = fbody//int2str_pad(s,2)//p%ext ! shared mem parallel rec
                call b%recvol%rec(p%stk, p, b%a, b%se, s, mul=p%mul, wmat=wmat)
                call b%recvol%clip(b%vol)
                call b%vol%write(p%vols(s), del_if_exists=.true.)
            endif
            deallocate(fbody)
        end do
        write(*,'(a)') "GENERATED VOLUMES: recvol*.ext"
        if( p%l_distr_exec )then
            fnr = get_fileunit()
            open(unit=fnr, FILE='JOB_FINISHED_'//int2str_pad(p%part,p%numlen),&
            &STATUS='REPLACE', action='WRITE', iostat=file_stat)
            call fopen_err( 'In: simple_recvol', file_stat )
            close(unit=fnr)
        endif
    end subroutine exec_rec
    
    subroutine exec_eorec( b, p, cline, fbody_in, wmat )
        class(build),               intent(inout) :: b
        class(params),              intent(inout) :: p
        class(cmdline),             intent(inout) :: cline    
        character(len=*), optional, intent(in)    :: fbody_in
        real,             optional, intent(in)    :: wmat(:,:)
        character(len=:), allocatable :: fbody, fname
        integer                       :: s, fri, toi, fnr, file_stat, ldim(3)
        real                          :: smpd
        ! rebuild b%vol according to box size (beacuse it is otherwise boxmatch)
        call b%vol%new([p%box,p%box,p%box], p%smpd, p%imgkind)
        if( cline%defined('state') )then ! setting iteration from/to state
            fri = p%state
            toi = p%state
        else
            fri = 1
            toi = p%nstates
        endif
        do s=fri,toi ! state loop
            if( debug ) write(*,*) 'processing state: ', s
            if( p%l_distr_exec )then ! embarrasingly parallel exec
                if( present(fbody_in) )then
                    allocate(fbody, source=trim(adjustl(fbody_in))//'_state')
                else
                    allocate(fbody, source='recvol_state')
                endif
                call b%eorecvol%eorec(p%stk, p, b%a, b%se, s, b%vol, mul=p%mul, part=p%part, fbody=fbody, wmat=wmat)
            else
                if( present(fbody_in) )then
                    allocate( fbody, source=trim(adjustl(fbody_in))//'_state' )
                else
                    allocate( fbody, source='recvol_state' )
                endif                
                call b%eorecvol%eorec(p%stk, p, b%a, b%se, s, b%vol, mul=p%mul, wmat=wmat)
                allocate(fname, source=fbody//int2str_pad(s,2)//p%ext)
                call b%vol%write(fname, del_if_exists=.true.)
                deallocate(fname)
            endif
            deallocate(fbody)
        end do
        if( p%l_distr_exec )then
            fnr = get_fileunit()
            open(unit=fnr, FILE='JOB_FINISHED_'//int2str_pad(p%part,p%numlen),&
            &STATUS='REPLACE', action='WRITE', iostat=file_stat)
            call fopen_err( 'In: simple_recvol', file_stat )
            close(unit=fnr)
        endif   
        write(*,'(a,1x,a)') "GENERATED VOLUMES: recvol*.ext"
    end subroutine exec_eorec
    
end module simple_rec_master
