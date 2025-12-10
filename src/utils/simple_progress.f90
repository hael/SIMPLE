! job progress estimation
module simple_progress
include 'simple_lib.f08'
implicit none

contains

    subroutine progressfile_init()
        integer      :: progress_fhandle, ok
        type(string) :: fname
        fname = '.progress'
        if( file_exists(fname) ) call del_file(fname)
        call fopen(progress_fhandle,file=fname, status='new', iostat=ok)
        write(progress_fhandle, '(RD,F4.2)') 0.0
        call fclose(progress_fhandle)
        call fname%kill
    end subroutine progressfile_init
    
    subroutine progressfile_init_parts(nparts)
        integer, intent(in) :: nparts
        integer             :: ipart
        do ipart = 1, nparts
            call progressfile_init_part(ipart)
        enddo
    end subroutine progressfile_init_parts
    
    subroutine progressfile_init_part(part)
        integer, intent(in) :: part
        integer             :: progress_fhandle, ok
        type(string) :: fname
        fname = '.progress_' // trim(adjustl(int2str(part)))
        if( file_exists(fname) ) call del_file(fname)
        call fopen(progress_fhandle,file=fname, status='new', iostat=ok)
        write(progress_fhandle, '(RD,F4.2)') 0.0
        call fclose(progress_fhandle)
        call fname%kill
    end subroutine progressfile_init_part
    
    subroutine progressfile_update(progress)
        real, intent(in) :: progress
        integer          :: progress_fhandle, ok
        type(string) :: fname
        fname = '.progress'
        if( file_exists(fname) ) call del_file(fname)
        call fopen(progress_fhandle,file=fname, status='new', iostat=ok)
        if(progress > 1.0) then
            write(progress_fhandle, '(RD,F4.2)') 1.0
        else if(progress < 0.0) then
            write(progress_fhandle, '(RD,F4.2)') 0.0
        else
            write(progress_fhandle, '(RD,F4.2)') progress
        endif
        call fclose(progress_fhandle)
        call fname%kill
    end subroutine progressfile_update
    
    subroutine progressfile_update_part(part, progress)
        real,    intent(in) :: progress
        integer, intent(in) :: part
        integer             :: progress_fhandle, ok
        type(string)        :: fname
        fname = '.progress_' // trim(adjustl(int2str(part)))
        if( file_exists(fname) ) call del_file(fname)
        call fopen(progress_fhandle,file=fname, status='new', iostat=ok)
        if(progress > 1.0) then
            write(progress_fhandle, '(RD,F4.2)') 1.0
        else if(progress < 0.0) then
            write(progress_fhandle, '(RD,F4.2)') 0.0
        else
            write(progress_fhandle, '(RD,F4.2)') progress
        endif
        call fclose(progress_fhandle)
        call fname%kill
    end subroutine progressfile_update_part
    
    subroutine progressfile_complete_parts(nparts)
        integer, intent(in) :: nparts
        integer             :: ipart
        type(string)        :: fname
        fname = '.progress_' // trim(adjustl(int2str(ipart)))
        do ipart = 1, nparts
            if( file_exists(fname) ) call del_file(fname)
        enddo
        call progressfile_init()
        call progressfile_update(1.0)
        call fname%kill
    end subroutine progressfile_complete_parts
    
    function progress_estimate_preprocess_stream(n_imported, n_added) result(progress)
        integer, intent(in) :: n_imported, n_added
        real progress
        progress = real(n_imported) / real(n_added)
    end function progress_estimate_preprocess_stream

    function progress_estimate_2D(iter, overlap, overlap_lim, fracsrch, fracsrch_lim, distinpl, distinpl_lim ) result(progress)
        real, intent(in) :: iter, overlap, overlap_lim, fracsrch, fracsrch_lim, distinpl, distinpl_lim
        real             :: overlap_sc, fracsrch_sc, overlap_contrib, fracsrch_contrib, progress
        overlap_contrib = 0.25
        fracsrch_contrib = 0.75
        progress = 0.0
        overlap_sc = 0.0
        fracsrch_sc = 0.0
        ! Test overlap
        if(overlap > 0.0 .AND. overlap_lim > 0.0) then
            overlap_sc = overlap / overlap_lim
            if(overlap_sc > 1.0) then
                overlap_sc = overlap_contrib
            else
                overlap_sc = overlap_contrib * overlap_sc
            end if
        endif
        ! Test frac_srch
        if(fracsrch > 0.0 .AND. fracsrch_lim > 0.0) then
            fracsrch_sc = fracsrch / fracsrch_lim
            if(fracsrch_sc > 1.0) then
                fracsrch_sc = fracsrch_contrib
            else
                fracsrch_sc = fracsrch_contrib * fracsrch_sc
            end if
        endif
        progress = overlap_sc + fracsrch_sc
        ! If set, distinpl drives progress
        if(distinpl > 0.0 .AND. distinpl_lim > 0.0) then
            progress = distinpl_lim / distinpl
        endif
        ! Apply scalings depending on iteration
        if(iter == 1.0) then
            progress = 0.05
        else if(iter == 2.0) then
            progress = 0.1
        else if (iter <= 10) then
            progress = progress / 2
        end if
    end function progress_estimate_2D
    
    subroutine lastfoundfile_update()
        integer       :: lastfound_fhandle, ok
        type(string)  :: fname
        character(8)  :: date
        character(10) :: time
        character(5)  :: zone
        integer,dimension(8) :: values
        ! using keyword arguments
        call date_and_time(date,time,zone,values)
        call date_and_time(DATE=date,ZONE=zone)
        call date_and_time(TIME=time)
        call date_and_time(VALUES=values)
        fname = '.lastfound'
        if( file_exists(fname) ) call del_file(fname)
        call fopen(lastfound_fhandle,file=fname, status='new', iostat=ok)
        write(lastfound_fhandle, '(I4,A,I2.2,A,I2.2,A,I2.2,A,I2.2)') values(1), '/', values(2), '/', values(3), '_', values(5), ':', values(6)
        call fclose(lastfound_fhandle)
    end subroutine lastfoundfile_update

end module simple_progress
