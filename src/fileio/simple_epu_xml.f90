!@descr: the beam shift of an EPU metadata file (BeamShift's a:_x and a:_y), read without an XML library
! EPU writes one XML document per exposure. SIMPLE needs only the beam shift from it, which sits anywhere
! in the document as
!   <BeamShift xmlns:a="..."><a:_x>value</a:_x><a:_y>value</a:_y></BeamShift>
! The first BeamShift element is read, and in it the first a:_x and a:_y, as the DOM lookup this replaces did.
module simple_epu_xml
implicit none
private

public :: read_epu_beamshift

contains

    !> shiftx, shifty: the a:_x and a:_y of the first BeamShift element of the EPU metadata file fname;
    !! l_ok is false, and both shifts zero, when the file cannot be read or holds no such element
    subroutine read_epu_beamshift( fname, shiftx, shifty, l_ok )
        character(len=*), intent(in)  :: fname
        real,             intent(out) :: shiftx, shifty
        logical,          intent(out) :: l_ok
        character(len=:), allocatable :: doc, beamshift
        logical :: l_okx, l_oky
        shiftx = 0.
        shifty = 0.
        l_ok   = .false.
        if( .not. read_text(fname, doc) ) return
        if( .not. element_content(doc, 'BeamShift', beamshift) ) return
        call element_real(beamshift, 'a:_x', shiftx, l_okx)
        call element_real(beamshift, 'a:_y', shifty, l_oky)
        l_ok = l_okx .and. l_oky
        if( .not. l_ok )then
            shiftx = 0.
            shifty = 0.
        endif
    end subroutine read_epu_beamshift

    !> the whole file as one string
    logical function read_text( fname, txt )
        character(len=*),              intent(in)  :: fname
        character(len=:), allocatable, intent(out) :: txt
        integer :: funit, fsize, ios
        read_text = .false.
        inquire(file=fname, size=fsize)
        if( fsize <= 0 ) return
        open(newunit=funit, file=fname, access='stream', form='unformatted', status='old', action='read', &
            &iostat=ios)
        if( ios /= 0 ) return
        allocate(character(len=fsize) :: txt)
        read(funit, iostat=ios) txt
        close(funit)
        read_text = ios == 0
    end function read_text

    !> the text between the first <name> (or <name attributes>) start tag of txt and the next </name>
    logical function element_content( txt, name, content )
        character(len=*),              intent(in)  :: txt, name
        character(len=:), allocatable, intent(out) :: content
        character(len=*), parameter :: AFTER_NAME = '>/ '//achar(9)//achar(10)//achar(13)
        integer :: from, start, past_name, gt, endtag
        element_content = .false.
        from = 1
        do
            start = index(txt(from:), '<'//name)
            if( start == 0 ) return
            start     = from + start - 1
            past_name = start + len(name) + 1
            if( past_name > len(txt) ) return
            ! a longer name that begins with this one, such as <BeamShiftX>, is not a match
            if( scan(txt(past_name:past_name), AFTER_NAME) > 0 ) exit
            from = past_name
        end do
        gt = index(txt(past_name:), '>')
        if( gt == 0 ) return
        gt = past_name + gt - 1
        if( txt(gt-1:gt-1) == '/' ) return              ! <name/> holds nothing
        endtag = index(txt(gt+1:), '</'//name)
        if( endtag == 0 ) return
        content = txt(gt+1:gt+endtag-1)
        element_content = .true.
    end function element_content

    !> the number held by the first <name> element of txt
    subroutine element_real( txt, name, val, l_ok )
        character(len=*), intent(in)  :: txt, name
        real,             intent(out) :: val
        logical,          intent(out) :: l_ok
        character(len=:), allocatable :: content
        integer :: i, ios
        val  = 0.
        l_ok = element_content(txt, name, content)
        if( .not. l_ok ) return
        ! line breaks and tabs around the number read as blanks
        do i = 1, len(content)
            if( iachar(content(i:i)) < 32 ) content(i:i) = ' '
        end do
        l_ok = len_trim(content) > 0
        if( .not. l_ok ) return
        read(content, *, iostat=ios) val
        l_ok = ios == 0
        if( .not. l_ok ) val = 0.
    end subroutine element_real

end module simple_epu_xml
