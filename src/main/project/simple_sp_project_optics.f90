submodule(simple_sp_project) simple_sp_project_optics
implicit none
#include "simple_local_flags.inc"
contains

    module subroutine import_optics_map( self, mapfileprefix )
        class(sp_project), intent(inout) :: self
        class(string),     intent(in)    :: mapfileprefix
        type(nrtxtfile)      :: mapfile
        real,    allocatable :: map_entries(:,:)
        integer, allocatable :: mics_optics_map(:)
        real                 :: min_importind, max_importind
        integer              :: il, nl, imic, iptcl, importind, stkind
        if(self%os_mic%get_noris() .eq. 0)                return
        if(.not.file_exists(mapfileprefix//METADATA_EXT)) then
            write(logfhandle, '(A)') mapfileprefix%to_char()//METADATA_EXT // " does not exist"
            return
        endif
        if(.not.file_exists(mapfileprefix//TXT_EXT)) then
            write(logfhandle, '(A)') mapfileprefix%to_char()//TXT_EXT // " does not exist"
            return
        endif
        call mapfile%new(mapfileprefix//TXT_EXT, 1)
        if( mapfile%get_nrecs_per_line() /= 2 )then
            THROW_WARN('INVALID FORMAT FOR: '//mapfileprefix%to_char()//TXT_EXT)
            call mapfile%kill
        endif
        nl = mapfile%get_ndatalines()
        allocate(map_entries(2, nl))
        do il = 1,nl
            call mapfile%readNextDataLine(map_entries(:,il))
        enddo
        call mapfile%kill()
        call self%os_mic%minmax('importind', min_importind, max_importind)
        allocate(mics_optics_map(int(max_importind)), source=1)
        do il = 1,nl
            if(int(map_entries(1, il)) <= max_importind) mics_optics_map(int(map_entries(1, il))) = int(map_entries(2,il))
        enddo
        call self%os_mic%set_all2single('ogid', 1)
        do imic=1, self%os_mic%get_noris()
            importind = self%os_mic%get_int(imic, 'importind')
            call self%os_mic%set(imic, 'ogid', mics_optics_map(importind))
        enddo
        if(self%os_stk%get_noris() .eq. self%os_mic%get_noris()) then
            do imic=1, self%os_mic%get_noris()
                call self%os_stk%set(imic, 'ogid', self%os_mic%get(imic, 'ogid'))
            enddo
        endif
        if(self%os_ptcl2D%get_noris() > 0) then
            do iptcl=1, self%os_ptcl2D%get_noris()
                stkind = self%os_ptcl2D%get_int(iptcl, 'stkind')
                call self%os_ptcl2D%set(iptcl, 'ogid', self%os_stk%get(stkind, 'ogid'))
            end do
        end if
        call self%read_segment('optics', mapfileprefix//METADATA_EXT)
        if(allocated(map_entries))     deallocate(map_entries)
        if(allocated(mics_optics_map)) deallocate(mics_optics_map)
    end subroutine import_optics_map

end submodule simple_sp_project_optics
