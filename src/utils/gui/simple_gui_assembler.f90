!@descr: assembles gui metadata objects into json to be sent to nice
module simple_gui_assembler
use simple_string
use simple_gui_metadata_api
implicit none

type :: gui_assembler
  private
  type(json_core)           :: json
  type(json_value), pointer :: json_root
  logical                   :: init = .false.

contains

  procedure :: new
  procedure :: kill
  procedure :: to_string
  procedure :: is_associated
  procedure :: assemble_stream_preprocess

end type gui_assembler

contains

  subroutine new( self )
    class(gui_assembler),                 intent(inout) :: self
    if( self%init ) call self%kill()
    self%init = .true.
    call self%json%initialize(no_whitespace=.true., compact_reals=.true.)
    call self%json%create_object(self%json_root, '') 
  end subroutine new

  subroutine kill( self )
    class(gui_assembler), intent(inout) :: self
    call self%json%destroy(self%json_root)
    nullify(self%json_root)
  end subroutine kill

  subroutine assemble_stream_preprocess( self, meta_preprocess, meta_micrographs, meta_histograms )
    class(gui_assembler),                              intent(inout) :: self
    type(gui_metadata_stream_preprocess),              intent(inout) :: meta_preprocess
    type(gui_metadata_micrograph),        allocatable, intent(inout) :: meta_micrographs(:)
    type(gui_metadata_histogram),         allocatable, intent(inout) :: meta_histograms(:)
    !type(gui_metadata_histogram),         allocatable, intent(inout) :: meta_timeplots(:)
    type(json_value),                     pointer                    :: json_ptr, json_mics_ptr, json_hists_ptr
    integer                                                          :: i_mic, n_mics, i_hist, n_hists
    json_ptr => meta_preprocess%jsonise()
    if( associated(json_ptr) ) then
      call self%json%rename(json_ptr, 'preprocessing')
      ! add micrographs array if present
      if( allocated(meta_micrographs) ) then
        n_mics = size(meta_micrographs)
        if( n_mics > 0 ) then
          call self%json%create_array(json_mics_ptr, 'micrographs')
          do i_mic=1, n_mics
            call self%json%add(json_mics_ptr, meta_micrographs(i_mic)%jsonise())
          enddo
          call self%json%add(json_ptr, json_mics_ptr)
        endif
      endif
      ! add histograms array if present
      if( allocated(meta_histograms) ) then
        n_hists = size(meta_histograms)
        if( n_hists > 0 ) then
          call self%json%create_array(json_hists_ptr, 'histograms')
          call self%json%add(json_ptr, json_hists_ptr)
          do i_hist=1, n_hists
            call self%json%add(json_hists_ptr, meta_histograms(i_hist)%jsonise())
          enddo
        endif
      endif
      call self%json%add(self%json_root, json_ptr)
    endif
    ! cleanup
    nullify(json_ptr, json_mics_ptr, json_hists_ptr)
  end subroutine assemble_stream_preprocess

  function to_string( self ) result( str )
    class(gui_assembler),     intent(inout) :: self
    type(string)                            :: str
    character(kind=CK,len=:), allocatable   :: buffer
    call self%json%print_to_string(self%json_root, buffer)
    str = buffer
    if( allocated(buffer) ) deallocate(buffer)
  end function to_string

  function is_associated( self ) result( assoc )
    class(gui_assembler), intent(inout) :: self
    logical                             :: assoc
    assoc = associated(self%json_root)
  end function is_associated

end module simple_gui_assembler