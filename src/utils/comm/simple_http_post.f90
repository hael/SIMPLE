!@descr: libcurl-based HTTP POST client with response capture
!==============================================================================
! MODULE: simple_http_post
!
! PURPOSE:
!   Provides a thin Fortran wrapper around libcurl for issuing HTTP POST
!   requests and capturing the response body, HTTP status code, and
!   content-type header.
!
! TYPES:
!   http_response — holds the server's response: body, content-type, and
!                   HTTP status code.
!   http_post     — manages a persistent curl session for a fixed URL;
!                   call new() once, request() as many times as needed,
!                   then kill().
!
! USAGE:
!   type(http_post)     :: poster
!   type(http_response) :: resp
!   logical             :: ok
!   call poster%new(string('http://example.com/api'))
!   ok = poster%request(resp, request_str=string('{"key":"val"}'))
!   call poster%kill()
!
! DEPENDENCIES:
!   iso_c_binding (intrinsic), curl
!==============================================================================
module simple_http_post
  use, intrinsic :: iso_c_binding, only: c_ptr, c_loc, c_funloc, &
                                         c_associated, c_f_pointer, c_size_t
  use curl, only: c_f_str_ptr, curl_global_init, curl_global_cleanup,         &
                  curl_easy_init, curl_easy_setopt, curl_easy_perform,        &
                  curl_easy_getinfo, curl_easy_cleanup,                       &
                  CURL_GLOBAL_DEFAULT, CURLE_OK,                              &
                  CURLOPT_URL, CURLOPT_VERBOSE, CURLOPT_TIMEOUT,              &
                  CURLOPT_NOSIGNAL, CURLOPT_COPYPOSTFIELDS,                   &
                  CURLOPT_POSTFIELDSIZE,                                      &
                  CURLOPT_WRITEFUNCTION, CURLOPT_WRITEDATA,                   &
                  CURLINFO_RESPONSE_CODE, CURLINFO_CONTENT_TYPE,              &
                  CURLOPT_SSL_VERIFYPEER
  use simple_string, only: string
  use simple_error,  only: simple_exception

  implicit none

  public  :: http_post
  public  :: http_response
  private
#include "simple_local_flags.inc"

  logical, parameter :: VERIFY_PEER = .false.  ! set to .false. to disable SSL peer verification (not recommended)

  ! Holds the server's reply after a successful request.
  type :: http_response
    type(string) :: content       ! response body (accumulated across callbacks)
    type(string) :: content_type  ! Content-Type header value
    integer      :: code          ! HTTP status code (e.g. 200, 404)
  end type http_response

  ! Manages a reusable curl session bound to a single URL.
  type :: http_post
    private
    type(string) :: url
    type(c_ptr)  :: curl_ptr
    logical      :: l_failed      = .false.
    logical      :: l_initialized = .false.
  contains
    procedure :: new
    procedure :: kill
    procedure :: initialised
    procedure :: request
  end type http_post

contains

  ! Initialise global curl state and bind this object to url.
  ! Must be called before any request().
  subroutine new( self, url )
    class(http_post), intent(inout) :: self
    type(string),     intent(in)    :: url
    integer                         :: rc
    self%l_initialized = .true.
    self%url           = url
    rc = curl_global_init(CURL_GLOBAL_DEFAULT)
    if( rc /= CURLE_OK ) THROW_HARD('Error: failed to set curl_global_init')
  end subroutine new

  ! Release global curl resources and reset initialisation flag.
  subroutine kill( self )
    class(http_post), intent(inout) :: self
    call curl_global_cleanup()
    call self%url%kill()
    self%l_initialized = .false.
  end subroutine kill

  ! Returns .true. if new() has been called and kill() has not.
  function initialised( self ) result( l_initialized )
    class(http_post), intent(in) :: self
    logical                      :: l_initialized
    l_initialized = self%l_initialized
  end function initialised

  ! Issue an HTTP POST to self%url. If request_str is present it is sent as
  ! the raw POST body; otherwise a body-less POST is issued. On return,
  ! response holds the HTTP status code, content-type, and accumulated body.
  ! Returns .true. on success, .false. if any curl operation failed.
  function request( self, response, request_str ) result( l_success )
    class(http_post),             intent(inout) :: self
    type(http_response), target,  intent(inout) :: response
    type(string),        optional, intent(in)   :: request_str
    character(len=:),              allocatable  :: content_type
    character(len=:),              allocatable  :: request_body
    integer                                     :: rc
    logical                                     :: l_success
    l_success     = .true.
    self%l_failed = .false.
    response%code = 0
    call response%content%kill()
    call response%content_type%kill()
    ! Initialise a new easy-curl handle for this request
    self%curl_ptr = curl_easy_init()
    if( .not. c_associated(self%curl_ptr) ) THROW_HARD('Error: curl_easy_init() failed')
    ! Configure curl options
    rc = curl_easy_setopt(self%curl_ptr, CURLOPT_URL,           self%url%to_char())
    if( rc /= CURLE_OK ) THROW_HARD('Error: failed to set curl option CURLOPT_URL')
    rc = curl_easy_setopt(self%curl_ptr, CURLOPT_VERBOSE,       0)
    if( rc /= CURLE_OK ) THROW_HARD('Error: failed to set curl option CURLOPT_VERBOSE')
    rc = curl_easy_setopt(self%curl_ptr, CURLOPT_TIMEOUT,       20)
    if( rc /= CURLE_OK ) THROW_HARD('Error: failed to set curl option CURLOPT_TIMEOUT')
    rc = curl_easy_setopt(self%curl_ptr, CURLOPT_NOSIGNAL,      1)
    if( rc /= CURLE_OK ) THROW_HARD('Error: failed to set curl option CURLOPT_NOSIGNAL')
    if( .not. VERIFY_PEER ) then
      THROW_WARN('Warning: SSL peer verification is disabled')
      rc = curl_easy_setopt(self%curl_ptr, CURLOPT_SSL_VERIFYPEER, 0)
      if( rc /= CURLE_OK ) THROW_HARD('Error: failed to set curl option CURLOPT_SSL_VERIFYPEER')
    endif
    ! Attach POST body if provided
    if( present(request_str) ) then
      request_body = request_str%to_char()
      rc = curl_easy_setopt(self%curl_ptr, CURLOPT_POSTFIELDSIZE, len(request_body))
      if( rc /= CURLE_OK ) THROW_HARD('Error: failed to set curl option CURLOPT_POSTFIELDSIZE')
      rc = curl_easy_setopt(self%curl_ptr, CURLOPT_COPYPOSTFIELDS, request_body)
      if( rc /= CURLE_OK ) THROW_HARD('Error: failed to set curl option CURLOPT_COPYPOSTFIELDS')
    endif
    ! Register response-body callback
    rc = curl_easy_setopt(self%curl_ptr, CURLOPT_WRITEFUNCTION,  c_funloc(response_callback))
    if( rc /= CURLE_OK ) THROW_HARD('Error: failed to set curl option CURLOPT_WRITEFUNCTION')
    rc = curl_easy_setopt(self%curl_ptr, CURLOPT_WRITEDATA,      c_loc(response))
    if( rc /= CURLE_OK ) THROW_HARD('Error: failed to set curl option CURLOPT_WRITEDATA')
    ! Execute the request
    rc = curl_easy_perform(self%curl_ptr)
    if( rc == CURLE_OK ) then
      rc = curl_easy_getinfo(self%curl_ptr, CURLINFO_RESPONSE_CODE, response%code)
      if( rc /= CURLE_OK ) then
        THROW_WARN('Warning: failed to get CURLINFO_RESPONSE_CODE')
        self%l_failed = .true.
      endif
      rc = curl_easy_getinfo(self%curl_ptr, CURLINFO_CONTENT_TYPE, content_type)
      if( rc == CURLE_OK ) then
        response%content_type = content_type
      else
        THROW_WARN('Warning: failed to get CURLINFO_CONTENT_TYPE')
        self%l_failed = .true.
      endif
    else
      THROW_WARN('Warning: curl request failed')
      self%l_failed = .true.
    endif
    ! Clean up curl handle and temporary buffers
    call curl_easy_cleanup(self%curl_ptr)
    if( allocated(content_type) ) deallocate(content_type)
    if( allocated(request_body) ) deallocate(request_body)
    if( self%l_failed ) l_success = .false.
  end function request

  ! libcurl write callback: invoked for each received chunk of the response
  ! body. Accumulates chunks into response_ptr%content.
  ! ptr        — C pointer to the current chunk
  ! size       — element size in bytes
  ! nmemb      — number of elements in this chunk
  ! client_data — user-supplied C pointer (cast to http_response here)
  ! Returns size * nmemb on success, 0 to signal an error to curl.
  function response_callback( ptr, size, nmemb, client_data ) bind(c)
    type(c_ptr),            intent(in), value :: ptr
    integer(kind=c_size_t), intent(in), value :: size
    integer(kind=c_size_t), intent(in), value :: nmemb
    type(c_ptr),            intent(in), value :: client_data
    integer(kind=c_size_t)                    :: response_callback
    type(http_response),    pointer           :: response_ptr
    character(len=:),       allocatable       :: buf
    response_callback = int(0, kind=c_size_t)
    ! Guard against null C pointers before dereferencing
    if( .not. c_associated(ptr) )         return
    if( .not. c_associated(client_data) ) return
    ! Convert C pointers to Fortran pointers and copy the chunk
    call c_f_pointer(client_data, response_ptr)
    call c_f_str_ptr(ptr, buf, int(size * nmemb, kind=8))
    if( .not. allocated(buf) ) return
    ! Append chunk to the accumulated response body
    if( response_ptr%content%is_allocated() ) then
      response_ptr%content = response_ptr%content // buf
    else
      response_ptr%content = buf
    endif
    deallocate(buf)
    response_callback = size * nmemb
  end function response_callback

end module simple_http_post
