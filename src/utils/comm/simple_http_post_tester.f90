!@descr: unit tests for simple_http_post against a loopback HTTP server (no network beyond localhost)
! A listener thread of ipc_tcp_socket_server answers HTTP/1.1 on 127.0.0.1: a request with
! a body is a POST and is echoed back with 201, a request without a body carries no body
! (libcurl issues a GET when no POST fields are set) and gets a canned JSON document with
! 200. Pinned: lifecycle, the body reaching the server verbatim, status code, content and
! content type as delivered, reset of the response between requests, and a fast failure
! with no listener. The request tests need the listener thread, which the ipc suites skip
! on macOS/FreeBSD and Windows.
module simple_http_post_tester
use iso_c_binding
use unix,                         only: c_pthread_mutex_init, c_pthread_mutex_destroy, c_pthread_mutex_lock, &
                                        c_pthread_mutex_unlock, c_read, c_close
use simple_string,                only: string
use simple_string_utils,          only: lowercase
use simple_http_post,             only: http_post, http_response
use simple_test_utils,            only: assert_true, assert_false, assert_int, assert_char, assert_string_eq
use simple_ipc_tcp_socket_server, only: ipc_tcp_socket_server, listener_args, recv_msg, repl_msg
implicit none

public  :: run_all_http_post_tests
private
#include "simple_local_flags.inc"

integer,          parameter :: REQ_BUF  = 4096
integer,          parameter :: RESP_BUF = 1024
character(len=*), parameter :: CRLF     = achar(13)//achar(10)
character(len=*), parameter :: CANNED   = '{"id":1,"title":"loopback","userId":1}'
character(len=*), parameter :: CTYPE    = 'application/json; charset=utf-8'
character(len=*), parameter :: BODY1    = '{"title":"foo","body":"bar","userId":1}'
character(len=*), parameter :: BODY2    = '{"title":"second"}'
integer,          parameter :: UNUSED_PORT = 65000

contains

    subroutine run_all_http_post_tests()
        write(*,'(A)') '**** running all http post tests ****'
        call test_create_and_kill()
#if defined(_WIN32) || defined(__FreeBSD__)
        write(*,'(A)') '**** skipping the loopback request tests on this platform (listener thread) ****'
#else
        call test_request_no_body()
        call test_request_with_body()
        call test_response_reset()
        call test_request_without_listener()
#endif
        write(*,'(A)') '**** http post tests done ****'
    end subroutine run_all_http_post_tests

    ! ---- the loopback server ---------------------------------------------------

    subroutine init_listener_args( largs )
        type(listener_args), intent(inout) :: largs
        integer(c_int) :: rc
        largs%fd       = -1_c_int
        largs%ready    = 0_c_int
        largs%data_ptr = c_null_ptr
        rc = c_pthread_mutex_init(largs%mutex, c_null_ptr)
        call assert_int(0, int(rc), 'listener mutex init succeeds')
    end subroutine init_listener_args

    subroutine destroy_listener_args( largs )
        type(listener_args), intent(inout) :: largs
        integer(c_int) :: rc
        rc = c_pthread_mutex_destroy(largs%mutex)
        call assert_int(0, int(rc), 'listener mutex destroy succeeds')
    end subroutine destroy_listener_args

    function loopback_url( port ) result( url )
        integer, intent(in) :: port
        type(string) :: url
        character(len=64) :: buf
        write(buf,'(A,I0,A)') 'http://127.0.0.1:', port, '/posts'
        url = string(trim(buf))
    end function loopback_url

    !> Content-Length of an HTTP header block, 0 when absent
    integer function content_length( hdr )
        character(len=*), intent(in) :: hdr
        character(len=len(hdr)) :: lc
        integer :: i, j, ios
        content_length = 0
        lc = lowercase(hdr)
        i  = index(lc, 'content-length:')
        if( i == 0 ) return
        i = i + len('content-length:')
        j = index(lc(i:), CRLF)
        if( j == 0 ) return
        read(lc(i:i+j-2),*,iostat=ios) content_length
        if( ios /= 0 ) content_length = 0
    end function content_length

    !> listener thread: serve requests until the server's kill sentinel arrives
    subroutine http_listener( arg_ptr ) bind(c)
        type(c_ptr), value :: arg_ptr
        type(listener_args), pointer :: args
        character(kind=c_char, len=REQ_BUF),  target :: req, chunk
        character(kind=c_char, len=RESP_BUF), target :: resp
        character(len=16)      :: clen_str
        integer(c_int)         :: conn_fd, rc
        integer(c_size_t)      :: nr
        integer                :: nread, ntot, hdr_end, clen, body_from, nbody, nresp, nwrite
        logical                :: ok
        if( .not. c_associated(arg_ptr) ) return
        call c_f_pointer(arg_ptr, args)
        rc = c_pthread_mutex_lock(args%mutex)
        args%ready = 1_c_int
        rc = c_pthread_mutex_unlock(args%mutex)
        do
            req     = ''
            conn_fd = -1_c_int
            call recv_msg(args%fd, conn_fd, req, nread, ok, close_after_read=.false.)
            if( .not. ok ) exit
            ! anything that is not an HTTP request line is the kill sentinel of server%kill
            if( nread < 4 .or. (req(1:4) /= 'POST' .and. req(1:4) /= 'GET ') )then
                rc = c_close(conn_fd)
                exit
            endif
            ! read until the header block is complete and Content-Length bytes of body are in
            ntot = nread
            do
                hdr_end = index(req(1:ntot), CRLF//CRLF)
                if( hdr_end > 0 )then
                    clen = content_length(req(1:hdr_end+3))
                    if( ntot - (hdr_end + 3) >= clen ) exit
                endif
                if( ntot >= REQ_BUF ) exit
                chunk = ''
                nr = c_read(conn_fd, c_loc(chunk), int(REQ_BUF - ntot, c_size_t))
                if( nr == 0_c_size_t .or. nr > int(REQ_BUF, c_size_t) ) exit
                req(ntot+1:ntot+int(nr)) = chunk(1:int(nr))
                ntot = ntot + int(nr)
            enddo
            if( hdr_end > 0 )then
                body_from = hdr_end + 4
                nbody     = max(0, min(clen, ntot - body_from + 1))
            else
                body_from = ntot + 1
                nbody     = 0
            endif
            if( req(1:4) == 'POST' )then
                write(clen_str,'(I0)') nbody
                resp  = 'HTTP/1.1 201 Created'//CRLF//'Content-Type: '//CTYPE//CRLF// &
                       &'Content-Length: '//trim(clen_str)//CRLF//'Connection: close'//CRLF//CRLF// &
                       &req(body_from:body_from+nbody-1)
            else
                write(clen_str,'(I0)') len(CANNED)
                resp  = 'HTTP/1.1 200 OK'//CRLF//'Content-Type: '//CTYPE//CRLF// &
                       &'Content-Length: '//trim(clen_str)//CRLF//'Connection: close'//CRLF//CRLF//CANNED
            endif
            nresp = len_trim(resp)
            call repl_msg(conn_fd, resp(1:nresp), nwrite, ok)
        enddo
    end subroutine http_listener

    ! ---- tests -----------------------------------------------------------------

    subroutine test_create_and_kill()
        type(http_post) :: post
        write(*,'(A)') 'test_create_and_kill'
        call post%new(loopback_url(UNUSED_PORT))
        call assert_true(post%initialised(),       'post is initialised after new()')
        call post%kill()
        call assert_true(.not. post%initialised(), 'post is uninitialised after kill()')
    end subroutine test_create_and_kill

    !> without a body the request carries none and the server answers 200 with its document
    subroutine test_request_no_body()
        type(ipc_tcp_socket_server)  :: server
        type(listener_args), target  :: largs
        type(http_post)              :: post
        type(http_response)          :: response
        write(*,'(A)') 'test_request_no_body'
        call init_listener_args(largs)
        call server%new(c_funloc(http_listener), c_loc(largs))
        call post%new(loopback_url(server%get_port()))
        call assert_true(post%request(response),     'body-less request succeeds')
        call assert_int(200, response%code,           'body-less request returns HTTP 200')
        call assert_string_eq(CANNED, response%content, 'content is the server document, verbatim')
        call assert_string_eq(CTYPE,  response%content_type, 'content type is delivered as sent')
        call post%kill()
        call assert_true(.not. post%initialised(),   'post is uninitialised after kill()')
        call server%kill()
        call destroy_listener_args(largs)
    end subroutine test_request_no_body

    !> the body reaches the server verbatim and comes back with 201
    subroutine test_request_with_body()
        type(ipc_tcp_socket_server)  :: server
        type(listener_args), target  :: largs
        type(http_post)              :: post
        type(http_response)          :: response
        write(*,'(A)') 'test_request_with_body'
        call init_listener_args(largs)
        call server%new(c_funloc(http_listener), c_loc(largs))
        call post%new(loopback_url(server%get_port()))
        call assert_true(post%request(response, request_str=string(BODY1)), 'request with body succeeds')
        call assert_int(201, response%code,                  'request with body returns HTTP 201')
        call assert_string_eq(BODY1, response%content,       'the POST body arrived verbatim (echoed)')
        call assert_string_eq(CTYPE, response%content_type,  'content type is delivered as sent')
        call post%kill()
        call assert_true(.not. post%initialised(), 'post is uninitialised after kill()')
        call server%kill()
        call destroy_listener_args(largs)
    end subroutine test_request_with_body

    !> the response is reset between requests on the same object
    subroutine test_response_reset()
        type(ipc_tcp_socket_server)  :: server
        type(listener_args), target  :: largs
        type(http_post)              :: post
        type(http_response)          :: response
        write(*,'(A)') 'test_response_reset'
        call init_listener_args(largs)
        call server%new(c_funloc(http_listener), c_loc(largs))
        call post%new(loopback_url(server%get_port()))
        call assert_true(post%request(response, request_str=string(BODY1)), 'first request succeeds')
        call assert_string_eq(BODY1, response%content, 'first body echoed')
        call assert_true(post%request(response, request_str=string(BODY2)), 'second request succeeds')
        call assert_string_eq(BODY2, response%content, 'second body echoed alone: the response was reset')
        call assert_int(201, response%code, 'second request returns HTTP 201')
        call assert_true(post%request(response), 'third, body-less request succeeds')
        call assert_int(200, response%code, 'body-less request after POSTs returns HTTP 200')
        call assert_string_eq(CANNED, response%content, 'body-less response replaces the echoed body')
        call post%kill()
        call assert_true(.not. post%initialised(), 'post is uninitialised after kill()')
        call server%kill()
        call destroy_listener_args(largs)
    end subroutine test_response_reset

    !> no listener: the request fails at once (connection refused), no code, no content
    subroutine test_request_without_listener()
        type(http_post)     :: post
        type(http_response) :: response
        write(*,'(A)') 'test_request_without_listener'
        call post%new(loopback_url(UNUSED_PORT))
        call assert_false(post%request(response, request_str=string(BODY1)), 'request without a listener fails')
        call assert_int(0, response%code, 'no status code without a listener')
        call assert_false(response%content%is_allocated(), 'no content without a listener')
        call post%kill()
    end subroutine test_request_without_listener

end module simple_http_post_tester
