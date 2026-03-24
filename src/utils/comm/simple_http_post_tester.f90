!@descr: unit tests for simple_http_post (lifecycle, body-less POST, POST with body, response reset)
!==============================================================================
! MODULE: simple_http_post_tester
!
! PURPOSE:
!   Exercises the http_post and http_response types across four test cases:
!     1. test_create_and_kill      — object lifecycle (new/kill/initialised)
!     2. test_request_no_body      — body-less POST; checks HTTP 200, and
!                                    FNV-1a hashes of the response body and
!                                    Content-Type header
!     3. test_request_with_body    — POST with a JSON body; checks HTTP 201
!                                    and that both response fields are populated
!     4. test_response_reset       — calls request() twice on the same object;
!                                    verifies the response is fully cleared
!                                    between calls (identical hashes expected)
!
! ENTRY POINT:
!   run_all_http_post_tests — call this from the top-level test harness.
!
! EXTERNAL DEPENDENCY:
!   https://jsonplaceholder.typicode.com — requires network access at test time.
!
! DEPENDENCIES:
!   simple_string, simple_http_post, simple_test_utils
!==============================================================================
module simple_http_post_tester
  use simple_string,     only: string
  use simple_http_post,  only: http_post, http_response
  use simple_test_utils, only: assert_true, assert_int, assert_char

  implicit none

  public  :: run_all_http_post_tests
  private
#include "simple_local_flags.inc"

contains

  ! Run all http_post unit tests in order.
  subroutine run_all_http_post_tests()
    write(*,'(A)') '**** running all http post tests ****'
    call test_create_and_kill()
    call test_request_no_body()
    call test_request_with_body()
    call test_response_reset()
  end subroutine run_all_http_post_tests

  ! Verify that new() marks the object as initialised and kill() clears it.
  subroutine test_create_and_kill()
    type(http_post) :: post
    write(*,'(A)') 'test_create_and_kill'
    call post%new(string('https://jsonplaceholder.typicode.com/posts'))
    call assert_true(post%initialised(),       'post is initialised after new()')
    call post%kill()
    call assert_true(.not. post%initialised(), 'post is uninitialised after kill()')
  end subroutine test_create_and_kill

  ! Issue a body-less POST to JSONPlaceholder and check:
  !   - request() returns .true.
  !   - HTTP status code is 200
  !   - response body and Content-Type match expected FNV-1a hashes
  subroutine test_request_no_body()
    type(http_post)     :: post
    type(http_response) :: response
    type(string)        :: content_hash, type_hash
    write(*,'(A)') 'test_request_no_body'
    call post%new(string('https://jsonplaceholder.typicode.com/posts'))
    call assert_true(post%initialised(),                 'post is initialised')
    call assert_true(post%request(response),             'body-less request succeeds')
    call assert_int( response%code, 200,                 'body-less request returns HTTP 200')
    content_hash = response%content%to_fnv1a_hash64()
    type_hash    = response%content_type%to_fnv1a_hash64()
    call assert_char(content_hash%to_char(), 'FF0C62AB7B8AEFF4', 'response content hash matches'     )
    call assert_char(type_hash%to_char(),    'AE256DC329DC5E8C', 'response content_type hash matches')
    call post%kill()
    call assert_true(.not. post%initialised(),           'post is uninitialised after kill()')
  end subroutine test_request_no_body

  ! POST a minimal JSON body to JSONPlaceholder and check:
  !   - request() returns .true.
  !   - HTTP status code is 201 (resource created)
  !   - response body and Content-Type are non-empty
  subroutine test_request_with_body()
    type(http_post)     :: post
    type(http_response) :: response
    write(*,'(A)') 'test_request_with_body'
    call post%new(string('https://jsonplaceholder.typicode.com/posts'))
    call assert_true(post%initialised(), 'post is initialised')
    call assert_true( &
      post%request(response, request_str=string('{"title":"foo","body":"bar","userId":1}')), &
      'request with body succeeds')
    call assert_int( response%code, 201,                     'request with body returns HTTP 201')
    call assert_true(response%content%is_allocated(),        'response content is populated')
    call assert_true(response%content_type%is_allocated(),   'response content_type is populated')
    call post%kill()
    call assert_true(.not. post%initialised(), 'post is uninitialised after kill()')
  end subroutine test_request_with_body

  ! Call request() twice on the same object with the same endpoint. Verify
  ! that the response content is fully reset between calls: the FNV-1a hash
  ! of the body must be identical for both, proving the first call's data was
  ! not accumulated into the second.
  subroutine test_response_reset()
    type(http_post)     :: post
    type(http_response) :: response
    type(string)        :: hash1, hash2
    write(*,'(A)') 'test_response_reset'
    call post%new(string('https://jsonplaceholder.typicode.com/posts'))
    call assert_true(post%request(response), 'first request succeeds')
    hash1 = response%content%to_fnv1a_hash64()
    call assert_true(post%request(response), 'second request succeeds')
    hash2 = response%content%to_fnv1a_hash64()
    call assert_char(hash1%to_char(), hash2%to_char(), &
      'response content hash matches on repeat request (response was reset)')
    call post%kill()
    call assert_true(.not. post%initialised(), 'post is uninitialised after kill()')
  end subroutine test_response_reset

end module simple_http_post_tester
