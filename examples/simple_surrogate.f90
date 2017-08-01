module simple_surrogate
implicit none
public :: surrogate
private
! This stateless type serves only for purposes of extension by other types.
! Its role is to serve as a substitute for the child type when that type
! is inaccessible because of Fortran's prohibition against circular references
type, abstract :: surrogate
end type
end module