!*****************************************************************************************
!> author: Jacob Williams
!  date: 12/22/2015
!  license: BSD
!
!  Numeric kind definitions.
!
!@note The default real kind (`wp`) can be
!      changed using optional preprocessor flags.
!      This library was built with real kind:
#ifdef REAL32
!      `real(kind=real32)` [4 bytes]
#elif REAL64
!      `real(kind=real64)` [8 bytes]
#elif REAL128
!      `real(kind=real128)` [16 bytes]
#else
!      `real(kind=real64)` [8 bytes]
#endif

    module dvode_kinds_module

    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

    private

#ifdef REAL32
    integer,parameter,public :: dvode_wp = real32   !! real kind used by this module [4 bytes]
#elif REAL64
    integer,parameter,public :: dvode_wp = real64   !! real kind used by this module [8 bytes]
#elif REAL128
    integer,parameter,public :: dvode_wp = real128  !! real kind used by this module [16 bytes]
#else
    integer,parameter,public :: dvode_wp = real64   !! real kind used by this module [8 bytes]
#endif

    end module dvode_kinds_module
!*****************************************************************************************
