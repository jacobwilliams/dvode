!*****************************************************************************************
!> author: Jacob Williams
!  date: 12/22/2015
!  license: BSD
!
!  Numeric kind definitions.

    module dvode_kinds

    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

    private

    integer,parameter,public :: wp = real64  !! Using "double precision" real kinds

    public :: get_wp

    contains
!*****************************************************************************************

!*****************************************************************************************
    integer function get_wp()

    !! Just a function that returns the value of `wp`, the
    !! working precision kind value for real numbers.

    implicit none

    get_wp = wp

    end function get_wp
!*****************************************************************************************

    end module dvode_kinds
!*****************************************************************************************
