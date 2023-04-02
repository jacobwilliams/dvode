!*******************************************************************************
!> license: BSD
!
!  BLAS support routines for DVODE.
!  These have been refactored into modern Fortran.

    module dvode_blas_module
      
#ifndef HAS_BLAS

    use dvode_kinds_module, only: wp => dvode_wp

    implicit none

    private

    real(wp),parameter :: zero   = 0.0_wp
    real(wp),parameter :: one    = 1.0_wp
    real(wp),parameter :: two    = 2.0_wp
    real(wp),parameter :: four   = 4.0_wp
    real(wp),parameter :: ten    = 10.0_wp
    real(wp),parameter :: hun    = 100.0_wp

    public :: daxpy,dcopy,ddot,dnrm2,dscal,idamax

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  constant times a vector plus a vector.
!  uses unrolled loops for increments equal to one.
!
!### Author
!  jack dongarra, linpack, 3/11/78.

      subroutine daxpy(n,da,dx,incx,dy,incy)
      implicit none

      integer,intent(in) :: n
      real(wp),intent(in) :: da
      real(wp),intent(in)  :: dx(*) 
      integer,intent(in)  :: incx  
      real(wp),intent(inout) :: dy(*)
      integer,intent(in)  :: incy

      integer :: i , ix , iy , m , mp1

      if ( n<=0 ) return
      if ( da==zero ) return
      if ( incx==1 .and. incy==1 ) then

        ! code for both increments equal to 1

        ! clean-up loop

         m = mod(n,4)
         if ( m/=0 ) then
            do i = 1 , m
               dy(i) = dy(i) + da*dx(i)
            end do
            if ( n<4 ) return
         end if
         mp1 = m + 1
         do i = mp1 , n , 4
            dy(i) = dy(i) + da*dx(i)
            dy(i+1) = dy(i+1) + da*dx(i+1)
            dy(i+2) = dy(i+2) + da*dx(i+2)
            dy(i+3) = dy(i+3) + da*dx(i+3)
         end do

      else

         ! code for unequal increments or equal increments
         ! not equal to 1

         ix = 1
         iy = 1
         if ( incx<0 ) ix = (-n+1)*incx + 1
         if ( incy<0 ) iy = (-n+1)*incy + 1
         do i = 1 , n
            dy(iy) = dy(iy) + da*dx(ix)
            ix = ix + incx
            iy = iy + incy
         end do

      end if

      end subroutine daxpy
!*******************************************************************************

!*******************************************************************************
!>
!  copies a vector, x, to a vector, y.
!  uses unrolled loops for increments equal to one.
!
!### Author
!  jack dongarra, linpack, 3/11/78.

    subroutine dcopy(n,dx,incx,dy,incy)

      implicit none

      integer,intent(in) :: n 
      real(wp),intent(in) :: dx(*) 
      integer,intent(in) :: incx 
      real(wp),intent(inout) :: dy(*)
      integer,intent(in) :: incy 

      integer :: i , ix , iy , m , mp1 

      if ( n<=0 ) return
      if ( incx==1 .and. incy==1 ) then

         ! code for both increments equal to 1

         ! clean-up loop

         m = mod(n,7)
         if ( m/=0 ) then
            do i = 1 , m
               dy(i) = dx(i)
            end do
            if ( n<7 ) return
         end if
         mp1 = m + 1
         do i = mp1 , n , 7
            dy(i) = dx(i)
            dy(i+1) = dx(i+1)
            dy(i+2) = dx(i+2)
            dy(i+3) = dx(i+3)
            dy(i+4) = dx(i+4)
            dy(i+5) = dx(i+5)
            dy(i+6) = dx(i+6)
         end do

      else

         ! code for unequal increments or equal increments
         ! not equal to 1

         ix = 1
         iy = 1
         if ( incx<0 ) ix = (-n+1)*incx + 1
         if ( incy<0 ) iy = (-n+1)*incy + 1
         do i = 1 , n
            dy(iy) = dx(ix)
            ix = ix + incx
            iy = iy + incy
         end do

      end if

    end subroutine dcopy
!*******************************************************************************

!*******************************************************************************
!>
!  forms the dot product of two vectors.
!  uses unrolled loops for increments equal to one.
!
!### Author
!  jack dongarra, linpack, 3/11/78.

      real(wp) function ddot(n,dx,incx,dy,incy)

      implicit none

      integer,intent(in) :: n
      real(wp),intent(in) :: dx(*)
      integer,intent(in) :: incx
      real(wp),intent(in) :: dy(*)
      integer,intent(in) :: incy

      real(wp) :: dtemp
      integer :: i , ix , iy , m , mp1

      ddot = zero
      dtemp = zero
      if ( n<=0 ) return
      if ( incx==1 .and. incy==1 ) then

         ! code for both increments equal to 1

         ! clean-up loop

         m = mod(n,5)
         if ( m/=0 ) then
            do i = 1 , m
               dtemp = dtemp + dx(i)*dy(i)
            end do
            if ( n<5 ) then
               ddot = dtemp
               return
            end if
         end if
         mp1 = m + 1
         do i = mp1 , n , 5
            dtemp = dtemp + dx(i)*dy(i) + dx(i+1)*dy(i+1) + &
                    dx(i+2)*dy(i+2) + dx(i+3)*dy(i+3) + dx(i+4)*dy(i+4)
         end do
         ddot = dtemp

      else

         ! code for unequal increments or equal increments
         ! not equal to 1

         ix = 1
         iy = 1
         if ( incx<0 ) ix = (-n+1)*incx + 1
         if ( incy<0 ) iy = (-n+1)*incy + 1
         do i = 1 , n
            dtemp = dtemp + dx(ix)*dy(iy)
            ix = ix + incx
            iy = iy + incy
         end do
         ddot = dtemp

      end if

      end function ddot
!*******************************************************************************

!*******************************************************************************
!>
!  Function that returns the Euclidean norm
!  \( \sqrt{ \mathbf{x}^T \mathbf{x} } \) of a vector \( \mathbf{x} \).
!
!### Further details
!
!  * this version written on 25-october-1982.
!  * modified on 14-october-1993 to inline the call to dlassq.
!    sven hammarling, nag ltd.
!  * Converted to modern Fortran, Jacob Williams, Jan. 2016.
!
!@note Replaced original SLSQP routine with this one from
!      [BLAS](http://netlib.sandia.gov/blas/dnrm2.f).

    real(wp) function dnrm2(n,x,incx)

        implicit none

        integer,intent(in)               :: incx
        integer,intent(in)               :: n
        real(wp),dimension(*),intent(in) :: x

        real(wp) :: absxi , norm , scale , ssq
        integer :: ix

        if ( n<1 .or. incx<1 ) then
           norm = zero
        else if ( n==1 ) then
           norm = abs(x(1))
        else
           scale = zero
           ssq = one
           ! the following loop is equivalent to this call to the lapack
           ! auxiliary routine:
           ! call dlassq( n, x, incx, scale, ssq )
           do ix = 1 , 1 + (n-1)*incx , incx
              if ( x(ix)/=zero ) then
                 absxi = abs(x(ix))
                 if ( scale<absxi ) then
                    ssq = one + ssq*(scale/absxi)**2
                    scale = absxi
                 else
                    ssq = ssq + (absxi/scale)**2
                 end if
              end if
           end do
           norm = scale*sqrt(ssq)
        end if

        dnrm2 = norm

    end function dnrm2
!*******************************************************************************

!*******************************************************************************
!>
!  scales a vector by a constant.
!  uses unrolled loops for increment equal to one.
!
!### Author
!  jack dongarra, linpack, 3/11/78.

    subroutine dscal(n,da,dx,incx)

      implicit none

      integer,intent(in) :: n
      real(wp),intent(in) :: da 
      real(wp),intent(inout) :: dx(*)
      integer,intent(in) :: incx

      integer :: i , m , mp1 , nincx

      if ( n<=0 .or. incx<=0 ) return
      if ( incx==1 ) then

         ! code for increment equal to 1

         ! clean-up loop

         m = mod(n,5)
         if ( m/=0 ) then
            do i = 1 , m
               dx(i) = da*dx(i)
            end do
            if ( n<5 ) return
         end if
         mp1 = m + 1
         do i = mp1 , n , 5
            dx(i) = da*dx(i)
            dx(i+1) = da*dx(i+1)
            dx(i+2) = da*dx(i+2)
            dx(i+3) = da*dx(i+3)
            dx(i+4) = da*dx(i+4)
         end do
      else

         ! code for increment not equal to 1

         nincx = n*incx
         do i = 1 , nincx , incx
            dx(i) = da*dx(i)
         end do

      end if

      end subroutine dscal
!*******************************************************************************

!*******************************************************************************
!>
!  Find the smallest index of that component of a vector
!  having the maximum magnitude.
!
!### Description
!  Find smallest index of maximum magnitude of double precision DX.
!  IDAMAX = first I, I = 1 to N, to maximize ABS(DX(IX+(I-1)*INCX)),
!  where IX = 1 if INCX >= 0, else IX = 1+(1-N)*INCX.
!
!### Authors
!  * Lawson, C. L., (JPL)
!  * Hanson, R. J., (SNLA)
!  * Kincaid, D. R., (U. of Texas)
!  * Krogh, F. T., (JPL)
!
!### Reference
!  * C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!    Krogh, Basic linear algebra subprograms for Fortran
!    usage, Algorithm No. 539, Transactions on Mathematical
!    Software 5, 3 (September 1979), pp. 308-323.
!
!### History
!  * 791001  DATE WRITTEN
!  * 890531  Changed all specific intrinsics to generic.  (WRB)
!  * 890531  REVISION DATE from Version 3.2
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 900821  Modified to correct problem with a negative increment. (WRB)
!  * 920501  Reformatted the REFERENCES section.  (WRB)
!  * Jacob Williams, 2/21/2016, converted to modern Fortran.

    integer function idamax(n,dx,incx)

    implicit none

    integer,intent(in)               :: n     !! number of elements in input vector(s)
    real(wp),dimension(*),intent(in) :: dx    !! double precision vector with `N` elements
    integer,intent(in)               :: incx  !! storage spacing between elements of `DX`

    real(wp) :: dmax, xmag
    integer :: i, ix

    idamax = 0
    if ( n<=0 ) return
    idamax = 1
    if ( n==1 ) return

    if ( incx==1 ) then

        ! code for increments equal to 1.

        dmax = abs(dx(1))
        do i = 2 , n
            xmag = abs(dx(i))
            if ( xmag>dmax ) then
                idamax = i
                dmax = xmag
            endif
        enddo

    else

        ! code for increments not equal to 1.

        ix = 1
        if ( incx<0 ) ix = (-n+1)*incx + 1
        dmax = abs(dx(ix))
        ix = ix + incx
        do i = 2 , n
            xmag = abs(dx(ix))
            if ( xmag>dmax ) then
                idamax = i
                dmax = xmag
            endif
            ix = ix + incx
        enddo

    end if

    end function idamax
!*******************************************************************************

#endif

!*******************************************************************************
    end module dvode_blas_module
!*******************************************************************************
