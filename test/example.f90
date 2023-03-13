!*****************************************************************************************
!>
!  Example problem from the DVODE comment block.
! 
!  the following is a simple example problem, with the coding
!  needed for its solution by dvode.  the problem is from chemical
!  kinetics, and consists of the following three rate equations:
!```
!  dy1/dt = -.04*y1 + 1.e4*y2*y3
!  dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2
!  dy3/dt = 3.e7*y2**2
!```
!  on the interval from t = 0.0 to t = 4.e10, with initial conditions
!  y1 = 1.0, y2 = y3 = 0.  the problem is stiff.
! 
!  the following coding solves this problem with dvode, using mf = 21
!  and printing results at t = .4, 4., ..., 4.e10.  it uses
!  itol = 2 and atol much smaller for y2 than y1 or y3 because
!  y2 has much smaller values.
!  at the end of the run, statistical quantities of interest are
!  printed. (see optional output in the full description below.)
!  to generate fortran source code, replace c in column 1 with a blank
!  in the coding below.

    program dvode_example

    use dvode_module
    use dvode_kinds_module, only: wp => dvode_wp
    use iso_fortran_env, only: output_unit

    implicit none

    real(wp) :: atol(3), rtol(1), rwork(67), t, tout, y(3)
    integer :: neq,iwork(33),itol,itask,istate,iopt,lrw,liw,mf,iout
    type(dvode_t) :: solver

    neq = 3
    y(1) = 1.0_wp
    y(2) = 0.0_wp
    y(3) = 0.0_wp
    t = 0.0_wp
    tout = 0.4_wp
    itol = 2
    rtol = 1.0e-4_wp
    atol(1) = 1.0e-8_wp
    atol(2) = 1.0e-14_wp
    atol(3) = 1.0e-6_wp
    itask = 1
    istate = 1
    iopt = 0
    lrw = 67
    liw = 33
    mf = 21

    call solver%initialize(f=fex, jac=jex)
    do iout = 1,12
      call solver%solve(neq,y,t,tout,itol,rtol,atol,itask,istate,&
                        iopt,rwork,lrw,iwork,liw,mf)
      write(output_unit,'(a,d12.4,a,3d14.6)') &
                ' at t =',t,'   y =',y(1),y(2),y(3)
      if (istate < 0) then
        write(output_unit,'(///a,i3)') ' error halt: istate =',istate
        stop
      end if
      tout = tout*10.0_wp
    end do

    write(output_unit,'(/a,i4,a,i4,a,i4,*(a,i4/))') &
            ' no. steps =', iwork(11),&
            '   no. f-s =', iwork(12),&
            '   no. j-s =', iwork(13),&
            '   no. lu-s =',iwork(19),&
            '  no. nonlinear iterations =', iwork(20),&
            '  no. nonlinear convergence failures =', iwork(21),&
            '  no. error test failures =', iwork(22)
    
    contains

    subroutine fex (me, neq, t, y, ydot)

    class(dvode_t),intent(inout) :: me
    integer :: neq
    real(wp) :: t
    real(wp) :: y(neq)
    real(wp) :: ydot(neq)

    ydot(1) = -0.04_wp*y(1) + 1.0e4_wp*y(2)*y(3)
    ydot(3) = 3.0e7_wp*y(2)*y(2)
    ydot(2) = -ydot(1) - ydot(3)

    end subroutine fex

    subroutine jex (me, neq, t, y, ml, mu, pd, nrowpd)

    class(dvode_t),intent(inout) :: me
    integer :: neq
    integer :: nrowpd
    real(wp) :: t, y(neq), pd(nrowpd,neq)
    integer :: ml, mu

    pd(1,1) = -0.04_wp
    pd(1,2) = 1.0e4_wp*y(3)
    pd(1,3) = 1.0e4_wp*y(2)
    pd(2,1) = 0.04_wp
    pd(2,3) = -pd(1,3)
    pd(3,2) = 6.0e7_wp*y(2)
    pd(2,2) = -pd(1,2) - pd(3,2)
    
    end subroutine jex

    end program dvode_example

! the following output was obtained from the above program on a
! cray-1 computer with the cft compiler.
!
! at t =  4.0000e-01   y =  9.851680e-01  3.386314e-05  1.479817e-02
! at t =  4.0000e+00   y =  9.055255e-01  2.240539e-05  9.445214e-02
! at t =  4.0000e+01   y =  7.158108e-01  9.184883e-06  2.841800e-01
! at t =  4.0000e+02   y =  4.505032e-01  3.222940e-06  5.494936e-01
! at t =  4.0000e+03   y =  1.832053e-01  8.942690e-07  8.167938e-01
! at t =  4.0000e+04   y =  3.898560e-02  1.621875e-07  9.610142e-01
! at t =  4.0000e+05   y =  4.935882e-03  1.984013e-08  9.950641e-01
! at t =  4.0000e+06   y =  5.166183e-04  2.067528e-09  9.994834e-01
! at t =  4.0000e+07   y =  5.201214e-05  2.080593e-10  9.999480e-01
! at t =  4.0000e+08   y =  5.213149e-06  2.085271e-11  9.999948e-01
! at t =  4.0000e+09   y =  5.183495e-07  2.073399e-12  9.999995e-01
! at t =  4.0000e+10   y =  5.450996e-08  2.180399e-13  9.999999e-01
!
! no. steps = 595   no. f-s = 832   no. j-s =  13   no. lu-s = 112
!  no. nonlinear iterations = 831
!  no. nonlinear convergence failures =   0
!  no. error test failures =  22