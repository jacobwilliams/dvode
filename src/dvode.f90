module dvode_module

    use dvode_kinds
    use blas_module
    use linpack_module
    use iso_fortran_env,    only: output_unit

    implicit none

    private

    real(wp),parameter :: epmach = epsilon(1.0_wp) !! machine epsilon
    integer,parameter :: iumach  = output_unit     !! standard output unit number

    public :: dvode

contains


!-----------------------------------------------------------------------
! dvode: variable-coefficient ordinary differential equation solver,
! with fixed-leading-coefficient implementation.
! this version is in double precision.
!
! dvode solves the initial value problem for stiff or nonstiff
! systems of first order odes,
!     dy/dt = f(t,y) ,  or, in component form,
!     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(neq)) (i = 1,...,neq).
! dvode is a package based on the episode and episodeb packages, and
! on the odepack user interface standard, with minor modifications.
!-----------------------------------------------------------------------
! authors:
!               peter n. brown and alan c. hindmarsh
!               center for applied scientific computing, l-561
!               lawrence livermore national laboratory
!               livermore, ca 94551
! and
!               george d. byrne
!               illinois institute of technology
!               chicago, il 60616
!-----------------------------------------------------------------------
! references:
!
! 1. p. n. brown, g. d. byrne, and a. c. hindmarsh, "vode: a variable
!    coefficient ode solver," siam j. sci. stat. comput., 10 (1989),
!    pp. 1038-1051.  also, llnl report ucrl-98412, june 1988.
! 2. g. d. byrne and a. c. hindmarsh, "a polyalgorithm for the
!    numerical solution of ordinary differential equations,"
!    acm trans. math. software, 1 (1975), pp. 71-96.
! 3. a. c. hindmarsh and g. d. byrne, "episode: an effective package
!    for the integration of systems of ordinary differential
!    equations," llnl report ucid-30112, rev. 1, april 1977.
! 4. g. d. byrne and a. c. hindmarsh, "episodeb: an experimental
!    package for the integration of systems of ordinary differential
!    equations with banded jacobians," llnl report ucid-30132, april
!    1976.
! 5. a. c. hindmarsh, "odepack, a systematized collection of ode
!    solvers," in scientific computing, r. s. stepleman et al., eds.,
!    north-holland, amsterdam, 1983, pp. 55-64.
! 6. k. r. jackson and r. sacks-davis, "an alternative implementation
!    of variable step-size multistep formulas for stiff odes," acm
!    trans. math. software, 6 (1980), pp. 295-318.
!-----------------------------------------------------------------------
! summary of usage.
!
! communication between the user and the dvode package, for normal
! situations, is summarized here.  this summary describes only a subset
! of the full set of options available.  see the full description for
! details, including optional communication, nonstandard options,
! and instructions for special situations.  see also the example
! problem (with program and output) following this summary.
!
! a. first provide a subroutine of the form:
!           subroutine f (neq, t, y, ydot, rpar, ipar)
!           double precision t, y(neq), ydot(neq), rpar
! which supplies the vector function f by loading ydot(i) with f(i).
!
! b. next determine (or guess) whether or not the problem is stiff.
! stiffness occurs when the jacobian matrix df/dy has an eigenvalue
! whose real part is negative and large in magnitude, compared to the
! reciprocal of the t span of interest.  if the problem is nonstiff,
! use a method flag mf = 10.  if it is stiff, there are four standard
! choices for mf (21, 22, 24, 25), and dvode requires the jacobian
! matrix in some form.  in these cases (mf .gt. 0), dvode will use a
! saved copy of the jacobian matrix.  if this is undesirable because of
! storage limitations, set mf to the corresponding negative value
! (-21, -22, -24, -25).  (see full description of mf below.)
! the jacobian matrix is regarded either as full (mf = 21 or 22),
! or banded (mf = 24 or 25).  in the banded case, dvode requires two
! half-bandwidth parameters ml and mu.  these are, respectively, the
! widths of the lower and upper parts of the band, excluding the main
! diagonal.  thus the band consists of the locations (i,j) with
! i-ml .le. j .le. i+mu, and the full bandwidth is ml+mu+1.
!
! c. if the problem is stiff, you are encouraged to supply the jacobian
! directly (mf = 21 or 24), but if this is not feasible, dvode will
! compute it internally by difference quotients (mf = 22 or 25).
! if you are supplying the jacobian, provide a subroutine of the form:
!           subroutine jac (neq, t, y, ml, mu, pd, nrowpd, rpar, ipar)
!           double precision t, y(neq), pd(nrowpd,neq), rpar
! which supplies df/dy by loading pd as follows:
!     for a full jacobian (mf = 21), load pd(i,j) with df(i)/dy(j),
! the partial derivative of f(i) with respect to y(j).  (ignore the
! ml and mu arguments in this case.)
!     for a banded jacobian (mf = 24), load pd(i-j+mu+1,j) with
! df(i)/dy(j), i.e. load the diagonal lines of df/dy into the rows of
! pd from the top down.
!     in either case, only nonzero elements need be loaded.
!
! d. write a main program which calls subroutine dvode once for
! each point at which answers are desired.  this should also provide
! for possible use of logical unit 6 for output of error messages
! by dvode.  on the first call to dvode, supply arguments as follows:
! f      = name of subroutine for right-hand side vector f.
!          this name must be declared external in calling program.
! neq    = number of first order odes.
! y      = array of initial values, of length neq.
! t      = the initial value of the independent variable.
! tout   = first point where output is desired (.ne. t).
! itol   = 1 or 2 according as atol (below) is a scalar or array.
! rtol   = relative tolerance parameter (scalar).
! atol   = absolute tolerance parameter (scalar or array).
!          the estimated local error in y(i) will be controlled so as
!          to be roughly less (in magnitude) than
!             ewt(i) = rtol*abs(y(i)) + atol     if itol = 1, or
!             ewt(i) = rtol*abs(y(i)) + atol(i)  if itol = 2.
!          thus the local error test passes if, in each component,
!          either the absolute error is less than atol (or atol(i)),
!          or the relative error is less than rtol.
!          use rtol = 0.0 for pure absolute error control, and
!          use atol = 0.0 (or atol(i) = 0.0) for pure relative error
!          control.  caution: actual (global) errors may exceed these
!          local tolerances, so choose them conservatively.
! itask  = 1 for normal computation of output values of y at t = tout.
! istate = integer flag (input and output).  set istate = 1.
! iopt   = 0 to indicate no optional input used.
! rwork  = real work array of length at least:
!             20 + 16*neq                      for mf = 10,
!             22 +  9*neq + 2*neq**2           for mf = 21 or 22,
!             22 + 11*neq + (3*ml + 2*mu)*neq  for mf = 24 or 25.
! lrw    = declared length of rwork (in user's dimension statement).
! iwork  = integer work array of length at least:
!             30        for mf = 10,
!             30 + neq  for mf = 21, 22, 24, or 25.
!          if mf = 24 or 25, input in iwork(1),iwork(2) the lower
!          and upper half-bandwidths ml,mu.
! liw    = declared length of iwork (in user's dimension statement).
! jac    = name of subroutine for jacobian matrix (mf = 21 or 24).
!          if used, this name must be declared external in calling
!          program.  if not used, pass a dummy name.
! mf     = method flag.  standard values are:
!          10 for nonstiff (adams) method, no jacobian used.
!          21 for stiff (bdf) method, user-supplied full jacobian.
!          22 for stiff method, internally generated full jacobian.
!          24 for stiff method, user-supplied banded jacobian.
!          25 for stiff method, internally generated banded jacobian.
! rpar,ipar = user-defined real and integer arrays passed to f and jac.
! note that the main program must declare arrays y, rwork, iwork,
! and possibly atol, rpar, and ipar.
!
! e. the output from the first call (or any call) is:
!      y = array of computed values of y(t) vector.
!      t = corresponding value of independent variable (normally tout).
! istate = 2  if dvode was successful, negative otherwise.
!          -1 means excess work done on this call. (perhaps wrong mf.)
!          -2 means excess accuracy requested. (tolerances too small.)
!          -3 means illegal input detected. (see printed message.)
!          -4 means repeated error test failures. (check all input.)
!          -5 means repeated convergence failures. (perhaps bad
!             jacobian supplied or wrong choice of mf or tolerances.)
!          -6 means error weight became zero during problem. (solution
!             component i vanished, and atol or atol(i) = 0.)
!
! f. to continue the integration after a successful return, simply
! reset tout and call dvode again.  no other parameters need be reset.
!
!-----------------------------------------------------------------------
! example problem
!
! the following is a simple example problem, with the coding
! needed for its solution by dvode.  the problem is from chemical
! kinetics, and consists of the following three rate equations:
!     dy1/dt = -.04*y1 + 1.e4*y2*y3
!     dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2
!     dy3/dt = 3.e7*y2**2
! on the interval from t = 0.0 to t = 4.e10, with initial conditions
! y1 = 1.0, y2 = y3 = 0.  the problem is stiff.
!
! the following coding solves this problem with dvode, using mf = 21
! and printing results at t = .4, 4., ..., 4.e10.  it uses
! itol = 2 and atol much smaller for y2 than y1 or y3 because
! y2 has much smaller values.
! at the end of the run, statistical quantities of interest are
! printed. (see optional output in the full description below.)
! to generate fortran source code, replace c in column 1 with a blank
! in the coding below.
!
!     external fex, jex
!     double precision atol, rpar, rtol, rwork, t, tout, y
!     dimension y(3), atol(3), rwork(67), iwork(33)
!     neq = 3
!     y(1) = 1.0d0
!     y(2) = 0.0d0
!     y(3) = 0.0d0
!     t = 0.0d0
!     tout = 0.4d0
!     itol = 2
!     rtol = 1.d-4
!     atol(1) = 1.d-8
!     atol(2) = 1.d-14
!     atol(3) = 1.d-6
!     itask = 1
!     istate = 1
!     iopt = 0
!     lrw = 67
!     liw = 33
!     mf = 21
!     do 40 iout = 1,12
!       call dvode(fex,neq,y,t,tout,itol,rtol,atol,itask,istate,
!    1            iopt,rwork,lrw,iwork,liw,jex,mf,rpar,ipar)
!       write(6,20)t,y(1),y(2),y(3)
! 20    format(' at t =',d12.4,'   y =',3d14.6)
!       if (istate .lt. 0) go to 80
! 40    tout = tout*10.
!     write(6,60) iwork(11),iwork(12),iwork(13),iwork(19),
!    1            iwork(20),iwork(21),iwork(22)
! 60  format(/' no. steps =',i4,'   no. f-s =',i4,
!    1       '   no. j-s =',i4,'   no. lu-s =',i4/
!    2       '  no. nonlinear iterations =',i4/
!    3       '  no. nonlinear convergence failures =',i4/
!    4       '  no. error test failures =',i4/)
!     stop
! 80  write(6,90)istate
! 90  format(///' error halt: istate =',i3)
!     stop
!     end
!
!     subroutine fex (neq, t, y, ydot, rpar, ipar)
!     double precision rpar, t, y, ydot
!     dimension y(neq), ydot(neq)
!     ydot(1) = -.04d0*y(1) + 1.d4*y(2)*y(3)
!     ydot(3) = 3.d7*y(2)*y(2)
!     ydot(2) = -ydot(1) - ydot(3)
!     return
!     end
!
!     subroutine jex (neq, t, y, ml, mu, pd, nrpd, rpar, ipar)
!     double precision pd, rpar, t, y
!     dimension y(neq), pd(nrpd,neq)
!     pd(1,1) = -.04d0
!     pd(1,2) = 1.d4*y(3)
!     pd(1,3) = 1.d4*y(2)
!     pd(2,1) = .04d0
!     pd(2,3) = -pd(1,3)
!     pd(3,2) = 6.d7*y(2)
!     pd(2,2) = -pd(1,2) - pd(3,2)
!     return
!     end
!
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
!-----------------------------------------------------------------------
! full description of user interface to dvode.
!
! the user interface to dvode consists of the following parts.
!
! i.   the call sequence to subroutine dvode, which is a driver
!      routine for the solver.  this includes descriptions of both
!      the call sequence arguments and of user-supplied routines.
!      following these descriptions is
!        * a description of optional input available through the
!          call sequence,
!        * a description of optional output (in the work arrays), and
!        * instructions for interrupting and restarting a solution.
!
! ii.  descriptions of other routines in the dvode package that may be
!      (optionally) called by the user.  these provide the ability to
!      alter error message handling, save and restore the internal
!      common, and obtain specified derivatives of the solution y(t).
!
! iii. descriptions of common blocks to be declared in overlay
!      or similar environments.
!
! iv.  description of two routines in the dvode package, either of
!      which the user may replace with his own version, if desired.
!      these relate to the measurement of errors.
!
!-----------------------------------------------------------------------
! part i.  call sequence.
!
! the call sequence parameters used for input only are
!     f, neq, tout, itol, rtol, atol, itask, iopt, lrw, liw, jac, mf,
! and those used for both input and output are
!     y, t, istate.
! the work arrays rwork and iwork are also used for conditional and
! optional input and optional output.  (the term output here refers
! to the return from subroutine dvode to the user's calling program.)
!
! the legality of input parameters will be thoroughly checked on the
! initial call for the problem, but not checked thereafter unless a
! change in input parameters is flagged by istate = 3 in the input.
!
! the descriptions of the call arguments are as follows.
!
! f      = the name of the user-supplied subroutine defining the
!          ode system.  the system must be put in the first-order
!          form dy/dt = f(t,y), where f is a vector-valued function
!          of the scalar t and the vector y.  subroutine f is to
!          compute the function f.  it is to have the form
!               subroutine f (neq, t, y, ydot, rpar, ipar)
!               double precision t, y(neq), ydot(neq), rpar
!          where neq, t, and y are input, and the array ydot = f(t,y)
!          is output.  y and ydot are arrays of length neq.
!          subroutine f should not alter y(1),...,y(neq).
!          f must be declared external in the calling program.
!
!          subroutine f may access user-defined real and integer
!          work arrays rpar and ipar, which are to be dimensioned
!          in the main program.
!
!          if quantities computed in the f routine are needed
!          externally to dvode, an extra call to f should be made
!          for this purpose, for consistent and accurate results.
!          if only the derivative dy/dt is needed, use dvindy instead.
!
! neq    = the size of the ode system (number of first order
!          ordinary differential equations).  used only for input.
!          neq may not be increased during the problem, but
!          can be decreased (with istate = 3 in the input).
!
! y      = a real array for the vector of dependent variables, of
!          length neq or more.  used for both input and output on the
!          first call (istate = 1), and only for output on other calls.
!          on the first call, y must contain the vector of initial
!          values.  in the output, y contains the computed solution
!          evaluated at t.  if desired, the y array may be used
!          for other purposes between calls to the solver.
!
!          this array is passed as the y argument in all calls to
!          f and jac.
!
! t      = the independent variable.  in the input, t is used only on
!          the first call, as the initial point of the integration.
!          in the output, after each call, t is the value at which a
!          computed solution y is evaluated (usually the same as tout).
!          on an error return, t is the farthest point reached.
!
! tout   = the next value of t at which a computed solution is desired.
!          used only for input.
!
!          when starting the problem (istate = 1), tout may be equal
!          to t for one call, then should .ne. t for the next call.
!          for the initial t, an input value of tout .ne. t is used
!          in order to determine the direction of the integration
!          (i.e. the algebraic sign of the step sizes) and the rough
!          scale of the problem.  integration in either direction
!          (forward or backward in t) is permitted.
!
!          if itask = 2 or 5 (one-step modes), tout is ignored after
!          the first call (i.e. the first call with tout .ne. t).
!          otherwise, tout is required on every call.
!
!          if itask = 1, 3, or 4, the values of tout need not be
!          monotone, but a value of tout which backs up is limited
!          to the current internal t interval, whose endpoints are
!          tcur - hu and tcur.  (see optional output, below, for
!          tcur and hu.)
!
! itol   = an indicator for the type of error control.  see
!          description below under atol.  used only for input.
!
! rtol   = a relative error tolerance parameter, either a scalar or
!          an array of length neq.  see description below under atol.
!          input only.
!
! atol   = an absolute error tolerance parameter, either a scalar or
!          an array of length neq.  input only.
!
!          the input parameters itol, rtol, and atol determine
!          the error control performed by the solver.  the solver will
!          control the vector e = (e(i)) of estimated local errors
!          in y, according to an inequality of the form
!                      rms-norm of ( e(i)/ewt(i) )   .le.   1,
!          where       ewt(i) = rtol(i)*abs(y(i)) + atol(i),
!          and the rms-norm (root-mean-square norm) here is
!          rms-norm(v) = sqrt(sum v(i)**2 / neq).  here ewt = (ewt(i))
!          is a vector of weights which must always be positive, and
!          the values of rtol and atol should all be non-negative.
!          the following table gives the types (scalar/array) of
!          rtol and atol, and the corresponding form of ewt(i).
!
!             itol    rtol       atol          ewt(i)
!              1     scalar     scalar     rtol*abs(y(i)) + atol
!              2     scalar     array      rtol*abs(y(i)) + atol(i)
!              3     array      scalar     rtol(i)*abs(y(i)) + atol
!              4     array      array      rtol(i)*abs(y(i)) + atol(i)
!
!          when either of these parameters is a scalar, it need not
!          be dimensioned in the user's calling program.
!
!          if none of the above choices (with itol, rtol, and atol
!          fixed throughout the problem) is suitable, more general
!          error controls can be obtained by substituting
!          user-supplied routines for the setting of ewt and/or for
!          the norm calculation.  see part iv below.
!
!          if global errors are to be estimated by making a repeated
!          run on the same problem with smaller tolerances, then all
!          components of rtol and atol (i.e. of ewt) should be scaled
!          down uniformly.
!
! itask  = an index specifying the task to be performed.
!          input only.  itask has the following values and meanings.
!          1  means normal computation of output values of y(t) at
!             t = tout (by overshooting and interpolating).
!          2  means take one step only and return.
!          3  means stop at the first internal mesh point at or
!             beyond t = tout and return.
!          4  means normal computation of output values of y(t) at
!             t = tout but without overshooting t = tcrit.
!             tcrit must be input as rwork(1).  tcrit may be equal to
!             or beyond tout, but not behind it in the direction of
!             integration.  this option is useful if the problem
!             has a singularity at or beyond t = tcrit.
!          5  means take one step, without passing tcrit, and return.
!             tcrit must be input as rwork(1).
!
!          note:  if itask = 4 or 5 and the solver reaches tcrit
!          (within roundoff), it will return t = tcrit (exactly) to
!          indicate this (unless itask = 4 and tout comes before tcrit,
!          in which case answers at t = tout are returned first).
!
! istate = an index used for input and output to specify the
!          the state of the calculation.
!
!          in the input, the values of istate are as follows.
!          1  means this is the first call for the problem
!             (initializations will be done).  see note below.
!          2  means this is not the first call, and the calculation
!             is to continue normally, with no change in any input
!             parameters except possibly tout and itask.
!             (if itol, rtol, and/or atol are changed between calls
!             with istate = 2, the new values will be used but not
!             tested for legality.)
!          3  means this is not the first call, and the
!             calculation is to continue normally, but with
!             a change in input parameters other than
!             tout and itask.  changes are allowed in
!             neq, itol, rtol, atol, iopt, lrw, liw, mf, ml, mu,
!             and any of the optional input except h0.
!             (see iwork description for ml and mu.)
!          note:  a preliminary call with tout = t is not counted
!          as a first call here, as no initialization or checking of
!          input is done.  (such a call is sometimes useful to include
!          the initial conditions in the output.)
!          thus the first call for which tout .ne. t requires
!          istate = 1 in the input.
!
!          in the output, istate has the following values and meanings.
!           1  means nothing was done, as tout was equal to t with
!              istate = 1 in the input.
!           2  means the integration was performed successfully.
!          -1  means an excessive amount of work (more than mxstep
!              steps) was done on this call, before completing the
!              requested task, but the integration was otherwise
!              successful as far as t.  (mxstep is an optional input
!              and is normally 500.)  to continue, the user may
!              simply reset istate to a value .gt. 1 and call again.
!              (the excess work step counter will be reset to 0.)
!              in addition, the user may increase mxstep to avoid
!              this error return.  (see optional input below.)
!          -2  means too much accuracy was requested for the precision
!              of the machine being used.  this was detected before
!              completing the requested task, but the integration
!              was successful as far as t.  to continue, the tolerance
!              parameters must be reset, and istate must be set
!              to 3.  the optional output tolsf may be used for this
!              purpose.  (note: if this condition is detected before
!              taking any steps, then an illegal input return
!              (istate = -3) occurs instead.)
!          -3  means illegal input was detected, before taking any
!              integration steps.  see written message for details.
!              note:  if the solver detects an infinite loop of calls
!              to the solver with illegal input, it will cause
!              the run to stop.
!          -4  means there were repeated error test failures on
!              one attempted step, before completing the requested
!              task, but the integration was successful as far as t.
!              the problem may have a singularity, or the input
!              may be inappropriate.
!          -5  means there were repeated convergence test failures on
!              one attempted step, before completing the requested
!              task, but the integration was successful as far as t.
!              this may be caused by an inaccurate jacobian matrix,
!              if one is being used.
!          -6  means ewt(i) became zero for some i during the
!              integration.  pure relative error control (atol(i)=0.0)
!              was requested on a variable which has now vanished.
!              the integration was successful as far as t.
!
!          note:  since the normal output value of istate is 2,
!          it does not need to be reset for normal continuation.
!          also, since a negative input value of istate will be
!          regarded as illegal, a negative output value requires the
!          user to change it, and possibly other input, before
!          calling the solver again.
!
! iopt   = an integer flag to specify whether or not any optional
!          input is being used on this call.  input only.
!          the optional input is listed separately below.
!          iopt = 0 means no optional input is being used.
!                   default values will be used in all cases.
!          iopt = 1 means optional input is being used.
!
! rwork  = a real working array (double precision).
!          the length of rwork must be at least
!             20 + nyh*(maxord + 1) + 3*neq + lwm    where
!          nyh    = the initial value of neq,
!          maxord = 12 (if meth = 1) or 5 (if meth = 2) (unless a
!                   smaller value is given as an optional input),
!          lwm = length of work space for matrix-related data:
!          lwm = 0             if miter = 0,
!          lwm = 2*neq**2 + 2  if miter = 1 or 2, and mf.gt.0,
!          lwm = neq**2 + 2    if miter = 1 or 2, and mf.lt.0,
!          lwm = neq + 2       if miter = 3,
!          lwm = (3*ml+2*mu+2)*neq + 2 if miter = 4 or 5, and mf.gt.0,
!          lwm = (2*ml+mu+1)*neq + 2   if miter = 4 or 5, and mf.lt.0.
!          (see the mf description for meth and miter.)
!          thus if maxord has its default value and neq is constant,
!          this length is:
!             20 + 16*neq                    for mf = 10,
!             22 + 16*neq + 2*neq**2         for mf = 11 or 12,
!             22 + 16*neq + neq**2           for mf = -11 or -12,
!             22 + 17*neq                    for mf = 13,
!             22 + 18*neq + (3*ml+2*mu)*neq  for mf = 14 or 15,
!             22 + 17*neq + (2*ml+mu)*neq    for mf = -14 or -15,
!             20 +  9*neq                    for mf = 20,
!             22 +  9*neq + 2*neq**2         for mf = 21 or 22,
!             22 +  9*neq + neq**2           for mf = -21 or -22,
!             22 + 10*neq                    for mf = 23,
!             22 + 11*neq + (3*ml+2*mu)*neq  for mf = 24 or 25.
!             22 + 10*neq + (2*ml+mu)*neq    for mf = -24 or -25.
!          the first 20 words of rwork are reserved for conditional
!          and optional input and optional output.
!
!          the following word in rwork is a conditional input:
!            rwork(1) = tcrit = critical value of t which the solver
!                       is not to overshoot.  required if itask is
!                       4 or 5, and ignored otherwise.  (see itask.)
!
! lrw    = the length of the array rwork, as declared by the user.
!          (this will be checked by the solver.)
!
! iwork  = an integer work array.  the length of iwork must be at least
!             30        if miter = 0 or 3 (mf = 10, 13, 20, 23), or
!             30 + neq  otherwise (abs(mf) = 11,12,14,15,21,22,24,25).
!          the first 30 words of iwork are reserved for conditional and
!          optional input and optional output.
!
!          the following 2 words in iwork are conditional input:
!            iwork(1) = ml     these are the lower and upper
!            iwork(2) = mu     half-bandwidths, respectively, of the
!                       banded jacobian, excluding the main diagonal.
!                       the band is defined by the matrix locations
!                       (i,j) with i-ml .le. j .le. i+mu.  ml and mu
!                       must satisfy  0 .le.  ml,mu  .le. neq-1.
!                       these are required if miter is 4 or 5, and
!                       ignored otherwise.  ml and mu may in fact be
!                       the band parameters for a matrix to which
!                       df/dy is only approximately equal.
!
! liw    = the length of the array iwork, as declared by the user.
!          (this will be checked by the solver.)
!
! note:  the work arrays must not be altered between calls to dvode
! for the same problem, except possibly for the conditional and
! optional input, and except for the last 3*neq words of rwork.
! the latter space is used for internal scratch space, and so is
! available for use by the user outside dvode between calls, if
! desired (but not for use by f or jac).
!
! jac    = the name of the user-supplied routine (miter = 1 or 4) to
!          compute the jacobian matrix, df/dy, as a function of
!          the scalar t and the vector y.  it is to have the form
!               subroutine jac (neq, t, y, ml, mu, pd, nrowpd,
!                               rpar, ipar)
!               double precision t, y(neq), pd(nrowpd,neq), rpar
!          where neq, t, y, ml, mu, and nrowpd are input and the array
!          pd is to be loaded with partial derivatives (elements of the
!          jacobian matrix) in the output.  pd must be given a first
!          dimension of nrowpd.  t and y have the same meaning as in
!          subroutine f.
!               in the full matrix case (miter = 1), ml and mu are
!          ignored, and the jacobian is to be loaded into pd in
!          columnwise manner, with df(i)/dy(j) loaded into pd(i,j).
!               in the band matrix case (miter = 4), the elements
!          within the band are to be loaded into pd in columnwise
!          manner, with diagonal lines of df/dy loaded into the rows
!          of pd. thus df(i)/dy(j) is to be loaded into pd(i-j+mu+1,j).
!          ml and mu are the half-bandwidth parameters. (see iwork).
!          the locations in pd in the two triangular areas which
!          correspond to nonexistent matrix elements can be ignored
!          or loaded arbitrarily, as they are overwritten by dvode.
!               jac need not provide df/dy exactly.  a crude
!          approximation (possibly with a smaller bandwidth) will do.
!               in either case, pd is preset to zero by the solver,
!          so that only the nonzero elements need be loaded by jac.
!          each call to jac is preceded by a call to f with the same
!          arguments neq, t, and y.  thus to gain some efficiency,
!          intermediate quantities shared by both calculations may be
!          saved in a user common block by f and not recomputed by jac,
!          if desired.  also, jac may alter the y array, if desired.
!          jac must be declared external in the calling program.
!               subroutine jac may access user-defined real and integer
!          work arrays, rpar and ipar, whose dimensions are set by the
!          user in the main program.
!
! mf     = the method flag.  used only for input.  the legal values of
!          mf are 10, 11, 12, 13, 14, 15, 20, 21, 22, 23, 24, 25,
!          -11, -12, -14, -15, -21, -22, -24, -25.
!          mf is a signed two-digit integer, mf = jsv*(10*meth + miter).
!          jsv = sign(mf) indicates the jacobian-saving strategy:
!            jsv =  1 means a copy of the jacobian is saved for reuse
!                     in the corrector iteration algorithm.
!            jsv = -1 means a copy of the jacobian is not saved
!                     (valid only for miter = 1, 2, 4, or 5).
!          meth indicates the basic linear multistep method:
!            meth = 1 means the implicit adams method.
!            meth = 2 means the method based on backward
!                     differentiation formulas (bdf-s).
!          miter indicates the corrector iteration method:
!            miter = 0 means functional iteration (no jacobian matrix
!                      is involved).
!            miter = 1 means chord iteration with a user-supplied
!                      full (neq by neq) jacobian.
!            miter = 2 means chord iteration with an internally
!                      generated (difference quotient) full jacobian
!                      (using neq extra calls to f per df/dy value).
!            miter = 3 means chord iteration with an internally
!                      generated diagonal jacobian approximation
!                      (using 1 extra call to f per df/dy evaluation).
!            miter = 4 means chord iteration with a user-supplied
!                      banded jacobian.
!            miter = 5 means chord iteration with an internally
!                      generated banded jacobian (using ml+mu+1 extra
!                      calls to f per df/dy evaluation).
!          if miter = 1 or 4, the user must supply a subroutine jac
!          (the name is arbitrary) as described above under jac.
!          for other values of miter, a dummy argument can be used.
!
! rpar     user-specified array used to communicate real parameters
!          to user-supplied subroutines.  if rpar is a vector, then
!          it must be dimensioned in the user's main program.  if it
!          is unused or it is a scalar, then it need not be
!          dimensioned.
!
! ipar     user-specified array used to communicate integer parameter
!          to user-supplied subroutines.  the comments on dimensioning
!          rpar apply to ipar.
!-----------------------------------------------------------------------
! optional input.
!
! the following is a list of the optional input provided for in the
! call sequence.  (see also part ii.)  for each such input variable,
! this table lists its name as used in this documentation, its
! location in the call sequence, its meaning, and the default value.
! the use of any of this input requires iopt = 1, and in that
! case all of this input is examined.  a value of zero for any
! of these optional input variables will cause the default value to be
! used.  thus to use a subset of the optional input, simply preload
! locations 5 to 10 in rwork and iwork to 0.0 and 0 respectively, and
! then set those of interest to nonzero values.
!
! name    location      meaning and default value
!
! h0      rwork(5)  the step size to be attempted on the first step.
!                   the default value is determined by the solver.
!
! hmax    rwork(6)  the maximum absolute step size allowed.
!                   the default value is infinite.
!
! hmin    rwork(7)  the minimum absolute step size allowed.
!                   the default value is 0.  (this lower bound is not
!                   enforced on the final step before reaching tcrit
!                   when itask = 4 or 5.)
!
! maxord  iwork(5)  the maximum order to be allowed.  the default
!                   value is 12 if meth = 1, and 5 if meth = 2.
!                   if maxord exceeds the default value, it will
!                   be reduced to the default value.
!                   if maxord is changed during the problem, it may
!                   cause the current order to be reduced.
!
! mxstep  iwork(6)  maximum number of (internally defined) steps
!                   allowed during one call to the solver.
!                   the default value is 500.
!
! mxhnil  iwork(7)  maximum number of messages printed (per problem)
!                   warning that t + h = t on a step (h = step size).
!                   this must be positive to result in a non-default
!                   value.  the default value is 10.
!
!-----------------------------------------------------------------------
! optional output.
!
! as optional additional output from dvode, the variables listed
! below are quantities related to the performance of dvode
! which are available to the user.  these are communicated by way of
! the work arrays, but also have internal mnemonic names as shown.
! except where stated otherwise, all of this output is defined
! on any successful return from dvode, and on any return with
! istate = -1, -2, -4, -5, or -6.  on an illegal input return
! (istate = -3), they will be unchanged from their existing values
! (if any), except possibly for tolsf, lenrw, and leniw.
! on any error return, output relevant to the error will be defined,
! as noted below.
!
! name    location      meaning
!
! hu      rwork(11) the step size in t last used (successfully).
!
! hcur    rwork(12) the step size to be attempted on the next step.
!
! tcur    rwork(13) the current value of the independent variable
!                   which the solver has actually reached, i.e. the
!                   current internal mesh point in t.  in the output,
!                   tcur will always be at least as far from the
!                   initial value of t as the current argument t,
!                   but may be farther (if interpolation was done).
!
! tolsf   rwork(14) a tolerance scale factor, greater than 1.0,
!                   computed when a request for too much accuracy was
!                   detected (istate = -3 if detected at the start of
!                   the problem, istate = -2 otherwise).  if itol is
!                   left unaltered but rtol and atol are uniformly
!                   scaled up by a factor of tolsf for the next call,
!                   then the solver is deemed likely to succeed.
!                   (the user may also ignore tolsf and alter the
!                   tolerance parameters in any other way appropriate.)
!
! nst     iwork(11) the number of steps taken for the problem so far.
!
! nfe     iwork(12) the number of f evaluations for the problem so far.
!
! nje     iwork(13) the number of jacobian evaluations so far.
!
! nqu     iwork(14) the method order last used (successfully).
!
! nqcur   iwork(15) the order to be attempted on the next step.
!
! imxer   iwork(16) the index of the component of largest magnitude in
!                   the weighted local error vector ( e(i)/ewt(i) ),
!                   on an error return with istate = -4 or -5.
!
! lenrw   iwork(17) the length of rwork actually required.
!                   this is defined on normal returns and on an illegal
!                   input return for insufficient storage.
!
! leniw   iwork(18) the length of iwork actually required.
!                   this is defined on normal returns and on an illegal
!                   input return for insufficient storage.
!
! nlu     iwork(19) the number of matrix lu decompositions so far.
!
! nni     iwork(20) the number of nonlinear (newton) iterations so far.
!
! ncfn    iwork(21) the number of convergence failures of the nonlinear
!                   solver so far.
!
! netf    iwork(22) the number of error test failures of the integrator
!                   so far.
!
! the following two arrays are segments of the rwork array which
! may also be of interest to the user as optional output.
! for each array, the table below gives its internal name,
! its base address in rwork, and its description.
!
! name    base address      description
!
! yh      21             the nordsieck history array, of size nyh by
!                        (nqcur + 1), where nyh is the initial value
!                        of neq.  for j = 0,1,...,nqcur, column j+1
!                        of yh contains hcur**j/factorial(j) times
!                        the j-th derivative of the interpolating
!                        polynomial currently representing the
!                        solution, evaluated at t = tcur.
!
! acor     lenrw-neq+1   array of size neq used for the accumulated
!                        corrections on each step, scaled in the output
!                        to represent the estimated local error in y
!                        on the last step.  this is the vector e in
!                        the description of the error control.  it is
!                        defined only on a successful return from dvode.
!
!-----------------------------------------------------------------------
! interrupting and restarting
!
! if the integration of a given problem by dvode is to be
! interrrupted and then later continued, such as when restarting
! an interrupted run or alternating between two or more ode problems,
! the user should save, following the return from the last dvode call
! prior to the interruption, the contents of the call sequence
! variables and internal common blocks, and later restore these
! values before the next dvode call for that problem.  to save
! and restore the common blocks, use subroutine dvsrco, as
! described below in part ii.
!
! in addition, if non-default values for either lun or mflag are
! desired, an extra call to xsetun and/or xsetf should be made just
! before continuing the integration.  see part ii below for details.
!
!-----------------------------------------------------------------------
! part ii.  other routines callable.
!
! the following are optional calls which the user may make to
! gain additional capabilities in conjunction with dvode.
! (the routines xsetun and xsetf are designed to conform to the
! slatec error handling package.)
!
!     form of call                  function
!  call xsetun(lun)           set the logical unit number, lun, for
!                             output of messages from dvode, if
!                             the default is not desired.
!                             the default value of lun is 6.
!
!  call xsetf(mflag)          set a flag to control the printing of
!                             messages by dvode.
!                             mflag = 0 means do not print. (danger:
!                             this risks losing valuable information.)
!                             mflag = 1 means print (the default).
!
!                             either of the above calls may be made at
!                             any time and will take effect immediately.
!
!  call dvsrco(rsav,isav,job) saves and restores the contents of
!                             the internal common blocks used by
!                             dvode. (see part iii below.)
!                             rsav must be a real array of length 49
!                             or more, and isav must be an integer
!                             array of length 40 or more.
!                             job=1 means save common into rsav/isav.
!                             job=2 means restore common from rsav/isav.
!                                dvsrco is useful if one is
!                             interrupting a run and restarting
!                             later, or alternating between two or
!                             more problems solved with dvode.
!
!  call dvindy(,,,,,)         provide derivatives of y, of various
!        (see below.)         orders, at a specified point t, if
!                             desired.  it may be called only after
!                             a successful return from dvode.
!
! the detailed instructions for using dvindy are as follows.
! the form of the call is:
!
!  call dvindy (t, k, rwork(21), nyh, dky, iflag)
!
! the input parameters are:
!
! t         = value of independent variable where answers are desired
!             (normally the same as the t last returned by dvode).
!             for valid results, t must lie between tcur - hu and tcur.
!             (see optional output for tcur and hu.)
! k         = integer order of the derivative desired.  k must satisfy
!             0 .le. k .le. nqcur, where nqcur is the current order
!             (see optional output).  the capability corresponding
!             to k = 0, i.e. computing y(t), is already provided
!             by dvode directly.  since nqcur .ge. 1, the first
!             derivative dy/dt is always available with dvindy.
! rwork(21) = the base address of the history array yh.
! nyh       = column length of yh, equal to the initial value of neq.
!
! the output parameters are:
!
! dky       = a real array of length neq containing the computed value
!             of the k-th derivative of y(t).
! iflag     = integer flag, returned as 0 if k and t were legal,
!             -1 if k was illegal, and -2 if t was illegal.
!             on an error return, a message is also written.
!-----------------------------------------------------------------------
! part iii.  common blocks.
! if dvode is to be used in an overlay situation, the user
! must declare, in the primary overlay, the variables in:
!   (1) the call sequence to dvode,
!   (2) the two internal common blocks
!         /dvod01/  of length  81  (48 double precision words
!                         followed by 33 integer words),
!         /dvod02/  of length  9  (1 double precision word
!                         followed by 8 integer words),
!
! if dvode is used on a system in which the contents of internal
! common blocks are not preserved between calls, the user should
! declare the above two common blocks in his main program to insure
! that their contents are preserved.
!
!-----------------------------------------------------------------------
! part iv.  optionally replaceable solver routines.
!
! below are descriptions of two routines in the dvode package which
! relate to the measurement of errors.  either routine can be
! replaced by a user-supplied version, if desired.  however, since such
! a replacement may have a major impact on performance, it should be
! done only when absolutely necessary, and only with great caution.
! (note: the means by which the package version of a routine is
! superseded by the user's version may be system-dependent.)
!
! (a) dewset.
! the following subroutine is called just before each internal
! integration step, and sets the array of error weights, ewt, as
! described under itol/rtol/atol above:
!     subroutine dewset (neq, itol, rtol, atol, ycur, ewt)
! where neq, itol, rtol, and atol are as in the dvode call sequence,
! ycur contains the current dependent variable vector, and
! ewt is the array of weights set by dewset.
!
! if the user supplies this subroutine, it must return in ewt(i)
! (i = 1,...,neq) a positive quantity suitable for comparison with
! errors in y(i).  the ewt array returned by dewset is passed to the
! dvnorm routine (see below.), and also used by dvode in the computation
! of the optional output imxer, the diagonal jacobian approximation,
! and the increments for difference quotient jacobians.
!
! in the user-supplied version of dewset, it may be desirable to use
! the current values of derivatives of y.  derivatives up to order nq
! are available from the history array yh, described above under
! optional output.  in dewset, yh is identical to the ycur array,
! extended to nq + 1 columns with a column length of nyh and scale
! factors of h**j/factorial(j).  on the first call for the problem,
! given by nst = 0, nq is 1 and h is temporarily set to 1.0.
! nyh is the initial value of neq.  the quantities nq, h, and nst
! can be obtained by including in dewset the statements:
!     double precision rvod, h, hu
!     common /dvod01/ rvod(48), ivod(33)
!     common /dvod02/ hu, ncfn, netf, nfe, nje, nlu, nni, nqu, nst
!     nq = ivod(28)
!     h = rvod(21)
! thus, for example, the current value of dy/dt can be obtained as
! ycur(nyh+i)/h  (i=1,...,neq)  (and the division by h is
! unnecessary when nst = 0).
!
! (b) dvnorm.
! the following is a real function routine which computes the weighted
! root-mean-square norm of a vector v:
!     d = dvnorm (n, v, w)
! where:
!   n = the length of the vector,
!   v = real array of length n containing the vector,
!   w = real array of length n containing weights,
!   d = sqrt( (1/n) * sum(v(i)*w(i))**2 ).
! dvnorm is called with n = neq and with w(i) = 1.0/ewt(i), where
! ewt is as set by subroutine dewset.
!
! if the user supplies this function, it should return a non-negative
! value of dvnorm suitable for use in the error control in dvode.
! none of the arguments should be altered by dvnorm.
! for example, a user-supplied dvnorm routine might:
!   -substitute a max-norm of (v(i)*w(i)) for the rms-norm, or
!   -ignore some components of v in the norm, with the effect of
!    suppressing the error control on those components of y.
!-----------------------------------------------------------------------
! revision history (yyyymmdd)
!  19890615  date written.  initial release.
!  19890922  added interrupt/restart ability, minor changes throughout.
!  19910228  minor revisions in line format,  prologue, etc.
!  19920227  modifications by d. pang:
!            (1) applied subgennam to get generic intrinsic names.
!            (2) changed intrinsic names to generic in comments.
!            (3) added *deck lines before each routine.
!  19920721  names of routines and labeled common blocks changed, so as
!            to be unique in combined single/double precision code (ach).
!  19920722  minor revisions to prologue (ach).
!  19920831  conversion to double precision done (ach).
!  19921106  fixed minor bug: etaq,etaqm1 in dvstep save statement (ach).
!  19921118  changed lunsav/mflgsv to ixsav (ach).
!  19941222  removed mf overwrite; attached sign to h in estimated second
!            deriv. in dvhin; misc. comment changes throughout (ach).
!  19970515  minor corrections to comments in prologue, dvjac (ach).
!  19981111  corrected block b by adding final line, go to 200 (ach).
!  20020430  various upgrades (ach): use odepack error handler package.
!            replaced d1mach by dumach.  various changes to main
!            prologue and other routine prologues.
!-----------------------------------------------------------------------
! other routines in the dvode package.
!
! in addition to subroutine dvode, the dvode package includes the
! following subroutines and function routines:
!  dvhin     computes an approximate step size for the initial step.
!  dvindy    computes an interpolated value of the y vector at t = tout.
!  dvstep    is the core integrator, which does one step of the
!            integration and the associated error control.
!  dvset     sets all method coefficients and test constants.
!  dvnlsd    solves the underlying nonlinear system -- the corrector.
!  dvjac     computes and preprocesses the jacobian matrix j = df/dy
!            and the newton iteration matrix p = i - (h/l1)*j.
!  dvsol     manages solution of linear system in chord iteration.
!  dvjust    adjusts the history array on a change of order.
!  dewset    sets the error weight vector ewt before each step.
!  dvnorm    computes the weighted r.m.s. norm of a vector.
!  dvsrco    is a user-callable routine to save and restore
!            the contents of the internal common blocks.
!  dacopy    is a routine to copy one two-dimensional array to another.
!  dgefa and dgesl   are routines from linpack for solving full
!            systems of linear algebraic equations.
!  dgbfa and dgbsl   are routines from linpack for solving banded
!            linear systems.
!  daxpy, dscal, and dcopy are basic linear algebra modules (blas).
!  dumach    sets the unit roundoff of the machine.
!  xerrwd, xsetun, xsetf, ixsav, and iumach handle the printing of all
!            error messages and warnings.  xerrwd is machine-dependent.
! note:  dvnorm, dumach, ixsav, and iumach are function routines.
! all the others are subroutines.
!
!-----------------------------------------------------------------------

subroutine dvode(f,neq,y,t,tout,itol,rtol,atol,itask,istate,iopt, &
                 rwork,lrw,iwork,liw,jac,mf,rpar,ipar)
implicit none
external f , jac
double precision y , t , tout , rtol , atol , rwork , rpar
integer neq , itol , itask , istate , iopt , lrw , iwork , liw ,  &
        mf , ipar
dimension y(*) , rtol(*) , atol(*) , rwork(lrw) , iwork(liw) ,    &
          rpar(*) , ipar(*)

!
! type declarations for labeled common block dvod01 --------------------
!
      double precision acnrm , ccmxj , conp , crate , drc , el , eta ,  &
                       etamax , h , hmin , hmxi , hnew , hscal , prl1 , &
                       rc , rl1 , tau , tq , tn , uround
      integer icf , init , ipup , jcur , jstart , jsv , kflag , kuth ,  &
              l , lmax , lyh , lewt , lacor , lsavf , lwm , liwm ,      &
              locjs , maxord , meth , miter , msbj , mxhnil , mxstep ,  &
              n , newh , newq , nhnil , nq , nqnyh , nqwait , nslj ,    &
              nslp , nyh
!
! type declarations for labeled common block dvod02 --------------------
!
      double precision hu
      integer ncfn , netf , nfe , nje , nlu , nni , nqu , nst
!
! type declarations for local variables --------------------------------
!
!      external dvnlsd
      logical ihit
      double precision atoli , big , ewti , four , h0 , hmax , hmx ,    &
                       hun , one , pt2 , rh , rtoli , size , tcrit ,    &
                       tnext , tolsf , tp , two , zero
      integer i , ier , iflag , imxer , jco , kgo , leniw , lenj ,      &
              lenp , lenrw , lenwm , lf0 , mband , mfa , ml , mord ,    &
              mu , mxhnl0 , mxstp0 , niter , nslast
      character(len=80) msg
!
! type declaration for function subroutines called ---------------------
!
!      double precision dumach , dvnorm
!
      dimension mord(2)
!-----------------------------------------------------------------------
! the following fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to dvode.
!-----------------------------------------------------------------------
      save mord , mxhnl0 , mxstp0
      save zero , one , two , four , pt2 , hun
!-----------------------------------------------------------------------
! the following internal common blocks contain variables which are
! communicated between subroutines in the dvode package, or which are
! to be saved between calls to dvode.
! in each block, real variables precede integers.
! the block /dvod01/ appears in subroutines dvode, dvindy, dvstep,
! dvset, dvnlsd, dvjac, dvsol, dvjust and dvsrco.
! the block /dvod02/ appears in subroutines dvode, dvindy, dvstep,
! dvnlsd, dvjac, and dvsrco.
!
! the variables stored in the internal common blocks are as follows:
!
! acnrm  = weighted r.m.s. norm of accumulated correction vectors.
! ccmxj  = threshhold on drc for updating the jacobian. (see drc.)
! conp   = the saved value of tq(5).
! crate  = estimated corrector convergence rate constant.
! drc    = relative change in h*rl1 since last dvjac call.
! el     = real array of integration coefficients.  see dvset.
! eta    = saved tentative ratio of new to old h.
! etamax = saved maximum value of eta to be allowed.
! h      = the step size.
! hmin   = the minimum absolute value of the step size h to be used.
! hmxi   = inverse of the maximum absolute value of h to be used.
!          hmxi = 0.0 is allowed and corresponds to an infinite hmax.
! hnew   = the step size to be attempted on the next step.
! hscal  = stepsize in scaling of yh array.
! prl1   = the saved value of rl1.
! rc     = ratio of current h*rl1 to value on last dvjac call.
! rl1    = the reciprocal of the coefficient el(1).
! tau    = real vector of past nq step sizes, length 13.
! tq     = a real vector of length 5 in which dvset stores constants
!          used for the convergence test, the error test, and the
!          selection of h at a new order.
! tn     = the independent variable, updated on each step taken.
! uround = the machine unit roundoff.  the smallest positive real number
!          such that  1.0 + uround .ne. 1.0
! icf    = integer flag for convergence failure in dvnlsd:
!            0 means no failures.
!            1 means convergence failure with out of date jacobian
!                   (recoverable error).
!            2 means convergence failure with current jacobian or
!                   singular matrix (unrecoverable error).
! init   = saved integer flag indicating whether initialization of the
!          problem has been done (init = 1) or not.
! ipup   = saved flag to signal updating of newton matrix.
! jcur   = output flag from dvjac showing jacobian status:
!            jcur = 0 means j is not current.
!            jcur = 1 means j is current.
! jstart = integer flag used as input to dvstep:
!            0  means perform the first step.
!            1  means take a new step continuing from the last.
!            -1 means take the next step with a new value of maxord,
!                  hmin, hmxi, n, meth, miter, and/or matrix parameters.
!          on return, dvstep sets jstart = 1.
! jsv    = integer flag for jacobian saving, = sign(mf).
! kflag  = a completion code from dvstep with the following meanings:
!               0      the step was succesful.
!              -1      the requested error could not be achieved.
!              -2      corrector convergence could not be achieved.
!              -3, -4  fatal error in vnls (can not occur here).
! kuth   = input flag to dvstep showing whether h was reduced by the
!          driver.  kuth = 1 if h was reduced, = 0 otherwise.
! l      = integer variable, nq + 1, current order plus one.
! lmax   = maxord + 1 (used for dimensioning).
! locjs  = a pointer to the saved jacobian, whose storage starts at
!          wm(locjs), if jsv = 1.
! lyh, lewt, lacor, lsavf, lwm, liwm = saved integer pointers
!          to segments of rwork and iwork.
! maxord = the maximum order of integration method to be allowed.
! meth/miter = the method flags.  see mf.
! msbj   = the maximum number of steps between j evaluations, = 50.
! mxhnil = saved value of optional input mxhnil.
! mxstep = saved value of optional input mxstep.
! n      = the number of first-order odes, = neq.
! newh   = saved integer to flag change of h.
! newq   = the method order to be used on the next step.
! nhnil  = saved counter for occurrences of t + h = t.
! nq     = integer variable, the current integration method order.
! nqnyh  = saved value of nq*nyh.
! nqwait = a counter controlling the frequency of order changes.
!          an order change is about to be considered if nqwait = 1.
! nslj   = the number of steps taken as of the last jacobian update.
! nslp   = saved value of nst as of last newton matrix update.
! nyh    = saved value of the initial value of neq.
! hu     = the step size in t last used.
! ncfn   = number of nonlinear convergence failures so far.
! netf   = the number of error test failures of the integrator so far.
! nfe    = the number of f evaluations for the problem so far.
! nje    = the number of jacobian evaluations so far.
! nlu    = the number of matrix lu decompositions so far.
! nni    = number of nonlinear iterations so far.
! nqu    = the method order last used.
! nst    = the number of steps taken for the problem so far.
!-----------------------------------------------------------------------
      common /dvod01/ acnrm , ccmxj , conp , crate , drc , el(13) ,     &
                      eta , etamax , h , hmin , hmxi , hnew , hscal ,   &
                      prl1 , rc , rl1 , tau(13) , tq(5) , tn , uround , &
                      icf , init , ipup , jcur , jstart , jsv , kflag , &
                      kuth , l , lmax , lyh , lewt , lacor , lsavf ,    &
                      lwm , liwm , locjs , maxord , meth , miter ,      &
                      msbj , mxhnil , mxstep , n , newh , newq , nhnil ,&
                      nq , nqnyh , nqwait , nslj , nslp , nyh
      common /dvod02/ hu , ncfn , netf , nfe , nje , nlu , nni , nqu ,  &
                      nst
!
      data mord(1)/12/ , mord(2)/5/ , mxstp0/500/ , mxhnl0/10/
      data zero/0.0d0/ , one/1.0d0/ , two/2.0d0/ , four/4.0d0/ ,        &
           pt2/0.2d0/ , hun/100.0d0/
!-----------------------------------------------------------------------
! block a.
! this code block is executed on every call.
! it tests istate and itask for legality and branches appropriately.
! if istate .gt. 1 but the flag init shows that initialization has
! not yet been done, an error return occurs.
! if istate = 1 and tout = t, return immediately.
!-----------------------------------------------------------------------
      if ( istate<1 .or. istate>3 ) then
!-----------------------------------------------------------------------
! block i.
! the following block handles all error returns due to illegal input
! (istate = -3), as detected before calling the core integrator.
! first the error message routine is called.   if the illegal input
! is a negative istate, the run is aborted (apparent infinite loop).
!-----------------------------------------------------------------------
         msg = 'dvode--  istate (=i1) illegal '
         call xerrwd(msg,30,1,1,1,istate,0,0,zero,zero)
         if ( istate>=0 ) goto 1500
!
         msg = 'dvode--  run aborted:  apparent infinite loop     '
         call xerrwd(msg,50,303,2,0,0,0,0,zero,zero)
         goto 99999
      else
         if ( itask<1 .or. itask>5 ) then
            msg = 'dvode--  itask (=i1) illegal  '
            call xerrwd(msg,30,2,1,1,itask,0,0,zero,zero)
            goto 1500
         else
            if ( istate==1 ) then
               init = 0
               if ( tout==t ) return
            elseif ( init/=1 ) then
               msg =                                                    &
          'dvode--  istate (=i1) .gt. 1 but dvode not initialized      '
               call xerrwd(msg,60,3,1,1,istate,0,0,zero,zero)
               goto 1500
            elseif ( istate==2 ) then
               goto 50
            endif
!-----------------------------------------------------------------------
! block b.
! the next code block is executed for the initial call (istate = 1),
! or for a continuation call with parameter changes (istate = 3).
! it contains checking of all input and various initializations.
!
! first check legality of the non-optional input neq, itol, iopt,
! mf, ml, and mu.
!-----------------------------------------------------------------------
            if ( neq<=0 ) then
               msg = 'dvode--  neq (=i1) .lt. 1     '
               call xerrwd(msg,30,4,1,1,neq,0,0,zero,zero)
               goto 1500
            else
               if ( istate/=1 ) then
                  if ( neq>n ) then
                     msg =                                              &
                    'dvode--  istate = 3 and neq increased (i1 to i2)  '
                     call xerrwd(msg,50,5,1,2,n,neq,0,zero,zero)
                     goto 1500
                  endif
               endif
               n = neq
               if ( itol<1 .or. itol>4 ) then
                  msg = 'dvode--  itol (=i1) illegal   '
                  call xerrwd(msg,30,6,1,1,itol,0,0,zero,zero)
                  goto 1500
               elseif ( iopt<0 .or. iopt>1 ) then
                  msg = 'dvode--  iopt (=i1) illegal   '
                  call xerrwd(msg,30,7,1,1,iopt,0,0,zero,zero)
                  goto 1500
               else
                  jsv = sign(1,mf)
                  mfa = abs(mf)
                  meth = mfa/10
                  miter = mfa - 10*meth
                  if ( meth<1 .or. meth>2 ) goto 800
                  if ( miter<0 .or. miter>5 ) goto 800
                  if ( miter>3 ) then
                     ml = iwork(1)
                     mu = iwork(2)
                     if ( ml<0 .or. ml>=n ) then
                        msg =                                           &
                    'dvode--  ml (=i1) illegal:  .lt.0 or .ge.neq (=i2)'
                        call xerrwd(msg,50,9,1,2,ml,neq,0,zero,zero)
                        goto 1500
                     elseif ( mu<0 .or. mu>=n ) then
                        msg =                                           &
                    'dvode--  mu (=i1) illegal:  .lt.0 or .ge.neq (=i2)'
                        call xerrwd(msg,50,10,1,2,mu,neq,0,zero,zero)
                        goto 1500
                     endif
                  endif
! next process and check the optional input. ---------------------------
                  if ( iopt==1 ) then
                     maxord = iwork(5)
                     if ( maxord<0 ) then
                        msg = 'dvode--  maxord (=i1) .lt. 0  '
                        call xerrwd(msg,30,11,1,1,maxord,0,0,zero,zero)
                        goto 1500
                     else
                        if ( maxord==0 ) maxord = 100
                        maxord = min(maxord,mord(meth))
                        mxstep = iwork(6)
                        if ( mxstep<0 ) then
                           msg = 'dvode--  mxstep (=i1) .lt. 0  '
                           call xerrwd(msg,30,12,1,1,mxstep,0,0,zero,   &
                                       zero)
                           goto 1500
                        else
                           if ( mxstep==0 ) mxstep = mxstp0
                           mxhnil = iwork(7)
                           if ( mxhnil<0 ) then
                              msg = 'dvode--  mxhnil (=i1) .lt. 0  '
                              call xerrwd(msg,30,13,1,1,mxhnil,0,0,zero,&
                                 zero)
                              goto 1500
                           else
                              if ( mxhnil==0 ) mxhnil = mxhnl0
                              if ( istate==1 ) then
                                 h0 = rwork(5)
                                 if ( (tout-t)*h0<zero ) then
                                    msg =                               &
                              'dvode--  tout (=r1) behind t (=r2)      '
                                    call xerrwd(msg,40,14,1,0,0,0,2,    &
                                       tout,t)
                                    msg =                               &
                    '      integration direction is given by h0 (=r1)  '
                                    call xerrwd(msg,50,14,1,0,0,0,1,h0, &
                                       zero)
                                    goto 1500
                                 endif
                              endif
                              hmax = rwork(6)
                              if ( hmax<zero ) then
                                 msg = 'dvode--  hmax (=r1) .lt. 0.0  '
                                 call xerrwd(msg,30,15,1,0,0,0,1,hmax,  &
                                    zero)
                                 goto 1500
                              else
                                 hmxi = zero
                                 if ( hmax>zero ) hmxi = one/hmax
                                 hmin = rwork(7)
                                 if ( hmin<zero ) then
                                    msg =                               &
                                       'dvode--  hmin (=r1) .lt. 0.0  '
                                    call xerrwd(msg,30,16,1,0,0,0,1,    &
                                       hmin,zero)
                                    goto 1500
                                 endif
                              endif
                           endif
                        endif
                     endif
                  else
                     maxord = mord(meth)
                     mxstep = mxstp0
                     mxhnil = mxhnl0
                     if ( istate==1 ) h0 = zero
                     hmxi = zero
                     hmin = zero
                  endif
!-----------------------------------------------------------------------
! set work array pointers and check lengths lrw and liw.
! pointers to segments of rwork and iwork are named by prefixing l to
! the name of the segment.  e.g., the segment yh starts at rwork(lyh).
! segments of rwork (in order) are denoted  yh, wm, ewt, savf, acor.
! within wm, locjs is the location of the saved jacobian (jsv .gt. 0).
!-----------------------------------------------------------------------
                  lyh = 21
                  if ( istate==1 ) nyh = n
                  lwm = lyh + (maxord+1)*nyh
                  jco = max(0,jsv)
                  if ( miter==0 ) lenwm = 0
                  if ( miter==1 .or. miter==2 ) then
                     lenwm = 2 + (1+jco)*n*n
                     locjs = n*n + 3
                  endif
                  if ( miter==3 ) lenwm = 2 + n
                  if ( miter==4 .or. miter==5 ) then
                     mband = ml + mu + 1
                     lenp = (mband+ml)*n
                     lenj = mband*n
                     lenwm = 2 + lenp + jco*lenj
                     locjs = lenp + 3
                  endif
                  lewt = lwm + lenwm
                  lsavf = lewt + n
                  lacor = lsavf + n
                  lenrw = lacor + n - 1
                  iwork(17) = lenrw
                  liwm = 1
                  leniw = 30 + n
                  if ( miter==0 .or. miter==3 ) leniw = 30
                  iwork(18) = leniw
                  if ( lenrw>lrw ) then
                     msg =                                              &
          'dvode--  rwork length needed, lenrw (=i1), exceeds lrw (=i2)'
                     call xerrwd(msg,60,17,1,2,lenrw,lrw,0,zero,zero)
                     goto 1500
                  elseif ( leniw>liw ) then
                     msg =                                              &
          'dvode--  iwork length needed, leniw (=i1), exceeds liw (=i2)'
                     call xerrwd(msg,60,18,1,2,leniw,liw,0,zero,zero)
                     goto 1500
                  else
! check rtol and atol for legality. ------------------------------------
                     rtoli = rtol(1)
                     atoli = atol(1)
                     do i = 1 , n
                        if ( itol>=3 ) rtoli = rtol(i)
                        if ( itol==2 .or. itol==4 ) atoli = atol(i)
                        if ( rtoli<zero ) goto 900
                        if ( atoli<zero ) goto 1000
                     enddo
                     if ( istate==1 ) then
!-----------------------------------------------------------------------
! block c.
! the next block is for the initial call only (istate = 1).
! it contains all remaining initializations, the initial call to f,
! and the calculation of the initial step size.
! the error weights in ewt are inverted after being loaded.
!-----------------------------------------------------------------------
                        uround = epmach
                        tn = t
                        if ( itask==4 .or. itask==5 ) then
                           tcrit = rwork(1)
                           if ( (tcrit-tout)*(tout-t)<zero ) goto 1300
                           if ( h0/=zero .and. (t+h0-tcrit)*h0>zero )   &
                                h0 = tcrit - t
                        endif
                        jstart = 0
                        if ( miter>0 ) rwork(lwm) = sqrt(uround)
                        ccmxj = pt2
                        msbj = 50
                        nhnil = 0
                        nst = 0
                        nje = 0
                        nni = 0
                        ncfn = 0
                        netf = 0
                        nlu = 0
                        nslj = 0
                        nslast = 0
                        hu = zero
                        nqu = 0
! initial call to f.  (lf0 points to yh(*,2).) -------------------------
                        lf0 = lyh + nyh
                        call f(n,t,y,rwork(lf0),rpar,ipar)
                        nfe = 1
! load the initial value vector in yh. ---------------------------------
                        call dcopy(n,y,1,rwork(lyh),1)
! load and invert the ewt array.  (h is temporarily set to 1.0.) -------
                        nq = 1
                        h = one
                        call dewset(n,itol,rtol,atol,rwork(lyh),        &
                                    rwork(lewt))
                        do i = 1 , n
                           if ( rwork(i+lewt-1)<=zero ) goto 1100
                           rwork(i+lewt-1) = one/rwork(i+lewt-1)
                        enddo
                        if ( h0==zero ) then
! call dvhin to set initial step size h0 to be attempted. --------------
                           call dvhin(n,t,rwork(lyh),rwork(lf0),f,rpar, &
                                      ipar,tout,uround,rwork(lewt),itol,&
                                      atol,y,rwork(lacor),h0,niter,ier)
                           nfe = nfe + niter
                           if ( ier/=0 ) then
                              msg =                                     &
          'dvode--  tout (=r1) too close to t(=r2) to start integration'
                              call xerrwd(msg,60,22,1,0,0,0,2,tout,t)
                              goto 1500
                           endif
                        endif
! adjust h0 if necessary to meet hmax bound. ---------------------------
                        rh = abs(h0)*hmxi
                        if ( rh>one ) h0 = h0/rh
! load h with h0 and scale yh(*,2) by h0. ------------------------------
                        h = h0
                        call dscal(n,h0,rwork(lf0),1)
                        goto 200
                     else
! if istate = 3, set flag to signal parameter changes to dvstep. -------
                        jstart = -1
! maxord was reduced below nq.  copy yh(*,maxord+2) into savf. ---------
                        if ( nq>maxord )                                &
                             call dcopy(n,rwork(lwm),1,rwork(lsavf),1)
! reload wm(1) = rwork(lwm), since lwm may have changed. ---------------
                        if ( miter>0 ) rwork(lwm) = sqrt(uround)
                     endif
                  endif
               endif
            endif
         endif
!-----------------------------------------------------------------------
! block d.
! the next code block is for continuation calls only (istate = 2 or 3)
! and is to check stop conditions before taking a step.
!-----------------------------------------------------------------------
 50      nslast = nst
         kuth = 0
         select case (itask)
         case (2)
            goto 100
         case (3)
            tp = tn - hu*(one+hun*uround)
            if ( (tp-tout)*h>zero ) then
               msg =                                                    &
          'dvode--  itask = i1 and tout (=r1) behind tcur - hu (= r2)  '
               call xerrwd(msg,60,23,1,1,itask,0,2,tout,tp)
               goto 1500
            else
               if ( (tn-tout)*h>=zero ) goto 300
               goto 100
            endif
         case (4)
            tcrit = rwork(1)
            if ( (tn-tcrit)*h>zero ) goto 1200
            if ( (tcrit-tout)*h<zero ) goto 1300
            if ( (tn-tout)*h>=zero ) then
               call dvindy(tout,0,rwork(lyh),nyh,y,iflag)
               if ( iflag/=0 ) goto 1400
               t = tout
               goto 400
            endif
         case (5)
            tcrit = rwork(1)
            if ( (tn-tcrit)*h>zero ) goto 1200
         case default
            if ( (tn-tout)*h<zero ) goto 100
            call dvindy(tout,0,rwork(lyh),nyh,y,iflag)
            if ( iflag/=0 ) goto 1400
            t = tout
            goto 400
         end select
         hmx = abs(tn) + abs(h)
         ihit = abs(tn-tcrit)<=hun*uround*hmx
         if ( ihit ) goto 300
         tnext = tn + hnew*(one+four*uround)
         if ( (tnext-tcrit)*h>zero ) then
            h = (tcrit-tn)*(one-four*uround)
            kuth = 1
         endif
      endif
!-----------------------------------------------------------------------
! block e.
! the next block is normally executed for all calls and contains
! the call to the one-step core integrator dvstep.
!
! this is a looping point for the integration steps.
!
! first check for too many steps being taken, update ewt (if not at
! start of problem), check for too much accuracy being requested, and
! check for h below the roundoff level in t.
!-----------------------------------------------------------------------
 100  if ( (nst-nslast)>=mxstep ) then
!-----------------------------------------------------------------------
! block h.
! the following block handles all unsuccessful returns other than
! those for illegal input.  first the error message routine is called.
! if there was an error test or convergence test failure, imxer is set.
! then y is loaded from yh, and t is set to tn.
! the optional output is loaded into the work arrays before returning.
!-----------------------------------------------------------------------
! the maximum number of steps was taken before reaching tout. ----------
         msg = 'dvode--  at current t (=r1), mxstep (=i1) steps   '
         call xerrwd(msg,50,201,1,0,0,0,0,zero,zero)
         msg = '      taken on this call before reaching tout     '
         call xerrwd(msg,50,201,1,1,mxstep,0,1,tn,zero)
         istate = -1
         goto 700
      else
         call dewset(n,itol,rtol,atol,rwork(lyh),rwork(lewt))
         do i = 1 , n
            if ( rwork(i+lewt-1)<=zero ) goto 500
            rwork(i+lewt-1) = one/rwork(i+lewt-1)
         enddo
      endif
 200  tolsf = uround*dvnorm(n,rwork(lyh),rwork(lewt))
      if ( tolsf<=one ) then
         if ( (tn+h)==tn ) then
            nhnil = nhnil + 1
            if ( nhnil<=mxhnil ) then
               msg =                                                    &
                    'dvode--  warning: internal t (=r1) and h (=r2) are'
               call xerrwd(msg,50,101,1,0,0,0,0,zero,zero)
               msg =                                                    &
          '      such that in the machine, t + h = t on the next step  '
               call xerrwd(msg,60,101,1,0,0,0,0,zero,zero)
               msg =                                                    &
                    '      (h = step size). solver will continue anyway'
               call xerrwd(msg,50,101,1,0,0,0,2,tn,h)
               if ( nhnil>=mxhnil ) then
                  msg =                                                 &
                    'dvode--  above warning has been issued i1 times.  '
                  call xerrwd(msg,50,102,1,0,0,0,0,zero,zero)
                  msg =                                                 &
                    '      it will not be issued again for this problem'
                  call xerrwd(msg,50,102,1,1,mxhnil,0,0,zero,zero)
               endif
            endif
         endif
!-----------------------------------------------------------------------
! call dvstep (y, yh, nyh, yh, ewt, savf, vsav, acor,
!              wm, iwm, f, jac, f, dvnlsd, rpar, ipar)
!-----------------------------------------------------------------------
         call dvstep(y,rwork(lyh),nyh,rwork(lyh),rwork(lewt),           &
                     rwork(lsavf),y,rwork(lacor),rwork(lwm),iwork(liwm),&
                     f,jac,f,dvnlsd,rpar,ipar)
         kgo = 1 - kflag
! branch on kflag.  note: in this version, kflag can not be set to -3.
!  kflag .eq. 0,   -1,  -2
         select case (kgo)
         case (2)
! kflag = -1.  error test failed repeatedly or with abs(h) = hmin. -----
            msg = 'dvode--  at t(=r1) and step size h(=r2), the error'
            call xerrwd(msg,50,204,1,0,0,0,0,zero,zero)
            msg = '      test failed repeatedly or with abs(h) = hmin'
            call xerrwd(msg,50,204,1,0,0,0,2,tn,h)
            istate = -4
            goto 600
         case (3)
! kflag = -2.  convergence failed repeatedly or with abs(h) = hmin. ----
            msg = 'dvode--  at t (=r1) and step size h (=r2), the    '
            call xerrwd(msg,50,205,1,0,0,0,0,zero,zero)
            msg = '      corrector convergence failed repeatedly     '
            call xerrwd(msg,50,205,1,0,0,0,0,zero,zero)
            msg = '      or with abs(h) = hmin   '
            call xerrwd(msg,30,205,1,0,0,0,2,tn,h)
            istate = -5
            goto 600
         case default
!-----------------------------------------------------------------------
! block f.
! the following block handles the case of a successful return from the
! core integrator (kflag = 0).  test for stop conditions.
!-----------------------------------------------------------------------
            init = 1
            kuth = 0
            select case (itask)
            case (2)
            case (3)
! itask = 3.  jump to exit if tout was reached. ------------------------
               if ( (tn-tout)*h<zero ) goto 100
            case (4)
! itask = 4.  see if tout or tcrit was reached.  adjust h if necessary.
               if ( (tn-tout)*h<zero ) then
                  hmx = abs(tn) + abs(h)
                  ihit = abs(tn-tcrit)<=hun*uround*hmx
                  if ( .not.(ihit) ) then
                     tnext = tn + hnew*(one+four*uround)
                     if ( (tnext-tcrit)*h>zero ) then
                        h = (tcrit-tn)*(one-four*uround)
                        kuth = 1
                     endif
                     goto 100
                  endif
               else
                  call dvindy(tout,0,rwork(lyh),nyh,y,iflag)
                  t = tout
                  goto 400
               endif
            case (5)
! itask = 5.  see if tcrit was reached and jump to exit. ---------------
               hmx = abs(tn) + abs(h)
               ihit = abs(tn-tcrit)<=hun*uround*hmx
            case default
! itask = 1.  if tout has been reached, interpolate. -------------------
               if ( (tn-tout)*h<zero ) goto 100
               call dvindy(tout,0,rwork(lyh),nyh,y,iflag)
               t = tout
               goto 400
            end select
         end select
      else
         tolsf = tolsf*two
         if ( nst==0 ) then
            msg = 'dvode--  at start of problem, too much accuracy   '
            call xerrwd(msg,50,26,1,0,0,0,0,zero,zero)
            msg =                                                       &
          '      requested for precision of machine:   see tolsf (=r1) '
            call xerrwd(msg,60,26,1,0,0,0,1,tolsf,zero)
            rwork(14) = tolsf
            goto 1500
         else
! too much accuracy requested for machine precision. -------------------
            msg = 'dvode--  at t (=r1), too much accuracy requested  '
            call xerrwd(msg,50,203,1,0,0,0,0,zero,zero)
            msg = '      for precision of machine:   see tolsf (=r2) '
            call xerrwd(msg,50,203,1,0,0,0,2,tn,tolsf)
            rwork(14) = tolsf
            istate = -2
            goto 700
         endif
      endif
!-----------------------------------------------------------------------
! block g.
! the following block handles all successful returns from dvode.
! if itask .ne. 1, y is loaded from yh and t is set accordingly.
! istate is set to 2, and the optional output is loaded into the work
! arrays before returning.
!-----------------------------------------------------------------------
 300  call dcopy(n,rwork(lyh),1,y,1)
      t = tn
      if ( itask==4 .or. itask==5 ) then
         if ( ihit ) t = tcrit
      endif
 400  istate = 2
      rwork(11) = hu
      rwork(12) = hnew
      rwork(13) = tn
      iwork(11) = nst
      iwork(12) = nfe
      iwork(13) = nje
      iwork(14) = nqu
      iwork(15) = newq
      iwork(19) = nlu
      iwork(20) = nni
      iwork(21) = ncfn
      iwork(22) = netf
      return
! ewt(i) .le. 0.0 for some i (not at start of problem). ----------------
 500  ewti = rwork(lewt+i-1)
      msg = 'dvode--  at t (=r1), ewt(i1) has become r2 .le. 0.'
      call xerrwd(msg,50,202,1,1,i,0,2,tn,ewti)
      istate = -6
      goto 700
! compute imxer if relevant. -------------------------------------------
 600  big = zero
      imxer = 1
      do i = 1 , n
         size = abs(rwork(i+lacor-1)*rwork(i+lewt-1))
         if ( big<size ) then
            big = size
            imxer = i
         endif
      enddo
      iwork(16) = imxer
! set y vector, t, and optional output. --------------------------------
 700  call dcopy(n,rwork(lyh),1,y,1)
      t = tn
      rwork(11) = hu
      rwork(12) = h
      rwork(13) = tn
      iwork(11) = nst
      iwork(12) = nfe
      iwork(13) = nje
      iwork(14) = nqu
      iwork(15) = nq
      iwork(19) = nlu
      iwork(20) = nni
      iwork(21) = ncfn
      iwork(22) = netf
      return
 800  msg = 'dvode--  mf (=i1) illegal     '
      call xerrwd(msg,30,8,1,1,mf,0,0,zero,zero)
      goto 1500
 900  msg = 'dvode--  rtol(i1) is r1 .lt. 0.0        '
      call xerrwd(msg,40,19,1,1,i,0,1,rtoli,zero)
      goto 1500
 1000 msg = 'dvode--  atol(i1) is r1 .lt. 0.0        '
      call xerrwd(msg,40,20,1,1,i,0,1,atoli,zero)
      goto 1500
 1100 ewti = rwork(lewt+i-1)
      msg = 'dvode--  ewt(i1) is r1 .le. 0.0         '
      call xerrwd(msg,40,21,1,1,i,0,1,ewti,zero)
      goto 1500
 1200 msg =                                                             &
          'dvode--  itask = 4 or 5 and tcrit (=r1) behind tcur (=r2)   '
      call xerrwd(msg,60,24,1,0,0,0,2,tcrit,tn)
      goto 1500
 1300 msg =                                                             &
          'dvode--  itask = 4 or 5 and tcrit (=r1) behind tout (=r2)   '
      call xerrwd(msg,60,25,1,0,0,0,2,tcrit,tout)
      goto 1500
 1400 msg =                                                             &
          'dvode--  trouble from dvindy.  itask = i1, tout = r1.       '
      call xerrwd(msg,60,27,1,1,itask,0,1,tout,zero)
!
 1500 istate = -3
      return

99999 end subroutine dvode

      subroutine dvhin(n,t0,y0,ydot,f,rpar,ipar,tout,uround,ewt,itol,   &
                       atol,y,temp,h0,niter,ier)
      implicit none
      external f
      double precision t0 , y0 , ydot , rpar , tout , uround , ewt ,    &
                       atol , y , temp , h0
      integer n , ipar , itol , niter , ier
      dimension y0(*) , ydot(*) , ewt(*) , atol(*) , y(*) , temp(*) ,   &
                rpar(*) , ipar(*)
!-----------------------------------------------------------------------
! call sequence input -- n, t0, y0, ydot, f, rpar, ipar, tout, uround,
!                        ewt, itol, atol, y, temp
! call sequence output -- h0, niter, ier
! common block variables accessed -- none
!
! subroutines called by dvhin:  f
! function routines called by dvhi: dvnorm
!-----------------------------------------------------------------------
! this routine computes the step size, h0, to be attempted on the
! first step, when the user has not supplied a value for this.
!
! first we check that tout - t0 differs significantly from zero.  then
! an iteration is done to approximate the initial second derivative
! and this is used to define h from w.r.m.s.norm(h**2 * yddot / 2) = 1.
! a bias factor of 1/2 is applied to the resulting h.
! the sign of h0 is inferred from the initial values of tout and t0.
!
! communication with dvhin is done with the following variables:
!
! n      = size of ode system, input.
! t0     = initial value of independent variable, input.
! y0     = vector of initial conditions, input.
! ydot   = vector of initial first derivatives, input.
! f      = name of subroutine for right-hand side f(t,y), input.
! rpar, ipar = dummy names for user's real and integer work arrays.
! tout   = first output value of independent variable
! uround = machine unit roundoff
! ewt, itol, atol = error weights and tolerance parameters
!                   as described in the driver routine, input.
! y, temp = work arrays of length n.
! h0     = step size to be attempted, output.
! niter  = number of iterations (and of f evaluations) to compute h0,
!          output.
! ier    = the error flag, returned with the value
!          ier = 0  if no trouble occurred, or
!          ier = -1 if tout and t0 are considered too close to proceed.
!-----------------------------------------------------------------------
!
! type declarations for local variables --------------------------------
!
      double precision afi , atoli , delyi , h , half , hg , hlb ,      &
                       hnew , hrat , hub , hun , pt1 , t1 , tdist ,     &
                       tround , two , yddnrm
      integer i , iter
!
! type declaration for function subroutines called ---------------------
!
!      double precision dvnorm
!-----------------------------------------------------------------------
! the following fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
      save half , hun , pt1 , two
      data half/0.5d0/ , hun/100.0d0/ , pt1/0.1d0/ , two/2.0d0/
!
      niter = 0
      tdist = abs(tout-t0)
      tround = uround*max(abs(t0),abs(tout))
      if ( tdist<two*tround ) then
! error return for tout - t0 too small. --------------------------------
         ier = -1
         goto 99999
      else
!
! set a lower bound on h based on the roundoff level in t0 and tout. ---
         hlb = hun*tround
! set an upper bound on h based on tout-t0 and the initial y and ydot. -
         hub = pt1*tdist
         atoli = atol(1)
         do i = 1 , n
            if ( itol==2 .or. itol==4 ) atoli = atol(i)
            delyi = pt1*abs(y0(i)) + atoli
            afi = abs(ydot(i))
            if ( afi*hub>delyi ) hub = delyi/afi
         enddo
!
! set initial guess for h as geometric mean of upper and lower bounds. -
         iter = 0
         hg = sqrt(hlb*hub)
! if the bounds have crossed, exit with the mean value. ----------------
         if ( hub<hlb ) then
            h0 = hg
            goto 200
         endif
!
! looping point for iteration. -----------------------------------------
! estimate the second derivative as a difference quotient in f. --------
 50      h = sign(hg,tout-t0)
         t1 = t0 + h
         do i = 1 , n
            y(i) = y0(i) + h*ydot(i)
         enddo
         call f(n,t1,y,temp,rpar,ipar)
         do i = 1 , n
            temp(i) = (temp(i)-ydot(i))/h
         enddo
         yddnrm = dvnorm(n,temp,ewt)
! get the corresponding new value of h. --------------------------------
         if ( yddnrm*hub*hub>two ) then
            hnew = sqrt(two/yddnrm)
         else
            hnew = sqrt(hg*hub)
         endif
         iter = iter + 1
!-----------------------------------------------------------------------
! test the stopping conditions.
! stop if the new and previous h values differ by a factor of .lt. 2.
! stop if four iterations have been done.  also, stop with previous h
! if hnew/hg .gt. 2 after first iteration, as this probably means that
! the second derivative value is bad because of cancellation error.
!-----------------------------------------------------------------------
         if ( iter<4 ) then
            hrat = hnew/hg
            if ( (hrat<=half) .or. (hrat>=two) ) then
               if ( (iter>=2) .and. (hnew>two*hg) ) then
                  hnew = hg
                  goto 100
               endif
               hg = hnew
               goto 50
            endif
         endif
!
! iteration done.  apply bounds, bias factor, and sign.  then exit. ----
 100     h0 = hnew*half
         if ( h0<hlb ) h0 = hlb
         if ( h0>hub ) h0 = hub
      endif
 200  h0 = sign(h0,tout-t0)
      niter = iter
      ier = 0
      return

99999 end subroutine dvhin

      subroutine dvindy(t,k,yh,ldyh,dky,iflag)
      implicit none
      double precision t , yh , dky
      integer k , ldyh , iflag
      dimension yh(ldyh,*) , dky(*)
!-----------------------------------------------------------------------
! call sequence input -- t, k, yh, ldyh
! call sequence output -- dky, iflag
! common block variables accessed:
!     /dvod01/ --  h, tn, uround, l, n, nq
!     /dvod02/ --  hu
!
! subroutines called by dvindy: dscal, xerrwd
! function routines called by dvindy: none
!-----------------------------------------------------------------------
! dvindy computes interpolated values of the k-th derivative of the
! dependent variable vector y, and stores it in dky.  this routine
! is called within the package with k = 0 and t = tout, but may
! also be called by the user for any k up to the current order.
! (see detailed instructions in the usage documentation.)
!-----------------------------------------------------------------------
! the computed values in dky are gotten by interpolation using the
! nordsieck history array yh.  this array corresponds uniquely to a
! vector-valued polynomial of degree nqcur or less, and dky is set
! to the k-th derivative of this polynomial at t.
! the formula for dky is:
!              q
!  dky(i)  =  sum  c(j,k) * (t - tn)**(j-k) * h**(-j) * yh(i,j+1)
!             j=k
! where  c(j,k) = j*(j-1)*...*(j-k+1), q = nqcur, tn = tcur, h = hcur.
! the quantities  nq = nqcur, l = nq+1, n, tn, and h are
! communicated by common.  the above sum is done in reverse order.
! iflag is returned negative if either k or t is out of bounds.
!
! discussion above and comments in driver explain all variables.
!-----------------------------------------------------------------------
!
! type declarations for labeled common block dvod01 --------------------
!
      double precision acnrm , ccmxj , conp , crate , drc , el , eta ,  &
                       etamax , h , hmin , hmxi , hnew , hscal , prl1 , &
                       rc , rl1 , tau , tq , tn , uround
      integer icf , init , ipup , jcur , jstart , jsv , kflag , kuth ,  &
              l , lmax , lyh , lewt , lacor , lsavf , lwm , liwm ,      &
              locjs , maxord , meth , miter , msbj , mxhnil , mxstep ,  &
              n , newh , newq , nhnil , nq , nqnyh , nqwait , nslj ,    &
              nslp , nyh
!
! type declarations for labeled common block dvod02 --------------------
!
      double precision hu
      integer ncfn , netf , nfe , nje , nlu , nni , nqu , nst
!
! type declarations for local variables --------------------------------
!
      double precision c , hun , r , s , tfuzz , tn1 , tp , zero
      integer i , ic , j , jb , jb2 , jj , jj1 , jp1
      character(len=80) msg
!-----------------------------------------------------------------------
! the following fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
      save hun , zero
!
      common /dvod01/ acnrm , ccmxj , conp , crate , drc , el(13) ,     &
                      eta , etamax , h , hmin , hmxi , hnew , hscal ,   &
                      prl1 , rc , rl1 , tau(13) , tq(5) , tn , uround , &
                      icf , init , ipup , jcur , jstart , jsv , kflag , &
                      kuth , l , lmax , lyh , lewt , lacor , lsavf ,    &
                      lwm , liwm , locjs , maxord , meth , miter ,      &
                      msbj , mxhnil , mxstep , n , newh , newq , nhnil ,&
                      nq , nqnyh , nqwait , nslj , nslp , nyh
      common /dvod02/ hu , ncfn , netf , nfe , nje , nlu , nni , nqu ,  &
                      nst
!
      data hun/100.0d0/ , zero/0.0d0/
!
      iflag = 0
      if ( k<0 .or. k>nq ) then
!
         msg = 'dvindy-- k (=i1) illegal      '
         call xerrwd(msg,30,51,1,1,k,0,0,zero,zero)
         iflag = -1
         return
      else
         tfuzz = hun*uround*sign(abs(tn)+abs(hu),hu)
         tp = tn - hu - tfuzz
         tn1 = tn + tfuzz
         if ( (t-tp)*(t-tn1)>zero ) then
            msg = 'dvindy-- t (=r1) illegal      '
            call xerrwd(msg,30,52,1,0,0,0,1,t,zero)
            msg =                                                       &
          '      t not in interval tcur - hu (= r1) to tcur (=r2)      '
            call xerrwd(msg,60,52,1,0,0,0,2,tp,tn)
            iflag = -2
            goto 99999
         else
!
            s = (t-tn)/h
            ic = 1
            if ( k/=0 ) then
               jj1 = l - k
               do jj = jj1 , nq
                  ic = ic*jj
               enddo
            endif
            c = real(ic)
            do i = 1 , n
               dky(i) = c*yh(i,l)
            enddo
            if ( k/=nq ) then
               jb2 = nq - k
               do jb = 1 , jb2
                  j = nq - jb
                  jp1 = j + 1
                  ic = 1
                  if ( k/=0 ) then
                     jj1 = jp1 - k
                     do jj = jj1 , j
                        ic = ic*jj
                     enddo
                  endif
                  c = real(ic)
                  do i = 1 , n
                     dky(i) = c*yh(i,jp1) + s*dky(i)
                  enddo
               enddo
               if ( k==0 ) return
            endif
         endif
      endif
      r = h**(-k)
      call dscal(n,r,dky,1)
      return

99999 end subroutine dvindy

      subroutine dvstep(y,yh,ldyh,yh1,ewt,savf,vsav,acor,wm,iwm,f,jac,  &
                        psol,vnls,rpar,ipar)
      implicit none
      external f , jac , psol , vnls
      double precision y , yh , yh1 , ewt , savf , vsav , acor , wm ,   &
                       rpar
      integer ldyh , iwm , ipar
      dimension y(*) , yh(ldyh,*) , yh1(*) , ewt(*) , savf(*) , vsav(*) &
                , acor(*) , wm(*) , iwm(*) , rpar(*) , ipar(*)
!-----------------------------------------------------------------------
! call sequence input -- y, yh, ldyh, yh1, ewt, savf, vsav,
!                        acor, wm, iwm, f, jac, psol, vnls, rpar, ipar
! call sequence output -- yh, acor, wm, iwm
! common block variables accessed:
!     /dvod01/  acnrm, el(13), h, hmin, hmxi, hnew, hscal, rc, tau(13),
!               tq(5), tn, jcur, jstart, kflag, kuth,
!               l, lmax, maxord, n, newq, nq, nqwait
!     /dvod02/  hu, ncfn, netf, nfe, nqu, nst
!
! subroutines called by dvstep: f, daxpy, dcopy, dscal,
!                               dvjust, vnls, dvset
! function routines called by dvstep: dvnorm
!-----------------------------------------------------------------------
! dvstep performs one step of the integration of an initial value
! problem for a system of ordinary differential equations.
! dvstep calls subroutine vnls for the solution of the nonlinear system
! arising in the time step.  thus it is independent of the problem
! jacobian structure and the type of nonlinear system solution method.
! dvstep returns a completion flag kflag (in common).
! a return with kflag = -1 or -2 means either abs(h) = hmin or 10
! consecutive failures occurred.  on a return with kflag negative,
! the values of tn and the yh array are as of the beginning of the last
! step, and h is the last step size attempted.
!
! communication with dvstep is done with the following variables:
!
! y      = an array of length n used for the dependent variable vector.
! yh     = an ldyh by lmax array containing the dependent variables
!          and their approximate scaled derivatives, where
!          lmax = maxord + 1.  yh(i,j+1) contains the approximate
!          j-th derivative of y(i), scaled by h**j/factorial(j)
!          (j = 0,1,...,nq).  on entry for the first step, the first
!          two columns of yh must be set from the initial values.
! ldyh   = a constant integer .ge. n, the first dimension of yh.
!          n is the number of odes in the system.
! yh1    = a one-dimensional array occupying the same space as yh.
! ewt    = an array of length n containing multiplicative weights
!          for local error measurements.  local errors in y(i) are
!          compared to 1.0/ewt(i) in various error tests.
! savf   = an array of working storage, of length n.
!          also used for input of yh(*,maxord+2) when jstart = -1
!          and maxord .lt. the current order nq.
! vsav   = a work array of length n passed to subroutine vnls.
! acor   = a work array of length n, used for the accumulated
!          corrections.  on a successful return, acor(i) contains
!          the estimated one-step local error in y(i).
! wm,iwm = real and integer work arrays associated with matrix
!          operations in vnls.
! f      = dummy name for the user supplied subroutine for f.
! jac    = dummy name for the user supplied jacobian subroutine.
! psol   = dummy name for the subroutine passed to vnls, for
!          possible use there.
! vnls   = dummy name for the nonlinear system solving subroutine,
!          whose real name is dependent on the method used.
! rpar, ipar = dummy names for user's real and integer work arrays.
!-----------------------------------------------------------------------
!
! type declarations for labeled common block dvod01 --------------------
!
      double precision acnrm , ccmxj , conp , crate , drc , el , eta ,  &
                       etamax , h , hmin , hmxi , hnew , hscal , prl1 , &
                       rc , rl1 , tau , tq , tn , uround
      integer icf , init , ipup , jcur , jstart , jsv , kflag , kuth ,  &
              l , lmax , lyh , lewt , lacor , lsavf , lwm , liwm ,      &
              locjs , maxord , meth , miter , msbj , mxhnil , mxstep ,  &
              n , newh , newq , nhnil , nq , nqnyh , nqwait , nslj ,    &
              nslp , nyh
!
! type declarations for labeled common block dvod02 --------------------
!
      double precision hu
      integer ncfn , netf , nfe , nje , nlu , nni , nqu , nst
!
! type declarations for local variables --------------------------------
!
      double precision addon , bias1 , bias2 , bias3 , cnquot , ddn ,   &
                       dsm , dup , etacf , etamin , etamx1 , etamx2 ,   &
                       etamx3 , etamxf , etaq , etaqm1 , etaqp1 ,       &
                       flotl , one , onepsm , r , thresh , told , zero
      integer i , i1 , i2 , iback , j , jb , kfc , kfh , mxncf , ncf ,  &
              nflag
!
! type declaration for function subroutines called ---------------------
!
!      double precision dvnorm
!-----------------------------------------------------------------------
! the following fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
      save addon , bias1 , bias2 , bias3 , etacf , etamin , etamx1 ,    &
         etamx2 , etamx3 , etamxf , etaq , etaqm1 , kfc , kfh , mxncf , &
         onepsm , thresh , one , zero
!-----------------------------------------------------------------------
      common /dvod01/ acnrm , ccmxj , conp , crate , drc , el(13) ,     &
                      eta , etamax , h , hmin , hmxi , hnew , hscal ,   &
                      prl1 , rc , rl1 , tau(13) , tq(5) , tn , uround , &
                      icf , init , ipup , jcur , jstart , jsv , kflag , &
                      kuth , l , lmax , lyh , lewt , lacor , lsavf ,    &
                      lwm , liwm , locjs , maxord , meth , miter ,      &
                      msbj , mxhnil , mxstep , n , newh , newq , nhnil ,&
                      nq , nqnyh , nqwait , nslj , nslp , nyh
      common /dvod02/ hu , ncfn , netf , nfe , nje , nlu , nni , nqu ,  &
                      nst
!
      data kfc/ - 3/ , kfh/ - 7/ , mxncf/10/
      data addon/1.0d-6/ , bias1/6.0d0/ , bias2/6.0d0/ , bias3/10.0d0/ ,&
           etacf/0.25d0/ , etamin/0.1d0/ , etamxf/0.2d0/ ,              &
           etamx1/1.0d4/ , etamx2/10.0d0/ , etamx3/10.0d0/ ,            &
           onepsm/1.00001d0/ , thresh/1.5d0/
      data one/1.0d0/ , zero/0.0d0/
!
      kflag = 0
      told = tn
      ncf = 0
      jcur = 0
      nflag = 0
      if ( jstart<=0 ) then
         if ( jstart==-1 ) goto 200
!-----------------------------------------------------------------------
! on the first call, the order is set to 1, and other variables are
! initialized.  etamax is the maximum ratio by which h can be increased
! in a single step.  it is normally 10, but is larger during the
! first step to compensate for the small initial h.  if a failure
! occurs (in corrector convergence or error test), etamax is set to 1
! for the next increase.
!-----------------------------------------------------------------------
         lmax = maxord + 1
         nq = 1
         l = 2
         nqnyh = nq*ldyh
         tau(1) = h
         prl1 = one
         rc = zero
         etamax = etamx1
         nqwait = 2
         hscal = h
         goto 400
!-----------------------------------------------------------------------
! take preliminary actions on a normal continuation step (jstart.gt.0).
! if the driver changed h, then eta must be reset and newh set to 1.
! if a change of order was dictated on the previous step, then
! it is done here and appropriate adjustments in the history are made.
! on an order decrease, the history array is adjusted by dvjust.
! on an order increase, the history array is augmented by a column.
! on a change of step size h, the history array yh is rescaled.
!-----------------------------------------------------------------------
      elseif ( kuth==1 ) then
         eta = min(eta,h/hscal)
         newh = 1
      endif
 100  if ( newh==0 ) goto 400
      if ( newq==nq ) goto 300
      if ( newq<nq ) then
         call dvjust(yh,ldyh,-1)
         nq = newq
         l = nq + 1
         nqwait = l
         goto 300
      endif
      if ( newq>nq ) then
         call dvjust(yh,ldyh,1)
         nq = newq
         l = nq + 1
         nqwait = l
         goto 300
      endif
!-----------------------------------------------------------------------
! the following block handles preliminaries needed when jstart = -1.
! if n was reduced, zero out part of yh to avoid undefined references.
! if maxord was reduced to a value less than the tentative order newq,
! then nq is set to maxord, and a new h ratio eta is chosen.
! otherwise, we take the same preliminary actions as for jstart .gt. 0.
! in any case, nqwait is reset to l = nq + 1 to prevent further
! changes in order for that many steps.
! the new h ratio eta is limited by the input h if kuth = 1,
! by hmin if kuth = 0, and by hmxi in any case.
! finally, the history array yh is rescaled.
!-----------------------------------------------------------------------
 200  lmax = maxord + 1
      if ( n/=ldyh ) then
         i1 = 1 + (newq+1)*ldyh
         i2 = (maxord+1)*ldyh
         if ( i1<=i2 ) then
            do i = i1 , i2
               yh1(i) = zero
            enddo
         endif
      endif
      if ( newq>maxord ) then
         flotl = real(lmax)
         if ( maxord<nq-1 ) then
            ddn = dvnorm(n,savf,ewt)/tq(1)
            eta = one/((bias1*ddn)**(one/flotl)+addon)
         endif
         if ( maxord==nq .and. newq==nq+1 ) eta = etaq
         if ( maxord==nq-1 .and. newq==nq+1 ) then
            eta = etaqm1
            call dvjust(yh,ldyh,-1)
         endif
         if ( maxord==nq-1 .and. newq==nq ) then
            ddn = dvnorm(n,savf,ewt)/tq(1)
            eta = one/((bias1*ddn)**(one/flotl)+addon)
            call dvjust(yh,ldyh,-1)
         endif
         eta = min(eta,one)
         nq = maxord
         l = lmax
      endif
      if ( kuth==1 ) eta = min(eta,abs(h/hscal))
      if ( kuth==0 ) eta = max(eta,hmin/abs(hscal))
      eta = eta/max(one,abs(hscal)*hmxi*eta)
      newh = 1
      nqwait = l
      if ( newq<=maxord ) goto 100
! rescale the history array for a change in h by a factor of eta. ------
 300  r = one
      do j = 2 , l
         r = r*eta
         call dscal(n,r,yh(1,j),1)
      enddo
      h = hscal*eta
      hscal = h
      rc = rc*eta
      nqnyh = nq*ldyh
!-----------------------------------------------------------------------
! this section computes the predicted values by effectively
! multiplying the yh array by the pascal triangle matrix.
! dvset is called to calculate all integration coefficients.
! rc is the ratio of new to old values of the coefficient h/el(2)=h/l1.
!-----------------------------------------------------------------------
 400  tn = tn + h
      i1 = nqnyh + 1
      do jb = 1 , nq
         i1 = i1 - ldyh
         do i = i1 , nqnyh
            yh1(i) = yh1(i) + yh1(i+ldyh)
         enddo
      enddo
      call dvset
      rl1 = one/el(2)
      rc = rc*(rl1/prl1)
      prl1 = rl1
!
! call the nonlinear system solver. ------------------------------------
!
      call vnls(y,yh,ldyh,vsav,savf,ewt,acor,iwm,wm,f,jac,psol,nflag,   &
                rpar,ipar)
!
      if ( nflag==0 ) then
!-----------------------------------------------------------------------
! the corrector has converged (nflag = 0).  the local error test is
! made and control passes to statement 500 if it fails.
!-----------------------------------------------------------------------
         dsm = acnrm/tq(2)
         if ( dsm>one ) then
!-----------------------------------------------------------------------
! the error test failed.  kflag keeps track of multiple failures.
! restore tn and the yh array to their previous values, and prepare
! to try the step again.  compute the optimum step size for the
! same order.  after repeated failures, h is forced to decrease
! more rapidly.
!-----------------------------------------------------------------------
            kflag = kflag - 1
            netf = netf + 1
            nflag = -2
            tn = told
            i1 = nqnyh + 1
            do jb = 1 , nq
               i1 = i1 - ldyh
               do i = i1 , nqnyh
                  yh1(i) = yh1(i) - yh1(i+ldyh)
               enddo
            enddo
            if ( abs(h)<=hmin*onepsm ) then
!-----------------------------------------------------------------------
! all returns are made through this section.
! on a successful return, etamax is reset and acor is scaled.
!-----------------------------------------------------------------------
               kflag = -1
               goto 600
            else
               etamax = one
               if ( kflag>kfc ) then
! compute ratio of new h to current h at the current order. ------------
                  flotl = real(l)
                  eta = one/((bias2*dsm)**(one/flotl)+addon)
                  eta = max(eta,hmin/abs(h),etamin)
                  if ( (kflag<=-2) .and. (eta>etamxf) ) eta = etamxf
                  goto 300
!-----------------------------------------------------------------------
! control reaches this section if 3 or more consecutive failures
! have occurred.  it is assumed that the elements of the yh array
! have accumulated errors of the wrong order.  the order is reduced
! by one, if possible.  then h is reduced by a factor of 0.1 and
! the step is retried.  after a total of 7 consecutive failures,
! an exit is taken with kflag = -1.
!-----------------------------------------------------------------------
               elseif ( kflag==kfh ) then
                  kflag = -1
                  goto 600
               elseif ( nq==1 ) then
                  eta = max(etamin,hmin/abs(h))
                  h = h*eta
                  hscal = h
                  tau(1) = h
                  call f(n,tn,y,savf,rpar,ipar)
                  nfe = nfe + 1
                  do i = 1 , n
                     yh(i,2) = h*savf(i)
                  enddo
                  nqwait = 10
                  goto 400
               else
                  eta = max(etamin,hmin/abs(h))
                  call dvjust(yh,ldyh,-1)
                  l = nq
                  nq = nq - 1
                  nqwait = l
                  goto 300
               endif
            endif
         else
!-----------------------------------------------------------------------
! after a successful step, update the yh and tau arrays and decrement
! nqwait.  if nqwait is then 1 and nq .lt. maxord, then acor is saved
! for use in a possible order increase on the next step.
! if etamax = 1 (a failure occurred this step), keep nqwait .ge. 2.
!-----------------------------------------------------------------------
            kflag = 0
            nst = nst + 1
            hu = h
            nqu = nq
            do iback = 1 , nq
               i = l - iback
               tau(i+1) = tau(i)
            enddo
            tau(1) = h
            do j = 1 , l
               call daxpy(n,el(j),acor,1,yh(1,j),1)
            enddo
            nqwait = nqwait - 1
            if ( (l/=lmax) .and. (nqwait==1) ) then
               call dcopy(n,acor,1,yh(1,lmax),1)
               conp = tq(5)
            endif
            if ( etamax/=one ) then
!-----------------------------------------------------------------------
! if nqwait = 0, an increase or decrease in order by one is considered.
! factors etaq, etaqm1, etaqp1 are computed by which h could
! be multiplied at order q, q-1, or q+1, respectively.
! the largest of these is determined, and the new order and
! step size set accordingly.
! a change of h or nq is made only if h increases by at least a
! factor of thresh.  if an order change is considered and rejected,
! then nqwait is set to 2 (reconsider it after 2 steps).
!-----------------------------------------------------------------------
! compute ratio of new h to current h at the current order. ------------
               flotl = real(l)
               etaq = one/((bias2*dsm)**(one/flotl)+addon)
               if ( nqwait==0 ) then
                  nqwait = 2
                  etaqm1 = zero
                  if ( nq/=1 ) then
! compute ratio of new h to current h at the current order less one. ---
                     ddn = dvnorm(n,yh(1,l),ewt)/tq(1)
                     etaqm1 = one/((bias1*ddn)**(one/(flotl-one))+addon)
                  endif
                  etaqp1 = zero
                  if ( l/=lmax ) then
! compute ratio of new h to current h at current order plus one. -------
                     cnquot = (tq(5)/conp)*(h/tau(2))**l
                     do i = 1 , n
                        savf(i) = acor(i) - cnquot*yh(i,lmax)
                     enddo
                     dup = dvnorm(n,savf,ewt)/tq(3)
                     etaqp1 = one/((bias3*dup)**(one/(flotl+one))+addon)
                  endif
                  if ( etaq<etaqp1 ) then
                     if ( etaqp1<=etaqm1 ) goto 420
                     eta = etaqp1
                     newq = nq + 1
                     call dcopy(n,acor,1,yh(1,lmax),1)
                     goto 450
                  elseif ( etaq<etaqm1 ) then
                     goto 420
                  endif
               endif
               eta = etaq
               newq = nq
               goto 450
            else
               if ( nqwait<2 ) nqwait = 2
               newq = nq
               newh = 0
               eta = one
               hnew = h
               goto 500
            endif
 420        eta = etaqm1
            newq = nq - 1
         endif
! test tentative new h against thresh, etamax, and hmxi, then exit. ----
 450     if ( eta<thresh .or. etamax==one ) then
            newq = nq
            newh = 0
            eta = one
            hnew = h
         else
            eta = min(eta,etamax)
            eta = eta/max(one,abs(h)*hmxi*eta)
            newh = 1
            hnew = h*eta
         endif
      else
!-----------------------------------------------------------------------
! the vnls routine failed to achieve convergence (nflag .ne. 0).
! the yh array is retracted to its values before prediction.
! the step size h is reduced and the step is retried, if possible.
! otherwise, an error exit is taken.
!-----------------------------------------------------------------------
         ncf = ncf + 1
         ncfn = ncfn + 1
         etamax = one
         tn = told
         i1 = nqnyh + 1
         do jb = 1 , nq
            i1 = i1 - ldyh
            do i = i1 , nqnyh
               yh1(i) = yh1(i) - yh1(i+ldyh)
            enddo
         enddo
         if ( nflag<-1 ) then
            if ( nflag==-2 ) kflag = -3
            if ( nflag==-3 ) kflag = -4
            goto 600
         elseif ( abs(h)<=hmin*onepsm ) then
            kflag = -2
            goto 600
         elseif ( ncf==mxncf ) then
            kflag = -2
            goto 600
         else
            eta = etacf
            eta = max(eta,hmin/abs(h))
            nflag = -1
            goto 300
         endif
      endif
 500  etamax = etamx3
      if ( nst<=10 ) etamax = etamx2
      r = one/tq(2)
      call dscal(n,r,acor,1)
 600  jstart = 1

      end subroutine dvstep

      subroutine dvset()
      implicit none
!-----------------------------------------------------------------------
! call sequence communication: none
! common block variables accessed:
!     /dvod01/ -- el(13), h, tau(13), tq(5), l(= nq + 1),
!                 meth, nq, nqwait
!
! subroutines called by dvset: none
! function routines called by dvset: none
!-----------------------------------------------------------------------
! dvset is called by dvstep and sets coefficients for use there.
!
! for each order nq, the coefficients in el are calculated by use of
!  the generating polynomial lambda(x), with coefficients el(i).
!      lambda(x) = el(1) + el(2)*x + ... + el(nq+1)*(x**nq).
! for the backward differentiation formulas,
!                                     nq-1
!      lambda(x) = (1 + x/xi*(nq)) * product (1 + x/xi(i) ) .
!                                     i = 1
! for the adams formulas,
!                              nq-1
!      (d/dx) lambda(x) = c * product (1 + x/xi(i) ) ,
!                              i = 1
!      lambda(-1) = 0,    lambda(0) = 1,
! where c is a normalization constant.
! in both cases, xi(i) is defined by
!      h*xi(i) = t sub n  -  t sub (n-i)
!              = h + tau(1) + tau(2) + ... tau(i-1).
!
!
! in addition to variables described previously, communication
! with dvset uses the following:
!   tau    = a vector of length 13 containing the past nq values
!            of h.
!   el     = a vector of length 13 in which vset stores the
!            coefficients for the corrector formula.
!   tq     = a vector of length 5 in which vset stores constants
!            used for the convergence test, the error test, and the
!            selection of h at a new order.
!   meth   = the basic method indicator.
!   nq     = the current order.
!   l      = nq + 1, the length of the vector stored in el, and
!            the number of columns of the yh array being used.
!   nqwait = a counter controlling the frequency of order changes.
!            an order change is about to be considered if nqwait = 1.
!-----------------------------------------------------------------------
!
! type declarations for labeled common block dvod01 --------------------
!
      double precision acnrm , ccmxj , conp , crate , drc , el , eta ,  &
                       etamax , h , hmin , hmxi , hnew , hscal , prl1 , &
                       rc , rl1 , tau , tq , tn , uround
      integer icf , init , ipup , jcur , jstart , jsv , kflag , kuth ,  &
              l , lmax , lyh , lewt , lacor , lsavf , lwm , liwm ,      &
              locjs , maxord , meth , miter , msbj , mxhnil , mxstep ,  &
              n , newh , newq , nhnil , nq , nqnyh , nqwait , nslj ,    &
              nslp , nyh
!
! type declarations for local variables --------------------------------
!
      double precision ahatn0 , alph0 , cnqm1 , cortes , csum , elp ,   &
                       em , em0 , floti , flotl , flotnq , hsum , one , &
                       rxi , rxis , s , six , t1 , t2 , t3 , t4 , t5 ,  &
                       t6 , two , xi , zero
      integer i , iback , j , jp1 , nqm1 , nqm2
!
      dimension em(13)
!-----------------------------------------------------------------------
! the following fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
      save cortes , one , six , two , zero
!
      common /dvod01/ acnrm , ccmxj , conp , crate , drc , el(13) ,     &
                      eta , etamax , h , hmin , hmxi , hnew , hscal ,   &
                      prl1 , rc , rl1 , tau(13) , tq(5) , tn , uround , &
                      icf , init , ipup , jcur , jstart , jsv , kflag , &
                      kuth , l , lmax , lyh , lewt , lacor , lsavf ,    &
                      lwm , liwm , locjs , maxord , meth , miter ,      &
                      msbj , mxhnil , mxstep , n , newh , newq , nhnil ,&
                      nq , nqnyh , nqwait , nslj , nslp , nyh
!
      data cortes/0.1d0/
      data one/1.0d0/ , six/6.0d0/ , two/2.0d0/ , zero/0.0d0/
!
      flotl = real(l)
      nqm1 = nq - 1
      nqm2 = nq - 2
      if ( meth==2 ) then
!
! set coefficients for bdf methods. ------------------------------------
         do i = 3 , l
            el(i) = zero
         enddo
         el(1) = one
         el(2) = one
         alph0 = -one
         ahatn0 = -one
         hsum = h
         rxi = one
         rxis = one
         if ( nq/=1 ) then
            do j = 1 , nqm2
! in el, construct coefficients of (1+x/xi(1))*...*(1+x/xi(j+1)). ------
               hsum = hsum + tau(j)
               rxi = h/hsum
               jp1 = j + 1
               alph0 = alph0 - one/real(jp1)
               do iback = 1 , jp1
                  i = (j+3) - iback
                  el(i) = el(i) + el(i-1)*rxi
               enddo
            enddo
            alph0 = alph0 - one/real(nq)
            rxis = -el(2) - alph0
            hsum = hsum + tau(nqm1)
            rxi = h/hsum
            ahatn0 = -el(2) - rxi
            do iback = 1 , nq
               i = (nq+2) - iback
               el(i) = el(i) + el(i-1)*rxis
            enddo
         endif
         t1 = one - ahatn0 + alph0
         t2 = one + real(nq)*t1
         tq(2) = abs(alph0*t2/t1)
         tq(5) = abs(t2/(el(l)*rxi/rxis))
         if ( nqwait==1 ) then
            cnqm1 = rxis/el(l)
            t3 = alph0 + one/real(nq)
            t4 = ahatn0 + rxi
            elp = t3/(one-t4+t3)
            tq(1) = abs(elp/cnqm1)
            hsum = hsum + tau(nq)
            rxi = h/hsum
            t5 = alph0 - one/real(nq+1)
            t6 = ahatn0 - rxi
            elp = t2/(one-t6+t5)
            tq(3) = abs(elp*rxi*(flotl+one)*t5)
         endif
!
! set coefficients for adams methods. ----------------------------------
      elseif ( nq/=1 ) then
         hsum = h
         em(1) = one
         flotnq = flotl - one
         do i = 2 , l
            em(i) = zero
         enddo
         do j = 1 , nqm1
            if ( (j==nqm1) .and. (nqwait==1) ) then
               s = one
               csum = zero
               do i = 1 , nqm1
                  csum = csum + s*em(i)/real(i+1)
                  s = -s
               enddo
               tq(1) = em(nqm1)/(flotnq*csum)
            endif
            rxi = h/hsum
            do iback = 1 , j
               i = (j+2) - iback
               em(i) = em(i) + em(i-1)*rxi
            enddo
            hsum = hsum + tau(j)
         enddo
! compute integral from -1 to 0 of polynomial and of x times it. -------
         s = one
         em0 = zero
         csum = zero
         do i = 1 , nq
            floti = real(i)
            em0 = em0 + s*em(i)/floti
            csum = csum + s*em(i)/(floti+one)
            s = -s
         enddo
! in el, form coefficients of normalized integrated polynomial. --------
         s = one/em0
         el(1) = one
         do i = 1 , nq
            el(i+1) = s*em(i)/real(i)
         enddo
         xi = hsum/h
         tq(2) = xi*em0/csum
         tq(5) = xi/el(l)
         if ( nqwait==1 ) then
! for higher order control constant, multiply polynomial by 1+x/xi(q). -
            rxi = one/xi
            do iback = 1 , nq
               i = (l+1) - iback
               em(i) = em(i) + em(i-1)*rxi
            enddo
! compute integral of polynomial. --------------------------------------
            s = one
            csum = zero
            do i = 1 , l
               csum = csum + s*em(i)/real(i+1)
               s = -s
            enddo
            tq(3) = flotl*em0/csum
         endif
      else
         el(1) = one
         el(2) = one
         tq(1) = one
         tq(2) = two
         tq(3) = six*tq(2)
         tq(5) = one
      endif
      tq(4) = cortes*tq(2)

      end subroutine dvset

      subroutine dvjust(yh,ldyh,iord)
      implicit none
      double precision yh
      integer ldyh , iord
      dimension yh(ldyh,*)
!-----------------------------------------------------------------------
! call sequence input -- yh, ldyh, iord
! call sequence output -- yh
! common block input -- nq, meth, lmax, hscal, tau(13), n
! common block variables accessed:
!     /dvod01/ -- hscal, tau(13), lmax, meth, n, nq,
!
! subroutines called by dvjust: daxpy
! function routines called by dvjust: none
!-----------------------------------------------------------------------
! this subroutine adjusts the yh array on reduction of order,
! and also when the order is increased for the stiff option (meth = 2).
! communication with dvjust uses the following:
! iord  = an integer flag used when meth = 2 to indicate an order
!         increase (iord = +1) or an order decrease (iord = -1).
! hscal = step size h used in scaling of nordsieck array yh.
!         (if iord = +1, dvjust assumes that hscal = tau(1).)
! see references 1 and 2 for details.
!-----------------------------------------------------------------------
!
! type declarations for labeled common block dvod01 --------------------
!
      double precision acnrm , ccmxj , conp , crate , drc , el , eta ,  &
                       etamax , h , hmin , hmxi , hnew , hscal , prl1 , &
                       rc , rl1 , tau , tq , tn , uround
      integer icf , init , ipup , jcur , jstart , jsv , kflag , kuth ,  &
              l , lmax , lyh , lewt , lacor , lsavf , lwm , liwm ,      &
              locjs , maxord , meth , miter , msbj , mxhnil , mxstep ,  &
              n , newh , newq , nhnil , nq , nqnyh , nqwait , nslj ,    &
              nslp , nyh
!
! type declarations for local variables --------------------------------
!
      double precision alph0 , alph1 , hsum , one , prod , t1 , xi ,    &
                       xiold , zero
      integer i , iback , j , jp1 , lp1 , nqm1 , nqm2 , nqp1
!-----------------------------------------------------------------------
! the following fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
      save one , zero
!
      common /dvod01/ acnrm , ccmxj , conp , crate , drc , el(13) ,     &
                      eta , etamax , h , hmin , hmxi , hnew , hscal ,   &
                      prl1 , rc , rl1 , tau(13) , tq(5) , tn , uround , &
                      icf , init , ipup , jcur , jstart , jsv , kflag , &
                      kuth , l , lmax , lyh , lewt , lacor , lsavf ,    &
                      lwm , liwm , locjs , maxord , meth , miter ,      &
                      msbj , mxhnil , mxstep , n , newh , newq , nhnil ,&
                      nq , nqnyh , nqwait , nslj , nslp , nyh
!
      data one/1.0d0/ , zero/0.0d0/
!
      if ( (nq==2) .and. (iord/=1) ) return
      nqm1 = nq - 1
      nqm2 = nq - 2
      if ( meth==2 ) then
!-----------------------------------------------------------------------
! stiff option...
! check to see if the order is being increased or decreased.
!-----------------------------------------------------------------------
         if ( iord==1 ) then
! order increase. ------------------------------------------------------
            do j = 1 , lmax
               el(j) = zero
            enddo
            el(3) = one
            alph0 = -one
            alph1 = one
            prod = one
            xiold = one
            hsum = hscal
            if ( nq/=1 ) then
               do j = 1 , nqm1
! construct coefficients of x*x*(x+xi(1))*...*(x+xi(j)). ---------------
                  jp1 = j + 1
                  hsum = hsum + tau(jp1)
                  xi = hsum/hscal
                  prod = prod*xi
                  alph0 = alph0 - one/real(jp1)
                  alph1 = alph1 + one/xi
                  do iback = 1 , jp1
                     i = (j+4) - iback
                     el(i) = el(i)*xiold + el(i-1)
                  enddo
                  xiold = xi
               enddo
            endif
            t1 = (-alph0-alph1)/prod
! load column l + 1 in yh array. ---------------------------------------
            lp1 = l + 1
            do i = 1 , n
               yh(i,lp1) = t1*yh(i,lmax)
            enddo
! add correction terms to yh array. ------------------------------------
            nqp1 = nq + 1
            do j = 3 , nqp1
               call daxpy(n,el(j),yh(1,lp1),1,yh(1,j),1)
            enddo
         else
! order decrease. ------------------------------------------------------
            do j = 1 , lmax
               el(j) = zero
            enddo
            el(3) = one
            hsum = zero
            do j = 1 , nqm2
! construct coefficients of x*x*(x+xi(1))*...*(x+xi(j)). ---------------
               hsum = hsum + tau(j)
               xi = hsum/hscal
               jp1 = j + 1
               do iback = 1 , jp1
                  i = (j+4) - iback
                  el(i) = el(i)*xi + el(i-1)
               enddo
            enddo
! subtract correction terms from yh array. -----------------------------
            do j = 3 , nq
               do i = 1 , n
                  yh(i,j) = yh(i,j) - yh(i,l)*el(j)
               enddo
            enddo
            return
         endif
!-----------------------------------------------------------------------
! nonstiff option...
! check to see if the order is being increased or decreased.
!-----------------------------------------------------------------------
      elseif ( iord==1 ) then
! order increase. ------------------------------------------------------
! zero out next column in yh array. ------------------------------------
         lp1 = l + 1
         do i = 1 , n
            yh(i,lp1) = zero
         enddo
         return
      else
! order decrease. ------------------------------------------------------
         do j = 1 , lmax
            el(j) = zero
         enddo
         el(2) = one
         hsum = zero
         do j = 1 , nqm2
! construct coefficients of x*(x+xi(1))*...*(x+xi(j)). -----------------
            hsum = hsum + tau(j)
            xi = hsum/hscal
            jp1 = j + 1
            do iback = 1 , jp1
               i = (j+3) - iback
               el(i) = el(i)*xi + el(i-1)
            enddo
         enddo
! construct coefficients of integrated polynomial. ---------------------
         do j = 2 , nqm1
            el(j+1) = real(nq)*el(j)/real(j)
         enddo
! subtract correction terms from yh array. -----------------------------
         do j = 3 , nq
            do i = 1 , n
               yh(i,j) = yh(i,j) - yh(i,l)*el(j)
            enddo
         enddo
         return
      endif

      end subroutine dvjust

      subroutine dvnlsd(y,yh,ldyh,vsav,savf,ewt,acor,iwm,wm,f,jac,pdum, &
                        nflag,rpar,ipar)
      implicit none
      external f , jac , pdum
      double precision y , yh , vsav , savf , ewt , acor , wm , rpar
      integer ldyh , iwm , nflag , ipar
      dimension y(*) , yh(ldyh,*) , vsav(*) , savf(*) , ewt(*) , acor(*)&
                , iwm(*) , wm(*) , rpar(*) , ipar(*)
!-----------------------------------------------------------------------
! call sequence input -- y, yh, ldyh, savf, ewt, acor, iwm, wm,
!                        f, jac, nflag, rpar, ipar
! call sequence output -- yh, acor, wm, iwm, nflag
! common block variables accessed:
!     /dvod01/ acnrm, crate, drc, h, rc, rl1, tq(5), tn, icf,
!                jcur, meth, miter, n, nslp
!     /dvod02/ hu, ncfn, netf, nfe, nje, nlu, nni, nqu, nst
!
! subroutines called by dvnlsd: f, daxpy, dcopy, dscal, dvjac, dvsol
! function routines called by dvnlsd: dvnorm
!-----------------------------------------------------------------------
! subroutine dvnlsd is a nonlinear system solver, which uses functional
! iteration or a chord (modified newton) method.  for the chord method
! direct linear algebraic system solvers are used.  subroutine dvnlsd
! then handles the corrector phase of this integration package.
!
! communication with dvnlsd is done with the following variables. (for
! more details, please see the comments in the driver subroutine.)
!
! y          = the dependent variable, a vector of length n, input.
! yh         = the nordsieck (taylor) array, ldyh by lmax, input
!              and output.  on input, it contains predicted values.
! ldyh       = a constant .ge. n, the first dimension of yh, input.
! vsav       = unused work array.
! savf       = a work array of length n.
! ewt        = an error weight vector of length n, input.
! acor       = a work array of length n, used for the accumulated
!              corrections to the predicted y vector.
! wm,iwm     = real and integer work arrays associated with matrix
!              operations in chord iteration (miter .ne. 0).
! f          = dummy name for user supplied routine for f.
! jac        = dummy name for user supplied jacobian routine.
! pdum       = unused dummy subroutine name.  included for uniformity
!              over collection of integrators.
! nflag      = input/output flag, with values and meanings as follows:
!              input
!                  0 first call for this time step.
!                 -1 convergence failure in previous call to dvnlsd.
!                 -2 error test failure in dvstep.
!              output
!                  0 successful completion of nonlinear solver.
!                 -1 convergence failure or singular matrix.
!                 -2 unrecoverable error in matrix preprocessing
!                    (cannot occur here).
!                 -3 unrecoverable error in solution (cannot occur
!                    here).
! rpar, ipar = dummy names for user's real and integer work arrays.
!
! ipup       = own variable flag with values and meanings as follows:
!              0,            do not update the newton matrix.
!              miter .ne. 0, update newton matrix, because it is the
!                            initial step, order was changed, the error
!                            test failed, or an update is indicated by
!                            the scalar rc or step counter nst.
!
! for more details, see comments in driver subroutine.
!-----------------------------------------------------------------------
! type declarations for labeled common block dvod01 --------------------
!
      double precision acnrm , ccmxj , conp , crate , drc , el , eta ,  &
                       etamax , h , hmin , hmxi , hnew , hscal , prl1 , &
                       rc , rl1 , tau , tq , tn , uround
      integer icf , init , ipup , jcur , jstart , jsv , kflag , kuth ,  &
              l , lmax , lyh , lewt , lacor , lsavf , lwm , liwm ,      &
              locjs , maxord , meth , miter , msbj , mxhnil , mxstep ,  &
              n , newh , newq , nhnil , nq , nqnyh , nqwait , nslj ,    &
              nslp , nyh
!
! type declarations for labeled common block dvod02 --------------------
!
      double precision hu
      integer ncfn , netf , nfe , nje , nlu , nni , nqu , nst
!
! type declarations for local variables --------------------------------
!
      double precision ccmax , crdown , cscale , dcon , del , delp ,    &
                       one , rdiv , two , zero
      integer i , ierpj , iersl , m , maxcor , msbp
!
! type declaration for function subroutines called ---------------------
!
!      double precision dvnorm
!-----------------------------------------------------------------------
! the following fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
      save ccmax , crdown , maxcor , msbp , rdiv , one , two , zero
!
      common /dvod01/ acnrm , ccmxj , conp , crate , drc , el(13) ,     &
                      eta , etamax , h , hmin , hmxi , hnew , hscal ,   &
                      prl1 , rc , rl1 , tau(13) , tq(5) , tn , uround , &
                      icf , init , ipup , jcur , jstart , jsv , kflag , &
                      kuth , l , lmax , lyh , lewt , lacor , lsavf ,    &
                      lwm , liwm , locjs , maxord , meth , miter ,      &
                      msbj , mxhnil , mxstep , n , newh , newq , nhnil ,&
                      nq , nqnyh , nqwait , nslj , nslp , nyh
      common /dvod02/ hu , ncfn , netf , nfe , nje , nlu , nni , nqu ,  &
                      nst
!
      data ccmax/0.3d0/ , crdown/0.3d0/ , maxcor/3/ , msbp/20/ ,        &
           rdiv/2.0d0/
      data one/1.0d0/ , two/2.0d0/ , zero/0.0d0/
!-----------------------------------------------------------------------
! on the first step, on a change of method order, or after a
! nonlinear convergence failure with nflag = -2, set ipup = miter
! to force a jacobian update when miter .ne. 0.
!-----------------------------------------------------------------------
      if ( jstart==0 ) nslp = 0
      if ( nflag==0 ) icf = 0
      if ( nflag==-2 ) ipup = miter
      if ( (jstart==0) .or. (jstart==-1) ) ipup = miter
! if this is functional iteration, set crate .eq. 1 and drop to 220
      if ( miter==0 ) then
         crate = one
         goto 100
      endif
!-----------------------------------------------------------------------
! rc is the ratio of new to old values of the coefficient h/el(2)=h/l1.
! when rc differs from 1 by more than ccmax, ipup is set to miter
! to force dvjac to be called, if a jacobian is involved.
! in any case, dvjac is called at least every msbp steps.
!-----------------------------------------------------------------------
      drc = abs(rc-one)
      if ( drc>ccmax .or. nst>=nslp+msbp ) ipup = miter
!-----------------------------------------------------------------------
! up to maxcor corrector iterations are taken.  a convergence test is
! made on the r.m.s. norm of each correction, weighted by the error
! weight vector ewt.  the sum of the corrections is accumulated in the
! vector acor(i).  the yh array is not altered in the corrector loop.
!-----------------------------------------------------------------------
 100  m = 0
      delp = zero
      call dcopy(n,yh(1,1),1,y,1)
      call f(n,tn,y,savf,rpar,ipar)
      nfe = nfe + 1
      if ( ipup>0 ) then
!-----------------------------------------------------------------------
! if indicated, the matrix p = i - h*rl1*j is reevaluated and
! preprocessed before starting the corrector iteration.  ipup is set
! to 0 as an indicator that this has been done.
!-----------------------------------------------------------------------
         call dvjac(y,yh,ldyh,ewt,acor,savf,wm,iwm,f,jac,ierpj,rpar,    &
                    ipar)
         ipup = 0
         rc = one
         drc = zero
         crate = one
         nslp = nst
! if matrix is singular, take error return to force cut in step size. --
         if ( ierpj/=0 ) goto 400
      endif
      do i = 1 , n
         acor(i) = zero
      enddo
! this is a looping point for the corrector iteration. -----------------
 200  if ( miter/=0 ) then
!-----------------------------------------------------------------------
! in the case of the chord method, compute the corrector error,
! and solve the linear system with that as right-hand side and
! p as coefficient matrix.  the correction is scaled by the factor
! 2/(1+rc) to account for changes in h*rl1 since the last dvjac call.
!-----------------------------------------------------------------------
         do i = 1 , n
            y(i) = (rl1*h)*savf(i) - (rl1*yh(i,2)+acor(i))
         enddo
         call dvsol(wm,iwm,y,iersl)
         nni = nni + 1
         if ( iersl>0 ) goto 300
         if ( meth==2 .and. rc/=one ) then
            cscale = two/(one+rc)
            call dscal(n,cscale,y,1)
         endif
         del = dvnorm(n,y,ewt)
         call daxpy(n,one,y,1,acor,1)
         do i = 1 , n
            y(i) = yh(i,1) + acor(i)
         enddo
      else
!-----------------------------------------------------------------------
! in the case of functional iteration, update y directly from
! the result of the last function evaluation.
!-----------------------------------------------------------------------
         do i = 1 , n
            savf(i) = rl1*(h*savf(i)-yh(i,2))
         enddo
         do i = 1 , n
            y(i) = savf(i) - acor(i)
         enddo
         del = dvnorm(n,y,ewt)
         do i = 1 , n
            y(i) = yh(i,1) + savf(i)
         enddo
         call dcopy(n,savf,1,acor,1)
      endif
!-----------------------------------------------------------------------
! test for convergence.  if m .gt. 0, an estimate of the convergence
! rate constant is stored in crate, and this is used in the test.
!-----------------------------------------------------------------------
      if ( m/=0 ) crate = max(crdown*crate,del/delp)
      dcon = del*min(one,crate)/tq(4)
      if ( dcon<=one ) then
!
! return for successful step. ------------------------------------------
         nflag = 0
         jcur = 0
         icf = 0
         if ( m==0 ) acnrm = del
         if ( m>0 ) acnrm = dvnorm(n,acor,ewt)
         goto 99999
      else
         m = m + 1
         if ( m/=maxcor ) then
            if ( m<2 .or. del<=rdiv*delp ) then
               delp = del
               call f(n,tn,y,savf,rpar,ipar)
               nfe = nfe + 1
               goto 200
            endif
         endif
      endif
!
 300  if ( miter/=0 .and. jcur/=1 ) then
         icf = 1
         ipup = miter
         goto 100
      endif
!
 400  nflag = -1
      icf = 2
      ipup = miter
      return

99999 end subroutine dvnlsd

      subroutine dvjac(y,yh,ldyh,ewt,ftem,savf,wm,iwm,f,jac,ierpj,rpar, &
                       ipar)
      implicit none
      external f , jac
      double precision y , yh , ewt , ftem , savf , wm , rpar
      integer ldyh , iwm , ierpj , ipar
      dimension y(*) , yh(ldyh,*) , ewt(*) , ftem(*) , savf(*) , wm(*) ,&
                iwm(*) , rpar(*) , ipar(*)
!-----------------------------------------------------------------------
! call sequence input -- y, yh, ldyh, ewt, ftem, savf, wm, iwm,
!                        f, jac, rpar, ipar
! call sequence output -- wm, iwm, ierpj
! common block variables accessed:
!     /dvod01/  ccmxj, drc, h, rl1, tn, uround, icf, jcur, locjs,
!               miter, msbj, n, nslj
!     /dvod02/  nfe, nst, nje, nlu
!
! subroutines called by dvjac: f, jac, dacopy, dcopy, dgbfa, dgefa,
!                              dscal
! function routines called by dvjac: dvnorm
!-----------------------------------------------------------------------
! dvjac is called by dvnlsd to compute and process the matrix
! p = i - h*rl1*j , where j is an approximation to the jacobian.
! here j is computed by the user-supplied routine jac if
! miter = 1 or 4, or by finite differencing if miter = 2, 3, or 5.
! if miter = 3, a diagonal approximation to j is used.
! if jsv = -1, j is computed from scratch in all cases.
! if jsv = 1 and miter = 1, 2, 4, or 5, and if the saved value of j is
! considered acceptable, then p is constructed from the saved j.
! j is stored in wm and replaced by p.  if miter .ne. 3, p is then
! subjected to lu decomposition in preparation for later solution
! of linear systems with p as coefficient matrix. this is done
! by dgefa if miter = 1 or 2, and by dgbfa if miter = 4 or 5.
!
! communication with dvjac is done with the following variables.  (for
! more details, please see the comments in the driver subroutine.)
! y          = vector containing predicted values on entry.
! yh         = the nordsieck array, an ldyh by lmax array, input.
! ldyh       = a constant .ge. n, the first dimension of yh, input.
! ewt        = an error weight vector of length n.
! savf       = array containing f evaluated at predicted y, input.
! wm         = real work space for matrices.  in the output, it contains
!              the inverse diagonal matrix if miter = 3 and the lu
!              decomposition of p if miter is 1, 2 , 4, or 5.
!              storage of matrix elements starts at wm(3).
!              storage of the saved jacobian starts at wm(locjs).
!              wm also contains the following matrix-related data:
!              wm(1) = sqrt(uround), used in numerical jacobian step.
!              wm(2) = h*rl1, saved for later use if miter = 3.
! iwm        = integer work space containing pivot information,
!              starting at iwm(31), if miter is 1, 2, 4, or 5.
!              iwm also contains band parameters ml = iwm(1) and
!              mu = iwm(2) if miter is 4 or 5.
! f          = dummy name for the user supplied subroutine for f.
! jac        = dummy name for the user supplied jacobian subroutine.
! rpar, ipar = dummy names for user's real and integer work arrays.
! rl1        = 1/el(2) (input).
! ierpj      = output error flag,  = 0 if no trouble, 1 if the p
!              matrix is found to be singular.
! jcur       = output flag to indicate whether the jacobian matrix
!              (or approximation) is now current.
!              jcur = 0 means j is not current.
!              jcur = 1 means j is current.
!-----------------------------------------------------------------------
!
! type declarations for labeled common block dvod01 --------------------
!
      double precision acnrm , ccmxj , conp , crate , drc , el , eta ,  &
                       etamax , h , hmin , hmxi , hnew , hscal , prl1 , &
                       rc , rl1 , tau , tq , tn , uround
      integer icf , init , ipup , jcur , jstart , jsv , kflag , kuth ,  &
              l , lmax , lyh , lewt , lacor , lsavf , lwm , liwm ,      &
              locjs , maxord , meth , miter , msbj , mxhnil , mxstep ,  &
              n , newh , newq , nhnil , nq , nqnyh , nqwait , nslj ,    &
              nslp , nyh
!
! type declarations for labeled common block dvod02 --------------------
!
      double precision hu
      integer ncfn , netf , nfe , nje , nlu , nni , nqu , nst
!
! type declarations for local variables --------------------------------
!
      double precision con , di , fac , hrl1 , one , pt1 , r , r0 ,     &
                       srur , thou , yi , yj , yjj , zero
      integer i , i1 , i2 , ier , ii , j , j1 , jj , jok , lenp , mba , &
              mband , meb1 , meband , ml , ml3 , mu , np1
!
! type declaration for function subroutines called ---------------------
!
!      double precision dvnorm
!-----------------------------------------------------------------------
! the following fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this subroutine.
!-----------------------------------------------------------------------
      save one , pt1 , thou , zero
!-----------------------------------------------------------------------
      common /dvod01/ acnrm , ccmxj , conp , crate , drc , el(13) ,     &
                      eta , etamax , h , hmin , hmxi , hnew , hscal ,   &
                      prl1 , rc , rl1 , tau(13) , tq(5) , tn , uround , &
                      icf , init , ipup , jcur , jstart , jsv , kflag , &
                      kuth , l , lmax , lyh , lewt , lacor , lsavf ,    &
                      lwm , liwm , locjs , maxord , meth , miter ,      &
                      msbj , mxhnil , mxstep , n , newh , newq , nhnil ,&
                      nq , nqnyh , nqwait , nslj , nslp , nyh
      common /dvod02/ hu , ncfn , netf , nfe , nje , nlu , nni , nqu ,  &
                      nst
!
      data one/1.0d0/ , thou/1000.0d0/ , zero/0.0d0/ , pt1/0.1d0/
!
      ierpj = 0
      hrl1 = h*rl1
! see whether j should be evaluated (jok = -1) or not (jok = 1). -------
      jok = jsv
      if ( jsv==1 ) then
         if ( nst==0 .or. nst>nslj+msbj ) jok = -1
         if ( icf==1 .and. drc<ccmxj ) jok = -1
         if ( icf==2 ) jok = -1
      endif
! end of setting jok. --------------------------------------------------
!
      if ( jok==-1 .and. miter==1 ) then
! if jok = -1 and miter = 1, call jac to evaluate jacobian. ------------
         nje = nje + 1
         nslj = nst
         jcur = 1
         lenp = n*n
         do i = 1 , lenp
            wm(i+2) = zero
         enddo
         call jac(n,tn,y,0,0,wm(3),n,rpar,ipar)
         if ( jsv==1 ) call dcopy(lenp,wm(3),1,wm(locjs),1)
      endif
!
      if ( jok==-1 .and. miter==2 ) then
! if miter = 2, make n calls to f to approximate the jacobian. ---------
         nje = nje + 1
         nslj = nst
         jcur = 1
         fac = dvnorm(n,savf,ewt)
         r0 = thou*abs(h)*uround*real(n)*fac
         if ( r0==zero ) r0 = one
         srur = wm(1)
         j1 = 2
         do j = 1 , n
            yj = y(j)
            r = max(srur*abs(yj),r0/ewt(j))
            y(j) = y(j) + r
            fac = one/r
            call f(n,tn,y,ftem,rpar,ipar)
            do i = 1 , n
               wm(i+j1) = (ftem(i)-savf(i))*fac
            enddo
            y(j) = yj
            j1 = j1 + n
         enddo
         nfe = nfe + n
         lenp = n*n
         if ( jsv==1 ) call dcopy(lenp,wm(3),1,wm(locjs),1)
      endif
!
      if ( jok==1 .and. (miter==1 .or. miter==2) ) then
         jcur = 0
         lenp = n*n
         call dcopy(lenp,wm(locjs),1,wm(3),1)
      endif
!
      if ( miter==1 .or. miter==2 ) then
! multiply jacobian by scalar, add identity, and do lu decomposition. --
         con = -hrl1
         call dscal(lenp,con,wm(3),1)
         j = 3
         np1 = n + 1
         do i = 1 , n
            wm(j) = wm(j) + one
            j = j + np1
         enddo
         nlu = nlu + 1
         call dgefa(wm(3),n,n,iwm(31),ier)
         if ( ier/=0 ) ierpj = 1
         return
      endif
! end of code block for miter = 1 or 2. --------------------------------
!
      if ( miter==3 ) then
! if miter = 3, construct a diagonal approximation to j and p. ---------
         nje = nje + 1
         jcur = 1
         wm(2) = hrl1
         r = rl1*pt1
         do i = 1 , n
            y(i) = y(i) + r*(h*savf(i)-yh(i,2))
         enddo
         call f(n,tn,y,wm(3),rpar,ipar)
         nfe = nfe + 1
         do i = 1 , n
            r0 = h*savf(i) - yh(i,2)
            di = pt1*r0 - h*(wm(i+2)-savf(i))
            wm(i+2) = one
            if ( abs(r0)>=uround/ewt(i) ) then
               if ( abs(di)==zero ) goto 50
               wm(i+2) = pt1*r0/di
            endif
         enddo
         return
 50      ierpj = 1
         return
      endif
! end of code block for miter = 3. -------------------------------------
!
! set constants for miter = 4 or 5. ------------------------------------
      ml = iwm(1)
      mu = iwm(2)
      ml3 = ml + 3
      mband = ml + mu + 1
      meband = mband + ml
      lenp = meband*n
!
      if ( jok==-1 .and. miter==4 ) then
! if jok = -1 and miter = 4, call jac to evaluate jacobian. ------------
         nje = nje + 1
         nslj = nst
         jcur = 1
         do i = 1 , lenp
            wm(i+2) = zero
         enddo
         call jac(n,tn,y,ml,mu,wm(ml3),meband,rpar,ipar)
         if ( jsv==1 ) call dacopy(mband,n,wm(ml3),meband,wm(locjs),    &
                                   mband)
      endif
!
      if ( jok==-1 .and. miter==5 ) then
! if miter = 5, make ml+mu+1 calls to f to approximate the jacobian. ---
         nje = nje + 1
         nslj = nst
         jcur = 1
         mba = min(mband,n)
         meb1 = meband - 1
         srur = wm(1)
         fac = dvnorm(n,savf,ewt)
         r0 = thou*abs(h)*uround*real(n)*fac
         if ( r0==zero ) r0 = one
         do j = 1 , mba
            do i = j , n , mband
               yi = y(i)
               r = max(srur*abs(yi),r0/ewt(i))
               y(i) = y(i) + r
            enddo
            call f(n,tn,y,ftem,rpar,ipar)
            do jj = j , n , mband
               y(jj) = yh(jj,1)
               yjj = y(jj)
               r = max(srur*abs(yjj),r0/ewt(jj))
               fac = one/r
               i1 = max(jj-mu,1)
               i2 = min(jj+ml,n)
               ii = jj*meb1 - ml + 2
               do i = i1 , i2
                  wm(ii+i) = (ftem(i)-savf(i))*fac
               enddo
            enddo
         enddo
         nfe = nfe + mba
         if ( jsv==1 ) call dacopy(mband,n,wm(ml3),meband,wm(locjs),    &
                                   mband)
      endif
!
      if ( jok==1 ) then
         jcur = 0
         call dacopy(mband,n,wm(locjs),mband,wm(ml3),meband)
      endif
!
! multiply jacobian by scalar, add identity, and do lu decomposition.
      con = -hrl1
      call dscal(lenp,con,wm(3),1)
      ii = mband + 2
      do i = 1 , n
         wm(ii) = wm(ii) + one
         ii = ii + meband
      enddo
      nlu = nlu + 1
      call dgbfa(wm(3),meband,n,ml,mu,iwm(31),ier)
      if ( ier/=0 ) ierpj = 1
! end of code block for miter = 4 or 5. --------------------------------
!

      end subroutine dvjac

      subroutine dacopy(nrow,ncol,a,nrowa,b,nrowb)
      implicit none
      double precision a , b
      integer nrow , ncol , nrowa , nrowb
      dimension a(nrowa,ncol) , b(nrowb,ncol)
!-----------------------------------------------------------------------
! call sequence input -- nrow, ncol, a, nrowa, nrowb
! call sequence output -- b
! common block variables accessed -- none
!
! subroutines called by dacopy: dcopy
! function routines called by dacopy: none
!-----------------------------------------------------------------------
! this routine copies one rectangular array, a, to another, b,
! where a and b may have different row dimensions, nrowa and nrowb.
! the data copied consists of nrow rows and ncol columns.
!-----------------------------------------------------------------------
      integer ic
!
      do ic = 1 , ncol
         call dcopy(nrow,a(1,ic),1,b(1,ic),1)
      enddo
!

      end subroutine dacopy

      subroutine dvsol(wm,iwm,x,iersl)
      implicit none
      double precision wm , x
      integer iwm , iersl
      dimension wm(*) , iwm(*) , x(*)
!-----------------------------------------------------------------------
! call sequence input -- wm, iwm, x
! call sequence output -- x, iersl
! common block variables accessed:
!     /dvod01/ -- h, rl1, miter, n
!
! subroutines called by dvsol: dgesl, dgbsl
! function routines called by dvsol: none
!-----------------------------------------------------------------------
! this routine manages the solution of the linear system arising from
! a chord iteration.  it is called if miter .ne. 0.
! if miter is 1 or 2, it calls dgesl to accomplish this.
! if miter = 3 it updates the coefficient h*rl1 in the diagonal
! matrix, and then computes the solution.
! if miter is 4 or 5, it calls dgbsl.
! communication with dvsol uses the following variables:
! wm    = real work space containing the inverse diagonal matrix if
!         miter = 3 and the lu decomposition of the matrix otherwise.
!         storage of matrix elements starts at wm(3).
!         wm also contains the following matrix-related data:
!         wm(1) = sqrt(uround) (not used here),
!         wm(2) = hrl1, the previous value of h*rl1, used if miter = 3.
! iwm   = integer work space containing pivot information, starting at
!         iwm(31), if miter is 1, 2, 4, or 5.  iwm also contains band
!         parameters ml = iwm(1) and mu = iwm(2) if miter is 4 or 5.
! x     = the right-hand side vector on input, and the solution vector
!         on output, of length n.
! iersl = output flag.  iersl = 0 if no trouble occurred.
!         iersl = 1 if a singular matrix arose with miter = 3.
!-----------------------------------------------------------------------
!
! type declarations for labeled common block dvod01 --------------------
!
      double precision acnrm , ccmxj , conp , crate , drc , el , eta ,  &
                       etamax , h , hmin , hmxi , hnew , hscal , prl1 , &
                       rc , rl1 , tau , tq , tn , uround
      integer icf , init , ipup , jcur , jstart , jsv , kflag , kuth ,  &
              l , lmax , lyh , lewt , lacor , lsavf , lwm , liwm ,      &
              locjs , maxord , meth , miter , msbj , mxhnil , mxstep ,  &
              n , newh , newq , nhnil , nq , nqnyh , nqwait , nslj ,    &
              nslp , nyh
!
! type declarations for local variables --------------------------------
!
      integer i , meband , ml , mu
      double precision di , hrl1 , one , phrl1 , r , zero
!-----------------------------------------------------------------------
! the following fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
      save one , zero
!
      common /dvod01/ acnrm , ccmxj , conp , crate , drc , el(13) ,     &
                      eta , etamax , h , hmin , hmxi , hnew , hscal ,   &
                      prl1 , rc , rl1 , tau(13) , tq(5) , tn , uround , &
                      icf , init , ipup , jcur , jstart , jsv , kflag , &
                      kuth , l , lmax , lyh , lewt , lacor , lsavf ,    &
                      lwm , liwm , locjs , maxord , meth , miter ,      &
                      msbj , mxhnil , mxstep , n , newh , newq , nhnil ,&
                      nq , nqnyh , nqwait , nslj , nslp , nyh
!
      data one/1.0d0/ , zero/0.0d0/
!
      iersl = 0
      select case (miter)
      case (3)
!
         phrl1 = wm(2)
         hrl1 = h*rl1
         wm(2) = hrl1
         if ( hrl1/=phrl1 ) then
            r = hrl1/phrl1
            do i = 1 , n
               di = one - r*(one-one/wm(i+2))
               if ( abs(di)==zero ) goto 100
               wm(i+2) = one/di
            enddo
         endif
!
         do i = 1 , n
            x(i) = wm(i+2)*x(i)
         enddo
         return
      case (4,5)
!
         ml = iwm(1)
         mu = iwm(2)
         meband = 2*ml + mu + 1
         call dgbsl(wm(3),meband,n,ml,mu,iwm(31),x,0)
         goto 99999
      case default
         call dgesl(wm(3),n,n,iwm(31),x,0)
         return
      end select
 100  iersl = 1
      return

99999 end subroutine dvsol

      subroutine dvsrco(rsav,isav,job)
      implicit none
      double precision rsav
      integer isav , job
      dimension rsav(*) , isav(*)
!-----------------------------------------------------------------------
! call sequence input -- rsav, isav, job
! call sequence output -- rsav, isav
! common block variables accessed -- all of /dvod01/ and /dvod02/
!
! subroutines/functions called by dvsrco: none
!-----------------------------------------------------------------------
! this routine saves or restores (depending on job) the contents of the
! common blocks dvod01 and dvod02, which are used internally by dvode.
!
! rsav = real array of length 49 or more.
! isav = integer array of length 41 or more.
! job  = flag indicating to save or restore the common blocks:
!        job  = 1 if common is to be saved (written to rsav/isav).
!        job  = 2 if common is to be restored (read from rsav/isav).
!        a call with job = 2 presumes a prior call with job = 1.
!-----------------------------------------------------------------------
      double precision rvod1 , rvod2
      integer ivod1 , ivod2
      integer i , leniv1 , leniv2 , lenrv1 , lenrv2
!-----------------------------------------------------------------------
! the following fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
      save lenrv1 , leniv1 , lenrv2 , leniv2
!
      common /dvod01/ rvod1(48) , ivod1(33)
      common /dvod02/ rvod2(1) , ivod2(8)
      data lenrv1/48/ , leniv1/33/ , lenrv2/1/ , leniv2/8/
!
      if ( job==2 ) then
!
         do i = 1 , lenrv1
            rvod1(i) = rsav(i)
         enddo
         do i = 1 , lenrv2
            rvod2(i) = rsav(lenrv1+i)
         enddo
!
         do i = 1 , leniv1
            ivod1(i) = isav(i)
         enddo
         do i = 1 , leniv2
            ivod2(i) = isav(leniv1+i)
         enddo
         goto 99999
      endif
      do i = 1 , lenrv1
         rsav(i) = rvod1(i)
      enddo
      do i = 1 , lenrv2
         rsav(lenrv1+i) = rvod2(i)
      enddo
!
      do i = 1 , leniv1
         isav(i) = ivod1(i)
      enddo
      do i = 1 , leniv2
         isav(leniv1+i) = ivod2(i)
      enddo
!
      return
!

99999 end subroutine dvsrco

      subroutine dewset(n,itol,rtol,atol,ycur,ewt)
      implicit none
!***begin prologue  dewset
!***subsidiary
!***purpose  set error weight vector.
!***type      double precision (sewset-s, dewset-d)
!***author  hindmarsh, alan c., (llnl)
!***description
!
!  this subroutine sets the error weight vector ewt according to
!      ewt(i) = rtol(i)*abs(ycur(i)) + atol(i),  i = 1,...,n,
!  with the subscript on rtol and/or atol possibly replaced by 1 above,
!  depending on the value of itol.
!
!***see also  dlsode
!***routines called  (none)
!***revision history  (yymmdd)
!   791129  date written
!   890501  modified prologue to slatec/ldoc format.  (fnf)
!   890503  minor cosmetic changes.  (fnf)
!   930809  renamed to allow single/double precision versions. (ach)
!***end prologue  dewset
!**end
      integer n , itol
      integer i
      double precision rtol , atol , ycur , ewt
      dimension rtol(*) , atol(*) , ycur(n) , ewt(n)
!

      select case (itol)
      case (2)
         do i = 1 , n
            ewt(i) = rtol(1)*abs(ycur(i)) + atol(i)
         enddo
         return
      case (3)
         do i = 1 , n
            ewt(i) = rtol(i)*abs(ycur(i)) + atol(1)
         enddo
         return
      case (4)
         do i = 1 , n
            ewt(i) = rtol(i)*abs(ycur(i)) + atol(i)
         enddo
         goto 99999
      case default
      end select
      do i = 1 , n
         ewt(i) = rtol(1)*abs(ycur(i)) + atol(1)
      enddo
      return

99999 end subroutine dewset

      double precision function dvnorm(n,v,w)
      implicit none
!***begin prologue  dvnorm
!***subsidiary
!***purpose  weighted root-mean-square vector norm.
!***type      double precision (svnorm-s, dvnorm-d)
!***author  hindmarsh, alan c., (llnl)
!***description
!
!  this function routine computes the weighted root-mean-square norm
!  of the vector of length n contained in the array v, with weights
!  contained in the array w of length n:
!    dvnorm = sqrt( (1/n) * sum( v(i)*w(i) )**2 )
!
!***see also  dlsode
!***routines called  (none)
!***revision history  (yymmdd)
!   791129  date written
!   890501  modified prologue to slatec/ldoc format.  (fnf)
!   890503  minor cosmetic changes.  (fnf)
!   930809  renamed to allow single/double precision versions. (ach)
!***end prologue  dvnorm
!**end
      integer n , i
      double precision v , w , sum
      dimension v(n) , w(n)

      sum = 0.0d0
      do i = 1 , n
         sum = sum + (v(i)*w(i))**2
      enddo
      dvnorm = sqrt(sum/n)

      end function dvnorm

      subroutine xerrwd(msg,nmes,nerr,level,ni,i1,i2,nr,r1,r2)
      implicit none
!***begin prologue  xerrwd
!***subsidiary
!***purpose  write error message with values.
!***category  r3c
!***type      double precision (xerrwv-s, xerrwd-d)
!***author  hindmarsh, alan c., (llnl)
!***description
!
!  subroutines xerrwd, xsetf, xsetun, and the function routine ixsav,
!  as given here, constitute a simplified version of the slatec error
!  handling package.
!
!  all arguments are input arguments.
!
!  msg    = the message (character array).
!  nmes   = the length of msg (number of characters).
!  nerr   = the error number (not used).
!  level  = the error level..
!           0 or 1 means recoverable (control returns to caller).
!           2 means fatal (run is aborted--see note below).
!  ni     = number of integers (0, 1, or 2) to be printed with message.
!  i1,i2  = integers to be printed, depending on ni.
!  nr     = number of reals (0, 1, or 2) to be printed with message.
!  r1,r2  = reals to be printed, depending on nr.
!
!  note..  this routine is machine-dependent and specialized for use
!  in limited context, in the following ways..
!  1. the argument msg is assumed to be of type character, and
!     the message is printed with a format of (1x,a).
!  2. the message is assumed to take only one line.
!     multi-line messages are generated by repeated calls.
!  3. if level = 2, control passes to the statement   stop
!     to abort the run.  this statement may be machine-dependent.
!  4. r1 and r2 are assumed to be in double precision and are printed
!     in d21.13 format.
!
!***routines called  ixsav
!***revision history  (yymmdd)
!   920831  date written
!   921118  replaced mflgsv/lunsav by ixsav. (ach)
!   930329  modified prologue to slatec format. (fnf)
!   930407  changed msg from character*1 array to variable. (fnf)
!   930922  minor cosmetic change. (fnf)
!***end prologue  xerrwd
!
!*internal notes:
!
! for a different default logical unit number, ixsav (or a subsidiary
! routine that it calls) will need to be modified.
! for a different run-abort command, change the statement following
! statement 100 at the end.
!-----------------------------------------------------------------------
! subroutines called by xerrwd.. none
! function routine called by xerrwd.. ixsav
!-----------------------------------------------------------------------
!**end
!
!  declare arguments.
!
      double precision r1 , r2
      integer nmes , nerr , level , ni , i1 , i2 , nr
      character*(*) msg
!
!  declare local variables.
!
      integer lunit , mesflg
!
!  get logical unit number and message print flag.
!

      lunit = ixsav(1,0,.false.)
      mesflg = ixsav(2,0,.false.)
      if ( mesflg/=0 ) then
!
!  write the message.
!
         write (lunit,99001) msg
99001    format (1x,a)
         if ( ni==1 ) write (lunit,99002) i1
99002    format (6x,'in above message,  i1 =',i10)
         if ( ni==2 ) write (lunit,99003) i1 , i2
99003    format (6x,'in above message,  i1 =',i10,3x,'i2 =',i10)
         if ( nr==1 ) write (lunit,99004) r1
99004    format (6x,'in above message,  r1 =',d21.13)
         if ( nr==2 ) write (lunit,99005) r1 , r2
99005    format (6x,'in above,  r1 =',d21.13,3x,'r2 =',d21.13)
      endif
!
!  abort the run if level = 2.
!
      if ( level/=2 ) return
      stop

      end subroutine xerrwd

      subroutine xsetf(mflag)
      implicit none
!***begin prologue  xsetf
!***purpose  reset the error print control flag.
!***category  r3a
!***type      all (xsetf-a)
!***keywords  error control
!***author  hindmarsh, alan c., (llnl)
!***description
!
!   xsetf sets the error print control flag to mflag:
!      mflag=1 means print all messages (the default).
!      mflag=0 means no printing.
!
!***see also  xerrwd, xerrwv
!***references  (none)
!***routines called  ixsav
!***revision history  (yymmdd)
!   921118  date written
!   930329  added slatec format prologue. (fnf)
!   930407  corrected see also section. (fnf)
!   930922  made user-callable, and other cosmetic changes. (fnf)
!***end prologue  xsetf
!
! subroutines called by xsetf.. none
! function routine called by xsetf.. ixsav
!-----------------------------------------------------------------------
!**end
      integer mflag , junk
!

      if ( mflag==0 .or. mflag==1 ) junk = ixsav(2,mflag,.true.)

      end subroutine xsetf

      subroutine xsetun(lun)
      implicit none
!***begin prologue  xsetun
!***purpose  reset the logical unit number for error messages.
!***category  r3b
!***type      all (xsetun-a)
!***keywords  error control
!***description
!
!   xsetun sets the logical unit number for error messages to lun.
!
!***author  hindmarsh, alan c., (llnl)
!***see also  xerrwd, xerrwv
!***references  (none)
!***routines called  ixsav
!***revision history  (yymmdd)
!   921118  date written
!   930329  added slatec format prologue. (fnf)
!   930407  corrected see also section. (fnf)
!   930922  made user-callable, and other cosmetic changes. (fnf)
!***end prologue  xsetun
!
! subroutines called by xsetun.. none
! function routine called by xsetun.. ixsav
!-----------------------------------------------------------------------
!**end
      integer lun , junk
!

      if ( lun>0 ) junk = ixsav(1,lun,.true.)

      end subroutine xsetun

      integer function ixsav(ipar,ivalue,iset)
      implicit none
!***begin prologue  ixsav
!***subsidiary
!***purpose  save and recall error message control parameters.
!***category  r3c
!***type      all (ixsav-a)
!***author  hindmarsh, alan c., (llnl)
!***description
!
!  ixsav saves and recalls one of two error message parameters:
!    lunit, the logical unit number to which messages are printed, and
!    mesflg, the message print flag.
!  this is a modification of the slatec library routine j4save.
!
!  saved local variables..
!   lunit  = logical unit number for messages.  the default is obtained
!            by a call to iumach (may be machine-dependent).
!   mesflg = print control flag..
!            1 means print all messages (the default).
!            0 means no printing.
!
!  on input..
!    ipar   = parameter indicator (1 for lunit, 2 for mesflg).
!    ivalue = the value to be set for the parameter, if iset = .true.
!    iset   = logical flag to indicate whether to read or write.
!             if iset = .true., the parameter will be given
!             the value ivalue.  if iset = .false., the parameter
!             will be unchanged, and ivalue is a dummy argument.
!
!  on return..
!    ixsav = the (old) value of the parameter.
!
!***see also  xerrwd, xerrwv
!***routines called  iumach
!***revision history  (yymmdd)
!   921118  date written
!   930329  modified prologue to slatec format. (fnf)
!   930915  added iumach call to get default output unit.  (ach)
!   930922  minor cosmetic changes. (fnf)
!   010425  type declaration for iumach added. (ach)
!***end prologue  ixsav
!
! subroutines called by ixsav.. none
! function routine called by ixsav.. iumach
!-----------------------------------------------------------------------
!**end
      logical iset
      integer ipar , ivalue
!-----------------------------------------------------------------------
      integer lunit , mesflg
!-----------------------------------------------------------------------
! the following fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this routine.
!-----------------------------------------------------------------------
      save lunit , mesflg
      data lunit/ - 1/ , mesflg/1/
!

      if ( ipar==1 ) then
         if ( lunit==-1 ) lunit = iumach
         ixsav = lunit
         if ( iset ) lunit = ivalue
      endif
!
      if ( ipar==2 ) then
         ixsav = mesflg
         if ( iset ) mesflg = ivalue
      endif
!

      end function ixsav

end module dvode_module
