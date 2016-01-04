module dvode_module

    implicit none

    private

    public :: dvode

contains

!DECK DVODE
      SUBROUTINE DVODE(F,Neq,Y,T,Tout,Itol,Rtol,Atol,Itask,Istate,Iopt, &
                       Rwork,Lrw,Iwork,Liw,JAC,Mf,Rpar,Ipar)
      IMPLICIT NONE
      EXTERNAL F , JAC
      DOUBLE PRECISION Y , T , Tout , Rtol , Atol , Rwork , Rpar
      INTEGER Neq , Itol , Itask , Istate , Iopt , Lrw , Iwork , Liw ,  &
              Mf , Ipar
      DIMENSION Y(*) , Rtol(*) , Atol(*) , Rwork(Lrw) , Iwork(Liw) ,    &
                Rpar(*) , Ipar(*)
!-----------------------------------------------------------------------
! DVODE: Variable-coefficient Ordinary Differential Equation solver,
! with fixed-leading-coefficient implementation.
! This version is in double precision.
!
! DVODE solves the initial value problem for stiff or nonstiff
! systems of first order ODEs,
!     dy/dt = f(t,y) ,  or, in component form,
!     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(NEQ)) (i = 1,...,NEQ).
! DVODE is a package based on the EPISODE and EPISODEB packages, and
! on the ODEPACK user interface standard, with minor modifications.
!-----------------------------------------------------------------------
! Authors:
!               Peter N. Brown and Alan C. Hindmarsh
!               Center for Applied Scientific Computing, L-561
!               Lawrence Livermore National Laboratory
!               Livermore, CA 94551
! and
!               George D. Byrne
!               Illinois Institute of Technology
!               Chicago, IL 60616
!-----------------------------------------------------------------------
! References:
!
! 1. P. N. Brown, G. D. Byrne, and A. C. Hindmarsh, "VODE: A Variable
!    Coefficient ODE Solver," SIAM J. Sci. Stat. Comput., 10 (1989),
!    pp. 1038-1051.  Also, LLNL Report UCRL-98412, June 1988.
! 2. G. D. Byrne and A. C. Hindmarsh, "A Polyalgorithm for the
!    Numerical Solution of Ordinary Differential Equations,"
!    ACM Trans. Math. Software, 1 (1975), pp. 71-96.
! 3. A. C. Hindmarsh and G. D. Byrne, "EPISODE: An Effective Package
!    for the Integration of Systems of Ordinary Differential
!    Equations," LLNL Report UCID-30112, Rev. 1, April 1977.
! 4. G. D. Byrne and A. C. Hindmarsh, "EPISODEB: An Experimental
!    Package for the Integration of Systems of Ordinary Differential
!    Equations with Banded Jacobians," LLNL Report UCID-30132, April
!    1976.
! 5. A. C. Hindmarsh, "ODEPACK, a Systematized Collection of ODE
!    Solvers," in Scientific Computing, R. S. Stepleman et al., eds.,
!    North-Holland, Amsterdam, 1983, pp. 55-64.
! 6. K. R. Jackson and R. Sacks-Davis, "An Alternative Implementation
!    of Variable Step-Size Multistep Formulas for Stiff ODEs," ACM
!    Trans. Math. Software, 6 (1980), pp. 295-318.
!-----------------------------------------------------------------------
! Summary of usage.
!
! Communication between the user and the DVODE package, for normal
! situations, is summarized here.  This summary describes only a subset
! of the full set of options available.  See the full description for
! details, including optional communication, nonstandard options,
! and instructions for special situations.  See also the example
! problem (with program and output) following this summary.
!
! A. First provide a subroutine of the form:
!           SUBROUTINE F (NEQ, T, Y, YDOT, RPAR, IPAR)
!           DOUBLE PRECISION T, Y(NEQ), YDOT(NEQ), RPAR
! which supplies the vector function f by loading YDOT(i) with f(i).
!
! B. Next determine (or guess) whether or not the problem is stiff.
! Stiffness occurs when the Jacobian matrix df/dy has an eigenvalue
! whose real part is negative and large in magnitude, compared to the
! reciprocal of the t span of interest.  If the problem is nonstiff,
! use a method flag MF = 10.  If it is stiff, there are four standard
! choices for MF (21, 22, 24, 25), and DVODE requires the Jacobian
! matrix in some form.  In these cases (MF .gt. 0), DVODE will use a
! saved copy of the Jacobian matrix.  If this is undesirable because of
! storage limitations, set MF to the corresponding negative value
! (-21, -22, -24, -25).  (See full description of MF below.)
! The Jacobian matrix is regarded either as full (MF = 21 or 22),
! or banded (MF = 24 or 25).  In the banded case, DVODE requires two
! half-bandwidth parameters ML and MU.  These are, respectively, the
! widths of the lower and upper parts of the band, excluding the main
! diagonal.  Thus the band consists of the locations (i,j) with
! i-ML .le. j .le. i+MU, and the full bandwidth is ML+MU+1.
!
! C. If the problem is stiff, you are encouraged to supply the Jacobian
! directly (MF = 21 or 24), but if this is not feasible, DVODE will
! compute it internally by difference quotients (MF = 22 or 25).
! If you are supplying the Jacobian, provide a subroutine of the form:
!           SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD, RPAR, IPAR)
!           DOUBLE PRECISION T, Y(NEQ), PD(NROWPD,NEQ), RPAR
! which supplies df/dy by loading PD as follows:
!     For a full Jacobian (MF = 21), load PD(i,j) with df(i)/dy(j),
! the partial derivative of f(i) with respect to y(j).  (Ignore the
! ML and MU arguments in this case.)
!     For a banded Jacobian (MF = 24), load PD(i-j+MU+1,j) with
! df(i)/dy(j), i.e. load the diagonal lines of df/dy into the rows of
! PD from the top down.
!     In either case, only nonzero elements need be loaded.
!
! D. Write a main program which calls subroutine DVODE once for
! each point at which answers are desired.  This should also provide
! for possible use of logical unit 6 for output of error messages
! by DVODE.  On the first call to DVODE, supply arguments as follows:
! F      = Name of subroutine for right-hand side vector f.
!          This name must be declared external in calling program.
! NEQ    = Number of first order ODEs.
! Y      = Array of initial values, of length NEQ.
! T      = The initial value of the independent variable.
! TOUT   = First point where output is desired (.ne. T).
! ITOL   = 1 or 2 according as ATOL (below) is a scalar or array.
! RTOL   = Relative tolerance parameter (scalar).
! ATOL   = Absolute tolerance parameter (scalar or array).
!          The estimated local error in Y(i) will be controlled so as
!          to be roughly less (in magnitude) than
!             EWT(i) = RTOL*abs(Y(i)) + ATOL     if ITOL = 1, or
!             EWT(i) = RTOL*abs(Y(i)) + ATOL(i)  if ITOL = 2.
!          Thus the local error test passes if, in each component,
!          either the absolute error is less than ATOL (or ATOL(i)),
!          or the relative error is less than RTOL.
!          Use RTOL = 0.0 for pure absolute error control, and
!          use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error
!          control.  Caution: Actual (global) errors may exceed these
!          local tolerances, so choose them conservatively.
! ITASK  = 1 for normal computation of output values of Y at t = TOUT.
! ISTATE = Integer flag (input and output).  Set ISTATE = 1.
! IOPT   = 0 to indicate no optional input used.
! RWORK  = Real work array of length at least:
!             20 + 16*NEQ                      for MF = 10,
!             22 +  9*NEQ + 2*NEQ**2           for MF = 21 or 22,
!             22 + 11*NEQ + (3*ML + 2*MU)*NEQ  for MF = 24 or 25.
! LRW    = Declared length of RWORK (in user's DIMENSION statement).
! IWORK  = Integer work array of length at least:
!             30        for MF = 10,
!             30 + NEQ  for MF = 21, 22, 24, or 25.
!          If MF = 24 or 25, input in IWORK(1),IWORK(2) the lower
!          and upper half-bandwidths ML,MU.
! LIW    = Declared length of IWORK (in user's DIMENSION statement).
! JAC    = Name of subroutine for Jacobian matrix (MF = 21 or 24).
!          If used, this name must be declared external in calling
!          program.  If not used, pass a dummy name.
! MF     = Method flag.  Standard values are:
!          10 for nonstiff (Adams) method, no Jacobian used.
!          21 for stiff (BDF) method, user-supplied full Jacobian.
!          22 for stiff method, internally generated full Jacobian.
!          24 for stiff method, user-supplied banded Jacobian.
!          25 for stiff method, internally generated banded Jacobian.
! RPAR,IPAR = user-defined real and integer arrays passed to F and JAC.
! Note that the main program must declare arrays Y, RWORK, IWORK,
! and possibly ATOL, RPAR, and IPAR.
!
! E. The output from the first call (or any call) is:
!      Y = Array of computed values of y(t) vector.
!      T = Corresponding value of independent variable (normally TOUT).
! ISTATE = 2  if DVODE was successful, negative otherwise.
!          -1 means excess work done on this call. (Perhaps wrong MF.)
!          -2 means excess accuracy requested. (Tolerances too small.)
!          -3 means illegal input detected. (See printed message.)
!          -4 means repeated error test failures. (Check all input.)
!          -5 means repeated convergence failures. (Perhaps bad
!             Jacobian supplied or wrong choice of MF or tolerances.)
!          -6 means error weight became zero during problem. (Solution
!             component i vanished, and ATOL or ATOL(i) = 0.)
!
! F. To continue the integration after a successful return, simply
! reset TOUT and call DVODE again.  No other parameters need be reset.
!
!-----------------------------------------------------------------------
! EXAMPLE PROBLEM
!
! The following is a simple example problem, with the coding
! needed for its solution by DVODE.  The problem is from chemical
! kinetics, and consists of the following three rate equations:
!     dy1/dt = -.04*y1 + 1.e4*y2*y3
!     dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2
!     dy3/dt = 3.e7*y2**2
! on the interval from t = 0.0 to t = 4.e10, with initial conditions
! y1 = 1.0, y2 = y3 = 0.  The problem is stiff.
!
! The following coding solves this problem with DVODE, using MF = 21
! and printing results at t = .4, 4., ..., 4.e10.  It uses
! ITOL = 2 and ATOL much smaller for y2 than y1 or y3 because
! y2 has much smaller values.
! At the end of the run, statistical quantities of interest are
! printed. (See optional output in the full description below.)
! To generate Fortran source code, replace C in column 1 with a blank
! in the coding below.
!
!     EXTERNAL FEX, JEX
!     DOUBLE PRECISION ATOL, RPAR, RTOL, RWORK, T, TOUT, Y
!     DIMENSION Y(3), ATOL(3), RWORK(67), IWORK(33)
!     NEQ = 3
!     Y(1) = 1.0D0
!     Y(2) = 0.0D0
!     Y(3) = 0.0D0
!     T = 0.0D0
!     TOUT = 0.4D0
!     ITOL = 2
!     RTOL = 1.D-4
!     ATOL(1) = 1.D-8
!     ATOL(2) = 1.D-14
!     ATOL(3) = 1.D-6
!     ITASK = 1
!     ISTATE = 1
!     IOPT = 0
!     LRW = 67
!     LIW = 33
!     MF = 21
!     DO 40 IOUT = 1,12
!       CALL DVODE(FEX,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,
!    1            IOPT,RWORK,LRW,IWORK,LIW,JEX,MF,RPAR,IPAR)
!       WRITE(6,20)T,Y(1),Y(2),Y(3)
! 20    FORMAT(' At t =',D12.4,'   y =',3D14.6)
!       IF (ISTATE .LT. 0) GO TO 80
! 40    TOUT = TOUT*10.
!     WRITE(6,60) IWORK(11),IWORK(12),IWORK(13),IWORK(19),
!    1            IWORK(20),IWORK(21),IWORK(22)
! 60  FORMAT(/' No. steps =',I4,'   No. f-s =',I4,
!    1       '   No. J-s =',I4,'   No. LU-s =',I4/
!    2       '  No. nonlinear iterations =',I4/
!    3       '  No. nonlinear convergence failures =',I4/
!    4       '  No. error test failures =',I4/)
!     STOP
! 80  WRITE(6,90)ISTATE
! 90  FORMAT(///' Error halt: ISTATE =',I3)
!     STOP
!     END
!
!     SUBROUTINE FEX (NEQ, T, Y, YDOT, RPAR, IPAR)
!     DOUBLE PRECISION RPAR, T, Y, YDOT
!     DIMENSION Y(NEQ), YDOT(NEQ)
!     YDOT(1) = -.04D0*Y(1) + 1.D4*Y(2)*Y(3)
!     YDOT(3) = 3.D7*Y(2)*Y(2)
!     YDOT(2) = -YDOT(1) - YDOT(3)
!     RETURN
!     END
!
!     SUBROUTINE JEX (NEQ, T, Y, ML, MU, PD, NRPD, RPAR, IPAR)
!     DOUBLE PRECISION PD, RPAR, T, Y
!     DIMENSION Y(NEQ), PD(NRPD,NEQ)
!     PD(1,1) = -.04D0
!     PD(1,2) = 1.D4*Y(3)
!     PD(1,3) = 1.D4*Y(2)
!     PD(2,1) = .04D0
!     PD(2,3) = -PD(1,3)
!     PD(3,2) = 6.D7*Y(2)
!     PD(2,2) = -PD(1,2) - PD(3,2)
!     RETURN
!     END
!
! The following output was obtained from the above program on a
! Cray-1 computer with the CFT compiler.
!
! At t =  4.0000e-01   y =  9.851680e-01  3.386314e-05  1.479817e-02
! At t =  4.0000e+00   y =  9.055255e-01  2.240539e-05  9.445214e-02
! At t =  4.0000e+01   y =  7.158108e-01  9.184883e-06  2.841800e-01
! At t =  4.0000e+02   y =  4.505032e-01  3.222940e-06  5.494936e-01
! At t =  4.0000e+03   y =  1.832053e-01  8.942690e-07  8.167938e-01
! At t =  4.0000e+04   y =  3.898560e-02  1.621875e-07  9.610142e-01
! At t =  4.0000e+05   y =  4.935882e-03  1.984013e-08  9.950641e-01
! At t =  4.0000e+06   y =  5.166183e-04  2.067528e-09  9.994834e-01
! At t =  4.0000e+07   y =  5.201214e-05  2.080593e-10  9.999480e-01
! At t =  4.0000e+08   y =  5.213149e-06  2.085271e-11  9.999948e-01
! At t =  4.0000e+09   y =  5.183495e-07  2.073399e-12  9.999995e-01
! At t =  4.0000e+10   y =  5.450996e-08  2.180399e-13  9.999999e-01
!
! No. steps = 595   No. f-s = 832   No. J-s =  13   No. LU-s = 112
!  No. nonlinear iterations = 831
!  No. nonlinear convergence failures =   0
!  No. error test failures =  22
!-----------------------------------------------------------------------
! Full description of user interface to DVODE.
!
! The user interface to DVODE consists of the following parts.
!
! i.   The call sequence to subroutine DVODE, which is a driver
!      routine for the solver.  This includes descriptions of both
!      the call sequence arguments and of user-supplied routines.
!      Following these descriptions is
!        * a description of optional input available through the
!          call sequence,
!        * a description of optional output (in the work arrays), and
!        * instructions for interrupting and restarting a solution.
!
! ii.  Descriptions of other routines in the DVODE package that may be
!      (optionally) called by the user.  These provide the ability to
!      alter error message handling, save and restore the internal
!      COMMON, and obtain specified derivatives of the solution y(t).
!
! iii. Descriptions of COMMON blocks to be declared in overlay
!      or similar environments.
!
! iv.  Description of two routines in the DVODE package, either of
!      which the user may replace with his own version, if desired.
!      these relate to the measurement of errors.
!
!-----------------------------------------------------------------------
! Part i.  Call Sequence.
!
! The call sequence parameters used for input only are
!     F, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK, IOPT, LRW, LIW, JAC, MF,
! and those used for both input and output are
!     Y, T, ISTATE.
! The work arrays RWORK and IWORK are also used for conditional and
! optional input and optional output.  (The term output here refers
! to the return from subroutine DVODE to the user's calling program.)
!
! The legality of input parameters will be thoroughly checked on the
! initial call for the problem, but not checked thereafter unless a
! change in input parameters is flagged by ISTATE = 3 in the input.
!
! The descriptions of the call arguments are as follows.
!
! F      = The name of the user-supplied subroutine defining the
!          ODE system.  The system must be put in the first-order
!          form dy/dt = f(t,y), where f is a vector-valued function
!          of the scalar t and the vector y.  Subroutine F is to
!          compute the function f.  It is to have the form
!               SUBROUTINE F (NEQ, T, Y, YDOT, RPAR, IPAR)
!               DOUBLE PRECISION T, Y(NEQ), YDOT(NEQ), RPAR
!          where NEQ, T, and Y are input, and the array YDOT = f(t,y)
!          is output.  Y and YDOT are arrays of length NEQ.
!          Subroutine F should not alter Y(1),...,Y(NEQ).
!          F must be declared EXTERNAL in the calling program.
!
!          Subroutine F may access user-defined real and integer
!          work arrays RPAR and IPAR, which are to be dimensioned
!          in the main program.
!
!          If quantities computed in the F routine are needed
!          externally to DVODE, an extra call to F should be made
!          for this purpose, for consistent and accurate results.
!          If only the derivative dy/dt is needed, use DVINDY instead.
!
! NEQ    = The size of the ODE system (number of first order
!          ordinary differential equations).  Used only for input.
!          NEQ may not be increased during the problem, but
!          can be decreased (with ISTATE = 3 in the input).
!
! Y      = A real array for the vector of dependent variables, of
!          length NEQ or more.  Used for both input and output on the
!          first call (ISTATE = 1), and only for output on other calls.
!          On the first call, Y must contain the vector of initial
!          values.  In the output, Y contains the computed solution
!          evaluated at T.  If desired, the Y array may be used
!          for other purposes between calls to the solver.
!
!          This array is passed as the Y argument in all calls to
!          F and JAC.
!
! T      = The independent variable.  In the input, T is used only on
!          the first call, as the initial point of the integration.
!          In the output, after each call, T is the value at which a
!          computed solution Y is evaluated (usually the same as TOUT).
!          On an error return, T is the farthest point reached.
!
! TOUT   = The next value of t at which a computed solution is desired.
!          Used only for input.
!
!          When starting the problem (ISTATE = 1), TOUT may be equal
!          to T for one call, then should .ne. T for the next call.
!          For the initial T, an input value of TOUT .ne. T is used
!          in order to determine the direction of the integration
!          (i.e. the algebraic sign of the step sizes) and the rough
!          scale of the problem.  Integration in either direction
!          (forward or backward in t) is permitted.
!
!          If ITASK = 2 or 5 (one-step modes), TOUT is ignored after
!          the first call (i.e. the first call with TOUT .ne. T).
!          Otherwise, TOUT is required on every call.
!
!          If ITASK = 1, 3, or 4, the values of TOUT need not be
!          monotone, but a value of TOUT which backs up is limited
!          to the current internal t interval, whose endpoints are
!          TCUR - HU and TCUR.  (See optional output, below, for
!          TCUR and HU.)
!
! ITOL   = An indicator for the type of error control.  See
!          description below under ATOL.  Used only for input.
!
! RTOL   = A relative error tolerance parameter, either a scalar or
!          an array of length NEQ.  See description below under ATOL.
!          Input only.
!
! ATOL   = An absolute error tolerance parameter, either a scalar or
!          an array of length NEQ.  Input only.
!
!          The input parameters ITOL, RTOL, and ATOL determine
!          the error control performed by the solver.  The solver will
!          control the vector e = (e(i)) of estimated local errors
!          in Y, according to an inequality of the form
!                      rms-norm of ( e(i)/EWT(i) )   .le.   1,
!          where       EWT(i) = RTOL(i)*abs(Y(i)) + ATOL(i),
!          and the rms-norm (root-mean-square norm) here is
!          rms-norm(v) = sqrt(sum v(i)**2 / NEQ).  Here EWT = (EWT(i))
!          is a vector of weights which must always be positive, and
!          the values of RTOL and ATOL should all be non-negative.
!          The following table gives the types (scalar/array) of
!          RTOL and ATOL, and the corresponding form of EWT(i).
!
!             ITOL    RTOL       ATOL          EWT(i)
!              1     scalar     scalar     RTOL*ABS(Y(i)) + ATOL
!              2     scalar     array      RTOL*ABS(Y(i)) + ATOL(i)
!              3     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL
!              4     array      array      RTOL(i)*ABS(Y(i)) + ATOL(i)
!
!          When either of these parameters is a scalar, it need not
!          be dimensioned in the user's calling program.
!
!          If none of the above choices (with ITOL, RTOL, and ATOL
!          fixed throughout the problem) is suitable, more general
!          error controls can be obtained by substituting
!          user-supplied routines for the setting of EWT and/or for
!          the norm calculation.  See Part iv below.
!
!          If global errors are to be estimated by making a repeated
!          run on the same problem with smaller tolerances, then all
!          components of RTOL and ATOL (i.e. of EWT) should be scaled
!          down uniformly.
!
! ITASK  = An index specifying the task to be performed.
!          Input only.  ITASK has the following values and meanings.
!          1  means normal computation of output values of y(t) at
!             t = TOUT (by overshooting and interpolating).
!          2  means take one step only and return.
!          3  means stop at the first internal mesh point at or
!             beyond t = TOUT and return.
!          4  means normal computation of output values of y(t) at
!             t = TOUT but without overshooting t = TCRIT.
!             TCRIT must be input as RWORK(1).  TCRIT may be equal to
!             or beyond TOUT, but not behind it in the direction of
!             integration.  This option is useful if the problem
!             has a singularity at or beyond t = TCRIT.
!          5  means take one step, without passing TCRIT, and return.
!             TCRIT must be input as RWORK(1).
!
!          Note:  If ITASK = 4 or 5 and the solver reaches TCRIT
!          (within roundoff), it will return T = TCRIT (exactly) to
!          indicate this (unless ITASK = 4 and TOUT comes before TCRIT,
!          in which case answers at T = TOUT are returned first).
!
! ISTATE = an index used for input and output to specify the
!          the state of the calculation.
!
!          In the input, the values of ISTATE are as follows.
!          1  means this is the first call for the problem
!             (initializations will be done).  See note below.
!          2  means this is not the first call, and the calculation
!             is to continue normally, with no change in any input
!             parameters except possibly TOUT and ITASK.
!             (If ITOL, RTOL, and/or ATOL are changed between calls
!             with ISTATE = 2, the new values will be used but not
!             tested for legality.)
!          3  means this is not the first call, and the
!             calculation is to continue normally, but with
!             a change in input parameters other than
!             TOUT and ITASK.  Changes are allowed in
!             NEQ, ITOL, RTOL, ATOL, IOPT, LRW, LIW, MF, ML, MU,
!             and any of the optional input except H0.
!             (See IWORK description for ML and MU.)
!          Note:  A preliminary call with TOUT = T is not counted
!          as a first call here, as no initialization or checking of
!          input is done.  (Such a call is sometimes useful to include
!          the initial conditions in the output.)
!          Thus the first call for which TOUT .ne. T requires
!          ISTATE = 1 in the input.
!
!          In the output, ISTATE has the following values and meanings.
!           1  means nothing was done, as TOUT was equal to T with
!              ISTATE = 1 in the input.
!           2  means the integration was performed successfully.
!          -1  means an excessive amount of work (more than MXSTEP
!              steps) was done on this call, before completing the
!              requested task, but the integration was otherwise
!              successful as far as T.  (MXSTEP is an optional input
!              and is normally 500.)  To continue, the user may
!              simply reset ISTATE to a value .gt. 1 and call again.
!              (The excess work step counter will be reset to 0.)
!              In addition, the user may increase MXSTEP to avoid
!              this error return.  (See optional input below.)
!          -2  means too much accuracy was requested for the precision
!              of the machine being used.  This was detected before
!              completing the requested task, but the integration
!              was successful as far as T.  To continue, the tolerance
!              parameters must be reset, and ISTATE must be set
!              to 3.  The optional output TOLSF may be used for this
!              purpose.  (Note: If this condition is detected before
!              taking any steps, then an illegal input return
!              (ISTATE = -3) occurs instead.)
!          -3  means illegal input was detected, before taking any
!              integration steps.  See written message for details.
!              Note:  If the solver detects an infinite loop of calls
!              to the solver with illegal input, it will cause
!              the run to stop.
!          -4  means there were repeated error test failures on
!              one attempted step, before completing the requested
!              task, but the integration was successful as far as T.
!              The problem may have a singularity, or the input
!              may be inappropriate.
!          -5  means there were repeated convergence test failures on
!              one attempted step, before completing the requested
!              task, but the integration was successful as far as T.
!              This may be caused by an inaccurate Jacobian matrix,
!              if one is being used.
!          -6  means EWT(i) became zero for some i during the
!              integration.  Pure relative error control (ATOL(i)=0.0)
!              was requested on a variable which has now vanished.
!              The integration was successful as far as T.
!
!          Note:  Since the normal output value of ISTATE is 2,
!          it does not need to be reset for normal continuation.
!          Also, since a negative input value of ISTATE will be
!          regarded as illegal, a negative output value requires the
!          user to change it, and possibly other input, before
!          calling the solver again.
!
! IOPT   = An integer flag to specify whether or not any optional
!          input is being used on this call.  Input only.
!          The optional input is listed separately below.
!          IOPT = 0 means no optional input is being used.
!                   Default values will be used in all cases.
!          IOPT = 1 means optional input is being used.
!
! RWORK  = A real working array (double precision).
!          The length of RWORK must be at least
!             20 + NYH*(MAXORD + 1) + 3*NEQ + LWM    where
!          NYH    = the initial value of NEQ,
!          MAXORD = 12 (if METH = 1) or 5 (if METH = 2) (unless a
!                   smaller value is given as an optional input),
!          LWM = length of work space for matrix-related data:
!          LWM = 0             if MITER = 0,
!          LWM = 2*NEQ**2 + 2  if MITER = 1 or 2, and MF.gt.0,
!          LWM = NEQ**2 + 2    if MITER = 1 or 2, and MF.lt.0,
!          LWM = NEQ + 2       if MITER = 3,
!          LWM = (3*ML+2*MU+2)*NEQ + 2 if MITER = 4 or 5, and MF.gt.0,
!          LWM = (2*ML+MU+1)*NEQ + 2   if MITER = 4 or 5, and MF.lt.0.
!          (See the MF description for METH and MITER.)
!          Thus if MAXORD has its default value and NEQ is constant,
!          this length is:
!             20 + 16*NEQ                    for MF = 10,
!             22 + 16*NEQ + 2*NEQ**2         for MF = 11 or 12,
!             22 + 16*NEQ + NEQ**2           for MF = -11 or -12,
!             22 + 17*NEQ                    for MF = 13,
!             22 + 18*NEQ + (3*ML+2*MU)*NEQ  for MF = 14 or 15,
!             22 + 17*NEQ + (2*ML+MU)*NEQ    for MF = -14 or -15,
!             20 +  9*NEQ                    for MF = 20,
!             22 +  9*NEQ + 2*NEQ**2         for MF = 21 or 22,
!             22 +  9*NEQ + NEQ**2           for MF = -21 or -22,
!             22 + 10*NEQ                    for MF = 23,
!             22 + 11*NEQ + (3*ML+2*MU)*NEQ  for MF = 24 or 25.
!             22 + 10*NEQ + (2*ML+MU)*NEQ    for MF = -24 or -25.
!          The first 20 words of RWORK are reserved for conditional
!          and optional input and optional output.
!
!          The following word in RWORK is a conditional input:
!            RWORK(1) = TCRIT = critical value of t which the solver
!                       is not to overshoot.  Required if ITASK is
!                       4 or 5, and ignored otherwise.  (See ITASK.)
!
! LRW    = The length of the array RWORK, as declared by the user.
!          (This will be checked by the solver.)
!
! IWORK  = An integer work array.  The length of IWORK must be at least
!             30        if MITER = 0 or 3 (MF = 10, 13, 20, 23), or
!             30 + NEQ  otherwise (abs(MF) = 11,12,14,15,21,22,24,25).
!          The first 30 words of IWORK are reserved for conditional and
!          optional input and optional output.
!
!          The following 2 words in IWORK are conditional input:
!            IWORK(1) = ML     These are the lower and upper
!            IWORK(2) = MU     half-bandwidths, respectively, of the
!                       banded Jacobian, excluding the main diagonal.
!                       The band is defined by the matrix locations
!                       (i,j) with i-ML .le. j .le. i+MU.  ML and MU
!                       must satisfy  0 .le.  ML,MU  .le. NEQ-1.
!                       These are required if MITER is 4 or 5, and
!                       ignored otherwise.  ML and MU may in fact be
!                       the band parameters for a matrix to which
!                       df/dy is only approximately equal.
!
! LIW    = the length of the array IWORK, as declared by the user.
!          (This will be checked by the solver.)
!
! Note:  The work arrays must not be altered between calls to DVODE
! for the same problem, except possibly for the conditional and
! optional input, and except for the last 3*NEQ words of RWORK.
! The latter space is used for internal scratch space, and so is
! available for use by the user outside DVODE between calls, if
! desired (but not for use by F or JAC).
!
! JAC    = The name of the user-supplied routine (MITER = 1 or 4) to
!          compute the Jacobian matrix, df/dy, as a function of
!          the scalar t and the vector y.  It is to have the form
!               SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD,
!                               RPAR, IPAR)
!               DOUBLE PRECISION T, Y(NEQ), PD(NROWPD,NEQ), RPAR
!          where NEQ, T, Y, ML, MU, and NROWPD are input and the array
!          PD is to be loaded with partial derivatives (elements of the
!          Jacobian matrix) in the output.  PD must be given a first
!          dimension of NROWPD.  T and Y have the same meaning as in
!          Subroutine F.
!               In the full matrix case (MITER = 1), ML and MU are
!          ignored, and the Jacobian is to be loaded into PD in
!          columnwise manner, with df(i)/dy(j) loaded into PD(i,j).
!               In the band matrix case (MITER = 4), the elements
!          within the band are to be loaded into PD in columnwise
!          manner, with diagonal lines of df/dy loaded into the rows
!          of PD. Thus df(i)/dy(j) is to be loaded into PD(i-j+MU+1,j).
!          ML and MU are the half-bandwidth parameters. (See IWORK).
!          The locations in PD in the two triangular areas which
!          correspond to nonexistent matrix elements can be ignored
!          or loaded arbitrarily, as they are overwritten by DVODE.
!               JAC need not provide df/dy exactly.  A crude
!          approximation (possibly with a smaller bandwidth) will do.
!               In either case, PD is preset to zero by the solver,
!          so that only the nonzero elements need be loaded by JAC.
!          Each call to JAC is preceded by a call to F with the same
!          arguments NEQ, T, and Y.  Thus to gain some efficiency,
!          intermediate quantities shared by both calculations may be
!          saved in a user COMMON block by F and not recomputed by JAC,
!          if desired.  Also, JAC may alter the Y array, if desired.
!          JAC must be declared external in the calling program.
!               Subroutine JAC may access user-defined real and integer
!          work arrays, RPAR and IPAR, whose dimensions are set by the
!          user in the main program.
!
! MF     = The method flag.  Used only for input.  The legal values of
!          MF are 10, 11, 12, 13, 14, 15, 20, 21, 22, 23, 24, 25,
!          -11, -12, -14, -15, -21, -22, -24, -25.
!          MF is a signed two-digit integer, MF = JSV*(10*METH + MITER).
!          JSV = SIGN(MF) indicates the Jacobian-saving strategy:
!            JSV =  1 means a copy of the Jacobian is saved for reuse
!                     in the corrector iteration algorithm.
!            JSV = -1 means a copy of the Jacobian is not saved
!                     (valid only for MITER = 1, 2, 4, or 5).
!          METH indicates the basic linear multistep method:
!            METH = 1 means the implicit Adams method.
!            METH = 2 means the method based on backward
!                     differentiation formulas (BDF-s).
!          MITER indicates the corrector iteration method:
!            MITER = 0 means functional iteration (no Jacobian matrix
!                      is involved).
!            MITER = 1 means chord iteration with a user-supplied
!                      full (NEQ by NEQ) Jacobian.
!            MITER = 2 means chord iteration with an internally
!                      generated (difference quotient) full Jacobian
!                      (using NEQ extra calls to F per df/dy value).
!            MITER = 3 means chord iteration with an internally
!                      generated diagonal Jacobian approximation
!                      (using 1 extra call to F per df/dy evaluation).
!            MITER = 4 means chord iteration with a user-supplied
!                      banded Jacobian.
!            MITER = 5 means chord iteration with an internally
!                      generated banded Jacobian (using ML+MU+1 extra
!                      calls to F per df/dy evaluation).
!          If MITER = 1 or 4, the user must supply a subroutine JAC
!          (the name is arbitrary) as described above under JAC.
!          For other values of MITER, a dummy argument can be used.
!
! RPAR     User-specified array used to communicate real parameters
!          to user-supplied subroutines.  If RPAR is a vector, then
!          it must be dimensioned in the user's main program.  If it
!          is unused or it is a scalar, then it need not be
!          dimensioned.
!
! IPAR     User-specified array used to communicate integer parameter
!          to user-supplied subroutines.  The comments on dimensioning
!          RPAR apply to IPAR.
!-----------------------------------------------------------------------
! Optional Input.
!
! The following is a list of the optional input provided for in the
! call sequence.  (See also Part ii.)  For each such input variable,
! this table lists its name as used in this documentation, its
! location in the call sequence, its meaning, and the default value.
! The use of any of this input requires IOPT = 1, and in that
! case all of this input is examined.  A value of zero for any
! of these optional input variables will cause the default value to be
! used.  Thus to use a subset of the optional input, simply preload
! locations 5 to 10 in RWORK and IWORK to 0.0 and 0 respectively, and
! then set those of interest to nonzero values.
!
! NAME    LOCATION      MEANING AND DEFAULT VALUE
!
! H0      RWORK(5)  The step size to be attempted on the first step.
!                   The default value is determined by the solver.
!
! HMAX    RWORK(6)  The maximum absolute step size allowed.
!                   The default value is infinite.
!
! HMIN    RWORK(7)  The minimum absolute step size allowed.
!                   The default value is 0.  (This lower bound is not
!                   enforced on the final step before reaching TCRIT
!                   when ITASK = 4 or 5.)
!
! MAXORD  IWORK(5)  The maximum order to be allowed.  The default
!                   value is 12 if METH = 1, and 5 if METH = 2.
!                   If MAXORD exceeds the default value, it will
!                   be reduced to the default value.
!                   If MAXORD is changed during the problem, it may
!                   cause the current order to be reduced.
!
! MXSTEP  IWORK(6)  Maximum number of (internally defined) steps
!                   allowed during one call to the solver.
!                   The default value is 500.
!
! MXHNIL  IWORK(7)  Maximum number of messages printed (per problem)
!                   warning that T + H = T on a step (H = step size).
!                   This must be positive to result in a non-default
!                   value.  The default value is 10.
!
!-----------------------------------------------------------------------
! Optional Output.
!
! As optional additional output from DVODE, the variables listed
! below are quantities related to the performance of DVODE
! which are available to the user.  These are communicated by way of
! the work arrays, but also have internal mnemonic names as shown.
! Except where stated otherwise, all of this output is defined
! on any successful return from DVODE, and on any return with
! ISTATE = -1, -2, -4, -5, or -6.  On an illegal input return
! (ISTATE = -3), they will be unchanged from their existing values
! (if any), except possibly for TOLSF, LENRW, and LENIW.
! On any error return, output relevant to the error will be defined,
! as noted below.
!
! NAME    LOCATION      MEANING
!
! HU      RWORK(11) The step size in t last used (successfully).
!
! HCUR    RWORK(12) The step size to be attempted on the next step.
!
! TCUR    RWORK(13) The current value of the independent variable
!                   which the solver has actually reached, i.e. the
!                   current internal mesh point in t.  In the output,
!                   TCUR will always be at least as far from the
!                   initial value of t as the current argument T,
!                   but may be farther (if interpolation was done).
!
! TOLSF   RWORK(14) A tolerance scale factor, greater than 1.0,
!                   computed when a request for too much accuracy was
!                   detected (ISTATE = -3 if detected at the start of
!                   the problem, ISTATE = -2 otherwise).  If ITOL is
!                   left unaltered but RTOL and ATOL are uniformly
!                   scaled up by a factor of TOLSF for the next call,
!                   then the solver is deemed likely to succeed.
!                   (The user may also ignore TOLSF and alter the
!                   tolerance parameters in any other way appropriate.)
!
! NST     IWORK(11) The number of steps taken for the problem so far.
!
! NFE     IWORK(12) The number of f evaluations for the problem so far.
!
! NJE     IWORK(13) The number of Jacobian evaluations so far.
!
! NQU     IWORK(14) The method order last used (successfully).
!
! NQCUR   IWORK(15) The order to be attempted on the next step.
!
! IMXER   IWORK(16) The index of the component of largest magnitude in
!                   the weighted local error vector ( e(i)/EWT(i) ),
!                   on an error return with ISTATE = -4 or -5.
!
! LENRW   IWORK(17) The length of RWORK actually required.
!                   This is defined on normal returns and on an illegal
!                   input return for insufficient storage.
!
! LENIW   IWORK(18) The length of IWORK actually required.
!                   This is defined on normal returns and on an illegal
!                   input return for insufficient storage.
!
! NLU     IWORK(19) The number of matrix LU decompositions so far.
!
! NNI     IWORK(20) The number of nonlinear (Newton) iterations so far.
!
! NCFN    IWORK(21) The number of convergence failures of the nonlinear
!                   solver so far.
!
! NETF    IWORK(22) The number of error test failures of the integrator
!                   so far.
!
! The following two arrays are segments of the RWORK array which
! may also be of interest to the user as optional output.
! For each array, the table below gives its internal name,
! its base address in RWORK, and its description.
!
! NAME    BASE ADDRESS      DESCRIPTION
!
! YH      21             The Nordsieck history array, of size NYH by
!                        (NQCUR + 1), where NYH is the initial value
!                        of NEQ.  For j = 0,1,...,NQCUR, column j+1
!                        of YH contains HCUR**j/factorial(j) times
!                        the j-th derivative of the interpolating
!                        polynomial currently representing the
!                        solution, evaluated at t = TCUR.
!
! ACOR     LENRW-NEQ+1   Array of size NEQ used for the accumulated
!                        corrections on each step, scaled in the output
!                        to represent the estimated local error in Y
!                        on the last step.  This is the vector e in
!                        the description of the error control.  It is
!                        defined only on a successful return from DVODE.
!
!-----------------------------------------------------------------------
! Interrupting and Restarting
!
! If the integration of a given problem by DVODE is to be
! interrrupted and then later continued, such as when restarting
! an interrupted run or alternating between two or more ODE problems,
! the user should save, following the return from the last DVODE call
! prior to the interruption, the contents of the call sequence
! variables and internal COMMON blocks, and later restore these
! values before the next DVODE call for that problem.  To save
! and restore the COMMON blocks, use subroutine DVSRCO, as
! described below in part ii.
!
! In addition, if non-default values for either LUN or MFLAG are
! desired, an extra call to XSETUN and/or XSETF should be made just
! before continuing the integration.  See Part ii below for details.
!
!-----------------------------------------------------------------------
! Part ii.  Other Routines Callable.
!
! The following are optional calls which the user may make to
! gain additional capabilities in conjunction with DVODE.
! (The routines XSETUN and XSETF are designed to conform to the
! SLATEC error handling package.)
!
!     FORM OF CALL                  FUNCTION
!  CALL XSETUN(LUN)           Set the logical unit number, LUN, for
!                             output of messages from DVODE, if
!                             the default is not desired.
!                             The default value of LUN is 6.
!
!  CALL XSETF(MFLAG)          Set a flag to control the printing of
!                             messages by DVODE.
!                             MFLAG = 0 means do not print. (Danger:
!                             This risks losing valuable information.)
!                             MFLAG = 1 means print (the default).
!
!                             Either of the above calls may be made at
!                             any time and will take effect immediately.
!
!  CALL DVSRCO(RSAV,ISAV,JOB) Saves and restores the contents of
!                             the internal COMMON blocks used by
!                             DVODE. (See Part iii below.)
!                             RSAV must be a real array of length 49
!                             or more, and ISAV must be an integer
!                             array of length 40 or more.
!                             JOB=1 means save COMMON into RSAV/ISAV.
!                             JOB=2 means restore COMMON from RSAV/ISAV.
!                                DVSRCO is useful if one is
!                             interrupting a run and restarting
!                             later, or alternating between two or
!                             more problems solved with DVODE.
!
!  CALL DVINDY(,,,,,)         Provide derivatives of y, of various
!        (See below.)         orders, at a specified point T, if
!                             desired.  It may be called only after
!                             a successful return from DVODE.
!
! The detailed instructions for using DVINDY are as follows.
! The form of the call is:
!
!  CALL DVINDY (T, K, RWORK(21), NYH, DKY, IFLAG)
!
! The input parameters are:
!
! T         = Value of independent variable where answers are desired
!             (normally the same as the T last returned by DVODE).
!             For valid results, T must lie between TCUR - HU and TCUR.
!             (See optional output for TCUR and HU.)
! K         = Integer order of the derivative desired.  K must satisfy
!             0 .le. K .le. NQCUR, where NQCUR is the current order
!             (see optional output).  The capability corresponding
!             to K = 0, i.e. computing y(T), is already provided
!             by DVODE directly.  Since NQCUR .ge. 1, the first
!             derivative dy/dt is always available with DVINDY.
! RWORK(21) = The base address of the history array YH.
! NYH       = Column length of YH, equal to the initial value of NEQ.
!
! The output parameters are:
!
! DKY       = A real array of length NEQ containing the computed value
!             of the K-th derivative of y(t).
! IFLAG     = Integer flag, returned as 0 if K and T were legal,
!             -1 if K was illegal, and -2 if T was illegal.
!             On an error return, a message is also written.
!-----------------------------------------------------------------------
! Part iii.  COMMON Blocks.
! If DVODE is to be used in an overlay situation, the user
! must declare, in the primary overlay, the variables in:
!   (1) the call sequence to DVODE,
!   (2) the two internal COMMON blocks
!         /DVOD01/  of length  81  (48 double precision words
!                         followed by 33 integer words),
!         /DVOD02/  of length  9  (1 double precision word
!                         followed by 8 integer words),
!
! If DVODE is used on a system in which the contents of internal
! COMMON blocks are not preserved between calls, the user should
! declare the above two COMMON blocks in his main program to insure
! that their contents are preserved.
!
!-----------------------------------------------------------------------
! Part iv.  Optionally Replaceable Solver Routines.
!
! Below are descriptions of two routines in the DVODE package which
! relate to the measurement of errors.  Either routine can be
! replaced by a user-supplied version, if desired.  However, since such
! a replacement may have a major impact on performance, it should be
! done only when absolutely necessary, and only with great caution.
! (Note: The means by which the package version of a routine is
! superseded by the user's version may be system-dependent.)
!
! (a) DEWSET.
! The following subroutine is called just before each internal
! integration step, and sets the array of error weights, EWT, as
! described under ITOL/RTOL/ATOL above:
!     SUBROUTINE DEWSET (NEQ, ITOL, RTOL, ATOL, YCUR, EWT)
! where NEQ, ITOL, RTOL, and ATOL are as in the DVODE call sequence,
! YCUR contains the current dependent variable vector, and
! EWT is the array of weights set by DEWSET.
!
! If the user supplies this subroutine, it must return in EWT(i)
! (i = 1,...,NEQ) a positive quantity suitable for comparison with
! errors in Y(i).  The EWT array returned by DEWSET is passed to the
! DVNORM routine (See below.), and also used by DVODE in the computation
! of the optional output IMXER, the diagonal Jacobian approximation,
! and the increments for difference quotient Jacobians.
!
! In the user-supplied version of DEWSET, it may be desirable to use
! the current values of derivatives of y.  Derivatives up to order NQ
! are available from the history array YH, described above under
! Optional Output.  In DEWSET, YH is identical to the YCUR array,
! extended to NQ + 1 columns with a column length of NYH and scale
! factors of h**j/factorial(j).  On the first call for the problem,
! given by NST = 0, NQ is 1 and H is temporarily set to 1.0.
! NYH is the initial value of NEQ.  The quantities NQ, H, and NST
! can be obtained by including in DEWSET the statements:
!     DOUBLE PRECISION RVOD, H, HU
!     COMMON /DVOD01/ RVOD(48), IVOD(33)
!     COMMON /DVOD02/ HU, NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
!     NQ = IVOD(28)
!     H = RVOD(21)
! Thus, for example, the current value of dy/dt can be obtained as
! YCUR(NYH+i)/H  (i=1,...,NEQ)  (and the division by H is
! unnecessary when NST = 0).
!
! (b) DVNORM.
! The following is a real function routine which computes the weighted
! root-mean-square norm of a vector v:
!     D = DVNORM (N, V, W)
! where:
!   N = the length of the vector,
!   V = real array of length N containing the vector,
!   W = real array of length N containing weights,
!   D = sqrt( (1/N) * sum(V(i)*W(i))**2 ).
! DVNORM is called with N = NEQ and with W(i) = 1.0/EWT(i), where
! EWT is as set by subroutine DEWSET.
!
! If the user supplies this function, it should return a non-negative
! value of DVNORM suitable for use in the error control in DVODE.
! None of the arguments should be altered by DVNORM.
! For example, a user-supplied DVNORM routine might:
!   -substitute a max-norm of (V(i)*W(i)) for the rms-norm, or
!   -ignore some components of V in the norm, with the effect of
!    suppressing the error control on those components of Y.
!-----------------------------------------------------------------------
! REVISION HISTORY (YYYYMMDD)
!  19890615  Date Written.  Initial release.
!  19890922  Added interrupt/restart ability, minor changes throughout.
!  19910228  Minor revisions in line format,  prologue, etc.
!  19920227  Modifications by D. Pang:
!            (1) Applied subgennam to get generic intrinsic names.
!            (2) Changed intrinsic names to generic in comments.
!            (3) Added *DECK lines before each routine.
!  19920721  Names of routines and labeled Common blocks changed, so as
!            to be unique in combined single/double precision code (ACH).
!  19920722  Minor revisions to prologue (ACH).
!  19920831  Conversion to double precision done (ACH).
!  19921106  Fixed minor bug: ETAQ,ETAQM1 in DVSTEP SAVE statement (ACH).
!  19921118  Changed LUNSAV/MFLGSV to IXSAV (ACH).
!  19941222  Removed MF overwrite; attached sign to H in estimated second
!            deriv. in DVHIN; misc. comment changes throughout (ACH).
!  19970515  Minor corrections to comments in prologue, DVJAC (ACH).
!  19981111  Corrected Block B by adding final line, GO TO 200 (ACH).
!  20020430  Various upgrades (ACH): Use ODEPACK error handler package.
!            Replaced D1MACH by DUMACH.  Various changes to main
!            prologue and other routine prologues.
!-----------------------------------------------------------------------
! Other Routines in the DVODE Package.
!
! In addition to subroutine DVODE, the DVODE package includes the
! following subroutines and function routines:
!  DVHIN     computes an approximate step size for the initial step.
!  DVINDY    computes an interpolated value of the y vector at t = TOUT.
!  DVSTEP    is the core integrator, which does one step of the
!            integration and the associated error control.
!  DVSET     sets all method coefficients and test constants.
!  DVNLSD    solves the underlying nonlinear system -- the corrector.
!  DVJAC     computes and preprocesses the Jacobian matrix J = df/dy
!            and the Newton iteration matrix P = I - (h/l1)*J.
!  DVSOL     manages solution of linear system in chord iteration.
!  DVJUST    adjusts the history array on a change of order.
!  DEWSET    sets the error weight vector EWT before each step.
!  DVNORM    computes the weighted r.m.s. norm of a vector.
!  DVSRCO    is a user-callable routine to save and restore
!            the contents of the internal COMMON blocks.
!  DACOPY    is a routine to copy one two-dimensional array to another.
!  DGEFA and DGESL   are routines from LINPACK for solving full
!            systems of linear algebraic equations.
!  DGBFA and DGBSL   are routines from LINPACK for solving banded
!            linear systems.
!  DAXPY, DSCAL, and DCOPY are basic linear algebra modules (BLAS).
!  DUMACH    sets the unit roundoff of the machine.
!  XERRWD, XSETUN, XSETF, IXSAV, and IUMACH handle the printing of all
!            error messages and warnings.  XERRWD is machine-dependent.
! Note:  DVNORM, DUMACH, IXSAV, and IUMACH are function routines.
! All the others are subroutines.
!
!-----------------------------------------------------------------------
!
! Type declarations for labeled COMMON block DVOD01 --------------------
!
      DOUBLE PRECISION ACNrm , CCMxj , CONp , CRAte , DRC , EL , ETA ,  &
                       ETAmax , H , HMIn , HMXi , HNEw , HSCal , PRL1 , &
                       RC , RL1 , TAU , TQ , TN , UROund
      INTEGER ICF , INIt , IPUp , JCUr , JSTart , JSV , KFLag , KUTh ,  &
              L , LMAx , LYH , LEWt , LACor , LSAvf , LWM , LIWm ,      &
              LOCjs , MAXord , METh , MITer , MSBj , MXHnil , MXStep ,  &
              N , NEWh , NEWq , NHNil , NQ , NQNyh , NQWait , NSLj ,    &
              NSLp , NYH
!
! Type declarations for labeled COMMON block DVOD02 --------------------
!
      DOUBLE PRECISION HU
      INTEGER NCFn , NETf , NFE , NJE , NLU , NNI , NQU , NST
!
! Type declarations for local variables --------------------------------
!
!      EXTERNAL DVNLSD
      LOGICAL ihit
      DOUBLE PRECISION atoli , big , ewti , four , h0 , hmax , hmx ,    &
                       hun , one , pt2 , rh , rtoli , size , tcrit ,    &
                       tnext , tolsf , tp , two , zero
      INTEGER i , ier , iflag , imxer , jco , kgo , leniw , lenj ,      &
              lenp , lenrw , lenwm , lf0 , mband , mfa , ml , mord ,    &
              mu , mxhnl0 , mxstp0 , niter , nslast
      CHARACTER*80 msg
!
! Type declaration for function subroutines called ---------------------
!
!      DOUBLE PRECISION DUMACH , DVNORM
!
      DIMENSION mord(2)
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to DVODE.
!-----------------------------------------------------------------------
      SAVE mord , mxhnl0 , mxstp0
      SAVE zero , one , two , four , pt2 , hun
!-----------------------------------------------------------------------
! The following internal COMMON blocks contain variables which are
! communicated between subroutines in the DVODE package, or which are
! to be saved between calls to DVODE.
! In each block, real variables precede integers.
! The block /DVOD01/ appears in subroutines DVODE, DVINDY, DVSTEP,
! DVSET, DVNLSD, DVJAC, DVSOL, DVJUST and DVSRCO.
! The block /DVOD02/ appears in subroutines DVODE, DVINDY, DVSTEP,
! DVNLSD, DVJAC, and DVSRCO.
!
! The variables stored in the internal COMMON blocks are as follows:
!
! ACNRM  = Weighted r.m.s. norm of accumulated correction vectors.
! CCMXJ  = Threshhold on DRC for updating the Jacobian. (See DRC.)
! CONP   = The saved value of TQ(5).
! CRATE  = Estimated corrector convergence rate constant.
! DRC    = Relative change in H*RL1 since last DVJAC call.
! EL     = Real array of integration coefficients.  See DVSET.
! ETA    = Saved tentative ratio of new to old H.
! ETAMAX = Saved maximum value of ETA to be allowed.
! H      = The step size.
! HMIN   = The minimum absolute value of the step size H to be used.
! HMXI   = Inverse of the maximum absolute value of H to be used.
!          HMXI = 0.0 is allowed and corresponds to an infinite HMAX.
! HNEW   = The step size to be attempted on the next step.
! HSCAL  = Stepsize in scaling of YH array.
! PRL1   = The saved value of RL1.
! RC     = Ratio of current H*RL1 to value on last DVJAC call.
! RL1    = The reciprocal of the coefficient EL(1).
! TAU    = Real vector of past NQ step sizes, length 13.
! TQ     = A real vector of length 5 in which DVSET stores constants
!          used for the convergence test, the error test, and the
!          selection of H at a new order.
! TN     = The independent variable, updated on each step taken.
! UROUND = The machine unit roundoff.  The smallest positive real number
!          such that  1.0 + UROUND .ne. 1.0
! ICF    = Integer flag for convergence failure in DVNLSD:
!            0 means no failures.
!            1 means convergence failure with out of date Jacobian
!                   (recoverable error).
!            2 means convergence failure with current Jacobian or
!                   singular matrix (unrecoverable error).
! INIT   = Saved integer flag indicating whether initialization of the
!          problem has been done (INIT = 1) or not.
! IPUP   = Saved flag to signal updating of Newton matrix.
! JCUR   = Output flag from DVJAC showing Jacobian status:
!            JCUR = 0 means J is not current.
!            JCUR = 1 means J is current.
! JSTART = Integer flag used as input to DVSTEP:
!            0  means perform the first step.
!            1  means take a new step continuing from the last.
!            -1 means take the next step with a new value of MAXORD,
!                  HMIN, HMXI, N, METH, MITER, and/or matrix parameters.
!          On return, DVSTEP sets JSTART = 1.
! JSV    = Integer flag for Jacobian saving, = sign(MF).
! KFLAG  = A completion code from DVSTEP with the following meanings:
!               0      the step was succesful.
!              -1      the requested error could not be achieved.
!              -2      corrector convergence could not be achieved.
!              -3, -4  fatal error in VNLS (can not occur here).
! KUTH   = Input flag to DVSTEP showing whether H was reduced by the
!          driver.  KUTH = 1 if H was reduced, = 0 otherwise.
! L      = Integer variable, NQ + 1, current order plus one.
! LMAX   = MAXORD + 1 (used for dimensioning).
! LOCJS  = A pointer to the saved Jacobian, whose storage starts at
!          WM(LOCJS), if JSV = 1.
! LYH, LEWT, LACOR, LSAVF, LWM, LIWM = Saved integer pointers
!          to segments of RWORK and IWORK.
! MAXORD = The maximum order of integration method to be allowed.
! METH/MITER = The method flags.  See MF.
! MSBJ   = The maximum number of steps between J evaluations, = 50.
! MXHNIL = Saved value of optional input MXHNIL.
! MXSTEP = Saved value of optional input MXSTEP.
! N      = The number of first-order ODEs, = NEQ.
! NEWH   = Saved integer to flag change of H.
! NEWQ   = The method order to be used on the next step.
! NHNIL  = Saved counter for occurrences of T + H = T.
! NQ     = Integer variable, the current integration method order.
! NQNYH  = Saved value of NQ*NYH.
! NQWAIT = A counter controlling the frequency of order changes.
!          An order change is about to be considered if NQWAIT = 1.
! NSLJ   = The number of steps taken as of the last Jacobian update.
! NSLP   = Saved value of NST as of last Newton matrix update.
! NYH    = Saved value of the initial value of NEQ.
! HU     = The step size in t last used.
! NCFN   = Number of nonlinear convergence failures so far.
! NETF   = The number of error test failures of the integrator so far.
! NFE    = The number of f evaluations for the problem so far.
! NJE    = The number of Jacobian evaluations so far.
! NLU    = The number of matrix LU decompositions so far.
! NNI    = Number of nonlinear iterations so far.
! NQU    = The method order last used.
! NST    = The number of steps taken for the problem so far.
!-----------------------------------------------------------------------
      COMMON /DVOD01/ ACNrm , CCMxj , CONp , CRAte , DRC , EL(13) ,     &
                      ETA , ETAmax , H , HMIn , HMXi , HNEw , HSCal ,   &
                      PRL1 , RC , RL1 , TAU(13) , TQ(5) , TN , UROund , &
                      ICF , INIt , IPUp , JCUr , JSTart , JSV , KFLag , &
                      KUTh , L , LMAx , LYH , LEWt , LACor , LSAvf ,    &
                      LWM , LIWm , LOCjs , MAXord , METh , MITer ,      &
                      MSBj , MXHnil , MXStep , N , NEWh , NEWq , NHNil ,&
                      NQ , NQNyh , NQWait , NSLj , NSLp , NYH
      COMMON /DVOD02/ HU , NCFn , NETf , NFE , NJE , NLU , NNI , NQU ,  &
                      NST
!
      DATA mord(1)/12/ , mord(2)/5/ , mxstp0/500/ , mxhnl0/10/
      DATA zero/0.0D0/ , one/1.0D0/ , two/2.0D0/ , four/4.0D0/ ,        &
           pt2/0.2D0/ , hun/100.0D0/
!-----------------------------------------------------------------------
! Block A.
! This code block is executed on every call.
! It tests ISTATE and ITASK for legality and branches appropriately.
! If ISTATE .gt. 1 but the flag INIT shows that initialization has
! not yet been done, an error return occurs.
! If ISTATE = 1 and TOUT = T, return immediately.
!-----------------------------------------------------------------------
      IF ( Istate<1 .OR. Istate>3 ) THEN
!-----------------------------------------------------------------------
! Block I.
! The following block handles all error returns due to illegal input
! (ISTATE = -3), as detected before calling the core integrator.
! First the error message routine is called.   If the illegal input
! is a negative ISTATE, the run is aborted (apparent infinite loop).
!-----------------------------------------------------------------------
         msg = 'DVODE--  ISTATE (=I1) illegal '
         CALL XERRWD(msg,30,1,1,1,Istate,0,0,zero,zero)
         IF ( Istate>=0 ) GOTO 1500
!
         msg = 'DVODE--  Run aborted:  apparent infinite loop     '
         CALL XERRWD(msg,50,303,2,0,0,0,0,zero,zero)
         GOTO 99999
      ELSE
         IF ( Itask<1 .OR. Itask>5 ) THEN
            msg = 'DVODE--  ITASK (=I1) illegal  '
            CALL XERRWD(msg,30,2,1,1,Itask,0,0,zero,zero)
            GOTO 1500
         ELSE
            IF ( Istate==1 ) THEN
               INIt = 0
               IF ( Tout==T ) RETURN
            ELSEIF ( INIt/=1 ) THEN
               msg =                                                    &
          'DVODE--  ISTATE (=I1) .gt. 1 but DVODE not initialized      '
               CALL XERRWD(msg,60,3,1,1,Istate,0,0,zero,zero)
               GOTO 1500
            ELSEIF ( Istate==2 ) THEN
               GOTO 50
            ENDIF
!-----------------------------------------------------------------------
! Block B.
! The next code block is executed for the initial call (ISTATE = 1),
! or for a continuation call with parameter changes (ISTATE = 3).
! It contains checking of all input and various initializations.
!
! First check legality of the non-optional input NEQ, ITOL, IOPT,
! MF, ML, and MU.
!-----------------------------------------------------------------------
            IF ( Neq<=0 ) THEN
               msg = 'DVODE--  NEQ (=I1) .lt. 1     '
               CALL XERRWD(msg,30,4,1,1,Neq,0,0,zero,zero)
               GOTO 1500
            ELSE
               IF ( Istate/=1 ) THEN
                  IF ( Neq>N ) THEN
                     msg =                                              &
                    'DVODE--  ISTATE = 3 and NEQ increased (I1 to I2)  '
                     CALL XERRWD(msg,50,5,1,2,N,Neq,0,zero,zero)
                     GOTO 1500
                  ENDIF
               ENDIF
               N = Neq
               IF ( Itol<1 .OR. Itol>4 ) THEN
                  msg = 'DVODE--  ITOL (=I1) illegal   '
                  CALL XERRWD(msg,30,6,1,1,Itol,0,0,zero,zero)
                  GOTO 1500
               ELSEIF ( Iopt<0 .OR. Iopt>1 ) THEN
                  msg = 'DVODE--  IOPT (=I1) illegal   '
                  CALL XERRWD(msg,30,7,1,1,Iopt,0,0,zero,zero)
                  GOTO 1500
               ELSE
                  JSV = SIGN(1,Mf)
                  mfa = ABS(Mf)
                  METh = mfa/10
                  MITer = mfa - 10*METh
                  IF ( METh<1 .OR. METh>2 ) GOTO 800
                  IF ( MITer<0 .OR. MITer>5 ) GOTO 800
                  IF ( MITer>3 ) THEN
                     ml = Iwork(1)
                     mu = Iwork(2)
                     IF ( ml<0 .OR. ml>=N ) THEN
                        msg =                                           &
                    'DVODE--  ML (=I1) illegal:  .lt.0 or .ge.NEQ (=I2)'
                        CALL XERRWD(msg,50,9,1,2,ml,Neq,0,zero,zero)
                        GOTO 1500
                     ELSEIF ( mu<0 .OR. mu>=N ) THEN
                        msg =                                           &
                    'DVODE--  MU (=I1) illegal:  .lt.0 or .ge.NEQ (=I2)'
                        CALL XERRWD(msg,50,10,1,2,mu,Neq,0,zero,zero)
                        GOTO 1500
                     ENDIF
                  ENDIF
! Next process and check the optional input. ---------------------------
                  IF ( Iopt==1 ) THEN
                     MAXord = Iwork(5)
                     IF ( MAXord<0 ) THEN
                        msg = 'DVODE--  MAXORD (=I1) .lt. 0  '
                        CALL XERRWD(msg,30,11,1,1,MAXord,0,0,zero,zero)
                        GOTO 1500
                     ELSE
                        IF ( MAXord==0 ) MAXord = 100
                        MAXord = MIN(MAXord,mord(METh))
                        MXStep = Iwork(6)
                        IF ( MXStep<0 ) THEN
                           msg = 'DVODE--  MXSTEP (=I1) .lt. 0  '
                           CALL XERRWD(msg,30,12,1,1,MXStep,0,0,zero,   &
                                       zero)
                           GOTO 1500
                        ELSE
                           IF ( MXStep==0 ) MXStep = mxstp0
                           MXHnil = Iwork(7)
                           IF ( MXHnil<0 ) THEN
                              msg = 'DVODE--  MXHNIL (=I1) .lt. 0  '
                              CALL XERRWD(msg,30,13,1,1,MXHnil,0,0,zero,&
                                 zero)
                              GOTO 1500
                           ELSE
                              IF ( MXHnil==0 ) MXHnil = mxhnl0
                              IF ( Istate==1 ) THEN
                                 h0 = Rwork(5)
                                 IF ( (Tout-T)*h0<zero ) THEN
                                    msg =                               &
                              'DVODE--  TOUT (=R1) behind T (=R2)      '
                                    CALL XERRWD(msg,40,14,1,0,0,0,2,    &
                                       Tout,T)
                                    msg =                               &
                    '      integration direction is given by H0 (=R1)  '
                                    CALL XERRWD(msg,50,14,1,0,0,0,1,h0, &
                                       zero)
                                    GOTO 1500
                                 ENDIF
                              ENDIF
                              hmax = Rwork(6)
                              IF ( hmax<zero ) THEN
                                 msg = 'DVODE--  HMAX (=R1) .lt. 0.0  '
                                 CALL XERRWD(msg,30,15,1,0,0,0,1,hmax,  &
                                    zero)
                                 GOTO 1500
                              ELSE
                                 HMXi = zero
                                 IF ( hmax>zero ) HMXi = one/hmax
                                 HMIn = Rwork(7)
                                 IF ( HMIn<zero ) THEN
                                    msg =                               &
                                       'DVODE--  HMIN (=R1) .lt. 0.0  '
                                    CALL XERRWD(msg,30,16,1,0,0,0,1,    &
                                       HMIn,zero)
                                    GOTO 1500
                                 ENDIF
                              ENDIF
                           ENDIF
                        ENDIF
                     ENDIF
                  ELSE
                     MAXord = mord(METh)
                     MXStep = mxstp0
                     MXHnil = mxhnl0
                     IF ( Istate==1 ) h0 = zero
                     HMXi = zero
                     HMIn = zero
                  ENDIF
!-----------------------------------------------------------------------
! Set work array pointers and check lengths LRW and LIW.
! Pointers to segments of RWORK and IWORK are named by prefixing L to
! the name of the segment.  E.g., the segment YH starts at RWORK(LYH).
! Segments of RWORK (in order) are denoted  YH, WM, EWT, SAVF, ACOR.
! Within WM, LOCJS is the location of the saved Jacobian (JSV .gt. 0).
!-----------------------------------------------------------------------
                  LYH = 21
                  IF ( Istate==1 ) NYH = N
                  LWM = LYH + (MAXord+1)*NYH
                  jco = MAX(0,JSV)
                  IF ( MITer==0 ) lenwm = 0
                  IF ( MITer==1 .OR. MITer==2 ) THEN
                     lenwm = 2 + (1+jco)*N*N
                     LOCjs = N*N + 3
                  ENDIF
                  IF ( MITer==3 ) lenwm = 2 + N
                  IF ( MITer==4 .OR. MITer==5 ) THEN
                     mband = ml + mu + 1
                     lenp = (mband+ml)*N
                     lenj = mband*N
                     lenwm = 2 + lenp + jco*lenj
                     LOCjs = lenp + 3
                  ENDIF
                  LEWt = LWM + lenwm
                  LSAvf = LEWt + N
                  LACor = LSAvf + N
                  lenrw = LACor + N - 1
                  Iwork(17) = lenrw
                  LIWm = 1
                  leniw = 30 + N
                  IF ( MITer==0 .OR. MITer==3 ) leniw = 30
                  Iwork(18) = leniw
                  IF ( lenrw>Lrw ) THEN
                     msg =                                              &
          'DVODE--  RWORK length needed, LENRW (=I1), exceeds LRW (=I2)'
                     CALL XERRWD(msg,60,17,1,2,lenrw,Lrw,0,zero,zero)
                     GOTO 1500
                  ELSEIF ( leniw>Liw ) THEN
                     msg =                                              &
          'DVODE--  IWORK length needed, LENIW (=I1), exceeds LIW (=I2)'
                     CALL XERRWD(msg,60,18,1,2,leniw,Liw,0,zero,zero)
                     GOTO 1500
                  ELSE
! Check RTOL and ATOL for legality. ------------------------------------
                     rtoli = Rtol(1)
                     atoli = Atol(1)
                     DO i = 1 , N
                        IF ( Itol>=3 ) rtoli = Rtol(i)
                        IF ( Itol==2 .OR. Itol==4 ) atoli = Atol(i)
                        IF ( rtoli<zero ) GOTO 900
                        IF ( atoli<zero ) GOTO 1000
                     ENDDO
                     IF ( Istate==1 ) THEN
!-----------------------------------------------------------------------
! Block C.
! The next block is for the initial call only (ISTATE = 1).
! It contains all remaining initializations, the initial call to F,
! and the calculation of the initial step size.
! The error weights in EWT are inverted after being loaded.
!-----------------------------------------------------------------------
                        UROund = DUMACH()
                        TN = T
                        IF ( Itask==4 .OR. Itask==5 ) THEN
                           tcrit = Rwork(1)
                           IF ( (tcrit-Tout)*(Tout-T)<zero ) GOTO 1300
                           IF ( h0/=zero .AND. (T+h0-tcrit)*h0>zero )   &
                                h0 = tcrit - T
                        ENDIF
                        JSTart = 0
                        IF ( MITer>0 ) Rwork(LWM) = SQRT(UROund)
                        CCMxj = pt2
                        MSBj = 50
                        NHNil = 0
                        NST = 0
                        NJE = 0
                        NNI = 0
                        NCFn = 0
                        NETf = 0
                        NLU = 0
                        NSLj = 0
                        nslast = 0
                        HU = zero
                        NQU = 0
! Initial call to F.  (LF0 points to YH(*,2).) -------------------------
                        lf0 = LYH + NYH
                        CALL F(N,T,Y,Rwork(lf0),Rpar,Ipar)
                        NFE = 1
! Load the initial value vector in YH. ---------------------------------
                        CALL DCOPY(N,Y,1,Rwork(LYH),1)
! Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
                        NQ = 1
                        H = one
                        CALL DEWSET(N,Itol,Rtol,Atol,Rwork(LYH),        &
                                    Rwork(LEWt))
                        DO i = 1 , N
                           IF ( Rwork(i+LEWt-1)<=zero ) GOTO 1100
                           Rwork(i+LEWt-1) = one/Rwork(i+LEWt-1)
                        ENDDO
                        IF ( h0==zero ) THEN
! Call DVHIN to set initial step size H0 to be attempted. --------------
                           CALL DVHIN(N,T,Rwork(LYH),Rwork(lf0),F,Rpar, &
                                      Ipar,Tout,UROund,Rwork(LEWt),Itol,&
                                      Atol,Y,Rwork(LACor),h0,niter,ier)
                           NFE = NFE + niter
                           IF ( ier/=0 ) THEN
                              msg =                                     &
          'DVODE--  TOUT (=R1) too close to T(=R2) to start integration'
                              CALL XERRWD(msg,60,22,1,0,0,0,2,Tout,T)
                              GOTO 1500
                           ENDIF
                        ENDIF
! Adjust H0 if necessary to meet HMAX bound. ---------------------------
                        rh = ABS(h0)*HMXi
                        IF ( rh>one ) h0 = h0/rh
! Load H with H0 and scale YH(*,2) by H0. ------------------------------
                        H = h0
                        CALL DSCAL(N,h0,Rwork(lf0),1)
                        GOTO 200
                     ELSE
! If ISTATE = 3, set flag to signal parameter changes to DVSTEP. -------
                        JSTart = -1
! MAXORD was reduced below NQ.  Copy YH(*,MAXORD+2) into SAVF. ---------
                        IF ( NQ>MAXord )                                &
                             CALL DCOPY(N,Rwork(LWM),1,Rwork(LSAvf),1)
! Reload WM(1) = RWORK(LWM), since LWM may have changed. ---------------
                        IF ( MITer>0 ) Rwork(LWM) = SQRT(UROund)
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
!-----------------------------------------------------------------------
! Block D.
! The next code block is for continuation calls only (ISTATE = 2 or 3)
! and is to check stop conditions before taking a step.
!-----------------------------------------------------------------------
 50      nslast = NST
         KUTh = 0
         SELECT CASE (Itask)
         CASE (2)
            GOTO 100
         CASE (3)
            tp = TN - HU*(one+hun*UROund)
            IF ( (tp-Tout)*H>zero ) THEN
               msg =                                                    &
          'DVODE--  ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)  '
               CALL XERRWD(msg,60,23,1,1,Itask,0,2,Tout,tp)
               GOTO 1500
            ELSE
               IF ( (TN-Tout)*H>=zero ) GOTO 300
               GOTO 100
            ENDIF
         CASE (4)
            tcrit = Rwork(1)
            IF ( (TN-tcrit)*H>zero ) GOTO 1200
            IF ( (tcrit-Tout)*H<zero ) GOTO 1300
            IF ( (TN-Tout)*H>=zero ) THEN
               CALL DVINDY(Tout,0,Rwork(LYH),NYH,Y,iflag)
               IF ( iflag/=0 ) GOTO 1400
               T = Tout
               GOTO 400
            ENDIF
         CASE (5)
            tcrit = Rwork(1)
            IF ( (TN-tcrit)*H>zero ) GOTO 1200
         CASE DEFAULT
            IF ( (TN-Tout)*H<zero ) GOTO 100
            CALL DVINDY(Tout,0,Rwork(LYH),NYH,Y,iflag)
            IF ( iflag/=0 ) GOTO 1400
            T = Tout
            GOTO 400
         END SELECT
         hmx = ABS(TN) + ABS(H)
         ihit = ABS(TN-tcrit)<=hun*UROund*hmx
         IF ( ihit ) GOTO 300
         tnext = TN + HNEw*(one+four*UROund)
         IF ( (tnext-tcrit)*H>zero ) THEN
            H = (tcrit-TN)*(one-four*UROund)
            KUTh = 1
         ENDIF
      ENDIF
!-----------------------------------------------------------------------
! Block E.
! The next block is normally executed for all calls and contains
! the call to the one-step core integrator DVSTEP.
!
! This is a looping point for the integration steps.
!
! First check for too many steps being taken, update EWT (if not at
! start of problem), check for too much accuracy being requested, and
! check for H below the roundoff level in T.
!-----------------------------------------------------------------------
 100  IF ( (NST-nslast)>=MXStep ) THEN
!-----------------------------------------------------------------------
! Block H.
! The following block handles all unsuccessful returns other than
! those for illegal input.  First the error message routine is called.
! if there was an error test or convergence test failure, IMXER is set.
! Then Y is loaded from YH, and T is set to TN.
! The optional output is loaded into the work arrays before returning.
!-----------------------------------------------------------------------
! The maximum number of steps was taken before reaching TOUT. ----------
         msg = 'DVODE--  At current T (=R1), MXSTEP (=I1) steps   '
         CALL XERRWD(msg,50,201,1,0,0,0,0,zero,zero)
         msg = '      taken on this call before reaching TOUT     '
         CALL XERRWD(msg,50,201,1,1,MXStep,0,1,TN,zero)
         Istate = -1
         GOTO 700
      ELSE
         CALL DEWSET(N,Itol,Rtol,Atol,Rwork(LYH),Rwork(LEWt))
         DO i = 1 , N
            IF ( Rwork(i+LEWt-1)<=zero ) GOTO 500
            Rwork(i+LEWt-1) = one/Rwork(i+LEWt-1)
         ENDDO
      ENDIF
 200  tolsf = UROund*DVNORM(N,Rwork(LYH),Rwork(LEWt))
      IF ( tolsf<=one ) THEN
         IF ( (TN+H)==TN ) THEN
            NHNil = NHNil + 1
            IF ( NHNil<=MXHnil ) THEN
               msg =                                                    &
                    'DVODE--  Warning: internal T (=R1) and H (=R2) are'
               CALL XERRWD(msg,50,101,1,0,0,0,0,zero,zero)
               msg =                                                    &
          '      such that in the machine, T + H = T on the next step  '
               CALL XERRWD(msg,60,101,1,0,0,0,0,zero,zero)
               msg =                                                    &
                    '      (H = step size). solver will continue anyway'
               CALL XERRWD(msg,50,101,1,0,0,0,2,TN,H)
               IF ( NHNil>=MXHnil ) THEN
                  msg =                                                 &
                    'DVODE--  Above warning has been issued I1 times.  '
                  CALL XERRWD(msg,50,102,1,0,0,0,0,zero,zero)
                  msg =                                                 &
                    '      it will not be issued again for this problem'
                  CALL XERRWD(msg,50,102,1,1,MXHnil,0,0,zero,zero)
               ENDIF
            ENDIF
         ENDIF
!-----------------------------------------------------------------------
! CALL DVSTEP (Y, YH, NYH, YH, EWT, SAVF, VSAV, ACOR,
!              WM, IWM, F, JAC, F, DVNLSD, RPAR, IPAR)
!-----------------------------------------------------------------------
         CALL DVSTEP(Y,Rwork(LYH),NYH,Rwork(LYH),Rwork(LEWt),           &
                     Rwork(LSAvf),Y,Rwork(LACor),Rwork(LWM),Iwork(LIWm),&
                     F,JAC,F,DVNLSD,Rpar,Ipar)
         kgo = 1 - KFLag
! Branch on KFLAG.  Note: In this version, KFLAG can not be set to -3.
!  KFLAG .eq. 0,   -1,  -2
         SELECT CASE (kgo)
         CASE (2)
! KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
            msg = 'DVODE--  At T(=R1) and step size H(=R2), the error'
            CALL XERRWD(msg,50,204,1,0,0,0,0,zero,zero)
            msg = '      test failed repeatedly or with abs(H) = HMIN'
            CALL XERRWD(msg,50,204,1,0,0,0,2,TN,H)
            Istate = -4
            GOTO 600
         CASE (3)
! KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ----
            msg = 'DVODE--  At T (=R1) and step size H (=R2), the    '
            CALL XERRWD(msg,50,205,1,0,0,0,0,zero,zero)
            msg = '      corrector convergence failed repeatedly     '
            CALL XERRWD(msg,50,205,1,0,0,0,0,zero,zero)
            msg = '      or with abs(H) = HMIN   '
            CALL XERRWD(msg,30,205,1,0,0,0,2,TN,H)
            Istate = -5
            GOTO 600
         CASE DEFAULT
!-----------------------------------------------------------------------
! Block F.
! The following block handles the case of a successful return from the
! core integrator (KFLAG = 0).  Test for stop conditions.
!-----------------------------------------------------------------------
            INIt = 1
            KUTh = 0
            SELECT CASE (Itask)
            CASE (2)
            CASE (3)
! ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
               IF ( (TN-Tout)*H<zero ) GOTO 100
            CASE (4)
! ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary.
               IF ( (TN-Tout)*H<zero ) THEN
                  hmx = ABS(TN) + ABS(H)
                  ihit = ABS(TN-tcrit)<=hun*UROund*hmx
                  IF ( .NOT.(ihit) ) THEN
                     tnext = TN + HNEw*(one+four*UROund)
                     IF ( (tnext-tcrit)*H>zero ) THEN
                        H = (tcrit-TN)*(one-four*UROund)
                        KUTh = 1
                     ENDIF
                     GOTO 100
                  ENDIF
               ELSE
                  CALL DVINDY(Tout,0,Rwork(LYH),NYH,Y,iflag)
                  T = Tout
                  GOTO 400
               ENDIF
            CASE (5)
! ITASK = 5.  See if TCRIT was reached and jump to exit. ---------------
               hmx = ABS(TN) + ABS(H)
               ihit = ABS(TN-tcrit)<=hun*UROund*hmx
            CASE DEFAULT
! ITASK = 1.  If TOUT has been reached, interpolate. -------------------
               IF ( (TN-Tout)*H<zero ) GOTO 100
               CALL DVINDY(Tout,0,Rwork(LYH),NYH,Y,iflag)
               T = Tout
               GOTO 400
            END SELECT
         END SELECT
      ELSE
         tolsf = tolsf*two
         IF ( NST==0 ) THEN
            msg = 'DVODE--  At start of problem, too much accuracy   '
            CALL XERRWD(msg,50,26,1,0,0,0,0,zero,zero)
            msg =                                                       &
          '      requested for precision of machine:   see TOLSF (=R1) '
            CALL XERRWD(msg,60,26,1,0,0,0,1,tolsf,zero)
            Rwork(14) = tolsf
            GOTO 1500
         ELSE
! Too much accuracy requested for machine precision. -------------------
            msg = 'DVODE--  At T (=R1), too much accuracy requested  '
            CALL XERRWD(msg,50,203,1,0,0,0,0,zero,zero)
            msg = '      for precision of machine:   see TOLSF (=R2) '
            CALL XERRWD(msg,50,203,1,0,0,0,2,TN,tolsf)
            Rwork(14) = tolsf
            Istate = -2
            GOTO 700
         ENDIF
      ENDIF
!-----------------------------------------------------------------------
! Block G.
! The following block handles all successful returns from DVODE.
! If ITASK .ne. 1, Y is loaded from YH and T is set accordingly.
! ISTATE is set to 2, and the optional output is loaded into the work
! arrays before returning.
!-----------------------------------------------------------------------
 300  CALL DCOPY(N,Rwork(LYH),1,Y,1)
      T = TN
      IF ( Itask==4 .OR. Itask==5 ) THEN
         IF ( ihit ) T = tcrit
      ENDIF
 400  Istate = 2
      Rwork(11) = HU
      Rwork(12) = HNEw
      Rwork(13) = TN
      Iwork(11) = NST
      Iwork(12) = NFE
      Iwork(13) = NJE
      Iwork(14) = NQU
      Iwork(15) = NEWq
      Iwork(19) = NLU
      Iwork(20) = NNI
      Iwork(21) = NCFn
      Iwork(22) = NETf
      RETURN
! EWT(i) .le. 0.0 for some i (not at start of problem). ----------------
 500  ewti = Rwork(LEWt+i-1)
      msg = 'DVODE--  At T (=R1), EWT(I1) has become R2 .le. 0.'
      CALL XERRWD(msg,50,202,1,1,i,0,2,TN,ewti)
      Istate = -6
      GOTO 700
! Compute IMXER if relevant. -------------------------------------------
 600  big = zero
      imxer = 1
      DO i = 1 , N
         size = ABS(Rwork(i+LACor-1)*Rwork(i+LEWt-1))
         IF ( big<size ) THEN
            big = size
            imxer = i
         ENDIF
      ENDDO
      Iwork(16) = imxer
! Set Y vector, T, and optional output. --------------------------------
 700  CALL DCOPY(N,Rwork(LYH),1,Y,1)
      T = TN
      Rwork(11) = HU
      Rwork(12) = H
      Rwork(13) = TN
      Iwork(11) = NST
      Iwork(12) = NFE
      Iwork(13) = NJE
      Iwork(14) = NQU
      Iwork(15) = NQ
      Iwork(19) = NLU
      Iwork(20) = NNI
      Iwork(21) = NCFn
      Iwork(22) = NETf
      RETURN
 800  msg = 'DVODE--  MF (=I1) illegal     '
      CALL XERRWD(msg,30,8,1,1,Mf,0,0,zero,zero)
      GOTO 1500
 900  msg = 'DVODE--  RTOL(I1) is R1 .lt. 0.0        '
      CALL XERRWD(msg,40,19,1,1,i,0,1,rtoli,zero)
      GOTO 1500
 1000 msg = 'DVODE--  ATOL(I1) is R1 .lt. 0.0        '
      CALL XERRWD(msg,40,20,1,1,i,0,1,atoli,zero)
      GOTO 1500
 1100 ewti = Rwork(LEWt+i-1)
      msg = 'DVODE--  EWT(I1) is R1 .le. 0.0         '
      CALL XERRWD(msg,40,21,1,1,i,0,1,ewti,zero)
      GOTO 1500
 1200 msg =                                                             &
          'DVODE--  ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)   '
      CALL XERRWD(msg,60,24,1,0,0,0,2,tcrit,TN)
      GOTO 1500
 1300 msg =                                                             &
          'DVODE--  ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)   '
      CALL XERRWD(msg,60,25,1,0,0,0,2,tcrit,Tout)
      GOTO 1500
 1400 msg =                                                             &
          'DVODE--  Trouble from DVINDY.  ITASK = I1, TOUT = R1.       '
      CALL XERRWD(msg,60,27,1,1,Itask,0,1,Tout,zero)
!
 1500 Istate = -3
      RETURN
!----------------------- End of Subroutine DVODE -----------------------
99999 END Subroutine DVODE
!DECK DVHIN
      SUBROUTINE DVHIN(N,T0,Y0,Ydot,F,Rpar,Ipar,Tout,Uround,Ewt,Itol,   &
                       Atol,Y,Temp,H0,Niter,Ier)
      IMPLICIT NONE
      EXTERNAL F
      DOUBLE PRECISION T0 , Y0 , Ydot , Rpar , Tout , Uround , Ewt ,    &
                       Atol , Y , Temp , H0
      INTEGER N , Ipar , Itol , Niter , Ier
      DIMENSION Y0(*) , Ydot(*) , Ewt(*) , Atol(*) , Y(*) , Temp(*) ,   &
                Rpar(*) , Ipar(*)
!-----------------------------------------------------------------------
! Call sequence input -- N, T0, Y0, YDOT, F, RPAR, IPAR, TOUT, UROUND,
!                        EWT, ITOL, ATOL, Y, TEMP
! Call sequence output -- H0, NITER, IER
! COMMON block variables accessed -- None
!
! Subroutines called by DVHIN:  F
! Function routines called by DVHI: DVNORM
!-----------------------------------------------------------------------
! This routine computes the step size, H0, to be attempted on the
! first step, when the user has not supplied a value for this.
!
! First we check that TOUT - T0 differs significantly from zero.  Then
! an iteration is done to approximate the initial second derivative
! and this is used to define h from w.r.m.s.norm(h**2 * yddot / 2) = 1.
! A bias factor of 1/2 is applied to the resulting h.
! The sign of H0 is inferred from the initial values of TOUT and T0.
!
! Communication with DVHIN is done with the following variables:
!
! N      = Size of ODE system, input.
! T0     = Initial value of independent variable, input.
! Y0     = Vector of initial conditions, input.
! YDOT   = Vector of initial first derivatives, input.
! F      = Name of subroutine for right-hand side f(t,y), input.
! RPAR, IPAR = Dummy names for user's real and integer work arrays.
! TOUT   = First output value of independent variable
! UROUND = Machine unit roundoff
! EWT, ITOL, ATOL = Error weights and tolerance parameters
!                   as described in the driver routine, input.
! Y, TEMP = Work arrays of length N.
! H0     = Step size to be attempted, output.
! NITER  = Number of iterations (and of f evaluations) to compute H0,
!          output.
! IER    = The error flag, returned with the value
!          IER = 0  if no trouble occurred, or
!          IER = -1 if TOUT and T0 are considered too close to proceed.
!-----------------------------------------------------------------------
!
! Type declarations for local variables --------------------------------
!
      DOUBLE PRECISION afi , atoli , delyi , h , half , hg , hlb ,      &
                       hnew , hrat , hub , hun , pt1 , t1 , tdist ,     &
                       tround , two , yddnrm
      INTEGER i , iter
!
! Type declaration for function subroutines called ---------------------
!
!      DOUBLE PRECISION DVNORM
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
      SAVE half , hun , pt1 , two
      DATA half/0.5D0/ , hun/100.0D0/ , pt1/0.1D0/ , two/2.0D0/
!
      Niter = 0
      tdist = ABS(Tout-T0)
      tround = Uround*MAX(ABS(T0),ABS(Tout))
      IF ( tdist<two*tround ) THEN
! Error return for TOUT - T0 too small. --------------------------------
         Ier = -1
         GOTO 99999
      ELSE
!
! Set a lower bound on h based on the roundoff level in T0 and TOUT. ---
         hlb = hun*tround
! Set an upper bound on h based on TOUT-T0 and the initial Y and YDOT. -
         hub = pt1*tdist
         atoli = Atol(1)
         DO i = 1 , N
            IF ( Itol==2 .OR. Itol==4 ) atoli = Atol(i)
            delyi = pt1*ABS(Y0(i)) + atoli
            afi = ABS(Ydot(i))
            IF ( afi*hub>delyi ) hub = delyi/afi
         ENDDO
!
! Set initial guess for h as geometric mean of upper and lower bounds. -
         iter = 0
         hg = SQRT(hlb*hub)
! If the bounds have crossed, exit with the mean value. ----------------
         IF ( hub<hlb ) THEN
            H0 = hg
            GOTO 200
         ENDIF
!
! Looping point for iteration. -----------------------------------------
! Estimate the second derivative as a difference quotient in f. --------
 50      h = SIGN(hg,Tout-T0)
         t1 = T0 + h
         DO i = 1 , N
            Y(i) = Y0(i) + h*Ydot(i)
         ENDDO
         CALL F(N,t1,Y,Temp,Rpar,Ipar)
         DO i = 1 , N
            Temp(i) = (Temp(i)-Ydot(i))/h
         ENDDO
         yddnrm = DVNORM(N,Temp,Ewt)
! Get the corresponding new value of h. --------------------------------
         IF ( yddnrm*hub*hub>two ) THEN
            hnew = SQRT(two/yddnrm)
         ELSE
            hnew = SQRT(hg*hub)
         ENDIF
         iter = iter + 1
!-----------------------------------------------------------------------
! Test the stopping conditions.
! Stop if the new and previous h values differ by a factor of .lt. 2.
! Stop if four iterations have been done.  Also, stop with previous h
! if HNEW/HG .gt. 2 after first iteration, as this probably means that
! the second derivative value is bad because of cancellation error.
!-----------------------------------------------------------------------
         IF ( iter<4 ) THEN
            hrat = hnew/hg
            IF ( (hrat<=half) .OR. (hrat>=two) ) THEN
               IF ( (iter>=2) .AND. (hnew>two*hg) ) THEN
                  hnew = hg
                  GOTO 100
               ENDIF
               hg = hnew
               GOTO 50
            ENDIF
         ENDIF
!
! Iteration done.  Apply bounds, bias factor, and sign.  Then exit. ----
 100     H0 = hnew*half
         IF ( H0<hlb ) H0 = hlb
         IF ( H0>hub ) H0 = hub
      ENDIF
 200  H0 = SIGN(H0,Tout-T0)
      Niter = iter
      Ier = 0
      RETURN
!----------------------- End of Subroutine DVHIN -----------------------
99999 END Subroutine DVHIN
!DECK DVINDY
      SUBROUTINE DVINDY(T,K,Yh,Ldyh,Dky,Iflag)
      IMPLICIT NONE
      DOUBLE PRECISION T , Yh , Dky
      INTEGER K , Ldyh , Iflag
      DIMENSION Yh(Ldyh,*) , Dky(*)
!-----------------------------------------------------------------------
! Call sequence input -- T, K, YH, LDYH
! Call sequence output -- DKY, IFLAG
! COMMON block variables accessed:
!     /DVOD01/ --  H, TN, UROUND, L, N, NQ
!     /DVOD02/ --  HU
!
! Subroutines called by DVINDY: DSCAL, XERRWD
! Function routines called by DVINDY: None
!-----------------------------------------------------------------------
! DVINDY computes interpolated values of the K-th derivative of the
! dependent variable vector y, and stores it in DKY.  This routine
! is called within the package with K = 0 and T = TOUT, but may
! also be called by the user for any K up to the current order.
! (See detailed instructions in the usage documentation.)
!-----------------------------------------------------------------------
! The computed values in DKY are gotten by interpolation using the
! Nordsieck history array YH.  This array corresponds uniquely to a
! vector-valued polynomial of degree NQCUR or less, and DKY is set
! to the K-th derivative of this polynomial at T.
! The formula for DKY is:
!              q
!  DKY(i)  =  sum  c(j,K) * (T - TN)**(j-K) * H**(-j) * YH(i,j+1)
!             j=K
! where  c(j,K) = j*(j-1)*...*(j-K+1), q = NQCUR, TN = TCUR, H = HCUR.
! The quantities  NQ = NQCUR, L = NQ+1, N, TN, and H are
! communicated by COMMON.  The above sum is done in reverse order.
! IFLAG is returned negative if either K or T is out of bounds.
!
! Discussion above and comments in driver explain all variables.
!-----------------------------------------------------------------------
!
! Type declarations for labeled COMMON block DVOD01 --------------------
!
      DOUBLE PRECISION ACNrm , CCMxj , CONp , CRAte , DRC , EL , ETA ,  &
                       ETAmax , H , HMIn , HMXi , HNEw , HSCal , PRL1 , &
                       RC , RL1 , TAU , TQ , TN , UROund
      INTEGER ICF , INIt , IPUp , JCUr , JSTart , JSV , KFLag , KUTh ,  &
              L , LMAx , LYH , LEWt , LACor , LSAvf , LWM , LIWm ,      &
              LOCjs , MAXord , METh , MITer , MSBj , MXHnil , MXStep ,  &
              N , NEWh , NEWq , NHNil , NQ , NQNyh , NQWait , NSLj ,    &
              NSLp , NYH
!
! Type declarations for labeled COMMON block DVOD02 --------------------
!
      DOUBLE PRECISION HU
      INTEGER NCFn , NETf , NFE , NJE , NLU , NNI , NQU , NST
!
! Type declarations for local variables --------------------------------
!
      DOUBLE PRECISION c , hun , r , s , tfuzz , tn1 , tp , zero
      INTEGER i , ic , j , jb , jb2 , jj , jj1 , jp1
      CHARACTER*80 msg
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
      SAVE hun , zero
!
      COMMON /DVOD01/ ACNrm , CCMxj , CONp , CRAte , DRC , EL(13) ,     &
                      ETA , ETAmax , H , HMIn , HMXi , HNEw , HSCal ,   &
                      PRL1 , RC , RL1 , TAU(13) , TQ(5) , TN , UROund , &
                      ICF , INIt , IPUp , JCUr , JSTart , JSV , KFLag , &
                      KUTh , L , LMAx , LYH , LEWt , LACor , LSAvf ,    &
                      LWM , LIWm , LOCjs , MAXord , METh , MITer ,      &
                      MSBj , MXHnil , MXStep , N , NEWh , NEWq , NHNil ,&
                      NQ , NQNyh , NQWait , NSLj , NSLp , NYH
      COMMON /DVOD02/ HU , NCFn , NETf , NFE , NJE , NLU , NNI , NQU ,  &
                      NST
!
      DATA hun/100.0D0/ , zero/0.0D0/
!
      Iflag = 0
      IF ( K<0 .OR. K>NQ ) THEN
!
         msg = 'DVINDY-- K (=I1) illegal      '
         CALL XERRWD(msg,30,51,1,1,K,0,0,zero,zero)
         Iflag = -1
         RETURN
      ELSE
         tfuzz = hun*UROund*SIGN(ABS(TN)+ABS(HU),HU)
         tp = TN - HU - tfuzz
         tn1 = TN + tfuzz
         IF ( (T-tp)*(T-tn1)>zero ) THEN
            msg = 'DVINDY-- T (=R1) illegal      '
            CALL XERRWD(msg,30,52,1,0,0,0,1,T,zero)
            msg =                                                       &
          '      T not in interval TCUR - HU (= R1) to TCUR (=R2)      '
            CALL XERRWD(msg,60,52,1,0,0,0,2,tp,TN)
            Iflag = -2
            GOTO 99999
         ELSE
!
            s = (T-TN)/H
            ic = 1
            IF ( K/=0 ) THEN
               jj1 = L - K
               DO jj = jj1 , NQ
                  ic = ic*jj
               ENDDO
            ENDIF
            c = REAL(ic)
            DO i = 1 , N
               Dky(i) = c*Yh(i,L)
            ENDDO
            IF ( K/=NQ ) THEN
               jb2 = NQ - K
               DO jb = 1 , jb2
                  j = NQ - jb
                  jp1 = j + 1
                  ic = 1
                  IF ( K/=0 ) THEN
                     jj1 = jp1 - K
                     DO jj = jj1 , j
                        ic = ic*jj
                     ENDDO
                  ENDIF
                  c = REAL(ic)
                  DO i = 1 , N
                     Dky(i) = c*Yh(i,jp1) + s*Dky(i)
                  ENDDO
               ENDDO
               IF ( K==0 ) RETURN
            ENDIF
         ENDIF
      ENDIF
      r = H**(-K)
      CALL DSCAL(N,r,Dky,1)
      RETURN
!----------------------- End of Subroutine DVINDY ----------------------
99999 END Subroutine DVINDY
!DECK DVSTEP
      SUBROUTINE DVSTEP(Y,Yh,Ldyh,Yh1,Ewt,Savf,Vsav,Acor,Wm,Iwm,F,JAC,  &
                        PSOL,VNLS,Rpar,Ipar)
      IMPLICIT NONE
      EXTERNAL F , JAC , PSOL , VNLS
      DOUBLE PRECISION Y , Yh , Yh1 , Ewt , Savf , Vsav , Acor , Wm ,   &
                       Rpar
      INTEGER Ldyh , Iwm , Ipar
      DIMENSION Y(*) , Yh(Ldyh,*) , Yh1(*) , Ewt(*) , Savf(*) , Vsav(*) &
                , Acor(*) , Wm(*) , Iwm(*) , Rpar(*) , Ipar(*)
!-----------------------------------------------------------------------
! Call sequence input -- Y, YH, LDYH, YH1, EWT, SAVF, VSAV,
!                        ACOR, WM, IWM, F, JAC, PSOL, VNLS, RPAR, IPAR
! Call sequence output -- YH, ACOR, WM, IWM
! COMMON block variables accessed:
!     /DVOD01/  ACNRM, EL(13), H, HMIN, HMXI, HNEW, HSCAL, RC, TAU(13),
!               TQ(5), TN, JCUR, JSTART, KFLAG, KUTH,
!               L, LMAX, MAXORD, N, NEWQ, NQ, NQWAIT
!     /DVOD02/  HU, NCFN, NETF, NFE, NQU, NST
!
! Subroutines called by DVSTEP: F, DAXPY, DCOPY, DSCAL,
!                               DVJUST, VNLS, DVSET
! Function routines called by DVSTEP: DVNORM
!-----------------------------------------------------------------------
! DVSTEP performs one step of the integration of an initial value
! problem for a system of ordinary differential equations.
! DVSTEP calls subroutine VNLS for the solution of the nonlinear system
! arising in the time step.  Thus it is independent of the problem
! Jacobian structure and the type of nonlinear system solution method.
! DVSTEP returns a completion flag KFLAG (in COMMON).
! A return with KFLAG = -1 or -2 means either ABS(H) = HMIN or 10
! consecutive failures occurred.  On a return with KFLAG negative,
! the values of TN and the YH array are as of the beginning of the last
! step, and H is the last step size attempted.
!
! Communication with DVSTEP is done with the following variables:
!
! Y      = An array of length N used for the dependent variable vector.
! YH     = An LDYH by LMAX array containing the dependent variables
!          and their approximate scaled derivatives, where
!          LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate
!          j-th derivative of y(i), scaled by H**j/factorial(j)
!          (j = 0,1,...,NQ).  On entry for the first step, the first
!          two columns of YH must be set from the initial values.
! LDYH   = A constant integer .ge. N, the first dimension of YH.
!          N is the number of ODEs in the system.
! YH1    = A one-dimensional array occupying the same space as YH.
! EWT    = An array of length N containing multiplicative weights
!          for local error measurements.  Local errors in y(i) are
!          compared to 1.0/EWT(i) in various error tests.
! SAVF   = An array of working storage, of length N.
!          also used for input of YH(*,MAXORD+2) when JSTART = -1
!          and MAXORD .lt. the current order NQ.
! VSAV   = A work array of length N passed to subroutine VNLS.
! ACOR   = A work array of length N, used for the accumulated
!          corrections.  On a successful return, ACOR(i) contains
!          the estimated one-step local error in y(i).
! WM,IWM = Real and integer work arrays associated with matrix
!          operations in VNLS.
! F      = Dummy name for the user supplied subroutine for f.
! JAC    = Dummy name for the user supplied Jacobian subroutine.
! PSOL   = Dummy name for the subroutine passed to VNLS, for
!          possible use there.
! VNLS   = Dummy name for the nonlinear system solving subroutine,
!          whose real name is dependent on the method used.
! RPAR, IPAR = Dummy names for user's real and integer work arrays.
!-----------------------------------------------------------------------
!
! Type declarations for labeled COMMON block DVOD01 --------------------
!
      DOUBLE PRECISION ACNrm , CCMxj , CONp , CRAte , DRC , EL , ETA ,  &
                       ETAmax , H , HMIn , HMXi , HNEw , HSCal , PRL1 , &
                       RC , RL1 , TAU , TQ , TN , UROund
      INTEGER ICF , INIt , IPUp , JCUr , JSTart , JSV , KFLag , KUTh ,  &
              L , LMAx , LYH , LEWt , LACor , LSAvf , LWM , LIWm ,      &
              LOCjs , MAXord , METh , MITer , MSBj , MXHnil , MXStep ,  &
              N , NEWh , NEWq , NHNil , NQ , NQNyh , NQWait , NSLj ,    &
              NSLp , NYH
!
! Type declarations for labeled COMMON block DVOD02 --------------------
!
      DOUBLE PRECISION HU
      INTEGER NCFn , NETf , NFE , NJE , NLU , NNI , NQU , NST
!
! Type declarations for local variables --------------------------------
!
      DOUBLE PRECISION addon , bias1 , bias2 , bias3 , cnquot , ddn ,   &
                       dsm , dup , etacf , etamin , etamx1 , etamx2 ,   &
                       etamx3 , etamxf , etaq , etaqm1 , etaqp1 ,       &
                       flotl , one , onepsm , r , thresh , told , zero
      INTEGER i , i1 , i2 , iback , j , jb , kfc , kfh , mxncf , ncf ,  &
              nflag
!
! Type declaration for function subroutines called ---------------------
!
!      DOUBLE PRECISION DVNORM
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
      SAVE addon , bias1 , bias2 , bias3 , etacf , etamin , etamx1 ,    &
         etamx2 , etamx3 , etamxf , etaq , etaqm1 , kfc , kfh , mxncf , &
         onepsm , thresh , one , zero
!-----------------------------------------------------------------------
      COMMON /DVOD01/ ACNrm , CCMxj , CONp , CRAte , DRC , EL(13) ,     &
                      ETA , ETAmax , H , HMIn , HMXi , HNEw , HSCal ,   &
                      PRL1 , RC , RL1 , TAU(13) , TQ(5) , TN , UROund , &
                      ICF , INIt , IPUp , JCUr , JSTart , JSV , KFLag , &
                      KUTh , L , LMAx , LYH , LEWt , LACor , LSAvf ,    &
                      LWM , LIWm , LOCjs , MAXord , METh , MITer ,      &
                      MSBj , MXHnil , MXStep , N , NEWh , NEWq , NHNil ,&
                      NQ , NQNyh , NQWait , NSLj , NSLp , NYH
      COMMON /DVOD02/ HU , NCFn , NETf , NFE , NJE , NLU , NNI , NQU ,  &
                      NST
!
      DATA kfc/ - 3/ , kfh/ - 7/ , mxncf/10/
      DATA addon/1.0D-6/ , bias1/6.0D0/ , bias2/6.0D0/ , bias3/10.0D0/ ,&
           etacf/0.25D0/ , etamin/0.1D0/ , etamxf/0.2D0/ ,              &
           etamx1/1.0D4/ , etamx2/10.0D0/ , etamx3/10.0D0/ ,            &
           onepsm/1.00001D0/ , thresh/1.5D0/
      DATA one/1.0D0/ , zero/0.0D0/
!
      KFLag = 0
      told = TN
      ncf = 0
      JCUr = 0
      nflag = 0
      IF ( JSTart<=0 ) THEN
         IF ( JSTart==-1 ) GOTO 200
!-----------------------------------------------------------------------
! On the first call, the order is set to 1, and other variables are
! initialized.  ETAMAX is the maximum ratio by which H can be increased
! in a single step.  It is normally 10, but is larger during the
! first step to compensate for the small initial H.  If a failure
! occurs (in corrector convergence or error test), ETAMAX is set to 1
! for the next increase.
!-----------------------------------------------------------------------
         LMAx = MAXord + 1
         NQ = 1
         L = 2
         NQNyh = NQ*Ldyh
         TAU(1) = H
         PRL1 = one
         RC = zero
         ETAmax = etamx1
         NQWait = 2
         HSCal = H
         GOTO 400
!-----------------------------------------------------------------------
! Take preliminary actions on a normal continuation step (JSTART.GT.0).
! If the driver changed H, then ETA must be reset and NEWH set to 1.
! If a change of order was dictated on the previous step, then
! it is done here and appropriate adjustments in the history are made.
! On an order decrease, the history array is adjusted by DVJUST.
! On an order increase, the history array is augmented by a column.
! On a change of step size H, the history array YH is rescaled.
!-----------------------------------------------------------------------
      ELSEIF ( KUTh==1 ) THEN
         ETA = MIN(ETA,H/HSCal)
         NEWh = 1
      ENDIF
 100  IF ( NEWh==0 ) GOTO 400
      IF ( NEWq==NQ ) GOTO 300
      IF ( NEWq<NQ ) THEN
         CALL DVJUST(Yh,Ldyh,-1)
         NQ = NEWq
         L = NQ + 1
         NQWait = L
         GOTO 300
      ENDIF
      IF ( NEWq>NQ ) THEN
         CALL DVJUST(Yh,Ldyh,1)
         NQ = NEWq
         L = NQ + 1
         NQWait = L
         GOTO 300
      ENDIF
!-----------------------------------------------------------------------
! The following block handles preliminaries needed when JSTART = -1.
! If N was reduced, zero out part of YH to avoid undefined references.
! If MAXORD was reduced to a value less than the tentative order NEWQ,
! then NQ is set to MAXORD, and a new H ratio ETA is chosen.
! Otherwise, we take the same preliminary actions as for JSTART .gt. 0.
! In any case, NQWAIT is reset to L = NQ + 1 to prevent further
! changes in order for that many steps.
! The new H ratio ETA is limited by the input H if KUTH = 1,
! by HMIN if KUTH = 0, and by HMXI in any case.
! Finally, the history array YH is rescaled.
!-----------------------------------------------------------------------
 200  LMAx = MAXord + 1
      IF ( N/=Ldyh ) THEN
         i1 = 1 + (NEWq+1)*Ldyh
         i2 = (MAXord+1)*Ldyh
         IF ( i1<=i2 ) THEN
            DO i = i1 , i2
               Yh1(i) = zero
            ENDDO
         ENDIF
      ENDIF
      IF ( NEWq>MAXord ) THEN
         flotl = REAL(LMAx)
         IF ( MAXord<NQ-1 ) THEN
            ddn = DVNORM(N,Savf,Ewt)/TQ(1)
            ETA = one/((bias1*ddn)**(one/flotl)+addon)
         ENDIF
         IF ( MAXord==NQ .AND. NEWq==NQ+1 ) ETA = etaq
         IF ( MAXord==NQ-1 .AND. NEWq==NQ+1 ) THEN
            ETA = etaqm1
            CALL DVJUST(Yh,Ldyh,-1)
         ENDIF
         IF ( MAXord==NQ-1 .AND. NEWq==NQ ) THEN
            ddn = DVNORM(N,Savf,Ewt)/TQ(1)
            ETA = one/((bias1*ddn)**(one/flotl)+addon)
            CALL DVJUST(Yh,Ldyh,-1)
         ENDIF
         ETA = MIN(ETA,one)
         NQ = MAXord
         L = LMAx
      ENDIF
      IF ( KUTh==1 ) ETA = MIN(ETA,ABS(H/HSCal))
      IF ( KUTh==0 ) ETA = MAX(ETA,HMIn/ABS(HSCal))
      ETA = ETA/MAX(one,ABS(HSCal)*HMXi*ETA)
      NEWh = 1
      NQWait = L
      IF ( NEWq<=MAXord ) GOTO 100
! Rescale the history array for a change in H by a factor of ETA. ------
 300  r = one
      DO j = 2 , L
         r = r*ETA
         CALL DSCAL(N,r,Yh(1,j),1)
      ENDDO
      H = HSCal*ETA
      HSCal = H
      RC = RC*ETA
      NQNyh = NQ*Ldyh
!-----------------------------------------------------------------------
! This section computes the predicted values by effectively
! multiplying the YH array by the Pascal triangle matrix.
! DVSET is called to calculate all integration coefficients.
! RC is the ratio of new to old values of the coefficient H/EL(2)=h/l1.
!-----------------------------------------------------------------------
 400  TN = TN + H
      i1 = NQNyh + 1
      DO jb = 1 , NQ
         i1 = i1 - Ldyh
         DO i = i1 , NQNyh
            Yh1(i) = Yh1(i) + Yh1(i+Ldyh)
         ENDDO
      ENDDO
      CALL DVSET
      RL1 = one/EL(2)
      RC = RC*(RL1/PRL1)
      PRL1 = RL1
!
! Call the nonlinear system solver. ------------------------------------
!
      CALL VNLS(Y,Yh,Ldyh,Vsav,Savf,Ewt,Acor,Iwm,Wm,F,JAC,PSOL,nflag,   &
                Rpar,Ipar)
!
      IF ( nflag==0 ) THEN
!-----------------------------------------------------------------------
! The corrector has converged (NFLAG = 0).  The local error test is
! made and control passes to statement 500 if it fails.
!-----------------------------------------------------------------------
         dsm = ACNrm/TQ(2)
         IF ( dsm>one ) THEN
!-----------------------------------------------------------------------
! The error test failed.  KFLAG keeps track of multiple failures.
! Restore TN and the YH array to their previous values, and prepare
! to try the step again.  Compute the optimum step size for the
! same order.  After repeated failures, H is forced to decrease
! more rapidly.
!-----------------------------------------------------------------------
            KFLag = KFLag - 1
            NETf = NETf + 1
            nflag = -2
            TN = told
            i1 = NQNyh + 1
            DO jb = 1 , NQ
               i1 = i1 - Ldyh
               DO i = i1 , NQNyh
                  Yh1(i) = Yh1(i) - Yh1(i+Ldyh)
               ENDDO
            ENDDO
            IF ( ABS(H)<=HMIn*onepsm ) THEN
!-----------------------------------------------------------------------
! All returns are made through this section.
! On a successful return, ETAMAX is reset and ACOR is scaled.
!-----------------------------------------------------------------------
               KFLag = -1
               GOTO 600
            ELSE
               ETAmax = one
               IF ( KFLag>kfc ) THEN
! Compute ratio of new H to current H at the current order. ------------
                  flotl = REAL(L)
                  ETA = one/((bias2*dsm)**(one/flotl)+addon)
                  ETA = MAX(ETA,HMIn/ABS(H),etamin)
                  IF ( (KFLag<=-2) .AND. (ETA>etamxf) ) ETA = etamxf
                  GOTO 300
!-----------------------------------------------------------------------
! Control reaches this section if 3 or more consecutive failures
! have occurred.  It is assumed that the elements of the YH array
! have accumulated errors of the wrong order.  The order is reduced
! by one, if possible.  Then H is reduced by a factor of 0.1 and
! the step is retried.  After a total of 7 consecutive failures,
! an exit is taken with KFLAG = -1.
!-----------------------------------------------------------------------
               ELSEIF ( KFLag==kfh ) THEN
                  KFLag = -1
                  GOTO 600
               ELSEIF ( NQ==1 ) THEN
                  ETA = MAX(etamin,HMIn/ABS(H))
                  H = H*ETA
                  HSCal = H
                  TAU(1) = H
                  CALL F(N,TN,Y,Savf,Rpar,Ipar)
                  NFE = NFE + 1
                  DO i = 1 , N
                     Yh(i,2) = H*Savf(i)
                  ENDDO
                  NQWait = 10
                  GOTO 400
               ELSE
                  ETA = MAX(etamin,HMIn/ABS(H))
                  CALL DVJUST(Yh,Ldyh,-1)
                  L = NQ
                  NQ = NQ - 1
                  NQWait = L
                  GOTO 300
               ENDIF
            ENDIF
         ELSE
!-----------------------------------------------------------------------
! After a successful step, update the YH and TAU arrays and decrement
! NQWAIT.  If NQWAIT is then 1 and NQ .lt. MAXORD, then ACOR is saved
! for use in a possible order increase on the next step.
! If ETAMAX = 1 (a failure occurred this step), keep NQWAIT .ge. 2.
!-----------------------------------------------------------------------
            KFLag = 0
            NST = NST + 1
            HU = H
            NQU = NQ
            DO iback = 1 , NQ
               i = L - iback
               TAU(i+1) = TAU(i)
            ENDDO
            TAU(1) = H
            DO j = 1 , L
               CALL DAXPY(N,EL(j),Acor,1,Yh(1,j),1)
            ENDDO
            NQWait = NQWait - 1
            IF ( (L/=LMAx) .AND. (NQWait==1) ) THEN
               CALL DCOPY(N,Acor,1,Yh(1,LMAx),1)
               CONp = TQ(5)
            ENDIF
            IF ( ETAmax/=one ) THEN
!-----------------------------------------------------------------------
! If NQWAIT = 0, an increase or decrease in order by one is considered.
! Factors ETAQ, ETAQM1, ETAQP1 are computed by which H could
! be multiplied at order q, q-1, or q+1, respectively.
! The largest of these is determined, and the new order and
! step size set accordingly.
! A change of H or NQ is made only if H increases by at least a
! factor of THRESH.  If an order change is considered and rejected,
! then NQWAIT is set to 2 (reconsider it after 2 steps).
!-----------------------------------------------------------------------
! Compute ratio of new H to current H at the current order. ------------
               flotl = REAL(L)
               etaq = one/((bias2*dsm)**(one/flotl)+addon)
               IF ( NQWait==0 ) THEN
                  NQWait = 2
                  etaqm1 = zero
                  IF ( NQ/=1 ) THEN
! Compute ratio of new H to current H at the current order less one. ---
                     ddn = DVNORM(N,Yh(1,L),Ewt)/TQ(1)
                     etaqm1 = one/((bias1*ddn)**(one/(flotl-one))+addon)
                  ENDIF
                  etaqp1 = zero
                  IF ( L/=LMAx ) THEN
! Compute ratio of new H to current H at current order plus one. -------
                     cnquot = (TQ(5)/CONp)*(H/TAU(2))**L
                     DO i = 1 , N
                        Savf(i) = Acor(i) - cnquot*Yh(i,LMAx)
                     ENDDO
                     dup = DVNORM(N,Savf,Ewt)/TQ(3)
                     etaqp1 = one/((bias3*dup)**(one/(flotl+one))+addon)
                  ENDIF
                  IF ( etaq<etaqp1 ) THEN
                     IF ( etaqp1<=etaqm1 ) GOTO 420
                     ETA = etaqp1
                     NEWq = NQ + 1
                     CALL DCOPY(N,Acor,1,Yh(1,LMAx),1)
                     GOTO 450
                  ELSEIF ( etaq<etaqm1 ) THEN
                     GOTO 420
                  ENDIF
               ENDIF
               ETA = etaq
               NEWq = NQ
               GOTO 450
            ELSE
               IF ( NQWait<2 ) NQWait = 2
               NEWq = NQ
               NEWh = 0
               ETA = one
               HNEw = H
               GOTO 500
            ENDIF
 420        ETA = etaqm1
            NEWq = NQ - 1
         ENDIF
! Test tentative new H against THRESH, ETAMAX, and HMXI, then exit. ----
 450     IF ( ETA<thresh .OR. ETAmax==one ) THEN
            NEWq = NQ
            NEWh = 0
            ETA = one
            HNEw = H
         ELSE
            ETA = MIN(ETA,ETAmax)
            ETA = ETA/MAX(one,ABS(H)*HMXi*ETA)
            NEWh = 1
            HNEw = H*ETA
         ENDIF
      ELSE
!-----------------------------------------------------------------------
! The VNLS routine failed to achieve convergence (NFLAG .NE. 0).
! The YH array is retracted to its values before prediction.
! The step size H is reduced and the step is retried, if possible.
! Otherwise, an error exit is taken.
!-----------------------------------------------------------------------
         ncf = ncf + 1
         NCFn = NCFn + 1
         ETAmax = one
         TN = told
         i1 = NQNyh + 1
         DO jb = 1 , NQ
            i1 = i1 - Ldyh
            DO i = i1 , NQNyh
               Yh1(i) = Yh1(i) - Yh1(i+Ldyh)
            ENDDO
         ENDDO
         IF ( nflag<-1 ) THEN
            IF ( nflag==-2 ) KFLag = -3
            IF ( nflag==-3 ) KFLag = -4
            GOTO 600
         ELSEIF ( ABS(H)<=HMIn*onepsm ) THEN
            KFLag = -2
            GOTO 600
         ELSEIF ( ncf==mxncf ) THEN
            KFLag = -2
            GOTO 600
         ELSE
            ETA = etacf
            ETA = MAX(ETA,HMIn/ABS(H))
            nflag = -1
            GOTO 300
         ENDIF
      ENDIF
 500  ETAmax = etamx3
      IF ( NST<=10 ) ETAmax = etamx2
      r = one/TQ(2)
      CALL DSCAL(N,r,Acor,1)
 600  JSTart = 1
!----------------------- End of Subroutine DVSTEP ----------------------
      END Subroutine DVSTEP
!DECK DVSET
      SUBROUTINE DVSET()
      IMPLICIT NONE
!-----------------------------------------------------------------------
! Call sequence communication: None
! COMMON block variables accessed:
!     /DVOD01/ -- EL(13), H, TAU(13), TQ(5), L(= NQ + 1),
!                 METH, NQ, NQWAIT
!
! Subroutines called by DVSET: None
! Function routines called by DVSET: None
!-----------------------------------------------------------------------
! DVSET is called by DVSTEP and sets coefficients for use there.
!
! For each order NQ, the coefficients in EL are calculated by use of
!  the generating polynomial lambda(x), with coefficients EL(i).
!      lambda(x) = EL(1) + EL(2)*x + ... + EL(NQ+1)*(x**NQ).
! For the backward differentiation formulas,
!                                     NQ-1
!      lambda(x) = (1 + x/xi*(NQ)) * product (1 + x/xi(i) ) .
!                                     i = 1
! For the Adams formulas,
!                              NQ-1
!      (d/dx) lambda(x) = c * product (1 + x/xi(i) ) ,
!                              i = 1
!      lambda(-1) = 0,    lambda(0) = 1,
! where c is a normalization constant.
! In both cases, xi(i) is defined by
!      H*xi(i) = t sub n  -  t sub (n-i)
!              = H + TAU(1) + TAU(2) + ... TAU(i-1).
!
!
! In addition to variables described previously, communication
! with DVSET uses the following:
!   TAU    = A vector of length 13 containing the past NQ values
!            of H.
!   EL     = A vector of length 13 in which vset stores the
!            coefficients for the corrector formula.
!   TQ     = A vector of length 5 in which vset stores constants
!            used for the convergence test, the error test, and the
!            selection of H at a new order.
!   METH   = The basic method indicator.
!   NQ     = The current order.
!   L      = NQ + 1, the length of the vector stored in EL, and
!            the number of columns of the YH array being used.
!   NQWAIT = A counter controlling the frequency of order changes.
!            An order change is about to be considered if NQWAIT = 1.
!-----------------------------------------------------------------------
!
! Type declarations for labeled COMMON block DVOD01 --------------------
!
      DOUBLE PRECISION ACNrm , CCMxj , CONp , CRAte , DRC , EL , ETA ,  &
                       ETAmax , H , HMIn , HMXi , HNEw , HSCal , PRL1 , &
                       RC , RL1 , TAU , TQ , TN , UROund
      INTEGER ICF , INIt , IPUp , JCUr , JSTart , JSV , KFLag , KUTh ,  &
              L , LMAx , LYH , LEWt , LACor , LSAvf , LWM , LIWm ,      &
              LOCjs , MAXord , METh , MITer , MSBj , MXHnil , MXStep ,  &
              N , NEWh , NEWq , NHNil , NQ , NQNyh , NQWait , NSLj ,    &
              NSLp , NYH
!
! Type declarations for local variables --------------------------------
!
      DOUBLE PRECISION ahatn0 , alph0 , cnqm1 , cortes , csum , elp ,   &
                       em , em0 , floti , flotl , flotnq , hsum , one , &
                       rxi , rxis , s , six , t1 , t2 , t3 , t4 , t5 ,  &
                       t6 , two , xi , zero
      INTEGER i , iback , j , jp1 , nqm1 , nqm2
!
      DIMENSION em(13)
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
      SAVE cortes , one , six , two , zero
!
      COMMON /DVOD01/ ACNrm , CCMxj , CONp , CRAte , DRC , EL(13) ,     &
                      ETA , ETAmax , H , HMIn , HMXi , HNEw , HSCal ,   &
                      PRL1 , RC , RL1 , TAU(13) , TQ(5) , TN , UROund , &
                      ICF , INIt , IPUp , JCUr , JSTart , JSV , KFLag , &
                      KUTh , L , LMAx , LYH , LEWt , LACor , LSAvf ,    &
                      LWM , LIWm , LOCjs , MAXord , METh , MITer ,      &
                      MSBj , MXHnil , MXStep , N , NEWh , NEWq , NHNil ,&
                      NQ , NQNyh , NQWait , NSLj , NSLp , NYH
!
      DATA cortes/0.1D0/
      DATA one/1.0D0/ , six/6.0D0/ , two/2.0D0/ , zero/0.0D0/
!
      flotl = REAL(L)
      nqm1 = NQ - 1
      nqm2 = NQ - 2
      IF ( METh==2 ) THEN
!
! Set coefficients for BDF methods. ------------------------------------
         DO i = 3 , L
            EL(i) = zero
         ENDDO
         EL(1) = one
         EL(2) = one
         alph0 = -one
         ahatn0 = -one
         hsum = H
         rxi = one
         rxis = one
         IF ( NQ/=1 ) THEN
            DO j = 1 , nqm2
! In EL, construct coefficients of (1+x/xi(1))*...*(1+x/xi(j+1)). ------
               hsum = hsum + TAU(j)
               rxi = H/hsum
               jp1 = j + 1
               alph0 = alph0 - one/REAL(jp1)
               DO iback = 1 , jp1
                  i = (j+3) - iback
                  EL(i) = EL(i) + EL(i-1)*rxi
               ENDDO
            ENDDO
            alph0 = alph0 - one/REAL(NQ)
            rxis = -EL(2) - alph0
            hsum = hsum + TAU(nqm1)
            rxi = H/hsum
            ahatn0 = -EL(2) - rxi
            DO iback = 1 , NQ
               i = (NQ+2) - iback
               EL(i) = EL(i) + EL(i-1)*rxis
            ENDDO
         ENDIF
         t1 = one - ahatn0 + alph0
         t2 = one + REAL(NQ)*t1
         TQ(2) = ABS(alph0*t2/t1)
         TQ(5) = ABS(t2/(EL(L)*rxi/rxis))
         IF ( NQWait==1 ) THEN
            cnqm1 = rxis/EL(L)
            t3 = alph0 + one/REAL(NQ)
            t4 = ahatn0 + rxi
            elp = t3/(one-t4+t3)
            TQ(1) = ABS(elp/cnqm1)
            hsum = hsum + TAU(NQ)
            rxi = H/hsum
            t5 = alph0 - one/REAL(NQ+1)
            t6 = ahatn0 - rxi
            elp = t2/(one-t6+t5)
            TQ(3) = ABS(elp*rxi*(flotl+one)*t5)
         ENDIF
!
! Set coefficients for Adams methods. ----------------------------------
      ELSEIF ( NQ/=1 ) THEN
         hsum = H
         em(1) = one
         flotnq = flotl - one
         DO i = 2 , L
            em(i) = zero
         ENDDO
         DO j = 1 , nqm1
            IF ( (j==nqm1) .AND. (NQWait==1) ) THEN
               s = one
               csum = zero
               DO i = 1 , nqm1
                  csum = csum + s*em(i)/REAL(i+1)
                  s = -s
               ENDDO
               TQ(1) = em(nqm1)/(flotnq*csum)
            ENDIF
            rxi = H/hsum
            DO iback = 1 , j
               i = (j+2) - iback
               em(i) = em(i) + em(i-1)*rxi
            ENDDO
            hsum = hsum + TAU(j)
         ENDDO
! Compute integral from -1 to 0 of polynomial and of x times it. -------
         s = one
         em0 = zero
         csum = zero
         DO i = 1 , NQ
            floti = REAL(i)
            em0 = em0 + s*em(i)/floti
            csum = csum + s*em(i)/(floti+one)
            s = -s
         ENDDO
! In EL, form coefficients of normalized integrated polynomial. --------
         s = one/em0
         EL(1) = one
         DO i = 1 , NQ
            EL(i+1) = s*em(i)/REAL(i)
         ENDDO
         xi = hsum/H
         TQ(2) = xi*em0/csum
         TQ(5) = xi/EL(L)
         IF ( NQWait==1 ) THEN
! For higher order control constant, multiply polynomial by 1+x/xi(q). -
            rxi = one/xi
            DO iback = 1 , NQ
               i = (L+1) - iback
               em(i) = em(i) + em(i-1)*rxi
            ENDDO
! Compute integral of polynomial. --------------------------------------
            s = one
            csum = zero
            DO i = 1 , L
               csum = csum + s*em(i)/REAL(i+1)
               s = -s
            ENDDO
            TQ(3) = flotl*em0/csum
         ENDIF
      ELSE
         EL(1) = one
         EL(2) = one
         TQ(1) = one
         TQ(2) = two
         TQ(3) = six*TQ(2)
         TQ(5) = one
      ENDIF
      TQ(4) = cortes*TQ(2)
!----------------------- End of Subroutine DVSET -----------------------
      END Subroutine DVSET
!DECK DVJUST
      SUBROUTINE DVJUST(Yh,Ldyh,Iord)
      IMPLICIT NONE
      DOUBLE PRECISION Yh
      INTEGER Ldyh , Iord
      DIMENSION Yh(Ldyh,*)
!-----------------------------------------------------------------------
! Call sequence input -- YH, LDYH, IORD
! Call sequence output -- YH
! COMMON block input -- NQ, METH, LMAX, HSCAL, TAU(13), N
! COMMON block variables accessed:
!     /DVOD01/ -- HSCAL, TAU(13), LMAX, METH, N, NQ,
!
! Subroutines called by DVJUST: DAXPY
! Function routines called by DVJUST: None
!-----------------------------------------------------------------------
! This subroutine adjusts the YH array on reduction of order,
! and also when the order is increased for the stiff option (METH = 2).
! Communication with DVJUST uses the following:
! IORD  = An integer flag used when METH = 2 to indicate an order
!         increase (IORD = +1) or an order decrease (IORD = -1).
! HSCAL = Step size H used in scaling of Nordsieck array YH.
!         (If IORD = +1, DVJUST assumes that HSCAL = TAU(1).)
! See References 1 and 2 for details.
!-----------------------------------------------------------------------
!
! Type declarations for labeled COMMON block DVOD01 --------------------
!
      DOUBLE PRECISION ACNrm , CCMxj , CONp , CRAte , DRC , EL , ETA ,  &
                       ETAmax , H , HMIn , HMXi , HNEw , HSCal , PRL1 , &
                       RC , RL1 , TAU , TQ , TN , UROund
      INTEGER ICF , INIt , IPUp , JCUr , JSTart , JSV , KFLag , KUTh ,  &
              L , LMAx , LYH , LEWt , LACor , LSAvf , LWM , LIWm ,      &
              LOCjs , MAXord , METh , MITer , MSBj , MXHnil , MXStep ,  &
              N , NEWh , NEWq , NHNil , NQ , NQNyh , NQWait , NSLj ,    &
              NSLp , NYH
!
! Type declarations for local variables --------------------------------
!
      DOUBLE PRECISION alph0 , alph1 , hsum , one , prod , t1 , xi ,    &
                       xiold , zero
      INTEGER i , iback , j , jp1 , lp1 , nqm1 , nqm2 , nqp1
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
      SAVE one , zero
!
      COMMON /DVOD01/ ACNrm , CCMxj , CONp , CRAte , DRC , EL(13) ,     &
                      ETA , ETAmax , H , HMIn , HMXi , HNEw , HSCal ,   &
                      PRL1 , RC , RL1 , TAU(13) , TQ(5) , TN , UROund , &
                      ICF , INIt , IPUp , JCUr , JSTart , JSV , KFLag , &
                      KUTh , L , LMAx , LYH , LEWt , LACor , LSAvf ,    &
                      LWM , LIWm , LOCjs , MAXord , METh , MITer ,      &
                      MSBj , MXHnil , MXStep , N , NEWh , NEWq , NHNil ,&
                      NQ , NQNyh , NQWait , NSLj , NSLp , NYH
!
      DATA one/1.0D0/ , zero/0.0D0/
!
      IF ( (NQ==2) .AND. (Iord/=1) ) RETURN
      nqm1 = NQ - 1
      nqm2 = NQ - 2
      IF ( METh==2 ) THEN
!-----------------------------------------------------------------------
! Stiff option...
! Check to see if the order is being increased or decreased.
!-----------------------------------------------------------------------
         IF ( Iord==1 ) THEN
! Order increase. ------------------------------------------------------
            DO j = 1 , LMAx
               EL(j) = zero
            ENDDO
            EL(3) = one
            alph0 = -one
            alph1 = one
            prod = one
            xiold = one
            hsum = HSCal
            IF ( NQ/=1 ) THEN
               DO j = 1 , nqm1
! Construct coefficients of x*x*(x+xi(1))*...*(x+xi(j)). ---------------
                  jp1 = j + 1
                  hsum = hsum + TAU(jp1)
                  xi = hsum/HSCal
                  prod = prod*xi
                  alph0 = alph0 - one/REAL(jp1)
                  alph1 = alph1 + one/xi
                  DO iback = 1 , jp1
                     i = (j+4) - iback
                     EL(i) = EL(i)*xiold + EL(i-1)
                  ENDDO
                  xiold = xi
               ENDDO
            ENDIF
            t1 = (-alph0-alph1)/prod
! Load column L + 1 in YH array. ---------------------------------------
            lp1 = L + 1
            DO i = 1 , N
               Yh(i,lp1) = t1*Yh(i,LMAx)
            ENDDO
! Add correction terms to YH array. ------------------------------------
            nqp1 = NQ + 1
            DO j = 3 , nqp1
               CALL DAXPY(N,EL(j),Yh(1,lp1),1,Yh(1,j),1)
            ENDDO
         ELSE
! Order decrease. ------------------------------------------------------
            DO j = 1 , LMAx
               EL(j) = zero
            ENDDO
            EL(3) = one
            hsum = zero
            DO j = 1 , nqm2
! Construct coefficients of x*x*(x+xi(1))*...*(x+xi(j)). ---------------
               hsum = hsum + TAU(j)
               xi = hsum/HSCal
               jp1 = j + 1
               DO iback = 1 , jp1
                  i = (j+4) - iback
                  EL(i) = EL(i)*xi + EL(i-1)
               ENDDO
            ENDDO
! Subtract correction terms from YH array. -----------------------------
            DO j = 3 , NQ
               DO i = 1 , N
                  Yh(i,j) = Yh(i,j) - Yh(i,L)*EL(j)
               ENDDO
            ENDDO
            RETURN
         ENDIF
!-----------------------------------------------------------------------
! Nonstiff option...
! Check to see if the order is being increased or decreased.
!-----------------------------------------------------------------------
      ELSEIF ( Iord==1 ) THEN
! Order increase. ------------------------------------------------------
! Zero out next column in YH array. ------------------------------------
         lp1 = L + 1
         DO i = 1 , N
            Yh(i,lp1) = zero
         ENDDO
         RETURN
      ELSE
! Order decrease. ------------------------------------------------------
         DO j = 1 , LMAx
            EL(j) = zero
         ENDDO
         EL(2) = one
         hsum = zero
         DO j = 1 , nqm2
! Construct coefficients of x*(x+xi(1))*...*(x+xi(j)). -----------------
            hsum = hsum + TAU(j)
            xi = hsum/HSCal
            jp1 = j + 1
            DO iback = 1 , jp1
               i = (j+3) - iback
               EL(i) = EL(i)*xi + EL(i-1)
            ENDDO
         ENDDO
! Construct coefficients of integrated polynomial. ---------------------
         DO j = 2 , nqm1
            EL(j+1) = REAL(NQ)*EL(j)/REAL(j)
         ENDDO
! Subtract correction terms from YH array. -----------------------------
         DO j = 3 , NQ
            DO i = 1 , N
               Yh(i,j) = Yh(i,j) - Yh(i,L)*EL(j)
            ENDDO
         ENDDO
         RETURN
      ENDIF
!----------------------- End of Subroutine DVJUST ----------------------
      END Subroutine DVJUST
!DECK DVNLSD
      SUBROUTINE DVNLSD(Y,Yh,Ldyh,Vsav,Savf,Ewt,Acor,Iwm,Wm,F,JAC,PDUM, &
                        Nflag,Rpar,Ipar)
      IMPLICIT NONE
      EXTERNAL F , JAC , PDUM
      DOUBLE PRECISION Y , Yh , Vsav , Savf , Ewt , Acor , Wm , Rpar
      INTEGER Ldyh , Iwm , Nflag , Ipar
      DIMENSION Y(*) , Yh(Ldyh,*) , Vsav(*) , Savf(*) , Ewt(*) , Acor(*)&
                , Iwm(*) , Wm(*) , Rpar(*) , Ipar(*)
!-----------------------------------------------------------------------
! Call sequence input -- Y, YH, LDYH, SAVF, EWT, ACOR, IWM, WM,
!                        F, JAC, NFLAG, RPAR, IPAR
! Call sequence output -- YH, ACOR, WM, IWM, NFLAG
! COMMON block variables accessed:
!     /DVOD01/ ACNRM, CRATE, DRC, H, RC, RL1, TQ(5), TN, ICF,
!                JCUR, METH, MITER, N, NSLP
!     /DVOD02/ HU, NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
!
! Subroutines called by DVNLSD: F, DAXPY, DCOPY, DSCAL, DVJAC, DVSOL
! Function routines called by DVNLSD: DVNORM
!-----------------------------------------------------------------------
! Subroutine DVNLSD is a nonlinear system solver, which uses functional
! iteration or a chord (modified Newton) method.  For the chord method
! direct linear algebraic system solvers are used.  Subroutine DVNLSD
! then handles the corrector phase of this integration package.
!
! Communication with DVNLSD is done with the following variables. (For
! more details, please see the comments in the driver subroutine.)
!
! Y          = The dependent variable, a vector of length N, input.
! YH         = The Nordsieck (Taylor) array, LDYH by LMAX, input
!              and output.  On input, it contains predicted values.
! LDYH       = A constant .ge. N, the first dimension of YH, input.
! VSAV       = Unused work array.
! SAVF       = A work array of length N.
! EWT        = An error weight vector of length N, input.
! ACOR       = A work array of length N, used for the accumulated
!              corrections to the predicted y vector.
! WM,IWM     = Real and integer work arrays associated with matrix
!              operations in chord iteration (MITER .ne. 0).
! F          = Dummy name for user supplied routine for f.
! JAC        = Dummy name for user supplied Jacobian routine.
! PDUM       = Unused dummy subroutine name.  Included for uniformity
!              over collection of integrators.
! NFLAG      = Input/output flag, with values and meanings as follows:
!              INPUT
!                  0 first call for this time step.
!                 -1 convergence failure in previous call to DVNLSD.
!                 -2 error test failure in DVSTEP.
!              OUTPUT
!                  0 successful completion of nonlinear solver.
!                 -1 convergence failure or singular matrix.
!                 -2 unrecoverable error in matrix preprocessing
!                    (cannot occur here).
!                 -3 unrecoverable error in solution (cannot occur
!                    here).
! RPAR, IPAR = Dummy names for user's real and integer work arrays.
!
! IPUP       = Own variable flag with values and meanings as follows:
!              0,            do not update the Newton matrix.
!              MITER .ne. 0, update Newton matrix, because it is the
!                            initial step, order was changed, the error
!                            test failed, or an update is indicated by
!                            the scalar RC or step counter NST.
!
! For more details, see comments in driver subroutine.
!-----------------------------------------------------------------------
! Type declarations for labeled COMMON block DVOD01 --------------------
!
      DOUBLE PRECISION ACNrm , CCMxj , CONp , CRAte , DRC , EL , ETA ,  &
                       ETAmax , H , HMIn , HMXi , HNEw , HSCal , PRL1 , &
                       RC , RL1 , TAU , TQ , TN , UROund
      INTEGER ICF , INIt , IPUp , JCUr , JSTart , JSV , KFLag , KUTh ,  &
              L , LMAx , LYH , LEWt , LACor , LSAvf , LWM , LIWm ,      &
              LOCjs , MAXord , METh , MITer , MSBj , MXHnil , MXStep ,  &
              N , NEWh , NEWq , NHNil , NQ , NQNyh , NQWait , NSLj ,    &
              NSLp , NYH
!
! Type declarations for labeled COMMON block DVOD02 --------------------
!
      DOUBLE PRECISION HU
      INTEGER NCFn , NETf , NFE , NJE , NLU , NNI , NQU , NST
!
! Type declarations for local variables --------------------------------
!
      DOUBLE PRECISION ccmax , crdown , cscale , dcon , del , delp ,    &
                       one , rdiv , two , zero
      INTEGER i , ierpj , iersl , m , maxcor , msbp
!
! Type declaration for function subroutines called ---------------------
!
!      DOUBLE PRECISION DVNORM
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
      SAVE ccmax , crdown , maxcor , msbp , rdiv , one , two , zero
!
      COMMON /DVOD01/ ACNrm , CCMxj , CONp , CRAte , DRC , EL(13) ,     &
                      ETA , ETAmax , H , HMIn , HMXi , HNEw , HSCal ,   &
                      PRL1 , RC , RL1 , TAU(13) , TQ(5) , TN , UROund , &
                      ICF , INIt , IPUp , JCUr , JSTart , JSV , KFLag , &
                      KUTh , L , LMAx , LYH , LEWt , LACor , LSAvf ,    &
                      LWM , LIWm , LOCjs , MAXord , METh , MITer ,      &
                      MSBj , MXHnil , MXStep , N , NEWh , NEWq , NHNil ,&
                      NQ , NQNyh , NQWait , NSLj , NSLp , NYH
      COMMON /DVOD02/ HU , NCFn , NETf , NFE , NJE , NLU , NNI , NQU ,  &
                      NST
!
      DATA ccmax/0.3D0/ , crdown/0.3D0/ , maxcor/3/ , msbp/20/ ,        &
           rdiv/2.0D0/
      DATA one/1.0D0/ , two/2.0D0/ , zero/0.0D0/
!-----------------------------------------------------------------------
! On the first step, on a change of method order, or after a
! nonlinear convergence failure with NFLAG = -2, set IPUP = MITER
! to force a Jacobian update when MITER .ne. 0.
!-----------------------------------------------------------------------
      IF ( JSTart==0 ) NSLp = 0
      IF ( Nflag==0 ) ICF = 0
      IF ( Nflag==-2 ) IPUp = MITer
      IF ( (JSTart==0) .OR. (JSTart==-1) ) IPUp = MITer
! If this is functional iteration, set CRATE .eq. 1 and drop to 220
      IF ( MITer==0 ) THEN
         CRAte = one
         GOTO 100
      ENDIF
!-----------------------------------------------------------------------
! RC is the ratio of new to old values of the coefficient H/EL(2)=h/l1.
! When RC differs from 1 by more than CCMAX, IPUP is set to MITER
! to force DVJAC to be called, if a Jacobian is involved.
! In any case, DVJAC is called at least every MSBP steps.
!-----------------------------------------------------------------------
      DRC = ABS(RC-one)
      IF ( DRC>ccmax .OR. NST>=NSLp+msbp ) IPUp = MITer
!-----------------------------------------------------------------------
! Up to MAXCOR corrector iterations are taken.  A convergence test is
! made on the r.m.s. norm of each correction, weighted by the error
! weight vector EWT.  The sum of the corrections is accumulated in the
! vector ACOR(i).  The YH array is not altered in the corrector loop.
!-----------------------------------------------------------------------
 100  m = 0
      delp = zero
      CALL DCOPY(N,Yh(1,1),1,Y,1)
      CALL F(N,TN,Y,Savf,Rpar,Ipar)
      NFE = NFE + 1
      IF ( IPUp>0 ) THEN
!-----------------------------------------------------------------------
! If indicated, the matrix P = I - h*rl1*J is reevaluated and
! preprocessed before starting the corrector iteration.  IPUP is set
! to 0 as an indicator that this has been done.
!-----------------------------------------------------------------------
         CALL DVJAC(Y,Yh,Ldyh,Ewt,Acor,Savf,Wm,Iwm,F,JAC,ierpj,Rpar,    &
                    Ipar)
         IPUp = 0
         RC = one
         DRC = zero
         CRAte = one
         NSLp = NST
! If matrix is singular, take error return to force cut in step size. --
         IF ( ierpj/=0 ) GOTO 400
      ENDIF
      DO i = 1 , N
         Acor(i) = zero
      ENDDO
! This is a looping point for the corrector iteration. -----------------
 200  IF ( MITer/=0 ) THEN
!-----------------------------------------------------------------------
! In the case of the chord method, compute the corrector error,
! and solve the linear system with that as right-hand side and
! P as coefficient matrix.  The correction is scaled by the factor
! 2/(1+RC) to account for changes in h*rl1 since the last DVJAC call.
!-----------------------------------------------------------------------
         DO i = 1 , N
            Y(i) = (RL1*H)*Savf(i) - (RL1*Yh(i,2)+Acor(i))
         ENDDO
         CALL DVSOL(Wm,Iwm,Y,iersl)
         NNI = NNI + 1
         IF ( iersl>0 ) GOTO 300
         IF ( METh==2 .AND. RC/=one ) THEN
            cscale = two/(one+RC)
            CALL DSCAL(N,cscale,Y,1)
         ENDIF
         del = DVNORM(N,Y,Ewt)
         CALL DAXPY(N,one,Y,1,Acor,1)
         DO i = 1 , N
            Y(i) = Yh(i,1) + Acor(i)
         ENDDO
      ELSE
!-----------------------------------------------------------------------
! In the case of functional iteration, update Y directly from
! the result of the last function evaluation.
!-----------------------------------------------------------------------
         DO i = 1 , N
            Savf(i) = RL1*(H*Savf(i)-Yh(i,2))
         ENDDO
         DO i = 1 , N
            Y(i) = Savf(i) - Acor(i)
         ENDDO
         del = DVNORM(N,Y,Ewt)
         DO i = 1 , N
            Y(i) = Yh(i,1) + Savf(i)
         ENDDO
         CALL DCOPY(N,Savf,1,Acor,1)
      ENDIF
!-----------------------------------------------------------------------
! Test for convergence.  If M .gt. 0, an estimate of the convergence
! rate constant is stored in CRATE, and this is used in the test.
!-----------------------------------------------------------------------
      IF ( m/=0 ) CRAte = MAX(crdown*CRAte,del/delp)
      dcon = del*MIN(one,CRAte)/TQ(4)
      IF ( dcon<=one ) THEN
!
! Return for successful step. ------------------------------------------
         Nflag = 0
         JCUr = 0
         ICF = 0
         IF ( m==0 ) ACNrm = del
         IF ( m>0 ) ACNrm = DVNORM(N,Acor,Ewt)
         GOTO 99999
      ELSE
         m = m + 1
         IF ( m/=maxcor ) THEN
            IF ( m<2 .OR. del<=rdiv*delp ) THEN
               delp = del
               CALL F(N,TN,Y,Savf,Rpar,Ipar)
               NFE = NFE + 1
               GOTO 200
            ENDIF
         ENDIF
      ENDIF
!
 300  IF ( MITer/=0 .AND. JCUr/=1 ) THEN
         ICF = 1
         IPUp = MITer
         GOTO 100
      ENDIF
!
 400  Nflag = -1
      ICF = 2
      IPUp = MITer
      RETURN
!----------------------- End of Subroutine DVNLSD ----------------------
99999 END Subroutine DVNLSD
!DECK DVJAC
      SUBROUTINE DVJAC(Y,Yh,Ldyh,Ewt,Ftem,Savf,Wm,Iwm,F,JAC,Ierpj,Rpar, &
                       Ipar)
      IMPLICIT NONE
      EXTERNAL F , JAC
      DOUBLE PRECISION Y , Yh , Ewt , Ftem , Savf , Wm , Rpar
      INTEGER Ldyh , Iwm , Ierpj , Ipar
      DIMENSION Y(*) , Yh(Ldyh,*) , Ewt(*) , Ftem(*) , Savf(*) , Wm(*) ,&
                Iwm(*) , Rpar(*) , Ipar(*)
!-----------------------------------------------------------------------
! Call sequence input -- Y, YH, LDYH, EWT, FTEM, SAVF, WM, IWM,
!                        F, JAC, RPAR, IPAR
! Call sequence output -- WM, IWM, IERPJ
! COMMON block variables accessed:
!     /DVOD01/  CCMXJ, DRC, H, RL1, TN, UROUND, ICF, JCUR, LOCJS,
!               MITER, MSBJ, N, NSLJ
!     /DVOD02/  NFE, NST, NJE, NLU
!
! Subroutines called by DVJAC: F, JAC, DACOPY, DCOPY, DGBFA, DGEFA,
!                              DSCAL
! Function routines called by DVJAC: DVNORM
!-----------------------------------------------------------------------
! DVJAC is called by DVNLSD to compute and process the matrix
! P = I - h*rl1*J , where J is an approximation to the Jacobian.
! Here J is computed by the user-supplied routine JAC if
! MITER = 1 or 4, or by finite differencing if MITER = 2, 3, or 5.
! If MITER = 3, a diagonal approximation to J is used.
! If JSV = -1, J is computed from scratch in all cases.
! If JSV = 1 and MITER = 1, 2, 4, or 5, and if the saved value of J is
! considered acceptable, then P is constructed from the saved J.
! J is stored in wm and replaced by P.  If MITER .ne. 3, P is then
! subjected to LU decomposition in preparation for later solution
! of linear systems with P as coefficient matrix. This is done
! by DGEFA if MITER = 1 or 2, and by DGBFA if MITER = 4 or 5.
!
! Communication with DVJAC is done with the following variables.  (For
! more details, please see the comments in the driver subroutine.)
! Y          = Vector containing predicted values on entry.
! YH         = The Nordsieck array, an LDYH by LMAX array, input.
! LDYH       = A constant .ge. N, the first dimension of YH, input.
! EWT        = An error weight vector of length N.
! SAVF       = Array containing f evaluated at predicted y, input.
! WM         = Real work space for matrices.  In the output, it containS
!              the inverse diagonal matrix if MITER = 3 and the LU
!              decomposition of P if MITER is 1, 2 , 4, or 5.
!              Storage of matrix elements starts at WM(3).
!              Storage of the saved Jacobian starts at WM(LOCJS).
!              WM also contains the following matrix-related data:
!              WM(1) = SQRT(UROUND), used in numerical Jacobian step.
!              WM(2) = H*RL1, saved for later use if MITER = 3.
! IWM        = Integer work space containing pivot information,
!              starting at IWM(31), if MITER is 1, 2, 4, or 5.
!              IWM also contains band parameters ML = IWM(1) and
!              MU = IWM(2) if MITER is 4 or 5.
! F          = Dummy name for the user supplied subroutine for f.
! JAC        = Dummy name for the user supplied Jacobian subroutine.
! RPAR, IPAR = Dummy names for user's real and integer work arrays.
! RL1        = 1/EL(2) (input).
! IERPJ      = Output error flag,  = 0 if no trouble, 1 if the P
!              matrix is found to be singular.
! JCUR       = Output flag to indicate whether the Jacobian matrix
!              (or approximation) is now current.
!              JCUR = 0 means J is not current.
!              JCUR = 1 means J is current.
!-----------------------------------------------------------------------
!
! Type declarations for labeled COMMON block DVOD01 --------------------
!
      DOUBLE PRECISION ACNrm , CCMxj , CONp , CRAte , DRC , EL , ETA ,  &
                       ETAmax , H , HMIn , HMXi , HNEw , HSCal , PRL1 , &
                       RC , RL1 , TAU , TQ , TN , UROund
      INTEGER ICF , INIt , IPUp , JCUr , JSTart , JSV , KFLag , KUTh ,  &
              L , LMAx , LYH , LEWt , LACor , LSAvf , LWM , LIWm ,      &
              LOCjs , MAXord , METh , MITer , MSBj , MXHnil , MXStep ,  &
              N , NEWh , NEWq , NHNil , NQ , NQNyh , NQWait , NSLj ,    &
              NSLp , NYH
!
! Type declarations for labeled COMMON block DVOD02 --------------------
!
      DOUBLE PRECISION HU
      INTEGER NCFn , NETf , NFE , NJE , NLU , NNI , NQU , NST
!
! Type declarations for local variables --------------------------------
!
      DOUBLE PRECISION con , di , fac , hrl1 , one , pt1 , r , r0 ,     &
                       srur , thou , yi , yj , yjj , zero
      INTEGER i , i1 , i2 , ier , ii , j , j1 , jj , jok , lenp , mba , &
              mband , meb1 , meband , ml , ml3 , mu , np1
!
! Type declaration for function subroutines called ---------------------
!
!      DOUBLE PRECISION DVNORM
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this subroutine.
!-----------------------------------------------------------------------
      SAVE one , pt1 , thou , zero
!-----------------------------------------------------------------------
      COMMON /DVOD01/ ACNrm , CCMxj , CONp , CRAte , DRC , EL(13) ,     &
                      ETA , ETAmax , H , HMIn , HMXi , HNEw , HSCal ,   &
                      PRL1 , RC , RL1 , TAU(13) , TQ(5) , TN , UROund , &
                      ICF , INIt , IPUp , JCUr , JSTart , JSV , KFLag , &
                      KUTh , L , LMAx , LYH , LEWt , LACor , LSAvf ,    &
                      LWM , LIWm , LOCjs , MAXord , METh , MITer ,      &
                      MSBj , MXHnil , MXStep , N , NEWh , NEWq , NHNil ,&
                      NQ , NQNyh , NQWait , NSLj , NSLp , NYH
      COMMON /DVOD02/ HU , NCFn , NETf , NFE , NJE , NLU , NNI , NQU ,  &
                      NST
!
      DATA one/1.0D0/ , thou/1000.0D0/ , zero/0.0D0/ , pt1/0.1D0/
!
      Ierpj = 0
      hrl1 = H*RL1
! See whether J should be evaluated (JOK = -1) or not (JOK = 1). -------
      jok = JSV
      IF ( JSV==1 ) THEN
         IF ( NST==0 .OR. NST>NSLj+MSBj ) jok = -1
         IF ( ICF==1 .AND. DRC<CCMxj ) jok = -1
         IF ( ICF==2 ) jok = -1
      ENDIF
! End of setting JOK. --------------------------------------------------
!
      IF ( jok==-1 .AND. MITer==1 ) THEN
! If JOK = -1 and MITER = 1, call JAC to evaluate Jacobian. ------------
         NJE = NJE + 1
         NSLj = NST
         JCUr = 1
         lenp = N*N
         DO i = 1 , lenp
            Wm(i+2) = zero
         ENDDO
         CALL JAC(N,TN,Y,0,0,Wm(3),N,Rpar,Ipar)
         IF ( JSV==1 ) CALL DCOPY(lenp,Wm(3),1,Wm(LOCjs),1)
      ENDIF
!
      IF ( jok==-1 .AND. MITer==2 ) THEN
! If MITER = 2, make N calls to F to approximate the Jacobian. ---------
         NJE = NJE + 1
         NSLj = NST
         JCUr = 1
         fac = DVNORM(N,Savf,Ewt)
         r0 = thou*ABS(H)*UROund*REAL(N)*fac
         IF ( r0==zero ) r0 = one
         srur = Wm(1)
         j1 = 2
         DO j = 1 , N
            yj = Y(j)
            r = MAX(srur*ABS(yj),r0/Ewt(j))
            Y(j) = Y(j) + r
            fac = one/r
            CALL F(N,TN,Y,Ftem,Rpar,Ipar)
            DO i = 1 , N
               Wm(i+j1) = (Ftem(i)-Savf(i))*fac
            ENDDO
            Y(j) = yj
            j1 = j1 + N
         ENDDO
         NFE = NFE + N
         lenp = N*N
         IF ( JSV==1 ) CALL DCOPY(lenp,Wm(3),1,Wm(LOCjs),1)
      ENDIF
!
      IF ( jok==1 .AND. (MITer==1 .OR. MITer==2) ) THEN
         JCUr = 0
         lenp = N*N
         CALL DCOPY(lenp,Wm(LOCjs),1,Wm(3),1)
      ENDIF
!
      IF ( MITer==1 .OR. MITer==2 ) THEN
! Multiply Jacobian by scalar, add identity, and do LU decomposition. --
         con = -hrl1
         CALL DSCAL(lenp,con,Wm(3),1)
         j = 3
         np1 = N + 1
         DO i = 1 , N
            Wm(j) = Wm(j) + one
            j = j + np1
         ENDDO
         NLU = NLU + 1
         CALL DGEFA(Wm(3),N,N,Iwm(31),ier)
         IF ( ier/=0 ) Ierpj = 1
         RETURN
      ENDIF
! End of code block for MITER = 1 or 2. --------------------------------
!
      IF ( MITer==3 ) THEN
! If MITER = 3, construct a diagonal approximation to J and P. ---------
         NJE = NJE + 1
         JCUr = 1
         Wm(2) = hrl1
         r = RL1*pt1
         DO i = 1 , N
            Y(i) = Y(i) + r*(H*Savf(i)-Yh(i,2))
         ENDDO
         CALL F(N,TN,Y,Wm(3),Rpar,Ipar)
         NFE = NFE + 1
         DO i = 1 , N
            r0 = H*Savf(i) - Yh(i,2)
            di = pt1*r0 - H*(Wm(i+2)-Savf(i))
            Wm(i+2) = one
            IF ( ABS(r0)>=UROund/Ewt(i) ) THEN
               IF ( ABS(di)==zero ) GOTO 50
               Wm(i+2) = pt1*r0/di
            ENDIF
         ENDDO
         RETURN
 50      Ierpj = 1
         RETURN
      ENDIF
! End of code block for MITER = 3. -------------------------------------
!
! Set constants for MITER = 4 or 5. ------------------------------------
      ml = Iwm(1)
      mu = Iwm(2)
      ml3 = ml + 3
      mband = ml + mu + 1
      meband = mband + ml
      lenp = meband*N
!
      IF ( jok==-1 .AND. MITer==4 ) THEN
! If JOK = -1 and MITER = 4, call JAC to evaluate Jacobian. ------------
         NJE = NJE + 1
         NSLj = NST
         JCUr = 1
         DO i = 1 , lenp
            Wm(i+2) = zero
         ENDDO
         CALL JAC(N,TN,Y,ml,mu,Wm(ml3),meband,Rpar,Ipar)
         IF ( JSV==1 ) CALL DACOPY(mband,N,Wm(ml3),meband,Wm(LOCjs),    &
                                   mband)
      ENDIF
!
      IF ( jok==-1 .AND. MITer==5 ) THEN
! If MITER = 5, make ML+MU+1 calls to F to approximate the Jacobian. ---
         NJE = NJE + 1
         NSLj = NST
         JCUr = 1
         mba = MIN(mband,N)
         meb1 = meband - 1
         srur = Wm(1)
         fac = DVNORM(N,Savf,Ewt)
         r0 = thou*ABS(H)*UROund*REAL(N)*fac
         IF ( r0==zero ) r0 = one
         DO j = 1 , mba
            DO i = j , N , mband
               yi = Y(i)
               r = MAX(srur*ABS(yi),r0/Ewt(i))
               Y(i) = Y(i) + r
            ENDDO
            CALL F(N,TN,Y,Ftem,Rpar,Ipar)
            DO jj = j , N , mband
               Y(jj) = Yh(jj,1)
               yjj = Y(jj)
               r = MAX(srur*ABS(yjj),r0/Ewt(jj))
               fac = one/r
               i1 = MAX(jj-mu,1)
               i2 = MIN(jj+ml,N)
               ii = jj*meb1 - ml + 2
               DO i = i1 , i2
                  Wm(ii+i) = (Ftem(i)-Savf(i))*fac
               ENDDO
            ENDDO
         ENDDO
         NFE = NFE + mba
         IF ( JSV==1 ) CALL DACOPY(mband,N,Wm(ml3),meband,Wm(LOCjs),    &
                                   mband)
      ENDIF
!
      IF ( jok==1 ) THEN
         JCUr = 0
         CALL DACOPY(mband,N,Wm(LOCjs),mband,Wm(ml3),meband)
      ENDIF
!
! Multiply Jacobian by scalar, add identity, and do LU decomposition.
      con = -hrl1
      CALL DSCAL(lenp,con,Wm(3),1)
      ii = mband + 2
      DO i = 1 , N
         Wm(ii) = Wm(ii) + one
         ii = ii + meband
      ENDDO
      NLU = NLU + 1
      CALL DGBFA(Wm(3),meband,N,ml,mu,Iwm(31),ier)
      IF ( ier/=0 ) Ierpj = 1
! End of code block for MITER = 4 or 5. --------------------------------
!
!----------------------- End of Subroutine DVJAC -----------------------
      END Subroutine DVJAC
!DECK DACOPY
      SUBROUTINE DACOPY(Nrow,Ncol,A,Nrowa,B,Nrowb)
      IMPLICIT NONE
      DOUBLE PRECISION A , B
      INTEGER Nrow , Ncol , Nrowa , Nrowb
      DIMENSION A(Nrowa,Ncol) , B(Nrowb,Ncol)
!-----------------------------------------------------------------------
! Call sequence input -- NROW, NCOL, A, NROWA, NROWB
! Call sequence output -- B
! COMMON block variables accessed -- None
!
! Subroutines called by DACOPY: DCOPY
! Function routines called by DACOPY: None
!-----------------------------------------------------------------------
! This routine copies one rectangular array, A, to another, B,
! where A and B may have different row dimensions, NROWA and NROWB.
! The data copied consists of NROW rows and NCOL columns.
!-----------------------------------------------------------------------
      INTEGER ic
!
      DO ic = 1 , Ncol
         CALL DCOPY(Nrow,A(1,ic),1,B(1,ic),1)
      ENDDO
!
!----------------------- End of Subroutine DACOPY ----------------------
      END Subroutine DACOPY
!DECK DVSOL
      SUBROUTINE DVSOL(Wm,Iwm,X,Iersl)
      IMPLICIT NONE
      DOUBLE PRECISION Wm , X
      INTEGER Iwm , Iersl
      DIMENSION Wm(*) , Iwm(*) , X(*)
!-----------------------------------------------------------------------
! Call sequence input -- WM, IWM, X
! Call sequence output -- X, IERSL
! COMMON block variables accessed:
!     /DVOD01/ -- H, RL1, MITER, N
!
! Subroutines called by DVSOL: DGESL, DGBSL
! Function routines called by DVSOL: None
!-----------------------------------------------------------------------
! This routine manages the solution of the linear system arising from
! a chord iteration.  It is called if MITER .ne. 0.
! If MITER is 1 or 2, it calls DGESL to accomplish this.
! If MITER = 3 it updates the coefficient H*RL1 in the diagonal
! matrix, and then computes the solution.
! If MITER is 4 or 5, it calls DGBSL.
! Communication with DVSOL uses the following variables:
! WM    = Real work space containing the inverse diagonal matrix if
!         MITER = 3 and the LU decomposition of the matrix otherwise.
!         Storage of matrix elements starts at WM(3).
!         WM also contains the following matrix-related data:
!         WM(1) = SQRT(UROUND) (not used here),
!         WM(2) = HRL1, the previous value of H*RL1, used if MITER = 3.
! IWM   = Integer work space containing pivot information, starting at
!         IWM(31), if MITER is 1, 2, 4, or 5.  IWM also contains band
!         parameters ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
! X     = The right-hand side vector on input, and the solution vector
!         on output, of length N.
! IERSL = Output flag.  IERSL = 0 if no trouble occurred.
!         IERSL = 1 if a singular matrix arose with MITER = 3.
!-----------------------------------------------------------------------
!
! Type declarations for labeled COMMON block DVOD01 --------------------
!
      DOUBLE PRECISION ACNrm , CCMxj , CONp , CRAte , DRC , EL , ETA ,  &
                       ETAmax , H , HMIn , HMXi , HNEw , HSCal , PRL1 , &
                       RC , RL1 , TAU , TQ , TN , UROund
      INTEGER ICF , INIt , IPUp , JCUr , JSTart , JSV , KFLag , KUTh ,  &
              L , LMAx , LYH , LEWt , LACor , LSAvf , LWM , LIWm ,      &
              LOCjs , MAXord , METh , MITer , MSBj , MXHnil , MXStep ,  &
              N , NEWh , NEWq , NHNil , NQ , NQNyh , NQWait , NSLj ,    &
              NSLp , NYH
!
! Type declarations for local variables --------------------------------
!
      INTEGER i , meband , ml , mu
      DOUBLE PRECISION di , hrl1 , one , phrl1 , r , zero
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
      SAVE one , zero
!
      COMMON /DVOD01/ ACNrm , CCMxj , CONp , CRAte , DRC , EL(13) ,     &
                      ETA , ETAmax , H , HMIn , HMXi , HNEw , HSCal ,   &
                      PRL1 , RC , RL1 , TAU(13) , TQ(5) , TN , UROund , &
                      ICF , INIt , IPUp , JCUr , JSTart , JSV , KFLag , &
                      KUTh , L , LMAx , LYH , LEWt , LACor , LSAvf ,    &
                      LWM , LIWm , LOCjs , MAXord , METh , MITer ,      &
                      MSBj , MXHnil , MXStep , N , NEWh , NEWq , NHNil ,&
                      NQ , NQNyh , NQWait , NSLj , NSLp , NYH
!
      DATA one/1.0D0/ , zero/0.0D0/
!
      Iersl = 0
      SELECT CASE (MITer)
      CASE (3)
!
         phrl1 = Wm(2)
         hrl1 = H*RL1
         Wm(2) = hrl1
         IF ( hrl1/=phrl1 ) THEN
            r = hrl1/phrl1
            DO i = 1 , N
               di = one - r*(one-one/Wm(i+2))
               IF ( ABS(di)==zero ) GOTO 100
               Wm(i+2) = one/di
            ENDDO
         ENDIF
!
         DO i = 1 , N
            X(i) = Wm(i+2)*X(i)
         ENDDO
         RETURN
      CASE (4,5)
!
         ml = Iwm(1)
         mu = Iwm(2)
         meband = 2*ml + mu + 1
         CALL DGBSL(Wm(3),meband,N,ml,mu,Iwm(31),X,0)
         GOTO 99999
      CASE DEFAULT
         CALL DGESL(Wm(3),N,N,Iwm(31),X,0)
         RETURN
      END SELECT
 100  Iersl = 1
      RETURN
!----------------------- End of Subroutine DVSOL -----------------------
99999 END Subroutine DVSOL
!DECK DVSRCO
      SUBROUTINE DVSRCO(Rsav,Isav,Job)
      IMPLICIT NONE
      DOUBLE PRECISION Rsav
      INTEGER Isav , Job
      DIMENSION Rsav(*) , Isav(*)
!-----------------------------------------------------------------------
! Call sequence input -- RSAV, ISAV, JOB
! Call sequence output -- RSAV, ISAV
! COMMON block variables accessed -- All of /DVOD01/ and /DVOD02/
!
! Subroutines/functions called by DVSRCO: None
!-----------------------------------------------------------------------
! This routine saves or restores (depending on JOB) the contents of the
! COMMON blocks DVOD01 and DVOD02, which are used internally by DVODE.
!
! RSAV = real array of length 49 or more.
! ISAV = integer array of length 41 or more.
! JOB  = flag indicating to save or restore the COMMON blocks:
!        JOB  = 1 if COMMON is to be saved (written to RSAV/ISAV).
!        JOB  = 2 if COMMON is to be restored (read from RSAV/ISAV).
!        A call with JOB = 2 presumes a prior call with JOB = 1.
!-----------------------------------------------------------------------
      DOUBLE PRECISION RVOd1 , RVOd2
      INTEGER IVOd1 , IVOd2
      INTEGER i , leniv1 , leniv2 , lenrv1 , lenrv2
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
      SAVE lenrv1 , leniv1 , lenrv2 , leniv2
!
      COMMON /DVOD01/ RVOd1(48) , IVOd1(33)
      COMMON /DVOD02/ RVOd2(1) , IVOd2(8)
      DATA lenrv1/48/ , leniv1/33/ , lenrv2/1/ , leniv2/8/
!
      IF ( Job==2 ) THEN
!
         DO i = 1 , lenrv1
            RVOd1(i) = Rsav(i)
         ENDDO
         DO i = 1 , lenrv2
            RVOd2(i) = Rsav(lenrv1+i)
         ENDDO
!
         DO i = 1 , leniv1
            IVOd1(i) = Isav(i)
         ENDDO
         DO i = 1 , leniv2
            IVOd2(i) = Isav(leniv1+i)
         ENDDO
         GOTO 99999
      ENDIF
      DO i = 1 , lenrv1
         Rsav(i) = RVOd1(i)
      ENDDO
      DO i = 1 , lenrv2
         Rsav(lenrv1+i) = RVOd2(i)
      ENDDO
!
      DO i = 1 , leniv1
         Isav(i) = IVOd1(i)
      ENDDO
      DO i = 1 , leniv2
         Isav(leniv1+i) = IVOd2(i)
      ENDDO
!
      RETURN
!
!----------------------- End of Subroutine DVSRCO ----------------------
99999 END Subroutine DVSRCO
!DECK DEWSET
      SUBROUTINE DEWSET(N,Itol,Rtol,Atol,Ycur,Ewt)
      IMPLICIT NONE
!***BEGIN PROLOGUE  DEWSET
!***SUBSIDIARY
!***PURPOSE  Set error weight vector.
!***TYPE      DOUBLE PRECISION (SEWSET-S, DEWSET-D)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!  This subroutine sets the error weight vector EWT according to
!      EWT(i) = RTOL(i)*ABS(YCUR(i)) + ATOL(i),  i = 1,...,N,
!  with the subscript on RTOL and/or ATOL possibly replaced by 1 above,
!  depending on the value of ITOL.
!
!***SEE ALSO  DLSODE
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791129  DATE WRITTEN
!   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
!   890503  Minor cosmetic changes.  (FNF)
!   930809  Renamed to allow single/double precision versions. (ACH)
!***END PROLOGUE  DEWSET
!**End
      INTEGER N , Itol
      INTEGER i
      DOUBLE PRECISION Rtol , Atol , Ycur , Ewt
      DIMENSION Rtol(*) , Atol(*) , Ycur(N) , Ewt(N)
!
!***FIRST EXECUTABLE STATEMENT  DEWSET
      SELECT CASE (Itol)
      CASE (2)
         DO i = 1 , N
            Ewt(i) = Rtol(1)*ABS(Ycur(i)) + Atol(i)
         ENDDO
         RETURN
      CASE (3)
         DO i = 1 , N
            Ewt(i) = Rtol(i)*ABS(Ycur(i)) + Atol(1)
         ENDDO
         RETURN
      CASE (4)
         DO i = 1 , N
            Ewt(i) = Rtol(i)*ABS(Ycur(i)) + Atol(i)
         ENDDO
         GOTO 99999
      CASE DEFAULT
      END SELECT
      DO i = 1 , N
         Ewt(i) = Rtol(1)*ABS(Ycur(i)) + Atol(1)
      ENDDO
      RETURN
!----------------------- END OF SUBROUTINE DEWSET ----------------------
99999 END SUBROUTINE DEWSET
!DECK DVNORM
      DOUBLE PRECISION FUNCTION DVNORM(N,V,W)
      IMPLICIT NONE
!***BEGIN PROLOGUE  DVNORM
!***SUBSIDIARY
!***PURPOSE  Weighted root-mean-square vector norm.
!***TYPE      DOUBLE PRECISION (SVNORM-S, DVNORM-D)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!  This function routine computes the weighted root-mean-square norm
!  of the vector of length N contained in the array V, with weights
!  contained in the array W of length N:
!    DVNORM = SQRT( (1/N) * SUM( V(i)*W(i) )**2 )
!
!***SEE ALSO  DLSODE
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791129  DATE WRITTEN
!   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
!   890503  Minor cosmetic changes.  (FNF)
!   930809  Renamed to allow single/double precision versions. (ACH)
!***END PROLOGUE  DVNORM
!**End
      INTEGER N , i
      DOUBLE PRECISION V , W , sum
      DIMENSION V(N) , W(N)
!
!***FIRST EXECUTABLE STATEMENT  DVNORM
      sum = 0.0D0
      DO i = 1 , N
         sum = sum + (V(i)*W(i))**2
      ENDDO
      DVNORM = SQRT(sum/N)
!----------------------- END OF FUNCTION DVNORM ------------------------
      END FUNCTION DVNORM
!DECK XERRWD
      SUBROUTINE XERRWD(Msg,Nmes,Nerr,Level,Ni,I1,I2,Nr,R1,R2)
      IMPLICIT NONE
!***BEGIN PROLOGUE  XERRWD
!***SUBSIDIARY
!***PURPOSE  Write error message with values.
!***CATEGORY  R3C
!***TYPE      DOUBLE PRECISION (XERRWV-S, XERRWD-D)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!  Subroutines XERRWD, XSETF, XSETUN, and the function routine IXSAV,
!  as given here, constitute a simplified version of the SLATEC error
!  handling package.
!
!  All arguments are input arguments.
!
!  MSG    = The message (character array).
!  NMES   = The length of MSG (number of characters).
!  NERR   = The error number (not used).
!  LEVEL  = The error level..
!           0 or 1 means recoverable (control returns to caller).
!           2 means fatal (run is aborted--see note below).
!  NI     = Number of integers (0, 1, or 2) to be printed with message.
!  I1,I2  = Integers to be printed, depending on NI.
!  NR     = Number of reals (0, 1, or 2) to be printed with message.
!  R1,R2  = Reals to be printed, depending on NR.
!
!  Note..  this routine is machine-dependent and specialized for use
!  in limited context, in the following ways..
!  1. The argument MSG is assumed to be of type CHARACTER, and
!     the message is printed with a format of (1X,A).
!  2. The message is assumed to take only one line.
!     Multi-line messages are generated by repeated calls.
!  3. If LEVEL = 2, control passes to the statement   STOP
!     to abort the run.  This statement may be machine-dependent.
!  4. R1 and R2 are assumed to be in double precision and are printed
!     in D21.13 format.
!
!***ROUTINES CALLED  IXSAV
!***REVISION HISTORY  (YYMMDD)
!   920831  DATE WRITTEN
!   921118  Replaced MFLGSV/LUNSAV by IXSAV. (ACH)
!   930329  Modified prologue to SLATEC format. (FNF)
!   930407  Changed MSG from CHARACTER*1 array to variable. (FNF)
!   930922  Minor cosmetic change. (FNF)
!***END PROLOGUE  XERRWD
!
!*Internal Notes:
!
! For a different default logical unit number, IXSAV (or a subsidiary
! routine that it calls) will need to be modified.
! For a different run-abort command, change the statement following
! statement 100 at the end.
!-----------------------------------------------------------------------
! Subroutines called by XERRWD.. None
! Function routine called by XERRWD.. IXSAV
!-----------------------------------------------------------------------
!**End
!
!  Declare arguments.
!
      DOUBLE PRECISION R1 , R2
      INTEGER Nmes , Nerr , Level , Ni , I1 , I2 , Nr
      CHARACTER*(*) Msg
!
!  Declare local variables.
!
      INTEGER lunit , mesflg
!
!  Get logical unit number and message print flag.
!
!***FIRST EXECUTABLE STATEMENT  XERRWD
      lunit = IXSAV(1,0,.FALSE.)
      mesflg = IXSAV(2,0,.FALSE.)
      IF ( mesflg/=0 ) THEN
!
!  Write the message.
!
         WRITE (lunit,99001) Msg
99001    FORMAT (1X,A)
         IF ( Ni==1 ) WRITE (lunit,99002) I1
99002    FORMAT (6X,'In above message,  I1 =',I10)
         IF ( Ni==2 ) WRITE (lunit,99003) I1 , I2
99003    FORMAT (6X,'In above message,  I1 =',I10,3X,'I2 =',I10)
         IF ( Nr==1 ) WRITE (lunit,99004) R1
99004    FORMAT (6X,'In above message,  R1 =',D21.13)
         IF ( Nr==2 ) WRITE (lunit,99005) R1 , R2
99005    FORMAT (6X,'In above,  R1 =',D21.13,3X,'R2 =',D21.13)
      ENDIF
!
!  Abort the run if LEVEL = 2.
!
      IF ( Level/=2 ) RETURN
      STOP
!----------------------- End of Subroutine XERRWD ----------------------
      END Subroutine XERRWD
!DECK XSETF
      SUBROUTINE XSETF(Mflag)
      IMPLICIT NONE
!***BEGIN PROLOGUE  XSETF
!***PURPOSE  Reset the error print control flag.
!***CATEGORY  R3A
!***TYPE      ALL (XSETF-A)
!***KEYWORDS  ERROR CONTROL
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!   XSETF sets the error print control flag to MFLAG:
!      MFLAG=1 means print all messages (the default).
!      MFLAG=0 means no printing.
!
!***SEE ALSO  XERRWD, XERRWV
!***REFERENCES  (NONE)
!***ROUTINES CALLED  IXSAV
!***REVISION HISTORY  (YYMMDD)
!   921118  DATE WRITTEN
!   930329  Added SLATEC format prologue. (FNF)
!   930407  Corrected SEE ALSO section. (FNF)
!   930922  Made user-callable, and other cosmetic changes. (FNF)
!***END PROLOGUE  XSETF
!
! Subroutines called by XSETF.. None
! Function routine called by XSETF.. IXSAV
!-----------------------------------------------------------------------
!**End
      INTEGER Mflag , junk
!
!***FIRST EXECUTABLE STATEMENT  XSETF
      IF ( Mflag==0 .OR. Mflag==1 ) junk = IXSAV(2,Mflag,.TRUE.)
!----------------------- End of Subroutine XSETF -----------------------
      END Subroutine XSETF
!DECK XSETUN
      SUBROUTINE XSETUN(Lun)
      IMPLICIT NONE
!***BEGIN PROLOGUE  XSETUN
!***PURPOSE  Reset the logical unit number for error messages.
!***CATEGORY  R3B
!***TYPE      ALL (XSETUN-A)
!***KEYWORDS  ERROR CONTROL
!***DESCRIPTION
!
!   XSETUN sets the logical unit number for error messages to LUN.
!
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***SEE ALSO  XERRWD, XERRWV
!***REFERENCES  (NONE)
!***ROUTINES CALLED  IXSAV
!***REVISION HISTORY  (YYMMDD)
!   921118  DATE WRITTEN
!   930329  Added SLATEC format prologue. (FNF)
!   930407  Corrected SEE ALSO section. (FNF)
!   930922  Made user-callable, and other cosmetic changes. (FNF)
!***END PROLOGUE  XSETUN
!
! Subroutines called by XSETUN.. None
! Function routine called by XSETUN.. IXSAV
!-----------------------------------------------------------------------
!**End
      INTEGER Lun , junk
!
!***FIRST EXECUTABLE STATEMENT  XSETUN
      IF ( Lun>0 ) junk = IXSAV(1,Lun,.TRUE.)
!----------------------- End of Subroutine XSETUN ----------------------
      END Subroutine XSETUN
!DECK IXSAV
      INTEGER FUNCTION IXSAV(Ipar,Ivalue,Iset)
      IMPLICIT NONE
!***BEGIN PROLOGUE  IXSAV
!***SUBSIDIARY
!***PURPOSE  Save and recall error message control parameters.
!***CATEGORY  R3C
!***TYPE      ALL (IXSAV-A)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!  IXSAV saves and recalls one of two error message parameters:
!    LUNIT, the logical unit number to which messages are printed, and
!    MESFLG, the message print flag.
!  This is a modification of the SLATEC library routine J4SAVE.
!
!  Saved local variables..
!   LUNIT  = Logical unit number for messages.  The default is obtained
!            by a call to IUMACH (may be machine-dependent).
!   MESFLG = Print control flag..
!            1 means print all messages (the default).
!            0 means no printing.
!
!  On input..
!    IPAR   = Parameter indicator (1 for LUNIT, 2 for MESFLG).
!    IVALUE = The value to be set for the parameter, if ISET = .TRUE.
!    ISET   = Logical flag to indicate whether to read or write.
!             If ISET = .TRUE., the parameter will be given
!             the value IVALUE.  If ISET = .FALSE., the parameter
!             will be unchanged, and IVALUE is a dummy argument.
!
!  On return..
!    IXSAV = The (old) value of the parameter.
!
!***SEE ALSO  XERRWD, XERRWV
!***ROUTINES CALLED  IUMACH
!***REVISION HISTORY  (YYMMDD)
!   921118  DATE WRITTEN
!   930329  Modified prologue to SLATEC format. (FNF)
!   930915  Added IUMACH call to get default output unit.  (ACH)
!   930922  Minor cosmetic changes. (FNF)
!   010425  Type declaration for IUMACH added. (ACH)
!***END PROLOGUE  IXSAV
!
! Subroutines called by IXSAV.. None
! Function routine called by IXSAV.. IUMACH
!-----------------------------------------------------------------------
!**End
      LOGICAL Iset
      INTEGER Ipar , Ivalue
!-----------------------------------------------------------------------
      INTEGER lunit , mesflg
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this routine.
!-----------------------------------------------------------------------
      SAVE lunit , mesflg
      DATA lunit/ - 1/ , mesflg/1/
!
!***FIRST EXECUTABLE STATEMENT  IXSAV
      IF ( Ipar==1 ) THEN
         IF ( lunit==-1 ) lunit = IUMACH()
         IXSAV = lunit
         IF ( Iset ) lunit = Ivalue
      ENDIF
!
      IF ( Ipar==2 ) THEN
         IXSAV = mesflg
         IF ( Iset ) mesflg = Ivalue
      ENDIF
!
!----------------------- End of Function IXSAV -------------------------
      END Function IXSAV
!DECK IUMACH
      INTEGER FUNCTION IUMACH()
      IMPLICIT NONE
!***BEGIN PROLOGUE  IUMACH
!***PURPOSE  Provide standard output unit number.
!***CATEGORY  R1
!***TYPE      INTEGER (IUMACH-I)
!***KEYWORDS  MACHINE CONSTANTS
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
! *Usage:
!        INTEGER  LOUT, IUMACH
!        LOUT = IUMACH()
!
! *Function Return Values:
!     LOUT : the standard logical unit for Fortran output.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   930915  DATE WRITTEN
!   930922  Made user-callable, and other cosmetic changes. (FNF)
!***END PROLOGUE  IUMACH
!
!*Internal Notes:
!  The built-in value of 6 is standard on a wide range of Fortran
!  systems.  This may be machine-dependent.
!**End
!***FIRST EXECUTABLE STATEMENT  IUMACH
      IUMACH = 6
!
!----------------------- End of Function IUMACH ------------------------
      END Function IUMACH
!DECK DUMACH
      DOUBLE PRECISION FUNCTION DUMACH()
      IMPLICIT NONE
!***BEGIN PROLOGUE  DUMACH
!***PURPOSE  Compute the unit roundoff of the machine.
!***CATEGORY  R1
!***TYPE      DOUBLE PRECISION (RUMACH-S, DUMACH-D)
!***KEYWORDS  MACHINE CONSTANTS
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
! *Usage:
!        DOUBLE PRECISION  A, DUMACH
!        A = DUMACH()
!
! *Function Return Values:
!     A : the unit roundoff of the machine.
!
! *Description:
!     The unit roundoff is defined as the smallest positive machine
!     number u such that  1.0 + u .ne. 1.0.  This is computed by DUMACH
!     in a machine-independent manner.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  DUMSUM
!***REVISION HISTORY  (YYYYMMDD)
!   19930216  DATE WRITTEN
!   19930818  Added SLATEC-format prologue.  (FNF)
!   20030707  Added DUMSUM to force normal storage of COMP.  (ACH)
!***END PROLOGUE  DUMACH
!
      DOUBLE PRECISION u , comp
!***FIRST EXECUTABLE STATEMENT  DUMACH
      u = 1.0D0
 100  u = u*0.5D0
      CALL DUMSUM(1.0D0,u,comp)
      IF ( comp/=1.0D0 ) GOTO 100
      DUMACH = u*2.0D0
!----------------------- End of Function DUMACH ------------------------
      END Function DUMACH

      SUBROUTINE DUMSUM(A,B,C)
      IMPLICIT NONE
!     Routine to force normal storing of A + B, for DUMACH.
      DOUBLE PRECISION A , B , C
      C = A + B
      END SUBROUTINE DUMSUM

!DECK DGEFA
      SUBROUTINE DGEFA(A,Lda,N,Ipvt,Info)
      IMPLICIT NONE
!***BEGIN PROLOGUE  DGEFA
!***PURPOSE  Factor a matrix using Gaussian elimination.
!***CATEGORY  D2A1
!***TYPE      DOUBLE PRECISION (SGEFA-S, DGEFA-D, CGEFA-C)
!***KEYWORDS  GENERAL MATRIX, LINEAR ALGEBRA, LINPACK,
!             MATRIX FACTORIZATION
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     DGEFA factors a double precision matrix by Gaussian elimination.
!
!     DGEFA is usually called by DGECO, but it can be called
!     directly with a saving in time if  RCOND  is not needed.
!     (Time for DGECO) = (1 + 9/N)*(Time for DGEFA) .
!
!     On Entry
!
!        A       DOUBLE PRECISION(LDA, N)
!                the matrix to be factored.
!
!        LDA     INTEGER
!                the leading dimension of the array  A .
!
!        N       INTEGER
!                the order of the matrix  A .
!
!     On Return
!
!        A       an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                The factorization can be written  A = L*U  where
!                L  is a product of permutation and unit lower
!                triangular matrices and  U  is upper triangular.
!
!        IPVT    INTEGER(N)
!                an integer vector of pivot indices.
!
!        INFO    INTEGER
!                = 0  normal value.
!                = K  if  U(K,K) .EQ. 0.0 .  This is not an error
!                     condition for this subroutine, but it does
!                     indicate that DGESL or DGEDI will divide by zero
!                     if called.  Use  RCOND  in DGECO for a reliable
!                     indication of singularity.
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  DAXPY, DSCAL, IDAMAX
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DGEFA
      INTEGER Lda , N , Ipvt(*) , Info
      DOUBLE PRECISION A(Lda,*)
!
      DOUBLE PRECISION t
      INTEGER j , k , kp1 , l , nm1
!
!     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
!
!***FIRST EXECUTABLE STATEMENT  DGEFA
      Info = 0
      nm1 = N - 1
      IF ( nm1>=1 ) THEN
         DO k = 1 , nm1
            kp1 = k + 1
!
!        FIND L = PIVOT INDEX
!
            l = IDAMAX(N-k+1,A(k,k),1) + k - 1
            Ipvt(k) = l
!
!        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
!
            IF ( A(l,k)==0.0D0 ) THEN
               Info = k
            ELSE
!
!           INTERCHANGE IF NECESSARY
!
               IF ( l/=k ) THEN
                  t = A(l,k)
                  A(l,k) = A(k,k)
                  A(k,k) = t
               ENDIF
!
!           COMPUTE MULTIPLIERS
!
               t = -1.0D0/A(k,k)
               CALL DSCAL(N-k,t,A(k+1,k),1)
!
!           ROW ELIMINATION WITH COLUMN INDEXING
!
               DO j = kp1 , N
                  t = A(l,j)
                  IF ( l/=k ) THEN
                     A(l,j) = A(k,j)
                     A(k,j) = t
                  ENDIF
                  CALL DAXPY(N-k,t,A(k+1,k),1,A(k+1,j),1)
               ENDDO
            ENDIF
         ENDDO
      ENDIF
      Ipvt(N) = N
      IF ( A(N,N)==0.0D0 ) Info = N
      END SUBROUTINE DGEFA

!DECK DGESL
      SUBROUTINE DGESL(A,Lda,N,Ipvt,B,Job)
      IMPLICIT NONE
!***BEGIN PROLOGUE  DGESL
!***PURPOSE  Solve the real system A*X=B or TRANS(A)*X=B using the
!            factors computed by DGECO or DGEFA.
!***CATEGORY  D2A1
!***TYPE      DOUBLE PRECISION (SGESL-S, DGESL-D, CGESL-C)
!***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     DGESL solves the double precision system
!     A * X = B  or  TRANS(A) * X = B
!     using the factors computed by DGECO or DGEFA.
!
!     On Entry
!
!        A       DOUBLE PRECISION(LDA, N)
!                the output from DGECO or DGEFA.
!
!        LDA     INTEGER
!                the leading dimension of the array  A .
!
!        N       INTEGER
!                the order of the matrix  A .
!
!        IPVT    INTEGER(N)
!                the pivot vector from DGECO or DGEFA.
!
!        B       DOUBLE PRECISION(N)
!                the right hand side vector.
!
!        JOB     INTEGER
!                = 0         to solve  A*X = B ,
!                = nonzero   to solve  TRANS(A)*X = B  where
!                            TRANS(A)  is the transpose.
!
!     On Return
!
!        B       the solution vector  X .
!
!     Error Condition
!
!        A division by zero will occur if the input factor contains a
!        zero on the diagonal.  Technically this indicates singularity
!        but it is often caused by improper arguments or improper
!        setting of LDA .  It will not occur if the subroutines are
!        called correctly and if DGECO has set RCOND .GT. 0.0
!        or DGEFA has set INFO .EQ. 0 .
!
!     To compute  INVERSE(A) * C  where  C  is a matrix
!     with  P  columns
!           CALL DGECO(A,LDA,N,IPVT,RCOND,Z)
!           IF (RCOND is too small) GO TO ...
!           DO 10 J = 1, P
!              CALL DGESL(A,LDA,N,IPVT,C(1,J),0)
!        10 CONTINUE
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  DAXPY, DDOT
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DGESL
      INTEGER Lda , N , Ipvt(*) , Job
      DOUBLE PRECISION A(Lda,*) , B(*)
!
      DOUBLE PRECISION t
      INTEGER k , kb , l , nm1
!***FIRST EXECUTABLE STATEMENT  DGESL
      nm1 = N - 1
      IF ( Job/=0 ) THEN
!
!        JOB = NONZERO, SOLVE  TRANS(A) * X = B
!        FIRST SOLVE  TRANS(U)*Y = B
!
         DO k = 1 , N
            t = DDOT(k-1,A(1,k),1,B(1),1)
            B(k) = (B(k)-t)/A(k,k)
         ENDDO
!
!        NOW SOLVE TRANS(L)*X = Y
!
         IF ( nm1>=1 ) THEN
            DO kb = 1 , nm1
               k = N - kb
               B(k) = B(k) + DDOT(N-k,A(k+1,k),1,B(k+1),1)
               l = Ipvt(k)
               IF ( l/=k ) THEN
                  t = B(l)
                  B(l) = B(k)
                  B(k) = t
               ENDIF
            ENDDO
         ENDIF
      ELSE
!
!        JOB = 0 , SOLVE  A * X = B
!        FIRST SOLVE  L*Y = B
!
         IF ( nm1>=1 ) THEN
            DO k = 1 , nm1
               l = Ipvt(k)
               t = B(l)
               IF ( l/=k ) THEN
                  B(l) = B(k)
                  B(k) = t
               ENDIF
               CALL DAXPY(N-k,t,A(k+1,k),1,B(k+1),1)
            ENDDO
         ENDIF
!
!        NOW SOLVE  U*X = Y
!
         DO kb = 1 , N
            k = N + 1 - kb
            B(k) = B(k)/A(k,k)
            t = -B(k)
            CALL DAXPY(k-1,t,A(1,k),1,B(1),1)
         ENDDO
      ENDIF
      END SUBROUTINE DGESL
!DECK DGBFA
      SUBROUTINE DGBFA(Abd,Lda,N,Ml,Mu,Ipvt,Info)
      IMPLICIT NONE
!***BEGIN PROLOGUE  DGBFA
!***PURPOSE  Factor a band matrix using Gaussian elimination.
!***CATEGORY  D2A2
!***TYPE      DOUBLE PRECISION (SGBFA-S, DGBFA-D, CGBFA-C)
!***KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     DGBFA factors a double precision band matrix by elimination.
!
!     DGBFA is usually called by DGBCO, but it can be called
!     directly with a saving in time if  RCOND  is not needed.
!
!     On Entry
!
!        ABD     DOUBLE PRECISION(LDA, N)
!                contains the matrix in band storage.  The columns
!                of the matrix are stored in the columns of  ABD  and
!                the diagonals of the matrix are stored in rows
!                ML+1 through 2*ML+MU+1 of  ABD .
!                See the comments below for details.
!
!        LDA     INTEGER
!                the leading dimension of the array  ABD .
!                LDA must be .GE. 2*ML + MU + 1 .
!
!        N       INTEGER
!                the order of the original matrix.
!
!        ML      INTEGER
!                number of diagonals below the main diagonal.
!                0 .LE. ML .LT.  N .
!
!        MU      INTEGER
!                number of diagonals above the main diagonal.
!                0 .LE. MU .LT.  N .
!                More efficient if  ML .LE. MU .
!     On Return
!
!        ABD     an upper triangular matrix in band storage and
!                the multipliers which were used to obtain it.
!                The factorization can be written  A = L*U  where
!                L  is a product of permutation and unit lower
!                triangular matrices and  U  is upper triangular.
!
!        IPVT    INTEGER(N)
!                an integer vector of pivot indices.
!
!        INFO    INTEGER
!                = 0  normal value.
!                = K  if  U(K,K) .EQ. 0.0 .  This is not an error
!                     condition for this subroutine, but it does
!                     indicate that DGBSL will divide by zero if
!                     called.  Use  RCOND  in DGBCO for a reliable
!                     indication of singularity.
!
!     Band Storage
!
!           If  A  is a band matrix, the following program segment
!           will set up the input.
!
!                   ML = (band width below the diagonal)
!                   MU = (band width above the diagonal)
!                   M = ML + MU + 1
!                   DO 20 J = 1, N
!                      I1 = MAX(1, J-MU)
!                      I2 = MIN(N, J+ML)
!                      DO 10 I = I1, I2
!                         K = I - J + M
!                         ABD(K,J) = A(I,J)
!                10    CONTINUE
!                20 CONTINUE
!
!           This uses rows  ML+1  through  2*ML+MU+1  of  ABD .
!           In addition, the first  ML  rows in  ABD  are used for
!           elements generated during the triangularization.
!           The total number of rows needed in  ABD  is  2*ML+MU+1 .
!           The  ML+MU by ML+MU  upper left triangle and the
!           ML by ML  lower right triangle are not referenced.
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  DAXPY, DSCAL, IDAMAX
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DGBFA
      INTEGER Lda , N , Ml , Mu , Ipvt(*) , Info
      DOUBLE PRECISION Abd(Lda,*)
!
      DOUBLE PRECISION t
      INTEGER i , i0 , j , ju , jz , j0 , j1 , k , kp1 , l ,   &
              lm , m , mm , nm1
!
!***FIRST EXECUTABLE STATEMENT  DGBFA
      m = Ml + Mu + 1
      Info = 0
!
!     ZERO INITIAL FILL-IN COLUMNS
!
      j0 = Mu + 2
      j1 = MIN(N,m) - 1
      IF ( j1>=j0 ) THEN
         DO jz = j0 , j1
            i0 = m + 1 - jz
            DO i = i0 , Ml
               Abd(i,jz) = 0.0D0
            ENDDO
         ENDDO
      ENDIF
      jz = j1
      ju = 0
!
!     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
!
      nm1 = N - 1
      IF ( nm1>=1 ) THEN
         DO k = 1 , nm1
            kp1 = k + 1
!
!        ZERO NEXT FILL-IN COLUMN
!
            jz = jz + 1
            IF ( jz<=N ) THEN
               IF ( Ml>=1 ) THEN
                  DO i = 1 , Ml
                     Abd(i,jz) = 0.0D0
                  ENDDO
               ENDIF
            ENDIF
!
!        FIND L = PIVOT INDEX
!
            lm = MIN(Ml,N-k)
            l = IDAMAX(lm+1,Abd(m,k),1) + m - 1
            Ipvt(k) = l + k - m
!
!        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
!
            IF ( Abd(l,k)==0.0D0 ) THEN
               Info = k
            ELSE
!
!           INTERCHANGE IF NECESSARY
!
               IF ( l/=m ) THEN
                  t = Abd(l,k)
                  Abd(l,k) = Abd(m,k)
                  Abd(m,k) = t
               ENDIF
!
!           COMPUTE MULTIPLIERS
!
               t = -1.0D0/Abd(m,k)
               CALL DSCAL(lm,t,Abd(m+1,k),1)
!
!           ROW ELIMINATION WITH COLUMN INDEXING
!
               ju = MIN(MAX(ju,Mu+Ipvt(k)),N)
               mm = m
               IF ( ju>=kp1 ) THEN
                  DO j = kp1 , ju
                     l = l - 1
                     mm = mm - 1
                     t = Abd(l,j)
                     IF ( l/=mm ) THEN
                        Abd(l,j) = Abd(mm,j)
                        Abd(mm,j) = t
                     ENDIF
                     CALL DAXPY(lm,t,Abd(m+1,k),1,Abd(mm+1,j),1)
                  ENDDO
               ENDIF
            ENDIF
         ENDDO
      ENDIF
      Ipvt(N) = N
      IF ( Abd(m,N)==0.0D0 ) Info = N
      END SUBROUTINE DGBFA
!DECK DGBSL
      SUBROUTINE DGBSL(Abd,Lda,N,Ml,Mu,Ipvt,B,Job)
      IMPLICIT NONE
!***BEGIN PROLOGUE  DGBSL
!***PURPOSE  Solve the real band system A*X=B or TRANS(A)*X=B using
!            the factors computed by DGBCO or DGBFA.
!***CATEGORY  D2A2
!***TYPE      DOUBLE PRECISION (SGBSL-S, DGBSL-D, CGBSL-C)
!***KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     DGBSL solves the double precision band system
!     A * X = B  or  TRANS(A) * X = B
!     using the factors computed by DGBCO or DGBFA.
!
!     On Entry
!
!        ABD     DOUBLE PRECISION(LDA, N)
!                the output from DGBCO or DGBFA.
!
!        LDA     INTEGER
!                the leading dimension of the array  ABD .
!
!        N       INTEGER
!                the order of the original matrix.
!
!        ML      INTEGER
!                number of diagonals below the main diagonal.
!
!        MU      INTEGER
!                number of diagonals above the main diagonal.
!
!        IPVT    INTEGER(N)
!                the pivot vector from DGBCO or DGBFA.
!
!        B       DOUBLE PRECISION(N)
!                the right hand side vector.
!
!        JOB     INTEGER
!                = 0         to solve  A*X = B ,
!                = nonzero   to solve  TRANS(A)*X = B , where
!                            TRANS(A)  is the transpose.
!
!     On Return
!
!        B       the solution vector  X .
!
!     Error Condition
!
!        A division by zero will occur if the input factor contains a
!        zero on the diagonal.  Technically this indicates singularity
!        but it is often caused by improper arguments or improper
!        setting of LDA .  It will not occur if the subroutines are
!        called correctly and if DGBCO has set RCOND .GT. 0.0
!        or DGBFA has set INFO .EQ. 0 .
!
!     To compute  INVERSE(A) * C  where  C  is a matrix
!     with  P  columns
!           CALL DGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)
!           IF (RCOND is too small) GO TO ...
!           DO 10 J = 1, P
!              CALL DGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0)
!        10 CONTINUE
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  DAXPY, DDOT
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DGBSL
      INTEGER Lda , N , Ml , Mu , Ipvt(*) , Job
      DOUBLE PRECISION Abd(Lda,*) , B(*)
!
      DOUBLE PRECISION t
      INTEGER k , kb , l , la , lb , lm , m , nm1
!***FIRST EXECUTABLE STATEMENT  DGBSL
      m = Mu + Ml + 1
      nm1 = N - 1
      IF ( Job/=0 ) THEN
!
!        JOB = NONZERO, SOLVE  TRANS(A) * X = B
!        FIRST SOLVE  TRANS(U)*Y = B
!
         DO k = 1 , N
            lm = MIN(k,m) - 1
            la = m - lm
            lb = k - lm
            t = DDOT(lm,Abd(la,k),1,B(lb),1)
            B(k) = (B(k)-t)/Abd(m,k)
         ENDDO
!
!        NOW SOLVE TRANS(L)*X = Y
!
         IF ( Ml/=0 ) THEN
            IF ( nm1>=1 ) THEN
               DO kb = 1 , nm1
                  k = N - kb
                  lm = MIN(Ml,N-k)
                  B(k) = B(k) + DDOT(lm,Abd(m+1,k),1,B(k+1),1)
                  l = Ipvt(k)
                  IF ( l/=k ) THEN
                     t = B(l)
                     B(l) = B(k)
                     B(k) = t
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
      ELSE
!
!        JOB = 0 , SOLVE  A * X = B
!        FIRST SOLVE L*Y = B
!
         IF ( Ml/=0 ) THEN
            IF ( nm1>=1 ) THEN
               DO k = 1 , nm1
                  lm = MIN(Ml,N-k)
                  l = Ipvt(k)
                  t = B(l)
                  IF ( l/=k ) THEN
                     B(l) = B(k)
                     B(k) = t
                  ENDIF
                  CALL DAXPY(lm,t,Abd(m+1,k),1,B(k+1),1)
               ENDDO
            ENDIF
         ENDIF
!
!        NOW SOLVE  U*X = Y
!
         DO kb = 1 , N
            k = N + 1 - kb
            B(k) = B(k)/Abd(m,k)
            lm = MIN(k,m) - 1
            la = m - lm
            lb = k - lm
            t = -B(k)
            CALL DAXPY(lm,t,Abd(la,k),1,B(lb),1)
         ENDDO
      ENDIF
      END SUBROUTINE DGBSL
!DECK DAXPY
      SUBROUTINE DAXPY(N,Da,Dx,Incx,Dy,Incy)
      IMPLICIT NONE
      INTEGER i , Incx , Incy , ix , iy , m , mp1 , N , ns
!***BEGIN PROLOGUE  DAXPY
!***PURPOSE  Compute a constant times a vector plus a vector.
!***CATEGORY  D1A7
!***TYPE      DOUBLE PRECISION (SAXPY-S, DAXPY-D, CAXPY-C)
!***KEYWORDS  BLAS, LINEAR ALGEBRA, TRIAD, VECTOR
!***AUTHOR  Lawson, C. L., (JPL)
!           Hanson, R. J., (SNLA)
!           Kincaid, D. R., (U. of Texas)
!           Krogh, F. T., (JPL)
!***DESCRIPTION
!
!                B L A S  Subprogram
!    Description of Parameters
!
!     --Input--
!        N  number of elements in input vector(s)
!       DA  double precision scalar multiplier
!       DX  double precision vector with N elements
!     INCX  storage spacing between elements of DX
!       DY  double precision vector with N elements
!     INCY  storage spacing between elements of DY
!
!     --Output--
!       DY  double precision result (unchanged if N .LE. 0)
!
!     Overwrite double precision DY with double precision DA*DX + DY.
!     For I = 0 to N-1, replace  DY(LY+I*INCY) with DA*DX(LX+I*INCX) +
!       DY(LY+I*INCY),
!     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
!     defined in a similar way using INCY.
!
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!                 Krogh, Basic linear algebra subprograms for Fortran
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DAXPY
      DOUBLE PRECISION Dx(*) , Dy(*) , Da
!***FIRST EXECUTABLE STATEMENT  DAXPY
      IF ( N<=0 .OR. Da==0.0D0 ) RETURN
      IF ( Incx==Incy ) THEN
         IF ( Incx<1 ) THEN
         ELSEIF ( Incx==1 ) THEN
!
!     Code for both increments equal to 1.
!
!     Clean-up loop so remaining vector length is a multiple of 4.
!
            m = MOD(N,4)
            IF ( m/=0 ) THEN
               DO i = 1 , m
                  Dy(i) = Dy(i) + Da*Dx(i)
               ENDDO
               IF ( N<4 ) RETURN
            ENDIF
            GOTO 100
         ELSE
!
!     Code for equal, positive, non-unit increments.
!
            ns = N*Incx
            DO i = 1 , ns , Incx
               Dy(i) = Da*Dx(i) + Dy(i)
            ENDDO
            GOTO 99999
         ENDIF
      ENDIF
!
!     Code for unequal or nonpositive increments.
!
      ix = 1
      iy = 1
      IF ( Incx<0 ) ix = (-N+1)*Incx + 1
      IF ( Incy<0 ) iy = (-N+1)*Incy + 1
      DO i = 1 , N
         Dy(iy) = Dy(iy) + Da*Dx(ix)
         ix = ix + Incx
         iy = iy + Incy
      ENDDO
      RETURN
 100  mp1 = m + 1
      DO i = mp1 , N , 4
         Dy(i) = Dy(i) + Da*Dx(i)
         Dy(i+1) = Dy(i+1) + Da*Dx(i+1)
         Dy(i+2) = Dy(i+2) + Da*Dx(i+2)
         Dy(i+3) = Dy(i+3) + Da*Dx(i+3)
      ENDDO
      RETURN
99999 END SUBROUTINE DAXPY
!DECK DCOPY
      SUBROUTINE DCOPY(N,Dx,Incx,Dy,Incy)
      IMPLICIT NONE
      INTEGER i , Incx , Incy , ix , iy , m , mp1 , N , ns
!***BEGIN PROLOGUE  DCOPY
!***PURPOSE  Copy a vector.
!***CATEGORY  D1A5
!***TYPE      DOUBLE PRECISION (SCOPY-S, DCOPY-D, CCOPY-C, ICOPY-I)
!***KEYWORDS  BLAS, COPY, LINEAR ALGEBRA, VECTOR
!***AUTHOR  Lawson, C. L., (JPL)
!           Hanson, R. J., (SNLA)
!           Kincaid, D. R., (U. of Texas)
!           Krogh, F. T., (JPL)
!***DESCRIPTION
!
!                B L A S  Subprogram
!    Description of Parameters
!
!     --Input--
!        N  number of elements in input vector(s)
!       DX  double precision vector with N elements
!     INCX  storage spacing between elements of DX
!       DY  double precision vector with N elements
!     INCY  storage spacing between elements of DY
!
!     --Output--
!       DY  copy of vector DX (unchanged if N .LE. 0)
!
!     Copy double precision DX to double precision DY.
!     For I = 0 to N-1, copy DX(LX+I*INCX) to DY(LY+I*INCY),
!     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
!     defined in a similar way using INCY.
!
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!                 Krogh, Basic linear algebra subprograms for Fortran
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DCOPY
      DOUBLE PRECISION Dx(*) , Dy(*)
!***FIRST EXECUTABLE STATEMENT  DCOPY
      IF ( N<=0 ) RETURN
      IF ( Incx==Incy ) THEN
         IF ( Incx<1 ) THEN
         ELSEIF ( Incx==1 ) THEN
!
!     Code for both increments equal to 1.
!
!     Clean-up loop so remaining vector length is a multiple of 7.
!
            m = MOD(N,7)
            IF ( m/=0 ) THEN
               DO i = 1 , m
                  Dy(i) = Dx(i)
               ENDDO
               IF ( N<7 ) RETURN
            ENDIF
            GOTO 100
         ELSE
!
!     Code for equal, positive, non-unit increments.
!
            ns = N*Incx
            DO i = 1 , ns , Incx
               Dy(i) = Dx(i)
            ENDDO
            GOTO 99999
         ENDIF
      ENDIF
!
!     Code for unequal or nonpositive increments.
!
      ix = 1
      iy = 1
      IF ( Incx<0 ) ix = (-N+1)*Incx + 1
      IF ( Incy<0 ) iy = (-N+1)*Incy + 1
      DO i = 1 , N
         Dy(iy) = Dx(ix)
         ix = ix + Incx
         iy = iy + Incy
      ENDDO
      RETURN
 100  mp1 = m + 1
      DO i = mp1 , N , 7
         Dy(i) = Dx(i)
         Dy(i+1) = Dx(i+1)
         Dy(i+2) = Dx(i+2)
         Dy(i+3) = Dx(i+3)
         Dy(i+4) = Dx(i+4)
         Dy(i+5) = Dx(i+5)
         Dy(i+6) = Dx(i+6)
      ENDDO
      RETURN
99999 END SUBROUTINE DCOPY
!DECK DDOT
      DOUBLE PRECISION FUNCTION DDOT(N,Dx,Incx,Dy,Incy)
      IMPLICIT NONE
      INTEGER i , Incx , Incy , ix , iy , m , mp1 , N , ns
!***BEGIN PROLOGUE  DDOT
!***PURPOSE  Compute the inner product of two vectors.
!***CATEGORY  D1A4
!***TYPE      DOUBLE PRECISION (SDOT-S, DDOT-D, CDOTU-C)
!***KEYWORDS  BLAS, INNER PRODUCT, LINEAR ALGEBRA, VECTOR
!***AUTHOR  Lawson, C. L., (JPL)
!           Hanson, R. J., (SNLA)
!           Kincaid, D. R., (U. of Texas)
!           Krogh, F. T., (JPL)
!***DESCRIPTION
!
!                B L A S  Subprogram
!    Description of Parameters
!
!     --Input--
!        N  number of elements in input vector(s)
!       DX  double precision vector with N elements
!     INCX  storage spacing between elements of DX
!       DY  double precision vector with N elements
!     INCY  storage spacing between elements of DY
!
!     --Output--
!     DDOT  double precision dot product (zero if N .LE. 0)
!
!     Returns the dot product of double precision DX and DY.
!     DDOT = sum for I = 0 to N-1 of  DX(LX+I*INCX) * DY(LY+I*INCY),
!     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
!     defined in a similar way using INCY.
!
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!                 Krogh, Basic linear algebra subprograms for Fortran
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DDOT
      DOUBLE PRECISION Dx(*) , Dy(*)
!***FIRST EXECUTABLE STATEMENT  DDOT
      DDOT = 0.0D0
      IF ( N<=0 ) RETURN
      IF ( Incx==Incy ) THEN
         IF ( Incx<1 ) THEN
         ELSEIF ( Incx==1 ) THEN
!
!     Code for both increments equal to 1.
!
!     Clean-up loop so remaining vector length is a multiple of 5.
!
            m = MOD(N,5)
            IF ( m/=0 ) THEN
               DO i = 1 , m
                  DDOT = DDOT + Dx(i)*Dy(i)
               ENDDO
               IF ( N<5 ) RETURN
            ENDIF
            GOTO 100
         ELSE
!
!     Code for equal, positive, non-unit increments.
!
            ns = N*Incx
            DO i = 1 , ns , Incx
               DDOT = DDOT + Dx(i)*Dy(i)
            ENDDO
            GOTO 99999
         ENDIF
      ENDIF
!
!     Code for unequal or nonpositive increments.
!
      ix = 1
      iy = 1
      IF ( Incx<0 ) ix = (-N+1)*Incx + 1
      IF ( Incy<0 ) iy = (-N+1)*Incy + 1
      DO i = 1 , N
         DDOT = DDOT + Dx(ix)*Dy(iy)
         ix = ix + Incx
         iy = iy + Incy
      ENDDO
      RETURN
 100  mp1 = m + 1
      DO i = mp1 , N , 5
         DDOT = DDOT + Dx(i)*Dy(i) + Dx(i+1)*Dy(i+1) + Dx(i+2)*Dy(i+2)  &
                + Dx(i+3)*Dy(i+3) + Dx(i+4)*Dy(i+4)
      ENDDO
      RETURN
99999 END FUNCTION DDOT
!DECK DNRM2
      DOUBLE PRECISION FUNCTION DNRM2(N,Dx,Incx)
      IMPLICIT NONE
      INTEGER i , Incx , j , N , nn
!***BEGIN PROLOGUE  DNRM2
!***PURPOSE  Compute the Euclidean length (L2 norm) of a vector.
!***CATEGORY  D1A3B
!***TYPE      DOUBLE PRECISION (SNRM2-S, DNRM2-D, SCNRM2-C)
!***KEYWORDS  BLAS, EUCLIDEAN LENGTH, EUCLIDEAN NORM, L2,
!             LINEAR ALGEBRA, UNITARY, VECTOR
!***AUTHOR  Lawson, C. L., (JPL)
!           Hanson, R. J., (SNLA)
!           Kincaid, D. R., (U. of Texas)
!           Krogh, F. T., (JPL)
!***DESCRIPTION
!
!                B L A S  Subprogram
!    Description of parameters
!
!     --Input--
!        N  number of elements in input vector(s)
!       DX  double precision vector with N elements
!     INCX  storage spacing between elements of DX
!
!     --Output--
!    DNRM2  double precision result (zero if N .LE. 0)
!
!     Euclidean norm of the N-vector stored in DX with storage
!     increment INCX.
!     If N .LE. 0, return with result = 0.
!     If N .GE. 1, then INCX must be .GE. 1
!
!     Four phase method using two built-in constants that are
!     hopefully applicable to all machines.
!         CUTLO = maximum of  SQRT(U/EPS)  over all known machines.
!         CUTHI = minimum of  SQRT(V)      over all known machines.
!     where
!         EPS = smallest no. such that EPS + 1. .GT. 1.
!         U   = smallest positive no.   (underflow limit)
!         V   = largest  no.            (overflow  limit)
!
!     Brief outline of algorithm.
!
!     Phase 1 scans zero components.
!     move to phase 2 when a component is nonzero and .LE. CUTLO
!     move to phase 3 when a component is .GT. CUTLO
!     move to phase 4 when a component is .GE. CUTHI/M
!     where M = N for X() real and M = 2*N for complex.
!
!     Values for CUTLO and CUTHI.
!     From the environmental parameters listed in the IMSL converter
!     document the limiting values are as follows:
!     CUTLO, S.P.   U/EPS = 2**(-102) for  Honeywell.  Close seconds are
!                   Univac and DEC at 2**(-103)
!                   Thus CUTLO = 2**(-51) = 4.44089E-16
!     CUTHI, S.P.   V = 2**127 for Univac, Honeywell, and DEC.
!                   Thus CUTHI = 2**(63.5) = 1.30438E19
!     CUTLO, D.P.   U/EPS = 2**(-67) for Honeywell and DEC.
!                   Thus CUTLO = 2**(-33.5) = 8.23181D-11
!     CUTHI, D.P.   same as S.P.  CUTHI = 1.30438D19
!     DATA CUTLO, CUTHI /8.232D-11,  1.304D19/
!     DATA CUTLO, CUTHI /4.441E-16,  1.304E19/
!
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!                 Krogh, Basic linear algebra subprograms for Fortran
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DNRM2
      INTEGER next
      DOUBLE PRECISION Dx(*) , cutlo , cuthi , hitest , sum , xmax ,    &
                       zero , one
      SAVE cutlo , cuthi , zero , one
      DATA zero , one/0.0D0 , 1.0D0/
!
      DATA cutlo , cuthi/8.232D-11 , 1.304D19/
!***FIRST EXECUTABLE STATEMENT  DNRM2
      IF ( N>0 ) THEN
!
         ASSIGN 200 TO next
         sum = zero
         nn = N*Incx
!
!                                                 BEGIN MAIN LOOP
!
         i = 1
      ELSE
         DNRM2 = zero
         GOTO 99999
      ENDIF
 100  GOTO next
 200  IF ( ABS(Dx(i))>cutlo ) GOTO 800
      ASSIGN 300 TO next
      xmax = zero
!
!                        PHASE 1.  SUM IS ZERO
!
 300  IF ( Dx(i)==zero ) GOTO 900
      IF ( ABS(Dx(i))>cutlo ) GOTO 800
!
!                                PREPARE FOR PHASE 2.
!
      ASSIGN 600 TO next
      GOTO 500
!
!                                PREPARE FOR PHASE 4.
!
 400  i = j
      ASSIGN 700 TO next
      sum = (sum/Dx(i))/Dx(i)
 500  xmax = ABS(Dx(i))
!
      sum = sum + (Dx(i)/xmax)**2
      GOTO 900
!
!                   PHASE 2.  SUM IS SMALL.
!                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
!
 600  IF ( ABS(Dx(i))>cutlo ) THEN
!
!                  PREPARE FOR PHASE 3.
!
         sum = (sum*xmax)*xmax
         GOTO 800
      ENDIF
!
!                     COMMON CODE FOR PHASES 2 AND 4.
!                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
!
 700  IF ( ABS(Dx(i))<=xmax ) THEN
         sum = sum + (Dx(i)/xmax)**2
      ELSE
         sum = one + sum*(xmax/Dx(i))**2
         xmax = ABS(Dx(i))
      ENDIF
      GOTO 900
!
!     FOR REAL OR D.P. SET HITEST = CUTHI/N
!     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
!
 800  hitest = cuthi/N
!
!                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
!
      DO j = i , nn , Incx
         IF ( ABS(Dx(j))>=hitest ) GOTO 400
         sum = sum + Dx(j)**2
      ENDDO
      DNRM2 = SQRT(sum)
      GOTO 99999
!
 900  i = i + Incx
      IF ( i<=nn ) GOTO 100
!
!              END OF MAIN LOOP.
!
!              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
!
      DNRM2 = xmax*SQRT(sum)
99999 END FUNCTION DNRM2
!DECK DSCAL
      SUBROUTINE DSCAL(N,Da,Dx,Incx)
      IMPLICIT NONE
!***BEGIN PROLOGUE  DSCAL
!***PURPOSE  Multiply a vector by a constant.
!***CATEGORY  D1A6
!***TYPE      DOUBLE PRECISION (SSCAL-S, DSCAL-D, CSCAL-C)
!***KEYWORDS  BLAS, LINEAR ALGEBRA, SCALE, VECTOR
!***AUTHOR  Lawson, C. L., (JPL)
!           Hanson, R. J., (SNLA)
!           Kincaid, D. R., (U. of Texas)
!           Krogh, F. T., (JPL)
!***DESCRIPTION
!
!                B L A S  Subprogram
!    Description of Parameters
!
!     --Input--
!        N  number of elements in input vector(s)
!       DA  double precision scale factor
!       DX  double precision vector with N elements
!     INCX  storage spacing between elements of DX
!
!     --Output--
!       DX  double precision result (unchanged if N.LE.0)
!
!     Replace double precision DX by double precision DA*DX.
!     For I = 0 to N-1, replace DX(IX+I*INCX) with  DA * DX(IX+I*INCX),
!     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
!
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!                 Krogh, Basic linear algebra subprograms for Fortran
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900821  Modified to correct problem with a negative increment.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DSCAL
      DOUBLE PRECISION Da , Dx(*)
      INTEGER i , Incx , ix , m , mp1 , N
!***FIRST EXECUTABLE STATEMENT  DSCAL
      IF ( N<=0 ) RETURN
      IF ( Incx==1 ) THEN
!
!     Code for increment equal to 1.
!
!     Clean-up loop so remaining vector length is a multiple of 5.
!
         m = MOD(N,5)
         IF ( m/=0 ) THEN
            DO i = 1 , m
               Dx(i) = Da*Dx(i)
            ENDDO
            IF ( N<5 ) RETURN
         ENDIF
         mp1 = m + 1
         DO i = mp1 , N , 5
            Dx(i) = Da*Dx(i)
            Dx(i+1) = Da*Dx(i+1)
            Dx(i+2) = Da*Dx(i+2)
            Dx(i+3) = Da*Dx(i+3)
            Dx(i+4) = Da*Dx(i+4)
         ENDDO
      ELSE
!
!     Code for increment not equal to 1.
!
         ix = 1
         IF ( Incx<0 ) ix = (-N+1)*Incx + 1
         DO i = 1 , N
            Dx(ix) = Da*Dx(ix)
            ix = ix + Incx
         ENDDO
         RETURN
      ENDIF
      END SUBROUTINE DSCAL
!DECK IDAMAX
      INTEGER FUNCTION IDAMAX(N,Dx,Incx)
      IMPLICIT NONE
!***BEGIN PROLOGUE  IDAMAX
!***PURPOSE  Find the smallest index of that component of a vector
!            having the maximum magnitude.
!***CATEGORY  D1A2
!***TYPE      DOUBLE PRECISION (ISAMAX-S, IDAMAX-D, ICAMAX-C)
!***KEYWORDS  BLAS, LINEAR ALGEBRA, MAXIMUM COMPONENT, VECTOR
!***AUTHOR  Lawson, C. L., (JPL)
!           Hanson, R. J., (SNLA)
!           Kincaid, D. R., (U. of Texas)
!           Krogh, F. T., (JPL)
!***DESCRIPTION
!
!                B L A S  Subprogram
!    Description of Parameters
!
!     --Input--
!        N  number of elements in input vector(s)
!       DX  double precision vector with N elements
!     INCX  storage spacing between elements of DX
!
!     --Output--
!   IDAMAX  smallest index (zero if N .LE. 0)
!
!     Find smallest index of maximum magnitude of double precision DX.
!     IDAMAX = first I, I = 1 to N, to maximize ABS(DX(IX+(I-1)*INCX)),
!     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
!
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!                 Krogh, Basic linear algebra subprograms for Fortran
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900821  Modified to correct problem with a negative increment.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  IDAMAX
      DOUBLE PRECISION Dx(*) , dmax , xmag
      INTEGER i , Incx , ix , N
!***FIRST EXECUTABLE STATEMENT  IDAMAX
      IDAMAX = 0
      IF ( N<=0 ) RETURN
      IDAMAX = 1
      IF ( N==1 ) RETURN
!
      IF ( Incx==1 ) THEN
!
!     Code for increments equal to 1.
!
         dmax = ABS(Dx(1))
         DO i = 2 , N
            xmag = ABS(Dx(i))
            IF ( xmag>dmax ) THEN
               IDAMAX = i
               dmax = xmag
            ENDIF
         ENDDO
         GOTO 99999
      ENDIF
!
!     Code for increments not equal to 1.
!
      ix = 1
      IF ( Incx<0 ) ix = (-N+1)*Incx + 1
      dmax = ABS(Dx(ix))
      ix = ix + Incx
      DO i = 2 , N
         xmag = ABS(Dx(ix))
         IF ( xmag>dmax ) THEN
            IDAMAX = i
            dmax = xmag
         ENDIF
         ix = ix + Incx
      ENDDO
      RETURN
99999 END FUNCTION IDAMAX

end module dvode_module
