!*******************************************************************************
!>
!  Module for the refactored DVODE.
!
!### See also
!  * The original DVODE.f (1989-2002) can be found
!    [here](https://computing.llnl.gov/projects/odepack/software)
!
!### History
!  * Translated into modern Fortran by Jacob Williams.
!
!@todo Replace GOTO statements

module dvode_module

    use dvode_blas_module
    use dvode_linpack_module
    use dvode_kinds_module, only: wp => dvode_wp
    use iso_fortran_env,    only: output_unit

    implicit none

    private

    real(wp),parameter :: zero = 0.0_wp
    real(wp),parameter :: epmach = epsilon(1.0_wp) !! machine epsilon
    integer,parameter :: iumach  = output_unit     !! standard output unit number
    integer,dimension(2),parameter :: mord = [12,5]

   type,public :: dvode_data_t

      !! Internal variables for DVODE.
      !! variables which are communicated between subroutines in the
      !! dvode package, or which are to be saved between calls to dvode.

      private

      ! type declarations formerly in common block dvod01
      real(wp) :: acnrm   = zero !! weighted r.m.s. norm of accumulated correction vectors.
      real(wp) :: ccmxj   = zero !! threshhold on drc for updating the jacobian. (see drc.)
      real(wp) :: conp    = zero !! the saved value of tq(5).
      real(wp) :: crate   = zero !! estimated corrector convergence rate constant.
      real(wp) :: drc     = zero !! relative change in h*rl1 since last dvjac call.
      real(wp) :: el(13)  = zero !! real array of integration coefficients.  see dvset.
      real(wp) :: eta     = zero !! saved tentative ratio of new to old h.
      real(wp) :: etamax  = zero !! saved maximum value of eta to be allowed.
      real(wp) :: h       = zero !! the step size.
      real(wp) :: hmin    = zero !! the minimum absolute value of the step size h to be used.
      real(wp) :: hmxi    = zero !! inverse of the maximum absolute value of h to be used.
                                 !! hmxi = 0.0 is allowed and corresponds to an infinite hmax.
      real(wp) :: hnew    = zero !! the step size to be attempted on the next step.
      real(wp) :: hscal   = zero !! step size `h` used in scaling of nordsieck array `yh` in [[dvjust]].
                                 !! (if `iord = +1`, [[dvjust]] assumes that `hscal = tau(1)`.)
                                 !! see references 1 and 2 for details.
      real(wp) :: prl1    = zero !! the saved value of rl1.
      real(wp) :: rc      = zero !! ratio of current h*rl1 to value on last dvjac call.
      real(wp) :: rl1     = zero !! the reciprocal of the coefficient el(1).
      real(wp) :: tau(13) = zero !! real vector of past nq step sizes, length 13.
      real(wp) :: tq(5)   = zero !! a real vector of length 5 in which dvset stores constants
                                 !! used for the convergence test, the error test, and the
                                 !! selection of h at a new order.
      real(wp) :: tn      = zero !! the independent variable, updated on each step taken.
      real(wp) :: uround  = zero !! the machine unit roundoff.  the smallest positive real number
                                 !! such that  1.0 + uround /= 1.0
      integer :: icf      = 0 !! integer flag for convergence failure in dvnlsd:
                              !!
                              !!  * 0 means no failures.
                              !!  * 1 means convergence failure with out of date jacobian
                              !!    (recoverable error).
                              !!  * 2 means convergence failure with current jacobian or
                              !!    singular matrix (unrecoverable error).
      integer :: init     = 0 !! saved integer flag indicating whether initialization of the
                              !! problem has been done (init = 1) or not.
      integer :: ipup     = 0 !! saved flag to signal updating of newton matrix.
                              !!
                              !! In [[dvnlsd]]: own variable flag with values and meanings as follows:
                              !!
                              !!  * 0, do not update the newton matrix.
                              !!  * `miter /= 0`, update newton matrix, because it is the
                              !!    initial step, order was changed, the error
                              !!    test failed, or an update is indicated by
                              !!    the scalar `rc` or step counter `nst`.
      integer :: jcur     = 0 !! output flag from [[dvjac]] showing jacobian status:
                              !!
                              !!  * jcur = 0 means j is not current.
                              !!  * jcur = 1 means j is current.
      integer :: jstart   = 0 !! integer flag used as input to dvstep:
                              !!
                              !!  * 0  means perform the first step.
                              !!  * 1  means take a new step continuing from the last.
                              !!  * -1 means take the next step with a new value of maxord,
                              !!    hmin, hmxi, n, meth, miter, and/or matrix parameters.
                              !!
                              !! on return, dvstep sets jstart = 1.
      integer :: jsv      = 0 !! integer flag for jacobian saving, = sign(mf).
      integer :: kflag    = 0 !! a completion code from dvstep with the following meanings:
                              !!
                              !!  *  0      the step was succesful.
                              !!  * -1      the requested error could not be achieved.
                              !!  * -2      corrector convergence could not be achieved.
                              !!  * -3, -4  fatal error in vnls (can not occur here).
      integer :: kuth     = 0 !! input flag to dvstep showing whether h was reduced by the
                              !! driver.  kuth = 1 if h was reduced, = 0 otherwise.
      integer :: l        = 0 !! integer variable, nq + 1, current order plus one.
      integer :: lmax     = 0 !! maxord + 1 (used for dimensioning).
      integer :: lyh      = 0 !! saved integer pointers to segments of rwork and iwork.
      integer :: lewt     = 0 !! saved integer pointers to segments of rwork and iwork.
      integer :: lacor    = 0 !! saved integer pointers to segments of rwork and iwork.
      integer :: lsavf    = 0 !! saved integer pointers to segments of rwork and iwork.
      integer :: lwm      = 0 !! saved integer pointers to segments of rwork and iwork.
      integer :: liwm     = 0 !! saved integer pointers to segments of rwork and iwork.
      integer :: locjs    = 0 !! a pointer to the saved jacobian, whose storage starts at
                              !! wm(locjs), if jsv = 1.
      integer :: maxord   = 0 !! the maximum order of integration method to be allowed.
      integer :: meth     = 0 !! the method flags.  see mf.
      integer :: miter    = 0 !! the method flags.  see mf.
      integer :: msbj     = 0 !! the maximum number of steps between j evaluations, = 50.
      integer :: mxhnil   = 0 !! saved value of optional input mxhnil.
      integer :: mxstep   = 0 !! saved value of optional input mxstep.
      integer :: n        = 0 !! the number of first-order odes, = neq.
      integer :: newh     = 0 !! saved integer to flag change of h.
      integer :: newq     = 0 !! the method order to be used on the next step.
      integer :: nhnil    = 0 !! saved counter for occurrences of t + h = t.
      integer :: nq       = 0 !! integer variable, the current integration method order.
      integer :: nqnyh    = 0 !! saved value of nq*nyh.
      integer :: nqwait   = 0 !! a counter controlling the frequency of order changes.
                              !! an order change is about to be considered if nqwait = 1.
      integer :: nslj     = 0 !! the number of steps taken as of the last jacobian update.
      integer :: nslp     = 0 !! saved value of nst as of last newton matrix update.
      integer :: nyh      = 0 !! saved value of the initial value of neq.

      ! type declarations formerly in common block dvod02
      real(wp) :: hu  = zero !! the step size in t last used.
      integer :: ncfn = 0 !! number of nonlinear convergence failures so far.
      integer :: netf = 0 !! the number of error test failures of the integrator so far.
      integer :: nfe  = 0 !! the number of f evaluations for the problem so far.
      integer :: nje  = 0 !! the number of jacobian evaluations so far.
      integer :: nlu  = 0 !! the number of matrix lu decompositions so far.
      integer :: nni  = 0 !! number of nonlinear iterations so far.
      integer :: nqu  = 0 !! the method order last used.
      integer :: nst  = 0 !! the number of steps taken for the problem so far.

   end type  dvode_data_t

!*****************************************************************************************
!>
!  Main DVODE class.
!
!### Summary of usage.
!
! communication between the user and the dvode package, for normal
! situations, is summarized here.  this summary describes only a subset
! of the full set of options available.  see the full description under
! [[dvode]] for details, including optional communication, nonstandard
! options,  and instructions for special situations.
!
! 1. first provide a subroutine `f` of the form [[f_func]]:
!    which supplies the vector function f by loading ydot(i) with f(i).
!
! 2. next determine (or guess) whether or not the problem is stiff.
!    stiffness occurs when the jacobian matrix df/dy has an eigenvalue
!    whose real part is negative and large in magnitude, compared to the
!    reciprocal of the t span of interest.  if the problem is nonstiff,
!    use a method flag mf = 10.  if it is stiff, there are four standard
!    choices for mf (21, 22, 24, 25), and dvode requires the jacobian
!    matrix in some form.  in these cases (mf > 0), dvode will use a
!    saved copy of the jacobian matrix.  if this is undesirable because of
!    storage limitations, set mf to the corresponding negative value
!    (-21, -22, -24, -25).  (see full description of mf below.)
!    the jacobian matrix is regarded either as full (mf = 21 or 22),
!    or banded (mf = 24 or 25).  in the banded case, dvode requires two
!    half-bandwidth parameters ml and mu.  these are, respectively, the
!    widths of the lower and upper parts of the band, excluding the main
!    diagonal.  thus the band consists of the locations (i,j) with
!    i-ml <= j <= i+mu, and the full bandwidth is ml+mu+1.
!
! 3. if the problem is stiff, you are encouraged to supply the jacobian
!    directly (mf = 21 or 24), but if this is not feasible, dvode will
!    compute it internally by difference quotients (mf = 22 or 25).
!    if you are supplying the jacobian, provide a subroutine of the form [[f_jac]]
!    which supplies df/dy by loading pd as follows (in either case, only nonzero elements need be loaded):
!
!      * for a full jacobian (mf = 21), load pd(i,j) with df(i)/dy(j),
!        the partial derivative of f(i) with respect to y(j).  (ignore the
!        ml and mu arguments in this case.)
!      * for a banded jacobian (mf = 24), load pd(i-j+mu+1,j) with
!        df(i)/dy(j), i.e. load the diagonal lines of df/dy into the rows of
!        pd from the top down.
!
! 4. write a main program which calls subroutine dvode once for
!    each point at which answers are desired.  this should also provide
!    for possible use of logical unit 6 for output of error messages
!    by dvode.  on the first call to dvode, supply arguments as follows:
!
!      * f      = name of subroutine for right-hand side vector f.
!      * jac    = name of subroutine for jacobian matrix (mf = 21 or 24).
!      * neq    = number of first order odes.
!      * y      = array of initial values, of length neq.
!      * t      = the initial value of the independent variable.
!      * tout   = first point where output is desired (/= t).
!      * itol   = 1 or 2 according as atol (below) is a scalar or array.
!      * rtol   = relative tolerance parameter.
!      * atol   = absolute tolerance parameter.
!      * itask  = 1 for normal computation of output values of y at t = tout.
!      * istate = integer flag (input and output).  set istate = 1.
!      * iopt   = 0 to indicate no optional input used.
!      * rwork  = real work array
!      * lrw    = declared length of rwork.
!      * iwork  = integer work array.
!      * liw    = declared length of iwork (in user's dimension statement).
!      * mf     = method flag.
!
! 5. the output from the first call (or any call) is:
!    * y = array of computed values of y(t) vector.
!    * t = corresponding value of independent variable (normally tout).
!    * istate = 2  if dvode was successful, negative otherwise.
!        * -1 means excess work done on this call. (perhaps wrong mf.)
!        * -2 means excess accuracy requested. (tolerances too small.)
!        * -3 means illegal input detected. (see printed message.)
!        * -4 means repeated error test failures. (check all input.)
!        * -5 means repeated convergence failures. (perhaps bad
!             jacobian supplied or wrong choice of mf or tolerances.)
!        * -6 means error weight became zero during problem. (solution
!             component i vanished, and atol or atol(i) = 0.)
!
! 6. to continue the integration after a successful return, simply
!    reset tout and call [[dvode]] again.  no other parameters need be reset.
!
!### Other routines callable
!
! the following are optional calls which the user may make to
! gain additional capabilities in conjunction with dvode.
! (the routines [[xsetun]] and [[xsetf]] are designed to conform to the
! slatec error handling package.)
!
!  * [[xsetun]]
!  * [[xsetf]]
!  * [[dvsrco]]
!  * [[dvindy]]

   type,public :: dvode_t

      private

      type(dvode_data_t) :: dat !! internal data (formerly in common blocks)

      ! formerly save variables
      real(wp) :: etaq   = zero
      real(wp) :: etaqm1 = zero
      integer :: lunit = -1 !! logical unit number for messages.  the default is obtained
                            !! by a call to iumach (may be machine-dependent).
      integer :: mesflg = 1 !! print control flag:
                            !!
                            !!  * 1 means print all messages (the default).
                            !!  * 0 means no printing.

      procedure(f_func),pointer :: f => null()
      procedure(f_jac),pointer :: jac => null()
      procedure(f_dewset),pointer :: dewset => null()
      procedure(f_dvnorm),pointer :: dvnorm => null()

      contains

      private

      procedure,public :: initialize !! this routine must be called first before using the solver.
      procedure,public :: solve => dvode !! main solver routine.
      procedure,public :: xsetun !! set the logical unit number, lun, for
                                 !! output of messages from dvode, if
                                 !! the default is not desired.
                                 !! the default value of lun is 6.
                                 !! this call may be made at any time and will take effect immediately.
      procedure,public :: xsetf !! set a flag to control the printing of
                                !! messages by dvode.
                                !! mflag = 0 means do not print. (danger:
                                !! this risks losing valuable information.)
                                !! mflag = 1 means print (the default).
                                !! this call may be made at any time and will take effect immediately.
      procedure,public :: dvsrco !! saves and restores the contents of
                                 !! the internal variables used by dvode.
                                 !! [[dvsrco]] is useful if one is
                                 !! interrupting a run and restarting
                                 !! later, or alternating between two or
                                 !! more problems solved with [[dvode]].
      procedure,public :: dvindy !! provide derivatives of y, of various
                                 !! orders, at a specified point t, if
                                 !! desired.  it may be called only after
                                 !! a successful return from dvode.

      procedure :: dvhin
      procedure :: dvstep
      procedure :: dvset
      procedure :: dvjust
      procedure :: dvnlsd
      procedure :: dvjac
      procedure :: dvsol
      procedure :: ixsav
      procedure :: xerrwd

    end type dvode_t

    abstract interface
      subroutine f_func (me, neq, t, y, ydot)
         import :: dvode_t, wp
         implicit none
         class(dvode_t),intent(inout) :: me
         integer :: neq
         real(wp) :: t
         real(wp) :: y(neq)
         real(wp) :: ydot(neq)
      end subroutine f_func
      subroutine f_jac (me, neq, t, y, ml, mu, pd, nrowpd)
         import :: dvode_t, wp
         implicit none
         class(dvode_t),intent(inout) :: me
         integer :: neq
         integer :: nrowpd
         real(wp) :: t, y(neq), pd(nrowpd,neq)
         integer :: ml, mu
      end subroutine f_jac
      subroutine f_dewset(me,n,itol,rtol,atol,ycur,ewt)
         import :: dvode_t, wp
         implicit none
         class(dvode_t),intent(inout) :: me
         integer,intent(in) :: n
         integer,intent(in) :: itol
         real(wp),intent(in) :: rtol(*)
         real(wp),intent(in) :: atol(*)
         real(wp),intent(in) :: ycur(n)
         real(wp),intent(out) :: ewt(n)
      end subroutine f_dewset
      real(wp) function f_dvnorm(me,n,v,w)
         import :: dvode_t, wp
         implicit none
         class(dvode_t),intent(inout) :: me
         integer,intent(in) :: n
         real(wp),intent(in) :: v(n)
         real(wp),intent(in) :: w(n)
      end function f_dvnorm
   end interface

   public :: dvode

contains

!*****************************************************************************************
!>
!  Set the function pointers.
!
!### optionally replaceable solver routines
!
!  The two routines (`dewset` and `dvnorm`) relate to the measurement of errors.
!  either routine can be replaced by a user-supplied version, if desired.
!  however, since such a replacement may have a major impact on performance,
!  it should be done only when absolutely necessary, and only with great caution.
!  (note: the means by which the package version of a routine is
!  superseded by the user's version may be system-dependent.)
!  If not changed by the users, the defaults are used:
!
!   * [[dewset_default]]
!   * [[dvnorm_default]]

   subroutine initialize(me,f,jac,dewset,dvnorm)

   implicit none

   class(dvode_t),intent(inout) :: me
   procedure(f_func) :: f !! the name of the user-supplied subroutine defining the
                          !! ode system.  the system must be put in the first-order
                          !! form dy/dt = f(t,y), where f is a vector-valued function
                          !! of the scalar t and the vector y.  subroutine f is to
                          !! compute the function f.  it is to have the form
                          !!      subroutine f (neq, t, y, ydot)
                          !!      real(wp) t, y(neq), ydot(neq)
                          !! where neq, t, and y are input, and the array ydot = f(t,y)
                          !! is output.  y and ydot are arrays of length neq.
                          !! subroutine f should not alter y(1),...,y(neq).
                          !! f must be declared external in the calling program.
                          !!
                          !! if quantities computed in the f routine are needed
                          !! externally to dvode, an extra call to f should be made
                          !! for this purpose, for consistent and accurate results.
                          !! if only the derivative dy/dt is needed, use [[dvindy]] instead.
   procedure(f_jac),optional :: jac !! the name of the user-supplied routine (miter = 1 or 4) to
                                    !! compute the jacobian matrix, df/dy, as a function of
                                    !! the scalar t and the vector y.  it is to have the form [[f_jac]],
                                    !! where neq, t, y, ml, mu, and nrowpd are input and the array
                                    !! pd is to be loaded with partial derivatives (elements of the
                                    !! jacobian matrix) in the output.  pd must be given a first
                                    !! dimension of nrowpd.  t and y have the same meaning as in
                                    !! subroutine f.
                                    !!
                                    !!  * in the full matrix case (miter = 1), ml and mu are
                                    !!    ignored, and the jacobian is to be loaded into pd in
                                    !!    columnwise manner, with df(i)/dy(j) loaded into pd(i,j).
                                    !!
                                    !!  * in the band matrix case (miter = 4), the elements
                                    !!    within the band are to be loaded into pd in columnwise
                                    !!    manner, with diagonal lines of df/dy loaded into the rows
                                    !!    of pd. thus df(i)/dy(j) is to be loaded into pd(i-j+mu+1,j).
                                    !!    ml and mu are the half-bandwidth parameters. (see iwork).
                                    !!    the locations in pd in the two triangular areas which
                                    !!    correspond to nonexistent matrix elements can be ignored
                                    !!    or loaded arbitrarily, as they are overwritten by dvode.
                                    !!
                                    !! jac need not provide df/dy exactly.  a crude
                                    !! approximation (possibly with a smaller bandwidth) will do.
                                    !!
                                    !! in either case, pd is preset to zero by the solver,
                                    !! so that only the nonzero elements need be loaded by jac.
                                    !! each call to jac is preceded by a call to f with the same
                                    !! arguments neq, t, and y.  thus to gain some efficiency,
                                    !! intermediate quantities shared by both calculations may be
                                    !! saved in a user common block by f and not recomputed by jac,
                                    !! if desired.  also, jac may alter the y array, if desired.
                                    !! jac must be declared external in the calling program.
   procedure(f_dewset),optional :: dewset !! the following subroutine is called just before each internal
                                          !! integration step, and sets the array of error weights, ewt, as
                                          !! described under itol/rtol/atol above:
                                          !! `subroutine dewset (neq, itol, rtol, atol, ycur, ewt)`
                                          !! where neq, itol, rtol, and atol are as in the [[dvode]] call sequence,
                                          !! ycur contains the current dependent variable vector, and
                                          !! ewt is the array of weights set by dewset.
                                          !!
                                          !! if the user supplies this subroutine, it must return in ewt(i)
                                          !! (i = 1,...,neq) a positive quantity suitable for comparison with
                                          !! errors in y(i).  the ewt array returned by dewset is passed to the
                                          !! dvnorm routine (see below.), and also used by dvode in the computation
                                          !! of the optional output imxer, the diagonal jacobian approximation,
                                          !! and the increments for difference quotient jacobians.
                                          !!
                                          !! in the user-supplied version of dewset, it may be desirable to use
                                          !! the current values of derivatives of y.  derivatives up to order nq
                                          !! are available from the history array yh, described above under
                                          !! optional output.  in dewset, yh is identical to the ycur array,
                                          !! extended to nq + 1 columns with a column length of nyh and scale
                                          !! factors of h**j/factorial(j).  on the first call for the problem,
                                          !! given by nst = 0, nq is 1 and h is temporarily set to 1.0.
                                          !! nyh is the initial value of neq.  the quantities nq, h, and nst
                                          !! can be obtained by including in dewset the statements:
                                          !!```fortran
                                          !!  type(dvode_data_t) :: sav
                                          !!  call me%dvsrco(sav,job=1)
                                          !!  nq = sav%nq
                                          !!  h = sav%h
                                          !!  nst = sav%nst
                                          !!```
                                          !! thus, for example, the current value of dy/dt can be obtained as
                                          !! `ycur(nyh+i)/h  (i=1,...,neq)`  (and the division by h is
                                          !! unnecessary when nst = 0).
   procedure(f_dvnorm),optional :: dvnorm !! the following is a real function routine which computes the weighted
                                          !! root-mean-square norm of a vector v:
                                          !! `d = dvnorm (n, v, w)`
                                          !! where:
                                          !!
                                          !!   * n = the length of the vector,
                                          !!   * v = real array of length n containing the vector,
                                          !!   * w = real array of length n containing weights,
                                          !!   * d = `sqrt( (1/n) * sum(v(i)*w(i))**2 )`.
                                          !!
                                          !! dvnorm is called with n = neq and with w(i) = 1.0/ewt(i), where
                                          !! ewt is as set by subroutine dewset.
                                          !!
                                          !! if the user supplies this function, it should return a non-negative
                                          !! value of dvnorm suitable for use in the error control in dvode.
                                          !! none of the arguments should be altered by dvnorm.
                                          !! for example, a user-supplied dvnorm routine might:
                                          !!
                                          !!   * substitute a max-norm of (v(i)*w(i)) for the rms-norm, or
                                          !!   * ignore some components of v in the norm, with the effect of
                                          !!     suppressing the error control on those components of y.

   me%f => f

   if (present(jac)) then
      me%jac => jac
   else
      me%jac => null()
   end if
   if (present(dewset)) then
      me%dewset => dewset
   else
      me%dewset => dewset_default
   end if
   if (present(dvnorm)) then
      me%dvnorm => dvnorm
   else
      me%dvnorm => dvnorm_default
   end if

   end subroutine initialize
!*****************************************************************************************

!*****************************************************************************************
!>
!  dvode: variable-coefficient ordinary differential equation solver,
!  with fixed-leading-coefficient implementation.
!
!  dvode solves the initial value problem for stiff or nonstiff
!  systems of first order odes,
!```
!  dy/dt = f(t,y) ,  or, in component form,
!  dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(neq)) (i = 1,...,neq).
!```
!  dvode is a package based on the `episode` and `episodeb` packages, and
!  on the `odepack` user interface standard, with minor modifications.
!
!### Optional input.
!
! the work arrays `rwork` and `iwork` are also used for conditional and
! optional input and optional output.  (the term output here refers
! to the return from subroutine dvode to the user's calling program.)
!
! the following is a list of the optional input provided for in the
! call sequence.  (see also part ii.)  for each such input variable,
! this table lists its name as used in this documentation, its
! location in the call sequence, its meaning, and the default value.
! the use of any of this input requires `iopt = 1`, and in that
! case all of this input is examined.  a value of zero for any
! of these optional input variables will cause the default value to be
! used.  thus to use a subset of the optional input, simply preload
! locations 5 to 10 in rwork and iwork to 0.0 and 0 respectively, and
! then set those of interest to nonzero values.
!```
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
!```
!### Optional output
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
!```
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
!```
! the following two arrays are segments of the rwork array which
! may also be of interest to the user as optional output.
! for each array, the table below gives its internal name,
! its base address in rwork, and its description.
!```
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
!```
!
!### Interrupting and restarting
!
! if the integration of a given problem by dvode is to be
! interrrupted and then later continued, such as when restarting
! an interrupted run or alternating between two or more ode problems,
! the user should save, following the return from the last dvode call
! prior to the interruption, the contents of the call sequence
! variables and internal variables, and later restore these
! values before the next dvode call for that problem.  to save
! and restore the variables, use subroutine [[dvsrco]], as
! described below in part ii.
!
! in addition, if non-default values for either lun or mflag are
! desired, an extra call to [[xsetun]] and/or [[xsetf]] should be made just
! before continuing the integration.  see part ii below for details.
!
!### Authors:
!  * peter n. brown and alan c. hindmarsh,
!    center for applied scientific computing, l-561,
!    lawrence livermore national laboratory,
!    livermore, ca 94551
!  * george d. byrne,
!    illinois institute of technology,
!    chicago, il 60616
!
!### References
!
!  1. p. n. brown, g. d. byrne, and a. c. hindmarsh, "vode: a variable
!     coefficient ode solver," siam j. sci. stat. comput., 10 (1989),
!     pp. 1038-1051.  also, llnl report ucrl-98412, june 1988.
!  2. g. d. byrne and a. c. hindmarsh, "a polyalgorithm for the
!     numerical solution of ordinary differential equations,"
!     acm trans. math. software, 1 (1975), pp. 71-96.
!  3. a. c. hindmarsh and g. d. byrne, "episode: an effective package
!     for the integration of systems of ordinary differential
!     equations," llnl report ucid-30112, rev. 1, april 1977.
!  4. g. d. byrne and a. c. hindmarsh, "episodeb: an experimental
!     package for the integration of systems of ordinary differential
!     equations with banded jacobians," llnl report ucid-30132, april
!     1976.
!  5. a. c. hindmarsh, "odepack, a systematized collection of ode
!     solvers," in scientific computing, r. s. stepleman et al., eds.,
!     north-holland, amsterdam, 1983, pp. 55-64.
!  6. k. r. jackson and r. sacks-davis, "an alternative implementation
!     of variable step-size multistep formulas for stiff odes," acm
!     trans. math. software, 6 (1980), pp. 295-318.
!
!### Revision history
!  * 19890615  date written.  initial release.
!  * 19890922  added interrupt/restart ability, minor changes throughout.
!  * 19910228  minor revisions in line format,  prologue, etc.
!  * 19920227  modifications by d. pang:
!              (1) applied subgennam to get generic intrinsic names.
!              (2) changed intrinsic names to generic in comments.
!              (3) added *deck lines before each routine.
!  * 19920721  names of routines and labeled common blocks changed, so as
!              to be unique in combined single/real(wp) code (ach).
!  * 19920722  minor revisions to prologue (ach).
!  * 19920831  conversion to real(wp) done (ach).
!  * 19921106  fixed minor bug: etaq,etaqm1 in dvstep save statement (ach).
!  * 19921118  changed lunsav/mflgsv to ixsav (ach).
!  * 19941222  removed mf overwrite; attached sign to h in estimated second
!              deriv. in dvhin; misc. comment changes throughout (ach).
!  * 19970515  minor corrections to comments in prologue, dvjac (ach).
!  * 19981111  corrected block b by adding final line, go to 200 (ach).
!  * 20020430  various upgrades (ach): use odepack error handler package.
!              replaced d1mach by dumach.  various changes to main
!              prologue and other routine prologues.
!
!@note the legality of input parameters will be thoroughly checked on the
!      initial call for the problem, but not checked thereafter unless a
!      change in input parameters is flagged by istate = 3 in the input.
!
!@note  the work arrays must not be altered between calls to dvode
!       for the same problem, except possibly for the conditional and
!       optional input, and except for the last 3*neq words of rwork.
!       the latter space is used for internal scratch space, and so is
!       available for use by the user outside dvode between calls, if
!       desired (but not for use by f or jac).

   subroutine dvode(me,neq,y,t,tout,itol,rtol,atol,itask,istate,iopt, &
                    rwork,lrw,iwork,liw,mf)

      class(dvode_t),intent(inout) :: me
      real(wp),intent(inout) :: y(*) !! a real array for the vector of dependent variables, of
                                     !! length `neq` or more.  used for both input and output on the
                                     !! first call (`istate = 1`), and only for output on other calls.
                                     !! on the first call, `y` must contain the vector of initial
                                     !! values.  in the output, `y` contains the computed solution
                                     !! evaluated at `t`.  if desired, the `y` array may be used
                                     !! for other purposes between calls to the solver.
                                     !!
                                     !! this array is passed as the `y` argument in all calls to
                                     !! `f` and `jac`.
      real(wp),intent(inout) :: t !! the independent variable.  in the input, `t` is used only on
                                  !! the first call, as the initial point of the integration.
                                  !! in the output, after each call, `t` is the value at which a
                                  !! computed solution `y` is evaluated (usually the same as `tout`).
                                  !! on an error return, `t` is the farthest point reached.
      real(wp),intent(in) :: tout !! the next value of `t` at which a computed solution is desired.
                                  !!
                                  !! when starting the problem (`istate = 1`), `tout` may be equal
                                  !! to `t` for one call, then should `/= t` for the next call.
                                  !! for the initial `t`, an input value of `tout /= t` is used
                                  !! in order to determine the direction of the integration
                                  !! (i.e. the algebraic sign of the step sizes) and the rough
                                  !! scale of the problem.  integration in either direction
                                  !! (forward or backward in `t`) is permitted.
                                  !!
                                  !! if `itask = 2 or 5` (one-step modes), `tout` is ignored after
                                  !! the first call (i.e. the first call with `tout /= t`).
                                  !! otherwise, `tout` is required on every call.
                                  !!
                                  !! if `itask = 1, 3, or 4`, the values of `tout` need not be
                                  !! monotone, but a value of `tout` which backs up is limited
                                  !! to the current internal `t` interval, whose endpoints are
                                  !! `tcur - hu` and `tcur`.  (see optional output, below, for
                                  !! `tcur` and `hu`.)
      real(wp),intent(in) :: rtol(*) !! a relative error tolerance parameter, either a scalar or
                                     !! an array of length `neq`.  see description under `atol`.
      real(wp),intent(in) :: atol(*) !! an absolute error tolerance parameter, either a scalar or
                                     !! an array of length `neq`.
                                     !!
                                     !! the input parameters itol, rtol, and atol determine
                                     !! the error control performed by the solver.  the solver will
                                     !! control the vector e = (e(i)) of estimated local errors
                                     !! in y, according to an inequality of the form
                                     !!```
                                     !!             rms-norm of ( e(i)/ewt(i) )   <=   1,
                                     !! where       ewt(i) = rtol(i)*abs(y(i)) + atol(i),
                                     !!```
                                     !! and the rms-norm (root-mean-square norm) here is
                                     !! rms-norm(v) = sqrt(sum v(i)**2 / neq).  here ewt = (ewt(i))
                                     !! is a vector of weights which must always be positive, and
                                     !! the values of rtol and atol should all be non-negative.
                                     !! the following table gives the types (scalar/array) of
                                     !! rtol and atol, and the corresponding form of ewt(i).
                                     !!
                                     !!```
                                     !!    itol    rtol       atol          ewt(i)
                                     !!     1     scalar     scalar     rtol*abs(y(i)) + atol
                                     !!     2     scalar     array      rtol*abs(y(i)) + atol(i)
                                     !!     3     array      scalar     rtol(i)*abs(y(i)) + atol
                                     !!     4     array      array      rtol(i)*abs(y(i)) + atol(i)
                                     !!```
                                     !!
                                     !! when either of these parameters is a scalar, it need not
                                     !! be dimensioned in the user's calling program.
                                     !!
                                     !! if none of the above choices (with itol, rtol, and atol
                                     !! fixed throughout the problem) is suitable, more general
                                     !! error controls can be obtained by substituting
                                     !! user-supplied routines for the setting of ewt and/or for
                                     !! the norm calculation.  see part iv below.
                                     !!
                                     !! if global errors are to be estimated by making a repeated
                                     !! run on the same problem with smaller tolerances, then all
                                     !! components of rtol and atol (i.e. of ewt) should be scaled
                                     !! down uniformly.
                                     !!
                                     !! use rtol = 0.0 for pure absolute error control, and
                                     !! use atol = 0.0 (or atol(i) = 0.0) for pure relative error
                                     !! control.  caution: actual (global) errors may exceed these
                                     !! local tolerances, so choose them conservatively.
      integer,intent(in) :: lrw !! the length of the array rwork, as declared by the user.
                                !! (this will be checked by the solver.)
      real(wp) :: rwork(lrw) !! a real working array.
                             !! the length of `rwork` must be at least
                             !! `20 + nyh*(maxord + 1) + 3*neq + lwm` where:
                             !!
                             !!  * nyh    = the initial value of neq,
                             !!  * maxord = 12 (if meth = 1) or 5 (if meth = 2) (unless a
                             !!             smaller value is given as an optional input),
                             !!  * lwm = length of work space for matrix-related data:
                             !!
                             !!      * `lwm = 0                    ` if miter = 0,
                             !!      * `lwm = 2*neq**2 + 2         ` if miter = 1 or 2, and mf>0,
                             !!      * `lwm = neq**2 + 2           ` if miter = 1 or 2, and mf<0,
                             !!      * `lwm = neq + 2              ` if miter = 3,
                             !!      * `lwm = (3*ml+2*mu+2)*neq + 2` if miter = 4 or 5, and mf>0,
                             !!      * `lwm = (2*ml+mu+1)*neq + 2  ` if miter = 4 or 5, and mf<0.
                             !!
                             !! (see the mf description for meth and miter.)
                             !! thus if maxord has its default value and neq is constant,
                             !! this length is:
                             !!
                             !!  * `20 + 16*neq                  `  for mf = 10,
                             !!  * `22 + 16*neq + 2*neq**2       `  for mf = 11 or 12,
                             !!  * `22 + 16*neq + neq**2         `  for mf = -11 or -12,
                             !!  * `22 + 17*neq                  `  for mf = 13,
                             !!  * `22 + 18*neq + (3*ml+2*mu)*neq`  for mf = 14 or 15,
                             !!  * `22 + 17*neq + (2*ml+mu)*neq  `  for mf = -14 or -15,
                             !!  * `20 +  9*neq                  `  for mf = 20,
                             !!  * `22 +  9*neq + 2*neq**2       `  for mf = 21 or 22,
                             !!  * `22 +  9*neq + neq**2         `  for mf = -21 or -22,
                             !!  * `22 + 10*neq                  `  for mf = 23,
                             !!  * `22 + 11*neq + (3*ml+2*mu)*neq`  for mf = 24 or 25.
                             !!  * `22 + 10*neq + (2*ml+mu)*neq  `  for mf = -24 or -25.
                             !!
                             !! the first 20 words of rwork are reserved for conditional
                             !! and optional input and optional output.
                             !!
                             !! the following word in rwork is a conditional input:
                             !!   rwork(1) = tcrit = critical value of t which the solver
                             !!              is not to overshoot.  required if itask is
                             !!              4 or 5, and ignored otherwise.  (see itask.)




      integer,intent(in) :: neq !! the size of the ode system (number of first order
                                !! ordinary differential equations).
                                !! `neq` may not be increased during the problem, but
                                !! can be decreased (with `istate = 3` in the input).
      integer,intent(in) :: itol !! an indicator for the type of error control.
                                 !! 1 or 2 according as `atol` is a scalar or array.
                                 !! see description under `atol`.
      integer,intent(in) :: itask !! an index specifying the task to be performed.
                                  !! input only.  itask has the following values and meanings:
                                  !!
                                  !!  * 1  means normal computation of output values of y(t) at
                                  !!       t = tout (by overshooting and interpolating).
                                  !! * 2  means take one step only and return.
                                  !! * 3  means stop at the first internal mesh point at or
                                  !!      beyond t = tout and return.
                                  !! * 4  means normal computation of output values of y(t) at
                                  !!      t = tout but without overshooting t = tcrit.
                                  !!      tcrit must be input as rwork(1).  tcrit may be equal to
                                  !!      or beyond tout, but not behind it in the direction of
                                  !!      integration.  this option is useful if the problem
                                  !!      has a singularity at or beyond t = tcrit.
                                  !! * 5  means take one step, without passing tcrit, and return.
                                  !!      tcrit must be input as rwork(1).
                                  !!
                                  !! note:  if itask = 4 or 5 and the solver reaches tcrit
                                  !! (within roundoff), it will return t = tcrit (exactly) to
                                  !! indicate this (unless itask = 4 and tout comes before tcrit,
                                  !! in which case answers at t = tout are returned first).
      integer,intent(inout) :: istate !! an index used for input and output to specify the
                                      !! the state of the calculation.
                                      !!
                                      !! in the input, the values of istate are as follows:
                                      !!
                                      !!  1.  means this is the first call for the problem
                                      !!      (initializations will be done).  see note below.
                                      !!  2.  means this is not the first call, and the calculation
                                      !!      is to continue normally, with no change in any input
                                      !!      parameters except possibly tout and itask.
                                      !!      (if itol, rtol, and/or atol are changed between calls
                                      !!      with istate = 2, the new values will be used but not
                                      !!      tested for legality.)
                                      !!  3.  means this is not the first call, and the
                                      !!      calculation is to continue normally, but with
                                      !!      a change in input parameters other than
                                      !!      tout and itask.  changes are allowed in
                                      !!      neq, itol, rtol, atol, iopt, lrw, liw, mf, ml, mu,
                                      !!      and any of the optional input except h0.
                                      !!      (see iwork description for ml and mu.)
                                      !!
                                      !! note:  a preliminary call with tout = t is not counted
                                      !! as a first call here, as no initialization or checking of
                                      !! input is done.  (such a call is sometimes useful to include
                                      !! the initial conditions in the output.)
                                      !! thus the first call for which tout /= t requires
                                      !! istate = 1 in the input.
                                      !!
                                      !! in the output, istate has the following values and meanings:
                                      !!
                                      !!  * 1  means nothing was done, as tout was equal to t with
                                      !!      istate = 1 in the input.
                                      !!  * 2  means the integration was performed successfully.
                                      !!  * -1  means an excessive amount of work (more than mxstep
                                      !!       steps) was done on this call, before completing the
                                      !!       requested task, but the integration was otherwise
                                      !!       successful as far as t.  (mxstep is an optional input
                                      !!       and is normally 500.)  to continue, the user may
                                      !!       simply reset istate to a value > 1 and call again.
                                      !!       (the excess work step counter will be reset to 0.)
                                      !!       in addition, the user may increase mxstep to avoid
                                      !!       this error return.  (see optional input below.)
                                      !!  * -2  means too much accuracy was requested for the precision
                                      !!       of the machine being used.  this was detected before
                                      !!       completing the requested task, but the integration
                                      !!       was successful as far as t.  to continue, the tolerance
                                      !!       parameters must be reset, and istate must be set
                                      !!       to 3.  the optional output tolsf may be used for this
                                      !!       purpose.  (note: if this condition is detected before
                                      !!       taking any steps, then an illegal input return
                                      !!       (istate = -3) occurs instead.)
                                      !!  * -3  means illegal input was detected, before taking any
                                      !!       integration steps.  see written message for details.
                                      !!       note:  if the solver detects an infinite loop of calls
                                      !!       to the solver with illegal input, it will cause
                                      !!       the run to stop.
                                      !!  * -4  means there were repeated error test failures on
                                      !!       one attempted step, before completing the requested
                                      !!       task, but the integration was successful as far as t.
                                      !!       the problem may have a singularity, or the input
                                      !!       may be inappropriate.
                                      !!  * -5  means there were repeated convergence test failures on
                                      !!       one attempted step, before completing the requested
                                      !!       task, but the integration was successful as far as t.
                                      !!       this may be caused by an inaccurate jacobian matrix,
                                      !!       if one is being used.
                                      !!  * -6  means ewt(i) became zero for some i during the
                                      !!       integration.  pure relative error control (atol(i)=0.0)
                                      !!       was requested on a variable which has now vanished.
                                      !!       the integration was successful as far as t.
                                      !!
                                      !! note:  since the normal output value of istate is 2,
                                      !! it does not need to be reset for normal continuation.
                                      !! also, since a negative input value of istate will be
                                      !! regarded as illegal, a negative output value requires the
                                      !! user to change it, and possibly other input, before
                                      !! calling the solver again.

      integer,intent(in) :: iopt !! an integer flag to specify whether or not any optional
                                 !! input is being used on this call.
                                 !! the optional input is listed separately below.
                                 !!
                                 !!  * iopt = 0 means no optional input is being used.
                                 !!    default values will be used in all cases.
                                 !!  * iopt = 1 means optional input is being used.

      integer,intent(in) :: liw !! the length of the array iwork, as declared by the user.
                                !! (this will be checked by the solver.)
      integer :: iwork(liw) !! an integer work array.  the length of iwork must be at least:
                            !!
                            !!  * `30`        if miter = 0 or 3 (mf = 10, 13, 20, 23), or
                            !!  * `30 + neq`  otherwise (abs(mf) = 11,12,14,15,21,22,24,25).
                            !!
                            !! the first 30 words of `iwork` are reserved for conditional and
                            !! optional input and optional output.
                            !!
                            !! the following 2 words in `iwork` are conditional input:
                            !!
                            !!  * `iwork(1) = ml`
                            !!  * `iwork(2) = mu`
                            !!
                            !! these are the lower and upper
                            !! half-bandwidths, respectively, of the
                            !! banded jacobian, excluding the main diagonal.
                            !! the band is defined by the matrix locations
                            !! `(i,j)` with `i-ml <= j <= i+mu`.  `ml` and `mu`
                            !! must satisfy  `0 <=  ml,mu  <= neq-1`.
                            !! these are required if `miter` is 4 or 5, and
                            !! ignored otherwise.  `ml` and `mu` may in fact be
                            !! the band parameters for a matrix to which
                            !! `df/dy` is only approximately equal.

      integer,intent(in) :: mf !! method flag.  standard values are:
                               !!
                               !!  * 10 for nonstiff (adams) method, no jacobian used.
                               !!  * 21 for stiff (bdf) method, user-supplied full jacobian.
                               !!  * 22 for stiff method, internally generated full jacobian.
                               !!  * 24 for stiff method, user-supplied banded jacobian.
                               !!  * 25 for stiff method, internally generated banded jacobian.
                               !!
                               !! the complete set of legal values of
                               !! mf are 10, 11, 12, 13, 14, 15, 20, 21, 22, 23, 24, 25,
                               !! -11, -12, -14, -15, -21, -22, -24, -25.
                               !!
                               !! `mf` is a signed two-digit integer, `mf = jsv*(10*meth + miter)`
                               !!
                               !! `jsv = sign(mf)` indicates the jacobian-saving strategy:
                               !!
                               !!  * jsv =  1 means a copy of the jacobian is saved for reuse
                               !!           in the corrector iteration algorithm.
                               !!  * jsv = -1 means a copy of the jacobian is not saved
                               !!          (valid only for miter = 1, 2, 4, or 5).
                               !!
                               !! `meth` indicates the basic linear multistep method:
                               !!
                               !!  * meth = 1 means the implicit adams method.
                               !!  * meth = 2 means the method based on backward
                               !!           differentiation formulas (bdf-s).
                               !!
                               !! `miter` indicates the corrector iteration method:
                               !!
                               !!  * miter = 0 means functional iteration (no jacobian matrix
                               !!    is involved).
                               !!  * miter = 1 means chord iteration with a user-supplied
                               !!    full (neq by neq) jacobian.
                               !!  * miter = 2 means chord iteration with an internally
                               !!    generated (difference quotient) full jacobian
                               !!    (using neq extra calls to f per df/dy value).
                               !!  * miter = 3 means chord iteration with an internally
                               !!    generated diagonal jacobian approximation
                               !!    (using 1 extra call to f per df/dy evaluation).
                               !!  * miter = 4 means chord iteration with a user-supplied
                               !!    banded jacobian.
                               !!  * miter = 5 means chord iteration with an internally
                               !!    generated banded jacobian (using ml+mu+1 extra
                               !!    calls to f per df/dy evaluation).
                               !!
                               !! if miter = 1 or 4, the user must supply a subroutine `jac`
                               !! (the name is arbitrary) as described above under `jac`.
                               !! for other values of miter, a dummy argument can be used.

      logical :: ihit
      real(wp) :: atoli , big , ewti , h0 , hmax , hmx , &
                  rh , rtoli , size , tcrit , &
                  tnext , tolsf , tp
      integer :: i , ier , iflag , imxer , jco , kgo , leniw , lenj , &
                 lenp , lenrw , lenwm , lf0 , mband , mfa , ml , &
                 mu , niter , nslast
      character(len=80) :: msg

      integer,parameter :: mxhnl0 = 10
      integer,parameter :: mxstp0 = 500
      real(wp),parameter :: one  = 1.0_wp
      real(wp),parameter :: two  = 2.0_wp
      real(wp),parameter :: four = 4.0_wp
      real(wp),parameter :: pt2  = 0.2_wp
      real(wp),parameter :: hun  = 100.0_wp

      !-----------------------------------------------------------------------
      ! block a.
      ! this code block is executed on every call.
      ! it tests istate and itask for legality and branches appropriately.
      ! if istate > 1 but the flag init shows that initialization has
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
         call me%xerrwd(msg,30,1,1,1,istate,0,0,zero,zero)
         if ( istate>=0 ) goto 1500
         msg = 'dvode--  run aborted:  apparent infinite loop     '
         call me%xerrwd(msg,50,303,2,0,0,0,0,zero,zero)
         return
      else
         if ( itask<1 .or. itask>5 ) then
            msg = 'dvode--  itask (=i1) illegal  '
            call me%xerrwd(msg,30,2,1,1,itask,0,0,zero,zero)
            goto 1500
         else
            if ( istate==1 ) then
               me%dat%init = 0
               if ( tout==t ) return
            else if ( me%dat%init/=1 ) then
               msg = 'dvode--  istate (=i1) > 1 but dvode not initialized      '
               call me%xerrwd(msg,60,3,1,1,istate,0,0,zero,zero)
               goto 1500
            else if ( istate==2 ) then
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
               msg = 'dvode--  neq (=i1) < 1     '
               call me%xerrwd(msg,30,4,1,1,neq,0,0,zero,zero)
               goto 1500
            else
               if ( istate/=1 ) then
                  if ( neq>me%dat%n ) then
                     msg = 'dvode--  istate = 3 and neq increased (i1 to i2)  '
                     call me%xerrwd(msg,50,5,1,2,me%dat%n,neq,0,zero,zero)
                     goto 1500
                  endif
               endif
               me%dat%n = neq
               if ( itol<1 .or. itol>4 ) then
                  msg = 'dvode--  itol (=i1) illegal   '
                  call me%xerrwd(msg,30,6,1,1,itol,0,0,zero,zero)
                  goto 1500
               else if ( iopt<0 .or. iopt>1 ) then
                  msg = 'dvode--  iopt (=i1) illegal   '
                  call me%xerrwd(msg,30,7,1,1,iopt,0,0,zero,zero)
                  goto 1500
               else
                  me%dat%jsv = sign(1,mf)
                  mfa = abs(mf)
                  me%dat%meth = mfa/10
                  me%dat%miter = mfa - 10*me%dat%meth
                  if ( me%dat%meth<1 .or. me%dat%meth>2 ) goto 800
                  if ( me%dat%miter<0 .or. me%dat%miter>5 ) goto 800
                  if ( me%dat%miter>3 ) then
                     ml = iwork(1)
                     mu = iwork(2)
                     if ( ml<0 .or. ml>=me%dat%n ) then
                        msg = 'dvode--  ml (=i1) illegal:  <0 or >=neq (=i2)'
                        call me%xerrwd(msg,50,9,1,2,ml,neq,0,zero,zero)
                        goto 1500
                     else if ( mu<0 .or. mu>=me%dat%n ) then
                        msg = 'dvode--  mu (=i1) illegal:  <0 or >=neq (=i2)'
                        call me%xerrwd(msg,50,10,1,2,mu,neq,0,zero,zero)
                        goto 1500
                     endif
                  endif
                  ! next process and check the optional input. ---------------------------
                  if ( iopt==1 ) then
                     me%dat%maxord = iwork(5)
                     if ( me%dat%maxord<0 ) then
                        msg = 'dvode--  maxord (=i1) < 0  '
                        call me%xerrwd(msg,30,11,1,1,me%dat%maxord,0,0,zero,zero)
                        goto 1500
                     else
                        if ( me%dat%maxord==0 ) me%dat%maxord = 100
                        me%dat%maxord = min(me%dat%maxord,mord(me%dat%meth))
                        me%dat%mxstep = iwork(6)
                        if ( me%dat%mxstep<0 ) then
                           msg = 'dvode--  mxstep (=i1) < 0  '
                           call me%xerrwd(msg,30,12,1,1,me%dat%mxstep,0,0,zero,zero)
                           goto 1500
                        else
                           if ( me%dat%mxstep==0 ) me%dat%mxstep = mxstp0
                           me%dat%mxhnil = iwork(7)
                           if ( me%dat%mxhnil<0 ) then
                              msg = 'dvode--  mxhnil (=i1) < 0  '
                              call me%xerrwd(msg,30,13,1,1,me%dat%mxhnil,0,0,zero,zero)
                              goto 1500
                           else
                              if ( me%dat%mxhnil==0 ) me%dat%mxhnil = mxhnl0
                              if ( istate==1 ) then
                                 h0 = rwork(5)
                                 if ( (tout-t)*h0<zero ) then
                                    msg = 'dvode--  tout (=r1) behind t (=r2)      '
                                    call me%xerrwd(msg,40,14,1,0,0,0,2,tout,t)
                                    msg = '      integration direction is given by h0 (=r1)  '
                                    call me%xerrwd(msg,50,14,1,0,0,0,1,h0,zero)
                                    goto 1500
                                 endif
                              endif
                              hmax = rwork(6)
                              if ( hmax<zero ) then
                                 msg = 'dvode--  hmax (=r1) < 0.0  '
                                 call me%xerrwd(msg,30,15,1,0,0,0,1,hmax,zero)
                                 goto 1500
                              else
                                 me%dat%hmxi = zero
                                 if ( hmax>zero ) me%dat%hmxi = one/hmax
                                 me%dat%hmin = rwork(7)
                                 if ( me%dat%hmin<zero ) then
                                    msg = 'dvode--  hmin (=r1) < 0.0  '
                                    call me%xerrwd(msg,30,16,1,0,0,0,1,me%dat%hmin,zero)
                                    goto 1500
                                 endif
                              endif
                           endif
                        endif
                     endif
                  else
                     me%dat%maxord = mord(me%dat%meth)
                     me%dat%mxstep = mxstp0
                     me%dat%mxhnil = mxhnl0
                     if ( istate==1 ) h0 = zero
                     me%dat%hmxi = zero
                     me%dat%hmin = zero
                  endif
                  !-----------------------------------------------------------------------
                  ! set work array pointers and check lengths lrw and liw.
                  ! pointers to segments of rwork and iwork are named by prefixing l to
                  ! the name of the segment.  e.g., the segment yh starts at rwork(lyh).
                  ! segments of rwork (in order) are denoted  yh, wm, ewt, savf, acor.
                  ! within wm, locjs is the location of the saved jacobian (jsv > 0).
                  !-----------------------------------------------------------------------
                  me%dat%lyh = 21
                  if ( istate==1 ) me%dat%nyh = me%dat%n
                  me%dat%lwm = me%dat%lyh + (me%dat%maxord+1)*me%dat%nyh
                  jco = max(0,me%dat%jsv)
                  if ( me%dat%miter==0 ) lenwm = 0
                  if ( me%dat%miter==1 .or. me%dat%miter==2 ) then
                     lenwm = 2 + (1+jco)*me%dat%n*me%dat%n
                     me%dat%locjs = me%dat%n*me%dat%n + 3
                  endif
                  if ( me%dat%miter==3 ) lenwm = 2 + me%dat%n
                  if ( me%dat%miter==4 .or. me%dat%miter==5 ) then
                     mband = ml + mu + 1
                     lenp = (mband+ml)*me%dat%n
                     lenj = mband*me%dat%n
                     lenwm = 2 + lenp + jco*lenj
                     me%dat%locjs = lenp + 3
                  endif
                  me%dat%lewt = me%dat%lwm + lenwm
                  me%dat%lsavf = me%dat%lewt + me%dat%n
                  me%dat%lacor = me%dat%lsavf + me%dat%n
                  lenrw = me%dat%lacor + me%dat%n - 1
                  iwork(17) = lenrw
                  me%dat%liwm = 1
                  leniw = 30 + me%dat%n
                  if ( me%dat%miter==0 .or. me%dat%miter==3 ) leniw = 30
                  iwork(18) = leniw
                  if ( lenrw>lrw ) then
                     msg = 'dvode--  rwork length needed, lenrw (=i1), exceeds lrw (=i2)'
                     call me%xerrwd(msg,60,17,1,2,lenrw,lrw,0,zero,zero)
                     goto 1500
                  else if ( leniw>liw ) then
                     msg = 'dvode--  iwork length needed, leniw (=i1), exceeds liw (=i2)'
                     call me%xerrwd(msg,60,18,1,2,leniw,liw,0,zero,zero)
                     goto 1500
                  else
                     ! check rtol and atol for legality. ------------------------------------
                     rtoli = rtol(1)
                     atoli = atol(1)
                     do i = 1 , me%dat%n
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
                        me%dat%uround = epmach
                        me%dat%tn = t
                        if ( itask==4 .or. itask==5 ) then
                           tcrit = rwork(1)
                           if ( (tcrit-tout)*(tout-t)<zero ) goto 1300
                           if ( h0/=zero .and. (t+h0-tcrit)*h0>zero ) &
                                h0 = tcrit - t
                        endif
                        me%dat%jstart = 0
                        if ( me%dat%miter>0 ) rwork(me%dat%lwm) = sqrt(me%dat%uround)
                        me%dat%ccmxj = pt2
                        me%dat%msbj = 50
                        me%dat%nhnil = 0
                        me%dat%nst = 0
                        me%dat%nje = 0
                        me%dat%nni = 0
                        me%dat%ncfn = 0
                        me%dat%netf = 0
                        me%dat%nlu = 0
                        me%dat%nslj = 0
                        nslast = 0
                        me%dat%hu = zero
                        me%dat%nqu = 0
                        ! initial call to f.  (lf0 points to yh(*,2).) -------------------------
                        lf0 = me%dat%lyh + me%dat%nyh
                        call me%f(me%dat%n,t,y(1:me%dat%n),rwork(lf0))
                        me%dat%nfe = 1
                        ! load the initial value vector in yh. ---------------------------------
                        call dcopy(me%dat%n,y,1,rwork(me%dat%lyh),1)
                        ! load and invert the ewt array.  (h is temporarily set to 1.0.) -------
                        me%dat%nq = 1
                        me%dat%h = one
                        call me%dewset(me%dat%n,itol,rtol(1:me%dat%n),atol(1:me%dat%n),&
                                       rwork(me%dat%lyh), rwork(me%dat%lewt))
                        do i = 1 , me%dat%n
                           if ( rwork(i+me%dat%lewt-1)<=zero ) goto 1100
                           rwork(i+me%dat%lewt-1) = one/rwork(i+me%dat%lewt-1)
                        enddo
                        if ( h0==zero ) then
                           ! call dvhin to set initial step size h0 to be attempted. --------------
                           call me%dvhin(me%dat%n,t,rwork(me%dat%lyh),rwork(lf0), &
                                         tout,me%dat%uround,rwork(me%dat%lewt),itol,&
                                         atol,y,rwork(me%dat%lacor),h0,niter,ier)
                           me%dat%nfe = me%dat%nfe + niter
                           if ( ier/=0 ) then
                              msg = &
                              'dvode--  tout (=r1) too close to t(=r2) to start integration'
                              call me%xerrwd(msg,60,22,1,0,0,0,2,tout,t)
                              goto 1500
                           endif
                        endif
                        ! adjust h0 if necessary to meet hmax bound. ---------------------------
                        rh = abs(h0)*me%dat%hmxi
                        if ( rh>one ) h0 = h0/rh
                        ! load h with h0 and scale yh(*,2) by h0. ------------------------------
                        me%dat%h = h0
                        call dscal(me%dat%n,h0,rwork(lf0),1)
                        goto 200
                     else
                        ! if istate = 3, set flag to signal parameter changes to dvstep. -------
                        me%dat%jstart = -1
                        ! maxord was reduced below nq.  copy yh(*,maxord+2) into savf. ---------
                        if ( me%dat%nq>me%dat%maxord ) &
                             call dcopy(me%dat%n,rwork(me%dat%lwm),1,rwork(me%dat%lsavf),1)
                        ! reload wm(1) = rwork(lwm), since lwm may have changed. ---------------
                        if ( me%dat%miter>0 ) rwork(me%dat%lwm) = sqrt(me%dat%uround)
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
 50      nslast = me%dat%nst
         me%dat%kuth = 0
         select case (itask)
         case (2)
            goto 100
         case (3)
            tp = me%dat%tn - me%dat%hu*(one+hun*me%dat%uround)
            if ( (tp-tout)*me%dat%h>zero ) then
               msg = 'dvode--  itask = i1 and tout (=r1) behind tcur - hu (= r2)  '
               call me%xerrwd(msg,60,23,1,1,itask,0,2,tout,tp)
               goto 1500
            else
               if ( (me%dat%tn-tout)*me%dat%h>=zero ) goto 300
               goto 100
            endif
         case (4)
            tcrit = rwork(1)
            if ( (me%dat%tn-tcrit)*me%dat%h>zero ) goto 1200
            if ( (tcrit-tout)*me%dat%h<zero ) goto 1300
            if ( (me%dat%tn-tout)*me%dat%h>=zero ) then
               call me%dvindy(tout,0,rwork(me%dat%lyh),me%dat%nyh,y,iflag)
               if ( iflag/=0 ) goto 1400
               t = tout
               goto 400
            endif
         case (5)
            tcrit = rwork(1)
            if ( (me%dat%tn-tcrit)*me%dat%h>zero ) goto 1200
         case default
            if ( (me%dat%tn-tout)*me%dat%h<zero ) goto 100
            call me%dvindy(tout,0,rwork(me%dat%lyh),me%dat%nyh,y,iflag)
            if ( iflag/=0 ) goto 1400
            t = tout
            goto 400
         end select
         hmx = abs(me%dat%tn) + abs(me%dat%h)
         ihit = abs(me%dat%tn-tcrit)<=hun*me%dat%uround*hmx
         if ( ihit ) goto 300
         tnext = me%dat%tn + me%dat%hnew*(one+four*me%dat%uround)
         if ( (tnext-tcrit)*me%dat%h>zero ) then
            me%dat%h = (tcrit-me%dat%tn)*(one-four*me%dat%uround)
            me%dat%kuth = 1
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
 100  if ( (me%dat%nst-nslast)>=me%dat%mxstep ) then
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
         call me%xerrwd(msg,50,201,1,0,0,0,0,zero,zero)
         msg = '      taken on this call before reaching tout     '
         call me%xerrwd(msg,50,201,1,1,me%dat%mxstep,0,1,me%dat%tn,zero)
         istate = -1
         goto 700
      else
         call me%dewset(me%dat%n,itol,rtol(1:me%dat%n),atol(1:me%dat%n),&
                        rwork(me%dat%lyh),rwork(me%dat%lewt))
         do i = 1 , me%dat%n
            if ( rwork(i+me%dat%lewt-1)<=zero ) goto 500
            rwork(i+me%dat%lewt-1) = one/rwork(i+me%dat%lewt-1)
         enddo
      endif
 200  tolsf = me%dat%uround*me%dvnorm(me%dat%n,rwork(me%dat%lyh),rwork(me%dat%lewt))
      if ( tolsf<=one ) then
         if ( (me%dat%tn+me%dat%h)==me%dat%tn ) then
            me%dat%nhnil = me%dat%nhnil + 1
            if ( me%dat%nhnil<=me%dat%mxhnil ) then
               msg = 'dvode--  warning: internal t (=r1) and h (=r2) are'
               call me%xerrwd(msg,50,101,1,0,0,0,0,zero,zero)
               msg = '      such that in the machine, t + h = t on the next step  '
               call me%xerrwd(msg,60,101,1,0,0,0,0,zero,zero)
               msg = '      (h = step size). solver will continue anyway'
               call me%xerrwd(msg,50,101,1,0,0,0,2,me%dat%tn,me%dat%h)
               if ( me%dat%nhnil>=me%dat%mxhnil ) then
                  msg = 'dvode--  above warning has been issued i1 times.  '
                  call me%xerrwd(msg,50,102,1,0,0,0,0,zero,zero)
                  msg = '      it will not be issued again for this problem'
                  call me%xerrwd(msg,50,102,1,1,me%dat%mxhnil,0,0,zero,zero)
               endif
            endif
         endif
         call me%dvstep(y,&                     ! y
                        rwork(me%dat%lyh),&     ! yh
                        me%dat%nyh,&            ! nyh
                        rwork(me%dat%lyh),&     ! yh
                        rwork(me%dat%lewt), &   ! ewt
                        rwork(me%dat%lsavf),&   ! savf
                        y,&                     ! vsav
                        rwork(me%dat%lacor),&   ! acor
                        rwork(me%dat%lwm),&     ! wm
                        iwork(me%dat%liwm))     ! iwm
         kgo = 1 - me%dat%kflag
         ! branch on kflag.  note: in this version, kflag can not be set to -3.
         !  kflag == 0,   -1,  -2
         select case (kgo)
         case (2)
            ! kflag = -1.  error test failed repeatedly or with abs(h) = hmin. -----
            msg = 'dvode--  at t(=r1) and step size h(=r2), the error'
            call me%xerrwd(msg,50,204,1,0,0,0,0,zero,zero)
            msg = '      test failed repeatedly or with abs(h) = hmin'
            call me%xerrwd(msg,50,204,1,0,0,0,2,me%dat%tn,me%dat%h)
            istate = -4
            goto 600
         case (3)
            ! kflag = -2.  convergence failed repeatedly or with abs(h) = hmin. ----
            msg = 'dvode--  at t (=r1) and step size h (=r2), the    '
            call me%xerrwd(msg,50,205,1,0,0,0,0,zero,zero)
            msg = '      corrector convergence failed repeatedly     '
            call me%xerrwd(msg,50,205,1,0,0,0,0,zero,zero)
            msg = '      or with abs(h) = hmin   '
            call me%xerrwd(msg,30,205,1,0,0,0,2,me%dat%tn,me%dat%h)
            istate = -5
            goto 600
         case default
            !-----------------------------------------------------------------------
            ! block f.
            ! the following block handles the case of a successful return from the
            ! core integrator (kflag = 0).  test for stop conditions.
            !-----------------------------------------------------------------------
            me%dat%init = 1
            me%dat%kuth = 0
            select case (itask)
            case (2)
            case (3)
               ! itask = 3.  jump to exit if tout was reached. ------------------------
               if ( (me%dat%tn-tout)*me%dat%h<zero ) goto 100
            case (4)
               ! itask = 4.  see if tout or tcrit was reached.  adjust h if necessary.
               if ( (me%dat%tn-tout)*me%dat%h<zero ) then
                  hmx = abs(me%dat%tn) + abs(me%dat%h)
                  ihit = abs(me%dat%tn-tcrit)<=hun*me%dat%uround*hmx
                  if ( .not.(ihit) ) then
                     tnext = me%dat%tn + me%dat%hnew*(one+four*me%dat%uround)
                     if ( (tnext-tcrit)*me%dat%h>zero ) then
                        me%dat%h = (tcrit-me%dat%tn)*(one-four*me%dat%uround)
                        me%dat%kuth = 1
                     endif
                     goto 100
                  endif
               else
                  call me%dvindy(tout,0,rwork(me%dat%lyh),me%dat%nyh,y,iflag)
                  t = tout
                  goto 400
               endif
            case (5)
               ! itask = 5.  see if tcrit was reached and jump to exit. ---------------
               hmx = abs(me%dat%tn) + abs(me%dat%h)
               ihit = abs(me%dat%tn-tcrit)<=hun*me%dat%uround*hmx
            case default
               ! itask = 1.  if tout has been reached, interpolate. -------------------
               if ( (me%dat%tn-tout)*me%dat%h<zero ) goto 100
               call me%dvindy(tout,0,rwork(me%dat%lyh),me%dat%nyh,y,iflag)
               t = tout
               goto 400
            end select
         end select
      else
         tolsf = tolsf*two
         if ( me%dat%nst==0 ) then
            msg = 'dvode--  at start of problem, too much accuracy   '
            call me%xerrwd(msg,50,26,1,0,0,0,0,zero,zero)
            msg = '      requested for precision of machine:   see tolsf (=r1) '
            call me%xerrwd(msg,60,26,1,0,0,0,1,tolsf,zero)
            rwork(14) = tolsf
            goto 1500
         else
            ! too much accuracy requested for machine precision. -------------------
            msg = 'dvode--  at t (=r1), too much accuracy requested  '
            call me%xerrwd(msg,50,203,1,0,0,0,0,zero,zero)
            msg = '      for precision of machine:   see tolsf (=r2) '
            call me%xerrwd(msg,50,203,1,0,0,0,2,me%dat%tn,tolsf)
            rwork(14) = tolsf
            istate = -2
            goto 700
         endif
      endif
      !-----------------------------------------------------------------------
      ! block g.
      ! the following block handles all successful returns from dvode.
      ! if itask /= 1, y is loaded from yh and t is set accordingly.
      ! istate is set to 2, and the optional output is loaded into the work
      ! arrays before returning.
      !-----------------------------------------------------------------------
 300  call dcopy(me%dat%n,rwork(me%dat%lyh),1,y,1)
      t = me%dat%tn
      if ( itask==4 .or. itask==5 ) then
         if ( ihit ) t = tcrit
      endif
 400  istate = 2
      rwork(11) = me%dat%hu
      rwork(12) = me%dat%hnew
      rwork(13) = me%dat%tn
      iwork(11) = me%dat%nst
      iwork(12) = me%dat%nfe
      iwork(13) = me%dat%nje
      iwork(14) = me%dat%nqu
      iwork(15) = me%dat%newq
      iwork(19) = me%dat%nlu
      iwork(20) = me%dat%nni
      iwork(21) = me%dat%ncfn
      iwork(22) = me%dat%netf
      return

      ! ewt(i) <= 0.0 for some i (not at start of problem). ----------------
 500  ewti = rwork(me%dat%lewt+i-1)
      msg = 'dvode--  at t (=r1), ewt(i1) has become r2 <= 0.'
      call me%xerrwd(msg,50,202,1,1,i,0,2,me%dat%tn,ewti)
      istate = -6
      goto 700

      ! compute imxer if relevant. -------------------------------------------
 600  big = zero
      imxer = 1
      do i = 1 , me%dat%n
         size = abs(rwork(i+me%dat%lacor-1)*rwork(i+me%dat%lewt-1))
         if ( big<size ) then
            big = size
            imxer = i
         endif
      enddo
      iwork(16) = imxer
      ! set y vector, t, and optional output. --------------------------------
 700  call dcopy(me%dat%n,rwork(me%dat%lyh),1,y,1)
      t = me%dat%tn
      rwork(11) = me%dat%hu
      rwork(12) = me%dat%h
      rwork(13) = me%dat%tn
      iwork(11) = me%dat%nst
      iwork(12) = me%dat%nfe
      iwork(13) = me%dat%nje
      iwork(14) = me%dat%nqu
      iwork(15) = me%dat%nq
      iwork(19) = me%dat%nlu
      iwork(20) = me%dat%nni
      iwork(21) = me%dat%ncfn
      iwork(22) = me%dat%netf

      return

 800  msg = 'dvode--  mf (=i1) illegal     '
      call me%xerrwd(msg,30,8,1,1,mf,0,0,zero,zero)
      goto 1500
 900  msg = 'dvode--  rtol(i1) is r1 < 0.0        '
      call me%xerrwd(msg,40,19,1,1,i,0,1,rtoli,zero)
      goto 1500
 1000 msg = 'dvode--  atol(i1) is r1 < 0.0        '
      call me%xerrwd(msg,40,20,1,1,i,0,1,atoli,zero)
      goto 1500
 1100 ewti = rwork(me%dat%lewt+i-1)
      msg = 'dvode--  ewt(i1) is r1 <= 0.0         '
      call me%xerrwd(msg,40,21,1,1,i,0,1,ewti,zero)
      goto 1500
 1200 msg = 'dvode--  itask = 4 or 5 and tcrit (=r1) behind tcur (=r2)   '
      call me%xerrwd(msg,60,24,1,0,0,0,2,tcrit,me%dat%tn)
      goto 1500
 1300 msg = 'dvode--  itask = 4 or 5 and tcrit (=r1) behind tout (=r2)   '
      call me%xerrwd(msg,60,25,1,0,0,0,2,tcrit,tout)
      goto 1500
 1400 msg = 'dvode--  trouble from dvindy.  itask = i1, tout = r1.       '
      call me%xerrwd(msg,60,27,1,1,itask,0,1,tout,zero)

 1500 istate = -3

   end subroutine dvode

!*****************************************************************************************
!>
!  this routine computes the step size, h0, to be attempted on the
!  first step, when the user has not supplied a value for this.
!
!  first we check that tout - t0 differs significantly from zero.  then
!  an iteration is done to approximate the initial second derivative
!  and this is used to define h from w.r.m.s.norm(h**2 * yddot / 2) = 1.
!  a bias factor of 1/2 is applied to the resulting h.
!  the sign of h0 is inferred from the initial values of tout and t0.

   subroutine dvhin(me,n,t0,y0,ydot,tout,uround,ewt,itol,&
                    atol,y,temp,h0,niter,ier)

      class(dvode_t),intent(inout) :: me
      real(wp),intent(in) :: t0 !! initial value of independent variable
      real(wp),intent(in) :: y0(*) !! vector of initial conditions
      real(wp),intent(in) :: ydot(*) !! vector of initial first derivatives
      real(wp),intent(in) :: tout !! first output value of independent variable
      real(wp),intent(in) :: uround !! machine unit roundoff
      real(wp),intent(in) :: ewt(*) !! error weights and tolerance parameters as described in the driver routine
      real(wp),intent(in) :: atol(*) !! error weights and tolerance parameters as described in the driver routine
      real(wp),intent(inout) :: y(n) !! work array of length n
      real(wp),intent(inout) :: temp(n) !! work array of length n
      integer,intent(in) :: n !! size of ode system
      integer,intent(in) :: itol !! error weights and tolerance parameters as described in the driver routine
      real(wp),intent(out) :: h0 !! step size to be attempted
      integer,intent(out) :: niter !! number of iterations (and of f evaluations) to compute h0
      integer,intent(out) :: ier !! the error flag, returned with the value:
                                 !!
                                 !!  * ier = 0  if no trouble occurred, or
                                 !!  * ier = -1 if tout and t0 are considered too close to proceed.

      real(wp) :: afi , atoli , delyi , h , hg , hlb , &
                  hnew , hrat , hub , t1 , tdist , &
                  tround , yddnrm
      integer :: i , iter

      real(wp),parameter :: half = 0.5_wp
      real(wp),parameter :: hun = 100.0_wp
      real(wp),parameter :: pt1 = 0.1_wp
      real(wp),parameter :: two = 2.0_wp

      main : block

         niter = 0
         tdist = abs(tout-t0)
         tround = uround*max(abs(t0),abs(tout))
         if ( tdist<two*tround ) then
            ! error return for tout - t0 too small. --------------------------------
            ier = -1
            return
         else
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
            ! set initial guess for h as geometric mean of upper and lower bounds. -
            iter = 0
            hg = sqrt(hlb*hub)
            ! if the bounds have crossed, exit with the mean value. ----------------
            if ( hub<hlb ) then
               h0 = hg
               exit main
            endif
            do
               ! looping point for iteration. -----------------------------------------
               ! estimate the second derivative as a difference quotient in f. --------
               h = sign(hg,tout-t0)
               t1 = t0 + h
               do i = 1 , n
                  y(i) = y0(i) + h*ydot(i)
               enddo
               call me%f(n,t1,y,temp)
               do i = 1 , n
                  temp(i) = (temp(i)-ydot(i))/h
               enddo
               yddnrm = me%dvnorm(n,temp,ewt(1:n))
               ! get the corresponding new value of h. --------------------------------
               if ( yddnrm*hub*hub>two ) then
                  hnew = sqrt(two/yddnrm)
               else
                  hnew = sqrt(hg*hub)
               endif
               iter = iter + 1
               !-----------------------------------------------------------------------
               ! test the stopping conditions.
               ! stop if the new and previous h values differ by a factor of < 2.
               ! stop if four iterations have been done.  also, stop with previous h
               ! if hnew/hg > 2 after first iteration, as this probably means that
               ! the second derivative value is bad because of cancellation error.
               !-----------------------------------------------------------------------
               if ( iter<4 ) then
                  hrat = hnew/hg
                  if ( (hrat<=half) .or. (hrat>=two) ) then
                     if ( (iter>=2) .and. (hnew>two*hg) ) then
                        hnew = hg
                        exit
                     endif
                     hg = hnew
                     cycle
                  endif
               endif
               exit
            end do
            ! iteration done.  apply bounds, bias factor, and sign.  then exit. ----
            h0 = hnew*half
            if ( h0<hlb ) h0 = hlb
            if ( h0>hub ) h0 = hub
         endif

      end block main

      h0 = sign(h0,tout-t0)
      niter = iter
      ier = 0

   end subroutine dvhin

!*****************************************************************************************
!>
!  dvindy computes interpolated values of the k-th derivative of the
!  dependent variable vector y, and stores it in dky.  this routine
!  is called within the package with k = 0 and t = tout, but may
!  also be called by the user for any k up to the current order.
!  (see detailed instructions in the usage documentation.)
!
!```
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
!```

      subroutine dvindy(me,t,k,yh,ldyh,dky,iflag)

      class(dvode_t),intent(inout) :: me
      real(wp),intent(in) :: t !! value of independent variable where answers are desired
                               !! (normally the same as the t last returned by [[dvode]]).
                               !! for valid results, t must lie between tcur - hu and tcur.
                               !! (see optional output for tcur and hu.)
      integer,intent(in) :: k !! integer order of the derivative desired.  k must satisfy
                              !! 0 <= k <= nqcur, where nqcur is the current order
                              !! (see optional output).  the capability corresponding
                              !! to k = 0, i.e. computing y(t), is already provided
                              !! by [[dvode]] directly.  since nqcur >= 1, the first
                              !! derivative dy/dt is always available with [[dvindy]].
      integer,intent(in) :: ldyh !! column length of yh, equal to the initial value of neq.
      real(wp),intent(in) :: yh(ldyh,*) !! the history array yh
      real(wp),intent(out) :: dky(*) !! a real array of length neq containing the computed value
                                     !! of the k-th derivative of y(t).
      integer,intent(out) :: iflag !! integer flag, returned as 0 if k and t were legal,
                                   !! -1 if k was illegal, and -2 if t was illegal.
                                   !! on an error return, a message is also written.

      real(wp) :: c , r , s , tfuzz , tn1 , tp
      integer :: i , ic , j , jb , jb2 , jj , jj1 , jp1
      character(len=80) :: msg

      real(wp),parameter :: hun = 100.0_wp

      iflag = 0
      if ( k<0 .or. k>me%dat%nq ) then
         msg = 'dvindy-- k (=i1) illegal      '
         call me%xerrwd(msg,30,51,1,1,k,0,0,zero,zero)
         iflag = -1
         return
      else
         tfuzz = hun*me%dat%uround*sign(abs(me%dat%tn)+abs(me%dat%hu),me%dat%hu)
         tp = me%dat%tn - me%dat%hu - tfuzz
         tn1 = me%dat%tn + tfuzz
         if ( (t-tp)*(t-tn1)>zero ) then
            msg = 'dvindy-- t (=r1) illegal      '
            call me%xerrwd(msg,30,52,1,0,0,0,1,t,zero)
            msg =  '      t not in interval tcur - hu (= r1) to tcur (=r2)      '
            call me%xerrwd(msg,60,52,1,0,0,0,2,tp,me%dat%tn)
            iflag = -2
            return
         else
            s = (t-me%dat%tn)/me%dat%h
            ic = 1
            if ( k/=0 ) then
               jj1 = me%dat%l - k
               do jj = jj1 , me%dat%nq
                  ic = ic*jj
               enddo
            endif
            c = real(ic,wp)
            do i = 1 , me%dat%n
               dky(i) = c*yh(i,me%dat%l)
            enddo
            if ( k/=me%dat%nq ) then
               jb2 = me%dat%nq - k
               do jb = 1 , jb2
                  j = me%dat%nq - jb
                  jp1 = j + 1
                  ic = 1
                  if ( k/=0 ) then
                     jj1 = jp1 - k
                     do jj = jj1 , j
                        ic = ic*jj
                     enddo
                  endif
                  c = real(ic,wp)
                  do i = 1 , me%dat%n
                     dky(i) = c*yh(i,jp1) + s*dky(i)
                  enddo
               enddo
               if ( k==0 ) return
            endif
         endif
      endif
      r = me%dat%h**(-k)
      call dscal(me%dat%n,r,dky,1)

   end subroutine dvindy

!*****************************************************************************************
!>
!  dvstep performs one step of the integration of an initial value
!  problem for a system of ordinary differential equations.
!
!  dvstep calls subroutine [[dvnlsd]] for the solution of the nonlinear system
!  arising in the time step.  thus it is independent of the problem
!  jacobian structure and the type of nonlinear system solution method.
!  dvstep returns a completion flag kflag (in common).
!  a return with kflag = -1 or -2 means either abs(h) = hmin or 10
!  consecutive failures occurred.  on a return with kflag negative,
!  the values of tn and the yh array are as of the beginning of the last
!  step, and h is the last step size attempted.

   subroutine dvstep(me,y,yh,ldyh,yh1,ewt,savf,vsav,acor,wm,iwm)

      class(dvode_t),intent(inout) :: me
      integer,intent(in) :: ldyh !! a constant integer >= n, the first dimension of yh.
                                 !! n is the number of odes in the system.
      real(wp),intent(inout) :: y(*) !! an array of length n used for the dependent variable vector.
      real(wp),intent(inout) :: yh(ldyh,*) !! an ldyh by lmax array containing the dependent variables
                                           !! and their approximate scaled derivatives, where
                                           !! lmax = maxord + 1.  yh(i,j+1) contains the approximate
                                           !! j-th derivative of y(i), scaled by h**j/factorial(j)
                                           !! (j = 0,1,...,nq).  on entry for the first step, the first
                                           !! two columns of yh must be set from the initial values.
      real(wp),intent(inout) :: yh1(*) !! a one-dimensional array occupying the same space as yh.
      real(wp),intent(in) :: ewt(*) !! an array of length n containing multiplicative weights
                                    !! for local error measurements.  local errors in y(i) are
                                    !! compared to 1.0/ewt(i) in various error tests.
      real(wp) :: savf(*) !! an array of working storage, of length n.
                          !! also used for input of yh(*,maxord+2) when jstart = -1
                          !! and maxord < the current order nq.
      real(wp) :: vsav(*) !! a work array of length n passed to subroutine [[dvnlsd]].
      real(wp) :: acor(*) !! a work array of length n, used for the accumulated
                          !! corrections.  on a successful return, acor(i) contains
                          !! the estimated one-step local error in y(i).
      real(wp) :: wm(*) !! real work array associated with matrix
                        !! operations in [[dvnlsd]].
      integer :: iwm(*) !! integer work array associated with matrix
                        !! operations in [[dvnlsd]].

      real(wp) :: cnquot , ddn , dsm , dup , etaqp1 , &
                  flotl , r , told
      integer :: i , i1 , i2 , iback , j , jb , ncf , nflag

      integer,parameter :: kfc = -3
      integer,parameter :: kfh = -7
      integer,parameter :: mxncf = 10
      real(wp),parameter :: one = 1.0_wp
      real(wp),parameter :: addon = 1.0e-6_wp
      real(wp),parameter :: bias1 = 6.0_wp
      real(wp),parameter :: bias2 = 6.0_wp
      real(wp),parameter :: bias3 = 10.0_wp
      real(wp),parameter :: etacf = 0.25_wp
      real(wp),parameter :: etamin = 0.1_wp
      real(wp),parameter :: etamxf = 0.2_wp
      real(wp),parameter :: etamx1 = 10000.0_wp
      real(wp),parameter :: etamx2 = 10.0_wp
      real(wp),parameter :: etamx3 = 10.0_wp
      real(wp),parameter :: onepsm = 1.00001_wp
      real(wp),parameter :: thresh = 1.5_wp

      me%dat%kflag = 0
      told = me%dat%tn
      ncf = 0
      me%dat%jcur = 0
      nflag = 0

      if ( me%dat%jstart<=0 ) then
         if ( me%dat%jstart==-1 ) then
            !-----------------------------------------------------------------------
            ! the following block handles preliminaries needed when jstart = -1.
            ! if n was reduced, zero out part of yh to avoid undefined references.
            ! if maxord was reduced to a value less than the tentative order newq,
            ! then nq is set to maxord, and a new h ratio eta is chosen.
            ! otherwise, we take the same preliminary actions as for jstart > 0.
            ! in any case, nqwait is reset to l = nq + 1 to prevent further
            ! changes in order for that many steps.
            ! the new h ratio eta is limited by the input h if kuth = 1,
            ! by hmin if kuth = 0, and by hmxi in any case.
            ! finally, the history array yh is rescaled.
            !-----------------------------------------------------------------------
            me%dat%lmax = me%dat%maxord + 1
            if ( me%dat%n/=ldyh ) then
               i1 = 1 + (me%dat%newq+1)*ldyh
               i2 = (me%dat%maxord+1)*ldyh
               if ( i1<=i2 ) then
                  do i = i1 , i2
                     yh1(i) = zero
                  enddo
               endif
            endif
            if ( me%dat%newq>me%dat%maxord ) then
               flotl = real(me%dat%lmax,wp)
               if ( me%dat%maxord<me%dat%nq-1 ) then
                  ddn = me%dvnorm(me%dat%n,savf(1:me%dat%n),ewt(1:me%dat%n))/me%dat%tq(1)
                  me%dat%eta = one/((bias1*ddn)**(one/flotl)+addon)
               endif
               if ( me%dat%maxord==me%dat%nq .and. me%dat%newq==me%dat%nq+1 ) me%dat%eta = me%etaq
               if ( me%dat%maxord==me%dat%nq-1 .and. me%dat%newq==me%dat%nq+1 ) then
                  me%dat%eta = me%etaqm1
                  call me%dvjust(yh,ldyh,-1)
               endif
               if ( me%dat%maxord==me%dat%nq-1 .and. me%dat%newq==me%dat%nq ) then
                  ddn = me%dvnorm(me%dat%n,savf(1:me%dat%n),ewt(1:me%dat%n))/me%dat%tq(1)
                  me%dat%eta = one/((bias1*ddn)**(one/flotl)+addon)
                  call me%dvjust(yh,ldyh,-1)
               endif
               me%dat%eta = min(me%dat%eta,one)
               me%dat%nq = me%dat%maxord
               me%dat%l = me%dat%lmax
            endif
            if ( me%dat%kuth==1 ) me%dat%eta = min(me%dat%eta,abs(me%dat%h/me%dat%hscal))
            if ( me%dat%kuth==0 ) me%dat%eta = max(me%dat%eta,me%dat%hmin/abs(me%dat%hscal))
            me%dat%eta = me%dat%eta/max(one,abs(me%dat%hscal)*me%dat%hmxi*me%dat%eta)
            me%dat%newh = 1
            me%dat%nqwait = me%dat%l
            if ( me%dat%newq>me%dat%maxord ) goto 300
         else
            !-----------------------------------------------------------------------
            ! on the first call, the order is set to 1, and other variables are
            ! initialized.  etamax is the maximum ratio by which h can be increased
            ! in a single step.  it is normally 10, but is larger during the
            ! first step to compensate for the small initial h.  if a failure
            ! occurs (in corrector convergence or error test), etamax is set to 1
            ! for the next increase.
            !-----------------------------------------------------------------------
            me%dat%lmax = me%dat%maxord + 1
            me%dat%nq = 1
            me%dat%l = 2
            me%dat%nqnyh = me%dat%nq*ldyh
            me%dat%tau(1) = me%dat%h
            me%dat%prl1 = one
            me%dat%rc = zero
            me%dat%etamax = etamx1
            me%dat%nqwait = 2
            me%dat%hscal = me%dat%h
            goto 400
         end if

      else if ( me%dat%kuth==1 ) then
         !-----------------------------------------------------------------------
         ! take preliminary actions on a normal continuation step (jstart>0).
         ! if the driver changed h, then eta must be reset and newh set to 1.
         ! if a change of order was dictated on the previous step, then
         ! it is done here and appropriate adjustments in the history are made.
         ! on an order decrease, the history array is adjusted by dvjust.
         ! on an order increase, the history array is augmented by a column.
         ! on a change of step size h, the history array yh is rescaled.
         !-----------------------------------------------------------------------
         me%dat%eta = min(me%dat%eta,me%dat%h/me%dat%hscal)
         me%dat%newh = 1
      endif

      if ( me%dat%newh==0 ) goto 400

      if ( me%dat%newq<me%dat%nq ) then
         call me%dvjust(yh,ldyh,-1)
         me%dat%nq = me%dat%newq
         me%dat%l = me%dat%nq + 1
         me%dat%nqwait = me%dat%l
      else if ( me%dat%newq>me%dat%nq ) then
         call me%dvjust(yh,ldyh,1)
         me%dat%nq = me%dat%newq
         me%dat%l = me%dat%nq + 1
         me%dat%nqwait = me%dat%l
      endif

      ! rescale the history array for a change in h by a factor of eta. ------
 300  r = one
      do j = 2 , me%dat%l
         r = r*me%dat%eta
         call dscal(me%dat%n,r,yh(1,j),1)
      enddo
      me%dat%h = me%dat%hscal*me%dat%eta
      me%dat%hscal = me%dat%h
      me%dat%rc = me%dat%rc*me%dat%eta
      me%dat%nqnyh = me%dat%nq*ldyh

      !-----------------------------------------------------------------------
      ! this section computes the predicted values by effectively
      ! multiplying the yh array by the pascal triangle matrix.
      ! dvset is called to calculate all integration coefficients.
      ! rc is the ratio of new to old values of the coefficient h/el(2)=h/l1.
      !-----------------------------------------------------------------------
 400  me%dat%tn = me%dat%tn + me%dat%h
      i1 = me%dat%nqnyh + 1
      do jb = 1 , me%dat%nq
         i1 = i1 - ldyh
         do i = i1 , me%dat%nqnyh
            yh1(i) = yh1(i) + yh1(i+ldyh)
         enddo
      enddo
      call me%dvset()
      me%dat%rl1 = one/me%dat%el(2)
      me%dat%rc = me%dat%rc*(me%dat%rl1/me%dat%prl1)
      me%dat%prl1 = me%dat%rl1

      ! call the nonlinear system solver. ------------------------------------
      call me%dvnlsd(y,yh,ldyh,vsav,savf,ewt,acor,iwm,wm,nflag)

      if ( nflag==0 ) then
         !-----------------------------------------------------------------------
         ! the corrector has converged (nflag = 0).  the local error test is
         ! made and control passes to statement 500 if it fails.
         !-----------------------------------------------------------------------
         dsm = me%dat%acnrm/me%dat%tq(2)
         if ( dsm>one ) then
            !-----------------------------------------------------------------------
            ! the error test failed.  kflag keeps track of multiple failures.
            ! restore tn and the yh array to their previous values, and prepare
            ! to try the step again.  compute the optimum step size for the
            ! same order.  after repeated failures, h is forced to decrease
            ! more rapidly.
            !-----------------------------------------------------------------------
            me%dat%kflag = me%dat%kflag - 1
            me%dat%netf = me%dat%netf + 1
            nflag = -2
            me%dat%tn = told
            i1 = me%dat%nqnyh + 1
            do jb = 1 , me%dat%nq
               i1 = i1 - ldyh
               do i = i1 , me%dat%nqnyh
                  yh1(i) = yh1(i) - yh1(i+ldyh)
               enddo
            enddo
            if ( abs(me%dat%h)<=me%dat%hmin*onepsm ) then
               !-----------------------------------------------------------------------
               ! all returns are made through this section.
               ! on a successful return, etamax is reset and acor is scaled.
               !-----------------------------------------------------------------------
               me%dat%kflag = -1
               goto 600
            else
               me%dat%etamax = one
               if ( me%dat%kflag>kfc ) then
                  ! compute ratio of new h to current h at the current order. ------------
                  flotl = real(me%dat%l,wp)
                  me%dat%eta = one/((bias2*dsm)**(one/flotl)+addon)
                  me%dat%eta = max(me%dat%eta,me%dat%hmin/abs(me%dat%h),etamin)
                  if ( (me%dat%kflag<=-2) .and. (me%dat%eta>etamxf) ) me%dat%eta = etamxf
                  goto 300
                  !-----------------------------------------------------------------------
                  ! control reaches this section if 3 or more consecutive failures
                  ! have occurred.  it is assumed that the elements of the yh array
                  ! have accumulated errors of the wrong order.  the order is reduced
                  ! by one, if possible.  then h is reduced by a factor of 0.1 and
                  ! the step is retried.  after a total of 7 consecutive failures,
                  ! an exit is taken with kflag = -1.
                  !-----------------------------------------------------------------------
               else if ( me%dat%kflag==kfh ) then
                  me%dat%kflag = -1
                  goto 600
               else if ( me%dat%nq==1 ) then
                  me%dat%eta = max(etamin,me%dat%hmin/abs(me%dat%h))
                  me%dat%h = me%dat%h*me%dat%eta
                  me%dat%hscal = me%dat%h
                  me%dat%tau(1) = me%dat%h
                  call me%f(me%dat%n,me%dat%tn,y(1:me%dat%n),savf(1:me%dat%n))
                  me%dat%nfe = me%dat%nfe + 1
                  do i = 1 , me%dat%n
                     yh(i,2) = me%dat%h*savf(i)
                  enddo
                  me%dat%nqwait = 10
                  goto 400
               else
                  me%dat%eta = max(etamin,me%dat%hmin/abs(me%dat%h))
                  call me%dvjust(yh,ldyh,-1)
                  me%dat%l = me%dat%nq
                  me%dat%nq = me%dat%nq - 1
                  me%dat%nqwait = me%dat%l
                  goto 300
               endif
            endif
         else
            !-----------------------------------------------------------------------
            ! after a successful step, update the yh and tau arrays and decrement
            ! nqwait.  if nqwait is then 1 and nq < maxord, then acor is saved
            ! for use in a possible order increase on the next step.
            ! if etamax = 1 (a failure occurred this step), keep nqwait >= 2.
            !-----------------------------------------------------------------------
            me%dat%kflag = 0
            me%dat%nst = me%dat%nst + 1
            me%dat%hu = me%dat%h
            me%dat%nqu = me%dat%nq
            do iback = 1 , me%dat%nq
               i = me%dat%l - iback
               me%dat%tau(i+1) = me%dat%tau(i)
            enddo
            me%dat%tau(1) = me%dat%h
            do j = 1 , me%dat%l
               call daxpy(me%dat%n,me%dat%el(j),acor,1,yh(1,j),1)
            enddo
            me%dat%nqwait = me%dat%nqwait - 1
            if ( (me%dat%l/=me%dat%lmax) .and. (me%dat%nqwait==1) ) then
               call dcopy(me%dat%n,acor,1,yh(1,me%dat%lmax),1)
               me%dat%conp = me%dat%tq(5)
            endif
            if ( me%dat%etamax/=one ) then
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
               flotl = real(me%dat%l,wp)
               me%etaq = one/((bias2*dsm)**(one/flotl)+addon)
               if ( me%dat%nqwait==0 ) then
                  me%dat%nqwait = 2
                  me%etaqm1 = zero
                  if ( me%dat%nq/=1 ) then
                     ! compute ratio of new h to current h at the current order less one. ---
                     ddn = me%dvnorm(me%dat%n,yh(1,me%dat%l),ewt(1:me%dat%n))/me%dat%tq(1)
                     me%etaqm1 = one/((bias1*ddn)**(one/(flotl-one))+addon)
                  endif
                  etaqp1 = zero
                  if ( me%dat%l/=me%dat%lmax ) then
                     ! compute ratio of new h to current h at current order plus one. -------
                     cnquot = (me%dat%tq(5)/me%dat%conp)*(me%dat%h/me%dat%tau(2))**me%dat%l
                     do i = 1 , me%dat%n
                        savf(i) = acor(i) - cnquot*yh(i,me%dat%lmax)
                     enddo
                     dup = me%dvnorm(me%dat%n,savf(1:me%dat%n),ewt(1:me%dat%n))/me%dat%tq(3)
                     etaqp1 = one/((bias3*dup)**(one/(flotl+one))+addon)
                  endif
                  if ( me%etaq<etaqp1 ) then
                     if ( etaqp1<=me%etaqm1 ) goto 420
                     me%dat%eta = etaqp1
                     me%dat%newq = me%dat%nq + 1
                     call dcopy(me%dat%n,acor,1,yh(1,me%dat%lmax),1)
                     goto 450
                  else if ( me%etaq<me%etaqm1 ) then
                     goto 420
                  endif
               endif
               me%dat%eta = me%etaq
               me%dat%newq = me%dat%nq
               goto 450
            else
               if ( me%dat%nqwait<2 ) me%dat%nqwait = 2
               me%dat%newq = me%dat%nq
               me%dat%newh = 0
               me%dat%eta = one
               me%dat%hnew = me%dat%h
               goto 500
            endif
 420        me%dat%eta = me%etaqm1
            me%dat%newq = me%dat%nq - 1
         endif

         ! test tentative new h against thresh, etamax, and hmxi, then exit. ----
 450     if ( me%dat%eta<thresh .or. me%dat%etamax==one ) then
            me%dat%newq = me%dat%nq
            me%dat%newh = 0
            me%dat%eta = one
            me%dat%hnew = me%dat%h
         else
            me%dat%eta = min(me%dat%eta,me%dat%etamax)
            me%dat%eta = me%dat%eta/max(one,abs(me%dat%h)*me%dat%hmxi*me%dat%eta)
            me%dat%newh = 1
            me%dat%hnew = me%dat%h*me%dat%eta
         endif
      else
         !-----------------------------------------------------------------------
         ! the vnls routine failed to achieve convergence (nflag /= 0).
         ! the yh array is retracted to its values before prediction.
         ! the step size h is reduced and the step is retried, if possible.
         ! otherwise, an error exit is taken.
         !-----------------------------------------------------------------------
         ncf = ncf + 1
         me%dat%ncfn = me%dat%ncfn + 1
         me%dat%etamax = one
         me%dat%tn = told
         i1 = me%dat%nqnyh + 1
         do jb = 1 , me%dat%nq
            i1 = i1 - ldyh
            do i = i1 , me%dat%nqnyh
               yh1(i) = yh1(i) - yh1(i+ldyh)
            enddo
         enddo
         if ( nflag<-1 ) then
            if ( nflag==-2 ) me%dat%kflag = -3
            if ( nflag==-3 ) me%dat%kflag = -4
            goto 600
         else if ( abs(me%dat%h)<=me%dat%hmin*onepsm ) then
            me%dat%kflag = -2
            goto 600
         else if ( ncf==mxncf ) then
            me%dat%kflag = -2
            goto 600
         else
            me%dat%eta = etacf
            me%dat%eta = max(me%dat%eta,me%dat%hmin/abs(me%dat%h))
            nflag = -1
            goto 300
         endif
      endif

 500  me%dat%etamax = etamx3
      if ( me%dat%nst<=10 ) me%dat%etamax = etamx2
      r = one/me%dat%tq(2)
      call dscal(me%dat%n,r,acor,1)

 600  me%dat%jstart = 1

   end subroutine dvstep

!*****************************************************************************************
!>
!  dvset is called by [[dvstep]] and sets coefficients for use there.
!
!  for each order `nq`, the coefficients in `el` are calculated by use of
!  the generating polynomial `lambda(x)`, with coefficients `el(i)`.
!
!  $$ \lambda(x) = \mathrm{el}(1) + \mathrm{el}(2) x + \cdots + \mathrm{el}(n_q+1)(x^n_q) $$
!
!  for the backward differentiation formulas,
!
!  $$ \lambda(x) = \left(1 + x/x_i \cdot (n_q) \right) \prod_{i=1}^{n_q-1} \left(1 + x/x_i(i) \right) $$
!
! for the adams formulas,
!
!  $$ \frac{d}{dx} \lambda(x) = c  \prod_{i=1}^{n_q-1} \left(1 + x/x_i(i) \right)  $$
!  $$ \lambda(-1) = 0 , \lambda(0) = 1 $$
!
!  where `c` is a normalization constant.
!  in both cases, \(x_i(i)\) is defined by
!
!  $$ \begin{align}
!    h x_i(i) &= t_n  -  t_{n-i} \\
!             &= h + \tau(1) + \tau(2) + \cdots + \tau(i-1)
!     \end{align} $$
!
!  in addition to variables described previously, communication
!  with [[dvset]] uses the following class variables:
!
!   * `tau`    = a vector of length 13 containing the past `nq` values
!                of h.
!   * `el`     = a vector of length 13 in which vset stores the
!                coefficients for the corrector formula.
!   * `tq`     = a vector of length 5 in which vset stores constants
!                used for the convergence test, the error test, and the
!                selection of `h` at a new order.
!   * `meth`   = the basic method indicator.
!   * `nq`     = the current order.
!   * `l`      = `nq + 1`, the length of the vector stored in el, and
!                the number of columns of the `yh` array being used.
!   * `nqwait` = a counter controlling the frequency of order changes.
!                an order change is about to be considered if `nqwait = 1`.

   subroutine dvset(me)

      class(dvode_t),intent(inout) :: me

      real(wp) :: ahatn0 , alph0 , cnqm1 , csum , elp , &
                  em(13) , em0 , floti , flotl , flotnq , hsum , &
                  rxi , rxis , s , t1 , t2 , t3 , t4 , t5 ,  &
                  t6 , xi
      integer :: i , iback , j , jp1 , nqm1 , nqm2

      real(wp),parameter :: cortes = 0.1_wp
      real(wp),parameter :: one = 1.0_wp
      real(wp),parameter :: six = 6.0_wp
      real(wp),parameter :: two = 2.0_wp

      flotl = real(me%dat%l,wp)
      nqm1 = me%dat%nq - 1
      nqm2 = me%dat%nq - 2
      if ( me%dat%meth==2 ) then
         ! set coefficients for bdf methods. ------------------------------------
         do i = 3 , me%dat%l
            me%dat%el(i) = zero
         enddo
         me%dat%el(1) = one
         me%dat%el(2) = one
         alph0 = -one
         ahatn0 = -one
         hsum = me%dat%h
         rxi = one
         rxis = one
         if ( me%dat%nq/=1 ) then
            do j = 1 , nqm2
               ! in el, construct coefficients of (1+x/xi(1))*...*(1+x/xi(j+1)). ------
               hsum = hsum + me%dat%tau(j)
               rxi = me%dat%h/hsum
               jp1 = j + 1
               alph0 = alph0 - one/real(jp1,wp)
               do iback = 1 , jp1
                  i = (j+3) - iback
                  me%dat%el(i) = me%dat%el(i) + me%dat%el(i-1)*rxi
               enddo
            enddo
            alph0 = alph0 - one/real(me%dat%nq,wp)
            rxis = -me%dat%el(2) - alph0
            hsum = hsum + me%dat%tau(nqm1)
            rxi = me%dat%h/hsum
            ahatn0 = -me%dat%el(2) - rxi
            do iback = 1 , me%dat%nq
               i = (me%dat%nq+2) - iback
               me%dat%el(i) = me%dat%el(i) + me%dat%el(i-1)*rxis
            enddo
         endif
         t1 = one - ahatn0 + alph0
         t2 = one + real(me%dat%nq,wp)*t1
         me%dat%tq(2) = abs(alph0*t2/t1)
         me%dat%tq(5) = abs(t2/(me%dat%el(me%dat%l)*rxi/rxis))
         if ( me%dat%nqwait==1 ) then
            cnqm1 = rxis/me%dat%el(me%dat%l)
            t3 = alph0 + one/real(me%dat%nq,wp)
            t4 = ahatn0 + rxi
            elp = t3/(one-t4+t3)
            me%dat%tq(1) = abs(elp/cnqm1)
            hsum = hsum + me%dat%tau(me%dat%nq)
            rxi = me%dat%h/hsum
            t5 = alph0 - one/real(me%dat%nq+1,wp)
            t6 = ahatn0 - rxi
            elp = t2/(one-t6+t5)
            me%dat%tq(3) = abs(elp*rxi*(flotl+one)*t5)
         endif
      ! set coefficients for adams methods. ----------------------------------
      else if ( me%dat%nq/=1 ) then
         hsum = me%dat%h
         em(1) = one
         flotnq = flotl - one
         do i = 2 , me%dat%l
            em(i) = zero
         enddo
         do j = 1 , nqm1
            if ( (j==nqm1) .and. (me%dat%nqwait==1) ) then
               s = one
               csum = zero
               do i = 1 , nqm1
                  csum = csum + s*em(i)/real(i+1,wp)
                  s = -s
               enddo
               me%dat%tq(1) = em(nqm1)/(flotnq*csum)
            endif
            rxi = me%dat%h/hsum
            do iback = 1 , j
               i = (j+2) - iback
               em(i) = em(i) + em(i-1)*rxi
            enddo
            hsum = hsum + me%dat%tau(j)
         enddo
         ! compute integral from -1 to 0 of polynomial and of x times it. -------
         s = one
         em0 = zero
         csum = zero
         do i = 1 , me%dat%nq
            floti = real(i,wp)
            em0 = em0 + s*em(i)/floti
            csum = csum + s*em(i)/(floti+one)
            s = -s
         enddo
         ! in el, form coefficients of normalized integrated polynomial. --------
         s = one/em0
         me%dat%el(1) = one
         do i = 1 , me%dat%nq
            me%dat%el(i+1) = s*em(i)/real(i,wp)
         enddo
         xi = hsum/me%dat%h
         me%dat%tq(2) = xi*em0/csum
         me%dat%tq(5) = xi/me%dat%el(me%dat%l)
         if ( me%dat%nqwait==1 ) then
            ! for higher order control constant, multiply polynomial by 1+x/xi(q). -
            rxi = one/xi
            do iback = 1 , me%dat%nq
               i = (me%dat%l+1) - iback
               em(i) = em(i) + em(i-1)*rxi
            enddo
            ! compute integral of polynomial. --------------------------------------
            s = one
            csum = zero
            do i = 1 , me%dat%l
               csum = csum + s*em(i)/real(i+1,wp)
               s = -s
            enddo
            me%dat%tq(3) = flotl*em0/csum
         endif
      else
         me%dat%el(1) = one
         me%dat%el(2) = one
         me%dat%tq(1) = one
         me%dat%tq(2) = two
         me%dat%tq(3) = six*me%dat%tq(2)
         me%dat%tq(5) = one
      endif
      me%dat%tq(4) = cortes*me%dat%tq(2)

   end subroutine dvset

!*****************************************************************************************
!>
!  this subroutine adjusts the `yh` array on reduction of order,
!  and also when the order is increased for the stiff option (meth = 2).

   subroutine dvjust(me,yh,ldyh,iord)

      class(dvode_t),intent(inout) :: me
      integer,intent(in) :: ldyh !! leading dimension of `yh`
      integer,intent(in) :: iord !! an integer flag used when meth = 2 to indicate an order
                                 !! increase (iord = +1) or an order decrease (iord = -1).
      real(wp),intent(inout) :: yh(ldyh,*)

      real(wp) :: alph0 , alph1 , hsum , prod , t1 , xi , xiold
      integer :: i , iback , j , jp1 , lp1 , nqm1 , nqm2 , nqp1

      real(wp),parameter :: one = 1.0_wp

      if ( (me%dat%nq==2) .and. (iord/=1) ) return
      nqm1 = me%dat%nq - 1
      nqm2 = me%dat%nq - 2
      if ( me%dat%meth==2 ) then
         !-----------------------------------------------------------------------
         ! stiff option...
         ! check to see if the order is being increased or decreased.
         !-----------------------------------------------------------------------
         if ( iord==1 ) then
            ! order increase. ------------------------------------------------------
            do j = 1 , me%dat%lmax
               me%dat%el(j) = zero
            enddo
            me%dat%el(3) = one
            alph0 = -one
            alph1 = one
            prod = one
            xiold = one
            hsum = me%dat%hscal
            if ( me%dat%nq/=1 ) then
               do j = 1 , nqm1
                  ! construct coefficients of x*x*(x+xi(1))*...*(x+xi(j)). ---------------
                  jp1 = j + 1
                  hsum = hsum + me%dat%tau(jp1)
                  xi = hsum/me%dat%hscal
                  prod = prod*xi
                  alph0 = alph0 - one/real(jp1,wp)
                  alph1 = alph1 + one/xi
                  do iback = 1 , jp1
                     i = (j+4) - iback
                     me%dat%el(i) = me%dat%el(i)*xiold + me%dat%el(i-1)
                  enddo
                  xiold = xi
               enddo
            endif
            t1 = (-alph0-alph1)/prod
            ! load column l + 1 in yh array. ---------------------------------------
            lp1 = me%dat%l + 1
            do i = 1 , me%dat%n
               yh(i,lp1) = t1*yh(i,me%dat%lmax)
            enddo
            ! add correction terms to yh array. ------------------------------------
            nqp1 = me%dat%nq + 1
            do j = 3 , nqp1
               call daxpy(me%dat%n,me%dat%el(j),yh(1,lp1),1,yh(1,j),1)
            enddo
         else
            ! order decrease. ------------------------------------------------------
            do j = 1 , me%dat%lmax
               me%dat%el(j) = zero
            enddo
            me%dat%el(3) = one
            hsum = zero
            do j = 1 , nqm2
               ! construct coefficients of x*x*(x+xi(1))*...*(x+xi(j)). ---------------
               hsum = hsum + me%dat%tau(j)
               xi = hsum/me%dat%hscal
               jp1 = j + 1
               do iback = 1 , jp1
                  i = (j+4) - iback
                  me%dat%el(i) = me%dat%el(i)*xi + me%dat%el(i-1)
               enddo
            enddo
            ! subtract correction terms from yh array. -----------------------------
            do j = 3 , me%dat%nq
               do i = 1 , me%dat%n
                  yh(i,j) = yh(i,j) - yh(i,me%dat%l)*me%dat%el(j)
               enddo
            enddo
            return
         endif
      !-----------------------------------------------------------------------
      ! nonstiff option...
      ! check to see if the order is being increased or decreased.
      !-----------------------------------------------------------------------
      else if ( iord==1 ) then
         ! order increase. ------------------------------------------------------
         ! zero out next column in yh array. ------------------------------------
         lp1 = me%dat%l + 1
         do i = 1 , me%dat%n
            yh(i,lp1) = zero
         enddo
         return
      else
         ! order decrease. ------------------------------------------------------
         do j = 1 , me%dat%lmax
            me%dat%el(j) = zero
         enddo
         me%dat%el(2) = one
         hsum = zero
         do j = 1 , nqm2
            ! construct coefficients of x*(x+xi(1))*...*(x+xi(j)). -----------------
            hsum = hsum + me%dat%tau(j)
            xi = hsum/me%dat%hscal
            jp1 = j + 1
            do iback = 1 , jp1
               i = (j+3) - iback
               me%dat%el(i) = me%dat%el(i)*xi + me%dat%el(i-1)
            enddo
         enddo
         ! construct coefficients of integrated polynomial. ---------------------
         do j = 2 , nqm1
            me%dat%el(j+1) = real(me%dat%nq,wp)*me%dat%el(j)/real(j,wp)
         enddo
         ! subtract correction terms from yh array. -----------------------------
         do j = 3 , me%dat%nq
            do i = 1 , me%dat%n
               yh(i,j) = yh(i,j) - yh(i,me%dat%l)*me%dat%el(j)
            enddo
         enddo
         return
      endif

   end subroutine dvjust

!*****************************************************************************************
!>
!  subroutine dvnlsd is a nonlinear system solver, which uses functional
!  iteration or a chord (modified newton) method.  for the chord method
!  direct linear algebraic system solvers are used.  subroutine dvnlsd
!  then handles the corrector phase of this integration package.

   subroutine dvnlsd(me,y,yh,ldyh,vsav,savf,ewt,acor,iwm,wm,nflag)

      class(dvode_t),intent(inout) :: me
      real(wp),intent(inout) :: y(*) !! the dependent variable, a vector of length n
      integer,intent(in) :: ldyh !! a constant >= n, the first dimension of yh
      real(wp),intent(inout) :: yh(ldyh,*) !! the nordsieck (taylor) array, ldyh by lmax, input
                                           !! and output.  on input, it contains predicted values.
      real(wp) :: vsav(*) !! unused work array.
      real(wp) :: savf(*) !! a work array of length n.
      real(wp),intent(in) :: ewt(*) !! an error weight vector of length n
      real(wp),intent(inout) :: acor(*) !! a work array of length n, used for the accumulated
                                        !! corrections to the predicted y vector.
      real(wp),intent(inout) :: wm(*) !! real work array associated with matrix
                                      !! operations in chord iteration (miter /= 0).
      integer,intent(inout) :: iwm(*) !! integer work array associated with matrix
                                      !! operations in chord iteration (miter /= 0).
      integer,intent(inout) :: nflag !! input/output flag, with values and meanings as follows:
                                     !!
                                     !! * **input:**
                                     !!     *  0 first call for this time step.
                                     !!     * -1 convergence failure in previous call to dvnlsd.
                                     !!     * -2 error test failure in dvstep.
                                     !!
                                     !! * **output:**
                                     !!     *  0 successful completion of nonlinear solver.
                                     !!     * -1 convergence failure or singular matrix.
                                     !!     * -2 unrecoverable error in matrix preprocessing (cannot occur here).
                                     !!     * -3 unrecoverable error in solution (cannot occur here).

      real(wp) :: cscale , dcon , del , delp
      integer :: i , ierpj , iersl , m

      integer,parameter :: maxcor = 3
      integer,parameter :: msbp = 20
      real(wp),parameter :: ccmax = 0.3_wp
      real(wp),parameter :: crdown = 0.3_wp
      real(wp),parameter :: rdiv = 2.0_wp
      real(wp),parameter :: one = 1.0_wp
      real(wp),parameter :: two = 2.0_wp

      !-----------------------------------------------------------------------
      ! on the first step, on a change of method order, or after a
      ! nonlinear convergence failure with nflag = -2, set ipup = miter
      ! to force a jacobian update when miter /= 0.
      !-----------------------------------------------------------------------
      if ( me%dat%jstart==0 ) me%dat%nslp = 0
      if ( nflag==0 ) me%dat%icf = 0
      if ( nflag==-2 ) me%dat%ipup = me%dat%miter
      if ( (me%dat%jstart==0) .or. (me%dat%jstart==-1) ) me%dat%ipup = me%dat%miter
      ! if this is functional iteration, set crate == 1 and drop to 220
      if ( me%dat%miter==0 ) then
         me%dat%crate = one
      else
         !-----------------------------------------------------------------------
         ! rc is the ratio of new to old values of the coefficient h/el(2)=h/l1.
         ! when rc differs from 1 by more than ccmax, ipup is set to miter
         ! to force dvjac to be called, if a jacobian is involved.
         ! in any case, dvjac is called at least every msbp steps.
         !-----------------------------------------------------------------------
         me%dat%drc = abs(me%dat%rc-one)
         if ( me%dat%drc>ccmax .or. me%dat%nst>=me%dat%nslp+msbp ) me%dat%ipup = me%dat%miter
      end if

      corrector : do
         !-----------------------------------------------------------------------
         ! up to maxcor corrector iterations are taken.  a convergence test is
         ! made on the r.m.s. norm of each correction, weighted by the error
         ! weight vector ewt.  the sum of the corrections is accumulated in the
         ! vector acor(i).  the yh array is not altered in the corrector loop.
         !-----------------------------------------------------------------------
         m = 0
         delp = zero
         call dcopy(me%dat%n,yh(1,1),1,y,1)
         call me%f(me%dat%n,me%dat%tn,y(1:me%dat%n),savf(1:me%dat%n))
         me%dat%nfe = me%dat%nfe + 1
         if ( me%dat%ipup>0 ) then
            !-----------------------------------------------------------------------
            ! if indicated, the matrix p = i - h*rl1*j is reevaluated and
            ! preprocessed before starting the corrector iteration.  ipup is set
            ! to 0 as an indicator that this has been done.
            !-----------------------------------------------------------------------
            call me%dvjac(y,yh,ldyh,ewt,acor,savf,wm,iwm,ierpj)
            me%dat%ipup = 0
            me%dat%rc = one
            me%dat%drc = zero
            me%dat%crate = one
            me%dat%nslp = me%dat%nst
            ! if matrix is singular, take error return to force cut in step size. --
            if ( ierpj/=0 ) exit corrector
         endif
         do i = 1 , me%dat%n
            acor(i) = zero
         enddo

         do
            ! this is a looping point for the corrector iteration. -----------------
            if ( me%dat%miter/=0 ) then
               !-----------------------------------------------------------------------
               ! in the case of the chord method, compute the corrector error,
               ! and solve the linear system with that as right-hand side and
               ! p as coefficient matrix.  the correction is scaled by the factor
               ! 2/(1+rc) to account for changes in h*rl1 since the last dvjac call.
               !-----------------------------------------------------------------------
               do i = 1 , me%dat%n
                  y(i) = (me%dat%rl1*me%dat%h)*savf(i) - (me%dat%rl1*yh(i,2)+acor(i))
               enddo
               call me%dvsol(wm,iwm,y,iersl)
               me%dat%nni = me%dat%nni + 1
               if ( iersl>0 ) exit
               if ( me%dat%meth==2 .and. me%dat%rc/=one ) then
                  cscale = two/(one+me%dat%rc)
                  call dscal(me%dat%n,cscale,y,1)
               endif
               del = me%dvnorm(me%dat%n,y(1:me%dat%n),ewt(1:me%dat%n))
               call daxpy(me%dat%n,one,y,1,acor,1)
               do i = 1 , me%dat%n
                  y(i) = yh(i,1) + acor(i)
               enddo
            else
               !-----------------------------------------------------------------------
               ! in the case of functional iteration, update y directly from
               ! the result of the last function evaluation.
               !-----------------------------------------------------------------------
               do i = 1 , me%dat%n
                  savf(i) = me%dat%rl1*(me%dat%h*savf(i)-yh(i,2))
               enddo
               do i = 1 , me%dat%n
                  y(i) = savf(i) - acor(i)
               enddo
               del = me%dvnorm(me%dat%n,y(1:me%dat%n),ewt(1:me%dat%n))
               do i = 1 , me%dat%n
                  y(i) = yh(i,1) + savf(i)
               enddo
               call dcopy(me%dat%n,savf,1,acor,1)
            endif
            !-----------------------------------------------------------------------
            ! test for convergence.  if m > 0, an estimate of the convergence
            ! rate constant is stored in crate, and this is used in the test.
            !-----------------------------------------------------------------------
            if ( m/=0 ) me%dat%crate = max(crdown*me%dat%crate,del/delp)
            dcon = del*min(one,me%dat%crate)/me%dat%tq(4)
            if ( dcon<=one ) then
               ! return for successful step. ------------------------------------------
               nflag = 0
               me%dat%jcur = 0
               me%dat%icf = 0
               if ( m==0 ) me%dat%acnrm = del
               if ( m>0 ) me%dat%acnrm = me%dvnorm(me%dat%n,acor(1:me%dat%n),ewt(1:me%dat%n))
               return
            else
               m = m + 1
               if ( m/=maxcor ) then
                  if ( m<2 .or. del<=rdiv*delp ) then
                     delp = del
                     call me%f(me%dat%n,me%dat%tn,y(1:me%dat%n),savf(1:me%dat%n))
                     me%dat%nfe = me%dat%nfe + 1
                     cycle
                  endif
               endif
               exit
            endif
         end do

         if ( me%dat%miter/=0 .and. me%dat%jcur/=1 ) then
            me%dat%icf = 1
            me%dat%ipup = me%dat%miter
         else
            exit corrector
         endif

      end do corrector

      nflag = -1
      me%dat%icf = 2
      me%dat%ipup = me%dat%miter

   end subroutine dvnlsd

!*****************************************************************************************
!>
!  dvjac is called by [[dvnlsd]] to compute and process the matrix
!  `p = i - h*rl1*j` , where j is an approximation to the jacobian.
!
!  here j is computed by the user-supplied routine `jac` if
!  miter = 1 or 4, or by finite differencing if miter = 2, 3, or 5.
!  if miter = 3, a diagonal approximation to j is used.
!  if jsv = -1, j is computed from scratch in all cases.
!  if jsv = 1 and miter = 1, 2, 4, or 5, and if the saved value of j is
!  considered acceptable, then p is constructed from the saved j.
!  j is stored in wm and replaced by p.  if miter /= 3, p is then
!  subjected to lu decomposition in preparation for later solution
!  of linear systems with p as coefficient matrix. this is done
!  by [[dgefa]] if miter = 1 or 2, and by dgbfa if miter = 4 or 5.

   subroutine dvjac(me,y,yh,ldyh,ewt,ftem,savf,wm,iwm,ierpj)

      class(dvode_t),intent(inout) :: me
      real(wp),intent(inout) :: y(*) !! vector containing predicted values on entry.
      integer,intent(in) :: ldyh !! a constant >= n, the first dimension of `yh`.
      real(wp),intent(in) :: yh(ldyh,*) !! the nordsieck array, an `ldyh` by `lmax` array.
      real(wp),intent(in) :: ewt(*) !! an error weight vector of length `n`.
      real(wp),intent(in) :: ftem(*)
      real(wp),intent(in) :: savf(*) !! array containing f evaluated at predicted y.
      real(wp),intent(inout) :: wm(*) !! real work space for matrices.  in the output, it contains
                                      !! the inverse diagonal matrix if miter = 3 and the lu
                                      !! decomposition of p if miter is 1, 2 , 4, or 5.
                                      !! storage of matrix elements starts at wm(3).
                                      !! storage of the saved jacobian starts at wm(locjs).
                                      !! wm also contains the following matrix-related data:
                                      !! wm(1) = sqrt(uround), used in numerical jacobian step.
                                      !! wm(2) = h*rl1, saved for later use if miter = 3.
      integer,intent(inout) :: iwm(*) !! integer work space containing pivot information,
                                      !! starting at iwm(31), if miter is 1, 2, 4, or 5.
                                      !! iwm also contains band parameters ml = iwm(1) and
                                      !! mu = iwm(2) if miter is 4 or 5.
      integer,intent(out) :: ierpj !! output error flag,  = 0 if no trouble, 1 if the p
                                   !! matrix is found to be singular.

      real(wp) :: con , di , fac , hrl1 , r , r0 , &
                  srur , yi , yj , yjj
      integer :: i , i1 , i2 , ier , ii , j , j1 , jj , jok , lenp , mba , &
                 mband , meb1 , meband , ml , ml3 , mu , np1

      real(wp),parameter :: one = 1.0_wp
      real(wp),parameter :: thou = 1000.0_wp
      real(wp),parameter :: pt1 = 0.1_wp

      ierpj = 0
      hrl1 = me%dat%h*me%dat%rl1 ! rl1 = 1/el(2)
      ! see whether j should be evaluated (jok = -1) or not (jok = 1). -------
      jok = me%dat%jsv
      if ( me%dat%jsv==1 ) then
         if ( me%dat%nst==0 .or. me%dat%nst>me%dat%nslj+me%dat%msbj ) jok = -1
         if ( me%dat%icf==1 .and. me%dat%drc<me%dat%ccmxj ) jok = -1
         if ( me%dat%icf==2 ) jok = -1
      endif
      ! end of setting jok. --------------------------------------------------

      if ( jok==-1 .and. me%dat%miter==1 ) then
         ! if jok = -1 and miter = 1, call me%jac to evaluate jacobian. ------------
         me%dat%nje = me%dat%nje + 1
         me%dat%nslj = me%dat%nst
         me%dat%jcur = 1
         lenp = me%dat%n*me%dat%n
         do i = 1 , lenp
            wm(i+2) = zero
         enddo
         call me%jac(me%dat%n,me%dat%tn,y(1:me%dat%n),0,0,wm(3),me%dat%n)
         if ( me%dat%jsv==1 ) call dcopy(lenp,wm(3),1,wm(me%dat%locjs),1)
      endif

      if ( jok==-1 .and. me%dat%miter==2 ) then
         ! if miter = 2, make n calls to f to approximate the jacobian. ---------
         me%dat%nje = me%dat%nje + 1
         me%dat%nslj = me%dat%nst
         me%dat%jcur = 1
         fac = me%dvnorm(me%dat%n,savf(1:me%dat%n),ewt(1:me%dat%n))
         r0 = thou*abs(me%dat%h)*me%dat%uround*real(me%dat%n,wp)*fac
         if ( r0==zero ) r0 = one
         srur = wm(1)
         j1 = 2
         do j = 1 , me%dat%n
            yj = y(j)
            r = max(srur*abs(yj),r0/ewt(j))
            y(j) = y(j) + r
            fac = one/r
            call me%f(me%dat%n,me%dat%tn,y(1:me%dat%n),ftem(1:me%dat%n))
            do i = 1 , me%dat%n
               wm(i+j1) = (ftem(i)-savf(i))*fac
            enddo
            y(j) = yj
            j1 = j1 + me%dat%n
         enddo
         me%dat%nfe = me%dat%nfe + me%dat%n
         lenp = me%dat%n*me%dat%n
         if ( me%dat%jsv==1 ) call dcopy(lenp,wm(3),1,wm(me%dat%locjs),1)
      endif

      if ( jok==1 .and. (me%dat%miter==1 .or. me%dat%miter==2) ) then
         me%dat%jcur = 0
         lenp = me%dat%n*me%dat%n
         call dcopy(lenp,wm(me%dat%locjs),1,wm(3),1)
      endif

      if ( me%dat%miter==1 .or. me%dat%miter==2 ) then
         ! multiply jacobian by scalar, add identity, and do lu decomposition. --
         con = -hrl1
         call dscal(lenp,con,wm(3),1)
         j = 3
         np1 = me%dat%n + 1
         do i = 1 , me%dat%n
            wm(j) = wm(j) + one
            j = j + np1
         enddo
         me%dat%nlu = me%dat%nlu + 1
         call dgefa(wm(3),me%dat%n,me%dat%n,iwm(31),ier)
         if ( ier/=0 ) ierpj = 1
         return
      endif
      ! end of code block for miter = 1 or 2. --------------------------------

      if ( me%dat%miter==3 ) then
         ! if miter = 3, construct a diagonal approximation to j and p. ---------
         me%dat%nje = me%dat%nje + 1
         me%dat%jcur = 1
         wm(2) = hrl1
         r = me%dat%rl1*pt1
         do i = 1 , me%dat%n
            y(i) = y(i) + r*(me%dat%h*savf(i)-yh(i,2))
         enddo
         call me%f(me%dat%n,me%dat%tn,y(1:me%dat%n),wm(3))
         me%dat%nfe = me%dat%nfe + 1
         do i = 1 , me%dat%n
            r0 = me%dat%h*savf(i) - yh(i,2)
            di = pt1*r0 - me%dat%h*(wm(i+2)-savf(i))
            wm(i+2) = one
            if ( abs(r0)>=me%dat%uround/ewt(i) ) then
               if ( abs(di)==zero ) then
                  ierpj = 1
                  return
               end if
               wm(i+2) = pt1*r0/di
            endif
         enddo
         return
      endif
      ! end of code block for miter = 3. -------------------------------------

      ! set constants for miter = 4 or 5. ------------------------------------
      ml = iwm(1)
      mu = iwm(2)
      ml3 = ml + 3
      mband = ml + mu + 1
      meband = mband + ml
      lenp = meband*me%dat%n

      if ( jok==-1 .and. me%dat%miter==4 ) then
         ! if jok = -1 and miter = 4, call me%jac to evaluate jacobian. ------------
         me%dat%nje = me%dat%nje + 1
         me%dat%nslj = me%dat%nst
         me%dat%jcur = 1
         do i = 1 , lenp
            wm(i+2) = zero
         enddo
         call me%jac(me%dat%n,me%dat%tn,y(1:me%dat%n),ml,mu,wm(ml3),meband)
         if ( me%dat%jsv==1 ) call dacopy(mband,me%dat%n,wm(ml3),meband,wm(me%dat%locjs),mband)
      endif

      if ( jok==-1 .and. me%dat%miter==5 ) then
         ! if miter = 5, make ml+mu+1 calls to f to approximate the jacobian. ---
         me%dat%nje = me%dat%nje + 1
         me%dat%nslj = me%dat%nst
         me%dat%jcur = 1
         mba = min(mband,me%dat%n)
         meb1 = meband - 1
         srur = wm(1)
         fac = me%dvnorm(me%dat%n,savf(1:me%dat%n),ewt(1:me%dat%n))
         r0 = thou*abs(me%dat%h)*me%dat%uround*real(me%dat%n,wp)*fac
         if ( r0==zero ) r0 = one
         do j = 1 , mba
            do i = j , me%dat%n , mband
               yi = y(i)
               r = max(srur*abs(yi),r0/ewt(i))
               y(i) = y(i) + r
            enddo
            call me%f(me%dat%n,me%dat%tn,y(1:me%dat%n),ftem(1:me%dat%n))
            do jj = j , me%dat%n , mband
               y(jj) = yh(jj,1)
               yjj = y(jj)
               r = max(srur*abs(yjj),r0/ewt(jj))
               fac = one/r
               i1 = max(jj-mu,1)
               i2 = min(jj+ml,me%dat%n)
               ii = jj*meb1 - ml + 2
               do i = i1 , i2
                  wm(ii+i) = (ftem(i)-savf(i))*fac
               enddo
            enddo
         enddo
         me%dat%nfe = me%dat%nfe + mba
         if ( me%dat%jsv==1 ) call dacopy(mband,me%dat%n,wm(ml3),meband,wm(me%dat%locjs),mband)
      endif

      if ( jok==1 ) then
         me%dat%jcur = 0
         call dacopy(mband,me%dat%n,wm(me%dat%locjs),mband,wm(ml3),meband)
      endif

      ! multiply jacobian by scalar, add identity, and do lu decomposition.
      con = -hrl1
      call dscal(lenp,con,wm(3),1)
      ii = mband + 2
      do i = 1 , me%dat%n
         wm(ii) = wm(ii) + one
         ii = ii + meband
      enddo
      me%dat%nlu = me%dat%nlu + 1
      call dgbfa(wm(3),meband,me%dat%n,ml,mu,iwm(31),ier)
      if ( ier/=0 ) ierpj = 1
      ! end of code block for miter = 4 or 5. --------------------------------

   end subroutine dvjac

!*****************************************************************************************
!>
! this routine copies one rectangular array, `a`, to another, `b`,
! where `a` and `b` may have different row dimensions, `nrowa` and `nrowb`.
! the data copied consists of `nrow` rows and `ncol` columns.

   subroutine dacopy(nrow,ncol,a,nrowa,b,nrowb)

      integer,intent(in) :: nrow
      integer,intent(in) :: ncol
      integer,intent(in) :: nrowa
      integer,intent(in) :: nrowb
      real(wp),intent(in) :: a(nrowa,ncol)
      real(wp),intent(out) :: b(nrowb,ncol)

      integer :: ic

      do ic = 1 , ncol
         call dcopy(nrow,a(1,ic),1,b(1,ic),1)
      enddo

   end subroutine dacopy

!*****************************************************************************************
!>
!  This routine manages the solution of the linear system arising from
!  a chord iteration.  it is called if `miter /= 0`:
!
!   * if miter is 1 or 2, it calls [[dgesl]] to accomplish this.
!   * if miter = 3 it updates the coefficient `h*rl1` in the diagonal
!     matrix, and then computes the solution.
!   * if miter is 4 or 5, it calls [[dgbsl]].

   subroutine dvsol(me,wm,iwm,x,iersl)

      class(dvode_t),intent(inout) :: me
      real(wp),intent(inout) :: wm(*) !! real work space containing the inverse diagonal matrix if
                                      !! miter = 3 and the lu decomposition of the matrix otherwise.
                                      !! storage of matrix elements starts at wm(3).
                                      !! wm also contains the following matrix-related data:
                                      !! wm(1) = sqrt(uround) (not used here),
                                      !! wm(2) = hrl1, the previous value of h*rl1, used if miter = 3.
      real(wp),intent(inout) :: x(*) !! the right-hand side vector on input,
                                     !! and the solution vector
                                     !! on output, of length n.
      integer :: iwm(*) !! integer work space containing pivot information, starting at
                        !! iwm(31), if miter is 1, 2, 4, or 5.  iwm also contains band
                        !! parameters ml = iwm(1) and mu = iwm(2) if miter is 4 or 5.
      integer,intent(out) :: iersl !! output flag.  iersl = 0 if no trouble occurred.
                                   !! iersl = 1 if a singular matrix arose with miter = 3.

      integer :: i , meband , ml , mu
      real(wp) :: di , hrl1 , phrl1 , r

      real(wp),parameter :: one = 1.0_wp

      iersl = 0
      select case (me%dat%miter)
      case (3)
         phrl1 = wm(2)
         hrl1 = me%dat%h*me%dat%rl1
         wm(2) = hrl1
         if ( hrl1/=phrl1 ) then
            r = hrl1/phrl1
            do i = 1 , me%dat%n
               di = one - r*(one-one/wm(i+2))
               if ( abs(di)==zero ) then
                  iersl = 1
                  return
               end if
               wm(i+2) = one/di
            enddo
         endif
         do i = 1 , me%dat%n
            x(i) = wm(i+2)*x(i)
         enddo
      case (4,5)
         ml = iwm(1)
         mu = iwm(2)
         meband = 2*ml + mu + 1
         call dgbsl(wm(3),meband,me%dat%n,ml,mu,iwm(31),x,0)
      case default
         call dgesl(wm(3),me%dat%n,me%dat%n,iwm(31),x,0)
      end select

   end subroutine dvsol

!*****************************************************************************************
!>
!  this routine saves or restores (depending on `job`) the contents of the
!  dvode internal variables.

   subroutine dvsrco(me,sav,job)

      class(dvode_t),intent(inout) :: me
      type(dvode_data_t),intent(inout) :: sav
      integer,intent(in) :: job !! flag indicating to save or restore the data:
                                !!
                                !! * `job = 1` if common is to be saved (written to `sav`).
                                !! * `job = 2` if common is to be restored (read from `sav`).
                                !!
                                !! a call with job = 2 presumes a prior call with job = 1.

      select case (job)
      case(1); sav = me%dat
      case(2); me%dat = sav
      case default
         error stop 'invalid input to dvsrco'
      end select

   end subroutine dvsrco

!*****************************************************************************************
!>
!  Set error weight vector.
!
!  this subroutine sets the error weight vector ewt according to
!```
!      ewt(i) = rtol(i)*abs(ycur(i)) + atol(i),  i = 1,...,n,
!```
!  with the subscript on rtol and/or atol possibly replaced by 1 above,
!  depending on the value of itol.
!
!### Author
!  * hindmarsh, alan c., (llnl)
!
!### Revision history
!  * 791129  date written
!  * 890501  modified prologue to slatec/ldoc format.  (fnf)
!  * 890503  minor cosmetic changes.  (fnf)
!  * 930809  renamed to allow single/real(wp) versions. (ach)

   subroutine dewset_default(me,n,itol,rtol,atol,ycur,ewt)

      class(dvode_t),intent(inout) :: me
      integer,intent(in) :: n
      integer,intent(in) :: itol
      real(wp),intent(in) :: rtol(*)
      real(wp),intent(in) :: atol(*)
      real(wp),intent(in) :: ycur(n)
      real(wp),intent(out) :: ewt(n)

      integer :: i

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
         return
      case default
      end select

      do i = 1 , n
         ewt(i) = rtol(1)*abs(ycur(i)) + atol(1)
      enddo

   end subroutine dewset_default

!*****************************************************************************************
!>
!  weighted root-mean-square vector norm.
!
!  this function routine computes the weighted root-mean-square norm
!  of the vector of length n contained in the array v, with weights
!  contained in the array w of length n:
!```
!    dvnorm = sqrt( (1/n) * sum( v(i)*w(i) )**2 )
!```
!
!### Author
!  * hindmarsh, alan c., (llnl)
!
!### Revision history
!  * 791129  date written
!  * 890501  modified prologue to slatec/ldoc format.  (fnf)
!  * 890503  minor cosmetic changes.  (fnf)
!  * 930809  renamed to allow single/real(wp) versions. (ach)

      real(wp) function dvnorm_default(me,n,v,w)

      class(dvode_t),intent(inout) :: me
      integer,intent(in) :: n
      real(wp),intent(in) :: v(n)
      real(wp),intent(in) :: w(n)

      integer :: i
      real(wp) :: sum

      sum = zero
      do i = 1 , n
         sum = sum + (v(i)*w(i))**2
      enddo
      dvnorm_default = sqrt(sum/n)

      end function dvnorm_default

!*****************************************************************************************
!>
!  write error message with values.
!
!  subroutines [[xerrwd]], [[xsetf]], [[xsetun]], and the function routine [[ixsav]],
!  as given here, constitute a simplified version of the slatec error
!  handling package.
!
!### Note
!
!  This routine is machine-dependent and specialized for use
!  in limited context, in the following ways:
!
!   1. the argument msg is assumed to be of type character, and
!      the message is printed with a format of (1x,a).
!   2. the message is assumed to take only one line.
!      multi-line messages are generated by repeated calls.
!   3. if level = 2, control passes to the statement   stop
!      to abort the run.  this statement may be machine-dependent.
!   4. r1 and r2 are assumed to be in real(wp) and are printed
!      in d21.13 format.
!
!### Author
!  * hindmarsh, alan c., (llnl)
!
!### Revision history
!  * 920831  date written
!  * 921118  replaced mflgsv/lunsav by ixsav. (ach)
!  * 930329  modified prologue to slatec format. (fnf)
!  * 930407  changed msg from character*1 array to variable. (fnf)
!  * 930922  minor cosmetic change. (fnf)
!
!### internal notes:
!
!  for a different default logical unit number, [[ixsav]] (or a subsidiary
!  routine that it calls) will need to be modified.
!  for a different run-abort command, change the statement at the end.

   subroutine xerrwd(me,msg,nmes,nerr,level,ni,i1,i2,nr,r1,r2)

      class(dvode_t),intent(inout) :: me
      character(len=*),intent(in) :: msg !! the message (character array).
      integer,intent(in) :: nmes !! the length of msg (number of characters).
      integer,intent(in) :: nerr !! the error number (not used).
      integer,intent(in) :: level !! the error level:
                                  !!
                                  !!  * 0 or 1 means recoverable (control returns to caller).
                                  !!  * 2 means fatal (run is aborted).
      integer,intent(in) :: ni !! number of integers (0, 1, or 2) to be printed with message.
      integer,intent(in) :: i1 !! integer to be printed, depending on `ni`.
      integer,intent(in) :: i2 !! integer to be printed, depending on `ni`.
      integer,intent(in) :: nr !! number of reals (0, 1, or 2) to be printed with message.
      real(wp),intent(in) :: r1 !! real to be printed, depending on `nr`.
      real(wp),intent(in) :: r2 !! real to be printed, depending on `nr`.

      integer :: lunit , mesflg

      ! get logical unit number and message print flag.
      lunit = me%ixsav(1,0,.false.)
      mesflg = me%ixsav(2,0,.false.)

      if ( mesflg/=0 ) then
         ! write the message.
         write (lunit,'(A)') trim(msg)
         select case (ni)
         case(1)
            write (lunit,'(6X,A,I10)') 'in above message,  i1 =', i1
         case(2)
            write (lunit,'(6X,A,I10,3X,A,I10)') 'in above message,  i1 =', i1, 'i2 =', i2
         end select
         select case (nr)
         case(1)
            write (lunit,'(6X,A,D21.13)') 'in above message,  r1 =', r1
         case (2)
            write (lunit,'(6X,A,D21.13,3X,A,D21.13)') 'in above message,  r1 =', r1, 'r2 =', r2
         end select
      endif

      ! abort the run if level = 2.
      if ( level==2 ) stop

   end subroutine xerrwd

!*****************************************************************************************
!>
!  reset the error print control flag.
!
!  xsetf sets the error print control flag to mflag:
!
!  * mflag=1 means print all messages (the default).
!  * mflag=0 means no printing.
!
!### Author
!  * hindmarsh, alan c., (llnl)
!
!### Revision history
!  * 921118  date written
!  * 930329  added slatec format prologue. (fnf)
!  * 930407  corrected see also section. (fnf)
!  * 930922  made user-callable, and other cosmetic changes. (fnf)

   subroutine xsetf(me,mflag)

      class(dvode_t),intent(inout) :: me
      integer :: mflag

      integer :: junk

      if ( mflag==0 .or. mflag==1 ) junk = me%ixsav(2,mflag,.true.)

   end subroutine xsetf

!*****************************************************************************************
!>
!  reset the logical unit number for error messages.
!
!  xsetun sets the logical unit number for error messages to lun.
!
!### Author
!  * hindmarsh, alan c., (llnl)
!
!### Revision history
!  * 921118  date written
!  * 930329  added slatec format prologue. (fnf)
!  * 930407  corrected see also section. (fnf)
!  * 930922  made user-callable, and other cosmetic changes. (fnf)

   subroutine xsetun(me,lun)

      class(dvode_t),intent(inout) :: me
      integer :: lun

      integer :: junk

      if ( lun>0 ) junk = me%ixsav(1,lun,.true.)

   end subroutine xsetun

!*****************************************************************************************
!>
!  save and recall error message control parameters.
!
!  [[ixsav]] saves and recalls one of two error message parameters:
!
!  * `lunit`, the logical unit number to which messages are printed, and
!  * `mesflg`, the message print flag.
!
!  this is a modification of the slatec library routine `j4save`.
!
!  on return:
!
!  * `ixsav` = the (old) value of the parameter.
!
!### See also
!  * [[xerrwd]]
!  * [[xerrwv]]
!
!### Author
!  * hindmarsh, alan c., (llnl)
!
!### Revision history
!  * 921118  date written
!  * 930329  modified prologue to slatec format. (fnf)
!  * 930915  added iumach call to get default output unit.  (ach)
!  * 930922  minor cosmetic changes. (fnf)
!  * 010425  type declaration for iumach added. (ach)

   integer function ixsav(me,ipar,ivalue,iset)

      class(dvode_t),intent(inout) :: me
      integer,intent(in) :: ipar !! parameter indicator (1 for `lunit`, 2 for `mesflg`).
      integer :: ivalue !! the value to be set for the parameter, if iset = .true.
      logical,intent(in) :: iset !! logical flag to indicate whether to read or write.
                                 !! if `iset = .true.`, the parameter will be given
                                 !! the value `ivalue`.  if `iset = .false.`, the parameter
                                 !! will be unchanged, and `ivalue` is a dummy argument.

      select case (ipar)
      case ( 1 )
         if ( me%lunit==-1 ) me%lunit = iumach
         ixsav = me%lunit
         if ( iset ) me%lunit = ivalue
      case ( 2 )
         ixsav = me%mesflg
         if ( iset ) me%mesflg = ivalue
      end select

   end function ixsav

end module dvode_module