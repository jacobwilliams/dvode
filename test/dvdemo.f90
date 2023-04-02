!*****************************************************************************************
!>
!  Demonstration program for the DVODE package.
!  This is the version of 30 April 2002.
!
!  This version is in double precision, suitable for short wordlength
!  computers, PCs, work-stations, etc.
!
!  The package is used to solve two simple problems,
!  one with a full Jacobian, the other with a banded Jacobian,
!  with all 12 of the appropriate values of MF in each case.
!  If the errors are too large, or other difficulty occurs,
!  a warning message is printed.  All output is on unit LOUT = 6.

      program dvdemo
      use dvode_module, wp => dvode_wp
      use iso_fortran_env, only: output_unit

      implicit none

      type(dvode_t) :: solver
      INTEGER :: i , iopar , iopt , iout , istate , itask , &
                 itol , iwork , jsv , leniw , lenrw , liw , lrw ,   &
                 mband , meth , mf , miter , ml , mu , ncfn , neq , nerr , &
                 netf , nfe , nfea , nje , nlu , nni , nout , nqu , nst
      real(wp) :: atol(1) , er , erm , ero , hu , &
                  rtol(1) , rwork , t , tout , y
      DIMENSION y(25) , rwork(847) , iwork(55)

      real(wp),parameter :: tout1 = 1.39283880203_wp
      real(wp),parameter :: dtout = 2.214773875_wp
      integer,parameter :: lout = output_unit
!
      rwork = 0
      iwork = 0
      nerr = 0
      itol = 1
      rtol = 0.0_wp
      atol = 1.0e-6_wp
      lrw = 847
      liw = 55
      iopt = 0
!
! First problem
!
      neq = 2
      nout = 4
      call solver%initialize(f=f1, jac=jac1)
      WRITE (lout,99001) neq , itol , rtol , atol
99001 FORMAT (' Demonstration program for DVODE package'//              &
              ' Problem 1:   Van der Pol oscillator:'/                  &
              '   xdotdot - 3*(1 - x**2)*xdot + x = 0,',                &
              '   x(0) = 2, xdot(0) = 0'/'   NEQ =',I2/'   ITOL =',I3,  &
              '   RTOL =',D10.1,'   ATOL =',D10.1)

      DO jsv = 1 , -1 , -2
         DO meth = 1 , 2
            DO miter = 0 , 3
               IF ( jsv>=0 .OR. miter/=0 ) THEN
                  IF ( jsv>=0 .OR. miter/=3 ) THEN
                     mf = jsv*(10*meth+miter)
                     WRITE (lout,99002) mf
99002                FORMAT (//70('*')//' Solution with MF =',I4//6X,   &
                                        't',15X,'x',15X,'xdot',7X,'NQ',6X,'H'/)
                     t = 0.0_wp
                     y(1) = 2.0_wp
                     y(2) = 0.0_wp
                     itask = 1
                     istate = 1
                     tout = tout1
                     ero = 0.0_wp
                     DO iout = 1 , nout
                        CALL solver%solve(neq,y,t,tout,itol,rtol,atol,  &
                                          itask,istate,iopt,rwork,lrw,iwork, &
                                          liw,mf)
                        hu = rwork(11)
                        nqu = iwork(14)
                        WRITE (lout,99003) t , y(1) , y(2) , nqu , hu
99003                   FORMAT (1X,D15.5,D16.5,D14.3,I5,D14.3)
                        IF ( istate<0 ) GOTO 2
                        iopar = iout - 2*(iout/2)
                        IF ( iopar==0 ) THEN
                           er = ABS(y(1))/atol(1)
                           ero = MAX(ero,er)
                           IF ( er>=10000.0_wp ) THEN
                              WRITE (lout,99008)
                              nerr = nerr + 1
                           ENDIF
                        ENDIF
                        tout = tout + dtout
                     ENDDO
 2                   IF ( istate<0 ) nerr = nerr + 1
                     nst = iwork(11)
                     nfe = iwork(12)
                     nje = iwork(13)
                     nlu = iwork(19)
                     lenrw = iwork(17)
                     leniw = iwork(18)
                     nni = iwork(20)
                     ncfn = iwork(21)
                     netf = iwork(22)
                     nfea = nfe
                     IF ( miter==2 ) nfea = nfe - neq*nje
                     IF ( miter==3 ) nfea = nfe - nje
                     WRITE (lout,99009) lenrw , leniw , nst , nfe ,     &
                                      & nfea , nje , nlu , nni , ncfn , &
                                      & netf , ero
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO
!
! Second problem

      call solver%initialize(f=f2, jac=jac2)
      neq = 25
      ml = 5
      mu = 0
      iwork(1) = ml
      iwork(2) = mu
      mband = ml + mu + 1
      nout = 5
      WRITE (lout,99004) neq , ml , mu , itol , rtol , atol
99004 FORMAT (//70('*')//' Demonstration program for DVODE package'//   &
              ' Problem 2: ydot = A * y , where',                       &
              ' A is a banded lower triangular matrix'/                 &
              '   derived from 2-D advection PDE'/'   NEQ =',I3,        &
              '   ML =',I2,'   MU =',I2/'   ITOL =',I3,'   RTOL =',     &
              D10.1,'   ATOL =',D10.1)
!
      DO jsv = 1 , -1 , -2
         DO meth = 1 , 2
            DO miter = 0 , 5
               IF ( miter/=1 .AND. miter/=2 ) THEN
                  IF ( jsv>=0 .OR. miter/=0 ) THEN
                     IF ( jsv>=0 .OR. miter/=3 ) THEN
                        mf = jsv*(10*meth+miter)
                        WRITE (lout,99005) mf
99005                   FORMAT (//70('*')//' Solution with MF =',I4//6X,&
                               &'t',13X,'Max.Err.',5X,'NQ',6X,'H'/)
                        t = 0.0_wp
                        DO i = 2 , neq
                           y(i) = 0.0_wp
                        ENDDO
                        y(1) = 1.0_wp
                        itask = 1
                        istate = 1
                        tout = 0.01_wp
                        ero = 0.0D0
                        DO iout = 1 , nout
                           CALL solver%solve(neq,y,t,tout,itol,rtol, &
                                    & atol,itask,istate,iopt,rwork,  &
                                    & lrw,iwork,liw,mf)
                           CALL EDIT2(y,t,erm)
                           hu = rwork(11)
                           nqu = iwork(14)
                           WRITE (lout,99006) t , erm , nqu , hu
99006                      FORMAT (1X,D15.5,D14.3,I5,D14.3)
                           IF ( istate<0 ) GOTO 4
                           er = erm/atol(1)
                           ero = MAX(ero,er)
                           IF ( er>1000.0_wp ) THEN
                              WRITE (lout,99008)
                              nerr = nerr + 1
                           ENDIF
                           tout = tout*10.0_wp
                        ENDDO
 4                      IF ( istate<0 ) nerr = nerr + 1
                        nst = iwork(11)
                        nfe = iwork(12)
                        nje = iwork(13)
                        nlu = iwork(19)
                        lenrw = iwork(17)
                        leniw = iwork(18)
                        nni = iwork(20)
                        ncfn = iwork(21)
                        netf = iwork(22)
                        nfea = nfe
                        IF ( miter==5 ) nfea = nfe - mband*nje
                        IF ( miter==3 ) nfea = nfe - nje
                        WRITE (lout,99009) lenrw , leniw , nst , nfe ,  &
                             & nfea , nje , nlu , nni , ncfn , netf ,   &
                             & ero
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO
!
      WRITE (lout,99007) nerr
99007 FORMAT (//70('*')//' Number of errors encountered =',I3)
99008 FORMAT (//' Warning: Error exceeds 10000 * tolerance'//)
99009 FORMAT (//' Final statistics for this run:'/' RWORK size =',I4,   &
             &'   IWORK size =',I4/' Number of steps =',                &
             &I5/' Number of f-s   =',I5/' (excluding J-s) =',          &
             &I5/' Number of J-s   =',I5/' Number of LU-s  =',          &
             &I5/' Number of nonlinear iterations =',                   &
             &I5/' Number of nonlinear convergence failures =',         &
             &I5/' Number of error test failures =',                    &
             &I5/' Error overrun =',D10.2)

      contains

      subroutine f1(me,neq,t,y,ydot)

      class(dvode_t),intent(inout) :: me
      integer :: neq
      real(wp) :: t , y(neq) , ydot(neq)

      ydot(1) = y(2)
      ydot(2) = 3.0_wp*(1.0_wp-y(1)*y(1))*y(2) - y(1)

      end subroutine f1

      subroutine jac1(me,neq,t,y,ml,mu,pd,nrowpd)

      class(dvode_t),intent(inout) :: me
      integer :: neq , ml , mu , nrowpd
      real(wp) :: t , y(neq) , pd(nrowpd,neq)

      pd(1,1) = 0.0_wp
      pd(1,2) = 1.0_wp
      pd(2,1) = -6.0_wp*y(1)*y(2) - 1.0_wp
      pd(2,2) = 3.0_wp*(1.0_wp-y(1)*y(1))

      end subroutine jac1

      subroutine f2(me,neq,t,y,ydot)

      class(dvode_t),intent(inout) :: me
      integer :: neq , i , j , k 
      real(wp) :: t , y(neq) , ydot(neq) , d

      real(wp),parameter :: alph1 = 1.0_wp 
      real(wp),parameter :: alph2 = 1.0_wp 
      integer,parameter :: ng = 5

      do j = 1 , ng
         do i = 1 , ng
            k = i + (j-1)*ng
            d = -2.0_wp*y(k)
            if ( i/=1 ) d = d + y(k-1)*alph1
            if ( j/=1 ) d = d + y(k-ng)*alph2
            ydot(k) = d
         enddo
      enddo

      end subroutine f2

      subroutine jac2(me,neq,t,y,ml,mu,pd,nrowpd)

      class(dvode_t),intent(inout) :: me
      integer :: neq , ml , mu , nrowpd , j , mband , mu1 , mu2
      real(wp) :: t , y(neq) , pd(nrowpd,neq)  

      real(wp),parameter :: alph1 = 1.0_wp 
      real(wp),parameter :: alph2 = 1.0_wp 
      integer,parameter :: ng = 5

      mband = ml + mu + 1
      mu1 = mu + 1
      mu2 = mu + 2
      do j = 1 , neq
         pd(mu1,j) = -2.0_wp
         pd(mu2,j) = alph1
         pd(mband,j) = alph2
      enddo
      do j = ng , neq , ng
         pd(mu2,j) = 0.0_wp
      enddo

      end subroutine jac2

      subroutine edit2(y,t,erm)
      implicit none
      integer :: i , j , k 
      real(wp) :: y(25) , t , erm , a1 , a2 , er , ex , yt

      real(wp),parameter :: alph1 = 1.0_wp 
      real(wp),parameter :: alph2 = 1.0_wp 
      integer,parameter :: ng = 5

      erm = 0.0_wp
      if ( t==0.0_wp ) return
      ex = 0.0_wp
      if ( t<=30.0_wp ) ex = exp(-2.0_wp*t)
      a2 = 1.0_wp
      do j = 1 , ng
         a1 = 1.0_wp
         do i = 1 , ng
            k = i + (j-1)*ng
            yt = t**(i+j-2)*ex*a1*a2
            er = abs(y(k)-yt)
            erm = max(erm,er)
            a1 = a1*alph1/i
         enddo
         a2 = a2*alph2/j
      enddo
      end subroutine edit2

   end program dvdemo
