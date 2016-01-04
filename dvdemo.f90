!DECK DVDEMO
      PROGRAM DVDEMO
      USE DVODE_MODULE
      IMPLICIT NONE
!-----------------------------------------------------------------------
! Demonstration program for the DVODE package.
! This is the version of 30 April 2002.
!
! This version is in double precision, suitable for short wordlength
! computers, PCs, work-stations, etc.
!
! The package is used to solve two simple problems,
! one with a full Jacobian, the other with a banded Jacobian,
! with all 12 of the appropriate values of MF in each case.
! If the errors are too large, or other difficulty occurs,
! a warning message is printed.  All output is on unit LOUT = 6.
!-----------------------------------------------------------------------
!
      EXTERNAL F1 , JAC1 , F2 , JAC2
      INTEGER i , iopar , iopt , iout , ipar(1) , istate , itask ,      &
            & itol , iwork , jsv , leniw , lenrw , liw , lout , lrw ,   &
            & mband , meth , mf , miter , ml , mu , ncfn , neq , nerr , &
            & netf , nfe , nfea , nje , nlu , nni , nout , nqu , nst
      DOUBLE PRECISION atol , dtout , er , erm , ero , hu , rpar(1) ,   &
                     & rtol , rwork , t , tout , tout1 , y
      DIMENSION y(25) , rwork(847) , iwork(55)
      DATA lout/6/ , tout1/1.39283880203D0/ , dtout/2.214773875D0/
!
      nerr = 0
      itol = 1
      rtol = 0.0D0
      atol = 1.0D-6
      lrw = 847
      liw = 55
      iopt = 0
!
! First problem
!
      neq = 2
      nout = 4
      WRITE (lout,99001) neq , itol , rtol , atol
99001 FORMAT (' Demonstration program for DVODE package'//              &
             &' Problem 1:   Van der Pol oscillator:'/                  &
             &'   xdotdot - 3*(1 - x**2)*xdot + x = 0,',                &
             &'   x(0) = 2, xdot(0) = 0'/'   NEQ =',I2/'   ITOL =',I3,  &
             &'   RTOL =',D10.1,'   ATOL =',D10.1)
!
      DO jsv = 1 , -1 , -2
         DO meth = 1 , 2
            DO miter = 0 , 3
               IF ( jsv>=0 .OR. miter/=0 ) THEN
                  IF ( jsv>=0 .OR. miter/=3 ) THEN
                     mf = jsv*(10*meth+miter)
                     WRITE (lout,99002) mf
99002                FORMAT (//70('*')//' Solution with MF =',I4//6X,   &
                           & 't',15X,'x',15X,'xdot',7X,'NQ',6X,'H'/)
                     t = 0.0D0
                     y(1) = 2.0D0
                     y(2) = 0.0D0
                     itask = 1
                     istate = 1
                     tout = tout1
                     ero = 0.0D0
                     DO iout = 1 , nout
                        CALL DVODE(F1,neq,y,t,tout,itol,[rtol],[atol],  &
                                 & itask,istate,iopt,rwork,lrw,iwork,   &
                                 & liw,JAC1,mf,rpar,ipar)
                        hu = rwork(11)
                        nqu = iwork(14)
                        WRITE (lout,99003) t , y(1) , y(2) , nqu , hu
99003                   FORMAT (1X,D15.5,D16.5,D14.3,I5,D14.3)
                        IF ( istate<0 ) GOTO 2
                        iopar = iout - 2*(iout/2)
                        IF ( iopar==0 ) THEN
                           er = ABS(y(1))/atol
                           ero = MAX(ero,er)
                           IF ( er>=10000.0D0 ) THEN
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
!
      neq = 25
      ml = 5
      mu = 0
      iwork(1) = ml
      iwork(2) = mu
      mband = ml + mu + 1
      nout = 5
      WRITE (lout,99004) neq , ml , mu , itol , rtol , atol
99004 FORMAT (//70('*')//' Demonstration program for DVODE package'//   &
             &' Problem 2: ydot = A * y , where',                       &
             &' A is a banded lower triangular matrix'/                 &
             &'   derived from 2-D advection PDE'/'   NEQ =',I3,        &
             &'   ML =',I2,'   MU =',I2/'   ITOL =',I3,'   RTOL =',     &
            & D10.1,'   ATOL =',D10.1)
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
                        t = 0.0D0
                        DO i = 2 , neq
                           y(i) = 0.0D0
                        ENDDO
                        y(1) = 1.0D0
                        itask = 1
                        istate = 1
                        tout = 0.01D0
                        ero = 0.0D0
                        DO iout = 1 , nout
                           CALL DVODE(F2,neq,y,t,tout,itol,[rtol],      &
                                    & [atol],itask,istate,iopt,rwork,   &
                                    & lrw,iwork,liw,JAC2,mf,rpar,ipar)
                           CALL EDIT2(y,t,erm)
                           hu = rwork(11)
                           nqu = iwork(14)
                           WRITE (lout,99006) t , erm , nqu , hu
99006                      FORMAT (1X,D15.5,D14.3,I5,D14.3)
                           IF ( istate<0 ) GOTO 4
                           er = erm/atol
                           ero = MAX(ero,er)
                           IF ( er>1000.0D0 ) THEN
                              WRITE (lout,99008)
                              nerr = nerr + 1
                           ENDIF
                           tout = tout*10.0D0
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
      END PROGRAM DVDEMO

      SUBROUTINE F1(Neq,T,Y,Ydot,Rpar,Ipar)
      IMPLICIT NONE
      INTEGER Neq , Ipar
      DOUBLE PRECISION T , Y , Ydot , Rpar
      DIMENSION Y(2) , Ydot(2) , Rpar(*) , Ipar(*)
      Ydot(1) = Y(2)
      Ydot(2) = 3.0D0*(1.0D0-Y(1)*Y(1))*Y(2) - Y(1)
      END SUBROUTINE F1

      SUBROUTINE JAC1(Neq,T,Y,Ml,Mu,Pd,Nrowpd,Rpar,Ipar)
      IMPLICIT NONE
      INTEGER Neq , Ml , Mu , Nrowpd , Ipar
      DOUBLE PRECISION T , Y , Pd , Rpar
      DIMENSION Y(2) , Pd(Nrowpd,2) , Rpar(*) , Ipar(*)
      Pd(1,1) = 0.0D0
      Pd(1,2) = 1.0D0
      Pd(2,1) = -6.0D0*Y(1)*Y(2) - 1.0D0
      Pd(2,2) = 3.0D0*(1.0D0-Y(1)*Y(1))
      END SUBROUTINE JAC1

      SUBROUTINE F2(Neq,T,Y,Ydot,Rpar,Ipar)
      IMPLICIT NONE
      INTEGER Neq , i , Ipar , j , k , ng
      DOUBLE PRECISION T , Y , Ydot , Rpar , alph1 , alph2 , d
      DIMENSION Y(Neq) , Ydot(Neq) , Rpar(*) , Ipar(*)
      DATA alph1/1.0D0/ , alph2/1.0D0/ , ng/5/
      DO j = 1 , ng
         DO i = 1 , ng
            k = i + (j-1)*ng
            d = -2.0D0*Y(k)
            IF ( i/=1 ) d = d + Y(k-1)*alph1
            IF ( j/=1 ) d = d + Y(k-ng)*alph2
            Ydot(k) = d
         ENDDO
      ENDDO
      END SUBROUTINE F2

      SUBROUTINE JAC2(Neq,T,Y,Ml,Mu,Pd,Nrowpd,Rpar,Ipar)
      IMPLICIT NONE
      INTEGER Neq , Ml , Mu , Nrowpd , Ipar , j , mband , mu1 , mu2 , ng
      DOUBLE PRECISION T , Y , Pd , Rpar , alph1 , alph2
      DIMENSION Y(Neq) , Pd(Nrowpd,Neq) , Rpar(*) , Ipar(*)
      DATA alph1/1.0D0/ , alph2/1.0D0/ , ng/5/
      mband = Ml + Mu + 1
      mu1 = Mu + 1
      mu2 = Mu + 2
      DO j = 1 , Neq
         Pd(mu1,j) = -2.0D0
         Pd(mu2,j) = alph1
         Pd(mband,j) = alph2
      ENDDO
      DO j = ng , Neq , ng
         Pd(mu2,j) = 0.0D0
      ENDDO
      END SUBROUTINE JAC2

      SUBROUTINE EDIT2(Y,T,Erm)
      IMPLICIT NONE
      INTEGER i , j , k , ng
      DOUBLE PRECISION Y , T , Erm , alph1 , alph2 , a1 , a2 , er , ex ,&
                     & yt
      DIMENSION Y(25)
      DATA alph1/1.0D0/ , alph2/1.0D0/ , ng/5/
      Erm = 0.0D0
      IF ( T==0.0D0 ) RETURN
      ex = 0.0D0
      IF ( T<=30.0D0 ) ex = EXP(-2.0D0*T)
      a2 = 1.0D0
      DO j = 1 , ng
         a1 = 1.0D0
         DO i = 1 , ng
            k = i + (j-1)*ng
            yt = T**(i+j-2)*ex*a1*a2
            er = ABS(Y(k)-yt)
            Erm = MAX(Erm,er)
            a1 = a1*alph1/i
         ENDDO
         a2 = a2*alph2/j
      ENDDO
      END SUBROUTINE EDIT2
