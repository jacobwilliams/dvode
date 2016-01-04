*DECK DVDEMO
C-----------------------------------------------------------------------
C Demonstration program for the DVODE package.
C This is the version of 30 April 2002.
C
C This version is in double precision, suitable for short wordlength
C computers, PCs, work-stations, etc.
C
C The package is used to solve two simple problems,
C one with a full Jacobian, the other with a banded Jacobian,
C with all 12 of the appropriate values of MF in each case.
C If the errors are too large, or other difficulty occurs,
C a warning message is printed.  All output is on unit LOUT = 6.
C-----------------------------------------------------------------------
C
      EXTERNAL F1, JAC1, F2, JAC2
      INTEGER I, IOPAR, IOPT, IOUT, IPAR, ISTATE, ITASK, ITOL, IWORK,
     1   JSV, LENIW, LENRW, LIW, LOUT, LRW, MBAND, METH, MF, MITER,
     2   ML, MU, NCFN, NEQ, NERR, NETF, NFE, NFEA, NJE, NLU, NNI, NOUT,
     3   NQU, NST
      DOUBLE PRECISION ATOL, DTOUT, ER, ERM, ERO, HU, RPAR, RTOL, RWORK,
     1   T, TOUT, TOUT1, Y
      DIMENSION Y(25), RWORK(847), IWORK(55)
      DATA LOUT/6/, TOUT1/1.39283880203D0/, DTOUT/2.214773875D0/
C
      NERR = 0
      ITOL = 1
      RTOL = 0.0D0
      ATOL = 1.0D-6
      LRW = 847
      LIW = 55
      IOPT = 0
C
C First problem
C
      NEQ = 2
      NOUT = 4
      WRITE (LOUT,110) NEQ,ITOL,RTOL,ATOL
 110  FORMAT(' Demonstration program for DVODE package'//
     1  ' Problem 1:   Van der Pol oscillator:'/
     2  '   xdotdot - 3*(1 - x**2)*xdot + x = 0,' ,
     3  '   x(0) = 2, xdot(0) = 0'/
     4  '   NEQ =',I2/
     5  '   ITOL =',I3,'   RTOL =',D10.1,'   ATOL =',D10.1)
C
      DO 195 JSV = 1,-1,-2
      DO 192 METH = 1,2
      DO 190 MITER = 0,3
      IF (JSV .LT. 0 .AND. MITER .EQ. 0) GO TO 190
      IF (JSV .LT. 0 .AND. MITER .EQ. 3) GO TO 190
      MF = JSV*(10*METH + MITER)
      WRITE (LOUT,120) MF
 120  FORMAT(//70('*')//' Solution with MF =',I4//
     1  6X,'t',15X,'x',15X,'xdot',7X,'NQ',6X,'H'/)
      T = 0.0D0
      Y(1) = 2.0D0
      Y(2) = 0.0D0
      ITASK = 1
      ISTATE = 1
      TOUT = TOUT1
      ERO = 0.0D0
      DO 170 IOUT = 1,NOUT
        CALL DVODE(F1,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,
     1     IOPT,RWORK,LRW,IWORK,LIW,JAC1,MF,RPAR,IPAR)
        HU = RWORK(11)
        NQU = IWORK(14)
        WRITE (LOUT,140) T,Y(1),Y(2),NQU,HU
 140    FORMAT(1X,D15.5,D16.5,D14.3,I5,D14.3)
        IF (ISTATE .LT. 0) GO TO 175
        IOPAR = IOUT - 2*(IOUT/2)
        IF (IOPAR .NE. 0) GO TO 170
        ER = ABS(Y(1))/ATOL
        ERO = MAX(ERO,ER)
        IF (ER .LT. 10000.0D0) GO TO 170
        WRITE (LOUT,150)
 150    FORMAT(//' Warning: Error exceeds 10000 * tolerance'//)
        NERR = NERR + 1
 170    TOUT = TOUT + DTOUT
 175  CONTINUE
      IF (ISTATE .LT. 0) NERR = NERR + 1
      NST = IWORK(11)
      NFE = IWORK(12)
      NJE = IWORK(13)
      NLU = IWORK(19)
      LENRW = IWORK(17)
      LENIW = IWORK(18)
      NNI = IWORK(20)
      NCFN = IWORK(21)
      NETF = IWORK(22)
      NFEA = NFE
      IF (MITER .EQ. 2) NFEA = NFE - NEQ*NJE
      IF (MITER .EQ. 3) NFEA = NFE - NJE
      WRITE (LOUT,180) LENRW,LENIW,NST,NFE,NFEA,NJE,NLU,NNI,NCFN,
     1                 NETF,ERO
 180  FORMAT(//' Final statistics for this run:'/
     1  ' RWORK size =',I4,'   IWORK size =',I4/
     2  ' Number of steps =',I5/
     3  ' Number of f-s   =',I5/
     4  ' (excluding J-s) =',I5/
     5  ' Number of J-s   =',I5/
     6  ' Number of LU-s  =',I5/
     7  ' Number of nonlinear iterations =', I5/
     8  ' Number of nonlinear convergence failures =', I5/
     9  ' Number of error test failures =', I5/
     1  ' Error overrun =',D10.2)
 190  CONTINUE
 192  CONTINUE
 195  CONTINUE
C
C Second problem
C
      NEQ = 25
      ML = 5
      MU = 0
      IWORK(1) = ML
      IWORK(2) = MU
      MBAND = ML + MU + 1
      NOUT = 5
      WRITE (LOUT,210) NEQ,ML,MU,ITOL,RTOL,ATOL
 210  FORMAT(//70('*')//' Demonstration program for DVODE package'//
     1  ' Problem 2: ydot = A * y , where',
     2  ' A is a banded lower triangular matrix'/
     3  '   derived from 2-D advection PDE'/
     4  '   NEQ =',I3,'   ML =',I2,'   MU =',I2/
     5  '   ITOL =',I3,'   RTOL =',D10.1,'   ATOL =',D10.1)
C
      DO 295 JSV = 1,-1,-2
      DO 292 METH = 1,2
      DO 290 MITER = 0,5
      IF (MITER .EQ. 1 .OR. MITER .EQ. 2) GO TO 290
      IF (JSV .LT. 0 .AND. MITER .EQ. 0) GO TO 290
      IF (JSV .LT. 0 .AND. MITER .EQ. 3) GO TO 290
      MF = JSV*(10*METH + MITER)
      WRITE (LOUT,220) MF
 220  FORMAT(//70('*')//' Solution with MF =',I4//
     1  6X,'t',13X,'Max.Err.',5X,'NQ',6X,'H'/)
      T = 0.0D0
      DO 230 I = 2,NEQ
 230    Y(I) = 0.0D0
      Y(1) = 1.0D0
      ITASK = 1
      ISTATE = 1
      TOUT = 0.01D0
      ERO = 0.0D0
      DO 270 IOUT = 1,NOUT
        CALL DVODE(F2,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,
     1     IOPT,RWORK,LRW,IWORK,LIW,JAC2,MF,RPAR,IPAR)
        CALL EDIT2(Y,T,ERM)
        HU = RWORK(11)
        NQU = IWORK(14)
        WRITE (LOUT,240) T,ERM,NQU,HU
 240    FORMAT(1X,D15.5,D14.3,I5,D14.3)
        IF (ISTATE .LT. 0) GO TO 275
        ER = ERM/ATOL
        ERO = MAX(ERO,ER)
        IF (ER .LE. 1000.0D0) GO TO 270
        WRITE (LOUT,150)
        NERR = NERR + 1
 270    TOUT = TOUT*10.0D0
 275  CONTINUE
      IF (ISTATE .LT. 0) NERR = NERR + 1
      NST = IWORK(11)
      NFE = IWORK(12)
      NJE = IWORK(13)
      NLU = IWORK(19)
      LENRW = IWORK(17)
      LENIW = IWORK(18)
      NNI = IWORK(20)
      NCFN = IWORK(21)
      NETF = IWORK(22)
      NFEA = NFE
      IF (MITER .EQ. 5) NFEA = NFE - MBAND*NJE
      IF (MITER .EQ. 3) NFEA = NFE - NJE
      WRITE (LOUT,180) LENRW,LENIW,NST,NFE,NFEA,NJE,NLU,NNI,NCFN,
     1                 NETF,ERO
 290  CONTINUE
 292  CONTINUE
 295  CONTINUE
C
      WRITE (LOUT,300) NERR
 300  FORMAT(//70('*')//' Number of errors encountered =',I3)
      STOP
      END
      SUBROUTINE F1 (NEQ, T, Y, YDOT, RPAR, IPAR)
      INTEGER NEQ, IPAR
      DOUBLE PRECISION T, Y, YDOT, RPAR
      DIMENSION Y(2), YDOT(2), RPAR(*), IPAR(*)
      YDOT(1) = Y(2)
      YDOT(2) = 3.0D0*(1.0D0 - Y(1)*Y(1))*Y(2) - Y(1)
      RETURN
      END
      SUBROUTINE JAC1 (NEQ, T, Y, ML, MU, PD, NROWPD, RPAR, IPAR)
      INTEGER NEQ, ML, MU, NROWPD, IPAR
      DOUBLE PRECISION T, Y, PD, RPAR
      DIMENSION Y(2), PD(NROWPD,2), RPAR(*), IPAR(*)
      PD(1,1) = 0.0D0
      PD(1,2) = 1.0D0
      PD(2,1) = -6.0D0*Y(1)*Y(2) - 1.0D0
      PD(2,2) = 3.0D0*(1.0D0 - Y(1)*Y(1))
      RETURN
      END
      SUBROUTINE F2 (NEQ, T, Y, YDOT, RPAR, IPAR)
      INTEGER NEQ, I, IPAR, J, K, NG
      DOUBLE PRECISION T, Y, YDOT, RPAR, ALPH1, ALPH2, D
      DIMENSION Y(NEQ), YDOT(NEQ), RPAR(*), IPAR(*)
      DATA ALPH1/1.0D0/, ALPH2/1.0D0/, NG/5/
      DO 10 J = 1,NG
      DO 10 I = 1,NG
        K = I + (J - 1)*NG
        D = -2.0D0*Y(K)
        IF (I .NE. 1) D = D + Y(K-1)*ALPH1
        IF (J .NE. 1) D = D + Y(K-NG)*ALPH2
 10     YDOT(K) = D
      RETURN
      END
      SUBROUTINE JAC2 (NEQ, T, Y, ML, MU, PD, NROWPD, RPAR, IPAR)
      INTEGER NEQ, ML, MU, NROWPD, IPAR, J, MBAND, MU1, MU2, NG
      DOUBLE PRECISION T, Y, PD, RPAR, ALPH1, ALPH2
      DIMENSION Y(NEQ), PD(NROWPD,NEQ), RPAR(*), IPAR(*)
      DATA ALPH1/1.0D0/, ALPH2/1.0D0/, NG/5/
      MBAND = ML + MU + 1
      MU1 = MU + 1
      MU2 = MU + 2
      DO 10 J = 1,NEQ
        PD(MU1,J) = -2.0D0
        PD(MU2,J) = ALPH1
 10     PD(MBAND,J) = ALPH2
      DO 20 J = NG,NEQ,NG
 20     PD(MU2,J) = 0.0D0
      RETURN
      END
      SUBROUTINE EDIT2 (Y, T, ERM)
      INTEGER I, J, K, NG
      DOUBLE PRECISION Y, T, ERM, ALPH1, ALPH2, A1, A2, ER, EX, YT
      DIMENSION Y(25)
      DATA ALPH1/1.0D0/, ALPH2/1.0D0/, NG/5/
      ERM = 0.0D0
      IF (T .EQ. 0.0D0) RETURN
      EX = 0.0D0
      IF (T .LE. 30.0D0) EX = EXP(-2.0D0*T)
      A2 = 1.0D0
      DO 60 J = 1,NG
        A1 = 1.0D0
        DO 50 I = 1,NG
          K = I + (J - 1)*NG
          YT = T**(I+J-2)*EX*A1*A2
          ER = ABS(Y(K)-YT)
          ERM = MAX(ERM,ER)
          A1 = A1*ALPH1/I
 50       CONTINUE
        A2 = A2*ALPH2/J
 60     CONTINUE
      RETURN
      END
