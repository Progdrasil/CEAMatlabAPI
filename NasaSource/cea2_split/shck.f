      SUBROUTINE SHCK
C***********************************************************************
C PRIMARY ROUTINE FOR SHOCK PROBLEMS.
C***********************************************************************
      IMPLICIT NONE
      INCLUDE 'cea.inc'
C LOCAL VARIABLES
      CHARACTER*1 cr12,cr52
      INTEGER i,iof,it1,it2,itr,j,n
      LOGICAL refl,seql,srefl
      REAL*8 ax,axx,b2,cormax,gg,hs,m2m1(NCOL),mis(13),mu12rt,p1,p21,
     &       p21l,p2p1(NCOL),pmn,rho12,rho52,rrho(NCOL),sg(78),t1,t21,
     &       t21l,t2t1(NCOL),ttmax,u1u2(NCOL),uis(13),utwo(NCOL),uu,wmx,
     &       ww
      REAL*8 DABS,DEXP,DLOG,DMIN1
      SAVE ax,axx,b2,cormax,cr12,cr52,gg,hs,i,iof,it1,it2,itr,j,m2m1,
     &  mis,mu12rt,n,p1,p21,p21l,p2p1,pmn,refl,rho12,rho52,rrho,seql,sg,
     &  srefl,t1,t21,t21l,t2t1,ttmax,u1u2,uis,utwo,uu,wmx,ww
C
      IF ( Trace.EQ.0. ) Trace = 5.E-9
      Tp = .TRUE.
      Cpmix = 0.
      srefl = .FALSE.
      IF ( .NOT.Short ) THEN
        WRITE (IOOUT,99001)
        WRITE (IOOUT,99002) Incdeq,Refleq,Incdfz,Reflfz
      ENDIF
      IF ( Refleq.OR.Reflfz ) srefl = .TRUE.
      seql = Incdeq
      IF ( T(1).EQ.0. ) T(1) = Rtemp(1)
      DO i = 1,Nsk
        uis(i) = U1(i)
        mis(i) = Mach1(i)
        IF ( Mach1(i).EQ.0.0.AND.U1(i).EQ.0.0 ) GOTO 100
      ENDDO
 100  IF ( Nsk.GT.NCOL ) THEN
        WRITE (IOOUT,99003) NCOL
        Nsk = NCOL
      ENDIF
      IF ( .NOT.Short ) THEN
        WRITE (IOOUT,99004) (U1(i),i=1,Nsk)
        WRITE (IOOUT,99005) (Mach1(i),i=1,Nsk)
      ENDIF
      iof = 0
 200  iof = iof + 1
      Oxfl = Oxf(iof)
      CALL NEWOF
      Incdeq = seql
 300  refl = .FALSE.
      it2 = 2
      it1 = 1
      Pp = P(1)
      Tt = T(1)
      IF ( .NOT.Incdeq ) THEN
C FROZEN
        DO n = 1,Nsk
          Dlvtp(n) = 1.
          Dlvpt(n) = -1.
        ENDDO
      ENDIF
      DO Npt = 1,Nsk
        Ppp(Npt) = P(Npt)
        Ttt(Npt) = T(Npt)
        IF ( Npt.GT.1 ) THEN
          IF ( Ppp(Npt).EQ.0. ) Ppp(Npt) = Ppp(Npt-1)
          IF ( Ttt(Npt).EQ.0. ) Ttt(Npt) = Ttt(Npt-1)
          Ssum(Npt) = Ssum(Npt-1)
          Hsum(Npt) = Hsum(Npt-1)
          IF ( Ttt(Npt).EQ.Tt.AND.Ppp(Npt).EQ.Pp ) GOTO 350
        ENDIF
        Pp = Ppp(Npt)
        Tt = Ttt(Npt)
        IF ( Tt.GE.Tg(1)*.8D0 ) THEN
          CALL HCALC
          Hsum(Npt) = Hsub0
        ELSE
          WRITE (IOOUT,99016) Tt,Npt
          GOTO 1000
        ENDIF
 350    IF ( Cpmix.NE.0. ) Gamma1 = Cpmix/(Cpmix-1./Wmix)
        A1 = (Rr*Gamma1*Tt/Wmix)**.5
        IF ( U1(Npt).EQ.0. ) U1(Npt) = A1*Mach1(Npt)
        IF ( Mach1(Npt).EQ.0. ) Mach1(Npt) = U1(Npt)/A1
        Wm(Npt) = Wmix
        Cpr(Npt) = Cpmix
        Gammas(Npt) = Gamma1
        Vlm(Npt) = Rr*Tt/(Wmix*Pp)
      ENDDO
      Npt = Nsk
C OUTPUT--1ST CONDITION
      WRITE (IOOUT,99006)
      IF ( .NOT.Incdeq ) THEN
        WRITE (IOOUT,99008)
      ELSE
        WRITE (IOOUT,99007)
      ENDIF
      Eql = .FALSE.
      CALL OUT1
      WRITE (IOOUT,99009)
      Fmt(4) = '13'
      Fmt(5) = ' '
      Fmt(7) = '4,'
      WRITE (IOOUT,Fmt) 'MACH NUMBER1   ',(Mach1(j),j=1,Npt)
      Fmt(7) = '2,'
      WRITE (IOOUT,Fmt) 'U1, M/SEC      ',(U1(j),j=1,Npt)
      CALL OUT2
C BEGIN CALCULATIONS FOR 2ND CONDITION
      IF ( Incdeq ) Eql = .TRUE.
      Npt = 1
 400  Gamma1 = Gammas(Npt)
      uu = U1(Npt)
      wmx = Wm(Npt)
      p1 = Ppp(Npt)
      t1 = Ttt(Npt)
      hs = Hsum(Npt)
      IF ( refl ) uu = u1u2(Npt)
      mu12rt = wmx*uu**2/(Rr*t1)
      IF ( refl ) THEN
C REFLECTED--SUBSCRIPTS 2=1, 5=2, P52=P21
        t21 = 2.
        b2 = (-1.-mu12rt-t21)/2.
        p21 = -b2 + SQRT(b2**2-t21)
      ELSE
        p21 = (2.*Gamma1*Mach1(Npt)**2-Gamma1+1.)/(Gamma1+1.)
C THE FOLLOWING IMPROVED FORMULATION FOR THE INITIAL ESTIMATE FOR THE
C 2ND CONDITION WAS MADE AND TESTED BY S. GORDON 7/10/89.
        IF ( .NOT.Eql ) THEN
          t21 = p21*(2./Mach1(Npt)**2+Gamma1-1.)/(Gamma1+1.)
        ELSE
          Pp = p21*p1
          Tp = .FALSE.
          Hp = .TRUE.
          Hsub0 = hs + uu**2/(2.*Rr)
          CALL EQLBRM
          t21 = Ttt(Npt)/t1
          Hp = .FALSE.
          Tp = .TRUE.
        ENDIF
      ENDIF
      p21l = DLOG(p21)
      ttmax = 1.05*Tg(4)/t1
      t21 = DMIN1(t21,ttmax)
      t21l = DLOG(t21)
      itr = 1
 500  IF ( Shkdbg ) WRITE (IOOUT,99010) itr,it2,it1,p21,it2,it1,t21,
     &     rho52
      Tt = t21*t1
      Pp = p21*p1
      IF ( .NOT.Eql ) THEN
C FROZEN
        Tln = DLOG(Tt)
        IF ( .NOT.Incdeq ) THEN
          CALL HCALC
          IF ( Tt.EQ.0. ) GOTO 600
          Hsum(Npt) = Hsub0
          Cpr(Npt) = Cpmix
        ELSE
          CALL CPHS
          Cpr(Npt) = Cpsum
          Hsum(Npt) = 0.
          DO j = 1,Ng
            Hsum(Npt) = Hsum(Npt) + H0(j)*En(j,Npt)
          ENDDO
          Hsum(Npt) = Hsum(Npt)*Tt
        ENDIF
      ELSE
        CALL EQLBRM
        IF ( Tt.EQ.0. ) GOTO 800
      ENDIF
      rho12 = wmx*t21/(Wm(Npt)*p21)
      gg = rho12*mu12rt
      rho52 = 1./rho12
      IF ( refl ) gg = -mu12rt*rho52/(rho52-1.)**2
      G(1,1) = -gg*Dlvpt(Npt) - p21
      G(1,2) = -gg*Dlvtp(Npt)
      G(1,3) = p21 - 1. + gg - mu12rt
      IF ( refl ) G(1,3) = p21 - 1. + gg*(rho52-1.)
      gg = gg*t1/wmx
      IF ( .NOT.refl ) gg = gg*rho12
      G(2,1) = -gg*Dlvpt(Npt) + Tt*(Dlvtp(Npt)-1.)/Wm(Npt)
      G(2,2) = -gg*Dlvtp(Npt) - Tt*Cpr(Npt)
      gg = 1. - rho12**2
      IF ( refl ) gg = (rho52+1.)/(rho52-1.)
      G(2,3) = Hsum(Npt) - hs - uu**2*gg/(2.*Rr)
      X(3) = G(1,1)*G(2,2) - G(1,2)*G(2,1)
      X(1) = (G(1,3)*G(2,2)-G(2,3)*G(1,2))/X(3)
      X(2) = (G(1,1)*G(2,3)-G(2,1)*G(1,3))/X(3)
      IF ( Shkdbg ) THEN
        WRITE (IOOUT,99011) G(1,1),G(1,2),G(1,3)
        WRITE (IOOUT,99011) G(2,1),G(2,2),G(2,3)
        WRITE (IOOUT,99012) X(1),X(2)
        WRITE (IOOUT,99013) Hsum(Npt),hs,uu,uu*rho12
      ENDIF
      ax = DABS(X(1))
      axx = DABS(X(2))
      IF ( axx.GT.ax ) ax = axx
      IF ( ax.GE..00005 ) THEN
        cormax = .40546511
        IF ( itr.GT.4 ) cormax = .22314355
        IF ( itr.GT.12 ) cormax = .09531018
        IF ( itr.GT.20 ) cormax = .04879016
        ax = ax/cormax
        IF ( ax.GT.1. ) THEN
          X(1) = X(1)/ax
          X(2) = X(2)/ax
        ENDIF
        p21l = p21l + X(1)
        t21l = t21l + X(2)
        p21 = DEXP(p21l)
        t21 = DEXP(t21l)
        IF ( Shkdbg ) WRITE (IOOUT,99014) cormax,X(1),X(2)
        IF ( itr.NE.1.OR.t21.LT.ttmax ) THEN
          itr = itr + 1
          IF ( itr.LT.61 ) GOTO 500
          WRITE (IOOUT,99015) U1(Npt)
        ELSE
          Tt = 0.
          Npt = Npt - 1
          GOTO 700
        ENDIF
      ENDIF
C CONVERGED OR TOOK 60 ITERATIONS WITHOUT CONVERGING.
C STORE RESULTS.
 600  rrho(Npt) = rho52
      m2m1(Npt) = Wm(Npt)/wmx
      p2p1(Npt) = p21
      t2t1(Npt) = t21
      utwo(Npt) = uu*rho12
      u1u2(Npt) = uu - utwo(Npt)
      IF ( Tt.GE.Tg(1)*.8D0.AND.Tt.LE.Tg(4)*1.1D0 ) THEN
        IF ( .NOT.Eql ) THEN
C FROZEN
          Ppp(Npt) = Pp
          Ttt(Npt) = Tt
          Gammas(Npt) = Cpr(Npt)/(Cpr(Npt)-1./wmx)
          Vlm(Npt) = Rr*Tt/(wmx*Pp)
          IF ( Incdeq ) THEN
            Ssum(Npt) = 0.
            DO j = 1,Ngc
              pmn = Pp*wmx*En(j,Npt)
              IF ( En(j,Npt).GT.0. ) Ssum(Npt) = Ssum(Npt) + En(j,Npt)
     &             *(S(j)-DLOG(pmn))
            ENDDO
          ENDIF
        ENDIF
        GOTO 900
      ENDIF
 700  WRITE (IOOUT,99016) Tt,Npt
      Tt = 0.
 800  IF ( Npt.LT.1 ) GOTO 1000
      Nsk = Npt
 900  IF ( Trnspt ) CALL TRANP
      Isv = 0
      IF ( Npt.LT.Nsk ) Isv = Npt
      IF ( Npt.EQ.1 ) Isv = -1
      Npt = Npt + 1
      IF ( Eql ) CALL SETEN
      IF ( Npt.LE.Nsk ) GOTO 400
      Npt = Nsk
      IF ( refl ) THEN
        IF ( .NOT.Eql ) WRITE (IOOUT,99020)
        IF ( Eql ) WRITE (IOOUT,99021)
        cr12 = '2'
        cr52 = '5'
      ELSE
        IF ( .NOT.Eql ) WRITE (IOOUT,99018)
        IF ( Eql ) WRITE (IOOUT,99019)
        cr12 = '1'
        cr52 = '2'
      ENDIF
      Fmt(7) = '2,'
      WRITE (IOOUT,Fmt) 'U'//cr52//', M/SEC      ',(utwo(j),j=1,Npt)
      CALL OUT2
      IF ( Trnspt ) CALL OUT4
      WRITE (IOOUT,99017)
      Fmt(7) = '3,'
      WRITE (IOOUT,Fmt) 'P'//cr52//'/P'//cr12//'           ',
     &              (p2p1(j),j=1,Npt)
      WRITE (IOOUT,Fmt) 'T'//cr52//'/T'//cr12//'           ',
     &              (t2t1(j),j=1,Npt)
      Fmt(7) = '4,'
      WRITE (IOOUT,Fmt) 'M'//cr52//'/M'//cr12//'           ',
     &              (m2m1(j),j=1,Npt)
      WRITE (IOOUT,Fmt) 'RHO'//cr52//'/RHO'//cr12//'       ',
     &              (rrho(j),j=1,Npt)
      Fmt(7) = '2,'
      IF ( .NOT.refl ) WRITE (IOOUT,Fmt) 'V2, M/SEC      ',(u1u2(j),
     &               j=1,Npt)
      IF ( refl ) WRITE (IOOUT,Fmt) 'U5+V2,M/SEC    ',(u1u2(j),j=1,Npt)
      IF ( .NOT.Eql ) THEN
C WRITE FROZEN MOLE (OR MASS) FRACTIONS
        Fmt(7) = '5,'
        IF ( .NOT.Incdeq ) THEN
          IF ( Massf ) THEN
            WRITE (IOOUT,99022) 'MASS'
          ELSE
            WRITE (IOOUT,99022) 'MOLE'
            ww = wmx
          ENDIF
          DO n = 1,Nreac
            j = Jray(n)
            IF ( Massf ) ww = Mw(j)
            WRITE (IOOUT,99023) Prod(j),(En(j,i)*ww,i=1,Npt)
          ENDDO
        ELSE
          Eql = .TRUE.
          CALL OUT3
          Eql = .FALSE.
        ENDIF
      ELSE
        CALL OUT3
      ENDIF
      Iplt = MIN(Iplt+Npt,500)
      IF ( srefl ) THEN
        IF ( .NOT.refl ) THEN
          refl = .TRUE.
          it2 = 5
          it1 = 2
          Eql = .TRUE.
          IF ( Reflfz ) THEN
            Eql = .FALSE.
            IF ( Refleq ) THEN
              j = 0
              DO i = 1,Npt
                j = j + 1
                sg(j) = u1u2(i)
                j = j + 1
                sg(j) = Wm(i)
                j = j + 1
                sg(j) = Ppp(i)
                j = j + 1
                sg(j) = Ttt(i)
                j = j + 1
                sg(j) = Hsum(i)
                j = j + 1
                sg(j) = Gammas(i)
              ENDDO
            ENDIF
          ENDIF
          Npt = 1
          GOTO 400
        ELSEIF ( .NOT.Eql.AND.Refleq ) THEN
          j = 1
          DO i = 1,Npt
            u1u2(i) = sg(j)
            Wm(i) = sg(j+1)
            Ppp(i) = sg(j+2)
            Ttt(i) = sg(j+3)
            Hsum(i) = sg(j+4)
            Gammas(i) = sg(j+5)
            j = j + 6
          ENDDO
          Eql = .TRUE.
          Npt = 1
          GOTO 400
        ENDIF
      ENDIF
      IF ( Incdeq.AND.Incdfz ) THEN
        Incdeq = .FALSE.
        Eql = .FALSE.
        GOTO 300
      ELSEIF ( iof.GE.Nof ) THEN
        Tp = .FALSE.
        DO n = 1,Nreac
          Rtemp(n) = T(1)
        ENDDO
      ELSE
        DO i = 1,Nsk
          U1(i) = uis(i)
          Mach1(i) = mis(i)
        ENDDO
        GOTO 200
      ENDIF
 1000 RETURN
99001 FORMAT (/'   *** INPUT FOR SHOCK PROBLEMS ***')
99002 FORMAT (/' INCDEQ =',L2,'   REFLEQ =',L2,'   INCDFZ =',L2,
     &        '    REFLFZ =',L2)
99003 FORMAT (/' WARNING!!  ONLY ',I2,' u1 OR mach1 VALUES ALLOWED ',
     &        '(SHCK)')
99004 FORMAT (/1p,' U1 =   ',5E13.6,/(8X,5E13.6))
99005 FORMAT (/1p,' MACH1 =',5E13.6,/(8X,5E13.6))
99006 FORMAT (////25X,'SHOCK WAVE PARAMETERS ASSUMING')
99007 FORMAT (/,16X,' EQUILIBRIUM COMPOSITION FOR INCIDENT SHOCKED',
     &        ' CONDITIONS'//)
99008 FORMAT (/,17X,' FROZEN COMPOSITION FOR INCIDENT SHOCKED',
     &        ' CONDITI1ONS'//)
99009 FORMAT (/' INITIAL GAS (1)')
99010 FORMAT (/' ITR NO.=',I3,3X,'P',I1,'/P',I1,' =',F9.4,3X,'T',I1,
     &        '/T',I1,' =',F9.4,'   RHO2/RHO1 =',F9.6)
99011 FORMAT (/' G(I,J)  ',3E15.8)
99012 FORMAT (/' X       ',2E15.8)
99013 FORMAT (/' HSUM HS UU U2 ',4E15.8)
99014 FORMAT (/' MAX.COR.=',e13.6,' X(1)=',e13.6,' X(2)=',e13.6)
99015 FORMAT (/6x,' WARNING!!  NO CONVERGENCE FOR u1=',F8.1,
     &        /'  ANSWERS NOT RELIABLE, SOLUTION MAY NOT EXIST (SHCK)')
99016 FORMAT (/' TEMPERATURE=',E12.4,' IS OUT OF EXTENDED RANGE ',
     &        'FOR POINT',I5,' (SHCK)')
99017 FORMAT ()
99018 FORMAT (/' SHOCKED GAS (2)--INCIDENT--FROZEN')
99019 FORMAT (/' SHOCKED GAS (2)--INCIDENT--EQUILIBRIUM')
99020 FORMAT (/' SHOCKED GAS (5)--REFLECTED--FROZEN')
99021 FORMAT (/' SHOCKED GAS (5)--REFLECTED--EQUILIBRIUM')
99022 FORMAT (/1x,A4,' FRACTIONS'/)
99023 FORMAT (' ',A16,F8.5,12F9.5)
      END
