      SUBROUTINE RKTOUT
C***********************************************************************
C SPECIAL OUTPUT FOR ROCKET PROBLEMS.
C***********************************************************************
      IMPLICIT NONE
      INCLUDE 'cea.inc'
C LOCAL VARIABLES
      CHARACTER*4 exit(11)
      CHARACTER*15 fi,fiv,fr,z(4)
      INTEGER i,i23,i46,i57,i68,i79,ione,ixfr,ixfz,j,k,line,ln,mae,mcf,
     &        misp,mivac,mmach,mppf,mppj,mxx(8),nex
      INTEGER INDEX
      REAL*8 agv,aw,gc,tem,tra,vaci(NCOL),ww
      SAVE agv,aw,fi,fiv,fr,gc,i,i23,i46,i57,i68,i79,ione,ixfr,ixfz,j,k,
     &  line,ln,mae,mcf,misp,mivac,mmach,mppf,mppj,mxx,nex,tem,tra,vaci,
     &  ww,z
C
      EQUIVALENCE (mxx(1),mppf)
      EQUIVALENCE (mxx(2),mppj)
      EQUIVALENCE (mxx(3),mmach)
      EQUIVALENCE (mxx(4),mae)
      EQUIVALENCE (mxx(5),mcf)
      EQUIVALENCE (mxx(6),mivac)
      EQUIVALENCE (mxx(7),misp)
      DATA exit/11*'EXIT'/
      IF ( .NOT.Eql ) THEN
        WRITE (IOOUT,99004)
        IF ( Nfz.GT.1 ) WRITE (IOOUT,99005) Nfz
      ELSE
        WRITE (IOOUT,99001)
        IF ( Iopt.NE.0 ) WRITE (IOOUT,99002)
        IF ( Iopt.EQ.0 ) WRITE (IOOUT,99003)
      ENDIF
      IF ( Ttt(1).EQ.T(It) ) WRITE (IOOUT,99006)
      tem = Ppp(1)*14.696006D0/1.01325D0
      WRITE (IOOUT,99009) 'Pin',tem
      i23 = 2
      IF ( Iopt.GT.0 ) THEN
        IF ( Iopt.EQ.1 ) WRITE (IOOUT,99007) Subar(1),App(2)
        IF ( Iopt.EQ.2 ) WRITE (IOOUT,99008) Ma,App(2)
        i23 = 3
      ENDIF
      CALL OUT1
      Fmt(4) = Fmt(6)
      nex = Npt - 2
      IF ( Page1 ) THEN
        ione = 0
        i46 = 4
        i57 = 5
        i68 = 6
        i79 = 7
      ELSE
        ione = i23
      ENDIF
C PRESSURE RATIOS
      IF ( Iopt.EQ.0 ) THEN
        WRITE (IOOUT,99011) (exit(i),i=1,nex)
        CALL VARFMT(App)
        WRITE (IOOUT,Fmt) 'Pinf/P         ',(App(j),j=1,Npt)
      ELSE
        nex = nex - 1
        WRITE (IOOUT,99010) (exit(i),i=1,nex)
        X(1) = 1.D0
        DO i = 2,Npt
          X(i) = Ppp(1)/Ppp(i)
        ENDDO
        CALL VARFMT(X)
        WRITE (IOOUT,Fmt) 'Pinj/P         ',(X(i),i=1,Npt)
      ENDIF
      CALL OUT2
      DO i = 1,8
        mxx(i) = 0
      ENDDO
      DO 100 i = 1,Nplt
        ixfz = INDEX(Pltvar(i)(2:),'fz')
        ixfr = INDEX(Pltvar(i)(2:),'fr')
        IF ( ixfz.NE.0.OR.ixfr.NE.0 ) THEN
          IF ( Eql ) GOTO 100
        ELSEIF ( .NOT.Eql ) THEN
          GOTO 100
        ENDIF
        IF ( Pltvar(i)(:4).EQ.'pi/p'.OR.Pltvar(i)(:3).EQ.'pip' ) THEN
          IF ( Iopt.EQ.0 ) mppf = i
          IF ( Iopt.NE.0 ) mppj = i
        ELSEIF ( Pltvar(i)(:4).EQ.'mach' ) THEN
          mmach = i
        ELSEIF ( Pltvar(i)(:2).EQ.'ae' ) THEN
          mae = i
        ELSEIF ( Pltvar(i)(:2).EQ.'cf' ) THEN
          mcf = i
        ELSEIF ( Pltvar(i)(:4).EQ.'ivac' ) THEN
          mivac = i
        ELSEIF ( Pltvar(i)(:3).EQ.'isp' ) THEN
          misp = i
        ENDIF
 100  CONTINUE
      IF ( Siunit ) THEN
        agv = 1.
        gc = 1.
        fr = 'CSTAR, M/SEC'
        fiv = 'Ivac, M/SEC'
        fi = 'Isp, M/SEC'
      ELSE
        gc = 32.174
        agv = 9.80665
        fr = 'CSTAR, FT/SEC'
        fiv = 'Ivac,LB-SEC/LB'
        fi = 'Isp, LB-SEC/LB'
      ENDIF
      DO k = 2,Npt
        Spim(k) = (2.*Rr*(Hsum(1)-Hsum(k)))**.5/agv
C AW IS THE LEFT SIDE OF EQ.(6.12) IN RP-1311,PT I.
        aw = Rr*Ttt(k)/(Ppp(k)*Wm(k)*Spim(k)*agv**2)
        IF ( k.EQ.i23 ) THEN
          IF ( Iopt.EQ.0 ) Cstr = gc*Ppp(1)*aw
          IF ( Iopt.NE.0 ) Cstr = gc*Ppp(1)/App(2)*aw
        ENDIF
        vaci(k) = Spim(k) + Ppp(k)*aw
        Vmoc(k) = 0.
        IF ( Sonvel(k).NE.0. ) Vmoc(k) = Spim(k)*agv/Sonvel(k)
      ENDDO
C MACH NUMBER
      Vmoc(1) = 0.
      IF ( Gammas(i23).EQ.0. ) Vmoc(i23) = 0.
      Fmt(7) = '3,'
      WRITE (IOOUT,Fmt) 'MACH NUMBER    ',(Vmoc(j),j=1,Npt)
      IF ( Trnspt ) CALL OUT4
      WRITE (IOOUT,99013)
C AREA RATIO
      Fmt(4) = '9x,'
      Fmt(i46) = '9x,'
      CALL VARFMT(Aeat)
      Fmt(5) = ' '
      Fmt(i57) = ' '
      WRITE (IOOUT,Fmt) 'Ae/At          ',(Aeat(j),j=2,Npt)
C C*
      Fmt(i57) = '13'
      Fmt(i68) = Fmt(i68+2)
      Fmt(i79) = '1,'
      WRITE (IOOUT,Fmt) fr,(Cstr,j=2,Npt)
C CF - THRUST COEFICIENT
      Fmt(i79) = '4,'
      DO i = 2,Npt
        X(i) = gc*Spim(i)/Cstr
      ENDDO
      WRITE (IOOUT,Fmt) 'CF             ',(X(j),j=2,Npt)
C VACUUM IMPULSE
      Fmt(i57) = '13'
      Fmt(i79) = '1,'
      WRITE (IOOUT,Fmt) fiv,(vaci(j),j=2,Npt)
C SPECIFIC IMPULSE
      WRITE (IOOUT,Fmt) fi,(Spim(j),j=2,Npt)
      IF ( Nplt.GT.0 ) THEN
        Spim(1) = 0
        Aeat(1) = 0
        Vmoc(1) = 0
        vaci(1) = 0
        X(1) = 0
        Spim(1) = 0
        DO i = ione + 1,Npt
          IF ( mppj.GT.0 ) Pltout(i+Iplt-ione,mppj) = Ppp(1)/Ppp(i)
          IF ( mppf.GT.0 ) Pltout(i+Iplt-ione,mppf) = App(i)
          IF ( mmach.GT.0 ) Pltout(i+Iplt-ione,mmach) = Vmoc(i)
          IF ( mae.GT.0 ) Pltout(i+Iplt-ione,mae) = Aeat(i)
          IF ( mcf.GT.0 ) Pltout(i+Iplt-ione,mcf) = X(i)
          IF ( mivac.GT.0 ) Pltout(i+Iplt-ione,mivac) = vaci(i)
          IF ( misp.GT.0 ) Pltout(i+Iplt-ione,misp) = Spim(i)
        ENDDO
      ENDIF
      WRITE (IOOUT,99012)
      Fmt(4) = ' '
      Fmt(5) = '13'
      Fmt(7) = '5,'
      IF ( Iopt.NE.0 ) THEN
        Fmt(i46) = Fmt(8)
        Fmt(i57) = Fmt(9)
      ENDIF
      IF ( .NOT.Eql ) THEN
        IF ( Massf ) THEN
          WRITE (IOOUT,99014) 'MASS'
        ELSE
          WRITE (IOOUT,99014) 'MOLE'
          ww = 1.D0/Totn(Nfz)
        ENDIF
C MOLE (OR MASS) FRACTIONS - FROZEN
        tra = 5.E-6
        IF ( Trace.NE.0. ) tra = Trace
        line = 0
        DO k = 1,Ngc
          IF ( Massf ) ww = Mw(k)
          X(line+1) = En(k,Nfz)*ww
          IF ( X(line+1).GE.tra ) THEN
            line = line + 1
            z(line) = Prod(k)
          ENDIF
          IF ( line.EQ.3.OR.k.EQ.Ngc ) THEN
            IF ( line.EQ.0 ) GOTO 200
            WRITE (IOOUT,99015) (z(ln),X(ln),ln=1,line)
            line = 0
          ENDIF
        ENDDO
      ENDIF
 200  CALL OUT3
      RETURN
99001 FORMAT (/////13x,' THEORETICAL ROCKET PERFORMANCE ASSUMING',
     &        ' EQUILIBRIUM')
99002 FORMAT (/11x,' COMPOSITION DURING EXPANSION FROM FINITE AREA',
     &        ' COMBUSTOR')
99003 FORMAT (/10x,' COMPOSITION DURING EXPANSION FROM INFINITE AREA',
     &        ' COMBUSTOR')
99004 FORMAT (/////10x,' THEORETICAL ROCKET PERFORMANCE ASSUMING FROZEN'
     &        ,' COMPOSITION')
99005 FORMAT (33X,'AFTER POINT',I2)
99006 FORMAT (25X,'AT AN ASSIGNED TEMPERATURE  ')
99007 FORMAT (' Ac/At =',F8.4,6x,'Pinj/Pinf =',F10.6)
99008 FORMAT (' MDOT/Ac =',F10.3,' (KG/S)/M**2',6x,'Pinj/Pinf =',F10.6)
99009 FORMAT (/1x,A3,' =',F8.1,' PSIA')
99010 FORMAT (/,17X,'INJECTOR  COMB END  THROAT',10(5X,A4))
99011 FORMAT (/17X,'CHAMBER   THROAT',11(5X,A4))
99012 FORMAT ()
99013 FORMAT (/' PERFORMANCE PARAMETERS'/)
99014 FORMAT (1x,A4,' FRACTIONS'/)
99015 FORMAT (1X,3(A15,F8.5,3X))
      END
