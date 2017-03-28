      SUBROUTINE CPHS
C***********************************************************************
C CALCULATES THERMODYNAMIC PROPERTIES FOR INDIVIDUAL SPECIES
C***********************************************************************
      IMPLICIT NONE
      INCLUDE 'cea.inc'
C LOCAL VARIABLES
      REAL*8 cx(7),hcx(7),scx(7)
      INTEGER i,ij,j,jj,k
      SAVE i,ij,j,jj,k,scx
C
      DATA cx/2*0.,1.D0,.5D0,.6666666666666667D0,.75D0,.8D0/
      DATA hcx(3)/1.D0/
      k = 1
      IF ( Tt.GT.Tg(2) ) k = 2
      IF ( Tt.GT.Tg(3) ) k = 3
      cx(2) = 1.D0/Tt
      cx(1) = cx(2)**2
      scx(3) = Tln
      scx(2) = -cx(2)
      hcx(2) = Tln*cx(2)
      hcx(1) = -cx(1)
      scx(1) = hcx(1)*.5D0
      DO i = 4,7
        hcx(i) = cx(i)*Tt
        scx(i) = cx(i-1)*Tt
      ENDDO
      DO j = 1,Ng
        H0(j) = 0.D0
        S(j) = 0.D0
      ENDDO
      DO i = 7,4, - 1
        DO j = 1,Ng
          S(j) = (S(j)+Coef(j,i,k))*scx(i)
          H0(j) = (H0(j)+Coef(j,i,k))*hcx(i)
        ENDDO
      ENDDO
      DO i = 1,3
        DO j = 1,Ng
          S(j) = S(j) + Coef(j,i,k)*scx(i)
          H0(j) = H0(j) + Coef(j,i,k)*hcx(i)
        ENDDO
      ENDDO
      DO j = 1,Ng
        S(j) = S(j) + Coef(j,9,k)
        H0(j) = H0(j) + Coef(j,8,k)*cx(2)
      ENDDO
      IF ( .NOT.Tp.OR.Convg ) THEN
        DO j = 1,Ng
          Cp(j) = 0.D0
        ENDDO
        DO i = 7,4, - 1
          DO j = 1,Ng
            Cp(j) = (Cp(j)+Coef(j,i,k))*Tt
          ENDDO
        ENDDO
        DO i = 1,3
          DO j = 1,Ng
            Cp(j) = Cp(j) + Coef(j,i,k)*cx(i)
          ENDDO
        ENDDO
      ENDIF
      IF ( Npr.NE.0.AND.k.NE.3.AND.Ng.NE.Ngc ) THEN
        DO ij = 1,Npr
          j = Jcond(ij)
          jj = Jcond(ij) - Ng
          Cp(j) = 0.D0
          H0(j) = 0.D0
          S(j) = 0.D0
          DO i = 7,4, - 1
            S(j) = (S(j)+Cft(jj,i))*scx(i)
            H0(j) = (H0(j)+Cft(jj,i))*hcx(i)
            Cp(j) = (Cp(j)+Cft(jj,i))*Tt
          ENDDO
          DO i = 1,3
            S(j) = S(j) + Cft(jj,i)*scx(i)
            H0(j) = H0(j) + Cft(jj,i)*hcx(i)
            Cp(j) = Cp(j) + Cft(jj,i)*cx(i)
          ENDDO
          S(j) = S(j) + Cft(jj,9)
          H0(j) = H0(j) + Cft(jj,8)*cx(2)
        ENDDO
      ENDIF
      GOTO 99999
      ENTRY ALLCON
      DO jj = 1,Nc
        j = jj + Ng
        Cp(j) = 0.D0
        H0(j) = 0.D0
        S(j) = 0.D0
        DO i = 7,4, - 1
          S(j) = (S(j)+Cft(jj,i))*scx(i)
          H0(j) = (H0(j)+Cft(jj,i))*hcx(i)
          Cp(j) = (Cp(j)+Cft(jj,i))*Tt
        ENDDO
        DO i = 1,3
          S(j) = S(j) + Cft(jj,i)*scx(i)
          H0(j) = H0(j) + Cft(jj,i)*hcx(i)
          Cp(j) = Cp(j) + Cft(jj,i)*cx(i)
        ENDDO
        S(j) = S(j) + Cft(jj,9)
        H0(j) = H0(j) + Cft(jj,8)*cx(2)
      ENDDO
99999 END
      SUBROUTINE DETON
C***********************************************************************
C CHAPMAN-JOUGUET DETONATIONS.
C***********************************************************************
      IMPLICIT NONE
      INCLUDE 'cea.inc'
C LOCAL VARIABLES
      CHARACTER*15 fdv,fg1,fh1,fhs1,fm1,fmm1,fpp1,frr1,ft1,ftt1
      CHARACTER*3 unit
      INTEGER i,ii,iof,itr,j,mdv,mgam,mh,mmach,mp,mson,mt,mxx(8)
      INTEGER INDEX
      REAL*8 a11,a12,a21,a22,alam,alfa,amm,b1,b2,cpl(NCOL),d,gam,
     &       gm1(NCOL),h1(NCOL),p1,pp1,pub(NCOL),rk,rr1,rrho(NCOL),t1,
     &       tem,tt1,tub(NCOL),ud,x1,x2
      SAVE a11,a12,a21,a22,alam,alfa,amm,b1,b2,cpl,d,gam,gm1,h1,i,ii,
     &  iof,itr,j,mdv,mgam,mh,mmach,mp,mson,mt,mxx,p1,pp1,pub,rk,rr1,
     &  rrho,t1,tem,tt1,tub,ud,unit,x1,x2
C
      EQUIVALENCE (mxx(1),mp)
      EQUIVALENCE (mxx(2),mt)
      EQUIVALENCE (mxx(3),mgam)
      EQUIVALENCE (mxx(4),mh)
      EQUIVALENCE (mxx(5),mdv)
      EQUIVALENCE (mxx(6),mson)
      EQUIVALENCE (mxx(7),mmach)
      DATA ft1/'T1, K'/,fh1/'H1, CAL/G'/,fhs1/'H1, KJ/KG'/,
     &     fm1/'M1, (1/n) '/,fg1/'GAMMA1'/,fpp1/'P/P1'/,ftt1/'T/T1'/,
     &     fmm1/'M/M1'/,frr1/'RHO/RHO1'/,fdv/'DET VEL,M/SEC'/
      iof = 0
      Eql = .TRUE.
      IF ( T(1).EQ.0. ) THEN
        T(1) = Rtemp(1)
        Nt = 1
      ENDIF
 100  Tt = T(1)
      iof = iof + 1
      Oxfl = Oxf(iof)
      CALL NEWOF
C BEGIN T LOOP.
      DO It = 1,Nt
        t1 = T(It)
C BEGIN P LOOP.
        DO Ip = 1,Np
          p1 = P(Ip)
          Tt = t1
          Pp = p1
          CALL HCALC
          IF ( Tt.EQ.0. ) RETURN
          IF ( Detdbg ) CALL OUT1
          h1(Npt) = Hsub0*R
          tub(Npt) = t1
          pub(Npt) = p1
          cpl(Npt) = Cpmix*R
          itr = 0
          Tt = 3800.
          pp1 = 15.
          Pp = pp1*p1
C CALCULATE ENTHALPY FOR INITIAL ESTIMATE OF T2(TT AFTER EQLBRM)
          Hsub0 = h1(Npt)/R + .75*t1*pp1/Wmix
          Tp = .FALSE.
          Hp = .TRUE.
          CALL EQLBRM
          Hsub0 = h1(Npt)/R
          Hp = .FALSE.
          IF ( Tt.NE.0. ) THEN
            gam = Gammas(Npt)
            tt1 = Tt/t1
            ii = 0
            tem = tt1 - .75*pp1/(Cpr(Npt)*Wmix)
            amm = Wm(Npt)/Wmix
            IF ( Detdbg ) WRITE (IOOUT,99001) Tt
C LOOP FOR IMPROVING T2/T1 AND P2/P1 INITIAL ESTIMATE.
            DO ii = 1,3
              alfa = amm/tt1
              pp1 = (1.+gam)*(1.+(1.-4.*gam*alfa/(1.+gam)**2)**.5)
     &              /(2.*gam*alfa)
              rk = pp1*alfa
              tt1 = tem + .5*pp1*gam*(rk*rk-1.)/(Wmix*Cpr(Npt)*rk)
              IF ( Detdbg ) WRITE (IOOUT,99002) ii,pp1,tt1
            ENDDO
            Tp = .TRUE.
            Tt = t1*tt1
            rr1 = pp1*amm/tt1
C BEGIN MAIN ITERATION LOOP.
 110        itr = itr + 1
            Pp = p1*pp1
            CALL EQLBRM
            IF ( Npt.EQ.0 ) GOTO 200
            IF ( Tt.NE.0. ) THEN
              gam = Gammas(Npt)
              amm = Wm(Npt)/Wmix
              rr1 = pp1*amm/tt1
              a11 = 1./pp1 + gam*rr1*Dlvpt(Npt)
              a12 = gam*rr1*Dlvtp(Npt)
              a21 = .5*gam*(rr1**2-1.-Dlvpt(Npt)*(1.+rr1**2))
     &              + Dlvtp(Npt) - 1.
              a22 = -.5*gam*Dlvtp(Npt)*(rr1**2+1.) - Wm(Npt)*Cpr(Npt)
              b1 = 1./pp1 - 1. + gam*(rr1-1.)
              b2 = Wm(Npt)*(Hsum(Npt)-h1(Npt)/R)
     &             /Tt - .5*gam*(rr1*rr1-1.)
              d = a11*a22 - a12*a21
              x1 = (a22*b1-a12*b2)/d
              x2 = (a11*b2-a21*b1)/d
              alam = 1.
              tem = x1
              IF ( tem.LT.0. ) tem = -tem
              IF ( x2.GT.tem ) tem = x2
              IF ( -x2.GT.tem ) tem = -x2
              IF ( tem.GT.0.4054652 ) alam = .4054652/tem
              pp1 = pp1*EXP(x1*alam)
              tt1 = tt1*EXP(x2*alam)
              Tt = t1*tt1
              ud = rr1*(Rr*gam*Tt/Wm(Npt))**.5
              IF ( Detdbg ) WRITE (IOOUT,99003) itr,pp1,tt1,rr1,x1,x2
C CONVERGENCE TEST
              IF ( itr.LT.8.AND.tem.GT.0.5E-04 ) GOTO 110
              IF ( itr.LT.8 ) THEN
                rrho(Npt) = rr1
                IF ( cpl(Npt).EQ.0. ) THEN
                  gm1(Npt) = 0.
                  Vmoc(Npt) = 0.
                ELSE
                  gm1(Npt) = cpl(Npt)/(cpl(Npt)-R/Wmix)
                  Vmoc(Npt) = ud/(Rr*gm1(Npt)*t1/Wmix)**.5
                ENDIF
              ELSE
                WRITE (IOOUT,99004)
                Npt = Npt - 1
                Tt = 0.
              ENDIF
              IF ( Trnspt ) CALL TRANP
              Isv = 0
              IF ( Ip.NE.Np.OR.It.NE.Nt.AND.Tt.NE.0. ) THEN
                Isv = Npt
                IF ( Npt.NE.NCOL ) GOTO 120
              ENDIF
            ENDIF
C OUTPUT
            WRITE (IOOUT,99005)
            CALL OUT1
C SET MXX ARRAY FOR PLOTTING PARAMETERS
            DO i = 1,8
              mxx(i) = 0
            ENDDO
            DO i = 1,Nplt
              IF ( INDEX(Pltvar(i)(2:),'1').NE.0 ) THEN
                IF ( Pltvar(i)(:3).EQ.'son' ) THEN
                  mson = i
                ELSEIF ( Pltvar(i)(:3).EQ.'gam' ) THEN
                  mgam = i
                ELSEIF ( Pltvar(i)(:1).EQ.'h' ) THEN
                  mh = i
                ELSEIF ( Pltvar(i)(:1).EQ.'t' ) THEN
                  mt = i
                ELSEIF ( Pltvar(i)(:1).EQ.'p' ) THEN
                  mp = i
                ENDIF
              ELSEIF ( INDEX(Pltvar(i),'vel').NE.0 ) THEN
                mdv = i
              ELSEIF ( INDEX(Pltvar(i),'mach').NE.0 ) THEN
                mmach = i
              ENDIF
            ENDDO
            WRITE (IOOUT,99006)
            Fmt(4) = '13'
            Fmt(5) = ' '
            Fmt(7) = '4,'
            DO i = 1,Npt
              IF ( Siunit ) THEN
                V(i) = pub(i)
                unit = 'BAR'
              ELSE
                V(i) = pub(i)/1.01325D0
                unit = 'ATM'
              ENDIF
              IF ( mp.GT.0 ) Pltout(i+Iplt,mp) = V(i)
            ENDDO
            WRITE (IOOUT,Fmt) 'P1, '//unit//'        ',(V(j),j=1,Npt)
            Fmt(7) = '2,'
            WRITE (IOOUT,Fmt) ft1,(tub(j),j=1,Npt)
            IF ( .NOT.Siunit ) WRITE (IOOUT,Fmt) fh1,(h1(j),j=1,Npt)
            IF ( Siunit ) WRITE (IOOUT,Fmt) fhs1,(h1(j),j=1,Npt)
            DO i = 1,Npt
              V(i) = Wmix
              Sonvel(i) = (Rr*gm1(i)*tub(i)/Wmix)**.5
            ENDDO
            Fmt(7) = '3,'
            WRITE (IOOUT,Fmt) fm1,(V(j),j=1,Npt)
            Fmt(7) = '4,'
            WRITE (IOOUT,Fmt) fg1,(gm1(j),j=1,Npt)
            Fmt(7) = '1,'
            WRITE (IOOUT,Fmt) 'SON VEL1,M/SEC ',(Sonvel(j),j=1,Npt)
            IF ( Nplt.GT.0 ) THEN
              DO i = 1,Npt
                IF ( mt.GT.0 ) Pltout(i+Iplt,mt) = tub(i)
                IF ( mgam.GT.0 ) Pltout(i+Iplt,mgam) = gm1(i)
                IF ( mh.GT.0 ) Pltout(i+Iplt,mh) = h1(i)
                IF ( mson.GT.0 ) Pltout(i+Iplt,mson) = Sonvel(i)
              ENDDO
            ENDIF
            WRITE (IOOUT,99007)
            Fmt(4) = Fmt(6)
            CALL OUT2
            IF ( Trnspt ) CALL OUT4
            WRITE (IOOUT,99008)
            Fmt(7) = '3,'
            DO i = 1,Npt
              V(i) = Ppp(i)/pub(i)
              Pcp(i) = Ttt(i)/tub(i)
              Sonvel(i) = Sonvel(i)*rrho(i)
              IF ( mmach.GT.0 ) Pltout(i+Iplt,mmach) = Vmoc(i)
              IF ( mdv.GT.0 ) Pltout(i+Iplt,mdv) = Sonvel(i)
            ENDDO
            WRITE (IOOUT,Fmt) fpp1,(V(j),j=1,Npt)
            WRITE (IOOUT,Fmt) ftt1,(Pcp(j),j=1,Npt)
            DO i = 1,Npt
              V(i) = Wm(i)/Wmix
            ENDDO
            Fmt(7) = '4,'
            WRITE (IOOUT,Fmt) fmm1,(V(j),j=1,Npt)
            WRITE (IOOUT,Fmt) frr1,(rrho(j),j=1,Npt)
            WRITE (IOOUT,Fmt) 'DET MACH NUMBER',(Vmoc(j),j=1,Npt)
            Fmt(7) = '1,'
            WRITE (IOOUT,Fmt) fdv,(Sonvel(j),j=1,Npt)
            Eql = .TRUE.
            CALL OUT3
            Iplt = MIN(Iplt+Npt,500)
            IF ( Isv.EQ.0.AND.iof.EQ.Nof ) GOTO 200
            IF ( Np.EQ.1.AND.Nt.EQ.1 ) GOTO 100
            WRITE (IOOUT,99009)
            Npt = 0
 120        Npt = Npt + 1
            IF ( Isv.EQ.1 ) Isv = -1
            CALL SETEN
          ENDIF
        ENDDO
      ENDDO
c     Iplt = MIN(Iplt+Npt,500)
      Iplt = MIN(Iplt+Npt-1,500)
      IF ( iof.LT.Nof ) GOTO 100
 200  Tp = .FALSE.
      RETURN
99001 FORMAT (/' T EST.=',F8.2/11X,'P/P1',17X,'T/T1')
99002 FORMAT (I5,2E20.8)
99003 FORMAT (/' ITER =',I2,5X,'P/P1 =',E15.8,/7X,'T/T1 =',E15.8,5X,
     &        'RHO/RHO1 =',E15.8,/7X,'DEL LN P/P1 =',E15.8,5X,
     &        'DEL LN T/T1 =',E15.8)
99004 FORMAT (/
     &        ' CONSERVATION EQNS NOT SATISFIED IN 8 ITERATIONS (DETON)'
     &        )
99005 FORMAT (//,21X,'DETONATION PROPERTIES OF AN IDEAL REACTING GAS')
99006 FORMAT (/' UNBURNED GAS'/)
99007 FORMAT (/' BURNED GAS'/)
99008 FORMAT (/' DETONATION PARAMETERS'/)
99009 FORMAT (///)
      END
