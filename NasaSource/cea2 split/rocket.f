      SUBROUTINE ROCKET
C***********************************************************************
C EXECUTIVE ROUTINE FOR ROCKET PROBLEMS.
C***********************************************************************
      IMPLICIT NONE
      INCLUDE 'cea.inc'
C LOCAL VARIABLES
      INTEGER i,i01,i12,iof,iplt1,iplte,ipp,isub,isup1,isupsv,itnum,
     &        itrot,nar,nipp,niter,nn,npr1,nptth
      LOGICAL done,seql,thi
      REAL*8 a1l,acatsv,aeatl,appl,aratio,asq,b1,c1,check,cprf,dd,dh,
     &       dlnp,dlnpe,dlt,dp,eln,mat,msq,p1,pa,pcpa,pcplt,pinf,pinj,
     &       pinjas,pjrat,ppa,pr,pracat,prat,pratsv,pvg,test,tmelt,usq
      REAL*8 DABS,DLOG,DMAX1,DSQRT
      SAVE acatsv,aeatl,appl,aratio,asq,check,cprf,dd,dh,dlnp,dlnpe,dlt,
     &  done,dp,eln,i,i01,i12,iof,iplt1,iplte,ipp,isub,isup1,isupsv,
     &  itnum,itrot,mat,msq,nar,nipp,niter,nn,npr1,nptth,p1,pcpa,pcplt,
     &  pinf,pinj,pinjas,pjrat,ppa,pr,pracat,prat,pratsv,pvg,seql,test,
     &  thi,tmelt,usq
C
      DATA a1l/ - 1.26505/,b1/1.0257/,c1/ - 1.2318/,pa/1.E05/
      iplte = Iplt
      isup1 = 1
      App(1) = 1.
      Iopt = 0
      Npp = Npp + 2
      nn = Npp
      i01 = 0
      i12 = 1
      nipp = 1
      nptth = 2
      IF ( Fac ) THEN
        Eql = .TRUE.
        Npp = Npp + 1
        IF ( Acat.NE.0. ) THEN
          Iopt = 1
        ELSEIF ( Ma.NE.0. ) THEN
          Iopt = 2
        ELSE
          WRITE (IOOUT,99001)
          Tt = 0.
          GOTO 1400
        ENDIF
        i01 = 1
        i12 = 2
        nipp = 2
        nptth = 3
        DO i = Nsub,1, - 1
          Subar(i+1) = Subar(i)
        ENDDO
        Nsub = Nsub + 1
        IF ( Iopt.NE.1 ) THEN
          IF ( Acat.EQ.0. ) Acat = 2.
        ENDIF
        Subar(1) = Acat
      ELSEIF ( .NOT.Eql.AND.Nfz.GT.1.AND.Nsub.GT.0 ) THEN
        Nsub = 0
        WRITE (IOOUT,99023)
      ENDIF
      nn = nn + Nsub + Nsup
      IF ( Nfz.GT.2.AND.nn.GT.NCOL-2 ) THEN
        WRITE (IOOUT,99002) NCOL - 2
        Nfz = 1
        Froz = .FALSE.
      ENDIF
      seql = Eql
      iof = 0
      Tt = Tcest
      Pp = P(1)
      App(i12) = 1.
C LOOP FOR EACH O/F
 100  It = 1
      iof = iof + 1
      Oxfl = Oxf(iof)
      IF ( T(1).NE.0. ) THEN
        Tp = .TRUE.
      ELSE
        Hp = .TRUE.
      ENDIF
      Sp = .FALSE.
      CALL NEWOF
      IF ( T(1).NE.0. ) Tt = T(1)
C LOOP FOR CHAMBER PRESSURES
 200  DO Ip = 1,Np
        itnum = 0
        Area = .FALSE.
        IF ( T(1).EQ.0. ) Hp = .TRUE.
        IF ( T(1).NE.0. ) Tp = .TRUE.
        Sp = .FALSE.
        Eql = .TRUE.
        isub = 1
        Isup = 1
        Pp = P(Ip)
        pinf = Pp
        ipp = 1
        itrot = 3
        isupsv = 1
        niter = 1
        Page1 = .TRUE.
        iplt1 = iplte
        Iplt = iplte
        done = .FALSE.
C LOOP FOR OUTPUT COLUMNS
 250    nar = Npt
        IF ( Eql ) THEN
          CALL EQLBRM
          IF ( Npt.EQ.Nfz ) cprf = Cpsum
        ELSE
          CALL FROZEN
        ENDIF
C TT = 0 IF NO CONVERGENCE
        IF ( Tt.NE.0. ) THEN
C TEST FOR FINITE AREA COMBUSTOR
          IF ( .NOT.Fac ) GOTO 400
          pinjas = P(Ip)*pa
          pinj = pinjas
          IF ( Npt.LE.2 ) THEN
            IF ( Npt.EQ.1.AND.Trnspt ) CALL TRANP
            IF ( Npt.EQ.2 ) pinf = Ppp(2)
          ENDIF
          IF ( Npt.NE.1 ) GOTO 400
C INITIAL ESTIMATE FOR PC (AND ACAT IF NOT ASSIGNED)
          DO i = 1,4
            prat = (b1+c1*Acat)/(1.+a1l*Acat)
            ppa = pinj*prat
            IF ( Iopt.EQ.1 ) GOTO 260
            Acat = ppa/(Ma*2350.)
            IF ( Acat.GE.1. ) THEN
              pratsv = prat
              IF ( Debugf ) THEN
                IF ( i.LE.1 ) WRITE (IOOUT,99004)
                WRITE (IOOUT,99005) i,ppa,Acat
              ENDIF
            ELSE
              WRITE (IOOUT,99003) Ma
              Tt = 0.
              GOTO 1400
            ENDIF
          ENDDO
          Subar(1) = Acat
 260      Pp = ppa/pa
          App(1) = Pp/Ppp(1)
          GOTO 1100
        ELSE
          IF ( Npt.LT.1 ) GOTO 1400
          IF ( .NOT.Area ) GOTO 600
          Npt = nar - 1
          Isup = Nsup + 2
          Isv = 0
          itnum = 0
          GOTO 950
        ENDIF
 300    Hp = .TRUE.
        Sp = .FALSE.
        niter = niter + 1
        Isv = 0
        Npt = 2
        ipp = 2
        CALL SETEN
        GOTO 250
 350    done = .TRUE.
        App(1) = Ppp(2)/Ppp(1)
        Area = .FALSE.
        IF ( Nsub.GT.1 ) isub = 2
        Isv = 4
        Npt = 2
        ipp = MIN(4,Npp)
        CALL SETEN
        Cpr(2) = Cpr(4)
        Dlvpt(2) = Dlvpt(4)
        Dlvtp(2) = Dlvtp(4)
        Gammas(2) = Gammas(4)
        Hsum(2) = Hsum(4)
        Ppp(2) = Ppp(4)
        App(2) = Ppp(1)/pinf
        Ssum(2) = Ssum(4)
        Totn(2) = Totn(4)
        Ttt(2) = Ttt(4)
        Vlm(2) = Vlm(4)
        Wm(2) = Wm(4)
        IF ( .NOT.Short ) WRITE (IOOUT,99009)
        GOTO 600
C INITIALIZE FOR THROAT
 400    IF ( ipp.GT.nipp ) THEN
          usq = 2.*(Hsum(1)-Hsum(Npt))*Rr
          IF ( ipp.GT.nptth ) GOTO 600
C THROAT
          IF ( .NOT.thi ) THEN
            Vv = Vlm(nptth)
            pvg = Pp*Vv*Gammas(nptth)
            IF ( pvg.EQ.0. ) THEN
              WRITE (IOOUT,99010)
              GOTO 550
            ELSE
              msq = usq/pvg
              IF ( Debug(1).OR.Debug(2) ) WRITE (IOOUT,99011) usq,pvg
              dh = DABS(msq-1.D0)
              IF ( dh.LE.0.4D-4 ) GOTO 550
              IF ( itrot.GT.0 ) THEN
                p1 = Pp
                IF ( Jsol.NE.0 ) THEN
                  tmelt = Tt
                  Pp = Pp*(1.D0+msq*Gammas(nptth))/(Gammas(nptth)+1.D0)
                ELSEIF ( tmelt.EQ.0. ) THEN
                  Pp = Pp*(1.D0+msq*Gammas(nptth))/(Gammas(nptth)+1.D0)
                ELSE
                  WRITE (IOOUT,99012)
                  dlt = DLOG(tmelt/Tt)
                  dd = dlt*Cpr(nptth)/(Enn*Dlvtp(nptth))
                  Pp = Pp*EXP(dd)
                  App(nptth) = P(Ip)/Pp
                  IF ( Fac ) App(nptth) = pinf/Pp
                  IF ( Eql.AND..NOT.Short ) WRITE (IOOUT,99013) 
     &                            App(nptth)
                  thi = .TRUE.
                  GOTO 250
                ENDIF
                GOTO 500
              ELSEIF ( itrot.LT.0 ) THEN
                IF ( itrot.LT.-19 ) THEN
                  WRITE (IOOUT,99010)
                  GOTO 550
                ELSE
                  IF ( Npr.NE.npr1 ) GOTO 550
                  Pp = Pp - dp
                  GOTO 500
                ENDIF
              ELSEIF ( Npr.EQ.npr1 ) THEN
                WRITE (IOOUT,99010)
                GOTO 550
              ELSE
                dp = DABS(Pp-p1)/20.
                Pp = DMAX1(Pp,p1)
                WRITE (IOOUT,99012)
                Pp = Pp - dp
                GOTO 500
              ENDIF
            ENDIF
          ELSE
            Gammas(nptth) = 0.
            GOTO 550
          ENDIF
        ELSE
          IF ( .NOT.Fac.AND.Trnspt ) CALL TRANP
          IF ( Npt.EQ.Nfz ) Eql = seql
          Tp = .FALSE.
          Hp = .FALSE.
          Sp = .TRUE.
          S0 = Ssum(i12)
        ENDIF
 450    tmelt = 0.
        itrot = 3
        thi = .FALSE.
        App(nptth) = ((Gammas(i12)+1.)/2.)
     &               **(Gammas(i12)/(Gammas(i12)-1.))
        IF ( Eql.AND..NOT.Short ) WRITE (IOOUT,99013) App(nptth)
        Pp = pinf/App(nptth)
        Isv = -i12
        GOTO 1200
 500    npr1 = Npr
        App(nptth) = P(Ip)/Pp
        IF ( Fac ) App(nptth) = pinf/Pp
        IF ( Eql.AND..NOT.Short ) WRITE (IOOUT,99013) App(nptth)
        itrot = itrot - 1
        GOTO 250
 550    Awt = Enn*Tt/(Pp*usq**.5)
        pcplt = DLOG(App(nptth))
 600    Isv = 0
        Aeat(Npt) = Enn*Ttt(Npt)/(Pp*usq**.5*Awt)
        IF ( Tt.EQ.0. ) GOTO 1150
        IF ( Area ) GOTO 750
        IF ( Trnspt.AND.(.NOT.Fac.OR.done.OR.Npt.GT.2) ) CALL TRANP
        IF ( Npt.EQ.Nfz ) Eql = seql
        IF ( Fac ) THEN
          IF ( Npt.EQ.nptth ) THEN
            Area = .TRUE.
            GOTO 750
          ELSEIF ( Npt.EQ.2.AND.done ) THEN
            Npt = 3
C  The following statement was corrected 1/30/2004.  Only fac parameters 
C    after combustion were affected--generally extra or missing points.
C  (remove) IF ( ipp.LE.Npp ) ipp = ipp - 1
            IF ( ipp.LT.Npp.OR.npp.EQ.4 ) ipp = ipp - 1
          ENDIF
        ENDIF
 650    IF ( ipp.LT.Npp ) GOTO 1100
 700    IF ( Nsub.EQ.i01.AND.Nsup.EQ.0 ) GOTO 1150
        Area = .TRUE.
C PCP ESTIMATES FOR AREA RATIOS
 750    IF ( itnum.EQ.0 ) THEN
          dlnp = 1.
          itnum = 1
          aratio = Subar(isub)
          IF ( (.NOT.Fac.OR.done).AND.Nsub.LE.i01 ) aratio = Supar(Isup)
          IF ( .NOT.Eql.AND.Nfz.GE.3 ) THEN
            IF ( aratio.LE.Aeat(Nfz) ) THEN
              WRITE (IOOUT,99014) Nfz
              GOTO 1050
            ENDIF
          ENDIF
          IF (aratio .LT. 1.d0 ) THEN
            WRITE (IOOUT,99025) 
            GOTO 1050
          ENDIF
          eln = DLOG(aratio)
          IF ( Fac ) THEN
            IF ( .NOT.done ) GOTO 800
          ENDIF
          IF ( Nsub.LE.i01 ) THEN
            IF ( Nfz.EQ.ipp ) isupsv = Isup
            IF ( Supar(Isup).LT.2. ) THEN
              appl = DSQRT(eln*(1.535d0+3.294d0*eln)) + pcplt
              GOTO 1100
            ELSE
              IF ( Isup.GT.isup1.AND.Supar(Isup-1).GE.2. ) GOTO 850
              appl = Gammas(nptth) + eln*1.4
              GOTO 1100
            ENDIF
          ENDIF
C TEST FOR CONVERGENCE ON AREA RATIO.
        ELSEIF ( Gammas(Npt).GT.0. ) THEN
          check = .00004
          IF ( Debug(Npt) ) WRITE (IOOUT,99016) itnum,aratio,Aeat(Npt),
     &                             App(Npt),dlnp
          IF ( DABS(Aeat(Npt)-aratio)/aratio.LE.check ) GOTO 900
          IF ( ABS(dlnp).LT..00004 ) GOTO 900
          aeatl = DLOG(Aeat(Npt))
          itnum = itnum + 1
          IF ( itnum.GT.10 ) THEN
            WRITE (IOOUT,99017) aratio
            GOTO 900
          ELSE
C IMPROVED PCP ESTIMATES.
            asq = Gammas(Npt)*Enn*Rr*Tt
            dlnpe = Gammas(Npt)*usq/(usq-asq)
            GOTO 850
          ENDIF
        ELSE
          WRITE (IOOUT,99015)
          Npt = Npt - 1
          IF ( Nsub.LE.0 ) isup1 = 100
          IF ( Nsub.LT.0. ) Nsup = Isup - 1
          IF ( Nsub.GT.0 ) Nsub = isub - 1
          GOTO 1000
        ENDIF
 800    appl = pcplt/(Subar(isub)+(10.587*eln**2+9.454)*eln)
        IF ( aratio.LT.1.09 ) appl = .9*appl
        IF ( aratio.GT.10. ) appl = appl/aratio
        IF ( isub.GT.1.OR.Npt.EQ.NCOL ) GOTO 1100
        GOTO 1200
 850    dlnp = dlnpe*eln - dlnpe*aeatl
        appl = appl + dlnp
        IF ( itnum.EQ.1 ) GOTO 1100
        IF ( appl.LT.0. ) appl = .000001
        App(Npt) = EXP(appl)
        Pp = pinf/App(Npt)
        GOTO 250
C CONVERGENCE HAS BEEN REACHED FOR ASSIGNED AREA RATIO
 900    Aeat(Npt) = aratio
        IF ( Fac ) THEN
          IF ( .NOT.done ) THEN
            IF ( Iopt.EQ.1 ) THEN
C OPTION 1 FOR FINITE AREA COMBUSTOR. INPUT IS ASSIGNED INJECTOR
C PRESSURE AND CONTRACTION RATIO. IMPROVED ESTIMATE FOR PC
              Area = .FALSE.
              itnum = 0
              ppa = Ppp(Npt)*pa
              pinj = ppa + 1.D05*usq/Vlm(Npt)
              test = (pinj-pinjas)/pinjas
              pcpa = pinf*pa
              IF ( Debugf ) THEN
                WRITE (IOOUT,99006)
                WRITE (IOOUT,99007) niter,test,pinjas,pinj,pcpa,ppa,
     &                          acatsv,Acat
              ENDIF
              IF ( ABS(test).LT.0.00002 ) GOTO 350
              prat = pinjas/pinj
              Pp = pinf*prat
              GOTO 300
            ELSEIF ( Iopt.EQ.2 ) THEN
C OPTION 2 FOR FINITE AREA COMBUSTOR. INPUT IS ASSIGNED INJECTOR
C PRESSURE AND MASS FLOW PER UNIT AREA. IMPROVED ESTIMATE FOR PC
C AND ACAT
              acatsv = Acat
              pratsv = prat
              Area = .FALSE.
              itnum = 0
              ppa = Ppp(4)*pa
              pinj = ppa + 1.D05*usq/Vlm(4)
              mat = pa/(Awt*Rr)
              Acat = mat/Ma
              prat = (b1+c1*Acat)/(1.+a1l*Acat)
              test = (pinj-pinjas)/pinjas
              pcpa = pinf*pa
              IF ( Debugf ) THEN
                WRITE (IOOUT,99006)
                WRITE (IOOUT,99007) niter,test,pinjas,pinj,pcpa,ppa,
     &                          acatsv,Acat
              ENDIF
              IF ( ABS(test).LT.0.00002 ) GOTO 350
              pjrat = pinj/pinjas
              Pp = pinf
              DO i = 1,2
                pracat = pratsv/prat
                pr = pjrat*pracat
                Pp = Pp/pr
                pcpa = Pp*pa
                Acat = Acat/pr
                Subar(1) = Acat
                pratsv = prat
                pjrat = 1.
                prat = (b1+c1*Acat)/(1.+a1l*Acat)
                IF ( Debugf ) WRITE (IOOUT,99008) pcpa,Acat,pjrat,pracat
              ENDDO
              GOTO 300
            ENDIF
          ENDIF
        ENDIF
 950    IF ( Trnspt ) CALL TRANP
        IF ( Npt.EQ.Nfz ) Eql = seql
 1000   itnum = 0
        IF ( Nsub.GT.i01 ) THEN
          isub = isub + 1
          IF ( isub.LE.Nsub ) GOTO 750
          isub = 1
          Nsub = -Nsub
          IF ( Isup.LE.Nsup ) GOTO 750
          Area = .FALSE.
          GOTO 1150
        ENDIF
 1050   Isup = Isup + 1
        itnum = 0
        IF ( Isup.LE.Nsup ) GOTO 750
        Isup = isupsv
        Area = .FALSE.
        GOTO 1150
C TEST FOR OUTPUT -- SCHEDULES COMPLETE OR NPT=NCOL
 1100   Isv = Npt
        IF ( Npt.NE.NCOL ) GOTO 1200
 1150   IF ( .NOT.Eql ) THEN
          IF ( Nfz.LE.1 ) THEN
            Cpr(Nfz) = cprf
            Gammas(Nfz) = cprf/(cprf-1./Wm(Nfz))
          ENDIF
        ENDIF
        CALL RKTOUT
        Iplt = Iplt + Npt
        IF ( .NOT.Page1 ) THEN
          Iplt = Iplt - 2
          IF ( Iopt.NE.0 ) Iplt = Iplt - 1
          Iplt = MIN(Iplt,500)
        ELSE
          Page1 = .FALSE.
        ENDIF
        iplte = MAX(iplte,Iplt)
        dlnp = 1.
        IF ( Tt.EQ.0. ) Area = .FALSE.
        IF ( .NOT.Eql.AND.Tt.EQ.0. ) WRITE (IOOUT,99018)
        IF ( Isv.EQ.0 ) THEN
C PCP, SUBAR, AND SUPAR SCHEDULES COMPLETED
          IF ( Nsub.LT.0 ) Nsub = -Nsub
          IF ( .NOT.Froz.OR..NOT.Eql ) GOTO 1300
C SET UP FOR FROZEN.
          IF ( Eql ) Iplt = iplt1
          Eql = .FALSE.
          Page1 = .TRUE.
          CALL SETEN
          Tt = Ttt(Nfz)
          ipp = Nfz
          IF ( Nfz.EQ.Npt ) GOTO 1150
          Npt = Nfz
          Enn = 1./Wm(Nfz)
          IF ( Nfz.EQ.1 ) GOTO 450
          IF ( Nsub.GT.0 ) THEN
            Nsub = -Nsub
            WRITE (IOOUT,99023)
          ENDIF
          IF ( App(Nfz).LT.App(nptth) ) THEN
            WRITE (IOOUT,99024)
          ELSE
            IF ( Nfz.LT.Npp ) GOTO 1200
            GOTO 700
          ENDIF
          GOTO 1300
        ELSE
          IF ( Eql ) WRITE (IOOUT,99019)
          Npt = nptth
        ENDIF
C SET INDICES AND ESTIMATES FOR NEXT POINT.
 1200   Npt = Npt + 1
        IF ( Eql.OR.(Isv.EQ.-i12.AND..NOT.seql) ) THEN
C THE FOLLOWING STATEMENT WAS ADDED TO TAKE CARE OF A SITUATION
C WHERE EQLBRM WENT SINGULAR WHEN STARTING FROM ESTIMATES WHERE
C BOTH SOLID AND LIQUID WERE INCLUDED.  JULY 27, 1990.
          IF ( Jliq.NE.0.AND.Isv.GT.0 ) Isv = 0
          CALL SETEN
        ENDIF
 1250   ipp = ipp + 1
        IF ( Npt.GT.nptth ) THEN
          IF ( Area ) THEN
            App(Npt) = EXP(appl)
          ELSE
            App(Npt) = Pcp(ipp-nptth)
            IF ( Fac ) App(Npt) = App(Npt)*pinf/Ppp(1)
            IF ( .NOT.Eql.AND.App(Npt).LT.App(Nfz) ) THEN
              WRITE (IOOUT,99020) Nfz
              GOTO 1250
            ENDIF
          ENDIF
          Pp = pinf/App(Npt)
          IF ( Fac ) THEN
            IF ( Area ) THEN
              IF ( isub.LE.Nsub.AND.isub.GT.i01.AND.aratio.GE.Aeat(2) )
     &             THEN
                WRITE (IOOUT,99021) aratio,Aeat(2)
                Npt = Npt - 1
                GOTO 1000
              ENDIF
            ELSEIF ( Npt.GT.nptth.AND.Pcp(ipp-3).LT.Ppp(1)/Ppp(2) ) THEN
              WRITE (IOOUT,99022) Pcp(ipp-3),Ppp(1)/Ppp(2)
              Npt = Npt - 1
              GOTO 650
            ENDIF
          ENDIF
        ENDIF
        GOTO 250
 1300   Npt = 1
C CHECK FOR COMPLETED SCHEDULES -
C 1) CHAMBER PRESSURES(IP = NP)
C 2) CHAMBER TEMPERATURES(IT = NT)
C 3) O/F VALUES(IOF = NOF)
        IF ( Ip.EQ.Np.AND.It.EQ.Nt.AND.iof.EQ.Nof ) GOTO 1400
        WRITE (IOOUT,99019)
        CALL SETEN
        Tt = Ttt(i12)
      ENDDO
      IF ( It.LT.Nt ) THEN
        It = It + 1
        Tt = T(It)
        GOTO 200
      ELSEIF ( iof.LT.Nof ) THEN
        GOTO 100
      ENDIF
 1400 Iplt = MAX(Iplt,iplte)
      RETURN
99001 FORMAT (/' FATAL ERROR!! EITHER mdot OR ac/at MISSING ',
     &        'FOR fac PROBLEM (ROCKET)')
99002 FORMAT (/' WARNING!!  nfz NOT ALLOWED TO BE > 2 IF THE TOTAL',/,
     &        ' NUMBER OF POINTS IS >',i3,' (ROCKET)')
99003 FORMAT (/' INPUT VALUE OF mdot/a =',F12.3,' IS TOO LARGE.'/
     &        ' GIVES CONTRACTION RATIO ESTIMATE LESS THAN 1 (ROCKET)')
99004 FORMAT (/'  ITERATION',9X,'PC',7X,'CONTRACTION RATIO')
99005 FORMAT (5X,I2,7X,F12.2,3X,F12.6)
99006 FORMAT (' ITER',3X,'TEST',3X,'ASSIGNED PINJ',1x,'CALC PINJ',5X,
     &        'PC',7X,'P AT ACAT',3X,'PREV ACAT',2X,'ACAT')
99007 FORMAT (I3,F10.6,1x,4F12.2,2F9.5)
99008 FORMAT (' NEW PC = ',F10.2,2X,'NEW ACAT = ',F9.6,2X,'PJRAT =',
     &        F10.7,' PRACAT =',F10.7)
99009 FORMAT (' END OF CHAMBER ITERATIONS')
99010 FORMAT (/' WARNING!!  DIFFICULTY IN LOCATING THROAT (ROCKET)')
99011 FORMAT (/' USQ=',E15.8,5X,'PVG=',E15.8)
99012 FORMAT (/' WARNING!!  DISCONTINUITY AT THE THROAT (ROCKET)')
99013 FORMAT (' Pinf/Pt =',F9.6)
99014 FORMAT (/,' WARNING!! FOR FROZEN PERFORMANCE, POINTS WERE OMITTED'
     &        ,' WHERE THE ASSIGNED',/,' SUPERSONIC AREA RATIOS WERE ',
     &        'LESS THAN THE VALUE AT POINT nfz =',I3,' (ROCKET)')
99015 FORMAT (/' WARNING!!  AREA RATIO CALCULATION CANNOT BE DONE ',
     &        'BECAUSE GAMMAs',/,' CALCULATION IMPOSSIBLE. (ROCKET)')
99016 FORMAT (/' ITER=',I2,2X,'ASSIGNED AE/AT=',F14.7,3X,'AE/AT=',F14.7,
     &        /,2X,'PC/P=',F14.7,2X,'DELTA LN PCP=',F14.7)
99017 FORMAT (/' WARNING!!  DID NOT CONVERGE FOR AREA RATIO =',F10.5,
     &        ' (ROCKET)')
99018 FORMAT (/' WARNING!!  CALCULATIONS WERE STOPPED BECAUSE NEXT ',
     &        'POINT IS MORE',/,' THAN 50 K BELOW THE TEMPERATURE',
     &        ' RANGE OF A CONDENSED SPECIES (ROCKET)')
99019 FORMAT (////)
99020 FORMAT (/,' WARNING!! FOR FROZEN PERFORMANCE, POINTS WERE OMITTED'
     &        ,' WHERE THE ASSIGNED',/,
     &        ' PRESSURE RATIOS WERE LESS THAN ',
     &        'THE VALUE AT POINT nfz =',I3,' (ROCKET)')
99021 FORMAT (/' WARNING!!  ASSIGNED subae/at =',f10.5,' IS NOT ',
     &        'PERMITTED TO BE GREATER'/' THAN ac/at =',f9.5,
     &        '.  POINT OMITTED (ROCKET)')
99022 FORMAT (/' WARNING!!  ASSIGNED pip =',F10.5,
     &        ' IS NOT PERMITTED'/' TO BE LESS THAN  Pinj/Pc =',f9.5,
     &        '. POINT OMITTED',' (ROCKET)')
99023 FORMAT (/' WARNING!!  FOR FROZEN PERFORMANCE, SUBSONIC AREA ',/,
     &       ' RATIOS WERE OMITTED SINCE nfz IS GREATER THAN 1 (ROCKET)'
     &       )
99024 FORMAT (/' WARNING!!  FREEZING IS NOT ALLOWED AT A SUBSONIC ',
     &        'PRESSURE RATIO FOR nfz GREATER'/' THAN 1. FROZEN ',
     &        'PERFORMANCE CALCULATIONS WERE OMITTED (ROCKET)')
99025 FORMAT (/' AN ASSIGNED AREA RATIO IS < 1 (ROCKET)' )
      END
