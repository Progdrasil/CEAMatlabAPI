      SUBROUTINE TRANIN
C***********************************************************************
C BRINGS IN AND SORTS OUT INPUT FOR TRANSPORT CALCULATIONS
C***********************************************************************
      IMPLICIT NONE
      INCLUDE 'cea.inc'
C LOCAL VARIABLES
      INTEGER i,ii,inds(MAXTR),ir,j,jtape(2),k,k1,k2,kt,kvc,l,loop,m,nms
      LOGICAL change,elc1,elc2,ion1,ion2,setx
      REAL*8 coeff,debye,ekt,enel,enmin,ionic,lamda,omega,prop,qc,ratio,
     &       stcf(MAXTR,MAXTR),stcoef(MAXTR),te,testen,testot,total,
     &       trc(6,3,2),wmols(MAXTR),wmred,xsel,xss(MAXTR)
      REAL*8 DABS,DEXP,DLOG,DSQRT
      SAVE change,coeff,debye,ekt,elc1,elc2,enel,enmin,i,ii,inds,ion1,
     &  ion2,ionic,ir,j,jtape,k,k1,k2,kt,kvc,l,lamda,loop,m,nms,omega,
     &  prop,qc,ratio,setx,stcf,stcoef,te,testen,testot,total,trc,wmols,
     &  wmred,xsel,xss
C
      IF ( .NOT.Eql ) THEN
        IF ( .NOT.Shock ) THEN
          IF ( .NOT.setx ) THEN
            setx = .TRUE.
            Nm = nms
            DO i = 1,Nm
              Xs(i) = xss(i)
              Wmol(i) = wmols(i)
              Ind(i) = inds(i)
            ENDDO
          ENDIF
          GOTO 300
        ELSEIF ( .NOT.Incdeq ) THEN
          IF ( Npt.LE.1 ) THEN
            Nm = Nreac
            DO i = 1,Nm
              j = Jray(i)
              Ind(i) = j
              Wmol(i) = Mw(j)
              Xs(i) = En(j,1)*Wm(1)
            ENDDO
          ENDIF
          GOTO 300
        ENDIF
      ENDIF
C PICK OUT IMPORTANT SPECIES
      Nm = 0
      total = 0.D0
      enmin = 1.0D-11/Wm(Npt)
      testot = 0.999999999D0/Wm(Npt)
      DO i = 1,Lsave
        j = Jcm(i)
        IF ( En(j,Npt).LE.0.D0.AND.j.LE.Ngc ) THEN
          IF ( (Enln(j)-Ennl+25.328436D0).GT.0.D0 ) En(j,Npt)
     &         = DEXP(Enln(j))
        ENDIF
        Nm = Nm + 1
        Ind(Nm) = j
        total = total + En(j,Npt)
        IF ( Mw(j).LT.1.0D0 ) enel = En(j,Npt)
        En(j,Npt) = -En(j,Npt)
      ENDDO
      testen = 1.D0/(Ng*Wm(Npt))
      loop = 0
 100  IF ( total.LE.testot.AND.loop.LE.Ng ) THEN
        loop = loop + 1
        testen = testen/10.
        DO j = 1,Ng
          IF ( En(j,Npt).GE.testen ) THEN
            IF ( Nm.GE.MAXTR ) THEN
              WRITE (IOOUT,99001) Nm,Npt
              GOTO 200
            ELSE
              total = total + En(j,Npt)
              Nm = Nm + 1
              Ind(Nm) = j
              En(j,Npt) = -En(j,Npt)
            ENDIF
          ENDIF
        ENDDO
        IF ( testen.GT.enmin ) GOTO 100
      ENDIF
C CALCULATE MOLE FRACTIONS FROM THE EN(J,NPT)
 200  DO j = 1,Ng
        En(j,Npt) = DABS(En(j,Npt))
      ENDDO
      DO i = 1,Nm
        j = Ind(i)
        Wmol(i) = Mw(j)
        Xs(i) = En(j,Npt)/total
      ENDDO
      IF ( Npt.EQ.Nfz ) THEN
        nms = Nm
        DO i = 1,Nm
          xss(i) = Xs(i)
          wmols(i) = Wmol(i)
          inds(i) = Ind(i)
        ENDDO
        setx = .FALSE.
      ENDIF
C REWRITE REACTIONS TO ELIMINATE TRACE ELEMENTS
      Nr = Nm - Lsave
      IF ( Nr.NE.0 ) THEN
        DO k = 1,MAXTR
          DO m = 1,MAXTR
            Stc(k,m) = 0.0D0
          ENDDO
        ENDDO
        k = 1
        DO i = Lsave + 1,Nm
          Stc(k,i) = -1.0D0
          j = Ind(i)
          DO m = 1,Lsave
            Stc(k,m) = A(m,j)
          ENDDO
          k = k + 1
        ENDDO
        DO i = 1,Nm
          IF ( Xs(i).LT.1.0D-10 ) THEN
            m = 1
            change = .FALSE.
            DO 210 j = 1,Nr
              coeff = Stc(j,i)
              IF ( ABS(coeff).GT.1.0D-05 ) THEN
                IF ( .NOT.change ) THEN
                  change = .TRUE.
                  DO k = 1,Nm
                    stcoef(k) = Stc(j,k)/coeff
                  ENDDO
                  GOTO 210
                ELSE
                  DO k = 1,Nm
                    Stc(j,k) = (Stc(j,k)/coeff) - stcoef(k)
                  ENDDO
                ENDIF
              ENDIF
              DO k = 1,Nm
                stcf(m,k) = Stc(j,k)
              ENDDO
              m = m + 1
 210        CONTINUE
            DO ii = 1,Nm
              DO j = 1,Nr
                Stc(j,ii) = stcf(j,ii)
              ENDDO
            ENDDO
            Nr = m - 1
          ENDIF
        ENDDO
      ENDIF
C FIND TRANSPORT DATA FOR IMPORTANT INTERACTIONS
 300  DO i = 1,Nm
        Con(i) = 0.0
        DO j = 1,Nm
          Eta(i,j) = 0.0
        ENDDO
      ENDDO
      REWIND IOSCH
      DO 400 ir = 1,Ntape
        READ (IOSCH) jtape,trc
        DO 350 k = 1,2
          DO i = 1,Nm
            j = Ind(i)
            IF ( j.EQ.jtape(k) ) THEN
              l = i
              IF ( k.EQ.2 ) THEN
                kvc = 1
 302            kt = 1
                IF ( trc(2,1,kvc).NE.0.E0 ) THEN
                  IF ( trc(2,2,kvc).NE.0.E0 ) THEN
                    IF ( Tt.GT.trc(2,1,kvc) ) kt = 2
                    IF ( trc(2,3,kvc).NE.0. ) THEN
                      IF ( Tt.GT.trc(2,2,kvc) ) kt = 3
                    ENDIF
                  ENDIF
                  prop = EXP(trc(6,kt,kvc)
     &                   +(trc(5,kt,kvc)/Tt+trc(4,kt,kvc))
     &                   /Tt+trc(3,kt,kvc)*Tln)
                  IF ( kvc.EQ.2 ) THEN
                    Con(l) = prop
                    GOTO 400
                  ELSE
                    Eta(l,m) = prop
                    IF ( l.NE.m ) Eta(m,l) = Eta(l,m)
                  ENDIF
                ELSEIF ( kvc.EQ.2 ) THEN
                  GOTO 400
                ENDIF
                kvc = 2
                GOTO 302
              ELSE
                m = i
                GOTO 350
              ENDIF
            ENDIF
          ENDDO
          GOTO 400
 350    CONTINUE
 400  CONTINUE
C MAKE ESTIMATES FOR MISSING DATA
C
C INCLUDES ION CROSS SECTION ESTIMATES
C ESTIMATES FOR  E-ION, ION-ION, E-NEUTRAL, ION-NEUTRAL
C DEBYE SHIELDING WITH IONIC CUTOFF DISTANCE
      IF ( Ions ) THEN
        te = Tt/1000.D0
        ekt = 4.8032D0**2/(Boltz*te)
        qc = 100.D0*(ekt**2)
        xsel = enel/total
        IF ( xsel.LT.1.0D-12 ) xsel = 1.0D-12
        debye = ((22.5D0/Pi)*(Rr/Avgdr*100.D0)*(te/xsel))/ekt**3
        ionic = ((810.D0/(4.0D0*Pi))*(Rr/Avgdr*100D0)*(te/xsel))
     &          **(2.0/3.0)/ekt**2
        lamda = DSQRT(debye+ionic)
        lamda = MAX(lamda,2.71828183D0)
      ENDIF
      DO i = 1,Nm
        k = Ind(i)
        Cprr(i) = Cp(k)
        IF ( .NOT.(Ions.AND.(DABS(A(Nlm,k)).EQ.1.D0).AND.
     &       (Eta(i,i).EQ.0.D0)) ) THEN
          IF ( Eta(i,i).EQ.0.D0 ) THEN
            omega = DLOG(50.D0*Wmol(i)**4.6/Tt**1.4)
            omega = MAX(omega,1.D0)
            Eta(i,i) = Viscns*DSQRT(Wmol(i)*Tt)/omega
          ENDIF
          IF ( Con(i).EQ.0.D0 ) Con(i) = Eta(i,i)
     &         *Rr*(.00375D0+.00132D0*(Cprr(i)-2.5D0))/Wmol(i)
        ENDIF
      ENDDO
      DO i = 1,Nm
        DO 450 j = i,Nm
          ion1 = .FALSE.
          ion2 = .FALSE.
          elc1 = .FALSE.
          elc2 = .FALSE.
          omega = 0.0
          IF ( Eta(i,j).EQ.0. ) Eta(i,j) = Eta(j,i)
          IF ( Eta(j,i).EQ.0. ) Eta(j,i) = Eta(i,j)
          IF ( Eta(i,j).EQ.0. ) THEN
            IF ( Ions ) THEN
C ESTIMATE FOR IONS
              k1 = Ind(i)
              k2 = Ind(j)
              IF ( ABS(A(Nlm,k1)).EQ.1.0 ) ion1 = .TRUE.
              IF ( ABS(A(Nlm,k2)).EQ.1.0 ) ion2 = .TRUE.
              IF ( Wmol(i).LT.1.0 ) elc1 = .TRUE.
              IF ( Wmol(j).LT.1.0 ) elc2 = .TRUE.
              IF ( ion1.AND.ion2 ) omega = 1.36D0*qc*DLOG(lamda)
              IF ( (ion1.AND.elc2).OR.(ion2.AND.elc1) )
     &             omega = 1.29D0*qc*DLOG(lamda)
              IF ( (ion1.AND..NOT.ion2).OR.(ion2.AND..NOT.ion1) )
     &             omega = EXP(6.776-0.4*Tln)
              IF ( omega.NE.0. ) THEN
                wmred = DSQRT(2.0*Tt*Wmol(i)*Wmol(j)/(Wmol(i)+Wmol(j)))
                Eta(i,j) = Viscns*wmred*Pi/omega
                Eta(j,i) = Eta(i,j)
                IF ( i.EQ.j ) THEN
                  Cprr(i) = Cp(k1)
                  Con(i) = Eta(i,i)
     &                     *Rr*(.00375D0+.00132D0*(Cprr(i)-2.5D0))
     &                     /Wmol(i)
                ENDIF
                GOTO 450
              ENDIF
            ENDIF
C ESTIMATE FOR UNLIKE INTERACTIONS FROM RIGID SPHERE ANALOGY
            ratio = DSQRT(Wmol(j)/Wmol(i))
            Eta(i,j) = 5.656854D0*Eta(i,i)
     &                 *SQRT(Wmol(j)/(Wmol(i)+Wmol(j)))
            Eta(i,j) = Eta(i,j)/(1.D0+DSQRT(ratio*Eta(i,i)/Eta(j,j)))**2
            Eta(j,i) = Eta(i,j)
          ENDIF
 450    CONTINUE
      ENDDO
      RETURN
99001 FORMAT (/' WARNING!!  MAXIMUM ALLOWED NO. OF SPECIES',I3,
     &        ' WAS USED IN ',
     &        /' TRANSPORT PROPERTY CALCULATIONS FOR POINT',I3,
     &        '(TRANIN))')
      END
