      SUBROUTINE TRANP
C***********************************************************************
C CALCULATES GAS TRANSPORT PROPERTIES
C
C   NUMBER OF GASEOUS SPECIES = NM   (MAXIMUM MAXTR)
C   NUMBER OF CHEMICAL REACTIONS = NR (NM - NLM)
C   ARRAY OF STOICHIOMETRIC COEFFICIENTS = STC
C***********************************************************************
      IMPLICIT NONE
      INCLUDE 'cea.inc'
C LOCAL VARIABLES
      INTEGER i,i1,j,jj,k,m,mm,nlmm,nmm
      REAL*8 cpreac,delh(MAXTR),gmat(MAXMAT,MAXMAT+1),phi(MAXTR,MAXTR),
     &       psi(MAXTR,MAXTR),reacon,rtpd(MAXTR,MAXTR),stx(MAXTR),
     &       stxij(MAXTR,MAXTR),sumc,sumv,wtmol,xskm(MAXTR,MAXTR)
      REAL*8 DABS
      SAVE cpreac,delh,gmat,i,i1,j,jj,k,m,mm,nlmm,nmm,phi,psi,reacon,
     &  rtpd,stx,stxij,sumc,sumv,wtmol,xskm
C
      CALL TRANIN
C CALCULATE VISCOSITY AND FROZEN THERMAL CONDUCTIVITY
      nmm = Nm - 1
      DO i = 1,Nm
        rtpd(i,i) = 0.D0
        phi(i,i) = 1.D0
        psi(i,i) = 1.D0
      ENDDO
      Confro(Npt) = 0.D0
      Vis(Npt) = 0.D0
      DO i = 1,nmm
        i1 = i + 1
CDIR$ IVDEP
        DO j = i1,Nm
          sumc = 2.D0/(Eta(i,j)*(Wmol(i)+Wmol(j)))
          phi(i,j) = sumc*Wmol(j)*Eta(i,i)
          phi(j,i) = sumc*Wmol(i)*Eta(j,j)
          sumc = (Wmol(i)+Wmol(j))**2
          psi(i,j) = phi(i,j)
     &               *(1.D0+2.41D0*(Wmol(i)-Wmol(j))*(Wmol(i)-.142D0*
     &               Wmol(j))/sumc)
          psi(j,i) = phi(j,i)
     &               *(1.D0+2.41D0*(Wmol(j)-Wmol(i))*(Wmol(j)-.142D0*
     &               Wmol(i))/sumc)
        ENDDO
      ENDDO
      DO i = 1,Nm
        sumc = 0.D0
        sumv = 0.D0
        DO j = 1,Nm
          sumc = sumc + psi(i,j)*Xs(j)
          sumv = sumv + phi(i,j)*Xs(j)
        ENDDO
        Vis(Npt) = Vis(Npt) + Eta(i,i)*Xs(i)/sumv
        Confro(Npt) = Confro(Npt) + Con(i)*Xs(i)/sumc
      ENDDO
      IF ( Eql.AND.Nr.GT.0 ) THEN
C CALCULATE REACTION HEAT CAPACITY AND THERMAL CONDUCTIVITY
        m = Nr + 1
        DO i = 1,Nr
          delh(i) = 0.0D0
          DO k = 1,Lsave
            j = Jcm(k)
            delh(i) = Stc(i,k)*H0(j) + delh(i)
          ENDDO
          nlmm = Lsave + 1
          DO k = nlmm,Nm
            j = Ind(k)
            delh(i) = Stc(i,k)*H0(j) + delh(i)
          ENDDO
          G(i,m) = delh(i)
        ENDDO
        DO i = 1,MAXTR
          DO j = 1,MAXTR
            IF ( DABS(Stc(i,j)).LT.1.0D-6 ) Stc(i,j) = 0.0D0
          ENDDO
        ENDDO
        jj = Nm - 1
        DO k = 1,jj
          mm = k + 1
          DO m = mm,Nm
            rtpd(k,m) = Wmol(k)*Wmol(m)/(1.1*Eta(k,m)*(Wmol(k)+Wmol(m)))
            xskm(k,m) = Xs(k)*Xs(m)
            xskm(m,k) = xskm(k,m)
            rtpd(m,k) = rtpd(k,m)
          ENDDO
        ENDDO
        DO i = 1,Nr
          DO j = i,Nr
            G(i,j) = 0.0D0
            gmat(i,j) = 0.0D0
          ENDDO
        ENDDO
        DO k = 1,jj
          mm = k + 1
          DO m = mm,Nm
            IF ( Xs(k).GE.1.0D-10.AND.Xs(m).GE.1.0D-10 ) THEN
              DO j = 1,Nr
                IF ( (Stc(j,k).EQ.0.D0).AND.(Stc(j,m).EQ.0.D0) ) stx(j)
     &               = 0.D0
                IF ( (Stc(j,k).NE.0.D0).OR.(Stc(j,m).NE.0.D0) ) stx(j)
     &               = Xs(m)*Stc(j,k) - Xs(k)*Stc(j,m)
              ENDDO
              DO i = 1,Nr
                DO j = i,Nr
                  stxij(i,j) = stx(i)*stx(j)/xskm(k,m)
                  G(i,j) = G(i,j) + stxij(i,j)
                  gmat(i,j) = gmat(i,j) + stxij(i,j)*rtpd(k,m)
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDDO
        m = 1 + Nr
        DO i = 1,Nr
CDIR$ IVDEP
          DO j = i,Nr
            G(j,i) = G(i,j)
          ENDDO
          G(i,m) = delh(i)
        ENDDO
        Imat = Nr
        CALL GAUSS
        cpreac = 0.D0
        DO i = 1,Nr
          G(i,m) = delh(i)
          cpreac = cpreac + R*delh(i)*X(i)
CDIR$ IVDEP
          DO j = i,Nr
            G(i,j) = gmat(i,j)
            G(j,i) = G(i,j)
          ENDDO
        ENDDO
        CALL GAUSS
        reacon = 0.D0
        DO i = 1,Nr
          reacon = reacon + R*delh(i)*X(i)
        ENDDO
        reacon = .6D0*reacon
      ELSE
        cpreac = 0.D0
        reacon = 0.D0
      ENDIF
C CALCULATE OTHER ANSWERS
      Cpfro(Npt) = 0.D0
      wtmol = 0.D0
      DO i = 1,Nm
        Cpfro(Npt) = Cpfro(Npt) + Xs(i)*Cprr(i)
        wtmol = wtmol + Xs(i)*Wmol(i)
      ENDDO
      Cpfro(Npt) = Cpfro(Npt)*R/wtmol
      Confro(Npt) = Confro(Npt)/1000.D0
      IF ( .NOT.Siunit ) Confro(Npt) = Confro(Npt)/4.184D0
      Vis(Npt) = Vis(Npt)/1000.D0
      Prfro(Npt) = Vis(Npt)*Cpfro(Npt)/Confro(Npt)
      IF ( Eql ) THEN
        cpreac = cpreac/wtmol
        reacon = reacon/1000.D0
        Cpeql(Npt) = cpreac + Cpfro(Npt)
        Coneql(Npt) = Confro(Npt) + reacon
        Preql(Npt) = Vis(Npt)*Cpeql(Npt)/Coneql(Npt)
      ENDIF
      END
