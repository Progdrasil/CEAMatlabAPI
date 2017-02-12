      SUBROUTINE MATRIX
C***********************************************************************
C SET UP ITERATION OR DERIVATIVE MATRIX.
C***********************************************************************
      IMPLICIT NONE
      INCLUDE 'cea.inc'
C LOCAL VARIABLES
      INTEGER i,iq,iq2,iq3,isym,j,k,kk,kmat
      REAL*8 energyl,f,h,ss,sss,term,term1
      SAVE energyl,f,h,i,iq,iq2,iq3,isym,j,k,kk,kmat,ss,sss,term,term1
C
      iq = Nlm + Npr
      Iq1 = iq + 1
      iq2 = Iq1 + 1
      iq3 = iq2 + 1
      kmat = iq3
      IF ( .NOT.Convg.AND.Tp ) kmat = iq2
      Imat = kmat - 1
C CLEAR MATRIX STORAGES TO ZERO
      DO i = 1,Imat
        DO k = 1,kmat
          G(i,k) = 0.0D0
        ENDDO
      ENDDO
      G(iq2,Iq1) = 0.D0
      sss = 0.D0
      Hsum(Npt) = 0.D0
C BEGIN SET-UP OF ITERATION OR DERIVATIVE MATRIX
      DO j = 1,Ng
        Mu(j) = H0(j) - S(j) + Enln(j) + Tm
        IF ( En(j,Npt).NE.0.D0 ) THEN
          h = H0(j)*En(j,Npt)
          f = Mu(j)*En(j,Npt)
          ss = h - f
          term1 = h
          IF ( kmat.EQ.iq2 ) term1 = f
          DO i = 1,Nlm
            IF ( A(i,j).NE.0. ) THEN
              term = A(i,j)*En(j,Npt)
              DO k = i,Nlm
                G(i,k) = G(i,k) + A(k,j)*term
              ENDDO
              G(i,Iq1) = G(i,Iq1) + term
              G(i,iq2) = G(i,iq2) + A(i,j)*term1
              IF ( .NOT.(Convg.OR.Tp) ) THEN
                G(i,iq3) = G(i,iq3) + A(i,j)*f
                IF ( Sp ) G(iq2,i) = G(iq2,i) + A(i,j)*ss
              ENDIF
            ENDIF
          ENDDO
          IF ( kmat.NE.iq2 ) THEN
            IF ( Convg.OR.Hp ) THEN
              G(iq2,iq2) = G(iq2,iq2) + H0(j)*h
              IF ( .NOT.Convg ) THEN
                G(iq2,iq3) = G(iq2,iq3) + H0(j)*f
                G(Iq1,iq3) = G(Iq1,iq3) + f
              ENDIF
            ELSE
              G(iq2,Iq1) = G(iq2,Iq1) + ss
              G(iq2,iq2) = G(iq2,iq2) + H0(j)*ss
              G(iq2,iq3) = G(iq2,iq3) + Mu(j)*ss
              G(Iq1,iq3) = G(Iq1,iq3) + f
            ENDIF
          ENDIF
          G(Iq1,iq2) = G(Iq1,iq2) + term1
        ENDIF
      ENDDO
C CONDENSED SPECIES
      IF ( Npr.NE.0 ) THEN
        DO k = 1,Npr
          j = Jcond(k)
          kk = Nlm + k
          Mu(j) = H0(j) - S(j)
          DO i = 1,Nlm
            G(i,kk) = A(i,j)
            G(i,kmat) = G(i,kmat) - A(i,j)*En(j,Npt)
          ENDDO
          G(kk,iq2) = H0(j)
          G(kk,kmat) = Mu(j)
          Hsum(Npt) = Hsum(Npt) + H0(j)*En(j,Npt)
          IF ( Sp ) THEN
            sss = sss + S(j)*En(j,Npt)
            G(iq2,kk) = S(j)
          ENDIF
        ENDDO
      ENDIF
      sss = sss + G(iq2,Iq1)
      Hsum(Npt) = Hsum(Npt) + G(Iq1,iq2)
      G(Iq1,Iq1) = Sumn - Enn
C REFLECT SYMMETRIC PORTIONS OF THE MATRIX
      isym = Iq1
      IF ( Hp.OR.Convg ) isym = iq2
      DO i = 1,isym
CDIR$ IVDEP
        DO j = i,isym
          G(j,i) = G(i,j)
        ENDDO
      ENDDO
C COMPLETE THE RIGHT HAND SIDE
      IF ( .NOT.Convg ) THEN
        DO i = 1,Nlm
          G(i,kmat) = G(i,kmat) + B0(i) - G(i,Iq1)
        ENDDO
        G(Iq1,kmat) = G(Iq1,kmat) + Enn - Sumn
C COMPLETE ENERGY ROW AND TEMPERATURE COLUMN
        IF ( kmat.NE.iq2 ) THEN
          IF ( Sp ) energyl = S0 + Enn - Sumn - sss
          IF ( Hp ) energyl = Hsub0/Tt - Hsum(Npt)
          G(iq2,iq3) = G(iq2,iq3) + energyl
          G(iq2,iq2) = G(iq2,iq2) + Cpsum
        ENDIF
      ELSE
        IF ( Pderiv ) THEN
C PDERIV = .TRUE.-- SET UP MATRIX TO SOLVE FOR DLVPT
          G(Iq1,iq2) = Enn
          DO i = 1,iq
            G(i,iq2) = G(i,Iq1)
          ENDDO
        ENDIF
        G(iq2,iq2) = G(iq2,iq2) + Cpsum
      ENDIF
      IF ( Vol.AND..NOT.Convg ) THEN
C CONSTANT VOLUME MATRIX
        IF ( kmat.EQ.iq2 ) THEN
          DO i = 1,iq
            G(i,Iq1) = G(i,iq2)
          ENDDO
        ELSE
CDIR$ IVDEP
          DO i = 1,iq
            G(Iq1,i) = G(iq2,i) - G(Iq1,i)
            G(i,Iq1) = G(i,iq2) - G(i,Iq1)
            G(i,iq2) = G(i,iq3)
          ENDDO
          G(Iq1,Iq1) = G(iq2,iq2) - G(Iq1,iq2) - G(iq2,Iq1)
          G(Iq1,iq2) = G(iq2,iq3) - G(Iq1,iq3)
          IF ( Hp ) G(Iq1,iq2) = G(Iq1,iq2) + Enn
        ENDIF
        kmat = Imat
        Imat = Imat - 1
      ENDIF
      END
