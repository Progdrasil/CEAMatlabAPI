      SUBROUTINE SETEN
C***********************************************************************
C USE COMPOSITIONS FROM PREVIOUS POINT AS INITIAL ESTIMATES FOR
C CURRENT POINT NPT.  IF -
C  ISV>0  USE COMPOSITIONS FROM POINT ISV.
C  ISV<0  SAVE COMPOSITIONS FROM POINT -ISV FOR POSSIBLE LATER USE.
C         ALSO USE COMPOSITIONS FROM POINT -ISV FOR NPT.
C  ISV=0  USE COMPOSITIONS SAVED WHEN ISV<0.
C***********************************************************************
      IMPLICIT NONE
      INCLUDE 'cea.inc'
C LOCAL VARIABLES
      INTEGER j,lsav
      REAL*8 tsave
      REAL*8 DEXP
      SAVE j,lsav,tsave
C
      IF ( Isv.LT.0 ) THEN
C FIRST T--SAVE COMPOSITIONS FOR FUTURE POINTS WITH THIS T
        Isv = -Isv
        tsave = Ttt(Isv)
        Ensave = Enn
        Enlsav = Ennl
        lsav = Lsave
        DO j = 1,Ng
          Sln(j) = Enln(j)
        ENDDO
        DO j = 1,Ng
          En(j,Npt) = En(j,Isv)
        ENDDO
        Npr = 0
        DO j = Ngp1,Ngc
          Sln(j) = En(j,Isv)
          En(j,Npt) = Sln(j)
          IF ( Jliq.EQ.j ) THEN
            En(Jsol,Npt) = En(Jsol,Isv) + En(Jliq,Isv)
            En(Jliq,Npt) = 0.
            Jsol = 0
            Jliq = 0
            tsave = tsave - 5.
            Tt = tsave
            Sln(j) = 0.
          ELSEIF ( En(j,Npt).GT.0. ) THEN
            Npr = Npr + 1
            Jcond(Npr) = j
          ENDIF
        ENDDO
      ELSEIF ( Isv.EQ.0 ) THEN
C NEXT POINT FIRST T IN SCHEDULE, USE PREVIOUS COMPOSITIONS FOR THIS T
        Jsol = 0
        Jliq = 0
        Enn = Ensave
        Ennl = Enlsav
        Lsave = lsav
        Npr = 0
        DO j = Ngp1,Ngc
          En(j,Npt) = Sln(j)
          IF ( En(j,Npt).GT.0.D0 ) THEN
            Npr = Npr + 1
            Jcond(Npr) = j
          ENDIF
        ENDDO
        DO j = 1,Ng
          En(j,Npt) = 0.
          Enln(j) = Sln(j)
          IF ( Sln(j).NE.0. ) THEN
            IF ( (Enln(j)-Ennl+18.5).GT.0. ) En(j,Npt) = DEXP(Enln(j))
          ENDIF
        ENDDO
        IF ( .NOT.Tp ) Tt = tsave
        Sumn = Enn
      ELSEIF ( Isv.GT.0 ) THEN
C USE COMPOSITIONS FROM PREVIOUS POINT
        DO j = 1,Ngc
          En(j,Npt) = En(j,Isv)
        ENDDO
      ENDIF
      END
