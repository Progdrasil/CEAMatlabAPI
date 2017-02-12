      SUBROUTINE FROZEN
C***********************************************************************
C CALCULATE PROPERTIES WITH FROZEN COMPOSITION AT ASSIGNED ENTROPY
C AND PRESSURE.  CALLED FROM ROCKET.
C***********************************************************************
      IMPLICIT NONE
      INCLUDE 'cea.inc'
C LOCAL VARIABLES
      INTEGER i,inc,iter,j,k,nnn
      REAL*8 DABS,DEXP,DLOG
      REAL*8 dlnt,dlpm
      SAVE dlnt,dlpm,i,inc,iter,j,k,nnn
C
      Convg = .FALSE.
      Tln = DLOG(Tt)
      dlpm = DLOG(Pp*Wm(Nfz))
      nnn = Npt
      Npt = Nfz
      DO j = 1,Ng
        IF ( En(j,Nfz).NE.0.D0 ) Deln(j) = -(DLOG(En(j,Nfz))+dlpm)
      ENDDO
      DO iter = 1,8
        Ssum(nnn) = 0.D0
        Cpsum = 0.D0
        CALL CPHS
        DO j = 1,Ng
          Cpsum = Cpsum + En(j,Nfz)*Cp(j)
          Ssum(nnn) = Ssum(nnn) + En(j,Nfz)*(S(j)+Deln(j))
        ENDDO
        IF ( Npr.NE.0 ) THEN
          DO k = 1,Npr
            j = Jcond(k)
            Cpsum = Cpsum + En(j,Nfz)*Cp(j)
            Ssum(nnn) = Ssum(nnn) + En(j,Nfz)*S(j)
          ENDDO
        ENDIF
        IF ( Convg ) THEN
          Npt = nnn
          Hsum(Npt) = 0.D0
          DO j = 1,Ngc
            Hsum(Npt) = Hsum(Npt) + En(j,Nfz)*H0(j)
          ENDDO
          Hsum(Npt) = Hsum(Npt)*Tt
          Ttt(Npt) = Tt
          Gammas(Npt) = Cpsum/(Cpsum-1./Wm(Nfz))
          Vlm(Npt) = Rr*Tt/(Wm(Nfz)*Pp)
          Wm(Npt) = Wm(Nfz)
          Dlvpt(Npt) = -1.
          Dlvtp(Npt) = 1.
          Totn(Npt) = Totn(Nfz)
          Ppp(Npt) = Pp
          Cpr(Npt) = Cpsum
          IF ( Tt.GE.(Tg(1)*.8D0) ) THEN
            DO i = Ngp1,Ngc
              IF ( En(i,Nfz).NE.0. ) THEN
                inc = i - Ng
                IF ( Tt.LT.(Temp(1,inc)-50.).OR.Tt.GT.(Temp(2,inc)+50.)
     &               ) GOTO 100
              ENDIF
            ENDDO
            GOTO 200
          ENDIF
          GOTO 100
        ELSE
          dlnt = (Ssum(Nfz)-Ssum(nnn))/Cpsum
          Tln = Tln + dlnt
          IF ( DABS(dlnt).LT.0.5D-4 ) Convg = .TRUE.
          Tt = DEXP(Tln)
        ENDIF
      ENDDO
      WRITE (IOOUT,99001)
 100  Tt = 0.
      Npt = Npt - 1
 200  RETURN
99001 FORMAT (/' FROZEN DID NOT CONVERGE IN 8 ITERATIONS (FROZEN)')
      END
