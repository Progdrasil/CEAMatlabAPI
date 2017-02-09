      SUBROUTINE THERMP
C***********************************************************************
C ASSIGNED THERMODYNAMIC STATES.  HP,SP,TP,UV,SV, AND TV PROBLEMS.
C***********************************************************************
      IMPLICIT NONE
      INCLUDE 'cea.inc'
C LOCAL VARIABLES
      INTEGER iof
      LOGICAL uv,tv,sv
      SAVE iof
C
      EQUIVALENCE (Hp,Uv)
      EQUIVALENCE (Tp,Tv)
      EQUIVALENCE (Sp,Sv)
      Eql = .TRUE.
      DO 100 iof = 1,Nof
        Oxfl = Oxf(iof)
        CALL NEWOF
C SET ASSIGNED P OR VOLUME
        DO Ip = 1,Np
          Pp = P(Ip)
C SET ASSIGNED T
          DO It = 1,Nt
            Vv = V(Ip)
            Tt = T(It)
            CALL EQLBRM
            IF ( Npt.EQ.0 ) GOTO 200
            IF ( Trnspt.AND.Tt.NE.0. ) CALL TRANP
            Isv = 0
            IF ( Ip.NE.Np.OR.It.NE.Nt.AND.Tt.NE.0. ) THEN
              Isv = Npt
              IF ( Npt.NE.NCOL ) GOTO 10
            ENDIF
            IF ( .NOT.Hp ) WRITE (IOOUT,99001)
            IF ( Hp ) WRITE (IOOUT,99002)
            IF ( .NOT.Vol ) THEN
              IF ( Hp ) WRITE (IOOUT,99006)
              IF ( Tp ) WRITE (IOOUT,99007)
              IF ( Sp ) WRITE (IOOUT,99008)
            ELSE
              IF ( Uv ) WRITE (IOOUT,99003)
              IF ( Tv ) WRITE (IOOUT,99004)
              IF ( Sv ) WRITE (IOOUT,99005)
            ENDIF
            CALL OUT1
            WRITE (IOOUT,99009)
            CALL OUT2
            IF ( Trnspt ) CALL OUT4
            CALL OUT3
            Iplt = MIN(Iplt+Npt,500)
            IF ( Isv.EQ.0.AND.iof.EQ.Nof ) GOTO 200
            WRITE (IOOUT,99010)
            Npt = 0
 10         Npt = Npt + 1
            IF ( .NOT.Tp.AND.Tt.NE.0. ) T(1) = Tt
            IF ( Nt.EQ.1.AND.Np.EQ.1 ) GOTO 100
            IF ( Ip.EQ.1.AND.It.EQ.1 ) Isv = -Isv
            IF ( Nt.NE.1 ) THEN
              IF ( It.EQ.Nt.OR.Tt.EQ.0. ) Isv = 0
            ENDIF
            CALL SETEN
          ENDDO
        ENDDO
 100  CONTINUE
 200  RETURN
99001 FORMAT (////15X,'THERMODYNAMIC EQUILIBRIUM PROPERTIES AT ASSIGNED'
     &        )
99002 FORMAT (////9X,
     &     'THERMODYNAMIC EQUILIBRIUM COMBUSTION PROPERTIES AT ASSIGNED'
     &     )
99003 FORMAT (/36X,' VOLUME'/)
99004 FORMAT (/28X,'TEMPERATURE AND VOLUME'/)
99005 FORMAT (/30X,'ENTROPY AND VOLUME'/)
99006 FORMAT (/34X,' PRESSURES'/)
99007 FORMAT (/27X,'TEMPERATURE AND PRESSURE'/)
99008 FORMAT (/29X,'ENTROPY AND PRESSURE'/)
99009 FORMAT (/' THERMODYNAMIC PROPERTIES'/)
99010 FORMAT (////)
      END
