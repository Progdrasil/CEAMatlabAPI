      SUBROUTINE NEWOF
C***********************************************************************
C CALCULATE NEW VALUES OF B0 AND HSUB0 FOR NEW OF RATIO
C***********************************************************************
      IMPLICIT NONE
      INCLUDE 'cea.inc'
C LOCAL VARIABLES
      INTEGER i,j
      REAL*8 assval,bigb,bratio,dbi,smalb,tem,v1,v2
      REAL*8 DABS,DLOG
      SAVE assval,bigb,bratio,dbi,i,j,smalb,tem,v1,v2
C
      IF ( .NOT.Short ) WRITE (IOOUT,99001) Oxfl
      Eqrat = 0.
      tem = Oxfl + 1.
      v2 = (Oxfl*Vmin(1)+Vmin(2))/tem
      v1 = (Oxfl*Vpls(1)+Vpls(2))/tem
      IF ( v2.NE.0. ) Eqrat = DABS(v1/v2)
      DO i = 1,Nlm
        B0(i) = (Oxfl*B0p(i,1)+B0p(i,2))/tem
        dbi = DABS(B0(i))
        IF ( i.EQ.1 ) THEN
          bigb = dbi
          smalb = dbi
        ELSEIF ( dbi.NE.0. ) THEN
          IF ( dbi.LT.smalb ) smalb = dbi
          IF ( dbi.GT.bigb ) bigb = dbi
        ENDIF
      ENDDO
      Bcheck = bigb*.000001D0
C CALCUALTE MOLECULAR WEIGHT OF TOTAL REACTANT, WMIX.
      IF ( Am(1).NE.0.0.AND.Am(2).NE.0.0 ) THEN
        Wmix = (Oxfl+1.)*Am(1)*Am(2)/(Am(1)+Oxfl*Am(2))
      ELSE
        Wmix = Am(2)
        IF ( Am(2).EQ.0.0 ) Wmix = Am(1)
      ENDIF
      Npt = 1
C IF ASSIGNED U OR H NOT GIVEN IN PROB DATA, INITIAL HSUB0 = 1.D30
      IF ( Size.EQ.0. ) assval = Hsub0
      IF ( assval.GE.1.D30 ) Hsub0 = (Oxfl*Hpp(1)+Hpp(2))/tem
C NOTE THAT "BRATIO" IS "BRATIO" IN SEC 3.2 IN RP-1311.
      bratio = smalb/bigb
      Size = 18.420681D0
      IF ( bratio.LT.1.D-5 ) Size = DLOG(1000.D0/bratio)
      Jsol = 0
      Jliq = 0
      IF ( .NOT.Short ) THEN
        WRITE (IOOUT,99002)
        IF ( Vol ) WRITE (IOOUT,99003)
        IF ( .NOT.Vol ) WRITE (IOOUT,99004)
        WRITE (IOOUT,99005) Hpp(2),Hpp(1),Hsub0
        WRITE (IOOUT,99006)
      ENDIF
      DO i = 1,Nlm
        j = Jcm(i)
        IF ( .NOT.Short ) WRITE (IOOUT,99007) Prod(j),B0p(i,2),B0p(i,1),
     &                           B0(i)
      ENDDO
      RETURN
99001 FORMAT (/' O/F = ',F10.6)
99002 FORMAT (/,23X,'EFFECTIVE FUEL',5X,'EFFECTIVE OXIDANT',8X,
     &        'MIXTURE')
99003 FORMAT (' INTERNAL ENERGY',11X,'u(2)/R',14X,'u(1)/R',14X,'u0/R')
99004 FORMAT (' ENTHALPY',18X,'h(2)/R',14X,'h(1)/R',15X,'h0/R')
99005 FORMAT (' (KG-MOL)(K)/KG',4X,E18.8,2E20.8)
99006 FORMAT (/' KG-FORM.WT./KG',13X,'bi(2)',15X,'bi(1)',15X,'b0i')
99007 FORMAT (1X,A16,3E20.8)
      END
