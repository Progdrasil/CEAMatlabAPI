      SUBROUTINE REACT
C***********************************************************************
C READ AND PROCESS REACTANT RECORDS.  CALLED FROM SUBROUTINE INPUT.
C***********************************************************************
      IMPLICIT NONE
      INCLUDE 'cea.inc'
C LOCAL VARIABLES
      CHARACTER*6 date
      CHARACTER*2 el(5)
      CHARACTER*15 sub
      INTEGER i,icf,ifaz,ifrmla,itot,j,jj,k,kk,kr,l,n,nall,nint,nj,
     &        ntgas,ntot
      LOGICAL fuel,rcoefs,wdone(2),hok
      REAL*8 bb(5),dat(35),dift,eform,pcwt,rcf(9,3),rm,t1,t2
      REAL*8 DABS,DLOG,t1save,t2save
      SAVE bb,dat,date,dift,eform,el,fuel,i,icf,ifaz,ifrmla,itot,j,jj,k,
     &  kk,kr,l,n,nall,nint,nj,ntgas,ntot,pcwt,rcf,rcoefs,rm,sub,t1,t2,
     &  wdone
C
      DO k = 1,2
        wdone(k) = .FALSE.
        Wp(k) = 0.
        Hpp(k) = 0.
        Vpls(k) = 0.
        Vmin(k) = 0.
        Am(k) = 0.
        Rh(k) = 0.
        DO j = 1,MAXEL
          Elmt(j) = ' '
          B0p(j,k) = 0.
        ENDDO
      ENDDO
      DO i = 1,MAXEL
        dat(i) = 0.
      ENDDO
C IF OXIDANT, KR = 1
C IF FUEL, KR = 2
      DO n = 1,Nreac
        hok = .false.
        t1save = 20000.d0
        t2save = 0.d0
        rcoefs = .TRUE.
        IF ( Energy(n).EQ.'lib'.OR.Rnum(n,1).EQ.0. ) THEN
          Tt = Rtemp(n)
          REWIND IOTHM
          READ (IOTHM) Tg,ntgas,ntot,nall
          DO 20 itot = 1,nall
            IF ( itot.LE.ntot ) THEN
              icf = 3
              IF ( itot.GT.ntgas ) icf = 1
              READ (IOTHM) sub,nint,date,(el(j),bb(j),j=1,5),ifaz,t1,t2,
     &                     rm,((rcf(i,j),i=1,9),j=1,icf)
            ELSE
              READ (IOTHM) sub,nint,date,(el(j),bb(j),j=1,5),ifaz,t1,t2,
     &                     rm,eform
              IF ( nint.GT.0 ) READ (IOTHM) ((rcf(i,j),i=1,9),j=1,nint)
            ENDIF
            IF ( sub.EQ.Rname(n).OR.sub.EQ.'*'//Rname(n) ) THEN
              IF ( nint.EQ.0 ) THEN
                rcoefs = .FALSE.
                hok = .true.
                Enth(n) = eform*1000.D0/Rr
                IF ( Tt.EQ.0 ) THEN
                  Tt = t1
                  Rtemp(n) = t1
                ELSE
                  dift = DABS(Tt-t1)
                  IF ( dift.GT.01d0 ) THEN
                    IF ( dift.GT.10.d0 ) THEN
                      WRITE (IOOUT,99001) Rname(n),t1,Tt
                      Nlm = 0
                      hok = .false.
                      GOTO 200
                    ELSE
                      WRITE (IOOUT,99002) Rname(n),t1,Tt
                      Tt = t1
                      Rtemp(n) = t1
                    ENDIF
                  ENDIF
                ENDIF
              ELSE
                if (ifaz.LE.0 ) then
                  t1save = min(t1save,.8d0*tg(1))
                  t2save = max(t2save,1.2d0*t2)
                else  
                  t1save = min(t1save,t1-.001d0)
                  t2save = max(t2save,t2+.001d0)
                endif
                if ( t1save .lt.Tt .and. t2save.gt.Tt ) hok = .true.
              ENDIF
              DO j = 1,5
                IF ( bb(j).EQ.0. ) GOTO 5
                Nfla(n) = j
                Ratom(n,j) = el(j)
                Rnum(n,j) = bb(j)
              ENDDO
 5            IF ( Tt.EQ.0. ) THEN
                IF ( .NOT.Hp ) GOTO 50
                WRITE (IOOUT,99004) n
                Nlm = 0
                GOTO 200
              ENDIF
              IF ( rcoefs.and.hok ) THEN
                Tln = DLOG(Tt)
                l = 1
                IF ( ifaz.LE.0 ) THEN
                  IF ( Tt.GT.Tg(2) ) l = 2
                  IF ( Tt.GT.Tg(3) ) l = 3
                ENDIF
                Enth(n) = (((((rcf(7,l)/5.D0)*Tt+rcf(6,l)/4.D0)*Tt+rcf(5
     &                    ,l)/3.D0)*Tt+rcf(4,l)/2.D0)*Tt+rcf(3,l))
     &                    *Tt - rcf(1,l)/Tt + rcf(2,l)*Tln + rcf(8,l)
                IF ( Vol.AND.ifaz.LE.0 ) Enth(n) = Enth(n) - Tt
              ENDIF
              if (hok) GOTO 50
            ENDIF
 20       CONTINUE
          if (.not.hok) then
            WRITE (IOOUT,99010) Tt,Rname(n),t1save,t2save
            Energy(n) = ' '
            Nlm = 0
            goto 200
          endif
        ENDIF
 50     ifrmla = Nfla(n)
        IF ( Fox(n)(:1).EQ.'f' ) THEN
          fuel = .TRUE.
          kr = 2
          Fox(n) = 'FUEL'
        ELSEIF ( Fox(n)(:4).EQ.'name' ) THEN
          fuel = .TRUE.
          kr = 2
          Fox(n) = 'NAME'
        ELSE
          kr = 1
          Fox(n) = 'OXIDANT'
        ENDIF
        DO j = 1,MAXEL
          dat(j) = 0.
        ENDDO
C STORE ATOMIC SYMBOLS IN ELMT ARRAY.
C CALCULATE MOLECULAR WEIGHT.
C TEMPORARILY STORE ATOMIC VALENCE IN X.
        rm = 0.D0
        DO 100 jj = 1,ifrmla
          DO j = 1,MAXEL
            nj = j
            IF ( Elmt(j).EQ.' ' ) GOTO 60
            IF ( Ratom(n,jj).EQ.Elmt(j) ) GOTO 80
          ENDDO
 60       Nlm = nj
          Elmt(j) = Ratom(n,jj)
 80       DO kk = 1,100
            IF ( Symbol(kk).EQ.Ratom(n,jj) ) THEN
              rm = rm + Rnum(n,jj)*Atmwt(kk)
              Atwt(j) = Atmwt(kk)
              X(j) = Valnce(kk)
              dat(j) = dat(j) + Rnum(n,jj)
              GOTO 100
            ENDIF
          ENDDO
          WRITE (IOOUT,99005) Ratom(n,jj)
          Nlm = 0
          GOTO 200
 100    CONTINUE
        IF ( Pecwt(n).LT.0. ) THEN
          Pecwt(n) = 0.
          IF ( .NOT.Moles.AND..NOT.wdone(kr) ) THEN
            wdone(kr) = .TRUE.
            Pecwt(n) = 100.
            WRITE (IOOUT,99006) n
          ELSE
            WRITE (IOOUT,99007) n
            Nlm = 0
            GOTO 200
          ENDIF
        ENDIF
C ADD CONTRIBUTIONS TO WP(K), HPP(K), AM(K), AND B0P(I,K)
        IF ( Pecwt(n).GT.0. ) wdone(kr) = .TRUE.
        pcwt = Pecwt(n)
        IF ( Moles ) pcwt = pcwt*rm
        Wp(kr) = Wp(kr) + pcwt
        IF ( rm.LE.0.D0 ) THEN
          Nlm = 0
          GOTO 200
        ELSE
          Hpp(kr) = Hpp(kr) + Enth(n)*pcwt/rm
          Am(kr) = Am(kr) + pcwt/rm
          IF ( Dens(n).NE.0. ) THEN
            Rh(kr) = Rh(kr) + pcwt/Dens(n)
          ELSE
            Rh(1) = 0.
            Rh(2) = 0.
          ENDIF
          DO j = 1,Nlm
            B0p(j,kr) = dat(j)*pcwt/rm + B0p(j,kr)
          ENDDO
          Rmw(n) = rm
        ENDIF
      ENDDO
      IF ( .NOT.fuel ) THEN
C 100 PERCENT OXIDANT, SWITCH INDICES
        DO n = 1,Nreac
          Fox(n) = ' '
        ENDDO
        Wp(2) = Wp(1)
        Wp(1) = 0.
        Hpp(2) = Hpp(1)
        Am(2) = Am(1)
        Am(1) = 0.
        DO j = 1,Nlm
          B0p(j,2) = B0p(j,1)
        ENDDO
      ENDIF
      IF ( Nlm.NE.0 ) THEN
C NORMALIZE HPP(KKR),AM(KR),B0P(I,KR), AND PECWT(N).
C CALCULATE V+(KR), AND V-(KR)
        DO kr = 1,2
          IF ( Wp(kr).NE.0. ) THEN
            Hpp(kr) = Hpp(kr)/Wp(kr)
            Am(kr) = Wp(kr)/Am(kr)
            IF ( Rh(kr).NE.0. ) Rh(kr) = Wp(kr)/Rh(kr)
            DO j = 1,Nlm
              B0p(j,kr) = B0p(j,kr)/Wp(kr)
              IF ( X(j).LT.0. ) Vmin(kr) = Vmin(kr) + B0p(j,kr)*X(j)
              IF ( X(j).GT.0. ) Vpls(kr) = Vpls(kr) + B0p(j,kr)*X(j)
            ENDDO
            IF ( .NOT.Moles ) THEN
              DO n = 1,Nreac
                IF ( Fox(n)(:1).NE.'O'.OR.kr.NE.2 ) THEN
                  IF ( Fox(n)(:1).EQ.'O'.OR.kr.NE.1 ) Pecwt(n)
     &                 = Pecwt(n)/Wp(kr)
                ENDIF
              ENDDO
            ENDIF
          ENDIF
        ENDDO
        IF ( .NOT.Short ) THEN
          IF ( Moles ) THEN
            WRITE (IOOUT,99008) ' MOLES '
          ELSE
            WRITE (IOOUT,99008) 'WT.FRAC'
          ENDIF
          DO n = 1,Nreac
            WRITE (IOOUT,99009) Fox(n),Rname(n),Pecwt(n),Enth(n),
     &           Rtemp(n),Dens(n),(Ratom(n,i),Rnum(n,i),i=1,Nfla(n)) 
          ENDDO
        ENDIF
      ENDIF
 200  RETURN
99001 FORMAT (/' REACTANT ',A15,'HAS BEEN DEFINED FOR THE TEMPERATURE', 
     &  F8.2,'K ONLY.'/' YOUR TEMPERATURE ASSIGNMENT',F8.2,
     &  ' IS MORE THAN 10 K FROM THIS VALUE. (REACT)')
99002 FORMAT (/' NOTE! REACTANT ',A15,'HAS BEEN DEFINED FOR ',
     &  'TEMPERATURE',F8.2,'K ONLY.'/' YOUR TEMPERATURE ASSIGNMENT',
     &  F8.2,' IS NOT = BUT <10 K FROM THIS VALUE. (REACT)')
99003 FORMAT (/' NOTE: ',A15,' IS EITHER NOT IN thermo.lib OR THE',
     &        ' TEMPERATURE ',/,
     &        ' IS OUT OF RANGE FOR THIS SPECIES (REACT)')
99004 FORMAT (/' TEMPERATURE MISSING FOR REACTANT NO.',I2,'(REACT)')
99005 FORMAT (/1x,a2,' NOT FOUND IN BLOCKDATA (REACT)')
99006 FORMAT (/' WARNING!!  AMOUNT MISSING FOR REACTANT',I3,'.',
     &        /' PROGRAM SETS WEIGHT PERCENT = 100. (REACT)')
99007 FORMAT (/' AMOUNT MISSING FOR REACTANT NO.',I2,'(REACT)')
99008 FORMAT (/4x,'REACTANT',10x,A7,3X,'(ENERGY/R),K',3X,
     &        'TEMP,K  DENSITY'/,8x,'EXPLODED FORMULA')
99009 FORMAT (1x,a1,': ',a15,f10.6,e15.6,f9.2,f8.4,/8x,5(2x,a2,f8.5))
99010 FORMAT (/' YOUR ASSIGNED TEMPERATURE',F8.2,'K FOR ',A15,/,
     & 'IS OUTSIDE ITS TEMPERATURE RANGE',F8.2,' TO',F9.2,'K (REACT)')
      END
