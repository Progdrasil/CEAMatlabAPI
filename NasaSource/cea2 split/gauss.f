      SUBROUTINE GAUSS
C***********************************************************************
C SOLVE ANY LINEAR SET OF UP TO MAXMAT EQUATIONS
C NUMBER OF EQUATIONS = IMAT
C***********************************************************************
      IMPLICIT NONE
      INCLUDE 'cea.inc'
C LOCAL VARIABLES
      INTEGER i,imatp1,j,k,nn,nnp1
      REAL*8 bigno,coefx(50),tmp
      REAL*8 DABS,DMAX1
      SAVE coefx,i,imatp1,j,k,nn,nnp1,tmp
C
      DATA bigno/1.E+25/
C BEGIN ELIMINATION OF NNTH VARIABLE
      imatp1 = Imat + 1
      DO nn = 1,Imat
        IF ( nn.NE.Imat ) THEN
C SEARCH FOR MAXIMUM COEFFICIENT IN EACH ROW
          nnp1 = nn + 1
          DO i = nn,Imat
            coefx(i) = bigno
            IF ( G(i,nn).NE.0. ) THEN
              coefx(i) = 0.
              DO j = nnp1,imatp1
                coefx(i) = DMAX1(coefx(i),DABS(G(i,j)))
              ENDDO
              tmp = DABS(G(i,nn))
              IF ( bigno*tmp.GT.coefx(i) ) THEN
                coefx(i) = coefx(i)/tmp
              ELSE
                coefx(i) = bigno
              ENDIF
            ENDIF
          ENDDO
C LOCATE ROW WITH SMALLEST MAXIMUM COEFFICIENT
          tmp = bigno
          i = 0
          DO j = nn,Imat
            IF ( coefx(j).LT.tmp ) THEN
              tmp = coefx(j)
              i = j
            ENDIF
          ENDDO
          IF ( i.EQ.0 ) THEN
            Msing = nn
            GOTO 99999
C INDEX I LOCATES EQUATION TO BE USED FOR ELIMINATING THE NTH
C VARIABLE FROM THE REMAINING EQUATIONS
C INTERCHANGE EQUATIONS I AND NN
          ELSEIF ( nn.NE.i ) THEN
            DO j = nn,imatp1
              tmp = G(i,j)
              G(i,j) = G(nn,j)
              G(nn,j) = tmp
            ENDDO
          ENDIF
        ELSEIF ( G(nn,nn).EQ.0 ) THEN
          Msing = nn
          GOTO 99999
        ENDIF
C DIVIDE NTH ROW BY NTH DIAGONAL ELEMENT AND ELIMINATE THE NTH
C VARIABLE FROM THE REMAINING EQUATIONS
        k = nn + 1
        tmp = G(nn,nn)
        IF ( tmp.EQ.0. ) THEN
          Msing = nn
          GOTO 99999
        ELSE
          DO j = k,imatp1
            G(nn,j) = G(nn,j)/tmp
          ENDDO
          IF ( k.NE.imatp1 ) THEN
            DO i = k,Imat
CDIR$ IVDEP
              DO j = k,imatp1
                G(i,j) = G(i,j) - G(i,nn)*G(nn,j)
              ENDDO
            ENDDO
          ENDIF
        ENDIF
      ENDDO
C BACKSOLVE FOR THE VARIABLES
      k = Imat
 100  j = k + 1
      X(k) = 0.0D0
      tmp = 0.0
      IF ( Imat.GE.j ) THEN
        DO i = j,Imat
          tmp = tmp + G(k,i)*X(i)
        ENDDO
      ENDIF
      X(k) = G(k,imatp1) - tmp
      k = k - 1
      IF ( k.NE.0 ) GOTO 100
99999 END
      SUBROUTINE HCALC
C***********************************************************************
C CALCULATE PROPERTIES FOR TOTAL REACTANT USING THERMO DATA FOR
C ONE OR MORE REACTANTS. USED ONLY FOR SHOCK AND DETON PROBLEMS.
C***********************************************************************
      IMPLICIT NONE
      INCLUDE 'cea.inc'
C LOCAL VARIABLES
      CHARACTER*6 date(MAXNGC)
      CHARACTER*2 el(5)
      CHARACTER*15 sub
      INTEGER i,icf,ifaz,itot,j,k,l,m,n,nall,nint,ntgas,ntot
      REAL*8 bb(5),enj,er,sj,t1,t2,tem,thermo(9,3),tsave
      REAL*8 DLOG
      SAVE bb,date,el,enj,er,i,icf,ifaz,itot,j,k,l,m,n,nall,nint,ntgas,
     &  ntot,sj,sub,t1,t2,tem,thermo,tsave
C
      tsave = Tt
      Tm = 0.
      IF ( Pp.GT.0. ) Tm = DLOG(Pp*Wmix)
      Ssum(Npt) = 0.
      Hpp(1) = 0.
      Hpp(2) = 0.
      Hsub0 = 0.
      Cpmix = 0.
      tem = (1.+Oxfl)
C LOOP ON REACTANTS.
C IF OXIDANT, K=1
C IF FUEL, K=2
      Nspr = Nspx
      DO n = 1,Nreac
        k = 2
        IF ( Fox(n)(:1).EQ.'O'.OR.Fox(n)(:1).EQ.'o' ) k = 1
        IF ( Tt.EQ.0. ) Tt = Rtemp(n)
        j = Jray(n)
        IF ( j.EQ.0 ) THEN
C SEARCH FOR REACTANT IN STORED THERMO SPECIES. STORE INDEX IN JRAY(N).
          ifaz = 0
          DO j = 1,Ngc
            IF ( Rname(n).EQ.Prod(j).OR.'*'//Rname(n).EQ.Prod(j) ) THEN
              Jray(n) = j
              IF ( j.GT.Ng ) THEN
                WRITE (IOOUT,99001)
                GOTO 20
              ENDIF
              GOTO 50
            ENDIF
          ENDDO
C SEARCH THERMO.LIB FOR SPECIES.
          REWIND IOTHM
          READ (IOTHM) Tg,ntgas,ntot,nall
          Nspr = Nspr + 1
          DO itot = 1,nall
            IF ( itot.LE.ntot ) THEN
              icf = 3
              IF ( itot.GT.ntgas ) icf = 1
              READ (IOTHM) sub,nint,date(Nspr),(el(j),bb(j),j=1,5),ifaz,
     &                     t1,t2,Mw(Nspr),((thermo(l,m),l=1,9),m=1,icf)
            ELSE
              READ (IOTHM) sub,nint,date(Nspr),(el(j),bb(j),j=1,5),ifaz,
     &                     t1,t2,Mw(Nspr),er
              IF ( nint.NE.0 ) THEN
                READ (IOTHM) ((thermo(i,j),i=1,9),j=1,nint)
                icf = nint
              ENDIF
            ENDIF
            IF ( sub.EQ.Rname(n).OR.sub.EQ.'*'//Rname(n) ) THEN
              IF ( ifaz.LE.0.AND.nint.GT.0 ) THEN
                DO j = 1,5
                  IF ( bb(j).EQ.0. ) GOTO 2
                  Nfla(n) = j
                  Ratom(n,j) = el(j)
                  Rnum(n,j) = bb(j)
                ENDDO
 2              Jray(n) = Nspr
                j = Nspr
                DO l = 1,icf
                  DO m = 1,9
                    Coef(j,m,l) = thermo(m,l)
                  ENDDO
                ENDDO
                GOTO 50
              ELSE
                IF ( ifaz.GT.0 ) WRITE (IOOUT,99001)
                IF ( nint.EQ.0 ) WRITE (IOOUT,99002) Rname(n)
                GOTO 20
              ENDIF
            ENDIF
          ENDDO
          Nspr = Nspr - 1
          WRITE (IOOUT,99003) Rname(n)
          Energy(n) = ' '
 20       Tt = 0.
          Cpmix = 0.
          GOTO 100
        ENDIF
C CALCULATE EN FOR REACTANT AND CALCULATE PROPERTIES.
 50     IF ( Moles ) enj = Pecwt(n)/Wp(k)
        IF ( .NOT.Moles ) enj = Pecwt(n)/Rmw(n)
        enj = enj/tem
        IF ( k.EQ.1 ) enj = enj*Oxfl
        Tln = DLOG(Tt)
        En(j,Npt) = enj
        l = 1
        IF ( ifaz.LE.0 ) THEN
          IF ( Tt.GT.Tg(2) ) l = 2
          IF ( Tt.GT.Tg(3).AND.ifaz.LT.0 ) l = 3
        ENDIF
        S(j) = ((((Coef(j,7,l)/4.)*Tt+Coef(j,6,l)/3.)*Tt+Coef(j,5,l)/2.)
     &         *Tt+Coef(j,4,l))*Tt - (Coef(j,1,l)*.5D0/Tt+Coef(j,2,l))
     &         /Tt + Coef(j,3,l)*Tln + Coef(j,9,l)
        H0(j) = ((((Coef(j,7,l)/5.)*Tt+Coef(j,6,l)/4.)*Tt+Coef(j,5,l)/3.
     &          )*Tt+Coef(j,4,l)/2.)
     &          *Tt - (Coef(j,1,l)/Tt-Coef(j,2,l)*Tln-Coef(j,8,l))
     &          /Tt + Coef(j,3,l)
        Cp(j) = (((Coef(j,7,l)*Tt+Coef(j,6,l))*Tt+Coef(j,5,l))
     &          *Tt+Coef(j,4,l))*Tt + (Coef(j,1,l)/Tt+Coef(j,2,l))
     &          /Tt + Coef(j,3,l)
        IF ( H0(j).GT.-.01.AND.H0(j).LT..01 ) H0(j) = 0.
C ADD CONTRIBUTION TO CP, H, AND S OF TOTAL REACTANT.
        Cpmix = Cpmix + Cp(j)*enj
C FOR CONDENSED SPECIES:  SJ = S(J)
        sj = S(j) - DLOG(enj) - Tm
        Ssum(Npt) = Ssum(Npt) + enj*sj
        er = H0(j)*enj*Tt
        Hsub0 = Hsub0 + er
        Hpp(k) = Hpp(k) + er
      ENDDO
      IF ( tsave.NE.0. ) Tt = tsave
 100  RETURN
99001 FORMAT (/' REACTANTS MUST BE GASEOUS FOR THIS PROBLEM (HCALC)')
99002 FORMAT (/' COEFFICIENTS FOR ',A15,' ARE NOT AVAILABLE (HCALC)')
99003 FORMAT (/' ERROR IN DATA FOR ',A15,' CHECK NAME AND TEMPERATURE',
     &        ' RANGE IN',/,' thermo.inp (HCALC)')
      END
