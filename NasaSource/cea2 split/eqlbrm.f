      SUBROUTINE EQLBRM
C***********************************************************************
C CALCULATE EQUILIBRIUM COMPOSITION AND PROPERTIES.
C***********************************************************************
      IMPLICIT NONE
      INCLUDE 'cea.inc'
C LOCAL VARIABLES
      CHARACTER*12 ae,cmp(MAXEL)
      CHARACTER*16 amb
      LOGICAL cpcalc,i2many,newcom,reduce
      INTEGER i,il,ilamb,ilamb1,inc,ipr,iq2,iter,ix,ixsing,iz,j,ja,jb,
     &        jbx,jc,jcondi,jcons,jdelg,jex,jj,jkg,jneg,jsw,k,kc,kg,kk,
     &        kmat,kneg,l,lc,lcs(MAXEL),le,lelim,lk,ll,lncvg,ls,lsing,
     &        lz,maxitn,ncvg,njc,nn,numb
      INTEGER IABS
      REAL*8 aa,ambda,ambda1,bigen,bigneg,delg,dlnt,dpie,ensol,esize,
     &       gap,gasfrc,pie,pisave(MAXMAT-2),siz9,sizeg,smalno,smnol,
     &       sum,sum1,szgj,tem,tmelt,tsize,ween,xi,xln,xsize,xx(MAXMAT)
      REAL*8 DABS,DEXP,DLOG,DMAX1
      SAVE aa,ae,amb,ambda,ambda1,bigen,bigneg,cmp,cpcalc,delg,dlnt,
     &  dpie,ensol,esize,gap,gasfrc,i,i2many,il,ilamb,ilamb1,inc,ipr,
     &  iq2,iter,ix,ixsing,iz,j,ja,jb,jbx,jc,jcondi,jcons,jdelg,jex,jj,
     &  jkg,jneg,jsw,k,kc,kg,kk,kmat,kneg,l,lc,lcs,le,lelim,lk,ll,lncvg,
     &  ls,lsing,lz,maxitn,ncvg,newcom,njc,nn,numb,pie,pisave,reduce,
     &  siz9,sizeg,sum,sum1,szgj,tem,tmelt,tsize,ween,xi,xln,xsize,xx
C
      DATA smalno/1.E-6/,smnol/ - 13.815511/
      ixsing = 0
      lsing = 0
      jsw = 0
      jdelg = 0
      maxitn = 50
      ncvg = 0
      lncvg = 3*Nlm
      reduce = .FALSE.
      siz9 = Size - 9.2103404D0
      tsize = Size
      xsize = Size + 6.90775528D0
      IF ( Trace.NE.0. ) THEN
        maxitn = maxitn + Ngc/2
        xsize = -DLOG(Trace)
        IF ( xsize.LT.Size ) xsize = Size + .1
      ENDIF
      IF ( xsize.GT.80. ) xsize = 80.D0
      esize = MIN(80.D0,xsize+6.90775528D0)
      jcons = 0
      pie = 0.
      i2many = .FALSE.
      Pderiv = .FALSE.
      Convg = .FALSE.
      numb = 0
      cpcalc = .TRUE.
      IF ( Tp ) cpcalc = .FALSE.
      IF ( Tt.NE.0.D0 ) THEN
        IF ( Npr.EQ.0.OR.(Tt.NE.T(1).AND..NOT.Tp) ) GOTO 400
        k = 1
      ELSE
        Tt = 3800.D0
        IF ( Npr.EQ.0 ) GOTO 400
        k = 1
      ENDIF
 100  j = Jcond(k)
      jc = j - Ng
      kg = -Ifz(jc)
      DO i = 1,9
        kg = kg + 1
        kc = jc + kg
        IF ( Tt.LE.Temp(2,kc) ) THEN
          IF ( kg.NE.0 ) THEN
            Jcond(k) = j + kg
            En(j+kg,Npt) = En(j,Npt)
            En(j,Npt) = 0.
            IF ( Prod(j).NE.Prod(j+kg).AND..NOT.Short ) 
     &           WRITE (IOOUT,99023) Prod(j),Prod(j+kg)
          ENDIF
          GOTO 300
        ELSEIF ( kc.GE.Nc.OR.Ifz(kc+1).LE.Ifz(kc) ) THEN
          GOTO 200
        ENDIF
      ENDDO
 200  IF ( .NOT.Tp ) THEN
        Tt = Temp(2,kc) - 10.D0
        k = 1
        GOTO 100
      ENDIF
      WRITE (IOOUT,99028) Prod(j)
      En(j,Npt) = 0.D0
      Enln(j) = 0.D0
      Deln(j) = 0.D0
      DO i = k,Npr
        Jcond(i) = Jcond(i+1)
      ENDDO
      Npr = Npr - 1
 300  k = k + 1
      IF ( k.LE.Npr ) GOTO 100
 400  Tln = DLOG(Tt)
      IF ( Vol ) Pp = Rr*Enn*Tt/Vv
      CALL CPHS
      Tm = DLOG(Pp/Enn)
      le = Nlm
      IF ( Lsave.NE.0.AND.Nlm.NE.Lsave ) THEN
        tem = EXP(-tsize)
        DO i = Lsave + 1,Nlm
          DO j = 1,Ng
            IF ( A(i,j).NE.0. ) THEN
              En(j,Npt) = tem
              Enln(j) = -tsize
            ENDIF
          ENDDO
        ENDDO
      ENDIF
      ls = Nlm
      lelim = 0
      lz = ls
      IF ( Ions ) lz = ls - 1
      IF ( Npt.EQ.1.AND..NOT.Shock.AND..NOT.Short ) WRITE (IOOUT,99001)
     &     (Elmt(i),i=1,Nlm)
      IF ( Debug(Npt) ) THEN
        DO i = 1,Nlm
          cmp(i) = Elmt(i)
        ENDDO
      ENDIF
C BEGIN ITERATION
 500  IF ( cpcalc ) THEN
        Cpsum = 0.D0
        DO j = 1,Ng
          Cpsum = Cpsum + En(j,Npt)*Cp(j)
        ENDDO
        IF ( Npr.NE.0 ) THEN
          DO k = 1,Npr
            j = Jcond(k)
            Cpsum = Cpsum + En(j,Npt)*Cp(j)
          ENDDO
          cpcalc = .FALSE.
        ENDIF
      ENDIF
      numb = numb + 1
      CALL MATRIX
      iq2 = Iq1 + 1
      IF ( Convg ) Imat = Imat - 1
      IF ( Debug(Npt) ) THEN
        IF ( .NOT.Convg ) THEN
          WRITE (IOOUT,99004) numb
        ELSE
          IF ( .NOT.Pderiv ) WRITE (IOOUT,99002)
          IF ( Pderiv ) WRITE (IOOUT,99003)
        ENDIF
        kmat = Imat + 1
        DO i = 1,Imat
          WRITE (IOOUT,99006) (G(i,k),k=1,kmat)
        ENDDO
      ENDIF
      Msing = 0
      CALL GAUSS
      IF ( Msing.EQ.0 ) THEN
        IF ( Debug(Npt) ) THEN
          WRITE (IOOUT,99005) (cmp(k),k=1,le)
          WRITE (IOOUT,99006) (X(i),i=1,Imat)
        ENDIF
        IF ( .NOT.Convg ) THEN
C OBTAIN CORRECTIONS TO THE ESTIMATES
          IF ( Vol ) X(iq2) = X(Iq1)
          IF ( Tp ) X(iq2) = 0.
          dlnt = X(iq2)
          sum = X(Iq1)
          IF ( Vol ) THEN
            X(Iq1) = 0.
            sum = -dlnt
          ENDIF
          DO 520 j = 1,Ng
            IF ( lelim.NE.0 ) THEN
              Deln(j) = 0.
              DO i = lelim,ls
                IF ( A(i,j).NE.0. ) GOTO 520
              ENDDO
            ENDIF
            Deln(j) = -Mu(j) + H0(j)*dlnt + sum
            DO k = 1,Nlm
              Deln(j) = Deln(j) + A(k,j)*X(k)
            ENDDO
            IF ( pie.NE.0. ) Deln(j) = Deln(j) + A(ls,j)*pie
 520      CONTINUE
          IF ( Npr.NE.0 ) THEN
            DO k = 1,Npr
              j = Jcond(k)
              kk = Nlm + k
              Deln(j) = X(kk)
            ENDDO
          ENDIF
C CALCULATE CONTROL FACTOR,AMBDA
          ambda = 1.D0
          ambda1 = 1.D0
          ilamb = 0
          ilamb1 = 0
          sum = DMAX1(DABS(X(Iq1)),DABS(dlnt))
          sum = sum*5.
          DO j = 1,Ng
            IF ( Deln(j).GT.0. ) THEN
              IF ( (Enln(j)-Ennl+Size).LE.0. ) THEN
                sum1 = DABS(Deln(j)-X(Iq1))
                IF ( sum1.GE.siz9 ) THEN
                  sum1 = DABS(-9.2103404D0-Enln(j)+Ennl)/sum1
                  IF ( sum1.LT.ambda1 ) THEN
                    ambda1 = sum1
                    ilamb1 = j
                  ENDIF
                ENDIF
              ELSEIF ( Deln(j).GT.sum ) THEN
                sum = Deln(j)
                ilamb = j
              ENDIF
            ENDIF
          ENDDO
          IF ( sum.GT.2.D0 ) ambda = 2.D0/sum
          IF ( ambda1.LE.ambda ) THEN
            ambda = ambda1
            ilamb = ilamb1
          ENDIF
          IF ( Debug(Npt) ) THEN
C INTERMEDIATE OUTPUT
            WRITE (IOOUT,99011) Tt,Enn,Ennl,Pp,Tm,ambda
            IF ( ambda.NE.1.D0 ) THEN
              amb = 'ENN'
              IF ( DABS(X(iq2)).GT.DABS(X(Iq1)) ) amb = 'TEMP'
              IF ( ilamb.NE.0 ) amb = Prod(ilamb)
              WRITE (IOOUT,99012) amb
            ENDIF
            IF ( Vol ) WRITE (IOOUT,99013) Vv*.001D0
            WRITE (IOOUT,99014)
            DO j = 1,Ngc
              WRITE (IOOUT,99015) Prod(j),En(j,Npt),Enln(j),Deln(j),
     &                       H0(j),S(j),H0(j) - S(j),Mu(j)
            ENDDO
          ENDIF
C APPLY CORRECTIONS TO ESTIMATES
          Totn(Npt) = 0.D0
          DO j = 1,Ng
            Enln(j) = Enln(j) + ambda*Deln(j)
          ENDDO
          DO 540 j = 1,Ng
            En(j,Npt) = 0.
            IF ( lelim.NE.0 ) THEN
              DO i = lelim,ls
                IF ( A(i,j).NE.0. ) GOTO 540
              ENDDO
            ENDIF
            IF ( (Enln(j)-Ennl+tsize).GT.0. ) THEN
              En(j,Npt) = DEXP(Enln(j))
              Totn(Npt) = Totn(Npt) + En(j,Npt)
            ENDIF
 540      CONTINUE
          IF ( Ions.AND.Elmt(Nlm).EQ.'E' ) THEN
            DO j = 1,Ng
              IF ( A(ls,j).NE.0..AND.En(j,Npt).EQ.0. ) THEN
                IF ( (Enln(j)-Ennl+esize).GT.0. ) THEN
                  En(j,Npt) = DEXP(Enln(j))
                  Totn(Npt) = Totn(Npt) + En(j,Npt)
                ENDIF
              ENDIF
            ENDDO
          ENDIF
          Sumn = Totn(Npt)
          IF ( Npr.NE.0 ) THEN
            DO k = 1,Npr
              j = Jcond(k)
              En(j,Npt) = En(j,Npt) + ambda*Deln(j)
              Totn(Npt) = Totn(Npt) + En(j,Npt)
            ENDDO
          ENDIF
          IF ( .NOT.Tp ) THEN
            Tln = Tln + ambda*dlnt
            Tt = DEXP(Tln)
            cpcalc = .TRUE.
            CALL CPHS
          ENDIF
          IF ( Vol ) THEN
            Enn = Sumn
            Ennl = DLOG(Enn)
            IF ( Vol ) Pp = Rr*Tt*Enn/Vv
          ELSE
            Ennl = Ennl + ambda*X(Iq1)
            Enn = DEXP(Ennl)
          ENDIF
          Tm = DLOG(Pp/Enn)
          IF ( Elmt(Nlm).EQ.'E' ) THEN
C CHECK ON REMOVING IONS
            DO j = 1,Ngc
              IF ( A(Nlm,j).NE.0. ) THEN
                IF ( En(j,Npt).GT.0. ) GOTO 560
              ENDIF
            ENDDO
            pie = X(Nlm)
            lelim = Nlm
            Nlm = Nlm - 1
            GOTO 500
          ENDIF
C TEST FOR CONVERGENCE
 560      IF ( numb.GT.maxitn ) THEN
            WRITE (IOOUT,99019) maxitn,Npt
            IF ( Nc.EQ.0.OR.i2many ) GOTO 1500
            i2many = .TRUE.
            IF ( .NOT.Hp.OR.Npt.NE.1.OR.Tt.GT.100. ) THEN
              IF ( Npr.NE.1.OR.Enn.GT.1.E-4 ) GOTO 1500
C HIGH TEMPERATURE, INCLUDED CONDENSED CONDITION
              WRITE (IOOUT,99020)
              Enn = .1
              Ennl = -2.3025851
              Sumn = Enn
              xi = Ng
              xi = Enn/xi
              xln = DLOG(xi)
              DO j = 1,Ng
                En(j,Npt) = xi
                Enln(j) = xln
              ENDDO
              j = Jcond(1)
              k = 1
              GOTO 1000
            ELSE
              WRITE (IOOUT,99008)
              GOTO 1500
            ENDIF
          ELSE
            sum = (X(Iq1)*Enn/Totn(Npt))
            IF ( DABS(sum).GT.0.5E-5 ) GOTO 500
            DO j = 1,Ng
              IF ( DABS(Deln(j))*En(j,Npt)/Totn(Npt).GT.0.5D-5 )
     &             GOTO 500
            ENDDO
            IF ( DABS(dlnt).GT.1.D-04 ) GOTO 500
            IF ( Npr.NE.0 ) THEN
              DO k = 1,Npr
                j = Jcond(k)
                IF ( DABS(Deln(j)/Totn(Npt)).GT.0.5D-5 ) GOTO 500
                IF ( En(j,Npt).LT.0. ) GOTO 700
              ENDDO
            ENDIF
            le = Nlm
            DO i = 1,Nlm
              IF ( DABS(B0(i)).GE.1.D-06 ) THEN
                sum = 0.
                DO j = 1,Ngc
                  sum = sum + En(j,Npt)*A(i,j)
                ENDDO
                IF ( DABS(B0(i)-sum).GT.Bcheck ) GOTO 500
              ENDIF
            ENDDO
            IF ( Trace.NE.0. ) THEN
              tsize = xsize
              tem = 1.
              IF ( numb.NE.1 ) THEN
                lk = lz
                IF ( Nlm.LT.lz ) lk = Nlm
                DO i = 1,lk
                  IF ( i.NE.lsing ) THEN
                    tem = 0.
                    IF ( X(i).NE.0. ) THEN
                      tem = DABS((pisave(i)-X(i))/X(i))
                      IF ( tem.GT..001 ) GOTO 565
                    ENDIF
                  ENDIF
                ENDDO
              ENDIF
 565          DO i = 1,Nlm
                pisave(i) = X(i)
              ENDDO
              IF ( tem.GT..001 ) GOTO 500
              IF ( Ions ) THEN
C CHECK ON ELECTRON BALANCE
                iter = 1
                IF ( pie.NE.0. ) THEN
                  le = Nlm + 1
                  X(le) = pie
                ENDIF
 566            sum1 = 0.
                sum = 0.
                pie = X(le)
                DO j = 1,Ng
                  IF ( A(ls,j).NE.0. ) THEN
                    En(j,Npt) = 0.
                    tem = 0.
                    IF ( Enln(j).GT.-87. ) tem = DEXP(Enln(j))
                    IF ( (Enln(j)-Ennl+tsize).GT.0..AND.Elmt(Nlm)
     &                   .EQ.'E' ) THEN
                      pie = 0.
                      En(j,Npt) = tem
                    ENDIF
                    aa = A(ls,j)*tem
                    sum = sum + aa
                    sum1 = sum1 + aa*A(ls,j)
                  ENDIF
                ENDDO
                IF ( sum1.NE.0. ) THEN
                  dpie = -sum/sum1
                  DO j = 1,Ng
                    IF ( A(ls,j).NE.0. ) Enln(j) = Enln(j) + A(ls,j)
     &                   *dpie
                  ENDDO
                  IF ( Debug(Npt) ) WRITE (IOOUT,99016) iter,dpie
                  IF ( DABS(dpie).GT..0001 ) THEN
                    X(le) = X(le) + dpie
                    iter = iter + 1
                    IF ( iter.LE.80 ) GOTO 566
                    WRITE (IOOUT,99017)
                    GOTO 1500
                  ELSEIF ( Elmt(Nlm).EQ.'E'.AND.pie.NE.0. ) THEN
                    Nlm = Nlm - 1
                    newcom = .TRUE.
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ELSEIF ( .NOT.Pderiv ) THEN
C TEMPERATURE DERIVATIVES--CONVG=T, PDERIV=F
          Dlvtp(Npt) = 1. - X(Iq1)
          Cpr(Npt) = G(iq2,iq2)
          DO j = 1,Iq1
            Cpr(Npt) = Cpr(Npt) - G(iq2,j)*X(j)
          ENDDO
C PRESSURE DERIVATIVE--CONVG=T, PDERIV=T
          Pderiv = .TRUE.
          GOTO 500
        ELSE
          Dlvpt(Npt) = -1. + X(Iq1)
          IF ( Jliq.EQ.0 ) THEN
            Gammas(Npt) = -1./(Dlvpt(Npt)+(Dlvtp(Npt)**2)*Enn/Cpr(Npt))
          ELSE
            En(Jsol,Npt) = ensol
            Hsum(Npt) = Hsum(Npt) + En(Jliq,Npt)*(H0(Jliq)-H0(Jsol))
            Gammas(Npt) = -1./Dlvpt(Npt)
            Npr = Npr + 1
            Jcond(Npr) = Jliq
          ENDIF
          GOTO 1400
        ENDIF
C SINGULAR MATRIX
      ELSE
        IF ( Convg ) THEN
          WRITE (IOOUT,99007)
          Dlvpt(Npt) = -1.
          Dlvtp(Npt) = 1.
          Cpr(Npt) = Cpsum
          Gammas(Npt) = -1./(Dlvpt(Npt)+(Dlvtp(Npt)**2)*Enn/Cpr(Npt))
          GOTO 1400
        ELSE
          WRITE (IOOUT,99009) numb,Msing
          lsing = Msing
          ixsing = ixsing + 1
          IF ( ixsing.LE.8 ) THEN
            xsize = 80.
            tsize = xsize
            IF ( Msing.GT.Nlm.AND.numb.LT.1.AND.Npr.GT.1.AND.
     &           jdelg.GT.0 ) THEN
              ween = 1000.
              j = 0
              DO 570 i = 1,Npr
                jcondi = Jcond(i)
                IF ( jcondi.NE.jdelg ) THEN
                  DO ll = 1,Nlm
                    IF ( A(ll,jdelg).NE.0.AND.A(ll,jcondi).NE.0. ) THEN
                      IF ( En(jcondi,Npt).LE.ween ) THEN
                        ween = En(jcondi,Npt)
                        j = jcondi
                        k = i
                      ENDIF
                      GOTO 570
                    ENDIF
                  ENDDO
                ENDIF
 570          CONTINUE
              IF ( j.GT.0 ) THEN
                WRITE (IOOUT,99020)
                GOTO 1000
              ENDIF
            ELSEIF ( .NOT.Hp.OR.Npt.NE.1.OR.Nc.EQ.0.OR.Tt.GT.100. ) THEN
              IF ( ixsing.GE.3 ) THEN
                IF ( Msing.LT.Iq1 ) THEN
                  IF ( reduce.AND.Msing.LE.Nlm ) THEN
                    IF ( Nlm.LT.lelim ) GOTO 1500
                    WRITE (IOOUT,99010) Npt,Elmt(Nlm)
                    Nlm = Nlm - 1
                    GOTO 500
                  ELSEIF ( Msing.LE.Nlm ) THEN
C FIND NEW COMPONENTS
                    IF ( .NOT.Ions ) GOTO 1100
                    IF ( Elmt(Nlm).NE.'E' ) GOTO 1100
                    DO j = 1,Ng
                      IF ( A(Nlm,j).NE.0. ) En(j,Npt) = 0.D0
                    ENDDO
                    pie = X(Nlm)
                    Nlm = Nlm - 1
                    IF ( Msing.GT.Nlm ) GOTO 500
                    GOTO 1100
                  ELSE
C REMOVE CONDENSED SPECIES TO CORRECT SINGULARITY
                    k = Msing - Nlm
                    j = Jcond(k)
                    IF ( j.NE.jcons ) THEN
                      jcons = j
                      GOTO 1000
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
              DO 575 jj = 1,Ng
                IF ( Ions ) THEN
                  IF ( Elmt(Nlm).NE.'E' ) THEN
                    IF ( A(ls,jj).NE.0. ) GOTO 575
                  ENDIF
                ENDIF
                IF ( En(jj,Npt).EQ.0. ) THEN
                  En(jj,Npt) = smalno
                  Enln(jj) = smnol
                ENDIF
 575          CONTINUE
              GOTO 500
            ELSE
              WRITE (IOOUT,99008)
            ENDIF
          ENDIF
        ENDIF
        GOTO 1500
      ENDIF
C CALCULATE ENTROPY, CHECK ON DELTA S FOR SP PROBLEMS
 600  Ssum(Npt) = 0.
      DO j = 1,Ng
        Ssum(Npt) = Ssum(Npt) + En(j,Npt)*(S(j)-Enln(j)-Tm)
      ENDDO
      IF ( Npr.GT.0 ) THEN
        DO k = 1,Npr
          j = Jcond(k)
          Ssum(Npt) = Ssum(Npt) + En(j,Npt)*S(j)
        ENDDO
      ENDIF
      IF ( .NOT.Sp ) THEN
        Convg = .TRUE.
      ELSE
        tem = Ssum(Npt) - S0
        IF ( DABS(tem).GT..0005 ) GOTO 500
        IF ( Debug(Npt) ) WRITE (IOOUT,99018) tem
        Convg = .TRUE.
      ENDIF
C CONVERGENCE TESTS ARE SATISFIED, TEST CONDENSED SPECIES.
 700  ncvg = ncvg + 1
      IF ( ncvg.GT.lncvg ) THEN
C ERROR, SET TT=0
        WRITE (IOOUT,99034) lncvg
        GOTO 1500
      ELSE
        IF ( .NOT.Shock ) THEN
          DO il = 1,le
            xx(il) = X(il)
          ENDDO
          IF ( .NOT.Short ) THEN
            IF ( newcom ) WRITE (IOOUT,99021) (cmp(k),k=1,le)
            WRITE (IOOUT,99022) Npt,numb,Tt,(xx(il),il=1,le)
          ENDIF
          IF ( .NOT.Tp.AND.Npr.EQ.0.AND.Tt.LE.Tg(1)*.2D0 ) THEN
            WRITE (IOOUT,99008)
            GOTO 1500
          ENDIF
          newcom = .FALSE.
        ENDIF
        IF ( Npr.NE.0 ) THEN
          bigneg = 0.
          jneg = 0
          DO k = 1,Npr
            j = Jcond(k)
            IF ( En(j,Npt)*Cp(j).LE.bigneg ) THEN
              bigneg = En(j,Npt)*Cp(j)
              jneg = j
              kneg = k
            ENDIF
          ENDDO
          IF ( jneg.NE.0 ) THEN
            j = jneg
            k = kneg
            IF ( j.EQ.Jsol.OR.j.EQ.Jliq ) THEN
              Jsol = 0
              Jliq = 0
            ENDIF
            GOTO 1000
          ENDIF
        ENDIF
        IF ( Ngc.NE.Ng.OR.Tp ) THEN
          Ng = Ngc
          CALL CPHS
          Ng = Ngp1 - 1
          cpcalc = .TRUE.
          IF ( Ngc.EQ.Ng ) GOTO 750
          CALL ALLCON
          IF ( Npr.NE.0.AND..NOT.Tp ) THEN
            gap = 50.
            DO 710 ipr = 1,Npr
              j = Jcond(ipr)
              IF ( j.NE.Jsol.AND.j.NE.Jliq ) THEN
                inc = j - Ng
                kg = -Ifz(inc)
                DO iz = 1,20
                  kg = kg + 1
                  kc = inc + kg
                  IF ( Tt.LE.Temp(2,kc) ) THEN
                    IF ( kg.NE.0 ) THEN
                      jkg = j + kg
                      IF ( IABS(kg).GT.1.OR.Prod(j).EQ.Prod(jkg) )
     &                     GOTO 740
                      IF ( jkg.EQ.jsw ) GOTO 720
                      IF ( Tt.LT.Temp(1,inc)-gap.OR.Tt.GT.Temp(2,inc)
     &                     +gap ) GOTO 740
                      GOTO 720
                    ENDIF
                    GOTO 710
                  ELSEIF ( Ifz(kc+1).LE.Ifz(kc) ) THEN
                    GOTO 710
                  ENDIF
                ENDDO
                IF ( Tt.GT.Temp(2,kc)*1.2D0 ) GOTO 1000
              ENDIF
 710        CONTINUE
          ENDIF
          sizeg = 0.
          szgj = 0.
          DO inc = 1,Nc
            j = inc + Ng
            IF ( Debug(Npt) ) WRITE (IOOUT,99024) Prod(j),Temp(1,inc),
     &                               Temp(2,inc),En(j,Npt)
            IF ( En(j,Npt).LE.0. ) THEN
              IF ( Tt.GT.Temp(1,inc).OR.Temp(1,inc).EQ.Tg(1) ) THEN
                IF ( Tt.LE.Temp(2,inc) ) THEN
                  sum = 0.
                  DO i = 1,Nlm
                    sum = sum + A(i,j)*X(i)
                  ENDDO
                  delg = (H0(j)-S(j)-sum)/Mw(j)
                  IF ( delg.LT.sizeg.AND.delg.LT.0. ) THEN
                    IF ( j.NE.jcons ) THEN
                      sizeg = delg
                      jdelg = j
                    ELSE
                      szgj = delg
                    ENDIF
                    ipr = ipr - 1
                  ENDIF
                  IF ( Debug(Npt) ) WRITE (IOOUT,99025) delg,sizeg
                ENDIF
              ENDIF
            ENDIF
          ENDDO
          IF ( sizeg.EQ.0..AND.szgj.EQ.0. ) GOTO 750
          IF ( sizeg.NE.0. ) THEN
            j = jdelg
            GOTO 800
          ELSE
            WRITE (IOOUT,99026) Prod(jcons)
            GOTO 1500
          ENDIF
 720      kk = MAX(0,kg)
          tmelt = Temp(kk+1,inc)
          Tt = tmelt
          Tln = DLOG(Tt)
          Jsol = MIN(j,jkg)
          Jliq = Jsol + 1
          En(jkg,Npt) = .5D0*En(j,Npt)
          En(j,Npt) = En(jkg,Npt)
          j = jkg
          GOTO 800
C WRONG PHASE INCLUDED FOR T INTERVAL, SWITCH EN
 740      En(jkg,Npt) = En(j,Npt)
          Jcond(ipr) = jkg
          En(j,Npt) = 0.
          jsw = j
          IF ( Prod(j).NE.Prod(jkg).AND..NOT.Short ) WRITE (IOOUT,99023)
     &         Prod(j),Prod(jkg)
          j = jkg
          GOTO 900
        ENDIF
C CONVERGED WITH NO CONDENSED CHANGES.  IF BOTH SOLID & LIQ PRESENT,
C TEMPORARILY REMOVE LIQ TO PREVENT SINGULAR DERIVATIVE MATRIX.
 750    Sumn = Enn
        IF ( Jsol.NE.0 ) THEN
          ensol = En(Jsol,Npt)
          En(Jsol,Npt) = En(Jsol,Npt) + En(Jliq,Npt)
          Dlvtp(Npt) = 0.
          Cpr(Npt) = 0.
          Gammas(Npt) = 0.
          Pderiv = .TRUE.
          DO k = 1,Npr
            IF ( Jcond(k).EQ.Jliq ) GOTO 760
          ENDDO
 760      DO i = k,Npr
            Jcond(i) = Jcond(i+1)
          ENDDO
          Npr = Npr - 1
        ENDIF
        GOTO 500
      ENDIF
C ADD CONDENSED SPECIES
 800  Npr = Npr + 1
      i = Npr
      DO ix = 2,Npr
        Jcond(i) = Jcond(i-1)
        i = i - 1
      ENDDO
      Jcond(1) = j
      IF ( .NOT.Short ) WRITE (IOOUT,99027) Prod(j)
 900  inc = j - Ng
      Convg = .FALSE.
      IF ( Tp ) cpcalc = .FALSE.
      numb = -1
      GOTO 500
C REMOVE CONDENSED SPECIES
 1000 En(j,Npt) = 0.D0
      Deln(j) = 0.D0
      Enln(j) = 0.D0
      DO i = k,Npr
        Jcond(i) = Jcond(i+1)
      ENDDO
      IF ( .NOT.Short ) WRITE (IOOUT,99028) Prod(j)
      Npr = Npr - 1
      DO i = 1,Nlm
        IF ( cmp(i).EQ.Prod(j) ) THEN
          numb = -1
          Convg = .FALSE.
          IF ( Tp ) cpcalc = .FALSE.
          GOTO 1100
        ENDIF
      ENDDO
      GOTO 900
 1100 newcom = .FALSE.
      nn = Nlm
      IF ( Elmt(Nlm).EQ.'E' ) nn = Nlm - 1
C FIND ORDER OF SPECIES FOR COMPONENTS - BIGGEST TO SMALLEST
      njc = 0
      DO lc = 1,nn
        lcs(lc) = 0
      ENDDO
 1200 bigen = -1.D-35
      DO j = 1,Ng
        IF ( En(j,Npt).GT.bigen ) THEN
          IF ( .NOT.Ions.OR.A(ls,j).EQ.0. ) THEN
            bigen = En(j,Npt)
            jbx = j
          ENDIF
        ENDIF
      ENDDO
      IF ( bigen.GT.0. ) THEN
        DO 1250 lc = 1,nn
          IF ( jbx.EQ.0 ) jbx = Jx(lc)
          IF ( A(lc,jbx).GT.smalno ) THEN
            IF ( njc.NE.0 ) THEN
              DO 1205 i = 1,njc
                l = lcs(i)
                IF ( l.EQ.lc ) GOTO 1250
                IF ( l.EQ.0 ) GOTO 1210
                j = Jcm(l)
                DO l = 1,nn
                  IF ( A(l,jbx).NE.A(l,j) ) GOTO 1205
                ENDDO
                GOTO 1250
 1205         CONTINUE
            ENDIF
 1210       DO i = 1,nn
              IF ( i.NE.lc ) THEN
                jex = Jx(i)
                IF ( DABS(A(lc,jbx)*A(i,jex)-A(lc,jex)*A(i,jbx))
     &               .LE.smalno ) GOTO 1250
              ENDIF
            ENDDO
            njc = njc + 1
            IF ( jbx.NE.Jcm(lc) ) newcom = .TRUE.
            Jcm(lc) = jbx
            lcs(njc) = lc
            GOTO 1300
          ENDIF
 1250   CONTINUE
 1300   En(jbx,Npt) = -En(jbx,Npt)
        IF ( njc.LT.nn ) GOTO 1200
      ENDIF
      DO j = 1,Ng
        En(j,Npt) = DABS(En(j,Npt))
      ENDDO
      IF ( newcom ) THEN
C SWITCH COMPONENTS
        DO lc = 1,nn
          jb = Jcm(lc)
          IF ( A(lc,jb).EQ.0. ) THEN
            jb = Jx(lc)
            Jcm(lc) = jb
          ENDIF
          tem = A(lc,jb)
          IF ( tem.NE.0. ) THEN
            pisave(lc) = H0(jb) - S(jb)
            IF ( jb.LE.Ng ) pisave(lc) = pisave(lc) + Enln(jb) + Tm
            cmp(lc) = Prod(jb)
C CALCULATE NEW COEFFICIENTS
            IF ( tem.NE.1. ) THEN
              B0(lc) = B0(lc)/tem
              B0p(lc,1) = B0p(lc,1)/tem
              B0p(lc,2) = B0p(lc,2)/tem
              DO j = 1,Nspx
                A(lc,j) = A(lc,j)/tem
              ENDDO
            ENDIF
            DO i = 1,nn
              IF ( A(i,jb).NE.0..AND.i.NE.lc ) THEN
                tem = A(i,jb)
                DO j = 1,Nspx
                  A(i,j) = A(i,j) - A(lc,j)*tem
                  IF ( DABS(A(i,j)).LT.1.E-5 ) A(i,j) = 0.
                ENDDO
                B0(i) = B0(i) - B0(lc)*tem
                B0p(i,1) = B0p(i,1) - B0p(lc,1)*tem
                B0p(i,2) = B0p(i,2) - B0p(lc,2)*tem
              ENDIF
            ENDDO
          ENDIF
        ENDDO
        IF ( Debug(Npt) ) THEN
          WRITE (IOOUT,99029)
          WRITE (IOOUT,99030) (cmp(k),k=1,nn)
        ENDIF
      ENDIF
      IF ( Msing.NE.0 ) THEN
C SWITCH ORDER OF MSING AND NLM COMPONENTS
        reduce = .TRUE.
        lelim = Nlm
        lsing = Nlm
        IF ( Msing.NE.Nlm ) THEN
          DO j = 1,Nspx
            aa = A(Msing,j)
            A(Msing,j) = A(Nlm,j)
            A(Nlm,j) = aa
          ENDDO
          ja = Jcm(Msing)
          Jcm(Msing) = Jcm(Nlm)
          Jcm(Nlm) = ja
          ae = cmp(Msing)
          cmp(Msing) = cmp(Nlm)
          cmp(Nlm) = ae
          ae = Elmt(Msing)
          Elmt(Msing) = Elmt(Nlm)
          Elmt(Nlm) = ae
          ja = Jx(Msing)
          Jx(Msing) = Jx(Nlm)
          Jx(Nlm) = ja
          aa = Atwt(Msing)
          Atwt(Msing) = Atwt(Nlm)
          Atwt(Nlm) = aa
          aa = B0(Msing)
          B0(Msing) = B0(Nlm)
          B0(Nlm) = aa
          aa = pisave(Msing)
          pisave(Msing) = pisave(Nlm)
          pisave(Nlm) = aa
          DO i = 1,2
            aa = B0p(Msing,i)
            B0p(Msing,i) = B0p(Nlm,i)
            B0p(Nlm,i) = aa
          ENDDO
        ENDIF
      ELSEIF ( .NOT.newcom.AND.Trace.EQ.0. ) THEN
        GOTO 600
      ENDIF
      Msing = 0
      tsize = xsize
      GOTO 500
 1400 Ttt(Npt) = Tt
      Ppp(Npt) = Pp
      Vlm(Npt) = Rr*Enn*Tt/Pp
      Hsum(Npt) = Hsum(Npt)*Tt
      Wm(Npt) = 1./Enn
      gasfrc = Enn/Totn(Npt)
      IF ( gasfrc.LT..0001 ) WRITE (IOOUT,99031) Npt,gasfrc
      IF ( Trace.NE.0. ) THEN
        DO 1450 j = 1,Ng
          IF ( lelim.NE.0 ) THEN
            DO i = lelim,ls
              IF ( A(i,j).NE.0. ) GOTO 1450
            ENDDO
          ENDIF
          IF ( Enln(j).GT.-87. ) En(j,Npt) = DEXP(Enln(j))
 1450   CONTINUE
      ENDIF
      IF ( Debug(Npt) ) WRITE (IOOUT,99032) Npt,Pp,Tt,Hsum(Npt),
     &                      Ssum(Npt),Wm(Npt),Cpr(Npt),Dlvpt(Npt),
     &                            Dlvtp(Npt),Gammas(Npt),Vlm(Npt)
      IF ( Tt.GE.Tg(1).AND.Tt.LE.Tg(4) ) GOTO 1600
      IF ( Shock ) GOTO 1600
      WRITE (IOOUT,99033) Tt,Npt
      IF ( Tt.GE.Tg(1)*.8D0.AND.Tt.LE.Tg(4)*1.1D0 ) GOTO 1600
      Npt = Npt + 1
 1500 Tt = 0.
      Npt = Npt - 1
      WRITE (IOOUT,99035) Npt
 1600 Lsave = Nlm
      Nlm = ls
      IF ( Npr.GT.0 ) Gonly = .FALSE.
      RETURN
99001 FORMAT (/' POINT ITN',6X,'T',10X,4(A4,8X)/(18X,5(A4,8X)))
99002 FORMAT (/' T DERIV MATRIX')
99003 FORMAT (/' P DERIV MATRIX')
99004 FORMAT (/' ITERATION',I3,6X,'MATRIX ')
99005 FORMAT (/' SOLUTION VECTOR',/,6x,5A15/8X,5A15)
99006 FORMAT (3X,5E15.6)
99007 FORMAT (/' DERIVATIVE MATRIX SINGULAR (EQLBRM)')
99008 FORMAT (/' LOW TEMPERATURE IMPLIES A CONDENSED SPECIES SHOULD HA',
     &        'VE BEEN INSERTED,',
     &        /' RESTART WITH insert DATASET (EQLBRM)')
99009 FORMAT (/' SINGULAR MATRIX, ITERATION',I3,'  VARIABLE',I3,
     &        '(EQLBRM)')
99010 FORMAT (/' WARNING!! POINT',I3,
     &        ' USES A REDUCED SET OF COMPONENTS',/
     &       ' SPECIES CONTAINING THE ELIMINATED COMPONENT ARE OMITTED.'
     &       ,/
     &   ' IT MAY BE NECESSARY TO RERUN WITH INSERTED CONDENSED SPECIES'
     &   ,/' CONTAINING COMPONENT ',A8,'(EQLBRM)')
99011 FORMAT (/' T=',E15.8,' ENN=',E15.8,' ENNL=',E15.8,' PP=',E15.8,
     &        /' LN P/N=',E15.8,' AMBDA=',E15.8)
99012 FORMAT (/' AMBDA SET BY ',A16)
99013 FORMAT (' VOLUME=',E15.8,'CC/G')
99014 FORMAT (/24X,'Nj',12X,'LN Nj',8X,'DEL LN Nj',6X,'H0j/RT',/,41X,
     &        'S0j/R',10X,' G0j/RT',8X,' Gj/RT')
99015 FORMAT (1X,A16,4E15.6,/35x,3E15.6)
99016 FORMAT (/' ELECTRON BALANCE ITER NO. =',I4,'  DELTA PI =',E14.7)
99017 FORMAT (/' DID NOT CONVERGE ON ELECTRON BALANCE (EQLBRM)')
99018 FORMAT (/' DELTA S/R =',E15.8)
99019 FORMAT (/,I4,' ITERATIONS DID NOT SATISFY CONVERGENCE',/,15x,
     &        ' REQUIREMENTS FOR THE POINT',I5,' (EQLBRM)')
99020 FORMAT (/' TRY REMOVING CONDENSED SPECIES (EQLBRM)')
99021 FORMAT (/' POINT ITN',6X,'T',10X,4A12/(18X,5A12))
99022 FORMAT (I4,I5,5F12.3,/(12X,5F12.3))
99023 FORMAT (' PHASE CHANGE, REPLACE ',A16,'WITH ',A16)
99024 FORMAT (/1X,A15,2F10.3,3X,E15.7)
99025 FORMAT (' [G0j-SUM(Aij*PIi)]/Mj =',E15.7,9X,'MAX NEG DELTA G =',
     &        E15.7)
99026 FORMAT (/' REINSERTION OF ',A16,' LIKELY TO CAUSE SINGULARITY,',
     &        '(EQLBRM)')
99027 FORMAT (' ADD ',A16)
99028 FORMAT (' REMOVE ',A16)
99029 FORMAT (/' NEW COMPONENTS')
99030 FORMAT (/2x,6A12)
99031 FORMAT (/' WARNING!  RESULTS MAY BE WRONG FOR POINT',I3,' DUE TO',
     &        /' LOW MOLE FRACTION OF GASES (',E15.8,') (EQLBRM)')
99032 FORMAT (/' POINT=',I3,3X,'P=',E13.6,3X,'T=',E13.6,/3X,'H/R=',
     &        E13.6,3X,'S/R=',E13.6,/3X,'M=',E13.6,3X,'CP/R=',E13.6,3X,
     &        'DLVPT=',E13.6,/3X,'DLVTP=',E13.6,3X,'GAMMA(S)=',E13.6,3X,
     &        'V=',E13.6)
99033 FORMAT (' THE TEMPERATURE=',E12.4,' IS OUT OF RANGE FOR POINT',I5,
     &        '(EQLBRM)')
99034 FORMAT (/,I3,' CONVERGENCES FAILED TO ESTABLISH SET OF CONDENSED',
     &        ' SPECIES (EQLBRM)')
99035 FORMAT (/' CALCULATIONS STOPPED AFTER POINT',I3,'(EQLBRM)')
      END
