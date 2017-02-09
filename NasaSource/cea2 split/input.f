      SUBROUTINE INPUT(Readok,Caseok,Ensert)
C***********************************************************************
C DECIPHER KEYWORDS, LITERAL VARIABLES, & NUMERICAL VARIABLES IN INPUT.
C***********************************************************************
      IMPLICIT NONE
      INCLUDE 'cea.inc'
C DUMMY ARGUMENTS
      LOGICAL Caseok,Readok
      CHARACTER*15 Ensert(20)
C LOCAL VARIABLES
      CHARACTER*15 cin(MAXNGC),cx15
      CHARACTER*4 code,cx4
      CHARACTER*1 cx1
      CHARACTER*2 cx2
      CHARACTER*3 cx3
      CHARACTER*26 lc,uc
      LOGICAL eqrats,incd,phi,pltdat,reacts,refl
      INTEGER i,ifrmla,ii,in,iv,ix,j,jj,k,lcin(MAXNGC),ncin,nmix
      INTEGER INDEX
      REAL*8 denmtr,dpin(MAXNGC),eratio,hr,mix(MAXNGC),ur,xyz
      REAL*8 DABS,DMIN1,DSQRT
C     gives the save atribute to the variables, so you can access them outsite this subroutine.
      SAVE cin,code,cx1,cx15,cx2,cx3,cx4,denmtr,dpin,eqrats,eratio,hr,i,
     &  ifrmla,ii,in,incd,iv,ix,j,jj,k,lcin,mix,ncin,nmix,phi,pltdat,
     &  reacts,refl,ur,xyz
C
      DATA uc/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/ !initialize uppercase
      DATA lc/'abcdefghijklmnopqrstuvwxyz'/ !initialize lowercase
      WRITE (IOOUT,99001) !ask for input file
      Caseok = .TRUE.
      Nonly = 0
      Nomit = 0
      Nsert = 0
      reacts = .FALSE.
      Trace = 0
      Short = .FALSE.
      Massf = .FALSE.
      DO i = 1,NCOL
        Debug(i) = .FALSE.
      ENDDO
      Nplt = 0
      Siunit = .TRUE.
      pltdat = .FALSE.
C CALL INFREE TO READ DATASET
 100  CALL INFREE(Readok,cin,ncin,lcin,dpin)
      IF ( .NOT.Readok ) GOTO 400
      code = cin(1)
      IF ( code.NE.'    ' ) THEN
C STORE PRODUCT NAMES FROM 'ONLY' DATASET
        IF ( code.EQ.'only' ) THEN
          Nonly = MIN(MAXNGC,ncin-1)
          DO i = 1,Nonly
            Prod(i) = cin(i+1)
          ENDDO
C STORE CONDENSED PRODUCT NAMES FROM 'INSERT' DATASET
        ELSEIF ( code.EQ.'inse' ) THEN
          Nsert = MIN(20,ncin-1)
          DO i = 1,Nsert
            Ensert(i) = cin(i+1)
          ENDDO
C STORE PRODUCT NAMES FROM 'OMIT' DATASET
        ELSEIF ( code.EQ.'omit' ) THEN
C CHECK OMIT DATASET
          Nomit = MIN(MAXNGC,ncin-1)
          DO i = 1,Nomit
            Omit(i) = cin(i+1)
          ENDDO
C KEYWORD 'THER' READ
C CALL UTHERM TO CONVERT FORMATTED THERMODYNAMIC DATA
        ELSEIF ( code.EQ.'ther' ) THEN
          Newr = .TRUE.
          REWIND IOTHM
          CALL UTHERM(Readok)
          IF ( .NOT.Readok ) THEN
            WRITE (IOOUT,99025)
            GOTO 400
          ENDIF
C KEYWORD 'TRAN' READ
C CALL UTRAN TO CONVERT FORMATTED TRANSPORT PROPERTIES
        ELSEIF ( code.EQ.'tran' ) THEN
          CALL UTRAN(Readok)
          IF ( .NOT.Readok ) THEN
            WRITE (IOOUT,99025)
            GOTO 400
          ENDIF
C PROCESS 'OUTP' DATASET.
        ELSEIF ( code.EQ.'outp' ) THEN
          DO 120 i = 2,ncin
            IF ( lcin(i).LT.0 ) THEN
              cx2 = cin(i)(1:2)
              cx3 = cin(i)(1:3)
              cx4 = cin(i)(1:4)
              IF ( cx3.EQ.'cal' ) THEN
                Siunit = .FALSE.
              ELSEIF ( cx4.EQ.'tran'.OR.cx3.EQ.'trn' ) THEN
                Trnspt = .TRUE.
              ELSEIF ( cx4.EQ.'trac' ) THEN
                Trace = dpin(i+1)
              ELSEIF ( cin(i)(:5).EQ.'short' ) THEN
                Short = .TRUE.
              ELSEIF ( cin(i)(:5).EQ.'massf' ) THEN
                Massf = .TRUE.
              ELSEIF ( cx3.EQ.'deb'.OR.cx3.EQ.'dbg' ) THEN
                DO j = i + 1,ncin
                  IF ( lcin(j).NE.i ) GOTO 120
                  k = dpin(j)
                  IF ( k.LE.NCOL ) Debug(k) = .TRUE.
                  lcin(j) = 0
                ENDDO
              ELSEIF ( cx2.EQ.'si' ) THEN
                Siunit = .TRUE.
              ELSEIF ( pltdat.AND.Nplt.LT.20 ) THEN
                Nplt = Nplt + 1
                Pltvar(Nplt) = cin(i)
              ELSEIF ( cx2.EQ.'pl' ) THEN
                pltdat = .TRUE.
              ELSE
                WRITE (IOOUT,99002) cin(i)
              ENDIF
            ENDIF
 120      CONTINUE
C SORT AND STORE DATA FROM 'REAC' DATASET.
        ELSEIF ( code.EQ.'reac' ) THEN
          reacts = .TRUE.
          Moles = .FALSE.
          Nreac = 0
          DO i = 1,MAXR
            Pecwt(i) = -1.
          ENDDO
          i = 1
 140      i = i + 1
          IF ( i.LE.ncin ) THEN
            IF ( lcin(i).NE.0 ) THEN
              IF ( lcin(i).GT.0 ) THEN
                WRITE (IOOUT,99003) cin(i)
                GOTO 140
              ENDIF
              cx15 = cin(i)
              cx1 = cx15(:1)
              cx2 = cx15(:2)
              cx3 = cx15(:3)
              cx4 = cx15(:4)
C NEW REACTANT
              IF ( cx2.NE.'na'.AND.cx2.NE.'ox'.AND.cx2.NE.'fu' ) THEN
C LOOK FOR PERCENTS
                IF ( cx1.EQ.'m'.OR.cx1.EQ.'w' ) THEN
                  IF ( lcin(i+1).GT.0 ) THEN
                    i = i + 1
                    Pecwt(Nreac) = dpin(i)
                  ELSE
                    Caseok = .FALSE.
                    WRITE (IOOUT,99004)
                  ENDIF
                  IF ( cx1.EQ.'m'.AND.Nreac.EQ.1 ) Moles = .TRUE.
                  IF ( cx1.EQ.'m'.AND..NOT.Moles.OR.cx1.EQ.'w'.AND.
     &                 Moles ) THEN
                    Caseok = .FALSE.
                    WRITE (IOOUT,99005)
                  ENDIF
                  GOTO 140
                ENDIF
C LOOK FOR TEMPERATURES
                IF ( cx1.EQ.'t' ) THEN
                  IF ( lcin(i+1).GT.0 ) THEN
                    i = i + 1
                    Rtemp(Nreac) = dpin(i)
                    IF ( lcin(i-1).LT.1 ) THEN
                      IF ( INDEX(cx15,'r').GT.0 ) Rtemp(Nreac)
     &                     = Rtemp(Nreac)/1.8D0
                      IF ( INDEX(cx15,'c').GT.0 ) Rtemp(Nreac)
     &                     = Rtemp(Nreac) + 273.15D0
                      IF ( INDEX(cx15,'f').GT.0 ) Rtemp(Nreac)
     &                     = (Rtemp(Nreac)-32.D0)/1.8D0 + 273.15D0
                    ENDIF
                  ELSE
                    WRITE (IOOUT,99006)
                    Caseok = .FALSE.
                  ENDIF
                  GOTO 140
                ENDIF
C LOOK FOR ENTHALPY
                IF ( cx1.EQ.'h'.OR.cx1.EQ.'u' ) THEN
                  Energy(Nreac) = cx15
                  IF ( lcin(i+1).GT.0 ) THEN
                    i = i + 1
                    Enth(Nreac) = dpin(i)*1000.D0/Rr
                    IF ( INDEX(cin(i-1),'c').GT.0 ) Enth(Nreac)
     &                   = Enth(Nreac)*4.184D0
                    IF ( INDEX(cin(i-1),'k').GT.0 ) Enth(Nreac)
     &                   = Enth(Nreac)*1000.D0
                  ENDIF
                  GOTO 140
                ENDIF
C LOOK FOR DENSITY
                IF ( cx3.EQ.'rho'.OR.cx3.EQ.'den' ) THEN
                  IF ( lcin(i+1).GT.0 ) THEN
                    i = i + 1
                    Dens(Nreac) = dpin(i)
                    IF ( INDEX(cx15,'kg').GT.0 ) Dens(Nreac)
     &                   = Dens(Nreac)/1000.D0
                  ENDIF
                  GOTO 140
                ENDIF
C CHECK FOR CHEMICAL SYMBOLS IN EXPLODED FORMULA
                IF ( (lcin(i).EQ.-1.OR.lcin(i).EQ.-2).AND.INDEX(uc,cx1)
     &               .GT.0 ) THEN
                  Energy(Nreac) = ' '
                  ifrmla = ifrmla + 1
                  Nfla(Nreac) = ifrmla
                  IF ( lcin(i).EQ.-2 ) THEN
                    ix = INDEX(lc,cx2(2:2))
                    IF ( ix.GT.0 ) cx2(2:2) = uc(ix:ix)
                  ENDIF
                  Ratom(Nreac,ifrmla) = cx2
                  IF ( lcin(i+1).EQ.i ) THEN
                    Rnum(Nreac,ifrmla) = dpin(i+1)
                  ELSE
                    Rnum(Nreac,ifrmla) = 1.
                  ENDIF
                  i = i + 1
                  GOTO 140
                ENDIF
                WRITE (IOOUT,99007) cin(i)
              ELSE
                Nreac = MIN(Nreac+1,MAXR)
                Fox(Nreac) = cx15
                i = i + 1
                IF ( lcin(i).LT.0 ) Rname(Nreac) = cin(i)
                ifrmla = 0
                Nfla(Nreac) = 0
                Energy(Nreac) = 'lib'
                Enth(Nreac) = 0.
                Jray(Nreac) = 0
                Pecwt(Nreac) = -1.
                Rnum(Nreac,1) = 0.
                Rmw(Nreac) = 0.
                Rtemp(Nreac) = 0.
              ENDIF
            ENDIF
            GOTO 140
          ENDIF
C SORT AND STORE INPUT FROM 'PROB' DATASET
        ELSEIF ( code.EQ.'prob' ) THEN
          Case = ' '
          DO i = 1,MAXPV
            P(i) = 0.
            V(i) = 0.
          ENDDO
          DO i = 1,MAXT
            T(i) = 0.
          ENDDO
          P(1) = 1.
          Trace = 0.
          Lsave = 0
          R = Rr/4184.D0
          S0 = 0.
          hr = 1.D30
          ur = 1.D30
          Tp = .FALSE.
          Hp = .FALSE.
          Sp = .FALSE.
          Rkt = .FALSE.
          Shock = .FALSE.
          Detn = .FALSE.
          Vol = .FALSE.
          Ions = .FALSE.
          Eql = .FALSE.
          Froz = .FALSE.
          Fac = .FALSE.
          Debugf = .FALSE.
          Acat = 0.
          Ma = 0.
          Nfz = 1
          Nsub = 0
          Nsup = 0
          Npp = 0
          Tcest = 3800.
          DO i = 1,NCOL
            Pcp(i) = 0.
            Pcp(i+NCOL) = 0.
            Supar(i) = 0.
            Subar(i) = 0.
            Mach1(i) = 0.
            U1(i) = 0.
          ENDDO
          Gamma1 = 0.
          phi = .FALSE.
          eqrats = .FALSE.
          incd = .FALSE.
          refl = .FALSE.
          Shkdbg = .FALSE.
          Incdeq = .FALSE.
          Incdfz = .FALSE.
          Refleq = .FALSE.
          Reflfz = .FALSE.
          Np = 0
          Nt = 1
          Trnspt = .FALSE.
C PROCESS LITERAL VARIABLES IN 'PROB' DATASET THAT DO NOT HAVE
C ASSOCIATED NUMERICAL DATA.
          DO 160 i = 2,ncin
            IF ( lcin(i).LT.0 ) THEN
              DO j = i + 1,ncin
                IF ( lcin(j).EQ.i ) GOTO 160
              ENDDO
              cx15 = cin(i)
              cx2 = cx15(:2)
              cx3 = cx15(:3)
              cx4 = cx15(:4)
              IF ( cx4.EQ.'case' ) THEN
                Case = cin(i+1)
                lcin(i+1) = 0
              ELSEIF ( cx2.EQ.'tp'.OR.cx2.EQ.'pt' ) THEN
                Tp = .TRUE.
              ELSEIF ( cx2.EQ.'hp'.OR.cx2.EQ.'ph' ) THEN
                Hp = .TRUE.
              ELSEIF ( cx2.EQ.'sp'.OR.cx2.EQ.'ps' ) THEN
                Sp = .TRUE.
              ELSEIF ( cx2.EQ.'sv'.OR.cx2.EQ.'vs' ) THEN
                Sp = .TRUE.
                Vol = .TRUE.
              ELSEIF ( cx2.EQ.'uv'.OR.cx2.EQ.'vu' ) THEN
                Hp = .TRUE.
                Vol = .TRUE.
              ELSEIF ( cx2.EQ.'tv'.OR.cx2.EQ.'vt' ) THEN
                Tp = .TRUE.
                Vol = .TRUE.
              ELSEIF ( cx2.EQ.'ro'.OR.cx3.EQ.'rkt' ) THEN
                Rkt = .TRUE.
              ELSEIF ( cx3.EQ.'dbg'.OR.cx3.EQ.'deb' ) THEN
                Debugf = .TRUE.
                Shkdbg = .TRUE.
                Detdbg = .TRUE.
              ELSEIF ( cx3.EQ.'fac' ) THEN
                Rkt = .TRUE.
                Eql = .TRUE.
                Fac = .TRUE.
                Froz = .FALSE.
              ELSEIF ( cx2.EQ.'eq' ) THEN
                Eql = .TRUE.
              ELSEIF ( cx2.EQ.'fr'.OR.cx2.EQ.'fz' ) THEN
                Froz = .TRUE.
              ELSEIF ( cx2.EQ.'sh' ) THEN
                Shock = .TRUE.
              ELSEIF ( cx3.EQ.'inc' ) THEN
                Shock = .TRUE.
                incd = .TRUE.
                IF ( INDEX(cx15,'eq').GT.0 ) Eql = .TRUE.
                IF ( INDEX(cx15,'fr').GT.0 ) Froz = .TRUE.
                IF ( INDEX(cx15,'fz').GT.0 ) Froz = .TRUE.
              ELSEIF ( cx3.EQ.'ref' ) THEN
                Shock = .TRUE.
                refl = .TRUE.
                IF ( INDEX(cx15,'eq').GT.0 ) Eql = .TRUE.
                IF ( INDEX(cx15,'fr').GT.0 ) Froz = .TRUE.
                IF ( INDEX(cx15,'fz').GT.0 ) Froz = .TRUE.
              ELSEIF ( cx3.EQ.'det' ) THEN
                Detn = .TRUE.
              ELSEIF ( cx4.EQ.'ions' ) THEN
                Ions = .TRUE.
              ELSE
                WRITE (IOOUT,99002) cx15
              ENDIF
              lcin(i) = 0
            ENDIF
 160      CONTINUE
          iv = 2
          Nof = 0
          GOTO 200
        ELSEIF ( code(1:3).EQ.'end' ) THEN
          IF ( Shock ) THEN
            IF ( incd.AND.Froz ) Incdfz = .TRUE.
            IF ( incd.AND.Eql ) Incdeq = .TRUE.
            IF ( refl.AND.Froz ) Reflfz = .TRUE.
            IF ( refl.AND.Eql ) Refleq = .TRUE.
          ENDIF
          Hsub0 = DMIN1(hr,ur)
          Size = 0.
          IF ( hr.GT..9D30 ) hr = 0.D0
          IF ( ur.GT..9D30 ) ur = 0.D0
          IF ( Trnspt ) Viscns = .3125*DSQRT(1.E5*Boltz/(Pi*Avgdr))
          IF ( Siunit ) R = Rr/1000.
          IF ( Detn.OR.Shock ) Newr = .TRUE.
          IF ( .NOT.Short ) THEN
            WRITE (IOOUT,99008) Tp,(Hp.AND..NOT.Vol),Sp,(Tp.AND.Vol),
     &                      (Hp.AND.Vol),(Sp.AND.Vol),Detn,Shock,refl,
     &                      incd,Rkt,Froz,Eql,Ions,Siunit,Debugf,Shkdbg,
     &                      Detdbg,Trnspt
            IF ( T(1).GT.0. ) WRITE (IOOUT,99009) (T(jj),jj=1,Nt)
            WRITE (IOOUT,99010) Trace,S0,hr,ur
            IF ( Np.GT.0.AND.Vol ) WRITE (IOOUT,99011)
     &           (V(jj)*1.D-05,jj=1,Np)
          ENDIF
          IF ( Rkt ) THEN
            IF ( Nt.EQ.0 ) Hp = .TRUE.
            IF ( .NOT.Short ) THEN
              WRITE (IOOUT,99012) (P(jj),jj=1,Np)
              WRITE (IOOUT,99013) (Pcp(jj),jj=1,Npp)
              WRITE (IOOUT,99014) (Subar(i),i=1,Nsub)
              WRITE (IOOUT,99015) (Supar(i),i=1,Nsup)
              WRITE (IOOUT,99016) Nfz,Ma,Acat
            ENDIF
          ELSE
            IF ( .NOT.Vol.AND..NOT.Short ) WRITE (IOOUT,99017)
     &           (P(jj),jj=1,Np)
          ENDIF
          IF ( reacts ) CALL REACT
          IF ( Nreac.EQ.0.OR.Nlm.LE.0 ) THEN
            WRITE (IOOUT,99018)
            Caseok = .FALSE.
            WRITE (IOOUT,99025)
            GOTO 400
          ENDIF
          IF ( Nof.EQ.0 ) THEN
            Nof = 1
            Oxf(1) = 0.
            IF ( Wp(2).GT.0. ) THEN
              Oxf(1) = Wp(1)/Wp(2)
            ELSE
              Caseok = .FALSE.
              WRITE (IOOUT,99004)
              WRITE (IOOUT,99025)
              GOTO 400
            ENDIF
          ELSEIF ( phi.OR.eqrats ) THEN
            DO i = 1,Nof
              eratio = Oxf(i)
              IF ( eqrats ) THEN
                xyz = -eratio*Vmin(2) - Vpls(2)
                denmtr = eratio*Vmin(1) + Vpls(1)
              ELSE
                xyz = -Vmin(2) - Vpls(2)
                denmtr = eratio*(Vmin(1)+Vpls(1))
              ENDIF
              IF ( DABS(denmtr).LT.1.D-30 ) THEN
                Caseok = .FALSE.
                WRITE (IOOUT,99019) eratio
                WRITE (IOOUT,99025)
                GOTO 400
              ENDIF
              Oxf(i) = xyz/denmtr
            ENDDO
          ENDIF
          IF ( .NOT.Sp.AND..NOT.Tp.AND..NOT.Hp.AND..NOT.Rkt.AND.
     &         .NOT.Detn.AND..NOT.Shock ) THEN
            Caseok = .FALSE.
            WRITE (IOOUT,99020)
          ELSEIF ( Tp.AND.T(1).LE.0. ) THEN
            Caseok = .FALSE.
            WRITE (IOOUT,99021)
          ELSEIF ( Np.LE.0 ) THEN
            Caseok = .FALSE.
            WRITE (IOOUT,99022)
          ENDIF
          IF ( .NOT.(Caseok.AND.Nlm.GT.0) ) WRITE (IOOUT,99025)
          GOTO 400
        ELSE
          WRITE (IOOUT,99023)
        ENDIF
      ENDIF
      GOTO 100
C PROCESS NUMERICAL DATA FOLLOWING 'PROB' LITERALS
 200  in = 0
      nmix = 0
      ii = iv
      DO i = ii,ncin
        iv = i
        IF ( lcin(i).NE.0 ) THEN
          IF ( lcin(i).LT.0 ) THEN
            IF ( in.GT.0 ) GOTO 300
            in = i
          ELSE
            IF ( lcin(i).NE.in ) GOTO 300
            nmix = nmix + 1
            mix(nmix) = dpin(i)
            lcin(i) = 0
          ENDIF
        ENDIF
      ENDDO
 300  IF ( nmix.LE.0 ) THEN
        IF ( iv.LT.ncin ) GOTO 200
        GOTO 100
      ENDIF
      cx15 = cin(in)
      cx1 = cx15(:1)
      cx2 = cx15(:2)
      cx3 = cx15(:3)
      cx4 = cx15(:4)
      IF ( cx1.EQ.'t' ) THEN
        Nt = nmix
        IF ( nmix.GT.MAXMIX ) THEN
          Nt = MAXMIX
          WRITE (IOOUT,99024) 't',Nt
        ENDIF
        DO i = 1,Nt
          IF ( cx4.NE.'tces' ) THEN
            T(i) = mix(i)
            IF ( lcin(in).LT.-1 ) THEN
              IF ( INDEX(cx15,'r').GT.0 ) T(i) = T(i)/1.8D0
              IF ( INDEX(cx15,'c').GT.0 ) T(i) = T(i) + 273.15D0
              IF ( INDEX(cx15,'f').GT.0 ) T(i) = (T(i)-32.D0)
     &             /1.8D0 + 273.15D0
            ENDIF
          ENDIF
        ENDDO
      ELSEIF ( (cx2.EQ.'pc'.OR.cx2.EQ.'pi').AND.INDEX(cx15(3:15),'p')
     &         .GT.0.AND.INDEX(cx15,'psi').EQ.0 ) THEN
        Npp = nmix
        IF ( nmix.GT.2*NCOL ) THEN
          Npp = 2*NCOL
          WRITE (IOOUT,99024) 'pcp',Npp
        ENDIF
        DO i = 1,Npp
          Pcp(i) = mix(i)
        ENDDO
      ELSEIF ( cx1.EQ.'p'.AND.cx3.NE.'phi' ) THEN
        Np = nmix
        IF ( nmix.GT.MAXPV ) THEN
          Np = MAXPV
          WRITE (IOOUT,99024) 'p',Np
        ENDIF
        DO 350 i = 1,Np
          P(i) = mix(i)
          IF ( INDEX(cx15,'psi').NE.0 ) THEN
            P(i) = P(i)/14.696006D0
          ELSEIF ( INDEX(cx15,'mmh').NE.0 ) THEN
            P(i) = P(i)/760.D0
          ELSEIF ( INDEX(cx15,'atm').EQ.0 ) THEN
            GOTO 350
          ENDIF
          P(i) = P(i)*1.01325D0
 350    CONTINUE
      ELSEIF ( cx3.EQ.'rho' ) THEN
        xyz = 1.D02
        IF ( INDEX(cx15,'kg').NE.0 ) xyz = 1.D05
        Np = nmix
        IF ( nmix.GT.MAXPV ) THEN
          Np = MAXPV
          WRITE (IOOUT,99024) 'rho',Np
        ENDIF
        DO i = 1,Np
          V(i) = xyz/mix(i)
        ENDDO
      ELSEIF ( cx1.EQ.'v' ) THEN
        xyz = 1.D02
        IF ( INDEX(cx15,'kg').NE.0 ) xyz = 1.D05
        Np = nmix
        IF ( nmix.GT.MAXPV ) THEN
          Np = MAXPV
          WRITE (IOOUT,99024) 'v',Np
        ENDIF
        DO i = 1,Np
          V(i) = mix(i)*xyz
        ENDDO
      ELSEIF ( cx3.EQ.'nfz'.OR.cx3.EQ.'nfr' ) THEN
        Nfz = mix(1)
        Froz = .TRUE.
      ELSEIF ( cx4.EQ.'tces' ) THEN
        Tcest = mix(1)
      ELSEIF ( cx4.EQ.'trac' ) THEN
        Trace = mix(1)
      ELSEIF ( cx3.EQ.'s/r' ) THEN
        S0 = mix(1)
      ELSEIF ( cx3.EQ.'u/r'.OR.cx2.EQ.'ur' ) THEN
        ur = mix(1)
      ELSEIF ( cx3.EQ.'h/r'.OR.cx2.EQ.'hr' ) THEN
        hr = mix(1)
      ELSEIF ( cx2.EQ.'u1' ) THEN
        Nsk = nmix
        IF ( nmix.GT.NCOL ) THEN
          Nsk = NCOL
          WRITE (IOOUT,99024) 'u1',Nsk
        ENDIF
        DO i = 1,Nsk
          U1(i) = mix(i)
        ENDDO
      ELSEIF ( cx4.EQ.'mach' ) THEN
        Nsk = nmix
        IF ( nmix.GT.NCOL ) THEN
          Nsk = NCOL
          WRITE (IOOUT,99024) 'mach1',Nsk
        ENDIF
        DO i = 1,Nsk
          Mach1(i) = mix(i)
        ENDDO
      ELSEIF ( cx3.EQ.'sub' ) THEN
        Nsub = nmix
        IF ( nmix.GT.13 ) THEN
          Nsub = 13
          WRITE (IOOUT,99024) 'subar',Nsub
        ENDIF
        DO i = 1,Nsub
          Subar(i) = mix(i)
        ENDDO
      ELSEIF ( cx3.EQ.'sup' ) THEN
        Nsup = nmix
        IF ( nmix.GT.13 ) THEN
          Nsup = 13
          WRITE (IOOUT,99024) 'supar',Nsup
        ENDIF
        DO i = 1,Nsup
          Supar(i) = mix(i)
        ENDDO
      ELSEIF ( cx2.EQ.'ac' ) THEN
        Acat = mix(1)
      ELSEIF ( cx4.EQ.'mdot'.OR.cx2.EQ.'ma' ) THEN
        Ma = mix(1)
      ELSEIF ( cx4.EQ.'case' ) THEN
        Case = cin(in+1)
        lcin(in+1) = 0
      ELSEIF ( Nof.EQ.0.AND.
     &         (cx3.EQ.'phi'.OR.cx3.EQ.'o/f'.OR.cx3.EQ.'f/a'.OR.
     &         cx2.EQ.'%f'.OR.cx1.EQ.'r') ) THEN
        Nof = nmix
        IF ( nmix.GT.MAXMIX ) THEN
          Nof = MAXMIX
          WRITE (IOOUT,99024) 'o/f',Nof
        ENDIF
        DO k = 1,Nof
          Oxf(k) = mix(k)
        ENDDO
        IF ( cx3.EQ.'phi' ) THEN
          phi = .TRUE.
        ELSEIF ( cx1.EQ.'r' ) THEN
          eqrats = .TRUE.
        ELSEIF ( cx3.EQ.'f/a' ) THEN
          DO k = 1,Nof
            IF ( Oxf(k).GT.0. ) Oxf(k) = 1./Oxf(k)
          ENDDO
        ELSEIF ( cx4.EQ.'%fue' ) THEN
          DO k = 1,Nof
            IF ( Oxf(k).GT.0. ) Oxf(k) = (100.-Oxf(k))/Oxf(k)
          ENDDO
        ENDIF
      ELSE
        WRITE (IOOUT,99002) cx15
      ENDIF
      IF ( iv.GE.ncin ) GOTO 100
      GOTO 200
 400  RETURN
99001 FORMAT (/,/)
99002 FORMAT ('  WARNING!!  DID NOT RECOGNIZE ',A15,' (INPUT)'/)
99003 FORMAT (/' WARNING!!  LITERAL EXPECTED FOR ',A15,'(INPUT)')
99004 FORMAT (/' REACTANT AMOUNT MISSING (INPUT)')
99005 FORMAT (/' MOLES AND WEIGHT PERCENTS SHOULD NOT BE MIXED (INPUT)')
99006 FORMAT (/' REACTANT TEMPERATURE MISSING (INPUT) ')
99007 FORMAT (/' WARNING!! ',A15,' NOT RECOGNIZED (INPUT)')
99008 FORMAT (/' OPTIONS: TP=',L1,'  HP=',L1,'  SP=',L1,'  TV=',L1,
     &        '  UV=',L1,'  SV=',L1,'  DETN=',L1,'  SHOCK=',L1,
     &        '  REFL=',L1,'  INCD=',L1,/' RKT=',L1,'  FROZ=',L1,
     &        '  EQL=',L1,'  IONS=',L1,'  SIUNIT=',L1,'  DEBUGF=',L1,
     &        '  SHKDBG=',L1,'  DETDBG=',L1,'  TRNSPT=',L1)
99009 FORMAT (/' T,K =',7F11.4)
99010 FORMAT (/1p,' TRACE=',E9.2,'  S/R=',E13.6,'  H/R=',E13.6,'  U/R=',
     &        E13.6)
99011 FORMAT (/' SPECIFIC VOLUME,M**3/KG =',1p,(4E14.7))
99012 FORMAT (/' Pc,BAR =',7F13.6)
99013 FORMAT (/' Pc/P =',9F11.4)
99014 FORMAT (/' SUBSONIC AREA RATIOS =',(5F11.4))
99015 FORMAT (/' SUPERSONIC AREA RATIOS =',(5F11.4))
99016 FORMAT (/' NFZ=',i3,1p,'  Mdot/Ac=',e13.6,'  Ac/At=',e13.6)
99017 FORMAT (/' P,BAR =',7F13.6)
99018 FORMAT (/' ERROR IN REACTANTS DATASET (INPUT)')
99019 FORMAT (/' UNABLE TO PROCESS EQUIVALENCE RATIO =',E11.4,'(INPUT)')
99020 FORMAT (/' TYPE OF PROBLEM NOT SPECIFIED (INPUT)')
99021 FORMAT (/' ASSIGNED VALUES OF TEMPERATURE ARE MISSING IN prob',
     &        ' DATASET (INPUT)')
99022 FORMAT (/' ASSIGNED PRESSURE (OR DENSITY) MISSING IN prob',
     &        ' DATASET (INPUT)')
99023 FORMAT (/' WARNING!!  A KEYWORD IS MISSING (INPUT)')
99024 FORMAT (/' NOTE!! MAXIMUM NUMBER OF ASSIGNED ',A5,' VALUES IS',I3,
     &        ' (INPUT)',/)
99025 FORMAT (/' FATAL ERROR IN DATASET (INPUT)')
      END
