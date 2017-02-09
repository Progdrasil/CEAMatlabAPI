      SUBROUTINE OUT1
C***********************************************************************
C OUT1 WRITES REACTANT AND FUEL-OXIDANT RATIO INFORMATION.
C ENTRY OUT2 WRITES THERMODYNAMIC PROPERTIES.
C ENTRY OUT3 WRITES MOLE FRACTIONS.
C ENTRY OUT4 WRITES TRANSPORT PROPERTIES.
C
C NOTE - ROCKET, SHOCK, AND DETON PROBLEMS HAVE ADDITIONAL OUTPUT.
C***********************************************************************
      IMPLICIT NONE
      INCLUDE 'cea.inc'
C LOCAL VARIABLES
      CHARACTER*15 fc,fgi,fh,fp,frh,fs,fu
      CHARACTER*4 mamo
      INTEGER i,im,ione,j,k,kin,m,mcond,mcondf,mcp,mdvp,mdvt,meq,mfa,
     &        mg,mgam,mh,mie,mm,mmw,mof,mp,mpf,mph,mpn,mpnf,mrho,ms,
     &        mson,mt,mvis,mxx(24),n,notuse
      INTEGER INDEX
      LOGICAL kok
      REAL*8 pfactor,pfuel,phi,rho,tem,tra,vnum
      SAVE fc,fgi,fh,fp,frh,fs,fu,i,im,ione,j,k,kin,kok,m,
     &  mamo,mcond,mcondf,mcp,mdvp,mdvt,meq,mfa,mg,mgam,mh,mie,mm,mmw,
     &  mof,mp,mpf,mph,mpn,mpnf,mrho,ms,mson,mt,mvis,mxx,n,notuse,
     &  pfactor,pfuel,phi,rho,tem,tra,vnum
C
      EQUIVALENCE (mxx(1),mp)
      EQUIVALENCE (mxx(2),mt)
      EQUIVALENCE (mxx(3),mrho)
      EQUIVALENCE (mxx(4),mh)
      EQUIVALENCE (mxx(5),mie)
      EQUIVALENCE (mxx(6),mg)
      EQUIVALENCE (mxx(7),ms)
      EQUIVALENCE (mxx(8),mm)
      EQUIVALENCE (mxx(9),mcp)
      EQUIVALENCE (mxx(10),mgam)
      EQUIVALENCE (mxx(11),mson)
      EQUIVALENCE (mxx(12),mcond)
      EQUIVALENCE (mxx(13),mvis)
      EQUIVALENCE (mxx(14),mpn)
      EQUIVALENCE (mxx(15),mpf)
      EQUIVALENCE (mxx(16),mof)
      EQUIVALENCE (mxx(17),mph)
      EQUIVALENCE (mxx(18),meq)
      EQUIVALENCE (mxx(19),mfa)
      EQUIVALENCE (mxx(20),mmw)
      EQUIVALENCE (mxx(21),mdvt)
      EQUIVALENCE (mxx(22),mdvp)
      EQUIVALENCE (mxx(23),mcondf)
      EQUIVALENCE (mxx(24),mpnf)
      WRITE (IOOUT,99001) Case
      IF ( Moles ) THEN
        WRITE (IOOUT,99002) '   MOLES   '
        IF ( .NOT.Siunit ) WRITE (IOOUT,99003)
        IF ( Siunit ) WRITE (IOOUT,99004)
      ELSE
        WRITE (IOOUT,99002) 'WT FRACTION'
        IF ( .NOT.Siunit ) WRITE (IOOUT,99005)
        IF ( Siunit ) WRITE (IOOUT,99006)
      ENDIF
      DO n = 1,Nreac
        WRITE (IOOUT,99007) Fox(n),Rname(n),Pecwt(n),Enth(n)*R,Rtemp(n)
      ENDDO
      phi = 0.
      tem = (Vpls(1)+Vmin(1))*Oxfl
      IF ( ABS(tem).GE.1.D-3 ) phi = -(Vmin(2)+Vpls(2))/tem
      IF ( Fox(1).EQ.'NAME' ) THEN
        pfuel = 0.
      ELSE
        pfuel = 100.D0/(1.D0+Oxfl)
      ENDIF
      IF ( Rh(1).NE.0..OR.Rh(2).NE.0. ) THEN
        IF ( Rh(1).EQ.0..OR.Rh(2).EQ.0. ) THEN
          rho = MAX(Rh(1),Rh(2))
        ELSE
          rho = (Oxfl+1.)*Rh(1)*Rh(2)/(Rh(1)+Oxfl*Rh(2))
        ENDIF
        IF ( Siunit ) THEN
          rho = rho*1000.D0
          WRITE (IOOUT,99021) rho
        ELSE
          WRITE (IOOUT,99022) rho
        ENDIF
      ENDIF
      WRITE (IOOUT,99008) Oxfl,pfuel,Eqrat,phi
      RETURN
C***********************************************************************
      ENTRY OUT2
      ione = 0
      IF ( Rkt.AND..NOT.Page1 ) THEN
        ione = 2
        IF ( Iopt.NE.0 ) ione = 3
      ENDIF
C SET MXX ARRAY FOR PLOTTING PARAMETERS
      DO i = 1,24
        mxx(i) = 0
      ENDDO
      DO 100 i = 1,Nplt
        IF ( INDEX(Pltvar(i)(2:),'1').EQ.0 ) THEN
          IF ( INDEX(Pltvar(i)(1:),'dlnt').NE.0 ) THEN
            mdvt = i
          ELSEIF ( INDEX(Pltvar(i)(1:),'dlnp').NE.0 ) THEN
            mdvp = i
          ELSEIF ( Pltvar(i)(:4).EQ.'pran' ) THEN
            IF ( INDEX(Pltvar(i)(3:),'fz').NE.0 .OR.
     &           INDEX(Pltvar(i)(3:),'fr').NE.0 ) THEN
              mpnf = i
            ELSE
              mpn = i
            ENDIF
          ELSEIF ( Pltvar(i)(:4).EQ.'cond' ) THEN
            IF ( INDEX(Pltvar(i)(3:),'fz').NE.0 .OR.
     &           INDEX(Pltvar(i)(3:),'fr').NE.0 ) THEN
              mcondf = i
            ELSE
              mcond = i
            ENDIF
          ELSEIF ( Pltvar(i)(:3).EQ.'phi' ) THEN
            mph = i
          ELSEIF ( Pltvar(i)(:2).EQ.'p ' ) THEN
            mp = i
          ELSEIF ( Pltvar(i)(:1).EQ.'t' ) THEN
            mt = i
          ELSEIF ( Pltvar(i)(:3).EQ.'rho' ) THEN
            mrho = i
          ELSEIF ( Pltvar(i)(:1).EQ.'h' ) THEN
            mh = i
          ELSEIF ( Pltvar(i)(:1).EQ.'u' ) THEN
            mie = i
          ELSEIF ( Pltvar(i)(:3).EQ.'gam' ) THEN
            mgam = i
          ELSEIF ( Pltvar(i)(:3).EQ.'son' ) THEN
            mson = i
          ELSEIF ( Pltvar(i)(:2).EQ.'g ' ) THEN
            mg = i
          ELSEIF ( Pltvar(i)(:2).EQ.'s ' ) THEN
            ms = i
          ELSEIF ( Pltvar(i)(:1).EQ.'m'.AND.Pltvar(i)(:2).NE.'ma' ) THEN
            IF ( .NOT.Gonly.AND.Pltvar(i)(:2).EQ.'mw' ) THEN
              mmw = i
            ELSE
              mm = i
            ENDIF
          ELSEIF ( Pltvar(i)(:2).EQ.'cp' ) THEN
            mcp = i
          ELSEIF ( Pltvar(i)(:3).EQ.'vis' ) THEN
            mvis = i
          ELSEIF ( Pltvar(i)(:3).EQ.'o/f' ) THEN
            mof = i
          ELSEIF ( Pltvar(i)(:2).EQ.'%f' ) THEN
            mpf = i
          ELSEIF ( Pltvar(i)(:3).EQ.'f/a' ) THEN
            mfa = i
          ELSEIF ( Pltvar(i)(:1).EQ.'r' ) THEN
            meq = i
          ENDIF
        ENDIF
 100  CONTINUE
      DO i = Iplt + 1,Iplt + Npt
        IF ( mof.GT.0 ) Pltout(i,mof) = Oxfl
        IF ( mpf.GT.0 ) Pltout(i,mpf) = pfuel
        IF ( mph.GT.0 ) Pltout(i,mph) = phi
        IF ( mfa.GT.0 ) Pltout(i,mfa) = 1.D0/Oxfl
        IF ( meq.GT.0 ) Pltout(i,meq) = Eqrat
      ENDDO
      IF ( Siunit ) THEN
        pfactor = 1.D0
        fp = 'P, BAR'
        vnum = 1.D05
        frh = 'RHO, KG/CU M'
        fh = 'H, KJ/KG'
        fu = 'U, KJ/KG'
        fgi = 'G, KJ/KG'
        fs = 'S, KJ/(KG)(K)'
        fc = 'Cp, KJ/(KG)(K)'
      ELSE
        pfactor = 1.D0/1.01325D0
        fp = 'P, ATM'
        vnum = 100.D0
        frh = 'RHO, G/CC'
        fh = 'H, CAL/G'
        fu = 'U, CAL/G'
        fgi = 'G, CAL/G'
        fs = 'S, CAL/(G)(K)'
        fc = 'Cp, CAL/(G)(K)'
      ENDIF
      Fmt(4) = Fmt(6)
C PRESSURE
      CALL VARFMT(Ppp)
      DO i = 1,Npt
        X(i) = Ppp(i)*pfactor
        IF ( Nplt.NE.0.AND.i.GT.ione ) THEN
          IF ( mp.GT.0 ) Pltout(i+Iplt-ione,mp) = X(i)
          IF ( mt.GT.0 ) Pltout(i+Iplt-ione,mt) = Ttt(i)
        ENDIF
      ENDDO
      WRITE (IOOUT,Fmt) fp,(X(j),j=1,Npt)
C TEMPERATURE
      Fmt(4) = '13'
      Fmt(5) = ' '
      Fmt(7) = '2,'
      WRITE (IOOUT,Fmt) 'T, K            ',(Ttt(j),j=1,Npt)
C DENSITY
      DO i = 1,Npt
        IF ( Vlm(i).NE.0. ) X(i) = vnum/Vlm(i)
        IF ( Nplt.NE.0.AND.i.GT.ione.AND.mrho.GT.0 )
     &       Pltout(i+Iplt-ione,mrho) = X(i)
      ENDDO
      CALL EFMT(Fmt(4),frh,X)
C ENTHALPY
      DO i = 1,Npt
        X(i) = Hsum(i)*R
        IF ( Nplt.NE.0.AND.i.GT.ione.AND.mh.GT.0 )
     &       Pltout(i+Iplt-ione,mh) = X(i)
      ENDDO
      Fmt(4) = Fmt(6)
      CALL VARFMT(X)
      WRITE (IOOUT,Fmt) fh,(X(j),j=1,Npt)
C INTERNAL ENERGY
      DO i = 1,Npt
        X(i) = (Hsum(i)-Ppp(i)*Vlm(i)/Rr)*R
        IF ( Nplt.NE.0.AND.i.GT.ione.AND.mie.GT.0 )
     &       Pltout(i+Iplt-ione,mie) = X(i)
      ENDDO
      CALL VARFMT(X)
      WRITE (IOOUT,Fmt) fu,(X(j),j=1,Npt)
C GIBBS ENERGY
      DO i = 1,Npt
        X(i) = (Hsum(i)-Ttt(i)*Ssum(i))*R
        IF ( Nplt.NE.0.AND.i.GT.ione ) THEN
          IF ( mg.GT.0 ) Pltout(i+Iplt-ione,mg) = X(i)
          IF ( mm.GT.0 ) Pltout(i+Iplt-ione,mm) = Wm(i)
          IF ( mmw.GT.0 ) Pltout(i+Iplt-ione,mmw) = 1.D0/Totn(i)
          IF ( ms.GT.0 ) Pltout(i+Iplt-ione,ms) = Ssum(i)*R
          IF ( mcp.GT.0 ) Pltout(i+Iplt-ione,mcp) = Cpr(i)*R
          IF ( mgam.GT.0 ) Pltout(i+Iplt-ione,mgam) = Gammas(i)
          IF ( mdvt.GT.0 ) Pltout(i+Iplt-ione,mdvt) = Dlvtp(i)
          IF ( mdvp.GT.0 ) Pltout(i+Iplt-ione,mdvp) = Dlvpt(i)
        ENDIF
      ENDDO
      CALL VARFMT(X)
      WRITE (IOOUT,Fmt) fgi,(X(j),j=1,Npt)
C ENTROPY
      Fmt(4) = '13'
      Fmt(5) = ' '
      Fmt(7) = '4,'
      WRITE (IOOUT,Fmt) fs,(Ssum(j)*R,j=1,Npt)
      WRITE (IOOUT,99009)
C MOLECULAR WEIGHT
      Fmt(7) = '3,'
      WRITE (IOOUT,Fmt) 'M, (1/n)        ',(Wm(j),j=1,Npt)
      IF ( .NOT.Gonly ) WRITE (IOOUT,Fmt) 'MW, MOL WT      ',
     &                                (1.D0/Totn(j),j=1,Npt)
C (DLV/DLP)T
      Fmt(7) = '5,'
      IF ( Eql ) WRITE (IOOUT,Fmt) '(dLV/dLP)t      ',(Dlvpt(j),j=1,Npt)
C (DLV/DLT)P
      Fmt(7) = '4,'
      IF ( Eql ) WRITE (IOOUT,Fmt) '(dLV/dLT)p      ',(Dlvtp(j),j=1,Npt)
C HEAT CAPACITY
      WRITE (IOOUT,Fmt) fc,(Cpr(j)*R,j=1,Npt)
C GAMMA(S)
      Fmt(7) = '4,'
      WRITE (IOOUT,Fmt) 'GAMMAs          ',(Gammas(j),j=1,Npt)
C SONIC VELOCITY
      Fmt(7) = '1,'
      DO i = 1,Npt
        Sonvel(i) = (Rr*Gammas(i)*Ttt(i)/Wm(i))**.5
        IF ( Nplt.NE.0.AND.i.GT.ione.AND.mson.GT.0 )
     &       Pltout(i+Iplt-ione,mson) = Sonvel(i)
      ENDDO
      WRITE (IOOUT,Fmt) 'SON VEL,M/SEC   ',(Sonvel(j),j=1,Npt)
      RETURN
C***********************************************************************
      ENTRY OUT3
      tra = 5.D-6
      IF ( Trace.NE.0. ) tra = Trace
C MASS OR MOLE FRACTIONS 
      IF ( Massf ) THEN
        mamo = 'MASS'
      ELSE
        mamo = 'MOLE'
      ENDIF
      IF ( Eql ) THEN
        WRITE (IOOUT,99010) mamo
        notuse = 0
        DO k = 1,Ngc
          kok = .TRUE.
          IF ( k.GT.Ng.AND.k.LT.Ngc.AND.Prod(k).EQ.Prod(k+1) ) THEN
            kok = .FALSE.
            im = 0
            GOTO 120
          ENDIF
          DO m = 1,Nplt
            im = 0
            IF ( Pltvar(m).EQ.Prod(k).OR.'*'//Pltvar(m).EQ.Prod(k) )
     &           THEN
              im = m
              GOTO 120
            ENDIF
          ENDDO
 120      kin = 0
          DO i = 1,Npt
            IF ( Massf ) THEN
              tem = Mw(k)
            ELSE
              tem = 1.D0/Totn(i)
            ENDIF
            IF ( k.LE.Ng ) THEN
              X(i) = En(k,i)*tem
            ELSE
              IF ( Prod(k).NE.Prod(k-1) ) X(i) = 0.D0
              IF ( En(k,i).GT.0.D0 ) X(i) = En(k,i)*tem
            ENDIF
            IF ( Nplt.NE.0.AND.i.GT.ione.AND.im.GT.0 )
     &           Pltout(i+Iplt-ione,im) = X(i)
            IF ( kok.AND.X(i).GE.tra ) kin = 1
          ENDDO
          IF ( kin.EQ.1 ) THEN
            IF ( Trace.EQ.0. ) THEN
              WRITE (IOOUT,99011) Prod(k),(X(i),i=1,Npt)
            ELSE
              CALL EFMT(Fmt(4),Prod(k),X)
            ENDIF
            IF ( Prod(k).EQ.Omit(notuse) ) notuse = notuse - 1
          ELSEIF ( Prod(k).NE.Prod(k-1) ) THEN
            notuse = notuse + 1
            Omit(notuse) = Prod(k)
          ENDIF
        ENDDO
      ENDIF
      WRITE (IOOUT,99012) Tg(4)
      IF ( .NOT.Short ) THEN
        WRITE (IOOUT,99013) mamo,tra
        WRITE (IOOUT,99014) (Omit(i),i=1,notuse)
      ENDIF
      IF ( .NOT.Moles ) WRITE (IOOUT,99015)
      GOTO 200
C***********************************************************************
      ENTRY OUT4
      WRITE (IOOUT,99009)
      WRITE (IOOUT,99016)
      IF ( Siunit ) THEN
        WRITE (IOOUT,99018)
      ELSE
        WRITE (IOOUT,99017)
      ENDIF
C TRANSPORT PROPERTIES
      Fmt(4) = Fmt(6)
      IF ( Nplt.GT.0 ) THEN
        DO i = 1,Npt
          IF ( i.GT.ione ) THEN
            IF ( mvis.GT.0 ) Pltout(i+Iplt-ione,mvis) = Vis(i)
            IF ( mcond.GT.0 ) Pltout(i+Iplt-ione,mcond) = Coneql(i)
            IF ( mpn.GT.0 ) Pltout(i+Iplt-ione,mpn) = Preql(i)
            IF ( mcondf.GT.0 ) Pltout(i+Iplt-ione,mcondf) = Confro(i)
            IF ( mpnf.GT.0 ) Pltout(i+Iplt-ione,mpnf) = Prfro(i)
          ENDIF
        ENDDO
      ENDIF
      CALL VARFMT(Vis)
      WRITE (IOOUT,Fmt) 'VISC,MILLIPOISE',(Vis(j),j=1,Npt)
      Fmt(4) = '13'
      Fmt(5) = ' '
      Fmt(7) = '4,'
      IF ( Eql ) THEN
        WRITE (IOOUT,99019)
C SPECIFIC HEAT
        WRITE (IOOUT,Fmt) fc,(Cpeql(j),j=1,Npt)
C CONDUCTIVITY
        WRITE (IOOUT,Fmt) 'CONDUCTIVITY    ',(Coneql(j),j=1,Npt)
C PRANDTL NUMBER
        WRITE (IOOUT,Fmt) 'PRANDTL NUMBER  ',(Preql(j),j=1,Npt)
      ENDIF
      WRITE (IOOUT,99020)
C SPECIFIC HEAT
      WRITE (IOOUT,Fmt) fc,(Cpfro(j),j=1,Npt)
C CONDUCTIVITY
      WRITE (IOOUT,Fmt) 'CONDUCTIVITY    ',(Confro(j),j=1,Npt)
C PRANDTL NUMBER
      WRITE (IOOUT,Fmt) 'PRANDTL NUMBER  ',(Prfro(j),j=1,Npt)
 200  RETURN
99001 FORMAT (' CASE = ',a15)
99002 FORMAT (/13X,'REACTANT',20x,a11,'      ENERGY',6x,'TEMP')
99003 FORMAT (57X,' CAL/MOL ',6x,'K')
99004 FORMAT (57X,'KJ/KG-MOL',6x,'K')
99005 FORMAT (42X,'(SEE NOTE)      CAL/MOL       K  ')
99006 FORMAT (42X,'(SEE NOTE)     KJ/KG-MOL      K  ')
99007 FORMAT (1x,a8,4x,a15,11x,f12.7,f14.3,f11.3)
99008 FORMAT (/' O/F=',F11.5,2X,'%FUEL=',F10.6,2X,'R,EQ.RATIO=',F9.6,2X,
     &        'PHI,EQ.RATIO=',F9.6)
99009 FORMAT ()
99010 FORMAT (/1x,A4,' FRACTIONS'/)
99011 FORMAT (1x,A15,F9.5,12F9.5)
99012 FORMAT (/'  * THERMODYNAMIC PROPERTIES FITTED TO',F7.0,'K')
99013 FORMAT (/'    PRODUCTS WHICH WERE CONSIDERED BUT WHOSE ',A4,
     &        ' FRACTIONS',/'    WERE LESS THAN',1PE13.6,
     &        ' FOR ALL ASSIGNED CONDITIONS'/)
99014 FORMAT (5(1x,A15))
99015 FORMAT (/' NOTE. WEIGHT FRACTION OF FUEL IN TOTAL FUELS AND OF',
     &        ' OXIDANT IN TOTAL OXIDANTS')
99016 FORMAT (' TRANSPORT PROPERTIES (GASES ONLY)')
99017 FORMAT ('   CONDUCTIVITY IN UNITS OF MILLICALORIES/(CM)(K)(SEC)'/)
99018 FORMAT ('   CONDUCTIVITY IN UNITS OF MILLIWATTS/(CM)(K)'/)
99019 FORMAT (/'  WITH EQUILIBRIUM REACTIONS'/)
99020 FORMAT (/'  WITH FROZEN REACTIONS'/)
99021 FORMAT (/' REACTANT DENSITY=',F8.2,' KG/CU M')
99022 FORMAT (/' REACTANT DENSITY=',F8.4,' G/CC')
      END
