      SUBROUTINE SEARCH
C***********************************************************************
C SEARCH THERMO.LIB FOR THERMO DATA FOR SPECIES TO BE CONSIDERED.
C***********************************************************************
      IMPLICIT NONE
      INCLUDE 'cea.inc'
C LOCAL VARIABLES
      CHARACTER*16 bin(2,40),pure(6),spece(2)
      CHARACTER*6 date(MAXNGC)
      CHARACTER*2 el(5)
      CHARACTER*15 sub
      INTEGER i,i5,ifaz,ii,ir,itot,j,jj(2),jk,k,lineb,nall,ne,nint,
     &        npure,nrec,ntgas,ntot
      REAL*8 b(5),t1,t2,thermo(9,3),trdata(36)
      SAVE b,bin,date,el,i,i5,ifaz,ii,ir,itot,j,jj,jk,k,lineb,nall,ne,
     &  nint,npure,nrec,ntgas,ntot,pure,spece,sub,t1,t2,thermo,trdata
C
      Nc = 0
      ne = 0
      DO i = 1,Nlm
        Jx(i) = 0
      ENDDO
      DO j = 1,MAXNGC
        S(j) = 0.
        H0(j) = 0.
        Deln(j) = 0.
        DO i = 1,Nlm
          A(i,j) = 0.
        ENDDO
      ENDDO
C READ TEMPERATURE RANGES FOR COEFFICIENTS OF GASEOUS SPECIES.
C SOME DEFINITIONS:
C   NTGAS = NUMBER OF GASEOUS SPECIES IN THERMO.LIB.
C   NTOT =  NTGAS PLUS NUMBER OF TEMPERATURE INTERVALS FOR CONDENSED.
C   NALL =  NTOT PLUS THE NUMBER OF REACTANT SPECIES IN THERMO.LIB.
C   NG =    NUMBER OF GASES WITH STORED COEFFICIENTS.
C   NC =    NUMBER OF CONDENSED INTERVALS WITH STORED COEFFICIENTS.
C   NGC =    NG + NC
C   THDATE = DATE READ FROM THERMO.INP FILE
      REWIND IOTHM
      READ (IOTHM) Tg,ntgas,ntot,nall,Thdate
      Ngc = 1
      Nc = 1
C BEGIN LOOP FOR READING SPECIES DATA FROM THERMO.LIB.
      DO 200 itot = 1,ntot
        IF ( itot.GT.ntgas ) THEN
          READ (IOTHM) sub,nint,date(Ngc),(el(j),b(j),j=1,5),Ifz(Nc),
     &                 Temp(1,Nc),Temp(2,Nc),Mw(Ngc),(Cft(Nc,k),k=1,9)
        ELSE
          READ (IOTHM) sub,nint,date(Ngc),(el(j),b(j),j=1,5),ifaz,t1,t2,
     &                 Mw(Ngc),thermo
        ENDIF
        IF ( Nonly.NE.0 ) THEN
          i = 1
 20       IF ( Prod(i).NE.sub.AND.'*'//Prod(i).NE.sub ) THEN
            i = i + 1
            IF ( i.LE.Nonly ) GOTO 20
            GOTO 200
          ELSE
            IF ( sub.EQ.Prod(Ngc-1) ) THEN
              Nonly = Nonly + 1
              DO k = Nonly,i + 1, - 1
                Prod(k) = Prod(k-1)
              ENDDO
            ELSE
              Prod(i) = Prod(Ngc)
            ENDIF
            Prod(Ngc) = sub
          ENDIF
        ELSEIF ( Nomit.NE.0 ) THEN
          DO i = 1,Nomit
            IF ( Omit(i).EQ.sub.OR.'*'//Omit(i).EQ.sub ) GOTO 200
          ENDDO
        ENDIF
        DO 50 k = 1,5
          IF ( b(k).EQ.0. ) GOTO 100
          DO i = 1,Nlm
            IF ( Elmt(i).EQ.el(k) ) THEN
              A(i,Ngc) = b(k)
              GOTO 50
            ENDIF
          ENDDO
          DO j = 1,Nlm
            A(j,Ngc) = 0.
          ENDDO
          GOTO 200
 50     CONTINUE
 100    Prod(Ngc) = sub
        IF ( itot.GT.ntgas ) THEN
          Nc = Nc + 1
          IF ( Nc.GT.MAXNC ) GOTO 400
        ELSE
          Ng = Ngc
          IF ( Ng.GT.MAXNG ) GOTO 400
          DO i = 1,3
            DO j = 1,9
              Coef(Ng,j,i) = thermo(j,i)
            ENDDO
          ENDDO
C IF SPECIES IS AN ATOMIC GAS, STORE INDEX IN JX
          IF ( b(2).EQ.0..AND.b(1).EQ.1. ) THEN
            DO i = 1,Nlm
              IF ( Elmt(i).EQ.el(1) ) THEN
                ne = ne + 1
                Jx(i) = Ngc
                Jcm(i) = Ngc
                GOTO 150
              ENDIF
            ENDDO
          ENDIF
        ENDIF
 150    Ngc = Ngc + 1
        IF ( Ngc.GT.MAXNGC ) GOTO 400
 200  CONTINUE
C FINISHED READING THERMO DATA FROM I/O UNIT IOTHM.
      Ifz(Nc) = 0
      Nc = Nc - 1
      Ngc = Ngc - 1
      Ngp1 = Ng + 1
      IF ( Ngc.LT.Nonly ) THEN
        DO k = Ngc + 1,Nonly
          WRITE (IOOUT,99001) Prod(k)
        ENDDO
      ENDIF
C FIND MISSING ELEMENTS (IF ANY) FOR COMPONENTS
      Nspx = Ngc
      IF ( ne.LT.Nlm ) THEN
        DO i = 1,Nlm
          IF ( Nspx.GT.MAXNGC ) GOTO 400
          IF ( Jx(i).EQ.0 ) THEN
            Nspx = Nspx + 1
            DO k = 1,Nlm
              A(k,Nspx) = 0.
            ENDDO
            A(i,Nspx) = 1.
            Prod(Nspx) = Elmt(i)
            DO k = 1,100
              IF ( Elmt(i).EQ.Symbol(k) ) THEN
                Mw(Nspx) = Atmwt(k)
                Atwt(i) = Atmwt(k)
                Cp(Nspx) = 2.5D0
                GOTO 210
              ENDIF
            ENDDO
 210        Jx(i) = Nspx
            Jcm(i) = Nspx
          ENDIF
        ENDDO
      ENDIF
C ARE ALL ELEMENTS IN PRODUCT SPECIES?
      DO 300 i = 1,Nlm
        DO j = 1,Ngc
          IF ( A(i,j).NE.0. ) GOTO 300
          ii = i
        ENDDO
        WRITE (IOOUT,99002) Elmt(ii)
        Ngc = 0
        GOTO 600
 300  CONTINUE
C WRITE POSSIBLE PRODUCT LIST
      IF ( .NOT.Short ) THEN
        WRITE (IOOUT,99003) Thdate
        DO i = 1,Ngc,3
          i5 = i + 2
          IF ( Ngc.LT.i5 ) i5 = Ngc
          WRITE (IOOUT,99004) (date(j),Prod(j),j=i,i5)
        ENDDO
      ENDIF
      GOTO 600
 400  WRITE (IOOUT,99005)
      Ngc = 0
      GOTO 600
C SEARCH FOR TRANSPORT PROPERTIES FOR THIS CHEMICAL SYSTEM
      ENTRY READTR
      REWIND IOTRN
      REWIND IOSCH
      Ntape = 0
      npure = 0
      lineb = 1
      IF ( .NOT.Short ) WRITE (IOOUT,99006)
      READ (IOTRN) nrec
      DO ir = 1,nrec
        READ (IOTRN) spece,trdata
        k = 1
 450    DO j = 1,Ng
          IF ( spece(k).EQ.Prod(j).OR.'*'//spece(k).EQ.Prod(j) ) THEN
            jj(k) = j
            IF ( k.EQ.2 ) THEN
C STORE NAMES FOR BINARIES IN BIN ARRAY.
              DO k = 1,2
                bin(k,lineb) = spece(k)
              ENDDO
              lineb = lineb + 1
              GOTO 500
            ELSE
              jj(2) = j
              IF ( spece(2).EQ.' ' ) THEN
C WRITE NAMES FOR PURE SPECIES.
                npure = npure + 1
                pure(npure) = spece(1)
                GOTO 500
              ELSE
                k = 2
                GOTO 450
              ENDIF
            ENDIF
          ENDIF
        ENDDO
        GOTO 550
 500    WRITE (IOSCH) jj,trdata
        Ntape = Ntape + 1
 550    IF ( npure.NE.0.AND.(npure.GE.6.OR.ir.GE.nrec) ) THEN
          IF ( .NOT.Short ) WRITE (IOOUT,99007) (pure(jk),jk=1,npure)
          npure = 0
        ENDIF
      ENDDO
      lineb = lineb - 1
      IF ( .NOT.Short ) THEN
        WRITE (IOOUT,99008)
        DO j = 1,lineb
          WRITE (IOOUT,99009) (bin(i,j),i=1,2)
        ENDDO
      ENDIF
      WRITE (IOOUT,99010)
 600  RETURN
99001 FORMAT (/' WARNING!!  ',A15,' NOT A PRODUCT IN thermo.lib FILE ',
     &        '(SEARCH)')
99002 FORMAT (/' PRODUCT SPECIES CONTAINING THE ELEMENT',A3,' MISSING',
     &        //,13x,'FATAL ERROR (SEARCH)')
99003 FORMAT (/2x,'SPECIES BEING CONSIDERED IN THIS SYSTEM',
     &        /' (CONDENSED PHASE MAY HAVE NAME LISTED SEVERAL TIMES)',
     &        /'  LAST thermo.inp UPDATE: ',A10,/)
99004 FORMAT (3(2X,A6,2X,A15))
99005 FORMAT (/' INSUFFICIENT STORAGE FOR PRODUCTS-SEE RP-1311,',
     &        /'   PART 2, PAGE 39. (SEARCH)')
99006 FORMAT (/' SPECIES WITH TRANSPORT PROPERTIES'//8X,'PURE SPECIES'/)
99007 FORMAT (4(2x,A16))
99008 FORMAT (/'     BINARY INTERACTIONS'/)
99009 FORMAT (5X,2A16)
99010 FORMAT ()
      END
