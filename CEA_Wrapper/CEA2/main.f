C***********************************************************************
C                     P R O G R A M      C E A 2
C
C             CHEMICAL EQULIBRIUM WITH APPLICATIONS         5/21/04
C***********************************************************************
C 11/01/94 - CORRECTION IN EQLBRM  TEST FOR CONVERGENCE ON CONDENSED
C            CHANGED FROM: DABS(DELN(J))<.5D-05
C                      TO: DABS(DELN(J)/TOTN(NPT))<.5D-05
C 06/03/96 - INFREE CHANGED TO BE MORE COMPATIABLE WITH VARIOUS FORTRAN
C            COMPILERS AND TO ACCEPT NUMERICS WITH UP TO 24 CHARACTERS.
C 07/05/96 - CORRECTIONS MADE FOR COMPATIBILITY WITH OTHER COMPILERS.
C 08/29/96 - CORRECTION IN EQLBRM TO CORRECT FOR CONDENSED SPECIES
C            TEST AT TEMPERATURE INTERVAL OVERLAP POINTS.
C 09/17/96 - ADD 2 MORE VARIABLES FOR PLOTING, THE VOLUME PARTIAL
C            DERIVATIVES: (DLNV/DLNT)P AND (DLNV/DLNP)T AS DEFINED
C            IN EQUATIONS (2.50) AND (2.51) IN RP-1311,PT1. THE PLOT
C            LIST VARIABLES IN DATASET 'OUTPUT' HAVE EMBEDDED STRINGS
C            'DLNT' AND 'DLNP' RESPECTIVELY (E.G. 'DLNV/DLNT' OR
C            '(DLNV/DLNP)T'
C 10/18/96 - CORRECT INFREE SO THAT NUMERICAL INPUT ARRAYS MAY BE SPLIT
C            FROM RECORD TO RECORD.
C 12/11/96 - INPUT: ELIMINATE 'PSI' POSSIBILITY IN PCP CHECK.
C 12/12/96 - INPUT: WRITE MESSAGE WHEN VALUES OF P OR T ARE MISSING.
C 01/31/97 - MAIN: REVERSE ORDER FOR INSERTING CONDENSED SO THAT HIGHEST
C            T INTERVAL IS INSERTED.
C 01/31/97 - EQLBRM: CHECK ON INITIAL T INTERVALS FOR NON-TP PROBLEMS
C            WHEN A T ESTIMATE IS GIVEN. IF NO T EST., USE 3800K.
C 01/31/97 - THERMP: REMOVE INITIAL ESTIMATE OF 3800K WHEN NO T IS GIVEN
C 03/18/97 - INPUT: INSURE THAT INPUT VARIABLES P,T,MACH,U1,V,RHO,SUBAR,
C            SUPAR,O/F,AND PCP DO NOT EXCEED STORAGE.  WRITE MESSAGE IF
C            TOO MANY VALUES ARE LISTED IN THE INPUT FILE.
C 03/18/97 - SEARCH: IFZ(NC+1) SET TO 0. CAUSED PROBLEMS WITH CONDENSED.
C 03/20/97 - MAIN: CHECK FOR 100%FUEL WHEN OXIDANT IS GIVEN AND OMITTING
C            THE OXIDANT CHANGES THE CHEMICAL SYSTEM. FATAL ERROR MESSAGE.
C 03/21/97 - MAIN: FOR NON-TP CASES AND FIRST POINT AND NO T ESTIMATE IS
C            GIVEN, AND CONDENSED SPECIES ARE INSERTED, ESTIMATE T TO BE
C            THE LOWEST T MAXIMUM OF THE INSERTED SPECIES MINUS 10K..
C 04/22/97 - HCALC AND DETON: CEA WAS ALLOWING ALLOWING EXECUTION OF
C            DETON PROBLEMS WITH CONDENSED REACTANTS.
C 04/28/97 - EQLBRM: CORRECT PHASES OF CONDENSED PRODUCTS WERE NOT BEING
C            CHECKED AT THE BEGINING OF THE ROUTINE FOR TP PROBLEMS.
C 08/11/97 - SEARCH: PRINT MESSAGE FOR EXCEEDING MAXNGC,MAXNG, OR MAXNC.
C 09/04/97 - MAIN,SEARCH: CHECK TO SEE IF PRODUCTS CONTAIN ALL ELEMENTS.
C 10/09/97 - REACT,HCALC: CHANGE ERROR MESSAGE FOR REACTANT DATA NOT
C            FOUND IN THERMO.LIB.
C 12/22/97 - INTRODUCE VARIABLE ARRAY SIZES FOR INPUT VARIABLES P,T,V,
C            (CORRECTED AND REACTANT MIXTURES. THESE MAXIMUM VALUES ARE
C             IN 1/29/98)  PARAMETER STATEMENTS
C             MAXPV - MAX NUMBER OF PRESSURES OR DENSITIES (DEFAULT 26)
C             MAXT - MAX NUMBER OF TEMPERATURES (DEFAULT 51)
C             MAXMIX - MAX NUMBER OF MIXTURE VALUES (DEFAULT 52)
C 02/12/98 - UTHERM: CHECK FOR NUMBER OF T INTERVALS (NT) BEING > 3 FOR
C            GASES (IFAZ=0).  IF SO, EXIT WITH MESSAGE.
C 03/02/98 - INCREASE DIMENSIONS OF VARIABLES FOR PLOT DUMP - PLTVAR,
C            PLTOUT,MAX NPLT, AND FORMAT STATEMENT IN MAIN.
C 03/10/98 - CORRECTED 3/2/98 CHANGE.
C 05/20/98 - UTHERM & HCALC. CORRECTED 'AIR' CALCULATIONS ABOVE 6000K
C            IN HCALC BY USING COEFS FOR 1000-6000. CORRECT MOLAR AMTS
C            IN UTHERM.
C 12/01/98 - INFREE & INPUT. CHANGED ICH ARRAY.
C 02/12/99 - UTHERM & SEARCH. ADDED DATE TO THERMO DATA.
C            VARFMT. ELIMINATE "," BEFORE ")" IN FORMAT FMT.
C 02/16/99 - INFREE. USE "CTRL I" FOR TAB.
C 04/16/99 - NOTE: THIS VERSION CONTAINS SOME CLEANUP CHANGES NOT
C            INCLUDED IN PREVIOUS VERSIONS.
C 04/16/99 - UTHERM: WHEN SETTING 'AIR' COMPOSITION, REMOVE CHECK ON C
C            ALSO INSERT WEIGHT WITH MORE FIGURES FOR E-.
C 04/16/99 - THERMO.INP CHANGED SO THAT ALL COEFFICIENTS LISTED IN ORDER:
C             -2  -1 0 1 2 3 4 WHETHER ALL COEFFICIENTS EXIST OR NOT.
C 05/26/99 - NEW ATOMIC WEIGHTS IN BLOCKDATA. CARBON NOW 12.0107.
C 06/11/99 - EQLBRM: MANY SMALL CHANGES THAT DO NOT AFFECT CALCULATIONS.
C            (CORRECTED 7/16/99)
C 07/01/99 - COMMON INPT AND NEWOF: BRATIO IS REPLACED BY BCHECK WHERE
C            BCHECK = (BIGGEST B0I)*.000001 AND USED IN EQN.(3.6A),RP-1311.
C            BRATIO AS DESCRIBED IN SEC.3.2 IS USED IN NEWOF ONLY.
C 08/30/99 - CHANGE I/O UNITS FOR INPUT AND OUTPUT FROM 5 & 6 TO 7 & 8.
C 11/05/99 - REWORD '99000 FORM...' IN SEARCH & '99003 FORM...' IN REACT.
C 06/13/00 - INPUT. BYPASS TCEST VARIABLE IN T(I) LOOP. CEA WAS CONVERTING
C            VALUE FROM CENTIGRADE TO KELVIN AND ASSIGNING TEMPERATURE.
C 09/15/00 - MAIN. PRINT MESSAGE IF INSERT SPECIES IS NOT FOUND.
C 10/17/00 - CHANGE MAX NUMBER OF ELEMENTS IN A REACTANT FROM 6 TO 12.
C 01/16/01 - EQLBRM: PRINT MESSAGE IF TT DRIVEN DOWN TO < .2*TG(1).
C 01/16/01 - INPUT: HP IN PRINTOUT CORRECTED TO F WHEN UV = T.
C 04/23/01 - TEST FOR THERMO PROPERTY EXTRAPOLATION FOR PURE SPECIES:
C            ( TT.GE.TG(1)*.8D0.AND.TT.LE.TG(4)*1.1D0) WHERE TG(1)= 200K.
C 08/06/01 - UTHERM: CORRECT LABELLING FOR CONDENSED REACTANTS.
C 05/30/02 - ROCKET: In "IF ( ipp.LT.Npp ) ipp = ipp-1", change LT to LE
C 06/21/02 - UPDATED FORTRAN
C 08/19/02 - INPUT: DISCONTINUE RUN FOR ILLEGAL EQUIVALENCE RATIO.
C 09/20/02 - REACT: FOR REACTANT TEMPS OUT OF RANGE, PRINT T RANGE.
C 10/18/02 - BLOCKDATA: New atomic weights inc. N=14.0067, S=32.065,...
C  1/17/03 - cea.inc: Change MAXTR from 30 to 40 for large systems.
C  3/31/03 - Declare I/O units 7 & 8 to be IOINP & IOOUT in cea.inc  
C  4/03/03 - ROCKET: Added check for assigned ae/at < 1.
C 10/22/03 - OUT2: Plot check on p,s,and g needs check on first 2 letters.
C 11/14/03 - OUT: Add cond...fz and pran...fz to .plt for all problems.
C 11/18/03 - INPUT: If exploded formula, Energy = ' ',not 'lib'.
C 11/25/03 - DETON: Extra line for .plt data with multiple DETON calls.
C 01/30/04 - ROCKET: To correct points printed after fac, change 
C            IF(ipp.LE.Npp)ipp=ipp-1 to IF(ipp.LT.Npp.OR.Npp.EQ.4) ...
C 02/05/04 - Chg. numbered ENDDO's to CONTINUE's for Watson compiler.
C 05/21/04 - Added labels to columns in .plt file
C***********************************************************************
      IMPLICIT NONE
      INCLUDE 'cea.inc'
C LOCAL VARIABLES
      CHARACTER*15 ensert(20)
      CHARACTER*200 infile,ofile
      CHARACTER*196 prefix
      LOGICAL caseok,ex,readok
      INTEGER i,inc,iof,j,ln,n
      INTEGER INDEX
      REAL*8 xi,xln
      REAL*8 DLOG
      SAVE caseok,ensert,ex,i,inc,infile,iof,j,ln,n,ofile,prefix,readok, !saves variables, but where?
     &  xi,xln
C
C     WRITE (*,99001) ! Ask for input file name
C     READ (*,99002) prefix ! Get input file name from user
      prefix = 'wrapper'
      ln = INDEX(prefix,' ') - 1 
      infile = prefix(1:ln)//'.inp' !add file suffix for input file
      ofile = prefix(1:ln)//'.dat'  !add file suffix for output file
      Pfile = prefix(1:ln)//'.plt'  !add file suffix for... i dont know
      INQUIRE (FILE=infile,EXIST=ex)! Get information on opened files
      IF ( .NOT.ex ) THEN !if file does not exist then
        PRINT *,infile,' DOES NOT EXIST' !tell user the file doesnt exist
        GOTO 400  !stop the program
      ENDIF
      OPEN (IOINP,FILE=infile,STATUS='old',FORM='formatted')! open input file given by user
      OPEN (IOOUT,FILE=ofile,STATUS='unknown',FORM='formatted') !create and open output file
      OPEN (IOSCH,STATUS='scratch',FORM='unformatted')! FIND OUT
      OPEN (IOTHM,FILE='thermo.lib',FORM='unformatted')!open termodynamics library 
      OPEN (IOTRN,FILE='trans.lib',FORM='unformatted')!Open trans library (WHAT IS TRANS LIBRARY?)
      WRITE (IOOUT,99006) !make a seperation line in output
      WRITE (IOOUT,99007) !write authors and name of application in output
      WRITE (IOOUT,99006) !make another seperation line in output
      readok = .TRUE.   ! Was able to read input file
      Newr = .FALSE.  
 100  Iplt = 0
      Nplt = 0
      CALL INPUT(readok,caseok,ensert) !Subroutine defined at line 2118
      IF ( caseok.AND.readok ) THEN
        DO iof = 1,Nof
          IF ( Oxf(iof).EQ.0..AND.B0p(1,1).NE.0. ) THEN
            DO i = 1,Nlm
              IF ( B0p(i,1).EQ.0..OR.B0p(i,2).EQ.0. ) THEN
                WRITE (IOOUT,99008)
                GOTO 200
              ENDIF
            ENDDO
          ENDIF
        ENDDO
        IF ( Ions ) THEN
          IF ( Elmt(Nlm).NE.'E' ) THEN
            Nlm = Nlm + 1
            Elmt(Nlm) = 'E'
            B0p(Nlm,1) = 0.
            B0p(Nlm,2) = 0.
          ENDIF
        ELSEIF ( Elmt(Nlm).EQ.'E' ) THEN
          Nlm = Nlm - 1
        ENDIF
        DO n = 1,Nreac
          Jray(n) = 0
        ENDDO
        CALL SEARCH
        IF ( Ngc.EQ.0 ) GOTO 300
        Newr = .FALSE.
        IF ( Trnspt ) CALL READTR
C INITIAL ESTIMATES
        Npr = 0
        Gonly = .TRUE.
        Enn = .1D0
        Ennl = -2.3025851
        Sumn = Enn
        xi = Ng
        IF ( xi.EQ.0. ) xi = 1.
        xi = Enn/xi
        xln = DLOG(xi)
        DO inc = 1,Nc
          j = Ng + inc
          En(j,1) = 0.D0
          Enln(j) = 0.D0
        ENDDO
        DO j = 1,Ng
          En(j,1) = xi
          Enln(j) = xln
        ENDDO
        IF ( Nc.NE.0.AND.Nsert.NE.0 ) THEN
          DO 120 i = 1,Nsert
            DO j = Ngc,Ngp1, - 1
              IF ( Prod(j).EQ.ensert(i) ) THEN
                Npr = Npr + 1
                Jcond(Npr) = j
                IF ( .NOT.Short ) WRITE (IOOUT,99003) Prod(j)
                GOTO 120
              ENDIF
            ENDDO
            WRITE (IOOUT,99004) ensert(i)
 120      CONTINUE
        ENDIF
        IF ( Rkt ) THEN
          CALL ROCKET
        ELSEIF ( Tp.OR.Hp.OR.Sp ) THEN
          CALL THERMP
        ELSEIF ( Detn ) THEN
          CALL DETON
        ELSEIF ( Shock ) THEN
          CALL SHCK
        ENDIF
        IF ( Nplt.GT.0 ) THEN
          OPEN (IOPLT,FILE=Pfile,FORM='formatted')
          WRITE (IOPLT,99009) (Pltvar(j),j=1,Nplt)
          DO i = 1,Iplt
            WRITE (IOPLT,99005) (Pltout(i,j),j=1,Nplt)
          ENDDO
          WRITE (IOPLT,99009) (Pltvar(j),j=1,Nplt)
        ENDIF
      ENDIF
 200  IF ( readok ) GOTO 100
 300  CLOSE (IOINP)
      CLOSE (IOOUT)
      CLOSE (IOSCH)
      CLOSE (IOTHM)
      CLOSE (IOTRN)
      CLOSE (IOPLT)
 400  STOP
99001 FORMAT (//' ENTER INPUT FILE NAME WITHOUT .inp EXTENSION.'/ 
     &        '   THE OUTPUT FILES FOR LISTING AND PLOTTING WILL HAVE',/
     &       ' THE SAME NAME WITH EXTENSIONS .out AND .plt RESPECTIVELY'
     &       //)
99002 FORMAT (a)
99003 FORMAT (1X,A16,'INSERTED')
99004 FORMAT (/' WARNING!!!',A16,'NOT FOUND FOR INSERTION')
99005 FORMAT (1x,1p,20E12.4)
99006 FORMAT (/' ***************************************************',
     &        '****************************')
99007 FORMAT (/,9x,'NASA-GLENN CHEMICAL EQUILIBRIUM PROGRAM CEA2,',
     &        ' MAY 21, 2004',/19x,'BY  BONNIE MCBRIDE', 
     &        ' AND SANFORD GORDON',/5x,
     &        ' REFS: NASA RP-1311, PART I, 1994',
     &        ' AND NASA RP-1311, PART II, 1996')
99008 FORMAT (/,'OXIDANT NOT PERMITTED WHEN SPECIFYING 100% FUEL(main)')
99009 FORMAT ('#',2x,20A12)
      END
