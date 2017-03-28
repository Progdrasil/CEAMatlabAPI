      SUBROUTINE INFREE(Readok,Cin,Ncin,Lcin,Dpin)
C***********************************************************************
C FREE-FORM READ FOR CEA.  READS AND DECIPHERS DATA FOR ONE DATASET.
C
C DEFINITIONS:
C   CH1  - INDIVIDUAL CHARACTERS IN RECORD, MAXIMUM 132.
C   NCH1 - COLUMN NUMBER FOR THE LAST NON-BLANK CHARACTER IN RECORD.
C   NCIN - NUMBER OF VARIABLES IN DATASET.
C   CIN  - CHARACTER STRINGS IN DATASET. MAXIMUM 15 CHARACTERS.
C   LCIN - NEG. LENGTH OF LITERALS.  FOR NUMERICS, INDEX OF PREVIOUS
C          LITERAL.  ZERO FOR UNACCEPTIBLE VARIABLES.  VARIABLE
C          FOLLOWING "CASE" IS ALWAYS ASSUMED TO BE LITERAL.
C   NB   - NUMBER OF DELIMITERS IN STRING.
C   NX   - NUMBER OF CHARACTERS IN STRING.
C   DPIN - NUMERICS AS DOUBLE PRECISION VARIABLE.
C   CNUM - CHARACTER STRING REPRESENTING DATASET NUMBERS. MAXIMUM 24
C          CHARACTERS.
C***********************************************************************
      IMPLICIT NONE
      INCLUDE 'cea.inc'
C DUMMY ARGUMENTS
      CHARACTER*15 Cin(MAXNGC)
      INTEGER Ncin
      INTEGER Lcin(MAXNGC)
      LOGICAL Readok
      REAL*8 Dpin(MAXNGC)
C LOCAL VARIABLES
      CHARACTER*1 ch1(132),cx,nums(13)
      CHARACTER*24 cnum
      CHARACTER*3 fmtl(3)
      CHARACTER*2 numg(24)
      CHARACTER*4 w1
      INTEGER i,ich1,j,kcin,nb,nch1,nx
C
      DATA fmtl/'(g','16','.0)'/
      DATA nums/'+','-','0','1','2','3','4','5','6','7','8','9','.'/
      DATA numg/'1','2','3','4','5','6','7','8','9','10','11','12','13',
     &     '14','15','16','17','18','19','20','21','22','23','24'/
      Ncin = 1
      Lcin(1) = 0
      kcin = 0
      Dpin(1) = 0.D0
 100  nb = 1
      nx = 0
      cnum = ' '
      Cin(Ncin) = ' '
      ch1(1) = ' '
      nch1 = 1
C READ CHARACTERS, ONE AT A TIME of input file
      READ (IOINP,99001,END=500,ERR=500) ch1
C FIND FIRST AND LAST NON-BLANK CHARACTER
      DO i = 132,1, - 1
        nch1 = i
        IF ( ch1(i).NE.' '.AND.ch1(i).NE.'	' ) GOTO 200
      ENDDO
 200  DO i = 1,nch1
        ich1 = i
        IF ( ch1(i).NE.' '.AND.ch1(i).NE.'	' ) GOTO 300
      ENDDO
 300  IF ( nch1.EQ.1.OR.ch1(ich1).EQ.'#'.OR.ch1(ich1).EQ.'!' ) THEN
        WRITE (IOOUT,99002) (ch1(i),i=1,nch1)
        GOTO 100
      ENDIF
      w1 = ch1(ich1)//ch1(ich1+1)//ch1(ich1+2)//ch1(ich1+3)
C IS STRING A KEYWORD SIGNALLING START OR END OF DATASET?
      IF ( w1.EQ.'ther'.OR.w1.EQ.'tran'.OR.w1.EQ.'prob'.OR.
     &     w1.EQ.'reac'.OR.w1.EQ.'outp'.OR.w1.EQ.'omit'.OR.
     &     w1.EQ.'only'.OR.w1.EQ.'inse'.OR.w1(1:3).EQ.'end' ) THEN
        IF ( Ncin.EQ.1 ) THEN
          Cin(Ncin) = w1
          IF ( w1(1:3).EQ.'end'.OR.w1.EQ.'ther'.OR.w1.EQ.'tran' ) THEN
            WRITE (IOOUT,99002) (ch1(i),i=1,nch1)
            RETURN
          ENDIF
          ich1 = ich1 + 4
          nx = 4
          Lcin(1) = -4
        ELSE
C KEYWORD READ FOR NEXT DATASET. END PROCESSING
          BACKSPACE IOINP
          IF ( nx.EQ.0 ) Ncin = Ncin - 1
          RETURN
        ENDIF
      ELSEIF ( Ncin.EQ.1 ) THEN
        WRITE (IOOUT,99003)
        GOTO 500
      ENDIF
      WRITE (IOOUT,99002) (ch1(i),i=1,nch1)
      DO 400 i = ich1,nch1
        cx = ch1(i)
C LOOK FOR DELIMITER STRINGS
        IF ( cx.EQ.','.AND.(Lcin(Ncin).GT.0.OR.nx.EQ.0) ) cx = ' '
        IF ( cx.EQ.'='.AND.(Lcin(Ncin).LT.0.OR.nx.EQ.0) ) cx = ' '
        IF ( cx.NE.' '.AND.cx.NE.'	' ) THEN
C LOOK FOR CHARACTER STRINGS
          nx = nx + 1
          IF ( Ncin.GT.1 ) THEN
            cnum(nx:nx) = cx
            IF ( nx.LE.15 ) Cin(Ncin) = cnum
            IF ( nx.EQ.1 ) THEN
C IS THIS A NUMERIC?
              DO j = 1,13
                IF ( ch1(i).EQ.nums(j) ) THEN
                  Lcin(Ncin) = kcin
                  GOTO 310
                ENDIF
              ENDDO
              Lcin(Ncin) = -1
              kcin = Ncin
            ELSEIF ( Lcin(Ncin).LT.0 ) THEN
              Lcin(Ncin) = -nx
            ENDIF
 310        nb = 1
          ENDIF
          IF ( i.LT.nch1.OR.Lcin(Ncin).LT.0 ) GOTO 400
        ENDIF
        IF ( nb.EQ.1..AND.nx.GT.0 ) THEN
          IF ( Ncin.GT.0.AND.Lcin(Ncin).GT.0 ) THEN
C CONVERT NUMERIC CHARACTER STRINGS TO REAL*8 VARIABLES (DPIN)
            fmtl(2) = numg(MIN(24,nx))
C INTERNAL READ TO CONVERT TO NUMERIC
            READ (cnum,fmtl,ERR=320) Dpin(Ncin)
          ENDIF
          GOTO 340
 320      IF ( Cin(Ncin-1)(:4).NE.'case' ) WRITE (IOOUT,99004) Cin(i)
          Lcin(Ncin) = 0
 340      Ncin = Ncin + 1
          Cin(Ncin) = ' '
          Lcin(Ncin) = 0
          Dpin(Ncin) = 0.D0
          nx = 0
          cnum = ' '
        ENDIF
        nb = nb + 1
 400  CONTINUE
      IF ( nx.GT.0 ) THEN
        Ncin = Ncin + 1
        Lcin(Ncin) = 0
        Dpin(Ncin) = 0.D0
      ENDIF
      GOTO 100
 500  Readok = .FALSE.
      RETURN
99001 FORMAT (132A1)
99002 FORMAT (1x,80A1)
99003 FORMAT (/' FATAL ERROR IN INPUT FORMAT (INFREE)')
99004 FORMAT (/' WARNING!!  UNACCEPTABLE NUMBER ',A15,' (INFREE)')
      END
