      SUBROUTINE SFINT(ICODE, LCODE, KPOWER, IORDER, LIMIT,
     *  MAXLBP, LISTBP, LENBP, IFAULT)
C
C         AS R82 (REMARK ON AS65) APPL.STATIST. (1990), VOL.39, NO.1
C
C         Expands a structure formula into a list of binary integers
C
      INTEGER IFAULT,IORDER,KPOWER,LCODE,LENBP,LIMIT,MAXLBP
      INTEGER ICODE(LCODE),LISTBP(MAXLBP)
      INTEGER I,IBP,IBS,IJB,ILS,IORD,IORJ,IP,IPRONO,IREL,IRS,
     1        J,JB,JBP,JC,JLFAC,JLSUM,JLSYM,JMFAC,JMSUM,JMSYM,JOP,JP,
     2        JRFAC,JRSUM,JRSYM,JTFAC,JTSUM,JTSYM,K,LLSUM,M,NPOWER,NS
      LOGICAL LCHECK,RCHECK
C
C         Set NBI so that 2**NBI-1 is a large legal positive integer
C
      INTEGER NBI
      PARAMETER (NBI=31)
C
C         Set codes for K-symbols and L according to number of symbols
C         (PNSYM) and number N of precedence classes (PNPREC), also
C         codes for symbols
C
      INTEGER PDOT,PK0,PK2NP1,PK2NP2,PL,PLBRKT,PMINSL,PMINST,PMINUS,
     1        PNPREC,PNSYM,PPLUS,PRBRKT,PSLASH,PSTAR,PSTAR2
      PARAMETER (PSTAR2=1,PDOT=2,PSLASH=3,PSTAR=4,PPLUS=5,PMINUS=6,
     1 PMINSL=7,PMINST=8,PLBRKT=9,PRBRKT=10, PNPREC=5,PNSYM=10,
     2 PK0=PNSYM+1,PK2NP1=PK0+2*PNPREC+1,PK2NP2=PK2NP1+1,PL=PK2NP2+1)
C
C         Initialise compact precedence table and set auxiliary
C         quantities
C
      INTEGER KILLEG,KOPEN,KCON,KCLOSE
      PARAMETER (KILLEG=0,KOPEN=1,KCON=2,KCLOSE=3)
      INTEGER IPRTAB(5,5), IENTRY(PNSYM+2*PNPREC+4),
     1        IPREC(PNSYM+2*PNPREC+4)
C
      DATA IPRTAB /2*0,-1,-1,1, -3,3,3*0, 2*0,3*1, -3,3,3*0, 2*3,3*0/
      DATA IENTRY /8*4, 3, 2, 13*1, 5/
      DATA IPREC /2, 4, 6, 8, 4*10, 12, 11,
     1      0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13/
C
C     LITERAL =   **  .  /  *  +  - -/ -*   (  )    FACTOR
C     SYMBOL  =   O2 O4 O6 O8 10 10 10 10   (  )    K0 K1 ... K12(=K2N+2)   L
C     CODE    =    1  2  3  4  5  6  7  8   9 10    11 12 ...  23          24
C     IENTRY  =    4  4  4  4  4  4  4  4   3  2     1  1 ...   1           5
C     IPREC   =    2  4  6  8 10 10 10 10  12 11     0  1 ...  12          13
C
C-------------------INITIALISATION--------------------------------------
C         IORD=highest order of retained terms
C         JC=pointer to next input symbol in icode
C         JB=pointer to top of sum stack LISTBP(1...)
C         JP=pointer to top of symbol stack LISTBP(...MAXLBP) (top.down)
C         containing 3-integer cells (symbol code,pointer to SUM,FAC)
C
      IFAULT = 0
      LENBP = 0
      IORD = LIMIT
      IF (IORD.LE.0) IORD = NBI
      JC = 0
      IF (LCODE) 1001, 999, 10
   10 JB = 0
      JP = MAXLBP + 1
C
C------------------ SYNTAX ANALYSER -----------------------------------
C-------- Stack limit symbol L and/or next ICODE symbol
C         Stack -input code for factor
C
  100 I = PL
      IF (JC.LT.1.OR.JC.GE.LCODE) GOTO 115
  110 JC = JC + 1
      I = ICODE(JC)
      IF (I.EQ.0.OR.I.GT.PNSYM) GOTO 1004
  115 JP = JP - 3
      IF (JB.GE.JP) GOTO 1002
      LISTBP(JP + 2) = -I
      LISTBP(JP + 1) = 0
      IF (I.LT.0) I = PK0
      LISTBP(JP) = I
      IF (JC.EQ.0) GOTO 110
C
C-------- SEARCH FOR PRODUCTION
C         If closing relation between top two stacked symbols scan for
C         matching opening relation to delimit possible rightpart
C
  120 IP = JP
      IBS = 0
  130 IP = IP + 3
C
C         IREL=relation between leftsymbol at IP and rightsymbol at IP-3
C
      ILS = LISTBP(IP)
      IRS = LISTBP(IP - 3)
      IREL = IPRTAB(IENTRY(ILS), IENTRY(IRS))
      IF (IREL.GE.0) GOTO 140
      IREL = -IREL
      M = IPREC(ILS) - IPREC(IRS)
      IF (IREL.EQ.KOPEN) M = 1 - M
      IF (M.EQ.0) IREL = KCON
      IF (M.GT.0) IREL = KILLEG
C
C         Branch on IBS=0/1 for terminal/earlier relation
C
  140 IF (IBS.EQ.1) GOTO 170
C
C         Terminal relation
C
      IF (IREL.EQ.KILLEG) GOTO 1005
      IBS = 1
      IF (IREL - KCLOSE) 100, 130, 100
C
C         Earlier relation
C
  170 IF (IREL.NE.KOPEN) GOTO 130
C
C-------- IDENTIFY PRODUCTION
C         Leftpart  rightpart                  Semantic action  IPRONO
C         K1        K0                         load             1
C         K1        ( K2N+1 )                  copy             2
C         KI        KI OI KI-1 I=2,4,...,2N    dyadic operation 3
C         KI        KI-1       I=1,2,...,2N+1         last two productions
C         K2N+2   L K2N+1 L                           are shortcircuited
C         Identify production, if any, with rightpart matching the NS
C         stacked symbols between the opening and closing relations
C         Copy terminal + rightpart cells
C
      NS = (IP - JP - 3) /3
      IF (NS.NE.1.AND.NS.NE.3) GOTO 1005
      JTSYM = LISTBP(JP)
      JTSUM = LISTBP(JP + 1)
      JTFAC = LISTBP(JP + 2)
      JRSYM = LISTBP(JP + 3)
      JRSUM = LISTBP(JP + 4)
      JRFAC = LISTBP(JP + 5)
      IPRONO = 4
      IF (NS.EQ.3) GOTO 180
      IF (JRSYM.EQ.PK0) IPRONO = 1
      GOTO 190
  180 JMSYM = LISTBP(JP + 6)
      JMSUM = LISTBP(JP + 7)
      JMFAC = LISTBP(JP + 8)
      JLSYM = LISTBP(JP + 9)
      JLSUM = LISTBP(JP + 10)
      JLFAC = LISTBP(JP + 11)
      IF (JMSYM.EQ.PK2NP1.AND.JLSYM.EQ.PLBRKT.AND.JRSYM.EQ.PRBRKT)
     *  IPRONO = 2
      I = IPREC(JMSYM)
      IF (IENTRY(JMSYM).EQ.4.AND.
     *    JLSYM.EQ.PK0 + I.AND.JRSYM.EQ.PK0 + I - 1) IPRONO = 3
C
C         Branch to semantics section
C
  190 GOTO (500, 600, 700, 1005), IPRONO
C
C-------- UPDATE STACK
C         Set JRSYM=K-symbol in leftpart. anticipate and shortcircuit
C         chain of 1-symbol productions by promoting JRSYM until it
C         concatenates with a stacked neighbouring symbol
C
  200 JRSYM = PK0 + MIN(IPREC(LISTBP(IP)) - 1, IPREC(LISTBP(JP)))
C
C         If leftpart=K2N+2 (formula) set LENBP and return
C         else stack terminal + leftpart cells
C
      IF (JRSYM.LT.PK2NP2) GOTO 210
      LENBP = JB
      GOTO 999
  210 JP = IP - 6
      LISTBP(JP) = JTSYM
      LISTBP(JP + 1) = JTSUM
      LISTBP(JP + 2) = JTFAC
      LISTBP(JP + 3) = JRSYM
      LISTBP(JP + 4) = JRSUM
      LISTBP(JP + 5) = JRFAC
      GOTO 120
C
C------------------ SEMANTICS ------------------------------------------
C-------- LOAD
C
  500 IF (JRFAC.GT.KPOWER) THEN
        JRFAC = - MIN(IORD, JRFAC - KPOWER)
      ELSE
        IF (JRFAC.GT.NBI) GOTO 1003
        JRFAC = 2 ** (JRFAC - 1)
      ENDIF
      JRSUM = JB
      JB = JB + 1
      IF (JB.GE.JP) GOTO 1002
      LISTBP(JB) = JRFAC
      GOTO 200
C
C-------- COPY
  600 JRSYM = JMSYM
      JRSUM = JMSUM
      JRFAC = JMFAC
      GOTO 200
C
C-------- DYADIC OPERATION
C
  700 JOP = JMSYM
      IJB = JB
      LLSUM = JRSUM
      RCHECK = (JOP.EQ.PSTAR2).OR.(JOP.EQ.PSTAR)
      LCHECK = (JOP.EQ.PSLASH).OR.(JOP.EQ.PSTAR)
C
C           **   .    /    *    +    -    -/   -*
      GOTO (702, 730, 710, 730, 780, 720, 720, 720), JOP
  702 NPOWER = - LISTBP(JRSUM + 1)
      IF(NPOWER.LE.0) GOTO 1005
      JB = JRSUM
      JRFAC = JLFAC
      IF (NPOWER.EQ.1) GOTO 775
      JRSUM = JLSUM
      IJB = LLSUM
      GOTO 730
  710 JB = JRSUM
      IBP = JLFAC
      I = JRSUM
      GOTO 740
  720 JB = JLSUM
C
C         Outer loop over left operand
C
  730 I = JLSUM
  735 I = I + 1
      IF (I.GT.LLSUM) GOTO 770
      IBP = LISTBP(I)
      IF(IBP.LT.0) GOTO 1005
C
C         Inner loop over right operand
C
  740 J = JRSUM
  745 J = J + 1
      IF (J.GT.IJB) GOTO 760
      JBP = LISTBP(J)
      IF(JBP.LT.0) GOTO 1005
      IORJ = IOR(IBP, JBP)
      IF ((JOP.EQ.PMINUS).OR.(JOP.EQ.PMINSL).OR.(JOP.EQ.PMINST))GOTO 750
C
C         ** . / and * only
C
      IF ((LCHECK.AND.(IORJ.EQ.IBP)).OR.
     1    (RCHECK.AND.(IORJ.EQ.JBP)).OR.
     2    (ICNT(IORJ).GT.IORD)) GOTO 745
      JB = JB + 1
      IF (JB.GE.JP) GOTO 1002
      LISTBP(JB) = IORJ
      GOTO 745
C
C         - -/ and -* only
C
  750 IF (IORJ.NE.IBP) GOTO 745
      K = PMINUS
      IF (IBP.EQ.JBP) K = PMINSL
      IF (JOP - K) 735, 745, 735
C
C         End of inner loop
C         Copy integer not deleted in - -/ or -* operation
C
  760 IF ((JOP.NE.PMINUS).AND.(JOP.NE.PMINSL).AND.(JOP.NE.PMINST))
     1  GOTO 735
      JB = JB + 1
      LISTBP(JB) = IBP
      GOTO 735
C
C         End of outer loop
C         For .  /  and * reorder sum and copy omitting repeats
C
  770 IF ((JOP.EQ.PMINUS).OR.(JOP.EQ.PMINSL).OR.(JOP.EQ.PMINST))GOTO 790
      IF (JOP.NE.PSTAR2) GOTO 775
      NPOWER = NPOWER - 1
      IF(NPOWER.EQ.1) GO TO 775
      J = IJB
      GO TO 776
  775 IF (JOP.NE.PDOT) IJB = JLSUM
      J = JLSUM
  776 CALL BISORT(LISTBP, IJB, JB, IORDER)
      IBP = -1
      DO 779 I = IJB + 1, JB
      IF (LISTBP(I).EQ.IBP) GOTO 779
      IBP = LISTBP(I)
      J = J + 1
      LISTBP(J) = IBP
  779 CONTINUE
      JB = J
      IF (JOP.NE.PSTAR2) GO TO 790
      IF (NPOWER.EQ.1) GO TO 790
C         For ** form next power
      JRSUM = IJB
      IJB = JB
      GO TO 730
C
C         + merge right with left sum
C
  780 JB = JRSUM
      DO 784 J = JRSUM + 1, IJB
      JBP = LISTBP(J)
      DO 782 I = JLSUM + 1, LLSUM
      IF(JBP.EQ.LISTBP(I)) GO TO 784
  782 CONTINUE
      JB = JB + 1
      LISTBP(JB) = JBP
  784 CONTINUE
C
C         Set JRSUM, JRFAC
C
  790 JRSUM = JLSUM
      JRFAC = IOR(JLFAC, JRFAC)
      IF ((JOP.EQ.PMINUS).OR.(JOP.EQ.PMINSL).OR.(JOP.EQ.PMINST))
     1  JRFAC = JLFAC
      GOTO 200
C
C------------------ RETURN AND ERROR SETTING ---------------------------
C
 1005 IFAULT = IFAULT + 1
 1004 IFAULT = IFAULT + 1
 1003 IFAULT = IFAULT + 1
 1002 IFAULT = IFAULT + 1
 1001 IFAULT = IFAULT + 1
      LENBP = JC
  999 RETURN
      END
