      SUBROUTINE HYBRID(NPARAM,PARAM,N,T0,RHO,NT,BL,BU,NP,POINTS,LOW,
     +                  IFAULT)
C     .. Scalar Arguments ..
      DOUBLE PRECISION RHO,T0
      INTEGER IFAULT,N,NP,NPARAM,NT
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION BL(NPARAM),BU(NPARAM),LOW(NPARAM+1),
     +                 PARAM(NPARAM),POINTS(NPARAM,N)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION BST,DIFF,NEWOBJ,OBJF,ONE,P,T,U,ZERO
      INTEGER I,I1,I2
C     ..
C     .. External Functions ..
      DOUBLE PRECISION RANDOM
      EXTERNAL RANDOM
C     ..
C     .. External Subroutines ..
      EXTERNAL FUNCT1,TRADNL
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DEXP
C
C      Define the constants ZERO and ONE and set
C      initial values for BST and IFAULT
C
C     .. Data statements ..
      DATA ZERO,ONE,BST/0.0D0,1.0D0,1.0D+10/
C     ..
      IFAULT = 0
C
C      ******************************************
C      *** Start of Annealing component to    ***
C      *** generate a set of starting points. ***
C      ******************************************
C
C      Initialise the temperature
C
      T = T0
C
C      Choose an initial point
C
      DO 10 I = 1,NPARAM
          PARAM(I) = RANDOM(BL(I),BU(I))
   10 CONTINUE
C
C      Calculate the function value for the initial point
C
      CALL FUNCT1(NPARAM,PARAM,OBJF)
C
C      Begin the loop over I1 until Nt Temp. reductions.
C
      DO 90 I1 = 1,NT
C
C      Set the counter for the No of points accepted to one
C
          NP = 1
C
C      Begin the loop over I2 until N points have been considered.
C
          DO 80 I2 = 1,N
C
C      Record the present point
C
              DO 20 I = 1,NPARAM
                  POINTS(I,NP+1) = PARAM(I)
   20         CONTINUE
C
C      Select a new point
C
              DO 30 I = 1,NPARAM
                  PARAM(I) = RANDOM(BL(I),BU(I))
   30         CONTINUE
C
C      Calculate the function value of this new point
C
              CALL FUNCT1(NPARAM,PARAM,NEWOBJ)
C
C      Calculate the difference
C
              DIFF = (NEWOBJ-OBJF)
C
C      If the new point is better than the old one then increase the
C      Np counter, and record it as a possible starting point. Move
C      to this new point.
C
              IF (DIFF.LT.ZERO) THEN
                  NP = NP + 1
                  DO 40 I = 1,NPARAM
                      POINTS(I,NP) = PARAM(I)
   40             CONTINUE
                  OBJF = NEWOBJ
                  IF (NEWOBJ.LT.BST) THEN
                      BST = NEWOBJ
                      DO 50 I = 1,NPARAM
                          POINTS(I,1) = PARAM(I)
   50                 CONTINUE
                  END IF
C
C      Otherwise use the Metropolis Criterion
C
              ELSE
C
C      P is the comparison statistic. U is a standard Uniform variate.
C
                  P = DEXP(-DIFF/T)
                  U = RANDOM(ZERO,ONE)
C
C      If U < P then accept this point. Increase the Np counter,
C      Record this point as a possible starting point, and move
C      to this new point.
C
                  IF (U.LT.P) THEN
                      NP = NP + 1
                      OBJF = NEWOBJ
                      DO 60 I = 1,NPARAM
                          POINTS(I,NP) = PARAM(I)
   60                 CONTINUE
C
C      If U >= P then reject this point and return to the old point.
C
                  ELSE
                      DO 70 I = 1,NPARAM
                          PARAM(I) = POINTS(I,NP+1)
   70                 CONTINUE
                  END IF
C
C      Both alternatives for the comparison of the new and old point
C      Have been considered.
C
              END IF
C
C      Continue with the loop until N new points have been considered
C
   80     CONTINUE
C
C      Set the error indicator to 1 if the algorithm appears to have
C      converged.
C
          IF (NP.EQ.1) THEN
              IFAULT = I1
          END IF
C
C      Once N points have been considered, lower the Temperature.
C
          T = T*RHO
C
C      Continue with the loop until the Temp. has been reduced Nt times
C
   90 CONTINUE
C
C       *****************************************
C       *** End of Annealing component, NP    ***
C       *** starting points are now stored in ***
C       *** the array POINTS.                 ***
C       *****************************************
C
C       Set an initially high value of LOW(NPARAM+1)
C
      LOW(NPARAM+1) = 1.0D+10
C
C       *******************************************
C       *** Start of the traditional component, ***
C       *** which loops over the NP starting    ***
C       *** points calling the traditional      ***
C       *** minimisation routine for each.      ***
C       *******************************************
C
C       Begin the loop over the Np starting points.
C
      DO 130 I1 = 1,NP
C
C       Set the PARAM array to be the values of the I1st start point
C
          DO 110 I2 = 1,NPARAM
              PARAM(I2) = POINTS(I2,I1)
  110     CONTINUE
C
C       Call traditional minimisation routine
C
          CALL TRADNL(NPARAM,PARAM,BL,BU,OBJF)
C
C       Check to see if this starting point gives the
C       lowest function value to date
C
          IF (OBJF.LT.LOW(NPARAM+1)) THEN
              DO 120 I = 1,NPARAM
                  LOW(I) = PARAM(I)
  120         CONTINUE
              LOW(NPARAM+1) = OBJF
          END IF
C
C       End the Np loop
C
  130 CONTINUE
C
C       **************************************
C       *** End of the second, traditional ***
C       *** component.                     ***
C       **************************************
C
C       Return to calling routine
C
      RETURN
C
C       End of subroutine HYBRID
C
      END

