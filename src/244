{
         Algorithm AS244  Appl. Statist. (1989) Vol. 38, No. 1  }

{ The procedure HiModel can be used to
     (1) Check decomposability for a hierarchical loglinear model;
     (2) Obtain the sets of variables onto which the model is collapsible;
     (3) Find the factor formula of MLEs for the model;
     (4) Find the reduced formula of the likelihood ratio test statistic. }


{ **********************************************************************

                          DATA STRUCTURES

  ********************************************************************** }

CONST  MaxVar=50;
       MaxGenerator=50;
       MaxSubGeClass=20;
       MaxInt=51;

TYPE   IntSet   = SET OF 1..MaxVar;
       SetArr   = ARRAY[1..MaxGenerator] OF IntSet;
       SetArr2  = ARRAY[1..MaxSubGeClass,1..MaxGenerator] OF IntSet;
       IntArr   = ARRAY [1..MaxSubGeClass] OF 1..MaxGenerator;
       IntGC    = 0..MaxGenerator;
       IntSubGC = 0..MaxSubGeClass;
       NonnegInt= 0..MaxInt;


{ **********************************************************************

                          PRINCIPAL PROCEDURE

  ********************************************************************** }

PROCEDURE HiModel(NumModels:                                       NonnegInt;
           M1SizeGeClass, M2SizeGeClass:                           IntGC;
           M1GeClass, M2GeClass:                                   SetArr;
       VAR M1NumNumer, M2NumNumer, M1NumDenom, M2NumDenom,
           LRNumNumer, LRNumDenom:                                 IntGC;
       VAR M1NumSubGC, M2NumSubGC, LRNumGCNumer, LRNumGCDenom:     IntSubGC;
       VAR Fault:                                                  NonnegInt;
       VAR M1Numerator, M2Numerator, M1Denominator,
           M2Denominator, LRNumerator, LRDenominator:              SetArr;
       VAR M1SubGeClass, M2SubGeClass, LRGCNumer, LRGCDenom:       SetArr2;
       VAR M1SizeSubGC, M2SizeSubGC, LRSizeNumerGC, LRSizeDenomGC: IntArr;
       VAR VariableSet:                                            IntSet);
 
VAR i,j,k,m,n: NonnegInt; VarSetModel: IntSet;
 
PROCEDURE Error(VAR Size: IntGC; MaxSize: IntGC; ErrorNumber: NonnegInt;
                    AnotherError: Boolean);
BEGIN  IF (Size>MaxSize) OR AnotherError THEN
         BEGIN  Fault:=ErrorNumber;  Size:=MaxSize
         END
END; { of Error }
 
{ Reduce the fraction: NumeratorSet/DenominatorSet. }
PROCEDURE Reduce(VAR NumNumerSet,NumDenomSet: IntGC;
                 VAR NumeratorSet, DenominatorSet: SetArr);
VAR i,j,k,m: IntGC;
BEGIN  k:=0;
 FOR i:=1 TO NumNumerSet DO
  BEGIN  j:=1;
   WHILE (NumeratorSet[i]<>DenominatorSet[j]) AND
         (j<=NumDenomSet) DO j:=j+1;
   IF j>NumDenomSet THEN
    BEGIN  k:=k+1;  NumeratorSet[k]:=NumeratorSet[i]
    END
   ELSE
    BEGIN
     FOR m:=j+1 TO NumDenomSet DO
      DenominatorSet[m-1]:=DenominatorSet[m];
     NumDenomSet:=NumDenomSet-1
    END
  END;
 NumNumerSet:=k
END; { of Reduce }
 
{ Check collapsibility and decomposability. }
PROCEDURE CollDeco (SizeGeClass: IntGC; GeClass: SetArr; VAR NumNumer,
  NumDenom: IntGC; VAR NumSubGC: IntSubGC; VAR Numerator,Denominator:
  SetArr; VAR SubGeClass: SetArr2; VAR SizeSubGC: IntArr);
VAR i: IntGC;
 
{ Check whether the generating class is decomposable or not. }
PROCEDURE Deco(VAR SizeGeClass: IntGC; VAR GeClass: SetArr);
VAR i,j: IntGC;  OneClass: IntSet;
    GeClassUnchange,GeneratorUnchange,NoSubGenerator: Boolean;
BEGIN
 REPEAT  GeClassUnchange:=true;
{ Delete variables which exist in only a single generator. }
  REPEAT  GeneratorUnchange:=true;
   FOR i:=1 TO SizeGeClass DO
    BEGIN  OneClass:=GeClass[i];
     FOR j:=1 TO SizeGeClass DO
      IF i<>j THEN OneClass:=OneClass-GeClass[j];
     IF OneClass<>[] THEN
      BEGIN
       NumNumer:=NumNumer+1;  Error(NumNumer,MaxGenerator,3,false);
       Numerator[NumNumer]:=GeClass[i];
       IF GeClass[i]-OneClass<>[] THEN
        BEGIN  NumDenom:=NumDenom+1;
         Error(NumDenom,MaxGenerator,4,false);
         Denominator[NumDenom]:=GeClass[i]-OneClass
        END;
        GeClass[i]:=GeClass[i]-OneClass;  GeneratorUnchange:=false
      END
    END
  UNTIL GeneratorUnchange;
{ Delete generators which are contained in others. }
  REPEAT  NoSubGenerator:=true;
   FOR i:=1 TO SizeGeClass DO
    FOR j:=1 TO SizeGeClass DO
     IF (i<>j) AND (GeClass[i]<=GeClass[j]) AND (GeClass[i]<>[]) THEN
      BEGIN  GeClass[i]:=[];
       NoSubGenerator:=false; GeClassUnchange:=false
      END;
{ Delete the null generators. }
   j:=0;
   FOR i:=1 TO SizeGeClass DO
    IF GeClass[i]<>[] THEN
     BEGIN j:=j+1;  GeClass[j]:=GeClass[i]
     END;
    SizeGeClass:=j
  UNTIL NoSubGenerator OR (SizeGeClass=0)
 UNTIL GeClassUnchange OR (SizeGeClass=0)
END; { of Deco }
 
{ Decompose generating class. }
PROCEDURE DecompGC(VAR SizeGeClass: IntGC; VAR GeClass: SetArr);
VAR i,j,k,InAGeneratorsNum,InBGeneratorsNum: IntGC; GeClassA: SetArr;
    BoundaryOfAB,VarInA,VarInB: IntSet;
    Researched: ARRAY [1..MaxGenerator] OF Boolean;
    DecomposedOrCannotDecompose,NoVarAddedInA: Boolean;
BEGIN
 Deco(SizeGeClass,GeClass);
 IF SizeGeClass<>0 THEN
  BEGIN  i:=0;  DecomposedOrCannotDecompose:=false;
{ Find variables in A. }
   REPEAT
    IF i=0 THEN BoundaryOfAB:=[] ELSE BoundaryOfAB:=GeClass[i];
    InAGeneratorsNum:=1;
    IF i<>1 THEN k:=1 ELSE k:=2;
    GeClassA[1]:=GeClass[k]; Researched[k]:=true; VarInA:=GeClassA[1];
    FOR j:=1 TO SizeGeClass DO IF j<>k THEN Researched[j]:=false;
    REPEAT  NoVarAddedInA:=true;
     FOR j:=1 TO SizeGeClass DO
      IF (i<>j) AND NOT Researched[j] AND
         (VarInA*GeClass[j]-BoundaryOfAB<>[]) THEN
       BEGIN NoVarAddedInA:=false; InAGeneratorsNum:=InAGeneratorsNum+1;
        GeClassA[InAGeneratorsNum]:=GeClass[j];  Researched[j]:=true;
        VarInA:=VarInA+GeClass[j]
       END
    UNTIL NoVarAddedInA;
    IF i<>0 THEN InBGeneratorsNum:=SizeGeClass-InAGeneratorsNum-1
            ELSE InBGeneratorsNum:=SizeGeClass-InAGeneratorsNum;
    VarInB:=[];
    IF InBGeneratorsNum<>0 THEN
{ The generating class can be decomposed by GeClass[i]. }
     BEGIN  j:=0;  DecomposedOrCannotDecompose:=true;
{ Find variables in B. }
      FOR k:=1 TO SizeGeClass DO
       IF (k<>i) AND NOT Researched[k] THEN
        BEGIN  j:=j+1; VarInB:=VarInB+GeClass[k]; GeClass[j]:=GeClass[k]
        END;
      IF i<>0 THEN
       BEGIN  InAGeneratorsNum:=InAGeneratorsNum+1;
        GeClassA[InAGeneratorsNum]:=BoundaryOfAB;
        InBGeneratorsNum:=InBGeneratorsNum+1;
        GeClass[InBGeneratorsNum]:=BoundaryOfAB*VarInB;
        NumDenom:=NumDenom+1;  Error(NumDenom,MaxGenerator,4,false);
        Denominator[NumDenom]:=GeClass[InBGeneratorsNum]
       END;
{ Decompose recursively A and B. }
      DecompGC(InAGeneratorsNum,GeClassA);
      DecompGC(InBGeneratorsNum,GeClass)
     END
    ELSE
{ The generating class cannot be decomposed by GeClass[i]. }
    IF i=SizeGeClass THEN
{ The generating class cannot be decomposed by any generator. }
     BEGIN  DecomposedOrCannotDecompose:=true; NumSubGC:=NumSubGC+1;
      IF NumSubGC>MaxSubGeClass THEN
       BEGIN  Fault:=5;  NumSubGC:=MaxSubGeClass
       END;
      SizeSubGC[NumSubGC]:=SizeGeClass;
      FOR k:=1 TO SizeGeClass DO SubGeClass[NumSubGC,k]:=GeClass[k]
     END
    ELSE i:=i+1
   UNTIL DecomposedOrCannotDecompose
  END
END; { of DecompGC }

BEGIN { of CollDeco }
 Error(SizeGeClass,MaxGenerator,1,SizeGeClass<1);  VarSetModel:=[];
 IF Fault=0 THEN
  FOR i:=1 TO SizeGeClass DO
   BEGIN
    IF NOT([1..MaxVar]>=GeClass[i]) OR (GeClass[i]=[]) THEN Fault:=2;
    VarSetModel:=VarSetModel+GeClass[i]
   END;
 NumNumer:=0;  NumDenom:=0;  NumSubGC:=0;
 IF Fault=0 THEN
  BEGIN
   DecompGC(SizeGeClass,GeClass);
   Reduce(NumNumer,NumDenom,Numerator,Denominator)
  END
END; { of CollDeco }
 
{ Check whether the sub-generating classes in models 1 and 2 are same. }
FUNCTION ClassesNotSame: Boolean;
VAR k: IntGC; m: NonnegInt; NotSame: Boolean;
BEGIN  k:=0;  NotSame:=M1SizeSubGC[i]<>LRSizeDenomGC[j];
 WHILE NOT NotSame AND (k<M1SizeSubGC[i]) DO
  BEGIN  k:=k+1;  m:=1;
   WHILE (m<=LRSizeDenomGC[j]) AND
         (M1SubGeClass[i,k]<>LRGCDenom[j,m]) DO m:=m+1;
   NotSame:=NotSame OR (m>LRSizeDenomGC[j])
  END;
 ClassesNotSame:=NotSame
END; { of ClassesNotSame }
 
BEGIN { of HiModel }
 M1NumNumer:=0; M1NumDenom:=0; M1NumSubGC:=0; VariableSet:=[];
 M2NumNumer:=0; M2NumDenom:=0; M2NumSubGC:=0; Fault:=0;
 LRNumNumer:=0; LRNumDenom:=0; LRNumGCNumer:=0; LRNumGCDenom:=0;
 IF (NumModels<1) OR (NumModels>2) THEN Fault:=6
 ELSE
  BEGIN
   CollDeco(M1SizeGeClass,M1GeClass,M1NumNumer,M1NumDenom,M1NumSubGC,
            M1Numerator,M1Denominator,M1SubGeClass,M1SizeSubGC);
   VariableSet:=VarSetModel;
   IF NumModels>1 THEN
    BEGIN
     CollDeco(M2SizeGeClass,M2GeClass,M2NumNumer,M2NumDenom,M2NumSubGC,
              M2Numerator,M2Denominator,M2SubGeClass,M2SizeSubGC);
     IF VariableSet<>VarSetModel THEN Fault:=9;
{ Find the likelihood ratio: M1/M2. }
     LRNumNumer:=M1NumNumer+M2NumDenom;
     Error(LRNumNumer,MaxGenerator,7,false);
     FOR i:=1 TO LRNumNumer DO
      IF i<=M1NumNumer THEN LRNumerator[i]:=M1Numerator[i]
         ELSE LRNumerator[i]:=M2Denominator[i-M1NumNumer];
     LRNumDenom:=M2NumNumer+M1NumDenom;
     Error(LRNumDenom,MaxGenerator,8,false);
     FOR i:=1 TO LRNumDenom DO
      IF i<=M1NumDenom THEN LRDenominator[i]:=M1Denominator[i]
         ELSE LRDenominator[i]:=M2Numerator[i-M1NumDenom];
     Reduce(LRNumNumer,LRNumDenom,LRNumerator,LRDenominator);
{ Reduce sub-generating classes of M1 and M2. }
     LRNumGCDenom:=M2NumSubGC;
     FOR m:=1 TO LRNumGCDenom DO
      BEGIN  LRSizeDenomGC[m]:=M2SizeSubGC[m];
       FOR n:=1 TO LRSizeDenomGC[m] DO LRGCDenom[m,n]:=M2SubGeClass[m,n]
      END;
     k:=0;
     FOR i:=1 TO M1NumSubGC DO
      BEGIN  j:=1;
       WHILE (j<=LRNumGCDenom) AND ClassesNotSame DO j:=j+1;
       IF j>LRNumGCDenom THEN
        BEGIN  k:=k+1;  LRSizeNumerGC[k]:=M1SizeSubGC[i];
         FOR m:=1 TO LRSizeNumerGC[k] DO
          LRGCNumer[k,m]:=M1SubGeClass[i,m]
        END
       ELSE
        BEGIN
         FOR m:=j+1 TO LRNumGCDenom DO
          BEGIN  LRSizeDenomGC[m-1]:=LRSizeDenomGC[m];
           FOR n:=1 TO LRSizeDenomGC[m-1] DO
            LRGCDenom[m-1,n]:=LRGCDenom[m,n]
          END;
         LRNumGCDenom:=LRNumGCDenom-1
        END
      END;
     LRNumGCNumer:=k
    END
  END
END; { of HiModel }
