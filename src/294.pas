PROCEDURE Error(VAR Size: IntGC; MaxSize: IntGC; ErrorNumber: NonnegInt;
                AnotherError: Boolean; VAR Fault: NonnegInt);

{ALGORITHM AS 294.1 APPL.STATIST. (1994) VOL.43, N0.3}

BEGIN  IF (Size>MaxSize) OR AnotherError THEN
        BEGIN  Fault:=ErrorNumber;  Size:=MaxSize
        END
END; {of Error}

PROCEDURE Reduce(VAR NumNumerSet,NumDenomSet:              IntGC;
    VAR NumeratorSet,NumerSetOC,DenominatorSet,DenomSetOC: SetArr;
        AnObsConfig:                                       IntSet);

{ALGORITHM AS 294.2 APPL.STATIST. (1994) VOL.43, N0.3}

{Reduce the fraction: NumeratorSet/DenominatorSet}

VAR i,j,k,m: IntGC;  Reducible: Boolean;
BEGIN  k:=0;
 FOR i:=1 TO NumNumerSet DO
  BEGIN
   IF AnObsConfig>=NumeratorSet[i] THEN
    BEGIN  j:=1;
     WHILE (j<=NumDenomSet) AND (NOT(AnObsConfig>=DenominatorSet[j]) OR
           (NumerSetOC[i]<>DenomSetOC[j]) OR
           (NumeratorSet[i]<>DenominatorSet[j])) DO j:=j+1;
     Reducible:= j <= NumDenomSet;
     IF Reducible THEN
      BEGIN
       FOR m:=j+1 TO NumDenomSet DO
        BEGIN DenominatorSet[m-1]:=DenominatorSet[m];
         DenomSetOC[m-1]:=DenomSetOC[m]
        END;
       NumDenomSet:=NumDenomSet-1
      END
    END;
   IF NOT Reducible OR NOT(AnObsConfig>=NumeratorSet[i]) THEN
    BEGIN k:=k+1; NumeratorSet[k]:=NumeratorSet[i];
     NumerSetOC[k]:=NumerSetOC[i]
    END
  END; {of FOR i}
 NumNumerSet:=k
END; {of Reduce}

PROCEDURE MissMLE (SizeGeClass, NumObsConfig:               IntGC;
      GeClass, ObsConfig:                                   SetArr;
  VAR NumNumer, NumDenom,NumSubGC:                          IntGC;
  VAR Numerator,Denominator,NumerOC,DenomOC,SGClassOC:      SetArr;
  VAR SubGeClass:                                           SetArr2;
  VAR SizeSubGC:                                            IntArr;
  VAR VarSetModel:                                          IntSet;
  VAR Fault:                                                NonnegInt);

{ALGORITHM AS 294.3 APPL.STATIST. (1994) VOL.43, N0.3}

{Find the factor formula of MLEs}

VAR i,j,k,m,OrigNumNumer: IntGC;  Incomplete: Boolean;
    Variables,NumerConfig,DenomConfig,OCInModel,Block,BlockOC: IntSet;

PROCEDURE Nest(Block:IntSet; VAR BlockOC: IntSet);

{Factorizing by nested pattern of observations}

VAR j,m: NonnegInt; OCj,OCUnion,Union,Intersect: IntSet;
    Unchange,Separa: Boolean;
BEGIN
 IF BlockOC<>[] THEN
{Find the union of inside observing configurations, OCUnion}
 BEGIN  Intersect:=Block; Union:=[]; OCUnion:=[];
  REPEAT  Unchange:=True;
   FOR j:=1 TO NumObsConfig DO
    IF j IN BlockOC THEN
     BEGIN OCj:=ObsConfig[j]*Block;
      IF Intersect>=OCj THEN
       BEGIN  Intersect:=OCj; OCUnion:=OCj; Union:=[j]
       END
      ELSE
       BEGIN Intersect:=Intersect*OCj;
        IF NOT(j IN Union) AND (NOT(OCj>=OCUnion) OR (OCj=OCUnion)) THEN
         BEGIN  OCUnion:=OCUnion+OCj; Union:=Union+[j]; Unchange:=false
         END
       END
     END
  UNTIL  Unchange;
  IF NOT (OCUnion>=Block) THEN
   BEGIN
{Check whether OCUnion is contained in a generator}
    IF i<=OrigNumNumer THEN Separa:=True
    ELSE
     BEGIN  m:=i-OrigNumNumer; j:=0;
      REPEAT j:=j+1; Separa:=OCUnion<=SubGeClass[m,j]
      UNTIL Separa OR (j>=SizeSubGC[m])
     END;
    IF Separa THEN
{Factorize the MLEs by nested pattern}
     BEGIN  NumNumer:=NumNumer+1;
      Error(NumNumer,MaxGenerator,3,false,Fault);
      Numerator[NumNumer]:=OCUnion; NumerOC[NumNUmer]:=BlockOC;
      NumDenom:=NumDenom+1; Error(NumDenom,MaxGenerator,4,false,Fault);
      Denominator[NumDenom]:=OCUnion; DenomOC[NumDenom]:=BlockOC-Union;
      BlockOC:=BlockOC-Union; Nest(Block,BlockOC)
     END
   END {of IF NOT (OCUnion>=Block) THEN}
 END {of IF BlockOC<>[] THEN}
END; {of Nest}

{Check Whether the block A is separable}
FUNCTION Separable(OtherConditionHolds:Boolean;
   VarsInA,Generator:IntSet; VAR Bound,OCInModel:IntSet): Boolean;
VAR k: NonnegInt; Unchange,Separa: Boolean;
BEGIN
 IF NOT OtherConditionHolds THEN Separable:=false
 ELSE
  IF Bound=[] THEN
{For the cases that the blocks A and B do not intersect:}
   BEGIN  Separable:=True; NumerConfig:=OCInModel;
    FOR k:=1 TO NumObsConfig DO
     IF k IN OCInModel THEN
      BEGIN  IF ObsConfig[k]*(Variables-VarsInA)=[]
             THEN OCInModel:=OCInModel-[k];
       IF ObsConfig[k]*VarsInA=[] THEN NumerConfig:=NumerConfig-[k]
      END
   END
  ELSE
{For the cases that the blocks A and B intersect:}
   BEGIN
{Adjust the Bound according to the observing configurations}
    REPEAT   Unchange:=true;
     FOR k:=1 TO NumObsConfig DO
      IF k IN OCInModel THEN
       IF (ObsConfig[k]*(VarsInA-Bound)<>[]) AND
          (ObsConfig[k]*(Variables-VarsInA)<>[]) AND
          NOT (ObsConfig[k]>=Bound) THEN
        BEGIN  Bound:=Bound+ObsConfig[k]*VarsInA; Unchange:=false
        END
    UNTIL  Unchange;
    IF NOT (Generator>=Bound) OR (VarsInA-Bound=[]) THEN
     Separable:=false
    ELSE
     BEGIN
{Check whether the block A is separable as a submodel}
      k:=0;  NumerConfig:=[];  Separa:=true;
      REPEAT  k:=k+1;
       IF k IN OCInModel THEN
        BEGIN  Separa:=ObsConfig[k]*(VarsInA-Bound)=[];
         IF NOT Separa THEN
          BEGIN  Separa:= Bound<=ObsConfig[k];
           IF Separa THEN NumerConfig:=NumerConfig+[k]
          END
        END
      UNTIL (k>=NumObsConfig) OR NOT Separa;
      IF Separa THEN DenomConfig:=NumerConfig
      ELSE
{Check whether the block A is separable as a main model}
       BEGIN  k:=0;  Separa:=true;
        REPEAT  k:=k+1;
         IF k IN OCInModel THEN
          Separa:= (VarsInA>=ObsConfig[k]*Variables)
                OR (Bound<=ObsConfig[k])
        UNTIL (k>=NumObsConfig) OR NOT Separa;
        IF Separa THEN
{Find sets of observing configurations of blocks A, B and denominator}
         BEGIN NumerConfig:=OCInModel; DenomConfig:=OCInModel;
          FOR k:=1 TO NumObsConfig DO
           IF (k IN OCInModel) AND(VarsInA>=ObsConfig[k]*Variables) THEN
            DenomConfig:=DenomConfig-[k];
          OCInModel:=DenomConfig
         END
       END; {of IF Separa ... ELSE}
      Separable:=Separa
     END {of IF NOT (Generator>=Bound)...ELSE}
   END {of IF Bound=[]...ELSE}
END; {of Separable}

{Check whether the generating class is decomposable or not}
PROCEDURE Deco(VAR SizeGeClass: IntGC; VAR GeClass: SetArr;
    VAR OCInModel:IntSet);
VAR i,j: IntGC;  Bound,OneClass: IntSet;
    GeClassUnchange,GeneratorUnchange,NoSubGenerator: Boolean;
BEGIN {of Deco}
 REPEAT  GeClassUnchange:=true;
{Delete variables which exist in only a single generator}
  REPEAT  GeneratorUnchange:=true;
   FOR i:=1 TO SizeGeClass DO
    BEGIN  OneClass:=GeClass[i];
     FOR j:=1 TO SizeGeClass DO
      IF i<>j THEN OneClass:=OneClass-GeClass[j];
     Bound:=GeClass[i]-OneClass;
     IF Separable(OneClass<>[],GeClass[i],GeClass[i],Bound,OCInModel)
      THEN
      BEGIN Variables:=Variables-GeClass[i]+Bound;
       NumNumer:=NumNumer+1; Error(NumNumer,MaxGenerator,3,false,Fault);
       Numerator[NumNumer]:=GeClass[i]; NumerOC[NumNumer]:=NumerConfig;
       IF Bound<>[] THEN
        BEGIN NumDenom:=NumDenom+1;
         Error(NumDenom,MaxGenerator,4,false,Fault);
         Denominator[NumDenom]:=Bound; DenomOC[NumDenom]:=DenomConfig
        END;
       GeClass[i]:=Bound;  GeneratorUnchange:=false
      END
    END
  UNTIL GeneratorUnchange;
{Delete generators which are contained in others}
  REPEAT  NoSubGenerator:=true;
   FOR i:=1 TO SizeGeClass DO
    FOR j:=1 TO SizeGeClass DO
     IF (i<>j) AND (GeClass[i]<=GeClass[j]) AND (GeClass[i]<>[]) THEN
      BEGIN GeClass[i]:=[];NoSubGenerator:=false;GeClassUnchange:=false
      END;
{Delete the null generators}
   j:=0;
   FOR i:=1 TO SizeGeClass DO
    IF GeClass[i]<>[] THEN
     BEGIN j:=j+1;  GeClass[j]:=GeClass[i]
     END;
   SizeGeClass:=j
  UNTIL NoSubGenerator OR (SizeGeClass=0)
 UNTIL GeClassUnchange OR (SizeGeClass=0)
END; {of Deco}

{Decompose generating class}
PROCEDURE DecompGC(VAR SizeGeClass: IntGC; VAR GeClass: SetArr;
    OCInModel: IntSet);
VAR i,j,k,m,InAGeneratorsNum,InBGeneratorsNum: IntGC; GeClassA: SetArr;
    BoundaryOfAB,Bound,VarInA,VarInB,T: IntSet;
    Researched: ARRAY [1..MaxGenerator] OF Boolean;
    DecomposedOrCannotDecompose,NoVarAddedInA: Boolean;
BEGIN  Variables:=[];
 FOR i:=1 TO SizeGeClass DO Variables:=Variables+GeClass[i];
 Deco(SizeGeClass,GeClass,OCInModel);
 IF SizeGeClass<>0 THEN
  BEGIN  i:=0;  DecomposedOrCannotDecompose:=false;
   REPEAT
{Set [] or GeClass[i] as the separator, BoundaryOfAB}
    IF i=0 THEN BoundaryOfAB:=[] ELSE BoundaryOfAB:=GeClass[i];  m:=0;
    REPEAT  m:=m+1;
     IF (m<>i) AND NOT ((i=0) AND (m<>1))  THEN
      BEGIN  InAGeneratorsNum:=1; GeClassA[1]:=GeClass[m];
       Researched[m]:=true; VarInA:=GeClassA[1];
       FOR j:=1 TO SizeGeClass DO IF j<>m THEN Researched[j]:=false;
{Find the set of variables in A, VarInA}
       REPEAT  NoVarAddedInA:=true;
        FOR j:=1 TO SizeGeClass DO
         IF (i<>j) AND NOT Researched[j] AND
            (VarInA*GeClass[j]-BoundaryOfAB<>[]) THEN
          BEGIN NoVarAddedInA:=false;  VarInA:=VarInA+GeClass[j];
           Researched[j]:=true; InAGeneratorsNum:=InAGeneratorsNum+1;
           GeClassA[InAGeneratorsNum]:=GeClass[j]
          END
       UNTIL NoVarAddedInA;
       IF i<>0 THEN InBGeneratorsNum:=SizeGeClass-InAGeneratorsNum-1
       ELSE  InBGeneratorsNum:=SizeGeClass-InAGeneratorsNum;
       VarInA:=VarInA+BoundaryOfAB;  VarInB:=[];
{Find the set of variables in B, VarInB}
       FOR k:=1 TO SizeGeClass DO
        IF (k<>i) AND NOT Researched[k] THEN VarInB:=VarInB+GeClass[k];
       Bound:=VarInA*VarInB;
       IF Separable(InBGeneratorsNum<>0,VarInA,BoundaryOfAB,
                    Bound,OCInModel) THEN
{The generating class can be decomposed by GeClass[i]}
        BEGIN  j:=0;  DecomposedOrCannotDecompose:=true;
         FOR k:=1 TO SizeGeClass DO
          IF (k<>i) AND NOT Researched[k] THEN
           BEGIN  j:=j+1; GeClass[j]:=GeClass[k]
           END;
         IF i<>0 THEN
          BEGIN  InAGeneratorsNum:=InAGeneratorsNum+1;
           GeClassA[InAGeneratorsNum]:=BoundaryOfAB;
           InBGeneratorsNum:=InBGeneratorsNum+1;
           GeClass[InBGeneratorsNum]:=Bound;
           NumDenom:=NumDenom+1;
           Error(NumDenom,MaxGenerator,4,false,Fault);
           Denominator[NumDenom]:=GeClass[InBGeneratorsNum];
           DenomOC[NumDenom]:=DenomConfig
          END;
{Decompose recursively A and B}
         DecompGC(InAGeneratorsNum,GeClassA,NumerConfig);
         DecompGC(InBGeneratorsNum,GeClass,OCInModel)
        END {of IF Separable(...) THEN}
      END {of IF (m<>i) AND NOT ((i=0) AND (m<>1)) THEN}
    UNTIL (m>=SizeGeClass) OR DecomposedOrCannotDecompose;
    IF NOT DecomposedOrCannotDecompose THEN
     BEGIN
{The generating class cannot be decomposed by GeClass[i]}
      IF i=SizeGeClass THEN
{The generating class cannot be decomposed by any generator}
       BEGIN  DecomposedOrCannotDecompose:=true; NumSubGC:=NumSubGC+1;
        Error(NumSubGC,MaxGenerator,5,false,Fault);
        SGClassOC[NumSubGC]:=OCInModel;SizeSubGC[NumSubGC]:=SizeGeClass;
        FOR k:=1 TO SizeGeClass DO SubGeClass[NumSubGC,k]:=GeClass[k]
       END
      ELSE i:=i+1
     END {of IF NOT Decomp... THEN}
   UNTIL DecomposedOrCannotDecompose
  END {of IF SizeGeClass<>0 THEN}
END; {of DecompGC}
BEGIN {of MissMLE}    Fault:=0;
 Error(SizeGeClass,MaxGenerator,1,SizeGeClass<1,Fault); VarSetModel:=[];
 IF (NumObsConfig<=0) OR (NumObsConfig>MaxGenerator) THEN Fault:=6;
 IF Fault=0 THEN
  FOR i:=1 TO SizeGeClass DO
   BEGIN
    IF NOT([1..MaxVar]>=GeClass[i]) OR (GeClass[i]=[]) THEN Fault:=2;
    VarSetModel:=VarSetModel+GeClass[i]
   END;
 IF Fault=0 THEN
  FOR i:=1 TO NumObsConfig DO
   IF NOT([1..MaxVar]>=ObsConfig[i])OR(ObsConfig[i]=[]) THEN Fault:=7;
 NumNumer:=0;  NumDenom:=0;  NumSubGC:=0;
 IF Fault=0 THEN
  BEGIN  OCInModel:=[1..NumObsConfig];
   DecompGC(SizeGeClass,GeClass,OCInModel);  OrigNumNumer:=NumNumer;
{Factorize MLEs according to nested pattern}
   FOR i:=1 TO OrigNumNumer+NumSubGC DO
    BEGIN  IF i<=OrigNumNumer THEN
            BEGIN Block:=Numerator[i]; Nest(Block,NumerOC[i])
            END
     ELSE
      BEGIN  m:=i-OrigNumNumer;  Block:=[];
       FOR j:=1 TO SizeSubGC[m] DO Block:=Block+SubGeClass[m,j];
       Nest(Block,SGClassOC[m])
      END
    END;
   Reduce(NumNumer,NumDenom,Numerator,NumerOC,Denominator,DenomOC,
          VarSetModel);
{Change the incomplete Numerator into SubGeClass}
   k:=0;
   FOR i:=1 TO NumNumer DO
    BEGIN  j:=1;
     Incomplete:=(j IN NumerOC[i]) AND NOT(ObsConfig[j]>=Numerator[i]);
     WHILE (j<=NumObsConfig) AND NOT Incomplete DO
      BEGIN  Incomplete:=(j IN NumerOC[i]) AND
                         NOT(ObsConfig[j]>=Numerator[i]);  j:=j+1
      END;
     IF Incomplete THEN
      BEGIN NumSubGC:=NumSubGC+1;
       Error(NumSubGC,MaxGenerator,5,false,Fault);
       SGClassOC[NumSubGC]:=NumerOC[i];  SizeSubGC[NumSubGC]:=1;
       SubGeClass[NumSubGC,1]:=Numerator[i]
      END
     ELSE
      BEGIN  k:=k+1;
       Numerator[k]:=Numerator[i]; NumerOC[k]:=NumerOC[i]
      END
    END; {of FOR i}
   NumNumer:=k
  END {of IF Fault=0 THEN}
END; {of MissMLE}

PROCEDURE MissHiMod(NumModels:                               NonnegInt;
      M1SizeGeClass,M2SizeGeClass,NumObsConfig:              IntGC;
      M1GeClass,M2GeClass,ObsConfig:                         SetArr;
  VAR LRNNumNumer,LRNNumDenom,LRNNumGeClass,
      LRDNumNumer,LRDNumDenom,LRDNumGeClass:                 IntArr;
  VAR LRNSizeGeClass,LRDSizeGeClass:                         IntArr2;
  VAR LRNNumer,LRNDenom,LRDNumer,LRDDenom,LRNOCNumer,LRNOCDenom,
      LRNOCGeClass,LRDOCNumer,LRDOCDenom,LRDOCGeClass:       SetArr2;
  VAR LRNGeClass,LRDGeClass:                                 SetArr3;
  VAR VariableSet:                                           IntSet;
  VAR Fault:                                                 NonnegInt);

{ALGORITHM AS 294.4 APPL.STATIST. (1994) VOL.43, N0.3}

{Find the reduced formula of LR statistic}

VAR i,j,k,m,n,p: NonnegInt;       VarSetModel,UnionN,Union:   IntSet;
    Reducible:   Boolean;         M1SubGeClass,M2SubGeClass:  SetArr2;
    M1NumNumer,M2NumNumer,M1NumDenom,M2NumDenom,
      M1NumSubGC,M2NumSubGC,NNumNumer,DNumNumer:              IntGC;
    M1SizeSubGC,M2SizeSubGC,NNumerInd,DNumerInd:              IntArr;
    M1Numerator,M2Numerator,M1Denominator,M2Denominator,
      M1NumerOC,M2NumerOC,M1DenomOC,M2DenomOC,M1SGClassOC,
      M2SGClassOC,NNumer,DNumer,NOCNumer,DOCNumer:            SetArr;
{Check whether the sub-generating classes in models 1 and 2 are same}
FUNCTION ClassesNotSame(h, i, j: NonnegInt): Boolean;
VAR k: IntGC;  m: NonnegInt;  NotSame: Boolean;
BEGIN  k:=0;  NotSame:=LRNSizeGeClass[h,i]<>LRDSizeGeClass[h,j];
 WHILE NOT NotSame AND (k<LRNSizeGeClass[h,i]) DO
  BEGIN  k:=k+1;  m:=1;
   WHILE (m<=LRDSizeGeClass[h,j]) AND
         (LRNGeClass[h,i,k]<>LRDGeClass[h,j,m]) DO m:=m+1;
   NotSame:=NotSame OR (m>LRDSizeGeClass[h,j])
  END;
 ClassesNotSame:=NotSame
END; {of ClassesNotSame}

{Find the union of generators}
PROCEDURE UnionSet(NumSet: IntGC;  Sets: SetArr;  VAR Union: IntSet);
VAR i: NonnegInt;
BEGIN  Union:=[];  FOR i:=1 TO NumSet DO  Union:=Union+Sets[i]
END;

{Set summed numerator and summed denominator of LR}
PROCEDURE Sum(ObsConfig: IntSet; LRNumNumer,LRNumGeClass: IntGC;
  VAR NumNumer: IntGC;  VAR LRNumer,LROCNumer,LROCGeClass,Numer,OCNumer:
   SetArr;  VAR LRSizeGeClass,NumerInd: IntArr; VAR LRGeClass: SetArr2);
VAR p:IntGC; j,k,m: NonnegInt;  Union,Unionj,UnionOther: IntSet;
BEGIN  UnionOther:=[];
 FOR k:=1 TO LRNumGeClass DO
  FOR m:=1 TO LRSizeGeClass[k] DO UnionOther:=UnionOther+LRGeClass[k,m];
 p:=0;
 FOR j:=1 TO LRNumNumer DO
  BEGIN Union:=UnionOther;
   FOR k:=1 TO LRNumNumer DO  IF k<>j THEN Union:=Union+LRNumer[k];
   IF LRNumer[j]*Union <= ObsConfig THEN
    BEGIN p:=p+1; Numer[p]:=LRNumer[j]*ObsConfig;
     OCNumer[p]:=LROCNumer[j]; NumerInd[p]:=j
    END
  END; {of FOR j}
 UnionOther:=[];
 FOR k:=1 TO LRNumNumer DO UnionOther:=UnionOther+LRNumer[k];
 FOR j:=1 TO LRNumGeClass DO
  BEGIN Union:=UnionOther; Unionj:=[];
   FOR k:=1 TO LRNumGeClass DO
    FOR m:=1 TO LRSizeGeClass[k] DO
     IF k=j THEN Unionj:=Unionj+LRGeClass[j,m]
     ELSE        Union :=Union +LRGeClass[k,m];
   IF Unionj*Union <= ObsConfig THEN
    BEGIN  m:=1;  Unionj:=Unionj*ObsConfig;
     WHILE (m<=LRSizeGeClass[j]) AND NOT(LRGeClass[j,m]>=Unionj) DO
      m:=m+1;
     IF m<=LRSizeGeClass[j] THEN
      BEGIN p:=p+1; Error(p,MaxGenerator,8,false,Fault);
       OCNumer[p]:=LROCGeClass[j]; Numer[p]:=Unionj;
       NumerInd[p]:=LRNumNumer+j
      END
    END
  END; {of FOR j}
 NumNumer:=p
END; {of Sum}

{Delete null numerator}
PROCEDURE Delete(VAR LRNumNumer,LRNumGeClass: IntGC;
  VAR LRNumer,LROCNumer,LROCGeClass: SetArr; VAR LRSizeGeClass: IntArr;
  VAR LRGeClass: SetArr2);
VAR j,k,m: NonnegInt;
BEGIN  k:=0;
 FOR j:=1 TO LRNumNumer DO  IF LRNumer[j] <> [] THEN
  BEGIN k:=k+1; LRNumer[k]:=LRNumer[j]; LROCNumer[k]:=LROCNumer[j]
  END;
 LRNumNumer:=k; k:=0;
 FOR j:=1 TO LRNumGeClass DO  IF LRSizeGeClass[j]<>0 THEN
  BEGIN k:=k+1; LRSizeGeClass[k]:=LRSizeGeClass[j];
   LROCGeClass[k]:=LROCGeClass[j];
   FOR m:=1 TO LRSizeGeClass[j] DO LRGeClass[k,m]:=LRGeClass[j,m]
  END;
 LRNumGeClass:=k
END; {of Delete}
BEGIN {of MissHiMod}
 M1NumNumer:=0; M1NumDenom:=0; M1NumSubGC:=0; VariableSet:=[];
 M2NumNumer:=0; M2NumDenom:=0; M2NumSubGC:=0; Fault:=0;
 IF NumModels <> 2 THEN Fault:=9
 ELSE
  BEGIN
   MissMLE(M1SizeGeClass,NumObsConfig,M1GeClass,ObsConfig,M1NumNumer,
     M1NumDenom,M1NumSubGC,M1Numerator,M1Denominator,M1NumerOC,
     M1DenomOC,M1SGClassOC,M1SubGeClass,M1SizeSubGC,VarSetModel,Fault);
   VariableSet:=VarSetModel;
   IF Fault=0 THEN
    BEGIN
     MissMLE(M2SizeGeClass,NumObsConfig,M2GeClass,ObsConfig,M2NumNumer,
      M2NumDenom,M2NumSubGC,M2Numerator,M2Denominator,M2NumerOC,
      M2DenomOC,M2SGClassOC,M2SubGeClass,M2SizeSubGC,VarSetModel,Fault);
     IF VariableSet<>VarSetModel THEN Fault:=10;
     IF Fault=0 THEN
{Find the likelihood ratio: M1/M2}
      FOR i:=1 TO NumObsConfig DO
       BEGIN  k:=0;
{Set the numerator of LR}
        FOR j:=1 TO M1NumNumer DO
         IF i IN M1NumerOC[j] THEN
          BEGIN  k:=k+1; LRNNumer[i,k]:=M1Numerator[j];
           LRNOCNumer[i,k]:=M1NumerOC[j]
          END;
        LRNNumNumer[i]:=k; k:=0;
        FOR j:=1 TO M1NumDenom DO
         IF i IN M1DenomOC[j] THEN
          BEGIN k:=k+1; LRNDenom[i,k]:=M1Denominator[j];
           LRNOCDenom[i,k]:=M1DenomOC[j]
          END;
        LRNNumDenom[i]:=k; k:=0;
        FOR j:=1 TO M1NumSubGC DO
         IF i IN M1SGClassOC[j] THEN
          BEGIN  k:=k+1; LRNSizeGeClass[i,k]:=M1SizeSubGC[j];
           LRNOCGeClass[i,k]:=M1SGClassOC[j];
           FOR m:=1 TO M1SizeSubGC[j] DO
            LRNGeClass[i,k,m]:=M1SubGeClass[j,m]
          END;
        LRNNumGeClass[i]:=k; k:=0;
{Set the denominator of LR}
        FOR j:=1 TO M2NumNumer DO
        IF i IN M2NumerOC[j] THEN
          BEGIN  k:=k+1; LRDNumer[i,k]:=M2Numerator[j];
           LRDOCNumer[i,k]:=M2NumerOC[j]
          END;
        LRDNumNumer[i]:=k; k:=0;
        FOR j:=1 TO M2NumDenom DO
         IF i IN M2DenomOC[j] THEN
          BEGIN k:=k+1; LRDDenom[i,k]:=M2Denominator[j];
           LRDOCDenom[i,k]:=M2DenomOC[j]
          END;
        LRDNumDenom[i]:=k; k:=0;
        FOR j:=1 TO M2NumSubGC DO
         IF i IN M2SGClassOC[j] THEN
          BEGIN  k:=k+1; LRDSizeGeClass[i,k]:=M2SizeSubGC[j];
           LRDOCGeClass[i,k]:=M2SGClassOC[j];
           FOR m:=1 TO M2SizeSubGC[j] DO
            LRDGeClass[i,k,m]:=M2SubGeClass[j,m]
          END;
        LRDNumGeClass[i]:=k;
{Reduce the numerator and denominator contained in ObsConfig[i]}
        Reduce(LRNNumNumer[i],LRDNumNumer[i],LRNNumer[i],LRNOCNumer[i],
               LRDNumer[i],LRDOCNumer[i],ObsConfig[i]);
        Reduce(LRNNumDenom[i],LRDNumDenom[i],LRNDenom[i],LRNOCDenom[i],
               LRDDenom[i],LRDOCDenom[i],ObsConfig[i]);
{Reduce sub-generating classes contained in ObsConfig[i]}
        k:=0;
        FOR m:=1 TO LRNNumGeClass[i] DO
         BEGIN  j:=1;
          UnionSet(LRNSizeGeClass[i,m],LRNGeClass[i,m],UnionN);
          IF ObsConfig[i] >= UnionN THEN
           BEGIN UnionSet(LRDSizeGeClass[i,j],LRDGeClass[i,j],Union);
            WHILE (j<=LRDNumGeClass[i]) AND (NOT(ObsConfig[i]>=Union) OR
                  (LRDOCGeClass[i,j]<>LRNOCGeClass[i,m]) OR
                  ClassesNotSame(i,m,j)) DO
             BEGIN j:=j+1;
              IF j<=LRDNumGeClass[i] THEN
               UnionSet(LRDSizeGeClass[i,j],LRDGeClass[i,j],Union)
             END;
            Reducible:= j <= LRDNumGeClass[i];
            IF Reducible THEN
             BEGIN
              FOR p:=j+1 TO LRDNumGeClass[i] DO
               BEGIN  LRDSizeGeClass[i,p-1]:=LRDSizeGeClass[i,p];
                LRDOCGeClass[i,p-1]:=LRDOCGeClass[i,p];
                FOR n:=1 TO LRDSizeGeClass[i,p-1] DO
                 LRDGeClass[i,p-1,n]:=LRDGeClass[i,p,n]
               END;
              LRDNumGeClass[i]:=LRDNumGeClass[i]-1
             END
           END;
          IF NOT Reducible OR NOT(ObsConfig[i] >= UnionN) THEN
           BEGIN k:=k+1; LRNSizeGeClass[i,k]:=LRNSizeGeClass[i,m];
            LRNOCGeClass[i,k]:=LRNOCGeClass[i,m];
            FOR n:=1 TO LRNSizeGeClass[i,k] DO
             LRNGeClass[i,k,n]:=LRNGeClass[i,m,n]
           END
         END; {of FOR m}
        LRNNumGeClass[i]:=k;
{Reduce summed numerator and denominator of LR}
        Sum(ObsConfig[i],LRNNumNumer[i],LRNNumGeClass[i],NNumNumer,
            LRNNumer[i],LRNOCNumer[i],LRNOCGeClass[i],NNumer,NOCNumer,
            LRNSizeGeClass[i],NNumerInd,LRNGeClass[i]);
        Sum(ObsConfig[i],LRDNumNumer[i],LRDNumGeClass[i],DNumNumer,
            LRDNumer[i],LRDOCNumer[i],LRDOCGeClass[i],DNumer,DOCNumer,
            LRDSizeGeClass[i],DNumerInd,LRDGeClass[i]);
        FOR p:=1 TO NNumNumer DO
         BEGIN  j:=1;
          WHILE (j<=DNumNumer) AND ((NOCNumer[p]<>DOCNumer[j]) OR
                (NNumer[p]<>DNumer[j])) DO j:=j+1;
          IF j<=DNumNumer THEN
           BEGIN n:=NNumerInd[p];
            IF n<=LRNNumNumer[i] THEN LRNNumer[i,n]:=[]
            ELSE LRNSizeGeClass[i,n-LRNNumNumer[i]]:=0;
            n:=DNumerInd[j];
            IF n<=LRDNumNumer[i] THEN LRDNumer[i,n]:=[]
            ELSE LRDSizeGeClass[i,n-LRDNumNumer[i]]:=0
           END
         END; {of FOR p}
        Delete(LRNNumNumer[i],LRNNumGeClass[i],LRNNumer[i],
         LRNOCNumer[i],LRNOCGeClass[i],LRNSizeGeClass[i],LRNGeClass[i]);
        Delete(LRDNumNumer[i],LRDNumGeClass[i],LRDNumer[i],
         LRDOCNumer[i],LRDOCGeClass[i],LRDSizeGeClass[i],LRDGeClass[i]);
{Reduce the whole LR if numerator is the same as denominator}
        Reducible:= (LRNNumNumer[i]=LRDNumNumer[i]) AND
                    (LRNNumDenom[i]=LRDNumDenom[i]) AND
                    (LRNNumGeClass[i]=LRDNumGeClass[i]);
        IF Reducible THEN
         BEGIN  j:=1;
          WHILE Reducible AND (j<=LRNNumNumer[i]) DO
           BEGIN  Reducible:=Reducible AND (LRNNumer[i,j]=LRDNumer[i,j])
                    AND (LRNOCNumer[i,j]=LRDOCNumer[i,j]);  j:=j+1
           END;
          j:=1;
          WHILE Reducible AND (j<=LRNNumDenom[i]) DO
           BEGIN  Reducible:=Reducible AND (LRNDenom[i,j]=LRDDenom[i,j])
                    AND (LRNOCDenom[i,j]=LRDOCDenom[i,j]);  j:=j+1
           END;
          j:=1;
          WHILE Reducible AND (j<=LRNNumGeClass[i]) DO
           BEGIN  Reducible:= Reducible
                    AND (LRNOCGeClass[i,j]=LRDOCGeClass[i,j]);  k:=1;
            WHILE ClassesNotSame(i,j,k) AND (k<=LRDNumGeClass[i]) DO
             k:=k+1;
            Reducible:=Reducible AND (k<=LRDNumGeClass[i]);  j:=j+1
           END;
          IF Reducible THEN
           BEGIN LRNNumNumer[i]:=0;LRDNumNumer[i]:=0;LRNNumDenom[i]:=0;
            LRDNumDenom[i]:=0; LRNNumGeClass[i]:=0; LRDNumGeClass[i]:=0
           END
         END {IF Reducible THEN}
       END {of FOR i}
    END
  END
END; {of MissHiMod}
