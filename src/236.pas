{ 
        Algorithm AS236  Appl. Statist. (1988) Vol. 37, No. 2
                                                                            }

{  To Enumerate all RxC tables having given row and column sums

   To Evaluate the conditional probability of each table under
   hypergeometric and multinomial models

   To Evaluate likelihood across all tables according to the multinomial
   model                                                                    }



{ ************************************************************************* }

{                          DATA STRUCTURES                                  }

{ ************************************************************************* }

  CONST
  maxtablesize=10;               { maximum dimensions of table }
  maxN=300;                      { maximum global total of table entries - the
                                  values 10 and 300 may be changed, if desired }

  TYPE
  intvalues=0..maxN;             { permissible table entry type }
  faultint=0..123;               { fault indicator type
                                        0,100,20,3,120,103,23,123  }
  dimension=1..maxtablesize;     { table dimension type }
  vector = ARRAY[1..maxtablesize] OF intvalues;{ row and column vector type }
  matrix = ARRAY[1..maxtablesize,1..maxtablesize] OF real; { matrix type }

  return = ^returnlist;          { pointer to function returns }
  head = ^headlist;              { pointer to NODE information }
  table  =  ^tablist;  {  pointer to NODULE.  All  tables  are stored in
                         in a single linked list of head and table records. }
  tablist = RECORD               { Each table/nodule record contains the
                                            value of the:-      }
              entry:intvalues;   { current table entry }
              loghyper:real;     { current contribution to hypergeometric }
              multi:real;        { current contribution to multinomial }
              rest: head ;       { a pointer to the rest of the family of
                                   tables built onto the current entry
                                   ie to the next node }
               next:table        { a pointer to the next possible entry value
                                   for the current position and associated
                                   family of tables ie to the next nodule }
            END;


  headlist = RECORD
               number:integer;   { number of subtables emanating from the
                                   current node }
               position:ARRAY[1..2] OF dimension;  { current table position }
               first:table       { pointer to the first cell value and
                                  associated tables ie to the first nodule }
             END;

  returnlist = RECORD
                 tables:head;    {  a function call returns a pointer to a  }
                 value:real      {  family of tables and the value of the
                                    likelihood function evaluated across these }
               END;

  VAR
  rowsum:vector;                 { required row sums for tables }
  colsum:vector;                 { required column sums for tables }
  rows,cols:dimension;           { number of rows and columns in table }
  pmatrix:matrix;                { proposed transition matrix }
  result:return;                 { variable to receive function return }
  loghypercons,multicons:real;   { constants for hypergeometric and
                                   multinomial probabilities }
  ifault:faultint;               { fault indicator }


{ ************************************************************************* }

{                        OUTPUT  PROCEDURES                                 }

{ ************************************************************************* }

PROCEDURE Report(r:return;i:faultint);

  TYPE
  tabmatrix = ARRAY[1..maxtablesize,1..maxtablesize] OF intvalues;

  PROCEDURE Printmat(tub:tabmatrix;prob1,prob2:real);

    VAR
    i,j:dimension;

    BEGIN                             { write out probabilities and current
                                                                      table }
    Writeln;
    Writeln('Hypergeometric probability = ',prob1);
    Writeln('Multinomial probability = ',prob2);
    FOR i:=1 TO rows DO
      BEGIN
      Writeln;
      Write('row ',i:1,' *  ');
      FOR j:=1 TO cols DO
        Write(tub[i,j],' ')
      END;

    Writeln
    END;

    PROCEDURE Printtables(v:real;example:head);

      CONST
      assumedzero=1.0e-15;

      VAR
      tab:tabmatrix;
      exam:table;
      i,j:dimension;
      hyperprob,multiprob:real;

      BEGIN
      WITH example^ DO
        BEGIN
        exam:=first;
        i:=position[1];
        j:=position[2]
        END;
      REPEAT
        WITH exam^ DO
          BEGIN
          tab[i,j]:=entry;            { copy current table entry into
                                         matrix for output }
          IF (i=rows) AND (j=cols) THEN { if table complete print out }
            BEGIN
            hyperprob:=Exp(loghyper+loghypercons);
            IF v>assumedzero THEN multiprob:=multi/v
            ELSE multiprob:=0.0;
            Printmat(tab,hyperprob,multiprob)
            END
          ELSE Printtables(v,rest)    { assemble rest of table through
                                        a recursive call to Printtables }
          END;
        exam:=exam^.next;
      UNTIL exam=NIL
      END;

  BEGIN {  Report  }
  Writeln;
  Writeln('IFAULT=',i:1);Writeln;
  IF i=0 THEN                         { if evaluation successful print out
                                             all tables and probabilities }
    BEGIN
    Writeln('Number of Tables Generated = ',r^.tables^.number);
    Writeln('Multinomial likelihood = ',r^.value*multicons);
    Printtables(r^.value,r^.tables)
    END
  END; {  Report  }



{ ************************************************************************* }

{                          PRINCIPAL FUNCTION                               }

{ ************************************************************************* }

FUNCTION Enum(rows,cols:dimension; rowsum,colsum:vector;
             pmatrix:matrix; VAR ifault:faultint):return;

  CONST
  oneminus=0.99999;
  oneplus=1.00001;

  TYPE
  dvector = ARRAY[0..maxN] OF real;   { real vector type }

  VAR
  logfact:dvector;
  sfault,tfault:faultint;
  i,j:dimension;
  rowtotal,coltotal:intvalues;
  total:real;

  FUNCTION max(a,b:integer):integer;
    BEGIN
    IF a>b THEN max:=a ELSE max:=b
    END; {  max  }

  FUNCTION min(a,b:integer):integer;
    BEGIN
    IF a<b THEN min:=a ELSE min:=b
    END; {  min  }

  FUNCTION Eval(n:intvalues;x:real):real;
    CONST assumedzero=1.0e-15;
    BEGIN
    IF n=0 THEN Eval:=1.0 ELSE
      IF x<assumedzero THEN Eval:=0.0 ELSE
        Eval:=Exp(n*Ln(x)-logfact[n])
    END; {  Eval  }

  PROCEDURE Evalconst(rowtotal:intvalues);
    VAR
    i:intvalues;
    logprod1,logprod2:real;
    BEGIN
    logfact[0]:=0.0;
    FOR i:=1 TO rowtotal DO
      logfact[i]:=logfact[i-1]+Ln(i);
    logprod1:=0.0;
    FOR i:=1 TO rows DO
      logprod1:=logprod1+logfact[rowsum[i]];
    logprod2:=0.0;

    FOR i:=1 TO cols DO
      logprod2:=logprod2+logfact[colsum[i]];
    loghypercons:=logprod1+logprod2-logfact[rowtotal];
    multicons:=Exp(logprod1)
    END; {  Evalconst  }


  FUNCTION RxC(rno,cno:dimension;rt,ct:intvalues;hp,mp:real;
                                  rowsum,colsum:vector):return;
    VAR
    count,val,upper,lower,newrt,newct:intvalues;
    curloghyper,curmulti,newhp,newmp,total:real;
    newrowsum:vector;
    newcolsum:vector;
    newhead:head;
    result,newresult:return;
    newtable,nextable:table;

    BEGIN  {  RxC  }

    New(newtable);
    newrowsum:=rowsum;
    newcolsum:=colsum;

    IF (rno<rows) AND (cno<cols) THEN {  current position rno,cno
                                                          in body of table }
      BEGIN
      nextable:=newtable;
      upper:=min(rowsum[rno],colsum[cno]);
      lower:=max(0,rowsum[rno]+colsum[cno]-ct);
      newct:= ct-colsum[cno];
      count:=0;
      total:=0.0;
      FOR val:= lower TO upper DO     { step 2 - consider all possible values
                                        for the current position }
        BEGIN
        WITH nextable^ DO
          BEGIN
          entry:=val;
                                      { evaluate contributions to likelihood }
          multi:=Eval(val,pmatrix[rno,cno]);
          loghyper:=-logfact[val];
          newrowsum[rno]:=rowsum[rno]-val;
          newcolsum[cno]:=colsum[cno]-val;
          newrt:=rt-val;
          newhp:=hp+loghyper;
          newmp:=mp*multi;
          New(result);
                                      { step 3 - obtain results for remainder
                                        of the table by recursive call to RxC }
          result:=RxC(rno,cno+1,newrt,newct,newhp,newmp,newrowsum,newcolsum);
          rest:=result^.tables;


          count:=count+rest^.number;
          total:=total+result^.value*multi;
          IF val<upper THEN
            BEGIN
            New(next);
            nextable:=next
            END
          ELSE next:=NIL
          END
        END
      END

    ELSE
      BEGIN
      New (result);

      IF rno<rows THEN                { current position in right-hand
                                                                     margin }
        BEGIN
        val:=rowsum[rno];             { only one possible value }
        newrowsum[rno]:=0;
        newcolsum[cno]:=colsum[cno]-val;
        newrt:=rt-val;
        newct:=newrt;
                                      { evaluate contributions to likelihood }
        curloghyper:=-logfact[val];
        curmulti:=Eval(val,pmatrix[rno,cno]);
        newhp:=hp+curloghyper;
        newmp:=mp*curmulti;
                                      { step 3 - obtain results for remainder
                                        of the table by recursive call to RxC }
        result:=RxC(rno+1,1,newrt,newct,newhp,newmp,newrowsum,newcolsum);
        count:=result^.tables^.number;
        total:=result^.value*curmulti
        END

      ELSE
        IF cno<cols THEN              { current position in last row }
          BEGIN
          val:=colsum[cno];           { only one possible value }
          newcolsum[cno]:=0;
          newrowsum[rno]:=rowsum[rno]-val;
          newrt:=rt-val;
          newct:=newrt;
                                      { evaluate contributions to likelihood }
          curloghyper:=-logfact[val];
          curmulti:=Eval(val,pmatrix[rno,cno]);
          newhp:=hp+curloghyper;
          newmp:=mp*curmulti;
                                      { step 3 - obtain results for remainder
                                        of the table by recursive call to RxC }
          result:=RxC(rno,cno+1,newrt,newct,newhp,newmp,newrowsum,newcolsum);
          count:=result^.tables^.number;

          total:=result^.value*curmulti
          END

        ELSE                          { final cell }
          BEGIN
          val:=rt;                    { only one possible value }
          result^.tables:=NIL;        { table complete }
          count:=1;
                                      { evaluate contributions to likelihood }
          curloghyper:=hp-logfact[val];
          curmulti:=Eval(val,pmatrix[rno,cno]);
          total:=curmulti;
          curmulti:=curmulti*mp
          END;

      WITH newtable^ DO
        BEGIN
        entry:=val;
        loghyper:=curloghyper;
        multi:=curmulti;
        rest:=result^.tables;
        next:=NIL
        END
      END;

    New(newhead);                     { assemble entries for current NODE }
    WITH newhead^ DO
      BEGIN
      number:=count;
      position[1]:=rno;
      position[2]:=cno;
      first:=newtable
      END;
    New(newresult);
    WITH newresult^ DO
      BEGIN
      tables:=newhead;
      value:=total
      END;
    RxC:=newresult
    END; {  RxC  }

  BEGIN {  Enum  }
                                      { check validity of input values }
  ifault:=0;
                                      { totals of row and column sums must be
                                                    equal }
  rowtotal:=0;
  coltotal:=0;
  FOR i:= 1 TO rows DO
    rowtotal:=rowtotal+rowsum[i];
  FOR i:=1 TO cols DO
    coltotal:=coltotal+colsum[i];

  IF rowtotal<>coltotal THEN ifault:=1;

                                      { pmatrix entries must be probabilities
                                          and rows must sum to one }
  sfault:=0;
  tfault:=0;
  FOR i:=1 TO rows DO
    BEGIN
    total:=0.0;
    FOR j:=1 TO cols DO
      BEGIN
      IF (pmatrix[i,j]<0.0) OR (pmatrix[i,j]>1.0) THEN sfault:=1;
      total:=total+pmatrix[i,j]
      END;
    IF (total<oneminus) OR (total>oneplus) THEN tfault:=1
    END;
  ifault:=100*ifault+20*sfault+3*tfault;

                                      { if valid input, evaluate tables and
                                                 probabilities }
  IF ifault=0 THEN
    BEGIN
    Evalconst(rowtotal);
    Enum:=RxC(1,1,rowtotal,coltotal,0.0,1.0,rowsum,colsum)
    END
  ELSE
    Enum:=NIL
  END;  {  Enum  }
