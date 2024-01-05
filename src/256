FUNCTION imhofbd(lambda:realarray; mult:posintarray; delta:realarray;
         noterms:units; eps1,eps2:real; VAR ifault:nonnegint):real;

{Algorithm AS 256.1  Appl. Statist. (1990) Vol.39, No.2}

{Pascal translation of Koerts and Abrahamse's procedure
 for evaluating Imhof's upper bound}

CONST pi=3.14159265;
VAR   i:units;
      count,hold,range,sum1,sum2:real;

BEGIN
IF (noterms < 1) OR (eps1 <= 0.0) OR (eps2 <= 0.0) THEN
   BEGIN
   ifault := 2; imhofbd := -2.0
   END
ELSE
   BEGIN
   count := 0.0;
   sum1 := 0.0; sum2 := 0.0;
   FOR i := 1 TO noterms DO
      BEGIN
      hold := abs(lambda[i]);
      IF hold > eps2 THEN
         BEGIN
         count := count+mult[i];
         sum1 := sum1+mult[i]*ln(hold);
         END;
      sum2 := sum2+delta[i]
      END;
   IF count < 0.9 THEN
      BEGIN
      ifault := 4; imhofbd := -4.0
      END
   ELSE
      BEGIN
      count := 0.5*count;
      sum1 := 0.5*sum1+ln(pi*count);
      range := exp(-(sum1+0.5*sum2+ln(eps1))/count);
      IF sum2=0.0 THEN
      range := range+5.0/count
      ELSE
         BEGIN
            REPEAT
            sum2 := 0.0;
            FOR i := 1 TO noterms DO
               BEGIN
               hold := sqr(range*lambda[i]);
               sum2 := sum2+delta[i]*hold/(1.0+hold)
               END;
            hold := exp(sum1+count*ln(range)+0.5*sum2);
            IF hold*eps1 <= 1.0 THEN range := range+5.0/count ELSE count := 0.0
            UNTIL count=0.0;
         END
      END;
   imhofbd := range
   END
END; {of imhofbd}

FUNCTION imhofint(lambda:realarray; mult:posintarray; delta:realarray;
         noterms:units; arg,bound,eps3:real; VAR ifault:nonnegint):real;

{Algorithm AS 256.2  Appl. Statist. (1990) Vol.39, No.2}

{Pascal translation of Koerts and Abrahamse's procedure
 for evaluating Imhof's integral by Simpson's Rule}

CONST maxit=14;
      pi=3.14159265;
VAR   i,j,n:nonnegint;
      eps4,hold,int1,int2,step,sum1,sum2,sum4:real;
      cgd:boolean;

FUNCTION imhoffn(u:real):real;

   {imhoffn evaluates Imhof's integrand}

   VAR i:units;
       hold2,hold3,rho,sum,theta:real;

   BEGIN
   rho := 0.0;
   sum := 0.0;
   theta := -u*arg;
   FOR i := 1 TO noterms DO
      BEGIN
      hold := u*lambda[i];
      hold2 := sqr(hold);
      hold3 := 1.0+hold2;
      theta := theta+mult[i]*arctan(hold)+delta[i]*hold/hold3;
      sum := sum+delta[i]*hold2/hold3;
      rho := rho+mult[i]*ln(hold3)
      END;
   imhoffn := sin(0.5*theta)/(u*exp(0.5*sum+0.25*rho))
   END; {of imhoffn}

BEGIN
IF (noterms < 1) OR (bound <= 0.0) OR (eps3 <= 0.0) THEN
   BEGIN
   ifault := 2; imhofint := -2.0
   END
ELSE
   BEGIN
   ifault := 5;
   n := 2; step := 0.5*bound;
   eps4 := 3.0*pi*eps3;
   sum1 := -arg;
   FOR i := 1 TO noterms DO
   sum1 := sum1+lambda[i]*(mult[i]+delta[i]);
   sum1 := 0.5*sum1+imhoffn(bound);
   sum2 := 0.0;
   sum4 := imhoffn(step);
   int2 := (sum1+4.0*sum4)*step;
   i := 0; cgd := false;
      REPEAT
      i := i+1; n := n+n;
      step := 0.5*step;
      int1 := int2;
      sum2 := sum2+sum4;
      int2 := sum1+sum2+sum2;
      sum4 := 0.0; j := 1;
         REPEAT
         sum4 := sum4+imhoffn(j*step);
         j := j+2
         UNTIL j > n;
      int2 := (int2+4.0*sum4)*step;
      IF i > 3 THEN
         BEGIN
         IF abs(int1-int2) < eps4 THEN
            BEGIN
            IF abs(int2) > 1.5*pi THEN
            ifault := 6 ELSE ifault := 0;
            cgd := true
            END
         END;
      UNTIL (i=maxit) OR cgd;
   imhofint := 0.5-int2/(3.0*pi)
   END
END; {of imhofint}

FUNCTION imhof(lambda:realarray; mult:posintarray; delta:realarray;
noterms:units; arg,bound,eps1,eps2,eps3:real; VAR ifault:nonnegint):real;

{Algorithm AS 256.3  Appl. Statist. (1990) Vol.39, No.2}

{Pascal translation of Koerts and Abrahamse's implementation of
 Imhof's procedure for evaluating the probability that a diagonal
 form in normal variables is less than a given value, arg}

VAR i:units;

BEGIN

{check parameters}
IF (noterms < 1) OR (eps1 <= 0.0) OR (eps2 <= 0.0) OR (eps3 <= 0.0) THEN
   BEGIN
   ifault := 2; imhof := -2.0
   END
ELSE
   BEGIN
   ifault := 1; i := 1;
   WHILE ifault=1 DO
      BEGIN
      i := i+1;
      IF i=noterms THEN ifault := 0;
      IF (mult[i] < 1) OR (delta[i] < 0.0) THEN
         BEGIN
         ifault := 3; imhof := -i
         END
      END;

   {main body}
   IF ifault=0 THEN
      BEGIN
      IF bound <= 0.0 THEN
      bound := imhofbd(lambda,mult,delta,noterms,eps1,eps2,ifault);
      IF ifault=0 THEN
      imhof := imhofint(lambda,mult,delta,noterms,arg,bound,eps3,ifault)
      ELSE imhof := -ifault
      END
   END
END; {of imhof}


FUNCTION imhotep(diag,subd:realarray; noterms:units;
nc,arg,bound,eps0,eps1,eps2,eps3:real; VAR ifault:nonnegint):real;

{Algorithm AS 256.4  Appl. Statist. (1990) Vol.39, No.2}

{imhotep uses Imhof's procedure with real arithmetic
 to evaluate the probability that a tridiagonal form
 in normal variables is less than a given value, arg}

VAR i:units;
    hold:real;
    delta:realarray;
    mult:posintarray;

BEGIN
IF (noterms < 1) OR (nc < 0.0) OR
   (eps0 <= 0.0) OR (eps1 <= 0.0) OR (eps2 <= 0.0) OR (eps3 <= 0.0) THEN
   BEGIN
   ifault := 2; imhotep := -2.0
   END
ELSE
   BEGIN
   ifault := 0; delta[noterms] := 1.0;
   IF noterms > 1 THEN
      BEGIN
      FOR i := noterms-1 DOWNTO 1 DO delta[i] := 0.0;
      tql3(noterms,diag,subd,delta,eps0,ifault)
      END;

   IF ifault <> 0 THEN imhotep := -ifault
   ELSE
      BEGIN
      FOR i := 1 TO noterms DO
         BEGIN
         hold := delta[i]; mult[i] := 1;
         delta[i] := nc*sqr(hold)
         END;
      imhotep := imhof(diag,mult,delta,noterms,arg,bound,eps1,eps2,eps3,ifault)
      END
   END
END; {of imhotep}

FUNCTION sneekbd(diag,subd:realarray; noterms:units;
nc,eps1,eps2:real; VAR ifault:nonnegint):real;

{Algorithm AS 256.5  Appl. Statist. (1990) Vol.39, No.2}

{sneekbd evaluates Imhof's upper bound
 by means of Sturm sequences}

CONST pi=3.14159265;
VAR   i,j:nonnegint;
      hold,sum1,sum2:real;

PROCEDURE sturm(eta:real; VAR det:real; VAR mult:nonnegint);

   {sturm evaluates the determinant of a tridiasonal matrix
    and the number of roots less than a given value, eta}

   VAR i:units;
       count:nonnegint;
       real0,real1,real2:real;

   BEGIN
   count := 0;
   real1 := 1.0; real2 := 0.0;
   FOR i := 1 TO noterms DO
      BEGIN
      real0 := real1*(diag[i]-eta)-real2*sqr(subd[i]);
      real2 := real1; real1 := real0;
      IF real1*real2 < 0.0 THEN count := count+1
      END;
   det := real0;
   mult := count
   END; {of sturm}

BEGIN
IF (noterms < 1) OR (nc < 0.0) OR (eps1 <= 0.0) OR (eps2 <= 0.0) THEN
   BEGIN
   ifault := 2; sneekbd := -2.0
   END
ELSE
   BEGIN
   {trap small eigenvalues}
   sturm(-eps2,hold,i);
   sturm( eps2,hold,j);
   IF i <> j THEN
      BEGIN
      ifault := 7; sneekbd := -7.0
      END
   ELSE
      BEGIN
      {compute bound}
      ifault := 0;
      sturm(0.0,hold,j);
      sum1 := 0.5*noterms;
      sum2 := sqrt(abs(hold))*pi*sum1*eps1;
      sneekbd := exp(-(ln(sum2)+0.5*nc)/sum1)+5.0/sum1
      END
   END
END; {of sneekbd}

FUNCTION sneekint(diag,subd:realarray; noterms:units;
nc,arg,bound,eps3:real; VAR ifault:nonnegint):real;

{Algorithm AS 256.6  Appl. Statist. (1990) Vol.39, No.2}

{sneekint uses complex arithmetic to
 evaluate Imhof's integral by Simpson's Rule}

CONST maxit=14;
      eta=0.0001;
      pi=3.14159265;
VAR   i,j,n:nonnegint;
      eps4,int1,int2,step,sum1,sum2,sum4:real;
      cgd:boolean;

FUNCTION sneekfn(u:real):real;

   {sneekfn evaluates Imhof's integrand}

   VAR sgn:-1..1;
       i:units;
       hold1,hold2,imag0,imag1,imag2,pies,real0,real1,real2:real;

   BEGIN

   {evaluate determinant}
   real1 := 1.0; imag1 := 0.0;
   real2 := 0.0; imag2 := 0.0;
   sgn := 1; pies := 0.0;
   FOR i := 1 TO noterms DO
      BEGIN
      hold1 := u*diag[i]; hold2 := sqr(u*subd[i]);
      real0 := real1+imag1*hold1+real2*hold2;
      imag0 := imag1-real1*hold1+imag2*hold2;
      real2 := real1; real1 := real0;
      imag2 := imag1; imag1 := imag0;
      IF sgn*real0 < 0.0 THEN
         BEGIN
         sgn := -sgn;
         IF sgn*imag0 < 0.0
             THEN pies := pies+pi
             ELSE pies := pies-pi
         END
      END;

   {evaluate exponent}
   hold1 := -nc/(sqr(real1)+sqr(imag1));
   real0 := nc+(real1*real2+imag1*imag2)*hold1;
   imag0 := (real1*imag2-imag1*real2)*hold1;

   {evaluate integrand}
   hold1 := imag1/real1;
   hold2 := sqr(hold1)+1.0;
   hold2 := sqrt(hold2)*abs(real1);
   hold2 := exp(-0.5*real0)/sqrt(hold2);
   hold1 := arctan(hold1)+pies+imag0+u*arg;
   sneekfn := sin(-0.5*hold1)*hold2/u
   END; {of sneekfn}

BEGIN
IF (noterms < 1) OR (nc < 0.0) OR (bound <= 0.0) OR (eps3 <= 0.0) THEN
   BEGIN
   ifault := 2; sneekint := -2.0
   END
ELSE
   BEGIN
   ifault := 5;
   n := 2; step := 0.5*bound;
   eps4 := 3.0*pi*eps3;
   sum1 := sneekfn(eta*bound)+sneekfn(bound);
   sum2 := 0.0;
   sum4 := sneekfn(step);
   int2 := (sum1+4.0*sum4)*step;
   i := 0; cgd := false;
      REPEAT
      i := i+1; n := n+n;
      step := 0.5*step;
      int1 := int2;
      sum2 := sum2+sum4;
      int2 := sum1+sum2+sum2;
      sum4 := 0.0; j := 1;
         REPEAT
         sum4 := sum4+sneekfn(j*step);
         j := j+2
         UNTIL j > n;
      int2 := (int2+4.0*sum4)*step;
      IF i > 3 THEN
         BEGIN
         IF abs(int1-int2) < eps4 THEN
            BEGIN
            IF abs(int2) > 1.5*pi THEN
            ifault := 6 ELSE ifault := 0;
            cgd := true
            END
         END
      UNTIL (i=maxit) OR cgd;
   sneekint := 0.5-int2/(3.0*pi)
   END
END; {of sneekint}

FUNCTION sneek(diag,subd:realarray; noterms:units;
nc,arg,bound,eps1,eps2,eps3:real; VAR ifault:nonnegint):real;

{Algorithm AS 256.7  Appl. Statist. (1990) Vol.39, No.2}

{sneek uses Imhof's procedure with complex arithmetic
 to evaluate the probability that a tridiagonal form
 in normal variables is less than a given value, arg}

BEGIN
IF (noterms < 1) OR (nc < 0.0) OR
   (eps1 <= 0.0) OR (eps2 <= 0.0) OR (eps3 <= 0.0) THEN
   BEGIN
   ifault := 2; sneek := -2.0
   END
ELSE
   BEGIN
   ifault := 0;
   IF bound <= 0.0 THEN
   bound := sneekbd(diag,subd,noterms,nc,eps1,eps2,ifault);
   IF ifault=0 THEN
   sneek := sneekint(diag,subd,noterms,nc,arg,bound,eps3,ifault)
   END
END; {of sneek}
