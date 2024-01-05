This version of AS204 is in Pascal.   It does not use AS 66 for the calculation
of normal probabilities, but one of Hastings' algorithms.

(* Program description: Journal of the Royal Statistical Society (Series C)
   Vol 33 No 3  1984, R.W. Farebrother,  Algorithm AS204, pp 332 - 339

   Use of Rubin's (1962) method to evaluate the
   expression  Pr[d k(i) K}(m(i),k(i)}) < c]   where k(i) and c are
   positive constants, and where K}(m(i),k(i)}) represents an independent
   chi-squared random variable with m(j) degrees of freedom and
   non-centrality parameter k(i)}.

     n =  number of chi-squared terms
     c =  Critical chi-squared value
     maxit = maximum number of iterations = 500 in Vol 33 No. 3 1984
     eps = degree of accuracy

     This includes an auxiliary algorithm, Function 26.2.17 from 'Handbook
     of Mathematical Functions', (1968), p 932, to obtain probabilities for
     the standard normal distribution.

     The program returns the cumulative probability value at the point c in
     the distribution.

*)
CONST
     n =  3;          (* number of terms in equation *)
     maxit = 500;     (* maximum number of iterations in eq. 3 of algorithm *)
     eps = 0.0001;    (* desired level of accuracy *)
     mode = 0.90625;  (* as in AS 204A  *)

VAR
    RUBEN, tol, beta, ao, aoinv   : Real;
    c, dnsty, z, rz               : Real;
    eps2, hold, hold2, sum, sum1  : Real;
    dans, lans, pans ,prbty, temp : Real;
    ifault, itemp, i, j, k, m     : Integer;
    lambda, delta                 : Array[1..n] of Real;
    mult                          : Array[1..n] of Integer;
    gamma, theta                  : Array[1..n] of Real;
    a, b                          : Array[1..maxit] of Real;
    exit, L                       : Boolean;

Procedure centnorm(VAR x,prob: Real);

const
  p  =  0.2316419;
  b1 =  0.319381530;
  b2 = -0.356563782;
  b3 =  1.781477937;
  b4 = -1.821255978;
  b5 =  1.330274429;

var
  constant, Z, t, Q        : real;
(*                         x       standardised value
                           Z       N(0,1) function value at x
                           Q       P(X > abs(x)) where X is N(0,1)
                        prob       P(-x < z < x)   *)

Begin
     constant := 1/sqrt(2*pi);
     If abs(x) < 4.5 then
     Begin
           if x < 0 then
                x := -x;
           Z := exp(-x*x/2)*constant;
           t := 1/(1+p*x);
           Q := ((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t*Z;
     end
     else
           Q := 0;
     prob := 1 - 2 * Q;
End;   (* centnorm  *)

BEGIN
     (* RUBEN evaluates the probability that a positive definite quadratic
        form in Normal variates is less than a given value *)

     c := 1.0;
     For i := 1 to n do
         delta[i] := 0.0;
     lambda[1] := 6;
     lambda[2] := 3;
     lambda[3] := 1;
     For i := 1 to n do
          mult[i] := 1;

     exit := false;
     L := false;
     If (n < 1) OR (c <= 0.0) OR (maxit < 1) OR (eps <= 0.0) then
     Begin
          RUBEN := -2.0;
          ifault := 2;
     End
     Else
     Begin
          tol := -200.0;
          (* preliminaries *)
          beta := lambda[1];
          sum := lambda[1];
          For i := 1 to n do
          Begin
               hold := lambda[i];
               If (hold <= 0.0) OR (mult[i] < 1) OR (delta[i] < 0.0) then
               Begin
                    RUBEN := -7.0;
                    ifault := -i;
                    exit := true;
                    i := n;
               End;
               If beta > hold then beta := hold;
               If sum < hold then sum := hold;
          End;
          If exit = false then
          Begin
               If mode > 0.0 then
               Begin
                    temp := mode * beta;
                    beta := temp;
               End
               Else
               Begin
                    temp := 2.0/(1.0/beta + 1.0/sum);
                    beta := temp;
               End;
               k := 0;
               sum := 1.0;
               sum1 := 0.0;
               for i := 1 to n do
               Begin
                    hold := beta/lambda[i];
                    gamma[i] := 1.0 - hold;
                    temp := 1;
                    for j := 1 to mult[i] do
                        temp := temp * hold;
                    sum := sum * temp;
                    sum1 := sum1 + delta[i];
                    k := k + mult[i];
                    theta[i] := 1.0;
               End;
               Writeln;
               ao := EXP(0.5 * (LN(sum) - sum1));
               If ao <= 0.0 then
               Begin
                    RUBEN := 0.0;
                    dnsty := 0.0;
                    ifault := 1;
               End
               Else
               Begin
                    z := c/beta;
                    (* evaluate probability and density of chi-squared on
                       k degrees of freedom.  The constant 0.22579135264473
                       is ln(sqrt(pi/2))  *)
                    itemp := (k DIV 2) * 2;
                    If k = itemp then
                    Begin
                         i := 2;
                         lans := -0.5 * z;
                         dans := EXP(lans);
                         pans := 1.0 - dans;
                    end
                    Else
                    Begin
                         i := 1;
                         lans := -0.5 * (z + LN(z)) - 0.22579135264473;
                         dans := EXP(lans);
                         rz := SQRT(z);
                         centnorm(rz,pans);
                         (* centnorm(x) will be a procedure or function to
                            evaluate the probability that a standard normal
                            variable lies in the interval [-x,x]
                            Algorithm AS66 recommended *)
                    End;
                    k := k - 2;
                    While i <= k do
                    Begin
                         If lans < tol then
                         Begin
                              lans := lans + LN(z/i);
                              dans := EXP(lans);
                         End
                         Else
                         Begin
                              temp := dans;
                              dans := temp * z/i;
                         End;
                         temp := pans;
                         pans := temp - dans;
                         i := i + 2;
                    End;
                    { evaluate successive terms of expansion  }
                    prbty := pans;
                    dnsty := dans;
                    eps2 := eps/ao;
                    aoinv := 1.0/ao;
                    sum := aoinv - 1.0;
                    For m := 1 to maxit do
                    Begin
                         sum1 := 0.0;
                         For i := 1 to n do
                         Begin
                              hold := theta[i];
                              theta[i] := hold * gamma[i];
                              hold2 := theta[i];
                              temp := hold2 * mult[i] + m * delta[i] * (hold - hold2);
                              sum1 := sum1 + temp;
                         End;
                         b[m] := 0.5 * sum1;
                         sum1 := b[m];
                         itemp := m - 1;
                         (* itemp = 0 when m = 1  !!!  *)
                         If itemp > 0 then
                         For i := itemp downto 1 do
                              sum1 := sum1 + b[i] * a[m-i];
                         a[m] := sum1/m;
                         sum1 := a[m];
                         k := k + 2;
                         If lans < tol then
                         Begin
                              lans := lans + LN(z/k);
                              dans := EXP(lans);
                         End
                         Else
                         Begin
                              temp := dans;
                              dans := temp * z/k;
                         End;
                         pans := pans - dans;
                         sum := sum - sum1;
                         temp := dnsty;
                         dnsty := temp + dans * sum1;
                         temp := sum1;
                         sum1 := pans * temp;
                         temp := prbty;
                         prbty := temp + sum1;
                         If prbty < -aoinv then
                         Begin
                              RUBEN := -3.0;
                              ifault := 3;
                              exit := true;
                              m := maxit;
                         End;
                         If exit = false then
                         Begin
                              temp := ABS(pans* sum);
                              If temp < eps2 then
                              Begin
                                   temp := ABS(sum1);
                                   If temp < eps2 then
                                   Begin
                                        ifault := 0;
                                        L := true;
                                        m := maxit;
                                   End;
                              end;
                         end;
                    End;  { m }
                    If L = false then ifault := 4;
                    If exit = false then
                    Begin
                         temp := dnsty;
                         dnsty := ao * temp/(beta + beta);
                         temp := prbty;
                         prbty := ao * temp;
                         If (prbty < 0.0) OR (prbty > 1.0) then
                            ifault := ifault + 5
                         Else
                             If dnsty < 0.0 then ifault := ifault + 6;
                         RUBEN := prbty;
                    end;
                end;
          End;
     End;
     Writeln('probability =  ',RUBEN:10:6);
     Writeln('fault code =  ',ifault:3);
END.


