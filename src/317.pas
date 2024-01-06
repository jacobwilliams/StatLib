program mixgood(input,output,data,results);
{
  Program to calculate the maximum likelihood estimates of a mixture of
  distributions and carry out the goodness of fit tests
}

const limit = 300;

type rarray = array[1..limit] of real;
     iarray = array[1..limit] of integer;
     darray = array[1..limit,1..limit] of real;

var a,k,m,c,sample_size,ifault: integer;
    reply: char;
    fullyclose: boolean;
    data: text; {text = file of char}
    results: text; {text = file of char}
    filename,filename1: string;
    alpha,mu,sigma: rarray;
    n,x: rarray;
    f: darray;

procedure mix;

{
  procedure to calculate the maximum likelihood estimates of a mixture of
  distributions (normal, exponential, poisson or binomial). The log-likelihood
  function and the number of iterations taken to satisfy the tolerance value
  are also calculated.

  Note that the following may depend on the compiler/operating system used:

  - the naming of "input" and "output" in the program name statement
  - the use of the statements "open", "close", "reset" and "rewrite" (and their
    parameters) to open and close the data and results files

  The format of the data file is:

      a, k, m, c                 type of distribution
				  a = 1: normal with K-S test
				      2: exponential with K-S test
				      3: poisson with K-S test
				      4: binomial with K-S test
				      5: normal with Chi-squared test
				      6: exponential with Chi-squared test
				      7: poisson with Chi-squared test
				      8: binomial with Chi-squared test)
				 k = number of component distributions
				 m = number of classes
				 c = number of iterations before each set
				     of printed results (c = 0 suppress
				     printing)
      x[1], n[1]                 first observed value, frequency
	:     :
      x[m], n[m]                 last observed value, frequency
				    (note that x[i+1] >= x[i])
      tol                        tolerance for first run
      alpha[1], mu[1], sigma[1]  proportion/mean/standard deviation of each
	  :       :        :        component distribution (sd is only
      alpha[k], mu[k], sigma[k]     required for the normal distribution)
      tol                        tolerances for subsequent runs (the final
       :                            tolerance should be <= 0.0)
      trials                     trials only needed for binomial mixture
      u[1]                       u[] are the upper class boundaries for the
	:                           Normal and Exponential mixtures only.
      u[m] or u[n]                  In the cases of the K-S test, these will
				    be all the x[1]..x[n] original values
				    arranged in ascending order
}

var i,j,counter: integer;
    tol,logl: real;

procedure putout(logl: real);

var j: integer;

begin
     for j := 1 to k
     do begin
	write(results,'Distribution  ');
	write(results,j:3,'    ',alpha[j]:10:6,mu[j]:10:6);

	if (a = 1) or (a = 5)
	then writeln(results,sigma[j]:10:6)
	else writeln(results);

	writeln(results)
	end;

     writeln(results,'Log likelihood ',logl:14:6);
     writeln(results)
end; {putout}

procedure mixture(tol: real; var logl: real; var counter: integer);

var i,j: integer;
    test: boolean;
    sumalpha,part,oldlogl: real;
    newalpha,newmu,newsigma,dt,nt,vt,g: rarray;
{    f: darray;}
    temp1,temp2: real;

begin
     {this algorithm calculates the maximum likelihood estimates of the
	parameters of a mixture of normal or exponential or poisson or
	binomial distributions. these parameters are the mixing proportions,
	the means and (in the normal distribution case) standard deviations.
	it also calculates the log-likelihood function and the number of
	iterations taken to satisfy the tolerance value.

	The returned values of ifault are:

	0 : no errors
	1 : distribution type (a) should be in the range 1 - 8
	2 : proportion (alpha) should be greater than 0.0 and less than 1.0
	3 : mean should be within the range of observed values (x)
	4 : standard deviation (sigma) should be greater than 0.0
	5 : observed values (x) should be in ascending order
	6 : number of observations should be at least twice the number of
	    distinct observed values (m)
	7 : an observed value (x) should be greater than or equal to 0.0
	8 : each mean should have a distinct value
	9 : each standard deviation (sigma) should have a distinct value if the
	    corresponding means have the same value}

     ifault := 0;

     if (a < 1) or (a > 8)
     then ifault := 1;

     if ifault = 0
     then begin

	  for i := 1 to (m-1)
	  do if (ifault = 0) and (x[i] > x[i+1])
	     then ifault := 5;

	  if (ifault = 0) and (sample_size < (2*m))
	  then ifault := 6;

	  if (ifault = 0)
	  then for i := 1 to m
	       do begin

		  if (ifault = 0) and (n[i] < 0)
		  then ifault := 6;

		  if (ifault = 0) and ((a <> 1) and (a <> 5))
		  then if x[i] < 0.0
		       then ifault := 7;

		  end;

	  oldlogl := 0.0;
	  counter := 0;
	  test := true;

	  while test and (ifault = 0)
	  do begin
	     counter := counter + 1;

	     for j := 1 to k
	     do begin

		if ((alpha[j] > 1.0) or (alpha[j] <= 0.0)) and (ifault = 0)
		then ifault := 2;

		if ((mu[j] > x[m]) or (mu[j] < x[1])) and (ifault = 0)
		then ifault := 3;

		if (ifault = 0) and ((a = 1) or (a = 5))
		then if sigma[j] <= 0.0
		     then ifault := 4;

		end;

	     if ifault = 0
	     then for i := 1 to (k-1)
		  do for j := (i+1) to k
		     do if (mu[i] = mu[j]) and (ifault = 0)
			then if (a = 1) or (a = 5)
			     then if sigma[i] = sigma[j]
				  then ifault := 9
			     else ifault := 8;

	     if ifault = 0
	     then begin
		  logl := 0.0;

		  for i := 1 to m
		  do begin
		     g[i] := 0.0;

		     for j := 1 to k
		     do begin
			{note that a**b = exp(b * log(a))}

			if (a = 1) or (a = 5)
			then begin
			     temp1 := (x[i]-mu[j])/sigma[j];
			     f[i,j] := exp(-((temp1*temp1)/2))/sigma[j]
			     end;

			if (a = 2) or (a = 6)
			then f[i,j] := (1/mu[j])*exp(-x[i]/mu[j]);

			if (a = 3) or (a = 7)
			then if x[i] = x[1]
			     then begin
				  temp1 := exp(x[i] * ln(mu[j]));
				  f[i,j] := exp(-mu[j])*temp1
				  end
			     else f[i,j] := f[i-1,j]*mu[j];

			if (a = 4) or (a = 8)
			then if x[i] = x[1]
			     then begin
				  temp1 := exp(x[m] * ln(1-mu[j]/x[m]));
				  temp2 := exp(x[i] * ln((mu[j])/(x[m]-mu[j])));
				  f[i,j] := temp1 * temp2
				  end
			     else f[i,j] := f[i-1,j]*(mu[j]/(x[m]-mu[j]));

			g[i] := g[i] + alpha[j]*f[i,j]
			end;

		     logl := logl + n[i]*ln(g[i])
		     end
		  end;

	     {calculation of the probability densities of the subpopulations
		which comprise the mixture and the log likelihood function}

	     if ifault = 0
	     then begin
		  test := false;
		  sumalpha := 0.0;

		  for j := 1 to k
		  do begin
		     nt[j] := 0;
		     dt[j] := 0;
		     vt[j] := 0;

		     for i := 1 to m
		     do begin
			part := f[i,j] * n[i] / g[i];
			dt[j] := dt[j] + part;
			nt[j] := nt[j] + part * x[i];

			if (a = 1) or (a = 5)
			then begin
			     temp1 := (x[i] - mu[j]);
			     vt[j] := vt[j] + part * temp1 * temp1
			     end

			end;

		{calculation of denominators and numerators of new estimates}
		     newmu[j] := nt[j] / dt[j];

		     if j <> k
		     then begin
			  newalpha[j] := alpha[j] * dt[j] / sample_size;
			  sumalpha := sumalpha + newalpha[j]
			  end
		     else newalpha[k] := 1 - sumalpha;

		     if (a = 1) or (a = 5)
		     then newsigma[j] := sqrt(vt[j] / dt[j]);

		     if abs(oldlogl - logl) > tol
		     then test := true;

		     oldlogl := logl;
		     alpha[j] := newalpha[j];
		     mu[j] := newmu[j];

		     if (a = 1) or (a = 5)
		     then sigma[j] := newsigma[j]
		     end;

		     if c > 0
		     then if (counter div c) * c = counter
			  then putout(logl)
		     end
		 end

	  end

end; {mixture}

begin {mix}

     {read the data file and record the data in the results file}
     readln(data,a,k,m,c);
     case a of
     1,5: begin
	  writeln(results,'Normal Mixture');
	  writeln(results,'**************')
	  end;

     2,6: begin
	  writeln(results,'Exponential Mixture');
	  writeln(results,'*******************')
	  end;

     3,7: begin
	  writeln(results,'Poisson Mixture');
	  writeln(results,'***************')
	  end;

     4,8: begin
	  writeln(results,'Binomial Mixture');
	  writeln(results,'****************')
	  end
     end;

     writeln(results);
     writeln(results,'Number of component distributions = ',k:4);
     writeln(results,'Number of classes = ',m:4);
     writeln(results);
     writeln(results,'   Value   Frequency');
     writeln(results);
     sample_size := 0;

     for i := 1 to m
     do begin
	readln(data,x[i],n[i]);
	writeln(results,x[i]:8:2,'  ',round(n[i]):10);
	sample_size := sample_size + round(n[i])
	end;

     writeln(results);
     write(results,'Sample size =  ');
     writeln(results,sample_size:8);
     writeln(results); write(results,'------------------------------');
     writeln(results,'------------------------------');writeln(results);
     write(results,'Tolerance = ');
     readln(data,tol);
     writeln(results,tol:10:8);
     writeln(results);
     writeln(results,'                       Initial estimates ');
     writeln(results);
     write(results,'                      Proportion     Mean ');

     if (a = 1) or (a = 5)
     then writeln(results,'     S.D. ')
     else writeln(results);

     writeln(results);

     for j := 1 to k
     do begin
	write(results,'Distribution  ');
	write(results,j:3,'    ');
	read(data,alpha[j],mu[j]);
	write(results,alpha[j]:10:6,mu[j]:10:6);

	if (a = 1) or (a = 5)
	then begin
	     readln(data,sigma[j]);
	     writeln(results,sigma[j]:10:6)
	     end
	else begin
	     readln(data);
	     writeln(results)
	     end;

	writeln(results)
	end;

     repeat
	if (c > 0)
	then begin
	     writeln(results,'                    Successive approximations of estimates');
	     writeln(results)
	     end;

	mixture(tol,logl,counter);

	writeln(results,'                    Calculated maximum likelihood estimates');
	writeln(results);

	putout(logl);

	write(results,'Number of iterations ');
	writeln(results,counter:5);
	writeln(results);

	if ifault <> 0
	then writeln(results,'Error ',ifault:3,' in mixture');

	readln(data,tol);
	if (tol > 0.0)
	then begin
	     write(results,'------------------------------');
	     writeln(results,'------------------------------');writeln(results);
	     writeln(results,'Input new tolerance or 0 to finish');
	     writeln(results);
	     write(results,'Tolerance = ');
	     writeln(results,tol:10:8);
	     writeln(results)
	     end
     until tol <= 0.0;

end; {mix}

procedure goodness;
{
 This algorithm carries out the Kolmogorov-Smirnov and the Chi-squared
 goodness-of-fit tests to mixtures of distributions of the same type. The
 distributions considered in the mixture are the normal, exponential,
 Poisson and binomial.
}

const
     ltone = 7.0;
     utzero = 18.47;
     zero = 0.0;
     half = 0.5;
     one = 1.0;
     con = 1.2;

var
   sum,res,highest1,highest2,normality,ks,ch: real;
   d,delta2,expf,expected1,u: rarray;
   percent,test,count,i,j,trials: integer;

procedure read_data1;

begin
     write(results,'------------------------------');
     writeln(results,'------------------------------');
     writeln(results);
     writeln(results,'Goodness of Fit');
     writeln(results,'***************');
     writeln(results);

     {read number of trials for binomial distribution only}
     if (a = 4) or (a = 8)
     then begin
	  read(data,trials);
	  writeln(results,'Number of trials = ',trials:3);
	  writeln(results)
	  end;
end; {read_data1}

procedure read_data2;

{reads in u-values for use in the binomial, exponential, poisson and
 normal distributions to calculate the expected probabilities or
 frequencies}

var
   t: integer;

begin
     if (a = 1) or (a = 2)
     then for t := 1 to sample_size
	  do read(data,u[t])
     else if (a = 5) or (a = 6)
     then for t := 1 to m
	  do read(data,u[t])
     else for t := 1 to m
	  do u[t] := x[t]
end; {read_data2}

procedure check;

{procedure checks that the sums of the proportions Alpha, is not outside
 the range from 0.0 to 1.0.}

var
   loop: integer;
   sum1: real;

begin
     sum1 := 0.0;

     for loop := 1 to k
     do sum1 := sum1 + alpha[loop];

     if (sum1 > 1.0) or (sum < 0.0)
     then begin
	  writeln(results);
	  writeln(results,'Ifault = 2 : The sum of the proportions should not be');
	  writeln(results,'greater than 1.0 or less than 0.0.');
	  writeln(results,'The sum of your proportions (Alpha) = ',sum1:3:2);
	  ifault := 1
	  end
end; {check}

procedure check1;

{this procedure checks that the sample values are arranged in ascending
 order, for the normal and exponential distributions, for calculation of
 their cumulative probabilities}

var
   loop: integer;

begin
     if (a = 1) or (a = 2)
     then begin
	  for loop := 1 to (sample_size - 1)
	  do begin
	     if (u[loop] > u[loop + 1])
	     then begin
		  writeln(results);
		  write(results,'Sample values have not been arranged in ascending order');
		  writeln(results);
		  ifault := 1
		  end
	     end;
	  end;
end; {check1}

procedure alnorm(var delta1: rarray; upper: boolean);

{adapted from algorithm AS66, calculates from minus infinity to u[i].
 This procedure is used in the calculation of the cumulative
 probabilities and frequencies for the normal distribution}

var
   z,y: real;
   up: boolean;

begin
     up := upper;
     z := ((u[i] - mu[j]) / sigma[j]);

     if (z < zero)
     then begin
	  up := not up;
	  z := -z;
	  end;

     y := half * (z * z);

     if (z > con)
     then delta1[j] := 0.398942280385 * exp(-y) /
		       (z - 3.8052e-8 + 1.00000615302 /
		       (z + 3.98064794e-4 + 1.98615381364 /
		       (z - 0.151679116635 + 5.29330324926 /
		       (z + 4.8385912808 - 15.1508972451 /
		       (z + 0.742380924027 + 30.789933034 /
		       (z + 3.99019417011))))))
     else delta1[j] := half - z * (0.398942280444 - 0.399903438504 * y /
		       (y + 5.75885480458 - 29.8213557808 /
		       (y + 2.62433121679 + 48.6959930692 /
		       (y + 5.92885724438))));

     if (up)
     then delta1[j] := one - delta1[j]
     else if (z > ltone) and (z <= utzero)
	  then delta1[j] := 0
	  else if (z > ltone) and (z > utzero)
	       then delta1[j] := 0
	       else if (up) and (z > utzero)
		    then begin
			 delta1[j] := 0;
			 delta1[j] := one - delta1[j];
			 end;
end; {alnorm}

procedure answer1(var sum: real);

begin
     sum := sum + (d[j] * alpha[j]);
end; {answer1}

procedure exponential(var sum: real);

{probability equation of the exponential distribution}

begin
     d[j] := alpha[j] * (1 - (exp(-u[i] / mu[j])));
     sum := sum + d[j];
end; {exponential}

procedure poiss1(var delta1: rarray);

{probabilty equation of the poisson distribution}

begin
     delta1[j] := alpha[j] * (exp(-mu[j]));
     count := 0;

     while (count < u[i])
     do begin
	count := count + 1;
	delta1[j] := (1 / count) * delta1[j] * mu[j];
	end;

     sum := sum + delta1[j];
end; {poiss1}

procedure binomial(var res: real);

{calculates (n) for use in the binomial distribution
	    (u)                (procedure binomial1)}

var
   loop,loop1: integer;
   start: real;

begin
     loop := 1;
     start := 1;

     if (u[i] = 0)
     then res := 1
     else begin
	  repeat
	  start := start * ((trials - loop) + 1);
	  loop := loop + 1;
	  until ((trials - loop) < count);

	  res := start;
	  loop := round(u[i]);

	  for loop1 := 1 to loop
	  do res := res * (1 / loop1);
	  end
end; {binomial}

procedure binomial1(var delta1: rarray; bu: rarray);

var
   loop: integer;

begin
     bu[j] := bu[j] / trials;
     delta1[j] := alpha[j];
     loop := 0;

     while (loop < u[i])
     do begin
	loop := loop + 1;
	delta1[j] := delta1[j] * bu[j];
	end;

     bu[j] := 1.0 - bu[j];
     loop := 0;

     while (loop < count)
     do begin
	loop := loop + 1;
	delta1[j] := delta1[j] * bu[j];
	end;

     sum := sum + res * (delta1[j]);
end; {binomial1}

procedure pdf;

{calculates first class interval for normal distribution from minus
 infinity and for exponential distribution from 0. Also the other class
 intervals are calculated by subtraction}

var
   loop,loop1: integer;

begin
     loop1 := 1;
     expected1[1] := expf[1];

     for loop := 2 to m
     do begin
	expected1[loop] := expf[loop1 + 1] - expf[loop1];
	loop1 := loop1 + 1;

	if (expected1[loop] < 0.0)
	then expected1[loop] := 0.0;
	end;
 end; {pdf}

procedure delist(var expected: rarray; expected1: rarray);

{calculates the last class interval for all distributions so that the
 area underneath the curve is equal to 1}

var
   loop,loop1: integer;
   sum2: real;

begin
     sum2 := 0.0;

     for loop1 := 1 to m
     do begin
	expected1[loop1] := expected1[loop1] * sample_size;
	sum2 := sum2 + expected1[loop1];
	end;

     if sum2 < sample_size
     then expected1[m] := sample_size - sum2 + expected1[m];

     for loop := 1 to m
     do begin
	expected[loop] := expected1[loop];

	if (a = 7) or (a = 8)
	then writeln(results,u[loop]:6:2,'      ',expected[loop]:14:2)
	else writeln(results,u[loop]:10:2,'        ',expected[loop]:14:2)

	end;
end; {delist}

procedure prob_data;

{reads in data cumulative probabilities for either the normal, poisson
 exponential and binomial distributions. The empirical values are
 calculated (depending if a = 1 or 2) 1:continuous or 2:discrete and
 stored in array observed}

var
   loop,upper: integer;
   sum: real;

begin
     if (a = 1) or (a = 2)
     then upper := sample_size
     else upper := m;

     for loop := 1 to upper
     do if (a = 1) or (a = 2)
	then n[loop] := (loop / sample_size)
     else begin
	  sum := sum + n[loop] / sample_size;
	  n[loop] := sum
	  end
end; {prob_data}

procedure test1;

var
   loop: integer;

begin
     if (a = 1) or (a = 2)
     then for loop := 1 to sample_size
	  do begin
	     d[loop] := n[loop] - expf[loop];
	     delta2[loop] := expf[loop] - (n[loop] - (1 / sample_size));
	     end
     else for loop := 1 to m
	  do begin
	     d[loop] := n[loop] - expf[loop];
	     delta2[loop] := expf[loop] - (n[loop] - (1 / m));
	     end
end; {test1}

procedure search;

{searches for highest positive numbers from delta1 and delta2 to
 calculate the maximum discrepancy}

var
   loop, upper: integer;

begin
     if (a = 1) or (a = 2)
     then upper := sample_size - 1
     else upper := m - 1;

     if d[1] >= 0
     then highest1 := d[1];

     if delta2[1] >= 0
     then highest2 := delta2[1];

     for loop := 1 to upper
     do begin
	if (d[loop + 1] >= 0) and (d[loop + 1] > highest1)
	then highest1 := d[loop + 1];

	if (delta2[loop + 1] >= 0) and (delta2[loop + 1] > highest2)
	then highest2 := delta2[loop + 1];
	end;
end; {search}

procedure search2;

{this procedure searches for the highest number from highest1 & highest2
 from procedure search}

begin
     if highest1 > highest2
     then ks := highest1
     else ks := highest2;
end; {search2}

procedure search1(var expected: rarray);

{this procedure searches for the expected frequency less than 5 and
 combines to the next expected frequency; which is also done to the
 corresponding observed frequency}

var
   loop1,loop2: integer;

begin
     count := 0;
     loop1 := 1;
     loop2 := 1;

     while (loop2 < m)
     do begin
	while (expected[loop1] < 5) and (loop2 < m)
	do begin
	   expected[loop1] := expected[loop1] + expected[loop2 + 1];
	   n[loop1] := n[loop1] + n[loop2 + 1];
	   loop2 := loop2 + 1;
	   end;

	if loop2 < m
	then begin
	     loop1 := loop1 + 1;
	     loop2 := loop2 + 1;
	     expected[loop1] := expected[loop2];
	     n[loop1] := n[loop2];
	     end;
	end;

     {checks if the last expected frequency is less than 5, if so the
      expected is added to the penultimate expected frequency}

     if (expected[loop1] < 5)
     then begin
	  expected[loop1 - 1] := expected[loop1 - 1] + expected[loop1];
	  n[loop1 - 1] := n[loop1 - 1] + n[loop1];
	  loop1 := loop1 - 1;
     end;

     count := loop1
end; {search1}

procedure chiq;

{calculates (O-E)/E for corresponding observed & expected frequencies
 read in by procedure read_data1}

var
   loop: integer;

begin
     ch := 0;

     for loop := 1 to count
     do begin
	d[loop] := (n[loop] - expf[loop]);
	d[loop] := (sqr(d[loop])) / expf[loop];
	ch := ch + d[loop];
	end;
end; {chiq}

procedure out;

{outputs the results of the chiq procedure}

var
   loop: integer;

begin
     writeln(results);
     writeln(results,'  Observed       Expected');
     writeln(results,'Frequency (O)  Frequency (E)   (O-E)**2/E');
     writeln(results);

     for loop := 1 to count
     do begin
	write(results,trunc(n[loop]):8,'       ',expf[loop]:8:2);
	write(results,'        ',d[loop]:8:2);
	writeln(results);
	end;

     writeln(results);
     writeln(results,'Calculated Chi-Squared value = ',ch:4:2);
end; {out}

procedure table;

var
   loop1,upper: integer;

begin
     writeln(results,'Fo     Fn  ');
     writeln(results,'       i/n     i/n - Fo       Fo - {i/n-1}  ');
     writeln(results,'***********************************************');
     writeln(results);

     if (a = 1) or (a = 2)
     then upper := sample_size
     else upper := m;

     for loop1 := 1 to upper
     do begin
	write(results,expf[loop1]:3:2,'     ',n[loop1]:3:2);
	write(results,'    ',d[loop1]:3:2,'      ',delta2[loop1]:3:2);
	writeln(results);
	end;

     writeln(results,'Calculated value of D = ',ks:6:4);
end; {table}

procedure mixfit;

{this procedure chooses the test statistic and the goodness-of-fit test}

var
   loop1,loop2: integer; {needed because the compiler will not allow
			  non-local variables i and j to be loop control
			  variables}

begin
     case a of
     1:   begin
	  test := 1;
{         writeln(results,'X Values    Cumulative Probabilities');
	  writeln(results);}

	  for loop1 := 1 to sample_size
	  do begin
	     i := loop1;

	     for loop2 := 1 to k
	     do begin
		j := loop2;
		alnorm(d,true);
		answer1(sum);
		end;

	     expf[i] := sum;
	     sum := 0;
{            writeln(results,u[i]:6:2,'      ',expf[i]:14:2);}
	     end;
{         writeln(results)}
	  end;

     2:   begin
	  test := 1;
{         writeln(results,'X Values    Cumulative Probabilities');
	  writeln(results);}

	  for loop1 := 1 to sample_size
	  do begin
	     i := loop1;

	     for loop2 := 1 to k
	     do begin
		j := loop2;
		exponential(sum);
		end;

	     expf[i] := sum;
	     sum := 0;
{            writeln(results,u[i]:6:2,'      ',expf[i]:14:2);}
	     end;
{         writeln(results)}
	  end;

     3:   begin
	  test := 1;
	  writeln(results,'X Values    Cumulative Probabilities');
	  writeln(results);

	  for loop1 := 1 to m
	  do begin
	     i := loop1;

	     for loop2 := 1 to k
	     do begin
		j := loop2;
		poiss1(d);
		end;

	     expf[i] := sum;
	     writeln(results,u[i]:6:2,'      ',expf[i]:14:2);
	     end;
	  writeln(results)
	  end;

     4:   begin
	  test := 1;
	  writeln(results,'X Values      Cumulative Probabilities');
	  writeln(results);

	  for loop1 := 1 to m
	  do begin
	     i := loop1;
	     count := (trials - round(u[i]));
	     binomial(res);

	     for loop2 := 1 to k
	     do begin
		j := loop2;
		binomial1(d,mu);
		end;

	     expf[i] := sum;
	     writeln(results,u[i]:6:2,'      ',expf[i]:14:2);
	     end;
	  writeln(results)
	  end;

     5:   begin
	  test := 2;
	  writeln(results,'Upper X Values    Expected Frequencies');
	  writeln(results);

	  for loop1 := 1 to m
	  do begin
	     i := loop1;

	     for loop2 := 1 to k
	     do begin
		j := loop2;
		alnorm(d,true);
		answer1(sum);
		end;

	     expf[i] := sum;
	     sum := 0;
	     end;
	  pdf;
	  delist(expf,expected1);
	  writeln(results)
	  end;

     6:   begin
	  test := 2;
	  writeln(results,'Upper X Values    Expected Frequencies');
	  writeln(results);

	  for loop1 := 1 to m
	  do begin
	     i := loop1;

	     for loop2 := 1 to k
	     do begin
		j := loop2;
		exponential(sum);
		end;

	     expf[i] := sum;
	     sum := 0;
	     end;
	  pdf;
	  delist(expf,expected1);
	  writeln(results)
	  end;

     7:   begin
	  test := 2;
	  writeln(results,'X Values    Expected Frequencies');
	  writeln(results);

	  for loop1 := 1 to m
	  do begin
	     i := loop1;

	     for loop2 := 1 to k
	     do begin
		j := loop2;
		poiss1(d);
		end;

	     expf[i] := sum;
	     sum := 0;
	     end;
	  delist(expf,expf);
	  writeln(results)
	  end;

     8:   begin
	  test := 2;
	  writeln(results,'U Values    Expected Frequencies');
	  writeln(results);

	  for loop1 := 1 to m
	  do begin
	     i := loop1;
	     count := (trials - round(u[i]));
	     binomial(res);

	     for loop2 := 1 to k
	     do begin
		j := loop2;
		binomial1(d,mu);
		end;

	     expf[i] := sum;
	     sum := 0;
	     end;
	  delist(expf,expf);
	  writeln(results)
	  end;
     end;

     case test of
     1:   begin
	  writeln(results,'Kolmogorov-Smirnov Test');
	  writeln(results,'***********************');
	  writeln(results);
	  prob_data;
	  test1;
	  highest1 := 0;
	  highest2 := 0;
	  search;
	  search2;
	  writeln(results,'Calculated value of D = ',ks:6:4);
	  end;

     2:   begin
	  writeln(results,'Chi-Squared Test');
	  writeln(results,'****************');
	  writeln(results);
	  search1(expf);
	  chiq;
	  out;
	  end;
     end;
end; {mixfit}

begin {goodness}
     {set initial values of i and j to satisfy the compiler}
     i := 0;
     j := 0;

     {call the procedures to read and check the data}
     ifault := 0;
     read_data1;
     read_data2;
     check;
     check1;

     if (ifault = 1)
     then begin
	  writeln('Data incorrect, please check errors printed in your output file');
	  end
     else begin
	  {produce the results}
	  sum := 0;
	  mixfit
	  end;

end; {goodness}

begin {main program - mixgood}
     {if fullyclose is true, all buffers are flushed when a file is closed}
     fullyclose := true;

     {open the input data and output results files}
     writeln('Enter mix input file name');
     readln(filename);
     assign(data,filename);

     writeln('Enter output file name');
     readln(filename1);
     assign(results,filename1);

     {initialise file}
     reset(data);
     rewrite(results);

     mix;
     goodness;

     close(data,fullyclose);
     close(results,fullyclose)

end. {mixgood}
