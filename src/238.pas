PROCEDURE karstor(n:nonnegint;VAR x,y,d:realarray;VAR alpha,beta,sad:real;
VAR ifault:nonnegint);

{Algorithm AS 238  Appl. Statist. (1988) Vol.37, No.3}

{Pascal translation of the original bubblesort variant of algorithm AS74
 for the L1 norm fitting of a straight line with real storage of integers
 Applied Statistics (1974),Vol 23,No 2}

CONST big=1000000;
     prec=0.000001;
VAR  a,b,i,j,it,p,q:nonnegint;
     sum,xa,ya,xi,yi:real;
     stop:boolean;
     t:ARRAY[nonnegint] of real; {pairwise slopes}

BEGIN
i:=1;xa:=x[1];
stop:=true;
WHILE (i<n) AND stop DO
   BEGIN
   i:=i+1;stop:=abs(x[i]-xa)<prec
   END;
IF stop THEN ifault:=1
ELSE
   BEGIN
   a:=1;b:=0;
   it:=0;ifault:=2;
   REPEAT
      it:=it+1;

      {calculate slopes}
      sum:=0.0;
      FOR i:=1 TO n DO
         BEGIN
         xa:=x[a];ya:=y[a];
         xi:=x[i]-xa;yi:=y[i]-ya;
         d[i]:=i;
         IF abs(xi)<prec THEN t[i]:=big
         ELSE t[i]:=yi/xi;
         sum:=sum+abs(xi)
         END;
      sum:=0.5*sum;

      {order slopes}
      FOR j:=n-1 DOWNTO 1 DO
      FOR i:=1 TO j DO
         BEGIN
         p:=trunc(d[i]+0.5);q:=trunc(d[i+1]+0.5);
         IF t[p]>t[q] THEN
            BEGIN
            d[i]:=q;d[i+1]:=p
            END
         END;

      {select best line}
      i:=0;xi:=0.0;
      REPEAT
         i:=i+1;p:=trunc(d[i]+0.5);
         xi:=xi+abs(x[p]-xa);
         IF xi>=sum THEN i:=n
      UNTIL i=n;

      IF p=b THEN
         BEGIN
         ifault:=0;it:=n
         END
      ELSE
         BEGIN
         b:=a;a:=p
         END
   UNTIL it=n;     

   {L1 estimates,norm and residuals}
   beta:=t[p];alpha:=y[a]-beta*x[a];
   sum:=0.0;
   FOR i:= 1 TO n DO
      BEGIN
      yi:=y[i]-alpha-beta*x[i];d[i]:=yi;
      sum:=sum+abs(yi)
      END;
   sad:=sum
   END
END;{of karstor}

   
PROCEDURE karst(n:nonnegint;VAR x,y,e:realarray;VAR alpha,beta,sad:real;
VAR ifault:nonnegint);

{Algorithm AS 238  Appl. Statist. (1988) Vol.37, No.3}

{partial insertionsort variant of algorithm AS74
for the L1 norm fitting of a straight line with a boolean trap
Applied Statistics (1974),Vol 23,No 2}

CONST big=10000000;
     prec=0.000001;
VAR  a,c,i,j,k:nonnegint;
     hold,sum,sum2,xa,ya:real;
     jump,stop:boolean;
     d:ARRAY[nonnegint] of nonnegint; {permuted indices}
     t:ARRAY[nonnegint] of real;      {pairwise slopes}
   old:ARRAY[nonnegint] of boolean;   {previous values of a}

BEGIN
i:=1;xa:=x[1];
stop:=true;
WHILE (i<n) AND stop DO
   BEGIN
   i:=i+1;stop:=abs(x[i]-xa)<prec
   END;
IF stop THEN ifault:=1
ELSE
   BEGIN
   a:=1;stop:=false;
   ifault:=0;
   FOR i:=2 TO n DO old[i]:=false;

   REPEAT

      {calculate slopes}      
      xa:=x[a];ya:=y[a];
      old[a]:=true;sum:=0.0;
      FOR i:=1 TO n DO
         BEGIN
         hold:=x[i]-xa;
         IF abs(hold)<prec
         THEN t[i]:=big
         ELSE t[i]:=(y[i]-ya)/hold;
         sum:=sum+abs(hold)
         END;
      sum2:=0.5*sum;

      {locate weighted median}
      k:=1;d[1]:=1;
      sum:=abs(x[1]-xa);
      FOR i:=2 TO n DO
         BEGIN k:=k+1;j:=k;
         sum:=sum+abs(x[i]-xa);hold:=t[i];
         REPEAT
            j:=j-1;c:=d[j];
            IF hold<t[c] THEN d[j+1]:=c
            ELSE
               BEGIN
               d[j+1]:=i;j:=0
               END
         UNTIL j<2;
         IF j=1 THEN d[1]:=i;
  
         jump:=false;
         REPEAT
            IF sum>sum2 THEN
               BEGIN
               hold:=sum-abs(x[d[k]]-xa);
               IF hold>sum2 THEN
                  BEGIN
                  sum:=hold;k:=k-1
                  END         
               ELSE jump:=true
               END
            ELSE jump:=true

         UNTIL jump                   
         END;

      c:=d[k];
      IF old[c] THEN stop:=true ELSE a:=c
   UNTIL stop;     

   {L1 estimates,norm and residuals}
   beta:=t[c];alpha:=y[a]-beta*x[a];
   sum:=0.0;
   FOR i:= 1 TO n DO
      BEGIN
      hold:=y[i]-alpha-beta*x[i];e[i]:=hold;
      sum:=sum+abs(hold)
      END;
   sad:=sum
   END
END;{of karst}


PROCEDURE simbarrob(n:nonnegint;VAR x,y,d:realarray;VAR alpha,beta,sad:real;
VAR ifault:nonnegint);

{Algorithm AS 238  Appl. Statist. (1988) Vol.37, No.3}

{Pascal translation of algorithm AS132
 for the L1 norm fitting of a straight line
 Applied Statistics (1978),Vol 27,No 3}

LABEL 160;
CONST big=10000000;
      acu=0.000001;
VAR  a1,a2,aaa,aaaa,ahalf,aone,bbb,bbbb,ddd,det,hold,ratio,sum,t,test,
     tot1,tot2,y1,y2:real;
     i,j,ibas1,ibas2,iin,iout,isave,iter:nonnegint;
     rho,zzz:-2..2;
     flag,more,stop:boolean;
     inext:ARRAY[nonnegint] of posnegint;

FUNCTION isign(i:posnegint;k:posnegint):posnegint;
{Pascal implementation of Fortran isign function}
BEGIN
IF k>=0 THEN isign:=i ELSE isign:=-i
END;{of isign}

BEGIN

{determine initial basis}
a1:=x[1];y1:=y[1];
i:=1;stop:=true;
WHILE (i<n) AND stop DO
   BEGIN
   i:=i+1;stop:=abs(a1-x[i])<acu
   END;
IF stop THEN ifault:=1
ELSE
   BEGIN
   ifault:=0;iter:=0;
   ahalf:=0.5+acu;aone:=ahalf+ahalf;
   a2:=x[i];y2:=y[i];
   ibas1:=1;ibas2:=i;
   det:=1.0/(a2-a1);
 
   {initial values of a,b,d}
   aaaa:=(a2*y1-a1*y2)*det;bbbb:=(y2-y1)*det;
   FOR i:=2 TO n DO
   IF y[i]>=(aaaa+bbbb*x[i]) THEN d[i]:=1.0 ELSE d[i]:=-1.0;
   d[ibas1]:= 0.0;d[ibas2]:=-1.0;
   tot1:=1.0;tot2:=x[ibas2];
   FOR i:=2 TO n DO
   BEGIN
      tot1:=tot1+d[i];tot2:=tot2+d[i]*x[i];
   END;
   t:=(a2*tot1-tot2)*det;flag:=abs(t)>aone; 
   
   {main iterative loop begins}
   REPEAT
   flag:= NOT flag;
   IF flag THEN
      BEGIN
      t:=(tot2-a1*tot1)*det;
      IF abs(t) < aone THEN GOTO 160;
      iout:=ibas2;x[iout]:=a1;
      aaa:=a1;bbb:=a2;
      END
   ELSE
      BEGIN
      IF iter=0 THEN det:=-det
      ELSE
         BEGIN
         t:=(tot2-a2*tot1)*det;
         IF abs(t) < aone THEN GOTO 160;
         END;
      aaa:=a2;bbb:=a1; 
      iout:=ibas1;
      END;{not flag}
   IF t>0.0 THEN rho:=1
   ELSE
   BEGIN
      rho:=-1;t:=-t;
      det:=-det 
   END;
   t:=0.5*t;

   {perform partial sort of ratios}
   ratio:=big;sum:=ahalf;
   inext[ibas1]:=ibas2;
   FOR i:=1 TO n DO
   BEGIN
      ddd:=(x[i]-aaa)*det;
      IF ddd*d[i]>acu THEN
      BEGIN
         test:=(y[i]-(aaaa+bbbb*x[i]))/ddd;
         IF test<ratio THEN
            BEGIN
            j:=ibas1;sum:=sum+abs(ddd);
            more:=true;stop:=false;
               REPEAT
               isave:=abs(inext[j]);
               IF test >= d[isave] THEN
                  BEGIN
                  more:=false;stop:=true
                  END
               ELSE
                  BEGIN
                  IF sum<t THEN stop:=true
                  ELSE
                     BEGIN
                     hold:=sum-abs((x[isave]-aaa)*det);
                     IF hold<t THEN stop:=true
                     ELSE
                        BEGIN
                        d[isave]:=isign(1,inext[j]);sum:=hold;
                        inext[j]:=inext[isave]
                        END{hold>=t}
                     END{sum>=t}
                  END;{test<d{isave]}
               UNTIL stop;

            IF more THEN
               REPEAT
               j:=isave;isave:=abs(inext[j]);
               UNTIL test>=d[isave];

            inext[i]:=inext[j];
            IF d[i]>=0.0 THEN inext[j]:=i ELSE inext[j]:=-i;
            d[i]:=test;
            IF sum>=t THEN
               BEGIN
               iin:=abs(inext[ibas1]);ratio:=d[iin]
               END
            END;{test<ratio}
        END;{ddd*d[i]>acu}
      END;{i loop}

   {update basis indicaTOrs}
   iin:=abs(inext[ibas1]);j:=iin;
   stop:=false;
      REPEAT
      isave:=abs(inext[j]);
      IF isave=ibas2 THEN stop:=true
      ELSE
         BEGIN
         zzz:=isign(1,inext[j]);d[isave]:=-zzz;
         zzz:=zzz+zzz;j:=isave;
         tot1:=tot1-zzz;tot2:=tot2-zzz*x[isave]
         END;
      UNTIL stop;

   zzz:=isign(1,inext[ibas1]);d[iout]:=-rho;
   tot1:=tot1-rho-zzz;tot2:=tot2-rho*bbb-zzz*x[iin];
   iter:=iter+1;
   IF flag THEN
      BEGIN
      x[ibas2]:=a2;ibas2:=iin;
      a2:=x[iin];y2:=y[iin];
      det:=1.0/(a1-a2);d[ibas2]:=-1.0;
      aaaa:=(a1*y2-a2*y1)*det;bbbb:=(y1-y2)*det
      END
   ELSE
      BEGIN
      ibas1:=iin;ibas1:=iin;
      a1:=x[iin];y1:=y[iin];
      det:=1.0/(a2-a1);d[ibas1]:=0.0;
      aaaa:=(a2*y1-a1*y2)*det;bbbb:=(y2-y1)*det
   END;
   UNTIL false;

   {calculate estimates,deviations,sad}
160:sum:=0.0;
   FOR i:=1 TO n DO
      BEGIN
      ddd:=y[i]-(aaaa+bbbb*x[i]);d[i]:=ddd;
      sum:=sum+abs(ddd)
   END;
   alpha:=aaaa;beta :=bbbb;
   sad:=sum;

   END
END;{of simbarrob}

