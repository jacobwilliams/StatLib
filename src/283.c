// The file contains the following:

// fishfun.c  /* C-source of Fisher's symmetry test function  + driver program */

// groupfun.c /* C source of Pitman's 2 sample test function + driver program */

// fishint.c /* exact (integer) version of fishfun.c */

// groupint.c /* exact (integer) version of groupfun.c */

// psydata /* data file of Rutter score data for fishfun */

// agedata /* data file of age data for groupfun */

// iagedat /* integer version of age data */

/* -----------------------------------------------------------------------*/
/* fishfun.c */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define MAXNUM 200 /* max. sample size */
float fisherprob(float *,int,float,float *, long *);

#define NAMLEN 40 /* max. length of file-name */
int compare(int *x,int *y)
	{
	return(*x-*y);
	}
main()
	{
	FILE *fp;
	char filename[NAMLEN]; /* array to hold filename */
	float a[MAXNUM],cv,p,ptol=.01; /* fractional error allowed on prob. */
	int scanval,n;
	long binsreq;

	puts("\nEnter input file name: ");
	scanf("%40s",filename); /* get file name (MAX 40 CHARS) */

	if ((fp = fopen(filename,"r")) == 0 ) /* open file */
		{
		puts ("Can't open input file"); /* unsuccessful */
		return;
		}

	n=0;
	while((scanval=fscanf(fp,"%f",&a[n])) != EOF && scanval != NULL )
		{
		printf("number %d is %f\n",n,a[n]);
		++n;
		}
	close(fp);
	printf("read in %d numbers\n",n);

	p=fisherprob(a,n,ptol,&cv,&binsreq);
	if(p<0.)
		{
		printf("estimating probability to a CV of %f requires %ld bins-sorry",ptol,binsreq);
		return;
		}
	printf("calculation used %ld histogram bins\n",binsreq);
	printf("1-tailed probability is %f\n",p);
	printf("Coefficient of variation on estimate of p is %f%%\n",cv); 
	/* can repeat with ptol less if cv is too large */
	}

float fisherprob(float a[],int n,float ptol,float *cv,long *binsreq)
/* This function returns the 1-tailed probability from
Fisher's symmetry test, the permutation, i.e. distribution-free
analogue of the paired t-test. */

	{
	int numberpos=0,numberneg=0, *positive, *x,smallnum,topbin,
	m, bin,i,ipossum=0,inegsum=0,s,n0,n0dash,ndash,ntop;
	float possum=0,negsum=0,scale,sum,p,factor,mu,sigsq,term1,term2,tail, *b;

	*binsreq=n;
	x=calloc(n,sizeof(int));
	positive=calloc(n,sizeof(int)); /* allocate work arrays */
	if(!x || !positive)return(-1.); /* and return -ve p-value if not possible */

	for(i=0,mu=0.,sigsq=0.;i<n;i++){
		sigsq+=a[i]*a[i]; /* accumulate permutation variance */
		if(a[i]>=0)
			{
			++numberpos; /* count number of +ve values */
			possum+=a[i]; /* and accumulate their sum */
			positive[i]=1; /* flag value as +ve */
			}
		else
				{
			++numberneg;
			negsum-=a[i];
			a[i]=-a[i];
			positive[i]=0;
			}
		}

	mu=(possum+negsum)/2.;
	sigsq/=4.; /* permutation mean and variance */
	sum=(possum<negsum)?possum:negsum; /* choose smaller sum as test statistic */
	m=(possum<negsum)?numberpos:numberneg;
	factor=sqrt((n-(mu-sum)*(mu-sum)/sigsq)/48.); /*for error calculation */
	term1=mu*(mu-sum)/(ptol*sigsq);
	term2=2./ptol;
	bin=factor*((term1>term2)?term1:term2);
	/* no. of bins needed for required error on p*/

	smallnum=(numberpos<numberneg)?numberpos:numberneg; 
	/* number in smaller group */
	if(smallnum>0)
		scale=bin/sum; /* scale sample values so that smaller
		                                                                                            sum is in bin number bin */
	else
		scale=1.; /* if only one sign present calculation is trivial. However,
		                                                                                            just go through the usual motions */
	for (i=0;i<n;++i){
		x[i]=a[i]*scale+0.5; /* rescale sample values. Adding 0.5 rounds to nearest integer rather than truncating*/
		if(positive[i])ipossum+=x[i];
		else
			inegsum+=x[i];
		}
	bin=(ipossum<inegsum)?ipossum:inegsum; /* rescale number of bins */

	*binsreq=bin+2;
	b=calloc(bin+2,sizeof(float)); /* allocate work array */
	if(!b)
		return(-1.); /* and return if not enough memory */

	qsort(x,n,sizeof(int),compare); /* sort sample values */
	for(i=0;i<n;++i)
		x[i]=(x[i]<bin+1)?x[i]:bin+1; /* any sample values
		                              greater than bin are put into the overflow bin */
	/*      for(i=0;i<=bin+1;i++)b[i]=0.0; */
	++b[0]; /* start histogram off - entries at zero and x0 */
	++b[x[0]];
	n0=x[0];
	for(s=1;s<n;s++)
		{
		n0dash=n0+x[s];
		if(n0dash>bin)
			{
			b[bin+1]*=2.; /* just double entries in overflow bin */
			for(ndash=bin;ndash>=bin-x[s]+1;ndash--)b[bin+1]+=b[ndash];
			}/* translate histogram */
		ntop=(bin>n0dash)?n0dash:bin;
		for(ndash=ntop;ndash>=x[s];ndash--)b[ndash]+=b[ndash-x[s]];
		n0=n0dash;
		}

	tail=0.; /* count entries in tail */
	for(i=0;i<=bin;++i)
		tail+=b[i];
	p=tail/(tail+b[bin+1]);
	*cv=100*b[bin]*factor/tail;
	free(b); /* free memory used for work arrays */
	free(x);
	free(positive);
	return(p);


	}

/* -----------------------------------------------------------------------*/
/* groupfun.c */

#include <stdio.h> /* for reading data file */
#include <math.h> /* for square root */
#include <stdlib.h>

#define MAXNUM 160 /* max. sample size */

#define LEN 40 /* max. length of filename */
int compare(int *x,int *y)
	{
	return(*x-*y);
	}
main()
	{
	FILE *fp;
	int n,firstgroup[MAXNUM],scanval;
	float a[MAXNUM],p,cv,ptol=.05; /* fractional error required on prob. */
	long binsreq;
	float pitmanprob(float *,int *, int, float, float *, long *);
	char filename[LEN]; /* array to hold filename */

	puts("\nEnter input file name: ");
	scanf("%40s",filename); /* get file name (MAX 40 CHARS) */

	if ((fp = fopen(filename,"r")) == 0 ) /* open file */
		{
		puts ("Can't open input file"); /* unsuccessful */
		return;
		}

	n=0;
	while((scanval=fscanf(fp,"%f%d",&a[n],&firstgroup[n])) != EOF && scanval != NULL)
		{
		printf("number %d is %f group %d\n",n,a[n],firstgroup[n]);
		++n;
		if(n>=MAXNUM)
			{
			printf("sample size too large. Maximum is %d\n",MAXNUM);
			return;
			}
		}
	close(fp);
	printf("read in %d numbers\n",n);

	p=pitmanprob(a,firstgroup,n,ptol,&cv,&binsreq); /* find tail probability */
	if(p<0.){
		printf("%ld array elements requested, which is too many.\n",binsreq);
		return;
		}
	else 	{
		printf("%ld array elements used\n",binsreq);
		printf("1-tailed probability is %f\n",p);
		printf("CV of estimated probability is %f %%",cv);
		}
	}

float pitmanprob(float a[],int firstgroup[],int n,float ptol,float *cv,long *binsreq)
/* returns 1-tail probability from Pitman's 2-sample permutation test,
the distribution-free analogue of the grouped t-test */

	{
	int numberone=0,numbertwo=0,
	smallnum,k,j,smallgroup,ntop,kbot,topbin,
	bin,i,isum=0,s,n0,n0dash,ndash,m, *dope, *x;
	float onesum=0,twosum=0,scale,sum,p,tail,term1,term2,
	minnum,maxnum,onemean,twomean,mu,sigsq,em,en,factor, *b;

	*binsreq=n;
	dope=calloc(n,sizeof(int));
	x=calloc(n,sizeof(int)); /* allocate work area */
	if(!dope || !x)return(-1.); /* return -ve p value if not possible */

	for(i=0;i<n;i++){
		if(i == 0){
			minnum=a[i]; /* find max. and min. sample values */
			maxnum=a[i];
			}
		else 	{
			minnum=(a[i]>minnum)?minnum:a[i];
			maxnum=(a[i]>maxnum)?a[i]:maxnum;
			}
		if(firstgroup[i] ==1)
			{
			++numberone; /* count numbers of in each group */
			onesum+=a[i]; /* and their sums */
			}
		else
				{
			++numbertwo;
			twosum+=a[i];
			}
		}

	onemean=onesum/numberone; /* find means for each group */
	twomean=twosum/numbertwo;
	smallgroup=(onemean>twomean)?2:1;
	for(i=0;i<n;i++)a[i]-=minnum; /* translate sample values to be >= 0 */
	smallnum=(numberone<numbertwo)?1:2;
	onesum-=numberone*minnum; /* and modify sums and maximum value */
	twosum-=numbertwo*minnum;
	maxnum-=minnum;

	if(smallgroup != smallnum)
		{ /* if necessary, reflect values so smaller group has smaller mean */
		for(i=0;i<n;i++)a[i]=maxnum-a[i];
		onesum=numberone*maxnum-onesum;
		twosum=numbertwo*maxnum-twosum;
		}
	sum=(smallnum == 1)?onesum:twosum;
	m=(numberone<numbertwo)?numberone:numbertwo; /* number in smaller group */

	em=m;
	en=n;
	for(mu=0.,sigsq=0.,i=0;i<n;i++){
		mu+=a[i];
		sigsq+=a[i]*a[i];
		for(j=i+1;j<n;j++)sigsq=sigsq-2.*a[i]*a[j]/(en-1.);
		}
	mu=(em/en)*mu;
	sigsq=(em/en)*(1.-(em/en))*sigsq; /* permutation mean and variance */
	factor=sqrt((em/en)*(1.-(em/en))*(en+(mu-sum)*(mu-sum)*(en+1.)/(en*sigsq))/12.);
	term1=mu*(mu-sum)/(ptol*sigsq);
	term2=2./ptol;
	bin=factor*((term1>term2)?term1:term2); /* number of histogram bins needed */

	scale=bin/sum; /* scale sample values so that sum of smaller sample
		                                                    is in bin number bin */
	for (i=0;i<n;i++){
		x[i]=a[i]*scale+0.5; /* rescaling sample values.
				                        Adding 0.5 rounds to nearest bin rather than truncating */
		if(firstgroup[i]==smallnum)isum+=x[i];
		}
	bin=isum; /* find exact histogram bin for start of reject H0 region */


	*binsreq=(long)(bin+2)*m; 
	b=calloc(*binsreq,sizeof(float)); /* allocate and zero workspace */
	if(!b) 
		return(-1.); /* not enough workspace */
	qsort(x,n,sizeof(int),compare); /* sort sample */
	for(i=0;i<n;i++)
		x[i]=(x[i]<bin+1)?x[i]:bin+1; /* any sample values greater than
		            the test statistic "bin" will cause b[>bin] to be occupied and are set
		            to bin+1, the overflow bin NOW rather than later so as not to cause integer
		            overflow */

	for(i=0;i<m;i++)dope[i]=i*(bin+2); /* offsets to mimick 2 dim. array */
	++b[dope[0]+x[0]]; /* start histogram off */
	/* for 2 dim. array code would be  ++b[0][x[0]];    */

	n0=x[0];
	for(s=1;s<n;s++)
		{ 
		kbot=(m+s-n>1)?m+s-n:1; 
		n0dash=n0+x[s]; /* translate each level by x[s] */
		for(k=(s>(m-1))?(m-1):s;k>=kbot;k--) /* only build up histograms
						                                    as far as mth level, which is the one we want */
			{
			if(n0dash>bin)for(ndash=bin-x[s]+1;ndash<=bin+1;ndash++)
				/* translate k-1 th level onto kth level...these terms go into overflow bin */
				b[dope[k]+bin+1]+=b[dope[k-1]+ndash]; 
			/* for 2 dim. array...    b[k][bin+1]+=b[k-1][ndash]; */
			ntop=(bin>n0dash)?n0dash:bin;
			for(ndash=x[s];ndash<=ntop;ndash++){
				/* these terms just translate along, but may go into overflow bin */
				b[dope[k]+ndash]+=b[dope[k-1]+ndash-x[s]];
				/* for 2 dim. array...  b[k][ndash]+=b[k-1][ndash-x[s]]; */
				}
			}

		++b[dope[0]+x[s]]; /* add new term to lowest level */
		/* for 2 dim. array...  ++b[0][x[s]]; */
		n0=n0dash;
		}

	tail=0.; /* find tail probability */
	for(i=0;i<=bin;i++)
		/* for 2 dim. array...   tail+=b[m-1][i]; */
		tail+=b[dope[m-1]+i];

	*cv=100.*b[dope[m-1]+bin]*factor/tail; /* find coefficient of variation */
	/* for 2d array...   *cv=100.*b[m-1][bin]*factor/tail; */
	free(x); /* free off workspace */
	free(dope);
	free(b);
	return(tail/(tail+b[dope[m-1]+bin+1]));
	}

/* -----------------------------------------------------------------------*/
/* fishint.c */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define MAXNUM 200 /* max. sample size */
float fishprob();

#define NAMLEN 40 /* max. length of file-name */
int comp(int *x,int *y)
{
return(*x-*y);
}
main()
      {
      FILE *fp;
      char filename[NAMLEN]; /* array to hold filename */
      float cv,p;
      int scanval,n,a[MAXNUM];
      long binsreq;

      puts("\nEnter input file name: ");
      scanf("%40s",filename); /* get file name (MAX 40 CHARS) */

      if ((fp = fopen(filename,"r")) == 0 ) /* open file */
            {
            puts ("Can't open input file"); /* unsuccessful */
            return;
            }

      n=0;
      while((scanval=fscanf(fp,"%d",&a[n])) != EOF && scanval != NULL )
            {
            printf("number %d is %d\n",n,a[n]);
            ++n;
            }
      close(fp);
      printf("read in %d numbers\n",n);

      p=fishprob(a,n,&binsreq);
      if(p<0.)
            {
            printf("calculating probability requires %ld bins-sorry",binsreq);
            return;
            }
      printf("calculation used %ld histogram bins\n",binsreq);
      printf("1-tailed probability is %f\n",p);
      }

float fishprob(int a[],int n,long *binsreq)
/* This function returns the 1-tailed probability from
Fisher's symmetry test, the permutation, i.e. distribution-free
analogue of the paired t-test. */

      {
      long bin,i,s,n0,n0dash,ndash,ntop;
      long possum=0,negsum=0,sum;
      float p,tail;
float *b;

      for(i=0;i<n;i++){
            if(a[i]>=0)
                  {
                  possum+=a[i]; /* and accumulate their sum */
                  }
            else
                        {
                  negsum-=a[i];
                  a[i]=-a[i];
                  }
            }

      sum=(possum<negsum)?possum:negsum; /* choose smaller sum as test statistic */
      bin=sum;
      /* no. of bins needed */
b=calloc(bin+2,sizeof(float));
      *binsreq=bin+2;
      if(!b)
            return(-1.);

      qsort(a,n,sizeof(int),comp);
      for(i=0;i<n;++i)
            a[i]=(a[i]<bin+1)?a[i]:bin+1; /* any sample values
                              greater than bin are put into the overflow bin */
      for(i=0;i<=bin+1;i++)b[i]=0.0;
      ++b[0]; /* start histogram off - entries at zero and x0 */
      ++b[a[0]];
      n0=a[0];
      for(s=1;s<n;s++)
            {
            n0dash=n0+a[s];
            if(n0dash>bin)
                  {
                  b[bin+1]=2*b[bin+1]; /* just double entries in overflow bin */
                  for(ndash=bin;ndash>=bin-a[s]+1;ndash--)b[bin+1]+=b[ndash];
                  }/* translate histogram */
            ntop=(bin>n0dash)?n0dash:bin;
            for(ndash=ntop;ndash>=a[s];ndash--)b[ndash]+=b[ndash-a[s]];
            n0=n0dash;
            }

      tail=0.; /* count entries in tail */
      for(i=0;i<=bin;++i)
            tail+=b[i];
      p=tail/(tail+b[bin+1]);
      free(b); /* free work space array */
      return(p);


      }

/* -----------------------------------------------------------------------*/
/* groupint.c */

#include <stdio.h> /* for reading data file */
#include <math.h> /* for square root */
#include <stdlib.h>

#define MAXNUM 200 /* max. sample size */

#define LEN 40 /* max. length of filename */
int comp(int *x,int *y)
{
return(*x-*y);
}
void main()
      {
      FILE *fp;
      int n,firstgroup[MAXNUM],scanval,x[MAXNUM];
      long binsreq;
      float p;
      float grouprob();
      char filename[LEN]; /* array to hold filename */

      puts("\nEnter input file name: ");
      scanf("%40s",filename); /* get file name (MAX 40 CHARS) */

      if ((fp = fopen(filename,"r")) == 0 ) /* open file */
            {
            puts ("Can't open input file"); /* unsuccessful */
            return;
            }

      n=0;
      while((scanval=fscanf(fp,"%d%d",&x[n],&firstgroup[n])) != EOF && scanval != NULL)
            {
            printf("number %d is %d group %d\n",n,x[n],firstgroup[n]);
            ++n;
            if(n>=MAXNUM)
                  {
                  printf("sample size too large. Maximum is %d\n",MAXNUM);
                  return;
                  }
            }
      close(fp);
      printf("read in %d numbers\n",n);

      p=grouprob(firstgroup,x,n,&binsreq); /* find tail probability */
      if(p<0.){
            printf("%ld array elements requested, which is too many.\n",binsreq);
            return;
            }
      else  {
            printf("%ld array elements used\n",binsreq);
            printf("1-tailed probability is %f\n",p);
            }
      }

float grouprob(int firstgroup[],int x[],int n,long *binsreq)
/* returns 1-tail probability from Pitman's 2-sample permutation test. */

      {
      int numberone=0,numbertwo=0,
      smallnum,k,j,smallgroup,ntop,kbot,topbin,
      bin,i,s,n0,n0dash,ndash,m,onesum=0,twosum=0,sum,minnum,maxnum;
      float p,tail,onemean,twomean;
int *dope;
float *b;
*binsreq=n;
dope=calloc(n,sizeof(int));
if(!dope)return(-1.);
      for(i=0;i<n;i++){
            if(i == 0){
                  minnum=x[i]; /* find max. and min. sample values */
                  maxnum=x[i];
                  }
            else  {
                  minnum=(x[i]>minnum)?minnum:x[i];
                  maxnum=(x[i]>maxnum)?x[i]:maxnum;
                  }
            if(firstgroup[i] ==1)
                  {
                  ++numberone; /* count numbers of in each group */
                  onesum+=x[i]; /* and their sums */
                  }
            else
                        {
                  ++numbertwo;
                  twosum+=x[i];
                  }
            }

      onemean=onesum/numberone; /* find means for each group */
      twomean=twosum/numbertwo;
      smallgroup=(onemean>twomean)?2:1;
      for(i=0;i<n;i++)x[i]-=minnum; /* translate sample values to be >= 0 */
      smallnum=(numberone<numbertwo)?1:2;
      onesum-=numberone*minnum; /* and modify sums and maximum value */
      twosum-=numbertwo*minnum;
      maxnum-=minnum;

      if(smallgroup != smallnum)
            { /* if necessary, reflect values so smaller group has smaller mean */
            for(i=0;i<n;i++)x[i]=maxnum-x[i];
            onesum=numberone*maxnum-onesum;
            twosum=numbertwo*maxnum-twosum;
            }
      sum=(smallnum == 1)?onesum:twosum;
      m=(numberone<numbertwo)?numberone:numbertwo; /* number in smaller group */

bin=sum;    

      qsort(x,n,sizeof(int),comp);
      for(i=0;i<n;i++)
            x[i]=(x[i]<bin+1)?x[i]:bin+1; /* any sample values greater than
            the test statistic "bin" will cause b[>bin] to be occupied and are set
            to bin+1, the overflow bin NOW rather than later so as not to cause integer
            overflow */
      *binsreq=(bin+2)*m;
b=calloc(*binsreq,sizeof(float));
      if(!b) 
            return(-1.); /* not enough workspace */
/*      for(i=0;i<*binsreq;i++)
            b[i]=0.0;  initialise work array */ 
      for(i=0;i<m;i++)dope[i]=i*(bin+2); /* offsets to mimic 2 dim. array */
      ++b[dope[0]+x[0]]; /* start histogram off */
      /* for 2 dim. array code would be  ++b[0][x[0]];    */

      n0=x[0];
      for(s=1;s<n;s++)
            { 
            kbot=(m+s-n>1)?m+s-n:1;
            n0dash=n0+x[s]; /* translate each level by x[s] */
            for(k=(s>(m-1))?(m-1):s;k>=kbot;k--) /* only build up histograms
                                    as far as mth level, which is the one we want */
                  {
                  if(n0dash>bin)for(ndash=bin-x[s]+1;ndash<=bin+1;ndash++)
                        /* translate k-1 th level onto kth level...these terms go into overflow bin */
                        b[dope[k]+bin+1]+=b[dope[k-1]+ndash]; 
                  /* for 2 dim. array...    b[k][bin+1]+=b[k-1][ndash]; */
                  ntop=(bin>n0dash)?n0dash:bin;
                  for(ndash=x[s];ndash<=ntop;ndash++){
                        /* these terms just translate along, but may go into overflow bin */
                        b[dope[k]+ndash]+=b[dope[k-1]+ndash-x[s]];
                        /* for 2 dim. array...  b[k][ndash]+=b[k-1][ndash-x[s]]; */
                        }
                  }

            ++b[dope[0]+x[s]]; /* add new term to lowest level */
            /* for 2 dim. array...  ++b[0][x[s]]; */
            n0=n0dash;
            }

      tail=0.; /* find tail probability */
      for(i=0;i<=bin;i++)
            /* for 2 dim. array...   tail+=b[m-1][i]; */
            tail+=b[dope[m-1]+i];

      return(tail/(tail+b[dope[m-1]+bin+1]));
      }

/* -----------------------------------------------------------------------*/
/* psydata */

-1 -2 1 -7 -5 0 3 0 -6 -1 -5 -1 -1 5 -1 -8 0 5 -2 -2 0 -4 -4 2
-1 -1 -1 4 0 1 -2 -2 -2 1 0 -4 -5 0 -5 -1 5 2 0

/* -----------------------------------------------------------------------*/
/* agedata */

1.5  2
2.5 2
1 2
6 2
7 2
2 2
9 2
8 1
5 1
6 1
5 1
6 1
7 1
4 1
9 2
6 2
9 2
7 2
4 2
4 2
5 1
5 1
4 1
7 1
9 1
6 1
6 1
6 1
5 2
3 2
1 2
3 2
5 1
4.5 1
6 1
4 1
6 1
4.5 1
5 1
9 1
5 2
3 2
3.5 2
5 2
7 2
3 2
3.5 2
7 2
3 2

/* -----------------------------------------------------------------------*/
/* iagedat */

3  2
5 2
2 2
12 2
14 2
4 2
18 2
16 1
10 1
12 1
10 1
12 1
14 1
8 1
18 2
12 2
18 2
14 2
8 2
8 2
10 1
10 1
8 1
14 1
18 1
12 1
12 1
12 1
10 2
6 2
2 2
6 2
10 1
9 1
12 1
8 1
12 1
9 1
10 1
18 1
10 2
6 2
7 2
10 2
14 2
6 2
7 2
14 2
6 2
