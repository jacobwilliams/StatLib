Unofficial mirror of http://lib.stat.cmu.edu/apstat/

## StatLib---Applied Statistics algorithms

The Royal Statistical Society has been publishing algorithms in its journal
Applied Statistics since 1968.  As of 1989, there are over 250 of them.
Most are in Fortran, though a few were in Algol, and some recent ones have
been in Pascal.   The book - Applied Statistics Algorithms by Griffiths, P.
and Hill, I.D., Ellis Horwood: Chichester (1985) contains translations of
several algorithms from Algol to Fortran.   A few of the other algorithms
have been supplied in Fortran translations, though a few are in Algol or
Pascal.   The index which follows indicates which algorithms are not in
Fortran.

The full source code is published in the journal.   Those available here
have been transcribed manually, mainly within CSIRO Division of Mathematics &
Statistics, or have been supplied directly by the authors or by the RSS
Algorithms Editor.   In some cases, later corrections or improvements
published in Applied Statistics, have been incorporated.

It is the policy of the editors of the algorithms section of Applied Statistics
that algorithms do not use double precision.   Many of the algorithms here have
been converted to double precision, though users should be careful to check
what precision is used.   Many of the algorithms require the use of other
algorithms, particularly functions for the gamma function and the normal
distribution function.   Where such functions or subroutines are required,
an appropriate comment has been added to the algorithm.

In a few cases, alternative algorithms from other sources have been added.
For instance, three algorithms are included in the file for as66
for calculating the area under the normal curve, and an alternative random
number generator is provided in the file for as183.   The Applied Statistics
algorithm for Nelder-Mead simplex minimization (AS 47) does not include the
fitting of a quadratic surface; the CSIRO/ Rothamsted implementation which
does this is included with AS 47.

Users must consult the original journal articles for details of the calling
arguments; they are not included with these algorithms.

*** Warning.   In some cases, there are different arguments (usually more)
in the Griffiths & Hill versions of these algorithms.   Users should check
this.

It has been assumed that the user will be using a Fortran-77 compiler, so
functions such as alog and amin1 have been converted to their generic forms
(log and min).   This simplifies conversion of code between single and
double precision.   Also all constants in the code, such as 1.0, 0.d0, etc.,
have been replaced with one, zero, etc., and these are defined in either data
or parameter statements.   Many of the algorithms have been entered in lower
case, which is not acceptable in Fortran, though most compilers accept
it.

Some of the algorithms need machine-dependent constants.   The user should
check this.   In most cases, these have been set for compilers which allow
a range of floating-point numbers from about 10**(-37) to 10**(+37), though
most modern compilers allow a much wider range in double precision.

No guarantee is given that the algorithms have been entered correctly, or
that they perform as claimed in the journal.

To obtain an algorithm, send an E-mail request of the form:
<listing>
  send index from apstat
  send 207 from apstat
</listing>
to statlib@lib.stat.cmu.edu

The Royal Statistical Society holds the copyright to these routines, but
has given its permission for their distribution provided that no fee is
charged.

The full collection of Applied Statistics algorithms is very large.
Please only request those algorithms which you need.  Requesting large
numbers of algorithms places a great strain on the StatLib system and
the underlying mail networks.

<HR>
<BLINK>
As of the end of 1997 the Applied Statistics journal does NOT accept algorithms
</BLINK>

### Listing of available Applied Statistics algorithms

No.  | Brief description (volume number/year in brackets)
---  | ---
[3](src/3.f) |        Student's t-distribution. (17/1968)  See also AS 27.
[5](src/5.f) |        The non-central t-distribution. (17/1968)  See also AS 243.
[6](src/6.f) |        Cholesky decomposition of a symmetric +ve definite matrix. (17/1968)
[7](src/7.f) |        Inversion of a symmetric matrix stored in triangular form. (17/1968)
[13](src/13.f) |       Minimum spanning tree. (18/1969)  See also AS 40.
[14](src/14.f) |       Printing the minimum spanning tree. (18/1969)
[15](src/15.f) |       Single-linkage cluster analysis. (18/1969)
[22](src/22.f) |       Calculate treatment effects in a complete factorial experiment for any numbers of levels of factors using the extended Yate's method. (19/1970)
[27](src/27.f) |       Upper tail area under Student's t-distribution. (19/1970) See also AS 3.
[30](src/30.f) |       Half-normal plotting. (19/1970)
[32](src/32.f) |       The incomplete gamma integral. (19/1970)   See also AS 239.
[34](src/34.f) |       Update inverse of symmetric banded matrix. (19/1970)
[38](src/38.f) |       Calculate R-squared for all possible regression subsets using the Gauss-Jordan method. (20/1971)  See also AS 268.
[40](src/40.f) |       Update a minimum spanning tree. (20/1971)
[41](src/41.f) |       Updates corrected sums of squares and products matrices. (20/1971) See also AS 240.
[45](src/45.f) |       Histogram plotting. (20/1971)
[46](src/46.f) |       Gram-Schmidt orthogonalization. (20/1971)
[47](src/47.f) |       Nelder & Mead simplex method of unconstrained minimization without requiring derivatives.  Does not include the quadratic surface fitting part of the algorithm.  A CSIRO/Rothamsted version of the algorithm, which does include fitting a quadratic surface, is also included. (20/1971) (Updated, 20/Dec/93)
[51](src/51.f) |       Log-linear fit for contingency tables. (21/1972)  See also AS 160.
[52](src/52.f) |       Calculation of sums of powers (up to 4) of deviations from the mean. (21/1972)
[53](src/53.f) |       Wishart variate generator. (21/1972)
[57](src/57.f) |       Printing multi-dimensional tables. (22/1973)
[58](src/58.f) |       Allocates observations to clusters to minimize within-cluster sum of squares. (22/1973)  See also AS 136.
[60](src/60.f) |       Eigenvalues/vectors of a symmetric matrix. (22/1973)
[62](src/62.f) |       Distribution of the Mann-Whitney U statistic. (22/1973)
[63](src/63.f) |       Incomplete beta function. (22/1973)  See also TOMS algorithm 708. TOMS algorithms are available from netlib.
[64](src/64.f) |       Inverse of the incomplete beta function ratio. (22/1973) The file here is actually the Griffiths & Hill version of AS 109.
[65](src/65.f) |       Expands structure formula to a list of binary integers.   This is actually remark R82 which replaces the original AS65. (39/1990)
[66](src/66.f) |       The normal distribution function.   Two other algorithms (not from Applied Statistics) have also been included. (22/1973)
[75](src/75.f) |       Algorithms for least-squares calculation using square-root free planar rotations (Morven Gentleman's package). (23/1974)
[76](src/76.f) |       An integral useful in calculating noncentral t and bivariate normal probabilities. (23/1974)
[77](src/77.f) |       Calculate exact null distribution of the largest root of a beta matrix. (23/1974)
[78](src/78.f) |       The mediancentre (i.e. the median in a multi-dimensional space). (23/1974)   See also AS 143.
[83](src/83.f) |       Complex discrete fast Fourier transform. (24/1975)
[84](src/84.f) |       Measures of multivariate skewness and kurtosis. (24/1975)
[88](src/88.f) |       Generate all nCr combinations by simulating nested Fortran DO-loops. (Jane Gentleman's routines). (24/1975)  See also AS 172.
[89](src/89.f) |       Tail probabilities for Spearman's rho. (24/1975)
[91](src/91.f) |       Percentage points of the chi-squared distribution. (24/1975)
[93](src/93.f) |       Calculates frequency distribution for the Ansari-Bradley test statistic. (25/1976)  A routine has been added to return the distribution function.
[95](src/95.f) |       Maximum likelihood estimation of scale and location parameters from grouped data.   User's distribution function. (25/1976)
[96](src/96.f) |       Finding `nice' scales for graphs. (25/1976)
[97](src/97.f) |       Real discrete fast Fourier transform.  Series length must be a power of 2.  (25/1976)  See also AS 117 and AS 176.
[99](src/99.f) |       Fitting Johnson curves by moments. (25/1976)
[100](src/100.f) |      Normal-Johnson and Johnson-Normal transformations. (25/1976)
[103](src/103.f) |      Psi or digamma function. (25/1976)
[107](src/107.f) |      Calculate operating characteristics and average sampling number for a general class of sequential sampling plans. (26/1977)
[108](src/108.f) |      Multiple linear regression minimizing the sum of absolute errors. (26/1977)  See also AS 238 (in Pascal).
[109](src/109.f) |      Inverse of the incomplete beta function. (26/1977)
[110](src/110.f) |      LP-Norm fit of straight line by extension of Schlossmacher. Updated 14.06.2013 (26/1977)
[111](src/111.f) |      Percentage points of the normal distribution. (26/1977)  See also AS 241.
[114](src/114.f) |      Compute the numerator of certain ordinal measures of association (Kendall's tau, Somer's d, Goodman and Kruskal's gamma) when the data are ordered categories. (26/1977)
[116](src/116.f) |      Calculate the tetrachoric correlation and its standard errors. (26/1977)
[117](src/117.f) |      Fast Fourier transform for series of any length using the CHIRP algorithm. (26/1977)
[121](src/121.f) |      Trigamma function. (27/1978)
[123](src/123.f) |      Distribution function of mixtures of beta distributions. (27/1978)
[125](src/125.f) |      Maximum likelihood estimation for censored exponential survival data with covariates. (27/1978)
[126](src/126.f) |      Distribution function of the range for the normal distribution. (27/1978)
[127](src/127.f) |      Generation of random orthogonal matrices. (27/1978)
[128](src/128.f) |      Computes approximate covariance matrix for normal order statistics. (27/1978)
[132](src/132.f) |      Simple regression minimizing the sum of absolute deviations. (27/1978)
[133](src/133.f) |      Finding the global maximum or minimum of a function of 1 variable. (27/1978)
[134](src/134.f) |      Generate random beta variates for alpha &lt 1 and beta &gt 1. (28/1979)
[135](src/135.f) |      Min-Max (L-infinity) estimates for linear multiple regression. (28/1979)
[136](src/136.f) |      A K-means clustering algorithm. (28/1979)
[138](src/138.f) |      Maximum likelihood estimates of the mean and standard deviation of the normal distribution with censored or confined observations. (28/1979)
[139](src/139.f) |      Maximum likelihood estimation in a linear model from confined and censored normal data. (28/1979)
[140](src/140.f) |      Clustering the nodes of a directed graph. (28/1979)
[141](src/141.f) |      Inversion of a symmetric matrix ignoring a specified row/column. (28/1979)
[142](src/142.f) |      Exact tests of significance in binary regression. (28/1979)
[143](src/143.f) |      Calculates the median centre. (28/1979)
[145](src/145.f) |      Exact distribution of the largest multinomial frequency. (28/1979)
[147](src/147.f) |      Incomplete gamma function. (29/1980)  See also AS 239.
[148](src/148.f) |      Removal of bias in the jackknife procedure. (This is actually ASR 62) (29/1980)
[149](src/149.f) |      Amalgamation of means in the case of simple ordering ('Up-and-Down Blocks' algorithm for isotonic regression). (29/1980)
[150](src/150.f) |      Computes estimate of spectrum of a point process using a centered moving average of the periodogram of the counting process. (29/1980)
[151](src/151.f) |      Smoothed spectral estimate for bivariate point processes. (29/1980)
[152](src/152.f) |      Cumulative hypergeometric probabilities.   (This is actually AS R77) (29/1980, revised in 38/1989)
[153](src/153.f) |      Distribution of weighted sum of squares of normal variables. (29/1980) Pan's procedure for the tail probabilities of the Durbin-Watson statistic.
[154](src/154.f) |      Exact maximum likelihood estimation of autoregressive-moving average models by Kalman filtering. (29/1980)  See also AS 182.
[155](src/155.f) |      Distribution function of a linear combination of non-central chi- squared random variables. (29/1980)  See also AS 204.   This is a Fortran translation supplied by the author.
[157](src/157.f) |      The runs-up and runs-down tests. (30/1981)
[158](src/158.f) |      Calculation of probabilities for inferences under simple order restrictions. (30/1981)  See also AS 198.
[159](src/159.f) |      Generate random 2-way table with given marginal totals. (30/1981)
[160](src/160.f) |      Partial and marginal association in multi-dimensional contingency tables. (30/1981)
[161](src/161.f) |      Critical regions of an unconditional non-randomized test of homogeneity in 2 x 2 contingency tables. (30/1981)
[162](src/162.f) |      Multivariate Conditional Logistic Analysis of  Stratum-matched Case-control Studies. (30/1981)   Includes a Fortran version of CACM algorithm 382 for generating all combinations of M out of N items.
[163](src/163.f) |      A Givens Algorithm for Moving from one Linear Model to  another without Going back to the Data. (30/1981)
[164](src/164.f) |      Least squares subject to linear constraints. (30/1981)
[165](src/165.f) |      Discriminant analysis of categorical data. (30/1981)
[166](src/166.f) |      Calculates the entanglement matrix for a specified design. (30/1981)
[167](src/167.f) |      Calculates efficiencies of estimation and their generalized inverse. (30/1981)
[168](src/168.f) |      Calculates 'neat' values for plotting scales. (30/1981)
[169](src/169.f) |      Produces scatter plots. (30/1981)
[170](src/170.f) |      Computation of probability and non-centrality parameter of a non- central chi-square distribution. (30/1981)
[171](src/171.f) |      Fisher's exact variance test for the Poisson distribution. (31/1982)
[172](src/172.f) |      Generates indices for simulated nested DO-loops.   Actually converts (either way) between a single index and a vector of subscripts. (31/1982)
[173](src/173.f) |      Generates design matrix for balanced factorial experiments. (31/1982)
[174](src/174.f) |      Multivariate rank sum test and median test. (31/1982)
[175](src/175.f) |      Cramer-Wold factorization of self-reciprocal polynomials as the product of two polynomials. (31/1982)
[176](src/176.f) |      Kernel density estimation using the fast Fourier transform. (31/1982) Also contains an alternative set of routines for density estimation.
[177](src/177.f) |      Expected values of normal order statistics. (31/1982)
[178](src/178.f) |      Gauss-Jordan sweep operator with multi-collinearity detection. (31/1982)
[179](src/179.f) |      Enumeration of all permutations of multi-sets with fixed repetition numbers. (31/1982)
[180](src/180.f) |      Linear rank estimate of the standard deviation after symmetric trimming. (31/1982)
181     | Withdrawn.  See R94 below
[182](src/182.f) |      Finite-sample prediction from ARIMA processes. (31/1982)  Uses AS 154.
[183](src/183.f) |      The Wichmann & Hill random number generator.   An alternative is also provided. (31/1982)
[184](src/184.f) |      Non-central studentized maximum and related multiple-t probabilities. (31/1982)
[185](src/185.f) |      Backward elimination procedure to find best-fitting log-linear models for contingency tables. (31/1982)  Uses AS 51.
[186](src/186.f) |      Discrete Fast Fourier Transform with data permutation.   Series length must be a power of 2. (31/1982)
[187](src/187.f) |      Derivatives of the incomplete gamma integral. (31/1982)
[188](src/188.f) |      Estimation of the order of dependence in sequences. (32/1983)
[189](src/189.f) |      Maximum likelihood estimation for the beta binomial distribution. (32/1983)
[190](src/190.f) |      Distribution function & its inverse, for the studentized range. (32/1983)
[191](src/191.f) |      Approximate likelihood calculation for ARMA and seasonal ARMA models. Includes a routine for the Banachiewicz (or modified Cholesky) factorization A = LDL'. (32/1983)
[192](src/192.f) |      Calculate approximate percentage points using Pearson curves and the first three or four moments. (32/1983)
[193](src/193.f) |      The Knuth spectral test for congruential random number generators. (32/1983)
[194](src/194.f) |      Test of goodness of fit of ARMA models. (32/1983)
[195](src/195.f) |      Multivariate normal probabilities for a rectangular region in N- dimensional space.  A version which calls IMSL routines is also included. (33/1984)  See AS 251 for a special case.
[196](src/196.f) |      Conditional multivariate logistic analysis of stratified case-control studies. (33/1984). Edited for errors by Bill Bardsley (bill.bardsley@man.ac.uk) on 18/08/2002.
[197](src/197.f) |      Likelihood function for an ARMA process. (33/1984)
[198](src/198.f) |      Calculation of level probabilities for order-restricted inference. (33/1984)
[199](src/199.f) |      Branch and bound algorithm to find the subset which maximizes a quadratic form. (33/1984)
[200](src/200.f) |      Approximate the sum of squares of normal score. (33/1984)
[201](src/201.f) |      Combine predictions about a statistic based on the orderings of a set of means with an F-test of differences between the means. (33/1984)
[202](src/202.f) |      Data-based nonparametric hazard estimation. (33/1984)
[203](src/203.f) |      Maximum likelihood estimation of mixtures of distributions (normal, exponential, Poisson and binomial). (33/1984)   See also AS221. This is a translation from Algol into Fortran.
[204](src/204.pas) |      Distribution of a sum of non-central chi-squared variables. (33/1984)   Translation from Algol into Pascal.
[205](src/205.f) |      Enumeration of all R x C contingency tables with given row and column totals, and calculation of hypergeometric probability for each table. (33/1984)
[206](src/206.f) |      Isotonic regression in two independent variables. (33/1984)
[207](src/207.f) |      Fit a generalized log-linear model. (33/1984)
[208](src/208.f) |      Fit multivariate logistic model by the method of moments. (34/1985)
[209](src/209.f) |      The distribution function of skewness and kurtosis. (34/1985)
[210](src/210.f) |      Fit 5-parameter Johnson SB curves by moments. (34/1985)
[211](src/211.f) |      The F-G diagonalization algorithm. (34/1985)  Modifications in ASR 71 & ASR 74 have been added to the file.
[212](src/212.f) |      Fitting the exponential curve by least squares. (34/1985)
[213](src/213.f) |      Generation of correlation matrices with specified eigenvalues. (34/1985)
[214](src/214.f) |      Calculation of Monte Carlo confidence intervals. (34/1985)
[215](src/215.f) |      Maximum likelihood estimation of the parameters of the generalized extreme-value distribution.  (34/1985). Updated on [15/Nov/99]
[216](src/216.f) |      Maximum likelihood fitting when model contains a linear part and auxiliary parameters. (34/1985)
[217](src/217.f) |      Computation of the Dip statistic to test for unimodality. (34/1985). Bug fixed by Ferenc Mechler (fmechler@med.cornell.edu). See [README file](src/217-readme.doc) for details. File changed to incorporate [this patch](src/217.patch) provided by Yong Lu [5/Sep/02] [4/Aug/05]
[218](src/218.f) |      Elements of the Fisher information matrix for the smallest extreme- value distribution, also allows for censored data. (35/1986)
[219](src/219.f) |      Height balanced trees. (35/1986)
[220](src/220.f) |      Operating characteristics of James-Stein and Efron-Morris estimation. (35/1986)
[221](src/221.f) |      Maximum likelihood estimation of a mixing distribution. (35/1986)
[222](src/222.f) |      Resistant smoothing using the fast Fourier transform. (36/1987)
[223](src/223.f) |      Optimum ridge parameter selection.  (36/1987)
[224](src/224.f) |      Combining component designs to form a design with several orthogonal blocking factors. (36/1987)
[225](src/225.f) |      Minimizing linear inequality-constrained Mahalanobis distances. (36/1987)
[226](src/226.f) |      Cumulative probabilities for the non-central beta distribution. (36/1987)
[227](src/227.f) |      Generate all possible N-bit binary codes. (36/1987)
[228](src/228.f) |      Finding I-projections (maximizing cross-entropy) subject to a finite set of linear inequality constraints. (36/1987)
[229](src/229.f) |      Quantile regression. (36/1987)
[230](src/230.f) |      Distribution of customers in M/G/m queues using an approximation due to Hokstad. (36/1987)
[231](src/231.pas) |      Distribution of a noncentral chi-squared variable with non-negative degrees of freedom. (36/1987) In Pascal.
[232](src/232.f) |      Computation of population and sample correlation and partial correlation matrices in MARMA(p,q) time series. (37/1988)
[233](src/233.f) |      Branch and bound algorithm for feature subset selection.  (37/1988)
[234](src/234.f) |      Approximating the percentage points of simple linear rank statistics with Cornish-Fisher expansions. (37/1988)
[235](src/235.f) |      Tally frequencies of distinct real values. (37/1988)
[236](src/236.pas) |      Recursive enumeration of RxC tables for exact likelihood evaluation. (37/1988) In Pascal.
[237](src/237.f) |      The corner method for identifying autoregressive moving average models. (37/1988)
[238](src/238.pas) |      Recursive procedure for the L1-norm fitting of a straight line. (37/1988) In Pascal.
[239](src/239.f) |      Incomplete gamma function. (37/1988)
[240](src/240.f) |      Updating the inverse of corrected sums of squares and products. (37/1988)
[241](src/241.f) |      Inverse normal - more accurate than AS 111. (37/1988)
[242](src/242.f) |      The exact likelihood of a vector autoregressive moving average model. (38/1989)
[243](src/243.f) |      Cumulative distribution function of the non-central t-distribution. (38/1989)
[244](src/244.pas) |      Decomposability and collapsibility for log-linear models. (38/1989) In Pascal.
[245](src/245.f) |      Log of the gamma function. (38/1989)   The file includes a slower but more accurate algorithm which is not an AS algorithm.
[246](src/246.f) |      Repeated measurements analysis of variance with unknown autoregressive parameter. (38/1989)
[247](src/247.pas) |      Updating the sufficient configurations for fitting ZPA (zero partial association) models to multi-dimensional contingency tables. (38/1989) In Pascal.
[248](src/248.f) |      Empirical distribution function goodness-of-fit tests for the uniform, normal and exponential distributions. (38/1989)
[249](src/249.f) |      Mean and covariance estimation for the truncated multivariate normal distribution. (38/1989)
[250](src/250.f) |      Test the equality of dispersion matrices using methods due to Puri & Sen (distribution-free) and Anderson (normal distribution). (38/1989)
[251](src/251.f) |      Multivariate normal probability integrals with product correlation structure. (38/1989)
[252](src/252.f) |      Generating classes for log-linear models. (39/1990)
[253](src/253.f) |      Maximum likelihood estimation of the RC(M) association model. (39/1990)
[254](src/254.f) |      Fit mixture of normal distributions to grouped & truncated data. (39/1990)
[255](src/255) |      Fitting two-way tables by means for rows, columns and cross-term. (39/1990) In Algol.
[256](src/256.pas) |      Distribution of a quadratic form in normal variables. (39/1990) In Pascal.
[257](src/257.f) |      Isotonic regression for umbrella orderings. (39/1990)
[258](src/258.f) |      Average run lengths for cumulative sum schemes. (39/1990)
[259](src/259.f) |      Extending confidence intervals by the Robbins-Monro search process. (39/1990)
[260](src/260.f) |      Distribution function of the square (R^2) of the sample multiple- correlation coefficient. (40/1991)
[261](src/261.f) |      Quantiles of the distribution of square (R^2) of the sample multiple- correlation coefficient. (40/1991)
[262](src/262.f) |      The Wei-Lachin two-sample test, the generalized Wilcoxon test and the log-rank test for incomplete multivariate observations with censoring. (40/1991)
[263](src/263.f) |      Construction of irredundant test sets. (40/1991)
[264](src/264.f) |      Printing of bit patterns. (40/1991)
[265](src/265.f) |      Waiting time distribution for the G/G/1 queue using the Fast Fourier transform. (40/1991)
[266](src/266.f) |      Maximum likelihood estimation of the parameters of the Dirichlet distribution. (40/1991)
[267](src/267.f) |      Calculate lower bound for the probability of correct selection of a subset containing the best populations from a set of populations. (40/1991)
[268](src/268.f) |      Calculates statistics for all possible subset regressions using an orthogonal decomposition. (40/1991)
[269](src/269.f) |      Calculates the Cornish-Fisher adjustment to distribution functions (from the normal distribution function) using higher cumulants. (41/1992)
[270](src/270.f) |      Maximum likelihood fitting of a 'key' density function, e.g. normal distribution, multiplied by a series in Hermite polynomials. (41/1992)
[271](src/271.f) |      Optimal (smallest mis-classification probability) joint classification procedure. (41/1992)  See also AS 276.
[272](src/272.f) |      Produce character Box-plots using supplied quartiles. (41/1992)
[273](src/273.f) |      Compare the least-squares fit of two subsets of regression variables using Spj0tvoll's method. (41/1992). Updated [04/01/03].
[274](src/274.f) |      Least squares algorithms to supplement those of Gentleman in AS 75. (41/1992)
[274-90](src/274.f90) |      Algorithm 274, translated into Fortran 90.
[275](src/275.f) |      Non-central chi-squared distribution function, number of degrees of freedom must be positive but not necessarily integer. (41/1992)
[276](src/276.f) |      Optimal combinatoric classification for two normal classes. (41/1992) See also AS 271.
[277](src/277.f) |      The Oja bivariate median. (41/1992)
[278](src/278.f) |      The distribution of quadratic forms of multivariate generalized Student variables. (41/1992)
[279](src/279.f) |      P-values for the generalized/alternative Durbin-Watson statistic with arbitrary lag using a Cholesky transformation method. (42/1993)
[280](src/280.f) |      The power function for Fisher's exact test. (42/1993)
[281](src/281.f) |      Scaling and rounding regression coefficients to force them to be integers. (42/1993)
[282](src/282.f) |      High breakdown regression and multivariate estimation. (42/1993)
[283](src/283.c) |      Rapid computation of the permutation paired and grouped t-tests. (42/1993)  In ANSI Standard C.
[284](src/284.f) |      Null distribution of a statistic for testing sphericity and additivity: a Jacobi polynomial expansion. (42/1993)
[285](src/285.f) |      Calculate multivariate normal probabilities using Monte Carlo methods for a general region defined by a user-supplied function. (42/1993)
[286](src/286.f) |      Estimation of parameters when there are errors in both the predictor variables and the dependent variable using least squares with a general model (includes implicit non-linear regression). (42/1993)
[287](src/287.f) |      Adaptive rejection sampling (to generate random variables) from log-concave density functions. (42/1993)
[288](src/288.f) |      P-value calculation for the generalized two-sample Smirnov tests. The tests are conditional on ties in the pooled sample. (43/1994)
[289](src/289.f) |      Convolves hypergeometric distributions generated by  several 2x2 tables (43/1994)
[290](src/290.f) |      Generate a grid of variance ratios for contour plotting of confidence regions for a pair of parameters in non-linear regression. (43/1994)
[291](src/291.f) |      Calculates overlap capability of subsequence SSEQ of length L. (43/1994)
[292](src/292.f) |      Computes Fisher information matrix elements for time (type I) or failure (type II) censored units from the smallest extreme value (sev), largest extreme value (lev), normal or logistic distribution.   (43/1994)
[293](src/293.f) |      Convolves conditional distributions generated by several 2xK tables (43/1994)
[294](src/294.pas) |      Pascal program for ....?
[295](src/295.f) |      Heuristic algorithm to pick N rows of X out of NCAND to maximize the determinant of X'X, using the Fedorov exchange algorithm. (43/1994)
[R94](src/R94.f) |      calculates Shapiro-Wilk normality test and P-value for sample sizes 3 &lt;= n &lt;= 5000 .  Handles censored or uncensored data. Corrects AS 181, which was found to be inaccurate for n > 50. (Uses function POLY in <a href="181">AS 181</a>) by Patrick Royston (44/1995)
[R96](src/R96.f) |     Interval for trapezoidal rule integration to limit error using moment generating function bound for weighted sum of chi-squared(1) random variables.
[297](src/297.f) |     Nonparametric regression using ultraspherical polynomials.
[298](src/298.f) |     Hybrid minimization routine using simulated annealing and any local minimizer supplied by the user.
[300](src/300.f) |     Efficient algorithm for determining a simulated percentile point with a fixed degree of accuracy.
[301](src/301.f) |     Computes the logarithms of F1, P(H) and P(H').
[302](src/302.f) |     Subroutine to compute a multiple isotonic regression for umbrella ordering with known peak.
[303](src/303.f) |     Generation of ordered multinomial frequencies. Given integers N (N >= 0) and K (K >= 1), generates all possible integer arrays NI of size K, whose elements are non-negative, non-decreasing and sum to N.
[304](src/304.f) |     Fisher's non-parametric randomization test for two small independent random samples.
[305](src/305.f) |     Compute the ARL of a CUSUM mean chart with linear drift.
[306](src/306.f) |     Calculation of a BIB product operation on any two matrices.
[307](src/307.f) |     Calculation of the simplicial depth and the halfspace depth
[310](src/310.f) |     Computes the cumulative distribution function of a non-central beta random variable.
[314](src/314.f) |     Inverts matrix with contents subject to modulo arithmetic.
[317](src/317.pas) |   Maximum Likelihood Estimation and Goodness-of-Fit Tests for Mixtures of Distributions. Submitted by Alan Miller (amiller@bigpond.net.au) [31/Oct/02]
[319](src/319.f)   |     A program to implement a quasi-newton method. Using new algorithm Varmet   July 1994
[319](src/319.f90) |     Fortran 90 version of 319.

<HR>
N.B. This library was maintained by Alan Miller, formerly of  CSIRO Division of Mathematics and Statistics, Clayton, Victoria 3169.
<br>

Fortran 90 translations of some of the above algorithms can be found
at https://jblevins.org/mirror/amiller/#apstat.

Recent algorithms are provided by "T.R.Hopkins" &lt;T.R.Hopkins@ukc.ac.uk&gt;.
The StatLib collection of Applied Statistics Algorithms is maintained
by Pantelis Vlachos, vlachos@stat.cmu.edu.

### Credit where credit is due
If you use an algorithm, dataset, or other information from StatLib,
please acknowledge both StatLib and the original contributor of the
material.

### Notes

 * Original site: Last modified: Thu Aug  4 09:57:59 EDT 2005 by <A HREF="http://lib.stat.cmu.edu/master/vlachos.html">Pantelis Vlachos</A>

### See also

* https://jblevins.org/mirror/amiller/#apstat 
* https://wp.csiro.au/alanmiller/apstat.html 