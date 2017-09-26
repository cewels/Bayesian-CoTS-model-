#include "Misc.h"

using namespace std;

/*double gamma_sample_v1(double a, double b, double c)
//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_SAMPLE samples the Gamma PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2004
//
//  Author:
//
//    Original FORTRAN77 version by Ahrens and U Dieter.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Joachim Ahrens, Ulrich Dieter,
//    Generating Gamma variates by a modified rejection technique,
//    Communications of the ACM,
//    Volume 25, Number 1, January 1982, pages 47-54.
//
//    Joachim Ahrens, Ulrich Dieter,
//    Computer methods for sampling from Gamma, Beta, Poisson and
//    binomial distributions,
//    Computing,
//    Volume 12, Number 3, September 1974, pages 223-246.
//
//    Joachim Ahrens, Klaus-Dieter Kohrt, Ulrich Dieter,<br>
//    Algorithm 599:
//    Sampling from Gamma and Poisson Distributions,<br>
//    ACM Transactions on Mathematical Software,<br>
//    Volume 9, Number 2, June 1983, pages 255-257.
//
//  Parameters:
//
//    Input, double A, B, C, the parameters of the PDF.
//    0.0 < B,
//    0.0 < C.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double GAMMA_SAMPLE, a sample of the PDF.
//
{
double a1 = 0.3333333;
double a2 = -0.2500030;
double a3 = 0.2000062;
double a4 = -0.1662921;
double a5 = 0.1423657;
double a6 = -0.1367177;
double a7 = 0.1233795;
double bcoef;
double co;
double d;
double e;
double e1 = 1.0;
double e2 = 0.4999897;
double e3 = 0.1668290;
double e4 = 0.0407753;
double e5 = 0.0102930;
double euler = 2.71828182845904;
double p;
double q;
double q0;
double q1 = 0.04166669;
double q2 = 0.02083148;
double q3 = 0.00801191;
double q4 = 0.00144121;
double q5 = -0.00007388;
double q6 = 0.00024511;
double q7 = 0.00024240;
double r;
double s;
double si;
double s2;
double t;
double u;
double v;
double w;
double x;

//
//  Allow C = 0.
//
if (c == 0.0)
{
x = a;
return x;
}
//
//  C < 1.
//
if (c < 1.0)
{
for (;;)
{
u = r8_uniform_01_sample();
t = 1.0 + c / euler;
p = u * t;

s = r8_exponential_01_sample();

if (p < 1.0)
{
x = exp(log(p) / c);
if (x <= s)
{
break;
}
}
else
{
x = -log((t - p) / c);
if ((1.0 - c) * log(x) <= s)
{
break;
}
}
}

x = a + b * x;
return x;
}
//
//  1 <= C.
//
else
{
s2 = c - 0.5;
s = sqrt(c - 0.5);
d = sqrt(32.0) - 12.0 * sqrt(c - 0.5);

t = r8_normal_01_sample();
x = (sqrt(c - 0.5) + 0.5 * t)*(sqrt(c - 0.5) + 0.5 * t);

if (0.0 <= t)
{
x = a + b * x;
return x;
}

u = r8_uniform_01_sample();

if (d * u <= t * t * t)
{
x = a + b * x;
return x;
}

r = 1.0 / c;

q0 = ((((((
q7   * r
+ q6) * r
+ q5) * r
+ q4) * r
+ q3) * r
+ q2) * r
+ q1) * r;

if (c <= 3.686)
{
bcoef = 0.463 + s - 0.178 * s2;
si = 1.235;
co = 0.195 / s - 0.079 + 0.016 * s;
}
else if (c <= 13.022)
{
bcoef = 1.654 + 0.0076 * s2;
si = 1.68 / s + 0.275;
co = 0.062 / s + 0.024;
}
else
{
bcoef = 1.77;
si = 0.75;
co = 0.1515 / s;
}

if (0.0 < sqrt(c - 0.5) + 0.5 * t)
{
v = 0.5 * t / s;

if (0.25 < r8_abs_v1(v))
{
q = q0 - s * t + 0.25 * t * t + 2.0 * s2 * log(1.0 + v);
}
else
{
q = q0 + 0.5 * t * t * ((((((
a7   * v
+ a6) * v
+ a5) * v
+ a4) * v
+ a3) * v
+ a2) * v
+ a1) * v;
}

if (log(1.0 - u) <= q)
{
x = a + b * x;
return x;
}
}

for (;;)
{
e = r8_exponential_01_sample();

u = r8_uniform_01_sample();

u = 2.0 * u - 1.0;
t = bcoef + r8_abs_v1(si * e) * r8_sign_v1(u);

if (-0.7187449 <= t)
{
v = 0.5 * t / s;

if (0.25 < r8_abs_v1(v))
{
q = q0 - s * t + 0.25 * t * t + 2.0 * s2 * log(1.0 + v);
}
else
{
q = q0 + 0.5 * t * t * ((((((
a7   * v
+ a6) * v
+ a5) * v
+ a4) * v
+ a3) * v
+ a2) * v
+ a1) * v;
}

if (0.0 < q)
{
if (0.5 < q)
{
w = exp(q) - 1.0;
}
else
{
w = ((((
e5   * q
+ e4) * q
+ e3) * q
+ e2) * q
+ e1) * q;
}

if (co * r8_abs_v1(u) <= w * exp(e - 0.5 * t * t))
{
x = a + b *(s + 0.5 * t)*(s + 0.5 * t);
return x;
}
}
}
}
}
}
void r8vec_unit_sum_v1(int n, double a[])

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIT_SUM normalizes an R8VEC to have unit sum.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input/output, double A[N], the vector to be normalized.
//    On output, the entries of A should have unit sum.  However, if
//    the input vector has zero sum, the routine halts.
//
{
double a_sum;
int i;

a_sum = 0.0;
for (i = 0; i < n; i++)
{
a_sum = a_sum + a[i];
}

if (a_sum == 0.0)
{
cout << "\n";
cout << "R8VEC_UNIT_SUM - Fatal error!\n";
cout << "  The vector entries sum to 0.\n";
exit(1);
}

for (i = 0; i < n; i++)
{
a[i] = a[i] / a_sum;
}

return;
}
double r8_sign_v1(double x)

//****************************************************************************80
//
//  Purpose:
//
//    R8_SIGN returns the sign of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the number whose sign is desired.
//
//    Output, double R8_SIGN, the sign of X.
//
{
if (x < 0.0)
{
return (-1.0);
}
else
{
return (1.0);
}
}
double r8_abs_v1(double x)

//****************************************************************************80
//
//  Purpose:
//
//    R8_ABS returns the absolute value of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the quantity whose absolute value is desired.
//
//    Output, double R8_ABS, the absolute value of X.
//
{
double value;

if (0.0 <= x)
{
value = x;
}
else
{
value = -x;
}
return value;
}

void dirichlet_sample_v1(int n, double a[], double *(&out))

//****************************************************************************80
//
//  Purpose:
//
//    DIRICHLET_SAMPLE samples the Dirichlet PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Jerry Banks, editor,
//    Handbook of Simulation,
//    Engineering and Management Press Books, 1998, page 169.
//
//  Parameters:
//
//    Input, int N, the number of components.
//
//    Input, double A(N), the probabilities for each component.
//    Each A(I) should be positive.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double DIRICHLET_SAMPLE[N], a sample of the PDF.  The entries
//    of X should sum to 1.
//
{
double a2;
double b2;
double c2;
int i;
double *x = new double[n];

a2 = 0.0;
b2 = 1.0;

for (i = 0; i < n; i++)
{
c2 = a[i];
x[i] = gamma_sample_v1(a2, b2, c2);
}
//
//  Rescale the vector to have unit sum.
//
r8vec_unit_sum_v1(n, x);

for (i = 0; i < n; i++)
{
out[i] = x[i];
}
delete[] x;
return;

}
*/
void dirichlet_sample(int n, double a[], double *(&out), base_generator_type& eng)

//****************************************************************************80
//
//  Purpose:
//
//    DIRICHLET_SAMPLE samples the Dirichlet PDF-based on gamma

//
//  Author:
//
//    me
//
//  Reference:
//
{
	int i;
	double sum;
	Context dist;

	Map<RowVectorXd> v(a, n);
	VectorXd v1(n);

	for (i = 0; i < n; i++)
	{
		v1(i) = dist.gamma_dis(v(i), eng);

	}

	sum = v1.sum();
	for (i = 0; i < n; i++)
	{
		out[i] = v1(i) / sum;
	}
	return;
}
void multi_boost(int n, double p[], int ncat, int *(&ix), base_generator_type& eng)

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_MULTINOMIAL_SAMPLE generates a multinomial random deviate.
//
// modified, using boost
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN77 version by Barry Brown, James Lovato.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Luc Devroye,
//    Non-Uniform Random Variate Generation,
//    Springer, 1986,
//    ISBN: 0387963057,
//    LC: QA274.D48.
//
//  Parameters:
//
//    Input, int N, the number of events, which will be
//    classified into one of the NCAT categories.
//
//    Input, double P[NCAT-1].  P(I) is the probability that an event
//    will be classified into category I.  Thus, each P(I) must be between 
//    0.0 and 1.0.  Only the first NCAT-1 values of P must be defined since 
//    P(NCAT) would be 1.0 minus the sum of the first NCAT-1 P's.
//
//    Input, int NCAT, the number of categories.
//
//    Output, int I4VEC_MULTINOMIAL_SAMPLE[NCAT], a random observation from 
//    the multinomial distribution.  All IX(i) will be nonnegative and their 
//    sum will be N.
//
{
	int i;
	int icat;
	int ntot;
	double prob;
	double ptot;
	Context dist;

	if (n < 0)
	{
		cerr << "\n";
		cerr << "I4VEC_MULTINOMIAL_SAMPLE - Fatal error!\n";
		cerr << "  N < 0\n";
		exit(1);
	}

	if (ncat <= 1)
	{
		cerr << "\n";
		cerr << "I4VEC_MULTINOMIAL_SAMPLE - Fatal error!\n";
		cerr << "  NCAT <= 1\n";
		exit(1);
	}

	for (i = 0; i < ncat - 1; i++)
	{
		if (p[i] < 0.0)
		{
			cerr << "\n";
			cerr << "I4VEC_MULTINOMIAL_SAMPLE - Fatal error!\n";
			cerr << "  Some P(i) < 0.\n";
			exit(1);
		}

		if (1.0 < p[i])
		{
			cerr << "\n";
			cerr << "I4VEC_MULTINOMIAL_SAMPLE - Fatal error!\n";
			cerr << "  Some 1 < P(i).\n";
			exit(1);
		}
	}

	ptot = 0.0;
	for (i = 0; i < ncat - 1; i++)
	{
		ptot = ptot + p[i];
	}

	if (1.00001 < ptot)
	{
		cerr << "\n";
		cerr << "I4VEC_MULTINOMIAL_SAMPLE - Fatal error!\n";
		cerr << "  1.0 < Sum of P().\n";
		exit(1);
	}
	//
	//  Initialize variables.
	//
	ntot = n;
	ptot = 1.0;

	for (i = 0; i < ncat; i++)
	{
		ix[i] = 0;
	}
	//
	//  Generate the observation.
	//


	for (icat = 0; icat < ncat - 1; icat++)
	{
		prob = p[icat] / ptot;
		ix[icat] = dist.binomial_dis(ntot, prob, eng);
		ntot = ntot - ix[icat];
		if (ntot <= 0)
		{
			return;
		}
		ptot = ptot - p[icat];
	}

	ix[ncat - 1] = ntot;

	return;
}

void multi_boost_1(int n, long double p[], int ncat, int *(&ix), base_generator_type& eng)

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_MULTINOMIAL_SAMPLE generates a multinomial random deviate.
//
// modified, using boost
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN77 version by Barry Brown, James Lovato.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Luc Devroye,
//    Non-Uniform Random Variate Generation,
//    Springer, 1986,
//    ISBN: 0387963057,
//    LC: QA274.D48.
//
//  Parameters:
//
//    Input, int N, the number of events, which will be
//    classified into one of the NCAT categories.
//
//    Input, double P[NCAT-1].  P(I) is the probability that an event
//    will be classified into category I.  Thus, each P(I) must be between 
//    0.0 and 1.0.  Only the first NCAT-1 values of P must be defined since 
//    P(NCAT) would be 1.0 minus the sum of the first NCAT-1 P's.
//
//    Input, int NCAT, the number of categories.
//
//    Output, int I4VEC_MULTINOMIAL_SAMPLE[NCAT], a random observation from 
//    the multinomial distribution.  All IX(i) will be nonnegative and their 
//    sum will be N.
//
{
	int i;
	int icat;
	int ntot;
	double prob;
	double ptot;
	Context dist;

	if (n < 0)
	{
		cerr << "\n";
		cerr << "I4VEC_MULTINOMIAL_SAMPLE - Fatal error!\n";
		cerr << "  N < 0\n";
		exit(1);
	}

	if (ncat <= 1)
	{
		cerr << "\n";
		cerr << "I4VEC_MULTINOMIAL_SAMPLE - Fatal error!\n";
		cerr << "  NCAT <= 1\n";
		exit(1);
	}

	for (i = 0; i < ncat - 1; i++)
	{
		if (p[i] < 0.0)
		{
			cerr << "\n";
			cerr << "I4VEC_MULTINOMIAL_SAMPLE - Fatal error!\n";
			cerr << "  Some P(i) < 0.\n";
			exit(1);
		}

		if (1.0 < p[i])
		{
			cerr << "\n";
			cerr << "I4VEC_MULTINOMIAL_SAMPLE - Fatal error!\n";
			cerr << "  Some 1 < P(i).\n";
			exit(1);
		}
	}

	ptot = 0.0;
	for (i = 0; i < ncat - 1; i++)
	{
		ptot = ptot + p[i];
	}

	if (1.00001 < ptot)
	{
		cerr << "\n";
		cerr << "I4VEC_MULTINOMIAL_SAMPLE - Fatal error!\n";
		cerr << "  1.0 < Sum of P().\n";
		exit(1);
	}
	//
	//  Initialize variables.
	//
	ntot = n;
	ptot = 1.0;

	for (i = 0; i < ncat; i++)
	{
		ix[i] = 0;
	}
	//
	//  Generate the observation.
	//


	for (icat = 0; icat < ncat - 1; icat++)
	{
		prob = p[icat] / ptot;
		ix[icat] = dist.binomial_dis(ntot, prob, eng);
		ntot = ntot - ix[icat];
		if (ntot <= 0)
		{
			return;
		}
		ptot = ptot - p[icat];
	}

	ix[ncat - 1] = ntot;

	return;
}

void multi_boost_long(long long n, double p[], int ncat, long long *(&ix), base_generator_type& eng)

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_MULTINOMIAL_SAMPLE generates a multinomial random deviate.
//
// modified, using boost
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN77 version by Barry Brown, James Lovato.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Luc Devroye,
//    Non-Uniform Random Variate Generation,
//    Springer, 1986,
//    ISBN: 0387963057,
//    LC: QA274.D48.
//
//  Parameters:
//
//    Input, int N, the number of events, which will be
//    classified into one of the NCAT categories.
//
//    Input, double P[NCAT-1].  P(I) is the probability that an event
//    will be classified into category I.  Thus, each P(I) must be between 
//    0.0 and 1.0.  Only the first NCAT-1 values of P must be defined since 
//    P(NCAT) would be 1.0 minus the sum of the first NCAT-1 P's.
//
//    Input, int NCAT, the number of categories.
//
//    Output, int I4VEC_MULTINOMIAL_SAMPLE[NCAT], a random observation from 
//    the multinomial distribution.  All IX(i) will be nonnegative and their 
//    sum will be N.
//
{
	int i;
	int icat;
	int ntot;
	double prob;
	double ptot;
	Context dist;

	if (n < 0)
	{
		cerr << "\n";
		cerr << "I4VEC_MULTINOMIAL_SAMPLE - Fatal error!\n";
		cerr << "  N < 0\n";
		exit(1);
	}

	if (ncat <= 1)
	{
		cerr << "\n";
		cerr << "I4VEC_MULTINOMIAL_SAMPLE - Fatal error!\n";
		cerr << "  NCAT <= 1\n";
		exit(1);
	}

	for (i = 0; i < ncat - 1; i++)
	{
		if (p[i] < 0.0)
		{
			cerr << "\n";
			cerr << "I4VEC_MULTINOMIAL_SAMPLE - Fatal error!\n";
			cerr << "  Some P(i) < 0.\n";
			exit(1);
		}

		if (1.0 < p[i])
		{
			cerr << "\n";
			cerr << "I4VEC_MULTINOMIAL_SAMPLE - Fatal error!\n";
			cerr << "  Some 1 < P(i).\n";
			exit(1);
		}
	}

	ptot = 0.0;
	for (i = 0; i < ncat - 1; i++)
	{
		ptot = ptot + p[i];
	}

	if (0.9999 < ptot)
	{
		cerr << "\n";
		cerr << "I4VEC_MULTINOMIAL_SAMPLE - Fatal error!\n";
		cerr << "  1.0 < Sum of P().\n";
		exit(1);
	}
	//
	//  Initialize variables.
	//
	ntot = n;
	ptot = 1.0;

	for (i = 0; i < ncat; i++)
	{
		ix[i] = 0;
	}
	//
	//  Generate the observation.
	//


	for (icat = 0; icat < ncat - 1; icat++)
	{
		prob = p[icat] / ptot;
		ix[icat] = dist.binomial_dis(ntot, prob, eng);
		ntot = ntot - ix[icat];
		if (ntot <= 0)
		{
			return;
		}
		ptot = ptot - p[icat];
	}

	ix[ncat - 1] = ntot;

	return;
}
/*
int sur_chlo(int *eggs, double *chlo){
// survival based only on the chlorophyll level
// input: number of eggs produced on a month at the location
// input: chlo is the chlorophyll of the month
// output: number of larvae survived at the end of 30 days
// no estimation- based on the equations in Fabricius et al 2010

int survived;
double sur;
double j;
double a, b;

dist.SetSeed((int)time(0) ^ (omp_get_thread_num() + 1));

a = dist.normal_dis(-0.153, 2.528);
b = dist.normal_dis(2.108, 1.701);
j = exp(a + b*(log2(*chlo)));
sur = j / (1 + j);
survived = dist.binomial_dis(*eggs, sur);
//cout << "Number of Larvae survived due to Chlorophyll level " << survived << endl;
return survived;
}
*/

long long sur_egg(long long *eggs, double lar_mor, base_generator_type& eng){

	long long survived;
	Context dist;
	survived = dist.binomial_dis_long(*eggs, pow((1 - lar_mor), 20), eng);

	return survived;
}

long long sur_chlo_sven(long long *eggs, double *chlo, int temp, double lar_mor, base_generator_type& eng){
	// survival rate is based only on the algae level and water temperature
	// instead of directly modifiy the survival, the effect of algae and water temperature merges with natural mortatlity
	// input: number of eggs produced on a month at the location
	// input: chlo is the chlorophyll of the month
	// input: water temperature, min=28 and max=30
	// output: number of larvae survived at the end of 30 days 
	// no estimation- based on the equations in Fabricius et al 2010

	long long survived;
	int w_temp;
	double algae;
	double l[3];
	double p[4];
	long long *bin = new long long[4];


	algae = (441.1 + 1336.5*(*chlo));
	w_temp = round(temp); //
	if (algae >= 1100){
		if (w_temp <= 28 && w_temp > 26){
			l[0] = -4.704 + 0.000179*algae;
			l[1] = -3.580 + 0.000241*algae;
			l[2] = -2.054 + 0.000268*algae;
			p[0] = exp(l[0]) / (1 + exp(l[0]));
			p[1] = exp(l[1]) / (1 + exp(l[1]));
			p[2] = exp(l[2]) / (1 + exp(l[2]));
			p[3] = 1 - p[2]; // ones not complete settlement after 24 days 
			multi_boost_long(*eggs, p, 4, bin, eng);
			survived = bin[1] * pow((1 - lar_mor), 10) + bin[2] * pow((1 - lar_mor), 17) + bin[3] * pow((1 - lar_mor), 24) + bin[3] * pow((1 - lar_mor), 40);
		}
		else if (w_temp == 29){
			l[0] = -4.704 + 0.000179*algae;
			l[1] = -3.580 + (0.000241 + 0.000232)*algae;
			l[2] = -2.054 + 0.000268*algae;
			p[0] = exp(l[0]) / (1 + exp(l[0]));
			p[1] = exp(l[1]) / (1 + exp(l[1]));
			p[2] = exp(l[2]) / (1 + exp(l[2]));
			p[3] = 1 - p[2]; // ones not complete settlement after 24 days 

			multi_boost_long(*eggs, p, 4, bin, eng);
			survived = bin[1] * pow((1 - lar_mor), 11) + bin[2] * pow((1 - lar_mor), 18) + bin[3] * pow((1 - lar_mor), 25) + bin[3] * pow((1 - lar_mor), 40);

		}
		else if (w_temp == 30 || w_temp == 31){
			l[0] = -4.704 + 0.000179*algae + 1.637;
			l[1] = -3.580 + (0.000241 + 0.000280)*algae;
			l[2] = -2.054 + 0.000268*algae;
			p[0] = exp(l[0]) / (1 + exp(l[0]));
			p[1] = exp(l[1]) / (1 + exp(l[1]));
			p[2] = exp(l[2]) / (1 + exp(l[2]));
			p[3] = 1 - p[2]; // ones not complete settlement after 24 days 
			multi_boost_long(*eggs, p, 4, bin, eng);
			survived = bin[1] * pow((1 - lar_mor), 11) + bin[2] * pow((1 - lar_mor), 18) + bin[3] * pow((1 - lar_mor), 25) + bin[3] * pow((1 - lar_mor), 40);
		}
		else {
			survived = 0;
		}

	}
	else{ // low nutrient no development
		if (w_temp > 32 && w_temp < 26){
			survived = 0;
		}
		else{ // low nutrient 
			survived = *eggs*(0.001)*pow((1 - lar_mor), 60);
		}
	}

	//cout << "Number of Larvae survived due to Chlorophyll level " << survived << endl;
	return survived;
}
long long nat_mor(long long *eggs, double sur, base_generator_type& eng){
	// Natural mortality of Larvae
	// Constant over time
	// input: eggs survived on month $t$ at reef $S$ after exposed to chlorophyll
	// output: Number of larvae survived after 30 days at reef $S$ after 
	// improvement: estimate mortality

	Context dist;


	long long out;
	if (sur < 0 || *eggs<0 ){
		cout << "problem nat_mor" << endl; 
		exit(0);
	}
	long long  mean = sur*(*eggs); 
	long long std = sqrt((1-sur)*sur*(*eggs));

	if (mean < 5000){
		out = dist.binomial_dis(*eggs, sur, eng);
		//cout << "bino "<< out << endl;
	}
	else{
		out = dist.normal_dis_long(mean, std, eng);
		//cout << out << " " << log(mean) << " " << log(var);
		//out = exp(out);
		//cout << "norm " << *eggs-out << endl;

	}
	if (out < 0){
		cout << out << endl;
		cout << "problem nat_mor" << endl;
		exit(0);
	}
	return out;
}
double degreeToRadian(double x)
{
	double degToRadFactor = PI / 180;
	return x*degToRadFactor;

}
double distance(double x1, double y1, double x2, double y2)
{// Haversine
	//km
	double a, c;
	double la1, la2, lo1, lo2, la_dif, lo_dif;
	double dis;
	la_dif = degreeToRadian(x2 - x1);
	lo_dif = degreeToRadian(y2 - y1);
	la1 = degreeToRadian(x1);
	la2 = degreeToRadian(x2);
	lo1 = degreeToRadian(y1);
	lo2 = degreeToRadian(y1);

	a = sin(la_dif / 2)*sin(la_dif / 2) + cos(la2)*cos(la1)*sin(lo_dif / 2)*sin(lo_dif / 2);
	c = 2 * atan2(sqrt(a), sqrt(1 - a));
	dis = 6371.01*c * 1000;
	return dis;
}

/*
void T_recru(struct reef_month *reef, int num_reef)
{
// temporary one- based on the distance between reef
// looping through the locations
int i, j;
double *dis = new double[num_reef];
double *prob = new double[num_reef];
int *lar = new int[num_reef];
int eggs;

for (i = 0; i < num_reef; i++){
double t_s = 0;
eggs = reef[i].eggs_out*(1 - T_12 - SELF_REC); // eggs is the number of eggs remain in the system
for (j = 0; j < num_reef; j++){
if (i != j){
dis[j] = 10000 / distance(reef[i].lat, reef[i].longi, reef[j].lat, reef[j].longi);
t_s += dis[j];
}
}

dis[i] = t_s*(SELF_REC / T_12); // allocation for sel-recruits
dirichlet_sample_v1(num_reef, dis, prob);
i4vec_multinomial_sample_v1(eggs, prob, num_reef, lar);
for (j = 0; j < num_reef; j++){
reef[j].larv_in += lar[j];
}
}
delete[] dis;
delete[] prob;
delete[] lar;
dis = nullptr;
prob = nullptr;
lar = nullptr;
return;
}
*/

void rec_2para(struct reef_month *reef, double rec[], double disp[],  int num_reef, base_generator_type& eng){
// two source of recruitment, self recruitment and recieving
// treating the recieveing from global, ie. sum of eggs recieved 
	// disp is actually proportion of (f)eggs produced by other reefs landed in this reef

	Context dist;
	int i;
	long long *disp_eggs = new long long [num_reef];
	long long tot=0;
	long long temp;

	for (i = 0; i < num_reef; i++){
		// self recruit
		if (reef[i].eggs_out != 0){
			reef[i].larv_in += dist.binomial_dis_long(reef[i].eggs_out, rec[i], eng);
			disp_eggs[i] = reef[i].eggs_out - reef[i].larv_in;
			tot += disp_eggs[i]; 
		}
	}
	for (i = 0; i < num_reef; i++){
		if (tot > 0){
			temp = dist.binomial_dis_long(tot - disp_eggs[i], disp[i], eng);
			reef[i].larv_in += temp;
			tot -= temp;
		}
	}
	return;
}
void rec_1para(struct reef_month *reef, double rec[], int num_reef, base_generator_type& eng){
	// proprotion recieved from the pool of eggs

	Context dist;
	int i;
	long long tot = 0;
	long long temp;

	for (i = 0; i < num_reef; i++){
		// total eggs in the system
		tot += reef[i].eggs_out;
	}
	for (i = 0; i < num_reef; i++){
		if (tot > 0){
			temp = dist.binomial_dis_long(tot, rec[i], eng);
			reef[i].larv_in = temp;
			tot -= temp;
		}
	}

	return;
}

void simple_rec(struct reef_month *reef, int num_reef, base_generator_type& eng)
{
	Context dist;
	int i, j;
	long long diff;

	for (i = 0; i < num_reef; i++){
		//sort out self-recruitment
		//cout << reef[i].REEF_ID << " LOC " << reef[i].eggs_out << " " << reef[i].larv_in << endl;
		if (reef[i].eggs_out != 0){
			reef[i].larv_in += dist.binomial_dis_long(reef[i].eggs_out, SELF_REC, eng);
			diff = reef[i].eggs_out - reef[i].larv_in;
			//cout << reef[i].REEF_ID << " " << reef[i].larv_in << " "  <<diff  << endl; 
			if (diff != 0){
				for (j = 0; j < num_reef; j++){
					//sort out self-recruitment
					if (i != j){
						reef[j].larv_in += dist.binomial_dis_long(diff, DISPER, eng);
						//cout << reef[j].REEF_ID << " " << reef[j].larv_in << endl;
					}
				}
			}
		}
	}
	return;
}
void recru_con(struct reef_month *reef, int num_reef, struct con *connect, int n_conn, base_generator_type& eng)
{
	// temporary one- based on the connectivity between reef
	// connect is the struct storing connectivity
	// looping through the locations
	// int n_conn is the length of connection records


	double *prob = new double[num_reef];
	string *sink = new string[num_reef];
	long long *lar = new long long[num_reef];
	int i, j, k, l;

	long long eggs;
	double sum;
	Context dist;

	for (i = 0; i < num_reef; i++){
		//sort out self-recruitment
		if (reef[i].eggs_out != 0){
			reef[i].larv_in += dist.binomial_dis_long(reef[i].eggs_out, SELF_REC, eng);
			eggs = reef[i].eggs_out - reef[i].larv_in; // eggs is the number of eggs remain in the system
			// find the connectivity of the current year; 

			if (eggs >1){
				sum = 0;
				k = 0;

				for (j = 0; j < n_conn; j++){
					if (reef[i].yr == connect[i].year[j] && reef[i].REEF_ID == connect[i].REEF_ID){
						prob[k] = connect[i].connect[j];
						sink[k] = connect[i].to[j];
						sum += prob[k];
						k++;
					}
				}

				if (k == 0){
					cout << "No connectivity data" << endl;
					exit(0);
				}
				if (k > num_reef){
					cout << "Too many connectivity data " << endl;
					exit(0);
				}

				prob[num_reef - 1] = 1 - sum; // last array stores the probability not recruit to 12 reefs

				for (l = 0; l < num_reef; l++){
					if (prob[l] < 0.000000000){
						cout << "hereh" << endl;
					}
				}
				multi_boost_long(eggs, prob, num_reef, lar, eng);
				// assign back to mth_reef

				for (j = 0; j < num_reef - 1; j++){
					for (k = 0; k < num_reef; k++){
						if (reef[k].REEF_ID == sink[j]){
							reef[k].larv_in += lar[j];
						}
					}
				}
			}
		}
	}
	return;
}
/*
void recru(int *eggs_num, struct parameters *param, int num_reef, double *con, int *out, int *abund2p){
/* recruitment to population, i.e 100% success of settlement
input:
1. Number of eggs survived from Chlorophyll and natural death
2. Parameters
3. num_reef: number of reefs in the system
4. abund2p: number of 2+ adult in each reef, array[num_reef]
5. con: array of connectivity from current reef to all other reefs (including itself) in the system and reef outside the systems
Output:
*out- number of larvae disperse from current reef to all reef + lost in space
*/
/*	double *V = new double[num_reef + 1];// store V
	double *r = new double[num_reef + 1]; // store r
	double *prob = new double[num_reef + 1]; // store p
	double tot = 0;
	double temp;
	int i;

	for (i = 0; i < num_reef; i++){
	temp = abund2p[i] / (*eggs_num*(1 + pow(abund2p[i], param->b)));
	V[i] = param->alpha*temp + (1 - param->alpha)*con[i];
	tot += temp;
	}
	if (tot > *eggs_num){
	cout << "Number of recruits is larger than number of eggs" << endl;
	exit(1);
	}
	V[num_reef + 1] = param->alpha*(1 - tot / (*eggs_num)) + (1 - param->alpha)*con[num_reef + 1];

	for (i = 0; i < num_reef; i++){
	r[i] = r8_normal_sample(V[i], param->sig_r);
	}

	dirichlet_sample_v1(num_reef, r, prob);
	i4vec_multinomial_sample_v1(*eggs_num, prob, num_reef, out);
	// next step assign out to reef struct
	delete[] V;
	delete[] r;
	delete[] prob;
	V = nullptr;
	r = nullptr;
	prob = nullptr;
	return;
	}
	*/
double trans_length(struct chlo_length *dat, int no_rec, double *chlo) {
	// length estimate from larvae to settlment
	// predicted length based on the chlorophyll level of the month
	// ATM- Deterministic, because the original data is fitted with "loess" function
	// Fabricius et al 2010
	// Because the monthly chlorophyll level is constant, all larvae produced in that month is assumed to have the same starting length
	// input cholorophyll and length datat, number of records in the data and chlorophyll level
	// output length (mm) - 2 decimal places
	double p_length = 0;
	int i;
	double min = 10, max = 0;
	// round to 2 decima
	double chlo_s = roundf(*chlo * 100) / 100;

	for (i = 0; i < no_rec; i++){
		if (dat[i].chro < min){ min = dat[i].chro; }
		if (dat[i].chro > max){ max = dat[i].chro; }
	}

	if (*chlo < min){
		p_length = 0.85;
	}
	else if (*chlo > max){
		p_length = 1.3;
	}
	else{
		for (i = 0; i < no_rec; i++){
			if (fabs(chlo_s - dat[i].chro) < 0.000001){// compare two double
				p_length = dat[i].p_len;
			}
		}
	}
	return p_length;
}

double growth_J(int in, base_generator_type& eng){
	//growth of juvenile COTS, i.e age <6 month old
	//Stage of eating coralline algae
	// assume growth is uniform from 1.5mm to 2.6 mm per month
	// reference Pratchett 2014
	// input: current length of creature (in)
	// output: Length of creature after a month of growth
	double out;
	//GROWTH_J_LO, GROWTh_J_HI 
	Context dist;
	out = in + dist.uniform_real_dis(GROWTH_J_LO, GROWTh_J_HI, eng);
	return out;
}
/*int growth_A(int den, double food, int age, struct parameters *param, int *in){
// Growth of COTs 6+ month old
// when diet changes from algae to coral
// growth rate is function
// abundance of density of COTs- number of COTs/area of reef (den)
// coral cover- only the live and corals in the "preferred group (food)
// age- age of COTs in month
// Input: three listed above + parameters + current length
// Output: length at the end of month

double out;
double lambda;

lambda = r8_normal_sample(param->b_0*(age) + param->b_1 *(den / food), param->sig_l);
out = *in *(1 + lambda*exp(-lambda));
return out;

}*/

/*int growth_A_post(int length){
	// post sexula matruity
	// growth in negligible- mean 0, stdev of 2
	int out;
	out = length + r8_normal_sample(0, 2);
	return out;
	}*/
/*double feed_rate(double den, double food, double k, double sigma, base_generator_type& eng){
	double out;
	double mean;
	Context dist;

	mean = food / den;
	out = dist.normal_dis(mean, sigma, eng);

	return out;
	}*/
double growth_A(double val, int length, double beta, double sigma, base_generator_type& eng){
	// Growth of COTs 6+ month old to 2 years
	// when diet changes from algae to coral
	// growth rate is function 
	// food: relative coral cover (from the first year)
	// density: relative density from 1st year
	// length - length of cots in month
	// Input: three listed above + parameters 
	// Output: length at the end of month

	double out;
	double k;
	Context dist;
	if (length <= 482){
		//	k = 0.5744 + 0.9077*log(length);
		k = 0.5744 + beta * val * log(length);
		//k = beta * (food*size / den) + 0.9077*log(length);
		out = dist.normal_dis(exp(k), sigma, eng);
	}
	else{
		out = length + dist.normal_dis(5, 2, eng);
	}
	return out;

}

double growth_A_withfeed(int length, double beta, double sigma, double food, base_generator_type& eng){
	// Growth of COTs 6+ month old to 2 years
	// when diet changes from algae to coral
	// growth rate is a step function, depending on FOOD availability
	//

	double out;
	double k;
	Context dist;

	if (food > 1){
		if (length <= 482){
			//	k = 0.5744 + 0.9077*log(length);
			k = 0.5744 + 0.9077 * log(length);
			//k = beta * (food*size / den) + 0.9077*log(length);
			out = dist.normal_dis(exp(k), 5, eng);
		}
		else{
			out = length + dist.normal_dis(5, 2, eng);
		}
	}
	else{
		k = 0.5744 + beta * log(length);
		out = dist.normal_dis(exp(k), sigma, eng);
	}
	return out;
}

double growth_A_withfeed_function(int length, double sigma, double *alpha, base_generator_type& eng){
	// Growth of COTs 6+ month old to 2 years
	// when diet changes from algae to coral
	// growth rate is a step function, depending on FOOD availability
	//
	double out;
	double k;
	Context dist;

	
	if (*alpha >= 0.907){ // coral cover increase
		if (length <= 482){
			//	k = 0.5744 + 0.9077*log(length);
			k = 0.5744 + 0.9077 * log(length);
			//k = beta * (food*size / den) + 0.9077*log(length);
			out = dist.normal_dis(exp(k), sigma, eng);
		}
		else{
			out = length + dist.normal_dis(5, 2, eng);
		}
	}
	else if (*alpha <= 0.87){// close to no growth if the food is drop by 5 % per month
		k = 0.5744 + 0.87 * log(length);
		//k = beta * (food*size / den) + 0.9077*log(length);
		out = dist.normal_dis(exp(k), sigma, eng);
	}
	else{// corval cover reduces
		k = 0.5744 + (*alpha)* log(length);
		out = dist.normal_dis(exp(k), sigma, eng);
	}
	return out;
}
/*double nat_mor_J(double *abun){
// Natural mortality for COT 6+mth to less than 24 months
// Based on MCCallum (1990)
// input:
// abundance of COTs on reef
double mor;
mor = NMJ_d + NMJ_a*NMJ_P / (*abun + NMJ_a / NMJ_s) + NMJ_K*(*abun);
return mor;
}*/
int cot_mor(double length, double min_sur, double max_sur, double k, double x0, base_generator_type& eng){

	// length based survival
	// P(Survival at the end of month) = (1-exp(-lambda* length))
	// Starting value of lambda is 0.1
	// input length, ouput survival
	// change MOR_LAMBDA to param->S_l later 

	double p_sur;
	int fate;
	Context dist;
	p_sur = min_sur + (max_sur / (1 + exp(-k*(length - x0))));
	//#pragma omp critical
	//cout << length << " " << p_sur << " " << k << " " << x0 << endl;
	if (p_sur <= 0){
		cout << length << endl;
		cout << p_sur << endl;
		exit(EXIT_FAILURE);
	}
	if (p_sur >= 0){
		fate = dist.bernoullie_dis(p_sur, eng);
	}
	// 1:alive;; 0- dead
	return fate;
}
long long relative_fec_op1(int length){
	// number of eggs produced by female
	// input: length of females
	long long f_eggs;
	double ln_wt;
	ln_wt = (-9.673964 + 2.929*log(length));
	f_eggs = exp(6.324359 + 1.439*ln_wt);
	// option 1: random draw from 40% to 100% of eggs are feterlized
	f_eggs = FER_DIS*f_eggs;

	return f_eggs;
}

void ini_length_struct(int *(&out), int pop_size, struct ini_length *in, int bin, double &ini_det, base_generator_type& rng)
{
	//draw length for individual of the initial population
	//Based on the 1986 cohort of Fisk(1992)
	//pop_size: population size
	//in: input, initial length structure
	//bin: number of bins in the initial length structure
	int i, ctr = 0, j;
	int *K = new int[bin]; // K stores numbers in each bin
	double *prob = new double[bin];
	Context dist;

	int t1;

	for (i = 0; i < bin; i++){
		prob[i] = in[i].prop;
	}
	multi_boost(pop_size, prob, bin, K, rng);
	for (i = 0; i < bin; i++){
		for (j = 0; j < K[i]; j++){
			do{
				t1 = dist.normal_dis(in[i].mean, in[i].std, rng);
			} while (t1<0);
			out[ctr] = t1; 
			if (t1 >= 300){
				ini_det++;
			}
			ctr++;
		}
	}
	delete[] K;
	delete[] prob;
	K = nullptr;
	prob = nullptr;
	return;
}

int length_age(int length, base_generator_type& eng){
	// predict age based on a given length
	int age;
	Context dist;

	if (length <= 30) age = dist.uniform_real_dis(6, 12, eng);
	else if (length > 30 && length <= 100) age = dist.uniform_real_dis(13, 24, eng);
	else if (length > 100 && length <= 250) age = dist.uniform_real_dis(25, 36, eng);
	else if (length > 250 && length <= 300) age = dist.uniform_real_dis(37, 48, eng);
	else if (length > 300 && length <= 400) age = dist.uniform_real_dis(49, 60, eng);
	else { age = dist.uniform_real_dis(61, MAX_ADULT_AGE, eng); }
	return age;
}

void ini_pop_distribution(struct reef_month *(&in), int no_loc, int N0) {
	// function: distribution the initial population according to weight of each reef
	// Dirichlet-multinomial 
	// input: weight for each reef and number of reefs 
	// output: number of COTS allocated to each reef

	base_generator_type rng(42);
	double *prop = new double[no_loc];
	int *n = new int[no_loc];
	int i;
	int ctr = 0;

	for (i = 0; i < no_loc; i++){
		prop[i] = in[i].rel;
	}

	multi_boost(N0, prop, no_loc, n, rng);

	for (i = 0; i < no_loc; i++)
	{
		in[i].N0 = n[i];
	}
	delete[] prop;
	//delete[] weight;
	delete[] n;
	prop = nullptr;
	//weight = nullptr;
	n = nullptr;
	return;
}

int convert_yr_mth_time(int s_yr, int s_mth, const int yr, const int mth){
	// convert year and month to simulation day
	// Why: Coral coverage and chlorophyll are input in year and month, here convert the Year_mth to sim_days
	// input: Year and month (array)
	// output: day

	int s_day;
	s_day = (yr * 12 + mth) - (s_yr * 12 + s_mth);
	return s_day;
}

void convert_coralcover_foo(struct coral_cov *in, struct reef_month *out, int no_loc, int sim_length, int size_cc_rec){
	// convert coral cover format to array of linkedlist

	int i, j, k;

	for (i = 0; i < no_loc; i++){
		for (j = 0; j < sim_length + 2; j++){
			out[i].s_day[j] = j;
		}
	}

	for (k = 0; k < size_cc_rec; k++){
		for (i = 0; i < no_loc; i++){
			for (j = 0; j < sim_length + 2; j++){
				if (in[k].REEF_ID == out[i].REEF_ID && in[k].s_day == out[i].s_day[j]){
					out[i].cov[j] = in[k].cov;
				}
			}
		}
	}
	return;
}

void format_chloro(struct chloro_ori *in, struct chlorophyll *out, int no_loc, int sim_length, int chloro_loc){
	// convert coral cover format to array of linkedlist

	int i, j, k, l;
	int no_yrs;
	no_yrs = sim_length / 12;

	for (k = 0; k < no_loc; k++){
		l = 0;
		for (i = 0; i < no_yrs; i++){
			for (j = 1; j <= 3; j++){ // Jan, Feb, March
				out[k].s_day[l] = i * 12 + j;
				l++;
			}
		}

	}

	for (k = 0; k < chloro_loc; k++){
		for (i = 0; i < no_loc; i++){
			for (j = 0; j < sim_length; j++){
				if (in[k].REEF_ID == out[i].REEF_ID && in[k].s_day == out[i].s_day[j]){
					out[i].chlo[j] = in[k].chlo;
				}
			}
		}
	}
	return;
}

void fomat_connect(struct connect_ori *in, struct con *(&out), int no_loc, int n_ori, int s_yrs){

	// s_yrs; number of simulation years
	// n_ori: length of the in records

	int i, j, k;

	for (i = 0; i < no_loc; i++){
		for (j = 0; j < n_ori; j++){
			if (in[j].from == out[i].REEF_ID){
				for (k = 0; k < (no_loc - 1)*(s_yrs + 1); k++){
					if (in[j].year == out[i].year[k] && in[j].to == out[i].to[k]){
						out[i].connect[k] = in[j].connect;
					}
				}
			}
		}
	}
	return;
}

void struct2mat(struct parameters *in, MatrixXd &out, int n_row){
	// convert params to matrix
	// n_row is the particle size
	int j;
	// w[1] is removed, because w[1]=1-w[0]
	for (j = 0; j < n_row; j++){
		out(j, 0) = in[j].alpha[0];
		out(j, 1) = in[j].alpha[1];
		out(j, 2) = in[j].alpha[2];
		out(j, 3) = in[j].beta;
		out(j, 4) = in[j].beta_sigma;
		out(j, 5) = in[j].detect[0];
		out(j, 6) = in[j].detect[1];
		out(j, 7) = in[j].detect[2];
		out(j, 8) = in[j].error_N;
		out(j, 9) = in[j].jv_sur;
		out(j, 10) = in[j].k;
		out(j, 11) = in[j].larv_mor;
		out(j, 12) = in[j].x0;
	}
	return;
}

void struct2array(struct parameters *in, VectorXd &out, int t){
	// convert params to a single array at give iteration, t
	out(0) = in[t].alpha[0];
	out(1) = in[t].alpha[1];
	out(2) = in[t].alpha[2];
	out(3) = in[t].beta;
	out(4) = in[t].beta_sigma;
	out(5) = in[t].detect[0];
	out(6) = in[t].detect[1];
	out(7) = in[t].detect[2];
	out(8) = in[t].error_N;
	out(9) = in[t].jv_sur;
	out(10) = in[t].k;
	out(11) = in[t].larv_mor;
	out(12) = in[t].x0;
	return;
}

void array2struct(struct parameters *&out, MatrixXd in, int t){
	// assign the array into params at the position t
	// t is the sequence of "PARTICLE"
	out[t].alpha[0] = in(0);
	out[t].alpha[1] = in(1);
	out[t].alpha[2] = in(2);
	out[t].beta = in(3);
	out[t].beta_sigma = in(4);
	out[t].detect[0]= in(5);
	out[t].detect[1] = in(6);
	out[t].detect[2] = in(7);
	out[t].error_N = in(8);
	out[t].jv_sur = in(9);
	out[t].k = in(10);
	out[t].larv_mor = in(11);
	out[t].x0 = in(12);
	return;
}
//***********************************************
/*
double *r8mat_copy_new(int m, int n, double a1[])

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_COPY_NEW copies one R8MAT to a "new" R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8's, which
//    may be stored as a vector in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 July 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A1[M*N], the matrix to be copied.
//
//    Output, double R8MAT_COPY_NEW[M*N], the copy of A1.
//
{
double *a2;
int i;
int j;

a2 = new double[m*n];

for (j = 0; j < n; j++)
{
for (i = 0; i < m; i++)
{
a2[i + j*m] = a1[i + j*m];
}
}
return a2;
}
*/
int check_array(MatrixXd a, int dim){

	//check if the parameter is less than zero
	int i, jj = 0;

	
	for (i = 0; i < dim; i++){
		if (a(i) < 0) {
			jj = 0;
			break;
		}
		else{
			jj = 1;
		}
	}

	if (a(0)>1 || a(1)>1 || a(2)>1 || a(9) > 1 || a(11) > 1){// checking larval and juv mortality
		jj = 0;
	}
	// check if detect is above1
	return jj;
}

double multivariate_normal_ln_pdf(VectorXd mean, MatrixXd cov, VectorXd val){

	double pi = 3.141592653589793;
	int k = mean.rows();
	double ln_pdf;


	VectorXd dif = val - mean;

	ln_pdf = (-0.5) *(k*log(2 * pi) + log(cov.determinant()) + dif.transpose()* cov.inverse() * dif);
	return ln_pdf;
}
double normal_pdf(double av, double sd, double rval)

//****************************************************************************80
//
//  Purpose:
//
//    R8_NORMAL_PDF evaluates the PDF of a normal distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN90 version by Guannan Zhang.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double AV, the mean value.
//
//    Input, double SD, the standard deviation.
//    0.0 < SD.
//
//    Input, double RVAL, the point where the PDF is evaluated.
//
//    Output, double R8_NORMAL_PDF, the value of the PDF at RVAL.
//
{
	double pi = 3.141592653589793;
	double rtemp;
	double value;

	if (sd <= 0.0)
	{
		cerr << "\n";
		cerr << "R8_NORMAL_PDF - Fatal error!\n";
		cerr << "  Standard deviation must be positive.\n";
		exit(1);
	}

	rtemp = (rval - av) * (rval - av) * 0.5 / (sd * sd);

	value = exp(-rtemp) / sd / sqrt(2.0 * pi);

	return value;
}

double lognormal_pdf_ln(double av, double sd, double rval){

	double ln_pdf;

	if (rval < 0){
		cerr << "\n";
		cerr << "Problem with lognormal, rval is less than 0" << endl;
		cout << "HHHHHHHHH " << rval << endl;
		exit(1);
	}

	ln_pdf = -(log(sd) + log(rval) + 0.9189385) - (log(rval) - av)*(log(rval) - av) / (2 * sd*sd);

	return ln_pdf; 

}

double unif_pdf_ln(double min, double max, double rval){

	double ln_pdf;

	if (rval < min || rval > max ){
		ln_pdf = 1;
	}
	else{
		ln_pdf = -log(max-min);
	}
	return ln_pdf;

}
double summary_stat_sd(struct reef_month *loc, int no_loc, struct reef_annual *in, int no_in){
	// loc: location reference
	// no_loc: number of reefs
	// in: struct store both predicted and observed values
	// no_in: number of in records
	// where records are NA, the value is "-9999"
	int i, j, k;
	double summary_stat = 0;

	for (i = 0; i < no_loc; i++)//location
	{
		// step 1: count number of relevant records at location [i]
		int N_rec = 0;

		for (j = 0; j < no_in; j++){
			//if (in[j].REEF_ID == loc[i].REEF_ID && in[j].observed >= 0){// all data
			if (in[j].REEF_ID == loc[i].REEF_ID && in[j].observed > 0){// does no include zero
				N_rec++;
			}
		}
		VectorXd diff(N_rec);

		k = 0;
		for (j = 0; j < no_in; j++){
			//if (in[j].REEF_ID == loc[i].REEF_ID && in[j].observed>=0){ // all data
			if (in[j].REEF_ID == loc[i].REEF_ID && in[j].observed > 0){
				// relative abudnace -
				diff(k) = in[j].pred_pop_size - in[j].observed;
				k++;
				//cout << in[j].year << endl;
			}
		}
		
		/*VectorXd x_bar = diff.array() - diff.mean();
		double temp = x_bar.adjoint()*x_bar;
		double sd = sqrt(temp / (double)(N_rec - 1));*/
		double temp = diff.adjoint()*diff;
		summary_stat += sqrt(temp /N_rec);

	}
	return summary_stat;
}

double zero_inf_p_pdf(double pz, int mean, int val){
	// calculate the pdf for the discrenpancy value
	double out;
	if (val == 0){
		out = pz + (1 - pz)*exp(-mean);
	}
	else{
		poisson dd(mean);
		out = pdf(dd, val)*(1 - pz);
	}
	return out; 
}


vector<size_t> sort_indexes(const vector<double> &v) {
	// really cool function- sorting out vector and keep the index
	// sorting from largest to smallest value
	// need to be used with
	/*	for (auto i : sort_indexes(rho)) {
		//params[i].rank = ctr; // to assign "rank" to the correct location 
		cout << i << " " << rho[i] << " " << params[i].rank << endl; // this is just to show values from largest to smallest; and the location of value (i) in the vector  
		ctr++;
	}*/
	// initialize original index locations
	vector<size_t> idx(v.size());
	iota(idx.begin(), idx.end(), 0);

	// sort indexes based on comparing values in v
	sort(idx.begin(), idx.end(),
		[&v](size_t i1, size_t i2) {return v[i1] > v[i2]; });

	return idx;
}

int check_ini_para(struct parameters *in, int p, int no_para, int no_reef){
	// n:number of reef
	// output of zero means need to redraw
	int i;
	int flag = 1;

	for (i = 0; i< no_reef; i++){
		if (in[p].alpha[i] <0 || in[p].alpha[i] >1 || in[p].detect[i] < 0){// checking larval and juv mortality
			flag = 0;
		}
	}

	if (in[p].beta < 0 || in[p].beta_sigma < 0 || in[p].error_N < 0 || in[p].error_N < 0 || in[p].jv_sur<0 || in[p].jv_sur>1 || in[p].k < 0 || in[p].larv_mor <0 || in[p].larv_mor>1)
	{
		flag = 0;
	}

	return flag;

}
