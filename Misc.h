#include "Eigen/Dense"
#include "Eigen/LU"
# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>
#include <omp.h>
#include "Dist.h"
#include "DefineParameters.h"
#include "struct.h"
using namespace Eigen;
using namespace std;
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/uuid/seed_rng.hpp>


typedef boost::random::lagged_fibonacci607 base_generator_type;
using boost::math::poisson;

void dirichlet_sample(int n, double a[], double *(&out), base_generator_type& eng);
void multi_boost(int n, double p[], int ncat, int *(&ix), base_generator_type& eng);
void multi_boost_1(int n, long double p[], int ncat, int *(&ix), base_generator_type& eng);
void multi_boost_long(long long n, double p[], int ncat, long long *(&ix), base_generator_type& eng);
long long sur_chlo_sven(long long *eggs, double *chlo, int temp, double lar_mor, base_generator_type& eng);
long long nat_mor(long long *eggs, double sur, base_generator_type& eng);
double degreeToRadian(double x);
double distance(double x1, double y1, double x2, double y2);
void recru_con(struct reef_month *reef, int num_reef, struct con *connect, int n_conn, base_generator_type& eng);
double trans_length(struct chlo_length *dat, int no_rec, double *chlo);
double growth_J(int in, base_generator_type& eng);
double growth_A(double val, int length, double beta, double sigma, base_generator_type& eng);
//double feed_rate(double den, double food, double k, double sigma, base_generator_type& eng);
double growth_A_withfeed(int length, double beta, double sigma, double food, base_generator_type& eng);
double growth_A_withfeed_function(int length, double sigma, double *alpha, base_generator_type& eng);
int cot_mor(double length, double min_sur, double max_sur, double k, double x0, base_generator_type& eng);
long long relative_fec_op1(int length);
void ini_length_struct(int *(&out), int pop_size, struct ini_length *in, int bin, double &ini_det, base_generator_type& rng);
int length_age(int length, base_generator_type& eng);
void ini_pop_distribution(struct reef_month *(&in), int no_loc, int N0);
int convert_yr_mth_time(int s_yr, int s_mth, const int yr, const int mth);
void convert_coralcover_foo(struct coral_cov *in, struct reef_month *out, int no_loc, int sim_length, int size_cc_rec);
void format_chloro(struct chloro_ori *in, struct chlorophyll *out, int no_loc, int sim_length, int chloro_loc);
void fomat_connect(struct connect_ori *in, struct con *(&out), int no_loc, int n_ori, int s_yrs);
void struct2mat(struct parameters *in, MatrixXd &out, int n_row);
void struct2array(struct parameters *in, VectorXd &out, int t);
void array2struct(struct parameters *&out, MatrixXd in, int t);
int check_array(MatrixXd a, int dim);
double multivariate_normal_ln_pdf(VectorXd mean, MatrixXd cov, VectorXd val);
double normal_pdf(double av, double sd, double rval);
void simple_rec(struct reef_month *reef, int num_reef, base_generator_type& eng);
long long sur_egg(long long *eggs, double lar_mor, base_generator_type& eng);
double summary_stat_sd(struct reef_month *loc, int no_loc, struct reef_annual *in, int no_in);
double lognormal_pdf_ln(double av, double sd, double rval);
double unif_pdf_ln(double min, double max, double rval);
void rec_2para(struct reef_month *reef, double rec[], double disp[], int num_reef, base_generator_type& eng);
void rec_1para(struct reef_month *reef, double rec[], int num_reef, base_generator_type& eng);
double zero_inf_p_pdf(double pz, int mean, int val);
vector<size_t> sort_indexes(const vector<double> &v);
int check_ini_para(struct parameters *in, int p, int no_para, int no_reef);