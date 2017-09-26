#ifndef structs_H
#define structs_H
// define all structs here

#define NUMVARPHI 2 //need to change this if we use mixture model
#include <string>
#include <vector>

struct adult{
	int age; // 1 months to max 8 year 
	int sex;  //1:female and 0: male
	double size; // mm
	int spawn; // indictor if the animal has spwaned, females only // 1: yes, 0: no and updated to 0 in March
	int bs_ind=0; // Body size indicator, if an animal's bosy size reduced in 10 consective month, it dies
	adult *next;
};

struct reef_month
{
	// dynamic, change monthly
	// length of the linklist is 13, number of reef
	// 
	int yr; 
	int mth;
	std::string REEF_ID;
	double lat; // latitude of reef
	double longi; // longitude of reef
	double rel=0; // initial relative abundance- weight of dirichlet 
	int ini_obs; // 1: no COTS were found in the first two years; 0: COTS were found on the first year
	std::string shelf; // I: inner M: Middle and O: outer
	int N0;
	double N0_det=0;
	//double size; // size of reef, KM^2
	long long eggs_out = 0; //number of eggs produced by the COTs on the reef that month
	long long larv_in = 0; // number of larvae recieved by the reef
	//int N6p[500]; // population of previous time point
	int s_day[1000]; // simulation day
	double cov[1000]; // coarl coverage corresponding to s_day
};

struct reef_annual{
	// length of the linklist is 13, number of reef
	int year;
	std::string REEF_ID;
	int mth;
	double pred_pop_size; // population size at Nov, before recruitment
	double observed=-99.0;
	int area_sampled;
};

struct parameters
{
	double larv_mor; // Larvae daily natural mortality
    double beta; //current (11/05) pre-sexual maturity growth scaling para
	double beta_sigma; // sigma for beta
	double detect[3]; // location specific
	double k;//current (11/05) survival -rate of increase
	double x0; // survival 
	double jv_sur; // juvenille monthly survival
	double alpha[3]; // recieveing parameter
	double error_N;
	int rank; // rank of the summary statistics - changes at each iteration
};

struct chlo_length{
	double chro;
	double p_len;
};

struct ini_length{
	double prop; // lower bound of class 
	double mean; // upper bound of class
	double std; // proporyion
};

struct coral_cov{
	std::string REEF_ID;
	int s_day;
	double cov; // coverage
};

struct chloro_ori{
	// store month chlorophyll
	// original
	std::string REEF_ID;
	int s_day;
	double chlo;
};
struct chlorophyll{
	// store monthly chlorophyll
	std::string REEF_ID;
	int s_day[1000];
	double chlo[1000];
};

struct connect_ori{
	// store orginal connectivity 
	// will be converted to usable format
	// connectivity is annual- not monthly
	std::string from; // sourcing reef
	int year; // simulation days
	std::string to; // sinking reefs
	double connect; // probability of connection
};

struct con{
	// store orginal connectivity 
	// will be converted to usable format
	// connectivity is annually, not monthly
	std::string REEF_ID; // sourcing reef
	int year[1008]; // simulation year, max capacity- 83 years
	std::string to[1008]; // sinking reefs
	double connect[1008]; // probability of connection
};

struct obs_abund{
	// observed data
	int year;
	int mth;
	std::string reef_ID;
	double CPUE; //mean poisson
	double pz; //probability of zeros
};
/*struct coral_comp{
	std::string reef_ID;
	double class1;
};*/

struct reef_size{
	std::string reef_ID;
	double size;
};

struct mth_out{
	int year;
	int mth;
	std::string reef_ID;
	double pop_size;
	int no_female;
	int no_male;
	int l15;
	int l1525;
	int l2540;
	int l40p;
};
struct w_temp{
	// no records on reef level
	int sim;
	int wTemp;
};
#endif
