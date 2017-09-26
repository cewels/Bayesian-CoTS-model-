// starting month is fixed to dec

#include<stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <omp.h>
#include "ANIMALS.h"
#include "Misc.h"
#include "DefineParameters.h"
#include "DefinePriors.h"
#include "eigenmvn.h"
#include <ctime>
#include "Dist.h"
#include<sys/types.h>
#include <random>
#include <thread>
#include <vector>
#include "ProcessHeader.h"
#include <boost/uuid/seed_rng.hpp>
#include <boost/random/mersenne_twister.hpp>

using namespace std;

//typedef boost::lagged_fibonacci607 base_generator_type;


const int get_int(const string& s) {
	stringstream ss(s);
	int ret;
	ss >> ret;
	return ret;
}

int count_Data(ifstream &filename)// struct weatherdata *enviroData)
{
	int numrec = 0;
	string lines;
	if (filename){
		getline(filename, lines); //get rid of the first line, 
		// ------------ read and store data into enviroData-------------
		while (getline(filename, lines))
		{
			numrec++;
		}

	}
	else{ cout << "Problem in opening file.\n"; exit(0); }

	return numrec;
}

double get_double(const string& s) {
	stringstream ss(s);
	double ret;
	ss >> ret;
	return  ret;
}

/*
void get_chlo_legth_dat(ifstream &filename, struct chlo_length *data){

int i = 0;
string lines;

if (filename){
getline(filename, lines);
i = 0;
while (getline(filename, lines)){
stringstream ss(lines);
vector<std::string> dat;
string field;
while (getline(ss, field, ',')) {
if (field.empty()){
field = "-9999";// missing value
}
dat.push_back(field);
}
data[i].chro = get_double(dat[0]);
data[i].p_len = roundf(get_double(dat[1]) * 100) / 100;;
vector<string>().swap(dat);
i++;
}
}
else{ cout << "Problem in opening CHLOROPHYLL BASED LENGTH file\n"; exit(0); }
}
*/

void get_old_params(ifstream &filename, struct parameters *params, vector<double> &rho){
	// update to mixute on 07/09
	int i = 0;
	string lines;

	if (filename){
		getline(filename, lines);
		i = 0;
		while (getline(filename, lines)){
			stringstream ss(lines);
			vector<std::string> dat;
			string field;
			while (getline(ss, field, '\t')) {
				if (field.empty()){
					cout << "There is missing value in param.txt" << endl;
					exit(0);
				}
				dat.push_back(field);
			}

			rho[i] = get_double(dat[1]);
			params[i].alpha[0] = get_double(dat[2]);
			params[i].alpha[1] = get_double(dat[3]);
			params[i].alpha[2] = get_double(dat[4]);
			params[i].beta = get_double(dat[5]);
			params[i].beta_sigma = get_double(dat[6]);
			params[i].detect[0] = get_double(dat[7]);
			params[i].detect[1] = get_double(dat[8]);
			params[i].detect[2] = get_double(dat[9]);
			params[i].error_N = get_double(dat[10]);
			params[i].jv_sur = get_double(dat[11]);
			params[i].k = get_double(dat[12]);
			params[i].larv_mor = get_double(dat[13]);
			params[i].x0 = get_double(dat[14]);
			vector<string>().swap(dat);
			i++;
		}
	}
	else{ cout << "Problem in opeing params.txt .\n"; exit(0); }
}


void get_ini_leng_dat(ifstream &filename, struct ini_length *ini_leng){
	// update to mixute on 07/09
	int i = 0;
	string lines;

	if (filename){
		getline(filename, lines);
		i = 0;
		while (getline(filename, lines)){
			stringstream ss(lines);
			vector<std::string> dat;
			string field;
			while (getline(ss, field, ' ')) {
				if (field.empty()){
					field = "-9999";// missing value
				}
				dat.push_back(field);
			}
			ini_leng[i].prop = get_double(dat[0]);
			ini_leng[i].mean = get_double(dat[1]);
			ini_leng[i].std = get_double(dat[2]);
			vector<string>().swap(dat);
			i++;
		}
	}
	else{ cout << "Problem in relative frequency among reefs file.\n"; exit(0); }
}
void get_ini_dist(ifstream &filename, struct reef_month *out){

	int i = 0;
	string lines;

	if (filename){
		getline(filename, lines);
		i = 0;
		while (getline(filename, lines)){
			stringstream ss(lines);
			vector<std::string> dat;
			string field;
			while (getline(ss, field, ' ')) {
				if (field.empty()){
					field = "-9999";// missing value
				}
				dat.push_back(field);
			}
			out[i].REEF_ID = dat[0];
			out[i].shelf = dat[1];
			out[i].lat = get_double(dat[2]);
			out[i].longi = get_double(dat[3]);
			out[i].rel = get_double(dat[4]);
			out[i].ini_obs = get_int(dat[5]);
			vector<string>().swap(dat);
			i++;
		}
	}
	else{ cout << "Problem in relative frequency among reefs file.\n"; exit(0); }

}
void get_coral_coverage(ifstream &filename, struct coral_cov *out, int s_yr, int s_mth){

	int i = 0;
	string lines;

	if (filename){
		getline(filename, lines);
		i = 0;
		while (getline(filename, lines)){
			stringstream ss(lines);
			vector<std::string> dat;
			string field;
			while (getline(ss, field, ' ')) {
				if (field.empty()){
					field = "-9999";// missing value
				}
				dat.push_back(field);
			}
			out[i].REEF_ID = dat[0];
			out[i].s_day = convert_yr_mth_time(s_yr, s_mth, get_int(dat[1]), get_int(dat[2]));
			out[i].cov = get_double(dat[3]);
			vector<string>().swap(dat);
			i++;
		}
	}
	else{ cout << "Problem in opening coral coverage.\n"; exit(0); }

}
void get_chlorophyll_dat(ifstream &filename, struct chloro_ori *out, int s_yr, int s_mth){

	int i = 0;
	string lines;

	if (filename){
		getline(filename, lines);
		i = 0;
		while (getline(filename, lines)){
			stringstream ss(lines);
			vector<std::string> dat;
			string field;
			while (getline(ss, field, ' ')) {
				if (field.empty()){
					field = "-9999";// missing value
				}
				dat.push_back(field);
			}
			out[i].REEF_ID = dat[0];
			out[i].s_day = convert_yr_mth_time(s_yr, s_mth, get_int(dat[1]), get_int(dat[2]));
			out[i].chlo = get_double(dat[3]);
			vector<string>().swap(dat);
			i++;
		}
	}
	else{ cout << "Problem in opening chlorophyll data.\n"; exit(0); }

}
void get_connect_dat(ifstream &filename, struct connect_ori *out){

	int i = 0;
	string lines;

	if (filename){
		getline(filename, lines);
		i = 0;
		while (getline(filename, lines)){
			stringstream ss(lines);
			vector<std::string> dat;
			string field;
			while (getline(ss, field, ' ')) {
				if (field.empty()){
					field = "-9999";// missing value
				}
				dat.push_back(field);
			}
			out[i].year = get_int(dat[0]);
			out[i].from = dat[1];
			out[i].to = dat[2];
			out[i].connect = get_double(dat[3]);
			vector<string>().swap(dat);
			i++;
		}
	}
	else{ cout << "Problem in opening connectivity data.\n"; exit(0); }

}
void get_obs_dat(ifstream &filename, struct obs_abund *out){

	int i = 0;
	string lines;

	if (filename){
		getline(filename, lines);
		i = 0;
		while (getline(filename, lines)){
			stringstream ss(lines);
			vector<std::string> dat;
			string field;
			while (getline(ss, field, ' ')) {
				if (field == "NA"){
					field = "-9999";// missing value
				}
				dat.push_back(field);
			}

			out[i].reef_ID = dat[0];
			out[i].year = get_int(dat[1]);
			out[i].mth = get_int(dat[2]);
			out[i].pz = get_double(dat[3]);
			out[i].CPUE = get_double(dat[4]);
			vector<string>().swap(dat);
			i++;
		}
	}
	else{ cout << "Problem in opening COTS abundance data.\n"; exit(0); }

}
void get_reef_size_dat(ifstream &filename, struct reef_size *out){

	int i = 0;
	string lines;

	if (filename){
		getline(filename, lines);
		i = 0;
		while (getline(filename, lines)){
			stringstream ss(lines);
			vector<std::string> dat;
			string field;
			while (getline(ss, field, '\t')) {
				if (field == "NA"){
					field = "-9999";// missing value
				}
				dat.push_back(field);
			}

			out[i].reef_ID = dat[1];
			out[i].size = get_double(dat[5]);

			vector<string>().swap(dat);
			i++;
		}
	}
	else{ cout << "Problem in opening REEF SIZE data.\n"; exit(0); }

}
void get_wtemp_dat(ifstream &filename, struct w_temp *out, int s_yr, int s_mth){

	int i = 0;
	string lines;

	if (filename){
		getline(filename, lines);
		i = 0;
		while (getline(filename, lines)){
			stringstream ss(lines);
			vector<std::string> dat;
			string field;
			while (getline(ss, field, ' ')) {
				if (field == "NA"){
					field = "-9999";// missing value
				}
				dat.push_back(field);
			}

			out[i].sim = convert_yr_mth_time(s_yr, s_mth, get_int(dat[0]), get_int(dat[1]));
			out[i].wTemp = get_int(dat[2]);
			vector<string>().swap(dat);
			i++;
		}
	}
	else{ cout << "Problem in opening WATER temperature data.\n"; exit(0); }

}
void usage()
{
	printf("Command line should have \n 1) simulation period (unit: Years) \n 2) starting population size  \n 3) Initial Sex ratio \n 4) Starting Year \n 5)Minimum tolerance value \n 6)Proportion of particles dropped per iteration \n 7) Number of total particles \n 8) Starting Tolerance value \n 9) Option (0:start new, 1:Continue) \n");
	exit(1);
}


int main(int argc, char **argv)
{
	ifstream aFile, bFile, cFile;
	ofstream myFile;
	int no_chl_rec, bins, bins_zero, coral_cov_rec;
	int i, p;
	int j, k, l;
	int chlo_rec;
	int n_conn; // number of connecticity records
	int no_loc;
	int SIM_LENGTH, INI_POP_SIZE, START_YR, SIM_YEARS;
	double INI_SEX_RAT;
	int DIM;

	double RATE;
	int PARTICLE_SIZE;
	double TOLERANCE_S, START_S;
	int OPTION;
	///---------------------------------------end very temp
	if (argc != 10) usage();
	SIM_YEARS = get_int(argv[1]);
	SIM_LENGTH = SIM_YEARS * 12;
	INI_POP_SIZE = get_int(argv[2]);
	INI_SEX_RAT = get_double(argv[3]);
	START_YR = get_int(argv[4]);
	TOLERANCE_S = get_double(argv[5]);// miniMUM TOLERANCE
	RATE = get_double(argv[6]);// Proportion of samples to be dropped
	PARTICLE_SIZE = get_int(argv[7]);// Number of particles to keep
	START_S = get_double(argv[8]);// miniMUM TOLERANCE
	OPTION = get_int(argv[9]); // 0: start from new; 1: provide the keep
	// store the parameters


	auto start = std::chrono::system_clock::now();


	// read in 
	//1) annual manta-tow data
	//2) AMPTO data 
	//-------------------Begin read in data--------------------------------------

	///////////////////////////////////////////
	/*	aFile.open("predlengthchro.txt", ios::in);
		no_chl_rec = count_Data(aFile);
		cout << "Nrow Chlorophyll Length data " << no_chl_rec << endl;
		aFile.clear();
		aFile.seekg(0);
		struct chlo_length *chlo_length;
		chlo_length = new struct chlo_length[no_chl_rec];
		get_chlo_legth_dat(aFile, &(*chlo_length));
		aFile.close();
		*/

	// for location where no COTS were seen on the first two years. 
	/*
	cFile.open("ini_size_stru_fisk_zero.txt", ios::in);
	bins_zero = count_Data(cFile);
	cout << "Number of Length bin: " << bins_zero << endl;
	cFile.clear();
	cFile.seekg(0);

	struct ini_length * ini_leng_zero;
	ini_leng_zero = new struct ini_length[bins_zero];
	get_ini_leng_dat(cFile, ini_leng_zero);
	cFile.close();
	*/

	//////////////////////////////////////////////
	bFile.open("initial_length_v2.txt", ios::in);
	bins = count_Data(bFile);
	cout << "Number of Length bin: " << bins << endl;
	bFile.clear();
	bFile.seekg(0);

	struct ini_length * ini_leng;
	ini_leng = new struct ini_length[bins];
	get_ini_leng_dat(bFile, ini_leng);
	bFile.close();

	// Read in Locations and relative abundance
	aFile.open("starting_wt.txt", ios::in);
	no_loc = count_Data(aFile);
	aFile.clear();
	aFile.seekg(0);
	struct reef_month *temp_mth_reef;
	temp_mth_reef = new struct reef_month[no_loc];
	get_ini_dist(aFile, temp_mth_reef);
	aFile.close();

	aFile.open("MonthlyCoralAdjCycle.txt", ios::in);
	coral_cov_rec = count_Data(aFile);
	aFile.clear();
	aFile.seekg(0);
	struct coral_cov *CC;
	CC = new struct coral_cov[coral_cov_rec];
	get_coral_coverage(aFile, CC, START_YR, START_MTH);
	aFile.close();

	// assign CC to mth_reef
	convert_coralcover_foo(CC, temp_mth_reef, no_loc, SIM_LENGTH, coral_cov_rec);
	delete[] CC;
	CC = nullptr;

	bFile.open("chlorophylldata_all.txt", ios::in);
	chlo_rec = count_Data(bFile);
	cout << "Chlorophyll data " << chlo_rec << endl;
	bFile.clear();
	bFile.seekg(0);
	struct chloro_ori *chlo_temp;
	chlo_temp = new struct chloro_ori[chlo_rec];
	get_chlorophyll_dat(bFile, chlo_temp, START_YR, START_MTH);
	bFile.close();

	struct chlorophyll *chlo_dat;
	chlo_dat = new struct chlorophyll[no_loc];

	for (i = 0; i < no_loc; i++){
		chlo_dat[i].REEF_ID = temp_mth_reef[i].REEF_ID;
	}

	format_chloro(chlo_temp, chlo_dat, no_loc, SIM_LENGTH, chlo_rec);

	delete[] chlo_temp;
	chlo_temp = nullptr;

	aFile.open("Connect_09.txt", ios::in);
	n_conn = count_Data(aFile);
	aFile.clear();
	aFile.seekg(0);
	struct connect_ori *con_temp = new struct connect_ori[n_conn];
	get_connect_dat(aFile, con_temp);
	aFile.close();

	// convert struct of connect_ori to 
	struct con *connect = new struct con[no_loc];

	int ctr1, yr;
	int s_yr;

	// assign soucing reef
	for (i = 0; i < no_loc; i++){
		s_yr = START_YR + 1; // recruitment starts in the following year
		connect[i].REEF_ID = temp_mth_reef[i].REEF_ID;
		ctr1 = 0;
		for (j = 0; j < (no_loc - 1)*(SIM_YEARS + 1); j++){
			connect[i].year[j] = s_yr;
			if (j == ((no_loc - 1) * (ctr1 + 1) - 1)){
				s_yr++;
				ctr1++;
			}

		}
		// assign sinking reefs
		l = 0;
		for (yr = 0; yr < SIM_YEARS + 1; yr++){
			for (k = 0; k < no_loc; k++){
				if (connect[i].REEF_ID != temp_mth_reef[k].REEF_ID){
					connect[i].to[l] = temp_mth_reef[k].REEF_ID;
					l++;
				}
			}
		}
	}
	fomat_connect(con_temp, connect, no_loc, n_conn, SIM_YEARS);

	delete[] con_temp;
	con_temp = nullptr;

	// read in reef size

	int nr;
	/*aFile.open("areaofReef.txt", ios::in);
	nr = count_Data(aFile);
	struct reef_size *size = new struct reef_size[nr];
	aFile.clear();
	aFile.seekg(0);
	get_reef_size_dat(aFile, size);
	aFile.close();

	for (i = 0; i < nr; i++){
	for (j = 0; j < no_loc; j++){
	if (temp_mth_reef[j].REEF_ID == size[i].reef_ID){
	temp_mth_reef[j].size = size[i].size;
	}
	}
	}

	delete[] size;
	size = nullptr;
	*/
	///---- read water temp data------
	bFile.open("Complete Pooled water temp.txt", ios::in);
	nr = count_Data(bFile);
	struct w_temp *wat_temp = new struct w_temp[nr];
	bFile.clear();
	bFile.seekg(0);
	get_wtemp_dat(bFile, wat_temp, START_YR, START_MTH);
	bFile.close();

	///// ---------------------- read the catch data

	int no_obs_rec;

	/*change to distribution on 08/09/2015
	aFile.open("LTMP_abun.txt", ios::in);
	no_obs_rec = count_Data(aFile);
	aFile.clear();
	aFile.seekg(0);
	struct obs_abund *obs_abundance = new struct obs_abund[no_obs_rec];
	get_obs_dat(aFile, obs_abundance);
	aFile.close();
	*/
	aFile.open("COT_abund_dis_cpue.txt", ios::in);
	no_obs_rec = count_Data(aFile);
	aFile.clear();
	aFile.seekg(0);
	struct obs_abund *obs_abundance = new struct obs_abund[no_obs_rec];
	get_obs_dat(aFile, obs_abundance);
	aFile.close();
	//-------------------------------end read in data--------------------------------

	///------ remove data

	// remove files
	remove("res_out.txt");
	remove("pop_str.txt");
	remove("rho_out.txt");
	remove("walltime.txt");
	remove("params.txt");
	remove("acc_rat.txt");
	//remove("annual.txt");


	/// assign universal initial structure
	// out put sort
	myFile.open("res_out.txt", ios::app);
	// just checking if struct2mat works
	myFile << "particle" << "\t" << "Loc" << "\t" << "Year" << "\t" << "Mth" << "\t" << "n_cots" << "\t" << "no_female" << "\t" << "no_male" << "\t" << "L15" << "\t" << "L1525" << "\t" << "L2540" << "\t" << "L40p" << endl;
	myFile.close();


	// Of a fixed number of COTS the follow functiond distribute COTS number according to the relative weight
	ini_pop_distribution(temp_mth_reef, no_loc, INI_POP_SIZE);

	myFile.open("params.txt", ios::app);
	//myFile << "ITI" << "\t" << "Par" << "\t" << "Summary_stat" << "\t" << "Beta" << "\t" << "detect" << "\t" << "jv_sur" << "\t" << "k" << "\t" << "larv_mor" << "\t" << "x0" << endl;
	myFile << "Par" << "\t" << "Summary_stat" << "\t";
	for (i = 0; i < no_loc; i++){
		myFile << "alpha_" << i << "\t";
	}
	myFile << "Beta" << "\t" << "Beta_sigma" << "\t";
	for (i = 0; i < no_loc; i++){
		myFile << "detect_" << i << "\t";
	}
	myFile << "Error" << "\t" << "jv_sur" << "\t" << "k" << "\t" << "larv_mor" << "\t" << "x0" << "\t" << endl;
	myFile.close();

	myFile.open("acc_rat.txt", ios::app);
	myFile << "ITI ratio i_acc R p_acc RN" << endl;
	myFile.close();

	///////////OMP////////////////////
	omp_set_dynamic(2);
	unsigned short max_thread = omp_get_max_threads();
	base_generator_type eng[max_thread];
	omp_set_num_threads(max_thread);



	struct parameters *params = new struct parameters[PARTICLE_SIZE];

	// add two elements to parameters

	// first runs
	Context dist;
	vector<double> rho(PARTICLE_SIZE);

	double e_max; // first epi_max
	int N_drop = (int)round(PARTICLE_SIZE* RATE); // numebr of partcles to be dropped
	struct parameters *param_keepers = new struct parameters[PARTICLE_SIZE - N_drop];

	if (OPTION == 0){

#pragma omp parallel for schedule(dynamic, 3) private(i, j, p) shared(rho,params, START_YR, no_loc, obs_abundance, PARTICLE_SIZE, bins, no_obs_rec,  temp_mth_reef,ini_leng, SIM_LENGTH, INI_SEX_RAT, eng) 
		for (p = 0; p < PARTICLE_SIZE; p++){

			double sum_stats = 100000;
			int flg = 1; // flag ==1 means population extinct
			int seed = boost::uuids::detail::seed_rng()();
			struct mth_out *out = new struct mth_out[SIM_LENGTH*no_loc]; // store the monthly records


			eng[omp_get_thread_num()].seed(seed);
			while ((sum_stats >= START_S) || (flg == 1)){
				//eng[omp_get_thread_num()].seed(static_cast<unsigned int>((seed++)*p + time(0))*max_thread*getpid());
				int inini;
				do{
					params[p].beta = dist.uniform_real_dis(BETA_MIN, BETA_MAX, eng[omp_get_thread_num()]);
					params[p].beta_sigma = dist.lognormal_dis(BETA_SIG_MEAN, BETA_SIG_STD, eng[omp_get_thread_num()]);
					params[p].jv_sur = dist.lognormal_dis(JUV_SUR_MEAN, JUV_SUR_STD, eng[omp_get_thread_num()]);
					params[p].k = dist.lognormal_dis(ADULT_MOR_K_MEAN, ADULT_MOR_K_STD, eng[omp_get_thread_num()]);
					params[p].larv_mor = 1 - dist.lognormal_dis(LAR_DAILY_MOR_MEAN, LAR_DAILY_MOR_STD, eng[omp_get_thread_num()]);
					params[p].x0 = dist.lognormal_dis(ADULT_MOR_X0_MEAN, ADULT_MOR_X0_STD, eng[omp_get_thread_num()]);
					params[p].error_N = dist.uniform_real_dis(VAR_ERROR_MIN, VAR_ERROR_MAX, eng[omp_get_thread_num()]);

					for (j = 0; j < no_loc; j++){
						params[p].detect[j] = dist.uniform_real_dis(DETECT_LO, DETECT_HI, eng[omp_get_thread_num()]);
						params[p].alpha[j] = dist.lognormal_dis(DIS_MEAN, DIS_STD, eng[omp_get_thread_num()]);
					}
					inini= check_ini_para(params, p, DIM, no_loc);

				} while (inini == 0);
 
					//----------------------------------call the process model here!!
					//setting up starting reef matrix
				struct reef_month *mth_reef = new struct reef_month[no_loc];

				for (i = 0; i < no_loc; i++){
					mth_reef[i] = temp_mth_reef[i];
					mth_reef[i].yr = START_YR;
					mth_reef[i].mth = START_MTH;
				}

				ProcessModel(&params[p], mth_reef, ini_leng, obs_abundance, bins, no_loc, SIM_LENGTH, INI_SEX_RAT, max_thread, no_obs_rec, sum_stats, flg, eng[omp_get_thread_num()], START_S, START_S, TOLERANCE_S, out);

				delete[]mth_reef;
				mth_reef = nullptr;

			}
			rho[p] = sum_stats;
			cout << p << " " << rho[p] << " flag " << flg << endl;

		}// end particles of the first iteration

		e_max = *max_element(rho.begin(), rho.end()); // first epi_max

		cout << "First e max " << e_max << endl;

		/*myFile.open("params.txt", ios::app);

		for (p = 0; p < PARTICLE_SIZE; p++){
		myFile << p << "\t" << rho[p] << "\t" << params[p].beta << "\t" << params[p].beta_sigma << "\t" << params[p].jv_sur << "\t" << params[p].k << "\t" << params[p].larv_mor << "\t" << params[p].x0 << "\t" << params[p].error_N << "\t";
		for (i = 0; i < no_loc; i++){
		myFile << params[p].alpha[i] << "\t" << params[p].detect[i] << "\t";
		}
		myFile << endl;
		}
		myFile.close();
		*/
	}

	else{ // read in the params[p] data

		aFile.open("params_old.txt", ios::in);
		get_old_params(aFile, params, rho);
		aFile.close();
		e_max = *max_element(rho.begin(), rho.end()); // first epi_max

		cout << "First e max " << e_max << endl;
	}

	double e_next;
	long long R, RN;
	R = INITIAL_R;
	long double p_acc;
	int tck;
	int u = 0;
	
	do{// et to e_max
		int ctr = 0;
		// ----- adaptive design
		RN = R* N_drop;
		int i_acc = 0;
		tck = 0;
		// --- adaptive
		// assign the rank to params
		for (auto i : sort_indexes(rho)) {
			params[i].rank = ctr;
			//cout << i << " " << rho[i] << " " << params[i].rank << endl; // this is just to show values from largest to smallest; and the location of value (i) in the vector  
			ctr++;
		}
		//for checking if the code works
		/*for (p = 0; p < PARTICLE_SIZE; p++){
		cout << rho[p] << " " << params[p].rank << endl;
		}*/
		// next step-replenish
		//----------------- keeppers

		int ctr_cc = 0;


		for (p = 0; p < PARTICLE_SIZE; p++){
			if (params[p].rank >= N_drop){
				param_keepers[ctr_cc] = params[p];
				ctr_cc++;
			}
			if (params[p].rank == N_drop){
				e_next = rho[p];
			}
		}

		cout << "next " << e_next << endl;
		/*
		myFile.open("params.txt", ios::app);

		for (i = 0; i < ctr_cc; i++){
		myFile << i << "\t" << param_keepers[i].rank << "\t" << param_keepers[i].beta << "\t" << param_keepers[i].beta_sigma << "\t" << param_keepers[i].jv_sur << "\t" << param_keepers[i].k << "\t" << param_keepers[i].larv_mor << "\t" << param_keepers[i].x0 << "\t" << param_keepers[i].error_N << "\t";
		for (j = 0; j < no_loc; j++){
		myFile << param_keepers[i].alpha[j] << "\t" << param_keepers[i].detect[j] << "\t";
		}
		myFile << endl;
		}
		myFile.close();
		*/
		// replensihing

		// parallel this part
		DIM = 13;

		// variance-covariance of

		MatrixXd S(DIM, DIM);
		MatrixXd temp_mat(PARTICLE_SIZE - N_drop, DIM);
		struct2mat(param_keepers, temp_mat, PARTICLE_SIZE - N_drop);

		MatrixXd x_bar = temp_mat.rowwise() - temp_mat.colwise().mean();
		S = 2 * (x_bar.adjoint() * x_bar) / double(temp_mat.rows());
		// end of mean and variance-cov matrix

		// code is fine to here
#pragma omp parallel for schedule(dynamic, 1) private(i,p) firstprivate(ctr_cc)shared(START_YR, no_loc, chlo_dat, obs_abundance, PARTICLE_SIZE, bins, no_obs_rec, temp_mth_reef,ini_leng, SIM_LENGTH, INI_SEX_RAT, DIM, S, N_drop, param_keepers, params, e_max, e_next, tck, RN,i_acc, rho, eng)
		for (p = 0; p < PARTICLE_SIZE; p++){

			double sum_stats;
			int flg = 0; // flg ==1 means population extinct
			struct parameters *temp_params = new struct parameters[1];
			struct mth_out *out = new struct mth_out[SIM_LENGTH*no_loc]; // store the monthly records

			//cout << "bf " << rho[p] << endl;
			//int seed = boost::uuids::detail::seed_rng()();

			//eng[omp_get_thread_num()].seed(seed);

			if (params[p].rank <= N_drop){// && tck< RN){ // replenish N_drop
				int flag = 0;
				int k = 0; 

				do{// line 2.5 to 2.12
#pragma omp atomic
					k++; 
#pragma omp atomic
					tck++;
					// step 1: select from keepers
					int oo;
					// q(.|theta)- line 2.4 
					// mean vector
					// draws - equal probability
					VectorXd mean(DIM);

					oo = dist.uniform_int_dis(0, ctr_cc - 1, eng[omp_get_thread_num()]);// ctr-1 is total number of particles kept
					//#pragma omp critical
					struct2array(param_keepers, mean, oo);

					// draw theta **

					EigenMultivariateNormal<double> normX(mean, S, boost::uuids::detail::seed_rng()());
					int flgg = 0;
					MatrixXd temp(DIM, 1);
					while (flgg == 0){
						temp = normX.samples(1);
						flgg = check_array(temp, DIM);
					}
					/*
#pragma omp critical
					{
					cout << "p " << p << " prop " << temp << endl;
					cout << "oo " << oo << "mean" << mean << endl;
					}
					*/

					//cout << "mean" <<  mean << "\t" << "temp" << temp << endl;
					// work out the MH ratio
					double MH;
					// order is extremely important 


					double top = lognormal_pdf_ln(DIS_MEAN, DIS_STD, temp(0)) + lognormal_pdf_ln(DIS_MEAN, DIS_STD, temp(1)) + lognormal_pdf_ln(DIS_MEAN, DIS_STD, temp(2)) + unif_pdf_ln(BETA_MIN, BETA_MAX, temp(3)) + lognormal_pdf_ln(BETA_SIG_MEAN, BETA_SIG_STD, temp(4)) - 3 * log(DETECT_HI - DETECT_LO) + unif_pdf_ln(VAR_ERROR_MIN, VAR_ERROR_MAX, temp(8)) + lognormal_pdf_ln(JUV_SUR_MEAN, JUV_SUR_STD, temp(9)) + lognormal_pdf_ln(ADULT_MOR_K_MEAN, ADULT_MOR_K_STD, temp(10)) + lognormal_pdf_ln(LAR_DAILY_MOR_MEAN, LAR_DAILY_MOR_STD, 1 - temp(11)) + lognormal_pdf_ln(ADULT_MOR_X0_MEAN, ADULT_MOR_X0_STD, temp(12)) + multivariate_normal_ln_pdf(temp, S, mean);

					double bot = lognormal_pdf_ln(DIS_MEAN, DIS_STD, mean(0)) + lognormal_pdf_ln(DIS_MEAN, DIS_STD, mean(1)) + lognormal_pdf_ln(DIS_MEAN, DIS_STD, mean(2)) + unif_pdf_ln(BETA_MIN, BETA_MAX, mean(3)) + lognormal_pdf_ln(BETA_SIG_MEAN, BETA_SIG_STD, mean(4)) - 3 * log(DETECT_HI - DETECT_LO) + unif_pdf_ln(VAR_ERROR_MIN, VAR_ERROR_MAX, mean(8)) + lognormal_pdf_ln(JUV_SUR_MEAN, JUV_SUR_STD, mean(9)) + lognormal_pdf_ln(ADULT_MOR_K_MEAN, ADULT_MOR_K_STD, mean(10)) + lognormal_pdf_ln(LAR_DAILY_MOR_MEAN, LAR_DAILY_MOR_STD, 1 - mean(11)) + lognormal_pdf_ln(ADULT_MOR_X0_MEAN, ADULT_MOR_X0_STD, mean(12)) + multivariate_normal_ln_pdf(mean, S, temp);

					MH = exp(top - bot);
					/*
#pragma omp critical
					{
					cout << " t " << top << " b " << bot;
					cout << " p " << p << " MH " << MH << endl;
					}
					*/

					if (MH > 1){ MH = 1; }// else MH is MH

					if (MH > dist.uniform_real_dis(0, 1, eng[omp_get_thread_num()])){
						array2struct(temp_params, temp, 0);
						// simulate data
						struct reef_month *mth_reef = new struct reef_month[no_loc];

						for (i = 0; i < no_loc; i++){
							mth_reef[i] = temp_mth_reef[i];
							mth_reef[i].yr = START_YR;
							mth_reef[i].mth = START_MTH;
						}
						/*
#pragma omp critical
						{
						cout << params[p].beta << " " << params[p].beta_sigma << " " << params[p].error_N << " " << params[p].jv_sur << " " << params[p].k << " " << params[p].larv_mor << " " << params[p].rank << " ";

						for (j = 0; j < no_loc; j++){
						cout << params[p].detect[j] << " " << params[p].alpha[j] << " ";
						}
						cout << endl;
						}
						*/
						ProcessModel(&temp_params[0], mth_reef, ini_leng, obs_abundance, bins, no_loc, SIM_LENGTH, INI_SEX_RAT, max_thread, no_obs_rec, sum_stats, flg, eng[omp_get_thread_num()], e_max, e_next, TOLERANCE_S, out);

						// if summary stat is less than E_Max
						if (sum_stats < e_next && flg == 0){
#pragma omp critical
							{
								cout << " p " << p << " sum_stat " << sum_stats << " flag " << flg << endl;
							}
							if (e_next < TOLERANCE_S){
#pragma omp critical
								{
									myFile.open("res_out.txt", ios::app);
									for (j = 0; j < no_loc*SIM_LENGTH; j++){
										myFile << p << "\t" << out[j].reef_ID << "\t" << out[j].year << "\t" << out[j].mth << "\t" << out[j].pop_size << "\t" << out[j].no_female << "\t" << out[j].no_male << "\t" << out[j].l15 << "\t" << out[j].l1525 << "\t" << out[j].l2540 << "\t" << out[j].l40p << endl;

									}
									myFile.close();

								}

							}
							flag = 1;
							array2struct(params, temp, p);
							rho[p] = sum_stats;
#pragma omp atomic
							i_acc++;
						}
						else{
							flag = 0;
						}
						delete[]mth_reef;
						mth_reef = nullptr;
					}
					else{ // MH else
						flag = 0;
					}
				}while(flag == 0 && k < R);

			} // end if: replenish N_drop

			delete[] temp_params;
			temp_params = nullptr;
			delete[] out;
			out = nullptr;
		} // end of p iteration- replenishing 
		cout << " icc " << i_acc << endl;

		p_acc = (double)i_acc / RN;
		R = log(0.1) / log(1 - p_acc);

		// update e_max;
		cout << "acceptance ratio " << p_acc << endl;
		cout << " R " << R << endl;
		cout << endl;
		u++;
		double uu = (double)i_acc / tck;
		myFile.open("acc_rat.txt", ios::app);
		myFile << u << " " << uu << " " << i_acc << " " << R << " " << p_acc << " "<< RN << endl;
		myFile.close();

		e_max = *max_element(rho.begin(), rho.end());
		cout << "e_max" << e_max << endl;

	} while (e_max > TOLERANCE_S);// && tck < RN);


	auto end = std::chrono::system_clock::now();
	auto elapsed = end - start;
	myFile.open("walltime.txt", ios::app);
	myFile << elapsed.count() << '\n';
	myFile.close();


	////-----------------
	myFile.open("params.txt", ios::app);

	for (p = 0; p < PARTICLE_SIZE; p++){
		myFile << p << "\t" << rho[p] << "\t";

		for (i = 0; i < no_loc; i++){
			myFile << params[p].alpha[i] << "\t";
		}

		myFile << params[p].beta << "\t" << params[p].beta_sigma << "\t";
		for (i = 0; i < no_loc; i++){
			myFile << params[p].detect[i] << "\t";
		}
		myFile << params[p].error_N << "\t" << params[p].jv_sur << "\t" << params[p].k << "\t" << params[p].larv_mor << "\t" << params[p].x0 << "\t" << endl;
	}
	myFile.close();

	cout << "Congratulations XXOO!!" << endl;

	delete[] ini_leng;
	delete[] temp_mth_reef;
	delete[] params;
	delete[] wat_temp;
	delete[] connect;
	delete[] chlo_dat;
	delete[] obs_abundance;
	delete[] param_keepers;

	param_keepers = nullptr;
	obs_abundance = nullptr;
	chlo_dat = nullptr;

	connect = nullptr;
	wat_temp = nullptr;
	params = nullptr;

	ini_leng = nullptr;

	return 0;
}
