// starting month is fixed to dec

#include<stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
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

#include <boost/random/mersenne_twister.hpp>

using namespace std;

typedef boost::lagged_fibonacci607  RandomGeneratorType;

void ProcessModel(struct parameters *params, struct reef_month *mth_reef, struct ini_length *ini_leng, struct obs_abund *obs_abundance, int bins, int no_loc, int SIM_LENGTH, double INI_SEX_RAT, int max_thread, int no_obs_rec, double &summary_stat, int &flag, base_generator_type& eng, int starting_S, double next_s, double min_TOL, struct mth_out *&out); // flag ==1 means population extinct