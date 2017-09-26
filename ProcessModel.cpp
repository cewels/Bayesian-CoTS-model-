#include "ProcessHeader.h"


void ProcessModel(struct parameters *params, struct reef_month *mth_reef, struct ini_length *ini_leng, struct obs_abund *obs_abundance, int bins, int no_loc, int SIM_LENGTH, double INI_SEX_RAT, int max_thread, int no_obs_rec, double &summary_stat, int &flag, base_generator_type& eng, int starting_S, double next_s, double min_TOL, struct mth_out *&out){
	//in: data
	// params[p]- one set only
	// mth_reef
	// no_loc: number of location
	// SIM_LENGTH
	// struct ini_leng: struct store initial population length structure
	// bin: number of length classes
	// INI_SEX_RAT: initial sex ratio
	// obs_abundance: observation data
	// no_obs_rec: number of observation records
	// out: summary statistics

	Context dist;
	int maxOfpop;
	int i, j;

	//base_generator_type eng;
	//eng.seed(static_cast<unsigned int>(++seed + time(0))*max_thread*getpid());

	LinkedList_COTS *COTS;
	COTS = new LinkedList_COTS[no_loc];
	struct reef_annual *annual = new struct reef_annual[SIM_LENGTH*no_loc];

	int counter = 0;

	// Set up initial COTs population

	for (i = 0; i < no_loc; i++){
		if (mth_reef[i].N0 >0){
			int *length = new int[mth_reef[i].N0];
			/* if the reef does not have COTS observed in the first two years
			if (mth_reef[i].ini_obs == 0){
			// initial length structure for reefs with initial observations ==0
			//ini_length_struct(length, mth_reef[i].N0, ini_leng_zero, bins_zero);
			}*/

			if (mth_reef[i].ini_obs == 1){
				ini_length_struct(length, mth_reef[i].N0, ini_leng, bins, mth_reef[i].N0_det, eng);
			}

			for (j = 0; j < mth_reef[i].N0; j++){
				COTS[i].insert_initial_pop(length[j], dist.bernoullie_dis(INI_SEX_RAT, eng), eng);
			}
			delete[] length;
			length = nullptr;
		}

	}

	// assign year and month to the initial reef structure
	//--------------------------------
	// Starting from time =0, looping through all years
	int x = 0;
	int tck = 0;
	int x1 = 1;
	int tck1 = 0;
	int nov = 11;

	int ctr = 0;
	double food;

	int tm;
	int tktk = 0;

	// cheating 
	long long *JAN = new long long[no_loc];
	long long *FEB = new long long[no_loc];
	long long *MAR = new long long[no_loc];

	for (tm = 0; tm < SIM_LENGTH; tm++){
		//#pragma omp critical
		//				cout << omp_get_thread_num() << " " << t << " " << p << " " << tm << endl; 
		RowVectorXi pop_check(no_loc);
		RowVectorXi pop_m(no_loc);
		RowVectorXi pop_f(no_loc);

		//cout << "S2" << endl;
		//-----------------------------------------------------------------------------------------------------
		//-----------------------------------/ JAN, FEB and MARCH only------------------------------------------------------------------

		if (tm == x1){
			//Natural mortality
			// Catharina's model, using the following. if I use the Sven's model, notural mortality is fused with chlorophyll
			//mth_reef[i].eggs_out = nat_mor(&mth_reef[i].eggs_out, sur_larvae);
			for (i = 0; i < no_loc; i++){
				/*c_rec = 0;
				while (chlo_dat[i].s_day[c_rec] != tm){
				//#pragma omp atomic
				c_rec++;  // to find the record corresponding to tm in the chlorophyll data
				}
				// chlorophyll - Fabricious
				//mth_reef[i].eggs_out = sur_chlo(&mth_reef[i].eggs_out, &chlo_dat[i].chlo[c_rec]);
				// Sven's model
				t_ctr = 0;
				while (wat_temp[t_ctr].sim != tm){
				//#pragma omp atomic
				t_ctr++;  // to find the record corresponding to tm in the chlorophyll data
				}*/

				if (mth_reef[i].eggs_out > 0){
					//14/07
					//mth_reef[i].eggs_out = sur_chlo_sven(&mth_reef[i].eggs_out, &chlo_dat[i].chlo[c_rec], wat_temp[t_ctr].wTemp, params[p].larv_mor, eng);
					mth_reef[i].eggs_out = nat_mor(&mth_reef[i].eggs_out, (1 - params->larv_mor), eng);

				}
				else{ mth_reef[i].eggs_out = 0; }

				// reset larvae in 
				mth_reef[i].larv_in = 0;
			}
			/* 14/07n_conn = (no_loc - 1)*(SIM_YEARS + 1);
			// think more about this
			recru_con(mth_reef, no_loc, connect, n_conn, eng);
			*/
			//simple_rec(mth_reef, no_loc, ]); // only for two locations
			//rec_2para(mth_reef, params->gamma, params->alpha, no_loc, eng);
			rec_1para(mth_reef, params->alpha, no_loc, eng);

			for (i = 0; i < no_loc; i++){
				if (tm == tck1 * 12 + 1){// jan
					JAN[i] = mth_reef[i].larv_in;
				}
				if (tm == tck1 * 12 + 2){
					FEB[i] = mth_reef[i].larv_in;
				}
				if (tm == tck1 * 12 + 3){
					MAR[i] = mth_reef[i].larv_in;
				}
			}

			if (x1 - tck1 * 12 == 3){ //MARCH
				//#pragma omp atomic
				tck1++;
				x1 = tck1 * 12 + 1;
				for (i = 0; i < no_loc; i++){
					COTS[i].reset_spawn(); // reset spawning
				}
			}
			else{
				//#pragma omp atomic
				x1++;
			}
		}

		///// JULY, AUGUST and SEP-------------------------------------------------


		if (tm == tktk * 12 + 7){
			// assume 10% survival per month
			for (i = 0; i < no_loc; i++){

				if (JAN[i] != 0){
					JAN[i] = dist.binomial_dis(JAN[i], pow(params->jv_sur, 4), eng);
					if (JAN[i] != 0){
						if (COTS[i].start != NULL){
							for (j = 0; j < JAN[i]; j++){
								COTS[i].updatepop(dist.uniform_int_dis(9, 20, eng), dist.bernoullie_dis(INI_SEX_RAT, eng));
							}
						}
						else{ // previous population die out
							for (j = 0; j < JAN[i]; j++){
								COTS[i].insert_new_pop(dist.uniform_int_dis(9, 20, eng), dist.bernoullie_dis(INI_SEX_RAT, eng));
							}
						}
					}
					JAN[i] = 0;
				}
			}
		}
		if (tm == tktk * 12 + 8){
			// assume 10% survival per month
			for (i = 0; i < no_loc; i++){
				if (FEB[i] != 0){
					FEB[i] = dist.binomial_dis(FEB[i], pow(params->jv_sur, 4), eng);
					if (FEB[i] != 0){
						if (COTS[i].start != NULL){
							for (j = 0; j < FEB[i]; j++){
								COTS[i].updatepop(dist.uniform_int_dis(9, 20, eng), dist.bernoullie_dis(INI_SEX_RAT, eng));
							}
						}
						else{ // previous population die out
							for (j = 0; j < FEB[i]; j++){
								COTS[i].insert_new_pop(dist.uniform_int_dis(9, 20, eng), dist.bernoullie_dis(INI_SEX_RAT, eng));
							}
						}
					}
					FEB[i] = 0;
				}
			}
		}
		if (tm == tktk * 12 + 9){
			// assume 10% survival per month
			for (i = 0; i < no_loc; i++){
				if (MAR[i] != 0){
					MAR[i] = dist.binomial_dis(MAR[i], pow(params->jv_sur, 4), eng);
					if (MAR[i] != 0){
						if (COTS[i].start != NULL){
							for (j = 0; j < MAR[i]; j++){
								COTS[i].updatepop(dist.uniform_int_dis(9, 20, eng), dist.bernoullie_dis(INI_SEX_RAT, eng));
							}
						}
						else{ // previous population die out
							for (j = 0; j < MAR[i]; j++){
								COTS[i].insert_new_pop(dist.uniform_int_dis(9, 20, eng), dist.bernoullie_dis(INI_SEX_RAT, eng));
							}
						}
						MAR[i] = 0;
					}
				}
			}
			tktk++;
		}


		//-------------------------------Every Month ----------------------------------------------------------------------

		for (i = 0; i < no_loc; i++){
			//cout << "ITI " << t << " par " << p << " step " << tm << " " << "Location " << mth_reef[i].REEF_ID << "\t" << "adult: " << COTS[i].count_COTS() << endl;
			//cout << "S3" << endl;
			int N_c = COTS[i].count_COTS();


			if (N_c != 0){
				// growth adult - older version
				// crude inclusion of coral
				if (tm == 0){
					food = 1;
				}
				else{
					food = (double)(mth_reef[i].cov[tm + 1] - mth_reef[i].cov[tm]) / (double)(mth_reef[i].cov[0]);
				}
				/* for simpler model
				food = 1;
				COTS[i].growth(food, params[p].beta, params[p].beta_sigma, eng);
				*/
				// step function
				//COTS[i].growth_feed_step(food, params[p].beta, params[p].beta_sigma, eng);

				double alpha = 0.907 + (food*params->beta / 0.5);

				COTS[i].growth_feed_function(&alpha, params->beta_sigma, eng);

				COTS[i].remove_starved_cots();

				// update age
				COTS[i].update_age();
				// mortality
				COTS[i].nat_mor_A(MAX_SUR, params->k, params->x0, params->jv_sur, eng);
			}
			//cout << " step " << tm << " " << "Location " << mth_reef[i].REEF_ID <<  "\t" << "adult: " << COTS[i].count_size10() << endl;

			//mth_reef[i].N6p[tm] = COTS[i].count_den_food();

			pop_check(i) = COTS[i].count_COTS();
			pop_m(i) = COTS[i].count_COTS_m();
			pop_f(i) = COTS[i].count_COTS_f();

		}

		int o;
		maxOfpop = pop_check.maxCoeff(&o);

		if (maxOfpop <= 2){
			if (pop_m(o) == 0 || pop_f(o) == 0){
				//cout << "Population is extinct" << endl;
				break;
			}
		}

		/////////--------------------end all months
		///*********************************Dec, Jan and Feb only
		int sex_m_m;
		if (tm == x){
			for (i = 0; i < no_loc; i++){
				// reset eggs
				sex_m_m = COTS[i].count_COTS_sex_mat_m();
				// only there is at least one male and one female on the reef, there will be spawning
				if (sex_m_m>0){ 
					mth_reef[i].eggs_out = COTS[i].spawn(eng);
				}
			}
			if (x - tck * 12 == 2){ //Feb
				//#pragma omp atomic
				tck++;
				x = tck * 12;
			}
			else{
				//#pragma omp atomic
				x++;
			}
		}
		for (i = 0; i < no_loc; i++){
			annual[ctr].year = mth_reef[i].yr;
			annual[ctr].mth = mth_reef[i].mth;
			annual[ctr].REEF_ID = mth_reef[i].REEF_ID;
			double sigg = abs(params->error_N);
			annual[ctr].pred_pop_size = COTS[i].count_24p_COTS(params->detect[i], sigg, mth_reef[i].N0_det, eng);
			if (next_s < min_TOL){
				out[counter].year = annual[ctr].year;
				out[counter].mth = annual[ctr].mth;
				out[counter].reef_ID = annual[ctr].REEF_ID;
				out[counter].pop_size = annual[ctr].pred_pop_size;
				out[counter].no_female = COTS[i].count_COTS_f();
				out[counter].no_male = COTS[i].count_COTS_m();
				out[counter].l15 = COTS[i].count_size15();
				out[counter].l1525 = COTS[i].count_size1525();
				out[counter].l2540 = COTS[i].count_size2540();
				out[counter].l40p = COTS[i].count_size40p();
				//#pragma omp atomic
				counter++;
			}
			ctr++;
		}

		for (i = 0; i < no_loc; i++){
			//#pragma omp critical
			mth_reef[i].mth++;
			// update year and month
			if (mth_reef[i].mth == 13){
				//#pragma omp atomic
				mth_reef[i].yr++;
				mth_reef[i].mth = 1;
			}

		}


	}// end of tm- COTs times

	if (maxOfpop <= 2){
		flag = 1;
		summary_stat = starting_S;
	}

	// compute summary statsitics 
	else{ // more than two cots remaining in the system 
		// step 1: merge annual and obs_abundance
		flag = 0;
		for (i = 0; i < ctr; i++){ /// ctr = number of simulation year* number of location
			for (j = 0; j < no_obs_rec; j++){
				if (annual[i].year == obs_abundance[j].year && annual[i].REEF_ID == obs_abundance[j].reef_ID && annual[i].mth == obs_abundance[j].mth && obs_abundance[j].CPUE >= 0){
					annual[i].observed = exp(obs_abundance[j].CPUE);
				}
			}
		}

		// step 2: calculating summary statistics

		summary_stat = summary_stat_sd(mth_reef, no_loc, annual, ctr);

	}

	/// The following code is only applicable during the forward simulation phase
	/*
	if (SIMM==1){
	myFile.open("res_out.txt", ios::app);
	for (j = 0; j < no_loc*SIM_LENGTH; j++){
	myFile << p << "\t" << out[j].reef_ID << "\t" << out[j].year << "\t" << out[j].mth << "\t" << out[j].pop_size << "\t" << out[j].no_female << "\t" << out[j].no_male << "\t" << out[j].l15 << "\t" << out[j].l1525 << "\t" << out[j].l2540 << "\t" << out[j].l40p << endl;

	}
	myFile.close();

	}
	*/
	///////-----------delete linked lists 
#pragma omp parallel for
	for (i = 0; i < no_loc; i++){
		COTS[i].~LinkedList_COTS();
	}
	delete[] annual;
	delete[] JAN;
	delete[] FEB;
	delete[] MAR;
	delete[] COTS;
	COTS = nullptr;
	JAN = nullptr;
	FEB = nullptr;
	MAR = nullptr;
	annual = nullptr;
	return; 
}
