#include "ANIMALS.h"
#include "Misc.h"
#include <numeric>


using namespace std;


LinkedList_COTS::LinkedList_COTS(){
	start = NULL;
	q = NULL;
}

LinkedList_COTS::~LinkedList_COTS()
{
	makeEmpty_COTS();
}

bool LinkedList_COTS::isEmpty_COTS()  // O(1)
{
	return start == 0;
}

void LinkedList_COTS::makeEmpty_COTS() // O(n)
{
	while (start){
		q = start->next;
		delete start;
		start = q;
	}
}
void LinkedList_COTS::print_COTS(){
	for (q = start; q != NULL; q = q->next)
	{
		cout << q->size << endl;
	}
	return;
}
int LinkedList_COTS::count_COTS()
{
	// counting every single COTS
	int i = 0;
	for (q = start; q != NULL; q = q->next)
		i++;

	return i;
}

int LinkedList_COTS::count_6p_COTS()
{
	// counting every single COTS
	int i = 0;
	for (q = start; q != NULL; q = q->next)
	{
		if (q->age >= 6){
			i++;
		}
	}
	return i;
}

double LinkedList_COTS::count_24p_COTS(double prob, double err, int N0, base_generator_type& eng)
{
	// counting COTS according to Ian's recommendation
	int no = 0;
	double out1;
	Context dist;
	for (q = start; q != NULL; q = q->next){
		if (q->size >= 300){
			no++;
		}
	}
	if (no != 0){
		out1 = dist.lognormal_dis((double)prob*(log(no)-log(N0)), err, eng);
	}
	else{ out1 = 0; }

	return out1;
}

int LinkedList_COTS::count_COTS_m()
{
	// counting number male COTS
	int i = 0;
	if (start != NULL){
		for (q = start; q != NULL; q = q->next){
			if (q->sex == 0){
				i++;
			}
		}
	}
	return i;
}
int LinkedList_COTS::count_COTS_f()
{
	// counting number of female COTS
	int i = 0;
	if (start != NULL){
		for (q = start; q != NULL; q = q->next){
			if (q->sex == 1){
				i++;
			}
		}
	}
	return i;
}
int LinkedList_COTS::count_COTS_sex_mat_m()
{
	// counting number of larvae in the system
	int i = 0;
	if (start != NULL){
		for (q = start; q != NULL; q = q->next){
			if (q->size > LENTH_AT_SEX_MAT && q->sex == 0){
				i++;
			}
		}
	}
	return i;
}
int LinkedList_COTS::count_COTS_sex_mat_f()
{
	// counting number of larvae in the system
	int i = 0;
	if (start != NULL){
		for (q = start; q != NULL; q = q->next){
			if (q->size > LENTH_AT_SEX_MAT && q->sex == 1){
				i++;
			}
		}
	}
	return i;
}
int LinkedList_COTS::count_size10(){
	int i = 0;
	if (start != NULL){
		for (q = start; q != NULL; q = q->next){
			if (q->size > 150){
				i++;
			}
		}
	}
	return i;
}
int LinkedList_COTS::count_size15(){
	int i = 0;
	if (start != NULL){
		for (q = start; q != NULL; q = q->next){
			if (q->size < 150){
				i++;
			}
		}
	}
	return i;
}
int LinkedList_COTS::count_size1525(){
	int i = 0;
	if (start != NULL){
		for (q = start; q != NULL; q = q->next){
			if (q->size >= 150 && q->size < 250){
				i++;
			}
		}
	}
	return i;
}
int LinkedList_COTS::count_size2540(){
	int i = 0;
	if (start != NULL){
		for (q = start; q != NULL; q = q->next){
			if (q->size >= 250 && q->size < 400){
				i++;
			}
		}
	}
	return i;
}
int LinkedList_COTS::count_size40p(){
	int i = 0;
	if (start != NULL){
		for (q = start; q != NULL; q = q->next){
			if (q->size >= 400){
				i++;
			}
		}
	}
	return i;
}

int LinkedList_COTS::count_den_food(){
	int i = 0, j=0;
	if (start != NULL){
		for (q = start; q != NULL; q = q->next){
			if (q->age > 6){
				if (q->size >= 280){
					i++;
				}
				else {
					j += q->size / 280;
				}
			}
		}
	}
	return i+j;

}

struct adult* LinkedList_COTS::insert_initial_pop(int length, int sex, base_generator_type& eng){
	// insert initial COTS population
	// input length
	// input lat and long - 
	// input 
	adult *COTS = new adult;// creat new node pointer 
	adult *tail = NULL;
	COTS->size = length;
	COTS->age = length_age(length, eng);
	COTS->sex = sex;
	COTS->spawn = 0;
	COTS->bs_ind = 0;
	COTS->next = start;
	tail = start;
	start = COTS;
	return start;
}

/*
void LinkedList_COTS::sex_mat(){

int bd = LENTH_AT_SEX_MAT;//r8_uniform_sample(100, 200); just fixed at
for (q = start; q != NULL; q = q->next){
if (q->size >= bd){
q->mat = 1;
}
else q->mat = 0;
}
return;
}
*/
void LinkedList_COTS::reset_spawn(){
	// reset spawning
	for (q = start; q != NULL; q = q->next){
		q->spawn = 0;
	}
	return;
}

long long LinkedList_COTS::spawn(base_generator_type& eng){
	// number of eggs per reef
	// female only spawn once per season 

	adult *q = start;
	long long eggs = 0;
	long long z1;

	int sp_indic;

	Context dist;

	while (q != NULL){

		if (q->sex == 1 && q->size >= LENTH_AT_SEX_MAT && q->spawn == 0 && q->bs_ind==0){
			// if animal is starved prior to spawning season, it will not spawn
			sp_indic = dist.bernoullie_dis(0.5, eng); //
			if (sp_indic == 1){
				z1 = relative_fec_op1(q->size);
				//eggs += poisson_sample_v1(z1);
				if (z1<0){
					cout << q->size << endl;
					cout << "LLLOOO" << endl;
					exit(EXIT_FAILURE);

				}
				if (z1 > 1000){
					eggs += dist.normal_dis(z1, sqrt(z1), eng);
				}
				else{
					eggs += dist.poisson_dis(z1, eng);
				}
				//eggs += simplerng_poisson_large(z1);
				q->spawn = 1;
			}
		}

		q = q->next;
	}
	return eggs*FER_DIS;
}
void LinkedList_COTS::update_age(){
	for (q = start; q != NULL; q = q->next){
		q->age++;
	}
	return;
}
void LinkedList_COTS::nat_mor_A(double max_sur, double k, double x0, double juv_sur, base_generator_type& eng){

	int S; //survival
	int i = 0;
	Context dist;
	double adj_max_sur = max_sur - juv_sur;

	if (start == NULL){
		//cout << "Adult population size is zero" << endl;
	}
	else{
		adult * trail = 0;
		q = start;
		while (q != NULL){
			
			if (q->age > MAX_ADULT_AGE){
				S = dist.bernoullie_dis(0.1, eng);
			}
			else if(q->age > 6 & q->age <= MAX_ADULT_AGE) {
				S = cot_mor(q->size, juv_sur, adj_max_sur, k, x0, eng);
			}
			else{ // age less than 6- will only be used right at the beginning
				if (juv_sur <=0){
					cout << "juv_sur is less than or equal to zero" << endl;
					exit(0);
				}
				S = dist.bernoullie_dis(juv_sur, eng);
			}
			if (S == 0){
				i++;
				if (start == q){
					start = q->next;
					delete q;
					trail = NULL;
					q = start;
				}
				// case 2: beyond the start
				else{
					trail->next = q->next;
					delete q;
					q = trail->next;
				}

			}
			else{
				trail = q;
				q = q->next;
			}
		}

	}
	//cout << "Number of 6mth+ COTS die from natural death: " << i << endl;
	return;
}
void LinkedList_COTS::remove_starved_cots(){

	int i = 0; 

	if (start == NULL){
		//cout << "Adult population size is zero" << endl;
	}
	else{
		adult * trail = 0;
		q = start;
		while (q != NULL){
			if (q->bs_ind >= 5){
				i++;
				if (start == q){
					start = q->next;
					delete q;
					trail = NULL;
					q = start;
				}
				// case 2: beyond the start
				else{
					trail->next = q->next;
					delete q;
					q = trail->next;
				}

			}
			else{
				trail = q;
				q = q->next;
			}
		}

	}
	//cout << "Number of 6mth+ COTS die from starvation: " << i << endl;
	return;
}
void LinkedList_COTS::growth(double val, double beta, double sigma, base_generator_type& eng){
	// growth
	// three phase
	// less than 6 months- coraline algae
	// post algae feeding stage

	for (q = start; q != NULL; q = q->next){
		
		if (q->age <= 6){
			q->size = growth_J(q->size, eng);
		}
		else{
			q->size = growth_A(val, q->size, beta, sigma, eng);

		}
	}
	return;
}

void LinkedList_COTS::growth_feed_step(double food, double beta, double sigma, base_generator_type& eng){
	// growth- step function, 
	// less than 6 months- coraline algae
	// post algae feeding stage

	double old;

	for (q = start; q != NULL; q = q->next){
		if (q->age <= 6){
			q->size = growth_J(q->size, eng);
		}

		else{
			old = q->size;
			q->size = growth_A_withfeed(q->size, beta, sigma, food, eng);

			if (q->size - old <= 0){ // body actually shinking
				if ((old - q->sex) / old > TOL_BODY_SHRINK){
					q->bs_ind = 5;
				}
				q->bs_ind++;
			}
			else{
				q->bs_ind = 0;  //
			}
		}
	}
	return;
}

void LinkedList_COTS::growth_feed_function(double *alpha, double sigma, base_generator_type& eng){
	// growth- step function, 
	// less than 6 months- coraline algae
	// post algae feeding stage

	double old;

	for (q = start; q != NULL; q = q->next){
		if (q->age <= 6){
			q->size = growth_J(q->size, eng);
		}

		else{
			old = q->size;
			q->size = growth_A_withfeed_function(q->size, sigma, alpha, eng);

			if (q->size - old <= 0){ // body actually shinking
				if ((old - q->sex) / old > TOL_BODY_SHRINK){
					q->bs_ind = 5;
				}
				q->bs_ind++;
			}
			else{
				q->bs_ind = 0;  //
			}
		}
	}
	return;
}
void LinkedList_COTS::updatepop(double size, int sex){
	adult *COTS = new adult;
	COTS->age = 6;
	//COTS->mat = 0;
	COTS->bs_ind = 0;
	COTS->spawn = 0;
	COTS->size = size;//trans_length(chdat, no_ch_rec, current_chlo);
	COTS->sex = sex;
	COTS->next = start;
	start = COTS;
}

struct adult* LinkedList_COTS::insert_new_pop(double size, int sex){
	// insert initial COTS population
	// input length
	// input lat and long - 
	// input 
	adult *COTS = new adult;// creat new node pointer 
	adult *tail = NULL;
	COTS->age = 6;
	//	COTS->mat = 0;
	COTS->spawn = 0;
	COTS->bs_ind = 0; 
	COTS->size = size;//trans_length(chdat, no_ch_rec, current_chlo);
	COTS->sex = sex;
	COTS->next = start;
	tail = start;
	start = COTS;
	return start;
}
