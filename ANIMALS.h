#ifndef LinkedListANIMALS_h
#define LinkedListANIMALS_h

#include<stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <iomanip>
#include "struct.h"
#include "Misc.h"
#include "DefineParameters.h"
#include "Dist.h"

class LinkedList_COTS
{
public:

	LinkedList_COTS();// set up default values
	~LinkedList_COTS();
	bool isEmpty_COTS();
	void makeEmpty_COTS();
	void print_COTS();
	int count_COTS();
	int count_6p_COTS();
	double count_24p_COTS(double prob, double err, int N0, base_generator_type& eng);
	int count_COTS_f();
	int count_COTS_m();
	int count_COTS_sex_mat_m();
	int count_COTS_sex_mat_f();
	int count_size10();
	int count_size15();
	int count_size1525();
	int count_size2540();
	int count_size40p();
	int count_den_food();
	struct adult* insert_initial_pop(int length, int sex, base_generator_type& eng);
	void sex_mat(); 
	void reset_spawn();
	long long spawn(base_generator_type& eng);
	void update_age();
	void remove_starved_cots();
	void nat_mor_A(double max_sur, double k, double x0, double juv_sur, base_generator_type& eng);
	void growth(double val, double beta, double sigma, base_generator_type& eng);
	void growth_feed_step(double food, double beta, double sigma, base_generator_type& eng);
	void growth_feed_function(double *alpha, double sigma, base_generator_type& eng);
	void updatepop(double size, int age);
	struct adult* insert_new_pop(double size, int sex);
	adult *start;
private:
	adult *q;
};

#endif 