#ifndef ODLS_SEQUENTIAL
#define ODLS_SEQUENTIAL

#include <iostream>
#include <map>
#include <vector>
#include <stdlib.h>   
#include <ctime> 
#include <stdio.h>
#include <time.h>
#include <chrono>
#include <set>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <iterator>

typedef std::vector<std::string> dls;

const unsigned LS_ORDER = 10;
const unsigned psuedotriple_char_arr_len = 3 * LS_ORDER * LS_ORDER;
const unsigned NUMBER_OF_COMB = 2445393; // 15953 72356
const int STOP_DUE_NO_DLS = -1;
const int STOP_DUE_LOW_LOCAL_BKV = -2;
//const double WAIT_FIRST_RESULTS_SECONDS = 1800;
//const double WAIT_FINAL_PROCESS_SECONDS = 36000;
const unsigned MAX_DIFF_VALUE_FROM_BKV = 7;
const unsigned NUM_OF_DLS_IN_ONE_CHECK = 1000000;

struct odls_pair
{
	dls dls_1;
	dls dls_2;
};

struct odls_pseudotriple
{
	dls dls_1;
	dls dls_2;
	dls dls_3;
	std::set<std::string> unique_orthogonal_cells;
};

struct fragment_data{
	unsigned orthogonal_value;
	double first_dls_generate_time;
	unsigned long long generated_DLS_count;
	double start_processing_time;
	double end_processing_time;
	short int result;
};

class odls_sequential
{
private:
	double dls_generate_start_time;
	double dls_generate_last_time;
	std::vector<odls_pair> odls_pair_vec;
public:
	odls_sequential();
	odls_pseudotriple best_one_dls_psudotriple;
	odls_pseudotriple best_all_dls_psudotriple;
	odls_pseudotriple dls_psudotriple;
	unsigned long long generated_DLS_count;
	double dls_total_time;
	double pseudotriples_total_time;
	double first_dls_generate_time;
	void readOdlsPairs(std::string known_podls_file_name);
	void makePseudotriple(odls_pair &orthogonal_pair, dls &new_dls, odls_pseudotriple &pseudotriple);
	int deterministicGeneratingDLS(unsigned fragment_index); // Alexey Zhuravlev function
	unsigned short processNewDLS(int fragment_index, unsigned short int *square);
	void odls_sequential::compareLocalRecordWithGlobal(int count, int local_max);
	int odls_sequential::generateDLS(unsigned fragment_index);
};

#endif