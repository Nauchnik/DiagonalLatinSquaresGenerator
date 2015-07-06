#ifdef _MPI
#include <mpi.h>
#endif

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

using namespace std;

typedef std::vector<std::string> dls;

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

void ControlProcess( int rank, int corecount );
void ComputeProcess( int rank, int corecount );
void ReadOdlsPairs( std::vector<odls_pair> &odls_pair_vec );
void MakePseudotriple( odls_pair &orthogonal_pair, dls &new_dls, odls_pseudotriple &pseudotriple );

void SetDefaultValues();
void CalculateFreeNumbers(int I, int J);
void Calculate();
int ChooseValue(int I, int J);
void Retrieve(int & I, int & J);
void RetrieveD(int & I, int & J, int & currentIndex);
void CleanForwardDiagonalBusyNumbers(int I, int J);
bool Forward(int & I, int & J);
void CleanBusyValues();
bool IsDiagonal();
void Print();

int numberOfInts = 8;
unsigned long long countOfCalc;
double limSecondsOneSquare = 100;
bool isTimeBreak = false;
int const numberOfIntsBase = 20;
map<int, int>freeNumbers[numberOfIntsBase][numberOfIntsBase];
map<int, int>busyNumbers[numberOfIntsBase][numberOfIntsBase];
map<int, int>diagonalBusyNumbers[numberOfIntsBase][numberOfIntsBase];
int valuesInTable[numberOfIntsBase][numberOfIntsBase];
odls_pseudotriple best_one_dls_psudotriple, best_all_dls_psudotriple, dls_psudotriple, best_total_pseudotriple;
bool isJustGeneratingDLS = false;

int main(int argc, char **argv)
{
#ifdef _DEBUG
	argc = 4;
	argv[1] = "10";
	argv[2] = "10";
	argv[3] = "10";
	//argv[4] = "pseudotriple_dls_10_template.cnf";
	//argv[5] = "-no_pairs";
#endif;
	int corecount = 0, rank = 1;
	isJustGeneratingDLS = true;
#ifdef _MPI
	MPI_Init( &argc, &argv );
	MPI_Comm_size( MPI_COMM_WORLD, &corecount );
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
#endif
	
	if ( argc < 4 ) {
		std::cerr << "Usage : LS order LS count lim_seconds [pseudotriple_template_cnf_name] -no_pairs";
		return 1;
	}
	
	numberOfInts = atoi(argv[1]);
	countOfCalc = atoi(argv[2]);
	limSecondsOneSquare = atof(argv[3]);
	bool isPairsUsing = true;
	std::string nopairs_str;
	if ( argc == 6 ) {
		std::cout << "argv[5] " << argv[5] << std::endl;
		nopairs_str = argv[5];
		if ( nopairs_str == "-no_pairs" )
			isPairsUsing = false;
	}
	
	std::cout << "isPairsUsing " << isPairsUsing << std::endl;  
	
	std::cout << "numberOfInts " << numberOfInts << std::endl;
	std::cout << "countOfCalc " << countOfCalc << std::endl;
	std::cout << "limSecondsOneSquare " << limSecondsOneSquare << std::endl;
	std::cout << endl;
	std::string str;
	std::vector<odls_pair> odls_pair_vec;

	if ( argc >= 5 ) {
		// checking founded SATisfying assignments mode, instead of pure DLS generating mode
		std::stringstream dls_pair_clauses_sstream, template_clauses_sstream, cells_restr_clause_sstream, tmp_sstream;
		std::string pseudotriple_template_cnf_name = argv[4];
		std::cout << "pseudotriple_template_cnf_name " << pseudotriple_template_cnf_name << std::endl;
		std::ifstream ifile( pseudotriple_template_cnf_name.c_str() );
		if ( !ifile.is_open() ) {
			std::cerr << pseudotriple_template_cnf_name << " not open" << std::endl;
			return 1;
		}
		std::vector<unsigned> cells_restr_var_numbers;
		unsigned uval;
		while ( std::getline( ifile, str ) ) {
			tmp_sstream.str(""); tmp_sstream.clear();
			tmp_sstream << str;
			str.erase( std::remove(str.begin(), str.end(), '\r'), str.end() );
			template_clauses_sstream << str << std::endl;
			if ( cells_restr_var_numbers.size() == 100 )
				continue;
			else
				cells_restr_var_numbers.clear();
			while ( tmp_sstream >> uval ) {
				if ( uval )
					cells_restr_var_numbers.push_back( uval );
			}
		}
		ifile.close();
		if ( cells_restr_var_numbers.size() != 100 ) {
			std::cerr << "cells_restr_var_numbers.size() != 100 ";
			exit(1);
		}
		ReadOdlsPairs( odls_pair_vec );
		std::string cur_pseudotriple_file_name; 
		std::ofstream cur_pseudotriple_file;
		unsigned cells_from = 60, cells_to = 70;
		unsigned pair_index = 0;

		if ( isPairsUsing ) {
			for ( auto &x : odls_pair_vec ) { // for every pair of dls make cnf for searching pseudotriple
				for ( unsigned i=0; i < x.dls_1.size(); i++ )
					for ( unsigned j=0; j < x.dls_1[i].size(); j++ )
						dls_pair_clauses_sstream << 100*i + 10*j + (x.dls_1[i][j]-48)+1 << " 0\n"; // char to int
				for ( unsigned i=0; i < x.dls_2.size(); i++ )
					for ( unsigned j=0; j < x.dls_2[i].size(); j++ ) {
						dls_pair_clauses_sstream << 1*1000 + 100*i + 10*j + (x.dls_2[i][j]-48)+1 << " 0";
						if ( ( i == x.dls_2.size() - 1 ) && ( j == x.dls_2[i].size() - 1 ) )
							dls_pair_clauses_sstream << " "; // for treengling
						else
							dls_pair_clauses_sstream << "\n";
					}
				for ( unsigned i=cells_from; i <= cells_to; i++) {
					cur_pseudotriple_file_name = "dls-pseudotriple_";
					tmp_sstream.clear(); tmp_sstream.str("");
					tmp_sstream << i;
					cur_pseudotriple_file_name += tmp_sstream.str();
					tmp_sstream.clear(); tmp_sstream.str("");
					cur_pseudotriple_file_name += "cells_pair";
					tmp_sstream << pair_index;
					cur_pseudotriple_file_name += tmp_sstream.str();
					cur_pseudotriple_file_name += ".cnf";
					cells_restr_clause_sstream << cells_restr_var_numbers[i-1] << " 0\n";
					cur_pseudotriple_file.open( cur_pseudotriple_file_name.c_str(), std::ios_base::out );
					cur_pseudotriple_file << template_clauses_sstream.str();
					cur_pseudotriple_file << cells_restr_clause_sstream.str();
					cur_pseudotriple_file << dls_pair_clauses_sstream.str();
					cur_pseudotriple_file.close();
					cells_restr_clause_sstream.str(""); cells_restr_clause_sstream.clear();
				}
				dls_pair_clauses_sstream.str(""); dls_pair_clauses_sstream.clear();
				pair_index++;
			}
		}
		else {
			// using no pairs
			for ( unsigned i=cells_from; i <= cells_to; i++) {
				cur_pseudotriple_file_name = "dls-pseudotriple_";
				tmp_sstream.clear(); tmp_sstream.str("");
				tmp_sstream << i;
				cur_pseudotriple_file_name += tmp_sstream.str();
				tmp_sstream.clear(); tmp_sstream.str("");
				cur_pseudotriple_file_name += ".cnf";
				cells_restr_clause_sstream << cells_restr_var_numbers[i-1] << " 0\n";
				cur_pseudotriple_file.open( cur_pseudotriple_file_name.c_str(), std::ios_base::out );
				cur_pseudotriple_file << template_clauses_sstream.str();
				cur_pseudotriple_file << cells_restr_clause_sstream.str();
				cur_pseudotriple_file << dls_pair_clauses_sstream.str();
				cur_pseudotriple_file.close();
				cells_restr_clause_sstream.str(""); cells_restr_clause_sstream.clear();
			}
		}
		
		/*
		// check solution
		stringstream sstream;
		ReadOdlsPairs( odls_pair_vec );
		std::string solutionfile_name = "out_treengeling_dls-pseudotriple_73cells_pair1.cnf";
		std::ifstream solutionfile( solutionfile_name.c_str(), std::ios_base::in );
		std::string str;
		dls new_dls;
		std::string dls_row;
		int val;
		if ( !solutionfile.is_open() ) {
			std::cerr << "solutionfile " << solutionfile_name << " not open" << std::endl;
			return 0;
		}
		
		while ( std::getline( solutionfile, str ) ) {
			if ( ( str[0] == 'v' ) && ( str[1] == ' ' ) ) {
				sstream << str.substr(2);
				while ( sstream >> val ) {
					if ( ( val >= 2001 ) && ( val <= 3000 ) ) {
						val = val % 10 ? (val % 10)-1 : 9;
						dls_row.push_back( '0' + val );
					}
					if ( dls_row.size() == 10 ) {
						new_dls.push_back( dls_row );
						std::cout << dls_row << std::endl;
						dls_row = "";
					}
				}
				sstream.clear(); sstream.str("");
			}
		}
		std::cout << std::endl;
		
		solutionfile.close();
		odls_pseudotriple pseudotriple;
		MakePseudotriple( odls_pair_vec[1], new_dls, pseudotriple );
		std::cout << "pseudotriple.unique_orthogonal_cells.size() " << pseudotriple.unique_orthogonal_cells.size() << std::endl;
		for ( auto &x : pseudotriple.unique_orthogonal_cells )
			std::cout << x << " ";
		 

		// check Brown pseudotriple
		std::cout << std::endl;
		MakePseudotriple( odls_pair_vec[17], odls_pair_vec[18].dls_1, pseudotriple );
		std::cout << "Brown pseudotriple.unique_orthogonal_cells.size() " << pseudotriple.unique_orthogonal_cells.size() << std::endl;
		for ( auto &x : pseudotriple.unique_orthogonal_cells )
			std::cout << x << " ";
		*/
		return 0;
	}
	
	if ( rank == 0 )
		ControlProcess( rank, corecount );
	else
		ComputeProcess( rank, corecount );
	
	return 0;
}

void ControlProcess( int rank, int corecount )
{
	std::cout << "ControlProcess()" << std::endl;
	std::chrono::high_resolution_clock::time_point t1, t2, finding_new_bkv_start_time, now_time;
	std::chrono::duration<double> time_span;
	
	char *psuedotriple_char_arr;
	unsigned psuedotriple_char_arr_len = 3 *numberOfInts * numberOfInts;
	//std::cout << "psuedotriple_char_arr_len " << psuedotriple_char_arr_len << std::endl;
	psuedotriple_char_arr = new char[psuedotriple_char_arr_len];
	std::ofstream ofile;
#ifdef _MPI
	MPI_Status status;
#endif
	dls new_dls;
	std::stringstream out_sstream;
	unsigned char_index;
	odls_pair cur_pair;
	cur_pair.dls_1.resize( numberOfInts );
	for( auto &x : cur_pair.dls_1 )
		x.resize( numberOfInts );
	cur_pair.dls_2.resize( numberOfInts );
	for( auto &x : cur_pair.dls_2 )
		x.resize( numberOfInts );
	new_dls.resize( numberOfInts );
	for( auto &x : new_dls )
		x.resize( numberOfInts );
	
	unsigned orthogonal_value_from_message = 0;
	t1 = std::chrono::high_resolution_clock::now();
	for ( ;; ) {
#ifdef _MPI
		//std::cout << "process " << rank << " before recieving" << std::endl;
		MPI_Recv( &orthogonal_value_from_message, 1, MPI_UNSIGNED, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status );
		//std::cout << "process " << rank << " recieved " << orthogonal_value_from_message << std::endl;
		MPI_Recv( psuedotriple_char_arr, psuedotriple_char_arr_len, MPI_CHAR, status.MPI_SOURCE, 0, MPI_COMM_WORLD, &status );		
		//std::cout << "process " << rank << " after recieving" << std::endl;
		/*std::cout << "psuedotriple_char_arr " << std::endl;
		for ( unsigned i = 0; i < psuedotriple_char_arr_len; i++ )
			std::cout << psuedotriple_char_arr[i];*/
		//std::cout << std::endl;
		out_sstream << "completed message from process " << status.MPI_SOURCE << std::endl;
		//std::cout << out_sstream.str();
#endif
		char_index = 0;
		for ( int i = 0; i < numberOfInts; i++ )
			for ( int j = 0; j < numberOfInts; j++ )
				cur_pair.dls_1[i][j] = psuedotriple_char_arr[char_index++];
		for ( int i = 0; i < numberOfInts; i++ )
			for ( int j = 0; j < numberOfInts; j++ )
				cur_pair.dls_2[i][j] = psuedotriple_char_arr[char_index++];
		for ( int i = 0; i < numberOfInts; i++ )	
			for ( int j = 0; j < numberOfInts; j++ )
				new_dls[i][j] = psuedotriple_char_arr[char_index++];
		MakePseudotriple( cur_pair, new_dls, dls_psudotriple );
		out_sstream << "dls_psudotriple.unique_orthogonal_cells.size() " << dls_psudotriple.unique_orthogonal_cells.size() << std::endl;
		if ( orthogonal_value_from_message != dls_psudotriple.unique_orthogonal_cells.size() ) {
			std::cerr << "error. orthogonal_value_from_message != dls_psudotriple.unique_orthogonal_cells.size()" << std::endl;
			std::cerr << orthogonal_value_from_message << " != " << dls_psudotriple.unique_orthogonal_cells.size() << std::endl;
#ifdef _MPI
			MPI_Abort( MPI_COMM_WORLD, 0 );
#endif
		}
		else {
			if ( dls_psudotriple.unique_orthogonal_cells.size() > best_total_pseudotriple.unique_orthogonal_cells.size() ) {
				best_total_pseudotriple = dls_psudotriple;
				t2 = std::chrono::high_resolution_clock::now();
				time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
				t1 = t2;
				out_sstream << std::endl << "new total_bkv " << best_total_pseudotriple.unique_orthogonal_cells.size() << std::endl;
				for ( auto &x : best_total_pseudotriple.dls_1 ) {
					for ( auto &y : x )
						out_sstream << y << " ";
					out_sstream << std::endl;
				}
				out_sstream << std::endl;
				for ( auto &x : best_total_pseudotriple.dls_2 ) {
					for ( auto &y : x )
						out_sstream << y << " ";
					out_sstream << std::endl;
				}
				out_sstream << std::endl;
				for ( auto &x : best_total_pseudotriple.dls_3 ) {
					for ( auto &y : x )
						out_sstream << y << " ";
					out_sstream << std::endl;
				}
				out_sstream << std::endl;
				for ( auto &x : best_total_pseudotriple.unique_orthogonal_cells )
					out_sstream << x << " ";
				out_sstream << std::endl;
				
				out_sstream << "time from previous BKV " << time_span.count() << std::endl << std::endl;
				ofile.open( "out", std::ios_base::app );
				ofile << out_sstream.str();
				ofile.close();
				out_sstream.clear(); out_sstream.str("");
			}
		}
	}
	delete[] psuedotriple_char_arr;
}

void ComputeProcess( int rank, int corecount )
{
	double sum_time = 0.0, min_time = 604800.0, max_time = 0.0;
	std::chrono::high_resolution_clock::time_point t1, t2, start_time, now_time;
	std::chrono::duration<double> time_span;
	double solving_time;

	// set of unique diagonal Latin squares
	std::set<std::string> dls_string_set;
	std::string cur_string_set;
	unsigned k;
	std::stringstream sstream;
	unsigned long long genereated_count = 0;
	unsigned best_first_pair_orthogonal_cells = 0, best_second_pair_orthogonal_cells = 0;
	dls new_dls;
	new_dls.resize( numberOfInts );
	start_time = std::chrono::high_resolution_clock::now();
	char *psuedotriple_char_arr;
	unsigned psuedotriple_char_arr_len = 3 *numberOfInts * numberOfInts;
	//std::cout << "psuedotriple_char_arr_len " << psuedotriple_char_arr_len << std::endl;
	psuedotriple_char_arr = new char[psuedotriple_char_arr_len];
	unsigned char_index;

	std::vector<odls_pair> odls_pair_vec;
	ReadOdlsPairs( odls_pair_vec );
	unsigned ortogonal_cells;

	// at first check known DLS's from file
	unsigned preprocess_bkv = 0;
	std::set<dls> dls_known;
	for ( auto &x : odls_pair_vec ) {
		dls_known.insert( x.dls_1 );
		dls_known.insert( x.dls_2 );
	}
	dls tmp_dls;
	odls_pseudotriple psudotriple;
	for (auto &x : odls_pair_vec) {
		for (auto &y : dls_known) {
			tmp_dls = y;
			if ((tmp_dls != x.dls_1) && (tmp_dls != x.dls_2)) {
				MakePseudotriple(x, tmp_dls, psudotriple);
				if (psudotriple.unique_orthogonal_cells.size() >= preprocess_bkv) {
					preprocess_bkv = psudotriple.unique_orthogonal_cells.size();
					std::cout << "preprocess_bkv " << preprocess_bkv << std::endl;
				}
				// check psuedotripy by Brown (2 pairs are orthogonal)
				if ((x.dls_1[0] == "0946175823") &&
					(x.dls_2[0] == "0851734692") &&
					(tmp_dls[0] == "0419827356")) {
					std::cout << "Brown pseudotriple BKV " << psudotriple.unique_orthogonal_cells.size() << std::endl;
					for (auto &x : psudotriple.unique_orthogonal_cells)
						std::cout << x << " ";
					std::cout << std::endl;
				}
			}
		}
	}
	
	std::cout << "preprocess_bkv based on knwon DLS from input file : " << preprocess_bkv << std::endl;

	for (unsigned long long i = 0; i < countOfCalc; i++) {
		t1 = std::chrono::high_resolution_clock::now();

		Calculate();

		t2 = std::chrono::high_resolution_clock::now();
		time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
		solving_time = time_span.count();
		
		if ( isTimeBreak ) { // if process was interrupted, start new process
			isTimeBreak = false;
			continue;
		}

		// here we have diagonal Latin square, let's add it to the set
		genereated_count++;
		k = 0;
		cur_string_set = "";
		
		for (int j1 = 0; j1 < numberOfInts; j1++) {
			for (int j2 = 0; j2 < numberOfInts; j2++)
				sstream << valuesInTable[j1][j2] - 1; // save values from 0 to 9
			new_dls[j1] = sstream.str();
			cur_string_set += sstream.str();
			sstream.clear(); sstream.str("");
		}
		//std::cout << cur_string_set << endl;
		//dls_string_set.insert( cur_string_set );
		
		/*std::cout << "new_dls" << std::endl;
		for ( auto &x : new_dls ) {
			for ( auto &y : x )
				std::cout << y << " ";
			std::cout << std::endl;
		}*/
		
		if (!isJustGeneratingDLS) {
			// make pseudotriple for every pair and choose one with maximum of orthogonal cells
			for( auto &x : odls_pair_vec ) {
				MakePseudotriple( x, new_dls, dls_psudotriple );
				if ( dls_psudotriple.unique_orthogonal_cells.size() > best_one_dls_psudotriple.unique_orthogonal_cells.size() )
					best_one_dls_psudotriple = dls_psudotriple;
				//std::cout << "cur_first_pair_orthogonal_cells " << cur_first_pair_orthogonal_cells << std::endl;
				//std::cout << "cur_second_pair_orthogonal_cells " << cur_second_pair_orthogonal_cells << std::endl;
				//std::cout << "cur_one_dls_psudotriple_orthogonal_cells " << cur_one_dls_psudotriple_orthogonal_cells << std::endl;
				//std::cout << "best_one_dls_psudotriple_orthogonal_cells " << best_one_dls_psudotriple_orthogonal_cells << std::endl;
			}
		
			if (best_one_dls_psudotriple.unique_orthogonal_cells.size() > best_all_dls_psudotriple.unique_orthogonal_cells.size()) {
				best_all_dls_psudotriple = best_one_dls_psudotriple;
				std::cout << "best_all_dls_psudotriple_orthogonal_cells " << best_all_dls_psudotriple.unique_orthogonal_cells.size() << std::endl;
				now_time = chrono::high_resolution_clock::now();
				time_span = std::chrono::duration_cast<std::chrono::duration<double>>(now_time - start_time);
				std::cout << "time from start " << time_span.count() << std::endl;
				std::cout << "genereated_count " << genereated_count << std::endl;
				char_index = 0;
				for (auto &x : best_all_dls_psudotriple.dls_1)
					for (auto &y : x)
						psuedotriple_char_arr[char_index++] = y;
				for (auto &x : best_all_dls_psudotriple.dls_2)
					for (auto &y : x)
						psuedotriple_char_arr[char_index++] = y;
				for (auto &x : best_all_dls_psudotriple.dls_3)
					for (auto &y : x)
						psuedotriple_char_arr[char_index++] = y;
				ortogonal_cells = best_all_dls_psudotriple.unique_orthogonal_cells.size();

	#ifdef _MPI
				//std::cout << "process " << rank << " before sending" << std::endl;
				MPI_Send( &ortogonal_cells, 1, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD );
				//std::cout << "process " << rank << " after sending" << std::endl;
				MPI_Send( psuedotriple_char_arr, psuedotriple_char_arr_len, MPI_CHAR, 0, 0, MPI_COMM_WORLD );
				//std::cout << "process " << rank << " after sending" << std::endl;
	#endif
			}
		}
		
		//Print();
		
		min_time = min_time < solving_time ? min_time : solving_time;
		max_time = max_time > solving_time ? max_time : solving_time;
		sum_time += solving_time;
		/*
		cout << "time in seconds : " << solving_time << endl;
		cout << "processed " << i+1 << " from " << countOfCalc << endl;
		cout << "generated " << genereated_count << endl;
		cout << "unique DLS : " << dls_string_set.size() << endl;
		cout << "current median time " << sum_time / double(i + 1) << endl;
		cout << "current min_time " << min_time << endl;
		cout << "current max_time " << max_time << endl;
		cout << endl;*/
	}
	std::cout << "genereated_count " << genereated_count << std::endl;

	delete[] psuedotriple_char_arr;
}

void MakePseudotriple( odls_pair &orthogonal_pair, dls &new_dls, odls_pseudotriple &pseudotriple )
{
	unsigned cur_first_pair_orthogonal_cells, cur_second_pair_orthogonal_cells;
	std::set<std::string> greece_latin_square1, greece_latin_square2; // for counting orthogonal cells between 2 DLS
	std::string cell_plus_cell;
	
	for ( unsigned j1=0; j1< new_dls.size(); j1++)
		for ( unsigned j2=0; j2 < new_dls[j1].size(); j2++) {
			cell_plus_cell = orthogonal_pair.dls_1[j1][j2];
			cell_plus_cell += new_dls[j1][j2]; 
			greece_latin_square1.insert( cell_plus_cell );
		}
	cur_first_pair_orthogonal_cells = greece_latin_square1.size();
	/*std::cout << "greece_latin_square1 " << std::endl;
	for ( auto &y : greece_latin_square1 )
		std::cout << y << " ";
	std::cout << std::endl;*/
	for ( unsigned j1=0; j1 < new_dls.size(); j1++)
		for ( unsigned j2=0; j2 < new_dls[j1].size(); j2++) {
			cell_plus_cell = orthogonal_pair.dls_2[j1][j2];
			cell_plus_cell += new_dls[j1][j2];
			greece_latin_square2.insert( cell_plus_cell );
		}
	cur_second_pair_orthogonal_cells = greece_latin_square2.size();
	pseudotriple.dls_1 = orthogonal_pair.dls_1;
	pseudotriple.dls_2 = orthogonal_pair.dls_2;
	pseudotriple.dls_3 = new_dls;
	pseudotriple.unique_orthogonal_cells.clear();
	std::set_intersection( greece_latin_square1.begin(),greece_latin_square1.end(),
		                   greece_latin_square2.begin(),greece_latin_square2.end(),
						   std::inserter(pseudotriple.unique_orthogonal_cells, pseudotriple.unique_orthogonal_cells.begin()));
}

void Calculate()
{
	SetDefaultValues();

	int stopedIndex = -10;
	int diagonalStopedIndex = -10;
	int currentIndex = -1;
	bool reachStoped = false;
	int i = 0;
	int j = 0;
	double cur_solving_time_sec;
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	std::chrono::high_resolution_clock::time_point t2;

	CalculateFreeNumbers(0, 0);
	int value;

	while (true) {
		//Вышли из проблемной зоны
		if (currentIndex > stopedIndex) {
			reachStoped = false;
			CleanBusyValues();
		}

		value = ChooseValue(i, j);

		////Не осталось свободного значения
		//if (value == -1 & !IsDiagonal())
		//{
		//	cout << "Hello" << endl;
		//	Retrieve(i, j);
		//
		//	RetrieveD(i, j, currentIndex);
		//
		//	Фиксация крайней точки в проблемной зоне
		//	
		//	if (!reachStoped)
		//	{
		//		reachStoped = true;
		//		stopedIndex = currentIndex;
		//	}
		//	currentIndex--;
		//}
		//
		//else

		if (value == -1 || !IsDiagonal())
		{
			//возвращаемся к предыдущей ячейке
			Retrieve(i, j);

			value = valuesInTable[i][j];
			busyNumbers[i][j].insert(pair<int, int>(value, value));
			
			valuesInTable[i][j] = 0;

			//Фиксация крайней точки в проблемной зоне
			if (!reachStoped) {
				reachStoped = true;
				stopedIndex = currentIndex;
			}
			currentIndex--;

		}
		else {
			currentIndex++;

			valuesInTable[i][j] = value;

			if (Forward(i, j))
				break;

		/*	if (i > 0 && i == j)
				diagonalBusyNumbers[i - 1][j - 1].clear();

			if (i > 0 && i == numberOfInts - j - 1)
				diagonalBusyNumbers[i - 1][j + 1].clear();*/

			busyNumbers[i][j].clear();
		}
		
		CalculateFreeNumbers(i, j);
		
		t2 = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
		cur_solving_time_sec = time_span.count();
		if ( cur_solving_time_sec >= limSecondsOneSquare ) {
			//std::cout << "limSecondsOneSquare " << limSecondsOneSquare << std::endl << std::endl;
			isTimeBreak = true;
			break;
		}
	}
}

void ReadOdlsPairs( std::vector<odls_pair> &odls_pair_vec )
{
	std::string odls_pairs_file_name = "ODLS_10_pairs.txt";
	std::ifstream odls_pairs_file( odls_pairs_file_name );
	if ( !odls_pairs_file.is_open() ) {
		std::cerr << odls_pairs_file_name << " is not open" << std::endl;
		exit(1);
	}
	std::string str, nonspace_str;
	odls_pair_vec.resize( 0 );
	odls_pair cur_odls_pair;
	while ( getline( odls_pairs_file, str ) ) {
		if ( ( str.size() <= 1 ) && ( cur_odls_pair.dls_1.size() > 0 ) ) { // if separate string
			odls_pair_vec.push_back( cur_odls_pair ); // add current odls pair
			cur_odls_pair.dls_1.clear();
			cur_odls_pair.dls_2.clear();
		}
		else { 
			// read two rows of 2 DLS for current pair
			nonspace_str = "";
			for ( auto &x : str )
				if ( ( x != ' ' ) && ( x != '\r' ) ) 
					nonspace_str += x;
			if ( nonspace_str.size() != 2*numberOfInts ) {
				std::cerr << "nonspace_str.size() != 2*numberOfInts" << std::endl;
				std::cerr << nonspace_str.size() << " != " << 2*numberOfInts << std::endl;
				exit(1);
			}
			cur_odls_pair.dls_1.push_back( nonspace_str.substr( 0, numberOfInts ) );
			cur_odls_pair.dls_2.push_back( nonspace_str.substr( numberOfInts, numberOfInts ) );
		}
	}
	odls_pairs_file.close();
	// if there is no separate string at the end of file, add last pair manually
	if ( cur_odls_pair.dls_1.size() != 0 )
		odls_pair_vec.push_back( cur_odls_pair ); // add current odls pair
	
	/*for ( auto &x : odls_pair_vec ) {
		for ( auto &y : x.dls_1 )
			std::cout << y << std::endl;
		std::cout << std::endl;
		for ( auto &z : x.dls_2 )
			std::cout << z << std::endl;
		std::cout << std::endl << std::endl;
	}*/
	
	std::set<std::string> greece_latin_square;
	string cell_plus_cell;
	// check every pair
	for ( auto &x : odls_pair_vec ) {
		for ( unsigned j1=0; j1 < x.dls_1.size(); j1++)
			for ( unsigned j2=0; j2 < x.dls_1[j1].size(); j2++) {
				cell_plus_cell =  x.dls_1[j1][j2];
				cell_plus_cell += x.dls_2[j1][j2];
				greece_latin_square.insert( cell_plus_cell );
			}
			//std::cout << "greece_latin_square.size() " << greece_latin_square.size() << std::endl;
			if ( greece_latin_square.size() != x.dls_1.size() * x.dls_1.size()  ) {
				std::cerr << "greece_latin_square.size() != x.dls_1.size() * x.dls_1.size() " << std::endl;
				std::cerr << greece_latin_square.size() << " != " << x.dls_1.size() * x.dls_1.size() << std::endl;
			}
	}
}

void SetDefaultValues()
{
	for (int i = 0; i < numberOfInts; i++)
		for (int j = 0; j < numberOfInts; j++)
			valuesInTable[i][j] = 0;
}

void CalculateFreeNumbers(int I, int J)
{
	for (int i = 1; i <= numberOfInts; i++)
		freeNumbers[I][J].insert(make_pair(i, i));

	int value;
	for (int i = 0; i < I; i++) {
		value = valuesInTable[i][J];
		freeNumbers[I][J].erase(value);
	}

	for (int j = 0; j < J; j++) {
		value = valuesInTable[I][j];
		freeNumbers[I][J].erase(value);
	}

	/*for (int i = 1; i <= numberOfInts; i++)
	{
		auto valueBusy = diagonalBusyNumbers[I][J].find(i);

		if (valueBusy != diagonalBusyNumbers[I][J].end())
		{
			freeNumbers[I][J].erase(valueBusy->second);
		}
	}*/

	for (int i = 1; i <= numberOfInts; i++) {
		auto valueBusy_it = busyNumbers[I][J].find(i);
		if (valueBusy_it != busyNumbers[I][J].end())
			freeNumbers[I][J].erase(valueBusy_it->second);
	}

	if (I != J && I != numberOfInts - J - 1)
		return;

	int valueBusy;
	if ( I == J)
		for (int i = 0; i < I; i++) {
			valueBusy = valuesInTable[i][i];
			freeNumbers[I][J].erase(valueBusy);
		}

	if (I == numberOfInts - J - 1)
		for (int i = 0; i < I; i++) {
			valueBusy = valuesInTable[i][numberOfInts - i - 1];
			freeNumbers[I][J].erase(valueBusy);
		}
}

int ChooseValue(int I, int J)
{
	vector <int> freeValues;

	for (int i = 1; i <= numberOfInts; i++) {
		auto value = freeNumbers[I][J].find(i);
		if (value != freeNumbers[I][J].end())
			freeValues.push_back(value->second);
	}

	if (freeValues.size() == 0)
		return -1;

	int k = rand() % freeValues.size();

	return freeValues[k];
}

//void RetrieveD(int & i, int & j, int & currentIndex)
//{
//	CleanForwardDiagonalBusyNumbers(i, j);
//	
//	if (i == j & i == numberOfInts - j - 1)
//		;
//
//	else
//	if (i == j)
//	{
//		do
//		{
//			busyNumbers[i][j].clear();
//
//			valuesInTable[i][j] = 0;
//
//			Retrieve(i, j);
//
//			currentIndex--;
//
//		} while (i != j);
//
//	}
//
//	else
//	if (i == numberOfInts - j - 1)
//	{
//		do
//		{
//			busyNumbers[i][j].clear();
//
//			valuesInTable[i][j] = 0;
//
//			Retrieve(i, j);
//
//			currentIndex--;
//
//		} while (i != numberOfInts - j - 1);
//
//	}
//		
//	int value = valuesInTable[i][j];
//
//	diagonalBusyNumbers[i][j].insert(pair<int, int>(value, value));
//
//	valuesInTable[i][j] = 0;
//
//	currentIndex--;
//	
//}

//void CleanForwardDiagonalBusyNumbers(int I, int J)
//{
//	for (int j = J; j < numberOfInts; j++)
//	{
//		diagonalBusyNumbers[I][j].clear();
//	}
//
//	for (int i = I + 1; i < numberOfInts; i++)
//	{
//		for (int j = 0; j < numberOfInts; j++)
//		{
//			diagonalBusyNumbers[i][j].clear();
//		}
//	}
//	
//}

void Retrieve(int & i, int & j)
{
	if (j > 0)
		j--;
	else {
		i--;
		j = numberOfInts - 1;
	}
}

bool Forward(int & i, int & j)
{
	if (j < numberOfInts - 1)
		j++;
	else {
		i++;
		j = 0;
	}
	return i == numberOfInts && j == 0;
}

void CleanBusyValues()
{
	for (int i = 0; i < numberOfInts; i++)
		for (int j = 0; j < numberOfInts; j++)
			busyNumbers[i][j].clear();
}

bool IsDiagonal()
{
	int a[numberOfIntsBase];

	for (int i = 0; i < numberOfInts; i++)
		a[i] = 0;

	int value;
	for (int i = 0; i < numberOfInts; i++) {
		value = valuesInTable[i][i];
		if (value == 0)
			break;
		a[value - 1]++;
	}

	for (int i = 0; i < numberOfInts; i++)
		if (a[i] > 1) return false;

	for (int i = 0; i < numberOfInts; i++)
		a[i] = 0;

	for (int i = 0; i < numberOfInts; i++) {
		value = valuesInTable[i][numberOfInts - i - 1];
		if (value == 0)
			break;
		a[value - 1]++;
	}

	for (int i = 0; i < numberOfInts; i++)
		if (a[i] > 1) return false;

	return true;
}

void Print()
{
	for (int i = 0; i < numberOfInts; i++) {
		for (int j = 0; j < numberOfInts; j++)
			printf("%3d", valuesInTable[i][j]);
		cout << endl;
	}
	cout << endl;
}