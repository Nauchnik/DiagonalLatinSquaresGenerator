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

#include "odls_sequential.h"

odls_pseudotriple best_total_pseudotriple;

void controlProcess(int rank, int corecount, odls_sequential odls_seq);
void computingProcess(int rank, int corecount, odls_sequential odls_seq);
void updateFragmentFile(std::vector<fragment_data> total_fragment_data, double mpi_start_time);

int main(int argc, char **argv)
{
#ifdef _DEBUG
	argc = 1;
	//argv[1] = "pseudotriple_dls_10_template.cnf";
	//argv[2] = "-no_pairs";
#endif;
	int corecount = 0, rank = 1;
	odls_sequential odls_seq;
	//isJustGeneratingDLS = true;
	
	MPI_Init( &argc, &argv );
	MPI_Comm_size( MPI_COMM_WORLD, &corecount );
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );

	std::cout << "corecount " << corecount << std::endl;

	if ( argc > 3 ) {
		std::cerr << "Usage : [pseudotriple_template_cnf_name] [-no_pairs]";
		return 1;
	}
	
	std::string str;
	
	if (argc >= 2) {
		std::string pseudotriple_template_cnf_name = argv[1];
		std::cout << "constructPseudotripleCNFs()" << std::endl;
		bool isPairsUsing = true;
		std::string nopairs_str;
		if (argc == 3) {
			std::cout << "argv[2] " << argv[2] << std::endl;
			nopairs_str = argv[2];
			if (nopairs_str == "-no_pairs")
				isPairsUsing = false;
		}
		std::cout << "isPairsUsing " << isPairsUsing << std::endl;
		odls_seq.constructPseudotripleCNFs(pseudotriple_template_cnf_name, isPairsUsing);
	}
	else {
		// MPI searching for DLS and constucting pseudotriples
		std::cout << "MPI searching for DLS and constucting pseudotriples" << std::endl;
		if (rank == 0)
			controlProcess(rank, corecount, odls_seq);
		else
			computingProcess(rank, corecount, odls_seq);
	}
	
	return 0;
}

void controlProcess(int rank, int corecount, odls_sequential odls_seq)
{
	double mpi_start_time = MPI_Wtime();
	std::cout << "ControlProcess()" << std::endl;
	std::chrono::high_resolution_clock::time_point t1, t2, finding_new_bkv_start_time, now_time, total_start_time;
	std::chrono::duration<double> time_span;
	total_start_time = std::chrono::high_resolution_clock::now();
	
	if (number_of_comb < corecount - 1) {
		std::cerr << "number_of_comb < corecount - 1" << std::endl;
		std::cerr << number_of_comb << " < " << corecount - 1 << std::endl;
		MPI_Abort(MPI_COMM_WORLD, 0);
	}
	
	MPI_Status status;
	MPI_Request request;
	mpi_start_time = MPI_Wtime();
	int fragment_index_to_send = corecount - 1;
	
	char psuedotriple_char_arr[psuedotriple_char_arr_len];
	std::ofstream ofile;
	std::vector<fragment_data> total_fragment_data;
	total_fragment_data.resize(number_of_comb);
	for (auto &x : total_fragment_data) {
		x.first_dls_generate_time = 0;
		x.generated_DLS_count = 0;
		x.orthogonal_value = 0;
		x.start_processing_time = 0.0;
		x.end_processing_time = 0.0;
		x.result = 0;
	}
	
	// send first tasks to all compute processes
	for (int i = 0; i < corecount - 1; i++) {
		MPI_Send(&i, 1, MPI_INT, i + 1, 0, MPI_COMM_WORLD);
		total_fragment_data[i].start_processing_time = MPI_Wtime();
	}
	
	dls new_dls;
	std::stringstream out_sstream;
	unsigned char_index;
	odls_pair cur_pair;
	cur_pair.dls_1.resize(LS_order);
	for( auto &x : cur_pair.dls_1 )
		x.resize(LS_order);
	cur_pair.dls_2.resize(LS_order);
	for( auto &x : cur_pair.dls_2 )
		x.resize(LS_order);
	new_dls.resize(LS_order);
	for( auto &x : new_dls )
		x.resize(LS_order);
	
	int orthogonal_value_from_computing_process = 0;
	int fragment_index_from_computing_process = 0;
	unsigned solved_tasks_count = 0;
	int current_bkv_neg = -100;
	int result_from_computing_process = 0;
	unsigned no_dls_stopped_count = 0, low_local_bkv_stopped_count = 0;
	t1 = std::chrono::high_resolution_clock::now();

	// start of receiving results and sending new tasks instead
	while (solved_tasks_count < number_of_comb) {
		//std::cout << "process " << rank << " before recieving" << std::endl;
		MPI_Recv( &fragment_index_from_computing_process, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
		MPI_Recv( &orthogonal_value_from_computing_process, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD, &status);
		
		// if processing of a task was interrupted, then send new task on a freed process
		if (( orthogonal_value_from_computing_process == STOP_DUE_NO_DLS) ||
			( orthogonal_value_from_computing_process == STOP_DUE_LOW_LOCAL_BKV))
		{
			result_from_computing_process = orthogonal_value_from_computing_process;

			// if a result from a computing process was received already, do nothing
			if (total_fragment_data[fragment_index_from_computing_process].result != 0)
				continue;
			
			solved_tasks_count++;
			out_sstream << "solved_tasks_count " << solved_tasks_count << std::endl;
			std::cout << "solved_tasks_count " << solved_tasks_count << std::endl;
			
			total_fragment_data[fragment_index_from_computing_process].end_processing_time = MPI_Wtime();
			total_fragment_data[fragment_index_from_computing_process].result = result_from_computing_process;
			
			no_dls_stopped_count = low_local_bkv_stopped_count = 0;
			for (auto &x : total_fragment_data) {
				if (x.result == STOP_DUE_NO_DLS)
					no_dls_stopped_count++;
				else if (x.result == STOP_DUE_LOW_LOCAL_BKV)
					low_local_bkv_stopped_count++;
			}
			out_sstream << "no_dls_stopped_count " << no_dls_stopped_count << std::endl;
			out_sstream << "low_local_bkv_stopped_count " << low_local_bkv_stopped_count << std::endl;
			
			if (solved_tasks_count != no_dls_stopped_count + low_local_bkv_stopped_count) {
				std::cerr << "solved_tasks_count != no_dls_stopped_count + low_local_bkv_stopped_count" << std::endl;
				std::cerr << solved_tasks_count << " != " << no_dls_stopped_count + low_local_bkv_stopped_count << std::endl;
				MPI_Abort(MPI_COMM_WORLD, 0);
			}
			
			// send new task if there are free ones
			if ( fragment_index_to_send < number_of_comb ) {
				MPI_Send(&fragment_index_to_send, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
				total_fragment_data[fragment_index_to_send].start_processing_time = MPI_Wtime();
				out_sstream << "fragment_index_to_send " << fragment_index_to_send << " was sent" << std::endl;
				out_sstream << std::endl;
				fragment_index_to_send++;
			}
			ofile.open("out", std::ios_base::app);
			ofile << out_sstream.str();
			ofile.close();
			out_sstream.clear(); out_sstream.str("");
			updateFragmentFile(total_fragment_data, mpi_start_time);
			continue;
		}
		else if ((orthogonal_value_from_computing_process > 0) && (orthogonal_value_from_computing_process <= LS_order*LS_order)) { // if new local BKV was received
			//std::cout << "process " << rank << " recieved " << orthogonal_value_from_message << std::endl;
			MPI_Recv(psuedotriple_char_arr, psuedotriple_char_arr_len, MPI_CHAR, status.MPI_SOURCE, 0, MPI_COMM_WORLD, &status);
			MPI_Recv(&odls_seq.first_dls_generate_time, 1, MPI_DOUBLE, status.MPI_SOURCE, 0, MPI_COMM_WORLD, &status);
			MPI_Recv(&odls_seq.generated_DLS_count, 1, MPI_UNSIGNED_LONG_LONG, status.MPI_SOURCE, 0, MPI_COMM_WORLD, &status);
			total_fragment_data[fragment_index_from_computing_process].first_dls_generate_time = odls_seq.first_dls_generate_time;
			total_fragment_data[fragment_index_from_computing_process].generated_DLS_count = odls_seq.generated_DLS_count;
			total_fragment_data[fragment_index_from_computing_process].orthogonal_value = orthogonal_value_from_computing_process;
			updateFragmentFile(total_fragment_data, mpi_start_time);
		}
		else {
			std::cerr << " incorrect orthogonal_value_from_computing_process value " << orthogonal_value_from_computing_process << std::endl;
			MPI_Abort(MPI_COMM_WORLD, 0);
		}
		
		char_index = 0;
		for (int i = 0; i < LS_order; i++)
			for (int j = 0; j < LS_order; j++)
				cur_pair.dls_1[i][j] = psuedotriple_char_arr[char_index++];
		for (int i = 0; i < LS_order; i++)
			for (int j = 0; j < LS_order; j++)
				cur_pair.dls_2[i][j] = psuedotriple_char_arr[char_index++];
		for (int i = 0; i < LS_order; i++)
			for (int j = 0; j < LS_order; j++)
				new_dls[i][j] = psuedotriple_char_arr[char_index++];
		odls_seq.makePseudotriple(cur_pair, new_dls, odls_seq.dls_psudotriple);
		out_sstream << "dls_psudotriple.unique_orthogonal_cells.size() " << odls_seq.dls_psudotriple.unique_orthogonal_cells.size() << std::endl;
		
		if ( orthogonal_value_from_computing_process != odls_seq.dls_psudotriple.unique_orthogonal_cells.size() ) {
			std::cerr << "error. orthogonal_value_from_message != dls_psudotriple.unique_orthogonal_cells.size()" << std::endl;
			std::cerr << orthogonal_value_from_computing_process << " != " << odls_seq.dls_psudotriple.unique_orthogonal_cells.size() << std::endl;
			MPI_Abort( MPI_COMM_WORLD, 0 );
		}
		
		if (odls_seq.dls_psudotriple.unique_orthogonal_cells.size() > best_total_pseudotriple.unique_orthogonal_cells.size()) {
			// new BKV found
			best_total_pseudotriple = odls_seq.dls_psudotriple;

			// send new BKV to all computing processes
			// send negative value to make difference between bkv and new search space index
			current_bkv_neg = -(int)(best_total_pseudotriple.unique_orthogonal_cells.size()); 
			std::cout << "current_bkv_neg " << current_bkv_neg << std::endl;
			for (int i = 0; i < corecount - 1; i++)
				MPI_Isend(&current_bkv_neg, 1, MPI_INT, i + 1, 0, MPI_COMM_WORLD, &request); 
			
			t2 = std::chrono::high_resolution_clock::now();
			time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
			t1 = t2;
			out_sstream << std::endl << "new total_bkv " << best_total_pseudotriple.unique_orthogonal_cells.size() << std::endl;
			out_sstream << "time from previous BKV " << time_span.count() << std::endl;
			time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - total_start_time);
			out_sstream << "time from start " << time_span.count() << std::endl << std::endl; 
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
			
			ofile.open( "out", std::ios_base::app );
			ofile << out_sstream.str();
			ofile.close();
			out_sstream.clear(); out_sstream.str("");
		}
		// check if it's time to stop program
		if (number_of_comb - solved_tasks_count < (unsigned)(corecount - 1)) {
			// there are some idle processes right now
			// check of every task was processed with enough time
			bool isTimeToInturrupt = true;
			unsigned final_tasks_count = 0;
			for (auto &x : total_fragment_data) {
				if (x.end_processing_time == 0.0) {
					final_tasks_count++;
					if (MPI_Wtime() - x.start_processing_time < WAIT_FINAL_PROCESS_SECONDS)
						isTimeToInturrupt = false;
				}
			}
			// if at least 1 process works less than limit then don't interrupt
			if (isTimeToInturrupt) {
				out_sstream << "*** stop criteria ***" << std::endl;
				out_sstream << "final_tasks_count " << final_tasks_count << std::endl;
				
				ofile.open("out", std::ios_base::app);
				ofile << out_sstream.str();
				ofile.close();
				out_sstream.clear(); out_sstream.str("");
			}
		}
	}
}

void computingProcess(int rank, int corecount, odls_sequential odls_seq)
{
	std::vector<odls_pair> odls_pair_vec;
	
	odls_seq.readOdlsPairs(odls_pair_vec);
	
	// check pseudotriples based on known DLS from pairs
	unsigned preprocess_bkv = 0;
	std::set<dls> dls_known;
	for ( auto &x : odls_pair_vec ) {
		dls_known.insert( x.dls_1 );
		dls_known.insert( x.dls_2 );
	}
	
	// on prerpocess make pseudotriples based on DLS from pairs 
	dls tmp_dls;
	odls_pseudotriple psudotriple;
	for (auto &x : odls_pair_vec) {
		for (auto &y : dls_known) {
			tmp_dls = y;
			if ((tmp_dls != x.dls_1) && (tmp_dls != x.dls_2)) {
				odls_seq.makePseudotriple(x, tmp_dls, psudotriple);
				if (psudotriple.unique_orthogonal_cells.size() >= preprocess_bkv) {
					preprocess_bkv = psudotriple.unique_orthogonal_cells.size();
					std::cout << "preprocess_bkv " << preprocess_bkv << std::endl;
				}
				// check Brown psuedotriple (2 pairs are orthogonal)
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
	std::cout << "preprocess_bkv based on known DLS from input file : " << preprocess_bkv << std::endl;
	
	MPI_Status status;
	int old_fragment_index;
	int result;
	int value_from_control_process;
	int bkv_from_control_process;
	int fragment_index = 0;
	bool isMessageSent;

	// repeat solving tasks from control process
	for (;;) {
		old_fragment_index = fragment_index;
		// firstly check messages with current BKV
		isMessageSent = false;
		do {
			MPI_Recv(&value_from_control_process, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
			if (value_from_control_process < 0) {
				bkv_from_control_process = abs(value_from_control_process);
				if ((unsigned)bkv_from_control_process > odls_seq.best_all_dls_psudotriple.unique_orthogonal_cells.size() + MAX_DIFF_VALUE_FROM_BKV) {
					result = STOP_DUE_LOW_LOCAL_BKV;
					if (!isMessageSent) {
						// result > 0, i.e. task solving was interrupted, ask for new task
						MPI_Send(&fragment_index, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
						MPI_Send(&result, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
						isMessageSent = true;
					}
				}
			}
		} while (value_from_control_process < 0);
		
		std::cout << "received fragment_index " << fragment_index << std::endl;
		fragment_index = value_from_control_process;
		if (fragment_index < old_fragment_index) {
			std::cerr << "fragment_index < old_fragment_index" << std::endl;
			std::cerr << fragment_index << " < " << old_fragment_index << std::endl;
			std::cerr << "may be insead of index we got new bkv" << std::endl;
			exit(1);
		}
		
		// TODO for Alexey: launch deterministic_generate_dls for short time - for checking conditions
		result = odls_seq.deterministicGeneratingDLS(odls_pair_vec, fragment_index);
		
		odls_seq.best_all_dls_psudotriple.unique_orthogonal_cells.clear();
		odls_seq.best_one_dls_psudotriple.unique_orthogonal_cells.clear();
		
		if (result == 0) { // whole subspace of the search space was processed
			std::cout << "fragment_index " << fragment_index << " processed all subspace" << std::endl;
			break;
		}
		
		std::cout << "deterministic_generate_dls() interrupted on rank " << rank << std::endl;
		
		if (result == STOP_DUE_NO_DLS)
			std::cout << "STOP_DUE_NO_DLS on rank " << rank << std::endl;
		else if (result == STOP_DUE_LOW_LOCAL_BKV)
			std::cout << "STOP_DUE_LOW_LOCAL_BKV on rank " << rank << std::endl;
		else {
			std::cerr << "incorrect result value " << result << std::endl;
			exit(1);
		}
		
		// result > 0, i.e. task solving was interrupted, ask for new task
		MPI_Send(&fragment_index, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		MPI_Send(&result, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}
}

// write to file current state for every fragment
void updateFragmentFile(std::vector<fragment_data> total_fragment_data, double mpi_start_time )
{
	std::ofstream fragment_file("fragment_file", std::ios_base::out);
	fragment_file << "total_data_from_fragment " << MPI_Wtime() - mpi_start_time << " seconds from start " << std::endl;
	fragment_file << "fragment_index first_dls_generate_time genereated_DLS_count orthogonal_value result time" << std::endl;
	unsigned k = 0;
	for (auto &x : total_fragment_data) {
		fragment_file << k++ << " " << x.first_dls_generate_time << " s " << x.generated_DLS_count << " " << x.orthogonal_value << " " << x.result << " ";
		if (x.result != 0)// if fragment was finished
			fragment_file << x.end_processing_time - x.start_processing_time << " s ";
		else
			fragment_file << "in progress";
		fragment_file << std::endl;
	}
	fragment_file.close();
	fragment_file.clear();
}