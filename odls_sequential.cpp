#ifdef _MPI
#include <mpi.h>
#endif

#include "odls_sequential.h"

odls_sequential::odls_sequential() :
	dls_generate_start_time (0),
	dls_generate_last_time(0),
	generated_DLS_count(0),
	dls_total_time(0),
	pseudotriples_total_time(0),
	first_dls_generate_time(0)
{
}

void odls_sequential::readOdlsPairs(std::string known_podls_file_name, std::vector<odls_pair> &odls_pair_vec)
{
	std::ifstream odls_pairs_file(known_podls_file_name.c_str());
	if (!odls_pairs_file.is_open()) {
		std::cerr << known_podls_file_name << " is not open" << std::endl;
		exit(1);
	}
	std::string str, nonspace_str;
	odls_pair_vec.resize(0);
	odls_pair cur_odls_pair;
	while (getline(odls_pairs_file, str)) {
		if ((str.size() <= 1) && (cur_odls_pair.dls_1.size() > 0)) { // if separate string
			odls_pair_vec.push_back(cur_odls_pair); // add current odls pair
			cur_odls_pair.dls_1.clear();
			cur_odls_pair.dls_2.clear();
		}
		else {
			// read two rows of 2 DLS for current pair
			nonspace_str = "";
			for (auto &x : str)
				if ((x != ' ') && (x != '\r'))
					nonspace_str += x;
			if (nonspace_str.size() != 2 * LS_ORDER) {
				std::cerr << "nonspace_str.size() != 2*LS_ORDER" << std::endl;
				std::cerr << nonspace_str.size() << " != " << 2 * LS_ORDER << std::endl;
				exit(1);
			}
			cur_odls_pair.dls_1.push_back(nonspace_str.substr(0, LS_ORDER));
			cur_odls_pair.dls_2.push_back(nonspace_str.substr(LS_ORDER, LS_ORDER));
		}
	}
	odls_pairs_file.close();
	// if there is no separate string at the end of file, add last pair manually
	if (cur_odls_pair.dls_1.size() != 0)
		odls_pair_vec.push_back(cur_odls_pair); // add current odls pair

	/*for ( auto &x : odls_pair_vec ) {
	for ( auto &y : x.dls_1 )
	std::cout << y << std::endl;
	std::cout << std::endl;
	for ( auto &z : x.dls_2 )
	std::cout << z << std::endl;
	std::cout << std::endl << std::endl;
	}*/

	std::set<std::string> greece_latin_square;
	std::string cell_plus_cell;
	// check every pair
	for (auto &x : odls_pair_vec) {
		for (unsigned j1 = 0; j1 < x.dls_1.size(); j1++)
			for (unsigned j2 = 0; j2 < x.dls_1[j1].size(); j2++) {
				cell_plus_cell = x.dls_1[j1][j2];
				cell_plus_cell += x.dls_2[j1][j2];
				greece_latin_square.insert(cell_plus_cell);
			}
		//std::cout << "greece_latin_square.size() " << greece_latin_square.size() << std::endl;
		if (greece_latin_square.size() != x.dls_1.size() * x.dls_1.size()) {
			std::cerr << "greece_latin_square.size() != x.dls_1.size() * x.dls_1.size() " << std::endl;
			std::cerr << greece_latin_square.size() << " != " << x.dls_1.size() * x.dls_1.size() << std::endl;
		}
	}
}

void odls_sequential::makePseudotriple(odls_pair &orthogonal_pair, dls &new_dls, odls_pseudotriple &pseudotriple)
{
	unsigned cur_first_pair_orthogonal_cells, cur_second_pair_orthogonal_cells;
	std::set<std::string> greece_latin_square1, greece_latin_square2; // for counting orthogonal cells between 2 DLS
	std::string cell_plus_cell;

	for (unsigned j1 = 0; j1< new_dls.size(); j1++)
		for (unsigned j2 = 0; j2 < new_dls[j1].size(); j2++) {
			cell_plus_cell = orthogonal_pair.dls_1[j1][j2];
			cell_plus_cell += new_dls[j1][j2];
			greece_latin_square1.insert(cell_plus_cell);
		}
	cur_first_pair_orthogonal_cells = greece_latin_square1.size();
	/*std::cout << "greece_latin_square1 " << std::endl;
	for ( auto &y : greece_latin_square1 )
	std::cout << y << " ";
	std::cout << std::endl;*/
	for (unsigned j1 = 0; j1 < new_dls.size(); j1++)
		for (unsigned j2 = 0; j2 < new_dls[j1].size(); j2++) {
			cell_plus_cell = orthogonal_pair.dls_2[j1][j2];
			cell_plus_cell += new_dls[j1][j2];
			greece_latin_square2.insert(cell_plus_cell);
		}
	cur_second_pair_orthogonal_cells = greece_latin_square2.size();
	pseudotriple.dls_1 = orthogonal_pair.dls_1;
	pseudotriple.dls_2 = orthogonal_pair.dls_2;
	pseudotriple.dls_3 = new_dls;
	pseudotriple.unique_orthogonal_cells.clear();
	std::set_intersection(greece_latin_square1.begin(), greece_latin_square1.end(),
		greece_latin_square2.begin(), greece_latin_square2.end(),
		std::inserter(pseudotriple.unique_orthogonal_cells, pseudotriple.unique_orthogonal_cells.begin()));
}

void odls_sequential::processNewDLS(std::vector<odls_pair> &odls_pair_vec, int fragment_index, unsigned short int *square)
{
#ifdef _MPI
	std::string cur_string_set;
	unsigned k;
	std::stringstream sstream;
	unsigned best_first_pair_orthogonal_cells = 0, best_second_pair_orthogonal_cells = 0;
	dls new_dls;
	new_dls.resize(LS_ORDER);
	char psuedotriple_char_arr[psuedotriple_char_arr_len];
	unsigned char_index;
	unsigned ortogonal_cells;

	dls_total_time += MPI_Wtime() - dls_generate_last_time;

	// time for generating first DLS
	if (generated_DLS_count == 0)
		first_dls_generate_time = MPI_Wtime() - dls_generate_last_time;

	dls_generate_last_time = MPI_Wtime();

	// here we have diagonal Latin square, let's add it to the set
	generated_DLS_count++;

	if (generated_DLS_count == 1)
		first_dls_generate_time = MPI_Wtime() - dls_generate_start_time;
	
	k = 0;
	cur_string_set = "";

	double prev_time = MPI_Wtime();

	for (int j1 = 0; j1 < LS_ORDER; j1++) {
		for (int j2 = 0; j2 < LS_ORDER; j2++)
			sstream << square[j1*LS_ORDER + j2];
		new_dls[j1] = sstream.str();
		cur_string_set += sstream.str();
		sstream.clear(); sstream.str("");
	}

	// make pseudotriple for every pair and choose one with maximum of orthogonal cells
	for (auto &x : odls_pair_vec) {
		makePseudotriple(x, new_dls, dls_psudotriple);
		if (dls_psudotriple.unique_orthogonal_cells.size() > best_one_dls_psudotriple.unique_orthogonal_cells.size())
			best_one_dls_psudotriple = dls_psudotriple;
		//std::cout << "cur_first_pair_orthogonal_cells " << cur_first_pair_orthogonal_cells << std::endl;
		//std::cout << "cur_second_pair_orthogonal_cells " << cur_second_pair_orthogonal_cells << std::endl;
		//std::cout << "cur_one_dls_psudotriple_orthogonal_cells " << cur_one_dls_psudotriple_orthogonal_cells << std::endl;
		//std::cout << "best_one_dls_psudotriple_orthogonal_cells " << best_one_dls_psudotriple_orthogonal_cells << std::endl;
	}

	if (best_one_dls_psudotriple.unique_orthogonal_cells.size() > best_all_dls_psudotriple.unique_orthogonal_cells.size()) {
		best_all_dls_psudotriple = best_one_dls_psudotriple;
		std::cout << "***" << std::endl;
		std::cout << "best_all_dls_psudotriple_orthogonal_cells " << best_all_dls_psudotriple.unique_orthogonal_cells.size() << std::endl;

		std::cout << "time from start " << MPI_Wtime() - dls_generate_start_time << std::endl;

		std::cout << "genereated_DLS_count " << generated_DLS_count << std::endl;
		std::cout << "dls_total_time " << dls_total_time << std::endl;
		std::cout << "pseudotriples_total_time " << pseudotriples_total_time << std::endl;
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

		MPI_Send(&fragment_index, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		MPI_Send(&ortogonal_cells, 1, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD);
		MPI_Send(psuedotriple_char_arr, psuedotriple_char_arr_len, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
		MPI_Send(&first_dls_generate_time, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		MPI_Send(&generated_DLS_count, 1, MPI_UNSIGNED_LONG_LONG, 0, 0, MPI_COMM_WORLD);
	}
	
	pseudotriples_total_time += MPI_Wtime() - prev_time;
#endif
}

// ДЛК находится в square[100]( а именнно square[0]-square[99]) в момент времени отмеченный коментарием на 774 строчке
int odls_sequential::deterministicGeneratingDLS(std::vector<odls_pair> &odls_pair_vec, unsigned fragment_index)
{
	return 0;
}

void compareLocalRecordWithGlobal(int count, int local_max)
{
	//получить глобальный рекорд
	//int global_max=твоя_функция_получения_глобального_максимума_ортогональности_среди_всех_областей();
		
	/*if(local_max>global_max)	
	{
		global_max=local_max;
		передача в глобальную область памяти нового рекорда
	}
	else if((local_max<=global_max)&&(local_max>=global_max-DIFFER_GLOBAL_LOCAL))
	{
		count=0;
	}
	else
	{
		return(local_max);
	}*/
}
