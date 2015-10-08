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

void odls_sequential::readOdlsPairs(std::vector<odls_pair> &odls_pair_vec)
{
	std::string odls_pairs_file_name = "ODLS_10_pairs.txt";
	std::ifstream odls_pairs_file(odls_pairs_file_name);
	if (!odls_pairs_file.is_open()) {
		std::cerr << odls_pairs_file_name << " is not open" << std::endl;
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
			if (nonspace_str.size() != 2 * LS_order) {
				std::cerr << "nonspace_str.size() != 2*LS_order" << std::endl;
				std::cerr << nonspace_str.size() << " != " << 2 * LS_order << std::endl;
				exit(1);
			}
			cur_odls_pair.dls_1.push_back(nonspace_str.substr(0, LS_order));
			cur_odls_pair.dls_2.push_back(nonspace_str.substr(LS_order, LS_order));
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
	new_dls.resize(LS_order);
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

	prev_time = MPI_Wtime();

	for (int j1 = 0; j1 < LS_order; j1++) {
		for (int j2 = 0; j2 < LS_order; j2++)
			sstream << square[j1*LS_order + j2];
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
#ifdef _MPI
	unsigned short int square[100] = { 0 };
	unsigned short int flag[100] = { 0 };
	unsigned short int start[15] = { 0 };
	unsigned short int end[15] = { 9 };

	int number_min = 0;
	int number_max = 0;
	int number_min_real = 0;
	int number_max_real = 0;

	long long int count = 0;
	int i = 0;
	int j = 0;
	int number_of_comb_in_one_part = 1;
	
#ifdef _MPI
	MPI_Status mpi_status;
	MPI_Request mpi_request;
#endif
	
	double elapsed_time;

	int iprobe_message = 0, message_size = 0;
	int bkv_from_control_process = 0;

	dls_generate_start_time = MPI_Wtime();
	dls_generate_last_time = MPI_Wtime();

	std::cout << "Start of generating DLS" << std::endl;

	generated_DLS_count = 0;
	dls_total_time = 0; 
	pseudotriples_total_time = 0;
	first_dls_generate_time = 0;
	best_all_dls_psudotriple.unique_orthogonal_cells.clear();
	
	/*float fparts = parts;
	fparts = number_of_comb / fparts;
	double fparts = parts;*/

	//округление в +//
	/*int inter;
	float finter;
	inter = fparts;
	finter = inter;
	if (fparts != finter)
	fparts = finter + 1;
	else
	fparts = finter;*/

	// генерация диапазона поиска
	for (square[0] = 0; square[0] <= 0; square[0]++) /*инициализцая 0 элемента*/
	{
		for (square[1] = 1; square[1] <= 1; square[1]++) /*инициализцая 1 элемента*/
		{
			for (square[2] = 2; ((square[2] <= 2) && (square[1] != square[0])); square[2]++)  /*инициализцая 2 элемента*/
			{
				if ((square[2] != square[1]) && (square[2] != square[0]))
				{
					flag[2] = 1;
				}
				for (square[3] = 3; ((square[3] <= 3) && (flag[2] == 1)); square[3]++)  /*инициализцая 3 элемента*/
				{
					if ((square[3] != square[2]) && (square[3] != square[1]) && (square[3] != square[0]))
					{
						flag[3] = 1;
					}
					for (square[4] = 4; ((square[4] <= 4) && (flag[3] == 1)); square[4]++) /*инициализцая 4 элемента*/
					{
						if ((square[4] != square[3]) && (square[4] != square[2]) && (square[4] != square[1]) && (square[4] != square[0]))
						{
							flag[4] = 1;
						}
						for (square[5] = 5; ((square[5] <= 5) && (flag[4] == 1)); square[5]++) /*инициализцая 5 элемента*/
						{
							if ((square[5] != square[4]) && (square[5] != square[3]) && (square[5] != square[2]) && (square[5] != square[1]) && (square[5] != square[0]))
							{
								flag[5] = 1;
							}
							for (square[6] = 6; ((square[6] <= 6) && (flag[5] == 1)); square[6]++) /*инициализцая 6 элемента*/
							{
								if ((square[6] != square[5]) && (square[6] != square[4]) && (square[6] != square[3]) && (square[6] != square[2]) && (square[6] != square[1]) && (square[6] != square[0]))
								{
									flag[6] = 1;
								}
								for (square[7] = 7; ((square[7] <= 7) && (flag[6] == 1)); square[7]++) /*инициализцая 7 элемента*/
								{
									if ((square[7] != square[6]) && (square[7] != square[5]) && (square[7] != square[4]) && (square[7] != square[3]) && (square[7] != square[2]) && (square[7] != square[1]) && (square[7] != square[0]))
									{
										flag[7] = 1;
									}
									for (square[8] = 8; ((square[8] <= 8) && (flag[7] == 1)); square[8]++) /*инициализцая 8 элемента*/
									{
										if ((square[8] != square[7]) && (square[8] != square[6]) && (square[8] != square[5]) && (square[8] != square[4]) && (square[8] != square[3]) && (square[8] != square[2]) && (square[8] != square[1]) && (square[8] != square[0]))
										{
											flag[8] = 1;
										}
										for (square[9] = 9; ((square[9] <= 9) && (flag[8] == 1)); square[9]++) /*инициализцая 9 элемента*/
										{
											if ((square[9] != square[8]) && (square[9] != square[7]) && (square[9] != square[6]) && (square[9] != square[5]) && (square[9] != square[4]) && (square[9] != square[3]) && (square[9] != square[2]) && (square[9] != square[1]) && (square[9] != square[0]))
											{
												flag[9] = 1;
											}
											for (square[10] = 0; ((square[10] < 10) && (flag[9] == 1)); square[10]++) /*инициализцая 10 элемента*/
											{
												for (square[11] = 0; ((square[11] < 10) && (square[10] != square[0])); square[11]++) /*инициализцая 11 элемента*/
												{
													if ((square[11] != square[10]) && (square[11] != square[1]) && (square[11] != square[0]))
													{
														flag[11] = 1;
													}
													for (square[12] = 0; ((square[12] < 10) && (flag[11] == 1)); square[12]++) /*инициализцая 12 элемента*/
													{
														if ((square[12] != square[11]) && (square[12] != square[10]) && (square[12] != square[2]))
														{
															flag[12] = 1;
														}
														for (square[13] = 0; ((square[13] < 10) && (flag[12] == 1)); square[13]++) /*инициализцая 13 элемента*/
														{
															if ((square[13] != square[12]) && (square[13] != square[11]) && (square[13] != square[10]) && (square[13] != square[3]))
															{
																flag[13] = 1;
															}
															for (square[14] = 0; ((square[14] < 10) && (flag[13] == 1)); square[14]++) /*инициализцая 14 элемента*/
															{
																if ((square[14] != square[13]) && (square[14] != square[12]) && (square[14] != square[11]) && (square[14] != square[10]) && (square[14] != square[4]))
																{
																	flag[14] = 1;

																	if (count == (number_of_comb_in_one_part*fragment_index))
																	{
																		for (j = 0; j<15; j++)
																		{
																			start[j] = square[j];
																		}
																	}

																	if ((count >= number_of_comb_in_one_part*fragment_index) && (count <= (number_of_comb_in_one_part*(fragment_index + 1) - 1)))
																	{
																		for (j = 0; j<15; j++)
																		{
																			end[j] = square[j];
																		}
																	}
																	count++;

																}
																flag[14] = 0;
															}
															flag[13] = 0;
														}
														flag[12] = 0;
													}
													flag[11] = 0;
												}

											}
											flag[9] = 0;
										}
										flag[8] = 0;
									}
									flag[7] = 0;
								}
								flag[6] = 0;
							}
							flag[5] = 0;
						}
						flag[4] = 0;
					}
					flag[3] = 0;
				}
				flag[2] = 0;
			}
		}
	}

	// Подсчет пороговых значений порогов
	number_min = (start[10] * 10000) + (start[11] * 1000) + (start[12] * 100) + (start[13] * 10) + start[14];
	number_max = (end[10] * 10000) + (end[11] * 1000) + (end[12] * 100) + (end[13] * 10) + end[14];

	//обнуление квадратов и флагов перед генерацией ДЛК
	for (i = 0; i<100; i++)
	{
		square[i] = 0;
		flag[i] = 0;
	}
	count = 0;

	//генерация ДЛК

	for (square[0] = 0; square[0] == 0; square[0]++) /*инициализцая 0 элемента*/
	{
		for (square[1] = 1; square[1] == 1; square[1]++) /*инициализцая 1 элемента*/
		{
			for (square[2] = 2; ((square[2] == 2) && (square[1] != square[0])); square[2]++)  /*инициализцая 2 элемента*/
			{
				if ((square[2] != square[1]) && (square[2] != square[0]))
				{
					flag[2] = 1;
				}
				for (square[3] = 3; ((square[3] == 3) && (flag[2] == 1)); square[3]++)  /*инициализцая 3 элемента*/
				{
					if ((square[3] != square[2]) && (square[3] != square[1]) && (square[3] != square[0]))
					{
						flag[3] = 1;
					}
					for (square[4] = 4; ((square[4] == 4) && (flag[3] == 1)); square[4]++) /*инициализцая 4 элемента*/
					{
						if ((square[4] != square[3]) && (square[4] != square[2]) && (square[4] != square[1]) && (square[4] != square[0]))
						{
							flag[4] = 1;
						}
						for (square[5] = 5; ((square[5] == 5) && (flag[4] == 1)); square[5]++) /*инициализцая 5 элемента*/
						{
							if ((square[5] != square[4]) && (square[5] != square[3]) && (square[5] != square[2]) && (square[5] != square[1]) && (square[5] != square[0]))
							{
								flag[5] = 1;
							}
							for (square[6] = 6; ((square[6] == 6) && (flag[5] == 1)); square[6]++) /*инициализцая 6 элемента*/
							{
								if ((square[6] != square[5]) && (square[6] != square[4]) && (square[6] != square[3]) && (square[6] != square[2]) && (square[6] != square[1]) && (square[6] != square[0]))
								{
									flag[6] = 1;
								}
								for (square[7] = 7; ((square[7] == 7) && (flag[6] == 1)); square[7]++) /*инициализцая 7 элемента*/
								{
									if ((square[7] != square[6]) && (square[7] != square[5]) && (square[7] != square[4]) && (square[7] != square[3]) && (square[7] != square[2]) && (square[7] != square[1]) && (square[7] != square[0]))
									{
										flag[7] = 1;
									}
									for (square[8] = 8; ((square[8] == 8) && (flag[7] == 1)); square[8]++) /*инициализцая 8 элемента*/
									{
										if ((square[8] != square[7]) && (square[8] != square[6]) && (square[8] != square[5]) && (square[8] != square[4]) && (square[8] != square[3]) && (square[8] != square[2]) && (square[8] != square[1]) && (square[8] != square[0]))
										{
											flag[8] = 1;
										}
										for (square[9] = 9; ((square[9] == 9) && (flag[8] == 1)); square[9]++) /*инициализцая 9 элемента*/
										{
											if ((square[9] != square[8]) && (square[9] != square[7]) && (square[9] != square[6]) && (square[9] != square[5]) && (square[9] != square[4]) && (square[9] != square[3]) && (square[9] != square[2]) && (square[9] != square[1]) && (square[9] != square[0]))
											{
												flag[9] = 1;
											}
											for (square[10] = 0; ((square[10] < 10) && (flag[9] == 1)); square[10]++) /*инициализцая 10 элемента*/
											{
												for (square[11] = 0; ((square[11] < 10) && (square[10] != square[0])); square[11]++) /*инициализцая 11 элемента*/
												{
													if ((square[11] != square[10]) && (square[11] != square[1]) && (square[11] != square[0]))
													{
														flag[11] = 1;
													}
													for (square[12] = 0; ((square[12] < 10) && (flag[11] == 1)); square[12]++) /*инициализцая 12 элемента*/
													{
														if ((square[12] != square[11]) && (square[12] != square[10]) && (square[12] != square[2]))
														{
															flag[12] = 1;
														}
														for (square[13] = 0; ((square[13] < 10) && (flag[12] == 1)); square[13]++) /*инициализцая 13 элемента*/
														{
															if ((square[13] != square[12]) && (square[13] != square[11]) && (square[13] != square[10]) && (square[13] != square[3]))
															{
																flag[13] = 1;
															}
															for (square[14] = 0; ((square[14] < 10) && (flag[13] == 1)); square[14]++) /*инициализцая 14 элемента*/
															{
																if ((square[14] != square[13]) && (square[14] != square[12]) && (square[14] != square[11]) && (square[14] != square[10]) && (square[14] != square[4]))
																{
																	flag[14] = 1;
																	number_min_real = (square[10] * 10000) + (square[11] * 1000) + (square[12] * 100) + (square[13] * 10) + square[14];
																	number_max_real = (square[10] * 10000) + (square[11] * 1000) + (square[12] * 100) + (square[13] * 10) + square[14];
																}
																for (square[15] = 0; ((square[15] < 10) && (flag[14] == 1) && (number_min_real >= number_min) && (number_max_real <= number_max)); square[15]++) /*инициализцая 15 элемента*/
																{
																	if ((square[15] != square[14]) && (square[15] != square[13]) && (square[15] != square[12]) && (square[15] != square[11]) && (square[15] != square[10]) && (square[15] != square[5]))
																	{
																		flag[15] = 1;
																	}
																	for (square[16] = 0; ((square[16] < 10) && (flag[15] == 1)); square[16]++) /*инициализцая 16 элемента*/
																	{
																		if ((square[16] != square[15]) && (square[16] != square[14]) && (square[16] != square[13]) && (square[16] != square[12]) && (square[16] != square[11]) && (square[16] != square[10]) && (square[16] != square[6]))
																		{
																			flag[16] = 1;
																		}
																		for (square[17] = 0; ((square[17] < 10) && (flag[16] == 1)); square[17]++) /*инициализцая 17 элемента*/
																		{
																			if ((square[17] != square[16]) && (square[17] != square[15]) && (square[17] != square[14]) && (square[17] != square[13]) && (square[17] != square[12]) && (square[17] != square[11]) && (square[17] != square[10]) && (square[17] != square[7]))
																			{
																				flag[17] = 1;
																			}
																			for (square[18] = 0; ((square[18] < 10) && (flag[17] == 1)); square[18]++) /*инициализцая 18 элемента*/
																			{
																				if ((square[18] != square[17]) && (square[18] != square[16]) && (square[18] != square[15]) && (square[18] != square[14]) && (square[18] != square[13]) && (square[18] != square[12]) && (square[18] != square[11]) && (square[18] != square[10]) && (square[18] != square[9]) && (square[18] != square[8]))
																				{
																					flag[18] = 1;
																				}
																				for (square[19] = 0; ((square[19] < 10) && (flag[18] == 1)); square[19]++) /*инициализцая 19 элемента*/
																				{
																					if ((square[19] != square[18]) && (square[19] != square[17]) && (square[19] != square[16]) && (square[19] != square[15]) && (square[19] != square[14]) && (square[19] != square[13]) && (square[19] != square[12]) && (square[19] != square[11]) && (square[19] != square[10]) && (square[19] != square[9]))
																					{
																						flag[19] = 1;
																					}
																					for (square[20] = 0; ((square[20] < 10) && (flag[19] == 1)); square[20]++) /*инициализцая 20 элемента*/
																					{
																						if ((square[20] != square[10]) && (square[20] != square[0]))
																						{
																							flag[20] = 1;
																						}
																						for (square[21] = 0; ((square[21] < 10) && (flag[20] == 1)); square[21]++) /*инициализцая 21 элемента*/
																						{
																							if ((square[21] != square[20]) && (square[21] != square[11]) && (square[21] != square[1]))
																							{
																								flag[21] = 1;
																							}
																							for (square[22] = 0; ((square[22] < 10) && (flag[21] == 1)); square[22]++) /*инициализцая 22 элемента*/
																							{
																								if ((square[22] != square[21]) && (square[22] != square[20]) && (square[22] != square[12]) && (square[22] != square[11]) && (square[22] != square[2]) && (square[22] != square[0]))
																								{
																									flag[22] = 1;
																								}
																								for (square[23] = 0; ((square[23] < 10) && (flag[22] == 1)); square[23]++) /*инициализцая 23 элемента*/
																								{
																									if ((square[23] != square[22]) && (square[23] != square[21]) && (square[23] != square[20]) && (square[23] != square[13]) && (square[23] != square[3]))
																									{
																										flag[23] = 1;
																									}
																									for (square[24] = 0; ((square[24] < 10) && (flag[23] == 1)); square[24]++) /*инициализцая 24 элемента*/
																									{
																										if ((square[24] != square[23]) && (square[24] != square[22]) && (square[24] != square[21]) && (square[24] != square[20]) && (square[24] != square[14]) && (square[24] != square[4]))
																										{
																											flag[24] = 1;
																										}
																										for (square[25] = 0; ((square[25] < 10) && (flag[24] == 1)); square[25]++) /*инициализцая 25 элемента*/
																										{
																											if ((square[25] != square[24]) && (square[25] != square[23]) && (square[25] != square[22]) && (square[25] != square[21]) && (square[25] != square[20]) && (square[25] != square[15]) && (square[25] != square[5]))
																											{
																												flag[25] = 1;
																											}
																											for (square[26] = 0; ((square[26] < 10) && (flag[25] == 1)); square[26]++) /*инициализцая 26 элемента*/
																											{
																												if ((square[26] != square[25]) && (square[26] != square[24]) && (square[26] != square[23]) && (square[26] != square[22]) && (square[26] != square[21]) && (square[26] != square[20]) && (square[26] != square[16]) && (square[26] != square[6]))
																												{
																													flag[26] = 1;
																												}
																												for (square[27] = 0; ((square[27] < 10) && (flag[26] == 1)); square[27]++) /*инициализцая 27 элемента*/
																												{
																													if ((square[27] != square[26]) && (square[27] != square[25]) && (square[27] != square[24]) && (square[27] != square[23]) && (square[27] != square[22]) && (square[27] != square[21]) && (square[27] != square[20]) && (square[27] != square[18]) && (square[27] != square[17]) && (square[27] != square[9]) && (square[27] != square[7]))
																													{
																														flag[27] = 1;
																													}
																													for (square[28] = 0; ((square[28] < 10) && (flag[27] == 1)); square[28]++) /*инициализцая 28 элемента*/
																													{
																														if ((square[28] != square[27]) && (square[28] != square[26]) && (square[28] != square[25]) && (square[28] != square[24]) && (square[28] != square[23]) && (square[28] != square[22]) && (square[28] != square[21]) && (square[28] != square[20]) && (square[28] != square[18]) && (square[28] != square[8]))
																														{
																															flag[28] = 1;
																														}
																														for (square[29] = 0; ((square[29] < 10) && (flag[28] == 1)); square[29]++) /*инициализцая 29 элемента*/
																														{
																															if ((square[29] != square[28]) && (square[29] != square[27]) && (square[29] != square[26]) && (square[29] != square[25]) && (square[29] != square[24]) && (square[29] != square[23]) && (square[29] != square[22]) && (square[29] != square[21]) && (square[29] != square[20]) && (square[29] != square[19]) && (square[29] != square[9]))
																															{
																																flag[29] = 1;
																															}
																															for (square[30] = 0; ((square[30] < 10) && (flag[29] == 1)); square[30]++) /*инициализцая 30 элемента*/
																															{
																																if ((square[30] != square[20]) && (square[30] != square[10]) && (square[30] != square[0]))
																																{
																																	flag[30] = 1;
																																}
																																for (square[31] = 0; ((square[31] < 10) && (flag[30] == 1)); square[31]++) /*инициализцая 31 элемента*/
																																{
																																	if ((square[31] != square[30]) && (square[31] != square[21]) && (square[31] != square[11]) && (square[31] != square[1]))
																																	{
																																		flag[31] = 1;
																																	}
																																	for (square[32] = 0; ((square[32] < 10) && (flag[31] == 1)); square[32]++) /*инициализцая 32 элемента*/
																																	{
																																		if ((square[32] != square[31]) && (square[32] != square[30]) && (square[32] != square[22]) && (square[32] != square[12]) && (square[32] != square[2]))
																																		{
																																			flag[32] = 1;
																																		}
																																		for (square[33] = 0; ((square[33] < 10) && (flag[32] == 1)); square[33]++) /*инициализцая 33 элемента*/
																																		{
																																			if ((square[33] != square[32]) && (square[33] != square[31]) && (square[33] != square[30]) && (square[33] != square[23]) && (square[33] != square[22]) && (square[33] != square[13]) && (square[33] != square[11]) && (square[33] != square[3]) && (square[33] != square[0]))
																																			{
																																				flag[33] = 1;
																																			}
																																			for (square[34] = 0; ((square[34] < 10) && (flag[33] == 1)); square[34]++) /*инициализцая 34 элемента*/
																																			{
																																				if ((square[34] != square[33]) && (square[34] != square[32]) && (square[34] != square[31]) && (square[34] != square[30]) && (square[34] != square[24]) && (square[34] != square[14]) && (square[34] != square[4]))
																																				{
																																					flag[34] = 1;
																																				}
																																				for (square[35] = 0; ((square[35] < 10) && (flag[34] == 1)); square[35]++) /*инициализцая 35 элемента*/
																																				{
																																					if ((square[35] != square[34]) && (square[35] != square[33]) && (square[35] != square[32]) && (square[35] != square[31]) && (square[35] != square[30]) && (square[35] != square[25]) && (square[35] != square[15]) && (square[35] != square[5]))
																																					{
																																						flag[35] = 1;
																																					}
																																					for (square[36] = 0; ((square[36] < 10) && (flag[35] == 1)); square[36]++) /*инициализцая 36 элемента*/
																																					{
																																						if ((square[36] != square[35]) && (square[36] != square[34]) && (square[36] != square[33]) && (square[36] != square[32]) && (square[36] != square[31]) && (square[36] != square[30]) && (square[36] != square[27]) && (square[36] != square[26]) && (square[36] != square[18]) && (square[36] != square[16]) && (square[36] != square[9]) && (square[36] != square[6]))
																																						{
																																							flag[36] = 1;
																																						}
																																						for (square[37] = 0; ((square[37] < 10) && (flag[36] == 1)); square[37]++) /*инициализцая 37 элемента*/
																																						{
																																							if ((square[37] != square[36]) && (square[37] != square[35]) && (square[37] != square[34]) && (square[37] != square[33]) && (square[37] != square[32]) && (square[37] != square[31]) && (square[37] != square[30]) && (square[37] != square[27]) && (square[37] != square[17]) && (square[37] != square[7]))
																																							{
																																								flag[37] = 1;
																																							}
																																							for (square[38] = 0; ((square[38] < 10) && (flag[37] == 1)); square[38]++) /*инициализцая 38 элемента*/
																																							{
																																								if ((square[38] != square[37]) && (square[38] != square[36]) && (square[38] != square[35]) && (square[38] != square[34]) && (square[38] != square[33]) && (square[38] != square[32]) && (square[38] != square[31]) && (square[38] != square[30]) && (square[38] != square[28]) && (square[38] != square[18]) && (square[38] != square[8]))
																																								{
																																									flag[38] = 1;
																																								}
																																								for (square[39] = 0; ((square[39] < 10) && (flag[38] == 1)); square[39]++) /*инициализцая 39 элемента*/
																																								{
																																									if ((square[39] != square[38]) && (square[39] != square[37]) && (square[39] != square[36]) && (square[39] != square[35]) && (square[39] != square[34]) && (square[39] != square[33]) && (square[39] != square[32]) && (square[39] != square[31]) && (square[39] != square[30]) && (square[39] != square[29]) && (square[39] != square[19]) && (square[39] != square[9]))
																																									{
																																										flag[39] = 1;
																																									}
																																									for (square[40] = 0; ((square[40] < 10) && (flag[39] == 1)); square[40]++) /*инициализцая 40 элемента*/
																																									{
																																										if ((square[40] != square[30]) && (square[40] != square[20]) && (square[40] != square[10]) && (square[40] != square[0]))
																																										{
																																											flag[40] = 1;
																																										}
																																										for (square[41] = 0; ((square[41] < 10) && (flag[40] == 1)); square[41]++) /*инициализцая 41 элемента*/
																																										{
																																											if ((square[41] != square[40]) && (square[41] != square[31]) && (square[41] != square[21]) && (square[41] != square[11]) && (square[41] != square[1]))
																																											{
																																												flag[41] = 1;
																																											}
																																											for (square[42] = 0; ((square[42] < 10) && (flag[41] == 1)); square[42]++) /*инициализцая 42 элемента*/
																																											{
																																												if ((square[42] != square[41]) && (square[42] != square[40]) && (square[42] != square[32]) && (square[42] != square[22]) && (square[42] != square[12]) && (square[42] != square[2]))
																																												{
																																													flag[42] = 1;
																																												}
																																												for (square[43] = 0; ((square[43] < 10) && (flag[42] == 1)); square[43]++) /*инициализцая 43 элемента*/
																																												{
																																													if ((square[43] != square[42]) && (square[43] != square[41]) && (square[43] != square[40]) && (square[43] != square[33]) && (square[43] != square[23]) && (square[43] != square[13]) && (square[43] != square[3]))
																																													{
																																														flag[43] = 1;
																																													}
																																													for (square[44] = 0; ((square[44] < 10) && (flag[43] == 1)); square[44]++) /*инициализцая 44 элемента*/
																																													{
																																														if ((square[44] != square[43]) && (square[44] != square[42]) && (square[44] != square[41]) && (square[44] != square[40]) && (square[44] != square[34]) && (square[44] != square[33]) && (square[44] != square[24]) && (square[44] != square[22]) && (square[44] != square[14]) && (square[44] != square[11]) && (square[44] != square[4]) && (square[44] != square[0]))
																																														{
																																															flag[44] = 1;
																																														}
																																														for (square[45] = 0; ((square[45] < 10) && (flag[44] == 1)); square[45]++) /*инициализцая 45 элемента*/
																																														{
																																															if ((square[45] != square[44]) && (square[45] != square[43]) && (square[45] != square[42]) && (square[45] != square[41]) && (square[45] != square[40]) && (square[45] != square[36]) && (square[45] != square[35]) && (square[45] != square[27]) && (square[45] != square[25]) && (square[45] != square[18]) && (square[45] != square[15]) && (square[45] != square[9]) && (square[45] != square[5]))
																																															{
																																																flag[45] = 1;
																																															}
																																															for (square[46] = 0; ((square[46] < 10) && (flag[45] == 1)); square[46]++) /*инициализцая 46 элемента*/
																																															{
																																																if ((square[46] != square[45]) && (square[46] != square[44]) && (square[46] != square[43]) && (square[46] != square[42]) && (square[46] != square[41]) && (square[46] != square[40]) && (square[46] != square[36]) && (square[46] != square[26]) && (square[46] != square[16]) && (square[46] != square[6]))
																																																{
																																																	flag[46] = 1;
																																																}
																																																for (square[47] = 0; ((square[47] < 10) && (flag[46] == 1)); square[47]++) /*инициализцая 47 элемента*/
																																																{
																																																	if ((square[47] != square[46]) && (square[47] != square[45]) && (square[47] != square[44]) && (square[47] != square[43]) && (square[47] != square[42]) && (square[47] != square[41]) && (square[47] != square[40]) && (square[47] != square[37]) && (square[47] != square[27]) && (square[47] != square[17]) && (square[47] != square[7]))
																																																	{
																																																		flag[47] = 1;
																																																	}
																																																	for (square[48] = 0; ((square[48] < 10) && (flag[47] == 1)); square[48]++) /*инициализцая 48 элемента*/
																																																	{
																																																		if ((square[48] != square[47]) && (square[48] != square[46]) && (square[48] != square[45]) && (square[48] != square[44]) && (square[48] != square[43]) && (square[48] != square[42]) && (square[48] != square[41]) && (square[48] != square[40]) && (square[48] != square[38]) && (square[48] != square[28]) && (square[48] != square[18]) && (square[48] != square[8]))
																																																		{
																																																			flag[48] = 1;
																																																		}
																																																		for (square[49] = 0; ((square[49] < 10) && (flag[48] == 1)); square[49]++) /*инициализцая 49 элемента*/
																																																		{
																																																			if ((square[49] != square[48]) && (square[49] != square[47]) && (square[49] != square[46]) && (square[49] != square[45]) && (square[49] != square[44]) && (square[49] != square[43]) && (square[49] != square[42]) && (square[49] != square[41]) && (square[49] != square[40]) && (square[49] != square[39]) && (square[49] != square[29]) && (square[49] != square[19]) && (square[49] != square[9]))
																																																			{
																																																				flag[49] = 1;
																																																			}
																																																			for (square[50] = 0; ((square[50] < 10) && (flag[49] == 1)); square[50]++) /*инициализцая 50 элемента*/
																																																			{
																																																				if ((square[50] != square[40]) && (square[50] != square[30]) && (square[50] != square[20]) && (square[50] != square[10]) && (square[50] != square[0]))
																																																				{
																																																					flag[50] = 1;
																																																				}
																																																				for (square[51] = 0; ((square[51] < 10) && (flag[50] == 1)); square[51]++) /*инициализцая 51 элемента*/
																																																				{
																																																					if ((square[51] != square[50]) && (square[51] != square[41]) && (square[51] != square[31]) && (square[51] != square[21]) && (square[51] != square[11]) && (square[51] != square[1]))
																																																					{
																																																						flag[51] = 1;
																																																					}
																																																					for (square[52] = 0; ((square[52] < 10) && (flag[51] == 1)); square[52]++) /*инициализцая 52 элемента*/
																																																					{
																																																						if ((square[52] != square[51]) && (square[52] != square[50]) && (square[52] != square[42]) && (square[52] != square[32]) && (square[52] != square[22]) && (square[52] != square[12]) && (square[52] != square[2]))
																																																						{
																																																							flag[52] = 1;
																																																						}
																																																						for (square[53] = 0; ((square[53] < 10) && (flag[52] == 1)); square[53]++) /*инициализцая 53 элемента*/
																																																						{
																																																							if ((square[53] != square[52]) && (square[53] != square[51]) && (square[53] != square[50]) && (square[53] != square[43]) && (square[53] != square[33]) && (square[53] != square[23]) && (square[53] != square[13]) && (square[53] != square[3]))
																																																							{
																																																								flag[53] = 1;
																																																							}
																																																							for (square[54] = 0; ((square[54] < 10) && (flag[53] == 1)); square[54]++) /*инициализцая 54 элемента*/
																																																							{
																																																								if ((square[54] != square[53]) && (square[54] != square[52]) && (square[54] != square[51]) && (square[54] != square[50]) && (square[54] != square[45]) && (square[54] != square[44]) && (square[54] != square[36]) && (square[54] != square[34]) && (square[54] != square[27]) && (square[54] != square[24]) && (square[54] != square[18]) && (square[54] != square[14]) && (square[54] != square[9]) && (square[54] != square[4]))
																																																								{
																																																									flag[54] = 1;
																																																								}
																																																								for (square[55] = 0; ((square[55] < 10) && (flag[54] == 1)); square[55]++) /*инициализцая 55 элемента*/
																																																								{
																																																									if ((square[55] != square[54]) && (square[55] != square[53]) && (square[55] != square[52]) && (square[55] != square[51]) && (square[55] != square[50]) && (square[55] != square[45]) && (square[55] != square[44]) && (square[55] != square[35]) && (square[55] != square[33]) && (square[55] != square[25]) && (square[55] != square[22]) && (square[55] != square[15]) && (square[55] != square[11]) && (square[55] != square[5]) && (square[55] != square[0]))
																																																									{
																																																										flag[55] = 1;
																																																									}
																																																									for (square[56] = 0; ((square[56] < 10) && (flag[55] == 1)); square[56]++) /*инициализцая 56 элемента*/
																																																									{
																																																										if ((square[56] != square[55]) && (square[56] != square[54]) && (square[56] != square[53]) && (square[56] != square[52]) && (square[56] != square[51]) && (square[56] != square[50]) && (square[56] != square[46]) && (square[56] != square[36]) && (square[56] != square[26]) && (square[56] != square[16]) && (square[56] != square[6]))
																																																										{
																																																											flag[56] = 1;
																																																										}
																																																										for (square[57] = 0; ((square[57] < 10) && (flag[56] == 1)); square[57]++) /*инициализцая 57 элемента*/
																																																										{
																																																											if ((square[57] != square[56]) && (square[57] != square[55]) && (square[57] != square[54]) && (square[57] != square[53]) && (square[57] != square[52]) && (square[57] != square[51]) && (square[57] != square[50]) && (square[57] != square[47]) && (square[57] != square[37]) && (square[57] != square[27]) && (square[57] != square[17]) && (square[57] != square[7]))
																																																											{
																																																												flag[57] = 1;
																																																											}
																																																											for (square[58] = 0; ((square[58] < 10) && (flag[57] == 1)); square[58]++) /*инициализцая 58 элемента*/
																																																											{
																																																												if ((square[58] != square[57]) && (square[58] != square[56]) && (square[58] != square[55]) && (square[58] != square[54]) && (square[58] != square[53]) && (square[58] != square[52]) && (square[58] != square[51]) && (square[58] != square[50]) && (square[58] != square[48]) && (square[58] != square[38]) && (square[58] != square[28]) && (square[58] != square[18]) && (square[58] != square[8]))
																																																												{
																																																													flag[58] = 1;
																																																												}
																																																												for (square[59] = 0; ((square[59] < 10) && (flag[58] == 1)); square[59]++) /*инициализцая 59 элемента*/
																																																												{
																																																													if ((square[59] != square[58]) && (square[59] != square[57]) && (square[59] != square[56]) && (square[59] != square[55]) && (square[59] != square[54]) && (square[59] != square[53]) && (square[59] != square[52]) && (square[59] != square[51]) && (square[59] != square[50]) && (square[59] != square[49]) && (square[59] != square[39]) && (square[59] != square[29]) && (square[59] != square[19]) && (square[59] != square[9]))
																																																													{
																																																														flag[59] = 1;
																																																													}
																																																													for (square[60] = 0; ((square[60] < 10) && (flag[59] == 1)); square[60]++) /*инициализцая 60 элемента*/
																																																													{
																																																														if ((square[60] != square[50]) && (square[60] != square[40]) && (square[60] != square[30]) && (square[60] != square[20]) && (square[60] != square[10]) && (square[60] != square[0]))
																																																														{
																																																															flag[60] = 1;
																																																														}
																																																														for (square[61] = 0; ((square[61] < 10) && (flag[60] == 1)); square[61]++) /*инициализцая 61 элемента*/
																																																														{
																																																															if ((square[61] != square[60]) && (square[61] != square[51]) && (square[61] != square[41]) && (square[61] != square[31]) && (square[61] != square[21]) && (square[61] != square[11]) && (square[61] != square[1]))
																																																															{
																																																																flag[61] = 1;
																																																															}
																																																															for (square[62] = 0; ((square[62] < 10) && (flag[61] == 1)); square[62]++) /*инициализцая 62 элемента*/
																																																															{
																																																																if ((square[62] != square[61]) && (square[62] != square[60]) && (square[62] != square[52]) && (square[62] != square[42]) && (square[62] != square[32]) && (square[62] != square[22]) && (square[62] != square[12]) && (square[62] != square[2]))
																																																																{
																																																																	flag[62] = 1;
																																																																}
																																																																for (square[63] = 0; ((square[63] < 10) && (flag[62] == 1)); square[63]++) /*инициализцая 63 элемента*/
																																																																{
																																																																	if ((square[63] != square[62]) && (square[63] != square[61]) && (square[63] != square[60]) && (square[63] != square[54]) && (square[63] != square[53]) && (square[63] != square[45]) && (square[63] != square[43]) && (square[63] != square[36]) && (square[63] != square[33]) && (square[63] != square[27]) && (square[63] != square[23]) && (square[63] != square[18]) && (square[63] != square[13]) && (square[63] != square[9]) && (square[63] != square[3]))
																																																																	{
																																																																		flag[63] = 1;
																																																																	}
																																																																	for (square[64] = 0; ((square[64] < 10) && (flag[63] == 1)); square[64]++) /*инициализцая 64 элемента*/
																																																																	{
																																																																		if ((square[64] != square[63]) && (square[64] != square[62]) && (square[64] != square[61]) && (square[64] != square[60]) && (square[64] != square[54]) && (square[64] != square[44]) && (square[64] != square[34]) && (square[64] != square[24]) && (square[64] != square[14]) && (square[64] != square[4]))
																																																																		{
																																																																			flag[64] = 1;
																																																																		}
																																																																		for (square[65] = 0; ((square[65] < 10) && (flag[64] == 1)); square[65]++) /*инициализцая 65 элемента*/
																																																																		{
																																																																			if ((square[65] != square[64]) && (square[65] != square[63]) && (square[65] != square[62]) && (square[65] != square[61]) && (square[65] != square[60]) && (square[65] != square[55]) && (square[65] != square[45]) && (square[65] != square[35]) && (square[65] != square[25]) && (square[65] != square[15]) && (square[65] != square[5]))
																																																																			{
																																																																				flag[65] = 1;
																																																																			}
																																																																			for (square[66] = 0; ((square[66] < 10) && (flag[65] == 1)); square[66]++) /*инициализцая 66 элемента*/
																																																																			{
																																																																				if ((square[66] != square[65]) && (square[66] != square[64]) && (square[66] != square[63]) && (square[66] != square[62]) && (square[66] != square[61]) && (square[66] != square[60]) && (square[66] != square[56]) && (square[66] != square[55]) && (square[66] != square[46]) && (square[66] != square[44]) && (square[66] != square[36]) && (square[66] != square[33]) && (square[66] != square[26]) && (square[66] != square[22]) && (square[66] != square[16]) && (square[66] != square[11]) && (square[66] != square[6]) && (square[66] != square[0]))
																																																																				{
																																																																					flag[66] = 1;
																																																																				}
																																																																				for (square[67] = 0; ((square[67] < 10) && (flag[66] == 1)); square[67]++) /*инициализцая 67 элемента*/
																																																																				{
																																																																					if ((square[67] != square[66]) && (square[67] != square[65]) && (square[67] != square[64]) && (square[67] != square[63]) && (square[67] != square[62]) && (square[67] != square[61]) && (square[67] != square[60]) && (square[67] != square[57]) && (square[67] != square[47]) && (square[67] != square[37]) && (square[67] != square[27]) && (square[67] != square[17]) && (square[67] != square[7]))
																																																																					{
																																																																						flag[67] = 1;
																																																																					}
																																																																					for (square[68] = 0; ((square[68] < 10) && (flag[67] == 1)); square[68]++) /*инициализцая 68 элемента*/
																																																																					{
																																																																						if ((square[68] != square[67]) && (square[68] != square[66]) && (square[68] != square[65]) && (square[68] != square[64]) && (square[68] != square[63]) && (square[68] != square[62]) && (square[68] != square[61]) && (square[68] != square[60]) && (square[68] != square[58]) && (square[68] != square[48]) && (square[68] != square[38]) && (square[68] != square[28]) && (square[68] != square[18]) && (square[68] != square[8]))
																																																																						{
																																																																							flag[68] = 1;
																																																																						}
																																																																						for (square[69] = 0; ((square[69] < 10) && (flag[68] == 1)); square[69]++) /*инициализцая 69 элемента*/
																																																																						{
																																																																							if ((square[69] != square[68]) && (square[69] != square[67]) && (square[69] != square[66]) && (square[69] != square[65]) && (square[69] != square[64]) && (square[69] != square[63]) && (square[69] != square[62]) && (square[69] != square[61]) && (square[69] != square[60]) && (square[69] != square[59]) && (square[69] != square[49]) && (square[69] != square[39]) && (square[69] != square[29]) && (square[69] != square[19]) && (square[69] != square[9]))
																																																																							{
																																																																								flag[69] = 1;
																																																																							}
																																																																							for (square[70] = 0; ((square[70] < 10) && (flag[69] == 1)); square[70]++) /*инициализцая 70 элемента*/
																																																																							{
																																																																								if ((square[70] != square[60]) && (square[70] != square[50]) && (square[70] != square[40]) && (square[70] != square[30]) && (square[70] != square[20]) && (square[70] != square[10]) && (square[70] != square[0]))
																																																																								{
																																																																									flag[70] = 1;
																																																																								}
																																																																								for (square[71] = 0; ((square[71] < 10) && (flag[70] == 1)); square[71]++) /*инициализцая 71 элемента*/
																																																																								{
																																																																									if ((square[71] != square[70]) && (square[71] != square[61]) && (square[71] != square[51]) && (square[71] != square[41]) && (square[71] != square[31]) && (square[71] != square[21]) && (square[71] != square[11]) && (square[71] != square[1]))
																																																																									{
																																																																										flag[71] = 1;
																																																																									}
																																																																									for (square[72] = 0; ((square[72] < 10) && (flag[71] == 1)); square[72]++) /*инициализцая 72 элемента*/
																																																																									{
																																																																										if ((square[72] != square[71]) && (square[72] != square[70]) && (square[72] != square[63]) && (square[72] != square[62]) && (square[72] != square[54]) && (square[72] != square[52]) && (square[72] != square[45]) && (square[72] != square[42]) && (square[72] != square[36]) && (square[72] != square[32]) && (square[72] != square[27]) && (square[72] != square[22]) && (square[72] != square[18]) && (square[72] != square[12]) && (square[72] != square[9]) && (square[72] != square[2]))
																																																																										{
																																																																											flag[72] = 1;
																																																																										}
																																																																										for (square[73] = 0; ((square[73] < 10) && (flag[72] == 1)); square[73]++) /*инициализцая 73 элемента*/
																																																																										{
																																																																											if ((square[73] != square[72]) && (square[73] != square[71]) && (square[73] != square[70]) && (square[73] != square[63]) && (square[73] != square[53]) && (square[73] != square[43]) && (square[73] != square[33]) && (square[73] != square[23]) && (square[73] != square[13]) && (square[73] != square[3]))
																																																																											{
																																																																												flag[73] = 1;
																																																																											}
																																																																											for (square[74] = 0; ((square[74] < 10) && (flag[73] == 1)); square[74]++) /*инициализцая 74 элемента*/
																																																																											{
																																																																												if ((square[74] != square[73]) && (square[74] != square[72]) && (square[74] != square[71]) && (square[74] != square[70]) && (square[74] != square[64]) && (square[74] != square[54]) && (square[74] != square[44]) && (square[74] != square[34]) && (square[74] != square[24]) && (square[74] != square[14]) && (square[74] != square[4]))
																																																																												{
																																																																													flag[74] = 1;
																																																																												}
																																																																												for (square[75] = 0; ((square[75] < 10) && (flag[74] == 1)); square[75]++) /*инициализцая 75 элемента*/
																																																																												{
																																																																													if ((square[75] != square[74]) && (square[75] != square[73]) && (square[75] != square[72]) && (square[75] != square[71]) && (square[75] != square[70]) && (square[75] != square[65]) && (square[75] != square[55]) && (square[75] != square[45]) && (square[75] != square[35]) && (square[75] != square[25]) && (square[75] != square[15]) && (square[75] != square[5]))
																																																																													{
																																																																														flag[75] = 1;
																																																																													}
																																																																													for (square[76] = 0; ((square[76] < 10) && (flag[75] == 1)); square[76]++) /*инициализцая 76 элемента*/
																																																																													{
																																																																														if ((square[76] != square[75]) && (square[76] != square[74]) && (square[76] != square[73]) && (square[76] != square[72]) && (square[76] != square[71]) && (square[76] != square[70]) && (square[76] != square[66]) && (square[76] != square[56]) && (square[76] != square[46]) && (square[76] != square[36]) && (square[76] != square[26]) && (square[76] != square[16]) && (square[76] != square[6]))
																																																																														{
																																																																															flag[76] = 1;
																																																																														}
																																																																														for (square[77] = 0; ((square[77] < 10) && (flag[76] == 1)); square[77]++) /*инициализцая 77 элемента*/
																																																																														{
																																																																															if ((square[77] != square[76]) && (square[77] != square[75]) && (square[77] != square[74]) && (square[77] != square[73]) && (square[77] != square[72]) && (square[77] != square[71]) && (square[77] != square[70]) && (square[77] != square[67]) && (square[77] != square[66]) && (square[77] != square[57]) && (square[77] != square[55]) && (square[77] != square[47]) && (square[77] != square[44]) && (square[77] != square[37]) && (square[77] != square[33]) && (square[77] != square[27]) && (square[77] != square[22]) && (square[77] != square[17]) && (square[77] != square[11]) && (square[77] != square[7]) && (square[77] != square[0]))
																																																																															{
																																																																																flag[77] = 1;
																																																																															}
																																																																															for (square[78] = 0; ((square[78] < 10) && (flag[77] == 1)); square[78]++) /*инициализцая 78 элемента*/
																																																																															{
																																																																																if ((square[78] != square[77]) && (square[78] != square[76]) && (square[78] != square[75]) && (square[78] != square[74]) && (square[78] != square[73]) && (square[78] != square[72]) && (square[78] != square[71]) && (square[78] != square[70]) && (square[78] != square[68]) && (square[78] != square[58]) && (square[78] != square[48]) && (square[78] != square[38]) && (square[78] != square[28]) && (square[78] != square[18]) && (square[78] != square[8]))
																																																																																{
																																																																																	flag[78] = 1;
																																																																																}
																																																																																for (square[79] = 0; ((square[79] < 10) && (flag[78] == 1)); square[79]++) /*инициализцая 79 элемента*/
																																																																																{
																																																																																	if ((square[79] != square[78]) && (square[79] != square[77]) && (square[79] != square[76]) && (square[79] != square[75]) && (square[79] != square[74]) && (square[79] != square[73]) && (square[79] != square[72]) && (square[79] != square[71]) && (square[79] != square[70]) && (square[79] != square[69]) && (square[79] != square[59]) && (square[79] != square[49]) && (square[79] != square[39]) && (square[79] != square[29]) && (square[79] != square[19]) && (square[79] != square[9]))
																																																																																	{
																																																																																		flag[79] = 1;
																																																																																	}
																																																																																	for (square[80] = 0; ((square[80] < 10) && (flag[79] == 1)); square[80]++) /*инициализцая 80 элемента*/
																																																																																	{
																																																																																		if ((square[80] != square[70]) && (square[80] != square[60]) && (square[80] != square[50]) && (square[80] != square[40]) && (square[80] != square[30]) && (square[80] != square[20]) && (square[80] != square[10]) && (square[80] != square[0]))
																																																																																		{
																																																																																			flag[80] = 1;
																																																																																		}
																																																																																		for (square[81] = 0; ((square[81] < 10) && (flag[80] == 1)); square[81]++) /*инициализцая 81 элемента*/
																																																																																		{
																																																																																			if ((square[81] != square[80]) && (square[81] != square[72]) && (square[81] != square[71]) && (square[81] != square[63]) && (square[81] != square[61]) && (square[81] != square[54]) && (square[81] != square[51]) && (square[81] != square[45]) && (square[81] != square[41]) && (square[81] != square[36]) && (square[81] != square[31]) && (square[81] != square[27]) && (square[81] != square[21]) && (square[81] != square[18]) && (square[81] != square[11]) && (square[81] != square[9]) && (square[81] != square[1]))
																																																																																			{
																																																																																				flag[81] = 1;
																																																																																			}
																																																																																			for (square[82] = 0; ((square[82] < 10) && (flag[81] == 1)); square[82]++) /*инициализцая 82 элемента*/
																																																																																			{
																																																																																				// wait some time, then check for interrupting conditions
																																																																																				elapsed_time = MPI_Wtime() - dls_generate_start_time;
																																																																																				
																																																																																				if (elapsed_time >= WAIT_FIRST_RESULTS_SECONDS) {
																																																																																					if (!generated_DLS_count)
																																																																																						return STOP_DUE_NO_DLS;
																																																																																					else {
																																																																																						iprobe_message = 0;
																																																																																						MPI_Iprobe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &iprobe_message, &mpi_status);
																																																																																						// if message from control process exists
																																																																																						if (iprobe_message) {
																																																																																							MPI_Irecv(&bkv_from_control_process, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &mpi_request);
																																																																																							bkv_from_control_process = abs(bkv_from_control_process);
																																																																																						}
																																																																																						// we could receive bkv_from_control_process later, check it
																																																																																						if ((unsigned)bkv_from_control_process > best_all_dls_psudotriple.unique_orthogonal_cells.size() + MAX_DIFF_VALUE_FROM_BKV)
																																																																																							return STOP_DUE_LOW_LOCAL_BKV;
																																																																																					}
																																																																																				}
																																																																																				if ((square[82] != square[81]) && (square[82] != square[80]) && (square[82] != square[72]) && (square[82] != square[62]) && (square[82] != square[52]) && (square[82] != square[42]) && (square[82] != square[32]) && (square[82] != square[22]) && (square[82] != square[12]) && (square[82] != square[2]))
																																																																																				{
																																																																																					flag[82] = 1;
																																																																																				}
																																																																																				for (square[83] = 0; ((square[83] < 10) && (flag[82] == 1)); square[83]++) /*инициализцая 83 элемента*/
																																																																																				{
																																																																																					if ((square[83] != square[82]) && (square[83] != square[81]) && (square[83] != square[80]) && (square[83] != square[73]) && (square[83] != square[63]) && (square[83] != square[53]) && (square[83] != square[43]) && (square[83] != square[33]) && (square[83] != square[23]) && (square[83] != square[13]) && (square[83] != square[3]))
																																																																																					{
																																																																																						flag[83] = 1;
																																																																																					}
																																																																																					for (square[84] = 0; ((square[84] < 10) && (flag[83] == 1)); square[84]++) /*инициализцая 84 элемента*/
																																																																																					{
																																																																																						if ((square[84] != square[83]) && (square[84] != square[82]) && (square[84] != square[81]) && (square[84] != square[80]) && (square[84] != square[74]) && (square[84] != square[64]) && (square[84] != square[54]) && (square[84] != square[44]) && (square[84] != square[34]) && (square[84] != square[24]) && (square[84] != square[14]) && (square[84] != square[4]))
																																																																																						{
																																																																																							flag[84] = 1;
																																																																																						}
																																																																																						for (square[85] = 0; ((square[85] < 10) && (flag[84] == 1)); square[85]++) /*инициализцая 85 элемента*/
																																																																																						{
																																																																																							if ((square[85] != square[84]) && (square[85] != square[83]) && (square[85] != square[82]) && (square[85] != square[81]) && (square[85] != square[80]) && (square[85] != square[75]) && (square[85] != square[65]) && (square[85] != square[55]) && (square[85] != square[45]) && (square[85] != square[35]) && (square[85] != square[25]) && (square[85] != square[15]) && (square[85] != square[5]))
																																																																																							{
																																																																																								flag[85] = 1;
																																																																																							}
																																																																																							for (square[86] = 0; ((square[86] < 10) && (flag[85] == 1)); square[86]++) /*инициализцая 86 элемента*/
																																																																																							{
																																																																																								if ((square[86] != square[85]) && (square[86] != square[84]) && (square[86] != square[83]) && (square[86] != square[82]) && (square[86] != square[81]) && (square[86] != square[80]) && (square[86] != square[76]) && (square[86] != square[66]) && (square[86] != square[56]) && (square[86] != square[46]) && (square[86] != square[36]) && (square[86] != square[26]) && (square[86] != square[16]) && (square[86] != square[6]))
																																																																																								{
																																																																																									flag[86] = 1;
																																																																																								}
																																																																																								for (square[87] = 0; ((square[87] < 10) && (flag[86] == 1)); square[87]++) /*инициализцая 87 элемента*/
																																																																																								{
																																																																																									if ((square[87] != square[86]) && (square[87] != square[85]) && (square[87] != square[84]) && (square[87] != square[83]) && (square[87] != square[82]) && (square[87] != square[81]) && (square[87] != square[80]) && (square[87] != square[77]) && (square[87] != square[67]) && (square[87] != square[57]) && (square[87] != square[47]) && (square[87] != square[37]) && (square[87] != square[27]) && (square[87] != square[17]) && (square[87] != square[7]))
																																																																																									{
																																																																																										flag[87] = 1;
																																																																																									}
																																																																																									for (square[88] = 0; ((square[88] < 10) && (flag[87] == 1)); square[88]++) /*инициализцая 88 элемента*/
																																																																																									{
																																																																																										if ((square[88] != square[87]) && (square[88] != square[86]) && (square[88] != square[85]) && (square[88] != square[84]) && (square[88] != square[83]) && (square[88] != square[82]) && (square[88] != square[81]) && (square[88] != square[80]) && (square[88] != square[78]) && (square[88] != square[77]) && (square[88] != square[68]) && (square[88] != square[66]) && (square[88] != square[58]) && (square[88] != square[55]) && (square[88] != square[48]) && (square[88] != square[44]) && (square[88] != square[38]) && (square[88] != square[33]) && (square[88] != square[28]) && (square[88] != square[22]) && (square[88] != square[18]) && (square[88] != square[11]) && (square[88] != square[8]) && (square[88] != square[0]))
																																																																																										{
																																																																																											flag[88] = 1;
																																																																																										}
																																																																																										for (square[89] = 0; ((square[89] < 10) && (flag[88] == 1)); square[89]++) /*инициализцая 89 элемента*/
																																																																																										{
																																																																																											if ((square[89] != square[88]) && (square[89] != square[87]) && (square[89] != square[86]) && (square[89] != square[85]) && (square[89] != square[84]) && (square[89] != square[83]) && (square[89] != square[82]) && (square[89] != square[81]) && (square[89] != square[80]) && (square[89] != square[79]) && (square[89] != square[69]) && (square[89] != square[59]) && (square[89] != square[49]) && (square[89] != square[39]) && (square[89] != square[29]) && (square[89] != square[19]) && (square[89] != square[9]))
																																																																																											{
																																																																																												flag[89] = 1;
																																																																																											}
																																																																																											for (square[90] = 0; ((square[90] < 10) && (flag[89] == 1)); square[90]++) /*инициализцая 90 элемента*/
																																																																																											{
																																																																																												if ((square[90] != square[81]) && (square[90] != square[80]) && (square[90] != square[72]) && (square[90] != square[70]) && (square[90] != square[63]) && (square[90] != square[60]) && (square[90] != square[54]) && (square[90] != square[50]) && (square[90] != square[45]) && (square[90] != square[40]) && (square[90] != square[36]) && (square[90] != square[30]) && (square[90] != square[27]) && (square[90] != square[20]) && (square[90] != square[18]) && (square[90] != square[10]) && (square[90] != square[9]) && (square[90] != square[0]))
																																																																																												{
																																																																																													flag[90] = 1;
																																																																																												}
																																																																																												for (square[91] = 0; ((square[91] < 10) && (flag[90] == 1)); square[91]++) /*инициализцая 91 элемента*/
																																																																																												{
																																																																																													if ((square[91] != square[90]) && (square[91] != square[81]) && (square[91] != square[71]) && (square[91] != square[61]) && (square[91] != square[51]) && (square[91] != square[41]) && (square[91] != square[31]) && (square[91] != square[21]) && (square[91] != square[11]) && (square[91] != square[1]))
																																																																																													{
																																																																																														flag[91] = 1;
																																																																																													}
																																																																																													for (square[92] = 0; ((square[92] < 10) && (flag[91] == 1)); square[92]++) /*инициализцая 92 элемента*/
																																																																																													{
																																																																																														if ((square[92] != square[91]) && (square[92] != square[90]) && (square[92] != square[82]) && (square[92] != square[72]) && (square[92] != square[62]) && (square[92] != square[52]) && (square[92] != square[42]) && (square[92] != square[32]) && (square[92] != square[22]) && (square[92] != square[12]) && (square[92] != square[2]))
																																																																																														{
																																																																																															flag[92] = 1;
																																																																																														}
																																																																																														for (square[93] = 0; ((square[93] < 10) && (flag[92] == 1)); square[93]++) /*инициализцая 93 элемента*/
																																																																																														{
																																																																																															if ((square[93] != square[92]) && (square[93] != square[91]) && (square[93] != square[90]) && (square[93] != square[83]) && (square[93] != square[73]) && (square[93] != square[63]) && (square[93] != square[53]) && (square[93] != square[43]) && (square[93] != square[33]) && (square[93] != square[23]) && (square[93] != square[13]) && (square[93] != square[3]))
																																																																																															{
																																																																																																flag[93] = 1;
																																																																																															}
																																																																																															for (square[94] = 0; ((square[94] < 10) && (flag[93] == 1)); square[94]++) /*инициализцая 94 элемента*/
																																																																																															{
																																																																																																if ((square[94] != square[93]) && (square[94] != square[92]) && (square[94] != square[91]) && (square[94] != square[90]) && (square[94] != square[84]) && (square[94] != square[74]) && (square[94] != square[64]) && (square[94] != square[54]) && (square[94] != square[44]) && (square[94] != square[34]) && (square[94] != square[24]) && (square[94] != square[14]) && (square[94] != square[4]))
																																																																																																{
																																																																																																	flag[94] = 1;
																																																																																																}
																																																																																																for (square[95] = 0; ((square[95] < 10) && (flag[94] == 1)); square[95]++) /*инициализцая 95 элемента*/
																																																																																																{
																																																																																																	if ((square[95] != square[94]) && (square[95] != square[93]) && (square[95] != square[92]) && (square[95] != square[91]) && (square[95] != square[90]) && (square[95] != square[85]) && (square[95] != square[75]) && (square[95] != square[65]) && (square[95] != square[55]) && (square[95] != square[45]) && (square[95] != square[35]) && (square[95] != square[25]) && (square[95] != square[15]) && (square[95] != square[5]))
																																																																																																	{
																																																																																																		flag[95] = 1;
																																																																																																	}
																																																																																																	for (square[96] = 0; ((square[96] < 10) && (flag[95] == 1)); square[96]++) /*инициализцая 96 элемента*/
																																																																																																	{
																																																																																																		if ((square[96] != square[95]) && (square[96] != square[94]) && (square[96] != square[93]) && (square[96] != square[92]) && (square[96] != square[91]) && (square[96] != square[90]) && (square[96] != square[86]) && (square[96] != square[76]) && (square[96] != square[66]) && (square[96] != square[56]) && (square[96] != square[46]) && (square[96] != square[36]) && (square[96] != square[26]) && (square[96] != square[16]) && (square[96] != square[6]))
																																																																																																		{
																																																																																																			flag[96] = 1;
																																																																																																		}
																																																																																																		for (square[97] = 0; ((square[97] < 10) && (flag[96] == 1)); square[97]++) /*инициализцая 97 элемента*/
																																																																																																		{
																																																																																																			if ((square[97] != square[96]) && (square[97] != square[95]) && (square[97] != square[94]) && (square[97] != square[93]) && (square[97] != square[92]) && (square[97] != square[91]) && (square[97] != square[90]) && (square[97] != square[87]) && (square[97] != square[77]) && (square[97] != square[67]) && (square[97] != square[57]) && (square[97] != square[47]) && (square[97] != square[37]) && (square[97] != square[27]) && (square[97] != square[17]) && (square[97] != square[7]))
																																																																																																			{
																																																																																																				flag[97] = 1;
																																																																																																			}
																																																																																																			for (square[98] = 0; ((square[98] < 10) && (flag[97] == 1)); square[98]++) /*инициализцая 98 элемента*/
																																																																																																			{
																																																																																																				if ((square[98] != square[97]) && (square[98] != square[96]) && (square[98] != square[95]) && (square[98] != square[94]) && (square[98] != square[93]) && (square[98] != square[92]) && (square[98] != square[91]) && (square[98] != square[90]) && (square[98] != square[88]) && (square[98] != square[78]) && (square[98] != square[68]) && (square[98] != square[58]) && (square[98] != square[48]) && (square[98] != square[38]) && (square[98] != square[28]) && (square[98] != square[18]) && (square[98] != square[8]))
																																																																																																				{
																																																																																																					flag[98] = 1;
																																																																																																				}
																																																																																																				for (square[99] = 0; ((square[99] < 10) && (flag[98] == 1)); square[99]++) /*инициализцая 99 элемента*/
																																																																																																				{
																																																																																																					if ((square[99] != square[98]) && (square[99] != square[97]) && (square[99] != square[96]) && (square[99] != square[95]) && (square[99] != square[94]) && (square[99] != square[93]) && (square[99] != square[92]) && (square[99] != square[91]) && (square[99] != square[90]) && (square[99] != square[89]) && (square[99] != square[88]) && (square[99] != square[79]) && (square[99] != square[77]) && (square[99] != square[69]) && (square[99] != square[66]) && (square[99] != square[59]) && (square[99] != square[55]) && (square[99] != square[49]) && (square[99] != square[44]) && (square[99] != square[39]) && (square[99] != square[33]) && (square[99] != square[29]) && (square[99] != square[22]) && (square[99] != square[19]) && (square[99] != square[11]) && (square[99] != square[9]) && (square[99] != square[0]))
																																																																																																					{
																																																																																																						/* ДЛК сгенерирован*/
																																																																																																						count++;
																																																																																																						/* В данный момент времени в square[0]-square[99] находится ДЛК */

																																																																																																						processNewDLS(odls_pair_vec, fragment_index, square);
																																																																																																					}
																																																																																																				}
																																																																																																				flag[98] = 0;
																																																																																																			}
																																																																																																			flag[97] = 0;
																																																																																																		}
																																																																																																		flag[96] = 0;
																																																																																																	}
																																																																																																	flag[95] = 0;
																																																																																																}
																																																																																																flag[94] = 0;
																																																																																															}
																																																																																															flag[93] = 0;
																																																																																														}
																																																																																														flag[92] = 0;
																																																																																													}
																																																																																													flag[91] = 0;
																																																																																												}
																																																																																												flag[90] = 0;
																																																																																											}
																																																																																											flag[89] = 0;
																																																																																										}
																																																																																										flag[88] = 0;
																																																																																									}
																																																																																									flag[87] = 0;
																																																																																								}
																																																																																								flag[86] = 0;
																																																																																							}
																																																																																							flag[85] = 0;
																																																																																						}
																																																																																						flag[84] = 0;
																																																																																					}
																																																																																					flag[83] = 0;
																																																																																				}
																																																																																				flag[82] = 0;
																																																																																			}
																																																																																			flag[81] = 0;
																																																																																		}
																																																																																		flag[80] = 0;
																																																																																	}
																																																																																	flag[79] = 0;
																																																																																}
																																																																																flag[78] = 0;
																																																																															}
																																																																															flag[77] = 0;
																																																																														}
																																																																														flag[76] = 0;
																																																																													}
																																																																													flag[75] = 0;
																																																																												}
																																																																												flag[74] = 0;
																																																																											}
																																																																											flag[73] = 0;
																																																																										}
																																																																										flag[72] = 0;
																																																																									}
																																																																									flag[71] = 0;
																																																																								}
																																																																								flag[70] = 0;
																																																																							}
																																																																							flag[69] = 0;
																																																																						}
																																																																						flag[68] = 0;
																																																																					}
																																																																					flag[67] = 0;
																																																																				}
																																																																				flag[66] = 0;
																																																																			}
																																																																			flag[65] = 0;
																																																																		}
																																																																		flag[64] = 0;
																																																																	}
																																																																	flag[63] = 0;
																																																																}
																																																																flag[62] = 0;
																																																															}
																																																															flag[61] = 0;
																																																														}
																																																														flag[60] = 0;
																																																													}
																																																													flag[59] = 0;
																																																												}
																																																												flag[58] = 0;
																																																											}
																																																											flag[57] = 0;
																																																										}
																																																										flag[56] = 0;
																																																									}
																																																									flag[55] = 0;
																																																								}
																																																								flag[54] = 0;
																																																							}
																																																							flag[53] = 0;
																																																						}
																																																						flag[52] = 0;
																																																					}
																																																					flag[51] = 0;
																																																				}
																																																				flag[50] = 0;
																																																			}
																																																			flag[49] = 0;
																																																		}
																																																		flag[48] = 0;
																																																	}
																																																	flag[47] = 0;
																																																}
																																																flag[46] = 0;
																																															}
																																															flag[45] = 0;
																																														}
																																														flag[44] = 0;
																																													}
																																													flag[43] = 0;
																																												}
																																												flag[42] = 0;
																																											}
																																											flag[41] = 0;
																																										}
																																										flag[40] = 0;
																																									}
																																									flag[39] = 0;
																																								}
																																								flag[38] = 0;
																																							}
																																							flag[37] = 0;
																																						}
																																						flag[36] = 0;
																																					}
																																					flag[35] = 0;
																																				}
																																				flag[34] = 0;
																																			}
																																			flag[33] = 0;
																																		}
																																		flag[32] = 0;
																																	}
																																	flag[31] = 0;
																																}
																																flag[30] = 0;
																															}
																															flag[29] = 0;
																														}
																														flag[28] = 0;
																													}
																													flag[27] = 0;
																												}
																												flag[26] = 0;
																											}
																											flag[25] = 0;
																										}
																										flag[24] = 0;
																									}
																									flag[23] = 0;
																								}
																								flag[22] = 0;
																							}
																							flag[21] = 0;
																						}
																						flag[20] = 0;
																					}
																					flag[19] = 0;
																				}
																				flag[18] = 0;
																			}
																			flag[17] = 0;
																		}
																		flag[16] = 0;
																	}
																	flag[15] = 0;
																}
																flag[14] = 0;
															}
															flag[13] = 0;
														}
														flag[12] = 0;
													}
													flag[11] = 0;
												}

											}
											flag[9] = 0;
										}
										flag[8] = 0;
									}
									flag[7] = 0;
								}
								flag[6] = 0;
							}
							flag[5] = 0;
						}
						flag[4] = 0;
					}
					flag[3] = 0;
				}
				flag[2] = 0;
			}
		}
	}
#endif
	return 0;
}