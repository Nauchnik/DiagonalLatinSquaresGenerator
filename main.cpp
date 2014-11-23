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

using namespace std;

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
double limSecondsOneSquare = 100;
bool isTimeBreak = false;
int const numberOfIntsBase = 20;
map<int, int>freeNumbers[numberOfIntsBase][numberOfIntsBase];
map<int, int>busyNumbers[numberOfIntsBase][numberOfIntsBase];
map<int, int>diagonalBusyNumbers[numberOfIntsBase][numberOfIntsBase];
int valuesInTable[numberOfIntsBase][numberOfIntsBase];

int main(int argc, char **argv)
{
	unsigned long long countOfCalc = 1;
	
	if ( argc == 1 ) {
		std::cerr << "Usage : LS order LS count lim_seconds";
		return 1;
	}

	if (argc > 1)
		numberOfInts = atoi(argv[1]);

	if (argc > 2)
		countOfCalc = atoi(argv[2]);

	if (argc > 3)
		limSecondsOneSquare = atof(argv[3]);

	//std::cout << "limSecondsOneSquare " << limSecondsOneSquare << std::endl;
	
	//srand((unsigned)time(0));

	cout << endl;
	double sum_time = 0.0, min_time = 604800.0, max_time = 0.0;
	std::chrono::high_resolution_clock::time_point t1, t2;
	std::chrono::duration<double> time_span;
	double solving_time;

	// set of unique diagonal Latin squares
	std::set<std::string> dls_string_set;
	std::string cur_string_set;
	cur_string_set.resize( numberOfInts * numberOfInts );
	unsigned k;
	std::stringstream sstream;
	unsigned long long genereated_count = 0;
	
	for (unsigned long long i = 0; i < countOfCalc; i++) {
		//time_t before = time(NULL);
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
		for (int j1 = 0; j1 < numberOfInts; j1++)
			for (int j2 = 0; j2 < numberOfInts; j2++)
				sstream << valuesInTable[j1][j2] - 1; // save values from 0 to 9
		cur_string_set = sstream.str();
		cout << cur_string_set << endl;
		dls_string_set.insert( cur_string_set );
		sstream.clear(); sstream.str("");
		
		//Print();
		
		min_time = min_time < solving_time ? min_time : solving_time;
		max_time = max_time > solving_time ? max_time : solving_time;
		sum_time += solving_time;

		cout << "time in seconds : " << solving_time << endl;
		cout << "processed " << i+1 << " from " << countOfCalc << endl;
		cout << "generated " << genereated_count << endl;
		cout << "unique DLS : " << dls_string_set.size() << endl;
		cout << "current median time " << sum_time / double(i + 1) << endl;
		cout << "current min_time " << min_time << endl;
		cout << "current max_time " << max_time << endl;
		cout << endl;
	}
	
	system("pause");

	return 0;
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