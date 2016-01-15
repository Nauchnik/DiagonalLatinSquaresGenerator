#include <stdio.h>

#define NUM_OF_PARTS 2445393
#define NUM_OF_DLK_IN_ONE_CHECK 1000000
#define DIFFER_GLOBAL_LOCAL 5

unsigned long long int generate_dlk(int parts, int part)
{
unsigned short int square[100] = { 0 };
unsigned short int stl[10][10] ={ 1 };
unsigned short int str[10][10] = { 1 };
unsigned short int diag[2][10] = { 1 };
unsigned short int flag[100] = { 0 };

unsigned short int start[100] = { 0 };
unsigned short int end[100] = { 9 };

unsigned long long int number_min=0;
unsigned long long int number_max=0;
unsigned long long int number_real=0;

unsigned short int result=0;
unsigned short int local_max=0;
unsigned short int global_max=0;

unsigned long long int count = 0;
int i=0;
int j=0;
int number_of_comb=NUM_OF_PARTS;
int number_of_comb_in_one_part=0;
float fparts=parts;
fparts=number_of_comb/fparts;


//округление в +//
int inter;
float finter;
inter=fparts;
finter=inter;
if (fparts!=finter)
fparts=finter+1;
else
fparts=finter;
number_of_comb_in_one_part=fparts;


//обнуление квадратов и флагов перед генерацией областей ДЛК
for (i=0;i<100;i++)
{
square[i]=0;
flag[i]=0;
}
count = 0;

for(i=0;i<10;i++)
{
	for(j=0;j<10;j++)
	{
		str[i][j]=1;
		stl[i][j]=1;
	}
}

for(i=0;i<2;i++)
{
	for(j=0;j<10;j++)
	{
		diag[i][j]=1;
	}
}

//поиск пороговых значений 
for (square[0] = 0; square[0] == 0; square[0]++) /*инициализцая 0 элемента*/
{
if(diag[0][square[0]] && str[0][square[0]] && stl[0][square[0]])
{
	str[0][square[0]]=0;
	stl[0][square[0]]=0;
	diag[0][square[0]]=0;
	flag[0] = 1;
}
for (square[1] = 1; (flag[0] && (square[1] == 1)); square[1]++) /*инициализцая 1 элемента*/
{
if(str[0][square[1]] && stl[1][square[1]])
{
	str[0][square[1]]=0;
	stl[1][square[1]]=0;
	flag[1] = 1;
}
for (square[2] = 2; (flag[1] && (square[2] == 2)); square[2]++)  /*инициализцая 2 элемента*/
{
if (str[0][square[2]] && stl[2][square[2]])
{
	str[0][square[2]]=0;
	stl[2][square[2]]=0;
	flag[2] = 1;
}
for (square[3] = 3; (flag[2] && (square[3] == 3)); square[3]++)  /*инициализцая 3 элемента*/
{
if (str[0][square[3]] && stl[3][square[3]])
{
	str[0][square[3]]=0;
	stl[3][square[3]]=0;
	flag[3] = 1;
}
for (square[4] = 4; (flag[3] && (square[4] == 4)); square[4]++) /*инициализцая 4 элемента*/
{
if (str[0][square[4]] && stl[4][square[4]])
{
	str[0][square[4]]=0;
	stl[4][square[4]]=0;
	flag[4] = 1;
}
for (square[5] = 5; (flag[4] && (square[5] == 5)); square[5]++) /*инициализцая 5 элемента*/
{
if (str[0][square[5]] && stl[5][square[5]])
{
	str[0][square[5]]=0;
	stl[5][square[5]]=0;
	flag[5] = 1;
}
for (square[6] = 6; (flag[5] && (square[6] == 6)); square[6]++) /*инициализцая 6 элемента*/
{
if (str[0][square[6]] && stl[6][square[6]])
{
	str[0][square[6]]=0;
	stl[6][square[6]]=0;
	flag[6] = 1;
}
for (square[7] = 7; (flag[6] && (square[7] == 7)); square[7]++) /*инициализцая 7 элемента*/
{
if (str[0][square[7]] && stl[7][square[7]])
{
	str[0][square[7]]=0;
	stl[7][square[7]]=0;
	flag[7] = 1;
}
for (square[8] = 8; (flag[7] && (square[8] == 8)); square[8]++) /*инициализцая 8 элемента*/
{
if (str[0][square[8]] && stl[8][square[8]])
{
	str[0][square[8]]=0;
	stl[8][square[8]]=0;
	flag[8] = 1;
}
for (square[9] = 9; (flag[8] && (square[9] == 9)); square[9]++) /*инициализцая 9 элемента*/
{
if (str[0][square[9]] && stl[9][square[9]] && diag[1][square[9]])
{
	str[0][square[9]]=0;
	stl[9][square[9]]=0;
	diag[1][square[9]]=0;
	flag[9] = 1;
}
for (square[11] = 0; (flag[9] && (square[11] < 10)); square[11]++) /*инициализцая 10 элемента*/
{
if (diag[0][square[11]] && stl[1][square[11]] && str[1][square[11]])
{
	str[1][square[11]]=0;
	stl[1][square[11]]=0;
	diag[0][square[11]]=0;
	flag[11] = 1;
}
for (square[18] = 0; (flag[11] && (square[18] < 10)); square[18]++) /*инициализцая 11 элемента*/
{
if (diag[1][square[18]] && stl[8][square[18]] && str[1][square[18]])
{
	str[1][square[18]]=0;
	stl[8][square[18]]=0;
	diag[1][square[18]]=0;
	flag[18] = 1;
}
for (square[22] = 0; (flag[18] && (square[22] < 10)); square[22]++) /*инициализцая 12 элемента*/
{
if (diag[0][square[22]] && stl[2][square[22]] && str[2][square[22]])
{
	str[2][square[22]]=0;
	stl[2][square[22]]=0;
	diag[0][square[22]]=0;
	flag[22] = 1;
}
for (square[27] = 0; (flag[22] && (square[27] < 10)); square[27]++) /*инициализцая 13 элемента*/
{
if (diag[1][square[27]] && stl[7][square[27]] && str[2][square[27]])
{
	str[2][square[27]]=0;
	stl[7][square[27]]=0;
	diag[1][square[27]]=0;
	flag[27] = 1;
}
for (square[33] = 0; (flag[27] && (square[33] < 10)); square[33]++) /*инициализцая 14 элемента*/
{
if (diag[0][square[33]] && stl[3][square[33]] && str[3][square[33]])
{
	str[3][square[33]]=0;
	stl[3][square[33]]=0;
	diag[0][square[33]]=0;
	flag[33] = 1;
}
for (square[36] = 0; (flag[33] && (square[36] < 10)); square[36]++) /*инициализцая 15 элемента*/
{
if (diag[1][square[36]] && stl[6][square[36]] && str[3][square[36]])
{
	str[3][square[36]]=0;
	stl[6][square[36]]=0;
	diag[1][square[36]]=0;
	flag[36] = 1;
}
for (square[44] = 0; (flag[36] && (square[44] < 10)); square[44]++) /*инициализцая 16 элемента*/
{
if (diag[0][square[44]] && stl[4][square[44]] && str[4][square[44]])
{
	str[4][square[44]]=0;
	stl[4][square[44]]=0;
	diag[0][square[44]]=0;
	flag[44] = 1;
}
for (square[45] = 0; (flag[44] && (square[45] < 10)); square[45]++) /*инициализцая 17 элемента*/
{
if (diag[1][square[45]] && stl[5][square[45]] && str[4][square[45]])
{
	str[4][square[45]]=0;
	stl[5][square[45]]=0;
	diag[1][square[45]]=0;
	flag[45] = 1;
}

if(flag[45])
{
	
if (count==(number_of_comb_in_one_part*part))
{
	for(j=0;j<100;j++)
	{
	start[j]=square[j];
	}
}

if((count>=number_of_comb_in_one_part*part) && (count <= (number_of_comb_in_one_part*(part+1)-1)))
{
for(j=0;j<100;j++)
{
end[j]=square[j];
}
}

count++;
}

if(flag[45])
{
	str[4][square[45]]=1;
	stl[5][square[45]]=1;
	diag[1][square[45]]=1;
	flag[45] = 0;
}
}
if(flag[44])
{
	str[4][square[44]]=1;
	stl[4][square[44]]=1;
	diag[0][square[44]]=1;
	flag[44] = 0;
}
}
if(flag[36])
{
	str[3][square[36]]=1;
	stl[6][square[36]]=1;
	diag[1][square[36]]=1;
	flag[36] = 0;
}
}
if(flag[33])
{
	str[3][square[33]]=1;
	stl[3][square[33]]=1;
	diag[0][square[33]]=1;
	flag[33] = 0;
}
}
if(flag[27])
{
	str[2][square[27]]=1;
	stl[7][square[27]]=1;
	diag[1][square[27]]=1;
	flag[27] = 0;
}
}
if(flag[22])
{
	str[2][square[22]]=1;
	stl[2][square[22]]=1;
	diag[0][square[22]]=1;
	flag[22] = 0;
}
}
if(flag[18])
{
	str[1][square[18]]=1;
	stl[8][square[18]]=1;
	diag[1][square[18]]=1;
	flag[18] = 0;
}
}
if(flag[11])
{
	str[1][square[11]]=1;
	stl[1][square[11]]=1;
	diag[0][square[11]]=1;
	flag[11] = 0;
}
}
if(flag[9])
{
	str[0][square[9]]=1;
	stl[9][square[9]]=1;
	diag[1][square[9]]=1;
	flag[9] = 0;
}
}
if(flag[8])
{
	str[0][square[8]]=1;
	stl[8][square[8]]=1;
	flag[8] = 0;
}
}
if(flag[7])
{
	str[0][square[7]]=1;
	stl[7][square[7]]=1;
	flag[7] = 0;
}
}
if(flag[6])
{
	str[0][square[6]]=1;
	stl[6][square[6]]=1;
	flag[6] = 0;
}
}
if(flag[5])
{
	str[0][square[5]]=1;
	stl[5][square[5]]=1;
	flag[5] = 0;
}
}
if(flag[4])
{
	str[0][square[4]]=1;
	stl[4][square[4]]=1;
	flag[4] = 0;
}
}
if(flag[3])
{
	str[0][square[3]]=1;
	stl[3][square[3]]=1;
	flag[3] = 0;
}
}
if(flag[2])
{
	str[0][square[2]]=1;
	stl[2][square[2]]=1;
	flag[2] = 0;
}
}
if(flag[1])
{
	str[0][square[1]]=1;
	stl[1][square[1]]=1;
	flag[1] = 0;
}
}
if(flag[0])
{
	str[0][square[0]]=1;
	stl[0][square[0]]=1;
	diag[0][square[0]]=1;
	flag[0] = 0;
}
}

//Подсчет пороговых значний порогов
number_min=(start[11]*1000000000)+(start[18]*100000000)+(start[22]*10000000)+(start[27]*1000000)+(start[33]*100000)+(start[36]*10000)+(start[44]*1000)+(start[45]*100); // +(start[18]*10)+start[19];
number_max=(end[11]*1000000000)+(end[18]*100000000)+(end[22]*10000000)+(end[27]*1000000)+(end[33]*100000)+(end[36]*10000)+(end[44]*1000)+(end[45]*100); //+(end[18]*10)+end[19];


//обнуление квадратов и флагов перед генерацией ДЛК
for (i=0;i<100;i++)
{
square[i]=0;
flag[i]=0;
}
count = 0;

for(i=0;i<10;i++)
{
	for(j=0;j<10;j++)
	{
		str[i][j]=1;
		stl[i][j]=1;
	}
}

for(i=0;i<2;i++)
{
	for(j=0;j<10;j++)
	{
		diag[i][j]=1;
	}
}


//генерация ДЛК
for (square[0] = 0; square[0] == 0; square[0]++) /*инициализцая 0 элемента*/
{
if(diag[0][square[0]] && str[0][square[0]] && stl[0][square[0]])
{
	str[0][square[0]]=0;
	stl[0][square[0]]=0;
	diag[0][square[0]]=0;
	flag[0] = 1;
}
for (square[1] = 1; (flag[0] && (square[1] == 1)); square[1]++) /*инициализцая 1 элемента*/
{
if(str[0][square[1]] && stl[1][square[1]])
{
	str[0][square[1]]=0;
	stl[1][square[1]]=0;
	flag[1] = 1;
}
for (square[2] = 2; (flag[1] && (square[2] == 2)); square[2]++)  /*инициализцая 2 элемента*/
{
if (str[0][square[2]] && stl[2][square[2]])
{
	str[0][square[2]]=0;
	stl[2][square[2]]=0;
	flag[2] = 1;
}
for (square[3] = 3; (flag[2] && (square[3] == 3)); square[3]++)  /*инициализцая 3 элемента*/
{
if (str[0][square[3]] && stl[3][square[3]])
{
	str[0][square[3]]=0;
	stl[3][square[3]]=0;
	flag[3] = 1;
}
for (square[4] = 4; (flag[3] && (square[4] == 4)); square[4]++) /*инициализцая 4 элемента*/
{
if (str[0][square[4]] && stl[4][square[4]])
{
	str[0][square[4]]=0;
	stl[4][square[4]]=0;
	flag[4] = 1;
}
for (square[5] = 5; (flag[4] && (square[5] == 5)); square[5]++) /*инициализцая 5 элемента*/
{
if (str[0][square[5]] && stl[5][square[5]])
{
	str[0][square[5]]=0;
	stl[5][square[5]]=0;
	flag[5] = 1;
}
for (square[6] = 6; (flag[5] && (square[6] == 6)); square[6]++) /*инициализцая 6 элемента*/
{
if (str[0][square[6]] && stl[6][square[6]])
{
	str[0][square[6]]=0;
	stl[6][square[6]]=0;
	flag[6] = 1;
}
for (square[7] = 7; (flag[6] && (square[7] == 7)); square[7]++) /*инициализцая 7 элемента*/
{
if (str[0][square[7]] && stl[7][square[7]])
{
	str[0][square[7]]=0;
	stl[7][square[7]]=0;
	flag[7] = 1;
}
for (square[8] = 8; (flag[7] && (square[8] == 8)); square[8]++) /*инициализцая 8 элемента*/
{
if (str[0][square[8]] && stl[8][square[8]])
{
	str[0][square[8]]=0;
	stl[8][square[8]]=0;
	flag[8] = 1;
}
for (square[9] = 9; (flag[8] && (square[9] == 9)); square[9]++) /*инициализцая 9 элемента*/
{
if (str[0][square[9]] && stl[9][square[9]] && diag[1][square[9]])
{
	str[0][square[9]]=0;
	stl[9][square[9]]=0;
	diag[1][square[9]]=0;
	flag[9] = 1;
}
for (square[11] = 0; (flag[9] && (square[11] < 10)); square[11]++) /*инициализцая 10 элемента*/
{
if (diag[0][square[11]] && stl[1][square[11]] && str[1][square[11]])
{
	str[1][square[11]]=0;
	stl[1][square[11]]=0;
	diag[0][square[11]]=0;
	flag[11] = 1;
}
for (square[18] = 0; (flag[11] && (square[18] < 10)); square[18]++) /*инициализцая 11 элемента*/
{
if (diag[1][square[18]] && stl[8][square[18]] && str[1][square[18]])
{
	str[1][square[18]]=0;
	stl[8][square[18]]=0;
	diag[1][square[18]]=0;
	flag[18] = 1;
}
for (square[22] = 0; (flag[18] && (square[22] < 10)); square[22]++) /*инициализцая 12 элемента*/
{
if (diag[0][square[22]] && stl[2][square[22]] && str[2][square[22]])
{
	str[2][square[22]]=0;
	stl[2][square[22]]=0;
	diag[0][square[22]]=0;
	flag[22] = 1;
}
for (square[27] = 0; (flag[22] && (square[27] < 10)); square[27]++) /*инициализцая 13 элемента*/
{
if (diag[1][square[27]] && stl[7][square[27]] && str[2][square[27]])
{
	str[2][square[27]]=0;
	stl[7][square[27]]=0;
	diag[1][square[27]]=0;
	flag[27] = 1;
}
for (square[33] = 0; (flag[27] && (square[33] < 10)); square[33]++) /*инициализцая 14 элемента*/
{
if (diag[0][square[33]] && stl[3][square[33]] && str[3][square[33]])
{
	str[3][square[33]]=0;
	stl[3][square[33]]=0;
	diag[0][square[33]]=0;
	flag[33] = 1;
}
for (square[36] = 0; (flag[33] && (square[36] < 10)); square[36]++) /*инициализцая 15 элемента*/
{
if (diag[1][square[36]] && stl[6][square[36]] && str[3][square[36]])
{
	str[3][square[36]]=0;
	stl[6][square[36]]=0;
	diag[1][square[36]]=0;
	flag[36] = 1;
}
for (square[44] = 0; (flag[36] && (square[44] < 10)); square[44]++) /*инициализцая 16 элемента*/
{
if (diag[0][square[44]] && stl[4][square[44]] && str[4][square[44]])
{
	str[4][square[44]]=0;
	stl[4][square[44]]=0;
	diag[0][square[44]]=0;
	flag[44] = 1;
}
for (square[45] = 0; (flag[44] && (square[45] < 10)); square[45]++) /*инициализцая 17 элемента*/
{
if (diag[1][square[45]] && stl[5][square[45]] && str[4][square[45]])
{
	str[4][square[45]]=0;
	stl[5][square[45]]=0;
	diag[1][square[45]]=0;
	flag[45] = 1;
	number_real=(square[11]*1000000000)+(square[18]*100000000)+(square[22]*10000000)+(square[27]*1000000)+(square[33]*100000)+(square[36]*10000)+(square[44]*1000)+(square[45]*100);
}
for (square[54] = 0; (flag[45] && (square[54] < 10) && (number_real>=number_min) && (number_real<=number_max)); square[54]++) /*инициализцая 18 элемента*/
{
if (diag[1][square[54]] && stl[4][square[54]] && str[5][square[54]])
{
	str[5][square[54]]=0;
	stl[4][square[54]]=0;
	diag[1][square[54]]=0;
	flag[54] = 1;
}
for (square[55] = 0; (flag[54] && (square[55] < 10)); square[55]++) /*инициализцая 19 элемента*/
{
if (diag[0][square[55]] && stl[5][square[55]] && str[5][square[55]])
{
	str[5][square[55]]=0;
	stl[5][square[55]]=0;
	diag[0][square[55]]=0;
	flag[55] = 1;
}
for (square[63] = 0; (flag[55] && (square[63] < 10)); square[63]++) /*инициализцая 20 элемента*/
{
if (diag[1][square[63]] && stl[3][square[63]] && str[6][square[63]])
{
	str[6][square[63]]=0;
	stl[3][square[63]]=0;
	diag[1][square[63]]=0;
	flag[63] = 1;
}
for (square[66] = 0; (flag[63] && (square[66] < 10)); square[66]++) /*инициализцая 21 элемента*/
{
if (diag[0][square[66]] && stl[6][square[66]] && str[6][square[66]])
{
	str[6][square[66]]=0;
	stl[6][square[66]]=0;
	diag[0][square[66]]=0;
	flag[66] = 1;
}
for (square[72] = 0; (flag[66] && (square[72] < 10)); square[72]++) /*инициализцая 22 элемента*/
{
if (diag[1][square[72]] && stl[2][square[72]] && str[7][square[72]])
{
	str[7][square[72]]=0;
	stl[2][square[72]]=0;
	diag[1][square[72]]=0;
	flag[72] = 1;
}
for (square[77] = 0; (flag[72] && (square[77] < 10)); square[77]++) /*инициализцая 23 элемента*/
{
if (diag[0][square[77]] && stl[7][square[77]] && str[7][square[77]])
{
	str[7][square[77]]=0;
	stl[7][square[77]]=0;
	diag[0][square[77]]=0;
	flag[77] = 1;
}
for (square[81] = 0; (flag[77] && (square[81] < 10)); square[81]++) /*инициализцая 24 элемента*/
{
if (diag[1][square[81]] && stl[1][square[81]] && str[8][square[81]])
{
	str[8][square[81]]=0;
	stl[1][square[81]]=0;
	diag[1][square[81]]=0;
	flag[81] = 1;
}
for (square[88] = 0; (flag[81] && (square[88] < 10)); square[88]++) /*инициализцая 25 элемента*/
{
if (diag[0][square[88]] && stl[8][square[88]] && str[8][square[88]])
{
	str[8][square[88]]=0;
	stl[8][square[88]]=0;
	diag[0][square[88]]=0;
	flag[88] = 1;
}
for (square[90] = 0; (flag[88] && (square[90] < 10)); square[90]++) /*инициализцая 26 элемента*/
{
if (diag[1][square[90]] && stl[0][square[90]] && str[9][square[90]])
{
	str[9][square[90]]=0;
	stl[0][square[90]]=0;
	diag[1][square[90]]=0;
	flag[90] = 1;
}
for (square[99] = 0; (flag[90] && (square[99] < 10)); square[99]++) /*инициализцая 27 элемента*/
{
if (diag[0][square[99]] && stl[9][square[99]] && str[9][square[99]])
{
	str[9][square[99]]=0;
	stl[9][square[99]]=0;
	diag[0][square[99]]=0;
	flag[99] = 1;
}
for (square[10] = 0; (flag[99] && (square[10] < 10)); square[10]++) /*инициализцая 28 элемента*/
{
if (stl[0][square[10]] && str[1][square[10]])
{
	str[1][square[10]]=0;
	stl[0][square[10]]=0;
	flag[10] = 1;
}
for (square[12] = 0; (flag[10] && (square[12] < 10)); square[12]++) /*инициализцая 29 элемента*/
{
if (str[1][square[12]] && stl[2][square[12]])
{
	str[1][square[12]]=0;
	stl[2][square[12]]=0;
	flag[12] = 1;
}
for (square[13] = 0; (flag[12] && (square[13] < 10)); square[13]++) /*инициализцая 30 элемента*/
{
if (str[1][square[13]] && stl[3][square[13]])
{
	str[1][square[13]]=0;
	stl[3][square[13]]=0;
	flag[13] = 1;
}
for (square[14] = 0; (flag[13] && (square[14] < 10)); square[14]++) /*инициализцая 14 элемента*/
{
if (str[1][square[14]] && stl[4][square[14]])
{
	str[1][square[14]]=0;
	stl[4][square[14]]=0;
	flag[14] = 1;
}
for (square[15] = 0; (flag[14] && (square[15] < 10)); square[15]++) /*инициализцая 15 элемента*/
{
if (str[1][square[15]] && stl[5][square[15]])
{
	str[1][square[15]]=0;
	stl[5][square[15]]=0;
	flag[15] = 1;
}
for (square[16] = 0; (flag[15] && (square[16] < 10)); square[16]++) /*инициализцая 16 элемента*/
{
if (str[1][square[16]] && stl[6][square[16]])
{
	str[1][square[16]]=0;
	stl[6][square[16]]=0;	
	flag[16] = 1;
}
for (square[17] = 0; (flag[16] && (square[17] < 10)); square[17]++) /*инициализцая 17 элемента*/
{
if (str[1][square[17]] && stl[7][square[17]])
{
	str[1][square[17]]=0;
	stl[7][square[17]]=0;
	flag[17] = 1;
}
for (square[19] = 0; (flag[17] && (square[19] < 10)); square[19]++) /*инициализцая 35 элемента*/
{
if (str[1][square[19]] && stl[9][square[19]])
{
	str[1][square[19]]=0;
	stl[9][square[19]]=0;
	flag[19] = 1;
}
for (square[20] = 0; (flag[19] && (square[20] < 10)); square[20]++) /*инициализцая 36 элемента*/
{
if (stl[0][square[20]] && str[2][square[20]])
{
	str[2][square[20]]=0;
	stl[0][square[20]]=0;
	flag[20] = 1;
}
for (square[21] = 0; (flag[20] && (square[21] < 10)); square[21]++) /*инициализцая 37 элемента*/
{
if (stl[1][square[21]] && str[2][square[21]])
{
	str[2][square[21]]=0;
	stl[1][square[21]]=0;
	flag[21] = 1;
}
for (square[23] = 0; (flag[21] && (square[23] < 10)); square[23]++) /*инициализцая 38 элемента*/
{
if (str[2][square[23]] && stl[3][square[23]])
{
	str[2][square[23]]=0;
	stl[3][square[23]]=0;
	flag[23] = 1;
}
for (square[24] = 0; (flag[23] && (square[24] < 10)); square[24]++) /*инициализцая 24 элемента*/
{
if (str[2][square[24]] && stl[4][square[24]])
{
	str[2][square[24]]=0;
	stl[4][square[24]]=0;
	flag[24] = 1;
}
for (square[25] = 0; (flag[24] && (square[25] < 10)); square[25]++) /*инициализцая 25 элемента*/
{
if (str[2][square[25]] && stl[5][square[25]])
{
	str[2][square[25]]=0;
	stl[5][square[25]]=0;
	flag[25] = 1;
}
for (square[26] = 0; (flag[25] && (square[26] < 10)); square[26]++) /*инициализцая 26 элемента*/
{
if (str[2][square[26]] && stl[6][square[26]])
{
	str[2][square[26]]=0;
	stl[6][square[26]]=0;
	flag[26] = 1;
}
for (square[28] = 0; (flag[26] && (square[28] < 10)); square[28]++) /*инициализцая 42 элемента*/
{
if (str[2][square[28]] && stl[8][square[28]])
{
	str[2][square[28]]=0;
	stl[8][square[28]]=0;
	flag[28] = 1;
}
for (square[29] = 0; (flag[28] && (square[29] < 10)); square[29]++) /*инициализцая 29 элемента*/
{
if (str[2][square[29]] && stl[9][square[29]])
{
	str[2][square[29]]=0;
	stl[9][square[29]]=0;
	flag[29] = 1;
}
for (square[30] = 0; (flag[29] && (square[30] < 10)); square[30]++) /*инициализцая 30 элемента*/
{
if (stl[0][square[30]] && str[3][square[30]])
{
	str[3][square[30]]=0;
	stl[0][square[30]]=0;
	flag[30] = 1;
}
for (square[31] = 0; (flag[30] && (square[31] < 10)); square[31]++) /*инициализцая 31 элемента*/
{
if (stl[1][square[31]] && str[3][square[31]])
{
	str[3][square[31]]=0;
	stl[1][square[31]]=0;
	flag[31] = 1;
}
for (square[32] = 0; (flag[31] && (square[32] < 10)); square[32]++) /*инициализцая 32 элемента*/
{
if (stl[2][square[32]] && str[3][square[32]])
{
	str[3][square[32]]=0;
	stl[2][square[32]]=0;
	flag[32] = 1;
}
for (square[34] = 0; (flag[32] && (square[34] < 10)); square[34]++) /*инициализцая 47 элемента*/
{
if (str[3][square[34]] && stl[4][square[34]])
{
	str[3][square[34]]=0;
	stl[4][square[34]]=0;
	flag[34] = 1;
}
for (square[35] = 0; (flag[34] && (square[35] < 10)); square[35]++) /*инициализцая 35 элемента*/
{
if (str[3][square[35]] && stl[5][square[35]])
{
	str[3][square[35]]=0;
	stl[5][square[35]]=0;
	flag[35] = 1;
}
for (square[37] = 0; (flag[35] && (square[37] < 10)); square[37]++) /*инициализцая 49 элемента*/
{
if (str[3][square[37]] && stl[7][square[37]])
{
	str[3][square[37]]=0;
	stl[7][square[37]]=0;
	flag[37] = 1;
}
for (square[38] = 0; (flag[37] && (square[38] < 10)); square[38]++) /*инициализцая 38 элемента*/
{
if (str[3][square[38]] && stl[8][square[38]])
{
	str[3][square[38]]=0;
	stl[8][square[38]]=0;
	flag[38] = 1;
}
for (square[39] = 0; (flag[38] && (square[39] < 10)); square[39]++) /*инициализцая 39 элемента*/
{
if (str[3][square[39]] && stl[9][square[39]])
{
	str[3][square[39]]=0;
	stl[9][square[39]]=0;
	flag[39] = 1;
}
for (square[40] = 0; (flag[39] && (square[40] < 10)); square[40]++) /*инициализцая 40 элемента*/
{
if (stl[0][square[40]] && str[4][square[40]])
{
	str[4][square[40]]=0;
	stl[0][square[40]]=0;
	flag[40] = 1;
}
for (square[41] = 0; (flag[40] && (square[41] < 10)); square[41]++) /*инициализцая 41 элемента*/
{
if (stl[1][square[41]] && str[4][square[41]])
{
	str[4][square[41]]=0;
	stl[1][square[41]]=0;
	flag[41] = 1;
}
for (square[42] = 0; (flag[41] && (square[42] < 10)); square[42]++) /*инициализцая 42 элемента*/
{
if (stl[2][square[42]] && str[4][square[42]])
{
	str[4][square[42]]=0;
	stl[2][square[42]]=0;
	flag[42] = 1;
}
for (square[43] = 0; (flag[42] && (square[43] < 10)); square[43]++) /*инициализцая 43 элемента*/
{
if (stl[3][square[43]] && str[4][square[43]])
{
	str[4][square[43]]=0;
	stl[3][square[43]]=0;
	flag[43] = 1;
}
for (square[46] = 0; (flag[43] && (square[46] < 10)); square[46]++) /*инициализцая 56 элемента*/
{
if (str[4][square[46]] && stl[6][square[46]])
{
	str[4][square[46]]=0;
	stl[6][square[46]]=0;
	flag[46] = 1;
}
for (square[47] = 0; (flag[46] && (square[47] < 10)); square[47]++) /*инициализцая 47 элемента*/
{
if (str[4][square[47]] && stl[7][square[47]])
{
	str[4][square[47]]=0;
	stl[7][square[47]]=0;
	flag[47] = 1;
}
for (square[48] = 0; (flag[47] && (square[48] < 10)); square[48]++) /*инициализцая 48 элемента*/
{
if (str[4][square[48]] && stl[8][square[48]])
{
	str[4][square[48]]=0;
	stl[8][square[48]]=0;
	flag[48] = 1;
}
for (square[49] = 0; (flag[48] && (square[49] < 10)); square[49]++) /*инициализцая 49 элемента*/
{
if (str[4][square[49]] && stl[9][square[49]])
{
	str[4][square[49]]=0;
	stl[9][square[49]]=0;
	flag[49] = 1;
}
for (square[50] = 0; (flag[49] && (square[50] < 10)); square[50]++) /*инициализцая 50 элемента*/
{
if (stl[0][square[50]] && str[5][square[50]])
{
	str[5][square[50]]=0;
	stl[0][square[50]]=0;
	flag[50] = 1;
}
for (square[51] = 0; (flag[50] && (square[51] < 10)); square[51]++) /*инициализцая 51 элемента*/
{
if (stl[1][square[51]] && str[5][square[51]])
{
	str[5][square[51]]=0;
	stl[1][square[51]]=0;
	flag[51] = 1;
}
for (square[52] = 0; (flag[51] && (square[52] < 10)); square[52]++) /*инициализцая 52 элемента*/
{
if (stl[2][square[52]] && str[5][square[52]])
{
	str[5][square[52]]=0;
	stl[2][square[52]]=0;
	flag[52] = 1;
}
for (square[53] = 0; (flag[52] && (square[53] < 10)); square[53]++) /*инициализцая 53 элемента*/
{
if (stl[3][square[53]] && str[5][square[53]])
{
	str[5][square[53]]=0;
	stl[3][square[53]]=0;
	flag[53] = 1;
}
for (square[56] = 0; (flag[53] && (square[56] < 10)); square[56]++) /*инициализцая 64 элемента*/
{
if (str[5][square[56]] && stl[6][square[56]])
{
	str[5][square[56]]=0;
	stl[6][square[56]]=0;
	flag[56] = 1;
}
for (square[57] = 0; (flag[56] && (square[57] < 10)); square[57]++) /*инициализцая 65 элемента*/
{
if (str[5][square[57]] && stl[7][square[57]])
{
	str[5][square[57]]=0;
	stl[7][square[57]]=0;
	flag[57] = 1;
}
for (square[58] = 0; (flag[57] && (square[58] < 10)); square[58]++) /*инициализцая 58 элемента*/
{
if (str[5][square[58]] && stl[8][square[58]])
{
	str[5][square[58]]=0;
	stl[8][square[58]]=0;
	flag[58] = 1;
}
for (square[59] = 0; (flag[58] && (square[59] < 10)); square[59]++) /*инициализцая 59 элемента*/
{
if (str[5][square[59]] && stl[9][square[59]])
{
	str[5][square[59]]=0;
	stl[9][square[59]]=0;
	flag[59] = 1;
}
for (square[60] = 0; (flag[59] && (square[60] < 10)); square[60]++) /*инициализцая 60 элемента*/
{
if (stl[0][square[60]] && str[6][square[60]])
{
	str[6][square[60]]=0;
	stl[0][square[60]]=0;
	flag[60] = 1;
}
for (square[61] = 0; (flag[60] && (square[61] < 10)); square[61]++) /*инициализцая 61 элемента*/
{
if (stl[1][square[61]] && str[6][square[61]])
{
	str[6][square[61]]=0;
	stl[1][square[61]]=0;
	flag[61] = 1;
}
for (square[62] = 0; (flag[61] && (square[62] < 10)); square[62]++) /*инициализцая 62 элемента*/
{
if (stl[2][square[62]] && str[6][square[62]])
{
	str[6][square[62]]=0;
	stl[2][square[62]]=0;
	flag[62] = 1;
}
for (square[64] = 0; (flag[62] && (square[64] < 10)); square[64]++) /*инициализцая 71 элемента*/
{
if (stl[4][square[64]] && str[6][square[64]])
{
	str[6][square[64]]=0;
	stl[4][square[64]]=0;
	flag[64] = 1;
}
for (square[65] = 0; (flag[64] && (square[65] < 10)); square[65]++) /*инициализцая 65 элемента*/
{
if (stl[5][square[65]] && str[6][square[65]])
{
	str[6][square[65]]=0;
	stl[5][square[65]]=0;
	flag[65] = 1;
}
for (square[67] = 0; (flag[65] && (square[67] < 10)); square[67]++) /*инициализцая 73 элемента*/
{
if (str[6][square[67]] && stl[7][square[67]])
{
	str[6][square[67]]=0;
	stl[7][square[67]]=0;
	flag[67] = 1;
}
for (square[68] = 0; (flag[67] && (square[68] < 10)); square[68]++) /*инициализцая 68 элемента*/
{
if (str[6][square[68]] && stl[8][square[68]])
{
	str[6][square[68]]=0;
	stl[8][square[68]]=0;
	flag[68] = 1;
}
for (square[69] = 0; (flag[68] && (square[69] < 10)); square[69]++) /*инициализцая 69 элемента*/
{
if (str[6][square[69]] && stl[9][square[69]])
{
	str[6][square[69]]=0;
	stl[9][square[69]]=0;
	flag[69] = 1;
}
for (square[70] = 0; (flag[69] && (square[70] < 10)); square[70]++) /*инициализцая 70 элемента*/
{
if (stl[0][square[70]] && str[7][square[70]])
{
	str[7][square[70]]=0;
	stl[0][square[70]]=0;
	flag[70] = 1;
}
for (square[71] = 0; (flag[70] && (square[71] < 10)); square[71]++) /*инициализцая 71 элемента*/
{
if (stl[1][square[71]] && str[7][square[71]])
{
	str[7][square[71]]=0;
	stl[1][square[71]]=0;
	flag[71] = 1;
}
for (square[73] = 0; (flag[71] && (square[73] < 10)); square[73]++) /*инициализцая 78 элемента*/
{
if (stl[3][square[73]] && str[7][square[73]])
{
	str[7][square[73]]=0;
	stl[3][square[73]]=0;
	flag[73] = 1;
}
for (square[74] = 0; (flag[73] && (square[74] < 10)); square[74]++) /*инициализцая 74 элемента*/
{
if (stl[4][square[74]] && str[7][square[74]])
{
	str[7][square[74]]=0;
	stl[4][square[74]]=0;
	flag[74] = 1;
}
for (square[75] = 0; (flag[74] && (square[75] < 10)); square[75]++) /*инициализцая 75 элемента*/
{
if (stl[5][square[75]] && str[7][square[75]])
{
	str[7][square[75]]=0;
	stl[5][square[75]]=0;
	flag[75] = 1;
}
for (square[76] = 0; (flag[75] && (square[76] < 10)); square[76]++) /*инициализцая 76 элемента*/
{
if (stl[6][square[76]] && str[7][square[76]])
{
	str[7][square[76]]=0;
	stl[6][square[76]]=0;
	flag[76] = 1;
}
for (square[78] = 0; (flag[76] && (square[78] < 10)); square[78]++) /*инициализцая 82 элемента*/
{
if (str[7][square[78]] && stl[8][square[78]])
{
	str[7][square[78]]=0;
	stl[8][square[78]]=0;
	flag[78] = 1;
}
for (square[79] = 0; (flag[78] && (square[79] < 10)); square[79]++) /*инициализцая 79 элемента*/
{
if (str[7][square[79]] && stl[9][square[79]])
{
	str[7][square[79]]=0;
	stl[9][square[79]]=0;
	flag[79] = 1;
}
for (square[80] = 0; (flag[79] && (square[80] < 10)); square[80]++) /*инициализцая 80 элемента*/
{
if (stl[0][square[80]] && str[8][square[80]])
{
	str[8][square[80]]=0;
	stl[0][square[80]]=0;
	flag[80] = 1;
}
for (square[82] = 0; (flag[80] && (square[82] < 10)); square[82]++) /*инициализцая 85 элемента*/
{
if (stl[2][square[82]] && str[8][square[82]])
{
	str[8][square[82]]=0;
	stl[2][square[82]]=0;
	flag[82] = 1;
}
for (square[83] = 0; (flag[82] && (square[83] < 10)); square[83]++) /*инициализцая 83 элемента*/
{
if (stl[3][square[83]] && str[8][square[83]])
{
	str[8][square[83]]=0;
	stl[3][square[83]]=0;
	flag[83] = 1;
}
for (square[84] = 0; (flag[83] && (square[84] < 10)); square[84]++) /*инициализцая 84 элемента*/
{
if (stl[4][square[84]] && str[8][square[84]])
{
	str[8][square[84]]=0;
	stl[4][square[84]]=0;
	flag[84] = 1;
}
for (square[85] = 0; (flag[84] && (square[85] < 10)); square[85]++) /*инициализцая 85 элемента*/
{
if (stl[5][square[85]] && str[8][square[85]])
{
	str[8][square[85]]=0;
	stl[5][square[85]]=0;
	flag[85] = 1;
}
for (square[86] = 0; (flag[85] && (square[86] < 10)); square[86]++) /*инициализцая 86 элемента*/
{
if (stl[6][square[86]] && str[8][square[86]])
{
	str[8][square[86]]=0;
	stl[6][square[86]]=0;
	flag[86] = 1;
}
for (square[87] = 0; (flag[86] && (square[87] < 10)); square[87]++) /*инициализцая 87 элемента*/
{
if (stl[7][square[87]] && str[8][square[87]])
{
	str[8][square[87]]=0;
	stl[7][square[87]]=0;
	flag[87] = 1;
}
for (square[89] = 0; (flag[87] && (square[89] < 10)); square[89]++) /*инициализцая 91 элемента*/
{
if (str[8][square[89]] && stl[9][square[89]])
{
	str[8][square[89]]=0;
	stl[9][square[89]]=0;
	flag[89] = 1;
}
for (square[91] = 0; (flag[89] && (square[91] < 10)); square[91]++) /*инициализцая 92 элемента*/
{
if (stl[1][square[91]] && str[9][square[91]])
{
	str[9][square[91]]=0;
	stl[1][square[91]]=0;
	flag[91] = 1;
}
for (square[92] = 0; (flag[91] && (square[92] < 10)); square[92]++) /*инициализцая 92 элемента*/
{
if (stl[2][square[92]] && str[9][square[92]])
{
	str[9][square[92]]=0;
	stl[2][square[92]]=0;
	flag[92] = 1;
}
for (square[93] = 0; (flag[92] && (square[93] < 10)); square[93]++) /*инициализцая 93 элемента*/
{
if (stl[3][square[93]] && str[9][square[93]])
{
	str[9][square[93]]=0;
	stl[3][square[93]]=0;
	flag[93] = 1;
}
for (square[94] = 0; (flag[93] && (square[94] < 10)); square[94]++) /*инициализцая 94 элемента*/
{
if (stl[4][square[94]] && str[9][square[94]])
{
	str[9][square[94]]=0;
	stl[4][square[94]]=0;
	flag[94] = 1;
}
for (square[95] = 0; (flag[94] && (square[95] < 10)); square[95]++) /*инициализцая 95 элемента*/
{
if (stl[5][square[95]] && str[9][square[95]])
{
	str[9][square[95]]=0;
	stl[5][square[95]]=0;
	flag[95] = 1;
}
for (square[96] = 0; (flag[95] && (square[96] < 10)); square[96]++) /*инициализцая 96 элемента*/
{
if (stl[6][square[96]] && str[9][square[96]])
{
	str[9][square[96]]=0;
	stl[6][square[96]]=0;
	flag[96] = 1;
}
for (square[97] = 0; (flag[96] && (square[97] < 10)); square[97]++) /*инициализцая 97 элемента*/
{ 
if (stl[7][square[97]] && str[9][square[97]])
{
	str[9][square[97]]=0;
	stl[7][square[97]]=0;
	flag[97] = 1;
}
for (square[98] = 0; (flag[97] && (square[98] < 10)); square[98]++) /*инициализцая 98 элемента*/
{
if (stl[8][square[98]] && str[9][square[98]])
{
	str[9][square[98]]=0;
	stl[8][square[98]]=0;
	flag[98] = 1;
}

if(flag[98])
{
	/* ДЛК сгенерирован*/
count++;

// тут расположить код проверки коэфициента ортогональности
result=твоя_функция_проверки();

if(result>local_max)
{
	local_max=result;
}

if(count==NUM_OF_DLK_IN_ONE_CHECK)
{
	//получить глобальный рекорд
	global_max=твоя_функция_получения_глобального_максимума_ортогональности_среди_всех_областей();
	
	if(local_max>global_max)	
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
	}
}

}


if(flag[98])
{
	str[9][square[98]]=1;
	stl[8][square[98]]=1;
	flag[98] = 0;
}
}

if(flag[97])
{
	str[9][square[97]]=1;
	stl[7][square[97]]=1;
	flag[97] = 0;
}
}
if(flag[96])
{
	str[9][square[96]]=1;
	stl[6][square[96]]=1;
	flag[96] = 0;
}
}
if(flag[95])
{
	str[9][square[95]]=1;
	stl[5][square[95]]=1;
	flag[95] = 0;
}
}
if(flag[94])
{
	str[9][square[94]]=1;
	stl[4][square[94]]=1;
	flag[94] = 0;
}
}
if(flag[93])
{
	str[9][square[93]]=1;
	stl[3][square[93]]=1;
	flag[93] = 0;
}
}
if(flag[92])
{
	str[9][square[92]]=1;
	stl[2][square[92]]=1;
	flag[92] = 0;
}
}
if(flag[91])
{
	str[9][square[91]]=1;
	stl[1][square[91]]=1;
	flag[91] = 0;
}
}
if(flag[89])
{
	str[8][square[89]]=1;
	stl[9][square[89]]=1;
	flag[89] = 0;
}
}
if(flag[87])
{
	str[8][square[87]]=1;
	stl[7][square[87]]=1;
	flag[87] = 0;
}
}
if(flag[86])
{
	str[8][square[86]]=1;
	stl[6][square[86]]=1;
	flag[86] = 0;
}
}
if(flag[85])
{
	str[8][square[85]]=1;
	stl[5][square[85]]=1;
	flag[85] = 0;
}
}
if(flag[84])
{
	str[8][square[84]]=1;
	stl[4][square[84]]=1;
	flag[84] = 0;
}
}
if(flag[83])
{
	str[8][square[83]]=1;
	stl[3][square[83]]=1;
	flag[83] = 0;
}
}
if(flag[82])
{
	str[8][square[82]]=1;
	stl[2][square[82]]=1;
	flag[82] = 0;
}
}
if(flag[80])
{
	str[8][square[80]]=1;
	stl[0][square[80]]=1;
	flag[80] = 0;
}
}
if(flag[79])
{
	str[7][square[79]]=1;
	stl[9][square[79]]=1;
	flag[79] = 0;
}
}
if(flag[78])
{
	str[7][square[78]]=1;
	stl[8][square[78]]=1;
	flag[78] = 0;
}
}
if(flag[76])
{
	str[7][square[76]]=1;
	stl[6][square[76]]=1;
	flag[76] = 0;
}
}
if(flag[75])
{
	str[7][square[75]]=1;
	stl[5][square[75]]=1;
	flag[75] = 0;
}
}
if(flag[74])
{
	str[7][square[74]]=1;
	stl[4][square[74]]=1;
	flag[74] = 0;
}
}
if(flag[73])
{
	str[7][square[73]]=1;
	stl[3][square[73]]=1;
	flag[73] = 0;
}
}
if(flag[71])
{
	str[7][square[71]]=1;
	stl[1][square[71]]=1;
	flag[71] = 0;
}
}
if(flag[70])
{
	str[7][square[70]]=1;
	stl[0][square[70]]=1;
	flag[70] = 0;
}
}
if(flag[69])
{
	str[6][square[69]]=1;
	stl[9][square[69]]=1;
	flag[69] = 0;
}
}
if(flag[68])
{
	str[6][square[68]]=1;
	stl[8][square[68]]=1;
	flag[68] = 0;
}
}
if(flag[67])
{
	str[6][square[67]]=1;
	stl[7][square[67]]=1;
	flag[67] = 0;
}
}
if(flag[65])
{
	str[6][square[65]]=1;
	stl[5][square[65]]=1;
	flag[65] = 0;
}
}
if(flag[64])
{
	str[6][square[64]]=1;
	stl[4][square[64]]=1;
	flag[64] = 0;
}
}
if(flag[62])
{
	str[6][square[62]]=1;
	stl[2][square[62]]=1;
	flag[62] = 0;
}
}
if(flag[61])
{
	str[6][square[61]]=1;
	stl[1][square[61]]=1;
	flag[61] = 0;
}
}
if(flag[60])
{
	str[6][square[60]]=1;
	stl[0][square[60]]=1;
	flag[60] = 0;
}
}
if(flag[59])
{
	str[5][square[59]]=1;
	stl[9][square[59]]=1;
	flag[59] = 0;
}
}
if(flag[58])
{
	str[5][square[58]]=1;
	stl[8][square[58]]=1;
	flag[58] = 0;
}
}
if(flag[57])
{
	str[5][square[57]]=1;
	stl[7][square[57]]=1;
	flag[57] = 0;
}
}
if(flag[56])
{
	str[5][square[56]]=1;
	stl[6][square[56]]=1;
	flag[56] = 0;
}
}
if(flag[53])
{
	str[5][square[53]]=1;
	stl[3][square[53]]=1;
	flag[53] = 0;
}
}
if(flag[52])
{
	str[5][square[52]]=1;
	stl[2][square[52]]=1;
	flag[52] = 0;
}
}
if(flag[51])
{
	str[5][square[51]]=1;
	stl[1][square[51]]=1;
	flag[51] = 0;
}
}
if(flag[50])
{
	str[5][square[50]]=1;
	stl[0][square[50]]=1;
	flag[50] = 0;
}
}
if(flag[49])
{
	str[4][square[49]]=1;
	stl[9][square[49]]=1;
	flag[49] = 0;
}
}
if(flag[48])
{
	str[4][square[48]]=1;
	stl[8][square[48]]=1;
	flag[48] = 0;
}
}
if(flag[47])
{
	str[4][square[47]]=1;
	stl[7][square[47]]=1;
	flag[47] = 0;
}
}
if(flag[46])
{
	str[4][square[46]]=1;
	stl[6][square[46]]=1;
	flag[46] = 0;
}
}
if(flag[43])
{
	str[4][square[43]]=1;
	stl[3][square[43]]=1;
	flag[43] = 0;
}
}
if(flag[42])
{
	str[4][square[42]]=1;
	stl[2][square[42]]=1;
	flag[42] = 0;
}
}
if(flag[41])
{
	str[4][square[41]]=1;
	stl[1][square[41]]=1;
	flag[41] = 0;
}
}
if(flag[40])
{
	str[4][square[40]]=1;
	stl[0][square[40]]=1;
	flag[40] = 0;
}
}
if(flag[39])
{
	str[3][square[39]]=1;
	stl[9][square[39]]=1;
	flag[39] = 0;
}
}
if(flag[38])
{
	str[3][square[38]]=1;
	stl[8][square[38]]=1;
	flag[38] = 0;
}
}
if(flag[37])
{
	str[3][square[37]]=1;
	stl[7][square[37]]=1;
	flag[37] = 0;
}
}
if(flag[35])
{
	str[3][square[35]]=1;
	stl[5][square[35]]=1;
	flag[35] = 0;
}
}
if(flag[34])
{
	str[3][square[34]]=1;
	stl[4][square[34]]=1;
	flag[34] = 0;
}
}
if(flag[32])
{
	str[3][square[32]]=1;
	stl[2][square[32]]=1;
	flag[32] = 0;
}
}
if(flag[31])
{
	str[3][square[31]]=1;
	stl[1][square[31]]=1;
	flag[31] = 0;
}
}
if(flag[30])
{
	str[3][square[30]]=1;
	stl[0][square[30]]=1;
	flag[30] = 0;
}
}
if(flag[29])
{
	str[2][square[29]]=1;
	stl[9][square[29]]=1;
	flag[29] = 0;
}
}
if(flag[28])
{
	str[2][square[28]]=1;
	stl[8][square[28]]=1;
	flag[28] = 0;
}
}
if(flag[26])
{
	str[2][square[26]]=1;
	stl[6][square[26]]=1;
	flag[26] = 0;
}
}
if(flag[25])
{
	str[2][square[25]]=1;
	stl[5][square[25]]=1;
	flag[25] = 0;
}
}
if(flag[24])
{
	str[2][square[24]]=1;
	stl[4][square[24]]=1;
	flag[24] = 0;
}
}
if(flag[23])
{
	str[2][square[23]]=1;
	stl[3][square[23]]=1;
	flag[23] = 0;
}
}
if(flag[21])
{
	str[2][square[21]]=1;
	stl[1][square[21]]=1;
	flag[21] = 0;
}
}
if(flag[20])
{
	str[2][square[20]]=1;
	stl[0][square[20]]=1;
	flag[20] = 0;
}
}
if(flag[19])
{
	str[1][square[19]]=1;
	stl[9][square[19]]=1;
	flag[19] = 0;
}
}
if(flag[17])
{
	str[1][square[17]]=1;
	stl[7][square[17]]=1;
	flag[17] = 0;
}
}
if(flag[16])
{
	str[1][square[16]]=1;
	stl[6][square[16]]=1;
	flag[16] = 0;
}
}
if(flag[15])
{
	str[1][square[15]]=1;
	stl[5][square[15]]=1;
	flag[15] = 0;
}
}
if(flag[14])
{
	str[1][square[14]]=1;
	stl[4][square[14]]=1;
	flag[14] = 0;
}
}
if(flag[13])
{
	str[1][square[13]]=1;
	stl[3][square[13]]=1;
	flag[13] = 0;
}
}
if(flag[12])
{
	str[1][square[12]]=1;
	stl[2][square[12]]=1;
	flag[12] = 0;
}
}
if(flag[10])
{
	str[1][square[10]]=1;
	stl[0][square[10]]=1;
	flag[10] = 0;
}
}
if(flag[99])
{
	str[9][square[99]]=1;
	stl[9][square[99]]=1;
	diag[0][square[99]]=1;
	flag[99] = 0;
}
}
if(flag[90])
{
	str[9][square[90]]=1;
	stl[0][square[90]]=1;
	diag[1][square[90]]=1;
	flag[90] = 0;
}
}
if(flag[88])
{
	str[8][square[88]]=1;
	stl[8][square[88]]=1;
	diag[0][square[88]]=1;
	flag[88] = 0;
}
}
if(flag[81])
{
	str[8][square[81]]=1;
	stl[1][square[81]]=1;
	diag[1][square[81]]=1;
	flag[81] = 0;
}
}
if(flag[77])
{
	str[7][square[77]]=1;
	stl[7][square[77]]=1;
	diag[0][square[77]]=1;
	flag[77] = 0;
}
}
if(flag[72])
{
	str[7][square[72]]=1;
	stl[2][square[72]]=1;
	diag[1][square[72]]=1;
	flag[72] = 0;
}
}
if(flag[66])
{
	str[6][square[66]]=1;
	stl[6][square[66]]=1;
	diag[0][square[66]]=1;
	flag[66] = 0;
}
}
if(flag[63])
{
	str[6][square[63]]=1;
	stl[3][square[63]]=1;
	diag[1][square[63]]=1;
	flag[63] = 0;
}
}
if(flag[55])
{
	str[5][square[55]]=1;
	stl[5][square[55]]=1;
	diag[0][square[55]]=1;
	flag[55] = 0;
}
}
if(flag[54])
{
	str[5][square[54]]=1;
	stl[4][square[54]]=1;
	diag[1][square[54]]=1;
	flag[54] = 0;
}
}
if(flag[45])
{
	str[4][square[45]]=1;
	stl[5][square[45]]=1;
	diag[1][square[45]]=1;
	flag[45] = 0;
}
}
if(flag[44])
{
	str[4][square[44]]=1;
	stl[4][square[44]]=1;
	diag[0][square[44]]=1;
	flag[44] = 0;
}
}
if(flag[36])
{
	str[3][square[36]]=1;
	stl[6][square[36]]=1;
	diag[1][square[36]]=1;
	flag[36] = 0;
}
}
if(flag[33])
{
	str[3][square[33]]=1;
	stl[3][square[33]]=1;
	diag[0][square[33]]=1;
	flag[33] = 0;
}
}
if(flag[27])
{
	str[2][square[27]]=1;
	stl[7][square[27]]=1;
	diag[1][square[27]]=1;
	flag[27] = 0;
}
}
if(flag[22])
{
	str[2][square[22]]=1;
	stl[2][square[22]]=1;
	diag[0][square[22]]=1;
	flag[22] = 0;
}
}
if(flag[18])
{
	str[1][square[18]]=1;
	stl[8][square[18]]=1;
	diag[1][square[18]]=1;
	flag[18] = 0;
}
}
if(flag[11])
{
	str[1][square[11]]=1;
	stl[1][square[11]]=1;
	diag[0][square[11]]=1;
	flag[11] = 0;
}
}
if(flag[9])
{
	str[0][square[9]]=1;
	stl[9][square[9]]=1;
	diag[1][square[9]]=1;
	flag[9] = 0;
}
}
if(flag[8])
{
	str[0][square[8]]=1;
	stl[8][square[8]]=1;
	flag[8] = 0;
}
}
if(flag[7])
{
	str[0][square[7]]=1;
	stl[7][square[7]]=1;
	flag[7] = 0;
}
}
if(flag[6])
{
	str[0][square[6]]=1;
	stl[6][square[6]]=1;
	flag[6] = 0;
}
}
if(flag[5])
{
	str[0][square[5]]=1;
	stl[5][square[5]]=1;
	flag[5] = 0;
}
}
if(flag[4])
{
	str[0][square[4]]=1;
	stl[4][square[4]]=1;
	flag[4] = 0;
}
}
if(flag[3])
{
	str[0][square[3]]=1;
	stl[3][square[3]]=1;
	flag[3] = 0;
}
}
if(flag[2])
{
	str[0][square[2]]=1;
	stl[2][square[2]]=1;
	flag[2] = 0;
}
}
if(flag[1])
{
	str[0][square[1]]=1;
	stl[1][square[1]]=1;
	flag[1] = 0;
}
}
if(flag[0])
{
	str[0][square[0]]=1;
	stl[0][square[0]]=1;
	diag[0][square[0]]=1;
	flag[0] = 0;
}
}

return (local_max);
}


int main()
{
	int i=0;
	generate_dlk(NUM_OF_PARTS,i);
}