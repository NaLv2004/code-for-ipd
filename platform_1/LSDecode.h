
#include <stdio.h>
#include <stdlib.h> 
#include <iostream>
#include <cmath>
#include <bitset>
#include <algorithm>
#include <complex>
#include <Eigen\Dense>
#include <Eigen\Core>
#include <unsupported\Eigen\KroneckerProduct>
#include <fstream>
using namespace std;

// int calcd(int CodeLength, int i, int* u, int **G);
int calcd(int* A, int K, int i, int* u, int** G, int* t);
int calcd2(int* A, int K, int i, int* u, int** G, int* t);
void calcdistance(float *y, float **d, int N);
int partition(float a[], int s, int t);
float findmink(float a[], int low, int high, int k);
void calbestworst(float** d, int K, float* best, float* worst, vector <int>* sinset);
// void LSDecode(float *y, int* A_Ac, int N, int K, int L, int **u, int **GG,int breakpoint,float *d, clock_t *T_start, clock_t *T_finish);
void LSDecode(int* A, float* y, int* A_Ac, int N, int K, int L, int** u, int** GG, vector <int>* sinset, int* DF, int** GGF, int breakpoint, float* d, clock_t* T_start, clock_t* T_finish, int* timesnumber);