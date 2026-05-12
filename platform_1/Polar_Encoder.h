#include <stdio.h>
#include <stdlib.h> 
#include <time.h> 
#include <iostream>
#include <cmath>
#include <vector>
#include <bitset>
#include <algorithm>
#include <xmmintrin.h>
#include <immintrin.h>
#include <smmintrin.h>
using namespace std;

void RM_codeconstruction(int N, int** GG, vector<int>& best_channel);
void B5G_codeconstruction(int N, int K, vector<int>& best_channel);
void GA_codeconstruction(int CodeLength, float sigma, vector<int>& best_channel);
void RM_GA_codeconstruction(int N, int K, int** GG, float sigma, vector<int>& best_channel);
void beta_codeconstruction(int N, vector<int>& best_channel);
float MyRecursiveFun(int i, int NN, float p);//����I(W)
//int* bitreorder(int* min, int len);
void PolarEncode_xor(int* uout, int* uin, int len);
void PolarEncodePartial(int** G, int* u_in, int* u_out, int N, int u_in_len);
double sampleNormal();
