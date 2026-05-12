#include "Polar_Encoder.h"
#include <xmmintrin.h>
#include <immintrin.h>
#include <smmintrin.h>

void decode(float* LLR, int* sum, int j, int* Ac, int a, int N, int K, int* Count);
//void print(int* LLR, int len);
void decode_list(float** LLR, int** sum, int j, int* Ac, int a, int N, int K, float* PM, int L, int m,
	float* LLR_in, float* W, int* index, int* better, int* worse, int* Count, int* Count_info, float* Qsumun);
//void decode_list(int** u1, float** LLR, int** sum, int j, int* Ac, int a, int N, int K, float* PM);

int cal_stage(int i, int n);