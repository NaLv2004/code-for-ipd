#include <vector>
#include <xmmintrin.h>
#include <immintrin.h>
#include <smmintrin.h>
#include "setting.h"

// void PolarDecodeIterVertical(float** oriLLRh, float** upLLR, int** umid, int* A_Acv, int* Av, int Nh, int Nv, int decodingi, int iter);
void PolarDecodeIterVertical(float** oriLLRh, float** upLLR, int** umid, int** G, int** uout_partial, int* A_Acv, int* Av, int Nh, int Nv, int decodingi, int iter);
int judge_rate(int* A_Ac, int len);
int judge_rate_L(int* A_Ac, int len);
void MITSOh(float** LLR, float* PM, int** sum, float** oriLLRh, int* indexpm, int iter, int decodingi, float** upLLR, float Qsumun);
void MITSOh_irregular(float** LLR, float* PM, int** sum, float** oriLLRh, int* indexpm, int iter, int decodingi, float** upLLR, float Qsumun, int* irregular_construction);
void MITSOv(float** LLR, float* PM, int** sum, float** oriLLRh, int* indexpm, int iter, int decodingi, float** upLLR, float Qsumun);
void Pyndiahh(float** LLR, int** sum, float** oriLLRh, int* indexpm, int iter, int decodingi, float** upLLR);
void Pyndiahh_getsum(float** LLR, int** sum, float** oriLLRh, int* indexpm, int iter, int decodingi, float** upLLR, double& sum_sdis);
void Pyndiahh_irregular(float** LLR, int** sum, float** oriLLRh, int* indexpm, int iter, int decodingi, float** upLLR, double& sum_sdis, int* irregular_construction);
void Pyndiahv(float** LLR, int** sum, float** oriLLRh, int* indexpm, int iter, int decodingi, float** upLLR);
void Pyndiahh_LSD(float** LLR, int** sum, float** oriLLRh, int* indexpm, int iter, int decodingi, float** upLLR);
void Pyndiahh_LSD_getsum(float** LLR, int** sum, float** oriLLRh, int* indexpm, int iter, int decodingi, float** upLLR, double& sdis_sum, int** GGh, int* Ach);
void Pyndiahv_adaptivebeta(float** LLR, int** sum, float** oriLLRh, int* indexpm, int iter, int decodingi, float** upLLR);
void Pyndiahh_adaptivebeta(float** LLR, int** sum, float** oriLLRh, int* indexpm, int iter, int decodingi, float** upLLR);
void Pyndiahh_LSD_adaptivebeta(float** LLR, int** sum, float** oriLLRh, int* indexpm, int iter, int decodingi, float** upLLR);
void Pyndiahh_bitflip(float** LLR, int** sum, float** oriLLRh, int* indexpm, int iter, int decodingi, float** upLLR, double& sdis_sum, int** tried_codewords_h, int** GGh, int* Ach, int* n_tried_patterns_h);
void Pyndiahv_bitflip(float** LLR, int** sum, float** oriLLRh, int* indexpm, int iter, int decodingi, float** upLLR, double& sdis_sum, int** tried_codewords_v, int** GGv, int* Acv, int* n_tried_patterns_v);
void my_sort_pyndiah(int* arr, int* idx_sort, int N);
void my_sort_pyndiah(float* arr, int* idx_sort, int N);
void f_function(float* a, int size);
void f_function_index(float* a, float* a_new, int size);
void g_function(float* a, int* b, int size);
void g_function_index(float* a, int* b, float* a_new, int size);
void combine(int* a, int* b, int size);
void combine_index(int* a, int* b, int* a_new, int size);
void hard_SIMD(float* a, int* b, int size);
float cal_sum(float* a, int len);
void replace_sum(int* a, int* b, int len);
void replace_LLR(float* a, float* b, int len);
void set_PM(float** a, float* b, int L, int len);
void stage_located(int L, int* index, float** LLR, int** sum, int** p, int Count, int N, vector<int> path_2, vector<int> path_0);
