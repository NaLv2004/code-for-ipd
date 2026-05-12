#pragma once
void listSphereDecode(int* u, float* y, int** G, int* A, int N, int L, int* n_paths);
void listSphereDecode_step(int* u, float* y, int** G, int* A, int N, int L, int* n_paths);
void LSD_decode_outall(int* u, int** paths_extend, float* y, int** G, int* A, int N, int L, int* n_paths);
bool judgeSame(int t, int* u, int* x, int** G, int N);
bool judgeSame_append(int t, int* u, int* x, int** G, int N, int bit_append);
void calc_npaths(int* A, int N, int* n_paths, int L);
void mySort(float* arr, int* idx_sort, int N);