#pragma once
# include<complex>
# include<stdio.h>
# include<iostream>
void modulation(int* input_bits, std::complex<float>* symbols_out, int mod_type, int nof_symbols);
void llr_mmse_bit(int Nt, float* u_sym_mmse, float* llr_bit, float sigma_sq, int ModType);
void llr_mmse_bit_exact(int Nt, float* mu_mmse, float* u_sym_mmse, float* llr_bit, float sigma_sq, int ModType);
void llr_kbest_bit(int Nt, float* u_sym_kbest, float* llr_bit, float sigma_sq, int ModType, int k_kbest, float** paths_kbest, float* ED);
void llr_ep_bit(int Nt, float* ep_extr_mean, float* ep_extr_var, float* llr_bit, int ModType);
void my_sort_ed(float* arr, int* idx_sort, int N);