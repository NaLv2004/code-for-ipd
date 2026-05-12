#pragma once
// void crc_list_extend(int** sub_codeword_list, float** edistance, int L, int Nv, int Nh);
void crc_list_extend(int** sub_codeword_list, int** whole_codeword_list, float** pm, int L, int Nv, int Nh, int L_large);
// void crc_list_extend_vertical(int** sub_codeword_list, int** whole_codeword_list, float** pm, int L, int Nv, int Nh, int L_large);
void crc_list_extend_vertical(int** sub_codeword_list, int** whole_codeword_list, vector<int>** intrlv_pattern, float** pm, float** y, float* H_real, int L, int Nv, int Nh, int Nr, int L_large, int ModType);
float calcChannelEuclidianDistance(int* x, float* y, float* H_real, std::vector<int>** intrlv_pattern, int idx_sym_col, int ModType, int Nh, int Nv, int Nr);
float calculate_norm(const float* y, const float* H, const float* x, int Nr, int Nt);
void crc_list_extend_vertical_singlestep(int** sub_codeword_list, int** whole_codeword_list, vector<int>** intrlv_pattern, float** pm, float** y, float* H_real,
  float* pm_renew, int L, int Nv, int Nh, int Nr, int L_mimo_blk, int ModType, int idx_sym_col);





void crc_list_extend_vertical_doublestage(int** sub_codeword_list, int** whole_codeword_list, vector<int>** intrlv_pattern, float** pm,
    float** y, float* H_real, int L, int Nv, int Nh, int Nr, int L_large, int L_mimo_blk, int ModType);