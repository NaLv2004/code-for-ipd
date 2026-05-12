#pragma once
#include <iostream>
#include <Eigen/Dense>
#include <complex>
#include "setting.h"
#include <vector>
#include <algorithm>
#include <omp.h>
using namespace Eigen;
// 假设 sampleNormal() 函数的定义
double sampleNormal_mimo();
MatrixXcf generateChannelMatrix(int Nr, int Nt);
void AWGN(int Nr, int Nt, std::complex<float>* x, std::complex<float>* y, const MatrixXcf& H, double n0);
// MatrixXf convertHToReal(const MatrixXcf& H);
void convertHToReal(const MatrixXcf& H, Matrix<float, 2 * Nr, 2 * Nt>& H_real);
void convertHToReal(const MatrixXcf& H, MatrixXf& H_real);
void convertYToReal(int Nr, std::complex<float>* y, float* y_real);
//void MMSEdetection(int Nr, int Nt, const MatrixXf& H_real, float* rx, float* y_hat, float sigma_n2);
//void MMSEdetection(const MatrixXf& H_real, const MatrixXf& identity, MatrixXf& received_signal, MatrixXf& mmse_matrix, MatrixXf& output_signal, MatrixXf& HTH, float* rx, float* y_hat, float sigma_n2);
//void MMSEdetection(const Matrix<double, 2 * Nr, 2 * Nt>& H_real, const Matrix<double, 2 * Nt, 2 * Nt>& identity, Matrix<float, 2 * Nr, 1>& received_signal, Matrix<float, 2 * Nt, 2 * Nr>& mmse_matrix, Matrix<float, 2 * Nt, 1>& output_signal, Matrix<float, 2 * Nt, 2 * Nt>& HTH, float* rx, float* y_hat, float sigma_n2);
void MMSEdetection(const Matrix<float, 2 * Nr, 2 * Nt>& H_real, const Matrix<float, 2 * Nt, 2 * Nt>& identity, Matrix<float, 2 * Nr, 1>& received_signal, Matrix<float, 2 * Nt, 2 * Nr>& mmse_matrix, Matrix<float, 2 * Nt, 1>& output_signal, Matrix<float, 2 * Nt, 2 * Nt>& HTH, float* rx, float* y_hat, float sigma_n2);
void harddec_bits(int L, float* llr_bits, int* dec_bits);
//void MMSEdetection(const MatrixXf H_real, const MatrixXf& identity, MatrixXf& received_signal, MatrixXf& mmse_matrix, MatrixXf& output_signal, MatrixXf& HTH, float* rx, float* y_hat, float sigma_n2);
void MMSEdetection(const MatrixXf H_real, const MatrixXf identity, MatrixXf received_signal, MatrixXf mmse_matrix, MatrixXf output_signal, MatrixXf HTH, float* rx, float* y_hat, float sigma_n2);
void generateIntrlvPattern(std::vector<int>& pattern, int len);
void generateIntrlvPatternInfoMapping(std::vector<int>& pattern, int nv, int kv, int* Av, int* Acv);
// void generateIntrlvPatternIrregular(std::vector<int>& pattern, int len, int* bad_channel_index, int* good_channel_index);
void generateIntrlvPatternIrregular(std::vector<int>& pattern, int len, int* bad_channel_index, int* good_channel_index, int cnt_bad_channels, int* use_short_codeword);
void generateIntrlvPatternIrregularArray(std::vector<int>& pattern, int len, int* bad_channel_index, int* good_channel_index, int cnt_bad_channels, int* use_short_codeword);
void randIntrlv(std::vector<int>& pattern, int* t_bits, int len);
void randDeintrlv(std::vector<int>& pattern, float* rx_llr, int len);
void randIntrlv(std::vector<int>& pattern, float* t_bits, int len);
void randIntrlvBlock(int** bits, int len_block, int height_block, int len_subblock, int interleave_rows);
void randDeIntrlvBlock(float** bits, int len_block, int height_block, int len_subblock, int interleave_rows);
void randIntrlvPartial(std::vector<int>& pattern, int* bits, int* info_pos, int N, int K_info);
void randIntrlvPartial(std::vector<int>& pattern, float* bits, int* info_pos, int N, int K_info);
void randDeIntrlvPartial(std::vector<int>& pattern, float* bits, int* info_pos, int N, int K_info);
struct Candidate;
std::vector<std::complex<double>> KBest(int k,
    const Eigen::MatrixXcd& R,
    const Eigen::VectorXcd& z_tilde,
    const std::vector<std::complex<double>>& sym);
// void KBestDetection(MatrixXf R, MatrixXf z, float* y_out, int k, int mod_type, int Nt);
void my_sort_kbest(float* arr, int* idx_sort, int N);
void KBestDetection(MatrixXf R, MatrixXf z, MatrixXf H, float* received_signal, float* ED, float* y_out, float** paths_extend, int k, int mod_type, int Nt, int Nr);
void KBestDetection_soft_JDD(MatrixXf R, MatrixXf z, MatrixXf H,
    float* y, float* ED, float* y_out, float* LLR_in, float** paths_extend,
    int k, int mod_type, int Nt, int Nr, int decoding_stage);
//void MMSEdetection(const MatrixXf& H_real, const MatrixXf& identity, MatrixXf& received_signal, MatrixXf& mmse_matrix, MatrixXf& output_signal, MatrixXf& HTH, float* rx, float* y_hat, float sigma_n2);
void MMSEdetection_noInverse(const MatrixXf H_real, const MatrixXf identity, MatrixXf received_signal, MatrixXf mmse_matrix, MatrixXf output_signal, MatrixXf HTH, MatrixXf HTHIinv, float* rx, float* y_hat, float sigma_n2);
void MMSEdetection_exact(const MatrixXf H_real, const MatrixXf identity, MatrixXf received_signal, MatrixXf mmse_matrix, MatrixXf output_signal, MatrixXf HTH, float* rx, float* y_hat, float* mu_mmse, float sigma_n2);
void EPD_detection_float(int flag, size_t Nt, size_t Nr, size_t iter, size_t mod_type,
    float es, float N0, float* HTH, float* HTY, int temp,
    float* extr_mean, float* extr_var,
    double* vec_alpha_in, double* vec_gamma_in, double* p_symbol_rearrange);
MatrixXcf generateChannelMatrixNonuniform(int Nr, int Nt, float** channel_scaling_factor);
void getChannelScalingFactor(int Nr, int Nt, int diff_gain, int* channel_index, float** channel_scale_factor, int predefined);
void precodeSVD(int Nr, int Nt, std::complex<float>* x, std::complex<float>* x_precode_out, const MatrixXcf& V);
void precodeSVD_preprocessing(int Nr, int Nt, std::complex<float>* y, std::complex<float>* y_preprocess_out, const MatrixXcf& U);
