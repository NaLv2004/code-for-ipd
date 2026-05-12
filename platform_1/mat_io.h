#pragma once
#include <iostream>
#include <Eigen/Dense>
#include <complex>
#include "setting.h"
#include <fstream>
void save_mat_txt(const Eigen::MatrixXf& mat, const std::string& file_name);
void save_array_txt(float SNR, float BER, float FER, const std::string& file_name);
void save_array_txt(int* arr, int len, const std::string& file_name);
void save_array_txt_commasep(int* arr, int len, const std::string& file_name);
void save_mat_txt(const float* mat, int cols, const std::string& file_name);
void display_array(int** arr, int rows, int cols);
void save_array_txt(float SNR, float FER, const std::string& file_name);