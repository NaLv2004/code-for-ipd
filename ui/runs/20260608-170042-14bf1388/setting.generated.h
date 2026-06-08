#pragma once
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
using namespace std;
//Coding Setting
const int dimensions = 2; // 1D or 2D
//1D
//const int N = 128;
//const int K = 64;    // pure info. bits
// 4096,3302

//const int N = 4096; // 628
//const int K = 3246;
//const string polarcon = "GA";  // RM/5G/beta/GA
//const int Nh = 64;
//const int Nv = 64;
//const int Kh = 8;
//const int Kv = 57;
//const int Nv_1D = 64;

//const int N = 2048; // 628
//const int K = 1635;
//const string polarcon = "GA";  // RM/5G/beta/GA
//const int Nh = 32;
//const int Nv = 64;
//const int Kh = 8;
//const int Kv = 57;
//const int Nv_1D = 32;

//const int N = 1024; // 628
//const int K = 823;
//const string polarcon = "GA";  // RM/5G/beta/GA
//const int Nh = 16;
//const int Nv = 64;
//const int Kh = 8;
//const int Kv = 57;
//const int Nv_1D = 16;
//2D
//const int Nh = 64;
//const int Nv = 64;
//const int Kh = 42;
//const int Kv = 42;

//const int Nh = 16;
//const int Nv = 32;
//const int Kh = 6;
//const int Kv = 26;

//const int Nh = 8;
//const int Nv = 32;
//const int Kh = 4;
//const int Kv = 4;

//const int N = 2048; // 628
//const int K = 1363;
//const string polarcon = "RM";  // RM/5G/beta/GA
//const int Nh = 64;
//const int Nv = 64;
//const int Kh = 57;
//const int Kv = 57;
//const int Nv_1D = 64;

//const int N = 512; // 628
//const int K = 417;
//const string polarcon = "GA";  // RM/5G/beta/GA
//const int Nh = 8;
//const int Nv = 64;
//const int Kh = 1;
//const int Kv = 1;
//const int Nv_1D = Nh;

//const int N = 256; // 628
//const int K = 214;
//const string polarcon = "GA";  // RM/5G/beta/GA
//const int Nh = 4;
//const int Nv = 64;
//const int Kh = 1;
//const int Kv = 1;
//const int Nv_1D = Nh;

//const int N = 1024; // 628
//const int K = 676;
//const string polarcon = "RM";  // RM/5G/beta/GA
//const int Nh = 64;
//const int Nv = 64;
//const int Kh = 57;
//const int Kv = 57;
//const int Nv_1D = Nh;

const int N = 1024; // 628
const int K = 567;
const string polarcon = "RM"; // RM/5G/beta/GA
const int Nh = 32;
const int Nv = 32;
const int Kh = 22;
const int Kv = 26;
const int Nv_1D = Nh;

const string polarconh = "RM-GA"; // RM/5G/beta/GA/RM-GA
const string polarconv = "RM"; // RM/5G/beta/5
//const int vertical_construction_fixed[Nv] = { 62,61,59,55,47,60,58,57,54,53,51,46,45,43,30,29,39,27,23,56,52,15,50,44,49,42,41,28,38,26,37,25,35,22,21,14,19,48,13,40,11,36,24,34,7,20,33,18,12,17,10,9,6,32,5,3,16,8,4,2,1,0,63,31 };
const bool sys = false;
const int softiter = 2; // prev: 3  prev (list-merge): 3
//MIMO 
const int Nt = 32; // 88  //86
const int Nr = 64; // 64   //16
const int trans_stream = 8;
const int ModType = 4;
const string MOD_DIM = "Temporal"; // Spatial or Temporal Modulation
const string MIMODET = "MMSE"; // "KBEST" or "MMSE" have been implemented. "KBEST" has only been implemented for 16QAM, spatial modulation.
// const int K_KBEST = 64;
const int K_KBEST = 256;
const int iter_MIMO = 7; // n of EP detector iterations
const bool flag_interleave = true; // whether interleave spatially
const bool flag_interleave_all = false;
const bool flag_interleave_partial = true; // true
const bool flag_interleave_spatial = true; // prev:true(1027,1103)
const bool flag_interleave_block = false;
const bool flag_interleave_inner = false;
const int interleave_rows = 64;
const float mimo_llr_thres = 0.3;
//const float alpha[8] = { 0,0.2,0.3,0.5,0.7,0.9,1,1 };
//const string softmethod = "Pyndiah"; //Pyndiah/MITSO
//Soft Iterations
const string softmethod = "Pyndiah"; // Pyndiah/MITSO
const float SOalpha = 0.5; // MITSO //prev 0.5
const int alphalen = 10;
// const float alpha[alphalen] = { 0.30, 0.30, 0.45, 0.50, 0.60, 0.65,0.50,0.60,0.8,0.85 }; // in use for scl
// const float alpha[alphalen] = {0, 0.30, 0.30, 0.450, 0.50,0.60, 0.65,0.7,0.75,0.8,0.85 };
// const float alpha[alphalen] = { 0.2, 0.3, 0.5, 0.7, 0.9, 1.0, 1.0 };
const float alpha[alphalen] = { 0.30, 0.30, 0.45, 0.50, 0.60, 0.65,0.7,0.75,0.8,0.85 }; // in use for scl*****
// const float alpha[alphalen] = { 0.30, 0.30, 0.45, 0.60, 0.70, 0.50,0.70,0.75,0.8,0.85 }; // in use for scl
// const float alpha[alphalen] = { 0.30, 0.30, 0.60, 0.50, 0.70, 0.60,0.8,0.75,0.8,0.85 }; // in use for scl
//const float alpha[alphalen] = { 0.20, 0.40, 0.45, 0.50, 0.60, 0.65,0.7,0.75,0.8,0.85 }; // in use for scl
// const float alpha[alphalen] = { 0,0.30,0.30,0.45,0.50,0.60,0.70,0.75,0.80,0.85 }; //in use for lsd
//const float alpha[alphalen] = { 0,0.30,0.30,0.45,0.50,0.60,0.70,0.75,0.80,0.85};
//const float alpha[alphalen] = { 0,0.30,0.30,0.45,0.50,0.60,0.65,0.70,0.75,0.80,0.85 };
const int usebeta = 0; // 1 for using beta; 0 for not using beta,
const int adaptive_beta = 0; // use adaptive betaf
const int betalen = 10;
//const float beta[betalen] = { 10,8,15,40,1,1,1,1 };
const float beta[betalen] = { 0.2,0.4,0.6,0.8,1,1,1,1,1,1 }; // 256QAM, 32 x 256 LSD T2
//const float beta[betalen] = { 1.0, 1.5, 2.0, 2.5, 3.0, 3, 3, 3, 3, 3 };
//(64,57)
const float a_beta[betalen] = { 0.60, 0.65, 0.70, 0.80, 0.90, 1.00, 1.00, 1.00 };
const float b_beta[betalen] = { 0.10, 0.15, 0.20, 0.35, 0.45, 0.55, 0.55, 0.55 };
const float k_beta[betalen] = { 1.00, 0.80, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10 };

//others
const int blk_cnt = 2;
const int CRCL = 0; // prev : 11
const int CRC_CHECK_ENABLE = 0;
const int PARITY_CHECK_ENABLE = 1;
const int CRC_VERTICAL = 0;
const int CRC_HORIZONTAL = 0;
const int L = 8;
const int L_WHOLE_CODEWORD = 1024;
const int L_FINAL_SEARCH = 200;
const int L_SUB_CODEWORD = 4; // 4 already seems good for SCL
const int L_MIMO_BLK = 64; // prev : 64
const string polarde = "SLD"; // DO NOT MODIFY! NO MATTER YOU ARE SIMULATING SLD OR SCL
const string polarde_array[6] = { "SLD","SCL","SCL","SCL","SCL","SCL" };
const int half_iter = 1; // half iteration at 2nd Pyndiah
const int bit_flip_pyndiah = 0;
const int double_stage = 0; // whether activate double-stage temporal decoding-
const int L_LSD = 32;
const int p = 10;
const int horizontal_ml = 1;
const int vertical_ml = 1;
const float max_llr_ratio = 3; // prev 2.0, 10.0
const int L_max = 8;
const int disable_max_ratio = false;
const int test_sc = (L==1);
const int test_cadsl = true;


//SISO or MIMO
const int atten = 2; // 1 for SISO; 2 for MIMO.  now only SISO
// const int atten_new = 2;

//SDD or JDD
const string system_archi = "SDD"; // JDD or SDD
const float beta_JDD = 0.4;
// Simulationsconst float Min_SNR = 23;
const float Min_SNR = 9; // last test (MIMO) : 10;
const float Max_SNR = 12;
const float SNR_step = 0.5;
const int errframe = 100;
const int record_frames = 100;
const int weirditer = 0;

const int restrict_llr = true; // prev: true
const int showllr = 0;
const int showllr_irr = 0;

const int irregular_coding = false;
const int num_time_constructions = 2;
const int K_time_array[num_time_constructions] = { 18, 18 }; // prev: 42,57 // prev: 16,26
const int border = 8; // 25 // prev: 21// prev: 15
const float amplify_ratio = 1;
const float amplify_ratio_array[5] = { 1.0,1.0,1.0,1.0,1.0 }; // 15
// const float amplify_ratio_array[5] = { 1.0,1,1.0,1.0,1.0 }; // 16
// const float amplify_ratio_array[5] = { 1.0,0.9,0.5,1.0,1.0 };// Tx16 Rx32
// const float amplify_ratio_array[5] = { 1.0,1.0,0.4,1.0,1.0 };// Tx16 Rx32
const bool irregular_iteration = true;
const float CLOG_RESTRICTION = 10000;
const int flag_information_mapping = false;
const int flag_separate_interleaving = false;


const int non_uniform_channel = 0;
const int USE_PLATFORM_CHANNEL = 0;
const int USE_ELAA_CHANNEL = 0;
const int N_FREQUENCY_SAMPLES = 1000; // platform : 624
const float diff_gain = 0;
const std::string MIMO_SVD_POWER_ALLOCATION_MODE = "SEPARATE";
const int mimo_svd_average_power_alloc = false;
const int mimo_waterfilling_power_alloc = false;
const int mimo_rate_matching_power_alloc = false;
const int MIMO_SVD_SEPARATE_AVERAGE_POWER_ALLOC = true;
const float MIMO_SVD_SEPARATE_POWER_RATIO = 0.4; // prev: 0.4 prev: 0.5(1103) // AWGN: 0.6
const float mimo_rate_matching_enlarge_ratio = 1.6; // prev : 1.6 // new prev 2.0
const int predefined_non_uniform = 1;
const int ENUMERATE_CODEWORD = 0;
// const int channel_index_predefined[Nv] = { 0,0,0,1,0,1,1,1,0,1,1,1,1,1,1,1,0,1,1,1,1,1,1,0,1,1,1,0,1,0,0,0};

//const int channel_index_predefined[Nv] = { 0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1 };
//
//const int channel_index_predefined[Nv] = { 0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1 };
// const int channel_index_predefined[Nv] = { 0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 };  // 25
// 
// 
//const int channel_index_predefined[Nv] = { 0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 };  // 24
//const int EXCLUSION_NUM = 4;
//const int exclusion_set[EXCLUSION_NUM] = { 28, 26, 22, 14 };
//const int PROTECTED_PARITY_NUM = 2;
//const int protected_parity_set[PROTECTED_PARITY_NUM] = { 0,1 };
//const int SPECIFIC_EXCLUSION = 1;

// const int channel_index_predefined[Nv] = { 0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 };  // 24
const int EXCLUSION_NUM = 4;
const int exclusion_set[EXCLUSION_NUM] = { 11,19,25,26 };
const int PROTECTED_PARITY_NUM = 2;
const int protected_parity_set[PROTECTED_PARITY_NUM] = { 0,4 };
const int SPECIFIC_EXCLUSION = 0;
const int OUTPUT_ERROR_PATTERN = 1;

// const int channel_index_predefined[Nv] = { 0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 };  // 23
// const int channel_index_predefined[Nv] = { 0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 };  // 21
// const int channel_index_predefined[Nv] = { 0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 };  // 20
// const int channel_index_predefined[Nv] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 };  // 18
// const int channel_index_predefined[Nv] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 };  // 16
// const int channel_index_predefined[Nv] = { 0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,1,1,1,1,1,1 };  // 15
// const int channel_index_predefined[Nv] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 };  // 15
// const int channel_index_predefined[Nv] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1 };  // 14
// const int channel_index_predefined[Nv] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1 };  // 13
// const int channel_index_predefined[32] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1 };  // 12
//const int channel_index_predefined[Nv] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1 };  // 11
// const int channel_index_predefined[Nv] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1 };  // 10
// const int channel_index_predefined[Nv] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1 };  // 9
// const int channel_index_predefined[Nv] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1 };  // 8
// const int channel_index_predefined[Nv] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1 };  // 7
// const int channel_index_predefined[Nv] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1 };  // 6
// const int channel_index_predefined[Nv] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1 };
const int channel_index_predefined[Nv] = { 0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1 }; // 6
// const int channel_index_predefined[Nv] = { 0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1 }; // 8
// const int channel_index_predefined[Nv] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1 }; // 2
// const int channel_index_predefined[Nv] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1 };  // 2
// const int channel_index_predefined[Nv] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1 };  // 1
// const int channel_index_predefined[Nv] = { 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
// const int channel_index_predefined[Nv] = { 0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1 };


// const int channel_index_predefined[Nv] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1 };
// const int channel_index_predefined[Nv] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 };
// const int channel_index_predefined[Nv] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 };
