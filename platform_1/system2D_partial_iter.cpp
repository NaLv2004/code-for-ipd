#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <fstream>
#include <vector>
#include <bitset>
#include <chrono> 
#include <complex>
#include <Eigen\Dense>
#include <Eigen\Core>
#include <unsupported\Eigen\KroneckerProduct>
#include "crc.h"
#include "printsetting.h"
#include "Polar_Encoder.h"
#include "Polar_Decoder.h"
#include "listSphereDecoding.h"
#include "Polar_Function.h"
#include "Polar_Constuction.h"
#include "modulation.h"
#include "setting.h"
#include "MIMO_Function.h"
#include "system_partial_iter.h"
#include <iomanip>
#include "mat_io.h"
using namespace Eigen;
using namespace std;


MatrixXf Fn(int N)
{
	int n = int(log(N) / log(2));
	MatrixXf m(2, 2);
	m(0, 0) = 1;
	m(1, 0) = 1;
	m(0, 1) = 0;
	m(1, 1) = 1;
	MatrixXf f = m;
	for (int i = 1; i < n; i++) {
		m = kroneckerProduct(m, f).eval();
	}
	return m;
}

void system2D_LSD_partial_iter(double* BER_array, double* FER_array)
{
	//freopen("result.txt", "w", stdout);
	printsetting();
	cout << "\nSNR	 " << "FER      BER" << endl;
	int N_Info = Kv * Kh - CRCL, Nmax;
	float CodeRate = (float)N_Info / (Nh * Nv);
	Nmax = max(Nh, Nv);
	// num of symbols per vertical codeword
	int Nsym = Nv / ModType;
	int symbol_length = ModType;
	// Initialization
	unsigned char* a_data = new unsigned char[Kv * Kh];// Information CodeWord
	int** u = new int* [Nv];// Original CodeWord u
	for (int i = 0; i < Nv; i++) { u[i] = new int[Nh]; }
	int** umid = new int* [Nv];// output CodeWord u
	for (int i = 0; i < Nv; i++) { umid[i] = new int[Nh]; }
	int** uout = new int* [Nv];// output CodeWord u
	for (int i = 0; i < Nv; i++) { uout[i] = new int[Nh]; }
	int** midx = new int* [Nv];//  Polar Encoded CodeWord x
	for (int i = 0; i < Nv; i++) { midx[i] = new int[Nh]; }
	int** x = new int* [Nv];//  Polar Encoded CodeWord x
	for (int i = 0; i < Nv; i++) { x[i] = new int[Nh]; }
	int** x_mmse_out = new int* [Nv]; // delete!
	for (int i = 0; i < Nv; i++) { x_mmse_out[i] = new int[Nh]; }
	std::complex<float>** x_mod = new std::complex<float>*[Nsym]; // modulated symbols x_mod // delete!
	for (int i = 0; i < Nsym; i++) { x_mod[i] = new std::complex<float>[Nh]; }
	std::complex<float>** y_mod = new std::complex<float>*[Nr]; // Output symbols after AWGN channels // delete!
	for (int i = 0; i < Nr; i++) { y_mod[i] = new std::complex<float>[Nh]; }
	float** y_mod_real = new float* [2 * Nr];
	for (int i = 0; i < 2 * Nr; i++) y_mod_real[i] = new float[Nh];
	float** y = new float* [Nv];// Output Signal After AWGN Channel
	for (int i = 0; i < Nv; i++) { y[i] = new float[Nh]; }
	float** y_mmse = new float* [2 * Nsym];
	for (int i = 0; i < 2 * Nsym; i++) { y_mmse[i] = new float[Nh]; }  // MMSE detector output // delete!
	float* y_mmse_tmp = new float[2 * Nsym];
	std::complex<float>* x_mod_tmp = new std::complex<float>[Nsym]; // modulated symbols of a single column // delete!
	std::complex<float>* y_mod_tmp = new std::complex<float>[Nr]; // received symbols of a single column // delete!
	float* y_mod_real_tmp = new float[2 * Nr];
	// float* y_mod_tmp = new float[Nr];  // delete!

	int* x_tmp = new int[Nv]; // delete!
	unsigned char* u_crc = new unsigned char[Kv * Kh];

	float** upLLR = new float* [Nv];
	for (int i = 0; i < Nv; i++) { upLLR[i] = new _MM_ALIGN16 float[Nh]; }
	//float** oriLLRv = new float* [Nh];
	//for (int i = 0; i < Nh; i++) { oriLLRv[i] = new _MM_ALIGN16 float[2 * Nv + 2]; }
	float** oriLLRh = new float* [Nv];
	for (int i = 0; i < Nv; i++) { oriLLRh[i] = new _MM_ALIGN16 float[Nh]; }
	float* oriLLR_tmp = new _MM_ALIGN16 float[Nv];  // delete!


	//*****************************************************************

	int** GGh = new int* [Nh];
	int** GGv = new int* [Nv];


	MatrixXf Gh = Fn(Nh);  // Generate Matrix of horizontal coding
	MatrixXf Gv = Fn(Nv);  // Generate Matrix of vertical coding
	MatrixXcf H(Nr, Nt); // Complex channel Matrix H: Nr x Nt

	// //MatrixX Definitions
	MatrixXf H_real(2 * Nr, 2 * Nt); // Real domain channel matrix
	MatrixXf received_signal(2 * Nr, 1);
	MatrixXf HTH(2 * Nt, 2 * Nt);
	MatrixXf Identity = MatrixXf::Identity(2 * Nt, 2 * Nt);
	MatrixXf output_signal(2 * Nt, 1);
	MatrixXf mmse_matrix(2 * Nt, 2 * Nr);

	// Matrix Definitions
	//Eigen::Matrix<float, 2 * Nr, 2 * Nt> H_real; // 实域信道矩阵
	//Eigen::Matrix<float, 2 * Nr, 1> received_signal; // 接收信号矩阵
	//Eigen::Matrix<float, 2 * Nt, 2 * Nt> HTH; // HTH矩阵
	//Eigen::Matrix<float, 2 * Nt, 2 * Nt> Identity = Eigen::Matrix<float, 2 * Nt, 2 * Nt>::Identity(); // 单位矩阵
	//Eigen::Matrix<float, 2 * Nt, 1> output_signal; // 输出信号矩阵
	//Eigen::Matrix<float, 2 * Nt, 2 * Nr> mmse_matrix; // MMSE矩阵


	for (int i = 0; i < Nh; i++) GGh[i] = new int[Nh];
	for (int i = 0; i < Nv; i++) GGv[i] = new int[Nv];
	for (int i = 0; i < Nh; i++)
		for (int j = 0; j < Nh; j++)
			GGh[i][j] = Gh(i, j);
	for (int i = 0; i < Nv; i++)
		for (int j = 0; j < Nv; j++)
			GGv[i][j] = Gv(i, j);
	//*****************************************************************
	vector<int> temph(Nh);
	vector<int> tempv(Nv);
	int* Ah = new int[Kh]; // Horizontal Info Bits
	int* Ach = new int[Nh - Kh]; // Horizontal Frozen Bits
	int* A_Ach = new int[Nh];
	int* Av = new int[Kv];
	int* Acv = new int[Nv - Kv];
	int* A_Acv = new int[Nv];

	//*****************************************************************
		// Main Function
	srand((unsigned int)time(0));
	int Num_Error, Num_Frame_Error;
	float Decoding_Duration;

	int Num_Error_new, Num_Frame_Error_new;
	float Decoding_Duration_new;

	clock_t T_start, T_finish;
	int SNR_count = 0;
	int anssc, anssd;
	for (double nSNR = Min_SNR; nSNR <= Max_SNR; nSNR += SNR_step)
	{
		// TEST: det err
		int n_det_err = 0; // n error bits at detector output
		int total_bit_cnt = 0;
		// TEST: det err
		int cntcnt = 0;
		//Initialization
		Num_Error = 0;
		Num_Frame_Error = 0;
		Decoding_Duration = 0;

		// double tmp = pow(10, nSNR / 10);
		float sigma_sq = 0;
		double tmp = 0;
		if (atten == 1)
		{
			tmp = pow(10, nSNR / 10);
			sigma_sq = (float)1 / (float)tmp / 2 / CodeRate;
		}
		if (atten == 2)
		{
			tmp = pow(10, nSNR / 10) * symbol_length * CodeRate * Nt;
			// tmp = pow(10, nSNR / 10) * symbol_length  * Nt;
			sigma_sq = (Nt * Nr) / tmp;
		}
		polar_codeconstruction(polarconh, Nh, Kh, GGh, (float)sqrt(sigma_sq), temph);
		polar_codeconstruction(polarconv, Nv, Kv, GGv, (float)sqrt(sigma_sq), tempv);

		// sort temph and tempv in ascending order
		for (int i = 0; i < Kh - 1; i++)
			for (int j = i + 1; j < Kh; j++)
				if (temph[i] > temph[j])
				{
					int ttt = temph[i];
					temph[i] = temph[j];
					temph[j] = ttt;
				}
		for (int i = 0; i < Kv - 1; i++)
			for (int j = i + 1; j < Kv; j++)
				if (tempv[i] > tempv[j])
				{
					int ttt = tempv[i];
					tempv[i] = tempv[j];
					tempv[j] = ttt;
				}

		// Decide the position of information bits and frozen bits
		for (int i = 0; i < Kh; i++) // Information Set
		{
			Ah[i] = temph[i];
			A_Ach[Ah[i]] = 1;
		}
		//	getchar();
		for (int i = 0; i < Nh - Kh; i++) // Frozen Bit
		{
			Ach[i] = temph[Kh + i];
			A_Ach[Ach[i]] = 0;
		}
		for (int i = 0; i < Kv; i++) // Information Set
		{
			Av[i] = tempv[i];
			A_Acv[Av[i]] = 1;
		}
		//	getchar();
		for (int i = 0; i < Nv - Kv; i++) // Frozen Bit
		{
			Acv[i] = tempv[Kv + i];
			A_Acv[Acv[i]] = 0;
		}
		// Frame Loop
		int passnum;
		int frame = 0;
		while (Num_Frame_Error < errframe)
		{
			frame++;
			//cout << frame << endl;
			int u_input[16] = { 0,     1     ,0,     1,
									0,     1,     0,     0,
									0 ,    1,     0,     0,
									1 ,    0,     0,     0 };
			for (int i = 0; i < N_Info; i++)
			{
				a_data[i] = (unsigned char)(rand() % 2);
				//		a_data[i] = u_input[i];
			}

			//cout << endl;
			if (CRCL) tx_append_crc(a_data, N_Info, CRCL);
			for (int i = 0; i < Nv; i++)
				for (int j = 0; j < Nh; j++)
					u[i][j] = 0;
#pragma omp parallel for
			for (int i = 0; i < Kv; i++)
			{
				for (int j = 0; j < Kh; j++)
				{
					u[Av[i]][Ah[j]] = (int)a_data[i * Kh + j];
					//cout << (int)u[Av[i]][Ah[j]] << " \n";
					//cout << Av[i] << "\n";
				}
				//cout << endl; 
			}
			//Polar Encoding
#pragma omp parallel for
			for (int i = 0; i < Nv; i++)
				PolarEncode_xor(midx[i], u[i], Nh);
#pragma omp parallel for
			for (int i = 0; i < Nh; i++)
			{
				int* tmpxv = new int[Nv]; int* tmpuv = new int[Nv];
				for (int j = 0; j < Nv; j++)
					tmpuv[j] = midx[j][i];
				PolarEncode_xor(tmpxv, tmpuv, Nv);
				for (int j = 0; j < Nv; j++)
					x[j][i] = tmpxv[j];
				delete[] tmpxv; delete[] tmpuv;
			}
			if (atten == 2)
			{
				// modulation
				// #pragma omp parallel for
				for (int i = 0; i < Nh; i++)
				{
					// Extracting a single column
					for (int j = 0; j < Nv; j++) { x_tmp[j] = x[j][i]; }
					modulation(x_tmp, x_mod_tmp, ModType, Nt);
					for (int j = 0; j < Nsym; j++) { x_mod[j][i] = x_mod_tmp[j]; }
					// test of transmitted bits
					/*if (i==0)
						for (int j = 0; j < Nv; j++)
						{
							if (!(j % ModType)) cout << endl;
							cout << x[j][i] << "  ";
						}*/
						// test of transmitted bits
				}
				// AWGN Channel (received y_mod )
				H = generateChannelMatrix(Nr, Nt);
				for (int i = 0; i < Nh; i++)
				{
					// Extracting a single column
					for (int j = 0; j < Nsym; j++) { x_mod_tmp[j] = x_mod[j][i]; }
					AWGN(Nr, Nt, x_mod_tmp, y_mod_tmp, H, sigma_sq);
					for (int j = 0; j < Nr; j++) { y_mod[j][i] = y_mod_tmp[j]; }
				}
				// Convert H and y to real domain
				convertHToReal(H, H_real);
				// TEST:H_real
				/*cout << "H_real";
				for (int i = 0; i < 2 * Nr; i++)
				{
					cout << endl;
					for (int j = 0; j < 2 * Nt; j++)
						printf("%.10f  ", H_real(i,j));
				}*/
				// TEST:H_real
				for (int i = 0; i < Nh; i++)
				{
					for (int j = 0; j < Nr; j++) { y_mod_tmp[j] = y_mod[j][i]; }
					convertYToReal(Nr, y_mod_tmp, y_mod_real_tmp);
					for (int j = 0; j < 2 * Nr; j++) { y_mod_real[j][i] = y_mod_real_tmp[j]; /*cout << y_mod_real_tmp[j] << endl;*/ }
				}
				// MMSE detection soft output
				for (int i = 0; i < Nh; i++)
				{
					for (int j = 0; j < 2 * Nr; j++) { y_mod_real_tmp[j] = y_mod_real[j][i]; }
					MMSEdetection(H_real, Identity, received_signal, mmse_matrix, output_signal, HTH, y_mod_real_tmp, y_mmse_tmp, sigma_sq);
					// TEST Detection
					/*if (i == 0)
					{
						save_mat_txt(H_real, "H_real.txt");
						save_mat_txt(received_signal, "rx.txt");
						save_mat_txt(mmse_matrix, "mat_mmse.txt");
						save_mat_txt(output_signal, "o_sym.txt");
						save_mat_txt(HTH, "HTH.txt");
					}*/
					//cout << "sigma2" << sigma_sq << endl;
					// TEST Detection
					for (int j = 0; j < 2 * Nt; j++) { y_mmse[j][i] = y_mmse_tmp[j]; /*cout << y_mmse_tmp[j] << endl;*/ }
				}

				// sym llr to bitwise llr
				// Detection is done for every column, that is, for every Nv bits
				// Interleaving is not currently considered
				for (int i = 0; i < Nh; i++)
				{
					// if (i == 0) { cout << "MMSE output LLRs"; }
					/*if (i == 0)
					{
						cout << "Detection Result" << endl;
					}
					if (i == 1)
					{
						cout << n_det_err << endl;
					}*/
					for (int j = 0; j < 2 * Nt; j++) { y_mmse_tmp[j] = y_mmse[j][i]; }
					llr_mmse_bit(Nt, y_mmse_tmp, oriLLR_tmp, sigma_sq, ModType);
					for (int j = 0; j < Nv; j++)
					{
						oriLLRh[j][i] = oriLLR_tmp[j];
						upLLR[j][i] = oriLLR_tmp[j];
						// TEST： det error
						x_mmse_out[j][i] = 0;
						if (oriLLR_tmp[j] < 0) { x_mmse_out[j][i] = 1; }
						if (x_mmse_out[j][i] != x[j][i]) { n_det_err++; }
						total_bit_cnt++;
						/*if (i == 0)
						{
							if (!(j % ModType)) cout << endl;
							cout << x_mmse_out[j][i] << " ";
						}
						cout << endl;*/
						// TEST: det error
					}
				}
			}
			if (atten == 1)
			{
				for (int i = 0; i < Nv; i++)
					for (int j = 0; j < Nh; j++)
					{
						// 0-1; 1-(-1)
						y[i][j] = (1 - 2 * x[i][j]) + (float)sqrt(sigma_sq) * (float)sampleNormal();
						oriLLRh[i][j] = 2 * y[i][j] / sigma_sq;
						//oriLLRv[j][i] = oriLLRh[i][j];
						upLLR[i][j] = oriLLRh[i][j];
					}
			}
			//Polar Decoding
			int blk_len = Nh / blk_cnt;
			for (int i_decod_blk = 0; i_decod_blk < blk_cnt; i_decod_blk++)
			{
				// column (spatial) decoding 
				int decodingi = 0;
				for (int iter = 1; iter <= softiter; iter++)
				{
					for (int i_decod_col = i_decod_blk * blk_len; i_decod_col < (i_decod_blk + 1) * blk_len; i_decod_col++)
					{
						// int iter = 1;
						decodingi = Nh - 1 - i_decod_col;
						//for SCL Decoding
						int Counth = 0, Countinfoh = 0;
						float Qsumunh = 0;
						int** sum = new int* [L];
						for (int i = 0; i < L; i++)  sum[i] = new _MM_ALIGN16 int[Nmax];
						float* LLR_in = new float[L];
						float* W = new float[2 * L];
						int* SCL_index = new int[L];
						int* better = new int[L];
						int* worse = new int[L];
						int* indexpm = new int[L];
						float* PM = new float[L];
						float** LLR = new float* [L];
						for (int i = 0; i < L; i++)  LLR[i] = new _MM_ALIGN16 float[2 * Nmax + 2];
						////////////////////////////////
						for (int j = 0; j < Nv; j++)
							for (int k = 0; k < L; k++)
								LLR[k][j] = upLLR[j][decodingi];
						//cout << LLR[0][j] << " ";

					//cout << endl;
						int index = 0;
						for (int k = 0; k < L; k++)  PM[k] = 0, indexpm[k] = k;
						//T_start = clock();
						decode_list(LLR, sum, (int)(log(Nv) / log(2)), A_Acv, 0, Nv, Kv, PM, L, (int)(log(L) / log(2)), LLR_in, W, SCL_index, better, worse, &Counth, &Countinfoh, &Qsumunh);
						//T_finish = clock();
						//Decoding_Duration += (T_finish - T_start);
						for (int i = 0; i < L - 1; i++)
							for (int j = i + 1; j < L; j++)
								if (PM[indexpm[i]] < PM[indexpm[j]])
								{
									int ttt = indexpm[i];
									indexpm[i] = indexpm[j];
									indexpm[j] = ttt;
								}
						//PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), N);
						index = indexpm[0];
						//}
						//calc distance s()

						if (softmethod == "Pyndiah") Pyndiahv(LLR, sum, oriLLRh, indexpm, iter * 2 - 2, decodingi, upLLR);
						if (softmethod == "MITSO") MITSOv(LLR, PM, sum, oriLLRh, indexpm, iter * 2 - 2, decodingi, upLLR, Qsumunh);

						////////delete SCL variables
						for (int i = 0; i < L; i++)
							delete[]LLR[i];
						delete[] LLR;
						for (int i = 0; i < L; i++)
							delete[]sum[i];
						delete[] sum;
						delete[] PM;
						delete[] LLR_in;
						delete[] W;
						delete[] SCL_index;
						delete[] better;
						delete[] worse;
						delete[] indexpm;

					}

					//row (temporal) decoding
					for (int decodingi = 0; decodingi < Nv; decodingi++)
					{
						//int iter = 0;
						//for SCL Decoding
						int Countv = 0, Countinfov = 0;
						float Qsumunv = 0;
						int** u1 = new int* [L];
						for (int i = 0; i < L; i++)  u1[i] = new _MM_ALIGN16 int[Nmax];
						int** u1_LSD = new int* [2 * L_LSD];
						for (int i = 0; i < 2 * L_LSD; i++) u1_LSD[i] = new int[Nmax];
						int** sum = new int* [L];
						for (int i = 0; i < L; i++)  sum[i] = new _MM_ALIGN16 int[Nmax];
						int** sum_LSD = new int* [2 * L_LSD];
						for (int i = 0; i < 2 * L_LSD; i++)  sum_LSD[i] = new int[Nmax];
						float* LLR_in = new float[L];
						int* u_decod = new int[Nh];
						int* n_paths = new int[Nh];
						float* W = new float[2 * L];
						int* SCL_index = new int[L];
						int* better = new int[L];
						int* worse = new int[L];
						int* indexpm = new int[L_LSD];
						float* PM = new float[L];
						float** LLR = new float* [L];
						for (int i = 0; i < L; i++)  LLR[i] = new _MM_ALIGN16 float[2 * Nmax + 2];

						for (int i = 0; i < L_LSD; i++) { indexpm[i] = i; }
						////////////////////////////////
						for (int j = 0; j < Nh; j++)
							for (int k = 0; k < L; k++)
								LLR[k][j] = upLLR[decodingi][j];
						//	cout << LLR[0][j] << " ";
						calc_npaths(A_Ach, Nh, n_paths, L_LSD);
						//	cout << endl;
						int index = 0;
						LSD_decode_outall(u_decod, u1_LSD, LLR[0], GGh, A_Ach, Nh, L_LSD, n_paths);
						// if (iter == softiter) PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), Nh);
						// for soft information exchange, inverse coding:
						for (int i = 0; i < L_LSD; i++) { PolarEncode_xor(sum_LSD[i], u1_LSD[0], Nh); }
						// if (iter == softiter) PolarEncode_xor(u1[0], sum_LSD[0], Nh);
						index = indexpm[0];
						// TEST
						/*cout << "SCL decoded bits" << endl;
						for (int i = 0; i < Nh; i++) { cout << sum[indexpm[0]][i] << "  "; }
						cout << endl;
						cout << "LSD decoded bits" << endl;
						for (int i = 0; i < Nh; i++) { cout << u1[0][i] << "  "; }
						cout << endl;*/
						// TEST
						//}
						//calc distance s()
						if (iter < softiter)
						{
							if (softmethod == "Pyndiah") Pyndiahh_LSD(LLR, sum_LSD, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR);
							if (softmethod == "MITSO") MITSOh(LLR, PM, sum_LSD, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR, Qsumunv);
							//	getchar();
						}
						else
							for (int j = 0; j < Kh; j++)
								umid[decodingi][Ah[j]] = u1_LSD[index][Ah[j]];
						////////delete SCL variables
						for (int i = 0; i < L; i++)
							delete[]u1[i];
						delete[] u1;
						for (int i = 0; i < 2 * L_LSD; i++)
							delete[]u1_LSD[i];
						delete[] u1_LSD;
						for (int i = 0; i < L; i++)
							delete[]LLR[i];
						delete[] LLR;
						for (int i = 0; i < L; i++)
							delete[]sum[i];
						delete[] sum;
						for (int i = 0; i < 2 * L_LSD; i++)
							delete[]sum_LSD[i];
						delete[] sum_LSD;
						delete[] PM;
						delete[] LLR_in;
						delete[] W;
						delete[] SCL_index;
						delete[] better;
						delete[] worse;
						delete[] indexpm;
						delete[] u_decod;
						delete[] n_paths;
					}
				}

			}
			// for (int decodingi = 0)
			int calcsum[Nv], calcu1[Nv];
			for (int i = 0; i < Kh; i++)
			{
				for (int j = 0; j < Nv; j++)
					calcsum[j] = umid[j][Ah[i]];
				PolarEncode_xor(calcu1, calcsum, Nv);
				for (int j = 0; j < Kv; j++)
					uout[Av[j]][Ah[i]] = calcu1[Av[j]];
			}
			//now without CRC
			int Temp_Error = Num_Error;
			for (int i = 0; i < Kv; i++)
			{
				for (int j = 0; j < Kh; j++)
				{
					//printf("%d ", (int)uout[Av[i]][Ah[j]]);
					Num_Error += abs(u[Av[i]][Ah[j]] - uout[Av[i]][Ah[j]]);
				}
				//cout << endl;
			}
			//getchar();
			if (Temp_Error < Num_Error) Num_Frame_Error++;
			if (frame % 10 == 0)
			{
				printf("\r<%d> SNR= %.2f FER= %d / %d = %f BER= %d / %d = %lf  Hard BER = %d / %d = %f",
					frame, nSNR, Num_Frame_Error, frame, (double)Num_Frame_Error / frame, Num_Error, (N_Info * frame), (double)Num_Error / N_Info / frame, n_det_err, total_bit_cnt, (double)n_det_err / total_bit_cnt);
				fflush(stdout);
			}
			//printf("%d %d", Num_Frame_Error, Num_Error);
			//getchar();
			//Decoding_Duration += (T_finish - T_start);
		}
		//Decoding_Duration = Decoding_Duration * 100 / CLOCKS_PER_SEC / frame;
		//Decoding_Duration = Decoding_Duration / CLOCKS_PER_SEC / Max_frame;

		//printf("%f  %f  %f  %f\n", nSNR, Decoding_Duration, (double)Num_Error / N_Info / Max_frame, (double)Num_Frame_Error / Max_frame);
		printf("%f %f %f\n", nSNR, (double)Num_Frame_Error / frame, (double)Num_Error / N_Info / frame);
		BER_array[SNR_count] = (double)Num_Error / N_Info / frame;
		FER_array[SNR_count] = (double)Num_Frame_Error / frame;
		SNR_count++;
	}
	//fclose(stdout);
	delete[] a_data;
	for (int i = 0; i < Nv; i++)
		delete[] u[i];
	delete[] u;
	for (int i = 0; i < Nv; i++)
		delete[] umid[i];
	delete[] umid;
	for (int i = 0; i < Nv; i++)
		delete[] uout[i];
	delete[] uout;
	for (int i = 0; i < Nv; i++)
		delete[] x[i];
	delete[] x;
	for (int i = 0; i < Nv; i++)
		delete[] x_mmse_out[i];
	delete[] x_mmse_out;
	for (int i = 0; i < Nv; i++)
		delete[] midx[i];
	delete[] midx;
	for (int i = 0; i < Nv; i++)
		delete[] y[i];
	delete[] y;
	delete[] u_crc;
	delete[] Ah;
	delete[] Ach;
	delete[] A_Ach;
	delete[] Av;
	delete[] Acv;
	delete[] A_Acv;
	for (int i = 0; i < Nv; i++)
		delete[] oriLLRh[i];
	delete[] oriLLRh;
	for (int i = 0; i < Nv; i++)
		delete[] upLLR[i];
	delete[] upLLR;
	for (int i = 0; i < Nh; i++)
		delete[] GGh[i];
	delete[] GGh;
	for (int i = 0; i < Nv; i++)
		delete[] GGv[i];
	delete[] GGv;
	for (int i = 0; i < Nsym; i++) {
		delete[] x_mod[i];
	}
	delete[] x_mod;
	for (int i = 0; i < Nr; i++) {
		delete[] y_mod[i];
	}
	delete[] y_mod;
	for (int i = 0; i < 2 * Nsym; i++) {
		delete[] y_mmse[i];
	}
	delete[] y_mmse;
	for (int i = 0; i < 2 * Nr; i++) {
		delete[] y_mod_real[i];
	}
	delete[] y_mod_real;
	delete[] x_mod_tmp;
	delete[] y_mod_tmp;
	delete[] y_mod_real_tmp;
	delete[] x_tmp;
	delete[] oriLLR_tmp;
}
