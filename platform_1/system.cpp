
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <random>
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
#include <iomanip>
#include <vector>
#include "mat_io.h"
#include"LSDecode.h"
#include <chrono>
#include "my_utils.h"
#include "IDD_Function.h"
#include "crc_list.h"
#include <mkl.h>
#include <random>
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
//ofstream y_out("y_out.txt");


void system1D(double* BER_array, double* FER_array)
{
	//freopen("result.txt", "w", stdout);
	printsetting();
	cout << "\nSNR	 " << "FER      BER" << endl;
	int N_Info = K - CRCL;
	float CodeRate = (float)N_Info / N;

	// Initialization
	unsigned char* a_data = new unsigned char[K];// Information CodeWord
	int* u = new int[N];// Original CodeWord u
	int* x = new int[N];// Polar Encoded CodeWord x
	float* y = new float[N];// Output Signal After AWGN Channel
	unsigned char* u_crc = new unsigned char[K];

	float** LLR = new float* [L];
	for (int i = 0; i < L; i++) { LLR[i] = new _MM_ALIGN16 float[2 * N + 2]; }
	int** sum = new int* [L];
	for (int i = 0; i < L; i++) { sum[i] = new _MM_ALIGN16 int[N]; }
	int** u1 = new int* [L];
	for (int i = 0; i < L; i++) { u1[i] = new _MM_ALIGN16 int[N]; }
	//*****************************************************************
	int** GG = new int* [N];

	MatrixXf G = Fn(N);
	for (int i = 0; i < N; i++)	GG[i] = new int[N];
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			GG[i][j] = G(i, j);

	//*****************************************************************
	vector<int> temp(N);

	int* A = new int[K];
	int* Ac = new int[N - K];
	int* A_Ac = new int[N];

	//for SCL Decoding
	float* LLR_in = new float[L];
	float* W = new float[2 * L];
	int* SCL_index = new int[L];
	int* better = new int[L];
	int* worse = new int[L];
	int* indexpm = new int[L];
	float* PM = new float[L];

	// for LSD Decoding
	int* n_paths = new int[N];
	//*****************************************************************
		// Main Function
		// srand((unsigned int)time(0));
	srand((unsigned int)0);
	int i, j;
	int Num_Error, Num_Frame_Error;
	float Decoding_Duration;

	int Num_Error_new, Num_Frame_Error_new;
	float Decoding_Duration_new;

	clock_t T_start, T_finish;
	int SNR_count = 0;
	int anssc, anssd;
	for (double nSNR = Min_SNR; nSNR <= Max_SNR; nSNR += SNR_step)
	{
		int cntcnt = 0;
		//Initialization
		Num_Error = 0;
		Num_Frame_Error = 0;
		Decoding_Duration = 0;

		double tmp = pow(10, nSNR / 10);
		float sigma_sq = (float)1 / (float)tmp / 2 / CodeRate;

		polar_codeconstruction(polarcon, N, K, GG, (float)sqrt(sigma_sq), temp);

		for (int i = 0; i < K - 1; i++)
			for (int j = i + 1; j < K; j++)
				if (temp[i] > temp[j])
				{
					int ttt = temp[i];
					temp[i] = temp[j];
					temp[j] = ttt;
				}


		for (int i = 0; i < K; i++) // Information Set
		{
			A[i] = temp[i];
			A_Ac[A[i]] = 1;
			//cout << "infoset[0][" << temp[i] << "]=1;\n";
			//cout << temp[i] << ",";
		}

		//	getchar();
		for (int i = 0; i < N - K; i++) // Frozen Bit
		{
			Ac[i] = temp[K + i];
			A_Ac[Ac[i]] = 0;
		}
		calc_npaths(A_Ac, N, n_paths, L_LSD);
		// TEST:
		/*for (int i = 0; i < N; i++)
		{
			cout << n_paths[i] << endl;
		}*/

		/*for (int i = 0; i < 127; i++) 	cout << A_Ac[i] << ",";*/
	//		cout << A_Ac[127] << endl;
	//	}
	//	getchar();

		// Frame Loop
		int passnum;
		int frame = 0;
		while (Num_Frame_Error < errframe)
		{
			frame++;

			for (i = 0; i < N_Info; i++)
			{
				// original info bits
				a_data[i] = (unsigned char)(rand() % 2);
				//	a_data[i] = u_input[i];
				//	cout << (int)a_data[i] << " ";
			}
			//cout << endl;
			if (CRCL) tx_append_crc(a_data, N_Info, CRCL);

			for (i = 0; i < K; i++)
			{
				u[A[i]] = (int)a_data[i];
				// cout << (int)A[i]<< " ";
			}
			//cout << endl;
			//getchar();
			for (i = 0; i < N - K; i++)//2048-1024
			{
				u[Ac[i]] = 0;
			}
			/*for (i = 0; i < K; i++)//2048-1024
			{
				cout << u[A[i]] << " ";
			}
			cout << "\n";
			getchar();*/
			//Polar Encoding
			PolarEncode_xor(x, u, N);

			// AWGN Channel
			for (i = 0; i < N; i++)
			{
				y[i] = (1 - 2 * x[i]) + (float)sqrt(sigma_sq) * (float)sampleNormal();
			}

			//Polar Decoding
			for (i = 0; i < N; i++)
			{
				for (j = 0; j < L; j++)
				{
					LLR[j][i] = 2 * y[i] / sigma_sq;
					//LLR[j][i] = -y[i];
				//	LLR[j][i] = LLR_input[i];
				}
				//u[i] = u_input[i];
			}
			//float* myLLR = new float[10];
			// LSD Decoding
			//TEST
			/*cout << endl;
			cout << "Codeword:" << endl;
			for (int i = 0; i < N; i++) { cout << u[i] << "  "; }
			cout << endl;*/
			//
			// listSphereDecode(u, LLR[0], GG, A_Ac, N, L_LSD, n_paths);
			int index = 0, C1D = 0, Count_info = 0;
			if (L == 1)
			{
				T_start = clock();
				decode(*LLR, *sum, (int)(log(N) / log(2)), A_Ac, 0, N, K, &C1D);
				PolarEncode_xor(*u1, *sum, N);
				T_finish = clock();
			}
			else
			{
				bool markcrc = false;
				for (int k = 0; k < L; k++)  PM[k] = 0;
				for (i = 0; i < L; i++)
					indexpm[i] = i;
				T_start = clock();
				float QQQ = 0;
				decode_list(LLR, sum, (int)(log(N) / log(2)), A_Ac, 0, N, K, PM, L, (int)(log(L) / log(2)), LLR_in, W, SCL_index, better, worse, &C1D, &Count_info, &QQQ);
				T_finish = clock();
				Decoding_Duration += (T_finish - T_start);
				for (i = 0; i < L - 1; i++)
					for (int j = i + 1; j < L; j++)
						if (PM[indexpm[i]] < PM[indexpm[j]])
						{
							int ttt = indexpm[i];
							indexpm[i] = indexpm[j];
							indexpm[j] = ttt;
						}
				PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), N);
				index = indexpm[0];
				if (CRCL > 0) {
					int A_cnt = 0;
					for (int i = 0; i < N; i++) {
						if (A_Ac[i] == 1) { u_crc[A_cnt] = (unsigned char)u1[index][i]; A_cnt++; }
					}
					if (rx_check_crc(u_crc, K, CRCL) == false) {
						int cnt = 1;
						while (cnt < L) {
							PolarEncode_xor(u1[indexpm[cnt]], sum[indexpm[cnt]], N);
							A_cnt = 0;
							for (int i = 0; i < N; i++) {
								if (A_Ac[i] == 1) { u_crc[A_cnt] = (unsigned char)u1[indexpm[cnt]][i]; A_cnt++; }
							}
							if (rx_check_crc(u_crc, K, CRCL) == true) {
								index = indexpm[cnt]; break;
							}
							cnt++;
						}
					}
				}
			}
			int Temp_Error = Num_Error;
			for (i = 0; i < K - CRCL; i++)
				Num_Error += abs(u1[index][A[i]] - u[A[i]]);
			if (Temp_Error < Num_Error) Num_Frame_Error++;
			//	cout << Num_Frame_Error;
			//	getchar();
			Decoding_Duration += (T_finish - T_start);
		}
		Decoding_Duration = Decoding_Duration * 100 / CLOCKS_PER_SEC / frame;
		//		Decoding_Duration = Decoding_Duration / CLOCKS_PER_SEC / Max_frame;

				//printf("%f  %f  %f  %f\n", nSNR, Decoding_Duration, (double)Num_Error / N_Info / Max_frame, (double)Num_Frame_Error / Max_frame);
		printf("%f %f %f\n", nSNR, (double)Num_Frame_Error / frame, (double)Num_Error / N_Info / frame);

		BER_array[SNR_count] = (double)Num_Error / N_Info / frame;
		FER_array[SNR_count] = (double)Num_Frame_Error / frame;
		SNR_count++;
	}
	fclose(stdout);
	delete[] a_data;
	delete[] u;
	delete[] x;
	delete[] y;
	delete[] u_crc;
	delete[] A;
	delete[] Ac;
	delete[] A_Ac;
	for (int i = 0; i < L; i++)
		delete[]LLR[i];
	delete[] LLR;
	for (int i = 0; i < L; i++)
		delete[]sum[i];
	delete[] sum;
	for (int i = 0; i < L; i++)
		delete[]u1[i];
	delete[] u1;
	delete[] PM;
	delete[] LLR_in;
	delete[] W;
	delete[] SCL_index;
	delete[] better;
	delete[] worse;
	for (int i = 0; i < N; i++)
		delete[]GG[i];
	delete[] GG;
	delete[] indexpm;
	delete[] n_paths;
}

void system1D_LSD(double* BER_array, double* FER_array)
{
	//freopen("result.txt", "w", stdout);
	printsetting();
	cout << "\nSNR	 " << "FER      BER" << endl;
	int N_Info = K - CRCL;
	float CodeRate = (float)N_Info / N;

	// Initialization
	unsigned char* a_data = new unsigned char[K];// Information CodeWord
	int* u = new int[N];// Original CodeWord u
	int* x = new int[N];// Polar Encoded CodeWord x
	float* y = new float[N];// Output Signal After AWGN Channel
	unsigned char* u_crc = new unsigned char[K];

	// error at each bit position
	int* nberr_everybit = new int[N];
	for (int i = 0; i < N; i++) { nberr_everybit[i] = 0; }
	float** LLR = new float* [L];
	for (int i = 0; i < L; i++) { LLR[i] = new _MM_ALIGN16 float[2 * N + 2]; }
	int** sum = new int* [L];
	for (int i = 0; i < L; i++) { sum[i] = new _MM_ALIGN16 int[N]; }
	int** u1 = new int* [L];
	for (int i = 0; i < L; i++) { u1[i] = new _MM_ALIGN16 int[N]; }
	//*****************************************************************
	int** GG = new int* [N];

	MatrixXf G = Fn(N);
	for (int i = 0; i < N; i++)	GG[i] = new int[N];
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			GG[i][j] = G(i, j);

	//*****************************************************************
	vector<int> temp(N);

	int* A = new int[K];
	int* Ac = new int[N - K];
	int* A_Ac = new int[N];

	//for SCL Decoding
	float* LLR_in = new float[L];
	float* W = new float[2 * L];
	int* SCL_index = new int[L];
	int* better = new int[L];
	int* worse = new int[L];
	int* indexpm = new int[L];
	float* PM = new float[L];

	// for LSD Decoding
	int* n_paths = new int[N];
	int* u_decod = new int[N];
	//*****************************************************************
		// Main Function
		// srand((unsigned int)time(0));
	srand((unsigned int)0);
	int i, j;
	int Num_Error, Num_Frame_Error;
	float Decoding_Duration;

	int Num_Error_new, Num_Frame_Error_new;
	float Decoding_Duration_new;

	clock_t T_start, T_finish;
	int SNR_count = 0;
	int anssc, anssd;
	for (double nSNR = Min_SNR; nSNR <= Max_SNR; nSNR += SNR_step)
	{
		for (int i = 0; i < N; i++) { nberr_everybit[i] = 0; }
		int cntcnt = 0;
		//Initialization
		Num_Error = 0;
		Num_Frame_Error = 0;
		Decoding_Duration = 0;

		double tmp = pow(10, nSNR / 10);
		float sigma_sq = (float)1 / (float)tmp / 2 / CodeRate;

		polar_codeconstruction(polarcon, N, K, GG, (float)sqrt(sigma_sq), temp);

		for (int i = 0; i < K - 1; i++)
			for (int j = i + 1; j < K; j++)
				if (temp[i] > temp[j])
				{
					int ttt = temp[i];
					temp[i] = temp[j];
					temp[j] = ttt;
				}


		for (int i = 0; i < K; i++) // Information Set
		{
			A[i] = temp[i];
			A_Ac[A[i]] = 1;
			//cout << "infoset[0][" << temp[i] << "]=1;\n";
			//cout << temp[i] << ",";
		}
		//	getchar();
		for (int i = 0; i < N - K; i++) // Frozen Bit
		{
			Ac[i] = temp[K + i];
			A_Ac[Ac[i]] = 0;
		}
		calc_npaths(A_Ac, N, n_paths, L_LSD);
		// TEST:
		/*for (int i = 0; i < N; i++)
		{
			cout << n_paths[i] << endl;
		}*/

		/*for (int i = 0; i < 127; i++) 	cout << A_Ac[i] << ",";*/
	//		cout << A_Ac[127] << endl;
	//	}
	//	getchar();

		// Frame Loop
		int passnum;
		int frame = 0;
		while (Num_Frame_Error < errframe)
		{
			frame++;
			bool correct = 1;
			int nberr_tmp = 0;
			for (i = 0; i < N_Info; i++)
			{
				// original info bits
				a_data[i] = (unsigned char)(rand() % 2);
				//	a_data[i] = u_input[i];
				//	cout << (int)a_data[i] << " ";
			}
			//cout << endl;
			if (CRCL) tx_append_crc(a_data, N_Info, CRCL);

			for (i = 0; i < K; i++)
			{
				u[A[i]] = (int)a_data[i];
				// cout << (int)A[i]<< " ";
			}
			//cout << endl;
			//getchar();
			for (i = 0; i < N - K; i++)//2048-1024
			{
				u[Ac[i]] = 0;
			}
			/*for (i = 0; i < K; i++)//2048-1024
			{
				cout << u[A[i]] << " ";
			}
			cout << "\n";
			getchar();*/
			//Polar Encoding
			PolarEncode_xor(x, u, N);

			// AWGN Channel
			for (i = 0; i < N; i++)
			{
				y[i] = (1 - 2 * x[i]) + (float)sqrt(sigma_sq) * (float)sampleNormal();
			}

			//Polar Decoding
			for (i = 0; i < N; i++)
			{
				for (j = 0; j < L; j++)
				{
					LLR[j][i] = 2 * y[i] / sigma_sq;
					//LLR[j][i] = -y[i];
				//	LLR[j][i] = LLR_input[i];
				}
				//u[i] = u_input[i];
			}
			//float* myLLR = new float[10];
			// LSD Decoding
			//TEST
			/*cout << endl;
			cout << "Codeword:" << endl;
			for (int i = 0; i < N; i++) { cout << u[i] << "  "; }
			cout << endl;*/
			//
			listSphereDecode_step(u_decod, LLR[0], GG, A_Ac, N, L_LSD, n_paths);
			//int index = 0, C1D = 0, Count_info = 0;
			for (int i = 0; i < N; i++)
			{
				if (u_decod[i] != u[i])
				{
					nberr_tmp++;
					correct = 0;
					nberr_everybit[i]++;
				}
			}
			Num_Error += nberr_tmp;
			if (!correct) Num_Frame_Error++;

			if (frame % 100 == 0 && frame > 0)
			{
				printf("\r<%d> SNR= %.2f FER= %d / %d = %f BER= %d / %d = %lf",
					frame, nSNR, Num_Frame_Error, frame, (double)Num_Frame_Error / frame, Num_Error, (N_Info * frame), (double)Num_Error / N_Info / frame);
				fflush(stdout);
			}

			// if(!(frame%10) && frame > 0){ printf("%f %f %f\n", nSNR, (double)Num_Frame_Error / frame, (double)Num_Error / N_Info / frame); }
		}
		// Decoding_Duration = Decoding_Duration * 100 / CLOCKS_PER_SEC / frame;
		//		Decoding_Duration = Decoding_Duration / CLOCKS_PER_SEC / Max_frame;

				//printf("%f  %f  %f  %f\n", nSNR, Decoding_Duration, (double)Num_Error / N_Info / Max_frame, (double)Num_Frame_Error / Max_frame);
		printf("\n");
		printf("%f %f %f\n", nSNR, (double)Num_Frame_Error / frame, (double)Num_Error / N_Info / frame);
		for (int i = 0; i < N; i++) { cout << "  " << nberr_everybit[i]; }
		cout << endl;
		BER_array[SNR_count] = (double)Num_Error / N_Info / frame;
		FER_array[SNR_count] = (double)Num_Frame_Error / frame;
		SNR_count++;
	}
	// fclose(stdout);
	delete[] a_data;
	delete[] u;
	delete[] x;
	delete[] y;
	delete[] u_crc;
	delete[] A;
	delete[] Ac;
	delete[] A_Ac;
	for (int i = 0; i < L; i++)
		delete[]LLR[i];
	delete[] LLR;
	for (int i = 0; i < L; i++)
		delete[]sum[i];
	delete[] sum;
	for (int i = 0; i < L; i++)
		delete[]u1[i];
	delete[] u1;
	delete[] PM;
	delete[] LLR_in;
	delete[] W;
	delete[] SCL_index;
	delete[] better;
	delete[] worse;
	for (int i = 0; i < N; i++)
		delete[]GG[i];
	delete[] GG;
	delete[] indexpm;
	delete[] n_paths;
	delete[] u_decod;
	delete[] nberr_everybit;
}


void system2D_LSD(double* BER_array, double* FER_array)
{
	std::ofstream file_total("total_timing.txt", std::ios::app);
	std::ofstream file_det("det_timing.txt", std::ios::app);
	std::ofstream file_dec("dec_timing.txt", std::ios::app);
	std::ofstream file_err_num("total_frames.txt", std::ios::app);
	auto det_end = std::chrono::system_clock::now();
	
	/*if (!file.is_open()) {
		std::cerr << "Error opening file: " << file_name << std::endl;
		return;
	}*/

	//file << setw(20) << SNR << setw(20) << BER << setw(20) << FER << "\n";

	//freopen("result.txt", "w", stdout);
	printsetting();
	cout << "\nSNR	 " << "FER      BER" << endl;
	int N_Info = Kv * Kh - CRCL, Nmax;
	//Interleaving pattern
	std::vector<int> intrlv_pattern(Nv);
	std::vector<int>* intrlv_pattern_2D = new vector<int>[Nh];
	for (int i = 0; i < Nh; i++) intrlv_pattern_2D[i].resize(Nv);
	float CodeRate = (float)N_Info / (Nh * Nv);
	Nmax = max(Nh, Nv);
	// num of symbols per vertical codeword
	int Nsym = 0;
	if (MOD_DIM == "Spatial") { Nsym = Nv / ModType; }
	else if (MOD_DIM == "Temporal") { Nsym = Nh / ModType; }
	int symbol_length = ModType;

	// Size of a mimo block is given by height(spatial) * length(temporal)
	int mimo_blk_height = 0;
	int mimo_blk_length = 0;
	if (MOD_DIM == "Spatial")
	{
		mimo_blk_height = Nsym;
		mimo_blk_length = Nh;
	}
	else if (MOD_DIM == "Temporal")
	{
		mimo_blk_height = Nv;
		mimo_blk_length = Nsym;
	}
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

	// Arrays for MIMO
	int** x_mmse_out = new int* [Nv]; // delete!
	for (int i = 0; i < Nv; i++) { x_mmse_out[i] = new int[Nh]; }

	//std::complex<float>** x_mod = new std::complex<float>* [Nsym]; // modulated symbols (spatial modulation) DELETE!
	//for (int i = 0; i < Nsym; i++) { x_mod[i] = new std::complex<float>[Nh]; }
	//std::complex<float>** x_mod_temporal = new std::complex<float>*[Nv]; // modulated symbols (temporal modulation) DELETE!
	//for (int i = 0; i < Nv; i++) x_mod_temporal[i] = new std::complex<float>[Nsym];
	//std::complex<float>** y_mod = new std::complex<float>* [Nr]; // Output symbols after AWGN channels // delete!
	//for (int i = 0; i < Nr; i++) { y_mod[i] = new std::complex<float>[Nh]; }
	//float** y_mod_real = new float* [2 * Nr];
	//for (int i = 0; i < 2*Nr; i++) y_mod_real[i] = new float[Nh];
	//float** y = new float* [Nv];// Output Signal After AWGN Channel
	//for (int i = 0; i < Nv; i++) { y[i] = new float[Nh]; }
	//float** y_mmse = new float* [2*Nsym];
	//for (int i = 0; i < 2*Nsym; i++) { y_mmse[i] = new float[Nh]; }  // MMSE detector output // delete!
	//float* y_mmse_tmp = new float[2*Nsym];
	//std::complex<float>* x_mod_tmp = new std::complex<float>[Nsym]; // modulated symbols of a single column // delete!
	//std::complex<float>* y_mod_tmp = new std::complex<float>[Nr]; // received symbols of a single column // delete!
	//float* y_mod_real_tmp = new float[2*Nr];
	//// float* y_mod_tmp = new float[Nr];  // delete!


	std::complex<float>** x_mod = new std::complex<float>*[mimo_blk_height]; // modulated symbols (spatial modulation) DELETE!
	for (int i = 0; i < mimo_blk_height; i++) { x_mod[i] = new std::complex<float>[mimo_blk_length]; }
	std::complex<float>** y_mod = new std::complex<float>*[Nr]; // Output symbols after AWGN channels // delete!
	for (int i = 0; i < Nr; i++) { y_mod[i] = new std::complex<float>[mimo_blk_length]; }
	float** y_mod_real = new float* [2 * Nr];
	for (int i = 0; i < 2 * Nr; i++) y_mod_real[i] = new float[mimo_blk_length];
	float** y = new float* [Nv];// Output Signal After AWGN Channel
	for (int i = 0; i < Nv; i++) { y[i] = new float[Nh]; }
	float** y_mmse = new float* [2 * mimo_blk_height];
	for (int i = 0; i < 2 * mimo_blk_height; i++) { y_mmse[i] = new float[mimo_blk_length]; }  // MMSE detector output // delete!
	float* y_mmse_tmp = new float[2 * mimo_blk_height];
	std::complex<float>* x_mod_tmp = new std::complex<float>[mimo_blk_height]; // modulated symbols of a single column // delete!
	std::complex<float>* y_mod_tmp = new std::complex<float>[Nr]; // received symbols of a single column // delete!
	float* y_mod_real_tmp = new float[2 * Nr];
	// float* y_mod_tmp = new float[Nr];  // delete!
	float** y_mod_real_temp_omp = new float* [mimo_blk_length]; // new delete!
	for (int i = 0; i < mimo_blk_length; i++) { y_mod_real_temp_omp[i] = new float[2 * Nr]; }
	float** y_mmse_temp_omp = new float* [mimo_blk_length]; // new delete!
	for (int i = 0; i < mimo_blk_length; i++) { y_mmse_temp_omp[i] = new float[2 * mimo_blk_height]; }


	int* x_tmp = new int[Nv]; // delete!

	// paths selected by k-best detection


	// CRC bits
	unsigned char* u_crc = new unsigned char[Kv * Kh];

	float** upLLR = new float* [Nv];
	for (int i = 0; i < Nv; i++) { upLLR[i] = new _MM_ALIGN16 float[Nh]; }
	//float** oriLLRv = new float* [Nh];
	//for (int i = 0; i < Nh; i++) { oriLLRv[i] = new _MM_ALIGN16 float[2 * Nv + 2]; }
	float** oriLLRh = new float* [Nv];
	for (int i = 0; i < Nv; i++) { oriLLRh[i] = new _MM_ALIGN16 float[Nh]; }
	int len_ori_llr_tmp = 0;
	if (MOD_DIM == "Spatial") { len_ori_llr_tmp = Nv; }
	else if (MOD_DIM == "Temporal") { len_ori_llr_tmp = ModType * Nt; }
	float* oriLLR_tmp = new _MM_ALIGN16 float[len_ori_llr_tmp];  // delete!
	//float** oriLLR_kbest = new float* [Nh];
	float** oriLLR_kbest = new float* [mimo_blk_length];
	for (int i = 0; i < Nh; i++) { oriLLR_kbest[i] = new float[len_ori_llr_tmp]; }
	double* sum_sdis_array = new double[Nv];
	for (int i = 0; i < Nv; i++) sum_sdis_array[i] = 0;

	//*****************************************************************

	int** GGh = new int* [Nh];
	int** GGv = new int* [Nv];


	MatrixXf Gh = Fn(Nh);  // Generate Matrix of horizontal coding
	MatrixXf Gv = Fn(Nv);  // Generate Matrix of vertical coding
	MatrixXcf H(Nr, Nt); // Complex channel Matrix H: Nr x Nt
	//Matrix<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic> H = MatrixXcf::Random(Nr, Nt); // Real domain channel matrix
	//// //MatrixX Definitions
	MatrixXf H_real(2 * Nr, 2 * Nt); // Real domain channel matrix
	MatrixXf received_signal(2 * Nr, 1);
	MatrixXf HTH(2 * Nt, 2 * Nt);
	MatrixXf Identity = MatrixXf::Identity(2 * Nt, 2 * Nt);
	MatrixXf output_signal(2 * Nt, 1);
	MatrixXf mmse_matrix(2 * Nt, 2 * Nr);
	MatrixXf HTHIinv(2 * Nt, 2 * Nt);

	// for KBest
	MatrixXf R(2 * Nt, 2 * Nt);
	MatrixXf Q(2 * Nr, 2 * Nt);
	MatrixXf z(2 * Nt, 1);



	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> H_real = MatrixXf::Random(2 * Nr, 2 * Nt); // Real domain channel matrix
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> received_signal = MatrixXf::Random(2 * Nr,1);
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> HTH = MatrixXf::Random(2 * Nt, 2 * Nt);
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> Identity = MatrixXf::Random(2 * Nt, 2 * Nt);
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> output_signal = MatrixXf::Random(2 * Nt, 1);
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> mmse_matrix = MatrixXf::Random(2 * Nt, 2 * Nr);


	// Matrix Definitions
	//Eigen::Matrix<float, 2 * Nr, 2 * Nt> H_real; // ʵ���ŵ�����
	//Eigen::Matrix<float, 2 * Nr, 1> received_signal; // �����źž���
	//Eigen::Matrix<float, 2 * Nt, 2 * Nt> HTH; // HTH����
	//Eigen::Matrix<float, 2 * Nt, 2 * Nt> Identity = Eigen::Matrix<float, 2 * Nt, 2 * Nt>::Identity(); // ��λ����
	//Eigen::Matrix<float, 2 * Nt, 1> output_signal; // ����źž���
	//Eigen::Matrix<float, 2 * Nt, 2 * Nr> mmse_matrix; // MMSE����


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
	// srand((unsigned int)time(0));
	std::random_device rd;              // ��ȡ���������
	std::mt19937 gen(rd());             // ʹ�����ӳ�ʼ��mt19937������
	std::uniform_int_distribution<> distrib(1, 1000000);
	int Num_Error, Num_Frame_Error;
	float Decoding_Duration;

	int Num_Error_new, Num_Frame_Error_new;
	float Decoding_Duration_new;

	clock_t T_start, T_finish;
	int SNR_count = 0;
	int anssc, anssd;
	for (double nSNR = Min_SNR; nSNR <= Max_SNR; nSNR += SNR_step)
	{
		for (int i = 0; i < Nv; i++) sum_sdis_array[i] = 0;
		double sdis_sum = 0;
		// TEST: det err
		int n_det_err = 0; // n error bits at detector output
		int n_det_err_mrb = 0; // error bits at detector output(MRBs)
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
		// test: Horizontal frozen bits

		/*cout << "Frozen Bits" <<endl;
		for (int i = 0; i < Nh - Kh; i++)
		{
			cout << setw(5) << Acv[i];
		}
		cout << endl;*/
		// test: Horizontal frozen bits
		// Frame Loop
		int passnum;
		int frame = 0;
		while (Num_Frame_Error < errframe)
		{
			/*test: processing time*/
			//auto start = std::chrono::system_clock::now();
			/*if (frame == 10) {
				cout << "AAAAAAAAAA" << endl;
			}*/
			frame++;
			// generate interleaving pattern
			if (flag_interleave)
			{
				generateIntrlvPattern(intrlv_pattern, Nv);
				for (int i = 0; i < Nh; i++) { generateIntrlvPattern(intrlv_pattern_2D[i], Nv); }
			}

			//cout << frame << endl;
			int u_input[16] = { 0,     1     ,0,     1,
									0,     1,     0,     0,
									0 ,    1,     0,     0,
									1 ,    0,     0,     0 };
			for (int i = 0; i < N_Info; i++)
			{
				a_data[i] = (unsigned char)(distrib(gen) % 2);
				//		a_data[i] = u_input[i];
			}

			//cout << endl;
			if (CRCL) tx_append_crc(a_data, N_Info, CRCL);
			for (int i = 0; i < Nv; i++)
				for (int j = 0; j < Nh; j++)
					u[i][j] = 0;
			//omp_set_num_threads(4);
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
			// test:original info
			//cout << "original info u:" << endl;
			////for (int j = 0; j < Nh; j++) { cout << u[Av[0]][j] << " "; }
			//for (int i = 0; i < Nv; i++)
			//{
			//	cout << endl;
			//	for (int j = 0; j < Nh; j++)
			//		cout << u[i][j];
			//}
			//cout << endl;
			// Polar Encoding
#pragma omp parallel for
			for (int i = 0; i < Nv; i++)
				PolarEncode_xor(midx[i], u[i], Nh);
			// test:horizontal encoding
			//cout << "u after horizontal encoding:" << endl;
			////for (int j = 0; j < Nh; j++) { cout << midx[Av[0]][j] << " "; }
			//for (int i = 0; i < Nv; i++)
			//{
			//	cout << endl;
			//	for (int j = 0; j < Nh; j++)
			//		cout << midx[i][j];
			//}
			//cout << endl;
			//cout << endl;
			// test:horizontal encoding
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
			/*for (int i = 0; i < 2 * Nr; i++)
				for (int j = 0; j < 2 * Nt; j++)
					H_real(i, j) = 0;*/
					/*MatrixXf H_real_test(2 * Nr, 2 * Nt);
					MatrixXf M = H_real_test.transpose();
					MatrixXf N = H_real_test;
					HTHIinv = M * N;*/
					//HTHIinv = M * H_real;
					// test:vertical encoding
					//cout << "u after vertical encoding:" << endl;
					////for (int j = 0; j < Nh; j++) { cout << x[Av[0]][j] << " "; }
					//for (int i = 0; i < Nv; i++)
					//{
					//	cout << endl;
					//	for (int j = 0; j < Nh; j++)
					//		cout << x[i][j];
					//}
					//cout << endl;
					//cout << endl;
					// test:vertical encoding
			if (atten == 2)
			{
				// modulation
				// #pragma omp parallel for
				// Interleaving
				if (flag_interleave)
				{
					for (int j = 0; j < Nh; j++)
					{
						for (int i = 0; i < Nv; i++) { x_tmp[i] = x[i][j]; }
						randIntrlv(intrlv_pattern_2D[j], x_tmp, Nv);
						for (int i = 0; i < Nv; i++) { x[i][j] = x_tmp[i]; }
					}
				}
				if (MOD_DIM == "Spatial")
				{
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
					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < 2 * Nr; j++) { y_mod_real_temp_omp[i][j] = y_mod_real[j][i]; }
					}
					// for KBest Detection
					Eigen::HouseholderQR<Eigen::MatrixXf> qr_Hreal(H_real);
					Q = qr_Hreal.householderQ();
					R = qr_Hreal.matrixQR().triangularView<Eigen::Upper>();
					//omp_set_num_threads(8);

#pragma omp parallel for private(z)
					for (int j = 0; j < Nh; j++)
					{
						MatrixXf received_signal_omp(2 * Nr, 1);
						int k_kbest_extend = pow(2, ModType / 2) * K_KBEST;
						float** paths_kbest = new float* [k_kbest_extend];
						for (int i = 0; i < k_kbest_extend; i++) paths_kbest[i] = new float[2 * Nt];
						float* ED = new float[K_KBEST];

						if (MIMODET == "MMSE") { MMSEdetection(H_real, Identity, received_signal, mmse_matrix, output_signal, HTH, y_mod_real_temp_omp[j], y_mmse_temp_omp[j], sigma_sq); }
						// ִ��һЩ����

						for (int i = 0; i < 2 * Nr; i++) { received_signal_omp(i, 0) = y_mod_real_temp_omp[j][i]; }
						z = Q.transpose() * received_signal_omp;
						if (MIMODET == "KBEST")
						{
							KBestDetection(R, z, H_real, y_mod_real_temp_omp[j], ED, y_mmse_temp_omp[j], paths_kbest, K_KBEST, ModType, Nt, Nr);
							llr_kbest_bit(Nt, y_mmse_temp_omp[j], oriLLR_kbest[j], sigma_sq, ModType, K_KBEST, paths_kbest, ED);
						}
						for (int i_path = 0; i_path < k_kbest_extend; i_path++) { delete[] paths_kbest[i_path]; }

						delete[] paths_kbest;
						delete[] ED;
					}

					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < 2 * Nt; j++) { y_mmse[j][i] = y_mmse_temp_omp[i][j]; }
					}
					// sym llr to bitwise llr
					// Detection is done for every column, that is, for every Nv bits
					// Interleaving is not currently considered
					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < 2 * Nt; j++) { y_mmse_tmp[j] = y_mmse[j][i]; }
						if (MIMODET == "MMSE") { llr_mmse_bit(Nt, y_mmse_tmp, oriLLR_tmp, sigma_sq, ModType); }
						if (MIMODET == "KBEST") { for (int j = 0; j < len_ori_llr_tmp; j++) { oriLLR_tmp[j] = oriLLR_kbest[i][j]; } }
						if (flag_interleave)
						{
							randDeintrlv(intrlv_pattern, oriLLR_tmp, Nv);
						}
						bool flagIsMax = false;
						for (int j = 0; j < Nv; j++)
						{
							if (!(j % 2))
							{
								if (abs(x_mmse_out[j][i]) > abs(x_mmse_out[j][i + 1])) { flagIsMax = true; }
								else { flagIsMax = false; }
							}
							else
							{
								if (abs(x_mmse_out[j][i]) > abs(x_mmse_out[j][i - 1])) { flagIsMax = true; }
								else { flagIsMax = false; }
							}
							oriLLRh[j][i] = oriLLR_tmp[j];
							upLLR[j][i] = oriLLR_tmp[j];
							// TEST�� det error
							x_mmse_out[j][i] = 0;
							if (oriLLR_tmp[j] < 0) { x_mmse_out[j][i] = 1; }
							if (x_mmse_out[j][i] != x[j][i]) {
								n_det_err++;
								// error at max-reliable bit. Warning: Only fit for 16QAM modulation
								if (flagIsMax) { n_det_err_mrb++; }
							}
							total_bit_cnt++;
						}
					}
					// test : detection TIME
					/*det_end = std::chrono::system_clock::now();
					auto duration_detection = std::chrono::duration_cast<std::chrono::microseconds>(det_end - start);
					file_det << duration_detection.count() << endl;*/
				}
				// Temporal modulation
				else if (MOD_DIM == "Temporal")
				{
					//modulation
					for (int i = 0; i < Nv; i++) // every row
					{
						modulation(x[i], x_mod[i], ModType, mimo_blk_length);
					}
					//AWGN channel (received y_mod)
					H = generateChannelMatrix(Nr, Nt);
					for (int j = 0; j < Nsym; j++)
					{
						// Extracting a single column
						for (int i = 0; i < Nv; i++) { x_mod_tmp[i] = x_mod[i][j]; }
						AWGN(Nr, Nt, x_mod_tmp, y_mod_tmp, H, sigma_sq);
						for (int i = 0; i < Nr; i++) { y_mod[i][j] = y_mod_tmp[i]; }
					}
					convertHToReal(H, H_real);

					for (int j = 0; j < Nsym; j++)
					{
						for (int i = 0; i < Nr; i++) { y_mod_tmp[i] = y_mod[i][j]; }
						convertYToReal(Nr, y_mod_tmp, y_mod_real_tmp);
						for (int i = 0; i < 2 * Nr; i++) { y_mod_real[i][j] = y_mod_real_tmp[i]; }
					}
					//MMSE detection done in paralell
					for (int i = 0; i < Nsym; i++)
					{
						for (int j = 0; j < 2 * Nr; j++) { y_mod_real_temp_omp[i][j] = y_mod_real[j][i]; }
					}
					Eigen::HouseholderQR<Eigen::MatrixXf> qr_Hreal(H_real);
					Q = qr_Hreal.householderQ();
					R = qr_Hreal.matrixQR().triangularView<Eigen::Upper>();

					//omp_set_num_threads(8);
					//cout << 1 << endl;
					//omp_set_num_threads(4);
					/*auto start = std::chrono::system_clock::now();*/
					/*MatrixXf H_real_trans(2*Nt, 2*Nr);
					for (int i = 0; i < 2 * Nt; i++)
						for (int j = 0; j < 2 * Nr; j++)
							H_real_trans(i, j) = H_real(j, i);
					HTHIinv = (H_real_trans * H_real + Identity).inverse();*/
# pragma omp parallel for private(z)
					for (int j = 0; j < Nsym; j++)

						//for (int x=0; x<Nsym; x++)
					{
						// int j = x % Nsym;
						MatrixXf received_signal_omp(2 * Nr, 1);
						int k_kbest_extend = pow(2, ModType / 2) * K_KBEST;
						float** paths_kbest = new float* [k_kbest_extend];
						for (int i = 0; i < k_kbest_extend; i++) { paths_kbest[i] = new float[2 * Nt]; }
						float* ED = new float[K_KBEST];

						if (MIMODET == "MMSE")
						{
							MMSEdetection(H_real, Identity, received_signal, mmse_matrix, output_signal, HTH, y_mod_real_temp_omp[j], y_mmse_temp_omp[j], sigma_sq);
							//MMSEdetection_noInverse(H_real, Identity, received_signal, mmse_matrix, output_signal, HTH, HTHIinv, y_mod_real_temp_omp[j], y_mmse_temp_omp[j], sigma_sq);
						}
						if (MIMODET == "KBEST")
						{
							for (int i = 0; i < 2 * Nr; i++) { received_signal_omp(i, 0) = y_mod_real_temp_omp[j][i]; }
							z = Q.transpose() * received_signal_omp;
							KBestDetection(R, z, H_real, y_mod_real_temp_omp[j], ED, y_mmse_temp_omp[j], paths_kbest, K_KBEST, ModType, Nt, Nr);
							llr_kbest_bit(Nt, y_mmse_temp_omp[j], oriLLR_kbest[j], sigma_sq, ModType, K_KBEST, paths_kbest, ED);
						}
						for (int i_path = 0; i_path < k_kbest_extend; i_path++) { delete[] paths_kbest[i_path]; }
						delete[] paths_kbest;
						delete[] ED;
					}
					// test: detection time
					/*auto detection_over = std::chrono::system_clock::now();
					auto duration_detection = std::chrono::duration_cast<std::chrono::microseconds>(detection_over - start);
					cout << endl;
					cout <<"Detection duration" << duration_detection.count() << "us" << endl;*/
					// test: detection time
					for (int i = 0; i < Nsym; i++)
					{
						for (int j = 0; j < 2 * Nt; j++) { y_mmse[j][i] = y_mmse_temp_omp[i][j]; }
					}
					//To bit domain output
					for (int j = 0; j < Nsym; j++)
					{
						for (int i = 0; i < 2 * Nt; i++) { y_mmse_tmp[i] = y_mmse[i][j]; }
						// With Spatial Modulation, Nt = Nv
						if (MIMODET == "MMSE") { llr_mmse_bit(Nt, y_mmse_tmp, oriLLR_tmp, sigma_sq, ModType); }
						if (MIMODET == "KBEST") { for (int i = 0; i < len_ori_llr_tmp; i++) { oriLLR_tmp[i] = oriLLR_kbest[j][i]; } }
						int i_bit = 0;
						int bit_col = 0;
						for (int j_inner = 0; j_inner < Nv; j_inner++)
						{
							for (int i_bit_col = 0; i_bit_col < ModType; i_bit_col++)
							{
								i_bit = ModType * j_inner + i_bit_col;
								bit_col = j * ModType + i_bit_col;
								oriLLRh[j_inner][bit_col] = oriLLR_tmp[i_bit];
								upLLR[j_inner][bit_col] = oriLLR_tmp[i_bit];
							}

						}
					}
					// test: LLR transfer time
					/*auto llr_to_bit_over = std::chrono::system_clock::now();
					auto duration_llr = std::chrono::duration_cast<std::chrono::microseconds>(llr_to_bit_over - detection_over);
					cout << endl;
					cout<< "LLR transfer duration" << duration_llr.count() << "us" << endl;*/
					// test: LLR transfer time
					if (flag_interleave)
					{
						for (int j = 0; j < Nh; j++)
						{
							for (int i = 0; i < Nv; i++) { oriLLR_tmp[i] = oriLLRh[i][j]; }
							randDeintrlv(intrlv_pattern_2D[j], oriLLR_tmp, Nv);
							for (int i = 0; i < Nv; i++)
							{
								oriLLRh[i][j] = oriLLR_tmp[i];
								upLLR[i][j] = oriLLR_tmp[i];
							}

						}
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
			for (int iter = 1; iter <= softiter; iter++)
			{
				//vertical decoding
		   // test: decoding time
			//auto decoding_start = std::chrono::system_clock::now();

			// test: decoding time
#pragma omp parallel for
				for (int decodingi = 0; decodingi < Nh; decodingi++) // 291.70kb
				{

					//for SCL Decoding
					int Counth = 0, Countinfoh = 0;
					float Qsumunh = 0;
					int** sum = new int* [L];
					//for (int i = 0; i < L; i++)  sum[i] = new  int[Nmax];
					int** u1 = new int* [L];
					for (int i = 0; i < L; i++)  u1[i] = new _MM_ALIGN16 int[Nmax];
					for (int i = 0; i < L; i++)  sum[i] = new _MM_ALIGN16 int[Nmax];
					//for (int i = 0; i < L; i++)  u1[i] = new  int[Nmax];

					float* LLR_in = new float[L];
					float* W = new float[2 * L];
					int* SCL_index = new int[L];
					int* better = new int[L];
					int* worse = new int[L];
					int* indexpm = new int[L];
					float* PM = new float[L];
					float** LLR = new float* [L];
					for (int i = 0; i < L; i++)  LLR[i] = new _MM_ALIGN16 float[2 * Nmax + 2];
					// for (int i = 0; i < L; i++)  LLR[i] = new float[2 * Nmax + 2];
					////////////////////////////////
					for (int j = 0; j < Nv; j++)
						for (int k = 0; k < L; k++)
							LLR[k][j] = upLLR[j][decodingi];
					//cout << LLR[0][j] << " ";

				//cout << endl;
					int index = 0;
					/*			if (L == 1)
								{
									//T_start = clock();
									Count = 0;
									decode(*LLR, *sum, (int)(log(Nv) / log(2)), A_Acv, flag, Nv, Kv);
									//PolarEncode_xor(*u1, *sum, N);
									//T_finish = clock();
								}
								else
								{*/
					for (int k = 0; k < L; k++) { PM[k] = 0; indexpm[k] = k; }
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
					if (adaptive_beta)
					{
						if (softmethod == "Pyndiah") Pyndiahv_adaptivebeta(LLR, sum, oriLLRh, indexpm, iter * 2 - 2, decodingi, upLLR);
						if (softmethod == "MITSO") MITSOv(LLR, PM, sum, oriLLRh, indexpm, iter * 2 - 2, decodingi, upLLR, Qsumunh);
					}
					else
					{
						if (softmethod == "Pyndiah") Pyndiahv(LLR, sum, oriLLRh, indexpm, iter * 2 - 2, decodingi, upLLR);
						if (softmethod == "MITSO") MITSOv(LLR, PM, sum, oriLLRh, indexpm, iter * 2 - 2, decodingi, upLLR, Qsumunh);
					}
					// half iteration Enabled
					if (half_iter && iter == softiter)
					{
						if (iter == softiter) PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), Nh);
						for (int j = 0; j < Kv; j++)
							umid[Av[j]][decodingi] = u1[index][Av[j]];
					}
					// half iteration
					////////delete SCL variables
					for (int i = 0; i < L; i++)
						delete[]LLR[i];
					delete[] LLR;
					for (int i = 0; i < L; i++)
						delete[]sum[i];
					for (int i = 0; i < L; i++)
						delete[]u1[i];
					delete[] u1;
					delete[] sum;
					delete[] PM;
					delete[] LLR_in;
					delete[] W;
					delete[] SCL_index;
					delete[] better;
					delete[] worse;
					delete[] indexpm;

				}
				//cout << "Horizon Start:\n";
				//horizontal decoding
				// #
				if (weirditer) { continue; }
				/*auto decoding_over = std::chrono::system_clock::now();
				auto duration_decoding = std::chrono::duration_cast<std::chrono::microseconds>(decoding_over - decoding_start);
				cout << endl;
				cout << "Decoding duration" << duration_decoding.count() << "us" << endl;*/

#pragma omp parallel for
				for (int decodingi = 0; decodingi < Nv; decodingi++)
				{
					// 293.70kb
					// cout << "KKKK" << endl;
					if (half_iter && (iter == softiter)) { break; }
					string polarde_now = polarde_array[iter - 1];
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

					/*			if (L == 1)
								{
									//T_start = clock();
									Count = 0;
									decode(*LLR, *sum, (int)(log(Nh) / log(2)), A_Ach, flag, Nh, Kh);
									if (iter == softiter) PolarEncode_xor(*u1, *sum, Nh);
									//T_finish = clock();
								}
								else
								{*/
								//for (int k = 0; k < L; k++)  PM[k] = 0, indexpm[k] = k; 
								////T_start = clock();
								//decode_list(LLR, sum, (int)(log(Nh) / log(2)), A_Ach, 0, Nh, Kh, PM, L, (int)(log(L) / log(2)), LLR_in, W, SCL_index, better, worse, &Countv, &Countinfov, &Qsumunv);
								////T_finish = clock();
								////Decoding_Duration += (T_finish - T_start);
								//for (int i = 0; i < L - 1; i++)
								//	for (int j = i + 1; j < L; j++)
								//		if (PM[indexpm[i]] < PM[indexpm[j]])
								//		{
								//			int ttt = indexpm[i];
								//			indexpm[i] = indexpm[j];
								//			indexpm[j] = ttt;
								//		}


							/*
							* Encoding: Info -> Horizontal Encoding -> A -> set to 0 -> B -> Horizontal Encoding -> Vertical Encoding
							* Decoding: Received Info -> Horizontal Decoding -> Vertical Decoding -> Horizontal Decoding
							*/
					if (polarde_now == "SLD")
					{
						LSD_decode_outall(u_decod, u1_LSD, LLR[0], GGh, A_Ach, Nh, L_LSD, n_paths);
						// if (iter == softiter) PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), Nh);
						// for soft information exchange, inverse coding:
						// THIS LINE is WRONG! BUT IT SEEMS BETTER THAN THE RIGHT VERSION
						for (int i = 0; i < L_LSD; i++) { PolarEncode_xor(sum_LSD[i], u1_LSD[i], Nh); }
						// for (int i = 0; i < L_LSD; i++) { PolarEncode_xor(sum_LSD[i], u1_LSD[i], Nh); }
						// if (iter == softiter) PolarEncode_xor(u1[0], sum_LSD[0], Nh);
						index = indexpm[0];
						// test: after horizontal decoding
						/*if (decodingi == Av[0])
						{
							cout << "Direct output of horizontal decoding:" << endl;
							for (int i = 0; i < Nh; i++) { cout << u1_LSD[0][i] << " "; }
							cout << endl;
							cout << "horizontal decoder out, encode once:" << endl;
							for (int i = 0; i < Nh; i++) { cout << sum_LSD[0][i] << " "; }
							cout << endl;
						}*/
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
							if (adaptive_beta)
							{
								if (softmethod == "Pyndiah") Pyndiahh_LSD_adaptivebeta(LLR, sum_LSD, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR);
								if (softmethod == "MITSO") MITSOh(LLR, PM, sum_LSD, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR, Qsumunv);
							}
							else
							{
								// cout << "IN LSD!!!" << endl;
								if (softmethod == "Pyndiah")
									Pyndiahh_LSD_getsum(LLR, sum_LSD, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR, sum_sdis_array[decodingi], GGh, Ach);
								//if (softmethod == "Pyndiah") Pyndiahh_LSD(LLR, sum_LSD, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR);
								if (softmethod == "MITSO") MITSOh(LLR, PM, sum_LSD, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR, Qsumunv);
							}
							//	getchar();
						}
						else
							for (int j = 0; j < Kh; j++)
								umid[decodingi][Ah[j]] = u1_LSD[index][Ah[j]];
					}
					if (polarde_now == "SCL")
					{
						for (int k = 0; k < L; k++)  PM[k] = 0, indexpm[k] = k;
						//T_start = clock();
						decode_list(LLR, sum, (int)(log(Nh) / log(2)), A_Ach, 0, Nh, Kh, PM, L, (int)(log(L) / log(2)), LLR_in, W, SCL_index, better, worse, &Countv, &Countinfov, &Qsumunv);
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
						if (iter == softiter) PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), Nh);
						index = indexpm[0];
						//}
						//calc distance s()
						if (iter < softiter)
						{
							if (adaptive_beta)
							{
								if (softmethod == "Pyndiah") Pyndiahh_adaptivebeta(LLR, sum, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR);
								if (softmethod == "MITSO") MITSOh(LLR, PM, sum, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR, Qsumunv);
								//	getchar();
							}
							else
							{
								//if (softmethod == "Pyndiah") Pyndiahh(LLR, sum, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR);
								if (softmethod == "Pyndiah") Pyndiahh_getsum(LLR, sum, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR, sum_sdis_array[decodingi]);
								if (softmethod == "MITSO") MITSOh(LLR, PM, sum, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR, Qsumunv);
								//	getchar();
							}
						}
						else
							for (int j = 0; j < Kh; j++)
								umid[decodingi][Ah[j]] = u1[index][Ah[j]];
					}

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
			int calcsum[Nv], calcu1[Nv];
			int calcu1h[Nh];
			// calcsum is from the direct output of the decoder
			// if output is the result of row decoding (meaning there is another column decoding to be done)
			if (half_iter)
			{
				for (int i = 0; i < Kv; i++)
				{
					PolarEncode_xor(calcu1h, umid[Av[i]], Nh);
					for (int j = 0; j < Kh; j++)
						uout[Av[i]][Ah[j]] = calcu1h[Ah[j]];
				}
			}
			// if output is the result of column decoding (meaning that there is another row decoding to be done)
			else
			{
				for (int i = 0; i < Kh; i++)
				{
					for (int j = 0; j < Nv; j++)
						calcsum[j] = umid[j][Ah[i]];
					PolarEncode_xor(calcu1, calcsum, Nv);
					for (int j = 0; j < Kv; j++)
						uout[Av[j]][Ah[i]] = calcu1[Av[j]];
				}
			}
			// test: final output
			/*cout << "final output" << endl;
			for (int k = 0; k < Nh; k++) { cout << uout[Av[0]][k] << " "; }
			cout << endl;*/
			// test: final output
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
				/*printf("\r<%d> SNR= %.2f FER= %d / %d = %f BER= %d / %d = %lf E-Det = %f, E-MRB = %f"  ,
					frame, nSNR, Num_Frame_Error, frame, (double)Num_Frame_Error / frame, Num_Error, (N_Info * frame), (double)Num_Error / N_Info / frame,
					(double)n_det_err/total_bit_cnt, (double)n_det_err_mrb/(total_bit_cnt/(ModType/2)));*/
				int L_sdis_sum = L_LSD;
				sdis_sum = 0;
				for (int i = 0; i < Nv; i++) sdis_sum += sum_sdis_array[i];
				if (polarde_array[0] == "SCL") L_sdis_sum = L;
				// cout << "sdis_sum_avg" << sdis_sum / frame / L_sdis_sum << endl;
				printf("\r<%d> SNR= %.2f FER= %d / %d = %f BER= %d / %d = %lf ,sdis_avg = %f",
					frame, nSNR, Num_Frame_Error, frame, (double)Num_Frame_Error / frame, Num_Error, (N_Info * frame), (double)Num_Error / N_Info / frame, (double)sdis_sum / frame / L_sdis_sum);
				fflush(stdout);
			}
			//printf("%d %d", Num_Frame_Error, Num_Error);
			//getchar();
			//Decoding_Duration += (T_finish - T_start);
			/*auto dec_end = std::chrono::system_clock::now();
			auto duration_decoding = std::chrono::duration_cast<std::chrono::microseconds>(dec_end - det_end);
			auto duration_total = std::chrono::duration_cast<std::chrono::microseconds>(dec_end - start);
			file_dec << duration_decoding.count() << endl;
			file_total << duration_total.count() << endl;*/
			if (frame % 16 == 0)
			{
				file_err_num << Num_Error << endl;
			}

		}
		//Decoding_Duration = Decoding_Duration * 100 / CLOCKS_PER_SEC / frame;
		//Decoding_Duration = Decoding_Duration / CLOCKS_PER_SEC / Max_frame;

		//printf("%f  %f  %f  %f\n", nSNR, Decoding_Duration, (double)Num_Error / N_Info / Max_frame, (double)Num_Frame_Error / Max_frame);
		printf("%f %f %f\n", nSNR, (double)Num_Frame_Error / frame, (double)Num_Error / N_Info / frame);
		BER_array[SNR_count] = (double)Num_Error / N_Info / frame;
		FER_array[SNR_count] = (double)Num_Frame_Error / frame;
		// save to .txt file
		string polarde_serial = "";
		for (int i = 0; i < softiter; i++)
			polarde_serial = polarde_serial + polarde_array[i] + string("_");
		string file_name = string("ResFinal\\2D_LSD_D_") + polarde_serial + MIMODET + string("T") + to_string(Nt) + string("R") + to_string(Nr)
			+ string("Q") + to_string(int(pow(2, ModType))) + string("L") + to_string(Nh) + string("K") + to_string(Kh)
			+ string("T") + to_string(softiter) + string("H") + to_string(half_iter) + string("UB") + to_string(usebeta) + string("I") + to_string(flag_interleave)
			+ string(".txt");
		save_array_txt(nSNR, BER_array[SNR_count], FER_array[SNR_count], file_name);
		int L_sdis_sum = L_LSD;
		sdis_sum = 0;
		for (int i = 0; i < Nv; i++) sdis_sum += sum_sdis_array[i];
		if (polarde_array[0] == "SCL") L_sdis_sum = L;
		cout << "sdis_sum_avg" << sdis_sum / frame / L_sdis_sum << endl;
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
	for (int i = 0; i < Nh; i++)
		delete[] oriLLR_kbest[i];
	delete[] oriLLR_kbest;
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
	for (int i = 0; i < mimo_blk_height; i++) {
		delete[] x_mod[i];
	}
	delete[] x_mod;
	for (int i = 0; i < Nr; i++) {
		delete[] y_mod[i];
	}
	delete[] y_mod;
	for (int i = 0; i < 2 * mimo_blk_height; i++) {
		delete[] y_mmse[i];
	}
	delete[] y_mmse;
	for (int i = 0; i < 2 * Nr; i++) {
		delete[] y_mod_real[i];
	}
	for (int i = 0; i < mimo_blk_length; i++)
	{
		delete[] y_mod_real_temp_omp[i];
	}
	delete[] y_mod_real_temp_omp;
	for (int i = 0; i < mimo_blk_length; i++)
	{
		delete[] y_mmse_temp_omp[i];
	}
	delete[] y_mmse_temp_omp;
	delete[] y_mod_real;
	delete[] x_mod_tmp;
	delete[] y_mod_tmp;
	delete[] y_mod_real_tmp;
	delete[] x_tmp;
	delete[] oriLLR_tmp;
	delete[] sum_sdis_array;
}



void system2D_LSD_epdet(double* BER_array, double* FER_array)
{
	std::ofstream file_total("total_timing.txt", std::ios::app);
	std::ofstream file_det("det_timing.txt", std::ios::app);
	std::ofstream file_dec("dec_timing.txt", std::ios::app);
	std::ofstream file_err_num("total_frames.txt", std::ios::app);
	// Lightweight progress log: one line every 10 newly-errored frames.
	// Filename encodes softiter and polarconh so multiple parallel runs do not collide.
	// Written to current working directory (caller is expected to launch the exe from
	// the desired output folder, e.g. result0512/run_<polarconh>_iter<softiter>/).
	std::string progress_log_name = std::string("progress_iter") + std::to_string(softiter) + "_" + polarconh + ".txt";
	std::ofstream file_progress(progress_log_name, std::ios::app);
	if (file_progress.is_open()) {
		file_progress << "# === run start: polarconh=" << polarconh << " polarconv=" << polarconv
			<< " softiter=" << softiter << " ===\n";
		file_progress << "# fields: SNR_dB  frames  err_frames  FER  err_bits  total_bits  BER\n";
		file_progress.flush();
	}
	auto det_end = std::chrono::system_clock::now();
	int Nt2 = 2 * Nt;
	int Nr2 = 2 * Nr;

	/*if (!file.is_open()) {
		std::cerr << "Error opening file: " << file_name << std::endl;
		return;
	}*/

	//file << setw(20) << SNR << setw(20) << BER << setw(20) << FER << "\n";

	//freopen("result.txt", "w", stdout);
	printsetting();
	cout << "\nSNR	 " << "FER      BER" << endl;
	int N_Info = Kv * Kh - CRCL, Nmax;
	//Interleaving pattern
	/*std::vector<int> intrlv_pattern(Nv);
	std::vector<int>* intrlv_pattern_2D = new vector<int>[Nh];
	for (int i = 0; i < Nh; i++) intrlv_pattern_2D[i].resize(Nv);*/
	int intrlv_pattern_num = Nh / ModType;
	std::vector<int> intrlv_pattern(Nv);
	std::vector<int>* intrlv_pattern_2D = new vector<int>[Nh];
	std::vector<int>* intrlv_pattern_2D_time = new vector<int>[Nv];
	std::vector<int>* intrlv_pattern_2D_inner = new vector<int>[Nh];
	std::vector<int>** intrlv_pattern_2D_time_partial = new vector<int> *[Nv];
	for (int i = 0; i < Nv; i++) intrlv_pattern_2D_time_partial[i] = new vector<int>[intrlv_pattern_num];
	for (int i = 0; i < Nv; i++)
	{
		for (int j = 0; j < intrlv_pattern_num; j++)
		{
			intrlv_pattern_2D_time_partial[i][j].resize(ModType);
		}
	}
	for (int i = 0; i < Nh; i++) intrlv_pattern_2D[i].resize(Nv);
	for (int i = 0; i < Nh; i++) intrlv_pattern_2D_inner[i].resize(Kv);
	for (int i = 0; i < Nv; i++) intrlv_pattern_2D_time[i].resize(Nh);

	float CodeRate = (float)N_Info / (Nh * Nv);
	Nmax = max(Nh, Nv);
	// num of symbols per vertical codeword
	int Nsym = 0;
	if (MOD_DIM == "Spatial") { Nsym = Nv / ModType; }
	else if (MOD_DIM == "Temporal") { Nsym = Nh / ModType; }
	int symbol_length = ModType;

	// Size of a mimo block is given by height(spatial) * length(temporal)
	int mimo_blk_height = 0;
	int mimo_blk_length = 0;
	if (MOD_DIM == "Spatial")
	{
		mimo_blk_height = Nsym;
		mimo_blk_length = Nh;
	}
	else if (MOD_DIM == "Temporal")
	{
		mimo_blk_height = Nv;
		mimo_blk_length = Nsym;
	}
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
	int** x_parity_check = new int* [Nv];
	for (int i = 0; i < Nh; i++) x_parity_check[i] = new int[Nh];

	// Arrays for MIMO
	int** x_mmse_out = new int* [Nv]; // delete!
	for (int i = 0; i < Nv; i++) { x_mmse_out[i] = new int[Nh]; }

	//std::complex<float>** x_mod = new std::complex<float>* [Nsym]; // modulated symbols (spatial modulation) DELETE!
	//for (int i = 0; i < Nsym; i++) { x_mod[i] = new std::complex<float>[Nh]; }
	//std::complex<float>** x_mod_temporal = new std::complex<float>*[Nv]; // modulated symbols (temporal modulation) DELETE!
	//for (int i = 0; i < Nv; i++) x_mod_temporal[i] = new std::complex<float>[Nsym];
	//std::complex<float>** y_mod = new std::complex<float>* [Nr]; // Output symbols after AWGN channels // delete!
	//for (int i = 0; i < Nr; i++) { y_mod[i] = new std::complex<float>[Nh]; }
	//float** y_mod_real = new float* [2 * Nr];
	//for (int i = 0; i < 2*Nr; i++) y_mod_real[i] = new float[Nh];
	//float** y = new float* [Nv];// Output Signal After AWGN Channel
	//for (int i = 0; i < Nv; i++) { y[i] = new float[Nh]; }
	//float** y_mmse = new float* [2*Nsym];
	//for (int i = 0; i < 2*Nsym; i++) { y_mmse[i] = new float[Nh]; }  // MMSE detector output // delete!
	//float* y_mmse_tmp = new float[2*Nsym];
	//std::complex<float>* x_mod_tmp = new std::complex<float>[Nsym]; // modulated symbols of a single column // delete!
	//std::complex<float>* y_mod_tmp = new std::complex<float>[Nr]; // received symbols of a single column // delete!
	//float* y_mod_real_tmp = new float[2*Nr];
	//// float* y_mod_tmp = new float[Nr];  // delete!


	std::complex<float>** x_mod = new std::complex<float>*[mimo_blk_height]; // modulated symbols (spatial modulation) DELETE!
	for (int i = 0; i < mimo_blk_height; i++) { x_mod[i] = new std::complex<float>[mimo_blk_length]; }
	std::complex<float>** y_mod = new std::complex<float>*[Nr]; // Output symbols after AWGN channels // delete!
	for (int i = 0; i < Nr; i++) { y_mod[i] = new std::complex<float>[mimo_blk_length]; }
	float** y_mod_real = new float* [2 * Nr];
	for (int i = 0; i < 2 * Nr; i++) y_mod_real[i] = new float[mimo_blk_length];
	float** y = new float* [Nv];// Output Signal After AWGN Channel
	for (int i = 0; i < Nv; i++) { y[i] = new float[Nh]; }
	float** y_mmse = new float* [2 * mimo_blk_height];
	for (int i = 0; i < 2 * mimo_blk_height; i++) { y_mmse[i] = new float[mimo_blk_length]; }  // MMSE detector output // delete!
	float* y_mmse_tmp = new float[2 * mimo_blk_height];
	std::complex<float>* x_mod_tmp = new std::complex<float>[mimo_blk_height]; // modulated symbols of a single column // delete!
	std::complex<float>* y_mod_tmp = new std::complex<float>[Nr]; // received symbols of a single column // delete!
	float* y_mod_real_tmp = new float[2 * Nr];
	// float* y_mod_tmp = new float[Nr];  // delete!
	float** y_mod_real_temp_omp = new float* [mimo_blk_length]; // new delete!
	for (int i = 0; i < mimo_blk_length; i++) { y_mod_real_temp_omp[i] = new float[2 * Nr]; }
	float** y_mmse_temp_omp = new float* [mimo_blk_length]; // new delete!
	for (int i = 0; i < mimo_blk_length; i++) { y_mmse_temp_omp[i] = new float[2 * mimo_blk_height]; }
	// ep
	float** ep_extr_mean_omp = new float* [mimo_blk_length];
	for (int i = 0; i < mimo_blk_length; i++) { ep_extr_mean_omp[i] = new float[2 * mimo_blk_height]; }
	float** ep_extr_var_omp = new float* [mimo_blk_length];
	for (int i = 0; i < mimo_blk_length; i++) { ep_extr_var_omp[i] = new float[2 * mimo_blk_height]; }
	int* x_tmp = new int[Nv]; // delete!

	// paths selected by k-best detection


	// CRC bits
	unsigned char* u_crc = new unsigned char[Kv * Kh];

	float** upLLR = new float* [Nv];
	for (int i = 0; i < Nv; i++) { upLLR[i] = new _MM_ALIGN16 float[Nh]; }
	//float** oriLLRv = new float* [Nh];
	//for (int i = 0; i < Nh; i++) { oriLLRv[i] = new _MM_ALIGN16 float[2 * Nv + 2]; }
	float** oriLLRh = new float* [Nv];
	for (int i = 0; i < Nv; i++) { oriLLRh[i] = new _MM_ALIGN16 float[Nh]; }
	int len_ori_llr_tmp = 0;
	if (MOD_DIM == "Spatial") { len_ori_llr_tmp = Nv; }
	else if (MOD_DIM == "Temporal") { len_ori_llr_tmp = ModType * Nt; }
	float* oriLLR_tmp = new _MM_ALIGN16 float[len_ori_llr_tmp];  // delete!
	//float** oriLLR_kbest = new float* [Nh];
	float** oriLLR_kbest = new float* [mimo_blk_length];
	for (int i = 0; i < Nh; i++) { oriLLR_kbest[i] = new float[len_ori_llr_tmp]; }
	double* sum_sdis_array = new double[Nv];
	for (int i = 0; i < Nv; i++) sum_sdis_array[i] = 0;
	int* n_tried_patterns_h = new int[softiter];
	for (int i = 0; i < softiter; i++) n_tried_patterns_h[i] = 0;
	int* n_tried_patterns_v = new int[softiter];
	for (int i = 0; i < softiter; i++) n_tried_patterns_v[i] = 0;

	//*****************************************************************
	// EP detector
	float* Hreal_array = new float[(2 * Nr) * (2 * Nt)]; // delete
	float* HTH_array = new float[(2 * Nt) * (2 * Nt)]; // delete
	float* HTY_array = new float[2 * Nt]; // delete 
	float** channel_scale_factor = new float* [Nr];
	for (int i = 0; i < Nr; i++) { channel_scale_factor[i] = new float[Nt]; for (int j = 0; j < Nt; j++) { channel_scale_factor[i][j] = 0.0; } }
		
	if (non_uniform_channel)
	{
		getChannelScalingFactor(Nr, Nt, diff_gain, nullptr, channel_scale_factor, 1);
		/*if (predefined_non_uniform) { getChannelScalingFactor(Nr, Nt, diff_gain, use_short_codeword, channel_scale_factor); }
		else { getChannelScalingFactor(Nr, Nt, diff_gain, use_short_codeword, channel_scale_factor); }*/
	}



	int** GGh = new int* [Nh];
	int** GGv = new int* [Nv];


	MatrixXf Gh = Fn(Nh);  // Generate Matrix of horizontal coding
	MatrixXf Gv = Fn(Nv);  // Generate Matrix of vertical coding
	MatrixXcf H(Nr, Nt); // Complex channel Matrix H: Nr x Nt
	//Matrix<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic> H = MatrixXcf::Random(Nr, Nt); // Real domain channel matrix
	//// //MatrixX Definitions
	MatrixXf H_real(2 * Nr, 2 * Nt); // Real domain channel matrix
	MatrixXf received_signal(2 * Nr, 1);
	MatrixXf HTH(2 * Nt, 2 * Nt);
	MatrixXf Identity = MatrixXf::Identity(2 * Nt, 2 * Nt);
	MatrixXf output_signal(2 * Nt, 1);
	MatrixXf mmse_matrix(2 * Nt, 2 * Nr);
	MatrixXf HTHIinv(2 * Nt, 2 * Nt);

	// for KBest
	MatrixXf R(2 * Nt, 2 * Nt);
	MatrixXf Q(2 * Nr, 2 * Nt);
	MatrixXf z(2 * Nt, 1);



	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> H_real = MatrixXf::Random(2 * Nr, 2 * Nt); // Real domain channel matrix
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> received_signal = MatrixXf::Random(2 * Nr,1);
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> HTH = MatrixXf::Random(2 * Nt, 2 * Nt);
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> Identity = MatrixXf::Random(2 * Nt, 2 * Nt);
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> output_signal = MatrixXf::Random(2 * Nt, 1);
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> mmse_matrix = MatrixXf::Random(2 * Nt, 2 * Nr);


	// Matrix Definitions
	//Eigen::Matrix<float, 2 * Nr, 2 * Nt> H_real; // ʵ���ŵ�����
	//Eigen::Matrix<float, 2 * Nr, 1> received_signal; // �����źž���
	//Eigen::Matrix<float, 2 * Nt, 2 * Nt> HTH; // HTH����
	//Eigen::Matrix<float, 2 * Nt, 2 * Nt> Identity = Eigen::Matrix<float, 2 * Nt, 2 * Nt>::Identity(); // ��λ����
	//Eigen::Matrix<float, 2 * Nt, 1> output_signal; // ����źž���
	//Eigen::Matrix<float, 2 * Nt, 2 * Nr> mmse_matrix; // MMSE����


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
	// srand((unsigned int)time(0));
	std::random_device rd;              // ��ȡ���������
	std::mt19937 gen(rd());             // ʹ�����ӳ�ʼ��mt19937������
	std::uniform_int_distribution<> distrib(1, 1000000);
	int Num_Error, Num_Frame_Error;
	float Decoding_Duration;

	int Num_Error_new, Num_Frame_Error_new;
	float Decoding_Duration_new;

	clock_t T_start, T_finish;
	int SNR_count = 0;
	int anssc, anssd;
	for (double nSNR = Min_SNR; nSNR <= Max_SNR; nSNR += SNR_step)
	{
		for (int i = 0; i < Nv; i++) sum_sdis_array[i] = 0;
		double sdis_sum = 0;
		// TEST: det err
		int n_det_err = 0; // n error bits at detector output
		int n_det_err_mrb = 0; // error bits at detector output(MRBs)
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


		// ENUMERATE CODEWORD
		while (ENUMERATE_CODEWORD) {
			const int N = Nh;
			const int K = Kh;
			const int n_irregular = 10;
			// int whole_Ah[26] = { 3, 5 6 7 9 10 11 12 13 14 15 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 };
			/*int idx_irr[n_irregular] = { 3,5,6,7,9,10,11,12,13,14,17,18,19,20,21,22,23,24,25,26 };*/
			int idx_irr[n_irregular] = { 3,5,6,7,9,10,11,12,13,14};
			// int idx_irr[n_irregular] = { 0,1,2,4,8,16,21,22,23,24,25,26,27,28,29,30,31 };
			int u[N];
			int cnt_nnz_pos[N];
			for (int i = 0; i < N; i++) { cnt_nnz_pos[i] = 0; }
			int codeword[N];
			for (int i = 0; i < Kh; i++) { cout << Ah[i] << " "; }
			// ��С�������ִ洢
			const int MAX_MIN_CODWORDS = 1E5; // �������1000����С��������
			int** min_codewords = new int* [MAX_MIN_CODWORDS];
			for (int i = 0; i < MAX_MIN_CODWORDS; i++) {
				min_codewords[i] = new int[N];
			}
			// int min_codewords[MAX_MIN_CODWORDS][N];
			int min_weight = N + 1; // ��ʼ��Ϊ����������+1
			int min_count = 0;

			// ö������2^K����Ϣ����
			long long total_combinations = 1LL << K; // 2^K

			for (long long info = 1; info < total_combinations; info++) {
				// ��ʼ����������������λ��Ϊ0��
				memset(u, 0, N * sizeof(int));

				// ������Ϣλ
				for (int i = 0; i < K; i++) {
					int bit = (info >> i) & 1; // ��ȡ��iλ
					u[Ah[i]] = bit;
				}

				// ��������
				PolarEncode_xor(codeword, u, N);

				// ���㵱ǰ��������
				// int current_weight = CalculateWeight(codeword, N);
				int weight = 0;
				for (int i = 0; i < N; i++) {
					weight += codeword[i];
				}
				int current_weight = weight;

				// ������С�������ּ���
				if (current_weight < min_weight) {
					min_weight = current_weight;
					min_count = 0;
					memcpy(min_codewords[min_count], codeword, N * sizeof(int));
					min_count++;
				}
				else if (current_weight == min_weight) {
					if (min_count < MAX_MIN_CODWORDS) {
						memcpy(min_codewords[min_count], codeword, N * sizeof(int));
						min_count++;
					}
					else {
						std::cerr << "���棺������С���ִ洢����" << std::endl;
					}
				}
			}

			// ������
			std::cout << "��С����: " << min_weight << std::endl;
			std::cout << "��С�������ָ���: " << min_count << std::endl;
			/*for (int i = 0; i < min_count; i++) {
				std::cout << "���� " << i + 1 << ": ";
				for (int j = 0; j < N; j++) {
					std::cout << min_codewords[i][j];
				}
				std::cout << std::endl;
			}*/
			int selected_min = 0;
			for (int i = 0; i < min_count; i++) {
				bool is_legal = true;
				for (int j = 0; j < n_irregular; j++) { if (min_codewords[i][idx_irr[j]] != 0) { is_legal = false; } }
				if (!is_legal) { continue; }
				for (int j = 0; j < N; j++) {
					if (min_codewords[i][j] != 0) { cnt_nnz_pos[j]++; }
				}
				selected_min++;
			}
			for (int i = 0; i < N; i++) { cout << i << ": " << cnt_nnz_pos[i] << endl; }
			cout << "Selected min codewords: " << selected_min << endl;
			// delete
			for (int i = 0; i < MAX_MIN_CODWORDS; i++) {
				delete[] min_codewords[i];
			}
			delete[] min_codewords;
		}

		// ENUMERATE CODEWORD

		/*for (int i = 0; i < Nh; i++)
		{
			cout << A_Ach[i] << ",";
		}
		cout << endl;*/
		// test: Horizontal frozen bits

		/*cout << "Frozen Bits" <<endl;
		for (int i = 0; i < Nh - Kh; i++)
		{
			cout << setw(5) << Acv[i];
		}
		cout << endl;*/
		// test: Horizontal frozen bits
		// Frame Loop
		int passnum;
		int frame = 0;

		// generate interleaving patterns
		if (flag_interleave)
		{
			string ofile_intrlv_name = string("intrlv_pattern_71_frame_49263.txt");
			string ofile_spatial_intrlv_name = string("intrlv_pattern_spatial71.txt") ;
			ifstream ofile_intrlv(ofile_intrlv_name);
			ifstream ofile_spatial_intrlv(ofile_spatial_intrlv_name);
			int file_in = 1;
			for (int i = 0; i < Nv; i++)
			{
				for (int j = 0; j < intrlv_pattern_num; j++)
				{

					for (int k = 0; k < ModType; k++)
						// ofile_intrlv >> file_in;
						ofile_intrlv >> intrlv_pattern_2D_time_partial[i][j][k];
					//cout << endl;
				}
			}
			ofile_intrlv.close();

			for (int i = 0; i < Nh; i++)
			{
				for (int j = 0; j < Nv; j++)
					ofile_spatial_intrlv>> intrlv_pattern_2D[i][j];
			}
			ofile_spatial_intrlv.close();
			/*for (int i = 0; i < Nh; i++)
			{
				for (int j = 0; j < Nv; j++)
					cout << intrlv_pattern_2D[i][j] << endl;
			}*/
			generateIntrlvPattern(intrlv_pattern, Nv);
			for (int i = 0; i < Nh; i++) { generateIntrlvPattern(intrlv_pattern_2D[i], Nv); }
			for (int i = 0; i < Nv; i++) { generateIntrlvPattern(intrlv_pattern_2D_time[i], Nh); }
			if (flag_interleave_partial)
			{
				for (int i = 0; i < Nv; i++)
				{
					for (int j = 0; j < intrlv_pattern_num; j++)
					{
						generateIntrlvPattern(intrlv_pattern_2D_time_partial[i][j], ModType);
						/*for (int k = 0; k < ModType; k++)
							cout << intrlv_pattern_2D_time_partial[i][j][k];
						cout << endl;*/
					}
				}
			}
			for (int i = 0; i < Nh; i++) { generateIntrlvPattern(intrlv_pattern_2D_inner[i], Kv); }
		}

		//ofstream ofile_intrlv("intrlv_pattern.txt");
		//for (int i=0; i<Nv; i++)
		//{
		//	for (int j = 0; j < intrlv_pattern_num; j+for+)
		//	{
		//		
		//		for (int k = 0; k < ModType; k++)
		//			ofile_intrlv << intrlv_pattern_2D_time_partial[i][j][k] << "   ";
		//		//cout << endl;
		//	}
		//}
		//ofile_intrlv.close();
		int frame_tmp = 0;
		int frame_delta = 0;
		int n_frame_error_pre = 0;
		while (Num_Frame_Error < errframe)
		{
			/*test: processing time*/
			//auto start = std::chrono::system_clock::now();
			/*if (frame == 10) {
				cout << "AAAAAAAAAA" << endl;
			}*/
			frame++;
			// generate interleaving pattern
			/*if (flag_interleave)
			{
				generateIntrlvPattern(intrlv_pattern, Nv);
				for (int i = 0; i < Nh; i++) { generateIntrlvPattern(intrlv_pattern_2D[i], Nv); }
				for (int i = 0; i < Nv; i++) { generateIntrlvPattern(intrlv_pattern_2D_time[i], Nh); }
				if (flag_interleave_partial)
				{
					for (int i = 0; i < Nv; i++)
					{
						for (int j = 0; j < intrlv_pattern_num; j++)
						{
							generateIntrlvPattern(intrlv_pattern_2D_time_partial[i][j], ModType);
						}
					}
				}
			}*/

			//cout << frame << endl;
			int u_input[16] = { 0,     1     ,0,     1,
									0,     1,     0,     0,
									0 ,    1,     0,     0,
									1 ,    0,     0,     0 };
			for (int i = 0; i < N_Info; i++)
			{
				a_data[i] = (unsigned char)(distrib(gen) % 2);
				//		a_data[i] = u_input[i];
			}

			//cout << endl;
			if (CRCL) tx_append_crc(a_data, N_Info, CRCL);
			for (int i = 0; i < Nv; i++)
				for (int j = 0; j < Nh; j++)
					u[i][j] = 0;
			//omp_set_num_threads(4);
//#pragma omp parallel for
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
			// test:original info
			//cout << "original info u:" << endl;
			////for (int j = 0; j < Nh; j++) { cout << u[Av[0]][j] << " "; }
			//for (int i = 0; i < Nv; i++)
			//{
			//	cout << endl;
			//	for (int j = 0; j < Nh; j++)
			//		cout << u[i][j];
			//}
			//cout << endl;
			// Polar Encoding
//#pragma omp parallel for

			// test: original info
			/*cout << "original info u:" << endl;
			display_array(u, Nv, Nh);*/
			// test: original info

			for (int i = 0; i < Nv; i++)
				PolarEncode_xor(midx[i], u[i], Nh);
			// test:horizontal encoding
			//cout << "u after horizontal encoding:" << endl;
			////for (int j = 0; j < Nh; j++) { cout << midx[Av[0]][j] << " "; }
			//for (int i = 0; i < Nv; i++)
			//{
			//	cout << endl;
			//	for (int j = 0; j < Nh; j++)
			//		cout << midx[i][j];
			//}
			//cout << endl;
			//cout << endl;
			// test:horizontal encoding
			// 
			// column-wise interleaving

			// test: row encoding
			/*cout << "u after horizontal encoding:" << endl;
			display_array(midx, Nv, Nh);*/
			// test: row encoding
			if (flag_interleave_inner)
			{
				int* tmp_to_interleave = new int[Nv];
				for (int i = 0; i < Nh; i++)
				{
					for (int j = 0; j < Nv; j++) { tmp_to_interleave[j] = midx[j][i]; }
					randIntrlvPartial(intrlv_pattern_2D_inner[i], tmp_to_interleave, Av, Nv, Kv);
					for (int j = 0; j < Nv; j++) { midx[j][i] = tmp_to_interleave[j]; }
				}
				delete[] tmp_to_interleave;
			}

			// test: innner interleaving
			/*cout << "u after inner interleaving:" << endl;
			display_array(midx, Nv, Nh);*/
			// test: innner interleaving
//#pragma omp parallel for
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

			// test: vertical encoding
			/*cout << "u after vertical encoding:" << endl;
			display_array(x, Nv, Nh);*/
			// test: vertical encoding
			/*for (int i = 0; i < 2 * Nr; i++)
				for (int j = 0; j < 2 * Nt; j++)
					H_real(i, j) = 0;*/
					/*MatrixXf H_real_test(2 * Nr, 2 * Nt);
					MatrixXf M = H_real_test.transpose();
					MatrixXf N = H_real_test;
					HTHIinv = M * N;*/
					//HTHIinv = M * H_real;
					// test:vertical encoding
					//cout << "u after vertical encoding:" << endl;
					////for (int j = 0; j < Nh; j++) { cout << x[Av[0]][j] << " "; }
					//for (int i = 0; i < Nv; i++)
					//{
					//	cout << endl;
					//	for (int j = 0; j < Nh; j++)
					//		cout << x[i][j];
					//}
					//cout << endl;
					//cout << endl;
					// test:vertical encoding
			if (atten == 2)
			{
				// modulation
				// #pragma omp parallel for
				// Interleaving
				if (flag_interleave)
				{
					if (flag_interleave_all )  // interleave temporally
					{
						for (int i = 0; i < Nv; i++)
						{
							randIntrlv(intrlv_pattern_2D_time[i], x[i], Nh);
						}
					}
					//// test: to be interleaved
					//cout << "x before interleaving:" << endl;
					//// print x line by line, add spaces between each 4 terms
					//for (int i = 0; i < Nv; i++)
					//{
					//	for (int j = 0; j < Nh; j++)
					//	{
					//		if (!(j % 4)) cout << " ";
					//		cout << x[i][j] << " ";
					//	}
					//	cout << endl;
					//}
					
					if (flag_interleave_partial)
					{
						for (int i = 0; i < Nv; i++)
						{
							for (int j = 0; j < intrlv_pattern_num; j++)
							{
								randIntrlv(intrlv_pattern_2D_time_partial[i][j], x[i]+j*ModType, ModType);
							}
						}
					}

					//// test: to be interleaved
					//cout << "x after interleaving:" << endl;
					//// print x line by line, add spaces between each 4 terms
					//for (int i = 0; i < Nv; i++)
					//{
					//	for (int j = 0; j < Nh; j++)
					//	{
					//		if (!(j % 4)) cout << " ";
					//		cout << x[i][j] << " ";
					//	}
					//	cout << endl;
					//}
					if ((!flag_interleave_block) && flag_interleave_spatial)
					{
						for (int j = 0; j < Nh; j++)
						{
							for (int i = 0; i < Nv; i++) { x_tmp[i] = x[i][j]; }
							randIntrlv(intrlv_pattern_2D[j], x_tmp, Nv);
							for (int i = 0; i < Nv; i++) { x[i][j] = x_tmp[i]; }
						}
					}


					if (flag_interleave_block)
					{
						randIntrlvBlock(x, Nh, Nv, ModType, interleave_rows);
					}
				}
				if (MOD_DIM == "Spatial")
				{
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
						/*for (int i = 0; i < mimo_blk_height; i++)
						{
							cout << x_mod_tmp[i] << endl;
						}
						cout << "111111" << endl;*/
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
					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < 2 * Nr; j++) { y_mod_real_temp_omp[i][j] = y_mod_real[j][i]; }
					}
					// for KBest Detection
					Eigen::HouseholderQR<Eigen::MatrixXf> qr_Hreal(H_real);
					Q = qr_Hreal.householderQ();
					R = qr_Hreal.matrixQR().triangularView<Eigen::Upper>();
					//omp_set_num_threads(8);

#pragma omp parallel for private(z)
					for (int j = 0; j < Nh; j++)
					{
						MatrixXf received_signal_omp(2 * Nr, 1);
						int k_kbest_extend = pow(2, ModType / 2) * K_KBEST;
						float** paths_kbest = new float* [k_kbest_extend];
						for (int i = 0; i < k_kbest_extend; i++) paths_kbest[i] = new float[2 * Nt];
						float* ED = new float[K_KBEST];

						if (MIMODET == "MMSE") { MMSEdetection(H_real, Identity, received_signal, mmse_matrix, output_signal, HTH, y_mod_real_temp_omp[j], y_mmse_temp_omp[j], sigma_sq); }
						// ִ��һЩ����

						for (int i = 0; i < 2 * Nr; i++) { received_signal_omp(i, 0) = y_mod_real_temp_omp[j][i]; }
						z = Q.transpose() * received_signal_omp;
						if (MIMODET == "KBEST")
						{
							KBestDetection(R, z, H_real, y_mod_real_temp_omp[j], ED, y_mmse_temp_omp[j], paths_kbest, K_KBEST, ModType, Nt, Nr);
							llr_kbest_bit(Nt, y_mmse_temp_omp[j], oriLLR_kbest[j], sigma_sq, ModType, K_KBEST, paths_kbest, ED);
						}
						for (int i_path = 0; i_path < k_kbest_extend; i_path++) { delete[] paths_kbest[i_path]; }

						delete[] paths_kbest;
						delete[] ED;
					}

					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < 2 * Nt; j++) { y_mmse[j][i] = y_mmse_temp_omp[i][j]; }
					}
					// sym llr to bitwise llr
					// Detection is done for every column, that is, for every Nv bits
					// Interleaving is not currently considered
					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < 2 * Nt; j++) { y_mmse_tmp[j] = y_mmse[j][i]; }
						if (MIMODET == "MMSE") { llr_mmse_bit(Nt, y_mmse_tmp, oriLLR_tmp, sigma_sq, ModType); }
						if (MIMODET == "KBEST") { for (int j = 0; j < len_ori_llr_tmp; j++) { oriLLR_tmp[j] = oriLLR_kbest[i][j]; } }
						if (flag_interleave)
						{
							randDeintrlv(intrlv_pattern, oriLLR_tmp, Nv);
						}
						bool flagIsMax = false;
						for (int j = 0; j < Nv; j++)
						{
							if (!(j % 2))
							{
								if (abs(x_mmse_out[j][i]) > abs(x_mmse_out[j][i + 1])) { flagIsMax = true; }
								else { flagIsMax = false; }
							}
							else
							{
								if (abs(x_mmse_out[j][i]) > abs(x_mmse_out[j][i - 1])) { flagIsMax = true; }
								else { flagIsMax = false; }
							}
							oriLLRh[j][i] = oriLLR_tmp[j];
							upLLR[j][i] = oriLLR_tmp[j];
							// TEST�� det error
							x_mmse_out[j][i] = 0;
							if (oriLLR_tmp[j] < 0) { x_mmse_out[j][i] = 1; }
							if (x_mmse_out[j][i] != x[j][i]) {
								n_det_err++;
								// error at max-reliable bit. Warning: Only fit for 16QAM modulation
								if (flagIsMax) { n_det_err_mrb++; }
							}
							total_bit_cnt++;
						}
					}
					// test : detection TIME
					/*det_end = std::chrono::system_clock::now();
					auto duration_detection = std::chrono::duration_cast<std::chrono::microseconds>(det_end - start);
					file_det << duration_detection.count() << endl;*/
				}
				// Temporal modulation
				else if (MOD_DIM == "Temporal")
				{
					//modulation
					//cout << "START TEST" << endl;
					for (int i = 0; i < Nv; i++) // every row
					{
						modulation(x[i], x_mod[i], ModType, mimo_blk_length);
					}
					/*for (int i = 0; i < mimo_blk_height; i++)
					{
						for (int j = 0; j < mimo_blk_length; j++)
							cout << x_mod[i][j] << endl;
					}
					cout << "111111" << endl;*/
					//AWGN channel (received y_mod)
					// test: modulated symbols
					/*cout << endl;
					cout << "modulated symbols" << endl;
					for (int i = 0; i < Nt; i++) cout << x_mod[i][0] << endl;
					cout << "modulated symbols" << endl;*/
					// test: modulated symbols
					// H = generateChannelMatrix(Nr, Nt);
					if (!non_uniform_channel)H = generateChannelMatrix(Nr, Nt);
					else H = generateChannelMatrixNonuniform(Nr, Nt, channel_scale_factor);
					for (int j = 0; j < Nsym; j++)
					{
						// Extracting a single column
						for (int i = 0; i < Nv; i++) { x_mod_tmp[i] = x_mod[i][j]; }
						AWGN(Nr, Nt, x_mod_tmp, y_mod_tmp, H, sigma_sq);
						for (int i = 0; i < Nr; i++) { y_mod[i][j] = y_mod_tmp[i]; }
					}
					convertHToReal(H, H_real);

					for (int j = 0; j < Nsym; j++)
					{
						for (int i = 0; i < Nr; i++) { y_mod_tmp[i] = y_mod[i][j]; }
						convertYToReal(Nr, y_mod_tmp, y_mod_real_tmp);
						for (int i = 0; i < 2 * Nr; i++) { y_mod_real[i][j] = y_mod_real_tmp[i]; }
					}
					//MMSE detection done in paralell
					for (int i = 0; i < Nsym; i++)
					{
						for (int j = 0; j < 2 * Nr; j++) { y_mod_real_temp_omp[i][j] = y_mod_real[j][i]; }
					}
					Eigen::HouseholderQR<Eigen::MatrixXf> qr_Hreal(H_real);
					Q = qr_Hreal.householderQ();
					R = qr_Hreal.matrixQR().triangularView<Eigen::Upper>();
					// for EP detection
					// copy H_real
					for (int i = 0; i < Nr2; i++)
					{
						for (int j = 0; j < Nt2; j++)
						{
							Hreal_array[i * Nt2 + j] = H_real(i, j);
						}
					}	
					// calculate HTH
					cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, Nt2, Nt2, Nr2, 1.0, Hreal_array, Nt2, Hreal_array, Nt2, 0.0, HTH_array, Nt2);
					//cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, Nt2, 1, Nr2, 1.0, ch_mtx_struct.H_real[cnt], Nt2, data_struct.rx_buffer_real[cnt], 1, 0.0, ch_mtx_struct.hty, 1);
		
					//omp_set_num_threads(8);
					//cout << 1 << endl;
					//omp_set_num_threads(4);
					/*auto start = std::chrono::system_clock::now();*/
					/*MatrixXf H_real_trans(2*Nt, 2*Nr);
					for (int i = 0; i < 2 * Nt; i++)
						for (int j = 0; j < 2 * Nr; j++)
							H_real_trans(i, j) = H_real(j, i);
					HTHIinv = (H_real_trans * H_real + Identity).inverse();*/
# pragma omp parallel for private(z)
					for (int j = 0; j < Nsym; j++)

						//for (int x=0; x<Nsym; x++)
					{
						// int j = x % Nsym;
						MatrixXf received_signal_omp(2 * Nr, 1);
						int k_kbest_extend = pow(2, ModType / 2) * K_KBEST;
						float** paths_kbest = new float* [k_kbest_extend];
						for (int i = 0; i < k_kbest_extend; i++) { paths_kbest[i] = new float[2 * Nt]; }
						float* ED = new float[K_KBEST];
						float* HTY_array_tmp = new float[2 * Nt];
						// Declared arrays for EP
						int flag = 1;
						int temp = 0;
						//float* extr_mean = new float[2 * Nr];
						//float* extr_var = new float[2 * Nr];
						double* vec_alpha_in = new double[2 * Nr];
						double* vec_gamma_in = new double[2 * Nr];
						double* p_symbol_rearrange = new double[Nt2 * pow(2,ModType/2)];
						for (int i = 0; i < Nt2 * pow(2,ModType/2); i++) { p_symbol_rearrange[i] = 1/(pow(2,ModType/2)); }
						for (int i = 0; i < Nr2; i++)
						{
							//extr_mean[i] = 0;
							//extr_var[i] = 0;
							vec_alpha_in[i] = 2;
							vec_gamma_in[i] = 0;
						}

						if (MIMODET == "MMSE")
						{
							MMSEdetection(H_real, Identity, received_signal, mmse_matrix, output_signal, HTH, y_mod_real_temp_omp[j], y_mmse_temp_omp[j], sigma_sq);
							//MMSEdetection_noInverse(H_real, Identity, received_signal, mmse_matrix, output_signal, HTH, HTHIinv, y_mod_real_temp_omp[j], y_mmse_temp_omp[j], sigma_sq);
						}
						if (MIMODET == "KBEST")
						{
							for (int i = 0; i < 2 * Nr; i++) { received_signal_omp(i, 0) = y_mod_real_temp_omp[j][i]; }
							z = Q.transpose() * received_signal_omp;
							KBestDetection(R, z, H_real, y_mod_real_temp_omp[j], ED, y_mmse_temp_omp[j], paths_kbest, K_KBEST, ModType, Nt, Nr);
							llr_kbest_bit(Nt, y_mmse_temp_omp[j], oriLLR_kbest[j], sigma_sq, ModType, K_KBEST, paths_kbest, ED);
						}
						if (MIMODET == "EP")
						{
							// test: save H_real, recv_signal
							/*save_mat_txt(y_mod_real_temp_omp[j], Nr2, "RxSig0128.txt");
							save_mat_txt(H_real, "Hreal0128.txt");*/
							// test: save H_real, recv_signal
							cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, Nt2, 1, Nr2, 1.0, Hreal_array, Nt2, y_mod_real_temp_omp[j], 1, 0.0, HTY_array_tmp, 1);
							//test:HTH_array
							/*for (int i = 0; i < Nt2; i++)
								for (int s = 0; s < Nt2; s++)
								{
									if (s == 0) std::cout << endl;
									std::cout << setw(10) << HTH_array[i * Nt2 + s];
								}
							cout << endl;*/
							// TEST:HTY
							/*for (int i = 0; i < Nr2; i++)
								received_signal_omp(i, 0) = y_mod_real_temp_omp[j][i];
							MatrixXf HTransY = H_real.transpose() * received_signal_omp;
							cout << "Eigen" << endl;
							for (int i = 0; i < Nt2; i++) cout << setw(10)<<HTransY(i, 0);
							cout << "mkl" << endl;
							for (int i = 0; i < Nt2; i++) cout << HTY_array[i];*/
							// TEST:HTY
							EPD_detection_float(flag, Nt, Nr, iter_MIMO, ModType, 1.0, sigma_sq, HTH_array, HTY_array_tmp, temp, ep_extr_mean_omp[j], ep_extr_var_omp[j], vec_alpha_in, vec_gamma_in, p_symbol_rearrange);
						}
						//test EP detection
						/*cout << "TEST:EP DETECTION" << endl;
						for (int i = 0; i < Nr2; i++) std::cout << extr_mean[i]<< endl;
						cout << "TEST:EP DETECTION" << endl;*/
						//test EP detection
						for (int i_path = 0; i_path < k_kbest_extend; i_path++) { delete[] paths_kbest[i_path]; }
						delete[] paths_kbest;
						delete[] ED;
						//delete[] extr_mean;
						//delete[] extr_var;
						delete[] vec_alpha_in;
						delete[] vec_gamma_in;
						delete[] p_symbol_rearrange;
						delete[] HTY_array_tmp;
					}
					// test: detection time
					/*auto detection_over = std::chrono::system_clock::now();
					auto duration_detection = std::chrono::duration_cast<std::chrono::microseconds>(detection_over - start);
					cout << endl;
					cout <<"Detection duration" << duration_detection.count() << "us" << endl;*/
					// test: detection time
					for (int i = 0; i < Nsym; i++)
					{
						for (int j = 0; j < 2 * Nt; j++) { y_mmse[j][i] = y_mmse_temp_omp[i][j]; }
					}
					//To bit domain output
					for (int j = 0; j < Nsym; j++)
					{
						for (int i = 0; i < 2 * Nt; i++) { y_mmse_tmp[i] = y_mmse[i][j]; }
						// With Spatial Modulation, Nt = Nv
						if (MIMODET == "MMSE") { llr_mmse_bit(Nt, y_mmse_tmp, oriLLR_tmp, sigma_sq, ModType); }
						if (MIMODET == "EP") { llr_ep_bit(Nt, ep_extr_mean_omp[j], ep_extr_var_omp[j], oriLLR_tmp, ModType); }
						if (MIMODET == "KBEST") { for (int i = 0; i < len_ori_llr_tmp; i++) { oriLLR_tmp[i] = oriLLR_kbest[j][i]; } }
						int i_bit = 0;
						int bit_col = 0;
						for (int j_inner = 0; j_inner < Nv; j_inner++)
						{
							for (int i_bit_col = 0; i_bit_col < ModType; i_bit_col++)
							{
								i_bit = ModType * j_inner + i_bit_col;
								bit_col = j * ModType + i_bit_col;
								oriLLRh[j_inner][bit_col] = oriLLR_tmp[i_bit];
								upLLR[j_inner][bit_col] = oriLLR_tmp[i_bit];
							}

						}
					}
					// test: LLR transfer time
					/*auto llr_to_bit_over = std::chrono::system_clock::now();
					auto duration_llr = std::chrono::duration_cast<std::chrono::microseconds>(llr_to_bit_over - detection_over);
					cout << endl;
					cout<< "LLR transfer duration" << duration_llr.count() << "us" << endl;*/
					// test: LLR transfer time
					if (flag_interleave)
					{
						if ((!flag_interleave_block) && flag_interleave_spatial)
						{
							for (int j = 0; j < Nh; j++)
							{
								for (int i = 0; i < Nv; i++) { oriLLR_tmp[i] = oriLLRh[i][j]; }
								randDeintrlv(intrlv_pattern_2D[j], oriLLR_tmp, Nv);
								for (int i = 0; i < Nv; i++)
								{
									oriLLRh[i][j] = oriLLR_tmp[i];
									upLLR[i][j] = oriLLR_tmp[i];
								}

							}
						}

						if (flag_interleave_block)
						{
							randDeIntrlvBlock(oriLLRh, Nh, Nv, ModType, interleave_rows);
							randDeIntrlvBlock(upLLR, Nh, Nv, ModType, interleave_rows);
						}

						if (flag_interleave_all) // Deinterleave in time domain
						{
							for (int j = 0; j < Nv; j++)
							{
								randDeintrlv(intrlv_pattern_2D_time[j], oriLLRh[j], Nh);
								randDeintrlv(intrlv_pattern_2D_time[j], upLLR[j], Nh);
							}
						}

						if (flag_interleave_partial)
						{
							for (int j = 0; j < Nv; j++)
							{
								for (int k = 0; k < intrlv_pattern_num; k++)
								{
									randDeintrlv(intrlv_pattern_2D_time_partial[j][k], oriLLRh[j] + k * ModType, ModType);
									randDeintrlv(intrlv_pattern_2D_time_partial[j][k], upLLR[j] + k * ModType, ModType);
								}
							}
						}

					}
					
					for (int i = 0; i < Nv; i++)
					{
						for (int j = 0; j < Nh; j++)
						{
							total_bit_cnt++;
							int bit_dec = (oriLLRh[i][j] < 0);
							if (bit_dec != x[i][j]) n_det_err++;
						}
					}

				}
				// test: detection errors
				
				// test: detection errors
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
			for (int iter = 1; iter <= softiter; iter++)
			{
				//vertical decoding
		   // test: decoding time
			//auto decoding_start = std::chrono::system_clock::now();

			// test: decoding time
#pragma omp parallel for
				for (int decodingi = 0; decodingi < Nh; decodingi++) // 291.70kb
				{

					//for SCL Decoding
					int Counth = 0, Countinfoh = 0;
					float Qsumunh = 0;
					int** sum = new int* [L];
					//for (int i = 0; i < L; i++)  sum[i] = new  int[Nmax];
					int** u1 = new int* [L];
					for (int i = 0; i < L; i++)  u1[i] = new _MM_ALIGN16 int[Nmax];
					for (int i = 0; i < L; i++)  sum[i] = new _MM_ALIGN16 int[Nmax];
					//for (int i = 0; i < L; i++)  u1[i] = new  int[Nmax];

					float* LLR_in = new float[L];
					float* W = new float[2 * L];
					int* SCL_index = new int[L];
					int* better = new int[L];
					int* worse = new int[L];
					int* indexpm = new int[L];
					float* PM = new float[L];
					float** LLR = new float* [L];
					for (int i = 0; i < L; i++)  LLR[i] = new _MM_ALIGN16 float[2 * Nmax + 2];
					// for (int i = 0; i < L; i++)  LLR[i] = new float[2 * Nmax + 2];
					////////////////////////////////
					for (int j = 0; j < Nv; j++)
						for (int k = 0; k < L; k++)
							LLR[k][j] = upLLR[j][decodingi];
					//cout << LLR[0][j] << " ";

				//cout << endl;
					int index = 0;
					/*			if (L == 1)
								{
									//T_start = clock();
									Count = 0;
									decode(*LLR, *sum, (int)(log(Nv) / log(2)), A_Acv, flag, Nv, Kv);
									//PolarEncode_xor(*u1, *sum, N);
									//T_finish = clock();
								}
								else
								{*/
					for (int k = 0; k < L; k++) { PM[k] = 0; indexpm[k] = k; }
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
					if (adaptive_beta)
					{
						if (softmethod == "Pyndiah") Pyndiahv_adaptivebeta(LLR, sum, oriLLRh, indexpm, iter * 2 - 2, decodingi, upLLR);
						if (softmethod == "MITSO") MITSOv(LLR, PM, sum, oriLLRh, indexpm, iter * 2 - 2, decodingi, upLLR, Qsumunh);
					}
					// iter * 2 - 2
					else
					{
						//cout << iter << endl;
						if (softmethod == "Pyndiah") Pyndiahv(LLR, sum, oriLLRh, indexpm, iter, decodingi, upLLR);
						if (softmethod == "MITSO") MITSOv(LLR, PM, sum, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunh);
					}
					// test: oriLLR and upLLR
					/*cout << endl;
					for (int i = 0; i < Nv; i++)
						for (int j = 0; j < Nh; j++)
							cout << oriLLRh[i][j] - upLLR[i][j] << "   ";
					cout << endl;*/
					/*cout << oriLLRh[0][0] << endl;
					cout << upLLR[0][0] << endl;*/
					// test: oriLLR and upLLR
					// half iteration Enabled
					if (half_iter && iter == softiter)
					{
						if (iter == softiter) PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), Nh);
						for (int j = 0; j < Kv; j++)
							umid[Av[j]][decodingi] = u1[index][Av[j]];
					}
					// half iteration
					////////delete SCL variables
					for (int i = 0; i < L; i++)
						delete[]LLR[i];
					delete[] LLR;
					for (int i = 0; i < L; i++)
						delete[]sum[i];
					for (int i = 0; i < L; i++)
						delete[]u1[i];
					delete[] u1;
					delete[] sum;
					delete[] PM;
					delete[] LLR_in;
					delete[] W;
					delete[] SCL_index;
					delete[] better;
					delete[] worse;
					delete[] indexpm;

				}
				//cout << "Horizon Start:\n";
				//horizontal decoding
				// #
				if (weirditer) { continue; }

				// test: LLR after horizontal decoding
				/*cout << "After horizontal decoding" << endl;
				for (int i = 0; i < Nv; i++)
				{
					for (int j = 0; j < Nh; j++)
					{
						cout << (oriLLRh[i][j] < 0) << " ";
					}
					cout << endl;
				}*/
				// test: LLR after horizontal decoding
				// inner deinterleaving
				if (flag_interleave_inner)
				{
					float* tmp_llr_to_interleave = new float[Nv];
					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < Nv; j++) { tmp_llr_to_interleave[j] = oriLLRh[j][i]; }
						randDeIntrlvPartial(intrlv_pattern_2D_inner[i], tmp_llr_to_interleave, Av, Nv, Kv);
						for (int j = 0; j < Nv; j++) { oriLLRh[j][i] = tmp_llr_to_interleave[j]; }
					}
					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < Nv; j++) { tmp_llr_to_interleave[j] = upLLR[j][i]; }
						randDeIntrlvPartial(intrlv_pattern_2D_inner[i], tmp_llr_to_interleave, Av, Nv, Kv);
						for (int j = 0; j < Nv; j++) { upLLR[j][i] = tmp_llr_to_interleave[j]; }
					}
					delete[] tmp_llr_to_interleave;
				}


				/*auto decoding_over = std::chrono::system_clock::now();
				auto duration_decoding = std::chrono::duration_cast<std::chrono::microseconds>(decoding_over - decoding_start);
				cout << endl;
				cout << "Decoding duration" << duration_decoding.count() << "us" << endl;*/

#pragma omp parallel for
				for (int decodingi = 0; decodingi < Nv; decodingi++)
				{
					// 293.70kb
					// cout << "KKKK" << endl;
					if (half_iter && (iter == softiter)) { break; }
					string polarde_now = polarde_array[iter - 1];
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

					/*			if (L == 1)
								{
									//T_start = clock();
									Count = 0;
									decode(*LLR, *sum, (int)(log(Nh) / log(2)), A_Ach, flag, Nh, Kh);
									if (iter == softiter) PolarEncode_xor(*u1, *sum, Nh);
									//T_finish = clock();
								}
								else
								{*/
								//for (int k = 0; k < L; k++)  PM[k] = 0, indexpm[k] = k; 
								////T_start = clock();
								//decode_list(LLR, sum, (int)(log(Nh) / log(2)), A_Ach, 0, Nh, Kh, PM, L, (int)(log(L) / log(2)), LLR_in, W, SCL_index, better, worse, &Countv, &Countinfov, &Qsumunv);
								////T_finish = clock();
								////Decoding_Duration += (T_finish - T_start);
								//for (int i = 0; i < L - 1; i++)
								//	for (int j = i + 1; j < L; j++)
								//		if (PM[indexpm[i]] < PM[indexpm[j]])
								//		{
								//			int ttt = indexpm[i];
								//			indexpm[i] = indexpm[j];
								//			indexpm[j] = ttt;
								//		}


							/*
							* Encoding: Info -> Horizontal Encoding -> A -> set to 0 -> B -> Horizontal Encoding -> Vertical Encoding
							* Decoding: Received Info -> Horizontal Decoding -> Vertical Decoding -> Horizontal Decoding
							*/
					if (polarde_now == "SLD")
					{
						LSD_decode_outall(u_decod, u1_LSD, LLR[0], GGh, A_Ach, Nh, L_LSD, n_paths);
						// if (iter == softiter) PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), Nh);
						// for soft information exchange, inverse coding:
						// THIS LINE is WRONG! BUT IT SEEMS BETTER THAN THE RIGHT VERSION
						for (int i = 0; i < L_LSD; i++) { PolarEncode_xor(sum_LSD[i], u1_LSD[i], Nh); }
						// for (int i = 0; i < L_LSD; i++) { PolarEncode_xor(sum_LSD[i], u1_LSD[i], Nh); }
						// if (iter == softiter) PolarEncode_xor(u1[0], sum_LSD[0], Nh);
						index = indexpm[0];
						// test: after horizontal decoding
						/*if (decodingi == Av[0])
						{
							cout << "Direct output of horizontal decoding:" << endl;
							for (int i = 0; i < Nh; i++) { cout << u1_LSD[0][i] << " "; }
							cout << endl;
							cout << "horizontal decoder out, encode once:" << endl;
							for (int i = 0; i < Nh; i++) { cout << sum_LSD[0][i] << " "; }
							cout << endl;
						}*/
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
						// iter*2-1
						if (iter < softiter)
						{
							if (adaptive_beta)
							{
								if (softmethod == "Pyndiah") Pyndiahh_LSD_adaptivebeta(LLR, sum_LSD, oriLLRh, indexpm, iter, decodingi, upLLR);
								if (softmethod == "MITSO") MITSOh(LLR, PM, sum_LSD, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunv);
							}
							else
							{
								// cout << "IN LSD!!!" << endl;
								if (softmethod == "Pyndiah")
									Pyndiahh_LSD_getsum(LLR, sum_LSD, oriLLRh, indexpm, iter, decodingi, upLLR, sum_sdis_array[decodingi], GGh, Ach);
								//if (softmethod == "Pyndiah") Pyndiahh_LSD(LLR, sum_LSD, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR);
								if (softmethod == "MITSO") MITSOh(LLR, PM, sum_LSD, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunv);
							}
							//	getchar();
						}
						else
							for (int j = 0; j < Kh; j++)
								umid[decodingi][Ah[j]] = u1_LSD[index][Ah[j]];
					}
					if (polarde_now == "SCL")
					{
						for (int k = 0; k < L; k++)  PM[k] = 0, indexpm[k] = k;
						//T_start = clock();
						decode_list(LLR, sum, (int)(log(Nh) / log(2)), A_Ach, 0, Nh, Kh, PM, L, (int)(log(L) / log(2)), LLR_in, W, SCL_index, better, worse, &Countv, &Countinfov, &Qsumunv);
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
						if (iter == softiter) PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), Nh);
						index = indexpm[0];
						//}
						//calc distance s()
						if (iter < softiter)
						{
							if (adaptive_beta)
							{
								if (softmethod == "Pyndiah") Pyndiahh_adaptivebeta(LLR, sum, oriLLRh, indexpm, iter, decodingi, upLLR);
								if (softmethod == "MITSO") MITSOh(LLR, PM, sum, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunv);
								//	getchar();
							}
							else
							{
								//if (softmethod == "Pyndiah") Pyndiahh(LLR, sum, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR);
								if (softmethod == "Pyndiah") Pyndiahh_getsum(LLR, sum, oriLLRh, indexpm, iter, decodingi, upLLR, sum_sdis_array[decodingi]);
								if (softmethod == "MITSO") MITSOh(LLR, PM, sum, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunv);
								//	getchar();
							}
						}
						else
							for (int j = 0; j < Kh; j++)
								umid[decodingi][Ah[j]] = u1[index][Ah[j]];
					}

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


				if (flag_interleave_inner)
				{
					float* tmp_llr_to_interleave = new float[Nv];
					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < Nv; j++) { tmp_llr_to_interleave[j] = oriLLRh[j][i]; }
						randIntrlvPartial(intrlv_pattern_2D_inner[i], tmp_llr_to_interleave, Av, Nv, Kv);
						for (int j = 0; j < Nv; j++) { oriLLRh[j][i] = tmp_llr_to_interleave[j]; }
					}
					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < Nv; j++) { tmp_llr_to_interleave[j] = upLLR[j][i]; }
						randIntrlvPartial(intrlv_pattern_2D_inner[i], tmp_llr_to_interleave, Av, Nv, Kv);
						for (int j = 0; j < Nv; j++) { upLLR[j][i] = tmp_llr_to_interleave[j]; }
					}
					delete[] tmp_llr_to_interleave;
				}
			}
			int calcsum[Nv], calcu1[Nv];
			int calcu1h[Nh];
			// calcsum is from the direct output of the decoder
			// if output is the result of row decoding (meaning there is another column decoding to be done)
			if (half_iter)
			{
				for (int i = 0; i < Kv; i++)
				{
					PolarEncode_xor(calcu1h, umid[Av[i]], Nh);
					for (int j = 0; j < Kh; j++)
						uout[Av[i]][Ah[j]] = calcu1h[Ah[j]];
				}
			}
			// if output is the result of column decoding (meaning that there is another row decoding to be done)
			else
			{
				for (int i = 0; i < Kh; i++)
				{
					for (int j = 0; j < Nv; j++)
						calcsum[j] = umid[j][Ah[i]];
					PolarEncode_xor(calcu1, calcsum, Nv);
					for (int j = 0; j < Kv; j++)
						uout[Av[j]][Ah[i]] = calcu1[Av[j]];
				}
			}
			// test: final output
			/*cout << "final output" << endl;
			for (int k = 0; k < Nh; k++) { cout << uout[Av[0]][k] << " "; }
			cout << endl;*/
			// test: final output
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
			if (Temp_Error < Num_Error) {
				Num_Frame_Error++; if (showllr_irr) {
					std::cout << endl << "ERR" << endl;
					for (int i = 0; i < Nh; i++) { std::cout << upLLR[3][i] << "  "; }
					std::cout << endl;
				}
			}
			if (frame % 10 == 0)
			{
				/*printf("\r<%d> SNR= %.2f FER= %d / %d = %f BER= %d / %d = %lf E-Det = %f, E-MRB = %f"  ,
					frame, nSNR, Num_Frame_Error, frame, (double)Num_Frame_Error / frame, Num_Error, (N_Info * frame), (double)Num_Error / N_Info / frame,
					(double)n_det_err/total_bit_cnt, (double)n_det_err_mrb/(total_bit_cnt/(ModType/2)));*/
				int L_sdis_sum = L_LSD;
				sdis_sum = 0;
				for (int i = 0; i < Nv; i++) sdis_sum += sum_sdis_array[i];
				if (polarde_array[0] == "SCL") L_sdis_sum = L;
				// cout << "sdis_sum_avg" << sdis_sum / frame / L_sdis_sum << endl;
				printf("\r<%d> SNR= %.2f FER= %d / %d = %f BER= %d / %d = %lf ,sdis_avg = %f, E-MRB = %f",
					frame, nSNR, Num_Frame_Error, frame, (double)Num_Frame_Error / frame, Num_Error, (N_Info * frame), (double)Num_Error / N_Info / frame, (double)sdis_sum / frame / L_sdis_sum, (double)n_det_err / total_bit_cnt);
				fflush(stdout);
			}

			// Lightweight progress log: write one line every 10 newly-errored frames per SNR point.
			// Uses a block-scope static counter; resets when a new SNR point starts (Num_Frame_Error==0).
			{
				static int s_last_logged_err = -1;
				if (Num_Frame_Error == 0) s_last_logged_err = -1;
				if (file_progress.is_open() && Num_Frame_Error > 0
					&& (Num_Frame_Error % 10 == 0)
					&& Num_Frame_Error != s_last_logged_err) {
					file_progress << nSNR << "\t" << frame << "\t" << Num_Frame_Error
						<< "\t" << ((double)Num_Frame_Error / frame)
						<< "\t" << Num_Error << "\t" << (N_Info * frame)
						<< "\t" << ((double)Num_Error / (double)N_Info / (double)frame) << "\n";
					file_progress.flush();
					s_last_logged_err = Num_Frame_Error;
				}
			}

			
			// output interleaving pattern
			//if (Num_Frame_Error % record_frames == 0 && Num_Frame_Error > 0 && Num_Frame_Error!=n_frame_error_pre)
			//{
			//	n_frame_error_pre = Num_Frame_Error;
			//	cout << endl;
			//	int n_file = Num_Frame_Error / record_frames;
			//	frame_delta = frame - frame_tmp;
			//	cout << "Delta Frames: " << frame_delta << endl;
			//	frame_tmp = frame;
			//	string ofile_intrlv_name = string("intrlv_pattern_") + to_string(n_file) + string("_frame_") + to_string(frame_delta) + string(".txt");
			//	string ofile_spatial_intrlv_name = string("intrlv_pattern_spatial") + to_string(n_file) + string(".txt");
			//	ofstream ofile_intrlv(ofile_intrlv_name);
			//	ofstream ofile_spatial_intrlv(ofile_spatial_intrlv_name);
			//	for (int i = 0; i < Nv; i++)
			//	{
			//		for (int j = 0; j < intrlv_pattern_num; j++)
			//		{

			//			for (int k = 0; k < ModType; k++)
			//				ofile_intrlv << intrlv_pattern_2D_time_partial[i][j][k] << " ";
			//			//cout << endl;
			//		}
			//	}
			//	ofile_intrlv.close();

			//	for (int i = 0; i < Nh; i++)
			//	{
			//		for (int j = 0; j < Nv; j++)
			//			ofile_spatial_intrlv << intrlv_pattern_2D[i][j] << " ";
			//	}
			//	ofile_spatial_intrlv.close();

			//	// re-generate interleaving pattern
			//	generateIntrlvPattern(intrlv_pattern, Nv);
			//	for (int i = 0; i < Nh; i++) { generateIntrlvPattern(intrlv_pattern_2D[i], Nv); }
			//	for (int i = 0; i < Nv; i++) { generateIntrlvPattern(intrlv_pattern_2D_time[i], Nh); }
			//	if (flag_interleave_partial)
			//	{
			//		for (int i = 0; i < Nv; i++)
			//		{
			//			for (int j = 0; j < intrlv_pattern_num; j++)
			//			{
			//				generateIntrlvPattern(intrlv_pattern_2D_time_partial[i][j], ModType);
			//				/*for (int k = 0; k < ModType; k++)
			//					cout << intrlv_pattern_2D_time_partial[i][j][k];
			//				cout << endl;*/
			//			}
			//		}
			//	}
			//	//Num_Frame_Error = 0;
			//}


			//printf("%d %d", Num_Frame_Error, Num_Error);
			//getchar();
			//Decoding_Duration += (T_finish - T_start);
			/*auto dec_end = std::chrono::system_clock::now();
			auto duration_decoding = std::chrono::duration_cast<std::chrono::microseconds>(dec_end - det_end);
			auto duration_total = std::chrono::duration_cast<std::chrono::microseconds>(dec_end - start);
			file_dec << duration_decoding.count() << endl;
			file_total << duration_total.count() << endl;*/
			if (frame % 16 == 0)
			{
				file_err_num << Num_Error << endl;
			}

		}
		//Decoding_Duration = Decoding_Duration * 100 / CLOCKS_PER_SEC / frame;
		//Decoding_Duration = Decoding_Duration / CLOCKS_PER_SEC / Max_frame;

		//printf("%f  %f  %f  %f\n", nSNR, Decoding_Duration, (double)Num_Error / N_Info / Max_frame, (double)Num_Frame_Error / Max_frame);
		printf("%f %f %f\n", nSNR, (double)Num_Frame_Error / frame, (double)Num_Error / N_Info / frame);
		BER_array[SNR_count] = (double)Num_Error / N_Info / frame;
		FER_array[SNR_count] = (double)Num_Frame_Error / frame;
		// save to .txt file
		string polarde_serial = "";
		for (int i = 0; i < softiter; i++)
			polarde_serial = polarde_serial + polarde_array[i] + string("_");
		string file_name = string("ResFinal\\Regular_") + polarde_serial + MIMODET + string("T") + to_string(Nt) + string("R") + to_string(Nr)
			+ string("Q") + to_string(int(pow(2, ModType))) + string("L") + to_string(Nh) + string("K") + to_string(Kh)
			+ string("T") + to_string(softiter) + string("H") + to_string(half_iter) + string("UB") + to_string(usebeta) + string("IS") + to_string(flag_interleave_spatial)
			+ string("DG") +to_string(diff_gain) + string(".txt");
		save_array_txt(nSNR,  FER_array[SNR_count], file_name);
		int L_sdis_sum = L_LSD;
		sdis_sum = 0;
		for (int i = 0; i < Nv; i++) sdis_sum += sum_sdis_array[i];
		if (polarde_array[0] == "SCL") L_sdis_sum = L;
		// cout << "sdis_sum_avg" << sdis_sum / frame / L_sdis_sum << endl;
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
	for (int i = 0; i < Nh; i++)
		delete[] oriLLR_kbest[i];
	delete[] oriLLR_kbest;
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
	for (int i = 0; i < mimo_blk_height; i++) {
		delete[] x_mod[i];
	}
	delete[] x_mod;
	for (int i = 0; i < Nr; i++) {
		delete[] y_mod[i];
	}
	delete[] y_mod;
	for (int i = 0; i < 2 * mimo_blk_height; i++) {
		delete[] y_mmse[i];
	}
	delete[] y_mmse;
	for (int i = 0; i < 2 * Nr; i++) {
		delete[] y_mod_real[i];
	}
	for (int i = 0; i < mimo_blk_length; i++)
	{
		delete[] y_mod_real_temp_omp[i];
	}
	delete[] y_mod_real_temp_omp;
	for (int i = 0; i < mimo_blk_length; i++)
	{
		delete[] y_mmse_temp_omp[i];
	}
	delete[] y_mmse_temp_omp;
	for (int i = 0; i < mimo_blk_length; i++)
	{
		delete[] ep_extr_mean_omp[i];
	}
	delete[] ep_extr_mean_omp;
	for (int i = 0; i < mimo_blk_length; i++)
	{
		delete[] ep_extr_var_omp[i];
	}
	delete[] ep_extr_var_omp;
	delete[] y_mod_real;
	delete[] x_mod_tmp;
	delete[] y_mod_tmp;
	delete[] y_mod_real_tmp;
	delete[] x_tmp;
	delete[] oriLLR_tmp;
	delete[] sum_sdis_array;
	delete[] Hreal_array;
	delete[] HTH_array;
	delete[] HTY_array;
	delete[] intrlv_pattern_2D_inner;
}




void system2D_LSD_epdet_with_CRC(double* BER_array, double* FER_array)
{
	std::ofstream file_crc_check_replacement_codeword("crc_replacement_codeword_log.txt", std::ios::app);
	std::ofstream file_crc_check_decoding("crc_check_decoding.txt", std::ios::app);
	std::ofstream file_total("total_timing.txt", std::ios::app);
	std::ofstream file_det("det_timing.txt", std::ios::app);
	std::ofstream file_dec("dec_timing.txt", std::ios::app);
	std::ofstream file_err_num("total_frames.txt", std::ios::app);
	auto det_end = std::chrono::system_clock::now();
	int Nt2 = 2 * Nt;
	int Nr2 = 2 * Nr;
	
	/*if (!file.is_open()) {
		std::cerr << "Error opening file: " << file_name << std::endl;
		return;
	}*/

	//file << setw(20) << SNR << setw(20) << BER << setw(20) << FER << "\n";

	//freopen("result.txt", "w", stdout);
	printsetting();
	cout << "\nSNR	 " << "FER      BER" << endl;
	int N_Info = Kv * Kh - CRCL, Nmax;
	//Interleaving pattern
	/*std::vector<int> intrlv_pattern(Nv);
	std::vector<int>* intrlv_pattern_2D = new vector<int>[Nh];
	for (int i = 0; i < Nh; i++) intrlv_pattern_2D[i].resize(Nv);*/
	int intrlv_pattern_num = Nh / ModType;
	std::vector<int> intrlv_pattern(Nv);
	std::vector<int>* intrlv_pattern_2D = new vector<int>[Nh];
	std::vector<int>* intrlv_pattern_2D_time = new vector<int>[Nv];
	std::vector<int>* intrlv_pattern_2D_inner = new vector<int>[Nh];
	std::vector<int>** intrlv_pattern_2D_time_partial = new vector<int> *[Nv];
	for (int i = 0; i < Nv; i++) intrlv_pattern_2D_time_partial[i] = new vector<int>[intrlv_pattern_num];
	for (int i = 0; i < Nv; i++)
	{
		for (int j = 0; j < intrlv_pattern_num; j++)
		{
			intrlv_pattern_2D_time_partial[i][j].resize(ModType);
		}
	}
	for (int i = 0; i < Nh; i++) intrlv_pattern_2D[i].resize(Nv);
	for (int i = 0; i < Nh; i++) intrlv_pattern_2D_inner[i].resize(Kv);
	for (int i = 0; i < Nv; i++) intrlv_pattern_2D_time[i].resize(Nh);

	float CodeRate = (float)N_Info / (Nh * Nv);
	Nmax = max(Nh, Nv);
	// num of symbols per vertical codeword
	int Nsym = 0;
	if (MOD_DIM == "Spatial") { Nsym = Nv / ModType; }
	else if (MOD_DIM == "Temporal") { Nsym = Nh / ModType; }
	int symbol_length = ModType;

	// Size of a mimo block is given by height(spatial) * length(temporal)
	int mimo_blk_height = 0;
	int mimo_blk_length = 0;
	if (MOD_DIM == "Spatial")
	{
		mimo_blk_height = Nsym;
		mimo_blk_length = Nh;
	}
	else if (MOD_DIM == "Temporal")
	{
		mimo_blk_height = Nv;
		mimo_blk_length = Nsym;
	}
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
	int** x_encoded_copy = new int* [Nv];
	for (int i = 0; i < Nv; i++) { x_encoded_copy[i] = new int[Nh]; }

	// Arrays for MIMO
	int** x_mmse_out = new int* [Nv]; // delete!
	for (int i = 0; i < Nv; i++) { x_mmse_out[i] = new int[Nh]; }
	int** x_parity_check = new int* [Nv];
	for (int i = 0; i < Nv; i++) { x_parity_check[i] = new int[Nh]; }

	//std::complex<float>** x_mod = new std::complex<float>* [Nsym]; // modulated symbols (spatial modulation) DELETE!
	//for (int i = 0; i < Nsym; i++) { x_mod[i] = new std::complex<float>[Nh]; }
	//std::complex<float>** x_mod_temporal = new std::complex<float>*[Nv]; // modulated symbols (temporal modulation) DELETE!
	//for (int i = 0; i < Nv; i++) x_mod_temporal[i] = new std::complex<float>[Nsym];
	//std::complex<float>** y_mod = new std::complex<float>* [Nr]; // Output symbols after AWGN channels // delete!
	//for (int i = 0; i < Nr; i++) { y_mod[i] = new std::complex<float>[Nh]; }
	//float** y_mod_real = new float* [2 * Nr];
	//for (int i = 0; i < 2*Nr; i++) y_mod_real[i] = new float[Nh];
	//float** y = new float* [Nv];// Output Signal After AWGN Channel
	//for (int i = 0; i < Nv; i++) { y[i] = new float[Nh]; }
	//float** y_mmse = new float* [2*Nsym];
	//for (int i = 0; i < 2*Nsym; i++) { y_mmse[i] = new float[Nh]; }  // MMSE detector output // delete!
	//float* y_mmse_tmp = new float[2*Nsym];
	//std::complex<float>* x_mod_tmp = new std::complex<float>[Nsym]; // modulated symbols of a single column // delete!
	//std::complex<float>* y_mod_tmp = new std::complex<float>[Nr]; // received symbols of a single column // delete!
	//float* y_mod_real_tmp = new float[2*Nr];
	//// float* y_mod_tmp = new float[Nr];  // delete!


	std::complex<float>** x_mod = new std::complex<float>*[mimo_blk_height]; // modulated symbols (spatial modulation) DELETE!
	for (int i = 0; i < mimo_blk_height; i++) { x_mod[i] = new std::complex<float>[mimo_blk_length]; }
	std::complex<float>** y_mod = new std::complex<float>*[Nr]; // Output symbols after AWGN channels // delete!
	for (int i = 0; i < Nr; i++) { y_mod[i] = new std::complex<float>[mimo_blk_length]; }
	float** y_mod_real = new float* [2 * Nr];
	for (int i = 0; i < 2 * Nr; i++) y_mod_real[i] = new float[mimo_blk_length];
	float** y = new float* [Nv];// Output Signal After AWGN Channel
	for (int i = 0; i < Nv; i++) { y[i] = new float[Nh]; }
	float** y_mmse = new float* [2 * mimo_blk_height];
	for (int i = 0; i < 2 * mimo_blk_height; i++) { y_mmse[i] = new float[mimo_blk_length]; }  // MMSE detector output // delete!
	float* y_mmse_tmp = new float[2 * mimo_blk_height];
	std::complex<float>* x_mod_tmp = new std::complex<float>[mimo_blk_height]; // modulated symbols of a single column // delete!
	std::complex<float>* y_mod_tmp = new std::complex<float>[Nr]; // received symbols of a single column // delete!
	float* y_mod_real_tmp = new float[2 * Nr];
	// float* y_mod_tmp = new float[Nr];  // delete!
	float** y_mod_real_temp_omp = new float* [mimo_blk_length]; // new delete!
	for (int i = 0; i < mimo_blk_length; i++) { y_mod_real_temp_omp[i] = new float[2 * Nr]; }
	float** y_mmse_temp_omp = new float* [mimo_blk_length]; // new delete!
	for (int i = 0; i < mimo_blk_length; i++) { y_mmse_temp_omp[i] = new float[2 * mimo_blk_height]; }
	// ep
	float** ep_extr_mean_omp = new float* [mimo_blk_length];
	for (int i = 0; i < mimo_blk_length; i++) { ep_extr_mean_omp[i] = new float[2 * mimo_blk_height]; }
	float** ep_extr_var_omp = new float* [mimo_blk_length];
	for (int i = 0; i < mimo_blk_length; i++) { ep_extr_var_omp[i] = new float[2 * mimo_blk_height]; }
	int* x_tmp = new int[Nv]; // delete!


	// paths selected by k-best detection
	// array for storing results of bit-flip
	int n_samples = pow(2, p);
	int*** tried_codewords_h = new int** [Nv];  // ��һά��Nv��ָ��

	for (int i = 0; i < Nv; ++i) {
		tried_codewords_h[i] = new int* [n_samples];  // �ڶ�ά��ÿ����һάԪ��ָ��n_samples��ָ��

		for (int j = 0; j < n_samples; ++j) {
			// ����ά��ÿ���ڶ�άԪ��ָ��Nh������
			tried_codewords_h[i][j] = new int[Nh]();  // ĩβ��()����ֵ��ʼ����0��
		}
	}

	int*** tried_codewords_v = new int** [Nh];  // ��һά��Nv��ָ��

	for (int i = 0; i < Nh; ++i) {
		tried_codewords_v[i] = new int* [n_samples];  // �ڶ�ά��ÿ����һάԪ��ָ��n_samples��ָ��

		for (int j = 0; j < n_samples; ++j) {
			// ����ά��ÿ���ڶ�άԪ��ָ��Nh������
			tried_codewords_v[i][j] = new int[Nv]();  // ĩβ��()����ֵ��ʼ����0��
		}
	}

	int* n_tried_patterns_h_all = new int[softiter];
	for (int i = 0; i < softiter; i++) n_tried_patterns_h_all[i] = 0;
	int* n_tried_patterns_v_all = new int[softiter];
	for (int i = 0; i < softiter; i++) n_tried_patterns_v_all[i] = 0;

	// CRC bits
	unsigned char* u_crc = new unsigned char[Kv * Kh];

	float** upLLR = new float* [Nv];
	for (int i = 0; i < Nv; i++) { upLLR[i] = new _MM_ALIGN16 float[Nh]; }
	//float** oriLLRv = new float* [Nh];
	//for (int i = 0; i < Nh; i++) { oriLLRv[i] = new _MM_ALIGN16 float[2 * Nv + 2]; }
	float** oriLLRh = new float* [Nv];
	for (int i = 0; i < Nv; i++) { oriLLRh[i] = new _MM_ALIGN16 float[Nh]; }
	int len_ori_llr_tmp = 0;
	if (MOD_DIM == "Spatial") { len_ori_llr_tmp = Nv; }
	else if (MOD_DIM == "Temporal") { len_ori_llr_tmp = ModType * Nt; }
	float* oriLLR_tmp = new _MM_ALIGN16 float[len_ori_llr_tmp];  // delete!
	//float** oriLLR_kbest = new float* [Nh];
	float** oriLLR_kbest = new float* [mimo_blk_length];
	for (int i = 0; i < Nh; i++) { oriLLR_kbest[i] = new float[len_ori_llr_tmp]; }
	double* sum_sdis_array = new double[Nv];
	for (int i = 0; i < Nv; i++) sum_sdis_array[i] = 0;

	//*****************************************************************
	// EP detector
	float* Hreal_array = new float[(2 * Nr) * (2 * Nt)]; // delete
	float* HTH_array = new float[(2 * Nt) * (2 * Nt)]; // delete
	float* HTY_array = new float[2 * Nt]; // delete 





	int** GGh = new int* [Nh];
	int** GGv = new int* [Nv];


	MatrixXf Gh = Fn(Nh);  // Generate Matrix of horizontal coding
	MatrixXf Gv = Fn(Nv);  // Generate Matrix of vertical coding
	MatrixXcf H(Nr, Nt); // Complex channel Matrix H: Nr x Nt
	//Matrix<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic> H = MatrixXcf::Random(Nr, Nt); // Real domain channel matrix
	//// //MatrixX Definitions
	MatrixXf H_real(2 * Nr, 2 * Nt); // Real domain channel matrix
	MatrixXf received_signal(2 * Nr, 1);
	MatrixXf HTH(2 * Nt, 2 * Nt);
	MatrixXf Identity = MatrixXf::Identity(2 * Nt, 2 * Nt);
	MatrixXf output_signal(2 * Nt, 1);
	MatrixXf mmse_matrix(2 * Nt, 2 * Nr);
	MatrixXf HTHIinv(2 * Nt, 2 * Nt);

	// for KBest
	MatrixXf R(2 * Nt, 2 * Nt);
	MatrixXf Q(2 * Nr, 2 * Nt);
	MatrixXf z(2 * Nt, 1);

	// for CRC
	int** whole_codeword_list = new int* [L_WHOLE_CODEWORD];
	for (int i = 0; i < L_WHOLE_CODEWORD; i++) { whole_codeword_list[i] = new int[Nh * Nv]; }

	int** whole_codeword_list_v = new int* [L_WHOLE_CODEWORD];
	for (int i = 0; i < L_WHOLE_CODEWORD; i++) { whole_codeword_list_v[i] = new int[Nh * Nv]; }

	int* n_err_codeword_v = new int[Nh];  // err of codeword at each column

	int** n_err_codeword_v_diffpos = new int*[Nh];
	for (int i = 0; i < Nh; i++) n_err_codeword_v_diffpos[i] = new int[L + 1];
	for (int i = 0; i < Nh; i++)
	{
		for (int j = 0; j < L + 1; j++) { n_err_codeword_v_diffpos[i][j] = 0; }
	}

	int* correct_idx_oneframe = new int[Nh];
	for (int i = 0; i < Nh; i++) { correct_idx_oneframe[i] = 0; }

	int* max_L_cnt = new int[L+1];
	for (int i = 0; i < L+1; i++) { max_L_cnt[i] = 0; }

	// for CRC
	
	int L_max = max(L_LSD, L);

	// stores the decoding result of Nv horizontal codewords
    int** sub_codeword_list = new int* [Nv];
	for (int i = 0; i < Nv; i++) { sub_codeword_list[i] = new int[L_max * Nh]; }
	int** sub_codeword_list_v = new int* [Nh];
	for (int i = 0; i < Nh; i++) { sub_codeword_list_v[i] = new int[L_max * Nv]; }
	float** edistance_list = new float* [Nv];
	for (int i = 0; i < Nv; i++) { edistance_list[i] = new float[L_max]; }
	float** edistance_list_v = new float* [Nh];
	for (int i = 0; i < Nh; i++) edistance_list_v[i] = new float[L_max];
	int** sub_encoded_codeword_list = new int* [Nv];
	for (int i = 0; i < Nv; i++) { sub_encoded_codeword_list[i] = new int[L_max * Nh]; }
	/*int** sub_encoded_codeword_list_v = new int* [Nh];
	for (int i = 0; i < Nh; i++) { sub_encoded_codeword_list_v[i] = new int[L_max * Nv]; }*/
	int** sub_encoded_codeword_list_v = new int* [Nh];
	for (int i = 0; i < Nh; i++) { sub_encoded_codeword_list_v[i] = new int[L_max * Nv]; }
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> H_real = MatrixXf::Random(2 * Nr, 2 * Nt); // Real domain channel matrix
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> received_signal = MatrixXf::Random(2 * Nr,1);
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> HTH = MatrixXf::Random(2 * Nt, 2 * Nt);
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> Identity = MatrixXf::Random(2 * Nt, 2 * Nt);
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> output_signal = MatrixXf::Random(2 * Nt, 1);
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> mmse_matrix = MatrixXf::Random(2 * Nt, 2 * Nr);


	// Matrix Definitions
	//Eigen::Matrix<float, 2 * Nr, 2 * Nt> H_real; // ʵ���ŵ�����
	//Eigen::Matrix<float, 2 * Nr, 1> received_signal; // �����źž���
	//Eigen::Matrix<float, 2 * Nt, 2 * Nt> HTH; // HTH����
	//Eigen::Matrix<float, 2 * Nt, 2 * Nt> Identity = Eigen::Matrix<float, 2 * Nt, 2 * Nt>::Identity(); // ��λ����
	//Eigen::Matrix<float, 2 * Nt, 1> output_signal; // ����źž���
	//Eigen::Matrix<float, 2 * Nt, 2 * Nr> mmse_matrix; // MMSE����


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
	// srand((unsigned int)time(0));
	std::random_device rd;              // ��ȡ���������
	std::mt19937 gen(rd());             // ʹ�����ӳ�ʼ��mt19937������
	std::uniform_int_distribution<> distrib(1, 1000000);
	int Num_Error, Num_Frame_Error;
	float Decoding_Duration;

	int Num_Error_new, Num_Frame_Error_new;
	float Decoding_Duration_new;

	clock_t T_start, T_finish;
	int SNR_count = 0;
	int anssc, anssd;
	for (double nSNR = Min_SNR; nSNR <= Max_SNR; nSNR += SNR_step)
	{
		for (int i = 0; i < softiter; i++) { n_tried_patterns_h_all[i] = 0; n_tried_patterns_v_all[i] = 0; }
		int* n_tried_patterns_h = new int[Nv];  
		int* n_tried_patterns_v = new int[Nh];
		for (int i = 0; i < Nv; i++) { n_tried_patterns_h[i] = 0;  }
		for (int i = 0; i < Nh; i++) { n_tried_patterns_v[i] = 0; }
		// crc check logs
		int n_err_frame_pass_crc = 0;  // err frames not detected by crc
		int n_err_frame_corrected_by_crc = 0; // err frames detected and corrected by crc
		int n_err_frame_uncorrected_by_crc = 0; // err frames detected but not corrected by crc
		int n_err_replacement_codeword_crc = 0;  // err of replacement codewords
		int n_err_pass_parity_check = 0; // n of err that passes parity check
		int n_err_frame_withno_replacement_codeword = 0; // n of err frames with no possible replacement codeword from the sub-lists
		int idx_replacement_codeword = 0; // idx of replacement codeword that passes crc check

		for (int i = 0; i < Nh; i++) { n_err_codeword_v[i] = 0; }

		for (int i = 0; i < Nh; i++)
		{
			for (int j = 0; j < L + 1; j++) { n_err_codeword_v_diffpos[i][j] = 0; }
		}
		for (int i = 0; i < L+1; i++) { max_L_cnt[i] = 0; }
		for (int i = 0; i < Nh; i++) { correct_idx_oneframe[i] = 0; }

		for (int i = 0; i < Nv; i++) sum_sdis_array[i] = 0;
		double sdis_sum = 0;
		// TEST: det err
		int n_det_err = 0; // n error bits at detector output
		int n_det_err_mrb = 0; // error bits at detector output(MRBs)
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

		for (int i = 0; i < Nh; i++)
		{
			cout << A_Ach[i] << ",";
		}
		cout << endl;
		// test: Horizontal frozen bits

		/*cout << "Frozen Bits" <<endl;
		for (int i = 0; i < Nh - Kh; i++)
		{
			cout << setw(5) << Acv[i];
		}
		cout << endl;*/
		// test: Horizontal frozen bits
		// Frame Loop
		int passnum;
		int frame = 0;

		// generate interleaving patterns
		if (flag_interleave)
		{
			string ofile_intrlv_name = string("intrlv_pattern_71_frame_49263.txt");
			string ofile_spatial_intrlv_name = string("intrlv_pattern_spatial71.txt");
			ifstream ofile_intrlv(ofile_intrlv_name);
			ifstream ofile_spatial_intrlv(ofile_spatial_intrlv_name);
			int file_in = 1;
			//for (int i = 0; i < Nv; i++)
			//{
			//	for (int j = 0; j < intrlv_pattern_num; j++)
			//	{

			//		for (int k = 0; k < ModType; k++)
			//			// ofile_intrlv >> file_in;
			//			ofile_intrlv >> intrlv_pattern_2D_time_partial[i][j][k];
			//		//cout << endl;
			//	}
			//}
			//ofile_intrlv.close();

			for (int i = 0; i < Nh; i++)
			{
				for (int j = 0; j < Nv; j++)
					ofile_spatial_intrlv >> intrlv_pattern_2D[i][j];
			}
			ofile_spatial_intrlv.close();
			/*for (int i = 0; i < Nh; i++)
			{
				for (int j = 0; j < Nv; j++)
					cout << intrlv_pattern_2D[i][j] << endl;
			}*/
			generateIntrlvPattern(intrlv_pattern, Nv);
			for (int i = 0; i < Nh; i++) { generateIntrlvPattern(intrlv_pattern_2D[i], Nv); }
			for (int i = 0; i < Nv; i++) { generateIntrlvPattern(intrlv_pattern_2D_time[i], Nh); }
			if (flag_interleave_partial)
			{
				for (int i = 0; i < Nv; i++)
				{
					for (int j = 0; j < intrlv_pattern_num; j++)
					{
						generateIntrlvPattern(intrlv_pattern_2D_time_partial[i][j], ModType);
						/*for (int k = 0; k < ModType; k++)
							cout << intrlv_pattern_2D_time_partial[i][j][k];
						cout << endl;*/
					}
				}
			}
			for (int i = 0; i < Nh; i++) { generateIntrlvPattern(intrlv_pattern_2D_inner[i], Kv); }
		}

		//ofstream ofile_intrlv("intrlv_pattern.txt");
		//for (int i=0; i<Nv; i++)
		//{
		//	for (int j = 0; j < intrlv_pattern_num; j+for+)
		//	{
		//		
		//		for (int k = 0; k < ModType; k++)
		//			ofile_intrlv << intrlv_pattern_2D_time_partial[i][j][k] << "   ";
		//		//cout << endl;
		//	}
		//}
		//ofile_intrlv.close();
		int frame_tmp = 0;
		int frame_delta = 0;
		int n_frame_error_pre = 0;
		while (Num_Frame_Error < errframe)
		{
			for (int i = 0; i < Nv; i++) { n_tried_patterns_h[i] = 0; }
			for (int i = 0; i < Nh; i++) { n_tried_patterns_v[i] = 0; }
			/*test: processing time*/
			//auto start = std::chrono::system_clock::now();
			/*if (frame == 10) {
				cout << "AAAAAAAAAA" << endl;
			}*/
			frame++;
			bool flag_pass_crc = 0;
			bool flag_pass_crc_undetected = 0;
			bool flag_pass_crc_corrected = 0;
			bool flag_found_replacement_codeword = 0;
			bool flag_pass_parity_check_in_crc = 0;

			generateIntrlvPattern(intrlv_pattern, Nv);
			for (int i = 0; i < Nh; i++) { generateIntrlvPattern(intrlv_pattern_2D[i], Nv); }
			for (int i = 0; i < Nv; i++) { generateIntrlvPattern(intrlv_pattern_2D_time[i], Nh); }
			if (flag_interleave_partial)
			{
				for (int i = 0; i < Nv; i++)
				{
					for (int j = 0; j < intrlv_pattern_num; j++)
					{
						generateIntrlvPattern(intrlv_pattern_2D_time_partial[i][j], ModType);
						/*for (int k = 0; k < ModType; k++)
							cout << intrlv_pattern_2D_time_partial[i][j][k];
						cout << endl;*/
					}
				}
			}
			for (int i = 0; i < Nh; i++) { generateIntrlvPattern(intrlv_pattern_2D_inner[i], Kv); }
			// generate interleaving pattern
			/*if (flag_interleave)
			{
				generateIntrlvPattern(intrlv_pattern, Nv);
				for (int i = 0; i < Nh; i++) { generateIntrlvPattern(intrlv_pattern_2D[i], Nv); }
				for (int i = 0; i < Nv; i++) { generateIntrlvPattern(intrlv_pattern_2D_time[i], Nh); }
				if (flag_interleave_partial)
				{
					for (int i = 0; i < Nv; i++)
					{
						for (int j = 0; j < intrlv_pattern_num; j++)
						{
							generateIntrlvPattern(intrlv_pattern_2D_time_partial[i][j], ModType);
						}
					}
				}
			}*/

			//cout << frame << endl;
			int u_input[16] = { 0,     1     ,0,     1,
									0,     1,     0,     0,
									0 ,    1,     0,     0,
									1 ,    0,     0,     0 };
			for (int i = 0; i < N_Info; i++)
			{
				a_data[i] = (unsigned char)(distrib(gen) % 2);
				//		a_data[i] = u_input[i];
			}

			//cout << endl;
			if (CRCL) tx_append_crc(a_data, N_Info, CRCL);
			for (int i = 0; i < Nv; i++)
				for (int j = 0; j < Nh; j++)
					u[i][j] = 0;
			//omp_set_num_threads(4);
//#pragma omp parallel for
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
			// test:original info
			//cout << "original info u:" << endl;
			////for (int j = 0; j < Nh; j++) { cout << u[Av[0]][j] << " "; }
			//for (int i = 0; i < Nv; i++)
			//{
			//	cout << endl;
			//	for (int j = 0; j < Nh; j++)
			//		cout << u[i][j];
			//}
			//cout << endl;
			// Polar Encoding
//#pragma omp parallel for

			// test: original info
			/*cout << "original info u:" << endl;
			display_array(u, Nv, Nh);*/
			// test: original info

			for (int i = 0; i < Nv; i++)
				PolarEncode_xor(midx[i], u[i], Nh);
			// test:horizontal encoding
			//cout << "u after horizontal encoding:" << endl;
			////for (int j = 0; j < Nh; j++) { cout << midx[Av[0]][j] << " "; }
			//for (int i = 0; i < Nv; i++)
			//{
			//	cout << endl;
			//	for (int j = 0; j < Nh; j++)
			//		cout << midx[i][j];
			//}
			//cout << endl;
			//cout << endl;
			// test:horizontal encoding
			// 
			// column-wise interleaving

			// test: row encoding
			/*cout << "u after horizontal encoding:" << endl;
			display_array(midx, Nv, Nh);*/
			// test: row encoding
			if (flag_interleave_inner)
			{
				int* tmp_to_interleave = new int[Nv];
				for (int i = 0; i < Nh; i++)
				{
					for (int j = 0; j < Nv; j++) { tmp_to_interleave[j] = midx[j][i]; }
					randIntrlvPartial(intrlv_pattern_2D_inner[i], tmp_to_interleave, Av, Nv, Kv);
					for (int j = 0; j < Nv; j++) { midx[j][i] = tmp_to_interleave[j]; }
				}
				delete[] tmp_to_interleave;
			}

			// test: innner interleaving
			/*cout << "u after inner interleaving:" << endl;
			display_array(midx, Nv, Nh);*/
			// test: innner interleaving
//#pragma omp parallel for
			for (int i = 0; i < Nh; i++)
			{
				int* tmpxv = new int[Nv]; int* tmpuv = new int[Nv];
				for (int j = 0; j < Nv; j++)
					tmpuv[j] = midx[j][i];
				PolarEncode_xor(tmpxv, tmpuv, Nv);
				for (int j = 0; j < Nv; j++)
				{
					x[j][i] = tmpxv[j];
					x_encoded_copy[j][i] = x[j][i];
				}
				delete[] tmpxv; delete[] tmpuv;
			}

			if (atten == 2)
			{
				// modulation
				// #pragma omp parallel for
				// Interleaving
				if (flag_interleave)
				{
					if (flag_interleave_all)  // interleave temporally
					{
						for (int i = 0; i < Nv; i++)
						{
							randIntrlv(intrlv_pattern_2D_time[i], x[i], Nh);
						}
					}
					//// test: to be interleaved
					//cout << "x before interleaving:" << endl;
					//// print x line by line, add spaces between each 4 terms
					//for (int i = 0; i < Nv; i++)
					//{
					//	for (int j = 0; j < Nh; j++)
					//	{
					//		if (!(j % 4)) cout << " ";
					//		cout << x[i][j] << " ";
					//	}
					//	cout << endl;
					//}

					if (flag_interleave_partial)
					{
						for (int i = 0; i < Nv; i++)
						{
							for (int j = 0; j < intrlv_pattern_num; j++)
							{
								randIntrlv(intrlv_pattern_2D_time_partial[i][j], x[i] + j * ModType, ModType);
							}
						}
					}

					//// test: to be interleaved
					//cout << "x after interleaving:" << endl;
					//// print x line by line, add spaces between each 4 terms
					//for (int i = 0; i < Nv; i++)
					//{
					//	for (int j = 0; j < Nh; j++)
					//	{
					//		if (!(j % 4)) cout << " ";
					//		cout << x[i][j] << " ";
					//	}
					//	cout << endl;
					//}
					if ((!flag_interleave_block) && flag_interleave_spatial)
					{
						for (int j = 0; j < Nh; j++)
						{
							for (int i = 0; i < Nv; i++) { x_tmp[i] = x[i][j]; }
							randIntrlv(intrlv_pattern_2D[j], x_tmp, Nv);
							for (int i = 0; i < Nv; i++) { x[i][j] = x_tmp[i]; }
						}
					}


					if (flag_interleave_block)
					{
						randIntrlvBlock(x, Nh, Nv, ModType, interleave_rows);
					}
				}
				if (MOD_DIM == "Spatial")
				{
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
						/*for (int i = 0; i < mimo_blk_height; i++)
						{
							cout << x_mod_tmp[i] << endl;
						}
						cout << "111111" << endl;*/
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
					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < 2 * Nr; j++) { y_mod_real_temp_omp[i][j] = y_mod_real[j][i]; }
					}
					// for KBest Detection
					Eigen::HouseholderQR<Eigen::MatrixXf> qr_Hreal(H_real);
					Q = qr_Hreal.householderQ();
					R = qr_Hreal.matrixQR().triangularView<Eigen::Upper>();
					//omp_set_num_threads(8);

#pragma omp parallel for private(z)
					for (int j = 0; j < Nh; j++)
					{
						MatrixXf received_signal_omp(2 * Nr, 1);
						int k_kbest_extend = pow(2, ModType / 2) * K_KBEST;
						float** paths_kbest = new float* [k_kbest_extend];
						for (int i = 0; i < k_kbest_extend; i++) paths_kbest[i] = new float[2 * Nt];
						float* ED = new float[K_KBEST];

						if (MIMODET == "MMSE") { MMSEdetection(H_real, Identity, received_signal, mmse_matrix, output_signal, HTH, y_mod_real_temp_omp[j], y_mmse_temp_omp[j], sigma_sq); }
						// ִ��һЩ����

						for (int i = 0; i < 2 * Nr; i++) { received_signal_omp(i, 0) = y_mod_real_temp_omp[j][i]; }
						z = Q.transpose() * received_signal_omp;
						if (MIMODET == "KBEST")
						{
							KBestDetection(R, z, H_real, y_mod_real_temp_omp[j], ED, y_mmse_temp_omp[j], paths_kbest, K_KBEST, ModType, Nt, Nr);
							llr_kbest_bit(Nt, y_mmse_temp_omp[j], oriLLR_kbest[j], sigma_sq, ModType, K_KBEST, paths_kbest, ED);
						}
						for (int i_path = 0; i_path < k_kbest_extend; i_path++) { delete[] paths_kbest[i_path]; }

						delete[] paths_kbest;
						delete[] ED;
					}

					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < 2 * Nt; j++) { y_mmse[j][i] = y_mmse_temp_omp[i][j]; }
					}
					// sym llr to bitwise llr
					// Detection is done for every column, that is, for every Nv bits
					// Interleaving is not currently considered
					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < 2 * Nt; j++) { y_mmse_tmp[j] = y_mmse[j][i]; }
						if (MIMODET == "MMSE") { llr_mmse_bit(Nt, y_mmse_tmp, oriLLR_tmp, sigma_sq, ModType); }
						if (MIMODET == "KBEST") { for (int j = 0; j < len_ori_llr_tmp; j++) { oriLLR_tmp[j] = oriLLR_kbest[i][j]; } }
						if (flag_interleave)
						{
							randDeintrlv(intrlv_pattern, oriLLR_tmp, Nv);
						}
						bool flagIsMax = false;
						for (int j = 0; j < Nv; j++)
						{
							if (!(j % 2))
							{
								if (abs(x_mmse_out[j][i]) > abs(x_mmse_out[j][i + 1])) { flagIsMax = true; }
								else { flagIsMax = false; }
							}
							else
							{
								if (abs(x_mmse_out[j][i]) > abs(x_mmse_out[j][i - 1])) { flagIsMax = true; }
								else { flagIsMax = false; }
							}
							oriLLRh[j][i] = oriLLR_tmp[j];
							upLLR[j][i] = oriLLR_tmp[j];
							// TEST�� det error
							x_mmse_out[j][i] = 0;
							if (oriLLR_tmp[j] < 0) { x_mmse_out[j][i] = 1; }
							if (x_mmse_out[j][i] != x[j][i]) {
								n_det_err++;
								// error at max-reliable bit. Warning: Only fit for 16QAM modulation
								if (flagIsMax) { n_det_err_mrb++; }
							}
							total_bit_cnt++;
						}
					}
					// test : detection TIME
					/*det_end = std::chrono::system_clock::now();
					auto duration_detection = std::chrono::duration_cast<std::chrono::microseconds>(det_end - start);
					file_det << duration_detection.count() << endl;*/
				}
				// Temporal modulation
				else if (MOD_DIM == "Temporal")
				{
					//modulation
					//cout << "START TEST" << endl;
					for (int i = 0; i < Nv; i++) // every row
					{
						modulation(x[i], x_mod[i], ModType, mimo_blk_length);
					}
					/*for (int i = 0; i < mimo_blk_height; i++)
					{
						for (int j = 0; j < mimo_blk_length; j++)
							cout << x_mod[i][j] << endl;
					}7
					cout << "111111" << endl;*/
					//AWGN channel (received y_mod)
					//// test: modulated symbols
					//cout << endl;
					//cout << "modulated symbols" << endl;
					//for (int i = 0; i < Nt; i++) cout << x_mod[i][0] << endl;
					//cout << "modulated symbols" << endl;
					//// test: modulated symbols
					H = generateChannelMatrix(Nr, Nt);
					for (int j = 0; j < Nsym; j++)
					{
						// Extracting a single column
						for (int i = 0; i < Nv; i++) { x_mod_tmp[i] = x_mod[i][j]; }

						AWGN(Nr, Nt, x_mod_tmp, y_mod_tmp, H, sigma_sq);
						for (int i = 0; i < Nr; i++) { y_mod[i][j] = y_mod_tmp[i]; }
					}
					convertHToReal(H, H_real);

					for (int j = 0; j < Nsym; j++)
					{
						for (int i = 0; i < Nr; i++) { y_mod_tmp[i] = y_mod[i][j]; }
						convertYToReal(Nr, y_mod_tmp, y_mod_real_tmp);
						for (int i = 0; i < 2 * Nr; i++) { y_mod_real[i][j] = y_mod_real_tmp[i]; }
					}
					//MMSE detection done in paralell
					for (int i = 0; i < Nsym; i++)
					{
						for (int j = 0; j < 2 * Nr; j++) { y_mod_real_temp_omp[i][j] = y_mod_real[j][i]; }
						// test : y_mod_real_temp_omp[0]
						/*if (i == 0)
						{
							cout << "y" << endl;
							for (int j = 0; j < 2 * Nr; j++) cout << y_mod_real_temp_omp[0][j] << endl;
							cout << "y" << endl;
						}*/
					}
					Eigen::HouseholderQR<Eigen::MatrixXf> qr_Hreal(H_real);
					Q = qr_Hreal.householderQ();
					R = qr_Hreal.matrixQR().triangularView<Eigen::Upper>();
					// for EP detection
					// copy H_real
					for (int i = 0; i < Nr2; i++)
					{
						for (int j = 0; j < Nt2; j++)
						{
							Hreal_array[i * Nt2 + j] = H_real(i, j);
						}
					}
					// calculate HTH
					cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, Nt2, Nt2, Nr2, 1.0, Hreal_array, Nt2, Hreal_array, Nt2, 0.0, HTH_array, Nt2);
					//cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, Nt2, 1, Nr2, 1.0, ch_mtx_struct.H_real[cnt], Nt2, data_struct.rx_buffer_real[cnt], 1, 0.0, ch_mtx_struct.hty, 1);

					//omp_set_num_threads(8);
					//cout << 1 << endl;
					//omp_set_num_threads(4);
					/*auto start = std::chrono::system_clock::now();*/
					/*MatrixXf H_real_trans(2*Nt, 2*Nr);
					for (int i = 0; i < 2 * Nt; i++)
						for (int j = 0; j < 2 * Nr; j++)
							H_real_trans(i, j) = H_real(j, i);
					HTHIinv = (H_real_trans * H_real + Identity).inverse();*/
# pragma omp parallel for private(z)
					for (int j = 0; j < Nsym; j++)

						//for (int x=0; x<Nsym; x++)
					{
						// int j = x % Nsym;
						MatrixXf received_signal_omp(2 * Nr, 1);
						int k_kbest_extend = pow(2, ModType / 2) * K_KBEST;
						float** paths_kbest = new float* [k_kbest_extend];
						for (int i = 0; i < k_kbest_extend; i++) { paths_kbest[i] = new float[2 * Nt]; }
						float* ED = new float[K_KBEST];
						float* HTY_array_tmp = new float[2 * Nt];
						// Declared arrays for EP
						int flag = 1;
						int temp = 0;
						//float* extr_mean = new float[2 * Nr];
						//float* extr_var = new float[2 * Nr];
						double* vec_alpha_in = new double[2 * Nr];
						double* vec_gamma_in = new double[2 * Nr];
						double* p_symbol_rearrange = new double[Nt2 * pow(2, ModType / 2)];
						for (int i = 0; i < Nt2 * pow(2, ModType / 2); i++) { p_symbol_rearrange[i] = 1 / (pow(2, ModType / 2)); }
						for (int i = 0; i < Nr2; i++)
						{
							//extr_mean[i] = 0;
							//extr_var[i] = 0;
							vec_alpha_in[i] = 2;
							vec_gamma_in[i] = 0;
						}

						if (MIMODET == "MMSE")
						{
							MMSEdetection(H_real, Identity, received_signal, mmse_matrix, output_signal, HTH, y_mod_real_temp_omp[j], y_mmse_temp_omp[j], sigma_sq);
							//MMSEdetection_noInverse(H_real, Identity, received_signal, mmse_matrix, output_signal, HTH, HTHIinv, y_mod_real_temp_omp[j], y_mmse_temp_omp[j], sigma_sq);
						}
						if (MIMODET == "KBEST")
						{
							for (int i = 0; i < 2 * Nr; i++) { received_signal_omp(i, 0) = y_mod_real_temp_omp[j][i]; }
							z = Q.transpose() * received_signal_omp;
							KBestDetection(R, z, H_real, y_mod_real_temp_omp[j], ED, y_mmse_temp_omp[j], paths_kbest, K_KBEST, ModType, Nt, Nr);
							llr_kbest_bit(Nt, y_mmse_temp_omp[j], oriLLR_kbest[j], sigma_sq, ModType, K_KBEST, paths_kbest, ED);
						}
						if (MIMODET == "EP")
						{
							// test: save H_real, recv_signal
							/*save_mat_txt(y_mod_real_temp_omp[j], Nr2, "RxSig0128.txt");
							save_mat_txt(H_real, "Hreal0128.txt");*/
							// test: save H_real, recv_signal
							cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, Nt2, 1, Nr2, 1.0, Hreal_array, Nt2, y_mod_real_temp_omp[j], 1, 0.0, HTY_array_tmp, 1);
							//test:HTH_array
							/*for (int i = 0; i < Nt2; i++)
								for (int s = 0; s < Nt2; s++)
								{
									if (s == 0) std::cout << endl;
									std::cout << setw(10) << HTH_array[i * Nt2 + s];
								}
							cout << endl;*/
							// TEST:HTY
							/*for (int i = 0; i < Nr2; i++)
								received_signal_omp(i, 0) = y_mod_real_temp_omp[j][i];
							MatrixXf HTransY = H_real.transpose() * received_signal_omp;
							cout << "Eigen" << endl;
							for (int i = 0; i < Nt2; i++) cout << setw(10)<<HTransY(i, 0);
							cout << "mkl" << endl;
							for (int i = 0; i < Nt2; i++) cout << HTY_array[i];*/
							// TEST:HTY
							EPD_detection_float(flag, Nt, Nr, iter_MIMO, ModType, 1.0, sigma_sq, HTH_array, HTY_array_tmp, temp, ep_extr_mean_omp[j], ep_extr_var_omp[j], vec_alpha_in, vec_gamma_in, p_symbol_rearrange);
						}
						//test EP detection
						/*cout << "TEST:EP DETECTION" << endl;
						for (int i = 0; i < Nr2; i++) std::cout << extr_mean[i]<< endl;
						cout << "TEST:EP DETECTION" << endl;*/
						//test EP detection
						for (int i_path = 0; i_path < k_kbest_extend; i_path++) { delete[] paths_kbest[i_path]; }
						delete[] paths_kbest;
						delete[] ED;
						//delete[] extr_mean;
						//delete[] extr_var;
						delete[] vec_alpha_in;
						delete[] vec_gamma_in;
						delete[] p_symbol_rearrange;
						delete[] HTY_array_tmp;
					}
					// test: detection time
					/*auto detection_over = std::chrono::system_clock::now();
					auto duration_detection = std::chrono::duration_cast<std::chrono::microseconds>(detection_over - start);
					cout << endl;
					cout <<"Detection duration" << duration_detection.count() << "us" << endl;*/
					// test: detection time
					for (int i = 0; i < Nsym; i++)
					{
						for (int j = 0; j < 2 * Nt; j++) { y_mmse[j][i] = y_mmse_temp_omp[i][j]; }
					}
					//To bit domain output
					for (int j = 0; j < Nsym; j++)
					{
						for (int i = 0; i < 2 * Nt; i++) { y_mmse_tmp[i] = y_mmse[i][j]; }
						// With Spatial Modulation, Nt = Nv
						if (MIMODET == "MMSE") { llr_mmse_bit(Nt, y_mmse_tmp, oriLLR_tmp, sigma_sq, ModType); }
						if (MIMODET == "EP") { llr_ep_bit(Nt, ep_extr_mean_omp[j], ep_extr_var_omp[j], oriLLR_tmp, ModType); }
						if (MIMODET == "KBEST") { for (int i = 0; i < len_ori_llr_tmp; i++) { oriLLR_tmp[i] = oriLLR_kbest[j][i]; } }
						int i_bit = 0;
						int bit_col = 0;
						for (int j_inner = 0; j_inner < Nv; j_inner++)
						{
							for (int i_bit_col = 0; i_bit_col < ModType; i_bit_col++)
							{
								i_bit = ModType * j_inner + i_bit_col;
								bit_col = j * ModType + i_bit_col;
								oriLLRh[j_inner][bit_col] = oriLLR_tmp[i_bit];
								upLLR[j_inner][bit_col] = oriLLR_tmp[i_bit];
							}

						}
					}
					// test: LLR transfer time
					/*auto llr_to_bit_over = std::chrono::system_clock::now();
					auto duration_llr = std::chrono::duration_cast<std::chrono::microseconds>(llr_to_bit_over - detection_over);
					cout << endl;
					cout<< "LLR transfer duration" << duration_llr.count() << "us" << endl;*/
					// test: LLR transfer time
					if (flag_interleave)
					{
						if ((!flag_interleave_block) && flag_interleave_spatial)
						{
							for (int j = 0; j < Nh; j++)
							{
								for (int i = 0; i < Nv; i++) { oriLLR_tmp[i] = oriLLRh[i][j]; }
								randDeintrlv(intrlv_pattern_2D[j], oriLLR_tmp, Nv);
								for (int i = 0; i < Nv; i++)
								{
									oriLLRh[i][j] = oriLLR_tmp[i];
									upLLR[i][j] = oriLLR_tmp[i];
								}

							}
						}

						if (flag_interleave_block)
						{
							randDeIntrlvBlock(oriLLRh, Nh, Nv, ModType, interleave_rows);
							randDeIntrlvBlock(upLLR, Nh, Nv, ModType, interleave_rows);
						}

						if (flag_interleave_all) // Deinterleave in time domain
						{
							for (int j = 0; j < Nv; j++)
							{
								randDeintrlv(intrlv_pattern_2D_time[j], oriLLRh[j], Nh);
								randDeintrlv(intrlv_pattern_2D_time[j], upLLR[j], Nh);
							}
						}

						if (flag_interleave_partial)
						{
							for (int j = 0; j < Nv; j++)
							{
								for (int k = 0; k < intrlv_pattern_num; k++)
								{
									randDeintrlv(intrlv_pattern_2D_time_partial[j][k], oriLLRh[j] + k * ModType, ModType);
									randDeintrlv(intrlv_pattern_2D_time_partial[j][k], upLLR[j] + k * ModType, ModType);
								}
							}
						}

					}

					for (int i = 0; i < Nv; i++)
					{
						for (int j = 0; j < Nh; j++)
						{
							total_bit_cnt++;
							int bit_dec = (oriLLRh[i][j] < 0);
							if (bit_dec != x[i][j]) n_det_err++;
						}
					}

				}
				// test: detection errors

				// test: detection errors
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
			for (int iter = 1; iter <= softiter; iter++)
			{
				for (int i = 0; i < Nv; i++) { n_tried_patterns_h[i] = 0; }
				for (int i = 0; i < Nh; i++) { n_tried_patterns_v[i] = 0; }
				//cout << '1' << endl;
				//vertical decoding
		   // test: decoding time
			//auto decoding_start = std::chrono::system_clock::now();

			// test: decoding time
#pragma omp parallel for
				for (int decodingi = 0; decodingi < Nh; decodingi++) // 291.70kb
				{

					//for SCL Decoding
					int Counth = 0, Countinfoh = 0;
					float Qsumunh = 0;
					int** sum = new int* [L];
					//for (int i = 0; i < L; i++)  sum[i] = new  int[Nmax];
					int** u1 = new int* [L];
					for (int i = 0; i < L; i++)  u1[i] = new _MM_ALIGN16 int[Nmax];
					for (int i = 0; i < L; i++)  sum[i] = new _MM_ALIGN16 int[Nmax];
					//for (int i = 0; i < L; i++)  u1[i] = new  int[Nmax];

					float* LLR_in = new float[L];
					float* W = new float[2 * L];
					int* SCL_index = new int[L];
					int* better = new int[L];
					int* worse = new int[L];
					int* indexpm = new int[L];
					float* PM = new float[L];
					float** LLR = new float* [L];
					for (int i = 0; i < L; i++)  LLR[i] = new _MM_ALIGN16 float[2 * Nmax + 2];
					// for (int i = 0; i < L; i++)  LLR[i] = new float[2 * Nmax + 2];
					////////////////////////////////
					for (int j = 0; j < Nv; j++)
						for (int k = 0; k < L; k++)
							LLR[k][j] = upLLR[j][decodingi];
					//cout << LLR[0][j] << " ";

				//cout << endl;
					int index = 0;
					/*			if (L == 1)
								{
									//T_start = clock();
									Count = 0;
									decode(*LLR, *sum, (int)(log(Nv) / log(2)), A_Acv, flag, Nv, Kv);
									//PolarEncode_xor(*u1, *sum, N);
									//T_finish = clock();
								}
								else
								{*/
					for (int k = 0; k < L; k++) { PM[k] = 0; indexpm[k] = k; }
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
					if (adaptive_beta)
					{
						if (softmethod == "Pyndiah") Pyndiahv_adaptivebeta(LLR, sum, oriLLRh, indexpm, iter * 2 - 2, decodingi, upLLR);
						if (softmethod == "MITSO") MITSOv(LLR, PM, sum, oriLLRh, indexpm, iter * 2 - 2, decodingi, upLLR, Qsumunh);
					}
					// iter * 2 - 2
					else
					{
						//cout << iter << endl;
						if (softmethod == "Pyndiah") {
							if (!bit_flip_pyndiah | (!vertical_ml)) { Pyndiahv(LLR, sum, oriLLRh, indexpm, iter, decodingi, upLLR); }
							else { Pyndiahv_bitflip(LLR, sum, oriLLRh, indexpm, iter, decodingi, upLLR, sum_sdis_array[decodingi], tried_codewords_v[decodingi], GGv, Acv, n_tried_patterns_v); }
						}
							
						if (softmethod == "MITSO") MITSOv(LLR, PM, sum, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunh);
					}
					// test: oriLLR and upLLR
					/*cout << endl;
					for (int i = 0; i < Nv; i++)
						for (int j = 0; j < Nh; j++)
							cout << oriLLRh[i][j] - upLLR[i][j] << "   ";
					cout << endl;*/
					/*cout << oriLLRh[0][0] << endl;
					cout << upLLR[0][0] << endl;*/
					// test: oriLLR and upLLR
					// half iteration Enabled
					if (half_iter && iter == softiter)
					{
						if (iter == softiter) PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), Nh);
						for (int j = 0; j < Kv; j++)
							umid[Av[j]][decodingi] = u1[index][Av[j]];
						for (int i = 0; i < L; i++)
						{
							int curr_index = indexpm[i];
							for (int j = 0; j < Nv; j++)
							{
								sub_encoded_codeword_list_v[decodingi][i * Nv + j] = sum[curr_index][j];
							}
							PolarEncode_xor(sub_codeword_list_v[decodingi] + i * Nv, sum[curr_index], Nv);
							edistance_list_v[decodingi][i] = PM[curr_index];
						}
						// check wheher an item of sum equals to x_encoded_copy[:][decodingi]
						int best_idx = indexpm[0];
						/*cout << endl;
						for (int k = 0; k < Nh; k++)
						{
							cout << sum[best_idx][k] << "  " << x_encoded_copy[k][decodingi] << endl;
						}*/


						int correct_codeword_idx = 0;
						int correct_idx_pre = -1;
						for (int i = 0; i < L; i++)
						{
							// cout << 1 << endl;
							int curr_index = indexpm[i];
							int found_correct = 0;
							for (int j = 0; j < Nv; j++)
							{
								if (sum[curr_index][j] != x_encoded_copy[j][decodingi]) { break; }
								if (j == Nv - 1) { correct_codeword_idx = curr_index; found_correct = 1; }
							}

							// display x_encoded_copy
				
							if (found_correct)
							{
								n_err_codeword_v[decodingi] += i;
								n_err_codeword_v_diffpos[decodingi][i] += 1;
								correct_idx_oneframe[decodingi] = i;
								break;
							}
							if (i == L - 1) {
								n_err_codeword_v[decodingi] += L;
								n_err_codeword_v_diffpos[decodingi][L] += 1;
								correct_idx_oneframe[decodingi] = L;
								break;
							}
						}
					}
					// copy to sub_encoded_codeword_list_v
					/*for (int i = 0; i < L; i++)
					{
						for (int j = 0; j < Nv; j++)
						{
							sub_encoded_codeword_list_v[i][i * Nv + j] = sum[indexpm[i]][j];
						}
					}*/
					// half iteration
					////////delete SCL variables
					for (int i = 0; i < L; i++)
						delete[]LLR[i];
					delete[] LLR;
					for (int i = 0; i < L; i++)
						delete[]sum[i];
					for (int i = 0; i < L; i++)
						delete[]u1[i];
					delete[] u1;
					delete[] sum;
					delete[] PM;
					delete[] LLR_in;
					delete[] W;
					delete[] SCL_index;
					delete[] better;
					delete[] worse;
					delete[] indexpm;

				}
				//cout << "Horizon Start:\n";
				//horizontal decoding
				// #
				if (weirditer) { continue; }

				// test: LLR after horizontal decoding
				/*cout << "After horizontal decoding" << endl;
				for (int i = 0; i < Nv; i++)
				{
					for (int j = 0; j < Nh; j++)
					{
						cout << (oriLLRh[i][j] < 0) << " ";
					}
					cout << endl;
				}*/
				// test: LLR after horizontal decoding
				// inner deinterleaving
				// randDeIntrlvPartial(intrlv_pattern_2D_inner[0],);
				if (flag_interleave_inner)
				{
					float* tmp_llr_to_interleave = new float[Nv];
					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < Nv; j++) { tmp_llr_to_interleave[j] = oriLLRh[j][i]; }
						randDeIntrlvPartial(intrlv_pattern_2D_inner[i], tmp_llr_to_interleave, Av, Nv, Kv);
						for (int j = 0; j < Nv; j++) { oriLLRh[j][i] = tmp_llr_to_interleave[j]; }
					}
					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < Nv; j++) { tmp_llr_to_interleave[j] = upLLR[j][i]; }
						randDeIntrlvPartial(intrlv_pattern_2D_inner[i], tmp_llr_to_interleave, Av, Nv, Kv);
						for (int j = 0; j < Nv; j++) { upLLR[j][i] = tmp_llr_to_interleave[j]; }
					}
					delete[] tmp_llr_to_interleave;
				}

				if ((iter == softiter) && half_iter)
				{
					int max_correct_idx = -1;
					for (int i = 0; i < Nh; i++)
					{
						if (correct_idx_oneframe[i] > max_correct_idx)
						{
							max_correct_idx = correct_idx_oneframe[i];
						}
					}
					max_L_cnt[max_correct_idx]++;
				}
				/*auto decoding_over = std::chrono::system_clock::now();
				auto duration_decoding = std::chrono::duration_cast<std::chrono::microseconds>(decoding_over - decoding_start);
				cout << endl;
				cout << "Decoding duration" << duration_decoding.count() << "us" << endl;*/

#pragma omp parallel for
				for (int decodingi = 0; decodingi < Nv; decodingi++)
				{
					// 293.70kb
					// cout << "KKKK" << endl;
					if (half_iter && (iter == softiter)) { break; }
					string polarde_now = polarde_array[iter - 1];
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

					/*			if (L == 1)
								{
									//T_start = clock();
									Count = 0;
									decode(*LLR, *sum, (int)(log(Nh) / log(2)), A_Ach, flag, Nh, Kh);
									if (iter == softiter) PolarEncode_xor(*u1, *sum, Nh);
									//T_finish = clock();
								}
								else
								{*/
								//for (int k = 0; k < L; k++)  PM[k] = 0, indexpm[k] = k; 
								////T_start = clock();
								//decode_list(LLR, sum, (int)(log(Nh) / log(2)), A_Ach, 0, Nh, Kh, PM, L, (int)(log(L) / log(2)), LLR_in, W, SCL_index, better, worse, &Countv, &Countinfov, &Qsumunv);
								////T_finish = clock();
								////Decoding_Duration += (T_finish - T_start);
								//for (int i = 0; i < L - 1; i++)
								//	for (int j = i + 1; j < L; j++)
								//		if (PM[indexpm[i]] < PM[indexpm[j]])
								//		{
								//			int ttt = indexpm[i];
								//			indexpm[i] = indexpm[j];
								//			indexpm[j] = ttt;
								//		}


							/*
							* Encoding: Info -> Horizontal Encoding -> A -> set to 0 -> B -> Horizontal Encoding -> Vertical Encoding
							* Decoding: Received Info -> Horizontal Decoding -> Vertical Decoding -> Horizontal Decoding
							*/
					if (polarde_now == "SLD")
					{
						LSD_decode_outall(u_decod, u1_LSD, LLR[0], GGh, A_Ach, Nh, L_LSD, n_paths);
						// if (iter == softiter) PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), Nh);
						// for soft information exchange, inverse coding:
						// THIS LINE is WRONG! BUT IT SEEMS BETTER THAN THE RIGHT VERSION
						for (int i = 0; i < L_LSD; i++) { PolarEncode_xor(sum_LSD[i], u1_LSD[i], Nh); }
						// for (int i = 0; i < L_LSD; i++) { PolarEncode_xor(sum_LSD[i], u1_LSD[i], Nh); }
						// if (iter == softiter) PolarEncode_xor(u1[0], sum_LSD[0], Nh);
						index = indexpm[0];
						// test: after horizontal decoding
						/*if (decodingi == Av[0])
						{
							cout << "Direct output of horizontal decoding:" << endl;
							for (int i = 0; i < Nh; i++) { cout << u1_LSD[0][i] << " "; }
							cout << endl;
							cout << "horizontal decoder out, encode once:" << endl;
							for (int i = 0; i < Nh; i++) { cout << sum_LSD[0][i] << " "; }
							cout << endl;
						}*/
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
						// iter*2-1
						if (iter < softiter)
						{
							if (adaptive_beta)
							{
								if (softmethod == "Pyndiah") Pyndiahh_LSD_adaptivebeta(LLR, sum_LSD, oriLLRh, indexpm, iter, decodingi, upLLR);
								if (softmethod == "MITSO") MITSOh(LLR, PM, sum_LSD, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunv);
							}
							else
							{
								// cout << "IN LSD!!!" << endl;
								if (softmethod == "Pyndiah")
									Pyndiahh_LSD_getsum(LLR, sum_LSD, oriLLRh, indexpm, iter, decodingi, upLLR, sum_sdis_array[decodingi], GGh, Ach);
								//if (softmethod == "Pyndiah") Pyndiahh_LSD(LLR, sum_LSD, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR);
								if (softmethod == "MITSO") MITSOh(LLR, PM, sum_LSD, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunv);
							}
							//	getchar();
						}
						else
							for (int j = 0; j < Kh; j++)
								umid[decodingi][Ah[j]] = u1_LSD[index][Ah[j]];
					}
					if (polarde_now == "SCL")
					{
						for (int k = 0; k < L; k++)  PM[k] = 0, indexpm[k] = k;
						//T_start = clock();
						decode_list(LLR, sum, (int)(log(Nh) / log(2)), A_Ach, 0, Nh, Kh, PM, L, (int)(log(L) / log(2)), LLR_in, W, SCL_index, better, worse, &Countv, &Countinfov, &Qsumunv);
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
						if (iter == softiter) PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), Nh);
						index = indexpm[0];
						for (int col = 0; col < Nh; col++) { x_parity_check[decodingi][col] = sum[indexpm[0]][col]; }
						// show sum[indexpm[0]]
						/*cout << "sum[indexpm[0]]" << endl;
						for (int i = 0; i < Nh; i++) { cout << sum[indexpm[0]][i] << "  "; }
						cout << "sum[indexpm[0]]" << endl;*/
						// show sum[indexpm[0]]
						//}
						//calc distance s()
						if (iter < softiter)
						{
							if (adaptive_beta)
							{
								if (softmethod == "Pyndiah") Pyndiahh_adaptivebeta(LLR, sum, oriLLRh, indexpm, iter, decodingi, upLLR);
								if (softmethod == "MITSO") MITSOh(LLR, PM, sum, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunv);
								//	getchar();
							}
							else
							{
								//if (softmethod == "Pyndiah") Pyndiahh(LLR, sum, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR);
								if (softmethod == "Pyndiah")
								{
									if (!bit_flip_pyndiah)  Pyndiahh_getsum(LLR, sum, oriLLRh, indexpm, iter, decodingi, upLLR, sum_sdis_array[decodingi]);
									else Pyndiahh_bitflip(LLR, sum, oriLLRh, indexpm, iter, decodingi, upLLR, sum_sdis_array[decodingi], tried_codewords_h[decodingi], GGh, Ach, n_tried_patterns_h);
									
								}
								if (softmethod == "MITSO") MITSOh(LLR, PM, sum, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunv);
								//	getchar();
							}
						}
						else
						{
							// prev
							/*for (int j = 0; j < Kh; j++)
								umid[decodingi][Ah[j]] = u1[index][Ah[j]];*/

							for (int j = 0; j < Nh; j++)
								umid[decodingi][j] = u1[index][j];
							// add result of current decoding to subcodeword_list
							for (int j = 0; j < L; j++)
							{
								int curr_index = indexpm[j];
								PolarEncode_xor(sub_codeword_list[decodingi] + j * Nh, sum[curr_index], Nh);
								for (int i = 0; i < Nh; i++) { sub_encoded_codeword_list[decodingi][j * Nh + i] = sum[curr_index][i]; }
								edistance_list[decodingi][j] = PM[curr_index];
							}

							if (test_cadsl && (!half_iter)) {
								int best_idx = indexpm[0];
								/*cout << endl;
								for (int k = 0; k < Nh; k++)
								{
									cout << sum[best_idx][k] << "  " << x_encoded_copy[k][decodingi] << endl;
								}*/


								int correct_codeword_idx = 0;
								int correct_idx_pre = -1;
								for (int i = 0; i < L; i++)
								{
									// cout << 1 << endl;
									int curr_index = indexpm[i];
									int found_correct = 0;
									for (int j = 0; j < Nv; j++)
									{
										if (sum[curr_index][j] != x_encoded_copy[decodingi][j]) { break; }
										if (j == Nv - 1) { correct_codeword_idx = curr_index; found_correct = 1; }
									}

									// display x_encoded_copy

									if (found_correct)
									{
										n_err_codeword_v[decodingi] += i;
										n_err_codeword_v_diffpos[decodingi][i] += 1;
										correct_idx_oneframe[decodingi] = i;
										break;
									}
									if (i == L - 1) {
										n_err_codeword_v[decodingi] += L;
										n_err_codeword_v_diffpos[decodingi][L] += 1;
										correct_idx_oneframe[decodingi] = L;
										break;
									}
								}
							}
						}
						     
					}

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

				if (iter == softiter && (!half_iter) && test_cadsl)
				{
					int max_correct_idx = -1;
					for (int i = 0; i < Nh; i++)
					{
						if (correct_idx_oneframe[i] > max_correct_idx)
						{
							max_correct_idx = correct_idx_oneframe[i];
						}
					}
					max_L_cnt[max_correct_idx]++;
				}

				// list extension after horizontal decoding
				//if (iter == softiter && CRCL>0)
				//{
				//	crc_list_extend(sub_codeword_list, whole_codeword_list, edistance_list, L, Nv, Nh, L_WHOLE_CODEWORD);
				//	/*for (int i_subcode = 0; i_subcode < Nv; i_subcode++)
				//	{
				//		cout << endl << "subcodeword in whole codeword list" << endl;
				//		for (int i_bit = 0; i_bit < Nh; i_bit++)
				//			cout << whole_codeword_list[0][i_subcode * Nh + i_bit];
				//		cout << endl << "subcodeword in decoded result" << endl;
				//		for (int i_bit = 0; i_bit < Nh; i_bit++)
				//			cout << umid[i_subcode][i_bit];
				//	}
				//	cout << "TEST ENDS" << endl;*/
				//}

				// test: double stage crc list extend
				/*if (iter == softiter)
				{
					crc_list_extend_vertical_doublestage(sub_codeword_list_v, whole_codeword_list_v, intrlv_pattern_2D_time_partial, edistance_list_v,
						y_mod_real_temp_omp, Hreal_array, L_SUB_CODEWORD, Nh, Nv, Nr, L_WHOLE_CODEWORD, L_MIMO_BLK, ModType);
				}*/
				// test: double stage crc list extend
				





				if (flag_interleave_inner)
				{
					float* tmp_llr_to_interleave = new float[Nv];
					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < Nv; j++) { tmp_llr_to_interleave[j] = oriLLRh[j][i]; }
						randIntrlvPartial(intrlv_pattern_2D_inner[i], tmp_llr_to_interleave, Av, Nv, Kv);
						for (int j = 0; j < Nv; j++) { oriLLRh[j][i] = tmp_llr_to_interleave[j]; }
					}
					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < Nv; j++) { tmp_llr_to_interleave[j] = upLLR[j][i]; }
						randIntrlvPartial(intrlv_pattern_2D_inner[i], tmp_llr_to_interleave, Av, Nv, Kv);
						for (int j = 0; j < Nv; j++) { upLLR[j][i] = tmp_llr_to_interleave[j]; }
					}
					delete[] tmp_llr_to_interleave;
				}
				for (int i = 0; i < Nv; i++) { n_tried_patterns_h_all[iter-1] += n_tried_patterns_h[i]; }
				for (int i = 0; i < Nh; i++) { n_tried_patterns_v_all[iter-1] += n_tried_patterns_v[i]; }
			}
			int calcsum[Nv], calcu1[Nv];
			int calcsumh[Nh], calcu1h[Nh];
			// calcsum is from the direct output of the decoder
			// if output is the result of row decoding (meaning there is another column decoding to be done)
			if (half_iter)
			{
				/*if (CRCL > 0)
				{
					crc_list_extend_vertical(sub_codeword_list_v, whole_codeword_list_v, intrlv_pattern_2D_time_partial, edistance_list_v, y_mod_real_temp_omp,
						Hreal_array, L, Nh, Nv, Nr, L_WHOLE_CODEWORD, ModType);
				}*/

				// test: sub_codeword_list_v
				/*cout << "SUB_CODEWORD_LIST_V" << endl;
				for (int i = 0; i < Nh; i++)
				{
					int* encoded = new int[Nv];
					PolarEncode_xor(encoded, sub_codeword_list_v[i], Nh);
					for (int j = 0; j < Nv; j++) sub_codeword_list_v[i][j] = encoded[j];
				}
				for (int i = 0; i < Nv; i++)
				{
					for (int j = 0; j < ModType; j++)
					{
						cout << sub_codeword_list_v[j][i] << "  ";
					}
					cout << endl;
				}
				cout << "SUB_CODEWORD_LIST_V" << endl;*/
				// test: sub_codeword_list_v
				/*crc_list_extend_vertical(sub_codeword_list_v, whole_codeword_list_v, intrlv_pattern_2D_time_partial, edistance_list_v, y_mod_real_temp_omp,
					Hreal_array, L, Nh, Nv, Nr, L_WHOLE_CODEWORD, ModType);*/

				/*for (int i = 0; i < Nv; i++)
				{
					for (int j = 0; j < Nh; j++) { calcsumh[j] = whole_codeword_list_v[0][j * Nv + i]; }
					PolarEncode_xor(calcu1h, calcsumh, Nh);
					for (int j = 0; j < Nh; j++)
						uout[i][j] = calcu1h[j];
				}*/
				/*for (int i = 0; i < Kv; i++)
				{
					PolarEncode_xor(calcu1h, umid[Av[i]], Nh);
					for (int j = 0; j < Kh; j++)
						uout[Av[i]][Ah[j]] = calcu1h[Ah[j]];
				}*/
				for (int i = 0; i < Nv; i++)
				{
					PolarEncode_xor(calcu1h, umid[i], Nh);
					for (int j = 0; j < Nh; j++)
						uout[i][j] = calcu1h[j];
				}

				if (CRCL > 0 && (CRC_CHECK_ENABLE))
				{
					flag_pass_crc = false;
					flag_found_replacement_codeword = false;



					for (int i = 0; i < Kv; i++)
					{
						for (int j = 0; j < Kh; j++)
							u_crc[i * Kh + j] = uout[Av[i]][Ah[j]];
					}
					if (rx_check_crc(u_crc, Kh * Kv, CRCL)) {
						flag_pass_crc = true;
					}


					// the i th horizontal codeword , j th col
					// whole_codeword_list_v[j*Nv+i];
					if (!flag_pass_crc)
					{
						//cout << endl << "NOT PASS" << endl;
						// float* H_real_array;
						// crc_list_extend_vertical(sub_codeword_list_v, whole_codeword_list_v, intrlv_pattern_2D_time_partial, edistance_list_v, y_mod_real_temp_omp,
							//Hreal_array, L, Nh, Nv, Nr, L_WHOLE_CODEWORD, ModType);
						// prev: 
						/*crc_list_extend_vertical(sub_codeword_list_v, whole_codeword_list_v, intrlv_pattern_2D_time_partial, edistance_list_v, y_mod_real_temp_omp,
								Hreal_array, L_SUB_CODEWORD, Nh, Nv, Nr, L_WHOLE_CODEWORD, ModType);*/
						crc_list_extend_vertical_doublestage(sub_codeword_list_v, whole_codeword_list_v, intrlv_pattern_2D_time_partial, edistance_list_v,
							y_mod_real_temp_omp, Hreal_array, L_SUB_CODEWORD, Nh, Nv, Nr, L_WHOLE_CODEWORD, L_MIMO_BLK, ModType);
						int i_code = 0;
						for (i_code = 0; i_code < L_FINAL_SEARCH; i_code++) // prev: i < L_WHOLE_CODEWORD
						{
							for (int i = 0; i < Kv; i++)
							{
								for (int j = 0; j < Nh; j++)
									calcsumh[j] = whole_codeword_list_v[i_code][j * Nv + Av[i]];
								PolarEncode_xor(calcu1h, calcsumh, Nh);
								for (int j = 0; j < Kh; j++) { uout[Av[i]][Ah[j]] = calcu1h[Ah[j]]; }
							}
							// temporal encoding to check whether the codewords is valid 
							for (int i = 0; i < Nv; i++)
							{
								flag_pass_parity_check_in_crc = true;
								int* whole_cwd_tmp = new int[Nh];
								int* whole_cwd_encoded_tmp = new int[Nh];
								for (int k = 0; k < Nh; k++) { whole_cwd_tmp[k] = 0; whole_cwd_encoded_tmp[k] = 0; }
								for (int j = 0; j < Nh; j++){ whole_cwd_tmp[j] = whole_codeword_list_v[i_code][j * Nv + i]; }
								PolarEncode_xor(whole_cwd_encoded_tmp, whole_cwd_tmp, Nh);
								//cout << endl;
								for (int j = 0; j < (Nh - Kh); j++) { /*cout << whole_cwd_encoded_tmp[Ach[j]]<<"  ";*/ if (whole_cwd_encoded_tmp[Ach[j]] != 0) flag_pass_parity_check_in_crc = false; }
							    delete[] whole_cwd_tmp;
								delete[] whole_cwd_encoded_tmp;
							}
							//for (int i = 0; i < Kh; i++)
							//{
							//	for (int j = 0; j < Nv; j++)
							//		calcsum[j] = whole_codeword_list[i_code][j * Nh + Ah[i]]; // Nv*Ah[i] + j
							//	PolarEncode_xor(calcu1, calcsum, Nv);
							//	for (int j = 0; j < Kv; j++)
							//		uout[Av[j]][Ah[i]] = calcu1[Av[j]];
							//}
							// copy from uout to u_crc
							for (int i = 0; i < Kv; i++)
							{
								for (int j = 0; j < Kh; j++)
									u_crc[i * Kh + j] = uout[Av[i]][Ah[j]];
							}
							if (rx_check_crc(u_crc, Kh * Kv, CRCL)) {
								// if (i_code == 0) { flag_pass_crc = true; }
								cout << "FLAG_PAR" << flag_pass_parity_check_in_crc << endl;
								if (flag_pass_parity_check_in_crc) { // if (true)
									cout << endl << "Found in " << i_code << endl; flag_found_replacement_codeword = true; idx_replacement_codeword = i_code; 
								}
								break;
							}
							//if (i_code == 0) { break; }
						}
						/*cout << endl << i_code << endl;
						if (i_code == L_WHOLE_CODEWORD) { 
							cout << "NOT FOUND" << endl; 
						}*/
					}
				}
			}
			// if output is the result of column decoding (meaning that there is another row decoding to be done)
			else
			{
				
					// prev
					/*for (int i = 0; i < Kh; i++)
					{
						for (int j = 0; j < Nv; j++)
							calcsum[j] = umid[j][Ah[i]];
						PolarEncode_xor(calcu1, calcsum, Nv);
						for (int j = 0; j < Kv; j++)
							uout[Av[j]][Ah[i]] = calcu1[Av[j]];
					}*/

					// the whole codeword
					for (int i = 0; i < Nh; i++)
					{
						//cout << endl;
						for (int j = 0; j < Nv; j++)
						{
							calcsum[j] = umid[j][i]; 
							//cout << umid[j][i] <<"  ";
						}
						for (int i_calcu1 = 0; i_calcu1 < Nv; i_calcu1++) { calcu1[i_calcu1] = 0; }
						PolarEncode_xor(calcu1, calcsum, Nv);
						for (int j = 0; j < Nv; j++)
							uout[j][i] = calcu1[j];
					}
					//cout << endl << "END TEST UMID" << endl;
					// output uout
					/*for (int i = 0; i < Nv; i++)
					{
						cout << endl;
						for (int j = 0; j < Nh; j++) { cout << uout[i][j] << "  "; }
					}*/
					// cout << endl;
				if (CRCL > 0)
				{
	                // conduct crc check on info bits of uout
					

					flag_pass_crc = false;
					flag_found_replacement_codeword = false;

					

					for (int i = 0; i < Kv; i++)
					{
						for (int j = 0; j < Kh; j++)
							u_crc[i * Kh + j] = uout[Av[i]][Ah[j]];
					}
					if (rx_check_crc(u_crc, Kh * Kv, CRCL)) {
						flag_pass_crc = true;
					}

					if (!flag_pass_crc)
					{
						// decide whether the row codewords can be possibly combined to form a correct codeword
						bool flag_found_in_all_subcodeword = true;
						bool flag_found_in_all_subcodeword_v = true;
						for (int i_row = 0; i_row < Nv; i_row++)
						{
							int i_sub = 0;
							for (i_sub = 0; i_sub < L; i_sub++)
							{
								bool is_correct = true;
								for (int i = 0; i < Nh; i++) { if (sub_encoded_codeword_list[i_row][i_sub * Nh + i] != x_encoded_copy[i_row][i]) is_correct = false; }
								if (is_correct) break;
							}
							if (i_sub == L) { flag_found_in_all_subcodeword = false;  }
							if (!flag_found_in_all_subcodeword) break;
						}

						for (int i_col = 0; i_col < Nh; i_col++)
						{
							int i_sub = 0;
							for (i_sub = 0; i_sub < L; i_sub++)
							{
								bool is_correct = true;
								for (int i = 0; i < Nv; i++) { if (sub_encoded_codeword_list[i_col][i_sub * Nv + i] != x_encoded_copy[i][i_col]) is_correct = false; }
								if (is_correct) break;
							}
							if (i_sub == L) { flag_found_in_all_subcodeword_v = false;   }
							if (!flag_found_in_all_subcodeword_v) break;
						}

						if ((!flag_found_in_all_subcodeword) && (!flag_found_in_all_subcodeword_v))
						{
							n_err_frame_withno_replacement_codeword++;
						}

						crc_list_extend(sub_codeword_list, whole_codeword_list, edistance_list, L, Nv, Nh, L_WHOLE_CODEWORD);
						for (int i_code = 0; i_code < L_WHOLE_CODEWORD; i_code++)
						{
							for (int i = 0; i < Kh; i++)
							{
								for (int j = 0; j < Nv; j++)
									calcsum[j] = whole_codeword_list[i_code][j * Nh + Ah[i]];
								PolarEncode_xor(calcu1, calcsum, Nv);
								for (int j = 0; j < Kv; j++)
									uout[Av[j]][Ah[i]] = calcu1[Av[j]];
							}
							// copy from uout to u_crc
							for (int i = 0; i < Kv; i++)
							{
								for (int j = 0; j < Kh; j++)
									u_crc[i * Kh + j] = uout[Av[i]][Ah[j]];
							}
							if (rx_check_crc(u_crc, Kh * Kv, CRCL)) {
								if (i_code == 0) { flag_pass_crc = true; }
								if (i_code > 0) { 
									cout << endl << "Found in " << i_code << endl; flag_found_replacement_codeword = true; idx_replacement_codeword = i_code;
									// row-wise parity check to column codewords

								}

								break;
							}
							//if (i_code == 0) { break; }
						}
					}
				}
			}
			// test: final output
			/*cout << "final output" << endl;
			for (int k = 0; k < Nh; k++) { cout << uout[Av[0]][k] << " "; }
			cout << endl;*/
			// test: final output
			//now without CRC
			int Temp_Error = Num_Error;
			for (int i = 0; i < Kv; i++)
			{
				for (int j = 0; j < Kh; j++)
				{
					//printf("%d ", (int)uout[Av[i]][Ah[j]]);
					Num_Error += abs(u[Av[i]][Ah[j]] - uout[Av[i]][Ah[j]]);
					/*if (!flag_pass_crc)
					{
						cout << abs(u[Av[i]][Ah[j]] - uout[Av[i]][Ah[j]]) << endl;
					}*/
				}
				//cout << endl;
			}


			// whether pass parity check for both row and column
			bool flag_pass_parity_check = true;

			//  check with x_parity_check
			for (int i = 0; i < Nh; i++)
			{
				//cout << endl;
				for (int j = 0; j < Nv; j++)
				{
					calcsum[j] = x_parity_check[j][i];
					//cout << umid[j][i] <<"  ";
				}
				for (int i_calcu1 = 0; i_calcu1 < Nv; i_calcu1++) { calcu1[i_calcu1] = 0; }
				PolarEncode_xor(calcu1, calcsum, Nv);
				// cout << endl;
				// for (int j = 0; j < Nv; j++) { cout << calcu1[j] << "  "; }
				for (int k = 0; k < Nv - Kv; k++) { if (calcu1[Acv[k]] != 0) flag_pass_parity_check = false;  }
			}
			


			/*for (int i = 0; i < Nv - Kv; i++)
			{
				for (int j = 0; j < Nh; j++) { if (uout[Acv[i]][j] != 0) { flag_pass_parity_check = false; break; } }
				if (!flag_pass_parity_check) break;
			}
			for (int i = 0; i<Nh - Kh; i++)
			{
				for (int j = 0; j < Nv; j++) { if (uout[j][Ach[i]] != 0) { flag_pass_parity_check = false; break; } }
				if (!flag_pass_parity_check) break;
			}*/

			// cout << "PARITY" << flag_pass_parity_check << endl;

			/*if (!flag_pass_crc) {
				cout << endl;
			}*/
			//getchar();
			if (flag_found_replacement_codeword)
			{   
				file_crc_check_replacement_codeword << endl;
				if (Temp_Error < Num_Error) { file_crc_check_replacement_codeword << "Index of Replacement Codeword: " << idx_replacement_codeword << "IsCorrect: " << 0 << endl; }
				if (Temp_Error == Num_Error) { file_crc_check_replacement_codeword << "Index of Replacement Codeword: " << idx_replacement_codeword << "IsCorrect: " << 1 << endl; }
			}
			if (Temp_Error != Num_Error) // prev: <
			{
				Num_Frame_Error++;
				if (flag_pass_crc) { n_err_frame_pass_crc++; }
				else { n_err_frame_uncorrected_by_crc++; }
				if (flag_found_replacement_codeword) {
					n_err_replacement_codeword_crc++;
				}
				if (flag_pass_parity_check) { n_err_pass_parity_check++; }
			}
			if (Temp_Error == Num_Error)
			{
				/*if (!flag_pass_parity_check) { 
					cout << endl << "FFFFFFFUCKKKKKK" << endl; 
					for (int i = 0; i < Nv; i++)
					{
						cout << endl;
						for (int j = 0; j < Nh; j++) { cout << uout[i][j] << "  "; }
					}
					cout << endl;
				}*/
				if (!flag_pass_crc) { n_err_frame_corrected_by_crc++; }
			}
			if (frame % 10 == 0)
			{
				/*printf("\r<%d> SNR= %.2f FER= %d / %d = %f BER= %d / %d = %lf E-Det = %f, E-MRB = %f"  ,
					frame, nSNR, Num_Frame_Error, frame, (double)Num_Frame_Error / frame, Num_Error, (N_Info * frame), (double)Num_Error / N_Info / frame,
					(double)n_det_err/total_bit_cnt, (double)n_det_err_mrb/(total_bit_cnt/(ModType/2)));*/
				int L_sdis_sum = L_LSD;
				sdis_sum = 0;
				for (int i = 0; i < Nv; i++) sdis_sum += sum_sdis_array[i];
				if (polarde_array[0] == "SCL") L_sdis_sum = L;
				// cout << "sdis_sum_avg" << sdis_sum / frame / L_sdis_sum << endl;
				printf("\r<%d> SNR= %.2f FER= %d / %d = %f BER= %d / %d = %lf ,sdis_avg = %f, E-MRB = %f",
					frame, nSNR, Num_Frame_Error, frame, (double)Num_Frame_Error / frame, Num_Error, (N_Info * frame), (double)Num_Error / N_Info / frame, (double)sdis_sum / frame / L_sdis_sum, (double)n_det_err / total_bit_cnt);
				fflush(stdout);
			}


			// output interleaving pattern
			//if (Num_Frame_Error % record_frames == 0 && Num_Frame_Error > 0 && Num_Frame_Error!=n_frame_error_pre)
			//{
			//	n_frame_error_pre = Num_Frame_Error;
			//	cout << endl;
			//	int n_file = Num_Frame_Error / record_frames;
			//	frame_delta = frame - frame_tmp;
			//	cout << "Delta Frames: " << frame_delta << endl;
			//	frame_tmp = frame;
			//	string ofile_intrlv_name = string("intrlv_pattern_") + to_string(n_file) + string("_frame_") + to_string(frame_delta) + string(".txt");
			//	string ofile_spatial_intrlv_name = string("intrlv_pattern_spatial") + to_string(n_file) + string(".txt");
			//	ofstream ofile_intrlv(ofile_intrlv_name);
			//	ofstream ofile_spatial_intrlv(ofile_spatial_intrlv_name);
			//	for (int i = 0; i < Nv; i++)
			//	{
			//		for (int j = 0; j < intrlv_pattern_num; j++)
			//		{

			//			for (int k = 0; k < ModType; k++)
			//				ofile_intrlv << intrlv_pattern_2D_time_partial[i][j][k] << " ";
			//			//cout << endl;
			//		}
			//	}
			//	ofile_intrlv.close();

			//	for (int i = 0; i < Nh; i++)
			//	{
			//		for (int j = 0; j < Nv; j++)
			//			ofile_spatial_intrlv << intrlv_pattern_2D[i][j] << " ";
			//	}
			//	ofile_spatial_intrlv.close();

			//	// re-generate interleaving pattern
			//	generateIntrlvPattern(intrlv_pattern, Nv);
			//	for (int i = 0; i < Nh; i++) { generateIntrlvPattern(intrlv_pattern_2D[i], Nv); }
			//	for (int i = 0; i < Nv; i++) { generateIntrlvPattern(intrlv_pattern_2D_time[i], Nh); }
			//	if (flag_interleave_partial)
			//	{
			//		for (int i = 0; i < Nv; i++)
			//		{
			//			for (int j = 0; j < intrlv_pattern_num; j++)
			//			{
			//				generateIntrlvPattern(intrlv_pattern_2D_time_partial[i][j], ModType);
			//				/*for (int k = 0; k < ModType; k++)
			//					cout << intrlv_pattern_2D_time_partial[i][j][k];
			//				cout << endl;*/
			//			}
			//		}
			//	}
			//	//Num_Frame_Error = 0;
			//}


			//printf("%d %d", Num_Frame_Error, Num_Error);
			//getchar();
			//Decoding_Duration += (T_finish - T_start);
			/*auto dec_end = std::chrono::system_clock::now();
			auto duration_decoding = std::chrono::duration_cast<std::chrono::microseconds>(dec_end - det_end);
			auto duration_total = std::chrono::duration_cast<std::chrono::microseconds>(dec_end - start);
			file_dec << duration_decoding.count() << endl;
			file_total << duration_total.count() << endl;*/
			if (frame % 16 == 0)
			{
				file_err_num << Num_Error << endl;
			}


			// logging crc output
			if ( (frame % 1000) == 0)
			{
				file_crc_check_decoding << endl;
				file_crc_check_decoding << "SNR: " << nSNR << endl;
				file_crc_check_decoding << "FER: " << (double)Num_Frame_Error / frame << endl;
				file_crc_check_decoding << "Total Frames: " << frame << endl;
				file_crc_check_decoding << "Error frames undetected by crc: " << n_err_frame_pass_crc << endl;
				file_crc_check_decoding << "Error frames corrected by crc: " << n_err_frame_corrected_by_crc << endl;
				file_crc_check_decoding << "Error frames failed to be corrected by crc: " << n_err_frame_uncorrected_by_crc << endl;
				file_crc_check_decoding << "Error of replacement codewords: " << n_err_replacement_codeword_crc << endl;
				file_crc_check_decoding << "Error frames that pass parity check: " << n_err_pass_parity_check << endl;
				file_crc_check_decoding << "Error frames with no possible replacement codewords: " << n_err_frame_withno_replacement_codeword << endl;
				file_crc_check_decoding << "Position of correct codeword within each list of vertical codewords" << endl;
				for (int i = 0; i < Nh; i++)
				{
					file_crc_check_decoding << n_err_codeword_v[i] << " ";
				}
				file_crc_check_decoding << endl<< "Num of correct codewords in each position" << endl;
				for (int i = 0; i < Nh; i++)
				{
					file_crc_check_decoding << i+1 << "th codeword: ";
					for (int j = 0; j < L + 1; j++)
						file_crc_check_decoding << n_err_codeword_v_diffpos[i][j] << "  ";
					file_crc_check_decoding << endl;
				}
				file_crc_check_decoding << endl << "Max correct codeword idx when the correct codeword is found" << endl;
				for (int i=0; i<L+1; i++)
					file_crc_check_decoding << max_L_cnt[i] << "  ";
				file_crc_check_decoding << endl;
			}

		}
		//Decoding_Duration = Decoding_Duration * 100 / CLOCKS_PER_SEC / frame;
		//Decoding_Duration = Decoding_Duration / CLOCKS_PER_SEC / Max_frame;

		//printf("%f  %f  %f  %f\n", nSNR, Decoding_Duration, (double)Num_Error / N_Info / Max_frame, (double)Num_Frame_Error / Max_frame);
		cout << endl;
		printf("%f %f %f\n", nSNR, (double)Num_Frame_Error / frame, (double)Num_Error / N_Info / frame);
		cout << endl;
		if (bit_flip_pyndiah)
		{
			for (int i = 0; i < softiter; i++)
			{
				cout << "# Tried Patterns (Horizontal):   " << setw(10) << (float)n_tried_patterns_h_all[i] / (float)(Nv * frame);
				cout << "# Tried Patterns (Vertical):   " << setw(10) << (float)n_tried_patterns_v_all[i] / (float)(Nh * frame);
				cout << endl;
			}
		}
		BER_array[SNR_count] = (double)Num_Error / N_Info / frame;
		FER_array[SNR_count] = (double)Num_Frame_Error / frame;
		// save to .txt file
		string polarde_serial = "";
		for (int i = 0; i < softiter; i++)
			polarde_serial = polarde_serial + polarde_array[i] + string("_");
		string file_name = string("ResFinal\\2D_LSD_D_") + polarde_serial + MIMODET + string("T") + to_string(Nt) + string("R") + to_string(Nr)
			+ string("Q") + to_string(int(pow(2, ModType))) + string("L") + to_string(Nh) + string("K") + to_string(Kh)
			+ string("T") + to_string(softiter) + string("H") + to_string(half_iter) + string("UB") + to_string(usebeta) + string("I") + to_string(flag_interleave)
			+ string(".txt");
		save_array_txt(nSNR, BER_array[SNR_count], FER_array[SNR_count], file_name);
		int L_sdis_sum = L_LSD;
		sdis_sum = 0;
		for (int i = 0; i < Nv; i++) sdis_sum += sum_sdis_array[i];
		if (polarde_array[0] == "SCL") L_sdis_sum = L;
		cout << "sdis_sum_avg" << sdis_sum / frame / L_sdis_sum << endl;
		SNR_count++;
		delete[] n_tried_patterns_h;
		delete[] n_tried_patterns_v;


		// logging crc details
		if (true)
		{
			file_crc_check_decoding << endl;
			file_crc_check_decoding << "FER: " << (double)Num_Frame_Error / frame << endl;
			file_crc_check_decoding << "Total Frames: " << frame << endl;
			file_crc_check_decoding << "Error frames undetected by crc: " << n_err_frame_pass_crc << endl;
			file_crc_check_decoding << "Error frames corrected by crc: " << n_err_frame_corrected_by_crc << endl;
			file_crc_check_decoding << "Error frames failed to be detected by crc: " << n_err_frame_uncorrected_by_crc << endl;
			file_crc_check_decoding << "Error of replacement codewords: " << n_err_replacement_codeword_crc << endl;
		}

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
		delete[] y[i];	for (int i = 0; i < Nh; i++)
		delete[] oriLLR_kbest[i];
	delete[] oriLLR_kbest;
	for (int i = 0; i < Nv; i++)
		delete sub_codeword_list[i];
	delete[] sub_codeword_list;
	for (int i = 0; i < L_WHOLE_CODEWORD; i++)
		delete whole_codeword_list[i];
	delete[] whole_codeword_list;
	for (int i=0; i<Nv; i++)
		delete[]edistance_list[i];
	delete[] edistance_list;
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
	for (int i = 0; i < mimo_blk_height; i++) {
		delete[] x_mod[i];
	}
	delete[] x_mod;
	for (int i = 0; i < Nr; i++) {
		delete[] y_mod[i];
	}
	delete[] y_mod;
	for (int i = 0; i < 2 * mimo_blk_height; i++) {
		delete[] y_mmse[i];
	}
	delete[] y_mmse;
	for (int i = 0; i < 2 * Nr; i++) {
		delete[] y_mod_real[i];
	}
	for (int i = 0; i < mimo_blk_length; i++)
	{
		delete[] y_mod_real_temp_omp[i];
	}
	delete[] y_mod_real_temp_omp;
	for (int i = 0; i < mimo_blk_length; i++)
	{
		delete[] y_mmse_temp_omp[i];
	}
	delete[] y_mmse_temp_omp;
	for (int i = 0; i < mimo_blk_length; i++)
	{
		delete[] ep_extr_mean_omp[i];
	}
	delete[] ep_extr_mean_omp;
	for (int i = 0; i < mimo_blk_length; i++)
	{
		delete[] ep_extr_var_omp[i];
	}
	for (int i = 0; i < Nv; ++i) {
		for (int j = 0; j < n_samples; ++j) {
			delete[] tried_codewords_h[i][j];  // �ͷŵ���ά����������
		}
		delete[] tried_codewords_h[i];  // �ͷŵڶ�ά��ָ������
	}
	delete[] tried_codewords_h;  // �ͷŵ�һά��ָ������
	for (int i = 0; i < Nh; ++i) {
		for (int j = 0; j < n_samples; ++j) {
			delete[] tried_codewords_v[i][j];  // �ͷŵ���ά����������
		}
		delete[] tried_codewords_v[i];  // �ͷŵڶ�ά��ָ������
	}
	delete[] tried_codewords_v;  // �ͷŵ�һά��ָ������
	delete[] ep_extr_var_omp;
	delete[] y_mod_real;
	delete[] x_mod_tmp;
	delete[] y_mod_tmp;
	delete[] y_mod_real_tmp;
	delete[] x_tmp;
	delete[] oriLLR_tmp;
	delete[] sum_sdis_array;
	delete[] Hreal_array;
	delete[] HTH_array;
	delete[] HTY_array;
	delete[] intrlv_pattern_2D_inner;
}








void system2D_LSD_epdet_with_CRC_vertical(double* BER_array, double* FER_array)
{
	std::ofstream file_crc_check_replacement_codeword("crc_replacement_codeword_log.txt", std::ios::app);
	std::ofstream file_crc_check_decoding("crc_check_decoding.txt", std::ios::app);
	std::ofstream file_total("total_timing.txt", std::ios::app);
	std::ofstream file_det("det_timing.txt", std::ios::app);
	std::ofstream file_dec("dec_timing.txt", std::ios::app);
	std::ofstream file_err_num("total_frames.txt", std::ios::app);
	auto det_end = std::chrono::system_clock::now();
	int Nt2 = 2 * Nt;
	int Nr2 = 2 * Nr;

	/*if (!file.is_open()) {
		std::cerr << "Error opening file: " << file_name << std::endl;
		return;
	}*/

	//file << setw(20) << SNR << setw(20) << BER << setw(20) << FER << "\n";

	//freopen("result.txt", "w", stdout);
	printsetting();
	cout << "\nSNR	 " << "FER      BER" << endl;
	int N_Info = Kv * Kh - CRCL, Nmax;
	//Interleaving pattern
	/*std::vector<int> intrlv_pattern(Nv);
	std::vector<int>* intrlv_pattern_2D = new vector<int>[Nh];
	for (int i = 0; i < Nh; i++) intrlv_pattern_2D[i].resize(Nv);*/
	int intrlv_pattern_num = Nh / ModType;
	std::vector<int> intrlv_pattern(Nv);
	std::vector<int>* intrlv_pattern_2D = new vector<int>[Nh];
	std::vector<int>* intrlv_pattern_2D_time = new vector<int>[Nv];
	std::vector<int>* intrlv_pattern_2D_inner = new vector<int>[Nh];
	std::vector<int>** intrlv_pattern_2D_time_partial = new vector<int> *[Nv];
	for (int i = 0; i < Nv; i++) intrlv_pattern_2D_time_partial[i] = new vector<int>[intrlv_pattern_num];
	for (int i = 0; i < Nv; i++)
	{
		for (int j = 0; j < intrlv_pattern_num; j++)
		{
			intrlv_pattern_2D_time_partial[i][j].resize(ModType);
		}
	}
	for (int i = 0; i < Nh; i++) intrlv_pattern_2D[i].resize(Nv);
	for (int i = 0; i < Nh; i++) intrlv_pattern_2D_inner[i].resize(Kv);
	for (int i = 0; i < Nv; i++) intrlv_pattern_2D_time[i].resize(Nh);

	float CodeRate = (float)N_Info / (Nh * Nv);
	Nmax = max(Nh, Nv);
	// num of symbols per vertical codeword
	int Nsym = 0;
	if (MOD_DIM == "Spatial") { Nsym = Nv / ModType; }
	else if (MOD_DIM == "Temporal") { Nsym = Nh / ModType; }
	int symbol_length = ModType;

	// Size of a mimo block is given by height(spatial) * length(temporal)
	int mimo_blk_height = 0;
	int mimo_blk_length = 0;
	if (MOD_DIM == "Spatial")
	{
		mimo_blk_height = Nsym;
		mimo_blk_length = Nh;
	}
	else if (MOD_DIM == "Temporal")
	{
		mimo_blk_height = Nv;
		mimo_blk_length = Nsym;
	}
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

	// Arrays for MIMO
	int** x_mmse_out = new int* [Nv]; // delete!
	for (int i = 0; i < Nv; i++) { x_mmse_out[i] = new int[Nh]; }
	int** x_parity_check = new int* [Nv];
	for (int i = 0; i < Nv; i++) { x_parity_check[i] = new int[Nh]; }

	//std::complex<float>** x_mod = new std::complex<float>* [Nsym]; // modulated symbols (spatial modulation) DELETE!
	//for (int i = 0; i < Nsym; i++) { x_mod[i] = new std::complex<float>[Nh]; }
	//std::complex<float>** x_mod_temporal = new std::complex<float>*[Nv]; // modulated symbols (temporal modulation) DELETE!
	//for (int i = 0; i < Nv; i++) x_mod_temporal[i] = new std::complex<float>[Nsym];
	//std::complex<float>** y_mod = new std::complex<float>* [Nr]; // Output symbols after AWGN channels // delete!
	//for (int i = 0; i < Nr; i++) { y_mod[i] = new std::complex<float>[Nh]; }
	//float** y_mod_real = new float* [2 * Nr];
	//for (int i = 0; i < 2*Nr; i++) y_mod_real[i] = new float[Nh];
	//float** y = new float* [Nv];// Output Signal After AWGN Channel
	//for (int i = 0; i < Nv; i++) { y[i] = new float[Nh]; }
	//float** y_mmse = new float* [2*Nsym];
	//for (int i = 0; i < 2*Nsym; i++) { y_mmse[i] = new float[Nh]; }  // MMSE detector output // delete!
	//float* y_mmse_tmp = new float[2*Nsym];
	//std::complex<float>* x_mod_tmp = new std::complex<float>[Nsym]; // modulated symbols of a single column // delete!
	//std::complex<float>* y_mod_tmp = new std::complex<float>[Nr]; // received symbols of a single column // delete!
	//float* y_mod_real_tmp = new float[2*Nr];
	//// float* y_mod_tmp = new float[Nr];  // delete!


	std::complex<float>** x_mod = new std::complex<float>*[mimo_blk_height]; // modulated symbols (spatial modulation) DELETE!
	for (int i = 0; i < mimo_blk_height; i++) { x_mod[i] = new std::complex<float>[mimo_blk_length]; }
	std::complex<float>** y_mod = new std::complex<float>*[Nr]; // Output symbols after AWGN channels // delete!
	for (int i = 0; i < Nr; i++) { y_mod[i] = new std::complex<float>[mimo_blk_length]; }
	float** y_mod_real = new float* [2 * Nr];
	for (int i = 0; i < 2 * Nr; i++) y_mod_real[i] = new float[mimo_blk_length];
	float** y = new float* [Nv];// Output Signal After AWGN Channel
	for (int i = 0; i < Nv; i++) { y[i] = new float[Nh]; }
	float** y_mmse = new float* [2 * mimo_blk_height];
	for (int i = 0; i < 2 * mimo_blk_height; i++) { y_mmse[i] = new float[mimo_blk_length]; }  // MMSE detector output // delete!
	float* y_mmse_tmp = new float[2 * mimo_blk_height];
	std::complex<float>* x_mod_tmp = new std::complex<float>[mimo_blk_height]; // modulated symbols of a single column // delete!
	std::complex<float>* y_mod_tmp = new std::complex<float>[Nr]; // received symbols of a single column // delete!
	float* y_mod_real_tmp = new float[2 * Nr];
	// float* y_mod_tmp = new float[Nr];  // delete!
	float** y_mod_real_temp_omp = new float* [mimo_blk_length]; // new delete!
	for (int i = 0; i < mimo_blk_length; i++) { y_mod_real_temp_omp[i] = new float[2 * Nr]; }
	float** y_mmse_temp_omp = new float* [mimo_blk_length]; // new delete!
	for (int i = 0; i < mimo_blk_length; i++) { y_mmse_temp_omp[i] = new float[2 * mimo_blk_height]; }
	// ep
	float** ep_extr_mean_omp = new float* [mimo_blk_length];
	for (int i = 0; i < mimo_blk_length; i++) { ep_extr_mean_omp[i] = new float[2 * mimo_blk_height]; }
	float** ep_extr_var_omp = new float* [mimo_blk_length];
	for (int i = 0; i < mimo_blk_length; i++) { ep_extr_var_omp[i] = new float[2 * mimo_blk_height]; }
	int* x_tmp = new int[Nv]; // delete!

	// paths selected by k-best detection


	// CRC bits
	unsigned char* u_crc = new unsigned char[Kv * Kh];

	float** upLLR = new float* [Nv];
	for (int i = 0; i < Nv; i++) { upLLR[i] = new _MM_ALIGN16 float[Nh]; }
	//float** oriLLRv = new float* [Nh];
	//for (int i = 0; i < Nh; i++) { oriLLRv[i] = new _MM_ALIGN16 float[2 * Nv + 2]; }
	float** oriLLRh = new float* [Nv];
	for (int i = 0; i < Nv; i++) { oriLLRh[i] = new _MM_ALIGN16 float[Nh]; }
	int len_ori_llr_tmp = 0;
	if (MOD_DIM == "Spatial") { len_ori_llr_tmp = Nv; }
	else if (MOD_DIM == "Temporal") { len_ori_llr_tmp = ModType * Nt; }
	float* oriLLR_tmp = new _MM_ALIGN16 float[len_ori_llr_tmp];  // delete!
	//float** oriLLR_kbest = new float* [Nh];
	float** oriLLR_kbest = new float* [mimo_blk_length];
	for (int i = 0; i < Nh; i++) { oriLLR_kbest[i] = new float[len_ori_llr_tmp]; }
	double* sum_sdis_array = new double[Nv];
	for (int i = 0; i < Nv; i++) sum_sdis_array[i] = 0;

	//*****************************************************************
	// EP detector
	float* Hreal_array = new float[(2 * Nr) * (2 * Nt)]; // delete
	float* HTH_array = new float[(2 * Nt) * (2 * Nt)]; // delete
	float* HTY_array = new float[2 * Nt]; // delete 





	int** GGh = new int* [Nh];
	int** GGv = new int* [Nv];


	MatrixXf Gh = Fn(Nh);  // Generate Matrix of horizontal coding
	MatrixXf Gv = Fn(Nv);  // Generate Matrix of vertical coding
	MatrixXcf H(Nr, Nt); // Complex channel Matrix H: Nr x Nt
	//Matrix<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic> H = MatrixXcf::Random(Nr, Nt); // Real domain channel matrix
	//// //MatrixX Definitions
	MatrixXf H_real(2 * Nr, 2 * Nt); // Real domain channel matrix
	MatrixXf received_signal(2 * Nr, 1);
	MatrixXf HTH(2 * Nt, 2 * Nt);
	MatrixXf Identity = MatrixXf::Identity(2 * Nt, 2 * Nt);
	MatrixXf output_signal(2 * Nt, 1);
	MatrixXf mmse_matrix(2 * Nt, 2 * Nr);
	MatrixXf HTHIinv(2 * Nt, 2 * Nt);

	// for KBest
	MatrixXf R(2 * Nt, 2 * Nt);
	MatrixXf Q(2 * Nr, 2 * Nt);
	MatrixXf z(2 * Nt, 1);

	// for CRC
	int** whole_codeword_list = new int* [L_WHOLE_CODEWORD];
	for (int i = 0; i < L_WHOLE_CODEWORD; i++) { whole_codeword_list[i] = new int[Nh * Nv]; }
	int L_max = max(L_LSD, L);

	// stores the decoding result of Nv horizontal codewords
	int** sub_codeword_list = new int* [Nv];
	for (int i = 0; i < Nv; i++) { sub_codeword_list[i] = new int[L_max * Nh]; }
	float** edistance_list = new float* [Nv];
	for (int i = 0; i < Nv; i++) { edistance_list[i] = new float[L_max]; }
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> H_real = MatrixXf::Random(2 * Nr, 2 * Nt); // Real domain channel matrix
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> received_signal = MatrixXf::Random(2 * Nr,1);
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> HTH = MatrixXf::Random(2 * Nt, 2 * Nt);
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> Identity = MatrixXf::Random(2 * Nt, 2 * Nt);
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> output_signal = MatrixXf::Random(2 * Nt, 1);
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> mmse_matrix = MatrixXf::Random(2 * Nt, 2 * Nr);


	// Matrix Definitions
	//Eigen::Matrix<float, 2 * Nr, 2 * Nt> H_real; // ʵ���ŵ�����
	//Eigen::Matrix<float, 2 * Nr, 1> received_signal; // �����źž���
	//Eigen::Matrix<float, 2 * Nt, 2 * Nt> HTH; // HTH����
	//Eigen::Matrix<float, 2 * Nt, 2 * Nt> Identity = Eigen::Matrix<float, 2 * Nt, 2 * Nt>::Identity(); // ��λ����
	//Eigen::Matrix<float, 2 * Nt, 1> output_signal; // ����źž���
	//Eigen::Matrix<float, 2 * Nt, 2 * Nr> mmse_matrix; // MMSE����


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
	// srand((unsigned int)time(0));
	std::random_device rd;              // ��ȡ���������
	std::mt19937 gen(rd());             // ʹ�����ӳ�ʼ��mt19937������
	std::uniform_int_distribution<> distrib(1, 1000000);
	int Num_Error, Num_Frame_Error;
	float Decoding_Duration;

	int Num_Error_new, Num_Frame_Error_new;
	float Decoding_Duration_new;

	clock_t T_start, T_finish;
	int SNR_count = 0;
	int anssc, anssd;
	for (double nSNR = Min_SNR; nSNR <= Max_SNR; nSNR += SNR_step)
	{
		// crc check logs
		int n_err_frame_pass_crc = 0;  // err frames not detected by crc
		int n_err_frame_corrected_by_crc = 0; // err frames detected and corrected by crc
		int n_err_frame_uncorrected_by_crc = 0; // err frames detected but not corrected by crc
		int n_err_replacement_codeword_crc = 0;  // err of replacement codewords
		int n_err_pass_parity_check = 0; // n of err that passes parity check
		int idx_replacement_codeword = 0; // idx of replacement codeword that passes crc check

		for (int i = 0; i < Nv; i++) sum_sdis_array[i] = 0;
		double sdis_sum = 0;
		// TEST: det err
		int n_det_err = 0; // n error bits at detector output
		int n_det_err_mrb = 0; // error bits at detector output(MRBs)
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

		for (int i = 0; i < Nh; i++)
		{
			cout << A_Ach[i] << ",";
		}
		cout << endl;
		// test: Horizontal frozen bits

		/*cout << "Frozen Bits" <<endl;
		for (int i = 0; i < Nh - Kh; i++)
		{
			cout << setw(5) << Acv[i];
		}
		cout << endl;*/
		// test: Horizontal frozen bits
		// Frame Loop
		int passnum;
		int frame = 0;

		// generate interleaving patterns
		if (flag_interleave)
		{
			string ofile_intrlv_name = string("intrlv_pattern_71_frame_49263.txt");
			string ofile_spatial_intrlv_name = string("intrlv_pattern_spatial71.txt");
			ifstream ofile_intrlv(ofile_intrlv_name);
			ifstream ofile_spatial_intrlv(ofile_spatial_intrlv_name);
			int file_in = 1;
			for (int i = 0; i < Nv; i++)
			{
				for (int j = 0; j < intrlv_pattern_num; j++)
				{

					for (int k = 0; k < ModType; k++)
						// ofile_intrlv >> file_in;
						ofile_intrlv >> intrlv_pattern_2D_time_partial[i][j][k];
					//cout << endl;
				}
			}
			ofile_intrlv.close();

			for (int i = 0; i < Nh; i++)
			{
				for (int j = 0; j < Nv; j++)
					ofile_spatial_intrlv >> intrlv_pattern_2D[i][j];
			}
			ofile_spatial_intrlv.close();
			/*for (int i = 0; i < Nh; i++)
			{
				for (int j = 0; j < Nv; j++)
					cout << intrlv_pattern_2D[i][j] << endl;
			}*/
			generateIntrlvPattern(intrlv_pattern, Nv);
			for (int i = 0; i < Nh; i++) { generateIntrlvPattern(intrlv_pattern_2D[i], Nv); }
			for (int i = 0; i < Nv; i++) { generateIntrlvPattern(intrlv_pattern_2D_time[i], Nh); }
			if (flag_interleave_partial)
			{
				for (int i = 0; i < Nv; i++)
				{
					for (int j = 0; j < intrlv_pattern_num; j++)
					{
						generateIntrlvPattern(intrlv_pattern_2D_time_partial[i][j], ModType);
						/*for (int k = 0; k < ModType; k++)
							cout << intrlv_pattern_2D_time_partial[i][j][k];
						cout << endl;*/
					}
				}
			}
			for (int i = 0; i < Nh; i++) { generateIntrlvPattern(intrlv_pattern_2D_inner[i], Kv); }
		}

		//ofstream ofile_intrlv("intrlv_pattern.txt");
		//for (int i=0; i<Nv; i++)
		//{
		//	for (int j = 0; j < intrlv_pattern_num; j+for+)
		//	{
		//		
		//		for (int k = 0; k < ModType; k++)
		//			ofile_intrlv << intrlv_pattern_2D_time_partial[i][j][k] << "   ";
		//		//cout << endl;
		//	}
		//}
		//ofile_intrlv.close();
		int frame_tmp = 0;
		int frame_delta = 0;
		int n_frame_error_pre = 0;
		while (Num_Frame_Error < errframe)
		{
			/*test: processing time*/
			//auto start = std::chrono::system_clock::now();
			/*if (frame == 10) {
				cout << "AAAAAAAAAA" << endl;
			}*/
			frame++;
			bool flag_pass_crc = 0;
			bool flag_pass_crc_undetected = 0;
			bool flag_pass_crc_corrected = 0;
			bool flag_found_replacement_codeword = 0;

			generateIntrlvPattern(intrlv_pattern, Nv);
			for (int i = 0; i < Nh; i++) { generateIntrlvPattern(intrlv_pattern_2D[i], Nv); }
			for (int i = 0; i < Nv; i++) { generateIntrlvPattern(intrlv_pattern_2D_time[i], Nh); }
			if (flag_interleave_partial)
			{
				for (int i = 0; i < Nv; i++)
				{
					for (int j = 0; j < intrlv_pattern_num; j++)
					{
						generateIntrlvPattern(intrlv_pattern_2D_time_partial[i][j], ModType);
						/*for (int k = 0; k < ModType; k++)
							cout << intrlv_pattern_2D_time_partial[i][j][k];
						cout << endl;*/
					}
				}
			}
			for (int i = 0; i < Nh; i++) { generateIntrlvPattern(intrlv_pattern_2D_inner[i], Kv); }
			// generate interleaving pattern
			/*if (flag_interleave)
			{
				generateIntrlvPattern(intrlv_pattern, Nv);
				for (int i = 0; i < Nh; i++) { generateIntrlvPattern(intrlv_pattern_2D[i], Nv); }
				for (int i = 0; i < Nv; i++) { generateIntrlvPattern(intrlv_pattern_2D_time[i], Nh); }
				if (flag_interleave_partial)
				{
					for (int i = 0; i < Nv; i++)
					{
						for (int j = 0; j < intrlv_pattern_num; j++)
						{
							generateIntrlvPattern(intrlv_pattern_2D_time_partial[i][j], ModType);
						}
					}
				}
			}*/

			//cout << frame << endl;
			int u_input[16] = { 0,     1     ,0,     1,
									0,     1,     0,     0,
									0 ,    1,     0,     0,
									1 ,    0,     0,     0 };
			for (int i = 0; i < N_Info; i++)
			{
				a_data[i] = (unsigned char)(distrib(gen) % 2);
				//		a_data[i] = u_input[i];
			}

			//cout << endl;
			if (CRCL) tx_append_crc(a_data, N_Info, CRCL);
			for (int i = 0; i < Nv; i++)
				for (int j = 0; j < Nh; j++)
					u[i][j] = 0;
			//omp_set_num_threads(4);
//#pragma omp parallel for
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
			// test:original info
			//cout << "original info u:" << endl;
			////for (int j = 0; j < Nh; j++) { cout << u[Av[0]][j] << " "; }
			//for (int i = 0; i < Nv; i++)
			//{
			//	cout << endl;
			//	for (int j = 0; j < Nh; j++)
			//		cout << u[i][j];
			//}
			//cout << endl;
			// Polar Encoding
//#pragma omp parallel for

			// test: original info
			/*cout << "original info u:" << endl;
			display_array(u, Nv, Nh);*/
			// test: original info

			for (int i = 0; i < Nv; i++)
				PolarEncode_xor(midx[i], u[i], Nh);
			// test:horizontal encoding
			//cout << "u after horizontal encoding:" << endl;
			////for (int j = 0; j < Nh; j++) { cout << midx[Av[0]][j] << " "; }
			//for (int i = 0; i < Nv; i++)
			//{
			//	cout << endl;
			//	for (int j = 0; j < Nh; j++)
			//		cout << midx[i][j];
			//}
			//cout << endl;
			//cout << endl;
			// test:horizontal encoding
			// 
			// column-wise interleaving

			// test: row encoding
			/*cout << "u after horizontal encoding:" << endl;
			display_array(midx, Nv, Nh);*/
			// test: row encoding
			if (flag_interleave_inner)
			{
				int* tmp_to_interleave = new int[Nv];
				for (int i = 0; i < Nh; i++)
				{
					for (int j = 0; j < Nv; j++) { tmp_to_interleave[j] = midx[j][i]; }
					randIntrlvPartial(intrlv_pattern_2D_inner[i], tmp_to_interleave, Av, Nv, Kv);
					for (int j = 0; j < Nv; j++) { midx[j][i] = tmp_to_interleave[j]; }
				}
				delete[] tmp_to_interleave;
			}

			// test: innner interleaving
			/*cout << "u after inner interleaving:" << endl;
			display_array(midx, Nv, Nh);*/
			// test: innner interleaving
//#pragma omp parallel for
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

			// show x[0] after horizontal decoding
			/*cout << "x[0]" << endl;
			for (int i = 0; i < Nh; i++) cout << x[0][i] << "  ";
			cout << endl;*/
			// show x[0] after horizontal decoding
			// test: vertical encoding
			/*cout << "u after vertical encoding:" << endl;
			display_array(x, Nv, Nh);*/
			// test: vertical encoding
			/*for (int i = 0; i < 2 * Nr; i++)
				for (int j = 0; j < 2 * Nt; j++)
					H_real(i, j) = 0;*/
					/*MatrixXf H_real_test(2 * Nr, 2 * Nt);
					MatrixXf M = H_real_test.transpose();
					MatrixXf N = H_real_test;
					HTHIinv = M * N;*/
					//HTHIinv = M * H_real;
					// test:vertical encoding
					//cout << "u after vertical encoding:" << endl;
					////for (int j = 0; j < Nh; j++) { cout << x[Av[0]][j] << " "; }
					//for (int i = 0; i < Nv; i++)
					//{
					//	cout << endl;
					//	for (int j = 0; j < Nh; j++)
					//		cout << x[i][j];
					//}
					//cout << endl;
					//cout << endl;
					// test:vertical encoding
			if (atten == 2)
			{
				// modulation
				// #pragma omp parallel for
				// Interleaving
				if (flag_interleave)
				{
					if (flag_interleave_all)  // interleave temporally
					{
						for (int i = 0; i < Nv; i++)
						{
							randIntrlv(intrlv_pattern_2D_time[i], x[i], Nh);
						}
					}
					//// test: to be interleaved
					//cout << "x before interleaving:" << endl;
					//// print x line by line, add spaces between each 4 terms
					//for (int i = 0; i < Nv; i++)
					//{
					//	for (int j = 0; j < Nh; j++)
					//	{
					//		if (!(j % 4)) cout << " ";
					//		cout << x[i][j] << " ";
					//	}
					//	cout << endl;
					//}

					if (flag_interleave_partial)
					{
						for (int i = 0; i < Nv; i++)
						{
							for (int j = 0; j < intrlv_pattern_num; j++)
							{
								randIntrlv(intrlv_pattern_2D_time_partial[i][j], x[i] + j * ModType, ModType);
							}
						}
					}

					//// test: to be interleaved
					//cout << "x after interleaving:" << endl;
					//// print x line by line, add spaces between each 4 terms
					//for (int i = 0; i < Nv; i++)
					//{
					//	for (int j = 0; j < Nh; j++)
					//	{
					//		if (!(j % 4)) cout << " ";
					//		cout << x[i][j] << " ";
					//	}
					//	cout << endl;
					//}
					if ((!flag_interleave_block) && flag_interleave_spatial)
					{
						for (int j = 0; j < Nh; j++)
						{
							for (int i = 0; i < Nv; i++) { x_tmp[i] = x[i][j]; }
							randIntrlv(intrlv_pattern_2D[j], x_tmp, Nv);
							for (int i = 0; i < Nv; i++) { x[i][j] = x_tmp[i]; }
						}
					}


					if (flag_interleave_block)
					{
						randIntrlvBlock(x, Nh, Nv, ModType, interleave_rows);
					}
				}
				if (MOD_DIM == "Spatial")
				{
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
						/*for (int i = 0; i < mimo_blk_height; i++)
						{
							cout << x_mod_tmp[i] << endl;
						}
						cout << "111111" << endl;*/
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
					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < 2 * Nr; j++) { y_mod_real_temp_omp[i][j] = y_mod_real[j][i]; }
					}
					// for KBest Detection
					Eigen::HouseholderQR<Eigen::MatrixXf> qr_Hreal(H_real);
					Q = qr_Hreal.householderQ();
					R = qr_Hreal.matrixQR().triangularView<Eigen::Upper>();
					//omp_set_num_threads(8);

#pragma omp parallel for private(z)
					for (int j = 0; j < Nh; j++)
					{
						MatrixXf received_signal_omp(2 * Nr, 1);
						int k_kbest_extend = pow(2, ModType / 2) * K_KBEST;
						float** paths_kbest = new float* [k_kbest_extend];
						for (int i = 0; i < k_kbest_extend; i++) paths_kbest[i] = new float[2 * Nt];
						float* ED = new float[K_KBEST];

						if (MIMODET == "MMSE") { MMSEdetection(H_real, Identity, received_signal, mmse_matrix, output_signal, HTH, y_mod_real_temp_omp[j], y_mmse_temp_omp[j], sigma_sq); }
						// ִ��һЩ����

						for (int i = 0; i < 2 * Nr; i++) { received_signal_omp(i, 0) = y_mod_real_temp_omp[j][i]; }
						z = Q.transpose() * received_signal_omp;
						if (MIMODET == "KBEST")
						{
							KBestDetection(R, z, H_real, y_mod_real_temp_omp[j], ED, y_mmse_temp_omp[j], paths_kbest, K_KBEST, ModType, Nt, Nr);
							llr_kbest_bit(Nt, y_mmse_temp_omp[j], oriLLR_kbest[j], sigma_sq, ModType, K_KBEST, paths_kbest, ED);
						}
						for (int i_path = 0; i_path < k_kbest_extend; i_path++) { delete[] paths_kbest[i_path]; }

						delete[] paths_kbest;
						delete[] ED;
					}

					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < 2 * Nt; j++) { y_mmse[j][i] = y_mmse_temp_omp[i][j]; }
					}
					// sym llr to bitwise llr
					// Detection is done for every column, that is, for every Nv bits
					// Interleaving is not currently considered
					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < 2 * Nt; j++) { y_mmse_tmp[j] = y_mmse[j][i]; }
						if (MIMODET == "MMSE") { llr_mmse_bit(Nt, y_mmse_tmp, oriLLR_tmp, sigma_sq, ModType); }
						if (MIMODET == "KBEST") { for (int j = 0; j < len_ori_llr_tmp; j++) { oriLLR_tmp[j] = oriLLR_kbest[i][j]; } }
						if (flag_interleave)
						{
							randDeintrlv(intrlv_pattern, oriLLR_tmp, Nv);
						}
						bool flagIsMax = false;
						for (int j = 0; j < Nv; j++)
						{
							if (!(j % 2))
							{
								if (abs(x_mmse_out[j][i]) > abs(x_mmse_out[j][i + 1])) { flagIsMax = true; }
								else { flagIsMax = false; }
							}
							else
							{
								if (abs(x_mmse_out[j][i]) > abs(x_mmse_out[j][i - 1])) { flagIsMax = true; }
								else { flagIsMax = false; }
							}
							oriLLRh[j][i] = oriLLR_tmp[j];
							upLLR[j][i] = oriLLR_tmp[j];
							// TEST�� det error
							x_mmse_out[j][i] = 0;
							if (oriLLR_tmp[j] < 0) { x_mmse_out[j][i] = 1; }
							if (x_mmse_out[j][i] != x[j][i]) {
								n_det_err++;
								// error at max-reliable bit. Warning: Only fit for 16QAM modulation
								if (flagIsMax) { n_det_err_mrb++; }
							}
							total_bit_cnt++;
						}
					}
					// test : detection TIME
					/*det_end = std::chrono::system_clock::now();
					auto duration_detection = std::chrono::duration_cast<std::chrono::microseconds>(det_end - start);
					file_det << duration_detection.count() << endl;*/
				}
				// Temporal modulation
				else if (MOD_DIM == "Temporal")
				{
					//modulation
					//cout << "START TEST" << endl;
					for (int i = 0; i < Nv; i++) // every row
					{
						modulation(x[i], x_mod[i], ModType, mimo_blk_length);
					}
					/*for (int i = 0; i < mimo_blk_height; i++)
					{
						for (int j = 0; j < mimo_blk_length; j++)
							cout << x_mod[i][j] << endl;
					}
					cout << "111111" << endl;*/
					//AWGN channel (received y_mod)
					// test: modulated symbols
					/*cout << endl;
					cout << "modulated symbols" << endl;
					for (int i = 0; i < Nt; i++) cout << x_mod[i][0] << endl;
					cout << "modulated symbols" << endl;*/
					// test: modulated symbols
					H = generateChannelMatrix(Nr, Nt);
					for (int j = 0; j < Nsym; j++)
					{
						// Extracting a single column
						for (int i = 0; i < Nv; i++) { x_mod_tmp[i] = x_mod[i][j]; }
						AWGN(Nr, Nt, x_mod_tmp, y_mod_tmp, H, sigma_sq);
						for (int i = 0; i < Nr; i++) { y_mod[i][j] = y_mod_tmp[i]; }
					}
					convertHToReal(H, H_real);

					for (int j = 0; j < Nsym; j++)
					{
						for (int i = 0; i < Nr; i++) { y_mod_tmp[i] = y_mod[i][j]; }
						convertYToReal(Nr, y_mod_tmp, y_mod_real_tmp);
						for (int i = 0; i < 2 * Nr; i++) { y_mod_real[i][j] = y_mod_real_tmp[i]; }
					}
					//MMSE detection done in paralell
					for (int i = 0; i < Nsym; i++)
					{
						for (int j = 0; j < 2 * Nr; j++) { y_mod_real_temp_omp[i][j] = y_mod_real[j][i]; }
					}
					Eigen::HouseholderQR<Eigen::MatrixXf> qr_Hreal(H_real);
					Q = qr_Hreal.householderQ();
					R = qr_Hreal.matrixQR().triangularView<Eigen::Upper>();
					// for EP detection
					// copy H_real
					for (int i = 0; i < Nr2; i++)
					{
						for (int j = 0; j < Nt2; j++)
						{
							Hreal_array[i * Nt2 + j] = H_real(i, j);
						}
					}
					// calculate HTH
					cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, Nt2, Nt2, Nr2, 1.0, Hreal_array, Nt2, Hreal_array, Nt2, 0.0, HTH_array, Nt2);
					//cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, Nt2, 1, Nr2, 1.0, ch_mtx_struct.H_real[cnt], Nt2, data_struct.rx_buffer_real[cnt], 1, 0.0, ch_mtx_struct.hty, 1);

					//omp_set_num_threads(8);
					//cout << 1 << endl;
					//omp_set_num_threads(4);
					/*auto start = std::chrono::system_clock::now();*/
					/*MatrixXf H_real_trans(2*Nt, 2*Nr);
					for (int i = 0; i < 2 * Nt; i++)
						for (int j = 0; j < 2 * Nr; j++)
							H_real_trans(i, j) = H_real(j, i);
					HTHIinv = (H_real_trans * H_real + Identity).inverse();*/
# pragma omp parallel for private(z)
					for (int j = 0; j < Nsym; j++)

						//for (int x=0; x<Nsym; x++)
					{
						// int j = x % Nsym;
						MatrixXf received_signal_omp(2 * Nr, 1);
						int k_kbest_extend = pow(2, ModType / 2) * K_KBEST;
						float** paths_kbest = new float* [k_kbest_extend];
						for (int i = 0; i < k_kbest_extend; i++) { paths_kbest[i] = new float[2 * Nt]; }
						float* ED = new float[K_KBEST];
						float* HTY_array_tmp = new float[2 * Nt];
						// Declared arrays for EP
						int flag = 1;
						int temp = 0;
						//float* extr_mean = new float[2 * Nr];
						//float* extr_var = new float[2 * Nr];
						double* vec_alpha_in = new double[2 * Nr];
						double* vec_gamma_in = new double[2 * Nr];
						double* p_symbol_rearrange = new double[Nt2 * pow(2, ModType / 2)];
						for (int i = 0; i < Nt2 * pow(2, ModType / 2); i++) { p_symbol_rearrange[i] = 1 / (pow(2, ModType / 2)); }
						for (int i = 0; i < Nr2; i++)
						{
							//extr_mean[i] = 0;
							//extr_var[i] = 0;
							vec_alpha_in[i] = 2;
							vec_gamma_in[i] = 0;
						}

						if (MIMODET == "MMSE")
						{
							MMSEdetection(H_real, Identity, received_signal, mmse_matrix, output_signal, HTH, y_mod_real_temp_omp[j], y_mmse_temp_omp[j], sigma_sq);
							//MMSEdetection_noInverse(H_real, Identity, received_signal, mmse_matrix, output_signal, HTH, HTHIinv, y_mod_real_temp_omp[j], y_mmse_temp_omp[j], sigma_sq);
						}
						if (MIMODET == "KBEST")
						{
							for (int i = 0; i < 2 * Nr; i++) { received_signal_omp(i, 0) = y_mod_real_temp_omp[j][i]; }
							z = Q.transpose() * received_signal_omp;
							KBestDetection(R, z, H_real, y_mod_real_temp_omp[j], ED, y_mmse_temp_omp[j], paths_kbest, K_KBEST, ModType, Nt, Nr);
							llr_kbest_bit(Nt, y_mmse_temp_omp[j], oriLLR_kbest[j], sigma_sq, ModType, K_KBEST, paths_kbest, ED);
						}
						if (MIMODET == "EP")
						{
							// test: save H_real, recv_signal
							/*save_mat_txt(y_mod_real_temp_omp[j], Nr2, "RxSig0128.txt");
							save_mat_txt(H_real, "Hreal0128.txt");*/
							// test: save H_real, recv_signal
							cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, Nt2, 1, Nr2, 1.0, Hreal_array, Nt2, y_mod_real_temp_omp[j], 1, 0.0, HTY_array_tmp, 1);
							//test:HTH_array
							/*for (int i = 0; i < Nt2; i++)
								for (int s = 0; s < Nt2; s++)
								{
									if (s == 0) std::cout << endl;
									std::cout << setw(10) << HTH_array[i * Nt2 + s];
								}
							cout << endl;*/
							// TEST:HTY
							/*for (int i = 0; i < Nr2; i++)
								received_signal_omp(i, 0) = y_mod_real_temp_omp[j][i];
							MatrixXf HTransY = H_real.transpose() * received_signal_omp;
							cout << "Eigen" << endl;
							for (int i = 0; i < Nt2; i++) cout << setw(10)<<HTransY(i, 0);
							cout << "mkl" << endl;
							for (int i = 0; i < Nt2; i++) cout << HTY_array[i];*/
							// TEST:HTY
							EPD_detection_float(flag, Nt, Nr, iter_MIMO, ModType, 1.0, sigma_sq, HTH_array, HTY_array_tmp, temp, ep_extr_mean_omp[j], ep_extr_var_omp[j], vec_alpha_in, vec_gamma_in, p_symbol_rearrange);
						}
						//test EP detection
						/*cout << "TEST:EP DETECTION" << endl;
						for (int i = 0; i < Nr2; i++) std::cout << extr_mean[i]<< endl;
						cout << "TEST:EP DETECTION" << endl;*/
						//test EP detection
						for (int i_path = 0; i_path < k_kbest_extend; i_path++) { delete[] paths_kbest[i_path]; }
						delete[] paths_kbest;
						delete[] ED;
						//delete[] extr_mean;
						//delete[] extr_var;
						delete[] vec_alpha_in;
						delete[] vec_gamma_in;
						delete[] p_symbol_rearrange;
						delete[] HTY_array_tmp;
					}
					// test: detection time
					/*auto detection_over = std::chrono::system_clock::now();
					auto duration_detection = std::chrono::duration_cast<std::chrono::microseconds>(detection_over - start);
					cout << endl;
					cout <<"Detection duration" << duration_detection.count() << "us" << endl;*/
					// test: detection time
					for (int i = 0; i < Nsym; i++)
					{
						for (int j = 0; j < 2 * Nt; j++) { y_mmse[j][i] = y_mmse_temp_omp[i][j]; }
					}
					//To bit domain output
					for (int j = 0; j < Nsym; j++)
					{
						for (int i = 0; i < 2 * Nt; i++) { y_mmse_tmp[i] = y_mmse[i][j]; }
						// With Spatial Modulation, Nt = Nv
						if (MIMODET == "MMSE") { llr_mmse_bit(Nt, y_mmse_tmp, oriLLR_tmp, sigma_sq, ModType); }
						if (MIMODET == "EP") { llr_ep_bit(Nt, ep_extr_mean_omp[j], ep_extr_var_omp[j], oriLLR_tmp, ModType); }
						if (MIMODET == "KBEST") { for (int i = 0; i < len_ori_llr_tmp; i++) { oriLLR_tmp[i] = oriLLR_kbest[j][i]; } }
						int i_bit = 0;
						int bit_col = 0;
						for (int j_inner = 0; j_inner < Nv; j_inner++)
						{
							for (int i_bit_col = 0; i_bit_col < ModType; i_bit_col++)
							{
								i_bit = ModType * j_inner + i_bit_col;
								bit_col = j * ModType + i_bit_col;
								oriLLRh[j_inner][bit_col] = oriLLR_tmp[i_bit];
								upLLR[j_inner][bit_col] = oriLLR_tmp[i_bit];
							}

						}
					}
					// test: LLR transfer time
					/*auto llr_to_bit_over = std::chrono::system_clock::now();
					auto duration_llr = std::chrono::duration_cast<std::chrono::microseconds>(llr_to_bit_over - detection_over);
					cout << endl;
					cout<< "LLR transfer duration" << duration_llr.count() << "us" << endl;*/
					// test: LLR transfer time
					if (flag_interleave)
					{
						if ((!flag_interleave_block) && flag_interleave_spatial)
						{
							for (int j = 0; j < Nh; j++)
							{
								for (int i = 0; i < Nv; i++) { oriLLR_tmp[i] = oriLLRh[i][j]; }
								randDeintrlv(intrlv_pattern_2D[j], oriLLR_tmp, Nv);
								for (int i = 0; i < Nv; i++)
								{
									oriLLRh[i][j] = oriLLR_tmp[i];
									upLLR[i][j] = oriLLR_tmp[i];
								}

							}
						}

						if (flag_interleave_block)
						{
							randDeIntrlvBlock(oriLLRh, Nh, Nv, ModType, interleave_rows);
							randDeIntrlvBlock(upLLR, Nh, Nv, ModType, interleave_rows);
						}

						if (flag_interleave_all) // Deinterleave in time domain
						{
							for (int j = 0; j < Nv; j++)
							{
								randDeintrlv(intrlv_pattern_2D_time[j], oriLLRh[j], Nh);
								randDeintrlv(intrlv_pattern_2D_time[j], upLLR[j], Nh);
							}
						}

						if (flag_interleave_partial)
						{
							for (int j = 0; j < Nv; j++)
							{
								for (int k = 0; k < intrlv_pattern_num; k++)
								{
									randDeintrlv(intrlv_pattern_2D_time_partial[j][k], oriLLRh[j] + k * ModType, ModType);
									randDeintrlv(intrlv_pattern_2D_time_partial[j][k], upLLR[j] + k * ModType, ModType);
								}
							}
						}

					}

					for (int i = 0; i < Nv; i++)
					{
						for (int j = 0; j < Nh; j++)
						{
							total_bit_cnt++;
							int bit_dec = (oriLLRh[i][j] < 0);
							if (bit_dec != x[i][j]) n_det_err++;
						}
					}

				}
				// test: detection errors

				// test: detection errors
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
			for (int iter = 1; iter <= softiter; iter++)
			{
				//vertical decoding
		   // test: decoding time
			//auto decoding_start = std::chrono::system_clock::now();

			// test: decoding time
#pragma omp parallel for
				for (int decodingi = 0; decodingi < Nh; decodingi++) // 291.70kb
				{

					//for SCL Decoding
					int Counth = 0, Countinfoh = 0;
					float Qsumunh = 0;
					int** sum = new int* [L];
					//for (int i = 0; i < L; i++)  sum[i] = new  int[Nmax];
					int** u1 = new int* [L];
					for (int i = 0; i < L; i++)  u1[i] = new _MM_ALIGN16 int[Nmax];
					for (int i = 0; i < L; i++)  sum[i] = new _MM_ALIGN16 int[Nmax];
					//for (int i = 0; i < L; i++)  u1[i] = new  int[Nmax];

					float* LLR_in = new float[L];
					float* W = new float[2 * L];
					int* SCL_index = new int[L];
					int* better = new int[L];
					int* worse = new int[L];
					int* indexpm = new int[L];
					float* PM = new float[L];
					float** LLR = new float* [L];
					for (int i = 0; i < L; i++)  LLR[i] = new _MM_ALIGN16 float[2 * Nmax + 2];
					// for (int i = 0; i < L; i++)  LLR[i] = new float[2 * Nmax + 2];
					////////////////////////////////
					for (int j = 0; j < Nv; j++)
						for (int k = 0; k < L; k++)
							LLR[k][j] = upLLR[j][decodingi];
					//cout << LLR[0][j] << " ";

				//cout << endl;
					int index = 0;
					/*			if (L == 1)
								{
									//T_start = clock();
									Count = 0;
									decode(*LLR, *sum, (int)(log(Nv) / log(2)), A_Acv, flag, Nv, Kv);
									//PolarEncode_xor(*u1, *sum, N);
									//T_finish = clock();
								}
								else
								{*/
					for (int k = 0; k < L; k++) { PM[k] = 0; indexpm[k] = k; }
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
					// for parity check
					if (PARITY_CHECK_ENABLE) {
						for (int row = 0; row < Nv; row++) { x_parity_check[row][decodingi] = sum[indexpm[0]][row]; }
					}
					//}
					//calc distance s()
					if (adaptive_beta)
					{
						if (softmethod == "Pyndiah") Pyndiahv_adaptivebeta(LLR, sum, oriLLRh, indexpm, iter * 2 - 2, decodingi, upLLR);
						if (softmethod == "MITSO") MITSOv(LLR, PM, sum, oriLLRh, indexpm, iter * 2 - 2, decodingi, upLLR, Qsumunh);
					}
					// iter * 2 - 2
					else
					{
						//cout << iter << endl;
						if (softmethod == "Pyndiah") Pyndiahv(LLR, sum, oriLLRh, indexpm, iter, decodingi, upLLR);
						if (softmethod == "MITSO") MITSOv(LLR, PM, sum, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunh);
					}
					// test: oriLLR and upLLR
					/*cout << endl;
					for (int i = 0; i < Nv; i++)
						for (int j = 0; j < Nh; j++)
							cout << oriLLRh[i][j] - upLLR[i][j] << "   ";
					cout << endl;*/
					/*cout << oriLLRh[0][0] << endl;
					cout << upLLR[0][0] << endl;*/
					// test: oriLLR and upLLR
					// half iteration Enabled
					if (half_iter && iter == softiter)
					{
						if (iter == softiter) PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), Nh);
						for (int j = 0; j < Kv; j++)
							umid[Av[j]][decodingi] = u1[index][Av[j]];
					}
					// half iteration
					////////delete SCL variables
					for (int i = 0; i < L; i++)
						delete[]LLR[i];
					delete[] LLR;
					for (int i = 0; i < L; i++)
						delete[]sum[i];
					for (int i = 0; i < L; i++)
						delete[]u1[i];
					delete[] u1;
					delete[] sum;
					delete[] PM;
					delete[] LLR_in;
					delete[] W;
					delete[] SCL_index;
					delete[] better;
					delete[] worse;
					delete[] indexpm;

				}
				//cout << "Horizon Start:\n";
				//horizontal decoding
				// #
				if (weirditer) { continue; }

				// test: LLR after horizontal decoding
				/*cout << "After horizontal decoding" << endl;
				for (int i = 0; i < Nv; i++)
				{
					for (int j = 0; j < Nh; j++)
					{
						cout << (oriLLRh[i][j] < 0) << " ";
					}
					cout << endl;
				}*/
				// test: LLR after horizontal decoding
				// inner deinterleaving
				if (flag_interleave_inner)
				{
					float* tmp_llr_to_interleave = new float[Nv];
					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < Nv; j++) { tmp_llr_to_interleave[j] = oriLLRh[j][i]; }
						randDeIntrlvPartial(intrlv_pattern_2D_inner[i], tmp_llr_to_interleave, Av, Nv, Kv);
						for (int j = 0; j < Nv; j++) { oriLLRh[j][i] = tmp_llr_to_interleave[j]; }
					}
					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < Nv; j++) { tmp_llr_to_interleave[j] = upLLR[j][i]; }
						randDeIntrlvPartial(intrlv_pattern_2D_inner[i], tmp_llr_to_interleave, Av, Nv, Kv);
						for (int j = 0; j < Nv; j++) { upLLR[j][i] = tmp_llr_to_interleave[j]; }
					}
					delete[] tmp_llr_to_interleave;
				}


				/*auto decoding_over = std::chrono::system_clock::now();
				auto duration_decoding = std::chrono::duration_cast<std::chrono::microseconds>(decoding_over - decoding_start);
				cout << endl;
				cout << "Decoding duration" << duration_decoding.count() << "us" << endl;*/

#pragma omp parallel for
				for (int decodingi = 0; decodingi < Nv; decodingi++)
				{
					// 293.70kb
					// cout << "KKKK" << endl;
					if (half_iter && (iter == softiter)) { break; }
					string polarde_now = polarde_array[iter - 1];
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

					/*			if (L == 1)
								{
									//T_start = clock();
									Count = 0;
									decode(*LLR, *sum, (int)(log(Nh) / log(2)), A_Ach, flag, Nh, Kh);
									if (iter == softiter) PolarEncode_xor(*u1, *sum, Nh);
									//T_finish = clock();
								}
								else
								{*/
								//for (int k = 0; k < L; k++)  PM[k] = 0, indexpm[k] = k; 
								////T_start = clock();
								//decode_list(LLR, sum, (int)(log(Nh) / log(2)), A_Ach, 0, Nh, Kh, PM, L, (int)(log(L) / log(2)), LLR_in, W, SCL_index, better, worse, &Countv, &Countinfov, &Qsumunv);
								////T_finish = clock();
								////Decoding_Duration += (T_finish - T_start);
								//for (int i = 0; i < L - 1; i++)
								//	for (int j = i + 1; j < L; j++)
								//		if (PM[indexpm[i]] < PM[indexpm[j]])
								//		{
								//			int ttt = indexpm[i];
								//			indexpm[i] = indexpm[j];
								//			indexpm[j] = ttt;
								//		}


							/*
							* Encoding: Info -> Horizontal Encoding -> A -> set to 0 -> B -> Horizontal Encoding -> Vertical Encoding
							* Decoding: Received Info -> Horizontal Decoding -> Vertical Decoding -> Horizontal Decoding
							*/
					if (polarde_now == "SLD")
					{
						LSD_decode_outall(u_decod, u1_LSD, LLR[0], GGh, A_Ach, Nh, L_LSD, n_paths);
						// if (iter == softiter) PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), Nh);
						// for soft information exchange, inverse coding:
						// THIS LINE is WRONG! BUT IT SEEMS BETTER THAN THE RIGHT VERSION
						for (int i = 0; i < L_LSD; i++) { PolarEncode_xor(sum_LSD[i], u1_LSD[i], Nh); }
						// for (int i = 0; i < L_LSD; i++) { PolarEncode_xor(sum_LSD[i], u1_LSD[i], Nh); }
						// if (iter == softiter) PolarEncode_xor(u1[0], sum_LSD[0], Nh);
						index = indexpm[0];
						// test: after horizontal decoding
						/*if (decodingi == Av[0])
						{
							cout << "Direct output of horizontal decoding:" << endl;
							for (int i = 0; i < Nh; i++) { cout << u1_LSD[0][i] << " "; }
							cout << endl;
							cout << "horizontal decoder out, encode once:" << endl;
							for (int i = 0; i < Nh; i++) { cout << sum_LSD[0][i] << " "; }
							cout << endl;
						}*/
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
						// iter*2-1
						if (iter < softiter)
						{
							if (adaptive_beta)
							{
								if (softmethod == "Pyndiah") Pyndiahh_LSD_adaptivebeta(LLR, sum_LSD, oriLLRh, indexpm, iter, decodingi, upLLR);
								if (softmethod == "MITSO") MITSOh(LLR, PM, sum_LSD, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunv);
							}
							else
							{
								// cout << "IN LSD!!!" << endl;
								if (softmethod == "Pyndiah")
									Pyndiahh_LSD_getsum(LLR, sum_LSD, oriLLRh, indexpm, iter, decodingi, upLLR, sum_sdis_array[decodingi], GGh, Ach);
								//if (softmethod == "Pyndiah") Pyndiahh_LSD(LLR, sum_LSD, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR);
								if (softmethod == "MITSO") MITSOh(LLR, PM, sum_LSD, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunv);
							}
							//	getchar();
						}
						else
							for (int j = 0; j < Kh; j++)
								umid[decodingi][Ah[j]] = u1_LSD[index][Ah[j]];
					}
					if (polarde_now == "SCL")
					{
						for (int k = 0; k < L; k++)  PM[k] = 0, indexpm[k] = k;
						//T_start = clock();
						decode_list(LLR, sum, (int)(log(Nh) / log(2)), A_Ach, 0, Nh, Kh, PM, L, (int)(log(L) / log(2)), LLR_in, W, SCL_index, better, worse, &Countv, &Countinfov, &Qsumunv);
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
						if (iter == softiter) PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), Nh);
						index = indexpm[0];
						for (int col = 0; col < Nh; col++) { x_parity_check[decodingi][col] = sum[indexpm[0]][col]; }
						// show sum[indexpm[0]]
						/*cout << "sum[indexpm[0]]" << endl;
						for (int i = 0; i < Nh; i++) { cout << sum[indexpm[0]][i] << "  "; }
						cout << "sum[indexpm[0]]" << endl;*/
						// show sum[indexpm[0]]
						//}
						//calc distance s()
						if (iter < softiter)
						{
							if (adaptive_beta)
							{
								if (softmethod == "Pyndiah") Pyndiahh_adaptivebeta(LLR, sum, oriLLRh, indexpm, iter, decodingi, upLLR);
								if (softmethod == "MITSO") MITSOh(LLR, PM, sum, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunv);
								//	getchar();
							}
							else
							{
								//if (softmethod == "Pyndiah") Pyndiahh(LLR, sum, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR);
								if (softmethod == "Pyndiah") Pyndiahh_getsum(LLR, sum, oriLLRh, indexpm, iter, decodingi, upLLR, sum_sdis_array[decodingi]);
								if (softmethod == "MITSO") MITSOh(LLR, PM, sum, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunv);
								//	getchar();
							}
						}
						else
						{
							// prev
							/*for (int j = 0; j < Kh; j++)
								umid[decodingi][Ah[j]] = u1[index][Ah[j]];*/

							for (int j = 0; j < Nh; j++)
								umid[decodingi][j] = u1[index][j];
							// add result of current decoding to subcodeword_list
							for (int j = 0; j < L; j++)
							{
								int curr_index = indexpm[j];
								PolarEncode_xor(sub_codeword_list[decodingi] + j * Nh, sum[curr_index], Nh);
								edistance_list[decodingi][j] = PM[curr_index];
							}
						}

					}

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

				// list extension after horizontal decoding
				//if (iter == softiter && CRCL>0)
				//{
				//	crc_list_extend(sub_codeword_list, whole_codeword_list, edistance_list, L, Nv, Nh, L_WHOLE_CODEWORD);
				//	/*for (int i_subcode = 0; i_subcode < Nv; i_subcode++)
				//	{
				//		cout << endl << "subcodeword in whole codeword list" << endl;
				//		for (int i_bit = 0; i_bit < Nh; i_bit++)
				//			cout << whole_codeword_list[0][i_subcode * Nh + i_bit];
				//		cout << endl << "subcodeword in decoded result" << endl;
				//		for (int i_bit = 0; i_bit < Nh; i_bit++)
				//			cout << umid[i_subcode][i_bit];
				//	}
				//	cout << "TEST ENDS" << endl;*/
				//}






				if (flag_interleave_inner)
				{
					float* tmp_llr_to_interleave = new float[Nv];
					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < Nv; j++) { tmp_llr_to_interleave[j] = oriLLRh[j][i]; }
						randIntrlvPartial(intrlv_pattern_2D_inner[i], tmp_llr_to_interleave, Av, Nv, Kv);
						for (int j = 0; j < Nv; j++) { oriLLRh[j][i] = tmp_llr_to_interleave[j]; }
					}
					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < Nv; j++) { tmp_llr_to_interleave[j] = upLLR[j][i]; }
						randIntrlvPartial(intrlv_pattern_2D_inner[i], tmp_llr_to_interleave, Av, Nv, Kv);
						for (int j = 0; j < Nv; j++) { upLLR[j][i] = tmp_llr_to_interleave[j]; }
					}
					delete[] tmp_llr_to_interleave;
				}
			}
			int calcsum[Nv], calcu1[Nv];
			int calcu1h[Nh];
			// calcsum is from the direct output of the decoder
			// if output is the result of row decoding (meaning there is another column decoding to be done)
			if (half_iter)
			{
				for (int i = 0; i < Kv; i++)
				{
					PolarEncode_xor(calcu1h, umid[Av[i]], Nh);
					for (int j = 0; j < Kh; j++)
						uout[Av[i]][Ah[j]] = calcu1h[Ah[j]];
				}
			}
			// if output is the result of column decoding (meaning that there is another row decoding to be done)
			else
			{

				// prev
				/*for (int i = 0; i < Kh; i++)
				{
					for (int j = 0; j < Nv; j++)
						calcsum[j] = umid[j][Ah[i]];
					PolarEncode_xor(calcu1, calcsum, Nv);
					for (int j = 0; j < Kv; j++)
						uout[Av[j]][Ah[i]] = calcu1[Av[j]];
				}*/

				// the whole codeword
				for (int i = 0; i < Nh; i++)
				{
					//cout << endl;
					for (int j = 0; j < Nv; j++)
					{
						calcsum[j] = umid[j][i];
						//cout << umid[j][i] <<"  ";
					}
					for (int i_calcu1 = 0; i_calcu1 < Nv; i_calcu1++) { calcu1[i_calcu1] = 0; }
					PolarEncode_xor(calcu1, calcsum, Nv);
					for (int j = 0; j < Nv; j++)
						uout[j][i] = calcu1[j];
				}
				//cout << endl << "END TEST UMID" << endl;
				// output uout
				/*for (int i = 0; i < Nv; i++)
				{
					cout << endl;
					for (int j = 0; j < Nh; j++) { cout << uout[i][j] << "  "; }
				}*/
				// cout << endl;
				if (CRCL > 0)
				{
					// conduct crc check on info bits of uout


					flag_pass_crc = false;
					flag_found_replacement_codeword = false;

					for (int i = 0; i < Kv; i++)
					{
						for (int j = 0; j < Kh; j++)
							u_crc[i * Kh + j] = uout[Av[i]][Ah[j]];
					}
					if (rx_check_crc(u_crc, Kh * Kv, CRCL)) {
						flag_pass_crc = true;
					}

					if (!flag_pass_crc)
					{
						crc_list_extend(sub_codeword_list, whole_codeword_list, edistance_list, L, Nv, Nh, L_WHOLE_CODEWORD);
						for (int i_code = 0; i_code < L_WHOLE_CODEWORD; i_code++)
						{
							for (int i = 0; i < Kh; i++)
							{
								for (int j = 0; j < Nv; j++)
									calcsum[j] = whole_codeword_list[i_code][j * Nh + Ah[i]];
								PolarEncode_xor(calcu1, calcsum, Nv);
								for (int j = 0; j < Kv; j++)
									uout[Av[j]][Ah[i]] = calcu1[Av[j]];
							}
							// copy from uout to u_crc
							for (int i = 0; i < Kv; i++)
							{
								for (int j = 0; j < Kh; j++)
									u_crc[i * Kh + j] = uout[Av[i]][Ah[j]];
							}
							if (rx_check_crc(u_crc, Kh * Kv, CRCL)) {
								if (i_code == 0) { flag_pass_crc = true; }
								if (i_code > 0) { cout << endl << "Found in " << i_code << endl; flag_found_replacement_codeword = true; idx_replacement_codeword = i_code; }
								break;
							}
							//if (i_code == 0) { break; }
						}
					}
				}
			}
			// test: final output
			/*cout << "final output" << endl;
			for (int k = 0; k < Nh; k++) { cout << uout[Av[0]][k] << " "; }
			cout << endl;*/
			// test: final output
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


			// whether pass parity check for both row and column
			bool flag_pass_parity_check = true;

			//  check with x_parity_check
			//  half_iter : result is column decoding
			if (half_iter) {
				int* row_cwd = new int[Nh];
				int* row_cwd_encoded = new int[Nh];
				for (int i = 0; i < Nv; i++) {
					for (int j = 0; j < Nh; j++) {
						row_cwd[j] = x_parity_check[i][j];
					}
					PolarEncode_xor(row_cwd_encoded, row_cwd, Nh);
					for (int k = 0; k < Nh - Kh; k++) { if (row_cwd_encoded[Ach[k]] != 0) flag_pass_parity_check = false; }
				}
				delete[] row_cwd;
				delete[] row_cwd_encoded;
			}
			else {
				for (int i = 0; i < Nh; i++)
				{
					//cout << endl;
					for (int j = 0; j < Nv; j++)
					{
						calcsum[j] = x_parity_check[j][i];
						//cout << umid[j][i] <<"  ";
					}
					for (int i_calcu1 = 0; i_calcu1 < Nv; i_calcu1++) { calcu1[i_calcu1] = 0; }
					PolarEncode_xor(calcu1, calcsum, Nv);
					// cout << endl;
					// for (int j = 0; j < Nv; j++) { cout << calcu1[j] << "  "; }
					for (int k = 0; k < Nv - Kv; k++) { if (calcu1[Acv[k]] != 0) flag_pass_parity_check = false; }
				}
			}



			/*for (int i = 0; i < Nv - Kv; i++)
			{
				for (int j = 0; j < Nh; j++) { if (uout[Acv[i]][j] != 0) { flag_pass_parity_check = false; break; } }
				if (!flag_pass_parity_check) break;
			}
			for (int i = 0; i<Nh - Kh; i++)
			{
				for (int j = 0; j < Nv; j++) { if (uout[j][Ach[i]] != 0) { flag_pass_parity_check = false; break; } }
				if (!flag_pass_parity_check) break;
			}*/

			// cout << "PARITY" << flag_pass_parity_check << endl;

			/*if (!flag_pass_crc) {
				cout << endl;
			}*/
			//getchar();
			if (flag_found_replacement_codeword)
			{
				file_crc_check_replacement_codeword << endl;
				if (Temp_Error < Num_Error) { file_crc_check_replacement_codeword << "Index of Replacement Codeword: " << idx_replacement_codeword << "IsCorrect: " << 0 << endl; }
				if (Temp_Error == Num_Error) { file_crc_check_replacement_codeword << "Index of Replacement Codeword: " << idx_replacement_codeword << "IsCorrect: " << 1 << endl; }
			}
			if (Temp_Error < Num_Error)
			{
				Num_Frame_Error++;
				if (flag_pass_crc) { n_err_frame_pass_crc++; }
				else { n_err_frame_uncorrected_by_crc++; }
				if (flag_found_replacement_codeword) {
					n_err_replacement_codeword_crc++;
				}
				if (flag_pass_parity_check) { n_err_pass_parity_check++; }
			}
			if (Temp_Error == Num_Error)
			{
				/*if (!flag_pass_parity_check) {
					cout << endl << "FFFFFFFUCKKKKKK" << endl;
					for (int i = 0; i < Nv; i++)
					{
						cout << endl;
						for (int j = 0; j < Nh; j++) { cout << uout[i][j] << "  "; }
					}
					cout << endl;
				}*/
				if (!flag_pass_crc) { n_err_frame_corrected_by_crc++; }
			}
			if (frame % 10 == 0)
			{
				/*printf("\r<%d> SNR= %.2f FER= %d / %d = %f BER= %d / %d = %lf E-Det = %f, E-MRB = %f"  ,
					frame, nSNR, Num_Frame_Error, frame, (double)Num_Frame_Error / frame, Num_Error, (N_Info * frame), (double)Num_Error / N_Info / frame,
					(double)n_det_err/total_bit_cnt, (double)n_det_err_mrb/(total_bit_cnt/(ModType/2)));*/
				int L_sdis_sum = L_LSD;
				sdis_sum = 0;
				for (int i = 0; i < Nv; i++) sdis_sum += sum_sdis_array[i];
				if (polarde_array[0] == "SCL") L_sdis_sum = L;
				// cout << "sdis_sum_avg" << sdis_sum / frame / L_sdis_sum << endl;
				printf("\r<%d> SNR= %.2f FER= %d / %d = %f BER= %d / %d = %lf ,sdis_avg = %f, E-MRB = %f",
					frame, nSNR, Num_Frame_Error, frame, (double)Num_Frame_Error / frame, Num_Error, (N_Info * frame), (double)Num_Error / N_Info / frame, (double)sdis_sum / frame / L_sdis_sum, (double)n_det_err / total_bit_cnt);
				fflush(stdout);
			}


			// output interleaving pattern
			//if (Num_Frame_Error % record_frames == 0 && Num_Frame_Error > 0 && Num_Frame_Error!=n_frame_error_pre)
			//{
			//	n_frame_error_pre = Num_Frame_Error;
			//	cout << endl;
			//	int n_file = Num_Frame_Error / record_frames;
			//	frame_delta = frame - frame_tmp;
			//	cout << "Delta Frames: " << frame_delta << endl;
			//	frame_tmp = frame;
			//	string ofile_intrlv_name = string("intrlv_pattern_") + to_string(n_file) + string("_frame_") + to_string(frame_delta) + string(".txt");
			//	string ofile_spatial_intrlv_name = string("intrlv_pattern_spatial") + to_string(n_file) + string(".txt");
			//	ofstream ofile_intrlv(ofile_intrlv_name);
			//	ofstream ofile_spatial_intrlv(ofile_spatial_intrlv_name);
			//	for (int i = 0; i < Nv; i++)
			//	{
			//		for (int j = 0; j < intrlv_pattern_num; j++)
			//		{

			//			for (int k = 0; k < ModType; k++)
			//				ofile_intrlv << intrlv_pattern_2D_time_partial[i][j][k] << " ";
			//			//cout << endl;
			//		}
			//	}
			//	ofile_intrlv.close();

			//	for (int i = 0; i < Nh; i++)
			//	{
			//		for (int j = 0; j < Nv; j++)
			//			ofile_spatial_intrlv << intrlv_pattern_2D[i][j] << " ";
			//	}
			//	ofile_spatial_intrlv.close();

			//	// re-generate interleaving pattern
			//	generateIntrlvPattern(intrlv_pattern, Nv);
			//	for (int i = 0; i < Nh; i++) { generateIntrlvPattern(intrlv_pattern_2D[i], Nv); }
			//	for (int i = 0; i < Nv; i++) { generateIntrlvPattern(intrlv_pattern_2D_time[i], Nh); }
			//	if (flag_interleave_partial)
			//	{
			//		for (int i = 0; i < Nv; i++)
			//		{
			//			for (int j = 0; j < intrlv_pattern_num; j++)
			//			{
			//				generateIntrlvPattern(intrlv_pattern_2D_time_partial[i][j], ModType);
			//				/*for (int k = 0; k < ModType; k++)
			//					cout << intrlv_pattern_2D_time_partial[i][j][k];
			//				cout << endl;*/
			//			}
			//		}
			//	}
			//	//Num_Frame_Error = 0;
			//}


			//printf("%d %d", Num_Frame_Error, Num_Error);
			//getchar();
			//Decoding_Duration += (T_finish - T_start);
			/*auto dec_end = std::chrono::system_clock::now();
			auto duration_decoding = std::chrono::duration_cast<std::chrono::microseconds>(dec_end - det_end);
			auto duration_total = std::chrono::duration_cast<std::chrono::microseconds>(dec_end - start);
			file_dec << duration_decoding.count() << endl;
			file_total << duration_total.count() << endl;*/
			if (frame % 16 == 0)
			{
				file_err_num << Num_Error << endl;
			}


			// logging crc output
			if ((frame % 1000) == 0)
			{
				file_crc_check_decoding << endl;
				file_crc_check_decoding << "FER: " << (double)Num_Frame_Error / frame << endl;
				file_crc_check_decoding << "Total Frames: " << frame << endl;
				file_crc_check_decoding << "Error frames undetected by crc: " << n_err_frame_pass_crc << endl;
				file_crc_check_decoding << "Error frames corrected by crc: " << n_err_frame_corrected_by_crc << endl;
				file_crc_check_decoding << "Error frames failed to be corrected by crc: " << n_err_frame_uncorrected_by_crc << endl;
				file_crc_check_decoding << "Error of replacement codewords: " << n_err_replacement_codeword_crc << endl;
				file_crc_check_decoding << "Error frames that pass parity check" << n_err_pass_parity_check << endl;
			}

		}
		//Decoding_Duration = Decoding_Duration * 100 / CLOCKS_PER_SEC / frame;
		//Decoding_Duration = Decoding_Duration / CLOCKS_PER_SEC / Max_frame;

		//printf("%f  %f  %f  %f\n", nSNR, Decoding_Duration, (double)Num_Error / N_Info / Max_frame, (double)Num_Frame_Error / Max_frame);
		printf("%f %f %f\n", nSNR, (double)Num_Frame_Error / frame, (double)Num_Error / N_Info / frame);
		BER_array[SNR_count] = (double)Num_Error / N_Info / frame;
		FER_array[SNR_count] = (double)Num_Frame_Error / frame;
		// save to .txt file
		string polarde_serial = "";
		for (int i = 0; i < softiter; i++)
			polarde_serial = polarde_serial + polarde_array[i] + string("_");
		string file_name = string("ResFinal\\2D_LSD_D_") + polarde_serial + MIMODET + string("T") + to_string(Nt) + string("R") + to_string(Nr)
			+ string("Q") + to_string(int(pow(2, ModType))) + string("L") + to_string(Nh) + string("K") + to_string(Kh)
			+ string("T") + to_string(softiter) + string("H") + to_string(half_iter) + string("UB") + to_string(usebeta) + string("I") + to_string(flag_interleave)
			+ string(".txt");
		save_array_txt(nSNR, BER_array[SNR_count], FER_array[SNR_count], file_name);
		int L_sdis_sum = L_LSD;
		sdis_sum = 0;
		for (int i = 0; i < Nv; i++) sdis_sum += sum_sdis_array[i];
		if (polarde_array[0] == "SCL") L_sdis_sum = L;
		cout << "sdis_sum_avg" << sdis_sum / frame / L_sdis_sum << endl;
		SNR_count++;


		// logging crc details
		if (true)
		{
			file_crc_check_decoding << endl;
			file_crc_check_decoding << "FER: " << (double)Num_Frame_Error / frame << endl;
			file_crc_check_decoding << "Total Frames: " << frame << endl;
			file_crc_check_decoding << "Error frames undetected by crc: " << n_err_frame_pass_crc << endl;
			file_crc_check_decoding << "Error frames corrected by crc: " << n_err_frame_corrected_by_crc << endl;
			file_crc_check_decoding << "Error frames failed to be detected by crc: " << n_err_frame_uncorrected_by_crc << endl;
			file_crc_check_decoding << "Error of replacement codewords: " << n_err_replacement_codeword_crc << endl;
		}

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
	for (int i = 0; i < Nh; i++)
		delete[] oriLLR_kbest[i];
	delete[] oriLLR_kbest;
	for (int i = 0; i < Nv; i++)
		delete sub_codeword_list[i];
	delete[] sub_codeword_list;
	for (int i = 0; i < L_WHOLE_CODEWORD; i++)
		delete whole_codeword_list[i];
	delete[] whole_codeword_list;
	for (int i = 0; i < Nv; i++)
		delete[]edistance_list[i];
	delete[] edistance_list;
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
	for (int i = 0; i < mimo_blk_height; i++) {
		delete[] x_mod[i];
	}
	delete[] x_mod;
	for (int i = 0; i < Nr; i++) {
		delete[] y_mod[i];
	}
	delete[] y_mod;
	for (int i = 0; i < 2 * mimo_blk_height; i++) {
		delete[] y_mmse[i];
	}
	delete[] y_mmse;
	for (int i = 0; i < 2 * Nr; i++) {
		delete[] y_mod_real[i];
	}
	for (int i = 0; i < mimo_blk_length; i++)
	{
		delete[] y_mod_real_temp_omp[i];
	}
	delete[] y_mod_real_temp_omp;
	for (int i = 0; i < mimo_blk_length; i++)
	{
		delete[] y_mmse_temp_omp[i];
	}
	delete[] y_mmse_temp_omp;
	for (int i = 0; i < mimo_blk_length; i++)
	{
		delete[] ep_extr_mean_omp[i];
	}
	delete[] ep_extr_mean_omp;
	for (int i = 0; i < mimo_blk_length; i++)
	{
		delete[] ep_extr_var_omp[i];
	}
	delete[] ep_extr_var_omp;
	delete[] y_mod_real;
	delete[] x_mod_tmp;
	delete[] y_mod_tmp;
	delete[] y_mod_real_tmp;
	delete[] x_tmp;
	delete[] oriLLR_tmp;
	delete[] sum_sdis_array;
	delete[] Hreal_array;
	delete[] HTH_array;
	delete[] HTY_array;
	delete[] intrlv_pattern_2D_inner;
}
void system2D_LSD_sys(double* BER_array, double* FER_array)
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
	//Eigen::Matrix<float, 2 * Nr, 2 * Nt> H_real; // ʵ���ŵ�����
	//Eigen::Matrix<float, 2 * Nr, 1> received_signal; // �����źž���
	//Eigen::Matrix<float, 2 * Nt, 2 * Nt> HTH; // HTH����
	//Eigen::Matrix<float, 2 * Nt, 2 * Nt> Identity = Eigen::Matrix<float, 2 * Nt, 2 * Nt>::Identity(); // ��λ����
	//Eigen::Matrix<float, 2 * Nt, 1> output_signal; // ����źž���
	//Eigen::Matrix<float, 2 * Nt, 2 * Nr> mmse_matrix; // MMSE����


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
			// test:original info
			/*cout << "original info u:" << endl;
			for (int j = 0; j < Nh; j++) { cout << u[Av[0]][j] << " "; }
			cout << endl;*/
			// Polar Encoding

#pragma omp parallel for
			for (int i = 0; i < Nv; i++)
				PolarEncode_xor(midx[i], u[i], Nh);
			// test:horizontal encoding
			/*cout << "u after horizontal encoding:" << endl;
			for (int j = 0; j < Nh; j++) { cout << midx[Av[0]][j] << " "; }
			cout << endl;*/
			// test:horizontal encoding
			for (int i = 0; i < Nv; i++)
			{
				int* code_tmp = new int[Nh];
				// TEST:
				/*if (i == Av[0])
				{
					cout << "After Stage 0 Encoding(Horizontal encoding, before setting to 0) \n";
					for (int i = 0; i < Nh; i++) { cout << " " << midx[Av[0]][i]; }
					cout << endl;
				}*/
				// TEST
				for (int j = 0; j < Nh - Kh; j++) { midx[i][Ach[j]] = 0; }
				// TEST:
				/*if (i == 0)
				{
					cout << "After Stage 1 Encoding(Setting to 0) \n";
					for (int i = 0; i < Nh; i++) { cout << " " << midx[Av[0]][i]; }
					cout << endl;
				}*/
				// TEST
				for (int j = 0; j < Nh; j++) { code_tmp[j] = midx[i][j]; }
				PolarEncode_xor(midx[i], code_tmp, Nh);
				delete[] code_tmp;
			}
			// TEST:
			/*cout << "After Stage 2 Encoding(Horizontal Encoding, after setting to 0) \n";
			for (int i = 0; i < Nh; i++) { cout << " " << midx[Av[0]][i]; }
			cout << endl;
			cout << "original info" << "\n";
			for (int i = 0; i < Nh; i++) { cout << " " << u[Av[0]][i]; }*/
			// TEST
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
			// test:vertical encoding
			/*cout << "u after vertical encoding:" << endl;
			for (int j = 0; j < Nh; j++) { cout << x[Av[0]][j] << " "; }
			cout << endl;*/
			// test:vertical encoding
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
						// TEST�� det error
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
			for (int iter = 1; iter <= softiter; iter++)
			{
				//	cout << "Vertical Start:\n";
					//vertical decoding

#pragma omp parallel for
				for (int decodingi = 0; decodingi < Nh; decodingi++)
				{

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
					/*			if (L == 1)
								{
									//T_start = clock();
									Count = 0;
									decode(*LLR, *sum, (int)(log(Nv) / log(2)), A_Acv, flag, Nv, Kv);
									//PolarEncode_xor(*u1, *sum, N);
									//T_finish = clock();
								}
								else
								{*/
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
				//cout << "Horizon Start:\n";
				//horizontal decoding
				// #pragma omp parallel for
				for (int decodingi = 0; decodingi < Nv; decodingi++)
				{
					string polarde_now = polarde_array[iter - 1];
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

					if (polarde_now == "SLD")
					{
						LSD_decode_outall(u_decod, u1_LSD, LLR[0], GGh, A_Ach, Nh, L_LSD, n_paths);
						// if (iter == softiter) PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), Nh);
						// for soft information exchange, inverse coding:
						// THIS LINE is WRONG! BUT IT SEEMS BETTER THAN THE RIGHT VERSION
						for (int i = 0; i < L_LSD; i++) { PolarEncode_xor(sum_LSD[i], u1_LSD[i], Nh); }
						// if (iter == softiter) PolarEncode_xor(u1[0], sum_LSD[0], Nh);
						index = indexpm[0];
						if (iter < softiter)
						{
							if (softmethod == "Pyndiah") Pyndiahh_LSD(LLR, sum_LSD, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR);
							if (softmethod == "MITSO") MITSOh(LLR, PM, sum_LSD, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR, Qsumunv);
							//	getchar();
						}
						else
							for (int j = 0; j < Kh; j++)
								umid[decodingi][Ah[j]] = u1_LSD[index][Ah[j]];
					}

					if (polarde_now == "SCL")
					{
						for (int k = 0; k < L; k++)  PM[k] = 0, indexpm[k] = k;
						//T_start = clock();
						decode_list(LLR, sum, (int)(log(Nh) / log(2)), A_Ach, 0, Nh, Kh, PM, L, (int)(log(L) / log(2)), LLR_in, W, SCL_index, better, worse, &Countv, &Countinfov, &Qsumunv);
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
						if (iter == softiter) PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), Nh);
						index = indexpm[0];
						//}
						//calc distance s()
						if (iter < softiter)
						{
							if (softmethod == "Pyndiah") Pyndiahh(LLR, sum, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR);
							if (softmethod == "MITSO") MITSOh(LLR, PM, sum, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR, Qsumunv);
							//	getchar();
						}
						else
							for (int j = 0; j < Kh; j++)
								umid[decodingi][Ah[j]] = u1[index][Ah[j]];
					}
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
			int calcsum[Nv], calcu1[Nv];
			// calcsum is from the direct output of the decoder
			for (int i = 0; i < Kh; i++)
			{
				for (int j = 0; j < Nv; j++)
					calcsum[j] = umid[j][Ah[i]];
				PolarEncode_xor(calcu1, calcsum, Nv);
				for (int j = 0; j < Kv; j++)
					uout[Av[j]][Ah[i]] = calcu1[Av[j]];
			}
			int calcsumh[Nh];
			for (int i = 0; i < Kv; i++)
			{
				for (int j = 0; j < Nh; j++) { calcsumh[j] = uout[Av[i]][j]; }
				PolarEncode_xor(uout[Av[i]], calcsumh, Nh);
			}
			// test: final output
			/*cout << "final output" << endl;
			for (int k = 0; k < Nh; k++) { cout << uout[Av[0]][k] << " "; }
			cout << endl;*/
			// test: final output
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
				printf("\r<%d> SNR= %.2f FER= %d / %d = %f BER= %d / %d = %lf ",
					frame, nSNR, Num_Frame_Error, frame, (double)Num_Frame_Error / frame, Num_Error, (N_Info * frame), (double)Num_Error / N_Info / frame);
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


void system2D(double* BER_array, double* FER_array)
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
	//Eigen::Matrix<float, 2 * Nr, 2 * Nt> H_real; // ʵ���ŵ�����
	//Eigen::Matrix<float, 2 * Nr, 1> received_signal; // �����źž���
	//Eigen::Matrix<float, 2 * Nt, 2 * Nt> HTH; // HTH����
	//Eigen::Matrix<float, 2 * Nt, 2 * Nt> Identity = Eigen::Matrix<float, 2 * Nt, 2 * Nt>::Identity(); // ��λ����
	//Eigen::Matrix<float, 2 * Nt, 1> output_signal; // ����źž���
	//Eigen::Matrix<float, 2 * Nt, 2 * Nr> mmse_matrix; // MMSE����


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
	// srand((unsigned int)time(0));
	std::random_device rd;              // ��ȡ���������
	std::mt19937 gen(rd());             // ʹ�����ӳ�ʼ��mt19937������


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
						// TEST�� det error
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
			for (int iter = 1; iter <= softiter; iter++)
			{
				//	cout << "Vertical Start:\n";
					//vertical decoding

#pragma omp parallel for
				for (int decodingi = 0; decodingi < Nh; decodingi++)
				{

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
					/*			if (L == 1)
								{
									//T_start = clock();
									Count = 0;
									decode(*LLR, *sum, (int)(log(Nv) / log(2)), A_Acv, flag, Nv, Kv);
									//PolarEncode_xor(*u1, *sum, N);
									//T_finish = clock();
								}
								else
								{*/
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
				//cout << "Horizon Start:\n";
				//horizontal decoding
#pragma omp parallel for
				for (int decodingi = 0; decodingi < Nv; decodingi++)
				{
					//for SCL Decoding
					int Countv = 0, Countinfov = 0;
					float Qsumunv = 0;
					int** u1 = new int* [L];
					for (int i = 0; i < L; i++)  u1[i] = new _MM_ALIGN16 int[Nmax];
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
					for (int j = 0; j < Nh; j++)
						for (int k = 0; k < L; k++)
							LLR[k][j] = upLLR[decodingi][j];
					//	cout << LLR[0][j] << " ";

				//	cout << endl;
					int index = 0;
					/*			if (L == 1)
								{
									//T_start = clock();
									Count = 0;
									decode(*LLR, *sum, (int)(log(Nh) / log(2)), A_Ach, flag, Nh, Kh);
									if (iter == softiter) PolarEncode_xor(*u1, *sum, Nh);
									//T_finish = clock();
								}
								else
								{*/
					for (int k = 0; k < L; k++)  PM[k] = 0, indexpm[k] = k;
					//T_start = clock();
					decode_list(LLR, sum, (int)(log(Nh) / log(2)), A_Ach, 0, Nh, Kh, PM, L, (int)(log(L) / log(2)), LLR_in, W, SCL_index, better, worse, &Countv, &Countinfov, &Qsumunv);
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
					if (iter == softiter) PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), Nh);
					index = indexpm[0];
					//}
					//calc distance s()
					if (iter < softiter)
					{
						if (softmethod == "Pyndiah") Pyndiahh(LLR, sum, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR);
						if (softmethod == "MITSO") MITSOh(LLR, PM, sum, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR, Qsumunv);
						//	getchar();
					}
					else
						for (int j = 0; j < Kh; j++)
							umid[decodingi][Ah[j]] = u1[index][Ah[j]];
					////////delete SCL variables
					for (int i = 0; i < L; i++)
						delete[]u1[i];
					delete[] u1;
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
			}
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
				/*printf("\r<%d> SNR= %.2f FER= %d / %d = %f BER= %d / %d = %lf  Hard BER = %d / %d = %f",
					frame, nSNR, Num_Frame_Error, frame, (double)Num_Frame_Error / frame, Num_Error, (N_Info * frame), (double)Num_Error / N_Info / frame, n_det_err, total_bit_cnt, (double)n_det_err / total_bit_cnt);
				fflush(stdout);*/
				//printf("\r<%d> SNR= %.2f FER= %d / %d = %f BER= %d / %d = %lf  Hard BER = %d / %d = %f",
				//	frame, nSNR, Num_Frame_Error, frame, (double)Num_Frame_Error / frame, Num_Error, (N_Info* frame), (double)Num_Error / N_Info / frame/*n_det_err, total_bit_cnt, (double)n_det_err / total_bit_cnt*/);
				//fflush(stdout);
				printf("\r<%d> SNR= %.2f FER= %d / %d = %f BER= %d / %d = %lf",
					frame, nSNR, Num_Frame_Error, frame, (double)Num_Frame_Error / frame, Num_Error, (N_Info * frame), (double)Num_Error / N_Info / frame/*n_det_err, total_bit_cnt, (double)n_det_err / total_bit_cnt*/);
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
	//Eigen::Matrix<float, 2 * Nr, 2 * Nt> H_real; // ʵ���ŵ�����
	//Eigen::Matrix<float, 2 * Nr, 1> received_signal; // �����źž���
	//Eigen::Matrix<float, 2 * Nt, 2 * Nt> HTH; // HTH����
	//Eigen::Matrix<float, 2 * Nt, 2 * Nt> Identity = Eigen::Matrix<float, 2 * Nt, 2 * Nt>::Identity(); // ��λ����
	//Eigen::Matrix<float, 2 * Nt, 1> output_signal; // ����źž���
	//Eigen::Matrix<float, 2 * Nt, 2 * Nr> mmse_matrix; // MMSE����


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
						// TEST�� det error
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
			//Polar Decoding for each blk
			int blk_len = Nh / blk_cnt;
			for (int i_decod_blk = 0; i_decod_blk < blk_cnt; i_decod_blk++)
			{
				// column (spatial) decoding 
				int decodingi = 0;
				for (int iter = 1; iter <= softiter; iter++)
				{
					for (int i_decod_col = 0; i_decod_col < (i_decod_blk + 1) * blk_len; i_decod_col++)
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


/*
* 2D system with double stage vertical decoding
*/
void system2D_LSD_doublestage(double* BER_array, double* FER_array)
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
	//Eigen::Matrix<float, 2 * Nr, 2 * Nt> H_real; // ʵ���ŵ�����
	//Eigen::Matrix<float, 2 * Nr, 1> received_signal; // �����źž���
	//Eigen::Matrix<float, 2 * Nt, 2 * Nt> HTH; // HTH����
	//Eigen::Matrix<float, 2 * Nt, 2 * Nt> Identity = Eigen::Matrix<float, 2 * Nt, 2 * Nt>::Identity(); // ��λ����
	//Eigen::Matrix<float, 2 * Nt, 1> output_signal; // ����źž���
	//Eigen::Matrix<float, 2 * Nt, 2 * Nr> mmse_matrix; // MMSE����


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
			// test:original info
			/*cout << "original info u:" << endl;
			for (int j = 0; j < Nh; j++) { cout << u[Av[0]][j] << " "; }*/
			//cout << endl;
			// Polar Encoding
#pragma omp parallel for
			for (int i = 0; i < Nv; i++)
				PolarEncode_xor(midx[i], u[i], Nh);
			// test:horizontal encoding
			/*cout << "u after horizontal encoding:" << endl;
			for (int j = 0; j < Nh; j++) { cout << midx[Av[0]][j] << " "; }
			cout << endl;*/
			// test:horizontal encoding
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
			// test:vertical encoding
			/*cout << "u after vertical encoding:" << endl;
			for (int j = 0; j < Nh; j++) { cout << x[Av[0]][j] << " "; }
			cout << endl;*/
			// test:vertical encoding
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
						// TEST�� det error
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
			for (int iter = 1; iter <= softiter; iter++)
			{
				//	cout << "Vertical Start:\n";
					//vertical decoding

#pragma omp parallel for
				for (int decodingi = 0; decodingi < Nh; decodingi++)
				{

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
					/*			if (L == 1)
								{
									//T_start = clock();
									Count = 0;
									decode(*LLR, *sum, (int)(log(Nv) / log(2)), A_Acv, flag, Nv, Kv);
									//PolarEncode_xor(*u1, *sum, N);
									//T_finish = clock();
								}
								else
								{*/
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
					if (adaptive_beta)
					{
						if (softmethod == "Pyndiah") Pyndiahv_adaptivebeta(LLR, sum, oriLLRh, indexpm, iter * 2 - 2, decodingi, upLLR);
						if (softmethod == "MITSO") MITSOv(LLR, PM, sum, oriLLRh, indexpm, iter * 2 - 2, decodingi, upLLR, Qsumunh);
					}
					else
					{
						if (softmethod == "Pyndiah") Pyndiahv(LLR, sum, oriLLRh, indexpm, iter * 2 - 2, decodingi, upLLR);
						if (softmethod == "MITSO") MITSOv(LLR, PM, sum, oriLLRh, indexpm, iter * 2 - 2, decodingi, upLLR, Qsumunh);
					}
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
				//cout << "Horizon Start:\n";
				//horizontal decoding
				// #pragma omp parallel for
				for (int decodingi = 0; decodingi < Nv; decodingi++)
				{
					string polarde_now = polarde_array[iter - 1];
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

					/*			if (L == 1)
								{
									//T_start = clock();
									Count = 0;
									decode(*LLR, *sum, (int)(log(Nh) / log(2)), A_Ach, flag, Nh, Kh);
									if (iter == softiter) PolarEncode_xor(*u1, *sum, Nh);
									//T_finish = clock();
								}
								else
								{*/
								//for (int k = 0; k < L; k++)  PM[k] = 0, indexpm[k] = k; 
								////T_start = clock();
								//decode_list(LLR, sum, (int)(log(Nh) / log(2)), A_Ach, 0, Nh, Kh, PM, L, (int)(log(L) / log(2)), LLR_in, W, SCL_index, better, worse, &Countv, &Countinfov, &Qsumunv);
								////T_finish = clock();
								////Decoding_Duration += (T_finish - T_start);
								//for (int i = 0; i < L - 1; i++)
								//	for (int j = i + 1; j < L; j++)
								//		if (PM[indexpm[i]] < PM[indexpm[j]])
								//		{
								//			int ttt = indexpm[i];
								//			indexpm[i] = indexpm[j];
								//			indexpm[j] = ttt;
								//		}


							/*
							* Encoding: Info -> Horizontal Encoding -> A -> set to 0 -> B -> Horizontal Encoding -> Vertical Encoding
							* Decoding: Received Info -> Horizontal Decoding -> Vertical Decoding -> Horizontal Decoding
							*/
					if (polarde_now == "SLD")
					{
						LSD_decode_outall(u_decod, u1_LSD, LLR[0], GGh, A_Ach, Nh, L_LSD, n_paths);
						// if (iter == softiter) PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), Nh);
						// for soft information exchange, inverse coding:
						// THIS LINE is WRONG! BUT IT SEEMS BETTER THAN THE RIGHT VERSION
						for (int i = 0; i < L_LSD; i++) { PolarEncode_xor(sum_LSD[i], u1_LSD[i], Nh); }
						// if (iter == softiter) PolarEncode_xor(u1[0], sum_LSD[0], Nh);
						index = indexpm[0];
						// test: after horizontal decoding
						/*if (decodingi == Av[0])
						{
							cout << "Direct output of horizontal decoding:" << endl;
							for (int i = 0; i < Nh; i++) { cout << u1_LSD[0][i] << " "; }
							cout << endl;
							cout << "horizontal decoder out, encode once:" << endl;
							for (int i = 0; i < Nh; i++) { cout << sum_LSD[0][i] << " "; }
							cout << endl;
						}*/
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
							if (adaptive_beta)
							{
								if (softmethod == "Pyndiah") Pyndiahh_LSD_adaptivebeta(LLR, sum_LSD, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR);
								if (softmethod == "MITSO") MITSOh(LLR, PM, sum_LSD, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR, Qsumunv);
							}
							else
							{
								if (softmethod == "Pyndiah") Pyndiahh_LSD(LLR, sum_LSD, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR);
								if (softmethod == "MITSO") MITSOh(LLR, PM, sum_LSD, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR, Qsumunv);
							}
							//	getchar();
						}
						else
							for (int j = 0; j < Kh; j++)
								umid[decodingi][Ah[j]] = u1_LSD[index][Ah[j]];
					}

					if (polarde_now == "SCL")
					{
						for (int k = 0; k < L; k++)  PM[k] = 0, indexpm[k] = k;
						//T_start = clock();
						decode_list(LLR, sum, (int)(log(Nh) / log(2)), A_Ach, 0, Nh, Kh, PM, L, (int)(log(L) / log(2)), LLR_in, W, SCL_index, better, worse, &Countv, &Countinfov, &Qsumunv);
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
						if (iter == softiter) PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), Nh);
						index = indexpm[0];
						//}
						//calc distance s()
						if (iter < softiter)
						{
							if (adaptive_beta)
							{
								if (softmethod == "Pyndiah") Pyndiahh_adaptivebeta(LLR, sum, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR);
								if (softmethod == "MITSO") MITSOh(LLR, PM, sum, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR, Qsumunv);
								//	getchar();
							}
							else
							{
								if (softmethod == "Pyndiah") Pyndiahh(LLR, sum, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR);
								if (softmethod == "MITSO") MITSOh(LLR, PM, sum, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR, Qsumunv);
								//	getchar();
							}
						}
						else
							for (int j = 0; j < Kh; j++)
								umid[decodingi][Ah[j]] = u1[index][Ah[j]];
					}

					// do another round of temporal decoding (SCL)
					if (double_stage)
					{
						for (int i = 0; i < L_LSD; i++) { indexpm[i] = i; }
						////////////////////////////////
						for (int j = 0; j < Nh; j++)
							for (int k = 0; k < L; k++)
								LLR[k][j] = upLLR[decodingi][j];
						for (int k = 0; k < L; k++)  PM[k] = 0, indexpm[k] = k;
						//T_start = clock();
						decode_list(LLR, sum, (int)(log(Nh) / log(2)), A_Ach, 0, Nh, Kh, PM, L, (int)(log(L) / log(2)), LLR_in, W, SCL_index, better, worse, &Countv, &Countinfov, &Qsumunv);
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
						if (iter == softiter) PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), Nh);
						index = indexpm[0];
						//}
						//calc distance s()
						if (iter < softiter)
						{
							if (adaptive_beta)
							{
								if (softmethod == "Pyndiah") Pyndiahh_adaptivebeta(LLR, sum, oriLLRh, indexpm, iter * 2, decodingi, upLLR);
								if (softmethod == "MITSO") MITSOh(LLR, PM, sum, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR, Qsumunv);
								//	getchar();
							}
							else
							{
								if (softmethod == "Pyndiah") Pyndiahh(LLR, sum, oriLLRh, indexpm, iter * 2, decodingi, upLLR);
								if (softmethod == "MITSO") MITSOh(LLR, PM, sum, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR, Qsumunv);
								//	getchar();
							}
						}
						else
							for (int j = 0; j < Kh; j++)
								umid[decodingi][Ah[j]] = u1[index][Ah[j]];
					}
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
			int calcsum[Nv], calcu1[Nv];
			// calcsum is from the direct output of the decoder
			for (int i = 0; i < Kh; i++)
			{
				for (int j = 0; j < Nv; j++)
					calcsum[j] = umid[j][Ah[i]];
				PolarEncode_xor(calcu1, calcsum, Nv);
				for (int j = 0; j < Kv; j++)
					uout[Av[j]][Ah[i]] = calcu1[Av[j]];
			}
			// test: final output
			/*cout << "final output" << endl;
			for (int k = 0; k < Nh; k++) { cout << uout[Av[0]][k] << " "; }
			cout << endl;*/
			// test: final output
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
				printf("\r<%d> SNR= %.2f FER= %d / %d = %f BER= %d / %d = %lf",
					frame, nSNR, Num_Frame_Error, frame, (double)Num_Frame_Error / frame, Num_Error, (N_Info * frame), (double)Num_Error / N_Info / frame);
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








void system1D_LSD_MIMO(double* BER_array, double* FER_array)
{
	float design_sigma = 0;
	//freopen("result.txt", "w", stdout);
	printsetting();
	cout << "\nSNR	 " << "FER      BER" << endl;
	int Nh = Nv_1D;
	int N_Info = K - CRCL;
	int Nmax = 0;
	//assert(N_Info == Nh * Nv);
	//Interleaving pattern
	std::vector<int> intrlv_pattern(N);
	float CodeRate = (float)N_Info / N;  // N_Info/N
	Nmax = max(Nh, Nv);
	// num of symbols per vertical codeword
	int Nsym = 0;
	if (MOD_DIM == "Spatial") { Nsym = Nv / ModType; }
	else if (MOD_DIM == "Temporal") { Nsym = Nh / ModType; }
	int symbol_length = ModType;

	// Size of a mimo block is given by height(spatial) * length(temporal)
	int mimo_blk_height = 0;
	int mimo_blk_length = 0;
	if (MOD_DIM == "Spatial")
	{
		mimo_blk_height = Nsym;
		mimo_blk_length = Nh;
	}
	else if (MOD_DIM == "Temporal")
	{
		mimo_blk_height = Nv;
		mimo_blk_length = Nsym;
	}
	// Initialization
	unsigned char* a_data = new unsigned char[K];// Information CodeWord
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

	// Arrays for MIMO
	int** x_mmse_out = new int* [Nv]; // delete!
	for (int i = 0; i < Nv; i++) { x_mmse_out[i] = new int[Nh]; }

	//std::complex<float>** x_mod = new std::complex<float>* [Nsym]; // modulated symbols (spatial modulation) DELETE!
	//for (int i = 0; i < Nsym; i++) { x_mod[i] = new std::complex<float>[Nh]; }
	//std::complex<float>** x_mod_temporal = new std::complex<float>*[Nv]; // modulated symbols (temporal modulation) DELETE!
	//for (int i = 0; i < Nv; i++) x_mod_temporal[i] = new std::complex<float>[Nsym];
	//std::complex<float>** y_mod = new std::complex<float>* [Nr]; // Output symbols after AWGN channels // delete!
	//for (int i = 0; i < Nr; i++) { y_mod[i] = new std::complex<float>[Nh]; }
	//float** y_mod_real = new float* [2 * Nr];
	//for (int i = 0; i < 2*Nr; i++) y_mod_real[i] = new float[Nh];
	//float** y = new float* [Nv];// Output Signal After AWGN Channel
	//for (int i = 0; i < Nv; i++) { y[i] = new float[Nh]; }
	//float** y_mmse = new float* [2*Nsym];
	//for (int i = 0; i < 2*Nsym; i++) { y_mmse[i] = new float[Nh]; }  // MMSE detector output // delete!
	//float* y_mmse_tmp = new float[2*Nsym];
	//std::complex<float>* x_mod_tmp = new std::complex<float>[Nsym]; // modulated symbols of a single column // delete!
	//std::complex<float>* y_mod_tmp = new std::complex<float>[Nr]; // received symbols of a single column // delete!
	//float* y_mod_real_tmp = new float[2*Nr];
	//// float* y_mod_tmp = new float[Nr];  // delete!


	std::complex<float>** x_mod = new std::complex<float>*[mimo_blk_height]; // modulated symbols (spatial modulation) DELETE!
	for (int i = 0; i < mimo_blk_height; i++) { x_mod[i] = new std::complex<float>[mimo_blk_length]; }
	std::complex<float>** y_mod = new std::complex<float>*[Nr]; // Output symbols after AWGN channels // delete!
	for (int i = 0; i < Nr; i++) { y_mod[i] = new std::complex<float>[mimo_blk_length]; }
	float** y_mod_real = new float* [2 * Nr];
	for (int i = 0; i < 2 * Nr; i++) y_mod_real[i] = new float[mimo_blk_length];
	float** y = new float* [Nv];// Output Signal After AWGN Channel
	for (int i = 0; i < Nv; i++) { y[i] = new float[Nh]; }
	float** y_mmse = new float* [2 * mimo_blk_height];
	for (int i = 0; i < 2 * mimo_blk_height; i++) { y_mmse[i] = new float[mimo_blk_length]; }  // MMSE detector output // delete!
	float* y_mmse_tmp = new float[2 * mimo_blk_height];
	std::complex<float>* x_mod_tmp = new std::complex<float>[mimo_blk_height]; // modulated symbols of a single column // delete!
	std::complex<float>* y_mod_tmp = new std::complex<float>[Nr]; // received symbols of a single column // delete!
	float* y_mod_real_tmp = new float[2 * Nr];
	// float* y_mod_tmp = new float[Nr];  // delete!
	float** y_mod_real_temp_omp = new float* [mimo_blk_length]; // new delete!
	for (int i = 0; i < mimo_blk_length; i++) { y_mod_real_temp_omp[i] = new float[2 * Nr]; }
	float** y_mmse_temp_omp = new float* [mimo_blk_length]; // new delete!
	for (int i = 0; i < mimo_blk_length; i++) { y_mmse_temp_omp[i] = new float[2 * mimo_blk_height]; }


	int* x_tmp = new int[Nv]; // delete!

	// CRC bits
	unsigned char* u_crc = new unsigned char[K];

	float** upLLR = new float* [Nv];
	for (int i = 0; i < Nv; i++) { upLLR[i] = new _MM_ALIGN16 float[Nh]; }
	//float** oriLLRv = new float* [Nh];
	//for (int i = 0; i < Nh; i++) { oriLLRv[i] = new _MM_ALIGN16 float[2 * Nv + 2]; }
	float** oriLLRh = new float* [Nv];
	for (int i = 0; i < Nv; i++) { oriLLRh[i] = new _MM_ALIGN16 float[Nh]; }
	int len_ori_llr_tmp = 0;
	if (MOD_DIM == "Spatial") { len_ori_llr_tmp = Nv; }
	else if (MOD_DIM == "Temporal") { len_ori_llr_tmp = ModType * Nt; }
	float* oriLLR_tmp = new _MM_ALIGN16 float[len_ori_llr_tmp];  // delete!


	//*****************************************************************

	int** GGh = new int* [Nh];
	int** GGv = new int* [Nv];
	int** GG = new int* [N];


	MatrixXf Gh = Fn(Nh);  // Generate Matrix of horizontal coding
	MatrixXf Gv = Fn(Nv);  // Generate Matrix of vertical coding
	MatrixXf G = Fn(N);
	MatrixXcf H(Nr, Nt); // Complex channel Matrix H: Nr x Nt
	//Matrix<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic> H = MatrixXcf::Random(Nr, Nt); // Real domain channel matrix
	//// //MatrixX Definitions
	MatrixXf H_real(2 * Nr, 2 * Nt); // Real domain channel matrix
	MatrixXf received_signal(2 * Nr, 1);
	MatrixXf HTH(2 * Nt, 2 * Nt);
	MatrixXf Identity = MatrixXf::Identity(2 * Nt, 2 * Nt);
	MatrixXf output_signal(2 * Nt, 1);
	MatrixXf mmse_matrix(2 * Nt, 2 * Nr);

	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> H_real = MatrixXf::Random(2 * Nr, 2 * Nt); // Real domain channel matrix
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> received_signal = MatrixXf::Random(2 * Nr,1);
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> HTH = MatrixXf::Random(2 * Nt, 2 * Nt);
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> Identity = MatrixXf::Random(2 * Nt, 2 * Nt);
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> output_signal = MatrixXf::Random(2 * Nt, 1);
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> mmse_matrix = MatrixXf::Random(2 * Nt, 2 * Nr);


	// Matrix Definitions
	//Eigen::Matrix<float, 2 * Nr, 2 * Nt> H_real; // ʵ���ŵ�����
	//Eigen::Matrix<float, 2 * Nr, 1> received_signal; // �����źž���
	//Eigen::Matrix<float, 2 * Nt, 2 * Nt> HTH; // HTH����
	//Eigen::Matrix<float, 2 * Nt, 2 * Nt> Identity = Eigen::Matrix<float, 2 * Nt, 2 * Nt>::Identity(); // ��λ����
	//Eigen::Matrix<float, 2 * Nt, 1> output_signal; // ����źž���
	//Eigen::Matrix<float, 2 * Nt, 2 * Nr> mmse_matrix; // MMSE����

	for (int i = 0; i < N; i++) GG[i] = new int[N];
	for (int i = 0; i < Nh; i++) GGh[i] = new int[Nh];
	for (int i = 0; i < Nv; i++) GGv[i] = new int[Nv];
	for (int i = 0; i < Nh; i++)
		for (int j = 0; j < Nh; j++)
			GGh[i][j] = Gh(i, j);
	for (int i = 0; i < Nv; i++)
		for (int j = 0; j < Nv; j++)
			GGv[i][j] = Gv(i, j);
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			GG[i][j] = G(i, j);
	//*****************************************************************
	vector<int> temph(Nh);
	vector<int> temp_1D(N);
	vector<int> tempv(Nv);
	int* Ah = new int[Kh]; // Horizontal Info Bits
	//int* Ach = new int[Nh - Kh]; // Horizontal Frozen Bits
	int* A_Ach = new int[Nh];
	int* Av = new int[Kv];
	int* Acv = new int[Nv - Kv];
	int* A_Acv = new int[Nv];

	// delete!
	int* A = new int[K];
	int* Ac = new int[N - K];
	int* A_Ac = new int[N];
	int* a_encode = new int[N];
	int* a_encode_tmp = new int[N];
	int* u_1D = new int[N];
	int* x_1D = new int[N];
	int* x_intrlv_1D = new int[N];
	int* u_decod_1D = new int[N];
	float* oriLLR_1D = new float[N];
	//for SCL Decoding
	float* LLR_in = new float[L];
	float* W = new float[2 * L];
	int* SCL_index = new int[L];
	int* better = new int[L];
	int* worse = new int[L];
	int* indexpm = new int[L];
	float* PM = new float[L];
	float** LLR = new float* [L];
	for (int i = 0; i < L; i++) { LLR[i] = new _MM_ALIGN16 float[2 * N + 2]; }
	int** sum = new int* [L];
	for (int i = 0; i < L; i++) { sum[i] = new _MM_ALIGN16 int[N]; }
	int** u1 = new int* [L];
	for (int i = 0; i < L; i++) { u1[i] = new _MM_ALIGN16 int[N]; }
	float** mu_mmse = new float*[mimo_blk_length];
	for (int i = 0; i < mimo_blk_length; i++) mu_mmse[i] = new float[2*Nt];

	// delete!

	//*****************************************************************
		// Main Function
	// srand((unsigned int)time(0));
	std::random_device rd;              // ��ȡ���������
	std::mt19937 gen(rd());             // ʹ�����ӳ�ʼ��mt19937������
	std::uniform_int_distribution<> distrib(1, 1000000);
	int Num_Error, Num_Frame_Error;
	float Decoding_Duration;

	int Num_Error_new, Num_Frame_Error_new;
	float Decoding_Duration_new;

	clock_t T_start, T_finish;
	int SNR_count = 0;
	int anssc, anssd;
	float design_snr = 11;
	if (atten == 2)
	{
		float tmp = 0;
		tmp = pow(10, design_snr / 10) * symbol_length * CodeRate * Nt;
		// tmp = pow(10, nSNR / 10) * symbol_length  * Nt;
		design_sigma = (Nt * Nr) / tmp;
	}
	for (double nSNR = Min_SNR; nSNR <= Max_SNR; nSNR += SNR_step)
	{
		// TEST: det err
		int n_det_err = 0; // n error bits at detector output
		int n_det_err_mrb = 0; // error bits at detector output(MRBs)
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
			//cout << tmp << endl;
			// tmp = pow(10, nSNR / 10) * symbol_length  * Nt;
			sigma_sq = (Nt * Nr) / tmp;
			//cout << sigma_sq << endl;
		}
		/*	polar_codeconstruction(polarconh, Nh, Kh, GGh, (float)sqrt(sigma_sq), temph);
			polar_codeconstruction(polarconv, Nv, Kv, GGv, (float)sqrt(sigma_sq), tempv);*/
		polar_codeconstruction(polarcon, N, K, GG, (float)sqrt(sigma_sq), temp_1D);

		// sort temph and tempv in ascending order
		for (int i = 0; i < K - 1; i++)
			for (int j = i + 1; j < K; j++)
				if (temp_1D[i] > temp_1D[j])
				{
					int ttt = temp_1D[i];
					temp_1D[i] = temp_1D[j];
					temp_1D[j] = ttt;
				}
				// Information set and frozen set
		for (int i = 0; i < K; i++) // Information Set
		{
			A[i] = temp_1D[i];
			A_Ac[A[i]] = 1;
		}
		//save_array_txt(A_Ac, N, "polarconstruction2048_0225.txt");
		// test A_Ac
		/*for (int i = 0; i < 1024; i++)
		{
			cout << A_Ac[i];
		}*/
		// test A_Ac
		//	getchar();
		for (int i = 0; i < N - K; i++) // Frozen Bit
		{
			Ac[i] = temp_1D[K + i];
			//cout << Ac[i] << endl;
			A_Ac[Ac[i]] = 0;
		}

		save_array_txt_commasep(A_Ac, N, "polarconstructionN128K95.txt");

		// test : construction from matlab
		// cout << "test for A, Ac and A_Ac" << endl;

        // polar construction form TY
		//std::ifstream ifile_construction("polarcons_tianyuN512K349.txt");
		//int cnt_tmp_A = 0;
		//int cnt_tmp_Ac = 0;
		//for (int i = 0; i < N; i++)
		//{
		//	ifile_construction >> A_Ac[i];
		//	//cout << A_Ac[i] << endl;
		//	if (A_Ac[i] == 1)
		//	{
		//		A[cnt_tmp_A] = i;
		//		//cout << A[cnt_tmp_A] << endl;
		//		cnt_tmp_A++;
		//	}
		//	else 
		//	{
		//		Ac[cnt_tmp_Ac] = i;
		//		cnt_tmp_Ac++;
		//	}
		//}

		// for (int i = 0; i < K; i++) cout << A[i] << endl;
		// cout << "test for A, Ac and A_Ac" << endl;
		// test : construction from matlab
		// Frame Loop
		int passnum;
		int frame = 0;
		while (Num_Frame_Error < errframe)
		{
			/*if (frame == 10) {
				cout << "AAAAAAAAAA" << endl;
			}*/
			frame++;
			// generate interleaving pattern
			if (flag_interleave) { generateIntrlvPattern(intrlv_pattern, N); }

			//cout << frame << endl;
			int u_input[16] = { 0,     1     ,0,     1,
									0,     1,     0,     0,
									0 ,    1,     0,     0,
									1 ,    0,     0,     0 };
			for (int i = 0; i < N_Info; i++)
			{
				a_data[i] = (unsigned char)(distrib(gen) % 2);
				//		a_data[i] = u_input[i];
			}
			
			//cout << endl;
			if (CRCL) tx_append_crc(a_data, N_Info, CRCL);

			// test: a_data
			/*for (int i = K-1; i < K; i++) {
				cout << a_data[i] << endl;
			}*/
			// test: a_data
			for (int i = 0; i < Nv; i++)
				for (int j = 0; j < Nh; j++)
					u[i][j] = 0;
			for (int i = 0; i < N; i++) { u_1D[i] = 0; }
			for (int i = 0; i < K; i++)
			{
				u_1D[A[i]] = a_data[i];
			}
			PolarEncode_xor(x_1D, u_1D, N);
				// test:vertical encoding
			for (int i = 0; i < N; i++) x_intrlv_1D[i] = x_1D[i];
			if (flag_interleave) {
				randIntrlv(intrlv_pattern, x_intrlv_1D, N);
			}
			//for (int j = 0; j < Nh; j++)
			//{
			//	for (int i = 0; i < Nv; i++)
			//	{
			//		oriLLR_1D[j * Nv + i] = oriLLRh[i][j];
			//		// test: oriLLR_1D
			//		// cout << oriLLR_1D[j * Nv + i] << endl;
			//		// test: oriLLR_1D

			//	}
			//}
			for (int i = 0; i < Nv; i++)
			{
				for (int j = 0; j < Nh; j++) {
					x[i][j] = x_intrlv_1D[j * Nv + i];
				}
			}
			if (atten == 2)
			{
				
				// Temporal modulation
			 if (MOD_DIM == "Temporal")
				{
					//redefine x_mod (Nv(Nt)*Nsym), x_mod_tmp(Nv*1), y_mod_tmp(Nv*1), y_mod(Nv * Nsym(Nh/ModType))
					// TEST: transmitted bits:
					/*cout << "Transmitted Bits" << endl;
					for (int i = 0; i < Nv; i++)
					{
						cout << x[i][0] << x[i][1] << x[i][2] << x[i][3] << endl;
					}*/
					//modulation
					for (int i = 0; i < Nv; i++) // every row
					{
						modulation(x[i], x_mod[i], ModType, mimo_blk_length);
						//cout << mimo_blk_height << endl;
					}
					//AWGN channel (received y_mod)
					H = generateChannelMatrix(Nr, Nt);
					for (int j = 0; j < Nsym; j++)
					{
						// Extracting a single column
						for (int i = 0; i < Nv; i++) { x_mod_tmp[i] = x_mod[i][j]; }
						AWGN(Nr, Nt, x_mod_tmp, y_mod_tmp, H, sigma_sq);
						for (int i = 0; i < Nr; i++) { y_mod[i][j] = y_mod_tmp[i]; }
					}
					// Test Modulation 
					/*cout << "Modulation Symbols" << endl;
					for (int i = 0; i < Nv; i++)
					{
						cout << x_mod[i][0] << endl;
					}*/
					// Test Modulation
					//Convert H and y to Real domain
					convertHToReal(H, H_real);
					for (int j = 0; j < Nsym; j++)
					{
						for (int i = 0; i < Nr; i++) { y_mod_tmp[i] = y_mod[i][j]; }
						convertYToReal(Nr, y_mod_tmp, y_mod_real_tmp);
						for (int i = 0; i < 2 * Nr; i++) { y_mod_real[i][j] = y_mod_real_tmp[i]; }
					}
					//MMSE detection done in paralell
					for (int i = 0; i < Nsym; i++)
					{
						for (int j = 0; j < 2 * Nr; j++) { y_mod_real_temp_omp[i][j] = y_mod_real[j][i]; }
					}
					# pragma omp parallel for
					for (int j = 0; j < Nsym; j++)
					{
						MMSEdetection(H_real, Identity, received_signal, mmse_matrix, output_signal, HTH, y_mod_real_temp_omp[j], y_mmse_temp_omp[j], sigma_sq);
						//MMSEdetection_exact(H_real, Identity, received_signal, mmse_matrix, output_signal, HTH, y_mod_real_temp_omp[j], y_mmse_temp_omp[j], mu_mmse[j], sigma_sq);
					}
					for (int i = 0; i < Nsym; i++)
					{
						for (int j = 0; j < 2 * Nt; j++) { y_mmse[j][i] = y_mmse_temp_omp[i][j]; }
					}
					// Test Detection
					/*cout << "Detected Symbols" << endl;
					for (int i = 0; i < 2*Nv; i++)
					{
						cout << y_mmse[i][0] << endl;
					}			*/
					// Test Detection
					//To bit domain output
					for (int j = 0; j < Nsym; j++)
					{
						// if (i == 0) { cout << "MMSE output LLRs"; }
						/*if (i == 0)
						{
							cout << "Detection Result" << endl;
						}*/
						/*if (j == 1)
						{
							cout << n_det_err << endl;
						}*/
						for (int i = 0; i < 2 * Nt; i++) { y_mmse_tmp[i] = y_mmse[i][j]; }
						// With Spatial Modulation, Nt = Nv
						// llr_mmse_bit_exact(Nt, mu_mmse[j], y_mmse_tmp, oriLLR_tmp, sigma_sq, ModType);
						llr_mmse_bit(Nt, y_mmse_tmp, oriLLR_tmp, sigma_sq, ModType);
						int i_bit = 0;
						int bit_col = 0;
						for (int j_inner = 0; j_inner < Nv; j_inner++)
						{
							for (int i_bit_col = 0; i_bit_col < ModType; i_bit_col++)
							{
								i_bit = ModType * j_inner + i_bit_col;
								bit_col = j * ModType + i_bit_col;
								oriLLRh[j_inner][bit_col] = oriLLR_tmp[i_bit];
								upLLR[j_inner][bit_col] = oriLLR_tmp[i_bit];
							}

						}

					}
				}
			}
			/*std::default_random_engine e;
			std::normal_distribution<double> norm(0, 1);*/
			if (atten == 1)
			{
				for (int i = 0; i < Nv; i++)
					for (int j = 0; j < Nh; j++)
					{
						// 0-1; 1-(-1)
						y[i][j] = (1 - 2 * x[i][j]) + (float)sqrt(sigma_sq) * sampleNormal();
						oriLLRh[i][j] = 2 * y[i][j] / sigma_sq;
						//oriLLRv[j][i] = oriLLRh[i][j];
						upLLR[i][j] = oriLLRh[i][j];
					}
				float noise_sum = 0;
				//ifstream ifile_llrmatlab("llr_matlab.txt");
				//for (int i = 0; i < N; i++)
				//{
				//	float norm_sample = sampleNormal();
				//	float y = (1 - 2 * x_1D[i]) + (float)sqrt(sigma_sq) * norm_sample;
				//	oriLLR_1D[i] = 2 * y / sigma_sq;
				//	noise_sum += norm_sample;//pow(abs(norm_sample),2);
				//	//ifile_llrmatlab >> oriLLR_1D[i];
				//}
				//cout << "noise_sum" << noise_sum << endl;
			}
			if (true)
			{
				for (int j = 0; j < Nh; j++)
				{
					for (int i = 0; i < Nv; i++)
					{
						oriLLR_1D[j * Nv + i] = oriLLRh[i][j];
						// test: oriLLR_1D
						// cout << oriLLR_1D[j * Nv + i] << endl;
						// test: oriLLR_1D
					}
				}
			}
			if (flag_interleave) { randDeintrlv(intrlv_pattern, oriLLR_1D, N); }

			// test: LLR output
			/*std::ofstream ofile_llr("llrout_1D_1.txt");
			for (int m = 0; m < N; m++) ofile_llr << oriLLR_1D[m]<<"  ";
			ofile_llr.close();
			std::ofstream ofile_bits("info_bits_u_1.txt");
			for (int m = 0; m < K; m++) ofile_bits << u_1D[A[m]]<<"   ";
			std::ofstream ofile_encoded_bits("encoded_bits_x_1.txt");
			for (int m = 0; m < N; m++) ofile_encoded_bits << x_1D[m] << "   ";
			ofile_bits.close();
			ofile_encoded_bits.close();*/
			// test: LLR output
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < L; j++)
				{
					LLR[j][i] = oriLLR_1D[i];
					//LLR[j][i] = -y[i];
				//	LLR[j][i] = LLR_input[i];
				}
				//u[i] = u_input[i];
			}
			 // test: LLR
			//for (int i = 0; i < N; i++)
			//{
			//	for (int j = 0; j < L; j++)
			//	{
			//		cout << LLR[j][i] << "   ";
			//	}
			//	//u[i] = u_input[i];
			//}
			// test: LLR
			//Polar Decoding
			/*

			Add Contents here!
			*/
			bool markcrc = false;
			for (int k = 0; k < L; k++)  PM[k] = 0;
			for (int i = 0; i < L; i++)
				indexpm[i] = i;
			//T_start = clock();
			float QQQ = 0;
			int index = 0, C1D = 0, Count_info = 0;
			decode_list(LLR, sum, (int)(log(N) / log(2)), A_Ac, 0, N, K, PM, L, (int)(log(L) / log(2)), LLR_in, W, SCL_index, better, worse, &C1D, &Count_info, &QQQ);
			// T_finish = clock();
			// Decoding_Duration += (T_finish - T_start);
			for (int i = 0; i < L - 1; i++)
				for (int j = i + 1; j < L; j++)
					if (PM[indexpm[i]] < PM[indexpm[j]])
					{
						int ttt = indexpm[i];
						indexpm[i] = indexpm[j];
						indexpm[j] = ttt;
					}
			// test:PM
			/*for (int i = 0; i < L; i++)
				cout << PM[i] << endl;*/
				// test:PM
			PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), N);
			index = indexpm[0]; 
			// test: decoder output
			/*std::ofstream ofile_dec_out("dec_out_u1.txt");
			for (int i = 0; i < K; i++) { ofile_dec_out << u1[index][A[i]] << "   "; }
			ofile_dec_out.close();*/
			// test: decoder output
			if (CRCL > 0) {
				int A_cnt = 0;
				for (int i = 0; i < N; i++) {
					if (A_Ac[i] == 1) { u_crc[A_cnt] = (unsigned char)u1[index][i]; A_cnt++; }
				}
				if (rx_check_crc(u_crc, K, CRCL) == false) {
					int cnt = 1;
					while (cnt < L) {
						PolarEncode_xor(u1[indexpm[cnt]], sum[indexpm[cnt]], N);
						A_cnt = 0;
						for (int i = 0; i < N; i++) {if (A_Ac[i] == 1) { u_crc[A_cnt] = (unsigned char)u1[indexpm[cnt]][i]; A_cnt++; }}
						//cout << rx_check_crc(u_crc, K, CRCL) << endl;
						if (rx_check_crc(u_crc, K, CRCL) == true) {
							index = indexpm[cnt]; break;
						}
						cnt++;
					}
				}
			}
			int calcsum[Nv], calcu1[Nv];
			// int* calcu1h=new int[Nh];
			// calcsum is from the direct output of the decoder
			// if output is the result of row decoding (meaning there is another column decoding to be done)

			// test: final output
			/*cout << "final output" << endl;
			for (int k = 0; k < Nh; k++) { cout << uout[Av[0]][k] << " "; }
			cout << endl;*/
			// test: final output
			//now without CRC
			int Temp_Error = Num_Error;
			for (int i = 0; i < K - CRCL; i++)
				Num_Error += abs(u1[index][A[i]] - u_1D[A[i]]);
			if (Temp_Error < Num_Error) Num_Frame_Error++;
			if (frame % 10 == 0)
			//if (true)
			{
				/*printf("\r<%d> SNR= %.2f FER= %d / %d = %f BER= %d / %d = %lf E-Det = %f, E-MRB = %f"  ,
					frame, nSNR, Num_Frame_Error, frame, (double)Num_Frame_Error / frame, Num_Error, (N_Info * frame), (double)Num_Error / N_Info / frame,
					(double)n_det_err/total_bit_cnt, (double)n_det_err_mrb/(total_bit_cnt/(ModType/2)));*/
				printf("\r<%d> SNR= %.2f FER= %d / %d = %f BER= %d / %d = %lf ",
					frame, nSNR, Num_Frame_Error, frame, (double)Num_Frame_Error / frame, Num_Error, (N_Info * frame), (double)Num_Error / N_Info / frame);
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
		// save to .txt file
		string polarde_serial = "";
		for (int i = 0; i < 1; i++)
			polarde_serial = polarde_serial + polarde_array[i] + string("_");
		string file_name = string("ResFinal\\1D_LSD_D_") + polarde_serial + string("T") + to_string(Nt) + string("R") + to_string(Nr)
			+ string("Q") + to_string(int(pow(2, ModType))) + string("L") + to_string(N) + string("K") + to_string(K)
			+ string("T") + to_string(softiter) + string("H") + to_string(half_iter) + string("UB") + to_string(usebeta)
			+ string(".txt");
		save_array_txt(nSNR, BER_array[SNR_count], FER_array[SNR_count], file_name);
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
	//delete[] Ach;
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
	for (int i = 0; i < N; i++)
		delete[] GG[i];
	delete[] GG;
	for (int i = 0; i < mimo_blk_height; i++) {
		delete[] x_mod[i];
	}
	delete[] x_mod;
	for (int i = 0; i < Nr; i++) {
		delete[] y_mod[i];
	}
	delete[] y_mod;
	for (int i = 0; i < 2 * mimo_blk_height; i++) {
		delete[] y_mmse[i];
	}
	delete[] y_mmse;
	for (int i = 0; i < 2 * Nr; i++) {
		delete[] y_mod_real[i];
	}
	for (int i = 0; i < mimo_blk_length; i++)
	{
		delete[] y_mod_real_temp_omp[i];
	}
	delete[] y_mod_real_temp_omp;
	for (int i = 0; i < mimo_blk_length; i++)
	{
		delete[] y_mmse_temp_omp[i];
	}
	delete[] y_mmse_temp_omp;
	delete[] y_mod_real;
	delete[] x_mod_tmp;
	delete[] y_mod_tmp;
	delete[] y_mod_real_tmp;
	delete[] x_tmp;
	delete[] oriLLR_tmp;

	delete[] A;
	delete[] Ac;
	delete[] A_Ac;
	delete[] a_encode;
	delete[] a_encode_tmp;
	delete[] u_1D;
	delete[] x_1D;
	delete[] u_decod_1D;
	delete[] oriLLR_1D;
	delete[] x_intrlv_1D;
	for (int i = 0; i < L; i++)
		delete[]LLR[i];
	delete[] LLR;
	for (int i = 0; i < L; i++)
		delete[]sum[i];
	delete[] sum;
	for (int i = 0; i < L; i++)
		delete[]u1[i];
	delete[] u1;
	delete[] PM;
	delete[] LLR_in;
	delete[] W;
	delete[] SCL_index;
	delete[] better;
	delete[] worse;
	for (int i = 0; i < mimo_blk_length; i++)
		delete[] mu_mmse[i];
	delete[] mu_mmse;
}



void system1D_LSD_MIMO_SVD(double* BER_array, double* FER_array)
{
	float design_sigma = 0;
	//freopen("result.txt", "w", stdout);
	printsetting();
	cout << "\nSNR	 " << "FER      BER" << endl;
	int Nh = Nv_1D;
	int N_Info = K - CRCL;
	int Nmax = 0;
	//assert(N_Info == Nh * Nv);
	//Interleaving pattern
	std::vector<int> intrlv_pattern(N);
	float CodeRate = (float)N_Info / N;  // N_Info/N
	Nmax = max(Nh, Nv);
	// num of symbols per vertical codeword
	int Nsym = 0;
	if (MOD_DIM == "Spatial") { Nsym = Nv / ModType; }
	else if (MOD_DIM == "Temporal") { Nsym = Nh / ModType; }
	int symbol_length = ModType;

	// Size of a mimo block is given by height(spatial) * length(temporal)
	int mimo_blk_height = 0;
	int mimo_blk_length = 0;
	if (MOD_DIM == "Spatial")
	{
		mimo_blk_height = Nsym;
		mimo_blk_length = Nh;
	}
	else if (MOD_DIM == "Temporal")
	{
		mimo_blk_height = Nv;
		mimo_blk_length = Nsym;
	}
	// Initialization
	unsigned char* a_data = new unsigned char[K];// Information CodeWord
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

	// Arrays for MIMO
	int** x_mmse_out = new int* [Nv]; // delete!
	for (int i = 0; i < Nv; i++) { x_mmse_out[i] = new int[Nh]; }

	//std::complex<float>** x_mod = new std::complex<float>* [Nsym]; // modulated symbols (spatial modulation) DELETE!
	//for (int i = 0; i < Nsym; i++) { x_mod[i] = new std::complex<float>[Nh]; }
	//std::complex<float>** x_mod_temporal = new std::complex<float>*[Nv]; // modulated symbols (temporal modulation) DELETE!
	//for (int i = 0; i < Nv; i++) x_mod_temporal[i] = new std::complex<float>[Nsym];
	//std::complex<float>** y_mod = new std::complex<float>* [Nr]; // Output symbols after AWGN channels // delete!
	//for (int i = 0; i < Nr; i++) { y_mod[i] = new std::complex<float>[Nh]; }
	//float** y_mod_real = new float* [2 * Nr];
	//for (int i = 0; i < 2*Nr; i++) y_mod_real[i] = new float[Nh];
	//float** y = new float* [Nv];// Output Signal After AWGN Channel
	//for (int i = 0; i < Nv; i++) { y[i] = new float[Nh]; }
	//float** y_mmse = new float* [2*Nsym];
	//for (int i = 0; i < 2*Nsym; i++) { y_mmse[i] = new float[Nh]; }  // MMSE detector output // delete!
	//float* y_mmse_tmp = new float[2*Nsym];
	//std::complex<float>* x_mod_tmp = new std::complex<float>[Nsym]; // modulated symbols of a single column // delete!
	//std::complex<float>* y_mod_tmp = new std::complex<float>[Nr]; // received symbols of a single column // delete!
	//float* y_mod_real_tmp = new float[2*Nr];
	//// float* y_mod_tmp = new float[Nr];  // delete!


	std::complex<float>** x_mod = new std::complex<float>*[mimo_blk_height]; // modulated symbols (spatial modulation) DELETE!
	for (int i = 0; i < mimo_blk_height; i++) { x_mod[i] = new std::complex<float>[mimo_blk_length]; }
	std::complex<float>** y_mod = new std::complex<float>*[Nr]; // Output symbols after AWGN channels // delete!
	for (int i = 0; i < Nr; i++) { y_mod[i] = new std::complex<float>[mimo_blk_length]; }
	float** y_mod_real = new float* [2 * Nr];
	for (int i = 0; i < 2 * Nr; i++) y_mod_real[i] = new float[mimo_blk_length];
	float** y = new float* [Nv];// Output Signal After AWGN Channel
	for (int i = 0; i < Nv; i++) { y[i] = new float[Nh]; }
	float** y_mmse = new float* [2 * mimo_blk_height];
	for (int i = 0; i < 2 * mimo_blk_height; i++) { y_mmse[i] = new float[mimo_blk_length]; }  // MMSE detector output // delete!
	float* y_mmse_tmp = new float[2 * mimo_blk_height];
	std::complex<float>* x_mod_tmp = new std::complex<float>[mimo_blk_height]; // modulated symbols of a single column // delete!
	std::complex<float>* y_mod_tmp = new std::complex<float>[Nr]; // received symbols of a single column // delete!
	float* y_mod_real_tmp = new float[2 * Nr];
	// float* y_mod_tmp = new float[Nr];  // delete!
	float** y_mod_real_temp_omp = new float* [mimo_blk_length]; // new delete!
	for (int i = 0; i < mimo_blk_length; i++) { y_mod_real_temp_omp[i] = new float[2 * Nr]; }
	float** y_mmse_temp_omp = new float* [mimo_blk_length]; // new delete!
	for (int i = 0; i < mimo_blk_length; i++) { y_mmse_temp_omp[i] = new float[2 * mimo_blk_height]; }


	int* x_tmp = new int[Nv]; // delete!

	// CRC bits
	unsigned char* u_crc = new unsigned char[K];

	float** upLLR = new float* [Nv];
	for (int i = 0; i < Nv; i++) { upLLR[i] = new _MM_ALIGN16 float[Nh]; }
	//float** oriLLRv = new float* [Nh];
	//for (int i = 0; i < Nh; i++) { oriLLRv[i] = new _MM_ALIGN16 float[2 * Nv + 2]; }
	float** oriLLRh = new float* [Nv];
	for (int i = 0; i < Nv; i++) { oriLLRh[i] = new _MM_ALIGN16 float[Nh]; }
	int len_ori_llr_tmp = 0;
	if (MOD_DIM == "Spatial") { len_ori_llr_tmp = Nv; }
	else if (MOD_DIM == "Temporal") { len_ori_llr_tmp = ModType * Nt; }
	float* oriLLR_tmp = new _MM_ALIGN16 float[len_ori_llr_tmp];  // delete!

	// for SVD precoding
	float* singular_values = new float[Nt];
	float* antenna_power = new float[Nt];
	std::complex<float>* x_precoded_tmp = new std::complex<float>[mimo_blk_height]; // precoded symbols of a single column // delete!
	std::complex<float>* y_precoded_tmp = new std::complex<float>[Nr];


	//*****************************************************************

	int** GGh = new int* [Nh];
	int** GGv = new int* [Nv];
	int** GG = new int* [N];


	MatrixXf Gh = Fn(Nh);  // Generate Matrix of horizontal coding
	MatrixXf Gv = Fn(Nv);  // Generate Matrix of vertical coding
	MatrixXf G = Fn(N);
	MatrixXcf H(Nr, Nt); // Complex channel Matrix H: Nr x Nt
	//Matrix<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic> H = MatrixXcf::Random(Nr, Nt); // Real domain channel matrix
	//// //MatrixX Definitions
	MatrixXf H_real(2 * Nr, 2 * Nt); // Real domain channel matrix
	MatrixXf received_signal(2 * Nr, 1);
	MatrixXf HTH(2 * Nt, 2 * Nt);
	MatrixXf Identity = MatrixXf::Identity(2 * Nt, 2 * Nt);
	MatrixXf output_signal(2 * Nt, 1);
	MatrixXf mmse_matrix(2 * Nt, 2 * Nr);

	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> H_real = MatrixXf::Random(2 * Nr, 2 * Nt); // Real domain channel matrix
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> received_signal = MatrixXf::Random(2 * Nr,1);
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> HTH = MatrixXf::Random(2 * Nt, 2 * Nt);
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> Identity = MatrixXf::Random(2 * Nt, 2 * Nt);
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> output_signal = MatrixXf::Random(2 * Nt, 1);
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> mmse_matrix = MatrixXf::Random(2 * Nt, 2 * Nr);


	// Matrix Definitions
	//Eigen::Matrix<float, 2 * Nr, 2 * Nt> H_real; // ʵ���ŵ�����
	//Eigen::Matrix<float, 2 * Nr, 1> received_signal; // �����źž���
	//Eigen::Matrix<float, 2 * Nt, 2 * Nt> HTH; // HTH����
	//Eigen::Matrix<float, 2 * Nt, 2 * Nt> Identity = Eigen::Matrix<float, 2 * Nt, 2 * Nt>::Identity(); // ��λ����
	//Eigen::Matrix<float, 2 * Nt, 1> output_signal; // ����źž���
	//Eigen::Matrix<float, 2 * Nt, 2 * Nr> mmse_matrix; // MMSE����

	for (int i = 0; i < N; i++) GG[i] = new int[N];
	for (int i = 0; i < Nh; i++) GGh[i] = new int[Nh];
	for (int i = 0; i < Nv; i++) GGv[i] = new int[Nv];
	for (int i = 0; i < Nh; i++)
		for (int j = 0; j < Nh; j++)
			GGh[i][j] = Gh(i, j);
	for (int i = 0; i < Nv; i++)
		for (int j = 0; j < Nv; j++)
			GGv[i][j] = Gv(i, j);
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			GG[i][j] = G(i, j);
	//*****************************************************************
	vector<int> temph(Nh);
	vector<int> temp_1D(N);
	vector<int> tempv(Nv);
	int* Ah = new int[Kh]; // Horizontal Info Bits
	//int* Ach = new int[Nh - Kh]; // Horizontal Frozen Bits
	int* A_Ach = new int[Nh];
	int* Av = new int[Kv];
	int* Acv = new int[Nv - Kv];
	int* A_Acv = new int[Nv];

	// delete!
	int* A = new int[K];
	int* Ac = new int[N - K];
	int* A_Ac = new int[N];
	int* a_encode = new int[N];
	int* a_encode_tmp = new int[N];
	int* u_1D = new int[N];
	int* x_1D = new int[N];
	int* x_intrlv_1D = new int[N];
	int* u_decod_1D = new int[N];
	float* oriLLR_1D = new float[N];
	//for SCL Decoding
	float* LLR_in = new float[L];
	float* W = new float[2 * L];
	int* SCL_index = new int[L];
	int* better = new int[L];
	int* worse = new int[L];
	int* indexpm = new int[L];
	float* PM = new float[L];
	float** LLR = new float* [L];
	for (int i = 0; i < L; i++) { LLR[i] = new _MM_ALIGN16 float[2 * N + 2]; }
	int** sum = new int* [L];
	for (int i = 0; i < L; i++) { sum[i] = new _MM_ALIGN16 int[N]; }
	int** u1 = new int* [L];
	for (int i = 0; i < L; i++) { u1[i] = new _MM_ALIGN16 int[N]; }
	float** mu_mmse = new float* [mimo_blk_length];
	for (int i = 0; i < mimo_blk_length; i++) mu_mmse[i] = new float[2 * Nt];

	// delete!

	//*****************************************************************
		// Main Function
	// srand((unsigned int)time(0));
	std::random_device rd;              // ��ȡ���������
	std::mt19937 gen(rd());             // ʹ�����ӳ�ʼ��mt19937������
	std::uniform_int_distribution<> distrib(1, 1000000);
	int Num_Error, Num_Frame_Error;
	float Decoding_Duration;

	int Num_Error_new, Num_Frame_Error_new;
	float Decoding_Duration_new;

	clock_t T_start, T_finish;
	int SNR_count = 0;
	int anssc, anssd;
	float design_snr = 11;
	if (atten == 2)
	{
		float tmp = 0;
		tmp = pow(10, design_snr / 10) * symbol_length * CodeRate * Nt;
		// tmp = pow(10, nSNR / 10) * symbol_length  * Nt;
		design_sigma = (Nt * Nr) / tmp;
	}
	for (double nSNR = Min_SNR; nSNR <= Max_SNR; nSNR += SNR_step)
	{
		// TEST: det err
		int n_det_err = 0; // n error bits at detector output
		int n_det_err_mrb = 0; // error bits at detector output(MRBs)
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
			//cout << tmp << endl;
			// tmp = pow(10, nSNR / 10) * symbol_length  * Nt;
			sigma_sq = (Nt * Nr) / tmp;
			//cout << sigma_sq << endl;
		}
		/*	polar_codeconstruction(polarconh, Nh, Kh, GGh, (float)sqrt(sigma_sq), temph);
			polar_codeconstruction(polarconv, Nv, Kv, GGv, (float)sqrt(sigma_sq), tempv);*/
		polar_codeconstruction(polarcon, N, K, GG, (float)sqrt(sigma_sq), temp_1D);

		// sort temph and tempv in ascending order
		for (int i = 0; i < K - 1; i++)
			for (int j = i + 1; j < K; j++)
				if (temp_1D[i] > temp_1D[j])
				{
					int ttt = temp_1D[i];
					temp_1D[i] = temp_1D[j];
					temp_1D[j] = ttt;
				}
		// Information set and frozen set
		for (int i = 0; i < K; i++) // Information Set
		{
			A[i] = temp_1D[i];
			A_Ac[A[i]] = 1;
		}
		//save_array_txt(A_Ac, N, "polarconstruction2048_0225.txt");
		// test A_Ac
		/*for (int i = 0; i < 1024; i++)
		{
			cout << A_Ac[i];
		}*/
		// test A_Ac
		//	getchar();
		for (int i = 0; i < N - K; i++) // Frozen Bit
		{
			Ac[i] = temp_1D[K + i];
			//cout << Ac[i] << endl;
			A_Ac[Ac[i]] = 0;
		}

		save_array_txt_commasep(A_Ac, N, "polarconstructionN128K95.txt");

		// test : construction from matlab
		// cout << "test for A, Ac and A_Ac" << endl;

		// polar construction form TY
		//std::ifstream ifile_construction("polarcons_tianyuN512K349.txt");
		//int cnt_tmp_A = 0;
		//int cnt_tmp_Ac = 0;
		//for (int i = 0; i < N; i++)
		//{
		//	ifile_construction >> A_Ac[i];
		//	//cout << A_Ac[i] << endl;
		//	if (A_Ac[i] == 1)
		//	{
		//		A[cnt_tmp_A] = i;
		//		//cout << A[cnt_tmp_A] << endl;
		//		cnt_tmp_A++;
		//	}
		//	else
		//	{
		//		Ac[cnt_tmp_Ac] = i;
		//		cnt_tmp_Ac++;
		//	}
		//}

		// for (int i = 0; i < K; i++) cout << A[i] << endl;
		// cout << "test for A, Ac and A_Ac" << endl;
		// test : construction from matlab
		// Frame Loop
		int passnum;
		int frame = 0;
		while (Num_Frame_Error < errframe)
		{
			/*if (frame == 10) {
				cout << "AAAAAAAAAA" << endl;
			}*/
			frame++;
			// generate interleaving pattern
			if (flag_interleave) { generateIntrlvPattern(intrlv_pattern, N); }

			//cout << frame << endl;
			int u_input[16] = { 0,     1     ,0,     1,
									0,     1,     0,     0,
									0 ,    1,     0,     0,
									1 ,    0,     0,     0 };
			for (int i = 0; i < N_Info; i++)
			{
				a_data[i] = (unsigned char)(distrib(gen) % 2);
				//		a_data[i] = u_input[i];
			}

			//cout << endl;
			if (CRCL) tx_append_crc(a_data, N_Info, CRCL);

			// test: a_data
			/*for (int i = K-1; i < K; i++) {
				cout << a_data[i] << endl;
			}*/
			// test: a_data
			for (int i = 0; i < Nv; i++)
				for (int j = 0; j < Nh; j++)
					u[i][j] = 0;
			for (int i = 0; i < N; i++) { u_1D[i] = 0; }
			for (int i = 0; i < K; i++)
			{
				u_1D[A[i]] = a_data[i];
			}
			PolarEncode_xor(x_1D, u_1D, N);
			// test:vertical encoding
			for (int i = 0; i < N; i++) x_intrlv_1D[i] = x_1D[i];
			if (flag_interleave) {
				randIntrlv(intrlv_pattern, x_intrlv_1D, N);
			}
			//for (int j = 0; j < Nh; j++)
			//{
			//	for (int i = 0; i < Nv; i++)
			//	{
			//		oriLLR_1D[j * Nv + i] = oriLLRh[i][j];
			//		// test: oriLLR_1D
			//		// cout << oriLLR_1D[j * Nv + i] << endl;
			//		// test: oriLLR_1D

			//	}
			//}
			for (int i = 0; i < Nv; i++)
			{
				for (int j = 0; j < Nh; j++) {
					x[i][j] = x_intrlv_1D[j * Nv + i];
				}
			}
			if (atten == 2)
			{

				// Temporal modulation
				if (MOD_DIM == "Temporal")
				{
					//redefine x_mod (Nv(Nt)*Nsym), x_mod_tmp(Nv*1), y_mod_tmp(Nv*1), y_mod(Nv * Nsym(Nh/ModType))
					// TEST: transmitted bits:
					/*cout << "Transmitted Bits" << endl;
					for (int i = 0; i < Nv; i++)
					{
						cout << x[i][0] << x[i][1] << x[i][2] << x[i][3] << endl;
					}*/
					//modulation
					for (int i = 0; i < Nv; i++) // every row
					{
						modulation(x[i], x_mod[i], ModType, mimo_blk_length);
						//cout << mimo_blk_height << endl;
					}
					//AWGN channel (received y_mod)
					H = generateChannelMatrix(Nr, Nt);
					JacobiSVD<MatrixXcf> svd(H, ComputeFullU | ComputeFullV);
					MatrixXcf U = svd.matrixU();
					MatrixXcf V = svd.matrixV();
					auto S = svd.singularValues();
					for (int i = 0; i < Nt; i++) { antenna_power[i] = 1.0; singular_values[i] = S[i]; }
					// SVD power allocation
					float* ev_tmp = new float[Nt];
					float sum_ev = 0;
					for (int i = 0; i < Nt; i++) {
						if (i < border) {
							ev_tmp[i] = 1 / (singular_values[i] * singular_values[i]);
						}
						else {
							ev_tmp[i] = 1 / (singular_values[i] * singular_values[i]);
						}
					}
					for (int i = 0; i < Nt; i++) {
						sum_ev += ev_tmp[i];
					}
					for (int i = 0; i < Nt; i++) {
						antenna_power[i] = Nt*ev_tmp[i] / sum_ev;
						// cout << antenna_power[i] << endl;
					}
					for (int i = 0; i < Nt; i++) { S[i] = S[i] * sqrt(antenna_power[i]); }
					for (int j = 0; j < Nsym; j++)
					{
						/*std::complex<float>* x_mod_tmp = new std::complex<float>[Nv];
						std::complex<float>* x_precoded_tmp = new std::complex<float>[Nt];
						std::complex<float>* y_mod_tmp = new std::complex<float>[Nr];
						std::complex<float>* y_precoded_tmp = new std::complex<float>[Nr];*/
						//// Extracting a single column
						for (int i = 0; i < Nv; i++) {
							x_mod_tmp[i] = x_mod[i][j] * sqrt(antenna_power[i]);
							// trans_energy += (x_mod_tmp[i].real() * x_mod_tmp[i].real() + x_mod_tmp[i].imag() * x_mod_tmp[i].imag());
						}
						precodeSVD(Nr, Nt, x_mod_tmp, x_precoded_tmp, V); // precoding: x_precode = V*x
						// for (int i = 0; i < Nt; i++) { x_precoded_tmp[i] = x_precoded_tmp[i] * sqrt(antenna_power[i]); }
						AWGN(Nr, Nt, x_precoded_tmp, y_mod_tmp, H, sigma_sq);
						precodeSVD_preprocessing(Nr, Nt, y_mod_tmp, y_precoded_tmp, U); // precoding: y_precode = U^H*y
						for (int i = 0; i < Nr; i++) { y_mod[i][j] = y_precoded_tmp[i]; }
						for (int i = 0; i < Nt; i++) { y_mod[i][j] = S[i] * y_mod[i][j] / (S[i] * S[i] + sigma_sq); } // MMSE detection with precoding
						//for (int i = 0; i < Nt; i++) { y_mod[i][j] = S[i] * UTy(i,j) / (S[i] * S[i] + sigma_sq); } // MMSE detection with precoding
						//delete[] x_mod_tmp; delete[] x_precoded_tmp; delete[] y_mod_tmp; delete[] y_precoded_tmp;
					}

					/*for (int j = 0; j < Nsym; j++) {
						for (int i = 0; i < Nv; i++) {
							std::cout << y_mod[i][j] <<"   " << x_mod[i][j] << std::endl;
						}
					}*/
					//for (int j = 0; j < Nsym; j++)
					//{
					//	// Extracting a single column
					//	for (int i = 0; i < Nv; i++) { x_mod_tmp[i] = x_mod[i][j]; }
					//	AWGN(Nr, Nt, x_mod_tmp, y_mod_tmp, H, sigma_sq);
					//	for (int i = 0; i < Nr; i++) { y_mod[i][j] = y_mod_tmp[i]; }
					//}
					// Test Modulation 
					/*cout << "Modulation Symbols" << endl;
					for (int i = 0; i < Nv; i++)
					{
						cout << x_mod[i][0] << endl;
					}*/
					// Test Modulation
					//Convert H and y to Real domain
					convertHToReal(H, H_real);
					for (int j = 0; j < Nsym; j++)
					{
						for (int i = 0; i < Nr; i++) { y_mod_tmp[i] = y_mod[i][j]; }
						convertYToReal(Nt, y_mod_tmp, y_mod_real_tmp);
						for (int i = 0; i < 2 * Nr; i++) { y_mod_real[i][j] = y_mod_real_tmp[i]; }
					}
					//MMSE detection done in paralell
					for (int i = 0; i < Nsym; i++)
					{
						for (int j = 0; j < 2 * Nr; j++) { y_mod_real_temp_omp[i][j] = y_mod_real[j][i]; }
					}
//# pragma omp parallel for
//					for (int j = 0; j < Nsym; j++)
//					{
//						MMSEdetection(H_real, Identity, received_signal, mmse_matrix, output_signal, HTH, y_mod_real_temp_omp[j], y_mmse_temp_omp[j], sigma_sq);
//						//MMSEdetection_exact(H_real, Identity, received_signal, mmse_matrix, output_signal, HTH, y_mod_real_temp_omp[j], y_mmse_temp_omp[j], mu_mmse[j], sigma_sq);
//					}
					for (int i = 0; i < Nsym; i++)
					{
						for (int j = 0; j < 2 * Nt; j++) { y_mmse[j][i] = y_mmse_temp_omp[i][j]; }
					}

					// std::cout << y_mod_real[0][0] << std::endl;
					// Test Detection
					/*cout << "Detected Symbols" << endl;
					for (int i = 0; i < 2*Nv; i++)
					{
						cout << y_mmse[i][0] << endl;
					}			*/
					// Test Detection
					//To bit domain output
					for (int j = 0; j < Nsym; j++)
					{
						// if (i == 0) { cout << "MMSE output LLRs"; }
						/*if (i == 0)
						{
							cout << "Detection Result" << endl;
						}*/
						/*if (j == 1)
						{
							cout << n_det_err << endl;
						}*/
						for (int i = 0; i < 2 * Nt; i++) { y_mmse_tmp[i] = y_mod_real[i][j]; }
						// With Spatial Modulation, Nt = Nv
						// llr_mmse_bit_exact(Nt, mu_mmse[j], y_mmse_tmp, oriLLR_tmp, sigma_sq, ModType);
						llr_mmse_bit(Nt, y_mmse_tmp, oriLLR_tmp, sigma_sq, ModType);
						int i_bit = 0;
						int bit_col = 0;
						for (int j_inner = 0; j_inner < Nv; j_inner++)
						{
							for (int i_bit_col = 0; i_bit_col < ModType; i_bit_col++)
							{
								i_bit = ModType * j_inner + i_bit_col;
								bit_col = j * ModType + i_bit_col;
								oriLLRh[j_inner][bit_col] = oriLLR_tmp[i_bit];
								upLLR[j_inner][bit_col] = oriLLR_tmp[i_bit];
							}

						}

					}
				}
			}

			/*int err_tmp = 0;
			for (int i = 0; i < Nv; i++) {
				for (int j = 0; j < Nh; j++) {
					if ((oriLLRh[i][j] > 0 && x[i][j] == 1) || (oriLLRh[i][j] < 0 && x[i][j] == 0)) { err_tmp++; }
					std::cout << oriLLRh[i][j] << "  " << x[i][j] << std::endl;
				}

			}
			std::cout << err_tmp << endl;*/


			/*std::default_random_engine e;
			std::normal_distribution<double> norm(0, 1);*/
			if (atten == 1)
			{
				for (int i = 0; i < Nv; i++)
					for (int j = 0; j < Nh; j++)
					{
						// 0-1; 1-(-1)
						y[i][j] = (1 - 2 * x[i][j]) + (float)sqrt(sigma_sq) * sampleNormal();
						oriLLRh[i][j] = 2 * y[i][j] / sigma_sq;
						//oriLLRv[j][i] = oriLLRh[i][j];
						upLLR[i][j] = oriLLRh[i][j];
					}
				float noise_sum = 0;
				//ifstream ifile_llrmatlab("llr_matlab.txt");
				//for (int i = 0; i < N; i++)
				//{
				//	float norm_sample = sampleNormal();
				//	float y = (1 - 2 * x_1D[i]) + (float)sqrt(sigma_sq) * norm_sample;
				//	oriLLR_1D[i] = 2 * y / sigma_sq;
				//	noise_sum += norm_sample;//pow(abs(norm_sample),2);
				//	//ifile_llrmatlab >> oriLLR_1D[i];
				//}
				//cout << "noise_sum" << noise_sum << endl;
			}
			if (true)
			{
				for (int j = 0; j < Nh; j++)
				{
					for (int i = 0; i < Nv; i++)
					{
						oriLLR_1D[j * Nv + i] = oriLLRh[i][j];
						// test: oriLLR_1D
						// cout << oriLLR_1D[j * Nv + i] << endl;
						// test: oriLLR_1D
					}
				}
			}
			if (flag_interleave) { randDeintrlv(intrlv_pattern, oriLLR_1D, N); }

			// test: LLR output
			/*std::ofstream ofile_llr("llrout_1D_1.txt");
			for (int m = 0; m < N; m++) ofile_llr << oriLLR_1D[m]<<"  ";
			ofile_llr.close();
			std::ofstream ofile_bits("info_bits_u_1.txt");
			for (int m = 0; m < K; m++) ofile_bits << u_1D[A[m]]<<"   ";
			std::ofstream ofile_encoded_bits("encoded_bits_x_1.txt");
			for (int m = 0; m < N; m++) ofile_encoded_bits << x_1D[m] << "   ";
			ofile_bits.close();
			ofile_encoded_bits.close();*/
			// test: LLR output
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < L; j++)
				{
					LLR[j][i] = oriLLR_1D[i];
					//LLR[j][i] = -y[i];
				//	LLR[j][i] = LLR_input[i];
				}
				//u[i] = u_input[i];
			}
			// test: LLR
		   //for (int i = 0; i < N; i++)
		   //{
		   //	for (int j = 0; j < L; j++)
		   //	{
		   //		cout << LLR[j][i] << "   ";
		   //	}
		   //	//u[i] = u_input[i];
		   //}
		   // test: LLR
		   //Polar Decoding
		   /*

		   Add Contents here!
		   */
			bool markcrc = false;
			for (int k = 0; k < L; k++)  PM[k] = 0;
			for (int i = 0; i < L; i++)
				indexpm[i] = i;
			//T_start = clock();
			float QQQ = 0;
			int index = 0, C1D = 0, Count_info = 0;
			decode_list(LLR, sum, (int)(log(N) / log(2)), A_Ac, 0, N, K, PM, L, (int)(log(L) / log(2)), LLR_in, W, SCL_index, better, worse, &C1D, &Count_info, &QQQ);
			// T_finish = clock();
			// Decoding_Duration += (T_finish - T_start);
			for (int i = 0; i < L - 1; i++)
				for (int j = i + 1; j < L; j++)
					if (PM[indexpm[i]] < PM[indexpm[j]])
					{
						int ttt = indexpm[i];
						indexpm[i] = indexpm[j];
						indexpm[j] = ttt;
					}
			// test:PM
			/*for (int i = 0; i < L; i++)
				cout << PM[i] << endl;*/
				// test:PM
			PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), N);
			index = indexpm[0];
			// test: decoder output
			/*std::ofstream ofile_dec_out("dec_out_u1.txt");
			for (int i = 0; i < K; i++) { ofile_dec_out << u1[index][A[i]] << "   "; }
			ofile_dec_out.close();*/
			// test: decoder output
			if (CRCL > 0) {
				int A_cnt = 0;
				for (int i = 0; i < N; i++) {
					if (A_Ac[i] == 1) { u_crc[A_cnt] = (unsigned char)u1[index][i]; A_cnt++; }
				}
				if (rx_check_crc(u_crc, K, CRCL) == false) {
					int cnt = 1;
					while (cnt < L) {
						PolarEncode_xor(u1[indexpm[cnt]], sum[indexpm[cnt]], N);
						A_cnt = 0;
						for (int i = 0; i < N; i++) { if (A_Ac[i] == 1) { u_crc[A_cnt] = (unsigned char)u1[indexpm[cnt]][i]; A_cnt++; } }
						//cout << rx_check_crc(u_crc, K, CRCL) << endl;
						if (rx_check_crc(u_crc, K, CRCL) == true) {
							index = indexpm[cnt]; break;
						}
						cnt++;
					}
				}
			}
			int calcsum[Nv], calcu1[Nv];
			// int* calcu1h=new int[Nh];
			// calcsum is from the direct output of the decoder
			// if output is the result of row decoding (meaning there is another column decoding to be done)

			// test: final output
			/*cout << "final output" << endl;
			for (int k = 0; k < Nh; k++) { cout << uout[Av[0]][k] << " "; }
			cout << endl;*/
			// test: final output
			//now without CRC
			int Temp_Error = Num_Error;
			for (int i = 0; i < K - CRCL; i++)
				Num_Error += abs(u1[index][A[i]] - u_1D[A[i]]);
			if (Temp_Error < Num_Error) Num_Frame_Error++;
			if (frame % 10 == 0)
				//if (true)
			{
				/*printf("\r<%d> SNR= %.2f FER= %d / %d = %f BER= %d / %d = %lf E-Det = %f, E-MRB = %f"  ,
					frame, nSNR, Num_Frame_Error, frame, (double)Num_Frame_Error / frame, Num_Error, (N_Info * frame), (double)Num_Error / N_Info / frame,
					(double)n_det_err/total_bit_cnt, (double)n_det_err_mrb/(total_bit_cnt/(ModType/2)));*/
				printf("\r<%d> SNR= %.2f FER= %d / %d = %f BER= %d / %d = %lf ",
					frame, nSNR, Num_Frame_Error, frame, (double)Num_Frame_Error / frame, Num_Error, (N_Info * frame), (double)Num_Error / N_Info / frame);
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
		// save to .txt file
		string polarde_serial = "";
		for (int i = 0; i < 1; i++)
			polarde_serial = polarde_serial + polarde_array[i] + string("_");
		string file_name = string("ResFinal\\1D_LSD_D_") + polarde_serial + string("T") + to_string(Nt) + string("R") + to_string(Nr)
			+ string("Q") + to_string(int(pow(2, ModType))) + string("L") + to_string(N) + string("K") + to_string(K)
			+ string("T") + to_string(softiter) + string("H") + to_string(half_iter) + string("UB") + to_string(usebeta)
			+ string(".txt");
		save_array_txt(nSNR, BER_array[SNR_count], FER_array[SNR_count], file_name);
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
	//delete[] Ach;
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
	for (int i = 0; i < N; i++)
		delete[] GG[i];
	delete[] GG;
	for (int i = 0; i < mimo_blk_height; i++) {
		delete[] x_mod[i];
	}
	delete[] x_mod;
	for (int i = 0; i < Nr; i++) {
		delete[] y_mod[i];
	}
	delete[] y_mod;
	for (int i = 0; i < 2 * mimo_blk_height; i++) {
		delete[] y_mmse[i];
	}
	delete[] y_mmse;
	for (int i = 0; i < 2 * Nr; i++) {
		delete[] y_mod_real[i];
	}
	for (int i = 0; i < mimo_blk_length; i++)
	{
		delete[] y_mod_real_temp_omp[i];
	}
	delete[] y_mod_real_temp_omp;
	for (int i = 0; i < mimo_blk_length; i++)
	{
		delete[] y_mmse_temp_omp[i];
	}
	delete[] y_mmse_temp_omp;
	delete[] y_mod_real;
	delete[] x_mod_tmp;
	delete[] y_mod_tmp;
	delete[] y_mod_real_tmp;
	delete[] x_tmp;
	delete[] oriLLR_tmp;

	delete[] A;
	delete[] Ac;
	delete[] A_Ac;
	delete[] a_encode;
	delete[] a_encode_tmp;
	delete[] u_1D;
	delete[] x_1D;
	delete[] u_decod_1D;
	delete[] oriLLR_1D;
	delete[] x_intrlv_1D;
	for (int i = 0; i < L; i++)
		delete[]LLR[i];
	delete[] LLR;
	for (int i = 0; i < L; i++)
		delete[]sum[i];
	delete[] sum;
	for (int i = 0; i < L; i++)
		delete[]u1[i];
	delete[] u1;
	delete[] PM;
	delete[] LLR_in;
	delete[] W;
	delete[] SCL_index;
	delete[] better;
	delete[] worse;
	for (int i = 0; i < mimo_blk_length; i++)
		delete[] mu_mmse[i];
	delete[] mu_mmse;
}




void system2D_JDD_copy(double* BER_array, double* FER_array)
{
	//freopen("result.txt", "w", stdout);
	printsetting();
	cout << "\nSNR	 " << "FER      BER" << endl;
	int N_Info = Kv * Kh - CRCL, Nmax;
	//Interleaving pattern
	std::vector<int> intrlv_pattern(Nv);
	float CodeRate = (float)N_Info / (Nh * Nv);
	Nmax = max(Nh, Nv);
	// num of symbols per vertical codeword
	int Nsym = 0;
	if (MOD_DIM == "Spatial") { Nsym = Nv / ModType; }
	else if (MOD_DIM == "Temporal") { Nsym = Nh / ModType; }
	int symbol_length = ModType;

	// Size of a mimo block is given by height(spatial) * length(temporal)
	int mimo_blk_height = 0;
	int mimo_blk_length = 0;
	if (MOD_DIM == "Spatial")
	{
		mimo_blk_height = Nsym;
		mimo_blk_length = Nh;
	}
	else if (MOD_DIM == "Temporal")
	{
		mimo_blk_height = Nv;
		mimo_blk_length = Nsym;
	}
	// Initialization
	//*****************************************************************

	int** GGh = new int* [Nh];
	int** GGv = new int* [Nv];


	MatrixXf Gh = Fn(Nh);  // Generate Matrix of horizontal coding
	MatrixXf Gv = Fn(Nv);  // Generate Matrix of vertical coding


	for (int i = 0; i < Nh; i++) GGh[i] = new int[Nh];
	for (int i = 0; i < Nv; i++) GGv[i] = new int[Nv];
	for (int i = 0; i < Nh; i++)
		for (int j = 0; j < Nh; j++)
			GGh[i][j] = Gh(i, j);
	for (int i = 0; i < Nv; i++)
		for (int j = 0; j < Nv; j++)
			GGv[i][j] = Gv(i, j);
	// test GGh
	/*for (int i = 0; i < Nv; i++)
	{
		cout << endl;
		for (int j = 0; j < Nv; j++)
			cout << setw(2) << GGh[i][j];
	}
	cout << endl;*/
	// test GGh
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
	// srand((unsigned int)time(0));
	std::random_device rd;              // ��ȡ���������
	std::mt19937 gen(rd());             // ʹ�����ӳ�ʼ��mt19937������
	std::uniform_int_distribution<> distrib(1, 1000000);
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
		int n_det_err_mrb = 0; // error bits at detector output(MRBs)
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
		// test: Horizontal frozen bits

		/*cout << "Frozen Bits" << endl;
		for (int i = 0; i < Nh - Kh; i++)
		{
			cout << setw(5) << Ach[i];
		}
		cout << endl;*/
		// test: Horizontal frozen bits
		// Frame Loop
		int passnum;
		int frame = 0;
		//MatrixXcf H(Nr, Nt); // Complex channel Matrix H: Nr x Nt
		////Matrix<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic> H = MatrixXcf::Random(Nr, Nt); // Real domain channel matrix
		////// //MatrixX Definitions
		//MatrixXf H_real(2 * Nr, 2 * Nt); // Real domain channel matrix
		//MatrixXf received_signal(2 * Nr, 1);
		//MatrixXf HTH(2 * Nt, 2 * Nt);
		//MatrixXf Identity = MatrixXf::Identity(2 * Nt, 2 * Nt);
		//MatrixXf output_signal(2 * Nt, 1);
		//MatrixXf mmse_matrix(2 * Nt, 2 * Nr);

		//// for KBest
		//MatrixXf R(2 * Nt, 2 * Nt);
		//MatrixXf Q(2 * Nr, 2 * Nt);
		//MatrixXf z(2 * Nr, 1);

		//cout << frame << endl;
		int u_input[16] = { 0,     1     ,0,     1,
								0,     1,     0,     0,
								0 ,    1,     0,     0,
								1 ,    0,     0,     0 };
		//while (Num_Frame_Error < errframe)
#pragma omp parallel for
		for (int sim_frame_cnt = 0; sim_frame_cnt < 10000000 ; sim_frame_cnt++)
		{ 
			frame++;
			if (Num_Frame_Error > errframe - 1) break;
			MatrixXcf H(Nr, Nt); // Complex channel Matrix H: Nr x Nt 
			//Matrix<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic> H = MatrixXcf::Random(Nr, Nt); // Real domain channel matrix
			//// //MatrixX Definitions
			MatrixXf H_real(2 * Nr, 2 * Nt); // Real domain channel matrix
			MatrixXf received_signal(2 * Nr, 1);
			MatrixXf HTH(2 * Nt, 2 * Nt);
			MatrixXf Identity = MatrixXf::Identity(2 * Nt, 2 * Nt);
			MatrixXf output_signal(2 * Nt, 1);
			MatrixXf mmse_matrix(2 * Nt, 2 * Nr);

			// for KBest
			MatrixXf R(2 * Nt, 2 * Nt);
			MatrixXf Q(2 * Nr, 2 * Nt);
			MatrixXf z(2 * Nr, 1);

			unsigned char* a_data = new unsigned char[Kv * Kh];// Information CodeWord ok
			int** u = new int* [Nv];// Original CodeWord u ok
			for (int i = 0; i < Nv; i++) { u[i] = new int[Nh]; }
			int** umid = new int* [Nv];// output CodeWord u ok
			for (int i = 0; i < Nv; i++) { umid[i] = new int[Nh]; }
			int** uout = new int* [Nv];// output CodeWord u ok
			for (int i = 0; i < Nv; i++) { uout[i] = new int[Nh]; }
			int** midx = new int* [Nv];//  Polar Encoded CodeWord x ok
			for (int i = 0; i < Nv; i++) { midx[i] = new int[Nh]; }
			int** x = new int* [Nv];//  Polar Encoded CodeWord x ok
			for (int i = 0; i < Nv; i++) { x[i] = new int[Nh]; }
			int** uout_vertical_partial = new int* [Nv]; // ok
			for (int i = 0; i < Nv; i++) uout_vertical_partial[i] = new int[Nh];

			// Arrays for MIMO
			int** x_mmse_out = new int* [Nv]; // delete! ok
			for (int i = 0; i < Nv; i++) { x_mmse_out[i] = new int[Nh]; }
			std::complex<float>** x_mod = new std::complex<float>*[mimo_blk_height]; //ok modulated symbols (spatial modulation) DELETE!
			for (int i = 0; i < mimo_blk_height; i++) { x_mod[i] = new std::complex<float>[mimo_blk_length]; } // ok
			std::complex<float>** y_mod = new std::complex<float>*[Nr]; // Output symbols after AWGN channels // delete!
			for (int i = 0; i < Nr; i++) { y_mod[i] = new std::complex<float>[mimo_blk_length]; }
			float** y_mod_real = new float* [2 * Nr];  // ok
			for (int i = 0; i < 2 * Nr; i++) y_mod_real[i] = new float[mimo_blk_length];
			float** y = new float* [Nv];// Output Signal After AWGN Channel
			for (int i = 0; i < Nv; i++) { y[i] = new float[Nh]; } //ok
			float** y_mmse = new float* [2 * mimo_blk_height]; // ok
			for (int i = 0; i < 2 * mimo_blk_height; i++) { y_mmse[i] = new float[mimo_blk_length]; }  // MMSE detector output // delete!
			float* y_mmse_tmp = new float[2 * mimo_blk_height]; // ok
			std::complex<float>* x_mod_tmp = new std::complex<float>[mimo_blk_height]; //ok modulated symbols of a single column // delete!
			std::complex<float>* y_mod_tmp = new std::complex<float>[Nr]; //ok received symbols of a single column // delete!
			float* y_mod_real_tmp = new float[2 * Nr]; //ok
			// float* y_mod_tmp = new float[Nr];  // delete!
			float** y_mod_real_temp_omp = new float* [mimo_blk_length]; // new delete! //ok
			for (int i = 0; i < mimo_blk_length; i++) { y_mod_real_temp_omp[i] = new float[2 * Nr]; }
			float** y_mmse_temp_omp = new float* [mimo_blk_length]; // new delete!
			for (int i = 0; i < mimo_blk_length; i++) { y_mmse_temp_omp[i] = new float[2 * mimo_blk_height]; }//ok


			int* x_tmp = new int[Nv]; // delete!//ok

			// paths selected by k-best detection


			// CRC bits
			unsigned char* u_crc = new unsigned char[Kv * Kh]; //ok

			float** upLLR = new float* [Nv];//ok
			for (int i = 0; i < Nv; i++) { upLLR[i] = new _MM_ALIGN16 float[Nh]; }
			//float** oriLLRv = new float* [Nh];
			//for (int i = 0; i < Nh; i++) { oriLLRv[i] = new _MM_ALIGN16 float[2 * Nv + 2]; }
			float** oriLLRh = new float* [Nv];//ok
			for (int i = 0; i < Nv; i++) { oriLLRh[i] = new _MM_ALIGN16 float[Nh]; }
			int len_ori_llr_tmp = 0;
			if (MOD_DIM == "Spatial") { len_ori_llr_tmp = Nv; }
			else if (MOD_DIM == "Temporal") { len_ori_llr_tmp = ModType * Nt; }
			float* oriLLR_tmp = new _MM_ALIGN16 float[len_ori_llr_tmp]; //ok // delete!
			float** oriLLR_kbest = new float* [Nh];
			for (int i = 0; i < Nh; i++) { oriLLR_kbest[i] = new float[len_ori_llr_tmp]; }
			//if (Num_Frame_Error > errframe) break;
			for (int i = 0; i < Nv; i++)
			{
				for (int j = 0; j < Nv; j++)
				{
					uout_vertical_partial[i][j] = 0;
				}
			}
			// generate interleaving pattern
			//cout << 0 << endl;
			if (flag_interleave) { generateIntrlvPattern(intrlv_pattern, Nv); }

			for (int i = 0; i < N_Info; i++)
			{
				a_data[i] = (unsigned char)(distrib(gen) % 2);
				//		a_data[i] = u_input[i];
			}

			//cout << endl;
			if (CRCL) tx_append_crc(a_data, N_Info, CRCL);
			for (int i = 0; i < Nv; i++)
				for (int j = 0; j < Nh; j++)
					u[i][j] = 0;
			//cout << 3 << endl;
//#pragma omp parallel for
			for (int i = 0; i < Kv; i++)
			{
				for (int j = 0; j < Kh; j++)
				{
					u[Av[i]][Ah[j]] = (int)a_data[i * Kh + j];
				}
			}
			//cout << 2 << endl;
			// test:original info
			//cout << "original info u:" << endl;
			////for (int j = 0; j < Nh; j++) { cout << u[Av[0]][j] << " "; }
			//for (int i = 0; i < Nv; i++)
			//{
			//	cout << endl;
			//	for (int j = 0; j < Nh; j++)
			//		cout << u[i][j];
			//}
			//cout << endl;
			// test: original info
			// Polar Encoding
			// Polar Encoding
//#pragma omp parallel for
			for (int i = 0; i < Nv; i++)
				PolarEncode_xor(midx[i], u[i], Nh);
			// test:horizontal encoding
			//cout << "u after horizontal encoding:" << endl;
			////for (int j = 0; j < Nh; j++) { cout << midx[Av[0]][j] << " "; }
			//for (int i = 0; i < Nv; i++)
			//{
			//	cout << endl;
			//	for (int j = 0; j < Nh; j++)
			//		cout << midx[i][j];
			//}
			//cout << endl;
			//cout << endl;
			// test:horizontal encoding
			//cout << 1 << endl;
//#pragma omp parallel for
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
			//cout << 4 << endl;
			// test:vertical encoding
			//cout << "u after vertical encoding:" << endl;
			////for (int j = 0; j < Nh; j++) { cout << x[Av[0]][j] << " "; }
			//for (int i = 0; i < Nv; i++)
			//{
			//	cout << endl;
			//	for (int j = 0; j < Nh; j++)
			//		cout << x[i][j];
			//}
			//cout << endl;
			//cout << endl;
			// test:vertical encoding

			// test:horizontal reverse-encoding
			/*cout << "after horizontal reverse-encoding" << endl;
			for (int i = 0; i < Nv; i++)
			{
				cout << endl;
				int* u_hori = new int[Nh];
				PolarEncode_xor(u_hori, x[i], Nh);
				for (int j = 0; j < Nh; j++) { cout << u_hori[j]; }
				delete[] u_hori;
			}
			cout << endl;*/
			// test:horizontal reverse-encoding


			if (atten == 2)
			{
				// modulation
				// #pragma omp parallel for
				// Interleaving
				if (flag_interleave)
				{
					for (int j = 0; j < Nh; j++)
					{
						for (int i = 0; i < Nv; i++) { x_tmp[i] = x[i][j]; }
						randIntrlv(intrlv_pattern, x_tmp, Nv);
						for (int i = 0; i < Nv; i++) { x[i][j] = x_tmp[i]; }
					}
				}

				// for JDD system, currently we only implement spatial modulation
				if (MOD_DIM == "Spatial")
				{
					for (int i = 0; i < Nh; i++)
					{
						// Extracting a single column
						for (int j = 0; j < Nv; j++) { x_tmp[j] = x[j][i]; }
						modulation(x_tmp, x_mod_tmp, ModType, Nt);
						for (int j = 0; j < Nsym; j++) { x_mod[j][i] = x_mod_tmp[j]; }
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
					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < Nr; j++) { y_mod_tmp[j] = y_mod[j][i]; }
						convertYToReal(Nr, y_mod_tmp, y_mod_real_tmp);
						for (int j = 0; j < 2 * Nr; j++) { y_mod_real[j][i] = y_mod_real_tmp[j]; /*cout << y_mod_real_tmp[j] << endl;*/ }
					}
					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < 2 * Nr; j++) { y_mod_real_temp_omp[i][j] = y_mod_real[j][i]; }
					}
					// for KBest Detection
					Eigen::HouseholderQR<Eigen::MatrixXf> qr_Hreal(H_real);
					Q = qr_Hreal.householderQ();
					R = qr_Hreal.matrixQR().triangularView<Eigen::Upper>();
					//#pragma omp parallel for private(z)
					// test for llr to detector
					/*float** llr_to_det_test = new float* [Nv];
					for (int i = 0; i < Nv; i++) llr_to_det_test[i] = new float[Nh];
					for (int i = 0; i < Nv; i++) { for (int j = 0; j < Nh; j++) llr_to_det_test[i][j] = 0; }*/
					// test for llr to detector
					//cout << 4 << endl;
					for (int k = 0; k < Nh; k++)
					{
						// cout << k << endl;
						int j = Nh - 1 - k;
						bool isJDDCheckpoint = my_find(Ach, j, Nh - Kh);
						MatrixXf received_signal_omp(2 * Nr, 1);
						int k_kbest_extend = pow(2, ModType / 2) * K_KBEST;
						float** paths_kbest = new float* [k_kbest_extend];
						for (int i = 0; i < k_kbest_extend; i++) paths_kbest[i] = new float[2 * Nt];
						float* ED = new float[K_KBEST];
						float* llr_to_detector = new float[Nv];
						float* llr_rearranged_tmp = new float[Nv];
						float* llr_ordered_tmp = new float[Nv];
						float** llr_rearranged_array = new float* [Nv];

						for (int i = 0; i < Nv; i++) { llr_rearranged_array[i] = new float[Nh]; }
						//cout << 200 << endl;
						if (MIMODET == "MMSE") { MMSEdetection(H_real, Identity, received_signal, mmse_matrix, output_signal, HTH, y_mod_real_temp_omp[j], y_mmse_temp_omp[j], sigma_sq); }
						for (int i = 0; i < 2 * Nr; i++) { received_signal_omp(i, 0) = y_mod_real_temp_omp[j][i]; }
						//MatrixXf Qtrans()
						z = Q.transpose() * received_signal_omp;
						//for (int s = 0; s < 32; s++) { cout << z(s, 0) << endl; }
						//cout << 200 << endl;
						if (MIMODET == "KBEST")
						{
							if (!isJDDCheckpoint) {
								KBestDetection(R, z, H_real, y_mod_real_temp_omp[j], ED, y_mmse_temp_omp[j], paths_kbest, K_KBEST, ModType, Nt, Nr);
							}
							// else if is JDD Checkpoint, activate KBEST_JDD
							else
							{
								// rearrange the LLRs from decoder
								for (int i_col = j + 1; i_col < Nh; i_col++)
								{
									for (int i_bit = 0; i_bit < Nv; i_bit++) { llr_ordered_tmp[i_bit] = upLLR[i_bit][i_col]; }
									rearrange_dec_det(llr_ordered_tmp, llr_rearranged_tmp, ModType, Nt);
									for (int i_bit = 0; i_bit < Nv; i_bit++) { llr_rearranged_array[i_bit][i_col] = llr_rearranged_tmp[i_bit]; }
								}
								// calculate the LLR to be passed to the detector
								calc_llr_idd(llr_to_detector, llr_rearranged_array, GGh, j, Nh, Nv);
								/*if (flag_interleave)
								{
									randIntrlv(intrlv_pattern, llr_to_detector, Nv);
								}*/
								// test: llr_idd
								//for (int j_row = 0; j_row < Nv; j_row++) { llr_to_det_test[j_row][j] = llr_to_detector[j_row]; }
								// test: llr_idd
								KBestDetection_soft_JDD(R, z, H_real, y_mod_real_temp_omp[j], ED, y_mmse_temp_omp[j], llr_to_detector, paths_kbest, K_KBEST, ModType, Nt, Nr, j);
							}
							llr_kbest_bit(Nt, y_mmse_temp_omp[j], oriLLR_kbest[j], sigma_sq, ModType, K_KBEST, paths_kbest, ED);
						}
						// for (int i_path = 0; i_path < k_kbest_extend; i_path++) { delete[] paths_kbest[i_path]; }

						// mmse: transform symbol LLR to bit LLR
						// kbest: copy from oriLLR_kbest to oriLLRh and execute vertical decoding
						//cout << 200 << endl;
						for (int i = 0; i < 2 * Nt; i++) { y_mmse[i][j] = y_mmse_temp_omp[j][i]; }
						for (int i = 0; i < 2 * Nt; i++) { y_mmse_tmp[i] = y_mmse[i][j]; }
						if (MIMODET == "MMSE") { llr_mmse_bit(Nt, y_mmse_tmp, oriLLR_tmp, sigma_sq, ModType); }
						if (MIMODET == "KBEST")
						{
							for (int i = 0; i < len_ori_llr_tmp; i++) { oriLLR_tmp[i] = oriLLR_kbest[j][i]; }
							// copy LLR
							for (int i = 0; i < Nv; i++)
							{
								oriLLRh[i][j] = oriLLR_tmp[i];
								upLLR[i][j] = oriLLR_tmp[i];
							}
							// First round of vertical polar decoding
							PolarDecodeIterVertical(oriLLRh, upLLR, umid, GGv, uout_vertical_partial, A_Acv, Av, Nh, Nv, j, 1);
							// test: oriLLRh
							/*int num_err_dec = 0;
							cout << "Dec from oriLLRh" << endl;
							for (int i = 0; i < Nv; i++)
							{
								cout << endl;
								for (int j = 0; j < Nh; j++)
								{
									cout << setw(2) << (oriLLRh[i][j] < 0);
									num_err_dec += abs(x[i][j] - (oriLLRh[i][j] < 0));
								}
							}
							cout << endl;
							cout << "Dec from oriLLRh" << endl << k << endl;
							cout << "Detector output error:  " << num_err_dec << endl;*/
							// test: oriLLRh
							// test: uout_vertical_partial
							/*cout << "uout_vertical_partial" << endl;
							cout << "k= " << k << endl;
							for (int i = 0; i < Nv; i++)
							{
								cout << endl;
								for (int j = 0; j < Nh; j++)
								{
									cout << uout_vertical_partial[i][j];
								}
							}*/
							// test: uout_vertical_partial
						}
						//cout << 2000 << endl;
						for (int i = 0; i < k_kbest_extend; i++) delete[] paths_kbest[i];
						delete[] paths_kbest;
						for (int i = 0; i < Nv; i++) delete[] llr_rearranged_array[i];
						delete[] llr_rearranged_array;
						delete[] ED;
						delete[] llr_rearranged_tmp;
						delete[] llr_ordered_tmp;
						delete[] llr_to_detector;
					}
					//cout << 5 << endl;
					// Temporal modulation
					// test for llr to detector
					/*cout << "IDD Dec" << endl;
					for (int i = 0; i < Nv; i++)
					{
						cout << endl;
						for (int j = 0; j < Nh; j++)
							cout << setw(2) << (llr_to_det_test[i][j]< 0);
					}
					cout << "IDD Dec" << endl;
					cout << "IDD LLR" << endl;
					for (int i = 0; i < Nv; i++)
					{
						cout << endl;
						for (int j = 0; j < Nh; j++)
							cout << setw(10) << llr_to_det_test[i][j] ;
					}
					cout << "IDD LLR" << endl;
					for (int i = 0; i < Nv; i++) delete[] llr_to_det_test[i];
					delete[] llr_to_det_test;*/
					// test for llr to detector

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
			//cout << 10000 << endl;
			for (int iter = 1; iter <= softiter; iter++)
			{
				//	cout << "Vertical Start:\n";
					//vertical decoding
				//cout << "Horizon Start:\n";
				//horizontal decoding
				// #
//#pragma omp parallel for
				for (int decodingi = 0; decodingi < Nh; decodingi++) // 291.70kb
				{
					if (iter == 1) break;
					//for SCL Decoding
					int Counth = 0, Countinfoh = 0;
					float Qsumunh = 0;
					int** sum = new int* [L];
					//for (int i = 0; i < L; i++)  sum[i] = new  int[Nmax];
					int** u1 = new int* [L];
					for (int i = 0; i < L; i++)  u1[i] = new _MM_ALIGN16 int[Nmax];
					for (int i = 0; i < L; i++)  sum[i] = new _MM_ALIGN16 int[Nmax];
					//for (int i = 0; i < L; i++)  u1[i] = new  int[Nmax];

					float* LLR_in = new float[L];
					float* W = new float[2 * L];
					int* SCL_index = new int[L];
					int* better = new int[L];
					int* worse = new int[L];
					int* indexpm = new int[L];
					float* PM = new float[L];
					float** LLR = new float* [L];
					for (int i = 0; i < L; i++)  LLR[i] = new _MM_ALIGN16 float[2 * Nmax + 2];
					// for (int i = 0; i < L; i++)  LLR[i] = new float[2 * Nmax + 2];
					////////////////////////////////
					for (int j = 0; j < Nv; j++)
						for (int k = 0; k < L; k++)
							LLR[k][j] = upLLR[j][decodingi];
					//cout << LLR[0][j] << " ";

				//cout << endl;
					int index = 0;
					/*			if (L == 1)
								{
									//T_start = clock();
									Count = 0;
									decode(*LLR, *sum, (int)(log(Nv) / log(2)), A_Acv, flag, Nv, Kv);
									//PolarEncode_xor(*u1, *sum, N);
									//T_finish = clock();
								}
								else
								{*/
					for (int k = 0; k < L; k++) { PM[k] = 0; indexpm[k] = k; }
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
					if (adaptive_beta)
					{
						if (softmethod == "Pyndiah") Pyndiahv_adaptivebeta(LLR, sum, oriLLRh, indexpm, iter * 2 - 2, decodingi, upLLR);
						if (softmethod == "MITSO") MITSOv(LLR, PM, sum, oriLLRh, indexpm, iter * 2 - 2, decodingi, upLLR, Qsumunh);
					}
					else
					{
						if (softmethod == "Pyndiah") Pyndiahv(LLR, sum, oriLLRh, indexpm, iter * 2 - 2, decodingi, upLLR);
						if (softmethod == "MITSO") MITSOv(LLR, PM, sum, oriLLRh, indexpm, iter * 2 - 2, decodingi, upLLR, Qsumunh);
					}
					// half iteration Enabled
					if (half_iter && iter == softiter)
					{
						if (iter == softiter) PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), Nh);
						for (int j = 0; j < Kv; j++)
							umid[Av[j]][decodingi] = u1[index][Av[j]];
					}
					// half iteration
					////////delete SCL variables
					for (int i = 0; i < L; i++)
						delete[]LLR[i];
					delete[] LLR;
					for (int i = 0; i < L; i++)
						delete[]sum[i];
					for (int i = 0; i < L; i++)
						delete[]u1[i];
					delete[] u1;
					delete[] sum;
					delete[] PM;
					delete[] LLR_in;
					delete[] W;
					delete[] SCL_index;
					delete[] better;
					delete[] worse;
					delete[] indexpm;

				}
				if (weirditer) { continue; }

				//#pragma omp parallel for
								//cout << 10000 << endl;
				for (int decodingi = 0; decodingi < Nv; decodingi++)
				{
					// 293.70kb
					// cout << "KKKK" << endl;
					if (half_iter && (iter == softiter)) { break; }
					string polarde_now = polarde_array[iter - 1];
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
					//cout << 10000 << endl;
					for (int i = 0; i < L_LSD; i++) { indexpm[i] = i; }
					////////////////////////////////
					for (int j = 0; j < Nh; j++)
						for (int k = 0; k < L; k++)
							LLR[k][j] = upLLR[decodingi][j];
					//	cout << LLR[0][j] << " ";
					calc_npaths(A_Ach, Nh, n_paths, L_LSD);
					//	cout << endl;
					int index = 0;
					if (polarde_now == "SCL")
					{
						for (int k = 0; k < L; k++)  PM[k] = 0, indexpm[k] = k;
						//T_start = clock();
						decode_list(LLR, sum, (int)(log(Nh) / log(2)), A_Ach, 0, Nh, Kh, PM, L, (int)(log(L) / log(2)), LLR_in, W, SCL_index, better, worse, &Countv, &Countinfov, &Qsumunv);
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
						if (iter == softiter) PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), Nh);
						index = indexpm[0];
						//}
						//calc distance s()
						if (iter < softiter)
						{
							if (adaptive_beta)
							{
								if (softmethod == "Pyndiah") Pyndiahh_adaptivebeta(LLR, sum, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR);
								if (softmethod == "MITSO") MITSOh(LLR, PM, sum, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR, Qsumunv);
								//	getchar();
							}
							else
							{
								if (softmethod == "Pyndiah") Pyndiahh(LLR, sum, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR);
								if (softmethod == "MITSO") MITSOh(LLR, PM, sum, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR, Qsumunv);
								//	getchar();
							}
						}
						else
							for (int j = 0; j < Kh; j++)
								umid[decodingi][Ah[j]] = u1[index][Ah[j]];
					}

					////////delete SCL variables
					//cout << 10000 << endl;
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
					//cout << 10000 << endl;
				}
			}
			//cout << 10000 << endl;
			int calcsum[Nv], calcu1[Nv];
			int calcu1h[Nh];
			// calcsum is from the direct output of the decoder
			// if output is the result of column decoding (meaning there is another row decoding to be done)
			if (half_iter)
			{
				for (int i = 0; i < Kv; i++)
				{
					PolarEncode_xor(calcu1h, umid[Av[i]], Nh);
					for (int j = 0; j < Kh; j++)
						uout[Av[i]][Ah[j]] = calcu1h[Ah[j]];
				}
			}
			// if output is the result of column decoding (meaning that there is another row decoding to be done)
			else
			{
				for (int i = 0; i < Kh; i++)
				{
					for (int j = 0; j < Nv; j++)
						calcsum[j] = umid[j][Ah[i]];
					PolarEncode_xor(calcu1, calcsum, Nv);
					for (int j = 0; j < Kv; j++)
						uout[Av[j]][Ah[i]] = calcu1[Av[j]];
				}
			}
			// test: final output
			/*cout << "final output" << endl;
			for (int k = 0; k < Nh; k++) { cout << uout[Av[0]][k] << " "; }
			cout << endl;*/
			// test: final output
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
			if (Temp_Error < Num_Error)
			{
				Num_Frame_Error++;
			}
			if (frame % 10 == 0)
			{
				/*printf("\r<%d> SNR= %.2f FER= %d / %d = %f BER= %d / %d = %lf E-Det = %f, E-MRB = %f"  ,
					frame, nSNR, Num_Frame_Error, frame, (double)Num_Frame_Error / frame, Num_Error, (N_Info * frame), (double)Num_Error / N_Info / frame,
					(double)n_det_err/total_bit_cnt, (double)n_det_err_mrb/(total_bit_cnt/(ModType/2)));*/
				printf("\r<%d> SNR= %.2f FER= %d / %d = %f BER= %d / %d = %lf ",
					frame, nSNR, Num_Frame_Error, frame, (double)Num_Frame_Error / frame, Num_Error, (N_Info * frame), (double)Num_Error / N_Info / frame);
				fflush(stdout);
			}
			//printf("%d %d", Num_Frame_Error, Num_Error);
			//getchar();
			//Decoding_Duration += (T_finish - T_start);

			// 391233.26kB
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
				delete[] uout_vertical_partial[i];
			delete[] uout_vertical_partial;
			for (int i = 0; i < Nv; i++)
				delete[] y[i];
			for (int i = 0; i < Nh; i++)
				delete[] oriLLR_kbest[i];
			delete[] oriLLR_kbest;
			delete[] y;
			delete[] u_crc;
			for (int i = 0; i < Nv; i++)
				delete[] oriLLRh[i];
			delete[] oriLLRh;
			for (int i = 0; i < Nv; i++)
				delete[] upLLR[i];
			delete[] upLLR;
			/*for (int i = 0; i < Nh; i++)
				delete[] GGh[i];
			delete[] GGh;
			for (int i = 0; i < Nv; i++)
				delete[] GGv[i];
			delete[] GGv;*/
			for (int i = 0; i < mimo_blk_height; i++) {
				delete[] x_mod[i];
			}
			delete[] x_mod;
			for (int i = 0; i < Nr; i++) {
				delete[] y_mod[i];
			}
			delete[] y_mod;
			for (int i = 0; i < 2 * mimo_blk_height; i++) {
				delete[] y_mmse[i];
			}
			delete[] y_mmse;
			for (int i = 0; i < 2 * Nr; i++) {
				delete[] y_mod_real[i];
			}
			for (int i = 0; i < mimo_blk_length; i++)
			{
				delete[] y_mod_real_temp_omp[i];
			}
			delete[] y_mod_real_temp_omp;
			for (int i = 0; i < mimo_blk_length; i++)
			{
				delete[] y_mmse_temp_omp[i];
			}
			/*for (int i = 0; i < Nh; i++)
				delete[] oriLLR_kbest[i];
			delete[] oriLLR_kbest;*/
			delete[] y_mmse_temp_omp;
			delete[] y_mod_real;
			delete[] x_mod_tmp;
			delete[] y_mod_tmp;
			delete[] y_mod_real_tmp;
			delete[] x_tmp;
			delete[] oriLLR_tmp;
			delete[] y_mmse_tmp;
		}
		//Decoding_Duration = Decoding_Duration * 100 / CLOCKS_PER_SEC / frame;
		//Decoding_Duration = Decoding_Duration / CLOCKS_PER_SEC / Max_frame;

		//printf("%f  %f  %f  %f\n", nSNR, Decoding_Duration, (double)Num_Error / N_Info / Max_frame, (double)Num_Frame_Error / Max_frame);
		printf("%f %f %f\n", nSNR, (double)Num_Frame_Error / frame, (double)Num_Error / N_Info / frame);
		BER_array[SNR_count] = (double)Num_Error / N_Info / frame;
		FER_array[SNR_count] = (double)Num_Frame_Error / frame;
		// save to .txt file
		string polarde_serial = "";
		for (int i = 0; i < softiter; i++)
			polarde_serial = polarde_serial + polarde_array[i] + string("_");
		string file_name = string("ResFinal\\2D_LSD_D_") + polarde_serial + MIMODET + string("T") + to_string(Nt) + string("R") + to_string(Nr)
			+ string("Q") + to_string(int(pow(2, ModType))) + string("L") + to_string(Nh) + string("K") + to_string(Kh)
			+ string("T") + to_string(softiter) + string("H") + to_string(half_iter) + string("UB") + to_string(usebeta)
			+ string(".txt");
		save_array_txt(nSNR, BER_array[SNR_count], FER_array[SNR_count], file_name);
		SNR_count++;
	}
	//fclose(stdout);

	delete[] Ah;
	delete[] Ach;
	delete[] A_Ach;
	delete[] Av;
	delete[] Acv;
	delete[] A_Acv;

}






void system2D_JDD(double* BER_array, double* FER_array)
{
	//freopen("result.txt", "w", stdout);
	printsetting();
	cout << "\nSNR	 " << "FER      BER" << endl;
	int N_Info = Kv * Kh - CRCL, Nmax;
	//Interleaving pattern
	std::vector<int> intrlv_pattern(Nv);
	float CodeRate = (float)N_Info / (Nh * Nv);
	Nmax = max(Nh, Nv);
	// num of symbols per vertical codeword
	int Nsym = 0;
	if (MOD_DIM == "Spatial") { Nsym = Nv / ModType; }
	else if (MOD_DIM == "Temporal") { Nsym = Nh / ModType; }
	int symbol_length = ModType;

	// Size of a mimo block is given by height(spatial) * length(temporal)
	int mimo_blk_height = 0;
	int mimo_blk_length = 0;
	if (MOD_DIM == "Spatial")
	{
		mimo_blk_height = Nsym;
		mimo_blk_length = Nh;
	}
	else if (MOD_DIM == "Temporal")
	{
		mimo_blk_height = Nv;
		mimo_blk_length = Nsym;
	}
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
	int** uout_vertical_partial = new int* [Nv];
	for (int i = 0; i < Nv; i++) uout_vertical_partial[i] = new int[Nh];

	// Arrays for MIMO
	int** x_mmse_out = new int* [Nv]; // delete!
	for (int i = 0; i < Nv; i++) { x_mmse_out[i] = new int[Nh]; }
	std::complex<float>** x_mod = new std::complex<float>*[mimo_blk_height]; // modulated symbols (spatial modulation) DELETE!
	for (int i = 0; i < mimo_blk_height; i++) { x_mod[i] = new std::complex<float>[mimo_blk_length]; }
	std::complex<float>** y_mod = new std::complex<float>*[Nr]; // Output symbols after AWGN channels // delete!
	for (int i = 0; i < Nr; i++) { y_mod[i] = new std::complex<float>[mimo_blk_length]; }
	float** y_mod_real = new float* [2 * Nr];
	for (int i = 0; i < 2 * Nr; i++) y_mod_real[i] = new float[mimo_blk_length];
	float** y = new float* [Nv];// Output Signal After AWGN Channel
	for (int i = 0; i < Nv; i++) { y[i] = new float[Nh]; }
	float** y_mmse = new float* [2 * mimo_blk_height];
	for (int i = 0; i < 2 * mimo_blk_height; i++) { y_mmse[i] = new float[mimo_blk_length]; }  // MMSE detector output // delete!
	float* y_mmse_tmp = new float[2 * mimo_blk_height];
	std::complex<float>* x_mod_tmp = new std::complex<float>[mimo_blk_height]; // modulated symbols of a single column // delete!
	std::complex<float>* y_mod_tmp = new std::complex<float>[Nr]; // received symbols of a single column // delete!
	float* y_mod_real_tmp = new float[2 * Nr];
	// float* y_mod_tmp = new float[Nr];  // delete!
	float** y_mod_real_temp_omp = new float* [mimo_blk_length]; // new delete!
	for (int i = 0; i < mimo_blk_length; i++) { y_mod_real_temp_omp[i] = new float[2 * Nr]; }
	float** y_mmse_temp_omp = new float* [mimo_blk_length]; // new delete!
	for (int i = 0; i < mimo_blk_length; i++) { y_mmse_temp_omp[i] = new float[2 * mimo_blk_height]; }


	int* x_tmp = new int[Nv]; // delete!

	// paths selected by k-best detection


	// CRC bits
	unsigned char* u_crc = new unsigned char[Kv * Kh];

	float** upLLR = new float* [Nv];
	for (int i = 0; i < Nv; i++) { upLLR[i] = new _MM_ALIGN16 float[Nh]; }
	//float** oriLLRv = new float* [Nh];
	//for (int i = 0; i < Nh; i++) { oriLLRv[i] = new _MM_ALIGN16 float[2 * Nv + 2]; }
	float** oriLLRh = new float* [Nv];
	for (int i = 0; i < Nv; i++) { oriLLRh[i] = new _MM_ALIGN16 float[Nh]; }
	int len_ori_llr_tmp = 0;
	if (MOD_DIM == "Spatial") { len_ori_llr_tmp = Nv; }
	else if (MOD_DIM == "Temporal") { len_ori_llr_tmp = ModType * Nt; }
	float* oriLLR_tmp = new _MM_ALIGN16 float[len_ori_llr_tmp];  // delete!
	float** oriLLR_kbest = new float* [Nh];
	for (int i = 0; i < Nh; i++) { oriLLR_kbest[i] = new float[len_ori_llr_tmp]; }

	//*****************************************************************

	int** GGh = new int* [Nh];
	int** GGv = new int* [Nv];


	MatrixXf Gh = Fn(Nh);  // Generate Matrix of horizontal coding
	MatrixXf Gv = Fn(Nv);  // Generate Matrix of vertical coding
	MatrixXcf H(Nr, Nt); // Complex channel Matrix H: Nr x Nt
	//Matrix<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic> H = MatrixXcf::Random(Nr, Nt); // Real domain channel matrix
	//// //MatrixX Definitions
	MatrixXf H_real(2 * Nr, 2 * Nt); // Real domain channel matrix
	MatrixXf received_signal(2 * Nr, 1);
	MatrixXf HTH(2 * Nt, 2 * Nt);
	MatrixXf Identity = MatrixXf::Identity(2 * Nt, 2 * Nt);
	MatrixXf output_signal(2 * Nt, 1);
	MatrixXf mmse_matrix(2 * Nt, 2 * Nr);

	// for KBest
	MatrixXf R(2 * Nt, 2 * Nt);
	MatrixXf Q(2 * Nr, 2 * Nt);
	MatrixXf z(2 * Nt, 1);

	for (int i = 0; i < Nh; i++) GGh[i] = new int[Nh];
	for (int i = 0; i < Nv; i++) GGv[i] = new int[Nv];
	for (int i = 0; i < Nh; i++)
		for (int j = 0; j < Nh; j++)
			GGh[i][j] = Gh(i, j);
	for (int i = 0; i < Nv; i++)
		for (int j = 0; j < Nv; j++)
			GGv[i][j] = Gv(i, j);
	// test GGh
	/*for (int i = 0; i < Nv; i++)
	{
		cout << endl;
		for (int j = 0; j < Nv; j++)
			cout << setw(2) << GGh[i][j];
	}
	cout << endl;*/
	// test GGh
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
	// srand((unsigned int)time(0));
	std::random_device rd;              // ��ȡ���������
	std::mt19937 gen(rd());             // ʹ�����ӳ�ʼ��mt19937������
	std::uniform_int_distribution<> distrib(1, 1000000);
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
		int n_det_err_mrb = 0; // error bits at detector output(MRBs)
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
		// test: Horizontal frozen bits

		/*cout << "Frozen Bits" << endl;
		for (int i = 0; i < Nh - Kh; i++)
		{
			cout << setw(5) << Ach[i];
		}
		cout << endl;*/
		// test: Horizontal frozen bits
		// Frame Loop
		int passnum;
		int frame = 0;
		while (Num_Frame_Error < errframe){
//#pragma omp parallel for
		/*for (int sim_frame_cnt = 0; sim_frame_cnt < 10000000; sim_frame_cnt++)
		{*/
			frame++;
			//if (Num_Frame_Error > errframe) break;
			for (int i = 0; i < Nv; i++)
			{
				for (int j = 0; j < Nh; j++)
				{
					uout_vertical_partial[i][j] = 0;
				}
			}
			// generate interleaving pattern
			if (flag_interleave) { generateIntrlvPattern(intrlv_pattern, Nv); }

			//cout << frame << endl;
			int u_input[16] = { 0,     1     ,0,     1,
									0,     1,     0,     0,
									0 ,    1,     0,     0,
									1 ,    0,     0,     0 };
			for (int i = 0; i < N_Info; i++)
			{
				a_data[i] = (unsigned char)(distrib(gen) % 2);
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
				}
			}
			// test:original info
			//cout << "original info u:" << endl;
			////for (int j = 0; j < Nh; j++) { cout << u[Av[0]][j] << " "; }
			//for (int i = 0; i < Nv; i++)
			//{
			//	cout << endl;
			//	for (int j = 0; j < Nh; j++)
			//		cout << u[i][j];
			//}
			//cout << endl;
			// test: original info
			// Polar Encoding
			// Polar Encoding
#pragma omp parallel for
			for (int i = 0; i < Nv; i++)
				PolarEncode_xor(midx[i], u[i], Nh);
			// test:horizontal encoding
			//cout << "u after horizontal encoding:" << endl;
			////for (int j = 0; j < Nh; j++) { cout << midx[Av[0]][j] << " "; }
			//for (int i = 0; i < Nv; i++)
			//{
			//	cout << endl;
			//	for (int j = 0; j < Nh; j++)
			//		cout << midx[i][j];
			//}
			//cout << endl;
			//cout << endl;
			// test:horizontal encoding

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
			// test:vertical encoding
			//cout << "u after vertical encoding:" << endl;
			////for (int j = 0; j < Nh; j++) { cout << x[Av[0]][j] << " "; }
			//for (int i = 0; i < Nv; i++)
			//{
			//	cout << endl;
			//	for (int j = 0; j < Nh; j++)
			//		cout << x[i][j];
			//}
			//cout << endl;
			//cout << endl;
			// test:vertical encoding

			// test:horizontal reverse-encoding
			/*cout << "after horizontal reverse-encoding" << endl;
			for (int i = 0; i < Nv; i++)
			{
				cout << endl;
				int* u_hori = new int[Nh];
				PolarEncode_xor(u_hori, x[i], Nh);
				for (int j = 0; j < Nh; j++) { cout << u_hori[j]; }
				delete[] u_hori;
			}
			cout << endl;*/
			// test:horizontal reverse-encoding


			if (atten == 2)
			{
				// modulation
				// #pragma omp parallel for
				// Interleaving
				if (flag_interleave)
				{
					for (int j = 0; j < Nh; j++)
					{
						for (int i = 0; i < Nv; i++) { x_tmp[i] = x[i][j]; }
						randIntrlv(intrlv_pattern, x_tmp, Nv);
						for (int i = 0; i < Nv; i++) { x[i][j] = x_tmp[i]; }
					}
				}

				// for JDD system, currently we only implement spatial modulation
				if (MOD_DIM == "Spatial")
				{
					for (int i = 0; i < Nh; i++)
					{
						// Extracting a single column
						for (int j = 0; j < Nv; j++) { x_tmp[j] = x[j][i]; }
						modulation(x_tmp, x_mod_tmp, ModType, Nt);
						for (int j = 0; j < Nsym; j++) { x_mod[j][i] = x_mod_tmp[j]; }
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
					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < Nr; j++) { y_mod_tmp[j] = y_mod[j][i]; }
						convertYToReal(Nr, y_mod_tmp, y_mod_real_tmp);
						for (int j = 0; j < 2 * Nr; j++) { y_mod_real[j][i] = y_mod_real_tmp[j]; /*cout << y_mod_real_tmp[j] << endl;*/ }
					}
					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < 2 * Nr; j++) { y_mod_real_temp_omp[i][j] = y_mod_real[j][i]; }
					}
					// for KBest Detection
					Eigen::HouseholderQR<Eigen::MatrixXf> qr_Hreal(H_real);
					Q = qr_Hreal.householderQ();
					R = qr_Hreal.matrixQR().triangularView<Eigen::Upper>();
					//#pragma omp parallel for private(z)
					// test for llr to detector
					/*float** llr_to_det_test = new float* [Nv];
					for (int i = 0; i < Nv; i++) llr_to_det_test[i] = new float[Nh];
					for (int i = 0; i < Nv; i++) { for (int j = 0; j < Nh; j++) llr_to_det_test[i][j] = 0; }*/
					// test for llr to detector
					for (int k = 0; k < Nh; k++)
					{
						// cout << k << endl;
						int j = Nh - 1 - k;
						bool isJDDCheckpoint = my_find(Ach, j, Nh - Kh);
						MatrixXf received_signal_omp(2 * Nr, 1);
						int k_kbest_extend = pow(2, ModType / 2) * K_KBEST;
						float** paths_kbest = new float* [k_kbest_extend];
						for (int i = 0; i < k_kbest_extend; i++) paths_kbest[i] = new float[2 * Nt];
						float* ED = new float[K_KBEST];
						float* llr_to_detector = new float[Nv];
						float* llr_rearranged_tmp = new float[Nv];
						float* llr_ordered_tmp = new float[Nv];
						float** llr_rearranged_array = new float* [Nv];
						for (int i = 0; i < Nv; i++) { llr_rearranged_array[i] = new float[Nh]; }

						if (MIMODET == "MMSE") { MMSEdetection(H_real, Identity, received_signal, mmse_matrix, output_signal, HTH, y_mod_real_temp_omp[j], y_mmse_temp_omp[j], sigma_sq); }
						for (int i = 0; i < 2 * Nr; i++) { received_signal_omp(i, 0) = y_mod_real_temp_omp[j][i]; }
						z = Q.transpose() * received_signal_omp;
						if (MIMODET == "KBEST")
						{
							if (!isJDDCheckpoint) { KBestDetection(R, z, H_real, y_mod_real_temp_omp[j], ED, y_mmse_temp_omp[j], paths_kbest, K_KBEST, ModType, Nt, Nr); }
							// else if is JDD Checkpoint, activate KBEST_JDD
							else
							{
								// rearrange the LLRs from decoder
								for (int i_col = j + 1; i_col < Nh; i_col++)
								{
									for (int i_bit = 0; i_bit < Nv; i_bit++) { llr_ordered_tmp[i_bit] = upLLR[i_bit][i_col]; }
									rearrange_dec_det(llr_ordered_tmp, llr_rearranged_tmp, ModType, Nt);
									for (int i_bit = 0; i_bit < Nv; i_bit++) { llr_rearranged_array[i_bit][i_col] = llr_rearranged_tmp[i_bit]; }
								}
								// calculate the LLR to be passed to the detector
								calc_llr_idd(llr_to_detector, llr_rearranged_array, GGh, j, Nh, Nv);
								/*if (flag_interleave)
								{
									randIntrlv(intrlv_pattern, llr_to_detector, Nv);
								}*/
								// test: llr_idd
								//for (int j_row = 0; j_row < Nv; j_row++) { llr_to_det_test[j_row][j] = llr_to_detector[j_row]; }
								// test: llr_idd
								KBestDetection_soft_JDD(R, z, H_real, y_mod_real_temp_omp[j], ED, y_mmse_temp_omp[j], llr_to_detector, paths_kbest, K_KBEST, ModType, Nt, Nr, j);
							}
							llr_kbest_bit(Nt, y_mmse_temp_omp[j], oriLLR_kbest[j], sigma_sq, ModType, K_KBEST, paths_kbest, ED);
						}
						// for (int i_path = 0; i_path < k_kbest_extend; i_path++) { delete[] paths_kbest[i_path]; }

						// mmse: transform symbol LLR to bit LLR
						// kbest: copy from oriLLR_kbest to oriLLRh and execute vertical decoding
						for (int i = 0; i < 2 * Nt; i++) { y_mmse[i][j] = y_mmse_temp_omp[j][i]; }
						for (int i = 0; i < 2 * Nt; i++) { y_mmse_tmp[i] = y_mmse[i][j]; }
						if (MIMODET == "MMSE") { llr_mmse_bit(Nt, y_mmse_tmp, oriLLR_tmp, sigma_sq, ModType); }
						if (MIMODET == "KBEST")
						{
							for (int i = 0; i < len_ori_llr_tmp; i++) { oriLLR_tmp[i] = oriLLR_kbest[j][i]; }
							// copy LLR
							for (int i = 0; i < Nv; i++)
							{
								oriLLRh[i][j] = oriLLR_tmp[i];
								upLLR[i][j] = oriLLR_tmp[i];
							}
							// First round of vertical polar decoding
							PolarDecodeIterVertical(oriLLRh, upLLR, umid, GGv, uout_vertical_partial, A_Acv, Av, Nh, Nv, j, 1);
							// test: oriLLRh
							/*int num_err_dec = 0;
							cout << "Dec from oriLLRh" << endl;
							for (int i = 0; i < Nv; i++)
							{
								cout << endl;
								for (int j = 0; j < Nh; j++)
								{
									cout << setw(2) << (oriLLRh[i][j] < 0);
									num_err_dec += abs(x[i][j] - (oriLLRh[i][j] < 0));
								}
							}
							cout << endl;
							cout << "Dec from oriLLRh" << endl << k << endl;
							cout << "Detector output error:  " << num_err_dec << endl;*/
							// test: oriLLRh
							// test: uout_vertical_partial
							/*cout << "uout_vertical_partial" << endl;
							cout << "k= " << k << endl;
							for (int i = 0; i < Nv; i++)
							{
								cout << endl;
								for (int j = 0; j < Nh; j++)
								{
									cout << uout_vertical_partial[i][j];
								}
							}*/
							// test: uout_vertical_partial
						}
						for (int i = 0; i < k_kbest_extend; i++) delete[] paths_kbest[i];
						delete[] paths_kbest;
						for (int i = 0; i < Nv; i++) delete[] llr_rearranged_array[i];
						delete[] llr_rearranged_array;
						delete[] ED;
						delete[] llr_rearranged_tmp;
						delete[] llr_ordered_tmp;
						delete[] llr_to_detector;
					}
					// Temporal modulation
					// test for llr to detector
					/*cout << "IDD Dec" << endl;
					for (int i = 0; i < Nv; i++)
					{
						cout << endl;
						for (int j = 0; j < Nh; j++)
							cout << setw(2) << (llr_to_det_test[i][j]< 0);
					}
					cout << "IDD Dec" << endl;
					cout << "IDD LLR" << endl;
					for (int i = 0; i < Nv; i++)
					{
						cout << endl;
						for (int j = 0; j < Nh; j++)
							cout << setw(10) << llr_to_det_test[i][j] ;
					}
					cout << "IDD LLR" << endl;
					for (int i = 0; i < Nv; i++) delete[] llr_to_det_test[i];
					delete[] llr_to_det_test;*/
					// test for llr to detector

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
			for (int iter = 1; iter <= softiter; iter++)
			{
				//	cout << "Vertical Start:\n";
					//vertical decoding
				//cout << "Horizon Start:\n";
				//horizontal decoding
				// #
//#pragma omp parallel for
				for (int decodingi = 0; decodingi < Nh; decodingi++) // 291.70kb
				{
					if (iter == 1) break;
					//for SCL Decoding
					int Counth = 0, Countinfoh = 0;
					float Qsumunh = 0;
					int** sum = new int* [L];
					//for (int i = 0; i < L; i++)  sum[i] = new  int[Nmax];
					int** u1 = new int* [L];
					for (int i = 0; i < L; i++)  u1[i] = new _MM_ALIGN16 int[Nmax];
					for (int i = 0; i < L; i++)  sum[i] = new _MM_ALIGN16 int[Nmax];
					//for (int i = 0; i < L; i++)  u1[i] = new  int[Nmax];

					float* LLR_in = new float[L];
					float* W = new float[2 * L];
					int* SCL_index = new int[L];
					int* better = new int[L];
					int* worse = new int[L];
					int* indexpm = new int[L];
					float* PM = new float[L];
					float** LLR = new float* [L];
					for (int i = 0; i < L; i++)  LLR[i] = new _MM_ALIGN16 float[2 * Nmax + 2];
					// for (int i = 0; i < L; i++)  LLR[i] = new float[2 * Nmax + 2];
					////////////////////////////////
					for (int j = 0; j < Nv; j++)
						for (int k = 0; k < L; k++)
							LLR[k][j] = upLLR[j][decodingi];
					//cout << LLR[0][j] << " ";

				//cout << endl;
					int index = 0;
					/*			if (L == 1)
								{
									//T_start = clock();
									Count = 0;
									decode(*LLR, *sum, (int)(log(Nv) / log(2)), A_Acv, flag, Nv, Kv);
									//PolarEncode_xor(*u1, *sum, N);
									//T_finish = clock();
								}
								else
								{*/
					for (int k = 0; k < L; k++) { PM[k] = 0; indexpm[k] = k; }
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
					if (adaptive_beta)
					{
						if (softmethod == "Pyndiah") Pyndiahv_adaptivebeta(LLR, sum, oriLLRh, indexpm, iter * 2 - 2, decodingi, upLLR);
						if (softmethod == "MITSO") MITSOv(LLR, PM, sum, oriLLRh, indexpm, iter * 2 - 2, decodingi, upLLR, Qsumunh);
					}
					else
					{
						if (softmethod == "Pyndiah") Pyndiahv(LLR, sum, oriLLRh, indexpm, iter * 2 - 2, decodingi, upLLR);
						if (softmethod == "MITSO") MITSOv(LLR, PM, sum, oriLLRh, indexpm, iter * 2 - 2, decodingi, upLLR, Qsumunh);
					}
					// half iteration Enabled
					if (half_iter && iter == softiter)
					{
						if (iter == softiter) PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), Nh);
						for (int j = 0; j < Kv; j++)
							umid[Av[j]][decodingi] = u1[index][Av[j]];
					}
					// half iteration
					////////delete SCL variables
					for (int i = 0; i < L; i++)
						delete[]LLR[i];
					delete[] LLR;
					for (int i = 0; i < L; i++)
						delete[]sum[i];
					for (int i = 0; i < L; i++)
						delete[]u1[i];
					delete[] u1;
					delete[] sum;
					delete[] PM;
					delete[] LLR_in;
					delete[] W;
					delete[] SCL_index;
					delete[] better;
					delete[] worse;
					delete[] indexpm;

				}
				if (weirditer) { continue; }

#pragma omp parallel for
				for (int decodingi = 0; decodingi < Nv; decodingi++)
				{
					// 293.70kb
					// cout << "KKKK" << endl;
					if (half_iter && (iter == softiter)) { break; }
					string polarde_now = polarde_array[iter - 1];
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
					if (polarde_now == "SCL")
					{
						for (int k = 0; k < L; k++)  PM[k] = 0, indexpm[k] = k;
						//T_start = clock();
						decode_list(LLR, sum, (int)(log(Nh) / log(2)), A_Ach, 0, Nh, Kh, PM, L, (int)(log(L) / log(2)), LLR_in, W, SCL_index, better, worse, &Countv, &Countinfov, &Qsumunv);
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
						if (iter == softiter) PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), Nh);
						index = indexpm[0];
						//}
						//calc distance s()
						if (iter < softiter)
						{
							if (adaptive_beta)
							{
								if (softmethod == "Pyndiah") Pyndiahh_adaptivebeta(LLR, sum, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR);
								if (softmethod == "MITSO") MITSOh(LLR, PM, sum, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR, Qsumunv);
								//	getchar();
							}
							else
							{
								if (softmethod == "Pyndiah") Pyndiahh(LLR, sum, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR);
								if (softmethod == "MITSO") MITSOh(LLR, PM, sum, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR, Qsumunv);
								//	getchar();
							}
						}
						else
							for (int j = 0; j < Kh; j++)
								umid[decodingi][Ah[j]] = u1[index][Ah[j]];
					}

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
			int calcsum[Nv], calcu1[Nv];
			int calcu1h[Nh];
			// calcsum is from the direct output of the decoder
			// if output is the result of column decoding (meaning there is another row decoding to be done)
			if (half_iter)
			{
				for (int i = 0; i < Kv; i++)
				{
					PolarEncode_xor(calcu1h, umid[Av[i]], Nh);
					for (int j = 0; j < Kh; j++)
						uout[Av[i]][Ah[j]] = calcu1h[Ah[j]];
				}
			}
			// if output is the result of column decoding (meaning that there is another row decoding to be done)
			else
			{
				for (int i = 0; i < Kh; i++)
				{
					for (int j = 0; j < Nv; j++)
						calcsum[j] = umid[j][Ah[i]];
					PolarEncode_xor(calcu1, calcsum, Nv);
					for (int j = 0; j < Kv; j++)
						uout[Av[j]][Ah[i]] = calcu1[Av[j]];
				}
			}
			// test: final output
			/*cout << "final output" << endl;
			for (int k = 0; k < Nh; k++) { cout << uout[Av[0]][k] << " "; }
			cout << endl;*/
			// test: final output
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
			if (Temp_Error < Num_Error)
			{
				Num_Frame_Error++;
			}
			if (frame % 10 == 0)
			{
				/*printf("\r<%d> SNR= %.2f FER= %d / %d = %f BER= %d / %d = %lf E-Det = %f, E-MRB = %f"  ,
					frame, nSNR, Num_Frame_Error, frame, (double)Num_Frame_Error / frame, Num_Error, (N_Info * frame), (double)Num_Error / N_Info / frame,
					(double)n_det_err/total_bit_cnt, (double)n_det_err_mrb/(total_bit_cnt/(ModType/2)));*/
				printf("\r<%d> SNR= %.2f FER= %d / %d = %f BER= %d / %d = %lf ",
					frame, nSNR, Num_Frame_Error, frame, (double)Num_Frame_Error / frame, Num_Error, (N_Info * frame), (double)Num_Error / N_Info / frame);
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
		// save to .txt file
		string polarde_serial = "";
		for (int i = 0; i < softiter; i++)
			polarde_serial = polarde_serial + polarde_array[i] + string("_");
		string file_name = string("ResFinal\\2D_LSD_D_") + polarde_serial + MIMODET + string("T") + to_string(Nt) + string("R") + to_string(Nr)
			+ string("Q") + to_string(int(pow(2, ModType))) + string("L") + to_string(Nh) + string("K") + to_string(Kh)
			+ string("T") + to_string(softiter) + string("H") + to_string(half_iter) + string("UB") + to_string(usebeta)
			+ string(".txt");
		save_array_txt(nSNR, BER_array[SNR_count], FER_array[SNR_count], file_name);
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
		delete[] uout_vertical_partial[i];
	delete[] uout_vertical_partial;
	for (int i = 0; i < Nv; i++)
		delete[] y[i];
	for (int i = 0; i < Nh; i++)
		delete[] oriLLR_kbest[i];
	delete[] oriLLR_kbest;
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
	for (int i = 0; i < mimo_blk_height; i++) {
		delete[] x_mod[i];
	}
	delete[] x_mod;
	for (int i = 0; i < Nr; i++) {
		delete[] y_mod[i];
	}
	delete[] y_mod;
	for (int i = 0; i < 2 * mimo_blk_height; i++) {
		delete[] y_mmse[i];
	}
	delete[] y_mmse;
	for (int i = 0; i < 2 * Nr; i++) {
		delete[] y_mod_real[i];
	}
	for (int i = 0; i < mimo_blk_length; i++)
	{
		delete[] y_mod_real_temp_omp[i];
	}
	delete[] y_mod_real_temp_omp;
	for (int i = 0; i < mimo_blk_length; i++)
	{
		delete[] y_mmse_temp_omp[i];
	}
	delete[] y_mmse_temp_omp;
	delete[] y_mod_real;
	delete[] x_mod_tmp;
	delete[] y_mod_tmp;
	delete[] y_mod_real_tmp;
	delete[] x_tmp;
	delete[] oriLLR_tmp;

}


void system2D_LSD_epdet_sys(double* BER_array, double* FER_array)
{
	std::ofstream file_total("total_timing.txt", std::ios::app);
	std::ofstream file_det("det_timing.txt", std::ios::app);
	std::ofstream file_dec("dec_timing.txt", std::ios::app);
	std::ofstream file_err_num("total_frames.txt", std::ios::app);
	auto det_end = std::chrono::system_clock::now();
	int Nt2 = 2 * Nt;
	int Nr2 = 2 * Nr;

	/*if (!file.is_open()) {
		std::cerr << "Error opening file: " << file_name << std::endl;
		return;
	}*/

	//file << setw(20) << SNR << setw(20) << BER << setw(20) << FER << "\n";

	//freopen("result.txt", "w", stdout);
	printsetting();
	cout << "\nSNR	 " << "FER      BER" << endl;
	int N_Info = Kv * Kh - CRCL, Nmax;
	//Interleaving pattern
	/*std::vector<int> intrlv_pattern(Nv);
	std::vector<int>* intrlv_pattern_2D = new vector<int>[Nh];
	for (int i = 0; i < Nh; i++) intrlv_pattern_2D[i].resize(Nv);*/
	int intrlv_pattern_num = Nh / ModType;
	std::vector<int> intrlv_pattern(Nv);
	std::vector<int>* intrlv_pattern_2D = new vector<int>[Nh];
	std::vector<int>* intrlv_pattern_2D_time = new vector<int>[Nv];
	std::vector<int>* intrlv_pattern_2D_inner = new vector<int>[Nh];
	std::vector<int>** intrlv_pattern_2D_time_partial = new vector<int> *[Nv];
	for (int i = 0; i < Nv; i++) intrlv_pattern_2D_time_partial[i] = new vector<int>[intrlv_pattern_num];
	for (int i = 0; i < Nv; i++)
	{
		for (int j = 0; j < intrlv_pattern_num; j++)
		{
			intrlv_pattern_2D_time_partial[i][j].resize(ModType);
		}
	}
	for (int i = 0; i < Nh; i++) intrlv_pattern_2D[i].resize(Nv);
	for (int i = 0; i < Nh; i++) intrlv_pattern_2D_inner[i].resize(Kv);
	for (int i = 0; i < Nv; i++) intrlv_pattern_2D_time[i].resize(Nh);

	float CodeRate = (float)N_Info / (Nh * Nv);
	Nmax = max(Nh, Nv);
	// num of symbols per vertical codeword
	int Nsym = 0;
	if (MOD_DIM == "Spatial") { Nsym = Nv / ModType; }
	else if (MOD_DIM == "Temporal") { Nsym = Nh / ModType; }
	int symbol_length = ModType;

	// Size of a mimo block is given by height(spatial) * length(temporal)
	int mimo_blk_height = 0;
	int mimo_blk_length = 0;
	if (MOD_DIM == "Spatial")
	{
		mimo_blk_height = Nsym;
		mimo_blk_length = Nh;
	}
	else if (MOD_DIM == "Temporal")
	{
		mimo_blk_height = Nv;
		mimo_blk_length = Nsym;
	}
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

	// Arrays for MIMO
	int** x_mmse_out = new int* [Nv]; // delete!
	for (int i = 0; i < Nv; i++) { x_mmse_out[i] = new int[Nh]; }

	//std::complex<float>** x_mod = new std::complex<float>* [Nsym]; // modulated symbols (spatial modulation) DELETE!
	//for (int i = 0; i < Nsym; i++) { x_mod[i] = new std::complex<float>[Nh]; }
	//std::complex<float>** x_mod_temporal = new std::complex<float>*[Nv]; // modulated symbols (temporal modulation) DELETE!
	//for (int i = 0; i < Nv; i++) x_mod_temporal[i] = new std::complex<float>[Nsym];
	//std::complex<float>** y_mod = new std::complex<float>* [Nr]; // Output symbols after AWGN channels // delete!
	//for (int i = 0; i < Nr; i++) { y_mod[i] = new std::complex<float>[Nh]; }
	//float** y_mod_real = new float* [2 * Nr];
	//for (int i = 0; i < 2*Nr; i++) y_mod_real[i] = new float[Nh];
	//float** y = new float* [Nv];// Output Signal After AWGN Channel
	//for (int i = 0; i < Nv; i++) { y[i] = new float[Nh]; }
	//float** y_mmse = new float* [2*Nsym];
	//for (int i = 0; i < 2*Nsym; i++) { y_mmse[i] = new float[Nh]; }  // MMSE detector output // delete!
	//float* y_mmse_tmp = new float[2*Nsym];
	//std::complex<float>* x_mod_tmp = new std::complex<float>[Nsym]; // modulated symbols of a single column // delete!
	//std::complex<float>* y_mod_tmp = new std::complex<float>[Nr]; // received symbols of a single column // delete!
	//float* y_mod_real_tmp = new float[2*Nr];
	//// float* y_mod_tmp = new float[Nr];  // delete!


	std::complex<float>** x_mod = new std::complex<float>*[mimo_blk_height]; // modulated symbols (spatial modulation) DELETE!
	for (int i = 0; i < mimo_blk_height; i++) { x_mod[i] = new std::complex<float>[mimo_blk_length]; }
	std::complex<float>** y_mod = new std::complex<float>*[Nr]; // Output symbols after AWGN channels // delete!
	for (int i = 0; i < Nr; i++) { y_mod[i] = new std::complex<float>[mimo_blk_length]; }
	float** y_mod_real = new float* [2 * Nr];
	for (int i = 0; i < 2 * Nr; i++) y_mod_real[i] = new float[mimo_blk_length];
	float** y = new float* [Nv];// Output Signal After AWGN Channel
	for (int i = 0; i < Nv; i++) { y[i] = new float[Nh]; }
	float** y_mmse = new float* [2 * mimo_blk_height];
	for (int i = 0; i < 2 * mimo_blk_height; i++) { y_mmse[i] = new float[mimo_blk_length]; }  // MMSE detector output // delete!
	float* y_mmse_tmp = new float[2 * mimo_blk_height];
	std::complex<float>* x_mod_tmp = new std::complex<float>[mimo_blk_height]; // modulated symbols of a single column // delete!
	std::complex<float>* y_mod_tmp = new std::complex<float>[Nr]; // received symbols of a single column // delete!
	float* y_mod_real_tmp = new float[2 * Nr];
	// float* y_mod_tmp = new float[Nr];  // delete!
	float** y_mod_real_temp_omp = new float* [mimo_blk_length]; // new delete!
	for (int i = 0; i < mimo_blk_length; i++) { y_mod_real_temp_omp[i] = new float[2 * Nr]; }
	float** y_mmse_temp_omp = new float* [mimo_blk_length]; // new delete!
	for (int i = 0; i < mimo_blk_length; i++) { y_mmse_temp_omp[i] = new float[2 * mimo_blk_height]; }
	// ep
	float** ep_extr_mean_omp = new float* [mimo_blk_length];
	for (int i = 0; i < mimo_blk_length; i++) { ep_extr_mean_omp[i] = new float[2 * mimo_blk_height]; }
	float** ep_extr_var_omp = new float* [mimo_blk_length];
	for (int i = 0; i < mimo_blk_length; i++) { ep_extr_var_omp[i] = new float[2 * mimo_blk_height]; }
	int* x_tmp = new int[Nv]; // delete!

	// paths selected by k-best detection


	// CRC bits
	unsigned char* u_crc = new unsigned char[Kv * Kh];

	float** upLLR = new float* [Nv];
	for (int i = 0; i < Nv; i++) { upLLR[i] = new _MM_ALIGN16 float[Nh]; }
	//float** oriLLRv = new float* [Nh];
	//for (int i = 0; i < Nh; i++) { oriLLRv[i] = new _MM_ALIGN16 float[2 * Nv + 2]; }
	float** oriLLRh = new float* [Nv];
	for (int i = 0; i < Nv; i++) { oriLLRh[i] = new _MM_ALIGN16 float[Nh]; }
	int len_ori_llr_tmp = 0;
	if (MOD_DIM == "Spatial") { len_ori_llr_tmp = Nv; }
	else if (MOD_DIM == "Temporal") { len_ori_llr_tmp = ModType * Nt; }
	float* oriLLR_tmp = new _MM_ALIGN16 float[len_ori_llr_tmp];  // delete!
	//float** oriLLR_kbest = new float* [Nh];
	float** oriLLR_kbest = new float* [mimo_blk_length];
	for (int i = 0; i < Nh; i++) { oriLLR_kbest[i] = new float[len_ori_llr_tmp]; }
	double* sum_sdis_array = new double[Nv];
	for (int i = 0; i < Nv; i++) sum_sdis_array[i] = 0;

	//*****************************************************************
	// EP detector
	float* Hreal_array = new float[(2 * Nr) * (2 * Nt)]; // delete
	float* HTH_array = new float[(2 * Nt) * (2 * Nt)]; // delete
	float* HTY_array = new float[2 * Nt]; // delete 





	int** GGh = new int* [Nh];
	int** GGv = new int* [Nv];


	MatrixXf Gh = Fn(Nh);  // Generate Matrix of horizontal coding
	MatrixXf Gv = Fn(Nv);  // Generate Matrix of vertical coding
	MatrixXcf H(Nr, Nt); // Complex channel Matrix H: Nr x Nt
	//Matrix<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic> H = MatrixXcf::Random(Nr, Nt); // Real domain channel matrix
	//// //MatrixX Definitions
	MatrixXf H_real(2 * Nr, 2 * Nt); // Real domain channel matrix
	MatrixXf received_signal(2 * Nr, 1);
	MatrixXf HTH(2 * Nt, 2 * Nt);
	MatrixXf Identity = MatrixXf::Identity(2 * Nt, 2 * Nt);
	MatrixXf output_signal(2 * Nt, 1);
	MatrixXf mmse_matrix(2 * Nt, 2 * Nr);
	MatrixXf HTHIinv(2 * Nt, 2 * Nt);

	// for KBest
	MatrixXf R(2 * Nt, 2 * Nt);
	MatrixXf Q(2 * Nr, 2 * Nt);
	MatrixXf z(2 * Nt, 1);



	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> H_real = MatrixXf::Random(2 * Nr, 2 * Nt); // Real domain channel matrix
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> received_signal = MatrixXf::Random(2 * Nr,1);
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> HTH = MatrixXf::Random(2 * Nt, 2 * Nt);
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> Identity = MatrixXf::Random(2 * Nt, 2 * Nt);
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> output_signal = MatrixXf::Random(2 * Nt, 1);
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> mmse_matrix = MatrixXf::Random(2 * Nt, 2 * Nr);


	// Matrix Definitions
	//Eigen::Matrix<float, 2 * Nr, 2 * Nt> H_real; // ʵ���ŵ�����
	//Eigen::Matrix<float, 2 * Nr, 1> received_signal; // �����źž���
	//Eigen::Matrix<float, 2 * Nt, 2 * Nt> HTH; // HTH����
	//Eigen::Matrix<float, 2 * Nt, 2 * Nt> Identity = Eigen::Matrix<float, 2 * Nt, 2 * Nt>::Identity(); // ��λ����
	//Eigen::Matrix<float, 2 * Nt, 1> output_signal; // ����źž���
	//Eigen::Matrix<float, 2 * Nt, 2 * Nr> mmse_matrix; // MMSE����


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
	// srand((unsigned int)time(0));
	std::random_device rd;              // ��ȡ���������
	std::mt19937 gen(rd());             // ʹ�����ӳ�ʼ��mt19937������
	std::uniform_int_distribution<> distrib(1, 1000000);
	int Num_Error, Num_Frame_Error;
	float Decoding_Duration;

	int Num_Error_new, Num_Frame_Error_new;
	float Decoding_Duration_new;

	clock_t T_start, T_finish;
	int SNR_count = 0;
	int anssc, anssd;
	for (double nSNR = Min_SNR; nSNR <= Max_SNR; nSNR += SNR_step)
	{
		for (int i = 0; i < Nv; i++) sum_sdis_array[i] = 0;
		double sdis_sum = 0;
		// TEST: det err
		int n_det_err = 0; // n error bits at detector output
		int n_det_err_mrb = 0; // error bits at detector output(MRBs)
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
		// test: Horizontal frozen bits

		/*cout << "Frozen Bits" <<endl;
		for (int i = 0; i < Nh - Kh; i++)
		{
			cout << setw(5) << Acv[i];
		}
		cout << endl;*/
		// test: Horizontal frozen bits
		// Frame Loop
		int passnum;
		int frame = 0;

		// generate interleaving patterns
		if (flag_interleave)
		{
			string ofile_intrlv_name = string("intrlv_pattern_71_frame_49263.txt");
			string ofile_spatial_intrlv_name = string("intrlv_pattern_spatial71.txt");
			ifstream ofile_intrlv(ofile_intrlv_name);
			ifstream ofile_spatial_intrlv(ofile_spatial_intrlv_name);
			int file_in = 1;
			for (int i = 0; i < Nv; i++)
			{
				for (int j = 0; j < intrlv_pattern_num; j++)
				{

					for (int k = 0; k < ModType; k++)
						// ofile_intrlv >> file_in;
						ofile_intrlv >> intrlv_pattern_2D_time_partial[i][j][k];
					//cout << endl;
				}
			}
			ofile_intrlv.close();

			for (int i = 0; i < Nh; i++)
			{
				for (int j = 0; j < Nv; j++)
					ofile_spatial_intrlv >> intrlv_pattern_2D[i][j];
			}
			ofile_spatial_intrlv.close();
			/*for (int i = 0; i < Nh; i++)
			{
				for (int j = 0; j < Nv; j++)
					cout << intrlv_pattern_2D[i][j] << endl;
			}*/
			generateIntrlvPattern(intrlv_pattern, Nv);
			for (int i = 0; i < Nh; i++) { generateIntrlvPattern(intrlv_pattern_2D[i], Nv); }
			for (int i = 0; i < Nv; i++) { generateIntrlvPattern(intrlv_pattern_2D_time[i], Nh); }
			if (flag_interleave_partial)
			{
				for (int i = 0; i < Nv; i++)
				{
					for (int j = 0; j < intrlv_pattern_num; j++)
					{
						generateIntrlvPattern(intrlv_pattern_2D_time_partial[i][j], ModType);
						/*for (int k = 0; k < ModType; k++)
							cout << intrlv_pattern_2D_time_partial[i][j][k];
						cout << endl;*/
					}
				}
			}
			for (int i = 0; i < Nh; i++) { generateIntrlvPattern(intrlv_pattern_2D_inner[i], Kv); }
		}

		//ofstream ofile_intrlv("intrlv_pattern.txt");
		//for (int i=0; i<Nv; i++)
		//{
		//	for (int j = 0; j < intrlv_pattern_num; j+for+)
		//	{
		//		
		//		for (int k = 0; k < ModType; k++)
		//			ofile_intrlv << intrlv_pattern_2D_time_partial[i][j][k] << "   ";
		//		//cout << endl;
		//	}
		//}
		//ofile_intrlv.close();
		int frame_tmp = 0;
		int frame_delta = 0;
		int n_frame_error_pre = 0;
		while (Num_Frame_Error < errframe)
		{
			/*test: processing time*/
			//auto start = std::chrono::system_clock::now();
			/*if (frame == 10) {
				cout << "AAAAAAAAAA" << endl;
			}*/
			frame++;
			// generate interleaving pattern
			/*if (flag_interleave)
			{
				generateIntrlvPattern(intrlv_pattern, Nv);
				for (int i = 0; i < Nh; i++) { generateIntrlvPattern(intrlv_pattern_2D[i], Nv); }
				for (int i = 0; i < Nv; i++) { generateIntrlvPattern(intrlv_pattern_2D_time[i], Nh); }
				if (flag_interleave_partial)
				{
					for (int i = 0; i < Nv; i++)
					{
						for (int j = 0; j < intrlv_pattern_num; j++)
						{
							generateIntrlvPattern(intrlv_pattern_2D_time_partial[i][j], ModType);
						}
					}
				}
			}*/

			//cout << frame << endl;
			int u_input[16] = { 0,     1     ,0,     1,
									0,     1,     0,     0,
									0 ,    1,     0,     0,
									1 ,    0,     0,     0 };
			for (int i = 0; i < N_Info; i++)
			{
				a_data[i] = (unsigned char)(distrib(gen) % 2);
				//		a_data[i] = u_input[i];
			}

			//cout << endl;
			if (CRCL) tx_append_crc(a_data, N_Info, CRCL);
			for (int i = 0; i < Nv; i++)
				for (int j = 0; j < Nh; j++)
					u[i][j] = 0;
			//omp_set_num_threads(4);
//#pragma omp parallel for
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
			// test:original info
			//cout << "original info u:" << endl;
			////for (int j = 0; j < Nh; j++) { cout << u[Av[0]][j] << " "; }
			//for (int i = 0; i < Nv; i++)
			//{
			//	cout << endl;
			//	for (int j = 0; j < Nh; j++)
			//		cout << u[i][j];
			//}
			//cout << endl;
			// Polar Encoding
//#pragma omp parallel for

			// test: original info
			/*cout << "original info u:" << endl;
			display_array(u, Nv, Nh);*/
			// test: original info

			for (int i = 0; i < Nv; i++)
				PolarEncode_xor(midx[i], u[i], Nh);
			// test:horizontal encoding
			//cout << "u after horizontal encoding:" << endl;
			////for (int j = 0; j < Nh; j++) { cout << midx[Av[0]][j] << " "; }
			//for (int i = 0; i < Nv; i++)
			//{
			//	cout << endl;
			//	for (int j = 0; j < Nh; j++)
			//		cout << midx[i][j];
			//}
			//cout << endl;
			//cout << endl;
			// test:horizontal encoding
			// 
			// column-wise interleaving

			// test: row encoding
			/*cout << "u after horizontal encoding:" << endl;
			display_array(midx, Nv, Nh);*/
			// test: row encoding
			if (flag_interleave_inner)
			{
				int* tmp_to_interleave = new int[Nv];
				for (int i = 0; i < Nh; i++)
				{
					for (int j = 0; j < Nv; j++) { tmp_to_interleave[j] = midx[j][i]; }
					randIntrlvPartial(intrlv_pattern_2D_inner[i], tmp_to_interleave, Av, Nv, Kv);
					for (int j = 0; j < Nv; j++) { midx[j][i] = tmp_to_interleave[j]; }
				}
				delete[] tmp_to_interleave;
			}

			// test: innner interleaving
			/*cout << "u after inner interleaving:" << endl;
			display_array(midx, Nv, Nh);*/
			// test: innner interleaving
			// 
			// systematic horizontal encoding
			for (int i = 0; i < Nv; i++)
			{
				int* code_tmp = new int[Nh];
				// TEST:
				/*if (i == Av[0])
				{
					cout << "After Stage 0 Encoding(Horizontal encoding, before setting to 0) \n";
					for (int i = 0; i < Nh; i++) { cout << " " << midx[Av[0]][i]; }
					cout << endl;
				}*/
				// TEST
				for (int j = 0; j < Nh - Kh; j++) { midx[i][Ach[j]] = 0; }
				// TEST:
				/*if (i == 0)
				{
					cout << "After Stage 1 Encoding(Setting to 0) \n";
					for (int i = 0; i < Nh; i++) { cout << " " << midx[Av[0]][i]; }
					cout << endl;
				}*/
				// TEST
				for (int j = 0; j < Nh; j++) { code_tmp[j] = midx[i][j]; }
				PolarEncode_xor(midx[i], code_tmp, Nh);
				delete[] code_tmp;
			}
//#pragma omp parallel for
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

			// test: vertical encoding
			/*cout << "u after vertical encoding:" << endl;
			display_array(x, Nv, Nh);*/
			// test: vertical encoding
			/*for (int i = 0; i < 2 * Nr; i++)
				for (int j = 0; j < 2 * Nt; j++)
					H_real(i, j) = 0;*/
					/*MatrixXf H_real_test(2 * Nr, 2 * Nt);
					MatrixXf M = H_real_test.transpose();
					MatrixXf N = H_real_test;
					HTHIinv = M * N;*/
					//HTHIinv = M * H_real;
					// test:vertical encoding
					//cout << "u after vertical encoding:" << endl;
					////for (int j = 0; j < Nh; j++) { cout << x[Av[0]][j] << " "; }
					//for (int i = 0; i < Nv; i++)
					//{
					//	cout << endl;
					//	for (int j = 0; j < Nh; j++)
					//		cout << x[i][j];
					//}
					//cout << endl;
					//cout << endl;
					// test:vertical encoding
			if (atten == 2)
			{
				// modulation
				// #pragma omp parallel for
				// Interleaving
				if (flag_interleave)
				{
					if (flag_interleave_all)  // interleave temporally
					{
						for (int i = 0; i < Nv; i++)
						{
							randIntrlv(intrlv_pattern_2D_time[i], x[i], Nh);
						}
					}
					//// test: to be interleaved
					//cout << "x before interleaving:" << endl;
					//// print x line by line, add spaces between each 4 terms
					//for (int i = 0; i < Nv; i++)
					//{
					//	for (int j = 0; j < Nh; j++)
					//	{
					//		if (!(j % 4)) cout << " ";
					//		cout << x[i][j] << " ";
					//	}
					//	cout << endl;
					//}

					if (flag_interleave_partial)
					{
						for (int i = 0; i < Nv; i++)
						{
							for (int j = 0; j < intrlv_pattern_num; j++)
							{
								randIntrlv(intrlv_pattern_2D_time_partial[i][j], x[i] + j * ModType, ModType);
							}
						}
					}

					//// test: to be interleaved
					//cout << "x after interleaving:" << endl;
					//// print x line by line, add spaces between each 4 terms
					//for (int i = 0; i < Nv; i++)
					//{
					//	for (int j = 0; j < Nh; j++)
					//	{
					//		if (!(j % 4)) cout << " ";
					//		cout << x[i][j] << " ";
					//	}
					//	cout << endl;
					//}
					if ((!flag_interleave_block) && flag_interleave_spatial)
					{
						for (int j = 0; j < Nh; j++)
						{
							for (int i = 0; i < Nv; i++) { x_tmp[i] = x[i][j]; }
							randIntrlv(intrlv_pattern_2D[j], x_tmp, Nv);
							for (int i = 0; i < Nv; i++) { x[i][j] = x_tmp[i]; }
						}
					}


					if (flag_interleave_block)
					{
						randIntrlvBlock(x, Nh, Nv, ModType, interleave_rows);
					}
				}
				if (MOD_DIM == "Spatial")
				{
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
						/*for (int i = 0; i < mimo_blk_height; i++)
						{
							cout << x_mod_tmp[i] << endl;
						}
						cout << "111111" << endl;*/
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
					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < 2 * Nr; j++) { y_mod_real_temp_omp[i][j] = y_mod_real[j][i]; }
					}
					// for KBest Detection
					Eigen::HouseholderQR<Eigen::MatrixXf> qr_Hreal(H_real);
					Q = qr_Hreal.householderQ();
					R = qr_Hreal.matrixQR().triangularView<Eigen::Upper>();
					//omp_set_num_threads(8);

#pragma omp parallel for private(z)
					for (int j = 0; j < Nh; j++)
					{
						MatrixXf received_signal_omp(2 * Nr, 1);
						int k_kbest_extend = pow(2, ModType / 2) * K_KBEST;
						float** paths_kbest = new float* [k_kbest_extend];
						for (int i = 0; i < k_kbest_extend; i++) paths_kbest[i] = new float[2 * Nt];
						float* ED = new float[K_KBEST];

						if (MIMODET == "MMSE") { MMSEdetection(H_real, Identity, received_signal, mmse_matrix, output_signal, HTH, y_mod_real_temp_omp[j], y_mmse_temp_omp[j], sigma_sq); }
						// ִ��һЩ����

						for (int i = 0; i < 2 * Nr; i++) { received_signal_omp(i, 0) = y_mod_real_temp_omp[j][i]; }
						z = Q.transpose() * received_signal_omp;
						if (MIMODET == "KBEST")
						{
							KBestDetection(R, z, H_real, y_mod_real_temp_omp[j], ED, y_mmse_temp_omp[j], paths_kbest, K_KBEST, ModType, Nt, Nr);
							llr_kbest_bit(Nt, y_mmse_temp_omp[j], oriLLR_kbest[j], sigma_sq, ModType, K_KBEST, paths_kbest, ED);
						}
						for (int i_path = 0; i_path < k_kbest_extend; i_path++) { delete[] paths_kbest[i_path]; }

						delete[] paths_kbest;
						delete[] ED;
					}

					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < 2 * Nt; j++) { y_mmse[j][i] = y_mmse_temp_omp[i][j]; }
					}
					// sym llr to bitwise llr
					// Detection is done for every column, that is, for every Nv bits
					// Interleaving is not currently considered
					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < 2 * Nt; j++) { y_mmse_tmp[j] = y_mmse[j][i]; }
						if (MIMODET == "MMSE") { llr_mmse_bit(Nt, y_mmse_tmp, oriLLR_tmp, sigma_sq, ModType); }
						if (MIMODET == "KBEST") { for (int j = 0; j < len_ori_llr_tmp; j++) { oriLLR_tmp[j] = oriLLR_kbest[i][j]; } }
						if (flag_interleave)
						{
							randDeintrlv(intrlv_pattern, oriLLR_tmp, Nv);
						}
						bool flagIsMax = false;
						for (int j = 0; j < Nv; j++)
						{
							if (!(j % 2))
							{
								if (abs(x_mmse_out[j][i]) > abs(x_mmse_out[j][i + 1])) { flagIsMax = true; }
								else { flagIsMax = false; }
							}
							else
							{
								if (abs(x_mmse_out[j][i]) > abs(x_mmse_out[j][i - 1])) { flagIsMax = true; }
								else { flagIsMax = false; }
							}
							oriLLRh[j][i] = oriLLR_tmp[j];
							upLLR[j][i] = oriLLR_tmp[j];
							// TEST�� det error
							x_mmse_out[j][i] = 0;
							if (oriLLR_tmp[j] < 0) { x_mmse_out[j][i] = 1; }
							if (x_mmse_out[j][i] != x[j][i]) {
								n_det_err++;
								// error at max-reliable bit. Warning: Only fit for 16QAM modulation
								if (flagIsMax) { n_det_err_mrb++; }
							}
							total_bit_cnt++;
						}
					}
					// test : detection TIME
					/*det_end = std::chrono::system_clock::now();
					auto duration_detection = std::chrono::duration_cast<std::chrono::microseconds>(det_end - start);
					file_det << duration_detection.count() << endl;*/
				}
				// Temporal modulation
				else if (MOD_DIM == "Temporal")
				{
					//modulation
					//cout << "START TEST" << endl;
					for (int i = 0; i < Nv; i++) // every row
					{
						modulation(x[i], x_mod[i], ModType, mimo_blk_length);
					}
					/*for (int i = 0; i < mimo_blk_height; i++)
					{
						for (int j = 0; j < mimo_blk_length; j++)
							cout << x_mod[i][j] << endl;
					}
					cout << "111111" << endl;*/
					//AWGN channel (received y_mod)
					// test: modulated symbols
					/*cout << endl;
					cout << "modulated symbols" << endl;
					for (int i = 0; i < Nt; i++) cout << x_mod[i][0] << endl;
					cout << "modulated symbols" << endl;*/
					// test: modulated symbols
					H = generateChannelMatrix(Nr, Nt);
					for (int j = 0; j < Nsym; j++)
					{
						// Extracting a single column
						for (int i = 0; i < Nv; i++) { x_mod_tmp[i] = x_mod[i][j]; }
						AWGN(Nr, Nt, x_mod_tmp, y_mod_tmp, H, sigma_sq);
						for (int i = 0; i < Nr; i++) { y_mod[i][j] = y_mod_tmp[i]; }
					}
					convertHToReal(H, H_real);

					for (int j = 0; j < Nsym; j++)
					{
						for (int i = 0; i < Nr; i++) { y_mod_tmp[i] = y_mod[i][j]; }
						convertYToReal(Nr, y_mod_tmp, y_mod_real_tmp);
						for (int i = 0; i < 2 * Nr; i++) { y_mod_real[i][j] = y_mod_real_tmp[i]; }
					}
					//MMSE detection done in paralell
					for (int i = 0; i < Nsym; i++)
					{
						for (int j = 0; j < 2 * Nr; j++) { y_mod_real_temp_omp[i][j] = y_mod_real[j][i]; }
					}
					Eigen::HouseholderQR<Eigen::MatrixXf> qr_Hreal(H_real);
					Q = qr_Hreal.householderQ();
					R = qr_Hreal.matrixQR().triangularView<Eigen::Upper>();
					// for EP detection
					// copy H_real
					for (int i = 0; i < Nr2; i++)
					{
						for (int j = 0; j < Nt2; j++)
						{
							Hreal_array[i * Nt2 + j] = H_real(i, j);
						}
					}
					// calculate HTH
					cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, Nt2, Nt2, Nr2, 1.0, Hreal_array, Nt2, Hreal_array, Nt2, 0.0, HTH_array, Nt2);
					//cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, Nt2, 1, Nr2, 1.0, ch_mtx_struct.H_real[cnt], Nt2, data_struct.rx_buffer_real[cnt], 1, 0.0, ch_mtx_struct.hty, 1);

					//omp_set_num_threads(8);
					//cout << 1 << endl;
					//omp_set_num_threads(4);
					/*auto start = std::chrono::system_clock::now();*/
					/*MatrixXf H_real_trans(2*Nt, 2*Nr);
					for (int i = 0; i < 2 * Nt; i++)
						for (int j = 0; j < 2 * Nr; j++)
							H_real_trans(i, j) = H_real(j, i);
					HTHIinv = (H_real_trans * H_real + Identity).inverse();*/
# pragma omp parallel for private(z)
					for (int j = 0; j < Nsym; j++)

						//for (int x=0; x<Nsym; x++)
					{
						// int j = x % Nsym;
						MatrixXf received_signal_omp(2 * Nr, 1);
						int k_kbest_extend = pow(2, ModType / 2) * K_KBEST;
						float** paths_kbest = new float* [k_kbest_extend];
						for (int i = 0; i < k_kbest_extend; i++) { paths_kbest[i] = new float[2 * Nt]; }
						float* ED = new float[K_KBEST];
						float* HTY_array_tmp = new float[2 * Nt];
						// Declared arrays for EP
						int flag = 1;
						int temp = 0;
						//float* extr_mean = new float[2 * Nr];
						//float* extr_var = new float[2 * Nr];
						double* vec_alpha_in = new double[2 * Nr];
						double* vec_gamma_in = new double[2 * Nr];
						double* p_symbol_rearrange = new double[Nt2 * pow(2, ModType / 2)];
						for (int i = 0; i < Nt2 * pow(2, ModType / 2); i++) { p_symbol_rearrange[i] = 1 / (pow(2, ModType / 2)); }
						for (int i = 0; i < Nr2; i++)
						{
							//extr_mean[i] = 0;
							//extr_var[i] = 0;
							vec_alpha_in[i] = 2;
							vec_gamma_in[i] = 0;
						}

						if (MIMODET == "MMSE")
						{
							MMSEdetection(H_real, Identity, received_signal, mmse_matrix, output_signal, HTH, y_mod_real_temp_omp[j], y_mmse_temp_omp[j], sigma_sq);
							//MMSEdetection_noInverse(H_real, Identity, received_signal, mmse_matrix, output_signal, HTH, HTHIinv, y_mod_real_temp_omp[j], y_mmse_temp_omp[j], sigma_sq);
						}
						if (MIMODET == "KBEST")
						{
							for (int i = 0; i < 2 * Nr; i++) { received_signal_omp(i, 0) = y_mod_real_temp_omp[j][i]; }
							z = Q.transpose() * received_signal_omp;
							KBestDetection(R, z, H_real, y_mod_real_temp_omp[j], ED, y_mmse_temp_omp[j], paths_kbest, K_KBEST, ModType, Nt, Nr);
							llr_kbest_bit(Nt, y_mmse_temp_omp[j], oriLLR_kbest[j], sigma_sq, ModType, K_KBEST, paths_kbest, ED);
						}
						if (MIMODET == "EP")
						{
							// test: save H_real, recv_signal
							/*save_mat_txt(y_mod_real_temp_omp[j], Nr2, "RxSig0128.txt");
							save_mat_txt(H_real, "Hreal0128.txt");*/
							// test: save H_real, recv_signal
							cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, Nt2, 1, Nr2, 1.0, Hreal_array, Nt2, y_mod_real_temp_omp[j], 1, 0.0, HTY_array_tmp, 1);
							//test:HTH_array
							/*for (int i = 0; i < Nt2; i++)
								for (int s = 0; s < Nt2; s++)
								{
									if (s == 0) std::cout << endl;
									std::cout << setw(10) << HTH_array[i * Nt2 + s];
								}
							cout << endl;*/
							// TEST:HTY
							/*for (int i = 0; i < Nr2; i++)
								received_signal_omp(i, 0) = y_mod_real_temp_omp[j][i];
							MatrixXf HTransY = H_real.transpose() * received_signal_omp;
							cout << "Eigen" << endl;
							for (int i = 0; i < Nt2; i++) cout << setw(10)<<HTransY(i, 0);
							cout << "mkl" << endl;
							for (int i = 0; i < Nt2; i++) cout << HTY_array[i];*/
							// TEST:HTY
							EPD_detection_float(flag, Nt, Nr, iter_MIMO, ModType, 1.0, sigma_sq, HTH_array, HTY_array_tmp, temp, ep_extr_mean_omp[j], ep_extr_var_omp[j], vec_alpha_in, vec_gamma_in, p_symbol_rearrange);
						}
						//test EP detection
						/*cout << "TEST:EP DETECTION" << endl;
						for (int i = 0; i < Nr2; i++) std::cout << extr_mean[i]<< endl;
						cout << "TEST:EP DETECTION" << endl;*/
						//test EP detection
						for (int i_path = 0; i_path < k_kbest_extend; i_path++) { delete[] paths_kbest[i_path]; }
						delete[] paths_kbest;
						delete[] ED;
						//delete[] extr_mean;
						//delete[] extr_var;
						delete[] vec_alpha_in;
						delete[] vec_gamma_in;
						delete[] p_symbol_rearrange;
						delete[] HTY_array_tmp;
					}
					// test: detection time
					/*auto detection_over = std::chrono::system_clock::now();
					auto duration_detection = std::chrono::duration_cast<std::chrono::microseconds>(detection_over - start);
					cout << endl;
					cout <<"Detection duration" << duration_detection.count() << "us" << endl;*/
					// test: detection time
					for (int i = 0; i < Nsym; i++)
					{
						for (int j = 0; j < 2 * Nt; j++) { y_mmse[j][i] = y_mmse_temp_omp[i][j]; }
					}
					//To bit domain output
					for (int j = 0; j < Nsym; j++)
					{
						for (int i = 0; i < 2 * Nt; i++) { y_mmse_tmp[i] = y_mmse[i][j]; }
						// With Spatial Modulation, Nt = Nv
						if (MIMODET == "MMSE") { llr_mmse_bit(Nt, y_mmse_tmp, oriLLR_tmp, sigma_sq, ModType); }
						if (MIMODET == "EP") { llr_ep_bit(Nt, ep_extr_mean_omp[j], ep_extr_var_omp[j], oriLLR_tmp, ModType); }
						if (MIMODET == "KBEST") { for (int i = 0; i < len_ori_llr_tmp; i++) { oriLLR_tmp[i] = oriLLR_kbest[j][i]; } }
						int i_bit = 0;
						int bit_col = 0;
						for (int j_inner = 0; j_inner < Nv; j_inner++)
						{
							for (int i_bit_col = 0; i_bit_col < ModType; i_bit_col++)
							{
								i_bit = ModType * j_inner + i_bit_col;
								bit_col = j * ModType + i_bit_col;
								oriLLRh[j_inner][bit_col] = oriLLR_tmp[i_bit];
								upLLR[j_inner][bit_col] = oriLLR_tmp[i_bit];
							}

						}
					}
					// test: LLR transfer time
					/*auto llr_to_bit_over = std::chrono::system_clock::now();
					auto duration_llr = std::chrono::duration_cast<std::chrono::microseconds>(llr_to_bit_over - detection_over);
					cout << endl;
					cout<< "LLR transfer duration" << duration_llr.count() << "us" << endl;*/
					// test: LLR transfer time
					if (flag_interleave)
					{
						if ((!flag_interleave_block) && flag_interleave_spatial)
						{
							for (int j = 0; j < Nh; j++)
							{
								for (int i = 0; i < Nv; i++) { oriLLR_tmp[i] = oriLLRh[i][j]; }
								randDeintrlv(intrlv_pattern_2D[j], oriLLR_tmp, Nv);
								for (int i = 0; i < Nv; i++)
								{
									oriLLRh[i][j] = oriLLR_tmp[i];
									upLLR[i][j] = oriLLR_tmp[i];
								}

							}
						}

						if (flag_interleave_block)
						{
							randDeIntrlvBlock(oriLLRh, Nh, Nv, ModType, interleave_rows);
							randDeIntrlvBlock(upLLR, Nh, Nv, ModType, interleave_rows);
						}

						if (flag_interleave_all) // Deinterleave in time domain
						{
							for (int j = 0; j < Nv; j++)
							{
								randDeintrlv(intrlv_pattern_2D_time[j], oriLLRh[j], Nh);
								randDeintrlv(intrlv_pattern_2D_time[j], upLLR[j], Nh);
							}
						}

						if (flag_interleave_partial)
						{
							for (int j = 0; j < Nv; j++)
							{
								for (int k = 0; k < intrlv_pattern_num; k++)
								{
									randDeintrlv(intrlv_pattern_2D_time_partial[j][k], oriLLRh[j] + k * ModType, ModType);
									randDeintrlv(intrlv_pattern_2D_time_partial[j][k], upLLR[j] + k * ModType, ModType);
								}
							}
						}

					}

					for (int i = 0; i < Nv; i++)
					{
						for (int j = 0; j < Nh; j++)
						{
							total_bit_cnt++;
							int bit_dec = (oriLLRh[i][j] < 0);
							if (bit_dec != x[i][j]) n_det_err++;
						}
					}

				}
				// test: detection errors

				// test: detection errors
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
			for (int iter = 1; iter <= softiter; iter++)
			{
				//vertical decoding
		   // test: decoding time
			//auto decoding_start = std::chrono::system_clock::now();

			// test: decoding time
#pragma omp parallel for
				for (int decodingi = 0; decodingi < Nh; decodingi++) // 291.70kb
				{

					//for SCL Decoding
					int Counth = 0, Countinfoh = 0;
					float Qsumunh = 0;
					int** sum = new int* [L];
					//for (int i = 0; i < L; i++)  sum[i] = new  int[Nmax];
					int** u1 = new int* [L];
					for (int i = 0; i < L; i++)  u1[i] = new _MM_ALIGN16 int[Nmax];
					for (int i = 0; i < L; i++)  sum[i] = new _MM_ALIGN16 int[Nmax];
					//for (int i = 0; i < L; i++)  u1[i] = new  int[Nmax];

					float* LLR_in = new float[L];
					float* W = new float[2 * L];
					int* SCL_index = new int[L];
					int* better = new int[L];
					int* worse = new int[L];
					int* indexpm = new int[L];
					float* PM = new float[L];
					float** LLR = new float* [L];
					for (int i = 0; i < L; i++)  LLR[i] = new _MM_ALIGN16 float[2 * Nmax + 2];
					// for (int i = 0; i < L; i++)  LLR[i] = new float[2 * Nmax + 2];
					////////////////////////////////
					for (int j = 0; j < Nv; j++)
						for (int k = 0; k < L; k++)
							LLR[k][j] = upLLR[j][decodingi];
					//cout << LLR[0][j] << " ";

				//cout << endl;
					int index = 0;
					/*			if (L == 1)
								{
									//T_start = clock();
									Count = 0;
									decode(*LLR, *sum, (int)(log(Nv) / log(2)), A_Acv, flag, Nv, Kv);
									//PolarEncode_xor(*u1, *sum, N);
									//T_finish = clock();
								}
								else
								{*/
					for (int k = 0; k < L; k++) { PM[k] = 0; indexpm[k] = k; }
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
					if (adaptive_beta)
					{
						if (softmethod == "Pyndiah") Pyndiahv_adaptivebeta(LLR, sum, oriLLRh, indexpm, iter * 2 - 2, decodingi, upLLR);
						if (softmethod == "MITSO") MITSOv(LLR, PM, sum, oriLLRh, indexpm, iter * 2 - 2, decodingi, upLLR, Qsumunh);
					}
					// iter * 2 - 2
					else
					{
						//cout << iter << endl;
						if (softmethod == "Pyndiah") Pyndiahv(LLR, sum, oriLLRh, indexpm, iter, decodingi, upLLR);
						if (softmethod == "MITSO") MITSOv(LLR, PM, sum, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunh);
					}
					// test: oriLLR and upLLR
					/*cout << endl;
					for (int i = 0; i < Nv; i++)
						for (int j = 0; j < Nh; j++)
							cout << oriLLRh[i][j] - upLLR[i][j] << "   ";
					cout << endl;*/
					/*cout << oriLLRh[0][0] << endl;
					cout << upLLR[0][0] << endl;*/
					// test: oriLLR and upLLR
					// half iteration Enabled
					if (half_iter && iter == softiter)
					{
						if (iter == softiter) PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), Nh);
						for (int j = 0; j < Kv; j++)
							umid[Av[j]][decodingi] = u1[index][Av[j]];
					}
					// half iteration
					////////delete SCL variables
					for (int i = 0; i < L; i++)
						delete[]LLR[i];
					delete[] LLR;
					for (int i = 0; i < L; i++)
						delete[]sum[i];
					for (int i = 0; i < L; i++)
						delete[]u1[i];
					delete[] u1;
					delete[] sum;
					delete[] PM;
					delete[] LLR_in;
					delete[] W;
					delete[] SCL_index;
					delete[] better;
					delete[] worse;
					delete[] indexpm;

				}
				//cout << "Horizon Start:\n";
				//horizontal decoding
				// #
				if (weirditer) { continue; }

				// test: LLR after horizontal decoding
				/*cout << "After horizontal decoding" << endl;
				for (int i = 0; i < Nv; i++)
				{
					for (int j = 0; j < Nh; j++)
					{
						cout << (oriLLRh[i][j] < 0) << " ";
					}
					cout << endl;
				}*/
				// test: LLR after horizontal decoding
				// inner deinterleaving
				if (flag_interleave_inner)
				{
					float* tmp_llr_to_interleave = new float[Nv];
					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < Nv; j++) { tmp_llr_to_interleave[j] = oriLLRh[j][i]; }
						randDeIntrlvPartial(intrlv_pattern_2D_inner[i], tmp_llr_to_interleave, Av, Nv, Kv);
						for (int j = 0; j < Nv; j++) { oriLLRh[j][i] = tmp_llr_to_interleave[j]; }
					}
					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < Nv; j++) { tmp_llr_to_interleave[j] = upLLR[j][i]; }
						randDeIntrlvPartial(intrlv_pattern_2D_inner[i], tmp_llr_to_interleave, Av, Nv, Kv);
						for (int j = 0; j < Nv; j++) { upLLR[j][i] = tmp_llr_to_interleave[j]; }
					}
					delete[] tmp_llr_to_interleave;
				}


				/*auto decoding_over = std::chrono::system_clock::now();
				auto duration_decoding = std::chrono::duration_cast<std::chrono::microseconds>(decoding_over - decoding_start);
				cout << endl;
				cout << "Decoding duration" << duration_decoding.count() << "us" << endl;*/

#pragma omp parallel for
				for (int decodingi = 0; decodingi < Nv; decodingi++)
				{
					// 293.70kb
					// cout << "KKKK" << endl;
					if (half_iter && (iter == softiter)) { break; }
					string polarde_now = polarde_array[iter - 1];
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

					/*			if (L == 1)
								{
									//T_start = clock();
									Count = 0;
									decode(*LLR, *sum, (int)(log(Nh) / log(2)), A_Ach, flag, Nh, Kh);
									if (iter == softiter) PolarEncode_xor(*u1, *sum, Nh);
									//T_finish = clock();
								}
								else
								{*/
								//for (int k = 0; k < L; k++)  PM[k] = 0, indexpm[k] = k; 
								////T_start = clock();
								//decode_list(LLR, sum, (int)(log(Nh) / log(2)), A_Ach, 0, Nh, Kh, PM, L, (int)(log(L) / log(2)), LLR_in, W, SCL_index, better, worse, &Countv, &Countinfov, &Qsumunv);
								////T_finish = clock();
								////Decoding_Duration += (T_finish - T_start);
								//for (int i = 0; i < L - 1; i++)
								//	for (int j = i + 1; j < L; j++)
								//		if (PM[indexpm[i]] < PM[indexpm[j]])
								//		{
								//			int ttt = indexpm[i];
								//			indexpm[i] = indexpm[j];
								//			indexpm[j] = ttt;
								//		}


							/*
							* Encoding: Info -> Horizontal Encoding -> A -> set to 0 -> B -> Horizontal Encoding -> Vertical Encoding
							* Decoding: Received Info -> Horizontal Decoding -> Vertical Decoding -> Horizontal Decoding
							*/
					if (polarde_now == "SLD")
					{
						LSD_decode_outall(u_decod, u1_LSD, LLR[0], GGh, A_Ach, Nh, L_LSD, n_paths);
						// if (iter == softiter) PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), Nh);
						// for soft information exchange, inverse coding:
						// THIS LINE is WRONG! BUT IT SEEMS BETTER THAN THE RIGHT VERSION
						for (int i = 0; i < L_LSD; i++) { PolarEncode_xor(sum_LSD[i], u1_LSD[i], Nh); }
						// for (int i = 0; i < L_LSD; i++) { PolarEncode_xor(sum_LSD[i], u1_LSD[i], Nh); }
						// if (iter == softiter) PolarEncode_xor(u1[0], sum_LSD[0], Nh);
						index = indexpm[0];
						// test: after horizontal decoding
						/*if (decodingi == Av[0])
						{
							cout << "Direct output of horizontal decoding:" << endl;
							for (int i = 0; i < Nh; i++) { cout << u1_LSD[0][i] << " "; }
							cout << endl;
							cout << "horizontal decoder out, encode once:" << endl;
							for (int i = 0; i < Nh; i++) { cout << sum_LSD[0][i] << " "; }
							cout << endl;
						}*/
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
						// iter*2-1
						if (iter < softiter)
						{
							if (adaptive_beta)
							{
								if (softmethod == "Pyndiah") Pyndiahh_LSD_adaptivebeta(LLR, sum_LSD, oriLLRh, indexpm, iter, decodingi, upLLR);
								if (softmethod == "MITSO") MITSOh(LLR, PM, sum_LSD, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunv);
							}
							else
							{
								// cout << "IN LSD!!!" << endl;
								if (softmethod == "Pyndiah")
									Pyndiahh_LSD_getsum(LLR, sum_LSD, oriLLRh, indexpm, iter, decodingi, upLLR, sum_sdis_array[decodingi], GGh, Ach);
								//if (softmethod == "Pyndiah") Pyndiahh_LSD(LLR, sum_LSD, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR);
								if (softmethod == "MITSO") MITSOh(LLR, PM, sum_LSD, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunv);
							}
							//	getchar();
						}
						else
							for (int j = 0; j < Kh; j++)
								umid[decodingi][Ah[j]] = u1_LSD[index][Ah[j]];
					}
					if (polarde_now == "SCL")
					{
						for (int k = 0; k < L; k++)  PM[k] = 0, indexpm[k] = k;
						//T_start = clock();
						decode_list(LLR, sum, (int)(log(Nh) / log(2)), A_Ach, 0, Nh, Kh, PM, L, (int)(log(L) / log(2)), LLR_in, W, SCL_index, better, worse, &Countv, &Countinfov, &Qsumunv);
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
						if (iter == softiter) PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), Nh);
						index = indexpm[0];
						//}
						//calc distance s()
						if (iter < softiter)
						{
							if (adaptive_beta)
							{
								if (softmethod == "Pyndiah") Pyndiahh_adaptivebeta(LLR, sum, oriLLRh, indexpm, iter, decodingi, upLLR);
								if (softmethod == "MITSO") MITSOh(LLR, PM, sum, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunv);
								//	getchar();
							}
							else
							{
								//if (softmethod == "Pyndiah") Pyndiahh(LLR, sum, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR);
								if (softmethod == "Pyndiah") Pyndiahh_getsum(LLR, sum, oriLLRh, indexpm, iter, decodingi, upLLR, sum_sdis_array[decodingi]);
								if (softmethod == "MITSO") MITSOh(LLR, PM, sum, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunv);
								//	getchar();
							}
						}
						else
							for (int j = 0; j < Kh; j++)
								umid[decodingi][Ah[j]] = u1[index][Ah[j]];
					}

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


				if (flag_interleave_inner)
				{
					float* tmp_llr_to_interleave = new float[Nv];
					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < Nv; j++) { tmp_llr_to_interleave[j] = oriLLRh[j][i]; }
						randIntrlvPartial(intrlv_pattern_2D_inner[i], tmp_llr_to_interleave, Av, Nv, Kv);
						for (int j = 0; j < Nv; j++) { oriLLRh[j][i] = tmp_llr_to_interleave[j]; }
					}
					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < Nv; j++) { tmp_llr_to_interleave[j] = upLLR[j][i]; }
						randIntrlvPartial(intrlv_pattern_2D_inner[i], tmp_llr_to_interleave, Av, Nv, Kv);
						for (int j = 0; j < Nv; j++) { upLLR[j][i] = tmp_llr_to_interleave[j]; }
					}
					delete[] tmp_llr_to_interleave;
				}
			}
			int calcsum[Nv], calcu1[Nv];
			int calcu1h[Nh];
			// calcsum is from the direct output of the decoder
			// if output is the result of row decoding (meaning there is another column decoding to be done)
			if (half_iter)
			{
				for (int i = 0; i < Kv; i++)
				{
					PolarEncode_xor(calcu1h, umid[Av[i]], Nh);
					for (int j = 0; j < Kh; j++)
						uout[Av[i]][Ah[j]] = calcu1h[Ah[j]];
				}
			}
			// if output is the result of column decoding (meaning that there is another row decoding to be done)
			else
			{
				for (int i = 0; i < Kh; i++)
				{
					for (int j = 0; j < Nv; j++)
						calcsum[j] = umid[j][Ah[i]];
					PolarEncode_xor(calcu1, calcsum, Nv);
					for (int j = 0; j < Kv; j++)
						uout[Av[j]][Ah[i]] = calcu1[Av[j]];
				}
			}

			int calcsumh[Nh];
			for (int i = 0; i < Kv; i++)
			{
				for (int j = 0; j < Nh; j++) { calcsumh[j] = uout[Av[i]][j]; }
				PolarEncode_xor(uout[Av[i]], calcsumh, Nh);
			}
			// test: final output
			/*cout << "final output" << endl;
			for (int k = 0; k < Nh; k++) { cout << uout[Av[0]][k] << " "; }
			cout << endl;*/
			// test: final output
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
				/*printf("\r<%d> SNR= %.2f FER= %d / %d = %f BER= %d / %d = %lf E-Det = %f, E-MRB = %f"  ,
					frame, nSNR, Num_Frame_Error, frame, (double)Num_Frame_Error / frame, Num_Error, (N_Info * frame), (double)Num_Error / N_Info / frame,
					(double)n_det_err/total_bit_cnt, (double)n_det_err_mrb/(total_bit_cnt/(ModType/2)));*/
				int L_sdis_sum = L_LSD;
				sdis_sum = 0;
				for (int i = 0; i < Nv; i++) sdis_sum += sum_sdis_array[i];
				if (polarde_array[0] == "SCL") L_sdis_sum = L;
				// cout << "sdis_sum_avg" << sdis_sum / frame / L_sdis_sum << endl;
				printf("\r<%d> SNR= %.2f FER= %d / %d = %f BER= %d / %d = %lf ,sdis_avg = %f, E-MRB = %f",
					frame, nSNR, Num_Frame_Error, frame, (double)Num_Frame_Error / frame, Num_Error, (N_Info * frame), (double)Num_Error / N_Info / frame, (double)sdis_sum / frame / L_sdis_sum, (double)n_det_err / total_bit_cnt);
				fflush(stdout);
			}


			// output interleaving pattern
			//if (Num_Frame_Error % record_frames == 0 && Num_Frame_Error > 0 && Num_Frame_Error!=n_frame_error_pre)
			//{
			//	n_frame_error_pre = Num_Frame_Error;
			//	cout << endl;
			//	int n_file = Num_Frame_Error / record_frames;
			//	frame_delta = frame - frame_tmp;
			//	cout << "Delta Frames: " << frame_delta << endl;
			//	frame_tmp = frame;
			//	string ofile_intrlv_name = string("intrlv_pattern_") + to_string(n_file) + string("_frame_") + to_string(frame_delta) + string(".txt");
			//	string ofile_spatial_intrlv_name = string("intrlv_pattern_spatial") + to_string(n_file) + string(".txt");
			//	ofstream ofile_intrlv(ofile_intrlv_name);
			//	ofstream ofile_spatial_intrlv(ofile_spatial_intrlv_name);
			//	for (int i = 0; i < Nv; i++)
			//	{
			//		for (int j = 0; j < intrlv_pattern_num; j++)
			//		{

			//			for (int k = 0; k < ModType; k++)
			//				ofile_intrlv << intrlv_pattern_2D_time_partial[i][j][k] << " ";
			//			//cout << endl;
			//		}
			//	}
			//	ofile_intrlv.close();

			//	for (int i = 0; i < Nh; i++)
			//	{
			//		for (int j = 0; j < Nv; j++)
			//			ofile_spatial_intrlv << intrlv_pattern_2D[i][j] << " ";
			//	}
			//	ofile_spatial_intrlv.close();

			//	// re-generate interleaving pattern
			//	generateIntrlvPattern(intrlv_pattern, Nv);
			//	for (int i = 0; i < Nh; i++) { generateIntrlvPattern(intrlv_pattern_2D[i], Nv); }
			//	for (int i = 0; i < Nv; i++) { generateIntrlvPattern(intrlv_pattern_2D_time[i], Nh); }
			//	if (flag_interleave_partial)
			//	{
			//		for (int i = 0; i < Nv; i++)
			//		{
			//			for (int j = 0; j < intrlv_pattern_num; j++)
			//			{
			//				generateIntrlvPattern(intrlv_pattern_2D_time_partial[i][j], ModType);
			//				/*for (int k = 0; k < ModType; k++)
			//					cout << intrlv_pattern_2D_time_partial[i][j][k];
			//				cout << endl;*/
			//			}
			//		}
			//	}
			//	//Num_Frame_Error = 0;
			//}


			//printf("%d %d", Num_Frame_Error, Num_Error);
			//getchar();
			//Decoding_Duration += (T_finish - T_start);
			/*auto dec_end = std::chrono::system_clock::now();
			auto duration_decoding = std::chrono::duration_cast<std::chrono::microseconds>(dec_end - det_end);
			auto duration_total = std::chrono::duration_cast<std::chrono::microseconds>(dec_end - start);
			file_dec << duration_decoding.count() << endl;
			file_total << duration_total.count() << endl;*/
			if (frame % 16 == 0)
			{
				file_err_num << Num_Error << endl;
			}

		}
		//Decoding_Duration = Decoding_Duration * 100 / CLOCKS_PER_SEC / frame;
		//Decoding_Duration = Decoding_Duration / CLOCKS_PER_SEC / Max_frame;

		//printf("%f  %f  %f  %f\n", nSNR, Decoding_Duration, (double)Num_Error / N_Info / Max_frame, (double)Num_Frame_Error / Max_frame);
		printf("%f %f %f\n", nSNR, (double)Num_Frame_Error / frame, (double)Num_Error / N_Info / frame);
		BER_array[SNR_count] = (double)Num_Error / N_Info / frame;
		FER_array[SNR_count] = (double)Num_Frame_Error / frame;
		// save to .txt file
		string polarde_serial = "";
		for (int i = 0; i < softiter; i++)
			polarde_serial = polarde_serial + polarde_array[i] + string("_");
		string file_name = string("ResFinal\\2D_LSD_D_") + polarde_serial + MIMODET + string("T") + to_string(Nt) + string("R") + to_string(Nr)
			+ string("Q") + to_string(int(pow(2, ModType))) + string("L") + to_string(Nh) + string("K") + to_string(Kh)
			+ string("T") + to_string(softiter) + string("H") + to_string(half_iter) + string("UB") + to_string(usebeta) + string("I") + to_string(flag_interleave)
			+ string(".txt");
		save_array_txt(nSNR, BER_array[SNR_count], FER_array[SNR_count], file_name);
		int L_sdis_sum = L_LSD;
		sdis_sum = 0;
		for (int i = 0; i < Nv; i++) sdis_sum += sum_sdis_array[i];
		if (polarde_array[0] == "SCL") L_sdis_sum = L;
		cout << "sdis_sum_avg" << sdis_sum / frame / L_sdis_sum << endl;
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
	for (int i = 0; i < Nh; i++)
		delete[] oriLLR_kbest[i];
	delete[] oriLLR_kbest;
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
	for (int i = 0; i < mimo_blk_height; i++) {
		delete[] x_mod[i];
	}
	delete[] x_mod;
	for (int i = 0; i < Nr; i++) {
		delete[] y_mod[i];
	}
	delete[] y_mod;
	for (int i = 0; i < 2 * mimo_blk_height; i++) {
		delete[] y_mmse[i];
	}
	delete[] y_mmse;
	for (int i = 0; i < 2 * Nr; i++) {
		delete[] y_mod_real[i];
	}
	for (int i = 0; i < mimo_blk_length; i++)
	{
		delete[] y_mod_real_temp_omp[i];
	}
	delete[] y_mod_real_temp_omp;
	for (int i = 0; i < mimo_blk_length; i++)
	{
		delete[] y_mmse_temp_omp[i];
	}
	delete[] y_mmse_temp_omp;
	for (int i = 0; i < mimo_blk_length; i++)
	{
		delete[] ep_extr_mean_omp[i];
	}
	delete[] ep_extr_mean_omp;
	for (int i = 0; i < mimo_blk_length; i++)
	{
		delete[] ep_extr_var_omp[i];
	}
	delete[] ep_extr_var_omp;
	delete[] y_mod_real;
	delete[] x_mod_tmp;
	delete[] y_mod_tmp;
	delete[] y_mod_real_tmp;
	delete[] x_tmp;
	delete[] oriLLR_tmp;
	delete[] sum_sdis_array;
	delete[] Hreal_array;
	delete[] HTH_array;
	delete[] HTY_array;
	delete[] intrlv_pattern_2D_inner;
}


void system2D_LSD_epdet_irregular(double* BER_array, double* FER_array)
{
	std::ofstream file_total("total_timing.txt", std::ios::app);
	std::ofstream file_det("det_timing.txt", std::ios::app);
	std::ofstream file_dec("dec_timing.txt", std::ios::app);
	std::ofstream file_err_num("total_frames.txt", std::ios::app);
	auto det_end = std::chrono::system_clock::now();
	int Nt2 = 2 * Nt;
	int Nr2 = 2 * Nr;
	int cnt_short_codewords = 0;
	for (int i = 0; i < Nv; i++) if (channel_index_predefined[i] == 1) { cnt_short_codewords++; }
	int* bad_channel_index = new int[cnt_short_codewords];
	int* good_channel_index = new int[Nv - cnt_short_codewords];
	int cnt_good = 0;
	int cnt_bad = 0;
	for (int i = 0; i < Nv; i++) { if (channel_index_predefined[i] == 0) { good_channel_index[cnt_good] = i; cnt_good++; } else { bad_channel_index[cnt_bad] = i; cnt_bad++; } }
	/*if (!file.is_open()) {
		std::cerr << "Error opening file: " << file_name << std::endl;
		return;
	}*/

	//file << setw(20) << SNR << setw(20) << BER << setw(20) << FER << "\n";

	//freopen("result.txt", "w", stdout);
	printsetting();
	cout << "\nSNR	 " << "FER      BER" << endl;
	int N_Info = Kv * Kh - CRCL, Nmax;
	if (irregular_coding) { N_Info = border * K_time_array[0] + (Kv - border) * K_time_array[1]  - CRCL; }
	//Interleaving pattern
	/*std::vector<int> intrlv_pattern(Nv);
	std::vector<int>* intrlv_pattern_2D = new vector<int>[Nh];
	for (int i = 0; i < Nh; i++) intrlv_pattern_2D[i].resize(Nv);*/
	int intrlv_pattern_num = Nh / ModType;
	std::vector<int> intrlv_pattern(Nv);
	std::vector<int>* intrlv_pattern_2D = new vector<int>[Nh];
	std::vector<int>* intrlv_pattern_2D_time = new vector<int>[Nv];
	std::vector<int>* intrlv_pattern_2D_inner = new vector<int>[Nh];
	std::vector<int>** intrlv_pattern_2D_time_partial = new vector<int> *[Nv];
	for (int i = 0; i < Nv; i++) intrlv_pattern_2D_time_partial[i] = new vector<int>[intrlv_pattern_num];
	for (int i = 0; i < Nv; i++)
	{
		for (int j = 0; j < intrlv_pattern_num; j++)
		{
			intrlv_pattern_2D_time_partial[i][j].resize(ModType);
		}
	}
	for (int i = 0; i < Nh; i++) intrlv_pattern_2D[i].resize(Nv);
	for (int i = 0; i < Nh; i++) intrlv_pattern_2D_inner[i].resize(Kv);
	for (int i = 0; i < Nv; i++) intrlv_pattern_2D_time[i].resize(Nh);

	float CodeRate = (float)(N_Info) / (Nh * Nv);
	Nmax = max(Nh, Nv);
	// num of symbols per vertical codeword
	int Nsym = 0;
	if (MOD_DIM == "Spatial") { Nsym = Nv / ModType; }
	else if (MOD_DIM == "Temporal") { Nsym = Nh / ModType; }
	int symbol_length = ModType;

	// Size of a mimo block is given by height(spatial) * length(temporal)
	int mimo_blk_height = 0;
	int mimo_blk_length = 0;
	if (MOD_DIM == "Spatial")
	{
		mimo_blk_height = Nsym;
		mimo_blk_length = Nh;
	}
	else if (MOD_DIM == "Temporal")
	{
		mimo_blk_height = Nv;
		mimo_blk_length = Nsym;
	}
	// Initialization
	unsigned char* a_data = new unsigned char[N_Info];// Information CodeWord
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
	int** x_parity_check = new int* [Nv];
	for (int i = 0; i < Nh; i++) x_parity_check[i] = new int[Nh];

	// Arrays for MIMO
	int** x_mmse_out = new int* [Nv]; // delete!
	for (int i = 0; i < Nv; i++) { x_mmse_out[i] = new int[Nh]; }

	//std::complex<float>** x_mod = new std::complex<float>* [Nsym]; // modulated symbols (spatial modulation) DELETE!
	//for (int i = 0; i < Nsym; i++) { x_mod[i] = new std::complex<float>[Nh]; }
	//std::complex<float>** x_mod_temporal = new std::complex<float>*[Nv]; // modulated symbols (temporal modulation) DELETE!
	//for (int i = 0; i < Nv; i++) x_mod_temporal[i] = new std::complex<float>[Nsym];
	//std::complex<float>** y_mod = new std::complex<float>* [Nr]; // Output symbols after AWGN channels // delete!
	//for (int i = 0; i < Nr; i++) { y_mod[i] = new std::complex<float>[Nh]; }
	//float** y_mod_real = new float* [2 * Nr];
	//for (int i = 0; i < 2*Nr; i++) y_mod_real[i] = new float[Nh];
	//float** y = new float* [Nv];// Output Signal After AWGN Channel
	//for (int i = 0; i < Nv; i++) { y[i] = new float[Nh]; }
	//float** y_mmse = new float* [2*Nsym];
	//for (int i = 0; i < 2*Nsym; i++) { y_mmse[i] = new float[Nh]; }  // MMSE detector output // delete!
	//float* y_mmse_tmp = new float[2*Nsym];
	//std::complex<float>* x_mod_tmp = new std::complex<float>[Nsym]; // modulated symbols of a single column // delete!
	//std::complex<float>* y_mod_tmp = new std::complex<float>[Nr]; // received symbols of a single column // delete!
	//float* y_mod_real_tmp = new float[2*Nr];
	//// float* y_mod_tmp = new float[Nr];  // delete!


	std::complex<float>** x_mod = new std::complex<float>*[mimo_blk_height]; // modulated symbols (spatial modulation) DELETE!
	for (int i = 0; i < mimo_blk_height; i++) { x_mod[i] = new std::complex<float>[mimo_blk_length]; }
	std::complex<float>** y_mod = new std::complex<float>*[Nr]; // Output symbols after AWGN channels // delete!
	for (int i = 0; i < Nr; i++) { y_mod[i] = new std::complex<float>[mimo_blk_length]; }
	float** y_mod_real = new float* [2 * Nr];
	for (int i = 0; i < 2 * Nr; i++) y_mod_real[i] = new float[mimo_blk_length];
	float** y = new float* [Nv];// Output Signal After AWGN Channel
	for (int i = 0; i < Nv; i++) { y[i] = new float[Nh]; }
	float** y_mmse = new float* [2 * mimo_blk_height];
	for (int i = 0; i < 2 * mimo_blk_height; i++) { y_mmse[i] = new float[mimo_blk_length]; }  // MMSE detector output // delete!
	float* y_mmse_tmp = new float[2 * mimo_blk_height];
	std::complex<float>* x_mod_tmp = new std::complex<float>[mimo_blk_height]; // modulated symbols of a single column // delete!
	std::complex<float>* y_mod_tmp = new std::complex<float>[Nr]; // received symbols of a single column // delete!
	float* y_mod_real_tmp = new float[2 * Nr];
	// float* y_mod_tmp = new float[Nr];  // delete!
	float** y_mod_real_temp_omp = new float* [mimo_blk_length]; // new delete!
	for (int i = 0; i < mimo_blk_length; i++) { y_mod_real_temp_omp[i] = new float[2 * Nr]; }
	float** y_mmse_temp_omp = new float* [mimo_blk_length]; // new delete!
	for (int i = 0; i < mimo_blk_length; i++) { y_mmse_temp_omp[i] = new float[2 * mimo_blk_height]; }
	// ep
	float** ep_extr_mean_omp = new float* [mimo_blk_length];
	for (int i = 0; i < mimo_blk_length; i++) { ep_extr_mean_omp[i] = new float[2 * mimo_blk_height]; }
	float** ep_extr_var_omp = new float* [mimo_blk_length];
	for (int i = 0; i < mimo_blk_length; i++) { ep_extr_var_omp[i] = new float[2 * mimo_blk_height]; }
	int* x_tmp = new int[Nv]; // delete!

	// paths selected by k-best detection


	// CRC bits
	unsigned char* u_crc = new unsigned char[Kv * Kh];

	float** upLLR = new float* [Nv];
	for (int i = 0; i < Nv; i++) { upLLR[i] = new _MM_ALIGN16 float[Nh]; }
	//float** oriLLRv = new float* [Nh];
	//for (int i = 0; i < Nh; i++) { oriLLRv[i] = new _MM_ALIGN16 float[2 * Nv + 2]; }
	float** oriLLRh = new float* [Nv];
	for (int i = 0; i < Nv; i++) { oriLLRh[i] = new _MM_ALIGN16 float[Nh]; }
	int len_ori_llr_tmp = 0;
	if (MOD_DIM == "Spatial") { len_ori_llr_tmp = Nv; }
	else if (MOD_DIM == "Temporal") { len_ori_llr_tmp = ModType * Nt; }
	float* oriLLR_tmp = new _MM_ALIGN16 float[len_ori_llr_tmp];  // delete!
	//float** oriLLR_kbest = new float* [Nh];
	float** oriLLR_kbest = new float* [mimo_blk_length];
	for (int i = 0; i < Nh; i++) { oriLLR_kbest[i] = new float[len_ori_llr_tmp]; }
	double* sum_sdis_array = new double[Nv];
	for (int i = 0; i < Nv; i++) sum_sdis_array[i] = 0;
	int* n_tried_patterns_h = new int[softiter];
	for (int i = 0; i < softiter; i++) n_tried_patterns_h[i] = 0;
	int* n_tried_patterns_v = new int[softiter];
	for (int i = 0; i < softiter; i++) n_tried_patterns_v[i] = 0;

	//*****************************************************************
	// EP detector
	float* Hreal_array = new float[(2 * Nr) * (2 * Nt)]; // delete
	float* HTH_array = new float[(2 * Nt) * (2 * Nt)]; // delete
	float* HTY_array = new float[2 * Nt]; // delete 





	int** GGh = new int* [Nh];
	int** GGv = new int* [Nv];


	MatrixXf Gh = Fn(Nh);  // Generate Matrix of horizontal coding
	MatrixXf Gv = Fn(Nv);  // Generate Matrix of vertical coding
	MatrixXcf H(Nr, Nt); // Complex channel Matrix H: Nr x Nt
	//Matrix<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic> H = MatrixXcf::Random(Nr, Nt); // Real domain channel matrix
	//// //MatrixX Definitions
	MatrixXf H_real(2 * Nr, 2 * Nt); // Real domain channel matrix
	MatrixXf received_signal(2 * Nr, 1);
	MatrixXf HTH(2 * Nt, 2 * Nt);
	MatrixXf Identity = MatrixXf::Identity(2 * Nt, 2 * Nt);
	MatrixXf output_signal(2 * Nt, 1);
	MatrixXf mmse_matrix(2 * Nt, 2 * Nr);
	MatrixXf HTHIinv(2 * Nt, 2 * Nt);

	// for KBest
	MatrixXf R(2 * Nt, 2 * Nt);
	MatrixXf Q(2 * Nr, 2 * Nt);
	MatrixXf z(2 * Nt, 1);



	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> H_real = MatrixXf::Random(2 * Nr, 2 * Nt); // Real domain channel matrix
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> received_signal = MatrixXf::Random(2 * Nr,1);
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> HTH = MatrixXf::Random(2 * Nt, 2 * Nt);
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> Identity = MatrixXf::Random(2 * Nt, 2 * Nt);
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> output_signal = MatrixXf::Random(2 * Nt, 1);
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> mmse_matrix = MatrixXf::Random(2 * Nt, 2 * Nr);


	// Matrix Definitions
	//Eigen::Matrix<float, 2 * Nr, 2 * Nt> H_real; // ʵ���ŵ�����
	//Eigen::Matrix<float, 2 * Nr, 1> received_signal; // �����źž���
	//Eigen::Matrix<float, 2 * Nt, 2 * Nt> HTH; // HTH����
	//Eigen::Matrix<float, 2 * Nt, 2 * Nt> Identity = Eigen::Matrix<float, 2 * Nt, 2 * Nt>::Identity(); // ��λ����
	//Eigen::Matrix<float, 2 * Nt, 1> output_signal; // ����źž���
	//Eigen::Matrix<float, 2 * Nt, 2 * Nr> mmse_matrix; // MMSE����


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
	vector<int>* construction_time_all = new vector<int>[Nv];
	for (int i = 0; i < Nv; i++) { construction_time_all[i].resize(Nh); }
	int* Ah = new int[Kh]; // Horizontal Info Bits
	int* Ach = new int[Nh - Kh]; // Horizontal Frozen Bits
	int* A_Ach = new int[Nh];
	int* Av = new int[Kv];
	int* Acv = new int[Nv - Kv];
	int* A_Acv = new int[Nv];
	int** A_time_all = new int* [Nv];
	for (int i = 0; i < Nv; i++) { A_time_all[i] = new int[Nh]; for (int j = 0; j < Nh; j++) A_time_all[i][j] = 0; }
	int** Ac_time_all = new int* [Nv];
	for (int i = 0; i < Nv; i++) { Ac_time_all[i] = new int[Nh]; for (int j = 0; j < Nh; j++) Ac_time_all[i][j] = 0;}
	int** A_Ac_time_all = new int* [Nv];
	for (int i = 0; i < Nv; i++) { A_Ac_time_all[i] = new int[Nh]; for (int j = 0; j < Nh; j++) A_Ac_time_all[i][j] = 0; }
	int* K_time_all = new int[Nv];
	for (int i = 0; i < Nv; i++) K_time_all[i] = 0;
	int* tempv_ordered = new int[Nv];
	int** all_error_positions = new int* [Nv];
	for (int i = 0; i < Nv; i++) {
		all_error_positions[i] = new int[Nh];
		for (int j = 0; j < Nh; j++) all_error_positions[i][j] = 0;
	}
	int* use_short_codeword = new int[Nv];
	for (int i = 0; i < Nv; i++) { use_short_codeword[i] = 0; }
	/*int* use_short_codeword = new int[Nv];
	for (int i = 0; i < Nv; i++) { use_short_codeword[i] = 0; }*/

	float** channel_scale_factor = new float* [Nr];
	for (int i = 0; i < Nr; i++) { channel_scale_factor[i] = new float[Nt]; for (int j = 0; j < Nt; j++) { channel_scale_factor[i][j] = 0.0; } }
	//*****************************************************************
		// Main Function
	// srand((unsigned int)time(0));
	std::random_device rd;              // ��ȡ���������
	std::mt19937 gen(rd());             // ʹ�����ӳ�ʼ��mt19937������
	std::uniform_int_distribution<> distrib(1, 1000000);
	int Num_Error, Num_Frame_Error;
	float Decoding_Duration;

	int Num_Error_new, Num_Frame_Error_new;
	float Decoding_Duration_new;

	clock_t T_start, T_finish;
	int SNR_count = 0;
	int anssc, anssd;
	for (double nSNR = Min_SNR; nSNR <= Max_SNR; nSNR += SNR_step)
	{
		for (int i = 0; i < Nv; i++) sum_sdis_array[i] = 0;
		double sdis_sum = 0;
		// TEST: det err
		int n_det_err = 0; // n error bits at detector output
		int n_det_err_mrb = 0; // error bits at detector output(MRBs)
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
		// for (int i = 0; i < Nv; i++) { tempv[i] = vertical_construction_fixed[i]; }
		// cout << endl;
		for (int i = 0; i < Nv; i++) { tempv_ordered[i] = tempv[i]; /*cout << tempv_ordered[i] << ",";*/ }

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

		/*for (int i = 0; i < Nh; i++)
		{
			cout << A_Ach[i] << ",";
		}
		cout << endl;*/
		// for (int i = 0; i < Kv; i++) { cout << Av[i] << endl; }
		for (int k = 0; k < Nv; k++)
		{
			bool bfind = false;
			int channel_index = 0;
			for (int i_find = 0; i_find < Kv; i_find++) { if (tempv_ordered[i_find] == k) { bfind = true;  channel_index = i_find; break; } }
			int k_time = K_time_array[1];
			int k_thres = Av[border];
			//cout << k_thres << endl;
			//cout << channel_index << endl;
			// if ((k < k_thres) && A_Acv[k]==1) { k_time = K_time_array[0]; }
			if (bfind && ((Kv -1 - channel_index) < border)) { k_time = K_time_array[0]; }
			K_time_all[k] = k_time;
			polar_codeconstruction(polarconh, Nh, k_time, GGh, (float)sqrt(sigma_sq), construction_time_all[k]);
			if (k_time == K_time_array[0]) { use_short_codeword[k] = 1; } // use long codeword
			else if (k_time == K_time_array[1]) { use_short_codeword[k] = 0; } // use short codeword
			for (int i = 0; i < k_time - 1; i++)
				for (int j = i + 1; j < k_time; j++)
					if (construction_time_all[k][i] > construction_time_all[k][j])
					{
						int ttt = construction_time_all[k][i];
						construction_time_all[k][i] = construction_time_all[k][j];
						construction_time_all[k][j] = ttt;
					}

			// Decide the position of information bits and frozen bits
			for (int i = 0; i < k_time; i++) // Information Set
			{
				A_time_all[k][i] = construction_time_all[k][i];
				A_Ac_time_all[k][A_time_all[k][i]] = 1;
			}
			//	getchar();
			for (int i = 0; i < Nh - k_time; i++) // Frozen Bit
			{
				Ac_time_all[k][i] = construction_time_all[k][k_time + i];
				A_Ac_time_all[k][Ac_time_all[k][i]] = 0;
			}
		}
		/*cout << endl;
		cout << "int channel_index_predefined[Nv] = {";
	    for (int i = 0; i < Nv; i++) { cout << use_short_codeword[i] << ","; }
		cout << "};" << endl;*/
		// Non-uniform channel
		vector<int> non_uniform_intrlv_pattern;
		non_uniform_intrlv_pattern.resize(Nv);
		if (non_uniform_channel)
		{
			getChannelScalingFactor(Nr, Nt, diff_gain, use_short_codeword, channel_scale_factor, predefined_non_uniform);
			if (predefined_non_uniform ) {
				generateIntrlvPatternIrregular(non_uniform_intrlv_pattern, Nv, bad_channel_index,
					good_channel_index, cnt_short_codewords, use_short_codeword);
				//// test
				/*cout << endl << "bad channel index: ";
				for (int i = 0; i < cnt_short_codewords; i++) { cout << bad_channel_index[i] << " "; }
				cout << endl << "good channel index: ";
				for (int i = 0; i < Nv-cnt_short_codewords; i++) { cout << good_channel_index[i] << " "; }
				cout << endl << "use short codeword: ";
			    for (int i = 0; i < Nv; i++) { cout << use_short_codeword[i] << " "; }
				cout << endl << "interleaving pattern: ";
				cout << "\n BAD" << endl;
				for (int i = 0; i < cnt_short_codewords; i++) { cout << bad_channel_index[i] << "(" << non_uniform_intrlv_pattern[bad_channel_index[i]] << ")" << "  "; }
				cout << "\n GOOD" << endl;
				for (int i = 0; i < Nv - cnt_short_codewords; i++) { cout << good_channel_index[i] << "(" << non_uniform_intrlv_pattern[good_channel_index[i]] << ")"<<"  "; }
				cout << endl;*/
			}
			/*if (predefined_non_uniform) { getChannelScalingFactor(Nr, Nt, diff_gain, use_short_codeword, channel_scale_factor); }
			else { getChannelScalingFactor(Nr, Nt, diff_gain, use_short_codeword, channel_scale_factor); }*/
		}
		// display A_Ac_time_all row by row
		/*for (int i = 0; i < Nv; i++)
		{
			for (int j = 0; j < Nh; j++)
			{
				cout << A_Ac_time_all[i][j] << ",";
			}
			cout << endl;
		}*/
		//// display K_time_all
		cout << endl;
		for (int i = 0; i < Nv; i++)
		{
			cout << K_time_all[i] << endl;
		}
		 //test: Horizontal frozen bits

		/*cout << "Frozen Bits" <<endl;
		for (int i = 0; i < Nv - Kv; i++)
		{
			cout << setw(5) << Acv[i];
		}
		cout << endl;*/
		// test: Horizontal frozen bits
		// Frame Loop
		int passnum;
		int frame = 0;

		// generate interleaving patterns
		if (flag_interleave)
		{
			string ofile_intrlv_name = string("intrlv_pattern_71_frame_49263.txt");
			string ofile_spatial_intrlv_name = string("intrlv_pattern_spatial71.txt");
			ifstream ofile_intrlv(ofile_intrlv_name);
			ifstream ofile_spatial_intrlv(ofile_spatial_intrlv_name);
			int file_in = 1;
			for (int i = 0; i < Nv; i++)
			{
				for (int j = 0; j < intrlv_pattern_num; j++)
				{

					for (int k = 0; k < ModType; k++)
						// ofile_intrlv >> file_in;
						ofile_intrlv >> intrlv_pattern_2D_time_partial[i][j][k];
					//cout << endl;
				}
			}
			ofile_intrlv.close();

			for (int i = 0; i < Nh; i++)
			{
				for (int j = 0; j < Nv; j++)
					ofile_spatial_intrlv >> intrlv_pattern_2D[i][j];
			}
			ofile_spatial_intrlv.close();
			/*for (int i = 0; i < Nh; i++)
			{
				for (int j = 0; j < Nv; j++)
					cout << intrlv_pattern_2D[i][j] << endl;
			}*/
			generateIntrlvPattern(intrlv_pattern, Nv);
			for (int i = 0; i < Nh; i++) { generateIntrlvPattern(intrlv_pattern_2D[i], Nv); }
			for (int i = 0; i < Nv; i++) { generateIntrlvPattern(intrlv_pattern_2D_time[i], Nh); }
			if (flag_interleave_partial)
			{
				for (int i = 0; i < Nv; i++)
				{
					for (int j = 0; j < intrlv_pattern_num; j++)
					{
						generateIntrlvPattern(intrlv_pattern_2D_time_partial[i][j], ModType);
						/*for (int k = 0; k < ModType; k++)
							cout << intrlv_pattern_2D_time_partial[i][j][k];
						cout << endl;*/
					}
				}
			}
			for (int i = 0; i < Nh; i++) { generateIntrlvPattern(intrlv_pattern_2D_inner[i], Kv); }
		}

		//ofstream ofile_intrlv("intrlv_pattern.txt");
		//for (int i=0; i<Nv; i++)
		//{
		//	for (int j = 0; j < intrlv_pattern_num; j+for+)
		//	{
		//		
		//		for (int k = 0; k < ModType; k++)
		//			ofile_intrlv << intrlv_pattern_2D_time_partial[i][j][k] << "   ";
		//		//cout << endl;
		//	}
		//}
		//ofile_intrlv.close();
		int frame_tmp = 0;
		int frame_delta = 0;
		int n_frame_error_pre = 0;
		while (Num_Frame_Error < errframe)
		{
			/*test: processing time*/
			//auto start = std::chrono::system_clock::now();
			/*if (frame == 10) {
				cout << "AAAAAAAAAA" << endl;
			}*/
			frame++;
			// generate interleaving pattern
			/*if (flag_interleave)
			{
				generateIntrlvPattern(intrlv_pattern, Nv);
				for (int i = 0; i < Nh; i++) { generateIntrlvPattern(intrlv_pattern_2D[i], Nv); }
				for (int i = 0; i < Nv; i++) { generateIntrlvPattern(intrlv_pattern_2D_time[i], Nh); }
				if (flag_interleave_partial)
				{
					for (int i = 0; i < Nv; i++)
					{
						for (int j = 0; j < intrlv_pattern_num; j++)
						{
							generateIntrlvPattern(intrlv_pattern_2D_time_partial[i][j], ModType);
						}
					}
				}
			}*/

			//cout << frame << endl;
			int u_input[16] = { 0,     1     ,0,     1,
									0,     1,     0,     0,
									0 ,    1,     0,     0,
									1 ,    0,     0,     0 };
			for (int i = 0; i < N_Info; i++)
			{
				a_data[i] = (unsigned char)(distrib(gen) % 2);
				//		a_data[i] = u_input[i];
			}

			// cout << N_Info << endl;
			if (CRCL) tx_append_crc(a_data, N_Info, CRCL);
			for (int i = 0; i < Nv; i++)
				for (int j = 0; j < Nh; j++)
				  u[i][j] = 0;
			//omp_set_num_threads(4);
//#pragma omp parallel for
			int info_count = 0;
			for (int i = 0; i < Kv; i++)
			{
				for (int j = 0; j < K_time_all[Av[i]]; j++)
				{
					u[Av[i]][A_time_all[Av[i]][j]] = (int)a_data[info_count];
					info_count++;
					//cout << (int)u[Av[i]][Ah[j]] << " \n";
					//cout << Av[i] << "\n";
				}
				//cout << endl; 
			}
			// test:original info
			//cout << "original info u:" << endl;
			////for (int j = 0; j < Nh; j++) { cout << u[Av[0]][j] << " "; }
			//for (int i = 0; i < Nv; i++)
			//{
			//	cout << endl;
			//	for (int j = 0; j < Nh; j++)
			//		cout << u[i][j];
			//}
			//cout << endl;
			// Polar Encoding
//#pragma omp parallel for

			// test: original info
			/*cout << "original info u:" << endl;
			display_array(u, Nv, Nh);*/
			// test: original info

			for (int i = 0; i < Nv; i++)
				PolarEncode_xor(midx[i], u[i], Nh);
			// test:horizontal encoding
			//cout << "u after horizontal encoding:" << endl;
			////for (int j = 0; j < Nh; j++) { cout << midx[Av[0]][j] << " "; }
			//for (int i = 0; i < Nv; i++)
			//{
			//	cout << endl;
			//	for (int j = 0; j < Nh; j++)
			//		cout << midx[i][j];
			//}
			//cout << endl;
			//cout << endl;
			// test:horizontal encoding
			// 
			// column-wise interleaving

			// test: row encoding
			/*cout << "u after horizontal encoding:" << endl;
			display_array(midx, Nv, Nh);*/
			// test: row encoding
			if (flag_interleave_inner)
			{
				int* tmp_to_interleave = new int[Nv];
				for (int i = 0; i < Nh; i++)
				{
					for (int j = 0; j < Nv; j++) { tmp_to_interleave[j] = midx[j][i]; }
					randIntrlvPartial(intrlv_pattern_2D_inner[i], tmp_to_interleave, Av, Nv, Kv);
					for (int j = 0; j < Nv; j++) { midx[j][i] = tmp_to_interleave[j]; }
				}
				delete[] tmp_to_interleave;
			}

			// test: innner interleaving
			/*cout << "u after inner interleaving:" << endl;
			display_array(midx, Nv, Nh);*/
			// test: innner interleaving
//#pragma omp parallel for
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

			// systematic spatial encoding
			// reset frozen bits to 0
			for (int i = 0; i < Nh; i++)
			{
				for (int j = 0; j < Nv - Kv; j++)
				{
					x[Acv[j]][i] = 0;
				}
			}
			// another round of spatial encoding
			for (int i = 0; i < Nh; i++)
			{
				int* tmpxv = new int[Nv]; int* tmpuv = new int[Nv];
				for (int j = 0; j < Nv; j++)
					tmpuv[j] = x[j][i];
				PolarEncode_xor(tmpxv, tmpuv, Nv);
				for (int j = 0; j < Nv; j++)
					x[j][i] = tmpxv[j];
				delete[] tmpxv; delete[] tmpuv;
			}



			// test: vertical encoding
			/*cout << "u after vertical encoding:" << endl;
			display_array(x, Nv, Nh);*/
			// test: vertical encoding
			/*for (int i = 0; i < 2 * Nr; i++)
				for (int j = 0; j < 2 * Nt; j++)
					H_real(i, j) = 0;*/
					/*MatrixXf H_real_test(2 * Nr, 2 * Nt);
					MatrixXf M = H_real_test.transpose();
					MatrixXf N = H_real_test;
					HTHIinv = M * N;*/
					//HTHIinv = M * H_real;
					// test:vertical encoding
					//cout << "u after vertical encoding:" << endl;
					////for (int j = 0; j < Nh; j++) { cout << x[Av[0]][j] << " "; }
					//for (int i = 0; i < Nv; i++)
					//{
					//	cout << endl;
					//	for (int j = 0; j < Nh; j++)
					//		cout << x[i][j];
					//}
					//cout << endl;
					//cout << endl;
					// test:vertical encoding
			if (atten == 2)
			{
				// modulation
				// #pragma omp parallel for
				// Interleaving
				if (flag_interleave)
				{
					if (flag_interleave_all)  // interleave temporally
					{
						for (int i = 0; i < Nv; i++)
						{
							randIntrlv(intrlv_pattern_2D_time[i], x[i], Nh);
						}
					}
					//// test: to be interleaved
					//cout << "x before interleaving:" << endl;
					//// print x line by line, add spaces between each 4 terms
					//for (int i = 0; i < Nv; i++)
					//{
					//	for (int j = 0; j < Nh; j++)
					//	{
					//		if (!(j % 4)) cout << " ";
					//		cout << x[i][j] << " ";
					//	}
					//	cout << endl;
					//}

					if (flag_interleave_partial)
					{
						for (int i = 0; i < Nv; i++)
						{
							for (int j = 0; j < intrlv_pattern_num; j++)
							{
								randIntrlv(intrlv_pattern_2D_time_partial[i][j], x[i] + j * ModType, ModType);
							}
						}
					}

					//// test: to be interleaved
					//cout << "x after interleaving:" << endl;
					//// print x line by line, add spaces between each 4 terms
					//for (int i = 0; i < Nv; i++)
					//{
					//	for (int j = 0; j < Nh; j++)
					//	{
					//		if (!(j % 4)) cout << " ";
					//		cout << x[i][j] << " ";
					//	}
					//	cout << endl;
					//}
					if ((!flag_interleave_block) && flag_interleave_spatial)
					{
						for (int j = 0; j < Nh; j++)
						{
							for (int i = 0; i < Nv; i++) { x_tmp[i] = x[i][j]; }
							randIntrlv(non_uniform_intrlv_pattern, x_tmp, Nv);
							for (int i = 0; i < Nv; i++) { x[i][j] = x_tmp[i]; }
						}
					}


					if (flag_interleave_block)
					{
						randIntrlvBlock(x, Nh, Nv, ModType, interleave_rows);
					}
				}
				if (MOD_DIM == "Spatial")
				{
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
						/*for (int i = 0; i < mimo_blk_height; i++)
						{
							cout << x_mod_tmp[i] << endl;
						}
						cout << "111111" << endl;*/
					}
					// AWGN Channel (received y_mod )
					if (!non_uniform_channel)H = generateChannelMatrix(Nr, Nt);
					else H = generateChannelMatrixNonuniform(Nr, Nt, channel_scale_factor);
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
					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < 2 * Nr; j++) { y_mod_real_temp_omp[i][j] = y_mod_real[j][i]; }
					}
					// for KBest Detection
					Eigen::HouseholderQR<Eigen::MatrixXf> qr_Hreal(H_real);
					Q = qr_Hreal.householderQ();
					R = qr_Hreal.matrixQR().triangularView<Eigen::Upper>();
					//omp_set_num_threads(8);

#pragma omp parallel for private(z)
					for (int j = 0; j < Nh; j++)
					{
						MatrixXf received_signal_omp(2 * Nr, 1);
						int k_kbest_extend = pow(2, ModType / 2) * K_KBEST;
						float** paths_kbest = new float* [k_kbest_extend];
						for (int i = 0; i < k_kbest_extend; i++) paths_kbest[i] = new float[2 * Nt];
						float* ED = new float[K_KBEST];

						if (MIMODET == "MMSE") { MMSEdetection(H_real, Identity, received_signal, mmse_matrix, output_signal, HTH, y_mod_real_temp_omp[j], y_mmse_temp_omp[j], sigma_sq); }
						// ִ��һЩ����

						for (int i = 0; i < 2 * Nr; i++) { received_signal_omp(i, 0) = y_mod_real_temp_omp[j][i]; }
						z = Q.transpose() * received_signal_omp;
						if (MIMODET == "KBEST")
						{
							KBestDetection(R, z, H_real, y_mod_real_temp_omp[j], ED, y_mmse_temp_omp[j], paths_kbest, K_KBEST, ModType, Nt, Nr);
							llr_kbest_bit(Nt, y_mmse_temp_omp[j], oriLLR_kbest[j], sigma_sq, ModType, K_KBEST, paths_kbest, ED);
						}
						for (int i_path = 0; i_path < k_kbest_extend; i_path++) { delete[] paths_kbest[i_path]; }

						delete[] paths_kbest;
						delete[] ED;
					}

					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < 2 * Nt; j++) { y_mmse[j][i] = y_mmse_temp_omp[i][j]; }
					}
					// sym llr to bitwise llr
					// Detection is done for every column, that is, for every Nv bits
					// Interleaving is not currently considered
					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < 2 * Nt; j++) { y_mmse_tmp[j] = y_mmse[j][i]; }
						if (MIMODET == "MMSE") { llr_mmse_bit(Nt, y_mmse_tmp, oriLLR_tmp, sigma_sq, ModType); }
						if (MIMODET == "KBEST") { for (int j = 0; j < len_ori_llr_tmp; j++) { oriLLR_tmp[j] = oriLLR_kbest[i][j]; } }
						if (flag_interleave)
						{
							randDeintrlv(intrlv_pattern, oriLLR_tmp, Nv);
						}
						bool flagIsMax = false;
						for (int j = 0; j < Nv; j++)
						{
							if (!(j % 2))
							{
								if (abs(x_mmse_out[j][i]) > abs(x_mmse_out[j][i + 1])) { flagIsMax = true; }
								else { flagIsMax = false; }
							}
							else
							{
								if (abs(x_mmse_out[j][i]) > abs(x_mmse_out[j][i - 1])) { flagIsMax = true; }
								else { flagIsMax = false; }
							}
							oriLLRh[j][i] = oriLLR_tmp[j];
							upLLR[j][i] = oriLLR_tmp[j];
							// TEST�� det error
							x_mmse_out[j][i] = 0;
							if (oriLLR_tmp[j] < 0) { x_mmse_out[j][i] = 1; }
							if (x_mmse_out[j][i] != x[j][i]) {
								n_det_err++;
								// error at max-reliable bit. Warning: Only fit for 16QAM modulation
								if (flagIsMax) { n_det_err_mrb++; }
							}
							total_bit_cnt++;
						}
					}
					// test : detection TIME
					/*det_end = std::chrono::system_clock::now();
					auto duration_detection = std::chrono::duration_cast<std::chrono::microseconds>(det_end - start);
					file_det << duration_detection.count() << endl;*/
				}
				// Temporal modulation
				else if (MOD_DIM == "Temporal")
				{
					//modulation
					//cout << "START TEST" << endl;
					for (int i = 0; i < Nv; i++) // every row
					{
						modulation(x[i], x_mod[i], ModType, mimo_blk_length);
					}
					/*for (int i = 0; i < mimo_blk_height; i++)
					{
						for (int j = 0; j < mimo_blk_length; j++)
							cout << x_mod[i][j] << endl;
					}
					cout << "111111" << endl;*/
					//AWGN channel (received y_mod)
					// test: modulated symbols
					/*cout << endl;
					cout << "modulated symbols" << endl;
					for (int i = 0; i < Nt; i++) cout << x_mod[i][0] << endl;
					cout << "modulated symbols" << endl;*/
					// test: modulated symbols
					// H = generateChannelMatrix(Nr, Nt);
					if (!non_uniform_channel)H = generateChannelMatrix(Nr, Nt);
					else H = generateChannelMatrixNonuniform(Nr, Nt, channel_scale_factor);
					for (int j = 0; j < Nsym; j++)
					{
						// Extracting a single column
						for (int i = 0; i < Nv; i++) { x_mod_tmp[i] = x_mod[i][j]; }
						AWGN(Nr, Nt, x_mod_tmp, y_mod_tmp, H, sigma_sq);
						for (int i = 0; i < Nr; i++) { y_mod[i][j] = y_mod_tmp[i]; }
					}
					convertHToReal(H, H_real);

					for (int j = 0; j < Nsym; j++)
					{
						for (int i = 0; i < Nr; i++) { y_mod_tmp[i] = y_mod[i][j]; }
						convertYToReal(Nr, y_mod_tmp, y_mod_real_tmp);
						for (int i = 0; i < 2 * Nr; i++) { y_mod_real[i][j] = y_mod_real_tmp[i]; }
					}
					//MMSE detection done in paralell
					for (int i = 0; i < Nsym; i++)
					{
						for (int j = 0; j < 2 * Nr; j++) { y_mod_real_temp_omp[i][j] = y_mod_real[j][i]; }
					}
					Eigen::HouseholderQR<Eigen::MatrixXf> qr_Hreal(H_real);
					Q = qr_Hreal.householderQ();
					R = qr_Hreal.matrixQR().triangularView<Eigen::Upper>();
					// for EP detection
					// copy H_real
					for (int i = 0; i < Nr2; i++)
					{
						for (int j = 0; j < Nt2; j++)
						{
							Hreal_array[i * Nt2 + j] = H_real(i, j);
						}
					}
					// calculate HTH
					cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, Nt2, Nt2, Nr2, 1.0, Hreal_array, Nt2, Hreal_array, Nt2, 0.0, HTH_array, Nt2);
					//cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, Nt2, 1, Nr2, 1.0, ch_mtx_struct.H_real[cnt], Nt2, data_struct.rx_buffer_real[cnt], 1, 0.0, ch_mtx_struct.hty, 1);

					//omp_set_num_threads(8);
					//cout << 1 << endl;
					//omp_set_num_threads(4);
					/*auto start = std::chrono::system_clock::now();*/
					/*MatrixXf H_real_trans(2*Nt, 2*Nr);
					for (int i = 0; i < 2 * Nt; i++)
						for (int j = 0; j < 2 * Nr; j++)
							H_real_trans(i, j) = H_real(j, i);
					HTHIinv = (H_real_trans * H_real + Identity).inverse();*/
# pragma omp parallel for private(z)
					for (int j = 0; j < Nsym; j++)

						//for (int x=0; x<Nsym; x++)
					{
						// int j = x % Nsym;
						MatrixXf received_signal_omp(2 * Nr, 1);
						int k_kbest_extend = pow(2, ModType / 2) * K_KBEST;
						float** paths_kbest = new float* [k_kbest_extend];
						for (int i = 0; i < k_kbest_extend; i++) { paths_kbest[i] = new float[2 * Nt]; }
						float* ED = new float[K_KBEST];
						float* HTY_array_tmp = new float[2 * Nt];
						// Declared arrays for EP
						int flag = 1;
						int temp = 0;
						//float* extr_mean = new float[2 * Nr];
						//float* extr_var = new float[2 * Nr];
						double* vec_alpha_in = new double[2 * Nr];
						double* vec_gamma_in = new double[2 * Nr];
						double* p_symbol_rearrange = new double[Nt2 * pow(2, ModType / 2)];
						for (int i = 0; i < Nt2 * pow(2, ModType / 2); i++) { p_symbol_rearrange[i] = 1 / (pow(2, ModType / 2)); }
						for (int i = 0; i < Nr2; i++)
						{
							//extr_mean[i] = 0;
							//extr_var[i] = 0;
							vec_alpha_in[i] = 2;
							vec_gamma_in[i] = 0;
						}

						if (MIMODET == "MMSE")
						{
							MMSEdetection(H_real, Identity, received_signal, mmse_matrix, output_signal, HTH, y_mod_real_temp_omp[j], y_mmse_temp_omp[j], sigma_sq);
							//MMSEdetection_noInverse(H_real, Identity, received_signal, mmse_matrix, output_signal, HTH, HTHIinv, y_mod_real_temp_omp[j], y_mmse_temp_omp[j], sigma_sq);
						}
						if (MIMODET == "KBEST")
						{
							for (int i = 0; i < 2 * Nr; i++) { received_signal_omp(i, 0) = y_mod_real_temp_omp[j][i]; }
							z = Q.transpose() * received_signal_omp;
							KBestDetection(R, z, H_real, y_mod_real_temp_omp[j], ED, y_mmse_temp_omp[j], paths_kbest, K_KBEST, ModType, Nt, Nr);
							llr_kbest_bit(Nt, y_mmse_temp_omp[j], oriLLR_kbest[j], sigma_sq, ModType, K_KBEST, paths_kbest, ED);
						}
						if (MIMODET == "EP")
						{
							// test: save H_real, recv_signal
							/*save_mat_txt(y_mod_real_temp_omp[j], Nr2, "RxSig0128.txt");
							save_mat_txt(H_real, "Hreal0128.txt");*/
							// test: save H_real, recv_signal
							cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, Nt2, 1, Nr2, 1.0, Hreal_array, Nt2, y_mod_real_temp_omp[j], 1, 0.0, HTY_array_tmp, 1);
							//test:HTH_array
							/*for (int i = 0; i < Nt2; i++)
								for (int s = 0; s < Nt2; s++)
								{
									if (s == 0) std::cout << endl;
									std::cout << setw(10) << HTH_array[i * Nt2 + s];
								}
							cout << endl;*/
							// TEST:HTY
							/*for (int i = 0; i < Nr2; i++)
								received_signal_omp(i, 0) = y_mod_real_temp_omp[j][i];
							MatrixXf HTransY = H_real.transpose() * received_signal_omp;
							cout << "Eigen" << endl;
							for (int i = 0; i < Nt2; i++) cout << setw(10)<<HTransY(i, 0);
							cout << "mkl" << endl;
							for (int i = 0; i < Nt2; i++) cout << HTY_array[i];*/
							// TEST:HTY
							EPD_detection_float(flag, Nt, Nr, iter_MIMO, ModType, 1.0, sigma_sq, HTH_array, HTY_array_tmp, temp, ep_extr_mean_omp[j], ep_extr_var_omp[j], vec_alpha_in, vec_gamma_in, p_symbol_rearrange);
						}
						//test EP detection
						/*cout << "TEST:EP DETECTION" << endl;
						for (int i = 0; i < Nr2; i++) std::cout << extr_mean[i]<< endl;
						cout << "TEST:EP DETECTION" << endl;*/
						//test EP detection
						for (int i_path = 0; i_path < k_kbest_extend; i_path++) { delete[] paths_kbest[i_path]; }
						delete[] paths_kbest;
						delete[] ED;
						//delete[] extr_mean;
						//delete[] extr_var;
						delete[] vec_alpha_in;
						delete[] vec_gamma_in;
						delete[] p_symbol_rearrange;
						delete[] HTY_array_tmp;
					}
					// test: detection time
					/*auto detection_over = std::chrono::system_clock::now();
					auto duration_detection = std::chrono::duration_cast<std::chrono::microseconds>(detection_over - start);
					cout << endl;
					cout <<"Detection duration" << duration_detection.count() << "us" << endl;*/
					// test: detection time
					for (int i = 0; i < Nsym; i++)
					{
						for (int j = 0; j < 2 * Nt; j++) { y_mmse[j][i] = y_mmse_temp_omp[i][j]; }
					}
					//To bit domain output
					for (int j = 0; j < Nsym; j++)
					{
						for (int i = 0; i < 2 * Nt; i++) { y_mmse_tmp[i] = y_mmse[i][j]; }
						// With Spatial Modulation, Nt = Nv
						if (MIMODET == "MMSE") { llr_mmse_bit(Nt, y_mmse_tmp, oriLLR_tmp, sigma_sq, ModType); }
						if (MIMODET == "EP") { llr_ep_bit(Nt, ep_extr_mean_omp[j], ep_extr_var_omp[j], oriLLR_tmp, ModType); }
						if (MIMODET == "KBEST") { for (int i = 0; i < len_ori_llr_tmp; i++) { oriLLR_tmp[i] = oriLLR_kbest[j][i]; } }
						int i_bit = 0;
						int bit_col = 0;
						for (int j_inner = 0; j_inner < Nv; j_inner++)
						{
							for (int i_bit_col = 0; i_bit_col < ModType; i_bit_col++)
							{
								i_bit = ModType * j_inner + i_bit_col;
								bit_col = j * ModType + i_bit_col;
								oriLLRh[j_inner][bit_col] = oriLLR_tmp[i_bit];
								upLLR[j_inner][bit_col] = oriLLR_tmp[i_bit];
							}

						}
					}
					// test: LLR transfer time
					/*auto llr_to_bit_over = std::chrono::system_clock::now();
					auto duration_llr = std::chrono::duration_cast<std::chrono::microseconds>(llr_to_bit_over - detection_over);
					cout << endl;
					cout<< "LLR transfer duration" << duration_llr.count() << "us" << endl;*/
					// test: LLR transfer time
					if (flag_interleave)
					{
						if ((!flag_interleave_block) && flag_interleave_spatial)
						{
							for (int j = 0; j < Nh; j++)
							{
								for (int i = 0; i < Nv; i++) { oriLLR_tmp[i] = oriLLRh[i][j]; }
								randDeintrlv(non_uniform_intrlv_pattern, oriLLR_tmp, Nv);
								for (int i = 0; i < Nv; i++)
								{
									oriLLRh[i][j] = oriLLR_tmp[i];
									upLLR[i][j] = oriLLR_tmp[i];
								}

							}
						}

						if (flag_interleave_block)
						{
							randDeIntrlvBlock(oriLLRh, Nh, Nv, ModType, interleave_rows);
							randDeIntrlvBlock(upLLR, Nh, Nv, ModType, interleave_rows);
						}

						if (flag_interleave_all) // Deinterleave in time domain
						{
							for (int j = 0; j < Nv; j++)
							{
								randDeintrlv(intrlv_pattern_2D_time[j], oriLLRh[j], Nh);
								randDeintrlv(intrlv_pattern_2D_time[j], upLLR[j], Nh);
							}
						}

						if (flag_interleave_partial)
						{
							for (int j = 0; j < Nv; j++)
							{
								for (int k = 0; k < intrlv_pattern_num; k++)
								{
									randDeintrlv(intrlv_pattern_2D_time_partial[j][k], oriLLRh[j] + k * ModType, ModType);
									randDeintrlv(intrlv_pattern_2D_time_partial[j][k], upLLR[j] + k * ModType, ModType);
								}
							}
						}

					}

					for (int i = 0; i < Nv; i++)
					{
						for (int j = 0; j < Nh; j++)
						{
							total_bit_cnt++;
							int bit_dec = (oriLLRh[i][j] < 0);
							if (bit_dec != x[i][j]) n_det_err++;
						}
					}

				}
				// test: detection errors

				// test: detection errors
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
			for (int iter = 1; iter <= softiter; iter++)
			{
				//vertical decoding
		   // test: decoding time
			//auto decoding_start = std::chrono::system_clock::now();

			// test: decoding time
#pragma omp parallel for
				for (int decodingi = 0; decodingi < Nh; decodingi++) // 291.70kb
				{

					//for SCL Decoding
					int Counth = 0, Countinfoh = 0;
					float Qsumunh = 0;
					int** sum = new int* [L];
					//for (int i = 0; i < L; i++)  sum[i] = new  int[Nmax];
					int** u1 = new int* [L];
					for (int i = 0; i < L; i++)  u1[i] = new _MM_ALIGN16 int[Nmax];
					for (int i = 0; i < L; i++)  sum[i] = new _MM_ALIGN16 int[Nmax];
					//for (int i = 0; i < L; i++)  u1[i] = new  int[Nmax];

					float* LLR_in = new float[L];
					float* W = new float[2 * L];
					int* SCL_index = new int[L];
					int* better = new int[L];
					int* worse = new int[L];
					int* indexpm = new int[L];
					float* PM = new float[L];
					float** LLR = new float* [L];
					for (int i = 0; i < L; i++)  LLR[i] = new _MM_ALIGN16 float[2 * Nmax + 2];
					// for (int i = 0; i < L; i++)  LLR[i] = new float[2 * Nmax + 2];
					////////////////////////////////
					for (int j = 0; j < Nv; j++)
						for (int k = 0; k < L; k++)
							LLR[k][j] = upLLR[j][decodingi];
					//cout << LLR[0][j] << " ";

				//cout << endl;
					int index = 0;
					/*			if (L == 1)
								{
									//T_start = clock();
									Count = 0;
									decode(*LLR, *sum, (int)(log(Nv) / log(2)), A_Acv, flag, Nv, Kv);
									//PolarEncode_xor(*u1, *sum, N);
									//T_finish = clock();
								}
								else
								{*/
					for (int k = 0; k < L; k++) { PM[k] = 0; indexpm[k] = k; }
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
					if (adaptive_beta)
					{
						if (softmethod == "Pyndiah") Pyndiahv_adaptivebeta(LLR, sum, oriLLRh, indexpm, iter * 2 - 2, decodingi, upLLR);
						if (softmethod == "MITSO") MITSOv(LLR, PM, sum, oriLLRh, indexpm, iter * 2 - 2, decodingi, upLLR, Qsumunh);
					}
					// iter * 2 - 2
					else
					{
						//cout << iter << endl;
						if (softmethod == "Pyndiah") Pyndiahv(LLR, sum, oriLLRh, indexpm, iter, decodingi, upLLR);
						if (softmethod == "MITSO") MITSOv(LLR, PM, sum, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunh);
					}
					// test: oriLLR and upLLR
					/*cout << endl;
					for (int i = 0; i < Nv; i++)
						for (int j = 0; j < Nh; j++)
							cout << oriLLRh[i][j] - upLLR[i][j] << "   ";
					cout << endl;*/
					/*cout << oriLLRh[0][0] << endl;
					cout << upLLR[0][0] << endl;*/
					// test: oriLLR and upLLR
					// half iteration Enabled
					// May Need Modification
					if (half_iter && iter == softiter)
					{
						if (iter == softiter) PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), Nh);
						for (int j = 0; j < Kv; j++)
							umid[Av[j]][decodingi] = u1[index][Av[j]];
					}
					// half iteration
					////////delete SCL variables
					for (int i = 0; i < L; i++)
						delete[]LLR[i];
					delete[] LLR;
					for (int i = 0; i < L; i++)
						delete[]sum[i];
					for (int i = 0; i < L; i++)
						delete[]u1[i];
					delete[] u1;
					delete[] sum;
					delete[] PM;
					delete[] LLR_in;
					delete[] W;
					delete[] SCL_index;
					delete[] better;
					delete[] worse;
					delete[] indexpm;

				}
				//cout << "Horizon Start:\n";
				//horizontal decoding
				// #
				if (weirditer) { continue; }

				// test: LLR after horizontal decoding
				/*cout << "After horizontal decoding" << endl;
				for (int i = 0; i < Nv; i++)
				{
					for (int j = 0; j < Nh; j++)
					{
						cout << (oriLLRh[i][j] < 0) << " ";
					}
					cout << endl;
				}*/
				// test: LLR after horizontal decoding
				// inner deinterleaving
				if (flag_interleave_inner)
				{
					float* tmp_llr_to_interleave = new float[Nv];
					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < Nv; j++) { tmp_llr_to_interleave[j] = oriLLRh[j][i]; }
						randDeIntrlvPartial(intrlv_pattern_2D_inner[i], tmp_llr_to_interleave, Av, Nv, Kv);
						for (int j = 0; j < Nv; j++) { oriLLRh[j][i] = tmp_llr_to_interleave[j]; }
					}
					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < Nv; j++) { tmp_llr_to_interleave[j] = upLLR[j][i]; }
						randDeIntrlvPartial(intrlv_pattern_2D_inner[i], tmp_llr_to_interleave, Av, Nv, Kv);
						for (int j = 0; j < Nv; j++) { upLLR[j][i] = tmp_llr_to_interleave[j]; }
					}
					delete[] tmp_llr_to_interleave;
				}


				/*auto decoding_over = std::chrono::system_clock::now();
				auto duration_decoding = std::chrono::duration_cast<std::chrono::microseconds>(decoding_over - decoding_start);
				cout << endl;
				cout << "Decoding duration" << duration_decoding.count() << "us" << endl;*/

#pragma omp parallel for
				for (int decodingi = 0; decodingi < Nv; decodingi++)
				{
					// 293.70kb
					// cout << "KKKK" << endl;
					if (half_iter && (iter == softiter)) { break; }
					string polarde_now = polarde_array[iter - 1];
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

					/*			if (L == 1)
								{
									//T_start = clock();
									Count = 0;
									decode(*LLR, *sum, (int)(log(Nh) / log(2)), A_Ach, flag, Nh, Kh);
									if (iter == softiter) PolarEncode_xor(*u1, *sum, Nh);
									//T_finish = clock();
								}
								else
								{*/
								//for (int k = 0; k < L; k++)  PM[k] = 0, indexpm[k] = k; 
								////T_start = clock();
								//decode_list(LLR, sum, (int)(log(Nh) / log(2)), A_Ach, 0, Nh, Kh, PM, L, (int)(log(L) / log(2)), LLR_in, W, SCL_index, better, worse, &Countv, &Countinfov, &Qsumunv);
								////T_finish = clock();
								////Decoding_Duration += (T_finish - T_start);
								//for (int i = 0; i < L - 1; i++)
								//	for (int j = i + 1; j < L; j++)
								//		if (PM[indexpm[i]] < PM[indexpm[j]])
								//		{
								//			int ttt = indexpm[i];
								//			indexpm[i] = indexpm[j];
								//			indexpm[j] = ttt;
								//		}


							/*
							* Encoding: Info -> Horizontal Encoding -> A -> set to 0 -> B -> Horizontal Encoding -> Vertical Encoding
							* Decoding: Received Info -> Horizontal Decoding -> Vertical Decoding -> Horizontal Decoding
							*/
					if (polarde_now == "SLD")
					{
						LSD_decode_outall(u_decod, u1_LSD, LLR[0], GGh, A_Ach, Nh, L_LSD, n_paths);
						// if (iter == softiter) PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), Nh);
						// for soft information exchange, inverse coding:
						// THIS LINE is WRONG! BUT IT SEEMS BETTER THAN THE RIGHT VERSION
						for (int i = 0; i < L_LSD; i++) { PolarEncode_xor(sum_LSD[i], u1_LSD[i], Nh); }
						// for (int i = 0; i < L_LSD; i++) { PolarEncode_xor(sum_LSD[i], u1_LSD[i], Nh); }
						// if (iter == softiter) PolarEncode_xor(u1[0], sum_LSD[0], Nh);
						index = indexpm[0];
						// test: after horizontal decoding
						/*if (decodingi == Av[0])
						{
							cout << "Direct output of horizontal decoding:" << endl;
							for (int i = 0; i < Nh; i++) { cout << u1_LSD[0][i] << " "; }
							cout << endl;
							cout << "horizontal decoder out, encode once:" << endl;
							for (int i = 0; i < Nh; i++) { cout << sum_LSD[0][i] << " "; }
							cout << endl;
						}*/
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
						// iter*2-1
						if (iter < softiter)
						{
							if (adaptive_beta)
							{
								if (softmethod == "Pyndiah") Pyndiahh_LSD_adaptivebeta(LLR, sum_LSD, oriLLRh, indexpm, iter, decodingi, upLLR);
								if (softmethod == "MITSO") MITSOh(LLR, PM, sum_LSD, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunv);
							}
							else if (irregular_iteration)
							{

							}
							else
							{
								// cout << "IN LSD!!!" << endl;
								if (softmethod == "Pyndiah")
									Pyndiahh_LSD_getsum(LLR, sum_LSD, oriLLRh, indexpm, iter, decodingi, upLLR, sum_sdis_array[decodingi], GGh, Ach);
								//if (softmethod == "Pyndiah") Pyndiahh_LSD(LLR, sum_LSD, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR);
								if (softmethod == "MITSO") MITSOh(LLR, PM, sum_LSD, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunv);
							}
							//	getchar();
						}
						else
							for (int j = 0; j < Kh; j++)
								umid[decodingi][Ah[j]] = u1_LSD[index][Ah[j]];
					}
					if (polarde_now == "SCL")
					{
						for (int k = 0; k < L; k++)  PM[k] = 0, indexpm[k] = k;
						//T_start = clock();
						// decode_list(LLR, sum, (int)(log(Nh) / log(2)), A_Ach, 0, Nh, Kh, PM, L, (int)(log(L) / log(2)), LLR_in, W, SCL_index, better, worse, &Countv, &Countinfov, &Qsumunv);
						decode_list(LLR, sum, (int)(log(Nh) / log(2)), A_Ac_time_all[decodingi], 0, Nh, K_time_all[decodingi], PM, L, (int)(log(L) / log(2)), LLR_in, W, SCL_index, better, worse, &Countv, &Countinfov, &Qsumunv);
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
						if (iter == softiter) PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), Nh);
						// since systematic encoding is used in space domain, we directly copy the best path (x) to u1 and encode on time domain when all temporal decoing is done
						/*if (iter == softiter)
						{
							for (int j = 0; j < Nh; j++) { u1[indexpm[0]][j] = sum[indexpm[0]][j]; }
						}*/
						index = indexpm[0];
						//}
						//calc distance s()
						if (iter < softiter)
						{
							if (adaptive_beta)
							{
								if (softmethod == "Pyndiah") Pyndiahh_adaptivebeta(LLR, sum, oriLLRh, indexpm, iter, decodingi, upLLR);
								if (softmethod == "MITSO") MITSOh(LLR, PM, sum, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunv);
								//	getchar();
							}
							else if (irregular_coding)
							{
								if (softmethod == "Pyndiah") Pyndiahh_irregular(LLR, sum, oriLLRh, indexpm, iter, decodingi, upLLR, sum_sdis_array[decodingi],use_short_codeword);
								if (softmethod == "MITSO") MITSOh_irregular(LLR, PM, sum, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunv,use_short_codeword);
							}
							else
							{
								//if (softmethod == "Pyndiah") Pyndiahh(LLR, sum, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR);
								if (softmethod == "Pyndiah") Pyndiahh_getsum(LLR, sum, oriLLRh, indexpm, iter, decodingi, upLLR, sum_sdis_array[decodingi]);
								if (softmethod == "MITSO") MITSOh(LLR, PM, sum, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunv);
								//	getchar();
							}
						}
						else
							for (int j = 0; j < Nh; j++)
								umid[decodingi][j] = u1[index][j];
					}

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


				if (flag_interleave_inner)
				{
					float* tmp_llr_to_interleave = new float[Nv];
					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < Nv; j++) { tmp_llr_to_interleave[j] = oriLLRh[j][i]; }
						randIntrlvPartial(intrlv_pattern_2D_inner[i], tmp_llr_to_interleave, Av, Nv, Kv);
						for (int j = 0; j < Nv; j++) { oriLLRh[j][i] = tmp_llr_to_interleave[j]; }
					}
					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < Nv; j++) { tmp_llr_to_interleave[j] = upLLR[j][i]; }
						randIntrlvPartial(intrlv_pattern_2D_inner[i], tmp_llr_to_interleave, Av, Nv, Kv);
						for (int j = 0; j < Nv; j++) { upLLR[j][i] = tmp_llr_to_interleave[j]; }
					}
					delete[] tmp_llr_to_interleave;
				}
			}
			int calcsum[Nv], calcu1[Nv];
			int calcu1h[Nh];
			// calcsum is from the direct output of the decoder
			// if output is the result of row decoding (meaning there is another column decoding to be done)
			if (half_iter)
			{
				for (int i = 0; i < Kv; i++)
				{
					PolarEncode_xor(calcu1h, umid[Av[i]], Nh);
					for (int j = 0; j < Kh; j++)
						uout[Av[i]][Ah[j]] = calcu1h[Ah[j]];
				}
			} 
			// if output is the result of column decoding (meaning that there is another row decoding to be done)
			else
			{
				for (int i = 0; i < Kh; i++)
				{
					for (int j = 0; j < Nv; j++)
						calcsum[j] = umid[j][Ah[i]];
					PolarEncode_xor(calcu1, calcsum, Nv);
					for (int j = 0; j < Kv; j++)
						uout[Av[j]][Ah[i]] = calcu1[Av[j]];
				}
			}
			// test: final output
			/*cout << "final output" << endl;
			for (int k = 0; k < Nh; k++) { cout << uout[Av[0]][k] << " "; }
			cout << endl;*/
			// test: final output
			//now without CRC
			int Temp_Error = Num_Error;
			for (int i = 0; i < Kv; i++)
			{
				//if ((i == 24) || (i == 56)) { continue; }
				for (int j = 0; j < K_time_all[Av[i]]; j++)
				{
					//printf("%d ", (int)uout[Av[i]][Ah[j]]);
					Num_Error += abs(u[Av[i]][A_time_all[Av[i]][j]] - umid[Av[i]][A_time_all[Av[i]][j]]);
					if (abs(u[Av[i]][A_time_all[Av[i]][j]] - umid[Av[i]][A_time_all[Av[i]][j]]) == 1)
					{
						all_error_positions[Av[i]][A_time_all[Av[i]][j]] += 1;
					}
				}
				//cout << endl;
			}
			//getchar();
			if (Temp_Error != Num_Error) Num_Frame_Error++;
			if (frame % 10 == 0)
			{
				/*printf("\r<%d> SNR= %.2f FER= %d / %d = %f BER= %d / %d = %lf E-Det = %f, E-MRB = %f"  ,
					frame, nSNR, Num_Frame_Error, frame, (double)Num_Frame_Error / frame, Num_Error, (N_Info * frame), (double)Num_Error / N_Info / frame,
					(double)n_det_err/total_bit_cnt, (double)n_det_err_mrb/(total_bit_cnt/(ModType/2)));*/
				int L_sdis_sum = L_LSD;
				sdis_sum = 0;
				for (int i = 0; i < Nv; i++) sdis_sum += sum_sdis_array[i];
				if (polarde_array[0] == "SCL") L_sdis_sum = L;
				// cout << "sdis_sum_avg" << sdis_sum / frame / L_sdis_sum << endl;
				printf("\r<%d> SNR= %.2f FER= %d / %d = %f BER= %d / %d = %lf ,sdis_avg = %f, E-MRB = %f",
					frame, nSNR, Num_Frame_Error, frame, (double)Num_Frame_Error / frame, Num_Error, (N_Info * frame), (double)Num_Error / N_Info / frame, (double)sdis_sum / frame / L_sdis_sum, (double)n_det_err / total_bit_cnt);
				fflush(stdout);
			}


			// output interleaving pattern
			//if (Num_Frame_Error % record_frames == 0 && Num_Frame_Error > 0 && Num_Frame_Error!=n_frame_error_pre)
			//{
			//	n_frame_error_pre = Num_Frame_Error;
			//	cout << endl;
			//	int n_file = Num_Frame_Error / record_frames;
			//	frame_delta = frame - frame_tmp;
			//	cout << "Delta Frames: " << frame_delta << endl;
			//	frame_tmp = frame;
			//	string ofile_intrlv_name = string("intrlv_pattern_") + to_string(n_file) + string("_frame_") + to_string(frame_delta) + string(".txt");
			//	string ofile_spatial_intrlv_name = string("intrlv_pattern_spatial") + to_string(n_file) + string(".txt");
			//	ofstream ofile_intrlv(ofile_intrlv_name);
			//	ofstream ofile_spatial_intrlv(ofile_spatial_intrlv_name);
			//	for (int i = 0; i < Nv; i++)
			//	{
			//		for (int j = 0; j < intrlv_pattern_num; j++)
			//		{

			//			for (int k = 0; k < ModType; k++)
			//				ofile_intrlv << intrlv_pattern_2D_time_partial[i][j][k] << " ";
			//			//cout << endl;
			//		}
			//	}
			//	ofile_intrlv.close();

			//	for (int i = 0; i < Nh; i++)
			//	{
			//		for (int j = 0; j < Nv; j++)
			//			ofile_spatial_intrlv << intrlv_pattern_2D[i][j] << " ";
			//	}
			//	ofile_spatial_intrlv.close();

			//	// re-generate interleaving pattern
			//	generateIntrlvPattern(intrlv_pattern, Nv);
			//	for (int i = 0; i < Nh; i++) { generateIntrlvPattern(intrlv_pattern_2D[i], Nv); }
			//	for (int i = 0; i < Nv; i++) { generateIntrlvPattern(intrlv_pattern_2D_time[i], Nh); }
			//	if (flag_interleave_partial)
			//	{
			//		for (int i = 0; i < Nv; i++)
			//		{
			//			for (int j = 0; j < intrlv_pattern_num; j++)
			//			{
			//				generateIntrlvPattern(intrlv_pattern_2D_time_partial[i][j], ModType);
			//				/*for (int k = 0; k < ModType; k++)
			//					cout << intrlv_pattern_2D_time_partial[i][j][k];
			//				cout << endl;*/
			//			}
			//		}
			//	}
			//	//Num_Frame_Error = 0;
			//}


			//printf("%d %d", Num_Frame_Error, Num_Error);
			//getchar();
			//Decoding_Duration += (T_finish - T_start);
			/*auto dec_end = std::chrono::system_clock::now();
			auto duration_decoding = std::chrono::duration_cast<std::chrono::microseconds>(dec_end - det_end);
			auto duration_total = std::chrono::duration_cast<std::chrono::microseconds>(dec_end - start);
			file_dec << duration_decoding.count() << endl;
			file_total << duration_total.count() << endl;*/
			if (frame % 16 == 0)
			{
				file_err_num << Num_Error << endl;
			}
			if (frame % 2000 == 0)
			{
				std::ofstream ofile_err("all_errors_log.txt");
				for (int i = 0; i < Nv; i++)
				{
					for (int j = 0; j < Nmax; j++)
					{
						ofile_err << all_error_positions[i][j] << " ";
					}
					ofile_err << endl;
				}
				ofile_err << endl;
				ofile_err.close();
			}

		}
		//Decoding_Duration = Decoding_Duration * 100 / CLOCKS_PER_SEC / frame;
		//Decoding_Duration = Decoding_Duration / CLOCKS_PER_SEC / Max_frame;

		//printf("%f  %f  %f  %f\n", nSNR, Decoding_Duration, (double)Num_Error / N_Info / Max_frame, (double)Num_Frame_Error / Max_frame);
		printf("%f %f %f\n", nSNR, (double)Num_Frame_Error / frame, (double)Num_Error / N_Info / frame);
		BER_array[SNR_count] = (double)Num_Error / N_Info / frame;
		FER_array[SNR_count] = (double)Num_Frame_Error / frame;
		// save to .txt file
		string polarde_serial = "";
		for (int i = 0; i < softiter; i++)
			polarde_serial = polarde_serial + polarde_array[i] + string("_");
		string file_name = string("ResFinal\\Irregular_") + polarde_serial + MIMODET + string("T") + to_string(Nt) + string("R") + to_string(Nr)
			+ string("Q") + to_string(int(pow(2, ModType))) + string("L") + to_string(Nh) + string("K") + to_string(Kh)
			+ string("T") + to_string(softiter) + string("H") + to_string(half_iter) + string("UB") + to_string(usebeta) + string("IS") + to_string(flag_interleave_spatial)
			+ string("DG") + to_string(diff_gain) + string(".txt");
		save_array_txt(nSNR, FER_array[SNR_count], file_name);
		int L_sdis_sum = L_LSD;
		sdis_sum = 0;
		for (int i = 0; i < Nv; i++) sdis_sum += sum_sdis_array[i];
		if (polarde_array[0] == "SCL") L_sdis_sum = L;
		// cout << "sdis_sum_avg" << sdis_sum / frame / L_sdis_sum << endl;
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
	for (int i = 0; i < Nh; i++)
		delete[] oriLLR_kbest[i];
	delete[] oriLLR_kbest;
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
	for (int i = 0; i < mimo_blk_height; i++) {
		delete[] x_mod[i];
	}
	delete[] x_mod;
	for (int i = 0; i < Nr; i++) {
		delete[] y_mod[i];
	}
	delete[] y_mod;
	for (int i = 0; i < 2 * mimo_blk_height; i++) {
		delete[] y_mmse[i];
	}
	delete[] y_mmse;
	for (int i = 0; i < 2 * Nr; i++) {
		delete[] y_mod_real[i];
	}
	for (int i = 0; i < mimo_blk_length; i++)
	{
		delete[] y_mod_real_temp_omp[i];
	}
	delete[] y_mod_real_temp_omp;
	for (int i = 0; i < mimo_blk_length; i++)
	{
		delete[] y_mmse_temp_omp[i];
	}
	delete[] y_mmse_temp_omp;
	for (int i = 0; i < mimo_blk_length; i++)
	{
		delete[] ep_extr_mean_omp[i];
	}
	delete[] ep_extr_mean_omp;
	for (int i = 0; i < mimo_blk_length; i++)
	{
		delete[] ep_extr_var_omp[i];
	}
	delete[] ep_extr_var_omp;
	delete[] y_mod_real;
	delete[] x_mod_tmp;
	delete[] y_mod_tmp;
	delete[] y_mod_real_tmp;
	delete[] x_tmp;
	delete[] oriLLR_tmp;
	delete[] sum_sdis_array;
	delete[] Hreal_array;
	delete[] HTH_array;
	delete[] HTY_array;
	delete[] intrlv_pattern_2D_inner;
}





void system2D_LSD_epdet_irregular_SVD_precode(double* BER_array, double* FER_array)
{
	std::ofstream file_total("total_timing.txt", std::ios::app);
	std::ofstream file_det("det_timing.txt", std::ios::app);
	std::ofstream file_dec("dec_timing.txt", std::ios::app);
	std::ofstream file_err_num("total_frames.txt", std::ios::app);
	std::ofstream file_error_pattern("error_patterns.txt", std::ios::app);
	auto det_end = std::chrono::system_clock::now();
	int Nt2 = 2 * Nt;
	int Nr2 = 2 * Nr;
	int cnt_short_codewords = 0;
	for (int i = 0; i < Nv; i++) if (channel_index_predefined[i] == 1) { cnt_short_codewords++; }
	int* bad_channel_index = new int[cnt_short_codewords];
	int* good_channel_index = new int[Nv - cnt_short_codewords];
	int cnt_good = 0;
	int cnt_bad = 0;
	for (int i = 0; i < Nv; i++) { if (channel_index_predefined[i] == 0) { good_channel_index[cnt_good] = i; cnt_good++; } else { bad_channel_index[cnt_bad] = i; cnt_bad++; } }
	/*if (!file.is_open()) {
		std::cerr << "Error opening file: " << file_name << std::endl;
		return;
	}*/

	//file << setw(20) << SNR << setw(20) << BER << setw(20) << FER << "\n";

	//freopen("result.txt", "w", stdout);
	printsetting();
	cout << "\nSNR	 " << "FER      BER" << endl;
	int N_Info = Kv * Kh - CRCL, Nmax;
	if (irregular_coding) { N_Info = border * K_time_array[0] + (Kv - border) * K_time_array[1] - CRCL; }
	//Interleaving pattern
	/*std::vector<int> intrlv_pattern(Nv);
	std::vector<int>* intrlv_pattern_2D = new vector<int>[Nh];
	for (int i = 0; i < Nh; i++) intrlv_pattern_2D[i].resize(Nv);*/
	int intrlv_pattern_num = Nh / ModType;
	std::vector<int> intrlv_pattern(Nv);
	std::vector<int> intrlv_pattern_info_mapping(Nv);
	std::vector<int>* intrlv_pattern_2D = new vector<int>[Nh];
	std::vector<int>* intrlv_pattern_2D_time = new vector<int>[Nv];
	std::vector<int>* intrlv_pattern_2D_inner = new vector<int>[Nh];
	std::vector<int>** intrlv_pattern_2D_time_partial = new vector<int> *[Nv];
	for (int i = 0; i < Nv; i++) intrlv_pattern_2D_time_partial[i] = new vector<int>[intrlv_pattern_num];
	for (int i = 0; i < Nv; i++)
	{
		for (int j = 0; j < intrlv_pattern_num; j++)
		{
			intrlv_pattern_2D_time_partial[i][j].resize(ModType);
		}
	}
	for (int i = 0; i < Nh; i++) intrlv_pattern_2D[i].resize(Nv);
	for (int i = 0; i < Nh; i++) intrlv_pattern_2D_inner[i].resize(Kv);
	for (int i = 0; i < Nv; i++) intrlv_pattern_2D_time[i].resize(Nh);

	float CodeRate = (float)(N_Info-2) / (Nh * Nv);
	Nmax = max(Nh, Nv);
	// num of symbols per vertical codeword
	int Nsym = 0;
	if (MOD_DIM == "Spatial") { Nsym = Nv / ModType; }
	else if (MOD_DIM == "Temporal") { Nsym = Nh / ModType; }
	int symbol_length = ModType;

	// Size of a mimo block is given by height(spatial) * length(temporal)
	int mimo_blk_height = 0;
	int mimo_blk_length = 0;
	if (MOD_DIM == "Spatial")
	{
		mimo_blk_height = Nsym;
		mimo_blk_length = Nh;
	}
	else if (MOD_DIM == "Temporal")
	{
		mimo_blk_height = Nv;
		mimo_blk_length = Nsym;
	}
	// Initialization
	unsigned char* a_data = new unsigned char[N_Info];// Information CodeWord
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
	int** x_parity_check = new int* [Nv];
	for (int i = 0; i < Nh; i++) x_parity_check[i] = new int[Nh];

	// Arrays for MIMO
	int** x_mmse_out = new int* [Nv]; // delete!
	for (int i = 0; i < Nv; i++) { x_mmse_out[i] = new int[Nh]; }

	//std::complex<float>** x_mod = new std::complex<float>* [Nsym]; // modulated symbols (spatial modulation) DELETE!
	//for (int i = 0; i < Nsym; i++) { x_mod[i] = new std::complex<float>[Nh]; }
	//std::complex<float>** x_mod_temporal = new std::complex<float>*[Nv]; // modulated symbols (temporal modulation) DELETE!
	//for (int i = 0; i < Nv; i++) x_mod_temporal[i] = new std::complex<float>[Nsym];
	//std::complex<float>** y_mod = new std::complex<float>* [Nr]; // Output symbols after AWGN channels // delete!
	//for (int i = 0; i < Nr; i++) { y_mod[i] = new std::complex<float>[Nh]; }
	//float** y_mod_real = new float* [2 * Nr];
	//for (int i = 0; i < 2*Nr; i++) y_mod_real[i] = new float[Nh];
	//float** y = new float* [Nv];// Output Signal After AWGN Channel
	//for (int i = 0; i < Nv; i++) { y[i] = new float[Nh]; }
	//float** y_mmse = new float* [2*Nsym];
	//for (int i = 0; i < 2*Nsym; i++) { y_mmse[i] = new float[Nh]; }  // MMSE detector output // delete!
	//float* y_mmse_tmp = new float[2*Nsym];
	//std::complex<float>* x_mod_tmp = new std::complex<float>[Nsym]; // modulated symbols of a single column // delete!
	//std::complex<float>* y_mod_tmp = new std::complex<float>[Nr]; // received symbols of a single column // delete!
	//float* y_mod_real_tmp = new float[2*Nr];
	//// float* y_mod_tmp = new float[Nr];  // delete!


	std::complex<float>** x_mod = new std::complex<float>*[mimo_blk_height]; // modulated symbols (spatial modulation) DELETE!
	for (int i = 0; i < mimo_blk_height; i++) { x_mod[i] = new std::complex<float>[mimo_blk_length]; }
	std::complex<float>** y_mod = new std::complex<float>*[Nr]; // Output symbols after AWGN channels // delete!
	for (int i = 0; i < Nr; i++) { y_mod[i] = new std::complex<float>[mimo_blk_length]; }
	float** y_mod_real = new float* [2 * Nr];
	for (int i = 0; i < 2 * Nr; i++) y_mod_real[i] = new float[mimo_blk_length];
	float** y = new float* [Nv];// Output Signal After AWGN Channel
	for (int i = 0; i < Nv; i++) { y[i] = new float[Nh]; }
	float** y_mmse = new float* [2 * mimo_blk_height];
	for (int i = 0; i < 2 * mimo_blk_height; i++) { y_mmse[i] = new float[mimo_blk_length]; }  // MMSE detector output // delete!
	float* y_mmse_tmp = new float[2 * mimo_blk_height];
	std::complex<float>* x_mod_tmp = new std::complex<float>[mimo_blk_height]; // modulated symbols of a single column // delete!
	std::complex<float>* y_mod_tmp = new std::complex<float>[Nr]; // received symbols of a single column // delete!
	std::complex<float>* x_precoded_tmp = new std::complex<float>[mimo_blk_height]; // precoded symbols of a single column // delete!
	std::complex<float>* y_precoded_tmp = new std::complex<float>[Nr];
	float* y_mod_real_tmp = new float[2 * Nr];
	// float* y_mod_tmp = new float[Nr];  // delete!
	float** y_mod_real_temp_omp = new float* [mimo_blk_length]; // new delete!
	for (int i = 0; i < mimo_blk_length; i++) { y_mod_real_temp_omp[i] = new float[2 * Nr]; }
	float** y_mmse_temp_omp = new float* [mimo_blk_length]; // new delete!
	for (int i = 0; i < mimo_blk_length; i++) { y_mmse_temp_omp[i] = new float[2 * mimo_blk_height]; }
	// ep
	float** ep_extr_mean_omp = new float* [mimo_blk_length];
	for (int i = 0; i < mimo_blk_length; i++) { ep_extr_mean_omp[i] = new float[2 * mimo_blk_height]; }
	float** ep_extr_var_omp = new float* [mimo_blk_length];
	for (int i = 0; i < mimo_blk_length; i++) { ep_extr_var_omp[i] = new float[2 * mimo_blk_height]; }
	int* x_tmp = new int[Nv]; // delete!

	// paths selected by k-best detection


	// CRC bits
	unsigned char* u_crc = new unsigned char[Kv * Kh];

	float** upLLR = new float* [Nv];
	for (int i = 0; i < Nv; i++) { upLLR[i] = new _MM_ALIGN16 float[Nh]; }
	//float** oriLLRv = new float* [Nh];
	//for (int i = 0; i < Nh; i++) { oriLLRv[i] = new _MM_ALIGN16 float[2 * Nv + 2]; }
	float** oriLLRh = new float* [Nv];
	for (int i = 0; i < Nv; i++) { oriLLRh[i] = new _MM_ALIGN16 float[Nh]; }
	int len_ori_llr_tmp = 0;
	if (MOD_DIM == "Spatial") { len_ori_llr_tmp = Nv; }
	else if (MOD_DIM == "Temporal") { len_ori_llr_tmp = ModType * Nt; }
	float* oriLLR_tmp = new _MM_ALIGN16 float[len_ori_llr_tmp];  // delete!
	//float** oriLLR_kbest = new float* [Nh];
	float** oriLLR_kbest = new float* [mimo_blk_length];
	for (int i = 0; i < Nh; i++) { oriLLR_kbest[i] = new float[len_ori_llr_tmp]; }
	double* sum_sdis_array = new double[Nv];
	for (int i = 0; i < Nv; i++) sum_sdis_array[i] = 0;
	int* n_tried_patterns_h = new int[softiter];
	for (int i = 0; i < softiter; i++) n_tried_patterns_h[i] = 0;
	int* n_tried_patterns_v = new int[softiter];
	for (int i = 0; i < softiter; i++) n_tried_patterns_v[i] = 0;

	//*****************************************************************
	// EP detector
	float* Hreal_array = new float[(2 * Nr) * (2 * Nt)]; // delete
	float* HTH_array = new float[(2 * Nt) * (2 * Nt)]; // delete
	float* HTY_array = new float[2 * Nt]; // delete 





	int** GGh = new int* [Nh];
	int** GGv = new int* [Nv];


	MatrixXf Gh = Fn(Nh);  // Generate Matrix of horizontal coding
	MatrixXf Gv = Fn(Nv);  // Generate Matrix of vertical coding
	MatrixXcf H(Nr, Nt); // Complex channel Matrix H: Nr x Nt
	//Matrix<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic> H = MatrixXcf::Random(Nr, Nt); // Real domain channel matrix
	//// //MatrixX Definitions
	MatrixXf H_real(2 * Nr, 2 * Nt); // Real domain channel matrix
	MatrixXf received_signal(2 * Nr, 1);
	MatrixXf HTH(2 * Nt, 2 * Nt);
	MatrixXf Identity = MatrixXf::Identity(2 * Nt, 2 * Nt);
	MatrixXf output_signal(2 * Nt, 1);
	MatrixXf mmse_matrix(2 * Nt, 2 * Nr);
	MatrixXf HTHIinv(2 * Nt, 2 * Nt);

	// for KBest
	MatrixXf R(2 * Nt, 2 * Nt);
	MatrixXf Q(2 * Nr, 2 * Nt);
	MatrixXf z(2 * Nt, 1);



	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> H_real = MatrixXf::Random(2 * Nr, 2 * Nt); // Real domain channel matrix
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> received_signal = MatrixXf::Random(2 * Nr,1);
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> HTH = MatrixXf::Random(2 * Nt, 2 * Nt);
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> Identity = MatrixXf::Random(2 * Nt, 2 * Nt);
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> output_signal = MatrixXf::Random(2 * Nt, 1);
	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> mmse_matrix = MatrixXf::Random(2 * Nt, 2 * Nr);


	// Matrix Definitions
	//Eigen::Matrix<float, 2 * Nr, 2 * Nt> H_real; // ʵ���ŵ�����
	//Eigen::Matrix<float, 2 * Nr, 1> received_signal; // �����źž���
	//Eigen::Matrix<float, 2 * Nt, 2 * Nt> HTH; // HTH����
	//Eigen::Matrix<float, 2 * Nt, 2 * Nt> Identity = Eigen::Matrix<float, 2 * Nt, 2 * Nt>::Identity(); // ��λ����
	//Eigen::Matrix<float, 2 * Nt, 1> output_signal; // ����źž���
	//Eigen::Matrix<float, 2 * Nt, 2 * Nr> mmse_matrix; // MMSE����


	for (int i = 0; i < Nh; i++) GGh[i] = new int[Nh];
	for (int i = 0; i < Nv; i++) GGv[i] = new int[Nv];
	for (int i = 0; i < Nh; i++)
		for (int j = 0; j < Nh; j++)
			GGh[i][j] = Gh(i, j);
	for (int i = 0; i < Nv; i++)
		for (int j = 0; j < Nv; j++)
			GGv[i][j] = Gv(i, j);

	// display GGh
	std::cout << endl;
	/*for (int i = 0; i < Nv; i++)
	{
		int nnz_cnt = 0;
		for (int j = 0; j < Nv; j++)
		{
			cout << GGv[i][j] << " ";
			nnz_cnt += GGv[i][j];	
		}
		cout << nnz_cnt << " " << endl;
	}*/

	//*****************************************************************
	vector<int> temph(Nh);
	vector<int> tempv(Nv);
	vector<int>* construction_time_all = new vector<int>[Nv];
	for (int i = 0; i < Nv; i++) { construction_time_all[i].resize(Nh); }
	int* Ah = new int[Kh]; // Horizontal Info Bits
	int* Ach = new int[Nh - Kh]; // Horizontal Frozen Bits
	int* A_Ach = new int[Nh];
	int* Av = new int[Kv];
	int* Acv = new int[Nv - Kv];
	int* A_Acv = new int[Nv];
	int** A_time_all = new int* [Nv];
	for (int i = 0; i < Nv; i++) { A_time_all[i] = new int[Nh]; for (int j = 0; j < Nh; j++) A_time_all[i][j] = 0; }
	int** Ac_time_all = new int* [Nv];
	for (int i = 0; i < Nv; i++) { Ac_time_all[i] = new int[Nh]; for (int j = 0; j < Nh; j++) Ac_time_all[i][j] = 0; }
	int** A_Ac_time_all = new int* [Nv];
	for (int i = 0; i < Nv; i++) { A_Ac_time_all[i] = new int[Nh]; for (int j = 0; j < Nh; j++) A_Ac_time_all[i][j] = 0; }
	int* K_time_all = new int[Nv];
	for (int i = 0; i < Nv; i++) K_time_all[i] = 0;
	int* tempv_ordered = new int[Nv];
	int** all_error_positions = new int* [Nv];
	for (int i = 0; i < Nv; i++) {
		all_error_positions[i] = new int[Nh];
		for (int j = 0; j < Nh; j++) all_error_positions[i][j] = 0;
	}
	int* use_short_codeword = new int[Nv];
	for (int i = 0; i < Nv; i++) { use_short_codeword[i] = 0; }
	/*int* use_short_codeword = new int[Nv];
	for (int i = 0; i < Nv; i++) { use_short_codeword[i] = 0; }*/

	float** channel_scale_factor = new float* [Nr];
	for (int i = 0; i < Nr; i++) { channel_scale_factor[i] = new float[Nt]; for (int j = 0; j < Nt; j++) { channel_scale_factor[i][j] = 0.0; } }
	float* singular_values_array = new float[Nt];
	for (int i = 0; i < Nt; i++) { singular_values_array[i] = 0.0; }
	float* antenna_power = new float[Nt];
	for (int i = 0; i < Nt; i++) { antenna_power[i] = 1.0; }
	float* antenna_power_single_ant = new float[Nv];
	for (int i = 0; i < Nv; i++) {
		antenna_power_single_ant[i] = 1.0;
	}

	// for test!
	int** u_sum = new int* [Nv];
	for (int i = 0; i < Nv; i++) { u_sum[i] = new int[Nh]; }
	for (int i = 0; i < Nv; i++) { for (int j = 0; j < Nh; j++) { u_sum[i][j] = 0; } }

	// for test
	int** x_encoder_output = new int* [Nv];
	for (int i = 0; i < Nv; i++) { x_encoder_output[i] = new int[Nh]; }
	for (int i = 0; i < Nv; i++) { for (int j = 0; j < Nh; j++) { x_encoder_output[i][j] = 0; } }

	int** x_decoder_output = new int* [Nv];
	for (int i = 0; i < Nv; i++) { x_decoder_output[i] = new int[Nh]; }
	for (int i = 0; i < Nv; i++) { for (int j = 0; j < Nh; j++) { x_decoder_output[i][j] = 0; } }

	int** x_encoder_sum = new int* [Nv];
	for (int i = 0; i < Nv; i++) { x_encoder_sum[i] = new int[Nh]; }
	for (int i = 0; i < Nv; i++) { for (int j = 0; j < Nh; j++) { x_encoder_sum[i][j] = 0; } }

	float** upLLR_pre = new float* [Nv];
	for (int i = 0; i < Nv; i++) { upLLR_pre[i] = new float[Nh]; }
	for (int i = 0; i < Nv; i++) { for (int j = 0; j < Nh; j++) { upLLR_pre[i][j] = 0.0; } }
	//*****************************************************************
	
		// Main Function
	// srand((unsigned int)time(0));
	std::random_device rd;              // ��ȡ���������
	std::mt19937 gen(rd());             // ʹ�����ӳ�ʼ��mt19937������
	std::uniform_int_distribution<> distrib(1, 1000000);
	int Num_Error, Num_Frame_Error;
	float Decoding_Duration;	

	int Num_Error_new, Num_Frame_Error_new;
	float Decoding_Duration_new;

	clock_t T_start, T_finish;
	int SNR_count = 0;
	int anssc, anssd;
	for (double nSNR = Min_SNR; nSNR <= Max_SNR; nSNR += SNR_step)
	{
		float trans_energy = 0.0;
		for (int i = 0; i < Nv; i++) sum_sdis_array[i] = 0;
		double sdis_sum = 0;
		// TEST: det err
		int n_det_err = 0; // n error bits at detector output
		int n_det_err_mrb = 0; // error bits at detector output(MRBs)
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
			// std::cout << sigma_sq << endl;
		}
		// cout << sigma_sq << endl;
		polar_codeconstruction(polarconh, Nh, Kh, GGh, (float)sqrt(sigma_sq), temph);
		polar_codeconstruction(polarconv, Nv, Kv, GGv, (float)sqrt(sigma_sq), tempv);
		// for (int i = 0; i < Nv; i++) { tempv[i] = vertical_construction_fixed[i]; }
		// cout << endl;
		for (int i = 0; i < Nv; i++) { tempv_ordered[i] = tempv[i];/* cout << tempv_ordered[i]+1 << ",";*/ }
		if (SPECIFIC_EXCLUSION)
		{
			for (int i = 0; i < EXCLUSION_NUM; i++)
			{
				int idx = exclusion_set[i];
				for (int j = 0; j < Nv; j++)
				{
					if (tempv_ordered[j] == idx)
					{
						// int tmp = tempv[j];
						tempv_ordered[j] = -1;
						break;
					}
				}
			}
		}
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

		/*for (int i = 0; i < Nh; i++)
		{
			cout << A_Ach[i] << ",";
		}
		cout << endl;*/
		// for (int i = 0; i < Kv; i++) { cout << Av[i] << endl; }
		for (int k = 0; k < Nv; k++)
		{
			bool bfind = false;
			int channel_index = 0;
			for (int i_find = 0; i_find < Kv; i_find++) { if (tempv_ordered[i_find] == k) { bfind = true;  channel_index = i_find; break; } }
			int k_time = K_time_array[1];
			int k_thres = Av[border];
			//cout << k_thres << endl;
			//cout << channel_index << endl;
			// if ((k < k_thres) && A_Acv[k]==1) { k_time = K_time_array[0]; }
			if (SPECIFIC_EXCLUSION) {
				bool bfind_extra = 0;
				for (int i_find = 0; i_find < PROTECTED_PARITY_NUM; i_find++) {
					if (protected_parity_set[i_find] == k) { bfind_extra = 1; break; }
				}
				if (bfind || bfind_extra) { k_time = K_time_array[0]; }
			}
			else {
				if (bfind && ((Kv - 1 - channel_index) < border)) { k_time = K_time_array[0]; }
			}
			K_time_all[k] = k_time;
			// std::cout << sigma_sq << endl;
			polar_codeconstruction(polarconh, Nh, k_time, GGh, (float)sqrt(sigma_sq), construction_time_all[k]);
			if (k_time == K_time_array[0]) { use_short_codeword[k] = 1; } // use long codeword
			else if (k_time == K_time_array[1]) { use_short_codeword[k] = 0; } // use short codeword
			for (int i = 0; i < k_time - 1; i++)
				for (int j = i + 1; j < k_time; j++)
					if (construction_time_all[k][i] > construction_time_all[k][j])
					{
						int ttt = construction_time_all[k][i];
						construction_time_all[k][i] = construction_time_all[k][j];
						construction_time_all[k][j] = ttt;
					}

			// Decide the position of information bits and frozen bits
			for (int i = 0; i < k_time; i++) // Information Set
			{
				A_time_all[k][i] = construction_time_all[k][i];
				A_Ac_time_all[k][A_time_all[k][i]] = 1;
			}
			//	getchar();
			for (int i = 0; i < Nh - k_time; i++) // Frozen Bit
			{
				Ac_time_all[k][i] = construction_time_all[k][k_time + i];
				A_Ac_time_all[k][Ac_time_all[k][i]] = 0;
			}
		}

		// for test
		/*cout << endl;
		cout << "int channel_index_predefined[Nv] = {";
		for (int i = 0; i < Nv; i++) { cout << use_short_codeword[i] << ","; }
		cout << "};" << endl;

		for (int i = 0; i < Nv; i++) {
			for (int j = 0; j < Nh; j++) { cout << A_Ac_time_all[i][j] << " "; }
			cout << K_time_all[i] << " ";
			cout << endl;
		}*/

		// Non-uniform channel
		vector<int> non_uniform_intrlv_pattern;
		vector<int>* non_uniform_intrlv_pattern_array = new vector<int>[Nh];
		non_uniform_intrlv_pattern.resize(Nv);
		for (int i = 0; i < Nh; i++) { non_uniform_intrlv_pattern_array[i].resize(Nv); }
		if (non_uniform_channel)
		{
			getChannelScalingFactor(Nr, Nt, diff_gain, use_short_codeword, channel_scale_factor, predefined_non_uniform);
			if (predefined_non_uniform && (K_time_array[0] != K_time_array[1])) {
				generateIntrlvPatternIrregular(non_uniform_intrlv_pattern, Nv, bad_channel_index,
					good_channel_index, cnt_short_codewords, use_short_codeword);
				if (flag_separate_interleaving) {
					for (int i = 0; i < Nh; i++) {
						generateIntrlvPatternIrregularArray(non_uniform_intrlv_pattern_array[i], Nv, bad_channel_index,
							good_channel_index, cnt_short_codewords, use_short_codeword);
					}
				}
				//// test
				/*cout << endl << "bad channel index: ";
				for (int i = 0; i < cnt_short_codewords; i++) { cout << bad_channel_index[i] << " "; }
				cout << endl << "good channel index: ";
				for (int i = 0; i < Nv-cnt_short_codewords; i++) { cout << good_channel_index[i] << " "; }
				cout << endl << "use short codeword: ";
				for (int i = 0; i < Nv; i++) { cout << use_short_codeword[i] << " "; }
				cout << endl << "interleaving pattern: ";
				cout << "\n BAD" << endl;
				for (int i = 0; i < cnt_short_codewords; i++) { cout << bad_channel_index[i] << "(" << non_uniform_intrlv_pattern_array[0][bad_channel_index[i]] << ")" << "  "; }
				cout << "\n GOOD" << endl;
				for (int i = 0; i < Nv - cnt_short_codewords; i++) { cout << good_channel_index[i] << "(" << non_uniform_intrlv_pattern_array[0][good_channel_index[i]] << ")"<<"  "; }
				cout << endl;*/
			}
			/*if (predefined_non_uniform) { getChannelScalingFactor(Nr, Nt, diff_gain, use_short_codeword, channel_scale_factor); }
			else { getChannelScalingFactor(Nr, Nt, diff_gain, use_short_codeword, channel_scale_factor); }*/
		}
		// display A_Ac_time_all row by row
		//cout << "A_Ac_time_all:" << endl;
		//for (int i = 0; i < Nv; i++)
		//{
		//	for (int j = 0; j < Nh; j++)
		//	{
		//		cout << A_Ac_time_all[i][j] << ",";
		//	}
		//	cout << endl;
		//}
		////// display K_time_all
		//cout << endl;
		//for (int i = 0; i < Nv; i++)
		//{
		//	cout << K_time_all[i] << endl;
		//}
		//test: Horizontal frozen bits

	   /*cout << "Frozen Bits" <<endl;
	   for (int i = 0; i < Nv - Kv; i++)
	   {
		   cout << setw(5) << Acv[i];
	   }
	   cout << endl;*/
	   // test: Horizontal frozen bits
	   // Frame Loop
		int passnum;
		int frame = 0;

		// generate interleaving patterns
		if (flag_interleave)
		{
			/*string ofile_intrlv_name = string("intrlv_pattern_71_frame_49263.txt");
			string ofile_spatial_intrlv_name = string("intrlv_pattern_spatial71.txt");
			ifstream ofile_intrlv(ofile_intrlv_name);
			ifstream ofile_spatial_intrlv(ofile_spatial_intrlv_name);*/
			int file_in = 1;
			//for (int i = 0; i < Nv; i++)
			//{
			//	for (int j = 0; j < intrlv_pattern_num; j++)
			//	{

			//		for (int k = 0; k < ModType; k++)
			//			// ofile_intrlv >> file_in;
			//			ofile_intrlv >> intrlv_pattern_2D_time_partial[i][j][k];
			//		//cout << endl;
			//	}
			//}
			//ofile_intrlv.close();

			//for (int i = 0; i < Nh; i++)
			//{
			//	for (int j = 0; j < Nv; j++)
			//		ofile_spatial_intrlv >> intrlv_pattern_2D[i][j];
			//}
			/*ofile_spatial_intrlv.close();*/
			/*for (int i = 0; i < Nh; i++)
			{
				for (int j = 0; j < Nv; j++)
					cout << intrlv_pattern_2D[i][j] << endl;
			}*/
			generateIntrlvPattern(intrlv_pattern, Nv);
			generateIntrlvPatternInfoMapping(intrlv_pattern_info_mapping, Nv, Kv, Av, Acv);

			// display intrlv_pattern_info_mapping
			// cout << "intrlv_pattern_info_mapping:" << endl;
			/*for (int i = 0; i < Nv; i++)
			{
				cout << i << "  "<< intrlv_pattern_info_mapping[i] << endl;
			}*/
			for (int i = 0; i < Nh; i++) { generateIntrlvPattern(intrlv_pattern_2D[i], Nv); }
			for (int i = 0; i < Nv; i++) { generateIntrlvPattern(intrlv_pattern_2D_time[i], Nh); }
			if (flag_interleave_partial)
			{
				for (int i = 0; i < Nv; i++)
				{
					for (int j = 0; j < intrlv_pattern_num; j++)
					{
						generateIntrlvPattern(intrlv_pattern_2D_time_partial[i][j], ModType);
						/*for (int k = 0; k < ModType; k++)
							cout << intrlv_pattern_2D_time_partial[i][j][k];
						cout << endl;*/
					}
				}
			}
			for (int i = 0; i < Nh; i++) { generateIntrlvPattern(intrlv_pattern_2D_inner[i], Kv); }
		}

		//ofstream ofile_intrlv("intrlv_pattern.txt");
		//for (int i=0; i<Nv; i++)
		//{
		//	for (int j = 0; j < intrlv_pattern_num; j+for+)
		//	{
		//		
		//		for (int k = 0; k < ModType; k++)
		//			ofile_intrlv << intrlv_pattern_2D_time_partial[i][j][k] << "   ";
		//		//cout << endl;
		//	}
		//}
		//ofile_intrlv.close();
		int frame_tmp = 0;
		int frame_delta = 0;
		int n_frame_error_pre = 0;
		while (Num_Frame_Error < errframe)
		{
			/*test: processing time*/
			//auto start = std::chrono::system_clock::now();
			/*if (frame == 10) {
				cout << "AAAAAAAAAA" << endl;
			}*/
			for (int i = 0; i < Nh; i++) { generateIntrlvPattern(intrlv_pattern_2D[i], Nv); }
			for (int i = 0; i < Nv; i++) { generateIntrlvPattern(intrlv_pattern_2D_time[i], Nh); }
			frame++;
			// generate interleaving pattern
			/*if (flag_interleave)
			{
				generateIntrlvPattern(intrlv_pattern, Nv);
				for (int i = 0; i < Nh; i++) { generateIntrlvPattern(intrlv_pattern_2D[i], Nv); }
				for (int i = 0; i < Nv; i++) { generateIntrlvPattern(intrlv_pattern_2D_time[i], Nh); }
				if (flag_interleave_partial)
				{
					for (int i = 0; i < Nv; i++)
					{
						for (int j = 0; j < intrlv_pattern_num; j++)
						{
							generateIntrlvPattern(intrlv_pattern_2D_time_partial[i][j], ModType);
						}
					}
				}
			}*/

			//cout << frame << endl;
			int u_input[16] = { 0,     1     ,0,     1,
									0,     1,     0,     0,
									0 ,    1,     0,     0,
									1 ,    0,     0,     0 };
			for (int i = 0; i < N_Info; i++)
			{
				 a_data[i] = (unsigned char)(distrib(gen) % 2);
				//		a_data[i] = u_input[i];
				// a_data[i] = 0.0;
			}

			// cout << N_Info << endl;
			if (CRCL) tx_append_crc(a_data, N_Info, CRCL);
			for (int i = 0; i < Nv; i++)
				for (int j = 0; j < Nh; j++)
					u[i][j] = 0;
			//omp_set_num_threads(4);
//#pragma omp parallel for
			int info_count = 0;
			for (int i = 0; i < Kv; i++)
			{
				for (int j = 0; j < K_time_all[Av[i]]; j++)
				{
					u[Av[i]][A_time_all[Av[i]][j]] = (int)a_data[info_count];
					info_count++;
					//cout << (int)u[Av[i]][Ah[j]] << " \n";
					//cout << Av[i] << "\n";
				}
				//cout << endl; 
			}
			// test:original info
			//cout << "original info u:" << endl;
			////for (int j = 0; j < Nh; j++) { cout << u[Av[0]][j] << " "; }
			//for (int i = 0; i < Nv; i++)
			//{
			//	cout << endl;
			//	for (int j = 0; j < Nh; j++)
			//		cout << u[i][j];
			//}
			//cout << endl;
			// Polar Encoding
//#pragma omp parallel for

			// test: original info
			/*cout << "original info u:" << endl;
			display_array(u, Nv, Nh);*/
			// test: original info

			for (int i = 0; i < Nv; i++)
				PolarEncode_xor(midx[i], u[i], Nh);
			// test:horizontal encoding
			//cout << "u after horizontal encoding:" << endl;
			////for (int j = 0; j < Nh; j++) { cout << midx[Av[0]][j] << " "; }
			//for (int i = 0; i < Nv; i++)
			//{
			//	cout << endl;
			//	for (int j = 0; j < Nh; j++)
			//		cout << midx[i][j];
			//}
			//cout << endl;
			//cout << endl;
			// test:horizontal encoding
			// 
			// column-wise interleaving

			// test: row encoding
			/*cout << "u after horizontal encoding:" << endl;
			display_array(midx, Nv, Nh);*/
			// test: row encoding
			if (flag_interleave_inner)
			{
				int* tmp_to_interleave = new int[Nv];
				for (int i = 0; i < Nh; i++)
				{
					for (int j = 0; j < Nv; j++) { tmp_to_interleave[j] = midx[j][i]; }
					randIntrlvPartial(intrlv_pattern_2D_inner[i], tmp_to_interleave, Av, Nv, Kv);
					for (int j = 0; j < Nv; j++) { midx[j][i] = tmp_to_interleave[j]; }
				}
				delete[] tmp_to_interleave;
			}

			// test: innner interleaving
			/*cout << "u after inner interleaving:" << endl;
			display_array(midx, Nv, Nh);*/
			// test: innner interleaving
//#pragma omp parallel for
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
			/*cout << "u after 1st spatial encoding" << endl;
			display_array(x, Nv, Nh);*/
			// systematic spatial encoding
			// reset frozen bits to 0
			for (int i = 0; i < Nh; i++)
			{
				for (int j = 0; j < Nv - Kv; j++)
				{
					x[Acv[j]][i] = 0;
				}
			}
			// another round of spatial encoding
			for (int i = 0; i < Nh; i++)
			{
				int* tmpxv = new int[Nv]; int* tmpuv = new int[Nv];
				for (int j = 0; j < Nv; j++)
					tmpuv[j] = x[j][i];
				PolarEncode_xor(tmpxv, tmpuv, Nv);
				for (int j = 0; j < Nv; j++)
					x[j][i] = tmpxv[j];
				delete[] tmpxv; delete[] tmpuv;
			}

			for (int i = 0; i < Nv; i++) {
				int* u_enc_tmp = new int[Nh];
				PolarEncode_xor(u_enc_tmp, x[i], Nh);
				for (int j = 0; j < Nh; j++) {
					u_sum[i][j] += u_enc_tmp[j];
				}
				delete[] u_enc_tmp;
			}


			if (OUTPUT_ERROR_PATTERN) {
				for (int i = 0; i < Nv; i++) {
					for (int j = 0; j < Nh; j++) {
						x_encoder_output[i][j] = x[i][j];
					}
				}
				int* test_enc = new int[Nh];
				for (int i = 0; i < Nv; i++)
				{
					PolarEncode_xor(test_enc, x_encoder_output[i], Nh);
					for (int j = 0; j < Nh; j++) { x_encoder_sum[i][j] += test_enc[j]; }
				}
				delete[] test_enc;
			}
			//if (frame > 20) {
			//	// display u_sum
			//	cout << "u_sum after frame " << frame << endl;
			//	for (int i = 0; i < Nv; i++) {
			//		for (int j = 0; j < Nh; j++) {
			//			cout << (u_sum[i][j]>0) << " ";
			//		}
			//		cout << endl;
			//	}
			//	cout << endl;
			//}
			/*cout << "u after 2nd spatial encoding" << endl;
			display_array(x, Nv, Nh);*/



			// test: vertical encoding
			/*cout << "u after vertical encoding:" << endl;
			display_array(x, Nv, Nh);*/
			// test: vertical encoding
			/*for (int i = 0; i < 2 * Nr; i++)
				for (int j = 0; j < 2 * Nt; j++)
					H_real(i, j) = 0;*/
					/*MatrixXf H_real_test(2 * Nr, 2 * Nt);
					MatrixXf M = H_real_test.transpose();
					MatrixXf N = H_real_test;
					HTHIinv = M * N;*/
					//HTHIinv = M * H_real;
					// test:vertical encoding
					//cout << "u after vertical encoding:" << endl;
					////for (int j = 0; j < Nh; j++) { cout << x[Av[0]][j] << " "; }
					//for (int i = 0; i < Nv; i++)
					//{
					//	cout << endl;
					//	for (int j = 0; j < Nh; j++)
					//		cout << x[i][j];
					//}
					//cout << endl;
					//cout << endl;
					// test:vertical encoding
			if (atten == 2)
			{
				// modulation
				// #pragma omp parallel for
				// Interleaving
				if (flag_interleave)
				{
					if (flag_interleave_all)  // interleave temporally
					{
						for (int i = 0; i < Nv; i++)
						{
							randIntrlv(intrlv_pattern_2D_time[i], x[i], Nh);
						}
					}

					if (flag_interleave_partial)
					{
						for (int i = 0; i < Nv; i++)
						{
							for (int j = 0; j < intrlv_pattern_num; j++)
							{
								randIntrlv(intrlv_pattern_2D_time_partial[i][j], x[i] + j * ModType, ModType);
							}
						}
					}
					if ((!flag_interleave_block) && flag_interleave_spatial)
					{
						for (int j = 0; j < Nh; j++)
						{
							for (int i = 0; i < Nv; i++) { x_tmp[i] = x[i][j]; }
							if (K_time_array[0] == K_time_array[1]) { 
								if (!flag_information_mapping)
								{
									randIntrlv(intrlv_pattern_2D[j], x_tmp, Nv);
								}
								else {
									randIntrlv(intrlv_pattern_info_mapping, x_tmp, Nv);
								}
							}
							else {
								if (!flag_separate_interleaving)
								{
									randIntrlv(non_uniform_intrlv_pattern, x_tmp, Nv);
								}
								else {
									randIntrlv(non_uniform_intrlv_pattern_array[j], x_tmp, Nv);
								}
							}
							for (int i = 0; i < Nv; i++) { x[i][j] = x_tmp[i]; }
						}
					}

					/*cout << "u after interleaving" << endl;
					display_array(x, Nv, Nh);*/

					if (flag_interleave_block)
					{
						randIntrlvBlock(x, Nh, Nv, ModType, interleave_rows);
					}
				}
				if (MOD_DIM == "Spatial")
				{
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
						/*for (int i = 0; i < mimo_blk_height; i++)
						{
							cout << x_mod_tmp[i] << endl;
						}
						cout << "111111" << endl;*/
					}
					// AWGN Channel (received y_mod )
					if (!non_uniform_channel)H = generateChannelMatrix(Nr, Nt);
					else H = generateChannelMatrixNonuniform(Nr, Nt, channel_scale_factor);
					JacobiSVD<MatrixXcf> svd(H, ComputeThinU | ComputeThinV);
					MatrixXcf U = svd.matrixU();
					MatrixXcf V = svd.matrixV();
					auto S = svd.singularValues();

					for (int i = 0; i < Nt; i++) { singular_values_array[i] = S[i]; }

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
						convertYToReal(Nt, y_mod_tmp, y_mod_real_tmp);
						for (int j = 0; j < 2 * Nr; j++) { y_mod_real[j][i] = y_mod_real_tmp[j]; /*cout << y_mod_real_tmp[j] << endl;*/ }
					}
					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < 2 * Nr; j++) { y_mod_real_temp_omp[i][j] = y_mod_real[j][i]; }
					}
					// for KBest Detection
					
					//omp_set_num_threads(8);

#pragma omp parallel for private(z)
					for (int j = 0; j < Nh; j++)
					{
						MatrixXf received_signal_omp(2 * Nr, 1);
						int k_kbest_extend = pow(2, ModType / 2) * K_KBEST;
						float** paths_kbest = new float* [k_kbest_extend];
						for (int i = 0; i < k_kbest_extend; i++) paths_kbest[i] = new float[2 * Nt];
						float* ED = new float[K_KBEST];

						if (MIMODET == "MMSE") { MMSEdetection(H_real, Identity, received_signal, mmse_matrix, output_signal, HTH, y_mod_real_temp_omp[j], y_mmse_temp_omp[j], sigma_sq); }
						// ִ��һЩ����

						
						if (MIMODET == "KBEST")

						{
							Eigen::HouseholderQR<Eigen::MatrixXf> qr_Hreal(H_real);
							Q = qr_Hreal.householderQ();
							R = qr_Hreal.matrixQR().triangularView<Eigen::Upper>();
							for (int i = 0; i < 2 * Nr; i++) { received_signal_omp(i, 0) = y_mod_real_temp_omp[j][i]; }
							z = Q.transpose() * received_signal_omp;
							KBestDetection(R, z, H_real, y_mod_real_temp_omp[j], ED, y_mmse_temp_omp[j], paths_kbest, K_KBEST, ModType, Nt, Nr);
							llr_kbest_bit(Nt, y_mmse_temp_omp[j], oriLLR_kbest[j], sigma_sq, ModType, K_KBEST, paths_kbest, ED);
						}
						for (int i_path = 0; i_path < k_kbest_extend; i_path++) { delete[] paths_kbest[i_path]; }

						delete[] paths_kbest;
						delete[] ED;
					}

					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < 2 * Nt; j++) { y_mmse[j][i] = y_mmse_temp_omp[i][j]; }
					}
					// sym llr to bitwise llr
					// Detection is done for every column, that is, for every Nv bits
					// Interleaving is not currently considered
					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < 2 * Nt; j++) { y_mmse_tmp[j] = y_mmse[j][i]; }
						if (MIMODET == "MMSE") { llr_mmse_bit(Nt, y_mmse_tmp, oriLLR_tmp, sigma_sq, ModType); }
						if (MIMODET == "KBEST") { for (int j = 0; j < len_ori_llr_tmp; j++) { oriLLR_tmp[j] = oriLLR_kbest[i][j]; } }
						if (flag_interleave)
						{
							randDeintrlv(intrlv_pattern, oriLLR_tmp, Nv);
						}
						bool flagIsMax = false;
						for (int j = 0; j < Nv; j++)
						{
							if (!(j % 2))
							{
								if (abs(x_mmse_out[j][i]) > abs(x_mmse_out[j][i + 1])) { flagIsMax = true; }
								else { flagIsMax = false; }
							}
							else
							{
								if (abs(x_mmse_out[j][i]) > abs(x_mmse_out[j][i - 1])) { flagIsMax = true; }
								else { flagIsMax = false; }
							}
							oriLLRh[j][i] = oriLLR_tmp[j];
							upLLR[j][i] = oriLLR_tmp[j];
							// TEST�� det error
							x_mmse_out[j][i] = 0;
							if (oriLLR_tmp[j] < 0) { x_mmse_out[j][i] = 1; }
							if (x_mmse_out[j][i] != x[j][i]) {
								n_det_err++;
								// error at max-reliable bit. Warning: Only fit for 16QAM modulation
								if (flagIsMax) { n_det_err_mrb++; }
							}
							total_bit_cnt++;
						}
					}
					// test : detection TIME
					/*det_end = std::chrono::system_clock::now();
					auto duration_detection = std::chrono::duration_cast<std::chrono::microseconds>(det_end - start);
					file_det << duration_detection.count() << endl;*/
				}
				// Temporal modulation
				else if (MOD_DIM == "Temporal")
				{
					
					//modulation
					//cout << "START TEST" << endl;
					for (int i = 0; i < Nv; i++) // every row
					{
						modulation(x[i], x_mod[i], ModType, mimo_blk_length);
					}
					/*for (int i = 0; i < mimo_blk_height; i++)
					{
						for (int j = 0; j < mimo_blk_length; j++)
							cout << x_mod[i][j] << endl;
					}
					cout << "111111" << endl;*/
					//AWGN channel (received y_mod)
					// test: modulated symbols
					/*cout << endl;
					cout << "modulated symbols" << endl;
					for (int i = 0; i < Nt; i++) cout << x_mod[i][0] << endl;
					cout << "modulated symbols" << endl;*/
					// test: modulated symbols
					// H = generateChannelMatrix(Nr, Nt);
					if (!non_uniform_channel)H = generateChannelMatrix(Nr, Nt);
					else H = generateChannelMatrixNonuniform(Nr, Nt, channel_scale_factor);
					JacobiSVD<MatrixXcf> svd(H, ComputeFullU | ComputeFullV);
					MatrixXcf U = svd.matrixU();
					// cout << U.rows() << " " << U.cols() << endl;
					MatrixXcf V = svd.matrixV();
					// cout << V.rows() << " " << V.cols() << endl;
					auto S = svd.singularValues();
					for (int i = 0; i < Nt; i++) { singular_values_array[i] = S[i]; }
					// if (nSNR == 10.0) { cout << "Temporal" << endl; }
					if (mimo_svd_average_power_alloc)
					{
						float my_sum = 0.0;
						if (K_time_array[0] == K_time_array[1]) {
							for (int i = 0; i < Nt; i++) { singular_values_array[i] = S[i]; }
							float sum_tmp = 0;
							for (int i = 0; i < Nt; i++) {
								sum_tmp += 1.0 / (singular_values_array[i] * singular_values_array[i]);
							}
							for (int i = 0; i < Nt; i++) {
								float sing_tmp = 1.0 / (singular_values_array[i] * singular_values_array[i]);
								antenna_power[i] = (float)Nt * sing_tmp / sum_tmp;
								my_sum += antenna_power[i];
								/*cout << antenna_power[i] << "  ";
								cout << sqrt(antenna_power[i]) * singular_values_array[i] << "  ";
								cout << endl;*/
							}
							// cout << "My_sum = " << my_sum << endl;
						}
						else {
							for (int i = 0; i < Nt; i++) { singular_values_array[i] = S[i]; }
							float sum_tmp = 0;
							int n_tmp = border;
							for (int i = 0; i < n_tmp; i++) {
								sum_tmp += 1.0 / (singular_values_array[i] * singular_values_array[i]);
							}
							for (int i = 0; i < n_tmp; i++) {
								float sing_tmp = 1.0 / (singular_values_array[i] * singular_values_array[i]);
								antenna_power[i] = (float)n_tmp * sing_tmp / sum_tmp;
								/*cout << antenna_power[i] << "  ";
								cout << sqrt(antenna_power[i]) * singular_values_array[i] << "  ";
								cout << endl;*/
							}
							sum_tmp = 0;
							n_tmp = Nt - border;
							for (int i = 0; i < n_tmp; i++) {
								sum_tmp += 1.0 / (singular_values_array[i+border] * singular_values_array[i+border]);
							}
							for (int i = 0; i < n_tmp; i++) {
								float sing_tmp = 1.0 / (singular_values_array[i+border] * singular_values_array[i+border]);
								antenna_power[i+border] = (float)n_tmp * sing_tmp / sum_tmp;
								/*cout << antenna_power[i] << "  ";
								cout << sqrt(antenna_power[i]) * singular_values_array[i] << "  ";
								cout << endl;*/
							}

						}
						for (int i = 0; i < Nt; i++) { S[i] = S[i] * sqrt(antenna_power[i]); }
					}
					else if (mimo_waterfilling_power_alloc) {
						for (int i = 0; i < Nt; i++) { singular_values_array[i] = S[i]; }
						float sum_tmp = 0;
						for (int i = 0; i < Nt; i++) {
							sum_tmp += (singular_values_array[i]);
						}
						for (int i = 0; i < Nt; i++) {
							float sing_tmp = singular_values_array[i];
							antenna_power[i] = (float)Nt * sing_tmp / sum_tmp;
							/*cout << antenna_power[i] << "  ";
							cout << sqrt(antenna_power[i]) * singular_values_array[i] << "  ";
							cout << endl;*/
						}
						for (int i = 0; i < Nt; i++) { S[i] = S[i] * sqrt(antenna_power[i]); }
					}
					else if (mimo_rate_matching_power_alloc) {
						if (!SPECIFIC_EXCLUSION) {
							for (int i = 0; i < Nt; i++) { singular_values_array[i] = S[i]; }
							float sum_tmp = 0;
							/*for (int i = 0; i < Nt; i++) {
								sum_tmp += (singular_values_array[i]);
							}*/
							sum_tmp = border * mimo_rate_matching_enlarge_ratio + (Nv - border) * 1.0;
							for (int i = 0; i < Nt; i++) {
								float power_tmp;
								if (i > Nt - border - 1) { power_tmp = mimo_rate_matching_enlarge_ratio; }
								else { power_tmp = 1.0; }
								// float sing_tmp = singular_values_array[i];
								antenna_power[i] = (float)Nt * power_tmp / sum_tmp;
								/*cout << antenna_power[i] << "  ";
								cout << sqrt(antenna_power[i]) * singular_values_array[i] << "  ";
								cout << endl;*/
							}
							tmp = 0;
							/*for (int i = 0; i < Nt; i++) { tmp += antenna_power[i]; }
							cout << tmp << endl;*/
							for (int i = 0; i < Nt; i++) { S[i] = S[i] * sqrt(antenna_power[i]); }
						}
						else {
							for (int i = 0; i < Nt; i++) { singular_values_array[i] = S[i]; }
							float sum_tmp = 0;
							/*for (int i = 0; i < Nt; i++) {
								sum_tmp += (singular_values_array[i]);
							}*/
							int border_exclude = 15;//border + PROTECTED_PARITY_NUM;
							sum_tmp = border_exclude * mimo_rate_matching_enlarge_ratio + (Nv - border_exclude) * 1.0;
							for (int i = 0; i < Nt; i++) {
								float power_tmp;
								if (i > Nt - border_exclude - 1) { power_tmp = mimo_rate_matching_enlarge_ratio; }
								else { power_tmp = 1.0; }
								// float sing_tmp = singular_values_array[i];
								antenna_power[i] = (float)Nt * power_tmp / sum_tmp;
								/*cout << antenna_power[i] << "  ";
								cout << sqrt(antenna_power[i]) * singular_values_array[i] << "  ";
								cout << endl;*/
							}
							tmp = 0;
							/*for (int i = 0; i < Nt; i++) { tmp += antenna_power[i]; }
							cout << tmp << endl;*/
							for (int i = 0; i < Nt; i++) { S[i] = S[i] * sqrt(antenna_power[i]); }
						}
					}
					else if (MIMO_SVD_SEPARATE_AVERAGE_POWER_ALLOC) {
						float my_sum = 0.0;
						for (int i = 0; i < Nt; i++) { singular_values_array[i] = S[i]; }
						float sum_tmp = 0;
						// int border = 15;
						for (int i = 0; i < Nt; i++) {
							if (i > Nt - border -1) {
								sum_tmp += 1.0 / (singular_values_array[i] * singular_values_array[i]) * MIMO_SVD_SEPARATE_POWER_RATIO;
							}
							else {
								sum_tmp += 1.0 / (singular_values_array[i] * singular_values_array[i]);
							}
						}
						for (int i = 0; i < Nt; i++) {
							float sing_tmp = 1.0 / (singular_values_array[i] * singular_values_array[i]);
							if (i > Nt - border - 1) {
								antenna_power[i] = (float)Nt * sing_tmp / sum_tmp * MIMO_SVD_SEPARATE_POWER_RATIO;
							}
							else {
								antenna_power[i] = (float)Nt * sing_tmp / sum_tmp;
							}
							my_sum += antenna_power[i];
							/*cout << antenna_power[i] << "  ";
							cout << sqrt(antenna_power[i]) * singular_values_array[i] << "  ";
							cout << endl;*/
						}
						for (int i = 0; i < Nt; i++) { S[i] = S[i] * sqrt(antenna_power[i]); }
						// cout << my_sum << endl;
					}
// # pragma omp parallel for
					for (int j = 0; j < Nsym; j++)
					{
						/*std::complex<float>* x_mod_tmp = new std::complex<float>[Nv];
						std::complex<float>* x_precoded_tmp = new std::complex<float>[Nt];
						std::complex<float>* y_mod_tmp = new std::complex<float>[Nr];
						std::complex<float>* y_precoded_tmp = new std::complex<float>[Nr];*/
						//// Extracting a single column
						for (int i = 0; i < Nv; i++) { x_mod_tmp[i] = x_mod[i][j] * sqrt(antenna_power[i]); 
						trans_energy += (x_mod_tmp[i].real()* x_mod_tmp[i].real()+ x_mod_tmp[i].imag()* x_mod_tmp[i].imag());
						}
						precodeSVD(Nr, Nt, x_mod_tmp, x_precoded_tmp, V); // precoding: x_precode = V*x
						// for (int i = 0; i < Nt; i++) { x_precoded_tmp[i] = x_precoded_tmp[i] * sqrt(antenna_power[i]); }
						AWGN(Nr, Nt, x_precoded_tmp, y_mod_tmp, H, sigma_sq);
						precodeSVD_preprocessing(Nr, Nt, y_mod_tmp, y_precoded_tmp, U); // precoding: y_precode = U^H*y
						for (int i = 0; i < Nr; i++) { y_mod[i][j] = y_precoded_tmp[i]; }
						for (int i = 0; i < Nt; i++) { y_mod[i][j] = S[i] * y_mod[i][j] / (S[i] * S[i]  +sigma_sq); } // MMSE detection with precoding
						//for (int i = 0; i < Nt; i++) { y_mod[i][j] = S[i] * UTy(i,j) / (S[i] * S[i] + sigma_sq); } // MMSE detection with precoding
						//delete[] x_mod_tmp; delete[] x_precoded_tmp; delete[] y_mod_tmp; delete[] y_precoded_tmp;
					}
					// for (int j = 0; j < Nt; j++) { std::cout << y_mod[j][0]*S[j](S[j]*S[j]+sigma_sq) << " "; } std::cout << std::endl;
					// for (int j = 0; j < Nt; j++) { std::cout << UTy(j,0) / S[j] << " "; } std::cout << std::endl;
					// for (int j = 0; j < Nt; j++) { std::cout << y_mod[j][0] << endl; }
					// convertHToReal(H, H_real);
					

					for (int j = 0; j < Nsym; j++)
					{
						for (int i = 0; i < Nr; i++) { y_mod_tmp[i] = y_mod[i][j]; }
						convertYToReal(Nt, y_mod_tmp, y_mod_real_tmp);
						for (int i = 0; i < 2 * Nr; i++) { y_mod_real[i][j] = y_mod_real_tmp[i]; }
					}
					//MMSE detection done in paralell
					for (int i = 0; i < Nsym; i++)
					{
						for (int j = 0; j < 2 * Nr; j++) { y_mod_real_temp_omp[i][j] = y_mod_real[j][i]; }
					}
					//Eigen::HouseholderQR<Eigen::MatrixXf> qr_Hreal(H_real);
					//Q = qr_Hreal.householderQ();
					//R = qr_Hreal.matrixQR().triangularView<Eigen::Upper>();
					//// for EP detection
					//// copy H_real
					//for (int i = 0; i < Nr2; i++)
					//{
					//	for (int j = 0; j < Nt2; j++)
					//	{
					//		Hreal_array[i * Nt2 + j] = H_real(i, j);
					//	}
					//}
					// calculate HTH
					// cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, Nt2, Nt2, Nr2, 1.0, Hreal_array, Nt2, Hreal_array, Nt2, 0.0, HTH_array, Nt2);
					//cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, Nt2, 1, Nr2, 1.0, ch_mtx_struct.H_real[cnt], Nt2, data_struct.rx_buffer_real[cnt], 1, 0.0, ch_mtx_struct.hty, 1);

					//omp_set_num_threads(8);
					//cout << 1 << endl;
					//omp_set_num_threads(4);
					/*auto start = std::chrono::system_clock::now();*/
					/*MatrixXf H_real_trans(2*Nt, 2*Nr);
					for (int i = 0; i < 2 * Nt; i++)
						for (int j = 0; j < 2 * Nr; j++)
							H_real_trans(i, j) = H_real(j, i);
					HTHIinv = (H_real_trans * H_real + Identity).inverse();*/
//# pragma omp parallel for private(z)
//					for (int j = 0; j < Nsym; j++)
//
//						//for (int x=0; x<Nsym; x++)
//					{
//						// int j = x % Nsym;
//						MatrixXf received_signal_omp(2 * Nr, 1);
//						int k_kbest_extend = pow(2, ModType / 2) * K_KBEST;
//						float** paths_kbest = new float* [k_kbest_extend];
//						for (int i = 0; i < k_kbest_extend; i++) { paths_kbest[i] = new float[2 * Nt]; }
//						float* ED = new float[K_KBEST];
//						float* HTY_array_tmp = new float[2 * Nt];
//						// Declared arrays for EP
//						int flag = 1;
//						int temp = 0;
//						//float* extr_mean = new float[2 * Nr];
//						//float* extr_var = new float[2 * Nr];
//						double* vec_alpha_in = new double[2 * Nr];
//						double* vec_gamma_in = new double[2 * Nr];
//						double* p_symbol_rearrange = new double[Nt2 * pow(2, ModType / 2)];
//						for (int i = 0; i < Nt2 * pow(2, ModType / 2); i++) { p_symbol_rearrange[i] = 1 / (pow(2, ModType / 2)); }
//						for (int i = 0; i < Nr2; i++)
//						{
//							//extr_mean[i] = 0;
//							//extr_var[i] = 0;
//							vec_alpha_in[i] = 2;
//							vec_gamma_in[i] = 0;
//						}
//
//						if (MIMODET == "MMSE")
//						{
//							MMSEdetection(H_real, Identity, received_signal, mmse_matrix, output_signal, HTH, y_mod_real_temp_omp[j], y_mmse_temp_omp[j], sigma_sq);
//							//MMSEdetection_noInverse(H_real, Identity, received_signal, mmse_matrix, output_signal, HTH, HTHIinv, y_mod_real_temp_omp[j], y_mmse_temp_omp[j], sigma_sq);
//						}
//						if (MIMODET == "KBEST")
//						{
//							for (int i = 0; i < 2 * Nr; i++) { received_signal_omp(i, 0) = y_mod_real_temp_omp[j][i]; }
//							z = Q.transpose() * received_signal_omp;
//							KBestDetection(R, z, H_real, y_mod_real_temp_omp[j], ED, y_mmse_temp_omp[j], paths_kbest, K_KBEST, ModType, Nt, Nr);
//							llr_kbest_bit(Nt, y_mmse_temp_omp[j], oriLLR_kbest[j], sigma_sq, ModType, K_KBEST, paths_kbest, ED);
//						}
//						if (MIMODET == "EP")
//						{
//							// test: save H_real, recv_signal
//							/*save_mat_txt(y_mod_real_temp_omp[j], Nr2, "RxSig0128.txt");
//							save_mat_txt(H_real, "Hreal0128.txt");*/
//							// test: save H_real, recv_signal
//							cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, Nt2, 1, Nr2, 1.0, Hreal_array, Nt2, y_mod_real_temp_omp[j], 1, 0.0, HTY_array_tmp, 1);
//							//test:HTH_array
//							/*for (int i = 0; i < Nt2; i++)
//								for (int s = 0; s < Nt2; s++)
//								{
//									if (s == 0) std::cout << endl;
//									std::cout << setw(10) << HTH_array[i * Nt2 + s];
//								}
//							cout << endl;*/
//							// TEST:HTY
//							/*for (int i = 0; i < Nr2; i++)
//								received_signal_omp(i, 0) = y_mod_real_temp_omp[j][i];
//							MatrixXf HTransY = H_real.transpose() * received_signal_omp;
//							cout << "Eigen" << endl;
//							for (int i = 0; i < Nt2; i++) cout << setw(10)<<HTransY(i, 0);
//							cout << "mkl" << endl;
//							for (int i = 0; i < Nt2; i++) cout << HTY_array[i];*/
//							// TEST:HTY
//							EPD_detection_float(flag, Nt, Nr, iter_MIMO, ModType, 1.0, sigma_sq, HTH_array, HTY_array_tmp, temp, ep_extr_mean_omp[j], ep_extr_var_omp[j], vec_alpha_in, vec_gamma_in, p_symbol_rearrange);
//						}
//						//test EP detection
//						/*cout << "TEST:EP DETECTION" << endl;
//						for (int i = 0; i < Nr2; i++) std::cout << extr_mean[i]<< endl;
//						cout << "TEST:EP DETECTION" << endl;*/
//						//test EP detection
//						for (int i_path = 0; i_path < k_kbest_extend; i_path++) { delete[] paths_kbest[i_path]; }
//						delete[] paths_kbest;
//						delete[] ED;
//						//delete[] extr_mean;
//						//delete[] extr_var;
//						delete[] vec_alpha_in;
//						delete[] vec_gamma_in;
//						delete[] p_symbol_rearrange;
//						delete[] HTY_array_tmp;
//					}
					// test: detection time
					/*auto detection_over = std::chrono::system_clock::now();
					auto duration_detection = std::chrono::duration_cast<std::chrono::microseconds>(detection_over - start);
					cout << endl;
					cout <<"Detection duration" << duration_detection.count() << "us" << endl;*/
					// test: detection time
					for (int i = 0; i < Nsym; i++)
					{
						for (int j = 0; j < 2 * Nt; j++) { y_mmse[j][i] = y_mmse_temp_omp[i][j]; }
					}
					//To bit domain output
					for (int j = 0; j < Nsym; j++)
					{
						for (int i = 0; i < 2 * Nt; i++) { y_mmse_tmp[i] = y_mod_real[i][j]; }
						// With Spatial Modulation, Nt = Nv
						if (MIMODET == "MMSE") { llr_mmse_bit(Nt, y_mmse_tmp, oriLLR_tmp, sigma_sq, ModType); }
						if (MIMODET == "EP") { llr_ep_bit(Nt, ep_extr_mean_omp[j], ep_extr_var_omp[j], oriLLR_tmp, ModType); }
						if (MIMODET == "KBEST") { for (int i = 0; i < len_ori_llr_tmp; i++) { oriLLR_tmp[i] = oriLLR_kbest[j][i]; } }
						int i_bit = 0;
						int bit_col = 0;
						for (int j_inner = 0; j_inner < Nv; j_inner++)
						{
							for (int i_bit_col = 0; i_bit_col < ModType; i_bit_col++)
							{
								i_bit = ModType * j_inner + i_bit_col;
								bit_col = j * ModType + i_bit_col;
								oriLLRh[j_inner][bit_col] = oriLLR_tmp[i_bit];
								upLLR[j_inner][bit_col] = oriLLR_tmp[i_bit];
							}

						}
					}
					// test: LLR transfer time
					/*auto llr_to_bit_over = std::chrono::system_clock::now();
					auto duration_llr = std::chrono::duration_cast<std::chrono::microseconds>(llr_to_bit_over - detection_over);
					cout << endl;
					cout<< "LLR transfer duration" << duration_llr.count() << "us" << endl;*/
					// test: LLR transfer time
					if (flag_interleave)
					{
						if ((!flag_interleave_block) && flag_interleave_spatial)
						{
							for (int j = 0; j < Nh; j++)
							{
								for (int i = 0; i < Nv; i++) { oriLLR_tmp[i] = oriLLRh[i][j]; }
								if (K_time_array[0] == K_time_array[1]) {
									if (!flag_information_mapping)  randDeintrlv(intrlv_pattern_2D[j], oriLLR_tmp, Nv);
									else { 
										randDeintrlv(intrlv_pattern_info_mapping, oriLLR_tmp, Nv); 
									}
								}
								else { 
									if (!flag_separate_interleaving) { randDeintrlv(non_uniform_intrlv_pattern, oriLLR_tmp, Nv); }
									else { randDeintrlv(non_uniform_intrlv_pattern_array[j], oriLLR_tmp, Nv); }
								}
								for (int i = 0; i < Nv; i++)
								{
									oriLLRh[i][j] = oriLLR_tmp[i];
									upLLR[i][j] = oriLLR_tmp[i];
								}

							}
						}

						if (flag_interleave_block)
						{
							randDeIntrlvBlock(oriLLRh, Nh, Nv, ModType, interleave_rows);
							randDeIntrlvBlock(upLLR, Nh, Nv, ModType, interleave_rows);
						}

						if (flag_interleave_all) // Deinterleave in time domain
						{
							for (int j = 0; j < Nv; j++)
							{
								randDeintrlv(intrlv_pattern_2D_time[j], oriLLRh[j], Nh);
								randDeintrlv(intrlv_pattern_2D_time[j], upLLR[j], Nh);
							}
						}

						if (flag_interleave_partial)
						{
							for (int j = 0; j < Nv; j++)
							{
								for (int k = 0; k < intrlv_pattern_num; k++)
								{
									randDeintrlv(intrlv_pattern_2D_time_partial[j][k], oriLLRh[j] + k * ModType, ModType);
									randDeintrlv(intrlv_pattern_2D_time_partial[j][k], upLLR[j] + k * ModType, ModType);
								}
							}
						}

					}

					for (int i = 0; i < Nv; i++)
					{
						for (int j = 0; j < Nh; j++)
						{
							total_bit_cnt++;
							int bit_dec = (oriLLRh[i][j] < 0);
							if (bit_dec != x[i][j]) n_det_err++;
						}
					}

				}
				// test: detection errors

				// test: detection errors
			}
			if (atten == 1)
			{
				if (flag_interleave)
				{
					if (flag_interleave_all)  // interleave temporally
					{
						for (int i = 0; i < Nv; i++)
						{
							randIntrlv(intrlv_pattern_2D_time[i], x[i], Nh);
						}
					}

					if (flag_interleave_partial)
					{
						for (int i = 0; i < Nv; i++)
						{
							for (int j = 0; j < intrlv_pattern_num; j++)
							{
								randIntrlv(intrlv_pattern_2D_time_partial[i][j], x[i] + j * ModType, ModType);
							}
						}
					}
					if ((!flag_interleave_block) && flag_interleave_spatial)
					{
						for (int j = 0; j < Nh; j++)
						{
							for (int i = 0; i < Nv; i++) { x_tmp[i] = x[i][j]; }
							if (true) {
								if (!flag_information_mapping)
								{
									randIntrlv(intrlv_pattern_2D[j], x_tmp, Nv);
								}
								else {
									randIntrlv(intrlv_pattern_info_mapping, x_tmp, Nv);
								}
							}
							else {
								if (!flag_separate_interleaving)
								{
									randIntrlv(non_uniform_intrlv_pattern, x_tmp, Nv);
								}
								else {
									randIntrlv(non_uniform_intrlv_pattern_array[j], x_tmp, Nv);
								}
							}
							for (int i = 0; i < Nv; i++) { x[i][j] = x_tmp[i]; }
						}
					}

					/*cout << "u after interleaving" << endl;
					display_array(x, Nv, Nh);*/

					if (flag_interleave_block)
					{
						randIntrlvBlock(x, Nh, Nv, ModType, interleave_rows);
					}
				}

				if (MIMO_SVD_SEPARATE_AVERAGE_POWER_ALLOC && (K_time_array[0]!=K_time_array[1])) {
					float low_rate_power = Nv * (MIMO_SVD_SEPARATE_POWER_RATIO) / (MIMO_SVD_SEPARATE_POWER_RATIO * border + 1.0 * (Nv - border));
					float high_rate_power = Nv * 1.0/ (MIMO_SVD_SEPARATE_POWER_RATIO * border + 1.0 * (Nv - border));
					for (int i = 0; i < Nv; i++) {
						if (K_time_all[i] == K_time_array[0]) { antenna_power_single_ant[i] = low_rate_power; }
						else { antenna_power_single_ant[i] = high_rate_power; }
					}
				}
				//for (int i = 0; i < Nv; i++) { std::cout << antenna_power_single_ant[i] << endl; }
				for (int i = 0; i < Nv; i++)
					for (int j = 0; j < Nh; j++)
					{
						// 0-1; 1-(-1)
						y[i][j] = sqrt(antenna_power_single_ant[i])*(1 - 2 * x[i][j]) + (float)sqrt(sigma_sq) * (float)sampleNormal();
						oriLLRh[i][j] = 2 * sqrt(antenna_power_single_ant[i]) * y[i][j] / sigma_sq;
						//oriLLRv[j][i] = oriLLRh[i][j];
						upLLR[i][j] = oriLLRh[i][j];
					}

				float* oriLLR_tmp_v = new float[Nv];
				float* oriLLR_tmp_v_de = new float[Nv];
				float* oriLLR_tmp_h = new float[Nh];
				if (flag_interleave) {
					if (flag_interleave_spatial) {
						for (int i = 0; i < Nh; i++) {
							for (int j = 0; j < Nv; j++) { oriLLR_tmp_v[j] = oriLLRh[j][i]; }
							randDeintrlv(intrlv_pattern_2D[i], oriLLR_tmp_v, Nv);
							for (int j = 0; j < Nv; j++) { oriLLRh[j][i] = oriLLR_tmp_v[j]; }
						}
					}
					if (flag_interleave_all) {
						for (int i = 0; i < Nv; i++) { randDeintrlv(intrlv_pattern_2D_time[i], oriLLRh[i], Nh); }
					}
				}
				delete[] oriLLR_tmp_v;
				delete[] oriLLR_tmp_h;
				delete[] oriLLR_tmp_v_de;
				for (int i = 0; i < Nv; i++) {
					for (int j = 0; j < Nh; j++) { upLLR[i][j] = oriLLRh[i][j]; }
				}
			}
			//Polar Decoding
			for (int iter = 1; iter <= softiter; iter++)
			{
				//vertical decoding
		   // test: decoding time
			//auto decoding_start = std::chrono::system_clock::now();

			// test: decoding time
#pragma omp parallel for
				for (int decodingi = 0; decodingi < Nh; decodingi++) // 291.70kb
				{

					//for SCL Decoding
					int Counth = 0, Countinfoh = 0;
					float Qsumunh = 0;
					int** sum = new int* [L];
					//for (int i = 0; i < L; i++)  sum[i] = new  int[Nmax];
					int** u1 = new int* [L];
					for (int i = 0; i < L; i++)  u1[i] = new _MM_ALIGN16 int[Nmax];
					for (int i = 0; i < L; i++)  sum[i] = new _MM_ALIGN16 int[Nmax];
					//for (int i = 0; i < L; i++)  u1[i] = new  int[Nmax];

					float* LLR_in = new float[L];
					float* W = new float[2 * L];
					int* SCL_index = new int[L];
					int* better = new int[L];
					int* worse = new int[L];
					int* indexpm = new int[L];
					float* PM = new float[L];
					float** LLR = new float* [L];
					for (int i = 0; i < L; i++)  LLR[i] = new _MM_ALIGN16 float[2 * Nmax + 2];
					// for (int i = 0; i < L; i++)  LLR[i] = new float[2 * Nmax + 2];
					////////////////////////////////
					for (int j = 0; j < Nv; j++)
						for (int k = 0; k < L; k++)
							LLR[k][j] = upLLR[j][decodingi];
					//cout << LLR[0][j] << " ";

				//cout << endl;
					int index = 0;
					/*			if (L == 1)
								{
									//T_start = clock();
									Count = 0;
									decode(*LLR, *sum, (int)(log(Nv) / log(2)), A_Acv, flag, Nv, Kv);
									//PolarEncode_xor(*u1, *sum, N);
									//T_finish = clock();
								}
								else
								{*/
					for (int k = 0; k < L; k++) { PM[k] = 0; indexpm[k] = k; }
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
					if (adaptive_beta)
					{
						if (softmethod == "Pyndiah") Pyndiahv_adaptivebeta(LLR, sum, oriLLRh, indexpm, iter * 2 - 2, decodingi, upLLR);
						if (softmethod == "MITSO") MITSOv(LLR, PM, sum, oriLLRh, indexpm, iter * 2 - 2, decodingi, upLLR, Qsumunh);
					}
					// iter * 2 - 2
					else
					{
						//cout << iter << endl;
						if (softmethod == "Pyndiah") Pyndiahv(LLR, sum, oriLLRh, indexpm, iter, decodingi, upLLR);
						if (softmethod == "MITSO") MITSOv(LLR, PM, sum, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunh);
					}
					// test: oriLLR and upLLR
					/*cout << endl;
					for (int i = 0; i < Nv; i++)
						for (int j = 0; j < Nh; j++)
							cout << oriLLRh[i][j] - upLLR[i][j] << "   ";
					cout << endl;*/
					/*cout << oriLLRh[0][0] << endl;
					cout << upLLR[0][0] << endl;*/
					// test: oriLLR and upLLR
					// half iteration Enabled
					// May Need Modification
					if (half_iter && iter == softiter)
					{
						if (iter == softiter) PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), Nh);
						for (int j = 0; j < Kv; j++)
							umid[Av[j]][decodingi] = u1[index][Av[j]];
					}
					// half iteration
					////////delete SCL variables
					for (int i = 0; i < L; i++)
						delete[]LLR[i];
					delete[] LLR;
					for (int i = 0; i < L; i++)
						delete[]sum[i];
					for (int i = 0; i < L; i++)
						delete[]u1[i];
					delete[] u1;
					delete[] sum;
					delete[] PM;
					delete[] LLR_in;
					delete[] W;
					delete[] SCL_index;
					delete[] better;
					delete[] worse;
					delete[] indexpm;

				}
				//cout << "Horizon Start:\n";
				//horizontal decoding
				// #
				if (weirditer) { continue; }

				// test: LLR after horizontal decoding
				/*cout << "After horizontal decoding" << endl;
				for (int i = 0; i < Nv; i++)
				{
					for (int j = 0; j < Nh; j++)
					{
						cout << (oriLLRh[i][j] < 0) << " ";
					}
					cout << endl;
				}*/
				// test: LLR after horizontal decoding
				// inner deinterleaving
				if (flag_interleave_inner)
				{
					float* tmp_llr_to_interleave = new float[Nv];
					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < Nv; j++) { tmp_llr_to_interleave[j] = oriLLRh[j][i]; }
						randDeIntrlvPartial(intrlv_pattern_2D_inner[i], tmp_llr_to_interleave, Av, Nv, Kv);
						for (int j = 0; j < Nv; j++) { oriLLRh[j][i] = tmp_llr_to_interleave[j]; }
					}
					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < Nv; j++) { tmp_llr_to_interleave[j] = upLLR[j][i]; }
						randDeIntrlvPartial(intrlv_pattern_2D_inner[i], tmp_llr_to_interleave, Av, Nv, Kv);
						for (int j = 0; j < Nv; j++) { upLLR[j][i] = tmp_llr_to_interleave[j]; }
					}
					delete[] tmp_llr_to_interleave;
				}


				/*auto decoding_over = std::chrono::system_clock::now();
				auto duration_decoding = std::chrono::duration_cast<std::chrono::microseconds>(decoding_over - decoding_start);
				cout << endl;
				cout << "Decoding duration" << duration_decoding.count() << "us" << endl;*/

#pragma omp parallel for
				for (int decodingi = 0; decodingi < Nv; decodingi++)
				{
					// 293.70kb
					// cout << "KKKK" << endl;
					if (half_iter && (iter == softiter)) { break; }
					string polarde_now = polarde_array[iter - 1];
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

					/*			if (L == 1)
								{
									//T_start = clock();
									Count = 0;
									decode(*LLR, *sum, (int)(log(Nh) / log(2)), A_Ach, flag, Nh, Kh);
									if (iter == softiter) PolarEncode_xor(*u1, *sum, Nh);
									//T_finish = clock();
								}
								else
								{*/
								//for (int k = 0; k < L; k++)  PM[k] = 0, indexpm[k] = k; 
								////T_start = clock();
								//decode_list(LLR, sum, (int)(log(Nh) / log(2)), A_Ach, 0, Nh, Kh, PM, L, (int)(log(L) / log(2)), LLR_in, W, SCL_index, better, worse, &Countv, &Countinfov, &Qsumunv);
								////T_finish = clock();
								////Decoding_Duration += (T_finish - T_start);
								//for (int i = 0; i < L - 1; i++)
								//	for (int j = i + 1; j < L; j++)
								//		if (PM[indexpm[i]] < PM[indexpm[j]])
								//		{
								//			int ttt = indexpm[i];
								//			indexpm[i] = indexpm[j];
								//			indexpm[j] = ttt;
								//		}


							/*
							* Encoding: Info -> Horizontal Encoding -> A -> set to 0 -> B -> Horizontal Encoding -> Vertical Encoding
							* Decoding: Received Info -> Horizontal Decoding -> Vertical Decoding -> Horizontal Decoding
							*/
					if (polarde_now == "SLD")
					{
						LSD_decode_outall(u_decod, u1_LSD, LLR[0], GGh, A_Ach, Nh, L_LSD, n_paths);
						// if (iter == softiter) PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), Nh);
						// for soft information exchange, inverse coding:
						// THIS LINE is WRONG! BUT IT SEEMS BETTER THAN THE RIGHT VERSION
						for (int i = 0; i < L_LSD; i++) { PolarEncode_xor(sum_LSD[i], u1_LSD[i], Nh); }
						// for (int i = 0; i < L_LSD; i++) { PolarEncode_xor(sum_LSD[i], u1_LSD[i], Nh); }
						// if (iter == softiter) PolarEncode_xor(u1[0], sum_LSD[0], Nh);
						index = indexpm[0];
						// test: after horizontal decoding
						/*if (decodingi == Av[0])
						{
							cout << "Direct output of horizontal decoding:" << endl;
							for (int i = 0; i < Nh; i++) { cout << u1_LSD[0][i] << " "; }
							cout << endl;
							cout << "horizontal decoder out, encode once:" << endl;
							for (int i = 0; i < Nh; i++) { cout << sum_LSD[0][i] << " "; }
							cout << endl;
						}*/
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
						// iter*2-1
						if (iter < softiter)
						{
							if (adaptive_beta)
							{
								if (softmethod == "Pyndiah") Pyndiahh_LSD_adaptivebeta(LLR, sum_LSD, oriLLRh, indexpm, iter, decodingi, upLLR);
								if (softmethod == "MITSO") MITSOh(LLR, PM, sum_LSD, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunv);
							}
							else if (irregular_iteration)
							{

							}
							else
							{
								// cout << "IN LSD!!!" << endl;
								if (softmethod == "Pyndiah")
									Pyndiahh_LSD_getsum(LLR, sum_LSD, oriLLRh, indexpm, iter, decodingi, upLLR, sum_sdis_array[decodingi], GGh, Ach);
								//if (softmethod == "Pyndiah") Pyndiahh_LSD(LLR, sum_LSD, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR);
								if (softmethod == "MITSO") MITSOh(LLR, PM, sum_LSD, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunv);
							}
							//	getchar();
						}
						else
							for (int j = 0; j < Kh; j++)
								umid[decodingi][Ah[j]] = u1_LSD[index][Ah[j]];
					}
					if (polarde_now == "SCL")
					{
						for (int k = 0; k < L; k++)  PM[k] = 0, indexpm[k] = k;
						//T_start = clock();
						// decode_list(LLR, sum, (int)(log(Nh) / log(2)), A_Ach, 0, Nh, Kh, PM, L, (int)(log(L) / log(2)), LLR_in, W, SCL_index, better, worse, &Countv, &Countinfov, &Qsumunv);
						decode_list(LLR, sum, (int)(log(Nh) / log(2)), A_Ac_time_all[decodingi], 0, Nh, K_time_all[decodingi], PM, L, (int)(log(L) / log(2)), LLR_in, W, SCL_index, better, worse, &Countv, &Countinfov, &Qsumunv);
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
						if (iter == softiter) { 
							PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), Nh); 
							if (OUTPUT_ERROR_PATTERN) {
								for (int i = 0; i < Nh; i++) { x_decoder_output[decodingi][i] = sum[indexpm[0]][i]; }
							}
						}
						// since systematic encoding is used in space domain, we directly copy the best path (x) to u1 and encode on time domain when all temporal decoing is done
						/*if (iter == softiter)
						{
							for (int j = 0; j < Nh; j++) { u1[indexpm[0]][j] = sum[indexpm[0]][j]; }
						}*/
						index = indexpm[0];
						//}
						//calc distance s()
						if (iter < softiter)
						{
							if (adaptive_beta)
							{
								if (softmethod == "Pyndiah") Pyndiahh_adaptivebeta(LLR, sum, oriLLRh, indexpm, iter, decodingi, upLLR);
								if (softmethod == "MITSO") MITSOh(LLR, PM, sum, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunv);
								//	getchar();
							}
							else if (irregular_coding)
							{
								if (softmethod == "Pyndiah") Pyndiahh_irregular(LLR, sum, oriLLRh, indexpm, iter, decodingi, upLLR, sum_sdis_array[decodingi], use_short_codeword);
								if (softmethod == "MITSO") MITSOh_irregular(LLR, PM, sum, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunv, use_short_codeword);
								if (iter == softiter - 1) {
									for (int i = 0; i < Nh; i++) { upLLR_pre[decodingi][i] = upLLR[decodingi][i]; }
								}
							}
							else
							{
								//if (softmethod == "Pyndiah") Pyndiahh(LLR, sum, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR);
								if (softmethod == "Pyndiah") Pyndiahh_getsum(LLR, sum, oriLLRh, indexpm, iter, decodingi, upLLR, sum_sdis_array[decodingi]);
								if (softmethod == "MITSO") MITSOh(LLR, PM, sum, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunv);
								//	getchar();
							}
						}
						else
							for (int j = 0; j < Nh; j++)
								umid[decodingi][j] = u1[index][j];
					}

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


				if (flag_interleave_inner)
				{
					float* tmp_llr_to_interleave = new float[Nv];
					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < Nv; j++) { tmp_llr_to_interleave[j] = oriLLRh[j][i]; }
						randIntrlvPartial(intrlv_pattern_2D_inner[i], tmp_llr_to_interleave, Av, Nv, Kv);
						for (int j = 0; j < Nv; j++) { oriLLRh[j][i] = tmp_llr_to_interleave[j]; }
					}
					for (int i = 0; i < Nh; i++)
					{
						for (int j = 0; j < Nv; j++) { tmp_llr_to_interleave[j] = upLLR[j][i]; }
						randIntrlvPartial(intrlv_pattern_2D_inner[i], tmp_llr_to_interleave, Av, Nv, Kv);
						for (int j = 0; j < Nv; j++) { upLLR[j][i] = tmp_llr_to_interleave[j]; }
					}
					delete[] tmp_llr_to_interleave;
				}
			}
			int calcsum[Nv], calcu1[Nv];
			int calcu1h[Nh];
			// calcsum is from the direct output of the decoder
			// if output is the result of row decoding (meaning there is another column decoding to be done)
			if (half_iter)
			{
				for (int i = 0; i < Kv; i++)
				{
					PolarEncode_xor(calcu1h, umid[Av[i]], Nh);
					for (int j = 0; j < Kh; j++)
						uout[Av[i]][Ah[j]] = calcu1h[Ah[j]];
				}
			}
			// if output is the result of column decoding (meaning that there is another row decoding to be done)
			else
			{
				for (int i = 0; i < Kh; i++)
				{
					for (int j = 0; j < Nv; j++)
						calcsum[j] = umid[j][Ah[i]];
					PolarEncode_xor(calcu1, calcsum, Nv);
					for (int j = 0; j < Kv; j++)
						uout[Av[j]][Ah[i]] = calcu1[Av[j]];
				}
			}
			// test: final output
			/*cout << "final output" << endl;
			for (int k = 0; k < Nh; k++) { cout << uout[Av[0]][k] << " "; }
			cout << endl;*/
			// test: final output
			//now without CRC
			int Temp_Error = Num_Error;
			for (int i = 0; i < Kv; i++)
			{
				//if ((i == 24) || (i == 56)) { continue; }
				for (int j = 0; j < K_time_all[Av[i]]; j++)
				{
					//printf("%d ", (int)uout[Av[i]][Ah[j]]);
					Num_Error += abs(u[Av[i]][A_time_all[Av[i]][j]] - umid[Av[i]][A_time_all[Av[i]][j]]);
					if (abs(u[Av[i]][A_time_all[Av[i]][j]] - umid[Av[i]][A_time_all[Av[i]][j]]) == 1)
					{
						all_error_positions[Av[i]][A_time_all[Av[i]][j]] += 1;
					}
				}
				//cout << endl;
			}
			//getchar();
			if (Temp_Error != Num_Error) {
				Num_Frame_Error++; 

				if (showllr_irr) {
					std::cout << endl << "ERR" << endl;
					for (int i = 0; i < Nh; i++) { std::cout << upLLR_pre[3][i] << "  "; }
					std::cout << endl;
					for (int i = 0; i < Nh; i++) { std::cout << upLLR[3][i] << "  "; }
					std::cout << endl;
				}
				if (OUTPUT_ERROR_PATTERN) {
					int** error_pattern_tmp = new int* [Nv];
					for (int i = 0; i < Nv; i++) { error_pattern_tmp[i] = new int[Nh]; }
					for (int i = 0; i < Nv; i++) { for (int j = 0; j < Nh; j++) { error_pattern_tmp[i][j] = 0; } }
					file_error_pattern << endl;
					for (int i = 0; i < Nv; i++) {
						for (int j = 0; j < Nh; j++) {
							int err = x_decoder_output[i][j] ^ x_encoder_output[i][j];
							file_error_pattern << err << " ";
							error_pattern_tmp[i][j] = err;
						}
						file_error_pattern << endl;
					}
					if (flag_interleave)
					{
						if (flag_interleave_all)  // interleave temporally
						{
							for (int i = 0; i < Nv; i++)
							{
								randIntrlv(intrlv_pattern_2D_time[i], error_pattern_tmp[i], Nh);
							}
						}

						if (flag_interleave_partial)
						{
							for (int i = 0; i < Nv; i++)
							{
								for (int j = 0; j < intrlv_pattern_num; j++)
								{
									randIntrlv(intrlv_pattern_2D_time_partial[i][j], error_pattern_tmp[i] + j * ModType, ModType);
								}
							}
						}
						if ((!flag_interleave_block) && flag_interleave_spatial)
						{
							for (int j = 0; j < Nh; j++)
							{
								for (int i = 0; i < Nv; i++) { x_tmp[i] = error_pattern_tmp[i][j]; }
								if (K_time_array[0] == K_time_array[1]) {
									if (!flag_information_mapping)
									{
										randIntrlv(intrlv_pattern_2D[j], x_tmp, Nv);
									}
									else {
										randIntrlv(intrlv_pattern_info_mapping, x_tmp, Nv);
									}
								}
								else {
									if (!flag_separate_interleaving)
									{
										randIntrlv(non_uniform_intrlv_pattern, x_tmp, Nv);
									}
									else {
										randIntrlv(non_uniform_intrlv_pattern_array[j], x_tmp, Nv);
									}
								}
								for (int i = 0; i < Nv; i++) { error_pattern_tmp[i][j] = x_tmp[i]; }
							}
						}
					}
					file_error_pattern << endl;
					float sum_sing = 0.0;
					int n_err_pos = 0.0;
					for (int i = 0; i < Nv; i++) {
						for (int j = 0; j < Nh; j++) {
							if (error_pattern_tmp[i][j]) {
								file_error_pattern << singular_values_array[i] * sqrt(antenna_power[i]) << "  ";
								sum_sing += singular_values_array[i];
								n_err_pos += 1;
							}
						}
					}
					file_error_pattern << endl << sum_sing / n_err_pos << endl;
					file_error_pattern << n_err_pos << endl;


				}

			}
			if (frame % 10 == 0)
			{
				/*printf("\r<%d> SNR= %.2f FER= %d / %d = %f BER= %d / %d = %lf E-Det = %f, E-MRB = %f"  ,
					frame, nSNR, Num_Frame_Error, frame, (double)Num_Frame_Error / frame, Num_Error, (N_Info * frame), (double)Num_Error / N_Info / frame,
					(double)n_det_err/total_bit_cnt, (double)n_det_err_mrb/(total_bit_cnt/(ModType/2)));*/
				int L_sdis_sum = L_LSD;
				sdis_sum = 0;
				for (int i = 0; i < Nv; i++) sdis_sum += sum_sdis_array[i];
				if (polarde_array[0] == "SCL") L_sdis_sum = L;
				// cout << "sdis_sum_avg" << sdis_sum / frame / L_sdis_sum << endl;
				printf("\r<%d> SNR= %.2f FER= %d / %d = %f BER= %d / %d = %lf ,sdis_avg = %f, E-MRB = %f",
					frame, nSNR, Num_Frame_Error, frame, (double)Num_Frame_Error / frame, Num_Error, (N_Info * frame), (double)Num_Error / N_Info / frame, (double)sdis_sum / frame / L_sdis_sum, (double)n_det_err / total_bit_cnt);
				fflush(stdout);

				//if (frame > 0) {
				//	// std::cout << endl << (float)trans_energy / (Nsym * Nv* frame) << endl;
				//	for (int i = 0; i < Nv; i++) {
				//		std::cout << endl;
				//		for (int j = 0; j < Nh; j++) {
				//			std::cout << (x_encoder_sum[i][j] > 0) << " ";
				//		}
				//	}
				//}
				//std::cout << "TEST ENDS" << endl;
			}


			// output interleaving pattern
			//if (Num_Frame_Error % record_frames == 0 && Num_Frame_Error > 0 && Num_Frame_Error!=n_frame_error_pre)
			//{
			//	n_frame_error_pre = Num_Frame_Error;
			//	cout << endl;
			//	int n_file = Num_Frame_Error / record_frames;
			//	frame_delta = frame - frame_tmp;
			//	cout << "Delta Frames: " << frame_delta << endl;
			//	frame_tmp = frame;
			//	string ofile_intrlv_name = string("intrlv_pattern_") + to_string(n_file) + string("_frame_") + to_string(frame_delta) + string(".txt");
			//	string ofile_spatial_intrlv_name = string("intrlv_pattern_spatial") + to_string(n_file) + string(".txt");
			//	ofstream ofile_intrlv(ofile_intrlv_name);
			//	ofstream ofile_spatial_intrlv(ofile_spatial_intrlv_name);
			//	for (int i = 0; i < Nv; i++)
			//	{
			//		for (int j = 0; j < intrlv_pattern_num; j++)
			//		{

			//			for (int k = 0; k < ModType; k++)
			//				ofile_intrlv << intrlv_pattern_2D_time_partial[i][j][k] << " ";
			//			//cout << endl;
			//		}
			//	}
			//	ofile_intrlv.close();

			//	for (int i = 0; i < Nh; i++)
			//	{
			//		for (int j = 0; j < Nv; j++)
			//			ofile_spatial_intrlv << intrlv_pattern_2D[i][j] << " ";
			//	}
			//	ofile_spatial_intrlv.close();

			//	// re-generate interleaving pattern
			//	generateIntrlvPattern(intrlv_pattern, Nv);
			//	for (int i = 0; i < Nh; i++) { generateIntrlvPattern(intrlv_pattern_2D[i], Nv); }
			//	for (int i = 0; i < Nv; i++) { generateIntrlvPattern(intrlv_pattern_2D_time[i], Nh); }
			//	if (flag_interleave_partial)
			//	{
			//		for (int i = 0; i < Nv; i++)
			//		{
			//			for (int j = 0; j < intrlv_pattern_num; j++)
			//			{
			//				generateIntrlvPattern(intrlv_pattern_2D_time_partial[i][j], ModType);
			//				/*for (int k = 0; k < ModType; k++)
			//					cout << intrlv_pattern_2D_time_partial[i][j][k];
			//				cout << endl;*/
			//			}
			//		}
			//	}
			//	//Num_Frame_Error = 0;
			//}


			//printf("%d %d", Num_Frame_Error, Num_Error);
			//getchar();
			//Decoding_Duration += (T_finish - T_start);
			/*auto dec_end = std::chrono::system_clock::now();
			auto duration_decoding = std::chrono::duration_cast<std::chrono::microseconds>(dec_end - det_end);
			auto duration_total = std::chrono::duration_cast<std::chrono::microseconds>(dec_end - start);
			file_dec << duration_decoding.count() << endl;
			file_total << duration_total.count() << endl;*/
			if (frame % 16 == 0)
			{
				file_err_num << Num_Error << endl;
			}
			if (frame % 500 == 0)
			{
				std::ofstream ofile_err("all_errors_log.txt");
				for (int i = 0; i < Nv; i++)
				{
					for (int j = 0; j < Nmax; j++)
					{
						ofile_err << all_error_positions[i][j] << " ";
					}
					ofile_err << endl;
				}
				ofile_err << endl;
				ofile_err.close();
			}

		}
		//Decoding_Duration = Decoding_Duration * 100 / CLOCKS_PER_SEC / frame;
		//Decoding_Duration = Decoding_Duration / CLOCKS_PER_SEC / Max_frame;

		//printf("%f  %f  %f  %f\n", nSNR, Decoding_Duration, (double)Num_Error / N_Info / Max_frame, (double)Num_Frame_Error / Max_frame);
		printf("%f %f %f\n", nSNR, (double)Num_Frame_Error / frame, (double)Num_Error / N_Info / frame);
		BER_array[SNR_count] = (double)Num_Error / N_Info / frame;
		FER_array[SNR_count] = (double)Num_Frame_Error / frame;
		// save to .txt file
		string polarde_serial = "";
		for (int i = 0; i < softiter; i++)
			polarde_serial = polarde_serial + polarde_array[i] + string("_");
		string file_name = string("ResIrr0922\\Irregular_") + string("Border") + to_string(border) + string("N") + to_string(Nh) + string("K0") + to_string(K_time_array[0]) + string("K1") + to_string(K_time_array[1])
			+ string("T") + to_string(softiter) + string("H") + to_string(half_iter) + string("UB") + to_string(usebeta) + string("IS") + to_string(flag_interleave_spatial)
			+ string("DG") + to_string(diff_gain) + string(".txt");
		save_array_txt(nSNR, FER_array[SNR_count], file_name);
		int L_sdis_sum = L_LSD;
		sdis_sum = 0;
		for (int i = 0; i < Nv; i++) sdis_sum += sum_sdis_array[i];
		if (polarde_array[0] == "SCL") L_sdis_sum = L;
		// cout << "sdis_sum_avg" << sdis_sum / frame / L_sdis_sum << endl;
		SNR_count++;
		delete[] non_uniform_intrlv_pattern_array;
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
	for (int i = 0; i < Nh; i++)
		delete[] oriLLR_kbest[i];
	delete[] oriLLR_kbest;
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
	for (int i = 0; i < mimo_blk_height; i++) {
		delete[] x_mod[i];
	}
	delete[] x_mod;
	for (int i = 0; i < Nr; i++) {
		delete[] y_mod[i];
	}
	delete[] y_mod;
	for (int i = 0; i < 2 * mimo_blk_height; i++) {
		delete[] y_mmse[i];
	}
	delete[] y_mmse;
	for (int i = 0; i < 2 * Nr; i++) {
		delete[] y_mod_real[i];
	}
	for (int i = 0; i < mimo_blk_length; i++)
	{
		delete[] y_mod_real_temp_omp[i];
	}
	delete[] y_mod_real_temp_omp;
	for (int i = 0; i < mimo_blk_length; i++)
	{
		delete[] y_mmse_temp_omp[i];
	}
	delete[] y_mmse_temp_omp;
	for (int i = 0; i < mimo_blk_length; i++)
	{
		delete[] ep_extr_mean_omp[i];
	}
	delete[] ep_extr_mean_omp;
	for (int i = 0; i < mimo_blk_length; i++)
	{
		delete[] ep_extr_var_omp[i];
	}
	delete[] ep_extr_var_omp;
	delete[] y_mod_real;
	delete[] x_mod_tmp;
	delete[] y_mod_tmp;
	delete[] y_mod_real_tmp;
	delete[] x_tmp;
	delete[] oriLLR_tmp;
	delete[] sum_sdis_array;
	delete[] Hreal_array;
	delete[] HTH_array;
	delete[] HTY_array;
	delete[] intrlv_pattern_2D_inner;
}





// puncturing
//void system2D_LSD_epdet_with_CRC_punctured(double* BER_array, double* FER_array)
//{
//	std::ofstream file_crc_check_replacement_codeword("crc_replacement_codeword_log.txt", std::ios::app);
//	std::ofstream file_crc_check_decoding("crc_check_decoding.txt", std::ios::app);
//	std::ofstream file_total("total_timing.txt", std::ios::app);
//	std::ofstream file_det("det_timing.txt", std::ios::app);
//	std::ofstream file_dec("dec_timing.txt", std::ios::app);
//	std::ofstream file_err_num("total_frames.txt", std::ios::app);
//	auto det_end = std::chrono::system_clock::now();
//	int Nt2 = 2 * Nt;
//	int Nr2 = 2 * Nr;
//
//	/*if (!file.is_open()) {
//		std::cerr << "Error opening file: " << file_name << std::endl;
//		return;
//	}*/
//
//	//file << setw(20) << SNR << setw(20) << BER << setw(20) << FER << "\n";
//
//	//freopen("result.txt", "w", stdout);
//	printsetting();
//	cout << "\nSNR	 " << "FER      BER" << endl;
//	int N_Info = Kv * Kh - CRCL, Nmax;
//	//Interleaving pattern
//	/*std::vector<int> intrlv_pattern(Nv);
//	std::vector<int>* intrlv_pattern_2D = new vector<int>[Nh];
//	for (int i = 0; i < Nh; i++) intrlv_pattern_2D[i].resize(Nv);*/
//	int intrlv_pattern_num = Nh / ModType;
//	std::vector<int> intrlv_pattern(Nv);
//	std::vector<int>* intrlv_pattern_2D = new vector<int>[Nh];
//	std::vector<int>* intrlv_pattern_2D_time = new vector<int>[Nv];
//	std::vector<int>* intrlv_pattern_2D_inner = new vector<int>[Nh];
//	std::vector<int>** intrlv_pattern_2D_time_partial = new vector<int> *[Nv];
//	for (int i = 0; i < Nv; i++) intrlv_pattern_2D_time_partial[i] = new vector<int>[intrlv_pattern_num];
//	for (int i = 0; i < Nv; i++)
//	{
//		for (int j = 0; j < intrlv_pattern_num; j++)
//		{
//			intrlv_pattern_2D_time_partial[i][j].resize(ModType);
//		}
//	}
//	for (int i = 0; i < Nh; i++) intrlv_pattern_2D[i].resize(Nv);
//	for (int i = 0; i < Nh; i++) intrlv_pattern_2D_inner[i].resize(Kv);
//	for (int i = 0; i < Nv; i++) intrlv_pattern_2D_time[i].resize(Nh);
//
//	float CodeRate = (float)N_Info / (Nh * Nv);
//	Nmax = max(Nh, Nv);
//	// num of symbols per vertical codeword
//	int Nsym = 0;
//	if (MOD_DIM == "Spatial") { Nsym = Nv / ModType; }
//	else if (MOD_DIM == "Temporal") { Nsym = Nh / ModType; }
//	int symbol_length = ModType;
//
//	// Size of a mimo block is given by height(spatial) * length(temporal)
//	int mimo_blk_height = 0;
//	int mimo_blk_length = 0;
//	if (MOD_DIM == "Spatial")
//	{
//		mimo_blk_height = Nsym;
//		mimo_blk_length = Nh;
//	}
//	else if (MOD_DIM == "Temporal")
//	{
//		mimo_blk_height = Nv;
//		mimo_blk_length = Nsym;
//	}
//	// Initialization
//	unsigned char* a_data = new unsigned char[Kv * Kh];// Information CodeWord
//	int** u = new int* [Nv];// Original CodeWord u
//	for (int i = 0; i < Nv; i++) { u[i] = new int[Nh]; }
//	int** umid = new int* [Nv];// output CodeWord u
//	for (int i = 0; i < Nv; i++) { umid[i] = new int[Nh]; }
//	int** uout = new int* [Nv];// output CodeWord u
//	for (int i = 0; i < Nv; i++) { uout[i] = new int[Nh]; }
//	int** midx = new int* [Nv];//  Polar Encoded CodeWord x
//	for (int i = 0; i < Nv; i++) { midx[i] = new int[Nh]; }
//	int** x = new int* [Nv];//  Polar Encoded CodeWord x
//	for (int i = 0; i < Nv; i++) { x[i] = new int[Nh]; }
//	int** x_encoded_copy = new int* [Nv];
//	for (int i = 0; i < Nv; i++) { x_encoded_copy[i] = new int[Nh]; }
//
//	// Arrays for MIMO
//	int** x_mmse_out = new int* [Nv]; // delete!
//	for (int i = 0; i < Nv; i++) { x_mmse_out[i] = new int[Nh]; }
//	int** x_parity_check = new int* [Nv];
//	for (int i = 0; i < Nv; i++) { x_parity_check[i] = new int[Nh]; }
//
//	//std::complex<float>** x_mod = new std::complex<float>* [Nsym]; // modulated symbols (spatial modulation) DELETE!
//	//for (int i = 0; i < Nsym; i++) { x_mod[i] = new std::complex<float>[Nh]; }
//	//std::complex<float>** x_mod_temporal = new std::complex<float>*[Nv]; // modulated symbols (temporal modulation) DELETE!
//	//for (int i = 0; i < Nv; i++) x_mod_temporal[i] = new std::complex<float>[Nsym];
//	//std::complex<float>** y_mod = new std::complex<float>* [Nr]; // Output symbols after AWGN channels // delete!
//	//for (int i = 0; i < Nr; i++) { y_mod[i] = new std::complex<float>[Nh]; }
//	//float** y_mod_real = new float* [2 * Nr];
//	//for (int i = 0; i < 2*Nr; i++) y_mod_real[i] = new float[Nh];
//	//float** y = new float* [Nv];// Output Signal After AWGN Channel
//	//for (int i = 0; i < Nv; i++) { y[i] = new float[Nh]; }
//	//float** y_mmse = new float* [2*Nsym];
//	//for (int i = 0; i < 2*Nsym; i++) { y_mmse[i] = new float[Nh]; }  // MMSE detector output // delete!
//	//float* y_mmse_tmp = new float[2*Nsym];
//	//std::complex<float>* x_mod_tmp = new std::complex<float>[Nsym]; // modulated symbols of a single column // delete!
//	//std::complex<float>* y_mod_tmp = new std::complex<float>[Nr]; // received symbols of a single column // delete!
//	//float* y_mod_real_tmp = new float[2*Nr];
//	//// float* y_mod_tmp = new float[Nr];  // delete!
//
//
//	std::complex<float>** x_mod = new std::complex<float>*[mimo_blk_height]; // modulated symbols (spatial modulation) DELETE!
//	for (int i = 0; i < mimo_blk_height; i++) { x_mod[i] = new std::complex<float>[mimo_blk_length]; }
//	std::complex<float>** y_mod = new std::complex<float>*[Nr]; // Output symbols after AWGN channels // delete!
//	for (int i = 0; i < Nr; i++) { y_mod[i] = new std::complex<float>[mimo_blk_length]; }
//	float** y_mod_real = new float* [2 * Nr];
//	for (int i = 0; i < 2 * Nr; i++) y_mod_real[i] = new float[mimo_blk_length];
//	float** y = new float* [Nv];// Output Signal After AWGN Channel
//	for (int i = 0; i < Nv; i++) { y[i] = new float[Nh]; }
//	float** y_mmse = new float* [2 * mimo_blk_height];
//	for (int i = 0; i < 2 * mimo_blk_height; i++) { y_mmse[i] = new float[mimo_blk_length]; }  // MMSE detector output // delete!
//	float* y_mmse_tmp = new float[2 * mimo_blk_height];
//	std::complex<float>* x_mod_tmp = new std::complex<float>[mimo_blk_height]; // modulated symbols of a single column // delete!
//	std::complex<float>* y_mod_tmp = new std::complex<float>[Nr]; // received symbols of a single column // delete!
//	float* y_mod_real_tmp = new float[2 * Nr];
//	// float* y_mod_tmp = new float[Nr];  // delete!
//	float** y_mod_real_temp_omp = new float* [mimo_blk_length]; // new delete!
//	for (int i = 0; i < mimo_blk_length; i++) { y_mod_real_temp_omp[i] = new float[2 * Nr]; }
//	float** y_mmse_temp_omp = new float* [mimo_blk_length]; // new delete!
//	for (int i = 0; i < mimo_blk_length; i++) { y_mmse_temp_omp[i] = new float[2 * mimo_blk_height]; }
//	// ep
//	float** ep_extr_mean_omp = new float* [mimo_blk_length];
//	for (int i = 0; i < mimo_blk_length; i++) { ep_extr_mean_omp[i] = new float[2 * mimo_blk_height]; }
//	float** ep_extr_var_omp = new float* [mimo_blk_length];
//	for (int i = 0; i < mimo_blk_length; i++) { ep_extr_var_omp[i] = new float[2 * mimo_blk_height]; }
//	int* x_tmp = new int[Nv]; // delete!
//
//
//	// paths selected by k-best detection
//	// array for storing results of bit-flip
//	int n_samples = pow(2, p);
//	int*** tried_codewords_h = new int** [Nv];  // ��һά��Nv��ָ��
//
//	for (int i = 0; i < Nv; ++i) {
//		tried_codewords_h[i] = new int* [n_samples];  // �ڶ�ά��ÿ����һάԪ��ָ��n_samples��ָ��
//
//		for (int j = 0; j < n_samples; ++j) {
//			// ����ά��ÿ���ڶ�άԪ��ָ��Nh������
//			tried_codewords_h[i][j] = new int[Nh]();  // ĩβ��()����ֵ��ʼ����0��
//		}
//	}
//
//	int*** tried_codewords_v = new int** [Nh];  // ��һά��Nv��ָ��
//
//	for (int i = 0; i < Nh; ++i) {
//		tried_codewords_v[i] = new int* [n_samples];  // �ڶ�ά��ÿ����һάԪ��ָ��n_samples��ָ��
//
//		for (int j = 0; j < n_samples; ++j) {
//			// ����ά��ÿ���ڶ�άԪ��ָ��Nh������
//			tried_codewords_v[i][j] = new int[Nv]();  // ĩβ��()����ֵ��ʼ����0��
//		}
//	}
//
//	int* n_tried_patterns_h_all = new int[softiter];
//	for (int i = 0; i < softiter; i++) n_tried_patterns_h_all[i] = 0;
//	int* n_tried_patterns_v_all = new int[softiter];
//	for (int i = 0; i < softiter; i++) n_tried_patterns_v_all[i] = 0;
//
//	// CRC bits
//	unsigned char* u_crc = new unsigned char[Kv * Kh];
//
//	float** upLLR = new float* [Nv];
//	for (int i = 0; i < Nv; i++) { upLLR[i] = new _MM_ALIGN16 float[Nh]; }
//	//float** oriLLRv = new float* [Nh];
//	//for (int i = 0; i < Nh; i++) { oriLLRv[i] = new _MM_ALIGN16 float[2 * Nv + 2]; }
//	float** oriLLRh = new float* [Nv];
//	for (int i = 0; i < Nv; i++) { oriLLRh[i] = new _MM_ALIGN16 float[Nh]; }
//	int len_ori_llr_tmp = 0;
//	if (MOD_DIM == "Spatial") { len_ori_llr_tmp = Nv; }
//	else if (MOD_DIM == "Temporal") { len_ori_llr_tmp = ModType * Nt; }
//	float* oriLLR_tmp = new _MM_ALIGN16 float[len_ori_llr_tmp];  // delete!
//	//float** oriLLR_kbest = new float* [Nh];
//	float** oriLLR_kbest = new float* [mimo_blk_length];
//	for (int i = 0; i < Nh; i++) { oriLLR_kbest[i] = new float[len_ori_llr_tmp]; }
//	double* sum_sdis_array = new double[Nv];
//	for (int i = 0; i < Nv; i++) sum_sdis_array[i] = 0;
//
//	//*****************************************************************
//	// EP detector
//	float* Hreal_array = new float[(2 * Nr) * (2 * Nt)]; // delete
//	float* HTH_array = new float[(2 * Nt) * (2 * Nt)]; // delete
//	float* HTY_array = new float[2 * Nt]; // delete 
//
//
//
//
//
//	int** GGh = new int* [Nh];
//	int** GGv = new int* [Nv];
//
//
//	MatrixXf Gh = Fn(Nh);  // Generate Matrix of horizontal coding
//	MatrixXf Gv = Fn(Nv);  // Generate Matrix of vertical coding
//	MatrixXcf H(Nr, Nt); // Complex channel Matrix H: Nr x Nt
//	//Matrix<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic> H = MatrixXcf::Random(Nr, Nt); // Real domain channel matrix
//	//// //MatrixX Definitions
//	MatrixXf H_real(2 * Nr, 2 * Nt); // Real domain channel matrix
//	MatrixXf received_signal(2 * Nr, 1);
//	MatrixXf HTH(2 * Nt, 2 * Nt);
//	MatrixXf Identity = MatrixXf::Identity(2 * Nt, 2 * Nt);
//	MatrixXf output_signal(2 * Nt, 1);
//	MatrixXf mmse_matrix(2 * Nt, 2 * Nr);
//	MatrixXf HTHIinv(2 * Nt, 2 * Nt);
//
//	// for KBest
//	MatrixXf R(2 * Nt, 2 * Nt);
//	MatrixXf Q(2 * Nr, 2 * Nt);
//	MatrixXf z(2 * Nt, 1);
//
//	// for CRC
//	int** whole_codeword_list = new int* [L_WHOLE_CODEWORD];
//	for (int i = 0; i < L_WHOLE_CODEWORD; i++) { whole_codeword_list[i] = new int[Nh * Nv]; }
//
//	int** whole_codeword_list_v = new int* [L_WHOLE_CODEWORD];
//	for (int i = 0; i < L_WHOLE_CODEWORD; i++) { whole_codeword_list_v[i] = new int[Nh * Nv]; }
//
//	int* n_err_codeword_v = new int[Nh];  // err of codeword at each column
//
//	int** n_err_codeword_v_diffpos = new int* [Nh];
//	for (int i = 0; i < Nh; i++) n_err_codeword_v_diffpos[i] = new int[L + 1];
//	for (int i = 0; i < Nh; i++)
//	{
//		for (int j = 0; j < L + 1; j++) { n_err_codeword_v_diffpos[i][j] = 0; }
//	}
//
//	int* correct_idx_oneframe = new int[Nh];
//	for (int i = 0; i < Nh; i++) { correct_idx_oneframe[i] = 0; }
//
//	int* max_L_cnt = new int[L + 1];
//	for (int i = 0; i < L + 1; i++) { max_L_cnt[i] = 0; }
//
//	// for CRC
//
//	int L_max = max(L_LSD, L);
//
//	// stores the decoding result of Nv horizontal codewords
//	int** sub_codeword_list = new int* [Nv];
//	for (int i = 0; i < Nv; i++) { sub_codeword_list[i] = new int[L_max * Nh]; }
//	int** sub_codeword_list_v = new int* [Nh];
//	for (int i = 0; i < Nh; i++) { sub_codeword_list_v[i] = new int[L_max * Nv]; }
//	float** edistance_list = new float* [Nv];
//	for (int i = 0; i < Nv; i++) { edistance_list[i] = new float[L_max]; }
//	float** edistance_list_v = new float* [Nh];
//	for (int i = 0; i < Nh; i++) edistance_list_v[i] = new float[L_max];
//	int** sub_encoded_codeword_list = new int* [Nv];
//	for (int i = 0; i < Nv; i++) { sub_encoded_codeword_list[i] = new int[L_max * Nh]; }
//	/*int** sub_encoded_codeword_list_v = new int* [Nh];
//	for (int i = 0; i < Nh; i++) { sub_encoded_codeword_list_v[i] = new int[L_max * Nv]; }*/
//	int** sub_encoded_codeword_list_v = new int* [Nh];
//	for (int i = 0; i < Nh; i++) { sub_encoded_codeword_list_v[i] = new int[L_max * Nv]; }
//	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> H_real = MatrixXf::Random(2 * Nr, 2 * Nt); // Real domain channel matrix
//	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> received_signal = MatrixXf::Random(2 * Nr,1);
//	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> HTH = MatrixXf::Random(2 * Nt, 2 * Nt);
//	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> Identity = MatrixXf::Random(2 * Nt, 2 * Nt);
//	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> output_signal = MatrixXf::Random(2 * Nt, 1);
//	//Matrix<float, Eigen::Dynamic, Eigen::Dynamic> mmse_matrix = MatrixXf::Random(2 * Nt, 2 * Nr);
//
//
//	// Matrix Definitions
//	//Eigen::Matrix<float, 2 * Nr, 2 * Nt> H_real; // ʵ���ŵ�����
//	//Eigen::Matrix<float, 2 * Nr, 1> received_signal; // �����źž���
//	//Eigen::Matrix<float, 2 * Nt, 2 * Nt> HTH; // HTH����
//	//Eigen::Matrix<float, 2 * Nt, 2 * Nt> Identity = Eigen::Matrix<float, 2 * Nt, 2 * Nt>::Identity(); // ��λ����
//	//Eigen::Matrix<float, 2 * Nt, 1> output_signal; // ����źž���
//	//Eigen::Matrix<float, 2 * Nt, 2 * Nr> mmse_matrix; // MMSE����
//
//
//	for (int i = 0; i < Nh; i++) GGh[i] = new int[Nh];
//	for (int i = 0; i < Nv; i++) GGv[i] = new int[Nv];
//	for (int i = 0; i < Nh; i++)
//		for (int j = 0; j < Nh; j++)
//			GGh[i][j] = Gh(i, j);
//	for (int i = 0; i < Nv; i++)
//		for (int j = 0; j < Nv; j++)
//			GGv[i][j] = Gv(i, j);
//	//*****************************************************************
//	vector<int> temph(Nh);
//	vector<int> tempv(Nv);
//	int* Ah = new int[Kh]; // Horizontal Info Bits
//	int* Ach = new int[Nh - Kh]; // Horizontal Frozen Bits
//	int* A_Ach = new int[Nh];
//	int* Av = new int[Kv];
//	int* Acv = new int[Nv - Kv];
//	int* A_Acv = new int[Nv];
//
//	//*****************************************************************
//		// Main Function
//	// srand((unsigned int)time(0));
//	std::random_device rd;              // ��ȡ���������
//	std::mt19937 gen(rd());             // ʹ�����ӳ�ʼ��mt19937������
//	std::uniform_int_distribution<> distrib(1, 1000000);
//	int Num_Error, Num_Frame_Error;
//	float Decoding_Duration;
//
//	int Num_Error_new, Num_Frame_Error_new;
//	float Decoding_Duration_new;
//
//	clock_t T_start, T_finish;
//	int SNR_count = 0;
//	int anssc, anssd;
//	for (double nSNR = Min_SNR; nSNR <= Max_SNR; nSNR += SNR_step)
//	{
//		for (int i = 0; i < softiter; i++) { n_tried_patterns_h_all[i] = 0; n_tried_patterns_v_all[i] = 0; }
//		int* n_tried_patterns_h = new int[Nv];
//		int* n_tried_patterns_v = new int[Nh];
//		for (int i = 0; i < Nv; i++) { n_tried_patterns_h[i] = 0; }
//		for (int i = 0; i < Nh; i++) { n_tried_patterns_v[i] = 0; }
//		// crc check logs
//		int n_err_frame_pass_crc = 0;  // err frames not detected by crc
//		int n_err_frame_corrected_by_crc = 0; // err frames detected and corrected by crc
//		int n_err_frame_uncorrected_by_crc = 0; // err frames detected but not corrected by crc
//		int n_err_replacement_codeword_crc = 0;  // err of replacement codewords
//		int n_err_pass_parity_check = 0; // n of err that passes parity check
//		int n_err_frame_withno_replacement_codeword = 0; // n of err frames with no possible replacement codeword from the sub-lists
//		int idx_replacement_codeword = 0; // idx of replacement codeword that passes crc check
//
//		for (int i = 0; i < Nh; i++) { n_err_codeword_v[i] = 0; }
//
//		for (int i = 0; i < Nh; i++)
//		{
//			for (int j = 0; j < L + 1; j++) { n_err_codeword_v_diffpos[i][j] = 0; }
//		}
//		for (int i = 0; i < L + 1; i++) { max_L_cnt[i] = 0; }
//		for (int i = 0; i < Nh; i++) { correct_idx_oneframe[i] = 0; }
//
//		for (int i = 0; i < Nv; i++) sum_sdis_array[i] = 0;
//		double sdis_sum = 0;
//		// TEST: det err
//		int n_det_err = 0; // n error bits at detector output
//		int n_det_err_mrb = 0; // error bits at detector output(MRBs)
//		int total_bit_cnt = 0;
//		// TEST: det err
//		int cntcnt = 0;
//		//Initialization
//		Num_Error = 0;
//		Num_Frame_Error = 0;
//		Decoding_Duration = 0;
//
//		// double tmp = pow(10, nSNR / 10);
//		float sigma_sq = 0;
//		double tmp = 0;
//		if (atten == 1)
//		{
//			tmp = pow(10, nSNR / 10);
//			sigma_sq = (float)1 / (float)tmp / 2 / CodeRate;
//		}
//		if (atten == 2)
//		{
//			tmp = pow(10, nSNR / 10) * symbol_length * CodeRate * Nt;
//			// tmp = pow(10, nSNR / 10) * symbol_length  * Nt;
//			sigma_sq = (Nt * Nr) / tmp;
//		}
//		polar_codeconstruction(polarconh, Nh, Kh, GGh, (float)sqrt(sigma_sq), temph);
//		polar_codeconstruction(polarconv, Nv, Kv, GGv, (float)sqrt(sigma_sq), tempv);
//
//		// sort temph and tempv in ascending order
//		for (int i = 0; i < Kh - 1; i++)
//			for (int j = i + 1; j < Kh; j++)
//				if (temph[i] > temph[j])
//				{
//					int ttt = temph[i];
//					temph[i] = temph[j];
//					temph[j] = ttt;
//				}
//		for (int i = 0; i < Kv - 1; i++)
//			for (int j = i + 1; j < Kv; j++)
//				if (tempv[i] > tempv[j])
//				{
//					int ttt = tempv[i];
//					tempv[i] = tempv[j];
//					tempv[j] = ttt;
//				}
//
//		// Decide the position of information bits and frozen bits
//		for (int i = 0; i < Kh; i++) // Information Set
//		{
//			Ah[i] = temph[i];
//			A_Ach[Ah[i]] = 1;
//		}
//		//	getchar();
//		for (int i = 0; i < Nh - Kh; i++) // Frozen Bit
//		{
//			Ach[i] = temph[Kh + i];
//			A_Ach[Ach[i]] = 0;
//		}
//		for (int i = 0; i < Kv; i++) // Information Set
//		{
//			Av[i] = tempv[i];
//			A_Acv[Av[i]] = 1;
//		}
//		//	getchar();
//		for (int i = 0; i < Nv - Kv; i++) // Frozen Bit
//		{
//			Acv[i] = tempv[Kv + i];
//			A_Acv[Acv[i]] = 0;
//		}
//
//		for (int i = 0; i < Nh; i++)
//		{
//			cout << A_Ach[i] << ",";
//		}
//		cout << endl;
//		// test: Horizontal frozen bits
//
//		/*cout << "Frozen Bits" <<endl;
//		for (int i = 0; i < Nh - Kh; i++)
//		{
//			cout << setw(5) << Acv[i];
//		}
//		cout << endl;*/
//		// test: Horizontal frozen bits
//		// Frame Loop
//		int passnum;
//		int frame = 0;
//
//		// generate interleaving patterns
//		if (flag_interleave)
//		{
//			string ofile_intrlv_name = string("intrlv_pattern_71_frame_49263.txt");
//			string ofile_spatial_intrlv_name = string("intrlv_pattern_spatial71.txt");
//			ifstream ofile_intrlv(ofile_intrlv_name);
//			ifstream ofile_spatial_intrlv(ofile_spatial_intrlv_name);
//			int file_in = 1;
//			//for (int i = 0; i < Nv; i++)
//			//{
//			//	for (int j = 0; j < intrlv_pattern_num; j++)
//			//	{
//
//			//		for (int k = 0; k < ModType; k++)
//			//			// ofile_intrlv >> file_in;
//			//			ofile_intrlv >> intrlv_pattern_2D_time_partial[i][j][k];
//			//		//cout << endl;
//			//	}
//			//}
//			//ofile_intrlv.close();
//
//			for (int i = 0; i < Nh; i++)
//			{
//				for (int j = 0; j < Nv; j++)
//					ofile_spatial_intrlv >> intrlv_pattern_2D[i][j];
//			}
//			ofile_spatial_intrlv.close();
//			/*for (int i = 0; i < Nh; i++)
//			{
//				for (int j = 0; j < Nv; j++)
//					cout << intrlv_pattern_2D[i][j] << endl;
//			}*/
//			generateIntrlvPattern(intrlv_pattern, Nv);
//			for (int i = 0; i < Nh; i++) { generateIntrlvPattern(intrlv_pattern_2D[i], Nv); }
//			for (int i = 0; i < Nv; i++) { generateIntrlvPattern(intrlv_pattern_2D_time[i], Nh); }
//			if (flag_interleave_partial)
//			{
//				for (int i = 0; i < Nv; i++)
//				{
//					for (int j = 0; j < intrlv_pattern_num; j++)
//					{
//						generateIntrlvPattern(intrlv_pattern_2D_time_partial[i][j], ModType);
//						/*for (int k = 0; k < ModType; k++)
//							cout << intrlv_pattern_2D_time_partial[i][j][k];
//						cout << endl;*/
//					}
//				}
//			}
//			for (int i = 0; i < Nh; i++) { generateIntrlvPattern(intrlv_pattern_2D_inner[i], Kv); }
//		}
//
//		//ofstream ofile_intrlv("intrlv_pattern.txt");
//		//for (int i=0; i<Nv; i++)
//		//{
//		//	for (int j = 0; j < intrlv_pattern_num; j+for+)
//		//	{
//		//		
//		//		for (int k = 0; k < ModType; k++)
//		//			ofile_intrlv << intrlv_pattern_2D_time_partial[i][j][k] << "   ";
//		//		//cout << endl;
//		//	}
//		//}
//		//ofile_intrlv.close();
//		int frame_tmp = 0;
//		int frame_delta = 0;
//		int n_frame_error_pre = 0;
//		while (Num_Frame_Error < errframe)
//		{
//			for (int i = 0; i < Nv; i++) { n_tried_patterns_h[i] = 0; }
//			for (int i = 0; i < Nh; i++) { n_tried_patterns_v[i] = 0; }
//			/*test: processing time*/
//			//auto start = std::chrono::system_clock::now();
//			/*if (frame == 10) {
//				cout << "AAAAAAAAAA" << endl;
//			}*/
//			frame++;
//			bool flag_pass_crc = 0;
//			bool flag_pass_crc_undetected = 0;
//			bool flag_pass_crc_corrected = 0;
//			bool flag_found_replacement_codeword = 0;
//			bool flag_pass_parity_check_in_crc = 0;
//
//			generateIntrlvPattern(intrlv_pattern, Nv);
//			for (int i = 0; i < Nh; i++) { generateIntrlvPattern(intrlv_pattern_2D[i], Nv); }
//			for (int i = 0; i < Nv; i++) { generateIntrlvPattern(intrlv_pattern_2D_time[i], Nh); }
//			if (flag_interleave_partial)
//			{
//				for (int i = 0; i < Nv; i++)
//				{
//					for (int j = 0; j < intrlv_pattern_num; j++)
//					{
//						generateIntrlvPattern(intrlv_pattern_2D_time_partial[i][j], ModType);
//						/*for (int k = 0; k < ModType; k++)
//							cout << intrlv_pattern_2D_time_partial[i][j][k];
//						cout << endl;*/
//					}
//				}
//			}
//			for (int i = 0; i < Nh; i++) { generateIntrlvPattern(intrlv_pattern_2D_inner[i], Kv); }
//			// generate interleaving pattern
//			/*if (flag_interleave)
//			{
//				generateIntrlvPattern(intrlv_pattern, Nv);
//				for (int i = 0; i < Nh; i++) { generateIntrlvPattern(intrlv_pattern_2D[i], Nv); }
//				for (int i = 0; i < Nv; i++) { generateIntrlvPattern(intrlv_pattern_2D_time[i], Nh); }
//				if (flag_interleave_partial)
//				{
//					for (int i = 0; i < Nv; i++)
//					{
//						for (int j = 0; j < intrlv_pattern_num; j++)
//						{
//							generateIntrlvPattern(intrlv_pattern_2D_time_partial[i][j], ModType);
//						}
//					}
//				}
//			}*/
//
//			//cout << frame << endl;
//			int u_input[16] = { 0,     1     ,0,     1,
//									0,     1,     0,     0,
//									0 ,    1,     0,     0,
//									1 ,    0,     0,     0 };
//			for (int i = 0; i < N_Info; i++)
//			{
//				a_data[i] = (unsigned char)(distrib(gen) % 2);
//				//		a_data[i] = u_input[i];
//			}
//
//			//cout << endl;
//			if (CRCL) tx_append_crc(a_data, N_Info, CRCL);
//			for (int i = 0; i < Nv; i++)
//				for (int j = 0; j < Nh; j++)
//					u[i][j] = 0;
//			//omp_set_num_threads(4);
////#pragma omp parallel for
//			for (int i = 0; i < Kv; i++)
//			{
//				for (int j = 0; j < Kh; j++)
//				{
//					u[Av[i]][Ah[j]] = (int)a_data[i * Kh + j];
//					//cout << (int)u[Av[i]][Ah[j]] << " \n";
//					//cout << Av[i] << "\n";
//				}
//				//cout << endl; 
//			}
//			// test:original info
//			//cout << "original info u:" << endl;
//			////for (int j = 0; j < Nh; j++) { cout << u[Av[0]][j] << " "; }
//			//for (int i = 0; i < Nv; i++)
//			//{
//			//	cout << endl;
//			//	for (int j = 0; j < Nh; j++)
//			//		cout << u[i][j];
//			//}
//			//cout << endl;
//			// Polar Encoding
////#pragma omp parallel for
//
//			// test: original info
//			/*cout << "original info u:" << endl;
//			display_array(u, Nv, Nh);*/
//			// test: original info
//
//			for (int i = 0; i < Nv; i++)
//				PolarEncode_xor(midx[i], u[i], Nh);
//			// test:horizontal encoding
//			//cout << "u after horizontal encoding:" << endl;
//			////for (int j = 0; j < Nh; j++) { cout << midx[Av[0]][j] << " "; }
//			//for (int i = 0; i < Nv; i++)
//			//{
//			//	cout << endl;
//			//	for (int j = 0; j < Nh; j++)
//			//		cout << midx[i][j];
//			//}
//			//cout << endl;
//			//cout << endl;
//			// test:horizontal encoding
//			// 
//			// column-wise interleaving
//
//			// test: row encoding
//			/*cout << "u after horizontal encoding:" << endl;
//			display_array(midx, Nv, Nh);*/
//			// test: row encoding
//			if (flag_interleave_inner)
//			{
//				int* tmp_to_interleave = new int[Nv];
//				for (int i = 0; i < Nh; i++)
//				{
//					for (int j = 0; j < Nv; j++) { tmp_to_interleave[j] = midx[j][i]; }
//					randIntrlvPartial(intrlv_pattern_2D_inner[i], tmp_to_interleave, Av, Nv, Kv);
//					for (int j = 0; j < Nv; j++) { midx[j][i] = tmp_to_interleave[j]; }
//				}
//				delete[] tmp_to_interleave;
//			}
//
//			// test: innner interleaving
//			/*cout << "u after inner interleaving:" << endl;
//			display_array(midx, Nv, Nh);*/
//			// test: innner interleaving
////#pragma omp parallel for
//			for (int i = 0; i < Nh; i++)
//			{
//				int* tmpxv = new int[Nv]; int* tmpuv = new int[Nv];
//				for (int j = 0; j < Nv; j++)
//					tmpuv[j] = midx[j][i];
//				PolarEncode_xor(tmpxv, tmpuv, Nv);
//				for (int j = 0; j < Nv; j++)
//				{
//					x[j][i] = tmpxv[j];
//					x_encoded_copy[j][i] = x[j][i];
//				}
//				delete[] tmpxv; delete[] tmpuv;
//			}
//
//			if (atten == 2)
//			{
//				// modulation
//				// #pragma omp parallel for
//				// Interleaving
//				if (flag_interleave)
//				{
//					if (flag_interleave_all)  // interleave temporally
//					{
//						for (int i = 0; i < Nv; i++)
//						{
//							randIntrlv(intrlv_pattern_2D_time[i], x[i], Nh);
//						}
//					}
//					//// test: to be interleaved
//					//cout << "x before interleaving:" << endl;
//					//// print x line by line, add spaces between each 4 terms
//					//for (int i = 0; i < Nv; i++)
//					//{
//					//	for (int j = 0; j < Nh; j++)
//					//	{
//					//		if (!(j % 4)) cout << " ";
//					//		cout << x[i][j] << " ";
//					//	}
//					//	cout << endl;
//					//}
//
//					if (flag_interleave_partial)
//					{
//						for (int i = 0; i < Nv; i++)
//						{
//							for (int j = 0; j < intrlv_pattern_num; j++)
//							{
//								randIntrlv(intrlv_pattern_2D_time_partial[i][j], x[i] + j * ModType, ModType);
//							}
//						}
//					}
//
//					//// test: to be interleaved
//					//cout << "x after interleaving:" << endl;
//					//// print x line by line, add spaces between each 4 terms
//					//for (int i = 0; i < Nv; i++)
//					//{
//					//	for (int j = 0; j < Nh; j++)
//					//	{
//					//		if (!(j % 4)) cout << " ";
//					//		cout << x[i][j] << " ";
//					//	}
//					//	cout << endl;
//					//}
//					if ((!flag_interleave_block) && flag_interleave_spatial)
//					{
//						for (int j = 0; j < Nh; j++)
//						{
//							for (int i = 0; i < Nv; i++) { x_tmp[i] = x[i][j]; }
//							randIntrlv(intrlv_pattern_2D[j], x_tmp, Nv);
//							for (int i = 0; i < Nv; i++) { x[i][j] = x_tmp[i]; }
//						}
//					}
//
//
//					if (flag_interleave_block)
//					{
//						randIntrlvBlock(x, Nh, Nv, ModType, interleave_rows);
//					}
//				}
//				if (MOD_DIM == "Spatial")
//				{
//					for (int i = 0; i < Nh; i++)
//					{
//						// Extracting a single column
//						for (int j = 0; j < Nv; j++) { x_tmp[j] = x[j][i]; }
//						modulation(x_tmp, x_mod_tmp, ModType, Nt);
//						for (int j = 0; j < Nsym; j++) { x_mod[j][i] = x_mod_tmp[j]; }
//						// test of transmitted bits
//						/*if (i==0)
//							for (int j = 0; j < Nv; j++)
//							{
//								if (!(j % ModType)) cout << endl;
//								cout << x[j][i] << "  ";
//							}*/
//							// test of transmitted bits
//						/*for (int i = 0; i < mimo_blk_height; i++)
//						{
//							cout << x_mod_tmp[i] << endl;
//						}
//						cout << "111111" << endl;*/
//					}
//					// AWGN Channel (received y_mod )
//					H = generateChannelMatrix(Nr, Nt);
//					for (int i = 0; i < Nh; i++)
//					{
//						// Extracting a single column
//						for (int j = 0; j < Nsym; j++) { x_mod_tmp[j] = x_mod[j][i]; }
//						AWGN(Nr, Nt, x_mod_tmp, y_mod_tmp, H, sigma_sq);
//						for (int j = 0; j < Nr; j++) { y_mod[j][i] = y_mod_tmp[j]; }
//					}
//					// Convert H and y to real domain
//					convertHToReal(H, H_real);
//					// TEST:H_real
//					/*cout << "H_real";
//					for (int i = 0; i < 2 * Nr; i++)
//					{
//						cout << endl;
//						for (int j = 0; j < 2 * Nt; j++)
//							printf("%.10f  ", H_real(i,j));
//					}*/
//					// TEST:H_real
//					for (int i = 0; i < Nh; i++)
//					{
//						for (int j = 0; j < Nr; j++) { y_mod_tmp[j] = y_mod[j][i]; }
//						convertYToReal(Nr, y_mod_tmp, y_mod_real_tmp);
//						for (int j = 0; j < 2 * Nr; j++) { y_mod_real[j][i] = y_mod_real_tmp[j]; /*cout << y_mod_real_tmp[j] << endl;*/ }
//					}
//					for (int i = 0; i < Nh; i++)
//					{
//						for (int j = 0; j < 2 * Nr; j++) { y_mod_real_temp_omp[i][j] = y_mod_real[j][i]; }
//					}
//					// for KBest Detection
//					Eigen::HouseholderQR<Eigen::MatrixXf> qr_Hreal(H_real);
//					Q = qr_Hreal.householderQ();
//					R = qr_Hreal.matrixQR().triangularView<Eigen::Upper>();
//					//omp_set_num_threads(8);
//
//#pragma omp parallel for private(z)
//					for (int j = 0; j < Nh; j++)
//					{
//						MatrixXf received_signal_omp(2 * Nr, 1);
//						int k_kbest_extend = pow(2, ModType / 2) * K_KBEST;
//						float** paths_kbest = new float* [k_kbest_extend];
//						for (int i = 0; i < k_kbest_extend; i++) paths_kbest[i] = new float[2 * Nt];
//						float* ED = new float[K_KBEST];
//
//						if (MIMODET == "MMSE") { MMSEdetection(H_real, Identity, received_signal, mmse_matrix, output_signal, HTH, y_mod_real_temp_omp[j], y_mmse_temp_omp[j], sigma_sq); }
//						// ִ��һЩ����
//
//						for (int i = 0; i < 2 * Nr; i++) { received_signal_omp(i, 0) = y_mod_real_temp_omp[j][i]; }
//						z = Q.transpose() * received_signal_omp;
//						if (MIMODET == "KBEST")
//						{
//							KBestDetection(R, z, H_real, y_mod_real_temp_omp[j], ED, y_mmse_temp_omp[j], paths_kbest, K_KBEST, ModType, Nt, Nr);
//							llr_kbest_bit(Nt, y_mmse_temp_omp[j], oriLLR_kbest[j], sigma_sq, ModType, K_KBEST, paths_kbest, ED);
//						}
//						for (int i_path = 0; i_path < k_kbest_extend; i_path++) { delete[] paths_kbest[i_path]; }
//
//						delete[] paths_kbest;
//						delete[] ED;
//					}
//
//					for (int i = 0; i < Nh; i++)
//					{
//						for (int j = 0; j < 2 * Nt; j++) { y_mmse[j][i] = y_mmse_temp_omp[i][j]; }
//					}
//					// sym llr to bitwise llr
//					// Detection is done for every column, that is, for every Nv bits
//					// Interleaving is not currently considered
//					for (int i = 0; i < Nh; i++)
//					{
//						for (int j = 0; j < 2 * Nt; j++) { y_mmse_tmp[j] = y_mmse[j][i]; }
//						if (MIMODET == "MMSE") { llr_mmse_bit(Nt, y_mmse_tmp, oriLLR_tmp, sigma_sq, ModType); }
//						if (MIMODET == "KBEST") { for (int j = 0; j < len_ori_llr_tmp; j++) { oriLLR_tmp[j] = oriLLR_kbest[i][j]; } }
//						if (flag_interleave)
//						{
//							randDeintrlv(intrlv_pattern, oriLLR_tmp, Nv);
//						}
//						bool flagIsMax = false;
//						for (int j = 0; j < Nv; j++)
//						{
//							if (!(j % 2))
//							{
//								if (abs(x_mmse_out[j][i]) > abs(x_mmse_out[j][i + 1])) { flagIsMax = true; }
//								else { flagIsMax = false; }
//							}
//							else
//							{
//								if (abs(x_mmse_out[j][i]) > abs(x_mmse_out[j][i - 1])) { flagIsMax = true; }
//								else { flagIsMax = false; }
//							}
//							oriLLRh[j][i] = oriLLR_tmp[j];
//							upLLR[j][i] = oriLLR_tmp[j];
//							// TEST�� det error
//							x_mmse_out[j][i] = 0;
//							if (oriLLR_tmp[j] < 0) { x_mmse_out[j][i] = 1; }
//							if (x_mmse_out[j][i] != x[j][i]) {
//								n_det_err++;
//								// error at max-reliable bit. Warning: Only fit for 16QAM modulation
//								if (flagIsMax) { n_det_err_mrb++; }
//							}
//							total_bit_cnt++;
//						}
//					}
//					// test : detection TIME
//					/*det_end = std::chrono::system_clock::now();
//					auto duration_detection = std::chrono::duration_cast<std::chrono::microseconds>(det_end - start);
//					file_det << duration_detection.count() << endl;*/
//				}
//				// Temporal modulation
//				else if (MOD_DIM == "Temporal")
//				{
//					//modulation
//					//cout << "START TEST" << endl;
//					for (int i = 0; i < Nv; i++) // every row
//					{
//						modulation(x[i], x_mod[i], ModType, mimo_blk_length);
//					}
//					/*for (int i = 0; i < mimo_blk_height; i++)
//					{
//						for (int j = 0; j < mimo_blk_length; j++)
//							cout << x_mod[i][j] << endl;
//					}7
//					cout << "111111" << endl;*/
//					//AWGN channel (received y_mod)
//					//// test: modulated symbols
//					//cout << endl;
//					//cout << "modulated symbols" << endl;
//					//for (int i = 0; i < Nt; i++) cout << x_mod[i][0] << endl;
//					//cout << "modulated symbols" << endl;
//					//// test: modulated symbols
//					H = generateChannelMatrix(Nr, Nt);
//					for (int j = 0; j < Nsym; j++)
//					{
//						// Extracting a single column
//						for (int i = 0; i < Nv; i++) { x_mod_tmp[i] = x_mod[i][j]; }
//
//						AWGN(Nr, Nt, x_mod_tmp, y_mod_tmp, H, sigma_sq);
//						for (int i = 0; i < Nr; i++) { y_mod[i][j] = y_mod_tmp[i]; }
//					}
//					convertHToReal(H, H_real);
//
//					for (int j = 0; j < Nsym; j++)
//					{
//						for (int i = 0; i < Nr; i++) { y_mod_tmp[i] = y_mod[i][j]; }
//						convertYToReal(Nr, y_mod_tmp, y_mod_real_tmp);
//						for (int i = 0; i < 2 * Nr; i++) { y_mod_real[i][j] = y_mod_real_tmp[i]; }
//					}
//					//MMSE detection done in paralell
//					for (int i = 0; i < Nsym; i++)
//					{
//						for (int j = 0; j < 2 * Nr; j++) { y_mod_real_temp_omp[i][j] = y_mod_real[j][i]; }
//						// test : y_mod_real_temp_omp[0]
//						/*if (i == 0)
//						{
//							cout << "y" << endl;
//							for (int j = 0; j < 2 * Nr; j++) cout << y_mod_real_temp_omp[0][j] << endl;
//							cout << "y" << endl;
//						}*/
//					}
//					Eigen::HouseholderQR<Eigen::MatrixXf> qr_Hreal(H_real);
//					Q = qr_Hreal.householderQ();
//					R = qr_Hreal.matrixQR().triangularView<Eigen::Upper>();
//					// for EP detection
//					// copy H_real
//					for (int i = 0; i < Nr2; i++)
//					{
//						for (int j = 0; j < Nt2; j++)
//						{
//							Hreal_array[i * Nt2 + j] = H_real(i, j);
//						}
//					}
//					// calculate HTH
//					cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, Nt2, Nt2, Nr2, 1.0, Hreal_array, Nt2, Hreal_array, Nt2, 0.0, HTH_array, Nt2);
//					//cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, Nt2, 1, Nr2, 1.0, ch_mtx_struct.H_real[cnt], Nt2, data_struct.rx_buffer_real[cnt], 1, 0.0, ch_mtx_struct.hty, 1);
//
//					//omp_set_num_threads(8);
//					//cout << 1 << endl;
//					//omp_set_num_threads(4);
//					/*auto start = std::chrono::system_clock::now();*/
//					/*MatrixXf H_real_trans(2*Nt, 2*Nr);
//					for (int i = 0; i < 2 * Nt; i++)
//						for (int j = 0; j < 2 * Nr; j++)
//							H_real_trans(i, j) = H_real(j, i);
//					HTHIinv = (H_real_trans * H_real + Identity).inverse();*/
//# pragma omp parallel for private(z)
//					for (int j = 0; j < Nsym; j++)
//
//						//for (int x=0; x<Nsym; x++)
//					{
//						// int j = x % Nsym;
//						MatrixXf received_signal_omp(2 * Nr, 1);
//						int k_kbest_extend = pow(2, ModType / 2) * K_KBEST;
//						float** paths_kbest = new float* [k_kbest_extend];
//						for (int i = 0; i < k_kbest_extend; i++) { paths_kbest[i] = new float[2 * Nt]; }
//						float* ED = new float[K_KBEST];
//						float* HTY_array_tmp = new float[2 * Nt];
//						// Declared arrays for EP
//						int flag = 1;
//						int temp = 0;
//						//float* extr_mean = new float[2 * Nr];
//						//float* extr_var = new float[2 * Nr];
//						double* vec_alpha_in = new double[2 * Nr];
//						double* vec_gamma_in = new double[2 * Nr];
//						double* p_symbol_rearrange = new double[Nt2 * pow(2, ModType / 2)];
//						for (int i = 0; i < Nt2 * pow(2, ModType / 2); i++) { p_symbol_rearrange[i] = 1 / (pow(2, ModType / 2)); }
//						for (int i = 0; i < Nr2; i++)
//						{
//							//extr_mean[i] = 0;
//							//extr_var[i] = 0;
//							vec_alpha_in[i] = 2;
//							vec_gamma_in[i] = 0;
//						}
//
//						if (MIMODET == "MMSE")
//						{
//							MMSEdetection(H_real, Identity, received_signal, mmse_matrix, output_signal, HTH, y_mod_real_temp_omp[j], y_mmse_temp_omp[j], sigma_sq);
//							//MMSEdetection_noInverse(H_real, Identity, received_signal, mmse_matrix, output_signal, HTH, HTHIinv, y_mod_real_temp_omp[j], y_mmse_temp_omp[j], sigma_sq);
//						}
//						if (MIMODET == "KBEST")
//						{
//							for (int i = 0; i < 2 * Nr; i++) { received_signal_omp(i, 0) = y_mod_real_temp_omp[j][i]; }
//							z = Q.transpose() * received_signal_omp;
//							KBestDetection(R, z, H_real, y_mod_real_temp_omp[j], ED, y_mmse_temp_omp[j], paths_kbest, K_KBEST, ModType, Nt, Nr);
//							llr_kbest_bit(Nt, y_mmse_temp_omp[j], oriLLR_kbest[j], sigma_sq, ModType, K_KBEST, paths_kbest, ED);
//						}
//						if (MIMODET == "EP")
//						{
//							// test: save H_real, recv_signal
//							/*save_mat_txt(y_mod_real_temp_omp[j], Nr2, "RxSig0128.txt");
//							save_mat_txt(H_real, "Hreal0128.txt");*/
//							// test: save H_real, recv_signal
//							cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, Nt2, 1, Nr2, 1.0, Hreal_array, Nt2, y_mod_real_temp_omp[j], 1, 0.0, HTY_array_tmp, 1);
//							//test:HTH_array
//							/*for (int i = 0; i < Nt2; i++)
//								for (int s = 0; s < Nt2; s++)
//								{
//									if (s == 0) std::cout << endl;
//									std::cout << setw(10) << HTH_array[i * Nt2 + s];
//								}
//							cout << endl;*/
//							// TEST:HTY
//							/*for (int i = 0; i < Nr2; i++)
//								received_signal_omp(i, 0) = y_mod_real_temp_omp[j][i];
//							MatrixXf HTransY = H_real.transpose() * received_signal_omp;
//							cout << "Eigen" << endl;
//							for (int i = 0; i < Nt2; i++) cout << setw(10)<<HTransY(i, 0);
//							cout << "mkl" << endl;
//							for (int i = 0; i < Nt2; i++) cout << HTY_array[i];*/
//							// TEST:HTY
//							EPD_detection_float(flag, Nt, Nr, iter_MIMO, ModType, 1.0, sigma_sq, HTH_array, HTY_array_tmp, temp, ep_extr_mean_omp[j], ep_extr_var_omp[j], vec_alpha_in, vec_gamma_in, p_symbol_rearrange);
//						}
//						//test EP detection
//						/*cout << "TEST:EP DETECTION" << endl;
//						for (int i = 0; i < Nr2; i++) std::cout << extr_mean[i]<< endl;
//						cout << "TEST:EP DETECTION" << endl;*/
//						//test EP detection
//						for (int i_path = 0; i_path < k_kbest_extend; i_path++) { delete[] paths_kbest[i_path]; }
//						delete[] paths_kbest;
//						delete[] ED;
//						//delete[] extr_mean;
//						//delete[] extr_var;
//						delete[] vec_alpha_in;
//						delete[] vec_gamma_in;
//						delete[] p_symbol_rearrange;
//						delete[] HTY_array_tmp;
//					}
//					// test: detection time
//					/*auto detection_over = std::chrono::system_clock::now();
//					auto duration_detection = std::chrono::duration_cast<std::chrono::microseconds>(detection_over - start);
//					cout << endl;
//					cout <<"Detection duration" << duration_detection.count() << "us" << endl;*/
//					// test: detection time
//					for (int i = 0; i < Nsym; i++)
//					{
//						for (int j = 0; j < 2 * Nt; j++) { y_mmse[j][i] = y_mmse_temp_omp[i][j]; }
//					}
//					//To bit domain output
//					for (int j = 0; j < Nsym; j++)
//					{
//						for (int i = 0; i < 2 * Nt; i++) { y_mmse_tmp[i] = y_mmse[i][j]; }
//						// With Spatial Modulation, Nt = Nv
//						if (MIMODET == "MMSE") { llr_mmse_bit(Nt, y_mmse_tmp, oriLLR_tmp, sigma_sq, ModType); }
//						if (MIMODET == "EP") { llr_ep_bit(Nt, ep_extr_mean_omp[j], ep_extr_var_omp[j], oriLLR_tmp, ModType); }
//						if (MIMODET == "KBEST") { for (int i = 0; i < len_ori_llr_tmp; i++) { oriLLR_tmp[i] = oriLLR_kbest[j][i]; } }
//						int i_bit = 0;
//						int bit_col = 0;
//						for (int j_inner = 0; j_inner < Nv; j_inner++)
//						{
//							for (int i_bit_col = 0; i_bit_col < ModType; i_bit_col++)
//							{
//								i_bit = ModType * j_inner + i_bit_col;
//								bit_col = j * ModType + i_bit_col;
//								oriLLRh[j_inner][bit_col] = oriLLR_tmp[i_bit];
//								upLLR[j_inner][bit_col] = oriLLR_tmp[i_bit];
//							}
//
//						}
//					}
//					// test: LLR transfer time
//					/*auto llr_to_bit_over = std::chrono::system_clock::now();
//					auto duration_llr = std::chrono::duration_cast<std::chrono::microseconds>(llr_to_bit_over - detection_over);
//					cout << endl;
//					cout<< "LLR transfer duration" << duration_llr.count() << "us" << endl;*/
//					// test: LLR transfer time
//					if (flag_interleave)
//					{
//						if ((!flag_interleave_block) && flag_interleave_spatial)
//						{
//							for (int j = 0; j < Nh; j++)
//							{
//								for (int i = 0; i < Nv; i++) { oriLLR_tmp[i] = oriLLRh[i][j]; }
//								randDeintrlv(intrlv_pattern_2D[j], oriLLR_tmp, Nv);
//								for (int i = 0; i < Nv; i++)
//								{
//									oriLLRh[i][j] = oriLLR_tmp[i];
//									upLLR[i][j] = oriLLR_tmp[i];
//								}
//
//							}
//						}
//
//						if (flag_interleave_block)
//						{
//							randDeIntrlvBlock(oriLLRh, Nh, Nv, ModType, interleave_rows);
//							randDeIntrlvBlock(upLLR, Nh, Nv, ModType, interleave_rows);
//						}
//
//						if (flag_interleave_all) // Deinterleave in time domain
//						{
//							for (int j = 0; j < Nv; j++)
//							{
//								randDeintrlv(intrlv_pattern_2D_time[j], oriLLRh[j], Nh);
//								randDeintrlv(intrlv_pattern_2D_time[j], upLLR[j], Nh);
//							}
//						}
//
//						if (flag_interleave_partial)
//						{
//							for (int j = 0; j < Nv; j++)
//							{
//								for (int k = 0; k < intrlv_pattern_num; k++)
//								{
//									randDeintrlv(intrlv_pattern_2D_time_partial[j][k], oriLLRh[j] + k * ModType, ModType);
//									randDeintrlv(intrlv_pattern_2D_time_partial[j][k], upLLR[j] + k * ModType, ModType);
//								}
//							}
//						}
//
//					}
//
//					for (int i = 0; i < Nv; i++)
//					{
//						for (int j = 0; j < Nh; j++)
//						{
//							total_bit_cnt++;
//							int bit_dec = (oriLLRh[i][j] < 0);
//							if (bit_dec != x[i][j]) n_det_err++;
//						}
//					}
//
//				}
//				// test: detection errors
//
//				// test: detection errors
//			}
//			if (atten == 1)
//			{
//				for (int i = 0; i < Nv; i++)
//					for (int j = 0; j < Nh; j++)
//					{
//						// 0-1; 1-(-1)
//						y[i][j] = (1 - 2 * x[i][j]) + (float)sqrt(sigma_sq) * (float)sampleNormal();
//						oriLLRh[i][j] = 2 * y[i][j] / sigma_sq;
//						//oriLLRv[j][i] = oriLLRh[i][j];
//						upLLR[i][j] = oriLLRh[i][j];
//					}
//			}
//			//Polar Decoding
//			for (int iter = 1; iter <= softiter; iter++)
//			{
//				for (int i = 0; i < Nv; i++) { n_tried_patterns_h[i] = 0; }
//				for (int i = 0; i < Nh; i++) { n_tried_patterns_v[i] = 0; }
//				//cout << '1' << endl;
//				//vertical decoding
//		   // test: decoding time
//			//auto decoding_start = std::chrono::system_clock::now();
//
//			// test: decoding time
//#pragma omp parallel for
//				for (int decodingi = 0; decodingi < Nh; decodingi++) // 291.70kb
//				{
//
//					//for SCL Decoding
//					int Counth = 0, Countinfoh = 0;
//					float Qsumunh = 0;
//					int** sum = new int* [L];
//					//for (int i = 0; i < L; i++)  sum[i] = new  int[Nmax];
//					int** u1 = new int* [L];
//					for (int i = 0; i < L; i++)  u1[i] = new _MM_ALIGN16 int[Nmax];
//					for (int i = 0; i < L; i++)  sum[i] = new _MM_ALIGN16 int[Nmax];
//					//for (int i = 0; i < L; i++)  u1[i] = new  int[Nmax];
//
//					float* LLR_in = new float[L];
//					float* W = new float[2 * L];
//					int* SCL_index = new int[L];
//					int* better = new int[L];
//					int* worse = new int[L];
//					int* indexpm = new int[L];
//					float* PM = new float[L];
//					float** LLR = new float* [L];
//					for (int i = 0; i < L; i++)  LLR[i] = new _MM_ALIGN16 float[2 * Nmax + 2];
//					// for (int i = 0; i < L; i++)  LLR[i] = new float[2 * Nmax + 2];
//					////////////////////////////////
//					for (int j = 0; j < Nv; j++)
//						for (int k = 0; k < L; k++)
//							LLR[k][j] = upLLR[j][decodingi];
//					//cout << LLR[0][j] << " ";
//
//				//cout << endl;
//					int index = 0;
//					/*			if (L == 1)
//								{
//									//T_start = clock();
//									Count = 0;
//									decode(*LLR, *sum, (int)(log(Nv) / log(2)), A_Acv, flag, Nv, Kv);
//									//PolarEncode_xor(*u1, *sum, N);
//									//T_finish = clock();
//								}
//								else
//								{*/
//					for (int k = 0; k < L; k++) { PM[k] = 0; indexpm[k] = k; }
//					//T_start = clock();
//					decode_list(LLR, sum, (int)(log(Nv) / log(2)), A_Acv, 0, Nv, Kv, PM, L, (int)(log(L) / log(2)), LLR_in, W, SCL_index, better, worse, &Counth, &Countinfoh, &Qsumunh);
//					//T_finish = clock();
//					//Decoding_Duration += (T_finish - T_start);
//					for (int i = 0; i < L - 1; i++)
//						for (int j = i + 1; j < L; j++)
//							if (PM[indexpm[i]] < PM[indexpm[j]])
//							{
//								int ttt = indexpm[i];
//								indexpm[i] = indexpm[j];
//								indexpm[j] = ttt;
//							}
//					//PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), N);
//					index = indexpm[0];
//					//}
//					//calc distance s()
//					if (adaptive_beta)
//					{
//						if (softmethod == "Pyndiah") Pyndiahv_adaptivebeta(LLR, sum, oriLLRh, indexpm, iter * 2 - 2, decodingi, upLLR);
//						if (softmethod == "MITSO") MITSOv(LLR, PM, sum, oriLLRh, indexpm, iter * 2 - 2, decodingi, upLLR, Qsumunh);
//					}
//					// iter * 2 - 2
//					else
//					{
//						//cout << iter << endl;
//						if (softmethod == "Pyndiah") {
//							if (!bit_flip_pyndiah | (!vertical_ml)) { Pyndiahv(LLR, sum, oriLLRh, indexpm, iter, decodingi, upLLR); }
//							else { Pyndiahv_bitflip(LLR, sum, oriLLRh, indexpm, iter, decodingi, upLLR, sum_sdis_array[decodingi], tried_codewords_v[decodingi], GGv, Acv, n_tried_patterns_v); }
//						}
//
//						if (softmethod == "MITSO") MITSOv(LLR, PM, sum, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunh);
//					}
//					// test: oriLLR and upLLR
//					/*cout << endl;
//					for (int i = 0; i < Nv; i++)
//						for (int j = 0; j < Nh; j++)
//							cout << oriLLRh[i][j] - upLLR[i][j] << "   ";
//					cout << endl;*/
//					/*cout << oriLLRh[0][0] << endl;
//					cout << upLLR[0][0] << endl;*/
//					// test: oriLLR and upLLR
//					// half iteration Enabled
//					if (half_iter && iter == softiter)
//					{
//						if (iter == softiter) PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), Nh);
//						for (int j = 0; j < Kv; j++)
//							umid[Av[j]][decodingi] = u1[index][Av[j]];
//						for (int i = 0; i < L; i++)
//						{
//							int curr_index = indexpm[i];
//							for (int j = 0; j < Nv; j++)
//							{
//								sub_encoded_codeword_list_v[decodingi][i * Nv + j] = sum[curr_index][j];
//							}
//							PolarEncode_xor(sub_codeword_list_v[decodingi] + i * Nv, sum[curr_index], Nv);
//							edistance_list_v[decodingi][i] = PM[curr_index];
//						}
//						// check wheher an item of sum equals to x_encoded_copy[:][decodingi]
//						int best_idx = indexpm[0];
//						/*cout << endl;
//						for (int k = 0; k < Nh; k++)
//						{
//							cout << sum[best_idx][k] << "  " << x_encoded_copy[k][decodingi] << endl;
//						}*/
//
//
//						int correct_codeword_idx = 0;
//						int correct_idx_pre = -1;
//						for (int i = 0; i < L; i++)
//						{
//							// cout << 1 << endl;
//							int curr_index = indexpm[i];
//							int found_correct = 0;
//							for (int j = 0; j < Nv; j++)
//							{
//								if (sum[curr_index][j] != x_encoded_copy[j][decodingi]) { break; }
//								if (j == Nv - 1) { correct_codeword_idx = curr_index; found_correct = 1; }
//							}
//
//							// display x_encoded_copy
//
//							if (found_correct)
//							{
//								n_err_codeword_v[decodingi] += i;
//								n_err_codeword_v_diffpos[decodingi][i] += 1;
//								correct_idx_oneframe[decodingi] = i;
//								break;
//							}
//							if (i == L - 1) {
//								n_err_codeword_v[decodingi] += L;
//								n_err_codeword_v_diffpos[decodingi][L] += 1;
//								correct_idx_oneframe[decodingi] = L;
//								break;
//							}
//						}
//					}
//					// copy to sub_encoded_codeword_list_v
//					/*for (int i = 0; i < L; i++)
//					{
//						for (int j = 0; j < Nv; j++)
//						{
//							sub_encoded_codeword_list_v[i][i * Nv + j] = sum[indexpm[i]][j];
//						}
//					}*/
//					// half iteration
//					////////delete SCL variables
//					for (int i = 0; i < L; i++)
//						delete[]LLR[i];
//					delete[] LLR;
//					for (int i = 0; i < L; i++)
//						delete[]sum[i];
//					for (int i = 0; i < L; i++)
//						delete[]u1[i];
//					delete[] u1;
//					delete[] sum;
//					delete[] PM;
//					delete[] LLR_in;
//					delete[] W;
//					delete[] SCL_index;
//					delete[] better;
//					delete[] worse;
//					delete[] indexpm;
//
//				}
//				//cout << "Horizon Start:\n";
//				//horizontal decoding
//				// #
//				if (weirditer) { continue; }
//
//				// test: LLR after horizontal decoding
//				/*cout << "After horizontal decoding" << endl;
//				for (int i = 0; i < Nv; i++)
//				{
//					for (int j = 0; j < Nh; j++)
//					{
//						cout << (oriLLRh[i][j] < 0) << " ";
//					}
//					cout << endl;
//				}*/
//				// test: LLR after horizontal decoding
//				// inner deinterleaving
//				// randDeIntrlvPartial(intrlv_pattern_2D_inner[0],);
//				if (flag_interleave_inner)
//				{
//					float* tmp_llr_to_interleave = new float[Nv];
//					for (int i = 0; i < Nh; i++)
//					{
//						for (int j = 0; j < Nv; j++) { tmp_llr_to_interleave[j] = oriLLRh[j][i]; }
//						randDeIntrlvPartial(intrlv_pattern_2D_inner[i], tmp_llr_to_interleave, Av, Nv, Kv);
//						for (int j = 0; j < Nv; j++) { oriLLRh[j][i] = tmp_llr_to_interleave[j]; }
//					}
//					for (int i = 0; i < Nh; i++)
//					{
//						for (int j = 0; j < Nv; j++) { tmp_llr_to_interleave[j] = upLLR[j][i]; }
//						randDeIntrlvPartial(intrlv_pattern_2D_inner[i], tmp_llr_to_interleave, Av, Nv, Kv);
//						for (int j = 0; j < Nv; j++) { upLLR[j][i] = tmp_llr_to_interleave[j]; }
//					}
//					delete[] tmp_llr_to_interleave;
//				}
//
//				if ((iter == softiter) && half_iter)
//				{
//					int max_correct_idx = -1;
//					for (int i = 0; i < Nh; i++)
//					{
//						if (correct_idx_oneframe[i] > max_correct_idx)
//						{
//							max_correct_idx = correct_idx_oneframe[i];
//						}
//					}
//					max_L_cnt[max_correct_idx]++;
//				}
//				/*auto decoding_over = std::chrono::system_clock::now();
//				auto duration_decoding = std::chrono::duration_cast<std::chrono::microseconds>(decoding_over - decoding_start);
//				cout << endl;
//				cout << "Decoding duration" << duration_decoding.count() << "us" << endl;*/
//
//#pragma omp parallel for
//				for (int decodingi = 0; decodingi < Nv; decodingi++)
//				{
//					// 293.70kb
//					// cout << "KKKK" << endl;
//					if (half_iter && (iter == softiter)) { break; }
//					string polarde_now = polarde_array[iter - 1];
//					//for SCL Decoding
//					int Countv = 0, Countinfov = 0;
//					float Qsumunv = 0;
//					int** u1 = new int* [L];
//					for (int i = 0; i < L; i++)  u1[i] = new _MM_ALIGN16 int[Nmax];
//					int** u1_LSD = new int* [2 * L_LSD];
//					for (int i = 0; i < 2 * L_LSD; i++) u1_LSD[i] = new int[Nmax];
//					int** sum = new int* [L];
//					for (int i = 0; i < L; i++)  sum[i] = new _MM_ALIGN16 int[Nmax];
//					int** sum_LSD = new int* [2 * L_LSD];
//					for (int i = 0; i < 2 * L_LSD; i++)  sum_LSD[i] = new int[Nmax];
//					float* LLR_in = new float[L];
//					int* u_decod = new int[Nh];
//					int* n_paths = new int[Nh];
//					float* W = new float[2 * L];
//					int* SCL_index = new int[L];
//					int* better = new int[L];
//					int* worse = new int[L];
//					int* indexpm = new int[L_LSD];
//					float* PM = new float[L];
//					float** LLR = new float* [L];
//					for (int i = 0; i < L; i++)  LLR[i] = new _MM_ALIGN16 float[2 * Nmax + 2];
//
//					for (int i = 0; i < L_LSD; i++) { indexpm[i] = i; }
//					////////////////////////////////
//					for (int j = 0; j < Nh; j++)
//						for (int k = 0; k < L; k++)
//							LLR[k][j] = upLLR[decodingi][j];
//					//	cout << LLR[0][j] << " ";
//					calc_npaths(A_Ach, Nh, n_paths, L_LSD);
//					//	cout << endl;
//					int index = 0;
//
//					/*			if (L == 1)
//								{
//									//T_start = clock();
//									Count = 0;
//									decode(*LLR, *sum, (int)(log(Nh) / log(2)), A_Ach, flag, Nh, Kh);
//									if (iter == softiter) PolarEncode_xor(*u1, *sum, Nh);
//									//T_finish = clock();
//								}
//								else
//								{*/
//								//for (int k = 0; k < L; k++)  PM[k] = 0, indexpm[k] = k; 
//								////T_start = clock();
//								//decode_list(LLR, sum, (int)(log(Nh) / log(2)), A_Ach, 0, Nh, Kh, PM, L, (int)(log(L) / log(2)), LLR_in, W, SCL_index, better, worse, &Countv, &Countinfov, &Qsumunv);
//								////T_finish = clock();
//								////Decoding_Duration += (T_finish - T_start);
//								//for (int i = 0; i < L - 1; i++)
//								//	for (int j = i + 1; j < L; j++)
//								//		if (PM[indexpm[i]] < PM[indexpm[j]])
//								//		{
//								//			int ttt = indexpm[i];
//								//			indexpm[i] = indexpm[j];
//								//			indexpm[j] = ttt;
//								//		}
//
//
//							/*
//							* Encoding: Info -> Horizontal Encoding -> A -> set to 0 -> B -> Horizontal Encoding -> Vertical Encoding
//							* Decoding: Received Info -> Horizontal Decoding -> Vertical Decoding -> Horizontal Decoding
//							*/
//					if (polarde_now == "SLD")
//					{
//						LSD_decode_outall(u_decod, u1_LSD, LLR[0], GGh, A_Ach, Nh, L_LSD, n_paths);
//						// if (iter == softiter) PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), Nh);
//						// for soft information exchange, inverse coding:
//						// THIS LINE is WRONG! BUT IT SEEMS BETTER THAN THE RIGHT VERSION
//						for (int i = 0; i < L_LSD; i++) { PolarEncode_xor(sum_LSD[i], u1_LSD[i], Nh); }
//						// for (int i = 0; i < L_LSD; i++) { PolarEncode_xor(sum_LSD[i], u1_LSD[i], Nh); }
//						// if (iter == softiter) PolarEncode_xor(u1[0], sum_LSD[0], Nh);
//						index = indexpm[0];
//						// test: after horizontal decoding
//						/*if (decodingi == Av[0])
//						{
//							cout << "Direct output of horizontal decoding:" << endl;
//							for (int i = 0; i < Nh; i++) { cout << u1_LSD[0][i] << " "; }
//							cout << endl;
//							cout << "horizontal decoder out, encode once:" << endl;
//							for (int i = 0; i < Nh; i++) { cout << sum_LSD[0][i] << " "; }
//							cout << endl;
//						}*/
//						// TEST
//						/*cout << "SCL decoded bits" << endl;
//						for (int i = 0; i < Nh; i++) { cout << sum[indexpm[0]][i] << "  "; }
//						cout << endl;
//						cout << "LSD decoded bits" << endl;
//						for (int i = 0; i < Nh; i++) { cout << u1[0][i] << "  "; }
//						cout << endl;*/
//						// TEST
//						//}
//						//calc distance s()
//						// iter*2-1
//						if (iter < softiter)
//						{
//							if (adaptive_beta)
//							{
//								if (softmethod == "Pyndiah") Pyndiahh_LSD_adaptivebeta(LLR, sum_LSD, oriLLRh, indexpm, iter, decodingi, upLLR);
//								if (softmethod == "MITSO") MITSOh(LLR, PM, sum_LSD, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunv);
//							}
//							else
//							{
//								// cout << "IN LSD!!!" << endl;
//								if (softmethod == "Pyndiah")
//									Pyndiahh_LSD_getsum(LLR, sum_LSD, oriLLRh, indexpm, iter, decodingi, upLLR, sum_sdis_array[decodingi], GGh, Ach);
//								//if (softmethod == "Pyndiah") Pyndiahh_LSD(LLR, sum_LSD, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR);
//								if (softmethod == "MITSO") MITSOh(LLR, PM, sum_LSD, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunv);
//							}
//							//	getchar();
//						}
//						else
//							for (int j = 0; j < Kh; j++)
//								umid[decodingi][Ah[j]] = u1_LSD[index][Ah[j]];
//					}
//					if (polarde_now == "SCL")
//					{
//						for (int k = 0; k < L; k++)  PM[k] = 0, indexpm[k] = k;
//						//T_start = clock();
//						decode_list(LLR, sum, (int)(log(Nh) / log(2)), A_Ach, 0, Nh, Kh, PM, L, (int)(log(L) / log(2)), LLR_in, W, SCL_index, better, worse, &Countv, &Countinfov, &Qsumunv);
//						//T_finish = clock();
//						//Decoding_Duration += (T_finish - T_start);
//						for (int i = 0; i < L - 1; i++)
//							for (int j = i + 1; j < L; j++)
//								if (PM[indexpm[i]] < PM[indexpm[j]])
//								{
//									int ttt = indexpm[i];
//									indexpm[i] = indexpm[j];
//									indexpm[j] = ttt;
//								}
//						if (iter == softiter) PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), Nh);
//						index = indexpm[0];
//						for (int col = 0; col < Nh; col++) { x_parity_check[decodingi][col] = sum[indexpm[0]][col]; }
//						// show sum[indexpm[0]]
//						/*cout << "sum[indexpm[0]]" << endl;
//						for (int i = 0; i < Nh; i++) { cout << sum[indexpm[0]][i] << "  "; }
//						cout << "sum[indexpm[0]]" << endl;*/
//						// show sum[indexpm[0]]
//						//}
//						//calc distance s()
//						if (iter < softiter)
//						{
//							if (adaptive_beta)
//							{
//								if (softmethod == "Pyndiah") Pyndiahh_adaptivebeta(LLR, sum, oriLLRh, indexpm, iter, decodingi, upLLR);
//								if (softmethod == "MITSO") MITSOh(LLR, PM, sum, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunv);
//								//	getchar();
//							}
//							else
//							{
//								//if (softmethod == "Pyndiah") Pyndiahh(LLR, sum, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR);
//								if (softmethod == "Pyndiah")
//								{
//									if (!bit_flip_pyndiah)  Pyndiahh_getsum(LLR, sum, oriLLRh, indexpm, iter, decodingi, upLLR, sum_sdis_array[decodingi]);
//									else Pyndiahh_bitflip(LLR, sum, oriLLRh, indexpm, iter, decodingi, upLLR, sum_sdis_array[decodingi], tried_codewords_h[decodingi], GGh, Ach, n_tried_patterns_h);
//
//								}
//								if (softmethod == "MITSO") MITSOh(LLR, PM, sum, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunv);
//								//	getchar();
//							}
//						}
//						else
//						{
//							// prev
//							/*for (int j = 0; j < Kh; j++)
//								umid[decodingi][Ah[j]] = u1[index][Ah[j]];*/
//
//							for (int j = 0; j < Nh; j++)
//								umid[decodingi][j] = u1[index][j];
//							// add result of current decoding to subcodeword_list
//							for (int j = 0; j < L; j++)
//							{
//								int curr_index = indexpm[j];
//								PolarEncode_xor(sub_codeword_list[decodingi] + j * Nh, sum[curr_index], Nh);
//								for (int i = 0; i < Nh; i++) { sub_encoded_codeword_list[decodingi][j * Nh + i] = sum[curr_index][i]; }
//								edistance_list[decodingi][j] = PM[curr_index];
//							}
//
//							if (test_cadsl && (!half_iter)) {
//								int best_idx = indexpm[0];
//								/*cout << endl;
//								for (int k = 0; k < Nh; k++)
//								{
//									cout << sum[best_idx][k] << "  " << x_encoded_copy[k][decodingi] << endl;
//								}*/
//
//
//								int correct_codeword_idx = 0;
//								int correct_idx_pre = -1;
//								for (int i = 0; i < L; i++)
//								{
//									// cout << 1 << endl;
//									int curr_index = indexpm[i];
//									int found_correct = 0;
//									for (int j = 0; j < Nv; j++)
//									{
//										if (sum[curr_index][j] != x_encoded_copy[decodingi][j]) { break; }
//										if (j == Nv - 1) { correct_codeword_idx = curr_index; found_correct = 1; }
//									}
//
//									// display x_encoded_copy
//
//									if (found_correct)
//									{
//										n_err_codeword_v[decodingi] += i;
//										n_err_codeword_v_diffpos[decodingi][i] += 1;
//										correct_idx_oneframe[decodingi] = i;
//										break;
//									}
//									if (i == L - 1) {
//										n_err_codeword_v[decodingi] += L;
//										n_err_codeword_v_diffpos[decodingi][L] += 1;
//										correct_idx_oneframe[decodingi] = L;
//										break;
//									}
//								}
//							}
//						}
//
//					}
//
//					////////delete SCL variables
//					for (int i = 0; i < L; i++)
//						delete[]u1[i];
//					delete[] u1;
//					for (int i = 0; i < 2 * L_LSD; i++)
//						delete[]u1_LSD[i];
//					delete[] u1_LSD;
//					for (int i = 0; i < L; i++)
//						delete[]LLR[i];
//					delete[] LLR;
//					for (int i = 0; i < L; i++)
//						delete[]sum[i];
//					delete[] sum;
//					for (int i = 0; i < 2 * L_LSD; i++)
//						delete[]sum_LSD[i];
//					delete[] sum_LSD;
//					delete[] PM;
//					delete[] LLR_in;
//					delete[] W;
//					delete[] SCL_index;
//					delete[] better;
//					delete[] worse;
//					delete[] indexpm;
//					delete[] u_decod;
//					delete[] n_paths;
//				}
//
//				if (iter == softiter && (!half_iter) && test_cadsl)
//				{
//					int max_correct_idx = -1;
//					for (int i = 0; i < Nh; i++)
//					{
//						if (correct_idx_oneframe[i] > max_correct_idx)
//						{
//							max_correct_idx = correct_idx_oneframe[i];
//						}
//					}
//					max_L_cnt[max_correct_idx]++;
//				}
//
//				// list extension after horizontal decoding
//				//if (iter == softiter && CRCL>0)
//				//{
//				//	crc_list_extend(sub_codeword_list, whole_codeword_list, edistance_list, L, Nv, Nh, L_WHOLE_CODEWORD);
//				//	/*for (int i_subcode = 0; i_subcode < Nv; i_subcode++)
//				//	{
//				//		cout << endl << "subcodeword in whole codeword list" << endl;
//				//		for (int i_bit = 0; i_bit < Nh; i_bit++)
//				//			cout << whole_codeword_list[0][i_subcode * Nh + i_bit];
//				//		cout << endl << "subcodeword in decoded result" << endl;
//				//		for (int i_bit = 0; i_bit < Nh; i_bit++)
//				//			cout << umid[i_subcode][i_bit];
//				//	}
//				//	cout << "TEST ENDS" << endl;*/
//				//}
//
//				// test: double stage crc list extend
//				/*if (iter == softiter)
//				{
//					crc_list_extend_vertical_doublestage(sub_codeword_list_v, whole_codeword_list_v, intrlv_pattern_2D_time_partial, edistance_list_v,
//						y_mod_real_temp_omp, Hreal_array, L_SUB_CODEWORD, Nh, Nv, Nr, L_WHOLE_CODEWORD, L_MIMO_BLK, ModType);
//				}*/
//				// test: double stage crc list extend
//
//
//
//
//
//
//				if (flag_interleave_inner)
//				{
//					float* tmp_llr_to_interleave = new float[Nv];
//					for (int i = 0; i < Nh; i++)
//					{
//						for (int j = 0; j < Nv; j++) { tmp_llr_to_interleave[j] = oriLLRh[j][i]; }
//						randIntrlvPartial(intrlv_pattern_2D_inner[i], tmp_llr_to_interleave, Av, Nv, Kv);
//						for (int j = 0; j < Nv; j++) { oriLLRh[j][i] = tmp_llr_to_interleave[j]; }
//					}
//					for (int i = 0; i < Nh; i++)
//					{
//						for (int j = 0; j < Nv; j++) { tmp_llr_to_interleave[j] = upLLR[j][i]; }
//						randIntrlvPartial(intrlv_pattern_2D_inner[i], tmp_llr_to_interleave, Av, Nv, Kv);
//						for (int j = 0; j < Nv; j++) { upLLR[j][i] = tmp_llr_to_interleave[j]; }
//					}
//					delete[] tmp_llr_to_interleave;
//				}
//				for (int i = 0; i < Nv; i++) { n_tried_patterns_h_all[iter - 1] += n_tried_patterns_h[i]; }
//				for (int i = 0; i < Nh; i++) { n_tried_patterns_v_all[iter - 1] += n_tried_patterns_v[i]; }
//			}
//			int calcsum[Nv], calcu1[Nv];
//			int calcsumh[Nh], calcu1h[Nh];
//			// calcsum is from the direct output of the decoder
//			// if output is the result of row decoding (meaning there is another column decoding to be done)
//			if (half_iter)
//			{
//				/*if (CRCL > 0)
//				{
//					crc_list_extend_vertical(sub_codeword_list_v, whole_codeword_list_v, intrlv_pattern_2D_time_partial, edistance_list_v, y_mod_real_temp_omp,
//						Hreal_array, L, Nh, Nv, Nr, L_WHOLE_CODEWORD, ModType);
//				}*/
//
//				// test: sub_codeword_list_v
//				/*cout << "SUB_CODEWORD_LIST_V" << endl;
//				for (int i = 0; i < Nh; i++)
//				{
//					int* encoded = new int[Nv];
//					PolarEncode_xor(encoded, sub_codeword_list_v[i], Nh);
//					for (int j = 0; j < Nv; j++) sub_codeword_list_v[i][j] = encoded[j];
//				}
//				for (int i = 0; i < Nv; i++)
//				{
//					for (int j = 0; j < ModType; j++)
//					{
//						cout << sub_codeword_list_v[j][i] << "  ";
//					}
//					cout << endl;
//				}
//				cout << "SUB_CODEWORD_LIST_V" << endl;*/
//				// test: sub_codeword_list_v
//				/*crc_list_extend_vertical(sub_codeword_list_v, whole_codeword_list_v, intrlv_pattern_2D_time_partial, edistance_list_v, y_mod_real_temp_omp,
//					Hreal_array, L, Nh, Nv, Nr, L_WHOLE_CODEWORD, ModType);*/
//
//					/*for (int i = 0; i < Nv; i++)
//					{
//						for (int j = 0; j < Nh; j++) { calcsumh[j] = whole_codeword_list_v[0][j * Nv + i]; }
//						PolarEncode_xor(calcu1h, calcsumh, Nh);
//						for (int j = 0; j < Nh; j++)
//							uout[i][j] = calcu1h[j];
//					}*/
//					/*for (int i = 0; i < Kv; i++)
//					{
//						PolarEncode_xor(calcu1h, umid[Av[i]], Nh);
//						for (int j = 0; j < Kh; j++)
//							uout[Av[i]][Ah[j]] = calcu1h[Ah[j]];
//					}*/
//				for (int i = 0; i < Nv; i++)
//				{
//					PolarEncode_xor(calcu1h, umid[i], Nh);
//					for (int j = 0; j < Nh; j++)
//						uout[i][j] = calcu1h[j];
//				}
//
//				if (CRCL > 0 && (CRC_CHECK_ENABLE))
//				{
//					flag_pass_crc = false;
//					flag_found_replacement_codeword = false;
//
//
//
//					for (int i = 0; i < Kv; i++)
//					{
//						for (int j = 0; j < Kh; j++)
//							u_crc[i * Kh + j] = uout[Av[i]][Ah[j]];
//					}
//					if (rx_check_crc(u_crc, Kh * Kv, CRCL)) {
//						flag_pass_crc = true;
//					}
//
//
//					// the i th horizontal codeword , j th col
//					// whole_codeword_list_v[j*Nv+i];
//					if (!flag_pass_crc)
//					{
//						//cout << endl << "NOT PASS" << endl;
//						// float* H_real_array;
//						// crc_list_extend_vertical(sub_codeword_list_v, whole_codeword_list_v, intrlv_pattern_2D_time_partial, edistance_list_v, y_mod_real_temp_omp,
//							//Hreal_array, L, Nh, Nv, Nr, L_WHOLE_CODEWORD, ModType);
//						// prev: 
//						/*crc_list_extend_vertical(sub_codeword_list_v, whole_codeword_list_v, intrlv_pattern_2D_time_partial, edistance_list_v, y_mod_real_temp_omp,
//								Hreal_array, L_SUB_CODEWORD, Nh, Nv, Nr, L_WHOLE_CODEWORD, ModType);*/
//						crc_list_extend_vertical_doublestage(sub_codeword_list_v, whole_codeword_list_v, intrlv_pattern_2D_time_partial, edistance_list_v,
//							y_mod_real_temp_omp, Hreal_array, L_SUB_CODEWORD, Nh, Nv, Nr, L_WHOLE_CODEWORD, L_MIMO_BLK, ModType);
//						int i_code = 0;
//						for (i_code = 0; i_code < L_FINAL_SEARCH; i_code++) // prev: i < L_WHOLE_CODEWORD
//						{
//							for (int i = 0; i < Kv; i++)
//							{
//								for (int j = 0; j < Nh; j++)
//									calcsumh[j] = whole_codeword_list_v[i_code][j * Nv + Av[i]];
//								PolarEncode_xor(calcu1h, calcsumh, Nh);
//								for (int j = 0; j < Kh; j++) { uout[Av[i]][Ah[j]] = calcu1h[Ah[j]]; }
//							}
//							// temporal encoding to check whether the codewords is valid 
//							for (int i = 0; i < Nv; i++)
//							{
//								flag_pass_parity_check_in_crc = true;
//								int* whole_cwd_tmp = new int[Nh];
//								int* whole_cwd_encoded_tmp = new int[Nh];
//								for (int k = 0; k < Nh; k++) { whole_cwd_tmp[k] = 0; whole_cwd_encoded_tmp[k] = 0; }
//								for (int j = 0; j < Nh; j++) { whole_cwd_tmp[j] = whole_codeword_list_v[i_code][j * Nv + i]; }
//								PolarEncode_xor(whole_cwd_encoded_tmp, whole_cwd_tmp, Nh);
//								//cout << endl;
//								for (int j = 0; j < (Nh - Kh); j++) { /*cout << whole_cwd_encoded_tmp[Ach[j]]<<"  ";*/ if (whole_cwd_encoded_tmp[Ach[j]] != 0) flag_pass_parity_check_in_crc = false; }
//								delete[] whole_cwd_tmp;
//								delete[] whole_cwd_encoded_tmp;
//							}
//							//for (int i = 0; i < Kh; i++)
//							//{
//							//	for (int j = 0; j < Nv; j++)
//							//		calcsum[j] = whole_codeword_list[i_code][j * Nh + Ah[i]]; // Nv*Ah[i] + j
//							//	PolarEncode_xor(calcu1, calcsum, Nv);
//							//	for (int j = 0; j < Kv; j++)
//							//		uout[Av[j]][Ah[i]] = calcu1[Av[j]];
//							//}
//							// copy from uout to u_crc
//							for (int i = 0; i < Kv; i++)
//							{
//								for (int j = 0; j < Kh; j++)
//									u_crc[i * Kh + j] = uout[Av[i]][Ah[j]];
//							}
//							if (rx_check_crc(u_crc, Kh * Kv, CRCL)) {
//								// if (i_code == 0) { flag_pass_crc = true; }
//								cout << "FLAG_PAR" << flag_pass_parity_check_in_crc << endl;
//								if (flag_pass_parity_check_in_crc) { // if (true)
//									cout << endl << "Found in " << i_code << endl; flag_found_replacement_codeword = true; idx_replacement_codeword = i_code;
//								}
//								break;
//							}
//							//if (i_code == 0) { break; }
//						}
//						/*cout << endl << i_code << endl;
//						if (i_code == L_WHOLE_CODEWORD) {
//							cout << "NOT FOUND" << endl;
//						}*/
//					}
//				}
//			}
//			// if output is the result of column decoding (meaning that there is another row decoding to be done)
//			else
//			{
//
//				// prev
//				/*for (int i = 0; i < Kh; i++)
//				{
//					for (int j = 0; j < Nv; j++)
//						calcsum[j] = umid[j][Ah[i]];
//					PolarEncode_xor(calcu1, calcsum, Nv);
//					for (int j = 0; j < Kv; j++)
//						uout[Av[j]][Ah[i]] = calcu1[Av[j]];
//				}*/
//
//				// the whole codeword
//				for (int i = 0; i < Nh; i++)
//				{
//					//cout << endl;
//					for (int j = 0; j < Nv; j++)
//					{
//						calcsum[j] = umid[j][i];
//						//cout << umid[j][i] <<"  ";
//					}
//					for (int i_calcu1 = 0; i_calcu1 < Nv; i_calcu1++) { calcu1[i_calcu1] = 0; }
//					PolarEncode_xor(calcu1, calcsum, Nv);
//					for (int j = 0; j < Nv; j++)
//						uout[j][i] = calcu1[j];
//				}
//				//cout << endl << "END TEST UMID" << endl;
//				// output uout
//				/*for (int i = 0; i < Nv; i++)
//				{
//					cout << endl;
//					for (int j = 0; j < Nh; j++) { cout << uout[i][j] << "  "; }
//				}*/
//				// cout << endl;
//				if (CRCL > 0)
//				{
//					// conduct crc check on info bits of uout
//
//
//					flag_pass_crc = false;
//					flag_found_replacement_codeword = false;
//
//
//
//					for (int i = 0; i < Kv; i++)
//					{
//						for (int j = 0; j < Kh; j++)
//							u_crc[i * Kh + j] = uout[Av[i]][Ah[j]];
//					}
//					if (rx_check_crc(u_crc, Kh * Kv, CRCL)) {
//						flag_pass_crc = true;
//					}
//
//					if (!flag_pass_crc)
//					{
//						// decide whether the row codewords can be possibly combined to form a correct codeword
//						bool flag_found_in_all_subcodeword = true;
//						bool flag_found_in_all_subcodeword_v = true;
//						for (int i_row = 0; i_row < Nv; i_row++)
//						{
//							int i_sub = 0;
//							for (i_sub = 0; i_sub < L; i_sub++)
//							{
//								bool is_correct = true;
//								for (int i = 0; i < Nh; i++) { if (sub_encoded_codeword_list[i_row][i_sub * Nh + i] != x_encoded_copy[i_row][i]) is_correct = false; }
//								if (is_correct) break;
//							}
//							if (i_sub == L) { flag_found_in_all_subcodeword = false; }
//							if (!flag_found_in_all_subcodeword) break;
//						}
//
//						for (int i_col = 0; i_col < Nh; i_col++)
//						{
//							int i_sub = 0;
//							for (i_sub = 0; i_sub < L; i_sub++)
//							{
//								bool is_correct = true;
//								for (int i = 0; i < Nv; i++) { if (sub_encoded_codeword_list[i_col][i_sub * Nv + i] != x_encoded_copy[i][i_col]) is_correct = false; }
//								if (is_correct) break;
//							}
//							if (i_sub == L) { flag_found_in_all_subcodeword_v = false; }
//							if (!flag_found_in_all_subcodeword_v) break;
//						}
//
//						if ((!flag_found_in_all_subcodeword) && (!flag_found_in_all_subcodeword_v))
//						{
//							n_err_frame_withno_replacement_codeword++;
//						}
//
//						crc_list_extend(sub_codeword_list, whole_codeword_list, edistance_list, L, Nv, Nh, L_WHOLE_CODEWORD);
//						for (int i_code = 0; i_code < L_WHOLE_CODEWORD; i_code++)
//						{
//							for (int i = 0; i < Kh; i++)
//							{
//								for (int j = 0; j < Nv; j++)
//									calcsum[j] = whole_codeword_list[i_code][j * Nh + Ah[i]];
//								PolarEncode_xor(calcu1, calcsum, Nv);
//								for (int j = 0; j < Kv; j++)
//									uout[Av[j]][Ah[i]] = calcu1[Av[j]];
//							}
//							// copy from uout to u_crc
//							for (int i = 0; i < Kv; i++)
//							{
//								for (int j = 0; j < Kh; j++)
//									u_crc[i * Kh + j] = uout[Av[i]][Ah[j]];
//							}
//							if (rx_check_crc(u_crc, Kh * Kv, CRCL)) {
//								if (i_code == 0) { flag_pass_crc = true; }
//								if (i_code > 0) {
//									cout << endl << "Found in " << i_code << endl; flag_found_replacement_codeword = true; idx_replacement_codeword = i_code;
//									// row-wise parity check to column codewords
//
//								}
//
//								break;
//							}
//							//if (i_code == 0) { break; }
//						}
//					}
//				}
//			}
//			// test: final output
//			/*cout << "final output" << endl;
//			for (int k = 0; k < Nh; k++) { cout << uout[Av[0]][k] << " "; }
//			cout << endl;*/
//			// test: final output
//			//now without CRC
//			int Temp_Error = Num_Error;
//			for (int i = 0; i < Kv; i++)
//			{
//				for (int j = 0; j < Kh; j++)
//				{
//					//printf("%d ", (int)uout[Av[i]][Ah[j]]);
//					Num_Error += abs(u[Av[i]][Ah[j]] - uout[Av[i]][Ah[j]]);
//					/*if (!flag_pass_crc)
//					{
//						cout << abs(u[Av[i]][Ah[j]] - uout[Av[i]][Ah[j]]) << endl;
//					}*/
//				}
//				//cout << endl;
//			}
//
//
//			// whether pass parity check for both row and column
//			bool flag_pass_parity_check = true;
//
//			//  check with x_parity_check
//			for (int i = 0; i < Nh; i++)
//			{
//				//cout << endl;
//				for (int j = 0; j < Nv; j++)
//				{
//					calcsum[j] = x_parity_check[j][i];
//					//cout << umid[j][i] <<"  ";
//				}
//				for (int i_calcu1 = 0; i_calcu1 < Nv; i_calcu1++) { calcu1[i_calcu1] = 0; }
//				PolarEncode_xor(calcu1, calcsum, Nv);
//				// cout << endl;
//				// for (int j = 0; j < Nv; j++) { cout << calcu1[j] << "  "; }
//				for (int k = 0; k < Nv - Kv; k++) { if (calcu1[Acv[k]] != 0) flag_pass_parity_check = false; }
//			}
//
//
//
//			/*for (int i = 0; i < Nv - Kv; i++)
//			{
//				for (int j = 0; j < Nh; j++) { if (uout[Acv[i]][j] != 0) { flag_pass_parity_check = false; break; } }
//				if (!flag_pass_parity_check) break;
//			}
//			for (int i = 0; i<Nh - Kh; i++)
//			{
//				for (int j = 0; j < Nv; j++) { if (uout[j][Ach[i]] != 0) { flag_pass_parity_check = false; break; } }
//				if (!flag_pass_parity_check) break;
//			}*/
//
//			// cout << "PARITY" << flag_pass_parity_check << endl;
//
//			/*if (!flag_pass_crc) {
//				cout << endl;
//			}*/
//			//getchar();
//			if (flag_found_replacement_codeword)
//			{
//				file_crc_check_replacement_codeword << endl;
//				if (Temp_Error < Num_Error) { file_crc_check_replacement_codeword << "Index of Replacement Codeword: " << idx_replacement_codeword << "IsCorrect: " << 0 << endl; }
//				if (Temp_Error == Num_Error) { file_crc_check_replacement_codeword << "Index of Replacement Codeword: " << idx_replacement_codeword << "IsCorrect: " << 1 << endl; }
//			}
//			if (Temp_Error != Num_Error) // prev: <
//			{
//				Num_Frame_Error++;
//				if (flag_pass_crc) { n_err_frame_pass_crc++; }
//				else { n_err_frame_uncorrected_by_crc++; }
//				if (flag_found_replacement_codeword) {
//					n_err_replacement_codeword_crc++;
//				}
//				if (flag_pass_parity_check) { n_err_pass_parity_check++; }
//			}
//			if (Temp_Error == Num_Error)
//			{
//				/*if (!flag_pass_parity_check) {
//					cout << endl << "FFFFFFFUCKKKKKK" << endl;
//					for (int i = 0; i < Nv; i++)
//					{
//						cout << endl;
//						for (int j = 0; j < Nh; j++) { cout << uout[i][j] << "  "; }
//					}
//					cout << endl;
//				}*/
//				if (!flag_pass_crc) { n_err_frame_corrected_by_crc++; }
//			}
//			if (frame % 10 == 0)
//			{
//				/*printf("\r<%d> SNR= %.2f FER= %d / %d = %f BER= %d / %d = %lf E-Det = %f, E-MRB = %f"  ,
//					frame, nSNR, Num_Frame_Error, frame, (double)Num_Frame_Error / frame, Num_Error, (N_Info * frame), (double)Num_Error / N_Info / frame,
//					(double)n_det_err/total_bit_cnt, (double)n_det_err_mrb/(total_bit_cnt/(ModType/2)));*/
//				int L_sdis_sum = L_LSD;
//				sdis_sum = 0;
//				for (int i = 0; i < Nv; i++) sdis_sum += sum_sdis_array[i];
//				if (polarde_array[0] == "SCL") L_sdis_sum = L;
//				// cout << "sdis_sum_avg" << sdis_sum / frame / L_sdis_sum << endl;
//				printf("\r<%d> SNR= %.2f FER= %d / %d = %f BER= %d / %d = %lf ,sdis_avg = %f, E-MRB = %f",
//					frame, nSNR, Num_Frame_Error, frame, (double)Num_Frame_Error / frame, Num_Error, (N_Info * frame), (double)Num_Error / N_Info / frame, (double)sdis_sum / frame / L_sdis_sum, (double)n_det_err / total_bit_cnt);
//				fflush(stdout);
//			}
//
//
//			// output interleaving pattern
//			//if (Num_Frame_Error % record_frames == 0 && Num_Frame_Error > 0 && Num_Frame_Error!=n_frame_error_pre)
//			//{
//			//	n_frame_error_pre = Num_Frame_Error;
//			//	cout << endl;
//			//	int n_file = Num_Frame_Error / record_frames;
//			//	frame_delta = frame - frame_tmp;
//			//	cout << "Delta Frames: " << frame_delta << endl;
//			//	frame_tmp = frame;
//			//	string ofile_intrlv_name = string("intrlv_pattern_") + to_string(n_file) + string("_frame_") + to_string(frame_delta) + string(".txt");
//			//	string ofile_spatial_intrlv_name = string("intrlv_pattern_spatial") + to_string(n_file) + string(".txt");
//			//	ofstream ofile_intrlv(ofile_intrlv_name);
//			//	ofstream ofile_spatial_intrlv(ofile_spatial_intrlv_name);
//			//	for (int i = 0; i < Nv; i++)
//			//	{
//			//		for (int j = 0; j < intrlv_pattern_num; j++)
//			//		{
//
//			//			for (int k = 0; k < ModType; k++)
//			//				ofile_intrlv << intrlv_pattern_2D_time_partial[i][j][k] << " ";
//			//			//cout << endl;
//			//		}
//			//	}
//			//	ofile_intrlv.close();
//
//			//	for (int i = 0; i < Nh; i++)
//			//	{
//			//		for (int j = 0; j < Nv; j++)
//			//			ofile_spatial_intrlv << intrlv_pattern_2D[i][j] << " ";
//			//	}
//			//	ofile_spatial_intrlv.close();
//
//			//	// re-generate interleaving pattern
//			//	generateIntrlvPattern(intrlv_pattern, Nv);
//			//	for (int i = 0; i < Nh; i++) { generateIntrlvPattern(intrlv_pattern_2D[i], Nv); }
//			//	for (int i = 0; i < Nv; i++) { generateIntrlvPattern(intrlv_pattern_2D_time[i], Nh); }
//			//	if (flag_interleave_partial)
//			//	{
//			//		for (int i = 0; i < Nv; i++)
//			//		{
//			//			for (int j = 0; j < intrlv_pattern_num; j++)
//			//			{
//			//				generateIntrlvPattern(intrlv_pattern_2D_time_partial[i][j], ModType);
//			//				/*for (int k = 0; k < ModType; k++)
//			//					cout << intrlv_pattern_2D_time_partial[i][j][k];
//			//				cout << endl;*/
//			//			}
//			//		}
//			//	}
//			//	//Num_Frame_Error = 0;
//			//}
//
//
//			//printf("%d %d", Num_Frame_Error, Num_Error);
//			//getchar();
//			//Decoding_Duration += (T_finish - T_start);
//			/*auto dec_end = std::chrono::system_clock::now();
//			auto duration_decoding = std::chrono::duration_cast<std::chrono::microseconds>(dec_end - det_end);
//			auto duration_total = std::chrono::duration_cast<std::chrono::microseconds>(dec_end - start);
//			file_dec << duration_decoding.count() << endl;
//			file_total << duration_total.count() << endl;*/
//			if (frame % 16 == 0)
//			{
//				file_err_num << Num_Error << endl;
//			}
//
//
//			// logging crc output
//			if ((frame % 1000) == 0)
//			{
//				file_crc_check_decoding << endl;
//				file_crc_check_decoding << "SNR: " << nSNR << endl;
//				file_crc_check_decoding << "FER: " << (double)Num_Frame_Error / frame << endl;
//				file_crc_check_decoding << "Total Frames: " << frame << endl;
//				file_crc_check_decoding << "Error frames undetected by crc: " << n_err_frame_pass_crc << endl;
//				file_crc_check_decoding << "Error frames corrected by crc: " << n_err_frame_corrected_by_crc << endl;
//				file_crc_check_decoding << "Error frames failed to be corrected by crc: " << n_err_frame_uncorrected_by_crc << endl;
//				file_crc_check_decoding << "Error of replacement codewords: " << n_err_replacement_codeword_crc << endl;
//				file_crc_check_decoding << "Error frames that pass parity check: " << n_err_pass_parity_check << endl;
//				file_crc_check_decoding << "Error frames with no possible replacement codewords: " << n_err_frame_withno_replacement_codeword << endl;
//				file_crc_check_decoding << "Position of correct codeword within each list of vertical codewords" << endl;
//				for (int i = 0; i < Nh; i++)
//				{
//					file_crc_check_decoding << n_err_codeword_v[i] << " ";
//				}
//				file_crc_check_decoding << endl << "Num of correct codewords in each position" << endl;
//				for (int i = 0; i < Nh; i++)
//				{
//					file_crc_check_decoding << i + 1 << "th codeword: ";
//					for (int j = 0; j < L + 1; j++)
//						file_crc_check_decoding << n_err_codeword_v_diffpos[i][j] << "  ";
//					file_crc_check_decoding << endl;
//				}
//				file_crc_check_decoding << endl << "Max correct codeword idx when the correct codeword is found" << endl;
//				for (int i = 0; i < L + 1; i++)
//					file_crc_check_decoding << max_L_cnt[i] << "  ";
//				file_crc_check_decoding << endl;
//			}
//
//		}
//		//Decoding_Duration = Decoding_Duration * 100 / CLOCKS_PER_SEC / frame;
//		//Decoding_Duration = Decoding_Duration / CLOCKS_PER_SEC / Max_frame;
//
//		//printf("%f  %f  %f  %f\n", nSNR, Decoding_Duration, (double)Num_Error / N_Info / Max_frame, (double)Num_Frame_Error / Max_frame);
//		cout << endl;
//		printf("%f %f %f\n", nSNR, (double)Num_Frame_Error / frame, (double)Num_Error / N_Info / frame);
//		cout << endl;
//		if (bit_flip_pyndiah)
//		{
//			for (int i = 0; i < softiter; i++)
//			{
//				cout << "# Tried Patterns (Horizontal):   " << setw(10) << (float)n_tried_patterns_h_all[i] / (float)(Nv * frame);
//				cout << "# Tried Patterns (Vertical):   " << setw(10) << (float)n_tried_patterns_v_all[i] / (float)(Nh * frame);
//				cout << endl;
//			}
//		}
//		BER_array[SNR_count] = (double)Num_Error / N_Info / frame;
//		FER_array[SNR_count] = (double)Num_Frame_Error / frame;
//		// save to .txt file
//		string polarde_serial = "";
//		for (int i = 0; i < softiter; i++)
//			polarde_serial = polarde_serial + polarde_array[i] + string("_");
//		string file_name = string("ResFinal\\2D_LSD_D_") + polarde_serial + MIMODET + string("T") + to_string(Nt) + string("R") + to_string(Nr)
//			+ string("Q") + to_string(int(pow(2, ModType))) + string("L") + to_string(Nh) + string("K") + to_string(Kh)
//			+ string("T") + to_string(softiter) + string("H") + to_string(half_iter) + string("UB") + to_string(usebeta) + string("I") + to_string(flag_interleave)
//			+ string(".txt");
//		save_array_txt(nSNR, BER_array[SNR_count], FER_array[SNR_count], file_name);
//		int L_sdis_sum = L_LSD;
//		sdis_sum = 0;
//		for (int i = 0; i < Nv; i++) sdis_sum += sum_sdis_array[i];
//		if (polarde_array[0] == "SCL") L_sdis_sum = L;
//		cout << "sdis_sum_avg" << sdis_sum / frame / L_sdis_sum << endl;
//		SNR_count++;
//		delete[] n_tried_patterns_h;
//		delete[] n_tried_patterns_v;
//
//
//		// logging crc details
//		if (true)
//		{
//			file_crc_check_decoding << endl;
//			file_crc_check_decoding << "FER: " << (double)Num_Frame_Error / frame << endl;
//			file_crc_check_decoding << "Total Frames: " << frame << endl;
//			file_crc_check_decoding << "Error frames undetected by crc: " << n_err_frame_pass_crc << endl;
//			file_crc_check_decoding << "Error frames corrected by crc: " << n_err_frame_corrected_by_crc << endl;
//			file_crc_check_decoding << "Error frames failed to be detected by crc: " << n_err_frame_uncorrected_by_crc << endl;
//			file_crc_check_decoding << "Error of replacement codewords: " << n_err_replacement_codeword_crc << endl;
//		}
//
//	}
//	//fclose(stdout);
//	delete[] a_data;
//	for (int i = 0; i < Nv; i++)
//		delete[] u[i];
//	delete[] u;
//	for (int i = 0; i < Nv; i++)
//		delete[] umid[i];
//	delete[] umid;
//	for (int i = 0; i < Nv; i++)
//		delete[] uout[i];
//	delete[] uout;
//	for (int i = 0; i < Nv; i++)
//		delete[] x[i];
//	delete[] x;
//	for (int i = 0; i < Nv; i++)
//		delete[] x_mmse_out[i];
//	delete[] x_mmse_out;
//	for (int i = 0; i < Nv; i++)
//		delete[] midx[i];
//	delete[] midx;
//	for (int i = 0; i < Nv; i++)
//		delete[] y[i];	for (int i = 0; i < Nh; i++)
//		delete[] oriLLR_kbest[i];
//	delete[] oriLLR_kbest;
//	for (int i = 0; i < Nv; i++)
//		delete sub_codeword_list[i];
//	delete[] sub_codeword_list;
//	for (int i = 0; i < L_WHOLE_CODEWORD; i++)
//		delete whole_codeword_list[i];
//	delete[] whole_codeword_list;
//	for (int i = 0; i < Nv; i++)
//		delete[]edistance_list[i];
//	delete[] edistance_list;
//	delete[] y;
//	delete[] u_crc;
//	delete[] Ah;
//	delete[] Ach;
//	delete[] A_Ach;
//	delete[] Av;
//	delete[] Acv;
//	delete[] A_Acv;
//	for (int i = 0; i < Nv; i++)
//		delete[] oriLLRh[i];
//	delete[] oriLLRh;
//	for (int i = 0; i < Nv; i++)
//		delete[] upLLR[i];
//	delete[] upLLR;
//	for (int i = 0; i < Nh; i++)
//		delete[] GGh[i];
//	delete[] GGh;
//	for (int i = 0; i < Nv; i++)
//		delete[] GGv[i];
//	delete[] GGv;
//	for (int i = 0; i < mimo_blk_height; i++) {
//		delete[] x_mod[i];
//	}
//	delete[] x_mod;
//	for (int i = 0; i < Nr; i++) {
//		delete[] y_mod[i];
//	}
//	delete[] y_mod;
//	for (int i = 0; i < 2 * mimo_blk_height; i++) {
//		delete[] y_mmse[i];
//	}
//	delete[] y_mmse;
//	for (int i = 0; i < 2 * Nr; i++) {
//		delete[] y_mod_real[i];
//	}
//	for (int i = 0; i < mimo_blk_length; i++)
//	{
//		delete[] y_mod_real_temp_omp[i];
//	}
//	delete[] y_mod_real_temp_omp;
//	for (int i = 0; i < mimo_blk_length; i++)
//	{
//		delete[] y_mmse_temp_omp[i];
//	}
//	delete[] y_mmse_temp_omp;
//	for (int i = 0; i < mimo_blk_length; i++)
//	{
//		delete[] ep_extr_mean_omp[i];
//	}
//	delete[] ep_extr_mean_omp;
//	for (int i = 0; i < mimo_blk_length; i++)
//	{
//		delete[] ep_extr_var_omp[i];
//	}
//	for (int i = 0; i < Nv; ++i) {
//		for (int j = 0; j < n_samples; ++j) {
//			delete[] tried_codewords_h[i][j];  // �ͷŵ���ά����������
//		}
//		delete[] tried_codewords_h[i];  // �ͷŵڶ�ά��ָ������
//	}
//	delete[] tried_codewords_h;  // �ͷŵ�һά��ָ������
//	for (int i = 0; i < Nh; ++i) {
//		for (int j = 0; j < n_samples; ++j) {
//			delete[] tried_codewords_v[i][j];  // �ͷŵ���ά����������
//		}
//		delete[] tried_codewords_v[i];  // �ͷŵڶ�ά��ָ������
//	}
//	delete[] tried_codewords_v;  // �ͷŵ�һά��ָ������
//	delete[] ep_extr_var_omp;
//	delete[] y_mod_real;
//	delete[] x_mod_tmp;
//	delete[] y_mod_tmp;
//	delete[] y_mod_real_tmp;
//	delete[] x_tmp;
//	delete[] oriLLR_tmp;
//	delete[] sum_sdis_array;
//	delete[] Hreal_array;
//	delete[] HTH_array;
//	delete[] HTY_array;
//	delete[] intrlv_pattern_2D_inner;
//}


// system function for GenAlg of Irregular Product Coding
//void system2D_LSD_epdet_irregular_SVD_precode_GenAlg(double* BER_array, double* FER_array, int* channel_index_predefined, int* low_rate_index)
//{
//	std::ofstream file_total("total_timing.txt", std::ios::app);
//	std::ofstream file_det("det_timing.txt", std::ios::app);
//	std::ofstream file_dec("dec_timing.txt", std::ios::app);
//	std::ofstream file_err_num("total_frames.txt", std::ios::app);
//	auto det_end = std::chrono::system_clock::now();
//	int Nt2 = 2 * Nt;
//	int Nr2 = 2 * Nr;
//	int cnt_short_codewords = 0;
//	for (int i = 0; i < Nv; i++) if (channel_index_predefined[i] == 1) { cnt_short_codewords++; }
//	int* bad_channel_index = new int[cnt_short_codewords];
//	int* good_channel_index = new int[Nv - cnt_short_codewords];
//	int cnt_good = 0;
//	int cnt_bad = 0;
//	for (int i = 0; i < Nv; i++) { if (channel_index_predefined[i] == 0) { good_channel_index[cnt_good] = i; cnt_good++; } else { bad_channel_index[cnt_bad] = i; cnt_bad++; } }
//	/*if (!file.is_open()) {
//		std::cerr << "Error opening file: " << file_name << std::endl;
//		return;
//	}*/
//
//	//file << setw(20) << SNR << setw(20) << BER << setw(20) << FER << "\n";
//
//	//freopen("result.txt", "w", stdout);
//	printsetting();
//	cout << "\nSNR	 " << "FER      BER" << endl;
//	int N_Info = Kv * Kh - CRCL, Nmax;
//	if (irregular_coding) { N_Info = border * K_time_array[0] + (Kv - border) * K_time_array[1] - CRCL; }
//	//Interleaving pattern
//	/*std::vector<int> intrlv_pattern(Nv);
//	std::vector<int>* intrlv_pattern_2D = new vector<int>[Nh];
//	for (int i = 0; i < Nh; i++) intrlv_pattern_2D[i].resize(Nv);*/
//	int intrlv_pattern_num = Nh / ModType;
//	std::vector<int> intrlv_pattern(Nv);
//	std::vector<int> intrlv_pattern_info_mapping(Nv);
//	std::vector<int>* intrlv_pattern_2D = new vector<int>[Nh];
//	std::vector<int>* intrlv_pattern_2D_time = new vector<int>[Nv];
//	std::vector<int>* intrlv_pattern_2D_inner = new vector<int>[Nh];
//	std::vector<int>** intrlv_pattern_2D_time_partial = new vector<int> *[Nv];
//	for (int i = 0; i < Nv; i++) intrlv_pattern_2D_time_partial[i] = new vector<int>[intrlv_pattern_num];
//	for (int i = 0; i < Nv; i++)
//	{
//		for (int j = 0; j < intrlv_pattern_num; j++)
//		{
//			intrlv_pattern_2D_time_partial[i][j].resize(ModType);
//		}
//	}
//	for (int i = 0; i < Nh; i++) intrlv_pattern_2D[i].resize(Nv);
//	for (int i = 0; i < Nh; i++) intrlv_pattern_2D_inner[i].resize(Kv);
//	for (int i = 0; i < Nv; i++) intrlv_pattern_2D_time[i].resize(Nh);
//
//	float CodeRate = (float)(N_Info) / (Nh * Nv);
//	Nmax = max(Nh, Nv);
//	// num of symbols per vertical codeword
//	int Nsym = 0;
//	if (MOD_DIM == "Spatial") { Nsym = Nv / ModType; }
//	else if (MOD_DIM == "Temporal") { Nsym = Nh / ModType; }
//	int symbol_length = ModType;
//
//	// Size of a mimo block is given by height(spatial) * length(temporal)
//	int mimo_blk_height = 0;
//	int mimo_blk_length = 0;
//	if (MOD_DIM == "Spatial")
//	{
//		mimo_blk_height = Nsym;
//		mimo_blk_length = Nh;
//	}
//	else if (MOD_DIM == "Temporal")
//	{
//		mimo_blk_height = Nv;
//		mimo_blk_length = Nsym;
//	}
//	// Initialization
//	unsigned char* a_data = new unsigned char[N_Info];// Information CodeWord
//	int** u = new int* [Nv];// Original CodeWord u
//	for (int i = 0; i < Nv; i++) { u[i] = new int[Nh]; }
//	int** umid = new int* [Nv];// output CodeWord u
//	for (int i = 0; i < Nv; i++) { umid[i] = new int[Nh]; }
//	int** uout = new int* [Nv];// output CodeWord u
//	for (int i = 0; i < Nv; i++) { uout[i] = new int[Nh]; }
//	int** midx = new int* [Nv];//  Polar Encoded CodeWord x
//	for (int i = 0; i < Nv; i++) { midx[i] = new int[Nh]; }
//	int** x = new int* [Nv];//  Polar Encoded CodeWord x
//	for (int i = 0; i < Nv; i++) { x[i] = new int[Nh]; }
//	int** x_parity_check = new int* [Nv];
//	for (int i = 0; i < Nh; i++) x_parity_check[i] = new int[Nh];
//
//	// Arrays for MIMO
//	int** x_mmse_out = new int* [Nv]; // delete!
//	for (int i = 0; i < Nv; i++) { x_mmse_out[i] = new int[Nh]; }
//	std::complex<float>** x_mod = new std::complex<float>*[mimo_blk_height]; // modulated symbols (spatial modulation) DELETE!
//	for (int i = 0; i < mimo_blk_height; i++) { x_mod[i] = new std::complex<float>[mimo_blk_length]; }
//	std::complex<float>** y_mod = new std::complex<float>*[Nr]; // Output symbols after AWGN channels // delete!
//	for (int i = 0; i < Nr; i++) { y_mod[i] = new std::complex<float>[mimo_blk_length]; }
//	float** y_mod_real = new float* [2 * Nr];
//	for (int i = 0; i < 2 * Nr; i++) y_mod_real[i] = new float[mimo_blk_length];
//	float** y = new float* [Nv];// Output Signal After AWGN Channel
//	for (int i = 0; i < Nv; i++) { y[i] = new float[Nh]; }
//	float** y_mmse = new float* [2 * mimo_blk_height];
//	for (int i = 0; i < 2 * mimo_blk_height; i++) { y_mmse[i] = new float[mimo_blk_length]; }  // MMSE detector output // delete!
//	float* y_mmse_tmp = new float[2 * mimo_blk_height];
//	std::complex<float>* x_mod_tmp = new std::complex<float>[mimo_blk_height]; // modulated symbols of a single column // delete!
//	std::complex<float>* y_mod_tmp = new std::complex<float>[Nr]; // received symbols of a single column // delete!
//	std::complex<float>* x_precoded_tmp = new std::complex<float>[mimo_blk_height]; // precoded symbols of a single column // delete!
//	std::complex<float>* y_precoded_tmp = new std::complex<float>[Nr];
//	float* y_mod_real_tmp = new float[2 * Nr];
//	// float* y_mod_tmp = new float[Nr];  // delete!
//	float** y_mod_real_temp_omp = new float* [mimo_blk_length]; // new delete!
//	for (int i = 0; i < mimo_blk_length; i++) { y_mod_real_temp_omp[i] = new float[2 * Nr]; }
//	float** y_mmse_temp_omp = new float* [mimo_blk_length]; // new delete!
//	for (int i = 0; i < mimo_blk_length; i++) { y_mmse_temp_omp[i] = new float[2 * mimo_blk_height]; }
//	// ep
//	float** ep_extr_mean_omp = new float* [mimo_blk_length];
//	for (int i = 0; i < mimo_blk_length; i++) { ep_extr_mean_omp[i] = new float[2 * mimo_blk_height]; }
//	float** ep_extr_var_omp = new float* [mimo_blk_length];
//	for (int i = 0; i < mimo_blk_length; i++) { ep_extr_var_omp[i] = new float[2 * mimo_blk_height]; }
//	int* x_tmp = new int[Nv]; // delete!
//
//	// paths selected by k-best detection
//
//
//	// CRC bits
//	unsigned char* u_crc = new unsigned char[Kv * Kh];
//
//	float** upLLR = new float* [Nv];
//	for (int i = 0; i < Nv; i++) { upLLR[i] = new _MM_ALIGN16 float[Nh]; }
//	//float** oriLLRv = new float* [Nh];
//	//for (int i = 0; i < Nh; i++) { oriLLRv[i] = new _MM_ALIGN16 float[2 * Nv + 2]; }
//	float** oriLLRh = new float* [Nv];
//	for (int i = 0; i < Nv; i++) { oriLLRh[i] = new _MM_ALIGN16 float[Nh]; }
//	int len_ori_llr_tmp = 0;
//	if (MOD_DIM == "Spatial") { len_ori_llr_tmp = Nv; }
//	else if (MOD_DIM == "Temporal") { len_ori_llr_tmp = ModType * Nt; }
//	float* oriLLR_tmp = new _MM_ALIGN16 float[len_ori_llr_tmp];  // delete!
//	//float** oriLLR_kbest = new float* [Nh];
//	float** oriLLR_kbest = new float* [mimo_blk_length];
//	for (int i = 0; i < Nh; i++) { oriLLR_kbest[i] = new float[len_ori_llr_tmp]; }
//	double* sum_sdis_array = new double[Nv];
//	for (int i = 0; i < Nv; i++) sum_sdis_array[i] = 0;
//	int* n_tried_patterns_h = new int[softiter];
//	for (int i = 0; i < softiter; i++) n_tried_patterns_h[i] = 0;
//	int* n_tried_patterns_v = new int[softiter];
//	for (int i = 0; i < softiter; i++) n_tried_patterns_v[i] = 0;
//
//	//*****************************************************************
//	// EP detector
//	float* Hreal_array = new float[(2 * Nr) * (2 * Nt)]; // delete
//	float* HTH_array = new float[(2 * Nt) * (2 * Nt)]; // delete
//	float* HTY_array = new float[2 * Nt]; // delete 
//
//
//
//
//
//	int** GGh = new int* [Nh];
//	int** GGv = new int* [Nv];
//
//
//	MatrixXf Gh = Fn(Nh);  // Generate Matrix of horizontal coding
//	MatrixXf Gv = Fn(Nv);  // Generate Matrix of vertical coding
//	MatrixXcf H(Nr, Nt); // Complex channel Matrix H: Nr x Nt
//	//Matrix<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic> H = MatrixXcf::Random(Nr, Nt); // Real domain channel matrix
//	//// //MatrixX Definitions
//	MatrixXf H_real(2 * Nr, 2 * Nt); // Real domain channel matrix
//	MatrixXf received_signal(2 * Nr, 1);
//	MatrixXf HTH(2 * Nt, 2 * Nt);
//	MatrixXf Identity = MatrixXf::Identity(2 * Nt, 2 * Nt);
//	MatrixXf output_signal(2 * Nt, 1);
//	MatrixXf mmse_matrix(2 * Nt, 2 * Nr);
//	MatrixXf HTHIinv(2 * Nt, 2 * Nt);
//
//	// for KBest
//	MatrixXf R(2 * Nt, 2 * Nt);
//	MatrixXf Q(2 * Nr, 2 * Nt);
//	MatrixXf z(2 * Nt, 1);
//
//	for (int i = 0; i < Nh; i++) GGh[i] = new int[Nh];
//	for (int i = 0; i < Nv; i++) GGv[i] = new int[Nv];
//	for (int i = 0; i < Nh; i++)
//		for (int j = 0; j < Nh; j++)
//			GGh[i][j] = Gh(i, j);
//	for (int i = 0; i < Nv; i++)
//		for (int j = 0; j < Nv; j++)
//			GGv[i][j] = Gv(i, j);
//
//	// display GGh
//	std::cout << endl;
//	/*for (int i = 0; i < Nv; i++)
//	{
//		int nnz_cnt = 0;
//		for (int j = 0; j < Nv; j++)
//		{
//			cout << GGv[i][j] << " ";
//			nnz_cnt += GGv[i][j];
//		}
//		cout << nnz_cnt << " " << endl;
//	}*/
//
//	//*****************************************************************
//	vector<int> temph(Nh);
//	vector<int> tempv(Nv);
//	vector<int>* construction_time_all = new vector<int>[Nv];
//	for (int i = 0; i < Nv; i++) { construction_time_all[i].resize(Nh); }
//	int* Ah = new int[Kh]; // Horizontal Info Bits
//	int* Ach = new int[Nh - Kh]; // Horizontal Frozen Bits
//	int* A_Ach = new int[Nh];
//	int* Av = new int[Kv];
//	int* Acv = new int[Nv - Kv];
//	int* A_Acv = new int[Nv];
//	int** A_time_all = new int* [Nv];
//	for (int i = 0; i < Nv; i++) { A_time_all[i] = new int[Nh]; for (int j = 0; j < Nh; j++) A_time_all[i][j] = 0; }
//	int** Ac_time_all = new int* [Nv];
//	for (int i = 0; i < Nv; i++) { Ac_time_all[i] = new int[Nh]; for (int j = 0; j < Nh; j++) Ac_time_all[i][j] = 0; }
//	int** A_Ac_time_all = new int* [Nv];
//	for (int i = 0; i < Nv; i++) { A_Ac_time_all[i] = new int[Nh]; for (int j = 0; j < Nh; j++) A_Ac_time_all[i][j] = 0; }
//	int* K_time_all = new int[Nv];
//	for (int i = 0; i < Nv; i++) K_time_all[i] = 0;
//	int* tempv_ordered = new int[Nv];
//	int** all_error_positions = new int* [Nv];
//	for (int i = 0; i < Nv; i++) {
//		all_error_positions[i] = new int[Nh];
//		for (int j = 0; j < Nh; j++) all_error_positions[i][j] = 0;
//	}
//	int* use_short_codeword = new int[Nv];
//	for (int i = 0; i < Nv; i++) { use_short_codeword[i] = 0; }
//	/*int* use_short_codeword = new int[Nv];
//	for (int i = 0; i < Nv; i++) { use_short_codeword[i] = 0; }*/
//
//	float** channel_scale_factor = new float* [Nr];
//	for (int i = 0; i < Nr; i++) { channel_scale_factor[i] = new float[Nt]; for (int j = 0; j < Nt; j++) { channel_scale_factor[i][j] = 0.0; } }
//	float* singular_values_array = new float[Nt];
//	for (int i = 0; i < Nt; i++) { singular_values_array[i] = 0.0; }
//	float* antenna_power = new float[Nt];
//	for (int i = 0; i < Nt; i++) { antenna_power[i] = 1.0; }
//	//*****************************************************************
//		// Main Function
//	// srand((unsigned int)time(0));
//	std::random_device rd;              // ��ȡ���������
//	std::mt19937 gen(rd());             // ʹ�����ӳ�ʼ��mt19937������
//	std::uniform_int_distribution<> distrib(1, 1000000);
//	int Num_Error, Num_Frame_Error;
//	float Decoding_Duration;
//
//	int Num_Error_new, Num_Frame_Error_new;
//	float Decoding_Duration_new;
//
//	clock_t T_start, T_finish;
//	int SNR_count = 0;
//	int anssc, anssd;
//	for (double nSNR = Min_SNR; nSNR <= Max_SNR; nSNR += SNR_step)
//	{
//		for (int i = 0; i < Nv; i++) sum_sdis_array[i] = 0;
//		double sdis_sum = 0;
//		// TEST: det err
//		int n_det_err = 0; // n error bits at detector output
//		int n_det_err_mrb = 0; // error bits at detector output(MRBs)
//		int total_bit_cnt = 0;
//		// TEST: det err
//		int cntcnt = 0;
//		//Initialization
//		Num_Error = 0;
//		Num_Frame_Error = 0;
//		Decoding_Duration = 0;
//
//		// double tmp = pow(10, nSNR / 10);
//		float sigma_sq = 0;
//		double tmp = 0;
//		if (atten == 1)
//		{
//			tmp = pow(10, nSNR / 10);
//			sigma_sq = (float)1 / (float)tmp / 2 / CodeRate;
//		}
//		if (atten == 2)
//		{
//			tmp = pow(10, nSNR / 10) * symbol_length * CodeRate * Nt;
//			// tmp = pow(10, nSNR / 10) * symbol_length  * Nt;
//			sigma_sq = (Nt * Nr) / tmp;
//		}
//		// cout << sigma_sq << endl;
//		polar_codeconstruction(polarconh, Nh, Kh, GGh, (float)sqrt(sigma_sq), temph);
//		polar_codeconstruction(polarconv, Nv, Kv, GGv, (float)sqrt(sigma_sq), tempv);
//		// for (int i = 0; i < Nv; i++) { tempv[i] = vertical_construction_fixed[i]; }
//		// cout << endl;
//		for (int i = 0; i < Nv; i++) { tempv_ordered[i] = tempv[i]; /*cout << tempv_ordered[i] << ",";*/ }
//
//		// sort temph and tempv in ascending order
//		for (int i = 0; i < Kh - 1; i++)
//			for (int j = i + 1; j < Kh; j++)
//				if (temph[i] > temph[j])
//				{
//					int ttt = temph[i];
//					temph[i] = temph[j];
//					temph[j] = ttt;
//				}
//		for (int i = 0; i < Kv - 1; i++)
//			for (int j = i + 1; j < Kv; j++)
//				if (tempv[i] > tempv[j])
//				{
//					int ttt = tempv[i];
//					tempv[i] = tempv[j];
//					tempv[j] = ttt;
//				}
//
//		// Decide the position of information bits and frozen bits
//		for (int i = 0; i < Kh; i++) // Information Set
//		{
//			Ah[i] = temph[i];
//			A_Ach[Ah[i]] = 1;
//		}
//		//	getchar();
//		for (int i = 0; i < Nh - Kh; i++) // Frozen Bit
//		{
//			Ach[i] = temph[Kh + i];
//			A_Ach[Ach[i]] = 0;
//		}
//		for (int i = 0; i < Kv; i++) // Information Set
//		{
//			Av[i] = tempv[i];
//			A_Acv[Av[i]] = 1;
//		}
//		//	getchar();
//		for (int i = 0; i < Nv - Kv; i++) // Frozen Bit
//		{
//			Acv[i] = tempv[Kv + i];
//			A_Acv[Acv[i]] = 0;
//		}
//
//		/*for (int i = 0; i < Nh; i++)
//		{
//			cout << A_Ach[i] << ",";
//		}
//		cout << endl;*/
//		// for (int i = 0; i < Kv; i++) { cout << Av[i] << endl; }
//		for (int k = 0; k < Nv; k++)
//		{
//			bool bfind = false;
//			int channel_index = 0;
//			for (int i_find = 0; i_find < Kv; i_find++) { if (tempv_ordered[i_find] == k) { bfind = true;  channel_index = i_find; break; } }
//			int k_time = K_time_array[1];
//			int k_thres = Av[border];
//			//cout << k_thres << endl;
//			//cout << channel_index << endl;
//			// if ((k < k_thres) && A_Acv[k]==1) { k_time = K_time_array[0]; }
//			if (bfind && ((Kv - 1 - channel_index) < border)) { k_time = K_time_array[0]; }
//			K_time_all[k] = k_time;
//			polar_codeconstruction(polarconh, Nh, k_time, GGh, (float)sqrt(sigma_sq), construction_time_all[k]);
//			if (k_time == K_time_array[0]) { use_short_codeword[k] = 1; } // use long codeword
//			else if (k_time == K_time_array[1]) { use_short_codeword[k] = 0; } // use short codeword
//			for (int i = 0; i < k_time - 1; i++)
//				for (int j = i + 1; j < k_time; j++)
//					if (construction_time_all[k][i] > construction_time_all[k][j])
//					{
//						int ttt = construction_time_all[k][i];
//						construction_time_all[k][i] = construction_time_all[k][j];
//						construction_time_all[k][j] = ttt;
//					}
//
//			// Decide the position of information bits and frozen bits
//			for (int i = 0; i < k_time; i++) // Information Set
//			{
//				A_time_all[k][i] = construction_time_all[k][i];
//				A_Ac_time_all[k][A_time_all[k][i]] = 1;
//			}
//			//	getchar();
//			for (int i = 0; i < Nh - k_time; i++) // Frozen Bit
//			{
//				Ac_time_all[k][i] = construction_time_all[k][k_time + i];
//				A_Ac_time_all[k][Ac_time_all[k][i]] = 0;
//			}
//		}
//		/*cout << endl;
//		cout << "int channel_index_predefined[Nv] = {";
//		for (int i = 0; i < Nv; i++) { cout << use_short_codeword[i] << ","; }
//		cout << "};" << endl;*/
//		// Non-uniform channel
//		vector<int> non_uniform_intrlv_pattern;
//		vector<int>* non_uniform_intrlv_pattern_array = new vector<int>[Nh];
//		non_uniform_intrlv_pattern.resize(Nv);
//		for (int i = 0; i < Nh; i++) { non_uniform_intrlv_pattern_array[i].resize(Nv); }
//		if (non_uniform_channel)
//		{
//			getChannelScalingFactor(Nr, Nt, diff_gain, use_short_codeword, channel_scale_factor, predefined_non_uniform);
//			if (predefined_non_uniform && (K_time_array[0] != K_time_array[1])) {
//				generateIntrlvPatternIrregular(non_uniform_intrlv_pattern, Nv, bad_channel_index,
//					good_channel_index, cnt_short_codewords, use_short_codeword);
//				if (flag_separate_interleaving) {
//					for (int i = 0; i < Nh; i++) {
//						generateIntrlvPatternIrregularArray(non_uniform_intrlv_pattern_array[i], Nv, bad_channel_index,
//							good_channel_index, cnt_short_codewords, use_short_codeword);
//					}
//				}
//				//// test
//				/*cout << endl << "bad channel index: ";
//				for (int i = 0; i < cnt_short_codewords; i++) { cout << bad_channel_index[i] << " "; }
//				cout << endl << "good channel index: ";
//				for (int i = 0; i < Nv-cnt_short_codewords; i++) { cout << good_channel_index[i] << " "; }
//				cout << endl << "use short codeword: ";
//				for (int i = 0; i < Nv; i++) { cout << use_short_codeword[i] << " "; }
//				cout << endl << "interleaving pattern: ";*/
//				/*cout << "\n BAD" << endl;
//				for (int i = 0; i < cnt_short_codewords; i++) { cout << bad_channel_index[i] << "(" << non_uniform_intrlv_pattern[bad_channel_index[i]] << ")" << "  "; }
//				cout << "\n GOOD" << endl;
//				for (int i = 0; i < Nv - cnt_short_codewords; i++) { cout << good_channel_index[i] << "(" << non_uniform_intrlv_pattern[good_channel_index[i]] << ")"<<"  "; }
//				cout << endl;*/
//			}
//			/*if (predefined_non_uniform) { getChannelScalingFactor(Nr, Nt, diff_gain, use_short_codeword, channel_scale_factor); }
//			else { getChannelScalingFactor(Nr, Nt, diff_gain, use_short_codeword, channel_scale_factor); }*/
//		}
//		// display A_Ac_time_all row by row
//		//cout << "A_Ac_time_all:" << endl;
//		//for (int i = 0; i < Nv; i++)
//		//{
//		//	for (int j = 0; j < Nh; j++)
//		//	{
//		//		cout << A_Ac_time_all[i][j] << ",";
//		//	}
//		//	cout << endl;
//		//}
//		////// display K_time_all
//		//cout << endl;
//		//for (int i = 0; i < Nv; i++)
//		//{
//		//	cout << K_time_all[i] << endl;
//		//}
//		//test: Horizontal frozen bits
//
//	   /*cout << "Frozen Bits" <<endl;
//	   for (int i = 0; i < Nv - Kv; i++)
//	   {
//		   cout << setw(5) << Acv[i];
//	   }
//	   cout << endl;*/
//	   // test: Horizontal frozen bits
//	   // Frame Loop
//		int passnum;
//		int frame = 0;
//
//		// generate interleaving patterns
//		if (flag_interleave)
//		{
//			/*string ofile_intrlv_name = string("intrlv_pattern_71_frame_49263.txt");
//			string ofile_spatial_intrlv_name = string("intrlv_pattern_spatial71.txt");
//			ifstream ofile_intrlv(ofile_intrlv_name);
//			ifstream ofile_spatial_intrlv(ofile_spatial_intrlv_name);*/
//			int file_in = 1;
//			//for (int i = 0; i < Nv; i++)
//			//{
//			//	for (int j = 0; j < intrlv_pattern_num; j++)
//			//	{
//
//			//		for (int k = 0; k < ModType; k++)
//			//			// ofile_intrlv >> file_in;
//			//			ofile_intrlv >> intrlv_pattern_2D_time_partial[i][j][k];
//			//		//cout << endl;
//			//	}
//			//}
//			//ofile_intrlv.close();
//
//			//for (int i = 0; i < Nh; i++)
//			//{
//			//	for (int j = 0; j < Nv; j++)
//			//		ofile_spatial_intrlv >> intrlv_pattern_2D[i][j];
//			//}
//			/*ofile_spatial_intrlv.close();*/
//			/*for (int i = 0; i < Nh; i++)
//			{
//				for (int j = 0; j < Nv; j++)
//					cout << intrlv_pattern_2D[i][j] << endl;
//			}*/
//			generateIntrlvPattern(intrlv_pattern, Nv);
//			generateIntrlvPatternInfoMapping(intrlv_pattern_info_mapping, Nv, Kv, Av, Acv);
//
//			// display intrlv_pattern_info_mapping
//			// cout << "intrlv_pattern_info_mapping:" << endl;
//			/*for (int i = 0; i < Nv; i++)
//			{
//				cout << i << "  "<< intrlv_pattern_info_mapping[i] << endl;
//			}*/
//			for (int i = 0; i < Nh; i++) { generateIntrlvPattern(intrlv_pattern_2D[i], Nv); }
//			for (int i = 0; i < Nv; i++) { generateIntrlvPattern(intrlv_pattern_2D_time[i], Nh); }
//			if (flag_interleave_partial)
//			{
//				for (int i = 0; i < Nv; i++)
//				{
//					for (int j = 0; j < intrlv_pattern_num; j++)
//					{
//						generateIntrlvPattern(intrlv_pattern_2D_time_partial[i][j], ModType);
//						/*for (int k = 0; k < ModType; k++)
//							cout << intrlv_pattern_2D_time_partial[i][j][k];
//						cout << endl;*/
//					}
//				}
//			}
//			for (int i = 0; i < Nh; i++) { generateIntrlvPattern(intrlv_pattern_2D_inner[i], Kv); }
//		}
//
//		//ofstream ofile_intrlv("intrlv_pattern.txt");
//		//for (int i=0; i<Nv; i++)
//		//{
//		//	for (int j = 0; j < intrlv_pattern_num; j+for+)
//		//	{
//		//		
//		//		for (int k = 0; k < ModType; k++)
//		//			ofile_intrlv << intrlv_pattern_2D_time_partial[i][j][k] << "   ";
//		//		//cout << endl;
//		//	}
//		//}
//		//ofile_intrlv.close();
//		int frame_tmp = 0;
//		int frame_delta = 0;
//		int n_frame_error_pre = 0;
//		while (Num_Frame_Error < errframe)
//		{
//			/*test: processing time*/
//			//auto start = std::chrono::system_clock::now();
//			/*if (frame == 10) {
//				cout << "AAAAAAAAAA" << endl;
//			}*/
//			for (int i = 0; i < Nh; i++) { generateIntrlvPattern(intrlv_pattern_2D[i], Nv); }
//			for (int i = 0; i < Nv; i++) { generateIntrlvPattern(intrlv_pattern_2D_time[i], Nh); }
//			frame++;
//			// generate interleaving pattern
//			/*if (flag_interleave)
//			{
//				generateIntrlvPattern(intrlv_pattern, Nv);
//				for (int i = 0; i < Nh; i++) { generateIntrlvPattern(intrlv_pattern_2D[i], Nv); }
//				for (int i = 0; i < Nv; i++) { generateIntrlvPattern(intrlv_pattern_2D_time[i], Nh); }
//				if (flag_interleave_partial)
//				{
//					for (int i = 0; i < Nv; i++)
//					{
//						for (int j = 0; j < intrlv_pattern_num; j++)
//						{
//							generateIntrlvPattern(intrlv_pattern_2D_time_partial[i][j], ModType);
//						}
//					}
//				}
//			}*/
//
//			//cout << frame << endl;
//			int u_input[16] = { 0,     1     ,0,     1,
//									0,     1,     0,     0,
//									0 ,    1,     0,     0,
//									1 ,    0,     0,     0 };
//			for (int i = 0; i < N_Info; i++)
//			{
//				a_data[i] = (unsigned char)(distrib(gen) % 2);
//				//		a_data[i] = u_input[i];
//			}
//
//			// cout << N_Info << endl;
//			if (CRCL) tx_append_crc(a_data, N_Info, CRCL);
//			for (int i = 0; i < Nv; i++)
//				for (int j = 0; j < Nh; j++)
//					u[i][j] = 0;
//			//omp_set_num_threads(4);
////#pragma omp parallel for
//			int info_count = 0;
//			for (int i = 0; i < Kv; i++)
//			{
//				for (int j = 0; j < K_time_all[Av[i]]; j++)
//				{
//					u[Av[i]][A_time_all[Av[i]][j]] = (int)a_data[info_count];
//					info_count++;
//					//cout << (int)u[Av[i]][Ah[j]] << " \n";
//					//cout << Av[i] << "\n";
//				}
//				//cout << endl; 
//			}
//			// test:original info
//			//cout << "original info u:" << endl;
//			////for (int j = 0; j < Nh; j++) { cout << u[Av[0]][j] << " "; }
//			//for (int i = 0; i < Nv; i++)
//			//{
//			//	cout << endl;
//			//	for (int j = 0; j < Nh; j++)
//			//		cout << u[i][j];
//			//}
//			//cout << endl;
//			// Polar Encoding
////#pragma omp parallel for
//
//			// test: original info
//			/*cout << "original info u:" << endl;
//			display_array(u, Nv, Nh);*/
//			// test: original info
//
//			for (int i = 0; i < Nv; i++)
//				PolarEncode_xor(midx[i], u[i], Nh);
//			// test:horizontal encoding
//			//cout << "u after horizontal encoding:" << endl;
//			////for (int j = 0; j < Nh; j++) { cout << midx[Av[0]][j] << " "; }
//			//for (int i = 0; i < Nv; i++)
//			//{
//			//	cout << endl;
//			//	for (int j = 0; j < Nh; j++)
//			//		cout << midx[i][j];
//			//}
//			//cout << endl;
//			//cout << endl;
//			// test:horizontal encoding
//			// 
//			// column-wise interleaving
//
//			// test: row encoding
//			/*cout << "u after horizontal encoding:" << endl;
//			display_array(midx, Nv, Nh);*/
//			// test: row encoding
//			if (flag_interleave_inner)
//			{
//				int* tmp_to_interleave = new int[Nv];
//				for (int i = 0; i < Nh; i++)
//				{
//					for (int j = 0; j < Nv; j++) { tmp_to_interleave[j] = midx[j][i]; }
//					randIntrlvPartial(intrlv_pattern_2D_inner[i], tmp_to_interleave, Av, Nv, Kv);
//					for (int j = 0; j < Nv; j++) { midx[j][i] = tmp_to_interleave[j]; }
//				}
//				delete[] tmp_to_interleave;
//			}
//
//			// test: innner interleaving
//			/*cout << "u after inner interleaving:" << endl;
//			display_array(midx, Nv, Nh);*/
//			// test: innner interleaving
////#pragma omp parallel for
//			for (int i = 0; i < Nh; i++)
//			{
//				int* tmpxv = new int[Nv]; int* tmpuv = new int[Nv];
//				for (int j = 0; j < Nv; j++)
//					tmpuv[j] = midx[j][i];
//				PolarEncode_xor(tmpxv, tmpuv, Nv);
//				for (int j = 0; j < Nv; j++)
//					x[j][i] = tmpxv[j];
//				delete[] tmpxv; delete[] tmpuv;
//			}
//			/*cout << "u after 1st spatial encoding" << endl;
//			display_array(x, Nv, Nh);*/
//			// systematic spatial encoding
//			// reset frozen bits to 0
//			for (int i = 0; i < Nh; i++)
//			{
//				for (int j = 0; j < Nv - Kv; j++)
//				{
//					x[Acv[j]][i] = 0;
//				}
//			}
//			// another round of spatial encoding
//			for (int i = 0; i < Nh; i++)
//			{
//				int* tmpxv = new int[Nv]; int* tmpuv = new int[Nv];
//				for (int j = 0; j < Nv; j++)
//					tmpuv[j] = x[j][i];
//				PolarEncode_xor(tmpxv, tmpuv, Nv);
//				for (int j = 0; j < Nv; j++)
//					x[j][i] = tmpxv[j];
//				delete[] tmpxv; delete[] tmpuv;
//			}
//			/*cout << "u after 2nd spatial encoding" << endl;
//			display_array(x, Nv, Nh);*/
//
//
//
//			// test: vertical encoding
//			/*cout << "u after vertical encoding:" << endl;
//			display_array(x, Nv, Nh);*/
//			// test: vertical encoding
//			/*for (int i = 0; i < 2 * Nr; i++)
//				for (int j = 0; j < 2 * Nt; j++)
//					H_real(i, j) = 0;*/
//					/*MatrixXf H_real_test(2 * Nr, 2 * Nt);
//					MatrixXf M = H_real_test.transpose();
//					MatrixXf N = H_real_test;
//					HTHIinv = M * N;*/
//					//HTHIinv = M * H_real;
//					// test:vertical encoding
//					//cout << "u after vertical encoding:" << endl;
//					////for (int j = 0; j < Nh; j++) { cout << x[Av[0]][j] << " "; }
//					//for (int i = 0; i < Nv; i++)
//					//{
//					//	cout << endl;
//					//	for (int j = 0; j < Nh; j++)
//					//		cout << x[i][j];
//					//}
//					//cout << endl;
//					//cout << endl;
//					// test:vertical encoding
//			if (atten == 2)
//			{
//				// modulation
//				// #pragma omp parallel for
//				// Interleaving
//				if (flag_interleave)
//				{
//					if (flag_interleave_all)  // interleave temporally
//					{
//						for (int i = 0; i < Nv; i++)
//						{
//							randIntrlv(intrlv_pattern_2D_time[i], x[i], Nh);
//						}
//					}
//					//// test: to be interleaved
//					//cout << "x before interleaving:" << endl;
//					//// print x line by line, add spaces between each 4 terms
//					//for (int i = 0; i < Nv; i++)
//					//{
//					//	for (int j = 0; j < Nh; j++)
//					//	{
//					//		if (!(j % 4)) cout << " ";
//					//		cout << x[i][j] << " ";
//					//	}
//					//	cout << endl;
//					//}
//
//					if (flag_interleave_partial)
//					{
//						for (int i = 0; i < Nv; i++)
//						{
//							for (int j = 0; j < intrlv_pattern_num; j++)
//							{
//								randIntrlv(intrlv_pattern_2D_time_partial[i][j], x[i] + j * ModType, ModType);
//							}
//						}
//					}
//
//					//// test: to be interleaved
//					//cout << "x after interleaving:" << endl;
//					//// print x line by line, add spaces between each 4 terms
//					//for (int i = 0; i < Nv; i++)
//					//{
//					//	for (int j = 0; j < Nh; j++)
//					//	{
//					//		if (!(j % 4)) cout << " ";
//					//		cout << x[i][j] << " ";
//					//	}
//					//	cout << endl;
//					//}
//					if ((!flag_interleave_block) && flag_interleave_spatial)
//					{
//						for (int j = 0; j < Nh; j++)
//						{
//							for (int i = 0; i < Nv; i++) { x_tmp[i] = x[i][j]; }
//							if (K_time_array[0] == K_time_array[1]) {
//								if (!flag_information_mapping)
//								{
//									randIntrlv(intrlv_pattern_2D[j], x_tmp, Nv);
//								}
//								else {
//									randIntrlv(intrlv_pattern_info_mapping, x_tmp, Nv);
//								}
//							}
//							else {
//								if (!flag_separate_interleaving)
//								{
//									randIntrlv(non_uniform_intrlv_pattern, x_tmp, Nv);
//								}
//								else {
//									randIntrlv(non_uniform_intrlv_pattern_array[j], x_tmp, Nv);
//								}
//							}
//							for (int i = 0; i < Nv; i++) { x[i][j] = x_tmp[i]; }
//						}
//					}
//
//					/*cout << "u after interleaving" << endl;
//					display_array(x, Nv, Nh);*/
//
//					if (flag_interleave_block)
//					{
//						randIntrlvBlock(x, Nh, Nv, ModType, interleave_rows);
//					}
//				}
//				if (MOD_DIM == "Spatial")
//				{
//					for (int i = 0; i < Nh; i++)
//					{
//						// Extracting a single column
//						for (int j = 0; j < Nv; j++) { x_tmp[j] = x[j][i]; }
//						modulation(x_tmp, x_mod_tmp, ModType, Nt);
//						for (int j = 0; j < Nsym; j++) { x_mod[j][i] = x_mod_tmp[j]; }
//						// test of transmitted bits
//						/*if (i==0)
//							for (int j = 0; j < Nv; j++)
//							{
//								if (!(j % ModType)) cout << endl;
//								cout << x[j][i] << "  ";
//							}*/
//							// test of transmitted bits
//						/*for (int i = 0; i < mimo_blk_height; i++)
//						{
//							cout << x_mod_tmp[i] << endl;
//						}
//						cout << "111111" << endl;*/
//					}
//					// AWGN Channel (received y_mod )
//					if (!non_uniform_channel)H = generateChannelMatrix(Nr, Nt);
//					else H = generateChannelMatrixNonuniform(Nr, Nt, channel_scale_factor);
//					JacobiSVD<MatrixXcf> svd(H, ComputeThinU | ComputeThinV);
//					MatrixXcf U = svd.matrixU();
//					MatrixXcf V = svd.matrixV();
//					auto S = svd.singularValues();
//
//					for (int i = 0; i < Nh; i++)
//					{
//						// Extracting a single column
//						for (int j = 0; j < Nsym; j++) { x_mod_tmp[j] = x_mod[j][i]; }
//						AWGN(Nr, Nt, x_mod_tmp, y_mod_tmp, H, sigma_sq);
//						for (int j = 0; j < Nr; j++) { y_mod[j][i] = y_mod_tmp[j]; }
//					}
//					// Convert H and y to real domain
//					convertHToReal(H, H_real);
//					// TEST:H_real
//					/*cout << "H_real";
//					for (int i = 0; i < 2 * Nr; i++)
//					{
//						cout << endl;
//						for (int j = 0; j < 2 * Nt; j++)
//							printf("%.10f  ", H_real(i,j));
//					}*/
//					// TEST:H_real
//					for (int i = 0; i < Nh; i++)
//					{
//						for (int j = 0; j < Nr; j++) { y_mod_tmp[j] = y_mod[j][i]; }
//						convertYToReal(Nt, y_mod_tmp, y_mod_real_tmp);
//						for (int j = 0; j < 2 * Nr; j++) { y_mod_real[j][i] = y_mod_real_tmp[j]; /*cout << y_mod_real_tmp[j] << endl;*/ }
//					}
//					for (int i = 0; i < Nh; i++)
//					{
//						for (int j = 0; j < 2 * Nr; j++) { y_mod_real_temp_omp[i][j] = y_mod_real[j][i]; }
//					}
//					// for KBest Detection
//
//					//omp_set_num_threads(8);
//
//#pragma omp parallel for private(z)
//					for (int j = 0; j < Nh; j++)
//					{
//						MatrixXf received_signal_omp(2 * Nr, 1);
//						int k_kbest_extend = pow(2, ModType / 2) * K_KBEST;
//						float** paths_kbest = new float* [k_kbest_extend];
//						for (int i = 0; i < k_kbest_extend; i++) paths_kbest[i] = new float[2 * Nt];
//						float* ED = new float[K_KBEST];
//
//						if (MIMODET == "MMSE") { MMSEdetection(H_real, Identity, received_signal, mmse_matrix, output_signal, HTH, y_mod_real_temp_omp[j], y_mmse_temp_omp[j], sigma_sq); }
//						// ִ��һЩ����
//
//
//						if (MIMODET == "KBEST")
//
//						{
//							Eigen::HouseholderQR<Eigen::MatrixXf> qr_Hreal(H_real);
//							Q = qr_Hreal.householderQ();
//							R = qr_Hreal.matrixQR().triangularView<Eigen::Upper>();
//							for (int i = 0; i < 2 * Nr; i++) { received_signal_omp(i, 0) = y_mod_real_temp_omp[j][i]; }
//							z = Q.transpose() * received_signal_omp;
//							KBestDetection(R, z, H_real, y_mod_real_temp_omp[j], ED, y_mmse_temp_omp[j], paths_kbest, K_KBEST, ModType, Nt, Nr);
//							llr_kbest_bit(Nt, y_mmse_temp_omp[j], oriLLR_kbest[j], sigma_sq, ModType, K_KBEST, paths_kbest, ED);
//						}
//						for (int i_path = 0; i_path < k_kbest_extend; i_path++) { delete[] paths_kbest[i_path]; }
//
//						delete[] paths_kbest;
//						delete[] ED;
//					}
//
//					for (int i = 0; i < Nh; i++)
//					{
//						for (int j = 0; j < 2 * Nt; j++) { y_mmse[j][i] = y_mmse_temp_omp[i][j]; }
//					}
//					// sym llr to bitwise llr
//					// Detection is done for every column, that is, for every Nv bits
//					// Interleaving is not currently considered
//					for (int i = 0; i < Nh; i++)
//					{
//						for (int j = 0; j < 2 * Nt; j++) { y_mmse_tmp[j] = y_mmse[j][i]; }
//						if (MIMODET == "MMSE") { llr_mmse_bit(Nt, y_mmse_tmp, oriLLR_tmp, sigma_sq, ModType); }
//						if (MIMODET == "KBEST") { for (int j = 0; j < len_ori_llr_tmp; j++) { oriLLR_tmp[j] = oriLLR_kbest[i][j]; } }
//						if (flag_interleave)
//						{
//							randDeintrlv(intrlv_pattern, oriLLR_tmp, Nv);
//						}
//						bool flagIsMax = false;
//						for (int j = 0; j < Nv; j++)
//						{
//							if (!(j % 2))
//							{
//								if (abs(x_mmse_out[j][i]) > abs(x_mmse_out[j][i + 1])) { flagIsMax = true; }
//								else { flagIsMax = false; }
//							}
//							else
//							{
//								if (abs(x_mmse_out[j][i]) > abs(x_mmse_out[j][i - 1])) { flagIsMax = true; }
//								else { flagIsMax = false; }
//							}
//							oriLLRh[j][i] = oriLLR_tmp[j];
//							upLLR[j][i] = oriLLR_tmp[j];
//							// TEST�� det error
//							x_mmse_out[j][i] = 0;
//							if (oriLLR_tmp[j] < 0) { x_mmse_out[j][i] = 1; }
//							if (x_mmse_out[j][i] != x[j][i]) {
//								n_det_err++;
//								// error at max-reliable bit. Warning: Only fit for 16QAM modulation
//								if (flagIsMax) { n_det_err_mrb++; }
//							}
//							total_bit_cnt++;
//						}
//					}
//					// test : detection TIME
//					/*det_end = std::chrono::system_clock::now();
//					auto duration_detection = std::chrono::duration_cast<std::chrono::microseconds>(det_end - start);
//					file_det << duration_detection.count() << endl;*/
//				}
//				// Temporal modulation
//				else if (MOD_DIM == "Temporal")
//				{
//
//					//modulation
//					//cout << "START TEST" << endl;
//					for (int i = 0; i < Nv; i++) // every row
//					{
//						modulation(x[i], x_mod[i], ModType, mimo_blk_length);
//					}
//					/*for (int i = 0; i < mimo_blk_height; i++)
//					{
//						for (int j = 0; j < mimo_blk_length; j++)
//							cout << x_mod[i][j] << endl;
//					}
//					cout << "111111" << endl;*/
//					//AWGN channel (received y_mod)
//					// test: modulated symbols
//					/*cout << endl;
//					cout << "modulated symbols" << endl;
//					for (int i = 0; i < Nt; i++) cout << x_mod[i][0] << endl;
//					cout << "modulated symbols" << endl;*/
//					// test: modulated symbols
//					// H = generateChannelMatrix(Nr, Nt);
//					if (!non_uniform_channel)H = generateChannelMatrix(Nr, Nt);
//					else H = generateChannelMatrixNonuniform(Nr, Nt, channel_scale_factor);
//					JacobiSVD<MatrixXcf> svd(H, ComputeFullU | ComputeFullV);
//					MatrixXcf U = svd.matrixU();
//					// cout << U.rows() << " " << U.cols() << endl;
//					MatrixXcf V = svd.matrixV();
//					// cout << V.rows() << " " << V.cols() << endl;
//					auto S = svd.singularValues();
//					// if (nSNR == 10.0) { cout << "Temporal" << endl; }
//					if (mimo_svd_average_power_alloc)
//					{
//						float my_sum = 0.0;
//						if (K_time_array[0] == K_time_array[1]) {
//							for (int i = 0; i < Nt; i++) { singular_values_array[i] = S[i]; }
//							float sum_tmp = 0;
//							for (int i = 0; i < Nt; i++) {
//								sum_tmp += 1.0 / (singular_values_array[i] * singular_values_array[i]);
//							}
//							for (int i = 0; i < Nt; i++) {
//								float sing_tmp = 1.0 / (singular_values_array[i] * singular_values_array[i]);
//								antenna_power[i] = (float)Nt * sing_tmp / sum_tmp;
//								my_sum += antenna_power[i];
//								/*cout << antenna_power[i] << "  ";
//								cout << sqrt(antenna_power[i]) * singular_values_array[i] << "  ";
//								cout << endl;*/
//							}
//							// cout << "My_sum = " << my_sum << endl;
//						}
//						else {
//							for (int i = 0; i < Nt; i++) { singular_values_array[i] = S[i]; }
//							float sum_tmp = 0;
//							int n_tmp = border;
//							for (int i = 0; i < n_tmp; i++) {
//								sum_tmp += 1.0 / (singular_values_array[i] * singular_values_array[i]);
//							}
//							for (int i = 0; i < n_tmp; i++) {
//								float sing_tmp = 1.0 / (singular_values_array[i] * singular_values_array[i]);
//								antenna_power[i] = (float)n_tmp * sing_tmp / sum_tmp;
//								/*cout << antenna_power[i] << "  ";
//								cout << sqrt(antenna_power[i]) * singular_values_array[i] << "  ";
//								cout << endl;*/
//							}
//							sum_tmp = 0;
//							n_tmp = Nt - border;
//							for (int i = 0; i < n_tmp; i++) {
//								sum_tmp += 1.0 / (singular_values_array[i + border] * singular_values_array[i + border]);
//							}
//							for (int i = 0; i < n_tmp; i++) {
//								float sing_tmp = 1.0 / (singular_values_array[i + border] * singular_values_array[i + border]);
//								antenna_power[i + border] = (float)n_tmp * sing_tmp / sum_tmp;
//								/*cout << antenna_power[i] << "  ";
//								cout << sqrt(antenna_power[i]) * singular_values_array[i] << "  ";
//								cout << endl;*/
//							}
//
//						}
//						for (int i = 0; i < Nt; i++) { S[i] = S[i] * sqrt(antenna_power[i]); }
//					}
//					else if (mimo_waterfilling_power_alloc) {
//						for (int i = 0; i < Nt; i++) { singular_values_array[i] = S[i]; }
//						float sum_tmp = 0;
//						for (int i = 0; i < Nt; i++) {
//							sum_tmp += (singular_values_array[i]);
//						}
//						for (int i = 0; i < Nt; i++) {
//							float sing_tmp = singular_values_array[i];
//							antenna_power[i] = (float)Nt * sing_tmp / sum_tmp;
//							/*cout << antenna_power[i] << "  ";
//							cout << sqrt(antenna_power[i]) * singular_values_array[i] << "  ";
//							cout << endl;*/
//						}
//						for (int i = 0; i < Nt; i++) { S[i] = S[i] * sqrt(antenna_power[i]); }
//					}
//					else if (mimo_rate_matching_power_alloc) {
//						for (int i = 0; i < Nt; i++) { singular_values_array[i] = S[i]; }
//						float sum_tmp = 0;
//						/*for (int i = 0; i < Nt; i++) {
//							sum_tmp += (singular_values_array[i]);
//						}*/
//						sum_tmp = border * mimo_rate_matching_enlarge_ratio + (Nv - border) * 1.0;
//						for (int i = 0; i < Nt; i++) {
//							float power_tmp;
//							if (i > Nt - border - 1) { power_tmp = mimo_rate_matching_enlarge_ratio; }
//							else { power_tmp = 1.0; }
//							// float sing_tmp = singular_values_array[i];
//							antenna_power[i] = (float)Nt * power_tmp / sum_tmp;
//							/*cout << antenna_power[i] << "  ";
//							cout << sqrt(antenna_power[i]) * singular_values_array[i] << "  ";
//							cout << endl;*/
//						}
//						tmp = 0;
//						/*for (int i = 0; i < Nt; i++) { tmp += antenna_power[i]; }
//						cout << tmp << endl;*/
//						for (int i = 0; i < Nt; i++) { S[i] = S[i] * sqrt(antenna_power[i]); }
//					}
//					// for (int i = 0; i < Nt; i++) { cout << S[i] << "  "; } cout << endl;
//					// cout << S[0] << endl;
//					// SVD: H = USV^T; 
//					//MatrixXcf y_mat(Nr, Nsym);
//					//MatrixXcf UTy(Nr, Nsym);
//					//MatrixXcf x_mat(Nt, Nsym);
//					//MatrixXcf Vx(Nt, 1);
//					//MatrixXcf noise(Nr, Nsym);
//					//for (int i = 0; i < Nv; i++)
//					//{
//					//	for (int j = 0; j < Nsym; j++)
//					//	{
//					//		x_mat(i, j) = x_mod[i][j];
//					//		{
//					//			float n_real = (float)sqrt(sigma_sq / 2) * (float)sampleNormal_mimo();
//					//			float n_imag = (float)sqrt(sigma_sq / 2) * (float)sampleNormal_mimo();
//					//			noise(i, j) = std::complex<float>(n_real, n_imag);
//					//		}
//					//	}
//					//}
//					//Vx = V * x_mat;
//					//y_mat = H * Vx + noise; // y = H*V*x + n
//					//UTy = U.conjugate().transpose() * y_mat; // y_precode = U^H*y
//					//
//					//test SVD
//					//for (int i = 0; i < Nv; i++) { cout << x_mod[i][0] << "  "; x_vec(i, 0) = x_mod[i][0]; } //cout << x_vec(i, 0) << "  "; }
//					//cout << endl;
//					//Vx = V * x_vec;
//					//y_vec = H * Vx;
//					//UTy = U.transpose() * y_vec;
//					//test SVD
//// # pragma omp parallel for
//					for (int j = 0; j < Nsym; j++)
//					{
//						/*std::complex<float>* x_mod_tmp = new std::complex<float>[Nv];
//						std::complex<float>* x_precoded_tmp = new std::complex<float>[Nt];
//						std::complex<float>* y_mod_tmp = new std::complex<float>[Nr];
//						std::complex<float>* y_precoded_tmp = new std::complex<float>[Nr];*/
//						//// Extracting a single column
//						for (int i = 0; i < Nv; i++) { x_mod_tmp[i] = x_mod[i][j] * sqrt(antenna_power[i]); }
//						precodeSVD(Nr, Nt, x_mod_tmp, x_precoded_tmp, V); // precoding: x_precode = V*x
//						// for (int i = 0; i < Nt; i++) { x_precoded_tmp[i] = x_precoded_tmp[i] * sqrt(antenna_power[i]); }
//						AWGN(Nr, Nt, x_precoded_tmp, y_mod_tmp, H, sigma_sq);
//						precodeSVD_preprocessing(Nr, Nt, y_mod_tmp, y_precoded_tmp, U); // precoding: y_precode = U^H*y
//						for (int i = 0; i < Nr; i++) { y_mod[i][j] = y_precoded_tmp[i]; }
//						for (int i = 0; i < Nt; i++) { y_mod[i][j] = S[i] * y_mod[i][j] / (S[i] * S[i] + sigma_sq); } // MMSE detection with precoding
//						//for (int i = 0; i < Nt; i++) { y_mod[i][j] = S[i] * UTy(i,j) / (S[i] * S[i] + sigma_sq); } // MMSE detection with precoding
//						//delete[] x_mod_tmp; delete[] x_precoded_tmp; delete[] y_mod_tmp; delete[] y_precoded_tmp;
//					}
//					// for (int j = 0; j < Nt; j++) { std::cout << y_mod[j][0]*S[j](S[j]*S[j]+sigma_sq) << " "; } std::cout << std::endl;
//					// for (int j = 0; j < Nt; j++) { std::cout << UTy(j,0) / S[j] << " "; } std::cout << std::endl;
//					// for (int j = 0; j < Nt; j++) { std::cout << y_mod[j][0] << endl; }
//					// convertHToReal(H, H_real);
//
//
//					for (int j = 0; j < Nsym; j++)
//					{
//						for (int i = 0; i < Nr; i++) { y_mod_tmp[i] = y_mod[i][j]; }
//						convertYToReal(Nt, y_mod_tmp, y_mod_real_tmp);
//						for (int i = 0; i < 2 * Nr; i++) { y_mod_real[i][j] = y_mod_real_tmp[i]; }
//					}
//					//MMSE detection done in paralell
//					for (int i = 0; i < Nsym; i++)
//					{
//						for (int j = 0; j < 2 * Nr; j++) { y_mod_real_temp_omp[i][j] = y_mod_real[j][i]; }
//					}
//					//Eigen::HouseholderQR<Eigen::MatrixXf> qr_Hreal(H_real);
//					//Q = qr_Hreal.householderQ();
//					//R = qr_Hreal.matrixQR().triangularView<Eigen::Upper>();
//					//// for EP detection
//					//// copy H_real
//					//for (int i = 0; i < Nr2; i++)
//					//{
//					//	for (int j = 0; j < Nt2; j++)
//					//	{
//					//		Hreal_array[i * Nt2 + j] = H_real(i, j);
//					//	}
//					//}
//					// calculate HTH
//					// cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, Nt2, Nt2, Nr2, 1.0, Hreal_array, Nt2, Hreal_array, Nt2, 0.0, HTH_array, Nt2);
//					//cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, Nt2, 1, Nr2, 1.0, ch_mtx_struct.H_real[cnt], Nt2, data_struct.rx_buffer_real[cnt], 1, 0.0, ch_mtx_struct.hty, 1);
//
//					//omp_set_num_threads(8);
//					//cout << 1 << endl;
//					//omp_set_num_threads(4);
//					/*auto start = std::chrono::system_clock::now();*/
//					/*MatrixXf H_real_trans(2*Nt, 2*Nr);
//					for (int i = 0; i < 2 * Nt; i++)
//						for (int j = 0; j < 2 * Nr; j++)
//							H_real_trans(i, j) = H_real(j, i);
//					HTHIinv = (H_real_trans * H_real + Identity).inverse();*/
//					//# pragma omp parallel for private(z)
//					//					for (int j = 0; j < Nsym; j++)
//					//
//					//						//for (int x=0; x<Nsym; x++)
//					//					{
//					//						// int j = x % Nsym;
//					//						MatrixXf received_signal_omp(2 * Nr, 1);
//					//						int k_kbest_extend = pow(2, ModType / 2) * K_KBEST;
//					//						float** paths_kbest = new float* [k_kbest_extend];
//					//						for (int i = 0; i < k_kbest_extend; i++) { paths_kbest[i] = new float[2 * Nt]; }
//					//						float* ED = new float[K_KBEST];
//					//						float* HTY_array_tmp = new float[2 * Nt];
//					//						// Declared arrays for EP
//					//						int flag = 1;
//					//						int temp = 0;
//					//						//float* extr_mean = new float[2 * Nr];
//					//						//float* extr_var = new float[2 * Nr];
//					//						double* vec_alpha_in = new double[2 * Nr];
//					//						double* vec_gamma_in = new double[2 * Nr];
//					//						double* p_symbol_rearrange = new double[Nt2 * pow(2, ModType / 2)];
//					//						for (int i = 0; i < Nt2 * pow(2, ModType / 2); i++) { p_symbol_rearrange[i] = 1 / (pow(2, ModType / 2)); }
//					//						for (int i = 0; i < Nr2; i++)
//					//						{
//					//							//extr_mean[i] = 0;
//					//							//extr_var[i] = 0;
//					//							vec_alpha_in[i] = 2;
//					//							vec_gamma_in[i] = 0;
//					//						}
//					//
//					//						if (MIMODET == "MMSE")
//					//						{
//					//							MMSEdetection(H_real, Identity, received_signal, mmse_matrix, output_signal, HTH, y_mod_real_temp_omp[j], y_mmse_temp_omp[j], sigma_sq);
//					//							//MMSEdetection_noInverse(H_real, Identity, received_signal, mmse_matrix, output_signal, HTH, HTHIinv, y_mod_real_temp_omp[j], y_mmse_temp_omp[j], sigma_sq);
//					//						}
//					//						if (MIMODET == "KBEST")
//					//						{
//					//							for (int i = 0; i < 2 * Nr; i++) { received_signal_omp(i, 0) = y_mod_real_temp_omp[j][i]; }
//					//							z = Q.transpose() * received_signal_omp;
//					//							KBestDetection(R, z, H_real, y_mod_real_temp_omp[j], ED, y_mmse_temp_omp[j], paths_kbest, K_KBEST, ModType, Nt, Nr);
//					//							llr_kbest_bit(Nt, y_mmse_temp_omp[j], oriLLR_kbest[j], sigma_sq, ModType, K_KBEST, paths_kbest, ED);
//					//						}
//					//						if (MIMODET == "EP")
//					//						{
//					//							// test: save H_real, recv_signal
//					//							/*save_mat_txt(y_mod_real_temp_omp[j], Nr2, "RxSig0128.txt");
//					//							save_mat_txt(H_real, "Hreal0128.txt");*/
//					//							// test: save H_real, recv_signal
//					//							cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, Nt2, 1, Nr2, 1.0, Hreal_array, Nt2, y_mod_real_temp_omp[j], 1, 0.0, HTY_array_tmp, 1);
//					//							//test:HTH_array
//					//							/*for (int i = 0; i < Nt2; i++)
//					//								for (int s = 0; s < Nt2; s++)
//					//								{
//					//									if (s == 0) std::cout << endl;
//					//									std::cout << setw(10) << HTH_array[i * Nt2 + s];
//					//								}
//					//							cout << endl;*/
//					//							// TEST:HTY
//					//							/*for (int i = 0; i < Nr2; i++)
//					//								received_signal_omp(i, 0) = y_mod_real_temp_omp[j][i];
//					//							MatrixXf HTransY = H_real.transpose() * received_signal_omp;
//					//							cout << "Eigen" << endl;
//					//							for (int i = 0; i < Nt2; i++) cout << setw(10)<<HTransY(i, 0);
//					//							cout << "mkl" << endl;
//					//							for (int i = 0; i < Nt2; i++) cout << HTY_array[i];*/
//					//							// TEST:HTY
//					//							EPD_detection_float(flag, Nt, Nr, iter_MIMO, ModType, 1.0, sigma_sq, HTH_array, HTY_array_tmp, temp, ep_extr_mean_omp[j], ep_extr_var_omp[j], vec_alpha_in, vec_gamma_in, p_symbol_rearrange);
//					//						}
//					//						//test EP detection
//					//						/*cout << "TEST:EP DETECTION" << endl;
//					//						for (int i = 0; i < Nr2; i++) std::cout << extr_mean[i]<< endl;
//					//						cout << "TEST:EP DETECTION" << endl;*/
//					//						//test EP detection
//					//						for (int i_path = 0; i_path < k_kbest_extend; i_path++) { delete[] paths_kbest[i_path]; }
//					//						delete[] paths_kbest;
//					//						delete[] ED;
//					//						//delete[] extr_mean;
//					//						//delete[] extr_var;
//					//						delete[] vec_alpha_in;
//					//						delete[] vec_gamma_in;
//					//						delete[] p_symbol_rearrange;
//					//						delete[] HTY_array_tmp;
//					//					}
//										// test: detection time
//										/*auto detection_over = std::chrono::system_clock::now();
//										auto duration_detection = std::chrono::duration_cast<std::chrono::microseconds>(detection_over - start);
//										cout << endl;
//										cout <<"Detection duration" << duration_detection.count() << "us" << endl;*/
//										// test: detection time
//					for (int i = 0; i < Nsym; i++)
//					{
//						for (int j = 0; j < 2 * Nt; j++) { y_mmse[j][i] = y_mmse_temp_omp[i][j]; }
//					}
//					//To bit domain output
//					for (int j = 0; j < Nsym; j++)
//					{
//						for (int i = 0; i < 2 * Nt; i++) { y_mmse_tmp[i] = y_mod_real[i][j]; }
//						// With Spatial Modulation, Nt = Nv
//						if (MIMODET == "MMSE") { llr_mmse_bit(Nt, y_mmse_tmp, oriLLR_tmp, sigma_sq, ModType); }
//						if (MIMODET == "EP") { llr_ep_bit(Nt, ep_extr_mean_omp[j], ep_extr_var_omp[j], oriLLR_tmp, ModType); }
//						if (MIMODET == "KBEST") { for (int i = 0; i < len_ori_llr_tmp; i++) { oriLLR_tmp[i] = oriLLR_kbest[j][i]; } }
//						int i_bit = 0;
//						int bit_col = 0;
//						for (int j_inner = 0; j_inner < Nv; j_inner++)
//						{
//							for (int i_bit_col = 0; i_bit_col < ModType; i_bit_col++)
//							{
//								i_bit = ModType * j_inner + i_bit_col;
//								bit_col = j * ModType + i_bit_col;
//								oriLLRh[j_inner][bit_col] = oriLLR_tmp[i_bit];
//								upLLR[j_inner][bit_col] = oriLLR_tmp[i_bit];
//							}
//
//						}
//					}
//					// test: LLR transfer time
//					/*auto llr_to_bit_over = std::chrono::system_clock::now();
//					auto duration_llr = std::chrono::duration_cast<std::chrono::microseconds>(llr_to_bit_over - detection_over);
//					cout << endl;
//					cout<< "LLR transfer duration" << duration_llr.count() << "us" << endl;*/
//					// test: LLR transfer time
//					if (flag_interleave)
//					{
//						if ((!flag_interleave_block) && flag_interleave_spatial)
//						{
//							for (int j = 0; j < Nh; j++)
//							{
//								for (int i = 0; i < Nv; i++) { oriLLR_tmp[i] = oriLLRh[i][j]; }
//								if (K_time_array[0] == K_time_array[1]) {
//									if (!flag_information_mapping)  randDeintrlv(intrlv_pattern_2D[j], oriLLR_tmp, Nv);
//									else {
//										randDeintrlv(intrlv_pattern_info_mapping, oriLLR_tmp, Nv);
//									}
//								}
//								else {
//									if (!flag_separate_interleaving) { randDeintrlv(non_uniform_intrlv_pattern, oriLLR_tmp, Nv); }
//									else { randDeintrlv(non_uniform_intrlv_pattern_array[j], oriLLR_tmp, Nv); }
//								}
//								for (int i = 0; i < Nv; i++)
//								{
//									oriLLRh[i][j] = oriLLR_tmp[i];
//									upLLR[i][j] = oriLLR_tmp[i];
//								}
//
//							}
//						}
//
//						if (flag_interleave_block)
//						{
//							randDeIntrlvBlock(oriLLRh, Nh, Nv, ModType, interleave_rows);
//							randDeIntrlvBlock(upLLR, Nh, Nv, ModType, interleave_rows);
//						}
//
//						if (flag_interleave_all) // Deinterleave in time domain
//						{
//							for (int j = 0; j < Nv; j++)
//							{
//								randDeintrlv(intrlv_pattern_2D_time[j], oriLLRh[j], Nh);
//								randDeintrlv(intrlv_pattern_2D_time[j], upLLR[j], Nh);
//							}
//						}
//
//						if (flag_interleave_partial)
//						{
//							for (int j = 0; j < Nv; j++)
//							{
//								for (int k = 0; k < intrlv_pattern_num; k++)
//								{
//									randDeintrlv(intrlv_pattern_2D_time_partial[j][k], oriLLRh[j] + k * ModType, ModType);
//									randDeintrlv(intrlv_pattern_2D_time_partial[j][k], upLLR[j] + k * ModType, ModType);
//								}
//							}
//						}
//
//					}
//
//					for (int i = 0; i < Nv; i++)
//					{
//						for (int j = 0; j < Nh; j++)
//						{
//							total_bit_cnt++;
//							int bit_dec = (oriLLRh[i][j] < 0);
//							if (bit_dec != x[i][j]) n_det_err++;
//						}
//					}
//
//				}
//				// test: detection errors
//
//				// test: detection errors
//			}
//			if (atten == 1)
//			{
//				for (int i = 0; i < Nv; i++)
//					for (int j = 0; j < Nh; j++)
//					{
//						// 0-1; 1-(-1)
//						y[i][j] = (1 - 2 * x[i][j]) + (float)sqrt(sigma_sq) * (float)sampleNormal();
//						oriLLRh[i][j] = 2 * y[i][j] / sigma_sq;
//						//oriLLRv[j][i] = oriLLRh[i][j];
//						upLLR[i][j] = oriLLRh[i][j];
//					}
//			}
//			//Polar Decoding
//			for (int iter = 1; iter <= softiter; iter++)
//			{
//				//vertical decoding
//		   // test: decoding time
//			//auto decoding_start = std::chrono::system_clock::now();
//
//			// test: decoding time
//#pragma omp parallel for
//				for (int decodingi = 0; decodingi < Nh; decodingi++) // 291.70kb
//				{
//
//					//for SCL Decoding
//					int Counth = 0, Countinfoh = 0;
//					float Qsumunh = 0;
//					int** sum = new int* [L];
//					//for (int i = 0; i < L; i++)  sum[i] = new  int[Nmax];
//					int** u1 = new int* [L];
//					for (int i = 0; i < L; i++)  u1[i] = new _MM_ALIGN16 int[Nmax];
//					for (int i = 0; i < L; i++)  sum[i] = new _MM_ALIGN16 int[Nmax];
//					//for (int i = 0; i < L; i++)  u1[i] = new  int[Nmax];
//
//					float* LLR_in = new float[L];
//					float* W = new float[2 * L];
//					int* SCL_index = new int[L];
//					int* better = new int[L];
//					int* worse = new int[L];
//					int* indexpm = new int[L];
//					float* PM = new float[L];
//					float** LLR = new float* [L];
//					for (int i = 0; i < L; i++)  LLR[i] = new _MM_ALIGN16 float[2 * Nmax + 2];
//					// for (int i = 0; i < L; i++)  LLR[i] = new float[2 * Nmax + 2];
//					////////////////////////////////
//					for (int j = 0; j < Nv; j++)
//						for (int k = 0; k < L; k++)
//							LLR[k][j] = upLLR[j][decodingi];
//					//cout << LLR[0][j] << " ";
//
//				//cout << endl;
//					int index = 0;
//					/*			if (L == 1)
//								{
//									//T_start = clock();
//									Count = 0;
//									decode(*LLR, *sum, (int)(log(Nv) / log(2)), A_Acv, flag, Nv, Kv);
//									//PolarEncode_xor(*u1, *sum, N);
//									//T_finish = clock();
//								}
//								else
//								{*/
//					for (int k = 0; k < L; k++) { PM[k] = 0; indexpm[k] = k; }
//					//T_start = clock();
//					decode_list(LLR, sum, (int)(log(Nv) / log(2)), A_Acv, 0, Nv, Kv, PM, L, (int)(log(L) / log(2)), LLR_in, W, SCL_index, better, worse, &Counth, &Countinfoh, &Qsumunh);
//					//T_finish = clock();
//					//Decoding_Duration += (T_finish - T_start);
//					for (int i = 0; i < L - 1; i++)
//						for (int j = i + 1; j < L; j++)
//							if (PM[indexpm[i]] < PM[indexpm[j]])
//							{
//								int ttt = indexpm[i];
//								indexpm[i] = indexpm[j];
//								indexpm[j] = ttt;
//							}
//					//PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), N);
//					index = indexpm[0];
//					//}
//					//calc distance s()
//					if (adaptive_beta)
//					{
//						if (softmethod == "Pyndiah") Pyndiahv_adaptivebeta(LLR, sum, oriLLRh, indexpm, iter * 2 - 2, decodingi, upLLR);
//						if (softmethod == "MITSO") MITSOv(LLR, PM, sum, oriLLRh, indexpm, iter * 2 - 2, decodingi, upLLR, Qsumunh);
//					}
//					// iter * 2 - 2
//					else
//					{
//						//cout << iter << endl;
//						if (softmethod == "Pyndiah") Pyndiahv(LLR, sum, oriLLRh, indexpm, iter, decodingi, upLLR);
//						if (softmethod == "MITSO") MITSOv(LLR, PM, sum, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunh);
//					}
//					// test: oriLLR and upLLR
//					/*cout << endl;
//					for (int i = 0; i < Nv; i++)
//						for (int j = 0; j < Nh; j++)
//							cout << oriLLRh[i][j] - upLLR[i][j] << "   ";
//					cout << endl;*/
//					/*cout << oriLLRh[0][0] << endl;
//					cout << upLLR[0][0] << endl;*/
//					// test: oriLLR and upLLR
//					// half iteration Enabled
//					// May Need Modification
//					if (half_iter && iter == softiter)
//					{
//						if (iter == softiter) PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), Nh);
//						for (int j = 0; j < Kv; j++)
//							umid[Av[j]][decodingi] = u1[index][Av[j]];
//					}
//					// half iteration
//					////////delete SCL variables
//					for (int i = 0; i < L; i++)
//						delete[]LLR[i];
//					delete[] LLR;
//					for (int i = 0; i < L; i++)
//						delete[]sum[i];
//					for (int i = 0; i < L; i++)
//						delete[]u1[i];
//					delete[] u1;
//					delete[] sum;
//					delete[] PM;
//					delete[] LLR_in;
//					delete[] W;
//					delete[] SCL_index;
//					delete[] better;
//					delete[] worse;
//					delete[] indexpm;
//
//				}
//				//cout << "Horizon Start:\n";
//				//horizontal decoding
//				// #
//				if (weirditer) { continue; }
//
//				// test: LLR after horizontal decoding
//				/*cout << "After horizontal decoding" << endl;
//				for (int i = 0; i < Nv; i++)
//				{
//					for (int j = 0; j < Nh; j++)
//					{
//						cout << (oriLLRh[i][j] < 0) << " ";
//					}
//					cout << endl;
//				}*/
//				// test: LLR after horizontal decoding
//				// inner deinterleaving
//				if (flag_interleave_inner)
//				{
//					float* tmp_llr_to_interleave = new float[Nv];
//					for (int i = 0; i < Nh; i++)
//					{
//						for (int j = 0; j < Nv; j++) { tmp_llr_to_interleave[j] = oriLLRh[j][i]; }
//						randDeIntrlvPartial(intrlv_pattern_2D_inner[i], tmp_llr_to_interleave, Av, Nv, Kv);
//						for (int j = 0; j < Nv; j++) { oriLLRh[j][i] = tmp_llr_to_interleave[j]; }
//					}
//					for (int i = 0; i < Nh; i++)
//					{
//						for (int j = 0; j < Nv; j++) { tmp_llr_to_interleave[j] = upLLR[j][i]; }
//						randDeIntrlvPartial(intrlv_pattern_2D_inner[i], tmp_llr_to_interleave, Av, Nv, Kv);
//						for (int j = 0; j < Nv; j++) { upLLR[j][i] = tmp_llr_to_interleave[j]; }
//					}
//					delete[] tmp_llr_to_interleave;
//				}
//
//
//				/*auto decoding_over = std::chrono::system_clock::now();
//				auto duration_decoding = std::chrono::duration_cast<std::chrono::microseconds>(decoding_over - decoding_start);
//				cout << endl;
//				cout << "Decoding duration" << duration_decoding.count() << "us" << endl;*/
//
//#pragma omp parallel for
//				for (int decodingi = 0; decodingi < Nv; decodingi++)
//				{
//					// 293.70kb
//					// cout << "KKKK" << endl;
//					if (half_iter && (iter == softiter)) { break; }
//					string polarde_now = polarde_array[iter - 1];
//					//for SCL Decoding
//					int Countv = 0, Countinfov = 0;
//					float Qsumunv = 0;
//					int** u1 = new int* [L];
//					for (int i = 0; i < L; i++)  u1[i] = new _MM_ALIGN16 int[Nmax];
//					int** u1_LSD = new int* [2 * L_LSD];
//					for (int i = 0; i < 2 * L_LSD; i++) u1_LSD[i] = new int[Nmax];
//					int** sum = new int* [L];
//					for (int i = 0; i < L; i++)  sum[i] = new _MM_ALIGN16 int[Nmax];
//					int** sum_LSD = new int* [2 * L_LSD];
//					for (int i = 0; i < 2 * L_LSD; i++)  sum_LSD[i] = new int[Nmax];
//					float* LLR_in = new float[L];
//					int* u_decod = new int[Nh];
//					int* n_paths = new int[Nh];
//					float* W = new float[2 * L];
//					int* SCL_index = new int[L];
//					int* better = new int[L];
//					int* worse = new int[L];
//					int* indexpm = new int[L_LSD];
//					float* PM = new float[L];
//					float** LLR = new float* [L];
//					for (int i = 0; i < L; i++)  LLR[i] = new _MM_ALIGN16 float[2 * Nmax + 2];
//
//					for (int i = 0; i < L_LSD; i++) { indexpm[i] = i; }
//					////////////////////////////////
//					for (int j = 0; j < Nh; j++)
//						for (int k = 0; k < L; k++)
//							LLR[k][j] = upLLR[decodingi][j];
//					//	cout << LLR[0][j] << " ";
//					calc_npaths(A_Ach, Nh, n_paths, L_LSD);
//					//	cout << endl;
//					int index = 0;
//
//					/*			if (L == 1)
//								{
//									//T_start = clock();
//									Count = 0;
//									decode(*LLR, *sum, (int)(log(Nh) / log(2)), A_Ach, flag, Nh, Kh);
//									if (iter == softiter) PolarEncode_xor(*u1, *sum, Nh);
//									//T_finish = clock();
//								}
//								else
//								{*/
//								//for (int k = 0; k < L; k++)  PM[k] = 0, indexpm[k] = k; 
//								////T_start = clock();
//								//decode_list(LLR, sum, (int)(log(Nh) / log(2)), A_Ach, 0, Nh, Kh, PM, L, (int)(log(L) / log(2)), LLR_in, W, SCL_index, better, worse, &Countv, &Countinfov, &Qsumunv);
//								////T_finish = clock();
//								////Decoding_Duration += (T_finish - T_start);
//								//for (int i = 0; i < L - 1; i++)
//								//	for (int j = i + 1; j < L; j++)
//								//		if (PM[indexpm[i]] < PM[indexpm[j]])
//								//		{
//								//			int ttt = indexpm[i];
//								//			indexpm[i] = indexpm[j];
//								//			indexpm[j] = ttt;
//								//		}
//
//
//							/*
//							* Encoding: Info -> Horizontal Encoding -> A -> set to 0 -> B -> Horizontal Encoding -> Vertical Encoding
//							* Decoding: Received Info -> Horizontal Decoding -> Vertical Decoding -> Horizontal Decoding
//							*/
//					if (polarde_now == "SLD")
//					{
//						LSD_decode_outall(u_decod, u1_LSD, LLR[0], GGh, A_Ach, Nh, L_LSD, n_paths);
//						// if (iter == softiter) PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), Nh);
//						// for soft information exchange, inverse coding:
//						// THIS LINE is WRONG! BUT IT SEEMS BETTER THAN THE RIGHT VERSION
//						for (int i = 0; i < L_LSD; i++) { PolarEncode_xor(sum_LSD[i], u1_LSD[i], Nh); }
//						// for (int i = 0; i < L_LSD; i++) { PolarEncode_xor(sum_LSD[i], u1_LSD[i], Nh); }
//						// if (iter == softiter) PolarEncode_xor(u1[0], sum_LSD[0], Nh);
//						index = indexpm[0];
//						// test: after horizontal decoding
//						/*if (decodingi == Av[0])
//						{
//							cout << "Direct output of horizontal decoding:" << endl;
//							for (int i = 0; i < Nh; i++) { cout << u1_LSD[0][i] << " "; }
//							cout << endl;
//							cout << "horizontal decoder out, encode once:" << endl;
//							for (int i = 0; i < Nh; i++) { cout << sum_LSD[0][i] << " "; }
//							cout << endl;
//						}*/
//						// TEST
//						/*cout << "SCL decoded bits" << endl;
//						for (int i = 0; i < Nh; i++) { cout << sum[indexpm[0]][i] << "  "; }
//						cout << endl;
//						cout << "LSD decoded bits" << endl;
//						for (int i = 0; i < Nh; i++) { cout << u1[0][i] << "  "; }
//						cout << endl;*/
//						// TEST
//						//}
//						//calc distance s()
//						// iter*2-1
//						if (iter < softiter)
//						{
//							if (adaptive_beta)
//							{
//								if (softmethod == "Pyndiah") Pyndiahh_LSD_adaptivebeta(LLR, sum_LSD, oriLLRh, indexpm, iter, decodingi, upLLR);
//								if (softmethod == "MITSO") MITSOh(LLR, PM, sum_LSD, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunv);
//							}
//							else if (irregular_iteration)
//							{
//
//							}
//							else
//							{
//								// cout << "IN LSD!!!" << endl;
//								if (softmethod == "Pyndiah")
//									Pyndiahh_LSD_getsum(LLR, sum_LSD, oriLLRh, indexpm, iter, decodingi, upLLR, sum_sdis_array[decodingi], GGh, Ach);
//								//if (softmethod == "Pyndiah") Pyndiahh_LSD(LLR, sum_LSD, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR);
//								if (softmethod == "MITSO") MITSOh(LLR, PM, sum_LSD, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunv);
//							}
//							//	getchar();
//						}
//						else
//							for (int j = 0; j < Kh; j++)
//								umid[decodingi][Ah[j]] = u1_LSD[index][Ah[j]];
//					}
//					if (polarde_now == "SCL")
//					{
//						for (int k = 0; k < L; k++)  PM[k] = 0, indexpm[k] = k;
//						//T_start = clock();
//						// decode_list(LLR, sum, (int)(log(Nh) / log(2)), A_Ach, 0, Nh, Kh, PM, L, (int)(log(L) / log(2)), LLR_in, W, SCL_index, better, worse, &Countv, &Countinfov, &Qsumunv);
//						decode_list(LLR, sum, (int)(log(Nh) / log(2)), A_Ac_time_all[decodingi], 0, Nh, K_time_all[decodingi], PM, L, (int)(log(L) / log(2)), LLR_in, W, SCL_index, better, worse, &Countv, &Countinfov, &Qsumunv);
//						//T_finish = clock();
//						//Decoding_Duration += (T_finish - T_start);
//						for (int i = 0; i < L - 1; i++)
//							for (int j = i + 1; j < L; j++)
//								if (PM[indexpm[i]] < PM[indexpm[j]])
//								{
//									int ttt = indexpm[i];
//									indexpm[i] = indexpm[j];
//									indexpm[j] = ttt;
//								}
//						if (iter == softiter) PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), Nh);
//						// since systematic encoding is used in space domain, we directly copy the best path (x) to u1 and encode on time domain when all temporal decoing is done
//						/*if (iter == softiter)
//						{
//							for (int j = 0; j < Nh; j++) { u1[indexpm[0]][j] = sum[indexpm[0]][j]; }
//						}*/
//						index = indexpm[0];
//						//}
//						//calc distance s()
//						if (iter < softiter)
//						{
//							if (adaptive_beta)
//							{
//								if (softmethod == "Pyndiah") Pyndiahh_adaptivebeta(LLR, sum, oriLLRh, indexpm, iter, decodingi, upLLR);
//								if (softmethod == "MITSO") MITSOh(LLR, PM, sum, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunv);
//								//	getchar();
//							}
//							else if (irregular_coding)
//							{
//								if (softmethod == "Pyndiah") Pyndiahh_irregular(LLR, sum, oriLLRh, indexpm, iter, decodingi, upLLR, sum_sdis_array[decodingi], use_short_codeword);
//								if (softmethod == "MITSO") MITSOh_irregular(LLR, PM, sum, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunv, use_short_codeword);
//							}
//							else
//							{
//								//if (softmethod == "Pyndiah") Pyndiahh(LLR, sum, oriLLRh, indexpm, iter * 2 - 1, decodingi, upLLR);
//								if (softmethod == "Pyndiah") Pyndiahh_getsum(LLR, sum, oriLLRh, indexpm, iter, decodingi, upLLR, sum_sdis_array[decodingi]);
//								if (softmethod == "MITSO") MITSOh(LLR, PM, sum, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunv);
//								//	getchar();
//							}
//						}
//						else
//							for (int j = 0; j < Nh; j++)
//								umid[decodingi][j] = u1[index][j];
//					}
//
//					////////delete SCL variables
//					for (int i = 0; i < L; i++)
//						delete[]u1[i];
//					delete[] u1;
//					for (int i = 0; i < 2 * L_LSD; i++)
//						delete[]u1_LSD[i];
//					delete[] u1_LSD;
//					for (int i = 0; i < L; i++)
//						delete[]LLR[i];
//					delete[] LLR;
//					for (int i = 0; i < L; i++)
//						delete[]sum[i];
//					delete[] sum;
//					for (int i = 0; i < 2 * L_LSD; i++)
//						delete[]sum_LSD[i];
//					delete[] sum_LSD;
//					delete[] PM;
//					delete[] LLR_in;
//					delete[] W;
//					delete[] SCL_index;
//					delete[] better;
//					delete[] worse;
//					delete[] indexpm;
//					delete[] u_decod;
//					delete[] n_paths;
//				}
//
//
//				if (flag_interleave_inner)
//				{
//					float* tmp_llr_to_interleave = new float[Nv];
//					for (int i = 0; i < Nh; i++)
//					{
//						for (int j = 0; j < Nv; j++) { tmp_llr_to_interleave[j] = oriLLRh[j][i]; }
//						randIntrlvPartial(intrlv_pattern_2D_inner[i], tmp_llr_to_interleave, Av, Nv, Kv);
//						for (int j = 0; j < Nv; j++) { oriLLRh[j][i] = tmp_llr_to_interleave[j]; }
//					}
//					for (int i = 0; i < Nh; i++)
//					{
//						for (int j = 0; j < Nv; j++) { tmp_llr_to_interleave[j] = upLLR[j][i]; }
//						randIntrlvPartial(intrlv_pattern_2D_inner[i], tmp_llr_to_interleave, Av, Nv, Kv);
//						for (int j = 0; j < Nv; j++) { upLLR[j][i] = tmp_llr_to_interleave[j]; }
//					}
//					delete[] tmp_llr_to_interleave;
//				}
//			}
//			int calcsum[Nv], calcu1[Nv];
//			int calcu1h[Nh];
//			// calcsum is from the direct output of the decoder
//			// if output is the result of row decoding (meaning there is another column decoding to be done)
//			if (half_iter)
//			{
//				for (int i = 0; i < Kv; i++)
//				{
//					PolarEncode_xor(calcu1h, umid[Av[i]], Nh);
//					for (int j = 0; j < Kh; j++)
//						uout[Av[i]][Ah[j]] = calcu1h[Ah[j]];
//				}
//			}
//			// if output is the result of column decoding (meaning that there is another row decoding to be done)
//			else
//			{
//				for (int i = 0; i < Kh; i++)
//				{
//					for (int j = 0; j < Nv; j++)
//						calcsum[j] = umid[j][Ah[i]];
//					PolarEncode_xor(calcu1, calcsum, Nv);
//					for (int j = 0; j < Kv; j++)
//						uout[Av[j]][Ah[i]] = calcu1[Av[j]];
//				}
//			}
//			// test: final output
//			/*cout << "final output" << endl;
//			for (int k = 0; k < Nh; k++) { cout << uout[Av[0]][k] << " "; }
//			cout << endl;*/
//			// test: final output
//			//now without CRC
//			int Temp_Error = Num_Error;
//			for (int i = 0; i < Kv; i++)
//			{
//				//if ((i == 24) || (i == 56)) { continue; }
//				for (int j = 0; j < K_time_all[Av[i]]; j++)
//				{
//					//printf("%d ", (int)uout[Av[i]][Ah[j]]);
//					Num_Error += abs(u[Av[i]][A_time_all[Av[i]][j]] - umid[Av[i]][A_time_all[Av[i]][j]]);
//					if (abs(u[Av[i]][A_time_all[Av[i]][j]] - umid[Av[i]][A_time_all[Av[i]][j]]) == 1)
//					{
//						all_error_positions[Av[i]][A_time_all[Av[i]][j]] += 1;
//					}
//				}
//				//cout << endl;
//			}
//			//getchar();
//			if (Temp_Error != Num_Error) Num_Frame_Error++;
//			if (frame % 10 == 0)
//			{
//				/*printf("\r<%d> SNR= %.2f FER= %d / %d = %f BER= %d / %d = %lf E-Det = %f, E-MRB = %f"  ,
//					frame, nSNR, Num_Frame_Error, frame, (double)Num_Frame_Error / frame, Num_Error, (N_Info * frame), (double)Num_Error / N_Info / frame,
//					(double)n_det_err/total_bit_cnt, (double)n_det_err_mrb/(total_bit_cnt/(ModType/2)));*/
//				int L_sdis_sum = L_LSD;
//				sdis_sum = 0;
//				for (int i = 0; i < Nv; i++) sdis_sum += sum_sdis_array[i];
//				if (polarde_array[0] == "SCL") L_sdis_sum = L;
//				// cout << "sdis_sum_avg" << sdis_sum / frame / L_sdis_sum << endl;
//				printf("\r<%d> SNR= %.2f FER= %d / %d = %f BER= %d / %d = %lf ,sdis_avg = %f, E-MRB = %f",
//					frame, nSNR, Num_Frame_Error, frame, (double)Num_Frame_Error / frame, Num_Error, (N_Info * frame), (double)Num_Error / N_Info / frame, (double)sdis_sum / frame / L_sdis_sum, (double)n_det_err / total_bit_cnt);
//				fflush(stdout);
//			}
//
//
//			// output interleaving pattern
//			//if (Num_Frame_Error % record_frames == 0 && Num_Frame_Error > 0 && Num_Frame_Error!=n_frame_error_pre)
//			//{
//			//	n_frame_error_pre = Num_Frame_Error;
//			//	cout << endl;
//			//	int n_file = Num_Frame_Error / record_frames;
//			//	frame_delta = frame - frame_tmp;
//			//	cout << "Delta Frames: " << frame_delta << endl;
//			//	frame_tmp = frame;
//			//	string ofile_intrlv_name = string("intrlv_pattern_") + to_string(n_file) + string("_frame_") + to_string(frame_delta) + string(".txt");
//			//	string ofile_spatial_intrlv_name = string("intrlv_pattern_spatial") + to_string(n_file) + string(".txt");
//			//	ofstream ofile_intrlv(ofile_intrlv_name);
//			//	ofstream ofile_spatial_intrlv(ofile_spatial_intrlv_name);
//			//	for (int i = 0; i < Nv; i++)
//			//	{
//			//		for (int j = 0; j < intrlv_pattern_num; j++)
//			//		{
//
//			//			for (int k = 0; k < ModType; k++)
//			//				ofile_intrlv << intrlv_pattern_2D_time_partial[i][j][k] << " ";
//			//			//cout << endl;
//			//		}
//			//	}
//			//	ofile_intrlv.close();
//
//			//	for (int i = 0; i < Nh; i++)
//			//	{
//			//		for (int j = 0; j < Nv; j++)
//			//			ofile_spatial_intrlv << intrlv_pattern_2D[i][j] << " ";
//			//	}
//			//	ofile_spatial_intrlv.close();
//
//			//	// re-generate interleaving pattern
//			//	generateIntrlvPattern(intrlv_pattern, Nv);
//			//	for (int i = 0; i < Nh; i++) { generateIntrlvPattern(intrlv_pattern_2D[i], Nv); }
//			//	for (int i = 0; i < Nv; i++) { generateIntrlvPattern(intrlv_pattern_2D_time[i], Nh); }
//			//	if (flag_interleave_partial)
//			//	{
//			//		for (int i = 0; i < Nv; i++)
//			//		{
//			//			for (int j = 0; j < intrlv_pattern_num; j++)
//			//			{
//			//				generateIntrlvPattern(intrlv_pattern_2D_time_partial[i][j], ModType);
//			//				/*for (int k = 0; k < ModType; k++)
//			//					cout << intrlv_pattern_2D_time_partial[i][j][k];
//			//				cout << endl;*/
//			//			}
//			//		}
//			//	}
//			//	//Num_Frame_Error = 0;
//			//}
//
//
//			//printf("%d %d", Num_Frame_Error, Num_Error);
//			//getchar();
//			//Decoding_Duration += (T_finish - T_start);
//			/*auto dec_end = std::chrono::system_clock::now();
//			auto duration_decoding = std::chrono::duration_cast<std::chrono::microseconds>(dec_end - det_end);
//			auto duration_total = std::chrono::duration_cast<std::chrono::microseconds>(dec_end - start);
//			file_dec << duration_decoding.count() << endl;
//			file_total << duration_total.count() << endl;*/
//			if (frame % 16 == 0)
//			{
//				file_err_num << Num_Error << endl;
//			}
//			if (frame % 2000 == 0)
//			{
//				std::ofstream ofile_err("all_errors_log.txt");
//				for (int i = 0; i < Nv; i++)
//				{
//					for (int j = 0; j < Nmax; j++)
//					{
//						ofile_err << all_error_positions[i][j] << " ";
//					}
//					ofile_err << endl;
//				}
//				ofile_err << endl;
//				ofile_err.close();
//			}
//
//		}
//		//Decoding_Duration = Decoding_Duration * 100 / CLOCKS_PER_SEC / frame;
//		//Decoding_Duration = Decoding_Duration / CLOCKS_PER_SEC / Max_frame;
//
//		//printf("%f  %f  %f  %f\n", nSNR, Decoding_Duration, (double)Num_Error / N_Info / Max_frame, (double)Num_Frame_Error / Max_frame);
//		printf("%f %f %f\n", nSNR, (double)Num_Frame_Error / frame, (double)Num_Error / N_Info / frame);
//		BER_array[SNR_count] = (double)Num_Error / N_Info / frame;
//		FER_array[SNR_count] = (double)Num_Frame_Error / frame;
//		// save to .txt file
//		string polarde_serial = "";
//		for (int i = 0; i < softiter; i++)
//			polarde_serial = polarde_serial + polarde_array[i] + string("_");
//		string file_name = string("ResFinal\\Irregular_") + string("Border") + to_string(border) + string("N") + to_string(Nh) + string("K0") + to_string(K_time_array[0]) + string("K1") + to_string(K_time_array[1])
//			+ string("T") + to_string(softiter) + string("H") + to_string(half_iter) + string("UB") + to_string(usebeta) + string("IS") + to_string(flag_interleave_spatial)
//			+ string("DG") + to_string(diff_gain) + string(".txt");
//		save_array_txt(nSNR, FER_array[SNR_count], file_name);
//		int L_sdis_sum = L_LSD;
//		sdis_sum = 0;
//		for (int i = 0; i < Nv; i++) sdis_sum += sum_sdis_array[i];
//		if (polarde_array[0] == "SCL") L_sdis_sum = L;
//		// cout << "sdis_sum_avg" << sdis_sum / frame / L_sdis_sum << endl;
//		SNR_count++;
//		delete[] non_uniform_intrlv_pattern_array;
//	}
//	//fclose(stdout);
//	delete[] a_data;
//	for (int i = 0; i < Nv; i++)
//		delete[] u[i];
//	delete[] u;
//	for (int i = 0; i < Nv; i++)
//		delete[] umid[i];
//	delete[] umid;
//	for (int i = 0; i < Nv; i++)
//		delete[] uout[i];
//	delete[] uout;
//	for (int i = 0; i < Nv; i++)
//		delete[] x[i];
//	delete[] x;
//	for (int i = 0; i < Nv; i++)
//		delete[] x_mmse_out[i];
//	delete[] x_mmse_out;
//	for (int i = 0; i < Nv; i++)
//		delete[] midx[i];
//	delete[] midx;
//	for (int i = 0; i < Nv; i++)
//		delete[] y[i];
//	for (int i = 0; i < Nh; i++)
//		delete[] oriLLR_kbest[i];
//	delete[] oriLLR_kbest;
//	delete[] y;
//	delete[] u_crc;
//	delete[] Ah;
//	delete[] Ach;
//	delete[] A_Ach;
//	delete[] Av;
//	delete[] Acv;
//	delete[] A_Acv;
//	for (int i = 0; i < Nv; i++)
//		delete[] oriLLRh[i];
//	delete[] oriLLRh;
//	for (int i = 0; i < Nv; i++)
//		delete[] upLLR[i];
//	delete[] upLLR;
//	for (int i = 0; i < Nh; i++)
//		delete[] GGh[i];
//	delete[] GGh;
//	for (int i = 0; i < Nv; i++)
//		delete[] GGv[i];
//	delete[] GGv;
//	for (int i = 0; i < mimo_blk_height; i++) {
//		delete[] x_mod[i];
//	}
//	delete[] x_mod;
//	for (int i = 0; i < Nr; i++) {
//		delete[] y_mod[i];
//	}
//	delete[] y_mod;
//	for (int i = 0; i < 2 * mimo_blk_height; i++) {
//		delete[] y_mmse[i];
//	}
//	delete[] y_mmse;
//	for (int i = 0; i < 2 * Nr; i++) {
//		delete[] y_mod_real[i];
//	}
//	for (int i = 0; i < mimo_blk_length; i++)
//	{
//		delete[] y_mod_real_temp_omp[i];
//	}
//	delete[] y_mod_real_temp_omp;
//	for (int i = 0; i < mimo_blk_length; i++)
//	{
//		delete[] y_mmse_temp_omp[i];
//	}
//	delete[] y_mmse_temp_omp;
//	for (int i = 0; i < mimo_blk_length; i++)
//	{
//		delete[] ep_extr_mean_omp[i];
//	}
//	delete[] ep_extr_mean_omp;
//	for (int i = 0; i < mimo_blk_length; i++)
//	{
//		delete[] ep_extr_var_omp[i];
//	}
//	delete[] ep_extr_var_omp;
//	delete[] y_mod_real;
//	delete[] x_mod_tmp;
//	delete[] y_mod_tmp;
//	delete[] y_mod_real_tmp;
//	delete[] x_tmp;
//	delete[] oriLLR_tmp;
//	delete[] sum_sdis_array;
//	delete[] Hreal_array;
//	delete[] HTH_array;
//	delete[] HTY_array;
//	delete[] intrlv_pattern_2D_inner;
//}


