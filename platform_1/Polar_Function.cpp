#include <xmmintrin.h>
#include <immintrin.h>
#include <smmintrin.h>
#include "Polar_Decoder.h"
#include <iostream>
#include <vector>
#include "setting.h"
#include <iomanip>
#include "Polar_Encoder.h"
using namespace std;
int sign(float x)
{
	if (x > 0) return 1;
	if (x < 0) return -1;
	return 0;
}

/*
Function: void PolarDecodeIterVertical(float** oriLLRh, float** upLLR, int Nh, int Nv, int decodingi, int L_SCL)
oriLLRh: the channel LLR, size = Nv x Nh;
upLLR: the previously updated LLR, also the LLR that is going to be updated by the function, size = Nv x Nh
umid: the intermidiate decoder output. It is modified only when
A_Acv: the info bits and frozen bits of the vertical polar code
decodingi: specifies which column is currently being decoded
*/

void MITSOv(float** LLR, float* PM, int** sum, float** oriLLRh, int* indexpm, int iter, int decodingi, float** upLLR, float Qsumun)
{
	for (int j = 0; j < Nv; j++)
	{
		// LLR to Prob
		float P0 = exp(LLR[0][j]) / (1 + exp(LLR[0][j]));
		float P1 = 1 / (1 + exp(LLR[0][j]));
		float num = 0, den = 0;
		for (int i = 0; i < L; i++)
		{
			//cout << PM[i] << endl;
			if (!sum[i][j]) num += exp(PM[i]);
			else den += exp(PM[i]);
		}
		//cout << num << " " << den << endl;
		num += Qsumun * P0;
		den += Qsumun * P1;
		//cout << num << " " << den << " "<< Qsumun<<endl;
		float clog = 0;
		if (den == 0) clog = CLOG_RESTRICTION;
		else
		{
			clog = log(num / den);
			if (clog > CLOG_RESTRICTION) clog = CLOG_RESTRICTION;
			if (clog < -CLOG_RESTRICTION) clog = -CLOG_RESTRICTION;
		}
		// Extrinsic
		upLLR[j][decodingi] = oriLLRh[j][decodingi] + SOalpha * (clog - LLR[0][j]);
	}
}
void MITSOh(float** LLR, float* PM, int** sum, float** oriLLRh, int* indexpm, int iter, int decodingi, float** upLLR, float Qsumun)
{
	for (int j = 0; j < Nh; j++)
	{
		float P0 = exp(LLR[0][j]) / (1 + exp(LLR[0][j]));
		float P1 = 1 / (1 + exp(LLR[0][j]));
		float num = 0, den = 0;
		for (int i = 0; i < L; i++)
		{
			//cout << PM[i] << endl;
			if (!sum[i][j]) num += exp(PM[i]);
			else den += exp(PM[i]);
		}
		num += Qsumun * P0;
		den += Qsumun * P1;
		float clog = 0;
		if (den == 0) clog = CLOG_RESTRICTION;
		else
		{
			clog = log(num / den);
			if (clog > CLOG_RESTRICTION) clog = CLOG_RESTRICTION;
			if (clog < -CLOG_RESTRICTION) clog = -CLOG_RESTRICTION;
		}
		upLLR[decodingi][j] = oriLLRh[decodingi][j] + SOalpha * (clog - LLR[0][j]);
	}
}
void MITSOh_irregular(float** LLR, float* PM, int** sum, float** oriLLRh, int* indexpm, int iter, int decodingi, float** upLLR, float Qsumun, int* irregular_construction)
{
	for (int j = 0; j < Nh; j++)
	{
		float P0 = exp(LLR[0][j]) / (1 + exp(LLR[0][j]));
		float P1 = 1 / (1 + exp(LLR[0][j]));
		float num = 0, den = 0;
		for (int i = 0; i < L; i++)
		{
			//cout << PM[i] << endl;
			if (!sum[i][j]) num += exp(PM[i]);
			else den += exp(PM[i]);
		}
		num += Qsumun * P0;
		den += Qsumun * P1;
		float clog = 0;
		// float clog = 0;
		if (den == 0) clog = CLOG_RESTRICTION;
		else
		{
			clog = log(num / den);
			if (clog > CLOG_RESTRICTION) clog = CLOG_RESTRICTION;
			if (clog < -CLOG_RESTRICTION) clog = -CLOG_RESTRICTION;
		}
		float r = 1.0;
		if (irregular_construction[decodingi] == 1) { r = amplify_ratio_array[iter-1]; }
		upLLR[decodingi][j] = oriLLRh[decodingi][j] + SOalpha * (r*clog - LLR[0][j]);
	}
}
void Pyndiahv(float** LLR, int** sum, float** oriLLRh, int* indexpm, int iter, int decodingi, float** upLLR)
{
	float sdis[L];
	for (int i = 0; i < L; i++)
	{
		sdis[i] = 0;
		for (int j = 0; j < Nv; j++)
			sdis[i] = sdis[i] + (1 - sign(LLR[0][j]) * (1 - 2 * sum[i][j])) * abs(LLR[0][j]) / 2;
	}
	//TEST:
	/*cout << "Pyndiah Dist" << endl;
	cout << endl;
	for (int i = 0; i < L; i++)
	{
		cout << sdis[i] << "   ";
	}
	cout << endl;
	cout << "Pyndiah" << endl;*/
	//TEST
	// Test Magnitude Restrict
	for (int i = 0; i < L; i++)
	{
		if (sdis[i] > 10000) sdis[i] = 10000;
		if (sdis[i] < -10000) sdis[i] = -10000;
	}
	if (showllr)
	{
		cout << "Pyndiah_SCL_v" << "iter=  " << iter << endl;
		for (int i = 0; i < L; i++) cout << setw(10) << sdis[i] << endl;
		cout << endl;
	}
	if (showllr) {
		cout << "original_llr" << endl;
		for (int j = 0; j < Nv; j++)
		{
			cout << setw(10) << upLLR[j][decodingi];
		}
		cout << endl;
		// test:upLLR
	}
	// Test Magnitude Restrict 
	// Test Original LLR
	// Test Original LLR
	float w[Nv];
	if (usebeta && restrict_llr)
	{
		for (int j = 0; j < Nv; j++)
		{
			if (upLLR[j][decodingi] < -beta[iter * 2 - 2])
				upLLR[j][decodingi] = -beta[iter * 2 - 2];
			if (upLLR[j][decodingi] > beta[iter * 2 - 2])
				upLLR[j][decodingi] = beta[iter * 2 - 2];
		}
	}
	for (int j = 0; j < Nv; j++)
	{
		int firstdiff = 0;
		for (int i = 1; i < L; i++)
			if (sum[indexpm[i]][j] != sum[indexpm[0]][j])
			{
				firstdiff = i;
				break;
			}
		if (firstdiff)
		{
			w[j] = (sdis[indexpm[firstdiff]] - sdis[indexpm[0]]) * (1 - 2 * sum[indexpm[0]][j]);
			if (usebeta && restrict_llr)
			{
				if (w[j] < -beta[iter * 2 - 2])
					w[j] = -beta[iter * 2 - 2];
				else if (w[j] > beta[iter * 2 - 2])
					w[j] = beta[iter * 2 - 2];
			}
		}
		else
		{
			if (usebeta) w[j] = beta[iter * 2 - 2] * (1 - 2 * sum[indexpm[0]][j]);
			else w[j] = (sdis[indexpm[L - 1]] - sdis[indexpm[0]]) * (1 - 2 * sum[indexpm[0]][j]);
		}
		/*upLLR[j][decodingi] = oriLLRh[j][decodingi] + alpha[iter * 2 - 2] * (w[j] - upLLR[j][decodingi] * (1.0));*/
		upLLR[j][decodingi] = oriLLRh[j][decodingi] + alpha[iter * 2 - 2] * (w[j] - upLLR[j][decodingi] * (1.0));
		// cout << upLLR[j][decodingi] << " ";
		if (usebeta && restrict_llr)
		{
			for (int j = 0; j < Nv; j++)
			{
				if (upLLR[j][decodingi] < -beta[iter * 2 - 2])
					upLLR[j][decodingi] = -beta[iter * 2 - 2];
				if (upLLR[j][decodingi] > beta[iter * 2 - 2])
					upLLR[j][decodingi] = beta[iter * 2 - 2];
			}
		}
	}
	// test:upLLR
	if (showllr) {
		cout << "upLLR_renew" << endl;
		for (int j = 0; j < Nv; j++)
		{
			cout << setw(10) << upLLR[j][decodingi];
		}
		cout << endl;
		cout << "w" << endl;
		for (int j = 0; j < Nv; j++)
		{
			cout << setw(10) << w[j];
		}
		cout << endl;
		// test:upLLR
	}
}

void Pyndiahh(float** LLR, int** sum, float** oriLLRh, int* indexpm, int iter, int decodingi, float** upLLR)
{
	float sdis[L];
	for (int i = 0; i < L; i++)
	{
		sdis[i] = 0;
		for (int j = 0; j < Nh; j++)
			sdis[i] = sdis[i] + (1 - sign(LLR[0][j]) * (1 - 2 * sum[i][j])) * abs(LLR[0][j]) / 2;
	}

	////TEST:
	//cout << "Pyndiah Dist" << endl;
	//cout << endl;
	//for (int i = 0; i < L_LSD; i++)
	//{
	//	cout << sdis[i] << "   ";
	//}
	//cout << endl;
	//cout << "Pyndiah" << endl;
	////TEST

	float w[Nh];
	for (int j = 0; j < Nh; j++)
	{
		int firstdiff = 0;
		for (int i = 1; i < L; i++)
			if (sum[indexpm[i]][j] != sum[indexpm[0]][j])
			{
				firstdiff = i;
				break;
			}
		if (firstdiff)
			w[j] = (sdis[indexpm[firstdiff]] - sdis[indexpm[0]]) * (1 - 2 * sum[indexpm[0]][j]);
		else
			if (usebeta) w[j] = beta[iter * 2 - 1] * (1 - 2 * sum[indexpm[0]][j]);
			else w[j] = (sdis[indexpm[L - 1]] - sdis[indexpm[0]]) * (1 - 2 * sum[indexpm[0]][j]);
		upLLR[decodingi][j] = oriLLRh[decodingi][j] + alpha[iter * 2 - 1] * (w[j] - upLLR[decodingi][j] * (1.0 - usebeta));
		//upLLR[j][decodingi] = oriLLRh[j][decodingi] + alpha[iter * 2 - 1] * (w[j]);
		//	cout << upLLR[decodingi][j] << " ";
	}
}

void Pyndiahh_getsum(float** LLR, int** sum, float** oriLLRh, int* indexpm, int iter, int decodingi, float** upLLR, double& sum_sdis)
{
	float sdis[L];
	for (int i = 0; i < L; i++)
	{
		sdis[i] = 0;
		for (int j = 0; j < Nh; j++)
			sdis[i] = sdis[i] + (1 - sign(LLR[0][j]) * (1 - 2 * sum[i][j])) * abs(LLR[0][j]) / 2;
	}

	for (int i = 0; i < L; i++) { sum_sdis += sdis[i]; }
	////TEST:
	//cout << "Pyndiah Dist" << endl;
	//cout << endl;
	//for (int i = 0; i < L_LSD; i++)
	//{
	//	cout << sdis[i] << "   ";
	//}
	//cout << endl;
	//cout << "Pyndiah" << endl;
	////TEST

	float w[Nh];
	for (int j = 0; j < Nh; j++)
	{
		int firstdiff = 0;
		for (int i = 1; i < L; i++)
			if (sum[indexpm[i]][j] != sum[indexpm[0]][j])
			{
				firstdiff = i;
				break;
			}
		if (firstdiff)
			w[j] = (sdis[indexpm[firstdiff]] - sdis[indexpm[0]]) * (1 - 2 * sum[indexpm[0]][j]);
		else
			if (usebeta) w[j] = beta[iter * 2 - 1] * (1 - 2 * sum[indexpm[0]][j]);
			else w[j] = (sdis[indexpm[L - 1]] - sdis[indexpm[0]]) * (1 - 2 * sum[indexpm[0]][j]);
		upLLR[decodingi][j] = oriLLRh[decodingi][j] + alpha[iter * 2 - 1] * (w[j] - upLLR[decodingi][j] * (1.0 - usebeta));
		//upLLR[j][decodingi] = oriLLRh[j][decodingi] + alpha[iter * 2 - 1] * (w[j]);
		//	cout << upLLR[decodingi][j] << " ";
	}
}

void Pyndiahh_LSD(float** LLR, int** sum, float** oriLLRh, int* indexpm, int iter, int decodingi, float** upLLR)
{
	float sdis[L_LSD];
	for (int i = 0; i < L_LSD; i++)
	{
		sdis[i] = 0;
		for (int j = 0; j < Nh; j++)
			sdis[i] = sdis[i] + (1 - sign(LLR[0][j]) * (1 - 2 * sum[i][j])) * abs(LLR[0][j]) / 2;
	}

	/*for (int i = 0; i < L_LSD; i++)
	{
		sdis_sum += sdis[i];
	}
	cout << sdis_sum << endl;*/
	//TEST:
	/*cout << "Pyndiah Dist" << endl;
	cout << endl;
	for (int i = 0; i < L_LSD; i++)
	{
		cout << sdis[i] << "   ";
	}
	cout << endl;
	cout << "Pyndiah" << endl;*/
	//TEST
	// Test Magnitude Restrict
	for (int i = 0; i < L_LSD; i++)
	{
		if (sdis[i] > 10000) sdis[i] = 10000;
		if (sdis[i] < -10000) sdis[i] = -10000;
	}
	/*if (showllr)
	{
		cout << "Pyndiah_LSD" << "iter=  " << iter << endl;
		for (int i = 0; i < L_LSD; i++) cout << sdis[i] << "   ";
		cout << endl;
	}*/
	// Test Magnitude Restrict 
	float w[Nh];
	if (usebeta && restrict_llr)
	{
		for (int j = 0; j < Nv; j++)
		{
			if (upLLR[j][decodingi] < -beta[iter * 2 - 2])
				upLLR[j][decodingi] = -beta[iter * 2 - 2];
			if (upLLR[j][decodingi] > beta[iter * 2 - 2])
				upLLR[j][decodingi] = beta[iter * 2 - 2];
		}
	}
	for (int j = 0; j < Nh; j++)
	{
		int firstdiff = 0;
		for (int i = 1; i < L_LSD; i++)
			if (sum[indexpm[i]][j] != sum[indexpm[0]][j])
			{
				firstdiff = i;
				break;
			}
		if (firstdiff)
		{
			w[j] = (sdis[indexpm[firstdiff]] - sdis[indexpm[0]]) * (1 - 2 * sum[indexpm[0]][j]);
			if (usebeta && restrict_llr)
			{
				if (w[j] < -beta[iter * 2 - 1])
					w[j] = -beta[iter * 2 - 1];
				else if (w[j] > beta[iter * 2 - 1])
					w[j] = beta[iter * 2 - 1];
			}
		}
		else
			if (usebeta) w[j] = beta[iter * 2 - 1] * (1 - 2 * sum[indexpm[0]][j]);
		//else w[j] = (sdis[indexpm[L_LSD - 1]] - sdis[indexpm[0]]) * (1 - 2 * sum[indexpm[0]][j]);
			else w[j] = (sdis[indexpm[L_LSD - 1]] - sdis[indexpm[0]]) * (1 - 2 * sum[indexpm[0]][j]);
		upLLR[decodingi][j] = oriLLRh[decodingi][j] + alpha[iter * 2 - 1] * (w[j] - upLLR[decodingi][j] * (1.0));
		//upLLR[j][decodingi] = oriLLRh[j][decodingi] + alpha[iter * 2 - 1] * (w[j]);

		if (usebeta && restrict_llr)
		{
			for (int j = 0; j < Nv; j++)
			{
				if (upLLR[j][decodingi] < -beta[iter * 2 - 1])
					upLLR[j][decodingi] = -beta[iter * 2 - 1];
				if (upLLR[j][decodingi] > beta[iter * 2 - 1])
					upLLR[j][decodingi] = beta[iter * 2 - 1];
			}
		}
		//upLLR[decodingi][j] = oriLLRh[decodingi][j] + 0.3 * (w[j] - upLLR[decodingi][j]);
		//	cout << upLLR[decodingi][j] << " ";
		// upLLR[decodingi][j] = oriLLRh[decodingi][j] + alpha[iter * 2 - 1] * (w[j] - upLLR[decodingi][j]);
	}
}

void my_sort_pyndiah(int* arr, int* idx_sort, int N) {
	for (int i = 0; i < N; ++i) {
		// cout << N << endl;
		idx_sort[i] = i;
	}
	std::sort(idx_sort, idx_sort + N, [&](int a, int b) {
		return arr[a] > arr[b];
		});

	float* sortedArr = new float[N];
	for (int i = 0; i < N; ++i) {
		sortedArr[i] = arr[idx_sort[i]];
	}


	for (int i = 0; i < N; ++i) {
		arr[i] = sortedArr[i];
	}

	delete[] sortedArr;
}


void my_sort_pyndiah(float* arr, int* idx_sort, int N) {
	for (int i = 0; i < N; ++i) {
		// cout << N << endl;
		idx_sort[i] = i;
	}
	std::sort(idx_sort, idx_sort + N, [&](int a, int b) {
		return arr[a] > arr[b];
		});

	float* sortedArr = new float[N];
	for (int i = 0; i < N; ++i) {
		sortedArr[i] = arr[idx_sort[i]];
	}


	for (int i = 0; i < N; ++i) {
		arr[i] = sortedArr[i];
	}

	delete[] sortedArr;
}

void Pyndiahh_LSD_getsum(float** LLR, int** sum, float** oriLLRh, int* indexpm, int iter, int decodingi, float** upLLR, double& sdis_sum, int** GGh, int* Ach)
{
	int p = 15;
	int n_samples = (int)pow(2, p);
	int num_first_diff = 0;
	int* first_diff_array = new int[Nh];
	int* idx_first_diff_sort = new int[Nh];
	float* sdis_sampled = new float[n_samples];
	int* valid_sample = new int[n_samples];
	int* bit_tmp = new int[Nh];
	int* bit_tmp_encoded = new int[Nh];
	for (int i = 0; i < n_samples; i++) { sdis_sampled[i] = 0; valid_sample[i] = 1; }
	for (int i = 0; i < Nh; i++) { idx_first_diff_sort[i] = -10 * Nh; }
	for (int i = 0; i < Nh; i++) { first_diff_array[i] = 0; bit_tmp[i] = 0; bit_tmp_encoded[i] = 0; }
	// decide bits with different decisions 
	for (int j = 0; j < Nh; j++)
	{
		int firstdiff = 0;
		for (int i = 1; i < L_LSD; i++)
			if (sum[indexpm[i]][j] != sum[indexpm[0]][j])
			{
				first_diff_array[j] = i;
				num_first_diff++;
				break;
			}
	}
	// arrange the positions of first error in ascending order to find the p bits with least certain decisions
	my_sort_pyndiah(first_diff_array, idx_first_diff_sort, Nh);
	float sdis[L_LSD];
	for (int i = 0; i < L_LSD; i++)
	{
		sdis[i] = 0;
		for (int j = 0; j < Nh; j++)
			sdis[i] = sdis[i] + (1 - sign(LLR[0][j]) * (1 - 2 * sum[i][j])) * abs(LLR[0][j]) / 2;
	}
	for (int i = 0; i < L_LSD; i++)
	{
		sdis_sum += sdis[i];
	}
	int num_firstdiff = 0;


	// test: LSD+Random Sampling
	// calculate the Euclidian distance of each sample.
	//int i_code = 0;
	//for (int i = 0; i < n_samples; i++)
	//{
	//	for (int ibit = 0; ibit < Nh; ibit++) {
	//		bit_tmp[ibit] = sum[indexpm[i_code]][ibit];
	//	}
	//	// bitset: 0 (LSB) to 19 (MSB)
	//	sdis_sampled[i] = sdis[indexpm[i_code]];
	//	std::bitset<20> candidate_bits(i+1);
	//	// cout << candidate_bits[0];
	//	for (int j = 0; j < p; j++)
	//	{
	//		int bit_flip_position = idx_first_diff_sort[j];
	//		bool bflip = (int)candidate_bits[j];
	//		if (!bflip) { continue; }
	//		else
	//		{
	//			// determined which bit is being flipped
	//			int flipped_bit = (sum[indexpm[i_code]][bit_flip_position]<0.5);
	//			bit_tmp[bit_flip_position] = flipped_bit;
	//			sdis_sampled[i] -= (1 - sign(LLR[0][bit_flip_position]) * (1 - 2 * sum[indexpm[i_code]][bit_flip_position])) * abs(LLR[0][bit_flip_position]) / 2;
	//			sdis_sampled[i] += (1 - sign(LLR[0][bit_flip_position]) * (1 - 2 * flipped_bit)) * abs(LLR[0][bit_flip_position]) / 2;
	//		}
	//	}
	//	PolarEncode_xor(bit_tmp_encoded, bit_tmp, Nh);
	//	for (int j = 0; j < Nh-Kh; j++)
	//	{
	//		if (bit_tmp_encoded[Ach[j]] != 0) { valid_sample[i] = 0; }
	//	}
	//}


	// cout << sdis_sum << endl;
	//TEST:
	/*cout << "Pyndiah Dist" << endl;
	cout << endl;
	for (int i = 0; i < L_LSD; i++)
	{
		cout << sdis[i] << "   ";
	}
	cout << endl;
	cout << "Pyndiah" << endl;

	cout << "Appended Dist" << endl;
	for (int i = 0; i < n_samples; i++)
	{
		if (valid_sample[i] && sdis_sampled[i]<sdis[L_LSD-1]) cout << sdis_sampled[i] << "     ";
	}
	cout << "Appended Dist" << endl;*/
	//TEST
	// Test Magnitude Restrict
	for (int i = 0; i < L_LSD; i++)
	{
		if (sdis[i] > 10000) sdis[i] = 10000;
		if (sdis[i] < -10000) sdis[i] = -10000;
	}
	/*if (showllr)
	{
		cout << "Pyndiah_LSD" << "iter=  " << iter << endl;
		for (int i = 0; i < L_LSD; i++) cout << sdis[i] << "   ";
		cout << endl;
	}*/
	// Test Magnitude Restrict 
	float w[Nh];
	if (usebeta)
	{
		for (int j = 0; j < Nv; j++)
		{
			if (upLLR[j][decodingi] < -beta[iter * 2 - 2])
				upLLR[j][decodingi] = -beta[iter * 2 - 2];
			if (upLLR[j][decodingi] > beta[iter * 2 - 2])
				upLLR[j][decodingi] = beta[iter * 2 - 2];
		}
	}
	for (int j = 0; j < Nh; j++)
	{
		int firstdiff = 0;
		for (int i = 1; i < L_LSD; i++)
			if (sum[indexpm[i]][j] != sum[indexpm[0]][j])
			{
				firstdiff = i;
				break;
			}
		if (firstdiff)
	    //if (false)
		{
			num_firstdiff++;
			w[j] = (sdis[indexpm[firstdiff]] - sdis[indexpm[0]]) * (1 - 2 * sum[indexpm[0]][j]);
			if (usebeta)
			{
				if (w[j] < -beta[iter * 2 - 1])
					w[j] = -beta[iter * 2 - 1];
				else if (w[j] > beta[iter * 2 - 1])
					w[j] = beta[iter * 2 - 1];
			}
		}
		else
			if (usebeta) w[j] = beta[iter * 2 - 1] * (1 - 2 * sum[indexpm[0]][j]);
		//else w[j] = (sdis[indexpm[L_LSD - 1]] - sdis[indexpm[0]]) * (1 - 2 * sum[indexpm[0]][j]);
			else w[j] = (sdis[indexpm[L_LSD - 1]] - sdis[indexpm[0]]) * (1 - 2 * sum[indexpm[0]][j]);
		upLLR[decodingi][j] = oriLLRh[decodingi][j] + alpha[iter * 2 - 1] * (w[j] - upLLR[decodingi][j] * (1.0));
		// test w, upLLR, oriLLR
		
		// test w, upLLR, oriLLR
		//upLLR[j][decodingi] = oriLLRh[j][decodingi] + alpha[iter * 2 - 1] * (w[j]);
		if (usebeta)
		{
			for (int j = 0; j < Nv; j++)
			{
				if (upLLR[j][decodingi] < -beta[iter * 2 - 1])
					upLLR[j][decodingi] = -beta[iter * 2 - 1];
				if (upLLR[j][decodingi] > beta[iter * 2 - 1])
					upLLR[j][decodingi] = beta[iter * 2 - 1];
			}
		}
		//upLLR[decodingi][j] = oriLLRh[decodingi][j] + 0.3 * (w[j] - upLLR[decodingi][j]);
		//	cout << upLLR[decodingi][j] << " ";
		// upLLR[decodingi][j] = oriLLRh[decodingi][j] + alpha[iter * 2 - 1] * (w[j] - upLLR[decodingi][j]);
	}
	// cout << "NUM_FIRSTDIFF:  " << num_firstdiff << endl;
	delete[] sdis_sampled;
	delete[] idx_first_diff_sort;
	delete[] first_diff_array;
	delete[] valid_sample;
	delete[] bit_tmp;
	delete[] bit_tmp_encoded;
}

void Pyndiahh_LSD_partial(float** LLR, int** sum, float** oriLLRh, int* indexpm, int iter, int decodingi, float** upLLR)
{
	float sdis[L_LSD];
	for (int i = 0; i < L_LSD; i++)
	{
		sdis[i] = 0;
		for (int j = 0; j < Nh; j++)
			sdis[i] = sdis[i] + (1 - sign(LLR[0][j]) * (1 - 2 * sum[i][j])) * abs(LLR[0][j]) / 2;
	}
	float w[Nh];
	for (int j = 0; j < Nh; j++)
	{
		int firstdiff = 0;
		for (int i = 1; i < L_LSD; i++)
			if (sum[indexpm[i]][j] != sum[indexpm[0]][j])
			{
				firstdiff = i;
				break;
			}
		if (firstdiff)
			w[j] = (sdis[indexpm[firstdiff]] - sdis[indexpm[0]]) * (1 - 2 * sum[indexpm[0]][j]);
		else
			if (usebeta) w[j] = beta[iter * 2 - 1] * (1 - 2 * sum[indexpm[0]][j]);
			else w[j] = (sdis[indexpm[L_LSD - 1]] - sdis[indexpm[0]]) * (1 - 2 * sum[indexpm[0]][j]);
		upLLR[decodingi][j] = oriLLRh[decodingi][j] + alpha[iter * 2 - 1] * (w[j] - upLLR[decodingi][j]);
		//	cout << upLLR[decodingi][j] << " ";
	}
}

int judge_rate(int* A_Ac, int len)
{
	bool temp0 = false;
	for (int k = len - 2; k >= 0; k--)
	{
		temp0 |= *(A_Ac + k);
		if (temp0) { break; return 4; }
	}
	if (temp0 == false)
	{
		if (temp0 || *(A_Ac + len - 1)) return 2; //repetition node
		else return 0;  //rate-0 node
	}
	else {
		bool temp1 = true;
		for (int k = 1; k < len; k++)
		{
			temp1 &= *(A_Ac + k);
			if (!temp1) { break; return 4; }
		}
		if (temp1 == true)
		{
			if (temp1 && *A_Ac) return 1; //rate-1 node
			else return 3; // SPC node
		}
		else return 4;
	}
}



// Use Adaptive Beta
void Pyndiahv_adaptivebeta(float** LLR, int** sum, float** oriLLRh, int* indexpm, int iter, int decodingi, float** upLLR)
{
	float sdis[L];
	for (int i = 0; i < L; i++)
	{
		sdis[i] = 0;
		for (int j = 0; j < Nv; j++)
			sdis[i] = sdis[i] + (1 - sign(LLR[0][j]) * (1 - 2 * sum[i][j])) * abs(LLR[0][j]) / 2;
	}
	float w[Nv];
	// adaptive beta
	for (int j = 0; j < Nv; j++)
	{
		int firstdiff = 0;
		for (int i = 1; i < L; i++)
			if (sum[indexpm[i]][j] != sum[indexpm[0]][j])
			{
				firstdiff = i;
				break;
			}
		if (firstdiff)
			w[j] = (sdis[indexpm[firstdiff]] - sdis[indexpm[0]]) * (1 - 2 * sum[indexpm[0]][j]);
		else
			if (usebeta)
			{
				if (adaptive_beta)
				{
					float beta_adapt = a_beta[iter * 2 - 2] - sdis[indexpm[0]] * k_beta[iter * 2 - 2];
					if (beta_adapt > b_beta[iter * 2 - 2]) {
						beta_adapt = b_beta[iter * 2 - 2];
					}
					w[j] = beta_adapt * (1 - 2 * sum[indexpm[0]][j]);
				}
				else
				{
					w[j] = beta[iter * 2 - 2] * (1 - 2 * sum[indexpm[0]][j]);
				}
			}
			else w[j] = (sdis[indexpm[L - 1]] - sdis[indexpm[0]]) * (1 - 2 * sum[indexpm[0]][j]);

		upLLR[j][decodingi] = oriLLRh[j][decodingi] + alpha[iter * 2 - 2] * (w[j] - upLLR[j][decodingi]);
		//cout << upLLR[j][decodingi] << " ";
	}
}

void Pyndiahh_adaptivebeta(float** LLR, int** sum, float** oriLLRh, int* indexpm, int iter, int decodingi, float** upLLR)
{
	float sdis[L];
	for (int i = 0; i < L; i++)
	{
		sdis[i] = 0;
		for (int j = 0; j < Nh; j++)
			sdis[i] = sdis[i] + (1 - sign(LLR[0][j]) * (1 - 2 * sum[i][j])) * abs(LLR[0][j]) / 2;
	}
	float w[Nh];
	for (int j = 0; j < Nh; j++)
	{
		int firstdiff = 0;
		for (int i = 1; i < L; i++)
			if (sum[indexpm[i]][j] != sum[indexpm[0]][j])
			{
				firstdiff = i;
				break;
			}
		if (firstdiff)
			w[j] = (sdis[indexpm[firstdiff]] - sdis[indexpm[0]]) * (1 - 2 * sum[indexpm[0]][j]);
		else
			if (usebeta)
			{
				if (adaptive_beta)
				{
					float beta_adapt = a_beta[iter * 2 - 1] - sdis[indexpm[0]] * k_beta[iter * 2 - 1];
					if (beta_adapt > b_beta[iter * 2 - 1]) {
						beta_adapt = b_beta[iter * 2 - 1];
					}
					w[j] = beta_adapt * (1 - 2 * sum[indexpm[0]][j]);
				}
				else
				{
					w[j] = beta[iter * 2 - 1] * (1 - 2 * sum[indexpm[0]][j]);
				}
			}
			else w[j] = (sdis[indexpm[L - 1]] - sdis[indexpm[0]]) * (1 - 2 * sum[indexpm[0]][j]);  //³ËÏµÊý
		upLLR[decodingi][j] = oriLLRh[decodingi][j] + alpha[iter * 2 - 1] * (w[j] - upLLR[decodingi][j]);
		//	cout << upLLR[decodingi][j] << " ";
	}
}

void Pyndiahh_LSD_adaptivebeta(float** LLR, int** sum, float** oriLLRh, int* indexpm, int iter, int decodingi, float** upLLR)
{
	float sdis[L_LSD];
	for (int i = 0; i < L_LSD; i++)
	{
		sdis[i] = 0;
		for (int j = 0; j < Nh; j++)
			sdis[i] = sdis[i] + (1 - sign(LLR[0][j]) * (1 - 2 * sum[i][j])) * abs(LLR[0][j]) / 2;
	}
	//TEST:
	/*cout << endl;
	for (int i = 0; i < L_LSD; i++)
	{
		cout << sdis[i] << "   ";
	}
	cout << endl;*/
	//TEST
	float w[Nh];
	for (int j = 0; j < Nh; j++)
	{
		int firstdiff = 0;
		for (int i = 1; i < L_LSD; i++)
			if (sum[indexpm[i]][j] != sum[indexpm[0]][j])
			{
				firstdiff = i;
				break;
			}
		if (firstdiff)
			w[j] = (sdis[indexpm[firstdiff]] - sdis[indexpm[0]]) * (1 - 2 * sum[indexpm[0]][j]);
		else
			if (usebeta)
			{
				if (adaptive_beta)
				{
					float beta_adapt = a_beta[iter * 2 - 1] - sdis[indexpm[0]] * k_beta[iter * 2 - 1];
					if (beta_adapt > b_beta[iter * 2 - 1]) {
						beta_adapt = b_beta[iter * 2 - 1];
					}
					w[j] = beta_adapt * (1 - 2 * sum[indexpm[0]][j]);
				}
				else
				{
					w[j] = beta[iter * 2 - 1] * (1 - 2 * sum[indexpm[0]][j]);
				}
			}
			else w[j] = (sdis[indexpm[L_LSD - 1]] - sdis[indexpm[0]]) * (1 - 2 * sum[indexpm[0]][j]);
		upLLR[decodingi][j] = oriLLRh[decodingi][j] + alpha[iter * 2 - 1] * (w[j] - upLLR[decodingi][j]);
		//	cout << upLLR[decodingi][j] << " ";
	}
}




int judge_rate_L(int* A_Ac, int len)
{
	bool temp0 = false;
	int cnt = 0;
	for (int k = len - 1; k >= 0; k--)
	{
		temp0 |= *(A_Ac + k);
		cnt += *(A_Ac + k);
	}
	if (temp0 == false)
	{
		return 0;  //rate-0 node
	}
	else {
		bool temp1 = true;
		for (int k = 0; k < len; k++)
		{
			temp1 &= *(A_Ac + k);
			if (!temp1) { break; }
		}
		if (temp1 == true)
		{
			return 1; // SPC node
		}
		else {
			if ((cnt <= 1) && len <= 16) return 2;
			else return 4;
		}
	}
}
void f_function(float* a, int size)
{
	if (size <= 4) {
		for (int i = 0; i < size; i += 4)
		{
			__m128 a_left = _mm_loadu_ps(a + i);//
			__m128 a_right = _mm_loadu_ps(a + i + size);
			__m128 SIGN_MASK = _mm_castsi128_ps(_mm_set1_epi32(0x80000000));
			__m128 sign = _mm_and_ps(_mm_xor_ps(a_left, a_right), SIGN_MASK);
			__m128 abs_a_left = _mm_andnot_ps(SIGN_MASK, a_left);
			__m128 abs_a_right = _mm_andnot_ps(SIGN_MASK, a_right);
			__m128 s = _mm_or_ps(sign, _mm_min_ps(abs_a_left, abs_a_right));
			_mm_storeu_ps(a + i + (size << 1), s);
		}
	}
	else {
		for (int i = 0; i < size; i += 8)
		{
			__m256 a_left = _mm256_loadu_ps(a + i);
			__m256 a_right = _mm256_loadu_ps(a + i + size);
			__m256 SIGN_MASK = _mm256_castsi256_ps(_mm256_set1_epi32(0x80000000));
			__m256 sign = _mm256_and_ps(_mm256_xor_ps(a_left, a_right), SIGN_MASK);
			__m256 abs_a_left = _mm256_andnot_ps(SIGN_MASK, a_left);
			__m256 abs_a_right = _mm256_andnot_ps(SIGN_MASK, a_right);
			__m256 s = _mm256_or_ps(sign, _mm256_min_ps(abs_a_left, abs_a_right));
			_mm256_storeu_ps(a + i + (size << 1), s);
		}
	}
}

void f_function_index(float* a, float* a_new, int size)
{
	if (size <= 4) {
		for (int i = 0; i < size; i += 4)
		{
			__m128 a_left = _mm_loadu_ps(a + i);//
			__m128 a_right = _mm_loadu_ps(a + i + size);
			__m128 SIGN_MASK = _mm_castsi128_ps(_mm_set1_epi32(0x80000000));
			__m128 sign = _mm_and_ps(_mm_xor_ps(a_left, a_right), SIGN_MASK);
			__m128 abs_a_left = _mm_andnot_ps(SIGN_MASK, a_left);
			__m128 abs_a_right = _mm_andnot_ps(SIGN_MASK, a_right);
			__m128 s = _mm_or_ps(sign, _mm_min_ps(abs_a_left, abs_a_right));
			_mm_storeu_ps(a_new + i + (size << 1), s);
		}
	}
	else {
		for (int i = 0; i < size; i += 8)
		{
			__m256 a_left = _mm256_loadu_ps(a + i);
			__m256 a_right = _mm256_loadu_ps(a + i + size);
			__m256 SIGN_MASK = _mm256_castsi256_ps(_mm256_set1_epi32(0x80000000));
			__m256 sign = _mm256_and_ps(_mm256_xor_ps(a_left, a_right), SIGN_MASK);
			__m256 abs_a_left = _mm256_andnot_ps(SIGN_MASK, a_left);
			__m256 abs_a_right = _mm256_andnot_ps(SIGN_MASK, a_right);
			__m256 s = _mm256_or_ps(sign, _mm256_min_ps(abs_a_left, abs_a_right));
			_mm256_storeu_ps(a_new + i + (size << 1), s);
		}
	}
}

void g_function(float* a, int* b, int size)
{
	if (size <= 4) {
		for (int i = 0; i < size; i += 4)
		{
			__m128 a_left = _mm_loadu_ps(a + i);
			__m128 a_right = _mm_loadu_ps(a + i + size);
			__m128 b_in = _mm_loadu_ps((float*)(b + i));
			__m128 MASK_0 = _mm_castsi128_ps(_mm_set1_epi32(0x00000000));
			__m128 SIGN_MASK = _mm_castsi128_ps(_mm_set1_epi32(0x80000000));
			__m128 sign = _mm_and_ps(_mm_cmpneq_ps(b_in, MASK_0), SIGN_MASK);
			__m128 s = _mm_add_ps(a_right, _mm_xor_ps(sign, a_left));
			_mm_storeu_ps(a + i + (size << 1), s);
		}
	}
	else {
		for (int i = 0; i < size; i += 8)
		{
			__m256 a_left = _mm256_loadu_ps(a + i);
			__m256 a_right = _mm256_loadu_ps(a + i + size);
			__m256 b_in = _mm256_loadu_ps((float*)(b + i));
			__m256 MASK_0 = _mm256_castsi256_ps(_mm256_set1_epi32(0x00000000));
			__m256 SIGN_MASK = _mm256_castsi256_ps(_mm256_set1_epi32(0x80000000));
			__m256 sign = _mm256_and_ps(_mm256_cmp_ps(b_in, MASK_0, 0), SIGN_MASK);
			__m256 s = _mm256_sub_ps(a_right, _mm256_xor_ps(sign, a_left));
			_mm256_storeu_ps(a + i + (size << 1), s);
		}
	}
}

void g_function_index(float* a, int* b, float* a_new, int size)
{
	if (size <= 4) {
		for (int i = 0; i < size; i += 4)
		{
			__m128 a_left = _mm_loadu_ps(a + i);
			__m128 a_right = _mm_loadu_ps(a + i + size);
			__m128 b_in = _mm_loadu_ps((float*)(b + i));
			__m128 MASK_0 = _mm_castsi128_ps(_mm_set1_epi32(0x00000000));
			__m128 SIGN_MASK = _mm_castsi128_ps(_mm_set1_epi32(0x80000000));
			__m128 sign = _mm_and_ps(_mm_cmpneq_ps(b_in, MASK_0), SIGN_MASK);
			__m128 s = _mm_add_ps(a_right, _mm_xor_ps(sign, a_left));
			_mm_storeu_ps(a_new + i + (size << 1), s);
		}
	}
	else {
		for (int i = 0; i < size; i += 8)
		{
			__m256 a_left = _mm256_loadu_ps(a + i);
			__m256 a_right = _mm256_loadu_ps(a + i + size);
			__m256 b_in = _mm256_loadu_ps((float*)(b + i));
			__m256 MASK_0 = _mm256_castsi256_ps(_mm256_set1_epi32(0x00000000));
			__m256 SIGN_MASK = _mm256_castsi256_ps(_mm256_set1_epi32(0x80000000));
			__m256 sign = _mm256_and_ps(_mm256_cmp_ps(b_in, MASK_0, 0), SIGN_MASK);
			__m256 s = _mm256_sub_ps(a_right, _mm256_xor_ps(sign, a_left));
			_mm256_storeu_ps(a_new + i + (size << 1), s);
		}
	}
}

void combine(int* a, int* b, int size)
{
	if (size <= 4) {
		__m128i* ptr_a_128 = (__m128i*) a;
		__m128i* ptr_b_128 = (__m128i*) b;
		for (int i = 0; i < size / 4; i++)
		{
			__m128i a_128 = _mm_loadu_si128(ptr_a_128 + i);
			__m128i b_128 = _mm_loadu_si128(ptr_b_128 + i);
			__m128i s_128 = _mm_xor_si128(a_128, b_128);
			_mm_storeu_si128(ptr_a_128 + i, s_128);
		}
	}
	else {
		__m256i* ptr_a_256 = (__m256i*) a;
		__m256i* ptr_b_256 = (__m256i*) b;
		for (int i = 0; i < size / 8; i++)
		{
			__m256i a_256 = _mm256_loadu_si256(ptr_a_256 + i);
			__m256i b_256 = _mm256_loadu_si256(ptr_b_256 + i);
			__m256i s_256 = _mm256_xor_si256(a_256, b_256);
			_mm256_storeu_si256(ptr_a_256 + i, s_256);
		}
	}
}
void combine_index(int* a, int* b, int* a_new, int size)
{
	if (size <= 4) {
		__m128i* ptr_a_128 = (__m128i*) a;
		__m128i* ptr_b_128 = (__m128i*) b;
		__m128i* ptr_a_new_128 = (__m128i*) a_new;
		for (int i = 0; i < size / 4; i++)
		{
			__m128i a_128 = _mm_loadu_si128(ptr_a_128 + i);
			__m128i b_128 = _mm_loadu_si128(ptr_b_128 + i);
			_mm_storeu_si128(ptr_a_new_128 + i + (size / 4), b_128);
			__m128i s_128 = _mm_xor_si128(a_128, b_128);
			_mm_storeu_si128(ptr_a_new_128 + i, s_128);
		}
	}
	else {
		__m256i* ptr_a_256 = (__m256i*) a;
		__m256i* ptr_b_256 = (__m256i*) b;
		__m256i* ptr_a_new_256 = (__m256i*) a_new;
		for (int i = 0; i < size / 8; i++)
		{
			__m256i a_256 = _mm256_loadu_si256(ptr_a_256 + i);
			__m256i b_256 = _mm256_loadu_si256(ptr_b_256 + i);
			_mm256_storeu_si256(ptr_a_new_256 + i + (size / 8), b_256);
			__m256i s_256 = _mm256_xor_si256(a_256, b_256);
			_mm256_storeu_si256(ptr_a_new_256 + i, s_256);
		}
	}
}
void hard_SIMD(float* a, int* b, int size)
{
	switch (size) {
	case 2:
	{
		__m128i* ptr_b_128 = (__m128i*) b;
		__m128 llr_128 = _mm_loadu_ps(a);
		__m128 s = _mm_cmpge_ps(_mm_setzero_ps(), llr_128);
		__m128i mask = _mm_set1_epi32(0x00000001);
		__m128i sum_128 = _mm_and_si128(mask, _mm_castps_si128(s));
		*b = sum_128.m128i_i32[0];
		*(b + 1) = sum_128.m128i_i32[1];
	}
	case 4:
	{
		__m128i* ptr_b_128 = (__m128i*) b;
		for (int i = 0; i < size / 4; i++)
		{
			__m128 llr_128 = _mm_loadu_ps(a + (i << 2));
			__m128 s = _mm_cmpge_ps(_mm_setzero_ps(), llr_128);
			__m128i mask = _mm_set1_epi32(0x00000001);
			__m128i sum_128 = _mm_and_si128(mask, _mm_castps_si128(s));
			_mm_storeu_si128(ptr_b_128 + i, sum_128);
		}

	}
	default:
	{
		__m256i* ptr_b_256 = (__m256i*) b;
		for (int i = 0; i < size / 8; i++)
		{
			__m256 llr_256 = _mm256_loadu_ps(a + (i << 3));
			__m256 s = _mm256_cmp_ps(_mm256_setzero_ps(), llr_256, 14);
			__m256i mask = _mm256_set1_epi32(0x00000001);
			__m256i sum_256 = _mm256_and_si256(mask, _mm256_castps_si256(s));
			_mm256_storeu_si256(ptr_b_256 + i, sum_256);
		}
	}
	}
}
float cal_sum(float* a, int len)
{
	float temp = 0;
	if (len <= 4) {
		int cnt = (len >> 2) - 1;
		__m128 s = _mm_loadu_ps(a);
		__m128 llr_128 = _mm_loadu_ps(a);
		s = _mm_add_ps(s, llr_128);
		s = _mm_hadd_ps(s, s);
		s = _mm_hadd_ps(s, s);
		temp += s.m128_f32[0];
	}
	else
	{
		int cnt = len >> 3;
		__m256 s = _mm256_loadu_ps(a);
		for (int i = 1; i < cnt; i++)
		{
			a += 8;
			__m256 llr_256 = _mm256_loadu_ps(a);
			s = _mm256_add_ps(s, llr_256);
		}
		s = _mm256_hadd_ps(s, s);
		s = _mm256_hadd_ps(s, s);
		temp += s.m256_f32[0];
		temp += s.m256_f32[5];
	}
	return temp;
}

void replace_sum(int* a, int* b, int len) {
	if (len < 8) {
		int cntBlock = len >> 2;
		int cntRem = len & 3;
		__m128i* p_new = (__m128i*)a;
		__m128i* p_old = (__m128i*)b;
		__m128i xidLoad;
		for (int i = 0; i < cntRem; i++)
			*(b + len - i - 1) = *(a + len - i - 1);
		for (int i = 0; i < cntBlock; i++)
		{
			xidLoad = _mm_loadu_si128(p_new);
			_mm_storeu_si128(p_old, xidLoad);
			p_old++; p_new++;
		}

	}
	else {
		int cntBlock = len >> 3;
		int cntRem = len & 7;
		__m256i* p_new = (__m256i*)a;
		__m256i* p_old = (__m256i*)b;
		__m256i xidLoad;
		for (int i = 0; i < cntRem; i++)
			*(b + len - i - 1) = *(a + len - i - 1);
		for (int i = 0; i < cntBlock; i++)
		{
			xidLoad = _mm256_loadu_si256(p_new);
			_mm256_storeu_si256(p_old, xidLoad);
			p_old++; p_new++;
		}

	}
}
void replace_LLR(float* a, float* b, int len) {
	if (len < 8) {
		int cntBlock = len >> 2;
		int cntRem = len & 3;
		__m128 xidLoad;
		for (int i = 0; i < cntBlock; i++)
		{
			xidLoad = _mm_loadu_ps(a);
			_mm_storeu_ps(b, xidLoad);
			a += 4; b += 4;
		}
		for (int i = 0; i < cntRem; i++)
			*(b - i - 1) = *(a - i - 1);
	}
	else {
		int cntBlock = len >> 3;
		int cntRem = len & 7;
		__m256 xidLoad;
		for (int i = 0; i < cntBlock; i++)
		{
			xidLoad = _mm256_loadu_ps(a);
			_mm256_storeu_ps(b, xidLoad);
			a += 8; b += 8;
		}
		for (int i = 0; i < cntRem; i++) {
			*(b - i - 1) = *(a - i - 1);
		}

	}
}
void set_PM(float** a, float* b, int L, int len)
{
	if (len == 2) {
		b[0] = *(*(a)+len);
		b[1] = *(*(a + 1) + len);
	}
	else if (len == 4) {
		__m128 xidLoad;
		__m128 ZERO = _mm_setzero_ps();
		for (int i = 0; i < L; i += 4)
		{
			xidLoad = _mm_setr_ps(*(*(a + i) + len), *(*(a + i + 1) + len), *(*(a + i + 2) + len), *(*(a + i + 3) + len));
			_mm_storeu_ps(b + i, _mm_add_ps(_mm_loadu_ps(b + i), _mm_min_ps(xidLoad, ZERO)));
		}
	}
	else
	{
		__m256 xidLoad;
		__m256 ZERO = _mm256_setzero_ps();
		for (int i = 0; i < L; i += 8)
		{
			xidLoad = _mm256_setr_ps(*(*(a + i) + len), *(*(a + i + 1) + len), *(*(a + i + 2) + len), *(*(a + i + 3) + len), *(*(a + i + 4) + len), *(*(a + i + 5) + len), *(*(a + i + 6) + len), *(*(a + i + 7) + len));
			_mm256_storeu_ps(b + i, _mm256_add_ps(_mm256_loadu_ps(b + i), _mm256_min_ps(xidLoad, ZERO)));
		}
	}
}
int index_posit(int a, int b)
{
	int c;
	int temp = 0;
	c = a ^ b;
	if (c == 0) return 0;
	else {
		while (c != 1) {
			c >>= 1;
			temp++;
		}
		return temp;
	}
}
int lower_posit(int a)
{
	int c = a + 1;
	int temp = 0;
	while ((c & 1) == 0) {
		c >>= 1;
		temp++;
	}
	return temp;
}
void stage_located(int L, int* index, float** LLR, int** sum, int** p, int Count, int N, vector<int> path_2, vector<int> path_0)
{
	for (int k = 0; k < path_0.size(); k++) {
		int cnt;
		int temp = (path_0[k] < path_2[k]) ? p[path_0[k]][path_2[k]] : p[path_2[k]][path_0[k]];
		//cout << Count << " " << temp << endl;

		int level = index_posit(Count, temp);
		//cout << level << endl;
		//getchar();
		int lower = lower_posit(Count);
		cnt = (((Count >> level) - 1) << level);
		replace_sum(*(sum + path_2[k]) + cnt, *(sum + path_0[k]) + cnt, Count - cnt + 1);

		cnt = N - (2 << level);
		int lower_cnt = N - (2 << lower_posit(Count));
		replace_LLR(*(LLR + path_2[k]) + cnt + N, *(LLR + path_0[k]) + cnt + N, N - cnt);
	}
	for (int k = 0; k < path_0.size(); k++)
	{

		for (int i = 0; i < L; i++)
		{
			if (i != path_2[k]) {
				if (i < path_0[k])
				{
					p[i][path_0[k]] = (i < path_2[k]) ? p[i][path_2[k]] : p[path_2[k]][i];
				}
				else if (i > path_0[k])
				{
					p[path_0[k]][i] = (i < path_2[k]) ? p[i][path_2[k]] : p[path_2[k]][i];
				}
			}
		}
		p[path_2[k]][path_0[k]] = Count;
		p[path_0[k]][path_2[k]] = Count;
		for (int k2 = 0; k2 < path_0.size(); k2++)
		{
			if (path_0[k] < path_0[k2])p[path_0[k]][path_0[k2]] = (path_2[k2] < path_2[k]) ? p[path_2[k2]][path_2[k]] : p[path_2[k]][path_2[k2]];
			else if (k > k2)p[path_0[k2]][path_0[k]] = (path_2[k2] < path_2[k]) ? p[path_2[k2]][path_2[k]] : p[path_2[k]][path_2[k2]];
		}
	}
}


void PolarDecodeIterVertical(float** oriLLRh, float** upLLR, int** umid, int** G, int** uout_partial, int* A_Acv, int* Av, int Nh, int Nv, int decodingi, int iter)
{
	int len_decoded = Nh - decodingi;
	//for SCL Decoding 
	int Nmax = max(Nh, Nv);
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

	// Vertical and horizontal reverse-encoding
	if (system_archi == "JDD")
	{
		// vertical reverse-encoding
		PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), Nh);
		for (int j = 0; j < Kv; j++)
			umid[Av[j]][decodingi] = u1[index][Av[j]];
		// horizontal reverse-encoding
		/*cout << "len_decoded" << len_decoded << endl;
		for (int i = 0; i < Nv; i++)
		{
			PolarEncodePartial(G, umid[i], uout_partial[i], Nv, len_decoded);
		}*/
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


void Pyndiahh_bitflip(float** LLR, int** sum, float** oriLLRh, int* indexpm, int iter, int decodingi, float** upLLR, double& sdis_sum, int** tried_codewords_h, int** GGh, int* Ach, int* n_tried_patterns_h)
{
	int n_samples = (int)pow(2, p);
	int num_first_diff = 0;
	float* reliability_metric = new float[Nh];
	int* first_diff_array = new int[Nh];
	int* idx_first_diff_sort = new int[Nh];
	float* sdis_sampled = new float[n_samples];
	int* valid_sample = new int[n_samples];
	int* bit_tmp = new int[Nh];
	int* bit_tmp_encoded = new int[Nh];
	vector<int> valid_samples_idx;
	//int** tried_codeword_set = new int[]
	for (int i = 0; i < n_samples; i++) { sdis_sampled[i] = 0; valid_sample[i] = 1; }
	for (int i = 0; i < Nh; i++) { idx_first_diff_sort[i] = -10 * Nh; }
	for (int i = 0; i < Nh; i++) { first_diff_array[i] = 0; bit_tmp[i] = 0; bit_tmp_encoded[i] = 0; }
	// decide bits with different decisions 
	for (int j = 0; j < Nh; j++)
	{
		int firstdiff = 0;
		for (int i = 1; i < L; i++)
			if (sum[indexpm[i]][j] != sum[indexpm[0]][j])
			{
				first_diff_array[j] = i;
				num_first_diff++;
				break;
			}
	}
	// use (2x-1)*llr to represent the reliability. Higher means less reliable
	for (int i = 0; i < Nh; i++) { reliability_metric[i] = (2*sum[indexpm[0]][i]-1) * LLR[0][i]; }
	// arrange the positions of first error in ascending order to find the p bits with least certain decisions
	// my_sort_pyndiah(first_diff_array, idx_first_diff_sort, Nh);
	/*cout << endl << "reliability_metric" << endl;
	for (int i = 0; i < Nh; i++) cout << reliability_metric[i] << "   ";*/
	my_sort_pyndiah(reliability_metric, idx_first_diff_sort, Nh);
	// test: reliability_metric and idx
	/*cout << endl << "Best Codeword" << endl;
	for (int i = 0; i < Nh; i++) cout << sum[indexpm[0]][i] << "   ";
	cout << endl << "idx" << endl;
	for (int i = 0; i < Nh; i++) cout << idx_first_diff_sort[i] << "   ";
	cout << endl;*/
	// test: reliability_metric and idx
	float sdis[L];
	for (int i = 0; i < L; i++)
	{
		sdis[i] = 0;
		for (int j = 0; j < Nh; j++)
			sdis[i] = sdis[i] + (1 - sign(LLR[0][j]) * (1 - 2 * sum[i][j])) * abs(LLR[0][j]) / 2;
	}
	for (int i = 0; i < L; i++)
	{
		sdis_sum += sdis[i];
	}
	int num_firstdiff = 0;
	

	// test: LSD+Random Sampling
	// calculate the Euclidian distance of each sample.
	int i_code = 0;
	int n_valid_samples = 0;
	for (int i = 0; i < n_samples; i++)
	{
		// cout << "i=" << i << endl;
		for (int ibit = 0; ibit < Nh; ibit++) {
			bit_tmp[ibit] = sum[indexpm[i_code]][ibit];
		}
		// bitset: 0 (LSB) to 19 (MSB)
		sdis_sampled[i] = sdis[indexpm[i_code]];
		std::bitset<20> candidate_bits(i+1);
		// cout << candidate_bits[0];
		for (int j = 0; j < p; j++)
		{
			int bit_flip_position = idx_first_diff_sort[j];
			bool bflip = (int)candidate_bits[j];
			//cout <<"BFLIPPPP" << bflip << endl;
			if (!bflip) { 
				continue; 
			}
			else
			{
				// determined which bit is being flipped
				// cout << bit_flip_position << endl;
				int flipped_bit = (sum[indexpm[i_code]][bit_flip_position]<0.5);
				
				bit_tmp[bit_flip_position] = flipped_bit;
				sdis_sampled[i] -= (1 - sign(LLR[0][bit_flip_position]) * (1 - 2 * sum[indexpm[i_code]][bit_flip_position])) * abs(LLR[0][bit_flip_position]) / 2;
				sdis_sampled[i] += (1 - sign(LLR[0][bit_flip_position]) * (1 - 2 * flipped_bit)) * abs(LLR[0][bit_flip_position]) / 2;
			}
		}
		if (!test_sc)
		{
			if (sdis_sampled[i] > max_llr_ratio * sdis[indexpm[L - 1]]) {
				valid_sample[i] = 0;
			}
		}
		// paths with the same distance are not valid
		for (int j = 0; j < L; j++)
		{
			if (abs(sdis_sampled[i] - sdis[j]) < 1E-5)
			{
				valid_sample[i] = 0;
			}
		}
		if (valid_sample[i]|disable_max_ratio)
		{
			PolarEncode_xor(bit_tmp_encoded, bit_tmp, Nh);
			bool is_cwd = true;
			for (int j = 0; j < Nh - Kh; j++)
			{
				if (bit_tmp_encoded[Ach[j]] != 0) { valid_sample[i] = 0; is_cwd = false; }
			}
			if (is_cwd) { n_tried_patterns_h[decodingi]++; }
		}
		
		if (valid_sample[i]) {
			n_valid_samples++;
			valid_samples_idx.push_back(i);
			for (int j = 0; j < Nh; j++)
				tried_codewords_h[i][j] = bit_tmp[j];
		}
	}

	int* indexpm_bitflip = new int[n_valid_samples+L];
	for (int i = 0; i < n_valid_samples + L; i++) { indexpm_bitflip[i] = i; }
	int** sum_bitflip = new int* [n_valid_samples + L];
	float* sdis_all = new float[n_valid_samples + L];
	for (int i = 0; i < (n_valid_samples + L); i++) { sum_bitflip[i] = new int[Nh]; }
	for (int i = 0; i < (n_valid_samples+L); i++)
	{
		for (int j = 0; j < Nh; j++) {
			if (i < L) {
				sum_bitflip[i][j] = sum[i][j];
				sdis_all[i] = sdis[i];
			}
			else {
				sum_bitflip[i][j] = tried_codewords_h[valid_samples_idx[i-L]][j];
				sdis_all[i] = sdis_sampled[valid_samples_idx[i - L]];
			}
		}
	}
	// arrange indexpm_bitflip
	for (int i = 0; i < (n_valid_samples+L-1); i++)
		for (int j = i + 1; j < (n_valid_samples+L); j++)
			if (sdis_all[indexpm_bitflip[i]] > sdis_all[indexpm_bitflip[j]])
			{
				int ttt = indexpm_bitflip[i];
				indexpm_bitflip[i] = indexpm_bitflip[j];
				indexpm_bitflip[j] = ttt;
			}
	// output indexpm_bitflip, sdis_all, sum_bitflip
	int test_mode = 0;
	// cout << "N_VALID_SAMPLES=" << (n_valid_samples + L) << endl;
	if (test_mode) {
		cout << "indexpm_bitflip" << endl;
		for (int i = 0; i < n_valid_samples + L; i++)
		{
			cout << indexpm_bitflip[i] << "   ";
		}
		cout << endl;
		cout << "sdis_all" << endl;
		for (int i = 0; i < n_valid_samples + L; i++)
		{
			cout << sdis_all[indexpm_bitflip[i]] << "   ";
		}
		cout << endl;
		cout << "sum_bitflip" << endl;
		for (int i = 0; i < n_valid_samples + L; i++)
		{
			for (int j = 0; j < Nh; j++)
			{
				cout << sum_bitflip[indexpm_bitflip[i]][j] << "   ";
			}
			cout << endl;
		}
		cout << endl;
	}
	// cout << sdis_sum << endl;
	////TEST:
	//cout << "Pyndiah Dist" << endl;
	//cout << endl;
	//for (int i = 0; i < L; i++)
	//{
	//	cout << sdis[i] << "   ";
	//}
	//cout << endl;
	//cout << "Pyndiah" << endl;

	//cout << "Appended Dist" << endl;
	//for (int i = 0; i < n_samples; i++)
	//{
	//	if (valid_sample[i] && sdis_sampled[i]<sdis[L-1]) cout << sdis_sampled[i] << "     ";
	//}
	//cout << "Appended Dist" << endl;
	////TEST
	//// Test Magnitude Restrict
	//for (int i = 0; i < L; i++)
	//{
	//	if (sdis[i] > 10000) sdis[i] = 10000;
	//	if (sdis[i] < -10000) sdis[i] = -10000;
	//}
	/*if (showllr)
	{
		cout << "Pyndiah_LSD" << "iter=  " << iter << endl;
		for (int i = 0; i < L; i++) cout << sdis[i] << "   ";
		cout << endl;
	}*/
	// Test Magnitude Restrict 
	// n of paths considered
	int L_CONSIDERED = L;
	if (L < 8) L_CONSIDERED = L + n_valid_samples;
	if (L_CONSIDERED > L_max) L_CONSIDERED = L_max;
	float w[Nh];
	if (usebeta)
	{
		for (int j = 0; j < Nv; j++)
		{
			if (upLLR[j][decodingi] < -beta[iter * 2 - 2])
				upLLR[j][decodingi] = -beta[iter * 2 - 2];
			if (upLLR[j][decodingi] > beta[iter * 2 - 2])
				upLLR[j][decodingi] = beta[iter * 2 - 2];
		}
	}
	for (int j = 0; j < Nh; j++)
	{
		int firstdiff = 0;
		for (int i = 1; i < L_CONSIDERED; i++)
			if (sum_bitflip[indexpm_bitflip[i]][j] != sum_bitflip[indexpm_bitflip[0]][j])
			{
				firstdiff = i;
				break;
			}
		if (firstdiff)
			//if (false)
		{
			num_firstdiff++;
			w[j] = (sdis_all[indexpm_bitflip[firstdiff]] - sdis_all[indexpm_bitflip[0]]) * (1 - 2 * sum_bitflip[indexpm_bitflip[0]][j]);
			if (usebeta)
			{
				if (w[j] < -beta[iter * 2 - 1])
					w[j] = -beta[iter * 2 - 1];
				else if (w[j] > beta[iter * 2 - 1])
					w[j] = beta[iter * 2 - 1];
			}
		}
		else
			if (usebeta) w[j] = beta[iter * 2 - 1] * (1 - 2 * sum_bitflip[indexpm_bitflip[0]][j]);
		//else w[j] = (sdis[indexpm[L - 1]] - sdis[indexpm[0]]) * (1 - 2 * sum[indexpm[0]][j]);
			else w[j] = (sdis_all[indexpm_bitflip[L_CONSIDERED - 1]] - sdis_all[indexpm_bitflip[0]]) * (1 - 2 * sum_bitflip[indexpm_bitflip[0]][j]);
		upLLR[decodingi][j] = oriLLRh[decodingi][j] + alpha[iter * 2 - 1] * (w[j] - upLLR[decodingi][j] * (1.0));
		// test w, upLLR, oriLLR

		// test w, upLLR, oriLLR
		//upLLR[j][decodingi] = oriLLRh[j][decodingi] + alpha[iter * 2 - 1] * (w[j]);
		if (usebeta)
		{
			for (int j = 0; j < Nv; j++)
			{
				if (upLLR[j][decodingi] < -beta[iter * 2 - 1])
					upLLR[j][decodingi] = -beta[iter * 2 - 1];
				if (upLLR[j][decodingi] > beta[iter * 2 - 1])
					upLLR[j][decodingi] = beta[iter * 2 - 1];
			}
		}
		//upLLR[decodingi][j] = oriLLRh[decodingi][j] + 0.3 * (w[j] - upLLR[decodingi][j]);
		//	cout << upLLR[decodingi][j] << " ";
		// upLLR[decodingi][j] = oriLLRh[decodingi][j] + alpha[iter * 2 - 1] * (w[j] - upLLR[decodingi][j]);
	}
	// cout << "NUM_FIRSTDIFF:  " << num_firstdiff << endl;
	delete[] sdis_sampled;
	delete[] idx_first_diff_sort;
	delete[] first_diff_array;
	delete[] valid_sample;
	delete[] bit_tmp;
	delete[] bit_tmp_encoded;
	for (int i = 0; i < (n_valid_samples + L); i++)
	{
		delete[] sum_bitflip[i];
	}
	// cout << n_valid_samples << endl;
	delete[] sum_bitflip;
	delete[] indexpm_bitflip;
	delete[] sdis_all;
	delete[] reliability_metric;
}



void Pyndiahv_bitflip(float** LLR, int** sum, float** oriLLRh, int* indexpm, int iter, int decodingi, float** upLLR, double& sdis_sum, int** tried_codewords_v, int** GGv, int* Acv, int* n_tried_patterns_v)
{
	/*int p_time = 9;
	int dull_var = 0;
	p_time = p_time + dull_var;*/
	int n_samples = (int)pow(2, p);
	int num_first_diff = 0;
	float* reliability_metric = new float[Nv];
	int* first_diff_array = new int[Nv];
	int* idx_first_diff_sort = new int[Nv];
	float* sdis_sampled = new float[n_samples];
	int* valid_sample = new int[n_samples];
	int* bit_tmp = new int[Nv];
	int* bit_tmp_encoded = new int[Nv];
	vector<int> valid_samples_idx;
	//int** tried_codeword_set = new int[]
	for (int i = 0; i < n_samples; i++) { sdis_sampled[i] = 0; valid_sample[i] = 1; }
	for (int i = 0; i < Nv; i++) { idx_first_diff_sort[i] = -10 * Nv; }
	for (int i = 0; i < Nv; i++) { first_diff_array[i] = 0; bit_tmp[i] = 0; bit_tmp_encoded[i] = 0; }
	// decide bits with different decisions 
	for (int j = 0; j < Nv; j++)
	{
		int firstdiff = 0;
		for (int i = 1; i < L; i++)
			if (sum[indexpm[i]][j] != sum[indexpm[0]][j])
			{
				first_diff_array[j] = i;
				num_first_diff++;
				break;
			}
	}
	// use (2x-1)*llr to represent the reliability. Higher means less reliable
	for (int i = 0; i < Nv; i++) { reliability_metric[i] = (2 * sum[indexpm[0]][i] - 1) * LLR[0][i]; }
	// arrange the positions of first error in ascending order to find the p bits with least certain decisions
	// my_sort_pyndiah(first_diff_array, idx_first_diff_sort, Nh);
	/*cout << endl << "reliability_metric" << endl;
	for (int i = 0; i < Nh; i++) cout << reliability_metric[i] << "   ";*/
	my_sort_pyndiah(reliability_metric, idx_first_diff_sort, Nv);
	// test: reliability_metric and idx
	/*cout << endl << "Best Codeword" << endl;
	for (int i = 0; i < Nh; i++) cout << sum[indexpm[0]][i] << "   ";
	cout << endl << "idx" << endl;
	for (int i = 0; i < Nh; i++) cout << idx_first_diff_sort[i] << "   ";
	cout << endl;*/
	// test: reliability_metric and idx
	float sdis[L];
	for (int i = 0; i < L; i++)
	{
		sdis[i] = 0;
		for (int j = 0; j < Nv; j++)
			sdis[i] = sdis[i] + (1 - sign(LLR[0][j]) * (1 - 2 * sum[i][j])) * abs(LLR[0][j]) / 2;
	}
	for (int i = 0; i < L; i++)
	{
		sdis_sum += sdis[i];
	}
	int num_firstdiff = 0;


	// test: LSD+Random Sampling
	// calculate the Euclidian distance of each sample.
	int i_code = 0;
	int n_valid_samples = 0;
	for (int i = 0; i < n_samples; i++)
	{
		// cout << "i=" << i << endl;
		for (int ibit = 0; ibit < Nv; ibit++) {
			bit_tmp[ibit] = sum[indexpm[i_code]][ibit];
		}
		// bitset: 0 (LSB) to 19 (MSB)
		sdis_sampled[i] = sdis[indexpm[i_code]];
		std::bitset<20> candidate_bits(i + 1);
		// cout << candidate_bits[0];
		for (int j = 0; j < p; j++)
		{
			int bit_flip_position = idx_first_diff_sort[j];
			bool bflip = (int)candidate_bits[j];
			//cout <<"BFLIPPPP" << bflip << endl;
			if (!bflip) {
				continue;
			}
			else
			{
				// determined which bit is being flipped
				// cout << bit_flip_position << endl;
				int flipped_bit = (sum[indexpm[i_code]][bit_flip_position] < 0.5);

				bit_tmp[bit_flip_position] = flipped_bit;
				sdis_sampled[i] -= (1 - sign(LLR[0][bit_flip_position]) * (1 - 2 * sum[indexpm[i_code]][bit_flip_position])) * abs(LLR[0][bit_flip_position]) / 2;
				sdis_sampled[i] += (1 - sign(LLR[0][bit_flip_position]) * (1 - 2 * flipped_bit)) * abs(LLR[0][bit_flip_position]) / 2;
			}
		}
		if (!test_sc)
		{
			if (sdis_sampled[i] > max_llr_ratio * sdis[indexpm[L - 1]]) {
				valid_sample[i] = 0;
			}
		}
		// paths with the same distance are not valid
		for (int j = 0; j < L; j++)
		{
			if (abs(sdis_sampled[i] - sdis[j]) < 1E-5)
			{
				valid_sample[i] = 0;
			}
		}
		if (valid_sample[i] | disable_max_ratio)
		{
			PolarEncode_xor(bit_tmp_encoded, bit_tmp, Nv);
			bool is_cwd = true;
			for (int j = 0; j < Nv - Kv; j++)
			{
				if (bit_tmp_encoded[Acv[j]] != 0) { valid_sample[i] = 0; is_cwd = false; }
			}
			if (is_cwd) { n_tried_patterns_v[decodingi]++; }
		}

		if (valid_sample[i]) {
			n_valid_samples++;
			valid_samples_idx.push_back(i);
			for (int j = 0; j < Nv; j++)
				tried_codewords_v[i][j] = bit_tmp[j];
		}
	}

	int* indexpm_bitflip = new int[n_valid_samples + L];
	for (int i = 0; i < n_valid_samples + L; i++) { indexpm_bitflip[i] = i; }
	int** sum_bitflip = new int* [n_valid_samples + L];
	float* sdis_all = new float[n_valid_samples + L];
	for (int i = 0; i < (n_valid_samples + L); i++) { sum_bitflip[i] = new int[Nv]; }
	for (int i = 0; i < (n_valid_samples + L); i++)
	{
		for (int j = 0; j < Nv; j++) {
			if (i < L) {
				sum_bitflip[i][j] = sum[i][j];
				sdis_all[i] = sdis[i];
			}
			else {
				sum_bitflip[i][j] = tried_codewords_v[valid_samples_idx[i - L]][j];
				sdis_all[i] = sdis_sampled[valid_samples_idx[i - L]];
			}
		}
	}
	// arrange indexpm_bitflip
	for (int i = 0; i < (n_valid_samples + L - 1); i++)
		for (int j = i + 1; j < (n_valid_samples + L); j++)
			if (sdis_all[indexpm_bitflip[i]] > sdis_all[indexpm_bitflip[j]])
			{
				int ttt = indexpm_bitflip[i];
				indexpm_bitflip[i] = indexpm_bitflip[j];
				indexpm_bitflip[j] = ttt;
			}
	// output indexpm_bitflip, sdis_all, sum_bitflip
	int test_mode = 0;
	// cout << "N_VALID_SAMPLES=" << (n_valid_samples + L) << endl;
	if (test_mode) {
		cout << "indexpm_bitflip" << endl;
		for (int i = 0; i < n_valid_samples + L; i++)
		{
			cout << indexpm_bitflip[i] << "   ";
		}
		cout << endl;
		cout << "sdis_all" << endl;
		for (int i = 0; i < n_valid_samples + L; i++)
		{
			cout << sdis_all[indexpm_bitflip[i]] << "   ";
		}
		cout << endl;
		cout << "sum_bitflip" << endl;
		for (int i = 0; i < n_valid_samples + L; i++)
		{
			for (int j = 0; j < Nv; j++)
			{
				cout << sum_bitflip[indexpm_bitflip[i]][j] << "   ";
			}
			cout << endl;
		}
		cout << endl;
	}
	// cout << sdis_sum << endl;
	////TEST:
	//cout << "Pyndiah Dist" << endl;
	//cout << endl;
	//for (int i = 0; i < L; i++)
	//{
	//	cout << sdis[i] << "   ";
	//}
	//cout << endl;
	//cout << "Pyndiah" << endl;

	//cout << "Appended Dist" << endl;
	//for (int i = 0; i < n_samples; i++)
	//{
	//	if (valid_sample[i] && sdis_sampled[i]<sdis[L-1]) cout << sdis_sampled[i] << "     ";
	//}
	//cout << "Appended Dist" << endl;
	////TEST
	//// Test Magnitude Restrict
	//for (int i = 0; i < L; i++)
	//{
	//	if (sdis[i] > 10000) sdis[i] = 10000;
	//	if (sdis[i] < -10000) sdis[i] = -10000;
	//}
	/*if (showllr)
	{
		cout << "Pyndiah_LSD" << "iter=  " << iter << endl;
		for (int i = 0; i < L; i++) cout << sdis[i] << "   ";
		cout << endl;
	}*/
	// Test Magnitude Restrict 
	// n of paths considered
	int L_CONSIDERED = L;
	if (L < 8) L_CONSIDERED = L + n_valid_samples;
	if (L_CONSIDERED > L_max) L_CONSIDERED = L_max;
	float w[Nh];
	if (usebeta)
	{
		for (int j = 0; j < Nv; j++)
		{
			if (upLLR[j][decodingi] < -beta[iter * 2 - 2])
				upLLR[j][decodingi] = -beta[iter * 2 - 2];
			if (upLLR[j][decodingi] > beta[iter * 2 - 2])
				upLLR[j][decodingi] = beta[iter * 2 - 2];
		}
	}
	for (int j = 0; j < Nv; j++)
	{
		int firstdiff = 0;
		for (int i = 1; i < L_CONSIDERED; i++)
			if (sum_bitflip[indexpm_bitflip[i]][j] != sum_bitflip[indexpm_bitflip[0]][j])
			{
				firstdiff = i;
				break;
			}
		if (firstdiff)
			//if (false)
		{
			num_firstdiff++;
			w[j] = (sdis_all[indexpm_bitflip[firstdiff]] - sdis_all[indexpm_bitflip[0]]) * (1 - 2 * sum_bitflip[indexpm_bitflip[0]][j]);
			if (usebeta)
			{
				if (w[j] < -beta[iter * 2 - 1])
					w[j] = -beta[iter * 2 - 1];
				else if (w[j] > beta[iter * 2 - 1])
					w[j] = beta[iter * 2 - 1];
			}
		}
		else
			if (usebeta) w[j] = beta[iter * 2 - 1] * (1 - 2 * sum_bitflip[indexpm_bitflip[0]][j]);
		//else w[j] = (sdis[indexpm[L - 1]] - sdis[indexpm[0]]) * (1 - 2 * sum[indexpm[0]][j]);
			else w[j] = (sdis_all[indexpm_bitflip[L_CONSIDERED - 1]] - sdis_all[indexpm_bitflip[0]]) * (1 - 2 * sum_bitflip[indexpm_bitflip[0]][j]);
		upLLR[j][decodingi] = oriLLRh[j][decodingi] + alpha[iter * 2 - 1] * (w[j] - upLLR[j][decodingi] * (1.0));
		// test w, upLLR, oriLLR
		// test w, upLLR, oriLLR
		//upLLR[j][decodingi] = oriLLRh[j][decodingi] + alpha[iter * 2 - 1] * (w[j]);
		if (usebeta)
		{
			for (int j = 0; j < Nv; j++)
			{
				if (upLLR[j][decodingi] < -beta[iter * 2 - 1])
					upLLR[j][decodingi] = -beta[iter * 2 - 1];
				if (upLLR[j][decodingi] > beta[iter * 2 - 1])
					upLLR[j][decodingi] = beta[iter * 2 - 1];
			}
		}
		//upLLR[decodingi][j] = oriLLRh[decodingi][j] + 0.3 * (w[j] - upLLR[decodingi][j]);
		//	cout << upLLR[decodingi][j] << " ";
		// upLLR[decodingi][j] = oriLLRh[decodingi][j] + alpha[iter * 2 - 1] * (w[j] - upLLR[decodingi][j]);
	}
	// cout << "NUM_FIRSTDIFF:  " << num_firstdiff << endl;
    
	delete[] sdis_sampled;
	delete[] idx_first_diff_sort;
	delete[] first_diff_array;
	delete[] valid_sample;
	delete[] bit_tmp;
	delete[] bit_tmp_encoded;
	for (int i = 0; i < (n_valid_samples + L); i++)
	{
		delete[] sum_bitflip[i];
	}
	// cout << n_valid_samples << endl;
	delete[] sum_bitflip;
	delete[] indexpm_bitflip;
	delete[] sdis_all;
	delete[] reliability_metric;
}

void Pyndiahh_irregular(float** LLR, int** sum, float** oriLLRh, int* indexpm, int iter, int decodingi, float** upLLR, double& sum_sdis, int* irregular_construction)
{
	float sdis[L];
	for (int i = 0; i < L; i++)
	{
		sdis[i] = 0;
		for (int j = 0; j < Nh; j++)
			sdis[i] = sdis[i] + (1 - sign(LLR[0][j]) * (1 - 2 * sum[i][j])) * abs(LLR[0][j]) / 2;
	}

	for (int i = 0; i < L; i++) { sum_sdis += sdis[i]; }
	////TEST:
	//cout << "Pyndiah Dist" << endl;
	//cout << endl;
	//for (int i = 0; i < L_LSD; i++)
	//{
	//	cout << sdis[i] << "   ";
	//}
	//cout << endl;
	//cout << "Pyndiah" << endl;
	////TEST
	bool show_tmp = ((decodingi == 3) && (showllr_irr));
	float upLLR_pre = 0.0;
	/*if (show_tmp) {
		std::cout << endl << setw(10) << iter;
	}*/

	float w[Nh];
	for (int j = 0; j < Nh; j++)
	{
		upLLR_pre = upLLR[decodingi][j];
		int firstdiff = 0;
		for (int i = 1; i < L; i++)
			if (sum[indexpm[i]][j] != sum[indexpm[0]][j])
			{
				firstdiff = i;
				break;
			}
		if (firstdiff)
			w[j] = (sdis[indexpm[firstdiff]] - sdis[indexpm[0]]) * (1 - 2 * sum[indexpm[0]][j]);
		else
			if (usebeta) w[j] = beta[iter * 2 - 1] * (1 - 2 * sum[indexpm[0]][j]);
			else w[j] = (sdis[indexpm[L - 1]] - sdis[indexpm[0]]) * (1 - 2 * sum[indexpm[0]][j]);
		if (irregular_construction[decodingi]==1)
			upLLR[decodingi][j] = oriLLRh[decodingi][j] + alpha[iter * 2 - 1]* (amplify_ratio_array[iter-1]*w[j] - upLLR[decodingi][j] * (1.0 - usebeta));
		else
			upLLR[decodingi][j] = oriLLRh[decodingi][j] + alpha[iter * 2 - 1] * (w[j] - upLLR[decodingi][j] * (1.0 - usebeta));
		//upLLR[j][decodingi] = oriLLRh[j][decodingi] + alpha[iter * 2 - 1] * (w[j]);
		//	cout << upLLR[decodingi][j] << " ";
		// for test; display w[j] and upLLR[decodingi][j];
		/*if (show_tmp) { std::cout << setw(20) << w[j] << " " << setw(20) << upLLR[decodingi][j]; }*/
	}
}