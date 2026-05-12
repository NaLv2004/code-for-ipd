# include "modulation.h"
# include <cmath>
# include <algorithm>
using namespace std;


float table_qpsk[2] = { -0.707107f, 0.707107f };
float table_16qam[4] = { -0.316228f, -0.948683f,  0.316228f,  0.948683f };
//float table_16qam[4]={}
float table_64qam[8] = { -0.462910f, -0.154303f, -0.771517f, -1.08012f,	0.462910f,  0.154303f,  0.771517f,  1.08012f };
float table_256qam[16] =
{ -0.383482f, -0.536875f, -0.230089f, -0.07669f,
-0.843661f, -0.690268f, -0.997054f, -1.15044f,
0.383482f, 0.536875f, 0.230089f, 0.07669f,
0.843661f, 0.690268f, 0.997054f, 1.15044f };

float s_bits_r2[2][1] = { {0}, {1} };
float s_bits_i2[2][1] = { {0}, {1} };
float s_bits_r4[4][2] = { {0,0},{0,1},{1,0},{1,1} };
float s_bits_i4[4][2] = { {0,0},{0,1},{1,0},{1,1} };
float s_bits_r6[8][3] = { {0,0,0},{0,0,1},{0,1,0},{0,1,1},{1,0,0},{1,0,1},{1,1,0},{1,1,1} };
float s_bits_i6[8][3] = { {0,0,0},{0,0,1},{0,1,0},{0,1,1},{1,0,0},{1,0,1},{1,1,0},{1,1,1} };
float s_bits_r8[16][4] = { {0,0,0,0},{0,0,0,1},{0,0,1,0},{0,0,1,1},{0,1,0,0},{0,1,0,1},{0,1,1,0},{0,1,1,1},{1,0,0,0},{1,0,0,1},{1,0,1,0},{1,0,1,1},{1,1,0,0},{1,1,0,1},{1,1,1,0},{1,1,1,1} };
float s_bits_i8[16][4] = { {0,0,0,0},{0,0,0,1},{0,0,1,0},{0,0,1,1},{0,1,0,0},{0,1,0,1},{0,1,1,0},{0,1,1,1},{1,0,0,0},{1,0,0,1},{1,0,1,0},{1,0,1,1},{1,1,0,0},{1,1,0,1},{1,1,1,0},{1,1,1,1} };

void modulation(int* input_bits, std::complex<float>* symbols_out, int mod_type, int nof_symbols) {
	int i, j;
	int tmp_i, tmp_q;
	int half_sym;

	half_sym = mod_type / 2;

	switch (mod_type)
	{
	case 2:
		for (i = 0; i < nof_symbols; i++)
		{
			tmp_i = 0;
			tmp_q = 0;

			for (j = 0; j < half_sym; j++)
			{
				tmp_i = tmp_i + (input_bits[i * mod_type + j] << (half_sym - j - 1));
				tmp_q = tmp_q + (input_bits[i * mod_type + half_sym + j] << (half_sym - j - 1));
			}

			float sym_real = table_qpsk[tmp_i];
			float sym_imag = table_qpsk[tmp_q];
			symbols_out[i] = std::complex<float>(sym_real, sym_imag);
		}
		break;
	case 4:
		for (i = 0; i < nof_symbols; i++)
		{
			// std::cout << nof_symbols << std::endl;
			// TEST:
			// cout << i << endl;
			// TEST
			tmp_i = 0;
			tmp_q = 0;
			for (j = 0; j < half_sym; j++)
			{
				// std::cout << input_bits[i * mod_type + j] << std::endl;
				tmp_i = tmp_i + (input_bits[i * mod_type + j] << (half_sym - j - 1));
				tmp_q = tmp_q + (input_bits[i * mod_type + half_sym + j] << (half_sym - j - 1));
			}
			/*std::cout << tmp_i << endl;
			std::cout << tmp_q << endl;*/
			//	symbols_out[i].real = table_16qam[tmp_i];
			//	symbols_out[i].imag = table_16qam[tmp_q];
			float sym_real = table_16qam[tmp_i];
			float sym_imag = table_16qam[tmp_q];
			/*cout << tmp_i << endl;
			cout << tmp_q << endl;*/
			symbols_out[i] = std::complex<float>(sym_real, sym_imag);
		}
		break;
	case 6:
		for (i = 0; i < nof_symbols; i++)
		{
			tmp_i = 0;
			tmp_q = 0;

			for (j = 0; j < half_sym; j++)
			{
				tmp_i = tmp_i + (input_bits[i * mod_type + j] << (half_sym - j - 1));
				tmp_q = tmp_q + (input_bits[i * mod_type + half_sym + j] << (half_sym - j - 1));
			}

			float sym_real = table_64qam[tmp_i];
			float sym_imag = table_64qam[tmp_q];
			symbols_out[i] = std::complex<float>(sym_real, sym_imag);

		}
		break;

	case 8:
		for (i = 0; i < nof_symbols; i++)
		{
			tmp_i = 0;
			tmp_q = 0;

			for (j = 0; j < half_sym; j++)
			{
				tmp_i = tmp_i + (input_bits[i * mod_type + j] << (half_sym - j - 1));
				tmp_q = tmp_q + (input_bits[i * mod_type + half_sym + j] << (half_sym - j - 1));
			}

			//	symbols_out[i].real = table_256qam[tmp_i];
			//	symbols_out[i].imag = table_256qam[tmp_q];
			float sym_real = table_256qam[tmp_i];
			float sym_imag = table_256qam[tmp_q];
			symbols_out[i] = std::complex<float>(sym_real, sym_imag);
		}
		break;
	default:
		printf("Invalid modulation type: 2:QPSK, 4:16-QAM \n");
	}
	//	printf("\n**********\n");
}

// sym output to soft bitwise llr output
void llr_mmse_bit(int Nt, float* u_sym_mmse, float* llr_bit, float sigma_sq, int ModType) {
	float prob1 = 0;
	float prob0 = 0;
	float norm_factor = 0;

	int Nt2 = 2 * Nt;
	int ModType_half = ModType / 2;
	int Q = pow(2, ModType_half);

	float* u_sym_reverse = new float[Nt2];
	float* sym_prob = new float[Q];   // prob of every half symbol in the alphabet


	for (int i = 0; i < Nt; i++)
	{
		u_sym_reverse[2 * i] = u_sym_mmse[i]; // real part
		u_sym_reverse[2 * i + 1] = u_sym_mmse[i + Nt]; // imag part
	}
	switch (ModType)
	{
	case 2:
		for (int i = 0; i < Nt2; i++)
		{
			// prob of every symbol in the alphabet tablee p(u) = exp[-(x-u_mmse)^2/(2*sigma_sq)]
			//norm_factor = 0;
			float sym_est = u_sym_reverse[i];
			// test: boundary decision
			llr_bit[i * ModType_half + 0] = -sym_est;
			llr_bit[i * ModType_half + 1] = -sym_est;
		}

		break;
	case 4:
		for (int i = 0; i < Nt2; i++)
		{
			// prob of every symbol in the alphabet tablee p(u) = exp[-(x-u_mmse)^2/(2*sigma_sq)]
			//norm_factor = 0;
			float sym_est = u_sym_reverse[i];
			// test: boundary decision
			llr_bit[i * ModType_half + 0] = -sym_est;
			llr_bit[i * ModType_half + 1] = -abs(sym_est) + 0.632455;
		}

		break;
	case 6:
		for (int i = 0; i < Nt2; i++)
		{
			// prob of every symbol in the alphabet tablee p(u) = exp[-(x-u_mmse)^2/(2*sigma_sq)]
			norm_factor = 0;
			float sym_est = u_sym_reverse[i];
			for (int i_sym = 0; i_sym < Q; i_sym++)
			{
				float sym_val = table_64qam[i_sym];
				float dist = abs(sym_val - sym_est);
				float dist2_nv = pow(dist, 2) / (2.0 * sigma_sq);
				sym_prob[i_sym] = exp(-dist2_nv);
				norm_factor += sym_prob[i_sym];
			}
			// normalization
			for (int i_sym = 0; i_sym < Q; i_sym++) { sym_prob[i_sym] /= norm_factor; }
			// p1 and p0
			for (int i_bit = 0; i_bit < ModType_half; i_bit++)
			{
				prob0 = 0;
				prob1 = 0;
				for (int i_sym = 0; i_sym < Q; i_sym++)
				{
					if (s_bits_r6[i_sym][i_bit] == 0) { prob0 += sym_prob[i_sym]; }
					if (s_bits_r6[i_sym][i_bit] == 1) { prob1 += sym_prob[i_sym]; }
				}
				// llr = log (p(0) / p(1)) 
				float logprob = log(prob0 / prob1);
				if (logprob > 10000) logprob = 10000;
				if (logprob < -10000) logprob = -10000;
				llr_bit[i * ModType_half + i_bit] = logprob;
			}
		}
		break;
	case 8:
		for (int i = 0; i < Nt2; i++)
		{
			// prob of every symbol in the alphabet tablee p(u) = exp[-(x-u_mmse)^2/(2*sigma_sq)]
			norm_factor = 0;
			float sym_est = u_sym_reverse[i];
			for (int i_sym = 0; i_sym < Q; i_sym++)
			{
				float sym_val = table_256qam[i_sym];
				float dist = abs(sym_val - sym_est);
				float dist2_nv = -pow(dist, 2) / (2.0 * sigma_sq);
				//sym_prob[i_sym] = exp(-dist2_nv);
				sym_prob[i_sym] = dist2_nv;
				norm_factor += sym_prob[i_sym];
			}
			// normalization
			// for (int i_sym = 0; i_sym < Q; i_sym++) { sym_prob[i_sym] /= norm_factor;  cout << sym_prob[i_sym] << "  "; }
			// p1 and p0
			for (int i_bit = 0; i_bit < ModType_half; i_bit++)
			{
				prob0 = 0;
				prob1 = 0;
				/*for (int i_sym = 0; i_sym < Q; i_sym++)
				{
					if (s_bits_r4[i_sym][i_bit] == 0) { prob0 += sym_prob[i_sym]; }
					if (s_bits_r4[i_sym][i_bit] == 1) { prob1 += sym_prob[i_sym]; }
				}*/
				float logprob = 0;
				float max_p0 = 0;
				float max_p1 = 0;
				int cnt_p0 = 0;
				int cnt_p1 = 0;
				for (int i_sym_prob = 0; i_sym_prob < Q; i_sym_prob++)
				{
					if (s_bits_r8[i_sym_prob][i_bit] == 0)
					{
						if (cnt_p0 == 0) { max_p0 = sym_prob[i_sym_prob]; }
						else
						{
							if (max_p0 < sym_prob[i_sym_prob])
								max_p0 = sym_prob[i_sym_prob];
						}
						cnt_p0++;
					}
					else
					{
						if (cnt_p1 == 0) { max_p1 = sym_prob[i_sym_prob]; }
						else
						{
							if (max_p1 < sym_prob[i_sym_prob])
								max_p1 = sym_prob[i_sym_prob];
						}
						cnt_p1++;
					}

				}
				logprob = max_p0 - max_p1;
				// llr = log (p(0) / p(1)) 
				// float logprob = log(prob0 / prob1);
				if (logprob > 10) logprob = 10;
				if (logprob < -10) logprob = -10;
				llr_bit[i * ModType_half + i_bit] = logprob;
			}
		}
		break;
	}
	delete[] u_sym_reverse;
	delete[] sym_prob;
}






void llr_mmse_bit_exact(int Nt, float* mu_mmse,  float* u_sym_mmse, float* llr_bit, float sigma_sq, int ModType) {
	float prob1 = 0;
	float prob0 = 0;
	float norm_factor = 0;
	float eps = 0.001; // sigma_mmse_min

	int Nt2 = 2 * Nt;
	int ModType_half = ModType / 2;
	int Q = pow(2, ModType_half);

	float* u_sym_reverse = new float[Nt2];
	float* sym_prob = new float[Q];   // prob of every half symbol in the alphabet


	for (int i = 0; i < Nt; i++)
	{
		u_sym_reverse[2 * i] = u_sym_mmse[i]; // real part
		u_sym_reverse[2 * i + 1] = u_sym_mmse[i + Nt]; // imag part
	}
	switch (ModType)
	{
	case 4:
		for (int i = 0; i < Nt2; i++)
		{
			float sigma_sq_mmse = mu_mmse[i] - mu_mmse[i] * mu_mmse[i];
			if (sigma_sq_mmse < eps) sigma_sq_mmse = eps;
			// prob of every symbol in the alphabet tablee p(u) = exp[-(x-u_mmse)^2/(2*sigma_sq)]
			//norm_factor = 0;
			float sym_est =  u_sym_reverse[i];
			for (int i_sym = 0; i_sym < Q; i_sym++)
			{
				float sym_val = table_16qam[i_sym];
				float dist = abs(mu_mmse[i]*sym_val - sym_est);
				//float dist2_nv = -pow(dist, 2) / (2.0 * sigma_sq);
				float dist2_nv = -pow(dist, 2) / (sigma_sq_mmse);
				//sym_prob[i_sym] = exp(-dist2_nv);
				sym_prob[i_sym] = dist2_nv;
				norm_factor += sym_prob[i_sym];
			}
			 //normalization
			// for (int i_sym = 0; i_sym < Q; i_sym++) { sym_prob[i_sym] /= norm_factor;  /*cout << sym_prob[i_sym] << "  ";*/ }
			 //p1 and p0
			for (int i_bit = 0; i_bit < ModType_half; i_bit++)
			{
				prob0 = 0;
				prob1 = 0;
				/*for (int i_sym = 0; i_sym < Q; i_sym++)
				{
					if (s_bits_r4[i_sym][i_bit] == 0) { prob0 += sym_prob[i_sym]; }
					if (s_bits_r4[i_sym][i_bit] == 1) { prob1 += sym_prob[i_sym]; }
				}*/
				float logprob = 0;
				if (i_bit == 0)
				{
					// cout << endl;
					logprob = std::max(sym_prob[0], sym_prob[1]) - std::max(sym_prob[2], sym_prob[3]);
					// cout << logprob << endl;
				}
				else if (i_bit == 1)
				{
					logprob = std::max(sym_prob[0], sym_prob[2]) - std::max(sym_prob[1], sym_prob[3]);
					//cout << logprob << endl;
				}
				// llr = log (p(0) / p(1)) 
				// float logprob = log(prob0 / prob1);
				// 16 QAM 32 x 256 : 0.8
				// 16 QAM 8 x 64 : 10 // 10
				// 1D: 10 or -10
				if (logprob > 100) logprob = 100;
				if (logprob < -100) logprob = -100;
				llr_bit[i * ModType_half + i_bit] = logprob;
			}
			// test: boundary decision
			/*llr_bit[i * ModType_half + 0] = -sym_est;
			llr_bit[i * ModType_half + 1] = -abs(sym_est) + 0.632455;*/
		}

		break;
	case 6:
		for (int i = 0; i < Nt2; i++)
		{
			// prob of every symbol in the alphabet tablee p(u) = exp[-(x-u_mmse)^2/(2*sigma_sq)]
			norm_factor = 0;
			float sym_est = u_sym_reverse[i];
			for (int i_sym = 0; i_sym < Q; i_sym++)
			{
				float sym_val = table_64qam[i_sym];
				float dist = abs(sym_val - sym_est);
				float dist2_nv = pow(dist, 2) / (2.0 * sigma_sq);
				sym_prob[i_sym] = exp(-dist2_nv);
				norm_factor += sym_prob[i_sym];
			}
			// normalization
			for (int i_sym = 0; i_sym < Q; i_sym++) { sym_prob[i_sym] /= norm_factor; }
			// p1 and p0
			for (int i_bit = 0; i_bit < ModType_half; i_bit++)
			{
				prob0 = 0;
				prob1 = 0;
				for (int i_sym = 0; i_sym < Q; i_sym++)
				{
					if (s_bits_r6[i_sym][i_bit] == 0) { prob0 += sym_prob[i_sym]; }
					if (s_bits_r6[i_sym][i_bit] == 1) { prob1 += sym_prob[i_sym]; }
				}
				// llr = log (p(0) / p(1)) 
				float logprob = log(prob0 / prob1);
				if (logprob > 10000) logprob = 10000;
				if (logprob < -10000) logprob = -10000;
				llr_bit[i * ModType_half + i_bit] = logprob;
			}
		}
		break;
	case 8:
		for (int i = 0; i < Nt2; i++)
		{
			// prob of every symbol in the alphabet tablee p(u) = exp[-(x-u_mmse)^2/(2*sigma_sq)]
			norm_factor = 0;
			float sym_est = u_sym_reverse[i];
			for (int i_sym = 0; i_sym < Q; i_sym++)
			{
				float sym_val = table_256qam[i_sym];
				float dist = abs(sym_val - sym_est);
				float dist2_nv = -pow(dist, 2) / (2.0 * sigma_sq);
				//sym_prob[i_sym] = exp(-dist2_nv);
				sym_prob[i_sym] = dist2_nv;
				norm_factor += sym_prob[i_sym];
			}
			// normalization
			// for (int i_sym = 0; i_sym < Q; i_sym++) { sym_prob[i_sym] /= norm_factor;  cout << sym_prob[i_sym] << "  "; }
			// p1 and p0
			for (int i_bit = 0; i_bit < ModType_half; i_bit++)
			{
				prob0 = 0;
				prob1 = 0;
				/*for (int i_sym = 0; i_sym < Q; i_sym++)
				{
					if (s_bits_r4[i_sym][i_bit] == 0) { prob0 += sym_prob[i_sym]; }
					if (s_bits_r4[i_sym][i_bit] == 1) { prob1 += sym_prob[i_sym]; }
				}*/
				float logprob = 0;
				float max_p0 = 0;
				float max_p1 = 0;
				int cnt_p0 = 0;
				int cnt_p1 = 0;
				for (int i_sym_prob = 0; i_sym_prob < Q; i_sym_prob++)
				{
					if (s_bits_r8[i_sym_prob][i_bit] == 0)
					{
						if (cnt_p0 == 0) { max_p0 = sym_prob[i_sym_prob]; }
						else
						{
							if (max_p0 < sym_prob[i_sym_prob])
								max_p0 = sym_prob[i_sym_prob];
						}
						cnt_p0++;
					}
					else
					{
						if (cnt_p1 == 0) { max_p1 = sym_prob[i_sym_prob]; }
						else
						{
							if (max_p1 < sym_prob[i_sym_prob])
								max_p1 = sym_prob[i_sym_prob];
						}
						cnt_p1++;
					}

				}
				logprob = max_p0 - max_p1;
				// llr = log (p(0) / p(1)) 
				// float logprob = log(prob0 / prob1);
				if (logprob > 10) logprob = 10;
				if (logprob < -10) logprob = -10;
				llr_bit[i * ModType_half + i_bit] = logprob;
			}
		}
		break;
	}
	delete[] u_sym_reverse;
	delete[] sym_prob;
}


// bit llr output for kbest detection
void llr_kbest_bit(int Nt, float* u_sym_kbest, float* llr_bit, float sigma_sq, int ModType, int k_kbest, float** paths_kbest, float* ED) {
	float prob1 = 0;
	float prob0 = 0;
	float norm_factor = 0;

	int Nt2 = 2 * Nt;
	int ModType_half = ModType / 2;
	int Q = pow(2, ModType_half);
	int* idx_ed = new int[k_kbest];
	int** bit_array = new int* [k_kbest];
	float* llr_bit_rearranged = new float[Nt2 * ModType_half];
	for (int i = 0; i < k_kbest; i++) { bit_array[i] = new int[Nt2 * ModType_half]; }
	for (int i = 0; i < k_kbest; i++)
	{
		int i_sym_in_set = 0;
		for (int j = 0; j < Nt2; j++)
		{
			i_sym_in_set = 0;
			for (int i_sym = 0; i_sym < Q; i_sym++)
			{
				if (paths_kbest[i][j] == table_16qam[i_sym])
				{
					i_sym_in_set = i_sym;
					break;
				}
			}
			for (int i_bit = 0; i_bit < ModType_half; i_bit++)
			{
				bit_array[i][j * ModType_half + i_bit] = s_bits_r4[i_sym_in_set][i_bit];
			}
		}
	}
	// test: bit array
	/*cout << "bit array" << endl;
	for (int i = 0; i < Nt2 * ModType_half; i++)
	{
		if (!(i % ModType_half)) { cout << endl; }
		cout << bit_array[0][i];
	}
	cout << endl;
	cout << "bit array" << endl;*/
	// test: bit array
	// Sort the EDs.
	my_sort_ed(ED, idx_ed, k_kbest);
	int best_idx = idx_ed[0];
	float dmin = ED[0];
	float dmin_diff = 0;
	float llr = 0;
	// Calculate the bitwise LLRs
	for (int i_bit = 0; i_bit < Nt2 * ModType_half; i_bit++)
	{
		int bit_dec = bit_array[best_idx][i_bit];
		int firstdiff = 0;
		for (int i_path = 0; i_path < k_kbest; i_path++)
		{
			if (bit_array[idx_ed[i_path]][i_bit] != bit_dec) { dmin_diff = ED[i_path]; firstdiff = 1; break; }
		}
		if (firstdiff)
		{
			// llr_bit_rearranged[i_bit] = (1 - 2 * bit_dec) * (dmin_diff - dmin);
			llr_bit_rearranged[i_bit] = (1 - 2 * bit_dec) * (dmin_diff - dmin) / (2 * sigma_sq);
		}
		else
		{
			// llr_bit_rearranged[i_bit] = (1 - 2 * bit_dec) * (ED[k_kbest - 1] - dmin);
			llr_bit_rearranged[i_bit] = (1 - 2 * bit_dec) * (ED[k_kbest - 1] - dmin) / (2 * sigma_sq);
		}
	}
	float* u_sym_reverse = new float[Nt2];
	float* sym_prob = new float[Q];   // prob of every half symbol in the alphabet
	for (int i = 0; i < Nt; i++)
	{
		u_sym_reverse[2 * i] = u_sym_kbest[i]; // real part
		u_sym_reverse[2 * i + 1] = u_sym_kbest[i + Nt]; // imag part
	}

	// Rearrange the bitwise llrs
	for (int i_sym = 0; i_sym < Nt; i_sym++)
	{
		for (int i_bit = 0; i_bit < ModType_half; i_bit++)
		{
			// real part
			llr_bit[i_sym * ModType + i_bit] = llr_bit_rearranged[i_sym * ModType_half + i_bit];
			// imaginary part
			llr_bit[i_sym * ModType + ModType_half + i_bit] = llr_bit_rearranged[Nt * ModType_half + i_sym * ModType_half + i_bit];
		}
	}
	// transform symbol sequence to bit sequence
	// test:llr_bit
	/*cout << "LLR_BIT" << endl;
	for (int i = 0; i < Nt2 * ModType_half; i++)
	{
		cout << "     " << llr_bit[i];
	}
	cout << "LLR_BIT" << endl;*/
	// test:llr_bit
	for (int i = 0; i < Nt2 * ModType_half; i++)
	{

		// test: for beta
		if (llr_bit[i] > 5) llr_bit[i] = 5;
		if (llr_bit[i] < -5) llr_bit[i] = -5;
		/*llr_bit[i] /= 30;
		if (llr_bit[i] > 3.5) llr_bit[i] = 3.5;
		if (llr_bit[i] < -3.5) llr_bit[i] = -3.5;*/
	}
	// test: bitwise LLRs
	/*cout << "Hard Decision from LLR" << endl;
	for (int i = 0; i < Nt * ModType; i++) { cout << (llr_bit_rearranged[i] < 0); }
	cout << endl;
	cout << "Hard Decision from LLR" << endl;*/
	// test: bitwise LLRs
	for (int i = 0; i < k_kbest; i++) { delete[] bit_array[i]; }
	delete[] idx_ed;
	delete[] bit_array;
	delete[] u_sym_reverse;
	delete[] sym_prob;
	delete[] llr_bit_rearranged;
}


// sym output to soft bitwise llr output
// u_sym_mmse: ep_extr_mean
// sigma_sq: ep_extr_var
// only 16QAM is hereby supported
void llr_ep_bit(int Nt, float* ep_extr_mean, float* ep_extr_var, float* llr_bit, int ModType) {
	float prob1 = 0;
	float prob0 = 0;
	float norm_factor = 0;

	int Nt2 = 2 * Nt;
	int ModType_half = ModType / 2;
	int Q = pow(2, ModType_half);

	float* u_sym_reverse = new float[Nt2];
	float* var_sym_reverse = new float[Nt2];
	float* sym_prob = new float[Q];   // prob of every half symbol in the alphabet


	for (int i = 0; i < Nt; i++)
	{
		u_sym_reverse[2 * i] = ep_extr_mean[i]; // real part
		u_sym_reverse[2 * i + 1] = ep_extr_mean[i + Nt]; // imag part
		var_sym_reverse[2 * i] = ep_extr_var[i]; // real part
		var_sym_reverse[2 * i + 1] = ep_extr_var[i + Nt]; // imag part
	}
	switch (ModType)
	{
	case 4:
		for (int i = 0; i < Nt2; i++)
		{
			// prob of every symbol in the alphabet tablee p(u) = exp[-(x-u_mmse)^2/(2*sigma_sq)]
			norm_factor = 0;
			float sym_est = u_sym_reverse[i];
			float sigma_sq = var_sym_reverse[i];
			for (int i_sym = 0; i_sym < Q; i_sym++)
			{
				float sym_val = table_16qam[i_sym];
				float dist = abs(sym_val - sym_est);
				float dist2_nv = -pow(dist, 2) / (2.0 * sigma_sq);
				//sym_prob[i_sym] = exp(-dist2_nv);
				sym_prob[i_sym] = dist2_nv;
				norm_factor += sym_prob[i_sym];
			}
			 //normalization
			//for (int i_sym = 0; i_sym < Q; i_sym++) { sym_prob[i_sym] /= norm_factor;  /*cout << sym_prob[i_sym] << "  ";*/ }
			 //p1 and p0
			int max_index = 0;
			//find closest symbol
			for (int i_sym_close = 1; i_sym_close < 4; i_sym_close++)
			{
				if (sym_prob[i_sym_close] > sym_prob[max_index]) max_index = i_sym_close;
			}
			for (int i_bit = 0; i_bit < ModType_half; i_bit++)
			{
				prob0 = 0;
				prob1 = 0;
				/*for (int i_sym = 0; i_sym < Q; i_sym++)
				{
					if (s_bits_r4[i_sym][i_bit] == 0) { prob0 += sym_prob[i_sym]; }
					if (s_bits_r4[i_sym][i_bit] == 1) { prob1 += sym_prob[i_sym]; }
				}*/
				float logprob = 0;
				if (i_bit == 0)
				{
					// cout << endl;
					logprob = std::max(sym_prob[0], sym_prob[1]) - std::max(sym_prob[2], sym_prob[3]);
					// cout << logprob << endl;
					//if (max_index == 0 || max_index == 1) { logprob = 10; }
					//else logprob = -10;
				}
				else if (i_bit == 1)
				{
					logprob = std::max(sym_prob[0], sym_prob[2]) - std::max(sym_prob[1], sym_prob[3]);
					//cout << logprob << endl;
					//if (max_index == 0 || max_index == 2) { logprob = 10; }
					//else logprob = -10;
				}
				// llr = log (p(0) / p(1)) 
				// float logprob = log(prob0 / prob1);
				// 16 QAM 32 x 256 : 0.8
				// 16 QAM 8 x 64 : 10 // 10
				// 1D: 10 or -10
				if (logprob > 10) logprob = 10;
				if (logprob < -10) logprob = -10;
				llr_bit[i * ModType_half + i_bit] = logprob;
			}
			// test: boundary decision
			/*llr_bit[i * ModType_half + 0] = -sym_est;
			llr_bit[i * ModType_half + 1] = -abs(sym_est) + 0.632455;*/
		}

		break;
	}
	delete[] u_sym_reverse;
	delete[] sym_prob;
	delete[] var_sym_reverse;
}


void my_sort_ed(float* arr, int* idx_sort, int N) {
	for (int i = 0; i < N; ++i) {
		// cout << N << endl;
		idx_sort[i] = i;
	}
	std::sort(idx_sort, idx_sort + N, [&](int a, int b) {
		return arr[a] < arr[b];
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

