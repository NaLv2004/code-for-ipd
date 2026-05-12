#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <stdio.h>
#include <stdlib.h> 
#include <time.h> 
#include <iostream>
#include <cmath>
#include <bitset>
#include <algorithm>
#include <vector>
#include <xmmintrin.h>
#include <immintrin.h>
#include <smmintrin.h>
using namespace std;


#define SIGN(n) (n==0?0:(n/abs(n)))
const int con5G128[128] = { 0,1,2,4,8,16,32,3,5,64,9,6,17,10,18,12,33,65,20,34,24,36,7,66,11,40,68,19,13,48,14,72,21,35,26,80,37,25,22,38,96,67,41,28,69,42,49,74,70,44,81,50,73,15,52,23,76,82,56,27,97,39,84,29,43,98,88,30,71,45,100,51,46,75,104,53,77,54,83,57,112,78,85,58,99,86,60,89,101,31,90,102,105,92,47,106,55,113,79,108,59,114,87,116,61,91,120,62,103,93,107,94,109,115,110,117,118,121,122,63,124,95,111,119,123,125,126,127 };
const int con5G1024[1024] = { 0,1,2,4,8,16,32,3,5,64,9,6,17,10,18,128,12,33,65,20,256,34,24,36,7,129,66,512,11,40,68,130,19,13,48,14,72,257,21,132,35,258,26,513,80,37,25,22,136,260,264,38,514,96,67,41,144,28,69,42,516,49,74,272,160,520,288,528,192,544,70,44,131,81,50,73,15,320,133,52,23,134,384,76,137,82,56,27,97,39,259,84,138,145,261,29,43,98,515,88,140,30,146,71,262,265,161,576,45,100,640,51,148,46,75,266,273,517,104,162,53,193,152,77,164,768,268,274,518,54,83,57,521,112,135,78,289,194,85,276,522,58,168,139,99,86,60,280,89,290,529,524,196,141,101,147,176,142,530,321,31,200,90,545,292,322,532,263,149,102,105,304,296,163,92,47,267,385,546,324,208,386,150,153,165,106,55,328,536,577,548,113,154,79,269,108,578,224,166,519,552,195,270,641,523,275,580,291,59,169,560,114,277,156,87,197,116,170,61,531,525,642,281,278,526,177,293,388,91,584,769,198,172,120,201,336,62,282,143,103,178,294,93,644,202,592,323,392,297,770,107,180,151,209,284,648,94,204,298,400,608,352,325,533,155,210,305,547,300,109,184,534,537,115,167,225,326,306,772,157,656,329,110,117,212,171,776,330,226,549,538,387,308,216,416,271,279,158,337,550,672,118,332,579,540,389,173,121,553,199,784,179,228,338,312,704,390,174,554,581,393,283,122,448,353,561,203,63,340,394,527,582,556,181,295,285,232,124,205,182,643,562,286,585,299,354,211,401,185,396,344,586,645,593,535,240,206,95,327,564,800,402,356,307,301,417,213,568,832,588,186,646,404,227,896,594,418,302,649,771,360,539,111,331,214,309,188,449,217,408,609,596,551,650,229,159,420,310,541,773,610,657,333,119,600,339,218,368,652,230,391,313,450,542,334,233,555,774,175,123,658,612,341,777,220,314,424,395,673,583,355,287,183,234,125,557,660,616,342,316,241,778,563,345,452,397,403,207,674,558,785,432,357,187,236,664,624,587,780,705,126,242,565,398,346,456,358,405,303,569,244,595,189,566,676,361,706,589,215,786,647,348,419,406,464,680,801,362,590,409,570,788,597,572,219,311,708,598,601,651,421,792,802,611,602,410,231,688,653,248,369,190,364,654,659,335,480,315,221,370,613,422,425,451,614,543,235,412,343,372,775,317,222,426,453,237,559,833,804,712,834,661,808,779,617,604,433,720,816,836,347,897,243,662,454,318,675,618,898,781,376,428,665,736,567,840,625,238,359,457,399,787,591,678,434,677,349,245,458,666,620,363,127,191,782,407,436,626,571,465,681,246,707,350,599,668,790,460,249,682,573,411,803,789,709,365,440,628,689,374,423,466,793,250,371,481,574,413,603,366,468,655,900,805,615,684,710,429,794,252,373,605,848,690,713,632,482,806,427,904,414,223,663,692,835,619,472,455,796,809,714,721,837,716,864,810,606,912,722,696,377,435,817,319,621,812,484,430,838,667,488,239,378,459,622,627,437,380,818,461,496,669,679,724,841,629,351,467,438,737,251,462,442,441,469,247,683,842,738,899,670,783,849,820,728,928,791,367,901,630,685,844,633,711,253,691,824,902,686,740,850,375,444,470,483,415,485,905,795,473,634,744,852,960,865,693,797,906,715,807,474,636,694,254,717,575,913,798,811,379,697,431,607,489,866,723,486,908,718,813,476,856,839,725,698,914,752,868,819,814,439,929,490,623,671,739,916,463,843,381,497,930,821,726,961,872,492,631,729,700,443,741,845,920,382,822,851,730,498,880,742,445,471,635,932,687,903,825,500,846,745,826,732,446,962,936,475,853,867,637,907,487,695,746,828,753,854,857,504,799,255,964,909,719,477,915,638,748,944,869,491,699,754,858,478,968,383,910,815,976,870,917,727,493,873,701,931,756,860,499,731,823,922,874,918,502,933,743,760,881,494,702,921,501,876,847,992,447,733,827,934,882,937,963,747,505,855,924,734,829,965,938,884,506,749,945,966,755,859,940,830,911,871,639,888,479,946,750,969,508,861,757,970,919,875,862,758,948,977,923,972,761,877,952,495,703,935,978,883,762,503,925,878,735,993,885,939,994,980,926,764,941,967,886,831,947,507,889,984,751,942,996,971,890,509,949,973,1000,892,950,863,759,1008,510,979,953,763,974,954,879,981,982,927,995,765,956,887,985,997,986,943,891,998,766,511,988,1001,951,1002,893,975,894,1009,955,1004,1010,957,983,958,987,1012,999,1016,767,989,1003,990,1005,959,1011,1013,895,1006,1014,1017,1018,991,1020,1007,1015,1019,1021,1022,1023 };


float MyRecursiveFun(int i, int NN, float p);//计算I(W)
void sub_sort(float* FilterArray1, int* Posit, int len);//从大到小排列并输出对应的下标
double phi(double t)
{
	if (t < 0.867861)
		return std::exp(0.0564 * t * t - 0.48560 * t);
	else // if(t >= phi_pivot)
		return std::exp(-0.4527 * std::pow(t, 0.8600) + 0.0218);
}

double phi_inv(double t)
{
	if (t > 0.6845772418)
		return 4.304964539 * (1 - sqrt(1 + 0.9567131408 * std::log(t)));
	else
		return std::pow(-2.208968 * std::log(t) + 0.0482, 1.1628);
}
void GA_codeconstruction(int CodeLength, float sigma, vector<int>& best_channel)
{
	int m = (int)(log(CodeLength) / log(2));
	std::vector<float> z(CodeLength, 0);
	const float alpha = -0.4527;
	const float beta = 0.0218;
	const float gamma = 0.8600;
	const float bisection_max = std::numeric_limits<float>::max();
	const float epsilon = 0.00000000001;

	for (unsigned i = 0; i < CodeLength; i++)
		best_channel[i] = i;

	//for (auto i = 0; i < std::exp2(m); i++)
	for (auto i = 0; i < pow(2.0, m); i++)
		z[i] = 2.0 / std::pow((float)sigma, 2.0);

	for (auto l = 1; l <= m; l++)
	{
		//auto o1 = (int)std::exp2(m - l + 1);
		//auto o2 = (int)std::exp2(m - l);
		auto o1 = (int)pow(2.0, m - l + 1);
		auto o2 = (int)pow(2.0, m - l);
		//for (auto t = 0; t < (int)std::exp2(l - 1); t++)
		for (auto t = 0; t < (int)pow(2.0, l - 1); t++)
		{
			float T = z[t * o1];

			z[t * o1] = phi_inv(1.0 - std::pow(1.0 - phi(T), 2.0));
			if (z[t * o1] == HUGE_VAL)
				z[t * o1] = T + M_LN2 / (alpha * gamma);

			z[t * o1 + o2] = 2.0 * T;
		}
	}
	std::sort(best_channel.begin(), best_channel.end(), [&](int i1, int i2) { return z[i1] > z[i2]; });
}

void B5G_codeconstruction(int N, int K, vector<int>& best_channel)
{
	int* setcon = new int[N];
	int cnt = 0;
	for (int i = 0; i < N; i++) setcon[i] = 0;
	if (N == 128)
		for (int i = 0; i < K; i++) setcon[con5G128[N - 1 - i]] = 1;
	else if (N == 1024)
		for (int i = 0; i < K; i++) setcon[con5G1024[N - 1 - i]] = 1;
	else {
		int end = 0;
		int i = 0;
		while (end < K) {
			if (con5G1024[1024 - 1 - i] < N) {
				setcon[con5G1024[1024 - 1 - i]] = 1;
				end++;
			}
			i++;
		}
	}
	for (int j = 0; j < N; j++)
		if (setcon[j]) best_channel[cnt++] = j;
	for (int j = 0; j < N; j++)
		if (!setcon[j]) best_channel[cnt++] = j;
}
float MyRecursiveFun(int i, int NN, float p)//p is epsilon of BSC
{
	/* i     the index of channel
	NN    the number of channel(Codelength)
	p     epsilon of BSC channel */
	float y = 0;
	if (i == 1 && NN == 1) y = 1 - p;
	else {
		if (i % 2 == 1) {
			float temp = MyRecursiveFun((i + 1) / 2, NN / 2, p);
			y = pow(temp, 2);
		}
		else {
			float temp = MyRecursiveFun(i / 2, NN / 2, p);
			y = 2 * temp - pow(temp, 2);
		}
	}
	return y;
}

void sub_sort(float* FilterArray1, int* Posit, int len)
{
	int i, j, count = 0;
	for (i = 0; i < len; i++)
		Posit[i] = 255;
	for (i = 0; i < len; i++)
	{
		count = 0;
		for (j = 0; j < len; j++)
			if (i != j)
				if (FilterArray1[i] < FilterArray1[j])
					count++;
		//如果此数大于所有的数，则count==14，Posit[14]记录下该数的下标
		while (Posit[count] != 255)count++;
		//此地已有数据了，则移动到下一个数据区，即count值和前面值一致，说明此数据和占据些地的数据值相等。
		Posit[count] = i;
	}
}
/*int* bitreorder(int* min, int len)
{
	int *a = new int[len];
	int *b = new int[len];
	int index = 0;
	for (int i = 0; i < len; i++) {
		a[i] = index;
		index++;
		bitset<n> bin(a[i]);
		bitset<n> bout;
		for (int j = 0; j < n; j++)
		{
			bout[n - j - 1] = bin[j];
		}
		b[i] = bout.to_ulong();
	}
	int* mout = new int[len];
	for (int i = 0; i < len; i++) {
		mout[i] = min[b[i]];
	}
	return mout;
}*/

void RM_codeconstruction(int N, int** GG, vector<int>& best_channel)
{
	int* discount = new int[N];
	for (int i = 0; i < N; i++)
	{
		discount[i] = 0;
		for (int j = 0; j < N; j++)
			if (GG[i][j])	discount[i]++;
	}
	for (int i = 0; i < N; i++)
	{
		int max = -1, index;
		for (int j = 0; j < N; j++)
			if (discount[j] > max)
			{
				max = discount[j];
				index = j;
			}
		discount[index] = -2;
		best_channel[i] = index;
	}
	delete[] discount;
	return;
}

void beta_codeconstruction(int N, vector<int>& best_channel)
{
	float* beta = new float[N];
	int tmp;
	for (int i = 0; i < N; i++)
	{
		best_channel[i] = i;
		int j = i;
		float k = 0;
		while (j > 0)
		{
			beta[i] += float(j % 2) * pow(2, k / 4);
			j = j / 2;
			k++;
		}
	}
	for (int i = 0; i < N - 1; i++)
		for (int j = i + 1; j < N; j++)
			if (beta[best_channel[i]] < beta[best_channel[j]])
			{
				tmp = best_channel[i];
				best_channel[i] = best_channel[j];
				best_channel[j] = tmp;
			}
	//for (int i = 0; i < N - 1; i++)
	//	if (beta[best_channel[i]] == beta[best_channel[i + 1]]) cout << "! ";
	//getchar();
	delete[] beta;
}


void PolarEncode_xor(int* uout, int* uin, int len)
{
	if (len == 1) return;
	else if (len == 2)
	{
		uout[0] = uin[0] + uin[1];
		uout[1] = uin[1];
	}
	else if (len == 4)
	{
		*(uout + 1) = *(uin + 1) ^ *(uin + 3);
		*(uout + 2) = *(uin + 2) ^ *(uin + 3);
		*(uout + 3) = *(uin + 3);
		*(uout) = *(uin) ^ *(uin + 1) ^ *(uout + 2);
	}
	else
	{
		for (int i = 0; i < len; i += 2)
		{
			uout[i] = uin[i] ^ uin[i + 1];
			uout[i + 1] = uin[i + 1];

		}
		for (int i = 0; i < len; i += 4)
		{
			for (int k = 0; k < 2; k++)
			{
				uout[i + k] ^= uout[i + k + 2];
			}
		}
		int temp = 8;
		/*int temp = 4;
		for (int i = 0; i < len; i += temp)
		{
			*(uout + i + 1) = *(uin + i + 1) ^ *(uin + i + 3);
			*(uout + i + 2) = *(uin + i + 2) ^ *(uin + i + 3);
			*(uout + i + 3) = *(uin + i + 3);
			*(uout + i) = *(uin + i) ^ *(uin + i + 1) ^ *(uout + i + 2);
		}
		temp <<= 1;*/

		__m128i* p = (__m128i*)uout;
		for (int i = 0; i < len; i += temp) {
			__m128i a_128i = _mm_loadu_si128(p);
			__m128i b_128i = _mm_loadu_si128(p + 1);
			_mm_storeu_si128(p, _mm_xor_si128(a_128i, b_128i));
			p += 2;
		}
		temp <<= 1;
		int cnt = 1;
		__m256i* q = (__m256i*)uout;
		while (len >= temp) {
			for (int i = 0; i < len; i += temp) {
				for (int j = 0; j < temp / 16; j++) {
					__m256i a_256i = _mm256_loadu_si256(q + j);
					__m256i b_256i = _mm256_loadu_si256(q + j + cnt);
					_mm256_storeu_si256(q + j, _mm256_xor_si256(a_256i, b_256i));
				}
				q += (cnt << 1);
			}
			temp <<= 1;
			cnt <<= 1;
			q = (__m256i*)uout;
		}
	}
}


// Encode a portion of the code
void PolarEncodePartial(int** G, int* u_in, int* u_out, int N, int u_in_len)
{
	for (int i = 0; i < u_in_len; i++)
	{
		int row_G_min = N - 1 - i;
		int col_G = N - 1 - i;
		//cout << col_G << endl;
		u_out[col_G] = 0;
		cout << row_G_min << endl;

		for (int row_G = N - 1; row_G >= row_G_min; row_G--)
		{

			u_out[col_G] ^= (u_in[row_G] & G[row_G][col_G]);
		}
	}
}


double sampleNormal() {
	double p = ((double)rand() / (RAND_MAX)) * 2 - 1;
	double v = ((double)rand() / (RAND_MAX)) * 2 - 1;
	double r = p * p + v * v;
	if (r == 0 || r > 1) return sampleNormal();
	double c = sqrt(-2 * log(r) / r);
	return p * c;
}

class PolarCodewords {
public:
	const int MAX_BITS = 32;
	int N;
	int K;
	int* A_Ac;
	int* A;
	int* Ac;
	int** test_patterns;
	int* test_pattern_indexes;

	PolarCodewords(int N, int K, int* A_Ac, int* A, int* Ac) {
		this->N = N;
		this->K = K;
		this->A_Ac = new int[N];
		this->A = new int[K];
		this->Ac = new int[N - K];
		for (int i = 0; i < N; i++) {
			this->A_Ac[i] = A_Ac[i];
		}
		for (int i = 0; i < K; i++) {
			this->A[i] = A[i];
		}
		for (int i = 0; i < (N - K); i++) {
			this->Ac[i] = Ac[i];
		}
	}

	void generateTestPatterns()
	{
		int n_test_patterns = pow(2, this->K);
		int n_test_pattern_idx = pow(2, this->N);
		this->test_patterns = new int* [n_test_patterns];
		int* u_to_encode = new int[this->N];
		for (int i = 0; i < N; i++) { u_to_encode[i] = 0; }
		for (int i = 0; i < n_test_patterns; i++) {
			test_patterns[i] = new int[this->N];
			for (int j = 0; j < this->N; j++) { this->test_patterns[i][j] = 0; }
		}
		this->test_pattern_indexes = new int[n_test_pattern_idx];
		for (int i = 0; i < n_test_pattern_idx; i++) {
			this->test_pattern_indexes[i] = 0;
		}
		// generate all possible polar codewords
		for (int i = 0; i < n_test_patterns; i++)
		{
			int* u_encoded = new int[this->N];
			for (int j = 0; j < this->N; j++) { u_encoded[j] = 0; }
			std::bitset<32> bits(i);
			for (int j = 0; j < this->K; j++) { u_to_encode[this->A[j]] = bits[j]; }
			PolarEncode_xor(u_encoded, u_to_encode, this->N);
			for (int j = 0; j < this->N; j++) { this->test_patterns[i][j] = u_encoded[j]; }
			int dec = 0;
			int weight = 1;
			for (int j = 0; j < this->N; j++) { dec += u_encoded[j] * weight; weight *= 2; }
			this->test_pattern_indexes[dec] = i;
		}
	}

	int getTestPattern(int test_pattern_id, int i_bit) {
		return this->test_patterns[test_pattern_indexes[test_pattern_id]][i_bit];
	}
};
