#include "Polar_Encoder.h"
#include "crc.h"
#include "polar_Function.h"
#include "SP_node.h"
#include <xmmintrin.h>
#include <immintrin.h>
#include <smmintrin.h>
#include <vector>
#include <fstream>
#include <ctime>
#include <functional>
#include <algorithm>
#include <queue>

#define hard(n) (n>=0?0:1)

ofstream fout("temp_output.txt");
ofstream wout("index_output.txt");
int* list_sub_sort2(float* FilterArray1, int len);

struct Comp {
	Comp(const vector<float> & v) : _v(v) {}
	bool operator ()(int a, int b) { return _v[a] > _v[b]; }
	const vector<float> & _v;
};
void decode(float* LLR, int* sum, int j, int* A_Ac, int a, int N, int K,int* Count)
{
	/*	LLR   N-length LLR as input
	sum   N-length beta value as output
	j     the number of stage, equaling to log2N
	A_Ac  the array which label the channel condition, frozen is 0, information is 1
	a     a variable which label the position the start of each array (alpha and beta)
	N     Codelength
	K     Inoformation bit length  */
	if (j == 0) {     // leaf node
		if (A_Ac[*Count] == 0) {
			*(sum + *Count) = 0;
		}
		else {
			*(sum + *Count) = hard(*(LLR + (N << 1) - 2));
		}
		*Count = *Count+1;
		return;
	}

	int node = 1 << j;
	int length = node >> 1;
	int start = (N << 1) - node;
	int last_start = (N << 1) - (node << 1);
	int type = 4;
	//if (j < 6) { //An adaptive operation, if j < 6, use Fast-SSC
	//	type = judge_rate(A_Ac + Count, node);
	//}
	//else type = 4;
	//type = 4;
	switch (type)
	{
	case 0:  //rate-0 node
	{
		memset(sum + *Count, 0, sizeof(int) * node);
		Count += node;
		break;
	}
	case 1:  //rate-1 node
	{
		hard_SIMD(LLR + last_start, sum + *Count, node);
		Count += node;
		break;
	}

	case 2:  //repetition node
	{
		float temp = cal_sum(LLR + last_start, node);
		if (temp >= 0) {
			memset(sum + *Count, 0, sizeof(int) * node);
		}
		else {
			for (int i = 0; i < node; i++)
			{
				*(sum + *Count + i) = 1;
			}
		}
		Count += node;
		break;
	}
	case 3:  //SPC node
	{
		int parity = 0;
		float temp = abs(*(LLR + last_start));
		int index = 0;
		for (int i = 0; i < node; i++)
		{
			if (*(LLR + last_start + i) >= 0)
			{
				*(sum + *Count + i) = 0;
				if (temp > *(LLR + last_start + i)) {
					temp = *(LLR + last_start + i); index = i;
				}
			}
			else
			{
				*(sum + *Count + i) = 1;
				parity ^= 1;
				if (temp > abs(*(LLR + last_start + i))) {
					temp = abs(*(LLR + last_start + i)); index = i;
				}
			}
		}
		*(sum + *Count + index) ^= parity;
		Count += node;
		break;
	}
	case 4:  //conventional SC operation
	{
		f_function(LLR + last_start, length);
		decode(LLR, sum, j - 1, A_Ac, a, N, K, Count);
		a = *Count >> (j - 1);

		g_function(LLR + last_start, sum + ((a - 1)*length), length);
		decode(LLR, sum, j - 1, A_Ac, a, N, K, Count);
		a = *Count >> (j - 1);
		int initial_posit = (a - 2) << (j - 1);
		if (length < 4) {
			for (int i = 0; i < length; i++)
			{
				*(sum + initial_posit + i) ^= *(sum + initial_posit + length + i);
			}
		}
		else {
			combine(sum + initial_posit, sum + initial_posit + length, length);
		}
		break;
	}
	default:cout << "error" << endl;
	}
	//printf("%d ",Count);
	//getchar();
	return;
}


int* list_sub_sort2(float* FilterArray1, int len)
{
	int i, j, cnt = 0;
	int* Posit = new int[len];
	for (i = 0; i < len; i++)
		Posit[i] = 255;
	for (i = 0; i < len; i++)
	{
		cnt = 0;
		for (j = 0; j < len; j++)
			if (i != j)
				if (FilterArray1[i]<FilterArray1[j])
					cnt++;
		//如果此数大于所有的数，则count==14，Posit[14]记录下该数的下标
		while (Posit[cnt] != 255)cnt++;
		//此地已有数据了，则移动到下一个数据区，即count值和前面值一致，说明此数据和占据些地的数据值相等。
		Posit[cnt] = i;
	}
	return Posit;
}
void sort_list(int* p, int* q, float* W, int L)
{
	int p_index = 0; int q_index = 0;
	for (int i = 0; i < L; i++) {
		for (int k = 1; k < L; k++)
		{
			if (W[p[p_index]] > W[p[k]]) {
				p_index = k;
			}
			if (W[q[q_index]] < W[q[k]]) {
				q_index = k;
			}
		}
		if (W[p[p_index]] < W[q[q_index]])
		{
			swap(p[p_index], q[q_index]);
			p_index = 0; q_index = 0;

		}
		else { break; }
	}
}


void decode_list(float** LLR, int** sum, int j, int* A_Ac, int a, int N, int K, float* PM, int L, int m,
	float* LLR_in, float* W, int* index ,int* better, int* worse,int* Count, int* Count_info, float* Qsumun)
{
	//if (sign == 1) return;
	/*	LLR   L*N LLR matrix as input
	sum   L*N beta matrix as output
	j     the number of stage, equaling to log2N
	A_Ac  the array which label the channel condition, frozen is 0, information is 1
	a     a variable which label the position the start of each array (alpha and beta)
	N     Codelength
	K     Information bit length
	PM    Path metric
	L     Number of paths
	m     Log2L
	p     The reference matrix in stage-located copy  */

	if (j == 0) {
		//getchar();
		//cout << "Count = " << Count << "\tis info. :" << A_Ac[Count] << endl;
		if (A_Ac[*Count] == 0)  //frozen bit condition
		{
		/*	for (int iii = 0; iii <= 2 * N - 2; iii++)
				cout << LLR[0][iii] << " ";
			cout << endl;
			cout << Count_info << endl;
			getchar();*/

			if (Count_info > 0) {
				for (int k = 0; k < L; k++) {
					if (LLR[k][2 * N - 2] < 0) { PM[k] = PM[k] + LLR[k][2 * N - 2]; }
				}
			}
			for (int k = 0; k < L; k++) {
				//cout << k << " " << Count << endl;
				sum[k][*Count] = 0;
			}
		}
		else  //information bit condition
		{
			*Count_info = *Count_info + 1;
			for (int k = 0; k < L; k++) {
				LLR_in[k] = *(*(LLR + k) + (N << 1) - 2);
			}
			if (*Count_info < m + 1)  // For the first 'm' information bits, do path expansion only
			{
				int valid = (L >> *Count_info);
				for (int k = 0; k < L; k += (valid << 1))
				{
					float PM_k = PM[k];
					for (int q = 0; q < valid; q++) {
						if (LLR_in[k] >= 0) {
							PM[k + valid + q] = PM_k - LLR_in[k];
						}
						else {
							PM[k + q] = PM_k + LLR_in[k];
							PM[k + valid + q] = PM_k;
						}
					}
				}
				for (int k = valid; k < L; k += (valid << 1))
				{
					for (int q = 0; q < valid; q++) {
						sum[k + q - valid][*Count] = 0;
						sum[k + q][*Count] = 1;
					}
				}
			}
			else if (*Count_info <= K) {  // The remaining information bits
				for (int k = 0; k < L; k++) {  //Path expansion and new path metric calculation
					if (LLR_in[k] >= 0) {
						W[(k << 1)] = PM[k];
						W[(k << 1) + 1] = PM[k] - LLR_in[k];
						better[k] = (k << 1);
						worse[k] = (k << 1) + 1;
					}
					else {
						W[(k << 1)] = PM[k] + LLR_in[k];
						W[(k << 1) + 1] = PM[k];
						worse[k] = (k << 1);
						better[k] = (k << 1) + 1;
					}
				}
				sort_list(better, worse, W, L);  // Path metric sortig
			
				for (int k = 0; k < L; k++) {
					index[k] = (better[k] >> 1);
				}
				float** tLLR = new float* [L];
				for (int i = 0; i < L; i++)  tLLR[i] = new float[N];
				int** tsum = new int* [L];
				for (int i = 0; i < L; i++)  tsum[i] = new int[N];
				for (int k = 0; k < L; k++) {  //path copy
					replace_sum(*(sum + index[k]), *(tsum + k), *Count);
					replace_LLR(*(LLR + index[k]) + N, *(tLLR + k), N);
				}
				float Qunvisit = 0;
				for (int k = 0; k < L; k++) {  //path copy
					replace_sum(*(tsum + k), *(sum + k), *Count);
					replace_LLR(*(tLLR + k), *(LLR + k) + N, N);
					sum[k][*Count] = better[k] % 2;
					PM[k] = W[better[k]];
					Qunvisit += exp(W[worse[k]]);
				}
				*Qsumun += pow(2, - (N - *Count - 1 - (K - *Count_info)))* Qunvisit;
				for (int i = 0; i < L; i++)
					delete[]tLLR[i];
				delete[] tLLR;
				for (int i = 0; i < L; i++)
					delete[]tsum[i];
				delete[] tsum;
			}
		}
		*Count= *Count+1;
		return;
	}
	int node = 1 << j;
	int length = node >> 1;
	int start = (N << 1) - node;
	int last_start = (N << 1) - (node << 1);

		for (int k = 0; k < L; k++) {
			f_function(*(LLR + k) + last_start, length);
		}

		decode_list(LLR, sum, j - 1, A_Ac, a, N, K, PM, L, m, LLR_in, W, index, better, worse, Count, Count_info, Qsumun);

		a = *Count >> (j - 1);

		for (int k = 0; k < L; k++) {
			g_function(*(LLR + k) + last_start, *(sum + k) + ((a - 1)*length), length);
		}

		decode_list(LLR, sum, j - 1, A_Ac, a, N, K, PM, L, m, LLR_in, W, index, better, worse, Count, Count_info, Qsumun);
	
		a = *Count >> (j - 1);
		int initial_posit = (a - 2) << (j - 1);
		for (int k = 0; k < L; k++) {
			if (length < 4) {
				for (int i = 0; i < length; i++)
				{
					*(*(sum + k) + initial_posit + i) ^= *(*(sum + k) + initial_posit + length + i);
				}
			}
			else {
				combine(*(sum + k) + initial_posit, *(sum + k) + initial_posit + length, length);
			}
		}
	
	return;
}
int cal_stage(int i, int n)
{
	int stage = 0;
	int N = pow(2, n);
	if (i == 0) stage = n - 1;
	else if (i < N) {
		char *str = new char[n];
		int j = 0;
		do {
			str[j++] = i % 2 + '0';
			i /= 2;
		} while (i);
		for (j = 0; j < n; j++)
		{
			if (str[j] == '1') {
				stage = j; break;
			}
		}
		delete[] str;
	}
	else stage = n - 1;

	return stage;
}
void cal_LLR(int s, float* LLR_new, float* LLR, int* psum, int fg, int N)
{
	int node = 1 << s;
	int length = node >> 1;
	int start = (N << 1) - node;
	int last_start = (N << 1) - (node << 1);
	if (fg == 0) {
		f_function_index(LLR+ last_start, LLR_new + last_start, length);
	}
	else {
		g_function_index(LLR + last_start, psum + N - node, LLR_new + last_start, length);
	}
}

void cal_sum(int s,int* psum_L, int* psum_R, int N)
{

}
