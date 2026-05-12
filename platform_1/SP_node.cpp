#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <bitset>
#include <vector>
#include <functional>
#include <queue>
#include "Polar_Encoder.h"
using namespace std;
#define hard(n) (n>=0?0:1)
struct Comp {
	Comp(const vector<float>& v) : _v(v) {}
	bool operator ()(int a, int b) { return _v[a] > _v[b]; }
	const vector<float>& _v;
};
void Polar_Encoder_bitset(bitset<16>* u_mat, int len) {
	switch (len) {
	case 2:
	{
		*(u_mat) ^= *(u_mat + 1);
		break;
	}
	case 4:
	{
		for (int i = 0; i < 4; i += 4)
		{
			*(u_mat + 2) ^= *(u_mat + 3);
			*u_mat ^= (*(u_mat + 1) ^ *(u_mat + 2));
			*(u_mat + 1) ^= *(u_mat + 3);
		}
		break;
	}
	case 8:
	{
		for (int i = 0; i < 8; i += 4)
		{
			*(u_mat + i + 2) ^= *(u_mat + i + 3);
			*(u_mat + i) ^= (*(u_mat + i + 1) ^ *(u_mat + i + 2));
			*(u_mat + i + 1) ^= *(u_mat + i + 3);
		}
		*(u_mat) ^= *(u_mat + 4);
		*(u_mat + 1) ^= *(u_mat + 5);
		*(u_mat + 2) ^= *(u_mat + 6);
		*(u_mat + 3) ^= *(u_mat + 7);
		break;
	}
	case 16:
	{
		for (int i = 0; i < 16; i += 4)
		{
			*(u_mat + i + 2) ^= *(u_mat + i + 3);
			*(u_mat + i) ^= (*(u_mat + i + 1) ^ *(u_mat + i + 2));
			*(u_mat + i + 1) ^= *(u_mat + i + 3);
		}
		for (int i = 0; i < 4; i++)
		{
			*(u_mat + i + 8) ^= *(u_mat + i + 12);
			*(u_mat + i) ^= (*(u_mat + i + 4) ^ *(u_mat + i + 8));
			*(u_mat + i + 4) ^= *(u_mat + i + 12);
		}
		break;
	}
	default: {cout << "error!" << endl; break; }
	}
}
void cal_NM_DMM(float* LLR, int Iv, int node, vector<int> A_pos, int** sum_all, float* NM_all, int qk, bitset<16> * u_mat, vector<float> NM, vector<int> res)
{
	int num = 1 << Iv;
	int i, j, k;
	memset(sum_all[0], 0, node*sizeof(float));
	for (i = 0; i < 16; i++)u_mat[i].reset();
	for (j = 0; j < num; j++) {
		bitset<16> bit(j);
		for (int i = 0; i < Iv; i++)
		{
			if (bit.test(i)) u_mat[A_pos[i]].set(j);
		}
	}
	Polar_Encoder_bitset(u_mat, node);
	for (i = 0; i < node; i++) {
		if (LLR[i] > 0)
		{
			for (j = 0; j < num; j++) { if (u_mat[i].test(j))NM[j] -= LLR[i]; }
		}
		else
		{
			for (j = 0; j < num; j++) { if (!u_mat[i].test(j))NM[j] += LLR[i]; }
		}
	}
	for (i = 0; i < num; i++) {
		res[i] = i;
	}
	if (qk < num) nth_element(res.begin(), res.begin() + qk, res.end(), Comp(NM));
	for (j = 0; j < qk; j++)
	{
		for (i = 0; i < node; i++)
		{
			sum_all[j][i] = u_mat[i].test(res[j]);
		}
		NM_all[j] = NM[res[j]];
	}
}