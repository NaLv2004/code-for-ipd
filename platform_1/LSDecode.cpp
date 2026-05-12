#include <stdio.h>
#include <stdlib.h> 
#include <iostream>
#include <cmath>
#include <bitset>
#include <algorithm>
#include <complex>
#include <time.h>
#include <fstream>
#include <vector>
#include "LSDecode.h"
using namespace std;
const int MAXIMUM = 1000000;
float minor = 0.000001;
int calcd2(int* A, int K, int i, int* u, int** G, int* t)
{
	int temp = 0;
	for (int j = i + 1; j < K; j++)
		if (G[A[j]][A[i]])
		{
			//	*t = *t + 1;
			temp = temp ^ (u[A[j]] & G[A[j]][A[i]]);
			//	cout << j << " "<<u[j] <<" "<< G[j][i] << endl;
		}
	return temp;
}
int calcd(int* A, int K, int i, int* u, int** G, int* t)
{
	int temp = 0;
	for (int j = i; j < K; j++)
		if (G[A[j]][A[i]])
		{
			//	*t = *t + 1;
			temp = temp ^ (u[A[j]] & G[A[j]][A[i]]);
			//	cout << j << " "<<u[j] <<" "<< G[j][i] << endl;
		}
	return temp;
}
/*
void calbestworst(float **d, int K, float *worst, vector <int> *sinset)
{
	worst[0] = 0;
	for (int j = 0;j < K;j++)
	{
		//best[j+1]=best[j];
		worst[j + 1] = worst[j];
		for (int i = 0;i < sinset[j].size();i++)
		{
			if (d[0][sinset[j][i]] < d[1][sinset[j][i]])
			{
				//best[j] = best[j]+d[0][sinset[j][i]];
				worst[j + 1] = worst[j + 1] + d[1][sinset[j][i]];
			}
			else
			{
				//best[j] = best[j]+d[1][sinset[j][i]];
				worst[j + 1] = worst[j + 1] + d[0][sinset[j][i]];
			}
		}
	}
}

void calcdistance(float *y, float **d, int N)
{
	float *y_slice = new float[N];
	for (int i = 0; i < N; i++)
		y_slice[i] = (1 + y[i]) / 2;

	for (int i = 0; i < N; i++)
	{
		d[0][i] = y_slice[i] * y_slice[i];
		d[1][i] = (y_slice[i] - 1) * (y_slice[i] - 1);
		if (d[0][i] < d[1][i])
		{
			d[1][i] = d[1][i] - d[0][i];
			d[0][i] = 0;
		}
		else
		{
			d[0][i] = d[0][i] - d[1][i];
			d[1][i] = 0;
		}
	}
	delete[]y_slice;
}*/


void calbestworst(float** d, int K, float* best, float* worst, vector <int>* sinset)
{
	best[0] = worst[0] = 0;
	for (int i = 0; i < sinset[0].size(); i++)
	{
		if (d[0][sinset[0][i]] < d[1][sinset[0][i]])
		{
			best[0] = best[0] + d[0][sinset[0][i]];
			worst[0] = worst[0] + d[1][sinset[0][i]];
		}
		else
		{
			best[0] = best[0] + d[1][sinset[0][i]];
			worst[0] = worst[0] + d[0][sinset[0][i]];
		}
	}
	for (int j = 1; j < K; j++)
	{
		best[j] = best[j - 1];
		worst[j] = worst[j - 1];
		for (int i = 0; i < sinset[j].size(); i++)
		{
			if (d[0][sinset[j][i]] < d[1][sinset[j][i]])
			{
				best[j] = best[j] + d[0][sinset[j][i]];
				worst[j] = worst[j] + d[1][sinset[j][i]];
			}
			else
			{
				best[j] = best[j] + d[1][sinset[j][i]];
				worst[j] = worst[j] + d[0][sinset[j][i]];
			}
		}
	}
}

void calcdistance(float* y, float** d, int N)
{
	float* y_slice = new float[N];
	for (int i = 0; i < N; i++)
	{
		y_slice[i] = (1 + y[i]) / 2;
	}
	//	d[0][0] = y_slice[0] * y_slice[0];
		//d[1][0] = (y_slice[0] - 1) * (y_slice[0] - 1);
	for (int i = 0; i < N; i++)
	{
		d[0][i] = y_slice[i] * y_slice[i];
		d[1][i] = (y_slice[i] - 1) * (y_slice[i] - 1);
		if (d[0][i] < d[1][i])
		{
			d[1][i] = d[1][i] - d[0][i];
			d[0][i] = 0;
		}
		else
		{
			d[0][i] = d[0][i] - d[1][i];
			d[1][i] = 0;
		}
	}
	delete[]y_slice;
}



int partition(float a[], int s, int t)
{
	int i, j;
	float x = a[s];
	i = s;
	j = t;
	do
	{
		while ((a[j] >= x) && i < j) j--;
		if (i < j) a[i++] = a[j];
		while ((a[i] <= x) && i < j) i++;
		if (i < j) a[j--] = a[i];
	} while (i < j);
	a[i] = x;
	return i;
}

float findmink(float a[], int low, int high, int k)
{
	int q;
	int index = -1;
	if (low < high)
	{
		q = partition(a, low, high);
		int len = q - low + 1;
		if (len == k)
			return a[q];
		if (len < k) return findmink(a, q + 1, high, k - len);
		return findmink(a, low, q - 1, k);
	}
	return a[high];
}



void LSDecode(int* A, float* y, int* A_Ac, int N, int K, int L, int** u, int** GG, vector <int>* sinset, int* DF, int** GGF, int breakpoint, float* d, clock_t* T_start, clock_t* T_finish, int* timesnumber)
//void LSDecode(float *y, int* A_Ac, int N, int K, int L, int **u, int **GG, int breakpoint, float *d)
{

	int* Temp = new int[N];
	int* mark = new int[L];
	float** distanceofy = new float* [2];// distance of y-(1 or 0)
	float* best = new float[K];
	float* worst = new float[K];
	distanceofy[0] = new float[N];
	distanceofy[1] = new float[N];
	calcdistance(y, distanceofy, N);
	calbestworst(distanceofy, K, best, worst, sinset);

	/*for (int i = 0; i < N; i++) cout << distanceofy[0][i] << " " << distanceofy[1][i] << "\n";
	cout<<endl;
	for (int i = 0; i < K; i++) cout << best[i] << " " << worst[i] << "\n";
	/*getchar();
	getchar();*/

	//***************************   calc the startoflist ***************
	int m = int(log(L) / log(2));
	int** startoflist = new int* [L];
	for (int i = 0; i < L; i++)
	{
		mark[i] = 0;
		startoflist[i] = new int[m];
	}
	for (int i = 0; i < m; i++) startoflist[0][i] = 0;
	for (int i = 1; i < L; i++)
	{
		startoflist[i][0] = startoflist[i - 1][0] + 1;
		int j = 0;
		while (startoflist[i][j] == 2)
		{
			startoflist[i][j] = 0;
			startoflist[i][j + 1] = startoflist[i - 1][j + 1] + 1;
			j++;
		}
		while (j < m - 1)
		{
			j++;
			startoflist[i][j] = startoflist[i - 1][j];
		}
	}
	//***************************  calc startoflist end   ***************
	//*******************************************************************
		//for (int i = 0; i < L; i++) u[i] = new int[N];
	float r = MAXIMUM;
	float* copyd = new float[2 * L];
	int Lc = 0, ans;
	float lastmin;
	*timesnumber = 0;
	//*T_start = clock();
	/*
	for (int i = K - 1; i > breakpoint; i--)
	{
		float d1, d2;
		u[0][sinset[i][0]] = 0;
		ans = calcd2(N, sinset[i][0], u[0], GG);
		//*timesnumber = *timesnumber + N - sinset[i][0]+1;
		d1 = distanceofy[ans][sinset[i][0]] + d[0];
		d2 = distanceofy[ans ^ 1][sinset[i][0]] + d[0];
		for (int kk = 1; kk < sinset[i].size(); kk++)
		{
			int ts = 0;
			for (int jj = 1; jj <= DF[sinset[i][kk]]; jj++)
				ts = ts + u[0][GGF[sinset[i][kk]][jj]];
			//*timesnumber = *timesnumber + DF[sinset[i][kk]];
			d1 = d1 + distanceofy[ts % 2][sinset[i][kk]];
			d2 = d2 + distanceofy[(ts+1) % 2][sinset[i][kk]];
		}
		d[0]=d1;
		if (d2 < d1)
		{
			u[0][sinset[i][0]] = 1;
			d[0] = d2;
		}
	}
	r = d[0];
	d[0] = 0;*/
	//float rr=r;

	//cout<<r<<endl;
	 //r=MAXIMUM;
//*******************************************************************
	float bound = r;
	*T_start = clock();
	//int mmmmm=0;
	for (int i = K - 1; i > breakpoint; i--)
	{
		/*	for (int j = 0; j < L; j++)
			{
				for (int k = i + 1;k < K;k++)
					cout << u[j][sinset[k][0]] << "　";
				cout << "PM = " << d[j]<<endl;
			}
			getchar();*/
		if (K - i <= m)
		{
			lastmin = MAXIMUM;
			for (int j = 0; j < L; j++)
			{
				u[j][sinset[i][0]] = startoflist[j][K - i - 1];
				d[j] = d[j] + distanceofy[calcd(A, K, i, u[j], GG, timesnumber)][sinset[i][0]];
				*timesnumber = *timesnumber + sinset[i].size();
				for (int kk = 1; kk < sinset[i].size(); kk++)
				{
					int ts = 0;
					for (int jj = 1; jj <= DF[sinset[i][kk]]; jj++)
						ts = ts + u[j][GGF[sinset[i][kk]][jj]];
					//	*timesnumber = *timesnumber + DF[sinset[i][kk]];
					d[j] = d[j] + distanceofy[ts % 2][sinset[i][kk]];
				}
				if (d[j] < lastmin) lastmin = d[j];
			}
			//	cout << i << " " << *timesnumber << endl;
			if (i > 0 && (lastmin + worst[i - 1]) < bound) bound = lastmin + worst[i - 1];
			if (K - i == m) Lc = L;
		}
		else
		{
			float d1, d2;
			int ans = 0;
			lastmin = MAXIMUM;
			for (int j = 0; j < L; j++)
				if (mark[j] != 1)
				{
					if (mark[j] == 2) mark[j] = 0;
					u[j][sinset[i][0]] = 0;
					ans = calcd2(A, K, i, u[j], GG, timesnumber);
					*timesnumber = *timesnumber + sinset[i].size();
					d1 = distanceofy[ans][sinset[i][0]] + d[j];
					d2 = distanceofy[ans ^ 1][sinset[i][0]] + d[j];
					//cout << "ans= " << ans << " " << d1 << " " << d2 << " " << sinset[i].size() << endl;
					for (int kk = 1; kk < sinset[i].size(); kk++)
					{
						int ts = 0;
						//cout << "frozen= " << sinset[i][kk] << " ";
						for (int jj = 1; jj <= DF[sinset[i][kk]]; jj++)
						{
							ts = ts + u[j][GGF[sinset[i][kk]][jj]];
							//		cout << "frozen impact= " << GGF[sinset[i][kk]][jj] << "    ";
						}
						//*timesnumber = *timesnumber + DF[sinset[i][kk]];
						d1 = d1 + distanceofy[ts % 2][sinset[i][kk]];
						d2 = d2 + distanceofy[(ts + 1) % 2][sinset[i][kk]];
						//cout << d1 << "　" << d2 << endl;
					}
					if (d1 < lastmin) lastmin = d1;
					if (d2 < lastmin) lastmin = d2;
					copyd[j] = d1;
					copyd[j + L] = d2;
					d[j] = d1;
					d[j + L] = d2;
				}
				else
				{
					d[j] = MAXIMUM;
					d[j + L] = MAXIMUM;
					copyd[j] = MAXIMUM;
					copyd[j + L] = MAXIMUM;
				}
			if (i > 0 && (lastmin + worst[i - 1]) < bound) bound = lastmin + worst[i - 1];
			int Ldone = 2 * Lc;
			if (Ldone > L) Ldone = L;

			/*cout << "Middle " << endl;
			for (int j = 0; j < 2*L; j++)
			{
			cout << "PM = " << d[j] << endl;
}
				cout << *timesnumber << endl;
getchar();*/
			float Lmin = findmink(d, 0, 2 * L - 1, Ldone);

			//		if (i>0 && Lmin+best[i-1] > bound)
			//		{
					/*	cout << i << "\n";
						cout << Lmin + best[i - 1] << "\n";
						cout << bound << "\n";*/
						//	getchar();
					//		Lmin = bound - best[i - 1]+0.00001;		
						//		cout << "！！！！！！\n";

							//mmmmm=1;

				//		}
						/*	cout<<bound<<endl;
							for (int j = 0; j < 2*L; j++)
							{
								cout << copyd[j] << " ";
							}
							cout << "\n";
								cout << Lmin<<"\n";
						if (mmmmm)						getchar();*/
			int j = 0, k = 0;

			//		cout << Lmin << "\n"; 
			Lc = 0;
			while (j < L && Lc < Ldone)
			{
				if (copyd[j] <= Lmin && copyd[j + L] <= Lmin)
				{
					u[j][sinset[i][0]] = 0;
					d[j] = copyd[j];
					while (k < L && (copyd[k] <= Lmin || copyd[k + L] <= Lmin)) k++;
					if (k != L && Lc + 1 < Ldone)
					{
						for (int p = i + 1; p < K; p++) u[k][sinset[p][0]] = u[j][sinset[p][0]];
						u[k][sinset[i][0]] = 1;
						d[k] = copyd[j + L];
						mark[k] = 2;
						k++;
						Lc++;
						//*expensionnum = *expensionnum + 1;
					}
					j++;
					Lc++;
				}
				else if (copyd[j] <= Lmin)
				{
					u[j][sinset[i][0]] = 0;
					d[j] = copyd[j];
					j++;
					Lc++;
				}
				else if (copyd[j + L] <= Lmin)
				{
					u[j][sinset[i][0]] = 1;
					d[j] = copyd[j + L];
					j++;
					Lc++;
				}
				else
				{
					if (mark[j] != 2)
					{
						d[j] = MAXIMUM;
						mark[j] = 1;
					}
					j++;
				}
			}
			/*	while (j < L)
				{
					if (copyd[j] <= Lmin)
					{
						u[j][sinset[i][0]] = 0;
						d[j] = copyd[j];
						j++;
					}
					else if (copyd[j + L] <= Lmin)
						{
							u[j][sinset[i][0]] = 1;
							d[j] = copyd[j+L];
							j++;
						}
					else
					{
						while (copyd[k] > Lmin || copyd[k + L] > Lmin) k++;
						for (int p = i + 1; p < K; p++) u[j][sinset[p][0]] = u[k][sinset[p][0]];
						u[j][sinset[i][0]] = 1;
						d[j] = copyd[k + L];
						k++;
						j++;
					}
				}*/


				/*if (Lc<L)
				{
					cout<<Lc<<" "<<i<<" !!!!!";
					getchar();
				}*/
		}

		/*for (int j = 0; j < L; j++)
		{
			cout << d[j] << " ";
					//	cout << mark[j] << " ";
		}		cout << "\n";
		for (int j = 0; j < L; j++)
		{
		//	cout << d[j] << " ";
						cout << mark[j] << " ";
		}
		cout << "\n";
		cout << i <<" "<< Lc<<" "<<bound<<" done!\n";
		if (mmmmm==1) getchar();*/
	}
	*T_finish = clock();
	/*cout << min << "\n";
	float s = 0;
	for (int i = N - 1; i >= 0; i--)
		s = s + distanceofy[calcd(N, i, uu, GG)][i];
	cout << s << "\n";
	getchar();*/
	/*float mmin=MAXIMUM;
	for (int j = 0; j < L; j++)
		if (d[j]<mmin) mmin=d[j];
	if (mmin>rr)
		{cout << "FFFFF\n";
	getchar();}*/
	for (int i = 0; i < 2; i++)
		delete[]distanceofy[i];
	delete[] distanceofy;
	for (int i = 0; i < L; i++)
		delete[]startoflist[i];
	delete[] startoflist;
	delete[] copyd;
	//delete[] DFi;
	delete[] Temp;
	delete[] mark;
	delete[] best;
	delete[] worst;
}





















//#include <stdio.h>
//#include <stdlib.h> 
//#include <iostream>
//#include <cmath>
//#include <bitset>
//#include <algorithm>
//#include <complex>
//#include <time.h>
//#include <fstream>
//
//using namespace std;
//
//
//int calcd2(int CodeLength, int i, int* u, int **G)
//{
//	int temp = 0;
//	for (int j = i+1; j < CodeLength; j++)
//		temp = temp ^ (u[j] & G[j][i]);
//	return temp;
//}
//int calcd(int CodeLength, int i, int* u, int **G)
//{
//	int temp = 0;
//	for (int j = i; j < CodeLength; j++)
//		temp = temp ^ (u[j] & G[j][i]);
//	return temp;
//}
//
//void calcdistance(float *y, float **d,int N)
//{
//	float *y_slice= new float[N];
//	for (int i = 0; i < N; i++)
//	{
//		y_slice[i] = (1 + y[i]) / 2;
//	}
//	for (int i = 0; i < N; i++)
//	{
//		d[0][i] = y_slice[i] * y_slice[i];
//		d[1][i] = (y_slice[i]-1) * (y_slice[i]-1);
//	}
//	delete []y_slice;
//}
//
//int partition(float a[], int s, int t)
//{
//	int i, j;
//	float x = a[s];
//	i = s;
//	j = t;
//	do
//	{
//		while ((a[j] >= x) && i < j) j--;
//		if (i < j) a[i++] = a[j];
//		while ((a[i] <= x) && i < j) i++;
//		if (i < j) a[j--] = a[i];
//	} while (i < j);
//	a[i] = x;
//	return i;
//}
//
//float findmink(float a[], int low, int high, int k)
//{
//	int q;
//	int index = -1;
//	if (low < high)
//	{
//		q = partition(a, low, high);
//		int len = q - low + 1;
//		if (len == k)
//			return a[q];
//		if (len < k) return findmink(a,q+1,high,k-len);
//		return findmink(a,low,q-1,k);
//	}
//	return a[high];
//}
//
//void LSDecode(float *y, int* A_Ac, int N, int K,int L,int **u,int **GG,int breakpoint,float *d, clock_t *T_start, clock_t *T_finish)
//{
//
//	float **distanceofy = new float*[2];// distance of y-(1 or 0)
//	distanceofy[0] = new float[N];
//	distanceofy[1] = new float[N];
//	calcdistance(y, distanceofy,N);
//	//for (int i = 0; i < N; i++) cout << distanceofy[0][i] << " " << distanceofy[1][i] << "\n";
////***************************   calc the startoflist ***************
//	int m = int(log(L) / log(2));
//	int **startoflist = new int *[L];
//	for (int i=0;i<L;i++)	startoflist[i] = new int[m];
//	for (int i = 0; i < m; i++) startoflist[0][i] = 0;
//	for (int i = 1; i < L; i++)
//	{
//		startoflist[i][0] = startoflist[i - 1][0] + 1;
//		int j = 0;
//		while (startoflist[i][j] == 2)
//		{
//			startoflist[i][j] = 0;
//			startoflist[i][j + 1] = startoflist[i - 1][j+1] + 1;
//			j++;
//		}
//		while (j<m-1)
//		{
//			j++;
//			startoflist[i][j] = startoflist[i - 1][j];
//		}
//	}
////***************************  calc startoflist end   ***************
//
//	//int **u = new int *[L];
//	//for (int i = 0; i < L; i++) u[i] = new int[N];
//
//	float *copyd = new float[2 * L];
//	int randbit = 0;
//	*T_start = clock();
//	for (int i = N - 1; i > breakpoint; i--)
//	{
//		if (A_Ac[i] == 1)
//		{
//			randbit++;
//			if (randbit<=m)
//				for (int j = 0; j < L; j++)
//				{
//					u[j][i] = startoflist[j][randbit-1];
//					d[j] = d[j] + distanceofy[calcd(N, i,u[j], GG)][i];
//				}
//			else
//			{
//				float d1, d2;
//				int ans = 0;
//				for (int j = 0; j < L;j++)
//				{
//					u[j][i] = 0;
//					ans = calcd2(N, i, u[j], GG);
//					d1 = distanceofy[ans][i] + d[j];
//					d2 = distanceofy[ans ^ (GG[i][i])][i] + d[j];
//					copyd[j] = d1;
//					copyd[j + L] = d2;
//					d[j] = d1;
//					d[j + L] = d2;
//				}
//				/*for (int j = 0; j < 2*L; j++)
//				{
//					cout << d[j] << " ";
//				}
//				cout << "\n";
//				getchar();*/
//				float Lmin = findmink(d, 0, 2*L-1, L);
//				int j = 0, k = 0;
//				while (j < L)
//				{
//					if (copyd[j] <= Lmin)
//					{
//						u[j][i] = 0;
//						d[j] = copyd[j];
//						j++;
//					}
//					else if (copyd[j + L] <= Lmin)
//						{
//							u[j][i] = 1;
//							d[j] = copyd[j+L];
//							j++;
//						}
//					else
//					{
//						while (copyd[k] > Lmin || copyd[k + L] > Lmin) k++;
//						for (int p = i + 1; p < N; p++) u[j][p] = u[k][p];
//						u[j][i] = 1;
//						d[j] = copyd[k + L];
//						k++;
//						j++;
//					}
//				}
//			}
//		}
//		else
//		{
//			for (int j = 0; j < L; j++)
//			{
//				u[j][i] = 0;
//				d[j] = distanceofy[calcd(N, i,u[j], GG)][i] + d[j];
//			}
//		}
//		/*for (int j = 0; j < L; j++)
//		{
//			cout << d[j] << " ";
//		}
//		cout << "\n";
//
//		getchar();*/
//	//	cout << i <<"done!\n";
//	}
//	*T_finish = clock();
//	/*cout << min << "\n";
//	float s = 0;
//	for (int i = N - 1; i >= 0; i--)
//		s = s + distanceofy[calcd(N, i, uu, GG)][i];
//	cout << s << "\n";
//	getchar();*/
//
//	for (int i = 0; i <2; i++)
//		delete[]distanceofy[i];
//	delete[] distanceofy;
//	for (int i = 0; i <L; i++)
//		delete[]startoflist[i];
//	delete[] startoflist;
//
//	delete[] copyd;
//}
