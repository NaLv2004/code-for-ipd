// Perform 1D/2D encoding and decoding.
//1D polar  copyright(c) 2017 by Yifei Shen - LEADS - SEU
//2D system copyright(c) 2024 by Huayi Zhou, Jiayan Xu - LEADS - SEU

//edit setting.h to change setting
//outputs in result.txt (see system.cpp)

#include "incluall.h"
#include<iomanip>
// #include "system_new.cpp"

#define _CRTDBG_MAP_ALLOC
#define EIGEN_USE_MKL_ALL
#include <stdlib.h>
#include <crtdbg.h>

using namespace std;
//beta[betalen] = { 0.2,0.2,0.6,0.8,1,1,1,1 };
// test interleaving
vector<int> pattern(5);
int A[5] = { 7,8,9,5,3 };
float B[5] = { 0.7,0.8,0.9,0.5,0.3 };



//#ifndef _DEBUG
//#define _DEBUG
//#endif
//
//#ifdef _DEBUG
//#undef new
//#define new new(_NORMAL_BLOCK, __FILE__, __LINE__)
//#endif

// test interleaving
int main()
{
	float MMM = 1.1;
	float PPP = MMM * 4 / 3 ;
	cout << PPP << endl;
	int** test_arr = new int*[5];
	for (int i = 0; i < 5; i++)
	{
		test_arr[i] = new int[5];
		for (int j = 0; j < 5; j++)
		{
			test_arr[i][j] = i*100 + j;
		}
	}

	int* p = test_arr[2] + 1;
	cout << p[2] << endl;
	//_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
	generateIntrlvPattern(pattern, 5);
	for (int i = 0; i < 5; i++) { cout << B[i]; }
	cout << endl;
	randIntrlv(pattern, B, 5);
	for (int i = 0; i < 5; i++) { cout << B[i]; }
	cout << endl;
	for (int i = 0; i < 5; i++) { cout << pattern[i]; }
	cout << endl;
	randDeintrlv(pattern, B, 5);
	for (int i = 0; i < 5; i++) { cout << B[i] << "    "; }
	cout << endl;
	/*double* BER_array112 = new double[10];
	#pragma omp parallel for
	for (int i = 0; i < 10; i++)
	{
		double* BER_array11 = new double[10];
		for (int j=0;j<10;j++)
			BER_array11[i] += j+i;
		BER_array112[i] = BER_array11[i];
	}
	for (int i = 0; i < 10; i++)
		cout << BER_array112[i] << endl;
	getchar();*/
	//freopen("sysresult.txt", "w", stdout);
	int SNR_len = (int)(Max_SNR - Min_SNR) / SNR_step;
	double * BER_array = new double[SNR_len];
	double * FER_array = new double[SNR_len];
	//printsetting();
	if (system_archi == "JDD") system2D_JDD_copy(BER_array, FER_array);
	if (dimensions==  1 && polarde=="SLD") system1D_LSD_MIMO(BER_array, FER_array);
	if (dimensions == 1 && polarde == "SCL") system1D(BER_array, FER_array);
	if (dimensions == 2 && polarde == "SCL") system2D(BER_array, FER_array);
	if (dimensions == 2 && polarde == "SLD" && sys==true) system2D_LSD_epdet_sys(BER_array, FER_array);
	if (dimensions == 2 && polarde == "SLD" && sys == false && MIMODET!="EP" && irregular_coding) system2D_LSD_epdet_irregular_SVD_precode(BER_array, FER_array);
	if (dimensions == 2 && polarde == "SLD" && sys == false && MIMODET != "EP" && (!irregular_coding)) system2D_LSD_epdet(BER_array, FER_array);
	if (dimensions == 2 && polarde == "SLD" && sys == false && MIMODET == "EP") system2D_LSD_epdet(BER_array, FER_array);
	// if (dimensions == 2 && polarde == "SLD" && sys == false && double_stage == 1) system2D_LSD_doublestage(BER_array, FER_array);
	// std::cout.clear();
	std::cout << "SNR	" << "FER  " << endl;
	for (int i = 0; i < SNR_len+1; i++) {
		std::cout << Min_SNR + i * SNR_step << "\t";
		std::cout << FER_array[i] << "\n";
	}
	cout << endl; cout << endl;
	cout << "SNR     " << "BER  " << endl;
	for (int i = 0; i < SNR_len + 1; i++) {
		cout << Min_SNR + i * SNR_step << "\t";
		cout << BER_array[i]<<"\n";
	}
	cout << endl;

	delete[] BER_array;
	delete[] FER_array;
	//fclose(stdout);
	//getchar();
	//getchar();
	//int* MMM = new int[100];
    //_CrtDumpMemoryLeaks();
	return 0;
}


// single simulation BLER-EbN0
//int main() {
//	sPhyConfig* PhyConfig = new sPhyConfig;
//	PhyConfig->nN = N;
//	PhyConfig->nK = K;
//	PhyConfig->polarCon = polarcon;
//	PhyConfig->nNh = Nh;
//	PhyConfig->nNv = Nv;
//	PhyConfig->nKh = Kh;
//	PhyConfig->nKv = Kv;
//	PhyConfig->nLowRateCodewords = border;
//	PhyConfig->nModType = ModType;
//	PhyConfig->nNt = Nt;
//	PhyConfig->nNr = Nr;
//	PhyConfig->nStream = trans_stream;
//	PhyConfig->bFlagInterleave = flag_interleave;
//	PhyConfig->bFlagInterleaveSpatial = flag_interleave_spatial;
//	PhyConfig->bFlagInterleaveTemporal = flag_interleave_all;
//	PhyConfig->bIrregularCoding = irregular_coding;
//	PhyConfig->powerAllocationScheme = MIMO_SVD_POWER_ALLOCATION_MODE;
//	PhyConfig->powerAllocationRatio = MIMO_SVD_SEPARATE_POWER_RATIO;
//	PhyConfig->nCRCL = CRCL;
//	PhyConfig->minSNR = 7;
//	PhyConfig->maxSNR = 300;
//	PhyConfig->stepSNR = 0.5;
//	PhyConfig->L = L;
//	PhyConfig->nMaxFrameError = errframe;
//	PhyConfig->pKTime = new int[num_time_constructions];
//	PhyConfig->polarConh = polarconh;
//	PhyConfig->polarConV = polarconv;
//	for (int i = 0; i < num_time_constructions; i++) { PhyConfig->pKTime[i] = K_time_array[i]; }
//	int lenSNRArr = round((Max_SNR - Min_SNR) / SNR_step) + 1;
//	float* BER = new float[lenSNRArr];
//	float* FER = new float[lenSNRArr];
//	printsetting();
//	/*system2D_new(PhyConfig, BER, FER);
//	system2D_new(PhyConfig, BER, FER);*/
//	system1D_new(PhyConfig, BER, FER);
//	delete[] BER;
//	delete[] FER;
//}

// transmission rate simulation: rate-EbN0 for target BLER
//int main() {
//	int numOpt = 1E4;
//	sPhyConfig* PhyConfig = new sPhyConfig;
//	soptParams* OptParam = new soptParams[numOpt];
//	float error = 0.3;
//	PhyConfig->nNh = Nh;
//	PhyConfig->nNv = Nv;
//	PhyConfig->nKh = Kh;
//	PhyConfig->nKv = Kv;
//	PhyConfig->nLowRateCodewords = border;
//	PhyConfig->nModType = ModType;
//	PhyConfig->nNt = Nt;
//	PhyConfig->nNr = Nr;
//	PhyConfig->nStream = trans_stream;
//	PhyConfig->bFlagInterleave = flag_interleave;
//	PhyConfig->bFlagInterleaveSpatial = flag_interleave_spatial;
//	PhyConfig->bFlagInterleaveTemporal = flag_interleave_all;
//	PhyConfig->bIrregularCoding = irregular_coding;
//	PhyConfig->powerAllocationScheme = MIMO_SVD_POWER_ALLOCATION_MODE;
//	PhyConfig->nCRCL = CRCL;
//	PhyConfig->minSNR = Min_SNR;
//	PhyConfig->maxSNR = Max_SNR;
//	PhyConfig->stepSNR = SNR_step;
//	PhyConfig->L = L;
//	PhyConfig->nMaxFrameError = errframe;
//	PhyConfig->pKTime = new int[num_time_constructions];
//	PhyConfig->polarConh = polarconh;
//	PhyConfig->polarConV = polarconv;
//	PhyConfig->nTimeConstructions = 2;
//	int minS = 15; // 22,32
//	int maxS = 29;
//	int nOpt = 0;
//	int border_max_1 = 23;
//	int border_max_2 = 18;
//	int K_time_array_1[2] = { 26, 31 };
//	int K_time_array_2[2] = { 16, 26 };
//	float targetErrorRate = 1E-2;
//	float designSNR = 9.5;
//	PhyConfig->minSNR = designSNR;
//	PhyConfig->maxSNR = designSNR;
//	ofstream result_file("OptimizationResult.txt");
//	std::cout << std::endl << setw(10) << "Idx Opt" << setw(10) << "Kv" << setw(10) << "Klow" << setw(10) << "Khigh" << setw(10)
//		<< "nLowRate" << setw(10) << "S" << setw(10) << "Rtrans" << setw(10) << "FER" <<std::endl;
//	// all possible stream count
//	for (int nSelStream = maxS; nSelStream >= minS; nSelStream--) {
//		float* FER = new float[1];
//		float* BER = new float[1];
//		int isErrorSatisfy = 0;
//		PhyConfig->nStream = nSelStream;
//		// (32, 26) * (26 + 31);
//		for (int i_border_1 = 0; i_border_1 < border_max_1; i_border_1+=1) {
//			if (i_border_1 == 0) {  // equal rates
//				PhyConfig->nKv = 26;
//				PhyConfig->pKTime[0] = 31;
//				PhyConfig->pKTime[1] = 31;
//				PhyConfig->powerAllocationScheme = "UNIFORM";
//				PhyConfig->bFlagInterleaveSpatial = true;
//				PhyConfig->bFlagInterleaveTemporal = true;
//				float designEbN0 = SNR2EbN0(designSNR, PhyConfig);
//				PhyConfig->minSNR = designEbN0;
//				PhyConfig->maxSNR = designEbN0;
//				system2D_new(PhyConfig, BER, FER);
//				float nowFER = FER[0];
//				if (nowFER < targetErrorRate*(1+error)) {
//					isErrorSatisfy = 1;
//				}
//			}
//			/*if (isErrorSatisfy) {
//				break;
//			}*/
//			else if (! isErrorSatisfy) { // unequal rates
//				PhyConfig->nKv = 26;
//				PhyConfig->pKTime[0] = K_time_array_1[0];
//				PhyConfig->pKTime[1] = K_time_array_1[1];
//				PhyConfig->nLowRateCodewords = i_border_1;
//				PhyConfig->powerAllocationScheme = "SEPARATE";
//				PhyConfig->bFlagInterleaveSpatial = false;
//				PhyConfig->bFlagInterleaveTemporal = true;
//				float designEbN0 = SNR2EbN0(designSNR, PhyConfig);
//				PhyConfig->minSNR = designEbN0;
//				PhyConfig->maxSNR = designEbN0;
//				system2D_new(PhyConfig, BER, FER);
//				float nowFER = FER[0];
//				// std::cout << endl << "nowFER" << nowFER << "TEP" << targetErrorRate << std::endl;
//				if (nowFER < targetErrorRate*(1+error)) {
//					// std::cout << "YYYYYYY" << std::endl;
//					isErrorSatisfy = 1;
//				}
//			}
//			setOptParam(PhyConfig, OptParam + nOpt, FER[0]);
//			std::cout << std::endl << setw(10) << nOpt << setw(10) << OptParam[nOpt].nKv << setw(10) << OptParam[nOpt].pKTime[0] << setw(10) << OptParam[nOpt].pKTime[1] << setw(10)
//				<< OptParam[nOpt].border << setw(10) << OptParam[nOpt].nStream << setw(10) << OptParam[nOpt].transmissionRate << setw(10) << FER[0] << std::endl;
//			result_file << std::endl << setw(10) << nOpt << setw(10) << OptParam[nOpt].nKv << setw(10) << OptParam[nOpt].pKTime[0] << setw(10) << OptParam[nOpt].pKTime[1] << setw(10)
//				<< OptParam[nOpt].border << setw(10) << OptParam[nOpt].nStream << setw(10) << OptParam[nOpt].transmissionRate << setw(10) << FER[0] << std::endl;
//			nOpt++;
//			if (isErrorSatisfy) {
//				// std::cout << "YYYYYYY" << std::endl;
//				break;
//			}
//		}
//
//		if (!isErrorSatisfy) {
//			for (int i_border_1 = 0; i_border_1 < border_max_2; i_border_1 += 1) {
//				if (i_border_1 == 0 && (!isErrorSatisfy)) {  // equal rates
//					PhyConfig->nKv = 26;
//					PhyConfig->pKTime[0] = 26;
//					PhyConfig->pKTime[1] = 26;
//					PhyConfig->powerAllocationScheme = "UNIFORM";
//					PhyConfig->bFlagInterleaveSpatial = true;
//					PhyConfig->bFlagInterleaveTemporal = true;
//					float designEbN0 = SNR2EbN0(designSNR, PhyConfig);
//					PhyConfig->minSNR = designEbN0;
//					PhyConfig->maxSNR = designEbN0;
//					system2D_new(PhyConfig, BER, FER);
//					float nowFER = FER[0];
//					if (nowFER < targetErrorRate*(1+error)) {
//						isErrorSatisfy = 1;
//					}
//				}
//				else if (!isErrorSatisfy) { // unequal rates
//					PhyConfig->nKv = 26;
//					PhyConfig->pKTime[0] = K_time_array_2[0];
//					PhyConfig->pKTime[1] = K_time_array_2[1];
//					PhyConfig->nLowRateCodewords = i_border_1;
//					PhyConfig->powerAllocationScheme = "SEPARATE";
//					PhyConfig->bFlagInterleaveSpatial = false;
//					PhyConfig->bFlagInterleaveTemporal = true;
//
//					/*setOptParam(PhyConfig, OptParam + nOpt);
//					std::cout << std::endl << setw(10) << nOpt << setw(10) << OptParam[nOpt].nKv << setw(10) << OptParam[nOpt].pKTime[0] << setw(10) << OptParam[nOpt].pKTime[1] << setw(10)
//						<< OptParam[nOpt].border << setw(10) << OptParam[nOpt].nStream << setw(10) << OptParam[nOpt].transmissionRate << setw(10) << FER[0] << std::endl;*/
//					float designEbN0 = SNR2EbN0(designSNR, PhyConfig);
//					PhyConfig->minSNR = designEbN0;
//					PhyConfig->maxSNR = designEbN0;
//					system2D_new(PhyConfig, BER, FER);
//					float nowFER = FER[0];
//					if (nowFER < targetErrorRate*(1+error)) {
//						isErrorSatisfy = 1;
//					}
//				}
//				setOptParam(PhyConfig, OptParam + nOpt, FER[0]);
//				std::cout << std::endl << setw(10) << nOpt << setw(10) << OptParam[nOpt].nKv << setw(10) << OptParam[nOpt].pKTime[0] << setw(10) << OptParam[nOpt].pKTime[1] << setw(10)
//					<< OptParam[nOpt].border << setw(10) << OptParam[nOpt].nStream << setw(10) << OptParam[nOpt].transmissionRate << setw(10) << FER[0] << std::endl;
//				result_file << std::endl << setw(10) << nOpt << setw(10) << OptParam[nOpt].nKv << setw(10) << OptParam[nOpt].pKTime[0] << setw(10) << OptParam[nOpt].pKTime[1] << setw(10)
//					<< OptParam[nOpt].border << setw(10) << OptParam[nOpt].nStream << setw(10) << OptParam[nOpt].transmissionRate << setw(10) << FER[0] << std::endl;
//				nOpt++;
//				if (isErrorSatisfy) {
//					break;
//				}
//			}
//		}
//	}
//	float rate_max = 0;
//	for (int iOpt = 0; iOpt < nOpt; iOpt++) {
//		if ((OptParam[iOpt].transmissionRate > rate_max) && (OptParam[iOpt].errorRate < targetErrorRate)) {
//			result_file << "Optimal Rate and Stream Selection" << std::endl;
//			result_file << std::endl << setw(10) << iOpt << setw(10) << OptParam[iOpt].nKv << setw(10) << OptParam[iOpt].pKTime[0] << setw(10) << OptParam[iOpt].pKTime[1] << setw(10)
//				<< OptParam[iOpt].border << setw(10) << OptParam[iOpt].nStream << setw(10) << OptParam[iOpt].transmissionRate << setw(10) << OptParam[iOpt].errorRate << std::endl;
//		}
//	}
//	result_file.close();
//}
//
//

// EbN0-rate for target BLER
//int main() {
//	std::ofstream ofile("EbN0VsCodeRate.txt");
//	float targetErrorRate = 1E-3;
//	float* codeRateArray = new float[100];
//	float* EbN0Array = new float[100];
//	int Kv_array_all[30] = { 0 };
//	int Kh_low_array_all[30] = { 0 };
//	int Kh_high_array_all[30] = { 0 };
//	// int is_irregular_all[30] = { 0 };
//	for (int i = 0; i < 30; i++) { Kh_low_array_all[i] = 16; Kh_high_array_all[i] = 26; }
//	Kh_low_array_all[0] = 26;
//	Kh_low_array_all[15] = 26;
//	for (int i = 0; i < 15; i++) { Kv_array_all[i] = 15;  }
//	for (int i = 15; i < 30; i++) { Kv_array_all[i] = 11; }
//	float minSNRAll = 5.0;
//	float maxSNRAll = 15.0;
//	float stepSNRAll = 0.2;
//	sPhyConfig* PhyConfig = new sPhyConfig;
//	PhyConfig->nNh = Nh;
//	PhyConfig->nNv = Nv;
//	PhyConfig->nKh = Kh;
//	PhyConfig->nKv = Kv;
//	PhyConfig->nLowRateCodewords = border;
//	PhyConfig->nModType = ModType;
//	PhyConfig->nNt = Nt;
//	PhyConfig->nNr = Nr;
//	PhyConfig->nStream = trans_stream;
//	PhyConfig->bFlagInterleave = flag_interleave;
//	PhyConfig->bFlagInterleaveSpatial = flag_interleave_spatial;
//	PhyConfig->bFlagInterleaveTemporal = flag_interleave_all;
//	PhyConfig->bIrregularCoding = irregular_coding;
//	PhyConfig->powerAllocationScheme = MIMO_SVD_POWER_ALLOCATION_MODE;
//	PhyConfig->nCRCL = CRCL;
//	PhyConfig->minSNR = 5.0;
//	PhyConfig->maxSNR = 12.0;
//	PhyConfig->stepSNR = SNR_step;
//	PhyConfig->L = L;
//	PhyConfig->nMaxFrameError = errframe;
//	PhyConfig->pKTime = new int[num_time_constructions];
//	PhyConfig->polarConh = polarconh;
//	PhyConfig->polarConV = polarconv;
//	int nLowRate = 0;
//	for (int i = 0; i < num_time_constructions; i++) { PhyConfig->pKTime[i] = K_time_array[i]; }
//	for (int i_rate = 16; i_rate < 25; i_rate+=2) {
//		float* FER = new float[1];
//		float* BER = new float[1];
//		nLowRate+=2;
//		PhyConfig->nKv = Kv_array_all[i_rate];
//		PhyConfig->pKTime[0] = Kh_low_array_all[i_rate];
//		PhyConfig->pKTime[1] = Kh_high_array_all[i_rate];
//		if (PhyConfig->pKTime[0] == PhyConfig->pKTime[1]) {
//			PhyConfig->nLowRateCodewords = 2;
//			PhyConfig->bFlagInterleaveSpatial = true;
//			PhyConfig->powerAllocationScheme = "UNIFORM";
//			nLowRate = 0;
//		}
//		else {
//			PhyConfig->nLowRateCodewords = nLowRate;
//			PhyConfig->bFlagInterleaveSpatial = false;
//			PhyConfig->powerAllocationScheme = "SEPARATE";
//			// nLowRate = 0;
//		}
//		// float nowCodeRate = ()
//		int nInfo = PhyConfig->pKTime[0] * PhyConfig->nLowRateCodewords + PhyConfig->pKTime[1] * (PhyConfig->nKv - PhyConfig->nLowRateCodewords) - CRCL;
//		int codeLength = PhyConfig->nNh * PhyConfig->nNv;
//		float nowCodeRate = (float)nInfo / (float)codeLength;
//		for (float nowEbN0 = minSNRAll; nowEbN0 < maxSNRAll; nowEbN0 += stepSNRAll) {
//			PhyConfig->minSNR = nowEbN0;
//			PhyConfig->maxSNR = nowEbN0;
//			system2D_new(PhyConfig, BER, FER);
//			if ((FER[0] < targetErrorRate * 1.3) && (FER[0] > 0)) {
//				EbN0Array[i_rate] = nowEbN0;
//				codeRateArray[i_rate] = nowCodeRate;
//				ofile << nowCodeRate << "  " << nowEbN0 << std::endl;
//				std::cout << std::endl << nowCodeRate << "  " << nowEbN0 << std::endl;
//				break;
//			}
//		}
//		delete[] FER;
//		delete[] BER;
//	}
//	ofile.close();
//	std::cout << std::endl;
//	for (int i = 0; i < 25; i++) {
//		std::cout << std::endl << codeRateArray[i] << "  " << EbN0Array[i] << std::endl;
//	}
//	/*int lenSNRArr = round((Max_SNR - Min_SNR) / SNR_step) + 1;
//	float* BER = new float[lenSNRArr];
//	float* FER = new float[lenSNRArr];
//	printsetting();
//	system2D_new(PhyConfig, BER, FER);
//	delete[] BER;
//	delete[] FER;*/
//	delete[] codeRateArray;
//	delete[] EbN0Array;
//}




// transmission rate for 1-D
//int main() {
//	int numOpt = 1E4;
//	int numDesignSNR = 20;
//	int* nOptArray = new int[20];
//	for (int i = 0; i < 20; i++) { nOptArray[i] = 0; }
//	sPhyConfig* PhyConfig = new sPhyConfig;
//	soptParams** OptParam = new soptParams*[numDesignSNR];
//	for (int i = 0; i < numDesignSNR; i++) { OptParam[i] = new soptParams[numDesignSNR]; }
//	float error = 0.3;
//	int minK = 300;
//	int maxK = 980;
//	int stepK = 20;
//	PhyConfig->nNh = Nh;
//	PhyConfig->nNv = Nv;
//	PhyConfig->nKh = Kh;
//	PhyConfig->nKv = Kv;
//	PhyConfig->nLowRateCodewords = border;
//	PhyConfig->nModType = ModType;
//	PhyConfig->nNt = Nt;
//	PhyConfig->nNr = Nr;
//	PhyConfig->nStream = trans_stream;
//	PhyConfig->bFlagInterleave = flag_interleave;
//	PhyConfig->bFlagInterleaveSpatial = flag_interleave_spatial;
//	PhyConfig->bFlagInterleaveTemporal = flag_interleave_all;
//	PhyConfig->bIrregularCoding = irregular_coding;
//	PhyConfig->powerAllocationScheme = MIMO_SVD_POWER_ALLOCATION_MODE;
//	PhyConfig->nCRCL = CRCL;
//	PhyConfig->minSNR = Min_SNR;
//	PhyConfig->maxSNR = Max_SNR;
//	PhyConfig->stepSNR = SNR_step;
//	PhyConfig->L = L;
//	PhyConfig->nMaxFrameError = errframe;
//	PhyConfig->pKTime = new int[num_time_constructions];
//	PhyConfig->polarConh = polarconh;
//	PhyConfig->polarConV = polarconv;
//	PhyConfig->nTimeConstructions = 2;
//	int minS = 20; // 22,32
//	int maxS = 32;
//	int nOpt = -1;
//	int nDesignSNR = 0;
//	int border_max_1 = 23;
//	int border_max_2 = 18;
//	int K_time_array_1[2] = { 26, 31 };
//	int K_time_array_2[2] = { 16, 26 };
//	float targetErrorRate = 1E-2;
//	float designSNR_min = 7.0;
//	float designSNR_max = 10.0;
//	float designSNR_step = 0.5;
//	/*PhyConfig->minSNR = designSNR;
//	PhyConfig->maxSNR = designSNR;*/
//	ofstream result_file("OptimizationResult.txt",std::ios::app);
//	/*std::cout << std::endl << setw(10) << "Idx Opt" << setw(10) << "Kv" << setw(10) << "Klow" << setw(10) << "Khigh" << setw(10)
//		<< "nLowRate" << setw(10) << "S" << setw(10) << "Rtrans" << setw(10) << "FER" <<std::endl;*/
//	// all possible stream count
//	for (float designSNR = designSNR_min; designSNR <= designSNR_max; designSNR += designSNR_step) {
//		nOpt = -1;
//		for (int nSelStream = maxS; nSelStream >= minS; nSelStream--) {
//			float* FER = new float[1];
//			float* BER = new float[1];
//			int isErrorSatisfy = 0;
//			PhyConfig->nStream = nSelStream;
//			// (32, 26) * (26 + 31);
//			for (int nowK = maxK; nowK >= minK; nowK -= stepK) {
//				PhyConfig->nK = nowK + PhyConfig->nCRCL;
//				PhyConfig->nN = N;
//				PhyConfig->minSNR = SNR2EbN0_1D(designSNR, PhyConfig);
//				PhyConfig->maxSNR = SNR2EbN0_1D(designSNR, PhyConfig);
//				system1D_new(PhyConfig, BER, FER);
//				if (FER[0] < 1.3 * targetErrorRate) {
//					isErrorSatisfy = 1;
//				}
//				nOpt++;
//				if (isErrorSatisfy) {
//					OptParam[nDesignSNR][nOpt].codeRate = (float)nowK / (float)PhyConfig->nN;
//					OptParam[nDesignSNR][nOpt].errorRate = FER[0];
//					OptParam[nDesignSNR][nOpt].nStream = nSelStream;
//					OptParam[nDesignSNR][nOpt].transmissionRate = nSelStream * PhyConfig->nModType * OptParam[nDesignSNR][nOpt].codeRate;
//					std::cout << std::endl <<setw(10)<< designSNR<< setw(10) << nSelStream << setw(10) << OptParam[nDesignSNR][nOpt].codeRate << setw(10) << OptParam[nDesignSNR][nOpt].transmissionRate << setw(10) << FER[0] << std::endl;
//					result_file << std::endl << setw(10) << designSNR << setw(10) << nSelStream << setw(10) << OptParam[nDesignSNR][nOpt].codeRate << setw(10) << OptParam[nDesignSNR][nOpt].transmissionRate << setw(10) << FER[0] << std::endl;
//					if (isErrorSatisfy)
//					{
//						break;
//					}
//				}
//			}
//			delete[] FER;
//			delete[] BER;
//
//		}
//		nOptArray[nDesignSNR] = nOpt + 1;
//		nDesignSNR++;
//	}
//	float rate_max = 0;
//	/*for (int iSNR = 0; iSNR < nDesignSNR; iSNR++) {
//		result_file << std::endl << "SNR=" << designSNR_min + iSNR * designSNR_step << std::endl;
//		for (int iOpt = 0; iOpt < nOptArray[iSNR]-1; iOpt++) {
//			if ((OptParam[iSNR][iOpt].transmissionRate > rate_max) && (OptParam[iSNR][iOpt].errorRate < targetErrorRate)) {
//				result_file << "Optimal Rate and Stream Selection" << std::endl;
//				result_file << std::endl << setw(10) << OptParam[iSNR][iOpt].nStream << setw(10) << OptParam[iSNR][iOpt].codeRate << setw(10) << OptParam[iSNR][iOpt].transmissionRate << std::endl;
//			}
//		}
//	}*/
//	result_file.close();
//}


struct intCombinations {
	int x;
	int y;
	int prod;
};

int findCombinations(int prod_min, int prod_max,
	int x_min, int x_max,
	int y_min, int y_max,
	vector<intCombinations>& comb_vec)
{
	vector<intCombinations> temp_results;
	for (int x = x_min; x <= x_max; ++x) {
		for (int y = y_min; y <= y_max; ++y) {
			long long product = static_cast<long long>(x) * y;
			if (product >= prod_min && product <= prod_max) {
				intCombinations combo;
				combo.x = x;
				combo.y = y;
				combo.prod = static_cast<int>(product);
				temp_results.push_back(combo);
			}
		}
	}
	std::sort(temp_results.begin(), temp_results.end(),
		[](const intCombinations& a, const intCombinations& b) {
			return a.prod > b.prod; // '>'  µœ÷Ωµ–Ú≈≈¡–
		}
	);
	comb_vec = temp_results;
	return comb_vec.size();
}

// transmission rate simulation for 2-D regular-TPC
//int main() {
//	int numOpt = 1E4;
//	int numDesignSNR = 20;
//	int* nOptArray = new int[20];
//	for (int i = 0; i < 20; i++) { nOptArray[i] = 0; }
//	sPhyConfig* PhyConfig = new sPhyConfig;
//	soptParams** OptParam = new soptParams * [numDesignSNR];
//	for (int i = 0; i < numDesignSNR; i++) { OptParam[i] = new soptParams[numDesignSNR]; }
//	float error = 0.3;
//	int minK = 500;
//	int maxK = 906;
//	int stepK = 20;
//	PhyConfig->nNh = Nh;
//	PhyConfig->nNv = Nv;
//	PhyConfig->nKh = Kh;
//	PhyConfig->nKv = Kv;
//	PhyConfig->nLowRateCodewords = border;
//	PhyConfig->nModType = ModType;
//	PhyConfig->nNt = Nt;
//	PhyConfig->nNr = Nr;
//	PhyConfig->nStream = trans_stream;
//	PhyConfig->bFlagInterleave = flag_interleave;
//	PhyConfig->bFlagInterleaveSpatial = flag_interleave_spatial;
//	PhyConfig->bFlagInterleaveTemporal = flag_interleave_all;
//	PhyConfig->bIrregularCoding = irregular_coding;
//	PhyConfig->powerAllocationScheme = MIMO_SVD_POWER_ALLOCATION_MODE;
//	PhyConfig->nCRCL = CRCL;
//	PhyConfig->minSNR = Min_SNR;
//	PhyConfig->maxSNR = Max_SNR;
//	PhyConfig->stepSNR = SNR_step;
//	PhyConfig->L = L;
//	PhyConfig->nMaxFrameError = errframe;
//	PhyConfig->pKTime = new int[num_time_constructions];
//	PhyConfig->polarConh = polarconh;
//	PhyConfig->polarConV = polarconv;
//	PhyConfig->nTimeConstructions = 2;
//	int minS = 20; // 22,32
//	int maxS = 26;
//	int nOpt = -1;
//	int nDesignSNR = 0;
//	int border_max_1 = 23;
//	int border_max_2 = 18;
//	int K_time_array_1[2] = { 26, 31 };
//	int K_time_array_2[2] = { 16, 26 };
//	float targetErrorRate = 1E-2;
//	float designSNR_min = 9.0;
//	float designSNR_max = 30.0;
//	float designSNR_step = 0.5;
//	vector<intCombinations> rateComb;
//	rateComb.resize(2 * PhyConfig->nNh * PhyConfig->nNv);
//	int nRateComb = findCombinations(minK, maxK, 1, PhyConfig->nNh-1, 1, PhyConfig->nNv-1, rateComb);
//	std::cout << setw(10) << "Kv" << setw(10) << "Kh" << setw(10) << "NInfo" << std::endl;
//	for (int i = 0; i < nRateComb; i++) {
//		std::cout << setw(10) << rateComb[i].y << setw(10) << rateComb[i].x << setw(10) << rateComb[i].prod << std::endl;
//	}
//	int RM_set_h[4] = { 6,16,26,31 };
//	int RM_set_v[4] = { 6,16,26,31 };
//	/*PhyConfig->minSNR = designSNR;
//	PhyConfig->maxSNR = designSNR;*/
//	ofstream result_file("OptimizationResult.txt", std::ios::app);
//	result_file << "Regular Polar-TPC" << std::endl;
//	/*std::cout << std::endl << setw(10) << "Idx Opt" << setw(10) << "Kv" << setw(10) << "Klow" << setw(10) << "Khigh" << setw(10)
//		<< "nLowRate" << setw(10) << "S" << setw(10) << "Rtrans" << setw(10) << "FER" <<std::endl;*/
//		// all possible stream count
//	for (float designSNR = designSNR_min; designSNR <= designSNR_max; designSNR += designSNR_step) {
//		nOpt = -1;
//		for (int nSelStream = maxS; nSelStream >= minS; nSelStream--) {
//			float* FER = new float[1];
//			float* BER = new float[1];
//			int isErrorSatisfy = 0;
//			PhyConfig->nStream = nSelStream;
//			// (32, 26) * (26 + 31);
//			for (int iRate = 0; iRate < nRateComb; iRate++) {
//				int nowK = rateComb[iRate].prod;
//				PhyConfig->nKh = rateComb[iRate].x;
//				PhyConfig->nKv = rateComb[iRate].y;
//				PhyConfig->bFlagInterleaveSpatial = true;
//				PhyConfig->bFlagInterleaveTemporal = true;
//				PhyConfig->powerAllocationScheme = "UNIFORM";
//				PhyConfig->polarConh = "5G";
//				PhyConfig->polarConV = "5G";
//				for (int i = 0; i < 4; i++) { if (PhyConfig->nKh == RM_set_h[i]) { PhyConfig->polarConh = "RM"; } }
//				for (int i = 0; i < 4; i++) { if (PhyConfig->nKv == RM_set_v[i]) { PhyConfig->polarConV = "RM"; } }
//				PhyConfig->pKTime[0] = PhyConfig->nKh;
//				PhyConfig->pKTime[1] = PhyConfig->nKh;
//				PhyConfig->minSNR = SNR2EbN0(designSNR, PhyConfig);
//				PhyConfig->maxSNR = SNR2EbN0(designSNR, PhyConfig);
//				system2D_new(PhyConfig, BER, FER);
//				if (FER[0] < 1.3 * targetErrorRate) {
//					isErrorSatisfy = 1;
//				}
//				nOpt++;
//				if (isErrorSatisfy) {
//					OptParam[nDesignSNR][nOpt].codeRate = (float)nowK / ((float)PhyConfig->nNh * (float)PhyConfig->nNv);
//					OptParam[nDesignSNR][nOpt].errorRate = FER[0];
//					OptParam[nDesignSNR][nOpt].nStream = nSelStream;
//					OptParam[nDesignSNR][nOpt].nKh = PhyConfig->nKh;
//					OptParam[nDesignSNR][nOpt].nKv = PhyConfig->nKv;
//					OptParam[nDesignSNR][nOpt].transmissionRate = nSelStream * PhyConfig->nModType * OptParam[nDesignSNR][nOpt].codeRate;
//					std::cout << std::endl << setw(10) << designSNR << setw(10) << nSelStream << setw(10) << OptParam[nDesignSNR][nOpt].codeRate << setw(10) << OptParam[nDesignSNR][nOpt].transmissionRate << setw(10) << FER[0] << " ("<<PhyConfig->nKv<<","<<PhyConfig->nKh<<")"<< std::endl;
//					result_file << std::endl << setw(10) << designSNR << setw(10) << nSelStream << setw(10) << OptParam[nDesignSNR][nOpt].codeRate << setw(10) << OptParam[nDesignSNR][nOpt].transmissionRate << setw(10) << FER[0] << " (" << PhyConfig->nKv << "," << PhyConfig->nKh << ")" << std::endl;;
//					if (isErrorSatisfy)
//					{
//						break;
//					}
//				}
//			}
//			delete[] FER;
//			delete[] BER;
//
//		}
//		nOptArray[nDesignSNR] = nOpt + 1;
//		nDesignSNR++;
//	}
//	float rate_max = 0;
//	/*for (int iSNR = 0; iSNR < nDesignSNR; iSNR++) {
//		result_file << std::endl << "SNR=" << designSNR_min + iSNR * designSNR_step << std::endl;
//		for (int iOpt = 0; iOpt < nOptArray[iSNR]-1; iOpt++) {
//			if ((OptParam[iSNR][iOpt].transmissionRate > rate_max) && (OptParam[iSNR][iOpt].errorRate < targetErrorRate)) {
//				result_file << "Optimal Rate and Stream Selection" << std::endl;
//				result_file << std::endl << setw(10) << OptParam[iSNR][iOpt].nStream << setw(10) << OptParam[iSNR][iOpt].codeRate << setw(10) << OptParam[iSNR][iOpt].transmissionRate << std::endl;
//			}
//		}
//	}*/
//	result_file.close();
//}


// code construction -> 
//int main() {
//	int numOpt = 3E4;
//	int numDesignSNR = 20;
//	float minEbN0 = 4;
//	float maxEbN0 = 20;
//	float stepEbN0 = 0.25;
//	int* nOptArray = new int[20];
//	for (int i = 0; i < 20; i++) { nOptArray[i] = 0; }
//	sPhyConfig* PhyConfig = new sPhyConfig;
//	soptParams** OptParam = new soptParams * [numDesignSNR];
//	for (int i = 0; i < numDesignSNR; i++) { OptParam[i] = new soptParams[numDesignSNR]; }
//	float error = 0.3;
//	int minK = 512;
//	int maxK = 666;
//	int stepK = 20;
//	PhyConfig->nNh = Nh;
//	PhyConfig->nNv = Nv;
//	PhyConfig->nKh = Kh;
//	PhyConfig->nKv = Kv;
//	PhyConfig->nLowRateCodewords = border;
//	PhyConfig->nModType = ModType;
//	PhyConfig->nNt = Nt;
//	PhyConfig->nNr = Nr;
//	PhyConfig->nStream = trans_stream;
//	PhyConfig->bFlagInterleave = flag_interleave;
//	PhyConfig->bFlagInterleaveSpatial = flag_interleave_spatial;
//	PhyConfig->bFlagInterleaveTemporal = flag_interleave_all;
//	PhyConfig->bIrregularCoding = irregular_coding;
//	PhyConfig->powerAllocationScheme = MIMO_SVD_POWER_ALLOCATION_MODE;
//	PhyConfig->nCRCL = CRCL;
//	PhyConfig->minSNR = Min_SNR;
//	PhyConfig->maxSNR = Max_SNR;
//	PhyConfig->stepSNR = SNR_step;
//	PhyConfig->L = L;
//	PhyConfig->nMaxFrameError = errframe;
//	PhyConfig->pKTime = new int[num_time_constructions];
//	PhyConfig->polarConh = polarconh;
//	PhyConfig->polarConV = polarconv;
//	PhyConfig->nTimeConstructions = 2;
//	int minS = 20; // 22,32
//	int maxS = 26;
//	int nOpt = -1;
//	int nDesignSNR = 0;
//	int border_max_1 = 23;
//	int border_max_2 = 18;
//	int K_time_array_1[2] = { 26, 31 };
//	int K_time_array_2[2] = { 16, 26 };
//	float targetErrorRate = 1E-3;
//	float designSNR_min = 9.0;
//	float designSNR_max = 30.0;
//	float designSNR_step = 0.5;
//	vector<intCombinations> rateComb;
//	rateComb.resize(2 * PhyConfig->nNh * PhyConfig->nNv);
//	int nRateComb = findCombinations(minK, maxK, 1, PhyConfig->nNh - 1, 1, PhyConfig->nNv - 1, rateComb);
//	std::cout << setw(10) << "Kv" << setw(10) << "Kh" << setw(10) << "NInfo" << std::endl;
//	for (int i = 0; i < nRateComb; i++) {
//		std::cout << setw(10) << rateComb[i].y << setw(10) << rateComb[i].x << setw(10) << rateComb[i].prod << std::endl;
//	}
//	std::cout << "Totally " << nRateComb << " code rates to simulate" << std::endl;
//	/*int RM_set_h[4] = { 6,16,26,31 };
//	int RM_set_v[4] = { 6,16,26,31 };*/
//	int RM_set_h[4] = { 6,16,26,31 };
//	int RM_set_v[4] = { 6,16,26,31 };
//	/*PhyConfig->minSNR = designSNR;
//	PhyConfig->maxSNR = designSNR;*/
//	ofstream result_file("AllRateRegularResult1024.txt", std::ios::app);
//	// result_file << "Regular Polar-TPC" << std::endl;
//	float* FER = new float[1];
//	float* BER = new float[1];
//	int isErrorSatisfy = 0;
//	PhyConfig->nStream = Nt;
//	// float nowEbN0 = minEbN0; nowEbN0 <= maxEbN0; nowEbN0++
//	// (32, 26) * (26 + 31);
//	for (int iRate = 0; iRate < nRateComb; iRate++) {
//		isErrorSatisfy = 0;
//		for (float nowEbN0 = minEbN0; nowEbN0 <= maxEbN0; nowEbN0+=stepEbN0) {
//			int nowK = rateComb[iRate].prod;
//			PhyConfig->nKh = rateComb[iRate].x;
//			PhyConfig->nKv = rateComb[iRate].y;
//			if (abs(PhyConfig->nKh-PhyConfig->nKv)>10||(PhyConfig->nKh>26)||(PhyConfig->nKv>26)) { break; }
//			PhyConfig->bFlagInterleaveSpatial = true;
//			PhyConfig->bFlagInterleaveTemporal = true;
//			PhyConfig->powerAllocationScheme = "UNIFORM";
//			PhyConfig->polarConh = "5G";
//			PhyConfig->polarConV = "5G";
//			for (int i = 0; i < 4; i++) { if (PhyConfig->nKh == RM_set_h[i]) { PhyConfig->polarConh = "RM"; } }
//			for (int i = 0; i < 4; i++) { if (PhyConfig->nKv == RM_set_v[i]) { PhyConfig->polarConV = "RM"; } }
//			PhyConfig->pKTime[0] = PhyConfig->nKh;
//			PhyConfig->pKTime[1] = PhyConfig->nKh;
//			PhyConfig->minSNR = nowEbN0;
//			PhyConfig->maxSNR = nowEbN0;
//			system2D_new(PhyConfig, BER, FER);
//			if (FER[0] < 1.3 * targetErrorRate) {
//				isErrorSatisfy = 1;
//			}
//			if (isErrorSatisfy) {
//				nOpt++;
//				OptParam[nDesignSNR][nOpt].codeRate = (float)nowK / ((float)PhyConfig->nNh * (float)PhyConfig->nNv);
//				OptParam[nDesignSNR][nOpt].errorRate = FER[0];
//				OptParam[nDesignSNR][nOpt].nStream = PhyConfig->nStream;
//				OptParam[nDesignSNR][nOpt].nKh = PhyConfig->nKh;
//				OptParam[nDesignSNR][nOpt].nKv = PhyConfig->nKv;
//				// OptParam[nDesignSNR][nOpt].transmissionRate = nSelStream * PhyConfig->nModType * OptParam[nDesignSNR][nOpt].codeRate;
//				// std::cout << std::endl << setw(10) << designSNR << setw(10) << nSelStream << setw(10) << OptParam[nDesignSNR][nOpt].codeRate << setw(10) << OptParam[nDesignSNR][nOpt].transmissionRate << setw(10) << FER[0] << " (" << PhyConfig->nKv << "," << PhyConfig->nKh << ")" << std::endl;
//				result_file << setw(10) << nowEbN0 << setw(10) << PhyConfig->nKv <<setw(10) <<PhyConfig->nKh<< std::endl;
//				std::cout << setw(10) << nowEbN0 << setw(10) << PhyConfig->nKv << setw(10) << PhyConfig->nKh << std::endl;
//				if (isErrorSatisfy)
//				{
//					break;
//				}
//			}
//		}
//		
//	}
//	delete[] FER;
//	delete[] BER;
//	float rate_max = 0;
//	result_file.close();
//}
//


// Simulation: BLER v.s. Alpha (scaling factor for power allocation)
//int main() {
//	int numOpt = 3E4;
//	int numDesignSNR = 20;
//	float minAlpha = 0.0;
//	float maxAlpha = 1.0;
//	float stepAlpha = 0.05;
//	float minEbN0 = 8.20;
//	float maxEbN0 = 8.20;
//	float stepEbN0 = 0.25;
//	int* nOptArray = new int[20];
//	for (int i = 0; i < 20; i++) { nOptArray[i] = 0; }
//	sPhyConfig* PhyConfig = new sPhyConfig;
//	soptParams** OptParam = new soptParams * [numDesignSNR];
//	for (int i = 0; i < numDesignSNR; i++) { OptParam[i] = new soptParams[numDesignSNR]; }
//	float error = 0.3;
//	int minK = 1;
//	int maxK = 351;
//	int stepK = 20;
//	PhyConfig->nNh = Nh;
//	PhyConfig->nNv = Nv;
//	PhyConfig->nKh = Kh;
//	PhyConfig->nKv = Kv;
//	PhyConfig->nLowRateCodewords = border;
//	PhyConfig->nModType = ModType;
//	PhyConfig->nNt = Nt;
//	PhyConfig->nNr = Nr;
//	PhyConfig->nStream = trans_stream;
//	PhyConfig->bFlagInterleave = flag_interleave;
//	PhyConfig->bFlagInterleaveSpatial = flag_interleave_spatial;
//	PhyConfig->bFlagInterleaveTemporal = flag_interleave_all;
//	PhyConfig->bIrregularCoding = irregular_coding;
//	PhyConfig->powerAllocationScheme = MIMO_SVD_POWER_ALLOCATION_MODE;
//	PhyConfig->nCRCL = CRCL;
//	PhyConfig->minSNR = Min_SNR;
//	PhyConfig->maxSNR = Max_SNR;
//	PhyConfig->stepSNR = SNR_step;
//	PhyConfig->L = L;
//	PhyConfig->nMaxFrameError = errframe;
//	PhyConfig->pKTime = new int[num_time_constructions];
//	PhyConfig->polarConh = polarconh;
//	PhyConfig->polarConV = polarconv;
//	PhyConfig->nTimeConstructions = 2;
//	int minS = 20; // 22,32
//	int maxS = 26;
//	int nOpt = -1;
//	int nDesignSNR = 0;
//	int border_max_1 = 23;
//	int border_max_2 = 18;
//	int K_time_array_1[2] = { 26, 31 };
//	int K_time_array_2[2] = { 16, 26 };
//	float targetErrorRate = 1E-2;
//	float designSNR_min = 9.0;
//	float designSNR_max = 30.0;
//	float designSNR_step = 0.5;
//	vector<intCombinations> rateComb;
//	rateComb.resize(2 * PhyConfig->nNh * PhyConfig->nNv);
//	int nRateComb = findCombinations(minK, maxK, 1, PhyConfig->nNh - 1, 1, PhyConfig->nNv - 1, rateComb);
//	std::cout << setw(10) << "Kv" << setw(10) << "Kh" << setw(10) << "NInfo" << std::endl;
//	for (int i = 0; i < nRateComb; i++) {
//		std::cout << setw(10) << rateComb[i].y << setw(10) << rateComb[i].x << setw(10) << rateComb[i].prod << std::endl;
//	}
//	int RM_set_h[4] = { 6,16,26,31 };
//	int RM_set_v[4] = { 6,16,26,31 };
//	/*PhyConfig->minSNR = designSNR;
//	PhyConfig->maxSNR = designSNR;*/
//	ofstream result_file("AllPowerAllocResult.txt", std::ios::app);
//	// result_file << "Regular Polar-TPC" << std::endl;
//	float* FER = new float[1];
//	float* BER = new float[1];
//	int isErrorSatisfy = 0;
//	PhyConfig->nStream = 31;
//	PhyConfig->nLowRateCodewords = border;
//	// float nowEbN0 = minEbN0; nowEbN0 <= maxEbN0; nowEbN0++
//	// (32, 26) * (26 + 31);
//	result_file << K_time_array[0] << "   " << K_time_array[1] << "   " << PhyConfig->nLowRateCodewords << std::endl;
//	std::cout << K_time_array[0] << "   " << K_time_array[1] << "   " << PhyConfig->nLowRateCodewords << std::endl;
//	for (float nowAlpha = minAlpha; nowAlpha <= maxAlpha; nowAlpha+=stepAlpha) {
//		isErrorSatisfy = 0;
//		for (float nowEbN0 = minEbN0; nowEbN0 <= maxEbN0; nowEbN0 += stepEbN0) {
//			PhyConfig->nKh = Kh;
//			PhyConfig->nKv = Kv;
//			PhyConfig->bFlagInterleaveSpatial = false;
//			PhyConfig->bFlagInterleaveTemporal = true;
//			PhyConfig->powerAllocationScheme = "SEPARATE";
//			PhyConfig->polarConh = "RM";
//			PhyConfig->polarConV = "GA";
//			PhyConfig->powerAllocationRatio = nowAlpha;
//			/*for (int i = 0; i < 4; i++) { if (PhyConfig->nKh == RM_set_h[i]) { PhyConfig->polarConh = "RM"; } }
//			for (int i = 0; i < 4; i++) { if (PhyConfig->nKv == RM_set_v[i]) { PhyConfig->polarConV = "RM"; } }*/
//			/*PhyConfig->pKTime[0] = PhyConfig->nKh;
//			PhyConfig->pKTime[1] = PhyConfig->nKh;*/
//			for (int i = 0; i < 2; i++) { PhyConfig->pKTime[i] = K_time_array[i]; }
//			PhyConfig->minSNR = nowEbN0;
//			PhyConfig->maxSNR = nowEbN0;
//			system2D_new(PhyConfig, BER, FER);
//			result_file << nowAlpha << "    " << FER[0] << std::endl;
//			std::cout << std::endl  << nowAlpha << "    " << FER[0] << std::endl;
//		}
//
//	}
//	delete[] FER;
//	delete[] BER;
//	float rate_max = 0;
//	result_file.close();
//}

