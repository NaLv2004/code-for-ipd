#pragma once
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
#include "system_new.h"

//using pInt = std::unique_ptr<int[]>;
//using pFloat = std::unique_ptr<float[]>;
//using pChar = std::unique_ptr<unsigned char[]>;
//using pCfloat = std::unique_ptr<std::complex<float>[]>;
//// using pVector = std::unique_ptr<std::vector[]>

Eigen::MatrixXf Fn(int N)
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

// 辅助函数：求最大公约数 (使用辗转相除法/欧几里得算法)
int gcd(int a, int b) {
	while (b != 0) {
		int temp = a % b;
		a = b;
		b = temp;
	}
	return a;
}

// 主函数：求最大公倍数
int lcm(int a, int b) {
	if (a == 0 || b == 0) return 0; // 任何数与 0 的最大公倍数设为 0
	return std::abs(a * b) / gcd(a, b);
}

void polar_code_construction_new(sPhyConfig* PhyConfig, float noiseVariance, pInt& A_H, pInt& Ac_H, pInt& A_Ac_H, pInt& A_V, pInt& Ac_V, pInt& A_Ac_V, pInt& kTimeAll, pInt& isLowRateCodeword,
	int** GGh, int** GGv) {
	vector<int> temph;
	vector<int> tempv;
	vector<int> tempv_ordered;
	vector<int>* construction_time_all = new vector<int>[PhyConfig->nNv];
	for (int i = 0; i < PhyConfig->nNv; i++) { construction_time_all[i].resize(PhyConfig->nNh); }
	temph.resize(PhyConfig->nNh);
	tempv.resize(PhyConfig->nNv);
	tempv_ordered.resize(PhyConfig->nNv);
	polar_codeconstruction(PhyConfig->polarConh, PhyConfig->nNh, PhyConfig->nKh, GGh, (float)sqrt(noiseVariance), temph);
	polar_codeconstruction(PhyConfig->polarConV, PhyConfig->nNv, PhyConfig->nKv, GGv, (float)sqrt(noiseVariance), tempv);
	// store the ordered reliability of vertical frozen bits
	for (int i = 0; i < PhyConfig->nNv; i++) { tempv_ordered[i] = tempv[i]; }
	// sort the first K elements of tempv and temph
	for (int i = 0; i < PhyConfig->nKh - 1; i++)
		for (int j = i + 1; j < PhyConfig->nKh; j++)
			if (temph[i] > temph[j])
			{
				int ttt = temph[i];
				temph[i] = temph[j];
				temph[j] = ttt;
			}
	for (int i = 0; i < PhyConfig->nKv - 1; i++)
		for (int j = i + 1; j < PhyConfig->nKv; j++)
			if (tempv[i] > tempv[j])
			{
				int ttt = tempv[i];
				tempv[i] = tempv[j];
				tempv[j] = ttt;
			}
	// fill in the vertical information and frozen set
	for (int i = 0; i < PhyConfig->nKv; i++) // Information Set
	{
		A_V[i] = tempv[i];
		A_Ac_V[A_V[i]] = 1;
	}
	//	getchar();
	for (int i = 0; i < PhyConfig->nNv - PhyConfig->nKv; i++) // Frozen Bit
	{
		Ac_V[i] = tempv[PhyConfig->nKv + i];
		A_Ac_V[Ac_V[i]] = 0;
	}
	// horizontal construction
	for (int k = 0; k < PhyConfig->nNv; k++)
	{
		bool bfind = false;
		int channel_index = 0;
		for (int i_find = 0; i_find < PhyConfig->nKv; i_find++) { if (tempv_ordered[i_find] == k) { bfind = true;  channel_index = i_find; break; } }
		int k_time = PhyConfig->pKTime[1];
	    if (bfind && ((PhyConfig->nKv - 1 - channel_index) < PhyConfig->nLowRateCodewords)) { k_time = PhyConfig->pKTime[0]; }  //prev: Kv-1
		kTimeAll[k] = k_time;
		// std::cout << noiseVariance << endl;
		polar_codeconstruction(PhyConfig->polarConh, PhyConfig->nNh, k_time, GGh, (float)sqrt(noiseVariance), construction_time_all[k]);
		if (k_time == PhyConfig->pKTime[0]) { isLowRateCodeword[k] = 1; } // use long codeword
		else if (k_time == PhyConfig->pKTime[1]) { isLowRateCodeword[k] = 0; } // use short codeword
		for (int i = 0; i < k_time - 1; i++)
			for (int j = i + 1; j < k_time; j++)
				if (construction_time_all[k][i] > construction_time_all[k][j])
				{
					int ttt = construction_time_all[k][i];
					construction_time_all[k][i] = construction_time_all[k][j];
					construction_time_all[k][j] = ttt;
				}
		// for (int i = 0; i < PhyConfig->nNh; i++) { std::cout << construction_time_all[k][i] << "  "; }
		// Decide the position of information bits and frozen bits
		for (int i = 0; i < k_time; i++) // Information Set
		{
			A_H[k*PhyConfig->nNh+i] = construction_time_all[k][i];
			A_Ac_H[k * PhyConfig->nNh+A_H[k * PhyConfig->nNh+i]] = 1;
		}
		//	getchar();
		for (int i = 0; i < PhyConfig->nNh - k_time; i++) // Frozen Bit
		{
			Ac_H[k * PhyConfig->nNh+i] = construction_time_all[k][k_time + i];
			A_Ac_H[k * PhyConfig->nNh+Ac_H[k * PhyConfig->nNh+i]] = 0;
		}
	}
	delete[] construction_time_all;
}


void polar_code_construction_new_1D(sPhyConfig* PhyConfig, float noiseVariance, pInt& A, pInt& A_Ac, pInt& Ac, int** GG) {
	vector<int> construction;
	construction.resize(PhyConfig->nN);
	polar_codeconstruction(PhyConfig->polarCon, PhyConfig->nN, PhyConfig->nK, GG, (float)sqrt(noiseVariance), construction);
	for (int i = 0; i < PhyConfig->nK - 1; i++)
		for (int j = i + 1; j < PhyConfig->nK; j++)
			if (construction[i] > construction[j])
			{
				int ttt = construction[i];
				construction[i] = construction[j];
				construction[j] = ttt;
			}
	for (int i = 0; i < PhyConfig->nK; i++) {
		A[i] = construction[i];
		A_Ac[A[i]] = 1;
	}
	// for (int i = 0; i < PhyConfig->nK; i++) { std::cout << A[i] << std::endl; }
	for (int i = 0; i < PhyConfig->nN - PhyConfig->nK; i++) {
		Ac[i] = construction[PhyConfig->nK + i];
		A_Ac[Ac[i]] = 0;
	}
}

void generateLargeInfoBlock(sPhyConfig* PhyConfig, pInt& idxComponentCodewordInBlock, pInt& idxLowRateCodeword, pInt& idxHighRateCodeword, pInt& isLowRateCodeword, pInt& uInfoBlock, pInt& originalInfoBlock, pInt& xEncodedBlock,
pInt& A_H, pInt& Ac_H, pInt& A_Ac_H, pInt& A_V, pInt& A_Ac_V, pInt& Ac_V, pInt& kTimeAll, pInt& nBadChannelsMIMOBlock, int nComponentInBlock, int nWholeCodewordInBlock, int nLenBlock, int nLenInfoBlock, int nInfo, int codeLength, int& nBadChannels,
	std::mt19937& gen, std::uniform_int_distribution<> distrib
) {
	pInt idxBadChannelInBlock(new int[nComponentInBlock]); // channels allocated to low rate codewords
	pInt idxGoodChannelInBlock(new int[nComponentInBlock]); //  channels allocated to high rate codewords
// # pragma omp parallel for
	for (int i = 0; i < nWholeCodewordInBlock; i++) {
		pInt xHorizontalEncode(new int[PhyConfig->nNh* PhyConfig->nNv]);
		pInt verticalEncodeTmp(new int[PhyConfig->nNv]);
		pInt verticalEncodeTmpMid(new int[PhyConfig->nNv]);
		for (int j = 0; j < PhyConfig->nNh * PhyConfig->nNv; j++) { xHorizontalEncode[j] = 0; }
		for (int j = 0; j < PhyConfig->nNv; j++) { verticalEncodeTmp[j] = 0; verticalEncodeTmpMid[j] = 0; }
		// Generate info bits
		for (int j = 0; j < nInfo; j++) { originalInfoBlock[i * codeLength  + j] = distrib(gen) % 2; }
		int i_bit_ori = 0;
		for (int j = 0; j < PhyConfig->nKv; j++) {
			for (int i_bit = 0; i_bit < kTimeAll[A_V[j]]; i_bit++) { uInfoBlock[i * codeLength + PhyConfig->nNh * A_V[j] + A_H[A_V[j]*PhyConfig->nNh+i_bit]] = originalInfoBlock[i * codeLength + i_bit_ori]; 
			// std::cout << i * codeLength + PhyConfig->nNh * A_V[j] + A_H[A_V[j] * PhyConfig->nNh + i_bit] << endl;
			i_bit_ori++; 
			}
			// std::cout << endl;
		}
		// Encoding
		// horizontal Encoding
		for (int j = 0; j < PhyConfig->nKv; j++) {
			int i_row = A_V[j];
			PolarEncode_xor(xHorizontalEncode.get() + i_row * PhyConfig->nNh, uInfoBlock.get() + i * codeLength + i_row * PhyConfig->nNh, PhyConfig->nNh);
		}
		// Vertical Encoding 
		// int ttttt = 0;
		for (int j = 0; j < PhyConfig->nNh; j++) {
			for (int i_row = 0; i_row < PhyConfig->nNv; i_row++) { verticalEncodeTmp[i_row] = xHorizontalEncode[i_row * PhyConfig->nNh + j]; }
			PolarEncode_xor(verticalEncodeTmpMid.get(), verticalEncodeTmp.get(), PhyConfig->nNv);
			for (int i_frozen = 0; i_frozen < PhyConfig->nNv - PhyConfig->nKv; i_frozen++) { verticalEncodeTmpMid[Ac_V[i_frozen]] = 0; }
			PolarEncode_xor(verticalEncodeTmp.get(), verticalEncodeTmpMid.get(), PhyConfig->nNv);
			for (int i_row = 0; i_row < PhyConfig->nNv; i_row++) { xHorizontalEncode[i_row * PhyConfig->nNh + j] = verticalEncodeTmp[i_row]; }
		}
		// fill to long block
		for (int j = 0; j < codeLength; j++) { xEncodedBlock[i * codeLength + j] = xHorizontalEncode[j];/* std::cout << i * codeLength + j << endl;*/ }
	}
	int nGoodChannels = PhyConfig->nStream - nBadChannels;
	int nMIMOInBlock = round(nComponentInBlock / PhyConfig->nStream);
	int nCodewordtoBadChannel = nMIMOInBlock * nBadChannels;
	int nCodewordtoGoodChannel = nMIMOInBlock * (PhyConfig->nStream - nBadChannels);
	int nHighRateCodewordInBlock = (PhyConfig->nNv - PhyConfig->nLowRateCodewords) * nWholeCodewordInBlock;
	int nLowRateCodewordInBlock = (PhyConfig->nLowRateCodewords) * nWholeCodewordInBlock;
	int nDiff = nLowRateCodewordInBlock - nCodewordtoBadChannel;
	for (int i = 0; i < nMIMOInBlock; i++) { nBadChannelsMIMOBlock[i] = nBadChannels; }
	if (PhyConfig->pKTime[0] != PhyConfig->pKTime[1]) {
		int iComponentInBlock = 0;
		for (int i = 0; i < nMIMOInBlock; i++) {
			if (i < nDiff) { nBadChannelsMIMOBlock[i] = nBadChannels + 1; for (int j = 0; j < nBadChannels + 1; j++) { idxBadChannelInBlock[iComponentInBlock] = i * PhyConfig->nStream + PhyConfig->nStream - 1 - j; iComponentInBlock++; } }
			else { nBadChannelsMIMOBlock[i] = nBadChannels; for (int j = 0; j < nBadChannels; j++) { idxBadChannelInBlock[iComponentInBlock] = i * PhyConfig->nStream + PhyConfig->nStream - 1 - j; iComponentInBlock++; } }
		}
		iComponentInBlock = 0;
		for (int i = 0; i < nMIMOInBlock; i++) {
			if (i < nDiff) { for (int j = 0; j < nGoodChannels - 1; j++) { idxGoodChannelInBlock[iComponentInBlock] = i * PhyConfig->nStream + j; iComponentInBlock++; } }
			else { for (int j = 0; j < nGoodChannels; j++) { idxGoodChannelInBlock[iComponentInBlock] = i * PhyConfig->nStream + j; iComponentInBlock++; } }
		}
		int iLow = 0;
		int iHigh = 0;
		for (int i = 0; i < nWholeCodewordInBlock; i++) {
			for (int j = 0; j < PhyConfig->nNv; j++) {
				if (isLowRateCodeword[j] == 1) { idxComponentCodewordInBlock[i * PhyConfig->nNv + j] = idxBadChannelInBlock[iLow]; iLow++; }
				else { idxComponentCodewordInBlock[i * PhyConfig->nNv + j] = idxGoodChannelInBlock[iHigh]; iHigh++; }
			}
		}
		// test
		/*std::cout << std::endl;
		for (int i = 0; i < 32; i++) std::cout << idxComponentCodewordInBlock[i]<< "  ";
		std::cout << std::endl;*/
	}
	else {
		for (int i = 0; i < nComponentInBlock; i++) { idxComponentCodewordInBlock[i] = i; }
	}
}

void randIntrlvWholeBlock(sPhyConfig* PhyConfig, pInt& xEncodedBlock, vector<vector<int>>& intrlvPatternTime , vector<vector<int>>& intrlvPatternSpace,  vector<vector<int>> intrlvPatternSpaceSep, int nWholeCodewordInBlock, int codeLength) {
	if (PhyConfig->bFlagInterleaveSpatial) {
		for (int i = 0; i < nWholeCodewordInBlock; i++) {
			int* intrlv_tmp = new int[PhyConfig->nNv];
			for (int j = 0; j < PhyConfig->nNh; j++) {
				int offset = i * codeLength;
				for (int k = 0; k < PhyConfig->nNv; k++) { intrlv_tmp[k] = xEncodedBlock[offset + k * PhyConfig->nNh + j]; }
				if (PhyConfig->pKTime[0] == PhyConfig->pKTime[1]) { randIntrlv(intrlvPatternSpace[j], intrlv_tmp, PhyConfig->nNv); }
				else { randIntrlv(intrlvPatternSpaceSep[j], intrlv_tmp, PhyConfig->nNv); }
				for (int k = 0; k < PhyConfig->nNv; k++) { xEncodedBlock[offset + k * PhyConfig->nNh + j] = intrlv_tmp[k]; }
			}
			delete[] intrlv_tmp;
		}
	}
# pragma omp parallel for 
	for (int i = 0; i < nWholeCodewordInBlock; i++) {
		for (int j = 0; j < PhyConfig->nNv; j++) {
			randIntrlv(intrlvPatternTime[j], xEncodedBlock.get() + i * codeLength + j * PhyConfig->nNh, PhyConfig->nNh);
		}
}
}


void randIntrlvWholeBlock1D(sPhyConfig* PhyConfig, pInt& xEncodedBlock, vector<int> intrlvPattern1D ,int nWholeCodewordInBlock, int codeLength) {
	for (int i = 0; i < nWholeCodewordInBlock; i++) {
		randIntrlv(intrlvPattern1D, xEncodedBlock.get() + i * codeLength, PhyConfig->nN);
	}
}

void codewordToAntennaMapping(sPhyConfig* PhyConfig, pInt& idxComponentCodewordInBlock, pInt& xTransmittedBlock , pInt& xEncodedBlock, int nComponentCodewordInBlock) {
	for (int i = 0; i < nComponentCodewordInBlock; i++) {
		for (int j = 0; j < PhyConfig->nNh; j++) {
			xTransmittedBlock[idxComponentCodewordInBlock[i] * PhyConfig->nNh + j] = xEncodedBlock[i * PhyConfig->nNh + j];
			// std::cout << idxComponentCodewordInBlock[i] * PhyConfig->nNh + j << std::endl;
			// std::cout << idxComponentCodewordInBlock[i] * PhyConfig->nNh + j << "   " << i * PhyConfig->nNh + j << std::endl;
		}
	}
}

void modulationBlock(sPhyConfig* PhyConfig, pInt& xTransmittedBlock, pCfloat& xModulationBlock, int nComponentCodewordInBlock) {
	for (int i = 0; i < nComponentCodewordInBlock; i++) {
		// std::cout << i << std::endl;
		modulation(xTransmittedBlock.get() + i * PhyConfig->nNh, xModulationBlock.get() + i * (PhyConfig->nNh / PhyConfig->nModType), PhyConfig->nModType, PhyConfig->nNh / PhyConfig->nModType);
	}
}

void generateChannelsBlock(sPhyConfig* PhyConfig, std::vector<MatrixXcf>& channelMatrixArray, int nMIMOInBlock, std::vector<MatrixXcf>& channelMatrixArrayPlatform, std::mt19937& gen, std::uniform_int_distribution<> distrib) {
// # pragma omp parallel for
	if (PhyConfig->use_platform_channel || PhyConfig->use_elaa_channel)
	{
		for (int i = 0; i < nMIMOInBlock; i++) {
			int idx = distrib(gen) % PhyConfig->n_frequency_samples;
			channelMatrixArray[i] = channelMatrixArrayPlatform[idx];
		}
    }
	else {
		for (int i = 0; i < nMIMOInBlock; i++) {
			channelMatrixArray[i] = generateChannelMatrix(PhyConfig->nNr, PhyConfig->nNt);
		}
	}
}

void generateChannelsBlock(sPhyConfig* PhyConfig, std::vector<MatrixXcf>& channelMatrixArray, int nMIMOInBlock) {
	for (int i = 0; i < nMIMOInBlock; i++) {
		channelMatrixArray[i] = generateChannelMatrix(PhyConfig->nNr, PhyConfig->nNt);
	}
}

void MMSEDetectionCfloat(pCfloat& y_recv, Eigen::MatrixXcf H, double noiseVariance, sPhyConfig* PhyConfig) {
	Eigen::MatrixXcf y_recv_vec(PhyConfig->nNr, 1);
	Eigen::MatrixXcf identityMatrix(PhyConfig->nNt, PhyConfig->nNt);
	identityMatrix.setIdentity();
	for (int i = 0; i < PhyConfig->nNr; i++) { y_recv_vec(i, 0) = y_recv[i]; }
	Eigen::MatrixXcf y_demod = (H.adjoint() * H+noiseVariance*identityMatrix).inverse() * H.adjoint()* y_recv_vec;
	for (int i = 0; i < PhyConfig->nNt; i++) { y_recv[i] = y_demod(i, 0); }
}

void transmissionMIMOBlock(sPhyConfig* PhyConfig, std::vector<MatrixXcf>& channelMatrixArray, int nMIMOInBlock, int nSymbols, int lenMIMOBlock, pCfloat& xTransmittedBlock, pFloat& llrBlock , pInt& nBadChannelsMIMOBlock, float noiseVariance) {
	// llrBlock is of nMIMOInBlock * nComponentCodewordInMIMO * Nh
	int lenMIMOBlockSym = lenMIMOBlock / PhyConfig->nModType;
	if (PhyConfig->powerAllocationScheme == "LONG_TERM") {
		vector<Eigen::MatrixXcf> vecU;
		vector<Eigen::MatrixXcf> vecV;
		vector<Eigen::JacobiSVD<Eigen::MatrixXcf>> vecSV;
		vecU.resize(nMIMOInBlock); vecV.resize(nMIMOInBlock); vecSV.resize(nMIMOInBlock);
		pFloat antennaPowerArray(new float[PhyConfig->nNt * nMIMOInBlock]);
		pFloat svArray2D(new float[PhyConfig->nNt * nMIMOInBlock]);
		float sumSVTmp = 0.0;
		float r = 1.0;
		for (int i = 0; i < PhyConfig->nNt * nMIMOInBlock; i++) { antennaPowerArray[i] = 0.0; }
# pragma omp parallel for
		for (int i_block = 0; i_block < nMIMOInBlock; i_block++) {
			JacobiSVD<MatrixXcf> svd(channelMatrixArray[i_block], ComputeFullU | ComputeFullV);
			vecU[i_block] = svd.matrixU();
			vecV[i_block] = svd.matrixV();
			auto S = svd.singularValues();
			for (int i = 0; i < PhyConfig->nNt; i++) { 
				if (i >= PhyConfig->nStream) { r = 0.0; }
				else {
					if (i > PhyConfig->nStream - nBadChannelsMIMOBlock[i_block] - 1) { r = PhyConfig->powerAllocationRatio; }
					else { r = 1.0; }
				}
				svArray2D[i_block * PhyConfig->nNt + i] = S[i]; sumSVTmp += r * 1.0 / (S[i] * S[i]);
			}
		}
		float sum_power = 0.0;
		for (int i_block = 0; i_block < nMIMOInBlock; i_block++) {
			for (int i = 0; i < PhyConfig->nNt; i++) {
				if (i >= PhyConfig->nStream) { r = 0.0; }
				else {
					if (i > PhyConfig->nStream - nBadChannelsMIMOBlock[i_block] - 1) { r = PhyConfig->powerAllocationRatio; }
					else { r = 1.0; }
				}
				antennaPowerArray[i_block * PhyConfig->nNt + i] = PhyConfig->nNt*(float)nMIMOInBlock*(1.0 / (svArray2D[i_block * PhyConfig->nNt+i] * svArray2D[i_block * PhyConfig->nNt+i]) * r) / sumSVTmp;
				sum_power += antennaPowerArray[i_block * PhyConfig->nNt + i];
			}
		}
		// std::cout << sum_power << std::endl;
		//std::cout << nMIMOInBlock << std::endl;
# pragma omp parallel for
		for (int i_block = 0; i_block < nMIMOInBlock; i_block++) {
			MatrixXcf U = vecU[i_block];
			MatrixXcf V = vecV[i_block];
			pFloat S(new float[PhyConfig->nNt]);
			for (int i = 0; i < PhyConfig->nNt; i++) {
				S[i] = svArray2D[i_block * PhyConfig->nNt + i];
			}
			double noiseVariance_tmp = noiseVariance;
			pCfloat x_trans_tmp(new std::complex<float>[PhyConfig->nNt]);
			// pCfloat x_trans_ori(new std::complex<float>[PhyConfig->nNt]); // for test
			pCfloat x_precoded_tmp(new std::complex<float>[PhyConfig->nNt]);
			pCfloat y_recv_tmp(new std::complex<float>[PhyConfig->nNr]);
			pCfloat y_preprocessed(new std::complex<float>[PhyConfig->nNr]);
			pFloat singular_values_array(new float[PhyConfig->nNt]);
			for (int i = 0; i < PhyConfig->nNt; i++) { x_trans_tmp[i] = 0.0; }
			for (int i = 0; i < PhyConfig->nNr; i++) { y_recv_tmp[i] = 0.0; }
			pFloat antennaPower(new float[PhyConfig->nNt]);
			for (int i = 0; i < PhyConfig->nStream; i++) { antennaPower[i] = antennaPowerArray[i_block*PhyConfig->nNt+i]; }
			for (int i = 0; i < PhyConfig->nStream; i++) { singular_values_array[i] = S[i]; }
			for (int i = PhyConfig->nStream; i < PhyConfig->nNt; i++) { antennaPower[i] = 0.0; }
			for (int i = 0; i < PhyConfig->nNt; i++) { S[i] = S[i] * sqrt(antennaPower[i]); }
			for (int i_sym_col = 0; i_sym_col < nSymbols; i_sym_col++) {
				for (int i = 0; i < PhyConfig->nStream; i++) {
					x_trans_tmp[i] = sqrt(antennaPower[i]) * xTransmittedBlock[i_block * lenMIMOBlockSym + i * nSymbols + i_sym_col];
					// x_trans_ori[i] =  xTransmittedBlock[i_block * lenMIMOBlockSym + i * nSymbols + i_sym_col];
				}
				// test: transmitted symbols 
				precodeSVD(PhyConfig->nNr, PhyConfig->nNt, x_trans_tmp.get(), x_precoded_tmp.get(), V);
				AWGN(PhyConfig->nNr, PhyConfig->nNt, x_precoded_tmp.get(), y_recv_tmp.get(), channelMatrixArray[i_block], noiseVariance_tmp);
				precodeSVD_preprocessing(PhyConfig->nNr, PhyConfig->nNt, y_recv_tmp.get(), y_preprocessed.get(), U);
				for (int i = 0; i < PhyConfig->nStream; i++) { y_recv_tmp[i] = S[i] * y_preprocessed[i] / (S[i] * S[i] + noiseVariance); }
				// for (int i = 0; i < PhyConfig->nStream; i++) { std::cout << x_trans_ori[i] << "  " << y_recv_tmp[i] << endl; }
				for (int i = 0; i < PhyConfig->nStream; i++) {
					for (int j = 0; j < PhyConfig->nModType; j++) {
						// int col = i_sym_col*PhyConfig->nModType+j
						// llrBlock[lenMIMOBlock*i_Block + i*PhyConfig->nNh+col]
						int col = i_sym_col * PhyConfig->nModType + j;
						// int col = i_sym_col * PhyConfig->nModType + j;
						if (j == 0) { llrBlock[lenMIMOBlock * i_block + i * PhyConfig->nNh + col] = -y_recv_tmp[i].real(); }
						else if (j == 1) { llrBlock[lenMIMOBlock * i_block + i * PhyConfig->nNh + col] = -abs(y_recv_tmp[i].real()) + 0.632455; }
						else if (j == 2) { llrBlock[lenMIMOBlock * i_block + i * PhyConfig->nNh + col] = -y_recv_tmp[i].imag(); }
						else if (j == 3) { llrBlock[lenMIMOBlock * i_block + i * PhyConfig->nNh + col] = -abs(y_recv_tmp[i].imag()) + 0.632455; }
					}
				}
			}
		}
	}
	else {
# pragma omp parallel for 
		for (int i_block = 0; i_block < nMIMOInBlock; i_block++) {
			double noiseVariance_tmp = noiseVariance;
			JacobiSVD<MatrixXcf> svd(channelMatrixArray[i_block], ComputeFullU | ComputeFullV);
			MatrixXcf U = svd.matrixU();
			MatrixXcf V = svd.matrixV();
			pCfloat x_trans_tmp(new std::complex<float>[PhyConfig->nNt]);
			// pCfloat x_trans_ori(new std::complex<float>[PhyConfig->nNt]); // for test
			pCfloat x_precoded_tmp(new std::complex<float>[PhyConfig->nNt]);
			pCfloat y_recv_tmp(new std::complex<float>[PhyConfig->nNr]);
			pCfloat y_preprocessed(new std::complex<float>[PhyConfig->nNr]);
			pFloat singular_values_array(new float[PhyConfig->nNt]);
			auto S = svd.singularValues();
			for (int i = 0; i < PhyConfig->nNt; i++) { x_trans_tmp[i] = 0.0; }
			for (int i = 0; i < PhyConfig->nNr; i++) { y_recv_tmp[i] = 0.0; }
			pFloat antennaPower(new float[PhyConfig->nNt]);
			for (int i = 0; i < PhyConfig->nStream; i++) { antennaPower[i] = (float)PhyConfig->nNt / PhyConfig->nStream; }
			for (int i = 0; i < PhyConfig->nStream; i++) { singular_values_array[i] = S[i]; }
			if (PhyConfig->powerAllocationScheme == "SEPARATE") {
				float my_sum = 0.0;
				// for (int i = 0; i < Nt; i++) { singular_values_array[i] = S[i]; }
				float sum_tmp = 0;
				// int border = 15;
				int border_tmp = nBadChannelsMIMOBlock[i_block];
				// std::cout << std::endl << border_tmp << std::endl;
				for (int i = 0; i < PhyConfig->nStream; i++) {
					if (i > PhyConfig->nStream - border_tmp - 1) {
						sum_tmp += 1.0 / (singular_values_array[i] * singular_values_array[i]) * PhyConfig->powerAllocationRatio;
					}
					else {
						sum_tmp += 1.0 / (singular_values_array[i] * singular_values_array[i]);
					}
				}
				for (int i = 0; i < PhyConfig->nStream; i++) {
					float sing_tmp = 1.0 / (singular_values_array[i] * singular_values_array[i]);
					if (i > PhyConfig->nStream - border_tmp - 1) {
						antennaPower[i] = (float)Nt * sing_tmp / sum_tmp * PhyConfig->powerAllocationRatio;
					}
					else {
						antennaPower[i] = (float)Nt * sing_tmp / sum_tmp;
					}
					my_sum += antennaPower[i];
					/*cout << antenna_power[i] << "  ";
					cout << sqrt(antenna_power[i]) * singular_values_array[i] << "  ";
					cout << endl;*/
				}
				// std::cout<< endl << my_sum << endl;
				for (int i = PhyConfig->nStream; i < PhyConfig->nNt; i++) { antennaPower[i] = 0.0; }
			}
			for (int i = PhyConfig->nStream; i < PhyConfig->nNt; i++) { antennaPower[i] = 0.0; }
			for (int i = 0; i < PhyConfig->nNt; i++) { S[i] = S[i] * sqrt(antennaPower[i]); /*std::cout << S[i] * S[i] << std::endl;*/ }
			for (int i_sym_col = 0; i_sym_col < nSymbols; i_sym_col++) {
				for (int i = 0; i < PhyConfig->nStream; i++) {
					x_trans_tmp[i] = sqrt(antennaPower[i]) * xTransmittedBlock[i_block * lenMIMOBlockSym + i * nSymbols + i_sym_col];
					// x_trans_ori[i] =  xTransmittedBlock[i_block * lenMIMOBlockSym + i * nSymbols + i_sym_col];
				}
				// test: transmitted symbols 
				precodeSVD(PhyConfig->nNr, PhyConfig->nNt, x_trans_tmp.get(), x_precoded_tmp.get(), V);
				AWGN(PhyConfig->nNr, PhyConfig->nNt, x_precoded_tmp.get(), y_recv_tmp.get(), channelMatrixArray[i_block], noiseVariance_tmp);
				precodeSVD_preprocessing(PhyConfig->nNr, PhyConfig->nNt, y_recv_tmp.get(), y_preprocessed.get(), U);
				for (int i = 0; i < PhyConfig->nStream; i++) { y_recv_tmp[i] = S[i] * y_preprocessed[i] / (S[i] * S[i] + noiseVariance); }
				// for (int i = 0; i < PhyConfig->nStream; i++) { std::cout << x_trans_ori[i] << "  " << y_recv_tmp[i] << endl; }
				for (int i = 0; i < PhyConfig->nStream; i++) {
					for (int j = 0; j < PhyConfig->nModType; j++) {
						// int col = i_sym_col*PhyConfig->nModType+j
						// llrBlock[lenMIMOBlock*i_Block + i*PhyConfig->nNh+col]
						int col = i_sym_col * PhyConfig->nModType + j;
						// int col = i_sym_col * PhyConfig->nModType + j;
						if (PhyConfig->nModType == 4) {
							if (j == 0) { llrBlock[lenMIMOBlock * i_block + i * PhyConfig->nNh + col] = -y_recv_tmp[i].real(); }
							else if (j == 1) { llrBlock[lenMIMOBlock * i_block + i * PhyConfig->nNh + col] = -abs(y_recv_tmp[i].real()) + 0.632455; }
							else if (j == 2) { llrBlock[lenMIMOBlock * i_block + i * PhyConfig->nNh + col] = -y_recv_tmp[i].imag(); }
							else if (j == 3) { llrBlock[lenMIMOBlock * i_block + i * PhyConfig->nNh + col] = -abs(y_recv_tmp[i].imag()) + 0.632455; }
						}
						else if (PhyConfig->nModType == 2) {
							if (j == 0) { llrBlock[lenMIMOBlock * i_block + i * PhyConfig->nNh + col] = -y_recv_tmp[i].real(); }
							else if (j == 1) { llrBlock[lenMIMOBlock * i_block + i * PhyConfig->nNh + col] = -y_recv_tmp[i].imag(); }
						}
					}
				}
			}
		}
	}
}


void transmissionMIMOBlockNoPrecode(sPhyConfig* PhyConfig, std::vector<MatrixXcf>& channelMatrixArray, int nMIMOInBlock, int nSymbols, int lenMIMOBlock, pCfloat& xTransmittedBlock, pFloat& llrBlock, pInt& nBadChannelsMIMOBlock, float noiseVariance) {
	// llrBlock is of nMIMOInBlock * nComponentCodewordInMIMO * Nh
	int lenMIMOBlockSym = lenMIMOBlock / PhyConfig->nModType;
	if (false) { ; }
	else {
# pragma omp parallel for 
		for (int i_block = 0; i_block < nMIMOInBlock; i_block++) {
			double noiseVariance_tmp = noiseVariance;
			/*JacobiSVD<MatrixXcf> svd(channelMatrixArray[i_block], ComputeFullU | ComputeFullV);
			MatrixXcf U = svd.matrixU();
			MatrixXcf V = svd.matrixV();*/
			pCfloat x_trans_tmp(new std::complex<float>[PhyConfig->nNt]);
			// pCfloat x_trans_ori(new std::complex<float>[PhyConfig->nNt]); // for test
			pCfloat x_precoded_tmp(new std::complex<float>[PhyConfig->nNt]);
			pCfloat y_recv_tmp(new std::complex<float>[PhyConfig->nNr]);
			pCfloat y_preprocessed(new std::complex<float>[PhyConfig->nNr]);
			pFloat singular_values_array(new float[PhyConfig->nNt]);
			// auto S = svd.singularValues();
			for (int i = 0; i < PhyConfig->nNt; i++) { x_trans_tmp[i] = 0.0; }
			for (int i = 0; i < PhyConfig->nNr; i++) { y_recv_tmp[i] = 0.0; }
			pFloat antennaPower(new float[PhyConfig->nNt]);
			for (int i = 0; i < PhyConfig->nStream; i++) { antennaPower[i] = (float)PhyConfig->nNt / PhyConfig->nStream; }
			for (int i = PhyConfig->nStream; i < PhyConfig->nNt; i++) { antennaPower[i] = 0.0; }
			// for (int i = 0; i < PhyConfig->nNt; i++) { S[i] = S[i] * sqrt(antennaPower[i]); /*std::cout << S[i] * S[i] << std::endl;*/ }
			for (int i_sym_col = 0; i_sym_col < nSymbols; i_sym_col++) {
				for (int i = 0; i < PhyConfig->nStream; i++) {
					x_trans_tmp[i] = sqrt(antennaPower[i]) * xTransmittedBlock[i_block * lenMIMOBlockSym + i * nSymbols + i_sym_col];
					// x_trans_ori[i] =  xTransmittedBlock[i_block * lenMIMOBlockSym + i * nSymbols + i_sym_col];
				}
				// test: transmitted symbols 
				// precodeSVD(PhyConfig->nNr, PhyConfig->nNt, x_trans_tmp.get(), x_precoded_tmp.get(), V);
				AWGN(PhyConfig->nNr, PhyConfig->nNt, x_trans_tmp.get(), y_recv_tmp.get(), channelMatrixArray[i_block], noiseVariance_tmp);
				/*std::cout << x_trans_tmp[0] << std::endl;
				std::cout << noiseVariance << std::endl;*/
				MMSEDetectionCfloat(y_recv_tmp, channelMatrixArray[i_block], noiseVariance, PhyConfig);
				// std::cout << channelMatrixArray[0](0, 0) << std::endl;
				// for (int kkk = 0; kkk < PhyConfig->nNr; kkk++) { std::cout << y_recv_tmp[kkk] << std::endl; }
				// precodeSVD_preprocessing(PhyConfig->nNr, PhyConfig->nNt, y_recv_tmp.get(), y_preprocessed.get(), U);
				// for (int i = 0; i < PhyConfig->nStream; i++) { y_recv_tmp[i] = S[i] * y_preprocessed[i] / (S[i] * S[i] + noiseVariance); }
				// for (int i = 0; i < PhyConfig->nStream; i++) { std::cout << x_trans_ori[i] << "  " << y_recv_tmp[i] << endl; }
				for (int i = 0; i < PhyConfig->nStream; i++) {
					for (int j = 0; j < PhyConfig->nModType; j++) {
						// int col = i_sym_col*PhyConfig->nModType+j
						// llrBlock[lenMIMOBlock*i_Block + i*PhyConfig->nNh+col]
						int col = i_sym_col * PhyConfig->nModType + j;
						// int col = i_sym_col * PhyConfig->nModType + j;
						if (PhyConfig->nModType == 4) {
							if (j == 0) { llrBlock[lenMIMOBlock * i_block + i * PhyConfig->nNh + col] = -y_recv_tmp[i].real(); }
							else if (j == 1) { llrBlock[lenMIMOBlock * i_block + i * PhyConfig->nNh + col] = -abs(y_recv_tmp[i].real()) + 0.632455; }
							else if (j == 2) { llrBlock[lenMIMOBlock * i_block + i * PhyConfig->nNh + col] = -y_recv_tmp[i].imag(); }
							else if (j == 3) { llrBlock[lenMIMOBlock * i_block + i * PhyConfig->nNh + col] = -abs(y_recv_tmp[i].imag()) + 0.632455; }
						}
						else if (PhyConfig->nModType == 2) {
							if (j == 0) { llrBlock[lenMIMOBlock * i_block + i * PhyConfig->nNh + col] = -y_recv_tmp[i].real(); }
							else if (j == 1) { llrBlock[lenMIMOBlock * i_block + i * PhyConfig->nNh + col] = -y_recv_tmp[i].imag(); }
						}
					}
				}
			}
		}
	}
}

void deMapDeIntrlvBlock(sPhyConfig* PhyConfig, pFloat& llrBlock, pFloat& llrDecoderBlock, pInt& idxComponentCodewordInBlock, vector<vector<int>> intrlvPatternTime, vector<vector<int>>& intrlvPatternSpace, vector<vector<int>>& intrlvPatternSpaceSep, int& nComponentCodewordInBlock) {
	for (int i = 0; i < nComponentCodewordInBlock; i++) {
		for (int j = 0; j < PhyConfig->nNh; j++) {
			llrDecoderBlock[i * PhyConfig->nNh + j] = llrBlock[idxComponentCodewordInBlock[i] * PhyConfig->nNh + j];
		}
		randDeintrlv(intrlvPatternTime[i % (PhyConfig->nNv)], llrDecoderBlock.get() + i * PhyConfig->nNh, PhyConfig->nNh); 
	}

	if (PhyConfig->bFlagInterleaveSpatial) {
		for (int i = 0; i < nComponentCodewordInBlock / PhyConfig->nNv; i++) {
			float* intrlv_tmp = new float[PhyConfig->nNv];
			for (int j = 0; j < PhyConfig->nNh; j++) {
				int offset = i * PhyConfig->nNh * PhyConfig->nNv;
				for (int k = 0; k < PhyConfig->nNv; k++) { intrlv_tmp[k] = llrDecoderBlock[offset + k * PhyConfig->nNh + j]; }
				if (PhyConfig->pKTime[0] == PhyConfig->pKTime[1]) { randDeintrlv(intrlvPatternSpace[j], intrlv_tmp, PhyConfig->nNv); }
				else { randDeintrlv(intrlvPatternSpaceSep[j], intrlv_tmp, PhyConfig->nNv); }
				for (int k = 0; k < PhyConfig->nNv; k++) { llrDecoderBlock[offset + k * PhyConfig->nNh + j] = intrlv_tmp[k]; }
			}
			delete[] intrlv_tmp;
		}
	}
}


void deMapDeIntrlvBlock1D(sPhyConfig* PhyConfig, pFloat& llrBlock, pFloat& llrDecoderBlock, pInt& idxComponentCodewordInBlock, vector<int>& intrlvPattern1D, int& nComponentCodewordInBlock, int nWholeCodewordInBlock) {
	for (int i = 0; i < nComponentCodewordInBlock; i++) {
		for (int j = 0; j < PhyConfig->nNh; j++) {
			llrDecoderBlock[i * PhyConfig->nNh + j] = llrBlock[idxComponentCodewordInBlock[i] * PhyConfig->nNh + j];
		}
		// randDeintrlv(intrlvPatternTime[i % (PhyConfig->nNv)], llrDecoderBlock.get() + i * PhyConfig->nNh, PhyConfig->nNh);
	}
	for (int i = 0; i < nWholeCodewordInBlock; i++) {
		randDeintrlv(intrlvPattern1D, llrDecoderBlock.get() + i * PhyConfig->nN, PhyConfig->nN);
	}
	/*if (PhyConfig->bFlagInterleaveSpatial) {
		for (int i = 0; i < nComponentCodewordInBlock / PhyConfig->nNv; i++) {
			float* intrlv_tmp = new float[PhyConfig->nNv];
			for (int j = 0; j < PhyConfig->nNh; j++) {
				int offset = i * PhyConfig->nNh * PhyConfig->nNv;
				for (int k = 0; k < PhyConfig->nNv; k++) { intrlv_tmp[k] = llrDecoderBlock[offset + k * PhyConfig->nNh + j]; }
				if (PhyConfig->pKTime[0] == PhyConfig->pKTime[1]) { randDeintrlv(intrlvPatternSpace[j], intrlv_tmp, PhyConfig->nNv); }
				else { randDeintrlv(intrlvPatternSpaceSep[j], intrlv_tmp, PhyConfig->nNv); }
				for (int k = 0; k < PhyConfig->nNv; k++) { llrDecoderBlock[offset + k * PhyConfig->nNh + j] = intrlv_tmp[k]; }
			}
			delete[] intrlv_tmp;
		}
	}*/
}

void generateIntrlvPatternIrrNew(sPhyConfig* PhyConfig, std::vector<int>& pattern, int len, pInt& isLowRateCodeword) {
	int* lr_codeword_index = new int[PhyConfig->nLowRateCodewords];
	int* hr_codeword_index = new int[len - PhyConfig->nLowRateCodewords];
	int cnt_lr = 0;
	int cnt_hr = 0;
	for (int i = 0; i < len; i++) {
		if (isLowRateCodeword[i] == 1) { lr_codeword_index[cnt_lr] = i; cnt_lr++; }
		else { hr_codeword_index[cnt_hr] = i; cnt_hr++; }
	}
	/*std::cout << cnt_lr;
	std::cout << */
	vector<int> intrlv_pattern_lr;
	vector<int> intrlv_pattern_hr;
	intrlv_pattern_lr.resize(PhyConfig->nLowRateCodewords);
	intrlv_pattern_hr.resize(PhyConfig->nNv - PhyConfig->nLowRateCodewords);
	generateIntrlvPattern(intrlv_pattern_lr, PhyConfig->nLowRateCodewords);
	generateIntrlvPattern(intrlv_pattern_hr, PhyConfig->nNv - PhyConfig->nLowRateCodewords);
	randIntrlv(intrlv_pattern_lr, lr_codeword_index, PhyConfig->nLowRateCodewords);
	randIntrlv(intrlv_pattern_hr, hr_codeword_index, len - PhyConfig->nLowRateCodewords);
	cnt_lr = 0;
	cnt_hr = 0;
	// for (int i = 0; i < len - PhyConfig->nLowRateCodewords; i++) { std::cout << hr_codeword_index[i] << endl; }
	for (int i = 0; i < len; i++) {
		if (isLowRateCodeword[i] == 1) { pattern[i] = lr_codeword_index[cnt_lr]; cnt_lr++; }
		else { pattern[i] = hr_codeword_index[cnt_hr]; cnt_hr++; }
	}
	delete[] lr_codeword_index;
	delete[] hr_codeword_index;
}

float SNR2EbN0(float SNR, sPhyConfig* PhyConfig) {
	/*tmp = pow(10, nowEbN0 / 10) * symbolLength * codeRate * PhyConfig->nNt;
	noiseVariance = (PhyConfig->nNt * PhyConfig->nNr) / tmp;*/
	int nInfo = PhyConfig->pKTime[0] * PhyConfig->nLowRateCodewords + PhyConfig->pKTime[1] * (PhyConfig->nKv - PhyConfig->nLowRateCodewords) - CRCL;
	int codeLength = PhyConfig->nNh * PhyConfig->nNv;
	int iEbN0 = 0;
	int symbolLength = PhyConfig->nModType;
	int nBadChannels = 0;
	float codeRate = (float)nInfo / codeLength;
	float linearSNR = pow(10, SNR / 10.0);
	float noiseVariance = PhyConfig->nNt / linearSNR;
	float linearEbN0 = (PhyConfig->nNt * PhyConfig->nNr) / noiseVariance / (symbolLength * codeRate * PhyConfig->nNt);
	float dBEbN0 = 10.0 * log10(linearEbN0);
	return dBEbN0;
}


float SNR2EbN0_1D(float SNR, sPhyConfig* PhyConfig) {
	/*tmp = pow(10, nowEbN0 / 10) * symbolLength * codeRate * PhyConfig->nNt;
	noiseVariance = (PhyConfig->nNt * PhyConfig->nNr) / tmp;*/
	int nInfo = PhyConfig->nK - CRCL;
	int codeLength = PhyConfig->nN;
	int iEbN0 = 0;
	int symbolLength = PhyConfig->nModType;
	int nBadChannels = 0;
	float codeRate = (float)nInfo / codeLength;
	float linearSNR = pow(10, SNR / 10.0);
	float noiseVariance = PhyConfig->nNt / linearSNR;
	float linearEbN0 = (PhyConfig->nNt * PhyConfig->nNr) / noiseVariance / (symbolLength * codeRate * PhyConfig->nNt);
	float dBEbN0 = 10.0 * log10(linearEbN0);
	return dBEbN0;
}



void generateLargeInfoBlock1D(sPhyConfig* PhyConfig, pInt& idxComponentCodewordInBlock, pInt& idxLowRateCodeword, pInt& idxHighRateCodeword, pInt& isLowRateCodeword, pInt& uInfoBlock, pInt& originalInfoBlock, pInt& xEncodedBlock,
	pInt& A_H, pInt& Ac_H, pInt& A_Ac_H, pInt& A_V, pInt& A_Ac_V, pInt& Ac_V, pInt& kTimeAll, pInt& nBadChannelsMIMOBlock, int nComponentInBlock, int nWholeCodewordInBlock, int nLenBlock, int nLenInfoBlock, int nInfo, int codeLength, int& nBadChannels,
	std::mt19937& gen, std::uniform_int_distribution<> distrib
) {
	pInt idxBadChannelInBlock(new int[nComponentInBlock]); // channels allocated to low rate codewords
	pInt idxGoodChannelInBlock(new int[nComponentInBlock]); //  channels allocated to high rate codewords
// # pragma omp parallel for
	for (int i = 0; i < nWholeCodewordInBlock; i++) {
		unsigned char* infoWithCRC = new unsigned char[PhyConfig->nN];
		for (int j = 0; j < PhyConfig->nN; j++) { infoWithCRC[j] = (unsigned char)0; }
		pInt xHorizontalEncode(new int[PhyConfig->nNh * PhyConfig->nNv]);
		pInt verticalEncodeTmp(new int[PhyConfig->nNv]);
		pInt verticalEncodeTmpMid(new int[PhyConfig->nNv]);
		for (int j = 0; j < PhyConfig->nNh * PhyConfig->nNv; j++) { xHorizontalEncode[j] = 0; }
		for (int j = 0; j < PhyConfig->nNv; j++) { verticalEncodeTmp[j] = 0; verticalEncodeTmpMid[j] = 0; }
		// Generate info bits
		/*for (int j = 0; j < nInfo; j++) { originalInfoBlock[i * codeLength + j] = distrib(gen) % 2; }*/
		for (int j = 0; j < nInfo; j++) { infoWithCRC[j] = (unsigned char) distrib(gen) % 2; }
		if (PhyConfig->nCRCL > 0) { 
			tx_append_crc(infoWithCRC, nInfo, PhyConfig->nCRCL);
		}
		for (int j = 0; j < PhyConfig->nK; j++) { originalInfoBlock[i * codeLength + j] = infoWithCRC[j]; } // prev: nInfo
		// std::cout << rx_check_crc(originalInfoBlock.get() + i * codeLength, PhyConfig->nK, PhyConfig->nCRCL);
		int i_bit_ori = 0;
		//for (int j = 0; j < PhyConfig->nKv; j++) {
		//	for (int i_bit = 0; i_bit < kTimeAll[A_V[j]]; i_bit++) {
		//		uInfoBlock[i * codeLength + PhyConfig->nNh * A_V[j] + A_H[A_V[j] * PhyConfig->nNh + i_bit]] = originalInfoBlock[i * codeLength + i_bit_ori];
		//		// std::cout << i * codeLength + PhyConfig->nNh * A_V[j] + A_H[A_V[j] * PhyConfig->nNh + i_bit] << endl;
		//		i_bit_ori++;
		//	}
		//	// std::cout << endl;
		//}
		for (i_bit_ori = 0; i_bit_ori < PhyConfig->nK; i_bit_ori++) {
			uInfoBlock[i * codeLength + A_H[i_bit_ori]] = originalInfoBlock[i * codeLength + i_bit_ori];// A_H here is actually A.
		}
		// Encoding
		PolarEncode_xor(xHorizontalEncode.get(), uInfoBlock.get() + i * codeLength, PhyConfig->nN);
		// horizontal Encoding
		//for (int j = 0; j < PhyConfig->nKv; j++) {
		//	int i_row = A_V[j];
		//	PolarEncode_xor(xHorizontalEncode.get() + i_row * PhyConfig->nNh, uInfoBlock.get() + i * codeLength + i_row * PhyConfig->nNh, PhyConfig->nNh);
		//}
		//// Vertical Encoding 
		//// int ttttt = 0;
		//for (int j = 0; j < PhyConfig->nNh; j++) {
		//	for (int i_row = 0; i_row < PhyConfig->nNv; i_row++) { verticalEncodeTmp[i_row] = xHorizontalEncode[i_row * PhyConfig->nNh + j]; }
		//	PolarEncode_xor(verticalEncodeTmpMid.get(), verticalEncodeTmp.get(), PhyConfig->nNv);
		//	for (int i_frozen = 0; i_frozen < PhyConfig->nNv - PhyConfig->nKv; i_frozen++) { verticalEncodeTmpMid[Ac_V[i_frozen]] = 0; }
		//	PolarEncode_xor(verticalEncodeTmp.get(), verticalEncodeTmpMid.get(), PhyConfig->nNv);
		//	for (int i_row = 0; i_row < PhyConfig->nNv; i_row++) { xHorizontalEncode[i_row * PhyConfig->nNh + j] = verticalEncodeTmp[i_row]; }
		//}
		// fill to long block
		for (int j = 0; j < codeLength; j++) { xEncodedBlock[i * codeLength + j] = xHorizontalEncode[j];/* std::cout << i * codeLength + j << endl;*/ }
		delete[] infoWithCRC;
	}
    for (int i = 0; i < nComponentInBlock; i++) { idxComponentCodewordInBlock[i] = i; }
}


/**
 * @brief MIMO信道功率分配函数 (基于ZF近似和最小化Worst-case Error Probability)
 *
 * @param PhyConfig 系统配置结构体指针
 * @param channelMatrixArray 信道矩阵数组 (输入原始H, 输出等效H)
 * @param nComponentMaxD 具有maxD距离的子码数量 (映射到最后几根天线)
 * @param minD 较小的最小汉明距离
 * @param maxD 较大的最小汉明距离
 * @param P 总发射功率
 */
void mimoChannelPowerAlloc(sPhyConfig* PhyConfig, std::vector<Eigen::MatrixXcf>& channelMatrixArray,
	int nComponentMaxD, float minD, float maxD, float P) {

	int nNt = PhyConfig->nNt;
	int nNr = PhyConfig->nNr;

	// 安全检查：确保天线配置支持ZF (Nr >= Nt)
	if (nNr < nNt) {
		std::cerr << "Error: ZF approximation requires Nr >= Nt." << std::endl;
		return;
	}

	// 遍历每一个频点/样本的信道矩阵
	for (int k = 0; k < PhyConfig->n_frequency_samples; ++k) {
		Eigen::MatrixXcf& H = channelMatrixArray[k];
		Eigen::MatrixXcf Gram = H.adjoint() * H;
		Eigen::MatrixXcf InvGram = Gram.inverse();
		std::vector<float> weights(nNt);
		float sum_weights = 0.0f;
		int boundaryIdx = nNt - nComponentMaxD;

		for (int i = 0; i < nNt; ++i) {
			float noise_amplification = InvGram(i, i).real();
			float d_i = (i < boundaryIdx) ? (float)minD : (float)maxD;
			weights[i] = noise_amplification / d_i;
			sum_weights += weights[i];
		}
		std::vector<float> allocated_powers(nNt); // 用于存储计算出的功率以便打印验证
		float current_P_sum = 0.0f;
		for (int i = 0; i < nNt; ++i) {
			float p_i = P * (weights[i] / sum_weights);
			allocated_powers[i] = p_i;
			current_P_sum += p_i;
			float scale_factor = std::sqrt(p_i);
			H.col(i) *= scale_factor;
		}
		/*if (k == 0) {
			std::cout << "--- Power Allocation Verification for Channel Sample 0 ---" << std::endl;
			std::cout << "Total Power Constraint P: " << P << std::endl;
			std::cout << "Calculated Power Sum:     " << current_P_sum << std::endl;
			std::cout << "Power Distribution (per antenna):" << std::endl;

			for (int i = 0; i < nNt; ++i) {
				std::string distType = (i < boundaryIdx) ? "minD" : "maxD";
				std::cout << "  Antenna " << i << " (" << distType << "): "
					<< allocated_powers[i] << std::endl;
			}
			std::cout << "--------------------------------------------------------" << std::endl;
		}*/
	}
}

void mimoChannelPowerAllocSeparate(sPhyConfig* PhyConfig, std::vector<Eigen::MatrixXcf>& channelMatrixArray,
	int nComponentMaxD, float minD, float maxD, float P, float noiseVariance) {

	int nNt = PhyConfig->nNt;
	int nNr = PhyConfig->nNr;

	// 安全检查：确保天线配置支持ZF (Nr >= Nt)
	if (nNr < nNt) {
		std::cerr << "Error: ZF approximation requires Nr >= Nt." << std::endl;
		return;
	}
	Eigen::MatrixXcf IdentityMat = Eigen::MatrixXcf::Identity(nNt, nNt);
	for (int i = 0; i < nNt; i++) { for (int j = 0; j < nNt; j++) { IdentityMat(i, j) *= std::complex<float>(0.0*noiseVariance, 0.0); } }
	// 遍历每一个频点/样本的信道矩阵
	for (int k = 0; k < PhyConfig->n_frequency_samples; ++k) {
		Eigen::MatrixXcf& H = channelMatrixArray[k];
		Eigen::MatrixXcf Gram = H.adjoint() * H;
		Eigen::MatrixXcf InvGram = (Gram+IdentityMat).inverse();

		std::vector<float> weights(nNt);
		std::vector<float> gains(nNt);
		float sum_weights = 0.0f;
		int boundaryIdx = nNt - nComponentMaxD;
		float sum_snr_high = 0;
		float sum_snr_low = 0;
		float avg_snr_high = 0;
		float avg_snr_low = 0;
		// std::cout << std::endl;
		for (int i = 0; i < nNt; i++) {
			gains[i] = 1.0 / (abs(InvGram(i, i).real())/**noiseVariance*/);
			// std::cout << gains[i] << std::endl;
		}
		
		
		for (int i = 0; i < nNt; i++) { if (i < boundaryIdx) { sum_snr_high += 1.0/gains[i]; } else { sum_snr_low += 1.0/gains[i]; } }
		avg_snr_high = 1.0 / sum_snr_high * (float)boundaryIdx;
		avg_snr_low = 1.0 / sum_snr_low * (float)nComponentMaxD;
		/*for (int i = 0; i < nNt; i++) { if (i < boundaryIdx) { sum_snr_high += gains[i]; } else { sum_snr_low += gains[i]; } }
		avg_snr_high = sum_snr_high / (float)boundaryIdx;
		avg_snr_low = sum_snr_low / (float)nComponentMaxD;*/
		// std::cout << std::endl;
		for (int i = 0; i < nNt; ++i) {
			// float noise_amplification = InvGram(i, i).real();
			float d_i = (i < boundaryIdx) ? (float)minD : (float)maxD;
			float snr_i = (i < boundaryIdx) ? (float)avg_snr_high : (float)avg_snr_low;
			weights[i] = 1.0 / (snr_i * d_i);
			// std::cout << snr_i << "   "<< d_i<<"   " << weights[i] << std::endl;
			sum_weights += weights[i];
			// std::cout << snr_i << std::endl;
		}
		std::vector<float> allocated_powers(nNt); // 用于存储计算出的功率以便打印验证
		float current_P_sum = 0.0f;
		for (int i = 0; i < nNt; ++i) {
			float p_i = P * (weights[i] / sum_weights);
			allocated_powers[i] = p_i;
			current_P_sum += p_i;
			float scale_factor = std::sqrt(p_i);
			H.col(i) *= scale_factor;
		}
	}
}
void system2D_new(sPhyConfig* PhyConfig, float* BER, float* FER) {
	int nInfo = PhyConfig->pKTime[0] * PhyConfig->nLowRateCodewords + PhyConfig->pKTime[1] * (PhyConfig->nKv - PhyConfig->nLowRateCodewords)-CRCL;
	int codeLength = PhyConfig->nNh * PhyConfig->nNv;
	int iEbN0 = 0;
	int symbolLength = PhyConfig->nModType;
	int nBadChannels = 0;
	float codeRate = (float)nInfo / codeLength;
	float noiseVariance = 0;
	std::random_device rd;
	std::mt19937 gen(rd());            
	std::uniform_int_distribution<> distrib(1, 1000000);
	for (int i = 1; i <= PhyConfig->nStream; i++) {
		if ((float)i / (float)PhyConfig->nStream > (float)PhyConfig->nLowRateCodewords / (float)PhyConfig->nNv) { nBadChannels = i - 1; break; }
	}

	// Array definitions
	pChar originalInfo(new unsigned char[nInfo]);
	pInt originalInfoCRC(new int[PhyConfig->nKv * PhyConfig->nKh]);
	pInt uInfo(new int[PhyConfig->nKv * PhyConfig->nKh]);
	pInt uInfoArray(new int[PhyConfig->nNv * PhyConfig->nNh]);
	pInt xEncodedArray(new int[PhyConfig->nNv * PhyConfig->nNh]);
	pInt A_H(new int[PhyConfig->nNv * PhyConfig->nNh]);
	pInt Ac_H(new int[PhyConfig->nNv * PhyConfig->nNh]);
	pInt A_Ac_H(new int[PhyConfig->nNv * PhyConfig->nNh]);
	pInt A_V(new int[PhyConfig->nNv]);
	pInt Ac_V(new int[PhyConfig->nNv]);
	pInt A_Ac_V(new int[PhyConfig->nNv]);
	pInt kTimeAll(new int[PhyConfig->nNv]);
	pInt isLowRateCodeword(new int[PhyConfig->nNv]);
	pInt idxLowRateCodeword(new int[PhyConfig->nNv]);
	pInt idxHighRateCodeword(new int[PhyConfig->nNv]);
	pInt errorPosArray(new int[PhyConfig->nNv * PhyConfig->nNh]);
	for (int i = 0; i < PhyConfig->nNh * PhyConfig->nNv; i++) { errorPosArray[i] = 0; }
	vector<vector<int>> intrlvPatternTime;
	vector<vector<int>> intrlvPatternSpace;
	vector<vector<int>> intrlvPatternSpaceSep;
	// vector<vector<int>> intrlvPatternSpaceSepH;
	vector<MatrixXcf> channelMatrixArray;
	vector<MatrixXcf> channelMatrixArrayPlatform;
	channelMatrixArrayPlatform.resize(PhyConfig->n_frequency_samples);
	intrlvPatternTime.resize(PhyConfig->nNv);
	for (int i = 0; i < PhyConfig->nNv; i++) { intrlvPatternTime[i].resize(PhyConfig->nNh); }
	for (int i = 0; i < PhyConfig->nNv; i++) {
		generateIntrlvPattern(intrlvPatternTime[i], PhyConfig->nNh);
	}
	intrlvPatternSpace.resize(PhyConfig->nNh);
	for (int i = 0; i < PhyConfig->nNh; i++) { intrlvPatternSpace[i].resize(PhyConfig->nNv); }
	for (int i = 0; i < PhyConfig->nNh; i++) {
		generateIntrlvPattern(intrlvPatternSpace[i], PhyConfig->nNv);
	}
	intrlvPatternSpaceSep.resize(PhyConfig->nNh);
	for (int i = 0; i < PhyConfig->nNh; i++) { intrlvPatternSpaceSep[i].resize(PhyConfig->nNv); }
	/*intrlvPatternSpaceSepL.resize(PhyConfig->nNh);
	for (int j = 0; j < PhyConfig->nNh; j++) { intrlvPatternSpaceSepL[j].resize(PhyConfig->nNv); }*/
	// for (int i = 0; i < PhyConfig->nNv; i++) { intrlvPatternTime[i].resize(PhyConfig->nNv); }
	int** GGh = new int* [PhyConfig->nNh];
	for (int i = 0; i < PhyConfig->nNh; i++) { GGh[i] = new int[PhyConfig->nNh]; }
	int** GGv = new int* [PhyConfig->nNv];
	for (int i = 0; i < Nv; i++) { GGv[i] = new int[PhyConfig->nNv]; }

	// Matrix definitions
	//Eigen::MatrixXf Gh = Fn(PhyConfig->nNh);  // Generate Matrix of horizontal coding
	//Eigen::MatrixXf Gv = Fn(PhyConfig->nNv);  // Generate Matrix of vertical coding
	Eigen::MatrixXcf H(PhyConfig->nNr, PhyConfig->nNt); // Complex channel Matrix H: Nr x Nt
	Eigen::MatrixXf H_real(2 * PhyConfig->nNr, 2 * PhyConfig->nNt); // Real domain channel matrix
	Eigen::MatrixXf receivedSignal(2 * PhyConfig->nNr, 1);
	Eigen::MatrixXf HTH(2 * PhyConfig->nNt, 2 * PhyConfig->nNt);
	Eigen::MatrixXf Identity = MatrixXf::Identity(2 * PhyConfig->nNt, 2 * PhyConfig->nNt);
	Eigen::MatrixXf outputSignal(2 * PhyConfig->nNt, 1);
	Eigen::MatrixXf mmseMatrix(2 * PhyConfig->nNt, 2 * PhyConfig->nNr);
	Eigen::MatrixXf HTHIInv(2 * PhyConfig->nNt, 2 * PhyConfig->nNt);
	std::ifstream platformChannelFile("platform_channel.txt");
	std::ifstream elaaChannelFile("H_eff_power_normalized.txt");
	std::ofstream elaaErrorPattern("error_pattern_elaa.txt", std::ios::app);
	int sorted[8] = { 5,3,1,0,7,2,4,6 };
	if ((PhyConfig->use_platform_channel)||(PhyConfig->use_elaa_channel)) {
		for (int i = 0; i < PhyConfig->n_frequency_samples; i++) {
			MatrixXcf H_tmp(PhyConfig->nNr, PhyConfig->nNt);
			channelMatrixArrayPlatform[i] = H_tmp;
			for (int j = 0; j < PhyConfig->nNr; j++) {
				for (int k = 0; k < PhyConfig->nNt; k++) {
					std::complex<float> h_tmp;
					if (PhyConfig->use_elaa_channel) { elaaChannelFile >> h_tmp; }
					else { platformChannelFile >> h_tmp; }
					channelMatrixArrayPlatform[i](k,j) = h_tmp;
					// std::cout << h_tmp << endl;
				}
			}
		}
	}
	if (PhyConfig->use_elaa_channel) {
		vector<float> gain_vec(PhyConfig->nNt);
		for (int i = 0; i < PhyConfig->nNt; i++) {
			float abs_tmp = 0;
			// int idx = 10;
			for (int idx = 0; idx < 1; idx++) {
				for (int j = 0; j < PhyConfig->nNr; j++) {
					abs_tmp = abs_tmp + channelMatrixArrayPlatform[idx](j, i).imag() * channelMatrixArrayPlatform[idx](j, i).imag() + channelMatrixArrayPlatform[idx](j, i).real() * channelMatrixArrayPlatform[idx](j, i).real();
				}
			}
			gain_vec[i] = abs_tmp;
		}
		std::cout << std::endl;
		for (int i = 0; i < PhyConfig->nNt; i++) { std::cout << gain_vec[i] << std::endl; }
	}
	// TEST
	// TEST
	// loop begins (EbN0)
	for (float nowEbN0 = PhyConfig->minSNR; nowEbN0 <= PhyConfig->maxSNR; nowEbN0=nowEbN0+PhyConfig->stepSNR) {
		vector<MatrixXcf> channelMatrixArrayPlatformCopy;
		channelMatrixArrayPlatformCopy.resize(PhyConfig->n_frequency_samples);
		channelMatrixArrayPlatformCopy = channelMatrixArrayPlatform;
		// calculate noise variance
		int nBitError = 0;
		int nFrameError = 0;
		int nTransFrame = 0;
		float tmp;
		if (!PhyConfig->bMIMO)
		{
			tmp = pow(10, nowEbN0 / 10);
			noiseVariance = (float)1 / (float)tmp / 2 / codeRate;
		}
		else
		{
			tmp = pow(10, nowEbN0 / 10) * symbolLength * codeRate * PhyConfig->nNt;
			noiseVariance = (PhyConfig->nNt * PhyConfig->nNr) / tmp; 
		}
		if (PhyConfig->use_elaa_channel && (PhyConfig->powerAllocationScheme == "SEPARATE")) {
			// mimoChannelPowerAllocSeparate(PhyConfig, channelMatrixArrayPlatformCopy, 4, 1, 3.0, PhyConfig->nNt, 0.1*noiseVariance);// prev: 0.001(0.1),3.0
			// mimoChannelPowerAllocSeparate(PhyConfig, channelMatrixArrayPlatformCopy, 4, 1, 3.0, PhyConfig->nNt, 0.1 * noiseVariance);// prev: 0.001(0.1),3.0
		}
		// polar code construction
		Eigen::MatrixXf Gv = Fn(PhyConfig->nNv);
		Eigen::MatrixXf Gh = Fn(PhyConfig->nNh);
		for (int i = 0; i < PhyConfig->nNv; i++) {
			for (int j = 0; j < PhyConfig->nNv; j++) {
				GGv[i][j] = Gv(i,j);
			}
		}

		for (int i = 0; i < PhyConfig->nNh; i++) {
			for (int j = 0; j < PhyConfig->nNh; j++) {
				GGh[i][j] = Gh(i, j);
			}
		}
		polar_code_construction_new(PhyConfig, noiseVariance, A_H, Ac_H, A_Ac_H, A_V, Ac_V, A_Ac_V, kTimeAll, isLowRateCodeword,
			GGh, GGv);
		// TEST
		if (PhyConfig->pKTime[0] != PhyConfig->pKTime[1] && PhyConfig->bFlagInterleaveSpatial) {
			for (int i = 0; i < PhyConfig->nNh; i++) {
				generateIntrlvPatternIrrNew(PhyConfig, intrlvPatternSpaceSep[i], PhyConfig->nNv, isLowRateCodeword);
				// for (int j = 0; j < PhyConfig->nNh; j++) { std::cout << intrlvPatternSpaceSep[i][j] << "  "; }
			}
		}
		for (int i = 0; i < PhyConfig->nNv; i++) {
			int nLowTmp = 0;
			int nHighTmp = 0;
			if (isLowRateCodeword[i]) { idxLowRateCodeword[nLowTmp] = i; nLowTmp++; }
			else { idxHighRateCodeword[nLowTmp] = i; nHighTmp++; }
		}
		// transmission begins
		int nComponentCodewordInBlock = lcm(PhyConfig->nStream, PhyConfig->nNv); // n of horizontal component codewords in a single block
		int nWholeCodewordInBlock = round(nComponentCodewordInBlock / PhyConfig->nNv); // n of entire codewords in a single block
		int lenMIMOBlock = PhyConfig->nStream * PhyConfig->nNh;
		int nLenBlock = nWholeCodewordInBlock * codeLength;  // n of bits in the entire block
		int nLenInfoBlock = nWholeCodewordInBlock * nInfo;
		int nMIMOInBlock = round(nComponentCodewordInBlock / PhyConfig->nStream);
		int nSymbols = PhyConfig->nNh / PhyConfig->nModType;
		pInt idxComponentCodewordInBlock(new int[nWholeCodewordInBlock * PhyConfig->nNv]); // index of horizontal codeword in whole block
		pInt uInfoBlock(new int[nLenBlock]);
		pInt uInfoBlockEsti(new int[nLenBlock]);
		pInt originalInfoBlock(new int[nLenBlock]);
		pInt xEncodedBlock(new int[nLenBlock]);
		pInt xEncodedBlockNoIntrlv(new int[nLenBlock]);
		pInt nBadChannelsMIMOBlock(new int[nMIMOInBlock]); // n of bad channels identified in each MIMO block;
		// std::cout << nLenBlock << endl;
		pInt xTransmittedBlock(new int[nLenBlock]);
		pCfloat xModulationBlock(new std::complex<float>[nLenBlock / PhyConfig->nModType]);
		pFloat llrBlock(new float[nLenBlock]);
		pFloat llrDecoderBlock(new float[nLenBlock]);
		pFloat upLLRBlock(new float[nLenBlock]);
		channelMatrixArray.resize(nMIMOInBlock);
		for (int i = 0; i < nMIMOInBlock; i++) {
			Eigen::MatrixXcf H_tmp(PhyConfig->nNr, PhyConfig->nNt);
			channelMatrixArray[i] = H_tmp;
		}
		for (int i = 0; i < nComponentCodewordInBlock; i++) {
			for (int j = 0; j < PhyConfig->nNh; j++) {
				uInfoBlock[i * PhyConfig->nNh + j] = 0;
				originalInfoBlock[i * PhyConfig->nNh + j] = 0;
				xEncodedBlock[i * PhyConfig->nNh + j] = 0;
			}
		}
		for (int i = 0; i < nLenBlock; i++) { llrBlock[i] = 0.0; }
		int nowFrame = 0;
		while (nFrameError < PhyConfig->nMaxFrameError)
		{
			if (PhyConfig->bEarlyStop) {
				// std::cout << nowFrame << " " << PhyConfig->frameThres << std::endl;
				if (nowFrame > PhyConfig->frameThres) {
					float nowFERTmp = (float)nFrameError / (float)nowFrame;
					// std::cout << "fer" << nowFERTmp << PhyConfig->ferThres std::endl;
					if (nowFERTmp < PhyConfig->ferThres) {
						FER[0] = nowFERTmp;
						break;
					}
				}
			}
			generateLargeInfoBlock(
				PhyConfig,
				idxComponentCodewordInBlock,
				idxLowRateCodeword,
				idxHighRateCodeword,
				isLowRateCodeword,
				uInfoBlock,
				originalInfoBlock,
				xEncodedBlock,
				A_H,
				Ac_H,
				A_Ac_H,
				A_V,
				A_Ac_V,
				Ac_V,
				kTimeAll,
				nBadChannelsMIMOBlock,
				nComponentCodewordInBlock,
				nWholeCodewordInBlock,
				nLenBlock,
				nLenInfoBlock,
				nInfo,
				codeLength,
				nBadChannels,
				gen,
				distrib
			);
			// for (int i = 0; i < nWholeCodewordInBlock * PhyConfig->nNv; i++) { std::cout << i << " " << idxComponentCodewordInBlock[i] << std::endl; }
			for (int i = 0; i < nLenBlock; i++) { xEncodedBlockNoIntrlv[i] = xEncodedBlock[i]; }
			randIntrlvWholeBlock(PhyConfig, xEncodedBlock, intrlvPatternTime,intrlvPatternSpace,intrlvPatternSpaceSep, nWholeCodewordInBlock, codeLength);
			codewordToAntennaMapping(PhyConfig, idxComponentCodewordInBlock, xTransmittedBlock, xEncodedBlock, nComponentCodewordInBlock);
			modulationBlock(PhyConfig, xTransmittedBlock, xModulationBlock, nComponentCodewordInBlock);
			if ((PhyConfig->use_platform_channel) || (PhyConfig->use_elaa_channel) ) {
				generateChannelsBlock(PhyConfig, channelMatrixArray, nMIMOInBlock, channelMatrixArrayPlatformCopy, gen, distrib);
			}
			else {
				generateChannelsBlock(PhyConfig, channelMatrixArray, nMIMOInBlock);
			}
			if (PhyConfig->use_elaa_channel) {
				transmissionMIMOBlockNoPrecode(PhyConfig, channelMatrixArray, nMIMOInBlock, nSymbols, lenMIMOBlock, xModulationBlock, llrBlock, nBadChannelsMIMOBlock, noiseVariance);
			}
			else {
				transmissionMIMOBlock(PhyConfig, channelMatrixArray, nMIMOInBlock, nSymbols, lenMIMOBlock, xModulationBlock, llrBlock, nBadChannelsMIMOBlock, noiseVariance);
			}
			deMapDeIntrlvBlock(PhyConfig, llrBlock, llrDecoderBlock, idxComponentCodewordInBlock, intrlvPatternTime, intrlvPatternSpace, intrlvPatternSpaceSep, nComponentCodewordInBlock);
			for (int i = 0; i < nLenBlock; i++) { upLLRBlock[i] = llrDecoderBlock[i]; }
			// Decoding
			int Nmax = max(PhyConfig->nNh, PhyConfig->nNv);
			for (int i_sub_frame = 0; i_sub_frame < nWholeCodewordInBlock; i_sub_frame++) {
				int frame_llr_offset = i_sub_frame * codeLength;
				float** oriLLRh = new float* [PhyConfig->nNv];
				for (int i = 0; i < Nv; i++) { oriLLRh[i] = new float[PhyConfig->nNh]; }
				float** upLLR = new float* [PhyConfig->nNv];
				for (int i = 0; i < Nv; i++) { upLLR[i] = new float[PhyConfig->nNh]; }
				double* sum_sdis_array = new double[PhyConfig->nNh];
				for (int i = 0; i < PhyConfig->nNv; i++) {
					for (int j = 0; j < PhyConfig->nNh; j++) {
						upLLR[i][j] = upLLRBlock[i_sub_frame * codeLength + i * PhyConfig->nNh + j];
						oriLLRh[i][j] = upLLR[i][j];
					}
				}
				for (int iter = 1; iter <= softiter; iter++)
				{
#pragma omp parallel for
					for (int decodingi = 0; decodingi < PhyConfig->nNh; decodingi++) // 291.70kb
					{

						//for SCL Decoding
						int Counth = 0, Countinfoh = 0;
						float Qsumunh = 0;
						int** sum = new int* [PhyConfig->L];
						//for (int i = 0; i < L; i++)  sum[i] = new  int[Nmax];
						int** u1 = new int* [PhyConfig->L];
						for (int i = 0; i < PhyConfig->L; i++)  u1[i] = new _MM_ALIGN16 int[Nmax];
						for (int i = 0; i < PhyConfig->L; i++)  sum[i] = new _MM_ALIGN16 int[Nmax];
						//for (int i = 0; i < L; i++)  u1[i] = new  int[Nmax];

						float* LLR_in = new float[PhyConfig->L];
						float* W = new float[2 * PhyConfig->L];
						int* SCL_index = new int[PhyConfig->L];
						int* better = new int[PhyConfig->L];
						int* worse = new int[PhyConfig->L];
						int* indexpm = new int[PhyConfig->L];
						float* PM = new float[PhyConfig->L];
						float** LLR = new float* [PhyConfig->L];
						for (int i = 0; i < PhyConfig->L; i++)  LLR[i] = new _MM_ALIGN16 float[2 * Nmax + 2];
						// for (int i = 0; i < L; i++)  LLR[i] = new float[2 * Nmax + 2];
						////////////////////////////////
						for (int j = 0; j < PhyConfig->nNv; j++)
							for (int k = 0; k < PhyConfig->L; k++)
								LLR[k][j] = upLLR[j][decodingi];
								// LLR[k][j] = upLLRBlock[frame_llr_offset + j * PhyConfig->nNh + decodingi];
						//cout << LLR[0][j] << " ";

					//cout << endl;
						int index = 0;
						for (int k = 0; k < PhyConfig->L; k++) { PM[k] = 0; indexpm[k] = k; }
						decode_list(LLR, sum, (int)(log(Nv) / log(2)), A_Ac_V.get(), 0, PhyConfig->nNv, PhyConfig->nKv, PM, L, (int)(log(PhyConfig->L) / log(2)), LLR_in, W, SCL_index, better, worse, &Counth, &Countinfoh, &Qsumunh);
						for (int i = 0; i < PhyConfig->L - 1; i++)
							for (int j = i + 1; j < PhyConfig->L; j++)
								if (PM[indexpm[i]] < PM[indexpm[j]])
								{
									int ttt = indexpm[i];
									indexpm[i] = indexpm[j];
									indexpm[j] = ttt;
								}
						index = indexpm[0];
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
						/*if (half_iter && iter == softiter)
						{
							if (iter == softiter) PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), Nh);
							for (int j = 0; j < Kv; j++)
								umid[A_V[j]][decodingi] = u1[index][A_V[j]];
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
#pragma omp parallel for
					for (int decodingi = 0; decodingi < PhyConfig->nNv; decodingi++)
					{
						// 293.70kb
						// cout << "KKKK" << endl;
						if (half_iter && (iter == softiter)) { break; }
						string polarde_now = polarde_array[iter - 1];
						//for SCL Decoding
						int Countv = 0, Countinfov = 0;
						float Qsumunv = 0;
						int** u1 = new int* [PhyConfig->L];
						for (int i = 0; i < PhyConfig->L; i++)  u1[i] = new _MM_ALIGN16 int[Nmax];
						int** u1_LSD = new int* [2 * L_LSD];
						for (int i = 0; i < 2 * L_LSD; i++) u1_LSD[i] = new int[Nmax];
						int** sum = new int* [PhyConfig->L];
						for (int i = 0; i < PhyConfig->L; i++)  sum[i] = new _MM_ALIGN16 int[Nmax];
						int** sum_LSD = new int* [2 * L_LSD];
						for (int i = 0; i < 2 * L_LSD; i++)  sum_LSD[i] = new int[Nmax];
						float* LLR_in = new float[PhyConfig->L];
						int* u_decod = new int[PhyConfig->nNh];
						int* n_paths = new int[PhyConfig->nNh];
						float* W = new float[2 * PhyConfig->L];
						int* SCL_index = new int[PhyConfig->L];
						int* better = new int[PhyConfig->L];
						int* worse = new int[PhyConfig->L];
						int* indexpm = new int[L_LSD];
						float* PM = new float[PhyConfig->L];
						float** LLR = new float* [PhyConfig->L];
						for (int i = 0; i < PhyConfig->L; i++)  LLR[i] = new _MM_ALIGN16 float[2 * Nmax + 2];

						for (int i = 0; i < L_LSD; i++) { indexpm[i] = i; }
						////////////////////////////////
						for (int j = 0; j < PhyConfig->nNh; j++)
							for (int k = 0; k < PhyConfig->L; k++)
								LLR[k][j] = upLLR[decodingi][j];
						//	cout << LLR[0][j] << " ";
						// calc_npaths(A_Ac_H.get(), Nh, n_paths, L_LSD);
						//	cout << endl;
						int index = 0;
						if (polarde_now == "SCL")
						{
							for (int k = 0; k < PhyConfig->L; k++)  PM[k] = 0, indexpm[k] = k;
							//T_start = clock();
							// decode_list(LLR, sum, (int)(log(Nh) / log(2)), A_Ach, 0, Nh, Kh, PM, L, (int)(log(L) / log(2)), LLR_in, W, SCL_index, better, worse, &Countv, &Countinfov, &Qsumunv);
							decode_list(LLR, sum, (int)(log(Nh) / log(2)), A_Ac_H.get() + decodingi * PhyConfig->nNh, 0, PhyConfig->nNh, kTimeAll[decodingi], PM, PhyConfig->L, (int)(log(PhyConfig->L) / log(2)), LLR_in, W, SCL_index, better, worse, &Countv, &Countinfov, &Qsumunv);
							//T_finish = clock();
							//Decoding_Duration += (T_finish - T_start);
							for (int i = 0; i < PhyConfig->L - 1; i++)
								for (int j = i + 1; j < PhyConfig->L; j++)
									if (PM[indexpm[i]] < PM[indexpm[j]])
									{
										int ttt = indexpm[i];
										indexpm[i] = indexpm[j];
										indexpm[j] = ttt;
									}
							if (iter == softiter) {
								PolarEncode_xor(*(u1 + indexpm[0]), *(sum + indexpm[0]), Nh);
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
									if (softmethod == "Pyndiah") Pyndiahh_irregular(LLR, sum, oriLLRh, indexpm, iter, decodingi, upLLR, sum_sdis_array[decodingi], isLowRateCodeword.get());
									if (softmethod == "MITSO") MITSOh_irregular(LLR, PM, sum, oriLLRh, indexpm, iter, decodingi, upLLR, Qsumunv, isLowRateCodeword.get());
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
								for (int j = 0; j < PhyConfig->nNh; j++)
									// umid[decodingi][j] = u1[index][j];
									uInfoBlockEsti[frame_llr_offset + decodingi * PhyConfig->nNh + j] = u1[index][j];
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
				for (int i = 0; i < PhyConfig->nNv; i++) { delete[] oriLLRh[i]; delete[] upLLR[i]; }
				delete[] oriLLRh; delete[] upLLR;
				delete[] sum_sdis_array;
			}
			for (int i_sub_frame = 0; i_sub_frame < nWholeCodewordInBlock; i_sub_frame++) {
				int singleFrameError = 0;
				for (int i_comp = 0; i_comp < PhyConfig->nKv; i_comp++) {
					for (int i_bit = 0; i_bit < kTimeAll[A_V[i_comp]]; i_bit++) {
						if (uInfoBlockEsti[i_sub_frame * codeLength + A_V[i_comp] * PhyConfig->nNh + A_H[A_V[i_comp] * PhyConfig->nNh + i_bit]] !=
							uInfoBlock[i_sub_frame * codeLength + A_V[i_comp] * PhyConfig->nNh + A_H[A_V[i_comp] * PhyConfig->nNh + i_bit]]) {
							singleFrameError++;
							errorPosArray[A_V[i_comp] * PhyConfig->nNh + A_H[A_V[i_comp] * PhyConfig->nNh + i_bit]]++;
						}
					}
				}
				if (singleFrameError > 0) { nBitError += singleFrameError; nFrameError++; }
				nowFrame++;
			}
			if ((nowFrame % nWholeCodewordInBlock) == 0) {
				printf("\r<%d> SNR= %.2f FER= %d / %d = %f BER= %d / %d = %lf",
					nowFrame, nowEbN0, nFrameError, nowFrame, (double)nFrameError / nowFrame, nBitError, (nInfo * nowFrame), (double)nBitError / nInfo / nowFrame);
				fflush(stdout);
			}
			if ((nowFrame % 1000) == 0) {
				elaaErrorPattern << std::endl;
				for (int i = 0; i < PhyConfig->nNv; i++) {
					for (int j = 0; j < PhyConfig->nNh; j++) {
						elaaErrorPattern << errorPosArray[i * PhyConfig->nNh + j] << "  ";
					}
					elaaErrorPattern << std::endl;
				}
			}
		}
		//Decoding_Duration = Decoding_Duration * 100 / CLOCKS_PER_SEC / frame;
		//Decoding_Duration = Decoding_Duration / CLOCKS_PER_SEC / Max_frame;
		std::cout << "  " << nowEbN0 << "  " << (double)nFrameError / nowFrame << "  " << (double)nBitError / nInfo / nowFrame << std::endl;
		FER[iEbN0] = (double)nFrameError / nowFrame;
	}
		
	iEbN0++;

	for (int i = 0; i < Nh; i++) { delete[] GGh[i]; } delete[] GGh;
	for (int i = 0; i < Nv; i++) { delete[] GGv[i]; } delete[] GGv;
}

void setOptParam(sPhyConfig* PhyConfig, soptParams* OptParams, float errorRate=0.0) {
	OptParams->border = PhyConfig->nLowRateCodewords;
	OptParams->nNh = PhyConfig->nNh;
	OptParams->nNv = PhyConfig->nNv;
	OptParams->nKv = PhyConfig->nKv;
	for (int i = 0; i < PhyConfig->nTimeConstructions; i++) { OptParams->pKTime[i] = PhyConfig->pKTime[i]; }
	OptParams->codeRate = (float)(OptParams->border * OptParams->pKTime[0] + (OptParams->nKv - OptParams->border) * OptParams->pKTime[1]) / (float)(OptParams->nNh * OptParams->nNv);
	OptParams->nStream = PhyConfig->nStream;
	OptParams->transmissionRate = OptParams->nStream * OptParams->codeRate * PhyConfig->nModType;
	OptParams->errorRate = errorRate;
}



void system1D_new(sPhyConfig* PhyConfig, float* BER, float* FER) {
	// int nInfo = PhyConfig->pKTime[0] * PhyConfig->nLowRateCodewords + PhyConfig->pKTime[1] * (PhyConfig->nKv - PhyConfig->nLowRateCodewords) - CRCL;
	int nInfo = PhyConfig->nK - PhyConfig->nCRCL;
	int codeLength = PhyConfig->nN;
	int iEbN0 = 0;
	int symbolLength = PhyConfig->nModType;
	int nBadChannels = 0;
	float codeRate = (float)nInfo / codeLength;
	float noiseVariance = 0;
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> distrib(1, 1000000);
	for (int i = 1; i <= PhyConfig->nStream; i++) {
		if ((float)i / (float)PhyConfig->nStream > (float)PhyConfig->nLowRateCodewords / (float)PhyConfig->nNv) { nBadChannels = i - 1; break; }
	}

	// Array definitions
	pChar originalInfo(new unsigned char[nInfo]);
	pInt originalInfoCRC(new int[PhyConfig->nKv * PhyConfig->nKh]);
	pInt uInfo(new int[PhyConfig->nKv * PhyConfig->nKh]);
	pInt uInfoArray(new int[PhyConfig->nNv * PhyConfig->nNh]);
	pInt xEncodedArray(new int[PhyConfig->nNv * PhyConfig->nNh]);
	pInt A_H(new int[PhyConfig->nNv * PhyConfig->nNh]);
	pInt Ac_H(new int[PhyConfig->nNv * PhyConfig->nNh]);
	pInt A_Ac_H(new int[PhyConfig->nNv * PhyConfig->nNh]);
	pInt A_V(new int[PhyConfig->nNv]);
	pInt Ac_V(new int[PhyConfig->nNv]);
	pInt A_Ac_V(new int[PhyConfig->nNv]);
	pInt A(new int[PhyConfig->nN]);
	pInt Ac(new int[PhyConfig->nN - PhyConfig->nK]);
	pInt A_Ac(new int[PhyConfig->nN]);
	pInt kTimeAll(new int[PhyConfig->nNv]);
	pInt isLowRateCodeword(new int[PhyConfig->nNv]);
	pInt idxLowRateCodeword(new int[PhyConfig->nNv]);
	pInt idxHighRateCodeword(new int[PhyConfig->nNv]);
	vector<int> intrlvPattern1D;
	intrlvPattern1D.resize(PhyConfig->nN);
	generateIntrlvPattern(intrlvPattern1D, PhyConfig->nN);
	/*vector<vector<int>> intrlvPatternTime;
	vector<vector<int>> intrlvPatternSpace;
	vector<vector<int>> intrlvPatternSpaceSep;*/
	// vector<vector<int>> intrlvPatternSpaceSepH;
	vector<MatrixXcf> channelMatrixArray;
	/*intrlvPatternTime.resize(PhyConfig->nNv);
	for (int i = 0; i < PhyConfig->nNv; i++) { intrlvPatternTime[i].resize(PhyConfig->nNh); }
	for (int i = 0; i < PhyConfig->nNv; i++) {
		generateIntrlvPattern(intrlvPatternTime[i], PhyConfig->nNh);
	}
	intrlvPatternSpace.resize(PhyConfig->nNh);
	for (int i = 0; i < PhyConfig->nNh; i++) { intrlvPatternSpace[i].resize(PhyConfig->nNv); }
	for (int i = 0; i < PhyConfig->nNh; i++) {
		generateIntrlvPattern(intrlvPatternSpace[i], PhyConfig->nNv);
	}
	intrlvPatternSpaceSep.resize(PhyConfig->nNh);
	for (int i = 0; i < PhyConfig->nNh; i++) { intrlvPatternSpaceSep[i].resize(PhyConfig->nNv); }*/
	/*intrlvPatternSpaceSepL.resize(PhyConfig->nNh);
	for (int j = 0; j < PhyConfig->nNh; j++) { intrlvPatternSpaceSepL[j].resize(PhyConfig->nNv); }*/
	// for (int i = 0; i < PhyConfig->nNv; i++) { intrlvPatternTime[i].resize(PhyConfig->nNv); }
	/*int** GGh = new int* [PhyConfig->nNh];
	for (int i = 0; i < PhyConfig->nNh; i++) { GGh[i] = new int[PhyConfig->nNh]; }
	int** GGv = new int* [PhyConfig->nNv];
	for (int i = 0; i < PhyConfig->nNv; i++) { GGv[i] = new int[PhyConfig->nNv]; }*/
	int** GG = new int* [PhyConfig->nN];

	for (int i = 0; i < PhyConfig->nN; i++) { GG[i] = new int[PhyConfig->nN]; }

	// Matrix definitions
	//Eigen::MatrixXf Gh = Fn(PhyConfig->nNh);  // Generate Matrix of horizontal coding
	//Eigen::MatrixXf Gv = Fn(PhyConfig->nNv);  // Generate Matrix of vertical coding
	Eigen::MatrixXcf H(PhyConfig->nNr, PhyConfig->nNt); // Complex channel Matrix H: Nr x Nt
	Eigen::MatrixXf H_real(2 * PhyConfig->nNr, 2 * PhyConfig->nNt); // Real domain channel matrix
	Eigen::MatrixXf receivedSignal(2 * PhyConfig->nNr, 1);
	Eigen::MatrixXf HTH(2 * PhyConfig->nNt, 2 * PhyConfig->nNt);
	Eigen::MatrixXf Identity = MatrixXf::Identity(2 * PhyConfig->nNt, 2 * PhyConfig->nNt);
	Eigen::MatrixXf outputSignal(2 * PhyConfig->nNt, 1);
	Eigen::MatrixXf mmseMatrix(2 * PhyConfig->nNt, 2 * PhyConfig->nNr);
	Eigen::MatrixXf HTHIInv(2 * PhyConfig->nNt, 2 * PhyConfig->nNt);
	// vector<MatrixXcf> channelMatrixArray;
	vector<MatrixXcf> channelMatrixArrayPlatform;
	channelMatrixArrayPlatform.resize(PhyConfig->n_frequency_samples);
	std::ifstream platformChannelFile("platform_channel.txt");
	std::ifstream elaaChannelFile("H_eff_power_normalized.txt");
	std::ofstream elaaErrorPattern("error_pattern_elaa.txt", std::ios::app);
	int sorted[8] = { 5,3,1,0,7,2,4,6 };
	if ((PhyConfig->use_platform_channel) || (PhyConfig->use_elaa_channel)) {
		for (int i = 0; i < PhyConfig->n_frequency_samples; i++) {
			MatrixXcf H_tmp(PhyConfig->nNr, PhyConfig->nNt);
			channelMatrixArrayPlatform[i] = H_tmp;
			for (int j = 0; j < PhyConfig->nNr; j++) {
				for (int k = 0; k < PhyConfig->nNt; k++) {
					std::complex<float> h_tmp;
					if (PhyConfig->use_elaa_channel) { elaaChannelFile >> h_tmp; }
					else { platformChannelFile >> h_tmp; }
					channelMatrixArrayPlatform[i](k, j) = h_tmp;
					// std::cout << h_tmp << endl;
				}
			}
		}
	}
	if (PhyConfig->use_elaa_channel) {
		vector<float> gain_vec(PhyConfig->nNt);
		for (int i = 0; i < PhyConfig->nNt; i++) {
			float abs_tmp = 0;
			// int idx = 10;
			for (int idx = 0; idx < 1; idx++) {
				for (int j = 0; j < PhyConfig->nNr; j++) {
					abs_tmp = abs_tmp + channelMatrixArrayPlatform[idx](j, i).imag() * channelMatrixArrayPlatform[idx](j, i).imag() + channelMatrixArrayPlatform[idx](j, i).real() * channelMatrixArrayPlatform[idx](j, i).real();
				}
			}
			gain_vec[i] = abs_tmp;
		}
		std::cout << std::endl;
		for (int i = 0; i < PhyConfig->nNt; i++) { std::cout << gain_vec[i] << std::endl; }
	}


	// loop begins (EbN0)
	for (float nowEbN0 = PhyConfig->minSNR; nowEbN0 <= PhyConfig->maxSNR; nowEbN0 = nowEbN0 + PhyConfig->stepSNR) {
		// calculate noise variance
		int nBitError = 0;
		int nFrameError = 0;
		int nTransFrame = 0;
		float tmp;
		if (!PhyConfig->bMIMO)
		{
			tmp = pow(10, nowEbN0 / 10);
			noiseVariance = (float)1 / (float)tmp / 2 / codeRate;
		}
		else
		{
			tmp = pow(10, nowEbN0 / 10) * symbolLength * codeRate * PhyConfig->nNt;
			noiseVariance = (PhyConfig->nNt * PhyConfig->nNr) / tmp;
		}
		// polar code construction
		Eigen::MatrixXf Gv = Fn(PhyConfig->nNv);
		Eigen::MatrixXf Gh = Fn(PhyConfig->nNh);
		// Eigen::MatrixXf G = Fn(PhyConfig->nN);
		std::ifstream  genMatFile("genMat1024.txt");
		for (int i = 0; i < PhyConfig->nN; i++) {
			for (int j = 0; j < PhyConfig->nN; j++) {
				genMatFile >> GG[i][j];
				if (GG[i][j] != 0 && GG[i][j] != 1) { std::cout << "ERROR" << std::endl; }
			}
		}
		
		genMatFile.close();
		/*for (int i = 0; i < PhyConfig->nNv; i++) {
			for (int j = 0; j < PhyConfig->nNv; j++) {
				GGv[i][j] = Gv(i, j);
			}
		}

		for (int i = 0; i < PhyConfig->nNh; i++) {
			for (int j = 0; j < PhyConfig->nNh; j++) {
				GGh[i][j] = Gh(i, j);
			}
		}*/
		/*polar_code_construction_new(PhyConfig, noiseVariance, A_H, Ac_H, A_Ac_H, A_V, Ac_V, A_Ac_V, kTimeAll, isLowRateCodeword,
			GGh, GGv);*/
		//if (PhyConfig->pKTime[0] != PhyConfig->pKTime[1] && PhyConfig->bFlagInterleaveSpatial) {
		//	for (int i = 0; i < PhyConfig->nNh; i++) {
		//		generateIntrlvPatternIrrNew(PhyConfig, intrlvPatternSpaceSep[i], PhyConfig->nNv, isLowRateCodeword);
		//		// for (int j = 0; j < PhyConfig->nNh; j++) { std::cout << intrlvPatternSpaceSep[i][j] << "  "; }
		//	}
		//}
		polar_code_construction_new_1D(PhyConfig, noiseVariance, A, A_Ac, Ac, GG);
		std::cout << std::endl;
		for (int i = 0; i < PhyConfig->nN; i++) { std::cout << A_Ac[i] << ","; }
		std::cout << std::endl;
		/*for (int i = 0; i < PhyConfig->nNv; i++) {
			int nLowTmp = 0;
			int nHighTmp = 0;
			if (isLowRateCodeword[i]) { idxLowRateCodeword[nLowTmp] = i; nLowTmp++; }
			else { idxHighRateCodeword[nLowTmp] = i; nHighTmp++; }
		}*/
		// transmission begins
		int nComponentCodewordInBlock = lcm(PhyConfig->nStream, PhyConfig->nNv); // n of horizontal component codewords in a single block
		int nWholeCodewordInBlock = round(nComponentCodewordInBlock / PhyConfig->nNv); // n of entire codewords in a single block
		int lenMIMOBlock = PhyConfig->nStream * PhyConfig->nNh;
		int nLenBlock = nWholeCodewordInBlock * codeLength;  // n of bits in the entire block
		int nLenInfoBlock = nWholeCodewordInBlock * nInfo;
		int nMIMOInBlock = round(nComponentCodewordInBlock / PhyConfig->nStream);
		int nSymbols = PhyConfig->nNh / PhyConfig->nModType;
		pInt idxComponentCodewordInBlock(new int[nWholeCodewordInBlock * PhyConfig->nNv]); // index of horizontal codeword in whole block
		pInt uInfoBlock(new int[nLenBlock]);
		pInt uInfoBlockEsti(new int[nLenBlock]);
		pInt originalInfoBlock(new int[nLenBlock]);
		pInt xEncodedBlock(new int[nLenBlock]);
		pInt xEncodedBlockNoIntrlv(new int[nLenBlock]);
		pInt nBadChannelsMIMOBlock(new int[nMIMOInBlock]); // n of bad channels identified in each MIMO block;
		// std::cout << nLenBlock << endl;
		pInt xTransmittedBlock(new int[nLenBlock]);
		pCfloat xModulationBlock(new std::complex<float>[nLenBlock / PhyConfig->nModType]);
		pFloat llrBlock(new float[nLenBlock]);
		pFloat llrDecoderBlock(new float[nLenBlock]);
		pFloat upLLRBlock(new float[nLenBlock]);
		channelMatrixArray.resize(nMIMOInBlock);
		for (int i = 0; i < nMIMOInBlock; i++) {
			Eigen::MatrixXcf H_tmp(PhyConfig->nNr, PhyConfig->nNt);
			channelMatrixArray[i] = H_tmp;
		}
		for (int i = 0; i < nComponentCodewordInBlock; i++) {
			for (int j = 0; j < PhyConfig->nNh; j++) {
				uInfoBlock[i * PhyConfig->nNh + j] = 0;
				originalInfoBlock[i * PhyConfig->nNh + j] = 0;
				xEncodedBlock[i * PhyConfig->nNh + j] = 0;
			}
		}
		for (int i = 0; i < nLenBlock; i++) { llrBlock[i] = 0.0; }
		int nowFrame = 0;
		while (nFrameError < PhyConfig->nMaxFrameError)
		{
			generateLargeInfoBlock1D(
				PhyConfig,
				idxComponentCodewordInBlock,
				idxLowRateCodeword,
				idxHighRateCodeword,
				isLowRateCodeword,
				uInfoBlock,
				originalInfoBlock,
				xEncodedBlock,
				A,
				Ac,
				A_Ac,
				A_V,
				A_Ac_V,
				Ac_V,
				kTimeAll,
				nBadChannelsMIMOBlock,
				nComponentCodewordInBlock,
				nWholeCodewordInBlock,
				nLenBlock,
				nLenInfoBlock,
				nInfo,
				codeLength,
				nBadChannels,
				gen,
				distrib
			);
			for (int i = 0; i < nLenBlock; i++) { xEncodedBlockNoIntrlv[i] = xEncodedBlock[i]; }
			// randIntrlvWholeBlock(PhyConfig, xEncodedBlock, intrlvPatternTime, intrlvPatternSpace, intrlvPatternSpaceSep, nWholeCodewordInBlock, codeLength);
			randIntrlvWholeBlock1D(PhyConfig, xEncodedBlock, intrlvPattern1D, nWholeCodewordInBlock, codeLength);
			// for (int i = 0; i < nLenBlock; i++) { std::cout << xEncodedBlock[i] << std::endl; }
			codewordToAntennaMapping(PhyConfig, idxComponentCodewordInBlock, xTransmittedBlock, xEncodedBlock, nComponentCodewordInBlock);
			// for (int i = 0; i < 16; i++) { std::cout << idxComponentCodewordInBlock[i] << std::endl; }
			/*for (int i = 0; i < 512; i++) { std::cout << xEncodedBlock[i] << std::endl; }
			for (int i = 0; i < 512; i++) { std::cout << xTransmittedBlock[i] << std::endl; }*/
			modulationBlock(PhyConfig, xTransmittedBlock, xModulationBlock, nComponentCodewordInBlock);
			// generateChannelsBlock(PhyConfig, channelMatrixArray, nMIMOInBlock);
			if ((PhyConfig->use_platform_channel) || (PhyConfig->use_elaa_channel)) {
				generateChannelsBlock(PhyConfig, channelMatrixArray, nMIMOInBlock, channelMatrixArrayPlatform, gen, distrib);
			}
			else {
				generateChannelsBlock(PhyConfig, channelMatrixArray, nMIMOInBlock);
			}
			if (PhyConfig->use_elaa_channel) {
				transmissionMIMOBlockNoPrecode(PhyConfig, channelMatrixArray, nMIMOInBlock, nSymbols, lenMIMOBlock, xModulationBlock, llrBlock, nBadChannelsMIMOBlock, noiseVariance);
			}
			else {
				transmissionMIMOBlock(PhyConfig, channelMatrixArray, nMIMOInBlock, nSymbols, lenMIMOBlock, xModulationBlock, llrBlock, nBadChannelsMIMOBlock, noiseVariance);
			}
			// transmissionMIMOBlock(PhyConfig, channelMatrixArray, nMIMOInBlock, nSymbols, lenMIMOBlock, xModulationBlock, llrBlock, nBadChannelsMIMOBlock, noiseVariance);
			deMapDeIntrlvBlock1D(PhyConfig, llrBlock, llrDecoderBlock, idxComponentCodewordInBlock, intrlvPattern1D, nComponentCodewordInBlock, nWholeCodewordInBlock);
			// deMapDeIntrlvBlock(PhyConfig, llrBlock, llrDecoderBlock, idxComponentCodewordInBlock, intrlvPatternTime, intrlvPatternSpace, intrlvPatternSpaceSep, nComponentCodewordInBlock);
			for (int i = 0; i < nLenBlock; i++) { upLLRBlock[i] = llrDecoderBlock[i]; }
			// Decoding
			int Nmax = PhyConfig->nN;
# pragma omp parallel for
			for (int i_sub_frame = 0; i_sub_frame < nWholeCodewordInBlock; i_sub_frame++) {
				float* llr_1D = new float[PhyConfig->nN];
				unsigned char* crc_check_cwd = new unsigned char[PhyConfig->nK];
				int* reencoded_cwd = new int[PhyConfig->nN];
				for (int i = 0; i < PhyConfig->nN; i++) { llr_1D[i] = llrDecoderBlock[i_sub_frame * codeLength + i]; }
				int Counth = 0, Countinfoh = 0;
				float Qsumunh = 0;
				int** sum = new int* [PhyConfig->L];
				//for (int i = 0; i < L; i++)  sum[i] = new  int[Nmax];
				int** u1 = new int* [PhyConfig->L];
				for (int i = 0; i < PhyConfig->L; i++)  u1[i] = new _MM_ALIGN16 int[Nmax];
				for (int i = 0; i < PhyConfig->L; i++)  sum[i] = new _MM_ALIGN16 int[Nmax];
				//for (int i = 0; i < L; i++)  u1[i] = new  int[Nmax];

				float* LLR_in = new float[PhyConfig->L];
				float* W = new float[2 * PhyConfig->L];
				int* SCL_index = new int[PhyConfig->L];
				int* better = new int[PhyConfig->L];
				int* worse = new int[PhyConfig->L];
				int* indexpm = new int[PhyConfig->L];
				float* PM = new float[PhyConfig->L];
				float** LLR = new float* [PhyConfig->L];
				for (int i = 0; i < PhyConfig->L; i++)  LLR[i] = new _MM_ALIGN16 float[2 * Nmax + 2];
				// for (int i = 0; i < L; i++)  LLR[i] = new float[2 * Nmax + 2];
				////////////////////////////////
				for (int j = 0; j < PhyConfig->nN; j++)
					for (int k = 0; k < PhyConfig->L; k++)
						LLR[k][j] = llr_1D[j];
				// LLR[k][j] = upLLRBlock[frame_llr_offset + j * PhyConfig->nNh + decodingi];
		//cout << LLR[0][j] << " ";

	//cout << endl;
				int index = 0;
				for (int k = 0; k < PhyConfig->L; k++) { PM[k] = 0; indexpm[k] = k; }
				decode_list(LLR, sum, (int)(log(PhyConfig->nN) / log(2)), A_Ac.get(), 0, PhyConfig->nN, PhyConfig->nK, PM, L, (int)(log(PhyConfig->L) / log(2)), LLR_in, W, SCL_index, better, worse, &Counth, &Countinfoh, &Qsumunh);
				for (int i = 0; i < PhyConfig->L - 1; i++)
					for (int j = i + 1; j < PhyConfig->L; j++)
						if (PM[indexpm[i]] < PM[indexpm[j]])
						{
							int ttt = indexpm[i];
							indexpm[i] = indexpm[j];
							indexpm[j] = ttt;
						}
				index = indexpm[0];
				bool isPass = 0;
				if (PhyConfig->nCRCL > 0) {
					for (int i = 0; i < PhyConfig->L; i++) {
						for (int kk = 0; kk < PhyConfig->nN; kk++) { reencoded_cwd[kk] = 0; }
						PolarEncode_xor(reencoded_cwd, sum[indexpm[i]], PhyConfig->nN);
						// for (int j = 0; j < PhyConfig->nK; j++) { crc_check_cwd[j] = (unsigned char)reencoded_cwd[A[j]]; }
						int A_cnt = 0;
						for (int j = 0; j < PhyConfig->nN; j++) {
							if (A_Ac[j] == 1) { 
								crc_check_cwd[A_cnt] = (unsigned char)reencoded_cwd[j];
								// std::cout << (int)crc_check_cwd[A_cnt] << originalInfoBlock[i*codeLength+A_cnt]  << std::endl;
								A_cnt++;
							}
							
						}
						if (rx_check_crc(crc_check_cwd, PhyConfig->nK, PhyConfig->nCRCL)) {
							index = indexpm[i];
							isPass = 1;
							// std::cout << index << std::endl;
							// if (i > 0) std::cout << "FOUND" << std::endl;
							break;
						}
					}
				}
				PolarEncode_xor(uInfoBlockEsti.get() + i_sub_frame * codeLength, sum[index], PhyConfig->nN);
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
				delete[] llr_1D;
				delete[] crc_check_cwd;
				delete[] reencoded_cwd;
			}

			for (int i_sub_frame = 0; i_sub_frame < nWholeCodewordInBlock; i_sub_frame++) {
				int singleFrameError = 0;
				for (int i = 0; i < nInfo; i++) { if (uInfoBlockEsti[i_sub_frame * codeLength + A[i]] != uInfoBlock[i_sub_frame * codeLength + A[i]]) { singleFrameError++; } }
				if (singleFrameError > 0) { nBitError += singleFrameError; nFrameError++; }
				nowFrame++;
			}
			if ((nowFrame % nWholeCodewordInBlock) == 0) {
				printf("\r<%d> SNR= %.2f FER= %d / %d = %f BER= %d / %d = %lf",
					nowFrame, nowEbN0, nFrameError, nowFrame, (double)nFrameError / nowFrame, nBitError, (nInfo * nowFrame), (double)nBitError / nInfo / nowFrame);
				fflush(stdout);
			}
		}
		//Decoding_Duration = Decoding_Duration * 100 / CLOCKS_PER_SEC / frame;
		//Decoding_Duration = Decoding_Duration / CLOCKS_PER_SEC / Max_frame;
		std::cout << "  " << nowEbN0 << "  " << (double)nFrameError / nowFrame << "  " << (double)nBitError / nInfo / nowFrame << std::endl;
		FER[iEbN0] = (double)nFrameError / nowFrame;
	}

	iEbN0++;
	/*std::cout << "deleted GGh" << std::endl;
	for (int i = 0; i < PhyConfig->nNh; i++) { delete[] GGh[i]; } delete[] GGh;
	
	for (int i = 0; i < PhyConfig->nNv; i++) { delete[] GGv[i]; } delete[] GGv;
	std::cout << "deleted GGv" << std::endl;*/
	for (int i = 0; i < PhyConfig->nN; i++) { delete[] GG[i]; } delete[] GG;
	// std::cout << "deleted GG" << std::endl;
}
