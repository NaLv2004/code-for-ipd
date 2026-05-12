#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <fstream>
#include <vector>
#include <bitset>
#include <chrono> 
#include <Eigen\Dense>
#include <Eigen\Core>
#include <unsupported\Eigen\KroneckerProduct>
#include "crc.h"
#include "Polar_Encoder.h"
#include "Polar_Decoder.h"
#include "Polar_Function.h"
#include "Polar_Constuction.h"
#include "setting.h"
#include <random>
using namespace Eigen;
using namespace std;
using pInt = std::unique_ptr<int[]>;
using pFloat = std::unique_ptr<float[]>;
using pChar = std::unique_ptr<unsigned char[]>;
using pCfloat = std::unique_ptr<std::complex<float>[]>;

struct sPhyConfig {
	int nK = 687;
	int nN = 1024;
	int nNv = 32;
	int nNh = 32;
	int nKv = 26;
	int nKh = 26;
	int nNt = 32;
	int nNr = 64;
	int nModType = 4;
	int nStream = 32;
	int nCRCL = 0;
	int nDecodingIter = 4;
	float minSNR = 7.0;
	float maxSNR = 9.0;
	float stepSNR = 0.5;
	int nTimeConstructions;
	int nLowRateCodewords = 6;
	int L = 8;
	int* pKTime = new int[2];
	int nMaxFrameError = 100;
	bool bFlagInterleave = 1;
	bool bFlagInterleaveTemporal = 1;
	bool bFlagInterleaveSpatial = 1;
	bool bMIMO = 1;
	bool bIrregularCoding = 1;
	std::string polarConh = "RM";
	std::string polarConV = "RM";
	std::string powerAllocationScheme =  "SEPARATE";
	std::string polarCon = "5G";
	float powerAllocationRatio = 0.4;
	float* alpha = nullptr;
	bool bEarlyStop = 0;
	int frameThres = 1000;
	float ferThres = 1E-3;
	int use_platform_channel = USE_PLATFORM_CHANNEL;
	int use_elaa_channel = USE_ELAA_CHANNEL;
	int n_frequency_samples = N_FREQUENCY_SAMPLES;
	//pFloat alpha(new float[nDecodingIter]);
};

struct soptParams {
	int nNh = 32;
	int nNv = 32;
	int nKh = 26;
	int nKv = 26;
	float codeRate = 1;
	int border = 6;
	int* pKTime = new int[2];
	int nStream = 32;
	float transmissionRate = 32;
	float errorRate = 0.0;
};

int gcd(int a, int b);
int lcm(int a, int b);
// struct sPhyConfig;
void polar_code_construction_new(sPhyConfig* PhyConfig, float noiseVariance, pInt& A_H, pInt& Ac_H, pInt& A_Ac_H, pInt& A_V, pInt& Ac_V, pInt& A_Ac_V, pInt& kTimeAll, pInt& isLowRateCodeword,
	int** GGh, int** GGv);
//void generateLargeInfoBlock(sPhyConfig* PhyConfig, pInt& idxComponentCodewordInBlock, pInt& idxLowRateCodeword, pInt& idxHighRateCodeword, pInt& isLowRateCodeword, pInt& uInfoBlock, pInt& originalInfoBlock, pInt& xEncodedBlock,
//	pInt& A_H, pInt& Ac_H, pInt& A_Ac_H, pInt& A_V, pInt& A_Ac_V, pInt& Ac_V, pInt& kTimeAll, int nComponentInBlock, int nWholeCodewordInBlock, int nLenBlock, int nLenInfoBlock, int nInfo, int codeLength, int& nBadChannels,
//	std::mt19937& gen, std::uniform_int_distribution<> distrib
//);
void generateLargeInfoBlock(sPhyConfig* PhyConfig, pInt& idxComponentCodewordInBlock, pInt& idxLowRateCodeword, pInt& idxHighRateCodeword, pInt& isLowRateCodeword, pInt& uInfoBlock, pInt& originalInfoBlock, pInt& xEncodedBlock,
	pInt& A_H, pInt& Ac_H, pInt& A_Ac_H, pInt& A_V, pInt& A_Ac_V, pInt& Ac_V, pInt& kTimeAll, pInt& nBadChannelsMIMOBlock, int nComponentInBlock, int nWholeCodewordInBlock, int nLenBlock, int nLenInfoBlock, int nInfo, int codeLength, int& nBadChannels,
	std::mt19937& gen, std::uniform_int_distribution<> distrib
);
// void randIntrlvWholeBlock(sPhyConfig* PhyConfig, pInt& xEncodedBlock, vector<vector<int>>& intrlvPatternTime, int nWholeCodewordInBlock, int codeLength);
void codewordToAntennaMapping(sPhyConfig* PhyConfig, pInt& idxComponentCodewordInBlock, pInt& xTransmittedBlock, pInt& xEncodedBlock, int nComponentCodewordInBlock);
void modulationBlock(sPhyConfig* PhyConfig, pInt& xTransmittedBlock, pCfloat& xModulationBlock, int nComponentCodewordInBlock);
void generateChannelsBlock(sPhyConfig* PhyConfig, std::vector<MatrixXcf>& channelMatrixArray, int nMIMOInBlock);
// void transmissionMIMOBlock(sPhyConfig* PhyConfig, std::vector<MatrixXcf>& channelMatrixArray, int nMIMOInBlock, int nSymbols, int lenMIMOBlock, pCfloat& xTransmittedBlock, pFloat& llrBlock, float noiseVariance);
void transmissionMIMOBlock(sPhyConfig* PhyConfig, std::vector<MatrixXcf>& channelMatrixArray, int nMIMOInBlock, int nSymbols, int lenMIMOBlock, pCfloat& xTransmittedBlock, pFloat& llrBlock, pInt& nBadChannelsMIMOBlock, float noiseVariance);
// void deMapDeIntrlvBlock(sPhyConfig* PhyConfig, pFloat& llrBlock, pFloat& llrDecoderBlock, pInt& idxComponentCodewordInBlock, vector<vector<int>> intrlvPatternTime, int& nComponentCodewordInBlock);
void system2D_new(sPhyConfig* PhyConfig, float* BER, float* FER);
void deMapDeIntrlvBlock(sPhyConfig* PhyConfig, pFloat& llrBlock, pFloat& llrDecoderBlock, pInt& idxComponentCodewordInBlock, vector<vector<int>> intrlvPatternTime, vector<vector<int>>& intrlvPatternSpace, vector<vector<int>>& intrlvPatternSpaceSep, int& nComponentCodewordInBlock);
void randIntrlvWholeBlock(sPhyConfig* PhyConfig, pInt& xEncodedBlock, vector<vector<int>>& intrlvPatternTime, vector<vector<int>>& intrlvPatternSpace, vector<vector<int>> intrlvPatternSpaceSep, int nWholeCodewordInBlock, int codeLength);
void generateIntrlvPatternIrrNew(sPhyConfig* PhyConfig, std::vector<int>& pattern, int len, pInt& isLowRateCodeword);
void setOptParam(sPhyConfig* PhyConfig, soptParams* OptParams, float errorRate);
float SNR2EbN0(float SNR, sPhyConfig* PhyConfig);
void polar_code_construction_new_1D(sPhyConfig* PhyConfig, float noiseVariance, pInt& A, pInt& A_Ac, pInt& Ac, int** GG);
void generateLargeInfoBlock1D(sPhyConfig* PhyConfig, pInt& idxComponentCodewordInBlock, pInt& idxLowRateCodeword, pInt& idxHighRateCodeword, pInt& isLowRateCodeword, pInt& uInfoBlock, pInt& originalInfoBlock, pInt& xEncodedBlock,
	pInt& A_H, pInt& Ac_H, pInt& A_Ac_H, pInt& A_V, pInt& A_Ac_V, pInt& Ac_V, pInt& kTimeAll, pInt& nBadChannelsMIMOBlock, int nComponentInBlock, int nWholeCodewordInBlock, int nLenBlock, int nLenInfoBlock, int nInfo, int codeLength, int& nBadChannels,
	std::mt19937& gen, std::uniform_int_distribution<> distrib
);
void randIntrlvWholeBlock1D(sPhyConfig* PhyConfig, pInt& xEncodedBlock, vector<int> intrlvPattern1D, int nWholeCodewordInBlock, int codeLength);
void deMapDeIntrlvBlock1D(sPhyConfig* PhyConfig, pFloat& llrBlock, pFloat& llrDecoderBlock, pInt& idxComponentCodewordInBlock, vector<int>& intrlvPattern1D, int& nComponentCodewordInBlock, int nWholeCodewordInBlock);
void system1D_new(sPhyConfig* PhyConfig, float* BER, float* FER);
float SNR2EbN0_1D(float SNR, sPhyConfig* PhyConfig);
void generateChannelsBlock(sPhyConfig* PhyConfig, std::vector<MatrixXcf>& channelMatrixArray, int nMIMOInBlock, std::vector<MatrixXcf>& channelMatrixArrayPlatform, std::mt19937& gen, std::uniform_int_distribution<> distrib);
void MMSEDetectionCfloat(pCfloat& y_recv, Eigen::MatrixXcf H, double noiseVariance, sPhyConfig* PhyConfig);
void transmissionMIMOBlockNoPrecode(sPhyConfig* PhyConfig, std::vector<MatrixXcf>& channelMatrixArray, int nMIMOInBlock, int nSymbols, int lenMIMOBlock, pCfloat& xTransmittedBlock, pFloat& llrBlock, pInt& nBadChannelsMIMOBlock, float noiseVariance);
void mimoChannelPowerAlloc(sPhyConfig* PhyConfig, std::vector<Eigen::MatrixXcf>& channelMatrixArray,
	int nComponentMaxD, float minD, float maxD, float P);
//void mimoChannelPowerAllocSeparate(sPhyConfig* PhyConfig, std::vector<Eigen::MatrixXcf>& channelMatrixArray,
//	int nComponentMaxD, float minD, float maxD, float P);
void mimoChannelPowerAllocSeparate(sPhyConfig* PhyConfig, std::vector<Eigen::MatrixXcf>& channelMatrixArray,
	int nComponentMaxD, float minD, float maxD, float P, float noiseVariance);