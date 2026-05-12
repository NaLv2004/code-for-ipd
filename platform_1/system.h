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
using namespace Eigen;
using namespace std;
void system1D(double* BER_array, double* FER_array);
void system2D(double* BER_array, double* FER_array);
void system1D_LSD(double* BER_array, double* FER_array);
void system2D_LSD(double* BER_array, double* FER_array);
void system2D_LSD_epdet(double* BER_array, double* FER_array);
void system2D_LSD_epdet_sys(double* BER_array, double* FER_array);
void system2D_LSD_partial_iter(double* BER_array, double* FER_array);
void system2D_LSD_sys(double* BER_array, double* FER_array);
void system2D_LSD_doublestage(double* BER_array, double* FER_array);
void system1D_LSD_MIMO(double* BER_array, double* FER_array);
void system1D_LSD_MIMO_SVD(double* BER_array, double* FER_array);
void system2D_JDD(double* BER_array, double* FER_array);
void system2D_JDD_copy(double* BER_array, double* FER_array);
void system2D_LSD_epdet_with_CRC(double* BER_array, double* FER_array);
void system2D_LSD_epdet_with_CRC_vertical(double* BER_array, double* FER_array);
void system2D_LSD_epdet_irregular(double* BER_array, double* FER_array);
void system2D_LSD_epdet_irregular_SVD_precode(double* BER_array, double* FER_array);
MatrixXf Fn(int N);