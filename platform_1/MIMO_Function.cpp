#include"MIMO_Function.h"
// #include"modulation.cpp"
#include<cmath>
#include<algorithm>
#include<vector>
#include<random>
#include<iterator>
#include<iomanip>
#include<omp.h>
#include"setting.h"
#include<chrono>
#include<mkl.h>
//using namespace std;

//#define cf_t	MKL_Complex8
#define PI 3.1415926

float table_qpsk_det[2] = { -0.707107f, 0.707107f };
float table_16qam_det[4] = { -0.316228f, -0.948683f,  0.316228f,  0.948683f };
float table_64qam_det[8] = { -0.462910f, -0.154303f, -0.771517f, -1.08012f,	0.462910f,  0.154303f,  0.771517f,  1.08012f };
float table_256qam_det[16] =
{ -0.383482f, -0.536875f, -0.230089f, -0.07669f,
-0.843661f, -0.690268f, -0.997054f, -1.15044f,
0.383482f, 0.536875f, 0.230089f, 0.07669f,
0.843661f, 0.690268f, 0.997054f, 1.15044f };

float s_bits_r2_det[2][1] = { {0}, {1} };
float s_bits_i2_det[2][1] = { {0}, {1} };
float s_bits_r4_det[4][2] = { {0,0},{0,1},{1,0},{1,1} };
float s_bits_i4_det[4][2] = { {0,0},{0,1},{1,0},{1,1} };
float s_bits_r6_det[8][3] = { {0,0,0},{0,0,1},{0,1,0},{0,1,1},{1,0,0},{1,0,1},{1,1,0},{1,1,1} };
float s_bits_i6_det[8][3] = { {0,0,0},{0,0,1},{0,1,0},{0,1,1},{1,0,0},{1,0,1},{1,1,0},{1,1,1} };
float s_bits_r8_det[16][4] = { {0,0,0,0},{0,0,0,1},{0,0,1,0},{0,0,1,1},{0,1,0,0},{0,1,0,1},{0,1,1,0},{0,1,1,1},{1,0,0,0},{1,0,0,1},{1,0,1,0},{1,0,1,1},{1,1,0,0},{1,1,0,1},{1,1,1,0},{1,1,1,1} };
float s_bits_i8_det[16][4] = { {0,0,0,0},{0,0,0,1},{0,0,1,0},{0,0,1,1},{0,1,0,0},{0,1,0,1},{0,1,1,0},{0,1,1,1},{1,0,0,0},{1,0,0,1},{1,0,1,0},{1,0,1,1},{1,1,0,0},{1,1,0,1},{1,1,1,0},{1,1,1,1} };

double sampleNormal_mimo() {
    double p = ((double)rand() / (RAND_MAX)) * 2 - 1;
    double v = ((double)rand() / (RAND_MAX)) * 2 - 1;
    double r = p * p + v * v;
    if (r == 0 || r > 1) return sampleNormal_mimo();
    double c = sqrt(-2 * log(r) / r);
    return p * c;
}

void inverseMatrix(float* pSrc, int dim)
{
    _int64* ivpv = new _int64[sizeof(int) * dim];
    float* pSrcBak = new float[dim * dim];  // LAPACKE_sgesv会覆盖A矩阵，因而将pSrc备份
    memset(pSrcBak, 0, sizeof(float) * dim * dim);
    for (int i = 0; i < dim; ++i)
    {
        // LAPACKE_sgesv函数计算AX=B，当B为单位矩阵时，X为inv(A)
        pSrcBak[i * (dim + 1)] = 1.0f;
    }

    size_t flag = LAPACKE_sgesv(LAPACK_ROW_MAJOR, dim, dim, pSrc, dim, ivpv, pSrcBak, dim);
    // 调用LAPACKE_sgesv后，会将inv(A)覆盖到X（即pDst）
    if (flag)
    {
        printf("Inverse is failed!\n");
    }

    memcpy(pSrc, pSrcBak, sizeof(float) * dim * dim);


    delete[] ivpv;
    ivpv = nullptr;

    delete[] pSrcBak;
    pSrcBak = nullptr;
}

// Interleaving and Deinterleaving
void generateIntrlvPattern(std::vector<int>& pattern, int len)
{
    for (int i = 0; i < len; i++) { pattern[i] = i; }
    std::mt19937 gen(std::random_device{}());
    std::shuffle(pattern.begin(), pattern.end(), gen);
}

void generateIntrlvPatternInfoMapping(std::vector<int>& pattern, int nv, int kv, int* Av, int* Acv)
{
    for (int i = 0; i < kv; i++) { pattern[i] = Av[i]; }
    for (int i = kv; i < nv; i++) { pattern[i] = Acv[i - kv]; }
}

void randIntrlv(std::vector<int>& pattern, int* t_bits, int len)
{
    int* copy_bits = new int[len];
    for (int i = 0; i < len; i++) { copy_bits[i] = t_bits[i]; }
    for (int i = 0; i < len; i++)
    {
        t_bits[i] = copy_bits[pattern[i]];
    }
    delete[] copy_bits;
}

void randIntrlv(std::vector<int>& pattern, float* t_bits, int len)
{
    float* copy_bits = new float[len];
    for (int i = 0; i < len; i++) { copy_bits[i] = t_bits[i]; }
    for (int i = 0; i < len; i++)
    {
        t_bits[i] = copy_bits[pattern[i]];
    }
    delete[] copy_bits;
}

void randDeintrlv(std::vector<int>& pattern, float* rx_llr, int len)
{
    float* copy_llr = new float[len];
    for (int i = 0; i < len; i++) { copy_llr[i] = rx_llr[i]; }
    for (int i = 0; i < len; i++) { rx_llr[pattern[i]] = copy_llr[i]; }
    delete[] copy_llr;
}

void generateIntrlvPatternIrregular(std::vector<int>& pattern, int len, int* bad_channel_index, int* good_channel_index, int cnt_bad_channels, int* use_short_codeword)
{
    int* lr_codeword_index = new int[cnt_bad_channels];
    int* hr_codeword_index = new int[len - cnt_bad_channels];
    int cnt_lr = 0;
    int cnt_hr = 0;
    for (int i = 0; i < len; i++) {
        if (use_short_codeword[i] == 1) { lr_codeword_index[cnt_lr] = i; cnt_lr++; }
        else { hr_codeword_index[cnt_hr] = i; cnt_hr++; }
    }
    cnt_lr = 0;
    cnt_hr = 0;
    for (int i = 0; i < cnt_bad_channels; i++)
    {
        pattern[bad_channel_index[i]] = lr_codeword_index[cnt_lr];
        cnt_lr++;
    }
    for (int i = 0; i < len-cnt_bad_channels; i++)
    {
        pattern[good_channel_index[i]] = hr_codeword_index[cnt_hr];
        cnt_hr++;
    }
    delete[] lr_codeword_index;
    delete[] hr_codeword_index;
}


void generateIntrlvPatternIrregularArray(std::vector<int>& pattern, int len, int* bad_channel_index, int* good_channel_index, int cnt_bad_channels, int* use_short_codeword)
{
    vector<int> intrlv_pattern_lr;
    vector<int> intrlv_pattern_hr;
    intrlv_pattern_lr.resize(cnt_bad_channels);
    intrlv_pattern_hr.resize(len - cnt_bad_channels);
    generateIntrlvPattern(intrlv_pattern_lr, cnt_bad_channels);
    generateIntrlvPattern(intrlv_pattern_hr, len - cnt_bad_channels);
    int* lr_codeword_index = new int[cnt_bad_channels];
    int* hr_codeword_index = new int[len - cnt_bad_channels];
    int cnt_lr = 0;
    int cnt_hr = 0;
    for (int i = 0; i < len; i++) {
        if (use_short_codeword[i] == 1) { lr_codeword_index[cnt_lr] = i; cnt_lr++; }
        else { hr_codeword_index[cnt_hr] = i; cnt_hr++; }
    }
    randIntrlv(intrlv_pattern_lr, lr_codeword_index, cnt_bad_channels);
    randIntrlv(intrlv_pattern_hr, hr_codeword_index, len - cnt_bad_channels);
    cnt_lr = 0;
    cnt_hr = 0;
    for (int i = 0; i < cnt_bad_channels; i++)
    {
        pattern[bad_channel_index[i]] = lr_codeword_index[cnt_lr];
        cnt_lr++;
    }
    for (int i = 0; i < len - cnt_bad_channels; i++)
    {
        pattern[good_channel_index[i]] = hr_codeword_index[cnt_hr];
        cnt_hr++;
    }
}



//void randIntrlvBlock(int** bits, int len_block, int height_block, int len_subblock, int interleave_rows)
//{
//    int n_subblock = len_block / len_subblock;
//    int len_array = len_subblock * height_block;
//    int* bits_before_interleave = new int[len_array];
//    int* bits_interleaved = new int[len_array];
//    for (int i_sub = 0; i_sub < n_subblock; i_sub++)
//    {
//
//    }
//    delete[] bits_before_interleave;
//    delete[] bits_interleaved;
//}

void randIntrlvBlock(int** bits, int len_block, int height_block, int len_subblock, int interleave_rows)
{
    int n_subblock = len_block / len_subblock;
    int len_array = len_subblock * height_block;
    int* bits_before_interleave = new int[len_array];
    int* bits_interleaved = new int[len_array];
    for (int i_sub = 0; i_sub < n_subblock; i_sub++)
    {
        // 提取当前子块到一维数组（按行优先）
        int index = 0;
        for (int r = 0; r < height_block; r++) {
            for (int c = 0; c < len_subblock; c++) {
                int original_col = i_sub * len_subblock + c;
                bits_before_interleave[index++] = bits[r][original_col];
            }
        }

        // 执行交织操作（按列读取interleave_rows行结构）
        int cols = len_array / interleave_rows;
        for (int j = 0; j < cols; j++) {
            for (int i = 0; i < interleave_rows; i++) {
                int original_idx = i * cols + j;
                int interleaved_idx = j * interleave_rows + i;
                bits_interleaved[interleaved_idx] = bits_before_interleave[original_idx];
            }
        }

        // 将交织后的数据写回原矩阵
        index = 0;
        for (int r = 0; r < height_block; r++) {
            for (int c = 0; c < len_subblock; c++) {
                bits[r][i_sub * len_subblock + c] = bits_interleaved[index++];
            }
        }
    }
    delete[] bits_before_interleave;
    delete[] bits_interleaved;
}


void randDeIntrlvBlock(float** bits, int len_block, int height_block, int len_subblock, int interleave_rows)
{
    int n_subblock = len_block / len_subblock;
    int len_array = len_subblock * height_block;
    float* bits_interleaved = new float[len_array];
    float* bits_deinterleaved = new float[len_array];

    for (int i_sub = 0; i_sub < n_subblock; i_sub++)
    {
        // 从当前子块提取交织后的数据（按行优先存储）
        int index = 0;
        for (int r = 0; r < height_block; r++) {
            for (int c = 0; c < len_subblock; c++) {
                int current_col = i_sub * len_subblock + c;
                bits_interleaved[index++] = bits[r][current_col];
            }
        }

        // 执行解交织操作
        int cols = len_array / interleave_rows;
        for (int i = 0; i < interleave_rows; i++) {
            for (int j = 0; j < cols; j++) {
                int src_idx = j * interleave_rows + i;
                int dst_idx = i * cols + j;
                bits_deinterleaved[dst_idx] = bits_interleaved[src_idx];
            }
        }

        // 将解交织数据写回原矩阵
        index = 0;
        for (int r = 0; r < height_block; r++) {
            for (int c = 0; c < len_subblock; c++) {
                bits[r][i_sub * len_subblock + c] = bits_deinterleaved[index++];
            }
        }
    }

    delete[] bits_interleaved;
    delete[] bits_deinterleaved;
}


void randIntrlvPartial(std::vector<int>& pattern, int* bits, int* info_pos, int N, int K_info)
{
    int* bits_interleaved = new int[N];
    for (int i = 0; i < K_info; i++)
    {
        bits_interleaved[info_pos[pattern[i]]] = bits[info_pos[i]];
    }
    for (int i = 0; i < K_info; i++)
    {
        bits[info_pos[i]] = bits_interleaved[info_pos[i]];
    }
    delete[] bits_interleaved;
}


void randIntrlvPartial(std::vector<int>& pattern, float* bits, int* info_pos, int N, int K_info)
{
    float* bits_interleaved = new float[N];
    for (int i = 0; i < K_info; i++)
    {
        bits_interleaved[info_pos[pattern[i]]] = bits[info_pos[i]];
    }
    for (int i = 0; i < K_info; i++)
    {
        bits[info_pos[i]] = bits_interleaved[info_pos[i]];
    }
    delete[] bits_interleaved;
}


void randDeIntrlvPartial(std::vector<int>& pattern, float* bits, int* info_pos, int N, int K_info)
{
    float* bits_deinterleaved = new float[N];
    for (int i = 0; i < K_info; i++)
    {
        bits_deinterleaved[info_pos[i]] = bits[info_pos[pattern[i]]];
    }
    for (int i = 0; i < K_info; i++)
    {
        bits[info_pos[i]] = bits_deinterleaved[info_pos[i]];
    }
    delete[] bits_deinterleaved;
}


// 生成信道矩阵
MatrixXcf generateChannelMatrix(int Nr, int Nt) {
    // 创建一个 Nr x Nt 的复数矩阵
    MatrixXcf H(Nr, Nt);

    // 填充信道矩阵 H
    for (int i = 0; i < Nr; i++) {
        for (int j = 0; j < Nt; j++) {
            // 从 sampleNormal() 生成实部和虚部
            float real = sqrt(0.5) * sampleNormal_mimo();  // 实部
            float imag = sqrt(0.5) * sampleNormal_mimo();  // 虚部

            // 将复数值赋给信道矩阵
            H(i, j) = std::complex<float>(real, imag);
            
        }
    }
    return H;
}


MatrixXcf generateChannelMatrixNonuniform(int Nr, int Nt, float** channel_scaling_factor) {
    // 创建一个 Nr x Nt 的复数矩阵
    MatrixXcf H(Nr, Nt);

    // 填充信道矩阵 H
    for (int i = 0; i < Nr; i++) {
        for (int j = 0; j < Nt; j++) {
            // 从 sampleNormal() 生成实部和虚部
            float real = sqrt(0.5) * sampleNormal_mimo()* channel_scaling_factor[i][j];  // 实部
            float imag = sqrt(0.5) * sampleNormal_mimo()* channel_scaling_factor[i][j];  // 虚部

            // 将复数值赋给信道矩阵
            H(i, j) = std::complex<float>(real, imag);

        }
    }
    return H;
}


void getChannelScalingFactor(int Nr, int Nt, int diff_gain, int* channel_index, float** channel_scale_factor, int predefined) {
    float tmp = diff_gain / 20.0;
    float diff_linear = pow(10.0, tmp);
    float sum = 0;
    for (int i = 0; i < Nr; i++)
    {
        for (int j = 0; j < Nt; j++)
        {
            if (!predefined)
            {
                if (channel_index[j] == 1) { channel_scale_factor[i][j] = 1.0; }
                else { channel_scale_factor[i][j] = diff_linear; }
                sum += channel_scale_factor[i][j];
            }
            else
            {
                if (channel_index_predefined[j] == 1) { channel_scale_factor[i][j] = 1.0; }
                else { channel_scale_factor[i][j] = diff_linear; }
                sum += channel_scale_factor[i][j];
            }
        }
    }
    float norm_sum = (float)Nr * (float)Nt;
    float norm_factor = norm_sum / sum;
    for (int i = 0; i < Nr; i++) { for (int j = 0; j < Nt; j++) { channel_scale_factor[i][j] *= norm_factor; /*cout << channel_scale_factor[i][j] << "  ";*/} }
   /* for (int i = 0; i < Nr; i++)
    {
        cout << endl;
        for (int j = 0; j < Nt; j++) { cout << channel_scale_factor[i][j] << "  "; }
    }
    cout << endl;*/
}

// 将复数矩阵 H 转换为实数域
void convertHToReal(const MatrixXcf& H, Matrix<float, 2 * Nr, 2 * Nt>& H_real) {
    int Nr = H.rows(); // 接收天线数量
    int Nt = H.cols(); // 发送天线数量,
    // 创建一个实数矩阵，行数为 2 * Nr，列数为 2 * Nt
   // MatrixXf H_real(2 * Nr, 2 * Nt);

    for (int i = 0; i < Nr; i++) {
        for (int j = 0; j < Nt; j++) {
            H_real(i, j) = H(i, j).real();
            H_real(i, j + Nt) = -H(i, j).imag();
            H_real(i + Nr, j) = H(i, j).imag();
            H_real(i + Nr, j + Nt) = H(i, j).real();
        }
    }
}

void convertHToReal(const MatrixXcf& H, MatrixXf& H_real) {
    int Nr = H.rows(); // 接收天线数量
    int Nt = H.cols(); // 发送天线数量,
    // 创建一个实数矩阵，行数为 2 * Nr，列数为 2 * Nt
   // MatrixXf H_real(2 * Nr, 2 * Nt);

    for (int i = 0; i < Nr; i++) {
        for (int j = 0; j < Nt; j++) {
            H_real(i, j) = H(i, j).real();
            H_real(i, j + Nt) = -H(i, j).imag();
            H_real(i + Nr, j) = H(i, j).imag();
            H_real(i + Nr, j + Nt) = H(i, j).real();
        }
    }
}

void AWGN(int Nr, int Nt, std::complex<float>* x, std::complex<float>* y, const MatrixXcf& H, double n0) {
    for (int i = 0; i < Nr; i++)
    {
        y[i] = 0;
        for (int j = 0; j < Nt; j++)
        {
            y[i] += x[j] * H(i, j);
        }
        float n_real = (float)sqrt(n0 / 2) * (float)sampleNormal_mimo();
        float n_imag = (float)sqrt(n0 / 2) * (float)sampleNormal_mimo();
        y[i] = y[i] + std::complex<float>(n_real, n_imag);
    }
}

void precodeSVD(int Nr, int Nt, std::complex<float>* x, std::complex<float>* x_precode_out, const MatrixXcf& V) {
    for (int i = 0; i < Nt; i++)
    {
        x_precode_out[i] = 0;
        for (int j = 0; j < Nt; j++)
        {
            x_precode_out[i] += x[j] * V(i, j);
        }
        /*float n_real = (float)sqrt(n0 / 2) * (float)sampleNormal_mimo();
        float n_imag = (float)sqrt(n0 / 2) * (float)sampleNormal_mimo();
        y[i] = y[i] + std::complex<float>(n_real, n_imag);*/
    }
}

void precodeSVD_preprocessing(int Nr, int Nt, std::complex<float>* y, std::complex<float>* y_preprocess_out, const MatrixXcf& U) {
    for (int i = 0; i < Nr; i++)
    {
        y_preprocess_out[i] = 0;
        for (int j = 0; j < Nr; j++)
        {
            std::complex<float> tmp = std::complex<float>(U(j, i).real(), -U(j, i).imag());
            y_preprocess_out[i] += y[j] * tmp;
        }
        /*float n_real = (float)sqrt(n0 / 2) * (float)sampleNormal_mimo();
        float n_imag = (float)sqrt(n0 / 2) * (float)sampleNormal_mimo();
        y[i] = y[i] + std::complex<float>(n_real, n_imag);*/
    }
}

void convertYToReal(int Nr, std::complex<float>* y, float* y_real) {
    for (int i = 0; i < Nr; i++)
    {
        y_real[i] = y[i].real();
        y_real[i + Nr] = y[i].imag();
    }
}

void MMSEdetection(const MatrixXf H_real, const MatrixXf& identity, MatrixXf& received_signal, MatrixXf& mmse_matrix, MatrixXf& output_signal, MatrixXf& HTH, float* rx, float* y_hat, float sigma_n2) {
    for (int i = 0; i < 2 * Nr; i++)
        received_signal(i, 0) = rx[i];
    HTH = H_real.transpose() * H_real;
    mmse_matrix = (HTH + sigma_n2 * identity).inverse() * H_real.transpose();
    output_signal = mmse_matrix * received_signal;
    for (int i = 0; i < 2 * Nt; ++i) {
        y_hat[i] = output_signal(i, 0);
    }
}

void MMSEdetection(const MatrixXf H_real, const MatrixXf identity, MatrixXf received_signal, MatrixXf mmse_matrix, MatrixXf output_signal, MatrixXf HTH, float* rx, float* y_hat, float sigma_n2) {
    //auto det_start = std::chrono::system_clock::now();
    for (int i = 0; i < 2 * Nr; i++)
        received_signal(i, 0) = rx[i];
    HTH = H_real.transpose() * H_real;
    mmse_matrix = (HTH + sigma_n2 * identity).inverse() * H_real.transpose();
    //mmse_matrix = (HTH + sigma_n2 * identity); //* H_real.transpose();
    output_signal = mmse_matrix * received_signal;
    for (int i = 0; i < 2 * Nt; ++i) {
        y_hat[i] = output_signal(i, 0);
    }
    //auto det_over = std::chrono::system_clock::now();
    //auto det_duration = std::chrono::duration_cast<std::chrono::microseconds>(det_over-det_start);
    //cout << det_duration.count() << endl;
   /* for (int i = 0; i < 2*Nt; i++)
        for (int s = 0; s < 2*Nt; s++)
        {
            if (s == 0) std::cout << endl;
            std::cout << setw(10) << HTH(i, s);
        }
    cout << endl;*/
}

void MMSEdetection_exact(const MatrixXf H_real, const MatrixXf identity, MatrixXf received_signal, MatrixXf mmse_matrix, MatrixXf output_signal, MatrixXf HTH, float* rx, float* y_hat, float* mu_mmse, float sigma_n2) {
    //auto det_start = std::chrono::system_clock::now();
    //float* mu_mmse = new float[2 * Nt];
    for (int i = 0; i < 2 * Nr; i++)
        received_signal(i, 0) = rx[i];
    HTH = H_real.transpose() * H_real;
    mmse_matrix = (HTH + sigma_n2 * identity).inverse() * H_real.transpose();
    //mmse_matrix = (HTH + sigma_n2 * identity); //* H_real.transpose();
    output_signal = mmse_matrix * received_signal;
    for (int i = 0; i < 2 * Nt; ++i) {
        y_hat[i] = output_signal(i, 0);
    }
    //auto det_over = std::chrono::system_clock::now();
    //auto det_duration = std::chrono::duration_cast<std::chrono::microseconds>(det_over-det_start);
    //cout << det_duration.count() << endl;
   /* for (int i = 0; i < 2*Nt; i++)
        for (int s = 0; s < 2*Nt; s++)
        {
            if (s == 0) std::cout << endl;
            std::cout << setw(10) << HTH(i, s);
        }
    cout << endl;*/
    // mmse_matrix: Nt * Nr (evry row has Nr elements)
    // H: Nr * Nt (every col has Nr elements)
    for (int i = 0; i < 2*Nt; i++)
    {
        mu_mmse[i] = 0;
        for (int j = 0; j < 2 * Nr; j++)
        {
            mu_mmse[i] += mmse_matrix(i, j) * H_real(j, i);
        }
       // cout << mu_mmse[i] << endl;
       // y_hat[i] /= mu_mmse[i];
    }
    //delete[] mu_mmse;
}

void MMSEdetection_noInverse(const MatrixXf H_real, const MatrixXf identity, MatrixXf received_signal, MatrixXf mmse_matrix, MatrixXf output_signal, MatrixXf HTH, MatrixXf HTHIinv, float* rx, float* y_hat, float sigma_n2) {
    //auto det_start = std::chrono::system_clock::now();
    for (int i = 0; i < 2 * Nr; i++)
        received_signal(i, 0) = rx[i];
    /**HTH = H_real.transpose() * H_real; */
    mmse_matrix = HTHIinv * H_real.transpose();
    //mmse_matrix = (HTH + sigma_n2 * identity); //* H_real.transpose();
    output_signal = mmse_matrix * received_signal;
    for (int i = 0; i < 2 * Nt; ++i) {
        y_hat[i] = output_signal(i, 0);
    }
    //auto det_over = std::chrono::system_clock::now();
    //auto det_duration = std::chrono::duration_cast<std::chrono::microseconds>(det_over-det_start);
    //cout << det_duration.count() << endl;
}

void MMSEdetection(const Matrix<float, 2 * Nr, 2 * Nt>& H_real, const Matrix<float, 2 * Nt, 2 * Nt>& identity, Matrix<float, 2 * Nr, 1>& received_signal, Matrix<float, 2 * Nt, 2 * Nr>& mmse_matrix, Matrix<float, 2 * Nt, 1>& output_signal, Matrix<float, 2 * Nt, 2 * Nt>& HTH, float* rx, float* y_hat, float sigma_n2) {
    for (int i = 0; i < 2 * Nr; i++)
        received_signal(i, 0) = rx[i];
    HTH = H_real.transpose() * H_real;
    // mmse_matrix = (HTH + sigma_n2 * identity).inverse() * H_real.transpose();
    mmse_matrix = (HTH + sigma_n2 * identity).colPivHouseholderQr().solve(H_real.transpose());
    output_signal = mmse_matrix * received_signal;
    for (int i = 0; i < 2 * Nt; ++i) {
        y_hat[i] = output_signal(i, 0);
    }
}

//void MMSEdetection(const MatrixXf H_real, const MatrixXf& identity, MatrixXf& received_signal, MatrixXf& mmse_matrix, MatrixXf& output_signal, MatrixXf& HTH, float* rx, float* y_hat, float sigma_n2) {
//    for (int i = 0; i < 2 * Nr; i++)
//        received_signal(i, 0) = rx[i];
//    HTH = H_real.transpose() * H_real;
//    mmse_matrix = (HTH + sigma_n2 * identity).inverse() * H_real.transpose();
//    output_signal = mmse_matrix * received_signal;
//    for (int i = 0; i < 2 * Nt; ++i) {
//        y_hat[i] = output_signal(i, 0);
//    }
//}

void harddec_bits(int L, float* llr_bits, int* dec_bits) {
    for (int i = 0; i < L; i++)
    {
        if (llr_bits[i] < 0) dec_bits[i] = 1;
        else dec_bits[i] = 0;
    }
}



// 定义候选路径结构体
struct Candidate {
    double cost;               // 累计成本
    Eigen::VectorXcd symbols;  // 符号路径（长度为 P）
};

// K-Best 算法实现
std::vector<std::complex<double>> KBest(int k,
    const Eigen::MatrixXcd& R,
    const Eigen::VectorXcd& z_tilde,
    const std::vector<std::complex<double>>& sym) {
    int P = z_tilde.size();       // 传输天线数
    int SYM_SIZE = sym.size();     // 符号集大小

    // 初始化候选路径向量
    std::vector<Candidate> candidates;
    candidates.reserve(SYM_SIZE);

    // 初始化阶段：处理最后一层
#pragma omp parallel
    {
        // 为每个线程准备一个本地候选路径列表，避免竞争
        std::vector<Candidate> local_candidates;
        local_candidates.reserve(SYM_SIZE / omp_get_num_threads() + 1);

#pragma omp for nowait
        for (int i = 0; i < SYM_SIZE; ++i) {
            // 计算成本 |z_tilde[P-1] - R(P-1, P-1) * sym_i|^2
            std::complex<double> temp = z_tilde[P - 1] - R(P - 1, P - 1) * sym[i];
            double cost = std::norm(temp);

            // 初始化符号路径
            Eigen::VectorXcd symbols = Eigen::VectorXcd::Zero(P);
            symbols(P - 1) = sym[i];

            // 添加到本地列表
            local_candidates.emplace_back(Candidate{ cost, symbols });
        }

        // 将本地候选路径添加到全局列表
#pragma omp critical
        {
            candidates.insert(candidates.end(), local_candidates.begin(), local_candidates.end());
        }
    }

    // 如果符号集大小超过 k，保留前 k 个最低成本的候选路径
    if (SYM_SIZE > k) {
        std::nth_element(candidates.begin(), candidates.begin() + k, candidates.end(),
            [](const Candidate& a, const Candidate& b) { return a.cost < b.cost; });
        candidates.resize(k);
    }

    // 迭代处理其余的每一层
    for (int level = P - 2; level >= 0; --level) {
        // 准备存储所有扩展后的候选路径
        std::vector<Candidate> all_candidates;
        all_candidates.reserve(k * SYM_SIZE);

        // 并行化扩展候选路径
#pragma omp parallel
        {
            // 每个线程的本地候选路径列表
            std::vector<Candidate> thread_candidates;
            thread_candidates.reserve(SYM_SIZE);

#pragma omp for nowait
            for (int c = 0; c < candidates.size(); ++c) {
                const Candidate& current = candidates[c];

                // 计算上三角部分的乘积 R(j, j+1:end) * x
                Eigen::VectorXcd R_upper = R.block(level, level + 1, 1, P - level - 1);
                Eigen::VectorXcd x_extension = current.symbols.segment(level + 1, P - level - 1);
                std::complex<double> uppTrigProd = R_upper.dot(x_extension);

                double Rjj = R(level, level).real(); // 假设 R 是实数或只取实部

                // 遍历所有可能的符号
                for (int s = 0; s < SYM_SIZE; ++s) {
                    // 计算对角线部分 R(j,j) * sym_i
                    std::complex<double> diagProd = R(level, level) * sym[s];

                    // 计算残差
                    std::complex<double> residual = z_tilde(level) - uppTrigProd - diagProd;

                    // 计算新增成本 |residual|^2
                    double cost_inc = std::norm(residual);

                    // 计算新的累计成本
                    double new_cost = current.cost + cost_inc;

                    // 创建新的符号路径
                    Eigen::VectorXcd new_symbols = current.symbols;
                    new_symbols(level) = sym[s];

                    // 添加到线程本地列表
                    thread_candidates.emplace_back(Candidate{ new_cost, new_symbols });
                }
            }

            // 将本地线程的候选路径添加到全局列表
#pragma omp critical
            {
                all_candidates.insert(all_candidates.end(), thread_candidates.begin(), thread_candidates.end());
            }
        }

        // 从所有扩展后的候选路径中选择前 k 个最低成本的路径
        if (all_candidates.size() > static_cast<size_t>(k)) {
            std::nth_element(all_candidates.begin(), all_candidates.begin() + k, all_candidates.end(),
                [](const Candidate& a, const Candidate& b) { return a.cost < b.cost; });
            all_candidates.resize(k);
        }

        // 更新当前的候选路径为新的 k 个最佳路径
        candidates = std::move(all_candidates);
    }

    // 在最终的候选路径中选择累计成本最低的路径作为估计结果
    int best_idx = 0;
    double min_cost = candidates[0].cost;
    for (int i = 1; i < candidates.size(); ++i) {
        if (candidates[i].cost < min_cost) {
            min_cost = candidates[i].cost;
            best_idx = i;
        }
    }

    // 将最佳路径转换为 std::vector<std::complex<double>>
    std::vector<std::complex<double>> x_hat(P);
    for (int i = 0; i < P; ++i) {
        x_hat[i] = candidates[best_idx].symbols(i);
    }

    return x_hat;
}

/*
* KBest MIMO detection
* Matrix R: from QR decomposition of the real domain channel
* Vector sym_set: the value for each symbol in the symbol set
* Vector z: z = Q^Hy
* The detection module generates k detection paths,
* each path with a length of Nt
*/
void KBestDetection(MatrixXf R, MatrixXf z, MatrixXf H, float* y, float* ED, float* y_out, float** paths_extend, int k, int mod_type, int Nt, int Nr)
{
    int n_half_sym = 2 * Nt;
    int mod_type_half = mod_type / 2;
    int n_symbols = pow(2, mod_type_half); // n of different symbols in the symbols set
    int k_extend = k * n_symbols;
    // Pointers
    // R matrix in array
    float** R_arr = new float* [n_half_sym];
    for (int i = 0; i < n_half_sym; i++) R_arr[i] = new float[n_half_sym];
    // z in array
    float* z_arr = new float[n_half_sym];
    // Candidate Paths
    //float** paths_extend = new float* [k_extend]; // delete!!
    //for (int i = 0; i < k_extend; i++) { paths_extend[i] = new float[2 * Nt]; } // delete!!
    float** paths_extend_copy = new float* [k_extend]; // delete!!
    for (int i = 0; i < k_extend; i++) { paths_extend_copy[i] = new float[2 * Nt]; } // delete!!
    // PED (Partial Euclidian Distance)
    float* PED = new float[k_extend]; // delete!!
    float* PED_copy = new float[k_extend];
    float* off_diagonal_product = new float[k_extend];
    float* diagonal_product = new float[k_extend];
    int* idx_sorted = new int[k_extend];
    float* Hx_hat = new float[2 * Nr];
    // Pointers
    //cout << 5000 << endl;
    //Initialization
    for (int i = 0; i < n_half_sym; i++)
        for (int j = 0; j < n_half_sym; j++)
            R_arr[i][j] = R(i, j);
    for (int i = 0; i < n_half_sym; i++) { z_arr[i] = z(i, 0); }
    int n_paths = 1;
    int n_paths_pre = 1;
    for (int i = 0; i < k_extend; i++)
        for (int j = 0; j < n_half_sym; j++)
            paths_extend[i][j] = 0;
    for (int i = 0; i < k_extend; i++)
        for (int j = 0; j < n_half_sym; j++)
            paths_extend_copy[i][j] = 0;
    for (int i = 0; i < k_extend; i++) { PED[i] = 0; }
    for (int i = 0; i < k_extend; i++) { off_diagonal_product[i] = 0; }
    for (int i = 0; i < k_extend; i++) { diagonal_product[i] = 0; }


    for (int i = 0; i < 2 * Nt; i++)
    {
        n_paths_pre = n_paths;
        n_paths = n_paths * n_symbols;
        if (n_paths > k) { n_paths = k; }
        for (int i_path_pre = 0; i_path_pre < n_paths_pre; i_path_pre++) { PED_copy[i_path_pre] = PED[i_path_pre]; }

        // off diagonal products
        for (int i_path_pre = 0; i_path_pre < n_paths_pre; i_path_pre++)
        {
            off_diagonal_product[i_path_pre] = 0;
            int row = n_half_sym - 1 - i;
            for (int col = n_half_sym - 1; col > n_half_sym - 1 - i; col--)
            {
                off_diagonal_product[i_path_pre] += R_arr[row][col] * paths_extend[i_path_pre][col];
            }
        }

        // diagonal products
        for (int j = 0; j < n_symbols; j++)
        {
            switch (mod_type)
            {
            case 4:
                diagonal_product[j] = table_16qam_det[j] * R_arr[n_half_sym - 1 - i][n_half_sym - 1 - i];
                break;

            case 6:
                diagonal_product[j] = table_64qam_det[j] * R_arr[n_half_sym - 1 - i][n_half_sym - 1 - i];
                break;

            case 8:
                diagonal_product[j] = table_256qam_det[j] * R_arr[n_half_sym - 1 - i][n_half_sym - 1 - i];
                break;
            }

        }

        int i_sym = 2 * Nt - 1 - i;

        // Calculate |Q^Hy - Rx|(PED) for each extended path.
        for (int i_path_pre = 0; i_path_pre < n_paths_pre; i_path_pre++)
        {
            for (int j = 0; j < n_symbols; j++)
            {
                PED[i_path_pre * n_symbols + j] = PED_copy[i_path_pre];
                PED[i_path_pre * n_symbols + j] += abs(z_arr[i_sym] - off_diagonal_product[i_path_pre]
                    - diagonal_product[j]) * abs(z_arr[i_sym] - off_diagonal_product[i_path_pre]
                        - diagonal_product[j]);
            }
        }

        // Sort the paths
        // After execution of this function, PED is already sorted.
        my_sort_kbest(PED, idx_sorted, n_paths);

        //// Copy the best PEDs
        //for (int i_path = 0; i_path < n_paths; i_path++)
        //{
        //    PED_copy[i_path] = PED[idx_sorted[i_path]];
        //}

        //// Copy back PEDs
        //for (int i_path = 0; i_path < n_paths; i_path++)
        //{
        //    PED[i_path] = PED_copy[i_path];
        //}

        // Copy paths 
        for (int i_path_pre = 0; i_path_pre < n_paths_pre; i_path_pre++)
            for (int j = 0; j < n_half_sym; j++)
                paths_extend_copy[i_path_pre][j] = paths_extend[i_path_pre][j];

        // Copy the best paths
        for (int i_survive_path = 0; i_survive_path < n_paths; i_survive_path++)
        {
            int i_original_path_extend = idx_sorted[i_survive_path];
            int i_original_path = i_original_path_extend / n_symbols;
            int i_symbol_in_set = i_original_path_extend % n_symbols;
            for (int j = 0; j < n_half_sym; j++)
                paths_extend[i_survive_path][j] = paths_extend_copy[i_original_path][j];
            switch (mod_type)
            {
            case 4:
                paths_extend[i_survive_path][n_half_sym - i - 1] = table_16qam_det[i_symbol_in_set];
                break;
            case 6:
                paths_extend[i_survive_path][n_half_sym - i - 1] = table_64qam_det[i_symbol_in_set];
                break;
            case 8:
                paths_extend[i_survive_path][n_half_sym - i - 1] = table_256qam_det[i_symbol_in_set];
                break;
            }

        }
        //// test: Paths
        //cout << "Paths   " << "Step=   " << i << endl;
        //cout << "n_paths= " << n_paths << endl;
        //for (int i_path = 0; i_path < n_paths; i_path++)
        //{
        //    for (int j = 0; j < n_half_sym; j++)
        //    {
        //        cout << setw(10) << paths_extend[i_path][j] << "   ";
        //    }
        //    cout << endl;
        //    cout << endl;
        //    // test: Paths
        //}
        //cout << "PED" << endl;
        //for (int j = 0; j < k_extend; j++) { cout << setw(15) << PED[j]; }
        //cout << endl;
    }

    // Output
    for (int i = 0; i < n_half_sym; i++)
    {
        y_out[i] = paths_extend[0][i];
    }


    // Rearrange the paths with real Euclidian distance |y-Hx|
    for (int i_path = 0; i_path < k; i_path++)
    {
        ED[i_path] = PED[i_path];
        /*for (int i_row = 0; i_row < 2 * Nr; i_row++)
        {
            Hx_hat[i_row] = 0;
            for (int i_col = 0; i_col < 2 * Nt; i_col++)
            {
                Hx_hat[i_row] += H(i_row, i_col) * paths_extend[i_path][i_col];
            }
            ED[i_path] += abs(y[i_row] - Hx_hat[i_row]) * abs(y[i_row] - Hx_hat[i_row]);
        }*/
    }

    //cout << 5000 << endl;
    // Delete Pointers
    for (int i = 0; i < n_half_sym; i++)
        delete[]R_arr[i];
    delete[]R_arr;

    /* for (int i = 0; i < k_extend; i++)
     {
         delete[] paths_extend[i];
     }*/

    for (int i = 0; i < k_extend; i++)
    {
        delete[] paths_extend_copy[i];
    }

    delete[]paths_extend_copy;
    delete[]z_arr;
    // delete []paths_extend;
    delete[]PED;
    delete[]PED_copy;
    delete[]off_diagonal_product;
    delete[]diagonal_product;
    delete[]idx_sorted;
    delete[] Hx_hat;
    // cout << 5000 << endl;
}



/*
* KBest MIMO detection
* Matrix R: from QR decomposition of the real domain channel matrix
* Vector sym_set: the value for each symbol in the symbol set
* Vector z: z = Q^Hy
* The detection module generates k detection paths,
* each path with a length of Nt
* bits_decoded_partial: the result of past detection and spatial decoding
* s_decoded_bits_partial : the spatially decoded bits
* GGv : vertical encoding matrix
* A_Acv: records the frozen and non-frozen bits
* k : n of kbest detection paths
* decodingi : decodingi+1 is the number of columns that have already been spatially decoded
*/
void KBestDetection_JDD(MatrixXf R, MatrixXf z, MatrixXf H,
    float* y, float* ED, float* y_out, float** paths_extend, int** s_decoded_bits_partial, int** GGv,
    int* A_Acv, int* A_Ach, int k, int mod_type, int Nt, int Nr, int Nv, int Nh, int decodingi)
{
    bool flag_is_frozen_horizontal = 0;  // every frozen bit is a check point
    int n_half_sym = 2 * Nt;
    int mod_type_half = mod_type / 2;
    int n_symbols = pow(2, mod_type_half); // n of different symbols in the symbols set
    int k_extend = k * n_symbols;
    // Pointers
    // R matrix in array
    int* flag_is_valid_path = new int[k_extend];
    int* flag_is_valid_path_small = new int[n_symbols];
    int* xor_res_partial = new int[Nv];
    float** R_arr = new float* [n_half_sym];
    for (int i = 0; i < n_half_sym; i++) R_arr[i] = new float[n_half_sym];
    // z in array
    float* z_arr = new float[n_half_sym];
    // Candidate Paths
    //float** paths_extend = new float* [k_extend]; // delete!!
    //for (int i = 0; i < k_extend; i++) { paths_extend[i] = new float[2 * Nt]; } // delete!!
    float** paths_extend_copy = new float* [k_extend]; // delete!!
    for (int i = 0; i < k_extend; i++) { paths_extend_copy[i] = new float[2 * Nt]; } // delete!!
    // PED (Partial Euclidian Distance)
    float* PED = new float[k_extend]; // delete!!
    float* PED_copy = new float[k_extend];
    float* off_diagonal_product = new float[k_extend];
    float* diagonal_product = new float[k_extend];
    int* idx_sorted = new int[k_extend];
    float* Hx_hat = new float[2 * Nr];
    // Pointers

    //Initialization
    for (int i = 0; i < n_half_sym; i++)
        for (int j = 0; j < n_half_sym; j++)
            R_arr[i][j] = R(i, j);
    for (int i = 0; i < n_half_sym; i++) { z_arr[i] = z(i, 0); }
    int n_paths = 1;
    int n_paths_pre = 1;
    for (int i = 0; i < k_extend; i++)
        for (int j = 0; j < n_half_sym; j++)
            paths_extend[i][j] = 0;
    for (int i = 0; i < k_extend; i++)
        for (int j = 0; j < n_half_sym; j++)
            paths_extend_copy[i][j] = 0;
    for (int i = 0; i < k_extend; i++) { PED[i] = 0; }
    for (int i = 0; i < k_extend; i++) { off_diagonal_product[i] = 0; }
    for (int i = 0; i < k_extend; i++) { diagonal_product[i] = 0; }

    // whether the algorithm should check the horizontally frozen bits
    flag_is_frozen_horizontal = !(A_Ach[Nh - 1 - decodingi]);
    // calculate xor_res_partial (horizontal reverse-encoding)
    for (int i_code = 0; i_code < Nv; i_code++)
    {
        xor_res_partial[i_code] = 0;
        int col_G = Nh - 1 - decodingi;
        int row_G_min = Nh - decodingi;
        for (int row_G = row_G_min; row_G < Nh; row_G++)
        {
            xor_res_partial[i_code] ^= (s_decoded_bits_partial[i_code][row_G] & GGv[row_G][col_G]);
        }
    }


    // G[Nh - decodingi - 1][col_G] & curr_decoded
    for (int i = 0; i < 2 * Nt; i++)
    {
        n_paths_pre = n_paths;
        n_paths = n_paths * n_symbols;
        if (n_paths > k) { n_paths = k; }
        for (int i_path_pre = 0; i_path_pre < n_paths_pre; i_path_pre++) { PED_copy[i_path_pre] = PED[i_path_pre]; }

        // off diagonal products
        for (int i_path_pre = 0; i_path_pre < n_paths_pre; i_path_pre++)
        {
            off_diagonal_product[i_path_pre] = 0;
            int row = n_half_sym - 1 - i;
            for (int col = n_half_sym - 1; col > n_half_sym - 1 - i; col--)
            {
                off_diagonal_product[i_path_pre] += R_arr[row][col] * paths_extend[i_path_pre][col];
            }
        }

        // diagonal products
        for (int j = 0; j < n_symbols; j++)
        {
            switch (mod_type)
            {
            case 4:
                diagonal_product[j] = table_16qam_det[j] * R_arr[n_half_sym - 1 - i][n_half_sym - 1 - i];
                break;

            case 6:
                diagonal_product[j] = table_64qam_det[j] * R_arr[n_half_sym - 1 - i][n_half_sym - 1 - i];
                break;

            case 8:
                diagonal_product[j] = table_256qam_det[j] * R_arr[n_half_sym - 1 - i][n_half_sym - 1 - i];
                break;
            }

        }

        int i_sym = 2 * Nt - 1 - i;

        // Calculate |Q^Hy - Rx|(PED) for each extended path.
        for (int i_path_pre = 0; i_path_pre < n_paths_pre; i_path_pre++)
        {
            for (int j = 0; j < n_symbols; j++)
            {
                PED[i_path_pre * n_symbols + j] = PED_copy[i_path_pre];
                PED[i_path_pre * n_symbols + j] += abs(z_arr[i_sym] - off_diagonal_product[i_path_pre]
                    - diagonal_product[j]) * abs(z_arr[i_sym] - off_diagonal_product[i_path_pre]
                        - diagonal_product[j]);
            }
        }


        // If the detector has arrived at a checkpoint
        if (flag_is_frozen_horizontal)
        {
            for (int j = 0; j < n_symbols; j++)
            {

            }
        }
        // Sort the paths
        // After execution of this function, PED is already sorted.
        my_sort_kbest(PED, idx_sorted, n_paths);

        //// Copy the best PEDs
        //for (int i_path = 0; i_path < n_paths; i_path++)
        //{
        //    PED_copy[i_path] = PED[idx_sorted[i_path]];
        //}

        //// Copy back PEDs
        //for (int i_path = 0; i_path < n_paths; i_path++)
        //{
        //    PED[i_path] = PED_copy[i_path];
        //}

        // Copy paths 
        for (int i_path_pre = 0; i_path_pre < n_paths_pre; i_path_pre++)
            for (int j = 0; j < n_half_sym; j++)
                paths_extend_copy[i_path_pre][j] = paths_extend[i_path_pre][j];

        // Copy the best paths
        for (int i_survive_path = 0; i_survive_path < n_paths; i_survive_path++)
        {
            int i_original_path_extend = idx_sorted[i_survive_path];
            int i_original_path = i_original_path_extend / n_symbols;
            int i_symbol_in_set = i_original_path_extend % n_symbols;
            for (int j = 0; j < n_half_sym; j++)
                paths_extend[i_survive_path][j] = paths_extend_copy[i_original_path][j];
            switch (mod_type)
            {
            case 4:
                paths_extend[i_survive_path][n_half_sym - i - 1] = table_16qam_det[i_symbol_in_set];
                break;
            case 6:
                paths_extend[i_survive_path][n_half_sym - i - 1] = table_64qam_det[i_symbol_in_set];
                break;
            case 8:
                paths_extend[i_survive_path][n_half_sym - i - 1] = table_256qam_det[i_symbol_in_set];
                break;
            }

        }
        //// test: Paths
        //cout << "Paths   " << "Step=   " << i << endl;
        //cout << "n_paths= " << n_paths << endl;
        //for (int i_path = 0; i_path < n_paths; i_path++)
        //{
        //    for (int j = 0; j < n_half_sym; j++)
        //    {
        //        cout << setw(10) << paths_extend[i_path][j] << "   ";
        //    }
        //    cout << endl;
        //    cout << endl;
        //    // test: Paths
        //}
        //cout << "PED" << endl;
        //for (int j = 0; j < k_extend; j++) { cout << setw(15) << PED[j]; }
        //cout << endl;
    }

    // Output
    for (int i = 0; i < n_half_sym; i++)
    {
        y_out[i] = paths_extend[0][i];
    }


    // Rearrange the paths with real Euclidian distance |y-Hx|
    for (int i_path = 0; i_path < k; i_path++)
    {
        ED[i_path] = 0;
        for (int i_row = 0; i_row < 2 * Nr; i_row++)
        {
            Hx_hat[i_row] = 0;
            for (int i_col = 0; i_col < 2 * Nt; i_col++)
            {
                Hx_hat[i_row] += H(i_row, i_col) * paths_extend[i_path][i_col];
            }
            ED[i_path] += abs(y[i_row] - Hx_hat[i_row]) * abs(y[i_row] - Hx_hat[i_row]);
        }
    }


    // Delete Pointers
    for (int i = 0; i < n_half_sym; i++)
        delete[]R_arr[i];
    delete[]R_arr;

    /* for (int i = 0; i < k_extend; i++)
     {
         delete[] paths_extend[i];
     }*/

    for (int i = 0; i < k_extend; i++)
    {
        delete[] paths_extend_copy[i];
    }

    delete[]paths_extend_copy;
    delete[]z_arr;
    // delete []paths_extend;
    delete[]PED;
    delete[]PED_copy;
    delete[]off_diagonal_product;
    delete[]diagonal_product;
    delete[]idx_sorted;
    delete[] Hx_hat;
    delete[] flag_is_valid_path;
    delete[] flag_is_valid_path_small;
    delete[] xor_res_partial;
}



void my_sort_kbest(float* arr, int* idx_sort, int N) {
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


// decide whether a symbol is valid when attached to a bit string
// xor_res_pre: the xor results of decoded bits: x(k) ^ x(k+1) ^ ... ^ x(N)
// idx_in_sym_set: index of symbol in the set
bool isValidSymbol(int idx_sym_in_set, int xor_res_pre, int mod_type, int decode_len, vector<int> frozen_bits_pos)
{
    int modtype_half = mod_type / 2;
    int* encoded_res = new int[modtype_half];
    for (int i_bit = 0; i_bit < modtype_half; i_bit++)
    {

    }
    delete[] encoded_res;
}









/*
* KBest MIMO detection
* Matrix R: from QR decomposition of the real domain channel
* Vector sym_set: the value for each symbol in the symbol set
* Vector z: z = Q^Hy
* LLR: ln(p(sum=0)/p(sum=1))
* The detection module generates k detection paths,
* each path with a length of Nt
*/
// KBEST Detection with soft JDD
void KBestDetection_soft_JDD(MatrixXf R, MatrixXf z, MatrixXf H,
    float* y, float* ED, float* y_out, float* LLR_in, float** paths_extend,
    int k, int mod_type, int Nt, int Nr, int decoding_stage)
{
    int n_half_sym = 2 * Nt;
    int mod_type_half = mod_type / 2;
    int n_symbols = pow(2, mod_type_half); // n of different symbols in the symbols set
    int k_extend = k * n_symbols;
    // Pointers
    // R matrix in array
    float** R_arr = new float* [n_half_sym];
    for (int i = 0; i < n_half_sym; i++) R_arr[i] = new float[n_half_sym];
    // z in array
    float* z_arr = new float[n_half_sym];
    // Candidate Paths
    //float** paths_extend = new float* [k_extend]; // delete!!
    //for (int i = 0; i < k_extend; i++) { paths_extend[i] = new float[2 * Nt]; } // delete!!
    float** paths_extend_copy = new float* [k_extend]; // delete!!
    for (int i = 0; i < k_extend; i++) { paths_extend_copy[i] = new float[2 * Nt]; } // delete!!
    // PED (Partial Euclidian Distance)
    float* PED = new float[k_extend]; // delete!!
    float* PED_copy = new float[k_extend];
    float* off_diagonal_product = new float[k_extend];
    float* diagonal_product = new float[k_extend];
    int* idx_sorted = new int[k_extend];
    float* Hx_hat = new float[2 * Nr];
    float* PM_dec = new float[n_symbols];
    // Pointers

    //Initialization
    for (int i = 0; i < n_half_sym; i++)
        for (int j = 0; j < n_half_sym; j++)
            R_arr[i][j] = R(i, j);
    for (int i = 0; i < n_half_sym; i++) { z_arr[i] = z(i, 0); }
    int n_paths = 1;
    int n_paths_pre = 1;
    for (int i = 0; i < k_extend; i++)
        for (int j = 0; j < n_half_sym; j++)
            paths_extend[i][j] = 0;
    for (int i = 0; i < k_extend; i++)
        for (int j = 0; j < n_half_sym; j++)
            paths_extend_copy[i][j] = 0;
    for (int i = 0; i < k_extend; i++) { PED[i] = 0; }
    for (int i = 0; i < k_extend; i++) { off_diagonal_product[i] = 0; }
    for (int i = 0; i < k_extend; i++) { diagonal_product[i] = 0; }
    for (int i = 0; i < n_symbols; i++) { PM_dec[i] = 0; }


    for (int i = 0; i < 2 * Nt; i++)
    {
        n_paths_pre = n_paths;
        n_paths = n_paths * n_symbols;
        if (n_paths > k) { n_paths = k; }
        for (int i_path_pre = 0; i_path_pre < n_paths_pre; i_path_pre++) { PED_copy[i_path_pre] = PED[i_path_pre]; }

        // off diagonal products
        for (int i_path_pre = 0; i_path_pre < n_paths_pre; i_path_pre++)
        {
            off_diagonal_product[i_path_pre] = 0;
            int row = n_half_sym - 1 - i;
            for (int col = n_half_sym - 1; col > n_half_sym - 1 - i; col--)
            {
                off_diagonal_product[i_path_pre] += R_arr[row][col] * paths_extend[i_path_pre][col];
            }
        }

        // diagonal products
        for (int j = 0; j < n_symbols; j++)
        {
            switch (mod_type)
            {
            case 4:
                diagonal_product[j] = table_16qam_det[j] * R_arr[n_half_sym - 1 - i][n_half_sym - 1 - i];
                break;

            case 6:
                diagonal_product[j] = table_64qam_det[j] * R_arr[n_half_sym - 1 - i][n_half_sym - 1 - i];
                break;

            case 8:
                diagonal_product[j] = table_256qam_det[j] * R_arr[n_half_sym - 1 - i][n_half_sym - 1 - i];
                break;
            }

        }
        int i_sym = 2 * Nt - 1 - i;
        // PM from decoder

        for (int j = 0; j < n_symbols; j++)
        {
            PM_dec[j] = 0;
            float PM_dec_onebit = 0;
            for (int i_bit_in_sym = 0; i_bit_in_sym < mod_type_half; i_bit_in_sym++)
            {
                PM_dec_onebit = LLR_in[i_sym * mod_type_half + i_bit_in_sym] * (1 - 2 * s_bits_r4_det[j][i_bit_in_sym]);
                PM_dec[j] -= PM_dec_onebit; //-=
            }
        }



        // Calculate |Q^Hy - Rx|(PED) for each extended path.
        for (int i_path_pre = 0; i_path_pre < n_paths_pre; i_path_pre++)
        {
            for (int j = 0; j < n_symbols; j++)
            {
                PED[i_path_pre * n_symbols + j] = PED_copy[i_path_pre];
                PED[i_path_pre * n_symbols + j] += abs(z_arr[i_sym] - off_diagonal_product[i_path_pre]
                    - diagonal_product[j]) * abs(z_arr[i_sym] - off_diagonal_product[i_path_pre]
                        - diagonal_product[j]);
                // PED[i_path_pre * n_symbols + j] += (beta_JDD * PM_dec[j]);//0.35 VERY GOOD MUCH BETTER THAN 0.4!  // 0.8, 13dB, L=32 may be improvements also 0.6
                 // float decod_stage_f = float(decoding_stage) / float(Nh);
                // float decod_stage_f = 1-sin(float(decoding_stage) / float(Nh)*3.1415926/2);
                PED[i_path_pre * n_symbols + j] += (beta_JDD * PM_dec[j]);
                // cout << "Original PED"<< abs(z_arr[i_sym] - off_diagonal_product[i_path_pre]
                    /*- diagonal_product[j]) * abs(z_arr[i_sym] - off_diagonal_product[i_path_pre]
                        - diagonal_product[j]) << endl;*/
                        // cout << "PM_dec" << PM_dec[j] << endl;
            }
        }

        // Sort the paths
        // After execution of this function, PED is already sorted.
        my_sort_kbest(PED, idx_sorted, n_paths);

        //// Copy the best PEDs
        //for (int i_path = 0; i_path < n_paths; i_path++)
        //{
        //    PED_copy[i_path] = PED[idx_sorted[i_path]];
        //}

        //// Copy back PEDs
        //for (int i_path = 0; i_path < n_paths; i_path++)
        //{
        //    PED[i_path] = PED_copy[i_path];
        //}

        // Copy paths 
        for (int i_path_pre = 0; i_path_pre < n_paths_pre; i_path_pre++)
            for (int j = 0; j < n_half_sym; j++)
                paths_extend_copy[i_path_pre][j] = paths_extend[i_path_pre][j];

        // Copy the best paths
        for (int i_survive_path = 0; i_survive_path < n_paths; i_survive_path++)
        {
            int i_original_path_extend = idx_sorted[i_survive_path];
            int i_original_path = i_original_path_extend / n_symbols;
            int i_symbol_in_set = i_original_path_extend % n_symbols;
            for (int j = 0; j < n_half_sym; j++)
                paths_extend[i_survive_path][j] = paths_extend_copy[i_original_path][j];
            switch (mod_type)
            {
            case 4:
                paths_extend[i_survive_path][n_half_sym - i - 1] = table_16qam_det[i_symbol_in_set];
                break;
            case 6:
                paths_extend[i_survive_path][n_half_sym - i - 1] = table_64qam_det[i_symbol_in_set];
                break;
            case 8:
                paths_extend[i_survive_path][n_half_sym - i - 1] = table_256qam_det[i_symbol_in_set];
                break;
            }

        }
        //// test: Paths
        //cout << "Paths   " << "Step=   " << i << endl;
        //cout << "n_paths= " << n_paths << endl;
        //for (int i_path = 0; i_path < n_paths; i_path++)
        //{
        //    for (int j = 0; j < n_half_sym; j++)
        //    {
        //        cout << setw(10) << paths_extend[i_path][j] << "   ";
        //    }
        //    cout << endl;
        //    cout << endl;
        //    // test: Paths
        //}
        //cout << "PED" << endl;
        //for (int j = 0; j < k_extend; j++) { cout << setw(15) << PED[j]; }
        //cout << endl;
    }

    // Output
    for (int i = 0; i < n_half_sym; i++)
    {
        y_out[i] = paths_extend[0][i];
    }


    // Rearrange the paths with real Euclidian distance |y-Hx|
    for (int i_path = 0; i_path < k; i_path++)
    {
        ED[i_path] = PED[i_path];
        /*for (int i_row = 0; i_row < 2 * Nr; i_row++)
        {
            Hx_hat[i_row] = 0;
            for (int i_col = 0; i_col < 2 * Nt; i_col++)
            {
                Hx_hat[i_row] += H(i_row, i_col) * paths_extend[i_path][i_col];
            }
            ED[i_path] += abs(y[i_row] - Hx_hat[i_row]) * abs(y[i_row] - Hx_hat[i_row]);
        }*/
    }


    // Delete Pointers
    for (int i = 0; i < n_half_sym; i++)
        delete[]R_arr[i];
    delete[]R_arr;

    /* for (int i = 0; i < k_extend; i++)
     {
         delete[] paths_extend[i];
     }*/

    for (int i = 0; i < k_extend; i++)
    {
        delete[] paths_extend_copy[i];
    }

    delete[]paths_extend_copy;
    delete[]z_arr;
    // delete []paths_extend;
    delete[]PED;
    delete[]PED_copy;
    delete[]off_diagonal_product;
    delete[]diagonal_product;
    delete[]idx_sorted;
    delete[] Hx_hat;
    delete[] PM_dec;
}



// V is Sigma 
// E is mu
void EPD_detection_float(int flag, size_t Nt, size_t Nr, size_t iter, size_t mod_type,
    float es, float N0, float* HTH, float* HTY, int temp,
    float* extr_mean, float* extr_var,
    double* vec_alpha_in, double* vec_gamma_in, double* p_symbol_rearrange)
{
    int Nt2 = Nt * 2;
    int Nr2 = Nr * 2;
    int order = pow(2, mod_type / 2);
    int ttt = order;
    int i, j, cnt;
    float beta = 0, efcelong = 1e-8, inita = 0.001;
    //float p_symbol_rearrange = pow(2, (-MODE_TYPE / 2));  p_symbol_rearrange is specified by information from the decoder
    //	printf("N0:%f\n", N0);
    // Initialization
    float*lambda = new float[2 * Nr];                           memset(lambda, 0, sizeof(float) * 2 * Nr);
    float*gamma = new float[2 * Nr];                            memset(gamma, 0, sizeof(float) * 2 * Nr);
    float*ma = new float[4 * Nr * Nr];                          memset(ma, 0, sizeof(float) * 4 * Nr * Nr);
    float*mb = new float[2 * Nr];                               memset(mb, 0, sizeof(float) * 2 * Nr);
    float*eye_ma = new float[4 * Nr * Nr];                      memset(eye_ma, 0, sizeof(float) * 4 * Nr * Nr);
    float*lambda_re = new float[2 * Nr];                        memset(lambda_re, 0, sizeof(float) * 2 * Nr);
    float* gamma_re = new float[2 * Nr];                         memset(gamma_re, 0, sizeof(float) * 2 * Nr);
    float* lambda_temp = new float[4 * Nr * Nr];                 memset(lambda_temp, 0, sizeof(float) * 4 * Nr * Nr);
    float* gamma_temp = new float[2 * Nr];                       memset(gamma_temp, 0, sizeof(float) * 2 * Nr);
    float* Bold_V = new float[4 * Nr * Nr];                      memset(Bold_V, 0, sizeof(float) * 4 * Nr * Nr);
    float* Bold_E = new float[2 * Nr];                           memset(Bold_E, 0, sizeof(float) * 2 * Nr);
    float* Vq = new float[2 * Nr];                               memset(Vq, 0, sizeof(float) * 2 * Nr);
    float* Eq = new float[2 * Nr];                               memset(Eq, 0, sizeof(float) * 2 * Nr);
    float* h = new float[2 * Nr];                                memset(h, 0, sizeof(float) * 2 * Nr);
    float* t = new float[2 * Nr];                                memset(t, 0, sizeof(float) * 2 * Nr);
    float* value = new float[2 * Nr * ttt];                      memset(value, 0, sizeof(float) * 2 * Nr * ttt);
    float* dex = new float[2 * Nr];                              memset(dex, 0, sizeof(float) * 2 * Nr);
    float* prob_muti = new float[2 * Nr * ttt];                  memset(prob_muti, 0, sizeof(float) * 2 * Nr * ttt);
    float* prob_out = new float[2 * Nr * ttt];                   memset(prob_out, 0, sizeof(float) * 2 * Nr * ttt);
    float* miu = new float[2 * Nr];                              memset(miu, 0, sizeof(float) * 2 * Nr);
    float* delta = new float[2 * Nr];                            memset(delta, 0, sizeof(float) * 2 * Nr);
    float* result = new float[2 * Nr];                           memset(result, 0, sizeof(float) * 2 * Nr);
    //cf_t* symb = new cf_t[Nr];                                  memset(symb, 0, sizeof(cf_t) * Nr);
    //cf_t* constellations = new cf_t[ttt * ttt];                 memset(constellations, 0, sizeof(cf_t) * ttt * ttt);
    float* consr = new float[ttt];                               memset(consr, 0, sizeof(float) * ttt);

    float* p_symbol = new float[2 * Nr * ttt];                   memset(p_symbol, 0, sizeof(float) * 2 * Nr * ttt);

    float* v = new float[2 * Nr];
    float* rou = new float[2 * Nr];
    float* miu_temp = new float[2 * Nr];
    float* t_new = new float[2 * Nr];

    float* var = new float[1 * 2 * Nr];
    float* mean = new float[1 * 2 * Nr];

    /*float* extr_mean = new float[1 * 2 * Nr];
    float* extr_var = new float[1 * 2 * Nr];*/
    float* extr_alpha = new float[1 * 2 * Nr];
    float* extr_gamma = new float[1 * 2 * Nr];

    // Initialization
    for (i = 0; i < Nt2; i++)
    {
        for (j = 0; j < order; j++)
        {
            p_symbol[i * order + j] = p_symbol_rearrange[i * order + j];
        }
    }
    for (i = 0; i < Nt2; i++)
    {
        lambda[i] = vec_alpha_in[i];
        gamma[i] = vec_gamma_in[i];
        gamma_temp[i] = vec_gamma_in[i];
    }
    for (i = 0; i < Nt2; i++)
    {
        //lambda[i] = inita;
        for (j = 0; j < Nt2; j++)
        {
            if (i == j)
            {
                eye_ma[i * 2 * Nt + j] = 1;
                lambda_temp[i * Nt2 + j] = lambda[i]; //modified 0408
            }

            else {
                lambda_temp[i * Nt2 + j] = 0;
                eye_ma[i * 2 * Nt + j] = 0;
            }
        }
    }
    /*******************?  ?  ?  *************************/
    for (i = 0; i < Nt2 * Nt2; i++)		ma[i] = 0;
    for (i = 0; i < Nt2; i++)			mb[i] = 0;
    for (i = 0; i < Nt2; i++)			gamma_re[i] = 0;
    //for (i = 0; i < Nt2; i++)			gamma_temp[i] = 0;
    //for (i = 0; i < Nt2; i++)			gamma[i] = 0;
    for (i = 0; i < Nt2; i++)			lambda_re[i] = 0;

    for (i = 0; i < Nt2; i++)			h[i] = 0;
    for (i = 0; i < Nt2; i++)			t[i] = 0;
    for (i = 0; i < Nt2; i++)			Vq[i] = 0;
    for (i = 0; i < Nt2; i++)			Eq[i] = 0;

    for (i = 0; i < Nt2; i++)			miu[i] = 0;
    for (i = 0; i < Nt2 * Nt2; i++)		Bold_V[i] = 0;
    for (i = 0; i < Nt2; i++)			Bold_E[i] = 0;

    for (i = 0; i < Nt2; i++)			delta[i] = 0;

    for (i = 0; i < Nt2 * order; i++)	value[i] = 0;
    for (i = 0; i < Nt2 * order; i++)	prob_muti[i] = 0;
    for (i = 0; i < Nt2 * order; i++)	prob_out[i] = 0;

    /**********************************************************/

    cblas_saxpy(Nt2 * Nt2, (1.0 / N0), HTH, 1, ma, 1);   // ma = ma+HTH/N0 => ma = HTH/N0
    cblas_saxpy(Nt2, (1.0 / N0), HTY, 1, mb, 1);   // mb = mb+HTY/N0 => mb = HTY/n0
    cblas_saxpy(Nt2 * Nt2, 1.0, ma, 1, lambda_temp, 1); // lambda_temp = lambda_temp + ma => lambda_temp = diag (Alpha) + HTH/N0
    inverseMatrix(lambda_temp, Nt2); // lambda_temp = lambda_temp^{-1} = (diag (Alpha) + HTH/N0)^-1

    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Nt2, Nt2, Nt2, 1.0, lambda_temp, Nt2, eye_ma, Nt2, 0.0, Bold_V, Nt2);
    // Bold_V = (diag (Alpha) + HTH/N0)^-1
    cblas_scopy(Nt2, gamma, 1, gamma_temp, 1); // gamma_temp = gamma
    cblas_saxpy(Nt2, 1.0, mb, 1, gamma_temp, 1); // gamma_temp = mb+gamma = HTY/n0 + gamma

    // mu = Sigma(H^Ty + Gamma)
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Nt2, 1, Nt2, 1.0, Bold_V, Nt2, gamma_temp, 1, 0.0, Bold_E, 1); // Bold_E = Bold_V*(HTY/n0 + gamma)

    // test: Bold_E and Bold_V
    /*cout << endl;
    cout << "Bold_V, T=0" << endl;
    for (int i = 0; i < Nr2; i++)
    {
        cout << Bold_V[i]<<"    ";
    }
    cout << endl;
    cout << "Bold_E, T=" << 0 << endl;
    for (int i = 0; i < Nt2; i++)
    {
        cout << Bold_E[i] << "    ";
    }
    cout << endl;*/


    for (cnt = 0; cnt < iter; cnt++)
    {
       // beta = fmin((exp(cnt / 1.5) / 10.0), 0.5);
        beta = 0.7;
        for (i = 0; i < Nt2; i++) { gamma_temp[i] = 0; }
        for (i = 0; i < Nt2 * Nt2; i++) { lambda_temp[i] = 0; }

        cblas_scopy(Nt2, lambda, 1, lambda_re, 1);
        cblas_scopy(Nt2, gamma, 1, gamma_re, 1);

        for (i = 0; i < Nt2; i++)
        {
            Vq[i] = Bold_V[i * Nt2 + i];
        }
        cblas_scopy(Nt2, Bold_E, 1, Eq, 1);

        //computing the deviation and average of the cavity marginal

        for (i = 0; i < Nt2; i++) h[i] = Vq[i] / (1.0 - Vq[i] * lambda[i]);
        for (i = 0; i < Nt2; i++) t[i] = h[i] * (Eq[i] / Vq[i] - gamma[i]);

        // ONLY SUPPORTS 16QAM!!
        for (i = 0; i < Nt2; i++)
        {
            double muti_sum = 0;
            double out_sum = 0;
            for (j = 0; j < order; j++)
            {
                value[i * order + j] = (-0.5 * pow((table_16qam_det[j] - t[i]), 2)) / (1.0 * h[i]);
                prob_muti[i * order + j] = fmax(exp(value[i * order + j]) / sqrtf(2 * PI * h[i]), 1e-15);
                muti_sum += prob_muti[i * order + j];
            }
            for (j = 0; j < order; j++)
            {
                prob_muti[i * order + j] = prob_muti[i * order + j] / muti_sum;

                prob_out[i * order + j] = fmax(prob_muti[i * order + j] * p_symbol_rearrange[i * order + j], 1e-15);
                out_sum += prob_out[i * order + j];
            }
            for (j = 0; j < order; j++)
            {
                prob_out[i * order + j] = prob_out[i * order + j] / out_sum;
            }
        }

        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, Nt2, 1, order, 1.0, prob_out, order, table_16qam_det, order, 0.0, miu, 1);

        /* delta */
        for (i = 0; i < Nt2; i++)
        {
            double delta_temp = 0;
            for (j = 0; j < order; j++)
            {
                delta_temp += (prob_out[i * order + j] * pow(abs(table_16qam_det[j] - miu[i]), 2));
            }
            delta[i] = fmax(delta_temp, 1e-8);
        }

        for (i = 0; i < Nt2; i++)
        {
            gamma[i] = (1 - beta) * gamma_re[i] + beta * (miu[i] / delta[i] - Eq[i] / Vq[i]);
            lambda[i] = (1 - beta) * lambda_re[i] + beta * (1.0 / delta[i] - 1.0 / Vq[i]);
        }

        // dealing with negative deviation
        for (i = 0; i < Nt2; i++)
        {
            if (lambda[i] <= 0)
            {
                gamma[i] = gamma_re[i];
                lambda[i] = lambda_re[i];
            }

        }
        for (i = 0; i < Nt2; i++)
        {
            for (j = 0; j < Nt2; j++)
            {
                if (i == j)
                {
                    lambda_temp[i * Nt2 + j] = lambda[i];
                }
                else
                {
                    lambda_temp[i * Nt2 + j] = 0;
                }
            }
        }

        cblas_saxpy(Nt2 * Nt2, 1.0, ma, 1, lambda_temp, 1);

        inverseMatrix(lambda_temp, Nt2);

        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Nt2, Nt2, Nt2, 1.0, lambda_temp, Nt2, eye_ma, Nt2, 0.0, Bold_V, Nt2);


        cblas_scopy(Nt2, gamma, 1, gamma_temp, 1);

        cblas_saxpy(Nt2, 1.0, mb, 1, gamma_temp, 1);

        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Nt2, 1, Nt2, 1.0, Bold_V, Nt2, gamma_temp, 1, 0.0, Bold_E, 1);

        /*cout << "TEST:EP DETECTION Bold_E" << endl;
        for (int i = 0; i < Nr2; i++) std::cout << Bold_E[i] << endl;
        cout << "TEST:EP DETECTION Bold_E" << endl;*/
        /*cout << endl;
        cout << "Bold_V, T="<<cnt << endl;
        for (int i = 0; i < Nr2; i++)
        {
            cout << Bold_V[i] << "    ";
        }
        cout << endl;

        cout << endl;
        cout << "Bold_E, T=" << cnt << endl;
        for (int i = 0; i < Nt2; i++)
        {
            cout << Bold_E[i] << "    ";
        }
        cout << endl;*/

    }
    cblas_scopy(Nt2, lambda, 1, lambda_re, 1);
    cblas_scopy(Nt2, gamma, 1, gamma_re, 1);

    for (i = 0; i < Nt2; i++)
    {
        Vq[i] = Bold_V[i * Nt2 + i];
    }

    cblas_scopy(Nt2, Bold_E, 1, Eq, 1);
    for (i = 0; i < Nt2; i++) h[i] = Vq[i] / (1 - Vq[i] * lambda[i]);
    for (i = 0; i < Nt2; i++) t[i] = h[i] * (Eq[i] / Vq[i] - gamma[i]);

    for (i = 0; i < Nt; i++)
    {
        extr_mean[temp * Nt + i] = t[i];
        extr_mean[temp * Nt + i + Nt * 1] = t[i + Nt];
        extr_var[temp * Nt + i] = h[i];
        extr_var[temp * Nt + i + Nt * 1] = h[i + Nt];
       /* extr_mean[temp * Nt + i] = Eq[i];
        extr_mean[temp * Nt + i + Nt * 1] = Eq[i + Nt];*/
    }
    // the function must reserve the extrinsic(cavity) obtained in the last turbo iteration
    // 删除指针并释放内存
    delete[] lambda;
    delete[] gamma;
    delete[] ma;
    delete[] mb;
    delete[] eye_ma;
    delete[] lambda_re;
    delete[] gamma_re;
    delete[] lambda_temp;
    delete[] gamma_temp;
    delete[] Bold_V;
    delete[] Bold_E;
    delete[] Vq;
    delete[] Eq;
    delete[] h;
    delete[] t;
    delete[] value;
    delete[] dex;
    delete[] prob_muti;
    delete[] prob_out;
    delete[] miu;
    delete[] delta;
    delete[] result;
    // 删除注释中指针的内存，如果有用到，取消注释
    // delete[] symb;
    // delete[] constellations;
    delete[] consr;
    delete[] p_symbol;
    delete[] v;
    delete[] rou;
    delete[] miu_temp;
    delete[] t_new;

    delete[] var;
    delete[] mean;
    // 如果用到这些指针，也释放内存
    // delete[] extr_mean;
    // delete[] extr_var;
    delete[] extr_alpha;
    delete[] extr_gamma;

}





