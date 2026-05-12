/*
form a large list based on the decoding result of the component codeword
@param in: sub_codeword_list: the candidate list for each sub-codeword; size is Nv*[L*Nh]
@param in: pm: the list of path metric of each sub-codeword; size is Nv*L
@param in: L: the size of each subcodeword_list
@param in: Nv: the number of sub-codewords
@param in: Nh: the length of each sub-codeword
@param in: L_large: the size of the final large list
@param out: whole_codeword_list: the final large list; size is L_large*[Nv*Nh]
*/

#include <stdlib.h>
#include <string.h>
#include <vector>
#include "modulation.h"
#include "MIMO_Function.h"
#include "mkl.h";
#include "crc_list.h"
#include "Polar_Encoder.h"


float table_qpsk_crc[2] = { -0.707107f, 0.707107f };
float table_16qam_crc[4] = { -0.316228f, -0.948683f,  0.316228f,  0.948683f };
//float table_16qam[4]={}
float table_64qam_crc[8] = { -0.462910f, -0.154303f, -0.771517f, -1.08012f,	0.462910f,  0.154303f,  0.771517f,  1.08012f };
float table_256qam_crc[16] =
{ -0.383482f, -0.536875f, -0.230089f, -0.07669f,
-0.843661f, -0.690268f, -0.997054f, -1.15044f,
0.383482f, 0.536875f, 0.230089f, 0.07669f,
0.843661f, 0.690268f, 0.997054f, 1.15044f };

typedef struct {
    int index;
    float pm;
} IndexPmPair;

static int compare_pairs(const void* a, const void* b) {
    const IndexPmPair* pa = (const IndexPmPair*)a;
    const IndexPmPair* pb = (const IndexPmPair*)b;
    if (pa->pm > pb->pm) return -1;
    if (pa->pm < pb->pm) return 1;
    return 0;
}


//void crc_list_extend_mimo_blk(int** sub_codeword_list_mimo_blk, )

void crc_list_extend(int** sub_codeword_list, int** whole_codeword_list, float** pm, int L, int Nv, int Nh, int L_large) {
    int** current_codewords = NULL;
    float* current_pm = NULL;
    int current_size = 0;

    for (int i = 0; i < Nv; i++) {
        int max_candidates = (current_size == 0) ? L : current_size * L;
        int** new_codewords = (int**)malloc(max_candidates * sizeof(int*));
        float* new_pm = (float*)malloc(max_candidates * sizeof(float));
        int new_candidate_count = 0;

        if (current_size == 0) {
            // 处理第一个子码
            for (int k = 0; k < L; k++) {
                int* codeword = (int*)malloc(Nv * Nh * sizeof(int));
                memcpy(codeword, sub_codeword_list[i] + k * Nh, Nh * sizeof(int));
                new_codewords[new_candidate_count] = codeword;
                new_pm[new_candidate_count] = pm[i][k];
                new_candidate_count++;
            }
        }
        else {
            // 处理后续子码
            for (int j = 0; j < current_size; j++) {
                for (int k = 0; k < L; k++) {
                    int* codeword = (int*)malloc(Nv * Nh * sizeof(int));
                    memcpy(codeword, current_codewords[j], i * Nh * sizeof(int));
                    memcpy(codeword + i * Nh, sub_codeword_list[i] + k * Nh, Nh * sizeof(int));
                    new_codewords[new_candidate_count] = codeword;
                    new_pm[new_candidate_count] = current_pm[j] + pm[i][k];
                    new_candidate_count++;
                }
            }
        }

        // 创建索引-PM对数组并排序
        IndexPmPair* pairs = (IndexPmPair*)malloc(new_candidate_count * sizeof(IndexPmPair));
        for (int m = 0; m < new_candidate_count; m++) {
            pairs[m].index = m;
            pairs[m].pm = new_pm[m];
        }
        qsort(pairs, new_candidate_count, sizeof(IndexPmPair), compare_pairs);

        // 确定保留数量
        int retain = (new_candidate_count < L_large) ? new_candidate_count : L_large;

        // 释放旧current数据
        if (current_codewords != NULL) {
            for (int j = 0; j < current_size; j++) free(current_codewords[j]);
            free(current_codewords);
            free(current_pm);
        }

        // 分配新current数据
        current_size = retain;
        current_codewords = (int**)malloc(retain * sizeof(int*));
        current_pm = (float*)malloc(retain * sizeof(float));

        // 填充新current数据
        for (int m = 0; m < retain; m++) {
            int idx = pairs[m].index;
            current_codewords[m] = new_codewords[idx];
            current_pm[m] = new_pm[idx];
        }

        // 释放未选中的码字和临时数组
        for (int m = retain; m < new_candidate_count; m++) {
            int idx = pairs[m].index;
            free(new_codewords[idx]);
        }
        free(new_codewords);
        free(new_pm);
        free(pairs);
    }

    // 复制结果到输出数组
    for (int m = 0; m < current_size; m++) {
        memcpy(whole_codeword_list[m], current_codewords[m], Nv * Nh * sizeof(int));
        free(current_codewords[m]);
    }

    // 释放current数据
    free(current_codewords);
    free(current_pm);
}



void crc_list_extend_vertical(int** sub_codeword_list, int** whole_codeword_list, vector<int>** intrlv_pattern, float** pm, float** y, float* H_real, int L, int Nv, int Nh, int Nr, int L_large, int ModType) {
    int** current_codewords = NULL;
    float* current_pm = NULL;
    int current_size = 0;
    float channel_euclidian_distance = 0;
    int idx_sym_col = 0;

    for (int i = 0; i < Nv; i++) {
        channel_euclidian_distance = 0;
        idx_sym_col = i / ModType;
        int max_candidates = (current_size == 0) ? L : current_size * L;
        int** new_codewords = (int**)malloc(max_candidates * sizeof(int*));
        float* new_pm = (float*)malloc(max_candidates * sizeof(float));
        int new_candidate_count = 0;

        if (current_size == 0) {
            // 处理第一个子码
            for (int k = 0; k < L; k++) {
                int* codeword = (int*)malloc(Nv * Nh * sizeof(int));
                memcpy(codeword, sub_codeword_list[i] + k * Nh, Nh * sizeof(int));
                new_codewords[new_candidate_count] = codeword;
                new_pm[new_candidate_count] = pm[i][k];
                new_candidate_count++;
            }
        }
        else {
            // 处理后续子码
            for (int j = 0; j < current_size; j++) {
                for (int k = 0; k < L; k++) {
                    int* codeword = (int*)malloc(Nv * Nh * sizeof(int));
                    // conduct vertical encoding on codeword

                    memcpy(codeword, current_codewords[j], i * Nh * sizeof(int));
                    memcpy(codeword + i * Nh, sub_codeword_list[i] + k * Nh, Nh * sizeof(int));
                    int mmm = codeword[0];
                    new_codewords[new_candidate_count] = codeword;
                    if ((i + 1) % ModType == 0) { 
                        int* codeword_encoded_v = (int*)malloc(Nv * Nh * sizeof(int));
                        for (int i_col = 0; i_col < (i+1); i_col++) { PolarEncode_xor(codeword_encoded_v + i_col * Nh, codeword + i_col * Nh, Nh); }
                        int offset = idx_sym_col * ModType * Nh;
                        channel_euclidian_distance = -calcChannelEuclidianDistance(codeword_encoded_v + offset, y[idx_sym_col], H_real, intrlv_pattern, idx_sym_col, ModType, Nh, Nv, Nr);
                        free(codeword_encoded_v);
                    }
                    new_pm[new_candidate_count] = current_pm[j] + pm[i][k] + channel_euclidian_distance;
                    // cout << channel_euclidian_distance << endl;
                    new_candidate_count++;
                    // free(codeword);
                }
            }
        }

        // 创建索引-PM对数组并排序
        IndexPmPair* pairs = (IndexPmPair*)malloc(new_candidate_count * sizeof(IndexPmPair));
        for (int m = 0; m < new_candidate_count; m++) {
            pairs[m].index = m;
            pairs[m].pm = new_pm[m];
        }
        qsort(pairs, new_candidate_count, sizeof(IndexPmPair), compare_pairs);

        // 确定保留数量
        int retain = (new_candidate_count < L_large) ? new_candidate_count : L_large;

        // 释放旧current数据
        if (current_codewords != NULL) {
            for (int j = 0; j < current_size; j++) free(current_codewords[j]);
            free(current_codewords);
            free(current_pm);
        }

        // 分配新current数据
        current_size = retain;
        current_codewords = (int**)malloc(retain * sizeof(int*));
        current_pm = (float*)malloc(retain * sizeof(float));

        // 填充新current数据
       //  cout << endl;
        for (int m = 0; m < retain; m++) {
            int idx = pairs[m].index;
            current_codewords[m] = new_codewords[idx];
            current_pm[m] = new_pm[idx];
            // cout << current_pm[m] << endl;
        }
        // cout << endl;

        // 释放未选中的码字和临时数组
        for (int m = retain; m < new_candidate_count; m++) {
            int idx = pairs[m].index;
            free(new_codewords[idx]);
        }
        free(new_codewords);
        free(new_pm);
        free(pairs);
    }
   /* for (int i = 0; i < current_size; i++)
    {
        cout << endl;
        for (int j = 0; j < Nh * Nv; j++) { cout << current_codewords[i][j] << "  "; }
        cout << endl;
    }*/
    // 复制结果到输出数组
    for (int m = 0; m < current_size; m++) {
        // cout << current_size << endl;
        // cout << m << endl;
        memcpy(whole_codeword_list[m], current_codewords[m], Nv * Nh * sizeof(int));
        free(current_codewords[m]);
    }

    // 释放current数据
    // check current_pm
    free(current_codewords);
    free(current_pm);
}






void crc_list_extend_vertical_second_stage(int** sub_codeword_list, int** whole_codeword_list, vector<int>** intrlv_pattern, float** pm, float** y, float* H_real, int L, int n_mimo_blk, int len_mimo_blk, int Nr, int L_large, int ModType) {
    // Nv to n_mimo_blk;
    // Nh to ModType*Nh; (len_mimo_blk)
    int** current_codewords = NULL;
    float* current_pm = NULL;
    int current_size = 0;
    float channel_euclidian_distance = 0;
    int idx_sym_col = 0;

    for (int i = 0; i < n_mimo_blk; i++) {
        channel_euclidian_distance = 0;
        idx_sym_col = i / ModType;
        int max_candidates = (current_size == 0) ? L : current_size * L;
        int** new_codewords = (int**)malloc(max_candidates * sizeof(int*));
        float* new_pm = (float*)malloc(max_candidates * sizeof(float));
        int new_candidate_count = 0;

        if (current_size == 0) {
            // 处理第一个子码
            for (int k = 0; k < L; k++) {
                int* codeword = (int*)malloc(n_mimo_blk * len_mimo_blk * sizeof(int));
                memcpy(codeword, sub_codeword_list[i] + k * len_mimo_blk, len_mimo_blk * sizeof(int));
                new_codewords[new_candidate_count] = codeword;
                new_pm[new_candidate_count] = pm[i][k];
                new_candidate_count++;
            }
        }
        else {
            // 处理后续子码
            for (int j = 0; j < current_size; j++) {
                for (int k = 0; k < L; k++) {
                    int* codeword = (int*)malloc(n_mimo_blk * len_mimo_blk * sizeof(int));
                    // conduct vertical encoding on codeword

                    memcpy(codeword, current_codewords[j], i * len_mimo_blk * sizeof(int));
                    memcpy(codeword + i * len_mimo_blk, sub_codeword_list[i] + k * len_mimo_blk, len_mimo_blk * sizeof(int));
                    int mmm = codeword[0];
                    new_codewords[new_candidate_count] = codeword;
                    /*if ((i + 1) % ModType == 0) {
                        int* codeword_encoded_v = (int*)malloc(Nv * Nh * sizeof(int));
                        for (int i_col = 0; i_col < (i + 1); i_col++) { PolarEncode_xor(codeword_encoded_v + i_col * Nh, codeword + i_col * Nh, Nh); }
                        int offset = idx_sym_col * ModType * Nh;
                        channel_euclidian_distance = -calcChannelEuclidianDistance(codeword_encoded_v + offset, y[idx_sym_col], H_real, intrlv_pattern, idx_sym_col, ModType, Nh, Nv, Nr);
                        free(codeword_encoded_v);
                    }*/
                    new_pm[new_candidate_count] = current_pm[j] + pm[i][k] + channel_euclidian_distance;
                    // cout << channel_euclidian_distance << endl;
                    new_candidate_count++;
                    // free(codeword);
                }
            }
        }

        // 创建索引-PM对数组并排序
        IndexPmPair* pairs = (IndexPmPair*)malloc(new_candidate_count * sizeof(IndexPmPair));
        for (int m = 0; m < new_candidate_count; m++) {
            pairs[m].index = m;
            pairs[m].pm = new_pm[m];
        }
        qsort(pairs, new_candidate_count, sizeof(IndexPmPair), compare_pairs);

        // 确定保留数量
        int retain = (new_candidate_count < L_large) ? new_candidate_count : L_large;

        // 释放旧current数据
        if (current_codewords != NULL) {
            for (int j = 0; j < current_size; j++) free(current_codewords[j]);
            free(current_codewords);
            free(current_pm);
        }

        // 分配新current数据
        current_size = retain;
        current_codewords = (int**)malloc(retain * sizeof(int*));
        current_pm = (float*)malloc(retain * sizeof(float));

        // 填充新current数据
       //  cout << endl;
        for (int m = 0; m < retain; m++) {
            int idx = pairs[m].index;
            current_codewords[m] = new_codewords[idx];
            current_pm[m] = new_pm[idx];
            // cout << current_pm[m] << endl;
        }
        // cout << endl;

        // 释放未选中的码字和临时数组
        for (int m = retain; m < new_candidate_count; m++) {
            int idx = pairs[m].index;
            free(new_codewords[idx]);
        }
        free(new_codewords);
        free(new_pm);
        free(pairs);
    }
    /* for (int i = 0; i < current_size; i++)
     {
         cout << endl;
         for (int j = 0; j < Nh * Nv; j++) { cout << current_codewords[i][j] << "  "; }
         cout << endl;
     }*/
     // 复制结果到输出数组
    for (int m = 0; m < current_size; m++) {
        // cout << current_size << endl;
        // cout << m << endl;
        memcpy(whole_codeword_list[m], current_codewords[m], Nv * Nh * sizeof(int));
        free(current_codewords[m]);
    }

    /*for (int m = 0; m < current_size; m++) {
        memcpy(whole_codeword_list[idx_sym_col] + m * ModType * Nh, current_codewords[m], ModType * Nh * sizeof(int));
    }*/

    // 释放current数据
    // check current_pm
    free(current_codewords);
    free(current_pm);
}

// double stage crc list extend
// first stage: search in every symbol column (this step generates a list of size L_mimo_blk)
// second stage: aggregate the search result in stage 1 (this step further generates a list of size L_large)
void crc_list_extend_vertical_doublestage(int** sub_codeword_list, int** whole_codeword_list, vector<int>** intrlv_pattern, float** pm,
    float** y, float* H_real, int L, int Nv, int Nh, int Nr, int L_large, int L_mimo_blk, int ModType) {
    //int** current_codewords = NULL;
    //float* current_pm = NULL;
    int current_size = 0;
    float channel_euclidian_distance = 0;
    int idx_sym_col = 0;

    // array for storing the codeword list for each MIMO block
    int n_mimo_blk = Nv / ModType;
    int** sub_codeword_list_mimo = (int**)malloc(n_mimo_blk * sizeof(int*));
    float** pm_renew = (float**)malloc(n_mimo_blk * sizeof(float*));


    for (int i = 0; i < n_mimo_blk; i++) { pm_renew[i] = (float*)malloc(L_mimo_blk * sizeof(float)); }
    for (int i = 0; i < n_mimo_blk; i++) { sub_codeword_list_mimo[i] = (int*)malloc(L_mimo_blk * Nh * ModType * sizeof(int)); }

    for (int idx_sym_col = 0; idx_sym_col < n_mimo_blk; idx_sym_col++)
    {
        crc_list_extend_vertical_singlestep(sub_codeword_list, sub_codeword_list_mimo, intrlv_pattern, pm, y,
            H_real, pm_renew[idx_sym_col], L, Nv, Nh, Nr, L_mimo_blk, ModType, idx_sym_col);
    }

    // test: sub_codeword_list and sub_codeword_mimo_list
    //cout << "SUB CODEWORD LIST" << endl;
    //for (int i = 0; i < Nv; i++)
    //    for (int j = 0; j < Nh; j++)
    //    {
    //        cout << sub_codeword_list[i][j] << "  ";
    //    }
    //cout << endl << "SUB CODEWORD LIST" << endl;
    //cout << "SUB CODEWORD MIMO LIST" << endl;
    //for (int i = 0; i < n_mimo_blk; i++)
    //{
    //    for (int j = 0; j < ModType * Nh; j++)
    //    {
    //        cout << sub_codeword_list_mimo[i][j] << "  ";
    //        /*for (int k = 0; k < Nh; k++) {
    //            cout << sub_codeword_list_mimo[i][j * Nh + k] << "  ";
    //        }*/
    //    }
    //}
    //cout << endl << "SUB CODEWORD MIMO LIST" << endl;
    // test: sub_codeword_list and sub_codeword_mimo_list

    // After this step, sub_codeword_list_mimo should retain L_mimo_blk decoding paths for the n_mimo_blk mimo blocks
    // The length of each retaned decoding paths is ModType*Nh
    // To visit each decodiing path, write: sub_codeword_list_mimo[i_mimo_blk][i_path*ModType*Nh+j_bit];
    // The first stage is done, now we need to do the second stage
    crc_list_extend_vertical_second_stage(sub_codeword_list_mimo, whole_codeword_list, intrlv_pattern, pm_renew, y, H_real, L_mimo_blk, n_mimo_blk, ModType * Nh,
        Nr, L_large, ModType);

    // test: whole_codeword_list
   /* cout << "WHOLE CODEWORD LIST" << endl;
    for (int i = 0; i < Nh * Nv; i++)
    {
        cout << whole_codeword_list[0][i] << "  ";
    }
    cout << endl << "WHOLE CODEWORD LIST" << endl;*/
    // test: whole_codeword_list 

    for (int i = 0; i < n_mimo_blk; i++) { free(sub_codeword_list_mimo[i]); }
    free(sub_codeword_list_mimo);
    for (int i = 0; i < n_mimo_blk; i++) { free(pm_renew[i]); }
    free(pm_renew);
    // free(current_codewords);
    // free(current_pm);
    // free(sub_codeword_list_mimo);
    // free(pm_renew);
}




void crc_list_extend_vertical_singlestep(int** sub_codeword_list, int** whole_codeword_list, vector<int>** intrlv_pattern, float** pm, float** y, float* H_real,
    float* pm_renew, int L, int Nv, int Nh, int Nr, int L_mimo_blk, int ModType, int idx_sym_col) {
    int** current_codewords = NULL;
    float* current_pm = NULL;
    int current_size = 0;
    int sub_codeword_offset = idx_sym_col * ModType;  // n of vertical codewords before the current codeword
    float channel_euclidian_distance = 0;
    //int idx_sym_col = 0;

    // array for storing the codeword list for each MIMO block
    int n_mimo_blk = Nv / ModType;
    int** sub_codeword_list_mimo = (int**)malloc(n_mimo_blk * sizeof(int*));
    for (int i = 0; i < n_mimo_blk; i++) { sub_codeword_list_mimo[i] = (int*)malloc(L_mimo_blk * Nh * ModType * sizeof(int)); }


    for (int i = 0; i < ModType; i++) {
        channel_euclidian_distance = 0;
        //idx_sym_col = i / ModType;
        int max_candidates = (current_size == 0) ? L : current_size * L;
        int** new_codewords = (int**)malloc(max_candidates * sizeof(int*));
        float* new_pm = (float*)malloc(max_candidates * sizeof(float));
        int new_candidate_count = 0;

        if (current_size == 0) {
            // 处理第一个子码
            for (int k = 0; k < L; k++) {
                int* codeword = (int*)malloc(Nv * Nh * sizeof(int));
                memcpy(codeword, sub_codeword_list[i + sub_codeword_offset] + k * Nh, Nh * sizeof(int));
                new_codewords[new_candidate_count] = codeword;
                new_pm[new_candidate_count] = pm[i + sub_codeword_offset][k];
                new_candidate_count++;
            }
        }
        else {
            // 处理后续子码
            for (int j = 0; j < current_size; j++) {
                for (int k = 0; k < L; k++) {
                    int* codeword = (int*)malloc(Nv * Nh * sizeof(int));
                    // conduct vertical encoding on codeword

                    memcpy(codeword, current_codewords[j], i * Nh * sizeof(int));
                    memcpy(codeword + i * Nh, sub_codeword_list[i + sub_codeword_offset] + k * Nh, Nh * sizeof(int));
                    int mmm = codeword[0];
                    new_codewords[new_candidate_count] = codeword;
                    if ((i + 1) % ModType == 0) {
                        int* codeword_encoded_v = (int*)malloc(Nv * Nh * sizeof(int));
                        for (int i_col = 0; i_col < (i + 1); i_col++) { PolarEncode_xor(codeword_encoded_v + i_col * Nh, codeword + i_col * Nh, Nh); }
                        // int offset = idx_sym_col * ModType * Nh;
                        int offset = 0;
                        channel_euclidian_distance = -calcChannelEuclidianDistance(codeword_encoded_v + offset, y[idx_sym_col], H_real, intrlv_pattern, idx_sym_col, ModType, Nh, Nv, Nr);
                        free(codeword_encoded_v);
                    }
                    new_pm[new_candidate_count] = current_pm[j] + pm[i + sub_codeword_offset][k] + channel_euclidian_distance;
                    // cout << channel_euclidian_distance << endl;
                    new_candidate_count++;
                    // free(codeword);
                }
            }
        }

        // 创建索引-PM对数组并排序
        IndexPmPair* pairs = (IndexPmPair*)malloc(new_candidate_count * sizeof(IndexPmPair));
        for (int m = 0; m < new_candidate_count; m++) {
            pairs[m].index = m;
            pairs[m].pm = new_pm[m];
        }
        qsort(pairs, new_candidate_count, sizeof(IndexPmPair), compare_pairs);

        // 确定保留数量
        int retain = (new_candidate_count < L_mimo_blk) ? new_candidate_count : L_mimo_blk;

        // 释放旧current数据
        if (current_codewords != NULL) {
            for (int j = 0; j < current_size; j++) free(current_codewords[j]);
            free(current_codewords);
            free(current_pm);
        }

        // 分配新current数据
        current_size = retain;
        current_codewords = (int**)malloc(retain * sizeof(int*));
        current_pm = (float*)malloc(retain * sizeof(float));

        // 填充新current数据
       //  cout << endl;
        for (int m = 0; m < retain; m++) {
            int idx = pairs[m].index;
            current_codewords[m] = new_codewords[idx];
            current_pm[m] = new_pm[idx];
            pm_renew[m] = new_pm[idx];
            // cout << current_pm[m] << endl;
        }
        // cout << endl;

        // 释放未选中的码字和临时数组
        for (int m = retain; m < new_candidate_count; m++) {
            int idx = pairs[m].index;
            free(new_codewords[idx]);
        }
        free(new_codewords);
        free(new_pm);
        free(pairs);
    }
    /* for (int i = 0; i < current_size; i++)
     {
         cout << endl;
         for (int j = 0; j < Nh * Nv; j++) { cout << current_codewords[i][j] << "  "; }
         cout << endl;
     }*/
     // 复制结果到输出数组
    //for (int m = 0; m < current_size; m++) {
    //    // cout << current_size << endl;
    //    // cout << m << endl;
    //    memcpy(whole_codeword_list[m], current_codewords[m], Nv * Nh * sizeof(int));
    //    free(current_codewords[m]);
    //}


    for (int m = 0; m < current_size; m++) {
        memcpy(whole_codeword_list[idx_sym_col] + m * ModType * Nh, current_codewords[m], ModType * Nh * sizeof(int));
    }


    // 释放current数据
    // check current_pm
    free(current_codewords);
    free(current_pm);
}




float calculate_norm(const float* y, const float* H, const float* x, int Nr, int Nt) {
    const MKL_INT M = 2 * Nr;   // 行数
    const MKL_INT N = 2 * Nt;   // 列数

    // 分配对齐内存用于存储y - Hx的结果
    float* diff = (float*)mkl_malloc(M * sizeof(float), 64);

    // 复制y到diff
    cblas_scopy(M, y, 1, diff, 1);

    // 计算diff = 1.0*diff + (-1.0)*H*x (即y - Hx)
    cblas_sgemv(CblasRowMajor, CblasNoTrans,
        M, N,           // 矩阵维度MxN
        -1.0f, H, N,    // -1.0*H，行主序矩阵的lda=N
        x, 1,           // x向量，步长1
        1.0f,           // 保持diff原有值（1.0*diff）
        diff, 1);       // 结果存储，步长1

    // 计算L2范数
    float norm = cblas_snrm2(M, diff, 1);

    // 释放内存
    mkl_free(diff);

    return norm;
}



// calculate |y-HM(x)|, where M(·) represents modulation
float calcChannelEuclidianDistance(int* x, float* y, float* H_real, std::vector<int>** intrlv_pattern, int idx_sym_col, int ModType, int Nh, int Nv, int Nr)
{
    int n_real_symbols = 2 * Nv;
    int len_real_symbol = ModType;
    int len_half_symbol = ModType / 2;
    float distance = 0;
    float* x_mod_real = new float[n_real_symbols];
    int* x_intrlv = new int[Nv * len_real_symbol]; // i*len_real_symbol+j
    for (int i = 0; i < Nv; i++) { for (int j = 0; j < len_real_symbol; j++) { x_intrlv[i * len_real_symbol + j] = x[j * Nv + i]; } }
    // map x (bits) to symbols
    // x before interleaving
    // i th row, j th col 
    //cout << "x before interleaving" << endl;
    //for (int i = 0; i < Nv; i++)
    //{
    //    cout << endl;
    //    for (int j = 0; j < len_real_symbol; j++) {
    //        cout << x_intrlv[i * len_real_symbol + j] << "  ";
    //        // cout << x_intrlv[j * Nv + i] << "  ";
    //    }
    //    cout << endl;
    //}
    //cout << "x before interleaving" << endl;
    // Interleaving
   
    for (int i = 0; i < Nv; i++) { randIntrlv(intrlv_pattern[i][idx_sym_col], x_intrlv + i * len_real_symbol, len_real_symbol); }
    /*cout << "x after interleaving" << endl;
    for (int i = 0; i < Nv; i++)
    {
        cout << endl;
        for (int j = 0; j < len_real_symbol; j++) {
            cout << x_intrlv[i * len_real_symbol + j] << "  ";
        }
        cout << endl;
    }
    cout << "x after interleaving" << endl;*/
    // Bits to symbols mapping
    for (int i = 0; i < Nv; i++)
    {
        int tmp_i = 0;
        int tmp_q = 0;
        for (int j = 0; j < len_half_symbol; j++) 
        {
            tmp_i = tmp_i + (x_intrlv[i * len_real_symbol + j] << (len_half_symbol - j - 1));
            tmp_q = tmp_q + (x_intrlv[i * len_real_symbol + len_half_symbol + j] << (len_half_symbol - j - 1));
        }
        x_mod_real[i] = table_16qam_crc[tmp_i];
        x_mod_real[Nv + i] = table_16qam_crc[tmp_q];

        // test: modulated symbols in crc
        /*if (i == 0)
        {
            cout << "modulated symbols in CRC" << endl;
            for (int j = 0; j < 2 * Nv; j++) { cout << x_mod_real[j] << endl; }
            cout << "modulated symbols in CRC" << endl;
        }*/
    }

    /*cout << "modulated symbols in CRC" << endl;
    for (int j = 0; j < 2 * Nv; j++) { cout << x_mod_real[j] << endl; }
    cout << "modulated symbols in CRC" << endl;*/

    // calculate x_mod_real
    distance = calculate_norm(y, H_real, x_mod_real, Nr, Nt);
    delete[] x_mod_real;
    delete[] x_intrlv;
    return distance;
}