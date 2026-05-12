# include<algorithm>
# include<iostream>
# include<iomanip>
# include "listSphereDecoding.h"
using namespace std;

/*
* function: listSphereDecode
* u: decoded bits
* y: bitwise LLR input
* A: info bit set
* N: code length
* L: list size
* PM: Surviving path metrics
* PM_extend: expanded path metrics of every decoding step (1 x 2L)
* x_hat: hard decision according to input llr y
*/
void listSphereDecode(int* u, float* y, int** G, int* A, int N, int L, int* n_paths)
{
    // TEST:
   /* cout << "Generation Matrix" << endl;
    for (int i = N - 10; i < N; i++)
    {
        cout << endl;
        for (int j = N - 10; j < N; j++)
            cout << "  " << G[i][j];
    }*/
    // TEST
    float* PM = new float[L];
    int* i_paths_sort = new int[2 * L];
    /*int** PM_extend = new int* [2];
    PM_extend[0] = new int[L];
    PM_extend[1] = new int[L];*/
    float* PM_extend = new float[2 * L];
    int** paths = new int* [L];
    for (int i = 0; i < L; i++) { paths[i] = new int[N]; }
    int** paths_extend = new int* [2 * L];
    for (int i = 0; i < 2 * L; i++) { paths_extend[i] = new int[N]; }
    int** paths_renew = new int* [L];
    for (int i = 0; i < L; i++) { paths_renew[i] = new int[N]; }
    int* x_hat = new int[N];

    // Initialization
    std::fill(PM, PM + L, 0);
    std::fill(PM_extend, PM_extend + 2*L, 0);
    for (int i = 0; i < L; i++) {
        std::fill(paths[i], paths[i] + N, 0);
    }
    for (int i = 0; i < 2 * L; i++) {
        std::fill(paths_extend[i], paths_extend[i] + N, 0);
    }
    for (int i = 0; i < N; i++) {
        if (y[i] > 0) { x_hat[i] = 0; }
        else { x_hat[i] = 1; }
    }
    // TEST:
    /*for (int i = 0; i < N; i++) cout << x_hat[i] << endl;*/
    // TEST
    // Decoding Steps
    for (int t = 0; t < N; t++)
    {
        // TEST
       /* if (t == 1)
        {
            cout << "TEST t-1" << endl;
            cout << paths_extend[0][127] << paths_extend[1][127] << endl;
        }
        if (t == 7)
        {
            cout << "TEST t-5" << endl;
        }*/
        // TEST
        int i_bit = N - 1 - t;
        // Add new paths
        // int n_paths_total = pow(2, t + 1); 
        int n_paths_survive = n_paths[t]; // surviving paths after the current decoding step
        // if (n_paths_total > L) { n_paths_total = 2 * L; n_paths_survive = L; }
        int n_paths_pre = 1;
        if (t > 0) n_paths_pre = n_paths[t - 1];
        /*for (int i_path = 0; i_path < n_paths_total; i_path++)
        {
            for (int i_bit = N - 1; i_bit >= 0; i_bit--)
            {

            }
        }*/
        // calculate all the path metrics
        for (int i_path_pre = 0; i_path_pre < n_paths_pre; i_path_pre++)
        {
            bool bSame0 = judgeSame_append(t, paths_extend[i_path_pre], x_hat, G, N, 0);
            // bool bSame1 = judgeSame_append(t, paths[i_path_pre], x_hat, G, N, 1);
            // current bit is frozen bit
            if (A[i_bit] == 0)
            {
                if (bSame0) { PM_extend[i_path_pre] = PM_extend[i_path_pre]; }
                else { PM_extend[i_path_pre] += abs(y[i_bit]); }
            }
            // current bit is information bit
            else
            {
                if (bSame0)
                {
                    PM_extend[i_path_pre] = PM_extend[i_path_pre];
                    PM_extend[i_path_pre + n_paths_pre] = PM_extend[i_path_pre] + abs(y[i_bit]);
                }
                else
                {
                    PM_extend[i_path_pre + n_paths_pre] = PM_extend[i_path_pre];
                    PM_extend[i_path_pre] = PM_extend[i_path_pre] + abs(y[i_bit]);
                }
            }
        }
        // sort the path metrics
        // Now PM_extend = [PM_renew(1 x L) 0(1 x L)]
        int n_paths_sort = 2 * n_paths_pre;
        if (A[i_bit] == 0) { n_paths_sort = n_paths_survive; }
        mySort(PM_extend, i_paths_sort, n_paths_sort);
        for (int i_path = n_paths_survive; i_path < 2 * L; i_path++) { PM_extend[i_path] = 0; }

        // Copy L paths with minimum PM from paths_extend to paths_renew
        for (int i_path = 0; i_path < n_paths_survive; i_path++)
        {
            int i_path_real = i_paths_sort[i_path];
            for (int k = i_bit; k < N; k++)
            {
                if (k == i_bit)
                {
                    if (i_path_real < n_paths_pre) { paths_renew[i_path][k] = 0; }
                    else { paths_renew[i_path][k] = 1; }
                }
                // else { paths_renew[i_path][k] = paths_extend[i_paths_sort[i_path]][k]; }
                else
                {
                    if (i_path_real < n_paths_pre) { paths_renew[i_path][k] = paths_extend[i_path_real][k]; }
                    else { paths_renew[i_path][k] = paths_extend[i_path_real-n_paths_pre][k]; }
                }
            }
        }

        // FOR TEST
        /*if (t == 1)
        { 
            cout << "i_paths_sort" << endl;
            cout << i_paths_sort[0] << i_paths_sort[1] << i_paths_sort[2] << i_paths_sort[3] <<endl;
            cout << "paths_renew" << endl;
            for (int i_path = 0; i_path < n_paths_survive; i_path++)
            {
                cout << endl;
                for (int k = i_bit; k < N; k++)
                {
                    cout << paths_renew[i_path][k];
                }
            }
            cout << endl;
            cout << "TEST END" << endl;
        }
        */
        // FOR TEST

        // Copy back to paths_extend
        for (int i_path = 0; i_path < n_paths_survive; i_path++)
            for (int k = i_bit; k < N; k++) { paths_extend[i_path][k] = paths_renew[i_path][k]; }

        // FOR TEST
        /*cout << "decoding_step:" << t << endl;
        cout << "n of surviving decoding paths:" << n_paths_survive << endl;
        cout << "Path Metrics" << endl;
        for (int i = 0; i < 2 * L; i++) { cout << setw(10) << PM_extend[i]; }
        cout << endl;
        cout << "Surviving Paths:" << endl;
        for (int i = 0; i < n_paths_survive; i++)
        {
            cout << endl;
            for (int j = 0; j < N; j++) cout << paths_extend[i][j];
        }
        cout << endl;*/
        // FOR TEST

    } // end for t
    // Decoder Output
    for (int i = 0; i < N; i++) { u[i] = paths_extend[0][i]; }
    delete[] PM;
    delete[] PM_extend;
    for (int i = 0; i < L; i++) { delete[] paths[i]; }
    delete[] paths;
    for (int i = 0; i < L; i++) { delete[] paths_renew[i]; }
    delete[] paths_renew;
    for (int i = 0; i < 2 * L; i++) { delete[] paths_extend[i]; }
    delete[] paths_extend;
    delete[] x_hat;
    delete[] i_paths_sort;
}

// err_step: 1 x N
// ouputs the decoding result step by step
void listSphereDecode_step(int* u, float* y, int** G, int* A, int N, int L, int* n_paths)
{
    // TEST:
   /* cout << "Generation Matrix" << endl;
    for (int i = N - 10; i < N; i++)
    {
        cout << endl;
        for (int j = N - 10; j < N; j++)
            cout << "  " << G[i][j];
    }*/
    // TEST
    float* PM = new float[L];
    int* i_paths_sort = new int[2 * L];
    /*int** PM_extend = new int* [2];
    PM_extend[0] = new int[L];
    PM_extend[1] = new int[L];*/
    float* PM_extend = new float[2 * L];
    int** paths = new int* [L];
    for (int i = 0; i < L; i++) { paths[i] = new int[N]; }
    int** paths_extend = new int* [2 * L];
    for (int i = 0; i < 2 * L; i++) { paths_extend[i] = new int[N]; }
    int** paths_renew = new int* [L];
    for (int i = 0; i < L; i++) { paths_renew[i] = new int[N]; }
    int* x_hat = new int[N];

    // Initialization
    std::fill(PM, PM + L, 0);
    std::fill(PM_extend, PM_extend + 2 * L, 0);
    for (int i = 0; i < L; i++) {
        std::fill(paths[i], paths[i] + N, 0);
    }
    for (int i = 0; i < 2 * L; i++) {
        std::fill(paths_extend[i], paths_extend[i] + N, 0);
    }
    for (int i = 0; i < N; i++) {
        if (y[i] > 0) { x_hat[i] = 0; }
        else { x_hat[i] = 1; }
    }
    // TEST:
    /*for (int i = 0; i < N; i++) cout << x_hat[i] << endl;*/
    // TEST
    // Decoding Steps
    for (int t = 0; t < N; t++)
    {
        // TEST
       /* if (t == 1)
        {
            cout << "TEST t-1" << endl;
            cout << paths_extend[0][127] << paths_extend[1][127] << endl;
        }
        if (t == 7)
        {
            cout << "TEST t-5" << endl;
        }*/
        // TEST
        int i_bit = N - 1 - t;
        // Add new paths
        // int n_paths_total = pow(2, t + 1); 
        int n_paths_survive = n_paths[t]; // surviving paths after the current decoding step
        // if (n_paths_total > L) { n_paths_total = 2 * L; n_paths_survive = L; }
        int n_paths_pre = 1;
        if (t > 0) n_paths_pre = n_paths[t - 1];
        /*for (int i_path = 0; i_path < n_paths_total; i_path++)
        {
            for (int i_bit = N - 1; i_bit >= 0; i_bit--)
            {

            }
        }*/
        // calculate all the path metrics
        for (int i_path_pre = 0; i_path_pre < n_paths_pre; i_path_pre++)
        {
            bool bSame0 = judgeSame_append(t, paths_extend[i_path_pre], x_hat, G, N, 0);
            // bool bSame1 = judgeSame_append(t, paths[i_path_pre], x_hat, G, N, 1);
            // current bit is frozen bit
            if (A[i_bit] == 0)
            {
                if (bSame0) { PM_extend[i_path_pre] = PM_extend[i_path_pre]; }
                else { PM_extend[i_path_pre] += abs(y[i_bit]); }
            }
            // current bit is information bit
            else
            {
                if (bSame0)
                {
                    PM_extend[i_path_pre] = PM_extend[i_path_pre];
                    PM_extend[i_path_pre + n_paths_pre] = PM_extend[i_path_pre] + abs(y[i_bit]);
                }
                else
                {
                    PM_extend[i_path_pre + n_paths_pre] = PM_extend[i_path_pre];
                    PM_extend[i_path_pre] = PM_extend[i_path_pre] + abs(y[i_bit]);
                }
            }
        }
        // sort the path metrics
        // Now PM_extend = [PM_renew(1 x L) 0(1 x L)]
        int n_paths_sort = 2 * n_paths_pre;
        if (A[i_bit] == 0) { n_paths_sort = n_paths_survive; }
        mySort(PM_extend, i_paths_sort, n_paths_sort);
        for (int i_path = n_paths_survive; i_path < 2 * L; i_path++) { PM_extend[i_path] = 0; }

        // Copy L paths with minimum PM from paths_extend to paths_renew
        for (int i_path = 0; i_path < n_paths_survive; i_path++)
        {
            int i_path_real = i_paths_sort[i_path];
            for (int k = i_bit; k < N; k++)
            {
                if (k == i_bit)
                {
                    if (i_path_real < n_paths_pre) { paths_renew[i_path][k] = 0; }
                    else { paths_renew[i_path][k] = 1; }
                }
                // else { paths_renew[i_path][k] = paths_extend[i_paths_sort[i_path]][k]; }
                else
                {
                    if (i_path_real < n_paths_pre) { paths_renew[i_path][k] = paths_extend[i_path_real][k]; }
                    else { paths_renew[i_path][k] = paths_extend[i_path_real - n_paths_pre][k]; }
                }
            }
        }

        // FOR TEST
        /*if (t == 1)
        {
            cout << "i_paths_sort" << endl;
            cout << i_paths_sort[0] << i_paths_sort[1] << i_paths_sort[2] << i_paths_sort[3] <<endl;
            cout << "paths_renew" << endl;
            for (int i_path = 0; i_path < n_paths_survive; i_path++)
            {
                cout << endl;
                for (int k = i_bit; k < N; k++)
                {
                    cout << paths_renew[i_path][k];
                }
            }
            cout << endl;
            cout << "TEST END" << endl;
        }
        */
        // FOR TEST

        // Copy back to paths_extend
        for (int i_path = 0; i_path < n_paths_survive; i_path++)
            for (int k = i_bit; k < N; k++) { paths_extend[i_path][k] = paths_renew[i_path][k]; }

        // FOR TEST
        /*cout << "decoding_step:" << t << endl;
        cout << "n of surviving decoding paths:" << n_paths_survive << endl;
        cout << "Path Metrics" << endl;
        for (int i = 0; i < 2 * L; i++) { cout << setw(10) << PM_extend[i]; }
        cout << endl;
        cout << "Surviving Paths:" << endl;
        for (int i = 0; i < n_paths_survive; i++)
        {
            cout << endl;
            for (int j = 0; j < N; j++) cout << paths_extend[i][j];
        }
        cout << endl;*/
        // FOR TEST
        /*if (!((t + 1) % 32))
        {
            u[i_bit] = paths_extend[0][i_bit];
            u[i_bit + 1] = paths_extend[0][i_bit + 1];
            u[i_bit + 2] = paths_extend[0][i_bit + 2];
            u[i_bit + 3] = paths_extend[0][i_bit + 3];
            u[i_bit + 4] = paths_extend[0][i_bit + 4];
            u[i_bit + 5] = paths_extend[0][i_bit + 5];
            u[i_bit + 6] = paths_extend[0][i_bit + 6];
            u[i_bit + 7] = paths_extend[0][i_bit + 7];
            u[i_bit + 8] = paths_extend[0][i_bit + 8];
            u[i_bit + 9] = paths_extend[0][i_bit + 9];
            u[i_bit + 10] = paths_extend[0][i_bit + 10];
            u[i_bit + 11] = paths_extend[0][i_bit + 11];
            u[i_bit + 12] = paths_extend[0][i_bit + 12];
            u[i_bit + 13] = paths_extend[0][i_bit + 13];
            u[i_bit + 14] = paths_extend[0][i_bit + 14];
            u[i_bit + 15] = paths_extend[0][i_bit + 15];
            u[i_bit + 16] = paths_extend[0][i_bit + 16];
            u[i_bit + 17] = paths_extend[0][i_bit + 17];
            u[i_bit + 18] = paths_extend[0][i_bit + 18];
            u[i_bit + 19] = paths_extend[0][i_bit + 19];
            u[i_bit + 20] = paths_extend[0][i_bit + 20];
            u[i_bit + 21] = paths_extend[0][i_bit + 21];
            u[i_bit + 22] = paths_extend[0][i_bit + 22];
            u[i_bit + 23] = paths_extend[0][i_bit + 23];
            u[i_bit + 24] = paths_extend[0][i_bit + 24];
            u[i_bit + 25] = paths_extend[0][i_bit + 25];
            u[i_bit + 26] = paths_extend[0][i_bit + 26];
            u[i_bit + 27] = paths_extend[0][i_bit + 27];
            u[i_bit + 28] = paths_extend[0][i_bit + 28];
            u[i_bit + 29] = paths_extend[0][i_bit + 29];
            u[i_bit + 30] = paths_extend[0][i_bit + 30];
            u[i_bit + 31] = paths_extend[0][i_bit + 31];
        }*/
        
    } // end for t
    // Decoder Output
    // for (int i = 0; i < N; i++) { u[i] = paths_extend[0][i]; }
    delete[] PM;
    delete[] PM_extend;
    for (int i = 0; i < L; i++) { delete[] paths[i]; }
    delete[] paths;
    for (int i = 0; i < L; i++) { delete[] paths_renew[i]; }
    delete[] paths_renew;
    for (int i = 0; i < 2 * L; i++) { delete[] paths_extend[i]; }
    delete[] paths_extend;
    delete[] x_hat;
    delete[] i_paths_sort;
}

/*
* function: LSD_Decode_outall
* u: decoded bits
* y: bitwise LLR input
* A: info bit set
* N: code length
* L: list size
* PM: Surviving path metrics
* PM_extend: expanded path metrics of every decoding step (1 x 2L)
* x_hat: hard decision according to input llr y
*/
void LSD_decode_outall(int* u, int** paths_extend, float* y, int** G, int* A, int N, int L, int* n_paths)
{
    // TEST:
   /* cout << "Generation Matrix" << endl;
    for (int i = N - 10; i < N; i++)
    {
        cout << endl;
        for (int j = N - 10; j < N; j++)
            cout << "  " << G[i][j];
    }*/
    // TEST
    float* PM = new float[L];
    int* i_paths_sort = new int[2 * L];
    /*int** PM_extend = new int* [2];
    PM_extend[0] = new int[L];
    PM_extend[1] = new int[L];*/
    float* PM_extend = new float[2 * L];
    int** paths = new int* [L];
    for (int i = 0; i < L; i++) { paths[i] = new int[N]; }
    /*int** paths_extend = new int* [2 * L];
    for (int i = 0; i < 2 * L; i++) { paths_extend[i] = new int[N]; }*/
    int** paths_renew = new int* [L];
    for (int i = 0; i < L; i++) { paths_renew[i] = new int[N]; }
    int* x_hat = new int[N];

    // Initialization
    std::fill(PM, PM + L, 0);
    std::fill(PM_extend, PM_extend + 2 * L, 0);
    for (int i = 0; i < L; i++) {
        std::fill(paths[i], paths[i] + N, 0);
    }
    for (int i = 0; i < 2 * L; i++) {
        std::fill(paths_extend[i], paths_extend[i] + N, 0);
    }
    for (int i = 0; i < N; i++) {
        if (y[i] > 0) { x_hat[i] = 0; }
        else { x_hat[i] = 1; }
    }
    // TEST:
    /*for (int i = 0; i < N; i++) cout << x_hat[i] << endl;*/
    // TEST
    // Decoding Steps
    for (int t = 0; t < N; t++)
    {
        // TEST
       /* if (t == 1)
        {
            cout << "TEST t-1" << endl;
            cout << paths_extend[0][127] << paths_extend[1][127] << endl;
        }
        if (t == 7)
        {
            cout << "TEST t-5" << endl;
        }*/
        // TEST
        int i_bit = N - 1 - t;
        // Add new paths
        // int n_paths_total = pow(2, t + 1); 
        int n_paths_survive = n_paths[t]; // surviving paths after the current decoding step
        // if (n_paths_total > L) { n_paths_total = 2 * L; n_paths_survive = L; }
        int n_paths_pre = 1;
        if (t > 0) n_paths_pre = n_paths[t - 1];
        /*for (int i_path = 0; i_path < n_paths_total; i_path++)
        {
            for (int i_bit = N - 1; i_bit >= 0; i_bit--)
            {

            }
        }*/
        // calculate all the path metrics
        for (int i_path_pre = 0; i_path_pre < n_paths_pre; i_path_pre++)
        {
            bool bSame0 = judgeSame_append(t, paths_extend[i_path_pre], x_hat, G, N, 0);
            // bool bSame1 = judgeSame_append(t, paths[i_path_pre], x_hat, G, N, 1);
            // current bit is frozen bit
            if (A[i_bit] == 0)
            {
                if (bSame0) { PM_extend[i_path_pre] = PM_extend[i_path_pre]; }
                else { PM_extend[i_path_pre] += abs(y[i_bit]); }
            }
            // current bit is information bit
            else
            {
                if (bSame0)
                {
                    PM_extend[i_path_pre] = PM_extend[i_path_pre];
                    PM_extend[i_path_pre + n_paths_pre] = PM_extend[i_path_pre] + abs(y[i_bit]);
                }
                else
                {
                    PM_extend[i_path_pre + n_paths_pre] = PM_extend[i_path_pre];
                    PM_extend[i_path_pre] = PM_extend[i_path_pre] + abs(y[i_bit]);
                }
            }
        }
        // sort the path metrics
        // Now PM_extend = [PM_renew(1 x L) 0(1 x L)]
        int n_paths_sort = 2 * n_paths_pre;
        if (A[i_bit] == 0) { n_paths_sort = n_paths_survive; }
        mySort(PM_extend, i_paths_sort, n_paths_sort);
        for (int i_path = n_paths_survive; i_path < 2 * L; i_path++) { PM_extend[i_path] = 0; }

        // Copy L paths with minimum PM from paths_extend to paths_renew
        for (int i_path = 0; i_path < n_paths_survive; i_path++)
        {
            int i_path_real = i_paths_sort[i_path];
            for (int k = i_bit; k < N; k++)
            {
                if (k == i_bit)
                {
                    if (i_path_real < n_paths_pre) { paths_renew[i_path][k] = 0; }
                    else { paths_renew[i_path][k] = 1; }
                }
                // else { paths_renew[i_path][k] = paths_extend[i_paths_sort[i_path]][k]; }
                else
                {
                    if (i_path_real < n_paths_pre) { paths_renew[i_path][k] = paths_extend[i_path_real][k]; }
                    else { paths_renew[i_path][k] = paths_extend[i_path_real - n_paths_pre][k]; }
                }
            }
        }

        // FOR TEST
        /*if (t == 1)
        {
            cout << "i_paths_sort" << endl;
            cout << i_paths_sort[0] << i_paths_sort[1] << i_paths_sort[2] << i_paths_sort[3] <<endl;
            cout << "paths_renew" << endl;
            for (int i_path = 0; i_path < n_paths_survive; i_path++)
            {
                cout << endl;
                for (int k = i_bit; k < N; k++)
                {
                    cout << paths_renew[i_path][k];
                }
            }
            cout << endl;
            cout << "TEST END" << endl;
        }
        */
        // FOR TEST

        // Copy back to paths_extend
        for (int i_path = 0; i_path < n_paths_survive; i_path++)
            for (int k = i_bit; k < N; k++) { paths_extend[i_path][k] = paths_renew[i_path][k]; }

        // FOR TEST
        /*if (t == N - 1) {
            cout << "decoding_step:" << t << endl;
            cout << "n of surviving decoding paths:" << n_paths_survive << endl;
            cout << "Path Metrics" << endl;
            for (int i = 0; i < 2 * L; i++) { cout << setw(10) << PM_extend[i]; }
            cout << endl;
            cout << "Surviving Paths:" << endl;
            for (int i = 0; i < n_paths_survive; i++)
            {
                cout << endl;
                for (int j = 0; j < N; j++) cout << paths_extend[i][j];
            }
            cout << endl;
        }*/
        // FOR TEST

    } // end for t
    // Decoder Output
    for (int i = 0; i < N; i++) { u[i] = paths_extend[0][i]; }
    // test PM
   /* cout << endl;
    for (int i = 0; i < L; i++)
        cout << PM_extend[i] << "    ";
    cout << endl;*/
    // test PM
    // test y
  /*  cout << endl;
    for (int i = 0; i < 32; i++)
        cout << y[i] << "    ";
        cout << endl;*/
    // test y
    delete[] PM;
    delete[] PM_extend;
    for (int i = 0; i < L; i++) { delete[] paths[i]; }
    delete[] paths;
    for (int i = 0; i < L; i++) { delete[] paths_renew[i]; }
    delete[] paths_renew;
    /*for (int i = 0; i < 2 * L; i++) { delete[] paths_extend[i]; }
    delete[] paths_extend;*/
    delete[] x_hat;
    delete[] i_paths_sort;
}






/*
* function: LSD_Decode_outall
* u: decoded bits
* y: bitwise LLR input
* A: info bit set
* N: code length
* L: list size
* PM: Surviving path metrics
* PM_extend: expanded path metrics of every decoding step (1 x 2L)
* x_hat: hard decision according to input llr y
*/
void StagedLSD_decode_outall(int* u, int** paths_extend, float* y, int** G, int* A, int N, int L, int* n_paths)
{
    // TEST:
   /* cout << "Generation Matrix" << endl;
    for (int i = N - 10; i < N; i++)
    {
        cout << endl;
        for (int j = N - 10; j < N; j++)
            cout << "  " << G[i][j];
    }*/
    // TEST
    float* PM = new float[L];
    int* i_paths_sort = new int[2 * L];
    /*int** PM_extend = new int* [2];
    PM_extend[0] = new int[L];
    PM_extend[1] = new int[L];*/
    float* PM_extend = new float[2 * L];
    int** paths = new int* [L];
    for (int i = 0; i < L; i++) { paths[i] = new int[N]; }
    /*int** paths_extend = new int* [2 * L];
    for (int i = 0; i < 2 * L; i++) { paths_extend[i] = new int[N]; }*/
    int** paths_renew = new int* [L];
    for (int i = 0; i < L; i++) { paths_renew[i] = new int[N]; }
    int* x_hat = new int[N];

    // Initialization
    std::fill(PM, PM + L, 0);
    std::fill(PM_extend, PM_extend + 2 * L, 0);
    for (int i = 0; i < L; i++) {
        std::fill(paths[i], paths[i] + N, 0);
    }
    for (int i = 0; i < 2 * L; i++) {
        std::fill(paths_extend[i], paths_extend[i] + N, 0);
    }
    for (int i = 0; i < N; i++) {
        if (y[i] > 0) { x_hat[i] = 0; }
        else { x_hat[i] = 1; }
    }
    // TEST:
    /*for (int i = 0; i < N; i++) cout << x_hat[i] << endl;*/
    // TEST
    // Decoding Steps
    for (int t = 0; t < N; t++)
    {
        // TEST
       /* if (t == 1)
        {
            cout << "TEST t-1" << endl;
            cout << paths_extend[0][127] << paths_extend[1][127] << endl;
        }
        if (t == 7)
        {
            cout << "TEST t-5" << endl;
        }*/
        // TEST
        int i_bit = N - 1 - t;
        // Add new paths
        // int n_paths_total = pow(2, t + 1); 
        int n_paths_survive = n_paths[t]; // surviving paths after the current decoding step
        // if (n_paths_total > L) { n_paths_total = 2 * L; n_paths_survive = L; }
        int n_paths_pre = 1;
        if (t > 0) n_paths_pre = n_paths[t - 1];
        /*for (int i_path = 0; i_path < n_paths_total; i_path++)
        {
            for (int i_bit = N - 1; i_bit >= 0; i_bit--)
            {

            }
        }*/
        // calculate all the path metrics
        for (int i_path_pre = 0; i_path_pre < n_paths_pre; i_path_pre++)
        {
            bool bSame0 = judgeSame_append(t, paths_extend[i_path_pre], x_hat, G, N, 0);
            // bool bSame1 = judgeSame_append(t, paths[i_path_pre], x_hat, G, N, 1);
            // current bit is frozen bit
            if (A[i_bit] == 0)
            {
                if (bSame0) { PM_extend[i_path_pre] = PM_extend[i_path_pre]; }
                else { PM_extend[i_path_pre] += abs(y[i_bit]); }
            }
            // current bit is information bit
            else
            {
                if (bSame0)
                {
                    PM_extend[i_path_pre] = PM_extend[i_path_pre];
                    PM_extend[i_path_pre + n_paths_pre] = PM_extend[i_path_pre] + abs(y[i_bit]);
                }
                else
                {
                    PM_extend[i_path_pre + n_paths_pre] = PM_extend[i_path_pre];
                    PM_extend[i_path_pre] = PM_extend[i_path_pre] + abs(y[i_bit]);
                }
            }
        }
        // sort the path metrics
        // Now PM_extend = [PM_renew(1 x L) 0(1 x L)]
        int n_paths_sort = 2 * n_paths_pre;
        if (A[i_bit] == 0) { n_paths_sort = n_paths_survive; }
        mySort(PM_extend, i_paths_sort, n_paths_sort);
        for (int i_path = n_paths_survive; i_path < 2 * L; i_path++) { PM_extend[i_path] = 0; }

        // Copy L paths with minimum PM from paths_extend to paths_renew
        for (int i_path = 0; i_path < n_paths_survive; i_path++)
        {
            int i_path_real = i_paths_sort[i_path];
            for (int k = i_bit; k < N; k++)
            {
                if (k == i_bit)
                {
                    if (i_path_real < n_paths_pre) { paths_renew[i_path][k] = 0; }
                    else { paths_renew[i_path][k] = 1; }
                }
                // else { paths_renew[i_path][k] = paths_extend[i_paths_sort[i_path]][k]; }
                else
                {
                    if (i_path_real < n_paths_pre) { paths_renew[i_path][k] = paths_extend[i_path_real][k]; }
                    else { paths_renew[i_path][k] = paths_extend[i_path_real - n_paths_pre][k]; }
                }
            }
        }

        // FOR TEST
        /*if (t == 1)
        {
            cout << "i_paths_sort" << endl;
            cout << i_paths_sort[0] << i_paths_sort[1] << i_paths_sort[2] << i_paths_sort[3] <<endl;
            cout << "paths_renew" << endl;
            for (int i_path = 0; i_path < n_paths_survive; i_path++)
            {
                cout << endl;
                for (int k = i_bit; k < N; k++)
                {
                    cout << paths_renew[i_path][k];
                }
            }
            cout << endl;
            cout << "TEST END" << endl;
        }
        */
        // FOR TEST

        // Copy back to paths_extend
        for (int i_path = 0; i_path < n_paths_survive; i_path++)
            for (int k = i_bit; k < N; k++) { paths_extend[i_path][k] = paths_renew[i_path][k]; }

        // FOR TEST
        /*if (t == N - 1) {
            cout << "decoding_step:" << t << endl;
            cout << "n of surviving decoding paths:" << n_paths_survive << endl;
            cout << "Path Metrics" << endl;
            for (int i = 0; i < 2 * L; i++) { cout << setw(10) << PM_extend[i]; }
            cout << endl;
            cout << "Surviving Paths:" << endl;
            for (int i = 0; i < n_paths_survive; i++)
            {
                cout << endl;
                for (int j = 0; j < N; j++) cout << paths_extend[i][j];
            }
            cout << endl;
        }*/
        // FOR TEST

    } // end for t
    // Decoder Output
    for (int i = 0; i < N; i++) { u[i] = paths_extend[0][i]; }
    delete[] PM;
    delete[] PM_extend;
    for (int i = 0; i < L; i++) { delete[] paths[i]; }
    delete[] paths;
    for (int i = 0; i < L; i++) { delete[] paths_renew[i]; }
    delete[] paths_renew;
    /*for (int i = 0; i < 2 * L; i++) { delete[] paths_extend[i]; }
    delete[] paths_extend;*/
    delete[] x_hat;
    delete[] i_paths_sort;
}












bool judgeSame(int t, int* u, int* x, int** G, int N)
{
    int x_decod = 0;
    bool bSame = 0;
    int j = N - t - 1;
    for (int k = j; k < N; k++)
    {
        x_decod = x_decod ^ (G[k][j] & u[k]);
    }
    if (x_decod == u[j]) { bSame = 1; }
    return bSame;
}

bool judgeSame_append(int t, int* u, int* x, int** G, int N, int bit_append)
{
    int x_decod = bit_append;
    bool bSame = 0;
    int j = N - t -1;
    for (int k = j; k < N; k++)
    {
        x_decod = x_decod ^ (G[k][j] & u[k]);
    }
    if (x_decod == x[j]) { bSame = 1; }
    return bSame;
}

void calc_npaths(int* A, int N, int* n_paths, int L)
{
    int n_path = 1;
    for (int t = 0; t < N; t++)
    {
        int i_bit = N - 1 - t;
        if (A[i_bit] == 1) { n_path *= 2; }
        if (n_path > L) { n_path = L; }
        n_paths[t] = n_path;
    }
}

void mySort(float* arr, int* idx_sort, int N) {
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

