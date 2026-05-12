#include "setting.h"
#include <string>
#include<iomanip>
using namespace std;
void printsetting()
{
	cout << "------------------------Simulation Params-----------------------------" << endl;
	cout << dimensions << "D Coding ";
	if (atten == 1) cout << "SISO"; else cout << "MIMO";
	cout << " Systems\t"; cout << "\tError Frame = " << errframe << "\n\n";
	cout << "------------------------Polar Cooding Params--------------------------" << endl;
	if (dimensions == 1) 	cout << "N = " << N << ", K = " << K << ", Polar Construction = " << polarcon << "\n\n";
	else
	{
		cout << "Nv = " << Nv << ", Nh = " << Nh << ", Kv = " << Kv << ", Kh = " << Kh << ", Polar Construction vertical = " << polarconv << ", Polar Construction horizon = " << polarconh << "\n\n";
		if (sys) { cout << "sys"; }
		else { cout << "non-sys"; }
		cout << endl;
		if (irregular_coding)
		{
			cout << "Irregular Coding" << endl;
			cout << "1st Component Code: " << K_time_array[0] << endl;
			cout << "2nd Component Code: " << K_time_array[1] << endl;
			cout << "Border: " << border << endl;
			if (mimo_rate_matching_power_alloc) {
				cout << "Rate matching power allocation: " << mimo_rate_matching_enlarge_ratio << endl;
			}
			if (MIMO_SVD_SEPARATE_AVERAGE_POWER_ALLOC) {
				cout << "Separate Average Power Allocation: " << MIMO_SVD_SEPARATE_POWER_RATIO << endl;
			}
			if (mimo_svd_average_power_alloc) {
				cout << "MIMO Average Power Allocation" << endl;
			}
		}
		cout << "------------------------Polar Decoding Params-------------------------" << endl;
		cout << "CRCL = " << CRCL << ", SCL_List = " << L << endl << endl;
		if (polarde == "SLD") { cout << "CRCL = " << CRCL << ", LSD_List = " << L_LSD << endl << endl; }
		if (dimensions == 2 && polarde == "SLD")
		{
			cout << "Temporal(Horizontal) Decoder";
			for (int i = 0; i < softiter; i++)
			{
				cout << "  " << polarde_array[i];
			}
			cout << endl;
		}
		if (CRCL) cout << "Whole Codeword List= " << L_WHOLE_CODEWORD << endl;
		cout << "------------------------2D Iteration Params----------------------------" << endl;
		cout << "Soft iterations = " << softiter << endl;

		if (softmethod == "Pyndiah")
		{
			cout << "Pyndiah\n";
			cout << "alpha[" << alphalen << "] = ";
			for (int i = 0; i < alphalen; i++)
				cout << alpha[i] << ", ";
			cout << endl;
			cout << "usebeta = " << usebeta << endl;
			if (usebeta)
			{
				if (adaptive_beta)
				{
					cout << "Adaptive Beta" << endl;
					cout << " a_beta[" << betalen << "] = ";
					for (int i = 0; i < betalen; i++)
						cout << a_beta[i] << ", ";
					cout << endl;
					cout << " b_beta[" << betalen << "] = ";
					for (int i = 0; i < betalen; i++)
						cout << b_beta[i] << ", ";
					cout << endl;
				}
				else
				{
					cout << " beta[" << betalen << "] = ";
					for (int i = 0; i < betalen; i++)
						cout << beta[i] << ", ";
					cout << endl;
				}

			}

		}
		if (softmethod == "MITSO")
		{
			cout << "MITSO\n";
			cout << "alpha = " << SOalpha << endl;
		}
		if (half_iter == 1) cout << "Final Decoding Iteration is Incomplete" << endl;
		if (!half_iter) cout << "Final Decoding Iteration is Complete" << endl;
		if (weirditer) cout << "Weird Iteration" << endl;
	}
	if (atten != 1)
	{
		cout << "------------------------MIMO Params------------------------------------" << endl;
		int Q = pow(2, ModType);
		cout << Q << "QAM" << endl;
		cout << "Tx:" << setw(5) << Nt << "   Rx:" << setw(5) << Nr << endl;
		cout << MOD_DIM << "  modulation" << endl;
		cout << "MIMO Detection Algorithm:  " << MIMODET << endl;
		if (MIMODET == "KBEST") { cout << "K of KBest=  " << K_KBEST << endl; }
		cout << "MIMO LLR THRES = " << mimo_llr_thres << endl;
		if (system_archi == "JDD")
		{
			cout << "JDD System" << endl;
			cout << "JDD Coefficient= " << beta_JDD << endl;
		}
		if (system_archi != "JDD") { cout << "SDD System" << endl; }
		cout << "Interleave=  " << flag_interleave << endl;
		cout << "Full Interleave= " << flag_interleave_all << endl;
		cout << "Partial S-T Interleave" << flag_interleave_partial << endl;
		if (non_uniform_channel) {
			cout << "Non-Uniform Channel" << endl;
			cout << "Diff Gain= " << diff_gain;
		}
	}
	
}