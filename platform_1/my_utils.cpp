# include "my_utils.h"
bool my_find(int* arr, int elem, int len)
{
	bool bfind = false;
	for (int i = 0; i < len; i++)
		if (elem == arr[i])
		{
			bfind = true;
			break;
		}
	return bfind;
}

// rearrange the output of the decoder for feedback to detector
void rearrange_dec_det(float* llr, float* llr_rearranged, int mod_type, int Nt)
{
	int Nt2 = 2 * Nt;
	int Nbits = Nt * mod_type;
	int mod_type_half = mod_type / 2;
	for (int i_sym = 0; i_sym < Nt; i_sym++)
	{
		for (int j = 0; j < mod_type; j++)
		{
			if (j < mod_type_half)
			{
				llr_rearranged[i_sym * mod_type_half + j] = llr[i_sym * mod_type + j];
			}
			else
			{
				llr_rearranged[(Nt + i_sym) * mod_type_half + j - mod_type_half] = llr[i_sym * mod_type + j];
			}
		}
	}
}
