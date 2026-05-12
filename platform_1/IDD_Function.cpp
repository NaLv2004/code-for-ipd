# include "IDD_Function.h"
# include <cmath>
void calc_llr_idd(float* llr_out, float** llr_in, int** GGh, int decoding_col_curr, int Nh, int Nv)
{
	for (int i = 0; i < Nv; i++) llr_out[i] = 0;
	int i_row_min = decoding_col_curr + 1;
	for (int i_code = 0; i_code < Nv; i_code++)
	{
		float p0 = 0;
		float p1 = 0;
		float p0_prod = 1;
		for (int i_bit = i_row_min; i_bit < Nh; i_bit++)
		{
			if (GGh[i_bit][decoding_col_curr]) 
			{ 
				p0_prod *= (1 - 2 /(1 + exp(llr_in[i_code][i_bit])));
			}
			// p0_prod *= GGh[decoding_col_curr][i_bit];
		}
		p0 = 0.5 * (1 + p0_prod);
		p1 = 1 - p0;
		float logP01 = log(p0 / p1);
		llr_out[i_code] = logP01;
	}
}

