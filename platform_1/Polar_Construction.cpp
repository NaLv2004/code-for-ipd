#include "Polar_Encoder.h"

void polar_codeconstruction(string method, int N, int K, int** GG, float sigma, vector<int>& best_channel)
{
	if (method == "RM") RM_codeconstruction(N, GG, best_channel);
	if (method == "beta") beta_codeconstruction(N, best_channel);
	if (method == "GA") GA_codeconstruction(N, sigma, best_channel);
	if (method == "5G") B5G_codeconstruction(N, K, best_channel);
	if (method == "RM-GA") RM_GA_codeconstruction(N, K, GG, sigma, best_channel);
}