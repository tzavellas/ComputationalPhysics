#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <TFile.h>
#include <TH1F.h>
#include <TRandom.h>

void partialSums(const double Lmin,const double Lmax, int N, double& meanSum, double& sigmaSum)
{
	TRandom3 ran(0);								// Random Number Generator
	meanSum = 0.;									// Variables that accumulate
	sigmaSum= 0.;									// Monte Carlo sums
	for(int i = 0; i < N; ++i)
	{
		double r = Lmin + (Lmax - Lmin)*ran.Rndm(); // Generate r in U(Lmin, Lmax)
		double fr = (r*r)*(1 + exp(r*r));			// Evaluate integrant at r
		meanSum += fr;								// Update sums
		sigmaSum += fr*fr;
	}
}

void StratifiedIntegration(const double Lmin, const double Lmax, const int N, const int strata, double& integ, double& err)
{
	double r_step = (Lmax - Lmin) / strata;
	double r_min = Lmin;
	double r_max = r_step;
	double meanSum = 0.;
	double sigmaSum = 0.;
	
	for (int i = 0; i < strata; ++i)
	{
		double ithMeanSum, ithSigmaSum;												// ith mean and sigma sum variables
		const int ithN = N/strata;													// each block gets equal amount of points
		partialSums(r_min, r_max, ithN, ithMeanSum, ithSigmaSum); 					// Evaluate mean and sigma^2 sums for the ith block
		meanSum += (ithMeanSum / ithN) * (r_max - r_min);							// Accumulate all block Integrals
		sigmaSum += sqrt(ithSigmaSum/ithN - (ithMeanSum/ithN)*(ithMeanSum/ithN));	// Accumulate all block variances
		r_min += r_step;															// Move to next block
		r_max += r_step;
	}
	integ = meanSum;																// Sum of all block integrals is the final integral
	err = sigmaSum * (r_max - r_min) / sqrt(N);										// Sum of all variances divided by N times (Lmax - Lmin) is the final error
}

void mc_MomentOfInertia()
{
	const double L = 4.;				// Rod length in cm
	const int N = 100000;				// Number of samples
	double momentInertia;				// Value of moment of intertia
	double momentInertiaErr;			// Value of moment of interia error
	
	// Stratified Integration in 1 part -> Crude Monte Carlo
	StratifiedIntegration(0., L, N, 1, momentInertia, momentInertiaErr);
	printf("\nCrude Monte Carlo (%d samples)\n", N);
	printf("momentInertia = %3.2e +- %4.2e\n", momentInertia, momentInertiaErr);
	
	// Stratified Integration in 4 strata
	int strata = 4.;
	StratifiedIntegration(0., L, N, strata, momentInertia, momentInertiaErr);
	printf("\nStratified Monte Carlo (%d strata)\n", strata);
	printf("momentInertia = %3.2e +- %4.2e\n", momentInertia, momentInertiaErr);
	
	// Stratified Integration in 16 strata
	strata = 16.;
	StratifiedIntegration(0., L, N, strata, momentInertia, momentInertiaErr);
	printf("\nStratified Monte Carlo (%d strata)\n", strata);
	printf("momentInertia = %3.2e +- %4.2e\n", momentInertia, momentInertiaErr);
	
	// Analytic Values for Comparison	
	const double M = 1.71975e7;
	printf("\nAnalytic Evaluation for Comparison\n", strata);
	printf("momentInertia = %3.2e\n", M);
}
