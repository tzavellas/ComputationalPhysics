#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <TFile.h>
#include <TH1F.h>
#include <TRandom.h>


void mc_rejection()
{
	TRandom3 ran(0);					// Random Number Generator
	const int N = 100000;				// Number of samples
	const double x_max = 2.;			// Max Limit for x variable
	const double y_max = 2.;			// Max Limit for y variable
	const double area = x_max * y_max;	// Max area
	std::vector<double> x, y;			// Store the values for plotting the points
	double mass, x_cm, y_cm, errorMass, errorXcm, errorYcm;
	double sumMass = 0, sumXcm = 0, sumYcm = 0, sumSigmaMass = 0, sumSigmaXcm = 0, sumSigmaYcm = 0;
	
	for(int i = 0; i < N; i++)
	{
		double xi = x_max * ran.Rndm();					// Generate x point in [0,2]
		double yi = y_max * ran.Rndm(); 				// Generate y point in [0,2]
		
		if( (yi + xi < 2.) && (xi*xi + yi*yi > 1.) ) // If pair "inside" object
		{
			x.push_back(xi);							// Store pair in vectors for later plot
			y.push_back(yi);

			double density = 1. + xi*xi + yi*yi; 		// surface density around random pair
			sumMass += density;							// Accumulate mass
			sumXcm += xi*density;						// Accumulate x_cm
			sumYcm += yi*density;						// Accumulate y_cm
			sumSigmaMass += density*density; 			// Accumulate mass variance
			sumSigmaXcm += (xi*density) * (xi*density);	// Accumulate x_cm variance
			sumSigmaYcm += (yi*density) * (yi*density);	// Accumulate y_cm variance
		}
	}
	
	delete gROOT->FindObject("c");									// Clean up previous drawings, if any
	delete gROOT->FindObject("g1");
	TCanvas* c = new TCanvas("c", "Rejection", 0, 0, 800, 600);		// Optional plot of the integrant
	TGraph* g1 = new TGraph(min(x.size(), y.size()), &x[0], &y[0]);	// and the actual points used for integration
	gStyle->SetOptStat("mr");
	g1->Draw("AP");
	
	mass = area * sumMass/N;							// Mass
	x_cm = (area * sumXcm/N)/mass; 						// x_cm
	y_cm = (area * sumYcm/N)/mass; 						// y_cm
	errorMass = mass*sqrt(1-x.size()/N)/sqrt(x.size());	// Mass error
	errorXcm = x_cm*sqrt(1-x.size()/N)/sqrt(x.size());	// x_cm error
	errorYcm = y_cm*sqrt(1-x.size()/N)/sqrt(x.size());	// y_cm error
	printf("\nMonte Carlo Rejection Method (%d samples)\n", N);
	printf("mass = %.3f +- %4.3f\n", mass, errorMass);
	printf("x_cm = %.3f +- %4.3f\n", x_cm, errorXcm);
	printf("y_cm = %.3f +- %4.3f\n", y_cm, errorYcm);
	
	// Analytic Values for Comparison	
	const double M = 3.48857;
	const double X = 0.840841;
	const double Y = 0.840841;
	printf("\nAnalytic Evaluation for Comparison\n");
	printf("mass = %.3f\n", M);
	printf("x_cm = %.3f\n", X);
	printf("y_cm = %.3f\n", Y);
}