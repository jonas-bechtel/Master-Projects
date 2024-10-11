#pragma once

#include <TApplication.h>  // For initializing the ROOT application
#include <TH1F.h>          // For creating 1D histograms
#include <TH3D.h>
#include <TCanvas.h>       // For drawing on canvas
#include <TRandom.h>       // For generating random numbers
#include <TSystem.h>
#include <TF1.h>
#include <TGraph.h>
#include <TRandom3.h>

#include "TUnuran.h"
#include "TUnuranContDist.h"
#include "TUnuranMultiContDist.h"
#include "TUnuranDiscrDist.h"
#include "TUnuranEmpDist.h"

#include <unuran.h>
#include <unuran_urng_rngstreams.h>

class ROOTStuff
{
public:
	ROOTStuff();
	void TestStuff(float mean);
	void SampleTest();
	void SampleTest3D();
	void SampleTestMCMC3D(float slider);
	TH3D* GenerateGoalHist();
	void updateCanvas();


private:
	// Initialize the ROOT application
	
	TCanvas* canvas1;
	TCanvas* canvas2;

};

