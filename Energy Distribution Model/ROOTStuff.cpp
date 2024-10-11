#include "ROOTStuff.h"

#include <string>
#include <iostream>
#include <random>
#include <ctime>
#include <chrono>
#include <algorithm>

ROOTStuff::ROOTStuff()
{
    // Create a canvas where the histogram will be drawn
    canvas1 = new TCanvas("canvas1", "Canvas 1", 800, 600);
    canvas2 = new TCanvas("canvas2", "Canvas 2", 800, 600);
    canvas1->Close();
    canvas2->Close();
}

void ROOTStuff::TestStuff(float mean)
{
    // Create a 1D histogram with 100 bins, range from 0 to 10
    TH1F* h1 = new TH1F("h1", "Random Numbers Histogram;X-axis;Y-axis", 100, 0, 3);

    TUnuran unr;

    TF1* pdfFunc = new TF1("pdf", "1", 0, 1);
    //pdfFunc->SetParameter(0, mean);
    //pdfFunc->SetParameter(1, mean);
    //pdfFunc->SetParameter(2, mean);
    
    pdfFunc->Draw();

    //1D case: create a distribution from two TF1 object pointers pdfFunc
    TUnuranContDist dist(pdfFunc);

    //initialize unuran passing the distribution and a string defining the method
    unr.Init(dist, "method=nrou");
    std::cout << "bla";
    //char buff[100];
    //snprintf(buff, sizeof(buff), "normal(%f,1)", mean);
    //std::string function = buff;
    //
    //if (!unr.Init("normal()", "method=arou")) {
    //    std::cout << "could not make the function sampler\n";
    //    return;
    //}

    // Fill the histogram with 10,000 random numbers following a Gaussian distribution
    for (int i = 0; i < 100; i++) {
        double x = unr.Sample();
        std::cout << x << std::endl;
        h1->Fill(x);  
    }

    // Draw the histogram on the canvas
    h1->Draw();
     
    
    // Cleanup
    //delete h1;
    //delete canvas;
}

void ROOTStuff::SampleTest()
{
    const int N = 1e6;

    const int n = 21;
    double data_x[n], data_y[n];
    
    for (int i = 0; i < n; i++) {
        data_x[i] = i * 0.05;
        data_y[i] = rand() % 10;
    }
    
    TGraph* graph = new TGraph(n, data_x, data_y);
    graph->Draw("AC*");

    TH1D* h1 = new TH1D("h1", "Random Numbers Histogram;X-axis;Y-axis", 30, 0, 1);

    std::default_random_engine generator(std::time(0));
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    int accepted_values = 0;

    //rejection sampling
    auto t_start = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < N; i++) {
        double x = distribution(generator);
        double y = distribution(generator) * *std::max_element(data_y, data_y + n);

        // check to accept or reject
        auto upper = std::upper_bound(data_x, data_x + n, x);
        int index_lower = (upper - data_x) - 1;

        if (upper == data_x + n)
        {
            std::cout << "not found for " << x << "\n";
            continue;
        }
        double x_upper = *upper;
        double x_lower = *(upper - 1);

        double y_upper = data_y[index_lower + 1];
        double y_lower = data_y[index_lower];

        double interpolated_value = (y_lower * (x_upper - x) + y_upper * (x - x_lower)) / (x_upper - x_lower);

        //std::cout << "x: " << x << "bounds: " << x_lower << " " << x_upper << " index of lower: " << (index) << "\n";
        //std::cout << "proposed sample: " << x << " function value: " << interpolated_value << " y: " << y << "\n";
        
        if (y <= interpolated_value)
        {
            // accept
            h1->Fill(x);
            accepted_values++;
        }


        //std::cout << x << std::endl;
        
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    double elapsed_time_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();

    std::cout << "rejection rate: " << (float)(N - accepted_values) / N << " time: " << elapsed_time_ms << "ms\n";
    std::cout << "Entries: " << h1->GetEntries() << " out of " << N << "\n";

    // Draw the histogram on the canvas
    h1->Scale(150.0 / h1->GetEntries());
    h1->Draw("HIST SAME");
}

void ROOTStuff::SampleTest3D()
{
    const int N = 1e5;
    int accepted_values = 0;

    TH3D* h1 = new TH3D("h1", "3D stuff", 10, -3, 3, 10, -3, 3, 10, -3, 3);
    TH3D* h2 = new TH3D("h2", "3D test", 30, -3, 3, 30, -3, 3, 30, -3, 3);

    // Create a random number generator
    TRandom3* rand = new TRandom3();

    // Number of entries to fill
    int nEntries = 1e4;

    // Fill the histogram with random Gaussian values
    for (int i = 0; i < nEntries; i++) {
        double x = rand->Gaus(0, 1); // Mean 0, Sigma 1
        double y = rand->Gaus(0, 1); // Mean 0, Sigma 1
        double z = rand->Gaus(0, 1); // Mean 0, Sigma 1
        h1->Fill(x, y, z);
    }
    //std::cout << h1->GetMaximum() << "\n";

    std::default_random_engine generator(std::time(0));
    std::uniform_real_distribution<double> distributionPoints(-3.0, 3.0);
    std::uniform_real_distribution<double> distributionValue(0.0, 1.0);

    //rejection sampling
    auto t_start = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < N; i++) {
        double x = distributionPoints(generator);
        double y = distributionPoints(generator);
        double z = distributionPoints(generator);
        double value = distributionValue(generator) * h1->GetMaximum();
        
        // check to accept or reject
        //std::cout << x << " " << y << " " << z << "\n";
        double interpolated_value = h1->Interpolate(x, y, z);
        //std::cout << interpolated_value << " drawn value: " << value << "\n";

        if (value <= interpolated_value)
        {
            // accept
            h2->Fill(x, y, z);
            accepted_values++;
            //std::cout << "accepted\n";
        }
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    double elapsed_time_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();

    std::cout << "rejection rate: " << (float)(N - accepted_values) / N << " time: " << elapsed_time_ms << "ms\n";
    std::cout << "Entries: " << h2->GetEntries() << " out of " << N << "\n";

    //h1->Draw("BOX2");
    h2->Draw("SAME BOX2 Z");
}

void ROOTStuff::SampleTestMCMC3D(float slider)
{
    const int N = 1e6;
    const int BurnIn = 1000;
    int BurnInCounter = 0;

    int nXBins = 30;
    int nYBins = 30;
    int nZBins = 30;

    float x_min = -3.0;
    float x_max = 3.0;
    float y_min = -3.0;
    float y_max = 3.0;
    float z_min = -3.0;
    float z_max = 3.0;

    TH3D* h1 = GenerateGoalHist();
    TH3D* h2 = new TH3D("h2", "Recreated Distribution", nXBins, x_min, x_max, nYBins, y_min, y_max, nZBins, z_min, z_max);
    float x_step = (x_max - x_min) / nXBins;

    float step_size = slider;

    std::default_random_engine generator(std::time(0));
    std::uniform_real_distribution<double> distributionStart(-3.0, 3.0);
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    std::normal_distribution<double> normalDistribution(0.0, step_size);

    // random start point
    double x = distributionStart(generator);
    double y = distributionStart(generator);
    double z = distributionStart(generator);

    auto t_start = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < N; i++)
    {
        // propose new sample
        double x_proposed = x + normalDistribution(generator);
        double y_proposed = y + normalDistribution(generator);
        double z_proposed = z + normalDistribution(generator);

        if (x_proposed > x_max - x_step or x_proposed < x_min + x_step)
        {
            //std::cout << "skipped\n";
            h2->Fill(x, y, z);
            continue;
        }

        // compute probabilities
        gErrorIgnoreLevel = 6001;
        double p_current = h1->Interpolate(x, y, z);
        double p_new = h1->Interpolate(x_proposed, y_proposed, z_proposed);

        // acceptance ratio
        double ratio = p_new / p_current;

        if (ratio >= 1 || distribution(generator) < ratio)
        {
            // Accept the new point
            x = x_proposed;
            y = y_proposed;
            z = z_proposed;
        }

        if (BurnInCounter > BurnIn)
        {
            h2->Fill(x, y, z);
        }
        BurnInCounter++;
       
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    double elapsed_time_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();

    std::cout << "time: " << elapsed_time_ms << "ms\n";

    TH1* h1_1D_z = h1->Project3D("z");
    h1_1D_z->SetLineColor(kBlue);
    TH1* h2_1D_z = h2->Project3D("z");
    h2_1D_z->SetLineColor(kRed);

    h1_1D_z->Scale(1 / h1_1D_z->GetMaximum());
    h2_1D_z->Scale(1 / h2_1D_z->GetMaximum());

    canvas1->cd();
    h1_1D_z->Draw("Hist");
    h2_1D_z->Draw("SAME Hist E");

    TH1* h1_1D_y = h1->Project3D("y");
    h1_1D_y->SetLineColor(kBlue);
    TH1* h2_1D_y = h2->Project3D("y");
    h2_1D_y->SetLineColor(kRed);

    h1_1D_y->Scale(1 / h1_1D_y->GetMaximum());
    h2_1D_y->Scale(1 / h2_1D_y->GetMaximum());

    canvas2->cd();
    h1_1D_y->Draw("Hist");
    h2_1D_y->Draw("SAME Hist E");
    //h2->Draw("BOX2");
}

TH3D* ROOTStuff::GenerateGoalHist()
{
    int nXBins = 10;
    int nYBins = 10;
    int nZBins = 10;

    TH3D* h1 = new TH3D("h1", "Goal Distribution", nXBins, -3, 3, nYBins, -3, 3, nZBins, -3, 3);

    // old random way
    //// Create a random number generator
    //TRandom3* rand = new TRandom3();
    //
    //// Number of entries to fill
    //int nEntries = 1e4;
    //
    //// Fill the histogram with random Gaussian values
    //for (int i = 0; i < nEntries; i++) {
    //    double x = rand->Uniform(-3.0, 3.0);
    //    double y = rand->Gaus(0, 0.5); // Mean 0, Sigma 1
    //    double z = rand->Gaus(0, 0.5); // Mean 0, Sigma 1
    //    h1->Fill(x, y, z);
    //}
    double sigma = 1;

    for (int i = 0; i < nXBins; i++) {
        for (int j = 0; j < nYBins; j++) {
            for (int k = 0; k < nZBins; k++) {
                double x = i; //-5 + 10.0 * i / (nx - 1);
                double y = -3 + 6.0 * j / (nYBins - 1);
                double z = -3 + 6.0 * k / (nZBins - 1);
                double value = exp(-(y * y + z * z) / 2.0 / sigma);  // Gaussian example
                h1->SetBinContent(i, j, k, value);
            }
        }
    }
    h1->GetXaxis()->SetTitle("X-Axis");
    h1->GetYaxis()->SetTitle("Y-Axis");
    h1->GetZaxis()->SetTitle("Z-Axis");
    return h1;
}



void ROOTStuff::updateCanvas()
{
    // Embed ROOT Canvas by manually updating it
    canvas1->cd();
    canvas1->Modified();
    canvas1->Update();

    canvas2->cd();
    canvas2->Modified();
    canvas2->Update();
    

   
}


