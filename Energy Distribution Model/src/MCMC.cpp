#include "MCMC.h"
#include "FileLoader.h"

#include <chrono>
#include <TAxis3D.h>
#include <TView.h>
#include <TRootCanvas.h>

MCMC::MCMC()
	: Module("MCMC")
{
}

void MCMC::SetTargetDistribution(TH3D* targetDist)
{
	if (targetDist)
	{
		this->targetDist = targetDist;
		this->targetDist->SetTitle("target density");
		axisRanges[0] = targetDist->GetXaxis()->GetXmin();
		axisRanges[1] = targetDist->GetXaxis()->GetXmax();
		axisRanges[2] = targetDist->GetYaxis()->GetXmin();
		axisRanges[3] = targetDist->GetYaxis()->GetXmax();
		axisRanges[4] = targetDist->GetZaxis()->GetXmin();
		axisRanges[5] = targetDist->GetZaxis()->GetXmax();
		PlotTargetDistribution();
	}
}

//void MCMC::GenerateTargetDistribution()
//{
//	int nXBins = numberBins[0];
//	int nYBins = numberBins[1];
//	int nZBins = numberBins[2];
//
//	delete targetDist;
//	targetDist = new TH3D("Goal Distribution", "Goal Distribution", nXBins, axisRanges[0], axisRanges[1], 
//													 nYBins, axisRanges[2], axisRanges[3],
//													 nZBins, axisRanges[4], axisRanges[5]);
//
//	double sigma = 1;
//	double A = 1;
//	double offset = 0;
//	double omega = 1 / 5.0;
//
//	for (int i = 1; i <= nXBins; i++) {
//		for (int j = 1; j <= nYBins; j++) {
//			for (int k = 1; k <= nZBins; k++) {
//				// Calculate the coordinates for this bin
//				double x = -3 + 6.0 * i / (nXBins - 1);   // Adjust x to range from -5 to 5
//				double y = -3 + 6.0 * j / (nYBins - 1);   // y range from -3 to 3
//				double z = -3 + 6.0 * k / (nZBins - 1);   // z range from -3 to 3
//
//				// Calculate the center of the tube in the y direction as a function of x
//				double y_center = A * sin(offset + omega * x * 3.14);  // Sinusoidal bend in the y direction
//
//				// Calculate the value using the Gaussian distribution centered at (y_center, z = 0)
//				double value = exp(-((y - y_center) * (y - y_center) + z * z) / (2.0 * sigma * sigma));
//
//				//double value = exp(-(y * y + z * z) / 2.0 / sigma);  // Gaussian example
//				targetDist->SetBinContent(i, j, k, value);
//			}
//		}
//	}
//	targetDist->GetXaxis()->SetTitle("X-Axis");
//	targetDist->GetYaxis()->SetTitle("Y-Axis");
//	targetDist->GetZaxis()->SetTitle("Z-Axis");
//}

std::vector<Point3D>& MCMC::GetSamples()
{
	return chain;
}

void MCMC::ShowUI()
{
	ImGui::InputInt("chain length", &chainLength);
	ImGui::InputInt("burn in", &burnIn);
	ImGui::InputInt("lag", &lag);
	if (ImGui::InputFloat3("Sigma of Proposal functions (x,y,z)", proposalSigma))
	{
		normalDistX = std::normal_distribution<double>(0.0, proposalSigma[0]);
		normalDistY = std::normal_distribution<double>(0.0, proposalSigma[1]);
		normalDistZ = std::normal_distribution<double>(0.0, proposalSigma[2]);
	}

	ImGui::Checkbox("change the seed", &changeSeed);
	ImGui::SameLine();
	ImGui::InputInt("Seed", &seed);
	RebinningFactorInput();

	if (ImGui::Button("generate chain"))
	{
		GenerateChain();

		PlotDistribution();
		PlotAutocorrelation();
		PlotProjections();
	}
	ImGui::SameLine();
	ImGui::Checkbox("async", &generateAsync);

	ImGui::LabelText("", "Took %.1f ms total. Interpolation took %.1f ms", totalTime, interpolationTime);
	ImGui::LabelText("", "Acceptance Rate: %.1f %%", acceptanceRate * 100);

}

void MCMC::GenerateChain()
{
	if (!targetDist)
	{
		std::cout << "no target histogram was set\n";
		return;
	}

	double nXBins = targetDist->GetXaxis()->GetNbins();
	double nYBins = targetDist->GetYaxis()->GetNbins();
	double nZBins = targetDist->GetZaxis()->GetNbins();

	delete m_distribution;
	m_distribution = new TH3D("generated Histogram", "generated Histogram", 
		nXBins, axisRanges[0], axisRanges[1],
		nYBins, axisRanges[2], axisRanges[3],
		nZBins, axisRanges[4], axisRanges[5]);

	m_distribution->SetXTitle("X-axis");
	m_distribution->SetYTitle("Y-axis");
	m_distribution->SetZTitle("Z-axis");

	chain.clear();
	chain.reserve(chainLength);
	futures.clear();
	m_distribution->Reset();
	interpolationTime = 0;

	if (changeSeed)
	{
		seed = std::time(0);
	}
	generator.seed(seed);

    auto t_start = std::chrono::high_resolution_clock::now();

	if (generateAsync)
	{
		// Launch asynchronous tasks
		for (int i = 0; i < numThreads; i++)
		{
			futures.push_back(std::async(std::launch::async, &MCMC::GenerateSubchain, this, i));
		}

		acceptanceRate = 0;

		// Wait for all tasks to complete
		for (auto& future : futures) {
			acceptanceRate += future.get();
		}

		// combine subchains
		for (int i = 0; i < numThreads; i++)
		{
			chain.insert(chain.end(), subchains[i].begin(), subchains[i].end());
		}

		acceptanceRate /= numThreads;

		// Fill graphs
		for (Point3D& point : chain)
		{
			m_distribution->Fill(point.x, point.y, point.z);
		}
	}
	else
	{
		// "random" start point
		double x = 0;	//uniformDist(generator);
		double y = 0;	//uniformDist(generator);
		double z = 0.3;	//uniformDist(generator);

		Point3D currentPoint(x, y, z);
		double currentValue = targetDist->Interpolate(x, y, z);

		int burnInCounter = 0;
		int lagCounter = 0;
		int acceptedValues = 0;
		int totalIterations = 0;

		while (chain.size() < chainLength)
		{
			totalIterations++;

			if (GenerateSample(currentPoint, currentValue, generator))
			{
				acceptedValues++;
			};

			if (burnInCounter < burnIn)
			{
				burnInCounter++;
				continue;
			}

			if (lagCounter == 0)
			{
				chain.push_back(currentPoint);
				m_distribution->Fill(currentPoint.x, currentPoint.y, currentPoint.z);
				lagCounter = lag;
			}
			else
			{
				lagCounter--;
			}

		}
		acceptanceRate = (float)acceptedValues / totalIterations;
	}

    auto t_end = std::chrono::high_resolution_clock::now();
    totalTime = std::chrono::duration<double, std::milli>(t_end - t_start).count();
}

float MCMC::GenerateSubchain(int index)
{
	thread_local std::mersenne_twister_engine<std::uint_fast64_t,
		64, 312, 156, 31,
		0xb5026f5aa96619e9, 29,
		0x5555555555555555, 17,
		0x71d67fffeda60000, 37,
		0xfff7eee000000000, 43,
		6364136223846793005> generatorLocal(seed + index);

	std::vector<Point3D>& subchain = subchains[index];
	int subchainLength = chainLength / numThreads;
	subchain.clear();
	subchain.reserve(subchainLength);

	// random start point
	double x = 0; //uniformDist(generator);
	double y = 0; //uniformDist(generator);
	double z = 0.3; //uniformDist(generator);

	Point3D currentPoint(x, y, z);
	double currentValue = targetDist->Interpolate(x, y, z);

	int burnInCounter = 0;
	int lagCounter = 0;
	int acceptedValues = 0;
	int totalIterations = 0;

	while (subchain.size() < subchainLength)
	{
		totalIterations++;

		if (GenerateSample(currentPoint, currentValue, generatorLocal))
		{
			acceptedValues++;
		};

		if (burnInCounter < burnIn)
		{
			burnInCounter++;
			continue;
		}

		if (lagCounter == 0)
		{
			subchain.push_back(currentPoint);
			lagCounter = lag;
		}
		else
		{
			lagCounter--;
		}
	}
	return (float)acceptedValues / totalIterations;
}

bool MCMC::GenerateSample(Point3D& current, double& currentValue, std::mersenne_twister_engine<std::uint_fast64_t,
	64, 312, 156, 31,
	0xb5026f5aa96619e9, 29,
	0x5555555555555555, 17,
	0x71d67fffeda60000, 37,
	0xfff7eee000000000, 43,
	6364136223846793005>& generator)
{
	// propose new sample
	double x_proposed = current.x + normalDistX(generator); //axisRanges[0] + uniformDist(generator) * (axisRanges[1] - axisRanges[0]); 
	double y_proposed = current.y + normalDistY(generator);
	double z_proposed = current.z + normalDistZ(generator);

	// check if point is outside histogram domain
	if (x_proposed < axisRanges[0] || x_proposed > axisRanges[1] ||
		y_proposed < axisRanges[2] || y_proposed > axisRanges[3] ||
		z_proposed < axisRanges[4] || z_proposed > axisRanges[5])
	{
		return false;
	}

	double x_nBins = targetDist->GetXaxis()->GetNbins();
	double y_nBins = targetDist->GetYaxis()->GetNbins();
	double z_nBins = targetDist->GetZaxis()->GetNbins();

	double x_modified = std::min(std::max(x_proposed, targetDist->GetXaxis()->GetBinCenter(1)), targetDist->GetXaxis()->GetBinCenter(x_nBins) - 1e-5);
	double y_modified = std::min(std::max(y_proposed, targetDist->GetYaxis()->GetBinCenter(1)), targetDist->GetYaxis()->GetBinCenter(y_nBins) - 1e-5);
	double z_modified = std::min(std::max(z_proposed, targetDist->GetZaxis()->GetBinCenter(1)), targetDist->GetZaxis()->GetBinCenter(z_nBins) - 1e-5);

	// compute probabilities
	auto t_int_start = std::chrono::high_resolution_clock::now();
	double p_new = targetDist->Interpolate(x_modified, y_modified, z_modified);

	//std::cout << "x: " << x_proposed << "y: " << y_proposed << "z: " << z_proposed << "\n";
	//std::cout << "x modified: " << x_proposed << std::min(std::max(x_proposed, targetDist->GetXaxis()->GetBinCenter(1)), targetDist->GetXaxis()->GetBinCenter(100) - 1e-5);
	auto t_int_end = std::chrono::high_resolution_clock::now();
	interpolationTime += std::chrono::duration<double, std::milli>(t_int_end - t_int_start).count();

	// acceptance ratio
	double ratio = p_new / currentValue;

	if (ratio >= 1 || uniformDist(generator) < ratio)
	{
		// Accept the new point
		current.x = x_proposed;
		current.y = y_proposed;
		current.z = z_proposed;
		currentValue = p_new;

		return true;
	}

	return false;
}

void MCMC::PlotTargetDistribution()
{
	if (!targetDist) return;
	if (targetDistSmall) delete targetDistSmall;

	m_mainCanvas->cd(1);
	targetDistSmall = (TH3D*)targetDist->Rebin3D(s_rebinningFactors[0],
												 s_rebinningFactors[1],
												 s_rebinningFactors[2], "target Distribution small");
	targetDistSmall->Draw("BOX2");
}

void MCMC::PlotAutocorrelation()
{
	if (chain.size() != chainLength)
	{
		std::cout << "the chain has the wrong size: " << chain.size() << "\n";
		return;
	}

	double means[3] = { 0, 0, 0 };
	double variances[3] = { 0, 0, 0 };

	// Compute means and variances
	for (const Point3D& point : chain) 
	{
		means[0] += point.x;
		means[1] += point.y;
		means[2] += point.z;
	}
	means[0] /= chainLength;
	means[1] /= chainLength;
	means[2] /= chainLength;

	// Compute variances
	for (const Point3D& point : chain) 
	{
		variances[0] += (point.x - means[0]) * (point.x - means[0]);
		variances[1] += (point.y - means[1]) * (point.y - means[1]);
		variances[2] += (point.z - means[2]) * (point.z - means[2]);
	}
	variances[0] /= chainLength;
	variances[1] /= chainLength;
	variances[2] /= chainLength;

	std::vector<double> autocorrX(chainLength);
	std::vector<double> autocorrY(chainLength);
	std::vector<double> autocorrZ(chainLength);

	// Compute autocorrelation
	for (int lag = 1; lag < 100; lag++) 
	{
		double sumX = 0;
		double sumY = 0;
		double sumZ = 0;

		for (int i = 0; i < chainLength - lag;  i++) 
		{
			Point3D point = chain[i];
			Point3D pointLag = chain[i + lag];
			sumX += (point.x - means[0]) * (pointLag.x - means[0]);
			sumY += (point.y - means[1]) * (pointLag.y - means[1]);
			sumZ += (point.z - means[2]) * (pointLag.z - means[2]);
		}
		autocorrX[lag] = sumX / (variances[0] * (chainLength - lag));
		autocorrY[lag] = sumY / (variances[1] * (chainLength - lag));
		autocorrZ[lag] = sumZ / (variances[2] * (chainLength - lag));
	}
	//std::cout << "x variance: " << variances[0] << " mean: " << means[0] << "\n";
	//std::cout << "y variance: " << variances[1] << " mean: " << means[1] << "\n";
	//std::cout << "z variance: " << variances[2] << " mean : " << means[2] << "\n";
		 
	m_secondCanvas->cd(4);
	TGraph* graphX = new TGraph(100, &autocorrX[0]);
	graphX->SetTitle("Autocorrelation X;Lag;Autocorrelation");
	graphX->Draw("AL");

	m_secondCanvas->cd(5);
	TGraph* graphY = new TGraph(100, &autocorrY[0]);
	graphY->SetTitle("Autocorrelation Y;Lag;Autocorrelation");
	graphY->Draw("AL");

	m_secondCanvas->cd(6);
	TGraph* graphZ = new TGraph(100, &autocorrZ[0]);
	graphZ->SetTitle("Autocorrelation Z;Lag;Autocorrelation");
	graphZ->Draw("AL");

}

void MCMC::PlotProjections()
{
	if (!(m_distribution || targetDist))
	{
		return;
	}
	delete projectionXTarget;
	delete projectionYTarget;
	delete projectionZTarget;
	delete projectionX;
	delete projectionY;
	delete projectionZ;

	m_secondCanvas->cd(1);
	projectionXTarget = targetDist->ProjectionX();
	projectionXTarget->Scale(1.0 / projectionXTarget->Integral());
	projectionX = m_distribution->ProjectionX();
	projectionX->SetLineColor(kRed);
	projectionX->Scale(1.0 / projectionX->Integral());
	projectionX->Draw("HIST E");
	projectionXTarget->Draw("HIST SAME");

	m_secondCanvas->cd(2);
	projectionYTarget = targetDist->ProjectionY();
	projectionYTarget->Scale(1.0 / projectionYTarget->Integral());
	projectionY = m_distribution->ProjectionY();
	projectionY->SetLineColor(kRed);
	projectionY->Scale(1.0 / projectionY->Integral());
	projectionY->Draw("HIST E");
	projectionYTarget->Draw("HIST SAME");

	m_secondCanvas->cd(3);
	projectionZTarget = targetDist->ProjectionZ();
	projectionZTarget->Scale(1.0 / projectionZTarget->Integral());
	projectionZ = m_distribution->ProjectionZ();
	projectionZ->SetLineColor(kRed);
	projectionZ->Scale(1.0 / projectionZ->Integral());
	projectionZ->Draw("HIST E");
	projectionZTarget->Draw("HIST SAME");
}


