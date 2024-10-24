#include <chrono>
#include <TAxis3D.h>
#include <TView.h>

#include "MCMC.h"
#include "ElectronBeam.h"
#include "IonBeam.h"
#include "EnergyDistributionModel.h"

MCMC::MCMC()
	: Distribution3D("MCMC")
{
}

void MCMC::SetupDistribution(std::filesystem::path file)
{
	IonBeam* ionBeam = (IonBeam*)Module::Get("Ion Beam");
	targetDist = ionBeam->MultiplyWithElectronDensities();

	targetDist->SetTitle("target density");
	axisRanges[0] = targetDist->GetXaxis()->GetXmin();
	axisRanges[1] = targetDist->GetXaxis()->GetXmax();
	axisRanges[2] = targetDist->GetYaxis()->GetXmin();
	axisRanges[3] = targetDist->GetYaxis()->GetXmax();
	axisRanges[4] = targetDist->GetZaxis()->GetXmin();
	axisRanges[5] = targetDist->GetZaxis()->GetXmax();
}

//void MCMC::SetTargetDistribution(TH3D* targetDist)
//{
//	if (targetDist)
//	{
//		this->targetDist = targetDist;
//		this->targetDist->SetTitle("target density");
//		axisRanges[0] = targetDist->GetXaxis()->GetXmin();
//		axisRanges[1] = targetDist->GetXaxis()->GetXmax();
//		axisRanges[2] = targetDist->GetYaxis()->GetXmin();
//		axisRanges[3] = targetDist->GetYaxis()->GetXmax();
//		axisRanges[4] = targetDist->GetZaxis()->GetXmin();
//		axisRanges[5] = targetDist->GetZaxis()->GetXmax();
//	}
//}

std::vector<Point3D>& MCMC::GetSamples()
{
	return chain;
}

MCMC_Parameters MCMC::GetParameter()
{
	return parameter;
}

void MCMC::SetParameter(MCMC_Parameters params)
{
	parameter = params;
}

void MCMC::ShowUI()
{
	ImGui::SetNextItemWidth(200.0f); ImGui::InputInt("chain length", &parameter.numberSamples);
	ImGui::SetNextItemWidth(200.0f); ImGui::InputInt("burn in", &parameter.burnIn);
	ImGui::SetNextItemWidth(200.0f); ImGui::InputInt("lag", &parameter.lag);

	ImGui::Checkbox("automatically set proposal sigma", &automaticProposalStd);
	ImGui::BeginDisabled(automaticProposalStd);
	ImGui::SetNextItemWidth(200.0f);
	if (ImGui::InputFloat3("Sigma of Proposal functions (x,y,z)", parameter.proposalSigma, "%.4f"))
	{
		normalDistX = std::normal_distribution<double>(0.0, parameter.proposalSigma[0]);
		normalDistY = std::normal_distribution<double>(0.0, parameter.proposalSigma[1]);
		normalDistZ = std::normal_distribution<double>(0.0, parameter.proposalSigma[2]);
	}
	ImGui::EndDisabled();

	ImGui::SetNextItemWidth(200.0f);
	ImGui::InputInt("Seed", &parameter.seed);
	ImGui::SameLine();
	ImGui::Checkbox("change the seed", &changeSeed);
	ImGui::SetNextItemWidth(200.0f);
	if (RebinningFactorInput())
	{
		IonBeam* ionBeam = (IonBeam*)Module::Get("Ion Beam");
		ElectronBeam* eBeam = (ElectronBeam*)Module::Get("Electron Beam");

		ionBeam->PlotDistribution();
		eBeam->PlotDistribution();
		PlotTargetDistribution();
	}

	if (ImGui::Button("generate chain"))
	{
		GenerateSamples();

		PlotDistribution();
		PlotAutocorrelation();
		PlotProjections();
	}
	ImGui::SameLine();
	ImGui::Checkbox("async", &generateAsync);

	ImGui::LabelText("", "Took %.1f ms total. Interpolation took %.1f ms", totalTime, interpolationTime);
	ImGui::LabelText("", "Acceptance Rate: %.1f %%", acceptanceRate * 100);
}

void MCMC::GenerateSamples()
{
	if (!targetDist)
	{
		std::cout << "no target histogram was set\n";
		return;
	}

	// use the std of the target distribution as sigmas for proposal functions
	if (automaticProposalStd)
	{
		parameter.proposalSigma[0] = (float)targetDist->GetStdDev(1);
		parameter.proposalSigma[1] = (float)targetDist->GetStdDev(2);
		parameter.proposalSigma[2] = (float)targetDist->GetStdDev(3);

		normalDistX = std::normal_distribution<double>(0.0, parameter.proposalSigma[0]);
		normalDistY = std::normal_distribution<double>(0.0, parameter.proposalSigma[1]);
		normalDistZ = std::normal_distribution<double>(0.0, parameter.proposalSigma[2]);
	}

	int nXBins = targetDist->GetXaxis()->GetNbins();
	int nYBins = targetDist->GetYaxis()->GetNbins();
	int nZBins = targetDist->GetZaxis()->GetNbins();

	delete m_distribution;
	m_distribution = new TH3D("generated Histogram", "generated Histogram", 
		nXBins, axisRanges[0], axisRanges[1],
		nYBins, axisRanges[2], axisRanges[3],
		nZBins, axisRanges[4], axisRanges[5]);

	m_distribution->SetXTitle("X-axis");
	m_distribution->SetYTitle("Y-axis");
	m_distribution->SetZTitle("Z-axis");

	chain.clear();
	chain.reserve(parameter.numberSamples);
	futures.clear();
	m_distribution->Reset();
	interpolationTime = 0;

	if (changeSeed)
	{
		parameter.seed = (int)std::time(0);
	}
	generator.seed(parameter.seed);

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

		while (chain.size() < parameter.numberSamples)
		{
			totalIterations++;

			if (GenerateSingleSample(currentPoint, currentValue, generator))
			{
				acceptedValues++;
			};

			if (burnInCounter < parameter.burnIn)
			{
				burnInCounter++;
				continue;
			}

			if (lagCounter == 0)
			{
				chain.push_back(currentPoint);
				m_distribution->Fill(currentPoint.x, currentPoint.y, currentPoint.z);
				lagCounter = parameter.lag;
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
		6364136223846793005> generatorLocal(parameter.seed + index);

	std::vector<Point3D>& subchain = subchains[index];
	int subchainLength = parameter.numberSamples / numThreads;
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

		if (GenerateSingleSample(currentPoint, currentValue, generatorLocal))
		{
			acceptedValues++;
		};

		if (burnInCounter < parameter.burnIn)
		{
			burnInCounter++;
			continue;
		}

		if (lagCounter == 0)
		{
			subchain.push_back(currentPoint);
			lagCounter = parameter.lag;
		}
		else
		{
			lagCounter--;
		}
	}
	return (float)acceptedValues / totalIterations;
}

bool MCMC::GenerateSingleSample(Point3D& current, double& currentValue, std::mersenne_twister_engine<std::uint_fast64_t,
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

	int x_nBins = targetDist->GetXaxis()->GetNbins();
	int y_nBins = targetDist->GetYaxis()->GetNbins();
	int z_nBins = targetDist->GetZaxis()->GetNbins();

	double x_modified = std::min(std::max(x_proposed, targetDist->GetXaxis()->GetBinCenter(1)), targetDist->GetXaxis()->GetBinCenter(x_nBins) - 1e-5);
	double y_modified = std::min(std::max(y_proposed, targetDist->GetYaxis()->GetBinCenter(1)), targetDist->GetYaxis()->GetBinCenter(y_nBins) - 1e-5);
	double z_modified = std::min(std::max(z_proposed, targetDist->GetZaxis()->GetBinCenter(1)), targetDist->GetZaxis()->GetBinCenter(z_nBins) - 1e-5);

	// compute probabilities
	auto t_int_start = std::chrono::high_resolution_clock::now();
	double p_new = targetDist->Interpolate(x_modified, y_modified, z_modified);

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
	if (chain.size() != parameter.numberSamples)
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
	means[0] /= parameter.numberSamples;
	means[1] /= parameter.numberSamples;
	means[2] /= parameter.numberSamples;

	// Compute variances
	for (const Point3D& point : chain) 
	{
		variances[0] += (point.x - means[0]) * (point.x - means[0]);
		variances[1] += (point.y - means[1]) * (point.y - means[1]);
		variances[2] += (point.z - means[2]) * (point.z - means[2]);
	}
	variances[0] /= parameter.numberSamples;
	variances[1] /= parameter.numberSamples;
	variances[2] /= parameter.numberSamples;

	std::vector<double> autocorrX(parameter.numberSamples);
	std::vector<double> autocorrY(parameter.numberSamples);
	std::vector<double> autocorrZ(parameter.numberSamples);

	// Compute autocorrelation
	for (int lag = 1; lag < 100; lag++) 
	{
		double sumX = 0;
		double sumY = 0;
		double sumZ = 0;

		for (int i = 0; i < parameter.numberSamples - lag;  i++)
		{
			Point3D point = chain[i];
			Point3D pointLag = chain[i + lag];
			sumX += (point.x - means[0]) * (pointLag.x - means[0]);
			sumY += (point.y - means[1]) * (pointLag.y - means[1]);
			sumZ += (point.z - means[2]) * (pointLag.z - means[2]);
		}
		autocorrX[lag] = sumX / (variances[0] * (parameter.numberSamples - lag));
		autocorrY[lag] = sumY / (variances[1] * (parameter.numberSamples - lag));
		autocorrZ[lag] = sumZ / (variances[2] * (parameter.numberSamples - lag));
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
	if (!m_distribution || !targetDist)
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

std::string MCMC_Parameters::String()
{
	std::string string = std::string(Form("# mcmc sampling parameter:\n")) +
						 std::string(Form("# number of samples: %d\n", numberSamples)) +
						 std::string(Form("# burn in: %d\n", burnIn)) +
						 std::string(Form("# lag: %d\n", lag)) +
						 std::string(Form("# proposal sigmas (x, y, z): %.4f, %.4f, %.3f m\n", proposalSigma[0], proposalSigma[1], proposalSigma[2])) +
						 std::string(Form("# seed: %d\n", seed));

	return string;
}
