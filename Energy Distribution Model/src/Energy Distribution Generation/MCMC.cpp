#include "pch.h"

#include "MCMC.h"
#include "ElectronBeam.h"
#include "IonBeam.h"
#include "EnergyDistributionManager.h"

MCMC_Window::MCMC_Window()
	: EnergyDistributionModule("MCMC", 1), m_parameters(activeDist.mcmcParameter)
{
	mcmc = this;

	if (generatorList.size() < numThreads)
	{
		std::cout << "list of generators is not large enough: size = " << generatorList.size() << " # threads = " << numThreads << std::endl;
	}
	for (unsigned int i = 0; i < numThreads; i++)
	{
		generatorList.at(i) = RNG_engine();
	}
}

void MCMC_Window::SetupDistribution(std::filesystem::path file)
{
	delete targetDist;
	targetDist = ionBeam->MultiplyWithElectronDensities();

	if (!targetDist)
	{
		std::cout << "could not multiply ion and electron beam" << std::endl;
		return;
	}	
}

void MCMC_Window::SetTargetDist(TH3D* target)
{
	delete targetDist;
	targetDist = target;

	if (!targetDist)
	{
		std::cout << "target dist was nullptr" << std::endl;
		return;
	}
}

std::vector<Point3D>& MCMC_Window::GetSamples()
{
	return chain;
}

void MCMC_Window::SetParameter(MCMC_Parameters params)
{
	m_parameters = params;
}

void MCMC_Window::ShowUI()
{
	if (ImGui::BeginChild("left side", ImVec2(300, -1), ImGuiChildFlags_ResizeX))
	{
		ImGui::PushItemWidth(90.0f);
		ShowList();
		ShowSettings();
		ImGui::Separator();
		ShowCanvasButtons();
		ImGui::PopItemWidth();
		ImGui::EndChild();
	}

	ImGui::SameLine();

	ShowMCMCDataPlots();

	ShowAutoCorrelationPlots();
}

void MCMC_Window::ShowList()
{
	if (ImGui::BeginListBox("##mcmc data list", ImVec2(-1, 270.0f)))
	{
		for (int i = 0; i < mcmcDataToLookAt.size(); i++)
		{
			ImGui::PushID(i);
			MCMC_Data& mcmcData = mcmcDataToLookAt.at(i);

			if (ImGui::Selectable(mcmcData.label.c_str(), i == selectedIndex, ImGuiSelectableFlags_AllowItemOverlap))
			{
				selectedIndex = i;
				SelectedItemChanged();
			}

			ImGui::SameLine();
			if (ImGui::SmallButton("x"))
			{
				RemoveMCMCDataFromList(i);
			}
			ImGui::PopID();
		}
		ImGui::EndListBox();
	}
}

void MCMC_Window::ShowSettings()
{
	ImGui::SetNextItemWidth(200.0f);
	if (ImGui::InputInt("chain length", m_parameters.numberSamples))
	{
		m_parameters.numberSamples.get() = ((m_parameters.numberSamples.get() + numThreads - 1) / numThreads) * numThreads;

	}
	ImGui::SetNextItemWidth(200.0f); ImGui::InputInt("burn in", m_parameters.burnIn);
	ImGui::SetNextItemWidth(200.0f); ImGui::InputInt("lag", m_parameters.lag);

	ImGui::Checkbox("automatically set proposal sigma", &automaticProposalStd);
	ImGui::BeginDisabled(automaticProposalStd);
	ImGui::SetNextItemWidth(200.0f);
	if (ImGui::InputFloat3("Sigma of Proposal functions (x,y,z)", m_parameters.proposalSigma, "%.4f"))
	{
		normalDistX = std::normal_distribution<double>(0.0, m_parameters.proposalSigma.get().x);
		normalDistY = std::normal_distribution<double>(0.0, m_parameters.proposalSigma.get().y);
		normalDistZ = std::normal_distribution<double>(0.0, m_parameters.proposalSigma.get().z);
	}
	ImGui::EndDisabled();
	ImGui::SetNextItemWidth(200.0f);
	ImGui::BeginDisabled(changeSeed);
	ImGui::InputInt("Seed", m_parameters.seed);
	ImGui::EndDisabled();
	ImGui::SameLine();
	ImGui::Checkbox("change the seed", &changeSeed);
	ImGui::SetNextItemWidth(200.0f);
	if (RebinningFactorInput())
	{
		ionBeam->PlotDistribution();
		eBeam->PlotDistribution();
		PlotTargetDistribution();
	}

	if (ImGui::Button("generate chain"))
	{
		SetupDistribution();
		GenerateSamples();

		if (m_distribution && targetDist)
		{
			MCMC_Data mcmcData;

			mcmcData.targetDist = (TH3D*)targetDist->Clone("target to look at");
			mcmcData.generatedDist = (TH3D*)m_distribution->Clone("generated to look at");

			mcmcData.targetDist->SetTitle("target to look at");
			mcmcData.generatedDist->SetTitle("generated to look at");

			mcmcData.FillData();
			mcmcData.label = eBeam->GetSelected()->label;

			AddMCMCDataToList(mcmcData);
		}
		UpdateAutocorrelationData();
	}
	ImGui::SameLine();
	ImGui::Checkbox("async", &generateAsync);
	ImGui::SameLine();
	ImGui::Checkbox("interpolate", &useInterpolation);

	ImGui::SetNextItemWidth(ImGui::GetContentRegionAvail().x);
	ImGui::LabelText("", "Took %.1f ms total. Interpolation took %.1f ms", totalTime, interpolationTime);
	ImGui::SetNextItemWidth(ImGui::GetContentRegionAvail().x);
	ImGui::LabelText("", "Acceptance Rate: %.1f %%", acceptanceRate * 100);

	ImGui::Checkbox("show autocorrelation", &showAutoCorrelationWindow);

	ImGui::SetNextItemWidth(200.0f);
	if (ImGui::SliderFloat("slice z", &SliceZ, -0.7f, 0.7f))
	{
		if (selectedIndex >= 0)
		{
			MCMC_Data& mcmc = mcmcDataToLookAt.at(selectedIndex);
			mcmc.targetSlice.FromTH3D(mcmc.targetDist, SliceZ);
			mcmc.generatedSlice.FromTH3D(mcmc.generatedDist, SliceZ);
		}
	}
}

void MCMC_Window::ShowMCMCDataPlots()
{
	if (ImPlot::BeginSubplots("##mcmc subplots", 2, 3, ImVec2(-1, -1), ImPlotSubplotFlags_ShareItems))
	{
		if (ImPlot::BeginPlot("Projection X"))
		{
			int i = 0;
			for (const MCMC_Data& mcmc_data : mcmcDataToLookAt)
			{
				std::string label = mcmc_data.label + " " + std::to_string(i++);
				ImPlot::PlotLine(label.c_str(), mcmc_data.xAxis.data(), mcmc_data.targetProjectionValuesX.data(), mcmc_data.xAxis.size(), ImPlotLineFlags_Segments);
				ImPlot::PlotLine(label.c_str(), mcmc_data.xAxis.data(), mcmc_data.generatedProjectionValuesX.data(), mcmc_data.xAxis.size());
			}
			ImPlot::EndPlot();
		}

		if (ImPlot::BeginPlot("Projection Y"))
		{
			int i = 0;
			for (const MCMC_Data& mcmc_data : mcmcDataToLookAt)
			{
				std::string label = mcmc_data.label + " " + std::to_string(i++);
				ImPlot::PlotLine(label.c_str(), mcmc_data.yAxis.data(), mcmc_data.targetProjectionValuesY.data(), mcmc_data.yAxis.size(), ImPlotLineFlags_Segments);
				ImPlot::PlotLine(label.c_str(), mcmc_data.yAxis.data(), mcmc_data.generatedProjectionValuesY.data(), mcmc_data.yAxis.size());
			}
			ImPlot::EndPlot();
		}

		if (ImPlot::BeginPlot("Projection Z"))
		{
			int i = 0;
			for (const MCMC_Data& mcmc_data : mcmcDataToLookAt)
			{
				std::string label = mcmc_data.label + " " + std::to_string(i++);
				ImPlot::PlotLine(label.c_str(), mcmc_data.zAxis.data(), mcmc_data.targetProjectionValuesZ.data(), mcmc_data.zAxis.size(), ImPlotLineFlags_Segments);
				ImPlot::PlotLine(label.c_str(), mcmc_data.zAxis.data(), mcmc_data.generatedProjectionValuesZ.data(), mcmc_data.zAxis.size());
			}
			ImPlot::EndPlot();
		}

		ImPlot::PushColormap(9);
		if (selectedIndex >= 0)
		{
			const MCMC_Data& mcmc_data = mcmcDataToLookAt.at(selectedIndex);
			const HeatMapData& targetSlice = mcmcDataToLookAt.at(selectedIndex).targetSlice;
			const HeatMapData& generatedSlice = mcmcDataToLookAt.at(selectedIndex).generatedSlice;
			std::string label = mcmc_data.label + " " + std::to_string(selectedIndex);

			if (ImPlot::BeginPlot("XY target Slice"))
			{
				ImPlot::SetupAxes("x", "y");
				ImPlot::PlotHeatmap(label.c_str(), targetSlice.values.data(), targetSlice.nRows,
					targetSlice.nCols, targetSlice.minValue, targetSlice.maxValue, nullptr,
					targetSlice.bottomLeft, targetSlice.topRight);
				ImPlot::EndPlot();
			}
			ImGui::SameLine();
			ImPlot::ColormapScale("##targetScale", targetSlice.minValue, targetSlice.maxValue, ImVec2(80, -1));

			ImGui::SameLine();
			if (ImPlot::BeginPlot("XY generated Slice"))
			{
				ImPlot::SetupAxes("x", "y");
				ImPlot::PlotHeatmap(label.c_str(), generatedSlice.values.data(), generatedSlice.nRows,
					generatedSlice.nCols, generatedSlice.minValue, generatedSlice.maxValue, nullptr,
					generatedSlice.bottomLeft, generatedSlice.topRight);
				ImPlot::EndPlot();
			}
			ImGui::SameLine();
			ImPlot::ColormapScale("##generatedScale", generatedSlice.minValue, generatedSlice.maxValue, ImVec2(80, -1));
		}
		ImPlot::PopColormap();

		ImPlot::EndSubplots();
	}
}

void MCMC_Window::ShowAutoCorrelationPlots()
{
	if (showAutoCorrelationWindow && ImGui::Begin("Autocorrelation window", &showAutoCorrelationWindow, ImGuiWindowFlags_NoDocking))
	{
		if (ImPlot::BeginSubplots("##autocorr subplots", 1, 3, ImVec2(-1, -1), ImPlotSubplotFlags_ShareItems))
		{
			if (ImPlot::BeginPlot("Autocorrelation X"))
			{
				ImPlot::PlotLine("##", lagValues.data(), autocorrX.data(), autocorrX.size());
				ImPlot::EndPlot();
			}

			if (ImPlot::BeginPlot("Autocorrelation Y"))
			{
				ImPlot::PlotLine("##", lagValues.data(), autocorrY.data(), autocorrY.size());
				ImPlot::EndPlot();
			}

			if (ImPlot::BeginPlot("Autocorrelation Z"))
			{
				ImPlot::PlotLine("##", lagValues.data(), autocorrZ.data(), autocorrZ.size());
				ImPlot::EndPlot();
			}
			ImPlot::EndSubplots();
		}

		ImGui::End();
	}
}

void MCMC_Window::GenerateSamples()
{
	if (!targetDist)
	{
		std::cout << "no target histogram was set\n";
		return;
	}

	// use the std of the target distribution as sigmas for proposal functions
	if (automaticProposalStd)
	{
		m_parameters.proposalSigma.get().x = (float)targetDist->GetStdDev(1);
		m_parameters.proposalSigma.get().y = (float)targetDist->GetStdDev(2);
		m_parameters.proposalSigma.get().z = (float)targetDist->GetStdDev(3);

		normalDistX = std::normal_distribution<double>(0.0, m_parameters.proposalSigma.get().x);
		normalDistY = std::normal_distribution<double>(0.0, m_parameters.proposalSigma.get().y);
		normalDistZ = std::normal_distribution<double>(0.0, m_parameters.proposalSigma.get().z);
	}

	delete m_distribution;
	m_distribution = (TH3D*)targetDist->Clone("generated Histogram");
	m_distribution->Reset();
	m_distribution->SetTitle("generated Histogram");

	chain.clear();
	chain.resize(m_parameters.numberSamples);
	futures.clear();
	futures.reserve(numThreads);
	m_distribution->Reset();
	interpolationTime = 0;

	if (changeSeed)
	{
		m_parameters.seed.set((int)std::time(0));
	}
	generatorList[0].seed(m_parameters.seed);

    auto t_start = std::chrono::high_resolution_clock::now();

	if (generateAsync)
	{
		// Launch asynchronous tasks
		for (unsigned int i = 0; i < numThreads; i++)
		{
			int length = m_parameters.numberSamples / numThreads;
			int offset = i * length;
			RNG_engine& generator = generatorList.at(i);
			generator.seed(m_parameters.seed.get() + i);
			
			futures.push_back(std::async(std::launch::async, &MCMC_Window::GenerateSubchain, this, length, offset, std::ref(generator)));
		}

		acceptanceRate = 0;

		// Wait for all tasks to complete
		for (auto& future : futures) 
		{
			acceptanceRate += future.get();
		}

		acceptanceRate /= numThreads;
	}
	else
	{
		acceptanceRate = GenerateSubchain(m_parameters.numberSamples, 0, generatorList[0]);
	}

	// Fill graphs
	for (const Point3D& point : chain)
	{
		m_distribution->Fill(point.x, point.y, point.z);
	}

    auto t_end = std::chrono::high_resolution_clock::now();
    totalTime = std::chrono::duration<double, std::milli>(t_end - t_start).count();
}

float MCMC_Window::GenerateSubchain(int length, int offset, RNG_engine& generator)
{
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
	int addedValues = 0;

	while (addedValues < length)
	{
		totalIterations++;

		if (GenerateSingleSample(currentPoint, currentValue, generator))
		{
			acceptedValues++;
		};

		if (burnInCounter < m_parameters.burnIn)
		{
			burnInCounter++;
			continue;
		}

		if (lagCounter == 0)
		{
			chain.at(offset + addedValues) = currentPoint;
			addedValues++;
			lagCounter = m_parameters.lag;
		}
		else
		{
			lagCounter--;
		}
	}
	return (float)acceptedValues / totalIterations;
}

bool MCMC_Window::GenerateSingleSample(Point3D& current, double& currentValue, RNG_engine& generator)
{
	// propose new sample
	double x_proposed = current.x + normalDistX(generator);
	double y_proposed = current.y + normalDistY(generator);
	double z_proposed = current.z + normalDistZ(generator);
	//double z_proposed = uniformDist(generator) * 0.7 ; 

	// check if point is outside histogram domain
	if (x_proposed < targetDist->GetXaxis()->GetXmin() || x_proposed > targetDist->GetXaxis()->GetXmax() ||
		y_proposed < targetDist->GetYaxis()->GetXmin() || y_proposed > targetDist->GetYaxis()->GetXmax() ||
		z_proposed < targetDist->GetZaxis()->GetXmin() || z_proposed > targetDist->GetZaxis()->GetXmax())
	{
		return false;
	}

	int x_nBins = targetDist->GetXaxis()->GetNbins();
	int y_nBins = targetDist->GetYaxis()->GetNbins();
	int z_nBins = targetDist->GetZaxis()->GetNbins();

	// compute probabilities
	auto t_int_start = std::chrono::high_resolution_clock::now();
	double p_new;
	if (useInterpolation)
	{
		double x_modified = std::min(std::max(x_proposed, targetDist->GetXaxis()->GetBinCenter(1)), targetDist->GetXaxis()->GetBinCenter(x_nBins) - 1e-5);
		double y_modified = std::min(std::max(y_proposed, targetDist->GetYaxis()->GetBinCenter(1)), targetDist->GetYaxis()->GetBinCenter(y_nBins) - 1e-5);
		double z_modified = std::min(std::max(z_proposed, targetDist->GetZaxis()->GetBinCenter(1)), targetDist->GetZaxis()->GetBinCenter(z_nBins) - 1e-5);

		p_new = targetDist->Interpolate(x_modified, y_modified, z_modified);
	}
	else
	{
		p_new = targetDist->GetBinContent(targetDist->FindBin(x_proposed, y_proposed, z_proposed));
	}

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

void MCMC_Window::SelectedItemChanged()
{
	MCMC_Data& newlySelected = mcmcDataToLookAt.at(selectedIndex);
	newlySelected.targetSlice.FromTH3D(newlySelected.targetDist, SliceZ);
	newlySelected.generatedSlice.FromTH3D(newlySelected.generatedDist, SliceZ);

	delete targetDist;
	delete m_distribution;
	targetDist = (TH3D*)newlySelected.targetDist->Clone("copied target");
	m_distribution = (TH3D*)newlySelected.generatedDist->Clone("copied generated");
	if (IsCanvasShown(m_mainCanvas))
	{
		PlotTargetDistribution();
		PlotDistribution();
	}
}

void MCMC_Window::AddMCMCDataToList(MCMC_Data& mcmcData)
{
	mcmcDataToLookAt.push_back(std::move(mcmcData));
	if (mcmcDataToLookAt.size() == 1)
	{
		selectedIndex = 0;
		SelectedItemChanged();
	}
}

void MCMC_Window::RemoveMCMCDataFromList(int index)
{
	mcmcDataToLookAt.erase(mcmcDataToLookAt.begin() + index);
	selectedIndex = std::min(selectedIndex, (int)mcmcDataToLookAt.size() - 1);
	if (selectedIndex >= 0)
	{
		SelectedItemChanged();
	}
}

void MCMC_Window::UpdateAutocorrelationData()
{
	// Compute means and variances
	for (const Point3D& point : chain)
	{
		means[0] += point.x;
		means[1] += point.y;
		means[2] += point.z;
	}
	means[0] /= m_parameters.numberSamples;
	means[1] /= m_parameters.numberSamples;
	means[2] /= m_parameters.numberSamples;

	// Compute variances
	for (const Point3D& point : chain)
	{
		variances[0] += (point.x - means[0]) * (point.x - means[0]);
		variances[1] += (point.y - means[1]) * (point.y - means[1]);
		variances[2] += (point.z - means[2]) * (point.z - means[2]);
	}
	variances[0] /= m_parameters.numberSamples;
	variances[1] /= m_parameters.numberSamples;
	variances[2] /= m_parameters.numberSamples;

	// Compute autocorrelation
	for (int lag = 0; lag < 100; lag++)
	{
		double sumX = 0;
		double sumY = 0;
		double sumZ = 0;

		for (int i = 0; i < m_parameters.numberSamples - lag; i++)
		{
			Point3D point = chain[i];
			Point3D pointLag = chain[i + lag];
			sumX += (point.x - means[0]) * (pointLag.x - means[0]);
			sumY += (point.y - means[1]) * (pointLag.y - means[1]);
			sumZ += (point.z - means[2]) * (pointLag.z - means[2]);
		}
		autocorrX[lag] = sumX / (variances[0] * (m_parameters.numberSamples - lag));
		autocorrY[lag] = sumY / (variances[1] * (m_parameters.numberSamples - lag));
		autocorrZ[lag] = sumZ / (variances[2] * (m_parameters.numberSamples - lag));

		lagValues[lag] = lag;
	}
}

void MCMC_Window::PlotTargetDistribution()
{
	if (!targetDist) return;

	m_mainCanvas->cd(1);
	delete targetDistSmall;
	targetDistSmall = (TH3D*)targetDist->Rebin3D(s_rebinningFactors[0],
												 s_rebinningFactors[1],
												 s_rebinningFactors[2], "target Distribution small");
	targetDistSmall->Draw("BOX2");
}

void MCMC_Data::FillData()
{
	xAxis.reserve(targetDist->GetNbinsX());
	yAxis.reserve(targetDist->GetNbinsY());
	zAxis.reserve(targetDist->GetNbinsZ());

	targetProjectionValuesX.reserve(targetDist->GetNbinsX());
	targetProjectionValuesY.reserve(targetDist->GetNbinsY());
	targetProjectionValuesZ.reserve(targetDist->GetNbinsZ());

	generatedProjectionValuesX.reserve(targetDist->GetNbinsX());
	generatedProjectionValuesY.reserve(targetDist->GetNbinsY());
	generatedProjectionValuesZ.reserve(targetDist->GetNbinsZ());

	for (int i = 1; i <= targetDist->GetNbinsX(); i++)
	{
		xAxis.push_back(targetDist->GetXaxis()->GetBinCenter(i));
	}
	for (int i = 1; i <= targetDist->GetNbinsY(); i++)
	{
		yAxis.push_back(targetDist->GetYaxis()->GetBinCenter(i));
	}
	for (int i = 1; i <= targetDist->GetNbinsZ(); i++)
	{
		zAxis.push_back(targetDist->GetZaxis()->GetBinCenter(i));
	}

	TH1D* targetProjectionX = targetDist->ProjectionX();
	TH1D* targetProjectionY = targetDist->ProjectionY();
	TH1D* targetProjectionZ = targetDist->ProjectionZ();
	targetProjectionX->Scale(1.0 / targetProjectionX->Integral());
	targetProjectionY->Scale(1.0 / targetProjectionY->Integral());
	targetProjectionZ->Scale(1.0 / targetProjectionZ->Integral());

	for (int i = 1; i <= targetProjectionX->GetNbinsX(); i++)
	{
		targetProjectionValuesX.push_back(targetProjectionX->GetBinContent(i));
	}
	for (int i = 1; i <= targetProjectionY->GetNbinsX(); i++)
	{
		targetProjectionValuesY.push_back(targetProjectionY->GetBinContent(i));
	}
	for (int i = 1; i <= targetProjectionZ->GetNbinsX(); i++)
	{
		targetProjectionValuesZ.push_back(targetProjectionZ->GetBinContent(i));
	}

	delete targetProjectionX;
	delete targetProjectionY;
	delete targetProjectionZ;

	TH1D* generatedProjectionX = generatedDist->ProjectionX();
	TH1D* generatedProjectionY = generatedDist->ProjectionY();
	TH1D* generatedProjectionZ = generatedDist->ProjectionZ();
	generatedProjectionX->Scale(1.0 / generatedProjectionX->Integral());
	generatedProjectionY->Scale(1.0 / generatedProjectionY->Integral());
	generatedProjectionZ->Scale(1.0 / generatedProjectionZ->Integral());

	for (int i = 1; i <= generatedProjectionX->GetNbinsX(); i++)
	{
		generatedProjectionValuesX.push_back(generatedProjectionX->GetBinContent(i));
	}
	for (int i = 1; i <= generatedProjectionY->GetNbinsX(); i++)
	{
		generatedProjectionValuesY.push_back(generatedProjectionY->GetBinContent(i));
	}
	for (int i = 1; i <= generatedProjectionZ->GetNbinsX(); i++)
	{
		generatedProjectionValuesZ.push_back(generatedProjectionZ->GetBinContent(i));
	}

	delete generatedProjectionX;
	delete generatedProjectionY;
	delete generatedProjectionZ;
}

MCMC_Data::MCMC_Data(MCMC_Data&& other) noexcept
{
	targetDist = other.targetDist;
	generatedDist = other.generatedDist;
	other.targetDist = nullptr;
	other.generatedDist = nullptr;

	xAxis = std::move(other.xAxis);
	yAxis = std::move(other.yAxis);
	zAxis = std::move(other.zAxis);

	targetProjectionValuesX = std::move(other.targetProjectionValuesX);
	targetProjectionValuesY = std::move(other.targetProjectionValuesY);
	targetProjectionValuesZ = std::move(other.targetProjectionValuesZ);

	generatedProjectionValuesX = std::move(other.generatedProjectionValuesX);
	generatedProjectionValuesY = std::move(other.generatedProjectionValuesY);
	generatedProjectionValuesZ = std::move(other.generatedProjectionValuesZ);

	targetSlice = std::move(other.targetSlice);
	generatedSlice = std::move(other.generatedSlice);

	label = std::move(other.label);
}

MCMC_Data& MCMC_Data::operator=(MCMC_Data&& other) noexcept
{
	if (this == &other) return *this;

	delete targetDist;
	delete generatedDist;

	targetDist = other.targetDist;
	generatedDist = other.generatedDist;
	other.targetDist = nullptr;
	other.generatedDist = nullptr;

	xAxis = std::move(other.xAxis);
	yAxis = std::move(other.yAxis);
	zAxis = std::move(other.zAxis);

	targetProjectionValuesX = std::move(other.targetProjectionValuesX);
	targetProjectionValuesY = std::move(other.targetProjectionValuesY);
	targetProjectionValuesZ = std::move(other.targetProjectionValuesZ);

	generatedProjectionValuesX = std::move(other.generatedProjectionValuesX);
	generatedProjectionValuesY = std::move(other.generatedProjectionValuesY);
	generatedProjectionValuesZ = std::move(other.generatedProjectionValuesZ);

	targetSlice = std::move(other.targetSlice);
	generatedSlice = std::move(other.generatedSlice);

	label = std::move(other.label);

	return *this;
}
