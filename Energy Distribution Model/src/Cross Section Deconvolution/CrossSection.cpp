#include "pch.h"

#include "Eigen/Dense"

#include "CrossSection.h"
#include "Constants.h"
#include "EnergyDistribution.h"
#include "RateCoefficient.h"
#include "EnergyDistributionSet.h"

#include "FileUtils.h"

CrossSection::CrossSection()
{

}

TH1D* CrossSection::GetHist()
{
	return hist;
}

std::string CrossSection::GetLabel() const
{
	return label;
}

void CrossSection::SetLabel(std::string newLabel)
{
	label = newLabel;
}

void CrossSection::SetupBinning(const CrossSectionBinningSettings& binSettings, const RateCoefficient& rc)
{
	// first edge needs to be 0
	std::vector<double> binEdges;
	double minEnergy = 0;    
	double maxEnergy = 100;  //just a guess, not fixed

	// binning like in the paper 
	if (binSettings.scheme == PaperBinning)
	{
		//EnergyDistribution& representativeEnergyDist = energyDistributionList.back();
		double kT_trans = 0.002; // representativeEnergyDist.eBeamParameter.transverse_kT;
		double kT_long = 0.00047; // representativeEnergyDist->eBeamParameter.longitudinal_kT; //eBeam->GetLongitudinal_kT(labEnergiesParameter.centerLabEnergy);

		binEdges.push_back(0);
		binEdges.push_back(kT_trans / 20);
		int binFactor = 1;

		for (int i = 1; binEdges[i] < kT_trans; i++)
		{
			binEdges.push_back(binFactor * 2 * binEdges.back());
		}
		while (binEdges.back() < maxEnergy)
		{
			double previousEdge = binEdges.back();
			double delta_E = sqrt(pow((kT_trans * log(2)), 2) + 16 * log(2) * kT_long * previousEdge);
			binEdges.push_back(previousEdge + binFactor * delta_E);

			//std::cout << "delta E " << delta_E << ", last edge:	" << previousEdge << ", ratio: " << previousEdge /delta_E << "\n";
			if (previousEdge / delta_E > binSettings.maxRatio) break;
		}
		double lastEdgeSoFar = binEdges.back();

		// add edges so rc points are in bin center
		for (size_t i = rc.detuningEnergies.size() - 1; i > 0; i--)
		{
			if (rc.detuningEnergies.at(i) < lastEdgeSoFar) continue;
			
			double Ed_i = rc.detuningEnergies.at(i);
			double Ed_after_i = rc.detuningEnergies.at(i - 1);
			//std::cout << "Ed[i] = " << Ed_i << " Ed[i + 1] = " << Ed_after_i << " i = " << i << std::endl;
			binEdges.push_back((Ed_i + Ed_after_i) / 2);
		}
		// add one last edge
		binEdges.push_back(binEdges.back() + 2 * (rc.detuningEnergies.front() - binEdges.back()));
	}
	// bin width increses by a constant factor
	if (binSettings.scheme == FactorBinning)
	{
		double min = minEnergy;
		double factor = TMath::Power((maxEnergy / min), (1.0 / binSettings.numberBins));

		binEdges.reserve(binSettings.numberBins + 1);
		binEdges.push_back(min);
		for (int i = 0; i < binSettings.numberBins; i++)
		{
			binEdges.push_back(binEdges[i] * factor);
		}
	}
	if (binSettings.scheme == PaperFactorMix)
	{
		//EnergyDistribution& representativeEnergyDist = energyDistributionList.back();
		double kT_trans = 0.002; // representativeEnergyDist.eBeamParameter.transverse_kT;
		//double kT_long = 0.00047; // representativeEnergyDist->eBeamParameter.longitudinal_kT; //eBeam->GetLongitudinal_kT(labEnergiesParameter.centerLabEnergy);

		binEdges.push_back(0);
		binEdges.push_back(kT_trans / 20);

		for (int i = 1; binEdges[i] < kT_trans; i++)
		{
			binEdges.push_back(2 * binEdges[i]);
			//std::cout << binEdges[i + 1] << "\n";
		}

		double factor = TMath::Power((maxEnergy / binEdges.back()), (1.0 / binSettings.numberBins));

		for (int i = 0; i < binSettings.numberBins; i++)
		{
			binEdges.push_back(binEdges.back() * factor);
		}
	}

	//for (double edge : binEdges)
	//{
	//	std::cout << edge << "\n";
	//}
	for (int i = 1; i < binEdges.size(); i++)
	{
		std::cout << "bin center: " << (binEdges[i] + binEdges[i - 1]) / 2 << "\tbin width:	" << (binEdges[i] - binEdges[i - 1]) << std::endl;
	}
	std::cout << "number cross section bins: " << binEdges.size() - 1 << "\n";

	hist = new TH1D("cross section fit", "cross section fit", binEdges.size() - 1, binEdges.data());
}

void CrossSection::SetupInitialGuess(const RateCoefficient& rc)
{
	if (!values.empty())
	{
		return;
	}
	values.clear();
	values.reserve(hist->GetNbinsX());
	errors.clear();
	errors.resize(hist->GetNbinsX());

	for (int i = 1; i <= hist->GetNbinsX(); i++)
	{
		double energy = hist->GetBinCenter(i);
		double velocity = TMath::Sqrt(2 * energy * TMath::Qe() / PhysicalConstants::electronMass);
		double alpha = rc.graph->Eval(energy);
		
		values.push_back(alpha / velocity);
		
		//initialGuess.push_back(alpha / velocity);
		energies.push_back(energy);
	}
}

void CrossSection::FitWithSVD(const RateCoefficient& rc, const EnergyDistributionSet& set)
{
	int n = set.distributions.size();
	int p = hist->GetNbinsX();

	Eigen::MatrixXd PsiMatrix(n, p);
	Eigen::VectorXd alphaVector(n);

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < p; j++)
		{
			// fill matrix and vector with 0 if p > n
			PsiMatrix(i, j) = set.distributions[i].psi[j];
		}
		alphaVector[i] = rc.value[i];
	}

	Eigen::JacobiSVD<Eigen::MatrixXd> svd(PsiMatrix, Eigen::ComputeThinU | Eigen::ComputeThinV);

	// Solve using the pseudoinverse (x = A^+ * b), where A^+ is the Moore-Penrose pseudoinverse
	Eigen::VectorXd result = svd.solve(alphaVector);

	for (int i = 0; i < values.size(); i++)
	{
		hist->SetBinContent(i + 1, result[i]);
		values[i] = result[i];
	}
}

void CrossSection::FitWithROOT(const RateCoefficient& rc, const EnergyDistributionSet& set)
{
	TF1 fitFunction("fit function", 
		[&rc, &set](double* x, double* p) 
		{
			return rc.ConvolveFit(x[0], p, set);
		}, 
		0, 100, hist->GetNbinsX());

	//ROOT::Math::MinimizerOptions::SetDefaultMinimizer("GSL");

	std::vector<double> squareRootValues;
	squareRootValues.reserve(values.size());
	for (double value : values)
	{
		squareRootValues.push_back(sqrt(value));
	}
	fitFunction.SetParameters(squareRootValues.data());
	
	rc.graph->Fit(&fitFunction, "RN");

	double* result = fitFunction.GetParameters();
	for (int i = 0; i < values.size(); i++)
	{
		hist->SetBinContent(i + 1, result[i] * result[i]);
		values[i] = result[i] * result[i];
	}
}

void CrossSection::FillWithOneOverE(int scale)
{
	// setup binning
	std::vector<double> binEdges;
	double minEnergy = 1e-4;
	double maxEnergy = 100;
	int numberBins = 1000;

	double factor = TMath::Power((maxEnergy / minEnergy), (1.0 / numberBins));

	binEdges.reserve(numberBins + 2);

	binEdges.push_back(0);
	binEdges.push_back(minEnergy);
	for (int i = 0; i < numberBins; i++)
	{
		binEdges.push_back(binEdges.back() * factor);
	}
	hist = new TH1D("cross section 1/E", "cross section 1/E", binEdges.size() - 1, binEdges.data());

	// fill hist and lists
	values.reserve(hist->GetNbinsX());
	energies.reserve(hist->GetNbinsX());

	for (int i = 1; i <= hist->GetNbinsX(); i++)
	{
		double energy = hist->GetBinCenter(i);
		double value = scale / energy;
		
		energies.push_back(energy);
		values.push_back(value);
		hist->SetBinContent(i, value);
	}
	label = scale + std::string(" over E cs");
}

void CrossSection::Deconvolve(const RateCoefficient& rc, EnergyDistributionSet& set, const FittingOptions& fitSettings, const CrossSectionBinningSettings& binSettings)
{
	if (rc.value.size() != set.distributions.size())
	{
		std::cout << "sizes of rate coefficients and energy distributions dont match: " <<
			rc.value.size() << " != " << set.distributions.size() << std::endl;
	}
	energyDistriubtionSetFolder = set.Label();
	mergedBeamRateCoefficientFile = rc.label;
	
	SetupBinning(binSettings, rc);

	set.CalculatePsisFromBinning(hist);
	SetupInitialGuess(rc);

	if (fitSettings.ROOT_fit)
	{
		FitWithROOT(rc, set);
	}
	else if (fitSettings.SVD_fit)
	{
		FitWithSVD(rc, set);
	}
}

void CrossSection::Plot(bool showMarkers) const
{
	if (showMarkers) ImPlot::SetNextMarkerStyle(ImPlotMarker_Square);
	ImPlot::PlotLine(label.c_str(), energies.data(), values.data(), values.size());
	ImPlot::PlotErrorBars(label.c_str(), energies.data(), values.data(), errors.data(), errors.size());
}

void CrossSection::Load(std::filesystem::path& file)
{
	std::ifstream infile(file);

	// Check if the file was successfully opened
	if (!infile.is_open())
	{
		std::cerr << "Error: Could not open the file " << file << std::endl;
		return;
	}

	std::string line;
	// skip first line
	std::getline(infile, line);
	while (std::getline(infile, line))
	{
		std::vector<std::string> tokens = FileUtils::SplitLine(line, "\t");
		energies.push_back(std::stod(tokens[0]));
		values.push_back(std::stod(tokens[1]));
		errors.push_back(std::stod(tokens[2]));
	}

	std::vector<double> binEdges = FileUtils::CalculateBinEdges(energies, false);

	hist = new TH1D("cross section fit", "cross section fit", binEdges.size() - 1, binEdges.data());

	for (int i = 1; i <= hist->GetNbinsX(); i++)
	{
		hist->SetBinContent(i, values.at(i - 1));
		hist->SetBinError(i, errors.at(i - 1));
	}
}

void CrossSection::Save() const
{
	// set the output filepath
	std::filesystem::path file = FileUtils::GetCrossSectionFolder() / (label + ".dat");

	// Create the directories if they don't exist
	if (!std::filesystem::exists(file.parent_path()))
	{
		std::filesystem::create_directories(file.parent_path());
	}

	std::ofstream outfile(file);

	if (!outfile.is_open())
	{
		std::cerr << "Error opening file" << std::endl;
		return;
	}
	outfile << std::setprecision(std::numeric_limits<double>::digits10 + 1) << std::fixed;

	outfile << "# Energy [eV]\tCross Section Value\terror\n";

	for (int i = 0; i < energies.size(); i++)
	{
		outfile << energies[i] << "\t" << values[i] << "\t" << errors[i] << "\n";
	}

	outfile.close();
}

void FittingOptions::ShowWindow(bool& show)
{
	if (!show)
	{
		return;
	}
	if (ImGui::Begin("Cross Section Fit settings", &show, ImGuiWindowFlags_NoDocking))
	{
		ImGui::Checkbox("ROOT fitting", &ROOT_fit);
		ImGui::SameLine();
		ImGui::Checkbox("SVD fitting", &SVD_fit);
		//ImGui::Checkbox("GD fitting", &GD_fit);
		//ImGui::InputInt("iterations", &iterations);
		//ImGui::InputDouble("learning rate", &learningRate);

	}
	ImGui::End();
}

void CrossSectionBinningSettings::ShowWindow(bool& show)
{
	if (!show)
	{
		return;
	}
	if (ImGui::Begin("Cross Section Binning settings", &show, ImGuiWindowFlags_NoDocking))
	{
		ImGui::Combo("binning options", (int*)&scheme, binningOptions, IM_ARRAYSIZE(binningOptions));

		if (scheme == FactorBinning || scheme == PaperFactorMix)
		{
			ImGui::SameLine();
			ImGui::SetNextItemWidth(100.0f);
			ImGui::InputInt("number bins", &numberBins);
		}
		if (scheme == FactorBinning)
		{
			//ImGui::SameLine();
			//ImGui::Checkbox("limit bin size", &limitBinSize);
			//ImGui::SameLine();
			//ImGui::BeginDisabled(!limitBinSize);
			//ImGui::SetNextItemWidth(100.0f);
			//ImGui::InputDouble("min bin size", &minBinSize, 0.0, 0.0, "%.1e");
			//ImGui::EndDisabled();
		}
		if (scheme == PaperBinning)
		{
			ImGui::SameLine();
			ImGui::InputInt("max ration", &maxRatio);
			//ImGui::InputInt("factor", &binFactor);
		}
	}
	ImGui::End();
}



//void CrossSectionManager::FitWithSVD()
//{
//	// solve Ax = b with SVD
//
//	EnergyDistributionManager* model = (EnergyDistributionManager*)Window::Get("Energy Distribution Manager");
//	std::vector<EnergyDistribution>& energyDistributionList = model->GetEnergyDistributions();
//
//	if (model->GetEnergyDistributions().empty())
//	{
//		std::cout << "no energy distributions\n";
//		return;
//	}
//
//	SetupFitCrossSectionHist();
//	CalculatePsis();
//
//	binValuesFit.clear();
//	binCentersFit.clear();
//	binValuesFit.reserve(crossSectionFit->GetNbinsX());
//	binCentersFit.reserve(crossSectionFit->GetNbinsX());
//
//	std::vector<int> parameterIndeces;
//	for (int i = 0; i < crossSectionFit->GetNbinsX(); i++)
//	{
//		parameterIndeces.push_back(i);
//	}
//	int nZeros = 10;
//
//	TMatrixD PsiMatrix;
//	TVectorD alphaVector;
//	int p = parameterIndeces.size(); // crossSectionFit->GetNbinsX();
//	std::vector<double> parameterResult;
//	parameterResult.resize(p);
//
//	while (nZeros > 0)
//	{
//		int n = energyDistributionList.size();
//		int p2 = parameterIndeces.size();
//
//		std::cout << "n: " << n << " p2: " << p2 << std::endl;
//		PsiMatrix.Clear();
//		alphaVector.Clear();
//
//		// matrix A is all the p Psis for all the n distributions, n >= p is required so the rest is filled with 0
//		PsiMatrix.ResizeTo(std::max(n, p2), p2);
//		// vector b with all the rate coeffictions
//		alphaVector.ResizeTo(std::max(n, p2));
//
//		// fill matrix and vector
//		for (int i = 0; i < std::max(n, p2); i++)
//		{
//			for (int j = 0; j < p2; j++)
//			{
//				// fill matrix and vector with 0 if p > n
//				if (i >= n)
//				{
//					PsiMatrix[i][j] = 0;
//				}
//				else
//				{
//					PsiMatrix[i][j] = energyDistributionList[i].psi[parameterIndeces[j]];
//				}
//			}
//			if (i >= n)
//			{
//				alphaVector[i] = 0;
//			}
//			else
//			{
//				alphaVector[i] = energyDistributionList[i].rateCoefficient;
//			}
//		}
//
//		alphaVector.Print();
//
//		TDecompSVD decomp(PsiMatrix);
//		// solve for x, result is put into input vector
//		decomp.Solve(alphaVector);
//
//		alphaVector.Print();
//
//		// remove all negative parameters
//		nZeros = 0;
//		for (int i = p2 - 1; i >= 0; i--)
//		{
//			if (alphaVector[parameterIndeces[i]] < 0)
//			{
//				parameterResult[parameterIndeces[i]] = 0;
//				parameterIndeces.erase(parameterIndeces.begin() + i);
//				nZeros++;
//			}
//			else
//			{
//				parameterResult[parameterIndeces[i]] = alphaVector[parameterIndeces[i]];
//			}
//		}
//		std::cout << "zeros: " << nZeros << std::endl;
//	}
//
//	FillFitPlots(parameterResult.data());
//}
