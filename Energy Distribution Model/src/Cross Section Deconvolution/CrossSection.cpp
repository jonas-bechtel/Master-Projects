#include "pch.h"

#include "CrossSection.h"
#include "Constants.h"
#include "EnergyDistribution.h"
#include "RateCoefficient.h"
#include "PlasmaRateCoefficient.h"
#include "EnergyDistributionSet.h"
#include "DeconvolutionWindow.h"

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
	double kT_trans = 0.002; // representativeEnergyDist.eBeamParameter.transverse_kT;
	double kT_long = 0.00047; // representativeEnergyDist->eBeamParameter.longitudinal_kT; //eBeam->GetLongitudinal_kT(labEnergiesParameter.centerLabEnergy);

	// first edge needs to be 0
	double minEnergy = 0;    
	double maxEnergy = 100;  //just a guess, not fixed
	double secondEdge = kT_trans / 20;

	Clear();

	std::vector<double> binEdges;

	binEdges.reserve(10);
	binEdges.push_back(minEnergy);
	binEdges.push_back(secondEdge);

	// first bins are always the same
	for (int i = 1; binEdges[i] < kT_trans; i++)
	{
		binEdges.push_back(2 * binEdges.back() * binSettings.binFactor);
	}

	// binning like in the paper 
	if (binSettings.scheme == PaperBinning)
	{
		int binFactor = 1;

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
		for (size_t i = 0; i < rc.detuningEnergies.size() - 1; i++)
		{
			if (rc.detuningEnergies.at(i) < lastEdgeSoFar) continue;

			double Ed_i = rc.detuningEnergies.at(i);
			double Ed_after_i = rc.detuningEnergies.at(i + 1);
			binEdges.push_back((Ed_i + Ed_after_i) / 2);
		}
		// add one last edge
		binEdges.push_back(binEdges.back() + 2 * (rc.detuningEnergies.back() - binEdges.back()));
	}

	// bin width increases by a constant factor
	if (binSettings.scheme == FactorBinning)
	{
		binEdges.reserve(binEdges.size() + binSettings.numberBins);
		double factor = TMath::Power((maxEnergy / binEdges.back()), (1.0 / (binSettings.numberBins)));

		for (int i = 0; i < binSettings.numberBins; i++)
		{
			binEdges.push_back(binEdges.back() * factor);
		}
	}

	// start with factor binning, then use rc points
	if (binSettings.scheme == PaperFactorMix)
	{
		binEdges.reserve(binEdges.size() + binSettings.numberBins);
		double factor = TMath::Power((binSettings.boundaryEnergy / binEdges.back()), (1.0 / (binSettings.numberBins)));

		for (int i = 0; i < binSettings.numberBins; i++)
		{
			binEdges.push_back(binEdges.back() * factor);
		}

		double lastEdgeSoFar = binEdges.back();

		// add edges so rc points are in bin center (rc is sorted in ascending order)
		for (size_t i = 0; i < rc.detuningEnergies.size() - 1; i++)
		{
			if (rc.detuningEnergies.at(i) < lastEdgeSoFar) continue;

			double Ed_i = rc.detuningEnergies.at(i);
			double Ed_after_i = rc.detuningEnergies.at(i + 1);
			binEdges.push_back((Ed_i + Ed_after_i) / 2);
			std::cout << binEdges.back() << std::endl;
		}
		// add one last edge
		binEdges.push_back(binEdges.back() + 2 * (rc.detuningEnergies.back() - binEdges.back()));
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

	values.resize(hist->GetNbinsX());
	errors.resize(hist->GetNbinsX());
	energies.resize(hist->GetNbinsX());

	for (int i = 1; i <= hist->GetNbinsX(); i++)
	{
		energies.at(i - 1) = hist->GetBinCenter(i);
	}
}

void CrossSection::SetInitialGuessValues(const RateCoefficient& rc)
{
	for (int i = 0; i < values.size(); i++)
	{
		double energy = hist->GetBinCenter(i + 1);
		double velocity = TMath::Sqrt(2 * energy * TMath::Qe() / PhysicalConstants::electronMass);
		double alpha = rc.graph->Eval(energy);
		
		values.at(i) = (alpha / velocity);
	}
}

double CrossSection::ConvolveFit(double Ed, double* csBins, const EnergyDistributionSet& set, bool squareCS,
	std::unordered_map<double, EnergyDistribution*>& map) const
{
	double sum = 0;
	const EnergyDistribution* distribution = nullptr;

	//check if Ed is in map
	if (map.find(Ed) != map.end())
	{
		distribution = map.at(Ed);
	}
	else
	{
		distribution = set.FindByEd(Ed);
		map.insert({ Ed, const_cast<EnergyDistribution*>(distribution) });
	}
	
	if (distribution == nullptr)
	{
		std::cout << "did not find distribution for detuning energy " << Ed << std::endl;
		return 0.0;
	}

	for (int i = 0; i < distribution->psi.size(); i++)
	{
		if (squareCS)
			sum += distribution->psi[i] * csBins[i] * csBins[i];

		else
			sum += distribution->psi[i] * csBins[i];
	}

	return sum;
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
		alphaVector[i] = rc.graph->GetPointY(i);
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

void CrossSection::FitWithROOT(const RateCoefficient& rc, const EnergyDistributionSet& set, const FittingOptions& fitSettings)
{
	// temporary cache of energy dists for deconvolution
	std::unordered_map<double, EnergyDistribution*> EdToDistMap;

	TF1 fitFunction("fit function", 
		[this, &set, &EdToDistMap](double* x, double* p)
		{
			return this->ConvolveFit(x[0], p, set, true, EdToDistMap);
		}, 
		0, 100, hist->GetNbinsX());

	//ROOT::Math::MinimizerOptions::SetDefaultMinimizer("GSL");

	std::vector<double> squareRootValues;
	squareRootValues.reserve(values.size());
	for (double value : values)
	{
		squareRootValues.push_back(sqrt(abs(value)));
	}
	fitFunction.SetParameters(squareRootValues.data());
	
	int numberFixedParameters = 0;
	if (fitSettings.fixParameters)
	{
		int fixedParamLowIndex = hist->FindBin(fitSettings.fixedParameterRange[0]) - 1;
		int fixedParamHighIndex = hist->FindBin(fitSettings.fixedParameterRange[1]) - 1;
		//std::cout << "fixed parameters: " << fixedParamLowIndex << " to " << fixedParamHighIndex << std::endl;
		//std::cout << "lower fixed parameter value:" << hist->GetBinCenter(fixedParamLowIndex + 1) << std::endl;
		//std::cout << "upper fixed parameter value:" << hist->GetBinCenter(fixedParamHighIndex + 1) << std::endl;
		for (size_t i = fixedParamLowIndex; i <= fixedParamHighIndex; i++)
		{
			int bla = std::min(i, squareRootValues.size() - 1);
			fitFunction.FixParameter(i, squareRootValues.at(bla));
			numberFixedParameters++;
		}
	}
	rc.graph->Fit(&fitFunction, "RNM");

	double* result = fitFunction.GetParameters();
	for (int i = 0; i < values.size(); i++)
	{
		hist->SetBinContent(i + 1, result[i] * result[i]);
		values[i] = result[i] * result[i];
	}

	//int N = rc.graph->GetN();
	//int p = fitFunction.GetNpar();
	//int DOF = N - (p - numberFixedParameters);
	//std::cout << "N: " << N << ", p_tot: " << p << ", p_fix: " << numberFixedParameters << ", DOF: " << DOF << std::endl;
	//double chi2 = fitFunction.GetChisquare();
	//
	//double reduced_chi2 = chi2 / DOF;
	//std::cout << "Reduced chi2: " << reduced_chi2 << std::endl;
}

void CrossSection::FitWithEigenNNLS(const RateCoefficient& rc, const EnergyDistributionSet& set, const FittingOptions& fitSettings)
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
		alphaVector[i] = rc.graph->GetPointY(i);
	}

	auto x = Eigen::NNLS(PsiMatrix, fitSettings.maxIterations, fitSettings.tolerance);
	auto res = x.solve(alphaVector);
	//std::cout << "number of iterations: " << x.iterations() << std::endl;
	//std::cout << res << std::endl;

	for (int i = 0; i < values.size(); i++)
	{
		hist->SetBinContent(i + 1, res[i]);
		values[i] = res[i];
	}
}

void CrossSection::ResetNonFixedParameters(const RateCoefficient& rc, const FittingOptions& fitSettings)
{
	int fixedParamLowIndex = hist->FindBin(fitSettings.fixedParameterRange[0]) - 1;
	int fixedParamHighIndex = hist->FindBin(fitSettings.fixedParameterRange[1]) - 1;
	for (int i = 0; i < values.size(); i++)
	{
		//std::cout << "value[" << i << "] = " << values.at(i) << std::endl;
		if (i < fixedParamLowIndex || i > fixedParamHighIndex)
		{
			double energy = hist->GetBinCenter(i + 1);
			double velocity = TMath::Sqrt(2 * energy * TMath::Qe() / PhysicalConstants::electronMass);
			double alpha = rc.graph->Eval(energy);

			values.at(i) = alpha / velocity;
		}
		//std::cout << "value[" << i << "] = " << values.at(i) << std::endl;
	}
}

void CrossSection::FillWithOneOverE(double scale)
{
	Clear();

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
		valueArray.push_back(value);
		hist->SetBinContent(i, value);
	}
	std::ostringstream oss;
	oss << std::scientific << std::setprecision(2) << scale;
	std::string scaleString = oss.str();
	label = scaleString + std::string(" over E cs");
}

void CrossSection::Deconvolve(RateCoefficient& rc, EnergyDistributionSet& set, const FittingOptions& fitSettings, const CrossSectionBinningSettings& binSettings)
{
	if (rc.value.size() != set.distributions.size())
	{
		std::cout << "sizes of rate coefficients and energy distributions dont match: " <<
			rc.value.size() << " != " << set.distributions.size() << std::endl;
		return;
	}
	energyDistriubtionSetFolder = set.Label();
	mergedBeamRateCoefficientFile = rc.label;
	
	SetupBinning(binSettings, rc);

	set.CalculatePsisFromBinning(hist);

	const int arraySize = fitSettings.errorIterations * hist->GetNbinsX();
	valueArray.clear();
	valueArray.resize(arraySize);
     
	for (int i = 0; i < fitSettings.errorIterations; i++)
	{
		SetInitialGuessValues(rc);

		if (fitSettings.ROOT_fit)
		{
			for (int i = 0; i < fitSettings.fit_iterations; i++)
			{
				FitWithROOT(rc, set, fitSettings);
			}
		}
		else if (fitSettings.SVD_fit)
		{
			FitWithSVD(rc, set);
		}
		else if (fitSettings.EigenNNLS_fit)
		{
			FitWithEigenNNLS(rc, set, fitSettings);
		}
		else if (fitSettings.NNLS_ROOT_fit)
		{
			FitWithEigenNNLS(rc, set, fitSettings);

			if (fitSettings.fixParameters)
			{
				ResetNonFixedParameters(rc, fitSettings);
			}

			for (int i = 0; i < fitSettings.fit_iterations; i++)
			{
				FitWithROOT(rc, set, fitSettings);
			}
		}

		for(int j = 0; j < values.size(); j++)
		{
			valueArray[i + j * fitSettings.errorIterations] = values[j];
		}

		rc.VaryGraphValues();
	}
	rc.ResetGraphValues();

	// Calculate mean and errors of error iterations
	for (int j = 0; j < values.size(); j++)
	{
		double mean = 0;
		double error = 0;

		for (int i = 0; i < fitSettings.errorIterations; i++)
		{
			double value = valueArray[i + j * fitSettings.errorIterations];
			mean += value;
		}
		mean /= fitSettings.errorIterations;

		for (int i = 0; i < fitSettings.errorIterations; i++)
		{
			error += pow(valueArray[i + j * fitSettings.errorIterations] - mean, 2);
		}
		error = sqrt(error / (fitSettings.errorIterations - 1));

		values[j] = mean;
		errors[j] = error;
		hist->SetBinContent(j + 1, mean);
		hist->SetBinError(j + 1, error);
	}

	// calculate chi2 of fit
	double chi2 = 0;
	for (int i = 0; i < rc.graph->GetN(); i++)
	{
		double fitValue = ConvolveFit(rc.graph->GetPointX(i), values.data(), set, false);
		double error = rc.graph->GetErrorY(i);
		//std::cout << "Point " << i << ": x = " << rc.graph->GetPointX(i) 
		//	<< ", fit value = " << fitValue << ", error = " << error << std::endl;
		if (error > 0)
		{
			double diff = rc.graph->GetPointY(i) - fitValue;
			chi2 += diff * diff / (error * error);
		}
		else
		{
			std::cout << "Error is zero for point " << i << ", skipping chi2 calculation for this point." << std::endl;
		}
	}
	double chi2_reduced = chi2 / (rc.graph->GetN() - hist->GetNbinsX());
	std::cout << "Chi2 reduced: " << chi2_reduced << std::endl;
}

void CrossSection::Plot(bool showMarkers) const
{
	if (showMarkers) ImPlot::SetNextMarkerStyle(ImPlotMarker_Square);
	ImPlot::PlotLine(label.c_str(), energies.data(), values.data(), values.size());
	ImPlot::PlotErrorBars(label.c_str(), energies.data(), values.data(), errors.data(), errors.size());
}

void CrossSection::Clear()
{
	delete hist;
	hist = nullptr;
	energies.clear();
	values.clear();
	errors.clear();
	valueArray.clear();
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

	Clear();

	std::string line;
	// skip first line
	std::getline(infile, line);
	while (std::getline(infile, line))
	{
		std::vector<std::string> tokens = FileUtils::SplitLine(line, "\t");
		energies.push_back(std::stod(tokens[0]));
		values.push_back(std::stod(tokens[1]));
		errors.push_back(std::stod(tokens[2]));

		valueArray.push_back(std::stod(tokens[1]));
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
	outfile << std::scientific << std::setprecision(6);

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
		ImGui::Checkbox("Eigen NNLS fitting", &EigenNNLS_fit);
		ImGui::Checkbox("NNLS ROOT combo", &NNLS_ROOT_fit);
		
		ImGui::InputInt("iterations", &fit_iterations);
		ImGui::InputInt("max iterations", &maxIterations);
		ImGui::InputDouble("tolerance", &tolerance, 0.0, 0.0, "%.1e");
		//ImGui::InputDouble("learning rate", &learningRate);
		ImGui::Checkbox("fix params", &fixParameters);
		ImGui::InputFloat2("fixed parameter range", fixedParameterRange, "%.4f");
		ImGui::InputInt("error iterations", &errorIterations);
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
		ImGui::SetNextItemWidth(150.0f);
		ImGui::Combo("binning options", (int*)&scheme, binningOptions, IM_ARRAYSIZE(binningOptions));
		ImGui::PushItemWidth(100.0f);
		ImGui::InputDouble("bin factor", &binFactor, 0.0, 0.0, "%.4f");
		if (scheme == FactorBinning || scheme == PaperFactorMix)
		{
			ImGui::InputInt("number bins", &numberBins);
		}
		if (scheme == PaperFactorMix)
		{
			ImGui::InputDouble("boundary energy", &boundaryEnergy, 0.0, 0.0, "%.4f");
		}
		if (scheme == PaperBinning)
		{
			ImGui::InputInt("max ratio", &maxRatio);
		}
		ImGui::PopItemWidth();
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
