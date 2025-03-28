#include "pch.h"
#include "CoolingForceValue.h"


#include "CoolingForceModel.h"

#include "Constants.h"
#include "FileUtils.h"

RNG_engine CoolingForceValue::generator = RNG_engine();

std::normal_distribution<double> CoolingForceValue::longitudinalNormalDistribution = std::normal_distribution<double>();
std::normal_distribution<double> CoolingForceValue::transverseNormalDistribution = std::normal_distribution<double>();


CoolingForceValue::CoolingForceValue()
{
	//std::cout << "calling cf value default constructor" << std::endl;
}

CoolingForceValue::~CoolingForceValue()
{
	//std::cout << "calling cf value destructor" << std::endl;
	delete forceX;
	delete forceY;
	delete forceZ;
	delete positionSamples;
}

CoolingForceValue::CoolingForceValue(CoolingForceValue&& other) noexcept
{
	//std::cout << "calling cf value move constructor" << std::endl;
	forceX = other.forceX;
	forceY = other.forceY;
	forceZ = other.forceZ;

	positionSamples = other.positionSamples;

	forceXSlice = std::move(other.forceXSlice);
	forceYSlice = std::move(other.forceYSlice);
	forceZSlice = std::move(other.forceZSlice);

	forceXIntegral = std::move(other.forceXIntegral);
	forceYIntegral = std::move(other.forceYIntegral);
	forceZIntegral = std::move(other.forceZIntegral);

	forceXValue = std::move(other.forceXValue);
	forceYValue = std::move(other.forceYValue);
	forceZValue = std::move(other.forceZValue);

	xAxis = std::move(other.xAxis);
	yAxis = std::move(other.yAxis);
	zAxis = std::move(other.zAxis);

	forceZProjectionX = std::move(other.forceZProjectionX);
	forceZProjectionY = std::move(other.forceZProjectionY);
	forceZProjectionZ = std::move(other.forceZProjectionZ);

	eBeamParameter = other.eBeamParameter;
	ionBeamParameter = other.ionBeamParameter;
	labEnergiesParameter = other.labEnergiesParameter;

	label = std::move(other.label);
	tags = std::move(other.tags);
	index = std::move(other.index);

	other.ResetDefaultValues();
}

CoolingForceValue& CoolingForceValue::operator=(CoolingForceValue&& other) noexcept
{
	//std::cout << "calling cf value move assignment" << std::endl;
	if (this == &other) return *this;

	forceX = other.forceX;
	forceY = other.forceY;
	forceZ = other.forceZ;

	positionSamples = other.positionSamples;

	forceXSlice = std::move(other.forceXSlice);
	forceYSlice = std::move(other.forceYSlice);
	forceZSlice = std::move(other.forceZSlice);

	forceXIntegral = std::move(other.forceXIntegral);
	forceYIntegral = std::move(other.forceYIntegral);
	forceZIntegral = std::move(other.forceZIntegral);

	forceXValue = std::move(other.forceXValue);
	forceYValue = std::move(other.forceYValue);
	forceZValue = std::move(other.forceZValue);

	xAxis = std::move(other.xAxis);
	yAxis = std::move(other.yAxis);
	zAxis = std::move(other.zAxis);

	forceZProjectionX = std::move(other.forceZProjectionX);
	forceZProjectionY = std::move(other.forceZProjectionY);
	forceZProjectionZ = std::move(other.forceZProjectionZ);

	eBeamParameter = other.eBeamParameter;
	ionBeamParameter = other.ionBeamParameter;
	labEnergiesParameter = other.labEnergiesParameter;

	label = std::move(other.label);
	tags = std::move(other.tags);
	index = std::move(other.index);

	other.ResetDefaultValues();

	return *this;
}

void CoolingForceValue::Calculate(std::filesystem::path descriptionFile, int index, bool onlyLongInLC)
{
	// get all necessary modules
	std::filesystem::path folder = descriptionFile.parent_path();

	// get 3 parameters: U drift tube, electron current, center E lab if index is in file
	std::array<float, 3> additionalParameter = FileUtils::GetParametersFromDescriptionFileAtIndex(descriptionFile, index);

	// if they are not found the index is not in the file
	if (!additionalParameter[0])
	{
		std::cout << "index " << index << " is not in the file " << descriptionFile.filename() << std::endl;
		return;
	}

	// set read electron current and center lab energy
	LabEnergy::SetDriftTubeVoltage(additionalParameter[0]);
	LabEnergy::SetCenterEnergy(additionalParameter[2]);
	ElectronBeam::SetElectronCurrent(additionalParameter[1]);
	ElectronBeam::CalculateEstimateLongkT();
	ElectronBeam::CalculateDetuningEnergy();
	ElectronBeam::CalculateDetuningVelocity();

	// full procedure to generate one energy distribution 
	// 1. setup necessary distributions
	std::filesystem::path densityfile = FileUtils::FindFileWithIndex(folder / "e-densities", index);
	if (densityfile.empty())
	{
		std::cout << "density file not found: " << densityfile.filename() << std::endl;
		return;
	}
	ElectronBeam::SetupDistribution(densityfile);

	std::filesystem::path energyfile = FileUtils::FindFileWithIndex(folder / "lab-energies", index);
	if (energyfile.empty())
	{
		std::cout << "energy file not found: " << energyfile.filename() << std::endl;
		return;
	}
	LabEnergy::SetupDistribution(energyfile);

	IonBeam::CreateFromReference(ElectronBeam::Get());

	// final setup of current curve
	CopyParameters();
	SetupLabel();
	SetupTags();

	// calculate cooling force
	std::vector<Point3D> samples = IonBeam::GeneratePositions();
	
	if (samples.empty())
	{
		std::cout << "sampling positions failed\n";
		return;
	}
	
	SetupHistograms(IonBeam::Get());
	
	for (const Point3D& point : samples)
	{
		double x = point.x;
		double y = point.y;
		double z = point.z;
	
		// calculate velocity magnitude from lab energy given as matrix (TH3D) from outside
		double labEnergy = LabEnergy::Get(x, y, z);
		double electronVelocityMagnitude = TMath::Sqrt(2 * labEnergy * TMath::Qe() / PhysicalConstants::electronMass);

		// determine direction of velocity based on beam trajectory function
		TVector3 longitudinalDirection = ElectronBeam::GetDirection(point.z);
		TVector3 transverseDirection = longitudinalDirection.Orthogonal();

		// add random values to velocity in transverse and longitudinal directions:
		// - calculate longitudinal kT, transverse kT is fixed
		double long_kT = ElectronBeam::GetLongitudinal_kT(labEnergy);
		double trans_kT = ElectronBeam::GetTransverse_kT();

		// - use kT to calculate sigmas of gaussians
		double longSigma = TMath::Sqrt(long_kT * TMath::Qe() / PhysicalConstants::electronMass);
		double transSigma = TMath::Sqrt(trans_kT * TMath::Qe() / PhysicalConstants::electronMass);

		// - sample from gaussians with these sigmas and add that to the electron velocity
		longitudinalNormalDistribution = std::normal_distribution<double>(0, longSigma);
		transverseNormalDistribution = std::normal_distribution<double>(0, transSigma);
		double longitudinalAddition = longitudinalNormalDistribution(generator);
		double transverseAdditionX = transverseNormalDistribution(generator);
		double transverseAdditionY = transverseNormalDistribution(generator);

		// we need a vector that is never in line with the longitudinalDirection
		TVector3 helpVector = TVector3(1, 0, 0);
		TVector3 transverseDirection1 = longitudinalDirection.Cross(helpVector);
		TVector3 transverseDirection2 = longitudinalDirection.Cross(transverseDirection1);

		TVector3 finalElectronVelocity = (electronVelocityMagnitude + longitudinalAddition) * longitudinalDirection
			+ transverseAdditionX * transverseDirection1
			+ transverseAdditionY * transverseDirection2;

		// calculate collision velocity vector and magnitude using a fixed ion beam velocity
		double ionVelocityMagnitude = TMath::Sqrt(2 * eBeamParameter.coolingEnergy * TMath::Qe() / PhysicalConstants::electronMass); // calc from cooling energy;
		TVector3 ionVelocityDirection = IonBeam::GetDirection();
		TVector3 ionVelocity = ionVelocityDirection * ionVelocityMagnitude;

		TVector3 collisionVelocity = ionVelocity - finalElectronVelocity;
		//double collisionVelocityMagnitude = collisionVelocity.Mag();

		double electronDensity = ElectronBeam::GetDensity(point);
		int ionCharge = IonBeam::GetCharge();
		TVector3 coolingforce = CoolingForceModel::CoolingForce(collisionVelocity, trans_kT, electronDensity, ionCharge, onlyLongInLC);
	
		positionSamples->Fill(x, y, z);
	
		int bin = forceX->FindBin(x, y, z);
		forceX->AddBinContent(bin, coolingforce.x());
		forceY->AddBinContent(bin, coolingforce.y());
		forceZ->AddBinContent(bin, coolingforce.z());
	}
	
	forceX->Divide(positionSamples);
	forceY->Divide(positionSamples);
	forceZ->Divide(positionSamples);
	
	FillData();
}

bool CoolingForceValue::ShowListItem(bool selected) const
{
	std::string labelTags = label;
	if (!tags.empty())
	{
		labelTags += "\n";
		labelTags += tags;
	}

	// Render each item as selectable
	bool clicked = false;
	if (ImGui::Selectable(labelTags.c_str(), selected, ImGuiSelectableFlags_AllowItemOverlap))
	{
		clicked = true;
	}

	if (ImGui::BeginItemTooltip())
	{
		ImGui::PushTextWrapPos(ImGui::GetFontSize() * 35.0f);
		ImGui::TextUnformatted(GetHeaderString().c_str());
		ImGui::PopTextWrapPos();
		ImGui::EndTooltip();
	}

	return clicked;
}

void CoolingForceValue::SetupHistograms(TH3D* reference)
{
	delete forceX;
	delete forceY;
	delete forceZ;
	delete positionSamples;

	forceX = (TH3D*)reference->Clone("cooling force X");
	forceY = (TH3D*)reference->Clone("cooling force Y");
	forceZ = (TH3D*)reference->Clone("cooling force Z");
	positionSamples = (TH3D*)reference->Clone("position samples");

	forceX->Reset();
	forceY->Reset();
	forceZ->Reset();
	positionSamples->Reset();

	forceX->SetTitle("cooling force X");
	forceY->SetTitle("cooling force Y");
	forceZ->SetTitle("cooling force Z");
	positionSamples->SetTitle("position samples");
}

void CoolingForceValue::FillData()
{
	forceXIntegral = CalculateIntegral(forceX);
	forceYIntegral = CalculateIntegral(forceY);
	forceZIntegral = CalculateIntegral(forceZ);

	forceXValue = forceXIntegral / 0.8;
	forceYValue = forceYIntegral / 0.8;
	forceZValue = forceZIntegral / 0.8;

	xAxis.reserve(forceZ->GetNbinsX());
	yAxis.reserve(forceZ->GetNbinsY());
	zAxis.reserve(forceZ->GetNbinsZ());

	forceZProjectionX.reserve(forceZ->GetNbinsX());
	forceZProjectionY.reserve(forceZ->GetNbinsY());
	forceZProjectionZ.reserve(forceZ->GetNbinsZ());

	for (int i = 1; i <= forceZ->GetNbinsX(); i++)
	{
		xAxis.push_back(forceZ->GetXaxis()->GetBinCenter(i));
	}
	for (int i = 1; i <= forceZ->GetNbinsY(); i++)
	{
		yAxis.push_back(forceZ->GetYaxis()->GetBinCenter(i));
	}
	for (int i = 1; i <= forceZ->GetNbinsZ(); i++)
	{
		zAxis.push_back(forceZ->GetZaxis()->GetBinCenter(i));
	}

	TH1D* projectionX = forceZ->ProjectionX();
	TH1D* projectionY = forceZ->ProjectionY();
	TH1D* projectionZ = forceZ->ProjectionZ();

	for (int i = 1; i <= projectionX->GetNbinsX(); i++)
	{
		forceZProjectionX.push_back(projectionX->GetBinContent(i));
	}
	for (int i = 1; i <= projectionY->GetNbinsX(); i++)
	{
		forceZProjectionY.push_back(projectionY->GetBinContent(i));
	}
	for (int i = 1; i <= projectionZ->GetNbinsX(); i++)
	{
		forceZProjectionZ.push_back(projectionZ->GetBinContent(i));
	}

	delete projectionX;
	delete projectionY;
	delete projectionZ;
}

double CoolingForceValue::CalculateIntegral(TH3D* hist)
{
	int nBinsZ = hist->GetZaxis()->GetNbins();
	double zMin = hist->GetZaxis()->GetXmin();
	double zMax = hist->GetZaxis()->GetXmax();

	// Create a 1D histogram for weighted averages along Z
	TH1D* H_weightedAvgZ = new TH1D("H_weightedAvgZ", "Weighted Average per Z", nBinsZ, zMin, zMax);

	// Loop over Z bins
	for (int iZ = 1; iZ <= nBinsZ; iZ++)
	{
		double weightedSum = 0.0;
		double weightSum = 0.0;

		// Loop over all (x, y) bins for this fixed Z bin
		for (int iX = 1; iX <= hist->GetXaxis()->GetNbins(); iX++)
		{
			for (int iY = 1; iY <= hist->GetYaxis()->GetNbins(); iY++)
			{
				int bin = hist->GetBin(iX, iY, iZ);

				double value = hist->GetBinContent(bin);
				double weight = positionSamples->GetBinContent(bin);

				weightedSum += value * weight;
				weightSum += weight;
			}
		}

		// Compute weighted average, handle division by zero
		double avg = (weightSum > 0) ? (weightedSum / weightSum) : 0.0;

		// Store result in the 1D histogram
		H_weightedAvgZ->SetBinContent(iZ, avg);
	}
	double result = H_weightedAvgZ->Integral("width");
	delete H_weightedAvgZ;

	return result;
}

void CoolingForceValue::Save(std::filesystem::path folder) const
{
	std::filesystem::path file = folder / (Filename() + ".asc");
	
	// Open a ROOT file for writing
	TFile outfile(file.string().c_str(), "RECREATE");
	if (!outfile.IsOpen())
	{
		std::cout << "error opening file: " << file << std::endl;
	}
	TNamed header("header", GetHeaderString());
	header.Write();

	// Write histograms into the file
	forceX->Write();
	forceY->Write();
	forceZ->Write();
	positionSamples->Write();

	outfile.Close();
}

void CoolingForceValue::Load(std::filesystem::path file)
{
	TFile infile(file.string().c_str(), "READ");

	// Read metadata
	TNamed* header = (TNamed*)infile.Get("header");
	if (header) 
	{
		eBeamParameter.fromString(header->GetTitle());
		ionBeamParameter.fromString(header->GetTitle());
		labEnergiesParameter.fromString(header->GetTitle());
	}
	SetupLabel();

	// Retrieve histograms
	forceX = (TH3D*)infile.Get("cooling force X");
	forceY = (TH3D*)infile.Get("cooling force Y");
	forceZ = (TH3D*)infile.Get("cooling force Z");
	positionSamples = (TH3D*)infile.Get("position samples");

	forceX->SetDirectory(0);
	forceY->SetDirectory(0);
	forceZ->SetDirectory(0);
	positionSamples->SetDirectory(0);

	// Check if histograms were loaded correctly
	if (!(forceX && forceY && forceZ && positionSamples)) 
	{
		std::cerr << "Error loading histograms!" << std::endl;
	}

	infile.Close();

	FillData();
}

void CoolingForceValue::CopyParameters()
{
	eBeamParameter = ElectronBeam::GetParameters();
	ionBeamParameter = IonBeam::GetParameters();
	labEnergiesParameter = LabEnergy::GetParameters();
}

void CoolingForceValue::SetupLabel()
{
	if (!eBeamParameter.densityFile.get().empty() && !labEnergiesParameter.energyFile.get().empty())
	{
		index = std::stoi(eBeamParameter.densityFile.get().filename().string().substr(0, 4));
	}

	label = Form("%d: U drift = %.2fV, v_d = %.1f", index, labEnergiesParameter.driftTubeVoltage.get(),
		eBeamParameter.detuningVelocity.get());
}

void CoolingForceValue::SetupTags()
{
	tags += ElectronBeam::GetTags();
	tags += LabEnergy::GetTags();
	tags += IonBeam::GetTags();
}

void CoolingForceValue::ResetDefaultValues()
{
	forceX = nullptr;
	forceY = nullptr;
	forceZ = nullptr;

	positionSamples = nullptr;

	forceXIntegral = 0.0;
	forceYIntegral = 0.0;
	forceZIntegral = 0.0;

	forceXValue = 0.0;
	forceYValue = 0.0;
	forceZValue = 0.0;

	label = "";
	tags = "";
	index = 0;
}

std::string CoolingForceValue::GetHeaderString() const
{
	std::string string =
		eBeamParameter.toString() +
		labEnergiesParameter.toString() +
		ionBeamParameter.toString(false);

	return string;
}

std::string CoolingForceValue::Filename() const
{
	std::ostringstream indexSS;
	indexSS << std::setw(4) << std::setfill('0') << index;

	std::string string = indexSS.str() + std::string(Form(" v_d %.0fmps", eBeamParameter.detuningVelocity.get()));

	return string;
}
