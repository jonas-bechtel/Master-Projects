#include "pch.h"
#include "CoolingForceValue.h"
#include "CoolingForceWindow.h"

#include "CoolingForceModel.h"

#include "Constants.h"
#include "FileUtils.h"
#include "HistUtils.h"
#include "Timer.h"

namespace CoolingForce
{
	bool Value::parallelForcePrecalculation = true;
	RNG_engine Value::generator = RNG_engine();

	std::normal_distribution<double> Value::longitudinalNormalDistribution = std::normal_distribution<double>();
	std::normal_distribution<double> Value::transverseNormalDistribution = std::normal_distribution<double>();


	Value::Value()
	{
		//std::cout << "calling cf value default constructor" << std::endl;
	}

	Value::~Value()
	{
		//std::cout << "calling cf value destructor" << std::endl;
		delete forceX;
		delete forceY;
		delete forceZ;
		delete positionSamples;
	}

	Value::Value(Value&& other) noexcept
	{
		//std::cout << "calling cf value move constructor" << std::endl;
		precalculatedForce = std::move(other.precalculatedForce);

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

	Value& Value::operator=(Value&& other) noexcept
	{
		//std::cout << "calling cf value move assignment" << std::endl;
		if (this == &other) return *this;

		precalculatedForce = std::move(other.precalculatedForce);

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

	void Value::CalculateOriginal(std::filesystem::path descriptionFile, int index, Model::Parameter params)
	{
		PrepareCalculation(descriptionFile, index);

		// calculate cooling force
		std::vector<Point3D> samples = IonBeam::GeneratePositions();

		if (samples.empty())
		{
			std::cout << "sampling positions failed\n";
			return;
		}

		for (const Point3D& point : samples)
		{
			double x = point.x;
			double y = point.y;
			double z = point.z;

			// calculate velocity magnitude from lab energy given as matrix (TH3D) from outside
			double labEnergy = LabEnergy::GetValue(x, y, z);
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

			Model::Parameter params;
			params.relativeVelocity = collisionVelocity;
			params.electronDensity = ElectronBeam::GetDensity(point);
			params.ionCharge = IonBeam::GetCharge();
			params.kT_long = long_kT;
			params.kT_trans = trans_kT;
			TVector3 coolingforce = Model::CoolingForce(params);

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

	void Value::CalculateHalfIntegrated(std::filesystem::path descriptionFile, int index, Model::Parameter params, bool interpolate)
	{
		Timer t;
		PrepareCalculation(descriptionFile, index);
		std::cout << "prep time: " << t.ElapsedMillis() << std::endl;
		t.Reset();
		PrecalculateForce(params);
		std::cout << "precalc time: " << t.ElapsedMillis() << std::endl;
		t.Reset();
		std::vector<Point3D> samples = IonBeam::GeneratePositions();
		std::cout << "generate samples time: " << t.ElapsedMillis() << std::endl;
		t.Reset();

		if (samples.empty())
		{
			std::cout << "sampling positions failed\n";
			return;
		}

		for (const Point3D& point : samples)
		{
			double x = point.x;
			double y = point.y;
			double z = point.z;

			TVector3 coolingforce;
			TH3D* preForceHist = precalculatedForce.GetHist();
			double xMin = preForceHist->GetXaxis()->GetBinCenter(1);
			double xMax = preForceHist->GetXaxis()->GetBinCenter(preForceHist->GetNbinsX());

			double yMin = preForceHist->GetYaxis()->GetBinCenter(1);
			double yMax = preForceHist->GetYaxis()->GetBinCenter(preForceHist->GetNbinsY());

			double zMin = preForceHist->GetZaxis()->GetBinCenter(1);
			double zMax = preForceHist->GetZaxis()->GetBinCenter(preForceHist->GetNbinsZ());

			if (x >= xMin && x <= xMax && y >= yMin && y <= yMax && abs(z) >= zMin && abs(z) <= zMax)
			{
				if (interpolate)
					coolingforce = { 0, 0, preForceHist->Interpolate(x, y, abs(z)) };
				else
					coolingforce = { 0, 0, preForceHist->GetBinContent(preForceHist->FindBin(x, y, abs(z))) };
			}
			else
			{
				coolingforce = { 0,0,0 };
			}

			positionSamples->Fill(x, y, z);

			int bin = forceX->FindBin(x, y, z);
			forceX->AddBinContent(bin, coolingforce.x());
			forceY->AddBinContent(bin, coolingforce.y());
			forceZ->AddBinContent(bin, coolingforce.z());
		}
		std::cout << "looping samples time: " << t.ElapsedMillis() << std::endl;
		t.Reset();

		forceX->Divide(positionSamples);
		forceY->Divide(positionSamples);
		forceZ->Divide(positionSamples);

		FillData();
		std::cout << "dividing and finish time: " << t.ElapsedMillis() << std::endl;
	}

	void Value::CalculateFullIntegrated(std::filesystem::path descriptionFile, int index, Model::Parameter params)
	{
		Timer t;
		PrepareCalculation(descriptionFile, index);
		std::cout << "prep time: " << t.ElapsedMillis() << std::endl;
		t.Reset();
		PrecalculateForce(params);
		std::cout << "precalc time: " << t.ElapsedMillis() << std::endl;
		t.Reset();
		TH3D* intermediate = (TH3D*)precalculatedForce.GetHist()->Clone("force times ion density");
		intermediate->Multiply(IonBeam::Get());

		forceZValue = intermediate->Integral("width") / 0.8;
		std::cout << "mult and integral time: " << t.ElapsedMillis() << std::endl;
	}

	void Value::CalculateFullIntegratedBetter(std::filesystem::path descriptionFile, int index, Model::Parameter params)
	{
		PrepareCalculation(descriptionFile, index);
		PrecalculateForce(params);

		float sigmaX = IonBeam::GetSigmaX();
		float sigmaY = IonBeam::GetSigmaY();

		float xLow = -5.0 * sigmaX;
		float xHigh = 5.0 * sigmaX;
		float yLow = -5.0 * sigmaY;
		float yHigh = 5.0 * sigmaY;
		float zLow = -0.7;
		float zHigh = 0.7;

		TF3 func("integral func", this, &Value::Integrand, xLow, xHigh, yLow, yHigh, zLow, zHigh, 0, 3);

		forceZValue = func.Integral(xLow, xHigh, yLow, yHigh, zLow, zHigh, 1.0e-7) / 0.8;
	}

	double Value::Integrand(double* position, double* params)
	{
		double x = position[0];
		double y = position[1];
		double z = position[2];

		double forceValue = HistUtils::GetValueAtPosition(precalculatedForce.GetHist(), { x,y,z });
		return IonBeam::GetValue({ x,y,z }) * forceValue;
	}

	bool Value::ShowListItem(bool selected) const
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

	bool Value::ShowParallelPrecalculationCheckbox()
	{
		return ImGui::Checkbox("parallel force calc", &parallelForcePrecalculation);
	}

	void Value::SetupHistogramsFromReference(TH3D* reference)
	{
		delete forceX;
		delete forceY;
		delete forceZ;
		delete positionSamples;
		//delete precalculatedForce;

		forceX = (TH3D*)reference->Clone("cooling force X");
		forceY = (TH3D*)reference->Clone("cooling force Y");
		forceZ = (TH3D*)reference->Clone("cooling force Z");
		positionSamples = (TH3D*)reference->Clone("position samples");
		TH3D* precalculatedForceHist = (TH3D*)reference->Clone("precalculated force");

		forceX->Reset();
		forceY->Reset();
		forceZ->Reset();
		positionSamples->Reset();
		precalculatedForceHist->Reset();

		forceX->SetTitle("cooling force X");
		forceY->SetTitle("cooling force Y");
		forceZ->SetTitle("cooling force Z");
		positionSamples->SetTitle("position samples");
		precalculatedForceHist->SetTitle("precalculated force");

		precalculatedForce = PlotBeamData(precalculatedForceHist);
	}

	void Value::FillData()
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

	double Value::CalculateIntegral(TH3D* hist)
	{
		int nBinsZ = hist->GetZaxis()->GetNbins();
		double zMin = hist->GetZaxis()->GetXmin();
		double zMax = hist->GetZaxis()->GetXmax();

		// Create a 1D histogram for weighted averages along Z
		TH1D* H_weightedAvgZ = new TH1D("H_weightedAvgZ", "Weighted Average per Z", nBinsZ, zMin, zMax);

		// Loop over Z bins
		for (int iZ = 0; iZ <= nBinsZ + 1; iZ++)
		{
			double weightedSum = 0.0;
			double weightSum = 0.0;

			// Loop over all (x, y) bins for this fixed Z bin
			for (int iX = 0; iX <= hist->GetXaxis()->GetNbins() + 1; iX++)
			{
				for (int iY = 0; iY <= hist->GetYaxis()->GetNbins() + 1; iY++)
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

	void Value::PlotPreForceSlize() const
	{
		precalculatedForce.PlotSlice();
	}

	void Value::UpdateSlice(float zValue)
	{
		precalculatedForce.UpdateSlice(zValue);
	}

	void Value::Save(std::filesystem::path folder) const
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

	void Value::Load(std::filesystem::path file)
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

	void Value::PrepareCalculation(std::filesystem::path descriptionFile, int index)
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

		SetupHistogramsFromReference(IonBeam::Get());
	}

	void Value::PrecalculateForce(Model::Parameter params)
	{
		TVector3 ionVelocity = IonBeam::GetVelocity();
		TH3D* hist = precalculatedForce.GetHist();

		int nBinsX = hist->GetXaxis()->GetNbins();
		int nBinsY = hist->GetYaxis()->GetNbins();
		int nBinsZ = hist->GetZaxis()->GetNbins();

		int totalBins = nBinsX * nBinsY * nBinsZ;
		std::vector<int> indices(totalBins);
		std::iota(indices.begin(), indices.end(), 0);

		if (parallelForcePrecalculation)
		{
			std::for_each_n(std::execution::par, indices.begin(), totalBins,
				[&](int idx)
				{
					int i = idx / (nBinsY * nBinsZ) + 1;
					int j = (idx / nBinsZ) % nBinsY + 1;
					int k = idx % nBinsZ + 1;

					double x = hist->GetXaxis()->GetBinCenter(i);
					double y = hist->GetYaxis()->GetBinCenter(j);
					double z = hist->GetZaxis()->GetBinCenter(k);
					if (z < 0) return;

					double labEnergy = LabEnergy::GetValue(x, y, z);
					TVector3 electronVelocity = ElectronBeam::GetVelocity(z, labEnergy);
					TVector3 relativeVelocity = ionVelocity - electronVelocity;

					Model::Parameter parameter = params;
					parameter.relativeVelocity = relativeVelocity;
					parameter.electronDensity = ElectronBeam::GetDensity({ x,y,z });
					parameter.kT_long = ElectronBeam::GetLongitudinal_kT(labEnergy);
					parameter.kT_trans = ElectronBeam::GetTransverse_kT();
					//std::cout << "outside: " << params.String() << std::endl;
					double value = Model::ForceZ(parameter);
					hist->SetBinContent(i, j, k, value);
					int bin = hist->FindBin(x, y, -z);
					hist->SetBinContent(bin, value);
				});
		}
		else
		{
			for (int i = 1; i <= hist->GetXaxis()->GetNbins(); i++)
			{
				for (int j = 1; j <= hist->GetYaxis()->GetNbins(); j++)
				{
					for (int k = 1; k <= hist->GetZaxis()->GetNbins(); k++)
					{
						double x = hist->GetXaxis()->GetBinCenter(i);
						double y = hist->GetYaxis()->GetBinCenter(j);
						double z = hist->GetZaxis()->GetBinCenter(k);
						if (z < 0) continue;

						double labEnergy = LabEnergy::GetValue(x, y, z);
						double electronVelocityMagnitude = TMath::Sqrt(2 * labEnergy * TMath::Qe() / PhysicalConstants::electronMass);
						TVector3 longitudinalDirection = ElectronBeam::GetDirection(z);
						TVector3 electronVelocity = electronVelocityMagnitude * longitudinalDirection;
						TVector3 relativeVelocity = ionVelocity - electronVelocity;

						params.relativeVelocity = relativeVelocity;
						params.electronDensity = ElectronBeam::GetDensity({ x,y,z });
						params.kT_long = ElectronBeam::GetLongitudinal_kT(labEnergy);
						params.kT_trans = ElectronBeam::GetTransverse_kT();

						double value = Model::ForceZ_Original(params);
						hist->SetBinContent(i, j, k, value);
						int bin = hist->FindBin(x, y, -z);
						hist->SetBinContent(bin, value);
					}
				}
			}
		}

		precalculatedForce.UpdateData();
		precalculatedForce.UpdateSlice(GetSliceValue());
	}

	void Value::CopyParameters()
	{
		eBeamParameter = ElectronBeam::GetParameters();
		ionBeamParameter = IonBeam::GetParameters();
		labEnergiesParameter = LabEnergy::GetParameters();
	}

	void Value::SetupLabel()
	{
		if (!eBeamParameter.densityFile.get().empty() && !labEnergiesParameter.energyFile.get().empty())
		{
			index = std::stoi(eBeamParameter.densityFile.get().filename().string().substr(0, 4));
		}

		label = Form("%d: U drift = %.2fV, v_d = %.1f", index, labEnergiesParameter.driftTubeVoltage.get(),
			eBeamParameter.detuningVelocity.get());
	}

	void Value::SetupTags()
	{
		tags += ElectronBeam::GetTags();
		tags += LabEnergy::GetTags();
		tags += IonBeam::GetTags();
	}

	void Value::ResetDefaultValues()
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

	std::string Value::GetHeaderString() const
	{
		std::string string =
			eBeamParameter.toString() +
			labEnergiesParameter.toString() +
			ionBeamParameter.toString(false);

		return string;
	}

	std::string Value::Filename() const
	{
		std::ostringstream indexSS;
		indexSS << std::setw(4) << std::setfill('0') << index;

		std::string string = indexSS.str() + std::string(Form(" v_d %.0fmps", eBeamParameter.detuningVelocity.get()));

		return string;
	}
}

