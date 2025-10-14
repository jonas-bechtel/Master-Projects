#include "pch.h"
#include "CoolingForceValue.h"
#include "CoolingForceWindow.h"

#include "CoolingForceModel.h"

#include "Constants.h"
#include "FileUtils.h"
#include "HistUtils.h"
#include "Timer.h"

#include "Application.h"

namespace CoolingForce
{
	//RNG_engine Value::generator = RNG_engine();
	//
	//std::normal_distribution<double> Value::longitudinalNormalDistribution = std::normal_distribution<double>();
	//std::normal_distribution<double> Value::transverseNormalDistribution = std::normal_distribution<double>();


	Value::Value()
	{
		//std::cout << "calling cf value default constructor" << std::endl;
	}

	Value::~Value()
	{
		//std::cout << "calling cf value destructor" << std::endl;
		
	}

	Value::Value(Value&& other) noexcept
	{
		//std::cout << "calling cf value move constructor" << std::endl;
		force3D = std::move(other.force3D);

		forceValue = other.forceValue;

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

		force3D = std::move(other.force3D);
		
		forceValue = other.forceValue;

		eBeamParameter = other.eBeamParameter;
		ionBeamParameter = other.ionBeamParameter;
		labEnergiesParameter = other.labEnergiesParameter;

		label = std::move(other.label);
		tags = std::move(other.tags);
		index = std::move(other.index);

		other.ResetDefaultValues();

		return *this;
	}

	void Value::Calculate(std::filesystem::path descriptionFile, int index, Model::Parameter params)
	{
		PrepareCalculation(descriptionFile, index);
		CalculateForce3D(params);

		float sigmaX = IonBeam::GetSigmaX();
		float sigmaY = IonBeam::GetSigmaY();

		float xLow = -5.0 * sigmaX;
		float xHigh = 5.0 * sigmaX;
		float yLow = -5.0 * sigmaY;
		float yHigh = 5.0 * sigmaY;
		float zLow = -0.7;
		float zHigh = 0.7;
		if (IonBeam::IsRangeLimited())
		{
			zLow = IonBeam::GetLimitedRange()[0];
			zHigh = IonBeam::GetLimitedRange()[1];
		}

		TF3 func("integral func", this, &Value::Integrand, xLow, xHigh, yLow, yHigh, zLow, zHigh, 0, 3);

		forceValue = func.Integral(xLow, xHigh, yLow, yHigh, zLow, zHigh, 1.0e-7) / CSR::coolerLength;
	}

	double Value::Integrand(double* position, double* params)
	{
		double x = position[0];
		double y = position[1];
		double z = position[2];

		double forceValue = HistUtils::GetValueAtPosition(force3D.GetHist(), { x,y,z });
		return IonBeam::GetValue(x, y, z) * forceValue;
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

		if (ImGui::BeginDragDropSource(ImGuiDragDropFlags_None))
		{
			ImGui::SetDragDropPayload("Cooling Force Value", this, sizeof(Value));
			ImGui::Text("drag to set params");
			ImGui::EndDragDropSource();
		}

		if (!Application::GetSettings().tooltipsDisabled && ImGui::BeginItemTooltip())
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

	bool Value::ShowCalcTransForceCheckbox()
	{
		return ImGui::Checkbox("calculate transverse force", &calculateTransverseForce);;
	}

	void Value::SetupHistogramsFromReference(TH3D* reference)
	{
		// read binning from reference hist
		int nBinsX = reference->GetNbinsX();
		double xMin = reference->GetXaxis()->GetXmin();
		double xMax = reference->GetXaxis()->GetXmax();

		int nBinsY = reference->GetNbinsY();
		double yMin = reference->GetYaxis()->GetXmin();
		double yMax = reference->GetYaxis()->GetXmax();

		int nBinsZ = reference->GetNbinsZ();
		double zMin = reference->GetZaxis()->GetXmin();
		double zMax = reference->GetZaxis()->GetXmax();

		// create histograms for force components
		TH3D* forceLongHist = new TH3D("forceLong", "force longitudinal", nBinsX, xMin, xMax, nBinsY, yMin, yMax, nBinsZ, zMin, zMax);
		
		// fill 3D hist data structures
		force3D = HistData3D(forceLongHist);
	}

	void Value::PlotForceSlice() const
	{
		force3D.PlotSlice();
	}

	void Value::UpdateSlice(float zValue)
	{
		force3D.UpdateSlice(zValue);
	}

	ElectronBeamParameters Value::GetElectronBeamParameters() const
	{
		return eBeamParameter;
	}

	IonBeamParameters Value::GetIonBeamParameters() const
	{
		return ionBeamParameter;
	}

	LabEnergyParameters Value::GetLabEnergyParameters() const
	{
		return labEnergiesParameter;
	}

	void Value::Save(std::filesystem::path folder) const
	{
		std::filesystem::path file = folder / (Filename() + ".root");

		// Open a ROOT file for writing
		TFile outfile(file.string().c_str(), "RECREATE");
		if (!outfile.IsOpen())
		{
			std::cout << "error opening file: " << file << std::endl;
		}
		TNamed header("header", GetHeaderString());
		header.Write();

		std::ostringstream value_ss;
		value_ss << std::scientific << forceValue;

		TNamed value("value", value_ss.str());
		value.Write();

		// Write histograms into the file
		force3D.GetHist()->Write();

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

		TNamed* value = (TNamed*)infile.Get("value");
		forceValue = std::stod(value->GetTitle());

		// Retrieve histograms
		force3D = HistData3D((TH3D*)infile.Get("forceLong"));

		// Check if histograms were loaded correctly
		if (!(force3D.GetHist()))
		{
			std::cerr << "Error loading histograms!" << std::endl;
		}
		else
		{
			force3D.GetHist()->SetDirectory(0);
		}

		infile.Close();
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

	void Value::CalculateForce3D(Model::Parameter params)
	{
		TVector3 ionVelocity = IonBeam::GetVelocity();
		TH3D* hist = force3D.GetHist();
		//TH3D* histTrans = forceTrans.GetHist();

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

					double value = Model::Force(parameter, calculateTransverseForce);
					
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

						double value = Model::Force(params, calculateTransverseForce);

						hist->SetBinContent(i, j, k, value);
						int bin = hist->FindBin(x, y, -z);
						hist->SetBinContent(bin, value);
					}
				}
			}
		}

		force3D.UpdateData();
		force3D.UpdateSlice(GetSliceValue());
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

		forceValue = 0.0;

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

