#include "pch.h"

#include "ElectronBeam.h"
#include "LabEnergies.h"
#include "Constants.h"
#include "FileUtils.h"
#include "ROOTcanvas.h"
#include "HeatMapData.h"
#include "HistUtils.h"

namespace ElectronBeam
{
	static ElectronBeamParameters parameter;

	// 3D Hist with main data
	static TH3D* beam = nullptr;

	// optional parameters
	static bool gaussianElectronBeam = false;
	static bool cylindricalElectronBeam = false;
	static bool noElectronBeamBend = false;
	static double electronBeamRadius = 0.05;
	static double electronBeamDensity = 1;
	static bool fixedLongitudinalTemperature = false;

	// parameters to increase histogram resolution by interpolation
	static bool increaseHist = false;
	static int factor = 3;

	static bool mirrorAroundZ = true;
	static bool cutOutZeros = true;

	// plotting data
	static ROOTCanvas* canvas = nullptr;

	static std::vector<PlotBeamData> plotBeams;
	static int selectedIndex = -1;

	static float SliceZ = 0.0f;

	// Trajectory values
	static float sliderZ = 0;
	static float sliderY = 0;

	void Init()
	{
		canvas = new ROOTCanvas("electron beam", "electron beam", 1200, 500);
		canvas->Divide(2, 1);
		PlotTrajectory();

		CalculateDetuningEnergy();
	}

	TH3D* Get()
	{
		return beam;
	}

	void SetupDistribution(std::filesystem::path densityfile)
	{
		delete beam;
		if (gaussianElectronBeam || cylindricalElectronBeam)
		{
			beam = GenerateElectronBeamDensity();
		}
		else
		{
			beam = LoadDensityFile(densityfile);
		}
	}

	void CalculateDetuningEnergy()
	{
		parameter.detuningEnergy = pow(sqrt(LabEnergy::GetCenterLabEnergy()) - sqrt(parameter.coolingEnergy), 2);
	}

	void CalculateDetuningVelocity()
	{
		double electronVel = TMath::Sqrt(2 * LabEnergy::GetCenterLabEnergy() * TMath::Qe() / PhysicalConstants::electronMass);
		double ionVel = TMath::Sqrt(2 * parameter.coolingEnergy * TMath::Qe() / PhysicalConstants::electronMass);
		parameter.detuningVelocity = ionVel - electronVel;
	}

	void CalculateEstimateLongkT()
	{
		parameter.longitudinal_kT_estimate = GetLongitudinal_kT(LabEnergy::GetCenterLabEnergy());
	}

	double GetVelocityMagnitude(double energy)
	{
		return TMath::Sqrt(2 * energy * TMath::Qe() / PhysicalConstants::electronMass);
	}

	double GetLongitudinal_kT(double labEnergy)
	{
		if (fixedLongitudinalTemperature)
		{
			return parameter.longitudinal_kT_estimate;
		}

		// intermediate unit is [J] final unit is [eV]
		double A = (1. + TMath::Power((parameter.expansionFactor - 1.) / parameter.expansionFactor, 2.))
			* (TMath::Power(TMath::K() * parameter.cathodeTemperature, 2.)) / (2. * TMath::Qe() * labEnergy);
		double B = 2.544008e-27;

		double C = parameter.LLR * TMath::Power(parameter.cathodeRadius, -2 / 3.) * TMath::Power(parameter.electronCurrent, 1 / 3.)
			* TMath::Power(parameter.extractionEnergy * TMath::Qe(), -1 / 6.) * parameter.extractionEnergy / labEnergy;

		double D = TMath::Power(parameter.sigmaLabEnergy * TMath::Qe(), 2.) / (2. * TMath::Qe() * labEnergy);
		return (A + B * C + D) / TMath::Qe();
	}

	double GetTransverse_kT()
	{
		return parameter.transverse_kT;
	}

	double GetCoolingEnergy()
	{
		return parameter.coolingEnergy;
	}

	double GetDensity(const Point3D& point)
	{
		return HistUtils::GetValueAtPosition(beam, point);
	}

	void SetElectronCurrent(double current)
	{
		parameter.electronCurrent = current;
	}

	std::string GetTags()
	{
		std::string tags = "";
		if (gaussianElectronBeam) tags += "e-gaus, ";
		if (cylindricalElectronBeam) tags += "e-cylinder, ";
		if (noElectronBeamBend) tags += "no bend, ";
		if (fixedLongitudinalTemperature) tags += "fixed kT||, ";

		return tags;
	}

	ElectronBeamParameters GetParameters()
	{
		return parameter;
	}

	std::string GetSelectedBeamLabel()
	{
		if (selectedIndex < 0) return nullptr;

		return plotBeams.at(selectedIndex).GetLabel();
	}

	TH3D* LoadDensityFile(std::filesystem::path file)
	{
		if (!file.empty())
		{
			TH3D* result = FileUtils::LoadMatrixFile(file);
			if(cutOutZeros) 
				result = CutZerosFromDistribution(result);

			if (mirrorAroundZ)
				result = MirrorDistributionAtZ(result);

			if (increaseHist)
				result = CreateLargeDistribution(result);
			
			result->SetTitle("electron distribution");
			result->SetName("electron distribution");
			parameter.densityFile.set(file);

			return result;
		}

		return nullptr;
	}

	PlotBeamData* GetSelected()
	{
		if (selectedIndex >= 0)
			return &plotBeams.at(selectedIndex);
		else
			return nullptr;
	}

	void SelectedItemChanged()
	{
		PlotBeamData& newlySelected = plotBeams.at(selectedIndex);
		newlySelected.UpdateSlice(SliceZ);

		if (canvas->IsShown()) newlySelected.Plot3D(canvas, 2);
	}

	void AddBeamToList(PlotBeamData& eBeam)
	{
		plotBeams.push_back(std::move(eBeam));
		if (plotBeams.size() == 1)
		{
			selectedIndex = 0;
			SelectedItemChanged();
		}
	}

	void RemoveBeamFromList(int index)
	{
		plotBeams.erase(plotBeams.begin() + index);
		selectedIndex = std::min(selectedIndex, (int)plotBeams.size() - 1);
		if (selectedIndex >= 0)
		{
			SelectedItemChanged();
		}
	}

	TVector3 GetDirection(double z)
	{
		if (noElectronBeamBend)
			return TVector3(0, 0, 1);

		double derivative = Derivative(z);
		TVector3 direction(0, derivative, 1);
		return direction.Unit();
	}

	TVector3 GetVelocity(double z, double energy)
	{
		return GetVelocityMagnitude(energy) * GetDirection(z);
	}

	TVector3 GetNormal(double z)
	{
		return GetDirection(z).Orthogonal();
	}

	TH3D* CutZerosFromDistribution(TH3D* input)
	{
		// Identify the non-zero bin ranges
		int minX = input->GetNbinsX(), maxX = 0;
		int minY = input->GetNbinsY(), maxY = 0;
		int minZ = 1;
		int maxZ = input->GetNbinsZ();

		// Loop over all bins to find the first and last non-zero bins in both x and y
		for (int i = 1; i <= input->GetNbinsX(); i++)
		{
			for (int j = 1; j <= input->GetNbinsY(); j++)
			{
				for (int k = 1; k <= input->GetNbinsZ(); k++)
				{
					if (input->GetBinContent(i, j, k) != 0)
					{
						if (i < minX) minX = i;
						if (i > maxX) maxX = i;
						if (j < minY) minY = j;
						if (j > maxY) maxY = j;
					}
				}
			}
		}

		// add/subtract 1 to keep one row of zeros in each direction for interpolation purposes
		minX -= 1;
		maxX += 1;
		minY -= 1;
		maxY += 1;

		int newBinsX = maxX - minX + 1;
		int newBinsY = maxY - minY + 1;
		int newBinsZ = input->GetNbinsZ();

		double xMin = input->GetXaxis()->GetBinLowEdge(minX);
		double xMax = input->GetXaxis()->GetBinUpEdge(maxX);
		double yMin = input->GetYaxis()->GetBinLowEdge(minY);
		double yMax = input->GetYaxis()->GetBinUpEdge(maxY);
		double zMin = input->GetZaxis()->GetBinLowEdge(minZ);
		double zMax = input->GetZaxis()->GetBinUpEdge(maxZ);

		std::string newName = std::string(input->GetName()) + "_cut";
		TH3D* output = new TH3D(newName.c_str(), input->GetTitle(),
			newBinsX, xMin, xMax,
			newBinsY, yMin, yMax,
			newBinsZ, zMin, zMax);

		output->GetXaxis()->SetTitle(input->GetXaxis()->GetTitle());
		output->GetYaxis()->SetTitle(input->GetYaxis()->GetTitle());
		output->GetZaxis()->SetTitle(input->GetZaxis()->GetTitle());

		// Step 4: Copy the relevant bin contents to the new histogram
		for (int i = minX; i <= maxX; i++)
		{
			for (int j = minY; j <= maxY; j++)
			{
				for (int k = minZ; k <= maxZ; k++)
				{
					double content = input->GetBinContent(i, j, k);
					output->SetBinContent(i - minX + 1, j - minY + 1, k, content); // Copy with adjusted indices
				}
			}
		}

		delete input;
		return output;
	}

	TH3D* MirrorDistributionAtZ(TH3D* input)
	{
		// Get bin counts
		int nBinsX = input->GetNbinsX();
		int nBinsY = input->GetNbinsY();
		int nBinsZ = input->GetNbinsZ();

		// Get axis ranges
		double xMin = input->GetXaxis()->GetXmin();
		double xMax = input->GetXaxis()->GetXmax();
		double yMin = input->GetYaxis()->GetXmin();
		double yMax = input->GetYaxis()->GetXmax();
		double zMin = -input->GetZaxis()->GetXmax(); // Extend into negative Z
		double zMax = input->GetZaxis()->GetXmax();

		// Create mirrored histogram with double Z range
		std::string newName = std::string(input->GetName()) + "_mirrored";
		TH3D* mirrored = new TH3D(newName.c_str(), input->GetTitle(),
			nBinsX, xMin, xMax,
			nBinsY, yMin, yMax,
			2 * nBinsZ - 1, zMin, zMax);

		mirrored->GetXaxis()->SetTitle(input->GetXaxis()->GetTitle());
		mirrored->GetYaxis()->SetTitle(input->GetYaxis()->GetTitle());
		mirrored->GetZaxis()->SetTitle(input->GetZaxis()->GetTitle());

		// Copy original data and mirror it
		for (int i = 1; i <= nBinsX; i++)
		{
			for (int j = 1; j <= nBinsY; j++)
			{
				for (int k = 1; k <= nBinsZ; k++)
				{
					double binContent = input->GetBinContent(i, j, k);
					double z = input->GetZaxis()->GetBinCenter(k);

					// Get bin indices in mirrored histogram
					int k_neg = mirrored->GetZaxis()->FindBin(-z);
					int k_pos = mirrored->GetZaxis()->FindBin(z);

					// Fill both positive and negative Z bins
					mirrored->SetBinContent(i, j, k_pos, binContent);
					mirrored->SetBinContent(i, j, k_neg, binContent);
				}
			}
		}

		delete input;
		return mirrored;
	}

	TH3D* CreateLargeDistribution(TH3D* input)
	{
		int numberBinsX = input->GetNbinsX() * factor;
		int numberBinsY = input->GetNbinsY() * factor;
		int numberBinsZ = input->GetNbinsZ() * factor;

		double firstBinCenterX = input->GetXaxis()->GetBinCenter(1);
		double firstBinCenterY = input->GetYaxis()->GetBinCenter(1);
		double firstBinCenterZ = input->GetZaxis()->GetBinCenter(1);

		double lastBinCenterX = input->GetXaxis()->GetBinCenter(input->GetNbinsX());
		double lastBinCenterY = input->GetYaxis()->GetBinCenter(input->GetNbinsY());
		double lastBinCenterZ = input->GetZaxis()->GetBinCenter(input->GetNbinsZ());

		std::string newName = std::string(input->GetName()) + "_large";
		TH3D* output = new TH3D(newName.c_str(), input->GetTitle(),
			numberBinsX, input->GetXaxis()->GetXmin(), input->GetXaxis()->GetXmax(),
			numberBinsY, input->GetYaxis()->GetXmin(), input->GetYaxis()->GetXmax(),
			numberBinsZ, input->GetZaxis()->GetXmin(), input->GetZaxis()->GetXmax());

		for (int i = 1; i <= numberBinsX; i++)
		{
			for (int j = 1; j <= numberBinsY; j++)
			{
				for (int k = 1; k <= numberBinsZ; k++)
				{
					double x = output->GetXaxis()->GetBinCenter(i);
					double y = output->GetYaxis()->GetBinCenter(j);
					double z = output->GetZaxis()->GetBinCenter(k);

					double x_modified = std::min(std::max(x, firstBinCenterX), lastBinCenterX - 1e-5);
					double y_modified = std::min(std::max(y, firstBinCenterY), lastBinCenterY - 1e-5);
					double z_modified = std::min(std::max(z, firstBinCenterZ), lastBinCenterZ - 1e-5);

					//std::cout << x << ", " << y << ", " << z << "\n";
					double value = input->Interpolate(x_modified, y_modified, z_modified);
					output->SetBinContent(i, j, k, value);
				}
			}
		}
		delete input;
		return output;
	}

	TH3D* GenerateElectronBeamDensity()
	{
		int nXBins = 100;
		int nYBins = 100;
		int nZBins = 200;

		double xmin = -0.08;
		double xmax = 0.08;
		double ymin = -0.08;
		double ymax = 0.08;
		double zmin = -0.7;
		double zmax = 0.7;

		TH3D* eBeam = new TH3D("generated densites", "generated densites", nXBins, xmin, xmax,
			nYBins, ymin, ymax,
			nZBins, zmin, zmax);

		for (int i = 1; i <= nXBins; i++)
		{
			for (int j = 1; j <= nYBins; j++)
			{
				for (int k = 1; k <= nZBins; k++)
				{
					// Calculate the coordinates for this bin
					double x = eBeam->GetXaxis()->GetBinCenter(i);
					double y = eBeam->GetYaxis()->GetBinCenter(j);
					double z = eBeam->GetZaxis()->GetBinCenter(k);

					double ymean = 0;
					double value = 0;
					if (!noElectronBeamBend)
					{
						ymean = Trajectory(z);
					}
					if (gaussianElectronBeam)
					{
						// Calculate the value using the Gaussian distribution centered at z = 0
						value = electronBeamDensity * exp(-(x * x + (y - ymean) * (y - ymean)) / (2.0 * pow(electronBeamRadius, 2)));
					}
					else if (cylindricalElectronBeam)
					{
						if (x * x + (y - ymean) * (y - ymean) <= pow(electronBeamRadius, 2))
						{
							// if inside the cylinder, set the value to an arbitrary constant value
							value = electronBeamDensity;
						}
					}
					eBeam->SetBinContent(i, j, k, value);
				}
			}
		}
		TH3D* newBeam = CutZerosFromDistribution(eBeam);
		newBeam->SetTitle("generated electron distribution");
		newBeam->SetName("generated electron distribution");
		return newBeam;
	}

	void PlotTrajectory()
	{
		canvas->cd(1);

		const int N = 1000;
		double zValues[N];
		double yValues[N];

		double z_min = -0.7;
		double z_max = 0.7;

		for (int i = 0; i < N; i++)
		{
			double z = z_min + i * (z_max - z_min) / (N - 1);
			zValues[i] = z;
			yValues[i] = Trajectory(z);
		}

		TGraph* graph = new TGraph(N, zValues, yValues);
		graph->SetTitle("Trajectory; z; y(z)");
		graph->Draw("AL");

		// Select a point on the curve where you want to draw the arrow
		Point3D selectedPoint(0, sliderY, sliderZ);
		double z0 = selectedPoint.z;
		double y0 = selectedPoint.y; //Trajectory(z0); 

		// Define the arrow length
		double arrowLength = 0.2;

		TVector3 tangent = GetDirection(selectedPoint.z);
		//std::cout << "y: " << tangent.y() << " z: " << tangent.z() << "\n";
		tangent *= arrowLength;// / (1 + 2000 * y0);

		// Draw the arrow pointing in the direction of the trajectory
		TArrow* tangentArrow = new TArrow(z0, y0, z0 + tangent.z(), y0 + tangent.y(), 0.02f, "|>");
		tangentArrow->SetLineColor(kRed);
		tangentArrow->SetLineWidth(2);
		tangentArrow->SetAngle(30); // Arrowhead angle
		tangentArrow->Draw();

		TVector3 normal = GetNormal(z0);
		normal *= arrowLength / 20;

		// Draw the normal vector arrow
		TArrow* normalArrow = new TArrow(z0, y0, z0 + normal.z(), y0 + normal.y(), 0.02f, "|>");
		normalArrow->SetLineColor(kBlue);
		normalArrow->SetLineWidth(2);
		normalArrow->SetAngle(30); // Arrowhead angle
		normalArrow->Draw();
	}

	// returns the y value as function of z
	double Trajectory(double z)
	{
		z = TMath::Abs(z);
		return pow(10, (-7.374 + 0.9 * z + 19.8 * z * z + 2.055 * z * z * z - 20.56 * z * z * z * z)) - pow(10, -7.374) + 1.0e-5;
	}

	double Derivative(double z)
	{
		double sign = TMath::Sign(1, z);
		z = TMath::Abs(z);
		return sign * pow(10, (0.9 * z + 19.8 * z * z + 2.055 * z * z * z - 20.56 * z * z * z * z - 7.374)) * (39.6 * z + 6.165 * z * z - 82.24 * z * z * z + 0.9) * log(10);
	}

	void ShowWindow()
	{
		if (ImGui::Begin("Electron Beam"))
		{
			if (ImGui::BeginChild("left side", ImVec2(300, -1), ImGuiChildFlags_ResizeX))
			{
				ShowList();

				if (ImGui::Button("Load e-density file"))
				{
					std::vector<std::filesystem::path> files = FileUtils::SelectFiles();
					for (const std::filesystem::path& file : files)
					{
							TH3D* hist = LoadDensityFile(file);
							std::cout << hist->GetName() << std::endl;
							PlotBeamData newBeam(hist);
							std::string index = file.filename().string().substr(0, 4);
							std::string label = file.parent_path().parent_path().filename().string() + ": index " + index;
							newBeam.SetLabel(label);

							AddBeamToList(newBeam);
					}
				}
				ImGui::SameLine();
				ImGui::BeginDisabled(!(cylindricalElectronBeam || gaussianElectronBeam));
				if (ImGui::Button("generate density"))
				{
					TH3D* hist = GenerateElectronBeamDensity();

					PlotBeamData newBeam(hist);
					std::string shape = gaussianElectronBeam ? "gaussian" : "cylindrical";
					std::string bend = noElectronBeamBend ? "no bend" : "";
					std::string radius = std::to_string(electronBeamRadius);
					newBeam.SetLabel(shape + " beam, radius " + radius + " " + bend);
					AddBeamToList(newBeam);
				}
				ImGui::EndDisabled();
				ImGui::SameLine();
				if (ImGui::Button("clear list"))
				{
					for (int i = plotBeams.size() - 1; i >= 0; i--)
					{
						RemoveBeamFromList(i);
					}
				}

				ImGui::SetNextItemWidth(200.0f);
				if (ImGui::SliderFloat("slice z", &SliceZ, -0.7f, 0.7f))
				{
					if (selectedIndex >= 0)
					{
						PlotBeamData& eBeam = plotBeams.at(selectedIndex);
						eBeam.UpdateSlice(SliceZ);
					}
				}

				ImGui::Separator();
				ShowParameterControls();
				ImGui::Separator();
				canvas->MakeShowHideButton();
				ImGui::SameLine();
				PlotBeamData::ShowRebinningFactorsInput();

				ImGui::SeparatorText("arrow sliders");
				ImGui::SetNextItemWidth(220.0f);
				if (ImGui::SliderFloat("z", &sliderZ, -0.7f, 0.7f))
				{
					PlotTrajectory();
				}
				ImGui::SetNextItemWidth(220.0f);
				if (ImGui::SliderFloat("y", &sliderY, 0, 0.05f, "%.4f"))
				{
					PlotTrajectory();
				}
			}
			ImGui::EndChild();

			ImGui::SameLine();

			ShowPlots();

		}
		ImGui::End();
		canvas->Render();
	}

	void ShowList()
	{
		if (ImGui::BeginListBox("##electron beam list", ImVec2(-1, 270.0f)))
		{
			for (int i = 0; i < plotBeams.size(); i++)
			{
				ImGui::PushID(i);
				PlotBeamData& eBeam = plotBeams.at(i);

				if (ImGui::Selectable(eBeam.GetLabel().c_str(), i == selectedIndex, ImGuiSelectableFlags_AllowItemOverlap))
				{
					selectedIndex = i;
					SelectedItemChanged();
				}

				ImGui::SameLine();
				if (ImGui::SmallButton("x"))
				{
					RemoveBeamFromList(i);
				}
				ImGui::PopID();
			}
			ImGui::EndListBox();
		}
	}

	void ShowParameterControls()
	{
		ImGui::BeginGroup();

		ImGui::PushItemWidth(90.0f);

		ImGui::Checkbox("use fixed longitudinal kT", &fixedLongitudinalTemperature);
		ImGui::BeginDisabled(!fixedLongitudinalTemperature);
		ImGui::InputDouble("longitudinal kT [eV]", parameter.longitudinal_kT_estimate);
		ImGui::EndDisabled();
		ImGui::InputDouble("cooling energy [eV]", parameter.coolingEnergy);
		ImGui::InputDouble("transverse kT [eV]", parameter.transverse_kT);
		ImGui::BeginDisabled(fixedLongitudinalTemperature);
		ImGui::InputDouble("cathode radius [m]", parameter.cathodeRadius);
		ImGui::InputDouble("cathode Temperature [K]", parameter.cathodeTemperature);
		ImGui::InputDouble("extraction energy [eV]", parameter.extractionEnergy);
		ImGui::InputDouble("expansion factor", parameter.expansionFactor);
		ImGui::InputDouble("LLR", parameter.LLR);
		ImGui::InputDouble("sigma lab energy [eV]", parameter.sigmaLabEnergy);
		ImGui::EndDisabled();

		ImGui::SeparatorText("Loading options");
		ImGui::Checkbox("increase bin number", &increaseHist);
		ImGui::SameLine();
		ImGui::BeginDisabled(!increaseHist);
		ImGui::InputInt("factor", &factor, 2);
		ImGui::EndDisabled();

		ImGui::Checkbox("mirror around z-axis", &mirrorAroundZ);
		ImGui::SameLine();
		ImGui::Checkbox("cut out zeros", &cutOutZeros);

		ImGui::SeparatorText("special beam shapes");
		if (ImGui::Checkbox("gaussian", &gaussianElectronBeam))
		{
			parameter.densityFile.get().clear();
			cylindricalElectronBeam = false;

			// remove bend if unticked
			if (!gaussianElectronBeam)
				noElectronBeamBend = false;
		}
		ImGui::SameLine();
		if (ImGui::Checkbox("cylindrical", &cylindricalElectronBeam))
		{
			parameter.densityFile.get().clear();
			gaussianElectronBeam = false;

			// remove bend if unticked
			if (!cylindricalElectronBeam)
				noElectronBeamBend = false;
		}
		ImGui::BeginDisabled(!(gaussianElectronBeam || cylindricalElectronBeam));
		ImGui::SameLine();
		ImGui::Checkbox("no bend", &noElectronBeamBend);
		ImGui::InputDouble("radius [m]", &electronBeamRadius);
		ImGui::InputDouble("density [1/m^3]", &electronBeamDensity, 0, 0, "%.2e");
		ImGui::EndDisabled();

		ImGui::PopItemWidth();
		ImGui::EndGroup();
	}

	void ShowPlots()
	{
		if (ImPlot::BeginSubplots("##eBeam subplots", 2, 3, ImVec2(-1, -1), ImPlotSubplotFlags_ShareItems))
		{
			if (ImPlot::BeginPlot("Projection X"))
			{
				ImPlot::SetupAxes("x", "normalised value");
				for (const PlotBeamData& eBeam : plotBeams)
				{
					eBeam.PlotProjectionX();
				}
				ImPlot::EndPlot();
			}

			if (ImPlot::BeginPlot("Projection Y"))
			{
				ImPlot::SetupAxes("y", "normalised value");
				for (const PlotBeamData& eBeam : plotBeams)
				{
					eBeam.PlotProjectionY();
				}
				ImPlot::EndPlot();
			}

			if (ImPlot::BeginPlot("Projection Z"))
			{
				ImPlot::SetupAxes("z", "normalised value");
				for (const PlotBeamData& eBeam : plotBeams)
				{
					eBeam.PlotProjectionZ();
				}
				ImPlot::EndPlot();
			}

			if (ImPlot::BeginPlot("Inside/Outside"))
			{
				ImPlot::SetupAxes("z", "value");
				for (const PlotBeamData& eBeam : plotBeams)
				{
					eBeam.PlotInsideOutsideValue();
				}
				ImPlot::EndPlot();
			}

			if (selectedIndex >= 0)
			{
				const PlotBeamData& sliceBeam = plotBeams.at(selectedIndex);
				sliceBeam.PlotSlice();
			}

			ImPlot::EndSubplots();
		}
	}
}

