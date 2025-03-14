#include "pch.h"

#include "ElectronBeam.h"
#include "Constants.h"
#include "FileHandler.h"

ElectronBeamWindow::ElectronBeamWindow()
	: EnergyDistributionModule("Electron Beam", 1), m_parameters(activeDist.eBeamParameter)
{
	eBeam = this;
	PlotTrajectory();

	CalculateDetuningEnergy();
}

void ElectronBeamWindow::SetupDistribution(std::filesystem::path densityfile)
{
	if (activeDist.simplifyParams.gaussianElectronBeam || activeDist.simplifyParams.cylindricalElectronBeam)
	{
		delete m_distribution;
		m_distribution = GenerateElectronBeamDensity();
	}
	else
	{
		LoadDensityFile(densityfile);
	}
}

void ElectronBeamWindow::CalculateDetuningEnergy()
{
	m_parameters.detuningEnergy = pow(sqrt(activeDist.labEnergiesParameter.centerLabEnergy)
		- sqrt(m_parameters.coolingEnergy), 2);
}

void ElectronBeamWindow::CalculateDetuningVelocity()
{
	double electronVel = TMath::Sqrt(2 * activeDist.labEnergiesParameter.centerLabEnergy * TMath::Qe() / PhysicalConstants::electronMass);
	double ionVel = TMath::Sqrt(2 * activeDist.eBeamParameter.coolingEnergy * TMath::Qe() / PhysicalConstants::electronMass);
	m_parameters.detuningVelocity = ionVel - electronVel;
}

void ElectronBeamWindow::LoadDensityFile(std::filesystem::path file)
{
	if (!file.empty())
	{
		delete m_distribution;
		m_distribution = FileHandler::GetInstance().LoadMatrixFile(file);
		m_distribution = CutZerosFromDistribution(m_distribution);
		if (mirrorAroundZ)
			m_distribution = MirrorDistributionAtZ(m_distribution);

		m_distribution->SetTitle("electron distribution");
		m_distribution->SetName("electron distribution");

		if(increaseHist)
			CreateLargeDistribution();

		m_parameters.densityFile.set(file);
	}
}

void ElectronBeamWindow::LoadToLookAt(std::filesystem::path file)
{
	if (!file.empty())
	{
		ElectronBeam newBeam;

		newBeam.fullHistogram = FileHandler::GetInstance().LoadMatrixFile(file);
		newBeam.fullHistogram = CutZerosFromDistribution(newBeam.fullHistogram);
		if (mirrorAroundZ)
			newBeam.fullHistogram = MirrorDistributionAtZ(newBeam.fullHistogram);

		newBeam.fullHistogram->SetTitle("electron beam");
		newBeam.fullHistogram->SetName("electron beam");

		newBeam.FillData();
		newBeam.label = file.parent_path().parent_path().filename().string() + ": index " + std::stoi(file.filename().string().substr(0, 4));

		AddBeamToList(newBeam);
	}
}

void ElectronBeamWindow::SelectedItemChanged()
{
	ElectronBeam& newlySelected = eBeamsToLookAt.at(selectedIndex);
	newlySelected.slice.FromTH3D(newlySelected.fullHistogram, SliceZ);

	delete m_distribution;
	m_distribution = (TH3D*)newlySelected.fullHistogram->Clone("copied e beam");
	if (IsCanvasShown(m_mainCanvas)) PlotDistribution();
}

void ElectronBeamWindow::AddBeamToList(ElectronBeam& eBeam)
{
	eBeamsToLookAt.push_back(std::move(eBeam));
	if (eBeamsToLookAt.size() == 1)
	{
		selectedIndex = 0;
		SelectedItemChanged();
	}
}

void ElectronBeamWindow::RemoveBeamFromList(int index)
{
	eBeamsToLookAt.erase(eBeamsToLookAt.begin() + index);
	selectedIndex = std::min(selectedIndex, (int)eBeamsToLookAt.size() - 1);
	if (selectedIndex >= 0)
	{
		SelectedItemChanged();
	}
}

TVector3 ElectronBeamWindow::GetDirection(double z)
{
	if (activeDist.simplifyParams.noElectronBeamBend)
		return TVector3(0, 0, 1);

	double derivative = Derivative(z);
	TVector3 direction(0, derivative, 1);
	return direction.Unit();
}

TVector3 ElectronBeamWindow::GetDirection(Point3D point)
{
	ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

	// Define the function to minimize
	ROOT::Math::Functor f([&](const double* z) { return DistancePointToTrajectoryOfZ(z[0], point); }, 1);
	minimizer->SetFunction(f);

	// Set initial value and step size for z
	double initial_z = point.z;
	minimizer->SetVariable(0, "z", initial_z, 0.1);

	minimizer->Minimize();

	// Get the result
	double resultingZ = minimizer->X()[0];

	return  GetDirection(resultingZ);
}

TVector3 ElectronBeamWindow::GetNormal(double z)
{
	return GetDirection(z).Orthogonal();
}

double ElectronBeamWindow::GetLongitudinal_kT(double labEnergy)
{
	if (activeDist.simplifyParams.fixedLongitudinalTemperature)
	{
		return activeDist.eBeamParameter.longitudinal_kT_estimate;
	}

	// intermediate unit is [J] final unit is [eV]
	double A = (1. + TMath::Power((m_parameters.expansionFactor - 1.) / m_parameters.expansionFactor, 2.))
		* (TMath::Power(TMath::K() * m_parameters.cathodeTemperature, 2.)) / (2. * TMath::Qe() * labEnergy);
	double B = 2.544008e-27;
	double C = m_parameters.LLR * TMath::Power(m_parameters.cathodeRadius, -2 / 3.) * TMath::Power(m_parameters.electronCurrent, 1 / 3.)
		* TMath::Power(m_parameters.extractionEnergy * TMath::Qe(), -1 / 6.) * m_parameters.extractionEnergy / labEnergy;
	double D = TMath::Power(m_parameters.sigmaLabEnergy * TMath::Qe(), 2.) / (2. * TMath::Qe() * labEnergy);
	return (A + B * C + D) / TMath::Qe();
}

double ElectronBeamWindow::GetTransverse_kT()
{
	return m_parameters.transverse_kT;
}

double ElectronBeamWindow::GetDensity(Point3D point)
{
	int numberBinsX = m_distribution->GetXaxis()->GetNbins();
	int numberBinsY = m_distribution->GetYaxis()->GetNbins();
	int numberBinsZ = m_distribution->GetZaxis()->GetNbins();

	if (point.x < m_distribution->GetXaxis()->GetBinCenter(1) || point.x > (m_distribution->GetXaxis()->GetBinCenter(numberBinsX) - 1e-4) ||
		point.y < m_distribution->GetYaxis()->GetBinCenter(1) || point.y > (m_distribution->GetYaxis()->GetBinCenter(numberBinsY) - 1e-4) ||
		point.z < m_distribution->GetZaxis()->GetBinCenter(1) || point.z > (m_distribution->GetZaxis()->GetBinCenter(numberBinsZ) - 1e-4))
	{
		return 0;
	}

	return m_distribution->Interpolate(point.x, point.y, point.z);
}

ElectronBeam* ElectronBeamWindow::GetSelected()
{
	if (selectedIndex < 0) return nullptr;

	return &eBeamsToLookAt.at(selectedIndex);
}

void ElectronBeamWindow::ShowUI()
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

	ShowElectronBeamPlots();
}

void ElectronBeamWindow::ShowList()
{
	if (ImGui::BeginListBox("##electron beam list", ImVec2(-1, 270.0f)))
	{
		for (int i = 0; i < eBeamsToLookAt.size(); i++)
		{
			ImGui::PushID(i);
			ElectronBeam& eBeam = eBeamsToLookAt.at(i);

			if (ImGui::Selectable(eBeam.label.c_str(), i == selectedIndex, ImGuiSelectableFlags_AllowItemOverlap))
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

void ElectronBeamWindow::ShowSettings()
{
	if (ImGui::Button("Load e-density file"))
	{
		std::vector<std::filesystem::path> files = FileHandler::GetInstance().SelectFiles();
		for (const std::filesystem::path& file : files)
		{
			LoadToLookAt(file);
		}
	}
	ImGui::SameLine();
	if (ImGui::Button("clear list"))
	{
		for (int i = eBeamsToLookAt.size() - 1; i >= 0; i--)
		{
			RemoveBeamFromList(i);
		}
	}
	ImGui::SetNextItemWidth(200.0f);
	if (ImGui::SliderFloat("slice z", &SliceZ, -0.7f, 0.7f))
	{
		if (selectedIndex >= 0)
		{
			ElectronBeam& eBeam = eBeamsToLookAt.at(selectedIndex);
			eBeam.slice.FromTH3D(eBeam.fullHistogram, SliceZ);
		}
	}
	ImGui::Checkbox("multiply bins", &increaseHist);
	ImGui::SameLine();
	ImGui::BeginDisabled(!increaseHist);
	ImGui::InputInt("factor", &factor, 2);
	ImGui::EndDisabled();

	ImGui::Checkbox("mirror around z-axis", &mirrorAroundZ);

	ImGui::SeparatorText("special beam shapes");
	if (ImGui::Checkbox("gaussian", activeDist.simplifyParams.gaussianElectronBeam))
	{
		m_parameters.densityFile.get().clear();
		activeDist.simplifyParams.cylindricalElectronBeam = false;
	}
	ImGui::SameLine();
	if (ImGui::Checkbox("cylindrical", activeDist.simplifyParams.cylindricalElectronBeam))
	{
		m_parameters.densityFile.get().clear();
		activeDist.simplifyParams.gaussianElectronBeam = false;
	}
	ImGui::BeginDisabled(!(activeDist.simplifyParams.gaussianElectronBeam || activeDist.simplifyParams.cylindricalElectronBeam));
	ImGui::InputDouble("radius [m]", activeDist.simplifyParams.electronBeamRadius);
	ImGui::SameLine();
	ImGui::Checkbox("no bend", activeDist.simplifyParams.noElectronBeamBend);
	if (ImGui::Button("generate density"))
	{
		ElectronBeam newBeam;
		newBeam.fullHistogram = GenerateElectronBeamDensity();
		
		newBeam.FillData();
		std::string shape = activeDist.simplifyParams.gaussianElectronBeam ? "gaussian" : "cylindrical";
		std::string bend = activeDist.simplifyParams.noElectronBeamBend ? "no bend" : "";
		std::string radius = std::to_string(activeDist.simplifyParams.electronBeamRadius);
		newBeam.label = shape + " beam, radius " + radius + " " + bend;
		AddBeamToList(newBeam);
	}
	ImGui::EndDisabled();

	ImGui::SeparatorText("main parameter");
	ImGui::Checkbox("use fixed longitudinal kT", activeDist.simplifyParams.fixedLongitudinalTemperature);
	ImGui::BeginDisabled(!activeDist.simplifyParams.fixedLongitudinalTemperature);
	ImGui::InputDouble("longitudinal kT [eV]", activeDist.eBeamParameter.longitudinal_kT_estimate);
	ImGui::EndDisabled();
	ImGui::InputDouble("cooling energy [eV]", m_parameters.coolingEnergy);	
	ImGui::InputDouble("transverse kT [eV]", m_parameters.transverse_kT);
	ImGui::BeginDisabled(activeDist.simplifyParams.fixedLongitudinalTemperature);
	ImGui::InputDouble("cathode radius [m]", m_parameters.cathodeRadius);
	ImGui::InputDouble("cathode Temperature [K]", m_parameters.cathodeTemperature);
	ImGui::InputDouble("extraction energy [eV]", m_parameters.extractionEnergy);
	ImGui::InputDouble("expansion factor", m_parameters.expansionFactor);
	ImGui::InputDouble("LLR", m_parameters.LLR);
	ImGui::InputDouble("sigma lab energy [eV]", m_parameters.sigmaLabEnergy);
	ImGui::EndDisabled();

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

void ElectronBeamWindow::ShowElectronBeamPlots()
{
	if (ImPlot::BeginSubplots("##eBeam subplots", 2, 3, ImVec2(-1, -1), ImPlotSubplotFlags_ShareItems))
	{
		if (ImPlot::BeginPlot("Projection X"))
		{
			for (const ElectronBeam& eBeam : eBeamsToLookAt)
			{
				ImPlot::PlotLine(eBeam.label.c_str(), eBeam.xAxis.data(), eBeam.projectionValuesX.data(), eBeam.xAxis.size());
			}
			ImPlot::EndPlot();
		}

		if (ImPlot::BeginPlot("Projection Y"))
		{
			for (const ElectronBeam& eBeam : eBeamsToLookAt)
			{
				ImPlot::PlotLine(eBeam.label.c_str(), eBeam.yAxis.data(), eBeam.projectionValuesY.data(), eBeam.yAxis.size());
			}
			ImPlot::EndPlot();
		}

		if (ImPlot::BeginPlot("Projection Z"))
		{
			for (const ElectronBeam& eBeam : eBeamsToLookAt)
			{
				ImPlot::PlotLine(eBeam.label.c_str(), eBeam.zAxis.data(), eBeam.projectionValuesZ.data(), eBeam.zAxis.size());
			}
			ImPlot::EndPlot();
		}

		ImPlot::PushColormap(9);
		if (selectedIndex >= 0)
		{
			const ElectronBeam& sliceBeam = eBeamsToLookAt.at(selectedIndex);

			if (ImPlot::BeginPlot("XY Slice"))
			{
				ImPlot::SetupAxes("x", "y");
				ImPlot::PlotHeatmap(sliceBeam.label.c_str(), sliceBeam.slice.values.data(), sliceBeam.slice.nRows,
					sliceBeam.slice.nCols, sliceBeam.slice.minValue, sliceBeam.slice.maxValue, nullptr,
					sliceBeam.slice.bottomLeft, sliceBeam.slice.topRight);
				ImPlot::EndPlot();
			}
			ImGui::SameLine();
			ImPlot::ColormapScale("##HeatScale", sliceBeam.slice.minValue, sliceBeam.slice.maxValue, ImVec2(100, -1));

		}
		ImPlot::PopColormap();

		ImPlot::EndSubplots();
	}
}

TH3D* ElectronBeamWindow::CutZerosFromDistribution(TH3D* input)
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

	int newBinsX = maxX - minX + 1;
	int newBinsY = maxY - minY + 1;
	int newBinsZ = input->GetNbinsZ();

	double xMin = input->GetXaxis()->GetBinLowEdge(minX);
	double xMax = input->GetXaxis()->GetBinUpEdge(maxX);
	double yMin = input->GetYaxis()->GetBinLowEdge(minY);
	double yMax = input->GetYaxis()->GetBinUpEdge(maxY);
	double zMin = input->GetZaxis()->GetXmin();
	double zMax = input->GetZaxis()->GetXmax();

	TH3D* output = new TH3D(input->GetName(), input->GetTitle(),
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

TH3D* ElectronBeamWindow::MirrorDistributionAtZ(TH3D* input)
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
	TH3D* mirrored = new TH3D(input->GetName(), input->GetTitle(),
		nBinsX, xMin, xMax,
		nBinsY, yMin, yMax,
		2 * nBinsZ, zMin, zMax);

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

void ElectronBeamWindow::CreateLargeDistribution()
{
	int numberBinsX = m_distribution->GetNbinsX() * factor;
	int numberBinsY = m_distribution->GetNbinsY() * factor;
	int numberBinsZ = m_distribution->GetNbinsZ() * factor;

	double firstBinCenterX = m_distribution->GetXaxis()->GetBinCenter(1);
	double firstBinCenterY = m_distribution->GetYaxis()->GetBinCenter(1);
	double firstBinCenterZ = m_distribution->GetZaxis()->GetBinCenter(1);

	double lastBinCenterX = m_distribution->GetXaxis()->GetBinCenter(m_distribution->GetNbinsX());
	double lastBinCenterY = m_distribution->GetYaxis()->GetBinCenter(m_distribution->GetNbinsY());
	double lastBinCenterZ = m_distribution->GetZaxis()->GetBinCenter(m_distribution->GetNbinsZ());

	TH3D* temp = new TH3D("electron distribution large", "electron distribution large",
		numberBinsX, m_distribution->GetXaxis()->GetXmin(), m_distribution->GetXaxis()->GetXmax(),
		numberBinsY, m_distribution->GetYaxis()->GetXmin(), m_distribution->GetYaxis()->GetXmax(),
		numberBinsZ, m_distribution->GetZaxis()->GetXmin(), m_distribution->GetZaxis()->GetXmax());

	for (int i = 1; i <= numberBinsX; i++) 
	{
		for (int j = 1; j <= numberBinsY; j++) 
		{
			for (int k = 1; k <= numberBinsZ; k++)
			{
				double x = temp->GetXaxis()->GetBinCenter(i);
				double y = temp->GetYaxis()->GetBinCenter(j);
				double z = temp->GetZaxis()->GetBinCenter(k);

				double x_modified = std::min(std::max(x, firstBinCenterX), lastBinCenterX - 1e-5);
				double y_modified = std::min(std::max(y, firstBinCenterY), lastBinCenterY - 1e-5);
				double z_modified = std::min(std::max(z, firstBinCenterZ), lastBinCenterZ - 1e-5);

				//std::cout << x << ", " << y << ", " << z << "\n";
				double value = m_distribution->Interpolate(x_modified, y_modified, z_modified);
				temp->SetBinContent(i, j, k, value);
			}
		}
	}
	delete m_distribution;
	m_distribution = temp;
}

TH3D* ElectronBeamWindow::GenerateElectronBeamDensity()
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
				if (!activeDist.simplifyParams.noElectronBeamBend)
				{
					ymean = Trajectory(z);
				}
				if (activeDist.simplifyParams.gaussianElectronBeam)
				{
					// Calculate the value using the Gaussian distribution centered at z = 0
					value = exp(-(x * x + (y - ymean) * (y - ymean)) / (2.0 * pow(activeDist.simplifyParams.electronBeamRadius, 2)));
				}
				else if (activeDist.simplifyParams.cylindricalElectronBeam)
				{
					if (x * x + (y - ymean) * (y - ymean) <= pow(activeDist.simplifyParams.electronBeamRadius, 2))
					{
						// if inside the cylinder, set the value to an arbitrary constant value
						value = 1;
					}
				}
				eBeam->SetBinContent(i, j, k, value);
			}
		}
	}
	TH3D* newBeam = CutZerosFromDistribution(eBeam);
	return newBeam;
}

void ElectronBeamWindow::PlotTrajectory()
{
	m_mainCanvas->cd(1);
	
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

	TVector3 tangent = GetDirection(selectedPoint);
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

void ElectronBeamWindow::Plot3DBeam(TH3D* eBeam)
{
	if (!IsCanvasShown(m_mainCanvas) || !eBeam || !m_mainCanvas) return;

	m_mainCanvas->cd(2);

	delete m_distributionSmall;
	m_distributionSmall = (TH3D*)eBeam->Rebin3D(s_rebinningFactors[0],
		s_rebinningFactors[1],
		s_rebinningFactors[2], "dist small");

	m_distributionSmall->GetXaxis()->SetTitle("x-axis");
	m_distributionSmall->GetYaxis()->SetTitle("y-axis");
	m_distributionSmall->GetZaxis()->SetTitle("z-axis");
	m_distributionSmall->Draw("BOX2 COLZ");

	
}

// returns the y value as function of z
double ElectronBeamWindow::Trajectory(double z) const
{
	z = TMath::Abs(z);
    return pow(10, (-7.374 + 0.9 * z + 19.8 * z * z + 2.055 * z * z * z - 20.56 * z * z * z * z)) - pow(10, -7.374) + 1.0e-5;
}

double ElectronBeamWindow::Derivative(double z) const
{
	double sign = TMath::Sign(1, z);
	z = TMath::Abs(z);
	return sign * pow(10, (0.9 * z + 19.8 * z * z + 2.055 * z * z * z - 20.56 * z * z * z * z - 7.374)) * (39.6 * z + 6.165 * z * z - 82.24 * z * z * z + 0.9) * log(10);
}

double ElectronBeamWindow::DistancePointToTrajectoryOfZ(double z, Point3D point)
{
	return TMath::Power(Trajectory(z) - point.y, 2) + TMath::Power(z - point.z, 2);
}

void ElectronBeam::FillData()
{
	xAxis.reserve(fullHistogram->GetNbinsX());
	yAxis.reserve(fullHistogram->GetNbinsY());
	zAxis.reserve(fullHistogram->GetNbinsZ());

	projectionValuesX.reserve(fullHistogram->GetNbinsX());
	projectionValuesY.reserve(fullHistogram->GetNbinsY());
	projectionValuesZ.reserve(fullHistogram->GetNbinsZ());

	for (int i = 1; i <= fullHistogram->GetNbinsX(); i++)
	{
		xAxis.push_back(fullHistogram->GetXaxis()->GetBinCenter(i));
	}
	for (int i = 1; i <= fullHistogram->GetNbinsY(); i++)
	{
		yAxis.push_back(fullHistogram->GetYaxis()->GetBinCenter(i));
	}
	for (int i = 1; i <= fullHistogram->GetNbinsZ(); i++)
	{
		zAxis.push_back(fullHistogram->GetZaxis()->GetBinCenter(i));
	}

	TH1D* projectionX = fullHistogram->ProjectionX();
	TH1D* projectionY = fullHistogram->ProjectionY();
	TH1D* projectionZ = fullHistogram->ProjectionZ();

	for (int i = 1; i <= projectionX->GetNbinsX(); i++)
	{
		projectionValuesX.push_back(projectionX->GetBinContent(i));
	}
	for (int i = 1; i <= projectionY->GetNbinsX(); i++)
	{
		projectionValuesY.push_back(projectionY->GetBinContent(i));
	}
	for (int i = 1; i <= projectionZ->GetNbinsX(); i++)
	{
		projectionValuesZ.push_back(projectionZ->GetBinContent(i));
	}

	delete projectionX;
	delete projectionY;
	delete projectionZ;
}

ElectronBeam::ElectronBeam(ElectronBeam&& other) noexcept
{
	fullHistogram = other.fullHistogram;
	other.fullHistogram = nullptr;

	xAxis = std::move(other.xAxis);
	yAxis = std::move(other.yAxis);
	zAxis = std::move(other.zAxis);

	projectionValuesX = std::move(other.projectionValuesX);
	projectionValuesY = std::move(other.projectionValuesY);
	projectionValuesZ = std::move(other.projectionValuesZ);

	slice = std::move(other.slice);

	label = std::move(other.label);
}

ElectronBeam& ElectronBeam::operator=(ElectronBeam&& other) noexcept
{
	if (this == &other) return *this;

	delete fullHistogram;
	fullHistogram = other.fullHistogram;
	other.fullHistogram = nullptr;

	xAxis = std::move(other.xAxis);
	yAxis = std::move(other.yAxis);
	zAxis = std::move(other.zAxis);

	projectionValuesX = std::move(other.projectionValuesX);
	projectionValuesY = std::move(other.projectionValuesY);
	projectionValuesZ = std::move(other.projectionValuesZ);

	slice = std::move(other.slice);

	label = std::move(other.label);

	return *this;
}
