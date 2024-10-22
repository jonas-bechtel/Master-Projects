#include "ElectronBeam.h"
#include "PhysicalConstants.h"
#include "FileHandler.h"

#include "imgui.h"

#include <math.h>

#include <TView.h>
#include <TAxis3D.h>
#include <TArrow.h>
#include <TMath.h>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include <TGraph.h>
#include <TVector3.h>
#include <TRootCanvas.h>

ElectronBeam::ElectronBeam()
	: Module("Electron Beam")
{
	PlotTrajectory();
}

ElectronBeamParameters ElectronBeam::GetParameter()
{
	return parameters;
}

TH3D* ElectronBeam::GetDistribution()
{
	if (parameters.hasGaussianShape)
	{
		return generatedBeamDensity;
	}
	return m_distribution;
}

void ElectronBeam::SetParameter(ElectronBeamParameters params)
{
	parameters = params;
}

void ElectronBeam::SetCurrent(double current)
{
	parameters.electronCurrent = current;
}

void ElectronBeam::LoadDensityFile(std::filesystem::path file)
{
	if (!file.empty())
	{
		delete m_distribution;
		m_distribution = FileHandler::GetInstance().LoadMatrixFile(file);
		CutZerosFromDistribution();
		m_distribution->SetTitle("electron distribution");
		m_distribution->SetName("electron distribution");

		if(increaseHist)
			CreateLargeDistribution();
		
		m_distribution->GetXaxis()->SetTitle("x-axis");
		m_distribution->GetYaxis()->SetTitle("y-axis");
		m_distribution->GetZaxis()->SetTitle("z-axis");

		parameters.densityFile = file;
	}
}

TVector3 ElectronBeam::GetDirection(double z)
{
	double derivative = Derivative(z);
	TVector3 direction(0, derivative, 1);
	return direction.Unit();
}

TVector3 ElectronBeam::GetDirection(Point3D point)
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

TVector3 ElectronBeam::GetNormal(double z)
{
	return GetDirection(z).Orthogonal();
}

double ElectronBeam::GetLongitudinal_kT(double labEnergy)
{
	if (parameters.hasFixedLongitudinalTemperature) 
	{
		return parameters.longitudinal_kT;
	}

	// intermediate unit is [J] final unit is [eV]
	double A = (1. + TMath::Power((parameters.expansionFactor - 1.) / parameters.expansionFactor, 2.)) 
		* (TMath::Power(TMath::K() * parameters.cathodeTemperature, 2.)) / (2. * TMath::Qe() * labEnergy);
	double B = 2.544008e-27;
	double C = parameters.LLR * TMath::Power(parameters.cathodeRadius, -2 / 3.) * TMath::Power(parameters.electronCurrent, 1 / 3.)
		* TMath::Power(parameters.extractionEnergy * TMath::Qe(), -1 / 6.) * parameters.extractionEnergy / labEnergy;
	double D = TMath::Power(parameters.sigmaLabEnergy * TMath::Qe(), 2.) / (2. * TMath::Qe() * labEnergy);
	return (A + B * C + D) / TMath::Qe();
}

double ElectronBeam::GetTransverse_kT()
{
	return parameters.transverse_kT;
}

void ElectronBeam::ShowUI()
{
	if (ImGui::Button("Load e-density file"))
	{
		std::filesystem::path file = FileHandler::GetInstance().OpenFileExplorer();
		LoadDensityFile(file);
		PlotDistribution();
		PlotProjections();

		IonBeam* ionBeam = (IonBeam*)Module::Get("Ion Beam");
		TH3D* result = ionBeam->MultiplyWithElectronDensities(m_distribution);
		ionBeam->PlotDistribution();

		MCMC* mcmc = (MCMC*)Module::Get("MCMC");
		mcmc->SetTargetDistribution(result);
		mcmc->PlotTargetDistribution();
	}
	ImGui::SameLine();
	ImGui::Checkbox("multiply bins", &increaseHist);
	ImGui::SameLine();
	ImGui::BeginDisabled(!increaseHist);
	ImGui::SetNextItemWidth(100.0f);
	ImGui::InputInt(" factor", &factor, 2);
	ImGui::EndDisabled();

	ImGui::Separator();
	if (ImGui::Checkbox("use gaussian beam shape", &parameters.hasGaussianShape))
	{
		parameters.densityFile.clear();
	}
	ImGui::BeginDisabled(!parameters.hasGaussianShape);
	ImGui::SetNextItemWidth(100.0f);
	ImGui::InputDouble("radius [m]", &parameters.radius);
	ImGui::SameLine();
	ImGui::Checkbox("exclude bend", &parameters.hasNoBending);
	if (ImGui::Button("generate density"))
	{
		GenerateElectronBeamDensity();
		PlotGeneratedDensities();
		PlotProjections();
		
		IonBeam* ionBeam = (IonBeam*)Module::Get("Ion Beam");
		TH3D* result = ionBeam->MultiplyWithElectronDensities(generatedBeamDensity);
		ionBeam->PlotDistribution();

		MCMC* mcmc = (MCMC*)Module::Get("MCMC");
		mcmc->SetTargetDistribution(result);
		mcmc->PlotTargetDistribution();
	}
	ImGui::EndDisabled();

	ImGui::Separator();
	ImGui::Checkbox("use fixed longitudinal kT", &parameters.hasFixedLongitudinalTemperature); ImGui::SetNextItemWidth(100.0f);
	ImGui::BeginDisabled(!parameters.hasFixedLongitudinalTemperature);
	ImGui::InputDouble("longitudinal kT [eV]", &parameters.longitudinal_kT);
	ImGui::EndDisabled();
	ImGui::Text("electron current: %e A", parameters.electronCurrent);			ImGui::SetNextItemWidth(100.0f);
	ImGui::InputDouble("cooling energy [eV]", &parameters.coolingEnergy);		ImGui::SetNextItemWidth(100.0f);
	ImGui::InputDouble("transverse kT [eV]", &parameters.transverse_kT);		ImGui::SetNextItemWidth(100.0f);
	ImGui::BeginDisabled(parameters.hasFixedLongitudinalTemperature);
	ImGui::InputDouble("cathode radius [m]", &parameters.cathodeRadius);		ImGui::SetNextItemWidth(100.0f);
	ImGui::InputDouble("cathode Temperature [K]", &parameters.cathodeTemperature); ImGui::SetNextItemWidth(100.0f);
	ImGui::InputDouble("extraction energy [eV]", &parameters.extractionEnergy); ImGui::SetNextItemWidth(100.0f);
	ImGui::InputDouble("expansion factor", &parameters.expansionFactor);		ImGui::SetNextItemWidth(100.0f);
	ImGui::InputDouble("LLR", &parameters.LLR);									ImGui::SetNextItemWidth(100.0f);
	ImGui::InputDouble("sigma lab energy [eV]", &parameters.sigmaLabEnergy);
	ImGui::EndDisabled();
	ImGui::Separator();

	if (ImGui::SliderFloat("z", &sliderZ, -0.7f, 0.7f))
	{
		PlotTrajectory();
	}
	if (ImGui::SliderFloat("y", &sliderY, 0, 0.05f, "%.4f"))
	{
		PlotTrajectory();
	}
}

void ElectronBeam::CutZerosFromDistribution()
{
	// Identify the non-zero bin ranges
	int minX = m_distribution->GetNbinsX(), maxX = 0;
	int minY = m_distribution->GetNbinsY(), maxY = 0;
	int minZ = 0;
	int maxZ = m_distribution->GetNbinsY();

	// Loop over all bins to find the first and last non-zero bins in both x and y
	for (int i = 1; i <= m_distribution->GetNbinsX(); i++)
	{
		for (int j = 1; j <= m_distribution->GetNbinsY(); j++)
		{
			for (int k = 1; k <= m_distribution->GetNbinsZ(); k++)
			{
				if (m_distribution->GetBinContent(i, j, k) != 0)
				{
					if (i < minX) minX = i;
					if (i > maxX) maxX = i;
					if (j < minY) minY = j;
					if (j > maxY) maxY = j;
				}
			}
			
		}
	}
	//std::cout << "minX: " << minX << "\n";
	//std::cout << "maxX: " << maxX << "\n";
	//std::cout << "minY: " << minY << "\n";
	//std::cout << "maxY: " << maxY << "\n";

	int newBinsX = maxX - minX + 1;
	int newBinsY = maxY - minY + 1;
	int newBinsZ = m_distribution->GetNbinsZ();

	//std::cout << "newBinsX: " << newBinsX << "\n";
	//std::cout << "newBinsY: " << newBinsY << "\n";

	double xMin = m_distribution->GetXaxis()->GetBinLowEdge(minX);
	double xMax = m_distribution->GetXaxis()->GetBinUpEdge(maxX);
	double yMin = m_distribution->GetYaxis()->GetBinLowEdge(minY);
	double yMax = m_distribution->GetYaxis()->GetBinUpEdge(maxY);
	double zMin = m_distribution->GetZaxis()->GetXmin();
	double zMax = m_distribution->GetZaxis()->GetXmax();

	//std::cout << "xMin: " << xMin << "\n";
	//std::cout << "xMax: " << xMax << "\n";

	TH3D* temp = new TH3D("temp", "temp",
		newBinsX, xMin, xMax,
		newBinsY, yMin, yMax,
		newBinsZ, zMin, zMax);

	// Step 4: Copy the relevant bin contents to the new histogram
	for (int i = minX; i <= maxX; i++) 
	{
		for (int j = minY; j <= maxY; j++) 
		{
			for (int k = minZ; k <= maxZ; k++)
			{
				double content = m_distribution->GetBinContent(i, j, k);
				temp->SetBinContent(i - minX + 1, j - minY + 1, k, content); // Copy with adjusted indices
			}
		}
	}
	//std::cout << "old bins\n";
	//for (int i = 1; i <= m_distribution->GetNbinsX(); i++)
	//{
	//	std::cout << m_distribution->GetXaxis()->GetBinCenter(i) << "\n";
	//}
	//std::cout << "new bins\n";
	//for (int i = 1; i <= temp->GetNbinsX(); i++)
	//{
	//	std::cout << temp->GetXaxis()->GetBinCenter(i) << "\n";
	//}

	delete m_distribution;
	m_distribution = temp;
}

void ElectronBeam::CreateLargeDistribution()
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

void ElectronBeam::GenerateElectronBeamDensity()
{
	int nXBins = 50;
	int nYBins = 50;
	int nZBins = 50;

	double xmin = -0.04;
	double xmax = 0.04;
	double ymin = -0.04;
	double ymax = 0.04;
	double zmin = 0.0;
	double zmax = 0.65;

	delete generatedBeamDensity;
	generatedBeamDensity = new TH3D("generated densites", "generated densites", nXBins, xmin, xmax,
																				nYBins, ymin, ymax,
																				nZBins, zmin, zmax);

	for (int i = 1; i <= nXBins; i++) {
		for (int j = 1; j <= nYBins; j++) {
			for (int k = 1; k <= nZBins; k++) {
				// Calculate the coordinates for this bin
				double x = generatedBeamDensity->GetXaxis()->GetBinCenter(i);
				double y = generatedBeamDensity->GetYaxis()->GetBinCenter(j);
				double z = generatedBeamDensity->GetZaxis()->GetBinCenter(k);

				double ymean = 0;
				if (!parameters.hasNoBending)
				{
					ymean = Trajectory(z);
				}
				// Calculate the value using the Gaussian distribution centered at z = 0
				double value = exp(-(x * x + (y - ymean) * (y - ymean)) / (2.0 * parameters.radius * parameters.radius));
				generatedBeamDensity->SetBinContent(i, j, k, value);
			}
		}
	}
}

void ElectronBeam::PlotGeneratedDensities()
{
	if (!generatedBeamDensity) return;
	if (generatedBeamDensitySmall) delete generatedBeamDensitySmall;

	m_mainCanvas->cd(2);
	generatedBeamDensitySmall = (TH3D*)generatedBeamDensity->Rebin3D(s_rebinningFactors[0],
																	 s_rebinningFactors[1],
																	 s_rebinningFactors[2], "generated densities small");

	generatedBeamDensitySmall->GetXaxis()->SetTitle("x-axis");
	generatedBeamDensitySmall->GetYaxis()->SetTitle("y-axis");
	generatedBeamDensitySmall->GetZaxis()->SetTitle("z-axis");
	generatedBeamDensitySmall->Draw("BOX2");
}

void ElectronBeam::PlotTrajectory()
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
	TArrow* tangentArrow = new TArrow(z0, y0, z0 + tangent.z(), y0 + tangent.y(), 0.02, "|>");
	tangentArrow->SetLineColor(kRed);
	tangentArrow->SetLineWidth(2);
	tangentArrow->SetAngle(30); // Arrowhead angle
	tangentArrow->Draw();

	TVector3 normal = GetNormal(z0);
	normal *= arrowLength / 20;

	// Draw the normal vector arrow
	TArrow* normalArrow = new TArrow(z0, y0, z0 + normal.z(), y0 + normal.y(), 0.02, "|>");
	normalArrow->SetLineColor(kBlue);
	normalArrow->SetLineWidth(2);
	normalArrow->SetAngle(30); // Arrowhead angle
	normalArrow->Draw();
}

void ElectronBeam::PlotProjections()
{
	if (m_distribution)
	{
		delete electronBeamProjectionX;
		delete electronBeamProjectionY;
		delete electronBeamProjectionZ;

		electronBeamProjectionX = m_distribution->ProjectionX();
		electronBeamProjectionY = m_distribution->ProjectionY();
		electronBeamProjectionZ = m_distribution->ProjectionZ();

		electronBeamProjectionX->Scale(1 / electronBeamProjectionX->GetMaximum());
		electronBeamProjectionY->Scale(1 / electronBeamProjectionY->GetMaximum());
		electronBeamProjectionZ->Scale(1 / electronBeamProjectionZ->GetMaximum());

		m_secondCanvas->cd(1);
		electronBeamProjectionX->Draw("hist");

		m_secondCanvas->cd(2);
		electronBeamProjectionY->Draw("hist");

		m_secondCanvas->cd(3);
		electronBeamProjectionZ->Draw("hist");
	}
	if (generatedBeamDensity)
	{
		std::cout << "hi\n";
		delete generatedBeamProjectionX;
		delete generatedBeamProjectionY;
		delete generatedBeamProjectionZ;

		generatedBeamProjectionX = generatedBeamDensity->ProjectionX();
		generatedBeamProjectionY = generatedBeamDensity->ProjectionY();
		generatedBeamProjectionZ = generatedBeamDensity->ProjectionZ();

		generatedBeamProjectionX->Scale(1 / generatedBeamProjectionX->GetMaximum());
		generatedBeamProjectionY->Scale(1 / generatedBeamProjectionY->GetMaximum());
		generatedBeamProjectionZ->Scale(1 / generatedBeamProjectionZ->GetMaximum());

		generatedBeamProjectionX->SetLineColor(kRed);
		generatedBeamProjectionY->SetLineColor(kRed);
		generatedBeamProjectionZ->SetLineColor(kRed);

		m_secondCanvas->cd(1);
		generatedBeamProjectionX->Draw("Same hist");

		m_secondCanvas->cd(2);
		generatedBeamProjectionY->Draw("Same hist");

		m_secondCanvas->cd(3);
		generatedBeamProjectionZ->Draw("Same hist");
	}

	

	//delete electronBeamProfileXY;
	//delete electronBeamProfileXZ;
	//delete electronBeamProfileYZ;

	

	//m_secondCanvas->cd(4);
	//electronBeamProfileXY = m_distribution->Project3DProfile("xy");
	//electronBeamProfileXY->Draw("COLZ");
	//
	//m_secondCanvas->cd(5);
	//electronBeamProfileXZ = m_distribution->Project3DProfile("xz");
	//electronBeamProfileXZ->Draw("COLZ");
	//
	//m_secondCanvas->cd(6);
	//electronBeamProfileYZ = m_distribution->Project3DProfile("yz");
	//electronBeamProfileYZ->Draw("COLZ");
}

// returns the y value as function of z
double ElectronBeam::Trajectory(double z)
{
	z = TMath::Abs(z);
    return pow(10, (-7.374 + 0.9 * z + 19.8 * z * z + 2.055 * z * z * z - 20.56 * z * z * z * z)) - pow(10, -7.374) + 1.0e-5;
}

double ElectronBeam::Derivative(double z)
{
	double sign = TMath::Sign(1, z);
	z = TMath::Abs(z);
	return sign * pow(10, (0.9 * z + 19.8 * z * z + 2.055 * z * z * z - 20.56 * z * z * z * z - 7.374)) * (39.6 * z + 6.165 * z * z - 82.24 * z * z * z + 0.9) * log(10);
}

double ElectronBeam::DistancePointToTrajectoryOfZ(double z, Point3D point)
{
	return TMath::Power(Trajectory(z) - point.y, 2) + TMath::Power(z - point.z, 2);
}

std::string ElectronBeamParameters::String()
{
	std::string string = std::string(Form("# electron beam parameter:\n"));
	if (hasGaussianShape)
	{
		string += std::string(Form("# using gaussian beam shape %s\n", (hasNoBending ? "excluding the bend" : "including the bend")));
		string += std::string(Form("# radius: %.4f m\n", radius));
	}
	else
	{
		string += std::string(Form("# density file: %s\n", densityFile.filename().string().c_str()));
	}
	if (hasFixedLongitudinalTemperature)
	{
		string += std::string(Form("# using fixed longitudial Temperature\n"));
		string += std::string(Form("# longitudinal kT: %.2e eV\n", longitudinal_kT));
	}
	else
	{
		string += std::string(Form("# electron current: %.2e A\n", electronCurrent));
		string += std::string(Form("# expansion factor: %.1f\n", expansionFactor));
		string += std::string(Form("# cathode radius: %.3e m\n", cathodeRadius));
		string += std::string(Form("# cathode tempperature: %.1f K\n", cathodeTemperature));
		string += std::string(Form("# LLR: %.1f\n", LLR));
		string += std::string(Form("# sigma lab energy: %.3f eV\n", sigmaLabEnergy));
		string += std::string(Form("# extraction energy: %.3f eV\n", extractionEnergy));
	}
						  
	string += std::string(Form("# cooling energy: %.3f eV\n", coolingEnergy));
	string += std::string(Form("# transverse kT: %.2e eV\n", transverse_kT));
	
	return string;
}
