#include "ElectronBeam.h"
#include "PhysicalConstants.h"
#include "FileLoader.h"

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
	: Module("electron beam")
{
	PlotTrajectory();
}

ElectronBeamParameters& ElectronBeam::GetParamters()
{
	return parameters;
}

void ElectronBeam::LoadDensityFile(std::filesystem::path file)
{
	if (!file.empty())
	{
		delete m_distribution;
		m_distribution = FileLoader::getInstance().LoadMatrixFile(file);
		m_distribution->SetTitle("electron m_distribution");
		m_distribution->SetName("electron m_distribution");
	}
}

bool ElectronBeam::HasDensityChanged()
{
	return densityChanged;
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
	densityChanged = false;
	if (ImGui::Button("Load e-density file"))
	{
		std::filesystem::path file = FileLoader::getInstance().openFileExplorer("data\\e-densities\\");
		LoadDensityFile(file);
		PlotDistribution();
		densityChanged = true;
	}

	ImGui::Text("electron current: %e A", parameters.electronCurrent);
	ImGui::InputDouble("cooling energy [eV]", &parameters.coolingEnergy);
	ImGui::InputDouble("transverse kT [eV]", &parameters.transverse_kT);
	ImGui::InputDouble("cathode radius [m]", &parameters.cathodeRadius);
	ImGui::InputDouble("cathode Temperature [K]", &parameters.cathodeTemperature);
	ImGui::InputDouble("extraction energy [eV]", &parameters.extractionEnergy);
	ImGui::InputDouble("expansion factor", &parameters.expansionFactor);
	ImGui::InputDouble("LLR", &parameters.LLR);
	ImGui::InputDouble("sigma lab energy [eV]", &parameters.sigmaLabEnergy);
	ImGui::Separator();

	if (ImGui::SliderFloat("z", &sliderZ, -0.7, 0.7))
	{
		PlotTrajectory();
	}
	if (ImGui::SliderFloat("y", &sliderY, 0, 0.05, "%.4f"))
	{
		PlotTrajectory();
	}
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

