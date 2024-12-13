
void CrossSectionManager::FitWithSVD()
{
	// solve Ax = b with SVD

	EnergyDistributionManager* model = (EnergyDistributionManager*)Window::Get("Energy Distribution Manager");
	std::vector<EnergyDistribution>& energyDistributionList = model->GetEnergyDistributions();

	if (model->GetEnergyDistributions().empty())
	{
		std::cout << "no energy distributions\n";
		return;
	}

	SetupFitCrossSectionHist();
	CalculatePsis();

	binValuesFit.clear();
	binCentersFit.clear();
	binValuesFit.reserve(crossSectionFit->GetNbinsX());
	binCentersFit.reserve(crossSectionFit->GetNbinsX());

	std::vector<int> parameterIndeces;
	for (int i = 0; i < crossSectionFit->GetNbinsX(); i++)
	{
		parameterIndeces.push_back(i);
	}
	int nZeros = 10;

	TMatrixD PsiMatrix;
	TVectorD alphaVector;
	int p = parameterIndeces.size(); // crossSectionFit->GetNbinsX();
	std::vector<double> parameterResult;
	parameterResult.resize(p);

	while (nZeros > 0)
	{
		int n = energyDistributionList.size();
		int p2 = parameterIndeces.size();

		std::cout << "n: " << n << " p2: " << p2 << std::endl;
		PsiMatrix.Clear();
		alphaVector.Clear();

		// matrix A is all the p Psis for all the n distributions, n >= p is required so the rest is filled with 0
		PsiMatrix.ResizeTo(std::max(n, p2), p2);
		// vector b with all the rate coeffictions
		alphaVector.ResizeTo(std::max(n, p2));

		// fill matrix and vector
		for (int i = 0; i < std::max(n, p2); i++)
		{
			for (int j = 0; j < p2; j++)
			{
				// fill matrix and vector with 0 if p > n
				if (i >= n)
				{
					PsiMatrix[i][j] = 0;
				}
				else
				{
					PsiMatrix[i][j] = energyDistributionList[i].psi[parameterIndeces[j]];
				}
			}
			if (i >= n)
			{
				alphaVector[i] = 0;
			}
			else
			{
				alphaVector[i] = energyDistributionList[i].rateCoefficient;
			}
		}

		alphaVector.Print();

		TDecompSVD decomp(PsiMatrix);
		// solve for x, result is put into input vector
		decomp.Solve(alphaVector);

		alphaVector.Print();

		// remove all negative parameters
		nZeros = 0;
		for (int i = p2 - 1; i >= 0; i--)
		{
			if (alphaVector[parameterIndeces[i]] < 0)
			{
				parameterResult[parameterIndeces[i]] = 0;
				parameterIndeces.erase(parameterIndeces.begin() + i);
				nZeros++;
			}
			else
			{
				parameterResult[parameterIndeces[i]] = alphaVector[parameterIndeces[i]];
			}
		}
		std::cout << "zeros: " << nZeros << std::endl;
	}

	FillFitPlots(parameterResult.data());
}

void CrossSectionManager::FitWithEigenSVD()
{
	EnergyDistributionManager* model = (EnergyDistributionManager*)Window::Get("Energy Distribution Manager");
	std::vector<EnergyDistribution>& energyDistributionList = model->GetEnergyDistributions();

	if (model->GetEnergyDistributions().empty())
	{
		std::cout << "no energy distributions\n";
		return;
	}

	SetupFitCrossSectionHist();
	CalculatePsis();

	binValuesFit.clear();
	binCentersFit.clear();
	binValuesFit.reserve(crossSectionFit->GetNbinsX());
	binCentersFit.reserve(crossSectionFit->GetNbinsX());

	int n = energyDistributionList.size();
	int p = crossSectionFit->GetNbinsX();

	Eigen::MatrixXd PsiMatrix(n, p);
	Eigen::VectorXd alphaVector(n);

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < p; j++)
		{
			// fill matrix and vector with 0 if p > n
			PsiMatrix(i, j) = energyDistributionList[i].psi[j];
		}
		alphaVector[i] = energyDistributionList[i].rateCoefficient;
	}

	Eigen::JacobiSVD<Eigen::MatrixXd> svd(PsiMatrix, Eigen::ComputeThinU | Eigen::ComputeThinV);

	// Solve using the pseudoinverse (x = A^+ * b), where A^+ is the Moore-Penrose pseudoinverse
	Eigen::VectorXd result = svd.solve(alphaVector);

	std::cout << result << std::endl;

	FillFitPlots(result.data());
}

// Project onto the non-negative orthant (set negative values to 0)
Eigen::VectorXd projectOntoNonNegative(const Eigen::VectorXd& x) {
	Eigen::VectorXd proj = x;
	proj = proj.cwiseMax(0);  // Set all negative values to 0
	return proj;
}

void CrossSectionManager::FitWithEigenGD()
{
	EnergyDistributionManager* model = (EnergyDistributionManager*)Window::Get("Energy Distribution Manager");
	std::vector<EnergyDistribution>& energyDistributionList = model->GetEnergyDistributions();

	if (model->GetEnergyDistributions().empty())
	{
		std::cout << "no energy distributions\n";
		return;
	}
	if (initialGuess.empty())
	{
		// create Fit cross section
		SetupFitCrossSectionHist();
		CalculatePsis();
		SetupInitialGuess();
	}
	else
	{
		initialGuess.clear();
		initialGuess = binValuesFit;
	}

	binValuesFit.clear();
	binCentersFit.clear();
	binValuesFit.reserve(crossSectionFit->GetNbinsX());
	binCentersFit.reserve(crossSectionFit->GetNbinsX());

	int n = energyDistributionList.size();
	int p = crossSectionFit->GetNbinsX();

	Eigen::MatrixXd PsiMatrix(n, p);
	Eigen::VectorXd alphaVector(n);

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < p; j++)
		{
			// fill matrix and vector with 0 if p > n
			PsiMatrix(i, j) = energyDistributionList[i].psi[j];
		}
		alphaVector[i] = energyDistributionList[i].rateCoefficient;
	}

	// Initial guess for x
	Eigen::Map<Eigen::VectorXd> x(initialGuess.data(), initialGuess.size());
	//Eigen::VectorXd x = Eigen::VectorXd::Zero(p);
	std::cout << x << std::endl;
	std::cout << alphaVector << std::endl;

	// Gradient descent parameters
	double factor = 1;
	Eigen::MatrixXd D = Eigen::MatrixXd::Zero(p - 2, p);
	for (int i = 0; i < p - 2; i++)
	{
		//factor = 10 * pow(i, 3);
		factor = 10000 * pow(crossSectionFit->GetBinCenter(i + 1), 1.5);
		D(i, i) = 1 * factor;       // x_i
		D(i, i + 1) = -2 * factor;  // -2 * x_{i+1}
		D(i, i + 2) = 1 * factor;   // x_{i+2}
	}

	// Precompute D^T D
	Eigen::MatrixXd DtD = D.transpose() * D;
	//Eigen::VectorXd lambdaVector(p);
	//for (int i = 0; i < p; i++)
	//{
	//	lambdaVector[i] = lambda * crossSectionFit->GetBinCenter(i);
	//}

	double tolerance = 1e-6;


	// Iterate using gradient descent
	for (int iter = 0; iter < iterations; ++iter) {
		// Compute the residual: r = A*x - b
		Eigen::VectorXd r = PsiMatrix * x - alphaVector;
		//std::cout << "residual: " << r << std::endl;
		// Compute the gradient: grad = A^T * r
		//std::cout << "test: " <<  lambdaVector * (DtD * x) << std::endl;
		//std::cout << "end: " << std::endl;
		//Eigen::VectorXd bla = lambdaVector.array() * (DtD * x).array();
		Eigen::VectorXd grad = PsiMatrix.transpose() * r + lambda * (DtD * x);
		//grad += 2 * lambda * x;
		//std::cout << "gradient: " << grad << std::endl;
		//std::cout << "gradient modified: " << grad.cwiseMin(0.001 / learningRate).cwiseMax(-0.001 / learningRate) << std::endl;
		// Update x using gradient descent step
		//double noise_scale = 0.00001 / (iter + 1);
		//Eigen::VectorXd noise = Eigen::VectorXd::Random(n);
		//std::cout << "random noise: " << noise << std::endl;
		x = x - learningRate * grad.cwiseMin(0.001 / learningRate).cwiseMax(-0.001 / learningRate);
		//+ noise * noise_scale);
	//std::cout << "new x: " << x << std::endl;
	// Project onto the non-negative orthant (i.e., enforce x >= 0)

		x = x.cwiseAbs();//projectOntoNonNegative(x);

		// Check for convergence (if the gradient is small enough, stop)
		if (grad.norm() < tolerance)
		{
			std::cout << "Converged in " << iter << " iterations." << std::endl;
			break;
		}
	}

	// Output the solution
	std::cout << "Solution x (with non-negativity constraints): " << std::endl << x << std::endl;
	FillFitPlots(x.data());
}

torch::Tensor CrossSectionManager::custom_loss(const torch::Tensor& x, const torch::Tensor& A, const torch::Tensor& b)
{
	// Compute the residual
	torch::Tensor residual = torch::matmul(A, x) - b;

	// Compute the residual loss (MSE)
	torch::Tensor loss_residual = torch::mean(torch::pow(residual, 2));

	// Add constraints
	// L2 regularization
	torch::Tensor l2_reg = torch::mean(torch::pow(x, 2));

	torch::Tensor negative_penalty = torch::relu(-x).sum();

	// Second derivative (smoothness constraint)
	if (x.size(0) >= 3) { // Ensure there are at least 3 elements to compute second differences
		torch::Tensor second_diff = x.slice(0, 2) - 2 * x.slice(0, 1, -1) + x.slice(0, 0, -2);
		torch::Tensor smoothness = torch::mean(torch::pow(second_diff, 2)); // Penalize curvature
		return loss_residual + 1000 * negative_penalty + l2regularisation * l2_reg + smoothRegularisation * smoothness;
	}
	else {
		// If x has fewer than 3 elements, return only the residual and L2 regularization
		return loss_residual + 1000 * negative_penalty + l2regularisation * l2_reg;
	}
}

void CrossSectionManager::FitWithTorch()
{
	EnergyDistributionManager* model = (EnergyDistributionManager*)Window::Get("Energy Distribution Manager");
	std::vector<EnergyDistribution>& energyDistributionList = model->GetEnergyDistributions();

	if (model->GetEnergyDistributions().empty())
	{
		std::cout << "no energy distributions\n";
		return;
	}
	if (initialGuess.empty())
	{
		// create Fit cross section
		SetupFitCrossSectionHist();
		CalculatePsis();
		SetupInitialGuess();
	}
	else
	{
		initialGuess.clear();
		initialGuess = binValuesFit;
	}

	binValuesFit.clear();
	binCentersFit.clear();
	binValuesFit.reserve(crossSectionFit->GetNbinsX());
	binCentersFit.reserve(crossSectionFit->GetNbinsX());

	int n = energyDistributionList.size();
	int p = crossSectionFit->GetNbinsX();

	torch::Tensor PsiMatrix = torch::zeros({ n, p }, torch::kDouble);
	torch::Tensor alphaVector = torch::zeros({ n }, torch::kDouble);

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < p; j++)
		{
			// fill matrix and vector with 0 if p > n
			PsiMatrix[i][j] = energyDistributionList[i].psi[j];
		}
		alphaVector[i] = energyDistributionList[i].rateCoefficient;
	}

	// Initial guess for x
	torch::Tensor x = torch::from_blob(initialGuess.data(), initialGuess.size(), torch::TensorOptions()
		.dtype(torch::kDouble)
		.requires_grad(true));

	//torch::Tensor x = torch::randn({ p }, torch::requires_grad(true));
	//std::cout << "matrix size: " << PsiMatrix.sizes() << std::endl;
	//std::cout << "alpha size: " << alphaVector.sizes() << std::endl;
	//std::cout << "parameter size: " << x.sizes() << std::endl;
	//std::cout << "psi matrix: \n" << PsiMatrix << std::endl;
	//std::cout << "alpha vector: \n" << alphaVector << std::endl;
	//std::cout << "parameters: \n" << x << std::endl;
	// 
	// Define an optimizer (e.g., SGD)
	//torch::optim::Adam optimizer({ x });
	//torch::optim::LBFGS optimizer({ x });
	torch::optim::SGD optimizer({ x }, torch::optim::SGDOptions(torchLearningRate).momentum(0.9));

	// Training loop to minimize the loss
	for (size_t epoch = 0; epoch < nEpochs; epoch++)
	{
		// Compute the loss
		torch::Tensor loss = custom_loss(x, PsiMatrix, alphaVector);

		// Backward pass and optimization step
		optimizer.zero_grad();
		loss.backward();
		optimizer.step();

		// Print progress
		if (epoch % (int)(nEpochs / 10) == 0 || epoch == nEpochs - 1)
		{
			std::cout << "Epoch [" << epoch + 1 << "], Loss: " << loss.item<float>() << std::endl;
		}
	}

	std::cout << "Optimized x:\n" << x << std::endl;
	FillFitPlots(x.data_ptr<double>());
}


void CrossSectionManager::FitCrossSectionHistogram()
{
	EnergyDistributionManager* model = (EnergyDistributionManager*)Window::Get("Energy Distribution Manager");
	if (model->GetEnergyDistributions().empty())
	{
		std::cout << "no energy distributions\n";
		return;
	}

	if (initialGuess.empty())
	{
		// create Fit cross section
		SetupFitCrossSectionHist();
		CalculatePsis();
		SetupInitialGuess();
	}
	else
	{
		initialGuess.clear();
		initialGuess = binValuesFit;
		for (double& value : initialGuess)
		{
			value = std::abs(value);
		}
	}
	TF1* fitFunction = new TF1("fit function", this, &CrossSection::FitFunction, 0, 99, crossSectionFit->GetNbinsX());

	if (limitParamRange)
	{
		for (int i = 0; i < crossSectionFit->GetNbinsX(); i++)
		{
			fitFunction->SetParLimits(i, 0, 1e30);
		}
	}

	fitFunction->SetParameters(initialGuess.data());


	binValuesFit.clear();
	binCentersFit.clear();
	binValuesFit.reserve(crossSectionFit->GetNbinsX());
	binCentersFit.reserve(crossSectionFit->GetNbinsX());

	rateCoefficients->Fit(fitFunction, "RN");

	double* parameter = fitFunction->GetParameters();
	FillFitPlots(parameter);

	fitFunction->Delete();
}

void CrossSectionManager::SetupInitialGuess()
{
	initialGuess.clear();

	for (int i = 1; i <= crossSectionFit->GetNbinsX(); i++)
	{
		double energy = crossSectionFit->GetBinCenter(i);
		double velocity = TMath::Sqrt(2 * energy * TMath::Qe() / PhysicalConstants::electronMass);
		double alpha = rateCoefficients->Eval(energy);
		initialGuess.push_back(alpha / velocity);
	}
}

void CrossSectionManager::SetupFitCrossSectionHist()
{
	EnergyDistributionManager* model = (EnergyDistributionManager*)Window::Get("Energy Distribution Manager");
	std::vector<EnergyDistribution>& energyDistributionList = model->GetEnergyDistributions();

	std::vector<double> binEdges;
	double maxEnergy = 93; //just a gues, not fixed
	double minEnergy = energyDistributionList.back().eBeamParameter.detuningEnergy / 10;

	// binning like in the paper 
	if (currentOption == PaperBinning)
	{
		EnergyDistribution& representativeEnergyDist = energyDistributionList.back();
		double kT_trans = representativeEnergyDist.eBeamParameter.transverse_kT;
		double kT_long = 0.00047; // representativeEnergyDist->eBeamParameter.longitudinal_kT; //eBeam->GetLongitudinal_kT(labEnergiesParameter.centerLabEnergy);

		binEdges.push_back(0);
		binEdges.push_back(kT_trans / 20);

		for (int i = 1; binEdges[i] < kT_trans; i++)
		{
			binEdges.push_back(binFactor * 2 * binEdges.back());
			//std::cout << binEdges.back() << "\n";
		}
		while (binEdges[binEdges.size() - 1] < maxEnergy)
		{
			double previousEdge = binEdges[binEdges.size() - 1];
			double delta_E = sqrt(pow((kT_trans * log(2)), 2) + 16 * log(2) * kT_long * previousEdge);
			binEdges.push_back(previousEdge + binFactor * delta_E);

			//std::cout << "delta E " << delta_E << "\n";
			//std::cout << binEdges[binEdges.size()] << "\n";
			//std::cout << (binEdges[binEdges.size()] < maxEnergy) << "\n";
		}
	}
	// constant binning
	if (currentOption == ConstantBinning)
	{
		binEdges.reserve(numberBins + 1);
		binEdges.push_back(0);
		double binWidth = maxEnergy / numberBins;
		for (int i = 0; i < numberBins; i++)
		{
			binEdges.push_back(binWidth * (i + 1));
		}
	}
	// bin width increses by a constant factor
	if (currentOption == FactorBinning)
	{
		if (limitBinSize)
		{
			double min = minEnergy; // std::max(energyRange[0], 1e-9f);
			double factor = TMath::Power((maxEnergy / min), (1.0 / numberBins));

			binEdges.reserve(numberBins + 1);
			binEdges.push_back(0);
			for (int i = 0; binEdges[binEdges.size() - 1] < maxEnergy; i++)
			{
				double nextBin = binEdges[i] * factor;
				double difference = std::max(nextBin - binEdges[i], minBinSize);
				binEdges.push_back(binEdges[i] + difference);
			}
		}
		else
		{
			float min = minEnergy;
			double factor = TMath::Power((maxEnergy / min), (1.0 / numberBins));

			binEdges.reserve(numberBins + 1);
			binEdges.push_back(min);
			for (int i = 0; i < numberBins; i++)
			{
				binEdges.push_back(binEdges[i] * factor);
			}
		}
	}
	if (currentOption == PaperFactorMix)
	{
		EnergyDistribution& representativeEnergyDist = energyDistributionList.back();
		double kT_trans = representativeEnergyDist.eBeamParameter.transverse_kT;
		//double kT_long = 0.00047; // representativeEnergyDist->eBeamParameter.longitudinal_kT; //eBeam->GetLongitudinal_kT(labEnergiesParameter.centerLabEnergy);

		binEdges.push_back(0);
		binEdges.push_back(kT_trans / 20);

		for (int i = 1; binEdges[i] < kT_trans; i++)
		{
			binEdges.push_back(2 * binEdges[i]);
			//std::cout << binEdges[i + 1] << "\n";
		}

		double factor = TMath::Power((maxEnergy / binEdges.back()), (1.0 / numberBins));

		for (int i = 0; i < numberBins; i++)
		{
			binEdges.push_back(binEdges.back() * factor);
		}
	}

	if (currentOption == Paper_FWHM)
	{

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

	crossSectionFit = new TH1D("cross section fit", "cross section fit", binEdges.size() - 1, binEdges.data());
}
