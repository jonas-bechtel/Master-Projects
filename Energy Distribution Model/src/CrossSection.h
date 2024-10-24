#pragma once
#include "Module.h"

class CrossSection : public Module
{
public:
	CrossSection();

private:
	void ShowUI() override;

private:
	TH1D* crossSection = nullptr;
};

