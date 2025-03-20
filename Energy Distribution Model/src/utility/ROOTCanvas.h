#pragma once

class ROOTCanvas : public TCanvas
{
public:
	ROOTCanvas(const char* name, const char* title, int width, int height);

	void MakeShowHideButton();
	bool IsShown();
	void Render();

private:
	void Show();
	void Hide();
	
private:
	bool m_shown = false;

	static TApplication app;
};

