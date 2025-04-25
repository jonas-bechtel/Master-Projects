#pragma once

namespace Application
{
	struct Settings
	{
		bool tooltipsDisabled = false;
		bool showImGuiDemoWindow = false;
		bool showImPlotDemoWindow = false;
	};

	void InitImGui();
	void Init();
	void Run();
	void ShowWindows();
	void ShutdownImGui();

	Settings& GetSettings();
}









//class Application
//{
//public:
//	void Run();
//	Application();
//	
//private:
//	void ShowWindows();
//
//private:
//	//TApplication app = TApplication("app", nullptr, nullptr);
//
//};

