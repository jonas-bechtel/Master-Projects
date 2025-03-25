#include "pch.h"
#include "ROOTCanvas.h"

TApplication ROOTCanvas::app = TApplication("app", nullptr, nullptr);

ROOTCanvas::ROOTCanvas(const char* name, const char* title, int width, int height)
	: TCanvas(name, title, width, height)
{
	Hide();
}

void ROOTCanvas::MakeShowHideButton()
{
	ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.0f, 0.2f, 0.5f, 1.0f));
	ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.0f, 0.3f, 0.6f, 1.0f));
	ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(0.0f, 0.2f, 0.4f, 1.0f));

	std::string label = (m_shown ? "Hide " : "Show ") + std::string(GetTitle());
	if(ImGui::Button(label.c_str()))
	{
		m_shown ? Hide() : Show();
	}
		
	ImGui::PopStyleColor(3);
}

bool ROOTCanvas::IsShown()
{
	return m_shown;
}

void ROOTCanvas::Show()
{
	if (m_shown) return;

	((TRootCanvas*)GetCanvasImp())->MapRaised();
	m_shown = true;
}

void ROOTCanvas::Hide()
{
	if (!m_shown) return;

	((TRootCanvas*)GetCanvasImp())->UnmapWindow();
	m_shown = false;
}

void ROOTCanvas::Render()
{
	cd();
	Modified();
	Update();
}
