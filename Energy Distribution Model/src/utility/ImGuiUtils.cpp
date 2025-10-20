#include "pch.h"
#include "ImGuiUtils.h"
#include <Application.h>

void ImGuiUtils::TextTooltip(std::string text)
{
	if (!Application::GetSettings().tooltipsDisabled && ImGui::BeginItemTooltip())
	{
		ImGui::PushTextWrapPos(ImGui::GetFontSize() * 15.0f);
		ImGui::Text(text.c_str());
		//ImGui::TextColored(ImVec4(0.8f, 0.3f, 0.3f, 1.0f), text.c_str());
		ImGui::PopTextWrapPos();
		ImGui::EndTooltip();
	}
}
