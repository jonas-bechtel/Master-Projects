#include "pch.h"

#include "RandomFunctions.h"

ImVec4 ModifyColor(const ImVec4& color, float factor) 
{
    ImVec4 newColor = { 0, 0, 0, 1 };
    newColor.x = std::clamp(color.x * factor, 0.0f, 1.0f);
    newColor.y = std::clamp(color.y * factor, 0.0f, 1.0f);
    newColor.z = std::clamp(color.z * factor, 0.0f, 1.0f);
    //std::cout << factor << std::endl;
    //std::cout << newColor.x << ", " << newColor.y << ", " << newColor.z << std::endl;

    return newColor;

}