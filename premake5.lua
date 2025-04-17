workspace "Energy Distribution Model"
	architecture "x86_64"
	startproject "Energy Distribution Model"

	configurations
	{
		"Debug",
		"Release",
	}

	flags
	{
		"MultiProcessorCompile"
	}

outputdir = "%{cfg.buildcfg}-%{cfg.system}-%{cfg.architecture}"

group "Dependencies"
	include "vendor/premake"
	include "Energy Distribution Model/vendor/tinyfiledialogs"
	include "Energy Distribution Model/vendor/imgui"
	include "Energy Distribution Model/vendor/implot"
	include "Energy Distribution Model/vendor/JSPEC"
group ""

group "Core"
	include "Energy Distribution Model"
group ""
