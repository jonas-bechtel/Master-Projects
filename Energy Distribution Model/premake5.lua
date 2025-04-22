ROOT_DIR = os.getenv("ROOT_DIR")
if ROOT_DIR == nil then
    print("ROOT_DIR environment variable is not set! Please make sure it is configured.")
end

project "Energy Distribution Model"
	kind "ConsoleApp"
	language "C++"
	cppdialect "C++17"
	staticruntime "off"

	targetdir ("%{wks.location}/bin/" .. outputdir .. "/%{prj.name}")
	objdir ("%{wks.location}/bin-int/" .. outputdir .. "/%{prj.name}")

	pchheader "pch.h"
	pchsource "src/pch.cpp"

	files	
	{
		"src/**.h",
		"src/**.cpp"
	}

	includedirs
	{
		"%{wks.location}/Energy Distribution Model/src",
		"%{wks.location}/Energy Distribution Model/src/Cooling Force",
		"%{wks.location}/Energy Distribution Model/src/Cross Section Deconvolution",
		"%{wks.location}/Energy Distribution Model/src/Energy Distribution Generation",
		"%{wks.location}/Energy Distribution Model/src/general",
		"%{wks.location}/Energy Distribution Model/src/utility",
		"%{wks.location}/Energy Distribution Model/vendor/imgui",
		"%{wks.location}/Energy Distribution Model/vendor/imgui/backends",
		"%{wks.location}/Energy Distribution Model/vendor/implot",
		"%{wks.location}/Energy Distribution Model/vendor/eigen",
		"%{wks.location}/Energy Distribution Model/vendor/tinyfiledialogs",
		"%{wks.location}/Energy Distribution Model/vendor/JSPEC/include",
		"%{ROOT_DIR}/include"
	}

	libdirs 
	{
		"%{ROOT_DIR}/lib"
	}

	links
	{
		"imgui",
		"implot",
		"tinyfiledialogs",
		"JSPEC",
		"Betacool",
		"libCore",
		"libRIO",
		"libHist",
		"libGpad",
		"libGraf",
		"libGraf3d",
		"libMatrix",
		"libMathCore",
		"libPhysics"
	}

	filter "system:windows"
		systemversion "latest"

		links
		{
			"d3d12",
			"d3dcompiler",
			"dxgi"
		}

		postbuildcommands
		{
			'cmd /c if exist vendor\\JSPEC\\lib\\*.dll xcopy /Q /Y /I vendor\\JSPEC\\lib\\*.dll "%{cfg.targetdir}" > nul'
		}


	filter "configurations:Debug"
		defines "_DEBUG"
		runtime "Debug"
		symbols "on"

	filter "configurations:Release"
		defines "NDEBUG"
		runtime "Release"
		optimize "on"
