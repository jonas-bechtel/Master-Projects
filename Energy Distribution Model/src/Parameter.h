#pragma once

#include <string>
#include <unordered_map>
#include <filesystem>
#include <variant> 

using Value = std::variant<int, double, std::filesystem::path, bool>;

struct ParameterValue
{
	std::string m_unit = "";
	Value m_value;

	ParameterValue() {};
	ParameterValue(Value value, std::string unit);

	operator int() const;
	operator bool() const;
	operator double() const;
	operator std::filesystem::path() const;
};

class Parameters
{
public:
	enum Params {};

	Parameters(std::string description);
	void Add(std::string name, Value value, std::string unit = "");
	//Value& Get(std::string name);

	ParameterValue& Get(int index);

	void Set(std::string name, Value value);

	std::string ToString();
	void FromString(std::string string);

private:
	std::string m_description;
	std::vector<ParameterValue> m_params;
	//std::unordered_map<std::string, ParameterValue> m_map;
};

class LabEnergyParameter : public Parameters
{
public:
	enum Params {CenterLabEnergy, EnergyFile, Count};

	//LabEnergyParameter(std::string description);
};


//------------------------------------------------------
#include <iostream>
#include <vector>
#include <variant>
#include <filesystem>
#include <typeinfo>

using Path = std::filesystem::path;
//using Value = std::variant<int, double, bool, Path>;

template<typename T>
struct ParameterValue
{
	union
	{
		int i;
		double d;
		bool b;
	};

	//T value;
	std::string string;

	ParameterValue(T val, const std::string& str)
		: value(val), string(str)
	{
	}

	T get()
	{
		return value;
	}
};

struct LabEnergyParameter
{
	ParameterValue<int> param1 = ParameterValue(5, "bla1");
	ParameterValue<double> param2 = ParameterValue(1.23, "bla2");
	ParameterValue<bool> param3 = ParameterValue(true, "bla3");

	LabEnergyParameter()
	{

	}
};

struct Parameter
{
	//template<typename T>
	//int& GetInt(int index)
	//{
	//	if (index < list.size())
	//	{
	//		if (std::holds_alternative<int>(list[index].value))
	//		{
	//			return std::get<int>(list[index].value);
	//		}
	//	}
	//	std::cerr << "could not get index " << index << " as an integer" << std::endl;
	//	int temp = 0;
	//	return temp;
	//}
	//
	//double& GetDouble(int index)
	//{
	//	if (std::holds_alternative<double>(list[index].value))
	//	{
	//		return std::get<double>(list[index].value);
	//	}
	//	std::cerr << "could not get index " << index << " as a double" << std::endl;
	//	double temp = 0;
	//	return temp;
	//}

protected:
	//std::vector<ParameterValue> list;
};

enum LabEnergyParameterNames
{
	Param1, Param2, NumberLabEnergyParameter
};

enum ElectronBeamParameterNames
{
	Current, Radius, NumberElectronBeamParameter
};


int main()
{

	LabEnergyParameter params;

	//ParameterValue<int> intParam(5, "Integer Parameter");


	std::cout << sizeof(ParameterValue<int>) << " " << sizeof(ParameterValue<Path>) << std::endl;

	//params.param1.i = 3;


	//int& param1 = params.GetInt(Param1);
	//
	std::cout << params.param1.get() << std::endl;
	//
	//param1 = 10;
	//
	//std::cout << params.GetDouble(Param1) << std::endl;
	//std::cout << params.GetDouble(Param2) << std::endl;

	return 1;
}

