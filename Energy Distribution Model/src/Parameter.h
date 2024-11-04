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
	Parameters(std::string description);
	void Add(std::string name, Value value, std::string unit = "");
	//Value& Get(std::string name);

	ParameterValue& Get(std::string name);

	void Set(std::string name, Value value);

	std::string ToString();
	void FromString(std::string string);

private:
	std::string m_description;
	std::unordered_map<std::string, ParameterValue> m_map;
};

class LabEnergyParameter : public Parameters
{
public:
	enum Params {CenterLabEnergy, EnergyFile, Count};

	//LabEnergyParameter(std::string description);
};

