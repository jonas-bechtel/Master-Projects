#include <iostream>
#include <type_traits>

#include "Parameter.h"

ParameterValue::ParameterValue(Value value, std::string unit)
	: m_value(value), m_unit(unit)
{
}

ParameterValue::operator int() const
{
	if (std::holds_alternative<int>(m_value))
	{
		return std::get<int>(m_value);
	}
	throw("wrong conversion");
}

ParameterValue::operator bool() const
{
	if (std::holds_alternative<bool>(m_value))
	{
		return std::get<bool>(m_value);
	}
	throw("wrong conversion");
}

ParameterValue::operator double() const
{
	if (std::holds_alternative<double>(m_value))
	{
		return std::get<double>(m_value);
	}
	throw("wrong conversion");
}

ParameterValue::operator std::filesystem::path() const
{
	if (std::holds_alternative<std::filesystem::path>(m_value))
	{
		return std::get<std::filesystem::path>(m_value);
	}
	throw("wrong conversion");
}

Parameters::Parameters(std::string description)
	: m_description(description)
{
}

void Parameters::Add(std::string name, Value value, std::string unit)
{
	m_map[name] = ParameterValue(value, unit);
}

//Value& Parameters::Get(std::string name)
//{
//	if (m_map.find(name) != m_map.end())
//	{
//		return m_map.at(name).m_value;
//	}
//	else
//	{
//		std::cout << name << " is not a Parameter." << std::endl;
//	}
//}

ParameterValue& Parameters::Get(std::string name)
{
	return m_map.at(name);
}

void Parameters::Set(std::string name, Value value)
{
	if (m_map.find(name) != m_map.end())
	{
		m_map.at(name).m_value = value;
	}
	else
	{
		std::cout << name << " is not a Parameter." << std::endl;
	}
}
