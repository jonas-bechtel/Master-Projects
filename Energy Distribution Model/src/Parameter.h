#pragma once

#include <string>
#include <iostream>
#include <sstream>
#include <filesystem>

using Path = std::filesystem::path;

struct float2
{
	float x, y;

	float2(float x_, float y_) : x(x_), y(y_) {}
	float2(std::string& str)
	{
		std::stringstream ss(str);
		std::string number;

		std::getline(ss, number, ',');
		x = std::stof(number);
		std::getline(ss, number, ',');
		y = std::stof(number);
	}
};

struct float3
{
	float x, y, z;

	float3(float x_, float y_, float z_) : x(x_), y(y_), z(z_) {}
	float3(std::string& str)
	{
		std::stringstream ss(str);
		std::string number;

		std::getline(ss, number, ',');
		x = std::stof(number);
		std::getline(ss, number, ',');
		y = std::stof(number);
		std::getline(ss, number, ',');
		z = std::stof(number);
	}
};

template<typename T>
struct ParameterValue
{
	enum Type { INT, DOUBLE, FLOAT_2, FLOAT_3, BOOL, PATH };

	// needs to be the size of the largest possible stored element
	using StorageType = std::aligned_storage_t<sizeof(Path), alignof(Path)>;

	int type = INT;
	StorageType value;
	std::string name;
	std::string format;

	ParameterValue(T val, const std::string& name, const std::string& format)
		: name(name), format(format)
	{
		new (&value) T(val);

		const std::type_info& typeInfo = typeid(T);
		if (typeInfo == typeid(int))
		{
			type = INT;
		}
		else if (typeInfo == typeid(double))
		{
			type = DOUBLE;
		}
		else if (typeInfo == typeid(float2))
		{
			type = FLOAT_2;
		}
		else if (typeInfo == typeid(float3))
		{
			type = FLOAT_3;
		}
		else if (typeInfo == typeid(bool))
		{
			type = BOOL;
		}
		else if (typeInfo == typeid(Path))
		{
			type = PATH;
		}
	}

	T& get()
	{
		return *reinterpret_cast<T*>(&value);
	}
	T* data()
	{
		return (T*)&value;
	}
	void set(T val)
	{
		reinterpret_cast<T*>(&value)->~T();

		// Construct the new value in the same storage
		new (&value) T(val);
	}

	operator T() 
	{
		return *(T*)(&value);
	}
	operator T* () 
	{
		return (T*)&value;
	}
	operator float* ()
	{
		return (float*)&value;
	}
	ParameterValue& operator=(const T& val)
	{
		set(val);
		return *this;
	}
};

struct Parameters
{
public:
	std::string toString()
	{
		std::string result = "# " + m_name + ":\n";
		char* pointer = (char*)this;

		for (int offset = m_dataStart; offset < GetSize(); offset += sizeof(ParameterValue<int>))
		{
			int type = *(int*)(pointer + offset) + offsetof(ParameterValue<int>, type);
			std::string name = *(std::string*)((pointer + offset) + offsetof(ParameterValue<int>, name));
			std::string format = *(std::string*)((pointer + offset) + offsetof(ParameterValue<int>, format));

			void* value = ((pointer + offset) + offsetof(ParameterValue<int>, value));

			result += "# " + name + ": ";
			switch (type)
			{
			case ParameterValue<int>::INT:
				result += createLine(format, *(int*)value);
				break;
			case ParameterValue<int>::DOUBLE:
				result += createLine(format, *(double*)value);
				break;
			case ParameterValue<int>::FLOAT_2:
				result += createLine(format, (*(float2*)value).x, (*(float2*)value).y);
				break;
			case ParameterValue<int>::FLOAT_3:
				result += createLine(format, (*(float3*)value).x, (*(float3*)value).y, (*(float3*)value).z);
				break;
			case ParameterValue<int>::BOOL:
				result += createLine(format, *(bool*)value);
				break;
			case ParameterValue<int>::PATH:
				result += createLine(format, (*(Path*)value).filename().string().c_str());
				break;
			}
			result += "\n";
		}

		return result;
	}
	void fromString(std::string header)
	{
		std::istringstream stream(header);
		std::string line;

		while (std::getline(stream, line))
		{
			line = line.substr(2);
			std::string name, value;

			// Create a stringstream for the line
			std::stringstream lineStream(line);

			// Read the name (everything before ':')
			std::getline(lineStream, name, ':');

			// Read the value (everything after ':')
			std::getline(lineStream, value);

			std::cout << "name: " << name << " value: " << value << std::endl;

			void* object = getParameterValue(name);
			int type = *(int*)(object)+offsetof(ParameterValue<int>, type);

			switch (type)
			{
			case ParameterValue<int>::INT:
				((ParameterValue<int>*)object)->set(std::stoi(value));
				break;
			case ParameterValue<double>::DOUBLE:
				((ParameterValue<double>*)object)->set(std::stod(value));
				break;
			case ParameterValue<int>::FLOAT_2:
				((ParameterValue<float2>*)object)->set(float2(value));
				break;
			case ParameterValue<int>::FLOAT_3:
				((ParameterValue<float3>*)object)->set(float3(value));
				break;
			case ParameterValue<int>::BOOL:
				((ParameterValue<bool>*)object)->set(std::stoi(value));
				break;
			case ParameterValue<int>::PATH:
				((ParameterValue<Path>*)object)->set(Path(value));
				break;
			}
		}
	}

private:
	template <typename T, typename... Args>
	std::string createLine(const std::string& format, const T& value1, const Args&... args)
	{
		char buffer[100];
		sprintf_s(buffer, format.c_str(), value1, args...);
		return std::string(buffer);
	}

	void* getParameterValue(const std::string& name)
	{
		char* pointer = (char*)this;

		for (int offset = m_dataStart; offset < GetSize(); offset += sizeof(ParameterValue<int>))
		{
			int type = *(int*)(pointer + offset) + offsetof(ParameterValue<int>, type);
			std::string ParamValueName = *(std::string*)((pointer + offset) + offsetof(ParameterValue<int>, name));
			if (ParamValueName == name)
			{
				return pointer + offset;
			}
		}
		return nullptr;
	}
	virtual int GetSize() = 0;

protected:
	std::string m_name = "";
	void setName(std::string&& name) { m_name = std::move(name); }

private:
	// offset in memory where the data of derived classes start, comes from vtable (8 bytes atm) and members of Parameters
	int m_dataStart = 8 + offsetof(Parameters, m_dataStart);
};



