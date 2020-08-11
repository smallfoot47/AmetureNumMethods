#pragma once
#include <iostream>
#include <typeinfo>

template <class T_t>
class only //hard type locking container
{
	T_t data;

public:
	
	only(const T_t& i_data)
	{
		data = i_data;
	}

	template <class L_t>
	only(const L_t& i_blocked_type)
	{
		try {
				std::string this_type = typeid(data).name();
				std::string that_type = typeid(i_blocked_type).name();

				std::string message = "ERR: ATTEMTPTING TO CHANGE only<" + this_type;
				message += "> to ";
				message += that_type + "(Illegal operation)\n";
				message += "SOL: Change the input to the function(s) that accept an only<";
				message += this_type + "> type or typecast ";
				message += that_type + " to ";
				message += this_type + "\n";

				throw std::runtime_error(message.c_str());
		}

		catch (std::exception e)
		{
			std::cout << e.what();
			exit(-69);
		}

	}

	inline operator T_t&()
	{
		return data;
	}
};

template <class D_t>
class flagged //weak access lock container
{
	D_t data;
	bool flag;

public:

	flagged(const D_t& i_data)
	{
		data = i_data;
	}

	flagged(const bool i_flag)
	{
		flag = i_flag;
	}

	operator D_t&()
	{
		return data;
	}

	operator bool&()
	{
		return flag;
	}
};

template <class D_t>
class locked //strong access lock container
{
	D_t data;
	bool flag;

public:

	locked(const D_t& i_data)
	{
		if (flag)
			data = i_data;
	}

	locked(const bool i_flag)
	{
		flag = i_flag;
	}

	operator D_t&()
	{
		if (flag)
			return data;

		return D_t(data);
	}

	operator bool&()
	{
		return flag;
	}
};

template <class T_1, class T_2 = T_1>
class dual //unique access locked heterotype container
{
	T_1 primary_state;
	T_2 secondary_state;
	bool state = true;

public:

	uint8_t SwitchState()
	{
		
		return (state = !state) + 1u;
	}

	uint8_t State()
	{
		return state + 1u;
	}

	void State(bool i_state)
	{
		state = i_state;
	}

	template <class T>
	dual(T& in_state)//setter
	{
		if (state)
			primary_state = (T_1)in_state;
		else
			secondary_state = (T_2)in_state;
	}

	template <class T>
	operator T&()//getter
	{
		if (state)
			return (T&)primary_state;
		return (T&)secondary_state;
	}
};