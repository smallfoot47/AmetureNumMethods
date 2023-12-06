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
