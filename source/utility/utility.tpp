/* Utility functions
 *
 * Provides useful functions for file handling
*/


/* NAMESPACE */
namespace analysis
{

	// read a single setting from a file
	template <typename Type>
	Type read_settings(std::string settings_file, std::string identifier)
	{
		std::ifstream ifs_settings;
		ifs_settings.open(settings_file.c_str());
		std::string dump; 
		Type result;
	
		// loop over the whole settings card and locate settings
		while (ifs_settings >> dump)
		{
			// find identifier
			if (dump.compare(identifier) == 0)
			{
				ifs_settings >> dump;
				if (dump.compare("=") == 0)
				{		
					ifs_settings >> result;
					break;
				}
			}
		}
		return result;
	}

	// read a list of settings from a file
	template <typename Type>
	std::vector<Type> read_settings_list(std::string settings_file, std::string identifier)
	{
		std::ifstream ifs_settings;
		ifs_settings.open(settings_file.c_str());
		std::string dump; std::stringstream stream;
	
		// loop over the whole settings card and locate settings
		while (ifs_settings >> dump)
		{
			// find identifier
			if (dump.compare(identifier) == 0)
			{
				ifs_settings >> dump;
				if (dump.compare("=") == 0)	
				{
					ifs_settings >> dump;
					stream << dump;
					break;
				}
			}
		}

		// extract vector of type from string stream
		std::vector<Type> result;
		char c; int index = 0; std::stringstream temp; 
		while (stream.get(c))
		{		
			if (c == '{')
				continue;
			if (c == ',' || c == '}') 
			{
				// reached end of entry => read and store it				
				Type entry; temp >> entry;
				temp.str(std::string()); temp.clear();
				result.push_back(entry);
				index++;				
			}
			else
				temp << c;
		}		
		return result;
	}

/* NAMESPACE */
}
