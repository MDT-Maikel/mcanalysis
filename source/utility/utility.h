/* Utility functions
 *
 * Provides useful functions for file handling
*/

#ifndef INC_UTILITY
#define INC_UTILITY

#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

#include <boost/filesystem.hpp>

#include "../../deps/gzstream/gzstream.h"

#include "../event/event.h"


/* NAMESPACE */
namespace analysis
{

	// returns a vector of paths pointing to all files in a directory matching the pattern
	std::vector<boost::filesystem::path> get_files(std::string dir, std::string type_pattern, bool recursive = false);

	// reads lhco events into a vector of events
	void read_lhco(std::vector<event*> & events, boost::filesystem::path file);

	// write lhco events into a file
	void write_lhco(const std::vector<event*> & events, boost::filesystem::path file);

	// read lhe events into a file
	void read_lhe(std::vector<event*> & events, boost::filesystem::path file);

	// write lhe events into a file
	void write_lhe(const std::vector<event*> & events, boost::filesystem::path file);

	// read a single setting from a file
	template <typename Type>
	Type read_settings(std::string settings_file, std::string identifier);

	// read a list of settings from a file
	template <typename Type>
	std::vector<Type> read_settings_list(std::string settings_file, std::string identifier);

/* NAMESPACE */
}

/* TEMPLATE IMPLEMENTATIONS */
#include "utility.tpp"

#endif
