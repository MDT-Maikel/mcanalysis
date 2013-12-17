/* Utility functions
 *
 * Provides useful functions for file handling
*/

#include "utility.h"


/* NAMESPACE */
namespace analysis
{
	
	// returns a vector of paths pointing to all files in a directory matching the pattern
	std::vector<boost::filesystem::path> get_files(std::string dir, std::string type_pattern, bool recursive)
	{
		// get all files in dir		
		boost::filesystem::path directory(dir);
		std::vector<boost::filesystem::path> files;
		copy(boost::filesystem::directory_iterator(directory), boost::filesystem::directory_iterator(), std::back_inserter(files));
		
		// also look into subfolders if recursive is true
		if (recursive)
		{
			unsigned int index = 0;
			while (index < files.size())
			{
				boost::filesystem::path dir = files[index];
				if (is_directory(dir))
				{
					// Add new files from this directory to the list and remove that directory.				
					std::vector<boost::filesystem::path> new_files;
					copy(boost::filesystem::directory_iterator(dir), boost::filesystem::directory_iterator(), std::back_inserter(new_files));
					files.erase(files.begin() + index);
					for (unsigned int i = 0; i < new_files.size(); i++)
						files.push_back(new_files[i]);
				}
				else
					index++;
			}
		}

		// remove all the files which do not match the file pattern
		for (int i = files.size() - 1; i >= 0 ; i--)
		{	
			std::string file = files[i].string();
			std::regex to_match("(.*)(" + type_pattern + ")");
			if (!std::regex_match(file, to_match))	
				files.erase(files.begin() + i);
		}

		// return the list of paths
		return files;
	}

	// reads lhco events into a vector of events
	void read_lhco(std::vector<event*> & events, boost::filesystem::path file)
	{
		// open the file with gzstream	
		std::string file_name = file.string();
		igzstream file_igz;
		file_igz.open(file_name.c_str(), std::ios::in);
		
		// variable to dump useless stuff into
		std::string dump;

		// find and remove header if existing
		bool found_header = false;
		while (file_igz)
		{
			file_igz >> dump;
			// locate end of header and remove it from the stream
			if (dump == "#</LesHouchesEvents>")
			{
				// remove additional information before event list
				file_igz >> dump >> dump >> dump >> dump >> dump >> dump;
				found_header = true;
				break;
			}
		}

		// reset input stream if header was not found
		igzstream file_igz_nohead;
		if (!found_header)
			file_igz_nohead.open(file_name.c_str(), std::ios::in);

		// loop over all events
		int num; 
		while ((found_header && file_igz >> num) || (!found_header && file_igz_nohead >> num))
		{
			// first entry already read, check if it is an event header
			if (num == 0)
			{
				if (found_header)
					file_igz >> dump >> dump;
				else
					file_igz_nohead >> dump >> dump;				
				continue;			
			}

			// it is not a header: store lhco into vector of events.
			event *ev = new event;
			while (true)
			{
				// read in lhco particle				
				lhco *p = new lhco;
				if (found_header)
					p->read(file_igz);
				else
					p->read(file_igz_nohead);

				ev->push_back(p);
				//p->write(std::cout);
				if (p->type() == particle::type_met)
				{
					break;
				}
				if (found_header)
					file_igz >> num;
				else
					file_igz_nohead >> num;
			}
			
			// push back into the event vector
			events.push_back(ev);
		}
	}

	// write lhco events into a file
	void write_lhco(const std::vector<event*> & events, boost::filesystem::path file)
	{
		// open the file with gzstream	
		std::string file_name = file.string();
		ogzstream file_ogz;
		file_ogz.open(file_name.c_str());
		
		// loop over all events
		for (unsigned int index = 0; index < events.size(); index++)
		{
			event *ev = events[index];
			
			// print the event header
			file_ogz << 0 << "\t" << index + 1 << "\t" << 0 << std::endl;
			
			// loop over all particles
			for (unsigned int i = 0; i < ev->size(); i++)
			{
				file_ogz << i + 1 << "\t";
				particle *p = (*ev)[i];
				p->write(file_ogz);
			}
		}
	}

	// reads lhe events into a vector of events
	void read_lhe(std::vector<event*> & events, boost::filesystem::path file)
	{
		// open the file with gzstream	
		std::string file_name = file.string();
		igzstream file_igz;
		file_igz.open(file_name.c_str(), std::ios::in);
		
		// variable to dump useless stuff into
		std::string dump;
		
		// loop over all events with a while loop 
		// and search for the <event> tag
		while (file_igz >> dump)
		{
			// found the start of an event
			if (dump == "<event>")
			{
				// get the number of particles in the event
				unsigned int nr_particles;
				file_igz >> nr_particles;
				
				// dump the five remaining useless variables
				file_igz >> dump >> dump >> dump >> dump >> dump;
				
				// make a new event and fill it
				event *ev = new event;
				for (unsigned int i = 0; i < nr_particles; i++)
				{
					lhe *p = new lhe;
					p->read(file_igz);
					ev->push_back(p);
				}

				// push back into the event vector
				events.push_back(ev);						
			}			
		}
	}

	// write lhe events into a file
	void write_lhe(const std::vector<event*> & events, boost::filesystem::path file)
	{
		// open the file with gzstream	
		std::string file_name = file.string();
		ogzstream file_ogz;
		file_ogz.open(file_name.c_str());
		
		// loop over all events
		for (unsigned int index = 0; index < events.size(); index++)
		{
			event *ev = events[index];
			
			// print the event header
			file_ogz << "<event>" << std::endl;
			file_ogz << ev->size() << "\t";
			file_ogz << "0 \t 0 \t 0 \t 0 \t 0" << std::endl;
			
			// loop over all particles
			for (unsigned int i = 0; i < ev->size(); i++)
			{
				particle *p = (*ev)[i];
				p->write(file_ogz);
			}
			file_ogz << "</event>" << std::endl;
		}
	}

/* NAMESPACE */
}
