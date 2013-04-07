#include "project_io.h"

namespace Geex
{
	const string ProjectIO::error_start_tag("Error: Bad format of project configuration file! Incorrect start tag.");
	const string ProjectIO::error_end_tag("Error: Bad format of project configuration file! Incorrect end tag.");
	const string ProjectIO::error_tag_dismatch("Error: Bad format of project configuration file! Tags do not match.");
	const string ProjectIO::error_prj_file_fail("Error: Failure to open project configuration file!");
	const string ProjectIO::error_pgn_file_fail("Error: Failure to open polygon file!");
	const string ProjectIO::error_mesh_file_fail("Error: Failure to open mesh file!");
	const string ProjectIO::error_texture_input("Error: No texture image input!\n");

	ProjectIO::ProjectIO(string prj_config_file)
	{
		has_density_input_ = false;
		load_project(prj_config_file);
	}

	void ProjectIO::load_project(const string& prj_config_file)
	{
		ifstream prj_file(prj_config_file.c_str(), std::ios_base::in);

		if (!prj_file)
			prompt_and_exit(error_prj_file_fail);

		attr_val["ProjectConfigFile"] = prj_config_file; // project configuration file path

		while (!prj_file.eof())
		{
			string name = get_start_token(prj_file);
			if (name.empty())
				break;
			string value("");
			char c;
			while ( prj_file.get(c) && c != '<'  ) // read the attribute value
				value += c;
			if (prj_file.eof())
				prompt_and_exit(error_end_tag);

			prj_file.putback(c);
			string end_name = get_end_token(prj_file);
			CaseInsensitiveTagCmp tag_match;
			if (tag_match(name, end_name) || tag_match(end_name, name))
				prompt_and_exit(error_tag_dismatch);

			// trim
			string::iterator it;
			for ( it = value.begin(); it != value.end() && std::isspace(*it); ++it);
				value.erase(value.begin(), it);
			string::reverse_iterator rit;
			for ( rit = value.rbegin(); rit != value.rend() && std::isspace(*rit); ++rit);
			value.erase(rit.base(), value.end());
			attr_val[name] = value;
		}
	}
	string ProjectIO::get_start_token(std::ifstream& fs)
	{
		string token;
		char c;
		if ( !(fs>>c) )
			return "";// no tag any more
		if ( c != '<' )
			prompt_and_exit(error_start_tag);
		while ( fs >> c && c != '>' )
			token += c;
		if ( c != '>' )
			prompt_and_exit(error_start_tag);
		return token;
	}
	string ProjectIO::get_end_token(std::ifstream& fs)
	{
		string token;
		char c;
		if ( !(fs >> c) || c != '<' )
			prompt_and_exit(error_end_tag);
		if ( !(fs >> c) || c != '/')
			prompt_and_exit(error_end_tag);
		while ( fs >> c && c != '>' )
			token += c;
		if ( c != '>' )
			prompt_and_exit(error_end_tag);
		return token;
	}

	void ProjectIO::debug_print()
	{
		for (TagLookupTable::const_iterator it = attr_val.begin(); it != attr_val.end(); it++)
			std::cout<<it->first<<": |"<<it->second<<'|'<<std::endl;
	}

	ProjectIO& ProjectIO::operator>>(vector<Polygon_2>& polygons)
	{
		TagLookupTable::const_iterator it = attr_val.find("PolygonFile");
		if ( it == attr_val.end() )
			prompt_and_exit(error_pgn_file_fail + " Not specified.");
		std::ifstream pgn_file(it->second.c_str());
		if (pgn_file.fail())
			prompt_and_exit(error_pgn_file_fail);
		polygons.clear();
		unsigned int nb_polygons = 0;
		pgn_file >> nb_polygons;
		polygons.reserve(nb_polygons);
		for (unsigned int i = 0; i < nb_polygons; i++)
		{	
			unsigned int nb_edges = 0;
			pgn_file >> nb_edges;
			polygons.push_back(Polygon_2());
			for (unsigned int j = 0; j < nb_edges; j++)
			{
				double x, y; 
				pgn_file >> x;
				pgn_file >> y;
				polygons.back().push_back(Point_2(x, y));
			}
		}

		return *this;
	}

	ProjectIO& ProjectIO::operator>>(std::vector<string>& texture_files)
	{
		texture_files.clear();
		TagLookupTable::const_iterator it = attr_val.find("TextureImageDir");
		if ( it == attr_val.end() )
		{
			std::cout<<error_texture_input<<std::endl;
			return;
		}
		std::vector<std::string> all_files;
		Geex::FileSystem::get_files(it->second, all_files);
		if (all_files.empty())
		{
			std::cout << "Error: The input directory \""<<it->second<<"\" does not exist!\n";
			return;
		}
		// filter out image file
		for (unsigned int i = 0; i < all_files.size(); i++)
			if (Geex::FileSystem::extension(all_files[i]) == "jpg")
				texture_files.push_back(all_files[i]);
		return *this;
	}
	ProjectIO& ProjectIO::operator>>(TriMesh& mesh)
	{
		TagLookupTable::const_iterator it = attr_val.find("MeshFile");
		if ( it == attr_val.end() )
			prompt_and_exit(error_mesh_file_fail+" Not specified.");
		mesh.load(it->second);

		// possibly load density(curvature) file
		it = attr_val.find("DensityFile");
		if ( it != attr_val.end() )
			has_density_input_ = mesh.load_density(it->second, 1.2);

		return *this;
	}

	string ProjectIO::attribute_value(const string& attribute_name)
	{
		TagLookupTable::const_iterator it = attr_val.find(attribute_name);
		if ( it == attr_val.end() )
			return "";
		else
			return it->second;
	}
}