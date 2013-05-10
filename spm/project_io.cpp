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
	const string ProjectIO::error_incomplete_file("Error: reading failed. incomplete or invalid file.");

	ProjectIO::ProjectIO(string prj_config_file)
	{
		has_density_input_ = false;
		load_project(prj_config_file);
		gamma = 1.0;
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
			trim(value);
			attr_val[name] = value;
		}
		//debug_print();
		check_project_validity();
	}
	
	void ProjectIO::trim(std::string& s)
	{
		string::iterator it;
		for ( it = s.begin(); it != s.end() && std::isspace(*it); ++it);
		s.erase(s.begin(), it);
		string::reverse_iterator rit;
		for ( rit = s.rbegin(); rit != s.rend() && std::isspace(*rit); ++rit);
		s.erase(rit.base(), s.end());	
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
	void ProjectIO::check_project_validity()
	{
		if (!mesh_polygon_coupled())
		{
			if (!texture_specified()) // polygons are in one single file, without texture
			{
				TagLookupTable::const_iterator it = attr_val.find("PolygonFile");
				if ( it == attr_val.end() )
					prompt_and_exit(error_pgn_file_fail + " Not specified.");
				single_tile_file = it->second;
			}
			else // polygons are in one directory, with texture
			{
				TagLookupTable::const_iterator it = attr_val.find("PolygonDir");
				if ( it == attr_val.end() )
					prompt_and_exit(error_pgn_file_fail + " Polygon directory not specified. ");
				single_tile_dir = it->second;
			}
			if (!multi_meshes_specified()) // one single mesh, one single file
			{
				TagLookupTable::const_iterator it = attr_val.find("MeshFile");
				if ( it == attr_val.end() )
					prompt_and_exit(error_mesh_file_fail+" Not specified.");
				single_mesh_file = it->second;
				it = attr_val.find("DensityFile");
				if (it != attr_val.end())
					single_density_file = it->second;
			}
			else // multiple meshes, multiple files in one single directory
			{
				TagLookupTable::const_iterator it = attr_val.find("MeshDir");
				if (it == attr_val.end() )
					prompt_and_exit("Mesh Directory not specified.");
				std::vector<std::string> all_files;
				Geex::FileSystem::get_files(it->second, all_files);
				if (all_files.empty())
					prompt_and_exit("Fatal error: program cannot load multiple meshes.");
				CaseInsensitiveTagCmp cmp;
				for (unsigned int i = 0; i < all_files.size(); i++)
				{
					std::string ext = Geex::FileSystem::extension(all_files[i]);
					if ( !cmp(ext, "obj") && !cmp("obj", ext) )
						multi_mesh_files.push_back(all_files[i]);
					else if ( !cmp(ext, "txt") && !cmp("txt", ext) )
						multi_density_files.push_back(all_files[i]);
				}
			}
		}
		else // one polygon set (with texture) is paired with one mesh inside one directory, several such pairs
		{
			TagLookupTable::const_iterator it = attr_val.find("MeshPolygonCoupleDir");
			if (it == attr_val.end())
				prompt_and_exit("No mesh-polygon pair (<MeshPolygonCoupleDir>) input specified!");
			// extract all the directories
			std::string::const_iterator semicolon_it,  copy_start_it = it->second.begin();
			mesh_tile_couple_dir.push_back(std::string());
			while ( (semicolon_it = std::find(copy_start_it, it->second.end(), ';')) != it->second.end() )
			{
				std::copy(copy_start_it, semicolon_it, std::back_insert_iterator<std::string>(mesh_tile_couple_dir.back()));
				trim(mesh_tile_couple_dir.back());
				copy_start_it = semicolon_it + 1;
				mesh_tile_couple_dir.push_back(std::string());
			}
			std::copy(copy_start_it, semicolon_it, std::back_insert_iterator<std::string>(mesh_tile_couple_dir.back()));
			trim(mesh_tile_couple_dir.back());
			// debug print
			std::cout<<"Mesh-polygon pair directories are: \n";
			for (unsigned int i = 0; i < mesh_tile_couple_dir.size(); i++)
				std::cout<<'|'<<mesh_tile_couple_dir[i]<<'|'<<std::endl;

			// get the file of all meshes
			for (unsigned int i = 0; i < mesh_tile_couple_dir.size(); i++)
			{
				std::vector<std::string> all_files;
				Geex::FileSystem::get_files(mesh_tile_couple_dir[i], all_files);
				if (all_files.empty())
					prompt_and_exit("Fatal error: program cannot load multiple meshes at the directory." + mesh_tile_couple_dir[i]);
				CaseInsensitiveTagCmp cmp;
				for (unsigned int i = 0; i < all_files.size(); i++)
				{
					std::string ext = Geex::FileSystem::extension(all_files[i]);
					if ( !cmp(ext, "obj") && !cmp("obj", ext) )
						multi_mesh_files.push_back(all_files[i]);
					else if ( !cmp(ext, "txt") && !cmp("txt", ext) )
						multi_density_files.push_back(all_files[i]);
				}
			}
		}
	}
	void ProjectIO::debug_print()
	{
		for (TagLookupTable::const_iterator it = attr_val.begin(); it != attr_val.end(); it++)
			std::cout<<it->first<<": |"<<it->second<<'|'<<std::endl;
	}

	void ProjectIO::read_polygons_from_file(std::vector<Ex_polygon_2>& res, const std::string& fn)
	{
		std::ifstream infile(fn.c_str());
		if (infile.fail())
			prompt_and_exit(error_pgn_file_fail);
		res.clear();
		unsigned int nb_polygons = 0;
		infile >> nb_polygons;
		res.reserve(nb_polygons);
		for (unsigned int i = 0; i < nb_polygons; i++)
		{	
			unsigned int nb_edges = 0;
			infile >> nb_edges;
			res.push_back(Ex_polygon_2());
			for (unsigned int j = 0; j < nb_edges; j++)
			{
				double x, y; 
				infile >> x;
				infile >> y;
				res.back().push_back(Point_2(x, y));
			}
		}
	}

	void ProjectIO::read_polygons_from_dir(std::vector<Ex_polygon_2>& res, const std::string& dir)
	{
		std::vector<std::string> all_files;
		Geex::FileSystem::get_files(dir, all_files);
		if (all_files.empty())
			prompt_and_exit(error_pgn_file_fail + " Directory " + dir +" does not exist.");
		std::vector<std::string> all_pgn_files;
		for (unsigned int i = 0; i < all_files.size(); i++)
			if (Geex::FileSystem::extension(all_files[i]) == "tply")
				all_pgn_files.push_back(all_files[i]);
		res.reserve(all_pgn_files.size());
		int texture_id = 0;
		for (unsigned int i = 0; i < all_pgn_files.size(); i++)
		{
			//std::cout<<"loading "<<all_pgn_files[i]<<std::endl;
			std::ifstream pgn_file(all_pgn_files[i].c_str());
			unsigned int nb_vert;
			pgn_file >> nb_vert;
			res.push_back(Ex_polygon_2());
			for (unsigned int j = 0; j < nb_vert; j++)
			{
				double x, y;
				double tx, ty;
				if (!(pgn_file >> x) )
					prompt_and_exit(error_incomplete_file); 
				if (!(pgn_file >> y) )
					prompt_and_exit(error_incomplete_file);
				res.back().push_back(Point_2(x, y));

				if (!(pgn_file >> tx))
					prompt_and_exit(error_incomplete_file);
				if (!(pgn_file >> ty))
					prompt_and_exit(error_incomplete_file);
				res.back().texture_coords.push_back(Point_2(tx, ty));
			}
			res.back().texture_id = texture_id;
			texture_id++;
		}
	}
	ProjectIO& ProjectIO::operator>>(vector<Ex_polygon_2>& polygons)
	{
		if (!texture_specified())
		{
			read_polygons_from_file(polygons, single_tile_file);
		}
		else // with texture
		{
			read_polygons_from_dir(polygons, single_tile_dir);
		}

		return *this;
	}

	ProjectIO& ProjectIO::operator>>(std::vector<std::vector<Ex_polygon_2>>& multi_polygon_set)
	{
		multi_polygon_set.reserve(mesh_tile_couple_dir.size());
		for (unsigned int i = 0; i < mesh_tile_couple_dir.size(); i++)
		{
			multi_polygon_set.push_back(std::vector<Ex_polygon_2>());
			read_polygons_from_dir(multi_polygon_set.back(), mesh_tile_couple_dir[i]);
		}
		return *this;
	}

	void ProjectIO::read_texture_files(std::vector<string>& texture_files)
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
		{
			if (Geex::FileSystem::extension(all_files[i]) == "jpg")
				texture_files.push_back(all_files[i]);
			else if (Geex::FileSystem::extension(all_files[i]) == "png")
				texture_files.push_back(all_files[i]);
		}
	}

	void ProjectIO::read_texture_files(std::vector<std::vector<std::string>>& texture_files_set)
	{
		std::for_each(texture_files_set.begin(), texture_files_set.end(), std::mem_fun_ref(&std::vector<std::string>::clear));
		texture_files_set.clear();
		texture_files_set.reserve(mesh_tile_couple_dir.size());
		for (unsigned int i = 0; i < mesh_tile_couple_dir.size(); i++)
		{
			
			std::vector<std::string> all_files;
			Geex::FileSystem::get_files(mesh_tile_couple_dir[i], all_files);
			if (all_files.empty())
				prompt_and_exit("The directory " + mesh_tile_couple_dir[i] + " does not exist. ");
			texture_files_set.push_back(std::vector<std::string>());
			for (unsigned int j = 0; j < all_files.size(); j++)
			{
				std::string ext = Geex::FileSystem::extension(all_files[j]);
				CaseInsensitiveTagCmp cmp;
				if ( !cmp(ext, "jpg") && !cmp("jpg", ext) )
					texture_files_set.back().push_back(all_files[j]);
				else if ( !cmp(ext, "png") && !cmp("png", ext) )
					texture_files_set.back().push_back(all_files[j]);
			}		
		}
	}
	
	ProjectIO& ProjectIO::operator>>(TriMesh& mesh)
	{
		if (single_mesh_file.empty())
			prompt_and_exit("No single mesh file input. Try multiple mesh files.");

		mesh.load(single_mesh_file);

		// possibly load density(curvature) file
		TagLookupTable::const_iterator gamma_it = attr_val.find("gamma");
		std::istringstream str_gamma;
		if (gamma_it != attr_val.end())
			str_gamma.str(gamma_it->second);
		else
			str_gamma.str("1.0");
		str_gamma >> gamma;
		std::cout<<"using gamma = "<<str_gamma.str()<<std::endl;
		has_density_input_ = mesh.load_density(single_density_file, gamma);
	
		return *this;
	}

	ProjectIO& ProjectIO::operator>>(std::vector<TriMesh>& multimesh)
	{
		if (multi_mesh_files.empty())
			prompt_and_exit("Mesh Directory not specified.");

		for (unsigned int i = 0; i < multi_mesh_files.size(); i++)
			//if (Geex::FileSystem::extension(multi_mesh_files[i]) == "obj")
			{
				multimesh.push_back(TriMesh());
				multimesh.back().load(multi_mesh_files[i]);
			}
		// load density
		TagLookupTable::const_iterator gamma_it = attr_val.find("gamma");
		if ( gamma_it != attr_val.end() )
		{
			std::istringstream str_gamma;
			str_gamma.str(gamma_it->second);
			str_gamma >> gamma;
			std::cout<<"using gamma = "<<gamma<<std::endl;
		}
		
		unsigned int idx = 0;
		bool load_density_ok = true;
		for (unsigned int i = 0; i < multi_density_files.size(); i++)
			//if (Geex::FileSystem::extension(all_mesh_files[i]) == "txt")
			{
				load_density_ok = load_density_ok && multimesh[idx].load_density(multi_density_files[i], gamma);
				idx++;
			}
		if (load_density_ok)
			has_density_input_ = true;
		else
			has_density_input_ = false;
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

	bool ProjectIO::read_physical_scales(double& min_scale, double& max_scale, int& levels)
	{
		TagLookupTable::const_iterator it = attr_val.find("PhysicalScales");
		if (it == attr_val.end())
			return false;
		std::istringstream is;
		is.str(it->second);
		is >> min_scale;
		is >> max_scale;
		is >> levels;
		return true;
	}

	bool ProjectIO::read_optimize_scales(double& min_scale, double& max_scale, int& levels)
	{
		TagLookupTable::const_iterator it = attr_val.find("OptimizationScales");
		if (it == attr_val.end())
			return false;
		std::istringstream is;
		is.str(it->second);
		is >> min_scale;
		is >> max_scale;
		is >> levels;
		return true;
	}
}