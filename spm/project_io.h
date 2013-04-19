#ifndef _PROJECT_IO_H_
#define  _PROJECT_IO_H_

/************************************************************************/
/*					Interaction with file system                        */
/*@Author: Wenchao Hu							                        */
/************************************************************************/

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <algorithm>
#include <functional>
#include <cctype>
#include <Geex/basics/file_system.h>
#include "spm_cgal.h"
#include "meshes.h"
#include "packing_object.h"

namespace Geex
{
	using std::string;
	using std::ifstream;
	using std::vector;
class ProjectIO
{

public:

	ProjectIO() {} 
	ProjectIO(string prj_config_file);

	/** access function **/
	string attribute_value(const string& attribute_name) ;
	bool has_density_input() const { return has_density_input_; }

	inline bool texture_specified() const 
	{
		TagLookupTable::const_iterator it = attr_val.find("WithTexture");
		CaseInsensitiveTagCmp cmp;
		return ( it != attr_val.end() && !cmp(it->second, "True") && !cmp("True", it->second) );
	}

	inline bool multi_meshes_specified() const
	{
		TagLookupTable::const_iterator it = attr_val.find("MultiMesh");
		CaseInsensitiveTagCmp cmp;
		return ( it != attr_val.end() && !cmp(it->second, "True") && !cmp("True", it->second) );
	}

	inline bool mesh_polygon_coupled() const
	{
		TagLookupTable::const_iterator it = attr_val.find("MeshPolygonCoupled");
		CaseInsensitiveTagCmp cmp;
		return ( it != attr_val.end() && !cmp(it->second, "True") && !cmp("True", it->second) );
	}

	const std::vector<std::string>& get_submesh_files() const { return multi_mesh_files; }

	/** input **/
	// load a whole project from a file, the same function as the constructor
	void load_project(const string& prj_configure_file);
	// read all texture files in one single directory
	void read_texture_files(std::vector<std::string>& texture_files);
	void read_texture_files(std::vector<std::vector<std::string>>& texture_files_set);
	// read triangle mesh
	ProjectIO& operator>>(TriMesh& mesh);
	ProjectIO& operator>>(std::vector<TriMesh>& multimesh);
	
	ProjectIO& operator>>(vector<Ex_polygon_2>& polygons); 
	ProjectIO& operator>>(std::vector<std::vector<Ex_polygon_2>>& multi_polygon_set);
	
	/** output **/
	void dump_results(const string& filename, const string& directory="");

	unsigned int nb_sub_pack_parts() const { return mesh_tile_couple_dir.size(); }

	/** debug **/
	void debug_print();

private:
	
	string get_start_token(std::ifstream& fs);
	string get_end_token(std::ifstream& fs);
	void trim(std::string& s);

	void read_polygons_from_file(std::vector<Ex_polygon_2>& res, const std::string& fn);
	void read_polygons_from_dir(std::vector<Ex_polygon_2>& res, const std::string& dir);

	void check_project_validity(); // check the validity of the project file and construct the basic file information
	//void read_mesh_from_file(TriMesh& mesh, const std::string& fn);
	
	struct CaseInsensitiveCmp : public std::binary_function<char, char, bool>
	{
		bool operator()(const char c0, const char c1) const
		{
			return std::tolower(c0) < std::tolower(c1);
		}
	};

	struct CaseInsensitiveTagCmp : public std::binary_function<string, string, bool>
	{
		bool operator()(const string& tag0, const string& tag1) const
		{
			return std::lexicographical_compare(tag0.begin(), tag0.end(), tag1.begin(), tag1.end(), CaseInsensitiveCmp());
		}
	};

private:
	typedef std::map<string, string, CaseInsensitiveTagCmp> TagLookupTable;
	TagLookupTable attr_val; // attribute-value pairs
	std::string single_mesh_file; // one single mesh
	std::string single_density_file; // density of one single mesh
	std::string single_tile_file; // all polygons are in one single file
	std::string single_tile_dir; // all polygons are in one directory with texture
	std::vector<std::string> multi_mesh_files; // multiple meshes
	std::vector<std::string> multi_density_files; // density file of multiple meshes
	//std::vector<std::string> tile_set_directories; // polygons are in different directories with texture
	std::vector<std::string> mesh_tile_couple_dir; // one polygon set and one mesh are in the same directory
	/** indicators of what are input **/
	bool has_density_input_; // each vertex has a weight, e.g. curvature

	/** i/o error message **/
	const static string error_start_tag;
	const static string error_end_tag;
	const static string error_tag_dismatch;
	const static string error_prj_file_fail;
	const static string error_pgn_file_fail;
	const static string error_mesh_file_fail;
	const static string error_texture_input;
	const static string error_texture_directory;
	const static string error_incomplete_file;
};

}
#endif