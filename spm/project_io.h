#ifndef _PROJECT_IO_H_
#define  _PROJECT_IO_H_

/************************************************************************/
/*					Interaction with file system                        */
/*@Author: Wenchao Hu							                        */
/************************************************************************/

#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>
#include <functional>
#include <cctype>
#include "spm_cgal.h"
#include "meshes.h"

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

	/** input **/
	// load a whole project from a file, the same function as the constructor
	void load_project(const string& prj_configure_file);
	// read 2D polygons
	ProjectIO& operator>>(vector<Polygon_2>& polygons);

	// read triangle mesh
	ProjectIO& operator>>(TriMesh& mesh);
	
	
	/** output **/
	void dump_results(const string& filename, const string& directory="");

	/** debug **/
	void debug_print();

private:
	
	string get_start_token(std::ifstream& fs);
	string get_end_token(std::ifstream& fs);
	
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

	/** indicators of what are input **/
	bool has_density_input_; // each vertex has a weight, e.g. curvature

	/** i/o error message **/
	const static string error_start_tag;
	const static string error_end_tag;
	const static string error_tag_dismatch;
	const static string error_prj_file_fail;
	const static string error_pgn_file_fail;
	const static string error_mesh_file_fail;
};

}
#endif