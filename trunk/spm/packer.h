#ifndef _PACKER_H_
#define _PACKER_H_

#include <vector>
#include <string>
#include <cfloat>
#include <cmath>
#include <set>
#include "project_io.h"
#include "packing_object.h"
#include "meshes.h"

namespace Geex
{
	class Packer
	{
	public:
		
		Packer();

		void get_bbox(real& x_min, real& y_min, real& z_min, real& x_max, real& y_max, real& z_max);

		/** initialization **/
		void load_project(const std::string& prj_config_file);

		void initialize(); // top-level initialization function. distribute polygons in some way
		

		/** access functions **/
		const TriMesh& mesh_domain() const { return mesh; }
		const vector<Packing_object>& get_tiles() const { return pack_objects; }

	private:

		void random_init_tiles(unsigned int nb_init_polygons);	// put polygons according to curvature (if any) and in a uniform way

	private:

		static Packer *instance_;
		
		ProjectIO pio;
		std::vector<Packing_object> pack_objects;
		std::vector<Polygon_2> pgn_lib;
		TriMesh mesh;
	};
}

#endif