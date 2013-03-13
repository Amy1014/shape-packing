#include "packer.h"

namespace Geex
{
	Packer* Packer::instance_ = nil;

	void Packer::get_bbox(real& x_min, real& y_min, real& z_min, real& x_max, real& y_max, real& z_max)
	{
		x_min = y_min = z_min = DBL_MAX;
		x_max = y_max = z_max = DBL_MIN;
		for(unsigned int i=0; i<mesh.size(); i++) 
			for(unsigned int j=0; j<3; j++) 
			{
				const vec3& p = mesh[i].vertex[j] ;
				x_min = gx_min(x_min, p.x); y_min = gx_min(y_min, p.y);	z_min = gx_min(z_min, p.z);
				x_max = gx_max(x_max, p.x);	y_max = gx_max(y_max, p.y); z_max = gx_max(z_max, p.z) ;
			}
	}

	Packer::Packer()
	{
		assert( instance_ == nil );
		instance_ = this;
		srand(1);
	}

	void Packer::load_project(const std::string& prj_config_file)
	{
		pio.load_project(prj_config_file);
		pio>>mesh>>pgn_lib;
		initialize();
		std::cout<<"Start computing RDT...\n";
		generate_RDT();
	}

	void Packer::initialize()
	{
		mesh.build_kdtree();
		// compute mesh area
		double mesh_area = 0.0;
		for (unsigned int i = 0; i < mesh.size(); i++)
			mesh_area += std::fabs(mesh[i].area());
		std::cout<<"Surface area: "<<mesh_area<<std::endl;

		// compute maximum and average polygon area in the library
		double max_pgn_area = DBL_MIN, mean_pgn_area = 0.0;
		for (unsigned int i = 0; i < pgn_lib.size(); i++)
		{
			double s = std::fabs(pgn_lib[i].area());
			max_pgn_area = std::max(max_pgn_area, s);
			mean_pgn_area += s;
		}
		mean_pgn_area /= pgn_lib.size();
		std::cout<<"Maximum polygon area: "<<max_pgn_area<<std::endl;

		// distribute polygons by different strategies
		// 1. random 
		unsigned int nb_init_polygons = mesh_area / mean_pgn_area;
		random_init_tiles(nb_init_polygons);
		
	}
	void Packer::random_init_tiles(unsigned int nb_init_polygons)
	{
		set<int> init_facets; // on which facets the location are
		if (pio.has_density_input()) // put according to density function
		{
			double total_weight = 0.0, res = 0.0;
			for (unsigned int i = 0; i < mesh.nb_vertices(); i++)
				total_weight += mesh.vertex(i).weight();
			for (unsigned int i = 0; i < mesh.nb_vertices(); i++)
			{
				double nf = nb_init_polygons*mesh.vertex(i).weight()/total_weight;
				int n(nf+res);
				if ( n >= 1 )
				{
					res = nf + res - n; // fractional part
					int nbput = std::min<int>(n, mesh.vertex(i).faces_.size());
					for (int j = 0; j < nbput; j++)
					{
						int idx = mesh.vertex(i).faces_[j];
						if ( init_facets.find(idx) != init_facets.end() )	continue;
						init_facets.insert(idx);
					}
				}
				else 
					res += nf;
			}
			// pick up the lost ones
			while ( init_facets.size() < nb_init_polygons )
			{
				int idx = ::rand()%mesh.size();
				while ( init_facets.find(idx) != init_facets.end() )
					idx = ::rand()%mesh.size();
				init_facets.insert(idx);
			}
		}
		else   // if no density specified, uniformly put	
		{
			for (unsigned int i = 0; i < nb_init_polygons; i++)
			{
				int idx = ::rand()%mesh.size();
				while ( init_facets.find(idx) != init_facets.end() )
					idx = ::rand()%mesh.size();
				init_facets.insert(idx);
			}
		}
		// choose polygons and distribute
		pack_objects.reserve(nb_init_polygons);
		unsigned int pgn_lib_idx = 0;
		for (set<int>::const_iterator it = init_facets.begin(); it != init_facets.end(); ++it)
		{
			const Facet& f = mesh[*it];
			vec3 gx_cent = (f.vertex[0] + f.vertex[1] + f.vertex[2])/3.0;
			//vec3 gx_normal = mesh[*it].normal();
			vec3 gx_normal = approx_normal(*it);
			const Polygon_2& pgn_2 = pgn_lib[pgn_lib_idx%pgn_lib.size()];
			// shrink factor
			double fa = std::fabs(f.area()), pa = std::fabs(pgn_2.area());
			double s = std::min(fa/pa, pa/fa);
			s = 0.1*std::sqrt(s);
			Polygon_2 init_polygon = CGAL::transform(Transformation_2(CGAL::SCALING, s), pgn_2);
			pack_objects.push_back(Packing_object(init_polygon, to_cgal_vec(gx_normal), to_cgal_pnt(gx_cent), s));
			pgn_lib_idx++;
		}
	}

	void Packer::generate_RDT()
	{
		rpvd.set_mesh(pio.attribute_value("MeshFile"));
		rpvd.set_trimesh(&mesh);
		rpvd.begin_insert();
		rpvd.insert_polygons(pack_objects.begin(), pack_objects.end(), 12);
		rpvd.insert_bounding_points(10);
		rpvd.end_insert();
		std::vector<Plane_3> tangent_planes;
		tangent_planes.reserve(pack_objects.size());
		for (unsigned int i = 0; i < pack_objects.size(); i++)
		{
			Point_3 c = pack_objects[i].centroid();
			Vector_3 n = pack_objects[i].norm();
			tangent_planes.push_back(Plane_3(c, n));
		}
		rpvd.compute_clipped_VD(tangent_planes);
	}

	vec3 Packer::approx_normal(unsigned int facet_idx)
	{
		const Facet& f = mesh[facet_idx];
		int vi0 = f.vertex_index[0], vi1 = f.vertex_index[1], vi2 = f.vertex_index[2];
		const MeshVertex& v0 = mesh.vertex(vi0);
		const MeshVertex& v1 = mesh.vertex(vi1);
		const MeshVertex& v2 = mesh.vertex(vi2);
		const vector<int>& neighbor0 = v0.faces_;
		const vector<int>& neighbor1 = v1.faces_;
		const vector<int>& neighbor2 = v2.faces_;
		vec3 avg_norm(0.0, 0.0, 0.0);
		if (!mesh.is_on_boundary(vi0) && !mesh.is_on_feature(vi0))
			for (unsigned int i = 0; i < neighbor0.size(); i++)
				avg_norm += mesh[neighbor0[i]].normal();
		if (!mesh.is_on_boundary(vi1) && !mesh.is_on_feature(vi1))
			for (unsigned int i = 0; i < neighbor1.size(); i++)
				avg_norm += mesh[neighbor1[i]].normal();
		if (!mesh.is_on_boundary(vi2) && !mesh.is_on_feature(vi2))
			for (unsigned int i = 0; i < neighbor2.size(); i++)
				avg_norm += mesh[neighbor2[i]].normal();
		if (avg_norm.x == 0.0 && avg_norm.y == 0.0 && avg_norm.z == 0.0)
			return f.normal();
		avg_norm /= avg_norm.length();
		return avg_norm;
	}
}