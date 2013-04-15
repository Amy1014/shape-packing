#include "packer.h"

namespace Geex
{
	Packer* Packer::instance_ = nil;

	double Packer::PI = 3.141592653589793;

	Packer::Packer()
	{
		assert( instance_ == nil );
		instance_ = this;
		//srand(1);
		srand(time(NULL));
		rot_lower_bd = -PI/6.0;
		rot_upper_bd = PI/6.0;

		frontier_edge_size = 1.0;
		hole_face_size = 1.0;

		stop_update_DT = false;

		epsilon = 0.15;
		match_weight = 0.0;

		samp_nb = 20;

		submesh_id = 0;
	}

	void Packer::load_project(const std::string& prj_config_file)
	{
		pio.load_project(prj_config_file);
		if (!pio.multi_meshes_specified())
		{
			pio>>mesh>>pgn_lib;
			initialize();
			std::cout<<"Start computing RDT...\n";
			generate_RDT();
			compute_clipped_VD();
		}
		else
		{
			pio>>pgn_lib;
			pio>>mesh_segments;
		}
		//rpvd.save_triangulation("rdt.obj");
	}

	void Packer::pack_next_submesh()
	{
		//static unsigned int mesh_id = 0;
		mesh.clear();
		pack_objects.clear();
		stop_update_DT = false;
		mesh = mesh_segments[submesh_id];
		initialize();
		std::cout<<"start "<<submesh_id<<" mesh\n";
		generate_RDT();
		compute_clipped_VD();
		//mesh_id++;
	}

	void Packer::write_to_results()
	{
		//static unsigned int mesh_id = 0;
		res_pack_objects.push_back(std::vector<Packing_object>());
		res_pack_objects.back().reserve(pack_objects.size());
		for (unsigned int i = 0; i < pack_objects.size(); i++)
			res_pack_objects.back().push_back(pack_objects[i]);
		submesh_id++;
		
	}

	void Packer::initialize()
	{
		mesh.build_kdtree();
		
		// compute mesh area
		mesh_area = 0.0;
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
		if (!pio.multi_meshes_specified())
			rpvd.set_mesh(pio.attribute_value("MeshFile"));
		else
			rpvd.set_mesh(pio.get_submesh_files().at(submesh_id));
		rpvd.set_trimesh(&mesh);

		// distribute polygons by different strategies
		// 1. random 
		unsigned int nb_init_polygons = mesh_area / mean_pgn_area;
		random_init_tiles(nb_init_polygons);
		std::cout<<pack_objects.size()<<" polygons are loaded initially.\n";
		// 2. multi-class or more (TODO)
	}
	void Packer::random_init_tiles(unsigned int nb_init_polygons)
	{
		set<int> init_facets; // on which facets the location are
		std::vector<std::pair<int, Point_3>> init_pos;
		init_pos.reserve(nb_init_polygons);
		if (pio.has_density_input()) // put according to density function
		{
			double total_weight = 0.0, res = 0.0;
			std::vector<double> region_area(mesh.nb_vertices());
			for (unsigned int i = 0; i < mesh.nb_vertices(); i++)
			{
				region_area[i] = mesh.surrounding_region_area(i);
				total_weight += mesh.vertex(i).weight()*region_area[i];
			}
			for (unsigned int i = 0; i < mesh.nb_vertices(); i++)
			{
				if (mesh.is_on_boundary(i) || mesh.is_on_feature(i))
					continue;
				double nf = nb_init_polygons*mesh.vertex(i).weight()*region_area[i]/total_weight;
				int n(nf+res);
				//int n(nf);
				//if ( n < 1)
				//	n = int(nf+0.5);
				//int n(nf+0.5);
				if ( n >= 1 )
				{
					res = nf + res - n; // fractional part
					int nbput = std::min<int>(n, mesh.vertex(i).faces_.size());
					if (nbput < n)
						std::cout<<"inadequate sampling!\n";
					for (int j = 0; j < nbput; j++)
					{
						int idx = mesh.vertex(i).faces_[j];
						Point_3 c = CGAL::centroid(to_cgal_pnt(mesh[idx].vertex[0]), to_cgal_pnt(mesh[idx].vertex[1]), to_cgal_pnt(mesh[idx].vertex[2]));
						init_pos.push_back(std::make_pair(idx, CGAL::midpoint(to_cgal_pnt(mesh.vertex(i).pos_), c)));
						//if ( init_facets.find(idx) != init_facets.end() )	continue;
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
				int vi0 = mesh[idx].vertex_index[0], vi1 = mesh[idx].vertex_index[1], vi2 = mesh[idx].vertex_index[2];
				if (mesh.is_on_boundary(vi0) || mesh.is_on_boundary(vi1) || mesh.is_on_boundary(vi2))
					continue;
				Point_3 c = CGAL::centroid(to_cgal_pnt(mesh[idx].vertex[0]), to_cgal_pnt(mesh[idx].vertex[1]), to_cgal_pnt(mesh[idx].vertex[2]));
				init_pos.push_back(std::make_pair(idx, c));
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
				Point_3 c = CGAL::centroid(to_cgal_pnt(mesh[idx].vertex[0]), to_cgal_pnt(mesh[idx].vertex[1]), to_cgal_pnt(mesh[idx].vertex[2]));
				init_pos.push_back(std::make_pair(idx, c));
			}
		}
		// choose polygons and distribute
		pack_objects.reserve(nb_init_polygons);
		unsigned int pgn_lib_idx = 0;
		for (unsigned int i = 0; i < init_pos.size(); i++)
		{
			const Facet& f = mesh[init_pos[i].first];
			vec3 gx_normal = approx_normal(init_pos[i].first);
			const Ex_polygon_2& pgn_2 = pgn_lib[pgn_lib_idx%pgn_lib.size()];
			// shrink factor
			double fa = std::fabs(f.area()), pa = std::fabs(pgn_2.area());
			double s = std::min(fa/pa, pa/fa);
			s = 0.05*std::sqrt(s);
			Polygon_2 init_polygon = CGAL::transform(Transformation_2(CGAL::SCALING, s), pgn_2);
			pack_objects.push_back(Packing_object(init_polygon, to_cgal_vec(gx_normal), init_pos[i].second, s));
			pack_objects.back().lib_idx = pgn_lib_idx%pgn_lib.size();
			pack_objects.back().facet_idx = init_pos[i].first;
			pack_objects.back().texture_id = pgn_2.texture_id;
			pack_objects.back().texture_coord.assign(pgn_2.texture_coords.begin(), pgn_2.texture_coords.end());
			pgn_lib_idx++;
		}
	}

#if 0
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
				double nf = (2-epsilon)*epsilon*nb_init_polygons*mesh.vertex(i).weight()/total_weight;
				//int n(nf+res);
				int n(nf);
				if ( n < 1)
					n = int(nf+0.5);
				//else
				//	n = int(nf);
				//int n(nf+0.5);
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
			const Ex_polygon_2& pgn_2 = pgn_lib[pgn_lib_idx%pgn_lib.size()];
			// shrink factor
			double fa = std::fabs(f.area()), pa = std::fabs(pgn_2.area());
			double s = std::min(fa/pa, pa/fa);
			s = 0.2*std::sqrt(s);
			Polygon_2 init_polygon = CGAL::transform(Transformation_2(CGAL::SCALING, s), pgn_2);
			pack_objects.push_back(Packing_object(init_polygon, to_cgal_vec(gx_normal), to_cgal_pnt(gx_cent), s));
			pack_objects.back().lib_idx = pgn_lib_idx%pgn_lib.size();
			pack_objects.back().facet_idx = *it;
			pack_objects.back().texture_id = pgn_2.texture_id;
			pack_objects.back().texture_coord.assign(pgn_2.texture_coords.begin(), pgn_2.texture_coords.end());
			pgn_lib_idx++;
		}
	}
#endif

	void Packer::compute_clipped_VD()
	{
		std::cout<<"Start computing clipped Voronoi region\n";
		std::vector<Plane_3> tangent_planes;
		tangent_planes.reserve(pack_objects.size());
		for (unsigned int i = 0; i < pack_objects.size(); i++)
		{
			Point_3 c = pack_objects[i].centroid();
			Vector_3 n = pack_objects[i].norm();
			tangent_planes.push_back(Plane_3(c, n));
		}
		rpvd.compute_clipped_VD(tangent_planes);
		std::cout<<"End computing.\n";
	}

	void Packer::generate_RDT()
	{
		rpvd.begin_insert();
		rpvd.insert_polygons(pack_objects.begin(), pack_objects.end(), samp_nb);
		rpvd.insert_bounding_points(40);
		CGAL::Timer t;
		t.start();
		rpvd.end_insert();
		t.stop();
		std::cout<<"time spent in RDT: "<<t.time()<<" seconds."<<std::endl;
		//compute_clipped_VD();
	}

	vec3 Packer::approx_normal(unsigned int facet_idx)
	{
		const Facet& f = mesh[facet_idx];
		vec3 original_normal = f.normal();
		//return original_normal;

		int vi0 = f.vertex_index[0], vi1 = f.vertex_index[1], vi2 = f.vertex_index[2];
		const MeshVertex& v0 = mesh.vertex(vi0);
		const MeshVertex& v1 = mesh.vertex(vi1);
		const MeshVertex& v2 = mesh.vertex(vi2);
		const vector<int>& neighbor0 = v0.faces_;
		const vector<int>& neighbor1 = v1.faces_;
		const vector<int>& neighbor2 = v2.faces_;
		vec3 avg_norm(0.0, 0.0, 0.0);
		const std::set<int>& high_cur_facets = mesh.getHighCurvatureFacets();
		if (!mesh.is_on_boundary(vi0) && !mesh.is_on_feature(vi0))
			for (unsigned int i = 0; i < neighbor0.size(); i++)
			{
				//if (high_cur_facets.find(neighbor0[i]) == high_cur_facets.end())
				vec3 n = mesh[neighbor0[i]].normal();
				if ( dot(n, original_normal) >= 0.5 )
					avg_norm += n;
			}
		if (!mesh.is_on_boundary(vi1) && !mesh.is_on_feature(vi1))
			for (unsigned int i = 0; i < neighbor1.size(); i++)
			{
				//if (high_cur_facets.find(neighbor1[i]) == high_cur_facets.end())
				vec3 n = mesh[neighbor1[i]].normal();
				if ( dot(n, original_normal) >= 0.5 )
					avg_norm += n;
			}
		if (!mesh.is_on_boundary(vi2) && !mesh.is_on_feature(vi2))
			for (unsigned int i = 0; i < neighbor2.size(); i++)
			{
				//if (high_cur_facets.find(neighbor2[i]) == high_cur_facets.end())
				vec3 n = mesh[neighbor2[i]].normal();
				if (dot(n, original_normal) >= 0.5)
					avg_norm += n;
			}
		if (avg_norm.x == 0.0 && avg_norm.y == 0.0 && avg_norm.z == 0.0)
			return f.normal();
		avg_norm /= (neighbor0.size()+neighbor1.size()+neighbor2.size());
		avg_norm /= avg_norm.length();
		return avg_norm;
	}

	Packer::Local_frame Packer::compute_local_frame(const Packing_object& tile)
	{
		Local_frame lf;
		lf.o = tile.centroid();
		lf.w = tile.norm(); // assume this vector has already been normalized
		Vector_3 u(lf.o, tile.vertex(0));
		lf.u = u/CGAL::sqrt(u.squared_length());
		lf.v = CGAL::cross_product(lf.w, lf.u);
		return lf;
	}
	Packer::Optimization_res Packer::optimize_one_polygon(unsigned int id, Local_frame& lf, Parameter& solution)
	{
		// choose a local frame
		lf = compute_local_frame(pack_objects[id]);

		// voronoi region
		std::vector<Segment_2> vr_region;
#ifdef _CILK_
		Containment& ctm = containments[id];
#else
		Containment& ctm = containment;
#endif
		// formulate optimization 
		ctm.clear();
		const RestrictedPolygonVoronoiDiagram::VertGroup& samp_pnts = rpvd.sample_points_group(id);
		if (!stop_update_DT)
			for (unsigned int i = 0; i < samp_pnts.size(); i++)
			{
				const vector<Point_3>& bisec_pnts = samp_pnts[i]->vd_vertices;
				Point_2 ref_point = lf.to_uv(samp_pnts[i]->point());
				if (bisec_pnts.size() == 0)
					continue;
				assert(bisec_pnts.size() != 1);
				for (unsigned int j = 0; j < bisec_pnts.size()-1; j++)
				{
					const Point_3& s = bisec_pnts[j];
					const Point_3& t = bisec_pnts[j+1];
					Point_2 s2 = lf.to_uv(s), t2 = lf.to_uv(t);
					ctm.const_list_.push_back(Constraint(lf.to_uv(samp_pnts[i]->point()), Segment_2(s2, t2), ref_point));
					//vr_region.push_back(Segment_2(s2, t2));
				}
				//const std::vector<double>& ws = samp_pnts[i]->weights;
				//for (unsigned int j = 0; j < bisec_pnts.size(); j++)
				//{
				//	containments[id].vd_vertices.push_back(lf.to_uv(bisec_pnts[j]));
				//	containments[id].weights.push_back(ws[j]);
				//}
			}
		else
			for (unsigned int i = 0; i < pack_objects[id].size(); i++)
			{
				for (unsigned int j = 0; j < samp_pnts.size(); j++)
				{
					const vector<Point_3>& bisec_pnts = samp_pnts[j]->vd_vertices;
					if (bisec_pnts.size() == 0)
						continue;
					assert(bisec_pnts.size() != 1);
					for (unsigned int k = 0; k < bisec_pnts.size()-1; k++)
					{
						const Point_3& s = bisec_pnts[k];
						const Point_3& t = bisec_pnts[k+1];
						Point_2 s2 = lf.to_uv(s), t2 = lf.to_uv(t);
						ctm.const_list_.push_back(Constraint(lf.to_uv(pack_objects[id].vertex(i)), Segment_2(s2, t2)));
						//vr_region.push_back(Segment_2(s2, t2));
					}
					//const std::vector<double>& ws = samp_pnts[j]->weights;
					//for (unsigned int k = 0; k < bisec_pnts.size(); k++)
					//{
					//	containments[id].vd_vertices.push_back(lf.to_uv(bisec_pnts[k]));
					//	containments[id].weights.push_back(ws[k]);
					//}
				}
			}
		// construct extended objective function
		//containments[id].region_area = 0.0;
		//for (unsigned int i = 0; i < vr_region.size(); i++)
		//	containments[id].region_area += std::fabs(CGAL::area(vr_region[0].source(), vr_region[i].source(), vr_region[i].target()));
		//containments[id].pgn_cent = lf.to_uv(lf.o);

		const unsigned int try_time_limit = 10;
		unsigned int try_times = 1;
		int knitro_res;
		while ( (knitro_res = KTR_optimize(&solution.k, &solution.theta, &solution.tx, &solution.ty, id)) && try_times < try_time_limit )
		{
			double r = Numeric::random_float64();
			solution = Parameter(1.0, (rot_upper_bd-rot_lower_bd)*r+rot_lower_bd, 0.0, 0.0);
			try_times++;
		}
		
		if (try_times == try_time_limit)
			return FAILED;
		else
		{
			//std::cout<<"k = "<<solution.k<<", theta = "<<solution.theta<<", tx = "<<solution.tx<<", ty = "<<solution.ty<<std::endl;
			return SUCCESS;
		}
	}

	Packer::Lloyd_res Packer::one_lloyd(bool enlarge, std::vector<Parameter>& solutions, std::vector<Local_frame>& lfs)
	{
		const double min_scalor = 1.05;
		std::vector<Optimization_res> opti_res(pack_objects.size());
#ifdef _CILK_
		containments.resize(pack_objects.size());
		cilk_for (unsigned int i = 0; i < pack_objects.size(); i++)
#else
		for (unsigned int i = 0; i < pack_objects.size(); i++)
#endif
		{
			opti_res[i] = optimize_one_polygon(i, lfs[i], solutions[i]);
			if (!enlarge)
				solutions[i].k = 1.0;
		}
		double mink = std::numeric_limits<double>::max();
		for (unsigned int i = 0; i < pack_objects.size(); i++)
		{
			if (opti_res[i] == SUCCESS)
				mink = std::min(solutions[i].k, mink);
			else
			{
				solutions[i].k = 1.0;
				solutions[i].theta = solutions[i].tx = solutions[i].ty = 0.0;
			}
		}
		if (mink >= min_scalor)
		{
			for (unsigned int i = 0; i < pack_objects.size(); i++)
				if (opti_res[i] == SUCCESS)
					solutions[i].k = mink;
		}
		std::cout<<"Minimum enlargement factor = "<<mink<<std::endl;
		std::vector<double> min_factors(solutions.size(), 1.0);// factors to constraint transformation
		if (!stop_update_DT)
		{
			constraint_transformation(solutions, lfs, min_factors);
		}
		//std::ofstream of("parameters.txt");
#ifdef _CILK_
		cilk_for (unsigned int i = 0; i < pack_objects.size(); i++)
#else
		for (unsigned int i = 0; i < pack_objects.size(); i++)
#endif
		{
			Parameter constrained_parameter = solutions[i]*min_factors[i];
			curv_constrained_transform(constrained_parameter, pack_objects[i].facet_idx, i);
			transform_one_polygon(i, lfs[i], constrained_parameter);
			solutions[i] -= constrained_parameter;// residual transformation
		}
		if (enlarge && mink < min_scalor)
			return NO_MORE_ENLARGEMENT;
		else
			return MORE_ENLARGEMENT;
	}

	void Packer::transform_one_polygon(unsigned int id, Local_frame& lf, Parameter& param)
	{
		Apply_transformation t(param, lf);
		std::transform(pack_objects[id].vertices_begin(), pack_objects[id].vertices_end(), pack_objects[id].vertices_begin(), t);
		pack_objects[id].factor *= param.k;
		RestrictedPolygonVoronoiDiagram::VertGroup& samp_pnts = rpvd.sample_points_group(id);
		std::transform(samp_pnts.begin(), samp_pnts.end(), samp_pnts.begin(), t);
		Point_3 c = pack_objects[id].centroid();
		vec3 dummyv, p;
		int fid;
		p = mesh.project_to_mesh(to_geex_pnt(c), dummyv, fid);
		vec3 n = approx_normal(fid);
		Transformation_3 same_t = pack_objects[id].align(to_cgal_vec(n), to_cgal_pnt(p));
		for (unsigned int j = 0; j < samp_pnts.size(); j++)
			samp_pnts[j]->point() = same_t(samp_pnts[j]->point());
		for (unsigned int j = 0; j < samp_pnts.size(); j++)
		{
			vec3 dummyv;
			Point_3 temp = samp_pnts[j]->point();
			samp_pnts[j]->mp = to_cgal_pnt(mesh.project_to_mesh(to_geex_pnt(temp), dummyv));
		}
		pack_objects[id].facet_idx = fid;
	}
	void Packer::lloyd(void (*post_action)(), bool enlarge)
	{
		static int times = 0;

		std::vector<Parameter> solutions(pack_objects.size());
		std::vector<Local_frame> local_frames(pack_objects.size());
		for (unsigned int i = 0; i < 1; i++)
		{
			std::cout<<"============ lloyd iteration "<<times++<<" ============\n";

			//rpvd.save_triangulation("pre_lloyd.obj");
			Lloyd_res res = one_lloyd(enlarge, solutions, local_frames);

			if (res == NO_MORE_ENLARGEMENT)
			{
				stop_update_DT = true;
				std::cout<<"Caution: Stopped updating iDT\n";
			}
			if (!stop_update_DT)
			{
				std::cout<<"Start iDT updating...\n";

				//std::ostringstream pre_fn(std::ostringstream::ate);
				//pre_fn.str("../imme_data/");
				//pre_fn<<"pre_rdt_"<<times<<".obj";
				//rpvd.save_triangulation(pre_fn.str());

				CGAL::Timer t;
				t.start();
				//rpvd.save_triangulation("pre_idt.obj");
				rpvd.iDT_update();
				//rpvd.save_triangulation("idt.obj");
				//generate_RDT();
				t.stop();
				std::cout<<"Time spent in iDT: "<<t.time()<<" seconds.\n";
				std::cout<<"End iDT updating\n";	

				if (post_action != NULL)
					post_action();

				// compensate for the constrained transformation in the last step
				//std::cout<<"Compensate for lost transformation...\n";
				//std::vector<double> min_factors(solutions.size());
				//for (unsigned int comp_times = 0; comp_times < 0; comp_times++)
				//{
				//	constraint_transformation(solutions, local_frames, min_factors);
				//	for (unsigned int j = 0; j < pack_objects.size(); j++)
				//	{
				//		Parameter constrained_param = solutions[j]*min_factors[j];
				//		transform_one_polygon(j, local_frames[j], constrained_param);
				//		solutions[j] -= constrained_param;
				//	}
				//	//rpvd.save_triangulation("pre_idt.obj");
				//	rpvd.iDT_update();
				//	//rpvd.save_triangulation("idt.obj");
				//	glut_viewer_redraw();
				//}
				compute_clipped_VD();
				//std::cout<<"End compensating.\n";

				//std::ostringstream fn(std::ostringstream::ate);
				//fn.str("../imme_data/");
				//fn<<"rdt_"<<times<<".obj";
				//rpvd.save_triangulation(fn.str());
			}

			if (post_action != NULL)
				post_action();
		}
	}
	
	void Packer::pack(void (*post_action)())
	{
		
		lloyd(post_action, true);

		// compute area coverage
#ifdef _CILK_
		cilk::reducer_opadd<double> sum_pgn_area(0.0);
		cilk_for (unsigned int i = 0; i < pack_objects.size(); i++)
			sum_pgn_area += pack_objects[i].area();
		std::cout<<"-- Area coverage ratio: "<<sum_pgn_area.get_value()/mesh_area<<std::endl;
#else
		double sum_pgn_area = 0.0;
		for (unsigned int i = 0; i < pack_objects.size(); i++)
			sum_pgn_area += pack_objects[i].area();
		std::cout<<"-- Area coverage ratio: "<<sum_pgn_area/mesh_area<<std::endl;
#endif

		
		if (post_action != NULL)
			post_action();
	}

	void Packer::constraint_transformation(vector<Parameter>& parameters, vector<Local_frame>& lfs, vector<double>& min_factors)
	{
		min_factors.assign(parameters.size(), 1.0);
		typedef RestrictedPolygonVoronoiDiagram RPVD;
		for (RPVD::Facet_iterator fit = rpvd.faces_begin(); fit != rpvd.faces_end(); ++fit)
		{
			RPVD::Halfedge_handle eh = fit->halfedge();
			Point_3 v0 = eh->vertex()->mp;
			int i0 = eh->vertex()->group_id;

			eh = eh->next();
			Point_3 v1 = eh->vertex()->mp;
			int i1 = eh->vertex()->group_id;

			eh = eh->next();
			Point_3 v2 = eh->vertex()->mp;
			int i2 = eh->vertex()->group_id;

// 			if (i0 < 0 || i1 < 0 || i2 < 0)
// 				continue;
			fit->n = CGAL::cross_product(Vector_3(v0, v1), Vector_3(v1, v2));
			fit->n = fit->n / CGAL::sqrt(fit->n.squared_length());
		}
		unsigned int N = 20, n = 1;
		bool stop = false;
		while (n < N && !stop)
		{
#ifdef _CILK_
			cilk_for (unsigned int i = 0; i < pack_objects.size(); i++)
#else
			for (unsigned int i = 0; i < pack_objects.size(); i++)
#endif
			{
				Parameter sp = parameters[i]*min_factors[i];
				Apply_transformation t(sp, lfs[i]);
				Packing_object future_pgn = pack_objects[i];
				std::transform(future_pgn.vertices_begin(), future_pgn.vertices_end(), future_pgn.vertices_begin(), t);
				RestrictedPolygonVoronoiDiagram::VertGroup& samp_pnts = rpvd.sample_points_group(i);
				//std::transform(samp_pnts.begin(), samp_pnts.end(), samp_pnts.begin(), t);
				std::vector<Point_3> future_samp_pnts;
				future_samp_pnts.reserve(samp_pnts.size());
				for (unsigned int j = 0; j < samp_pnts.size(); j++)
					future_samp_pnts.push_back(samp_pnts[j]->point());
				std::transform(future_samp_pnts.begin(), future_samp_pnts.end(), future_samp_pnts.begin(), t);
				Point_3 c = future_pgn.centroid();
				vec3 dummyv, p;
				int fid;
				p = mesh.project_to_mesh(to_geex_pnt(c), dummyv, fid);
				vec3 n = approx_normal(fid);
				Transformation_3 same_t = future_pgn.align(to_cgal_vec(n), to_cgal_pnt(p));
				for (unsigned int j = 0; j < future_samp_pnts.size(); j++)
				{
					future_samp_pnts[j] = same_t(future_samp_pnts[j]);
				}
				for (unsigned int j = 0; j < future_samp_pnts.size(); j++)
				{
					vec3 dummyv;
					samp_pnts[j]->mp = to_cgal_pnt(mesh.project_to_mesh(to_geex_pnt(future_samp_pnts[j]), dummyv));
				}
			}
			stop = true;
			std::vector<bool> flipped(parameters.size(), false);
			int nb_flipped = 0;
			for (RPVD::Facet_iterator fit = rpvd.faces_begin(); fit != rpvd.faces_end(); ++fit)
			{
				RPVD::Halfedge_handle eh = fit->halfedge();
				Point_3 v0 = eh->vertex()->mp;
				int i0 = eh->vertex()->group_id;

				eh = eh->next();
				Point_3 v1 = eh->vertex()->mp;
				int i1 = eh->vertex()->group_id;

				eh = eh->next();
				Point_3 v2 = eh->vertex()->mp;
				int i2 = eh->vertex()->group_id;
// 
// 				if (i0 < 0 || i1 < 0 || i2 < 0)
// 					continue;
				Vector_3 n = CGAL::cross_product(Vector_3(v0, v1), Vector_3(v1, v2));
				double nlen2 = n.squared_length();
				if (nlen2 == 0.0) //degenerate case
				{
					//std::cout<<"degenerate at <"<<i0<<", "<<i1<<", "<<i2<<">\n";
					continue;
				}
				n = n / CGAL::sqrt(nlen2);
				if (n*fit->n < 0.0)
				{
					nb_flipped++;
					if (i0 >= 0)	flipped[i0] = true;
					if (i1 >= 0)	flipped[i1] = true;
					if (i2 >= 0)	flipped[i2] = true;
					stop = false;
				}
			}
			//std::cout<<"number of flipped facets: "<<nb_flipped<<std::endl;
			n++;
			for (unsigned int j = 0; j < parameters.size() ; j++)
			{
				if (flipped[j])
				{
					if (/*n>= N-2 && n <= N*/n>=N-1) // last time
						min_factors[j] *= 0.0;
					else
						min_factors[j] *= 0.5;
				}
			}
		}
		std::cout<<"Try times for constraint: "<<n<<std::endl;
			
	}

	void Packer::detect_holes()
	{
		std::for_each(holes.begin(), holes.end(), std::mem_fun_ref(&Hole::clear));
		holes.clear();

		double avg_area = 0.0;
		for (unsigned int i = 0; i < pack_objects.size(); i++)
			avg_area += pack_objects[i].area();
		avg_area /= pack_objects.size();
		double front_edge_len2_threshold = frontier_edge_size*avg_area;
		double hole_facet_area_threshold = hole_face_size*avg_area;
		std::queue<Facet_iterator> seed_faces;
		// detect faces in holes
		for (Facet_iterator fit = rpvd.faces_begin(); fit != rpvd.faces_end(); ++fit)
		{
			Halfedge_handle eh = fit->halfedge();
			Vertex_handle v0 = eh->vertex();
			eh = eh->next();
			Vertex_handle v1 = eh->vertex();
			eh = eh->next();
			Vertex_handle v2 = eh->vertex();
			if (v0->group_id == v1->group_id && v0->group_id == v2->group_id && v1->group_id == v2->group_id)
				fit->vacant = false;
			else
				fit->vacant = true;
			if (fit->vacant)
			{
				double	d01 = CGAL::squared_distance(v0->mp, v1->mp),
						d02 = CGAL::squared_distance(v0->mp, v2->mp),
						d12 = CGAL::squared_distance(v1->mp, v2->mp);
				if (d01 >= front_edge_len2_threshold || d02 >= front_edge_len2_threshold || d12 >= front_edge_len2_threshold)
					seed_faces.push(fit);
			}
		}

		while (!seed_faces.empty())
		{
			Facet_iterator f = seed_faces.front();
			seed_faces.pop();
			if (!f->vacant)
				continue;
			f->vacant = false;
			std::queue<Halfedge_handle> active_front_edges;
			std::list<Facet_iterator> ahole;
			std::vector<Halfedge_handle> hole_bd;
			Halfedge_handle e = f->halfedge();
			active_front_edges.push(e);
			active_front_edges.push(e->next());
			active_front_edges.push(e->next()->next());
			ahole.push_back(f);
			while (!active_front_edges.empty())
			{
				Halfedge_handle curr_edge = active_front_edges.front();
				active_front_edges.pop();
				Halfedge_handle occur_edge = curr_edge->opposite();
				double len = CGAL::squared_distance(curr_edge->vertex()->mp, occur_edge->vertex()->mp);
				if (!occur_edge->is_border() && len >= front_edge_len2_threshold)
				{		
					Facet_handle neighbor_face = occur_edge->face();
					if (neighbor_face->vacant)
					{
						active_front_edges.push(occur_edge->prev());
						active_front_edges.push(occur_edge->next());
						neighbor_face->vacant = false;
						ahole.push_back(neighbor_face);
					}
					else
						hole_bd.push_back(curr_edge);
				}
				else
					hole_bd.push_back(curr_edge);
			}
			double ahole_area = 0.0;
			for (std::list<Facet_iterator>::const_iterator it = ahole.begin(); it != ahole.end(); ++it)
			{
				Halfedge_handle eh = (*it)->halfedge();
				Triangle_3 tri(eh->vertex()->mp, eh->next()->vertex()->mp, eh->next()->next()->vertex()->mp);
				ahole_area += CGAL::sqrt(tri.squared_area());
			}

			if (ahole_area >= hole_facet_area_threshold)
			{
				// convert to point representation
				holes.push_back(Hole());
				holes.back().reserve(hole_bd.size());
				for (unsigned int i = 0; i < hole_bd.size(); i++)
				{
					Vertex_handle src = hole_bd[i]->vertex();
					Vertex_handle tgt = hole_bd[i]->opposite()->vertex();
					holes.back().push_back(std::make_pair(src, tgt));
				}
				//holes.push_back(hole_bd);
			}
		}
		std::cout<<holes.size()<<" holes were detected. \n";
	}

	void Packer::replace()
	{
		std::cout<<"Start replacing...\n";
		CGAL::Timer replace_timer;
		replace_timer.start();
		//std::ofstream res("parameters.txt");
		//rpvd.save_triangulation("pre_replace.obj");
		
#ifdef _CILK_
		cilk_for (unsigned int i = 0; i < pack_objects.size(); i++)
#else
		for (unsigned int i = 0; i < pack_objects.size(); i++)
#endif
		{
			const RestrictedPolygonVoronoiDiagram::VertGroup& samp_pnts = rpvd.sample_points_group(i);

			Local_frame lf = compute_local_frame(pack_objects[i]);

			std::vector<Segment_2> region2d;

			bool contain_non_delaunay_facet = false;
			bool penetration = false;
			
			for (unsigned int j = 0; j < samp_pnts.size(); j++)
			{
				const std::vector<Point_3>& bisec_pnts = samp_pnts[j]->vd_vertices;
				if (bisec_pnts.size() >= 2)
				{
					for (unsigned int k = 0; k < bisec_pnts.size()-1; k++)
					{
						const Point_3& s = bisec_pnts[k];
						const Point_3& t = bisec_pnts[k+1];
						region2d.push_back(Segment_2(lf.to_uv(s), lf.to_uv(t)));
					}
				}
				contain_non_delaunay_facet |= samp_pnts[j]->contain_non_delaunay_facet;
				penetration |= samp_pnts[j]->penetration;
			}
			Polygon_matcher pm(region2d, 200);
			//double region_area = pm.get_model_area();
			//double tile_area = pack_objects[i].area();
			//if (tile_area/region_area < 0.84)
			//{
			std::priority_queue<Match_info_item<unsigned int>, std::vector<Match_info_item<unsigned int>>, Match_measure> match_res;
			for (unsigned int idx = 0; idx < pgn_lib.size(); idx++)
				match_res.push(pm.affine_match(pgn_lib[idx], idx, match_weight));
			
			// choose the result with the smallest match error now
			Match_info_item<unsigned int> matcher = match_res.top();

			const Ex_polygon_2& match_pgn = pgn_lib[matcher.val];
			Polygon_2 transformed_pgn2d = CGAL::transform(matcher.t, match_pgn);

			pack_objects[i].clear();
			
			for (unsigned int j = 0; j < transformed_pgn2d.size(); j++)
			{
				Point_3 p = lf.to_xy(transformed_pgn2d.vertex(j));
				pack_objects[i].push_back(p);
			}

			Point_3 c = pack_objects[i].centroid();
			vec3 v;
			int fid;
			vec3 prjp = mesh.project_to_mesh(to_geex_pnt(c), v, fid);

			v = approx_normal(fid);
			pack_objects[i].align(to_cgal_vec(v), to_cgal_pnt(prjp));

			double shrink_factor;
			if (contain_non_delaunay_facet || penetration)
				shrink_factor = 0.3;
			else
				shrink_factor = 0.6;

			Transformation_3 rescalor = Transformation_3(CGAL::TRANSLATION, Vector_3(CGAL::ORIGIN, to_cgal_pnt(prjp))) *
										( Transformation_3(CGAL::SCALING, shrink_factor) *
										Transformation_3(CGAL::TRANSLATION, Vector_3(to_cgal_pnt(prjp), CGAL::ORIGIN)) );	

			std::transform(pack_objects[i].vertices_begin(), pack_objects[i].vertices_end(), pack_objects[i].vertices_begin(), rescalor);
			//
			pack_objects[i].lib_idx = matcher.val;
			pack_objects[i].factor = matcher.scale*shrink_factor;
			pack_objects[i].facet_idx = fid;
			pack_objects[i].texture_coord.assign(match_pgn.texture_coords.begin(), match_pgn.texture_coords.end());
			pack_objects[i].texture_id = match_pgn.texture_id;
		}
		replace_timer.stop();
		// rebuild the restricted delaunay triangulation and voronoi cell
		generate_RDT();
		compute_clipped_VD();
		stop_update_DT = false;
		std::cout<<"End replacing. Computation time: "<< replace_timer.time()<<" seconds.\n";
	}

	void Packer::enlarge_one_polygon(unsigned int id, double f, double theta, double tx, double ty)
	{
		// choose a local frame
		Local_frame lf;
		lf.o = pack_objects[id].centroid();
		lf.w = pack_objects[id].norm(); // assume this vector has already been normalized
		Vector_3 u(lf.o, pack_objects[id].vertex(0));
		lf.u = u/CGAL::sqrt(u.squared_length());
		lf.v = CGAL::cross_product(lf.w, lf.u);

		Parameter p(f, theta, tx, ty);

		transform_one_polygon(id, lf, p);
	}

	void Packer::curv_constrained_transform(Parameter& para, int fid, unsigned int pgn_id)
	{
		// approximate curvature at the facet
		const Facet& f = mesh[fid];
		double v0_cur = mesh.curvature_at_vertex(f.vertex_index[0]),
				v1_cur = mesh.curvature_at_vertex(f.vertex_index[1]),
				v2_cur = mesh.curvature_at_vertex(f.vertex_index[2]);
		//std::cout<<v0_cur<<", "<<v1_cur<<", "<<v2_cur<<std::endl;
		//system("pause");
		double avg_cur = (v0_cur + v1_cur + v2_cur) / 3.0 /*std::max(v0_cur, v1_cur)*/;
		//avg_cur = std::max(avg_cur, v2_cur);
		if (avg_cur == 0.0)
		{
			//system("pause");
			return;
		}
		double r = 1.0 / avg_cur;
		double squared_translation = para.tx * para.tx + para.ty * para.ty;
		double threshold = (2.0 - epsilon) * epsilon;
		if (squared_translation > threshold * r * r)
		{
			double temp = std::sqrt(threshold)*r/std::sqrt(squared_translation);
			//std::cout<<"scaler = "<<temp<<std::endl;
			para.tx *= temp;
			para.ty *= temp;
		}
		//double original_area = pack_objects[pgn_id].area();
		double pgn_radius2 = std::numeric_limits<double>::min();
		double mean_radius = 0.0;
		Point_3 c = pack_objects[pgn_id].centroid();
		for (unsigned int i = 0; i < pack_objects[pgn_id].size(); i++)
		{
			mean_radius += CGAL::sqrt(CGAL::squared_distance(c, pack_objects[pgn_id].vertex(i)));
			pgn_radius2 = std::max(pgn_radius2, CGAL::squared_distance(c, pack_objects[pgn_id].vertex(i)));
		}
		mean_radius /= pack_objects[pgn_id].size();
		if ( para.k * para.k * mean_radius * mean_radius/*pgn_radius2*/ > threshold * r * r)
		{
			//double temp = std::sqrt(threshold)*r/(para.k*std::sqrt(pgn_radius2));
			double temp = std::sqrt(threshold)*r/(para.k*mean_radius);
			para.k = std::max(1.0, para.k*temp);
			//std::cout<<"Curvature constrains scale factor.\n";
		}
	}

	void Packer::fill_holes()
	{
		std::cout<<"Start filling holes...\n";
		for (unsigned int i = 0; i < holes.size(); i++)
		{
			pack_objects.push_back(Packing_object());
			fill_one_hole(holes[i], pack_objects.back());
		}
		generate_RDT();
		stop_update_DT = false;
		compute_clipped_VD();
		std::cout<<"End filling holes.\n";
	}

	void Packer::fill_one_hole(Hole& hl, Packing_object& filler)
	{
		std::vector<Segment_2> hole_bd_2d;
		std::vector<Point_3> hole_verts;
		hole_verts.reserve(hl.size()*2);
		for (unsigned int j = 0; j < hl.size(); j++)
		{
			//const Segment_3& he = hl[j];
			hole_verts.push_back(hl[j].first->mp);
			hole_verts.push_back(hl[j].second->mp);
		}
		Point_3 hole_cent = CGAL::centroid(hole_verts.begin(), hole_verts.end(), CGAL::Dimension_tag<0>());
		vec3 v, prj_cent;
		int fid;
		prj_cent = mesh.project_to_mesh(to_geex_pnt(hole_cent), v, fid);
		Local_frame lf;
		lf.o = to_cgal_pnt(prj_cent);
		lf.w = to_cgal_vec(v);
		cgal_vec_normalize(lf.w);
		Plane_3 prj_plane(lf.o, lf.w);
		lf.u = prj_plane.base1();
		cgal_vec_normalize(lf.u);
		lf.v = CGAL::cross_product(lf.w, lf.u);
		for (unsigned int j = 0; j < hl.size(); j++)
		{
			//const Segment_3& he = hl[j];
			Point_2 s = lf.to_uv( prj_plane.projection(hl[j].first->mp) );
			Point_2 t = lf.to_uv( prj_plane.projection(hl[j].second->mp) );
			hole_bd_2d.push_back( Segment_2(s, t) );
		}
		Polygon_matcher pm(hole_bd_2d, 200);
		std::priority_queue<Match_info_item<unsigned int>, std::vector<Match_info_item<unsigned int>>, Match_measure> match_res;
		for (unsigned int idx = 0; idx < pgn_lib.size(); idx++)
			match_res.push(pm.affine_match(pgn_lib[idx], idx));
		// choose the result with the smallest match error now
		Match_info_item<unsigned int> matcher = match_res.top();

		const Ex_polygon_2& match_pgn = pgn_lib[matcher.val];
		Polygon_2 transformed_pgn2d = CGAL::transform(matcher.t, match_pgn);

		filler.clear();

		for (unsigned int j = 0; j < transformed_pgn2d.size(); j++)
		{
			Point_3 p = lf.to_xy(transformed_pgn2d.vertex(j));
			filler.push_back(p);
		}

		filler.align(lf.w, lf.o);
		double shrink_factor = 0.5;

		Transformation_3 rescalor = Transformation_3(CGAL::TRANSLATION, Vector_3(CGAL::ORIGIN, lf.o)) *
			( Transformation_3(CGAL::SCALING, shrink_factor) *
			Transformation_3(CGAL::TRANSLATION, Vector_3(lf.o, CGAL::ORIGIN)) );	

		std::transform(filler.vertices_begin(), filler.vertices_end(), filler.vertices_begin(), rescalor);
		filler.factor = matcher.scale * shrink_factor;
		filler.facet_idx = fid;
		filler.lib_idx = matcher.val;
		filler.texture_coord.assign(match_pgn.texture_coords.begin(), match_pgn.texture_coords.end());
		filler.texture_id = match_pgn.texture_id;	
	}

	void Packer::remove_one_polygon(unsigned int id, Hole& hole)
	{
		// build the hole after removing the polygon
		const RestrictedPolygonVoronoiDiagram::VertGroup& samp_pnts = rpvd.sample_points_group(id);
		typedef RDT_data_structure::Halfedge_around_vertex_circulator Edge_circulator;
		hole.clear();
		//std::list<Halfedge_handle> halfedge_bd;
		for (unsigned int i = 0; i < samp_pnts.size(); i++)
		{
			Edge_circulator start_edge = samp_pnts[i]->vertex_begin();
			Edge_circulator current_edge = start_edge;
			//bool penetration = true;
			do 
			{
				Vertex_handle v_adj = current_edge->opposite()->vertex();
				if (v_adj->group_id == samp_pnts[i]->group_id)
				{
					//penetration = false;
					break;
				}
				++current_edge;
			} while (current_edge != start_edge);
			//if (penetration)
				//std::cout<<"Caution: penetration may happened.! Polygon id: "<<id<<std::endl;
			start_edge = current_edge;
			do 
			{
				if (current_edge->facet() != Facet_handle() )
				{
					//Halfedge_handle opposite_edge = current_edge->prev();
					//Vertex_handle v0 = opposite_edge->vertex(), v1 = opposite_edge->opposite()->vertex();
					Vertex_handle v0 = current_edge->prev()->vertex(), v1 = current_edge->next()->vertex();
					if (v0->group_id != samp_pnts[i]->group_id && v1->group_id != samp_pnts[i]->group_id)
					{
						//halfedge_bd.push_back(opposite_edge);
						//Point_3 src = opposite_edge->vertex()->mp, tgt = opposite_edge->opposite()->vertex()->mp;
						//Vertex_handle src = opposite_edge->vertex(), tgt = opposite_edge->opposite()->vertex();
						hole.push_back(std::make_pair(v0, v1));
					}
				}
				++current_edge;
			} while (current_edge != start_edge);
		}

		// remove the corresponding triangulation
		std::set<Facet_handle> removed_facets;
		for (unsigned int i = 0; i < samp_pnts.size(); i++)
		{
			Edge_circulator start_edge = samp_pnts[i]->vertex_begin();
			Edge_circulator current_edge = start_edge;	
			bool border_break = false;
			do 
			{
				if (current_edge->is_border() || current_edge->opposite()->is_border())
				{
					border_break = true;
					break;
				}
				Facet_handle f = current_edge->facet();
				removed_facets.insert(f);
				++current_edge;
			} while (current_edge != start_edge);
			if (border_break)
			{
				//std::cout<<"border broken\n";
				current_edge = start_edge;
				do 
				{
					--current_edge;
					if (current_edge->is_border() || current_edge->opposite()->is_border())
					{
						border_break = true;
						break;
					}
					Facet_handle f = current_edge->facet();
					removed_facets.insert(f);			
				} while (current_edge != start_edge);
			}
		}
		for (std::set<Facet_handle>::iterator it = removed_facets.begin(); it != removed_facets.end(); ++it)
			rpvd.erase_facet((*it)->halfedge());
		rpvd.delete_point_group(id);
	}

	void Packer::remove_polygons()
	{
		static unsigned int id = 0;
		//for (unsigned int id = 0; id < pack_objects.size(); id++) 
		{
			holes.push_back(Hole());
			remove_one_polygon(id, holes.back());
		}

		id++;
	}

	bool Packer::replace_one_polygon(unsigned int id, Hole& region)
	{
		fill_one_hole(region, pack_objects[id]);


		// place the replacing polygon at a better place through optimization
		Local_frame lf = compute_local_frame(pack_objects[id]);
		Polygon_2 pgn2d;
		Plane_3 pln(lf.o, lf.w);
#if 0
		for (unsigned int i = 0; i < pack_objects[id].size(); i++)
			pgn2d.push_back(lf.to_uv(pack_objects[id].vertex(i)));
		std::vector<Segment_2> hole2d;
		hole2d.reserve(region.size());
		for (unsigned int i = 0; i < region.size(); i++)
		{
			Point_3 tan_src = pln.projection(region[i].first->mp);
			Point_3 tan_tgt = pln.projection(region[i].second->mp);
			hole2d.push_back(Segment_2(lf.to_uv(tan_src), lf.to_uv(tan_tgt)));
		}
#ifdef	_CILK_
		Containment& ctm = containments[id];
#else
		Containment& ctm = containment;
#endif
		ctm.load_polygons(pgn2d, hole2d);
		const unsigned int try_time_limit = 10;
		unsigned int try_times = 1;
		int knitro_res;
		Parameter solution;
		while ( (knitro_res = KTR_optimize(&solution.k, &solution.theta, &solution.tx, &solution.ty, id)) && try_times < try_time_limit )
		{
			double r = Numeric::random_float64();
			solution = Parameter(1.0, (rot_upper_bd-rot_lower_bd)*r+rot_lower_bd, 0.0, 0.0);
			try_times++;
		}

		if (try_times != try_time_limit)
		{
			transform_one_polygon(id, lf, solution);
		}
		lf = compute_local_frame(pack_objects[id]);
		pln = Plane_3(lf.o, lf.w);
#endif

		// re-triangulate the region nearby
		// triangulate the polygon first	
		//CDT cdt;
		cdt.clear();
		// sample the polygon
		double perimeter = 0.0;
		std::vector<double> edge_len(pack_objects[id].size());
		unsigned int nb_polygon_samp = 0; // for debug
		for ( unsigned int i = 0; i < pack_objects[id].size(); i++ )
		{
			Segment_3 e = pack_objects[id].edge(i);
			edge_len[i] = CGAL::sqrt(e.squared_length());
			perimeter += edge_len[i];
		}
		for ( unsigned int i = 0; i < pack_objects[id].size(); i++ )
		{
			Segment_3 e = pack_objects[id].edge(i);
			unsigned int n = edge_len[i] / perimeter * samp_nb;
			Point_3 src = e.source(), tgt = e.target();
			for ( unsigned int j = 0; j < n+1; j++ )
			{
				Point_3 p = CGAL::ORIGIN + ( (n+1-j)*(src - CGAL::ORIGIN) + j*(tgt - CGAL::ORIGIN) )/(n+1);
				//p = e.supporting_line().projection(p);
				vec3 dv;
				vec3 pp = mesh.project_to_mesh(to_geex_pnt(p), dv);
				CDT::Vertex_handle vh = cdt.insert(lf.to_uv(pln.projection(p)));
				vh->geo_info = MyPoint(p, to_cgal_pnt(pp), id);
				vh->already_exist_in_rdt = false;
				nb_polygon_samp++;
			}
		}
		// hole boundary segments as constraints
		std::map<Vertex_handle, CDT::Vertex_handle> hole_bd_pnts;
		for (unsigned int i = 0; i < region.size(); i++)
		{
			Vertex_handle s = region[i].first, t = region[i].second;

			if ( hole_bd_pnts.find(s) == hole_bd_pnts.end() )
			{
				CDT::Vertex_handle vh = cdt.insert(lf.to_uv(pln.projection(s->point())));
				vh->geo_info = MyPoint(s->point(),s->mp, s->group_id);
				vh->already_exist_in_rdt = true;
				vh->rdt_handle = s;
				hole_bd_pnts.insert(std::make_pair(s, vh));
			}
			if ( hole_bd_pnts.find(t) == hole_bd_pnts.end() )
			{
				CDT::Vertex_handle vh = cdt.insert(lf.to_uv(pln.projection(t->point())));
				vh->geo_info = MyPoint(t->point(), t->mp, t->group_id);
				vh->already_exist_in_rdt = true;
				vh->rdt_handle = t;
				hole_bd_pnts.insert(std::make_pair(t, vh));
			}
		}
		// add constraints
		std::vector<Segment_2> hole_edges;
		hole_edges.reserve(region.size());
		for (unsigned int i = 0;i < region.size(); i++)
		{
			CDT::Vertex_handle sh = hole_bd_pnts[region[i].first];
			CDT::Vertex_handle th = hole_bd_pnts[region[i].second];
			hole_edges.push_back(Segment_2(sh->point(), th->point()));
			cdt.insert_constraint(sh, th);
		}
		
		for (CDT::All_faces_iterator fit = cdt.all_faces_begin(); fit != cdt.all_faces_end(); ++fit)
		{
			if (cdt.is_infinite(fit))
			{
				fit->inside_hole = false;
				continue;
			}
			CDT::Vertex_handle v[] = {fit->vertex(0), fit->vertex(1), fit->vertex(2)};
			Point_2 c = CGAL::centroid(v[0]->point(), v[1]->point(), v[2]->point());
			if (inside_polygon(c, hole_edges))
				fit->inside_hole = true;
			else
				fit->inside_hole = false;
		}
		if (cdt.number_of_vertices() != (region.size() + nb_polygon_samp))
		{
			std::cout<<"assertion on number relationship failed!\n";
			std::cout<<"number of vertices in cdt: "<<cdt.number_of_vertices()<<std::endl;
			std::cout<<"number of vertices in polygons: "<<region.size() + nb_polygon_samp<<std::endl;
			//std::cout<<"now shrink relevant polygons...\n";
			//if (cdt.number_of_vertices() > (region.size() + nb_polygon_samp)) // constraint intersection
			//{
			//	std::set<unsigned int> pene_pgn_id;
			//	for (unsigned int k = 0; k < region.size(); k++)
			//		for (unsigned int l = k+1; l < region.size(); l++)
			//		{
			//			Segment_2 ek(hole_bd_pnts[region[k].first]->point(), hole_bd_pnts[region[k].second]->point());
			//			Segment_2 el(hole_bd_pnts[region[l].first]->point(), hole_bd_pnts[region[l].second]->point());
			//			if (CGAL::do_intersect(ek, el))
			//			{
			//				pene_pgn_id.insert(region[k].first->group_id);
			//				pene_pgn_id.insert(region[k].second->group_id);
			//				pene_pgn_id.insert(region[l].first>group_id);
			//				pene_pgn_id.insert(region[l].second->group_id);
			//			}
			//		}
			//	// shrink the relevant polygons
			//	Parameter shrink_trans(0.95, 0.0, 0.0, 0.0);
			//	for (std::set<unsigned int>::const_iterator it = pene_pgn_id.begin(); it != pene_pgn_id.end(); ++it)
			//	{
			//		Local_frame lf = compute_local_frame(pack_objects[*it]);
			//		transform_one_polygon(*it, lf, shrink_trans);
			//	}
			//}
			return false;
			//system("pause");
			//exit(0);
		}
		else
		{
			rpvd.delegate(CDTtoRDT(cdt, mesh, hole_bd_pnts));
			return true;
		}
	}

	void Packer::ex_replace()
	{
		//static unsigned int id = 0;
		//replace_one_polygon(id, holes.back());
		//id++;
		generate_RDT();
		eliminate_penetration();
		unsigned int id;
		for ( id = 0; id < pack_objects.size(); id++)
		{
			holes.push_back(Hole());
			backup = pack_objects[id];
			remove_one_polygon(id, holes.back());
			if (!replace_one_polygon(id, holes.back()))
				break;
		}
		if (id < pack_objects.size())
		{
			pack_objects[id] = backup;
			generate_RDT();
		}

	}
	void Packer::eliminate_penetration()
	{
		std::cout<<"Start eliminating penetration\n";
		for (unsigned int i = 0; i < pack_objects.size(); i++)
		{
			// collect neighboring polygons
			const RestrictedPolygonVoronoiDiagram::VertGroup& samp_pnts = rpvd.sample_points_group(i);
			std::set<unsigned int> nid;
			for (unsigned int j = 0; j < samp_pnts.size(); j++)
			{
				Halfedge_handle e = samp_pnts[j]->halfedge();
				Halfedge_handle end = e;
				do
				{
					nid.insert(e->opposite()->vertex()->group_id);
					e = e->next()->opposite();
				} while ( e != end);
			}
			bool penetration_detected = true;
			while (penetration_detected)
			{
				penetration_detected = false;
				for (std::set<unsigned int>::const_iterator it = nid.begin(); it != nid.end(); ++it)
					if ( *it != i)
					{
						if (pair_penetration(i, *it))
						{
							std::cout<<"Penetration or overlap detected!\n";
							penetration_detected = true;
							break;
						}
					}
				if (penetration_detected)
				{
					Local_frame lf = compute_local_frame(pack_objects[i]);
					Parameter shrink(0.90, 0.0, 0.0, 0.0);
					transform_one_polygon(i, lf, shrink);
					rpvd.iDT_update();
				}				
			}
		}
		std::cout<<"Eliminating penetration ends.\n";
	}

	bool Packer::pair_penetration(unsigned int id0, unsigned int id1)
	{
		Local_frame lf = compute_local_frame(pack_objects[id0]);
		Polygon_2 pgn0_2d;
		std::vector<Point_2> pgn1_2d;
		const RestrictedPolygonVoronoiDiagram::VertGroup& samp_pnts = rpvd.sample_points_group(id1);
		pgn1_2d.reserve(samp_pnts.size());
		for (unsigned int i = 0; i < pack_objects[id0].size(); i++)
			pgn0_2d.push_back(lf.to_uv(pack_objects[id0].vertex(i)));
		Plane_3 pln(lf.o, lf.w);
		for (unsigned int i = 0; i < samp_pnts.size(); i++)
			pgn1_2d.push_back(lf.to_uv(pln.projection(samp_pnts[i]->point())));
		// check whether id1 overlaps with id0
		for (unsigned int i = 0; i < pgn1_2d.size(); i++)
		{
			if (!pgn0_2d.has_on_unbounded_side(pgn1_2d[i]))
				return true;
		}
		return false;
	}
	void Packer::save_curvature_and_area()
	{
		std::ofstream cur_area_file("cur_area.txt");
		for (unsigned int i = 0; i < pack_objects.size(); i++)
		{
			double pgn_area = pack_objects[i].area();
			int fid = pack_objects[i].facet_idx;
			double cur = (mesh.curvature_at_vertex(mesh[fid].vertex_index[0])
						+ mesh.curvature_at_vertex(mesh[fid].vertex_index[1])
						+ mesh.curvature_at_vertex(mesh[fid].vertex_index[2]) ) / 3.0;
			cur_area_file << pgn_area << "\t" << cur<<std::endl;
		}
	}
	int Packer::KTR_optimize(double* io_k, double* io_theta, double* io_t1, double* io_t2, unsigned int idx)
	{
		int  nStatus;
		/* variables that are passed to KNITRO */
		KTR_context *kc;
		int n, m, nnzJ, nnzH, objGoal, objType;
		int *cType, *jacIndexVars, *jacIndexCons;
		double obj, *x, *lambda;
		double *xLoBnds, *xUpBnds, *xInitial, *cLoBnds, *cUpBnds;
		int i, j, k; // convenience variables

		/*problem size and mem allocation */
		n = 4;
#ifdef _CILK_
		m = containments[idx].get_constaint_list_size();
#else
		m = containment.get_constaint_list_size();
#endif
		nnzJ = n*m;
		nnzH = 0;
		x = (double *)malloc(n * sizeof(double));
		lambda = (double *)malloc((m+n) * sizeof(double));

		xLoBnds      = (double *)malloc(n * sizeof(double));
		xUpBnds      = (double *)malloc(n * sizeof(double));
		xInitial     = (double *)malloc(n * sizeof(double));
		cType        = (int    *)malloc(m * sizeof(int));
		cLoBnds      = (double *)malloc(m * sizeof(double));
		cUpBnds      = (double *)malloc(m * sizeof(double));
		jacIndexVars = (int    *)malloc(nnzJ * sizeof(int));
		jacIndexCons = (int    *)malloc(nnzJ * sizeof(int));

		/* objective type */
		objType = KTR_OBJTYPE_LINEAR;
		objGoal = KTR_OBJGOAL_MAXIMIZE;
		/* bounds and constraints type */
		xLoBnds[0] = 0.0;
		xUpBnds[0] = KTR_INFBOUND;
		xLoBnds[1] = rot_lower_bd;
		xUpBnds[1] = rot_upper_bd;
		for (i = 2; i < n; i++) {
			xLoBnds[i] = -KTR_INFBOUND;
			xUpBnds[i] = KTR_INFBOUND;
		}
		for (j = 0; j < m; j++) {
			cType[j] = KTR_CONTYPE_GENERAL;
			cLoBnds[j] = 0.0;
			cUpBnds[j] = KTR_INFBOUND;
		}
		/* initial point */
		xInitial[0] = *io_k;
		xInitial[1] = *io_theta;
		xInitial[2] = *io_t1;
		xInitial[3] = *io_t2;
		/* sparsity pattern (here, of a full matrix) */
		k = 0;
		for (j = 0; j < m; j++) 
			for (i = 0; i < n; i++){
				jacIndexCons[k] = j;
				jacIndexVars[k] = i;
				k++;
			}
		/* create a KNITRO instance */
		kc = KTR_new();
		if (kc == NULL)
			return -1 ; // probably a license issue
		/* set options: automatic gradient and hessian matrix */
		if (KTR_set_int_param_by_name(kc, "gradopt", KTR_GRADOPT_EXACT) != 0)
			return -1 ;
		if (KTR_set_int_param_by_name(kc, "hessopt", KTR_HESSOPT_LBFGS) != 0)
			return -1 ;
		if (KTR_set_int_param_by_name(kc, "outlev", 0) != 0)
			return -1 ;
		if (KTR_set_int_param_by_name(kc, "maxit", 500) != 0)
			return -1 ;
		/* register the callback function */
		if (KTR_set_func_callback(kc, &callback) != 0)
			return -1 ;
		if (KTR_set_grad_callback(kc, &callback) != 0)
			return -1 ;
		/* pass the problem definition to KNITRO */
		nStatus = KTR_init_problem(kc, n, objGoal, objType, xLoBnds, xUpBnds, m, cType, cLoBnds, cUpBnds,
									nnzJ, jacIndexVars, jacIndexCons, nnzH, NULL, NULL, xInitial, NULL);
		/* free memory (KNITRO maintains its own copy) */
		free(xLoBnds);	free(xUpBnds);
		free(xInitial);
		free(cType);
		free(cLoBnds);	free(cUpBnds);
		free(jacIndexVars);	free(jacIndexCons);

#ifdef _CILK_
		nStatus = KTR_solve(kc, x, lambda, 0, &obj, NULL, NULL, NULL, NULL, NULL, &idx);
#else
		nStatus = KTR_solve(kc, x, lambda, 0, &obj, NULL, NULL, NULL, NULL, NULL, NULL);
#endif

		if(nStatus ==0)
		{	
			*io_k =  x[0];
			*io_theta = x[1];
			*io_t1 = x[2];
			*io_t2 = x[3];
		}

		/* delete the KNITRO instance and primal/dual solution */
		KTR_free(&kc);
		free(x);
		free(lambda);
		return nStatus;
	}
	int Packer::callback(const int evalRequestCode, const int n, const int m, const int nnzJ, const int nnzH,
						const double * const x,	const double * const lambda, double * const obj, double * const c,
						double * const objGrad,	double * const jac,	double * const hessian,	double * const hessVector, void * userParams)
	{
		Packer* p = instance() ;
		if (evalRequestCode == KTR_RC_EVALFC)
		{
			/* objective function */
			// x[0]- x[4]: k cos sin t1 t2
			/*************************** No weighted objective ****************************/
			*obj    = x[0];
			/**************************** Weighted objective ***********************/
			//*obj = p->containments[*(unsigned int*)userParams].compute_extended_object(x[0], x[2], x[3]);
#ifdef _CILK_
			p->containments[*(unsigned int*)userParams].compute_constraints(x, c, m);
#else
			p->containment.compute_constraints(x,c,m);
#endif
			return 0;
		}
		else if (evalRequestCode == KTR_RC_EVALGA)
		{
			/*************************** No weighted objective ****************************/
			objGrad[0] = 1.0;
			objGrad[1] = 0.0;
			objGrad[2] = 0.0;
			objGrad[3] = 0.0;
			/**************************** Weighted objective ***********************/
			//p->containments[*(unsigned int*)userParams].compute_extended_object_grad(x[0], x[2], x[3], &objGrad[0], &objGrad[2], &objGrad[3]);
			//objGrad[1] = 0.0;

#ifdef _CILK_
			p->containments[*(unsigned int*)userParams].compute_constraint_grads(x, jac);
#else
			p->containment.compute_constraint_grads(x, jac);
#endif
			return 0;
		}
		else 
		{
			printf ("Wrong evalRequestCode in callback function.\n");
			return -1;
		}
	}

	void Packer::get_bbox(real& x_min, real& y_min, real& z_min, real& x_max, real& y_max, real& z_max)
	{
		x_min = y_min = z_min = DBL_MAX;
		x_max = y_max = z_max = DBL_MIN;
		if (!pio.multi_meshes_specified())
			for(unsigned int i=0; i<mesh.size(); i++) 
				for(unsigned int j=0; j<3; j++) 
				{
					const vec3& p = mesh[i].vertex[j] ;
					x_min = gx_min(x_min, p.x); y_min = gx_min(y_min, p.y);	z_min = gx_min(z_min, p.z);
					x_max = gx_max(x_max, p.x);	y_max = gx_max(y_max, p.y); z_max = gx_max(z_max, p.z) ;
				}
		else
		{
			for (unsigned int i = 0; i < mesh_segments.size(); i++)
				for (unsigned int j = 0; j < mesh_segments[i].size(); j++)
					for (unsigned int k = 0; k < 3; k++)
					{
						const vec3& p = mesh_segments[i][j].vertex[k] ;
						x_min = gx_min(x_min, p.x); y_min = gx_min(y_min, p.y);	z_min = gx_min(z_min, p.z);
						x_max = gx_max(x_max, p.x);	y_max = gx_max(y_max, p.y); z_max = gx_max(z_max, p.z) ;
					}
		}
	}
}