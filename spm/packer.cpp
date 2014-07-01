#include "packer.h"

namespace Geex
{
	Packer* Packer::instance_ = nil;

	double Packer::PI = 3.141592653589793;

	Packer::Packer()
	{
		assert( instance_ == nil );
		instance_ = this;

		srand(time(NULL));
		rot_lower_bd = -PI/6.0;
		rot_upper_bd = PI/6.0;

		frontier_edge_size = 0.4;
		hole_face_size = 0.4;

		stop_update_DT = false;

		epsilon = 0.15;
		match_weight = 0.0;

		samp_nb = 20;

		sub_pack_id = 0;

		sync_opt = true;
		phony_upper_scale = 100.0;
		use_voronoi_cell_ = true;

		replace_factor = 0.8;

		vector_field = false;
		align_constraint = false;

		only_allow_scale = false;
	}

	void Packer::load_project(const std::string& prj_config_file)
	{
		pio.load_project(prj_config_file);
		std::cout<<"\""<<pio.attribute_value("SPMProject")<<"\" project loaded.\n";
		if (pio.mesh_polygon_coupled())
		{
			std::cout<<"reading pairs of polygons and meshes...\n";
			pio >> pgn_lib_sets;
			pio >> mesh_segments;
		}
		else
		{
			if (!pio.multi_meshes_specified())
			{
				pio>>mesh>>pgn_lib;
				/*std::cout<<"filter the polygons\n";
				std::map<double, Ex_polygon_2> area_polygon;
				for (size_t k = 0; k < pgn_lib.size(); k++)
					area_polygon[std::fabs(pgn_lib[k].area())] = pgn_lib[k];
				pgn_lib.clear();
				std::map<double, Ex_polygon_2>::iterator it = area_polygon.begin();
				std::advance(it, 80);
				for (; it != area_polygon.end(); ++it)
					pgn_lib.push_back(it->second);
				std::ofstream ofs("voronoi_207.poly");
				ofs << pgn_lib.size() << std::endl;
				for (size_t k = 0; k < pgn_lib.size(); k++)
				{
					ofs << pgn_lib[k].size() << std::endl;
					for (size_t l = 0; l < pgn_lib[k].size(); l++)
						ofs << pgn_lib[k].vertex(l).x()<<' '<<pgn_lib[k].vertex(l).y()<<std::endl;
				}*/
				if (!pio.attribute_value("VectorField").empty())
					mesh.load_vecfield(pio.attribute_value("VectorField"));
				vector_field = mesh.with_vector_field();
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
		}
		double min_val, max_val;
		int nb_levels;
		if (pio.read_physical_scales(min_val, max_val, nb_levels))
		{
			phy_disc_barr.set(min_val, max_val, nb_levels);
			std::cout<<"Physical: < "<<min_val<<", "<<max_val<<", "<<nb_levels<<'>'<<std::endl;
		}
		if (pio.read_optimize_scales(min_val, max_val, nb_levels))
		{
			opt_disc_barr.set(min_val, max_val, nb_levels);
			std::cout<<"Auxiliary: < "<<min_val<<", "<<max_val<<", "<<nb_levels<<'>'<<std::endl;
		}

		if (phy_disc_barr.nb_levels() != 0 && opt_disc_barr.nb_levels() != 0)
			disc_barr.set(phy_disc_barr, opt_disc_barr);
	}

	void Packer::pack_next_submesh() // start a new packing process
	{
		//static unsigned int mesh_id = 0;
		if (mesh_segments.size() > 0)
		{
			mesh.clear();
			pack_objects.clear();
			stop_update_DT = false;
			mesh = mesh_segments[sub_pack_id];
		}

		// if polygon is coupled with mesh
		if (pgn_lib_sets.size() > 0)
		{
			pgn_lib.clear();
			pgn_lib.reserve(pgn_lib_sets[sub_pack_id].size());
			for (unsigned int i = 0; i < pgn_lib_sets[sub_pack_id].size(); i++)
				pgn_lib.push_back(pgn_lib_sets[sub_pack_id][i]);
		}
		initialize();
		std::cout<<"start another packing process...\n";
		generate_RDT();
		compute_clipped_VD();
	}

	void Packer::save_sub_result()
	{
		//static unsigned int mesh_id = 0;
		if (pack_objects.size() == 0)
			return;
		res_pack_objects.push_back(std::vector<Packing_object>());
		res_pack_objects.back().reserve(pack_objects.size());
		for (unsigned int i = 0; i < pack_objects.size(); i++)
			res_pack_objects.back().push_back(pack_objects[i]);
		sub_pack_id++;	
	}

	void Packer::initialize()
	{
		mesh.build_kdtree();

		//build similarity matrix
		sim_mat.resize(pgn_lib.size());
		std::for_each( sim_mat.begin(), sim_mat.end(), std::bind2nd(std::mem_fun1_ref(&std::vector<double>::resize), pgn_lib.size()) );

//#ifdef _CILK_ 
//		cilk_for (unsigned int i = 0; i < pgn_lib.size(); i++)
//#else 
		for (unsigned int i = 0; i < pgn_lib.size(); i++)
//#endif 		
		{
			Polygon_matcher pm(pgn_lib[i]);
			for (unsigned int j = 0; j < pgn_lib.size(); j++)
			{
				if (j == i)
					sim_mat[i][j] = 0.0;
				else
				{
					Match_info_item<unsigned int> match_res = pm.affine_match(pgn_lib[j], j, 0.0);
					sim_mat[i][j] = match_res.error;
				}
			}
		}
		
		// compute mesh area
		mesh_area = 0.0;
		for (unsigned int i = 0; i < mesh.size(); i++)
			mesh_area += std::fabs(mesh[i].area());
		std::cout<<"Surface area: "<<mesh_area<<std::endl;
		double max_cur = std::numeric_limits<double>::min();
		double min_cur = std::numeric_limits<double>::max();
		for (unsigned int i = 0; i < mesh.size(); i++)
		{
			double cur = mesh.curvature_at_face(i);
			max_cur = std::max(max_cur, cur);
			min_cur = std::min(min_cur, cur);
		}
		std::cout<<"Maximum curvature = "<<max_cur<<", minimum curvature = "<<min_cur<<std::endl;
		
		// compute maximum and average polygon area in the library
		double max_pgn_area = -std::numeric_limits<double>::max(), mean_pgn_area = 0.0;
		for (unsigned int i = 0; i < pgn_lib.size(); i++)
		{
			double s = std::fabs(pgn_lib[i].area());
			max_pgn_area = std::max(max_pgn_area, s);
			mean_pgn_area += s;
		}
		mean_pgn_area /= pgn_lib.size();
		std::cout<<"Maximum polygon area: "<<max_pgn_area<<std::endl;
		std::cout<<"Mean polygon area: "<<mean_pgn_area<<std::endl;
		
		unsigned int nb_tiles = 0;
		if (!pio.attribute_value("NumberTiles").empty())
		{
			std::istringstream iss(pio.attribute_value("NumberTiles"));
			iss >> nb_tiles;
		}
		if (mesh_segments.size() == 0)
			rpvd.set_mesh(pio.attribute_value("MeshFile"));
		else
			rpvd.set_mesh(pio.get_submesh_files().at(sub_pack_id));
		rpvd.set_trimesh(&mesh);
		
		if (nb_tiles > 0)
		{
			double adjust_factor = mesh_area/nb_tiles/mean_pgn_area;
			mean_pgn_area *= adjust_factor;
			adjust_factor = std::sqrt(adjust_factor);
			Transformation_2 t(CGAL::SCALING, adjust_factor);
			for (unsigned int i = 0; i < pgn_lib.size(); i++)
				pgn_lib[i] *= adjust_factor;
		}
		// distribute polygons by different strategies
		// 1. random 
		unsigned int nb_init_polygons;
		float expected_coverage = 0.7f;
		if (!pio.attribute_value("ExpectedCoverage").empty())
		{
			std::istringstream iss(pio.attribute_value("ExpectedCoverage"));
			iss >> expected_coverage;
		}
		if (nb_tiles > 0)
			nb_init_polygons = nb_tiles;
		else
		{
			std::cout<<"Expected Coverage = "<<expected_coverage<<std::endl;
			nb_init_polygons = mesh_area / mean_pgn_area * expected_coverage;
		}
		std::cout<<nb_init_polygons<<" polygons are expected.\n";
		random_init_tiles(nb_init_polygons);
		std::cout<<pack_objects.size()<<" polygons are loaded initially.\n";

		activate_all();
	}
	void Packer::random_init_tiles(unsigned int nb_init_polygons)
	{
		//std::ofstream pos_ofs("init_pos.txt");
		//std::ofstream fid_ofs("init_fid.txt");
		//Init_point_generator *ipg = new CVT_point_generator(mesh, nb_init_polygons, pio.attribute_value("MeshFile"));

		/* for rigid size version, curvature sensitive initialization */
		std::vector<double> vtx_cur(mesh.nb_vertices());
		for (size_t i = 0; i < mesh.nb_vertices(); i++)
			vtx_cur[i] = -mesh.curvature_at_vertex(i);
		std::sort(vtx_cur.begin(), vtx_cur.end());
		std::map<double, size_t> area_idx;
		for (size_t i = 0; i < pgn_lib.size(); i++)
			area_idx[std::fabs(pgn_lib[i].area())] = i;
		/**************************************************************/

		Init_point_generator *ipg = new Density_point_generator(mesh, nb_init_polygons);
		// choose polygons and distribute
		unsigned int pgn_lib_idx = 0;
		unsigned int nb_pnts = ipg->nb_points();
		pack_objects.reserve(nb_pnts);
		//pos_ofs << nb_pnts << std::endl;
		//fid_ofs << nb_pnts << std::endl;
		for (unsigned int i = 0; i < nb_pnts; i++)
		{
			Point_3 pos = ipg->next_point();
			int f_idx;
			vec3 n;
			mesh.project_to_mesh(to_geex_pnt(pos), n, f_idx);

			/* for rigid size version, curvature sensitive initialization */
			double fcur = mesh.curvature_at_face(f_idx);
			std::vector<double>::iterator cur_it = std::lower_bound(vtx_cur.begin(), vtx_cur.end(), -fcur);
			size_t cur_rank = cur_it - vtx_cur.begin();
			std::map<double, size_t>::iterator it = area_idx.begin();
			std::advance(it, float(cur_rank)/vtx_cur.size()*pgn_lib.size());
			if (it == area_idx.end())	--it;
			pgn_lib_idx = it->second;
			/**************************************************************/			

			const Facet& f = mesh[f_idx];
			vec3 gx_normal = approx_normal(f_idx);
			const Ex_polygon_2& pgn_2 = pgn_lib[pgn_lib_idx%pgn_lib.size()];
			// shrink factor
			double fa = std::fabs(f.area()), pa = std::fabs(pgn_2.area());
			double s = std::min(fa/pa, pa/fa);
			s = 0.1*std::sqrt(s);
			Polygon_2 init_polygon = CGAL::transform(Transformation_2(CGAL::SCALING, s), pgn_2);
			pack_objects.push_back(Packing_object(init_polygon, to_cgal_vec(gx_normal), pos, s));
			pack_objects.back().lib_idx = pgn_lib_idx%pgn_lib.size();
			pack_objects.back().facet_idx = f_idx;
			pack_objects.back().texture_id = pgn_2.texture_id;
			pack_objects.back().texture_coord.assign(pgn_2.texture_coords.begin(), pgn_2.texture_coords.end());

			pgn_lib_idx++;

			//pos_ofs << pos.x() << ' ' << pos.y() << ' ' << pos.z() << std::endl;
			//fid_ofs << f_idx << std::endl;
		}
		delete ipg;
	}
	
	void Packer::adjust(double factor)
	{
		if (factor <= 0.0)
			return;
		for (unsigned int i = 0; i < pack_objects.size(); i++)
		{
			if (pack_objects[i].active)
			{
				Local_frame lf = pack_objects[i].local_frame();
				Parameter scale(factor, 0.0, 0.0, 0.0);
				transform_one_polygon(i, lf, scale);
			}
		}
		generate_RDT();
		compute_clipped_VD();

		if (sync_opt)
		{
			double min_size = std::numeric_limits<double>::max();
			for (unsigned int i = 0; i < pack_objects.size(); i++)
			{				
				//min_size = std::min(pack_objects[i].rel_factor(mesh.curvature_at_face(pack_objects[i].facet_idx)), min_size);
				min_size = std::min(min_size, pack_objects[i].factor);
			}
			disc_barr.set_current_barrier(min_size);
		}
	}

	void Packer::auto_adjust()
	{
		double min_scale = std::numeric_limits<double>::max();
		for (unsigned int i = 0; i < pack_objects.size(); i++)
			min_scale = std::min(min_scale, pack_objects[i].factor);
		for (unsigned int i = 0; i < pack_objects.size(); i++)
		{
			Local_frame lf = pack_objects[i].local_frame();
			Parameter scale(min_scale/pack_objects[i].factor, 0.0, 0.0, 0.0);
			transform_one_polygon(i, lf, scale);
		}
		generate_RDT();
		compute_clipped_VD();

		if (sync_opt)
		{
			double min_size = std::numeric_limits<double>::max();
			for (unsigned int i = 0; i < pack_objects.size(); i++)			
				//min_size = std::min(pack_objects[i].rel_factor(mesh.curvature_at_face(pack_objects[i].facet_idx)), min_size);
				min_size = std::min(min_size, pack_objects[i].factor);
			disc_barr.set_current_barrier(min_size);
		}

		//activate_all();
	}
	void Packer::compute_clipped_VD(bool approx)
	{
		//std::cout<<"Start computing clipped Voronoi region\n";
		std::vector<Plane_3> tangent_planes;
		tangent_planes.reserve(pack_objects.size());
		std::vector<Point_3> ref_pnts;
		ref_pnts.reserve(pack_objects.size());
		for (unsigned int i = 0; i < pack_objects.size(); i++)
		{
			Point_3 c = pack_objects[i].centroid();
			Vector_3 n = pack_objects[i].norm();
			tangent_planes.push_back(Plane_3(c, n));
			ref_pnts.push_back(c);
		}
		rpvd.compute_clipped_VD(tangent_planes, ref_pnts);
		rpvd.compute_midpoint_VD(tangent_planes, ref_pnts);
		//std::cout<<"End computing.\n";
	}


	void Packer::generate_RDT()
	{
		rpvd.begin_insert();
		rpvd.insert_polygons(pack_objects.begin(), pack_objects.end(), samp_nb);
		//rpvd.insert_polygons(pack_objects.begin(), pack_objects.end());
		rpvd.insert_bounding_points(10);
		CGAL::Timer t;
		t.start();
		rpvd.end_insert();
		t.stop();
		std::cout<<"time spent in RDT: "<<t.time()<<" seconds."<<std::endl;
		std::for_each(holes.begin(), holes.end(), std::mem_fun_ref(&Hole::clear));
		holes.clear();
	}

	vec3 Packer::approx_normal(unsigned int facet_idx)
	{
		const Facet& f = mesh[facet_idx];
		vec3 original_normal = f.normal();
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

	Packer::Optimization_res Packer::optimize_one_polygon(unsigned int id, Local_frame& lf, Parameter& solution)
	{
		// choose a local frame
		lf = pack_objects[id].local_frame();

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
		if (use_voronoi_cell_)
		{
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
					}
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
						}
					}
				}	
		}
		else
		{
   			for (unsigned int i = 0; i < samp_pnts.size(); i++)
   			{
   				const vector<Segment_3>& med_segs = samp_pnts[i]->med_segs;
   				Point_2 ref_point = lf.to_uv(samp_pnts[i]->point());
   				if (med_segs.size() == 0)
   					continue;
   				for (unsigned int j = 0; j < med_segs.size(); j++)
   				{
   					Segment_2 s2d(lf.to_uv(med_segs[j]));
   					ctm.const_list_.push_back(Constraint(/*lf.to_uv(samp_pnts[i]->point())*/ref_point, s2d, ref_point));
   				}
   			}
		}

 		if (vector_field)
 			ctm.align_angle = vecfield_align_angle(id);

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
		{
			solution.k = 1.0;
			solution.theta = solution.tx = solution.ty = 0.0;
			return FAILED;
		}
		else
			return SUCCESS;
	}

	Packer::Lloyd_res Packer::one_lloyd(bool enlarge, std::vector<Parameter>& solutions, std::vector<Local_frame>& lfs)
	{
  		if (vector_field && !align_constraint)
  			vecfield_align();

		double min_scalor;

		min_scalor = 1.05;
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
// 			if (solutions[i].k < 1.0)
// 			{
// 				solutions[i].k = 1.0;
// 				solutions[i].theta = solutions[i].tx = solutions[i].ty = 0.0;
// 			}
		}

		// do statistics on transformation
		double mink = std::numeric_limits<double>::max();
		for (unsigned int i = 0; i < pack_objects.size(); i++)
		{
			if (opti_res[i] == SUCCESS && pack_objects[i].active)
				mink = std::min(solutions[i].k, mink);
			else if (opti_res[i] == SUCCESS && !pack_objects[i].active)
				solutions[i].k = 1.0;
			else
			{
				solutions[i].k = 1.0;
				solutions[i].theta = solutions[i].tx = solutions[i].ty = 0.0;
			}
		}
		if (mink >= min_scalor)
		{
			for (unsigned int i = 0; i < pack_objects.size(); i++)
				if (opti_res[i] == SUCCESS && pack_objects[i].active)
					solutions[i].k = mink;
		}
		std::cout<<"Minimum enlargement factor = "<<mink<<std::endl;

		//if (!vector_field)
			constraint_transformation(solutions, lfs, false);
			constraint_transformation(solutions, lfs, true);

#ifdef _CILK_
		cilk_for (unsigned int i = 0; i < pack_objects.size(); i++)
#else
		for (unsigned int i = 0; i < pack_objects.size(); i++)
#endif
		{
			curv_constrained_transform(solutions[i], pack_objects[i].facet_idx, i);
			if (pack_objects[i].active)
				transform_one_polygon(i, lfs[i], solutions[i]);
		}

		// update facet normal
		for (RestrictedPolygonVoronoiDiagram::Facet_iterator fit = rpvd.faces_begin(); fit != rpvd.faces_end(); ++fit)
		{
			RestrictedPolygonVoronoiDiagram::Halfedge_handle eh = fit->halfedge();
			Point_3 v0 = eh->vertex()->mp;
			eh = eh->next();
			Point_3 v1 = eh->vertex()->mp;
			eh = eh->next();
			Point_3 v2 = eh->vertex()->mp;
			Vector_3 v01(v0, v1), v12(v1, v2);
			cgal_vec_normalize(v01);
			cgal_vec_normalize(v12);
			fit->n = CGAL::cross_product(v01, v12);
			fit->n = fit->n / CGAL::sqrt(fit->n.squared_length());
		}
		if (enlarge && mink < min_scalor)
			return NO_MORE_ENLARGEMENT;
		else
			return MORE_ENLARGEMENT;
	}

	bool Packer::discrete_one_lloyd(bool enlarge, std::vector<Parameter>& solutions, std::vector<Local_frame>& lfs, 
									double barrier_factor, double nxt_barrier)
	{
  		if (vector_field && !align_constraint)
  			vecfield_align();

		double min_factor = 1.005; // to determine convergence

		std::vector<Optimization_res> opti_res(pack_objects.size());
#ifdef _CILK_
		containments.resize(pack_objects.size());
		cilk_for (unsigned int i = 0; i < pack_objects.size(); i++)
#else
		for (unsigned int i = 0; i < pack_objects.size(); i++)
#endif
		{
			opti_res[i] = optimize_one_polygon(i, lfs[i], solutions[i]);
		}
		for (size_t i = 0; i < pack_objects.size(); i++)
			if (!pack_objects[i].active)
			{
				solutions[i].k = 1.0;
				solutions[i].theta = solutions[i].tx = solutions[i].ty = 0.0;
			}
		// check whether the tile has already been close to its voronoi region
		for (unsigned int i = 0; i < pack_objects.size(); i++)
		{
			if (!enlarge)
				solutions[i].k = 1.0;
		}
	
		double mink = std::numeric_limits<double>::max();
		int nb_failures = 0;
		// suppress the inactive and optimization failure ones
		//for (unsigned int i = 0; i < pack_objects.size(); i++)
		//{
		//	if (opti_res[i] == SUCCESS && pack_objects[i].active)
		//		mink = std::min(solutions[i].k, mink);
		//	else if (opti_res[i] == SUCCESS && !pack_objects[i].active)
		//		solutions[i].k = 1.0;
		//	else if (opti_res[i] != SUCCESS)
		//	{
		//		nb_failures++;
		//		solutions[i].k = 1.0;
		//		solutions[i].theta = solutions[i].tx = solutions[i].ty = 0.0;
		//	}
		//}

		//std::cout<<"number of failures: "<<nb_failures<<std::endl;
		
		if (!only_allow_scale)
			constraint_transformation(solutions, lfs, false);
		else
			for (unsigned int i = 0; i < solutions.size(); i++)
				solutions[i].theta = solutions[i].tx = solutions[i].ty = 0.0;
		//if (!vector_field)
		constraint_transformation(solutions, lfs, true);

//  #ifdef _CILK_
//  		cilk_for (unsigned int i = 0; i < pack_objects.size(); i++)
//  			curv_constrained_transform(solutions[i], pack_objects[i].facet_idx, i);
//  #else
//  		for (unsigned int i = 0; i < pack_objects.size(); i++)
//  			curv_constrained_transform(solutions[i], pack_objects[i].facet_idx, i);
//  #endif
		// after a series of filtering, compute the minimum scale factor
		for (unsigned int i = 0; i < pack_objects.size(); i++)
		{
			if (opti_res[i] == SUCCESS && pack_objects[i].active)
				mink = std::min(solutions[i].k, mink);
			//else if (opti_res[i] == SUCCESS && !pack_objects[i].active)
			//	solutions[i].k = 1.0;
		}
		for (size_t i = 0; i < pack_objects.size(); i++)
			if (opti_res[i] == SUCCESS && pack_objects[i].active)
				solutions[i].k = mink;
		std::cout<<"Minimum scale is "<<mink<<std::endl;
		std::cout<<"number of failures: "<<nb_failures<<std::endl;

		// check state of each tile
		std::vector<bool> has_growth_room(pack_objects.size());
		if (enlarge)
			for (unsigned int i = 0; i < pack_objects.size(); i++)
			{
				//pack_objects[i].reach_barrier = false;
				//double cur = mesh.curvature_at_face(pack_objects[i].facet_idx);
				//double rel_factor = pack_objects[i].rel_factor(cur);
				
				if (pack_objects[i].active)
				{
					//if (pack_objects[i].factor < nxt_barrier)
					if (pack_objects[i].factor*solutions[i].k > barrier_factor)
					{
						solutions[i].k = barrier_factor / pack_objects[i].factor; // restrict it at this barrier
					//if (rel_factor < nxt_barrier)
					//	solutions[i].k = barrier_factor / rel_factor;
					//else
					//	solutions[i].k = 1.0;
						has_growth_room[i] = false;
					}
					else
						has_growth_room[i] = true;
					//pack_objects[i].reach_barrier = true;
				}
				//else if (pack_objects[i].active && solutions[i].k < min_factor)
				//	has_growth_room[i] = false;
				else// if (pack_objects[i].active && solutions[i].k >= min_factor)
					has_growth_room[i] = false;
				//else 
				//	has_growth_room[i] = false;
			}		
		//else
		//	for (unsigned int i = 0; i < pack_objects.size(); i++)
		//		pack_objects[i].reach_barrier = true;
#ifdef _CILK_
		cilk_for (unsigned int i = 0; i < pack_objects.size(); i++)
#else
		for (unsigned int i = 0; i < pack_objects.size(); i++)
#endif
			if (pack_objects[i].active)
				transform_one_polygon(i, lfs[i], solutions[i]);


		bool stay_at_this_barrier = false;
		size_t nb_potential_growth = std::count(has_growth_room.begin(), has_growth_room.end(), true);
		float perc = nb_potential_growth/(float)pack_objects.size();
		std::cout<<"percentage of potential growth: "<<perc*100.0f<<'%'<<std::endl;
		//stay_at_this_barrier = (perc > 0.005); // if only few tiles can be enlarged
		stay_at_this_barrier = (perc > 0);
		if (!enlarge)
		{
			return true;
		}
		else
			return stay_at_this_barrier;
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

	void Packer::vecfield_align()
	{
		for (unsigned int i = 0; i < pack_objects.size(); i++)
		{
			if (!pack_objects[i].active)	continue;
			const Facet& f = mesh[pack_objects[i].facet_idx];
			const vec3 v0 = mesh.vector_field_at(f.vertex_index[0]);
			const vec3 v1 = mesh.vector_field_at(f.vertex_index[1]);
			const vec3 v2 = mesh.vector_field_at(f.vertex_index[2]);
			vec3 v = (v0 + v1 + v2) / 3.0;
 			Local_frame lf = pack_objects[i].local_frame();
			Vector_3 av(v.x, v.y, v.z);
			Line_3 axis;
			Vector_3 principle_dir;
			double longest_edge_len = -std::numeric_limits<double>::max();
			for (unsigned int j = 0 ; j < pack_objects[i].size(); j++)
			{
				double sl = pack_objects[i].edge(j).squared_length();
				if ( sl > longest_edge_len)
				{
					longest_edge_len = sl;
					principle_dir = pack_objects[i].edge(j).to_vector();
				}
			}
			Vector_2 principle_dir_2 = lf.to_uv(principle_dir);
			Vector_2 av_2 = lf.to_uv(av);
			cgal_vec_normalize(principle_dir_2);
			cgal_vec_normalize(av_2);
			double det = principle_dir_2.x()*av_2.y() - principle_dir_2.y()*av_2.x();
			double cos_theta = principle_dir_2*av_2;
			if (cos_theta > 1.0)	cos_theta = 1.0;
			if (cos_theta < -1.0)	cos_theta = -1.0;
			double theta = std::acos(cos_theta);
			if (det < 0.0)	theta = -theta;
			if (theta > PI/3 && theta <= PI/2)	theta = theta - PI/2;
			else if (theta > PI/2 && theta <= 2*PI/3)	theta = theta - PI/2;
			else if (theta >= -PI/2 && theta < - PI/3)	theta = PI/2 + theta;
			else if (theta >= -2*PI/3 && theta < -PI/2)	theta = PI/2 + theta;
			transform_one_polygon(i, lf, Parameter(1.0, theta, 0.0, 0.0));
		}
		generate_RDT();
		compute_clipped_VD(false);
	}

	double Packer::vecfield_align_angle(unsigned int idx)
	{
		const Facet& f = mesh[pack_objects[idx].facet_idx];
		const vec3 v0 = mesh.vector_field_at(f.vertex_index[0]);
		const vec3 v1 = mesh.vector_field_at(f.vertex_index[1]);
		const vec3 v2 = mesh.vector_field_at(f.vertex_index[2]);
		vec3 v = (v0 + v1 + v2) / 3.0;
		Local_frame lf = pack_objects[idx].local_frame();
		Vector_3 av(v.x, v.y, v.z);
		Line_3 axis;
		CGAL::linear_least_squares_fitting_3(pack_objects[idx].vertices_begin(), pack_objects[idx].vertices_end(), axis, CGAL::Dimension_tag<0>());
		Vector_3 principle_dir = axis.to_vector();
		Vector_2 principle_dir_2 = lf.to_uv(principle_dir);
		Vector_2 av_2 = lf.to_uv(av);
		cgal_vec_normalize(principle_dir_2);
		cgal_vec_normalize(av_2);
		double det = principle_dir_2.x()*av_2.y() - principle_dir_2.y()*av_2.x();
		double cos_theta = principle_dir_2*av_2;
		if (cos_theta > 1.0)	cos_theta = 1.0;
		if (cos_theta < -1.0)	cos_theta = -1.0;
		double theta = std::acos(cos_theta);
		if (det < 0.0)	theta = -theta;
		return theta;
	}
	void Packer::lloyd(void (*post_action)(), bool enlarge)
	{
		static int times = 0;

		for (unsigned int i = 0; i < 1; i++)
		{
			std::vector<Parameter> solutions(pack_objects.size());
			std::vector<Local_frame> local_frames(pack_objects.size());

			std::cout<<"============ lloyd iteration "<<times++<<" ============\n";

			Lloyd_res res = one_lloyd(enlarge, solutions, local_frames);
			if (res == NO_MORE_ENLARGEMENT)
			{
			//	stop_update_DT = true;
				std::cout<<"Caution: Stopped updating iDT\n";
			}
			//if (!stop_update_DT)
			//{
				std::cout<<"Start iDT updating...\n";
				CGAL::Timer t;
				t.start();
				//if (vector_field)
				//	generate_RDT();
				//else
					rpvd.iDT_update();
				t.stop();
				std::cout<<"Time spent in iDT: "<<t.time()<<" seconds.\n";
				std::cout<<"End iDT updating\n";	
				compute_clipped_VD();
			//}			
		
			if (post_action != NULL)
				post_action();
		}

	}

	bool Packer::discrete_lloyd(void (*post_action)(), bool enlarge)
	{
		static int times = 0;
		bool stay_at_this_barrier;
		for (unsigned int i = 0; i < 1; i++)
		{
			double current_barrier;
			if ( !disc_barr.get_current_barrier(current_barrier))
			{
				std::cout<<"The maximum barrier factor achieved. Do nothing.\n";
				return false;
			}

			std::vector<Parameter> solutions(pack_objects.size());
			std::vector<Local_frame> local_frames(pack_objects.size());

			std::cout<<"============ discrete lloyd iteration "<<times++<<" ============\n";
			
			double nxt_barrier;
			//if (!disc_barr.get_next_barrier(nxt_barrier))
			//{
				//nxt_barrier = 1.0e20;
			//	std::cout<<"====== The largest barrier has been achieved ======\n";
			//	break;
			//}
			std::cout<<"current = "<<current_barrier<<std::endl;
			stay_at_this_barrier = discrete_one_lloyd(enlarge, solutions, local_frames, current_barrier, nxt_barrier);
			
			if (!stay_at_this_barrier)  // go to the next barrier
			{
				disc_barr.go_to_next_barrier();
				std::cout<<"go to next barrier"<<std::endl;
			}

			// check whether to use approximate voronoi region
			//bool use_appox_VD = false;
			double min_factor = std::numeric_limits<double>::max();
			for (unsigned int j = 0; j < solutions.size(); j++)
				if (solutions[j].k > 1.0)
					min_factor = std::min(solutions[j].k, min_factor);		
			use_voronoi_cell_ = !(enlarge && (min_factor < 1.05));
			rpvd.iDT_update();

			compute_clipped_VD(false);

			//	compute_clipped_VD(true);

			if (post_action != NULL)
				post_action();
		}
		return true;
	}
	
	void Packer::interface_lloyd(void (*post_action)())
	{
		//if (!sync_opt)
		//	lloyd(post_action, false);
		//else
			discrete_lloyd(post_action, false);
	}
	bool Packer::pack(void (*post_action)())
	{	
		//if (!sync_opt)
		//	lloyd(post_action, true);
		//else
		//{
			//if (sync_opt && disc_barr.get_max() != phony_upper_scale)
			//	disc_barr.append(phony_upper_scale);
			bool not_finished = discrete_lloyd(post_action, true);
		//}
		print_area_coverage();	
		return not_finished;
	}

	void Packer::constraint_transformation(vector<Parameter>& parameters, vector<Local_frame>& lfs, bool constrain_scale)
	{
		std::vector<double> min_factors(parameters.size(), 1.0);
		typedef RestrictedPolygonVoronoiDiagram RPVD;
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
				Parameter sp = parameters[i]; // constrained transformation to try
				if (constrain_scale)
					sp.k = (sp.k - 1.0) * min_factors[i] + 1.0;
				else
					sp = parameters[i]*min_factors[i];
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
					//samp_pnts[j]->mp = future_samp_pnts[j];
				}
			}
			stop = true;
			std::vector<bool> flipped(parameters.size(), false);
			int nb_flipped = 0;
			for (RPVD::Facet_iterator fit = rpvd.faces_begin(); fit != rpvd.faces_end(); ++fit)
			{
				RPVD::Halfedge_handle eh = fit->halfedge();
				Point_3 v0 = eh->vertex()->mp;
				//Point_3 dv0 = eh->vertex()->point();
				int i0 = eh->vertex()->group_id;

				eh = eh->next();
				Point_3 v1 = eh->vertex()->mp;
				//Point_3 dv1 = eh->vertex()->point();
				int i1 = eh->vertex()->group_id;

				eh = eh->next();
				Point_3 v2 = eh->vertex()->mp;
				//Point_3 dv2 = eh->vertex()->point();
				int i2 = eh->vertex()->group_id;

				if (i0 == i1 && i0 == i2)
					continue;

				Vector_3 v01(v0, v1), v12(v1, v2);
				cgal_vec_normalize(v01);
				cgal_vec_normalize(v12);
				Vector_3 n = CGAL::cross_product(v01, v12);
				double nlen2 = n.squared_length();
 				if (nlen2 <= 1.0e-8) //degenerate case
 				{
 					nb_flipped++;
 					if (i0 >= 0)	flipped[i0] = true;
 					if (i1 >= 0)	flipped[i1] = true;
 					if (i2 >= 0)	flipped[i2] = true;
 					stop = false;
 					continue;
 				}
				cgal_vec_normalize(n);
				if (n*fit->n <= 1.0e-6)
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

		if (constrain_scale)
			for (unsigned int i = 0; i < parameters.size(); i++)
				parameters[i].k = (parameters[i].k - 1.0) * min_factors[i] + 1.0;
		else
			for (unsigned int i = 0; i < parameters.size(); i++)
				parameters[i] *= min_factors[i];
		std::cout<<"Trial times for constraint: "<<n<<std::endl;
			
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

		std::vector<Polygon_2> polygon_2d(pack_objects.size());
		std::vector<Local_frame> lf;
		lf.reserve(pack_objects.size());
		for (unsigned int i = 0; i < pack_objects.size(); i++)
		{
			Local_frame loc_frame = pack_objects[i].local_frame();
			lf.push_back(loc_frame);
			for (unsigned int j = 0; j < pack_objects[i].size(); j++)
				polygon_2d[i].push_back(loc_frame.to_uv(pack_objects[i].vertex(j)));
		}

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
			{
				if (v0->group_id < 0)
					fit->vacant = false;
				else
				{
					Point_3 c = CGAL::ORIGIN + (Vector_3(CGAL::ORIGIN, v0->point())+Vector_3(CGAL::ORIGIN, v1->point())+Vector_3(CGAL::ORIGIN, v2->point()))/3.0;
					int idx = v0->group_id;
					Point_2 c_2d = lf[idx].to_uv(c);
					if (polygon_2d[idx].has_on_unbounded_side(c_2d))
						fit->vacant = true;
					else
						fit->vacant = false;
				}
			}
			else
			{	
				fit->vacant = true;
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
		std::map< int, Match_info_item<unsigned int> > replace_map;
		for (size_t i = 0; i < holes.size(); i++)
		{
			// collect the tiles surrounding a hole
			std::set<int> nghb_tiles;
			for (size_t j = 0; j < holes[i].size(); j++)
			{
				nghb_tiles.insert(holes[i][j].first->group_id);
				nghb_tiles.insert(holes[i][j].second->group_id);
			}
			std::map<Match_info_item<unsigned int>, int, Match_measure> match_idx;
			for (std::set<int>::const_iterator it = nghb_tiles.begin(); it != nghb_tiles.end(); ++it)
			{
				if (*it < 0)	continue;
				// extract the surrounding space
				Hole hl;
				vacate_one_polygon(*it, hl);
				std::vector<Segment_2> hole_bd_2d;
				std::vector<Point_3> hole_verts;
				hole_verts.reserve(hl.size()*2);
				for (unsigned int j = 0; j < hl.size(); j++)
				{
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
					Point_2 s = lf.to_uv( prj_plane.projection(hl[j].first->point()) );
					Point_2 t = lf.to_uv( prj_plane.projection(hl[j].second->point()) );
					hole_bd_2d.push_back( Segment_2(s, t) );
				}
				Polygon_matcher pm(hole_bd_2d, 100);
				std::priority_queue<Match_info_item<unsigned int>, std::vector< Match_info_item<unsigned int> >, Match_measure> match_res;
				for (unsigned int idx = 0; idx < pgn_lib.size(); idx++)
					match_res.push(pm.affine_match(pgn_lib[idx], idx));
				// choose the result with the smallest match error now
				Match_info_item<unsigned int> matcher = match_res.top();
				match_idx[matcher] = *it;
			}
			// the best matching one
			std::map<Match_info_item<unsigned int>, int, Match_measure>::reverse_iterator best_match = match_idx.rbegin();
			replace_map[best_match->second] = best_match->first;
		}
		for (std::map< int, Match_info_item<unsigned int> >::iterator it = replace_map.begin(); it != replace_map.end(); ++it)
		{
			Packing_object& filler = pack_objects[it->first];
			Hole hl;
			vacate_one_polygon(it->first, hl);
			std::vector<Point_3> hole_verts;
			hole_verts.reserve(hl.size()*2);
			for (unsigned int j = 0; j < hl.size(); j++)
			{
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
			filler.clear();
			Match_info_item<unsigned int>& matcher = it->second;
			const Ex_polygon_2& match_pgn = pgn_lib[matcher.val];
			Polygon_2 transformed_pgn2d = CGAL::transform(matcher.t, match_pgn);
			for (unsigned int j = 0; j < transformed_pgn2d.size(); j++)
			{
				Point_3 p = lf.to_xy(transformed_pgn2d.vertex(j));
				filler.push_back(p);
			}
			filler.align(lf.w, lf.o);
			//std::transform(filler.vertices_begin(), filler.vertices_end(), filler.vertices_begin(), rescalor);
			filler.factor = 1.0;
			filler.facet_idx = fid;
			filler.lib_idx = matcher.val;
			filler.texture_coord.assign(match_pgn.texture_coords.begin(), match_pgn.texture_coords.end());
			filler.texture_id = match_pgn.texture_id;
		}
		// rebuild the restricted delaunay triangulation and voronoi cell
		generate_RDT();
		compute_clipped_VD();
		stop_update_DT = false;
		use_voronoi_cell_ = false;
		std::for_each(holes.begin(), holes.end(), std::mem_fun_ref(&Hole::clear));
		holes.clear();
		std::for_each(pack_objects.begin(), pack_objects.end(), std::mem_fun_ref(&Packing_object::activate));
		if (sync_opt)
		{
			double min_size = std::numeric_limits<double>::max();
			for (unsigned int i = 0; i < pack_objects.size(); i++)
				min_size = std::min(min_size, pack_objects[i].factor);
			disc_barr.set_current_barrier(min_size);
		}
		//std::cout<<"End replacing. Computation time: "<< replace_timer.time()<<" seconds.\n";
	}

	void Packer::curv_constrained_transform(Parameter& para, int fid, unsigned int pgn_id)
	{
		// approximate curvature at the facet
		const Facet& f = mesh[fid];
		double v0_cur = mesh.curvature_at_vertex(f.vertex_index[0]),
				v1_cur = mesh.curvature_at_vertex(f.vertex_index[1]),
				v2_cur = mesh.curvature_at_vertex(f.vertex_index[2]);
		double avg_cur = (v0_cur + v1_cur + v2_cur) / 3.0 /*std::max(v0_cur, v1_cur)*/;
		if (avg_cur == 0.0)
			return;
		double r = 1.0 / avg_cur;
		double squared_translation = para.tx * para.tx + para.ty * para.ty;
		double threshold = (2.0 - epsilon) * epsilon;
		if (squared_translation > threshold * r * r)
		{
			double temp = std::sqrt(threshold)*r/std::sqrt(squared_translation);
			para.tx *= temp;
			para.ty *= temp;
		}
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
			double temp = std::sqrt(threshold)*r/(para.k*mean_radius);
			para.k = std::max(1.0, para.k*temp);
		}
	}

	void Packer::fill_holes()
	{
		std::cout<<"Start filling holes...\n";
		for (unsigned int i = 0; i < pack_objects.size(); i++)
			pack_objects[i].active = false;
		size_t nb_filled = 0;
		for (unsigned int i = 0; i < holes.size(); i++)
		{
			pack_objects.push_back(Packing_object());
			std::set<int> hole_neighbors;
			for (size_t j = 0; j < holes[i].size(); j++)
			{
				if (holes[i][j].first->group_id >= 0)
					rpvd.get_neighbor_indices(holes[i][j].first->group_id, hole_neighbors);
				if (holes[i][j].second->group_id >= 0)
					rpvd.get_neighbor_indices(holes[i][j].second->group_id, hole_neighbors);
			}
			bool success = fill_one_hole(holes[i], pack_objects.back());
			if (success)
			{
				pack_objects.back().activate();
				nb_filled++;
				for (std::set<int>::const_iterator it = hole_neighbors.begin(); it != hole_neighbors.end(); ++it)
				{
					//std::cout<<*it<<' ';
					if (*it >= 0)
						pack_objects[*it].activate();
				}
				//std::cout<<std::endl;
			}
			else
				pack_objects.pop_back();
		}
		generate_RDT();
		stop_update_DT = false;
		use_voronoi_cell_ = false;
		compute_clipped_VD();
		for (size_t i = 0; i < pack_objects.size(); i++)
			if (pack_objects[i].active)
			{
				Local_frame lf = pack_objects[i].local_frame();
				Parameter scale(replace_factor, 0.0, 0.0, 0.0);
				transform_one_polygon(i, lf, scale);
			}
		std::cout<<"End filling holes: "<<nb_filled<<" holes were filled.\n";
		std::for_each(holes.begin(), holes.end(), std::mem_fun_ref(&Hole::clear));
		holes.clear();
		//std::for_each(pack_objects.begin(), pack_objects.end(), std::mem_fun_ref(&Packing_object::activate));
		if (sync_opt)
		{
			double min_size = std::numeric_limits<double>::max();
			for (unsigned int i = 0; i < pack_objects.size(); i++)
			{				
				//min_size = std::min(pack_objects[i].rel_factor(mesh.curvature_at_face(pack_objects[i].facet_idx)), min_size);
				min_size = std::min(min_size, pack_objects[i].factor);
			}
			disc_barr.set_current_barrier(min_size);
			//disc_barr.set(0.9, 1.0, 2);
		}
	}

	bool Packer::fill_one_hole(Hole& hl, Packing_object& filler)
	{
		std::vector<Segment_2> hole_bd_2d;
		std::vector<Point_3> hole_verts;
		hole_verts.reserve(hl.size()*2);
		for (unsigned int j = 0; j < hl.size(); j++)
		{
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
			Point_2 s = lf.to_uv( prj_plane.projection(hl[j].first->point()) );
			Point_2 t = lf.to_uv( prj_plane.projection(hl[j].second->point()) );
			hole_bd_2d.push_back( Segment_2(s, t) );
		}
		Polygon_matcher pm(hole_bd_2d, 200);
		std::priority_queue<Match_info_item<unsigned int>, std::vector< Match_info_item<unsigned int> >, Match_measure> match_res;
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
		double shrink_factor = 1.0;
		// if the tile is near the boundary, contract it more
		for (size_t j = 0; j < hl.size(); j++)
			if (hl[j].first->group_id < 0 || hl[j].second->group_id < 0)
			{
				shrink_factor = 0.6;
				break;
			}
		Transformation_3 rescalor = Transformation_3(CGAL::TRANSLATION, Vector_3(CGAL::ORIGIN, lf.o)) *
			( Transformation_3(CGAL::SCALING, shrink_factor) *
			Transformation_3(CGAL::TRANSLATION, Vector_3(lf.o, CGAL::ORIGIN)) );	

		std::transform(filler.vertices_begin(), filler.vertices_end(), filler.vertices_begin(), rescalor);
		filler.factor = matcher.scale * shrink_factor;
		filler.facet_idx = fid;
		filler.lib_idx = matcher.val;
		filler.texture_coord.assign(match_pgn.texture_coords.begin(), match_pgn.texture_coords.end());
		filler.texture_id = match_pgn.texture_id;	
		return true;
	}

	void Packer::con_replace()
	{
		std::cout<<"Start replacing...\n";
		CGAL::Timer replace_timer;
		replace_timer.start();
		
		/* for rigid size version, curvature sensitive initialization */
// 		std::vector<double> vtx_cur(mesh.nb_vertices());
// 		for (size_t i = 0; i < mesh.nb_vertices(); i++)
// 			vtx_cur[i] = -mesh.curvature_at_vertex(i);
// 		std::sort(vtx_cur.begin(), vtx_cur.end());
// 		std::map<double, size_t> area_idx;
// 		for (size_t i = 0; i < pgn_lib.size(); i++)
// 			area_idx[std::fabs(pgn_lib[i].area())] = i;
		/**************************************************************/

#ifdef _CILK_
		cilk_for (unsigned int i = 0; i < pack_objects.size(); i++)
#else
		for (unsigned int i = 0; i < pack_objects.size(); i++)
#endif
		{
			const RestrictedPolygonVoronoiDiagram::VertGroup& samp_pnts = rpvd.sample_points_group(i);
			// check whether near boundary
			bool near_boundary = false;
			for (unsigned int j = 0; j < samp_pnts.size(); j++)
			{
				RestrictedPolygonVoronoiDiagram::Vertex_const_handle v = samp_pnts[j]->halfedge()->opposite()->vertex();
				if (v->group_id < 0)
				{
					near_boundary = true;
					break;
				}
			}

			const std::vector<Point_3>& smoothed_region = rpvd.get_smoothed_voronoi_regions().at(i);

			Local_frame lf = pack_objects[i].local_frame();

			std::vector<Segment_2> region2d;
			
			for (unsigned int j = 0; j < smoothed_region.size(); j++)
			{
				Point_2 s = lf.to_uv(smoothed_region[j]);
				Point_2 t = lf.to_uv(smoothed_region[(j+1)%smoothed_region.size()]);
				region2d.push_back(Segment_2(s, t));
			}
			Polygon_matcher pm(region2d, 200);
			std::priority_queue<Match_info_item<unsigned int>, std::vector< Match_info_item<unsigned int> >, Match_measure> match_res;
			for (unsigned int idx = 0; idx < pgn_lib.size(); idx++)
			{
				Match_info_item<unsigned int> res = pm.affine_match(pgn_lib[idx], idx, match_weight);
				//res.error = (1-match_weight)*res.error + match_weight*sim_mat[pack_objects[i].lib_idx][idx];
				match_res.push(res);
			}

			// choose the result with the smallest match error now 

			double shrink_factor = replace_factor;

			/* for rigid size version, curvature sensitive initialization */
// 			double fcur = mesh.curvature_at_face(pack_objects[i].facet_idx);
// 			std::vector<double>::iterator cur_it = std::lower_bound(vtx_cur.begin(), vtx_cur.end(), -fcur);
// 			size_t cur_rank = cur_it - vtx_cur.begin();
// 			std::map<double, size_t>::iterator it = area_idx.begin();
// 			std::advance(it, float(cur_rank)/vtx_cur.size()*pgn_lib.size());
			/************************************************************************/


			Match_info_item<unsigned int> matcher = match_res.top();
			// recompute the right transformation
			matcher = pm.affine_match(pgn_lib[matcher.val], matcher.val, match_weight);
			if (sync_opt)
			{
				if ( matcher.scale * shrink_factor > disc_barr.get_max() )
					shrink_factor = disc_barr.get_max() / matcher.scale;
			}
			if (near_boundary)	
				shrink_factor *= 0.8;

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

			Transformation_3 rescalor = Transformation_3(CGAL::TRANSLATION, Vector_3(CGAL::ORIGIN, to_cgal_pnt(prjp))) *
				( Transformation_3(CGAL::SCALING, shrink_factor) *
				Transformation_3(CGAL::TRANSLATION, Vector_3(to_cgal_pnt(prjp), CGAL::ORIGIN)) );	

			std::transform(pack_objects[i].vertices_begin(), pack_objects[i].vertices_end(), pack_objects[i].vertices_begin(), rescalor);
			
			pack_objects[i].lib_idx = matcher.val;
			pack_objects[i].factor = matcher.scale*shrink_factor;
			pack_objects[i].facet_idx = fid;
			pack_objects[i].texture_coord.assign(match_pgn.texture_coords.begin(), match_pgn.texture_coords.end());
			pack_objects[i].texture_id = match_pgn.texture_id;
		}
		replace_timer.stop();
		std::cout<<"End replacing. Computation time: "<< replace_timer.time()<<" seconds.\n";
		// rebuild the restricted delaunay triangulation and voronoi cell
		generate_RDT();
		compute_clipped_VD();
		stop_update_DT = false;
		use_voronoi_cell_ = true;
		std::for_each(pack_objects.begin(), pack_objects.end(), std::mem_fun_ref(&Packing_object::activate));
		if (sync_opt)
		{
			double min_size = std::numeric_limits<double>::max();
			for (unsigned int i = 0; i < pack_objects.size(); i++)
			{				
				min_size = std::min(pack_objects[i].rel_factor(mesh.curvature_at_face(pack_objects[i].facet_idx)), min_size);
			}
			disc_barr.set_current_barrier(min_size);
		}
		
	}
	void Packer::swap()
	{
		CGAL::Timer swap_timer;
		swap_timer.start();
// 		std::vector<unsigned int> group_indices(pack_objects.size());
// 		for (unsigned int i = 0; i < pack_objects.size(); i++)
// 			group_indices[i] = i;
		// grouping
		std::vector<unsigned int> group_indices;
		// build the adjacency matrix
		std::vector< std::vector<bool> > adj_mat(pack_objects.size());
		for (unsigned int i = 0; i < pack_objects.size(); i++)
		{
			adj_mat[i].assign(pack_objects.size(), false);
			adj_mat[i][i] = true;
		}
		for (RestrictedPolygonVoronoiDiagram::Halfedge_const_iterator eit = rpvd.halfedges_begin(); eit != rpvd.halfedges_end(); ++eit)
		{
			RestrictedPolygonVoronoiDiagram::Vertex_const_handle v0 = eit->vertex(), v1 = eit->opposite()->vertex();
			if (v0->group_id >= 0 && v1->group_id >= 0)
			{
				if (v0->group_id != v1->group_id)
				{
					adj_mat[v0->group_id][v1->group_id] = true;
					adj_mat[v1->group_id][v0->group_id] = true;
				}
			}
		}
		// form the group
		for (unsigned int i = 0; i < pack_objects.size(); i++)
		{
			if (swapped_indices.find(i) != swapped_indices.end())	continue;
			bool adj_to_any = false;
			for (unsigned int j = 0; j < group_indices.size(); j++)
				if (adj_mat[group_indices[j]][i])
				{
					adj_to_any = true;
					break;
				}
			if (!adj_to_any)
			{
				group_indices.push_back(i);
				swapped_indices.insert(i);
			}
		}
		unsigned int nb_swapped = group_swap(group_indices);
		swap_timer.stop();
		std::cout<<"Swap time: "<<swap_timer.time()<<" seconds. Number of swapped tiles: "<<nb_swapped<<std::endl;
		generate_RDT();
		compute_clipped_VD();
		stop_update_DT = false;
		use_voronoi_cell_ = true;
		std::for_each(pack_objects.begin(), pack_objects.end(), std::mem_fun_ref(&Packing_object::activate));
		if (sync_opt)
		{
			double min_size = std::numeric_limits<double>::max();
			for (unsigned int i = 0; i < pack_objects.size(); i++)			
				min_size = std::min(pack_objects[i].rel_factor(mesh.curvature_at_face(pack_objects[i].facet_idx)), min_size);
			disc_barr.set_current_barrier(min_size);
		}	
		if (swapped_indices.size() == pack_objects.size()) // all swapped
		{
			std::cout<< "All swapped.\n";
			swapped_indices.clear();
		}
		swapped_this_time.assign(group_indices.begin(), group_indices.end());
	}
	void Packer::selective_swap()
	{
		// select the tiles with large surrounding area
		std::map< double, unsigned int, std::greater<double> > vacancy_index;
		for (unsigned int i = 0; i < pack_objects.size(); i++)
		{
			Hole h;
			vacate_one_polygon(i, h);
			// the area of the hole
			Point_3 c = pack_objects[i].centroid();
			double hole_area = 0.0;
			for (unsigned int j = 0; j < h.size(); j++)
			{
				Vector_3 area_vec = CGAL::cross_product(Vector_3(c, h[j].first->point()), Vector_3(c, h[j].second->point()));
				hole_area += CGAL::sqrt(area_vec.squared_length());
			}
			hole_area /= 2;
			// the area of neighboring space
			double vacant_area = hole_area - pack_objects[i].area();
			if (vacant_area < 0.0)	vacant_area = 0.0;
			vacancy_index[vacant_area] = i;
			pack_objects[i].active = false;
		}
		const unsigned int N = match_weight*pack_objects.size();
		std::map< double, unsigned int, std::greater<double> >::const_iterator it = vacancy_index.begin();
		std::vector<unsigned int> selected_tiles;
		selected_tiles.reserve(N);
		for (unsigned int j = 0; j < N; j++, ++it)
		{
			selected_tiles.push_back(it->second);
		}
		group_swap(selected_tiles);
		for (unsigned int i = 0; i < selected_tiles.size(); i++)
			pack_objects[selected_tiles[i]].active = true;
	}
	unsigned int Packer::group_swap(const std::vector<unsigned int>& indices)
	{
		if (indices.size() <= 1)	return 0;
		vicinity_holes.clear();
		mat m(indices.size()); 
		//std::vector< std::vector<Segment_2> > region2d(indices.size());
		std::vector< HoleEdges > region2d(indices.size());
		for (unsigned int i = 0; i < indices.size(); i++)
		{
			m[i].resize(indices.size());	
			
			Local_frame lf = pack_objects[indices[i]].local_frame();
			/////////////////////		Vicinity		////////////////////////////
			Hole h;
			vacate_one_polygon(indices[i], h);
 			region2d[i].reserve(h.size());
 			vicinity_holes.push_back(std::list<Segment_3>());
 			for (unsigned int j = 0; j < h.size(); j++)
 			{
 				Point_2 s = lf.to_uv(h[j].first->point());
 				Point_2 t = lf.to_uv(h[j].second->point());
 				//region2d[i].push_back(Segment_2(s, t));
				region2d[i].push_back(std::make_pair(vec2(s.x(), s.y()), vec2(t.x(), t.y())));
 				vicinity_holes.back().push_back(Segment_3(h[j].first->point(),h[j].second->point()));
 			}
			///////////////////////////////////////////////////////////////////////
			//////////////////		Hausdorff	///////////////////////
//  			std::vector<Vertex_handle> loop;
//  			sortHoleEdge(h, loop);
//  			for (unsigned int j = 0; j < loop.size(); j++)
//  			{
//  				region2d[i].push_back(lf.to_uv(loop[j]->point()));
//  			}
 			std::vector<double> hausdorff_errors(pgn_lib.size());
 			for (unsigned int j = 0; j < pgn_lib.size(); j++)
 			{
 				double tx, ty, k, theta;
 				hausdorff_errors[j] = matching_score(region2d[i], pgn_lib[j], tx, ty, k, theta);
 			}
 			for (unsigned int j = 0; j < indices.size(); j++)
 				m[i][j] = -hausdorff_errors[pack_objects[indices[j]].lib_idx];
			///////////////////////////////////////////////////////////
			//////////////////		Affine match		///////////////
// 			Polygon_matcher pm(region2d[i], 200);
// 			std::vector< Match_info_item<unsigned int> > match_res;
// 			match_res.reserve(pgn_lib.size());
// 			for (unsigned int j = 0; j < pgn_lib.size(); j++)
// 			{
// 				Match_info_item<unsigned int> res = pm.affine_match(pgn_lib[j], j, match_weight);
// 				match_res.push_back(res);
// 			}
// 			for (unsigned int j = 0; j < indices.size(); j++)
// 				m[i][j] = -match_res[pack_objects[indices[j]].lib_idx].error;
			///////////////////////////////////////////////////////////
		}
		iAuction auction(m);
		auction.MainAlgo();
		// begin swap
		unsigned int nb_swapped = 0;
		std::vector<unsigned int> pre_indices(indices.size());
		for (unsigned int i = 0; i < indices.size(); i++)
			pre_indices[i] = pack_objects[indices[i]].lib_idx;
		for (unsigned int i = 0; i < indices.size(); i++)
		{
			unsigned int pid = indices[i];
			const RestrictedPolygonVoronoiDiagram::VertGroup& samp_pnts = rpvd.sample_points_group(pid);
			bool near_boundary = false;
			for (unsigned int j = 0; j < samp_pnts.size(); j++)
			{
				RestrictedPolygonVoronoiDiagram::Vertex_const_handle v = samp_pnts[j]->halfedge()->opposite()->vertex();
				if (v->group_id < 0)
				{
					near_boundary = true;
					break;
				}
			}
			unsigned int idx = pre_indices[auction.GetAssignedCol(i)];
			if (idx != pack_objects[pid].lib_idx)		nb_swapped++;
			// choose the result with the smallest match error now 
			double shrink_factor = replace_factor;
//  		Polygon_matcher pm(region2d[i], 200);
//  		Match_info_item<unsigned int> matcher = pm.affine_match(pgn_lib[idx], idx, match_weight);
 			double tx, ty, k, theta;
 			matching_score(region2d[i], pgn_lib[idx], tx, ty, k, theta);
 			Match_info_item<unsigned int> matcher;
 			matcher.scale = k;
 			matcher.val = idx;
			if (sync_opt)
			{
				if ( matcher.scale * shrink_factor > disc_barr.get_max() )
					shrink_factor = disc_barr.get_max() / matcher.scale;
			}
			if (near_boundary)	
				shrink_factor *= 0.8;

			const Ex_polygon_2& match_pgn = pgn_lib[idx];
			//Polygon_2 transformed_pgn2d = CGAL::transform(matcher.t, match_pgn);
 			Polygon_2 transformed_pgn2d = match_pgn;
 			transform_poly(transformed_pgn2d, tx, ty, k, theta);

			Local_frame lf = pack_objects[pid].local_frame();

			pack_objects[pid].clear();

			for (unsigned int j = 0; j < transformed_pgn2d.size(); j++)
			{
				Point_3 p = lf.to_xy(transformed_pgn2d.vertex(j));
				pack_objects[pid].push_back(p);
			}

			Point_3 c = pack_objects[pid].centroid();
			vec3 v;
			int fid;
			vec3 prjp = mesh.project_to_mesh(to_geex_pnt(c), v, fid);

			v = approx_normal(fid);
			pack_objects[pid].align(to_cgal_vec(v), to_cgal_pnt(prjp));

			Transformation_3 rescalor = Transformation_3(CGAL::TRANSLATION, Vector_3(CGAL::ORIGIN, to_cgal_pnt(prjp))) *
				( Transformation_3(CGAL::SCALING, shrink_factor) *
				Transformation_3(CGAL::TRANSLATION, Vector_3(to_cgal_pnt(prjp), CGAL::ORIGIN)) );	

			std::transform(pack_objects[pid].vertices_begin(), pack_objects[pid].vertices_end(), pack_objects[pid].vertices_begin(), rescalor);

			pack_objects[pid].lib_idx = matcher.val;
			pack_objects[pid].factor = matcher.scale*shrink_factor;
			pack_objects[pid].facet_idx = fid;
			pack_objects[pid].texture_coord.assign(match_pgn.texture_coords.begin(), match_pgn.texture_coords.end());
			pack_objects[pid].texture_id = match_pgn.texture_id;
		}
		return nb_swapped;
	}
	void Packer::vacate_one_polygon(unsigned int id, Hole& hole)
	{
		const RestrictedPolygonVoronoiDiagram::VertGroup& samp_pnts = rpvd.sample_points_group(id);
		typedef RDT_data_structure::Halfedge_around_vertex_circulator Edge_circulator;
		hole.clear();
		for (unsigned int i = 0; i < samp_pnts.size(); i++)
		{
			Edge_circulator start_edge = samp_pnts[i]->vertex_begin();
			Edge_circulator current_edge = start_edge;
			do 
			{
				Vertex_handle v_adj = current_edge->opposite()->vertex();
				if (v_adj->group_id == samp_pnts[i]->group_id)
					break;
				++current_edge;
			} while (current_edge != start_edge);
			start_edge = current_edge;
			do 
			{
				if (current_edge->facet() != Facet_handle() )
				{
					Vertex_handle v0 = current_edge->prev()->vertex(), v1 = current_edge->next()->vertex();
					if (v0->group_id != samp_pnts[i]->group_id && v1->group_id != samp_pnts[i]->group_id)
						hole.push_back(std::make_pair(v0, v1));
				}
				++current_edge;
			} while (current_edge != start_edge);
		}
	}
	void Packer::report()
	{
		double max_factor = -std::numeric_limits<double>::max();
		double min_factor = std::numeric_limits<double>::max();
		double cur_max_factor = -std::numeric_limits<double>::max();
		double cur_min_factor = std::numeric_limits<double>::max();
		double a_min, c_min, a_max, c_max;
		unsigned int min_idx, max_idx;
		double mean_factor = 0.0;
		double cur_mean_factor = 0.0;
		//std::set<unsigned int> lib_indices;
		std::map<unsigned int, unsigned int> frequency;
		double sum_curvature = 0.0;
		for (unsigned int i = 0; i < pack_objects.size(); i++)
		{
			double temp = mesh.curvature_at_face(pack_objects[i].facet_idx);
			sum_curvature += 1.0/(temp*temp+1.0);
		}

		for (unsigned int i = 0; i < pack_objects.size(); i++)
		{
			double res_area = pack_objects[i].area();
			//double normalized_factor = /*std::pow(pack_objects[i].factor, 0.75) * */pow(avg_cur, pio.gamma_used())*res_area;
			double cur = mesh.curvature_at_face(pack_objects[i].facet_idx);
			//double normalized_factor = std::sqrt(res_area / (area_coverage*mesh_area/(sum_curvature*(cur*cur+1.0))));
			//double normalized_factor = (avg_cur*avg_cur+1.0)*pack_objects[i].factor*pack_objects[i].factor;
			double normalized_factor = pack_objects[i].rel_factor(cur);
			max_factor = std::max(max_factor, pack_objects[i].factor);
			min_factor = std::min(min_factor, pack_objects[i].factor);
			//cur_max_factor = std::max(cur_max_factor, normalized_factor);
			//cur_min_factor = std::min(cur_min_factor, normalized_factor);
			if (cur_max_factor < normalized_factor)
			{
				cur_max_factor = normalized_factor;
				c_max = cur;
				a_max = res_area;
				max_idx = i;
			}
			if (cur_min_factor > normalized_factor)
			{
				cur_min_factor = normalized_factor;
				c_min = cur;
				a_min = res_area;
				min_idx = i;
			}
			mean_factor += pack_objects[i].factor;
			cur_mean_factor += normalized_factor;
			//lib_indices.insert(pack_objects[i].lib_idx);
			frequency[pack_objects[i].lib_idx]++;
		}
		mean_factor /= pack_objects.size();
		cur_mean_factor /= pack_objects.size();
		std::cout<<"Absolute size: \n";
		std::cout<<"\t-- Maximum factor: "<<max_factor<<std::endl;
		std::cout<<"\t-- Minimum factor: "<<min_factor<<std::endl;
		std::cout<<"\t-- Mean factor: "<<mean_factor<<std::endl;
		//std::cout<<"Curvature relative size: \n";
		//std::cout<<"\t-- Maximum factor: "<<cur_max_factor<<" : "<<a_max<<", "<<c_max<<", "<<max_idx<<std::endl;
		//std::cout<<"\t-- Minimum factor: "<<cur_min_factor<<" : "<<a_min<<", "<<c_min<<", "<<min_idx<<std::endl;
		//std::cout<<"\t-- Mean factor: "<<cur_mean_factor<<std::endl;
		//std::cout<<"Polygon variety:\n";
		//std::cout<<"\t-- Number of polygon types from input: "<<pgn_lib.size()<<std::endl;
		//std::cout<<"\t-- Number of polygon types: "<<frequency.size()<<std::endl;
		//std::ofstream freq_file("usage-frequency.txt");
		//for (std::map<unsigned int, unsigned int>::const_iterator it = frequency.begin(); it != frequency.end(); ++it)
		//	freq_file<<it->first<<' '<<it->second<<std::endl;
	}

	void Packer::print_area_coverage()
	{
		// compute area coverage
#ifdef _CILK_
		cilk::reducer_opadd<double> sum_pgn_area(0.0);
		cilk_for (unsigned int i = 0; i < pack_objects.size(); i++)
			sum_pgn_area += pack_objects[i].area();
		area_coverage = sum_pgn_area.get_value()/mesh_area;
		std::cout<<"-- Area coverage ratio: "<<area_coverage<<std::endl;
#else
		double sum_pgn_area = 0.0;
		for (unsigned int i = 0; i < pack_objects.size(); i++)
			sum_pgn_area += pack_objects[i].area();
		area_coverage = sum_pgn_area/mesh_area;
		std::cout<<"-- Area coverage ratio: "<<area_coverage<<std::endl;
#endif
	}

	void Packer::save_tiles()
	{
		std::string obj_fn = pio.attribute_value("OutputDir") + "/tiles.obj";
		std::ofstream ofile(obj_fn.c_str());
		unsigned int vtx_idx = 1;
		std::vector<std::vector<unsigned int>> global_vtx_idx;
		global_vtx_idx.reserve(pack_objects.size());
		if (!pio.texture_specified())
		{
			// pseudo texture, this can be removed safely
			ofile << "mtllib "<<"mat.mtl\n";
			/////////////////////////////////////////////////////
			for (unsigned int i = 0; i < pack_objects.size(); i++)
			{
				global_vtx_idx.push_back(std::vector<unsigned int>());
				global_vtx_idx.back().reserve(pack_objects[i].size());
				for (unsigned int j = 0; j < pack_objects[i].size(); j++)
				{
					Point_3 v = pack_objects[i].vertex(j);
					ofile<<"v "<<v.x()<<' '<<v.y()<<' '<<v.z()<<std::endl;
					global_vtx_idx.back().push_back(vtx_idx);
					vtx_idx++;
				}
			}
			// pseudo texture, this can be removed safely
			for (unsigned int i = 0; i < pack_objects.size(); i++)
			{
				double xmin = std::numeric_limits<double>::max(), ymin = std::numeric_limits<double>::max();
				double xmax = std::numeric_limits<double>::min(), ymax = std::numeric_limits<double>::min();
				Local_frame lf = pack_objects[i].local_frame();
				Polygon_2 pgn2d;
				for (unsigned int j = 0; j < pack_objects[i].size(); j++)
					pgn2d.push_back(lf.to_uv(pack_objects[i].vertex(j)));
				for (unsigned int j = 0; j < pgn2d.size(); j++)
				{
					Point_2 v = pgn2d.vertex(j);
					xmin = std::min(v.x(), xmin);
					xmax = std::max(v.x(), xmax);
					ymin = std::min(v.y(), ymin);
					ymax = std::max(v.y(), ymax);
				}
				double w = std::max(xmax - xmin, ymax - ymin);
				double bd_area = w * w;
				//double pgn_area = std::fabs(pgn2d.area());
				int nx = 10, ny = 10;
				int idx = pack_objects[i].lib_idx%(nx*ny);
				int c = idx/nx, r = idx%ny;
				double f = std::sqrt(1.0/(nx*ny)/bd_area);
				double lbc_x = c*(1.0/nx), lbc_y = r*(1.0/ny);
				Point_2 lbc(lbc_x, lbc_y);
				Transformation_2 to_org(CGAL::TRANSLATION, Vector_2(Point_2(xmin, ymin), /*CGAL::ORIGIN*/lbc));
				Transformation_2 in_bd(CGAL::SCALING, 0.99*f);
				Transformation_2 t = in_bd*to_org;
				Polygon_2 tex_pgn = CGAL::transform(t, pgn2d);
				for (unsigned int j = 0; j < tex_pgn.size(); j++)
				{
					Point_2 v = tex_pgn.vertex(j);
					double tx = std::min(v.x(), 1.0);
					tx = std::max(v.x(), 0.0);
					double ty = std::min(v.y(), 1.0);
					ty = std::max(v.y(), 0.0);
					ofile<<"vt "<<tx<<' '<<ty<<std::endl;
				}
			}
			for (unsigned int i = 0; i < pack_objects.size(); i++)
			{
				Vector_3 n = pack_objects[i].norm();
				for (unsigned int j = 0; j < pack_objects[i].size(); j++)
					ofile<<"vn "<<n.x()<<' '<<n.y()<<' '<<n.z()<<std::endl;
			}
			/////////////////////////////////////////////////////
			for (unsigned int i = 0; i < pack_objects.size(); i++)
			{
				ofile<<"usemtl marble"<<std::endl;
				//for (unsigned int j = 1; j < pack_objects[i].size()-1; j++)
				//	ofile<<"f "<<global_vtx_idx[i][0]<<' '<<global_vtx_idx[i][j]<<' '<<global_vtx_idx[i][j+1]<<std::endl;
				ofile<<"f ";
				for (unsigned int j = 0; j < pack_objects[i].size(); j++)
					ofile<<global_vtx_idx[i][j]<<'/'<<global_vtx_idx[i][j]<<'/'<<global_vtx_idx[i][j]<<' ';
				ofile<<std::endl;
			}
		}
		else
		{
			ofile << "mtllib "<<"mat.mtl\n";
			for (unsigned int i = 0; i < pack_objects.size(); i++)
			{
				global_vtx_idx.push_back(std::vector<unsigned int>());
				global_vtx_idx.back().reserve(pack_objects[i].size());
				for (unsigned int j = 0; j < pack_objects[i].size(); j++)
				{
					Point_3 v = pack_objects[i].vertex(j);
					ofile<<"v "<<v.x()<<' '<<v.y()<<' '<<v.z()<<std::endl;
					global_vtx_idx.back().push_back(vtx_idx);
					vtx_idx++;
				}
			}
			// write texture coordinates
			for (unsigned int i = 0; i < pack_objects.size(); i++)
			{
				std::vector<Point_2>& text_coord = pack_objects[i].texture_coord;
				for (int j = 0; j < pack_objects[i].size(); j++)
				{
					ofile<<"vt "<<text_coord[j].x()<<' '<<1.0 - text_coord[j].y()<<std::endl;
				}
			}
			// write normals
			for (unsigned int i = 0; i < pack_objects.size(); i++)
			{
				Vector_3 n = pack_objects[i].norm();
				for (unsigned int j = 0; j < pack_objects[i].size(); j++)
					ofile<<"vn "<<n.x()<<' '<<n.y()<<' '<<n.z()<<std::endl;
			}
			std::vector<std::string> texture_files;
			pio.read_texture_files(texture_files);
			std::vector<std::string> texture_names;
			// pick out only the file name
			// save material file
			for (unsigned int i = 0; i < texture_files.size(); i++)
			{
				unsigned int last_slash_pos = texture_files[i].rfind('/');
				if (last_slash_pos == std::string::npos)
					last_slash_pos = texture_files[i].rfind('\\');
				if (last_slash_pos == std::string::npos)
					last_slash_pos = 0;
				std::string fn(texture_files[i], last_slash_pos+1);
				unsigned int dot_pos = fn.find('.');
				std::string texture_name(fn, 0, dot_pos);
				texture_names.push_back(texture_name);
			}
			// write faces
			for (unsigned int i = 0; i < pack_objects.size(); i++)
			{
				ofile<<"usemtl "<<texture_names[pack_objects[i].texture_id]<<std::endl;
				ofile<<"f ";
				for (unsigned int j = 0; j < pack_objects[i].size(); j++)
					ofile<<global_vtx_idx[i][j]<<'/'<<global_vtx_idx[i][j]<<'/'<<global_vtx_idx[i][j]<<' ';
				ofile<<std::endl;
			}
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
 		if (vector_field && !align_constraint)
 			n = 3;
 		else
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
 		if (vector_field && !align_constraint)
 		{
 			xLoBnds[0] = 0.0;
 			xUpBnds[0] = KTR_INFBOUND;
 			for (i = 1; i < n; i++) {
 				xLoBnds[i] = -KTR_INFBOUND;
 				xUpBnds[i] = KTR_INFBOUND;
 			}
 			/* initial point */
 			//xInitial[0] = *io_k;
			xInitial[0] = 0.9;
 			xInitial[1] = *io_t1;
 			xInitial[2] = *io_t2;
 		}
 		else
 		{
			xLoBnds[0] = 0.0;
			xUpBnds[0] = KTR_INFBOUND;
			xLoBnds[1] = rot_lower_bd;
			xUpBnds[1] = rot_upper_bd;
			for (i = 2; i < n; i++) {
				xLoBnds[i] = -KTR_INFBOUND;
				xUpBnds[i] = KTR_INFBOUND;
			}
			/* initial point */
			//xInitial[0] = *io_k;
			xInitial[0] = 0.9;
			xInitial[1] = *io_theta;
			xInitial[2] = *io_t1;
			xInitial[3] = *io_t2;
		}

		for (j = 0; j < m; j++) {
			cType[j] = KTR_CONTYPE_GENERAL;
			cLoBnds[j] = 0.0;
			cUpBnds[j] = KTR_INFBOUND;
		}

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
		if (vector_field)
		{
			if (KTR_set_func_callback(kc, vcallback) != 0)
				return -1 ;
			if (KTR_set_grad_callback(kc, vcallback) != 0)
				return -1 ;
		}
		else
		{
			if (KTR_set_func_callback(kc, callback) != 0)
				return -1 ;
			if (KTR_set_grad_callback(kc, callback) != 0)
				return -1 ;
		}

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

		if(nStatus == 0)
		{	
			*io_k =  x[0];
			if (vector_field && !align_constraint)
			{
 				*io_theta = 0.0;
 				*io_t1 = x[1];
 				*io_t2 = x[2];
			}
			else
			{
				*io_theta = x[1];
				*io_t1 = x[2];
				*io_t2 = x[3];
			}
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
			*obj    = x[0];
#ifdef _CILK_
			p->containments[*(unsigned int*)userParams].compute_constraints(x, c, m);
#else
			p->containment.compute_constraints(x,c,m);
#endif
			return 0;
		}
		else if (evalRequestCode == KTR_RC_EVALGA)
		{
			objGrad[0] = 1.0;
			objGrad[1] = 0.0;
			objGrad[2] = 0.0;
			objGrad[3] = 0.0;

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

	int Packer::vcallback(const int evalRequestCode, const int n, const int m, const int nnzJ, const int nnzH,
							const double * const x,	const double * const lambda, double * const obj, double * const c,
							double * const objGrad,	double * const jac,	double * const hessian,	double * const hessVector, void * userParams)
	{
		Packer* p = instance() ;
		double f = p->rot_upper_bd*p->rot_upper_bd;
		Containment *ctm;
#ifdef _CILK_
		ctm = &(p->containments[*(unsigned int*)userParams]);
#else
		ctm = &p->containment;
#endif
		if (evalRequestCode == KTR_RC_EVALFC)
		{
			if (!p->align_constraint)
			{
				*obj = x[0];
				ctm->vcompute_constraints(x, c, m);
			}
			else
			{
				//*obj = x[0] + (f*4)/((x[1] - ctm->align_angle)*(x[1] - ctm->align_angle) + f);
				*obj = x[0] + (f*2)/((x[1] - ctm->align_angle)*(x[1] - ctm->align_angle) + f) + (f*2)/((x[1] - PI + ctm->align_angle)*(x[1] - PI + ctm->align_angle) + f);
				ctm->compute_constraints(x, c, m);
			}
			return 0;
		}
		else if (evalRequestCode == KTR_RC_EVALGA)
		{
			if (p->align_constraint)
			{
				objGrad[0] = 1.0;
				//objGrad[1] = 0.0;
				double temp = (x[1] - ctm->align_angle)*(x[1] - ctm->align_angle) + f;
				temp *= temp;
				objGrad[1] = (f*2*2*(x[1] - ctm->align_angle))/temp;
				temp = (x[1] - PI + ctm->align_angle)*(x[1] - PI + ctm->align_angle) + f;
				temp *= temp;
				objGrad[1] += (f*2*2*(x[1] - PI + ctm->align_angle))/temp;
				objGrad[2] = 0.0;
				objGrad[3] = 0.0;
				ctm->compute_constraint_grads(x, jac);
			}
			else
			{
				objGrad[0] = 1.0;
				objGrad[1] = 0.0;
				objGrad[2] = 0.0;
				ctm->vcompute_constraint_grads(x, jac);
			}

			//ctm->vcompute_constraint_grads(x, jac);
			//ctm->compute_constraint_grads(x, jac);

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
		if (mesh_segments.size() == 0)
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

	void Packer::sortHoleEdge(const Hole& hole, std::vector<Vertex_handle>& sorted_loop)
	{
		std::map< Vertex_handle, std::vector<Vertex_handle> > vtx_nbs;

		for (int i=0; i<hole.size(); i++)
		{
			const std::pair<Vertex_handle, Vertex_handle>& edge = hole[i];
			vtx_nbs[edge.first].push_back(edge.second);
			vtx_nbs[edge.second].push_back(edge.first);
		}

		sorted_loop.push_back(hole[0].first);
		sorted_loop.push_back(hole[0].second);
		bool finish = false;
		while (!finish)
		{
			Vertex_handle vh = sorted_loop.back();
			Vertex_handle prev_vh = sorted_loop[sorted_loop.size()-2];
			for (std::vector<Vertex_handle>::iterator nb_itr = vtx_nbs[vh].begin(); nb_itr != vtx_nbs[vh].end(); nb_itr++)
			{
				if ((*nb_itr) != prev_vh)
				{
					if (*nb_itr == sorted_loop[0])
					{
						finish = true;
					}
					else
					{
						sorted_loop.push_back(*nb_itr);
					}

					break;
				}
			}
		}
	}
}