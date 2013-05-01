#include "packer.h"

namespace Geex
{
	Packer* Packer::instance_ = nil;

	double Packer::PI = 3.141592653589793;

	Packer::Packer() /*: phy_disc_barr(0.3, 1.2, 10), opt_disc_barr(0.2, 1.1 , 10), disc_barr(phy_disc_barr, opt_disc_barr)*/
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

		sub_pack_id = 0;

		sync_opt = true;
		phony_upper_scale = 100.0;
		use_voronoi_cell_ = true;

		replace_factor = 0.8;
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
		{
			disc_barr.set(phy_disc_barr, opt_disc_barr);
			max_scale = disc_barr.get_max();
			min_scale = disc_barr.get_min();
			levels = disc_barr.nb_levels();
		}
		else
		{
			max_scale = std::numeric_limits<double>::max();
			min_scale = std::numeric_limits<double>::min();
			levels = 0;
		}
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
			//pgn_lib = pgn_lib_sets[sub_pack_id];
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
		
		// compute mesh area
		mesh_area = 0.0;
		for (unsigned int i = 0; i < mesh.size(); i++)
			mesh_area += std::fabs(mesh[i].area());
		std::cout<<"Surface area: "<<mesh_area<<std::endl;

		// compute maximum and average polygon area in the library
		double max_pgn_area = std::numeric_limits<double>::min(), mean_pgn_area = 0.0;
		for (unsigned int i = 0; i < pgn_lib.size(); i++)
		{
			double s = std::fabs(pgn_lib[i].area());
			max_pgn_area = std::max(max_pgn_area, s);
			mean_pgn_area += s;
		}
		mean_pgn_area /= pgn_lib.size();
		std::cout<<"Maximum polygon area: "<<max_pgn_area<<std::endl;
		std::cout<<"Mean polygon area: "<<mean_pgn_area<<std::endl;

		if (mesh_segments.size() == 0)
			rpvd.set_mesh(pio.attribute_value("MeshFile"));
		else
			rpvd.set_mesh(pio.get_submesh_files().at(sub_pack_id));
		rpvd.set_trimesh(&mesh);

		// distribute polygons by different strategies
		// 1. random 
		unsigned int nb_init_polygons = mesh_area / mean_pgn_area * 0.85;
		std::cout<<nb_init_polygons<<" polygons are expected.\n";
		random_init_tiles(nb_init_polygons);
		std::cout<<pack_objects.size()<<" polygons are loaded initially.\n";

		std::for_each(pack_objects.begin(), pack_objects.end(), std::mem_fun_ref(&Packing_object::activate));
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
				if (mesh.near_boundary(i) || mesh.is_on_feature(i))
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
			while ( init_pos.size() < nb_init_polygons && init_facets.size() < mesh.size() )
			{
				int idx = ::rand()%mesh.size();
				while ( init_facets.find(idx) != init_facets.end() )
					idx = ::rand()%mesh.size();
				init_facets.insert(idx);
				int vi0 = mesh[idx].vertex_index[0], vi1 = mesh[idx].vertex_index[1], vi2 = mesh[idx].vertex_index[2];
				bool v_near_bd[] = { mesh.is_on_boundary(vi0), mesh.is_on_boundary(vi1), mesh.is_on_boundary(vi2)};
				if ( v_near_bd[0] && v_near_bd[1] || v_near_bd[1] && v_near_bd[2] || v_near_bd[2] && v_near_bd[0] )
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
		//if (approx)
		//	rpvd.compute_approx_VD(tangent_planes, ref_pnts);
		//else
		//	rpvd.compute_clipped_VD(tangent_planes, ref_pnts);
		rpvd.compute_clipped_VD(tangent_planes, ref_pnts);
		rpvd.compute_midpoint_VD(tangent_planes, ref_pnts);
		//std::cout<<"End computing.\n";
	}

	void Packer::compute_clipped_VD(std::vector<bool> use_approx)
	{
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
		rpvd.compute_mixed_VD(tangent_planes, ref_pnts, use_approx);
	}
	void Packer::generate_RDT()
	{
		rpvd.begin_insert();
		rpvd.insert_polygons(pack_objects.begin(), pack_objects.end(), samp_nb);
		//rpvd.insert_polygons(pack_objects.begin(), pack_objects.end());
		rpvd.insert_bounding_points(40);
		CGAL::Timer t;
		t.start();
		rpvd.end_insert();
		t.stop();
		std::cout<<"time spent in RDT: "<<t.time()<<" seconds."<<std::endl;
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
				assert(bisec_pnts.size() != 1);
				for (unsigned int j = 0; j < med_segs.size(); j++)
				{
					Segment_2 s2d(lf.to_uv(med_segs[j]));
					ctm.const_list_.push_back(Constraint(lf.to_uv(samp_pnts[i]->point()), s2d, ref_point));
				}
			}
/*			for (unsigned int i = 0; i < pack_objects[id].size(); i++)
			{
				for (unsigned int j = 0; j < samp_pnts.size(); j++)
				{
					const vector<Segment_3>& med_segs = samp_pnts[j]->med_segs;
					if (med_segs.size() == 0)
						continue;
					assert(bisec_pnts.size() != 1);
					for (unsigned int k = 0; k < med_segs.size(); k++)
					{
						ctm.const_list_.push_back(Constraint(lf.to_uv(pack_objects[id].vertex(i)), lf.to_uv(med_segs[k])));
					}
				}
			}*/	
		}
	

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
		double min_scalor;
		//if (discrete_scaling)
		//	min_scalor = 1.01;
		//else
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

		//if (!stop_update_DT)
		//{
			// first constrain rotation and translation
			//constraint_transformation(solutions, lfs, false);
			// second constrain enlargement
			//if (mink < min_scalor && enlarge)
			//	constraint_transformation(solutions, lfs, true);
			//else
				constraint_transformation(solutions, lfs, false);
				//if (discrete_scaling)
				//	constraint_transformation(solutions, lfs, true);
		//}
		//rpvd.save_triangulation("before_transform.obj");
		//std::ofstream of("parameters.txt");
#ifdef _CILK_
		cilk_for (unsigned int i = 0; i < pack_objects.size(); i++)
#else
		for (unsigned int i = 0; i < pack_objects.size(); i++)
#endif
		{
			curv_constrained_transform(solutions[i], pack_objects[i].facet_idx, i);
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
									double barrier_factor, double nxt_barrier/*, std::vector<bool>& approx_vd*/)
	{
		double min_factor = 1.05; // to determine convergence

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
		// check whether the tile has already been close to its voronoi region
		for (unsigned int i = 0; i < pack_objects.size(); i++)
		{
			//if (opti_res[i] == SUCCESS && solutions[i].k >= 1.065)
			//	approx_vd[i] = false;
			//else
			//	approx_vd[i] = true;

			if (!enlarge)
				solutions[i].k = 1.0;
		}
	
		double mink = std::numeric_limits<double>::max();
		//int nb_failures = 0;
		// suppress the inactive and optimization failure ones
		for (unsigned int i = 0; i < pack_objects.size(); i++)
		{
			if (opti_res[i] == SUCCESS && pack_objects[i].active)
				mink = std::min(solutions[i].k, mink);
			else if (opti_res[i] == SUCCESS && !pack_objects[i].active)
				solutions[i].k = 1.0;
			else if (opti_res[i] != SUCCESS)
			{
				//nb_failures++;
				solutions[i].k = 1.0;
				solutions[i].theta = solutions[i].tx = solutions[i].ty = 0.0;
			}
		}

		//std::cout<<"number of failures: "<<nb_failures<<std::endl;
		
		//if (mink >= min_factor)
		//{
		//	for (unsigned int i = 0; i < pack_objects.size(); i++)
		//		if (opti_res[i] == SUCCESS && pack_objects[i].active)
		//			solutions[i].k = mink;
		//}
		constraint_transformation(solutions, lfs, false);

		std::cout<<"Minimum scale is "<<mink<<std::endl;
		if (mink <= 1.1)
			constraint_transformation(solutions, lfs, true);

//#ifdef _CILK_
//		cilk_for (unsigned int i = 0; i < pack_objects.size(); i++)
//			curv_constrained_transform(solutions[i], pack_objects[i].facet_idx, i);
//#else
//		for (unsigned int i = 0; i < pack_objects.size(); i++)
//			curv_constrained_transform(solutions[i], pack_objects[i].facet_idx, i);
//#endif
		// check state of each tile
		std::vector<bool> has_growth_room(pack_objects.size());
		if (enlarge)
			for (unsigned int i = 0; i < pack_objects.size(); i++)
			{
				pack_objects[i].reach_barrier = false;
				if (pack_objects[i].active && pack_objects[i].factor*solutions[i].k > barrier_factor)
				{
					if (pack_objects[i].factor < nxt_barrier)
						solutions[i].k = barrier_factor / pack_objects[i].factor; // restrict it at this barrier
					else
						solutions[i].k = 1.0;
					has_growth_room[i] = false;
					pack_objects[i].reach_barrier = true;
				}
				else if (pack_objects[i].active && solutions[i].k < min_factor)
					has_growth_room[i] = false;
				else if (pack_objects[i].active && solutions[i].k >= min_factor)
					has_growth_room[i] = true;
				else 
					has_growth_room[i] = false;
			}		
		else
			for (unsigned int i = 0; i < pack_objects.size(); i++)
				pack_objects[i].reach_barrier = true;
#ifdef _CILK_
		cilk_for (unsigned int i = 0; i < pack_objects.size(); i++)
#else
		for (unsigned int i = 0; i < pack_objects.size(); i++)
#endif
		{	
			curv_constrained_transform(solutions[i], pack_objects[i].facet_idx, i);
			transform_one_polygon(i, lfs[i], solutions[i]);
		}
		
		bool stay_at_this_barrier = false;
		int nb_potential_growth = 0;
		for (unsigned int i = 0; i < has_growth_room.size(); i++)
		{
			//stay_at_this_barrier = stay_at_this_barrier || has_growth_room[i];
			if (has_growth_room[i])
				nb_potential_growth++;
		}
		double perc = nb_potential_growth/(double)pack_objects.size();
		std::cout<<"percentage of potential growth: "<<perc*100.0<<'%'<<std::endl;
		stay_at_this_barrier = (perc > 0.005); // if only few tiles can be enlarged
		if (!enlarge)
			return true;
		else
			return stay_at_this_barrier;

	}

	void Packer::split_large_tiles()
	{
		std::set<unsigned int> indices_to_split;
		for (unsigned int i = 0; i < pack_objects.size(); i++)
		{
			if ( disc_barr.beyond_range() && pack_objects[i].reach_barrier )
				indices_to_split.insert(i);
		}
		split(indices_to_split);
	}
	void Packer::split(std::set<unsigned int>& tile_indices)
	{
		std::cout<<"Start splitting...\n";
		if (tile_indices.size() == 0)
		{
			std::cout<<"Nothing to split.\n";
			return;
		}
		std::vector<Packing_object> new_tile_set;
		new_tile_set.reserve(tile_indices.size()*2 + pack_objects.size());
		unsigned int split_nb = 0;
		for (unsigned int i = 0; i < pack_objects.size(); i++)
		{
			if (tile_indices.find(i) != tile_indices.end()) // this tile needs splitting
			{
				split_nb++;
				std::vector<Point_3>& smoothed_vd_region = rpvd.get_smoothed_voronoi_regions().at(i);
				Point_3 c = CGAL::centroid(smoothed_vd_region.begin(), smoothed_vd_region.end(), CGAL::Dimension_tag<0>());
				Local_frame lf = pack_objects[i].local_frame();
				std::vector<Point_2> vd_region_2d(smoothed_vd_region.size());
				for (unsigned int j = 0; j < smoothed_vd_region.size(); j++)
					vd_region_2d[j] = lf.to_uv(smoothed_vd_region[j]);
				Point_2 c_2d = lf.to_uv(c);
				unsigned int min_idx;
				double min_dist = std::numeric_limits<double>::max();
				for (unsigned int j = 0; j < vd_region_2d.size(); j++)
				{
					Point_2 midp = CGAL::midpoint(vd_region_2d[j], vd_region_2d[(j+1)%vd_region_2d.size()]);
					double dist = CGAL::squared_distance(midp, c_2d);
					if (dist < min_dist)
					{
						min_idx = j;
						min_dist = dist;
					}
				}
				Point_2 end_pnt0 = CGAL::midpoint(vd_region_2d[min_idx], vd_region_2d[(min_idx+1)%vd_region_2d.size()]), end_pnt1(1.0e20, 1.0e20);
				Line_2 sep_ln(c_2d, end_pnt0);
				for (unsigned int j = 0; j < vd_region_2d.size(); j++)
				{
					if (j == min_idx)
						continue;
					Segment_2 edge(vd_region_2d[j], vd_region_2d[(j+1)%vd_region_2d.size()]);
					if (CGAL::do_intersect(sep_ln, edge))
					{
						CGAL::Object o = CGAL::intersection(sep_ln, edge);
						if ( const Point_2 *inter_pnt = CGAL::object_cast<Point_2>(&o))
							end_pnt1 = *inter_pnt;
					}
				}
				if (end_pnt1 == Point_2(1.0e20, 1.0e20))
					prompt_and_exit("No intersection with the opposite edge.");

				std::vector<Point_2> sep_pnts[2];
				for (unsigned int j = 0; j < vd_region_2d.size(); j++)
				{
					if (sep_ln.has_on_negative_side(vd_region_2d[j]))
						sep_pnts[0].push_back(vd_region_2d[j]);
					else
						sep_pnts[1].push_back(vd_region_2d[j]);
				}
				sep_pnts[0].push_back(end_pnt0);
				sep_pnts[0].push_back(end_pnt1);
				sep_pnts[1].push_back(end_pnt0);
				sep_pnts[1].push_back(end_pnt1);
				for (int k = 0; k < 2; k++)
				{
					std::vector<Point_2> con_region_pnts;
					CGAL::ch_graham_andrew(sep_pnts[k].begin(), sep_pnts[k].end(), std::back_insert_iterator<std::vector<Point_2>>(con_region_pnts));
					// fill the left
					std::vector<Segment_2> closed_region;
					closed_region.reserve(con_region_pnts.size());
					for (unsigned int j = 0; j < con_region_pnts.size(); j++)
						closed_region.push_back(Segment_2(con_region_pnts[j], con_region_pnts[(j+1)%con_region_pnts.size()]));
					Polygon_matcher pm(closed_region, 200);
					std::priority_queue<Match_info_item<unsigned int>, std::vector<Match_info_item<unsigned int>>, Match_measure> match_res;
					for (unsigned int idx = 0; idx < pgn_lib.size(); idx++)
						match_res.push(pm.affine_match(pgn_lib[idx], idx, 0.0));

					double shrink_factor = 0.6;
					Match_info_item<unsigned int> matcher = match_res.top();
					if ( matcher.scale * shrink_factor > disc_barr.get_max() )
						shrink_factor = disc_barr.get_max() / matcher.scale;
					const Ex_polygon_2& match_pgn = pgn_lib[matcher.val];
					Polygon_2 transformed_pgn2d = CGAL::transform(matcher.t, match_pgn);

					Packing_object new_tile;

					for (unsigned int j = 0; j < transformed_pgn2d.size(); j++)
					{
						Point_3 p = lf.to_xy(transformed_pgn2d.vertex(j));
						new_tile.push_back(p);
					}

					Point_3 new_c = new_tile.centroid();
					vec3 v;
					int fid;
					vec3 prjp = mesh.project_to_mesh(to_geex_pnt(new_c), v, fid);

					v = approx_normal(fid);
					new_tile.align(to_cgal_vec(v), to_cgal_pnt(prjp));

					Transformation_3 rescalor = Transformation_3(CGAL::TRANSLATION, Vector_3(CGAL::ORIGIN, to_cgal_pnt(prjp))) *
																( Transformation_3(CGAL::SCALING, shrink_factor) *
																	Transformation_3(CGAL::TRANSLATION, Vector_3(to_cgal_pnt(prjp), CGAL::ORIGIN)) );	

					std::transform(new_tile.vertices_begin(), new_tile.vertices_end(), new_tile.vertices_begin(), rescalor);
					//
					new_tile.lib_idx = matcher.val;
					new_tile.factor = matcher.scale*shrink_factor;
					new_tile.facet_idx = fid;
					new_tile.texture_coord.assign(match_pgn.texture_coords.begin(), match_pgn.texture_coords.end());
					new_tile.texture_id = match_pgn.texture_id;
					new_tile_set.push_back(new_tile);
				}
			}
			else
				new_tile_set.push_back(pack_objects[i]);
		}
		pack_objects.clear();
		pack_objects.reserve(new_tile_set.size());
		std::copy(new_tile_set.begin(), new_tile_set.end(), std::back_insert_iterator<std::vector<Packing_object>>(pack_objects));
		generate_RDT();
		compute_clipped_VD();
		stop_update_DT = false;
		std::for_each(pack_objects.begin(), pack_objects.end(), std::mem_fun_ref(&Packing_object::activate));
		double min_size = std::numeric_limits<double>::max();
		for (unsigned int i = 0; i < pack_objects.size(); i++)
			min_size = std::min(pack_objects[i].factor, min_size);
		disc_barr.set_current_barrier(min_size);
		std::cout<<"End splitting."<<split_nb<<" tiles are split.\n";
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
			//std::vector<bool> approx_vd(pack_objects.size());

			std::cout<<"============ discrete lloyd iteration "<<times++<<" ============\n";
			
			double nxt_barrier;
			if (!disc_barr.get_next_barrier(nxt_barrier))
				nxt_barrier = 1.0e20;
			stay_at_this_barrier = discrete_one_lloyd(enlarge, solutions, local_frames, current_barrier, nxt_barrier);
			
			if (!stay_at_this_barrier)  // go to the next barrier
			{
				// roll back to previous barrier for tiles that do not reach current barrier
				//if (!sync_opt)
				//	for (unsigned int j = 0; j < pack_objects.size(); j++)
				//	{
				//		if (pack_objects[j].active && !pack_objects[j].reach_barrier)
				//		{
				//			pack_objects[j].active = false;
				//			//double closest_phy_barrier;
				//			//phy_disc_barr.set_current_barrier(pack_objects[j].factor);
				//			//if (!phy_disc_barr.get_prev_barrier(closest_phy_barrier))
				//			//	continue;
				//			//Parameter p(closest_phy_barrier / pack_objects[j].factor, 0.0, 0.0, 0.0);
				//			//Local_frame lf = pack_objects[j].local_frame();
				//			//transform_one_polygon(j, lf, p);
				//		}
				//	}
				//else
				//{
				//	// check whether any chance to be larger
				//	bool all_reach_barrier = true;
				//	for (unsigned int j = 0; j < pack_objects.size(); j++)
				//		all_reach_barrier = all_reach_barrier && pack_objects[j].reach_barrier;
				//	if (all_reach_barrier)
				//		std::cout<<"Synchronized optimization should stop.\n";
				//}
				disc_barr.go_to_next_barrier();
				std::cout<<"go to next barrier\n";
			}
			else
			{
				bool no_reach_barrier = true;
				for (unsigned int j = 0; j < pack_objects.size(); j++)
					no_reach_barrier = no_reach_barrier && !pack_objects[j].reach_barrier;
				if (no_reach_barrier)
				{
					std::cout<<"No tile can be enlarged any more.\n";
					disc_barr.go_to_end(); // just to terminate this series of iterations using barrier
				}
			}

			// check whether to use approximate voronoi region
			//bool use_appox_VD = false;
			double min_factor = std::numeric_limits<double>::max();
			for (unsigned int j = 0; j < solutions.size(); j++)
				if (solutions[j].k > 1.0)
					min_factor = std::min(solutions[j].k, min_factor);		
			use_voronoi_cell_ = !(enlarge && (min_factor < 1.05));

			// update idt
			rpvd.iDT_update();
			//compute_clipped_VD(approx_vd);

			//if (!use_appox_VD)
				compute_clipped_VD(false);
			//else 
			//	compute_clipped_VD(true);

			if (post_action != NULL)
				post_action();
		}
		return !stay_at_this_barrier;
	}
	
	void Packer::interface_lloyd(void (*post_action)())
	{
		if (!sync_opt)
			lloyd(post_action, false);
		else
			discrete_lloyd(post_action, false);
	}
	void Packer::pack(void (*post_action)())
	{
		
		if (!sync_opt)
			lloyd(post_action, true);
		else
		{
			if (sync_opt && disc_barr.get_max() != phony_upper_scale)
				disc_barr.append(phony_upper_scale);
			if (discrete_lloyd(post_action, true))
			{
				//for (int i = 0; i < 1; i++)
					//discrete_lloyd(post_action, false);
			}
		}
		print_area_coverage();	
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
		
#ifdef _CILK_
		cilk_for (unsigned int i = 0; i < pack_objects.size(); i++)
#else
		for (unsigned int i = 0; i < pack_objects.size(); i++)
#endif
		{
			const RestrictedPolygonVoronoiDiagram::VertGroup& samp_pnts = rpvd.sample_points_group(i);

			Local_frame lf = pack_objects[i].local_frame();

			std::vector<Segment_2> region2d;

			//bool contain_non_delaunay_facet = false;
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
				//contain_non_delaunay_facet |= samp_pnts[j]->contain_non_delaunay_facet;
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
			if (penetration)
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
			//std::cout<<"Restricting size\n";
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
			Point_2 s = lf.to_uv( prj_plane.projection(hl[j].first->point()) );
			Point_2 t = lf.to_uv( prj_plane.projection(hl[j].second->point()) );
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
		double shrink_factor = 0.3;

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

	void Packer::con_replace()
	{
		std::cout<<"Start replacing...\n";
		CGAL::Timer replace_timer;
		replace_timer.start();

#ifdef _CILK_
		//__cilkrts_set_param("nworkers", "24");
		cilk_for (unsigned int i = 0; i < pack_objects.size(); i++)
#else
		for (unsigned int i = 0; i < pack_objects.size(); i++)
#endif
		{
			const RestrictedPolygonVoronoiDiagram::VertGroup& samp_pnts = rpvd.sample_points_group(i);

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
			std::priority_queue<Match_info_item<unsigned int>, std::vector<Match_info_item<unsigned int>>, Match_measure> match_res;
			for (unsigned int idx = 0; idx < pgn_lib.size(); idx++)
				match_res.push(pm.affine_match(pgn_lib[idx], idx, match_weight));

			// choose the result with the smallest match error now
			double shrink_factor = replace_factor;

			Match_info_item<unsigned int> matcher = match_res.top();
			if (sync_opt)
			{
				//Match_info_item<unsigned int> best_backup = matcher;
				//while ( matcher.scale*shrink_factor > discrete_factors.back() )
				//{
				//	match_res.pop();
				//	if (match_res.empty())
				//		break;
				//	matcher = match_res.top();
				//}	
				//if (match_res.empty())
				//{
				//	std::cout<<"No suitable matcher!\n";
				//	shrink_factor = discrete_factors.back() / best_backup.scale;
				//	matcher = best_backup;
				//}
				if ( matcher.scale * shrink_factor > disc_barr.get_max() )
					shrink_factor = disc_barr.get_max() / matcher.scale;
			}
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
		use_voronoi_cell_ = true;
		std::for_each(pack_objects.begin(), pack_objects.end(), std::mem_fun_ref(&Packing_object::activate));
		if (sync_opt)
		{
			double min_size = std::numeric_limits<double>::max();
			//current_factor = discrete_factors.size();
			for (unsigned int i = 0; i < pack_objects.size(); i++)
				min_size = std::min(pack_objects[i].factor, min_size);
			disc_barr.set_current_barrier(min_size);
		}
		std::cout<<"End replacing. Computation time: "<< replace_timer.time()<<" seconds.\n";
	}
	void Packer::report()
	{
		double max_factor = std::numeric_limits<double>::min();
		double min_factor = std::numeric_limits<double>::max();
		double cur_max_factor = std::numeric_limits<double>::min();
		double cur_min_factor = std::numeric_limits<double>::max();
		double mean_factor = 0.0;
		double cur_mean_factor = 0.0;
		std::set<unsigned int> lib_indices;
		for (unsigned int i = 0; i < pack_objects.size(); i++)
		{
			//double normalized_factor = pack_objects[i].factor*pgn_lib[pack_objects[i].lib_idx].normalize_factor;
			const Facet& f = mesh[pack_objects[i].facet_idx];
			double v0_cur = mesh.curvature_at_vertex(f.vertex_index[0]),
				   v1_cur = mesh.curvature_at_vertex(f.vertex_index[1]),
				   v2_cur = mesh.curvature_at_vertex(f.vertex_index[2]);
			double avg_cur = (v0_cur + v1_cur + v2_cur) / 3.0 ;
			double res_area = pack_objects[i].area();
			double normalized_factor = /*std::pow(pack_objects[i].factor, 0.75) * *//*pow(avg_cur, pio.gamma_used())*/avg_cur*avg_cur*res_area;
			max_factor = std::max(max_factor, pack_objects[i].factor);
			min_factor = std::min(min_factor, pack_objects[i].factor);
			cur_max_factor = std::max(cur_max_factor, normalized_factor);
			cur_min_factor = std::min(cur_min_factor, normalized_factor);
			//if (cur_min_factor > cur_max_factor)
			//{
			//	std::cout<<"curvature = "<<avg_cur<<", factor = "<<pack_objects[i].factor<<std::endl;
			//	std::cout<<"( "<<cur_max_factor<<','<<cur_min_factor<<")\n";
			//	system("pause");
			//}
			mean_factor += pack_objects[i].factor;
			cur_mean_factor += normalized_factor;
			lib_indices.insert(pack_objects[i].lib_idx);
		}
		mean_factor /= pack_objects.size();
		cur_mean_factor /= pack_objects.size();
		std::cout<<"Absolute size: \n";
		std::cout<<"\t-- Maximum factor: "<<max_factor<<std::endl;
		std::cout<<"\t-- Minimum factor: "<<min_factor<<std::endl;
		std::cout<<"\t-- Mean factor: "<<mean_factor<<std::endl;
		std::cout<<"Curvature relative size: \n";
		std::cout<<"\t-- Maximum curvature correction factor: "<<cur_max_factor<<std::endl;
		std::cout<<"\t-- Minimum curvature correction factor: "<<cur_min_factor<<std::endl;
		std::cout<<"\t-- Mean factor: "<<cur_mean_factor<<std::endl;
		std::cout<<"Polygon variety:\n";
		std::cout<<"\t-- Number of polygon types from input: "<<pgn_lib.size()<<std::endl;
		std::cout<<"\t-- Number of polygon types: "<<lib_indices.size()<<std::endl;
	}

	void Packer::print_area_coverage()
	{
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
	}
	void Packer::discretize_tiles()
	{
		// compute all scaling levels
		std::vector<double> discrete_scale(levels+1);
		std::cout<<"Discretization scales: < ";
		for (int i = 0; i < levels+1; i++)
		{
			discrete_scale[i] = ( (levels-i)*min_scale + i*max_scale ) / levels;
			std::cout<<discrete_scale[i]<<" ";
		}
		std::cout<<'>'<<std::endl;
		for (unsigned int i = 0; i < pack_objects.size(); i++)
		{
			std::vector<double>::iterator closest_scale_it = std::lower_bound(discrete_scale.begin(), discrete_scale.end(), pack_objects[i].factor);
			double closest_scale = 1.0;
			if (closest_scale_it == discrete_scale.begin())
				std::cout<<"====== Too small scale, no discrete scale found! ======\n";
			else if (*closest_scale_it != pack_objects[i].factor)
			{
				closest_scale = *(closest_scale_it - 1) / pack_objects[i].factor;
				//if (closest_scale > pack_objects[i].factor)
				//	system("pause");
				pack_objects[i] *= closest_scale;
			}
		}
		generate_RDT();
		compute_clipped_VD();
		stop_update_DT = false;
		print_area_coverage();
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
			for (unsigned int i = 0; i < pack_objects.size(); i++)
			{
				for (unsigned int j = 1; j < pack_objects[i].size()-1; j++)
					ofile<<"f "<<global_vtx_idx[i][0]<<' '<<global_vtx_idx[i][j]<<' '<<global_vtx_idx[i][j+1]<<std::endl;
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
}