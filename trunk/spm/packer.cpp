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

		rpvd.set_mesh(pio.attribute_value("MeshFile"));
		rpvd.set_trimesh(&mesh);

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

	void Packer::compute_clipped_VD()
	{
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

	void Packer::generate_RDT()
	{
		rpvd.begin_insert();
		rpvd.insert_polygons(pack_objects.begin(), pack_objects.end(), 12);
		rpvd.insert_bounding_points(10);
		rpvd.end_insert();
		std::cout<<"Start computing clipped Voronoi region\n";
		compute_clipped_VD();
		std::cout<<"End computing.\n";
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

	Packer::Optimization_res Packer::optimize_one_polygon(unsigned int id, Local_frame& lf, Parameter& solution)
	{
		// choose a local frame
		lf.o = pack_objects[id].centroid();
		lf.w = pack_objects[id].norm(); // assume this vector has already been normalized
		Vector_3 u(lf.o, pack_objects[id].vertex(0));
		lf.u = u/CGAL::sqrt(u.squared_length());
		lf.v = CGAL::cross_product(lf.w, lf.u);

		// formulate optimization 
#ifdef _CILK_
		containments[id].clear();
#else
		containment.clear();
#endif
		const RestrictedPolygonVoronoiDiagram::VertGroup& samp_pnts = rpvd.sample_points_group(id);
		for (unsigned int i = 0; i < samp_pnts.size(); i++)
		{
			const vector<Point_3>& bisec_pnts = samp_pnts[i]->vd_vertices;
			Point_2 ref_point = lf.to_uv(samp_pnts[i]->point_3());
			if (bisec_pnts.size() == 0)
				continue;
			assert(bisec_pnts.size() != 1);
			for (unsigned int j = 0; j < bisec_pnts.size()-1; j++)
			{
				const Point_3& s = bisec_pnts[j];
				const Point_3& t = bisec_pnts[j+1];
#ifdef _CILK_
				//containments[id].const_list_.push_back(Constraint(lf.to_uv(samp_pnts[i]->point_3()), lf.to_uv(Segment_3(s,t)), ref_point));
				containments[id].const_list_.push_back(Constraint(lf.to_uv(samp_pnts[i]->point_3()), lf.to_uv(Segment_3(s, t))));
#else
				containment.const_list_.push_back(Constraint(lf.to_uv(samp_pnts[i]->point_3()), lf.to_uv(Segment_3(s,t)), ref_point));
				//containment.const_list_.push_back(Constraint(lf.to_uv(samp_pnts[i]->point_3()), lf.to_uv(Segment_3(s, t))));
#endif
			}
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

	Packer::Lloyd_res Packer::lloyd(bool enlarge)
	{
		std::vector<Local_frame> lfs(pack_objects.size());
		std::vector<Parameter> solutions(pack_objects.size());
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
		if (mink >= 1.005)
		{
			for (unsigned int i = 0; i < pack_objects.size(); i++)
				if (opti_res[i] == SUCCESS)
					solutions[i].k = mink;
		}
		constraint_transformation(solutions, lfs);
#ifdef _CILK_
		cilk_for (unsigned int i = 0; i < pack_objects.size(); i++)
#else
		for (unsigned int i = 0; i < pack_objects.size(); i++)
#endif
		{
			Apply_transformation t(solutions[i], lfs[i]);
			std::transform(pack_objects[i].vertices_begin(), pack_objects[i].vertices_end(), pack_objects[i].vertices_begin(), t);
			//RestrictedPolygonVoronoiDiagram::VertGroup& samp_pnts = rpvd.sample_points_group(i);
			//std::transform(samp_pnts.begin(), samp_pnts.end(), samp_pnts.begin(), t);
			Point_3 c = pack_objects[i].centroid();
			vec3 dummyv, p;
			int fid;
			p = mesh.project_to_mesh(to_geex_pnt(c), dummyv, fid);
			vec3 n = approx_normal(fid);
			Transformation_3 same_t = pack_objects[i].align(to_cgal_vec(n), to_cgal_pnt(p));
			//for (unsigned int j = 0; j < samp_pnts.size(); j++)
			//{
				//vec3 dummyv;
				//Point_3 temp = same_t(samp_pnts[j]->point_3());
				//samp_pnts[j]->mp = to_cgal_pnt(mesh.project_to_mesh(to_geex_pnt(temp), dummyv));
			//}
		}
		for (RestrictedPolygonVoronoiDiagram::Face_iterator fit = rpvd.faces_begin(); fit != rpvd.faces_end(); ++fit)
		{
			int i0 = fit->vertex(0)->group_id, i1 = fit->vertex(1)->group_id, i2 = fit->vertex(2)->group_id;
			if (i0 < 0 || i1 < 0 || i2 < 0)
				continue;
			Point_3 p0 = fit->vertex(0)->mp, p1 = fit->vertex(1)->mp, p2 = fit->vertex(2)->mp;
			Vector_3 v01(p0, p1), v12(p1, p2);
			Vector_3 v = CGAL::cross_product(v01, v12);
			Vector_3 ov = fit->n;
			//std::cout<<v*ov<<std::endl;
			//system("pause");
			//if (v*ov < 0.0)
			//	system("pause");
		}
		if (enlarge && mink < 1.005)
			return NO_MORE_ENLARGEMENT;
		else
			return MORE_ENLARGEMENT;
		//generate_RDT();
	}

	void Packer::pack(void (*post_action)())
	{
		for (unsigned int i = 0; i < 1; i++)
		{
			lloyd(false);
			generate_RDT();
			//rpvd.iDT_update();
			//compute_clipped_VD();
			if (post_action != NULL)
				post_action();
		}
		//Lloyd_res r;
		//do 
		//{
		//	r = lloyd(true);
		//	generate_RDT();
		//	if (post_action != NULL)
		//		post_action();
		//} while ( r == MORE_ENLARGEMENT);
	}
	
	void Packer::rpack(void (*post_action)())
	{
		Lloyd_res r;
		r = lloyd(true);
		//generate_RDT();
		update_iDT();
		if (post_action != NULL)
			post_action();
	}

	void Packer::constraint_transformation(vector<Parameter>& parameters, vector<Local_frame>& lfs)
	{
		vector<double> min_factors(parameters.size(), 1.0); // minimum constraint control factors
		typedef RestrictedPolygonVoronoiDiagram RPVD;
		for (RPVD::Face_iterator fit = rpvd.faces_begin(); fit != rpvd.faces_end(); ++fit)
		{
			int i0 = fit->vertex(0)->group_id, i1 = fit->vertex(1)->group_id, i2 = fit->vertex(2)->group_id;
			if (i0 < 0 || i1 < 0 || i2 < 0) // involving boundary edges
				continue;
			// original normal
			//Vector_3 n01(fit->vertex(0)->point_3(), fit->vertex(1)->point_3());
			//Vector_3 n12(fit->vertex(1)->point_3(), fit->vertex(2)->point_3());
			Vector_3 n01(fit->vertex(0)->mp, fit->vertex(1)->mp);
			Vector_3 n12(fit->vertex(1)->mp, fit->vertex(2)->mp);
			Vector_3 on = CGAL::cross_product(n01, n12);
			fit->n = on;
			double m = 1.0;
			unsigned int N = 1000, n = 1;
			while (n < N)
			{
				Apply_transformation t0(parameters[i0]*m, lfs[i0]), t1(parameters[i1]*m, lfs[i1]), t2(parameters[i2]*m, lfs[i2]);
				vec3 dummyn;
				//Point_3 tp0 = t0(fit->vertex(0)->point_3());
				//Point_3 tp1 = t1(fit->vertex(1)->point_3());
				//Point_3 tp2 = t2(fit->vertex(2)->point_3());
				Point_3 tp0 = to_cgal_pnt(mesh.project_to_mesh(to_geex_pnt(t0(fit->vertex(0)->point_3())), dummyn));
				Point_3 tp1 = to_cgal_pnt(mesh.project_to_mesh(to_geex_pnt(t1(fit->vertex(1)->point_3())), dummyn));
				Point_3 tp2 = to_cgal_pnt(mesh.project_to_mesh(to_geex_pnt(t2(fit->vertex(2)->point_3())), dummyn));
				Vector_3 tn = CGAL::cross_product(Vector_3(tp0, tp1), Vector_3(tp1, tp2));
				if (tn*on > 0.0)
					break;
				m /= 2.0;
				n++;
			}
			if (n == N)	m = 0.0;
			min_factors[i0] = std::min(m, min_factors[i0]);
			min_factors[i1] = std::min(m, min_factors[i1]);
			min_factors[i2] = std::min(m, min_factors[i2]);
		}
		for (unsigned int i = 0; i < parameters.size(); i++)
		{
			if (min_factors[i] < 1.0)
				parameters[i] *= min_factors[i];
			if (min_factors[i] == 0.0)
				parameters[i].k = 1.0;// don't allow the polygon to change
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
}