#include "Init_point_generator.h"

namespace Geex
{

	inline double Init_point_generator::weight(unsigned int facet_idx)
	{
		double cur = mesh.curvature_at_face(facet_idx);
		const Facet& f = mesh[facet_idx];
		double area = f.area();
		return std::sqrt(area)*1.0*density(cur); 
	}

	/************************************************************************/
	/*					Uniformly sampling                                  */
	/************************************************************************/
	Point_3 Uniform_point_generator::next_point()
	{
		if (nb_pnt_fetched > nb_pnts)
		{
			std::cout<<"Too many times of fetching!!!\n";
			return Point_3(0.0, 0.0, 0.0);
		}
		
		int rand_idx = ::rand()%mesh.size();

		while (facets_used.find(rand_idx) != facets_used.end())
			rand_idx = ::rand()%mesh.size();
		facets_used.insert(rand_idx);

		const Facet& f = mesh[rand_idx];
		nb_pnt_fetched++;

		return CGAL::centroid(to_cgal_pnt(f.vertex[0]), to_cgal_pnt(f.vertex[1]), to_cgal_pnt(f.vertex[2]));
	}

	/************************************************************************/
	/*					Sampling based on density                           */
	/************************************************************************/
	Density_point_generator::Density_point_generator(const TriMesh& trimesh, unsigned int nb_init_pnts) : Init_point_generator(trimesh, nb_init_pnts)
	{
		pnts.reserve(nb_init_pnts);
		std::set<int> init_facets;

		double total_weight = 0.0, res = 0.0;
		std::vector<double> region_area(mesh.nb_vertices());
		for (unsigned int i = 0; i < mesh.nb_vertices(); i++)
		{
			double cur = mesh.curvature_at_vertex(i);
			//total_weight += std::sqrt(cur_size_map(cur));
			double neighbor_area = 0.0;
			for (unsigned int j = 0; j < mesh.vertex(i).faces_.size(); j++)
				neighbor_area += std::fabs(mesh[mesh.vertex(i).faces_.at(j)].area());
			region_area[i] = neighbor_area;
			//region_area[i] = std::pow(region_area[i], 0.75);
			//total_weight += density(cur)*region_area[i];
			total_weight += mesh.vertex(i).weight_;
		}
		for (unsigned int i = 0; i < mesh.nb_vertices(); i++)
		{
			if (mesh.near_boundary(i) || mesh.is_on_feature(i))
				continue;
			double cur = mesh.curvature_at_vertex(i);
			//double nf = nb_init_pnts*(density(cur)*region_area[i]/total_weight);
			//double nf = nb_init_polygons*(std::sqrt(cur_size_map(cur))/total_weight);
			double nf = nb_init_pnts*(mesh.vertex(i).weight_/total_weight);
			int n(nf+res);
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
					pnts.push_back(CGAL::midpoint(to_cgal_pnt(mesh.vertex(i).pos_), c));
					init_facets.insert(idx);
				}
			}
			else 
				res += nf;
		}
		// pick up the lost ones
		while ( pnts.size() < nb_init_pnts && init_facets.size() < mesh.size() )
		{
			int idx = ::rand()%mesh.size();
			while ( init_facets.find(idx) != init_facets.end() )
				idx = ::rand()%mesh.size();
			init_facets.insert(idx);
			int vi0 = mesh[idx].vertex_index[0], vi1 = mesh[idx].vertex_index[1], vi2 = mesh[idx].vertex_index[2];
			bool v_near_bd[] = {mesh.is_on_boundary(vi0), mesh.is_on_boundary(vi1), mesh.is_on_boundary(vi2)};
			if ( v_near_bd[0] /*&& v_near_bd[1]*/ || v_near_bd[1] /*&& v_near_bd[2]*/ || v_near_bd[2] /*&& v_near_bd[0]*/ )
				continue;
			Point_3 c = CGAL::centroid(to_cgal_pnt(mesh[idx].vertex[0]), to_cgal_pnt(mesh[idx].vertex[1]), to_cgal_pnt(mesh[idx].vertex[2]));
			pnts.push_back(c);
		}

	}

	Point_3 Density_point_generator::next_point()
	{
		if (nb_pnt_fetched > pnts.size())
		{
			std::cout<<"Too many times of fetching!!!\n";
			nb_pnt_fetched = 0;
		}
		return pnts[nb_pnt_fetched++];
	}

	/************************************************************************/
	/*						Sampling after CVT                              */
	/************************************************************************/

	CVT_point_generator::CVT_point_generator(const TriMesh& trimesh, unsigned int nb_init_pnts, const std::string& mesh_filename) 
											: Init_point_generator(trimesh, nb_init_pnts)
	{
// 		try
// 		{

		// generate initial points based on curvature
		pnts.reserve(nb_init_pnts);
		topo_mesh = new TopoPolyMesh;
		topo_mesh->load(mesh_filename);

		double total_weight = 0.0;
		std::vector<double> facet_weights;
		facet_weights.reserve(mesh.size());
		std::vector<double> facet_density;
		facet_density.reserve(mesh.size());
		for (unsigned int i = 0; i < mesh.size(); i++)
		{
			facet_weights.push_back(weight(i));
			facet_density.push_back(density(i));
			total_weight += facet_weights.back();
		}

		double res = 0.0;
		for (unsigned int i = 0; i < mesh.size(); i++)
		{
			const Facet& f = mesh[i];
			double nf = nb_init_pnts*facet_weights[i]/total_weight;
			int n(nf+res);
			if (n>=1)
			{
				res = nf + res - n;
				for (int j = 0; j < n; j++)
				{
					double alpha = Geex::Numeric::random_float64();
					double beta = Geex::Numeric::random_float64();
					vec3 p = alpha*f.vertex[0] + beta*f.vertex[1] + (1.0 - alpha - beta)*f.vertex[2];
					pnts.push_back(p);
				}
			}
			else
				res += nf;
		}
		while (pnts.size() < nb_init_pnts)
		{
			int rand_fidx = Geex::Numeric::random_int32()%mesh.size();
			const Facet& f = mesh[rand_fidx];
			double alpha = Geex::Numeric::random_float64();
			double beta = Geex::Numeric::random_float64();
			vec3 p = alpha*f.vertex[0] + beta*f.vertex[1] + (1.0 - alpha - beta)*f.vertex[2];
			pnts.push_back(p);
		}

		const int lloyd_times = 1;
		std::cout<<"Starting lloyd, do it "<<lloyd_times<<" times...\n";
		for (int i = 0; i < lloyd_times; i++)
		{
			Delaunay *del = Delaunay::create("CGAL");
			del->set_vertices(pnts);
			RVD rvd(del, topo_mesh);
			Compute_centroids cc(mesh, facet_density, pnts.size());
			rvd.for_each_facet(cc);
			cc.get_centroids(pnts);
			for (unsigned int j = 0; j < pnts.size(); j++)
			{
				vec3 n ;
				pnts[j] = mesh.project_to_mesh(pnts[j], n);
			}
			delete del;
		}
		std::cout<<"End lloyd.\n";
// 		}
// 		catch (...)
// 		{
// 			prompt_and_exit("Unknown exception happens.");
// 		}
	}

	Point_3 CVT_point_generator::next_point()
	{
		if (nb_pnt_fetched > pnts.size())
		{
			std::cout<<"Too many times of fetching!!!\n";
			nb_pnt_fetched = 0;
		}
		return to_cgal_pnt(pnts[nb_pnt_fetched++]);
	}

	unsigned int CVT_point_generator::nb_points() 
	{
		return pnts.size();
	}

	Compute_centroids::Compute_centroids(const TriMesh& _mesh, const std::vector<double>& _facets_density, unsigned int nb_seeds)
											: mesh(_mesh), facets_density(_facets_density)
	{
		denominator_sum.resize(nb_seeds, 0.0);
		numerator_sum.resize(nb_seeds, vec3(0.0, 0.0, 0.0));
	}

	void Compute_centroids::operator()(unsigned int vi, TopoPolyMesh *clipped_region) const
	{
		for (unsigned int f = 0; f < clipped_region->nb_facets(); f++)
		{
			std::vector<Point_3> verts;
			verts.reserve(clipped_region->facet_end(f) - clipped_region->facet_begin(f));
			for (unsigned int i = clipped_region->facet_begin(f); i < clipped_region->facet_end(f); i++)
				verts.push_back(to_cgal_pnt(clipped_region->vertex(i)));
			Point_3 facet_cent = CGAL::centroid(verts.begin(), verts.end(), CGAL::Dimension_tag<0>());
			// area of this facet
			double facet_area = 0.0;
			for (unsigned int i = 0; i < verts.size(); i++)
			{
				Vector_3 v0(facet_cent, verts[i]), v1(facet_cent, verts[(i+1)%verts.size()]);
				Vector_3 v = CGAL::cross_product(v0, v1);
				facet_area += CGAL::sqrt(v.squared_length());
			}
			facet_area /= 2.0;
			// get density
			vec3 n;
			int fidx;
			vec3 geex_cent = to_geex_pnt(facet_cent);
			mesh.project_to_mesh(geex_cent, n, fidx);
			double dens = facets_density[fidx];
			numerator_sum[vi] = numerator_sum[vi] + (dens * facet_area) * geex_cent;
			denominator_sum[vi] += dens * facet_area;
		}
	}

	void Compute_centroids::get_centroids(std::vector<vec3>& centroids) const
	{
		centroids.resize(numerator_sum.size());
		for (unsigned int i = 0; i < centroids.size(); i++)
			centroids[i] = numerator_sum[i] / denominator_sum[i];
	}


}