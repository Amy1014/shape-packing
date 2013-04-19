#pragma once

#include <Geex/graphics/geexob.h>
#include "packer.h"
#include "spm_graphics.h"

namespace Geex 
{
	class SPM : public Geexob, public Packer, public SPM_Graphics 
	{
	public:
		SPM() ;
		virtual void get_bbox(
			real& x_min, real& y_min, real& z_min,
			real& x_max, real& y_max, real& z_max
			) ;

		virtual void do_draw() ;
		void load_project(const std::string& prj_file)
		{
			Packer::load_project(prj_file);
			SPM_Graphics::load_graphics();
			contain_multi_pack_ = get_project_ioer().mesh_polygon_coupled();
		}

		void start_new_packing()
		{
			if (pack_process_id > 0 && pack_process_id <= get_project_ioer().nb_sub_pack_parts() )
				Packer::save_sub_result();
			if (pack_process_id >= 0 && pack_process_id < get_project_ioer().nb_sub_pack_parts())
			{
				Packer::pack_next_submesh();
				SPM_Graphics::load_next_texture_lib();
				redraw_triangulation();
				redraw_voronoi_cell();
				glut_viewer_redraw();
			}
			pack_process_id++;
		}

		bool contain_multi_packing()
		{
			return contain_multi_pack_;
		}

	private:
		bool contain_multi_pack_;
		unsigned int pack_process_id;
	} ;

}//namespace Geex
