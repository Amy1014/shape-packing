#pragma once

#include <Geex/graphics/geexob.h>
#include "packer.h"
#include "spm_graphics.h"

namespace Geex 
{
	class SPM : public Geexob, public Packer, public SPM_Graphics 
	{
	public:
		SPM();
		virtual void get_bbox(
			real& x_min, real& y_min, real& z_min,
			real& x_max, real& y_max, real& z_max
			) ;

		virtual void do_draw() ;
		void load_project(const std::string& prj_file)
		{
			Packer::load_project(prj_file);
			SPM_Graphics::load_graphics();
		}
	} ;

}//namespace Geex
