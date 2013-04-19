
#include "SPM.h"

namespace Geex
{
	SPM::SPM(): SPM_Graphics(this)
	{
		contain_multi_pack_ = false; 
		pack_process_id = 0;
	}
	void SPM::do_draw()
	{
		SPM_Graphics::draw();
	}
	void SPM::get_bbox(real& x_min, real& y_min, real& z_min,
					real& x_max, real& y_max, real& z_max) 
	{
		Packer::get_bbox(x_min, y_min, z_min, x_max, y_max, z_max);
	}
}