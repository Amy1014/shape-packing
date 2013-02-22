
#include "SPM.h"

namespace Geex
{
	SPM::SPM(): View(this)
	{

	}
	void SPM::do_draw()
	{
		View::draw();
	}
	void SPM::get_bbox(real& x_min, real& y_min, real& z_min,
					real& x_max, real& y_max, real& z_max) 
	{
		Doc::get_bbox(x_min, y_min, z_min, x_max, y_max, z_max);
	}
}