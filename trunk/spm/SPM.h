#pragma once

#include <Geex/graphics/geexob.h>
#include "Doc.h"
#include "View.h"

namespace Geex 
{
	class SPM : public Geexob, public Doc, public View 
	{
	public:
		SPM();
		virtual void get_bbox(
			real& x_min, real& y_min, real& z_min,
			real& x_max, real& y_max, real& z_max
			) ;

		virtual void do_draw() ;

	} ;

}//namespace Geex
