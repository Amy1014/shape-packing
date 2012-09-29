#include "P0.h"
#include <Geex/symbolics/const_pow.h>
namespace Geex {
   P0::P0() : Function(3,0,3){}
   void P0::eval(bool do_f, bool do_g, bool do_H) {
      if(do_f) {
{
double f0 = p(0);
f(0) = f0;
}
{
double f1 = p(1);
f(1) = f1;
}
{
double f2 = p(2);
f(2) = f2;
}
      }
      if(do_g) {
      }
   }
}
