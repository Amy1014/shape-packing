#include "P2P.h"
#include <Geex/symbolics/const_pow.h>
namespace Geex {
   P2P::P2P() : Function(3,3,0){}
   void P2P::eval(bool do_f, bool do_g, bool do_H) {
      if(do_f) {
{
double f0 = x(0);
f(0) = f0;
}
{
double f1 = x(1);
f(1) = f1;
}
{
double f2 = x(2);
f(2) = f2;
}
      }
      if(do_g) {
{
double g0_0 = 1.0;
g(0,0) = g0_0;
}
{
double g0_1 = 0.0;
g(0,1) = g0_1;
}
{
double g0_2 = 0.0;
g(0,2) = g0_2;
}
{
double g1_0 = 0.0;
g(1,0) = g1_0;
}
{
double g1_1 = 1.0;
g(1,1) = g1_1;
}
{
double g1_2 = 0.0;
g(1,2) = g1_2;
}
{
double g2_0 = 0.0;
g(2,0) = g2_0;
}
{
double g2_1 = 0.0;
g(2,1) = g2_1;
}
{
double g2_2 = 1.0;
g(2,2) = g2_2;
}
      }
   }
}
