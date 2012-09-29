#include "dist_3d.h"
#include <Geex/symbolics/const_pow.h>
namespace Geex {
   dist_3d::dist_3d() : Function(1,9,9){}
   void dist_3d::eval(bool do_f, bool do_g, bool do_H) {
      if(do_f) {
{
double tmp_2 = x(5);
double tmp_3 = x(8);
double tmp_4 = -tmp_3;
double tmp_5 =  tmp_4+tmp_2;
double tmp_7 = p(1);
double tmp_8 = x(6);
double tmp_9 = x(3);
double tmp_10 = -tmp_9;
double tmp_11 =  tmp_10+tmp_8;
double tmp_14 = x(7);
double tmp_15 = x(4);
double tmp_16 = -tmp_15;
double tmp_17 =  tmp_16+tmp_14;
double f0 =  ( tmp_7*tmp_17-p(2)*tmp_5+tmp_11*p(0))*tmp_11-( p(5)*tmp_5-tmp_7*tmp_11-p(4)*tmp_17)*tmp_17-( p(6)*tmp_11+p(7)*tmp_17-p(8)*tmp_5)*tmp_5;
f(0) = f0;
}
      }
      if(do_g) {
{
double g0_0 = 0.0;
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
double tmp_13 = x(5);
double tmp_14 = x(8);
double tmp_15 = -tmp_14;
double tmp_16 =  tmp_15+tmp_13;
double g0_3 =  -2.0*p(1)*( x(7)-x(4))+p(6)*tmp_16+-2.0*( x(6)-x(3))*p(0)+p(2)*tmp_16;
g(0,3) = g0_3;
}
{
double tmp_2 = x(5);
double tmp_3 = x(8);
double tmp_4 = -tmp_3;
double tmp_5 =  tmp_2+tmp_4;
double g0_4 =  p(5)*tmp_5+-2.0*p(1)*( x(6)-x(3))+-2.0*p(4)*( x(7)-x(4))+p(7)*tmp_5;
g(0,4) = g0_4;
}
{
double tmp_2 = x(7);
double tmp_3 = x(4);
double tmp_4 = -tmp_3;
double tmp_5 =  tmp_4+tmp_2;
double tmp_8 = x(6);
double tmp_9 = x(3);
double tmp_10 = -tmp_9;
double tmp_11 =  tmp_10+tmp_8;
double g0_5 = -tmp_11*p(6)-p(7)*tmp_5-p(2)*tmp_11-p(5)*tmp_5+2.0*p(8)*( x(5)-x(8));
g(0,5) = g0_5;
}
{
double tmp_13 = x(5);
double tmp_14 = x(8);
double tmp_15 = -tmp_14;
double tmp_16 =  tmp_15+tmp_13;
double g0_6 =  2.0*p(1)*( x(7)-x(4))-tmp_16*p(6)+2.0*( x(6)-x(3))*p(0)-tmp_16*p(2);
g(0,6) = g0_6;
}
{
double tmp_2 = x(5);
double tmp_3 = x(8);
double tmp_4 = -tmp_3;
double tmp_5 =  tmp_4+tmp_2;
double g0_7 = -p(5)*tmp_5+2.0*p(1)*( x(6)-x(3))+2.0*p(4)*( x(7)-x(4))-p(7)*tmp_5;
g(0,7) = g0_7;
}
{
double tmp_2 = x(7);
double tmp_3 = x(4);
double tmp_4 = -tmp_3;
double tmp_5 =  tmp_2+tmp_4;
double tmp_8 = x(6);
double tmp_9 = x(3);
double tmp_10 = -tmp_9;
double tmp_11 =  tmp_8+tmp_10;
double g0_8 =  tmp_11*p(6)+tmp_5*p(5)+tmp_11*p(2)+p(7)*tmp_5+-2.0*p(8)*( x(5)-x(8));
g(0,8) = g0_8;
}
      }
   }
}
