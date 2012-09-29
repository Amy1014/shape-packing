#include <Geex/symbolics/function.h>
namespace Geex {

   class dist_3d: public Function {
      public:
      dist_3d();
      virtual void eval(bool do_f, bool do_g, bool do_H) ;
   };
}
