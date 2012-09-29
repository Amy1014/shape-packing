#include <Geex/symbolics/function.h>
namespace Geex {

   class P4: public Function {
      public:
      P4();
      virtual void eval(bool do_f, bool do_g, bool do_H) ;
   };
}
