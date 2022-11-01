
// Based on https://gcc.gnu.org/wiki/Visibility
#if defined _WIN32 || defined __CYGWIN__
    #ifdef __GNUC__
        #define DLL_EXPORT __attribute__ ((dllexport))
    #else
        #define DLL_EXPORT __declspec(dllexport)
    #endif
#else
    #define DLL_EXPORT __attribute__ ((visibility ("default")))
#endif

#include <dolfin/common/Array.h>
#include <dolfin/math/basic.h>
#include <dolfin/mesh/SubDomain.h>
#include <Eigen/Dense>


// cmath functions
using std::cos;
using std::sin;
using std::tan;
using std::acos;
using std::asin;
using std::atan;
using std::atan2;
using std::cosh;
using std::sinh;
using std::tanh;
using std::exp;
using std::frexp;
using std::ldexp;
using std::log;
using std::log10;
using std::modf;
using std::pow;
using std::sqrt;
using std::ceil;
using std::fabs;
using std::floor;
using std::fmod;
using std::max;
using std::min;

const double pi = DOLFIN_PI;


namespace dolfin
{
  class dolfin_subdomain_b52c1a2e7ea15353a421996e34af254a : public SubDomain
  {
     public:
       

       dolfin_subdomain_b52c1a2e7ea15353a421996e34af254a()
          {
            
          }

       // Return true for points inside the sub domain
       bool inside(const Eigen::Ref<const Eigen::VectorXd> x, bool on_boundary) const final
       {
         return on_boundary && x[0]>1.5.0 && x[0]<9.5 && x[1]>1.0 && x[1]<5.0;
       }

       void set_property(std::string name, double value)
       {

       }

       double get_property(std::string name) const
       {

         return 0.0;
       }

  };
}

extern "C" DLL_EXPORT dolfin::SubDomain * create_dolfin_subdomain_b52c1a2e7ea15353a421996e34af254a()
{
  return new dolfin::dolfin_subdomain_b52c1a2e7ea15353a421996e34af254a;
}

