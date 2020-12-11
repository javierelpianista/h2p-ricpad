#include <sstream>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/mpc.hpp>
#include <ginac/ginac.h>

namespace mp = boost::multiprecision;
namespace gi = GiNaC;

namespace ricpad {
namespace conversions {

using mp::mpfr_float;
using mp::mpc_complex;
using std::cout;
using std::endl;

gi::ex mp_to_ex( const mpfr_float& );
gi::ex mp_to_ex( const mpc_complex& );

template <class T>
T ex_to_mp(const gi::ex& num) 
{
    std::ostringstream oss;
    oss << num;

    return T(oss.str());
};

}; // namespace conversions
}; // namespace ricpad
