#ifndef CONVERSIONS
#define CONVERSIONS

#include <sstream>
#include <boost/multiprecision/mpfr.hpp>
#include <ginac/ginac.h>

namespace mp = boost::multiprecision;
namespace gi = GiNaC;

namespace ricpad {
namespace conversions {

using mp::mpfr_float;
using std::cout;
using std::endl;

gi::ex mp_to_ex( const mpfr_float& );

template <typename T>
T ex_to_mp(const gi::ex& num) 
{
    std::ostringstream oss;
    oss << num;

    return T(oss.str());
};

template <typename T>
T assign_h(const mpfr_float &tol) {
    cout << "assign_h is not defined for the selected typename." << endl;
    return -1;
};

}; // namespace conversions
}; // namespace ricpad

#endif
