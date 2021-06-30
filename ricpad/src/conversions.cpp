#include <sstream>
#include <ricpad/conversions.hpp>
#include <boost/multiprecision/mpfr.hpp>

namespace mp = boost::multiprecision;
using mp::mpfr_float;

namespace ricpad::conversions {

gi::ex mp_to_ex ( const mpfr_float& num ) {
    std::ostringstream oss;
    oss << std::setprecision(num.precision()) << num;

    gi::parser reader;

    return reader(oss.str());
};

template <>
mpfr_float assign_h(const mpfr_float &tol) {
    return tol;
}
}
