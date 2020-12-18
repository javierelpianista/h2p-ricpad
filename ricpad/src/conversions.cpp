#include <sstream>
#include <ricpad/conversions.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/mpc.hpp>

namespace mp = boost::multiprecision;
using mp::mpfr_float;
using mp::mpc_complex;

namespace ricpad::conversions {

gi::ex mp_to_ex ( const mpfr_float& num ) {
    std::ostringstream oss;
    oss << std::setprecision(num.precision()) << num;

    gi::parser reader;

    return reader(oss.str());
};

gi::ex mp_to_ex ( const mpc_complex& num ) {
    std::ostringstream oss;
    oss << std::setprecision(num.precision()) << num.real() 
        << "+" << num.imag() << "*I";

    gi::parser reader;

    gi::ex ans = reader(oss.str());

    return reader(oss.str());
};

template <>
mpc_complex ex_to_mp<mpc_complex>( const gi::ex& numex ) {
    std::ostringstream oss;
    gi::numeric num = gi::ex_to<gi::numeric>(numex);

    oss << std::setprecision(gi::Digits) << 
        "(" << num.real() << "," << num.imag() << ")" << endl;

    mpc_complex ans(oss.str());

    return ans;
};

template <>
mpfr_float assign_h(const mpfr_float &tol) {
    return tol;
}

template <>
mpc_complex assign_h(const mpfr_float &tol) {
    return mp::sqrt(mpfr_float(2))*tol*mpc_complex("(1,1)");
}
};
