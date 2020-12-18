#ifndef REQ
#define REQ

#include <iostream>
#include <vector>

#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/mpc.hpp>

#include <ricpad/conversions.hpp>
#include <ricpad/solver.hpp>
#include <ricpad/hankdet.hpp>

#include <h2p.hpp>

namespace mp = boost::multiprecision;

using mp::mpfr_float;
using mp::mpc_complex;

template<typename num_t>
void R_eq(
    const int Dmin,
    const int Dmax,
    const num_t &U0, 
    const num_t &A0,
    const num_t &R0,
    // Use a scheme of automatic determination of numerical precision in each step
    const int d = 0,
    const bool auto_prec = true,
    mpfr_float tol = -1,
    num_t h = -1, 
    num_t h2 = -1, 
    const mpfr_float &prec_ini = -1
) {
    using ricpad::hankdet::hankdet;
    using ricpad::solver::NR_solve;
    using ricpad::conversions::assign_h;

    int ndigits;

    if ( auto_prec ) {
        if ( prec_ini == -1 ) {
            tol = mpfr_float("1E-10");
        } else {
            tol = prec_ini;
        }
        h = assign_h<num_t>(tol*tol);
        h2 = h*h;

        ndigits = (-mp::log10(tol)*8).convert_to<int>();
        num_t::default_precision(ndigits);
        mpfr_float::default_precision(ndigits);
        gi::Digits = ndigits;
    } else {
        ndigits = num_t::default_precision();
    }

    num_t U = U0, R = R0, A = A0;

    int D = Dmin, m = 0, s = 0;

    std::function<num_t(Eigen::Matrix<num_t,3,1>&)> 
        fun_m, fun_l, fun_d; 

    fun_m = [&D, &d, &s, &m]
        ( Eigen::Matrix<num_t,3,1>& param ) -> num_t {
            num_t &U = param[0];
            num_t &A = param[1];
            num_t &R = param[2];

            std::vector<num_t> pmv, qmv, coefs; 
            pmv = coefficients::pm<num_t>(2*D+d+1);
            qmv = coefficients::qm<num_t>(2*D+d+1, U, A, R, m);

            coefs = coefsm<num_t>(2*D+d-1, s, pmv, qmv);
            coefs.erase(coefs.begin(), coefs.begin()+d+1);
            num_t ans = hankdet<num_t>(D,coefs);

            return ans;
    };

    fun_l = [&D, &d, &m]
        ( Eigen::Matrix<num_t,3,1>& param ) -> num_t { 
            num_t &U = param[0];
            num_t &A = param[1];
            num_t &R = param[2];

            std::vector<num_t> plv, qlv, coefs;
            plv = coefficients::pl<num_t>(2*D+d+2);
            qlv = coefficients::ql<num_t>(2*D+d+2, U, A, R, m);

            coefs = coefsl<num_t>(2*D+d-1, U, A, R, m, plv, qlv);
            coefs.erase(coefs.begin(), coefs.begin()+d+2);
            num_t ans = hankdet<num_t>(D,coefs);

            return ans;
    };
    
    fun_d = [&D, &d, &m, &s, &h2, &fun_m, &fun_l]
        ( Eigen::Matrix<num_t,3,1>& param ) -> num_t {
            using ricpad::differentiate::differentiate;

            num_t &U = param[0];
            num_t &A = param[1];
            num_t &R = param[2];

            return (
                differentiate<num_t, 3>(fun_m, param, 1, h2) *
                  differentiate<num_t, 3>(fun_l, param, 2, h2) -
                differentiate<num_t, 3>(fun_m, param, 2, h2) * 
                  differentiate<num_t, 3>(fun_l, param, 1, h2) 
                  );
        };

    std::vector<decltype(fun_m)> F;

    F.push_back(fun_m);
    F.push_back(fun_l);
    F.push_back(fun_d);

    Eigen::Matrix<num_t,3,1> params, params_old;
    Eigen::Matrix<mpfr_float,3,1> dif_vec;
    mpfr_float dif;

    params(0) = U;
    params(1) = A;
    params(2) = R;

    std::ofstream U_stream, A_stream, R_stream;
    
    for ( ; D<=Dmax; D++ ) {
        params_old = params;
        params = NR_solve<num_t, mpfr_float, 3>(F, params, tol, h);
        U = mpfr_float(params(0));
        A = mpfr_float(params(1));
        R = mpfr_float(params(2));

        cout.precision(30);
        cout << " D = " << std::setw(2) << D << " ";
        cout << params(0).real() << " " << params(1).real() << " " << 
            params(2).real() << " " << std::setprecision(4) << dif << endl;

        U_stream.open("U.dat", std::ios::app);
        U_stream << std::setprecision(params(0).precision()) << params(0).real()
            << endl;
        U_stream.close();

        A_stream.open("A.dat", std::ios::app);
        A_stream << std::setprecision(params(1).precision()) << params(1).real()
            << endl;
        A_stream.close();

        R_stream.open("R.dat", std::ios::app);
        R_stream << std::setprecision(params(2).precision()) << params(2).real()
            << endl;
        R_stream.close();

        dif_vec = params - params_old;
        dif = 1;

        for ( int i = 0; i<3; i++ ) {
            dif = min(dif, mp::abs(dif_vec(i)));
        }

        if ( auto_prec ) {
            tol = min(tol, dif*1E-10);
            h = assign_h<num_t>(tol*tol);
            h2 = h*h;

            ndigits = (-mp::log10(tol)*16).convert_to<int>();
            num_t::default_precision(ndigits);
            mpfr_float::default_precision(ndigits);
            gi::Digits = ndigits;
        }
    }
}
#endif
