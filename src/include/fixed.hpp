#ifndef FIXED
#define FIXED

#include <iostream>
#include <vector>

#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/mpc.hpp>

#include <ricpad/conversions.hpp>
#include <ricpad/solver.hpp>
#include <ricpad/hankdet.hpp>

namespace mp = boost::multiprecision;

using mp::mpfr_float;
using mp::mpc_complex;

template<typename num_t>
void UA(
    const int Dmin,
    const int Dmax,
    const int m,
    const int s, 
    const num_t &U0, 
    const num_t &A0,
    const num_t &R,
    // Use a scheme of automatic determination of numerical precision in each step
    const int d = 0,
    mpfr_float tol = -1,
    num_t h = -1, 
    num_t h2 = -1,
    const bool use_E = false,
    const std::string output_file = "",
    const bool log_nr = false, 
    const int nr_max_iter = 100
    )

{
    using ricpad::hankdet::hankdet;
    using ricpad::solver::NR_solve;
    using ricpad::conversions::assign_h;

    num_t U = U0, A = A0, Up;

    if ( use_E ) {
        U = U + 1/R;
        // Add this to the variable U to print what we want
        Up = -1/R;
    } else {
        Up = 0;
    }

    int D = Dmin;

    std::function<num_t(Eigen::Matrix<num_t,2,1>&)> 
        fun_m, fun_l; 

    fun_m = [&D, &d, &s, &m, &R]
        ( Eigen::Matrix<num_t,2,1>& param ) -> num_t {
            num_t &U = param[0];
            num_t &A = param[1];

            std::vector<num_t> pmv, qmv, coefs; 
            pmv = coefficients::pm<num_t>(2*D+d+1);
            qmv = coefficients::qm<num_t>(2*D+d+1, U, A, R, m);

            coefs = coefsm<num_t>(2*D+d-1, s, pmv, qmv);
            coefs.erase(coefs.begin(), coefs.begin()+d+1);
            num_t ans = hankdet<num_t>(D,coefs);

            return ans;
    };

    fun_l = [&D, &d, &m, &R]
        ( Eigen::Matrix<num_t,2,1>& param ) -> num_t { 
            num_t &U = param[0];
            num_t &A = param[1];

            std::vector<num_t> plv, qlv, coefs;
            plv = coefficients::pl<num_t>(2*D+d+2);
            qlv = coefficients::ql<num_t>(2*D+d+2, U, A, R, m);

            coefs = coefsl<num_t>(2*D+d-1, U, A, R, m, plv, qlv);
            coefs.erase(coefs.begin(), coefs.begin()+d+2);
            num_t ans = hankdet<num_t>(D,coefs);

            return ans;
    };
    
    std::vector<decltype(fun_m)> F;

    F.push_back(fun_m);
    F.push_back(fun_l);

    Eigen::Matrix<num_t,2,1> params, params_old;
    Eigen::Matrix<mpfr_float,2,1> dif_vec;
    mpfr_float dif;

    params(0) = U;
    params(1) = A;

    std::ofstream out_stream;

    for ( ; D<=Dmax || Dmax == -1; D++ ) {
        params_old = params;

        if ( output_file != "" ) {
            out_stream.open(output_file, std::ios::app);
        }
        params = NR_solve<num_t, mpfr_float, 2>(F, params, tol, h, 
            nr_max_iter, log_nr, out_stream);

        cout.precision(30);
        cout << " D = " << std::setw(2) << D << " ";
        cout << params(0)+ Up<< " " << params(1)<< 
            " " << std::setprecision(4) << dif << endl;

        if ( output_file != "" ) {
            out_stream << " D = " 
                << std::left
                << std::setw(3) << D << " "
                << std::setw(params(0).precision() + 10)
                << std::setprecision(params(0).precision()) 
                << params(0) + Up << params(1) << endl; 

            out_stream.close();
        }

        dif_vec = params - params_old;
        dif = 1;

        for ( int i = 0; i<2; i++ ) {
            dif = min(dif, mp::abs(dif_vec(i)));
        }
    }
}
#endif
