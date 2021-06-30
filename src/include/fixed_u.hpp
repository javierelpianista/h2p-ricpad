#ifndef COUPLING
#define COUPLING 

#endif

#include <iostream>
#include <vector>

#include <boost/multiprecision/mpfr.hpp>

#include <ricpad/conversions.hpp>
#include <ricpad/solver.hpp>
#include <ricpad/hankdet.hpp>

namespace mp = boost::multiprecision;

using mp::mpfr_float;

template<typename num_t>
void U(
    const int Dmin,
    const int Dmax,
    const int m,
    const int s, 
    const num_t &U0, 
    const num_t &A0,
    const num_t &R,
    const int d = 0,
    mpfr_float tol = -1,
    num_t h = -1, 
    const bool use_E = false,
    const std::string output_file = "",
    const bool log_nr = false, 
    const int nr_max_iter = 100
    )

{
    using ricpad::hankdet::hankdet;
    using ricpad::solver::NR_solve;
    using ricpad::conversions::assign_h;

    num_t U = U0, A = A0, Aold, Up;

    if ( use_E ) {
        U = U + 1/R;
        // Add this to the variable U to print what we want
        Up = -1/R;
    } else {
        Up = 0;
    }

    int D = Dmin;

    std::function<num_t(Eigen::Matrix<num_t,1,1>&)> fun;

    fun = [&D, &d, &m, &A, &R]
        ( Eigen::Matrix<num_t,1,1>& param ) -> num_t { 
            num_t &U = param[0];

            std::vector<num_t> plv, qlv, coefs;
            plv = coefficients::pl<num_t>(4*D+d+2);
            qlv = coefficients::ql<num_t>(4*D+d+2, U, A, R, m);

            coefs = coefsl<num_t>(4*D+d-1, U, A, R, m, plv, qlv);
            coefs.erase(coefs.begin(), coefs.begin()+d+2);
            num_t ans = hankdet<num_t>(2*D,coefs);

            return ans;
    };
    
    std::vector<decltype(fun)> F;
    F.push_back(fun);

    std::ofstream out_stream;

    Eigen::Matrix<num_t,1,1> params, params_old;
    Eigen::Matrix<mpfr_float,1,1> dif_vec;
    mpfr_float dif;

    params(0) = U;

    for ( ; D<=Dmax || Dmax == -1; D++ ) {
        params_old = params;

        if ( output_file != "" ) {
            out_stream.open(output_file, std::ios::app);
        }
        params = NR_solve<num_t, mpfr_float, 1>(F, params, tol, h, 
            nr_max_iter, log_nr, out_stream);

        cout.precision(30);
        cout << " D = " << std::setw(2) << D << " ";
        cout << params(0) + Up << " " << std::setprecision(4) << dif << endl;

        if ( output_file != "" ) {
            out_stream << " D = " 
                << std::left
                << std::setw(3) << D << " "
                << std::setw(params(0).precision() + 10)
                << std::setprecision(params(0).precision()) 
                << params(0) + Up << endl;

            out_stream.close();
        }

        dif_vec = params - params_old;
        dif = 1;

        dif = mp::abs(dif_vec(0));
    }
}
