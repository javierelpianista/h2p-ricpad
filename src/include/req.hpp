#ifndef REQ
#define REQ

#include <iostream>
#include <vector>

#include <boost/multiprecision/mpfr.hpp>

#include <ricpad/conversions.hpp>
#include <ricpad/solver.hpp>
#include <ricpad/hankdet.hpp>

#include <h2p.hpp>

namespace mp = boost::multiprecision;

using mp::mpfr_float;

template<typename num_t>
void R_eq(
    const int Dmin,
    const int Dmax,
    const int m, 
    const int s, 
    const num_t &U0, 
    const num_t &A0,
    const num_t &R0,
    // Use a scheme of automatic determination of numerical precision in each step
    const int d = 0,
    mpfr_float tol = -1,
    num_t h = -1, 
    num_t h2 = -1, 
    const std::string output_file = "",
    const bool log_nr = false,
    const int nr_max_iter = 100,
    const bool double_D_lambda = false
) {
    using ricpad::hankdet::hankdet;
    using ricpad::solver::NR_solve;
    using ricpad::conversions::assign_h;

    int ndigits;

    num_t U = U0, R = R0, A = A0;

    int D = Dmin, doub;

    if ( double_D_lambda ) {
        doub = 2;
    } else {
        doub = 1;
    };

    cout << "s = " << s << endl;
    cout << "m = " << m << endl;
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

    fun_l = [&D, &d, &m, &doub]
        ( Eigen::Matrix<num_t,3,1>& param ) -> num_t { 
            num_t &U = param[0];
            num_t &A = param[1];
            num_t &R = param[2];

            std::vector<num_t> plv, qlv, coefs;
            plv = coefficients::pl<num_t>(doub*(2*D+1)+d);
            qlv = coefficients::ql<num_t>(doub*(2*D+1)+d, U, A, R, m);

            coefs = coefsl<num_t>(doub*2*D+d-1, U, A, R, m, plv, qlv);
            coefs.erase(coefs.begin(), coefs.begin()+d+2);
            num_t ans = hankdet<num_t>(doub*D,coefs);

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

    std::ofstream out_stream;

    for ( ; D<=Dmax || Dmax == -1; D++ ) {
        params_old = params;
        if ( output_file != "" ) {
            out_stream.open(output_file, std::ios::app);
        }
        params = NR_solve<num_t, mpfr_float, 3>(F, params, tol, h,
                nr_max_iter, log_nr, out_stream);

        U = mpfr_float(params(0));
        A = mpfr_float(params(1));
        R = mpfr_float(params(2));

        cout.precision(30);
        cout << " D = " << std::setw(2) << D << " ";
        cout << params(0).real() << " " << params(1).real() << " " << 
            params(2).real() << " " << std::setprecision(4) << dif << endl;

        if ( output_file != "" ) {
            out_stream << " D = " 
                << std::left
                << std::setw(3) << D << " "
                << std::setw(params(0).precision() + 10)
                << std::setprecision(params(0).precision()) 
                << params(0)  
                << std::setw(params(1).precision() + 10)
                << std::setprecision(params(1).precision()) 
                << params(1)  
                << std::setw(params(2).precision() + 10)
                << std::setprecision(params(2).precision()) 
                << params(2)  
                << endl;

            out_stream.close();
        }

        dif_vec = params - params_old;
        dif = 1;

        for ( int i = 0; i<3; i++ ) {
            dif = min(dif, mp::abs(dif_vec(i)));
        }

        for ( int i = 0; i<2; i++ ) {
            dif = min(dif, mp::abs(dif_vec(i)));
        }
    }
}
#endif
