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
    const bool use_E = false
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

            /*
            cout << "-----------------------------------------" << endl;
            cout  << "U" << " "<< U << endl;
            cout  << "A" << " "<< A << endl;
            cout  << "R" << " "<< R << endl;
            cout  << "Dm"<< " " << ans << endl;
            cout << "-----------------------------------------" << endl;
            cout << endl;
            */

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

            /*
            cout << "-----------------------------------------" << endl;
            cout << "U" <<  " " << U << endl;
            cout << "A" <<  " " << A << endl;
            cout << "R" <<  " " << R << endl;
            cout << "Dl" << " " <<  ans << endl;
            cout << "-----------------------------------------" << endl;
            cout << endl;
            */

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

    std::ofstream U_stream, A_stream, R_stream;

    for ( ; D<=Dmax; D++ ) {
        params_old = params;
        params = NR_solve<num_t, mpfr_float, 2>(F, params, tol, h);

        cout.precision(30);
        cout << " D = " << std::setw(2) << D << " ";
        cout << params(0).real() + Up.real() << " " << params(1).real() << 
            " " << std::setprecision(4) << dif << endl;

        U_stream.open("U.dat", std::ios::app);
        U_stream << std::setprecision(params(0).precision()) << 
            params(0).real() + Up.real() << endl;
        U_stream.close();

        A_stream.open("A.dat", std::ios::app);
        A_stream << std::setprecision(params(1).precision()) << params(1).real()
            << endl;
        A_stream.close();

        dif_vec = params - params_old;
        dif = 1;

        for ( int i = 0; i<2; i++ ) {
            dif = min(dif, mp::abs(dif_vec(i)));
        }
    }
}
#endif
