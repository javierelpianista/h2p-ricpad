#ifndef RICPAD_SOLVER
#define RICPAD_SOLVER

#include <iostream>

#include <boost/multiprecision/mpfr.hpp>
#include <boost/math/policies/error_handling.hpp>

#include <ricpad/differentiate.hpp>
#include <ricpad/hankdet.hpp>

namespace mp = boost::multiprecision;
using mp::mpfr_float;
using std::cout;
using std::endl;

namespace ricpad::solver { 
template <
    typename C, // complex number type
    typename R, // real number type
    int N       // Number of equations and unknowns
>
Eigen::Matrix<C,N,1> NR_solve( 
        const std::vector<std::function<C(Eigen::Matrix<C,N,1>&)>> &f, 
        Eigen::Matrix<C,N,1> x0,
        const R& tol, 
        const C& h,
        const int maxiter = 100,
        const bool log_iters = false,
        std::ostream &out_stream = cout
        ) 
{
    using ricpad::differentiate::differentiate;
    Eigen::Matrix<C, N, N> jacobian, inv_jacobian;
    Eigen::Matrix<C, N, 1> x(x0), xold;
    R desv = tol + 1;
    static const char* function = "ricpad::solver::NR_solver";

    int niter = 0;

    C val;

    while ( desv > tol ) {
        for ( int i = 0; i < x.size(); i++ ) {
            for ( int j = 0; j < x.size(); j++ ) {
                val = differentiate<C, N>(f[i], x, j, h);
                jacobian(i, j) = std::move(val);
            }
        }

        inv_jacobian = jacobian.inverse();
        
        Eigen::Matrix<C, N, 1> F;

        for ( int i = 0; i < N; i++ ) F(i) = f[i](x); 

        xold = x;
        x = x - inv_jacobian * F;

        if ( log_iters ) {
            out_stream << " << ";
            for ( int i = 0; i < N; i++ ) {
                out_stream << std::setw(x0(0).precision() + 10) 
                    << std::left 
                    << std::setprecision(x0(0).precision())
                    << x(i);
            }
            out_stream << endl;
        }

        desv = (x - xold).norm();

        if ( niter++ > maxiter ) {
            return boost::math::policies::raise_evaluation_error(
                function, 
                "Maximum number of Newton-Raphson iterations reached", x0, 
                boost::math::policies::policy<>());
        }
    }

    return x;
} 

} // namespace ricpad::solver

#endif
