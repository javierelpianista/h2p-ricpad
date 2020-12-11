#include <iostream>
#include <vector>
#include <fstream>

#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/mpc.hpp>
#include <boost/multiprecision/eigen.hpp>
#include <ginac/ginac.h>
#include <Eigen/Dense>

#include <ricpad/hankdet.hpp>
#include <ricpad/solver.hpp>
#include <ricpad/conversions.hpp>

using boost::multiprecision::mpfr_float;
using boost::multiprecision::mpc_complex;

using std::cout;
using std::endl;

namespace mp = boost::multiprecision;
namespace gi = GiNaC;

namespace coefficients {
template <typename T>
std::vector<T> pl( int N ) {
    std::vector<T> ans;
    ans.reserve(N+1);
    gi::symbol x("x");
    gi::ex pl = 2*(x+1)/(x*(x+2));

    pl = gi::series_to_poly(pl.series(x, N+3));

    for ( int i = 0; i <= N; i++ ) {
        ans.emplace_back(
                ricpad::conversions::ex_to_mp<T>(
                    gi::ex_to<gi::numeric>(
                        pl.coeff(x, i).evalf()
                        )
                    )
                );
    }

    return ans;
}

template <typename T>
std::vector<T> ql( 
        const int N, 
        const T& U, 
        const T& A, 
        const T& R, 
        const int m ) {

    gi::ex Ee, Ae, Re;
    T E = R/2 - U*R*R/2;

    Ee = ricpad::conversions::mp_to_ex(E);
    Ae = ricpad::conversions::mp_to_ex(A);
    Re = ricpad::conversions::mp_to_ex(R);

    std::vector<T> ans;
    ans.reserve(N+1);
    gi::symbol x("x");
    gi::ex ql = (2*Re*(x+1)+Ae-Ee*gi::power((x+1),2))/(x*(x+2)) 
                    - gi::power(m,2)/(gi::power(x,2)*gi::power((x+2),2));

    ql = gi::series_to_poly(ql.series(x, N+3));

    for ( int i = 0; i <= N; i++ ) {
        ans.emplace_back(
                ricpad::conversions::ex_to_mp<T>(
                    gi::ex_to<gi::numeric>(
                        ql.coeff(x, i).evalf()
                        )
                    )
                );
    }

    return ans;
}


template <typename T>
std::vector<T> pm(
        const int N
        )
{
    std::vector<T> ans;
    ans.reserve(N+1);

    for ( int i = 0; i <= N; i++ ) {
        ans.emplace_back(T(-2));
    };

    return ans;
}


// Return only the coefficients of Qm with even exponents (the odd ones are 0)
template <typename T>
std::vector<T> qm( 
        const int N, 
        const T& U, 
        const T& A, 
        const T& R, 
        const int m ) 
{
    std::vector<T> ans;
    ans.reserve(N+1);
    T E = R/2 - U*R*R/2;

    gi::symbol x("x");
    gi::ex qm, Ee, Ae;

    Ee = ricpad::conversions::mp_to_ex(E);
    Ae = ricpad::conversions::mp_to_ex(A);

    qm = (Ee*gi::power(x,2)-Ae)/(1-gi::power(x,2))
        -gi::power(m,2)/(gi::power(1-gi::power(x,2),2));
    qm = gi::series_to_poly(qm.series(x, 2*N+1));

    for ( int i = 0; i <= N; i++ ) {
        ans.emplace_back(
                ricpad::conversions::ex_to_mp<T>(
                    gi::ex_to<gi::numeric>(
                        qm.coeff(x, 2*i).evalf()
                        )
                    )
                );
    }

    return ans;
}
}

// Returns a vector with coefficients f[-1] ... f[N]
//
// Vectors p and q contain only the coefficients multiplying positive 
// (including 0) powers of x. The negative ones are considered explicitly
// in f[-1] == coefs[0].
template <class T>
std::vector<T> coefsl(
        const int N,
        const T& U,
        const T& A, 
        const T& R, 
        const int m, 
        // p[0], p[1], p[2], ...
        const std::vector<T>& p_coefs,
        // q[0], q[1], q[2], ...
        const std::vector<T>& q_coefs
        )
{
    T E = R/2 - U*R*R/2;
    std::vector<T> coefs;
    coefs.reserve(N+2);

    const int ma = abs(m);
    T sum;

    coefs.push_back((R + (A-E)/2 + ma*ma/2+ma/2)/(ma+1));

    for ( int j = 0; j <= N; j++ ) {
        sum = 0;
        for ( int k = 0; k <= j; k++ ) 
            // coefs[0] <==> f[-1], coefs[1] <==> f[0], etc.
            sum += (coefs[k]-p_coefs[k])*coefs[j-k];
        sum += q_coefs[j] + ma*p_coefs[j+1]/2;
        sum = sum/(ma + j + 2);
        coefs.push_back(sum);
    }
    
    return coefs;
}

// Returns a vector with coefficients f[0] ... f[N]. Here s can be either 
// 0 or 1, denoting even (odd) eigenfunctions.

template <class T>
std::vector<T> coefsm(
        const int N,
        const int s,
        // p[0], p[1], p[2], ...
        const std::vector<T>& p_coefs,
        // q[0], q[1], q[2], ...
        const std::vector<T>& q_coefs
        )
{
    std::vector<T> coefs;
    coefs.reserve(N+1);

    T sum;

    for ( int j = 0; j <= N; j++ ) {
        sum = 0;
        for ( int k = 0; k <= j-1; k++ ) 
            sum += (coefs[k]-p_coefs[k])*coefs[j-k-1];
        sum += q_coefs[j] + s*p_coefs[j+1]/2;
        sum = sum/(2*j + 2*s + 1);
        coefs.push_back(sum);
    }
    
    return coefs;
}

int main(int argc, char * argv[]) {
    using ricpad::hankdet::hankdet;
    using ricpad::solver::NR_solve;

    mpc_complex::default_precision(1200);
    gi::Digits = 1200;
    cout.precision(30);

    mpc_complex U("-0.6"), R("2.0"), A("1.0");
    mpc_complex h("(1E-300,1E-300)");

    mpfr_float tol("1E-150");

    int D = 2, d = 2, m = 0, s = 0;

    std::function<mpc_complex(Eigen::Matrix<mpc_complex,3,1>&)> 
        fun_m, fun_l, fun_d; 

    fun_m = [&D, &d, &s, &m]
        ( Eigen::Matrix<mpc_complex,3,1>& param ) -> mpc_complex {
            mpc_complex &U = param[0];
            mpc_complex &A = param[1];
            mpc_complex &R = param[2];

            std::vector<mpc_complex> pmv, qmv, coefs; 
            pmv = coefficients::pm<mpc_complex>(2*D+d+1);
            qmv = coefficients::qm<mpc_complex>(2*D+d+1, U, A, R, m);

            coefs = coefsm<mpc_complex>(2*D+d-1, s, pmv, qmv);
            coefs.erase(coefs.begin(), coefs.begin()+d+1);
            mpc_complex ans = hankdet<mpc_complex>(D,coefs);

            return ans;
    };

    fun_l = [&D, &d, &m]
        ( Eigen::Matrix<mpc_complex,3,1>& param ) -> mpc_complex { 
            mpc_complex &U = param[0];
            mpc_complex &A = param[1];
            mpc_complex &R = param[2];

            std::vector<mpc_complex> plv, qlv, coefs;
            plv = coefficients::pl<mpc_complex>(2*D+d+2);
            qlv = coefficients::ql<mpc_complex>(2*D+d+2, U, A, R, m);

            coefs = coefsl<mpc_complex>(2*D+d-1, U, A, R, m, plv, qlv);
            coefs.erase(coefs.begin(), coefs.begin()+d+2);
            mpc_complex ans = hankdet<mpc_complex>(D,coefs);

            return ans;
    };
    
    fun_d = [&D, &d, &m, &s, &fun_m, &fun_l]
        ( Eigen::Matrix<mpc_complex,3,1>& param ) -> mpc_complex {
            using ricpad::differentiate::differentiate;

            mpc_complex &U = param[0];
            mpc_complex &A = param[1];
            mpc_complex &R = param[2];

            mpc_complex h("(1E-600,1E-600)");

            return (
                differentiate<mpc_complex, 3>(fun_m, param, 1, h) *
                  differentiate<mpc_complex, 3>(fun_l, param, 2, h) -
                differentiate<mpc_complex, 3>(fun_m, param, 2, h) * 
                  differentiate<mpc_complex, 3>(fun_l, param, 1, h) 
                  );
        };

    std::vector<decltype(fun_m)> F;

    F.push_back(fun_m);
    F.push_back(fun_l);
    F.push_back(fun_d);

    Eigen::Matrix<mpc_complex,3,1> params;

    params(0) = U;
    params(1) = A;
    params(2) = R;

    std::ofstream U_stream, A_stream, R_stream;
    
    for ( ; D<=80; D++ ) {
        params = NR_solve<mpc_complex, mpfr_float, 3>(F, params, tol, h);
        U = params(0);
        A = params(1);
        R = params(2);

        cout.precision(30);
        cout << " D = " << std::setw(2) << D << " ";
        cout << params(0).real() << " " << params(1).real() << " " << 
            params(2).real() << endl;

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
    }

    return 0;
}
