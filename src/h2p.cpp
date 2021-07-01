#include <iostream>
#include <vector>
#include <fstream>

#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/eigen.hpp>
#include <ginac/ginac.h>
#include <Eigen/Dense>

#include <ricpad/hankdet.hpp>
#include <ricpad/solver.hpp>
#include <ricpad/conversions.hpp>

using boost::multiprecision::mpfr_float;

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

/*
*/
