#include <vector>
#include <boost/math/policies/error_handling.hpp>

#ifndef RICPAD_HANKDET
#define RICPAD_HANKDET

namespace ricpad::hankdet {
// Return the hankel determinant with D and d, for Problem problem and 
// parameter s, evaluated at point E0
template <class T>
T hankdet(
        const int D, 
        // Coefficients f[d+2]...f[2*D+d-2]. This vector is destroyed.
        std::vector<T>& coefs
        ) {
    using std::vector;

    vector<T> coefsm1, coefsm2;
    static const char* function = "ricpad::hankdet::hankdet<%1%>";

    // Check if we have enough coefficients
    if ( coefs.size() < 2*D-1 ) {
        return boost::math::policies::raise_evaluation_error(
            function,
            "Input vector coefs should contain at least %1% elements. ", 
            2*D-1,
            boost::math::policies::policy<>());
    }

    if ( D == 0 ) {
        return 1;
    } else {
        coefsm1 = std::vector<T>(2*D, T(1));

        if ( D == 1 ) return coefs[0];

        for ( int j = 2; j <= D; j++ ) {
            coefsm2 = (vector<T>&&)(coefsm1);
            coefsm1 = (vector<T>&&)(coefs);

            for ( int k = 0; k <= 2*(D-j); k++ ) {
                coefs.emplace_back(
                        (coefsm1[k]*coefsm1[k+2] - 
                        coefsm1[k+1]*coefsm1[k+1]) /
                        coefsm2[k+2]
                        );
            }
        }
    }

    return coefs[0];
}

} // namespace
#endif
