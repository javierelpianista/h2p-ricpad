#ifndef RICPAD_DIFFERENTIATE
#define RICPAD_DIFFERENTIATE

#include <Eigen/Dense>

namespace ricpad {
namespace differentiate {

template <typename T, int N>
T differentiate(
        // A function which takes an Eigen Matrix
        const std::function<T(Eigen::Matrix<T,N,1>&)>& f,
        // The variables contained in an Eigen matrix
        const Eigen::Matrix<T, N, 1>& x,
        // With respect to which variable we differentiate
        const int k,
        // Step size
        const T& h
      )
{
    Eigen::Matrix<T, N, 1> xp(x), xm(x);

    xp(k) += h;
    xm(k) -= h;

    T ans = f(xp) - f(xm);
    ans /= (2*h);

    return ans;
};

}; // namespace differentiate
}; // namespace ricpad
#endif
