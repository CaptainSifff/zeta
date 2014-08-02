#include <iostream>
#include <chrono>
#include <cmath>
#include <complex>
#include "riemann_zeta.tcc"
using namespace std;
int main(int argc, char **argv) {
  typedef double FPType;
  uint n = 5000;
//   cout.precision(14);
//     for(uint i = 0; i < n; ++i)
//     {
//       FPType x = 40.0* static_cast<FPType>(i)/n + 1.1;
// //      std::cout<<std::scientific<<x<<" "<<
//       std::tr1::__detail::__riemann_zeta(x)
//       ;//between 1 and 10 riemann_zeta_glob is called
// //      <<std::endl;
//     }
  std::cout<<mytr1::__detail::__hurwitz_zeta_glob(5.0, 2.0)<<std::endl;
  std::cout<<mytr1::__detail::hurwitz_zeta(5.0, std::complex<FPType>(2.0))<<std::endl;
    return 0;
}
