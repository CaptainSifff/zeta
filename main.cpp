#include <iostream>
#include <chrono>
//#include <tr1/cmath> //for the original implementation
#include <complex>
#include <cmath>
#include "riemann_zeta.tcc"
using namespace std;
int main(int argc, char **argv) {
  typedef double FPType;
  uint n = 5000;
  /*this part of code was for performance testing. the old implementation takes about 2.8s on my core2 and the new one 0.8s*/
//   cout.precision(14);
//     for(uint i = 0; i < n; ++i)
//     {
//       FPType x = 10.0* static_cast<FPType>(i)/n + 1.1;
// //      std::cout<<std::scientific<<x<<" "<<
//       mytr1::__detail::__riemann_zeta(x)
// //      std::tr1::__detail::__riemann_zeta(x)
//       ;//between 1 and 10 riemann_zeta_glob is called
// //      <<std::endl;
//     }

/*something that didn't work in the original implementation*/
  std::cout<<mytr1::__detail::__hurwitz_zeta_glob(5.1, 0.5)<<std::endl;
  std::cout<<mytr1::__detail::hurwitz_zeta(5.1, std::complex<FPType>(0.5))<<std::endl;
// std::cout<<mytr1::__detail::PolyLog_Exp(2.5, std::complex<FPType>(15,1.0))<<std::endl;
// for(uint k = 0; k < 32; ++k)
// {
//   std::cout<<"=======  "<<k<<"  =========="<<std::endl;
//   std::cout<<mytr1::__detail::PolyLog_Exp(4.0, std::complex<FPType>(0.0, 2.0*M_PI * k/32.0))<<std::endl;
// std::cout<<mytr1::__detail::PolyLog_Exp(-4.0, std::complex<FPType>(0.0, 2.0*M_PI * k/32.0))<<std::endl;
// std::cout<<mytr1::__detail::PolyLog_Exp(2.6, std::complex<FPType>(0.0, 2.0*M_PI * k/32.0))<<std::endl;
// std::cout<<mytr1::__detail::PolyLog_Exp(-2.6, std::complex<FPType>(0.0, 2.0*M_PI * k/32.0))<<std::endl;
// }
    return 0;
}
