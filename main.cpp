#include <iostream>
#include <chrono>
//#include <tr1/cmath> //for the original implementation
#include <complex>
#include <fstream>
#include <cmath>
#include "riemann_zeta.tcc"
using namespace std;
int main(int argc, char **argv) {
  typedef double FPType;
  uint n = 5000;
  /*this part of code was for performance testing. the old implementation takes about 2.8s on my core2 and the new one 0.8s*/
   cout.precision(14);
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
std::cout<<mytr1::__detail::PolyLog_Exp(2.5, std::complex<FPType>(15,1.0))<<std::endl;
// for(uint k = 0; k < 32; ++k)
// {
//   std::cout<<"=======  "<<k<<"  =========="<<std::endl;
//   std::cout<<mytr1::__detail::PolyLog_Exp(4.0, std::complex<FPType>(0.0, 2.0*M_PI * k/32.0))<<std::endl;
// std::cout<<std::scientific<<mytr1::__detail::PolyLog_Exp(-4.0, std::complex<FPType>(0.0, 2.0*M_PI * k/32.0))<<std::endl;
// std::cout<<mytr1::__detail::PolyLog_Exp(2.6, std::complex<FPType>(0.0, 2.0*M_PI * k/32.0))<<std::endl;
// std::cout<<mytr1::__detail::PolyLog_Exp(-2.6, std::complex<FPType>(0.0, 2.0*M_PI * k/32.0))<<std::endl;
// }
//std::cout<<mytr1::__detail::PolyLog_Exp(2.6, std::complex<FPType>(M_PI, M_PI))<<std::endl;
/*for(uint k = 0; k < 10; ++k)
{
std::cout<<mytr1::__detail::PolyLog_Exp(-4.0, std::complex<FPType>(-M_PI/2 - M_PI*4/20 , 0))<<std::endl;
std::cout<<mytr1::__detail::Poly_log_exp_negative_real_part(-4.0, std::complex<FPType>(-M_PI/2 - M_PI*4/20, 0))<<std::endl;
}*/
//std::cout<<mytr1::__detail::PolyLog_Exp_neg(-50.5, std::complex<FPType>(1.0, 1.0))<<std::endl;
//std::cout<<mytr1::__detail::PolyLog_Exp_neg_four(-4, std::complex<FPType>(1.0, 1.0))<<std::endl;
//std::cout<<mytr1::__detail::PolyLog_Exp_pos(2.3, std::complex<FPType>(1.0, 1.0))<<std::endl;
std::cout<<mytr1::__detail::PolyLog_Exp_asym(60.4, std::complex<FPType>(30.0, 0.0))<<std::endl;//Don't trust Mathematica for small s
// auto l = 2;
// auto p = std::atan(l);
// FPType alpha[] = {0.5, 1, 1.5, 4};
// std::ofstream data("el20.txt");
// for(uint a = 0; a < sizeof(alpha)/sizeof(FPType); ++a)
// {
//   for(int s = -1 ; s <= 1; s += 2)
//   {
// for(auto k = -M_PI; k < M_PI; k+= 0.002)
//   data<<k<<" "<<std::sqrt(1.0 + l*l)*real(std::exp(std::complex<FPType>(0, -s*p)) / (std::exp(std::complex<FPType>(0, k)) - std::exp(-alpha[a])))<<std::endl;
// data<<"&"<<std::endl;
//   }
// }
    return 0;
}
