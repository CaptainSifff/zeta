#ifndef HURWITZ_ZETA_H
#define HURWITZ_ZETA_H
#include "PolyLog.h"

template <typename FPType>
inline std::complex<FPType> hurwitz_zeta(const FPType s, std::complex<FPType> x)
    {
      if ( 
	((x.imag() >= 0) && ((x.real() >= 0.0) && (x.real() <  1.0) )) || 
	((x.imag() <  0) && ((x.real() >  0.0) && (x.real() <= 1.0) ))
      )
      {
	constexpr FPType tp = 2.0*M_PI;
	FPType t = 1.0-s;
	std::complex<FPType> lpe = PolyLog_Exp(t, std::complex<FPType>(0.0, tp)/*== 2 pi I */ * x );
	//FIXME: This prefactor is prone to overflow
	return std::tgamma(t)* std::pow(tp, -t)* ( std::exp(std::complex<FPType>(0.0, -M_PI/2.0 * t)) * lpe + std::exp(std::complex<FPType>(0.0, M_PI/2.0 * t)) * conj(lpe) );
      }
      else 
      {
	std::cout<<"domain not (yet) supported!!"<<std::endl;
	std::__throw_domain_error(__N("Bad argument to Hurwitz zeta"));
      }
    }
#endif