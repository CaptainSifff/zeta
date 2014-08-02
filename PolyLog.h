#ifndef POLYLOG_H
#define POLY_LOG_H

/**
 * A function to reliably compare two floating point numbers
 * @param a
 * @param b
 * @return returns true if a and b are equal to zero or differ only by max(a,b)*eps
 * */
template <typename FPType> 
bool fpequal(const FPType& a, const FPType& b)
{
  bool retval = true;
    if ((a != FPType(0)) || (b != FPType(0)))//looks mean, but is necessary that the next line has sense.
      retval = (std::fabs(a - b) < 5.0*std::max(std::fabs(a), std::fabs(b))*std::numeric_limits<FPType>::epsilon());
  return retval;
}

template <typename FPType>
inline std::complex<FPType> PolyLog_Exp(const unsigned int s, std::complex<FPType> w)
{//positive integer s
}

template <typename FPType>
inline std::complex<FPType> PolyLog_Exp_neg(const unsigned int s, std::complex<FPType> w)
{//negative integer s
}

template <typename FPType>
inline std::complex<FPType> PolyLog_Exp_neg(const FPType s, std::complex<FPType> w)
{//negative s
}

template <typename FPType>
inline std::complex<FPType> PolyLog_Exp_pos(const FPType s, std::complex<FPType> w)
{//positive s
}

template <typename FPType>
inline std::complex<FPType> PolyLog_Exp(const FPType s, std::complex<FPType> w)
    {
      /*reduce the imaginary part to the range where the series converges quickly*/
      if (fpequal<FPType>(std::rint(s), s))
      {//capture the cases of positive integer index
	int nu = static_cast<int> (lrint(s));
	if(0 == nu)
	{
	  std::complex<FPType> t = std::exp(w);
	  return t/(1.0 - t);
	}
	else if (1 == nu)
	  return -std::log(1.0 - std::exp(w));
	else if (nu > 1)
	{//FIXME: check for real or non-real argument. asymptotic expansions
	  while (w.imag() > M_PI) w = std::complex<FPType>(w.real(), w.imag() - 2.0*M_PI);//branch cuts??
	  while (w.imag() <= -M_PI) w = std::complex<FPType>(w.real(), w.imag() + 2.0*M_PI);
	  if (fpequal(arg(w), 0.0) || fpequal(arg(w), 2.0*M_PI))
	    return mytr1::__detail::__riemann_zeta(s);
	  else
	    return PolyLog_Exp(static_cast<uint>(nu) , w);
	}
	else
	  return PolyLog_Exp_neg(s, w);
      }
      else
      {//FIXME: check for real or non-real argument. asymptotic expansions
	  while (w.imag() > M_PI) w = std::complex<FPType>(w.real(), w.imag() - 2.0*M_PI);//branch cuts??
	  while (w.imag() <= -M_PI) w = std::complex<FPType>(w.real(), w.imag() + 2.0*M_PI);
	if (s < 0)
	  return PolyLog_Exp_neg(s, w);
	else
	{
	  if (fpequal(arg(w), 0.0) || fpequal(arg(w), 2.0*M_PI))
	    return mytr1::__detail::__riemann_zeta(s);
	  else
	    return PolyLog_Exp_pos(s, w);
	}
      }
    }
    
template <typename FPType>
inline std::complex<FPType> PolyLog(const FPType s, std::complex<FPType> w)
{
}

#endif