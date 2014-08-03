#ifndef POLYLOG_H
#define POLY_LOG_H

/**
 * A function to reliably compare two floating point numbers.
 * @param a
 * @param b
 * @return returns true if a and b are equal to zero or differ only by max(a,b)* 5* eps
 * */
template <typename FPType> 
bool fpequal(const FPType& a, const FPType& b)
{
  bool retval = true;
    if ((a != FPType(0)) || (b != FPType(0)))//looks mean, but is necessary that the next line has sense.
      retval = (std::fabs(a - b) < 5.0*std::max(std::fabs(a), std::fabs(b))*std::numeric_limits<FPType>::epsilon());
  return retval;
}

/** This function catches the cases of positive integer index s
 * @param s the index s
 * @param w
 */
template <typename FPType>
inline std::complex<FPType> PolyLog_Exp_pos(const unsigned int s, std::complex<FPType> w)
{//positive integer s
}

/** This function catches the cases of negative integer index s
 * @param s the index s
 * @param w
 */
template <typename FPType>
inline std::complex<FPType> PolyLog_Exp_neg(const unsigned int s, std::complex<FPType> w)
{//negative integer s
}

/** This function catches the cases of negative real index s
 * @param s the index s
 * @param w
 */
template <typename FPType>
inline std::complex<FPType> PolyLog_Exp_neg(const FPType s, std::complex<FPType> w)
{//negative s
}

/** This function catches the cases of positive real index s
 * @param s the index s
 * @param w
 */
template <typename FPType>
inline std::complex<FPType> PolyLog_Exp_pos(const FPType s, std::complex<FPType> w)
{//positive s
}

/** This is the Frontend function which calculates Li_s( e^w )
 * @param s the index s
 * @param w complex w
 */
template <typename FPType>
inline std::complex<FPType> PolyLog_Exp(const FPType s, std::complex<FPType> w)
    {
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
	  /*The reductions of the imaginary part yield the same the same results as Mathematica then.
	   * Necessary to improve the speed of convergence
	   */
	  while (w.imag() > M_PI) w = std::complex<FPType>(w.real(), w.imag() - 2.0*M_PI);
	  while (w.imag() <= -M_PI) w = std::complex<FPType>(w.real(), w.imag() + 2.0*M_PI);
// 	  if (fpequal(arg(w), 0.0) || fpequal(arg(w), 2.0*M_PI))
// 	    return mytr1::__detail::__riemann_zeta(s);
// 	  else
	    return PolyLog_Exp_pos(static_cast<uint>(nu) , w);
	}
	else
	  return PolyLog_Exp_neg(s, w);
      }
      else
      {//FIXME: check for real or non-real argument. asymptotic expansions
	  /*The reductions of the imaginary part yield the same the same results as Mathematica then.
	   * Necessary to improve the speed of convergence
	   */
	while (w.imag() > M_PI) w = std::complex<FPType>(w.real(), w.imag() - 2.0*M_PI);//branch cuts??
	while (w.imag() <= -M_PI) w = std::complex<FPType>(w.real(), w.imag() + 2.0*M_PI);
	if (s < 0)
	  return PolyLog_Exp_neg(s, w);
	else
	{
// 	  if (fpequal(arg(w), 0.0) || fpequal(arg(w), 2.0*M_PI))
// 	    return mytr1::__detail::__riemann_zeta(s);
// 	  else
	    return PolyLog_Exp_pos(s, w);
	}
      }
    }

/** This is the Frontend function which calculates Li_s( e^w )
 * @param s the index s
 * @param w real w
 */
template <typename FPType>
inline std::complex<FPType> PolyLog_Exp(const FPType s, FPType w)
    {
      if (fpequal<FPType>(std::rint(s), s))
      {//capture the cases of positive integer index
	int nu = static_cast<int> (lrint(s));
	if(0 == nu)
	{
	  FPType t = std::exp(w);
	  return t/(1.0 - t);
	}
	else if (1 == nu)
	  return -std::log(1.0 - std::exp(w));
	else if (nu > 1)
	{//FIXME: check for real or non-real argument. asymptotic expansions
// 	  if (fpequal(arg(w), 0.0) || fpequal(arg(w), 2.0*M_PI))
// 	    return mytr1::__detail::__riemann_zeta(s);
// 	  else
	    return PolyLog_Exp(static_cast<uint>(nu) , w);
	}
	else
	  return PolyLog_Exp_neg(s, w);
      }
      else
      {//FIXME: check for real or non-real argument. asymptotic expansions
	if (s < 0)
	  return PolyLog_Exp_neg(s, w);
	else
	{
// 	  if (fpequal(arg(w), 0.0) || fpequal(arg(w), 2.0*M_PI))
// 	    return mytr1::__detail::__riemann_zeta(s);
// 	  else
	    return PolyLog_Exp_pos(s, w);
	}
      }
    }

/** A function to implement the PolyLog in those cases where we can calculate it.
 * @param s The index s
 * @param w A complex w
 */
template <typename FPType>
inline std::complex<FPType> PolyLog(const FPType s, std::complex<FPType> w)
{
  if(fpequal(arg(w), 0.0))
    return PolyLog(s, real(w));
  else
  {
    	std::cout<<"domain not (yet) supported!!"<<std::endl;
	std::__throw_domain_error(__N("Bad argument to PolyLog"));
  }
}

/** A function to implement the PolyLog for two real arguments.
 * @param s The index s
 * @param w A real w
*/
template <typename FPType>
inline std::complex<FPType> PolyLog(const FPType s, FPType w)
{
  if (fpequal(w, 0.0)) return 0.0;//According to Mathematica
  if (w < 0)
  {//use the square formula to access negative values.
    FPType wp = -w;
    FPType x = std::log(wp);
    return PolyLog_Exp(s, 2.0 * x) * std::pow(2.0, 1.0-s) - PolyLog_Exp(s, x);
  }
  else
  {
    FPType x = std::log(w);
    return PolyLog_Exp(s, w);
  }
}

#endif