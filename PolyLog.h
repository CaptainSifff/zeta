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
      std::cout<<"Integer Series for positive s"<<std::endl;
      std::complex<FPType> res = mytr1::__detail::__riemann_zeta(static_cast<FPType>(s));//optimization possibility: s are positive integers
      std::complex<FPType> wpower = w;
      FPType fac = 1.0;
      FPType harmonicN = 1.0;//HarmonicNumber_1
      for (uint k = 1; k <= s-2; ++k)
      {
	res += wpower*fac*mytr1::__detail::__riemann_zeta(static_cast<FPType>(s - k));
	wpower *= w;
	FPType temp = 1.0/(1.0 + k);
	fac *= temp;
	harmonicN += temp;
      }
      //harmonicN now contains H_{s-1}
      //fac should be 1/(n-1)!
      res += (harmonicN - std::log(-w))*wpower*fac;
      wpower *= w;
      fac /= s;
      res -= wpower*fac/2.0;
      wpower *= w;
      //now comes the remainder of the series.
      const FPType tp = 2.0 * M_PI;
      const std::complex<FPType> pref = wpower/M_PI/tp;
      const unsigned int maxit = 200;
      unsigned int j = 1;
      bool terminate = false;
      fac /= (s+1.0);//(1/(n+1)!)
      res -= M_PI*M_PI/6.0*fac * pref; //subtract the zeroth order term.
      //remainder of series
      fac *= 3.0*2.0/(s + 2.0)/(s+3.0);
      std::complex<FPType> upfac = -(w/tp)*(w/tp);
      std::complex<FPType> w2 = upfac;
      while (!terminate)//assume uniform convergence
      {
        FPType rzarg = static_cast<FPType>(2*j+2);
        FPType rz = mytr1::__detail::__riemann_zeta(rzarg);
//        std::cout<<rz<<" "<<fac<<" "<<w2<<std::endl;
	std::complex<FPType> nextterm = (rz*fac)*w2;
	w2 *= upfac;
	fac *= rzarg/(rzarg + s) * (rzarg+1.0)/(rzarg + s + 1.0);
	++j;
	terminate = (fpequal( std::abs(res - pref*nextterm), std::abs(res) ) || (j > maxit));
	res -= pref * nextterm;
      }
      std::cout<<"Iterations in Integer Series: "<<j<<std::endl;
      return res;
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

/** This function implements the asymptotic series for the PolyLog.
 * It is given by 2 \sum \limits_{k=0}^\infty \zeta(2k) w^{s-2k}/Gamma(s-2k+1) -i \pi w^(s-1)/Gamma(s)
 * for Re(w) >> 1
 * Don't check this against Mathematica 8. 
 * For real x the imaginary part of the PolyLog is given by Im(Li_s(e^u)) = - \pi u^{s-1}/Gamma(s)
 * Check this relation for any benchmark that you use.
 * @param s the index s
 * @param w
 */
template <typename FPType>
inline std::complex<FPType> PolyLog_Exp_asym(const FPType s, std::complex<FPType> w)
{//asymptotic expansion
  std::cout<<"asymptotic expansions"<<std::endl;
  std::complex<FPType> wgamma = std::pow(w, s-1.0)/std::tgamma(s);/*wgamma = w^(s-1)/Gamma(s)*/
  std::complex<FPType> res = std::complex<FPType>(0.0, -M_PI)* wgamma;
  wgamma *= w/s;/*wgamma = w^s / Gamma(s+1)*/
  constexpr uint maxiter = 100;
  bool terminate = false;
  std::complex<FPType> oldterm = -0.5*wgamma; /*zeta(0) * w^s / Gamma(s+1)*/
  res += 2.0 * oldterm;
  std::complex<FPType> newterm;
  std::complex<FPType> wq = 1.0/(w*w);
  uint k = 1;
  while (!terminate)
  {
    wgamma *= wq * (s + 1.0 - 2*k) * (s + 2.0 - 2*k);
    newterm = mytr1::__detail::__riemann_zeta(static_cast<FPType> (2*k) ) * wgamma;
//    std::cout<<k<<" "<<newterm<<" "<< std::endl;
    if(std::abs(newterm) > std::abs(oldterm)) terminate = true;//termination due to failure of asymptotic expansion
    if(fpequal(std::abs(res + 2.0* newterm), std::abs(res))) terminate = true; // precision goal reached.
    if(k > maxiter) terminate = true;//stop the iteration somewhen
    if(!terminate)
    {
      res += 2.0*newterm;
      oldterm = newterm;
      ++k;
    }
  }
  std::cout<<"Iterations: "<<k<<std::endl;
  return res;
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
	{//FIXME: check for real or non-real argument.
	  /*The reductions of the imaginary part yield the same the same results as Mathematica then.
	   * Necessary to improve the speed of convergence
	   */
	  if(real(w) < 15.0)//arbitrary transition point...
	  {
	    while (w.imag() > M_PI) w = std::complex<FPType>(w.real(), w.imag() - 2.0*M_PI);
	    while (w.imag() <= -M_PI) w = std::complex<FPType>(w.real(), w.imag() + 2.0*M_PI);
	  if (w == 0.0)//catch the case of the PolyLog evaluated evaluated at e^0
	    return mytr1::__detail::__riemann_zeta(s);
	  else
	    return PolyLog_Exp_pos(static_cast<uint>(nu) , w);
	  }
	  else
	  {
	    //wikipedia says that this is required for Wood's formula
	    while (w.imag() > 0) w = std::complex<FPType>(w.real(), w.imag() - 2.0*M_PI);
	    while (w.imag() <= -2.0*M_PI) w = std::complex<FPType>(w.real(), w.imag() + 2.0*M_PI);
	    return PolyLog_Exp_asym(s,w);//FIXME: the series should terminate after a finite number of terms.
	  }
	}
	else
	  return PolyLog_Exp_neg(s, w);//no asymptotic expansion available...
      }
      else
      {//FIXME: check for real or non-real argument.
	if(real(w) < 15.0)//arbitrary transition point
	{
	  /*The reductions of the imaginary part yield the same the same results as Mathematica then.
	   * Necessary to improve the speed of convergence
	   */
	  while (w.imag() > M_PI) w = std::complex<FPType>(w.real(), w.imag() - 2.0*M_PI);//branch cuts??
	  while (w.imag() <= -M_PI) w = std::complex<FPType>(w.real(), w.imag() + 2.0*M_PI);
	  if (s < 0)
	    return PolyLog_Exp_neg(s, w);
	  else
	    {
	  if (w == 0.0)//catch the case of the PolyLog evaluated evaluated at e^0
	    return mytr1::__detail::__riemann_zeta(s);
	  else
	      return PolyLog_Exp_pos(s, w);
	}
	}
	else
	{
	  //wikipedia says that this is required for Wood's formula
	  while (w.imag() > 0) w = std::complex<FPType>(w.real(), w.imag() - 2.0*M_PI);
	  while (w.imag() <= -2.0*M_PI) w = std::complex<FPType>(w.real(), w.imag() + 2.0*M_PI);
	  return PolyLog_Exp_asym(s,w);
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