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
{   //positive integer s
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

/** This function catches the cases of negative real index s.
 * Theoretical convergence is present for |w| < 2*pi.
 * @param s the index s
 * @param w
 */
template <typename FPType>
inline std::complex<FPType> PolyLog_Exp_neg(const FPType s, std::complex<FPType> w)
{   //basic general loop, but s is a negative quantity here
    //TODO: optimize/fix for the case of negative Integers
    std::cout<<"Negative real s"<<std::endl;
    std::complex<FPType> res = std::tgamma(1-s)*std::pow(-w, s-1);
    constexpr FPType tp = 2.0 * M_PI;
    const std::complex<FPType> wup = w/tp;
    std::complex<FPType> w2 = wup;
    std::complex<FPType> pref = std::pow(tp, s)/M_PI;
    std::complex<FPType> gam = std::tgamma(1.0-s); //here we factor up the ratio of Gamma(1 - s + k)/k! . This ratio should be well behaved even for large k

    FPType sp, cp;
    sincos(M_PI/2.0 * s, &sp, &cp);
    //Here we add the expression that would result from ignoring the zeta function in the series.
    std::complex<FPType> expis(cp, sp);
    std::complex<FPType> p = tp - std::complex<FPType>(0.0, 1.0) * w;
    std::complex<FPType> q = tp + std::complex<FPType>(0.0, 1.0) * w;
    res += std::complex<FPType>(0.0, 1.0) * gam * (conj(expis) * std::pow(p, s-1.0) - expis *std::pow(q, s-1.0));//this can be optimized for real values of w
    //The above expression is the result of sum_k Gamma(1+k-s) /k! * sin(pi /2* (s-k)) * (w/2/pi)^k
    //Therefore we only need to sample values of zeta(n) on the real axis that really differ from one
    res += pref * sp * gam * (mytr1::__detail::__riemann_zeta(1-s) - 1.0);
    const unsigned int maxit = 200;
    unsigned int j = 1;
    bool terminate = false;
    gam*= (1.0 - s);
    while (!terminate)//assume uniform convergence
    {
        FPType rzarg = 1 + j - s;
        FPType rz = (mytr1::__detail::__riemann_zeta(rzarg) - 1.0);//only the difference to one is needed. FIXME: this expression underflows for rzarg > 50
        FPType sine;
        if(j & 1)//save the repeated recalculation of the sines
        {   /*odd*/
            sine = cp;
            if ( !((j-1)/ 2 & 1) )
                sine = -sine;
        }
        else
        {   /*even*/
            sine = sp;
            if((j/2) & 1)
                sine = -sine;
        }
        std::complex<FPType> nextterm =  w2 * gam * (sine * rz);
//	std::cout<<j<<" "<<nextterm<<" "<<rz<<" "<<std::endl;
        w2 *= wup;
        ++j;
        gam  *= rzarg/(j);//equal to 1/(j+1) since we have incremented j in the line above
        terminate = (fpequal( std::abs(res + pref*nextterm), std::abs(res) ) || (j > maxit)) && !fpequal(std::abs(rz), 0.0);/*this last check is necessary at integer s*/
        res += pref*nextterm;
    }
    std::cout<<"Iterations in PolyLogExp_neg: "<<j<<std::endl;
    return res;
}

/** This function catches the cases of negative integer index s which are multiples of four. In that case the sine occuring in the expansion 
 * occasionally takes on the value zero. We use that to provide an optimized series for p = -4n:
 * Li_p(e^w) = Gamma(1-p) * (-w)^{p-1} - A_p(w) - B_p(w)
 * with
 * A_p(w) = -2 / \sqrt{2 \pi} (-p)! / (2 \pi)^(-p/2) (1 + w^2/(4 pi^2))^{-1/2 + p/2} cos((1 - p) ArcTan(2 pi/ w))
 * and 
 * B_p(w) = - (2 pi)^p / pi * \sum \limits_{k = 0}^\infty \Gamma(2 + 2k - p)/ (2k+1)! (-1)^k (w/2/\pi)^(2k+1) (Zeta(2 + 2k - p) - 1.0)
 * suitable for |w| < 2 pi
 * @param n the index n = 4k
 * @param w
 */
template <typename FPType>
inline std::complex<FPType> PolyLog_Exp_neg_four(const int n, std::complex<FPType> w)
{
  std::cout<<"Negative integer s = -4k"<<std::endl;
  std::complex<FPType> res = std::tgamma(1-n)*std::pow(-w, n-1);
  constexpr FPType tp = 2.0 * M_PI;
  std::complex<FPType> wup = w/tp;
  std::complex<FPType> wq = wup*wup;
  FPType pref = std::pow(tp, n)/M_PI;
  res += std::tgamma(1-n)* pref * std::pow(1.0 + wq, -0.5 + n/2) * 
//  std::sin( static_cast<FPType>(1-n) * (M_PI/2-std::atan(1.0/wup)));//subtract  the expression A_p(w)
  //using the sine here yields better results...
  std::cos( static_cast<FPType>(1-n) * std::atan(1.0/wup));//subtract  the expression A_p(w)
  uint k = 0;
  bool terminate = false;
  uint maxit = 300;
  FPType gam = std::tgamma(2-n);
  while(!terminate)
  {
    std::complex<FPType> newterm = ( gam * (mytr1::__detail::__riemann_zeta(static_cast<FPType>(2 + 2*k - n)) - 1.0)) * wup;
    gam *= - static_cast<FPType>(2 + 2*k -n + 1) / (2*k + 1 + 2) * static_cast<FPType>(2 + 2*k -n) / (2 * k + 1 + 1);
    wup *= wq;
    terminate = (fpequal( std::abs(res - pref*newterm), std::abs(res) ) || (k > maxit));
    res -= pref*newterm;
    ++k;
  }
  std::cout<<"Iterations in the series for s = -4n : "<<k<<std::endl;
  return res;
}

/** This function catches the cases of negative integer index s
 * @param s the index s
 * @param w
 */
template <typename FPType>
inline std::complex<FPType> PolyLog_Exp_neg(const int s, std::complex<FPType> w)
{//negative integer s
  if ((-s)% 4 == 0 )//Divisible by four. Then the sine in the resulting series is occasionaly zero...
    return PolyLog_Exp_neg_four(s, w);
  else
    return PolyLog_Exp_neg(static_cast<FPType>(s), w);
}

/** This function catches the cases of positive real index s
 * @param s the index s
 * @param w
 */
template <typename FPType>
inline std::complex<FPType> PolyLog_Exp_pos(const FPType s, std::complex<FPType> w)
{   //positive s
    std::cout<<"Series for real positive s"<<std::endl;
    std::complex<FPType> res = mytr1::__detail::__riemann_zeta(s);
    std::complex<FPType> wpower = w;
    FPType fac = 1.0;
    uint m = static_cast<uint>(std::floor(s));
    for (uint k = 1; k <= m; ++k)
    {
        res += wpower*fac*mytr1::__detail::__riemann_zeta(static_cast<FPType>(s - k));
        wpower *= w;
        FPType temp = 1.0/(1.0 + k);
        fac *= temp;
    }
    //fac should be 1/(m+1)!
    res += std::tgamma(1-s)*std::pow(-w, s-1);
    const FPType tp = 2.0 * M_PI;
    const FPType pref = 2.0 * std::pow(tp, s-1);
    //now comes the remainder of the series
    const unsigned int maxit = 100;
    unsigned int j = 0;
    bool terminate = false;
    std::complex<FPType> wup = w/tp;
    std::complex<FPType> w2 = std::pow(wup, m+1);
    std::complex<FPType> gam = std::tgamma(2.0-s+m)*fac; //here we factor up the ratio of Gamma(1 - s + k)/k! . This ratio should be well behaved even for large k
    FPType sp, cp;
    sincos(M_PI/2.0 * s, &sp, &cp);
    while (!terminate)//assume uniform convergence
    {   //FIXME: optimize.
        int idx = m + 1 + j;
        FPType zetaarg = 1 + idx - s;
        FPType sine;
        if(idx & 1)//save the reperated calculation of the sines
        {   /*odd*/
            sine = cp;
            if ( !((idx-1)/ 2 & 1) )
                sine = -sine;
        }
        else
        {   /*even*/
            sine = sp;
            if((idx/2) & 1)
                sine = -sine;
        }
        std::complex<FPType> nextterm = (mytr1::__detail::__riemann_zeta(zetaarg) * sine * gam) * w2;
//	std::cout<<j<<" "<<nextterm<<" used Gamma = "<<gam<<std::endl;
        w2 *= wup;
        gam *= zetaarg/(1.0 + idx);
        ++j;
        terminate = (fpequal( std::abs(res + pref*nextterm), std::abs(res) ) || (j > maxit));
        res += pref * nextterm;
    }
    std::cout<<"Iterations in PolyLogExp_pos: "<<j<<std::endl;
    return res;
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
{   //asymptotic expansion
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

/**
 * Theoretical convergence for Re(w) < 0. Seems to beat the other expansions for Re(w) < -pi/2 - pi/5
 */
template <typename FPType>
std::complex<FPType> Poly_log_exp_negative_real_part(FPType s, std::complex<FPType> w)
{
  std::cout<<"negative real part series (exponential)"<<std::endl;
  std::complex<FPType> ew = std::exp(w);
  std::complex<FPType> up = ew;
  std::complex<FPType> res = ew;
  uint maxiter = 500;
  bool terminate = false;
  uint k = 2;
  while(!terminate)
  {
    ew *= up;
    std::complex<FPType> newterm = std::pow(k, -s) * ew;
    terminate = (fpequal(std::abs(res + newterm), std::abs(res))) || (k > maxiter);
    res += newterm;
    ++k;
  }
  std::cout<<"iterations in Poly_log_exp_negative_real_part: "<<k<<std::endl;
  return res;
}

/** This is the Frontend function which calculates Li_s( e^w )
 * @param s the index s
 * @param w complex w
 */
template <typename FPType>
inline std::complex<FPType> PolyLog_Exp(const FPType s, std::complex<FPType> w)
{
    if(fpequal(real(w), 0.0) && (fpequal(imag(w), 0.0) || fpequal(imag(w), 2.0*M_PI)))//catch the case of the PolyLog evaluated at e^0 FIXME: higher multiples of 2 pi
    {
        if (s > 1.0)
            return mytr1::__detail::__riemann_zeta(s);
        else
            return std::numeric_limits<FPType>::infinity();
    }
    if (fpequal<FPType>(std::rint(s), s))
    {   //capture the cases of positive integer index
        int nu = static_cast<int> (lrint(s));
        if(0 == nu)
        {
            std::complex<FPType> t = std::exp(w);
            return t/(1.0 - t);
        }
        else if (1 == nu)
            return -std::log(1.0 - std::exp(w));
        else if (nu > 1)
        {   //FIXME: check for real or non-real argument.
            if(real(w) < 15.0)//arbitrary transition point...
            {
	       /*The reductions of the imaginary part yield the same results as Mathematica.
               * Necessary to improve the speed of convergence
               */
                while (w.imag() > M_PI) w = std::complex<FPType>(w.real(), w.imag() - 2.0*M_PI);
                while (w.imag() <= -M_PI) w = std::complex<FPType>(w.real(), w.imag() + 2.0*M_PI);
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
	{
	  if (( ((-nu) & 1) == 0) && fpequal(real(w), 0.0))
	  {
	    FPType ip = imag(w);//get imaginary part
	    while (ip <= -M_PI) ip += 2.0*M_PI;
	    while (ip > M_PI) ip -= 2.0*M_PI;
	    if(fpequal(ip, M_PI))
	      return 0.0;//Li_{-n}(-1) + (-1)^n Li_{-n}(1/-1) = 0 
	    else
	    {
	      return PolyLog_Exp_neg(nu, std::complex<FPType>(w.real(), w.imag()/*why does using ip not work here??*/ ));//no asymptotic expansion available... check the reduction
	    }
	  }
	  else
//	  if ((-nu) & 1)
	    return PolyLog_Exp_neg(nu, w);//no asymptotic expansion available...
//	    else return 0.0;//Li_{-n}(1) + (-1)^n Li_{-n}(1) = 0 
	}
    }
    else
    {   //FIXME: check for real or non-real argument.
        if(real(w) < 15.0)//arbitrary transition point
        {
            /*The reductions of the imaginary part yield the same results as Mathematica then.
             * Necessary to improve the speed of convergence
             */
            while (w.imag() > M_PI) w = std::complex<FPType>(w.real(), w.imag() - 2.0*M_PI);//branch cuts??
            while (w.imag() <= -M_PI) w = std::complex<FPType>(w.real(), w.imag() + 2.0*M_PI);
            if (s < 0)
                return PolyLog_Exp_neg(s, w);
            else
                return PolyLog_Exp_pos(s, w);
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
    if(fpequal(real(w), 0.0) && (fpequal(imag(w), 0.0) || fpequal(imag(w), 2.0*M_PI)))//catch the case of the PolyLog evaluated evaluated at e^0
    {
        if (s > 1.0)
            return mytr1::__detail::__riemann_zeta(s);
        else
            return std::numeric_limits<FPType>::infinity();
    }
    if (fpequal<FPType>(std::rint(s), s))
    {   //capture the cases of positive integer index
        int nu = static_cast<int> (lrint(s));
        if(0 == nu)
        {
            FPType t = std::exp(w);
            return t/(1.0 - t);
        }
        else if (1 == nu)
            return -std::log(1.0 - std::exp(w));
        else if (nu > 1)
        {
            if(real(w) < 15.0)//arbitrary transition point...
                return PolyLog_Exp_pos(static_cast<uint>(nu) , w);
            else
                return PolyLog_Exp_asym(s,w);//FIXME: the series should terminate after a finite number of terms.
        }
        else//e^w with real w is always a positive quantity
	    return PolyLog_Exp_neg(s, w);//no asymptotic expansion available...
    }
    else
    {
        if(real(w) < 15.0)//arbitrary transition point
        {
            if (s < 0)
                return PolyLog_Exp_neg(s, w);
            else
                return PolyLog_Exp_pos(s, w);
        }
        else
            return PolyLog_Exp_asym(s,w);
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
    {   //use the square formula to access negative values.
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

/** A function to implement Dirichlet's Eta function
 * @param w A w
*/
template <typename FPType>
inline std::complex<FPType> Dirichlet_eta(std::complex<FPType> w)
{
    if(w.imag() == 0.0)
        return -PolyLog(w.real(), -1.0);
    else
    {   std::cout<<"Domain not (yet) supported!!"<<std::endl;
        std::__throw_domain_error(__N("Bad argument to Dirichlet Eta."));
    }
}

/** A function to implement Dirichlet's beta function
 * @param w A w
*/
template <typename FPType>
inline FPType Dirichlet_beta(std::complex<FPType> w)
{
    if(w.imag() == 0.0)
        return imag(PolyLog(w.real(), std::complex<FPType>(0.0, 1.0)));
    else
    {   std::cout<<"domain not (yet) supported!!"<<std::endl;
        std::__throw_domain_error(__N("Bad argument to Dirichlet Eta."));
    }
}

/** A function to implement Claussen's series Sl.
 * Notation and connection to polylog from wikipedia
 * @param w A  w
 * FIXME: Check the restriction to positive integers m
*/
template <typename FPType>
inline FPType Claussen_Sl(uint m, std::complex<FPType> w)
{
    std::complex<FPType> ple = PolyLog_Exp(m, std::complex<FPType>(0.0, 1.0) * w);
    if (m & 1)
        return imag(ple);
    else
        return real(ple);
}

/** A function to implement Claussen's series Cl
 * @param w A  w
 * FIXME: Check the restriction to positive integers m
*/
template <typename FPType>
inline FPType Claussen_Cl(uint m, std::complex<FPType> w)
{
    std::complex<FPType> ple = PolyLog_Exp(m, std::complex<FPType>(0.0, 1.0) * w);
    if (m & 1)
        return real(ple);
    else
        return imag(ple);
}
#endif
