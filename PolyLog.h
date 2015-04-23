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

/**
 * A function to calculate the values of zeta at even positive integers. For values smaller than thirty a table is used.
 * @param k 
 * @return zeta(k)
 */
template <typename FPType = double> 
inline FPType evenzeta(const uint& k)
{
  //the following constants were calculated with Mathematica 8
  constexpr FPType data[] = {
-0.50000000000000000000000000,
 1.6449340668482264364724152,
 1.0823232337111381915160037,
 1.0173430619844491397145179,
 1.0040773561979443393786852,
 1.0009945751278180853371460,
 1.0002460865533080482986380,
 1.0000612481350587048292585,
 1.0000152822594086518717326,
 1.0000038172932649998398565,
 1.0000009539620338727961132,
 1.0000002384505027277329900,
 1.0000000596081890512594796,
 1.0000000149015548283650412,
 1.0000000037253340247884571,   
  };
  constexpr auto maxk = 2 * sizeof(data)/sizeof(FPType);
  FPType retval;
  if (k < maxk)
    retval = data[k/2];
  else
    retval = mytr1::__detail::__riemann_zeta(static_cast<FPType>(k));
  return retval;
}

/** This function catches the cases of positive integer index s.
 * Li_s(e^w) = \sum_{k=0, k != s-1} \zeta(s-k) w^k/k! + (H_{s-1} - log(-w)) w^(s-1)/(s-1)!
 * The radius of convergence is |w| < 2 pi.
 * Note that this series involves a log(-x).
 * gcc and Mathematica differ in their implementation of \log(e^(i \pi)):
 * gcc: \log(e^(+- i * \pi)) = +- i \pi
 * whereas Mathematica doesn't preserve the sign in this case: \log(e^(+- i * \pi)) = +i \pi
 * @param s the index s
 * @param w
 */
template <typename FPType>
std::complex<FPType> PolyLog_Exp_pos(const unsigned int s, std::complex<FPType> w)
{   //positive integer s
//    std::cout<<"Integer Series for positive s - 1"<<std::endl;
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
        uint rzarg = 2*j+2;
//        FPType rz = mytr1::__detail::__riemann_zeta(rzarg);
	FPType rz = evenzeta<FPType>(rzarg);
//        std::cout<<rz<<" "<<fac<<" "<<w2<<std::endl;
        std::complex<FPType> nextterm = (rz*fac)*w2;
        w2 *= upfac;
        fac *= rzarg/static_cast<FPType>(rzarg + s) * (rzarg+1)/static_cast<FPType>(rzarg + s + 1);
        ++j;
        terminate = (fpequal( std::abs(res - pref*nextterm), std::abs(res) ) || (j > maxit));
        res -= pref * nextterm;
    }
    std::cout<<"Iterations in Integer Series: "<<j<<'\n';
    return res;
}

/** This function catches the cases of positive integer index s for real w.
 * This specialization is worthwhile to catch the differing behaviour of log(x).
 * Li_s(e^w) = \sum_{k=0, k != s-1} \zeta(s-k) w^k/k! + (H_{s-1} - log(-w)) w^(s-1)/(s-1)!
 * The radius of convergence is |w| < 2 pi.
 * Note that this series involves a log(-x).
 * The use of evenzeta yields a speedup of about 2.5
 * gcc and Mathematica differ in their implementation of \log(e^(i \pi)):
 * gcc: \log(e^(+- i * \pi)) = +- i \pi
 * whereas Mathematica doesn't preserve the sign in this case: \log(e^(+- i * \pi)) = +i \pi
 * @param s the index s
 * @param w
 */
template <typename FPType>
inline std::complex<FPType> PolyLog_Exp_pos(const unsigned int s, FPType w)
{   //positive integer s
//    std::cout<<"Integer Series for positive s - 2 "<<std::endl;
    FPType res = mytr1::__detail::__riemann_zeta(static_cast<FPType>(s));//optimization possibility: s are positive integers
    FPType wpower = w;
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
    std::complex<FPType> imagtemp = (harmonicN - std::log(std::complex<FPType>(-w, 0.0)))*wpower*fac;
    res += real(imagtemp);
//    res += (harmonicN - std::log(-w))*wpower*fac;
    wpower *= w;
    fac /= s;
    res -= wpower*fac/2.0;
    wpower *= w;
    //now comes the remainder of the series.
    const FPType tp = 2.0 * M_PI;
    const FPType pref = wpower/M_PI/tp;
    const unsigned int maxit = 200;
    unsigned int j = 1;
    bool terminate = false;
    fac /= (s+1.0);//(1/(n+1)!)
    res -= M_PI*M_PI/6.0*fac * pref; //subtract the zeroth order term.
    //remainder of series
    fac *= 3.0*2.0/(s + 2.0)/(s+3.0);
    FPType upfac = -(w/tp)*(w/tp);
    FPType w2 = upfac;
    while (!terminate)//assume uniform convergence
    {
        FPType rzarg = static_cast<FPType>(2*j+2);
//        FPType rz = mytr1::__detail::__riemann_zeta(rzarg);
	FPType rz = evenzeta<FPType>(rzarg);
//        std::cout<<rz<<" "<<fac<<" "<<w2<<std::endl;
        FPType nextterm = rz*fac*w2;
        w2 *= upfac;
        fac *= rzarg/(rzarg + s) * (rzarg+1.0)/(rzarg + s + 1.0);
        ++j;
        terminate = (fpequal( std::abs(res - pref*nextterm), std::abs(res) ) || (j > maxit));
        res -= pref * nextterm;
    }
//    std::cout<<"Iterations in Integer Series: "<<j<<'\n';
    return std::complex<FPType>(res, imag(imagtemp));
}

/** This function catches the cases of negative real index s.
 * Theoretical convergence is present for |w| < 2*pi.
 * We use an optimized version of
 * Li_s(e^w) = Gamma(1-s)*(-w)^(s-1) + (2*pi)^(-s)/pi * A_p(w)
 * A_p(w)= \sum_k Gamma(1+k-s)/k!*Sin(pi/2*(s-k))*(w/2/\pi)^k*zeta(1+k-s)
 * @param s the index s
 * @param w
 */
template <typename FPType>
inline std::complex<FPType> PolyLog_Exp_neg(const FPType s, std::complex<FPType> w)
{   //basic general loop, but s is a negative quantity here
    //FIXME Large s makes problems. The series should be rearrangeable so that we only need the ratio Gamma(1-s)/(2 pi)^s
    std::cout<<"Negative real s - 3"<<std::endl;
    FPType ls = std::lgamma(1.0-s);
    std::complex<FPType> res = std::exp(ls - (1.0-s) * std::log(-w));
    constexpr FPType tp = 2.0 * M_PI;
    const std::complex<FPType> wup = w/tp;
    std::complex<FPType> w2 = wup;
    std::complex<FPType> pref = 2.0 * std::pow(tp, - (1.0-s));
    //here we factor up the ratio of Gamma(1 - s + k)/k! .
    //This ratio should be well behaved even for large k in the series afterwards
    //Note that we have a problem for large s
    //Since s is negative we evaluate the Gamma Function on the positive real axis where it is real.
    FPType gam = std::exp(ls);

    FPType sp, cp;
    sincos(M_PI/2.0 * s, &sp, &cp);
    //Here we add the expression that would result from ignoring the zeta function in the series.
    std::complex<FPType> expis(cp, sp);
    std::complex<FPType> p = tp - std::complex<FPType>(0.0, 1.0) * w;
    std::complex<FPType> q = tp + std::complex<FPType>(0.0, 1.0) * w;
    res += std::complex<FPType>(0.0, 1.0) * gam * (conj(expis) * std::pow(p, s-1.0) - expis *std::pow(q, s-1.0));//this can be optimized for real values of w
    //The above expression is the result of sum_k Gamma(1+k-s) /k! * sin(pi /2* (s-k)) * (w/2/pi)^k
    //Therefore we only need to sample values of zeta(n) on the real axis that really differ from one
    res += pref * (sp * gam * (mytr1::__detail::__riemann_zeta(1.0-s) - 1.0));
    constexpr unsigned int maxit = 200;
    unsigned int j = 1;
    bool terminate = false;
    gam*= (1.0 - s);
    while (!terminate)//assume uniform convergence
    {
        FPType rzarg = (1.0 - s) + j;
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
        std::complex<FPType> nextterm =  w2 * (gam * sine * rz);
//	std::cout<<j<<" "<<nextterm<<" "<<rz<<" "<<std::endl;
        w2 *= wup;
        ++j;
        gam  *= rzarg/(j);//equal to 1/(j+1) since we have incremented j in the line above
        terminate = (fpequal( std::abs(res + pref*nextterm), std::abs(res) ) || (j > maxit));
        res += pref*nextterm;
    }
    std::cout<<"Iterations in PolyLogExp_neg: "<<j<<std::endl;
    return res;
}

/** This function catches the cases of negative integer index s which are multiples of two. In that case the sine occuring in the expansion 
 * occasionally takes on the value zero. We use that to provide an optimized series for p = 2n:
 * In the template parameter sigma we transport whether p = 4k (sigma = 1) or p = 4k + 2  (sigma = -1)
 * Li_p(e^w) = Gamma(1-p) * (-w)^{p-1} - A_p(w) - sigma * B_p(w)
 * with
 * A_p(w) = 2 (2\pi)^(p-1) (-p)! / (2 \pi)^(-p/2) (1 + w^2/(4 pi^2))^{-1/2 + p/2} cos((1 - p) ArcTan(2 pi/ w))
 * and 
 * B_p(w) = - 2 (2 pi)^(p-1) * \sum \limits_{k = 0}^\infty \Gamma(2 + 2k - p)/ (2k+1)! (-1)^k (w/2/\pi)^(2k+1) (Zeta(2 + 2k - p) - 1.0)
 * This is suitable for |w| < 2 pi
 * The original series is (This might be worthwhile if we use the already present table of the Bernoullis)
 * Li_p(e^w) = Gamma(1-p) * (-w)^{p-1} - sigma (2 pi)^p / pi * \sum \limits_{k = 0}^\infty \Gamma(2 + 2k - p)/ (2k+1)! (-1)^k (w/2/\pi)^(2k+1) Zeta(2 + 2k - p)
 * @param n the index n = 4k
 * @param w
 */
template <typename FPType, int sigma>
inline std::complex<FPType> PolyLog_Exp_neg_even(const uint n, std::complex<FPType> w)
{
//  std::cout<<"Negative even integer s = -2k , - 4"<<std::endl;
  const uint np = 1+n;
  FPType lnp = std::lgamma(np);
  std::complex<FPType> res = std::exp(lnp - FPType(np) * std::log(-w));
  constexpr FPType tp = 2.0 * M_PI;
  std::complex<FPType> wup = w/tp;
  std::complex<FPType> wq = wup*wup;
  FPType pref = 2.0 * std::pow(tp, -int(1 + n));
  //subtract  the expression A_p(w)
  res -= std::exp(lnp - 0.5*np*std::log( 1.0 + wq)) * pref * std::cos( static_cast<FPType>(np) * std::atan(1.0/wup));
  uint k = 0;
  bool terminate = false;
  constexpr uint maxit = 300;
  FPType gam = std::tgamma(2+n);
  if(sigma != 1)
    pref = -pref;
  while(!terminate)
  {
//    std::complex<FPType> newterm = ( gam * (mytr1::__detail::__riemann_zeta(static_cast<FPType>(2*k + 2 + n)) - 1.0)) * wup;
    std::complex<FPType> newterm = ( gam * (evenzeta<FPType>(2*k + 2 + n) - 1.0)) * wup;
    gam *= - static_cast<FPType>(2 * k + 2 + n + 1) / (2*k + 2 + 1) * static_cast<FPType>(2*k + 2 + n) / (2 * k + 1 + 1);
    wup *= wq;
    terminate = (fpequal( std::abs(res - pref*newterm), std::abs(res) ) || (k > maxit));
    res -= pref*newterm;
    ++k;
  }
//  std::cout<<"Iterations in the series for s = -4n : "<<k<<'\n';
  return res;
}

/** This function catches the cases of negative integer index s which are odd. In that case the sine occuring in the expansion 
 * occasionally takes on the value zero. We use that to provide an optimized series for p = 1 + 2k:
 * Int the template parameter sigma we transport whether p = 1 + 4k (sigma = 1) or p = 3 + 4k  (sigma = -1)
 * Li_p(e^w) = Gamma(1-p) * (-w)^{p-1} + sigma * A_p(w) - sigma * B_p(w)
 * with
 * A_p(w) = 2 (2\pi)^(p-1) * Gamma(1-p) (1 + w^2/(4 pi^2))^{-1/2 + p/2} cos((1 - p) ArcTan(2 pi/ w))
 * and 
 * B_p(w) = 2 (2 pi)^(p-1) * \sum \limits_{k = 0}^\infty \Gamma(1 + 2k - p)/ (2k)! (-w^2/4/\pi^2)^k (Zeta(1 + 2k - p) - 1.0)
 * This is suitable for |w| < 2 pi .
 * The use of evenzeta gives a speedup of about 50
 * The original series is (This might be worthwhile if we use the already present table of the Bernoullis)
 * Li_p(e^w) = Gamma(1-p) * (-w)^{p-1} - sigma *2*(2 pi)^(p-1) * \sum \limits_{k = 0}^\infty \Gamma(1 + 2k - p)/ (2k)! (-1)^k (w/2/\pi)^(2k) Zeta(1 + 2k - p)
 * @param n the index n = 4k
 * @param w
 */
template <typename FPType, int sigma>
inline std::complex<FPType> PolyLog_Exp_neg_odd(const uint n, std::complex<FPType> w)
{
//  std::cout<<"Negative odd integer s = -(1 + 2k), - 5"<<std::endl;
  const uint np = 1+n;
  FPType lnp = std::lgamma(np);
  std::complex<FPType> res = std::exp(lnp - FPType(np) * std::log(-w));
  constexpr FPType itp = 1.0/(2.0 * M_PI);
  std::complex<FPType> wq = -w * itp * w*itp;
  FPType pref = 2.0 * std::pow(itp, np);
  //subtract  the expression A_p(w)
  res += std::exp(lnp -0.5*np*std::log(1.0 - wq))* pref * std::cos( static_cast<FPType>(np) * std::atan(2.0 * M_PI/w));
  if(sigma != 1)
    pref = -pref;
  bool terminate = false;
  constexpr uint maxit = 300;
  FPType gam = std::exp(lnp);
  //zeroth order
  res -= pref * gam * (evenzeta<FPType>(np) - 1.0);
  uint k = 0;
  std::complex<FPType> wup = wq;
  while(!terminate)
  {
    uint zk = 2*k;
    gam *= static_cast<FPType>(zk + np)/(1 + zk) * static_cast<FPType>(1+zk + np) / (zk+2);
//    std::complex<FPType> newterm = ( gam * (mytr1::__detail::__riemann_zeta(static_cast<FPType>(zk + 2 + np)) - 1.0)) * wup;
    std::complex<FPType> newterm = ( gam * (evenzeta<FPType>(zk + 2 + np) - 1.0)) * wup;
    wup *= wq;
    terminate = (fpequal( std::abs(res - pref*newterm), std::abs(res) ) || (k > maxit));
    res -= pref*newterm;
    ++k;
  }
//  std::cout<<"Iterations in the series for s = -(1+2*k) : "<<k<<'\n';
  return res;
}

/** This function catches the cases of negative integer index s
 * @param s the index s
 * @param w
 */
template <typename FPType>
inline std::complex<FPType> PolyLog_Exp_neg(const int s, std::complex<FPType> w)
{//negative integer s
  const uint n = -s;
  switch(n%4)
  {
    case 0:
      return PolyLog_Exp_neg_even<FPType, 1>(n, w);
    case 1:
      return PolyLog_Exp_neg_odd<FPType, 1>(n, w);
    case 2:
      return PolyLog_Exp_neg_even<FPType, -1>(n, w);
    case 3:
      return PolyLog_Exp_neg_odd<FPType, -1>(n, w);
      break;
  }
}

/** This function catches the cases of positive real index s.
 * The defining series is
 * Li_s(e^w) = A_s(w) + B_s(w)+ \Gamma(1-s)(-w)^(s-1)
 * with
 * A_s(w) = \sum_{k=0}^{m} \zeta(s-k)w^k/k!
 * B_s(w) = \sum_{k=m+1}^\infty \sin(\pi/2(s-k)) \Gamma(1-s+k)\zeta(1-s+k) (w/2/\pi)^k/k!
 * @param s the index s
 * @param w
 */
template <typename FPType>
inline std::complex<FPType> PolyLog_Exp_pos(const FPType s, std::complex<FPType> w)
{   //positive s
    std::cout<<"Series for real positive s - 6"<<std::endl;
    std::complex<FPType> res = mytr1::__detail::__riemann_zeta(s);
    std::complex<FPType> wpower = w;
    FPType fac = 1.0;
    const uint m = static_cast<uint>(std::floor(s));
    for (uint k = 1; k <= m; ++k)
    {
        res += wpower*fac*mytr1::__detail::__riemann_zeta(static_cast<FPType>(s - k));
        wpower *= w;
        FPType temp = 1.0/(1.0 + k);
        fac *= temp;
    }
    //fac should be 1/(m+1)!
    //We revert here to the plain evaluation of the Gamma function instead of lgamma since we require evaluation at negative values.
    res += std::tgamma(1.0-s)*std::pow(-w, s-1.0);
    const FPType tp = 2.0 * M_PI;
    const FPType pref = 2.0 * std::pow(tp, s-1);
    //now comes the remainder of the series
    const unsigned int maxit = 100;
    unsigned int j = 0;
    bool terminate = false;
    std::complex<FPType> wup = w/tp;
    std::complex<FPType> w2 = std::pow(wup, m+1);
    //It is 1 < 2 - s + m < 2 => Gamma(2-s+m) will not overflow
    FPType gam = std::tgamma(2.0-s+m)*fac; //here we factor up the ratio of Gamma(1 - s + k)/k! . This ratio should be well behaved even for large k
    FPType sp, cp;
    sincos(M_PI/2.0 * s, &sp, &cp);
    while (!terminate)//assume uniform convergence
    {   //FIXME: optimize.
        int idx = m + 1 + j;
        FPType zetaarg = 1 + idx - s;
        FPType sine;
        if(idx & 1)//save the repeated calculation of the sines
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
 * For real u the imaginary part of the PolyLog is given by Im(Li_s(e^u)) = - \pi u^{s-1}/Gamma(s)
 * Check this relation for any benchmark that you use.
 * The use of evenzeta leads to a speedup of about 1000.
 * @param s the index s
 * @param w
 */
template <typename FPType>
inline std::complex<FPType> PolyLog_Exp_asym(const FPType s, std::complex<FPType> w)
{   //asymptotic expansion
    std::cout<<"asymptotic expansions , -7 "<<std::endl;
    std::complex<FPType> wgamma = std::exp((s-1.0)*std::log(w) - std::lgamma(s));/*wgamma = w^(s-1)/Gamma(s)*/
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
//        newterm = mytr1::__detail::__riemann_zeta(static_cast<FPType> (2*k) ) * wgamma;
	newterm = evenzeta<FPType>(2*k) * wgamma;
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
    std::cout<<"Iterations: "<<k<<'\n';
    return res;
}

/**
 * Theoretical convergence for Re(w) < 0. Seems to beat the other expansions for Re(w) < -pi/2 - pi/5.
 * Note that this is an implementation of the basic series:
 * Li_s(e^z) = \sum_{k=1} e^(k*z) * k^(-s)
 * @param s is an arbitrary type, Integer or float
 * @param w something with a negative real part
 */
template <typename PowerType, typename T>
inline T PolyLog_Exp_negative_real_part(PowerType s, T w)
{
  std::cout<<"negative real part series (exponential) - 8"<<std::endl;
  T ew = std::exp(w);
  const T up = ew;
  T res = ew;
  uint maxiter = 500;
  bool terminate = false;
  uint k = 2;
  while(!terminate)
  {
    ew *= up;
    T temp = std::pow(k, s);//This saves us a type conversion
    T newterm = ew / temp;
    terminate = (fpequal(std::abs(res + newterm), std::abs(res))) || (k > maxiter);
    res += newterm;
    ++k;
  }
  std::cout<<"iterations in PolyLog_Exp_negative_real_part: "<<k<<std::endl;
  return res;
}

/* This is the case where s is a positive integer.
 */
template <typename FPType>
inline std::complex<FPType> PolyLog_Exp_int_pos(const uint s, std::complex<FPType> w)
{
    FPType rw = w.real();
    FPType iw = w.imag();
    if(fpequal(rw, 0.0) && fpequal(std::remainder(iw, 2.0*M_PI), 0.0))
    {
        if (s > 1)
            return mytr1::__detail::__riemann_zeta(FPType(s));
        else
            return std::numeric_limits<FPType>::infinity();
    }
        if(0 == s)
        {
            std::complex<FPType> t = std::exp(w);
            return t/(1.0 - t);
        }
        else if (1 == s)
            return -std::log(1.0 - std::exp(w));
        else
        {
	    if(rw < -(M_PI/2.0 + M_PI/5.0)   )
	    {
	      //choose the exponentially converging series
	      return PolyLog_Exp_negative_real_part(s, w);
	    }
            //The transition point chosen here, is quite arbitrary and needs more testing.
            if(rw < 6.0)
            {
	       /*The reductions of the imaginary part yield the same results as Mathematica.
               * Necessary to improve the speed of convergence
               */
                while (w.imag() > M_PI) w = std::complex<FPType>(w.real(), w.imag() - 2.0*M_PI);
                while (w.imag() <= -M_PI) w = std::complex<FPType>(w.real(), w.imag() + 2.0*M_PI);
                return PolyLog_Exp_pos(s , w);
            }
            else
            {
                //wikipedia says that this is required for Wood's formula
                while (w.imag() > 0) w = std::complex<FPType>(w.real(), w.imag() - 2.0*M_PI);
                while (w.imag() <= -2.0*M_PI) w = std::complex<FPType>(w.real(), w.imag() + 2.0*M_PI);
                return PolyLog_Exp_asym(static_cast<FPType>(s), w);//FIXME: the series should terminate after a finite number of terms.
            }
        }
}

/* This is the case where s is a positive integer. And w is supposed to be a real
 */
template <typename FPType>
inline std::complex<FPType> PolyLog_Exp_int_pos(const uint s, FPType w)
{
    if(fpequal(w, 0.0))
    {
        if (s > 1)
            return mytr1::__detail::__riemann_zeta(FPType(s));
        else
            return std::numeric_limits<FPType>::infinity();
    }
        if(0 == s)
        {
            FPType t = std::exp(w);
            return t/(1.0 - t);
        }
        else if (1 == s)
            return -std::log(1.0 - std::exp(w));
        else
        {
	    if(w < -(M_PI/2.0 + M_PI/5.0)   )
	    {
	      //choose the exponentially converging series
	      return PolyLog_Exp_negative_real_part(s, std::complex<FPType>(w));
	    }
            //The transition point chosen here, is quite arbitrary and needs more testing.
            if(w < 6.0)
            {
                return PolyLog_Exp_pos(s, w);
            }
            else
            {
                return PolyLog_Exp_asym(static_cast<FPType>(s), std::complex<FPType>(w));//FIXME: the series should terminate after a finite number of terms.
            }
        }
}

/* This is the case where s is a negative integer.
 */
template <typename FPType>
inline std::complex<FPType> PolyLog_Exp_int_neg(const int s, std::complex<FPType> w)
{
    if (( ((-s) & 1) == 0) && fpequal(real(w), 0.0))
    {
      //Now s is odd and w on the unit-circle
      FPType iw = imag(w);//get imaginary part
	    FPType rem = std::remainder(iw, 2.0*M_PI);
	    if(fpequal(std::abs(rem), 0.5))
	    {
	      //Due to: Li_{-n}(-1) + (-1)^n Li_{-n}(1/-1) = 0 
	      return 0.0;
	    }
	    else
	    {
	      return PolyLog_Exp_neg(s, std::complex<FPType>(w.real(), rem));//no asymptotic expansion available... check the reduction
	    }
	  }
	  else
	  {
          if(real(w) < -(M_PI/2.0 + M_PI/5.0)   )//choose the exponentially converging series
          {
	    return PolyLog_Exp_negative_real_part(s, w);
          }
            if(real(w) < 6.0)//arbitrary transition point...
            {
	       /*The reductions of the imaginary part yield the same results as Mathematica.
               * Necessary to improve the speed of convergence
               */
                while (w.imag() > M_PI) w = std::complex<FPType>(w.real(), w.imag() - 2.0*M_PI);
                while (w.imag() <= -M_PI) w = std::complex<FPType>(w.real(), w.imag() + 2.0*M_PI);
                return PolyLog_Exp_neg(s , w);
            }
            else
            {
                //wikipedia says that this is required for Wood's formula
                while (w.imag() > 0) w = std::complex<FPType>(w.real(), w.imag() - 2.0*M_PI);
                while (w.imag() <= -2.0*M_PI) w = std::complex<FPType>(w.real(), w.imag() + 2.0*M_PI);
                return PolyLog_Exp_asym(static_cast<FPType>(s), w);//FIXME: the series should terminate after a finite number of terms.
            }
	  }
}

/* This is the case where s is a negative integer. and w is a real
 */
template <typename FPType>
inline std::complex<FPType> PolyLog_Exp_int_neg(const int s, FPType w)
{
  if(w < -(M_PI/2.0 + M_PI/5.0)   )//choose the exponentially converging series
  {
    return PolyLog_Exp_negative_real_part(s, std::complex<FPType>(w));
  }
  if (fpequal(w, 0.0)) return std::numeric_limits<FPType>::infinity();
  if(w < 6.0)//arbitrary transition point...
  {
    return PolyLog_Exp_neg(s , std::complex<FPType>(w));
    
  }
  else
  {
    return PolyLog_Exp_asym(static_cast<FPType>(s), std::complex<FPType>(w));//FIXME: the series should terminate after a finite number of terms.
  }
}

/* This is the case where s is a positive real value.
 */
template <typename FPType>
inline std::complex<FPType> PolyLog_Exp_real_pos(const FPType s, std::complex<FPType> w)
{
    FPType rw = w.real();
    FPType iw = w.imag();
    if(fpequal(rw, 0.0) && fpequal(std::remainder(iw, 2.0*M_PI), 0.0))
    {
        if (s > 1.0)
            return mytr1::__detail::__riemann_zeta(s);
        else
            return std::numeric_limits<FPType>::infinity();
    }
  if(rw < -(M_PI/2.0 + M_PI/5.0)   )//choose the exponentially converging series
  {
    return PolyLog_Exp_negative_real_part(s, w);
  }
  if(rw < 6.0)//arbitrary transition point
        {
            /*The reductions of the imaginary part yield the same results as Mathematica then.
             * Necessary to improve the speed of convergence
             */
            while (w.imag() > M_PI) w = std::complex<FPType>(w.real(), w.imag() - 2.0*M_PI);//branch cuts??
            while (w.imag() <= -M_PI) w = std::complex<FPType>(w.real(), w.imag() + 2.0*M_PI);
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

/* This is the case where s is a positive real value. and w is a plain real.
 */
template <typename FPType>
inline std::complex<FPType> PolyLog_Exp_real_pos(const FPType s, FPType w)
{
  if(fpequal(w, 0.0))
  {
      if (s > 1.0)
          return mytr1::__detail::__riemann_zeta(s);
      else
          return std::numeric_limits<FPType>::infinity();
  }
  if(w < -(M_PI/2.0 + M_PI/5.0)   )//choose the exponentially converging series
  {
    return PolyLog_Exp_negative_real_part(s, w);
  }
  if(w < 6.0)//arbitrary transition point
  {
    return PolyLog_Exp_pos(s, std::complex<FPType>(w));
  }
  else
  {
    return PolyLog_Exp_asym(s, std::complex<FPType>(w));
  }
}

/* This is the case where s is a negative real value.
 * Now we branch depending on the properties of w in the specific functions
 */
template <typename FPType>
inline std::complex<FPType> PolyLog_Exp_real_neg(const FPType s, std::complex<FPType> w)
{
  FPType rw = w.real();
  FPType iw = w.imag();
  if(rw < -(M_PI/2.0 + M_PI/5.0)   )//choose the exponentially converging series
  {
    return PolyLog_Exp_negative_real_part(s, w);
  }
  if(rw < 6)//arbitrary transition point
        {
            /*The reductions of the imaginary part yield the same results as Mathematica then.
             * Necessary to improve the speed of convergence
             */
            while (w.imag() > M_PI) w = std::complex<FPType>(w.real(), w.imag() - 2.0*M_PI);//branch cuts??
            while (w.imag() <= -M_PI) w = std::complex<FPType>(w.real(), w.imag() + 2.0*M_PI);
	    return PolyLog_Exp_neg(s, w);
        }
        else
        {
            //wikipedia says that this is required for Wood's formula
            while (w.imag() > 0) w = std::complex<FPType>(w.real(), w.imag() - 2.0*M_PI);
            while (w.imag() <= -2.0*M_PI) w = std::complex<FPType>(w.real(), w.imag() + 2.0*M_PI);
            return PolyLog_Exp_asym(s,w);
        }
}

/* This is the case where s is a negative real value.
 * Now we branch depending on the properties of w in the specific functions.
 * w is a real quantity
 */
template <typename FPType>
inline std::complex<FPType> PolyLog_Exp_real_neg(const FPType s, FPType w)
{
  if(w < -(M_PI/2.0 + M_PI/5.0)   )//choose the exponentially converging series
  {
    return PolyLog_Exp_negative_real_part(s, std::complex<FPType>(w));
  }
  if(w < 6)//arbitrary transition point
  {
    return PolyLog_Exp_neg(s, std::complex<FPType>(w));
  }
  else
  {
    return PolyLog_Exp_asym(s, std::complex<FPType>(w));
  }
}

/** This is the frontend function which calculates Li_s( e^w )
 * First we branch into different parts depending on the properties of s.
 * This function is the same irrespective of a real or complex w, hence the template parameter ArgType.
 * @param s the index s
 * @param w complex w
 * @return the value of Li_s(e^w).
 */
template <typename FPType, typename ArgType>
inline std::complex<FPType> PolyLog_Exp(const FPType s, ArgType w)
{
  if(s > 45.0)//cutoff chosen arbitrarily. Not much testing was involved
    return PolyLog_Exp_negative_real_part(s, w);
  std::complex<FPType> ret;
  if (fpequal<FPType>(std::rint(s), s))
    {
      //In this branch of the if statement, s is an integer
      int p = int(std::lrint(s));
      if(p > 0)
        ret = PolyLog_Exp_int_pos(p, w);
      else
	ret = PolyLog_Exp_int_neg(p, w);
    }
    else
    {
      if (s > 0)
      {
	ret = PolyLog_Exp_real_pos(s, w);
      }
      else
	ret = PolyLog_Exp_real_neg(s, w);
    }
    return ret;
}

/** A function to implement the PolyLog for two real arguments.
 * @param s The index s
 * @param x A real x
*/
template <typename FPType>
inline std::complex<FPType> PolyLog(const FPType s, FPType x)
{
    if (fpequal(x, 0.0)) return 0.0;//According to Mathematica
    if (x < 0)
    {   //use the square formula to access negative values.
        FPType xp = -x;
        FPType y = std::log(xp);
        return PolyLog_Exp(s, 2.0 * y) * std::pow(2.0, 1.0-s) - PolyLog_Exp(s, y);
    }
    else
    {
        FPType y = std::log(x);
        return PolyLog_Exp(s, y);
    }
}

/** A function to implement the PolyLog in those cases where we can calculate it.
 * @param s The index s
 * @param w A complex w
 */
template <typename FPType>
inline std::complex<FPType> PolyLog(const FPType s, std::complex<FPType> w)
{
  std::cout<<s<<" "<<w<<std::endl;
    if(fpequal(imag(w), 0.0))
        return PolyLog(s, real(w));
    else
        return PolyLog_Exp(s, std::log(w));
}

/** A function to implement Dirichlet's Eta function
 * @param w A w
*/
template <typename FPType>
inline std::complex<FPType> Dirichlet_eta(std::complex<FPType> w)
{
    if(fpequal(imag(w), 0.0))
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
    if(fpequal(imag(w), 0.0))
        return imag(PolyLog(w.real(), std::complex<FPType>(0.0, 1.0)));
    else
    {   std::cout<<"domain not (yet) supported!!"<<std::endl;
        std::__throw_domain_error(__N("Bad argument to Dirichlet Eta."));
    }
}

/** A function to implement Claussen's series Sl.
 * Notation and connection to polylog from wikipedia
 * @param w A  w
 * FIXME: Check the restriction to positive integers m.
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
 * FIXME: Check the restriction to positive integers m.
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
