#include "MUQ/Approximation/Polynomials/Monomial.h"

using namespace muq::Approximation;

Monomial::Monomial() : IndexedScalarBasis() {}

Monomial::~Monomial() {}

double Monomial::DerivativeEvaluate(int const polyOrder, int const derivOrder, double const x) const {

    if((derivOrder > polyOrder) || (polyOrder==0))
        return 0.0;

    double c = 1.0;
    for(int k=polyOrder; k>polyOrder-derivOrder; --k)
        c *= k;

    return c*std::pow(x, polyOrder-derivOrder);

}

double Monomial::BasisEvaluate(int const order, double const x) const {
    return std::pow(x, order);
}

REGISTER_SCALARBASIS_FAMILY(Monomial)

void Monomial::MonomialDivision(Eigen::VectorXd const& A,
                                            Eigen::VectorXd const& B,
                                            Eigen::VectorXd      & Q,
                                            Eigen::VectorXd      & R){
  const int m = A.size()-1;
  const int n = B.size()-1;

  assert(m>0);
  assert(n>0);

  const int mMn = m - n;
  int i;

  if ( mMn < 0 ){
    R = -1.0*A;
    return;
  }

  assert(n);

  // first, compute the quotient
  const double iB0 = 1.0/B[n];
  int nj;
  Q.resize(mMn+1);
  for ( i = 0; i <= mMn; ++ i )
  {
    nj = std::min<int>(i,n);//i > n ? n : i;
    Q[mMn-i] = A[m-i];
    for ( int j = 1; j <= nj; ++ j )
      Q[mMn-i] -= B[n-j] * Q[mMn - i + j];

    Q[mMn-i] *= iB0;
  }

  // now, compute the remainder
  R.resize(n);
  for ( i = 1; i <= n; ++ i )
  {
    R[i-1] = A[i - 1];
    nj = std::min<int>(i, mMn + 1);//mMn + 1 > i ? i : mMn + 1;
    for ( int j = 0; j < nj; ++ j )
      R[i-1] -= B[i - 1 - j] * Q[j];

  }
}


double Monomial::MonomialEvaluate(Eigen::VectorXd const& P, double x) {

  const int Psize = P.size();
  double val = P(Psize-1);
  for (int i = Psize-2; i >=0; --i)
    val = val * x + P(i);

  return val;
}

Eigen::VectorXd Monomial::MonomialRoots(Eigen::VectorXd const& Pin, double tol) {
  assert(tol>0);

  // remove the high order zeros if there are any
  int polyOrder=Pin.size()-1;
  while(Pin(polyOrder)==0){
    --polyOrder;
  }

  // check the inputs
  if(polyOrder==0)
    return Eigen::VectorXd();

  if(polyOrder==1){
    return Eigen::VectorXd::Constant(1,-1.0*Pin(0)/Pin(1));

  }else if(polyOrder==2){

    // use the quadratic equation to find the roots
    const double a = Pin(2);
    const double b = Pin(1);
    const double c = Pin(0);
    const double part = b*b-4*a*c;
    if(part<0){
      return Eigen::VectorXd();
    }else{
      Eigen::VectorXd output;
      double root1 = (-b-sqrt(part))/(2.0*a);
      double root2 = (-b+sqrt(part))/(2.0*a);
      if(abs(root1-root2)>1e-14){
        output.resize(2);
        output(0) = fmin(root1,root2);
        output(1) = fmax(root1,root2);
      }else{
        output.resize(1);
        output(0) = 0.5*(root1+root2);
      }
      return output;
    }
  }

  // Estimate the roots in general
  Eigen::VectorXd P = Pin.head(polyOrder+1);

  // each term in the sturm sequence is a polynomial, this vector stores the coefficients
  std::vector<Eigen::VectorXd> sturmSeq(polyOrder+1);
  sturmSeq.at(0) = P;

  // initialize the second term in the sturm sequence with the polynomial derivative
  sturmSeq.at(1) = Eigen::VectorXd(polyOrder);
  for(int i=polyOrder-1; i>=0; --i)
    sturmSeq.at(1)(i) = (i+1)*P(i+1);

  // now use polynomial division to fill in the rest of the sturm sequence with remainders
  Eigen::VectorXd Q;
  for(int i=2; i<polyOrder+1; ++i){
    MonomialDivision(sturmSeq.at(i-2), sturmSeq.at(i-1), Q, sturmSeq.at(i));
    sturmSeq.at(i) *= -1.0;
  }

  // compute the number of sign changes at -infty and +infty.  infSigns holds the signs of each polynomial in the Sturm sequence at +\infty
  std::vector<bool> infSigns(polyOrder+1); // 0 = negative, 1 = positive
  infSigns.at(0) = sturmSeq.at(0)(sturmSeq.at(0).size()-1)>=0; // the highest order polynomial will dominate at +\infty, so the sign is equal to the sign of the coefficient

  std::vector<bool> minusInfSigns(polyOrder+1); //minusInfSigns holds the signs of each polynomial in the Sturm sequence at -\infty
  minusInfSigns.at(0) = (sturmSeq.at(0).size()%2==0) ? !infSigns.at(0) : infSigns.at(0);

  int numSignChanges[] = {0,0};
  for(int i=1; i<polyOrder+1; ++i){

    // the sign at +infty is just the sign of the largest coefficient
    infSigns.at(i) = sturmSeq.at(i)(sturmSeq.at(i).size()-1)>=0;
    if(infSigns.at(i)^infSigns.at(i-1))
      ++numSignChanges[1];

    // the sign at -infty needs to take into account whether the power is even or odd
    bool newIsOdd = sturmSeq.at(i).size()%2==0;
    minusInfSigns.at(i) = newIsOdd ? !infSigns.at(i) : infSigns.at(i);

    if(minusInfSigns.at(i)^minusInfSigns.at(i-1))
      ++numSignChanges[0];
  }

  int numRealRoots = numSignChanges[0]-numSignChanges[1];
  if(numRealRoots<=0){
    std::cout << "ERROR: No roots exist for monomial coefficients:\n  " << P.transpose() << std::endl;
    std::cout << "Sturm sequence is given by:\n";
    for(int i=0; i<sturmSeq.size(); ++i)
      std::cout << sturmSeq.at(i).transpose() << "  -->  " << infSigns.at(i) << " , " << minusInfSigns.at(i) << std::endl;
    assert(numRealRoots>0);
  }


  // first, find a lower bound and upper bound that bound ALL the real roots
  double lb = -1e3; // must be negative to start with
  double ub = 1e3;  // must be positive to start with

  bool foundLB = false;
  int midSignChanges;
  while(!foundLB){

    double oldPolyVal = MonomialEvaluate(sturmSeq.at(0),lb);

    // there is a chance we landed on a root exactly! This is currently a hack to get things to work, we should instead store this location and use it later!
    if(abs(oldPolyVal)<1e-13){
      lb -= 1;
      oldPolyVal = MonomialEvaluate(sturmSeq.at(0),lb);
    }

    midSignChanges = 0;
    for(int j=1; j<sturmSeq.size(); ++j){
      double polyVal = MonomialEvaluate(sturmSeq.at(j),lb);
      if(std::signbit(oldPolyVal)!=std::signbit(polyVal))
        ++midSignChanges;
      oldPolyVal = polyVal;
    }

    if(midSignChanges==numSignChanges[0]){
      foundLB = true;
    }else{
      lb *= 2.0;
    }
  }


  bool foundUB = false;
  while(!foundUB){

    double oldPolyVal = MonomialEvaluate(sturmSeq.at(0),ub);

    // there is a chance we landed on a root exactly! This is currently a hack to get things to work, we should instead store this location and use it later!
    if(abs(oldPolyVal)<1e-13){
      ub += 1;
      oldPolyVal = MonomialEvaluate(sturmSeq.at(0),ub);
    }

    midSignChanges = 0;
    for(int j=1; j<sturmSeq.size(); ++j){
      double polyVal = MonomialEvaluate(sturmSeq.at(j),ub);
      if(std::signbit(oldPolyVal)!=std::signbit(polyVal))
        ++midSignChanges;
      oldPolyVal = polyVal;
    }

    if(midSignChanges==numSignChanges[1]){
      foundUB = true;
    }else{
      ub *= 2.0;
    }
  }

  double bound_gap = ub-lb;

  // now that we know how many roots there are, we need to bracket them so we can use bisection to solve
  Eigen::VectorXd lowerBounds = Eigen::VectorXd::Constant(numRealRoots,lb); // lowerBounds will hold the lower value of each interval holding a root
  lowerBounds(0) = lb;

  double localTol = bound_gap; // the width of the search interval

  std::vector<int> lowerSignChanges(numRealRoots);
  std::vector<int> upperSignChanges(numRealRoots);
  lowerSignChanges.at(0) = numSignChanges[0];
  upperSignChanges.at(0) = numSignChanges[1];

  int numIntervals = 1; // how many intervals have we found that contain exactly one root


  while(numIntervals<numRealRoots && localTol > tol){ // bisection width
    localTol *= 0.5; // half the search interval

    int nloc = numIntervals;
    for(int i=0; i<nloc; ++i){ // loop over the number of intervals that have been found, initially just 1
      double x = lowerBounds(i) + localTol;

      double oldPolyVal = MonomialEvaluate(sturmSeq.at(0),x);

      // there is a chance we landed on a root exactly! This is currently a hack to get things to work, we should instead store this location and use it later!
      if(abs(oldPolyVal)<1e-13){
        x -= tol;
        oldPolyVal = MonomialEvaluate(sturmSeq.at(0),x);
      }


      midSignChanges = 0;
      for(int j=1; j<sturmSeq.size(); ++j){
        double polyVal = MonomialEvaluate(sturmSeq.at(j),x);
        if(std::signbit(oldPolyVal)!=std::signbit(polyVal))
          ++midSignChanges;
        oldPolyVal = polyVal;
      }

      if(midSignChanges==lowerSignChanges.at(i)){

        // if the number of sign changes is the same as before, no roots are to the left of x and we can update the lower bound
        lowerBounds(i)=x;

      }else if(midSignChanges != upperSignChanges.at(i)){
        // if we make it in here, we know that there are roots to the left and right of x, so we need to split the interval

        // if there is a real root to the left of x, x will serve as the lower bound in the search interval for the next root
        lowerBounds(numIntervals) = x;
        upperSignChanges.at(numIntervals) = upperSignChanges.at(i);

        upperSignChanges.at(i) = midSignChanges;
        lowerSignChanges.at(numIntervals) =  midSignChanges;
        numIntervals++;
      }
    }
  }

  std::sort(lowerBounds.data(),lowerBounds.data()+lowerBounds.size());

  /////////////////////////////////////////////////////////////////////////////////////////
  // Now that we've bracketed the roots, let's try to find them with a bisection solver
  Eigen::VectorXd rootLocs(numRealRoots);
  for(int i=0; i<numRealRoots; ++i){

    // bisection method
    double upperPos = (i<numRealRoots-1) ? lowerBounds(i+1) : ub;
    double lowerPos = lowerBounds(i);

    //double lowerVal = MonomialEvaluate(sturmSeq.at(0),lowerPos);
    double upperVal = MonomialEvaluate(sturmSeq.at(0),upperPos);
    bool foundExact = false;
    while(upperPos-lowerPos>tol){

      double midPos = 0.5*(upperPos-lowerPos)+lowerPos;
      double midVal = MonomialEvaluate(sturmSeq.at(0),midPos);

      // did we find a root exactly?
      if(std::abs<double>(midVal)<1e-13){
        rootLocs(i) = midPos;
        foundExact = true;
        break;
      }

      // continue bisecting the region
      if( (midVal>0)^(upperVal>0) ){
        lowerPos = midPos;
        //lowerVal = midVal;
      }else{
        upperPos = midPos;
        upperVal = midVal;
      }
    }

    if(!foundExact)
      rootLocs(i) = 0.5*(upperPos-lowerPos)+lowerPos;
  }

  // at this point, we have intervals containing all of our roots
  return rootLocs;
}

Eigen::VectorXd Monomial::GetMonomialCoeffs(unsigned int polyOrder) const {
  Eigen::VectorXd coeffs = Eigen::VectorXd::Zero(polyOrder+1);
  coeffs(polyOrder) = 1.0;
  return coeffs;
}
