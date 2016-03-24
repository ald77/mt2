#include "mt2.hpp"

#include <cmath>

#include "polynomial.hpp"

using namespace std;

void MT2::SetTestMass(double test_mass){
  cached_momenta_ = false;
  cached_mt2_ = false;
  m_ = max(test_mass, 0.);
}

void MT2::SetMomenta(const TLorentzVector &visible_A,
                     const TLorentzVector &visible_B,
                     double invisible_px, double invisible_py){
  cached_momenta_ = false;
  cached_mt2_ = false;
  mA_ = max(visible_A.M(), 0.);
  xA_ = visible_A.Px();
  yA_ = visible_A.Py();
  mB_ = max(visible_B.M(), 0.);
  xB_ = visible_B.Px();
  yB_ = visible_B.Py();
  x_ = invisible_px;
  y_ = invisible_py;
  TransformMomenta();
}

void MT2::SetMomenta(double visible_mass_a, double visible_px_a, double visible_py_a,
                     double visible_mass_b, double visible_px_b, double visible_py_b,
                     double invisible_px, double invisible_py){
  cached_momenta_ = false;
  cached_mt2_ = false;
  mA_ = max(visible_mass_a, 0.);
  xA_ = visible_px_a;
  yA_ = visible_py_a;
  mB_ = max(visible_mass_b, 0.);
  xB_ = visible_px_b;
  yB_ = visible_py_b;
  x_ = invisible_px;
  y_ = invisible_py;
  TransformMomenta();
}

double MT2::GetMT2() const{
  ComputeMT2();
  return mt2_;
}

double MT2::GetTrialMT2(double invisible_px_a,
                        double invisible_py_a) const{
  Rotate(invisible_px_a, invisible_py_a, rotation_angle_);
  double invisible_px_b = x_ - invisible_px_a;
  double invisible_py_b = y_ - invisible_py_a;
  if(sides_swapped_){
    swap(invisible_px_a, invisible_px_b);
    swap(invisible_py_a, invisible_py_b);
  }
  double mta = ComputeMT(mA_, xA_, yA_, m_, invisible_px_a, invisible_py_a);
  double mtb = ComputeMT(mB_, xB_, yB_, m_, invisible_px_b, invisible_py_b);
  return mta>mtb ? mta : mtb;
}

void MT2::GetInvisibleMomenta(TLorentzVector &invisible_a,
                              TLorentzVector &invisible_b) const{
  ComputeMomenta();
  double invisible_px_a = 0.5*(x_+deltaX_);
  double invisible_py_a = 0.5*(y_+deltaY_);
  double invisible_px_b = 0.5*(x_-deltaX_);
  double invisible_py_b = 0.5*(y_-deltaY_);
  if(sides_swapped_){
    swap(invisible_px_a, invisible_px_b);
    swap(invisible_py_a, invisible_py_b);
  }
  Rotate(invisible_px_a, invisible_py_a, -rotation_angle_);
  Rotate(invisible_px_b, invisible_py_b, -rotation_angle_);
  invisible_a.SetXYZM(invisible_px_a, invisible_py_a, 0., m_);
  invisible_b.SetXYZM(invisible_px_b, invisible_py_b, 0., m_);
}

void MT2::GetInvisibleMomenta(double &invisible_px_a,
                              double &invisible_py_a,
                              double &invisible_px_b,
                              double &invisible_py_b) const{
  ComputeMomenta();
  invisible_px_a = 0.5*(x_+deltaX_);
  invisible_py_a = 0.5*(y_+deltaY_);
  invisible_px_b = 0.5*(x_-deltaX_);
  invisible_py_b = 0.5*(y_-deltaY_);
  if(sides_swapped_){
    swap(invisible_px_a, invisible_px_b);
    swap(invisible_py_a, invisible_py_b);
  }
  Rotate(invisible_px_a, invisible_py_a, -rotation_angle_);
  Rotate(invisible_px_b, invisible_py_b, -rotation_angle_);
}

double MT2::IsUnbalanced() const{
  ComputeMT2();
  return unbalanced_;
}

void MT2::TransformMomenta(){
  if(mA_ < mB_){
    swap(mA_, mB_);
    swap(xA_, xB_);
    swap(yA_, yB_);
    sides_swapped_ = true;
  }else{
    sides_swapped_ = false;
  }
  rotation_angle_ = -atan2(yA_, xA_);
  Rotate(xA_, yA_, rotation_angle_);
  Rotate(xB_, yB_, rotation_angle_);
  Rotate(x_, y_, rotation_angle_);
}

void MT2::Rotate(double &x, double &y, double phi){
  double pt = hypot(y, x);
  phi += atan2(y, x);
  x = pt * cos(phi);
  y = pt * sin(phi);
}

void MT2::ComputeMT2() const{
  if(cached_mt2_) return;

  double lower_bound, upper_bound;
  ComputeUnbalancedMT2(lower_bound, upper_bound);
  double middle = lower_bound + 0.5*(upper_bound - lower_bound);
  unbalanced_ = upper_bound-lower_bound;
  while(lower_bound < middle && middle < upper_bound){
    unsigned num_solutions = GetNumSolutions(middle);
    if(num_solutions > 1){
      upper_bound = middle;
    }else{
      lower_bound = middle;
    }
    middle = lower_bound + 0.5*(upper_bound - lower_bound);
  }
  
  mt2_ = upper_bound;
  cached_mt2_ = true;
}

void MT2::ComputeUnbalancedMT2(double &lower_bound, double &upper_bound) const{
  lower_bound = m_ + mA_;
  upper_bound = -1.;
  if(mA_ == 0. && m_ == 0. /* && mB_ == 0. implied since mA_ >= mB_*/){
    //Check for analytic solution in massless case
    double alpha = (x_*yB_ - xB_*y_)/(xA_*yB_ - xB_*yA_);
    double beta = (xA_*y_ - x_*yA_)/(xA_*yB_ - xB_*yA_);
    if(alpha >= 0. && beta >= 0.){
      deltaX_ = 2.*alpha*xA_ - x_;
      deltaY_ = 2.*alpha*yA_ - y_;
      lower_bound = 0.;
      upper_bound = 0.;
      cached_momenta_ = true;
    }
  }else if(mA_ > 0.){
    //Check for unbalanced case
    double xa_mina = m_/mA_ * xA_;
    double ya_mina = m_/mA_ * yA_;
    lower_bound = m_ + mA_;
    upper_bound = ComputeMT(mB_, xB_, yB_, m_, x_-xa_mina, y_-ya_mina);
    if(upper_bound <= lower_bound){
      upper_bound = lower_bound;
      deltaX_ = 2.*xa_mina - x_;
      deltaY_ = 2.*ya_mina - y_;
      cached_momenta_ = true;
    }
  }
  if(lower_bound > upper_bound){
    //Couldn't find analytic solution. Try to set good bounds...
    lower_bound = m_ + mA_;
    double mta = ComputeMT(mA_, xA_, yA_, m_, 0.5*x_, 0.5*y_);
    double mtb = ComputeMT(mB_, xB_, yB_, m_, 0.5*x_, 0.5*y_);
    upper_bound = max(mta, mtb);
    if(upper_bound < lower_bound){
      //Can only happen due to rounding errors in ComputeMT
      upper_bound = lower_bound;
    }
  }
}

double MT2::ComputeMT(double m1, double x1, double y1,
                      double m2, double x2, double y2){
  return sqrt(m1*m1+m2*m2+2.*(sqrt((m1*m1+x1*x1+y1*y1)*(m2*m2+x2*x2+y2*y2))-x1*x2-y1*y2));
}

void MT2::ComputeMomenta() const{
  if(cached_momenta_) return;
  ComputeMT2();

  //Find deltaX_ as root of quartic
  long double c4, c3, c2, c1, c0;
  GetCoefficients(mt2_, c4, c3, c2, c1, c0);
  double high = 10000.;
  double low = -high;
  unsigned num_solutions = NumRoots(Polynomial<long double>({c0, c1, c2, c3, c4}), low, high);
  while(num_solutions == 0 && high < numeric_limits<double>::max()){
    low*=2.;
    high*=2.;
    num_solutions = NumRoots(Polynomial<long double>({c0, c1, c2, c3, c4}), low, high);
  }
  double middle = low + 0.5*(high-low);
  while(low < middle && middle < high){
    num_solutions = NumRoots(Polynomial<long double>({c0, c1, c2, c3, c4}), low, middle);
    if(num_solutions == 0){
      low=middle;
    }else{
      high=middle;
    }
    middle = low + 0.5*(high-low);
  }
  deltaX_ = middle;

  //Get two possible deltaY from knowin mt of side A and deltaX_
  long double T2 = 4.*m_*m_ + x_*x_ + y_*y_;
  long double zetaA2 = xA_*x_ + yA_*y_ + mt2_*mt2_ - mA_*mA_ - m_*m_;
  long double ETA2 = mA_*mA_ + xA_*xA_ + yA_*yA_;
  long double XA2 = mA_*mA_ + xA_*xA_;
  long double YA2 = mA_*mA_ + yA_*yA_;
  long double kappaA2 = xA_*yA_;
  long double chiA3 = ETA2*x_ - xA_*zetaA2;
  long double gammaA3 = ETA2*y_ - yA_*zetaA2;
  long double lambdaA4 = ETA2*T2 - zetaA2*zetaA2;
  long double xiA4 = kappaA2*kappaA2 - XA2*YA2;
  long double rhoA5 = kappaA2*gammaA3 + XA2*chiA3;
  long double sigmaA6 = gammaA3*gammaA3 - XA2*lambdaA4;
  double dya1 = (kappaA2*deltaX_-gammaA3+sqrt(xiA4*deltaX_*deltaX_-2.*rhoA5*deltaX_+sigmaA6))/XA2;
  double dya2 = (kappaA2*deltaX_-gammaA3-sqrt(xiA4*deltaX_*deltaX_-2.*rhoA5*deltaX_+sigmaA6))/XA2;

  //Pick the one with the right mt2
  double mta1 = ComputeMT(mA_, xA_, yA_, m_, 0.5*(x_+deltaX_), 0.5*(y_+dya1));
  double mtb1 = ComputeMT(mB_, xB_, yB_, m_, 0.5*(x_-deltaX_), 0.5*(y_-dya1));
  double mta2 = ComputeMT(mA_, xA_, yA_, m_, 0.5*(x_+deltaX_), 0.5*(y_+dya2));
  double mtb2 = ComputeMT(mB_, xB_, yB_, m_, 0.5*(x_-deltaX_), 0.5*(y_-dya2));
  double mt21 = max(mta1, mtb1);
  double mt22 = max(mta2, mtb2);
  if(mt21 <= mt22){
    deltaY_ = dya1;
  }else{
    deltaY_ = dya2;
  }
  
  cached_momenta_ = true;
}

unsigned MT2::GetNumSolutions(double mt) const{
  long double c4, c3, c2, c1, c0;
  GetCoefficients(mt, c4, c3, c2, c1, c0);
  auto answer = NumRoots(Polynomial<long double>({c0, c1, c2, c3, c4}));
  return answer;
}

void MT2::GetCoefficients(double mt,
                          long double &c4,
                          long double &c3,
                          long double &c2,
                          long double &c1,
                          long double &c0) const{
  long double T2 = 4.*m_*m_ + x_*x_ + y_*y_;
  long double zetaA2 = xA_*x_ + yA_*y_ + mt*mt - mA_*mA_ - m_*m_;
  long double ETA2 = mA_*mA_ + xA_*xA_ + yA_*yA_;
  long double XA2 = mA_*mA_ + xA_*xA_;
  long double YA2 = mA_*mA_ + yA_*yA_;
  long double kappaA2 = xA_*yA_;
  long double chiA3 = ETA2*x_ - xA_*zetaA2;
  long double gammaA3 = ETA2*y_ - yA_*zetaA2;
  long double lambdaA4 = ETA2*T2 - zetaA2*zetaA2;
  long double xiA4 = kappaA2*kappaA2 - XA2*YA2;
  long double rhoA5 = kappaA2*gammaA3 + XA2*chiA3;
  long double sigmaA6 = gammaA3*gammaA3 - XA2*lambdaA4;
  long double zetaB2 = xB_*x_ + yB_*y_ + mt*mt - mB_*mB_ - m_*m_;
  long double ETB2 = mB_*mB_ + xB_*xB_ + yB_*yB_;
  long double XB2 = mB_*mB_ + xB_*xB_;
  long double YB2 = mB_*mB_ + yB_*yB_;
  long double kappaB2 = xB_*yB_;
  long double chiB3 = ETB2*x_ - xB_*zetaB2;
  long double gammaB3 = ETB2*y_ - yB_*zetaB2;
  long double lambdaB4 = ETB2*T2 - zetaB2*zetaB2;
  //long double xiB4 = kappaB2*kappaB2 - XB2*YB2;
  //long double rhoB5 = kappaB2*gammaB3 + XB2*chiB3;
  //long double sigmaB6 = gammaB3*gammaB3 - XB2*lambdaB4;
  long double D4 = XA2*kappaB2 - kappaA2*XB2;
  long double F5 = gammaA3*XB2 + XA2*gammaB3;
  long double G6 = XA2*XA2*YB2 - 2.*XA2*kappaA2*kappaB2 + kappaA2*kappaA2*XB2 + xiA4*XB2;
  long double H7 = XA2*gammaA3*kappaB2 - kappaA2*gammaA3*XB2 - rhoA5*XB2 - XA2*XA2*chiB3 - XA2*kappaA2*gammaB3;
  long double J8 = gammaA3*gammaA3*XB2 + sigmaA6*XB2 + 2.*XA2*gammaA3*gammaB3 + XA2*XA2*lambdaB4;
  c4 = G6*G6 - 4.*D4*D4*xiA4;
  c3 = 4.*(G6*H7 - 2.*xiA4*D4*F5 + 2.*rhoA5*D4*D4);
  c2 = 2.*(G6*J8 + 2.*H7*H7 - 2.*xiA4*F5*F5 + 8.*rhoA5*D4*F5 - 2.*sigmaA6*D4*D4);
  c1 = 4.*(H7*J8 - 2.*sigmaA6*D4*F5 + 2.*rhoA5*F5*F5);
  c0 = J8*J8 - 4.*sigmaA6*F5*F5;
}
