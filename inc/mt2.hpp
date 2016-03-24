#ifndef H_MT2
#define H_MT2

#include "TLorentzVector.h"

class MT2{
public:
  MT2() = default;
  MT2(const MT2 &) = default;
  MT2& operator=(const MT2 &) = default;
  MT2(MT2 &&) = default;
  MT2& operator=(MT2 &&) = default;
  ~MT2() = default;

  void SetTestMass(double test_mass);
  
  void SetMomenta(const TLorentzVector &visible_A,
                  const TLorentzVector &visible_B,
                  double invisible_px, double invisible_py);
  void SetMomenta(double visible_mass_a, double visible_px_a, double visible_py_a,
                  double visible_mass_b, double visible_px_b, double visible_py_b,
                  double invisible_px, double invisible_py);

  double GetMT2() const;
  double GetTrialMT2(double invisible_px_a,
                     double invisible_py_a) const;
  
  void GetInvisibleMomenta(TLorentzVector &invisible_a,
                           TLorentzVector &invisible_b) const;
  void GetInvisibleMomenta(double &invisible_px_a,
                           double &invisible_py_a,
                           double &invisible_px_b,
                           double &invisible_py_b) const;

  double IsUnbalanced() const;

private:
  double mA_, xA_, yA_;
  double mB_, xB_, yB_;
  double m_, x_, y_;

  mutable double deltaX_, deltaY_, mt2_, unbalanced_;
  mutable bool cached_momenta_, cached_mt2_;

  double rotation_angle_;
  bool sides_swapped_;

  void TransformMomenta();
  static void Rotate(double &x, double &y, double phi);

  void ComputeMT2() const;
  void ComputeUnbalancedMT2(double &lower_bound, double &upper_bound) const;
  static double ComputeMT(double m1, double x1, double y1,
                          double m2, double x2, double y2);
  void ComputeMomenta() const;
  
  unsigned GetNumSolutions(double mt) const;
  void GetCoefficients(double mt,
                       long double &c4,
                       long double &c3,
                       long double &c2,
                       long double &c1,
                       long double &c0) const;
};

#endif
