#ifndef H_POLYNOMIAL
#define H_POLYNOMIAL

#include <vector>

namespace{
  template<typename CoeffType>
    unsigned GetSignChanges(const std::vector<CoeffType> &coeffs){
    unsigned num_roots = 0;
    for(size_t i = 0; i < coeffs.size(); ++i){
      if(coeffs.at(i) == 0) continue;
      for(size_t j = i + 1; j < coeffs.size(); ++j){
        if(coeffs.at(j) == 0){
          continue;
        }else{
          if((coeffs.at(i) < 0 && coeffs.at(j) > 0)
             || (coeffs.at(i) > 0 && coeffs.at(j) < 0)) ++num_roots;
          break;
        }
      }
    }
    return num_roots;
  }
}

template<typename CoeffType>
class Polynomial{
public:
  Polynomial() = default;
  Polynomial(const Polynomial &) = default;
  Polynomial& operator=(const Polynomial &) = default;
  Polynomial(Polynomial &&) = default;
  Polynomial& operator=(Polynomial &&) = default;
  ~Polynomial() = default;

  Polynomial(const CoeffType &coeff):
    coeffs_(1, coeff){
    Shrink();
  }

  Polynomial(const std::vector<CoeffType> &coeffs):
    coeffs_(coeffs){
    Shrink();
  }

  void SetCoefficient(std::size_t degree, CoeffType coefficient){
    if(degree >= coeffs_.size() && coefficient != 0){
      coeffs_.resize(degree + 1,0);
      coeffs_.at(degree) = coefficient;
    }else if(degree < coeffs_.size()){
      coeffs_.at(degree) = coefficient;
      if(coefficient == 0 && degree+1 == coeffs_.size()){
        Shrink();
      }
    }
  }

  CoeffType GetCoefficient(std::size_t degree) const{
    if(degree >= coeffs_.size()) return 0;
    return coeffs_.at(degree);
  }

  CoeffType LeadingCoefficient() const{
    return coeffs_.size() ? coeffs_.back() : 0;
  }

  std::size_t Degree() const{
    return coeffs_.size() ? coeffs_.size()-1 : 0;
  }

  void TruncateToDegree(std::size_t degree){
    if(degree < Degree()){
      coeffs_.resize(degree+1,0);
      Shrink();
    }
  }

  bool IsZero() const{
    return coeffs_.size() == 0;
  }

  operator std::vector<CoeffType>() const{
    return coeffs_;
  }

  template<typename ArgType>
  auto operator()(const ArgType &x) const -> decltype(CoeffType(0)*x){
    if(coeffs_.size() == 0) return 0;
    decltype(coeffs_.at(0)*x) value = coeffs_.back();
    for(size_t degree = coeffs_.size()-2; degree < coeffs_.size(); --degree){
      value = value * x + coeffs_.at(degree);
    }
    return value;
  }

  Polynomial & operator+=(const Polynomial &p){
    if(p.coeffs_.size() > coeffs_.size()){
      coeffs_.resize(p.coeffs_.size(),0);
    }
    for(size_t i = 0; i < p.coeffs_.size(); ++i){
      coeffs_.at(i) += p.coeffs_.at(i);
    }
    Shrink();
    return *this;
  }
 
  Polynomial & operator-=(const Polynomial &p){
    if(p.coeffs_.size() > coeffs_.size()){
      coeffs_.resize(p.coeffs_.size(),0);
    }
    for(size_t i = 0; i < p.coeffs_.size(); ++i){
      coeffs_.at(i) -= p.coeffs_.at(i);
    }
    Shrink();
    return *this;
  }

  Polynomial & operator*=(const Polynomial &p){
    std::vector<CoeffType> output(Degree()+p.Degree()+1, 0);
    for(size_t deg_a = 0; deg_a < coeffs_.size(); ++deg_a){
      for(size_t deg_b = 0; deg_b < p.coeffs_.size(); ++deg_b){
        output.at(deg_a+deg_b) += coeffs_.at(deg_a)*p.coeffs_.at(deg_b);
      }
    }
    coeffs_ = output;
    Shrink();
    return *this;
  }

  Polynomial & operator/=(const Polynomial &p){
    Polynomial quo;
    Polynomial rem = *this;
    while(!rem.IsZero() && rem.Degree()>=p.Degree()){
      size_t before_deg = rem.Degree();
      CoeffType t = rem.LeadingCoefficient()/p.LeadingCoefficient();
      Polynomial term;
      term.SetCoefficient(rem.Degree()-p.Degree(), t);
      quo += term;
      rem -= term*p;
      rem.TruncateToDegree(before_deg-1);
    }
    *this = quo;
    Shrink();
    return *this;
  }
 
  Polynomial & operator%=(const Polynomial &p){
    Polynomial quo(0);
    Polynomial rem = *this;
    size_t pdeg = p.Degree();
    CoeffType pcoeff_inv = 1./p.LeadingCoefficient();
    while(!rem.IsZero() && rem.Degree()>=pdeg){
      size_t before_deg = rem.Degree();
      CoeffType t = rem.LeadingCoefficient()*pcoeff_inv;
      Polynomial term;
      term.SetCoefficient(before_deg-pdeg, t);
      quo += term;
      rem -= term*p;
      if(before_deg == 0){
        break;
      }else{
        rem.TruncateToDegree(before_deg-1);
      }
    }
    *this = rem;
    Shrink();
    return *this;
  }

  Polynomial operator-() const{
    Polynomial out = *this;
    for(auto &c: out.coeffs_){
      c = -c;
    }
    return out;
  }

  Polynomial operator+() const{
    return *this;
  }

  bool operator!=(const Polynomial &p){
    return coeffs_!=p.coeffs_;
  }
 
  bool operator==(const Polynomial &p){
    return coeffs_==p.coeffs_;
  }
 
private:
  std::vector<CoeffType> coeffs_;

  void Shrink(){
    if(coeffs_.size() == 0) return;

    bool is_zero = true;
    size_t degree = 0;
    for(size_t i = coeffs_.size()-1; i < coeffs_.size() && is_zero; --i){
      if(coeffs_.at(i) != 0){
        degree = i;
        is_zero = false;
      }
    }
    
    coeffs_.resize(is_zero ? 0 : degree + 1, 0);
  }
};

template<typename CoeffType>
Polynomial<CoeffType> operator+(Polynomial<CoeffType> a, Polynomial<CoeffType> b){
  return a+=b;
}

template<typename CoeffType>
Polynomial<CoeffType> operator-(Polynomial<CoeffType> a, Polynomial<CoeffType> b){
  return a-=b;
}

template<typename CoeffType>
Polynomial<CoeffType> operator*(Polynomial<CoeffType> a, Polynomial<CoeffType> b){
  return a*=b;
}

template<typename CoeffType>
Polynomial<CoeffType> operator/(Polynomial<CoeffType> a, Polynomial<CoeffType> b){
  return a/=b;
}

template<typename CoeffType>
Polynomial<CoeffType> operator%(Polynomial<CoeffType> a, Polynomial<CoeffType> b){
  return a%=b;
}

template<typename CoeffType>
std::ostream & operator<<(std::ostream &stream, const Polynomial<CoeffType> &p){
  stream << "Poly[" << p.Degree() << "](";
  for(size_t deg = 0; deg < p.Degree(); ++deg){
    stream << p.GetCoefficient(deg) << ',';
  }
  stream << p.GetCoefficient(p.Degree()) << ')';
  return stream;
}

template<typename CoeffType>
Polynomial<CoeffType> Derivative(const Polynomial<CoeffType> &poly, size_t n = 1){
  auto coeffs = static_cast<std::vector<CoeffType> >(poly);
  for(; n != 0; --n){
    for(size_t from_degree = 1; from_degree < coeffs.size(); ++from_degree){
      coeffs.at(from_degree-1) = from_degree * coeffs.at(from_degree);
    }
    coeffs.pop_back();
  }
  return coeffs;
}

template<typename CoeffType>
std::vector<Polynomial<CoeffType> > SturmSequence(const Polynomial<CoeffType> &poly){
  std::vector<Polynomial<CoeffType> > out(1, poly);
  if(poly.Degree() == 0) return out;
  out.push_back(Derivative(poly));
  size_t i = 0;
  auto negrem = -(out.at(i) % out.at(i+1));
  while(!negrem.IsZero()){
    //DBG(negrem);
    out.push_back(negrem);
    ++i;
    negrem  = -(out.at(i) % out.at(i+1));
  }
  return out;
}

template<typename CoeffType>
unsigned NumRoots(const Polynomial<CoeffType> &poly){
  auto sturm = SturmSequence(poly);
  std::vector<CoeffType> pos_seq, neg_seq;
  for(const auto &p: sturm){
    pos_seq.push_back(p.LeadingCoefficient());
    neg_seq.push_back(p.Degree()%2 == 0 ? p.LeadingCoefficient() : -p.LeadingCoefficient());
  }
  
  int changes_pos = GetSignChanges(pos_seq);
  int changes_neg = GetSignChanges(neg_seq);
  
  return changes_neg - changes_pos;
}

template<typename CoeffType, typename ArgType = double>
unsigned NumRoots(const Polynomial<CoeffType> &poly, const ArgType &low, const ArgType &high){
  auto sturm = SturmSequence(poly);
  std::vector<CoeffType> pos_seq, neg_seq;
  for(const auto &p: sturm){
    pos_seq.push_back(p(high));
    neg_seq.push_back(p(low));
  }
  int changes_pos = GetSignChanges(pos_seq);
  int changes_neg = GetSignChanges(neg_seq);
  
  return changes_neg - changes_pos;
}

#endif
