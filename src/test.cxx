// Generates simple MC of decays of pair produced particles
// Computes mt2 using our calculator
// Computes mt2 with the bisection calculator as a reference for comparison

#include "test.hpp"

#include <cmath>

#include <iostream>
#include <random>
#include <array>
#include <functional>

#include "TFile.h"
#include "TLorentzVector.h"
#include "TTree.h"

#include "mt2_bisect.hpp"

#include "mt2.hpp"
#include "timer.hpp"
#include "polynomial.hpp"

using namespace std;

namespace{
  //Top mass and number of random events to test
  double mt_ = 172.5;
  size_t nevents_ = 1e3;

  //For random number distribution
  mt19937_64 prng_ = InitializePRNG();
  bernoulli_distribution bernoulli_(7./8.);
  uniform_real_distribution<> uniform_;
  uniform_real_distribution<> unihalf_(0., 0.5);
  uniform_real_distribution<> uni2pi_(0., 2.*acos(-1.));
  exponential_distribution<> exponential_(0.01);
}

int main(){
  mt2_bisect::mt2 theirs;
  MT2 mine;

  TFile file("mt2.root", "recreate");
  file.cd();
  TLorentzVector p4A, p4a, p4B, p4b, test_a, test_b;
  double test_mass, their_mt2, my_mt2, test_mt2, unbalanced;
  TTree tree("tree","tree");
  tree.Branch("A", &p4A);
  tree.Branch("a", &p4a);
  tree.Branch("B", &p4B);
  tree.Branch("b", &p4b);
  tree.Branch("test_a", &test_a);
  tree.Branch("test_b", &test_b);
  tree.Branch("test_mass", &test_mass);
  tree.Branch("theirs", &their_mt2);
  tree.Branch("mine", &my_mt2);
  tree.Branch("test_mt2", &test_mt2);
  tree.Branch("unbalanced", &unbalanced);
  Timer timer(nevents_, 1.);
  timer.Start();
  for(size_t event = 0; event < nevents_; ++event){
    timer.Iterate();
    double mA = bernoulli_(prng_) ? unihalf_(prng_)*mt_ : 0.;
    double mB = bernoulli_(prng_) ? unihalf_(prng_)*mt_ : 0.;
    double ma = bernoulli_(prng_) ? unihalf_(prng_)*mt_ : 0.;
    double mb = bernoulli_(prng_) ? unihalf_(prng_)*mt_ : 0.;
    double pxtA = bernoulli_(prng_) ? exponential_(prng_) : 0.;
    double pytA = bernoulli_(prng_) ? exponential_(prng_) : 0.;
    double pztA = bernoulli_(prng_) ? exponential_(prng_) : 0.;
    double pxtB = bernoulli_(prng_) ? exponential_(prng_) : 0.;
    double pytB = bernoulli_(prng_) ? exponential_(prng_) : 0.;
    double pztB = bernoulli_(prng_) ? exponential_(prng_) : 0.;
    double costA = uniform_(prng_);
    double phiA = uni2pi_(prng_);
    double costB = uniform_(prng_);
    double phiB = uni2pi_(prng_);
    test_mass =  bernoulli_(prng_) ? exponential_(prng_) : 0.;

    double thetaA = acos(costA);
    double thetaa = acos(-1.)-thetaA;
    double etaA = -log(tan(0.5*thetaA));
    double etaa = -log(tan(0.5*thetaa));
    double sumA = mA+ma;
    double diffA = mA-ma;
    double pA = sqrt((mt_*mt_-sumA*sumA)*(mt_*mt_-diffA*diffA))/(2.*mt_);

    TLorentzVector p4TA;
    p4TA.SetXYZM(pxtA, pytA, pztA, mt_);
    p4A.SetPtEtaPhiM(pA, etaA, phiA, mA);
    p4a.SetPtEtaPhiM(pA, etaa, -phiA, ma);
    p4A.Boost(p4TA.BoostVector());
    p4a.Boost(p4TA.BoostVector());

    double thetaB = acos(costB);
    double thetab = acos(-1.)-thetaB;
    double etaB = -log(tan(0.5*thetaB));
    double etab = -log(tan(0.5*thetab));
    double sumB = mB+mb;
    double diffB = mB-mb;
    double pB = sqrt((mt_*mt_-sumB*sumB)*(mt_*mt_-diffB*diffB))/(2.*mt_);
    TLorentzVector p4TB;
    p4TB.SetXYZM(pxtB, pytB, pztB, mt_);
    p4B.SetPtEtaPhiM(pB, etaB, phiB, mB);
    p4b.SetPtEtaPhiM(pB, etab, -phiB, mb);
    p4B.Boost(p4TB.BoostVector());
    p4b.Boost(p4TB.BoostVector());

    TLorentzVector p4c = p4a+p4b;

    double pa0[3] = {p4A.M(), p4A.Px(), p4A.Py()};
    double pb0[3] = {p4B.M(), p4B.Px(), p4B.Py()};
    double pmiss0[3] = {0., p4c.Px(), p4c.Py()};
    theirs.set_mn(test_mass);
    theirs.set_momenta(pa0, pb0, pmiss0);
    mine.SetMomenta(p4A, p4B, p4c.Px(), p4c.Py());
    mine.SetTestMass(test_mass);
    their_mt2 = theirs.get_mt2();
    my_mt2 = mine.GetMT2();
    mine.GetInvisibleMomenta(test_a, test_b);
    test_mt2 = mine.GetTrialMT2(test_a.Px(), test_a.Py());
    unbalanced = mine.IsUnbalanced();
    tree.Fill();
  }
  tree.Write();
  file.Close();
}

mt19937_64 InitializePRNG(){
  array<int, 128> sd;
  random_device r;
  generate_n(sd.begin(), sd.size(), ref(r));
  seed_seq ss(begin(sd), end(sd));
  return mt19937_64(ss);
}
