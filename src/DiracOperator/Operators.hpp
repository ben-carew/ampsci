#pragma once
#include "DiracOperator/DiracOperator.hpp"
#include <cmath>
#include <functional>
#include <string>
#include <vector>

namespace DiracOperator {

//******************************************************************************
//! @brief General function of r, even scalar operator
//! @details Pass only grid to just get r, or either pass a lambda/function
//! [f(r)], or a number, n, (for r^n).
class RadialF final : public ScalarOperator {
  // = some function of r
  // Pass only grid to just get r, or
  // either pass a lambda/function [f(r)], or a number, n, (for r^n)
  // Explicitely even with rank 0, so ...
private:
  std::vector<double> fillVec(const Grid &gr,
                              const std::function<double(double)> &f) {
    std::vector<double> f_r;
    f_r.reserve(gr.num_points);
    for (auto r : gr.r)
      f_r.push_back(f(r));
    return f_r;
  }

public:
  RadialF(const Grid &rgrid, const std::function<double(double)> &f)
      : ScalarOperator(Parity::even, 1.0, fillVec(rgrid, f)) {}
  RadialF(const Grid &rgrid, const double n)
      : ScalarOperator(Parity::even, 1.0, fillVec(rgrid, [n](double r) {
                         return std::pow(r, n);
                       })) {}
  std::string name() const override { return "RadialFunction"; }
  std::string units() const override { return "au"; }
};

//******************************************************************************
//! Electric dipole operator: -|e|r = -er
/*! @details
\f[<a||d||b> = R (-1)^{ja+1/2} \sqrt{[ja][jb]} \, tjs(ja,jb,1, -1/2,1/2,0)\f]
\f[<a||d||b> = R <k_a||C^k||k_b>\f]
\f[R = -e \int r(f_a f_b + g_a g_b) \, dr\f]
*/
class E1 final : public TensorOperator {
public:
  E1(const Grid &gr) : TensorOperator(1, Parity::odd, -1.0, gr.r, 0) {}

  double angularF(const int ka, const int kb) const override {
    return Angular::Ck_kk(1, ka, kb);
  }
  std::string name() const override { return "E1"; }
  std::string units() const override { return "aB"; }
};

//******************************************************************************
//! E^k (electric multipole) operator
class Ek final : public TensorOperator {
public:
  Ek(const Grid &gr, const int k)
      : TensorOperator(k, Angular::evenQ(k) ? Parity::even : Parity::odd, -1.0,
                       gr.rpow(k), 0),
        m_k(k) {}
  double angularF(const int ka, const int kb) const override {
    return Angular::Ck_kk(m_k, ka, kb);
  }
  std::string name() const override {
    return std::string("E") + std::to_string(m_k);
  }
  std::string units() const override {
    return std::string("aB^") + std::to_string(m_k);
  }

private:
  int m_k;
};

//******************************************************************************
//! @brief Electric dipole operator, v-form:
//! \f$ \frac{ie}{\omega \alpha} \vec{\alpha}\f$
/*! @details
\f[  <a||d_v||b> =  R \f]
\f[ R = -\frac{2e}{\omega \alpha}
\int( f_ag_b <ka||s||-kb> - g_af_b <-ka||s||kb>)\,dr \f]
*/
class E1_vform final : public TensorOperator
// d_v = (ie/w alpha) v{alpha}   [v{a} = g0v{g}]\f$
// <a||dv||b> = -2e/(w alpha) Int[ fagb <ka||s||-kb> - gafb <-ka||s||kb>]
{
public:
  E1_vform(const double alpha, const double omega = 0.0)
      : TensorOperator(1, Parity::odd, -0.0, {}, 0, Realness::real, true),
        m_alpha(alpha) {
    updateFrequency(omega);
  }
  std::string name() const override { return "E1v"; }
  std::string units() const override { return "aB"; }

  double angularF(const int, const int) const override { return 1.0; }

private:
  virtual double angularCff(int, int) const override { return 0; }
  virtual double angularCgg(int, int) const override { return 0; }
  virtual double angularCfg(int ka, int kb) const override {
    return Angular::S_kk(ka, -kb);
  }
  virtual double angularCgf(int ka, int kb) const override {
    return -Angular::S_kk(-ka, kb);
  }

  void updateFrequency(const double omega) override {
    m_constant = -2.0 / (m_alpha * omega);
  }

private:
  double m_alpha; // (including var-alpha)
};

//******************************************************************************
//! Magnetic dipole operator: <a||M1||b>
/*! @details
\f[ <a||M1||b> = R (k_a + k_b) <-k_a||C^1||k_b> \f]
\f[R = \frac{3}{\alpha^2\omega} \int (f_ag_b+g_af_b) j_1(kr) \, dr\f]
\f$ k = \omega/c = \omega*\alpha \f$.
Check sign?
*/
class M1 final : public TensorOperator {
public:
  // XXX Check sign!
  M1(const Grid &gr, const double alpha, const double omega = 0.0)
      : TensorOperator(1, Parity::even, +0.0, {}, 0, Realness::real, true),
        m_gr(&gr),
        m_alpha(alpha) {
    updateFrequency(omega);
  }
  M1 &operator=(const M1 &) = default;
  M1(const M1 &) = default;
  ~M1() = default;
  std::string name() const override { return std::string("M1"); }
  std::string units() const override { return std::string("mu_B"); }

  double angularF(const int ka, const int kb) const override {
    return (ka + kb) * Angular::Ck_kk(1, -ka, kb);
  }
  double angularCff(int, int) const override { return 0.0; }
  double angularCgg(int, int) const override { return 0.0; }
  double angularCfg(int, int) const override { return 1.0; }
  double angularCgf(int, int) const override { return 1.0; }

  void updateFrequency(const double omega) override {
    // XXX Check sign!
    if (std::abs(omega) > 0) {
      m_constant = +3.0 / (m_alpha * m_alpha * omega);
      m_vec = SphericalBessel::fillBesselVec_kr(1, omega * m_alpha, m_gr->r);
    } else {
      m_constant = +1.0 / (m_alpha);
      m_vec.clear();
    }
  }

private:
  const Grid *const m_gr;
  const double m_alpha;
};

//******************************************************************************
//! @brief Magnetic hyperfine operator
//! @details Note: 'hfs_F' function (magnetization distribtuion) **includes**
//! the 1/r^2 (slightly different to definition in paper)
class Hyperfine final : public TensorOperator {

  using Func_R2_R = std::function<double(double, double)>; // save typing

public: // F(r) functions. NOTE: includes 1/r^2 !
  static inline auto sphericalBall_F() -> Func_R2_R {
    return [=](double r, double rN) {
      return (r > rN) ? 1. / (r * r) : r / (rN * rN * rN);
    };
  }
  static inline auto sphericalShell_F() -> Func_R2_R {
    return [=](double r, double rN) { return (r > rN) ? 1. / (r * r) : 0.0; };
  }
  static inline auto pointlike_F() -> Func_R2_R {
    return [=](double r, double) { return 1. / (r * r); };
  }

  //! 'Volotka' single-particle model
  static inline auto volotkaBW_F(double mu, double I_nuc, double l_pn, int gl)
      -> Func_R2_R
  // Function that returns generates + returns F_BW Bohr-Weiskopf
  // gl = 1 for proton, =0 for neutron. Double allowed for testing..
  // mu is in units of Nuclear Magneton!
  {
    const auto two_I = int(2 * I_nuc + 0.0001);
    const auto two_l = int(2 * l_pn + 0.0001);
    auto g_l = double(gl); // just safety
    auto gI = mu / I_nuc;
    auto K = (l_pn * (l_pn + 1.0) - (3. / 4.)) / (I_nuc * (I_nuc + 1.0));
    const double g_s = (2.0 * gI - g_l * (K + 1.0)) / (1.0 - K);
    std::cout << "BW using: gl=" << g_l << ", gs=" << g_s << ", l=" << l_pn
              << ", gI=" << gI << " (I=" << I_nuc << ")\n";
    const double factor =
        (two_I == two_l + 1)
            ? g_s * (1 - two_I) / (4.0 * (two_I + 2)) + g_l * 0.5 * (two_I - 1)
            : g_s * (3 + two_I) / (4.0 * (two_I + 2)) +
                  g_l * 0.5 * two_I * (two_I + 3) / (two_I + 2);
    if (two_I != two_l + 1 && two_I != two_l - 1) {
      std::cerr << "\nFAIL:59 in Hyperfine (volotkaBW_F):\n "
                   "we must have I = l +/- 1/2, but we have: I,l = "
                << I_nuc << "," << l_pn << "\n";
      return [](double, double) { return 0.0; };
    }
    return [=](double r, double rN) {
      return (r > rN) ? 1. / (r * r)
                      : (r / (rN * rN * rN)) *
                            (1.0 - (3.0 / mu) * std::log(r / rN) * factor);
    };
  }

  //! 'Volotka' single-particle model, for doubly-odd nuclei
  static inline auto doublyOddBW_F(double mut, double It, double mu1, double I1,
                                   double l1, int gl1, double I2, double l2)
      -> Func_R2_R
  // F(r) * g = 0.5 [ g1F1 + g2F2 + (g1F1 - g2F2) * K]
  // K = [I1(I1+1) - I2(I2+1)] / [I(I+1)]
  // return F(r) [divide by g]
  // volotkaBW_F(mu, I, l, g_l); //gl is 1 or 0
  // g2 : from: g = 0.5 [ g1 + g2 + (g1 - g2) * K]
  {
    auto K = (I1 * (I1 + 1.0) - I2 * (I2 + 1.0)) / (It * (It + 1.0));
    auto gt = mut / It;
    auto g1 = mu1 / I1;
    auto g2 = (g1 * (K + 1.0) - 2.0 * gt) / (K - 1.0);
    auto mu2 = g2 * I2;
    auto gl2 = (gl1 == 0) ? 1 : 0;
    auto F1 = volotkaBW_F(mu1, I1, l1, gl1);
    auto F2 = volotkaBW_F(mu2, I2, l2, gl2);
    return [=](double r, double rN) {
      return (0.5 / gt) * (g1 * F1(r, rN) + g2 * F2(r, rN) +
                           K * (g1 * F1(r, rN) - g2 * F2(r, rN)));
    };
  }

private: // helper
  static inline std::vector<double> RadialFunc(double rN, const Grid &rgrid,
                                               const Func_R2_R &hfs_F) {
    std::vector<double> rfunc;
    rfunc.reserve(rgrid.num_points);
    for (auto r : rgrid.r)
      rfunc.push_back(hfs_F(r, rN));
    return rfunc;
  }

public: // constructor
  Hyperfine(double muN, double IN, double rN, const Grid &rgrid,
            const Func_R2_R &hfs_F = sphericalBall_F())
      : TensorOperator(1, Parity::even, muN * PhysConst::muN_CGS_MHz / IN,
                       RadialFunc(rN, rgrid, hfs_F), 0),
        Inuc(IN) {}
  std::string name() const override { return "hfs"; }
  std::string units() const override { return "MHz"; }

  double angularF(const int ka, const int kb) const override {
    return (ka + kb) * Angular::Ck_kk(1, -ka, kb);
  }

  static double convertRMEtoA(const DiracSpinor &Fa, const DiracSpinor &Fb) {
    return 0.5 / Fa.jjp1() / Angular::Ck_kk(1, -Fa.k, Fb.k);
    // Correct for diag. Off diag? Prob not defined?
  }

  double hfsA(const DiracSpinor &Fa) {
    auto Raa = radialIntegral(Fa, Fa);
    return Raa * Fa.k / (Fa.jjp1());
    // nb: in MHz
  }

  static double hfsA(const TensorOperator *h, const DiracSpinor &Fa) {
    auto Raa = h->radialIntegral(Fa, Fa);
    return Raa * Fa.k / (Fa.jjp1());
  }

  double de_F(const DiracSpinor &Fa, double jF) {
    auto Ahfs = hfsA(Fa); // nb: in MHz
    return 0.5 * Ahfs * (jF * (jF + 1.0) - Fa.jjp1() - Inuc * (Inuc + 1.0));
  }

private:
  double Inuc;

private:
  virtual double angularCff(int, int) const override { return 0; }
  virtual double angularCgg(int, int) const override { return 0; }
  virtual double angularCfg(int, int) const override { return 1.0; }
  virtual double angularCgf(int, int) const override { return 1.0; }
};

//******************************************************************************
//! Radiative QED operator, electric part
class Hrad_el final : public ScalarOperator {
public:
  Hrad_el(const std::vector<double> &Hel)
      : ScalarOperator(Parity::even, 1.0, Hel, {1, 0, 0, 1}) {}
  std::string name() const override { return "Hrad_el"; }
  std::string units() const override { return "au"; }
};
//! Radiative QED operator, off-diagonal magnetic part
class Hrad_mag final : public ScalarOperator {
public:
  Hrad_mag(const std::vector<double> &Hmag)
      : ScalarOperator(Parity::even, -137.036, Hmag, {0, 1, 1, 0}) {}
  std::string name() const override { return "Hrad_mag"; }
  std::string units() const override { return "au"; }
};

//******************************************************************************
//! @brief Nuclear-spin independent PNC operator (Qw)
/*! @details
\f[
h_{PNCnsi} = \frac{G_F \, Q_W}{2\sqrt{2}} \, \rho(r) \, \gamma_5.
\f]
  - Check -ve sign?
  - Ouput is in units of (Qw * 1.e-11.) by default. To get (Qw/-N), multiply
  by (-N) [can go into optional 'factor']

Generates rho(r) according to fermi distribution, given c and t [c and t in
FERMI / femptometers].
*/
class PNCnsi final : public ScalarOperator {
public:
  PNCnsi(double c, double t, const Grid &rgrid, double factor = 1,
         const std::string &in_units = "iQw*e-11")
      : ScalarOperator(Parity::odd, factor * PhysConst::GFe11 / std::sqrt(8.0),
                       Nuclear::fermiNuclearDensity_tcN(t, c, 1, rgrid),
                       {0, 1, -1, 0}, 0, Realness::imaginary),
        m_unit(in_units) {}
  std::string name() const override { return "pnc-nsi"; }
  std::string units() const override { return m_unit; }

private:
  const std::string m_unit{"iQw*e-11"};
};

} // namespace DiracOperator
