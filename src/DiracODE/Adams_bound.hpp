#pragma once
#include "Adams_coefs.hpp"
#include <utility>
#include <vector>
class DiracSpinor;
class Grid;

/*!
@brief
Solve Dirac equation. Functions defined in "<DiracODE/DiracODE.hpp>"
@details
  - v is local potential (e.g., v = v_dir + v_nuc)
  - H_mag is off-diagonal magnetic form factor (QED radiative potential); enter
as '{}' to not use.
  - alpha: \f$\alpha = \lambda\alpha_0\f$ is the effective value of
fine-structure constant
*/
namespace DiracODE {

//******************************************************************************
//! @brief Solves bound-state problem for local potential (en < 0)
/*! @details
\f[ (H_0 + v - \epsilon_a)F_a = 0\f]
en0 is initial energy guess (must be reasonably good).
log_eps: log10(eps); eps is convergence target for energy.
*/
void boundState(DiracSpinor &Fa, const double en0, const std::vector<double> &v,
                const std::vector<double> &H_mag, const double alpha,
                int log_eps = 14, const DiracSpinor *const VxFa = nullptr,
                const DiracSpinor *const Fa0 = nullptr, double zion = 1);

//! @brief For given energy en, solves (local) DE with correct boundary
//! conditions at the origin
void regularAtOrigin(DiracSpinor &Fa, const double en,
                     const std::vector<double> &v,
                     const std::vector<double> &H_mag, const double alpha);

//! @brief For given energy en, solves (local) DE with correct boundary
//! conditions at infinity
void regularAtInfinity(DiracSpinor &Fa, const double en,
                       const std::vector<double> &v,
                       const std::vector<double> &H_mag, const double alpha);

//! @brief For given energy en (en > 0), solves (local) DE for continuum state
//! (with energy normalisation).
/*! @details
ext_grid is an 'extended grid' (see Grid class); needed since we need to solve
to large r to enforce boundary conditions (especially for large en). r_asym0 is
initial guess for 'asymptotic region'. ext_grid must extend past r_asym0.
*/
void solveContinuum(DiracSpinor &Fa, const double en,
                    const std::vector<double> &v, const Grid &ext_grid,
                    const double r_asym0, const double alpha,
                    const DiracSpinor *const VxFa = nullptr,
                    const DiracSpinor *const Fa0 = nullptr);

//******************************************************************************

//! @brief Solves inhomogeneous Dirac equation
/*! @details
\f[ (H_0 + v -\epsilon_a)F_a = S \f]
with `source' term, S. Solves for \f$\psi_\kappa\f$ with angular momentum kappa.
en = \f$\epsilon\f$ is given. Note sign of S.
Uses Green's method (see Method documentation).
*/
DiracSpinor solve_inhomog(const int kappa, const double en,
                          const std::vector<double> &v,
                          const std::vector<double> &H_mag, const double alpha,
                          const DiracSpinor &source);

//! @brief Solves inhomogeneous Dirac equation
/*! @details
As above. Overload to accept/overwrite solution to Fa. kappa is taken from Fa.
*/
void solve_inhomog(DiracSpinor &Fa, const double en,
                   const std::vector<double> &v,
                   const std::vector<double> &H_mag, const double alpha,
                   const DiracSpinor &source);

//! @brief Solves inhomogeneous Dirac equation
/*! @details
As above. Overload to accept/overwrite solution to Fa.
All these routines solve also for Fzero, Finf, which are solutions to
homogeneous equation  (H-en)Fa = 0 [reg @ origin, and infinity, respectively].
  - The first two throw these solutions away, the third keeps them (in some
cases they can be re-used)
  - These Spinors are solved internally and over-written, they don't need to be
solved first (i.e., they are out parameters, not in/out parameters)
*/
void solve_inhomog(DiracSpinor &Fa, DiracSpinor &Fzero, DiracSpinor &Finf,
                   const double en, const std::vector<double> &v,
                   const std::vector<double> &H_mag, const double alpha,
                   const DiracSpinor &source);

namespace Adams {

//******************************************************************************
// Parameters used for Adams-Moulton mehtod:
namespace Param {
constexpr int AMO = 7;            // Adams-Moulton (between 5 and 8)
constexpr AdamsCoefs<AMO> AMcoef; // Adamns-Moulton coeficients [defined .h]
constexpr double cALR = 550;      // 'assymptotically large r [kinda..]' (=800)
constexpr int max_its = 99;       // Max # attempts at converging [sove bs] (30)
constexpr double lfrac_de = 0.12; // 'large' energy variations (0.1 => 10%)
constexpr int d_ctp = 4;          // Num points past ctp +/- d_ctp.

// order of the expansion coeficients in 'inwardAM'  (15 orig.)
constexpr int nx = 15;
// convergance for expansion in `inwardAM'' (10^-8)
constexpr double nx_eps = 1.e-12;
// # of outdir runs [finds first Param::num_loops*AMO+1 points (3)]
constexpr int num_loops = 1;

// weighting function for meshing in/out solutions
// nb: must be positive, but i may be negative (?) [ctp - d_ctp]
constexpr auto weight = [](const int i) {
  return 1.0 / static_cast<double>(i * i + 1);
};

static_assert(Param::AMO >= 5 && Param::AMO <= 8,
              "\nFAIL 8 in Adams: parameter AMO must be between 5 and 8\n");

} // namespace Param

//******************************************************************************
class DiracMatrix {
  // Notation:
  // df = af - bg
  // dg = -cf + dg
public:
  DiracMatrix(const Grid &in_grid, const std::vector<double> &in_v,
              const int in_k, const double in_en, const double in_alpha,
              const std::vector<double> &Hmag = {},
              const DiracSpinor *const VxFa = nullptr,
              const DiracSpinor *const iFa0 = nullptr, double zion = 1);

  const Grid *const pgr;
  const std::vector<double> *const v;
  const std::vector<double> *const Hmag;
  const DiracSpinor *const VxFa;
  const DiracSpinor *const Fa0;
  const double zion = 1.0;
  const int k;
  const double en, alpha, cc;

  // update a and d for off-diag additional potential (magnetic form-fac, QED)
  double a(std::size_t i) const;
  double b(std::size_t i) const;
  double c(std::size_t i) const;
  double d(std::size_t i) const;
  std::tuple<double, double, double, double> abcd(std::size_t i) const;
  double dfdu(const std::vector<double> &f, const std::vector<double> &g,
              std::size_t i) const;
  double dgdu(const std::vector<double> &f, const std::vector<double> &g,
              std::size_t i) const;
  // Note: these are UN-SCALED
  double dfdu_X(std::size_t i) const;
  double dgdu_X(std::size_t i) const;
};

struct TrackEnGuess {
  // Number of times there was too many/too few nodes:
  int count_toomany = 0;
  int count_toofew = 0;
  // Upper and lower energy window before correct # nodes
  double high_en = 0.0;
  double low_en = 0.0;
};

// -----------------------------------------------------------------------------
int findPracticalInfinity(const double en, const std::vector<double> &v,
                          const std::vector<double> &r, const double alr);

int findClassicalTurningPoint(const double en, const std::vector<double> &v,
                              const int pinf, const int d_ctp);

void trialDiracSolution(std::vector<double> &f, std::vector<double> &g,
                        std::vector<double> &dg, const double en, const int ka,
                        const std::vector<double> &v,
                        const std::vector<double> &H_mag, const Grid &gr,
                        const int ctp, const int d_ctp, const int pinf,
                        const double alpha,
                        const DiracSpinor *const VxFa = nullptr,
                        const DiracSpinor *const Fa0 = nullptr,
                        double zion = 1);

int countNodes(const std::vector<double> &f, const int maxi);

// c++17: could use structured binding, move in/outs to returns
void largeEnergyChange(double *en, TrackEnGuess *sofar, double frac_de,
                       bool toomany_nodes);

double smallEnergyChangePT(const double en, const double anorm,
                           const std::vector<double> &f,
                           const std::vector<double> &dg, const int ctp,
                           const int d_ctp, const double alpha,
                           const TrackEnGuess &sofar);

void outwardAM(std::vector<double> &f, std::vector<double> &g,
               const DiracMatrix &Hd, const int final);

void inwardAM(std::vector<double> &f, std::vector<double> &g,
              const DiracMatrix &Hd, const int ctp, const int pinf);

void adamsMoulton(std::vector<double> &f, std::vector<double> &g,
                  const DiracMatrix &Hd, const int ni, const int nf);

void joinInOutSolutions(std::vector<double> &f, std::vector<double> &g,
                        std::vector<double> &dg,
                        const std::vector<double> &f_in,
                        const std::vector<double> &g_in, const int ctp,
                        const int d_ctp, const int pinf);

} // namespace Adams
} // namespace DiracODE
