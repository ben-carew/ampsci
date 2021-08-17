#include "Modules/ladder.hpp"
#include "Coulomb/Coulomb.hpp"
#include "Coulomb/QkTable.hpp"
#include "Coulomb/YkTable.hpp"
#include "IO/InputBlock.hpp"
#include "Wavefunction/Wavefunction.hpp"

namespace Module {

// Module for testing ladder diagram implementation
void ladder(const IO::InputBlock &input, const Wavefunction &wf) {
  std::cout << "\nLadder Module:\n\n";

  input.checkBlock2({{"option1", "example option 1 (double)"},
                     {"option2", "example option 2 (int)"}});

  // Example for retrieving input options:
  const auto val1 = input.get("option1", 0.0);
  const auto val2 = input.get("option2", 0);
  std::cout << "option1 = " << val1 << "\n";
  std::cout << "option2 = " << val2 << "\n";

  // "YkTable" calculates the Y^k_ab(r) "Hartree Screening"
  // functions These save much time when calculating Q^k
  // coeficients, since many different Qk coefs will re-use the same
  // y(r). Note: this is only used in 'qk.fill()'
  const Coulomb::YkTable yk(wf.rgrid, &wf.basis);

  // COnstruct an empty "Qk" table

  // QkTable assumes the "Qk" symmetry: {abcd} = cbad = adcb = cdab =
  // badc = bcda = dabc = dcba.
  Coulomb::QkTable qk;

  // WkTable assumes the "Wk" symmetry: {abcd} = badc = cdab = dcba
  // Coulomb::WkTable qk;

  // NkTable assumes no symmetry. Might need to use this to store L^k
  // Coulomb::NkTable qk;

  // Fill the Qk table will all possible (non-zero) Coulomb integrals
  std::cout << "Fill Qk table:\n";
  qk.fill(yk);

  // As an example, print out some Q^k values for core states:
  for (const auto &Fi : wf.core) {
    for (const auto &Fj : wf.core) {
      const auto [kmin, kmax] = Coulomb::k_minmax_Q(Fi, Fj, Fi, Fj);
      for (int k = kmin; k <= kmax; ++k) {
        // Print to screen. Note useful, just an example for how to use
        std::cout << "Q(" << k << "|" << Fi.shortSymbol() << ","
                  << Fj.shortSymbol() << "," << Fi.shortSymbol() << ","
                  << Fj.shortSymbol() << ") = " << qk.Q(k, Fi, Fj, Fi, Fj)
                  << "\n";
        // also have 'P' and 'W':
        // qk.P(k, Fi, Fj, Fi, Fj);
        // qk.W(k, Fi, Fj, Fi, Fj);
        // note: k_minmax_Q given min/max k for Q only
        // use: k_minmax_P to find non-zero k for 'P'
        // use: k_minmax_W to find non-zero k for 'W'
        // Get energies like: Fi.en()
        // Some more examples you might need:
        std::cout << "Energy: " << Fj.en() << "\n";
        std::cout << "n: " << Fj.n << "\n";
        std::cout << "kappa: " << Fj.k << "\n";
        std::cout << "2*j: " << Fj.twoj() << "\n";
      }
    }
  }

  // wf.core; // vector of core (Hartree Fock) states
  // wf.valence; // vector of valence (Hartree Fock) states
  // wf.basis; // vector of 'basis' states; include core and excited sates
  // together

  // To calculate Sigma, need to separate the core/excited basis states.
  // Either, use wf.isInCore(Fj) (returns true if 'Fj' is in the core)
  // Or, use wf.en_coreval_gap(), which:
  // => Fj.en < wf.en_coreval_gap() means Fj is in the core
  // May be easier to separate the basis into two new vectors, one for core, one
  // for excited.
}

} // namespace Module
