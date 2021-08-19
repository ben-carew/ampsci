#include "Modules/AKtests.hpp"
#include "DMionisation/AKF_akFunctions.hpp"
#include "IO/InputBlock.hpp"
#include "Physics/PhysConst_constants.hpp" // For GHz unit conversion
#include "Wavefunction/Wavefunction.hpp"

namespace Module {

void AKtests(const IO::InputBlock &input, const Wavefunction &wf) {

  // int calculateK_nk(const Wavefunction &wf, std::size_t nk, int max_L,
  //                   double dE,
  //                   std::vector<std::vector<std::vector<double>>> &jLqr_f,
  //                   std::vector<float> &K_nk, double Zeff = -1);

  // void sphericalBesselTable(
  //     std::vector<std::vector<std::vector<double>>> & jLqr_f, int max_L,
  //     const std::vector<double> &q_array, const std::vector<double> &r);

  // loop through the core states:
  std::cout << "Core: \n";
  for (const auto &Fc : wf.core) {
    std::cout << Fc.symbol() << " " << Fc.en() << "\n";
  }

  // loop through the basis states:
  std::cout << "Basis: \n";
  for (const auto &Fb : wf.basis) {
    std::cout << Fb.symbol() << " " << Fb.en() << "\n";
  }
}

} // namespace Module
