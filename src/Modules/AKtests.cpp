#include "Modules/AKtests.hpp"
#include "DMionisation/AKF_akFunctions.hpp"
#include "IO/InputBlock.hpp"
#include "Physics/PhysConst_constants.hpp" // For GHz unit conversion
#include "Wavefunction/Wavefunction.hpp"

const double AUMEV = PhysConst::c / PhysConst::m_e_MeV;

namespace Module {

void AKtests(const IO::InputBlock &input, const Wavefunction &wf) {
  std::cout << "\nAK tests module:\n\n";

  std::ofstream myfile;
  myfile.open("K_nk.txt");

  // // Some test outputs for practice
  // int value = input.get("testinput",0);
  // std::vector<int> v;
  // v.push_back(1);
  // AKF::addThirty(v);
  // for (int i=0; i<v.size(); i++){
  //   std::cout << "v["<<i<<"] = " << v[i] << ", ";}
  //   std::cout << "size = " << v.size() << '\n';
  // std::cout << "value = " << value << '\n';

  const int max_L = input.get("max_L", 0);
  const double dE = input.get("dE", 0);

  const size_t n = 100;
  const int q_min = 1;
  const int q_max = 1000;
  auto q_array = AKF::LogVect(q_min, q_max, n);
  std::vector<float> K_nk(q_array.size());

  // const std::vector<double> r{1,2,3,4,5};

  // int calculateK_nk(const Wavefunction &wf, std::size_t is, int max_L,
  //                   double dE,
  //                   std::vector<std::vector<std::vector<double>>> &jLqr_f,
  //                   std::vector<float> &K_nk, double Zeff = -1);

  // void sphericalBesselTable(
  //     std::vector<std::vector<std::vector<double>>> & jLqr_f, int max_L,
  //     const std::vector<double> &q_array, const std::vector<double> &r);
  const auto &r = wf.rgrid->r();

  const auto jLqr_f = AKF::sphericalBesselTable(max_L, q_array, r);

  // const auto &Jlq_r = jLqr_f[iL][iq];

  // loop through the core states:
  std::cout << "Core: \n";
  // for (const auto &Fc : wf.core) {
  //   std::cout << Fc.symbol() << " " << Fc.en() << "\n";
  // }

  for (std::size_t i = 0; i < wf.core.size(); ++i) {
    std::cout << wf.core[i].symbol() << " " << wf.core[i].en() << "\n";
    if (wf.core[i].n == 3 && wf.core[i].k == -1) {
      AKF::calculateK_nk(wf, i, max_L, dE, jLqr_f, K_nk);
    }
  }

  // std::cout << "Atomic Kernel K_nk" << '\n';
  for (std::size_t i = 0; i < q_array.size(); i++) {
    myfile << q_array[i] / AUMEV << ' ' << K_nk[i] << '\n';
  }
  myfile.close();

  // loop through the basis states:
  // std::cout << "Basis: \n";
  // for (const auto &Fb : wf.basis) {
  //   std::cout << Fb.symbol() << " " << Fb.en() << "\n";
  // }
} // AKtests
} // namespace Module
