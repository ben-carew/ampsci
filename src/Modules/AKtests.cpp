#include "Modules/AKtests.hpp"
#include "DMionisation/AKF_akFunctions.hpp"
#include "IO/InputBlock.hpp"
#include "Physics/PhysConst_constants.hpp" // For GHz unit conversion
#include "Wavefunction/Wavefunction.hpp"

const double AUMEV = PhysConst::c / PhysConst::m_e_MeV;

namespace Module {

void AKtests(const IO::InputBlock &input, const Wavefunction &wf) {
  std::cout << "\nAK tests module:\n";

  const int max_L = input.get("max_L", 0);
  const int core_n = input.get("core_n", 1);
  const int core_k = input.get("core_k", -1);
  const bool sum = input.get("sum", false);

  const int Npoints = input.get("Npoints", 1);
  const double dE = input.get("dE", 0);
  const double dEa = input.get("dEa", 0);
  const double dEb = input.get("dEb", 0);
  const double deltaE = (dEb - dEa) / Npoints;
  auto E_array = AKF::LinVect(dEa, dEb, Npoints);

  const int n = input.get("q_n", 0);
  const double q_min = 0.001;
  const double q_max = 1000;
  auto q_array = AKF::LogVect(q_min, q_max, n);
  size_t q_size = q_array.size();
  size_t k_size = wf.core.size();
  std::vector<float> K_nk(q_size);
  std::vector<float> K_sum(q_size);
  std::vector<float> K_bsum(q_size);
  std::vector<float> K_basis(q_size);
  std::vector<float> K_prime(q_size);
  const auto &r = wf.rgrid->r();

  const auto jLqr_f = AKF::sphericalBesselTable(max_L, q_array, r);

  if (sum == true) {
    // sum over all orbitals
    for (std::size_t j = 0; j < k_size; ++j) {
      for (std::size_t k = 0; k < E_array.size(); k++) {
        K_nk = AKF::calculateK_nk(wf, j, max_L, E_array[k], jLqr_f, q_size);
        for (std::size_t i = 0; i < q_size; i++) {
          K_prime[i] = K_prime[i] + K_nk[i] * deltaE;
        }
      }
      K_basis = AKF::basisK_nk(wf, j, max_L, dEa, dEb, jLqr_f, q_size);
      for (std::size_t i = 0; i < q_size; i++) {
        K_sum[i] = K_prime[i] + K_sum[i];
        K_bsum[i] = K_basis[i] + K_bsum[i];
      }
      std::fill(K_prime.begin(), K_prime.end(), 0);
    }
  } else {
    // calculate single orbitals
    for (std::size_t j = 0; j < k_size; ++j) {
      if (wf.core[j].n == core_n && wf.core[j].k == core_k) {
        for (std::size_t k = 0; k < E_array.size(); k++) {
          K_nk = AKF::calculateK_nk(wf, j, max_L, E_array[k], jLqr_f, q_size);
          for (std::size_t i = 0; i < q_size; i++) {
            K_prime[i] = K_prime[i] + K_nk[i] * deltaE;
          }
        }
        K_basis = AKF::basisK_nk(wf, j, max_L, dEa, dEb, jLqr_f, q_size);
      }
      // std::ofstream core_dat;
      // core_dat.open("j" + std::to_string(j) + ".txt");
      // for (std::size_t i = 0; i < q_array.size(); i++) {
      //   core_dat << q_array[i] / AUMEV << ' ' << K_prime[i] << ' ' <<
      //   K_basis[i]
      //            << '\n';
      // }
      // core_dat.close();
      // std::fill(K_nk.begin(), K_nk.end(), 0);
      // std::fill(K_prime.begin(), K_prime.end(), 0);
      // std::fill(K_basis.begin(), K_basis.end(), 0);
    }
  }

  // data ouput
  std::ofstream data1;
  data1.open("K_nk.txt");
  std::ofstream data2;
  data2.open("K_basis.txt");
  std::ofstream data3;
  data3.open("K_all.txt");
  if (sum == false) {
    for (std::size_t j = 0; j < q_array.size(); j++) {
      data1 << q_array[j] / AUMEV << ' ' << K_nk[j] << '\n';
      data2 << q_array[j] / AUMEV << ' ' << K_basis[j] << '\n';
      data3 << q_array[j] / AUMEV << ' ' << K_prime[j] << ' ' << K_basis[j]
            << '\n';
    }
  } else {
    for (std::size_t j = 0; j < q_array.size(); j++) {
      data1 << q_array[j] / AUMEV << ' ' << K_sum[j] << '\n';
      data2 << q_array[j] / AUMEV << ' ' << K_bsum[j] << '\n';
      data3 << q_array[j] / AUMEV << ' ' << K_sum[j] << ' ' << K_bsum[j]
            << '\n';
    }
  }
  data1.close();
  data2.close();
  data3.close();

  // // loop through the core states:
  // std::cout << "Core: \n";
  // for (std::size_t j = 0; j < k_size; ++j) {
  //   std::cout << wf.core[j].symbol() << " " << wf.core[j].en() << "\n";
  // }
  // // loop through the basis states :
  // std::cout << "Basis: \n";
  // for (const auto &Fb : wf.basis) {
  //   std::cout << Fb.symbol() << " " << Fb.en() << "\n";
  // }
} // AKtests
} // namespace Module
