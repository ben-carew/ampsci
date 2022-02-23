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
  const bool sum = input.get("sum", true);
  bool calc_j = input.get("calc_j", true);
  const bool do_finite_diff = input.get("do_finite_diff", false);
  const bool do_basis = input.get("do_basis", true);

  const int Npoints = input.get("Npoints", 1);
  const double dE = input.get("dE", 0);
  const double dEa = input.get("dEa", 0);
  const double dEb = input.get("dEb", 0);
  const double deltaE = (dEb - dEa) / Npoints;
  auto E_array = AKF::LinVect(dEa, dEb, Npoints);

  const int n = input.get("q_n", 200);
  const double q_min = 0.001;
  const double q_max = 1000;
  auto q_array = AKF::LogVect(q_min, q_max, n);
  const size_t q_size = q_array.size();
  const size_t k_size = wf.core.size();
  const auto &r = wf.rgrid->r();
  const size_t r_size = r.size();
  std::vector<float> K_nk(q_size);
  std::vector<float> K_sum(q_size);
  std::vector<float> K_bsum(q_size);
  std::vector<float> K_basis(q_size);
  std::vector<float> K_prime(q_size);

  std::vector<std::vector<std::vector<double>>> jLqr_f;

  if (calc_j == true) {
    jLqr_f = AKF::sphericalBesselTable(max_L, q_array, r);
    // std::ofstream test1;
    // test1.open("jLqr_before.txt");
    // for (int iq = 0; iq < q_size; iq++) {
    //   for (int ir = 0; ir < r_size; ir++) {
    //     test1 << ir * iq << ' ' << jLqr_f[0][iq][ir] << '\n';
    //   }
    // }
  } else {
    int L, iq, ir;
    char delim = ' ';
    jLqr_f.resize(max_L + 1, std::vector<std::vector<double>>(
                                 q_size, std::vector<double>(r_size)));
    std::vector<std::string> j_read;
    std::string line;
    std::ifstream data_in("jLqr.txt");
    if (data_in.is_open()) {
      for (int l = 1; l <= (max_L + 1) * q_size * r_size; l++) {
        getline(data_in, line, delim);
        L = std::stoi(line);
        getline(data_in, line, delim);
        iq = std::stoi(line);
        getline(data_in, line, delim);
        ir = std::stoi(line);
        getline(data_in, line);
        jLqr_f[L][iq][ir] = std::stod(line);
        /*std::cout << l << ' ' << L << ' ' << iq << ' ' << ir << ' '
                  << jLqr_f[L][iq][ir] << '\n';*/
      }
    }
    // std::ofstream test2;
    // test2.open("jLqr_after.txt");
    // for (int iq = 0; iq < q_size; iq++) {
    //   for (int ir = 0; ir < r_size; ir++) {
    //     test2 << ir * iq << ' ' << jLqr_f[0][iq][ir] << '\n';
    //   }
    // }
  }

  if (sum == true) {
    // sum over all orbitals
    for (std::size_t j = 0; j < k_size; ++j) {
      if (do_finite_diff == true) {
        for (std::size_t k = 0; k < E_array.size(); k++) {
          K_nk = AKF::calculateK_nk(wf, j, max_L, E_array[k], jLqr_f, q_size);
          for (std::size_t i = 0; i < q_size; i++) {
            K_prime[i] = K_prime[i] + K_nk[i] * deltaE;
          }
        }
        for (std::size_t i = 0; i < q_size; i++) {
          K_sum[i] = K_prime[i] + K_sum[i];
        }
      }
      if (do_basis == true) {
        K_basis = AKF::basisK_nk(wf, j, max_L, dEa, dEb, jLqr_f, q_size);
        for (std::size_t i = 0; i < q_size; i++) {
          K_bsum[i] = K_basis[i] + K_bsum[i];
        }
        std::fill(K_prime.begin(), K_prime.end(), 0);
      }
    }
  } else {
    // calculate single orbitals
    for (std::size_t j = 0; j < k_size; ++j) {
      if (wf.core[j].n == core_n && wf.core[j].k == core_k) {
        if (do_finite_diff == true) {
          for (std::size_t k = 0; k < E_array.size(); k++) {
            K_nk = AKF::calculateK_nk(wf, j, max_L, E_array[k], jLqr_f, q_size);
            for (std::size_t i = 0; i < q_size; i++) {
              K_prime[i] = K_prime[i] + K_nk[i] * deltaE;
            }
          }
        }
        if (do_basis == true) {
          K_basis = AKF::basisK_nk(wf, j, max_L, dEa, dEb, jLqr_f, q_size);
        }
      }
      // // Calculate all orbitals individually
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
  std::ofstream data;
  data.open("K_all.txt");
  if (sum == false) {
    for (std::size_t j = 0; j < q_array.size(); j++) {
      data << q_array[j] / AUMEV << ' ' << K_prime[j] << ' ' << K_basis[j]
           << '\n';
    }
  } else {
    for (std::size_t j = 0; j < q_array.size(); j++) {
      data << q_array[j] / AUMEV << ' ' << K_sum[j] << ' ' << K_bsum[j] << '\n';
    }
  }
  data.close();

} // AKtests
} // namespace Module
