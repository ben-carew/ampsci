#include "Modules/AKtests.hpp"
#include "DMionisation/AKF_akFunctions.hpp"
#include "IO/InputBlock.hpp"
#include "Maths/Interpolator.hpp"
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
  const bool calc_j = input.get("calc_j", true);
  const bool do_finite_diff = input.get("do_finite_diff", false);
  const bool do_basis = input.get("do_basis", true);
  const bool energy_dependent = input.get("energy_dependent", true);
  const bool dmX = input.get("dmX", true);

  const int Npoints = input.get("Npoints", 1);
  [[maybe_unused]] const double dE = input.get("dE", 0);
  double dEa = input.get("dEa", 0);
  double dEb = input.get("dEb", 0);
  const double deltaE = (dEb - dEa) / Npoints;
  auto E_array = AKF::LinVect(dEa, dEb, Npoints);

  const int n = input.get("q_n", 200);
  double q_min = 0.001;
  double q_max = 1000;
  auto q_array = AKF::LogVect(q_min, q_max, n);
  const size_t q_size = q_array.size();
  const size_t k_size = wf.core.size();
  const auto &r = wf.rgrid->r();
  const size_t r_size = r.size();
  std::vector<double> K_nk(q_size);
  std::vector<double> K_sum(q_size);
  std::vector<double> K_bsum(q_size);
  std::vector<double> K_basis(q_size);
  std::vector<double> K_basis_single(q_size);
  std::vector<double> K_prime(q_size);
  std::vector<double> K_en_dep(q_size);
  std::vector<double> K_en_single(q_size);
  std::vector<double> K_en_sum(q_size);
  std::vector<std::vector<double>> K_energy;
  std::vector<std::vector<double>> K_E_sum(q_size,
                                           std::vector<double>(Npoints));
  std::vector<std::vector<std::vector<double>>> K_enq(
      Npoints,
      std::vector<std::vector<double>>(k_size, std::vector<double>(q_size)));

  std::vector<std::vector<std::vector<double>>> jLqr_f;

  if (calc_j == true) {
    jLqr_f = AKF::sphericalBesselTable(max_L, q_array, r);
  } else {
    int L, iq, ir;
    char delim = ' ';
    jLqr_f.resize(max_L + 1, std::vector<std::vector<double>>(
                                 q_size, std::vector<double>(r_size)));
    std::string line;
    std::string L_check;
    std::string q_check;
    std::string r_check;
    std::ifstream data_in;
    data_in.open("jLqr.txt");
    if (data_in.is_open()) {
      getline(data_in, L_check, delim);
      getline(data_in, q_check, delim);
      getline(data_in, r_check, delim);
      if (max_L == std::stoi(L_check) && (int)q_size == std::stoi(q_check) &&
          (int)r_size == std::stoi(r_check)) {
        std::cout << "Reading jLqr..." << std::endl;
        for (int l = 1; l < (max_L + 1) * (int)(q_size) * (int)(r_size); l++) {
          getline(data_in, line, delim);
          L = std::stoi(line);
          getline(data_in, line, delim);
          iq = std::stoi(line);
          getline(data_in, line, delim);
          ir = std::stoi(line);
          getline(data_in, line);
          jLqr_f[L][iq][ir] = std::stod(line);
          // std::cout << l << ' ' << L << ' ' << iq << ' ' << ir << ' '
          //           << jLqr_f[L][iq][ir] << '\n';
        }
        std::cout << "jLqr read successful." << std::endl;
      } else {
        std::cout << "Error: jLqr is of wrong size, need to recalculate."
                  << std::endl;
        exit(1);
      }
    }
    data_in.close();
  }

  if (sum == true) {
    // sum over all orbitals
    for (std::size_t j = 0; j < k_size; ++j) {
      if (do_finite_diff == true) {
        std::cout << "Calculating finite difference sum" << std::endl;
        for (std::size_t k = 0; k < E_array.size(); k++) {
          K_nk = AKF::calculateK_nk(wf, j, max_L, E_array[k], jLqr_f, q_size);
          for (std::size_t i = 0; i < q_size; i++) {
            K_prime[i] = K_prime[i] + K_nk[i] * deltaE;
            // std::cout << K_prime.at(i) << std::endl;
          }
        }
        for (std::size_t i = 0; i < q_size; i++) {
          K_sum[i] = K_prime[i] + K_sum[i];
        }
      }
      if (do_basis == true) {
        std::cout << "Calculating with basis states" << std::endl;
        K_basis = AKF::basisK_nk(wf, j, max_L, dEa, dEb, jLqr_f, q_size);
        for (std::size_t i = 0; i < q_size; i++) {
          K_bsum[i] = K_basis[i] + K_bsum[i];
        }
      }
    }
    std::fill(K_prime.begin(), K_prime.end(), 0);
  } else {
    // calculate single orbitals
    for (std::size_t j = 0; j < k_size; ++j) {
      if (wf.core[j].n == core_n && wf.core[j].k == core_k) {
        if (do_finite_diff == true) {
          std::cout << "Calculating with finite difference" << std::endl;
          for (int k = 0; k < Npoints; k++) {
            K_nk =
                AKF::calculateK_nk(wf, j, max_L, E_array.at(k), jLqr_f, q_size);
            for (std::size_t i = 0; i < q_size; i++) {
              K_prime.at(i) += K_nk.at(i) * deltaE;
            }
          }
        }
        if (do_basis == true) {
          std::cout << "Calculating with basis states" << std::endl;
          K_basis_single =
              AKF::basisK_nk(wf, j, max_L, dEa, dEb, jLqr_f, q_size);
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

  if (energy_dependent == true) {
    // make K_basis dependent on E.
    std::cout << "Calculating with basis states" << std::endl;
    if (sum == true) {
      for (std::size_t j = 0; j < k_size; ++j) {
        // if (wf.core[j].n == core_n && wf.core[j].k == core_k) {
        K_energy = AKF::basisK_energy(wf, j, max_L, dEa, dEb, jLqr_f, q_array,
                                      Npoints);
        for (std::size_t iq = 0; iq < K_energy.size(); iq++) {
          for (std::size_t ie = 0; ie < K_energy.at(iq).size(); ie++) {
            K_E_sum.at(iq).at(ie) += K_energy.at(iq).at(ie);
            K_en_sum.at(iq) += K_E_sum.at(iq).at(ie);
          }
        }
      }
    } else {
      for (std::size_t j = 0; j < k_size; ++j) {
        if (wf.core[j].n == core_n && wf.core[j].k == core_k) {
          K_energy = AKF::basisK_energy(wf, j, max_L, dEa, dEb, jLqr_f, q_array,
                                        Npoints);
          for (std::size_t iq = 0; iq < q_size; iq++) {
            for (auto ie = 0; ie < Npoints; ie++) {
              K_en_single.at(iq) += K_energy.at(iq).at(ie) * std::abs(deltaE);
            }
          }
        }
      }
    }
  }

  if (dmX == true) {
    const std::string fname = "ak-Xe_test";
    std::vector<std::string> nklst(k_size);
    for (std::size_t j = 0; j < k_size; ++j) {
      nklst.at(j) = wf.core.at(j).symbol();
    }

    for (std::size_t j = 0; j < k_size; ++j) {
      K_energy =
          AKF::basisK_energy(wf, j, max_L, dEa, dEb, jLqr_f, q_array, Npoints);
      for (auto ie = 0; ie < Npoints; ie++) {
        for (std::size_t iq = 0; iq < q_size; iq++) {
          K_enq.at(ie).at(j).at(iq) = K_energy.at(iq).at(ie);
        }
      }
    }
    AKF::akReadWrite(fname, true, K_enq, nklst, q_min, q_max, dEa, dEb);
  }

  std::vector<double> K_wrt_E(Npoints, 0);
  for (int ie = 0; ie < Npoints; ie++) {
    for (std::size_t iq = 0; iq < 1; iq++) {
      //  K_wrt_E.at(ie) += K_E_sum.at(iq).at(ie) * q_array.at(iq);
      K_wrt_E.at(ie) +=
          K_energy.at(iq).at(ie) * q_array.at(iq) * q_array.at(iq);
    }
    // std::cout << K_wrt_E.at(ie) << std::endl;
  }

  std::ofstream data1;
  data1.open("K_Energy.txt");
  for (int j = 0; j < Npoints; j++) {
    data1 << E_array.at(j) << ' ' << std::abs(K_wrt_E.at(j)) << std::endl;
  }
  data1.close();

  // data ouput
  std::ofstream data2;
  data2.open("K_all.txt");
  for (std::size_t j = 0; j < q_array.size(); j++) {
    data2 << q_array[j] / AUMEV << ' ' << K_prime.at(j) << ' '
          << K_basis_single.at(j) << ' ' << std::abs(K_en_single.at(j)) << ' '
          << K_sum.at(j) << ' ' << K_bsum.at(j) << ' '
          << std::abs(K_en_sum.at(j)) << '\n';
  }
  data2.close();

  // // Calculate sum of differences between two methods
  // std::ofstream data2;
  // data2.open("difference.txt", std::ios_base::app);
  // const int bas_num = input.get("bas_num", 100);
  // double difference;
  // double sum_diff = 0;
  // double sum_K, sum_rel;
  // for (auto &vec_n : K_prime) {
  //   sum_K += vec_n;
  // }
  // for (std::size_t j = 0; j < q_size; j++) {
  //   difference = std::abs(K_prime[j] - K_basis[j]);
  //   sum_diff = sum_diff + difference;
  // }
  // sum_rel = sum_diff / sum_K;
  // data2 << bas_num << ' ' << sum_diff << ' ' << sum_rel << '\n';
  // data2.close();

  // Print basis energies and k values
  // for (size_t ii = 0; ii < wf.basis.size(); ii++) {
  //   std::cout << wf.basis.at(ii).k << " " << wf.basis.at(ii).en() <<
  //   std::endl;
  // }
} // namespace Module
} // namespace Module
