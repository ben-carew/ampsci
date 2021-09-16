#pragma once
#include <string>
#include <vector>
class Wavefunction;

namespace AKF {

double CLkk(int L, int ka, int kb);

void writeToTextFile(const std::string &fname,
                     const std::vector<std::vector<std::vector<float>>> &AK,
                     const std::vector<std::string> &nklst, double qmin,
                     double qmax, double demin, double demax);

int akReadWrite(const std::string &fname, bool write,
                std::vector<std::vector<std::vector<float>>> &AK,
                std::vector<std::string> &nklst, double &qmin, double &qmax,
                double &dEmin, double &dEmax);

std::vector<float>
calculateK_nk(const Wavefunction &wf, std::size_t nk, int max_L, double dE,
              const std::vector<std::vector<std::vector<double>>> &jLqr_f,
              std::size_t q_size);

void calculateK_nk(const Wavefunction &wf, std::size_t nk, int max_L, double dE,
                   const std::vector<std::vector<std::vector<double>>> &jLqr_f,
                   std::vector<float> &K_nk);

std::vector<float>
basisK_nk(const Wavefunction &wf, std::size_t is, int max_L, double dEa,
          double dEb,
          const std::vector<std::vector<std::vector<double>>> &jLqr_f,
          std::size_t q_size);

int calculateKpw_nk(const Wavefunction &wf, std::size_t nk, double dE,
                    std::vector<std::vector<double>> &jl_qr,
                    std::vector<float> &K_nk);

std::vector<std::vector<std::vector<double>>>
sphericalBesselTable(int max_L, const std::vector<double> &q_array,
                     const std::vector<double> &r);
void sphericalBesselTable(std::vector<std::vector<std::vector<double>>> &jLqr_f,
                          int max_L, const std::vector<double> &q_array,
                          const std::vector<double> &r);

std::vector<double> LogVect(double min, double max, int num_points);
std::vector<double> LinVect(double min, double max, int n);

} // namespace AKF
