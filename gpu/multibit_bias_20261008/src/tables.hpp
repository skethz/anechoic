// Per-step tables exactly as the V80 host builds them.
// The body of build_tables() is copied verbatim from fpga/v80_sca/src/host_sca_multi_v6.cpp (lines ~117-127 and geom()).
#pragma once
#include "sca_ref.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

namespace v80host {

inline std::vector<double> geom(double t0, double t1, int S) {
    std::vector<double> T(S); for (int t = 0; t < S; ++t) T[t] = t0 * std::pow(t1 / t0, double(t) / std::max(1, S - 1)); return T;
}

struct Schedule {
    int S = 560;
    double t0 = 30, t1 = 5, q = 4, lam = 0, jv = 0, kappa = 0;
    bool ramp = false;
};

inline sca::Tables build_tables(const Schedule& c) {
    const int S = c.S; const double t0 = c.t0, t1 = c.t1, q = c.q, lam = c.lam, jv = c.jv, kappa = c.kappa; const bool ramp = c.ramp;
    // ---- verbatim from host_sca_multi_v6.cpp ----
        std::vector<double> L(S, lam);
        if (ramp) for (int t = 0; t < S; ++t) if (t >= 0.7 * S) L[t] = lam * (S - 1 - t) / std::max(1.0, S - 1 - 0.7 * S);
        sca::Tables tb = sca::make_tables(geom(t0, t1, S), std::vector<double>(S, q), L, jv);
        if (kappa > 0) {  // TEC-T (same as host_sca.cpp): kconst[t] = kappa*T_{t-1}*ramp(t) in Q16.16, no popcount term
            const std::vector<double> Tg = geom(t0, t1, S);
            for (int t = 1; t < S; ++t) {
                const double rf = (ramp && t >= 0.7 * S) ? (S - 1 - t) / std::max(1.0, S - 1 - 0.7 * S) : 1.0;
                tb.kconst[t] = int32_t(std::llround(kappa * Tg[t - 1] * rf * 65536.0));
                tb.kcorr[t] = 0;
            }
        }
    // ---- end verbatim ----
    return tb;
}

// Parse the V80 host's schedule flags (--steps --t0 --t1 --q --lambda --ramp --tec-jv --tecT). Returns true if consumed.
inline bool parse_schedule_flag(Schedule& c, const std::string& k, int& a, int argc, char** argv) {
    if (k == "--ramp") { c.ramp = true; return true; }
    auto val = [&]() { if (a + 1 >= argc) throw std::runtime_error("missing value for " + k); return std::string(argv[++a]); };
    if (k == "--steps") c.S = std::stoi(val());
    else if (k == "--t0") c.t0 = std::stod(val());
    else if (k == "--t1") c.t1 = std::stod(val());
    else if (k == "--q") c.q = std::stod(val());
    else if (k == "--lambda") c.lam = std::stod(val());
    else if (k == "--tec-jv") c.jv = std::stod(val());
    else if (k == "--tecT") c.kappa = std::stod(val());
    else return false;
    return true;
}

}  // namespace v80host
