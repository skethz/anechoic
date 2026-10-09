// Checks that the engine's flag-built tables (tables.hpp::build_tables, the V80 host code used by sca_gpu_bias and
// check_ref_bias) equal mb_common.hpp::tables (scale = 1), which the ReAIM study's hw_check uses, for every schedule given
// on stdin as lines "t0 t1 S q lam jv kappa ramp". Prints one JSON line; exits 1 on any difference.
#include "mb_common.hpp"
#include "tables.hpp"
#include <iostream>
int main() {
    double t0, t1, q, lam, jv, kappa; int S, ramp; long n = 0, bad = 0;
    while (std::cin >> t0 >> t1 >> S >> q >> lam >> jv >> kappa >> ramp) {
        mb::Sched a; a.t0 = t0; a.t1 = t1; a.S = S; a.q = q; a.lam = lam; a.jv = jv; a.kappa = kappa; a.ramp = ramp != 0;
        v80host::Schedule b; b.t0 = t0; b.t1 = t1; b.S = S; b.q = q; b.lam = lam; b.jv = jv; b.kappa = kappa; b.ramp = ramp != 0;
        const sca::Tables x = mb::tables(a), y = v80host::build_tables(b);
        ++n; bad += !(x.fourT == y.fourT && x.q == y.q && x.kcorr == y.kcorr && x.kconst == y.kconst);
    }
    std::cout << "{\"schedules\":" << n << ",\"different\":" << bad << "}\n";
    return bad ? 1 : 0;
}
