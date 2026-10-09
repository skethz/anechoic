// Independent CPU checker for the bias kernel: re-runs sca::run_trial_bias(n, dense, bias, tables, seed, trial, trace = true)
// (fpga/v80_sca/src/sca_ref_bias.hpp, unchanged; -DSCA_LANES=256) for every trial written by sca_gpu_bias and compares final
// spins, sum_i s_i h_i (device), total flips and every per-step n_lin value bit for bit; also checks the trial's sum_bs and
// energy fields against the reference (E = -(sum_sh + sum_bs) / 2).
// Usage: check_ref_bias (--graph K2000.bin | --gset FILE | --jint8 FILE) [--bias FILE] --prefix OUT [schedule flags]
//                       [--seed S] [--threads T]
#include "sca_ref_bias.hpp"
#include "tables.hpp"
#include "jmat.hpp"
#include <atomic>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <thread>
#include <vector>

static int64_t field_of(const std::string& line, const std::string& key) {
    const std::string k = "\"" + key + "\":";
    const size_t p = line.find(k);
    if (p == std::string::npos) throw std::runtime_error("missing " + key);
    return std::stoll(line.substr(p + k.size()));
}

int main(int argc, char** argv) {
    try {
        v80host::Schedule sc; std::string jkind = "graph", jpath = "data/K2000.bin", bpath, prefix; uint64_t seed = 20261004; int nth = 32;
        for (int a = 1; a < argc; ++a) {
            std::string k = argv[a];
            if (v80host::parse_schedule_flag(sc, k, a, argc, argv)) continue;
            auto val = [&]() { if (a + 1 >= argc) throw std::runtime_error("missing value for " + k); return std::string(argv[++a]); };
            if (k == "--graph") { jkind = "graph"; jpath = val(); }
            else if (k == "--gset") { jkind = "gset"; jpath = val(); }
            else if (k == "--jint8") { jkind = "jint8"; jpath = val(); }
            else if (k == "--bias") bpath = val();
            else if (k == "--prefix") prefix = val(); else if (k == "--seed") seed = std::stoull(val());
            else if (k == "--threads") nth = std::stoi(val()); else throw std::runtime_error("unknown argument " + k);
        }
        const jmat::IntJ J = jmat::load(jkind, jpath);
        const std::vector<int32_t> bias = bpath.empty() ? std::vector<int32_t>() : jmat::load_bias(bpath, J.N);
        const int32_t* bptr = bpath.empty() ? nullptr : bias.data();
        const sca::Tables tb = v80host::build_tables(sc);
        const int S = sc.S;
        std::ifstream tr(prefix + ".trials.jsonl"), sp(prefix + ".spins.bin", std::ios::binary), tc(prefix + ".trace.bin", std::ios::binary);
        if (!tr || !sp || !tc) throw std::runtime_error("missing GPU output files (need trace run)");
        struct Rec { uint32_t id; int64_t sh, sbs, e, edev, flips; uint64_t spins[32]; std::vector<int32_t> trace; };
        std::vector<Rec> recs; std::string line;
        while (std::getline(tr, line)) {
            Rec r; r.id = uint32_t(field_of(line, "trial_id")); r.sh = field_of(line, "sum_sh"); r.sbs = field_of(line, "sum_bs");
            r.e = field_of(line, "energy"); r.edev = field_of(line, "device_energy"); r.flips = field_of(line, "flips");
            sp.read(reinterpret_cast<char*>(r.spins), sizeof r.spins);
            r.trace.resize(S); tc.read(reinterpret_cast<char*>(r.trace.data()), size_t(S) * 4);
            if (!sp || !tc) throw std::runtime_error("short binary file");
            recs.push_back(std::move(r));
        }
        std::atomic<int> next{0}, bad{0};
        std::vector<std::string> notes(recs.size());
        auto work = [&]() {
            for (int i; (i = next++) < int(recs.size());) {
                const Rec& r = recs[i];
                const sca::ResultBias ref = sca::run_trial_bias(J.N, J.dense.data(), bptr, tb, seed, r.id, true);
                const int64_t eref = -(ref.sum_sh + ref.sum_bs) / 2;
                bool ok = ref.sum_sh == r.sh && ref.sum_bs == r.sbs && eref == r.e && eref == r.edev && ref.flips == r.flips;
                for (int w = 0; w < sca::WORDS; ++w) ok &= ref.spins[w] == r.spins[w];
                int first_bad = -1;
                for (int t = 0; t < S; ++t) if (ref.n_lin[t] != r.trace[t]) { ok = false; if (first_bad < 0) first_bad = t; }
                if (!ok) {
                    ++bad;
                    std::ostringstream o; o << "trial " << r.id << ": ref sum_sh " << ref.sum_sh << " gpu " << r.sh << ", ref E " << eref << " gpu " << r.e << "/" << r.edev
                                            << ", ref flips " << ref.flips << " gpu " << r.flips << ", first n_lin mismatch step " << first_bad;
                    notes[i] = o.str();
                }
            }
        };
        std::vector<std::thread> th; for (int i = 0; i < nth; ++i) th.emplace_back(work); for (auto& t : th) t.join();
        int shown = 0;
        for (auto& n : notes) if (!n.empty() && shown++ < 20) std::cerr << n << "\n";
        std::cout << "{\"prefix\":\"" << prefix << "\",\"instance\":\"" << jpath << "\",\"bias\":\"" << bpath << "\",\"N\":" << J.N << ",\"trials\":" << recs.size()
                  << ",\"exact\":" << int(recs.size()) - bad.load() << ",\"mismatch\":" << bad.load() << ",\"steps\":" << S << ",\"seed\":" << seed << "}\n";
        return bad.load() ? 2 : 0;
    } catch (const std::exception& e) { std::cerr << "FAIL: " << e.what() << "\n"; return 1; }
}
