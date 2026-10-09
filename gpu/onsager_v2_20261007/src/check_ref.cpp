// Independent CPU checker: re-runs sca::run_trial(trace=true) (sca_ref.hpp, -DSCA_LANES=256) for every trial written by sca_gpu
// and compares final spins, cut, total flips and every per-step n_lin value bit for bit.
// Usage: check_ref --graph data/K2000.bin --prefix OUT [schedule flags as given to sca_gpu] [--seed S] [--threads T]
#include "sca_ref.hpp"
#include "tables.hpp"
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
        v80host::Schedule sc; std::string graph = "data/K2000.bin", prefix; uint64_t seed = 20261004; int nth = 32;
        for (int a = 1; a < argc; ++a) {
            std::string k = argv[a];
            if (v80host::parse_schedule_flag(sc, k, a, argc, argv)) continue;
            auto val = [&]() { if (a + 1 >= argc) throw std::runtime_error("missing value for " + k); return std::string(argv[++a]); };
            if (k == "--graph") graph = val(); else if (k == "--prefix") prefix = val(); else if (k == "--seed") seed = std::stoull(val());
            else if (k == "--threads") nth = std::stoi(val()); else throw std::runtime_error("unknown argument " + k);
        }
        const sca::Graph g = sca::load_graph(graph);
        if (g.sumw != -1040) throw std::runtime_error("unexpected graph");
        const sca::Tables tb = v80host::build_tables(sc);
        const int S = sc.S;
        std::ifstream tr(prefix + ".trials.jsonl"), sp(prefix + ".spins.bin", std::ios::binary), tc(prefix + ".trace.bin", std::ios::binary);
        if (!tr || !sp || !tc) throw std::runtime_error("missing GPU output files (need trace run)");
        struct Rec { uint32_t id; int64_t cut, dcut, flips; uint64_t spins[32]; std::vector<int32_t> trace; };
        std::vector<Rec> recs; std::string line;
        while (std::getline(tr, line)) {
            Rec r; r.id = uint32_t(field_of(line, "trial_id")); r.cut = field_of(line, "cut"); r.dcut = field_of(line, "device_cut");
            r.flips = field_of(line, "flips");
            sp.read(reinterpret_cast<char*>(r.spins), sizeof r.spins);
            r.trace.resize(S); tc.read(reinterpret_cast<char*>(r.trace.data()), size_t(S) * 4);
            if (!sp || !tc) throw std::runtime_error("short binary file");
            recs.push_back(std::move(r));
        }
        const std::vector<int8_t> dense = sca::dense_matrix(g);
        std::atomic<int> next{0}, bad{0};
        std::vector<std::string> notes(recs.size());
        auto work = [&]() {
            for (int i; (i = next++) < int(recs.size());) {
                const Rec& r = recs[i];
                const sca::Result ref = sca::run_trial(g, tb, seed, r.id, true, dense.data());
                bool ok = ref.cut == r.cut && ref.cut == r.dcut && ref.flips == r.flips;
                for (int w = 0; w < sca::WORDS; ++w) ok &= ref.spins[w] == r.spins[w];
                int first_bad = -1;
                for (int t = 0; t < S; ++t) if (ref.n_lin[t] != r.trace[t]) { ok = false; if (first_bad < 0) first_bad = t; }
                if (!ok) {
                    ++bad;
                    std::ostringstream o; o << "trial " << r.id << ": ref cut " << ref.cut << " gpu " << r.cut << "/" << r.dcut << ", ref flips " << ref.flips
                                            << " gpu " << r.flips << ", first n_lin mismatch step " << first_bad;
                    notes[i] = o.str();
                }
            }
        };
        std::vector<std::thread> th; for (int i = 0; i < nth; ++i) th.emplace_back(work); for (auto& t : th) t.join();
        for (auto& n : notes) if (!n.empty()) std::cerr << n << "\n";
        std::cout << "{\"prefix\":\"" << prefix << "\",\"trials\":" << recs.size() << ",\"exact\":" << int(recs.size()) - bad.load()
                  << ",\"mismatch\":" << bad.load() << ",\"steps\":" << S << ",\"seed\":" << seed << "}\n";
        return bad.load() ? 2 : 0;
    } catch (const std::exception& e) { std::cerr << "FAIL: " << e.what() << "\n"; return 1; }
}
