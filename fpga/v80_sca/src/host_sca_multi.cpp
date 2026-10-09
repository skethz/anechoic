// Physical runner for the v5.2 engine interface (packed tables) with 1..5 replicated engines.
// --verify: seeded trials with per-step trace, each checked bit for bit against sca_ref.hpp (build with -DSCA_LANES=256).
// --cohort: timed trials. Every round starts all engines, then polls them; each engine's passive counter gives its own
// start-to-done cycles and the round time is their maximum. Every returned state is rescored independently on the host.
// Control windows: BAR0 offset CONTROL + e*0x1000. HBM (BAR2): J 0x0, packed tables 0x80000, out_e 0x100000 + e*0x40000,
// trace_e 0x400000 + e*0x200000.
#include "sca_ref.hpp"
#include "xsca_v80_hw.h"
#include <ami.h>
#include <ami_device.h>
#include <ami_mem_access.h>
#include <algorithm>
#include <chrono>
#include <cctype>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#ifndef SCA_TIMER_ID
#error "define SCA_TIMER_ID / SCA_CYC_LO / SCA_CYC_HI to match prepare_ip.py"
#endif

namespace {
constexpr uint64_t CONTROL = 0x01100000, HBM = 0x004000000000;
constexpr uint32_t J_OFF = 0x000000, TAB_OFF = 0x080000, OUT_BASE = 0x100000, OUT_STRIDE = 0x40000;
constexpr uint32_t TRACE_BASE = 0x400000, TRACE_STRIDE = 0x200000, BAR2_SIZE = 0x01000000;
constexpr uint32_t TIMER_MAGIC = 0x53434131;  // "SCA1"
constexpr int MAX_E = 5;

void check(int status, const char* what) {
    if (status != AMI_STATUS_OK) throw std::runtime_error(std::string(what) + ": " + ami_get_last_error());
}

struct Device {
    ami_device* h = nullptr;
    Device() { check(ami_dev_find("01:00.0", &h), "Find V80"); }
    ~Device() { if (h) ami_dev_delete(&h); }
    std::string uuid() {
        char v[AMI_LOGIC_UUID_SIZE] = {}; check(ami_dev_read_uuid(h, v), "Read UUID");
        std::string r(v); while (!r.empty() && std::isspace((unsigned char)r.back())) r.pop_back(); return r;
    }
    uint32_t rd(int e, uint32_t off) { uint32_t v = 0; check(ami_mem_bar_read(h, 0, CONTROL + e * 0x1000 + off, &v), "Read control"); return v; }
    void wr(int e, uint32_t off, uint32_t v) { check(ami_mem_bar_write(h, 0, CONTROL + e * 0x1000 + off, v), "Write control"); }
    void wr64(int e, uint32_t off, uint64_t v) { wr(e, off, uint32_t(v)); wr(e, off + 4, uint32_t(v >> 32)); }
    std::vector<uint32_t> read(uint32_t off, size_t bytes) {
        if (bytes % 4 || uint64_t(off) + bytes > BAR2_SIZE) throw std::runtime_error("HBM read range");
        std::vector<uint32_t> r(bytes / 4);
        for (size_t w = 0; w < r.size(); w += 4096)
            check(ami_mem_bar_read_range(h, 2, off + w * 4, uint32_t(std::min<size_t>(4096, r.size() - w)), r.data() + w), "Read HBM");
        return r;
    }
    void write(uint32_t off, const void* data, size_t bytes) {
        if (bytes % 4 || uint64_t(off) + bytes > BAR2_SIZE) throw std::runtime_error("HBM write range");
        std::vector<uint32_t> w(bytes / 4); std::memcpy(w.data(), data, bytes);
        for (size_t i = 0; i < w.size(); i += 4096)
            check(ami_mem_bar_write_range(h, 2, off + i * 4, uint32_t(std::min<size_t>(4096, w.size() - i)), w.data() + i), "Write HBM");
        if (read(off, bytes) != w) throw std::runtime_error("HBM staging readback mismatch");
    }
    uint64_t cycles(int e) { return uint64_t(rd(e, SCA_CYC_LO)) | (uint64_t(rd(e, SCA_CYC_HI)) << 32); }
};

std::vector<double> geom(double t0, double t1, int S) {
    std::vector<double> T(S); for (int t = 0; t < S; ++t) T[t] = t0 * std::pow(t1 / t0, double(t) / std::max(1, S - 1)); return T;
}

int64_t scalar_cut(const sca::Graph& g, const uint64_t* bits) {  // independent scorer
    int64_t c = 0;
    for (int x = 0; x < sca::N; ++x) {
        const bool sx = (bits[x >> 6] >> (x & 63)) & 1;
        for (int y = x + 1; y < sca::N; ++y)
            if (sx != (((bits[y >> 6] >> (y & 63)) & 1) != 0)) c += -sca::jval(g, x, y);
    }
    return c;
}
}  // namespace

int main(int argc, char** argv) {
    try {
        std::string action, uuid, output;
        int S = 560, trials = 1, per_launch = 1, E = 1;
        double t0 = 30, t1 = 5, q = 4, lam = 0, jv = 0, mhz = 0;
        bool ramp = false;
        uint64_t seed = 20261004; uint32_t offset = 0;
        for (int a = 1; a < argc; ++a) {
            std::string k = argv[a];
            if (k == "--probe" || k == "--verify" || k == "--cohort") { action = k; continue; }
            if (k == "--ramp") { ramp = true; continue; }
            if (a + 1 == argc) throw std::runtime_error("missing value for " + k);
            std::string v = argv[++a];
            if (k == "--uuid") uuid = v; else if (k == "--output") output = v; else if (k == "--steps") S = std::stoi(v);
            else if (k == "--trials") trials = std::stoi(v); else if (k == "--per-launch") per_launch = std::stoi(v);
            else if (k == "--engines") E = std::stoi(v);
            else if (k == "--t0") t0 = std::stod(v); else if (k == "--t1") t1 = std::stod(v); else if (k == "--q") q = std::stod(v);
            else if (k == "--lambda") lam = std::stod(v); else if (k == "--tec-jv") jv = std::stod(v);
            else if (k == "--seed") seed = std::stoull(v); else if (k == "--offset") offset = uint32_t(std::stoul(v));
            else if (k == "--clock-mhz") mhz = std::stod(v);
            else throw std::runtime_error("unknown argument " + k);
        }
        Device dev;
        if (action == "--probe") { std::cout << "{\"logic_uuid\":\"" << dev.uuid() << "\",\"hardware_writes\":false}\n"; return 0; }
        const char* task = std::getenv("TASK_ROOT");
        if (!task) throw std::runtime_error("TASK_ROOT missing");
        if (uuid.empty() || dev.uuid() != uuid) throw std::runtime_error("expected --uuid does not match board");
        if (output.empty() || !(mhz > 0 && mhz < 1000) || S < 2 || S > 4096 || trials < 1 || per_launch < 1 || per_launch > 512 ||
            E < 1 || E > MAX_E || trials % (E * per_launch) != 0)
            throw std::runtime_error("need --output --clock-mhz, valid schedule, 1<=engines<=5, trials divisible by engines*per-launch");
        if (std::ifstream(output + ".jsonl").good()) throw std::runtime_error("refusing to overwrite " + output);
        sca::Graph g = sca::load_graph(std::string(task) + "/data/K2000.bin");
        if (g.sumw != -1040) throw std::runtime_error("unexpected graph");
        std::vector<double> L(S, lam);
        if (ramp) for (int t = 0; t < S; ++t) if (t >= 0.7 * S) L[t] = lam * (S - 1 - t) / std::max(1.0, S - 1 - 0.7 * S);
        sca::Tables tb = sca::make_tables(geom(t0, t1, S), std::vector<double>(S, q), L, jv);
        const bool verify = action == "--verify";
        if (verify && (trials > 64 || per_launch * S * 4 > int(TRACE_STRIDE))) throw std::runtime_error("verify: small runs only");

        check(ami_dev_request_access(dev.h), "Request AMI access");
        for (int e = 0; e < E; ++e) {
            if (dev.rd(e, SCA_TIMER_ID) != TIMER_MAGIC) throw std::runtime_error("timing monitor identity missing on engine " + std::to_string(e));
            if (!(dev.rd(e, 0) & 4)) throw std::runtime_error("engine not idle");
        }
        std::vector<uint64_t> J64(size_t(sca::NP) * sca::WORDS, 0);
        for (int x = 0; x < sca::N; ++x) for (int w = 0; w < sca::WORDS; ++w) J64[size_t(x) * sca::WORDS + w] = g.J[x][w];
        dev.write(J_OFF, J64.data(), J64.size() * 8);
        std::vector<uint64_t> tab(size_t(4) * S, 0);
        for (int t = 0; t < S; ++t) {
            tab[4 * t] = uint64_t(uint32_t(tb.fourT[t])) | (uint64_t(uint32_t(tb.q[t])) << 32);
            tab[4 * t + 1] = uint64_t(int64_t(tb.kconst[t]));
            tab[4 * t + 2] = uint64_t(tb.kcorr[t]);
        }
        dev.write(TAB_OFF, tab.data(), tab.size() * 8);
        for (int e = 0; e < E; ++e) {
            dev.wr64(e, XSCA_V80_CONTROL_ADDR_J64_DATA, HBM + J_OFF);
            dev.wr64(e, XSCA_V80_CONTROL_ADDR_TAB_IN_DATA, HBM + TAB_OFF);
            dev.wr64(e, XSCA_V80_CONTROL_ADDR_OUT_R_DATA, HBM + OUT_BASE + uint64_t(e) * OUT_STRIDE);
            dev.wr64(e, XSCA_V80_CONTROL_ADDR_TRACE_DATA, HBM + TRACE_BASE + uint64_t(e) * TRACE_STRIDE);
            dev.wr(e, XSCA_V80_CONTROL_ADDR_S_DATA, uint32_t(S));
            dev.wr64(e, XSCA_V80_CONTROL_ADDR_SEED_DATA, seed);
            dev.wr64(e, XSCA_V80_CONTROL_ADDR_SUMW_DATA, uint64_t(g.sumw));
            dev.wr(e, XSCA_V80_CONTROL_ADDR_DO_TRACE_DATA, verify ? 1u : 0u);
            dev.wr(e, 0x04, 0);
        }
        std::ofstream rec(output + ".jsonl"), states(output + ".spins.bin", std::ios::binary);
        if (!rec || !states) throw std::runtime_error("cannot create result files");
        rec << std::setprecision(12);
        int done = 0, fails = 0, round = 0;
        while (done < trials) {
            const int n = per_launch;
            std::vector<uint32_t> off(E);
            for (int e = 0; e < E; ++e) {
                off[e] = offset + uint32_t(done + e * n);
                dev.wr(e, XSCA_V80_CONTROL_ADDR_TRIALS_DATA, uint32_t(n));
                dev.wr(e, XSCA_V80_CONTROL_ADDR_TRIAL_OFFSET_DATA, off[e]);
                dev.wr(e, XSCA_V80_CONTROL_ADDR_LOAD_J_DATA, round == 0 ? 1u : 0u);
                if (!(dev.rd(e, 0) & 4)) throw std::runtime_error("engine not idle before launch");
            }
            const auto t_sub = std::chrono::steady_clock::now();
            for (int e = 0; e < E; ++e) dev.wr(e, 0, 1);
            for (int e = 0; e < E; ++e)
                while (!(dev.rd(e, 0) & 2))
                    if (std::chrono::steady_clock::now() - t_sub > std::chrono::minutes(5)) throw std::runtime_error("kernel timeout");
            const auto t_done = std::chrono::steady_clock::now();
            std::vector<uint64_t> cyc(E); uint64_t cmax = 0;
            for (int e = 0; e < E; ++e) { cyc[e] = dev.cycles(e); cmax = std::max(cmax, cyc[e]); }
            for (int e = 0; e < E; ++e) {
                auto outw = dev.read(OUT_BASE + e * OUT_STRIDE, size_t(n) * 36 * 8);
                std::vector<uint64_t> out(size_t(n) * 36); std::memcpy(out.data(), outw.data(), out.size() * 8);
                std::vector<int32_t> trace;
                if (verify) { auto tw = dev.read(TRACE_BASE + e * TRACE_STRIDE, size_t(n) * S * 4); trace.resize(tw.size()); std::memcpy(trace.data(), tw.data(), tw.size() * 4); }
                dev.wr(e, 0, 0x10);
                for (int i = 0; i < n; ++i) {
                    const uint64_t* o = &out[size_t(i) * 36];
                    const int64_t cut_dev = int64_t(o[32]), flips = int64_t(o[33]), cut_host = scalar_cut(g, o);
                    bool ok = cut_dev == cut_host && (o[31] >> 16) == 0, exact = true;
                    if (verify) {
                        sca::Result r = sca::run_trial(g, tb, seed, off[e] + uint32_t(i), true);
                        for (int w = 0; w < sca::WORDS; ++w) exact &= o[w] == r.spins[w];
                        exact &= r.cut == cut_dev && r.flips == flips;
                        for (int t = 0; t < S; ++t) exact &= trace[size_t(i) * S + t] == r.n_lin[t];
                        ok &= exact;
                    }
                    fails += !ok;
                    states.write(reinterpret_cast<const char*>(o), 32 * 8);
                    rec << "{\"trial_id\":" << off[e] + i << ",\"engine\":" << e << ",\"round\":" << round << ",\"seed\":" << seed
                        << ",\"steps\":" << S << ",\"t0\":" << t0 << ",\"t1\":" << t1 << ",\"q\":" << q << ",\"lambda\":" << lam
                        << ",\"ramp\":" << (ramp ? "true" : "false") << ",\"tec_jv\":" << jv << ",\"cut\":" << cut_host
                        << ",\"device_cut\":" << cut_dev << ",\"success\":" << (cut_host >= 33000 ? "true" : "false")
                        << ",\"flips\":" << flips << ",\"launch_trials\":" << n << ",\"launch_cycles\":" << cyc[e]
                        << ",\"round_max_cycles\":" << cmax << ",\"engines\":" << E
                        << ",\"launch_device_ms\":" << cyc[e] / (mhz * 1000.0) << ",\"round_device_ms\":" << cmax / (mhz * 1000.0)
                        << ",\"round_host_ms\":" << std::chrono::duration<double, std::milli>(t_done - t_sub).count()
                        << ",\"load_J\":" << (round == 0 ? "true" : "false")
                        << ",\"scores_agree\":" << (cut_dev == cut_host ? "true" : "false")
                        << ",\"exact_reference\":" << (verify ? (exact ? "\"PASS\"" : "\"FAIL\"") : "null") << "}\n";
                }
            }
            rec.flush(); states.flush();
            done += E * n; ++round;
            std::cerr << "completed " << done << "/" << trials << ", round_max_cycles=" << cmax << ", failures=" << fails << "\n";
        }
        std::cout << "{\"action\":\"" << action << "\",\"engines\":" << E << ",\"trials\":" << trials << ",\"failures\":" << fails << "}\n";
        return fails ? 2 : 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << "\n";
        return 1;
    }
}
