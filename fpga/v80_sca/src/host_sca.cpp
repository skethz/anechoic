// Physical runner for the V80 Onsager-SCA engine.
// --verify: seeded trials with per-step trace, each checked bit for bit against sca_ref.hpp.
// --cohort: timed trials (device cycles from the passive counter), with every returned state rescored independently.
// Register offsets come from the HLS-generated driver header (xsca_v80_hw.h); the timer offsets are set at compile time
// to match prepare_ip.py (SCA_TIMER_ID, SCA_CYC_LO, SCA_CYC_HI).
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
constexpr uint32_t J_OFF = 0x000000, FOURT_OFF = 0x080000, Q_OFF = 0x084000, KCONST_OFF = 0x088000, KCORR_OFF = 0x08C000;
constexpr uint32_t OUT_OFF = 0x100000, TRACE_OFF = 0x200000, BAR2_SIZE = 0x01000000;
constexpr uint32_t TIMER_MAGIC = 0x53434131;  // "SCA1"

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
    uint32_t rd(uint32_t off) { uint32_t v = 0; check(ami_mem_bar_read(h, 0, CONTROL + off, &v), "Read control"); return v; }
    void wr(uint32_t off, uint32_t v) { check(ami_mem_bar_write(h, 0, CONTROL + off, v), "Write control"); }
    void wr64(uint32_t off, uint64_t v) { wr(off, uint32_t(v)); wr(off + 4, uint32_t(v >> 32)); }
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
    uint64_t cycles() { return uint64_t(rd(SCA_CYC_LO)) | (uint64_t(rd(SCA_CYC_HI)) << 32); }
};

std::vector<double> geom(double t0, double t1, int S) {
    std::vector<double> T(S); for (int t = 0; t < S; ++t) T[t] = t0 * std::pow(t1 / t0, double(t) / std::max(1, S - 1)); return T;
}

// Independent scorer: scalar edge loop over unpacked spins (does not use the engine's field arithmetic).
int64_t scalar_cut(const sca::Graph& g, const uint64_t* bits) {
    int64_t c = 0;
    for (int x = 0; x < sca::N; ++x) {
        const bool sx = (bits[x >> 6] >> (x & 63)) & 1;
        for (int y = x + 1; y < sca::N; ++y)
            if (sx != (((bits[y >> 6] >> (y & 63)) & 1) != 0)) c += -sca::jval(g, x, y);  // w = -J
    }
    return c;
}
}  // namespace

int main(int argc, char** argv) {
    try {
        std::string action, uuid, output;
        int S = 560, trials = 1, per_launch = 1, do_trace = 0;
        double t0 = 30, t1 = 5, q = 4, lam = 0, jv = 0, mhz = 0, kappa = 0;
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
            else if (k == "--t0") t0 = std::stod(v); else if (k == "--t1") t1 = std::stod(v); else if (k == "--q") q = std::stod(v);
            else if (k == "--lambda") lam = std::stod(v); else if (k == "--tec-jv") jv = std::stod(v);
            else if (k == "--tecT") kappa = std::stod(v);
            else if (k == "--seed") seed = std::stoull(v); else if (k == "--offset") offset = uint32_t(std::stoul(v));
            else if (k == "--clock-mhz") mhz = std::stod(v);
            else throw std::runtime_error("unknown argument " + k);
        }
        Device dev;
        if (action == "--probe") {
            std::cout << "{\"logic_uuid\":\"" << dev.uuid() << "\",\"hardware_writes\":false}\n";
            return 0;
        }
        const char* task = std::getenv("TASK_ROOT");
        if (!task) throw std::runtime_error("TASK_ROOT missing");
        if (uuid.empty() || dev.uuid() != uuid) throw std::runtime_error("expected --uuid does not match board");
        if (output.empty() || !(mhz > 0 && mhz < 1000) || S < 2 || S > 4096 || trials < 1 || per_launch < 1 || per_launch > 1024)
            throw std::runtime_error("need --output --clock-mhz and valid schedule/trials");
        if (std::ifstream(output + ".jsonl").good()) throw std::runtime_error("refusing to overwrite " + output);
        sca::Graph g = sca::load_graph(std::string(task) + "/data/K2000.bin");
        if (g.sumw != -1040) throw std::runtime_error("unexpected graph");
        std::vector<double> L(S, lam);
        if (ramp) for (int t = 0; t < S; ++t) if (t >= 0.7 * S) L[t] = lam * (S - 1 - t) / std::max(1.0, S - 1 - 0.7 * S);
        sca::Tables tb = sca::make_tables(geom(t0, t1, S), std::vector<double>(S, q), L, jv);
        if (kappa > 0) {  // TEC-T: field - kappa*T_{t-1}*ramp(t)*s(t-1), i.e. kconst[t] = kappa*T_{t-1}*ramp(t) in Q16.16, no popcount term
            const std::vector<double> Tg = geom(t0, t1, S);
            for (int t = 1; t < S; ++t) {
                const double rf = (ramp && t >= 0.7 * S) ? (S - 1 - t) / std::max(1.0, S - 1 - 0.7 * S) : 1.0;
                tb.kconst[t] = int32_t(std::llround(kappa * Tg[t - 1] * rf * 65536.0));
                tb.kcorr[t] = 0;
            }
        }
        const bool verify = action == "--verify";
        do_trace = verify ? 1 : 0;
        if (verify && (trials > 64 || per_launch * S * 4 > int(BAR2_SIZE - TRACE_OFF))) throw std::runtime_error("verify: small runs only");

        check(ami_dev_request_access(dev.h), "Request AMI access");
        if (dev.rd(SCA_TIMER_ID) != TIMER_MAGIC) throw std::runtime_error("timing monitor identity missing");
        if (!(dev.rd(0) & 4)) throw std::runtime_error("kernel not idle");
        std::vector<uint64_t> J64(size_t(sca::NP) * sca::WORDS, 0);
        for (int x = 0; x < sca::N; ++x) for (int w = 0; w < sca::WORDS; ++w) J64[size_t(x) * sca::WORDS + w] = g.J[x][w];
        dev.write(J_OFF, J64.data(), J64.size() * 8);
        dev.write(FOURT_OFF, tb.fourT.data(), S * 4); dev.write(Q_OFF, tb.q.data(), S * 4);
        dev.write(KCONST_OFF, tb.kconst.data(), S * 4); dev.write(KCORR_OFF, tb.kcorr.data(), S * 8);
        dev.wr64(XSCA_V80_CONTROL_ADDR_J64_DATA, HBM + J_OFF);
        dev.wr64(XSCA_V80_CONTROL_ADDR_FOURT_IN_DATA, HBM + FOURT_OFF);
        dev.wr64(XSCA_V80_CONTROL_ADDR_Q_IN_DATA, HBM + Q_OFF);
        dev.wr64(XSCA_V80_CONTROL_ADDR_KCONST_IN_DATA, HBM + KCONST_OFF);
        dev.wr64(XSCA_V80_CONTROL_ADDR_KCORR_IN_DATA, HBM + KCORR_OFF);
        dev.wr64(XSCA_V80_CONTROL_ADDR_OUT_R_DATA, HBM + OUT_OFF);
        dev.wr64(XSCA_V80_CONTROL_ADDR_TRACE_DATA, HBM + TRACE_OFF);
        dev.wr(XSCA_V80_CONTROL_ADDR_S_DATA, uint32_t(S));
        dev.wr64(XSCA_V80_CONTROL_ADDR_SEED_DATA, seed);
        dev.wr64(XSCA_V80_CONTROL_ADDR_SUMW_DATA, uint64_t(g.sumw));
        dev.wr(XSCA_V80_CONTROL_ADDR_DO_TRACE_DATA, uint32_t(do_trace));
        dev.wr(0x04, 0);  // global interrupt off; host polls

        std::ofstream rec(output + ".jsonl"), states(output + ".spins.bin", std::ios::binary);
        if (!rec || !states) throw std::runtime_error("cannot create result files");
        rec << std::setprecision(12);
        bool first = true; int done = 0, fails = 0;
        while (done < trials) {
            const int n = std::min(per_launch, trials - done);
            const uint32_t off = offset + uint32_t(done);
            dev.wr(XSCA_V80_CONTROL_ADDR_TRIALS_DATA, uint32_t(n));
            dev.wr(XSCA_V80_CONTROL_ADDR_TRIAL_OFFSET_DATA, off);
            dev.wr(XSCA_V80_CONTROL_ADDR_LOAD_J_DATA, first ? 1u : 0u);
            if (!(dev.rd(0) & 4)) throw std::runtime_error("kernel not idle before launch");
            const auto t_sub = std::chrono::steady_clock::now();
            dev.wr(0, 1);
            while (!(dev.rd(0) & 2))
                if (std::chrono::steady_clock::now() - t_sub > std::chrono::minutes(5)) throw std::runtime_error("kernel timeout");
            const auto t_done = std::chrono::steady_clock::now();
            const uint64_t cyc = dev.cycles();
            auto outw = dev.read(OUT_OFF, size_t(n) * 36 * 8);
            std::vector<uint64_t> out(size_t(n) * 36); std::memcpy(out.data(), outw.data(), out.size() * 8);
            std::vector<int32_t> trace;
            if (verify) { auto tw = dev.read(TRACE_OFF, size_t(n) * S * 4); trace.resize(tw.size()); std::memcpy(trace.data(), tw.data(), tw.size() * 4); }
            dev.wr(0, 0x10);  // acknowledge done
            for (int i = 0; i < n; ++i) {
                const uint64_t* o = &out[size_t(i) * 36];
                const int64_t cut_dev = int64_t(o[32]), flips = int64_t(o[33]), cut_host = scalar_cut(g, o);
                bool ok = cut_dev == cut_host && (o[31] >> 16) == 0;
                bool exact = true;
                if (verify) {
                    sca::Result r = sca::run_trial(g, tb, seed, off + uint32_t(i), true);
                    for (int w = 0; w < sca::WORDS; ++w) exact &= o[w] == r.spins[w];
                    exact &= r.cut == cut_dev && r.flips == flips;
                    for (int t = 0; t < S; ++t) exact &= trace[size_t(i) * S + t] == r.n_lin[t];
                    ok &= exact;
                }
                fails += !ok;
                states.write(reinterpret_cast<const char*>(o), 32 * 8);
                rec << "{\"trial_id\":" << off + i << ",\"seed\":" << seed << ",\"steps\":" << S << ",\"t0\":" << t0 << ",\"t1\":" << t1
                    << ",\"q\":" << q << ",\"lambda\":" << lam << ",\"ramp\":" << (ramp ? "true" : "false") << ",\"tec_jv\":" << jv << ",\"tecT_kappa\":" << kappa
                    << ",\"cut\":" << cut_host << ",\"device_cut\":" << cut_dev << ",\"success\":" << (cut_host >= 33000 ? "true" : "false")
                    << ",\"flips\":" << flips << ",\"launch_trials\":" << n << ",\"launch_cycles\":" << cyc
                    << ",\"launch_device_ms\":" << cyc / (mhz * 1000.0) << ",\"launch_host_ms\":"
                    << std::chrono::duration<double, std::milli>(t_done - t_sub).count() << ",\"load_J\":" << (first ? "true" : "false")
                    << ",\"scores_agree\":" << (cut_dev == cut_host ? "true" : "false")
                    << ",\"exact_reference\":" << (verify ? (exact ? "\"PASS\"" : "\"FAIL\"") : "null") << "}\n";
            }
            rec.flush(); states.flush();
            first = false; done += n;
            std::cerr << "completed " << done << "/" << trials << ", launch_cycles=" << cyc << ", failures=" << fails << "\n";
        }
        std::cout << "{\"action\":\"" << action << "\",\"trials\":" << trials << ",\"failures\":" << fails << "}\n";
        return fails ? 2 : 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << "\n";
        return 1;
    }
}
