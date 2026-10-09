// Revision-4 runner (sca_core_mb_r4.v + sca_mover_mb3.cpp): copy of host_sca_multi_mb3.cpp (frozen for Amendment 3) with
// helpers from mb_common4.hpp (adds --matrix jint8:<J file>:<bias file>) and S up to 2^MB_TLOG (default 13: 8192).
// (host_sca_multi_mb3.cpp text follows) Revision-3 runner (sca_core_mb_r3.v + sca_mover_mb3.cpp): copy of host_sca_multi_mb.cpp (frozen) with
//  * runtime active spins n = the problem's size (1..2048; header H0[47:36] via load_J bits 27..16), so G-set N = 800 / 1000
//    run with the n < 2048 semantics of sca_ref_bias.hpp; output bits of spins >= n must be 0;
//  * per-spin biases when the image has them (-DMB_BIAS=1): bias mode (load_J bit 2) for problems with a linear term
//    (--matrix ising:<file> or rand:...:<n>:<bmax>); the 128-word bias segment follows the planes in HBM; output word 32
//    is then sum_i s_i h_i and the host reports the energy H = -(sum_sh + sum_bs) / 2 (--target-energy: success if H <= it);
//  * golden model sca_ref_bias.hpp run_trial_bias (unchanged shared file) for --verify;
//  * host checks of the engine's field width: max_y sum_x |J_xy| + max |b| < 2^(HB-1), HB = MB_HBX or 1 + clog2(2047 M + 1);
//  * tables at J_BYTES + 64 KiB (the bias segment sits at J_BYTES).
// (Original header follows.)
// Physical runner for the multi-bit engine (sca_core_mb.v + sca_mover_mb.cpp, MB_K coupling planes), 1..12 engines.
// Copy of host_sca_multi_v6_power.cpp (rev 2) with integer coupling matrices (mb_common.hpp); the frozen v6 hosts are
// unchanged. Differences:
//  * --matrix k2000:<file> | gset:<file> | rand:<seed>:<density>:<maxabs> (default k2000:$TASK_ROOT/data/K2000.bin);
//    the engine gets MB_K sign-magnitude planes; the golden model is sca_ref.hpp run_trial(..., dense), unchanged.
//  * --target <cut> success threshold (default 33000, the K2000 target); every state is rescored independently
//    (K2000: packed +-1 graph as before; G-set: raw edge list).
//  * --corr-scale <x|auto>: correction scale folded into lambda (Onsager) and kappa (TEC-T); auto = mean_i sum_j J_ij^2 / N
//    (MULTIBIT_SPEC). Default 1, which reproduces the v6 tables exactly (the K2000 regression).
//  * HBM layout (BAR2, 16 MiB): J at 0 (MB_K x 512 KiB), packed tables after it, then out_e and trace_e (see below).
// --verify: seeded trials with per-step trace, each checked bit for bit against the reference.
// --cohort: timed trials; --power-run <s>: sustained back-to-back load with no readback (board power only).
#include "mb_common4.hpp"
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
#ifndef MB_K
#error "define MB_K (coupling planes of the image)"
#endif
#ifndef MB_BIAS
#define MB_BIAS 0
#endif
#ifndef MB_HBX
#define MB_HBX 0
#endif
#ifndef MB_TLOG
#define MB_TLOG 13
#endif

namespace {
constexpr uint64_t CONTROL = 0x01100000, HBM = 0x004000000000;
constexpr int MAX_E = 12;
constexpr uint32_t BAR2_SIZE = 0x01000000;
constexpr uint32_t J_OFF = 0, J_BYTES = uint32_t(MB_K) * 0x80000;
constexpr uint32_t B_OFF = J_BYTES, TAB_OFF = J_BYTES + 0x10000, OUT_BASE = TAB_OFF + 0x40000, OUT_STRIDE = 0x40000;
constexpr int HB_ENGINE = MB_HBX > 0 ? MB_HBX : [] { int m = (1 << (MB_K - 1)) - 1, v = 2047 * m, b = 0; while ((1 << b) < v + 1) ++b; return 1 + b; }();
constexpr uint32_t TRACE_BASE = OUT_BASE + MAX_E * OUT_STRIDE + 0x40000;
constexpr uint32_t TRACE_STRIDE = ((BAR2_SIZE - TRACE_BASE) / MAX_E) & ~uint32_t(0xFFFF);
static_assert(TRACE_STRIDE >= 0x80000, "trace window too small");
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
}  // namespace

int main(int argc, char** argv) {
    try {
        std::string action, uuid, output, matrix, scale_arg = "1", bias_arg = "auto";
        int64_t target_energy = INT64_MIN;
        int S = 560, trials = 1, per_launch = 1, E = 1;
        double t0 = 30, t1 = 5, q = 4, lam = 0, jv = 0, mhz = 0, kappa = 0, power_s = 0;
        int64_t target = 33000;
        bool ramp = false, keep_tables = false;
        uint64_t seed = 20261004; uint32_t offset = 0;
        for (int a = 1; a < argc; ++a) {
            std::string k = argv[a];
            if (k == "--probe" || k == "--verify" || k == "--cohort") { action = k; continue; }
            if (k == "--ramp") { ramp = true; continue; }
            if (k == "--keep-tables") { keep_tables = true; continue; }
            if (a + 1 == argc) throw std::runtime_error("missing value for " + k);
            std::string v = argv[++a];
            if (k == "--uuid") uuid = v; else if (k == "--output") output = v; else if (k == "--steps") S = std::stoi(v);
            else if (k == "--trials") trials = std::stoi(v); else if (k == "--per-launch") per_launch = std::stoi(v);
            else if (k == "--engines") E = std::stoi(v);
            else if (k == "--t0") t0 = std::stod(v); else if (k == "--t1") t1 = std::stod(v); else if (k == "--q") q = std::stod(v);
            else if (k == "--lambda") lam = std::stod(v); else if (k == "--tec-jv") jv = std::stod(v);
            else if (k == "--tecT") kappa = std::stod(v);
            else if (k == "--seed") seed = std::stoull(v); else if (k == "--offset") offset = uint32_t(std::stoul(v));
            else if (k == "--clock-mhz") mhz = std::stod(v);
            else if (k == "--power-run") power_s = std::stod(v);
            else if (k == "--matrix") matrix = v; else if (k == "--target") target = std::stoll(v);
            else if (k == "--corr-scale") scale_arg = v;
            else if (k == "--target-energy") target_energy = std::stoll(v);
            else if (k == "--bias-mode") bias_arg = v;
            else throw std::runtime_error("unknown argument " + k);
        }
        Device dev;
        if (action == "--probe") { std::cout << "{\"logic_uuid\":\"" << dev.uuid() << "\",\"hardware_writes\":false}\n"; return 0; }
        const char* task = std::getenv("TASK_ROOT");
        if (!task) throw std::runtime_error("TASK_ROOT missing");
        if (uuid.empty() || dev.uuid() != uuid) throw std::runtime_error("expected --uuid does not match board");
        if (output.empty() || !(mhz > 0 && mhz < 1000) || S < 2 || S > (1 << MB_TLOG) || trials < 1 || per_launch < 1 || per_launch > 512 ||
            E < 1 || E > MAX_E || trials % (E * per_launch) != 0)
            throw std::runtime_error("need --output --clock-mhz, valid schedule, 1<=engines<=12, trials divisible by engines*per-launch");
        if (std::ifstream(output + ".jsonl").good()) throw std::runtime_error("refusing to overwrite " + output);
        if (matrix.empty()) matrix = std::string("k2000:") + task + "/data/K2000.bin";
        const mb4::Matrix mat = mb4::load(matrix);
        if (mat.from_packed && mat.sumw != -1040) throw std::runtime_error("unexpected graph");
        const sca::Graph g = mb4::graph_of(mat);
        const double scale = scale_arg == "auto" ? mb4::j2_scale(mat) : std::stod(scale_arg);
        const bool bmode = bias_arg == "auto" ? mat.has_bias : bias_arg == "1";
        if (bmode && !MB_BIAS) throw std::runtime_error("problem has a bias but the image has no bias memory (MB_BIAS=0)");
        if (!bmode) for (int i = 0; i < sca::NP; ++i) if (mat.bias[i] != 0) throw std::runtime_error("nonzero bias without bias mode");
        if (mb4::field_bound(mat) >= (int64_t(1) << (HB_ENGINE - 1))) throw std::runtime_error("fields exceed the engine's width HB");
        if (bmode && target_energy == INT64_MIN) throw std::runtime_error("bias problems need --target-energy");
        const uint32_t nflags = (uint32_t(mat.n) << 16) | (bmode ? 4u : 0u);   // r3: n and bias mode on every launch
        mb4::Sched sc; sc.t0 = t0; sc.t1 = t1; sc.S = S; sc.q = q; sc.lam = lam; sc.ramp = ramp; sc.jv = jv; sc.kappa = kappa;
        sc.scale = scale;
        const sca::Tables tb = mb4::tables(sc);
        mb4::check_tables(tb, MB_K);
        const std::vector<uint64_t> planes = mb4::planes(mat, MB_K);
        const std::vector<uint64_t> bwords = mb4::bias_beats(mat);   // 512 x 64 bit = 128 words of 256 bit
        const bool verify = action == "--verify";
        if (verify && (trials > 64 || uint64_t(per_launch) * S * 4 > TRACE_STRIDE)) throw std::runtime_error("verify: small runs only");
        if (uint64_t(per_launch) * 36 * 8 > OUT_STRIDE) throw std::runtime_error("per-launch output exceeds its window");

        check(ami_dev_request_access(dev.h), "Request AMI access");
        for (int e = 0; e < E; ++e) {
            if (dev.rd(e, SCA_TIMER_ID) != TIMER_MAGIC) throw std::runtime_error("timing monitor identity missing on engine " + std::to_string(e));
            if (!(dev.rd(e, 0) & 4)) throw std::runtime_error("engine not idle");
        }
        dev.write(J_OFF, planes.data(), planes.size() * 8);
        if (planes.size() * 8 != J_BYTES) throw std::runtime_error("plane size");
        if (bmode) dev.write(B_OFF, bwords.data(), bwords.size() * 8);
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
        if (power_s > 0) {   // back-to-back launches on all engines, results not read back (board-power measurement only)
            const auto ts = std::chrono::steady_clock::now();
            long trials_run = 0; double busy_s = 0, host_s = 0; int round = 0;
            double win_host = 0, win_dev = 0; int win_n = 0;
            auto t_win = ts;
            while (std::chrono::duration<double>(std::chrono::steady_clock::now() - ts).count() < power_s) {
                for (int e = 0; e < E; ++e) {
                    dev.wr(e, XSCA_V80_CONTROL_ADDR_TRIALS_DATA, uint32_t(per_launch));
                    dev.wr(e, XSCA_V80_CONTROL_ADDR_TRIAL_OFFSET_DATA, offset + uint32_t(trials_run + e * per_launch));
                    dev.wr(e, XSCA_V80_CONTROL_ADDR_LOAD_J_DATA, (round == 0 ? 1u : (keep_tables ? 2u : 0u)) | nflags);
                }
                const auto t_sub = std::chrono::steady_clock::now();
                for (int e = 0; e < E; ++e) dev.wr(e, 0, 1);
                for (int e = 0; e < E; ++e) while (!(dev.rd(e, 0) & 2)) {}
                const auto t_done = std::chrono::steady_clock::now();
                uint64_t cmax = 0;
                for (int e = 0; e < E; ++e) cmax = std::max(cmax, dev.cycles(e));
                if (round > 0) {
                    const double hs = std::chrono::duration<double>(t_done - t_sub).count(), ds = double(cmax) / (mhz * 1e6);
                    busy_s += ds; host_s += hs; win_host += hs; win_dev += ds; ++win_n;
                    if (std::chrono::duration<double>(t_done - t_win).count() > 2.0) {
                        std::cerr << std::setprecision(6) << "t=" << std::chrono::duration<double>(t_done - ts).count()
                                  << " rounds=" << win_n << " host_ms=" << 1e3 * win_host / win_n << " dev_ms=" << 1e3 * win_dev / win_n
                                  << " host/dev=" << win_host / win_dev << "\n";
                        win_host = win_dev = 0; win_n = 0; t_win = t_done;
                    }
                }
                for (int e = 0; e < E; ++e) dev.wr(e, 0, 0x10);
                trials_run += long(E) * per_launch; ++round;
            }
            const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - ts).count();
            std::cout << std::setprecision(9) << "{\"action\":\"--power-run\",\"engines\":" << E << ",\"per_launch\":" << per_launch
                      << ",\"rounds\":" << round << ",\"trials\":" << trials_run << ",\"seconds\":" << sec
                      << ",\"device_busy_seconds\":" << busy_s << ",\"device_busy_fraction\":" << busy_s / sec
                      << ",\"host_submit_to_done_seconds\":" << host_s << ",\"host_over_device\":" << host_s / busy_s
                      << ",\"matrix\":\"" << mat.source << "\",\"K\":" << MB_K << ",\"n\":" << mat.n
                      << ",\"bias_mode\":" << (bmode ? "true" : "false") << "}\n";
            return 0;
        }
        std::ofstream rec(output + ".jsonl"), states(output + ".spins.bin", std::ios::binary);
        if (!rec || !states) throw std::runtime_error("cannot create result files");
        rec << std::setprecision(12);
        char fp[17]; std::snprintf(fp, sizeof fp, "%016llx", (unsigned long long)mb4::fingerprint(mat));
        int done = 0, fails = 0, round = 0;
        while (done < trials) {
            const int n = per_launch;
            std::vector<uint32_t> off(E);
            for (int e = 0; e < E; ++e) {
                off[e] = offset + uint32_t(done + e * n);
                dev.wr(e, XSCA_V80_CONTROL_ADDR_TRIALS_DATA, uint32_t(n));
                dev.wr(e, XSCA_V80_CONTROL_ADDR_TRIAL_OFFSET_DATA, off[e]);
                dev.wr(e, XSCA_V80_CONTROL_ADDR_LOAD_J_DATA, (round == 0 ? 1u : (keep_tables ? 2u : 0u)) | nflags);
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
                    // r3: word 32 is the cut (Max-Cut) or sum_i s_i h_i (bias mode); both rescored from the spins
                    const int64_t w32 = int64_t(o[32]), flips = int64_t(o[33]);
                    int64_t cut_host = 0, sh_host = 0, bs_host = 0, energy = 0;
                    if (bmode) { const auto p = mb4::sh_bs_of(mat, o); sh_host = p.first; bs_host = p.second; energy = -(sh_host + bs_host) / 2; }
                    else cut_host = mb4::cut_of(mat, o);
                    const int64_t score_host = bmode ? sh_host : cut_host;
                    bool pad0 = true;   // bits of spins >= n are 0
                    for (int b = mat.n; b < sca::NP; ++b) pad0 &= ((o[b >> 6] >> (b & 63)) & 1) == 0;
                    bool ok = w32 == score_host && pad0, exact = true;
                    if (verify) {
                        const sca::ResultBias r = sca::run_trial_bias(mat.n, mat.J.data(), bmode ? mat.bias.data() : nullptr, tb, seed,
                                                                      off[e] + uint32_t(i), true);
                        for (int w = 0; w < sca::WORDS; ++w) exact &= o[w] == r.spins[w];
                        exact &= (bmode ? r.sum_sh : (g.sumw + r.sum_sh / 2) / 2) == w32 && r.flips == flips;
                        for (int t = 0; t < S; ++t) exact &= trace[size_t(i) * S + t] == r.n_lin[t];
                        ok &= exact;
                    }
                    fails += !ok;
                    states.write(reinterpret_cast<const char*>(o), 32 * 8);
                    rec << "{\"trial_id\":" << off[e] + i << ",\"engine\":" << e << ",\"round\":" << round << ",\"seed\":" << seed
                        << ",\"steps\":" << S << ",\"t0\":" << t0 << ",\"t1\":" << t1 << ",\"q\":" << q << ",\"lambda\":" << lam
                        << ",\"ramp\":" << (ramp ? "true" : "false") << ",\"tec_jv\":" << jv << ",\"tecT_kappa\":" << kappa
                        << ",\"corr_scale\":" << scale << ",\"K\":" << MB_K << ",\"n\":" << mat.n << ",\"bias_mode\":" << (bmode ? "true" : "false")
                        << ",\"matrix\":\"" << mat.source << "\",\"matrix_fnv\":\"" << fp << "\"";
                    if (bmode)
                        rec << ",\"target_energy\":" << target_energy << ",\"energy\":" << energy << ",\"sum_sh\":" << sh_host
                            << ",\"device_sum_sh\":" << w32 << ",\"success\":" << (energy <= target_energy ? "true" : "false");
                    else
                        rec << ",\"target\":" << target << ",\"cut\":" << cut_host
                            << ",\"device_cut\":" << w32 << ",\"success\":" << (cut_host >= target ? "true" : "false");
                    rec << ",\"padding_zero\":" << (pad0 ? "true" : "false")
                        << ",\"flips\":" << flips << ",\"launch_trials\":" << n << ",\"launch_cycles\":" << cyc[e]
                        << ",\"round_max_cycles\":" << cmax << ",\"engines\":" << E
                        << ",\"launch_device_ms\":" << cyc[e] / (mhz * 1000.0) << ",\"round_device_ms\":" << cmax / (mhz * 1000.0)
                        << ",\"round_host_ms\":" << std::chrono::duration<double, std::milli>(t_done - t_sub).count()
                        << ",\"load_J\":" << (round == 0 ? "true" : "false")
                        << ",\"scores_agree\":" << (w32 == score_host ? "true" : "false")
                        << ",\"exact_reference\":" << (verify ? (exact ? "\"PASS\"" : "\"FAIL\"") : "null") << "}\n";
                }
            }
            rec.flush(); states.flush();
            done += E * n; ++round;
            std::cerr << "completed " << done << "/" << trials << ", round_max_cycles=" << cmax << ", failures=" << fails << "\n";
        }
        std::cout << "{\"action\":\"" << action << "\",\"engines\":" << E << ",\"trials\":" << trials << ",\"failures\":" << fails
                  << ",\"matrix\":\"" << mat.source << "\",\"K\":" << MB_K << ",\"n\":" << mat.n << ",\"bias_mode\":" << (bmode ? "true" : "false")
                  << ",\"HB\":" << HB_ENGINE << "}\n";
        return fails ? 2 : 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << "\n";
        return 1;
    }
}
