// Integer coupling matrices for the multi-bit engine (fpga/v80_sca/MULTIBIT_SPEC.md).
// Sources: K2000.bin (SNOWGPU1 packed +-1, via sca::load_graph), G-set edge lists ("n m" then "i j w", 1-indexed; J = -w),
// and the dense int8 binary format JINT8 below. Everything is converted to the dense NP x NP int8 array that
// sca::run_trial(..., dense) consumes: dense[x*NP + y] = J_xy, zero diagonal, zero padding (rows/columns N..2047).
//
// JINT8 file: 8-byte magic "SCAJINT8", uint32 N, uint32 reserved (0), then N*N int8 row-major (J[0][0..N-1], J[1][..], ...).
// BIAS file (bias version, 8 October 2026): 8-byte magic "SCABIAS1", uint32 N, uint32 reserved (0), then N int32 b_i.
// Biases are limited to |b| <= 32767 (B_BITS = 16, MULTIBIT_SPEC "Bias (external field)"); padding biases are 0.
#pragma once
#include "sca_ref.hpp"
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace jmat {

struct IntJ {
    int N = 0;
    std::vector<int8_t> dense = std::vector<int8_t>(size_t(sca::NP) * sca::NP, 0);
    int64_t sumw = 0;   // sum_{i<j} w_ij with w = -J (cut = (sumw + sum s*h / 2) / 2)
    int maxabs = 0;     // max |J_xy|
    int64_t hmax = 0;   // max_y sum_x |J_xy| (field range bound)
    double mean_row_sq = 0;   // mean_i sum_j J_ij^2
    std::string source;
    std::vector<int32_t> row_ptr, col; std::vector<int8_t> val;   // sparse rows (all nonzeros), for host rescoring
    int8_t at(int x, int y) const { return dense[size_t(x) * sca::NP + y]; }
};

inline void finish(IntJ& m) {
    if (m.N < 1 || m.N > sca::NP) throw std::runtime_error("N out of range");
    int64_t s = 0, hm = 0; int mx = 0; double sq = 0;
    for (int x = 0; x < sca::NP; ++x) {
        int64_t row = 0;
        for (int y = 0; y < sca::NP; ++y) {
            const int v = m.at(x, y);
            if ((x >= m.N || y >= m.N || x == y) && v != 0) throw std::runtime_error("nonzero diagonal or padding coupling");
            if (v != m.at(y, x)) throw std::runtime_error("J not symmetric");
            if (x < y) s -= v;
            row += std::abs(v); mx = std::max(mx, std::abs(v)); sq += double(v) * v;
        }
        hm = std::max(hm, row);
    }
    m.sumw = s; m.maxabs = mx; m.hmax = hm; m.mean_row_sq = sq / m.N;
    m.row_ptr.assign(1, 0); m.col.clear(); m.val.clear();
    for (int x = 0; x < m.N; ++x) {
        for (int y = 0; y < m.N; ++y) if (m.at(x, y)) { m.col.push_back(y); m.val.push_back(m.at(x, y)); }
        m.row_ptr.push_back(int32_t(m.col.size()));
    }
}

inline IntJ from_k2000(const std::string& path) {
    const sca::Graph g = sca::load_graph(path);
    IntJ m; m.N = sca::N; m.dense = sca::dense_matrix(g); m.source = path; finish(m);
    if (m.sumw != g.sumw) throw std::runtime_error("sumw mismatch");
    return m;
}

inline IntJ from_gset(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open " + path);
    long n = 0, e = 0;
    if (!(f >> n >> e)) throw std::runtime_error("bad G-set header");
    IntJ m; m.N = int(n); m.source = path;
    if (n < 1 || n > sca::NP) throw std::runtime_error("G-set N out of range");
    for (long k = 0; k < e; ++k) {
        long i, j, w;
        if (!(f >> i >> j >> w)) throw std::runtime_error("short G-set file");
        if (i < 1 || j < 1 || i > n || j > n || i == j) throw std::runtime_error("bad G-set edge");
        if (w < -127 || w > 127) throw std::runtime_error("G-set weight out of int8 range");
        int8_t& a = m.dense[size_t(i - 1) * sca::NP + (j - 1)];
        int8_t& b = m.dense[size_t(j - 1) * sca::NP + (i - 1)];
        if (a != 0 || b != 0) throw std::runtime_error("duplicate G-set edge");
        a = b = int8_t(-w);   // J = -w
    }
    long extra; if (f >> extra) throw std::runtime_error("trailing data in G-set file");
    finish(m);
    return m;
}

inline IntJ from_jint8(const std::string& path) {
    FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) throw std::runtime_error("cannot open " + path);
    char magic[8]; uint32_t n = 0, rsv = 0;
    if (std::fread(magic, 1, 8, f) != 8 || std::string(magic, 8) != "SCAJINT8" || std::fread(&n, 4, 1, f) != 1 ||
        std::fread(&rsv, 4, 1, f) != 1 || n < 1 || n > uint32_t(sca::NP) || rsv != 0)
        throw std::runtime_error("bad JINT8 header");
    IntJ m; m.N = int(n); m.source = path;
    std::vector<int8_t> row(n);
    for (uint32_t x = 0; x < n; ++x) {
        if (std::fread(row.data(), 1, n, f) != n) throw std::runtime_error("short JINT8 file");
        for (uint32_t y = 0; y < n; ++y) m.dense[size_t(x) * sca::NP + y] = row[y];
    }
    char c; if (std::fread(&c, 1, 1, f) == 1) throw std::runtime_error("trailing data in JINT8 file");
    std::fclose(f);
    finish(m);
    return m;
}

inline void write_jint8(const IntJ& m, const std::string& path) {
    FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) throw std::runtime_error("cannot create " + path);
    const uint32_t n = uint32_t(m.N), rsv = 0;
    std::fwrite("SCAJINT8", 1, 8, f); std::fwrite(&n, 4, 1, f); std::fwrite(&rsv, 4, 1, f);
    for (int x = 0; x < m.N; ++x) std::fwrite(&m.dense[size_t(x) * sca::NP], 1, size_t(m.N), f);
    std::fclose(f);
}

// Load by explicit kind: --graph (K2000.bin), --gset (edge list), --jint8 (dense).
inline IntJ load(const std::string& kind, const std::string& path) {
    if (kind == "graph") return from_k2000(path);
    if (kind == "gset") return from_gset(path);
    if (kind == "jint8") return from_jint8(path);
    throw std::runtime_error("unknown J source kind " + kind);
}

// Reference graph object for sca::run_trial(..., dense): only sumw is read when dense != nullptr.
inline sca::Graph ref_graph(const IntJ& m) {
    if (m.N != sca::N) throw std::runtime_error("the golden reference (sca_ref.hpp) is compiled for N = 2000");
    sca::Graph g; g.sumw = m.sumw; return g;
}

// Independent scorer: cut = sum_{x<y, s_x != s_y} w_xy, w = -J (sparse rows: every nonzero coupling, each pair once).
inline bool spin_up(const uint64_t* bits, int i) { return (bits[i >> 6] >> (i & 63)) & 1; }
inline int64_t scalar_cut(const IntJ& m, const uint64_t* bits) {
    int64_t c = 0;
    for (int x = 0; x < m.N; ++x) {
        const bool sx = spin_up(bits, x);
        for (int k = m.row_ptr[x]; k < m.row_ptr[x + 1]; ++k)
            if (m.col[k] > x && sx != spin_up(bits, m.col[k])) c -= m.val[k];
    }
    return c;
}

// Independent energy: E = -sum_{i<j} J_ij s_i s_j - sum_i b_i s_i over the active spins.
inline int64_t energy(const IntJ& m, const int32_t* bias, const uint64_t* bits) {
    int64_t e = 0;
    for (int x = 0; x < m.N; ++x) {
        const int sx = spin_up(bits, x) ? 1 : -1;
        e -= int64_t(bias[x]) * sx;
        for (int k = m.row_ptr[x]; k < m.row_ptr[x + 1]; ++k)
            if (m.col[k] > x) e -= int64_t(m.val[k]) * sx * (spin_up(bits, m.col[k]) ? 1 : -1);
    }
    return e;
}

// Number of active spins whose returned field differs from a from-scratch recomputation h_i = b_i + sum_j J_ij s_j.
inline int64_t field_mismatches(const IntJ& m, const int32_t* bias, const uint64_t* bits, const int32_t* hdev) {
    int64_t bad = 0;
    for (int x = 0; x < m.N; ++x) {
        int64_t h = bias[x];
        for (int k = m.row_ptr[x]; k < m.row_ptr[x + 1]; ++k) h += int64_t(m.val[k]) * (spin_up(bits, m.col[k]) ? 1 : -1);
        bad += h != hdev[x];
    }
    return bad;
}

// Bias files (SCABIAS1). Returns 2048 entries (zero padding).
inline std::vector<int32_t> load_bias(const std::string& path, int N) {
    FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) throw std::runtime_error("cannot open " + path);
    char magic[8]; uint32_t n = 0, rsv = 0;
    if (std::fread(magic, 1, 8, f) != 8 || std::string(magic, 8) != "SCABIAS1" || std::fread(&n, 4, 1, f) != 1 ||
        std::fread(&rsv, 4, 1, f) != 1 || rsv != 0)
        throw std::runtime_error("bad bias header");
    if (int(n) != N) throw std::runtime_error("bias length differs from N");
    std::vector<int32_t> b(sca::NP, 0);
    if (std::fread(b.data(), 4, n, f) != n) throw std::runtime_error("short bias file");
    char c; if (std::fread(&c, 1, 1, f) == 1) throw std::runtime_error("trailing data in bias file");
    std::fclose(f);
    for (uint32_t i = 0; i < n; ++i) if (b[i] < -32767 || b[i] > 32767) throw std::runtime_error("bias outside |b| <= 32767 (B_BITS = 16)");
    return b;
}
inline void write_bias(const std::vector<int32_t>& b, int N, const std::string& path) {
    FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) throw std::runtime_error("cannot create " + path);
    const uint32_t n = uint32_t(N), rsv = 0;
    std::fwrite("SCABIAS1", 1, 8, f); std::fwrite(&n, 4, 1, f); std::fwrite(&rsv, 4, 1, f); std::fwrite(b.data(), 4, n, f);
    std::fclose(f);
}

}  // namespace jmat
