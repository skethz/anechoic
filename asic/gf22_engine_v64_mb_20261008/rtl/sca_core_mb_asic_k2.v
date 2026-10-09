// sca_core_mb_asic.v -- ASIC variant (GF 22FDX estimate) of the multi-bit V80 engine core (sca_core_mb_r2.v, K-bit
// sign-magnitude couplings, MULTIBIT_SPEC.md), with the external-field bias and a runtime number of active spins.
// Generated from the unmodified sca_core_mb_r2.v by scripts/make_mb_asic_rtl.py (exact, asserted replacements).
// Bit-exact at the stream interface with sca_ref.hpp run_trial(..., dense) (b = 0, n = N) and with sca_ref_bias.hpp
// run_trial_bias(n, dense, bias, ...) (see ../sim).
//
// Changes against sca_core_mb_r2.v:
//  * Memories (as the 1.25 GHz v6.4b ASIC). Coupling residue banks at the top level: plane p of residue class e (rows
//    x = 8a + e at address a) is CPL_M = 32 / CPL_G single-port macros of 256 x (64 CPL_G) bits (CPL_G consecutive
//    groups per macro); all 32 groups read the same row in the same cycle. The macro output takes the place of the FPGA
//    read register; the group registers it once more, so the latency is unchanged. A bank is enabled only for a load
//    write or a valid flip of its extractor; the group ignores the word of a port without a valid flip (V4 = 0).
//    Step table 4096 x 128 and n_lin trace 4096 x 12: single-port macros. The table read is a two-cycle path (capture
//    with tab_ld3, next-step registers one cycle later): cycle counts unchanged.
//  * Bias (MULTIBIT_SPEC "Bias (external field)"). Bias memory 256 x 128 (2048 x int16; b_y in word y >> 3, bits
//    16 (y & 7) +: 16), loaded as its own stream segment (header bit 35) right after the coupling planes. With header
//    bit 36 the fields start at h(0) = b + sum_x J_xy s_x(0): the fields are cleared in the first cycle of the seed
//    phase, the 256 bias words are read and written into the fields during the seed phase (finished long before the
//    INIT batch; the phase would wait for it, which never happens, checked in simulation), then the INIT batch adds
//    J s(0) as before. The per-step datapath and the cycles per step are unchanged.
//  * Runtime n. Header H1[43:32] = number of active spins n (0 means N): lanes with spin index >= n are padding
//    (never updated, not counted in n_lin or flips), output spin bits >= n are zero (sca_ref_bias.hpp semantics).
//  * Widths. Fields HB = 1 + clog2(2047 M + 2^15 + 1) bits (17 at K = 2 and K = 4: any 16-bit bias with any K-bit
//    matrix), z high part ZH = HB + 1 (exact for every representable field and any int32 C).
//  * Header bit 37 (raw_score): the score word is sum_i s_i h_i (biased fields) instead of the cut.
//
// ---------------- sca_core_mb_r2.v header (unchanged) ----------------
// sca_core_mb_r2.v -- revision 2 of sca_core_mb.v (7 October 2026, after the first board run). Only change: the cut
// computation. sca_core_mb.v (built as board_mb2_e12_250mhz) inherited from v6.4b `(x + {63'd0, x[63]}) >>> 1`; the
// unsigned concatenation makes the expression unsigned, so >>> shifts logically and the halving (truncation toward zero)
// is wrong when sum_i s_i h_i < 0, i.e. when cut < sumw / 2. That never happens on K2000 (sumw = -1040, cut >= 0), so v6.4b
// is unaffected; on G22 (sumw = 19990) the all-equal period-2 states (cut 0) returned 2^62. Spins, flips and n_lin were
// always exact. Here the sign operand is $signed, so the shift is arithmetic (same module name, interface and cycles).
//
// sca_core_mb.v -- V80 Onsager-SCA engine core with K-bit couplings (MULTIBIT_SPEC.md, 7 October 2026).
// Derived from v6.4b (sca_core.v, SHA256 32e1c2f8..., kept unchanged). Default K = 2: ternary couplings {-1, 0, +1}.
// Bit-exact with ../sca_ref.hpp run_trial(..., dense) built with SCA_LANES=256, for any integer J with |J| <= 2^(K-1) - 1,
// J symmetric, J_xx = 0 (rows/columns >= N are never used). Same step pipeline, latencies and cycles per step as v6.4b.
//
// Changes against v6.4b (everything else is the v6.4b text):
//  * Coupling store. Sign-magnitude bit planes: plane 0 = sign (1 iff J < 0), planes 1..K-1 = bits 0..K-2 of |J|.
//    v6.4b keeps four 1-bit copies of every row so that its 8 read ports (one per flip extractor) see all rows. But
//    extractor e only ever issues spins x with x mod 8 == e (o_x[2:0] = E, checked in simulation), so port e only reads
//    residue class e. Here each port has its own bank holding exactly those rows, once: bank (e, p) of group b holds
//    plane p of rows x = 8a + e at address a = x >> 3 (256 x 64 bits, bit j = plane_p[x][b + 32j]). K * 2048^2 bits per
//    engine (8 Mbit at K = 2; v6.4b: 16 Mbit). Residues 0..3 are block RAM (one 256 x 64 bank per plane), residues 4..7
//    UltraRAM (one URAM per plane pair: port A plane 2u, port B plane 2u + 1, through XPM as in v6.4b). The read
//    latency is unchanged: address register (ga) -> memory read register -> output register -> field update.
//  * Field update. With tau_e = sigma_e * (-1)^sg in {+1, -1}, a flip (x_e, sigma_e) adds 2 sigma_e J = 2 m tau_e to
//    field y, so per magnitude plane b the group needs a_b = #{e: m_b & tau_e = +1} and c_b = #{e: m_b}, and
//    pend = sum_b 2^b * mode * (2 a_b - c_b) (mode = 2 in DEC, 1 in INIT): two 8-input population counts per plane,
//    no multipliers. K = 1 would reduce to v6.4b's 2 pc - nvalid.
//  * Diagonal. J_xx = 0 (spec), so the stored field is the true field (v6.4b stores J_xx = -1, i.e. h - s, and adds 1
//    back in the decision lane); the lane computes s*h directly. The decisions and the score are unchanged.
//  * Widths. Fields HB = 1 + clog2(2047 M + 1) bits (12 / 15 / 19 at K = 2 / 4 / 8), pend K + 4 bits, s*h HB + 1
//    bits, score accumulators widened to match, z high part ZH = 16 bits for K <= 4 (exact while |C| < 2^29, checked by
//    the host) and HB + 2 bits above.
//  * Load. K * 32768 beats of couplings, plane-major (plane p, row x, columns 128c..128c+127 at beat 32768 p + 16 x + c).
//
// Stream protocol (input 128-bit beats, output 64-bit words):
//   in : H0 = {seed[63:0], 29'b0, keep_tables, do_trace, load_J, S[15:0], trials[15:0]}  H1 = {sumw[63:0], 32'b0, trial_offset[31:0]}
//        if load_J: K * 32768 beats of coupling planes (see above)
//        unless keep_tables: S beats of tables {kcorr[31:0], kconst[31:0], q[31:0], fourT[31:0]}  (kcorr must fit in int32)
//   out: per trial 32 spin words (bit = 1 iff spin = +1, spins >= N zero), cut (int64), flips (uint64),
//        then, if do_trace, S words {32'b0, n_lin}.
//
// ---------------- v6.4b notes (unchanged parts) ----------------
// v6.4: 8 flip extractors; extractor e drains groups e, e+8, e+16, e+24 (32-bit round words, bit p = 4m + q'); flipped
// spin x = {k, p, e}. v6.2: header bit 34 (keep_tables) skips the table stream. v6.1: replicated enables/constants,
// 16-bit split comparisons in the lanes, registered coupling-address mux, pipelined score path, input skid buffer.
// Organisation: 32 groups. Group b owns the 64 spins/fields y = b + 32j (j = 0..63), the 8 decision lanes l = b + 32m
// (m = 0..7; lane l handles spins l + 256k in round k, i.e. j = 8k + m), and its 64 columns of the coupling banks.
`timescale 1ns/1ps
`default_nettype none

module sca_core #(
    parameter integer N = 2000,
    parameter integer K = 2,      // coupling bit planes (sign + K-1 magnitude bits); |J| <= 2^(K-1) - 1
    parameter integer DRAIN = 2,  // cycles between "extractors quiet" and the next round command (>= 2; 1 fails in simulation)
    parameter integer CPL_G = 2   // ASIC: groups per coupling macro; macro = 256 words x (64 * CPL_G) bits
) (
    input  wire         clk,
    input  wire         rst_n,
    input  wire [127:0] s_axis_tdata,
    input  wire         s_axis_tvalid,
    output wire         s_axis_tready,
    output wire [63:0]  m_axis_tdata,
    output wire         m_axis_tvalid,
    input  wire         m_axis_tready
);
    localparam integer MJ = (1 << (K - 1)) - 1;            // largest |J|
    localparam integer HB = 1 + $clog2(2047 * MJ + 32768 + 1);   // field width with a 16-bit bias (17 at K = 2, 4)
    localparam integer SW = HB + 1;                         // s*h
    localparam integer ZH = HB + 1;                         // high part of z: exact for any field and any int32 C
    localparam integer CPL_W = 64 * CPL_G;                  // ASIC: coupling macro width
    localparam integer CPL_M = 32 / CPL_G;                  // ASIC: macros per residue-bank plane
    localparam [11:0] N12 = N;
    localparam integer SCG = SW + 6;                        // per-group score accumulator (64 terms)
    localparam integer PB = (K > 2) ? $clog2(K) : 1;        // plane index bits
    localparam integer JBEATS = K * 32768;                  // coupling beats per load

    localparam [3:0] ST_HDR0 = 4'd0, ST_HDR1 = 4'd1, ST_LOADJ = 4'd2, ST_LOADT = 4'd3, ST_TRIAL = 4'd4,
                     ST_SEED = 4'd5, ST_BATCH = 4'd6, ST_WAIT = 4'd7, ST_SCOREW = 4'd8, ST_OUT = 4'd9,
                     ST_LOADB = 4'd10;
    localparam [1:0] M_INIT = 2'd0, M_DEC = 2'd1, M_SCORE = 2'd2;
    localparam integer SEED_LAT = 15;   // sd_run (aligned with sd_xa) -> sca_mix output

    reg [3:0] st;
    reg [15:0] p_trials; reg [15:0] p_S; reg p_loadJ, p_trace, p_keepT; reg [63:0] p_seed; reg [31:0] p_toff;
    reg p_loadB, p_useB, p_raw; reg [11:0] p_n;      // bias segment, biased init, raw score; active spins
    reg signed [63:0] p_sumw;
    reg [15:0] tr;           // trial index within the launch
    reg [15:0] t_step;       // number of DEC batches started in this trial (= current step + 1)
    reg [15:0] cnt;          // load counter / generic counter
    reg [19:0] jcnt;         // MB: coupling-load beat counter

    // ------------------------------------------------------------------ v6.1: input skid buffer (registered tready)
    reg [127:0] in_d, sk_d; reg in_v, sk_v;
    wire core_ready = (st == ST_HDR0) | (st == ST_HDR1) | (st == ST_LOADJ) | (st == ST_LOADT) | (st == ST_LOADB);
    wire in_fire = in_v & core_ready;
    wire in_acc = s_axis_tvalid & ~sk_v;
    assign s_axis_tready = ~sk_v;
    always @(posedge clk) begin
        if (~rst_n) begin in_v <= 1'b0; sk_v <= 1'b0; end
        else if (~in_v | in_fire) begin
            if (sk_v) begin in_d <= sk_d; in_v <= 1'b1; sk_v <= 1'b0; end
            else begin in_v <= in_acc; if (in_acc) in_d <= s_axis_tdata; end
        end else if (in_acc) begin sk_v <= 1'b1; sk_d <= s_axis_tdata; end
    end

    // ------------------------------------------------------------------ coupling load (row assembly, one plane per row)
    // rowbuf shifts with every input beat (its content only matters in LOADJ, where every valid beat is consumed);
    // the write strobe is registered, and the groups take the completed plane row from rowbuf one cycle later.
    reg [2047:0] rowbuf;
    always @(posedge clk) if (in_v) rowbuf <= {in_d, rowbuf[2047:128]};
    reg j_we_r; reg [10:0] j_addr_r; reg [PB-1:0] j_pl_r;
    always @(posedge clk) begin
        j_we_r <= (st == ST_LOADJ) & in_v & (jcnt[3:0] == 4'd15);
        j_addr_r <= jcnt[14:4];
        j_pl_r <= jcnt[15 +: PB];
    end

    // ------------------------------------------------------------------ tables and trace (ASIC: single-port SRAM macros)
    reg tab_ld3, tab_ld2, tab_ld1, tab_ld4, tab_ld5;  // table read in flight: address with tab_ld1, data in tab_q with tab_ld4
    reg trace_rd; reg trace_rd_d;
    reg [11:0] tab_raddr; reg [127:0] tab_q;
    wire [127:0] tab_q1;                       // macro output = the FPGA read register tab_q1
    wire tab_we = (st == ST_LOADT) & in_v;
    mem_tab u_tab (.clk(clk), .ce(tab_we | tab_ld1), .we(tab_we), .addr(tab_we ? cnt[11:0] : tab_raddr),
                   .wdata(in_d), .rdata(tab_q1));
    // 1.25 GHz: capture the macro output two cycles after the read edge (multicycle path 2: the output holds until the
    // next read, which is >= 8 cycles later)
    always @(posedge clk) if (tab_ld3) tab_q <= tab_q1;
    reg trace_we; reg [11:0] trace_wa; reg [11:0] trace_wd;
    reg [11:0] trace_ra; wire [11:0] trace_q;  // macro output = the FPGA read register trace_q
    mem_trc u_trc (.clk(clk), .ce(trace_we | trace_rd), .we(trace_we), .addr(trace_we ? trace_wa : trace_ra),
                   .wdata(trace_wd), .rdata(trace_q));

    // ------------------------------------------------------------------ bias memory (256 x 128 = 2048 x int16)
    // Written in LOADB (word cnt = biases 8 cnt .. 8 cnt + 7); read during the seed phase (bl_run, words 0..255 in order).
    // bl_v2 / bl_a2 / bl_d2: registered copy of the read word, broadcast to the groups.
    reg bl_run; reg [7:0] bl_a; reg bl_v1, bl_v2; reg [7:0] bl_a1, bl_a2; reg [127:0] bl_d2;
    wire bias_we = (st == ST_LOADB) & in_v;
    wire [127:0] bias_q;
    mem_bias u_bias (.clk(clk), .ce(bias_we | bl_run), .we(bias_we), .addr(bias_we ? cnt[7:0] : bl_a),
                     .wdata(in_d), .rdata(bias_q));
    always @(posedge clk) begin
        bl_v1 <= rst_n & bl_run; bl_a1 <= bl_a;
        bl_v2 <= rst_n & bl_v1; bl_a2 <= bl_a1; bl_d2 <= bias_q;
    end
    wire bl_busy = bl_run | bl_v1 | bl_v2;

    // ------------------------------------------------------------------ seeder: splitmix64 x2 per lane, chain order
    reg sd_run; reg [8:0] sd_p;
    wire [7:0] sd_lane = {sd_p[2:0], sd_p[7:3]};          // chain position p = 8b + m  ->  lane l = b + 32m
    wire [31:0] trial_id = p_toff + {16'd0, tr};
    wire [63:0] sd_x0 = p_seed ^ ({32'd0, trial_id} << 20) ^ ({56'd0, sd_lane} << 48) ^ 64'h0000_0000_5CA0_F00D;
    reg [63:0] sd_xa, sd_xb;
    always @(posedge clk) begin
        sd_xa <= sd_x0 + 64'h9E37_79B9_7F4A_7C15;
        sd_xb <= sd_x0 + 64'h3C6E_F372_FE94_F82A;      // 2 * golden (mod 2^64)
    end
    wire [63:0] sd_a, sd_b;
    sca_mix_mb mixa (.clk(clk), .z(sd_xa), .o(sd_a));
    sca_mix_mb mixb (.clk(clk), .z(sd_xb), .o(sd_b));
    reg [SEED_LAT-1:0] sd_vpipe;
    always @(posedge clk) sd_vpipe <= rst_n ? {sd_vpipe[SEED_LAT-2:0], sd_run} : {SEED_LAT{1'b0}};
    wire sd_ov = sd_vpipe[SEED_LAT-1];
    wire [127:0] sd_state_raw = {sd_b, sd_a};
    wire [127:0] sd_state = (sd_state_raw == 128'd0) ? 128'd1 : sd_state_raw;
    // v6.4b: the shift strobe gets one more register (its replication tree was route-limited); the data is delayed to match
    reg sd_ov_r; reg seed_shift_c; reg [127:0] seed_d1, seed_d2, seed_d3;
    always @(posedge clk) begin
        sd_ov_r <= sd_ov;
        seed_shift_c <= sd_ov_r;
        seed_d1 <= sd_state;
        seed_d2 <= seed_d1;
        seed_d3 <= seed_d2;
    end

    // ------------------------------------------------------------------ broadcast control registers (v6.1: replicated)
    reg lc_v; reg [2:0] lc_k; reg [1:0] lc_mode;
    reg [26:0] k_fourT; reg [31:0] k_csame, k_cdiff;
    reg [15:0] k_Uh, k_Lh; reg [3:0] k_flags;
    reg h_clr_c, sp_set_c, out_load_c, out_shift_c, amode;

    // ------------------------------------------------------------------ groups
    wire [127:0] chain [0:32];
    assign chain[32] = seed_d3;
    wire [31:0] g_exv; wire [95:0] g_exk; wire [255:0] g_exb, g_exs;
    wire [223:0] g_nl; wire [32*SCG-1:0] g_sc; wire [63:0] g_ob;
    wire [87:0] ax_bus; wire [7:0] av_bus, as_bus;
    reg [11:0] n_c;                            // active spins, broadcast to the groups
    always @(posedge clk) n_c <= p_n;

    // ------------------------------------------------------------------ ASIC: coupling residue banks (see header)
    // Write data: group b's 64 columns of the completed plane row (bit j = plane_p[x][b + 32j]), registered like the
    // FPGA group's r_jdata; address / enables registered per macro like the FPGA group's ga / wsel.
    wire [2047:0] cw_next;
    genvar gb, gj, ge, gp, gi;
    generate for (gb = 0; gb < 32; gb = gb + 1) begin : cwm
        for (gj = 0; gj < 64; gj = gj + 1) begin : jcol
            assign cw_next[64 * gb + gj] = rowbuf[gb + 32 * gj];
        end
    end endgenerate
    reg [2047:0] cw_d;
    always @(posedge clk) cw_d <= cw_next;
    wire [2048 * 8 * K - 1:0] cq;              // bank (e, p), group b: cq[2048 * (K e + p) + 64 b +: 64]
    generate for (ge = 0; ge < 8; ge = ge + 1) begin : cbank
        for (gp = 0; gp < K; gp = gp + 1) begin : pl
            for (gi = 0; gi < CPL_M; gi = gi + 1) begin : cm
                reg [7:0] ca; reg cce, cwe;
                always @(posedge clk) begin
                    ca  <= j_we_r ? j_addr_r[10:3] : ax_bus[11 * ge + 3 +: 8];
                    cwe <= j_we_r & (j_addr_r[2:0] == ge) & (j_pl_r == gp);
                    cce <= j_we_r ? ((j_addr_r[2:0] == ge) & (j_pl_r == gp)) : av_bus[ge];
                end
                mem_cpl #(.W(CPL_W)) u (.clk(clk), .ce(cce), .we(cwe), .addr(ca),
                    .wdata(cw_d[CPL_W * gi +: CPL_W]), .rdata(cq[2048 * (K * ge + gp) + CPL_W * gi +: CPL_W]));
            end
        end
    end endgenerate

    generate for (gb = 0; gb < 32; gb = gb + 1) begin : grp
        wire [8 * 64 * K - 1:0] cqg;           // this group's slices, cqg[64 * (K e + p) +: 64]
        for (ge = 0; ge < 8; ge = ge + 1) begin : cge
            for (gp = 0; gp < K; gp = gp + 1) begin : cgp
                assign cqg[64 * (K * ge + gp) +: 64] = cq[2048 * (K * ge + gp) + 64 * gb +: 64];
            end
        end
        sca_group_mb #(.B(gb), .N(N), .K(K), .HB(HB), .SW(SW), .ZH(ZH), .SCG(SCG), .PB(PB)) u (
            .clk(clk), .lc_v(lc_v), .lc_k(lc_k), .lc_mode(lc_mode),
            .k_fourT(k_fourT), .k_csame(k_csame), .k_cdiff(k_cdiff), .k_Uh(k_Uh), .k_Lh(k_Lh), .k_flags(k_flags),
            .seed_shift(seed_shift_c), .seed_in(chain[gb + 1]), .seed_out(chain[gb]),
            .h_clr(h_clr_c), .sp_set(sp_set_c),
            .av(av_bus), .asg(as_bus), .amode(amode), .cqm(cqg), .nact(n_c),
            .bl_v(bl_v2), .bl_a(bl_a2), .bl_d(bl_d2),
            .ex_v(g_exv[gb]), .ex_k(g_exk[3 * gb +: 3]), .exb(g_exb[8 * gb +: 8]), .exs(g_exs[8 * gb +: 8]),
            .nl_acc(g_nl[7 * gb +: 7]), .sc_acc(g_sc[SCG * gb +: SCG]),
            .out_load(out_load_c), .out_shift(out_shift_c), .out_bits(g_ob[2 * gb +: 2]));
    end endgenerate

    // ------------------------------------------------------------------ extractors
    wire [7:0] ex_idle;
    wire ex_in_v = g_exv[0];
    wire [2:0] ex_in_k = g_exk[2:0];
    generate for (gb = 0; gb < 8; gb = gb + 1) begin : exu
        wire [31:0] wb, ws;
        genvar q, m;
        for (q = 0; q < 4; q = q + 1) begin : qq
            for (m = 0; m < 8; m = m + 1) begin : mm
                assign wb[4 * m + q] = g_exb[8 * (gb + 8 * q) + m];
                assign ws[4 * m + q] = g_exs[8 * (gb + 8 * q) + m];
            end
        end
        sca_extractor_mb #(.E(gb)) u (
            .clk(clk), .rst(~rst_n), .in_v(ex_in_v), .in_k(ex_in_k), .in_b(wb), .in_s(ws),
            .o_v(av_bus[gb]), .o_x(ax_bus[11 * gb +: 11]), .o_s(as_bus[gb]), .idle(ex_idle[gb]));
    end endgenerate

    // words received in the current batch (counted one cycle after the slot write, so "pending" already includes it)
    reg [3:0] words_seen; reg ex_in_v_d;
    reg batch_start;   // pulse: clears words_seen
    always @(posedge clk) begin
        ex_in_v_d <= ex_in_v;
        if (batch_start) words_seen <= 4'd0;
        else if (ex_in_v_d) words_seen <= words_seen + 4'd1;
    end
    wire quiet = (words_seen == 4'd8) & (&ex_idle) & ~ex_in_v & ~ex_in_v_d;

    // flips (DEC mode only)
    reg [31:0] flips_total;
    wire [3:0] nv_issue = av_bus[0] + av_bus[1] + av_bus[2] + av_bus[3] + av_bus[4] + av_bus[5] + av_bus[6] + av_bus[7];

    // ------------------------------------------------------------------ capture: n_lin (tap 7) and score (tap 9)
    // The round-7 command is registered at the end of cycle c7; group nl_acc is final from c7+8 and the pipelined group
    // score sc_acc from c7+10; cap_sr[k] is visible at c7+1+k.
    reg [9:0] cap_sr; reg [19:0] capm_sr;
    reg cap1, cap2; reg [1:0] cap1_m, cap2_m;
    reg [9:0] nl1 [0:3];
    reg [11:0] nl2;
    integer ii;
    always @(posedge clk) begin
        cap1 <= cap_sr[7] & (capm_sr[15:14] != M_SCORE); cap1_m <= capm_sr[15:14];
        for (ii = 0; ii < 4; ii = ii + 1)
            nl1[ii] <= g_nl[7 * (8 * ii + 0) +: 7] + g_nl[7 * (8 * ii + 1) +: 7] + g_nl[7 * (8 * ii + 2) +: 7] +
                       g_nl[7 * (8 * ii + 3) +: 7] + g_nl[7 * (8 * ii + 4) +: 7] + g_nl[7 * (8 * ii + 5) +: 7] +
                       g_nl[7 * (8 * ii + 6) +: 7] + g_nl[7 * (8 * ii + 7) +: 7];
        cap2 <= cap1; cap2_m <= cap1_m;
        nl2 <= nl1[0] + nl1[1] + nl1[2] + nl1[3];
    end
    // v6.1: score tree, one adder level per register stage (32 -> 16 -> 8 -> 4 -> 1); MB: widths follow SCG
    reg sA, sB, sC, sD; reg signed [SCG:0] scA [0:15]; reg signed [SCG+1:0] scB [0:7]; reg signed [SCG+2:0] scC [0:3];
    reg signed [SCG+4:0] sc_tot;
    always @(posedge clk) begin
        sA <= cap_sr[9] & (capm_sr[19:18] == M_SCORE);
        for (ii = 0; ii < 16; ii = ii + 1)
            scA[ii] <= $signed(g_sc[SCG * (2 * ii) +: SCG]) + $signed(g_sc[SCG * (2 * ii + 1) +: SCG]);
        sB <= sA; for (ii = 0; ii < 8; ii = ii + 1) scB[ii] <= scA[2 * ii] + scA[2 * ii + 1];
        sC <= sB; for (ii = 0; ii < 4; ii = ii + 1) scC[ii] <= scB[2 * ii] + scB[2 * ii + 1];
        sD <= sC; sc_tot <= scC[0] + scC[1] + scC[2] + scC[3];
    end

    // next-step table entry and the correction pipeline: corr = int32(((nlin * kcorr) >> 8) + kconst)
    reg [31:0] tn_fourT, tn_q, tn_kconst, tn_kcorr;
    reg [31:0] tn_U, tn_L;                     // v6.1: U = twoT - 1, L = 1 - twoT (int32), from the prefetched entry
    wire [31:0] tn_twoT = {1'b0, tn_fourT[31:1]};
    always @(posedge clk) begin
        tab_ld2 <= tab_ld1; tab_ld3 <= tab_ld2; tab_ld4 <= tab_ld3; tab_ld5 <= tab_ld4;
        if (tab_ld4) begin   // 1.25 GHz: one cycle later than the FPGA core (still >= 13 cycles before first use)
            tn_fourT <= tab_q[31:0]; tn_q <= tab_q[63:32]; tn_kconst <= tab_q[95:64]; tn_kcorr <= tab_q[127:96];
        end
        if (tab_ld5) begin tn_U <= tn_twoT - 32'd1; tn_L <= 32'd1 - tn_twoT; end
    end
    reg c_go1, c_go2, c_go3, c_go4;
    reg [11:0] c_nl; reg signed [31:0] c_kc;
    reg signed [44:0] c_m, c_p;
    reg [31:0] c_corr;
    reg nk_ready; reg [26:0] nk_fourT; reg [31:0] nk_csame, nk_cdiff; reg [15:0] nk_Uh, nk_Lh; reg [3:0] nk_flags;
    always @(posedge clk) begin
        c_go1 <= cap2;
        c_nl <= (cap2_m == M_INIT) ? 12'd0 : nl2;
        c_kc <= $signed(tn_kcorr);
        c_go2 <= c_go1; c_m <= $signed({1'b0, c_nl}) * c_kc;
        c_go3 <= c_go2; c_p <= c_m;
        c_go4 <= c_go3; c_corr <= c_p[39:8] + tn_kconst;        // (c_p >>> 8)[31:0] + kconst, int32 wrap
    end
    // low halves of Csame = q - corr and Cdiff = q + corr, compared with the low halves of U and L (unsigned)
    wire [15:0] cs_lo = tn_q[15:0] - c_corr[15:0], cd_lo = tn_q[15:0] + c_corr[15:0];

    // score -> cut = (sumw + shsum/2)/2 with C truncation toward zero
    reg cut_go1, cut_go2; reg signed [63:0] cut_half, cut_val;
    wire signed [63:0] sh64 = sc_tot;                    // sign-extended (signed source)
    wire signed [63:0] cut_sum = p_sumw + cut_half;
    always @(posedge clk) begin
        cut_go1 <= sD;
        cut_half <= (sh64 + $signed({63'd0, sh64[63]})) >>> 1;     // r2: signed operand, so >>> is arithmetic
        cut_go2 <= cut_go1;
        cut_val <= (cut_sum + $signed({63'd0, cut_sum[63]})) >>> 1;
    end

    // ------------------------------------------------------------------ output FIFO
    reg [63:0] of_mem [0:7]; reg [3:0] of_wp, of_rp; wire [3:0] of_cnt = of_wp - of_rp;
    reg of_push; reg [63:0] of_din;
    wire of_pop = m_axis_tvalid & m_axis_tready;
    assign m_axis_tvalid = (of_cnt != 4'd0);
    assign m_axis_tdata = of_mem[of_rp[2:0]];
    always @(posedge clk) begin
        if (~rst_n) begin of_wp <= 4'd0; of_rp <= 4'd0; end
        else begin
            if (of_push) begin of_mem[of_wp[2:0]] <= of_din; of_wp <= of_wp + 4'd1; end
            if (of_pop) of_rp <= of_rp + 4'd1;
        end
    end
    wire of_room = (of_cnt <= 4'd4);   // at most two pushes in flight

    // spin word assembly: bit p < 32 from group p (bit 0), bit p >= 32 from group p-32 (bit 1)
    wire [63:0] spin_word;
    generate for (gb = 0; gb < 32; gb = gb + 1) begin : sw
        assign spin_word[gb] = g_ob[2 * gb];
        assign spin_word[32 + gb] = g_ob[2 * gb + 1];
    end endgenerate

    // ------------------------------------------------------------------ main FSM
    reg [3:0] iss;            // rounds issued in the current batch
    reg [1:0] b_mode;         // mode of the current batch
    reg [3:0] dr;             // drain counter
    reg last_dec;             // v6.1: the batch in flight is the last before SCORE (precomputed)
    reg [7:0] oc;             // output sequencer
    reg [1:0] ocad;           // 3-cycle cadence for spin words
    reg [15:0] ot;            // trace word index
    // output spin-word mask: om_rem = n - 64 w for the next word w (set while waiting for out_load), om_mask = bits < om_rem
    reg signed [12:0] om_rem; reg [63:0] om_mask;
    always @(posedge clk)
        om_mask <= (om_rem >= 13'sd64) ? {64{1'b1}} : (om_rem <= 13'sd0) ? 64'd0 : ((64'd1 << om_rem[5:0]) - 64'd1);

    always @(posedge clk) begin
        // defaults
        lc_v <= 1'b0; batch_start <= 1'b0; tab_ld1 <= 1'b0; trace_we <= 1'b0;
        out_load_c <= 1'b0; out_shift_c <= 1'b0; of_push <= 1'b0; trace_rd <= 1'b0; sd_run <= 1'b0;
        cap_sr <= {cap_sr[8:0], 1'b0}; capm_sr <= {capm_sr[17:0], 2'b00};
        trace_rd_d <= trace_rd;
        if (c_go4) begin
            nk_csame <= tn_q - c_corr; nk_cdiff <= tn_q + c_corr;
            nk_fourT <= tn_fourT[26:0];
            nk_Uh <= tn_U[31:16]; nk_Lh <= tn_L[31:16];
            // {ge_diff, ge_same, le_diff, le_same}: low-half comparisons used when the high halves are equal
            nk_flags <= {cd_lo >= tn_L[15:0], cs_lo >= tn_L[15:0], cd_lo <= tn_U[15:0], cs_lo <= tn_U[15:0]};
            nk_ready <= 1'b1;
        end
        if (cap2 & (cap2_m == M_DEC)) begin trace_we <= 1'b1; trace_wa <= t_step[11:0] - 12'd1; trace_wd <= nl2; end
        if (amode & (st != ST_SEED)) flips_total <= flips_total + {28'd0, nv_issue};
        if (bl_run) begin bl_a <= bl_a + 8'd1; if (bl_a == 8'd255) bl_run <= 1'b0; end

        if (~rst_n) begin
            st <= ST_HDR0; cnt <= 16'd0; jcnt <= 20'd0; nk_ready <= 1'b0; cap_sr <= 10'd0; h_clr_c <= 1'b0; sp_set_c <= 1'b0;
            amode <= 1'b0; bl_run <= 1'b0;
        end else case (st)
        ST_HDR0: if (in_fire) begin
            p_trials <= in_d[15:0]; p_S <= in_d[31:16];
            p_loadJ <= in_d[32]; p_trace <= in_d[33]; p_keepT <= in_d[34]; p_seed <= in_d[127:64];
            p_loadB <= in_d[35]; p_useB <= in_d[36]; p_raw <= in_d[37];
            st <= ST_HDR1;
        end
        ST_HDR1: if (in_fire) begin
            p_toff <= in_d[31:0]; p_sumw <= in_d[127:64];
            p_n <= (in_d[43:32] == 12'd0) ? N12 : in_d[43:32];
            cnt <= 16'd0; jcnt <= 20'd0;
            st <= p_loadJ ? ST_LOADJ : p_loadB ? ST_LOADB : ((p_S != 16'd0 && !p_keepT) ? ST_LOADT : ST_TRIAL);
            tr <= 16'd0;
        end
        ST_LOADJ: if (in_fire) begin
            jcnt <= jcnt + 20'd1;
            if (jcnt == JBEATS - 1) begin
                jcnt <= 20'd0; st <= p_loadB ? ST_LOADB : ((p_S != 16'd0 && !p_keepT) ? ST_LOADT : ST_TRIAL);
            end
        end
        ST_LOADB: if (in_fire) begin
            cnt <= cnt + 16'd1;
            if (cnt == 16'd255) begin cnt <= 16'd0; st <= (p_S != 16'd0 && !p_keepT) ? ST_LOADT : ST_TRIAL; end
        end
        ST_LOADT: if (in_fire) begin
            cnt <= cnt + 16'd1;
            if (cnt == p_S - 16'd1) begin cnt <= 16'd0; st <= ST_TRIAL; end
        end
        ST_TRIAL: begin
            if (tr == p_trials) st <= ST_HDR0;
            else begin
                st <= ST_SEED; cnt <= 16'd0; sd_p <= 9'd0;
                h_clr_c <= 1'b1; sp_set_c <= 1'b1; flips_total <= 32'd0; amode <= 1'b0;
                bl_run <= p_useB; bl_a <= 8'd0;            // bias words are written into the fields during SEED
                tab_raddr <= 12'd0; tab_ld1 <= 1'b1;      // step 0 entry
                nk_ready <= 1'b0;
            end
        end
        ST_SEED: begin
            h_clr_c <= 1'b0;   // ASIC: one-cycle field clear (nothing updates the fields during SEED except the bias load)
            // issue 256 seeder inputs, then wait until the last state is in the chain (2 cycles after sd_ov)
            if (sd_p != 9'd256) begin sd_run <= 1'b1; sd_p <= sd_p + 9'd1; end
            if (seed_shift_c) cnt <= cnt + 16'd1;
            if (cnt == 16'd256 && ~seed_shift_c && ~bl_busy) begin
                // last shift happens in the groups at the end of the next cycle; the first INIT round reads the
                // RNG two cycles after its command, so starting the batch now is safe.
                h_clr_c <= 1'b0; sp_set_c <= 1'b0;
                b_mode <= M_INIT; iss <= 4'd0; batch_start <= 1'b1; amode <= 1'b0;
                last_dec <= (p_S == 16'd0);
                st <= ST_BATCH;
            end
        end
        ST_BATCH: begin
            lc_v <= 1'b1; lc_k <= iss[2:0]; lc_mode <= b_mode;
            iss <= iss + 4'd1;
            if (iss == 4'd7) begin
                cap_sr[0] <= 1'b1; capm_sr[1:0] <= b_mode;
                st <= (b_mode == M_SCORE) ? ST_SCOREW : ST_WAIT;
                dr <= 4'd0;
            end
        end
        ST_WAIT: begin
            // quiet = all 8 words received and all extractors idle: the last flip was issued at most one cycle ago,
            // its field update is visible 6 cycles after issue; the next round's fields are read 2 cycles after its
            // command -> command no earlier than quiet + 3 (DRAIN = 2 states here plus the BATCH state).
            if ((dr != 4'd0 || quiet) && dr != 4'd15) dr <= dr + 4'd1;
            if (dr >= DRAIN[3:0]) begin
                if (last_dec) begin
                    b_mode <= M_SCORE; iss <= 4'd0; batch_start <= 1'b1; st <= ST_BATCH;
                end else if (nk_ready) begin
                    // start DEC step t_step (0 after INIT)
                    if (b_mode == M_INIT) t_step <= 16'd1; else t_step <= t_step + 16'd1;
                    last_dec <= (b_mode == M_INIT) ? (p_S == 16'd1) : (t_step + 16'd1 == p_S);
                    k_fourT <= nk_fourT; k_csame <= nk_csame; k_cdiff <= nk_cdiff; k_Uh <= nk_Uh; k_Lh <= nk_Lh;
                    k_flags <= nk_flags; nk_ready <= 1'b0;
                    tab_raddr <= (b_mode == M_INIT) ? 12'd1 : t_step[11:0] + 12'd1; tab_ld1 <= 1'b1;
                    b_mode <= M_DEC; amode <= 1'b1; iss <= 4'd0; batch_start <= 1'b1; st <= ST_BATCH;
                end
            end
        end
        ST_SCOREW: begin
            if (cut_go2) begin st <= ST_OUT; oc <= 8'd0; ocad <= 2'd0; out_load_c <= 1'b1; end
        end
        ST_OUT: begin
            // oc: 0..2 wait for out_load; 3..34 spin words (3-cycle cadence); 35 cut; 36 flips; 37 trace; 38 done
            if (oc < 8'd3) begin oc <= oc + 8'd1; om_rem <= $signed({1'b0, p_n}); end
            else if (oc < 8'd35) begin
                if (ocad == 2'd0) begin
                    if (of_room) begin
                        of_push <= 1'b1;
                        of_din <= spin_word & om_mask;      // spins >= n are zero
                        out_shift_c <= 1'b1; ocad <= 2'd1; oc <= oc + 8'd1; om_rem <= om_rem - 13'sd64;
                    end
                end else ocad <= (ocad == 2'd2) ? 2'd0 : ocad + 2'd1;
            end else if (oc == 8'd35) begin
                if (of_room) begin of_push <= 1'b1; of_din <= p_raw ? sh64 : cut_val; oc <= oc + 8'd1; end
            end else if (oc == 8'd36) begin
                if (of_room) begin of_push <= 1'b1; of_din <= {32'd0, flips_total}; oc <= oc + 8'd1; ot <= 16'd0; end
            end else if (oc == 8'd37) begin
                if (~p_trace || ot == p_S) oc <= 8'd38;
                else if (~trace_rd & ~trace_rd_d & of_room) begin trace_ra <= ot[11:0]; trace_rd <= 1'b1; end
                if (trace_rd_d) begin of_push <= 1'b1; of_din <= {52'd0, trace_q}; ot <= ot + 16'd1; end
            end else begin
                tr <= tr + 16'd1; st <= ST_TRIAL;
            end
        end
        default: st <= ST_HDR0;
        endcase
    end

`ifdef SCA_SIM_CHECKS
    // single-port legality (a memory never sees a read and a write in the same cycle), the residue property the banks
    // rely on (extractor e only issues x with x mod 8 == e), and that the bias load never delays the INIT batch
    integer chk_err = 0, chk_e;
    always @(posedge clk) if (rst_n) begin
        if (j_we_r & (|av_bus)) begin chk_err = chk_err + 1; $display("CHECK: coupling load write during a flip read"); end
        if (tab_we & tab_ld1) begin chk_err = chk_err + 1; $display("CHECK: table write during a table read"); end
        if (trace_we & trace_rd) begin chk_err = chk_err + 1; $display("CHECK: trace write during a trace read"); end
        if (bias_we & bl_run) begin chk_err = chk_err + 1; $display("CHECK: bias write during a bias read"); end
        for (chk_e = 0; chk_e < 8; chk_e = chk_e + 1)
            if (av_bus[chk_e] && (ax_bus[11 * chk_e +: 3] != chk_e)) begin
                chk_err = chk_err + 1; $display("CHECK: residue violation, extractor %0d", chk_e);
            end
        if (st == ST_SEED && cnt == 16'd256 && ~seed_shift_c && bl_busy) begin
            chk_err = chk_err + 1; $display("CHECK: bias load extends the seed phase");
        end
    end
`endif
endmodule

// --------------------------------------------------------------------------------------------------------------------
// splitmix64 output mix, pipelined: o = mix(z) after 15 cycles (z^z>>30, *C1, z^z>>27, *C2, z^z>>31). (v6.4b sca_mix)
module sca_mix_mb (
    input  wire        clk,
    input  wire [63:0] z,
    output wire [63:0] o
);
    localparam integer ML = 6;   // registers after each 64x64 multiply (absorbed into the DSP cascade)
    reg [63:0] a1, a2, o_r;
    reg [63:0] m1 [0:ML-1];
    reg [63:0] m2 [0:ML-1];
    integer i;
    always @(posedge clk) begin
        a1 <= z ^ (z >> 30);
        m1[0] <= a1 * 64'hBF58_476D_1CE4_E5B9;
        for (i = 1; i < ML; i = i + 1) m1[i] <= m1[i - 1];
        a2 <= m1[ML - 1] ^ (m1[ML - 1] >> 27);
        m2[0] <= a2 * 64'h94D0_49BB_1331_11EB;
        for (i = 1; i < ML; i = i + 1) m2[i] <= m2[i - 1];
        o_r <= m2[ML - 1] ^ (m2[ML - 1] >> 31);
    end
    assign o = o_r;
endmodule

// --------------------------------------------------------------------------------------------------------------------
// Flip extractor e (0..7): slot k holds the 32 flip bits of round k (bit p = 4m + q' -> spin x = {k, p, e}) and new spins.
// (v6.4b sca_extractor.) Note o_x[2:0] = E: extractor e only ever issues spins with x mod 8 == e.
module sca_extractor_mb #(
    parameter integer E = 0
) (
    input  wire        clk,
    input  wire        rst,
    input  wire        in_v,
    input  wire [2:0]  in_k,
    input  wire [31:0] in_b,
    input  wire [31:0] in_s,
    output reg o_v,
    output reg [10:0] o_x,
    output reg o_s,
    output wire        idle
);
    reg [31:0] slot_b [0:7];
    reg [31:0] slot_s [0:7];
    reg [7:0] pend;
    reg [31:0] cur, curs; reg [2:0] curk;
    wire [31:0] lb = cur & (~cur + 32'd1);
    wire [31:0] rest = cur & ~lb;
    wire last = (rest == 32'd0);
    wire [2:0] kn = pend[0] ? 3'd0 : pend[1] ? 3'd1 : pend[2] ? 3'd2 : pend[3] ? 3'd3 :
                    pend[4] ? 3'd4 : pend[5] ? 3'd5 : pend[6] ? 3'd6 : 3'd7;
    wire take = last & (pend != 8'd0);
    function [4:0] enc32(input [31:0] oh);
        enc32 = {|(oh & 32'hFFFF_0000), |(oh & 32'hFF00_FF00), |(oh & 32'hF0F0_F0F0), |(oh & 32'hCCCC_CCCC), |(oh & 32'hAAAA_AAAA)};
    endfunction
    always @(posedge clk) begin
        o_v <= (cur != 32'd0);
        o_x <= {curk, enc32(lb), E[2:0]};
        o_s <= |(curs & lb);
        if (rst) begin
            cur <= 32'd0; pend <= 8'd0;
        end else begin
            if (!last) cur <= rest;
            else if (take) begin cur <= slot_b[kn]; curs <= slot_s[kn]; curk <= kn; end
            else cur <= 32'd0;
            pend <= (pend & ~(take ? (8'd1 << kn) : 8'd0)) | ((in_v && in_b != 32'd0) ? (8'd1 << in_k) : 8'd0);
        end
        if (in_v) begin slot_b[in_k] <= in_b; slot_s[in_k] <= in_s; end
    end
    assign idle = (cur == 32'd0) & (pend == 8'd0);
endmodule

// --------------------------------------------------------------------------------------------------------------------
// Group b: 64 spins/fields y = b + 32j, 8 decision lanes (lane m -> spins j = 8k + m), coupling residue banks b.
module sca_group_mb #(
    parameter integer B = 0,
    parameter integer N = 2000,
    parameter integer K = 2,
    parameter integer HB = 12,
    parameter integer SW = 13,
    parameter integer ZH = 16,
    parameter integer SCG = 19,
    parameter integer PB = 1
) (
    input  wire          clk,
    input  wire          lc_v,
    input  wire [2:0]    lc_k,
    input  wire [1:0]    lc_mode,
    input  wire [26:0]   k_fourT,
    input  wire [31:0]   k_csame,
    input  wire [31:0]   k_cdiff,
    input  wire [15:0]   k_Uh,
    input  wire [15:0]   k_Lh,
    input  wire [3:0]    k_flags,
    input  wire          seed_shift,
    input  wire [127:0]  seed_in,
    output wire [127:0]  seed_out,
    input  wire          h_clr,
    input  wire          sp_set,
    input  wire [7:0]    av,
    input  wire [7:0]    asg,
    input  wire          amode,
    input  wire [8*64*K-1:0] cqm,   // ASIC: this group's 64-bit slices of the 8 x K residue-bank plane outputs
    input  wire [11:0]   nact,      // runtime number of active spins
    input  wire          bl_v,      // bias word valid (seed phase)
    input  wire [7:0]    bl_a,      // bias word address: biases 8 bl_a .. 8 bl_a + 7
    input  wire [127:0]  bl_d,      // bias word (int16 x 8)
    output reg           ex_v,
    output reg  [2:0]    ex_k,
    output reg  [7:0]    exb,
    output reg  [7:0]    exs,
    output reg  [6:0]    nl_acc,
    output reg  [SCG-1:0] sc_acc,
    input  wire          out_load,
    input  wire          out_shift,
    output wire [1:0]    out_bits
);
    localparam [1:0] M_INIT = 2'd0, M_DEC = 2'd1, M_SCORE = 2'd2;
    localparam integer PW = K + 4;   // pend width: |pend| <= 16 * (2^(K-1) - 1)
    integer m, j;

    // ---- one register hop for every broadcast input (v6.1: replicated toward their loads)
    reg g_v; reg [2:0] g_k; reg [1:0] g_mode;
    reg [26:0] r_fourT; reg [31:0] r_csame, r_cdiff;
    reg [15:0] r_Uh, r_Lh; reg [3:0] r_flags;
    reg rng_en;                   // = g_v & (g_mode != SCORE), registered directly
    reg [7:0] gv, gs; reg gm;
    reg r_hclr; reg r_spset;
    reg r_seedsh; reg r_outload, r_outshift;
    localparam [4:0] B5 = B;
    reg [11:0] r_n; reg r_blv; reg [5:0] r_blj; reg [15:0] r_bld;
    always @(posedge clk) begin
        // runtime n and the bias word: field j of this group is y = B + 32 j, its bias is in word 4 j + B[4:3], slice B[2:0]
        r_n <= nact;
        r_blv <= bl_v & (bl_a[1:0] == B5[4:3]); r_blj <= bl_a[7:2]; r_bld <= bl_d[16 * (B % 8) +: 16];
    end
    always @(posedge clk) begin
        g_v <= lc_v; g_k <= lc_k; g_mode <= lc_mode;
        rng_en <= lc_v & (lc_mode != M_SCORE);
        r_fourT <= k_fourT; r_csame <= k_csame; r_cdiff <= k_cdiff; r_Uh <= k_Uh; r_Lh <= k_Lh; r_flags <= k_flags;
        gv <= av; gs <= asg; gm <= amode;
        r_hclr <= h_clr; r_spset <= sp_set; r_seedsh <= seed_shift; r_outload <= out_load; r_outshift <= out_shift;
    end

    // ---- state
    reg [HB-1:0] h [0:63];      // field (J_xx = 0, so the true field: b_y + sum_x J_xy s_x)
    reg [PW-1:0] pend [0:63];
    reg [63:0] s, sp;
    reg [127:0] rs [0:7];        // xoshiro128** states {s3, s2, s1, s0}
    assign seed_out = rs[0];
    assign out_bits = sp[1:0];

    function [31:0] rng_out(input [127:0] x);
        reg [31:0] a, r5;
        begin
            r5 = x[63:32] * 32'd5;
            a = {r5[24:0], r5[31:25]};
            rng_out = a * 32'd9;
        end
    endfunction
    function [127:0] rng_next(input [127:0] x);
        reg [31:0] s0, s1, s2, s3, t;
        begin
            s0 = x[31:0]; s1 = x[63:32]; s2 = x[95:64]; s3 = x[127:96];
            t = s1 << 9;
            s2 = s2 ^ s0; s3 = s3 ^ s1; s1 = s1 ^ s2; s0 = s0 ^ s3; s2 = s2 ^ t; s3 = {s3[20:0], s3[31:21]};
            rng_next = {s3, s2, s1, s0};
        end
    endfunction

    // ---- RNG chain / update (stage 0)
    always @(posedge clk) begin
        for (m = 0; m < 8; m = m + 1) begin
            if (r_seedsh) rs[m] <= (m == 7) ? seed_in : rs[(m + 1) & 7];
            else if (rng_en) rs[m] <= rng_next(rs[m]);
        end
    end

    // (MB: the field/spin index j = 8k + m is written as 8 * g_k + m instead of v6.4b's {g_k, m[2:0]}; identical in
    // Verilog and in synthesis, but the Vivado 2025.1 simulator evaluates a concatenated index built from a loop
    // variable as m alone, so v6.4b's text cannot be simulated with xsim.)
    // ---- lane pipeline. z = s*h * 2^16 + C (C = Csame or Cdiff, int32) is kept as {zh, zl}: zl = C[15:0] is a per-step
    // constant, so z < thr and -twoT < z < twoT reduce to ZH-bit signed comparisons of zh plus precomputed low-half
    // comparisons (two's-complement lexicographic order: high parts signed, low halves unsigned).
    reg v1, v2, v3, v4, v5; reg [2:0] k1, k2, k3, k4, k5; reg [1:0] md1, md2, md3, md4, md5;
    reg [HB-1:0] hs1 [0:7]; reg [7:0] si1, spi1; reg [15:0] u1 [0:7];
    reg signed [SW-1:0] sh2 [0:7]; reg [7:0] sel2, b2; reg signed [26:0] da2 [0:7]; reg signed [15:0] db2 [0:7];
    reg [ZH-1:0] zh3 [0:7]; reg signed [SW-1:0] sh3 [0:7]; reg [7:0] b3, sel3; reg signed [42:0] mreg [0:7];
    reg [ZH-1:0] zh4 [0:7]; reg [7:0] ltU4, eqU4, gtL4, eqL4, sel4; reg signed [SW-1:0] sh4 [0:7]; reg [7:0] b4;
    reg signed [42:0] preg [0:7];
    reg [7:0] fl5, nl5, b5, val5; reg signed [SW-1:0] sh5 [0:7];

    wire [7:0] val4;
    genvar gm2;
    generate for (gm2 = 0; gm2 < 8; gm2 = gm2 + 1) begin : vl
        assign val4[gm2] = (B + 32 * (8 * k4 + gm2)) < r_n;   // runtime n (spins >= n are padding)
    end endgenerate

    // sign-extended high parts (ZH bits) of the per-step constants
    wire signed [15:0] cs_hi16 = r_csame[31:16], cd_hi16 = r_cdiff[31:16], uh16 = r_Uh, lh16 = r_Lh;
    wire signed [ZH-1:0] cs_hi = cs_hi16, cd_hi = cd_hi16, uh_z = uh16, lh_z = lh16;

    reg signed [SW-1:0] sh_c; reg signed [ZH-1:0] shz_c; reg [15:0] zl_c, tl_c; reg signed [ZH-1:0] th_c;
    reg signed [10:0] th11_c; reg signed [26:0] thr_c; reg signed [HB:0] hx_c;
    reg lt_c, nlU_c, nlL_c;
    always @(posedge clk) begin
        // stage 0 -> 1
        v1 <= g_v; k1 <= g_k; md1 <= g_mode;
        for (m = 0; m < 8; m = m + 1) begin
            hs1[m] <= h[8 * g_k + m];
            si1[m] <= s[8 * g_k + m];
            spi1[m] <= sp[8 * g_k + m];
            u1[m] <= rng_out(rs[m]) >> 16;
        end
        // stage 1 -> 2: s*h (v6.4b: s*(h - s) + 1 with its stored diagonal -1; identical value)
        v2 <= v1; k2 <= k1; md2 <= md1;
        for (m = 0; m < 8; m = m + 1) begin
            hx_c = $signed({hs1[m][HB-1], hs1[m]});
            sh_c = si1[m] ? hx_c : -hx_c;
            sh2[m] <= sh_c;
            sel2[m] <= (si1[m] == spi1[m]);
            b2[m] <= u1[m][15];
            da2[m] <= $signed(r_fourT);
            db2[m] <= $signed({~u1[m][15], u1[m][14:0]});
        end
        // stage 2 -> 3: zh = s*h + C[31:16] (ZH-bit wrap = high part of the sum, since s*h*2^16 has a zero low half)
        v3 <= v2; k3 <= k2; md3 <= md2;
        for (m = 0; m < 8; m = m + 1) begin
            shz_c = sh2[m];
            zh3[m] <= shz_c + (sel2[m] ? cs_hi : cd_hi);
            sh3[m] <= sh2[m];
            b3[m] <= b2[m];
            sel3[m] <= sel2[m];
            mreg[m] <= da2[m] * db2[m];
        end
        // stage 3 -> 4: high-part comparisons with U = twoT - 1 and L = 1 - twoT
        v4 <= v3; k4 <= k3; md4 <= md3;
        for (m = 0; m < 8; m = m + 1) begin
            zh4[m] <= zh3[m];
            ltU4[m] <= $signed(zh3[m]) < uh_z;
            eqU4[m] <= $signed(zh3[m]) == uh_z;
            gtL4[m] <= $signed(zh3[m]) > lh_z;
            eqL4[m] <= $signed(zh3[m]) == lh_z;
            sel4[m] <= sel3[m];
            sh4[m] <= sh3[m];
            b4[m] <= b3[m];
            preg[m] <= mreg[m];
        end
        // stage 4 -> 5: thr = (prod >> 16) as int32 = {th, tl}; flip = z < thr; n_lin: L <= z <= U
        v5 <= v4; k5 <= k4; md5 <= md4;
        for (m = 0; m < 8; m = m + 1) begin
            thr_c = preg[m] >>> 16;
            th11_c = thr_c[26:16];
            th_c = th11_c;
            tl_c = thr_c[15:0];
            zl_c = sel4[m] ? r_csame[15:0] : r_cdiff[15:0];
            lt_c = ($signed(zh4[m]) < th_c) | (($signed(zh4[m]) == th_c) & (zl_c < tl_c));
            nlU_c = ltU4[m] | (eqU4[m] & (sel4[m] ? r_flags[0] : r_flags[1]));
            nlL_c = gtL4[m] | (eqL4[m] & (sel4[m] ? r_flags[2] : r_flags[3]));
            fl5[m] <= val4[m] & lt_c;
            nl5[m] <= val4[m] & nlU_c & nlL_c;
            val5[m] <= val4[m];
            b5[m] <= b4[m];
            sh5[m] <= sh4[m];
        end
    end

    // ---- stage 5: write-back, extraction words, n_lin accumulator; spin/output registers
    reg [3:0] nl_pc;
    always @(posedge clk) begin
        ex_v <= v5 & (md5 != M_SCORE);
        ex_k <= k5;
        if (g_v && g_k == 3'd0) nl_acc <= 7'd0;
        if (v5 && md5 == M_DEC) begin
            nl_pc = 4'd0;
            for (m = 0; m < 8; m = m + 1) begin
                sp[8 * k5 + m] <= s[8 * k5 + m];
                s[8 * k5 + m] <= s[8 * k5 + m] ^ fl5[m];
                exb[m] <= fl5[m];
                exs[m] <= ~s[8 * k5 + m];
                nl_pc = nl_pc + {3'd0, nl5[m]};
            end
            nl_acc <= nl_acc + {3'd0, nl_pc};
        end else if (v5 && md5 == M_INIT) begin
            for (m = 0; m < 8; m = m + 1) begin
                s[8 * k5 + m] <= val5[m] ? b5[m] : 1'b1;
                exb[m] <= val5[m];
                exs[m] <= b5[m];
            end
        end
        if (r_spset) sp <= {64{1'b1}};
        else if (r_outload) sp <= s;
        else if (r_outshift) sp <= {2'b00, sp[63:2]};
    end

    // ---- v6.1: score accumulation pipelined (pairs, quads, accumulate): sc_acc final 2 cycles after stage 5
    reg v6s, v7s; reg signed [SW:0] scp [0:3]; reg signed [SW+1:0] scq [0:1];
    always @(posedge clk) begin
        v6s <= v5 & (md5 == M_SCORE);
        for (m = 0; m < 4; m = m + 1)
            scp[m] <= (val5[2 * m] ? {sh5[2 * m][SW-1], sh5[2 * m]} : {(SW+1){1'b0}}) +
                      (val5[2 * m + 1] ? {sh5[2 * m + 1][SW-1], sh5[2 * m + 1]} : {(SW+1){1'b0}});
        v7s <= v6s;
        scq[0] <= scp[0] + scp[1]; scq[1] <= scp[2] + scp[3];
        if (g_v && g_k == 3'd0) sc_acc <= {SCG{1'b0}};
        else if (v7s) sc_acc <= sc_acc + {{(SCG-SW-2){scq[0][SW+1]}}, scq[0]} + {{(SCG-SW-2){scq[1][SW+1]}}, scq[1]};
    end

    // ---- couplings (ASIC): the residue-bank plane outputs come from the macros at the top level (the macro output takes
    // the place of the FPGA read register qr) and are registered once more here (= the FPGA output register q), so cq
    // has the FPGA core's latency: address register (ca) -> macro -> cq -> field update.
    reg [8 * 64 * K - 1:0] cq;
    always @(posedge clk) cq <= cqm;

    // ---- field update: control aligned with the coupling words (two register stages after ga)
    reg [7:0] c3v, c3s;
    reg [7:0] P4, Q4, V4;
    // v6.4b: the accumulate mode (init m=1 / decide m=2) only changes at batch boundaries, when no flip is in flight
    // (quiet + DRAIN), and the first flip of a batch reaches this stage >= 12 cycles after the change. So the mode used
    // here may lag by a few cycles: it comes through a two-level tree (gm -> 8 copies -> per-8-field copies).
    reg [7:0] md_a; reg [7:0] md_b;
    always @(posedge clk) begin
        c3v <= gv; c3s <= gs;
        md_a <= {8{gm}}; md_b <= md_a;
        P4 <= c3v & c3s; Q4 <= c3v & ~c3s; V4 <= c3v;
    end
    // per magnitude plane b: a = #{e: m_b & tau_e = +1} (tau_e = +1 iff sign bit 0 and new spin +1, or sign bit 1 and new
    // spin -1), c = #{e: m_b}; pend = sum_b 2^(b-1) * (mode ? 4a - 2c : 2a - c) in PW-bit two's complement (exact: the
    // result lies in [-2^(PW-1), 2^(PW-1)) and intermediate wrap is exact modulo 2^PW)
    reg [3:0] pa_c, pm_c; reg [PW-1:0] pt_c, pv_c;
    integer e, bb;
    always @(posedge clk) begin
        for (j = 0; j < 64; j = j + 1) begin
            pv_c = {PW{1'b0}};
            for (bb = 1; bb < K; bb = bb + 1) begin
                pa_c = 4'd0; pm_c = 4'd0;
                for (e = 0; e < 8; e = e + 1) begin
                    pa_c = pa_c + {3'd0, cq[64 * (K * e + bb) + j] & (cq[64 * (K * e) + j] ? Q4[e] : P4[e])};
                    pm_c = pm_c + {3'd0, cq[64 * (K * e + bb) + j] & V4[e]};
                end
                pt_c = md_b[j >> 3] ? (({{K{1'b0}}, pa_c} << 2) - ({{K{1'b0}}, pm_c} << 1))
                                    : (({{K{1'b0}}, pa_c} << 1) - {{K{1'b0}}, pm_c});
                pv_c = pv_c + (pt_c << (bb - 1));
            end
            pend[j] <= pv_c;
            h[j] <= (r_blv && r_blj == j) ? {{(HB-16){r_bld[15]}}, r_bld} :      // bias load (seed phase)
                    r_hclr ? {HB{1'b0}} : h[j] + {{(HB-PW){pend[j][PW-1]}}, pend[j]};
        end
    end
endmodule
`default_nettype wire
