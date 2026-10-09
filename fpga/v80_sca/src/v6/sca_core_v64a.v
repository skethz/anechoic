// sca_core.v -- V80 Onsager-SCA engine core, v6.4 (hand-written RTL, streaming I/O).
// v6.4: 8 flip extractors and 4 coupling copies (2 BRAM + 2 URAM, 8 read ports): extraction was 64% of step cycles with 4.
// Extractor e drains groups e, e+8, e+16, e+24 (32-bit round words, bit p = 4m + q'); flipped spin x = {k, p, e}.
// v6.2 (4 extractors) is kept in sca_core_v62.v.
// v6.2: header bit 34 (keep_tables) skips the table stream; the tables of the previous launch stay resident (v6.1 is kept in
// sca_core_v61.v). In multi-engine rounds the per-launch table stream was NoC-bound (~20 cycles/step on the V80).
// Bit-exact with ../sca_ref.hpp built with SCA_LANES=256 (same RNG streams, decisions, fields, flips, score, n_lin trace).
// v6.0 (first verified version) is kept in sca_core_v60.v; v6.1 changes timing structure only (see "v6.1" notes):
// replicated enables/constants, 16-bit split comparisons in the lanes, registered coupling-address mux, pipelined score
// path, input skid buffer, precomputed end-of-run flag. Cycles per step are unchanged.
//
// Why v6: the HLS engines (v4-v5.3) outline every pipelined loop into its own module and pass the register-resident state
// (fields, RNG states, spin words) in and out, so the state exists two or three times (355k FF per engine) and the
// loop-exit handshakes drive 32k-43k-fanout clock enables through BUFG_FABRIC (the limiting paths at 250-300 MHz).
// Here every state bit exists once and there is no AXI handshake anywhere in the datapath.
//
// Organisation: 32 groups. Group b owns the 64 spins/fields y = b + 32j (j = 0..63), the 8 decision lanes l = b + 32m
// (m = 0..7; lane l handles spins l + 256k in round k, i.e. j = 8k + m), and its 64 columns of both coupling copies
// (bank b, address x: bit j = J[x][b + 32j]). So the field update, the lane field multiplexers and the coupling banks
// are local to a group. Four flip extractors (e = 0..3) drain the flip bits of groups e, e+4, ..., e+28, one 64-bit word
// per decision round (bit p = 8m + q for group e + 4q, lane m), so extraction overlaps with the remaining rounds; the
// flipped spin is x = {k, p, e} = 256k + 32m + 4q + e.
//
// Stream protocol (input 128-bit beats, output 64-bit words):
//   in : H0 = {seed[63:0], 29'b0, keep_tables, do_trace, load_J, S[15:0], trials[15:0]}  H1 = {sumw[63:0], 32'b0, trial_offset[31:0]}
//        if load_J: 32768 beats of couplings, beat = row x, columns 128c..128c+127 (c = 0..15), bit = 1 iff J = +1
//        unless keep_tables: S beats of tables {kcorr[31:0], kconst[31:0], q[31:0], fourT[31:0]}  (kcorr must fit in int32)
//   out: per trial 32 spin words (bit = 1 iff spin = +1, spins >= N zero), cut (int64), flips (uint64),
//        then, if do_trace, S words {32'b0, n_lin}.
`timescale 1ns/1ps
`default_nettype none

module sca_core #(
    parameter integer N = 2000,
    parameter integer DRAIN = 2   // cycles between "extractors quiet" and the next round command (>= 2; 1 fails in simulation)
) (
    (* X_INTERFACE_INFO = "xilinx.com:signal:clock:1.0 clk CLK" *)
    (* X_INTERFACE_PARAMETER = "ASSOCIATED_BUSIF s_axis:m_axis, ASSOCIATED_RESET rst_n" *)
    input  wire         clk,
    (* X_INTERFACE_INFO = "xilinx.com:signal:reset:1.0 rst_n RST" *)
    (* X_INTERFACE_PARAMETER = "POLARITY ACTIVE_LOW" *)
    input  wire         rst_n,
    input  wire [127:0] s_axis_tdata,
    input  wire         s_axis_tvalid,
    output wire         s_axis_tready,
    output wire [63:0]  m_axis_tdata,
    output wire         m_axis_tvalid,
    input  wire         m_axis_tready
);
    localparam [3:0] ST_HDR0 = 4'd0, ST_HDR1 = 4'd1, ST_LOADJ = 4'd2, ST_LOADT = 4'd3, ST_TRIAL = 4'd4,
                     ST_SEED = 4'd5, ST_BATCH = 4'd6, ST_WAIT = 4'd7, ST_SCOREW = 4'd8, ST_OUT = 4'd9;
    localparam [1:0] M_INIT = 2'd0, M_DEC = 2'd1, M_SCORE = 2'd2;
    localparam integer SEED_LAT = 15;   // sd_run (aligned with sd_xa) -> sca_mix output

    reg [3:0] st;
    reg [15:0] p_trials; reg [15:0] p_S; reg p_loadJ, p_trace, p_keepT; reg [63:0] p_seed; reg [31:0] p_toff;
    reg signed [63:0] p_sumw;
    reg [15:0] tr;           // trial index within the launch
    reg [15:0] t_step;       // number of DEC batches started in this trial (= current step + 1)
    reg [15:0] cnt;          // load counter / generic counter

    // ------------------------------------------------------------------ v6.1: input skid buffer (registered tready)
    reg [127:0] in_d, sk_d; reg in_v, sk_v;
    wire core_ready = (st == ST_HDR0) | (st == ST_HDR1) | (st == ST_LOADJ) | (st == ST_LOADT);
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

    // ------------------------------------------------------------------ coupling load (row assembly)
    // rowbuf shifts with every input beat (its content only matters in LOADJ, where every valid beat is consumed);
    // the write strobe is registered, and the groups take the completed row from rowbuf one cycle later.
    reg [2047:0] rowbuf;
    always @(posedge clk) if (in_v) rowbuf <= {in_d, rowbuf[2047:128]};
    (* max_fanout = 8 *) reg j_we_r; (* max_fanout = 8 *) reg [10:0] j_addr_r;
    always @(posedge clk) begin
        j_we_r <= (st == ST_LOADJ) & in_v & (cnt[3:0] == 4'd15);
        j_addr_r <= cnt[14:4];
    end

    // ------------------------------------------------------------------ tables and trace
    reg [127:0] tab_mem [0:4095];
    reg [11:0] tab_raddr; reg [127:0] tab_q1, tab_q;
    always @(posedge clk) begin
        if ((st == ST_LOADT) & in_v) tab_mem[cnt[11:0]] <= in_d;
        tab_q1 <= tab_mem[tab_raddr];
        tab_q <= tab_q1;
    end
    reg [11:0] trace_mem [0:4095];
    reg trace_we; reg [11:0] trace_wa; reg [11:0] trace_wd;
    reg [11:0] trace_ra; reg [11:0] trace_q;
    always @(posedge clk) begin
        if (trace_we) trace_mem[trace_wa] <= trace_wd;
        trace_q <= trace_mem[trace_ra];
    end

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
    sca_mix mixa (.clk(clk), .z(sd_xa), .o(sd_a));
    sca_mix mixb (.clk(clk), .z(sd_xb), .o(sd_b));
    reg [SEED_LAT-1:0] sd_vpipe;
    always @(posedge clk) sd_vpipe <= rst_n ? {sd_vpipe[SEED_LAT-2:0], sd_run} : {SEED_LAT{1'b0}};
    wire sd_ov = sd_vpipe[SEED_LAT-1];
    wire [127:0] sd_state_raw = {sd_b, sd_a};
    wire [127:0] sd_state = (sd_state_raw == 128'd0) ? 128'd1 : sd_state_raw;
    (* max_fanout = 8 *) reg seed_shift_c; reg [127:0] seed_d1, seed_d2;
    always @(posedge clk) begin
        seed_shift_c <= sd_ov;
        seed_d1 <= sd_state;
        seed_d2 <= seed_d1;
    end

    // ------------------------------------------------------------------ broadcast control registers (v6.1: replicated)
    (* max_fanout = 8 *) reg lc_v; (* max_fanout = 8 *) reg [2:0] lc_k; (* max_fanout = 8 *) reg [1:0] lc_mode;
    (* max_fanout = 8 *) reg [26:0] k_fourT; (* max_fanout = 8 *) reg [31:0] k_csame, k_cdiff;
    (* max_fanout = 8 *) reg [15:0] k_Uh, k_Lh; (* max_fanout = 8 *) reg [3:0] k_flags;
    (* max_fanout = 8 *) reg h_clr_c, sp_set_c, out_load_c, out_shift_c, amode;

    // ------------------------------------------------------------------ groups
    wire [127:0] chain [0:32];
    assign chain[32] = seed_d2;
    wire [31:0] g_exv; wire [95:0] g_exk; wire [255:0] g_exb, g_exs;
    wire [223:0] g_nl; wire [607:0] g_sc; wire [63:0] g_ob;
    wire [87:0] ax_bus; wire [7:0] av_bus, as_bus;

    genvar gb;
    generate for (gb = 0; gb < 32; gb = gb + 1) begin : grp
        wire [63:0] jd;
        genvar gj;
        for (gj = 0; gj < 64; gj = gj + 1) begin : jcol
            assign jd[gj] = rowbuf[gb + 32 * gj];
        end
        sca_group #(.B(gb), .N(N)) u (
            .clk(clk), .lc_v(lc_v), .lc_k(lc_k), .lc_mode(lc_mode),
            .k_fourT(k_fourT), .k_csame(k_csame), .k_cdiff(k_cdiff), .k_Uh(k_Uh), .k_Lh(k_Lh), .k_flags(k_flags),
            .seed_shift(seed_shift_c), .seed_in(chain[gb + 1]), .seed_out(chain[gb]),
            .h_clr(h_clr_c), .sp_set(sp_set_c),
            .ax(ax_bus), .av(av_bus), .asg(as_bus), .amode(amode),
            .jwe(j_we_r), .jaddr(j_addr_r), .jdata(jd),
            .ex_v(g_exv[gb]), .ex_k(g_exk[3 * gb +: 3]), .exb(g_exb[8 * gb +: 8]), .exs(g_exs[8 * gb +: 8]),
            .nl_acc(g_nl[7 * gb +: 7]), .sc_acc(g_sc[19 * gb +: 19]),
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
        sca_extractor #(.E(gb)) u (
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
    // v6.1: score tree, one adder level per register stage (32 -> 16 -> 8 -> 4 -> 1)
    reg sA, sB, sC, sD; reg signed [19:0] scA [0:15]; reg signed [20:0] scB [0:7]; reg signed [21:0] scC [0:3];
    reg signed [23:0] sc_tot;
    always @(posedge clk) begin
        sA <= cap_sr[9] & (capm_sr[19:18] == M_SCORE);
        for (ii = 0; ii < 16; ii = ii + 1)
            scA[ii] <= $signed(g_sc[19 * (2 * ii) +: 19]) + $signed(g_sc[19 * (2 * ii + 1) +: 19]);
        sB <= sA; for (ii = 0; ii < 8; ii = ii + 1) scB[ii] <= scA[2 * ii] + scA[2 * ii + 1];
        sC <= sB; for (ii = 0; ii < 4; ii = ii + 1) scC[ii] <= scB[2 * ii] + scB[2 * ii + 1];
        sD <= sC; sc_tot <= scC[0] + scC[1] + scC[2] + scC[3];
    end

    // next-step table entry and the correction pipeline: corr = int32(((nlin * kcorr) >> 8) + kconst)
    reg [31:0] tn_fourT, tn_q, tn_kconst, tn_kcorr;
    reg [31:0] tn_U, tn_L;                     // v6.1: U = twoT - 1, L = 1 - twoT (int32), from the prefetched entry
    reg tab_ld3, tab_ld2, tab_ld1, tab_ld4;    // table read in flight: address with tab_ld1, data in tab_q with tab_ld3
    wire [31:0] tn_twoT = {1'b0, tn_fourT[31:1]};
    always @(posedge clk) begin
        tab_ld2 <= tab_ld1; tab_ld3 <= tab_ld2; tab_ld4 <= tab_ld3;
        if (tab_ld3) begin
            tn_fourT <= tab_q[31:0]; tn_q <= tab_q[63:32]; tn_kconst <= tab_q[95:64]; tn_kcorr <= tab_q[127:96];
        end
        if (tab_ld4) begin tn_U <= tn_twoT - 32'd1; tn_L <= 32'd1 - tn_twoT; end
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
    wire signed [63:0] sh64 = {{40{sc_tot[23]}}, sc_tot};
    wire signed [63:0] cut_sum = p_sumw + cut_half;
    always @(posedge clk) begin
        cut_go1 <= sD;
        cut_half <= (sh64 + {63'd0, sh64[63]}) >>> 1;
        cut_go2 <= cut_go1;
        cut_val <= (cut_sum + {63'd0, cut_sum[63]}) >>> 1;
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
    reg trace_rd; reg trace_rd_d;

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

        if (~rst_n) begin
            st <= ST_HDR0; cnt <= 16'd0; nk_ready <= 1'b0; cap_sr <= 10'd0; h_clr_c <= 1'b0; sp_set_c <= 1'b0;
            amode <= 1'b0;
        end else case (st)
        ST_HDR0: if (in_fire) begin
            p_trials <= in_d[15:0]; p_S <= in_d[31:16];
            p_loadJ <= in_d[32]; p_trace <= in_d[33]; p_keepT <= in_d[34]; p_seed <= in_d[127:64];
            st <= ST_HDR1;
        end
        ST_HDR1: if (in_fire) begin
            p_toff <= in_d[31:0]; p_sumw <= in_d[127:64];
            cnt <= 16'd0;
            st <= p_loadJ ? ST_LOADJ : ((p_S != 16'd0 && !p_keepT) ? ST_LOADT : ST_TRIAL);
            tr <= 16'd0;
        end
        ST_LOADJ: if (in_fire) begin
            cnt <= cnt + 16'd1;
            if (cnt == 16'd32767) begin cnt <= 16'd0; st <= (p_S != 16'd0 && !p_keepT) ? ST_LOADT : ST_TRIAL; end
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
                tab_raddr <= 12'd0; tab_ld1 <= 1'b1;      // step 0 entry
                nk_ready <= 1'b0;
            end
        end
        ST_SEED: begin
            // issue 256 seeder inputs, then wait until the last state is in the chain (2 cycles after sd_ov)
            if (sd_p != 9'd256) begin sd_run <= 1'b1; sd_p <= sd_p + 9'd1; end
            if (seed_shift_c) cnt <= cnt + 16'd1;
            if (cnt == 16'd256 && ~seed_shift_c) begin
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
            if (oc < 8'd3) oc <= oc + 8'd1;
            else if (oc < 8'd35) begin
                if (ocad == 2'd0) begin
                    if (of_room) begin
                        of_push <= 1'b1;
                        of_din <= (oc == 8'd34) ? (spin_word & 64'h0000_0000_0000_FFFF) : spin_word;
                        out_shift_c <= 1'b1; ocad <= 2'd1; oc <= oc + 8'd1;
                    end
                end else ocad <= (ocad == 2'd2) ? 2'd0 : ocad + 2'd1;
            end else if (oc == 8'd35) begin
                if (of_room) begin of_push <= 1'b1; of_din <= cut_val; oc <= oc + 8'd1; end
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
endmodule

// --------------------------------------------------------------------------------------------------------------------
// splitmix64 output mix, pipelined: o = mix(z) after 15 cycles (z^z>>30, *C1, z^z>>27, *C2, z^z>>31).
module sca_mix (
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
module sca_extractor #(
    parameter integer E = 0
) (
    input  wire        clk,
    input  wire        rst,
    input  wire        in_v,
    input  wire [2:0]  in_k,
    input  wire [31:0] in_b,
    input  wire [31:0] in_s,
    (* max_fanout = 8 *) output reg o_v,
    (* max_fanout = 8 *) output reg [10:0] o_x,
    (* max_fanout = 8 *) output reg o_s,
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
// Group b: 64 spins/fields y = b + 32j, 8 decision lanes (lane m -> spins j = 8k + m), coupling bank b (two copies).
module sca_group #(
    parameter integer B = 0,
    parameter integer N = 2000
) (
    input  wire         clk,
    input  wire         lc_v,
    input  wire [2:0]   lc_k,
    input  wire [1:0]   lc_mode,
    input  wire [26:0]  k_fourT,
    input  wire [31:0]  k_csame,
    input  wire [31:0]  k_cdiff,
    input  wire [15:0]  k_Uh,
    input  wire [15:0]  k_Lh,
    input  wire [3:0]   k_flags,
    input  wire         seed_shift,
    input  wire [127:0] seed_in,
    output wire [127:0] seed_out,
    input  wire         h_clr,
    input  wire         sp_set,
    input  wire [87:0]  ax,
    input  wire [7:0]   av,
    input  wire [7:0]   asg,
    input  wire         amode,
    input  wire         jwe,
    input  wire [10:0]  jaddr,
    input  wire [63:0]  jdata,
    output reg          ex_v,
    output reg  [2:0]   ex_k,
    output reg  [7:0]   exb,
    output reg  [7:0]   exs,
    output reg  [6:0]   nl_acc,
    output reg  [18:0]  sc_acc,
    input  wire         out_load,
    input  wire         out_shift,
    output wire [1:0]   out_bits
);
    localparam [1:0] M_INIT = 2'd0, M_DEC = 2'd1, M_SCORE = 2'd2;
    integer m, j;

    // ---- one register hop for every broadcast input (v6.1: replicated toward their loads)
    reg g_v; (* max_fanout = 16 *) reg [2:0] g_k; (* max_fanout = 16 *) reg [1:0] g_mode;
    (* max_fanout = 4 *) reg [26:0] r_fourT; (* max_fanout = 4 *) reg [31:0] r_csame, r_cdiff;
    (* max_fanout = 4 *) reg [15:0] r_Uh, r_Lh; (* max_fanout = 4 *) reg [3:0] r_flags;
    (* max_fanout = 32 *) reg rng_en;                   // = g_v & (g_mode != SCORE), registered directly
    (* max_fanout = 4 *) reg [10:0] ga0, ga2, ga4, ga6;   // port A addresses (coupling write address during the load)
    (* max_fanout = 4 *) reg [10:0] ga1, ga3, ga5, ga7;
    reg [7:0] gv, gs; reg gm;
    (* max_fanout = 64 *) reg r_hclr; (* max_fanout = 16 *) reg r_spset;
    (* max_fanout = 32 *) reg r_seedsh; (* max_fanout = 16 *) reg r_outload, r_outshift;
    (* max_fanout = 4 *) reg r_jwe; reg [63:0] r_jdata;
    always @(posedge clk) begin
        g_v <= lc_v; g_k <= lc_k; g_mode <= lc_mode;
        rng_en <= lc_v & (lc_mode != M_SCORE);
        r_fourT <= k_fourT; r_csame <= k_csame; r_cdiff <= k_cdiff; r_Uh <= k_Uh; r_Lh <= k_Lh; r_flags <= k_flags;
        ga0 <= jwe ? jaddr : ax[10:0]; ga1 <= ax[21:11]; ga2 <= jwe ? jaddr : ax[32:22]; ga3 <= ax[43:33];
        ga4 <= jwe ? jaddr : ax[54:44]; ga5 <= ax[65:55]; ga6 <= jwe ? jaddr : ax[76:66]; ga7 <= ax[87:77];
        gv <= av; gs <= asg; gm <= amode;
        r_hclr <= h_clr; r_spset <= sp_set; r_seedsh <= seed_shift; r_outload <= out_load; r_outshift <= out_shift;
        r_jwe <= jwe; r_jdata <= jdata;
    end

    // ---- state
    reg [11:0] h [0:63];        // stored field = true field - s (|value| <= 2000)
    reg signed [5:0] pend [0:63];
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

    // ---- lane pipeline. z = (s*h + 1) * 2^16 + C (C = Csame or Cdiff, int32) is kept as {zh, zl}: zl = C[15:0] is a
    // per-step constant, so z < thr and -twoT < z < twoT reduce to 16-bit signed comparisons of zh plus precomputed
    // low-half comparisons (two's-complement lexicographic order: high halves signed, low halves unsigned).
    reg v1, v2, v3, v4, v5; reg [2:0] k1, k2, k3, k4, k5; reg [1:0] md1, md2, md3, md4, md5;
    reg [11:0] hs1 [0:7]; reg [7:0] si1, spi1; reg [15:0] u1 [0:7];
    reg signed [12:0] sh2 [0:7]; reg [7:0] sel2, b2; reg signed [26:0] da2 [0:7]; reg signed [15:0] db2 [0:7];
    reg [15:0] zh3 [0:7]; reg signed [12:0] sh3 [0:7]; reg [7:0] b3, sel3; reg signed [42:0] mreg [0:7];
    reg [15:0] zh4 [0:7]; reg [7:0] ltU4, eqU4, gtL4, eqL4, sel4; reg signed [12:0] sh4 [0:7]; reg [7:0] b4;
    reg signed [42:0] preg [0:7];
    reg [7:0] fl5, nl5, b5, val5; reg signed [12:0] sh5 [0:7];

    wire [7:0] val4;
    genvar gm2;
    generate for (gm2 = 0; gm2 < 8; gm2 = gm2 + 1) begin : vl
        assign val4[gm2] = (B + 32 * (8 * k4 + gm2)) < N;
    end endgenerate

    reg signed [12:0] sh_c; reg [15:0] zh_c, zl_c, th_c, tl_c; reg signed [26:0] thr_c;
    reg lt_c, nlU_c, nlL_c;
    always @(posedge clk) begin
        // stage 0 -> 1
        v1 <= g_v; k1 <= g_k; md1 <= g_mode;
        for (m = 0; m < 8; m = m + 1) begin
            hs1[m] <= h[{g_k, m[2:0]}];
            si1[m] <= s[{g_k, m[2:0]}];
            spi1[m] <= sp[{g_k, m[2:0]}];
            u1[m] <= rng_out(rs[m]) >> 16;
        end
        // stage 1 -> 2
        v2 <= v1; k2 <= k1; md2 <= md1;
        for (m = 0; m < 8; m = m + 1) begin
            sh_c = si1[m] ? ($signed({hs1[m][11], hs1[m]}) + 13'sd1) : (13'sd1 - $signed({hs1[m][11], hs1[m]}));
            sh2[m] <= sh_c;
            sel2[m] <= (si1[m] == spi1[m]);
            b2[m] <= u1[m][15];
            da2[m] <= $signed(r_fourT);
            db2[m] <= $signed({~u1[m][15], u1[m][14:0]});
        end
        // stage 2 -> 3: zh = sh + C[31:16] (16-bit wrap = high half of the int32 sum, since sh*2^16 has a zero low half)
        v3 <= v2; k3 <= k2; md3 <= md2;
        for (m = 0; m < 8; m = m + 1) begin
            zh3[m] <= {{3{sh2[m][12]}}, sh2[m]} + (sel2[m] ? r_csame[31:16] : r_cdiff[31:16]);
            sh3[m] <= sh2[m];
            b3[m] <= b2[m];
            sel3[m] <= sel2[m];
            mreg[m] <= da2[m] * db2[m];
        end
        // stage 3 -> 4: high-half comparisons with U = twoT - 1 and L = 1 - twoT
        v4 <= v3; k4 <= k3; md4 <= md3;
        for (m = 0; m < 8; m = m + 1) begin
            zh4[m] <= zh3[m];
            ltU4[m] <= $signed(zh3[m]) < $signed(r_Uh);
            eqU4[m] <= zh3[m] == r_Uh;
            gtL4[m] <= $signed(zh3[m]) > $signed(r_Lh);
            eqL4[m] <= zh3[m] == r_Lh;
            sel4[m] <= sel3[m];
            sh4[m] <= sh3[m];
            b4[m] <= b3[m];
            preg[m] <= mreg[m];
        end
        // stage 4 -> 5: thr = (prod >> 16) as int32 = {th, tl}; flip = z < thr; n_lin: L <= z <= U
        v5 <= v4; k5 <= k4; md5 <= md4;
        for (m = 0; m < 8; m = m + 1) begin
            thr_c = preg[m] >>> 16;
            th_c = {{5{thr_c[26]}}, thr_c[26:16]};
            tl_c = thr_c[15:0];
            zl_c = sel4[m] ? r_csame[15:0] : r_cdiff[15:0];
            lt_c = ($signed(zh4[m]) < $signed(th_c)) | ((zh4[m] == th_c) & (zl_c < tl_c));
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
                sp[{k5, m[2:0]}] <= s[{k5, m[2:0]}];
                s[{k5, m[2:0]}] <= s[{k5, m[2:0]}] ^ fl5[m];
                exb[m] <= fl5[m];
                exs[m] <= ~s[{k5, m[2:0]}];
                nl_pc = nl_pc + {3'd0, nl5[m]};
            end
            nl_acc <= nl_acc + {3'd0, nl_pc};
        end else if (v5 && md5 == M_INIT) begin
            for (m = 0; m < 8; m = m + 1) begin
                s[{k5, m[2:0]}] <= val5[m] ? b5[m] : 1'b1;
                exb[m] <= val5[m];
                exs[m] <= b5[m];
            end
        end
        if (r_spset) sp <= {64{1'b1}};
        else if (r_outload) sp <= s;
        else if (r_outshift) sp <= {2'b00, sp[63:2]};
    end

    // ---- v6.1: score accumulation pipelined (pairs, quads, accumulate): sc_acc final 2 cycles after stage 5
    reg v6s, v7s; reg signed [13:0] scp [0:3]; reg signed [14:0] scq [0:1];
    always @(posedge clk) begin
        v6s <= v5 & (md5 == M_SCORE);
        for (m = 0; m < 4; m = m + 1)
            scp[m] <= (val5[2 * m] ? {sh5[2 * m][12], sh5[2 * m]} : 14'sd0) +
                      (val5[2 * m + 1] ? {sh5[2 * m + 1][12], sh5[2 * m + 1]} : 14'sd0);
        v7s <= v6s;
        scq[0] <= scp[0] + scp[1]; scq[1] <= scp[2] + scp[3];
        if (g_v && g_k == 3'd0) sc_acc <= 19'd0;
        else if (v7s) sc_acc <= sc_acc + {{4{scq[0][14]}}, scq[0]} + {{4{scq[1][14]}}, scq[1]};
    end

    // ---- couplings: two copies x two read ports; port A of each copy also writes during the load
    // (v6.1: the write/read address mux is in front of the ga0/ga2 registers, so every BRAM address pin is driven by a
    // register; jwe/jaddr arrive one cycle ahead of jdata's registered copy, matching r_jwe/r_jdata.)
    reg [63:0] mem0 [0:2047];
    reg [63:0] mem1 [0:2047];
    reg [63:0] q0a_r, q0b_r, q1a_r, q1b_r, q0a, q0b, q1a, q1b;
    always @(posedge clk) begin
        if (r_jwe) mem0[ga0] <= r_jdata;
        q0a_r <= mem0[ga0];
        q0a <= q0a_r;
    end
    always @(posedge clk) begin
        q0b_r <= mem0[ga1];
        q0b <= q0b_r;
    end
    always @(posedge clk) begin
        if (r_jwe) mem1[ga2] <= r_jdata;
        q1a_r <= mem1[ga2];
        q1a <= q1a_r;
    end
    always @(posedge clk) begin
        q1b_r <= mem1[ga3];
        q1b <= q1b_r;
    end
    // v6.4: copies 2 and 3 in UltraRAM, read latency 2 like the BRAM copies. Synthesis uses the XPM macro (URAM inference
    // from behavioural code was rejected: "invalid write mode"); simulation uses the equivalent behavioural model
    // (no_change write mode: port A output holds while writing, which only happens during the coupling load).
    wire [63:0] q2a, q2b, q3a, q3b;
`ifdef SYNTHESIS
    xpm_memory_tdpram #(
        .ADDR_WIDTH_A(11), .ADDR_WIDTH_B(11), .BYTE_WRITE_WIDTH_A(64), .BYTE_WRITE_WIDTH_B(64), .CLOCKING_MODE("common_clock"),
        .MEMORY_PRIMITIVE("ultra"), .MEMORY_SIZE(131072), .READ_DATA_WIDTH_A(64), .READ_DATA_WIDTH_B(64),
        .READ_LATENCY_A(2), .READ_LATENCY_B(2), .WRITE_DATA_WIDTH_A(64), .WRITE_DATA_WIDTH_B(64),
        .WRITE_MODE_A("no_change"), .WRITE_MODE_B("no_change"), .USE_MEM_INIT(0), .MEMORY_INIT_FILE("none"),
        .ECC_MODE("no_ecc"), .AUTO_SLEEP_TIME(0), .WAKEUP_TIME("disable_sleep"), .MESSAGE_CONTROL(0),
        .USE_EMBEDDED_CONSTRAINT(0), .SIM_ASSERT_CHK(0)
    ) u_mem2 (
        .clka(clk), .clkb(clk), .ena(1'b1), .enb(1'b1), .regcea(1'b1), .regceb(1'b1), .rsta(1'b0), .rstb(1'b0),
        .wea(r_jwe), .web(1'b0), .addra(ga4), .addrb(ga5), .dina(r_jdata), .dinb(64'd0), .douta(q2a), .doutb(q2b),
        .injectdbiterra(1'b0), .injectsbiterra(1'b0), .injectdbiterrb(1'b0), .injectsbiterrb(1'b0), .sleep(1'b0),
        .dbiterra(), .sbiterra(), .dbiterrb(), .sbiterrb());
    xpm_memory_tdpram #(
        .ADDR_WIDTH_A(11), .ADDR_WIDTH_B(11), .BYTE_WRITE_WIDTH_A(64), .BYTE_WRITE_WIDTH_B(64), .CLOCKING_MODE("common_clock"),
        .MEMORY_PRIMITIVE("ultra"), .MEMORY_SIZE(131072), .READ_DATA_WIDTH_A(64), .READ_DATA_WIDTH_B(64),
        .READ_LATENCY_A(2), .READ_LATENCY_B(2), .WRITE_DATA_WIDTH_A(64), .WRITE_DATA_WIDTH_B(64),
        .WRITE_MODE_A("no_change"), .WRITE_MODE_B("no_change"), .USE_MEM_INIT(0), .MEMORY_INIT_FILE("none"),
        .ECC_MODE("no_ecc"), .AUTO_SLEEP_TIME(0), .WAKEUP_TIME("disable_sleep"), .MESSAGE_CONTROL(0),
        .USE_EMBEDDED_CONSTRAINT(0), .SIM_ASSERT_CHK(0)
    ) u_mem3 (
        .clka(clk), .clkb(clk), .ena(1'b1), .enb(1'b1), .regcea(1'b1), .regceb(1'b1), .rsta(1'b0), .rstb(1'b0),
        .wea(r_jwe), .web(1'b0), .addra(ga6), .addrb(ga7), .dina(r_jdata), .dinb(64'd0), .douta(q3a), .doutb(q3b),
        .injectdbiterra(1'b0), .injectsbiterra(1'b0), .injectdbiterrb(1'b0), .injectsbiterrb(1'b0), .sleep(1'b0),
        .dbiterra(), .sbiterra(), .dbiterrb(), .sbiterrb());
`else
    reg [63:0] mem2 [0:2047];
    reg [63:0] mem3 [0:2047];
    reg [63:0] q2a_r, q2b_r, q3a_r, q3b_r, q2a_q, q2b_q, q3a_q, q3b_q;
    always @(posedge clk) begin
        if (r_jwe) mem2[ga4] <= r_jdata;
        else q2a_r <= mem2[ga4];
        q2a_q <= q2a_r;
    end
    always @(posedge clk) begin
        q2b_r <= mem2[ga5];
        q2b_q <= q2b_r;
    end
    always @(posedge clk) begin
        if (r_jwe) mem3[ga6] <= r_jdata;
        else q3a_r <= mem3[ga6];
        q3a_q <= q3a_r;
    end
    always @(posedge clk) begin
        q3b_r <= mem3[ga7];
        q3b_q <= q3b_r;
    end
    assign q2a = q2a_q; assign q2b = q2b_q; assign q3a = q3a_q; assign q3b = q3b_q;
`endif

    // ---- field update: control aligned with the coupling words (two register stages after ga*)
    reg [7:0] c3v, c3s; reg c3m;
    (* max_fanout = 16 *) reg [7:0] P4, Q4; (* max_fanout = 16 *) reg [4:0] Bc4; (* max_fanout = 16 *) reg dec4;
    always @(posedge clk) begin
        c3v <= gv; c3s <= gs; c3m <= gm;
        P4 <= c3v & c3s; Q4 <= c3v & ~c3s; dec4 <= c3m;
        Bc4 <= ({4'd0, c3v[0]} + {4'd0, c3v[1]} + {4'd0, c3v[2]} + {4'd0, c3v[3]} +
                {4'd0, c3v[4]} + {4'd0, c3v[5]} + {4'd0, c3v[6]} + {4'd0, c3v[7]}) << c3m;
    end
    // pend = m*(2*pc - nvalid) in [-16, 16] (6-bit two's complement; intermediate wrap is exact modulo 64)
    reg [3:0] pc_c; reg [5:0] pv_c;
    always @(posedge clk) begin
        for (j = 0; j < 64; j = j + 1) begin
            pc_c = {3'd0, q0a[j] ? P4[0] : Q4[0]} + {3'd0, q0b[j] ? P4[1] : Q4[1]} +
                   {3'd0, q1a[j] ? P4[2] : Q4[2]} + {3'd0, q1b[j] ? P4[3] : Q4[3]} +
                   {3'd0, q2a[j] ? P4[4] : Q4[4]} + {3'd0, q2b[j] ? P4[5] : Q4[5]} +
                   {3'd0, q3a[j] ? P4[6] : Q4[6]} + {3'd0, q3b[j] ? P4[7] : Q4[7]};
            pv_c = (dec4 ? {pc_c, 2'b00} : {1'b0, pc_c, 1'b0}) - {1'b0, Bc4};
            pend[j] <= pv_c;
            h[j] <= r_hclr ? 12'd0 : h[j] + {{6{pend[j][5]}}, pend[j]};
        end
    end
endmodule
`default_nettype wire
