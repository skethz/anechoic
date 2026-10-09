// sca_core.v -- V80 Onsager-SCA engine core, v6 (hand-written RTL, streaming I/O).
// Bit-exact with ../sca_ref.hpp built with SCA_LANES=256 (same RNG streams, decisions, fields, flips, score, n_lin trace).
//
// Why v6: the HLS engines (v4-v5.3) outline every pipelined loop into its own module and pass the register-resident state
// (fields, RNG states, spin words) in and out, so the state exists two or three times (355k FF per engine) and the
// loop-exit handshakes drive 32k-43k-fanout clock enables through BUFG_FABRIC (the limiting paths at 250-300 MHz).
// Here every state bit exists once, there is no AXI handshake anywhere in the datapath, and all wide broadcasts are
// one register hop into a group.
//
// Organisation: 32 groups. Group b owns the 64 spins/fields y = b + 32j (j = 0..63), the 8 decision lanes l = b + 32m
// (m = 0..7; lane l handles spins l + 256k in round k, i.e. j = 8k + m), and its 64 columns of both coupling copies
// (bank b, address x: bit j = J[x][b + 32j]). So the field update, the lane field multiplexers and the coupling banks
// are local to a group. Four flip extractors (e = 0..3) drain the flip bits of groups e, e+4, ..., e+28, one 64-bit word
// per decision round (bit p = 8m + q for group e + 4q, lane m), so extraction overlaps with the remaining rounds; the
// flipped spin is x = {k, p, e} = 256k + 32m + 4q + e.
//
// Stream protocol (input 128-bit beats, output 64-bit words):
//   in : H0 = {seed[63:0], 30'b0, do_trace, load_J, S[15:0], trials[15:0]}  H1 = {sumw[63:0], 32'b0, trial_offset[31:0]}
//        if load_J: 32768 beats of couplings, beat = row x, columns 128c..128c+127 (c = 0..15), bit = 1 iff J = +1
//        S beats of tables {kcorr[31:0], kconst[31:0], q[31:0], fourT[31:0]}  (kcorr must fit in int32)
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
    reg [15:0] p_trials; reg [15:0] p_S; reg p_loadJ, p_trace; reg [63:0] p_seed; reg [31:0] p_toff;
    reg signed [63:0] p_sumw;
    reg [15:0] tr;           // trial index within the launch
    reg [15:0] t_step;       // step index of the current DEC batch
    reg [15:0] cnt;          // load counter / generic counter

    assign s_axis_tready = (st == ST_HDR0) | (st == ST_HDR1) | (st == ST_LOADJ) | (st == ST_LOADT);
    wire in_fire = s_axis_tvalid & s_axis_tready;

    // ------------------------------------------------------------------ coupling load (row assembly)
    reg [2047:0] rowbuf;
    wire [2047:0] rowfull = {s_axis_tdata, rowbuf[2047:128]};
    wire j_we_c = (st == ST_LOADJ) & s_axis_tvalid & (cnt[3:0] == 4'd15);
    wire [10:0] j_addr_c = cnt[14:4];
    always @(posedge clk) if ((st == ST_LOADJ) & s_axis_tvalid) rowbuf <= rowfull;

    // ------------------------------------------------------------------ tables and trace
    reg [127:0] tab_mem [0:4095];
    reg [11:0] tab_raddr; reg [127:0] tab_q1, tab_q;
    always @(posedge clk) begin
        if ((st == ST_LOADT) & s_axis_tvalid) tab_mem[cnt[11:0]] <= s_axis_tdata;
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
    reg seed_shift_c; reg [127:0] seed_d1, seed_d2;
    always @(posedge clk) begin
        seed_shift_c <= sd_ov;
        seed_d1 <= sd_state;
        seed_d2 <= seed_d1;
    end

    // ------------------------------------------------------------------ broadcast control registers
    reg lc_v; reg [2:0] lc_k; reg [1:0] lc_mode;
    reg [26:0] k_fourT; reg [31:0] k_csame, k_cdiff; reg [32:0] k_twoTm1, k_twoTm2;
    reg h_clr_c, sp_set_c, out_load_c, out_shift_c;
    reg amode;

    // ------------------------------------------------------------------ groups
    wire [127:0] chain [0:32];
    assign chain[32] = seed_d2;
    wire [31:0] g_exv; wire [95:0] g_exk; wire [255:0] g_exb, g_exs;
    wire [223:0] g_nl; wire [607:0] g_sc; wire [63:0] g_ob;
    wire [43:0] ax_bus; wire [3:0] av_bus, as_bus;

    genvar gb;
    generate for (gb = 0; gb < 32; gb = gb + 1) begin : grp
        wire [63:0] jd;
        genvar gj;
        for (gj = 0; gj < 64; gj = gj + 1) begin : jcol
            assign jd[gj] = rowfull[gb + 32 * gj];
        end
        sca_group #(.B(gb), .N(N)) u (
            .clk(clk), .lc_v(lc_v), .lc_k(lc_k), .lc_mode(lc_mode),
            .k_fourT(k_fourT), .k_csame(k_csame), .k_cdiff(k_cdiff), .k_twoTm1(k_twoTm1), .k_twoTm2(k_twoTm2),
            .seed_shift(seed_shift_c), .seed_in(chain[gb + 1]), .seed_out(chain[gb]),
            .h_clr(h_clr_c), .sp_set(sp_set_c),
            .ax(ax_bus), .av(av_bus), .asg(as_bus), .amode(amode),
            .jwe(j_we_c), .jaddr(j_addr_c), .jdata(jd),
            .ex_v(g_exv[gb]), .ex_k(g_exk[3 * gb +: 3]), .exb(g_exb[8 * gb +: 8]), .exs(g_exs[8 * gb +: 8]),
            .nl_acc(g_nl[7 * gb +: 7]), .sc_acc(g_sc[19 * gb +: 19]),
            .out_load(out_load_c), .out_shift(out_shift_c), .out_bits(g_ob[2 * gb +: 2]));
    end endgenerate

    // ------------------------------------------------------------------ extractors
    wire [3:0] ex_idle;
    wire in_v = g_exv[0];
    wire [2:0] in_k = g_exk[2:0];
    generate for (gb = 0; gb < 4; gb = gb + 1) begin : exu
        wire [63:0] wb, ws;
        genvar q, m;
        for (q = 0; q < 8; q = q + 1) begin : qq
            for (m = 0; m < 8; m = m + 1) begin : mm
                assign wb[8 * m + q] = g_exb[8 * (gb + 4 * q) + m];
                assign ws[8 * m + q] = g_exs[8 * (gb + 4 * q) + m];
            end
        end
        sca_extractor #(.E(gb)) u (
            .clk(clk), .rst(~rst_n), .in_v(in_v), .in_k(in_k), .in_b(wb), .in_s(ws),
            .o_v(av_bus[gb]), .o_x(ax_bus[11 * gb +: 11]), .o_s(as_bus[gb]), .idle(ex_idle[gb]));
    end endgenerate

    // words received in the current batch (counted one cycle after the slot write, so "pending" already includes it)
    reg [3:0] words_seen; reg in_v_d;
    reg batch_start;   // pulse: clears words_seen
    always @(posedge clk) begin
        in_v_d <= in_v;
        if (batch_start) words_seen <= 4'd0;
        else if (in_v_d) words_seen <= words_seen + 4'd1;
    end
    wire quiet = (words_seen == 4'd8) & (&ex_idle) & ~in_v & ~in_v_d;

    // flips (DEC mode only)
    reg [31:0] flips_total;
    wire [2:0] nv_issue = av_bus[0] + av_bus[1] + av_bus[2] + av_bus[3];

    // ------------------------------------------------------------------ capture: n_lin / score reduction
    reg [7:0] cap_sr; reg [15:0] capm_sr;    // token + mode, 8 cycles after the round-7 command
    reg cap1, cap2; reg [1:0] cap1_m, cap2_m;
    reg [9:0] nl1 [0:3]; reg signed [21:0] sc1 [0:3];
    reg [11:0] nl2; reg signed [23:0] sc2;
    integer ii, gg;
    always @(posedge clk) begin
        cap1 <= cap_sr[7]; cap1_m <= capm_sr[15:14];
        for (ii = 0; ii < 4; ii = ii + 1) begin
            nl1[ii] <= g_nl[7 * (8 * ii + 0) +: 7] + g_nl[7 * (8 * ii + 1) +: 7] + g_nl[7 * (8 * ii + 2) +: 7] +
                       g_nl[7 * (8 * ii + 3) +: 7] + g_nl[7 * (8 * ii + 4) +: 7] + g_nl[7 * (8 * ii + 5) +: 7] +
                       g_nl[7 * (8 * ii + 6) +: 7] + g_nl[7 * (8 * ii + 7) +: 7];
            sc1[ii] <= $signed(g_sc[19 * (8 * ii + 0) +: 19]) + $signed(g_sc[19 * (8 * ii + 1) +: 19]) +
                       $signed(g_sc[19 * (8 * ii + 2) +: 19]) + $signed(g_sc[19 * (8 * ii + 3) +: 19]) +
                       $signed(g_sc[19 * (8 * ii + 4) +: 19]) + $signed(g_sc[19 * (8 * ii + 5) +: 19]) +
                       $signed(g_sc[19 * (8 * ii + 6) +: 19]) + $signed(g_sc[19 * (8 * ii + 7) +: 19]);
        end
        cap2 <= cap1; cap2_m <= cap1_m;
        nl2 <= nl1[0] + nl1[1] + nl1[2] + nl1[3];
        sc2 <= sc1[0] + sc1[1] + sc1[2] + sc1[3];
    end

    // next-step table entry and the correction pipeline: corr = int32(((nlin * kcorr) >> 8) + kconst)
    reg [31:0] tn_fourT, tn_q, tn_kconst, tn_kcorr;
    reg tab_ld3, tab_ld2, tab_ld1;   // table read in flight: address with tab_ld1, data in tab_q with tab_ld3
    always @(posedge clk) begin
        tab_ld2 <= tab_ld1; tab_ld3 <= tab_ld2;
        if (tab_ld3) begin
            tn_fourT <= tab_q[31:0]; tn_q <= tab_q[63:32]; tn_kconst <= tab_q[95:64]; tn_kcorr <= tab_q[127:96];
        end
    end
    // tab_ld1 is a one-cycle pulse issued with tab_raddr (set in the FSM); the data is in tab_q two cycles later
    reg c_go1, c_go2, c_go3, c_go4;
    reg [11:0] c_nl; reg signed [31:0] c_kc;
    reg signed [44:0] c_m, c_p;
    reg [31:0] c_corr;
    reg nk_ready; reg [26:0] nk_fourT; reg [31:0] nk_csame, nk_cdiff; reg [32:0] nk_twoTm1, nk_twoTm2;
    wire [31:0] tn_twoT = {1'b0, tn_fourT[31:1]};
    always @(posedge clk) begin
        c_go1 <= cap2 & (cap2_m != M_SCORE);
        c_nl <= (cap2_m == M_INIT) ? 12'd0 : nl2;
        c_kc <= $signed(tn_kcorr);
        c_go2 <= c_go1; c_m <= $signed({1'b0, c_nl}) * c_kc;
        c_go3 <= c_go2; c_p <= c_m;
        c_go4 <= c_go3; c_corr <= c_p[39:8] + tn_kconst;        // (c_p >>> 8)[31:0] + kconst, int32 wrap
    end

    // score -> cut = (sumw + shsum/2)/2 with C truncation toward zero
    reg cut_go1, cut_go2; reg signed [63:0] cut_half, cut_val;
    wire signed [63:0] sh64 = {{40{sc2[23]}}, sc2};
    wire signed [63:0] cut_sum = p_sumw + cut_half;
    always @(posedge clk) begin
        cut_go1 <= cap2 & (cap2_m == M_SCORE);
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
    reg [7:0] oc;             // output sequencer
    reg [1:0] ocad;           // 3-cycle cadence for spin words
    reg [15:0] ot;            // trace word index
    reg trace_rd; reg trace_rd_d;

    always @(posedge clk) begin
        // defaults
        lc_v <= 1'b0; batch_start <= 1'b0; tab_ld1 <= 1'b0; trace_we <= 1'b0;
        out_load_c <= 1'b0; out_shift_c <= 1'b0; of_push <= 1'b0; trace_rd <= 1'b0; sd_run <= 1'b0;
        cap_sr <= {cap_sr[6:0], 1'b0}; capm_sr <= {capm_sr[13:0], 2'b00};
        trace_rd_d <= trace_rd;
        if (c_go4) begin
            nk_csame <= tn_q - c_corr; nk_cdiff <= tn_q + c_corr;
            nk_fourT <= tn_fourT[26:0];
            nk_twoTm1 <= {1'b0, tn_twoT} - 33'd1;
            nk_twoTm2 <= {tn_twoT, 1'b0} - 33'd2;
            nk_ready <= 1'b1;
        end
        if (cap2 & (cap2_m == M_DEC)) begin trace_we <= 1'b1; trace_wa <= t_step[11:0] - 12'd1; trace_wd <= nl2; end
        if (amode & (st != ST_SEED)) flips_total <= flips_total + {29'd0, nv_issue};

        if (~rst_n) begin
            st <= ST_HDR0; cnt <= 16'd0; nk_ready <= 1'b0; cap_sr <= 8'd0; h_clr_c <= 1'b0; sp_set_c <= 1'b0;
            amode <= 1'b0;
        end else case (st)
        ST_HDR0: if (in_fire) begin
            p_trials <= s_axis_tdata[15:0]; p_S <= s_axis_tdata[31:16];
            p_loadJ <= s_axis_tdata[32]; p_trace <= s_axis_tdata[33]; p_seed <= s_axis_tdata[127:64];
            st <= ST_HDR1;
        end
        ST_HDR1: if (in_fire) begin
            p_toff <= s_axis_tdata[31:0]; p_sumw <= s_axis_tdata[127:64];
            cnt <= 16'd0;
            st <= p_loadJ ? ST_LOADJ : (p_S != 16'd0 ? ST_LOADT : ST_TRIAL);
            tr <= 16'd0;
        end
        ST_LOADJ: if (in_fire) begin
            cnt <= cnt + 16'd1;
            if (cnt == 16'd32767) begin cnt <= 16'd0; st <= (p_S != 16'd0) ? ST_LOADT : ST_TRIAL; end
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
            // command -> command no earlier than quiet + 3 (DRAIN).
            if ((dr != 4'd0 || quiet) && dr != 4'd15) dr <= dr + 4'd1;
            if (dr >= DRAIN[3:0]) begin
                if ((b_mode == M_INIT && p_S == 16'd0) || (b_mode == M_DEC && t_step == p_S)) begin
                    b_mode <= M_SCORE; iss <= 4'd0; batch_start <= 1'b1; st <= ST_BATCH;
                end else if (nk_ready) begin
                    // start DEC step t_step (0 after INIT)
                    if (b_mode == M_INIT) t_step <= 16'd1; else t_step <= t_step + 16'd1;
                    k_fourT <= nk_fourT; k_csame <= nk_csame; k_cdiff <= nk_cdiff; k_twoTm1 <= nk_twoTm1;
                    k_twoTm2 <= nk_twoTm2; nk_ready <= 1'b0;
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
// Flip extractor e: slot k holds the 64 flip bits of round k (bit p = 8m + q -> spin x = {k, p, e}) and the new spins.
module sca_extractor #(
    parameter integer E = 0
) (
    input  wire        clk,
    input  wire        rst,
    input  wire        in_v,
    input  wire [2:0]  in_k,
    input  wire [63:0] in_b,
    input  wire [63:0] in_s,
    output reg         o_v,
    output reg  [10:0] o_x,
    output reg         o_s,
    output wire        idle
);
    reg [63:0] slot_b [0:7];
    reg [63:0] slot_s [0:7];
    reg [7:0] pend;
    reg [63:0] cur, curs; reg [2:0] curk;
    wire [63:0] lb = cur & (~cur + 64'd1);
    wire [63:0] rest = cur & ~lb;
    wire last = (rest == 64'd0);
    wire [2:0] kn = pend[0] ? 3'd0 : pend[1] ? 3'd1 : pend[2] ? 3'd2 : pend[3] ? 3'd3 :
                    pend[4] ? 3'd4 : pend[5] ? 3'd5 : pend[6] ? 3'd6 : 3'd7;
    wire take = last & (pend != 8'd0);
    function [5:0] enc64(input [63:0] oh);
        enc64 = {|(oh & 64'hFFFF_FFFF_0000_0000), |(oh & 64'hFFFF_0000_FFFF_0000), |(oh & 64'hFF00_FF00_FF00_FF00),
                 |(oh & 64'hF0F0_F0F0_F0F0_F0F0), |(oh & 64'hCCCC_CCCC_CCCC_CCCC), |(oh & 64'hAAAA_AAAA_AAAA_AAAA)};
    endfunction
    always @(posedge clk) begin
        o_v <= (cur != 64'd0);
        o_x <= {curk, enc64(lb), E[1:0]};
        o_s <= |(curs & lb);
        if (rst) begin
            cur <= 64'd0; pend <= 8'd0;
        end else begin
            if (!last) cur <= rest;
            else if (take) begin cur <= slot_b[kn]; curs <= slot_s[kn]; curk <= kn; end
            else cur <= 64'd0;
            pend <= (pend & ~(take ? (8'd1 << kn) : 8'd0)) | ((in_v && in_b != 64'd0) ? (8'd1 << in_k) : 8'd0);
        end
        if (in_v) begin slot_b[in_k] <= in_b; slot_s[in_k] <= in_s; end
    end
    assign idle = (cur == 64'd0) & (pend == 8'd0);
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
    input  wire [32:0]  k_twoTm1,
    input  wire [32:0]  k_twoTm2,
    input  wire         seed_shift,
    input  wire [127:0] seed_in,
    output wire [127:0] seed_out,
    input  wire         h_clr,
    input  wire         sp_set,
    input  wire [43:0]  ax,
    input  wire [3:0]   av,
    input  wire [3:0]   asg,
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

    // ---- one register hop for every broadcast input
    reg g_v; reg [2:0] g_k; reg [1:0] g_mode;
    reg [26:0] r_fourT; reg [31:0] r_csame, r_cdiff; reg [32:0] r_twoTm1, r_twoTm2;
    reg [10:0] ga0, ga1, ga2, ga3; reg [3:0] gv, gs; reg gm;
    reg r_hclr, r_spset, r_seedsh, r_outload, r_outshift;
    reg r_jwe; reg [10:0] r_jaddr; reg [63:0] r_jdata;
    always @(posedge clk) begin
        g_v <= lc_v; g_k <= lc_k; g_mode <= lc_mode;
        r_fourT <= k_fourT; r_csame <= k_csame; r_cdiff <= k_cdiff; r_twoTm1 <= k_twoTm1; r_twoTm2 <= k_twoTm2;
        ga0 <= ax[10:0]; ga1 <= ax[21:11]; ga2 <= ax[32:22]; ga3 <= ax[43:33]; gv <= av; gs <= asg; gm <= amode;
        r_hclr <= h_clr; r_spset <= sp_set; r_seedsh <= seed_shift; r_outload <= out_load; r_outshift <= out_shift;
        r_jwe <= jwe; r_jaddr <= jaddr; r_jdata <= jdata;
    end

    // ---- state
    reg [11:0] h [0:63];        // stored field = true field - s (|value| <= 2000)
    reg signed [4:0] pend [0:63];
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
    wire run_rng = g_v & (g_mode != M_SCORE);
    always @(posedge clk) begin
        for (m = 0; m < 8; m = m + 1) begin
            if (r_seedsh) rs[m] <= (m == 7) ? seed_in : rs[(m + 1) & 7];
            else if (run_rng) rs[m] <= rng_next(rs[m]);
        end
    end

    // ---- lane pipeline
    reg v1, v2, v3, v4, v5; reg [2:0] k1, k2, k3, k4, k5; reg [1:0] md1, md2, md3, md4, md5;
    reg [11:0] hs1 [0:7]; reg [7:0] si1, spi1; reg [15:0] u1 [0:7];
    reg signed [12:0] sh2 [0:7]; reg [7:0] sel2, b2; reg signed [26:0] da2 [0:7]; reg signed [15:0] db2 [0:7];
    reg [31:0] z3 [0:7]; reg signed [12:0] sh3 [0:7]; reg [7:0] b3; reg signed [42:0] mreg [0:7];
    reg [31:0] z4 [0:7]; reg [32:0] zp4 [0:7]; reg signed [12:0] sh4 [0:7]; reg [7:0] b4; reg signed [42:0] preg [0:7];
    reg [7:0] fl5, nl5, b5, val5; reg signed [12:0] sh5 [0:7];

    wire [7:0] val4;
    genvar gm2;
    generate for (gm2 = 0; gm2 < 8; gm2 = gm2 + 1) begin : vl
        assign val4[gm2] = (B + 32 * (8 * k4 + gm2)) < N;
    end endgenerate

    reg signed [12:0] sh_c; reg [15:0] zhi_c; reg [15:0] zlo_c; reg signed [26:0] thr_c; reg [32:0] zp_c;
    reg lt_c, nlb_c;
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
        // stage 2 -> 3
        v3 <= v2; k3 <= k2; md3 <= md2;
        for (m = 0; m < 8; m = m + 1) begin
            zhi_c = {{3{sh2[m][12]}}, sh2[m]} + (sel2[m] ? r_csame[31:16] : r_cdiff[31:16]);
            zlo_c = sel2[m] ? r_csame[15:0] : r_cdiff[15:0];
            z3[m] <= {zhi_c, zlo_c};
            sh3[m] <= sh2[m];
            b3[m] <= b2[m];
            mreg[m] <= da2[m] * db2[m];
        end
        // stage 3 -> 4
        v4 <= v3; k4 <= k3; md4 <= md3;
        for (m = 0; m < 8; m = m + 1) begin
            z4[m] <= z3[m];
            zp4[m] <= {z3[m][31], z3[m]} + r_twoTm1;
            sh4[m] <= sh3[m];
            b4[m] <= b3[m];
            preg[m] <= mreg[m];
        end
        // stage 4 -> 5
        v5 <= v4; k5 <= k4; md5 <= md4;
        for (m = 0; m < 8; m = m + 1) begin
            thr_c = preg[m] >>> 16;
            lt_c = $signed({z4[m][31], z4[m]}) < $signed({{6{thr_c[26]}}, thr_c});
            nlb_c = ~zp4[m][32] & ($signed(zp4[m]) <= $signed(r_twoTm2));   // -twoT < z < twoT (exact for twoT = 0)
            fl5[m] <= val4[m] & lt_c;
            nl5[m] <= val4[m] & nlb_c;
            val5[m] <= val4[m];
            b5[m] <= b4[m];
            sh5[m] <= sh4[m];
        end
    end

    // ---- stage 5: write-back, extraction words, accumulators; spin/output registers
    reg [3:0] nl_pc; reg signed [18:0] sc_sum;
    always @(posedge clk) begin
        ex_v <= v5 & (md5 != M_SCORE);
        ex_k <= k5;
        if (g_v && g_k == 3'd0) begin nl_acc <= 7'd0; sc_acc <= 19'd0; end
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
        end else if (v5 && md5 == M_SCORE) begin
            sc_sum = 19'sd0;
            for (m = 0; m < 8; m = m + 1)
                if (val5[m]) sc_sum = sc_sum + {{6{sh5[m][12]}}, sh5[m]};
            sc_acc <= sc_acc + sc_sum;
        end
        if (r_spset) sp <= {64{1'b1}};
        else if (r_outload) sp <= s;
        else if (r_outshift) sp <= {2'b00, sp[63:2]};
    end

    // ---- couplings: two copies x two read ports; port A of each copy also writes during the load
    reg [63:0] mem0 [0:2047];
    reg [63:0] mem1 [0:2047];
    reg [63:0] q0a_r, q0b_r, q1a_r, q1b_r, q0a, q0b, q1a, q1b;
    wire [10:0] a0 = r_jwe ? r_jaddr : ga0;
    wire [10:0] a2 = r_jwe ? r_jaddr : ga2;
    always @(posedge clk) begin
        if (r_jwe) mem0[a0] <= r_jdata;
        q0a_r <= mem0[a0];
        q0a <= q0a_r;
    end
    always @(posedge clk) begin
        q0b_r <= mem0[ga1];
        q0b <= q0b_r;
    end
    always @(posedge clk) begin
        if (r_jwe) mem1[a2] <= r_jdata;
        q1a_r <= mem1[a2];
        q1a <= q1a_r;
    end
    always @(posedge clk) begin
        q1b_r <= mem1[ga3];
        q1b <= q1b_r;
    end

    // ---- field update: control aligned with the coupling words (two register stages after ga*)
    reg [3:0] c3v, c3s; reg c3m;
    reg [3:0] P4, Q4; reg [3:0] Bc4; reg dec4;
    always @(posedge clk) begin
        c3v <= gv; c3s <= gs; c3m <= gm;
        P4 <= c3v & c3s; Q4 <= c3v & ~c3s; dec4 <= c3m;
        Bc4 <= ({3'd0, c3v[0]} + {3'd0, c3v[1]} + {3'd0, c3v[2]} + {3'd0, c3v[3]}) << c3m;
    end
    reg [2:0] pc_c; reg [4:0] pv_c;
    always @(posedge clk) begin
        for (j = 0; j < 64; j = j + 1) begin
            pc_c = {2'd0, q0a[j] ? P4[0] : Q4[0]} + {2'd0, q0b[j] ? P4[1] : Q4[1]} +
                   {2'd0, q1a[j] ? P4[2] : Q4[2]} + {2'd0, q1b[j] ? P4[3] : Q4[3]};
            pv_c = (dec4 ? {pc_c, 2'b00} : {1'b0, pc_c, 1'b0}) - {1'b0, Bc4};
            pend[j] <= pv_c;
            h[j] <= r_hclr ? 12'd0 : h[j] + {{7{pend[j][4]}}, pend[j]};
        end
    end
endmodule
`default_nettype wire
