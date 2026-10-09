#!/usr/bin/env python3
"""Derive rtl/sca_core_mb_asic.v from the unmodified multi-bit FPGA core orig/sca_core_mb_r2.v (MULTIBIT_SPEC.md):
ASIC memories in the 1.25 GHz organisation of the v6.4b ASIC (single-copy residue banks as 256 x 128 single-port SRAM
macros, two-cycle step-table read), plus the external-field bias (2048 x int16 bias memory, biased field initialisation),
a runtime number of active spins n (sca_ref_bias.hpp run_trial_bias) and a raw-score output mode.
Every edit is an exact, asserted string replacement; the decision/field-update datapath is the sca_core_mb_r2.v text
except for the field width parameters and the field register's bias-load input.
Usage: make_mb_asic_rtl.py orig/sca_core_mb_r2.v rtl/sca_core_mb_asic.v"""
import re
import sys

src_path, dst_path = sys.argv[1], sys.argv[2]
s = open(src_path).read()


def rep(old, new, count=1):
    global s
    n = s.count(old)
    assert n == count, (n, old[:90])
    s = s.replace(old, new)


HEADER = """// sca_core_mb_asic.v -- ASIC variant (GF 22FDX estimate) of the multi-bit V80 engine core (sca_core_mb_r2.v, K-bit
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
"""
s = HEADER + s

# FPGA-only attributes (X_INTERFACE_*, max_fanout, dont_touch, ram_style)
s = re.sub(r'\(\*[^*]*\*\)\s*', '', s)

# ---------------------------------------------------------------- top: name, parameters, widths
rep("module sca_core_mb #(", "module sca_core #(")
rep("""    parameter integer DRAIN = 2   // cycles between "extractors quiet" and the next round command (>= 2; 1 fails in simulation)
) (""", """    parameter integer DRAIN = 2,  // cycles between "extractors quiet" and the next round command (>= 2; 1 fails in simulation)
    parameter integer CPL_G = 2   // ASIC: groups per coupling macro; macro = 256 words x (64 * CPL_G) bits
) (""")
rep("""    localparam integer HB = 1 + $clog2(2047 * MJ + 1);     // field width
    localparam integer SW = HB + 1;                         // s*h
    localparam integer ZH = (K <= 4) ? 16 : HB + 2;         // high part of z""",
    """    localparam integer HB = 1 + $clog2(2047 * MJ + 32768 + 1);   // field width with a 16-bit bias (17 at K = 2, 4)
    localparam integer SW = HB + 1;                         // s*h
    localparam integer ZH = HB + 1;                         // high part of z: exact for any field and any int32 C
    localparam integer CPL_W = 64 * CPL_G;                  // ASIC: coupling macro width
    localparam integer CPL_M = 32 / CPL_G;                  // ASIC: macros per residue-bank plane
    localparam [11:0] N12 = N;""")
rep("""                     ST_SEED = 4'd5, ST_BATCH = 4'd6, ST_WAIT = 4'd7, ST_SCOREW = 4'd8, ST_OUT = 4'd9;""",
    """                     ST_SEED = 4'd5, ST_BATCH = 4'd6, ST_WAIT = 4'd7, ST_SCOREW = 4'd8, ST_OUT = 4'd9,
                     ST_LOADB = 4'd10;""")
rep("""    reg [15:0] p_trials; reg [15:0] p_S; reg p_loadJ, p_trace, p_keepT; reg [63:0] p_seed; reg [31:0] p_toff;""",
    """    reg [15:0] p_trials; reg [15:0] p_S; reg p_loadJ, p_trace, p_keepT; reg [63:0] p_seed; reg [31:0] p_toff;
    reg p_loadB, p_useB, p_raw; reg [11:0] p_n;      // bias segment, biased init, raw score; active spins""")
rep("""    wire core_ready = (st == ST_HDR0) | (st == ST_HDR1) | (st == ST_LOADJ) | (st == ST_LOADT);""",
    """    wire core_ready = (st == ST_HDR0) | (st == ST_HDR1) | (st == ST_LOADJ) | (st == ST_LOADT) | (st == ST_LOADB);""")

# ---------------------------------------------------------------- top: tables and trace -> macros (1.25 GHz form), bias memory
rep("""    // ------------------------------------------------------------------ tables and trace
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
""", """    // ------------------------------------------------------------------ tables and trace (ASIC: single-port SRAM macros)
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
""")
rep("    reg tab_ld3, tab_ld2, tab_ld1, tab_ld4;    // table read in flight: address with tab_ld1, data in tab_q with tab_ld3\n    wire [31:0] tn_twoT",
    "    wire [31:0] tn_twoT")
rep("""        tab_ld2 <= tab_ld1; tab_ld3 <= tab_ld2; tab_ld4 <= tab_ld3;
        if (tab_ld3) begin
            tn_fourT <= tab_q[31:0]; tn_q <= tab_q[63:32]; tn_kconst <= tab_q[95:64]; tn_kcorr <= tab_q[127:96];
        end
        if (tab_ld4) begin tn_U <= tn_twoT - 32'd1; tn_L <= 32'd1 - tn_twoT; end""",
    """        tab_ld2 <= tab_ld1; tab_ld3 <= tab_ld2; tab_ld4 <= tab_ld3; tab_ld5 <= tab_ld4;
        if (tab_ld4) begin   // 1.25 GHz: one cycle later than the FPGA core (still >= 13 cycles before first use)
            tn_fourT <= tab_q[31:0]; tn_q <= tab_q[63:32]; tn_kconst <= tab_q[95:64]; tn_kcorr <= tab_q[127:96];
        end
        if (tab_ld5) begin tn_U <= tn_twoT - 32'd1; tn_L <= 32'd1 - tn_twoT; end""")
rep("    reg [15:0] ot;            // trace word index\n    reg trace_rd; reg trace_rd_d;\n",
    "    reg [15:0] ot;            // trace word index\n")

# ---------------------------------------------------------------- top: residue banks, group wiring
rep("""    wire [87:0] ax_bus; wire [7:0] av_bus, as_bus;

    genvar gb;
    generate for (gb = 0; gb < 32; gb = gb + 1) begin : grp
        wire [63:0] jd;
        genvar gj;
        for (gj = 0; gj < 64; gj = gj + 1) begin : jcol
            assign jd[gj] = rowbuf[gb + 32 * gj];
        end
        sca_group_mb #(.B(gb), .N(N), .K(K), .HB(HB), .SW(SW), .ZH(ZH), .SCG(SCG), .PB(PB)) u (""",
    """    wire [87:0] ax_bus; wire [7:0] av_bus, as_bus;
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
        sca_group_mb #(.B(gb), .N(N), .K(K), .HB(HB), .SW(SW), .ZH(ZH), .SCG(SCG), .PB(PB)) u (""")
rep("""            .ax(ax_bus), .av(av_bus), .asg(as_bus), .amode(amode),
            .jwe(j_we_r), .jaddr(j_addr_r), .jpl(j_pl_r), .jdata(jd),""",
    """            .av(av_bus), .asg(as_bus), .amode(amode), .cqm(cqg), .nact(n_c),
            .bl_v(bl_v2), .bl_a(bl_a2), .bl_d(bl_d2),""")

# ---------------------------------------------------------------- top: score word, output mask, FSM
rep("""                if (of_room) begin of_push <= 1'b1; of_din <= cut_val; oc <= oc + 8'd1; end""",
    """                if (of_room) begin of_push <= 1'b1; of_din <= p_raw ? sh64 : cut_val; oc <= oc + 8'd1; end""")
rep("""            if (oc < 8'd3) oc <= oc + 8'd1;
            else if (oc < 8'd35) begin
                if (ocad == 2'd0) begin
                    if (of_room) begin
                        of_push <= 1'b1;
                        of_din <= (oc == 8'd34) ? (spin_word & 64'h0000_0000_0000_FFFF) : spin_word;
                        out_shift_c <= 1'b1; ocad <= 2'd1; oc <= oc + 8'd1;
                    end""",
    """            if (oc < 8'd3) begin oc <= oc + 8'd1; om_rem <= $signed({1'b0, p_n}); end
            else if (oc < 8'd35) begin
                if (ocad == 2'd0) begin
                    if (of_room) begin
                        of_push <= 1'b1;
                        of_din <= spin_word & om_mask;      // spins >= n are zero
                        out_shift_c <= 1'b1; ocad <= 2'd1; oc <= oc + 8'd1; om_rem <= om_rem - 13'sd64;
                    end""")
rep("""    reg [15:0] ot;            // trace word index

    always @(posedge clk) begin
        // defaults""",
    """    reg [15:0] ot;            // trace word index
    // output spin-word mask: om_rem = n - 64 w for the next word w (set while waiting for out_load), om_mask = bits < om_rem
    reg signed [12:0] om_rem; reg [63:0] om_mask;
    always @(posedge clk)
        om_mask <= (om_rem >= 13'sd64) ? {64{1'b1}} : (om_rem <= 13'sd0) ? 64'd0 : ((64'd1 << om_rem[5:0]) - 64'd1);

    always @(posedge clk) begin
        // defaults""")
rep("""        if (amode & (st != ST_SEED)) flips_total <= flips_total + {28'd0, nv_issue};
""", """        if (amode & (st != ST_SEED)) flips_total <= flips_total + {28'd0, nv_issue};
        if (bl_run) begin bl_a <= bl_a + 8'd1; if (bl_a == 8'd255) bl_run <= 1'b0; end
""")
rep("""            st <= ST_HDR0; cnt <= 16'd0; jcnt <= 20'd0; nk_ready <= 1'b0; cap_sr <= 10'd0; h_clr_c <= 1'b0; sp_set_c <= 1'b0;
            amode <= 1'b0;""",
    """            st <= ST_HDR0; cnt <= 16'd0; jcnt <= 20'd0; nk_ready <= 1'b0; cap_sr <= 10'd0; h_clr_c <= 1'b0; sp_set_c <= 1'b0;
            amode <= 1'b0; bl_run <= 1'b0;""")
rep("""            p_loadJ <= in_d[32]; p_trace <= in_d[33]; p_keepT <= in_d[34]; p_seed <= in_d[127:64];""",
    """            p_loadJ <= in_d[32]; p_trace <= in_d[33]; p_keepT <= in_d[34]; p_seed <= in_d[127:64];
            p_loadB <= in_d[35]; p_useB <= in_d[36]; p_raw <= in_d[37];""")
rep("""            p_toff <= in_d[31:0]; p_sumw <= in_d[127:64];
            cnt <= 16'd0; jcnt <= 20'd0;
            st <= p_loadJ ? ST_LOADJ : ((p_S != 16'd0 && !p_keepT) ? ST_LOADT : ST_TRIAL);""",
    """            p_toff <= in_d[31:0]; p_sumw <= in_d[127:64];
            p_n <= (in_d[43:32] == 12'd0) ? N12 : in_d[43:32];
            cnt <= 16'd0; jcnt <= 20'd0;
            st <= p_loadJ ? ST_LOADJ : p_loadB ? ST_LOADB : ((p_S != 16'd0 && !p_keepT) ? ST_LOADT : ST_TRIAL);""")
rep("""            if (jcnt == JBEATS - 1) begin jcnt <= 20'd0; st <= (p_S != 16'd0 && !p_keepT) ? ST_LOADT : ST_TRIAL; end
        end""",
    """            if (jcnt == JBEATS - 1) begin
                jcnt <= 20'd0; st <= p_loadB ? ST_LOADB : ((p_S != 16'd0 && !p_keepT) ? ST_LOADT : ST_TRIAL);
            end
        end
        ST_LOADB: if (in_fire) begin
            cnt <= cnt + 16'd1;
            if (cnt == 16'd255) begin cnt <= 16'd0; st <= (p_S != 16'd0 && !p_keepT) ? ST_LOADT : ST_TRIAL; end
        end""")
rep("""                h_clr_c <= 1'b1; sp_set_c <= 1'b1; flips_total <= 32'd0; amode <= 1'b0;""",
    """                h_clr_c <= 1'b1; sp_set_c <= 1'b1; flips_total <= 32'd0; amode <= 1'b0;
                bl_run <= p_useB; bl_a <= 8'd0;            // bias words are written into the fields during SEED""")
rep("""        ST_SEED: begin
            // issue 256 seeder inputs, then wait until the last state is in the chain (2 cycles after sd_ov)""",
    """        ST_SEED: begin
            h_clr_c <= 1'b0;   // ASIC: one-cycle field clear (nothing updates the fields during SEED except the bias load)
            // issue 256 seeder inputs, then wait until the last state is in the chain (2 cycles after sd_ov)""")
rep("""            if (cnt == 16'd256 && ~seed_shift_c) begin""",
    """            if (cnt == 16'd256 && ~seed_shift_c && ~bl_busy) begin""")

# ---------------------------------------------------------------- simulation checks
rep("""        default: st <= ST_HDR0;
        endcase
    end
endmodule""", """        default: st <= ST_HDR0;
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
endmodule""")

# ---------------------------------------------------------------- group: ports and registers
rep("""    input  wire [87:0]   ax,
    input  wire [7:0]    av,
    input  wire [7:0]    asg,
    input  wire          amode,
    input  wire          jwe,
    input  wire [10:0]   jaddr,
    input  wire [PB-1:0] jpl,
    input  wire [63:0]   jdata,""",
    """    input  wire [7:0]    av,
    input  wire [7:0]    asg,
    input  wire          amode,
    input  wire [8*64*K-1:0] cqm,   // ASIC: this group's 64-bit slices of the 8 x K residue-bank plane outputs
    input  wire [11:0]   nact,      // runtime number of active spins
    input  wire          bl_v,      // bias word valid (seed phase)
    input  wire [7:0]    bl_a,      // bias word address: biases 8 bl_a .. 8 bl_a + 7
    input  wire [127:0]  bl_d,      // bias word (int16 x 8)""")
rep("""    reg [63:0] r_jdata;
    always @(posedge clk) begin""",
    """    localparam [4:0] B5 = B;
    reg [11:0] r_n; reg r_blv; reg [5:0] r_blj; reg [15:0] r_bld;
    always @(posedge clk) begin
        // runtime n and the bias word: field j of this group is y = B + 32 j, its bias is in word 4 j + B[4:3], slice B[2:0]
        r_n <= nact;
        r_blv <= bl_v & (bl_a[1:0] == B5[4:3]); r_blj <= bl_a[7:2]; r_bld <= bl_d[16 * (B % 8) +: 16];
    end
    always @(posedge clk) begin""")
rep("""        r_jdata <= jdata;
""", "")
rep("""        assign val4[gm2] = (B + 32 * (8 * k4 + gm2)) < N;""",
    """        assign val4[gm2] = (B + 32 * (8 * k4 + gm2)) < r_n;   // runtime n (spins >= n are padding)""")

# ---------------------------------------------------------------- group: residue banks -> registered macro outputs
i0 = s.index("    // ---- couplings: single-copy residue banks (see the file header).")
i1 = s.index("    // ---- field update: control aligned with the coupling words (two register stages after ga)")
s = s[:i0] + """    // ---- couplings (ASIC): the residue-bank plane outputs come from the macros at the top level (the macro output takes
    // the place of the FPGA read register qr) and are registered once more here (= the FPGA output register q), so cq
    // has the FPGA core's latency: address register (ca) -> macro -> cq -> field update.
    reg [8 * 64 * K - 1:0] cq;
    always @(posedge clk) cq <= cqm;

""" + s[i1:]
rep("""            h[j] <= r_hclr ? {HB{1'b0}} : h[j] + {{(HB-PW){pend[j][PW-1]}}, pend[j]};""",
    """            h[j] <= (r_blv && r_blj == j) ? {{(HB-16){r_bld[15]}}, r_bld} :      // bias load (seed phase)
                    r_hclr ? {HB{1'b0}} : h[j] + {{(HB-PW){pend[j][PW-1]}}, pend[j]};""")
rep("    reg [HB-1:0] h [0:63];      // field (J_xx = 0, so the true field; |value| <= 1999 * (2^(K-1) - 1))",
    "    reg [HB-1:0] h [0:63];      // field (J_xx = 0, so the true field: b_y + sum_x J_xy s_x)")
open(dst_path, 'w').write(s)
print('written', dst_path, len(s.splitlines()), 'lines')
