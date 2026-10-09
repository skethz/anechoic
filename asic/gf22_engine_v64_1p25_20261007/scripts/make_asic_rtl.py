#!/usr/bin/env python3
"""Derive rtl/sca_core_asic.v (ASIC memories) from the unmodified v6.4b RTL (orig/sca_core_v64b.v).
Every edit is an exact, asserted string replacement; all non-memory logic is left as is."""
import re
import sys

src_path, dst_path = sys.argv[1], sys.argv[2]
s = open(src_path).read()


def rep(old, new, count=1):
    global s
    n = s.count(old)
    assert n == count, (n, old[:80])
    s = s.replace(old, new)


HEADER = """// sca_core_asic.v -- ASIC variant (GF 22FDX estimate) of the V80 Onsager-SCA engine core v6.4b (sca_core.v, 8 extractors).
// Generated from the unmodified v6.4b source by scripts/make_asic_rtl.py; only the memories change (plus removal of the
// FPGA-only attributes). Cycle- and bit-exact with v6.4b at the stream interface (verified, see ../sim).
//
// Memory changes (v6.4b -> ASIC):
//  * Couplings. v6.4b keeps, per group, four copies (2 BRAM + 2 URAM) of a 2048 x 64 bank with two read ports each: port
//    e (e = 0..7) reads the row of extractor e. Extractor e drains groups e, e+8, e+16, e+24 and its flipped spin is
//    x = {k, p, e} (o_x[2:0] = E), so port e only ever reads rows x with x mod 8 == e. Here every port has its own
//    single-port residue bank holding exactly those rows at address x >> 3: bank e = 256 rows x 2048 bits (64 columns
//    per group), 4 Mbit per engine in total (v6.4b: 4 copies = 16 Mbit). All 32 groups read the same row in the same
//    cycle, so the slices of CPL_G consecutive groups share one macro of 256 x (64 * CPL_G) bits. The load writes row x
//    into bank x mod 8 at address x >> 3. Read latency is unchanged: address register (ga* in v6.4b, ca here) -> macro
//    (flow-through; it takes the place of the read register q*_r) -> second register q* in the group -> field update.
//    A bank is enabled only when its extractor issues a valid flip (or for a load write); the group ignores the coupling
//    word of a port without a valid flip (P4[e] = Q4[e] = 0), so the gated reads are exact.
//  * Table memory (4096 x 128) and n_lin trace (4096 x 12): single-port SRAM macros, reads enabled only when the data is
//    consumed (tab_ld1, trace_rd); writes never coincide with these reads (checked in simulation with SCA_SIM_CHECKS).
//  * Other arrays stay in flip-flops (output FIFO, extractor slots, group state, row buffer).
// Memory macros are instantiated through mem_cpl / mem_tab / mem_trc (behavioral models in sca_mem_behav.v; the
// compiler-macro wrappers used for synthesis exist only on the synthesis server).
//
// ---------------- v6.4b header (unchanged) ----------------
"""
s = HEADER + s

# FPGA-only attributes (X_INTERFACE_*, max_fanout, dont_touch)
s = re.sub(r'\(\*[^*]*\*\)\s*', '', s)

# ---- top: parameter for the macro grouping
rep("""    parameter integer DRAIN = 2   // cycles between "extractors quiet" and the next round command (>= 2; 1 fails in simulation)
) (""", """    parameter integer DRAIN = 2,  // cycles between "extractors quiet" and the next round command (>= 2; 1 fails in simulation)
    parameter integer CPL_G = 4   // ASIC: groups per coupling macro (1, 2 or 4); macro = 256 words x (64 * CPL_G) bits
) (""")
rep("""    localparam integer SEED_LAT = 15;   // sd_run (aligned with sd_xa) -> sca_mix output
""", """    localparam integer SEED_LAT = 15;   // sd_run (aligned with sd_xa) -> sca_mix output
    localparam integer CPL_W = 64 * CPL_G;   // ASIC: coupling macro width
    localparam integer CPL_M = 32 / CPL_G;   // ASIC: macros per residue bank
""")

# ---- top: tables and trace -> single-port macros
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
    reg tab_ld3, tab_ld2, tab_ld1, tab_ld4;    // table read in flight: address with tab_ld1, data in tab_q with tab_ld3
    reg trace_rd; reg trace_rd_d;
    reg [11:0] tab_raddr; reg [127:0] tab_q;
    wire [127:0] tab_q1;                       // macro output = v6.4b's read register tab_q1
    wire tab_we = (st == ST_LOADT) & in_v;
    mem_tab u_tab (.clk(clk), .ce(tab_we | tab_ld1), .we(tab_we), .addr(tab_we ? cnt[11:0] : tab_raddr),
                   .wdata(in_d), .rdata(tab_q1));
    always @(posedge clk) tab_q <= tab_q1;
    reg trace_we; reg [11:0] trace_wa; reg [11:0] trace_wd;
    reg [11:0] trace_ra; wire [11:0] trace_q;  // macro output = v6.4b's read register trace_q
    mem_trc u_trc (.clk(clk), .ce(trace_we | trace_rd), .we(trace_we), .addr(trace_we ? trace_wa : trace_ra),
                   .wdata(trace_wd), .rdata(trace_q));
""")
rep("    reg tab_ld3, tab_ld2, tab_ld1, tab_ld4;    // table read in flight: address with tab_ld1, data in tab_q with tab_ld3\n    wire [31:0] tn_twoT",
    "    wire [31:0] tn_twoT")
rep("    reg [15:0] ot;            // trace word index\n    reg trace_rd; reg trace_rd_d;\n",
    "    reg [15:0] ot;            // trace word index\n")

# ---- top: residue banks in front of the groups; groups get cq instead of the coupling address/load ports
rep("""    wire [87:0] ax_bus; wire [7:0] av_bus, as_bus;

    genvar gb;
    generate for (gb = 0; gb < 32; gb = gb + 1) begin : grp
        wire [63:0] jd;
        genvar gj;
        for (gj = 0; gj < 64; gj = gj + 1) begin : jcol
            assign jd[gj] = rowbuf[gb + 32 * gj];
        end
        sca_group #(.B(gb), .N(N)) u (""", """    wire [87:0] ax_bus; wire [7:0] av_bus, as_bus;

    // ------------------------------------------------------------------ ASIC: coupling residue banks (see header)
    // Write data: group b's 64 columns of the completed row (bit j = J[x][b + 32j]), registered like v6.4b's r_jdata.
    wire [2047:0] cw_next;
    genvar gb, gj;
    generate for (gb = 0; gb < 32; gb = gb + 1) begin : cwm
        for (gj = 0; gj < 64; gj = gj + 1) begin : jcol
            assign cw_next[64 * gb + gj] = rowbuf[gb + 32 * gj];
        end
    end endgenerate
    reg [2047:0] cw_d;
    always @(posedge clk) cw_d <= cw_next;
    wire [16383:0] cq;                         // bank e, group b: cq[2048 * e + 64 * b +: 64]
    genvar ge, gi;
    generate for (ge = 0; ge < 8; ge = ge + 1) begin : cbank
        for (gi = 0; gi < CPL_M; gi = gi + 1) begin : cm
            // per-macro copies of the address / enable registers (v6.4b: ga0..ga7 and r_jwe in every group)
            reg [7:0] ca; reg cce, cwe;
            always @(posedge clk) begin
                ca  <= j_we_r ? j_addr_r[10:3] : ax_bus[11 * ge + 3 +: 8];
                cwe <= j_we_r & (j_addr_r[2:0] == ge);
                cce <= j_we_r ? (j_addr_r[2:0] == ge) : av_bus[ge];
            end
            mem_cpl #(.W(CPL_W)) u (.clk(clk), .ce(cce), .we(cwe), .addr(ca),
                .wdata(cw_d[CPL_W * gi +: CPL_W]), .rdata(cq[2048 * ge + CPL_W * gi +: CPL_W]));
        end
    end endgenerate

    generate for (gb = 0; gb < 32; gb = gb + 1) begin : grp
        sca_group #(.B(gb), .N(N)) u (""")
rep("""            .ax(ax_bus), .av(av_bus), .asg(as_bus), .amode(amode),
            .jwe(j_we_r), .jaddr(j_addr_r), .jdata(jd),""", """            .av(av_bus), .asg(as_bus), .amode(amode),
            .cq({cq[2048 * 7 + 64 * gb +: 64], cq[2048 * 6 + 64 * gb +: 64], cq[2048 * 5 + 64 * gb +: 64],
                 cq[2048 * 4 + 64 * gb +: 64], cq[2048 * 3 + 64 * gb +: 64], cq[2048 * 2 + 64 * gb +: 64],
                 cq[2048 * 1 + 64 * gb +: 64], cq[64 * gb +: 64]}),""")

# ---- simulation checks (single-port legality)
rep("""        default: st <= ST_HDR0;
        endcase
    end
endmodule""", """        default: st <= ST_HDR0;
        endcase
    end

`ifdef SCA_SIM_CHECKS
    // single-port legality: a memory never sees a read request and a write in the same cycle
    integer chk_err = 0;
    always @(posedge clk) if (rst_n) begin
        if (j_we_r & (|av_bus)) begin chk_err = chk_err + 1; $display("CHECK: coupling load write during a flip read"); end
        if (tab_we & tab_ld1) begin chk_err = chk_err + 1; $display("CHECK: table write during a table read"); end
        if (trace_we & trace_rd) begin chk_err = chk_err + 1; $display("CHECK: trace write during a trace read"); end
    end
`endif
endmodule""")

# ---- group ports
rep("""    input  wire [87:0]  ax,
    input  wire [7:0]   av,""", """    input  wire [7:0]   av,""")
rep("""    input  wire         amode,
    input  wire         jwe,
    input  wire [10:0]  jaddr,
    input  wire [63:0]  jdata,""", """    input  wire         amode,
    input  wire [511:0] cq,          // ASIC: this group's 64-bit slice of the 8 residue-bank outputs {bank7, ..., bank0}""")
rep("""    reg [10:0] ga0, ga2, ga4, ga6;   // port A addresses (coupling write address during the load)
    reg [10:0] ga1, ga3, ga5, ga7;
""", "")
rep("""    reg r_jwe; reg [63:0] r_jdata;
""", "")
rep("""        ga0 <= jwe ? jaddr : ax[10:0]; ga1 <= ax[21:11]; ga2 <= jwe ? jaddr : ax[32:22]; ga3 <= ax[43:33];
        ga4 <= jwe ? jaddr : ax[54:44]; ga5 <= ax[65:55]; ga6 <= jwe ? jaddr : ax[76:66]; ga7 <= ax[87:77];
""", "")
rep("""        r_jwe <= jwe; r_jdata <= jdata;
""", "")

# ---- group coupling memories -> registered bank outputs
i0 = s.index("    // ---- couplings: two copies x two read ports; port A of each copy also writes during the load")
i1 = s.index("    // ---- field update: control aligned with the coupling words (two register stages after ga*)")
s = s[:i0] + """    // ---- couplings (ASIC): the residue-bank outputs (= v6.4b's read registers q*_r) are registered once more (= q*).
    // Port e (bank e) feeds P4[e]/Q4[e] exactly as v6.4b's q0a, q0b, q1a, q1b, q2a, q2b, q3a, q3b (ga0..ga7).
    reg [63:0] q0a, q0b, q1a, q1b, q2a, q2b, q3a, q3b;
    always @(posedge clk) begin
        q0a <= cq[63:0];    q0b <= cq[127:64];  q1a <= cq[191:128]; q1b <= cq[255:192];
        q2a <= cq[319:256]; q2b <= cq[383:320]; q3a <= cq[447:384]; q3b <= cq[511:448];
    end

""" + s[i1:]
rep("// Group b: 64 spins/fields y = b + 32j, 8 decision lanes (lane m -> spins j = 8k + m), coupling bank b (two copies).",
    "// Group b: 64 spins/fields y = b + 32j, 8 decision lanes (lane m -> spins j = 8k + m). ASIC: the coupling residue banks\n"
    "// are at the top level; cq carries this group's 64-bit slice of the eight bank outputs.")
open(dst_path, 'w').write(s)
print('written', dst_path, len(s.splitlines()), 'lines')
